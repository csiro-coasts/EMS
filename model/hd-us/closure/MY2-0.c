/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/closure/MY2-0.c
 *  
 *  Description:
 *  Implements Mellor-Yamada level 2 mixing scheme.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: MY2-0.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#if defined(OS_IRIX) || defined(OS_SOLARIS)
#define HAS_CUBE_ROOT
extern double cbrt(double);     /* Not always proto-typed */
#endif

double get_Lscale(geometry_t *window, window_t *windat, win_priv_t *wincon,
                  double Ri, double N2, double Lo, double q, double z,
                  double sl, double bl, double top, double bot, int kts,
                  int ktb, double z0);
int cg;

/* Constants */

#define KAPPA (0.4)             /* Von Karman's constant */
#define rho0 (1024.0)           /* standard density */
#define A1 (0.92)
#define B1 (16.6)
#define A2 (0.74)
#define B2 (10.1)
#define C1 (0.08)
#define G1 (1.0/3.0 - 2.0* A1 / B1)
#define G2 (B2/B1 + 6.0* A1 / B1)
#define BGC1 (B1*(G1 - C1))
#define Rcrit (G1/(G1+G2))
#define EPS (0.00001)
#define MOLK (0.000002)
#define MOLD (0.0001)

/*-------------------------------------------------------------------*/
/* Routine to initialize the Mellor-Yamada 2.0 mixing scheme         */
/*-------------------------------------------------------------------*/
void closure_MY2_init(geometry_t *geom, /* Sparse global geometry
                                           structure */
                      master_t *master  /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coordinates */

  master->q = d_alloc_1d(geom->sgsiz);
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    master->q[c] = EPS;
    master->Vz[c] = master->vz0;
    master->Kz[c] = master->kz0;
  }
}

/* END closure_MY2_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the Mellor-Yamada 2.0 closure scheme.                             */
/*-------------------------------------------------------------------*/
void closure_MY2(geometry_t *window,  /* Processing window */
                 window_t *windat,  /* Window data structure */
                 win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c3;                /* Sparse coodinate / counter */
  int cs;                       /* 2D sparse coordinate */
  int cb;                       /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */
  double sumq;
  double sumzq;
  double lscale;
  double lasymp;
  double l1, l2;
  double dz;
  double du1, du2;
  double drho;
  double du1dz;
  double du2dz;
  double drhodz;
  double GH;
  double GM;
  double Rgrad;
  double Rflux;
  double SH;
  double SM;
  double q3;
  double *vu1;
  double *vu2;

  vu1 = wincon->w1;
  vu2 = wincon->w2;
  memset(vu1, 0, window->sgsiz * sizeof(double));
  memset(vu2, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate gradients over the whole grid */
  /* Vertical u and v gradients */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    cb = window->bot_t[cc];
    zm1 = window->zm1[c];
    while (c < cb) {
      vu1[c] = windat->u[c] - windat->u[zm1];
      vu2[c] = windat->v[c] - windat->v[zm1];
      c = zm1;
      zm1 = window->zm1[c];
    }
  }

  for (cc = 1; cc <= window->b2_t; cc++) {
    c3 = c = window->nsur_t[cc];
    cs = window->m2d[c];
    cb = window->bot_t[cc];
    zm1 = window->zm1[c];

    /*---------------------------------------------------------------*/
    /* Compute q-moment length scale */
    sumq = 0;
    sumzq = 0;
    while (c < cb) {
      /* SIGMA : multiply by depth */
      dz = (window->cellz[c] - window->cellz[zm1]) * wincon->Ds[cs];
      sumq += wincon->q[c] * dz;
      sumzq += wincon->q[c] * fabs(window->cellz[c] * wincon->Ds[cs]) * dz;
      c = zm1;
      zm1 = window->zm1[c];
    }
    lasymp = 0.3 * sumzq / sumq;

    /*---------------------------------------------------------------*/
    /* Compute cell centred viscosity, Vz, and diffusion, Kz. Use */
    /* the minimum of surface,bottom and q-moment length scales.  */
    c = c3;
    zm1 = window->zm1[c];

    while (c < cb) {

      /* SIGMA : Multiply by depth */
      dz = (window->cellz[c] - window->cellz[zm1]) * wincon->Ds[cs];
      l1 =
        KAPPA * ((window->cellz[c] - window->cellz[cb]) * wincon->Ds[cs] +
                 wincon->z0[cs]);
      l2 =
        KAPPA * ((window->cellz[c3] - window->cellz[c]) * wincon->Ds[cs] +
                 wincon->zs);
      lscale = min(l1, l2);
      lscale = lscale / (1 + lscale / lasymp);
      if (windat->L)
        windat->L[c] = lscale;

      /* Vertical u1 gradient */
      du1 = vu1[c];

      /* Vertical u2 gradient */
      du2 = vu2[c];

      /* Vertical density gradient */
      drho = windat->dens_0[c] - windat->dens_0[zm1];
      du1dz = du1 / dz;
      du2dz = du2 / dz;
      drhodz = (drho / dz) / rho0;

      GH = -wincon->g * drhodz;
      GM = du1dz * du1dz + du2dz * du2dz;

      /* Shear production                                            */
      if (windat->s_prod)
	windat->s_prod[c] = windat->Vz[c] * (du1dz * du1dz + du2dz * du2dz);

      /* Buoyancy production                                         */
      if (windat->b_prod) windat->b_prod[c] = -windat->Kz[c] * GH;

      if (GM > 1.e-10) {
        Rgrad = GH / GM;
        Rflux =
          0.725 * (Rgrad + 0.186 -
                   sqrt(Rgrad * Rgrad - 0.316 * Rgrad + 0.0346));
        if (Rflux < Rcrit) {
          SH = 3 * A2 * (G1 - (G1 + G2) * Rflux) / (1 - Rflux);
          SM = SH * A1 * (BGC1 - (BGC1 + (6 * A1 + 3 * A2)) * Rflux) /
            (A2 * (B1 * G1 - (B1 * (G1 + G2) - 3 * A1) * Rflux));
          q3 =
            B1 * lscale * (windat->Vz[c] * GM +
                           windat->Kz[c] * wincon->g * drhodz);

          if (q3 > 0)
            wincon->q[c] = pow(q3, 1.0 / 3.0);
          else {
            hd_quit
              ("Vz: q3 must be greater than zero at %d (%d %d %d) : %f %f %f %f %f\n",
               c, window->s2i[c], window->s2j[c], window->s2k[c], windat->Vz[c],
               windat->Kz[c], lscale, drhodz, GM);
          }
          windat->Vz[c] = lscale * wincon->q[c] * SM;
          windat->Kz[c] = lscale * wincon->q[c] * SH;
        } else {
          windat->Vz[c] = wincon->vz0;
          windat->Kz[c] = wincon->kz0;
        }
      } else {
        windat->Vz[c] = wincon->vz0;
        windat->Kz[c] = wincon->kz0;
      }
      c = zm1;
      zm1 = window->zm1[c];
    }
  }
}

/* END closure_MY2()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate mixing coefficientsusing the Mellor Yamada   */
/* 2.0 turbulence closure scheme applicable to estuaries.            */
/* Output is Vz, Kz and q for every cell.                            */
/*-------------------------------------------------------------------*/
void closure_MY2_est(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c3;                /* Sparse coodinate / counter */
  int cs;                       /* 2D sparse coordinate */
  int cb;                       /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */
  int kts, ktb;
  double sumq;
  double sumzq;
  double lscale;
  double lasymp;
  double dz;
  double du1, du2;
  double drho;
  double du1dz;
  double du2dz;
  double drhodz;
  double GH;
  double GM;
  double Rgrad;
  double Rflux;
  double SH;
  double SM;
  double q3;
  double top_m, bot_m;
  double *vu1;
  double *vu2;

  vu1 = wincon->w1;
  vu2 = wincon->w2;
  memset(vu1, 0, window->sgsiz * sizeof(double));
  memset(vu2, 0, window->sgsiz * sizeof(double));
  /*-----------------------------------------------------------------*/
  /* Calculate gradients over the whole grid */
  /* Vertical u and v gradients */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    cb = window->bot_t[cc];
    zm1 = window->zm1[c];
    while (c < cb) {
      vu1[c] = windat->u[c] - windat->u[zm1];
      vu2[c] = windat->v[c] - windat->v[zm1];
      c = zm1;
      zm1 = window->zm1[c];
    }
  }

  for (cc = 1; cc <= window->b2_t; cc++) {
    c3 = c = window->nsur_t[cc];
    cs = window->m2d[c];
    cb = window->bot_t[cc];
    zm1 = window->zm1[c];

    /*---------------------------------------------------------------*/
    /* Compute q-moment length scale */
    sumq = 0;
    sumzq = 0;
    while (c < cb) {
      /* SIGMA : multiply by depth */
      dz = (window->cellz[c] - window->cellz[zm1]) * wincon->Ds[cs];
      sumq += wincon->q[c] * dz;
      sumzq += wincon->q[c] * fabs(window->cellz[c] * wincon->Ds[cs]) * dz;
      c = zm1;
      zm1 = window->zm1[c];
    }
    lasymp = 0.3 * sumzq / sumq;

    /*---------------------------------------------------------------*/
    /* Get the depths of the top and bottom of the pycnocline */
    mld(window, windat, wincon, c3, cb, &top_m, &bot_m, &kts, &ktb);

    /*---------------------------------------------------------------*/
    /* Compute cell centred viscosity, Vz, and diffusion, Kz. Use */
    /* the minimum of surface,bottom and q-moment length scales.  */
    c = c3;
    zm1 = window->zm1[c];

    while (c < cb) {
      cg = c;
      /* SIGMA : Multiply by depth */
      dz = (window->cellz[c] - window->cellz[zm1]) * wincon->Ds[cs];

      /* Vertical u1 gradient */
      du1 = vu1[c];

      /* Vertical u2 gradient */
      du2 = vu2[c];

      /* Vertical density gradient */
      drho = windat->dens_0[c] - windat->dens_0[zm1];

      du1dz = du1 / dz;
      du2dz = du2 / dz;
      drhodz = (drho / dz) / rho0;
      GH = -wincon->g * drhodz;
      GM = du1dz * du1dz + du2dz * du2dz;

      if (GM > 1.e-10) {
        Rgrad = GH / GM;
        Rflux =
          0.725 * (Rgrad + 0.186 -
                   sqrt(Rgrad * Rgrad - 0.316 * Rgrad + 0.0346));

        lscale =
          get_Lscale(window, windat, wincon, Rflux, GH, lasymp,
                     wincon->q[c], window->cellz[c], window->cellz[c3],
                     window->cellz[cb], top_m, bot_m, kts, ktb,
                     wincon->z0[cs]);
        if (windat->L)
          windat->L[c] = lscale;

        if (Rflux < Rcrit) {
          SH = 3 * A2 * (G1 - (G1 + G2) * Rflux) / (1 - Rflux);
          SM = SH * A1 * (BGC1 - (BGC1 + (6 * A1 + 3 * A2)) * Rflux) /
            (A2 * (B1 * G1 - (B1 * (G1 + G2) - 3 * A1) * Rflux));
          q3 =
            B1 * lscale * (windat->Vz[c] * GM +
                           windat->Kz[c] * wincon->g * drhodz);
          if (q3 > 0)
            wincon->q[c] = pow(q3, 1.0 / 3.0);
          else {
            hd_quit
              ("Vz: q3 must be greater than zero at %d (%d %d %d) : %f %f %f %f\n",
               c, geom->s2i[c], geom->s2j[c], geom->s2k[c], windat->Vz[c],
               windat->Kz[c], lscale, drhodz);
          }
          q3 = lscale * wincon->q[c];
          windat->Vz[c] = max(wincon->vz0, q3 * SM);
          windat->Kz[c] = max(wincon->kz0, q3 * SH);
        } else {
          windat->Vz[c] = wincon->vz0;
          windat->Kz[c] = wincon->kz0;
        }
      } else {
        windat->Vz[c] = wincon->vz0;
        windat->Kz[c] = wincon->kz0;
      }
      c = zm1;
      zm1 = window->zm1[c];
    }
  }
}

/* END closure_MY2_est()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the length scale for a water column which is       */
/* characterised by 3 regions; surface and bottom mixed layers and   */
/* a stably stratified interior layer. Follows the formulation of    */
/* Eifler and Schrimpf (1992) and Demirov et al (1998).              */
/*-------------------------------------------------------------------*/
double get_Lscale(geometry_t *window, /* Processing window */
                  window_t *windat, /* Window data structure */
                  win_priv_t *wincon, /* Window geometry / constants */
                  double Ri,    /* Flux Richardson number */
                  double N2,    /* Brunt-Vaisala frequency squared */
                  double Lo,    /* q length scale */
                  double q,     /* Turbulent kinetic energy */
                  double z,     /* Depth in the water column */
                  double sl,    /* Top of the water column */
                  double bl,    /* Bottom of the water column */
                  double top,   /* Top of the stratified layer */
                  double bot,   /* Bottom of the stratified layer */
                  int kts,      /* k index of the top of stratified layer */
                  int ktb,      /* k index of the bottom of stratified
                                   layer */
                  double z0     /* Bottom roughness length */
  )
{
  double e = 3.0;               /* Calibratable parameter */
  double c2 = 0.065;            /* Calibratable parameter */
  double Lmin = 0.01;           /* Minimum length scale in stratified area 
                                 */
  double N = 0.0;               /* Brunt-Vaisala frequency */
  double g = wincon->g;         /* Acceleration due to gravity */
  double zz1;                   /* Distance from the surface */
  double zz2;                   /* Distance from the bottom */
  double L, l1, l2;             /* Turbulence length scale */
  double drhodz;                /* Vertical density gradient */
  int zm1;                      /* Sparse coodinate at k-1 */
  int lof = 0;                  /* Flag to include q length scale (this */
  /* decreases L by approximately half).  */

  e = wincon->eparam;
  Lmin = wincon->Lmin;

  /* If Ri<0 (convective instability) do not correct the length */
  /* scale for stratification.  */
  if (Ri < 0.0)
    Ri = 0.0;

  /*-----------------------------------------------------------------*/
  /* One layer system. Assumes boundary layer overlap is one layer */
  if (top == sl && bot == bl) {
    double lt = sl - bl;        /* Thickness of the boundary layer */
    zz1 = sl - z;
    zz2 = z - bl;
    l1 =
      KAPPA * zz1 * pow(1.0 - Ri,
                        e) / (1.0 + KAPPA * zz1 / (c2 * lt)) + wincon->zs;
    l2 =
      KAPPA * zz2 * pow(1.0 - Ri,
                        e) / (1.0 + KAPPA * zz2 / (c2 * lt)) + z0;
    L = min(l1, l2);
    if (lof)
      L = L / (1 + L / Lo);
  }
  /*-----------------------------------------------------------------*/
  /* Three layer system */
  else {
    /* Unstably stratified, set L to the upper limit (Galperin et */
    /* al; 1988, eqn 26 and 28).  */
    if (N2 < 0.0) {
      L = q * sqrt(-0.0233 / N2);
    } else {
      /* L in the surface layer */
      if (z <= sl && z > top) {
        double lt = sl - top;
        zz1 = sl - z;
        L =
          KAPPA * zz1 * pow(1.0 - Ri,
                            e) / (1.0 + KAPPA * zz1 / (c2 * lt)) +
          wincon->zs;
        if (lof)
          L = L / (1 + L / Lo);
      }
      /* L in the bottom layer */
      else if (z < bot && z >= bl) {
        double lt = bot - bl;
        zz2 = z - bl;
        L =
          KAPPA * zz2 * pow(1.0 - Ri,
                            e) / (1.0 + KAPPA * zz2 / (c2 * lt)) + z0;
        if (lof)
          L = L / (1 + L / Lo);
      }
      /* L in the stratified region */
      else {
        double ci;
        double zz;
        double zos;
        double Ni;
        double qi;
        double dz;
        double mid = 0.5 * (top + bot);
        if (z > mid) {
          zz = sl - top;
          zos = wincon->zs;
          qi = wincon->q[kts];
          zm1 = window->zm1[kts];
          dz = window->cellz[kts] - window->cellz[zm1];
          drhodz = (windat->dens_0[kts] - windat->dens_0[zm1]) / dz;
          Ni = sqrt(-g * drhodz / rho0);
        } else if (z < mid) {
          zz = bot - bl;
          zos = z0;
          zm1 = window->zm1[ktb];
          qi = wincon->q[ktb];
          dz = window->cellz[ktb] - window->cellz[zm1];
          drhodz = (windat->dens_0[ktb] - windat->dens_0[zm1]) / dz;
          Ni = sqrt(-g * drhodz / rho0);
        } else {
          /* Takes into account the 2 layer system here : top == bot */
          zz = 0.5 * (sl - top + bot - bl);
          zos = 0.5 * (z0 + wincon->zs);
          zm1 = window->zm1[kts];
          qi = 0.5 * (wincon->q[kts] + wincon->q[ktb]);
          dz = window->cellz[kts] - window->cellz[zm1];
          drhodz = (windat->dens_0[kts] - windat->dens_0[zm1]) / dz;
          l1 = -g * drhodz / rho0;
          zm1 = window->zm1[ktb];
          dz = window->cellz[ktb] - window->cellz[zm1];
          drhodz = (windat->dens_0[ktb] - windat->dens_0[zm1]) / dz;
          /* Buoyancy frequency is the mean of that at the top and */
          /* bottom of the stratified layer.  */
          Ni = sqrt(0.5 * (-g * drhodz / rho0 + l1));
        }
        N = sqrt(N2);
        /* Note the relation q=sqrt(2k) holds; hence sqrt(k)=q/sqrt(2) */
        ci =
          Ni * (zos +
                KAPPA * zz * pow(1.0 - Ri,
                                 e) / (1.0 +
                                       KAPPA / c2)) / (qi / sqrt(2.0));
        /* Limit the maximum length scale according to Galperin et */
        /* al (1988), eqn 22.  */
        if (ci > 0.53)
          ci = 0.53;
        L = ci * q / (sqrt(2.0) * N);
        if (L < Lmin)
          L = Lmin;
      }
    }
  }

  /* Get the length scale based on the triangular profile and use */
  /* this as a maximum.  */
  l1 = KAPPA * (z - bl + z0);
  l2 = KAPPA * (sl - z + wincon->zs);
  l1 = min(l1 / (1 + l1 / Lo), l2 / (1 + l2 / Lo));
  if (L > l1)
    L = l1;
  return (L);
}

/* END Get_Lscale()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the top and bottom of the mixed layer based on    */
/* the turbulent kinetic energy. When the TKE (m2s-1) becomes less   */
/* than the threshold as distance is incremented from the top or     */
/* bottom, then the top or bottom of the mixed layer is              */
/* established. Note that the turbulence intensity, q=sqrt(2TKE).    */
/*-------------------------------------------------------------------*/
void mldk(geometry_t *window,   /* Processing window */
          window_t *windat,     /* Window data structure */
          win_priv_t *wincon,   /* Window geometry / constants */
          int cs,               /* z coordinate of surface */
          int cb,               /* z coordinate of bottom */
          double *tm,           /* Depth of the pycnocline surface */
          double *bm,           /* Depth of the pycnocline bottom */
          int *css,             /* k index of the pycnocline surface */
          int *cbb              /* k index of the pycnocline bottom */
  )
{
  int c;                        /* Sparse coordinate */
  double tke;                   /* Turbulent kinetic energy */
  double thr = 1e-5;            /* Density gradient threshold */

  /* Initialise the pycnocline depths */
  *tm = window->cellz[cs];
  *bm = window->cellz[cb];
  *css = cs;
  *cbb = cb;

  /*-----------------------------------------------------------------*/
  /* Get the depth of the top of the pycnocline */
  c = cs;
  while (c < cb) {
    tke = (wincon->q) ? 0.5 * wincon->q[c] * wincon->q[c] : windat->tke[c];
    if (tke < thr) {
      *tm = window->cellz[c];
      *css = c;
      break;
    }
    c = window->zm1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get the depth of the bottom of the pycnocline */
  c = cb;
  while (c != cs) {
    c = window->zp1[c];
    tke = (wincon->q) ? 0.5 * wincon->q[c] * wincon->q[c] : windat->tke[c];
    if (tke < thr) {
      *bm = window->cellz[c];
      *cbb = c;
      break;
    }
  }
}

/* END mldk()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the top and bottom of the mixed layer based on    */
/* the gradient of density. When the density gradient becomes        */
/* less than the threshold (density gradients are negative) as       */
/* distance is incremented from the top or bottom, then the top or   */
/* bottom of the mixed layer is established.                         */
/*-------------------------------------------------------------------*/
void mld(geometry_t *window,    /* Processing window */
         window_t *windat,      /* Window data structure */
         win_priv_t *wincon,    /* Window geometry / constants */
         int cs,                /* z coordinate of surface */
         int cb,                /* z coordinate of bottom */
         double *tm,            /* Depth of the pycnocline surface */
         double *bm,            /* Depth of the pycnocline bottom */
         int *css,              /* k index of the pycnocline surface */
         int *cbb               /* k index of the pycnocline bottom */
  )
{
  int c, zm1;                   /* Sparse coordinate */
  double dz;                    /* Layer thickness */
  double drhodz;                /* Density gradient */
  double thr = -0.01;           /* Density gradient threshold */
  int n, ns = 100;              /* Iteration safety number */
  /* -0.01 => generous, -0.003 conservative */

  /* Initialise the pycnocline depths */
  *tm = window->cellz[cs];
  *bm = window->cellz[cb];
  *css = cs;
  *cbb = cb;

  /*-----------------------------------------------------------------*/
  /* Get the depth of the top of the pycnocline */
  n = 0;
  c = cs;
  zm1 = window->zm1[c];
  while (c < cb && c != zm1) {
    dz = window->cellz[c] - window->cellz[zm1];
    drhodz = windat->dens_0[c] - windat->dens_0[zm1];
    drhodz = (drhodz / dz);
    if (drhodz < thr) {
      *tm = window->cellz[c];
      *css = c;
      break;
    }
    c = zm1;
    zm1 = window->zm1[c];
    n++;
    if (n > ns) {
      hd_warn("mld: Can't find top of pycnocline at (%d %d)\n", window->s2i[cs],window->s2j[cs]);
      break;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the depth of the bottom of the pycnocline */
  n = 0;
  c = cb;
  while (c != cs && c != window->zp1[c]) {
    zm1 = c;
    c = window->zp1[c];
    dz = window->cellz[c] - window->cellz[zm1];
    drhodz = windat->dens_0[c] - windat->dens_0[zm1];
    drhodz = (drhodz / dz);
    if (drhodz < thr) {
      *bm = window->cellz[c];
      *cbb = c;
      break;
    }
    n++;
    if (n > ns) {
      hd_warn("mld: Can't find bottom of pycnocline at (%d %d)\n", window->s2i[cs],window->s2j[cs]);
      break;
    }
  }
}

/* END mld()                                                         */
/*-------------------------------------------------------------------*/
