/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/closure/MY2-5.c
 *  
 *  Description:
 *  Implements Mellor-Yamada level 2.5 mixing scheme.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: MY2-5.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Constants                                                        */
#define KAPPA (0.4)             /* Von Karman's constant */
#define RHO_0 (1025.0)          /* Mean density */
#define A1 (0.92)
#define B1 (16.6)
#define A2 (0.74)
#define B2 (10.1)
#define C1 (0.08)
#define CM0 (0.5562)

/* Minimum turbulent kinetic energy - GOTM code                      */
#define MIN_TKE          (1e-8)
/* Minimum dissipation - from GOTM code                              */
#define MIN_DISS         (1e-12)

/* Minimum turbulent kinetic energy - Burchard et al (1998)          */
/*#define MIN_TKE          (7.6e-6)*/
/* Minimum dissipation - from CRS code                               */
/*#define MIN_DISS         (5e-10)*/

/* Minimum turbulent kinetic energy - Warner et al (2005)            */
/*#define MIN_TKE          (5e-6)*/
/* Minimum dissipation - Warner et al (2005) using kL=1e-8           */
/*#define MIN_DISS         (9.62e-7)*/

/* Maximum possible viscosity/diffusivity (for stability)            */
#define MAXVZ	 	 (1.0)
#define MAXKZ		 (1.0)
/* Stability function for turbulent diffusion; Mellor and Yamada     */
/* (1982), p862. Note sqrt(2).SQ=ck in  Burchard et al (1998).       */
#define SQ (0.2)

#define E1 (1.8)
#define E2 (1.33)
#define E3 (1.8)
#define CI (0.53)
#define r1 (3.0)
#define r2 (6.0)
#define r3 (9.0)
#define r4 (18.0)
#define r7 (2.0/3.0)
#define COEF1 (A2*(1.0-r2*A1/B1))
#define COEF2 (r1*A2*B2+r4*A1*A2)
#define COEF3 (A1*(1.0-r1*C1-r2*A1/B1))
#define COEF4 (r4*A1*A1+r3*A1*A2)
#define COEF5 (r3*A1*A2)


/*-------------------------------------------------------------------*/
/* Routine to initialize the Mellor-Yamada 2.0 mixing scheme         */
/*-------------------------------------------------------------------*/
void closure_MY2_5_init(geometry_t *geom, /* Sparse global geometry
                                             structure */
                        master_t *master  /* Master data structure */
  )
{
  int c, c2, cc;                /* Sparse coordinates */
  double l1, l2;                /* Distances to surface / bottom */

  /* Sanity check - make sure appropriate variables are present */
  if (!master->Q2 || !master->Q2L || !master->L || !master->Kq)
    hd_quit
      ("MY2.5 mixing scheme requires tracers called tki, tki_l, lscale and Kq\n");

  /* Set the minimum length scales */
  if (!master->min_tke)
    master->min_tke = MIN_TKE;
  if (!master->min_diss)
    master->min_diss = MIN_DISS;
  /* Minimum length scale from Burchard et al (1998) Eqn 11. */
  /* This approach is used in GOTM */
  if (!master->Lmin)
    master->Lmin = CM0 * CM0 * CM0 * pow(master->min_tke, 1.5) / 
      master->min_diss;

  /* In the MY2.5 mixing model, Q2 and Q2L are vertically mixed */
  /* using the eddy viscosity. For this version, diffusion is */
  /* handled in the code below, so tracer diffusion is turned off in */
  /* master_build() for these variables.  */

  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    master->Vz[c] = master->vz0;
    master->Kz[c] = master->kz0;
    master->Q2[c] = 2.0 * master->min_tke;
    /* Initialise L with a triangular function */
    l1 = KAPPA * ((geom->cellz[c] - geom->botz[c2]) * master->Ds[c2] +
                  master->z0[c2]);
    l2 = KAPPA * ((geom->cellz[c2] - geom->cellz[c]) * master->Ds[c2] +
                  master->zs);
    master->L[c] = min(l1, l2);
    master->Q2L[c] = master->Q2[c] * master->L[c];
    master->Kq[c] = SQ * sqrt(master->Q2[c]) * master->L[c];
  }
}

/* END closure_MY2_5_init()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the Mellor-Yamada 2.5 closure scheme.                             */
/*-------------------------------------------------------------------*/
void closure_MY2_5(geometry_t *window,  /* Processing window */
                   window_t *windat,  /* Window data structure */
                   win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c2, k;             /* Sparse coodinate / counter */
  int cs, ks;                   /* Surface sparse coordinate */
  int cb, kb;                   /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */
  int xp1, yp1;                 /* Cell to the west, north of cell c */
  int winsize = window->sgsiz;  /* 3D work array size */
  int nzsize = window->nz + 1;  /* 1D work array size */
  double SH;                    /* Stability function for Kz */
  double SM;                    /* Stability function for Vz */
  double LQ;                    /* L * Q */
  double GH;                    /* Richardson number for stability
                                   functions */
  double P;                     /* Shear production */
  double B;                     /* Buoyancy production */
  double CONST1;                /* Constant for boundary conditions */
  double plus, minus;           /* Source and sink terms */
  double diss;                  /* Dissipation */
  double dzdt;                  /* dz / dt */
  double v1, v2;                /* Stress components */
  double l1, l2;                /* Mixing lengths */
  double ustr;                  /* Surface or bottom friction velocity */
  double q2, q2l;               /* Face centered Q2 and Q2L */
  double mq2;                   /* Minimum Q2 */
  double mq2l;                  /* Minimum Q2L */
  double GM;                    /* Shear number */

  /* Set pointers.  */
  /* Note: the 3D work arrays wincon->w# could be used for the */
  /* dummy arrays below, but execution speed is considerably faster */
  /* when work array access is sequential in memory, hence the */
  /* mapping to a contiguous vertical 1D work array wincon->v#. The */
  /* code below contains a mixture of variables accessed in the */
  /* sparse array, c, and local array, k, with mapping achieved via */
  /* k=window->s2k[c].  */
  int *ctp = window->nsur_t;    /* Surface cells to process */
  int vcs = window->b2_t;       /* Last index of surface cells */
  double *C = wincon->v1;       /* Main diagonal */
  double *Cm1 = wincon->v2;     /* Lower diagonal */
  double *Cp1 = wincon->v3;     /* Upper diagonal */
  double *rhs = wincon->v4;     /* Right hand side */
  double *dzface = wincon->v5;  /* Cell thickness for faces */
  double *Q2_Splus = wincon->v6;  /* Source term for Q2 */
  double *Q2_Sminus = wincon->v7; /* Sink term for Q2 */
  double *Q2L_Splus = wincon->v8; /* Source term for Q2L */
  double *Q2L_Sminus = wincon->v9;  /* Sink term for Q2L */
  double *dz = wincon->v10;     /* Cell thickness in local coords */
  double *N2 = wincon->w1;      /* Buoyancy frequency */
  double *Kq = wincon->w2;      /* Mixing diffusivity, local coord */
  double *Q2 = wincon->w3;      /* Q2 in local coordinates */
  double *Q2L = wincon->w4;     /* Q2L in local coordinates */
  double *M2 = wincon->w5;      /* Shear frequency */
  double *rho = windat->dens_0;
  double *Kz = windat->Kz;
  double *Vz = windat->Vz;
  double *L = windat->L;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(N2, 0, winsize * sizeof(double));

  /* Multiplier for boundary conditions; eqn. 54 Mellor&Yamada(1982) */
  CONST1 = pow(B1, r7);
  /* Minimum turbulent intensity : k=q^2/2                           */
  mq2 = 2.0 * wincon->min_tke;
  mq2l = mq2 * wincon->Lmin;

  /*-----------------------------------------------------------------*/
  /* Loop therough the surface cells in this window */
  for (cc = 1; cc <= vcs; cc++) {

    cs = c = ctp[cc];           /* Set cs to the surface sparse coordinate 
                                 */
    c2 = window->m2d[c];        /* 2D sparse location corresponding to c */
    ks = window->s2k[cs];       /* Surface local coordinate */
    zm1 = window->zm1[c];

    /*---------------------------------------------------------------*/
    /* Single layer case (i.e. surface lies in the bottom layer) */
    if (zm1 == window->zm1[zm1]) {
      continue;
    }

    /*---------------------------------------------------------------*/
    /* Find the bottom sparse coordinate. The sparse coordinate, c, */
    /* lies in the sediment at the end of the loop.  */
    c = cs;
    while (c != zm1 && wincon->dz[c]) {
      dz[window->s2k[c]] = wincon->dz[c];
      c = zm1;
      zm1 = window->zm1[c];
    }
    cb = window->zp1[c];        /* Set cb to the bottom sparse coordinate */
    kb = window->s2k[cb];       /* Bottom local coordinate */

    /*---------------------------------------------------------------*/
    /* Map the sparse arrays to local coordinates.  */
    c = cs;
    for (k = ks; k >= kb; k--) {
      zm1 = window->zm1[c];
      dzface[k] = (k == kb) ? 0.5 * dz[k] : 0.5 * (dz[k - 1] + dz[k]);
      Kq[k] = windat->Kq[c];
      Q2[k] = windat->Q2[c];
      Q2L[k] = windat->Q2L[c];
      c = zm1;
    }
    memset(Q2_Splus, 0, nzsize * sizeof(double));
    memset(Q2_Sminus, 0, nzsize * sizeof(double));
    memset(Q2L_Splus, 0, nzsize * sizeof(double));
    memset(Q2L_Sminus, 0, nzsize * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Loop over the layer interfaces, calculating source/sink terms */
    /* in the Q2 and Q2l equations.                                  */
    c = cb;
    do {
      double drho;
      double du1, du2;

      zm1 = c;
      c = window->zp1[c];
      k = window->s2k[c];

      /* Vertical u1 gradient                                        */
      xp1 = window->xp1[c];
      du1 = 0.5 * (windat->u1[c] - windat->u1[zm1] +
                   windat->u1[xp1] - windat->u1[window->zm1[xp1]]);

      /* Vertical u2 gradient                                        */
      yp1 = window->yp1[c];
      du2 = 0.5 * (windat->u2[c] - windat->u2[zm1] +
                   windat->u2[yp1] - windat->u2[window->zm1[yp1]]);

      /* Density gradient                                            */
      drho = rho[c] - rho[zm1];

      /* Brunt Vaisala frequency (squared). Also used later to       */
      /* calculate stability functions.                              */
      N2[c] = -(wincon->g / RHO_0) * (drho / dzface[k]);

      /* Shear frequency squared.                                    */
      M2[c] = (du1 * du1 + du2 * du2) / (dzface[k] * dzface[k]);

      /* Shear production - Rodi (1984) Eqn 2.48, or Burchard et al  */
      /* (1998) Eqn 9, with additional correction for internal wave  */
      /* shear as discussed by Mellor (1989), p267.                  */
      P = Vz[c] * (M2[c] + 0.7 * N2[c]);
      if (windat->s_prod) windat->s_prod[c] = P;

      /* Buoyancy production - refs as above.                        */
      B = -Kz[c] * N2[c];
      if (windat->b_prod) windat->b_prod[c] = B;

      /* Interpolate Q2 and Q2L from layer centres to interfaces     */
      q2 = max(0.5 * (Q2[k - 1] + Q2[k]), mq2);
      q2l = max(0.5 * (Q2L[k - 1] + Q2L[k]), mq2l);

      /* Dissipation : loss term for Q2 (e=q^3/B1.L eqn 11 Burchard  */
      /* et al 1998 and eqn 25c & 12 Mellor & Yamada (1982)).        */
      /* Note : do not use a length scale as q2l / q2, since this    */
      /* will not have been limited according to                     */
      /* Galperin et al (1988).                                      */
      /* Note that loss terms are multiplied by q2(t+1)/q2(t) in the */
      /* implicit solution, e.g. Burchard et al (1998), p10547.      */
      diss = minus = max(q2 * sqrt(q2) / (B1 * L[c]), wincon->min_diss);

      /* Q2 source/sink, splitting positive and negative parts       */
      /* (Patankar 1980). Note : k=q^2/2 hence the Q2 equation is    */
      /* equal to the 2*k equation. Source and sink terms for Q2     */
      /* must therefore be multiplied by 2.0 (Mellor and Yamada,     */
      /* 1982 eqn 16). This cancels the 0.5 required to get the mean */
      /* between layers when sources/sinks are distributed to        */
      /* adjacent layers below.                                      */
      if ((P + B) > 0) {
        plus = P + B;
      } else {
        plus = P;
        minus -= B;
      }

      /* Distribute Q2 sources/sinks to adjacent layers              */
      Q2_Splus[k - 1] += plus;
      Q2_Splus[k] += plus;
      Q2_Sminus[k - 1] += minus;
      Q2_Sminus[k] += minus;

      /* Sources and sinks for Q2L.                                  */
      /* Loss term for Q2l : eqn 48 Mellor & Yamada (1982).          */
      l1 = windat->eta[c2] - window->gridz[c] * wincon->Ds[c2] + wincon->zs;
      l2 = window->gridz[c] * wincon->Ds[c2] - window->botz[c2] +
	wincon->z0[c2];
      /* Triangular prescribed length scale                          */
      l1 = KAPPA * min(l1, l2); 
      /* Parabolic prescribed length scale                           */
      /*l1 = KAPPA * l1 * l2 / (l1 + l2);*/
      minus =  diss * L[c] * (1.0 + E2 * (L[c] / l1) * (L[c] / l1));

      /* Q2L source/sink, splitting positive and negative parts      */
      /* (Patankar 1980).                                            */
      if ((P * E1 + B * E3) > 0) {
        plus = (P * E1 + B * E3) * L[c];
      } else {
        plus = P * E1 * L[c];
        minus -= B * E3 * L[c];
      }

      /* Distribute Q2L sources/sinks to adjacent layers             */
      Q2L_Splus[k - 1] += 0.5 * plus;
      Q2L_Splus[k] += 0.5 * plus;
      Q2L_Sminus[k - 1] += 0.5 * minus;
      Q2L_Sminus[k] += 0.5 * minus;

    } while (c != cs);

    /*---------------------------------------------------------------*/
    /* Get Q2 at the forward time step. This is solved using         */
    /* simplified Gaussian elimination (Thomas algorithm). First set */
    /* up the tri-diagonal matrix.                                   */
    /* Q2 bottom boundary condition - zero flux (Neumann condition)  */
    dzdt = dz[kb] / windat->dt;
    Cm1[kb] = 0.0;
    Cp1[kb] = -Kq[kb + 1] / dzface[kb + 1];
    C[kb] = dzdt - Cp1[kb] + dz[kb] * Q2_Sminus[kb] / max(Q2[kb], mq2);
    rhs[kb] = dzdt * Q2[kb] + dz[kb] * Q2_Splus[kb];

    /* Q2 surface boundary condition - zero flux condition           */
    dzdt = dz[ks] / windat->dt;
    Cm1[ks] = -Kq[ks] / dzface[ks];
    Cp1[ks] = 0.0;
    C[ks] = dzdt - Cm1[ks] + dz[ks] * Q2_Sminus[ks] / max(Q2[ks], mq2);
    rhs[ks] = dzdt * Q2[ks] + dz[ks] * Q2_Splus[ks];

    rhs[ks] = dzdt * Q2[ks] + dz[ks] * Q2_Splus[ks];

    /* Q2 mid-water column                                           */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / windat->dt;
      Cm1[k] = -Kq[k] / dzface[k];
      Cp1[k] = -Kq[k + 1] / dzface[k + 1];
      C[k] = dzdt - Cm1[k] - Cp1[k] + dz[k] * Q2_Sminus[k] / 
	max(Q2[k], mq2);
      rhs[k] = dzdt * Q2[k] + dz[k] * Q2_Splus[k];
    }

    /* Solve tridiagonal system for Q2                               */
    tridiagonal(window, cs, cb, C, Cm1, Cp1, rhs, windat->Q2, mq2);

    /* Set Q2 at the surface and bottom using a Dirichlet condition  */
    /* Mellor and Yamada (1982) eqn 54.                              */
    /* Surface. Note ustr=sqrt(sqrt(v1^2+v2^2)/rho) and              */
    /* Q2=CONST1*ustr*ustr.                                          */
    xp1 = window->xp1[c2];
    yp1 = window->yp1[c2];
    v1 = 0.5 * (windat->wind1[c] + windat->wind1[xp1]);
    v2 = 0.5 * (windat->wind2[c] + windat->wind2[yp1]);
    ustr = sqrt(v1 * v1 + v2 * v2) / RHO_0;
    windat->Q2[cs] = CONST1 * ustr;

    /* Bottom                                                        */
    xp1 = window->xp1[cb];
    yp1 = window->yp1[cb];
    v1 = 0.5 * (windat->u1[cb] + windat->u1[xp1]);
    v2 = 0.5 * (windat->u2[cb] + windat->u2[yp1]);
    ustr = wincon->Cd[c2] * (v1 * v1 + v2 * v2);
    windat->Q2[cb] = CONST1 * ustr;

    /*---------------------------------------------------------------*/
    /* Get Q2L at the forward time step.                             */
    /* Q2L bottom boundary condition - zero flux (Neumann condition) */
    dzdt = dz[kb] / windat->dt;
    C[kb] = dzdt - Cp1[kb] + dz[kb] * Q2L_Sminus[kb] / 
      max(Q2L[kb], mq2l);
    rhs[kb] = dzdt * Q2L[kb] + dz[kb] * Q2L_Splus[kb];

    /* Q2L surface boundary condition - zero flux condition          */
    dzdt = dz[ks] / windat->dt;
    C[ks] = dzdt - Cm1[ks] + dz[ks] * Q2L_Sminus[ks] / 
      max(Q2L[ks], mq2l);
    rhs[ks] = dzdt * Q2L[ks] + dz[ks] * Q2L_Splus[ks];

    /* Q2L mid-water column                                          */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / windat->dt;
      C[k] = dzdt - Cm1[k] - Cp1[k] + dz[k] * Q2L_Sminus[k] / 
	max(Q2L[k], mq2l);
      rhs[k] = dzdt * Q2L[k] + dz[k] * Q2L_Splus[k];
    }

    /* Solve tridiagonal system for Q2L                              */
    tridiagonal(window, cs, cb, C, Cm1, Cp1, rhs, windat->Q2L, mq2l);

    /* Set Q2L at the surface and bottom using a Dirichlet condition */
    /* Surface.                                                      */
    l1 = wincon->zs;
    if (wincon->waves & VERTMIX) l1 = windat->wave_amp[cs];
    L[cs] =
      KAPPA * (windat->eta[c2] - window->cellz[cs] * wincon->Ds[c2] + l1);
    windat->Q2L[cs] = windat->Q2[cs] * L[cs];

    /* Bottom                                                        */
    L[cb] =
      KAPPA * (window->cellz[cb] * wincon->Ds[c2] - window->botz[c2] +
               wincon->z0[c2]);
    windat->Q2L[cb] = windat->Q2[cb] * L[cb];

    /*---------------------------------------------------------------*/
    /* Get the mixing coefficients. Note: SM and SH limit to         */
    /* infinity when GH approaches 0.0288 (Mellor & Yamada, Fig 3)   */
    c = window->zm1[cs];
    while (c != cb) {
      zm1 = window->zm1[c];

      /* Diffusivity and viscosity correspond to layers faces, hence */
      /* Q2 and Q2L are interpolated from centers to faces.          */
      q2 = 0.5 * (windat->Q2[zm1] + windat->Q2[c]);
      q2l = 0.5 * (windat->Q2L[zm1] + windat->Q2L[c]);

      /* Limit the maximum length scale according to Galperin et al  */
      /* (1988), eqn 22.                                             */
      if (N2[c] > 0.0)
	L[c] = max(wincon->Lmin, min(q2l / q2, CI * sqrt(q2 / N2[c])));
      else
	L[c] = max(wincon->Lmin, q2l / q2);

      /* Get the Richardson number using the limits                  */
      GH = -N2[c] * L[c] * L[c] / q2;
      GM = L[c] * L[c] * M2[c] / q2;

      /* Get the stability functions (Galperin et al, 1988, eqn 24 & */
      /* 25. Note these are derived assuming P+B/diss=1.             */
      wincon->s_func(-2.0 * GH, 2.0 * GM, &SM, &SH);

      /* Stability functions accept aN and aM as arguments, and      */
      /* return cmu and cmu_d rather rthan SH and SM.                */
      SM /= sqrt(2.0);
      SH /= sqrt(2.0);

      /* Get the mixing coefficients */
      LQ = L[c] * sqrt(q2);
      windat->Kq[c] = min(MAXVZ, (LQ * SQ + windat->Kq[c]) * 0.5);
      Vz[c] = min(MAXVZ, max((LQ * SM + Vz[c]) * 0.5, wincon->vz0));
      Kz[c] = min(MAXKZ, max((LQ * SH + Kz[c]) * 0.5, wincon->kz0));
      c = zm1;
    }

    /* Set a no gradient condition at the surface and bottom layers  */
    Vz[cs] = Vz[window->zm1[cs]];
    Kz[cs] = Kz[window->zm1[cs]];
    windat->Kq[cs] = windat->Kq[window->zm1[cs]];
    Vz[cb] = Vz[window->zp1[cb]];
    Kz[cb] = Kz[window->zp1[cb]];
    windat->Kq[cb] = windat->Kq[window->zm1[cb]];
  }
}

/* END closure_MY2_5()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Solves a tridiagonal system of linear equations using simplified  */
/* Gaussian elimination (Thomas algorithm).                          */
/*-------------------------------------------------------------------*/
void tridiagonal(geometry_t *window,  /* Window geometry structure */
                 int cs, int cb,  /* Surface and bottom sparse coordinate */
                 double *C,     /* Main diagonal */
                 double *Cm1,   /* Lower diagonal */
                 double *Cp1,   /* Upper diagonal */
                 double *rhs,   /* Right hand side of equation */
                 double *var,   /* df_variable_t to solve */
                 double minvar  /* Minimum variable value */
  )
{
  win_priv_t *wincon = window->wincon;
  int c, k;                     /* Sparse coordinate at k */
  int ks = window->s2k[cs];     /* Surface local coordinate */
  int kb = window->s2k[cb];     /* Bottom local coordinate */
  double div;                   /* Dummy */
  double *sol = wincon->v11;    /* Work array */
  double *ud = wincon->v12;     /* Work array */

  /* Solve tridiagonal system for var */
  div = C[kb];
  sol[kb] = rhs[kb] / div;
  for (k = kb + 1; k <= ks; k++) {
    ud[k] = Cp1[k - 1] / div;
    div = C[k] - Cm1[k] * ud[k];
    if (div == 0.0)
	hd_quit_and_dump("tridiagonal : zero divisor @ (%d %d %d)\n", 
			 window->s2i[cs], window->s2j[cs], k);
    sol[k] = (rhs[k] - Cm1[k] * sol[k - 1]) / div;
  }
  for (k = ks - 1; k >= kb; k--) {
    sol[k] -= ud[k + 1] * sol[k + 1];
  }
  c = cb;
  for (k = kb; k <= ks; k++) {
    var[c] = max(sol[k], minvar);
    c = window->zp1[c];
  }
}

/* END tridiagonal()                                                 */
/*-------------------------------------------------------------------*/
