/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/closure/MY2-5.c
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
 *  $Id: MY2-5_HAR.c 5898 2018-08-23 02:07:21Z her127 $
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
#define C1S (C1)
#define C2  (0.7)
#define C2S (C2)
#define C3 (0.2)
#define CM0 (0.5562)
#define SEF (1.0)
#define GHMAX  (0.029)
#define GVMAX  (0.024)
#define GHMIN  (-0.28)
#define SQ2    (0.41)
#define SQ2L   (0.41)
#define ZRS    (0.1)

#define WAGE     (0.0)
#define WAVEFAC  (0.0)
#define ZRB      (0.1)
#define KMMIN    (0.0)
#define KHMIN    (0.0)
#define CKBG     (5.0e-3)
#define KBGCON   (5.0e-3)
#define EPS      (1.0e-30)
#define SMALL    (1.0e-10)
#define DOCEAN   (5.0e3)

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
#define E2 (1.0)
#define E3 (5.0)
#define E4 (1.33)
#define E5 (0.04)
#define E6 (5.0)
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

#define DS0   (A1*(1.0-3.0*C1S-6.0*A1/B1))
#define DM0   (A1*(1.0-3.0*C1-6.0*A1/B1))
#define DH0   (A2*(1.0-6.0*A1/B1))

/*-------------------------------------------------------------------*/
/* Routine to initialize the Mellor-Yamada 2.0 mixing scheme         */
/*-------------------------------------------------------------------*/
void closure_MY2_5_HAR_init(geometry_t *geom, /* Sparse global geometry
						 structure */
			    master_t *master  /* Master data structure */
  )
{
  int c, c2, cc;                /* Sparse coordinates */
  double l1, l2;                /* Distances to surface / bottom */

  /* Sanity check - make sure appropriate variables are present */
  if (!master->Q2 || !master->Q2L || !master->L || !master->Kq)
    hd_quit
      ("MY2.5 Harcourt mixing scheme requires tracers called tki, tki_l, lscale and Kq\n");

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

/* END closure_MY2_5_HAR_init()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the Mellor-Yamada 2.5 closure scheme.                             */
/*-------------------------------------------------------------------*/
void closure_MY2_5_HAR(geometry_t *window,  /* Processing window */
		       window_t *windat,  /* Window data structure */
		       win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c2, k;             /* Sparse coodinate / counter */
  int cs, ks;                   /* Surface sparse coordinate */
  int cb, kb;                   /* Bottom sparse coordinate */
  int zm1, zm2;                 /* Cell below cell c */
  int winsize = window->szc;    /* 3D work array size */
  int nzsize = window->nz + 1;  /* 1D work array size */
  double SH;                    /* Stability function for Kz */
  double SM;                    /* Stability function for Vz */
  double LQ;                    /* L * Q */
  double GM;                    /* Shear number */
  double GH;                    /* Richardson number for stability
                                   functions */
  double GS;                    /* Stokes forcing function */
  double GV;                    /* Stokes vortex forcing function */
  double P;                     /* Shear production */
  double B;                     /* Buoyancy production */
  double SVP;                   /* Stokes vortex production */
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
  double z0t;                   /* Surface roughness */
  double zrsw;                  /* Wave dependent surface roughness */
  double znrm;
  double zscl;
  double kbg;
  double dbmax, bfmax;
  double DS1, DS2, DM1, DM2, DM3, DH1, DH2, DH3, SS;
  double d1;                    /* Dummy */

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
  double *VP2 = wincon->w6;     /* Stokes vortex production */
  double *SP2 = wincon->w7;     /* Stokes shear production */
  double *rho = windat->dens_0;
  double *Kz = windat->Kz;
  double *Vz = windat->Vz;
  double *L = windat->L;
  double stkx[window->nz+1];
  double stky[window->nz+1];
  double fzs[window->nz+1];
  double *wu = wincon->d1;
  double *wv = wincon->d2;

  /*-----------------------------------------------------------------*/
  /* Get the east and north cell centered wind stress                */
  vel_cen(window, windat, wincon, windat->wind1, NULL, wu, wv, NULL, NULL, 1);

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(N2, 0, winsize * sizeof(double));
  memset(SP2, 0, winsize * sizeof(double));
  memset(VP2, 0, winsize * sizeof(double));

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
    /* Map the sparse arrays to local coordinates.                   */
    c = cs;
    for (k = ks; k >= kb; k--) {

      zm1 = window->zm1[c];
      dzface[k] = (k == kb) ? 0.5 * dz[k] : 0.5 * (dz[k - 1] + dz[k]);
      Kq[k] = windat->Kq[c];
      Kq[k] = windat->Kz[c];
      Q2[k] = windat->Q2[c];
      Q2L[k] = windat->Q2L[c];
      /* Stokes data                                                 */
      if (wincon->waves & STOKES_MIX) {
        if (wincon->waves & SPECTRAL) {
	  stkx[k] = windat->wave_stke1[c];
	  stky[k] = windat->wave_stke2[c];
	} else {
	  double depth = windat->eta[c2] - window->gridz[c];
	  double w = (windat->wave_period[c2]) ? 2.0 * PI / windat->wave_period[c2] : 0.0;
	  double kn = w * w / wincon->g;
	  double attn = exp(-2.0 * kn * depth);
	  stkx[k] = attn * windat->wave_ste1[c2];
	  stky[k] = attn * windat->wave_ste2[c2];
	}
      }
      c = zm1;
    }
    stkx[kb] = stkx[kb + 1];
    stky[kb] = stky[kb + 1];
    memset(Q2_Splus, 0, nzsize * sizeof(double));
    memset(Q2_Sminus, 0, nzsize * sizeof(double));
    memset(Q2L_Splus, 0, nzsize * sizeof(double));
    memset(Q2L_Sminus, 0, nzsize * sizeof(double));

    /* Surface roughness                                             */
    z0t = wincon->zs;
    /* z0t = max(wincon->zs, 1.6 * wincon->wave_hf * sqrt(windat->wave_amp[cs])); */
    /* z0t = wincon->wave_hf * windat->wave_amp[cs]; */
    zrsw = z0t;

    /*---------------------------------------------------------------*/
    /* Loop over the layer interfaces, calculating source/sink terms */
    /* in the Q2 and Q2l equations.                                  */
    znrm = zscl = 0.0;
    c = cb;
    do {
      double drho;
      double du1, du2;
      double su1, su2;

      zm1 = c;
      c = window->zp1[c];
      k = window->s2k[c];

      /* Vertical u1 gradient                                        */
      du1 = windat->u[c] - windat->u[zm1];

      /* Vertical u2 gradient                                        */
      du2 = windat->v[c] - windat->v[zm1];

      /* Density gradient                                            */
      drho = rho[c] - rho[zm1];

      /* Brunt Vaisala frequency (squared). Also used later to       */
      /* calculate stability functions.                              */
      N2[c] = -(wincon->g / RHO_0) * (drho / dzface[k]);

      /* Shear frequency squared.                                    */
      M2[c] = (du1 * du1 + du2 * du2) / (dzface[k] * dzface[k]);

      /* Stokes Vortex production (VPROD in SMCH) */
      if (wincon->waves & STOKES_MIX) {
        su1 = stkx[k] - stkx[k-1];
	su2 = stky[k] - stky[k-1];
	SP2[c] = (su1 * su1 + su2 * su2) / (dzface[k] * dzface[k]);
	VP2[c] = (su1 * du1 + su2 * du2) / (dzface[k] * dzface[k]);
      }

      /* Shear production - Rodi (1984) Eqn 2.48, or Burchard et al  */
      /* (1998) Eqn 9, with additional correction for internal wave  */
      /* shear as discussed by Mellor (1989), p267.                  */
      P = Vz[c] * (M2[c] + 0.7 * N2[c]);
      if (windat->s_prod) windat->s_prod[c] = P;

      /* Buoyancy production - refs as above.                        */
      B = -Kz[c] * N2[c];
      if (windat->b_prod) windat->b_prod[c] = B;

      /* Stokes-gradient momentum flux shear production */
      P += windat->Kq[c] * VP2[c];
      SVP = Vz[c] * VP2[c] + windat->Kq[c] * SP2[c];
      znrm += max(0.0, SVP);
      zscl += windat->L[c] * max(0.0, SVP);

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
      if ((P + B + SVP) > 0) {
        plus = P + B + SVP;
      } else {
        plus = P + SVP;
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
      /*minus =  diss * L[c] * (1.0 + E2 * (L[c] / l1) * (L[c] / l1));*/

      /* Waves                                                       */
      d1 = L[c] * (DOCEAN + ZRB + zrsw) / (KAPPA * (zrsw - window->gridz[c2]) * 
					   (DOCEAN + window->gridz[c] + ZRB));
      minus = E2 * diss * (1.0 + E4 * (d1 * d1));


      /* Q2L source/sink, splitting positive and negative parts      */
      /* (Patankar 1980).                                            */
      if ((P * E1 + B * E3 + E6 * SVP) > 0) {
        plus = (P * E1 + B * E3 + E6 * SVP) * L[c];
      } else {
	plus = (P * E1 + E6 * SVP) * L[c];
        minus -= B * E3 * L[c];
      }
      /*
      plus = E1 * P + E3 * B + E6 * SVP;
      */

      /* Distribute Q2L sources/sinks to adjacent layers             */
      Q2L_Splus[k - 1] += 0.5 * plus;
      Q2L_Splus[k] += 0.5 * plus;
      Q2L_Sminus[k - 1] += 0.5 * minus;
      Q2L_Sminus[k] += 0.5 * minus;

    } while (c != cs);

    zscl = (znrm) ? zscl / znrm : EPS;
    for (k = ks; k >= kb; k--) {
      double depth = windat->eta[c2] - window->gridz[c];
      fzs[k] = 1.0 + tanh(-0.25 * depth / zscl);
      c = window->zm1[c];
    } 
    c = cs;

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
    v1 = wu[c2];
    v2 = wv[c2] ;
    ustr = sqrt(v1 * v1 + v2 * v2) / RHO_0;
    /* Drichlet condition for wave breaking source                   */
    zm1 = window->zm1[cs];
    zm2 = window->zm1[zm1];
    d1 = 2.0 * windat->Q2[zm1]-windat->Q2[zm2];
    windat->Q2[cs] = max(max(d1, CONST1 * ustr), windat->Q2[cs]);

    /* Bottom                                                        */
    v1 = windat->u[cb];
    v2 = windat->v[cb];
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
    if (wincon->waves & VERTMIX) l1 = windat->wave_amp[c2];
    L[cs] =
      KAPPA * (0.5 *(windat->eta[c2] - window->gridz[cs] * wincon->Ds[c2]) + l1);
    /*L[cs] = 0.5 * (windat->Q2L[zm1] + windat->Q2L[c]) / 
      max(SMALL, 0.5 * (windat->Q2[zm1] + windat->Q2[c]));*/
    windat->Q2L[cs] = windat->Q2[cs] * L[cs];
    windat->Q2L[cs] = windat->Q2[cs] * zrsw;

    /* Limit Q2L at the surface                                      */
    windat->Q2L[cs] = max(zrsw * windat->Q2[cs], windat->Q2L[cs]);

    /* Bottom                                                        */
    L[cb] =
      KAPPA * (window->cellz[cb] * wincon->Ds[c2] - window->botz[c2] +
               wincon->z0[c2]);
    windat->Q2L[cb] = windat->Q2[cb] * L[cb];

    /*---------------------------------------------------------------*/
    /* Get the mixing coefficients. Note: SM and SH limit to         */
    /* infinity when GH approaches 0.0288 (Mellor & Yamada, Fig 3)   */
    c = window->zm1[cs];
    dbmax = L[cs];
    bfmax = 0.0;
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
      /*L[c] = q2l / max(q2, SMALL);*/
      if (N2[c] > 0.0) {
	double dbuoy = sqrt(q2 / N2[c]);
	double depth = windat->eta[c2] - window->gridz[c];
	if (L[c] < max(dbuoy, depth))
	  dbmax = max(dbmax, L[c]);
	if (dbuoy >= dbmax)
	  bfmax = max(N2[c], bfmax);
      }

      /*
      if (wincon->s_func == NULL) {
      */
	/* Get the Richardson number using the limits                */
	GH = -N2[c] * L[c] * L[c] / q2;
	GM = L[c] * L[c] * M2[c] / q2;
	GS = L[c] * L[c] * SP2[c] / q2;
	GV = L[c] * L[c] * VP2[c] / q2;

	d1 = (1.0 - fzs[window->s2k[c]]);
	GS *= (d1 * d1);
	GV *= d1;
	GH = max(GH, GHMIN);
	GV = max(GV, GHMIN);
	GH = min(GH, GHMAX);
	GV = min(GV, GVMAX);
	GS = min(GS, GVMAX);

	/* Get the stability functions (Galperin et al, 1988, eqn 24 */
	/* & 25. Note these are derived assuming P+B/diss=1.         */
	DS1 = 1.0 - 9.0 * A1 * (A2 * GH + A1 * GV);
	DS2 = -9.0 * A1 * A2 * C2S * GH;
	DM1 = 1.0 - 9.0 * A1 * (A2 * GH + 4.0 * A1 * GV);
	DM2 = 9.0 * A1 * (2.0 * A1 + A2 * (1.0 - C2)) * GH;
	DM3 = 27.0 * A1 * A1 * GS;
	DH1 = 1.0 - 3.0 * A2 * ((6.0 * A1 + B2 * (1.0 - C3)) * GH +
				3.0 * A2 * (1.0 - C2) * GV - 3.0 * A2 * C2S * GS);
	DH2 = 9.0 * A2 * GV * (2.0 * A1 + A2);
	DH3 = 9.0 * A2 * GS * (2.0 * A1 + A2);
	
	/* Stability functions.                                      */
	SH = (DH0 * DM1 * DS1 + DH2 * DM0 * DS1 + (DH3 * DM1 + DM3 * DH2) * DS0 ) /
	    (DH1 * DM1 * DS1 - DH2 * DM2 * DS1 - (DH3 * DM1 + DM3 * DH2) * DS2);
	SH = max(SH, 0.0);

	SS = (DS0 + DS2 * SH) / DS1;
	SS = max(SS, 0.0);

	SM = (DM0 + DM2 * SH + DM3 * SS) / DM1;
	SM = max(SM, 0.0);

	/* Get the mixing coefficients */
	LQ = L[c] * sqrt(q2);
	windat->Kq[c] = min(MAXVZ, (LQ * SS * (1.0 - fzs[window->s2k[c]]) + windat->Kq[c]) * 0.5);
	/*
      } else {
	*/
        /* Alternative stability functions                           */
	/*
        wincon->s_func(-2.0 * GH, 2.0 * GM, &SM, &SH);
	*/
	/* Stability functions accept aN and aM as arguments, and    */
	/* return cmu and cmu_d rather rthan SH and SM.              */
	/*
	SM /= sqrt(2.0);
	SH /= sqrt(2.0);
	LQ = L[c] * sqrt(q2);
	windat->Kq[c] = min(MAXVZ, (LQ * SQ + windat->Kq[c]) * 0.5);
      }
	*/
      /* Get the mixing coefficients */
      Vz[c] = min(MAXVZ, max((LQ * SM + Vz[c]) * 0.5, wincon->vz0));
      Kz[c] = min(MAXKZ, max((LQ * SH + Kz[c]) * 0.5, wincon->kz0));

      if (N2[c] / M2[c] >= 0.0) {
	d1 = N2[c] / (0.7 * M2[c]);
	kbg = max(0.0, 1.0 - d1 * d1);
	kbg = CKBG *kbg * kbg * kbg;
      } else
        kbg = KBGCON;

      if (N2[c] > bfmax) {
        Vz[c] = max(Vz[c], wincon->vz0 + kbg);
	Kz[c] = max(Vz[c], wincon->kz0 + kbg);
	windat->Kq[c] = max(windat->Kq[c], wincon->kz0 + kbg);
      }

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

/* END closure_MY2_5_HAR()                                           */
/*-------------------------------------------------------------------*/
