/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/closure/k-e-implicit.c
 *  
 *  Description:
 *  Implements K-epsilon mixing scheme. The code
 *  here is based primarily on equations described
 *  in the following references:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: k-e-implicit.c 5898 2018-08-23 02:07:21Z her127 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

double K_e(geometry_t *window, int cs, double wflux, double z0s);
double K_k(geometry_t *window, int cs, double wflux, double z0s);

/*-------------------------------------------------------------------*/
/* Constants                                                         */
#define VK		(VON_KAR)         /* Von Karman's constant */
#define KAPPA (0.4)             /* Von Karman's constant */

/* Minimum possible viscosity/diffusivity (akin to molecular diffusion) */
#define NU		(1.0e-5)

/* Maximum possible viscosity/diffusivity (for stability)            */
#define MAXVZ		(1.0)
#define MAXKZ		(1.0)

/* Burchard et al (1998) Table 1. Note that in CRS code, and also in */
/* Rodi (1984), the value of c_mu is given respectively as CM0**4,   */
/* and 0.09. These are almost equivalent.                            */
#define CM0		(0.5562)
#define CM0c            (CM0*CM0*CM0)

/* Burchard et al (1998) Table 1                                     */
#define C2E		(1.92)

/* Rodi (1984) p28 gives a formula for this, but here we use the     */
/* value from Burchard et al (1998) Table 1, which is almost         */
/* identical.                                                        */
#define C1E		(1.44)

/* Burchard et al (1998) Table 1                                     */
#define SIG_K    	(1.0)

/* Burchard et al (1998) gives an equation (22) for sigma_e, which   */
/* depends on the values above. The definition here would look like  */
/*     (VK*VK/((C2E-C1E)*CM0*CM0)                                    */
/* Here we just use their quoted  value from Table 1. Rodi (1984)    */
/* uses a different value  of 1.3.                                   */
#define SIG_E    	(1.08)

/* Minimum turbulent kinetic energy - Burchard et al (1998)          */
#define MIN_TKE          (7.6e-6)

/* Minimum dissipation - Warner et al (2005) Table 1 (also GOTM code)*/
/* define MIN_DISS          (1e-12)*/

/* Minimum dissipation - from CRS code                               */
#define MIN_DISS          (5e-10)

/* Mean density */
#define RHO_0		(1025.0)

/* Wave parameters */
#define a_e  (-17.78)  /* Umlauf and Burchard, 2003, Table 6 */
#define L_e  (0.025)   /* Umlauf and Burchard, 2003, Table 6 */

/*-------------------------------------------------------------------*/
/* Routine to initialize the implicit k-e mixing scheme              */
/*-------------------------------------------------------------------*/
void closure_k_e_implicit_init(geometry_t *geom,  /* Sparse global
                                                     geometry structure */
                               master_t *master /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coordinates */

  /* Sanity check - make sure appropriate variables are present */
  if (!master->tke || !master->diss)
    hd_quit
      ("k-epsilon mixing scheme requires tracers called tke and diss\n");

  /* In the k-e mixing model, TKE and dissipation are vertically */
  /* mixed using the eddy viscosity (divided by the Schmidt numbers */
  /* sigma_k and sigma_e respectively - see Burchard et al 1998 Eqns */
  /* 20 and 21). For this version, diffusion is handled in the code */
  /* below, so tracer diffusion is turned off in master_build() for */
  /* these variables.  */

  /* Read optional minimum values for tke and diffusion. If not */
  /* defined, set to default values defined above.  */
  if (!master->min_tke)
    master->min_tke = MIN_TKE;
  if (!master->min_diss)
    master->min_diss = MIN_DISS;
  if (!master->vz0)
    master->vz0 = NU;
  if (!master->kz0)
    master->kz0 = NU;

  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    master->Vz[c] = master->vz0;
    master->Kz[c] = master->kz0;
  }
}

/* END closure_k_e_implicit_init()                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the k-e closure scheme.                                           */
/*-------------------------------------------------------------------*/
void closure_k_e_implicit(geometry_t *window,  /* Window geometry    */
                          window_t *windat,    /* Window data        */
                          win_priv_t *wincon   /* Window constants   */
  )
{
  int c, c2, cc, k;             /* Cell coodinate / counter          */

  int cs, ks;                   /* 2D sparse coordinate */
  int cb, kb;                   /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */
  int winsize = window->szc;    /* 3D cell work array size           */
  int nzsize = window->nz + 1;  /* 1D work array size */
  double z0t, z0;
  double z0b;
  double dzdt;
  double stke;
  double sdiss;
  double sqrttke;
  double v1, wflux, wave_induced, diss0, b1;

  /*-----------------------------------------------------------------*/
  /* Note: the 3D work arrays wincon->w# could be used for the */
  /* dummy arrays below, but execution speed is considerably faster */
  /* when work array access is sequential in memory, hence the */
  /* mapping to a contiguous vertical 1D work array wincon->v#. The */
  /* code below contains a mixture of variables accessed in the */
  /* sparse array, c, and local array, k, with mapping achieved via */
  /* k=window->s2k[c].  */
  int *ctp = window->nsur_t;    /* Surface cells to process */
  int vcs = window->b2_t;      /* Last index of surface cells */
  double *C = wincon->v1;       /* Main diagonal */
  double *Cm1 = wincon->v2;     /* Lower diagonal */
  double *Cp1 = wincon->v3;     /* Upper diagonal */
  double *rhs = wincon->v4;     /* Right hand side */
  double *dzface = wincon->v5;  /* Cell thickness for faces */
  double *tke_Splus = wincon->v6; /* Source term for tke */
  double *tke_Sminus = wincon->v7;  /* Sink term for tke */
  double *diss_Splus = wincon->v8;  /* Source term for diss */
  double *diss_Sminus = wincon->v9; /* Sink term for diss */
  double *dz = wincon->v10;     /* Cell thickness in local coords */
  double *N2 = wincon->w1;      /* Buoyancy frequency */
  double *tke = wincon->w2;     /* tke in local coordinates */
  double *diss = wincon->w3;    /* diss in local coordinates */
  double *du1 = wincon->w4;     /* Vertical u1 gradient */
  double *du2 = wincon->w5;     /* Vertical u2 gradient */
  double *Kz = wincon->w6;      /* Diffusion in local coordinates */
  double *Vz = wincon->w7;      /* Viscosity in local coordinates */
  double *rho = windat->dens_0; /* Pointer to potential density */
  double *wu = wincon->d1;
  double *wv = wincon->d2;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(N2, 0, winsize * sizeof(double));
  memset(tke, 0, winsize * sizeof(double));
  memset(diss, 0, winsize * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Get the east and north cell centered wind stress                */
  vel_cen(window, windat, wincon, windat->wind1, NULL, wu, wv, NULL, NULL, 1);

  /*-----------------------------------------------------------------*/
  /* Calculate gradients over the whole grid */
  /* Vertical u and v gradients */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    cb = window->bot_t[cc];
    zm1 = window->zm1[c];
    while (c < cb) {
      du1[c] = windat->u[c] - windat->u[zm1];
      du2[c] = windat->v[c] - windat->v[zm1];
      c = zm1;
      zm1 = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the mixing coefficients over the whole grid */
  for (cc = 1; cc <= vcs; cc++) {
    cs = c = ctp[cc];           /* Set cs to the surface sparse coordinate 
                                 */
    c2 = window->m2d[c];        /* 2D sparse location corresponding to c */
    ks = window->s2k[cs];       /* Surface local coordinate */
    zm1 = window->zm1[c];

    /*---------------------------------------------------------------*/
    /* One wet layer : no mixing allowed but may need to alter TKE */
    /* and dissipation values, and the particle tracking code may */
    /* still require vertical diffusion values.  */
    if (zm1 == window->zm1[zm1]) {
      k_e_1layer(window, windat, wincon, wu, wv, c, c2);
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
    /* Get the wave contribution to tke if required. Note, wave      */
    /* is input into the water column, hence the negative sign (+ve  */
    /* fluxes are out of the water column).                          */
    wflux = 0.0;
    if (wincon->waves & VERTMIX) {
      if (wincon->waves & JONES) {
	wflux = sqrt(sqrt(wu[c2] * wu[c2] + wv[c2] * wv[c2]) / RHO_0);
	wflux = -wincon->wave_alpha * wflux * wflux * wflux;
      }
    }

    /*---------------------------------------------------------------*/
    /* Map the sparse arrays to local coordinates.  */
    c = cs;
    for (k = ks; k >= kb; k--) {
      zm1 = window->zm1[c];
      dzface[k] = (k == kb) ? 0.5 * dz[k] : 0.5 * (dz[k - 1] + dz[k]);
      tke[k] = windat->tke[c];
      diss[k] = windat->diss[c];
      Vz[k] = windat->Vz[c];
      Kz[k] = windat->Kz[c];
      c = zm1;
    }
    memset(tke_Splus, 0, nzsize * sizeof(double));
    memset(tke_Sminus, 0, nzsize * sizeof(double));
    memset(diss_Splus, 0, nzsize * sizeof(double));
    memset(diss_Sminus, 0, nzsize * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Loop over the layer interfaces, calculating source/sink terms */
    /* in the TKE and dissipation equations.  */
    c = cb;
    do {
      double M2;
      double drho;
      double P;
      double B;
      double c3e;
      double plus;
      double minus;


      zm1 = c;
      c = window->zp1[c];
      k = window->s2k[c];

      /* Density gradient */
      drho = rho[c] - rho[zm1];

      /* Constant in dissipation equation - Burchard et al (1998) */
      /* Section 5.2. Value is positive for unstable stratification, */
      /* negative for stable stratification.  */
      c3e = (drho > 0.0) ? 1.0 : -0.4;

      /* Brunt Vaisala frequency (squared). Also used later to */
      /* calculate stability functions.  */
      N2[c] = -(wincon->g / RHO_0) * (drho / dzface[k]);

      /* Shear frequency squared. Take care to avoid division by */
      /* zero if this is used this to calculate a Richardson number. */
      M2 = (du1[c] * du1[c] + du2[c] * du2[c]) / (dzface[k] * dzface[k]);

      /* Interpolate TKE and dissipation from layer centres to layer */
      /* interface.  */
      stke = max(0.5 * (tke[k - 1] + tke[k]), wincon->min_tke);
      sdiss = max(0.5 * (diss[k - 1] + diss[k]), wincon->min_diss);

      /* Ensure minimum level of dissipation - Burchard et al (1998) */
      /* Eqn 16.  */
      if ((N2[c] > 0.0) && (sdiss * sdiss < 0.045 * stke * stke * N2[c]))
        sdiss = max(sqrt(0.045 * stke * stke * N2[c]), wincon->min_diss);

      /* Set the Malek Ghantous wave parameters if required          */
      /* http://www.nonlin-processes-geophys.net/21/325/2014/npg-21-325-2014.html */
      /* Wave orbital method, Babanin & Haus, (2009), JPO, 39        */
      wave_induced = 0.0;
      if (wincon->waves & VERTMIX) {
	if (wincon->waves & WOM) {
	  double angfreq = 2.0 * PI / windat->wave_period[c2];
	  /* Deep water dispersion relation */
	  double wk = angfreq * angfreq / wincon->g; 
	  b1 = 0.0;
	  if (wincon->wave_b1 > 0.0) {
	    /* Set constant b1 for babanin dissipation                 */
	    /*  Ian Young's value                                      */
	    b1=0.0014; 
	    /* Alex Babanin's value, from Ardhuin's data (2009)        */
	    b1=0.002;  
	    b1 = wincon->wave_b1;
	  } else if (wincon->wave_b1 < 0) {
	    /* Set quadratic-in-steepness b1 according to Babanin      */
	    /* (pc 2013).                                              */
	    b1 = 5.0 * (wk * windat->wave_amp[c2] / 2.0);
	    b1 = b1 * b1;
	  }
	  wave_induced = b1 * wk * pow(angfreq * windat->wave_amp[c2] /
				       2.0 * exp(wk * window->cellz[c]), 3.0);
	}
      }

      /* Shear production - Rodi (1984) Eqn 2.48, or Burchard et al */
      /* (1998) Eqn 9, with additional correction for internal wave */
      /* shear as discussed by Mellor (1989), p267.  */
      P = Vz[k] * (M2 + 0.7 * N2[c]) + wave_induced;
      if (windat->s_prod) windat->s_prod[c] = P;

      /* Buoyancy production - refs as above */
      B = -Kz[k] * N2[c];
      if (windat->b_prod) windat->b_prod[c] = B;

      /* TKE source/sink, splitting positive and negative parts */
      /* (Patankar 1980).  */
      if ((P + B) > 0) {
        plus = P + B;
        minus = sdiss;
      } else {
        plus = P;
        minus = sdiss - B;
      }

      /* Distribute TKE sources to adjacent layers */
      tke_Splus[k - 1] += 0.5 * plus;
      tke_Splus[k] += 0.5 * plus;
      tke_Sminus[k - 1] += 0.5 * minus;
      tke_Sminus[k] += 0.5 * minus;

      /* Dissipation source/sink, splitting positive and negative */
      /* parts (Patankar 1980).  */
      if ((C1E * P + c3e * B) > 0) {
        plus = (sdiss / stke) * (C1E * P + c3e * B);
        minus = (sdiss / stke) * C2E * sdiss;
      } else {
        plus = (sdiss / stke) * C1E * P;
        minus = (sdiss / stke) * (C2E * sdiss - c3e * B);
      }

      /* Distribute dissipation sources to adjacent layers */
      diss_Splus[k - 1] += 0.5 * plus;
      diss_Splus[k] += 0.5 * plus;
      diss_Sminus[k - 1] += 0.5 * minus;
      diss_Sminus[k] += 0.5 * minus;

    } while (c != cs);

    /*---------------------------------------------------------------*/
    /* Now all the source and sink terms have been calculated: set */
    /* up tri-diagonal systems of equations for TKE and dissipation, */
    /* with appropriate surface and bottom boundary conditions.  */
    z0t = wincon->zs;
    if (wincon->waves & VERTMIX) {
      if (wincon->waves & JONES)
	z0t = wincon->wave_hf * windat->wave_amp[cs];
    }
    z0 = z0t;
    z0t = max(0.5 * dz[ks] + z0t, wincon->hmin);

    /*---------------------------------------------------------------*/
    /* TKE bottom boundary condition - zero flux (Neumann condition) */
    dzdt = dz[kb] / windat->dt;
    Cm1[kb] = 0.0;
    Cp1[kb] = -Vz[kb + 1] / (dzface[kb + 1] * SIG_K);
    C[kb] = dzdt - Cp1[kb] + dz[kb] * tke_Sminus[kb] /
      max(tke[kb], wincon->min_tke);
    rhs[kb] = dzdt * tke[kb] + dz[kb] * tke_Splus[kb];

    /* TKE surface boundary condition - zero flux */
    /* GOTM manual; Eq. 127
    wflux = K_k(window, cs, wflux, z0t); */
    /* Get the surface boundary condition by differentiating the 
    power law relation (k - K(z0-z)^a_e) and evaluating at z=0.
    v1 = K_e(window, cs, wflux, z0);
    wflux = -a_e * v1 * pow(z0, a_e - 1.0);
    */
    /* Jones & Monosmith Eq. 9 use wflux on the RHS. */

    dzdt = dz[ks] / windat->dt;
    Cm1[ks] = -Vz[ks] / (dzface[ks] * SIG_K);
    Cp1[ks] = 0.0;
    C[ks] = dzdt - Cm1[ks] + dz[ks] * tke_Sminus[ks] /
      max(tke[ks], wincon->min_tke);
    rhs[ks] = dzdt * tke[ks] + dz[ks] * tke_Splus[ks] - wflux;

    /* TKE mid-water column */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / windat->dt;
      Cm1[k] = -Vz[k] / (dzface[k] * SIG_K);
      Cp1[k] = -Vz[k + 1] / (dzface[k + 1] * SIG_K);
      C[k] = dzdt - Cm1[k] - Cp1[k] + dz[k] * tke_Sminus[k] /
        max(tke[k], wincon->min_tke);
      rhs[k] = dzdt * tke[k] + dz[k] * tke_Splus[k];
    }

    /* Solve tridiagonal system for TKE */
    tridiagonal(window, cs, cb, C, Cm1, Cp1, rhs, windat->tke,
                wincon->min_tke);

    /*---------------------------------------------------------------*/
    /* Dissipation bottom boundary condition - Dirichlet condition */
    /* specified from Burchard et al 1998 Eqns 18 and 11.  */
    z0b = max(0.5 * dz[kb], wincon->hmin);
    Cm1[kb] = 0.0;
    Cp1[kb] = 0.0;
    C[kb] = 1.0;
    rhs[kb] =
      max(CM0c * windat->tke[cb] * sqrt(windat->tke[cb]) /
          (z0b * KAPPA), wincon->min_diss);

    /* Dissipation surface boundary condition, Dirichlet condition. */
    /* For logarithmic boundary layers, see GOTM manual: Eq. 134.   */
    /* www.gotm.net/pages/documentation/manual/stable/pdf/a4.pdf    */
    /* For shear-free boundary layers, see GOTM manual: Eq. 136.    */
    /* www.gotm.net/pages/documentation/manual/stable/pdf/a4.pdf    */
    Cm1[ks] = 0.0;
    Cp1[ks] = 0.0;
    C[ks] = 1.0;

    if (wincon->waves & VERTMIX) {
      if (wincon->waves & JONES) {
	/* Shear-free boundary layers */
	v1 = K_e(window, cs, wflux, z0);
	diss0 = CM0c * v1 * sqrt(v1) * pow(z0t, 1.5*a_e-1.0) / L_e;
      }
    } else {
      /* Logarithmic boundary layers */
      diss0 = max(CM0c * windat->tke[cs] * sqrt(windat->tke[cs]) /
		  (z0t * KAPPA), wincon->min_diss);
    }
    rhs[ks] = max(diss0, wincon->min_diss);

    /* Dissipation mid-water column */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / windat->dt;
      Cm1[k] = -Vz[k] / (dzface[k] * SIG_E);
      Cp1[k] = -Vz[k + 1] / (dzface[k + 1] * SIG_E);
      C[k] = dzdt - Cm1[k] - Cp1[k] + dz[k] * diss_Sminus[k] /
        max(diss[k], wincon->min_diss);
      rhs[k] = dzdt * diss[k] + dz[k] * diss_Splus[k];
    }

    /* Solve tridiagonal system for dissipation */
    tridiagonal(window, cs, cb, C, Cm1, Cp1, rhs, windat->diss, wincon->min_diss);

    /*---------------------------------------------------------------*/
    /* Loop over the layers, calculating new diffusivity values.  */
    /* This can be done in two ways: */
    /* 1. At the layer interfaces directly, which requires averaging */
    /* the TKE and dissipation values from adjacent layers.  */
    /* 2. At the layer centres, which means that the diffusivity */
    /* values need to be interpolated to layer interfaces after the */
    /* event.  */
    /* Method 1 gives very good agreement with Burchard et al (1998) */
    /* test cases at all layers.  */
    /* Method 2 gives good agreement also, except for minor */
    /* variations near the top and bottom boundary.  */
    /* Use method 1 here.  */
    c = cb;
    do {
      double L, aN, aM, cmu, cmu_d, Lmax;
      zm1 = c;
      c = window->zp1[c];
      k = window->s2k[c];

      stke = 0.5 * (windat->tke[zm1] + windat->tke[c]);
      sqrttke = sqrt(stke);
      sdiss = 0.5 * (windat->diss[zm1] + windat->diss[c]);

      /* Turbulent length scale - Burchard et al (1998) Eqn 11.  */
      L = CM0 * CM0 * CM0 * stke * sqrttke / sdiss;
      Lmax = 0.56 * stke / N2[c];
      L = (N2[c] > 0) ? min(L, sqrt(Lmax)) : L;
      if (windat->L)
        windat->L[c] = L;

      /* Stability function parameter - Burchard et al (1998) Eqn 27 */
      aN = L * L * N2[c] / stke;
      aM = L * L * (du1[c] * du1[c] + du2[c] * du2[c]) / (dzface[k] * dzface[k] * stke);
      wincon->s_func(aN, aM, &cmu, &cmu_d);

      /* Updated viscosity and diffusivity */

      windat->Vz[c] = min(MAXVZ, max(cmu * sqrttke * L, wincon->vz0));
      windat->Kz[c] = min(MAXKZ, max(cmu_d * sqrttke * L, wincon->kz0));

      /* Wave enhanced mixing for BV nethod                          */
      if (wincon->waves & VERTMIX) {
	if (wincon->waves & BVM) {
	  double bv = get_bv_mono(window, windat, wincon, c, cs);
	  windat->Vz[c] += bv;
	  windat->Kz[c] += bv;
	}
      }
    } while (c != cs);

    /* Set Vz and Kz at bottom to be the same as at the next highest */
    /* interface. These values are never used in the mixing code, */
    /* but may be used by the particle tracking code.  */
    windat->Vz[cb] = windat->Vz[window->zp1[cb]];
    windat->Kz[cb] = windat->Kz[window->zp1[cb]];
  }
}

/* END closure_k_e_implicit()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* This routine deals with the case where there is only 1 wet layer. */
/*-------------------------------------------------------------------*/
void k_e_1layer(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon, /* Window geometry / constants */
		double *wu,
		double *wv,
                int c,          /* Sparse coordinate */
                int cs          /* 2D sparse coordinate */
  )
{

  double ustr_bot, ustr_surf, ustr;
  double z0;
  double tke;
  double diss;

  /* Surface stress */
  ustr_surf = sqrt(sqrt(wu[cs] * wu[cs] + wv[cs] * wv[cs]) / RHO_0);

  /* Bottom stress */
  ustr_bot = sqrt(wincon->Cd[cs] * (windat->u[c] * windat->u[c] + 
				    windat->v[c] * windat->v[c]));

  /* Pick maximum of surface and bottom friction velocity */
  ustr = max(ustr_surf, ustr_bot);

  /* z0 value ??? */
  z0 = max(0.5 * wincon->dz[c], wincon->zs);

  /* TKE value */
  tke = max(ustr * ustr / (CM0 * CM0), wincon->min_tke);

  /* Dissipation value */
  diss = max(ustr * ustr * ustr / (KAPPA * z0), wincon->min_diss);
  windat->tke[c] = tke;
  windat->diss[c] = diss;

  /* Vertical viscosity and diffusivity - Burchard et al (1998) */
  /* Eqns 11, 14, 30 */
  windat->Vz[c] = CM0 * CM0 * CM0 * CM0 * tke * tke / diss;
  windat->Kz[c] = 0.6985 * CM0 * CM0 * CM0 * tke * tke / diss;
}

/* END k_e_1layer()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute the constant K for shear-free boundary layers  */
/* (e.g. wave forced).                                               */
/*-------------------------------------------------------------------*/
double K_e(geometry_t *window, int cs, double wflux, double z0)
{
  window_t *windat = window->windat;
  double stke, sdiss, tkeL, cmu; 

  stke = 0.5 * (windat->tke[cs] + windat->tke[window->zm1[cs]]);
  sdiss = 0.5 * (windat->diss[cs] + windat->diss[window->zm1[cs]]);
  /* Turbulent length scale * sqrt(tke), Burchard et al (1998) Eq 11 */
  tkeL = CM0 * CM0 * CM0 * stke * stke / sdiss;
  /* Stability function                                              */
  cmu = windat->Vz[cs] / tkeL;
  /* Function K, Jones and Monosmith (2008), Eq. 19                  */
  /* Note: SIG_K is solved implicity on the LHS.                     */
  return(pow(-wflux / (cmu * a_e * L_e), 0.6667) / pow(z0, a_e));
}

/* END K_e()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute the Neumann surface condition for shear-free   */
/* boundary layers (e.g. wave forced) according to Eq. 127 in the    */
/* GOTM manual:                                                      */
/* http://www.gotm.net/pages/documentation/manual/stable/pdf/a4.pdf  */
 /* Substitute Eq. 128 into Eq. 127 reduces to wflux                 */
/*-------------------------------------------------------------------*/
double K_k(geometry_t *window, int cs, double wflux, double z0)
{
  window_t *windat = window->windat;
  double stke, sdiss, tkeL, cmu, K; 

  stke = 0.5 * (windat->tke[cs] + windat->tke[window->zm1[cs]]);
  sdiss = 0.5 * (windat->diss[cs] + windat->diss[window->zm1[cs]]);
  /* Turbulent length scale * sqrt(tke), Burchard et al (1998) Eq 11 */
  tkeL = CM0 * CM0 * CM0 * stke * stke / sdiss;
  /* Stability function                                              */
  cmu = windat->Vz[cs] / tkeL;
  /* Function K, Jones and Monosmith (2008), Eq. 19                  */
  /* Note: SIG_K is solved implicity on the LHS.                     */
  K = pow(-wflux / (cmu * a_e * L_e), 0.6667) / pow(z0, a_e);
  return(cmu * pow(K, 1.5) * L_e * a_e * pow(z0, 1.5 * a_e) / SIG_K);
}

/* END K_k()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the wave induced mixing coefficient, Bv.                 */
/* Qiao et al, (2009), Ocean Dynamics, 1339-1355.                    */
/* Qiao et al, (2013), JGR, 118, 4514 - 4524                         */
/*-------------------------------------------------------------------*/
double get_bv_mono(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon, /* Window constants           */
		   int c,              /* Sparse coordinate          */
		   int cs              /* 2D sparse coordinate       */
  )
{
  double bv, d1;
  double A = windat->wave_amp[cs];               /* Amplitude        */
  double w = 2.0 * PI / windat->wave_period[cs]; /* Angular freq     */
  double k = w * w / wincon->g;                  /* Wave number      */
  double cp = w / k;                             /* Phase velocity   */
  double stokes = cp * A * A * k * k;
  double depth = windat->eta[cs] - window->botz[cs];

  d1 = sinh(k * (depth + window->cellz[c])) / sinh(k * depth);
  bv = wincon->wave_alpha * A * stokes * pow(d1, 3.0);
  return(bv);
}

/* END get_bv_mono()                                                 */
/*-------------------------------------------------------------------*/
