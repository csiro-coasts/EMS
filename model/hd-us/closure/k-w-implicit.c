/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/closure/k-w-implicit.c
 *  
 *  Description:
 *  Implements K-omega mixing scheme. The code
 *  here is based primarily on equations described
 *  in the following references:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: k-w-implicit.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Constants                                                         */
#define VK		(VON_KAR)         /* Von Karman's constant */
#define KAPPA (0.4)             /* Von Karman's constant */

/* Minimum possible viscosity/diffusivity (akin to molecular diffusion) */
#define NU		(1.0e-5)

/* Maximum possible viscosity/diffusivity (for stability)            */
#define MAXVZ		(1.0)
#define MAXKZ		(1.0)

/* Umlauf et al (2003) Table 1                                       */
#define CM0		(0.5562)
#define CM0c            (CM0*CM0*CM0)

/* Umlauf et al (2003) Table 1                                       */
#define C1W             (0.52)

/* Umlauf et al (2003) Table 1                                       */
#define C2W             (0.8)

/* Umlauf et al (2003) Table 6, using KW88 (C3E=C3W+1)               */
#define C3W             (-0.642)

/* Umlauf et al (2003) Table 1                                       */
#define SIG_K    	(2.0)

/* Umlauf et al (2003) Table 1                                       */
#define SIG_W    	(2.0)

/* Minimum turbulent kinetic energy - Burchard et al (1998) Eqn 17   */
#define MIN_TKE          (7.6e-6)

/* Minimum dissipation - from CRS code                               */
/*#define MIN_DISS         (5e-10)*/
/* From Warner et al (2005), Table 1 and Umlauf et al eqn. 10        */
/* i.e. min_omega = 1e-12                                            */
#define MIN_DISS         (7.27e-19)

/* Mean density                                                      */
#define RHO_0		(1025.0)

/* Stability function constants, Canuto et al (2001) Eqn 18,19,20b   */
#define L1              (0.1070)
#define L2              (0.0032)
#define L3              (0.0864)
#define L4              (0.1200)
#define L5              (11.900)
#define L6              (0.4000)
#define L7              (0.0000)
#define L8              (0.4800)

#define s0              (1.5*L1*L5*L5)
#define s1              (-L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8)
#define s2              (-0.375*L1*(L6*L6-L7*L7))
#define s4              (2.*L5)
#define s5              (2.*L4)
#define s6              (2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7))
#define d0              (3.*L5*L5)
#define d1              (L5*(7.*L4+3.*L8))
#define d2              (L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7))
#define d3              (L4*(4.*L4+3.*L8))
#define d4              (L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2))
#define d5              (0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7))

/* Wave parameters */
#define a_w  (-2.53)  /* Jones and Monosmith, JGR, (2008, p5) */
#define L_w  (0.25)   /* Jones and Monosmith, JGR, (2008, p5) */

void k_w_1layer(geometry_t *window, window_t *windat, win_priv_t *wincon,
                double *wu, double *wv, int c, int cs);
void fcu_f(geometry_t *window, window_t *windat, win_priv_t *wincon,
           double *chik);
void fcw_f(geometry_t *window, window_t *windat, win_priv_t *wincon,
           double *chiw);
double K_w(geometry_t *window, int cs, double wflux, double N2, double fcu, double z0);

/*-------------------------------------------------------------------*/
/* Routine to initialize the implicit k-e mixing scheme              */
/*-------------------------------------------------------------------*/
void closure_k_w_implicit_init(geometry_t *geom,  /* Sparse global
                                                     geometry structure */
                               master_t *master /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coordinates */

  /* Sanity check - make sure appropriate variables are present */
  if (!master->tke || !master->omega)
    hd_quit
      ("k-omega mixing scheme requires tracers called tke and omega\n");

  /* In the k-w mixing model, TKE and omega are vertically */
  /* mixed using the eddy viscosity (divided by the Schmidt numbers */
  /* sigma_k and sigma_e respectively - see Umlauf et al 1998 Eqns */
  /* 7 and 8). For this version, diffusion is handled in the code */
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

/* END closure_k_w_implicit_init()                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the k-w closure scheme.                                           */
/*-------------------------------------------------------------------*/
void closure_k_w_implicit(geometry_t *window,  /* Window geometry    */
                          window_t *windat,    /* Window data        */
                          win_priv_t *wincon   /* Window constants   */
  )
{
  int c, c2, cc, k;             /* Sparse coodinate / counter */
  int cs, ks;                   /* 2D sparse coordinate */
  int cb, kb;                   /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */
  int winsize = window->szc;    /* 3D work array size */
  int nzsize = window->nz + 1;  /* 1D work array size */
  double z0t, z0;
  double z0b;
  double dzdt;
  double stke;
  double somega;
  double diss;
  double du1, du2;
  double omin;
  double min_omega;
  double wflux, diss0;
  double wave_induced, b1;

  /*-----------------------------------------------------------------*/
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
  double *tke_Splus = wincon->v6; /* Source term for tke */
  double *tke_Sminus = wincon->v7;  /* Sink term for tke */
  double *omega_Splus = wincon->v8; /* Source term for diss */
  double *omega_Sminus = wincon->v9;  /* Sink term for diss */
  double *dz = wincon->v10;     /* Cell thickness in local coords */
  double *N2 = wincon->w1;      /* Buoyancy frequency */
  double *M2 = wincon->w10;     /* Velocity shear */
  double *tke = wincon->w2;     /* tke in local coordinates */
  double *omega = wincon->w3;   /* diss in local coordinates */
  double *Kz = wincon->w4;      /* Diffusion in local coordinates */
  double *Vz = wincon->w5;      /* Viscosity in local coordinates */
  double *rho = windat->dens_0; /* Pointer to potential density */
  double *fcu = wincon->w6;     /* Parameter function fcu */
  double *fcw = wincon->w7;     /* Parameter function fcw */
  double *wu = wincon->d6;
  double *wv = wincon->d7;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(N2, 0, winsize * sizeof(double));
  memset(M2, 0, winsize * sizeof(double));

  /* Minimum dissipation - using Umlauf et al Eqn 10, MIN_TKE,       */
  /* MIN_DISS and fcu=1.                                             */
  min_omega = wincon->min_diss / (CM0c * CM0 * wincon->min_tke);

  /*-----------------------------------------------------------------*/
  /* Get the parameter function for the loss term of omega */
  fcu_f(window, windat, wincon, fcu);
  fcw_f(window, windat, wincon, fcw);

  /*-----------------------------------------------------------------*/
  /* Get the east and north cell centered wind stress                */
  vel_cen(window, windat, wincon, windat->wind1, wu, wv, NULL, NULL, 1);

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
      k_w_1layer(window, windat, wincon, wu, wv, c, cs);
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
      omega[k] = windat->omega[c];
      Vz[k] = windat->Vz[c];
      Kz[k] = windat->Kz[c];
      c = zm1;
    }

    memset(tke_Splus, 0, nzsize * sizeof(double));
    memset(tke_Sminus, 0, nzsize * sizeof(double));
    memset(omega_Splus, 0, nzsize * sizeof(double));
    memset(omega_Sminus, 0, nzsize * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Loop over the layer interfaces, calculating source/sink terms */
    /* in the TKE and dissipation equations.  */
    c = cb;
    do {
      double drho;
      double P;
      double B;
      double c3w;
      double plus;
      double minus;

      zm1 = c;
      c = window->zp1[c];
      k = window->s2k[c];
 
     /* Vertical u1 gradient */
      du1 = windat->u[c] - windat->u[zm1];

      /* Vertical u2 gradient */
      du2 = windat->v[c] - windat->v[zm1];

      /* Density gradient */
      drho = rho[c] - rho[zm1];

      /* Constant in dissipation equation, Umlauf et al (2003) Table */
      /* 6. Value is positive for unstable stratification, negative */
      /* for stable stratification. Note; this parameter must be */
      /* calculated using the stability functions of Canuto et al */
      /* (2001).  */
      c3w = (drho > 0.0) ? 1.0 : C3W;

      /* Brunt Vaisala frequency (squared). Also used later to */
      /* calculate stability functions.  */
      N2[c] = -(wincon->g / RHO_0) * (drho / dzface[k]);

      /* Shear frequency squared. Take care to avoid division by */
      /* zero if this is used this to calculate a Richardson number. */
      M2[c] = (du1 * du1 + du2 * du2) / (dzface[k] * dzface[k]);

      /* Interpolate TKE and dissipation from layer centres to layer */
      /* interface.  */
      stke = max(0.5 * (tke[k - 1] + tke[k]), wincon->min_tke);
      omin =
        max((N2[c] > 0) ? sqrt(N2[c]) / (fcu[c] * sqrt(0.56) * CM0) : 0,
            min_omega);
      somega = max(0.5 * (omega[k - 1] + omega[k]), omin);

      /* Dissipation, Umlauf et al (2003) Eqn 10 */
      diss = CM0 * CM0c * fcu[c] * stke * somega;

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
      P = Vz[k] * (M2[c] + 0.7 * N2[c]) + wave_induced;
      if (windat->s_prod) windat->s_prod[c] = P;

      /* Buoyancy production - refs as above */
      B = -Kz[k] * N2[c];
      if (windat->b_prod) windat->b_prod[c] = B;

      /* TKE source/sink, splitting positive and negative parts */
      /* (Patankar 1980).  */
      if ((P + B) > 0) {
        plus = P + B;
        minus = diss;
      } else {
        plus = P;
        minus = diss - B;
      }

      /* Distribute TKE sources to adjacent layers */
      tke_Splus[k - 1] += 0.5 * plus;
      tke_Splus[k] += 0.5 * plus;
      tke_Sminus[k - 1] += 0.5 * minus;
      tke_Sminus[k] += 0.5 * minus;

      /* Omega source/sink, splitting positive and negative parts */
      /* (Patankar 1980).  */
      if ((C1W * P + c3w * B) > 0) {
        plus = somega * (C1W * P + c3w * B) / stke;
        minus = somega * C2W * diss * fcw[c] / (fcu[c] * stke);
      } else {
        plus = somega * C1W * P / stke;
        minus = somega * (C2W * diss * fcw[c] / fcu[c] - c3w * B) / stke;
      }

      /* Distribute omega sources to adjacent layers */
      omega_Splus[k - 1] += 0.5 * plus;
      omega_Splus[k] += 0.5 * plus;
      omega_Sminus[k - 1] += 0.5 * minus;
      omega_Sminus[k] += 0.5 * minus;

    } while (c != cs);

    /*---------------------------------------------------------------*/
    /* Now all the source and sink terms have been calculated: set */
    /* up tri-diagonal systems of equations for TKE and omega, with */
    /* appropriate surface and bottom boundary conditions.  */

    /*---------------------------------------------------------------*/
    /* TKE bottom boundary condition - zero flux (Neumann condition) */
    dzdt = dz[kb] / windat->dt;
    Cm1[kb] = 0.0;
    Cp1[kb] = -Vz[kb + 1] / (dzface[kb + 1] * SIG_K);
    C[kb] = dzdt - Cp1[kb] + dz[kb] * tke_Sminus[kb] /
      max(tke[kb], wincon->min_tke);
    rhs[kb] = dzdt * tke[kb] + dz[kb] * tke_Splus[kb];

    /* TKE surface boundary condition - zero flux */
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
    /* Omega bottom boundary condition - Dirichlet condition */
    /* specified from www.gotm.net/code_html/node50.html (BCs for */
    /* the generic turbulence model).  */
    z0b = max(0.5 * dz[kb] + wincon->z0[c2], wincon->hmin);
    Cm1[kb] = 0.0;
    Cp1[kb] = 0.0;
    C[kb] = 1.0;
    rhs[kb] = max(sqrt(windat->tke[cb]) / (CM0 * z0b * KAPPA), min_omega);

    /* Omega surface boundary condition, Dirichlet condition.       */
    /* For logarithmic boundary layers, see GOTM manual:            */
    /* www.gotm.net/pages/documentation/manual/pdf/a4.pdf, Eq. 138. */
    /* For shear-free boundary layers, see GOTM manual:             */
    /* www.gotm.net/pages/documentation/manual/pdf/a4.pdf, Eq. 140. */
    z0t = wincon->zs;
    if (wincon->waves & VERTMIX) {
      if (wincon->waves & JONES)
	z0t = wincon->wave_hf * windat->wave_amp[cs];
    }
    z0 = z0t;
    z0t = max(0.5 * dz[ks] + z0t, wincon->hmin);

    Cm1[ks] = 0.0;
    Cp1[ks] = 0.0;
    C[ks] = 1.0;
    if (wincon->waves & VERTMIX) {
      if (wincon->waves & JONES) {
	/* Shear-free boundary layers */
	du1 = K_w(window, cs, wflux, N2[cs], fcu[cs], z0);
	diss0 = sqrt(du1) * pow(z0t, 0.5*a_w-1.0) / (CM0 * L_w);
      }
    } else {
      /* Logarithmic boundary layers */
      diss0 = sqrt(windat->tke[cs]) / (CM0 * z0t * KAPPA);
    }
    rhs[ks] = max(diss0, min_omega);

    /* Omega mid-water column */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / windat->dt;
      Cm1[k] = -Vz[k] / (dzface[k] * SIG_W);
      Cp1[k] = -Vz[k + 1] / (dzface[k + 1] * SIG_W);
      C[k] = dzdt - Cm1[k] - Cp1[k] + dz[k] * omega_Sminus[k] /
        max(omega[k], min_omega);
      rhs[k] = dzdt * omega[k] + dz[k] * omega_Splus[k];
    }

    /* Solve tridiagonal system for omega */
    tridiagonal(window, cs, cb, C, Cm1, Cp1, rhs, windat->omega, 0.0);

    /* Limit the turbulence frequency using Galperin et al (1988); */
    /* e.g. Burchard et al (1998) Eqn 16 with Umlauf et al (2003) */
    /* Eqn 10.  */
    for (c = cs; c <= cb; c = window->zm1[c]) {
      omin =
        max((N2[c] > 0) ? 0.045 * sqrt(N2[c]) / (fcu[c] * CM0c * CM0) : 0,
            min_omega);
      windat->omega[c] = max(omin, windat->omega[c]);
    }

    /*---------------------------------------------------------------*/
    /* Loop over the layers, calculating new diffusivity values at */
    /* the layer interfaces directly, which requires averaging the */
    /* TKE and omega values from adjacent layers.  */
    c = cb;
    do {

      double L, aN, aM, cmu, cmu_d, Lk, sqrttke, d, Lmax;
      double c3 = sqrt(2.0) / CM0c;

      zm1 = c;
      c = window->zp1[c];

      stke = 0.5 * (windat->tke[zm1] + windat->tke[c]);
      sqrttke = sqrt(stke);
      somega = 0.5 * (windat->omega[zm1] + windat->omega[c]);

      /* Turbulent length scale - Burchard et al (1998) Eqn 11 with */
      /* Umlauf et al (2003) Eqn 10. Limit the length scale using */
      /* Galperin et al (1988), Eqn 22.  */
      L = sqrttke / (CM0 * fcu[c] * somega);
      Lmax = sqrt(0.56 * stke / N2[c]);
      L = (N2[c] > 0) ? min(L, Lmax) : L;
      if (windat->L)
        windat->L[c] = L;

      /* Buoyancy and shear numbers, Burchard and Bolding (2001) Eqn 20 */
      Lk = L * L / stke;
      aN = N2[c] * Lk;
      aM = M2[c] * Lk;

      /* Note; the scaling c3 is required if mixing coefficients are */
      /* expressed in terms of tke and L (Canuto et al (2001) Eqn 16d) */
      /* rather than tke and diss.  */
      /* Stability functions accept aN and aM as arguments, and      */
      /* return cmu and cmu_d rather than SH and SM.                 */
      wincon->s_func(aN, aM, &cmu, &cmu_d);
      cmu = c3 * cmu;
      cmu_d = c3 * cmu_d;

      /* Updated viscosity and diffusivity */
      windat->Vz[c] = min(MAXVZ, max(cmu * sqrttke * L, wincon->vz0));
      windat->Kz[c] = min(MAXKZ, max(cmu_d * sqrttke * L, wincon->kz0));

      /* Wave enhanced mixing for BV nethod                          */
      if (wincon->waves & VERTMIX && wincon->waves & BVM) {
	double bv = get_bv_mono(window, windat, wincon, c, cs);
	windat->Vz[c] += bv;
	windat->Kz[c] += bv;
      }

    } while (c != cs); 

    /* Set Vz and Kz at bottom to be the same as at the next highest */
    /* interface. These values are never used in the mixing code, */
    /* but may be used by the particle tracking code.  */
    windat->Vz[cb] = windat->Vz[window->zp1[cb]];
    windat->Kz[cb] = windat->Kz[window->zp1[cb]];
  }
}

/* END closure_k_w_implicit()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* This routine deals with the case where there is only 1 wet layer. */
/*-------------------------------------------------------------------*/
void k_w_1layer(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon, /* Window geometry / constants */
		double *wu,
		double *wv,
                int c,          /* Sparse coordinate */
                int cs          /* 2D sparse coordinate */
  )
{
  double v1, v2;
  double ustr_bot, ustr_surf, ustr;
  double z0;
  double tke;
  double omega;
  double min_omega;

  /* Minimum dissipation - using Umlauf et al Eqn 10, MIN_TKE,       */
  /* MIN_DISS and fcu=1.                                             */
  min_omega = wincon->min_diss / (CM0c * CM0 * wincon->min_tke);

  ustr_surf = sqrt(sqrt(wu[cs] * wu[cs] + wv[cs] * wv[cs]) / RHO_0);

  /* Bottom stress */
  ustr_bot = sqrt(wincon->Cd[cs] * (windat->u[c] * windat->u[c] + 
				    windat->v[c] * windat->v[c]));

  /* Pick maximum of surface and bottom friction velocity */
  ustr = max(ustr_surf, ustr_bot);

  /* z0 value ??? */
  z0 = max(0.5 * wincon->dz[c], wincon->hmin);

  /* TKE value */
  tke = max(ustr * ustr / (CM0 * CM0), wincon->min_tke);

  /* Dissipation value */
  omega = max(ustr * ustr * ustr / (KAPPA * z0), min_omega);
  windat->tke[c] = tke;
  windat->omega[c] = omega;

  /* Vertical viscosity and diffusivity - Burchard et al (1998) */
  /* Eqns 11, 14, 30 */
  windat->Vz[c] = 2.0 * s0 * tke / (CM0c * CM0 * d0 * omega);
  windat->Kz[c] = 2.0 * s4 * tke / (CM0c * CM0 * d0 * omega);
}

/* END k_w_1layer()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the function fcu, Umlauf et al (2003) eqn 11 */
/*-------------------------------------------------------------------*/
void fcu_f(geometry_t *window,  /* Processing window */
           window_t *windat,    /* Window data structure */
           win_priv_t *wincon,  /* Window geometry / constants */
           double *chik         /* Return parameter */
  )
{
  int c, cs, cc, c1, c2;        /* Sparse coordinates */
  int ee, e, es;
  int p1, m1, p1m, m1m;         /* Directional maps */
  double fcu;                   /* Return value, chik>0 */
  double *tke;                  /* Turbulent kinetic energy */
  double *omega;                /* Turbulent frequency */
  double dzz;                   /* Cell spacing between layers */
  double chiks;                 /* chik squared */
  double *tke_g = wincon->w2;
  double *omega_g = wincon->w3;

  /* W88 model */
  if (!wincon->fcf) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      chik[c] = 1.0;
    }
    return;
  }

  not_implemented("k-w fcu_f function");
  hd_warn("Use W88 model!");
  return;

  /* Set pointers and initialise */
  tke = windat->tke;
  omega = windat->omega;

  /* Get the gradient of tke and omega at edges */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    chik[e] = (tke[c1] - tke[c2]) / window->h2au1[es];
    omega_g[e] = (omega[c1] - omega[c2]) / window->h2au1[es];
  }
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    tke_g[c] = 0.0;
    for (ee = 1; ee <= window->npe[cs]; ee++)
      tke_g[c] += chik[e] / (double)window->npe[cs];
  }
  memcpy(chik, omega_g, window->sze * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    omega_g[c] = 0.0;
    for (ee = 1; ee <= window->npe[cs]; ee++)
      omega_g[c] += chik[e] / (double)window->npe[cs];
  }

  /* Calculate the function fcu */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    p1 = window->xp1[c];
    m1 = window->zm1[c];
    p1m = window->zm1[p1];
    m1m = window->zm1[m1];
    chik[c] = 0.5 * (tke_g[c] + tke_g[m1]) + 0.5 * (omega_g[c] + omega_g[m1]);
    m1 = window->zm1[c];
    dzz = 0.5 * (wincon->dz[c] * wincon->dz[m1]);
    chik[c] += (tke[c] - tke[m1]) * (omega[c] - omega[m1]) / (dzz * dzz);
    chik[c] /= (omega[c] * omega[c] * omega[c]);
    chiks = chik[c] * chik[c];
    fcu = (1.0 + 680.0 * chiks) / (1.0 + 400.0 * chiks);
    chik[c] = (chik[c] <= 0.0) ? 1.0 : fcu;
    /* chik[c]=1.0; */
  }
}

/* END fcu()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the function fcw, Umlauf et al (2003) eqn 12 */
/*-------------------------------------------------------------------*/
void fcw_f(geometry_t *window,  /* Processing window */
           window_t *windat,    /* Window data structure */
           win_priv_t *wincon,  /* Window geometry / constants */
           double *chiw         /* Return parameter */
  )
{
  int c, cc, c2, cs, cb;        /* Sparse coordinates */
  int xp1, xm1, yp1, ym1;       /* Directional maps */
  int zp1, zm1, xpym1, xmyp1;   /* Directional maps */
  double *tke;                  /* Turbulent kinetic energy */
  double *omega;                /* Turbulent frequency */
  double *u1, *u2, *w;          /* Pointers to velocity */
  double *u1au2, *u2au1;        /* Velocity at opposite cell faces */
  double wtop, wbot;            /* w at zp1 and zm1 */
  double lat;                   /* Latitude */
  double om2, om3, om = 7.292e-5; /* Rotational components (local coords) */
  double w12, w13, w23;
  double s12, s13, s23;
  double dudx, dvdy, dwdz;      /* Velocity gradients */
  double dudy, dudz, dvdx, dvdz;  /* Velocity gradients */
  double dwdx, dwdy;            /* Velocity gradients */
  double dzz;                   /* Cell spacing between layers */

  /* Initialise */

  /* W88 model */
  if (!wincon->fcf) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      chiw[c] = 1.0;
    }
    return;
  }

  not_implemented("k-w fcw_f function");
  hd_warn("Use W88 model!");
  return;

  /* Set pointers and initialise */
  tke = windat->tke;
  omega = windat->omega;
  u1 = windat->u1;
  u2 = windat->u2;
  w = windat->w;
  u1au2 = wincon->w8;
  u2au1 = wincon->w9;
  memset(chiw, 0, window->sgsiz * sizeof(double));

  /* Get u1 at e2 faces */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    u1au2[c] = 0.25 * (u1[ym1] + u1[xpym1] + u1[c] + u1[xp1]);
  }

  /* Get u2 at e1 faces */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    u2au1[c] = 0.25 * (u2[xm1] + u2[xmyp1] + u2[c] + u2[yp1]);
  }

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = cs  = window->nsur_t[cc];
    c2 = window->m2d[c];
    cb = window->bot_t[cc];
    lat = asin(wincon->u1c5[c2] / (2.0 * om));
    om2 = om * cos(lat);
    om3 = om * sin(lat);

    while (c <= cb) {
      zm1 = window->zm1[c];
      zp1 = window->zp1[c];
      xp1 = window->xp1[c];
      xm1 = window->xm1[c];
      yp1 = window->yp1[c];
      ym1 = window->ym1[c];
      dzz = 0.5 * (wincon->dz[c] * wincon->dz[zm1]);
      
      /* All velocity gradients are defined at cell faces */
      /* Horizontal u1 gradient in e1 direction */
      dudx = (0.5 * (u1[xp1] - u1[c] + u1[window->zm1[xp1]] - u1[zm1])) /
	window->h1acell[c2];
      
      /* Horizontal u1 gradient in e2 direction */
      dudy =
	(0.5 *
	 (u1au2[yp1] - u1au2[c] + u1au2[window->zm1[yp1]] -
	  u1au2[zm1])) / window->h2acell[c2];
      
      /* Horizontal u2 gradient in e2 direction */
      dvdy = (u2[yp1] - u2[c] + u2[window->zm1[yp1]] - u2[zm1]) /
	window->h2acell[c2];
      
      /* Horizontal u2 gradient in e2 direction */
      dvdx =
	(0.5 *
	 (u2au1[xp1] - u2au1[c] + u2au1[window->zm1[xp1]] -
	  u2au1[zm1])) / window->h1acell[c2];

      /* Horizontal w gradient in e1 direction */
      dwdx = (w[xp1] - w[xm1]) / (2.0 * window->h1acell[c2]);

      /* Horizontal w gradient in e2 direction */
      dwdy = (w[yp1] - w[ym1]) / (2.0 * window->h2acell[c2]);

      /* Vertical u1 gradient */
      dudz =
	(0.5 * (u1[c] - u1[zm1] + u1[xp1] - u1[window->zm1[xp1]])) / dzz;

      /* Vertical u2 gradient */
      dvdz =
	(0.5 * (u2[c] - u2[zm1] + u2[yp1] - u2[window->zm1[yp1]])) / dzz;

      /* Vertical w gradient */
      wtop = (c == cs) ? windat->wtop[c2] : w[zp1];
      wbot = (c == cb) ? windat->wbot[c2] : w[zm1];
      dwdz = 0.5 * (wtop - wbot) / dzz;

      w12 = 0.5 * (dudy - dvdx) - om3;
      w13 = 0.5 * (dudz - dwdx) + om2;
      w23 = 0.5 * (dvdz - dwdy);

      s12 = 0.5 * (dudy + dvdx);
      s13 = 0.5 * (dudz + dwdx);
      s23 = 0.5 * (dvdz + dwdy);

      chiw[c] += (-w12 * w12 * dvdy + w12 * w23 * s13);
      chiw[c] += (-w13 * w12 * s23 - w13 * w13 * dwdz);
      chiw[c] += (-w12 * w12 * dudx - w12 * w13 * s23);
      chiw[c] += (w23 * w12 * s13 - w23 * w23 * dwdz);
      chiw[c] += (-w13 * w13 * dudx - w13 * w23 * s12);
      chiw[c] += (-w23 * w13 * s12 - w23 * w23 * dvdy);
      chiw[c] =
	fabs(chiw[c] / (pow(CM0, 12) * omega[c] * omega[c] * omega[c]));
      chiw[c] = (1.0 + 70.0 * chiw[c]) / (1.0 + 80.0 * chiw[c]);
      /* chiw[c]=1.0; */
      c = zm1;
    }
  }
}

/* END fcw()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute the constant K for shear-free boundary layers  */
/* (e.g. wave forced).                                               */
/*-------------------------------------------------------------------*/
double K_w(geometry_t *window, int cs, double wflux, double N2, double fcu, double z0)
{
  window_t *windat = window->windat;
  double stke, somega, tkeL, L, Lmax, cmu; 

  stke = 0.5 * (windat->tke[cs] + windat->tke[window->zm1[cs]]);
  somega = 0.5 * (windat->omega[cs] + windat->omega[window->zm1[cs]]);
  /* Turbulent length scale - Burchard et al (1998) Eqn 11 with */
  /* Umlauf et al (2003) Eqn 10.  Limit the length scale using */
  /* Galperin et al (1988), Eqn 22.  */
  L = sqrt(stke) / (CM0 * fcu * somega);
  Lmax = sqrt(0.56 * stke / N2);
  L = (N2 > 0) ? min(L, Lmax) : L;
  /* Stability function                                              */
  cmu = windat->Vz[cs] / sqrt(stke) * L;
  /* Function K, Kones and Monosmith (2008), Eq. 19                  */
  /* Note: SIG_K is solved implicity on the LHS.                     */
  return(pow(-wflux / (cmu * a_w * L_w), 0.6667) / pow(z0, a_w));
}

/* END K_w()                                                         */
/*-------------------------------------------------------------------*/
