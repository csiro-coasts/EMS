/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/closure/closure.c
 *  
 *  Description:
 *  Selects the vertical mixing scheme to be used.
 *  Note the table of pointers to routines which
 *  implement the various vertical mixing schemes
 *  in meco. These routines are found in the other
 *  source files in this directory.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: closure.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/* Declarations for routines in other files */
void closure_csanady_init(geometry_t *geom, master_t *master);
void closure_csanady(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);

void closure_MY2_init(geometry_t *geom, master_t *master);
void closure_MY2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void closure_MY2_est(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);

void closure_MY2_5_init(geometry_t *geom, master_t *master);
void closure_MY2_5(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);

void closure_MY2_5_HAR_init(geometry_t *geom, master_t *master);
void closure_MY2_5_HAR(geometry_t *window, window_t *windat,
		       win_priv_t *wincon);

void closure_constant_init(geometry_t *geom, master_t *master);
void closure_constant(geometry_t *window, window_t *windat,
                      win_priv_t *wincon);

void closure_k_e_implicit_init(geometry_t *geom, master_t *master);
void closure_k_e_implicit(geometry_t *window, window_t *windat,
                          win_priv_t *wincon);

void closure_k_w_implicit_init(geometry_t *geom, master_t *master);
void closure_k_w_implicit(geometry_t *window, window_t *windat,
                          win_priv_t *wincon);

/*-------------------------------------------------------------------*/
/* This routine reads the parameter file to determine the vertical   */
/* mixing scheme to be used, sets the appropriate pointers, and      */
/* calls the mixing scheme initialisation and calculation routines.  */
/*-------------------------------------------------------------------*/
void closure_init(parameters_t *params, /* Parameter data structure */
                  master_t *master  /* Master data structure */
  )
{
  /* Table of mixing scheme parameter names, long names, and */
  /* pointers to routines which implement them. Note that there may */
  /* may be more than one table entry per mixing scheme, to allow a */
  /* scheme to have several names, any of which may be used in the */
  /* the parameter file.  */
  struct {
    char *name;
    char *description;
    void (*init) (geometry_t *geom, master_t *master);
    void (*calc) (geometry_t *window, window_t *windat,
                  win_priv_t *wincon);
    void (*stab) (double aN, double aM, double *cmu, double *cmu_d);
  } mixing_list[] = {
    {
    "constant",
      "Constant mixing value",
      closure_constant_init, closure_constant, s_null}, {
    "csanady",
      "Csanady-Richardson mixing scheme",
      closure_csanady_init, closure_csanady, s_null}, {
    "MY2",
      "Mellor-Yamada level 2 mixing scheme",
    closure_MY2_init, closure_MY2, s_null}, {
    "MY2_5",
      "Mellor-Yamada level 2.5 mixing scheme",
      closure_MY2_5_init, closure_MY2_5, s_pom}, {
    "MY2_5_HAR",
      "Harcourt (2015) modified Mellor-Yamada level 2.5 mixing scheme",
      closure_MY2_5_HAR_init, closure_MY2_5_HAR, s_null}, {
    "harcourt",
      "Harcourt (2015) modified Mellor-Yamada level 2.5 mixing scheme",
      closure_MY2_5_HAR_init, closure_MY2_5_HAR, s_null}, {
    "mellor_yamada_2_0",
      "Mellor-Yamada level 2 mixing scheme",
      closure_MY2_init, closure_MY2, s_null}, {
    "mellor_yamada_2_5",
      "Mellor-Yamada level 2.5 mixing scheme",
      closure_MY2_5_init, closure_MY2_5, s_pom}, {
    "mellor_yamada_2_0_estuarine",
      "Mellor-Yamada level 2 mixing scheme : improved mixing length",
      closure_MY2_init, closure_MY2_est, s_null}, {
    "k-e-implicit",
      "k-epsilon mixing scheme, implicit",
      closure_k_e_implicit_init, closure_k_e_implicit, s_galperin}, {
    "k-e",
      "k-epsilon mixing scheme, implicit",
      closure_k_e_implicit_init, closure_k_e_implicit, s_galperin}, {
    "W88",
      "Wilcox (1988)",
      closure_k_w_implicit_init, closure_k_w_implicit, s_galperin}, {
    "k-w",
      "k-omega mixing scheme, implicit",
    closure_k_w_implicit_init, closure_k_w_implicit, s_canutoA}, {
      NULL, NULL, NULL, NULL, NULL}
  };

  void (*init) (geometry_t *geom, master_t *master) = NULL;
  int i, n;

  /*-----------------------------------------------------------------*/
  /* Read the parameter file to determine the mixing scheme used */
  prm_set_errfn(hd_silent_warn);
  if (params->mixsc != NULL) {
    /* Search for the mixing scheme name in the above table */
    for (i = 0; (init == NULL) && mixing_list[i].name; ++i) {
      if (strcasecmp(params->mixsc, mixing_list[i].name) == 0) {
        init = mixing_list[i].init;
	master->s_func = mixing_list[i].stab;
	master->calc_closure = mixing_list[i].calc;
      }
    }
  } else
    strcpy(params->mixsc, "");

  /* Set the stability function if required */
  if (params->s_func != NULL) {
    if (strcmp(params->s_func, "CANUTO_A") == 0)
      master->s_func = s_canutoA;
    if (strcmp(params->s_func, "CANUTO_B") == 0)
      master->s_func = s_canutoA;
    if (strcmp(params->s_func, "KANTHA&CLAYSON") == 0)
      master->s_func = s_kantha_clayson;
    if (strcmp(params->s_func, "GALPERIN") == 0)
      master->s_func = s_galperin;
    if (strcmp(params->s_func, "POM") == 0)
      master->s_func = s_pom;
    if (strcmp(params->s_func, "MUNK&ANDERSON") == 0)
      master->s_func = s_munk_and;
    if (strcmp(params->s_func, "EIFLER&SCHRIMPF") == 0)
      master->s_func = s_eifler_schrim;
    if (strcmp(params->s_func, "SCHUMANN&GERZ") == 0)
      master->s_func = s_schum_gerz;
  }

  /*-----------------------------------------------------------------*/
  /* If something went wrong, issue a warning and quit */
  if (init == NULL) {
    if (strlen(params->mixsc) != 0)
      hd_warn
        ("closure_init: The mixing scheme '%s' specified in the parameter file is unknown or unimplemented.\n",
         params->mixsc);
    else
      hd_warn("closure_init: No MIXING_SCHEME parameter was specified.\n");
    hd_warn("closure_init: The supported mixing schemes are :\n");
    hd_warn("closure_init: Identifier         Description\n");
    hd_warn("closure_init: ------------------------------\n");
    for (i = 0; mixing_list[i].name; ++i)
      hd_warn("closure_init: %-17s  %s\n", mixing_list[i].name,
              mixing_list[i].description);
    hd_quit
      ("closure_init: Please specify a valid vertical mixing scheme in the parameter file.\n");
  }


  if (strcasecmp(params->mixsc, "W88") == 0)
    master->fcf = params->fcf = 0.0;
  else
    master->fcf = params->fcf = 1.0;

  /*-----------------------------------------------------------------*/
  /* Call the initialisation routine for the chosen scheme */
  init(geom, master);
}

/* END closure_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set vertical mixing coefficients on all open           */
/* boundaries in a given window.                                     */
/*-------------------------------------------------------------------*/
void bdry_closure(geometry_t *window, /* Processing window */
                  window_t *windat, /* Window data structure */
                  win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int n;                        /* Counters */

  open_bdrys_t **open = window->open;

  /* Get the mean Kz if required */
  if (wincon->means & KZ_M && windat->Kzm) {
    int cc, c, cs;
    double t = windat->dtf;
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      windat->Kzm[c] = (windat->Kzm[c] * windat->meanc[cs] + 
		       windat->Kz[c] * t)  / (windat->meanc[cs] + t);
    }
  }

  /* Smooth Vz and Kz if required */
  smooth_closure(window, windat, wincon);

  for (n = 0; n < window->nobc; n++) {

    if (!(open[n]->bcond_Vz & NOTHIN)) {
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_t,
              open[n]->obc_t, open[n]->oi1_t, open[n]->oi2_t,
              open[n]->cyc_t, windat->Vz, NULL, NULL, open[n]->bcond_Vz, 
              windat->dt, NULL, NULL, 0, 0);

    }
    if (!(open[n]->bcond_Kz & NOTHIN)) {
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_t,
              open[n]->obc_t, open[n]->oi1_t, open[n]->oi2_t,
              open[n]->cyc_t, windat->Kz, NULL, NULL, open[n]->bcond_Kz, 
              windat->dt, NULL, NULL, 0, 0);
    }
  }
}

/* END bdry_closure()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* The condition Rflux < Rcrit in MY2.0 can create grid point noise  */
/* when the background coeficients are used (Rflux >= Rcrit). This   */
/* can create noise in the temperature solution if a heatflux is     */
/* imposed. Apply a Shuman filter to remove the noise in this  case. */
/* This routine may be applied to other turbulence closure schemes.  */
/* Shuman filters are described in Kowakik and Murty, p144.          */
/*-------------------------------------------------------------------*/
void smooth_closure(geometry_t *window,   /* Processing geometry     */
		    window_t *windat,     /* Window data             */
		    win_priv_t *wincon)   /* Window constants        */
{
  if (wincon->smooth_VzKz) {
    shuman_3d(window, windat->Kz, 0.5);
    shuman_3d(window, windat->Vz, 0.5);
  }
}

/* END closure_smooth()                                              */
/*-------------------------------------------------------------------*/


#define CM0		(0.5562)
#define CM0c            (CM0*CM0*CM0)
#define SQ2             (sqrt(2.0))
#define SMALL            (1e-10)

void s_null(double aN, double aM, double *cmu, double *cmu_d)
{
}

/*-------------------------------------------------------------------*/
/* Canuto et. al. (2001) model A stabilityy functions, Canuto et. al */
/* (2001) eq. 14c, 17, 20b.                                          */
/*-------------------------------------------------------------------*/
#define L1   (0.1070)
#define L2   (0.0032)
#define L3   (0.0864)
#define L4   (0.1200)
#define L5   (11.900)
#define L6   (0.4000)
#define L7   (0.0000)
#define L8   (0.4800)

#define s0   (1.5*L1*L5*L5)
#define s1   (-L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8)
#define s2   (-0.375*L1*(L6*L6-L7*L7))
#define s4   (2.*L5)
#define s5   (2.*L4)
#define s6   (2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7))
#define d0   (3.*L5*L5)
#define d1   (L5*(7.*L4+3.*L8))
#define d2   (L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7))
#define d3   (L4*(4.*L4+3.*L8))
#define d4   (L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2))
#define d5   (0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7))

#define GHmin_CA  (-0.28)
#define GH0_CA    (0.0329)
#define GHcrit_CA (0.03)
#define Cf        (4.0 / (CM0c * CM0c))

void s_canutoA(double aN, double aM, double *cmu, double *cmu_d)
{
  double d, GH, GMmin;

  /* Get max and min in terms of GH; GH = -0.5aN */
  GH = -0.5 * aN;

  /* Smooth GH near the maximum limit (Warner et. al. (2005) Eq. 33) */
  /*
  d = (GH - GHcrit_CA);
  GH = (GH - d * d) / (GH + GH0_CA - 2.0*GHcrit_CA);
  */
  /* Set minimim and maximum shear and buoyancy numbers; Warner et.  */
  /* al. 2005, Table 4 and eq. 40b (note aN=-2GH).                   */
  d = (d0/(2.0*Cf) - d1*GH + d3*2.0*Cf*GH*GH) / (d2 - d4*2.0*Cf*GH);
  aM = min(aM, 2.0*d);
  aN = -2.0 * max(GHmin_CA, min(GH, GH0_CA));

  /* Scale buoyancy number applicable for these stability functions; */
  /* Canuto uses 4tke^2/diss^2 as an argument. Note: tke^2/diss^2 =  */
  /* L^2/(CM0^6*tke), derived from Burchard et al (1998) Eqn 11.     */
  /* The factor 4 comes from the Canuto et al (2001) Eqn 14c in      */
  /* Eqn 17a (alternatively Eq 32, 37 and 39 in Warner et al (2005). */
  aN *= Cf;
  aM *= Cf;

  /* Stability function parameter - based on Canuto et al (2001)     */
  /* Eqn 17-19.                                                      */
  d =
    d0 + d1 * aN + d2 * aM + d3 * aN * aN + d4 * aN * aM +
    d5 * aM * aM;
  *cmu = SQ2 * ((s0 + s1 * aN + s2 * aM) / d);
  *cmu_d = SQ2 * ((s4 + s5 * aN + s6 * aM) / d);
}

/* END s_canutoA()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Canuto et. al. (2001) model A stabilityy functions, Canuto et. al */
/* (2001) eq. 14c, 17, 22c.                                          */
/*-------------------------------------------------------------------*/
#define Lb1  (0.1270)
#define Lb2  (0.00336)
#define Lb3  (0.0906)
#define Lb4  (0.101)
#define Lb5  (11.2)
#define Lb6  (0.4000)
#define Lb7  (0.0000)
#define Lb8  (0.3180)

#define S0   (1.5*Lb1*Lb5*Lb5)
#define S1   (-Lb4*(Lb6+Lb7)+2.*Lb4*Lb5*(Lb1-1./3.*Lb2-Lb3)+1.5*Lb1*Lb5*Lb8)
#define S2   (-0.375*Lb1*(Lb6*Lb6-Lb7*Lb7))
#define S4   (2.*Lb5)
#define S5   (2.*Lb4)
#define S6   (2./3.*Lb5*(3.*Lb3*Lb3-Lb2*Lb2)-0.5*Lb5*Lb1*(3.*Lb3-Lb2)+0.75*Lb1*(Lb6-Lb7))
#define D0   (3.*Lb5*Lb5)
#define D1   (Lb5*(7.*Lb4+3.*Lb8))
#define D2   (Lb5*Lb5*(3.*Lb3*Lb3-Lb2*Lb2)-0.75*(Lb6*Lb6-Lb7*Lb7))
#define D3   (Lb4*(4.*Lb4+3.*Lb8))
#define D4   (Lb4*(Lb2*Lb6-3.*Lb3*Lb7-Lb5*(Lb2*Lb2-Lb3*Lb3))+Lb5*Lb8*(3.*Lb3*Lb3-Lb2*Lb2))
#define D5   (0.25*(Lb2*Lb2-3*Lb3*Lb3)*(Lb6*Lb6-Lb7*Lb7))

#define GHmin_CB  (-0.28)
#define GH0_CB    (0.0444)
#define GHcrit_CB (0.0414)

void s_canutoB(double aN, double aM, double *cmu, double *cmu_d)
{
  double d, GH, GMmin;

  /* Get max and min in terms of GH; GH = -0.5aN */
  GH = -0.5 * aN;

  /* Smooth GH near the maximum limit (Warner et. al. (2005) Eq. 33) */
  /*
  d = (GH - GHcrit_CB);
  GH = (GH - d * d) / (GH + GH0_CB - 2.0*GHcrit_CB);
  */
  /* Set minimim and maximum shear and buoyancy numbers; Warner et.  */
  /* al. 2005, Table 4 and eq. 40b (note aN=-2GH).                   */
  d = (d0/(2.0*Cf) - d1*GH + d3*2.0*Cf*GH*GH) / (d2 - d4*2.0*Cf*GH);
  aM = min(aM, 2.0*d);
  aN = -2.0 * max(GHmin_CB, min(GH, GH0_CB));

  /* Scale buoyancy number applicable for these stability functions; */
  /* Canuto uses 4tke^2/diss^2 as an argument. Note: tke^2/diss^2 =  */
  /* L^2/(CM0^6*tke), derived from Burchard et al (1998) Eqn 11.     */
  /* The factor 4 comes from the Canuto et al (2001) Eqn 14c in      */
  /* Eqn 17a (alternatively Eq 32, 37 and 39 in Warner et al (2005). */
  aN *= Cf;
  aM *= Cf;

  /* Stability function parameter - based on Canuto et al (2001)     */
  /* Eqn 17-19.                                                      */
  d =
    D0 + D1 * aN + D2 * aM + D3 * aN * aN + D4 * aN * aM +
    D5 * aM * aM;
  *cmu = SQ2 * ((S0 + S1 * aN + S2 * aM) / d);
  *cmu_d = SQ2 * ((S4 + S5 * aN + S6 * aM) / d);
}

/* END s_canutoB()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Kantha and Clayson stability functions, Warner et al (2005) Eq 30 */
/* 31.                                                               */
/*-------------------------------------------------------------------*/
#define A1 (0.92)
#define B1 (16.6)
#define A2 (0.74)
#define B2 (10.1)
#define C1 (0.08)
#define C2 (0.7)
#define C3 (0.3)

#define r1 (3.0)
#define r2 (6.0)
#define r3 (9.0)
#define r4 (18.0)
#define r5 (-1.0/3.0)
#define r7 (2.0/3.0)
#define CO1 (A2*(1.0 - r2*A1/B1))
#define CO2 (r1*A2*(r2*A1+B2*(1.0-C3)))
#define CO3 (pow(B1,r5))
#define CO4 (r4*A1*A1+r3*A1*A2*(1-C2))
#define CO5 (r3*A1*A2)

#define GHmin_KC  (-0.28)
#define GH0_KC    (0.0233)
#define GHcrit_KC (0.02)

void s_kantha_clayson(double aN, double aM, double *cmu, double *cmu_d)
{
  double d, GH;

  /* Scale buoyancy number applicable for these stability functions; */
  /* uses GH = -L^2/2K as an argument (Warner 92005) Eq. 32.         */
  GH = -0.5 * aN;

  /* Smooth GH near the maximum limit (Warner et. al. (2005) Eq. 33) */
  d = (GH - GHcrit_KC);
  GH = (GH - d * d) / (GH + GH0_KC - 2.0*GHcrit_KC);

  /* Set minimim buoyancy numbers; Warner et. al. 2005, Table 3      */
  GH = max(GHmin_KC, min(GH, GH0_KC));

  *cmu_d = CO1/ (1.0 - CO2 * GH);
  *cmu = CO3 + *cmu * CO4 * GH;
  *cmu = SQ2 * (*cmu_d / (1.0 - CO5 * GH));
  *cmu_d *= SQ2;
}

/* END s_kantha_clayson()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Stability functions used in POM code; Mellor, 1992 (User's Guide) */
/* Eq. 14.2.                                                         */
/*-------------------------------------------------------------------*/
#define COEF1 (A2*(1.0-r2*A1/B1))
#define COEF2 (r1*A2*B2+r4*A1*A2)
#define COEF3 (A1*(1.0-r1*C1-r2*A1/B1))
#define COEF4 (r4*A1*A1+r3*A1*A2)
#define COEF5 (r3*A1*A2)

#define GHmin_MY  (-0.28)
#define GH0_MY    (0.0233)
#define GHcrit_MY (0.02)

/* Critical GH, Warner et al (2005) Eqn 33, Table 3                  */
#define CRIT_GH           (0.02)

void s_pom(double aN, double aM, double *cmu, double *cmu_d)
{
  double GH;

  /* Scale buoyancy number applicable for these stability functions; */
  /* uses GH = L^2/q^2 = -L^2/2k as an argument.                     */
  GH = -0.5 * aN;

  /* Maximum value of GH; Galperin et al (1988) Eqn 29, or Mellor &  */
  /* Yamada (1982), p859.                                            */
  /* Minimum value of GH; Galperin et al (1988) Eqn 30.              */
  GH = max(GHmin_MY, min(GH, GH0_MY));

  *cmu_d = COEF1 / (1.0 - COEF2 * GH);
  *cmu = COEF3 + *cmu * COEF4 * GH;
  *cmu = SQ2 * (*cmu_d / (1.0 - COEF5 * GH));
  *cmu_d *= SQ2;
}

/* END s_pom()                                                       */
/*-------------------------------------------------------------------*/

#define aN_min   (-0.0466)
#define aN_max   (0.56)

/*-------------------------------------------------------------------*/
/* Galperin et. al. (1988) stability functions, e.g. Burchard et. al */
/* (1998), Eq. 29, 30.                                               */
/*-------------------------------------------------------------------*/
void s_galperin(double aN, double aM, double *cmu, double *cmu_d)
{
  /* Limit the buoyancy number, Burchard et. al. (1998), p10,547     */
  aN = max(min(aN_max, aN), aN_min);

  *cmu = (CM0 + 2.182 * aN) / (1 + 20.4 * aN + 53.12 * aN * aN);
  *cmu_d = 0.6985 / (1.0 + 17.34 * aN);
}

/* END s_galperin()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Munk and Anderson (1948) stability functions, e.g. Burchard et al */
/* (1999). GOTM, a general ocean turbulence model. Theory,           */
/* implementation and test cases. Eq. 2.73 and 2.74.                 */
/*-------------------------------------------------------------------*/
#define Pr0    (0.74)      /* used in GOTM code */
#define Pr0a    (0.7143)
void s_munk_and(double aN, double aM, double *cmu, double *cmu_d)
{
  double Ri, Pr;

  /* Limit the buoyancy number, Burchard et. al. (1998), p10,547     */
  aN = max(min(aN_max, aN), aN_min);
  Ri = aN / (aM + SMALL);

  Pr = (Ri >= SMALL) ? Pr0 * pow(1.+3.33*Ri,1.5)/sqrt(1.+10.0*Ri) : Pr0;
  *cmu = CM0;
  *cmu_d = CM0 / Pr;
}

/* END s_munk_and()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Eifler and Schrimpf (1992) stablity functions, e.g. Burchard et   */
/* al (1999). GOTM, a general ocean turbulence model. Theory,        */
/* implementation and test cases. Eq. 2.76 and 2.77.                 */
/*-------------------------------------------------------------------*/
void s_eifler_schrim(double aN, double aM, double *cmu, double *cmu_d)
{
  double Ri, Rf, Pr;

  /* Limit the buoyancy number, Burchard et. al. (1998), p10,547     */
  aN = max(min(aN_max, aN), aN_min);
  Ri = 0.5 / Pr0 * aN / (aM + SMALL);

  Rf = sqrt(Ri * Ri + 1.0) - Ri;
  Rf = max(Ri*Ri, 2.0);
  Pr = max(min(1.0 /Pr0 * sqrt(Rf), 0.18), 2.0);

  *cmu = CM0;
  *cmu_d = CM0 * Pr;
}

/* END s_eifler_schrim()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Schumann and Gerz (1995) stablity functions.                      */
/*-------------------------------------------------------------------*/
void s_schum_gerz(double aN, double aM, double *cmu, double *cmu_d)
{
  double Ri, Pr;
  double limit = 3;

  /* Limit the buoyancy number, Burchard et. al. (1998), p10,547     */
  aN = max(min(aN_max, aN), aN_min);
  Ri = aN / (aM + SMALL);
  Pr = (Ri >= SMALL) ? Pr0 * exp(-Ri/(Pr0*0.25))+Ri/0.25 : Pr0;

  *cmu = CM0;
  *cmu_d = CM0 / min(limit, Pr);
}

/* END s_shum_gerz()                                                 */
/*-------------------------------------------------------------------*/
