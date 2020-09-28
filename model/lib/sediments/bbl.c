/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/sediment/bbl.c
 *  
 *  Description:
 *  Calculate bottom stresses.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: bbl.c 5955 2018-09-17 00:23:31Z mar644 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"

#define VK 0.4
#define EPS 0.00001

/* function prototypes */
static int madsen94(sed_column_t *sm, double ubr, double wr, double ucr,
                     double zr, double phiwc, double zo, int iverbose,
                     double *pustrc, double *pustrwm,
                     double *pustrr, double *pfwc, double *pzoa);

static int wripples(sed_column_t *sm, double ubr, double wr,
                     double zo);
static int swart74(sed_column_t *sm, double ubr, double wr, double ucr,
                    double zr, double zo, double zoc, double *pustrc, 
		    double *pustrwm,
                    double *pustrr, double *pfwc, double *pzoa);
static void reference_velocity(sed_column_t *sm, double *pzr, double *pui,
                               double *puj, double *pzo);
static void calc_phys_roughness(sed_column_t *sm);

void reef_scale_depth(sediment_t *sediment, sed_column_t *sm);

/* Functions for calculating enhanced bottom friction.
 * The madsen routines where kindly provided by Chris Sherwood.
 */

  void bbl(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;

    /* Compute the enchanced bottom friction coefficient. */
    double ub;
    double period;
    double dir;
    double u;
    double udir;
    double ui;
    double uj;
    double uval;
    double vval;
    double ustrc;
    double ustrwm;
    double ustrr;
    double fwc;
    double zoa;
    double zr;


    /* Waves: Get ub, period and direction */
    ub = max (sm->wave_ub, 0);
    period = sm->wave_period;
    if (period < 1e-21) period = 8.; // if period is not specified, set default 
    dir = sm->wave_dir;

   /* Ripples and physical roughness*/
    if (param->calc_ripples) {
      if (sm->coh_sed[sm->topk_sed] < 10.) {
	/*sandy seabed*/
        if (wripples(sm, ub, 2 * PI / period, sm->z0_skin))
	  i_set_error(sediment->hmodel, sm->col_number, LFATAL, "bbl:wripples: Error encountered in wripples().\n");
      }
      else {
        sm->hripples = param->bioriph;
        sm->lripples = param->bioripl;
      }
    } else {
      double physriph = param->physriph;
      double physripl = param->physripl;
      if (param->physriph_spv != NULL)
	physriph = param->physriph_spv[sm->col_number-1];
      if (param->physripl_spv != NULL)
	physripl = param->physripl_spv[sm->col_number-1];
      
      sm->hripples = max(param->bioriph, physriph);
      sm->lripples = max(param->bioripl, physripl);
    }

    /* Calculate total (grain + ripples) physical roughness */
    calc_phys_roughness(sm);

    /* Currents: magnitude and direction of the velocity 
       at the cell centre. */
    reference_velocity(sm, &zr, &ui, &uj, &sm->z0_phys);
    // velocity referenced to the earth coordinates (ie W-E, S-N axes)
    // theta is cell-centered angle between grid and earth coordinates 
    uval = ui * cos(sm->theta) - uj * sin(sm->theta);
    vval = ui * sin(sm->theta) + uj * cos(sm->theta);
    if ((uval == 0.0) && (vval == 0.0)) {
      u = 0.0;
      udir = dir;
    } else {
      u = sqrt(uval * uval + vval * vval);
      //      udir = fmod(450.0 - atan2(vval, uval) * 180.0 / PI, 180.0);
      // direction angle: from W-E counterclockwise
      //(note that the wave direction above is the direction to destination)
      // if the wave dir is the direction to origin, change atan2 sign)
      // if waves are referenced to S-N, add 90 to 360)
      udir = fmod(360.0 + atan2(vval, uval) * 180.0 / PI, 360.0);

    }

    /* Bottom friction and aparent roughness */
    if (param->bbl_nonlinear) {
      if (madsen94(sm, ub, 2 * PI / period, u, zr, dir - udir,
		   sm->z0_phys, 0, &ustrc, &ustrwm, &ustrr, &fwc, &zoa))
	i_set_error(sediment->hmodel, sm->col_number, LFATAL, "bbl:madsen94: Error encountered in madsen94().\n");
    }
    else {
      if (swart74(sm, ub, 2 * PI / period, u, zr, sm->z0_skin,
		  sm->z0_phys, &ustrc, &ustrwm, &ustrr, &fwc, &zoa))
	i_set_error(sediment->hmodel, sm->col_number, LFATAL, "bbl:swart74: Error encountered in swart74().\n");
    }
    /* Store total bottom friction velocity */
      sm->ustrcw = ustrr;
      sm->ustrcw_skin = ustrr;
}

/********************************************************/
static void reference_velocity(sed_column_t *sm, double *pzr, double *pui,
                               double *puj, double *pz0)
{
  int kb = sm->botk_wc;
  int kt = sm->topk_wc;
  double bot = sm->botz_wc;
  double top = sm->topz_wc;
  double zc;
  int kzc;
  /* Reference height */
  if (kb < kt) 
     zc = (sm->gridz_wc[kb + 1] - bot) / 2.;
  else
    zc = (top - bot) / 2.;
  kzc=kb;

  /*Fix for 3d-z-grid atrefacts: near-bottom cells with stagnant water, 
   occassionally produced by a 3d-z-grid over complex bathymetry, result 
   in zero friction and artificial accumulation of sediments in these cells. 
   To mitigate such grid-dependency, Mecosed takes velocity from the next 
   cell above the near bottom stagnant one. 
  */
  // u2bcc =  sm->u1_wc[kb]*sm->u1_wc[kb]+sm->u2_wc[kb]*sm->u2_wc[kb];
  // if(fabs(u2bcc) < 1e-9 && kb < kt) {

  // tmp fix: taking all ref velocity from the 2nd cell above the ground
  // comment out if-block below to reference velocities in the the near-bottom cell
  if(kb < kt) {
     zc = (sm->cellz_wc[kb + 1] - bot);
     kzc=kb+1;
  }

  /* If the reference height is less than 1m use log-profile to 
    interpolate velocities to 1m hight, so that in a typical application 
    zr exceeds the thickness of wave bbl */
    *pzr = zc;
    if (zc < 1.) *pzr = 1.;

  /*check that the cell centre is above the roughness hight */
  if (zc < *pz0)  zc = (*pz0) * 2;

  /* Reference velocity */
  *pui = sm->u1_wc[kzc] * log(*pzr / *pz0) / log(zc / *pz0);
  *puj = sm->u2_wc[kzc] * log(*pzr / *pz0) / log(zc / *pz0);

  return;
}

/** Wave-current friction factor.
  * (Equations 32 and 33 in Madsen, 1994).
  *
  * @param cmu Relative strength of currents from Eqn 27
  * @param cukw Relative roughness = cmu*ubr/(kn*wr)
  * @return Wave-current friction factor
  *
  * @author Chris Sherwood, CSIRO
  * @version 11 July 1997
  */
static double fwc94(double cmu, double cukw)
{
  if (cukw <= 100.)
    return (cmu * exp(7.02 * pow(cukw, -0.078) - 8.82));
  else                          /* ( cukw > 100. ) */
    return (cmu * exp(5.61 * pow(cukw, -0.109) - 7.30));
}

/** Grant-Madsen wave-current calculations from Madsen (1994).
 * 
 * This version supercedes madsen94 (it is better optimized) but the
 * argument list is differant.
 * 
 *  (Note that this formulation does not address angular difference
 *  between stress and velocity).
 *
 * @param ci Cell i index.
 * @param cj Cell j index.
 * @param ub Wave-orbital velocity outside wbl [m/s]
 * @param T Wave period [s]
 * @param ur Current velocity at height zr [m/s]
 * @param zr Reference height for current velocity [m]
 * @param phiwc Angle between currents and waves at zr (degrees)
 * @param zo Bottom roughness height [m]
 * @param iverbose Switch; when 1, extra output
 * @param pustrc Current shear velocity u*c [m/s]
 * @param pustrr Wave shear velocity u*w [m/s]
 * @param pphicw Angle between u*cw and waves at bed (degrees)
 * @param pustrcw Wave-current combined shear velocity u*cw [m/s]
 * @param pfwc Wave friction factor [ ]
 * @param pzoa Apparent bottom roughness [m]
 *
 * @author Chris Sherwood, CSIRO
 * @version 5 June 1997
 * MH 07/2012. Changed function type to int, returning 1 on fail.
 */
static int madsen94(sed_column_t *sm, double ubr, double wr, double ucr,
                     double zr, double phiwc, double zo, int iverbose,
                     double *pustrc, double *pustrwm,
                     double *pustrr, double *pfwc, double *pzoa)
{
#define MAXIT 20
#define LOWU (0.01)
  double rmu[MAXIT], Cmu[MAXIT], fwc[MAXIT], dwc[MAXIT];
  double ustrwm2[MAXIT], ustrr2[MAXIT], ustrc[MAXIT];
  double kN, cosphiwc, ustrr, lnzr, lndw, lnln, bigsqr, diff, zoa;
  int i, nit, errflg;
  int ci = sm->i;
  int cj = sm->j;

  kN = 30. * zo;
  zoa = zo;

  /* junk return values */
  *pustrc = 0.0;
  *pustrwm = 0.0;
  *pustrr = 0.0;
  *pfwc = .4;
  *pzoa = zoa;

  /* some data checks */
  if (wr <= 0.) {
    sedtag(LWARN,"sed:bbl:madsen94", "Bad value for ang. freq. in Madsen94 at (%d,%d): wr = %g", ci, cj, wr); 
    return(1);
  }
  else if (ubr < 0.) {
      sedtag(LWARN,"sed:bbl:madsen94", "Bad value for orbital vel. in Madsen94 at (%d,%d): ub = %g",ci, cj, ubr); 
      return(1);
  }
  else if (kN <= 0.) {
      sedtag(LWARN,"sed:bbl:madsen94", "Negative roughness in Madsen94 at (%d,%d): kN = %g", ci, cj, kN); 
      return(1);
  }

/* The reference height for current velocity zr should be located above the turbulent wave boundary layer. Here following Madsen 1994, we use a somewhat arbitrary criteria to check this condition (eg assuming that the wave bbl thickness ~ kN). If zr is too small then go to a linear bbl model swart74*/
  if (zr < 5. * kN) {
      if (zr < 5. * zo) { 
	/* The reference height zr should exceed the physical roughness zo*/
	  sedtag(LWARN,"sed:bbl:madsen94", "Low value for ref. level in Madsen94 at (%d,%d): zr = %g", ci, cj, zr); 
	  return(1);
      }
      if(swart74(sm,ubr,wr,ucr,zr,zo,zo,pustrc,pustrwm,pustrr,pfwc,pzoa))
	return(1);
      else
	return(0);
  }
  
  if (ubr <= LOWU) {
    if (ucr <= LOWU) {
      *pustrc = 0.;
      *pustrwm = 0.;
      *pustrr = 0.;
      *pzoa = zo;
      return(0);
    }
    ustrc[0] = ucr * VK / log(zr / zo);
    *pustrc = ustrc[0];
    *pustrwm = 0.;
    *pustrr = ustrc[0];
    *pzoa = zo;
    return(0);
  }

  cosphiwc = fabs(cos(phiwc * (PI / 180.)));
  rmu[0] = 0.;
  Cmu[0] = 1.;
  fwc[0] = fwc94(Cmu[0], (Cmu[0] * ubr / (kN * wr))); /* Eqn. 32 or 33 */
  ustrwm2[0] = 0.5 * fwc[0] * ubr * ubr;  /* Eqn. 29 */
  ustrr2[0] = Cmu[0] * ustrwm2[0];  /* Eqn. 26 */
  ustrr = sqrt(ustrr2[0]);
  dwc[0] = (Cmu[0] * ubr / (kN * wr)) >= 8. ? 2. * VK * ustrr / wr : kN;
  lnzr = log(zr / dwc[0]);
  lndw = log(dwc[0] / zo);
  lnln = lnzr / lndw;
  bigsqr =
    (-1. + sqrt(1 + ((4. * VK * lndw) / (lnzr * lnzr)) * ucr / ustrr));
  ustrc[0] = 0.5 * ustrr * lnln * bigsqr;
  errflg = 0;
  for (i = 1, nit = 1, diff = 1.; diff > 0.0005 && i < MAXIT; i++, nit++) {
    rmu[i] = ustrc[i - 1] * ustrc[i - 1] / ustrwm2[i - 1];  /* Eqn. 28 */
    Cmu[i] = sqrt(1. + 2. * rmu[i] * cosphiwc + rmu[i] * rmu[i]); /* Eqn
                                                                     27 */
    fwc[i] = fwc94(Cmu[i], (Cmu[i] * ubr / (kN * wr))); /* Eqn. 32 or 33 */
    ustrwm2[i] = 0.5 * fwc[i] * ubr * ubr;  /* Eqn. 29 */
    ustrr2[i] = Cmu[i] * ustrwm2[i];  /* Eqn. 26 */
    ustrr = sqrt(ustrr2[i]);
    dwc[i] = (Cmu[i] * ubr / (kN * wr)) >= 8. ? 2. * VK * ustrr / wr : kN;  /* Eqn. 
                                                                               36 
                                                                             */
    if (dwc[i] > 0.8 * zr) {
      errflg = 1;
      dwc[i] = 0.8 * zr;        /* clumsy fix for ill posed cases */
    }
    lnzr = log(zr / dwc[i]);
    lndw = log(dwc[i] / zo);
    lnln = lnzr / lndw;
    bigsqr =
      (-1. + sqrt(1 + ((4. * VK * lndw) / (lnzr * lnzr)) * ucr / ustrr));
    ustrc[i] = 0.5 * ustrr * lnln * bigsqr; /* Eqn. 38 */
    diff = fabs((fwc[i] - fwc[i - 1]) / fwc[i]);
  }
  if (errflg)
      sedtag(LWARN,"sed:bbl:madsen94", " dwc > 0.8*zr");

  *pustrwm = sqrt(ustrwm2[nit - 1]);
  *pustrc = ustrc[nit - 1];
  *pustrr = ustrr;
  *pfwc = fwc[nit - 1];
  zoa =
    exp(log(dwc[nit - 1]) -
        (ustrc[nit - 1] / ustrr) * log(dwc[nit - 1] / zo));
  *pzoa = zoa;
#ifdef HAVE_EMSLIB_MODULE
  if (is_log_enabled(LTRACE)) {
#endif
    for (i = 0; i < nit; i++) {
      sedtag(LTRACE,"sed:bbl:madsen94", " i=%2d rmu=%9.6f Cmu=%9.6f fwc=%9.6f dwc=%9.6f u*c=%9.4f u*wm=%9.4f u*r=%9.4f", i, rmu[i], Cmu[i], fwc[i], dwc[i], ustrc[i], sqrt(ustrwm2[i]), sqrt(ustrr2[i]));
    }
#ifdef HAVE_EMSLIB_MODULE
  }
#endif
  return(0);
}

/********************************************************/


/********************************************************/
/* MH 07/2012. Changed function type to int, returning 1 on fail. */

static int swart74(sed_column_t *sm, double ubr, double wr, double ucr,
                    double zr, double zo, double zoc, double *pustrc, 
		    double *pustrwm,
                    double *pustrr, double *pfwc, double *pzoa)
{
  double fwc;
  int i = sm->i;
  int j = sm->j;
  double kN = 30. * zo;

  /* junk return values */
  *pustrc = 0.0;
  *pustrwm = 0.0;
  *pustrr = 0.0;
  *pfwc = .4;
  *pzoa = zo;

  /* some data checks */
  if (wr <= 0.) {
    sedtag(LWARN,"sed:bbl:swart74", "Bad value for ang. freq. in swart74 at (%d,%d): wr = %g", i, j, wr); 
    return(1);
  }
  else if (ubr < 0.) {
    sedtag(LWARN,"sed:bbl:swart74", "Bad value for orbital vel. in swart74 at (%d,%d): ub = %g",i, j, ubr); 
    return(1);
  }
  else if (kN <= 0.) {
    fprintf(stderr,"swart74: Negative roughness in swart94 at (%d,%d): kN = %g", i, j, kN); 
    return(1);
  }

  if (zr < 5. * zoc) {
    sedtag(LWARN,"sed:bbl:swart74", "Low value for ref. level in Madsen94 at (%d,%d): zr = %g", i, j, zr); 
    return(1);
  }

  if (ubr <= LOWU) {
    if (ucr <= LOWU) {
      *pustrc = 0.;
      *pustrwm = 0.;
      *pustrr = 0.;
      return(0);
    }
    *pustrc = ucr * VK / log(zr / zoc);
    *pustrwm = 0.;
    *pustrr = ucr * VK / log(zr / zoc);
    return(0);
  }

  /* Wave Shear Velocity */
  fwc = exp(5.213 * pow(kN * wr / ubr, 0.194) - 5.977); /* wave friction
                                                           coeff */
  *pfwc = fwc;
  *pustrwm = sqrt(0.5 * fwc * ubr * ubr);
  /* Current Shear Velocity */
  *pustrc = ucr * VK / log(zr / zoc);
  /* Wave-Current Shear Velocity */
  *pustrr = (*pustrwm) + (*pustrc);
  return(0);
}

/********************************************************/

/* Black K., Oldman J. 1999, Marine Geology 162, pp 121-132. */
static int wripples(sed_column_t *sm, double ubr, double wr,
                     double zo)
{
  double teta;
  double hrip = 0.;
  double lrip = 0.; 
  double z;
  int i = sm->i;
  int j = sm->j;
  double kN = 30. * zo;
  double fwc;

  /* some data checks */
  if (wr <= 0.) {
    sedtag(LWARN,"sed:bbl:wripples", "Bad value for ang. freq. in wripples at (%d,%d): wr = %g",i, j, wr); 
    return(1);
  } else if (ubr < 0.) {
    sedtag(LWARN,"sed:bbl:wripples", "Bad value for orbital vel. in wripples at (%d,%d): ub = %g",i, j, ubr);
    return(1);
  }

  /* Start Ripples */
  fwc = exp(5.213 * pow(kN * wr / ubr, 0.194) - 5.977);
  teta = 0.5 * fwc * ubr * ubr / ((2.65 - 1.) * 9.8 * zo);
  z = (4. * teta * 1.e-6) / sqrt((2.65 - 1.) * 9.8 * zo * zo * zo);

  if (z <= 0.0016) {
    hrip = sm->hripples;
    lrip = sm->lripples;
  } else if (0.0016 < z && z < 0.012) {
    hrip = 0.018 * (ubr / wr) * pow(z, -0.50);
    lrip = 0.120 * (ubr / wr) * pow(z, -0.49);
  } else if (0.012 <= z && z < 0.18) {
    hrip = 0.0007 * (ubr / wr) * pow(z, -1.23);
    lrip = 0.0667 * (ubr / wr) * pow(z, -0.58);
  } else if (0.18 <= z) {
    hrip = 0.0;
    lrip = 0.0;
  }

  sm->hripples = hrip;
  sm->lripples = lrip;

  return(0);
}

/*************************/
static void calc_phys_roughness(sed_column_t *sm)
{
  double hrip = sm->hripples;
  double lrip = sm->lripples;
  double z0_skin = sm->z0_skin;

  /* Ripple Roughness */
  if (lrip > 0.001)
    sm->z0_phys = z0_skin + (28. * hrip * hrip / lrip) / 30.;
  else
    sm->z0_phys = z0_skin;

  return;
}


#ifdef __cplusplus
}
#endif
