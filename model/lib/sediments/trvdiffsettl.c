/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/sediment/trvdiffsettl.c
 *  
 *  Description:
 *  Calculate vertical diffusion and settling of tracer
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: trvdiffsettl.c 5955 2018-09-17 00:23:31Z mar644 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"

static int implicit_vadv_vdiff(double dt, int nz, double *dz, double *dzold,
                                double *por, double *porold, double *var,
                                double *dvar, double *Kz, double *w,
                                double botinflux, double topinflux,
                                double botoutfluxcoef,
                                double topoutfluxcoef, int kb, int kt,
                                double *Splus, double *Sminus, double *dzface);

void vdiff_sedim_wc(sediment_t *sediment, sed_column_t *sm,
		    sed_tracer_t *tracer,
                    sed_tracer_t *tracercar, double botinflux,
                    double topinflux, double botoutfluxcoef,
                    double topoutfluxcoef, double *dvar)
{
  sed_params_t *param = sediment->msparam;
  long k, kt, kb;
  double top, bot;
  double *Kzij;                 /* diffusion coefficients for 1 column */
  double *w;
  double *var = sm->tr_wc[tracer->n];
  double s = 1.;                /* Kz scale factor */
  double *dummy;
  int ncar = tracercar->n;

  /* Allocate temporary arrays */
  Kzij = d_alloc_1d(param->nz + 1);
  w = d_alloc_1d(param->nz + 1);
  dummy = d_alloc_1d(param->nz + 1);
  kb = sm->botk_wc;
  kt = sm->topk_wc;
  top = sm->topz_wc;
  bot = sm->botz_wc;
  /* nothing to do if dry */
  /* if( top <= bot ) */
  /* continue; */

  if ( tracer->diffuse == 0 )
   s=0;

  /* Get diffusion coefficients for this column */
  Kzij[kt + 1] = 0;
  Kzij[kt] = s * sm->Kz_wc[kt];
  for (k = kt - 1; k > kb; k--)
    Kzij[k] = s * sm->Kz_wc[k];
  Kzij[kb] = 0.;
  /* sediment velocity */
    for (k = kb; k <= kt + 1; k++)
      w[k] = -sm->gridvel_wc[k] + sm->svel_wc[ncar][k];

  for (k = kb; k <= kt; k++)
    dummy[k] = 1.;

  /* MH 07/2012: included instability handling */
  if (implicit_vadv_vdiff(param->dt, param->nz, sm->dz_wc, sm->dzold_wc,
			  dummy, dummy, var, dvar, Kzij, w,
			  botinflux, topinflux, botoutfluxcoef, topoutfluxcoef,
			  kb, kt, NULL, NULL,sm->dzface_wc))
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "trvdiffsettl:vdiff_sedim_wc: Error encountered in implicit_vadv_vdiff().\n");

  /* free temporary arrays */
  d_free_1d(Kzij);
  d_free_1d(w);
  d_free_1d(dummy);
}

/******************************************************************************/

void vdiff_dissolved_wc(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
			double botinflux, double topinflux,
			double botoutfluxcoef, double topoutfluxcoef,
			double *dvar)
{
  sed_params_t *param = sediment->msparam;
  long k, kt, kb;
  double top, bot;
  double *Kzij;                 /* diffusion coefficients for 1 column */
  double *w;
  double *var = sm->tr_wc[tracer->n];
  double s = 1.;                /* Kz scale factor */
  double *dummy;
  /* Allocate temporary arrays */
  Kzij = d_alloc_1d(param->nz + 1);
  w = d_alloc_1d(param->nz + 1);
  dummy = d_alloc_1d(param->nz + 1);
  kb = sm->botk_wc;
  kt = sm->topk_wc;
  top = sm->topz_wc;
  bot = sm->botz_wc;
  /* nothing to do if dry */
  /* if( top <= bot ) */
  /* continue; */
  if ( tracer->diffuse == 0 )
   s=0;

  /* Get diffusion coefficients for this column */
  Kzij[kt + 1] = 0;
  Kzij[kt] = s * sm->Kz_wc[kt];
  for (k = kt - 1; k > kb; k--)
    Kzij[k] = s * sm->Kz_wc[k];
  Kzij[kb] = 0.;
  for (k = kb; k <= kt + 1; k++)
    w[k] = sm->watvel_wc[k];
  for (k = kb; k <= kt; k++)
    dummy[k] = 1.;
  /* MH 07/2012: included instability handling */
  if (implicit_vadv_vdiff(param->dt, param->nz, sm->dz_wc, sm->dzold_wc,
			  sm->por_wc, sm->porold_wc, var, dvar, Kzij, w,
			  botinflux, topinflux, botoutfluxcoef, topoutfluxcoef,
			  kb, kt, NULL, NULL,sm->dzface_wc))
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "trvdiffsettl:vdiff_dissolved_wc: Error encountered in implicit_vadv_vdiff().\n");

  /* free temporary arrays */
  d_free_1d(Kzij);
  d_free_1d(w);
  d_free_1d(dummy);
}

/******************************************************************************/
void vdiff_sedim_sed(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                     double botinflux, double topinflux,
                     double botoutfluxcoef, double topoutfluxcoef,
                     double *dvar)
{
  sed_params_t *param = sediment->msparam;
  long k, kt, kb;
  double top, bot;
  double *Kzij;                 /* diffusion coefficients for 1 column */
  double *var;
  double *w;
  double s = 1.;                /* Kz scale factor */
  int n = tracer->n;
  double *dummy;
  /* Allocate temporary arrays */
  Kzij = d_alloc_1d(param->sednz + 1);
  w = d_alloc_1d(param->sednz + 1);
  var = d_alloc_1d(param->sednz + 1);
  dummy = d_alloc_1d(param->sednz + 1);
  kb = sm->botk_sed;
  kt = sm->topk_sed;
  top = sm->topz_sed;
  bot = sm->botz_sed;
  /* nothing to do if dry */
  /* if( top <= bot ) */
  /* continue; */
  /* Get diffusion coefficients and vertical velocities for this column */
  Kzij[kt + 1] = 0;
  Kzij[kt] = s * sm->partic_kz[kt];
  for (k = kt - 1; k > kb; k--)
    Kzij[k] = s * sm->partic_kz[k];
  Kzij[kb] = 0.;
  /* relative velocity of bottom sediment */
  w[kt + 1] = 0.;
  for (k = kb; k <= kt; k++)
    w[k] = -sm->gridvel_sed[k] + sm->svel_consolid[k];
  w[kb] = 0.;
  for (k = kb; k <= kt; k++)
    var[k] = sm->tr_sed[n][k];
  for (k = kb; k <= kt; k++)
    dummy[k] = 1.;
  /* MH 07/2012: included instability handling */
  if (implicit_vadv_vdiff(param->dt, param->sednz, sm->dz_sed, sm->dzold_sed,
			  dummy, dummy, var, dvar, Kzij, w,
			  botinflux, topinflux, botoutfluxcoef, topoutfluxcoef,
			  kb, kt, NULL, NULL,sm->dzface_sed))
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "trvdiffsettl:vdiff_sedim_sed: Error encountered in implicit_vadv_vdiff().\n");
  /* free temporary arrays */
  d_free_1d(Kzij);
  d_free_1d(w);
  d_free_1d(var);
  d_free_1d(dummy);
}

/******************************************************************************/
void vdiff_dissolved_sed(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                         double botinflux, double topinflux,
                         double botoutfluxcoef, double topoutfluxcoef,
                         double *dvar)
{
  sed_params_t *param = sediment->msparam;
  long k, kt, kb;
  double top, bot;
  double *Kzij;                 /* diffusion coefficients for 1 column */
  double *var;
  double *w;
  double s = 1.;                /* Kz scale factor */
  int n = tracer->n;
  /* Allocate temporary arrays */
  Kzij = d_alloc_1d(param->sednz + 1);
  w = d_alloc_1d(param->sednz + 1);
  var = d_alloc_1d(param->sednz + 1);
  kb = sm->botk_sed;
  kt = sm->topk_sed;
  top = sm->topz_sed;
  bot = sm->botz_sed;
  /* nothing to do if dry */
  /* if( top <= bot ) */
  /* continue; */
  /* Get diffusion coefficients and vertical velocities for this column */
  Kzij[kt + 1] = 0;
  Kzij[kt] = s * sm->dissol_kz[kt];
  for (k = kt - 1; k > kb; k--)
    Kzij[k] = s * sm->dissol_kz[k];
  Kzij[kb] = 0.;
  for (k = kb; k <= kt; k++)
    var[k] = sm->tr_sed[n][k];
  /* relative velocity of bed water */
  for (k = kb; k <= kt + 1; k++)
    w[k] = sm->watvel_sed[k];
  /* MH 07/2012: included instability handling */
  if (implicit_vadv_vdiff(param->dt, param->sednz, sm->dz_sed, sm->dzold_sed,
			  sm->por_sed, sm->porold_sed, var, dvar, Kzij, w,
			  botinflux, topinflux, botoutfluxcoef, topoutfluxcoef,
			  kb, kt, NULL, NULL,sm->dzface_sed))
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "trvdiffsettl:vdiff_dissolved_sed: Error encountered in implicit_vadv_vdiff().\n");
  /* free temporary arrays */
  d_free_1d(Kzij);
  d_free_1d(w);
  d_free_1d(var);

}

/*********************************************************************

    File:           FIvdiffsettl.c

    Created:        Wed Jul 14 18:03:11 EST 1993

    Author:

    Purpose:        A general routine to implement fully
		    implicit vertical mixing for a single
		    column (i,j) of any of
		    the model 3d variables. The calling
		    routine should itself deal with trivial
		    cases where there is only 1 wet layer.

    Arguments:      var    - pointer to the 1d array of input values
		    dvar   - output 1d array of changes (ie new value
			     is var+dvar)
		    Kz     - 1d array of mixing coefficients at
			     layer interfaces
		    fb     - flux at bottom (+ve up)
		    ft     - flux at water surface (+ve up)
		    kb     - bottom k index
		    kt     - top k index
		    Splus  - Positive part of source term (see below)
		    Sminus - Negative part of source term (see below)

    Returns:        void

    Revisions:
		    11/01/1999 SJW
		    Temporary arrays are now locally allocated.
		    Added parameters and code to allow sources/sinks.
		    Source linearisation and positivity treatment
		    according to Patankar 1980. Source is make up of
		    two parts:

		        S = Splus - Sminus

		    where Splus and Sminus are both positive.

                    07/2012 MH
                    Made the function type int, returning 1 if it fails.

    $Id: trvdiffsettl.c 5955 2018-09-17 00:23:31Z mar644 $

*********************************************************************/

static int implicit_vadv_vdiff(double dt, int nz, double *dz, double *dzold, double *por, double *porold, double *var,  /* Array of variable values */
  double *dvar, /* Array of changes in values */
  double *Kz, /* Diffusivity values */
  double *w, 
  double botinflux,  /* Flux out of bottom */
  double topinflux, /* Flux out of top */
  double botoutfluxcoef,  /* Flux out of bottom */
  double topoutfluxcoef,  /* Flux out of top  */
  int    kb, 
  int    kt, /* Vertical indices for water bottom and top */
  double *Splus,  /* Positive part of source term */
  double *Sminus,  /* Negative part of source term */
  double *dzface  /*UR added, calculate the dface ones, not for every tracer */
  )
{
  int k = 0;
  double dzdt, dzdtold;
  double div;

#if HAVE_ALLOCA
  double *Cm1 = (double *)alloca((nz + 1) * sizeof(double));
  double *C = (double *)alloca((nz + 1) * sizeof(double));
  double *Cp1 = (double *)alloca((nz + 1) * sizeof(double));
  double *rhs = (double *)alloca((nz + 1) * sizeof(double));
  double *sol = (double *)alloca((nz + 1) * sizeof(double));
  double *ud = (double *)alloca((nz + 1) * sizeof(double));
/*  double *dzface = (double *)alloca((nz + 1) * sizeof(double)); */
#else
  double *Cm1 = d_alloc_1d(nz + 1);
  double *C = d_alloc_1d(nz + 1);
  double *Cp1 = d_alloc_1d(nz + 1);
  double *rhs = d_alloc_1d(nz + 1);
  double *sol = d_alloc_1d(nz + 1);
  double *ud = d_alloc_1d(nz + 1);
/*  double *dzface = d_alloc_1d(nz + 1); */
#endif

  /* Single layer - the calling routine should deal with this */
  if (kt <= kb) {
    sedtag(LWARN,"sed:trvdiffsettl:implicit_vadv_vdiff"," less than 2 layers\n");
    return(1);
  }
  /* Set up tri-diagonal set of equations */
  /* Bottom layer */
  dzdt = por[kb] * dz[kb] / dt;
  dzdtold = porold[kb] * dzold[kb] / dt;

  Cm1[kb] = 0.0;
  Cp1[kb] = -Kz[kb + 1] / dzface[kb + 1];
  C[kb] = dzdt - Cp1[kb];
  rhs[kb] = dzdtold * var[kb] + botinflux;

  if (w[kb + 1] < 0.)
    Cp1[kb] += w[kb + 1];
  else
    C[kb] += w[kb + 1];

  C[kb] -= botoutfluxcoef;

  /* Mid-water layers */
  for (k = kb + 1; k < kt; k++) {
    dzdt = por[k] * dz[k] / dt;
    dzdtold = porold[k] * dzold[k] / dt;

    Cm1[k] = -Kz[k] / dzface[k];
    Cp1[k] = -Kz[k + 1] / dzface[k + 1];
    C[k] = dzdt - Cm1[k] - Cp1[k];
    rhs[k] = dzdtold * var[k];

    if (w[k + 1] < 0.)
      Cp1[k] += w[k + 1];
    else
      C[k] += w[k + 1];

    if (w[k] < 0.)
      C[k] -= w[k];
    else
      Cm1[k] -= w[k];
  }
  /* Surface layer */
  dzdt = por[kt] * dz[kt] / dt;
  dzdtold = porold[kt] * dzold[kt] / dt;

  Cm1[kt] = -Kz[kt] / dzface[kt];
  Cp1[kt] = 0.0;
  C[kt] = dzdt - Cm1[kt];
  rhs[kt] = dzdtold * var[kt] - topinflux;

  if (w[kt] < 0.)
    C[kt] -= w[kt];
  else
    Cm1[kt] -= w[kt];

  C[kt] += topoutfluxcoef;

  /* Add positive part of source terms, if specified */
  if (Splus)
    for (k = kb; k <= kt; k++)
      rhs[k] += dz[k] * Splus[k];
  /* Add negative part of source terms, if specified */
  if (Sminus)
    for (k = kb; k <= kt; k++)
      C[k] += dz[k] * Sminus[k] / var[k];

  /* Solve tridiagonal system */
  div = C[kb];
  sol[kb] = rhs[kb] / div;
  for (k = kb + 1; k <= kt; k++) {
    ud[k] = Cp1[k - 1] / div;
    div = C[k] - Cm1[k] * ud[k];
    if (div == 0.0) {
      sedtag(LWARN,"sed:trvdiffsettl:implicit_vadv_vdiff","(solvetri): zero divisor\n");
      return(1);
    }
    sol[k] = (rhs[k] - Cm1[k] * sol[k - 1]) / div;
  }
  dvar[kt] = sol[kt] - var[kt];
  for (k = kt - 1; k >= kb; k--) {
    sol[k] -= ud[k + 1] * sol[k + 1];
    dvar[k] = sol[k] - var[k];
  }

#if !HAVE_ALLOCA
  /* Free temporary storage */
  d_free_1d(Cm1);
  d_free_1d(C);
  d_free_1d(Cp1);
  d_free_1d(rhs);
  d_free_1d(sol);
  d_free_1d(ud);
#endif
  return(0);
}


#ifdef __cplusplus
}
#endif
