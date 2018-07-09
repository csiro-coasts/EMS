/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/sediment/vtransp.c
 *  
 *  Description:
 *  Calculate vertical advection-diffusion
 *  in water column and sediment bed
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vtransp.c 5848 2018-06-29 05:01:15Z riz008 $
 *
 */

#ifdef __cplusplus
#include "stdafx.h"
extern "C" {
#endif

#include "sediments.h"


void reef_scale_depth(sediment_t *sediment, sed_column_t *sm);

static void calc_active_dz(sediment_t *sediment, sed_column_t *sm);
static void cbfiltered(sediment_t *sediment, sed_column_t *sm);
static void biodiffkz(sediment_t *sediment, sed_column_t *sm);
static double bioprofile(char profile, double z, double scale);
static void update_svel_floc(sediment_t *sediment, sed_column_t *sm);
void cohsedcontent(sediment_t *sediment, sed_column_t *sm);
static void update_sedim_wc(sediment_t *sediment, sed_column_t *sm);
static void calc_svel_consolid(sediment_t *sediment, sed_column_t *sm);
static void update_sed_levels(sediment_t *sediment, sed_column_t *sm);
static void update_sedim_sed(sediment_t *sediment, sed_column_t *sm);
static void correct_wc(sediment_t *sediment, sed_column_t *sm);
static void update_dissolved_tr(sediment_t *sediment, sed_column_t *sm);
static void mix_tr_wcbed(sediment_t *sediment, sed_column_t *sm);
static void adsorbed_wc(sediment_t *sediment, sed_column_t *sm);
static void adsorbed_sed(sediment_t *sediment, sed_column_t *sm);
static void truncate_errors(sediment_t *sediment, sed_column_t *sm);
static void sorption(sediment_t *sediment, sed_column_t *sm);
static void sed_tracer_decay(sediment_t *sediment, sed_column_t *sm);
static double erdep_fine(sediment_t *sediment, sed_column_t *sm,sed_tracer_t *tracer);
static double calc_css(sediment_t *sediment, sed_column_t *sm, int k);
static double erdep_coarse(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer);
static void calculate_dzface_sed(sed_column_t *sm);
static void calculate_dzface_wc(sed_column_t *sm);

static void srf_flux_inout(sediment_t *sediment, sed_column_t *sm, 
			   int nt, double *srf_flux_in, double *srf_flux_out_coef);
  
double wgt_tophat(double x,double scale);
double wgt_linear(double x,double scale);
double wgt_parabolic(double x,double scale);
double wgt_gaussian(double x,double scale);

void sinterface_getsvel_custom(void* hmodel, int c, int n,
     char *name, double *svel_wc, int topk_wc, int botk_wc);

#if defined(ENABLE_REF_C_COEF)
static double ref_c_coef(sediment_t *sediment, sed_column_t *sm, int n);
#endif

void update_porosity_wc(sediment_t *sediment, sed_column_t *sm);
void update_porosity_sed(sediment_t *sediment, sed_column_t *sm);

extern void vdiff_sedim_wc(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                           sed_tracer_t *tracercar, double botinflux,
                           double topinflux, double botoutfluxcoef,
                           double topoutfluxcoef, double *dvar);
extern void vdiff_sedim_sed(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                            double botinflux, double topinflux,
                            double botoutfluxcoef, double topoutfluxcoef,
                            double *dvar);
extern void vdiff_dissolved_wc(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                               double botinflux, double topinflux,
                               double botoutfluxcoef,
                               double topoutfluxcoef, double *dvar);
extern void vdiff_dissolved_sed(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer,
                                double botinflux, double topinflux,
                                double botoutfluxcoef,
                                double topoutfluxcoef, double *dvar);


void vert_transport(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;

  if (param->sednz < 1)
    return;

  /* CALCULATE SEDIMENT FLUXES and DISTRIBUTION in WC */
  /* Calculate thickness of the upper most, sediment bed active layer. */
  calc_active_dz(sediment, sm);
  /* Calculate HF filtered sediment concentration in active layer */
  cbfiltered(sediment, sm);
  /* Calculate diffusion coefficient for bioturbation and bioirrigation */
  biodiffkz(sediment, sm);

  /* Find cohesive sediment content in sediment bed; 
     also estimates tss_wc to be used in the next routine NMY2013 */
  cohsedcontent(sediment, sm);
 
 /* Find settling velocities of flocculating particles in wc */
  update_svel_floc(sediment, sm);


  /*UR added  Calculate dz values at wc interfaces */
  calculate_dzface_wc(sm);
  /* Update sediment concentration in wc and find sediment fluxes:
     sm->erdepflux[n] */
  update_sedim_wc(sediment, sm);

  /* UPDATE SEDIMENT THICKNESS */
  /* Now, sediment fluxes at the wc-bottom interface are known and
     sediment thickness and sediment levels can be updated */
  /* First find consolidation velocity */
  calc_svel_consolid(sediment, sm);

  /* Then update sediment grid */
  update_sed_levels(sediment, sm);

  /*UR added  Calculate dz values at sed interfaces */
  calculate_dzface_sed(sm);
  /* CALCULATE SEDIMENT DISTRIBUTION in BOTTOM */
  update_sedim_sed(sediment, sm);

  /* Correct wc grid and sediment concentration in wc taking into account
     changed sediment depth */
  correct_wc(sediment, sm);

  if (param->geomorph > 0) {
    /*UR added  Calculate dz values at interfaces
     * if they gave changed in correct_wc*/
    calculate_dzface_sed(sm);
    calculate_dzface_wc(sm);
  }
  /* Update porosity in wc and sedbed and redistribute concentration if
     necessary */
  update_porosity_wc(sediment, sm);
  update_porosity_sed(sediment, sm);

  /* CALCULATE DISSOLVED TRACER */
  /* dissolved transport in wc and bottom */
  update_dissolved_tr(sediment, sm);
  /* mixing of the dissolved tracer across the wc_sedbed interface */
  mix_tr_wcbed(sediment, sm);

  /* TRACER on SEDIMENT */
  /* adsorbed concentrations in wc */
  adsorbed_wc(sediment, sm);
  /* adsorbed concentrations in bed */
  adsorbed_sed(sediment, sm);
  /* END TRANSPORT */

  /* MISC */
  /* truncate numericall errors */
  truncate_errors(sediment, sm);
  /* Sorption - desorption and decay */
  sorption(sediment, sm);
  sed_tracer_decay(sediment, sm);

  /* Print section */
  /*
     { int k; int tk_sed = sm->topk_sed; int bk_sed = sm->botk_sed; warn("
     water depth = %f", sm->depth_wc); warn(" sedim depth = %f",
     sm->depth_sedi); warn(" total depth = %f",
     sm->depth_wc+sm->depth_sedi); for (k=bk_sed; k<=tk_sed; k++) warn("
     por_sed[%d] = %f", k, sm->por_sed[k]); } */


}

/***************************************/
/* Calculate thickness of an active sediment layer */

static void calc_active_dz(sediment_t *sediment, sed_column_t *sm)
{
  /* sed_params_t *param = sediment->msparam; */
  // double hrip = max(sm->hripples, 0.005);
  sm->depth_sedi = sm->topz_sed - sm->botz_sed;

  // do not update dzactive 2011 (bgc integrator cant handle very thin layers) 
  // sm->dzactive = 0.25 * hrip * sm->depth_sedi / (hrip + sm->depth_sedi);

  if (sm->dzactive > sm->depth_sedi) {
    /* MH 07/2012: included instability handling */
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:calc_active_dz: Thickness of an active layer greater then the total sediment thickness");
    /*sedtag(LFATAL,"sed:vtransp:calc_active_dz"," Thickness of an active layer greater then the total sediment thickness");
      exit(1);*/
    /* END MH */
  }
  /* if (sm->depth_sedi < param->mindepth_sedi) warn("sediment bed almost
     totally eroded"); */
}

/***************************************************/
/* Calculate HF filtered concentrations in the top sediment bed layer.
   HF filtered sediment concentrations are used to estimate
   sediment resuspension, because the numerical scheme for the sediment
   concentration in water column is explicit with respect to the
   particles concentration in bed.
*/

static void cbfiltered(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int tk_sed = sm->topk_sed;
  int n;

  for (n = 0; n < param->ntr; n++) {

    /* 2011 tmp: Restore Markovian property (may cause stability problems)
       sm->cbfilt[n] = 0.5*(sm->tr_sed[n][tk_sed] + sm->tr_sed[n][tk_sed - 1]);
    */

    sm->cbfilt[n] = (sm->cbnm1[n] + sm->cbnm2[n] + sm->cbnm3[n] +
                     sm->cbnm4[n] + sm->tr_sed[n][tk_sed] +
                     sm->tr_sed[n][tk_sed - 1]) / 6.;
    sm->cbnm4[n] = sm->cbnm3[n];
    sm->cbnm3[n] = sm->cbnm2[n];
    sm->cbnm2[n] = sm->cbnm1[n];
    sm->cbnm1[n] = sm->cbfilt[n];
  }
}

/**************************************************/
/* Find the coefficients of the dissolved
   and particulate tracers diffusion in the sediment bed
*/
static void biodiffkz(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;
  int bk_sed = sm->botk_sed;
  double scale   = param->maxbiodepth;
  double biodens = param->biodens;
  double bi_dissol_kz = param->bi_dissol_kz;
  double bt_partic_kz = param->bt_partic_kz;
  double bi_dissol_kz_i = param->bi_dissol_kz_i;
  double bt_partic_kz_i = param->bt_partic_kz_i;
  double z;
  char profile = param->biosedprofile;
  
  /* Handle spatially varying cases */
  if (param->maxbiodepth_spv)
    scale = param->maxbiodepth_spv[sm->col_number-1];
  if (param->biodens_spv)
    biodens = param->biodens_spv[sm->col_number-1];
  if (param->bi_dissol_kz_spv)
    bi_dissol_kz = param->bi_dissol_kz_spv[sm->col_number-1];
  if (param->bi_dissol_kz_i_spv)
    bi_dissol_kz_i = param->bi_dissol_kz_i_spv[sm->col_number-1];
  if (param->bt_partic_kz_spv)
    bt_partic_kz = param->bt_partic_kz_spv[sm->col_number-1];
  if (param->bt_partic_kz_i_spv)
    bt_partic_kz_i = param->bt_partic_kz_i_spv[sm->col_number-1];

  for (k = bk_sed; k <= tk_sed + 1; k++) {
    /* z = fabs(sm->gridz_sed[k]); */
    z = fabs(sm->gridz_sed[tk_sed + 1] - sm->gridz_sed[k]);
    sm->dissol_kz[k] =
      bi_dissol_kz * biodens * bioprofile(profile, z, scale);
    sm->partic_kz[k] =
      bt_partic_kz * biodens * bioprofile(profile, z, scale);
  }
  /* sediment-water interface */
  sm->dissol_kz_i = bi_dissol_kz_i * biodens;
  sm->partic_kz_i = bt_partic_kz_i * biodens;

  if (param->reef_scale_depth >1e-11)
    reef_scale_depth(sediment, sm);// 2013

}

static double bioprofile(char profile, double z, double scale)
{
    double r = 1.;
  switch (profile) {
  case 'c':
    return (wgt_tophat(z, scale));
  case 'l':
    return (wgt_linear(z, scale));
  case 'p':
    return (wgt_parabolic(z, scale));
  case 'g':
    return (wgt_gaussian(z, scale));
  default:
      {
    sedtag(LFATAL,"sed:vtransp:bioprofile","Unknown profile (%c)", profile);
    exit(1);
      }
  }
  return r;
}


/*********************************************/
/* Calculate settling velocities of
   flocculating particles
*/
static void update_svel_floc(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n, n_salt;                        /* tracer index */
  double c_fine, c_salt;
  double a_coef = param->flocprm1;
  double b_coef = param->flocprm2;
  double rab1,rab2,rab3;

  /* tested in sed_optimization */
  n_salt=0;

  for (n = 0; n < param->ntr; n++)
  {
    sm->erflux[n] = 0.;
    /*no sediment flux at the water surface*/
    sm->svel_wc[n][tk_wc+1] = 0.;
  }



  switch (param->flocmode) {
    case 0:
    {
      for (k = bk_wc; k <= tk_wc; k++)
        sm->svel_floc[k] = 0.0;   /* no flocculation */
      break;
    }
    case 1:
    {
      for (k = bk_wc; k <= tk_wc; k++) {
        c_fine = 0.;
        /* get rtracer which contribute to c-fine */
        for (n = 0; n < param->n_vtransp_floc; n++) {
          c_fine += sm->tr_wc[param->vtransp_floc[n]][k];
        }

        c_salt = sm->tr_wc[n_salt][k];
        if (c_salt>0.1)
          sm->svel_floc[k] = -a_coef * pow(c_fine, b_coef); /* power law */
        else
          sm->svel_floc[k] = 0.;
      }
      break;
    }
    case 2:
    {
      // scale with salinity
      for (k = bk_wc; k <= tk_wc; k++) {
	c_salt = sm->tr_wc[n_salt][k];
	if (c_salt>0.1)
	  sm->svel_floc[k] = -a_coef *  (b_coef*c_salt/(1+b_coef*c_salt));
	else
	  sm->svel_floc[k] = 0.;
      }
 
      //     sedtag(LFATAL,"vtransp:update_svel_floc","Option FLOC_MODE = 2 is not available");
      // exit(1);
            
      break;
    }
  case 3:
    {
      // scale with Kz
    for (k = bk_wc; k <= tk_wc; k++) {
      c_salt = sm->tr_wc[n_salt][k];
      if (c_salt>0.1)
	sm->svel_floc[k] = -a_coef * sm->Kz_wc[k];
      else
	sm->svel_floc[k] = 0.;
    }

    // sedtag(LFATAL,"vtransp:update_svel_floc","Option FLOC_MODE = 3 is not available");
    // exit(1);     
      break;
    }
    case 4:
    {

 for (k = bk_wc; k <= tk_wc; k++) {
   c_fine = 0.;
   /* get rtracer which contribute to c-fine */
   for (n = 0; n < param->n_vtransp_floc; n++) {
     c_fine += sm->tr_wc[param->vtransp_floc[n]][k];
   }
   c_salt = sm->tr_wc[n_salt][k];
   if (c_salt>0.1)
     sm->svel_floc[k] = -a_coef*exp(0.945*c_fine - 0.105*c_fine*c_fine);
   else
     sm->svel_floc[k] = 0.;
 }

 // sedtag(LFATAL,"vtransp:update_svel_floc","Option FLOC_MODE = 4 is not available");  exit(1);
      break;
    }

    default:
      break;
  }

  /* Sediment net settling velocity */
  for (n = 0; n < param->n_vtransp_partic; n++)
  {
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_partic[n]];
    for (k = bk_wc; k <= tk_wc; k++) {
      if (tracer->floc)
        sm->svel_wc[param->vtransp_partic[n]][k] = max(-0.1, min(sm->svel_floc[k], tracer->svel[sm->col_number-1]));
      else
        sm->svel_wc[param->vtransp_partic[n]][k] = tracer->svel[sm->col_number-1];

      // NMY13 hindered settling
      if(k > bk_wc) {
	if(param->hindered_svel_patch) {
	  // adhoc formulation 
	  // imposes stronger constraint on conentration fields 
	  // suppresses spurious sources produced over the hard 
	  // substrate by semi-lagrange transport model;
	  rab1 = sm->tss_wc[k-1];
	  rab2 = 10.; 
	  sm->svel_wc[param->vtransp_partic[n]][k] *= max(0, (rab2 - rab1)/rab2);
	}
	else 
	  if(param->hindered_svel)
	    {
	      // by van Rijn "Principles of sed transport in river, 
	      // estuaries and coastal seas", Aqua Publications, 
	      // Amsterdam, 1993 (no sequential pagination) (page 3.16)
	      rab1=sm->tss_wc[k-1] / 2650.;
	      rab2=(1-2.15*rab1) * (1-0.75*pow(rab1,0.33));
	      rab3= max(0,rab2);
	      sm->svel_wc[param->vtransp_partic[n]][k] *= rab3; //vRijn	  
	    }	
      } //end if(k>bk_w) i.e. end of hindered settling section

    }
  }

  for (n = 0; n < param->ntr; n++) {
      sed_tracer_t *tr = &sediment->mstracers[n];
      sinterface_getsvel_custom(sediment->hmodel, sm->hd_col_index, n,
       tr->name, sm->svel_wc[n], sm->topk_wc, sm->botk_wc);
  }
}
/*end */
/*UR-OPT end */


/***************************************************/
/* Find cohesive sediment content
   in the sediment bed
*/
void cohsedcontent(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;
  int bk_sed = sm->botk_sed;
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  double c_fine, c_coarse, c_total;
  /* Find cohesive sediment content in bed */

  for (k = bk_sed; k <= tk_sed; k++) {
    sm->coh_sed[k] = 0.;
    sm->tss_sed[k]=0.;
    c_coarse = 0.;
    c_fine = 0.;

    for (n = 0; n < param->n_vtransp_nd_insed_partic; n++) {
      sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_nd_insed_partic[n]];
      if(tracer->u_scale==1) // if units are kg m-3 (ie non-organic particles)
	sm->tss_sed[k] += sm->tr_sed[param->vtransp_nd_insed_partic[n]][k];

      if (tracer->u_scale==1 && tracer->cohesive)
        c_fine += sm->tr_sed[param->vtransp_nd_insed_partic[n]][k];
      else
        c_coarse += sm->tr_sed[param->vtransp_nd_insed_partic[n]][k];
    }
    c_total = c_fine + c_coarse;
    if (c_total < MINVAL)
      sedtag(LWARN,"sed:vtransp:cohsedcontent"," Low sediment content in bed");
    // sm->coh_sed[k] = 100. * c_fine / (c_total+1.e-13);
    sm->coh_sed[k] = c_fine; //2010 coh_sed given by concentration raher than fraction
  }

for (k = bk_wc; k <= tk_wc; k++) {
    sm->coh_wc[k] = 0.;
    sm->tss_wc[k]=0.;
    c_coarse = 0.;
    c_fine = 0.;
    for (n = 0; n < param->n_vtransp_nd_insed_partic; n++) {
      sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_nd_insed_partic[n]];     
      if(tracer->u_scale==1) // if units are kg m-3 (ie non-organic particles)
	sm->tss_wc[k] += sm->tr_wc[param->vtransp_nd_insed_partic[n]][k];

      if (tracer->u_scale==1 && tracer->cohesive)
        c_fine += sm->tr_wc[param->vtransp_nd_insed_partic[n]][k];
      else
        c_coarse += sm->tr_wc[param->vtransp_nd_insed_partic[n]][k];
    }
    c_total = c_fine + c_coarse;  
    // sm->coh_wc[k] = 100. * c_fine / (c_total+1.e-13);
    sm->coh_wc[k] = c_fine; //2010 coh_wc given by concentration raher than fraction
  }

}
/* end  */
/*UR-OPT end */


/************************************************************************/
/* Advect, diffuse, settle and resuspend
   sediments in water column
*/
static void update_sedim_wc(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  double srf_flux_in, srf_flux_out_coef; 

  /* check for minimum of two sed columns moved to sed_optimization */
  /*
     Loop over all tracers to update sediment concentration in wc and find
     sediment fluxes sm->erdepflux[n], using explicit scheme for the
     erosion term and implicit representation for the settling and
     deposition processes. */

  for (n = 0; n < param->n_vtransp_nd_insed_partic; n++)
  {
    const double zeroval = 0.;
    int nt =  param->vtransp_nd_insed_partic[n];
    // sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_nd_insed_partic[n]];
    sed_tracer_t *tracer = &sediment->mstracers[nt];
    double depfluxcoef;
    double *dvar = d_alloc_1d(param->nz + 1);

    //  nt = param->vtransp_nd_insed_partic[n];

    // NMY 2013 surface flux
    srf_flux_inout(sediment, sm, nt, &srf_flux_in, &srf_flux_out_coef);
    /* Define erosion flux sm->erflux[n] (explicit term [kg m-2 s-1] )
       and coefficient for implicit deposition depfluxcoef (implicit
       term [m s-1] ) */
    if (tracer->cohesive)
      depfluxcoef = erdep_fine(sediment, sm, tracer);
    else
      depfluxcoef = erdep_coarse(sediment, sm, tracer); /* (kg m-2 s-1) */

    /* sediment mixing across wc-sediment bed */
    depfluxcoef += -sm->partic_kz_i /
        min(sm->dz_wc[bk_wc], sm->dz_sed[tk_sed]);
    sm->erflux[nt] += sm->partic_kz_i * sm->cbfilt[nt] /
        min(sm->dz_wc[bk_wc], sm->dz_sed[tk_sed]);

    /* check sediment thickness */
    if ((sm->topz_sed - sm->botz_sed) < param->mindepth_sedi) {
      sm->erflux[nt] = 0.0;
      depfluxcoef *= 0.01;
    }
    /* check erosion/deposition activators */
    if (!tracer->resuspend)
      sm->erflux[nt] = 0.0;
    if (!tracer->deposit)
      depfluxcoef = 0.0;
    // if particles buoyant  
    if (sm->svel_wc[nt][bk_wc] > 0.) {
      sm->erflux[nt] = 0.0;
      depfluxcoef = 0.0;
    }
    //2012 hardsub scaling
    sm->erflux[nt] *= sm->hardsub_scale;
    depfluxcoef *= sm->hardsub_scale;
    // scaling with depth
    if(sm->depth_wc > 500)  sm->erflux[nt] = 0;

    /* Do diffusion, settling, erosion and deposition */ 
    vdiff_sedim_wc(sediment, sm, tracer, tracer, sm->erflux[nt], srf_flux_in, depfluxcoef, srf_flux_out_coef, dvar);

    /* find deposition flux */
    sm->depflux[nt] =
        depfluxcoef * (sm->tr_wc[nt][bk_wc] + dvar[bk_wc]);
    /* find total flux = erosion + deposition */
    sm->erdepflux[nt] = sm->erflux[nt] + sm->depflux[nt];
    /* if there is not enough sediment in active layer reduce erosion
       rate and update concentrations and sediment flux */
    if (sm->tr_sed[nt][tk_sed] * sm->dz_sed[tk_sed] <
        sm->erdepflux[nt] * param->dt) {
      /* transfer this sediment from active layer to wc */
      sm->tr_wc[nt][bk_wc] =
          (sm->tr_wc[nt][bk_wc] * sm->dz_wc[bk_wc] +
          sm->tr_sed[nt][tk_sed] * sm->dz_sed[tk_sed]) /
          sm->dz_wc[bk_wc];

      vdiff_sedim_wc(sediment,sm,tracer,tracer,zeroval,srf_flux_in,zeroval, srf_flux_out_coef,dvar);

      sm->depflux[nt] = 0.0;
      sm->erflux[nt] =
          sm->tr_sed[nt][tk_sed] * sm->dz_sed[tk_sed] / param->dt;
      sm->erdepflux[nt] = sm->erflux[nt] + sm->depflux[nt];
    }

    // if the sediment thickness exceeds critical value and sediments settle 2010
    if(sm->depth_sedi > param->max_thick_sed) {
      if( sm->erdepflux[nt] < 0) {
	for (k = bk_wc; k <= tk_wc; k++)
	  dvar[k] = 0;
	sm->depflux[nt] = 0;
	sm->erflux[nt] = 0;
	sm->erdepflux[nt] = 0;
      }}


    /* update concentration */
    for (k = bk_wc; k <= tk_wc; k++)
      sm->tr_wc[nt][k] += dvar[k];

    d_free_1d(dvar);
  }
}
/* end */
/*UR-OPT end */


/**************************************************/
/* Calculate particles velocity in a
   swelling/consolidating sediment bed
*/
static void calc_svel_consolid(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* tk in sediment */
  double crate = param->consolrate;
  if (param->calcvol_sed && param->consolidate) {
    sm->svel_consolid[bk_sed] = 0.;

    for (k = bk_sed; k <= tk_sed; k++) {
      sm->svel_consolid[k + 1] = sm->svel_consolid[k] -
        (1. / (crate * 86400.)) * (sm->dz_sed[k] / (1. - sm->por_sed[k])) *
        (sm->por_sed[k] - param->finpor_sed);

    }
    /* Limit maximum sediment consolidation distance, in a single time
       step, to the half depth of initial sediment thickness */
    while (fabs(sm->svel_consolid[tk_sed + 1] * param->dt) >
           0.5 * (sm->topz_sed - sm->botz_sed) * param->minpor_sed) {
      for (k = bk_sed; k <= tk_sed + 1; k++) {
        sm->svel_consolid[k] *= 0.5;

      }
    }
  } else {
    for (k = bk_sed; k <= tk_sed + 1; k++)
      sm->svel_consolid[k] = 0.;
  }
  /* one water column cell */
  if (sm->topk_wc == sm->botk_wc) {
    /* MH 07/2012: included instability handling */
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:svel_consolid: At least two wc cells must be present.");
    /*sedtag(LFATAL,"sed:vtransp:svel_consolid"," At least two wc cells must be present.");
      exit(1);*/
    /* END MH */
  }

}

/***************************************************/
/* Update coordinate levels of sediment bed */
static void update_sed_levels(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* tk in sediment */
  int n;                        /* tracer index */
  double flowrate, voidratio, voidr;
  double sigmagridz, newgridz, dhinter;
  double gscale = 1.e-5;
  char etext[MAXSTRLEN];
  if (param->calcvol_sed) {
    /* find void ratio in sediment top layer */
    if (sm->por_sed[tk_sed] < 1.)
      voidratio = sm->por_sed[tk_sed] / (1. - sm->por_sed[tk_sed]);
    else {
      /* MH 07/2012: included instability handling */
      sprintf(etext, "sed:vtransp:update_sed_levels: Sediment volume in bed less/equal 0.0: por_sed = %f, i=%d, j=%d \n",
	      sm->por_sed[tk_sed], sm->i, sm->j);
      i_set_error(sediment->hmodel, sm->col_number, LFATAL, etext);
      /*sedtag(LFATAL,"sed:vtransp:update_sed_levels"," Sediment volume in bed less/equal 0.0: por_sed = %d, i=%d, j=%d \n",sm->por_sed[tk_sed], sm->i, sm->j);*/
      voidratio = 1.;
      /*exit(1);*/
      /* END MH */
    }
    /* find erosion/deposition rate term in equation for the sediment
       thickness */
    flowrate = 0.;
    for (n = 0; n < param->n_vtransp_por_sed; n++) {
      sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_por_sed[n]];
          voidr =
            (sm->erdepflux[param->vtransp_por_sed[n]] >
             0.0) ? voidratio : (tracer->b_dens / tracer->i_conc - 1.);
          flowrate += (1. + voidr) * sm->erdepflux[param->vtransp_por_sed[n]] / tracer->b_dens;
    }
  } else {
    flowrate = 0.;
  }

  /* find new coordinate of the sediment top level */

  newgridz = sm->gridz_sed[tk_sed + 1] +
    param->dt * (sm->svel_consolid[tk_sed + 1] - flowrate);
  sm->gridvel_sed[tk_sed + 1] =
    (newgridz - sm->gridz_sed[tk_sed + 1]) / param->dt;
  sm->gridz_sed[tk_sed + 1] = newgridz;
  /* find coordinate of the second level (bottom of an active layer) */
  newgridz = sm->gridz_sed[tk_sed + 1] - sm->dzactive;
  sm->gridvel_sed[tk_sed] = (newgridz - sm->gridz_sed[tk_sed]) / param->dt;
  sm->gridz_sed[tk_sed] = newgridz;

  /* find coordinates of internal layers */
  for (k = bk_sed; k <= tk_sed + 1; k++)
    sm->tmp_sed[k] = sm->gridz_sed[k];

  if (sm->gridvel_sed[tk_sed + 1] < 0.) {
    for (k = tk_sed; k > bk_sed; k--) {
      if ((sm->tmp_sed[k] - sm->tmp_sed[k - 1]) < param->minseddz)
        sm->tmp_sed[k - 1] = sm->tmp_sed[k] - param->minseddz;
    }
  } else {
    for (k = tk_sed; k > bk_sed; k--)
      sm->gridvel_sed[k - 1] = 0.;
  }

  dhinter = (sm->gridz_sed[tk_sed] - sm->gridz_sed[bk_sed]);
  for (k = tk_sed - 1; k > bk_sed; k--) {
    sigmagridz = sm->gridz_sed[tk_sed] + sm->sigma_sed[k] * dhinter;
    newgridz = sm->tmp_sed[k] + gscale * (sigmagridz - sm->tmp_sed[k]);
    sm->gridvel_sed[k] = (newgridz - sm->gridz_sed[k]) / param->dt;
    sm->gridz_sed[k] = newgridz;
  }

  if (sm->gridz_sed[bk_sed + 1] <= sm->gridz_sed[bk_sed]) {
    /* MH 07/2012: included instability handling */
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:update_sed_levels: Wrong allocation of the grid levels");
    /*sedtag(LFATAL,"sed:vtransp:update_sed_levels"," Wrong allocation of the grid levels");
      exit(1);*/
    /* END MH */
  }
  /* new2 end */

  sm->gridvel_sed[bk_sed] = 0.0;

/*UR-OPT start
  for (k = bk_sed; k <= tk_sed; k++)
    sm->cellz_sed[k] = (sm->gridz_sed[k] + sm->gridz_sed[k + 1]) / 2.;

  sm->botz_sed = sm->gridz_sed[sm->botk_sed];
  sm->topz_sed = sm->gridz_sed[sm->topk_sed + 1];
  sm->depth_sedi = sm->topz_sed - sm->botz_sed;

  for (k = bk_sed; k <= tk_sed; k++)
    sm->dzold_sed[k] = sm->dz_sed[k];

  if (2 > 1) {
    double cbot, ctop;
    cbot = sm->botz_sed;
    for (k = bk_sed; k < tk_sed; k++) {
      ctop = sm->gridz_sed[k + 1];
      sm->dz_sed[k] = ctop - cbot;
      cbot = ctop;
    }
    ctop = sm->topz_sed;
    sm->dz_sed[tk_sed] = ctop - cbot;

  }

  / *UR-OPT 4/2006 */
	/* merge the three loops, remove the unecessary conditional statement */
  sm->botz_sed = sm->gridz_sed[sm->botk_sed];
  sm->topz_sed = sm->gridz_sed[sm->topk_sed + 1];
  sm->depth_sedi = sm->topz_sed - sm->botz_sed;
  {
  	double cbot, ctop;
    cbot = sm->botz_sed;
    /* reduce exit condition by one and merge the three loops */
   	for (k = bk_sed; k < tk_sed; k++)
   	{
    	sm->cellz_sed[k] = (sm->gridz_sed[k] + sm->gridz_sed[k + 1]) / 2.;
    	sm->dzold_sed[k] = sm->dz_sed[k];
      ctop = sm->gridz_sed[k + 1];
      sm->dz_sed[k] = ctop - cbot;
      cbot = ctop;
    }
    /* make up for the reduction from above */
   	sm->cellz_sed[tk_sed] = (sm->gridz_sed[tk_sed] + sm->gridz_sed[tk_sed + 1]) / 2.;
   	sm->dzold_sed[tk_sed] = sm->dz_sed[tk_sed];
    ctop = sm->topz_sed;
    sm->dz_sed[tk_sed] = ctop - cbot;

  }

  /*UR-OPT end */

}
/*end */
/*UR-OPT end */


/* **************************************** */
/* Advect and diffuse sediment
   in benthic layers
*/
static void update_sedim_sed(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* bk in sediment */
  int n;                        /* index */

  /* Update sediment concentration in bottom */
  for (n = 0; n < param->n_vtransp_nd_insed_partic; n++) {
    double zeroval = 0.;
    double *dvar;
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_nd_insed_partic[n]];
    dvar = d_alloc_1d(param->sednz + 1);
    vdiff_sedim_sed(sediment, sm, tracer, zeroval, 
		    sm->erdepflux[param->vtransp_nd_insed_partic[n]],zeroval, zeroval, dvar);
    for (k = bk_sed; k <= tk_sed; k++)
      sm->tr_sed[param->vtransp_nd_insed_partic[n]][k] += dvar[k];
    d_free_1d(dvar);
  }
}
/* end */



/*********************************************/
/* Update coordinate levels
   in water column
*/
static void correct_wc(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_wc = sm->topk_wc;      /* tk in wc */
  int bk_wc = sm->botk_wc;      /* tk in wc */
  int tk_sed = sm->topk_sed;    /* tk in wc */
  int n;                        /* tracer index */
  double ctop, cbot;
  double Ho, Hn, Htotal, newgridz;
  /* Correct for changed water volume */

  /* update dz, gridz and cellz in wc */

  if (param->geomorph > 0) {
    Ho = sm->topz_wc - sm->botz_wc;
    Hn = Ho - sm->gridvel_sed[tk_sed + 1] * param->dt;
    sm->depth_wc = Hn;
    Htotal = Hn + sm->topz_sed - sm->botz_sed;
    for (k = bk_wc; k <= tk_wc; k++) {
      newgridz = sm->gridz_wc[k] * Hn / Ho;
      sm->gridvel_wc[k] = (newgridz - sm->gridz_wc[k]) / param->dt;
      sm->gridz_wc[k] = newgridz;
    }
    sm->gridvel_wc[tk_wc + 1] = 0.;

    for (k = bk_wc; k <= tk_wc; k++)
      sm->cellz_wc[k] = (sm->gridz_wc[k] + sm->gridz_wc[k + 1]) / 2.;
    sm->botz_wc = sm->gridz_wc[bk_wc];
    for (k = bk_wc; k <= tk_wc; k++)
      sm->dzold_wc[k] = sm->dz_wc[k];

    cbot = sm->botz_wc;
    for (k = bk_wc; k < tk_wc; k++) {
      ctop = sm->gridz_wc[k + 1];
      sm->dz_wc[k] = ctop - cbot;
      cbot = ctop;
    }
    ctop = sm->topz_wc;
    sm->dz_wc[tk_wc] = ctop - cbot;

    /* update sediment concentration in wc */
    for (k = bk_wc; k <= tk_wc; k++) {
      for (n = 0; n < param->n_vtransp_nd_insed_partic; n++) {
        sm->tr_wc[param->vtransp_nd_insed_partic[n]][k] =
            sm->tr_wc[param->vtransp_nd_insed_partic[n]][k] * sm->dzold_wc[k] / sm->dz_wc[k];
      }
    }
  } else {
    sm->depth_wc = sm->topz_wc - sm->botz_wc;
  }
}
/* end */
/*UR-OPT end */


/*******************************/
/* Update porosity
   in water column
*/
  void update_porosity_wc(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_wc = sm->topk_wc;      /* tk in wc */
  int bk_wc = sm->botk_wc;      /* bk in wc */
  int n;                        /* index */
  double sed_vol;
  /* new: update porosity in wc and redistribute concentration if
     necessary */

  for (k = bk_wc; k <= tk_wc; k++) {
    sed_vol = 0.0;
    for (n = 0; n < param->n_vtransp_por_wc; n++) {
      sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_por_wc[n]];
      if (fabs(sm->tr_wc[param->vtransp_por_wc[n]][k]) < MINVAL)
        sm->tr_wc[param->vtransp_por_wc[n]][k] = 0.;
      sed_vol += sm->tr_wc[param->vtransp_por_wc[n]][k] / tracer->b_dens;
    }

    sm->porold_wc[k] = sm->por_wc[k];
    sm->por_wc[k] = 1. - sed_vol;
    if (sm->por_wc[k] < param->minpor_wc) {
      double dm;
      if (k < tk_wc) {
        for (n = 0; n < param->n_vtransp_por_wc; n++) {
          dm =
            sm->tr_wc[param->vtransp_por_wc[n]][k] * (1. -
            (1. - param->minpor_wc) / sed_vol) * sm->dz_wc[k];
           sm->tr_wc[param->vtransp_por_wc[n]][k] *= (1. - param->minpor_wc) / sed_vol;
           sm->tr_wc[param->vtransp_por_wc[n]][k + 1] += dm / sm->dz_wc[k + 1];
           sm->por_wc[k] = param->minpor_wc;
        }
      } else
        sedtag(LWARN,"sed:vtransp:por_wc","Too much sediment in wc.");
    } else {
      if (sm->por_wc[k] > 1.) {
	/* MH 07/2012: included instability handling */
	i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:por_wc: Negative concentration of sediment in wc");
        /*sedtag(LFATAL,"sed:vtransp:por_wc"," Negative concentration of sediment in wc");
	  exit(1);*/
	/* END MH */
      }
    }
  }
}/*end update_porosity_wc*/
/*UR-OPT end */


/*********************************************************/
/* update porosity
   in sediment bed
*/
void update_porosity_sed(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sed */
  int bk_sed = sm->botk_sed;    /* bk in sed */
  int n,m;                        /* tracer index */
  double sed_vol;
  /* update porosity in sediment bed and redistribute sediment when
     necessary */
  for (k = bk_sed; k <= tk_sed; k++) {
    sed_vol = 0.0;
    for (n = 0; n < param->n_vtransp_por_sed; n++) {
      sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_por_sed[n]];
      if (fabs(sm->tr_sed[param->vtransp_por_sed[n]][k]) < MINVAL)
        sm->tr_sed[param->vtransp_por_sed[n]][k] = 0.;
      sed_vol += sm->tr_sed[param->vtransp_por_sed[n]][k] / tracer->b_dens;
    }
    sm->porold_sed[k] = sm->por_sed[k];
    sm->por_sed[k] = 1. - sed_vol;
    
    //NMY Oct11: if initial concentration of sediments is inconsistent with the init grid structure
    // rescale init concentrations and update porosity
    if ( param->nstep <= 1 && sm->por_sed[k] < param->minpor_sed) {
      for (n = 0; n < param->n_vtransp_por_sed; n++) {
	m = param->vtransp_por_sed[n];
	sm->tr_sed[m][k] *= (1-2.*param->minpor_sed)/sed_vol; 
      }
      sm->por_sed[k] = 2.*param->minpor_sed;
   }


    // this seem become redundant now
    if (sm->por_sed[k] < param->minpor_sed) {
      double dm;
      if (k < tk_sed) {
	// redistribute excess mass into the upper layer
        for (n = 0; n < param->n_vtransp_por_sed; n++) {
	  m = param->vtransp_por_sed[n];
          dm =
           sm->tr_sed[m][k] * (1. -
           (1. - param->minpor_sed) / sed_vol) * sm->dz_sed[k];
          sm->tr_sed[m][k] *= (1. - param->minpor_sed) / sed_vol;
          sm->tr_sed[m][k + 1] += dm / sm->dz_sed[k + 1];
          sm->por_sed[k] = param->minpor_sed;
        }
      } else {
        sedtag(LERROR,"sed:vtransp:update_porosity_sed"," Too much sediment in bottom.");
      }
    } else {
      if (sm->por_sed[k] > 1.){
	/* MH 07/2012: included instability handling */
	i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:pr_sed: Negative concentration of sediment in bottom");
        /*sedtag(LFATAL,"sed:vtransp:pr_sed"," Negative concentration of sediment in bottom");
	  exit(1);*/
	/* END MH */
      }
    }
  }
}
/*end update_porosity_sed */
/*UR-OPT end */

/**************************
sets wc surface influx and outflux_coefficient for a tracer nt
updates: srf_flux_in, srf_flux_out_coef 
*/
static void srf_flux_inout(sediment_t *sediment, sed_column_t *sm,
			   int nt, double *srf_flux_in, double *srf_flux_out_coef)
{
  sed_tracer_t *tr = &sediment->mstracers[nt];
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  double srf_flux = sm->tr_srf_flux[nt];

  // NMY 2013 surface flux
  if(srf_flux > 0 ) {  // outflux
    *srf_flux_in = 0;
    if (sm->tr_wc[nt][tk_wc] > 1e-21) 
      *srf_flux_out_coef = srf_flux / sm->tr_wc[nt][tk_wc];
    else
      *srf_flux_out_coef = 0;
  } 
  else                       // influx
    {
      *srf_flux_in = srf_flux;
      *srf_flux_out_coef = 0;	
    }
  
}

/********************************************/
/* Advect and diffuse dissolved tracer
   in water column and sediment bed
*/
static void update_dissolved_tr(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* bk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  double srf_flux_in, srf_flux_out_coef; 

  /* interstitial water velocity*porosity */
  /* do this test at the start once - if possible this should go into init */
  if (tk_wc == bk_wc) { /* one wc layer */
    /* MH 07/2012: included instability handling */
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:update_dissolved_tr: At least two wc cells must be present");
    /*sedtag(LFATAL,"sed:vtransp:update_dissolved_tr"," At least two wc cells must be present");
      exit(1);*/
    /* END MH */
  }



  sm->watvel_sed[bk_sed] = 0.;
  for (k = bk_sed; k <= tk_sed; k++)
    sm->watvel_sed[k + 1] = sm->watvel_sed[k] -
      (sm->por_sed[k] * sm->dz_sed[k] -
       sm->porold_sed[k] * sm->dzold_sed[k]) / param->dt;
  /* wc water velocity */
  sm->watvel_wc[tk_wc + 1] = 0.;
  for (k = tk_wc; k >= bk_wc; k--)
    sm->watvel_wc[k] = sm->watvel_wc[k + 1] +
      (sm->por_wc[k] * sm->dz_wc[k] -
       sm->porold_wc[k] * sm->dzold_wc[k]) / param->dt;
  sm->watvel_wc[bk_wc] = sm->watvel_sed[tk_sed + 1];

  if (sm->watvel_sed[tk_sed + 1] > 0.) { // if deposition event
    /* Update dissolved concentration first in bottom then in water */
    for (n = 0; n < param->n_vtransp_dissolv; n++) {
      double zeroval = 0.;
      double botinflux = 0.;
      double topoutfluxcoef = 0.;
      double *dvar;
      int nt = param->vtransp_dissolv[n];
      sed_tracer_t *tracer = &sediment->mstracers[nt];
      dvar = d_alloc_1d(param->sednz + 1);    
      topoutfluxcoef = sm->watvel_sed[tk_sed + 1];

      vdiff_dissolved_sed(sediment, sm, tracer, botinflux, zeroval,
        zeroval, topoutfluxcoef, dvar);

      for (k = bk_sed; k <= tk_sed; k++)
        sm->tr_sed[nt][k] += dvar[k];
      d_free_1d(dvar);

      /* explicit influx in water from sediment */
      // NMY 2013 surface flux
      srf_flux_inout(sediment, sm, nt, &srf_flux_in, &srf_flux_out_coef);
      botinflux = sm->watvel_sed[tk_sed + 1] * sm->tr_sed[nt][tk_sed];
      if (!param->geomorph && (strcmp("temp",tracer->name) == 0 ||
          strcmp("salt",tracer->name) == 0) )
        botinflux = 0.;     /* Cut temp and salt influx */
      dvar = d_alloc_1d(param->nz + 1);
      vdiff_dissolved_wc(sediment, sm, tracer, botinflux, srf_flux_in,
			 zeroval, srf_flux_out_coef, dvar);
      for (k = bk_wc; k <= tk_wc; k++)
        sm->tr_wc[nt][k] += dvar[k];
      d_free_1d(dvar);
    }
  } else { // if erosion event
    /* Update dissolved concentration first in water then in bottom */
    for (n = 0; n < param->n_vtransp_dissolv; n++) {
      double zeroval = 0.;
      double topinflux = 0.;
      double botoutfluxcoef = 0.;
      double *dvar;
      int nt = param->vtransp_dissolv[n];
      sed_tracer_t *tracer = &sediment->mstracers[nt];
      dvar = d_alloc_1d(param->nz + 1);
      // NMY 2013 surface flux
      // sets srf_flux_in, srf_flux_out_coef
      srf_flux_inout(sediment, sm, nt, &srf_flux_in, &srf_flux_out_coef);
      botoutfluxcoef = sm->watvel_sed[tk_sed + 1];
      if (!param->geomorph && (strcmp("temp",tracer->name) == 0 ||  strcmp("salt",tracer->name) == 0))
        botoutfluxcoef = 0.;  /* Cut temp and salt outflux */
      vdiff_dissolved_wc(sediment, sm, tracer, zeroval, srf_flux_in,
			 botoutfluxcoef, srf_flux_out_coef, dvar);

      for (k = bk_wc; k <= tk_wc; k++)
        sm->tr_wc[nt][k] += dvar[k];

      d_free_1d(dvar);

      /* explicit influx in sediment from water */
      topinflux = sm->watvel_sed[tk_sed + 1] * sm->tr_wc[nt][bk_wc];
      dvar = d_alloc_1d(param->sednz + 1);
      vdiff_dissolved_sed(sediment, sm, tracer, zeroval, topinflux,
                          zeroval, zeroval, dvar);
      for (k = bk_sed; k <= tk_sed; k++)
        sm->tr_sed[nt][k] += dvar[k];
      d_free_1d(dvar);

    }
  }
}
/* end update_dissolved_tr */
/*UR-OPT end */

/*********************************************/
/* Diffusion of the dissolved tracer
   across the water column and sediment bed
*/
static void mix_tr_wcbed(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  /* mixing of dissolved tracers across the "wc - sediment bed" interface */
  /* solid particles have already been mixed */
  for (n = 0; n < param->n_vtransp_dissolv; n++) {

    //sed_tracer_t *tracer = &sediment->mstracers[n];
    double k0;
    double A, B, Cw, D, R;
    double E, Cb, minz, m0, m1;
    double dt = param->dt;
    k0 = sm->dissol_kz_i;
    if (k0 > 0.) {
      m0 = sm->tr_wc[param->vtransp_dissolv[n]][bk_wc] * sm->dz_wc[bk_wc] * sm->por_wc[bk_wc];
      minz = min(sm->dz_wc[bk_wc], sm->dz_sed[tk_sed]);
      R = k0 / (sm->dz_wc[bk_wc] * sm->por_wc[bk_wc] * minz);
      D = k0 / (sm->dz_sed[tk_sed] * sm->por_sed[tk_sed] * minz);
      B = (sm->tr_wc[param->vtransp_dissolv[n]][bk_wc] * sm->dz_wc[bk_wc] * sm->por_wc[bk_wc] +
           sm->tr_sed[param->vtransp_dissolv[n]][tk_sed] * sm->dz_sed[tk_sed] *
           sm->por_sed[tk_sed])
        / (R + D);
      A =
        sm->tr_wc[param->vtransp_dissolv[n]][bk_wc] * sm->dz_wc[bk_wc] * sm->por_wc[bk_wc] -
        D * B;
      E = A * exp(-(R + D) * dt);
      Cw = (E + D * B) / (sm->dz_wc[bk_wc] * sm->por_wc[bk_wc]);
      Cb = (-E + R * B) / (sm->dz_sed[tk_sed] * sm->por_sed[tk_sed]);
      sm->tr_wc[param->vtransp_dissolv[n]][bk_wc] = Cw > 0. ? Cw : 0.;
      sm->tr_sed[param->vtransp_dissolv[n]][tk_sed] = Cb > 0. ? Cb : 0.;
      /* calculate diffusion flux */
      m1 = sm->tr_wc[param->vtransp_dissolv[n]][bk_wc] * sm->dz_wc[bk_wc] * sm->por_wc[bk_wc];
      sm->erdepflux[param->vtransp_dissolv[n]] += (m1 - m0) / dt;
    }
  }
}


/*********************************************/
/* Advect, diffuse, resuspend and settle
   adsorbed tracer in water column
*/
void adsorbed_wc(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* index */

  for (n = 0; n < param->n_vtransp_adsorb; n++) {
    const double zeroval = 0.;
    double *dvar, depfluxcoef;
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_adsorb[n]];
    int ns = tracer->carriernum;
    sed_tracer_t *tracercar = &sediment->mstracers[ns];
    dvar = d_alloc_1d(param->nz + 1);
    /* Define erosion flux sm->erflux[n] (explicit term [kg m-2 s-1] )
       and coefficient for implicit deposition depfluxcoef (implicit
       term [m s-1] ) */

    if (sm->erdepflux[ns] < 0.) {
      depfluxcoef = sm->erdepflux[ns] / (sm->tr_wc[ns][bk_wc] + MINVAL);
      sm->erflux[param->vtransp_adsorb[n]] = 0.;
    } else {
      depfluxcoef = 0.;     /* (kg m-2 s-1) */
      sm->erflux[param->vtransp_adsorb[n]] =
        sm->erdepflux[ns] * sm->tr_sed[param->vtransp_adsorb[n]][tk_sed] /
        (sm->tr_sed[ns][tk_sed] + MINVAL);
    }
    if ((sm->topz_sed - sm->botz_sed) < param->mindepth_sedi)
      sm->erflux[param->vtransp_adsorb[n]] = 0.0;

    /* Do diffusion, settling, erosion and deposition */
    vdiff_sedim_wc(sediment, sm, tracer, tracercar, sm->erflux[param->vtransp_adsorb[n]],
                     zeroval, depfluxcoef, zeroval, dvar);
    /* find deposition flux */
    sm->depflux[param->vtransp_adsorb[n]] = depfluxcoef *
        (sm->tr_wc[param->vtransp_adsorb[n]][bk_wc] + dvar[bk_wc]);
    /* find total flux = erosion + deposition */
    sm->erdepflux[param->vtransp_adsorb[n]] =
        sm->erflux[param->vtransp_adsorb[n]] + sm->depflux[param->vtransp_adsorb[n]];

    /* if there is not enough tracer in active layer reduce erosion
       rate and update concentrations and sediment flux */
    if (sm->tr_sed[param->vtransp_adsorb[n]][tk_sed] * sm->dz_sed[tk_sed] <
        sm->erdepflux[param->vtransp_adsorb[n]] * param->dt) {
      /* transfer sediment from active layer to wc */
      sm->tr_wc[param->vtransp_adsorb[n]][bk_wc] =
          (sm->tr_wc[param->vtransp_adsorb[n]][bk_wc] * sm->dz_wc[bk_wc] +
          sm->tr_sed[param->vtransp_adsorb[n]][tk_sed] * sm->dz_sed[tk_sed]) /
          sm->dz_wc[bk_wc];
      vdiff_sedim_wc(sediment, sm, tracer, tracercar, zeroval, zeroval,
                     zeroval, zeroval, dvar);
      sm->depflux[param->vtransp_adsorb[n]] = 0.0;
      sm->erflux[param->vtransp_adsorb[n]] =
        sm->tr_sed[param->vtransp_adsorb[n]][tk_sed] * sm->dz_sed[tk_sed] / param->dt;
      sm->erdepflux[param->vtransp_adsorb[n]] =
          sm->erflux[param->vtransp_adsorb[n]] + sm->depflux[param->vtransp_adsorb[n]];
    }
    /* update concentration */
    for (k = bk_wc; k <= tk_wc; k++)
      sm->tr_wc[param->vtransp_adsorb[n]][k] += dvar[k];
    d_free_1d(dvar);

  }
}

/**************************************************/
/* Advect and diffuse adsorbed tracer
   in sediment bed
*/
void adsorbed_sed(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* bk in sediment */
  int n;                        /* tracer index */

  for (n = 0; n < param->n_vtransp_adsorb; n++) {
    double zeroval = 0.;
    double *dvar;
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_adsorb[n]];
    dvar = d_alloc_1d(param->sednz + 1);
    vdiff_sedim_sed(sediment, sm, tracer, zeroval, sm->erdepflux[param->vtransp_adsorb[n]],
                    zeroval, zeroval, dvar);
    for (k = bk_sed; k <= tk_sed; k++) {
      sm->tr_sed[param->vtransp_adsorb[n]][k] += dvar[k];
    }
    d_free_1d(dvar);
  }
}
/* end */
/*UR-OPT end */


/*******************************************/
/* Cut numerical truncation errors */
void truncate_errors(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* bk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  for (n = 0; n < param->ntr; n++) {
    /* sed_tracer_t* tracer = &sediment->mstracers[n]; */
    for (k = bk_wc; k <= tk_wc; k++)
      if (fabs(sm->tr_wc[n][k]) < MINVAL)
        sm->tr_wc[n][k] = 0.;
    for (k = bk_sed; k <= tk_sed; k++)
      if (fabs(sm->tr_sed[n][k]) < MINVAL)
        sm->tr_sed[n][k] = 0.;
  }
}

/********************************************/
/* Calculate sorption/desorption exchange */
void sorption(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* tk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */

  for (n = 0; n < param->n_vtransp_adsorb; n++) {
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_adsorb[n]];
    int nd, ns;
    int np = param->vtransp_adsorb[n];
    nd = tracer->dissolvednum;
    ns = tracer->carriernum;
    if (nd >= 0) {          /* if there is sorption-desorption */
      sed_tracer_t *tracerd = &sediment->mstracers[nd];
      double Kd = tracer->adsorbkd;
      double rate = tracer->adsorbrate;
      double A, B, Cp, Cw, M, D, R;
      if (rate > 0.) {      /* non equilibrium sorption */
        double rate_wc, rate_sed;
        double dt = param->dt;
        rate = 1. / (24. * 3600. * rate);
        if (!tracerd->dissol){
	  /* MH 07/2012: included instability handling */
	  i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:sorption: desorbed traced must be dissolved");
          /*sedtag(LFATAL,"sed:vtransp:sorption"," desorbed traced must be dissolved");
	    exit(1);*/
	  /* END MH */
        }
        for (k = bk_wc; k <= tk_wc; k++) {
          R = (sm->tr_wc[ns][k] / sm->por_wc[k]) * Kd;
          B =
            (sm->tr_wc[np][k] +
             sm->tr_wc[nd][k] * sm->por_wc[k]) / (R + 1.);
          A = sm->tr_wc[nd][k] * sm->por_wc[k] - B;
          rate_wc = rate / (R + 1.);
          D = A * exp(-rate_wc * (R + 1) * dt);
          Cp = -D + R * B;
          Cw = (D + B) / sm->por_wc[k];
          sm->tr_wc[np][k] = Cp > 0. ? Cp : 0.;
          sm->tr_wc[nd][k] = Cw > 0. ? Cw : 0.;
        }
        for (k = bk_sed; k <= tk_sed; k++) {
          R = (sm->tr_sed[ns][k] / sm->por_sed[k]) * Kd;
          B =
            (sm->tr_sed[np][k] +
             sm->tr_sed[nd][k] * sm->por_sed[k]) / (R + 1.);
          A = sm->tr_sed[nd][k] * sm->por_sed[k] - B;
          rate_sed = rate / (R + 1.);
          D = A * exp(-rate_sed * (R + 1) * dt);
          Cp = -D + R * B;
          Cw = (D + B) / sm->por_sed[k];
          sm->tr_sed[np][k] = Cp > 0. ? Cp : 0.;
          sm->tr_sed[nd][k] = Cw > 0. ? Cw : 0.;


        }
      } else {              /* equilibrium sorption */

        for (k = bk_wc; k <= tk_wc; k++) {
          M = sm->tr_wc[np][k] + sm->tr_wc[nd][k] * sm->por_wc[k];
          R = sm->tr_wc[ns][k] * Kd;
          Cw = M / (R + sm->por_wc[k]);
          Cp = M - Cw * sm->por_wc[k];
          sm->tr_wc[np][k] = Cp > 0. ? Cp : 0.;
          sm->tr_wc[nd][k] = Cw > 0. ? Cw : 0.;

        }
        for (k = bk_sed; k <= tk_sed; k++) {
          M = sm->tr_sed[np][k] + sm->tr_sed[nd][k] * sm->por_sed[k];
          R = sm->tr_sed[ns][k] * Kd;
          Cw = M / (R + sm->por_sed[k]);
          Cp = M - Cw * sm->por_sed[k];
          sm->tr_sed[np][k] = Cp > 0. ? Cp : 0.;
          sm->tr_sed[nd][k] = Cw > 0. ? Cw : 0.;

        }
      }
    }
  }
}
/*end */
/*UR-OPT end  */


/**************************************/
/* Tracer decay */
void sed_tracer_decay(sediment_t *sediment, sed_column_t *sm)
{
  sed_params_t *param = sediment->msparam;
  int k;
  int tk_sed = sm->topk_sed;    /* tk in sediment */
  int bk_sed = sm->botk_sed;    /* bk in sediment */
  int tk_wc = sm->topk_wc;
  int bk_wc = sm->botk_wc;
  int n;                        /* tracer index */
  /* decay */
  for (n = 0; n < param->n_vtransp_decay; n++) {
    sed_tracer_t *tracer = &sediment->mstracers[param->vtransp_decay[n]];
    double rate = tracer->decay;
    double dt = param->dt;
    rate = 1. / (24. * 3600. * rate);
    for (k = bk_wc; k <= tk_wc; k++)
      sm->tr_wc[param->vtransp_decay[n]][k] *= exp(-rate * dt);
    for (k = bk_sed; k <= tk_sed; k++)
      sm->tr_sed[param->vtransp_decay[n]][k] *= exp(-rate * dt);
  }
}
/* end */
/*UR-OPT end */

/********************************************/
/* Calculate erosion rate on cohesive bed, in kg m-2 s-1 */

static double erdep_fine(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer)
{
  sed_params_t *param = sediment->msparam;
  int tk_sed = sm->topk_sed;
  int bk_wc = sm->botk_wc;
  double rho_w = 1025.;
  double bs;                    /* bottom shear stress */
  double nxsbs;                 /* normalized excess bottom shear stress */
  double depfluxcoef;           /* deposition term coefficient for
                                   implicit calculation */
  double r;
  double v;
  int n = tracer->n;
 // double flux_scaling = 1;
  char etext[MAXSTRLEN];


  /* if unints are not kg m-3, they are set to mg m-3 */
  //  if (strcmp(tracer->units,"kg m-3") != 0)
  //  flux_scaling=1.e6;

  /* Calculate critical shear stress of erosion for cohesive sediment */
  sm->css[tk_sed] = calc_css(sediment, sm, tk_sed);
  /* Calculate erosion term explicitly */
  bs = rho_w * sm->ustrcw_skin * sm->ustrcw_skin; /* skin friction */
  if (sm->css[tk_sed] <= 0.){
    /* MH 07/2012: included instability handling */
    sprintf(etext, "sed:vtransp:erdep_fine: Negative critical shear on cohesive sed.bed %f",sm->css[tk_sed]);
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, etext);
    /*sedtag(LFATAL,"sed:vtransp:erdep_fine"," Negative critical shear on cohesive sed.bed %f",sm->css[tk_sed]);
      exit(1);*/
    /* END MH */
  }
  nxsbs = (bs > sm->css[tk_sed]) ? (bs / sm->css[tk_sed] - 1.0) : 0.0;
/* normalised excess stress */

  /* find volumetric fraction of the given sediment class in bottom */
  //  r = (sm->cbfilt[n] / tracer->b_dens) / (1. - sm->por_sed[tk_sed]);

  r = sm->cbfilt[n]/sm->tss_sed[tk_sed];

  sm->erflux[n] = param->erflux_scale *
      r * 0.002 * sm->css[tk_sed] * nxsbs;  /* erosion flux, in kg m-2 s-1 */
  /* Define deposition term coefficient for implicit calculation */
  nxsbs = (bs < sm->css_dep) ? (1.0 - bs / sm->css_dep) : 0.0;  
  /* normalised excess stress */

  v = sm->svel_wc[n][bk_wc];
  depfluxcoef = v * nxsbs;

  return (depfluxcoef);
}

/**********************************************************************/

/* Routine to calculate critical shear stress of erosion.
 * The units of the value returned here are N m-2.
 */
static double calc_css(sediment_t *sediment, sed_column_t *sm, int k)
{
  sed_params_t *param = sediment->msparam;
  double vr;                    /* calculated void ratio */
  double rho_b, rho_d;          /* bulk (rho_b) and dry (rho_d) density of
                                   the sediment layer */
  double rho_w = 1025.;         /* water density */
  double rho_p = 2650.;         /* density of sediment particles */
/*  Govindaraju et al., 1999,              Mitcheler           Hwang and Mehta
* vr_i=17, rho_b = 1113, rho_d=144 kg/cub.m css = 0.47 N/m**2        0.53
* vr_i=50, rho_b=1056, rho_d=50,            css = 0.29               0.05
*   Torfs et al., 1996,
* vr_f=1.6, rho_b=1631, rho_d=1000,         css = 1.67                0.84
* vr_f=3, rho_b = 1419, rho_d=650 kg/cub.m  css = 1.23
*   Govindaraju et al., 1999,
* vr_f=5.4, rho_b=1271, rho_d= 406,         css = 0.9
*/
  vr = sm->por_sed[k] / (1. - sm->por_sed[k]);
  if (vr < 0.) {
    /* MH 07/2012: included instability handling */
    i_set_error(sediment->hmodel, sm->col_number, LFATAL, "sed:vtransp:calc_css: negative void ratio in sediment bed");
    /*sedtag(LFATAL,"sed:vtransp:calc_css"," negative void ratio in sediment bed");
      exit(1);*/
    /* END MH */
  }
  rho_b = (vr * rho_w + rho_p) / (1. + vr);
  rho_d = rho_p / (1. + vr);
/* Critical shear stress for erosion */
  switch (param->cssmode) {
  case 0:
    // return (param->css);
    return (sm->css[k]);
    break;
  case 1:
    /* (Hwang and Mehta 1989) */
    if ((rho_b = rho_b * 0.001) < 1.065)
      rho_b = 1.065;            /* transfer to g/cub.cm */
    return (0.883 * pow((rho_b - 1.065), 0.2) + 0.05);
    break;
  case 2:
    /* (Mitcheller et. all 1996) */
    return (0.015 * pow((rho_b - 1000.), 0.73));
    break;
  case 3:
    /* (Ockenden and Delo 1988)) */
    return (0.0012 * pow(rho_d, 1.2));
    break;
  default:
    break;
  case 4:
    return (sm->css[k]);
    break;
  }
  return (sm->css[k]);
}

/*********************************************************************/

/* Calculate erosion rate on non-cohesive sediment bed, in kg m-2 s-1 */

static double erdep_coarse(sediment_t *sediment, sed_column_t *sm, sed_tracer_t *tracer)
{

  int tk_sed = sm->topk_sed;    /* number of active sed layer */
  int bk_wc = sm->botk_wc;      /* number of wc bottom cell */
  double rho_w = 1025.;         /* water density */
  double bs;                    /* skin friction */
  double css_coarse;            /* critical shear stress */
  double nxsbs;                 /* normalized excess skin friction */
  double gamma0 = 0.002;        /* empirical coefficient */
  int n = tracer->n;            /* tracer number */
  double svel = tracer->svel[sm->col_number-1];   /* settling velocity */
  double ref_c_c;               /* multiplier transforming cell averaged
                                   concentration to ref conc */
  double depfluxcoef;           /* deposition flux coefficient for
                                   implicit calculations */

  bs = rho_w * sm->ustrcw_skin * sm->ustrcw_skin;
  css_coarse =
    max(sm->css[tk_sed],
        rho_w * 0.25 * svel * 0.25 * svel);
  nxsbs = (bs > css_coarse) ? (bs / css_coarse - 1.0) : 0.0;
  sm->ref_c_eq[n] = sm->cbfilt[n] * gamma0 * nxsbs / (1. + gamma0 * nxsbs);
  sm->erflux[n] = -svel * sm->ref_c_eq[n];

#if defined(ENABLE_REF_C_COEF)
  ref_c_c = ref_c_coef(sediment, sm, n);
#else
  ref_c_c = 1.;                 /* test */
#endif

  sm->ref_c[n] = ref_c_c * sm->tr_wc[n][bk_wc];
  depfluxcoef = svel * ref_c_c;
  return (depfluxcoef);

}

/* Calculate dz values at interfaces */
static void calculate_dzface_sed(sed_column_t *sm)
{
  int k;
  for (k = sm->botk_sed + 1; k <= sm->topk_sed; k++)
    sm->dzface_sed[k] = (sm->dz_sed[k - 1] + sm->dz_sed[k]) / 2.;
}


static void calculate_dzface_wc(sed_column_t *sm)
{
  int k;
  for (k = sm->botk_wc + 1; k <= sm->topk_wc; k++)
    sm->dzface_wc[k] = (sm->dz_wc[k - 1] + sm->dz_wc[k]) / 2.;
}


#if defined(ENABLE_REF_C_COEF)
/***************************************************************/

/* Define ( actual reference concentration / averaged concentration).
   Simplest equilibrium case, with zero vertical flux.
*/

double ref_c_coef(sediment_t *sediment, sed_column_t *sm, int n)
{
  sed_tracer_t *tracer = &sediment->mstracers[n];
  int bk_wc = sm->botk_wc;      /* number of the near bottom cell in water
                                   column */
  double v;                     /* settling velocity (defined as positive
                                   value here) */
  double zref;                  /* reference height */
  double gama;                  /* parameter: gama = v/diffcoef */
  double dzc = sm->dz_wc[bk_wc];  /* thickness of the near bottom cell (in
                                     water column) */
  double c0;                    /* concentration in the near bottom cell */
  double ref_c_c;               /* actual concentration at the reference
                                   height */
  double K0 = 10.e-6;
  double pp;
  double delta;
  double kapu;


  zref = 7. * 30. * sm->z0_skin;
  dzc = dzc;

  v = -sm->svel_wc[n][bk_wc];

  c0 = sm->tr_wc[n][bk_wc];
  kapu = (sm->Kz_wc[bk_wc + 1] - K0) / dzc;
  /* kapu = 0.4 * sm->ustrcw; */
  delta = (kapu < K0) ? 1.0 : (K0 / kapu);
  gama = v / K0;
  pp = -1.0 + delta * gama;
  if (fabs(pp) < 0.01) {
    ref_c_c = dzc / log((dzc + delta) / delta);
    ref_c_c = ref_c_c * (1. / (zref + delta));
  } else {
    double pdzc;
    double pzref;
    pdzc = pow(delta / (dzc + delta), pp);
    pzref = pow(delta / (zref + delta), pp);
    ref_c_c = (dzc * pp) / (1. - pdzc);
    ref_c_c = ref_c_c * (pzref / (zref + delta));
  }

  return (max(ref_c_c, 1.));
}
#endif


#ifdef __cplusplus
}
#endif
