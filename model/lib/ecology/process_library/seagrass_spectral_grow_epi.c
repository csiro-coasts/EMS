/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_spectral_grow_epi.c
 *  
 *  Description:
 *  
 *  Seagrass growth from spectrally-resolved light field, and from water column and porewater nitrate, ammonia and phosphorus.
 *
 *  Options: seagrass_spectral_grow_epi(Z|H|D|P)
 *  
 *  Z - Zostera (SG)
 *  H - Halophila (SGH)
 *  D - Deep Halophila (SGD)
 *  P - Posidonia (SGP)
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: seagrass_spectral_grow_epi.c 7199 2022-09-14 06:36:12Z bai155 $
 *
 */

/* 
 *  Model described in: Baird, M. E., M. P. Adams, R. C. Babcock, K. Oubelkheir, M. Mongin, 
 *                      K. A. Wild-Allen, J. Skerratt, B. J. Robson, K. Petrou, P. J. Ralph, 
 *                      K. R. O'Brien, A. B. Carter, J. C. Jarvis, M. A. Rasheed (2016) 
 *		        A biophysical representation of seagrass growth for application in a 
 *                      complex shallow-water biogeochemical model Ecol. Mod. 325: 13-27.  
 */

/* Notes: 1. Mass balance assumes that if there are any seagrass species, that 'Z' is one of them.

          2. Includes preferential NH4 uptake, multiple sediment layer roots, and minimium light requirement (MLR) and
             determined respiration. 

	  3. Additional to Baird et al., (2016): (1) preferential nutrient uptake by roots; (2) oxygen balance includes NO3.

	  4. Mass balance acheieved by putting NO3_sed into TN, and NO3_sed*48/12 into BOD.

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "constants.h"
#include "seagrass_spectral_grow_epi.h"

#define EPS_DIN 1.0e-12
#define EPS_DIP 1.0e-12
#define unitch 1000.0

int ginterface_getsedtopk(void *model, int b);
int ginterface_getsedbotk(void *model, int b);

typedef struct {
  int do_mb;                  /* flags */
  char species;
  
  /*
   * parameters
   */
  double umax_t0;
  double omega;
  double KN;
  double KP;
  double SGfrac;
  double SGtransrate;
  double SGrootdepth;
  double SGmlr;
  double SG_mL;
  double SGorient;

  int rootbot_i;
  int nsed;
  
  /*
   * epis
   */
  int SG_N_i;
  int SG_N_pr_i;
  int SG_N_gr_i;
  int EpiTN_i;
  int EpiTP_i;
  int EpiTC_i;
  int EpiBOD_i;

  int EpiOxy_pr_i;
  
  int KI_SG_i;
  int SGROOT_N_i;
  
  /*
   * tracers
   */
  int NH4_sed_i;
  int NO3_sed_i;

  // NEW STUFF BELOW FOR LEAF UPTAKE OF NUTRIENTS
  int NH4_wc_i;
  int NO3_wc_i;
  int DIP_wc_i;
  int DNO3_i;
  int DPO4_i;
  int DNH4_i;
  int ustrcw_skin_i;
  /*
   * The above stuff is just indices (integers)
   * NEW STUFF ABOVE FOR LEAF UPTAKE OF NUTRIENTS
   */

  int DIP_sed_i;
  int DIC_wc_i;
  int Oxygen_wc_i;
  int Oxy_pr_wc_i;

  /*
   * Multiple sediment layer indicies
   */
  int *NO3_multi_sed_i; 
  int *NH4_multi_sed_i;
  int *DIP_multi_sed_i;
  
  /*
   * common cell variables
   */
  int Tfactor_i;
  int umax_i;
  int resp_i;

} workspace;

void seagrass_spectral_grow_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    char* prm = p->prms->se[0]->s;

    int OFFSET_SED = tracers->n;
    int OFFSET_EPI = tracers->n * 2;

    int topk_sed = ginterface_getsedtopk(e->model, 0);
    int botk_sed = ginterface_getsedbotk(e->model, 0);

    ws->nsed = abs(topk_sed - botk_sed) + 1;

    eco_write_setup(e,"Seagrass growth uses %d sediment layers \n",ws->nsed);

    int k;
    e->use_multi_sed = 1;

    p->workspace = ws;

    ws->species = prm[0];
 
    switch (ws->species){

    case 'Z':   /* Zostera - standard terms */

      eco_write_setup(e,"Seagrass species Zostera chosen \n");

      ws->umax_t0 = get_parameter_value(e, "SGumax");
      ws->omega = get_parameter_value(e, "SGleafden");
      ws->KN = get_parameter_value(e, "SG_KN");
      ws->KP = get_parameter_value(e, "SG_KP");
      ws->SGfrac = get_parameter_value(e, "SGfrac");
      ws->SGtransrate = get_parameter_value(e, "SGtransrate");
      ws->SG_mL = get_parameter_value(e, "SG_mL");

      ws->SGrootdepth = try_parameter_value(e, "SGrootdepth");
      if (isnan(ws->SGrootdepth))
	ws->SGrootdepth = -0.15;

      ws->SGmlr = try_parameter_value(e, "SGmlr");
      if (isnan(ws->SGmlr))
	ws->SGmlr = 4.5;

      ws->SG_N_i = e->find_index(epis, "SG_N", e) + OFFSET_EPI;
      ws->SGROOT_N_i = e->find_index(epis, "SGROOT_N", e) + OFFSET_EPI;
      ws->umax_i = find_index_or_add(e->cv_cell, "SGumax", e);
      ws->resp_i = find_index_or_add(e->cv_cell, "SGresp", e);
      ws->rootbot_i = find_index_or_add(e->cv_cell, "SGrootbot", e);

      ws->SGorient = try_parameter_value(e, "SGorient");
      if (isnan(ws->SGorient))
	ws->SGorient = 1.0;

      /*non essential diagnostic tracers*/

      ws->SG_N_pr_i = e->try_index(epis, "SG_N_pr", e);
      if (ws->SG_N_pr_i > -1) 
	ws->SG_N_pr_i += OFFSET_EPI;

      ws->SG_N_gr_i = e->try_index(epis, "SG_N_gr", e);
      if (ws->SG_N_gr_i > -1) 
	ws->SG_N_gr_i += OFFSET_EPI;

      break;

    case 'H':   /* Halophila */

      /* to avoid changing code, SGH terms are assigned to SG */

      eco_write_setup(e,"Seagrass species Halophila chosen \n");

      ws->umax_t0 = get_parameter_value(e, "SGHumax");
      ws->omega = get_parameter_value(e, "SGHleafden");
      ws->KN = get_parameter_value(e, "SGH_KN");
      ws->KP = get_parameter_value(e, "SGH_KP");
      ws->SGfrac = get_parameter_value(e, "SGHfrac");
      ws->SGtransrate = get_parameter_value(e, "SGHtransrate");
      ws->SG_mL = get_parameter_value(e, "SGH_mL");

      ws->SGrootdepth = try_parameter_value(e, "SGHrootdepth");
      if (isnan(ws->SGrootdepth))
	ws->SGrootdepth = -0.08;

      ws->SGmlr = try_parameter_value(e, "SGHmlr");
      if (isnan(ws->SGmlr))
	ws->SGmlr = 2.8; // Note this is 2.0 in B2p0

      ws->SGorient = try_parameter_value(e, "SGHorient");
      if (isnan(ws->SGorient))
	ws->SGorient = 1.0;

      ws->SG_N_i = e->find_index(epis, "SGH_N", e) + OFFSET_EPI;
      ws->SGROOT_N_i = e->find_index(epis, "SGHROOT_N", e) + OFFSET_EPI;
      ws->umax_i = find_index_or_add(e->cv_cell, "SGHumax", e);
      ws->resp_i = find_index_or_add(e->cv_cell, "SGHresp", e);
      ws->rootbot_i = find_index_or_add(e->cv_cell, "SGHrootbot", e);

      /*non essential diagnostic tracers*/

      ws->SG_N_pr_i = e->try_index(epis, "SGH_N_pr", e);
      if (ws->SG_N_pr_i > -1) 
	ws->SG_N_pr_i += OFFSET_EPI;

      ws->SG_N_gr_i = e->try_index(epis, "SGH_N_gr", e);
      if (ws->SG_N_gr_i > -1) 
	ws->SG_N_gr_i += OFFSET_EPI;
      break;

    case 'D':   /* Deep */

      /* to avoid changing code, SGH terms are assigned to SG */

      eco_write_setup(e,"Deep seagrass species chosen \n");

      ws->umax_t0 = get_parameter_value(e, "SGDumax");
      ws->omega = get_parameter_value(e, "SGDleafden");
      ws->KN = get_parameter_value(e, "SGD_KN");
      ws->KP = get_parameter_value(e, "SGD_KP");
      ws->SGfrac = get_parameter_value(e, "SGDfrac");
      ws->SGtransrate = get_parameter_value(e, "SGDtransrate");
      ws->SG_mL = get_parameter_value(e, "SGD_mL");

      ws->SGrootdepth = try_parameter_value(e, "SGDrootdepth");
      if (isnan(ws->SGrootdepth))
	ws->SGrootdepth = -0.05;

      ws->SGmlr = try_parameter_value(e, "SGDmlr");
      if (isnan(ws->SGmlr))
	ws->SGmlr = 1.5;

      ws->SGorient = try_parameter_value(e, "SGDorient");
      if (isnan(ws->SGorient))
	ws->SGorient = 1.0;

      ws->SG_N_i = e->find_index(epis, "SGD_N", e) + OFFSET_EPI;
      ws->SGROOT_N_i = e->find_index(epis, "SGDROOT_N", e) + OFFSET_EPI;
      ws->umax_i = find_index_or_add(e->cv_cell, "SGDumax", e);
      ws->resp_i = find_index_or_add(e->cv_cell, "SGDresp", e);
      ws->rootbot_i = find_index_or_add(e->cv_cell, "SGDrootbot", e);

      /*non essential diagnostic tracers*/

      ws->SG_N_pr_i = e->try_index(epis, "SGD_N_pr", e);
      if (ws->SG_N_pr_i > -1) 
	ws->SG_N_pr_i += OFFSET_EPI;

      ws->SG_N_gr_i = e->try_index(epis, "SGD_N_gr", e);
      if (ws->SG_N_gr_i > -1) 
	ws->SG_N_gr_i += OFFSET_EPI;
      break;
    case 'P':   /* Posidonia */
 
      /* still to write posidonia code */
      break;
    
    default:
      e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected seagrass type: expected \"Zostera\", \"Halophila\", \"Deep\"or \"Posidonia\"\n", e->processfname, p->name, prm);

    }

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;
    /*
     * tracers
     */
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_sed_i = e->find_index(tracers, "NH4", e) + OFFSET_SED;
    ws->NO3_sed_i = e->find_index(tracers, "NO3", e) + OFFSET_SED;
    ws->DIP_sed_i = e->find_index(tracers, "DIP", e) + OFFSET_SED;
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);

    // NEW STUFF BELOW FOR LEAF UPTAKE OF NUTRIENTS
    ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);
    ws->DIP_wc_i = e->find_index(tracers, "DIP", e);
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e); 
    ws->DNH4_i = find_index(e->cv_cell, "DNH4", e);
    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;
  /*
    * The above stuff populates indices (integers)
    * NEW STUFF ABOVE FOR LEAF UPTAKE OF NUTRIENTS
    */

    /*
     * Allocate arrays for multiple sediment tracer values
     */
    ws->NO3_multi_sed_i = malloc((topk_sed-botk_sed+1)*sizeof(int));
    ws->NH4_multi_sed_i = malloc((topk_sed-botk_sed+1)*sizeof(int));
    ws->DIP_multi_sed_i = malloc((topk_sed-botk_sed+1)*sizeof(int));
    for (k = topk_sed; k >= botk_sed; k--){
      ws->NO3_multi_sed_i[k] = get_multi_sed_index(e, k, topk_sed, "NO3");
      ws->NH4_multi_sed_i[k] = get_multi_sed_index(e, k, topk_sed, "NH4");
      ws->DIP_multi_sed_i[k] = get_multi_sed_index(e, k, topk_sed, "DIP");
    }
   
    /*non essential diagnostic tracer*/

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);
    ws->EpiOxy_pr_i = e->try_index(epis, "EpiOxy_pr", e);
      if (ws->EpiOxy_pr_i > -1) 
	ws->EpiOxy_pr_i  += OFFSET_EPI;

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
}

void seagrass_spectral_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;

    /* In this routine KI_SG is generic, but species-specific cell variables 
       must be added for light_spectral_epi.c */

    switch (ws->species){

    case 'Z':   /* Zostera - standard terms */
      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SG", e);
      break;
    case 'H':   /* Halophila */
      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SGH", e);
      break;
    case 'P':   /* Posidonia */
      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SGP", e);
      break;
    case 'D':
      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SGD", e);
      break;
    }
}

void seagrass_spectral_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->NO3_multi_sed_i);
    free(ws->NH4_multi_sed_i);
    free(ws->DIP_multi_sed_i);
    free(ws);
}

void seagrass_spectral_grow_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    column *col = c->col;
    double* cv = c->cv;
    double* y = c->y;
    double* dz_multi_sed = c->dz_multi_sed;
    double* porosity_multi_sed = c->porosity_multi_sed;
    int k;

    double NO3_sed = 0.0;
    double NH4_sed = 0.0;
    double DIP_sed = 0.0;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    /* Multiply SG_N and SGROOT_N by unitch so that mass balance units are consistent */

    double SG_N = y[ws->SG_N_i] * unitch;
    double SGROOT_N = y[ws->SGROOT_N_i] * unitch;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    /* determine layers, max of sediment layers (i.e. surface). */
       
    cv[ws->rootbot_i] = col->topk_sed;

    double total_sed_depth = 0.0;

    for (k = col->topk_sed; k >= col->botk_sed; k--){
      total_sed_depth += dz_multi_sed[k];
      if (total_sed_depth < (-ws->SGrootdepth)){ 
	cv[ws->rootbot_i] = k-1;
      }
    }
    
    /* if roots go deeper than the sediment layers */
    
    if (cv[ws->rootbot_i] < 0)
      cv[ws->rootbot_i] = 0;

    /* mass balance - add sediment dissolved nutrients for one species, not including top layer */

    if (ws->species == 'Z'){
      for (k = col->topk_sed-1; k >= col->botk_sed; k--){
	NO3_sed += y[ws->NO3_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
	NH4_sed += y[ws->NH4_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
	DIP_sed += y[ws->DIP_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
      }
    }

    if (ws->do_mb) {
      y[ws->EpiTN_i] += SG_N + SGROOT_N + NO3_sed + NH4_sed;
      y[ws->EpiTP_i] += (SG_N + SGROOT_N)* atk_W_P + DIP_sed;
      y[ws->EpiTC_i] += (SG_N + SGROOT_N)* atk_W_C;
      y[ws->EpiBOD_i] += (SG_N + SGROOT_N)* atk_W_O - NO3_sed*48.0/14.01;
    }

    /* Calculated photorespiration term */

    /* Photorespiration at low biomass is ~4.5 mol photon m-2 d-1, 0.7 is the wavelength average absorbance */

    /* Factor of two accounts for mlr being a daily dose, and SG_mL also occuring over 24 h, while resp only 
       acts during daylight hours, so multiple by 2 to get to 24 hours */

    cv[ws->resp_i] = 2.0 * (ws->omega * 0.7 * ws->SGmlr * (30.0/5500.0) * 14.01 / 86400.0 * ws->SGorient - ws->SG_mL);
}

void seagrass_spectral_grow_epi_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  double* y = ia->y;
  double* y1 = ia->y1;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double* dz_multi_sed = c->dz_multi_sed;
  double dz_wc = c->dz_wc;
  double* porosity_multi_sed = c->porosity_multi_sed;
  column *col = c->col;
  ecology *e = col->e;
  int k;
  
  double SG_N = y[ws->SG_N_i];

  if (SG_N > 1e-12){
    
    double SGROOT_N = y[ws->SGROOT_N_i];
    double NO3_multi_sed[ws->nsed];
    double NH4_multi_sed[ws->nsed]; 
    double DIP_multi_sed[ws->nsed];
    
    double NO3_sed = 0.0;
    double NH4_sed = 0.0;
    double DIP_sed = 0.0;
    double NO3frac[ws->nsed];
    double NH4frac[ws->nsed];
    double DIPfrac[ws->nsed];
    double sumdz = 0.0;
    double sumdzNO3 = 0.0;double sumdzNH4 = 0.0;double sumdzDIP = 0.0;
    
    double NO3_water = y[ws->NO3_wc_i];
    double NH4_water = y[ws->NH4_wc_i];
    double DIP_water = y[ws->DIP_wc_i];
    
    int rootbot = cv[ws->rootbot_i];

    /* These are all the sediment layers */

    for (k = col->topk_sed; k >= rootbot; k--) {
      NO3_multi_sed[k] = max(y[ws->NO3_multi_sed_i[k]],EPS_DIN);
      NH4_multi_sed[k] = max(y[ws->NH4_multi_sed_i[k]],EPS_DIN);
      DIP_multi_sed[k] = max(y[ws->DIP_multi_sed_i[k]],EPS_DIP);
      sumdz += dz_multi_sed[k]*porosity_multi_sed[k];
      if (dz_multi_sed[k] > 0.001){
	sumdzNO3 += dz_multi_sed[k]*porosity_multi_sed[k]*NO3_multi_sed[k];
	sumdzNH4 += dz_multi_sed[k]*porosity_multi_sed[k]*NH4_multi_sed[k];
	sumdzDIP += dz_multi_sed[k]*porosity_multi_sed[k]*DIP_multi_sed[k];
      }
    }

    /* Calculate a single "effective" nutrient concentration. */

    for (k = col->topk_sed; k >= rootbot; k--) {
      NO3_sed += NO3_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdz;
      NH4_sed += NH4_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdz;
      DIP_sed += DIP_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdz;
    }
   
    /* Calculated fraction taken from each layer */

    for (k = col->topk_sed; k >= rootbot; k--) {
      if (dz_multi_sed[k] > 0.001){
	NO3frac[k] = NO3_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdzNO3;
	NH4frac[k] = NH4_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdzNH4;
	DIPfrac[k] = DIP_multi_sed[k]*dz_multi_sed[k]*porosity_multi_sed[k]/sumdzDIP;
      }else{
	NO3frac[k] = 0.0;
	NH4frac[k] = 0.0;
	DIPfrac[k] = 0.0; 
      }
    }

    double umax = cv[ws->umax_i]; /* s-1 */

    /* kI is calculated in light_spectral_uq_epi.c in mol photon m-2 s-1 and converted to turnover time */

    /* SG_N now in g N m-2, so need to multiply mgN2molN by 1000 */

    double kI = c->cv[ws->KI_SG_i] * atk_A_N / ( atk_A_I * SG_N * mgN2molN * 1000.0);

    kI = max(0.0,kI - cv[ws->resp_i]);  // s-1

    if (kI > 1.0e-9){  // avoid divide by zeros.

      /* CALCULATING kN FOR WATER (FOLLOWING MACROALGAE CODE) BELOW */
      double Sc[3];
      Sc[0] = 1.05e-6 / c->cv[ws->DNO3_i];
      Sc[1] = 1.05e-6 / c->cv[ws->DNH4_i];
      Sc[2] = 1.05e-6 / c->cv[ws->DPO4_i];
      
      double tau = y[ws->ustrcw_skin_i]*y[ws->ustrcw_skin_i];
      double S_NO3 = 2850.0*pow((2.0*tau),0.38)*pow(Sc[0],-0.6)/86400.0; /* m s-1 */
      double S_NH4 = 2850.0*pow((2.0*tau),0.38)*pow(Sc[1],-0.6)/86400.0; /* m s-1 */
      double S_DIP = 2850.0*pow((2.0*tau),0.38)*pow(Sc[2],-0.6)/86400.0; /* m s-1 */
      double SA = 1.0-exp(-SG_N * ws->omega); // dimensionless
      
      // The above doesn't take into account bending angle or absorbance (the latter being wavelength-dependent).
      
      // Supply in units of turnover of biomass //
      
      double kN_NO3_water = S_NO3 * SA * NO3_water / (SG_N * mgN2molN * 1000.0); /* s-1 */
      double kN_NH4_water = S_NH4 * SA * NH4_water / (SG_N * mgN2molN * 1000.0); /* s-1 */
      double kP_water = S_DIP * SA * DIP_water * (atk_A_N / atk_A_P ) / (SG_N * mgN2molN * 1000.0); /*  s-1 */
      
      /* determined what growth is, then work to attribute it to different sources */
      
      double maxNuptake = umax * NH4_sed / ( ws->KN + NH4_sed) + umax * NO3_sed / ( ws->KN + NO3_sed) + kN_NO3_water + kN_NH4_water;
      double maxPuptake = umax * DIP_sed / ( ws->KP + DIP_sed) + kP_water;
      
      double growthrate = min(umax,min(kI,min(maxNuptake,maxPuptake)));
      
      /* Assume sediment nutrients used first, with NH4 before NO3. */
      
      double kNH4_sed = min(growthrate,umax * NH4_sed / ( ws->KN + NH4_sed));
      double kNO3_sed = min(growthrate - kNH4_sed, umax * NO3_sed / ( ws->KN + NO3_sed));
      double kP_sed = min(growthrate,umax * DIP_sed / ( ws->KP + DIP_sed));
      
      double kN_sed = kNH4_sed + kNO3_sed;
      
      /* Limit water column uptake so it doesn't exceed maximum growth rate */
      
      kP_water = min(kP_water,growthrate - kP_sed);
      kN_NH4_water = min(kN_NH4_water, growthrate - kN_sed);
      kN_NO3_water = min(kN_NO3_water, growthrate - kN_sed - kN_NH4_water);
      
      // Law of the minimum, but including maximum growth rate.
      
      double kN_water = kN_NO3_water + kN_NH4_water;
      double kN = kN_sed + kN_water;
      double kP = kP_sed + kP_water;           
      double growth = SG_N * growthrate;
      
      if (kN > 1.0e-12){   // kN is s-1 

	y1[ws->SG_N_i] += growth;
	
	/* Now need to take out of the layers based on their thickness / conc. */
	
	/* PARTITION NUTRIENT UPTAKE BELOW */
	
	double NH4uptake_frac_from_water = kN_NH4_water/kN;
	double NO3uptake_frac_from_water = kN_NO3_water/kN;
	double DIPuptake_frac_from_water = kP_water/kP;
	
	// Assume even NO3 taken up in sediments releases O2 in water column - i.e. comes up xylem as NO3
	
	double Oxy_pr = (growth * unitch * atk_W_O + (growth * NO3uptake_frac_from_water + SG_N * kNO3_sed) * unitch * 48.0/14.01)/ dz_wc;
	
	/* PARTITION NUTRIENT UPTAKE ABOVE */
	
	for (k = col->topk_sed; k >= rootbot; k--) {
	  y1[ws->NH4_multi_sed_i[k]] -= unitch * SG_N * kNH4_sed * NH4frac[k] / dz_multi_sed[k] / porosity_multi_sed[k];
	  y1[ws->NO3_multi_sed_i[k]] -= unitch * SG_N * kNO3_sed * NO3frac[k] / dz_multi_sed[k] / porosity_multi_sed[k];
	  y1[ws->DIP_multi_sed_i[k]] -= unitch * growth * (1.0-DIPuptake_frac_from_water) * atk_W_P * DIPfrac[k] / dz_multi_sed[k] / porosity_multi_sed[k];
	}
	
	y1[ws->DIC_wc_i] -= growth * unitch * atk_W_C / dz_wc;
	
	/* ADDED DERIVATIVE FOR WC NUTRIENTS BELOW */

	y1[ws->NH4_wc_i] -= growth * NH4uptake_frac_from_water * unitch / dz_wc;
	y1[ws->NO3_wc_i] -= growth * NO3uptake_frac_from_water * unitch / dz_wc;
	y1[ws->DIP_wc_i] -= growth * atk_W_P * DIPuptake_frac_from_water * unitch / dz_wc;

	/* ADDED DERIVATIVE FOR WC NUTRIENTS ABOVE */
	
	y1[ws->Oxygen_wc_i] += Oxy_pr;
	
	if (ws-> SG_N_gr_i > -1)
	  y1[ws->SG_N_gr_i] += growthrate * SEC_PER_DAY; /* d-1 */
	if (ws-> SG_N_pr_i > -1)
	  y1[ws->SG_N_pr_i] += growth * SEC_PER_DAY * atk_W_C; /* gC m-2 d-1 */
	if (ws->Oxy_pr_wc_i > -1)
	  y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY ;  /*mgO m-3 d-1 */
	if (ws->EpiOxy_pr_i > -1)
	  y1[ws->EpiOxy_pr_i] += Oxy_pr * SEC_PER_DAY * dz_wc;  /*mgO m-2 d-1 */
	
      } // if growth > 0.0
    }  // if kI < 0.0
    
    /* Translocation of N between root and leaf - simplfy version in paper */ 
    
    double SGfrac = ws->SGfrac;
    double SGtransrate = ws->SGtransrate;
    
    double translocate_to_roots = (SGfrac * (SG_N+SGROOT_N) - SGROOT_N) * SGtransrate;
    
    y1[ws->SGROOT_N_i] += translocate_to_roots;
    y1[ws->SG_N_i] -= translocate_to_roots;
    
    }
  }
  
void seagrass_spectral_grow_epi_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    column *col = c->col;
    double* y = c->y;
    double* dz_multi_sed = c->dz_multi_sed;
    double* porosity_multi_sed = c->porosity_multi_sed;
    int k;

    double* cv = c->cv;

    double NO3_sed = 0.0;
    double NH4_sed = 0.0;
    double DIP_sed = 0.0;

    double SG_N = y[ws->SG_N_i] * unitch;
    double SGROOT_N = y[ws->SGROOT_N_i] * unitch;

    /* mass balance - only count sediment change for one seagrass species. */
    if (ws->species == 'Z'){
      for (k = col->topk_sed-1; k >= col->botk_sed; k--){
	NO3_sed += y[ws->NO3_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
	NH4_sed += y[ws->NH4_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
	DIP_sed += y[ws->DIP_multi_sed_i[k]]*dz_multi_sed[k]*porosity_multi_sed[k];
      }
    }

    y[ws->EpiTN_i] += SG_N + SGROOT_N + NO3_sed + NH4_sed;
    y[ws->EpiTP_i] += (SG_N + SGROOT_N)* atk_W_P + DIP_sed;
    y[ws->EpiTC_i] += (SG_N + SGROOT_N)* atk_W_C;
    y[ws->EpiBOD_i] += (SG_N + SGROOT_N)* atk_W_O - NO3_sed*48.0/14.01;
}
