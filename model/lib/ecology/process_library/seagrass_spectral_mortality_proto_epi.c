/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_spectral_mortality_proto_epi.c
 *  
 *  Description: Seagrass mortality including linear mortality, quadratic mortality 
 *               and shear stress mortality, affecting both above and below ground biomass.
 *  
 *  Options: seagrass_spectral_mortality_proto_epi(Z|H|D|P)
 *
 *  Z - Zostera (SG)
 *  H - Halophila (SGH)
 *  D - Deep Halophila (SGD)
 *  P - Posidonia (SGP)
 * 
 *  Model described in: Baird, M. E., M. P. Adams, R. C. Babcock, K. Oubelkheir, M. Mongin, 
 *                      K. A. Wild-Allen, J. Skerratt, B. J. Robson, K. Petrou, P. J. Ralph, 
 *                      K. R. O'Brien, A. B. Carter, J. C. Jarvis, M. A. Rasheed (2016) 
 *                      A biophysical representation of seagrass growth for application in a 
 *                      complex shallow-water biogeochemical model Ecol. Mod. 325: 13-27. 
 *    
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: seagrass_spectral_mortality_proto_epi.c 5942 2018-09-13 04:28:00Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "seagrass_spectral_mortality_proto_epi.h"

typedef struct {
  char species;
  
  /*
   * parameters
   */
  
  double mL_t0;
  double mL_ROOT_t0;
  double SGleafden;
  double SGfrac;
  double SG_mQ;
  
  double seed;

  double SG_tau_critical;
  double SG_tau_efold;
  
  /*
   * epis
   */

  int SGROOT_N_i;
  int SG_N_i;
  
  /*
   * tracers
   */
  int DetBL_N_sed_i;
  
  /*
   * common cell variables
   */
  int Tfactor_i;
  int mL_i;
  int mL_ROOT_i;

  int ustrcw_skin_i;
  int SG_shear_mort_i;


} workspace;

void seagrass_spectral_mortality_proto_epi_init(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = malloc(sizeof(workspace));
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;

  char* prm = p->prms->se[0]->s;
  
  int OFFSET_SED = tracers->n;
  int OFFSET_EPI = tracers->n * 2;
  
  p->workspace = ws;

  ws->species = prm[0];
  
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->DetBL_N_sed_i = e->find_index(tracers, "DetBL_N", e) + OFFSET_SED;
  
  // NEW STUFF FOR SHEAR STRESS-DEPENDENT MORTALITY
  ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

  ws->SG_mQ = 0.0;

  switch (ws->species){

  case 'Z': /* Zostera - standard terms */

    eco_write_setup(e,"Chosen Zostera in spectral mort photo. \n");

    ws->mL_ROOT_t0 = get_parameter_value(e, "SGROOT_mL");
    ws->mL_t0 = get_parameter_value(e, "SG_mL");
    ws->SGleafden = get_parameter_value(e, "SGleafden");
    ws->mL_i = find_index_or_add(e->cv_cell, "SG_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGROOT_mL", e);
    
    ws->seed = try_parameter_value(e, "SGseedfrac");
    if (isnan(ws->seed)){
      ws->seed = 0.0;
      eco_write_setup(e,"Code default: SGseedfrac = %e \n",ws->seed);
    }

    ws->SGfrac = get_parameter_value(e, "SGfrac");
    
    ws->SG_tau_critical = try_parameter_value(e, "SG_tau_critical");
    if (isnan(ws->SG_tau_critical)){
      ws->SG_tau_critical = 1.0;
       eco_write_setup(e,"Code default:  SG_tau_critical = %e \n",ws->SG_tau_critical);
    }

    ws->SG_tau_efold = try_parameter_value(e, "SG_tau_efold");
    if (isnan(ws->SG_tau_efold)){
      ws->SG_tau_efold = 0.5 * 24.0 * 3600.0;
      eco_write_setup(e,"Code default: SG_tau_efold = %e \n",ws->SG_tau_efold);
    }

    ws->SG_mQ = try_parameter_value(e, "SG_mQ");
    if (isnan(ws->SG_mQ)){
      ws->SG_mQ = 0.0;
      eco_write_setup(e,"Code default:  SG_mQ  = %e \n", ws->SG_mQ);
    }

    ws->SG_N_i = e->find_index(epis, "SG_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGROOT_N", e) + OFFSET_EPI;
    ws->SG_shear_mort_i = e->try_index(epis, "SG_shear_mort", e);
    if (ws->SG_shear_mort_i > -1)
      ws->SG_shear_mort_i += OFFSET_EPI;
    break;
    
  case 'H':  /* Halophila */

    eco_write_setup(e,"Chosen Halophila in spectral mort photo. \n");
    
    ws->mL_ROOT_t0 = get_parameter_value(e, "SGHROOT_mL");  
    ws->mL_t0 = get_parameter_value(e, "SGH_mL");
    ws->SGleafden = get_parameter_value(e, "SGHleafden");
    ws->mL_i = find_index_or_add(e->cv_cell, "SGH_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGHROOT_mL", e);
    
    ws->seed = try_parameter_value(e, "SGHseedfrac");
    if (isnan(ws->seed)){
      ws->seed = 0.0;
      eco_write_setup(e,"Code default: SGHseedfrac = %e \n",ws->seed);
    }
    ws->SGfrac = get_parameter_value(e, "SGHfrac");
    
    ws->SG_tau_critical = try_parameter_value(e, "SGH_tau_critical");
    if (isnan(ws->SG_tau_critical)){
      ws->SG_tau_critical = 1.0;
      eco_write_setup(e,"Code default:  SGH_tau_critical = %e \n",ws->SG_tau_critical);
    }

    ws->SG_tau_efold = try_parameter_value(e, "SGH_tau_efold");
    if (isnan(ws->SG_tau_efold)){
      ws->SG_tau_efold = 0.5 * 24.0 * 3600.0;
      eco_write_setup(e,"Code default: SGH_tau_efold = %e \n",ws->SG_tau_efold);
    }
    ws->SG_mQ = try_parameter_value(e, "SGH_mQ");
    if (isnan(ws->SG_mQ)){
      ws->SG_mQ = 0.0;
      eco_write_setup(e,"Code default:  SGH_mQ  = %e \n", ws->SG_mQ);
    }

    ws->SG_N_i = e->find_index(epis, "SGH_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGHROOT_N", e) + OFFSET_EPI;
    ws->SG_shear_mort_i = e->try_index(epis, "SGH_shear_mort", e);
    if (ws->SG_shear_mort_i > -1)
      ws->SG_shear_mort_i += OFFSET_EPI;
    break;

  case 'P':   /* Posidonia */

    eco_write_setup(e,"Chosen Posidonia in spectral mort photo. \n");

    ws->mL_ROOT_t0 = get_parameter_value(e, "SGPROOT_mL");  
    ws->mL_t0 = get_parameter_value(e, "SGP_mL");
    ws->SGleafden = get_parameter_value(e, "SGPleafden");
    ws->mL_i = find_index_or_add(e->cv_cell, "SGP_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGPROOT_mL", e);
    ws->seed = try_parameter_value(e, "SGPseedfrac");
    if (isnan(ws->seed))
      ws->seed = 0.0;
    ws->SGfrac = get_parameter_value(e, "SGPfrac");
    
    ws->SG_tau_critical = try_parameter_value(e, "SGP_tau_critical");
    if (isnan(ws->SG_tau_critical))
      ws->SG_tau_critical = 1.0;
    
    ws->SG_tau_efold = try_parameter_value(e, "SGP_tau_efold");
    if (isnan(ws->SG_tau_efold))
      ws->SG_tau_efold = 0.5*24.0*3600.0;

    ws->SG_mQ = try_parameter_value(e, "SGP_mQ");
    if (isnan(ws->SG_mQ))
      ws->SG_mQ = 0.0;
    
    ws->SG_N_i = e->find_index(epis, "SGP_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGPROOT_N", e) + OFFSET_EPI;
    ws->SG_shear_mort_i = e->try_index(epis, "SGP_shear_mort", e);
    if (ws->SG_shear_mort_i > -1)
      ws->SG_shear_mort_i += OFFSET_EPI;
    break;

  case 'D':  /* Deep - using Halophila values for the moment */

    eco_write_setup(e,"Chosen Deep in spectral mort photo. \n");

    ws->mL_ROOT_t0 = get_parameter_value(e, "SGDROOT_mL");  
    ws->mL_t0 = get_parameter_value(e, "SGD_mL");
    ws->SGleafden = get_parameter_value(e, "SGDleafden");
    ws->mL_i = find_index_or_add(e->cv_cell, "SGD_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGDROOT_mL", e);
    
    ws->seed = try_parameter_value(e, "SGDseedfrac");
    if (isnan(ws->seed)){
      ws->seed = 0.0;
      eco_write_setup(e,"Code default: SGDseedfrac = %e \n",ws->seed);
    }

    ws->SGfrac = get_parameter_value(e, "SGDfrac");
    
    ws->SG_tau_critical = try_parameter_value(e, "SGD_tau_critical");
    if (isnan(ws->SG_tau_critical)){
      ws->SG_tau_critical = 1.0;
      eco_write_setup(e,"Code default:  SGD_tau_critical = %e \n",ws->SG_tau_critical);
    }
    
    ws->SG_tau_efold = try_parameter_value(e, "SGD_tau_efold");
    if (isnan(ws->SG_tau_efold)){
      ws->SG_tau_efold = 0.5*24.0*3600.0;
      eco_write_setup(e,"Code default: SGD_tau_efold = %e \n",ws->SG_tau_efold);
    }

    ws->SG_mQ = try_parameter_value(e, "SGD_mQ");
    if (isnan(ws->SG_mQ)){
      ws->SG_mQ = (0.1/86400.0)/0.02;
      eco_write_setup(e,"Code default:  SGH_mQ  = %e \n", ws->SG_mQ);
    }

    /* overwrite default if recom for UQ group */

    if (process_present(e,PT_WC,"recom_extras")){
      ws->SG_mQ = (0.1/86400.0)/0.02;
      eco_write_setup(e,"Written over because in RECOM for UQ:  SGH_mQ  = %e \n", ws->SG_mQ);
    }

    ws->SG_N_i = e->find_index(epis, "SGD_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGDROOT_N", e) + OFFSET_EPI;
    ws->SG_shear_mort_i = e->try_index(epis, "SGD_shear_mort", e);
    if (ws->SG_shear_mort_i > -1)
      ws->SG_shear_mort_i += OFFSET_EPI;
    break;
  default:
    e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected seagrass type: expected \"Zostera\", \"Halophila\" or \"Posidonia\"\n", e->processfname, p->name, prm);
  }    
}

void seagrass_spectral_mortality_proto_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void seagrass_spectral_mortality_proto_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;
    cv[ws->mL_ROOT_i] = ws->mL_ROOT_t0 * Tfactor;

}

void seagrass_spectral_mortality_proto_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;

    /* leave seed population of 1 % cover */

    double seed = ws->seed / ws->SGleafden;

    double SG_N = max(0.0,y[ws->SG_N_i] - seed * (1.0 - ws->SGfrac));
    double SGROOT_N = max(0.0,y[ws->SGROOT_N_i] - seed * ws->SGfrac);
    double mortality = cv[ws->mL_i] * SG_N;
    double mort_root = cv[ws->mL_ROOT_i] * SGROOT_N;

    // ADDITIONAL CODE FOR SHEAR STRESS-DEPENDENT MORTALITY

    double tau = y[ws->ustrcw_skin_i]*y[ws->ustrcw_skin_i]*1000.0;  // density-specific shear stress x density.

    // max = 10 d-1, so as not to make the ODE system too stiff.

    double shear_stress_mort = min(2.0/86400,max((tau - ws->SG_tau_critical)/ws->SG_tau_critical, 0.0) * (1.0 / ws->SG_tau_efold));

    double shear_stress_mortality_A = shear_stress_mort * SG_N;
    double shear_stress_mortality_B = shear_stress_mort * SGROOT_N;

    // ADDITIONAL QUADRATIC TERM FOR DEEP SEAGRASS //

    mortality = mortality + ws->SG_mQ * SG_N * SG_N;   // quadratic term 
 
    // CODE BELOW CHANGED TO ACCOUNT FOR SHEAR STRESS-DEPENDENT MORTALITY
    y1[ws->SGROOT_N_i] -= (mort_root + shear_stress_mortality_B);
    y1[ws->SG_N_i] -= (mortality + shear_stress_mortality_A);
    y1[ws->DetBL_N_sed_i] += 1000.0 * (mortality + shear_stress_mortality_A + mort_root + shear_stress_mortality_B ) / c->dz_sed;

    if (ws->SG_shear_mort_i > -1){
      y1[ws->SG_shear_mort_i] += shear_stress_mortality_A * SEC_PER_DAY;
    }
}

void seagrass_spectral_mortality_proto_epi_postcalc(eprocess* p, void* pp)
{
}
