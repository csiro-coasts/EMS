/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/massbalance_epi.c
 *  
 *  Description:
 *  Mass balance for C, N, P and O in epibenthic.
 *
 *  1. Because seagrass take nutrients from multiple sediment layers, TN, TP and TC and BOD 
 *     include mass from lower sediment layers. 
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: massbalance_epi.c 6293 2019-08-12 02:50:06Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "utils.h"
#include "eprocess.h"
#include "column.h"
#include "cell.h"
#include "einterface.h"
#include "stringtable.h"
#include "constants.h"
#include "massbalance_epi.h"

/* definition of MASSBALANCE_EPS to common header file ecology/constants.h */

typedef struct {
  /*
   * tracers
   */
  int TN_wc_i;
  int TP_wc_i;
  int TC_wc_i;
  int BOD_wc_i;
  int COD_wc_i;
  int Nfix_wc_i;
  int Den_fl_wc_i;
  int Amm_fl_wc_i;

  int TN_sed_i;
  int TP_sed_i;
  int TC_sed_i;
  int BOD_sed_i;
  int COD_sed_i;
  int Den_fl_sed_i;
  int Amm_fl_sed_i;
  int TN_epi_i;
  int TP_epi_i;
  int TC_epi_i;
  int BOD_epi_i;
  int Gnet_i;
  int CO2_flux_i;
  int O2_flux_i;
  int Oxygen_wc_i;
  int Oxygen_sed_i;

  int NO3_wc_i;
  int NO3_sed_i;

  int mb_oxygen;
  
  /*
   * common cell variables
   */
  int TN_old_i;
  int TP_old_i;
  int TC_old_i;
  int TO_old_i;
} workspace;

void massbalance_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;
    workspace* ws = malloc(sizeof(workspace));

    int OFFSET_SED = e->ntr;
    int OFFSET_EPI = e->ntr * 2;

    p->workspace = ws;

    /*
     * tracers
     */
    ws->TN_wc_i = e->find_index(tracers, "TN", e);
    ws->TP_wc_i = e->find_index(tracers, "TP", e);
    ws->TC_wc_i = e->find_index(tracers, "TC", e);
    ws->BOD_wc_i = e->find_index(tracers, "BOD", e);
    ws->COD_wc_i = e->try_index(tracers, "COD", e);
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->Nfix_wc_i = e->try_index(tracers, "Nfix", e);
    ws->CO2_flux_i = e->try_index(tracers, "CO2_flux", e);
    ws->O2_flux_i = e->try_index(tracers, "O2_flux", e);

    ws->Den_fl_wc_i = e->find_index(tracers, "Den_fl", e);
    ws->Amm_fl_wc_i = e->try_index(tracers, "Amm_fl", e);

    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);

    ws->TN_sed_i = e->find_index(tracers, "TN", e) + OFFSET_SED;
    ws->TP_sed_i = e->find_index(tracers, "TP", e) + OFFSET_SED;
    ws->TC_sed_i = e->find_index(tracers, "TC", e) + OFFSET_SED;
    ws->BOD_sed_i = e->try_index(tracers, "BOD", e) + OFFSET_SED;
    ws->COD_sed_i = e->try_index(tracers, "COD", e) + OFFSET_SED;
    ws->Oxygen_sed_i = e->find_index(tracers, "Oxygen", e) + OFFSET_SED;

    ws->Den_fl_sed_i = e->find_index(tracers, "Den_fl", e) + OFFSET_SED;

    ws->NO3_sed_i = e->find_index(tracers, "NO3", e) + OFFSET_SED;

    ws->Amm_fl_sed_i = e->try_index(tracers, "Amm_fl", e);
    if (ws->Amm_fl_sed_i > -1)
	     ws->Amm_fl_sed_i += OFFSET_SED;
    /*
     * epis
     */

    ws->Gnet_i = e->try_index(epis, "Gnet", e);
    if (ws->Gnet_i > -1) 
	 ws->Gnet_i  += OFFSET_EPI;

    ws->TN_epi_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->TP_epi_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->TC_epi_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;

    ws->mb_oxygen = 0;

    if (ws->COD_wc_i > -1){
      ws->mb_oxygen = 1;
      ws->BOD_epi_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;
      ws->TO_old_i = find_index_or_add(e->cv_cell, "TO_old", e);
    }
    
    ws->TN_old_i = find_index_or_add(e->cv_cell, "TN_old", e);
    ws->TP_old_i = find_index_or_add(e->cv_cell, "TP_old", e);
    ws->TC_old_i = find_index_or_add(e->cv_cell, "TC_old", e);

    /*
     * set a flag indicating doing mass balance calculations
     */
    stringtable_add_ifabscent(e->cv_model, "massbalance_wc", -1);
    stringtable_add_ifabscent(e->cv_model, "massbalance_epi", -1);
    stringtable_add_ifabscent(e->cv_model, "massbalance_sed", -1);
    eco_write_setup(e,"\nMass balance in epibenthic to fractional difference of %e \n",MASSBALANCE_EPS*1000.0);
}

void massbalance_epi_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
  
  if (process_present(e,PT_EPI,"coral_spectral_grow_epi") && (ws->Gnet_i < 0)){
    emstag(LPANIC,"eco:massbalance_epi:init","Mass balance in epi not possible with corals calcifying and Gnet not in the tracer list. Put Gnet in tracer list or remove massbalance_epi.");
  }
  
  if (process_present(e,PT_EPI,"coral_spectral_carb_epi") && (ws->Gnet_i < 0)){
    emstag(LPANIC,"eco:massbalance_epi:init","Mass balance in epi not possible with corals calcifying and Gnet not in the tracer list. Put Gnet in tracer list or remove massbalance_epi.");
  }
}

void massbalance_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void massbalance_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
    double dz_wc = c->dz_wc;
    double dz_sed = c->dz_sed;
    double porosity = c->porosity;

    cv[ws->TN_old_i] = y[ws->TN_wc_i] * dz_wc + y[ws->TN_sed_i] * dz_sed + y[ws->TN_epi_i];
    cv[ws->TP_old_i] = y[ws->TP_wc_i] * dz_wc + y[ws->TP_sed_i] * dz_sed + y[ws->TP_epi_i];
    cv[ws->TC_old_i] = y[ws->TC_wc_i] * dz_wc + y[ws->TC_sed_i] * dz_sed + y[ws->TC_epi_i];

    if (ws->COD_wc_i > -1){

      cv[ws->TO_old_i] = (y[ws->Oxygen_wc_i] - y[ws->COD_wc_i] - y[ws->BOD_wc_i] + y[ws->NO3_wc_i] / 14.01 * 48.0 ) * dz_wc;
      
      cv[ws->TO_old_i] += (y[ws->Oxygen_sed_i] - y[ws->COD_sed_i] + y[ws->NO3_sed_i] / 14.01 * 48.0) * dz_sed * porosity - y[ws->BOD_sed_i] * dz_sed; 
      
      cv[ws->TO_old_i] += - y[ws->BOD_epi_i];

      y[ws->BOD_wc_i] = 0.0;
      y[ws->BOD_sed_i] = 0.0;
      y[ws->BOD_epi_i] = 0.0;
    }

    y[ws->TN_wc_i] = 0.0;
    y[ws->TP_wc_i] = 0.0;
    y[ws->TC_wc_i] = 0.0;
    y[ws->TN_sed_i] = 0.0;
    y[ws->TP_sed_i] = 0.0;
    y[ws->TC_sed_i] = 0.0;
    y[ws->TN_epi_i] = 0.0;
    y[ws->TP_epi_i] = 0.0;
    y[ws->TC_epi_i] = 0.0;    

}

void massbalance_epi_postcalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
    double dz_wc = c->dz_wc;
    double dz_sed = c->dz_sed;
    double Nfix_wc = (ws->Nfix_wc_i >= 0) ? y[ws->Nfix_wc_i] : 0.0;
    double CO2_flux = (ws->CO2_flux_i >= 0) ? y[ws->CO2_flux_i] : 0.0;

    double Gnet = (ws->Gnet_i >= 0) ? y[ws->Gnet_i] : 0.0;

    int ij[2];

    double porosity = c->porosity;

    double TO = (y[ws->Oxygen_wc_i] - y[ws->COD_wc_i] - y[ws->BOD_wc_i] + y[ws->NO3_wc_i] / 14.01 * 48.0 ) * dz_wc;

    TO += (y[ws->Oxygen_sed_i] - y[ws->COD_sed_i] + y[ws->NO3_sed_i] / 14.01 * 48.0) * dz_sed * porosity - y[ws->BOD_sed_i] * dz_sed; 

    TO += - y[ws->BOD_epi_i];

    double TN = (y[ws->TN_wc_i] - Nfix_wc) * dz_wc + y[ws->TN_sed_i] * dz_sed + y[ws->Den_fl_sed_i] / SEC_PER_DAY + y[ws->TN_epi_i];

    if (ws->Amm_fl_sed_i > -1)
      TN += (y[ws->Amm_fl_sed_i] + y[ws->Amm_fl_wc_i]) / SEC_PER_DAY;

    TN += y[ws->Den_fl_wc_i] / SEC_PER_DAY ;

    double TP = y[ws->TP_wc_i] * dz_wc + y[ws->TP_sed_i] * dz_sed + y[ws->TP_epi_i] + 0.0;
    double TC = y[ws->TC_wc_i] * dz_wc + y[ws->TC_sed_i] * dz_sed + y[ws->TC_epi_i] + Gnet;

    double eps = fabs(TN - cv[ws->TN_old_i]) / (TN + cv[ws->TN_old_i]);

    if (eps > MASSBALANCE_EPS*1000.0){
      einterface_get_ij(c->e->model,c->b,ij);
      e_warn("ecology: error: N balance (%e,%e) violation in epibenthic cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d <%u %u>\n", TN,cv[ws->TN_old_i], eps, e->nstep, c->nsubstep, c->col->b,(p->type == PT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
        }
    eps = fabs(TP - cv[ws->TP_old_i]) / (TP + cv[ws->TP_old_i]);

    if (eps > MASSBALANCE_EPS*1000.0){
      einterface_get_ij(c->e->model,c->b,ij);
      e_warn("ecology: error: P balance (%e,%e) violation in epibenthic cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d <%u %u>\n", TP,cv[ws->TP_old_i], eps, e->nstep, c->nsubstep, c->col->b,(p->type == PT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
    }
    if (CO2_flux == 0.0){

      /* Can't do mass balance on epi when there is only one water column layer as it does not know that it is having only called massbalance_wc once */
  
    eps = fabs(TC - cv[ws->TC_old_i]) / (TC + cv[ws->TC_old_i]);
      if (eps > MASSBALANCE_EPS*1000.0)
        e_warn("ecology: error: C balance (%e,%e) violation in epibenthic cell by %.3g, nstep = %d, nsubstep = %d, b = %d\n", TC,cv[ws->TC_old_i],eps, e->nstep, c->nsubstep, c->col->b);
      }

    if (ws->mb_oxygen){

      if (y[ws->O2_flux_i] == 0.0){

      /* because TO can be close to zero */

	eps = fabs(TO - cv[ws->TO_old_i]) / max(fabs(TO + cv[ws->TO_old_i]),8000.0);

	if (eps > MASSBALANCE_EPS*1000.0){
	  e_warn("ecology: error: O2 balance (%e,%e) violation in epibenthic cell by %.3g, nstep = %d, nsubstep = %d, b = %d\n", TO,cv[ws->TO_old_i],eps, e->nstep, c->nsubstep, c->col->b);
	}
      }
    }
}
