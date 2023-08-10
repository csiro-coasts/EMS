/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/massbalance_wc.c
 *  
 *  Description:
 *  Process implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: massbalance_wc.c 7253 2022-10-26 10:19:01Z bai155 $
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
#include "massbalance_wc.h"

/*UR definition of MASSBALANCE_EPI to common header file eco_constants.h
*/

typedef struct {
    /*
     * tracers
     */
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;
  int COD_i;
  int Oxygen_i;
  int Nfix_i;
  int CO2_flux_i;
  int O2_flux_i;
  int COD_flux_i;

  int NO3_i;
  int Den_fl_i;
  int Amm_fl_i;

    /*
     * common cell variables
     */

  int TN_old_i;
  int TP_old_i;
  int TC_old_i;
  int TO_old_i;

} workspace;

void massbalance_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * tracers
     */
    ws->BOD_i = -1;
    ws->COD_i = -1;
    ws->Nfix_i = -1;
    ws->CO2_flux_i = -1;
    ws->O2_flux_i = -1;
    ws->COD_flux_i = -1;
    ws->Den_fl_i = -1;
    ws->Amm_fl_i = -1;

    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->BOD_i = e->try_index(tracers, "BOD", e);
    ws->COD_i = e->try_index(tracers, "COD", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);

    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->Nfix_i = e->try_index(tracers, "Nfix", e);

    ws->CO2_flux_i = e->try_index(tracers, "CO2_flux", e);
    ws->O2_flux_i = e->try_index(tracers, "O2_flux", e);
    ws->COD_flux_i = e->try_index(tracers, "COD_flux", e);
    
    ws->Den_fl_i = e->find_index(tracers, "Den_fl", e);
    ws->Amm_fl_i = e->try_index(tracers, "Amm_fl", e);
   
    /*
     * common cell variables
     */
    ws->TN_old_i = find_index_or_add(e->cv_cell, "TN_old", e);
    ws->TP_old_i = find_index_or_add(e->cv_cell, "TP_old", e);
    ws->TC_old_i = find_index_or_add(e->cv_cell, "TC_old", e);

    ws->TO_old_i = -1;

    if (ws->COD_i > -1 && ws->O2_flux_i > -1 && ws->BOD_i > -1 && ws->COD_flux_i > -1){
      ws->TO_old_i = find_index_or_add(e->cv_cell, "TO_old", e);
    }else{
      eco_write_setup(e,"\nNot doing oxygen balance because one of COD, O2_flux, BOD or COD_flux is not in tracer list\n");
      }
    /*
     * set a flag indicating doing mass balance calculations
     */
    stringtable_add_ifabscent(e->cv_model, "massbalance_wc", -1);

    eco_write_setup(e,"\nMass balance in water column to fractional difference of %e \n",MASSBALANCE_EPS);

}

void massbalance_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void massbalance_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    /*
     * if this is a child cell, return, the parent should take care of
     * checking the mass balance
     */
    if (c->parent != NULL)
        return;

    cv[ws->TN_old_i] = y[ws->TN_i];
    cv[ws->TP_old_i] = y[ws->TP_i];
    cv[ws->TC_old_i] = y[ws->TC_i];

    if (ws->TO_old_i > -1){
      cv[ws->TO_old_i] = y[ws->Oxygen_i] - y[ws->BOD_i] - y[ws->COD_i] + y[ws->NO3_i] / 14.01 * 48.0;
      y[ws->BOD_i] = 0.0;
    }

    y[ws->TN_i] = 0.0;
    y[ws->TP_i] = 0.0;
    y[ws->TC_i] = 0.0;
   
}

void massbalance_wc_postcalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y;
    double TN, Nfix, TP, TC, TO, TO_old, eps, CO2_flux, O2_flux, COD_flux;

    /*
     * if this is a child cell, return, the parent should take care of
     * checking the mass balance
     */
    if (c->parent != NULL)
        return;

    y = c->y;
    TN = y[ws->TN_i] + y[ws->Den_fl_i] / SEC_PER_DAY / c->dz_wc;
    Nfix = (ws->Nfix_i >= 0) ? y[ws->Nfix_i] : 0.0;
    TP = y[ws->TP_i];
    TC = y[ws->TC_i];

    if (ws->Amm_fl_i > -1){
      TN += y[ws->Amm_fl_i] / SEC_PER_DAY / c->dz_wc;
    }

    CO2_flux = (ws->CO2_flux_i >= 0) ? y[ws->CO2_flux_i] : 0.0;

    eps = fabs(TN - cv[ws->TN_old_i] - Nfix) / (TN + cv[ws->TN_old_i] + Nfix);

    if (eps > MASSBALANCE_EPS)
     e_warn("ecology: error: N balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);

    eps = fabs(TP - cv[ws->TP_old_i]) / (TP + cv[ws->TP_old_i]);
 
    if (eps > MASSBALANCE_EPS)
      e_warn("ecology: error: P balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);
    
    eps = fabs(TC - cv[ws->TC_old_i] + CO2_flux * e->dt / c->dz_wc) / (TC + cv[ws->TC_old_i]);
    
    if (eps > MASSBALANCE_EPS)
      e_warn("ecology: error: C balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);

    if (ws->TO_old_i > -1){
	O2_flux = (ws->O2_flux_i >= 0) ? y[ws->O2_flux_i] : 0.0;
	COD_flux = (ws->COD_flux_i >= 0) ? y[ws->COD_flux_i] : 0.0;

          TO = y[ws->Oxygen_i] - y[ws->COD_i]  - y[ws->BOD_i];

	  // add nitrate into oxygen mass balance.

	  TO = TO + y[ws->NO3_i] / 14.01 * 48.0;

	  if (ws->Amm_fl_i > -1){
	    TO += y[ws->Amm_fl_i]/2.0 * 48.0 / 14.01 / SEC_PER_DAY / c->dz_wc;
	  }

	  /* because TO can be close to zero */

	  eps = fabs(TO - cv[ws->TO_old_i] + COD_flux + O2_flux / c->dz_wc) / max(fabs(TO + COD_flux + cv[ws->TO_old_i]),8000.0);

	  if (eps > MASSBALANCE_EPS)
	    e_warn("ecology: error: Oxygen - BOD - COD (%e,%e) imbalance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", TO, cv[ws->TO_old_i], eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);
    }
}









