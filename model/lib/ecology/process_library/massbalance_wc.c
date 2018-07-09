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
 *  $Id: massbalance_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->BOD_i = e->try_index(tracers, "BOD", e);
    ws->COD_i = e->try_index(tracers, "COD", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);

    ws->Nfix_i = e->try_index(tracers, "Nfix", e);

    ws->CO2_flux_i = e->try_index(tracers, "CO2_flux", e);
    ws->O2_flux_i = e->try_index(tracers, "O2_flux", e);
   
    /*
     * common cell variables
     */
    ws->TN_old_i = find_index_or_add(e->cv_cell, "TN_old", e);
    ws->TP_old_i = find_index_or_add(e->cv_cell, "TP_old", e);
    ws->TC_old_i = find_index_or_add(e->cv_cell, "TC_old", e);

    if (ws->COD_i > -1){
      ws->TO_old_i = find_index_or_add(e->cv_cell, "TO_old", e);
    }
    /*
     * set a flag indicating doing mass balance calculations
     */
    stringtable_add_ifabscent(e->cv_model, "massbalance_wc", -1);
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

    if (ws->COD_i > -1){
      cv[ws->TO_old_i] = y[ws->Oxygen_i] - y[ws->BOD_i] - y[ws->COD_i];
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
    double TN, Nfix, TP, TC, TO, TO_old, eps, CO2_flux, O2_flux;

    /*
     * if this is a child cell, return, the parent should take care of
     * checking the mass balance
     */
    if (c->parent != NULL)
        return;

    y = c->y;
    TN = y[ws->TN_i];
    Nfix = (ws->Nfix_i >= 0) ? y[ws->Nfix_i] : 0.0;
    TP = y[ws->TP_i];
    TC = y[ws->TC_i];

    CO2_flux = (ws->CO2_flux_i >= 0) ? y[ws->CO2_flux_i] : 0.0;

    eps = fabs(TN - cv[ws->TN_old_i] - Nfix) / (TN + cv[ws->TN_old_i] + Nfix);

    if (eps > MASSBALANCE_EPS)
     e->quitfn("ecology: error: N balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);

    eps = fabs(TP - cv[ws->TP_old_i]) / (TP + cv[ws->TP_old_i]);
 
    if (eps > MASSBALANCE_EPS)
      e->quitfn("ecology: error: P balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);
    
    eps = fabs(TC - cv[ws->TC_old_i] + CO2_flux * e->dt / c->dz_wc) / (TC + cv[ws->TC_old_i]);
    
    if (eps > MASSBALANCE_EPS)
      e->quitfn("ecology: error: C balance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);

    if (ws->COD_i > -1){
	O2_flux = (ws->O2_flux_i >= 0) ? y[ws->O2_flux_i] : 0.0;

          TO = y[ws->Oxygen_i] - y[ws->COD_i]  - y[ws->BOD_i];

	  /* because TO can be close to zero */

	  eps = fabs(TO - cv[ws->TO_old_i] + O2_flux / c->dz_wc) / max(fabs(TO + cv[ws->TO_old_i]),8000.0);
	  
	  if (eps > MASSBALANCE_EPS)
	    e->quitfn("ecology: error: Oxygen - BOD - COD (%e,%e) imbalance violation in water cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", TO, cv[ws->TO_old_i], eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);
      }
}









