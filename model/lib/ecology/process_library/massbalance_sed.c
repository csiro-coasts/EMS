/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/massbalance_sed.c
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
 *  $Id: massbalance_sed.c 6689 2021-03-24 00:54:38Z wil00y $
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
#include "massbalance_sed.h"

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
  int Den_fl_i;
  int Amm_fl_i;
  int NO3_i;
  
  /*
   * common cell variables
   */
  int TN_old_i;
  int TP_old_i;
  int TC_old_i;
  int TO_old_i;

} workspace;

void massbalance_sed_init(eprocess* p)
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
    ws->Den_fl_i = e->find_index(tracers, "Den_fl", e);
    ws->Amm_fl_i = e->try_index(tracers, "Amm_fl", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);

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

    stringtable_add_ifabscent(e->cv_model, "massbalance_sed", -1);

    eco_write_setup(e,"\n Mass balance in sediment to fractional difference of %e \n",MASSBALANCE_EPS*1000.0);
}

void massbalance_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void massbalance_sed_precalc(eprocess* p, void* pp)
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

    double porosity = c->porosity;

    cv[ws->TN_old_i] = y[ws->TN_i];
    cv[ws->TP_old_i] = y[ws->TP_i];
    cv[ws->TC_old_i] = y[ws->TC_i];

    if (ws->COD_i > -1){
      cv[ws->TO_old_i] = (y[ws->Oxygen_i] - y[ws->COD_i] + y[ws->NO3_i] / 14.01 * 48.0 ) * porosity - y[ws->BOD_i];

      y[ws->BOD_i] = 0.0;
    }

    y[ws->TN_i] = 0.0;
    y[ws->TP_i] = 0.0;
    y[ws->TC_i] = 0.0;

}

void massbalance_sed_postcalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y;
    double TN, TP, TC, TO, TO_old, eps;

    /*
     * if this is a child cell, return, the parent should take care of
     * checking the mass balance
     */
    if (c->parent != NULL)
        return;

    y = c->y;

    double porosity = c->porosity;

    TN = y[ws->TN_i] + y[ws->Den_fl_i] / SEC_PER_DAY / c->dz_sed;
    TP = y[ws->TP_i];
    TC = y[ws->TC_i];

    if (ws->Amm_fl_i > -1){
      TN += y[ws->Amm_fl_i] / SEC_PER_DAY / c->dz_sed ;
    }

    eps = fabs(TN - cv[ws->TN_old_i]) / (TN + cv[ws->TN_old_i]);
    if (eps > MASSBALANCE_EPS*1000.0)
        e_warn("ecology: error: N balance violation in sediment cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_sed);
    eps = fabs(TP - cv[ws->TP_old_i]) / (TP + cv[ws->TP_old_i]);

    if (eps > MASSBALANCE_EPS*1000.0)
        e_warn("ecology: error: P balance violation in sediment cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_sed);
    eps = fabs(TC - cv[ws->TC_old_i]) / (TC + cv[ws->TC_old_i]);

    if (eps > MASSBALANCE_EPS*1000.0)
        e_warn("ecology: error: C balance violation in sediment cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", eps, e->nstep, c->nsubstep, c->col->b, c->k_sed);

    if (ws->COD_i > -1){

      TO = (y[ws->Oxygen_i] - y[ws->COD_i]) * porosity - y[ws->BOD_i];

      // add dissolved nitrate into oxygen mass balance.

      TO = TO + y[ws->NO3_i] / 14.01 * 48.0 * porosity ;

      if (ws->Amm_fl_i > -1){
	TO += y[ws->Amm_fl_i]/2.0 * 48.0 / 14.01 / SEC_PER_DAY / c->dz_sed ;
      }

      /* because TO can be close to zero */

      eps = fabs(TO - cv[ws->TO_old_i]) / max(fabs(TO + cv[ws->TO_old_i]),8000.0);
     
      /* The balance for oxygen is Oxygen - Biological Ocean Demand  */ 
      
      if (eps > MASSBALANCE_EPS*1000.0)
       	e_warn("ecology: error: Oxygen - BOD - COD imbalance (%e,%e) violation in sediment cell by %.3g, nstep = %d, nsubstep = %d, b = %d, k = %d\n", TO, cv[ws->TO_old_i], eps, e->nstep, c->nsubstep, c->col->b, c->k_wc);
    }
}
