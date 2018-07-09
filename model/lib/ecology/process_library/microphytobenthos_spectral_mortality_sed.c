/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/microphytobenthos_spectral_mortality_sed.c
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
 *  $Id: microphytobenthos_spectral_mortality_sed.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "microphytobenthos_spectral_mortality_sed.h"

typedef struct {
  int do_mb;                        /* flag */
    
  /* parameters */
  
  double mQ_t0;
  double KO_aer;

  /* tracers */

  int NH4_i;
  int NH4_pr_i;
  int DIP_i;
  int DIC_i;
  int MPB_N_i;
  int MPB_NR_i;
  int MPB_I_i;
  int DetPL_N_i;
  int Oxygen_i;

  int MPB_PR_i;
  int MPB_Chl_i;

  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;
  int COD_i;

  /* common cell variables */

  int Tfactor_i;
  int mQ_i;
} workspace;

void microphytobenthos_spectral_mortality_sed_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /* parameters */

  ws->mQ_t0 = get_parameter_value(e, "MPB_mQ");
  ws->KO_aer = get_parameter_value(e, "KO_aer");

  /* tracers   */

  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->DIP_i = e->find_index(tracers, "DIP", e);
  ws->DIC_i = e->find_index(tracers, "DIC", e);
  ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
  ws->MPB_N_i = e->find_index(tracers, "MPB_N", e);
  ws->MPB_NR_i = e->find_index(tracers, "MPB_NR", e);
  ws->MPB_I_i = e->find_index(tracers, "MPB_I", e);
  ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);

  ws->TN_i = e->find_index(tracers, "TN", e);
  ws->TP_i = e->find_index(tracers, "TP", e);
  ws->TC_i = e->find_index(tracers, "TC", e);
  ws->BOD_i = e->find_index(tracers, "BOD", e);
  ws->COD_i = e->try_index(tracers, "COD", e);

  ws->MPB_Chl_i = e->find_index(tracers, "MPB_Chl", e);
  ws->MPB_PR_i = e->find_index(tracers, "MPB_PR", e);

  /* common cell variables  */

  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->mQ_i = find_index_or_add(e->cv_cell, "MPB_mQ", e);
}

void microphytobenthos_spectral_mortality_sed_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = (workspace*) p->workspace;

  ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
}

void microphytobenthos_spectral_mortality_sed_destroy(eprocess* p)
{
  free(p->workspace);
}

void microphytobenthos_spectral_mortality_sed_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  
  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
  
  cv[ws->mQ_i] = ws->mQ_t0 * Tfactor;

  if (ws->do_mb) {
    double* y = c->y;
    double MPB_N = y[ws->MPB_N_i];
    double MPB_NR = y[ws->MPB_NR_i];
    double MPB_PR = y[ws->MPB_PR_i];
    double MPB_I = y[ws->MPB_I_i];

    y[ws->TN_i] += MPB_N + MPB_NR;
    y[ws->TP_i] += MPB_N * red_W_P + MPB_PR;
    y[ws->TC_i] += MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01;
    y[ws->BOD_i] += (MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01)*C_O_W;
  }
}

void microphytobenthos_spectral_mortality_sed_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double porosity = c->porosity;
  double MPB_N = y[ws->MPB_N_i];
  double MPB_NR = y[ws->MPB_NR_i];
  double MPB_I = y[ws->MPB_I_i];
  double MPB_PR = y[ws->MPB_PR_i];
  double MPB_Chl = y[ws->MPB_Chl_i];
  double Oxygen = e_max(y[ws->Oxygen_i]);

  /* MPB mortality is quadratic for MPB_N. To be consistent with other components 
     (say Chl) it is MPB_N * cv[ws->mQ_i] * MPB_Chl   */

  double mortality = MPB_N * cv[ws->mQ_i];

/* MPB_N mortality at Redfield but MPB_NR is excess nutrients so they must go into individual pools */

  y1[ws->MPB_N_i] -= MPB_N * mortality;
  y1[ws->MPB_NR_i] -= MPB_NR * mortality;
  y1[ws->MPB_I_i] -= MPB_I * mortality;
  y1[ws->MPB_PR_i] -= MPB_PR * mortality;
  y1[ws->MPB_Chl_i] -= MPB_Chl * mortality;
  y1[ws->DetPL_N_i] += MPB_N * mortality;
  y1[ws->NH4_i] += MPB_NR * mortality / porosity ;
  if (ws->NH4_pr_i > -1)
    y1[ws->NH4_pr_i] += MPB_NR * mortality * SEC_PER_DAY * c->dz_sed * porosity;
  y1[ws->DIP_i] +=  MPB_PR * mortality / porosity;
  y1[ws->DIC_i] +=  MPB_I * mortality * 106.0/1060.0*12.01 / porosity;

  double Oxygen2 = Oxygen * Oxygen; 
  double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

  y1[ws->Oxygen_i] -= MPB_I * mortality * 106.0/1060.0*12.01 * C_O_W * sigmoid / porosity;

  if (ws->COD_i > -1){
    y1[ws->COD_i] += MPB_I * mortality * 106.0/1060.0*12.01 * C_O_W * (1.0-sigmoid)/ porosity;
  }
}

void microphytobenthos_spectral_mortality_sed_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

  double MPB_N = y[ws->MPB_N_i];
  double MPB_NR = y[ws->MPB_NR_i];
  double MPB_PR = y[ws->MPB_PR_i];
  double MPB_I = y[ws->MPB_I_i];

  y[ws->TN_i] += MPB_N + MPB_NR;
  y[ws->TP_i] += MPB_N * red_W_P + MPB_PR;
  y[ws->TC_i] += MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01;
  y[ws->BOD_i] += (MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01)*C_O_W;

}
