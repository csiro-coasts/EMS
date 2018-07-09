/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nitrification_wc.c
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
 *  $Id: nitrification_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <emslogger.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "nitrification_wc.h"


typedef struct {
  /*
   * parameters
   */
  double r_nit_wc_t0;
  double KO_Nit;
  
  /*
   * tracers
   */
  int NH4_i;
  int NH4_pr_i;
  int NO3_i;
  int Oxygen_i;
  int Oxy_pr_i;
  int COD_i;
  
  /*
   * common variables
   */
  int Tfactor_i;
  int r_nit_wc_i;
} workspace;

void nitrification_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*
   * parameters
   */
  ws->r_nit_wc_t0 = get_parameter_value(e, "r_nit_wc");
  ws->KO_Nit = get_parameter_value(e, "KO_Nit");
  /*
   * tracers
   */
  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->NO3_i = e->find_index(tracers, "NO3", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
  ws->COD_i = e->try_index(tracers, "COD", e);
  
  /*non essential diagnostic tracer*/
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
  ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
  
  /*
   * common variables
   */
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->r_nit_wc_i = find_index_or_add(e->cv_cell, "r_nit_wc", e);
}

void nitrification_wc_destroy(eprocess* p)
{
  free(p->workspace);
}

void nitrification_wc_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  
  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
  
  cv[ws->r_nit_wc_i] = ws->r_nit_wc_t0 * Tfactor;
}

void nitrification_wc_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double NH4 = y[ws->NH4_i];
  double Oxygen = e_max(y[ws->Oxygen_i]);
  double Nitrification = cv[ws->r_nit_wc_i] * NH4 * Oxygen / (ws->KO_Nit + Oxygen);
  
  y1[ws->NH4_i] -= Nitrification;
  y1[ws->NO3_i] += Nitrification;
  
  /*
   * KWA added : 2 moles of DO lost for every mole of N == 4.57 g (NIT_N_0)
   *             DO per g of N nitrified and not denitrified 
   */
  y1[ws->Oxygen_i] -= Nitrification * NIT_N_0 ;

  /* Need to account for oxygen */

  if (ws->COD_i > -1)
    y1[ws->COD_i] -= Nitrification * NIT_N_0 ;
  
  if (ws->NH4_pr_i> -1)
    y1[ws->NH4_pr_i] -= Nitrification * SEC_PER_DAY;
  
  /* porosity and layer thickness were erroneously in next line until March 2014 */
  
  if (ws->Oxy_pr_i> -1)
    y1[ws->Oxy_pr_i] -= Nitrification * NIT_N_0 * SEC_PER_DAY;  
}

void nitrification_wc_postcalc(eprocess* p, void* pp)
{
}
