/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_mortality_epi.c
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
 *  $Id: macroalgae_mortality_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "macroalgae_mortality_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
    /*
     * parameters
     */
  double mL_t0;
  double Lamda;
  double K;

    /*
     * epis
     */
    

    /*
     * tracers
     */
  int DetBL_N_wc_i;
  int SWF_N_i;
  int SWS_N_i;
  int NH4_i;
  
    /*
     * common cell variables
     */
    int Tfactor_i;
    int mL_i;
  
  double time;
  /*
   * 3D velocity
   */
  // double vel[3]; 
} workspace;

void macroalgae_mortality_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;




    p->workspace = ws;

    /*
     * parameters
     */
    ws->mL_t0 = get_parameter_value(e, "MA_mL");
    ws->Lamda = get_parameter_value(e, "MA_Lamda");
    ws->K = get_parameter_value(e, "MA_K");

    /*
     * epis
     */
   

    /*
     * tracers
     */
    ws->DetBL_N_wc_i = e->find_index(tracers, "DetBL_N", e);
    ws->SWF_N_i = e->find_index(tracers, "SWF_N", e);
    ws->SWS_N_i = e->find_index(tracers, "SWS_N", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    
    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, "MA_mL", e);
    ws->time = 0;
}

void macroalgae_mortality_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void macroalgae_mortality_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;
    ws->time = (c->e->nstep)*(c->e->dt)/ SEC_PER_DAY;
    // Get 3d velocity from the interface
    // einterface_3d_vel(c->col->model, c->b,  c->k_wc, ws->vel);
   
}

void macroalgae_mortality_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;


    double SWF_N = y[ws->SWF_N_i];
    double SWS_N = y[ws->SWS_N_i];

    double Lamda = ws->Lamda;
    double K = ws->K;
    double time = (ws->time);
    double age  = 105.0*(1-exp(-pow(time/Lamda,K)));
    double mortality = (K/Lamda)*pow(age/Lamda,K-1)*exp(-pow(age/Lamda,K))/SEC_PER_DAY;
    // double mortality = c->cv[ws->mL_i];
    
    y1[ws->SWF_N_i] -= mortality * SWF_N;
    y1[ws->SWS_N_i] -= mortality * SWS_N;
    y1[ws->DetBL_N_wc_i] += mortality * SWF_N * 1000;
    y1[ws->NH4_i] += mortality * SWS_N * 1000;
    
}

void macroalgae_mortality_wc_postcalc(eprocess* p, void* pp)
{





 
}
