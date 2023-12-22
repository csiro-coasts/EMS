/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/trichodesmium_mortality_photoacclimate.c
 *  
 *  Description:
 *  Quadratic mortality of Trichodesmium due to phages
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: trichodesmium_mortality_photoacclimate.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "trichodesmium_mortality_photoacclimate.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double Tricho_crit;
    double mQ_t0;
    double KO_aer;

    /*
     * tracers
     */
    int Tricho_N_i;
    int Tricho_NR_i;
    int Tricho_I_i;
    int Tricho_PR_i;
    int Tricho_Chl_i;
    int Oxygen_i;
    int COD_i;
    int Oxy_pr_i;
    int DIC_i;
    int DIP_i;
    int NH4_i;
    int NH4_pr_i;
    int DetPL_N_i;
    int TN_i;
    int TP_i;
    int TC_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mQ_i;
} workspace;

void trichodesmium_mortality_photoacclimate_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->mQ_t0 = get_parameter_value(e, "Tricho_mQ");
    ws->Tricho_crit = get_parameter_value(e, "Tricho_crit");
    ws->KO_aer = get_parameter_value(e, "KO_aer");

    /*
     * tracers
     */

    ws->Tricho_N_i = e->find_index(tracers, "Tricho_N", e);
    ws->Tricho_NR_i = e->find_index(tracers, "Tricho_NR", e);
    ws->Tricho_I_i = e->find_index(tracers, "Tricho_I", e);

    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->COD_i = e->try_index(tracers, "COD", e);
    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);

    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);

    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    ws->Tricho_PR_i = e->find_index(tracers, "Tricho_PR", e);
    ws->Tricho_Chl_i = e->find_index(tracers, "Tricho_Chl", e);

    /*
     * common cell variables
     */

    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mQ_i = find_index_or_add(e->cv_cell, "Tricho_mQ", e);
}

void trichodesmium_mortality_photoacclimate_postinit(eprocess* p)
{
}

void trichodesmium_mortality_photoacclimate_destroy(eprocess* p)
{
    free(p->workspace);
}

void trichodesmium_mortality_photoacclimate_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mQ_i] = ws->mQ_t0 * Tfactor;
}

void trichodesmium_mortality_photoacclimate_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Tricho_N = y[ws->Tricho_N_i];
    double Tricho_NR = y[ws->Tricho_NR_i];
    double Tricho_PR = y[ws->Tricho_PR_i];
    double Tricho_I = y[ws->Tricho_I_i];
    double Tricho_Chl = y[ws->Tricho_Chl_i];
    double Oxygen = e_max(y[ws->Oxygen_i]);
    double mortality = (Tricho_N > ws->Tricho_crit) ? (Tricho_N - ws->Tricho_crit) * (Tricho_N - ws->Tricho_crit) * cv[ws->mQ_i] : 0.0;
    double NH4release = Tricho_NR / Tricho_N * mortality;
    double energy_loss = Tricho_I /Tricho_N * mortality;
    double DIPrelease = Tricho_PR / Tricho_N * mortality;
    double Chla_loss = Tricho_Chl / Tricho_N * mortality;

    double Oxygen2 = Oxygen * Oxygen; 
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    double Oxy_pr = -energy_loss * 106.0/1060.0*12.01 * C_O_W * sigmoid;

    if (ws->COD_i > -1){
      y1[ws->COD_i] += energy_loss * 106.0/1060.0*12.01 * C_O_W *(1.0-sigmoid);
    }
    
    /* Structural material (Tricho_N) released at Redfield to DetPL, but all 
       reserves must go to individual pools (NR, PR), or lost (I, Chl)  */
    
    y1[ws->Tricho_N_i] -= mortality;
    y1[ws->Tricho_NR_i] -= NH4release;
    y1[ws->Tricho_I_i] -= energy_loss;
    y1[ws->Tricho_PR_i] -= DIPrelease;
    y1[ws->Tricho_Chl_i] -= Chla_loss;

    y1[ws->DetPL_N_i] += mortality ;
    y1[ws->NH4_i] += NH4release ;
    if (ws->NH4_pr_i > -1)
      y1[ws->NH4_pr_i] += NH4release * SEC_PER_DAY ;
    y1[ws->DIP_i] += DIPrelease ;
    y1[ws->DIC_i] += energy_loss * 106.0/1060.0*12.01;
    y1[ws->Oxygen_i] += Oxy_pr;
    if (ws->Oxy_pr_i > -1)
      y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY;
}

void trichodesmium_mortality_photoacclimate_postcalc(eprocess* p, void* pp)
{
} 

