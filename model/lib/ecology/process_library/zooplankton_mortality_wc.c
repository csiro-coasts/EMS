/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/zooplankton_mortality_wc.c
 *  
 *  Description: Zooplankton mortality based on a quadratic mortality rate coefficient.
 *
 *  Options: zooplankton_spectral_mortality_wc(small|large)
 *
 *  Small - State variables ZooS_*, parameters ZS*
 *  Large - State variables ZooL_*, parameters ZL*
 *
 *  WARNING: Mass balance check requires that if this processes is active, then so must be  
 *           phytoplankton_spectral_grow_wc
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: zooplankton_mortality_wc.c 5934 2018-09-12 03:34:13Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "zooplankton_mortality_wc.h"

typedef struct {
    /*
     * flag: 1 for large zooplankton, 0 for small zooplankton
     */
    int large;

    /*
     * parameters
     */
    double mQ_t0;
    double FDM;
    double KO_aer;

    /*
     * tracers
     */
    int Zoo_N_i;
    int Oxygen_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int DetPL_N_i;
    int Oxy_pr_i;
    int NH4_pr_i;
    int COD_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mQ_i;
} workspace;

void zooplankton_mortality_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    char* prm = p->prms->se[0]->s;
    int large;

    p->workspace = ws;

    if (toupper(prm[0]) == 'L')
        ws->large = 1;
    else if (toupper(prm[0]) == 'S')
        ws->large = 0;
    else
        e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected parameter: expected \"small\" or \"large\"\n", e->processfname, p->name, prm);
    large = ws->large;

    /*
     * parameters
     */
    ws->mQ_t0 = get_parameter_value(e, (large) ? "ZL_mQ" : "ZS_mQ");
    ws->FDM = get_parameter_value(e, (large) ? "ZL_FDM" : "ZS_FDM");
    ws->KO_aer = get_parameter_value(e, "KO_aer");

    /*
     * tracers
     */
    ws->Zoo_N_i = e->find_index(tracers, (large) ? "ZooL_N" : "ZooS_N", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->COD_i = e->find_index(tracers, "COD", e);

  /*non essential diagnostic variables*/

    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mQ_i = find_index_or_add(e->cv_cell, (large) ? "ZL_mQ" : "ZS_mQ", e);
}

void zooplankton_mortality_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void zooplankton_mortality_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mQ_i] = ws->mQ_t0 * Tfactor;
}

void zooplankton_mortality_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    double* cv = ((cell*) ia->media)->cv;

    double Zoo_N = y[ws->Zoo_N_i];
    double Oxygen = e_max(y[ws->Oxygen_i]);
    double mortality = Zoo_N * Zoo_N * cv[ws->mQ_i];
    double NH4release = mortality * (1.0 - ws->FDM);

    double Oxygen2 = Oxygen * Oxygen; 
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    double Oxy_pr = -NH4release * red_W_O * sigmoid;

    y1[ws->Zoo_N_i] -= mortality;
    y1[ws->NH4_i] += NH4release;
    y1[ws->DIP_i] += NH4release * red_W_P;
    y1[ws->DIC_i] += NH4release * red_W_C;
    y1[ws->DetPL_N_i] += mortality * ws->FDM;
    y1[ws->Oxygen_i] += Oxy_pr;

    if (ws->COD_i > -1)
      y1[ws->COD_i] += NH4release * red_W_O * (1.0-sigmoid);

 if (ws-> NH4_pr_i > -1)
    y1[ws->NH4_pr_i] += NH4release * SEC_PER_DAY;

 if (ws->  Oxy_pr_i> -1)
    y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY;
}

void zooplankton_mortality_wc_postcalc(eprocess* p, void* pp)
{
}
