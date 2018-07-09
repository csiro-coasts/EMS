/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/microphytobenthos_diel_mortality_sed.c
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
 *  $Id: microphytobenthos_diel_mortality_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "microphytobenthos_diel_mortality_sed.h"

typedef struct {
    /*
     * parameters
     */
    double mQ_t0;

    /*
     * tracers
     */
    int NH4_i;
    int NH4_pr_i;
    int DIP_i;
    int MPB_N_i;
    int MPB_NR_i;
    int MPB_I_i;
    int DetPL_N_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mQ_i;
} workspace;

void microphytobenthos_diel_mortality_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->mQ_t0 = get_parameter_value(e, "MPB_mQ");

    /*
     * tracers
     */
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->NH4_pr_i = e->find_index(tracers, "NH4_pr", e);
    ws->MPB_N_i = e->find_index(tracers, "MPB_N", e);
    ws->MPB_NR_i = e->find_index(tracers, "MPB_NR", e);
    ws->MPB_I_i = e->find_index(tracers, "MPB_I", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mQ_i = find_index_or_add(e->cv_cell, "MPB_mQ", e);
}

void microphytobenthos_diel_mortality_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void microphytobenthos_diel_mortality_sed_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mQ_i] = ws->mQ_t0 * Tfactor;
}

void microphytobenthos_diel_mortality_sed_calc(eprocess* p, void* pp)
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
    double mortality1 = MPB_N  * (MPB_N * cv[ws->mQ_i]);
    double mortality2 = MPB_NR * (MPB_N * cv[ws->mQ_i]);
    double mortality3 = MPB_I  * (MPB_N * cv[ws->mQ_i]);

/* MPB_N mortality at Redfield but MPB_NR is excess nutrients so they must go into individual pools */

    y1[ws->MPB_N_i] -= mortality1;
    y1[ws->MPB_NR_i] -= mortality2;
    y1[ws->MPB_I_i] -= mortality3;
    y1[ws->DetPL_N_i] += mortality1;
    y1[ws->NH4_i] += mortality2 / porosity ;
    y1[ws->NH4_pr_i] += mortality2 * SEC_PER_DAY * c->dz_sed * porosity;
    y1[ws->DIP_i] += mortality2 * red_W_P / porosity;
}

void microphytobenthos_diel_mortality_sed_postcalc(eprocess* p, void* pp)
{
}
