/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nodularia_mortality_wc.c
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
 *  $Id: nodularia_mortality_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "nodularia_mortality_wc.h"

typedef struct {
    /*
     * parameters
     */
    double mQ_t0;
    double mL_t0;

    /*
     * tracers
     */
    int Phy_N_i;
    int DetPL_N_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mL_i;
    int mQ_i;
} workspace;

void nodularia_mortality_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->mQ_t0 = get_parameter_value(e, "NO_mQ");
    ws->mL_t0 = get_parameter_value(e, "NO_mL");

    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, "PhyN_N", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, "NO_mL", e);
    ws->mQ_i = find_index_or_add(e->cv_cell, "NO_mQ", e);
}

void nodularia_mortality_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void nodularia_mortality_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mQ_i] = ws->mQ_t0 * Tfactor;
    cv[ws->mL_i] = ws->mL_t0 * Tfactor;
}

void nodularia_mortality_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Phy_N = y[ws->Phy_N_i];
    double mortality = cv[ws->mQ_i] * Phy_N * Phy_N + cv[ws->mL_i] * Phy_N;

    y1[ws->DetPL_N_i] += mortality;
    y1[ws->Phy_N_i] -= mortality;
}

void nodularia_mortality_wc_postcalc(eprocess* p, void* pp)
{
}
