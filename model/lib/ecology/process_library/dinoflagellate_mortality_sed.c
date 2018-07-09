/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dinoflagellate_mortality_sed.c
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
 *  $Id: dinoflagellate_mortality_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "dinoflagellate_mortality_sed.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double mL_t0;
    double NtoCHL;

    /*
     * tracers
     */
    int PhyD_N_i;
    int PhyD_C_i;
    int DIC_i;
    int DetPL_N_i;
    int Chl_a_i;
    int TN_i;
    int TP_i;
    int TC_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mL_i;
} workspace;

void dinoflagellate_mortality_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->mL_t0  = get_parameter_value(e, "PD_mL");
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");

    /*
     * tracers
     */
    ws->PhyD_N_i = e->find_index(tracers, "PhyD_N", e);
    ws->PhyD_C_i = e->find_index(tracers, "PhyD_C", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, "PD_mL", e);
}

void dinoflagellate_mortality_sed_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = (workspace*) p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
}

void dinoflagellate_mortality_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void dinoflagellate_mortality_sed_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;

    if (ws->do_mb) {
        double PhyD_N = y[ws->PhyD_N_i];
        double PhyD_C = y[ws->PhyD_C_i];

        y[ws->TN_i] += PhyD_N;
        y[ws->TP_i] += PhyD_N * red_W_P;
        y[ws->TC_i] += PhyD_C;
    }
}

void dinoflagellate_mortality_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double PhyD_N = y[ws->PhyD_N_i];
    double PhyD_C = y[ws->PhyD_C_i];
    double mL = cv[ws->mL_i];
    double mortality = PhyD_N * mL;
    double porosity = c->porosity;

    y1[ws->PhyD_N_i] -= mortality;
    y1[ws->DetPL_N_i] += mortality;
    y1[ws->DIC_i] += (PhyD_C - PhyD_N * red_W_C) * mL / porosity;
    y1[ws->PhyD_C_i] -= PhyD_C * mL;
}

void dinoflagellate_mortality_sed_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double PhyD_N = y[ws->PhyD_N_i];
    double PhyD_C = y[ws->PhyD_C_i];

    y[ws->Chl_a_i] += PhyD_N / ws->NtoCHL;
    y[ws->TN_i] += PhyD_N;
    y[ws->TP_i] += PhyD_N * red_W_P;
    y[ws->TC_i] += PhyD_C;
}
