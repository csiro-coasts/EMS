/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/phytoplankton_mortality_sed.c
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
 *  $Id: phytoplankton_mortality_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "phytoplankton_mortality_sed.h"

typedef struct {
    int do_mb;                  /* flag */
    /*
     * flag: 1 for large phytoplankton, 0 for small phytoplankton
     */
    int large;

    /*
     * parameters
     */
    double mL_t0;
    double NtoCHL;

    /*
     * tracers
     */
    int Phy_N_i;
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

void phytoplankton_mortality_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    char* prm = p->prms->se[0]->s;

    p->workspace = ws;

    if (toupper(prm[0]) == 'L')
        ws->large = 1;
    else if (toupper(prm[0]) == 'S')
        ws->large = 0;
    else
        e_quit("error: ecology: \"%s\": \"%s(%s)\": unexpected parameter: expected \"small\" or \"large\"\n", e->processfname, p->name, prm);

    /*
     * parameters
     */
    ws->mL_t0 = get_parameter_value(e, (ws->large) ? "PhyL_mL" : "PhyS_mL");
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");

    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, (ws->large) ? "PhyL_N" : "PhyS_N", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, (ws->large) ? "PhyL_mL" : "PhyS_mL", e);
}

void phytoplankton_mortality_sed_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = (workspace*) p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
}

void phytoplankton_mortality_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void phytoplankton_mortality_sed_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;

    if (ws->do_mb) {
        double* y = c->y;
        double Phy_N = y[ws->Phy_N_i];

        y[ws->TN_i] += Phy_N;
        y[ws->TP_i] += Phy_N * red_W_P;
        y[ws->TC_i] += Phy_N * red_W_C;
    }
}

void phytoplankton_mortality_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* cv = ((cell*) ia->media)->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Phy_N = y[ws->Phy_N_i];
    double mortality = Phy_N * cv[ws->mL_i];

    y1[ws->Phy_N_i] -= mortality;
    y1[ws->DetPL_N_i] += mortality;
}

void phytoplankton_mortality_sed_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double Phy_N = y[ws->Phy_N_i];

  /*  y[ws->Chl_a_i] += Phy_N / ws->NtoCHL;*/
    y[ws->TN_i] += Phy_N;
    y[ws->TP_i] += Phy_N * red_W_P;
    y[ws->TC_i] += Phy_N * red_W_C;
}
