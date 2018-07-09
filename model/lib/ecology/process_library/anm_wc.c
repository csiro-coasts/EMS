/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/anm_wc.c
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
 *  $Id: anm_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "column.h"
#include "cell.h"
#include "einterface.h"
#include "anm_wc.h"

#define TEMP_BACK 16.0          /* background temperature */
#define TEMP_DIFF 14.0          /* difference between ANM and background
                                 * temperature */

typedef struct {
    /*
     * parameters
     */
    double rB;
    double r_UVB;
    double r_UVB_irr;
    double r_ads_20L;
    double r_ads_CL;
    double r_ads_20R;
    double Kd_RES;
    double Kd_riv;
    /*
     * tracers
     */
    int salt_i;
    int light_i;
    int ldra_i;
    int lpra_i;
    int lpra_a_i;
    int rdra_i;
    int rpra_i;
    int kd_i;
    int temp_i;
} workspace;

void anm_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->rB = get_parameter_value(e, "rb");
    ws->r_UVB = get_parameter_value(e, "r_uvb");
    ws->r_UVB_irr = get_parameter_value(e, "r_uvb_irr");
    ws->r_ads_20L = get_parameter_value(e, "r_ads_20l");
    ws->r_ads_CL = get_parameter_value(e, "r_ads_cl");
    ws->r_ads_20R = get_parameter_value(e, "r_ads_20r");
    ws->Kd_RES = get_parameter_value(e, "kd_res");
    ws->Kd_riv = get_parameter_value(e, "kd_riv");
    /*
     * tracers
     */
    ws->salt_i = e->find_index(tracers, "salt", e);
    ws->light_i = e->find_index(tracers, "light", e);
    ws->ldra_i = e->find_index(tracers, "ldra", e);
    ws->lpra_i = e->find_index(tracers, "lpra", e);
    ws->lpra_a_i = e->find_index(tracers, "lpra_a", e);
    ws->rdra_i = e->find_index(tracers, "rdra", e);
    ws->rpra_i = e->find_index(tracers, "rpra", e);
    ws->kd_i = e->find_index(tracers, "kd", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
}

void anm_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void anm_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* y = c->y;

    double temp = y[ws->temp_i];
    double anm = (temp - TEMP_BACK) / TEMP_DIFF;

    if (c->k <= 0)
        y[ws->kd_i] = ws->Kd_riv + ws->Kd_RES * anm;
    else
        y[ws->kd_i] = FLT_MAX;
}

void anm_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* conc = ia->y;
    double* y1 = ia->y1;

    double S20 = conc[ws->salt_i] / 20.0;
    double Light = conc[ws->light_i];
    double LDRA = conc[ws->ldra_i];
    double LPRA = conc[ws->lpra_i];
    double LPRA_A = conc[ws->lpra_a_i];
    double RDRA = conc[ws->rdra_i];

    double r_UVB = ws->r_UVB * Light / ws->r_UVB_irr;
    double r_ads_R = ws->r_ads_20R * S20;
    double r_ads_L = ws->r_ads_CL + ws->r_ads_20L * S20;

    y1[ws->rdra_i] += -r_ads_R * RDRA;
    y1[ws->rpra_i] += r_ads_R * RDRA;
    y1[ws->ldra_i] += -(r_ads_L + ws->rB + r_UVB) * LDRA;
    y1[ws->lpra_i] += r_ads_L * LDRA - (ws->rB + r_UVB) * LPRA;
    y1[ws->lpra_a_i] += (r_ads_L + ws->rB) * LDRA - r_UVB * LPRA_A;
}
