/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/anm_epi.c
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
 *  $Id: anm_epi.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "anm_epi.h"

typedef struct {
    /*
     * parameters
     */
    double D_RES;
    double F_ads;
    double nu;
    /*
     * tracers
     */
    int wc_ads_ldra_i;
    int wc_ads_rdra_i;
    int wc_ldra_i;
    int wc_rdra_i;
    int sed_ldra_i;
    int sed_rdra_i;
} workspace;

void anm_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    int OFFSET_SED = tracers->n;

    p->workspace = ws;

    /*
     * parameters
     */
    ws->D_RES = get_parameter_value(e, "D_RES");
    ws->F_ads = get_parameter_value(e, "F_ads");
    ws->nu = get_parameter_value(e, "nu");
    /*
     * tracers
     */
    ws->wc_ads_ldra_i = e->find_index(tracers, "ads_LDRA", e);
    ws->wc_ads_rdra_i = e->find_index(tracers, "ads_RDRA", e);
    ws->wc_ldra_i = e->find_index(tracers, "LDRA", e);
    ws->wc_rdra_i = e->find_index(tracers, "RDRA", e);
    ws->sed_ldra_i = e->find_index(tracers, "LDRA", e) + OFFSET_SED;
    ws->sed_rdra_i = e->find_index(tracers, "RDRA", e) + OFFSET_SED;
}

void anm_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void anm_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* conc = ia->y;
    double* y1 = ia->y1;

    double dz_wc = c->dz_wc;
    double dz_sed = c->dz_sed;
    double porosity = c->porosity;
    double ustrcw = einterface_getustrcw(p->ecology->model, c->b);
    double ads_koeff = ws->D_RES * ws->F_ads * ustrcw / ws->nu;

    double wc_LDRA = conc[ws->wc_ldra_i];
    double wc_RDRA = conc[ws->wc_rdra_i];
    double ads_LDRA = ads_koeff * wc_LDRA;
    double ads_RDRA = ads_koeff * wc_RDRA;

    y1[ws->wc_ldra_i] -= ads_LDRA / dz_wc;
    y1[ws->wc_rdra_i] -= ads_RDRA / dz_wc;
    y1[ws->sed_ldra_i] += ads_LDRA / (dz_sed * porosity);
    y1[ws->sed_rdra_i] += ads_RDRA / (dz_sed * porosity);

    /*
     * diagnostics
     */
    y1[ws->wc_ads_ldra_i] += ads_LDRA;
    y1[ws->wc_ads_rdra_i] += ads_RDRA;
}
