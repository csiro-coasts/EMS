/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_grow_epi.c
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
 *  $Id: macroalgae_grow_epi.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "constants.h"
#include "macroalgae_grow_epi.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double umax_t0;
    double aA;
    double m;
    int n;
    double* grow;
    double* rate;
    double len;
    double Benth_resp;
    double cf;
    double Ub;
    double ks;

    /*
     * epis
     */
    int MA_N_i;
    int MA_N_pr_i;
    int MA_N_gr_i;
    int Epilight_i;
    int Epilightatt_i;
    int EpiTN_i;
    int EpiTP_i;
    int EpiTC_i;

    /*
     * tracers
     */
    int DIC_wc_i;
    int NH4_wc_i;
    int NO3_wc_i;
    int DIP_wc_i;
    int Oxygen_wc_i;
    int Oxy_pr_wc_i;

    /*
     * common cell variables
     */
    int DNO3_i;
    int DPO4_i;
    int Tfactor_i;
    int MAumax_i;
} workspace;

void macroalgae_grow_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "MAumax");
    ws->aA = get_parameter_value(e, "MAaA");
    ws->m = get_parameter_value(e, "MAm");
    ws->n = rint(get_parameter_value(e, "MAn"));
    {
        char* table = get_parameter_stringvalue(e, "MAtable");

        ws->len = extract_grow(ws->n, table, &ws->grow, &ws->rate);
    }
    ws->Benth_resp = get_parameter_value(e, "Benth_resp") * atk_A_I;
    ws->cf = get_parameter_value(e, "cf");
    ws->Ub = get_parameter_value(e, "Ub");
    ws->ks = get_parameter_value(e, "ks");

    /*
     * epis
     */
    ws->MA_N_i = e->find_index(epis, "MA_N", e) + OFFSET_EPI;
    ws->MA_N_pr_i = e->find_index(epis, "MA_N_pr", e) + OFFSET_EPI;
    ws->MA_N_gr_i = e->find_index(epis, "MA_N_gr", e) + OFFSET_EPI;
    ws->Epilight_i = e->find_index(epis, "Epilight", e) + OFFSET_EPI;
    ws->Epilightatt_i = e->find_index(epis, "Epilightatt", e) + OFFSET_EPI;
    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;

    /*
     * tracers
     */
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->Oxy_pr_wc_i = e->find_index(tracers, "Oxy_pr", e);
    ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);
    ws->DIP_wc_i = e->find_index(tracers, "DIP", e);

    /*
     * common variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->MAumax_i = find_index_or_add(e->cv_cell, "MAumax", e);
}

void macroalgae_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;
}

void macroalgae_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->grow);
    free(ws->rate);
    free(ws);
}

void macroalgae_grow_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MA_N = y[ws->MA_N_i];

    cv[ws->MAumax_i] = ws->umax_t0 * Tfactor;

    y[ws->Epilightatt_i] += MA_N * ws->aA;

    if (ws->do_mb) {
        y[ws->EpiTN_i] += MA_N;
        y[ws->EpiTP_i] += MA_N * atk_W_P;
        y[ws->EpiTC_i] += MA_N * atk_W_C;
    }
}

void macroalgae_grow_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double dz_wc = c->dz_wc;

    double MA_N = y[ws->MA_N_i];
    double NO3_wc = y[ws->NO3_wc_i];
    double NH4_wc = y[ws->NH4_wc_i];
    double din_wc = NO3_wc + NH4_wc;
    double DIN_wc = (din_wc > 0.0) ? din_wc : EPS_DIN;
    double dip_wc = y[ws->DIP_wc_i];
    double DIP_wc = (dip_wc > 0.0) ? dip_wc : EPS_DIP;
    double Epilight = y[ws->Epilight_i];
    double I = 2.77e18 * Epilight / AV;
    double kI = I * (1.0 - exp(-MA_N * ws->aA));
    double conc[2];
    double diff[2];

    conc[0] = DIN_wc * mgN2molN;
    conc[1] = DIP_wc * mgP2molP;
    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];

    {
        double umax = cv[ws->MAumax_i];
        double growthrate = benthgrow(umax, ws->cf, KINEMATIC_VISCOSITY, ws->Ub, ws->ks, diff, ws->m * MA_N, kI, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Benth_resp);

        double growth = MA_N * growthrate;
        double Oxy_pr = growth * atk_W_O / dz_wc;

        y1[ws->MA_N_i] += growth;
        y1[ws->MA_N_gr_i] += growthrate / umax;
        y1[ws->MA_N_pr_i] += growth * SEC_PER_DAY;
        y1[ws->NH4_wc_i] -= growth * NH4_wc / DIN_wc / dz_wc;
        y1[ws->NO3_wc_i] -= growth * NO3_wc / DIN_wc / dz_wc;
        y1[ws->DIP_wc_i] -= growth * atk_W_P / dz_wc;
        y1[ws->DIC_wc_i] -= growth * atk_W_C / dz_wc;
        y1[ws->Oxygen_wc_i] += Oxy_pr;
        y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY;
    }
}

void macroalgae_grow_epi_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double MA_N = y[ws->MA_N_i];

    y[ws->EpiTN_i] += MA_N;
    y[ws->EpiTP_i] += MA_N * atk_W_P;
    y[ws->EpiTC_i] += MA_N * atk_W_C;
}
