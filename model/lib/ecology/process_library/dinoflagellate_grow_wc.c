/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dinoflagellate_grow_wc.c
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
 *  $Id: dinoflagellate_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "dinoflagellate_grow_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double umax_t0;
    double aA;
    double psi;
    double Sh;
    double m;
    double IDF;
    int n;
    double NtoCHL;
    double DFCtoNvar;
    double Plank_resp;
    int len;
    double* grow;
    double* rate;

    /*
     * tracers
     */
    int Phy_N_i;
    int Phy_N_pr_i;
    int Phy_N_gr_i;
    int Phy_C_i;
    int NO3_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Light_i;
    int Chl_a_i;
    int Kd_i;
    int Oxy_pr_i;
    int temp_i;
    int salt_i;
    int TN_i;
    int TP_i;
    int TC_i;

    /*
     * common cell variables
     */
    int DNO3_i;
    int DPO4_i;
    int Tfactor_i;
    int umax_i;
} workspace;

void dinoflagellate_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "DFumax");
    ws->aA = try_parameter_value(e, "DFaA");
    if (isnan(ws->aA)) {
        double absorb = get_parameter_value(e, "DFabsorb");
        double rad = get_parameter_value(e, "DFrad");

        ws->aA = aa(rad, absorb);
    }
    ws->psi = try_parameter_value(e, "DFpsi");
    if (isnan(ws->psi))
        ws->psi = psi(get_parameter_value(e, "DFrad"));
    ws->Sh = get_parameter_value(e, "DFSh");
    ws->m = try_parameter_value(e, "DFm");
    if (isnan(ws->m))
        ws->m = PhyCellMass(get_parameter_value(e, "DFrad"));
    ws->IDF = get_parameter_value(e, "IDF");
    ws->n = rint(get_parameter_value(e, "DFn"));
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");
    ws->Plank_resp = get_parameter_value(e, "Plank_resp") * red_A_I;
    ws->len = extract_grow(ws->n, get_parameter_stringvalue(e, "DFtable"), &ws->grow, &ws->rate);
    ws->DFCtoNvar = get_parameter_value(e, "DFCtoNvar");

    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, "PhyD_N", e);
    ws->Phy_N_pr_i = e->find_index(tracers, "PhyD_N_pr", e);
    ws->Phy_N_gr_i = e->find_index(tracers, "PhyD_N_gr", e);
    ws->Phy_C_i = e->find_index(tracers, "PhyD_C", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->Kd_i = e->find_index(tracers, "Kd", e);
    ws->Oxy_pr_i = e->find_index(tracers, "Oxy_pr", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
    ws->salt_i = e->find_index(tracers, "salt", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    /*
     * common cell variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "DFumax", e);
}

void dinoflagellate_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = (workspace*) p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void dinoflagellate_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->grow);
    free(ws->rate);
    free(ws);
}

void dinoflagellate_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double Phy_N = y[ws->Phy_N_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    y[ws->Kd_i] += ws->aA * Phy_N * mgN2molN / ws->m / red_A_N;

    if (ws->do_mb) {
        double Phy_C = y[ws->Phy_C_i];

        y[ws->TN_i] += Phy_N;
        y[ws->TP_i] += Phy_N * red_W_P;
        y[ws->TC_i] += Phy_C;
    }
}

void dinoflagellate_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Phy_N = y[ws->Phy_N_i];
    double Phy_C = y[ws->Phy_C_i];
    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Light = y[ws->Light_i];
    double I = 2.77e18 * Light / AV;
    /*TODO KWA  confirm change*/
    /* double I = Light/E2W; */
    double diff[2];
    double conc[2];

    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];
    conc[0] = DIN * mgN2molN;
    conc[1] = DIP * mgP2molP;

    {
        double umax = cv[ws->umax_i];

        double growthrate_nut = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, ws->IDF, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
        double growthrate = growthrate_nut * (0.5 + 0.5 * tanh(10.0 * (1.0 - Phy_N * red_W_C / e_max(Phy_C))));
        double growth = Phy_N * growthrate;

        double cells = Phy_N * mgN2molN / ws->m / red_A_N;
        double photosynthesis = ws->aA * I * red_A_C * MW_Carb * cells * 1000.0 / red_A_I;
        double growth_C = photosynthesis * umax * Phy_C / e_max(photosynthesis + umax * Phy_C) * (0.5 + 0.5 * tanh(10.0 * (1.0 - Phy_C / e_max(Phy_N) / red_W_C / ws->DFCtoNvar)));

        y1[ws->Phy_N_i] += growth;
        y1[ws->NH4_i] -= growth * NH4 / DIN;
        y1[ws->NO3_i] -= growth * NO3 / DIN;
        y1[ws->DIP_i] -= growth * red_W_P;
        y1[ws->Phy_C_i] += growth_C;
        y1[ws->DIC_i] -= growth_C;
        y1[ws->Oxygen_i] += growth_C * red_W_O / red_W_C;
        y1[ws->Phy_N_pr_i] += growth * SEC_PER_DAY * c->dz_wc;
        y1[ws->Phy_N_gr_i] = growthrate / umax;
        y1[ws->Oxy_pr_i] += growth_C * red_W_O / red_W_C * SEC_PER_DAY;
    }
}

void dinoflagellate_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double Phy_N = y[ws->Phy_N_i];
    double Phy_C = y[ws->Phy_C_i];

    y[ws->Chl_a_i] += Phy_N / ws->NtoCHL;
    y[ws->TN_i] += Phy_N;
    y[ws->TP_i] += Phy_N * red_W_P;
    y[ws->TC_i] += Phy_C;
}
