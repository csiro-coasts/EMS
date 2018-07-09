/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nodularia_grow_wc.c
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
 *  $Id: nodularia_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "nodularia_grow_wc.h"

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
    int n;
    double NtoCHL;
    double NO_TG;               /* former TempgrowBG */
    double NO_SG;               /* former SaltgrowBG */
    double Plank_resp;
    int len;
    double* grow;
    double* rate;

    /*
     * tracers
     */
    int Phy_N_i;
    int NO3_i;
    int Nfix_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Light_i;
    int Phy_N_pr_i;
    int Phy_N_gr_i;
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

void nodularia_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "NOumax");
    ws->aA = try_parameter_value(e, "NOaA");
    if (isnan(ws->aA)) {
        double absorb = get_parameter_value(e, "NOabsorb");
        double rad = get_parameter_value(e, "NOrad");

        ws->aA = aa(rad, absorb);
    }
    ws->psi = try_parameter_value(e, "NOpsi");
    if (isnan(ws->psi))
        ws->psi = psi(get_parameter_value(e, "NOrad"));
    ws->Sh = get_parameter_value(e, "NOSh");
    ws->m = try_parameter_value(e, "NOm");
    if (isnan(ws->m))
        ws->m = PhyCellMass(get_parameter_value(e, "NOrad"));
    ws->n = rint(get_parameter_value(e, "NOn"));
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");
    ws->NO_TG = get_parameter_value(e, "NO_TG");
    ws->NO_SG = get_parameter_value(e, "NO_SG");
    ws->Plank_resp = get_parameter_value(e, "Plank_resp") * red_A_I;
    ws->len = extract_grow(ws->n, get_parameter_stringvalue(e, "NOtable"), &ws->grow, &ws->rate);

    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, "PhyN_N", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->Nfix_i = e->find_index(tracers, "Nfix", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Phy_N_pr_i = e->find_index(tracers, "PhyN_N_pr", e);
    ws->Phy_N_gr_i = e->find_index(tracers, "PhyN_N_gr", e);
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
    ws->umax_i = find_index_or_add(e->cv_cell, "NOumax", e);
}

void nodularia_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void nodularia_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->grow);
    free(ws->rate);
    free(ws);
}

void nodularia_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double temp = y[ws->temp_i];
    double salt = y[ws->salt_i];
    double Phy_N = y[ws->Phy_N_i];
    double tfact = exp(temp - ws->NO_TG);
    double sfact = exp(ws->NO_SG - salt);
    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;
    cv[ws->umax_i] *= tfact / (1.0 + tfact) * sfact / (1.0 + sfact);

    y[ws->Kd_i] += ws->aA * Phy_N * mgN2molN / ws->m / red_A_N;

    if (ws->do_mb) {
        y[ws->TN_i] += Phy_N;
        y[ws->TP_i] += Phy_N * red_W_P;
        y[ws->TC_i] += Phy_N * red_W_C;
    }
}

void nodularia_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Phy_N = y[ws->Phy_N_i];
    double NH4 = y[ws->NH4_i];
    /*double NO3 = y[ws->NO3_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN; */
    double DIN2 = (NH4 > 0.0) ? NH4 : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Light = y[ws->Light_i];
    double I = 2.77e18 * Light / AV;
    double diff[2];
    double conc[2];
    double concfix[2];

    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];
    /* conc[0] = DIN * mgN2molN;
     * UR correction according to J.Parslow 12/2004 - Nodularia can only take up NH4
     */
    conc[0] = DIN2 * mgN2molN;
    conc[1] = DIP * mgP2molP;
    concfix[0] = 100.0;
    concfix[1] = DIP * mgP2molP;

    {
        double umax = cv[ws->umax_i];

        double growthrate = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, I, concfix, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
        double growthrate_N = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, I, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
        /* UR is this appropriate
        if(growthrate < growthrate_N)
          growthrate = growthrate_N;
        */
        double growth = Phy_N * growthrate;
        double growth_N = Phy_N * growthrate_N;
        double uptake_N = e_min(growth, growth_N);

        y1[ws->Phy_N_i] += growth;
        y1[ws->NH4_i] -= uptake_N;
        y1[ws->Nfix_i] += growth - uptake_N;
        y1[ws->DIP_i] -= growth * red_W_P;
        y1[ws->DIC_i] -= growth * red_W_C;
        y1[ws->Oxygen_i] += growth * red_W_O;
        y1[ws->Phy_N_pr_i] += growth * SEC_PER_DAY * c->dz_wc;
        y1[ws->Phy_N_gr_i] = growthrate / umax;
        y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
    }
}

void nodularia_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double Phy_N = y[ws->Phy_N_i];

    y[ws->Chl_a_i] += Phy_N / ws->NtoCHL;
    y[ws->TN_i] += Phy_N;
    y[ws->TP_i] += Phy_N * red_W_P;
    y[ws->TC_i] += Phy_N * red_W_C;
}
