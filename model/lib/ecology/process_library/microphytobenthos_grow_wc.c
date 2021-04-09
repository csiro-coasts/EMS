/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/microphytobenthos_grow_wc.c
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
 *  $Id: microphytobenthos_grow_wc.c 6697 2021-03-24 01:10:38Z wil00y $
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
#include "microphytobenthos_grow_wc.h"

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
    double Plank_resp;
    int len;
    double* grow;
    double* rate;

    /*
     * tracers
     */
    int MPB_N_i;
    int NO3_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Light_i;
    int MPB_N_pr_i;
    int MPB_N_gr_i;
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

void microphytobenthos_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "MBumax");
    ws->aA = try_parameter_value(e, "MBaA");
    if (isnan(ws->aA)) {
        /*UR CHANGED to exit with some meaningful error message */
        double absorb = try_parameter_value(e, "MBabsorb");
        double rad = try_parameter_value(e, "MBrad");
        if(isnan(rad) || isnan(absorb))
          e->quitfn("If tracer 'MBaA' is not provided 'MBrad' and 'MBabsorb' are required to calculate microphytobenthos_grow_wc!");
        ws->aA = aa(rad, absorb);
    }
    ws->psi = try_parameter_value(e, "MBpsi");
    if (isnan(ws->psi))
    {
      if(isnan(try_parameter_value(e, "MBrad")))
        e->quitfn("Either tracer 'MBpsi' or 'MBrad' are required to calculate microphytobenthos_grow_wc!");
      ws->psi = psi(try_parameter_value(e, "MBrad"));
    }
    ws->Sh = get_parameter_value(e, "MBSh");
    ws->m = try_parameter_value(e, "MBm");
    if (isnan(ws->m))
    {
        /*UR CHANGED to exit with some meaningful error message */
        if(isnan(try_parameter_value(e, "MBrad")))
          e->quitfn("Either tracer 'MBm' or 'MBrad' are required to calculate microphytobenthos_grow_wc!");
        ws->m = PhyCellMass(try_parameter_value(e, "MBrad"));
    }
    ws->n = rint(get_parameter_value(e, "MBn"));
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");
    ws->Plank_resp = get_parameter_value(e, "Plank_resp") * red_A_I;
    {
        char* table = get_parameter_stringvalue(e, "MBtable");

        ws->len = extract_grow(ws->n, table, &ws->grow, &ws->rate);
    }

    /*
     * tracers
     */
    ws->MPB_N_i = e->find_index(tracers, "MPB_N", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->MPB_N_pr_i = e->find_index(tracers, "MPB_N_pr", e);
    ws->MPB_N_gr_i = e->find_index(tracers, "MPB_N_gr", e);
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
    ws->umax_i = find_index_or_add(e->cv_cell, "MBumax", e);
}

void microphytobenthos_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void microphytobenthos_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->grow);
    free(ws->rate);
    free(ws);
}

void microphytobenthos_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MPB_N = y[ws->MPB_N_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    y[ws->Kd_i] += ws->aA * MPB_N * mgN2molN / ws->m / red_A_N;

    if (ws->do_mb) {
        y[ws->TN_i] += MPB_N;
        y[ws->TP_i] += MPB_N * red_W_P;
        y[ws->TC_i] += MPB_N * red_W_C;
    }
}

void microphytobenthos_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double MPB_N = y[ws->MPB_N_i];
    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Light = y[ws->Light_i];
    double I = 2.77e18 * Light / AV;
    double diff[2];
    double conc[2];

    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];
    conc[0] = DIN * mgN2molN;
    conc[1] = DIP * mgP2molP;

    {
        double umax = cv[ws->umax_i];

        double growthrate = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, I, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
        double growth = MPB_N * growthrate;

        y1[ws->MPB_N_i] += growth;
        y1[ws->NH4_i] -= growth * NH4 / DIN;
        y1[ws->NO3_i] -= growth * NO3 / DIN;
        y1[ws->DIP_i] -= growth * red_W_P;
        y1[ws->DIC_i] -= growth * red_W_C;
        y1[ws->Oxygen_i] += growth * red_W_O;
        y1[ws->MPB_N_pr_i] += growth * SEC_PER_DAY * red_W_C;
        y1[ws->MPB_N_gr_i] = growthrate * SEC_PER_DAY;
        y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
    }
}

void microphytobenthos_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double MPB_N = y[ws->MPB_N_i];

    y[ws->Chl_a_i] += MPB_N / ws->NtoCHL;
    y[ws->TN_i] += MPB_N;
    y[ws->TP_i] += MPB_N * red_W_P;
    y[ws->TC_i] += MPB_N * red_W_C;
}
