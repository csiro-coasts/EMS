/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_grow_epi.c
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
 *  $Id: seagrass_grow_epi.c 6699 2021-03-24 01:12:35Z wil00y $
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
#include "seagrass_grow_epi.h"

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
    double Benth_resp;
    double MAaA;
    double* MAgrow;
    double* MArate;
    double MAlen;
    double KN;
    double KP;

    /*
     * epis
     */
    int SG_N_i;
    int SG_N_pr_i;
    int SG_N_gr_i;
    int MA_N_i;
    int Epilight_i;
    int Epilightatt_i;
    int EpiTN_i;
    int EpiTP_i;
    int EpiTC_i;

    /*
     * tracers
     */
    int NH4_sed_i;
    int NO3_sed_i;
    int DIP_sed_i;
    int DIC_sed_i;
    int Oxygen_wc_i;
    int Oxy_pr_wc_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int umax_i;
} workspace;

void seagrass_grow_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    int OFFSET_SED = tracers->n;
    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "SGumax");
    ws->aA = get_parameter_value(e, "SGaA");
    ws->m = get_parameter_value(e, "SGm");
    ws->Benth_resp = get_parameter_value(e, "Benth_resp") * atk_A_I;
    ws->MAaA = get_parameter_value(e, "MAaA");
    {
        int MAn = get_parameter_value(e, "MAn");
        char* MAtable = get_parameter_stringvalue(e, "MAtable");

        ws->MAlen = extract_grow(MAn, MAtable, &ws->MAgrow, &ws->MArate);
    }
    ws->KN = get_parameter_value(e, "SG_KN");
    ws->KP = get_parameter_value(e, "SG_KP");

    /*
     * epis
     */
    ws->SG_N_i = e->find_index(epis, "SG_N", e) + OFFSET_EPI;
    ws->SG_N_pr_i = e->find_index(epis, "SG_N_pr", e) + OFFSET_EPI;
    ws->SG_N_gr_i = e->find_index(epis, "SG_N_gr", e) + OFFSET_EPI;
    ws->MA_N_i = e->find_index(epis, "MA_N", e) + OFFSET_EPI;
    ws->Epilight_i = e->find_index(epis, "Epilight", e) + OFFSET_EPI;
    ws->Epilightatt_i = e->find_index(epis, "Epilightatt", e) + OFFSET_EPI;
    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;

    /*
     * tracers
     */
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_sed_i = e->find_index(tracers, "NH4", e) + OFFSET_SED;
    ws->NO3_sed_i = e->find_index(tracers, "NO3", e) + OFFSET_SED;
    ws->DIP_sed_i = e->find_index(tracers, "DIP", e) + OFFSET_SED;
    ws->DIC_sed_i = e->find_index(tracers, "DIC", e) + OFFSET_SED;



    /*non essential diagnosic tracer*/

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "SGumax", e);
}

void seagrass_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;
}

void seagrass_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->MAgrow);
    free(ws->MArate);
    free(ws);
}

void seagrass_grow_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double SG_N = y[ws->SG_N_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    y[ws->Epilightatt_i] += SG_N * ws->aA;

    if (ws->do_mb) {
        y[ws->EpiTN_i] += SG_N;
        y[ws->EpiTP_i] += SG_N * atk_W_P;
        y[ws->EpiTC_i] += SG_N * atk_W_C;
    }
    if(c->porosity == 0.0)
      p->ecology->quitfn("ecology:seagrass_grow_epi: porosity cannot be 0.0!");
}

void seagrass_grow_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double dz_sed = c->dz_sed;
    double dz_wc = c->dz_wc;
    double porosity = c->porosity;

    double NO3_sed = y[ws->NO3_sed_i];
    double NH4_sed = y[ws->NH4_sed_i];
    double din_sed = NO3_sed + NH4_sed;
    double DIN_sed = (din_sed > 0.0) ? din_sed : EPS_DIN;
    double dip_sed = y[ws->DIP_sed_i];
    double DIP_sed = (dip_sed > 0.0) ? dip_sed : EPS_DIP;
    double Epilight = y[ws->Epilight_i];
    double I = 2.77e18 * Epilight / AV;
    double SG_N = y[ws->SG_N_i];
    double MA_N = y[ws->MA_N_i];
    double umax = cv[ws->umax_i];
    double kI = I * (1.0 - exp(-SG_N * ws->aA)) * exp(-MA_N * ws->MAaA);
    double kN = DIN_sed * umax * atk_A_N * ws->m * SG_N / ws->KN;
    double kP = DIP_sed * umax * atk_A_P * ws->m * SG_N / ws->KP;
    double up[3];

    up[0] = kI;
    up[1] = kN;
    up[2] = kP;

    {
        double growthrate = CRgrowth(umax, up, ws->m * SG_N, ws->MAgrow, ws->MArate, ws->MAlen, 3, ws->Benth_resp);
        double growth = SG_N * growthrate;
        double Oxy_pr = growth * atk_W_O / dz_wc;

        y1[ws->SG_N_i] += growth;
        y1[ws->SG_N_gr_i] += growthrate * SEC_PER_DAY; /* d-1 */
        y1[ws->SG_N_pr_i] += growth * SEC_PER_DAY * atk_W_C; /* gC m-2 d-1 */;
        y1[ws->NH4_sed_i] -= growth * NH4_sed / DIN_sed / dz_sed / porosity;
        y1[ws->NO3_sed_i] -= growth * NO3_sed / DIN_sed / dz_sed / porosity;
        y1[ws->DIP_sed_i] -= growth * atk_W_P / dz_sed / porosity;
        y1[ws->DIC_sed_i] -= growth * atk_W_C / dz_sed / porosity;
        y1[ws->Oxygen_wc_i] += Oxy_pr;


	if (ws->  Oxy_pr_wc_i> -1)
	  y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY;
    }
}

void seagrass_grow_epi_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double SG_N = y[ws->SG_N_i];

    y[ws->EpiTN_i] += SG_N;
    y[ws->EpiTP_i] += SG_N * atk_W_P;
    y[ws->EpiTC_i] += SG_N * atk_W_C;
}
