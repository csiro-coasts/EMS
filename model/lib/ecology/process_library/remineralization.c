/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/remineralization.c
 *  
 *  Description:
 *  
 *  - Remineralisation processes affecting C, N, P and O budgets.
 *  - Works in both water column and sediments (porewaters).
 *  - Does mass balance terms for all particulate organic and inorganic nutrients.
 *  
 *  19/12/18 MEB - changed COD / DO interaction term.
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: remineralization.c 6450 2020-01-22 04:15:06Z wil00y $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "remineralization.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double r_DetPL_t0;
    double r_DetBL_t0;
    double r_RD_t0;
    double r_DOM_t0;
    double F_LD_RD;
    double F_LD_DOM;
    double F_RD_DOM;
    double KO_aer;
    double r_RD_NtoP;
    double r_DOM_NtoP;
    double tau_COD;

    /*
     * tracers
     */
    int DetPL_N_i;
    int DetBL_N_i;
    int DetR_N_i;
    int DetR_C_i;
    int DetR_P_i;
    int DOR_N_i;
    int DOR_C_i;
    int DOR_P_i;
    int NH4_i;
    int NH4_pr_i;
    int Oxygen_i;
    int Oxy_pr_i;
    int DIC_i;
    int DIP_i;
    int NO3_i;
    int TN_i;
    int TC_i;
    int TP_i;
    int BOD_i;
    int COD_i;

    /*
     * common variables
     */
    int Tfactor_i;
    int r_DetPL_i;
    int r_DetBL_i;
    int r_RD_i;
    int r_DOM_i;
    int NH4_remin_pr_i;         /* used in sediment cells only */
} workspace;

void remineralization_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->r_DetPL_t0 = get_parameter_value(e, "r_DetPL");
    ws->r_DetBL_t0 = get_parameter_value(e, "r_DetBL");
    ws->r_RD_t0 = get_parameter_value(e, "r_RD");
    ws->r_DOM_t0 = get_parameter_value(e, "r_DOM");
    ws->F_LD_RD = get_parameter_value(e, "F_LD_RD");
    ws->F_LD_DOM = get_parameter_value(e, "F_LD_DOM");
    if (ws->F_LD_RD + ws->F_LD_DOM > 1.0)
        e->quitfn("ecology: error: F_LD_RD + F_LD_DOM > 1\n");
    ws->F_RD_DOM = get_parameter_value(e, "F_RD_DOM");
    if (ws->F_RD_DOM > 1.0)
        e->quitfn("ecology: error: F_RD_DOM > 1\n");
    ws->KO_aer = get_parameter_value(e, "KO_aer");

    ws->r_RD_NtoP = try_parameter_value(e, "r_RD_NtoP");
    if (isnan(ws->r_RD_NtoP)){
      ws->r_RD_NtoP = 1.0;
      eco_write_setup(e,"Code default of r_RD_NtoP = %e \n",ws->r_RD_NtoP);
    }

    ws->r_DOM_NtoP = try_parameter_value(e, "r_DOM_NtoP");
    if (isnan(ws->r_DOM_NtoP)){
      ws->r_DOM_NtoP = 1.0;
      eco_write_setup(e,"Code default of r_DOM_NtoP = %e \n",ws->r_DOM_NtoP);
    }

    ws->tau_COD = try_parameter_value(e, "tau_COD");
    if (isnan(ws->tau_COD)){
      ws->tau_COD = 1.0/3600.0;
      eco_write_setup(e,"Code default of tau_COD = %e \n",ws->tau_COD);
    }

    /*
     * tracers
     */
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->DetBL_N_i = e->find_index(tracers, "DetBL_N", e);
    ws->DetR_N_i = e->find_index(tracers, "DetR_N", e);
    ws->DetR_C_i = e->find_index(tracers, "DetR_C", e);
    ws->DetR_P_i = e->find_index(tracers, "DetR_P", e);
    ws->DOR_N_i = e->find_index(tracers, "DOR_N", e);
    ws->DOR_C_i = e->find_index(tracers, "DOR_C", e);
    ws->DOR_P_i = e->find_index(tracers, "DOR_P", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);

    ws->COD_i = e->try_index(tracers, "COD", e);

    /*non essential diagnostic tracer*/
    
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->r_DetPL_i = find_index_or_add(e->cv_cell, "r_DetPL", e);
    ws->r_DetBL_i = find_index_or_add(e->cv_cell, "r_DetBL", e);
    ws->r_RD_i = find_index_or_add(e->cv_cell, "r_RD", e);
    ws->r_DOM_i = find_index_or_add(e->cv_cell, "r_DOM", e);
    ws->NH4_remin_pr_i = (p->type == PT_SED) ? find_index_or_add(e->cv_cell, "NH4_remin_pr", e) : -1;
}

void remineralization_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, (p->type == PT_WC) ? "massbalance_wc" : "massbalance_sed", e) >= 0) ? 1 : 0;
}

void remineralization_destroy(eprocess* p)
{
    free(p->workspace);
}

void remineralization_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double DOR_N = y[ws->DOR_N_i];
    double DetPL_N = y[ws->DetPL_N_i];
    double DetBL_N = y[ws->DetBL_N_i];
    double DetR_N = y[ws->DetR_N_i];

    cv[ws->r_DetPL_i] = Tfactor * ws->r_DetPL_t0;
    cv[ws->r_DetBL_i] = Tfactor * ws->r_DetBL_t0;
    cv[ws->r_RD_i] = Tfactor * ws->r_RD_t0;
    cv[ws->r_DOM_i] = Tfactor * ws->r_DOM_t0;

    if(c->porosity == 0.0)
      p->ecology->quitfn("ecology:remineralization: porosity cannot be 0.0!");

    if (ws->do_mb) {
        double porosity = (c->type == CT_SED) ? c->porosity : 1.0;
        double NH4 = y[ws->NH4_i];
        double NO3 = y[ws->NO3_i];
        double DIC = y[ws->DIC_i];
        double DIP = y[ws->DIP_i];
        double DOR_C = y[ws->DOR_C_i];
        double DetR_C = y[ws->DetR_C_i];
        double DetR_P = y[ws->DetR_P_i];
        double DOR_P = y[ws->DOR_P_i];
	double BOD = y[ws->BOD_i];

        y[ws->TN_i] += (NH4 + NO3 + DOR_N) * porosity + DetPL_N + DetBL_N + DetR_N;
        y[ws->TC_i] += (DIC + DOR_C) * porosity + DetPL_N * red_W_C + DetBL_N * atk_W_C + DetR_C;
        y[ws->TP_i] += (DIP + DOR_P) * porosity + DetPL_N * red_W_P + DetBL_N * atk_W_P + DetR_P;
	y[ws->BOD_i] += (DOR_C * C_O_W * porosity + DetPL_N * red_W_O + DetBL_N * atk_W_O + DetR_C * C_O_W);
    }

}

void remineralization_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double porosity = (c->type == CT_SED) ? c->porosity : 1.0;

    /*
     * multiple for "*_pr" tracers in sediment
     */
    double dz_pr = (c->type == CT_SED) ? c->dz_sed : 1.0;

    double DetPL_N = y[ws->DetPL_N_i];
    double DetBL_N = y[ws->DetBL_N_i];
    double DetR_C = y[ws->DetR_C_i];
    double DetR_N = y[ws->DetR_N_i];
    double DetR_P = y[ws->DetR_P_i];
    double DOR_N = y[ws->DOR_N_i];
    double DOR_C = y[ws->DOR_C_i];
    double DOR_P = y[ws->DOR_P_i];
    double Oxygen2 = e_max(y[ws->Oxygen_i]);

    Oxygen2 = Oxygen2 * Oxygen2; 

    /* limit detrital remineralisation to concentrations > 0.001 */
    double DetPL_N_break = 0.0;
    if (DetPL_N > 0.001){
      DetPL_N_break = cv[ws->r_DetPL_i] * DetPL_N;}

    double DetBL_N_break = 0.0;
    if (DetBL_N > 0.001){
      DetBL_N_break = cv[ws->r_DetBL_i] * DetBL_N;}

    double r_RD = cv[ws->r_RD_i];

    double DetR_C_break = 0.0;
    if (DetR_C > 0.001){
      DetR_C_break = r_RD * DetR_C;}

    double DetR_N_break = 0.0;
    if (DetR_N > 0.001){
      DetR_N_break = r_RD * DetR_N;}

    double DetR_P_break = 0.0;
    if (DetR_P > 0.001){
      DetR_P_break = ws->r_RD_NtoP * r_RD * DetR_P;}

    double DetPL_N_remin = DetPL_N_break * (1.0 - ws->F_LD_RD - ws->F_LD_DOM);
    double DetBL_N_remin = DetBL_N_break * (1.0 - ws->F_LD_RD - ws->F_LD_DOM);
    double DetR_N_remin = DetR_N_break * (1.0 - ws->F_RD_DOM);
    double DetR_C_remin = DetR_C_break * (1.0 - ws->F_RD_DOM);
    double DetR_P_remin = DetR_P_break * (1.0 - ws->F_RD_DOM);
    double r_DOM = cv[ws->r_DOM_i];
    double DO_N_remin = r_DOM * DOR_N * porosity;
    double DO_C_remin = r_DOM * DOR_C * porosity;
    double DO_P_remin = ws->r_DOM_NtoP * r_DOM * DOR_P * porosity;
    double NH4_pr = DetPL_N_remin + DetBL_N_remin + DetR_N_remin + DO_N_remin;

    /* split oxic and anoxic remineralisation */

    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    double remin_oxy = (DetPL_N_remin * red_W_O + DetBL_N_remin * atk_W_O + (DetR_C_remin + DO_C_remin) * C_O_W);

    double Oxy_pr = - remin_oxy * sigmoid;

    y1[ws->DetPL_N_i] -= DetPL_N_break;
    y1[ws->DetBL_N_i] -= DetBL_N_break;
    y1[ws->DetR_C_i] += DetPL_N_break * ws->F_LD_RD * red_W_C + DetBL_N_break * ws->F_LD_RD * atk_W_C - DetR_C_break;
    y1[ws->DetR_N_i] += DetPL_N_break * ws->F_LD_RD + DetBL_N_break * ws->F_LD_RD - DetR_N_break;
    y1[ws->DetR_P_i] += DetPL_N_break * ws->F_LD_RD * red_W_P + DetBL_N_break * ws->F_LD_RD * atk_W_P - DetR_P_break;
    y1[ws->DOR_N_i] += (DetPL_N_break * ws->F_LD_DOM + DetBL_N_break * ws->F_LD_DOM + DetR_N_break * ws->F_RD_DOM - DO_N_remin) / porosity;
    y1[ws->DOR_C_i] += (DetPL_N_break * ws->F_LD_DOM * red_W_C + DetBL_N_break * ws->F_LD_DOM * atk_W_C + DetR_C_break * ws->F_RD_DOM - DO_C_remin) / porosity;
    y1[ws->DOR_P_i] += (DetPL_N_break * ws->F_LD_DOM * red_W_P + DetBL_N_break * ws->F_LD_DOM * atk_W_P + DetR_P_break * ws->F_RD_DOM - DO_P_remin) / porosity;
    y1[ws->NH4_i] += NH4_pr / porosity;
    y1[ws->Oxygen_i] += Oxy_pr / porosity;
    y1[ws->DIC_i] += (DetPL_N_remin * red_W_C + DetBL_N_remin * atk_W_C + DetR_C_remin + DO_C_remin) / porosity;
    y1[ws->DIP_i] += (DetPL_N_remin * red_W_P + DetBL_N_remin * atk_W_P + DetR_P_remin + DO_P_remin) / porosity;

    if (ws->COD_i > -1){

      y1[ws->COD_i] += remin_oxy * (1.0-sigmoid) / porosity;

    /* Now consume oxygen with COD */

      // double tau_oxy = 1.0/3600.0;   // consume oxygen with COD within an hour.    
      // y1[ws->Oxygen_i] -= tau_oxy * y[ws->COD_i] * (y[ws->Oxygen_i] / 8000.0); 
      // y1[ws->COD_i] -= tau_oxy * y[ws->COD_i] * (y[ws->Oxygen_i] / 8000.0);

      double consume = ws->tau_COD * min(y[ws->COD_i],8000.0) * (y[ws->Oxygen_i] / 8000.0);
      y1[ws->Oxygen_i] -= consume;
      y1[ws->COD_i] -= consume;
    }

 if (ws-> NH4_pr_i > -1)
   y1[ws->NH4_pr_i] += NH4_pr * SEC_PER_DAY * dz_pr;

 if (ws-> Oxy_pr_i > -1)
   y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY * dz_pr;

 if (p->type == PT_SED)
        cv[ws->NH4_remin_pr_i] = NH4_pr;
}

void remineralization_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = ((cell*) pp);
    double* y = c->y;

    double porosity = (c->type == CT_SED) ? c->porosity : 1.0;

    double NH4 = y[ws->NH4_i];
    double NO3 = y[ws->NO3_i];
    double DetPL_N = y[ws->DetPL_N_i];
    double DetBL_N = y[ws->DetBL_N_i];
    double DetR_N = y[ws->DetR_N_i];
    double DIC = y[ws->DIC_i];
    double DIP = y[ws->DIP_i];
    double DetR_C = y[ws->DetR_C_i];
    double DetR_P = y[ws->DetR_P_i];
    double DOR_N = y[ws->DOR_N_i];
    double DOR_C = y[ws->DOR_C_i];
    double DOR_P = y[ws->DOR_P_i];

    y[ws->TN_i] += (NH4 + NO3 + DOR_N) * porosity + DetPL_N + DetBL_N + DetR_N;
    y[ws->TC_i] += (DIC + DOR_C) * porosity + DetPL_N * red_W_C + DetBL_N * atk_W_C + DetR_C;
    y[ws->TP_i] += (DIP + DOR_P) * porosity + DetPL_N * red_W_P + DetBL_N * atk_W_P + DetR_P;
    y[ws->BOD_i] += (DOR_C * C_O_W * porosity + DetPL_N * red_W_O + DetBL_N * atk_W_O + DetR_C * C_O_W);
}
