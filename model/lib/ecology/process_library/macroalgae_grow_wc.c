/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_grow_wc.c
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
 *  $Id: macroalgae_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "macroalgae_grow_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
    int do_mb;                  /* flag */
   
    

    /*
     * parameters
     */
    double umax;
    double Vmax_NH4;
    double Vmax_NO3;
    double Khalf_NH4;
    double Khalf_NO3;
    double aA;
    double m;
    int n;
    double MA_WC_resp;
    double h;
    double Qmax;
    double Qmin;
    double Kc;
    double Topt;
    double Tran;
    double Is;
 

    /*
     * tracers
     */
    int temp_i;
    int SWS_N_i;
    int SWF_N_i;
    int SWS_N_gr_i;
    int SWF_N_gr_i;
    int NO3_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Oxy_pr_i;
    int Light_i;
    int Chl_a_i;
    int Kd_i;
    int Salt_i;
    int TN_i;
    int TP_i;
    int TC_i;

    /*
     * common cell variables
     */
   
    int DNO3_i;
    int DPO4_i;
    int DNH4_i;
   
} workspace;

void macroalgae_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    char* prm = p->prms->se[0]->s;
    

    p->workspace = ws;

   
    /*
     * parameters
     */
    ws->umax = get_parameter_value(e, "MAumax");
    ws->aA = try_parameter_value(e, "MAaA");
    ws->m = try_parameter_value(e, "MAm");
    ws->n = rint(get_parameter_value(e,"MAn"));
    ws->h = get_parameter_value(e, "h");
    ws->Qmin = get_parameter_value(e, "Qmin");
    ws->Qmax = get_parameter_value(e, "Qmax");
    ws->Kc = get_parameter_value(e, "Kc");
    ws->Vmax_NH4 = get_parameter_value(e, "Vmax_NH4");
    ws->Vmax_NO3 = get_parameter_value(e, "Vmax_NO3");
    ws->Khalf_NH4 = get_parameter_value(e, "Khalf_NH4");
    ws->Khalf_NO3 = get_parameter_value(e, "Khalf_NO3");
    ws->Tran = get_parameter_value(e, "Tran");
    ws->Topt = get_parameter_value(e, "Topt");
    ws->Is = get_parameter_value(e, "Is");
    ws->MA_WC_resp = get_parameter_value(e, "MA_WC_resp") * atk_A_I;
   
    /*
     * tracers
     */
    ws->temp_i = e->find_index(tracers,"temp", e);
    ws->SWS_N_i = e->find_index(tracers,"SWS_N", e);
    ws->SWF_N_i = e->find_index(tracers,"SWF_N", e);
    ws->SWS_N_gr_i = e->find_index(tracers, "SWS_N_gr", e);
    ws->SWF_N_gr_i = e->find_index(tracers, "SWF_N_gr", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Kd_i = e->find_index(tracers, "Kd", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->Salt_i = e->find_index(tracers, "Salt", e);
    ws->Oxy_pr_i = e->find_index(tracers, "Oxy_pr", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
   
    /*
     * common cell variables
     */
   
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
   
}

void macroalgae_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void macroalgae_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;
    free(ws);
}

void macroalgae_grow_wc_precalc(eprocess* p, void* pp)
{
  
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
    
  
    double SWS_N = y[ws->SWS_N_i];
    double SWF_N = y[ws->SWF_N_i];
    
   
    y[ws->Kd_i] += ((ws->aA)*(SWF_N)) * 1000; //change units from gm to mg
     
   


    if (ws->do_mb) {
      y[ws->TN_i] += (SWS_N + SWF_N) * 1000;  //change units to mgN m-3
        y[ws->TP_i] += SWF_N * atk_W_P * 1000;
        y[ws->TC_i] += SWF_N * atk_W_C * 1000;
    }
}

void macroalgae_grow_wc_calc(eprocess* p, void* pp)
{
  //Here I put in my calculations so stuff from this cell into the next
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double h = c->dz_wc;
    double* y = ia->y;
    double* y1 = ia->y1;
   
    double SWS_N = y[ws->SWS_N_i];
    double SWF_N = y[ws->SWF_N_i];
    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    
    double Iz = y[ws->Light_i];
    double I =(2.77e18 * Iz / AV)*1e6;
    double diff[2];
    double conc[2];

    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];
    conc[0] = DIN * mgN2molN;
    conc[1] = DIP * mgP2molP;

    {
        double umax = ws->umax;
        double Qmax = ws->Qmax;
        double Qmin = ws->Qmin;
        double Kc = ws->Kc;
        double Topt = ws->Topt;
        double Tran = ws->Tran;
        double temp = y[ws->temp_i];
        double k = y[ws->Kd_i];
        double Is = ws->Is;
        double Vmax_NH4 = ws->Vmax_NH4;
        double Vmax_NO3 = ws->Vmax_NO3;
        double Khalf_NH4 = ws->Khalf_NH4;
        double Khalf_NO3 = ws->Khalf_NO3;
        double B = SWF_N/Qmin;
        double Q =  Qmin*(1.0 + (SWS_N/SWF_N));
        double up_NH4 = ((Vmax_NH4*NH4)/(Khalf_NH4 + NH4))*((Qmax - Q)/(Qmax - Qmin));
        double up_NO3 = ((Vmax_NO3*NO3)/(Khalf_NO3 + NO3))*((Qmax - Q)/(Qmax - Qmin));
        double Tlim = 1.0/(1.0 + exp(-((temp-Topt)/Tran)));
        double Qlim = (Q-Qmin)/(Q-Kc);
	double Ilim = (exp(1.0)/(k*h))*(exp(-(I*exp(-k*h))/Is)-exp(-I/Is));
        double dic = y[ws->DIC_i];
        double up_C = dic/(dic+1000);
        double growthrate = umax*Ilim*Tlim*Qlim*up_C;
        double growth = SWS_N * growthrate;
	// printf("up_C: %f %f \n", up_C, growth);


        y1[ws->SWS_N_i] += (up_NH4 + up_NO3) * B - growth;
        y1[ws->SWF_N_i] += growth;
        y1[ws->NH4_i] -=up_NH4 * B *1000;
        y1[ws->NO3_i] -= up_NO3 * B * 1000;
        y1[ws->DIP_i] -= growth * atk_W_P * 1000;
        y1[ws->DIC_i] -= growth * atk_W_C * 1000;
        y1[ws->Oxygen_i] += growth * atk_W_O *1000;
	// y1[ws->SWS_N_gr_i] = (up_NH4 + up_NO3) * SEC_PER_DAY;
	// y1[ws->SWF_N_gr_i] = growthrate * SEC_PER_DAY;

	// y1[ws->SWF_N_pr_i] += growth * SEC_PER_DAY * atk_W_C;
	// y1[ws->SWF_N_gr_i] = growthrate * SEC_PER_DAY;
      
        y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY * 1000;
    }
}

void macroalgae_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double SWS_N = y[ws->SWS_N_i];
    double SWF_N = y[ws->SWF_N_i];
    // y[ws->Chl_a_i] += SWF_N / ws->NtoCHL;
    y[ws->TN_i] += (SWS_N + SWF_N) * 1000;
    y[ws->TP_i] += SWF_N * atk_W_P * 1000;
    y[ws->TC_i] += SWF_N * atk_W_C * 1000;
}
