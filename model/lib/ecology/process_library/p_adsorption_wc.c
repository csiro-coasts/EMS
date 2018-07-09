/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/p_adsorption_wc.c
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
 *  $Id: p_adsorption_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "p_adsorption_wc.h"
#include "einterface.h"

// #define EFI_EPS 1.0e-4

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    
    double Pads_r_t0;
    double Pads_Kwc;
    double Pads_KO;
    double Pads_exp;

    /*
     * tracers
     */
    int PIP_i;
    int PIPI_i;
    int Oxygen_i;
    int DIP_i;
    int TP_i;
    int EFI_i;
    int Mud_i;

    /*
     * common variables
     */
    int Tfactor_i;
    int Pads_r_i;

} workspace;

void p_adsorption_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->Pads_r_t0 = get_parameter_value(e, "Pads_r");
    ws->Pads_Kwc = get_parameter_value(e, "Pads_Kwc");
    ws->Pads_KO = get_parameter_value(e, "Pads_KO");
    ws->Pads_exp = get_parameter_value(e, "Pads_exp");

    /*
     * tracers
     */
    ws->EFI_i = e->find_index(e->tracers, "EFI", e);
    ws->PIP_i = e->find_index(tracers, "PIP", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->TP_i = e->find_index(tracers, "TP", e);

    ws->Mud_i = e->try_index(e->tracers, "Mud", e);

    /* include PIPI in TP, even though no change in this processes */

    ws->PIPI_i = e->try_index(tracers, "PIPI", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->Pads_r_i = find_index_or_add(e->cv_cell, "Pads_r", e);
}

void p_adsorption_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    /* test for equal sinking rates of PIP and EFI components */

    double v1,v2,v3;
      
    /*
     * Not valid during a pre_build (RECOM)
     */

    if (!e->pre_build) {
      if (ws->Mud_i > -1){
	v1 = einterface_gettracersvel(e->model,"PIP");
	v2 = einterface_gettracersvel(e->model,"Mud");
	v3 = einterface_gettracersvel(e->model,"FineSed");
	
	if ((v1!=v2)||(v1!=v3)){
	  e->quitfn("Sinking rates of PIP %e, Mud %e, FineSed %e are not equal. Must be fixed in .prm or .tran file",v1,v2,v3);
	}
      }
    }

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void p_adsorption_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void p_adsorption_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor  = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double PIPI = 0.0;

    cv[ws->Pads_r_i] = ws->Pads_r_t0 * Tfactor;


    if (ws->do_mb) {

      if (ws->PIPI_i > -1)
        PIPI =  y[ws->PIPI_i];

        y[ws->TP_i] += y[ws->PIP_i] + PIPI;
    }
}

void p_adsorption_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Oxygen = e_max(y[ws->Oxygen_i]);
    double DIP = y[ws->DIP_i];
    double EFI = y[ws->EFI_i];
    double PIP = y[ws->PIP_i];
    // double npow = 1.0 / ws->Pads_exp;
    double Pads_r = cv[ws->Pads_r_i];
    
    /*    double adsorption = Pads_r * EFI * ws->Pads_Kwc * DIP * (Oxygen / (ws->Pads_KO + e_max(Oxygen)));
    
	  double net_desorp = Pads_r * EFI * ws->Pads_Kwc * (pow(e_max(PIP) / (e_max(EFI + EFI_EPS) * ws->Pads_Kwc), npow)) - adsorption; */
 
    double net_desorp  = Pads_r * (PIP / (e_max(EFI) * ws->Pads_Kwc) -  DIP * (Oxygen / (ws->Pads_KO + Oxygen)));

    y1[ws->DIP_i] += net_desorp;
    y1[ws->PIP_i] -= net_desorp;
}

void p_adsorption_wc_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = ((cell*) pp);
    double* y = c->y;

    double PIPI = 0.0;
 
    if (ws->PIPI_i > -1)
        PIPI =  y[ws->PIPI_i];

    y[ws->TP_i] += y[ws->PIP_i] + PIPI;
    
}

