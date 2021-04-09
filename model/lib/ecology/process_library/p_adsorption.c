/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/p_adsorption.c
 *  
 *  Description: Carries out phosphorus adsorption-desorption in both water column and sediments, 
 *               and considers two size classes of sediments (Dust and non-Dust)
 * 
 *  Immobilisation only occurs in the sediments.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: p_adsorption.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "p_adsorption.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    
    double Pads_r_t0;
    double Pads_K;
    double Pads_KO;
    double Pads_exp;
    double r_immob_PIP_t0;

    /*
     * tracers
     */
    int PIP_i;
    int PIP_Dust_i;
    int PIPI_i;
    int DIP_i;
    int Oxygen_i;
    int TP_i;
    int EFI_i;

    int Mud_i;
    int Mud_carbonate_i;
    int Dust_i;
  
    /*
     * common variables
     */
    int Tfactor_i;
    int Pads_r_i;
    int r_immob_PIP_i;

} workspace;

void p_adsorption_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */

    ws->Pads_r_t0 = try_parameter_value(e,"Pads_r");
    if (isnan(ws->Pads_r_t0)){
      ws->Pads_r_t0 = 4.000000e-02/86400.0;
      eco_write_setup(e,"Code default of Pads_r = %e \n",ws->Pads_r_t0);
    }
    
    ws->Pads_KO = try_parameter_value(e,"Pads_KO");
    if (isnan(ws->Pads_KO)){
      ws->Pads_KO = 2000.0;
      eco_write_setup(e,"Code default of Pads_KO = %e \n",ws->Pads_KO);
    }
    
    ws->Pads_exp = try_parameter_value(e,"Pads_exp");
    if (isnan(ws->Pads_exp)){
      ws->Pads_exp = 1.0;
      eco_write_setup(e,"Code default of Pads_exp = %e \n",ws->Pads_exp);
    }

    // ws->Pads_r_t0 = get_parameter_value(e, "Pads_r");
    // ws->Pads_KO = get_parameter_value(e, "Pads_KO");
    // ws->Pads_exp = get_parameter_value(e, "Pads_exp");
    // ws->r_immob_PIP_t0 = get_parameter_value(e, "r_immob_PIP");

    if (p->type == PT_SED){
      //ws->Pads_K = get_parameter_value(e, "Pads_Ksed");
      ws->Pads_K = try_parameter_value(e,"Pads_Ksed");
      if (isnan(ws->Pads_K)){
	ws->Pads_K = 74.0;
	eco_write_setup(e,"Code default of Pads_Ksed = %e \n",ws->Pads_K);
      }

      ws->r_immob_PIP_t0 = try_parameter_value(e,"r_immob_PIP");
      if (isnan(ws->r_immob_PIP_t0)){
	ws->r_immob_PIP_t0 = 0.0;
	eco_write_setup(e,"Code default of r_immob_PIP_t0 = %e \n",ws->r_immob_PIP_t0);
      }
      
    }else{
      // ws->Pads_K = get_parameter_value(e, "Pads_Kwc");
      ws->Pads_K = try_parameter_value(e,"Pads_Kwc");
      if (isnan(ws->Pads_K)){
	ws->Pads_K = 30.0;
	eco_write_setup(e,"Code default of Pads_Kwc = %e \n",ws->Pads_K);
      }
      ws->r_immob_PIP_t0 = 0.0;
    }
    /*
     * tracers
     */
    ws->EFI_i = e->find_index(e->tracers, "EFI", e);
    ws->Dust_i = e->find_index(e->tracers, "Dust", e);
    ws->PIP_i = e->find_index(tracers, "PIP", e);
    ws->PIPI_i = e->find_index(tracers, "PIPI", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->TP_i = e->find_index(tracers, "TP", e);

    ws->Mud_i = e->try_index(e->tracers, "Mud", e);
    ws->Mud_carbonate_i = e->try_index(e->tracers, "Mud-carbonate", e);

    ws->PIP_Dust_i = e->find_index(tracers, "PIP_Dust", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->Pads_r_i = find_index_or_add(e->cv_cell, "Pads_r", e);
    ws->r_immob_PIP_i = find_index_or_add(e->cv_cell, "r_immob_PIP", e);
}

void p_adsorption_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

/* test for equal sinking rates of PIP and EFI components */

    double v1,v2,v3,v4;
      
    /*
     * Not valid during a pre_build (RECOM)
     */

    if (!e->pre_build) {
      if (p->type == PT_WC){
	if (ws->Mud_i > -1){
	  v1 = einterface_gettracersvel(e->model,"PIP");
	  v2 = einterface_gettracersvel(e->model,"Mud");
	  v3 = einterface_gettracersvel(e->model,"FineSed");

	  eco_write_setup(e,"\np_adsorption: sinking rates of PIP %e, Mud %e, FineSed %e \n",v1,v2,v3);
	  
	  if ((v1!=v2)||(v1!=v3)){
	    e->quitfn("Sinking rates of PIP %e, Mud %e, FineSed %e are not equal. Must be fixed in .prm or .tran file",v1,v2,v3);
	  }
	}
	if (ws->Mud_carbonate_i > -1){
	  v1 = einterface_gettracersvel(e->model,"PIP");
	  v2 = einterface_gettracersvel(e->model,"Mud-carbonate");
	  v3 = einterface_gettracersvel(e->model,"Mud-mineral");
	  v4 = einterface_gettracersvel(e->model,"FineSed");

	  eco_write_setup(e,"\np_adsorption: sinking rates of PIP = %e, Mud-carbonate = %e, Mud-mineral = %e FineSed = %e \n",v1,v2,v3,v4);
	  
	  if ((v1!=v2)||(v1!=v3)||(v1!=v4)){
	    e->quitfn("Sinking rates of PIP %e, Mud-carbonate %e, Mud-mineral %e FineSed %e are not equal. Must be fixed in .prm or .tran file",v1,v2,v3,v4);
	  }
	}
	if (ws->PIP_Dust_i > -1){
	  v1 = einterface_gettracersvel(e->model,"PIP_Dust");
	  v2 = einterface_gettracersvel(e->model,"Dust");

	  eco_write_setup(e,"\np_adsorption: sinking rates of PIP_Dust = %e, Dust = %e \n",v1,v2);
	  
	  if (v1!=v2){
	    e->quitfn("Sinking rates of PIP_Dust %e, Dust %e,are not equal. Must be fixed in .prm or .tran file",v1,v2);
	  }
	}
      }
      if (p->type == PT_WC){
	ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
      }else{
	ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
      }
    }
}
void p_adsorption_destroy(eprocess* p)
{
    free(p->workspace);
}

void p_adsorption_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->Pads_r_i] = ws->Pads_r_t0 * Tfactor;
    cv[ws->r_immob_PIP_i] = ws->r_immob_PIP_t0 * Tfactor;

    if (ws->do_mb) {
        double* y = c->y;
        double PIP = y[ws->PIP_i];
	double PIP_Dust = y[ws->PIP_Dust_i];
        double PIPI = y[ws->PIPI_i];

        y[ws->TP_i] += PIP + PIPI + PIP_Dust;
    }
}

void p_adsorption_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double porosity = c->porosity;

    double Pads_r = cv[ws->Pads_r_i];
    double r_immob_PIP = cv[ws->r_immob_PIP_i];

    double Oxygen = e_max(y[ws->Oxygen_i]);
    double DIP = y[ws->DIP_i];
    double EFI = y[ws->EFI_i];
    double PIP = y[ws->PIP_i];

    // New code

    EFI = max(1.0e-6,y[ws->EFI_i]-y[ws->Dust_i]);
    double Dust = max(1.0e-6,y[ws->Dust_i]);
    double PIP_Dust = y[ws->PIP_Dust_i];

    //

    double net_desorp  = Pads_r * (PIP / (EFI * ws->Pads_K) -  DIP * (Oxygen / (ws->Pads_KO + Oxygen)));
       
    double PIP_immob = r_immob_PIP * PIP;    
    
    y1[ws->DIP_i] += (net_desorp) / porosity;
    y1[ws->PIP_i] += -net_desorp - PIP_immob;
    y1[ws->PIPI_i] += PIP_immob;

    // New code

    double net_desorp_Dust = Pads_r * ( PIP_Dust / (Dust * ws->Pads_K) - DIP * (Oxygen / (ws->Pads_KO + Oxygen)));
       
    double PIP_immob_Dust = r_immob_PIP * PIP_Dust;
    
    y1[ws->DIP_i] += (net_desorp_Dust) / porosity;
    y1[ws->PIP_Dust_i] += -net_desorp_Dust - PIP_immob_Dust;
    y1[ws->PIPI_i] += PIP_immob_Dust;

    //
}

void p_adsorption_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = ((cell*) pp);
    double* y = c->y;

    if (ws->do_mb) {
      double PIP = y[ws->PIP_i];
      double PIPI = y[ws->PIPI_i];
      double PIP_Dust = y[ws->PIP_Dust_i];
      y[ws->TP_i] += PIP + PIPI + PIP_Dust;
    }
}
