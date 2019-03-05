/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dinoflagellate_spectral_mortality_wc.c
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
 *  $Id: dinoflagellate_spectral_mortality_wc.c 6151 2019-03-05 02:00:45Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "dinoflagellate_spectral_mortality_wc.h"

typedef struct {
  int do_mb;                  /* flag */

    /*
     * parameters
     */
  double mL_t0;
  double KO_aer;
    /*
     * tracers
     */
  int Phy_N_i;
  int Phy_NR_i;
  int Phy_I_i;
  int DIC_i;
  int DIP_i;
  int NH4_i;
  int NH4_pr_i;
  int DetPL_N_i;
  int Oxygen_i; 

  int Phy_PR_i;
  int Phy_Chl_i;

  int TN_i;
  int TP_i;
  int TC_i;
  int COD_i;

    /*
     * common cell variables
     */
  int Tfactor_i;
  int mL_i;
} workspace;

void dinoflagellate_spectral_mortality_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->mL_t0 = get_parameter_value(e, "PD_mL");
    ws->KO_aer = get_parameter_value(e, "KO_aer");
    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, "PhyD_N", e);
    ws->Phy_NR_i = e->find_index(tracers, "PhyD_NR", e);
    ws->Phy_I_i = e->find_index(tracers, "PhyD_I", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);

    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->COD_i = e->try_index(tracers, "COD", e);

    ws->Phy_PR_i = e->find_index(tracers, "PhyD_PR", e);
    ws->Phy_Chl_i = e->find_index(tracers, "PhyD_Chl", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, "PD_mL", e);
}

void dinoflagellate_spectral_mortality_wc_postinit(eprocess* p)
{
}

void dinoflagellate_spectral_mortality_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void dinoflagellate_spectral_mortality_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;

}

void dinoflagellate_spectral_mortality_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double Phy_N = y[ws->Phy_N_i];
    double Phy_NR = y[ws->Phy_NR_i];
    double Phy_PR = y[ws->Phy_PR_i];
    double Phy_I = y[ws->Phy_I_i];
    double Phy_Chl = y[ws->Phy_Chl_i];
    double Oxygen = e_max(y[ws->Oxygen_i]);

    double mortality1 = Phy_N * cv[ws->mL_i];
    double mortality2 = Phy_NR * cv[ws->mL_i];
    double mortality3 = Phy_I * cv[ws->mL_i];
    double mortality4 = Phy_PR * cv[ws->mL_i];
    
    /* Structural material (Phy_N) released at Redfield to DetPL, but all 
       reserves must go to individual pools (NR, PR), or lost (I, Chl)  */
    
    y1[ws->Phy_N_i] -= mortality1;
    y1[ws->Phy_NR_i] -= mortality2;
    y1[ws->Phy_I_i] -= mortality3;
    y1[ws->Phy_PR_i] -= mortality4;
    y1[ws->Phy_Chl_i] -= Phy_Chl * cv[ws->mL_i];
    
    y1[ws->DetPL_N_i] += mortality1 ;
    y1[ws->NH4_i] += mortality2 ;
    if (ws->NH4_pr_i > -1)
      y1[ws->NH4_pr_i] += mortality2 * SEC_PER_DAY ;
    y1[ws->DIP_i] += mortality4 ;
    y1[ws->DIC_i] += mortality3 * 106.0/1060.0*12.01 ;
    
    double Oxygen2 = Oxygen * Oxygen; 
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    y1[ws->Oxygen_i] -= mortality3 * 106.0/1060.0*12.01 * C_O_W * sigmoid ;
    
    if (ws->COD_i > -1){
      y1[ws->COD_i] += mortality3 * 106.0/1060.0*12.01 * C_O_W * (1.0-sigmoid);
    }

}

void dinoflagellate_spectral_mortality_wc_postcalc(eprocess* p, void* pp)
{
}
