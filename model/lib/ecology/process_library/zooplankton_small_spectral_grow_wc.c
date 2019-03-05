/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/zooplankton_small_spectral_grow_wc.c
 *  
 *  Description: Zooplankton growth from grazing on small zooplankton. 
 *
 *  Grazing equations include a maximum encounter rate and maximum growth rate limit.
 *
 *  Baird, M. E., S. J. Walker, B. B. Wallace, I. T. Webster and J. S. Parslow (2003) 
 *        The use of mechanistic descriptions of algal growth and zooplankton grazing 
 *        in an estuarine eutrophication model. 
 *        Estuarine, Coastal and Shelf Science 56, 685-695.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: zooplankton_small_spectral_grow_wc.c 6049 2018-12-20 03:02:48Z bai155 $
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
#include "zooplankton_small_spectral_grow_wc.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double umax_t0;
    double rad;
    char* meth;
    double swim_t0;
    double TKEeps;
    double PSrad;
    double m;
    double E;
    double FDG;
    double KO_aer;

    /*
     * tracers
     */
    int ZooS_N_i;
    int PhyS_N_i;
    int PhyS_NR_i;
    int PhyS_I_i;
    int NH4_i;
    int DetPL_N_i;
    int DIC_i;
    int Oxygen_i;
    int DIP_i;
    int temp_i;
    int salt_i;
    int ZooS_N_gr_i;
    int ZooS_N_rm_i;
    int NH4_pr_i;
    int Oxy_pr_i;
    int TN_i;
    int TP_i;
    int TC_i;
    int BOD_i;
    int COD_i;

    int PhyS_PR_i;
    int PhyS_Chl_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int vis_i;
    int umax_i;
    int phi_PS_i;
} workspace;

void zooplankton_small_spectral_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "ZSumax");
    ws->rad = get_parameter_value(e, "ZSrad");
    ws->meth = get_parameter_stringvalue(e, "ZSmeth");
    ws->swim_t0 = get_parameter_value(e, "ZSswim");
    ws->TKEeps = get_parameter_value(e, "TKEeps");
    ws->PSrad = get_parameter_value(e, "PSrad");
    ws->m = try_parameter_value(e, "ZSm");
    if (isnan(ws->m)){
      ws->m = ZooCellMass(ws->rad);
      eco_write_setup(e,"Code default of ZSm %e \n",ws->m);
    }
    ws->E = get_parameter_value(e, "ZS_E");
    ws->FDG = get_parameter_value(e, "ZS_FDG");
    ws->KO_aer = get_parameter_value(e, "KO_aer");

    /*
     * tracers
     */
    ws->ZooS_N_i = e->find_index(tracers, "ZooS_N", e);
    ws->PhyS_N_i = e->find_index(tracers, "PhyS_N", e);
    ws->PhyS_NR_i = e->find_index(tracers, "PhyS_NR", e);
    ws->PhyS_I_i = e->find_index(tracers, "PhyS_I", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
    ws->salt_i = e->find_index(tracers, "salt", e);
 
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);
    ws->COD_i = e->try_index(tracers, "COD", e);

    ws->PhyS_PR_i = e->find_index(tracers, "PhyS_PR", e);
    ws->PhyS_Chl_i = e->find_index(tracers, "PhyS_Chl", e);

    /*
     * common cell variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->vis_i = find_index(e->cv_cell, "viscosity", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "ZSumax", e);
    ws->phi_PS_i = find_index_or_add(e->cv_cell, "phi_PS", e);

    /*non essential diagnostic variables*/

   ws->ZooS_N_gr_i = e->try_index(tracers, "ZooS_N_gr", e);
    ws->ZooS_N_rm_i = e->try_index(tracers, "ZooS_N_rm", e);
    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);

}



void zooplankton_small_spectral_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
}

void zooplankton_small_spectral_grow_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void zooplankton_small_spectral_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* y = c->y;
    double* cv = c->cv;
    double temp = y[ws->temp_i];
    double vis = cv[ws->vis_i];
    double swim = ws->swim_t0;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;
    swim *= Tfactor;

    cv[ws->phi_PS_i] = phi(ws->meth, ws->PSrad, 0.004 * pow(ws->PSrad, 0.26), 0.0, ws->rad, swim, 0.0, ws->TKEeps, vis, 1000.0, temp);

    if (ws->do_mb) {
        double ZooS_N = y[ws->ZooS_N_i];

        y[ws->TN_i] += ZooS_N;
        y[ws->TP_i] += ZooS_N * red_W_P;
        y[ws->TC_i] += ZooS_N * red_W_C;
        y[ws->BOD_i] += ZooS_N * red_W_O;
    }
}

void zooplankton_small_spectral_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;
    double ZooS_N = y[ws->ZooS_N_i];
    double PhyS_N = y[ws->PhyS_N_i];
    double PhyS_NR = y[ws->PhyS_NR_i];
    double PhyS_I = y[ws->PhyS_I_i];
    double cells = ZooS_N * mgN2molN / ws->m / red_A_N;
    double phi_PS = cv[ws->phi_PS_i];
    double max_enc = PhyS_N * phi_PS;
    double umax = cv[ws->umax_i];
    double max_ing = umax * ws->m * red_A_N / mgN2molN / ws->E;
    double graze = cells * e_min(max_ing, max_enc);
    double NH4_release = graze * (1.0 - ws->E) * (1.0 - ws->FDG);
    double growth = graze * ws->E;
    double DIC_release = NH4_release * red_W_C;
    double Oxygen = e_max(y[ws->Oxygen_i]);

    double PhyS_PR = y[ws->PhyS_PR_i];
    double PhyS_Chl = y[ws->PhyS_Chl_i];

    if (graze == 0.0)
        max_enc = 1.0;

    double gr = graze * phi_PS / max_enc;

    /* all grazed nutrient reserve must go into dissolved pools as it is not at Redfield */

    y1[ws->ZooS_N_i] += growth;
    y1[ws->PhyS_N_i] -= gr * PhyS_N ;
    y1[ws->PhyS_NR_i] -= gr * PhyS_NR ;
    y1[ws->PhyS_I_i] -= gr * PhyS_I;
    y1[ws->NH4_i] += NH4_release + gr * PhyS_NR;

    y1[ws->DIP_i] += NH4_release * red_W_P + gr * PhyS_PR;
    y1[ws->PhyS_PR_i] -= gr * PhyS_PR ;
    y1[ws->PhyS_Chl_i] -= gr * PhyS_Chl;

    DIC_release += gr * PhyS_I * 106.0/1060.0 * 12.01;

    y1[ws->DIC_i] += DIC_release;

    double Oxygen2 = Oxygen * Oxygen; 
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    double Oxy_pr = -DIC_release * red_W_O / red_W_C * sigmoid;
    
    y1[ws->Oxygen_i] += Oxy_pr;

    if (ws->COD_i > -1)
      y1[ws->COD_i] += DIC_release * red_W_O / red_W_C * (1.0 - sigmoid);

    y1[ws->DetPL_N_i] += graze * (1.0 - ws->E) * ws->FDG;

    if (ws->ZooS_N_rm_i > -1)
      y1[ws->ZooS_N_rm_i] += gr * PhyS_N * red_W_C * SEC_PER_DAY;
    
    if (ws->NH4_pr_i > -1)
      y1[ws->NH4_pr_i] += (NH4_release + gr * PhyS_NR) * SEC_PER_DAY;
    
    if (ws->Oxy_pr_i > -1)
      y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY;
    
    if (ws->ZooS_N_gr_i  > -1)
      y1[ws->ZooS_N_gr_i] += growth * SEC_PER_DAY;
    
    
    /*UR 29/3/2005 changed to provide absolute growth rate
     * as advised by JP
     y1[ws->ZooS_N_gr_i] += growth / umax;
    */
}

void zooplankton_small_spectral_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double ZooS_N = y[ws->ZooS_N_i];

    y[ws->TN_i] += ZooS_N;
    y[ws->TP_i] += ZooS_N * red_W_P;
    y[ws->TC_i] += ZooS_N * red_W_C;
    y[ws->BOD_i] += ZooS_N * red_W_O;

}
