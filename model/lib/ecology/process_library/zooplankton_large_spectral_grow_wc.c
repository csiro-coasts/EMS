/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/zooplankton_large_spectral_grow_wc.c
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
 *  $Id: zooplankton_large_spectral_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "zooplankton_large_spectral_grow_wc.h"

typedef struct {
  int do_mb;                  /* flag */
  int with_df;                /* flag */
  int with_mpb;               /* flag */
  int with_Tricho;            /* flag */
  
  /*
   * parameters
   */
  double umax_t0;
  double rad;
  char* meth;
  double swim_t0;
  double TKEeps;
  double MBrad;
  double PLrad;
  double DFrad;
  double Tricho_rad;
  double m;
  double E;
  double FDG;
  double KO_aer;
  double ZLdvmrate;
  double ZLpar;

  /*
   * tracers
   */
  int ZooL_N_i;
  int PhyL_N_i;
  int PhyL_NR_i;
  int PhyL_I_i;
  int PhyD_N_i;
  int PhyD_NR_i;
  int PhyD_I_i;
  int MPB_N_i;
  int MPB_NR_i;
  int MPB_I_i;
  int Tricho_N_i;
  int Tricho_NR_i;
  int Tricho_I_i;
  int NH4_i;
  int DetPL_N_i;
  int DIC_i;
  int Oxygen_i;
  int DIP_i;
  int temp_i;
  int salt_i;
  int ZooL_N_gr_i;
  int ZooL_N_rm_i;
  int NH4_pr_i;
  int Oxy_pr_i;
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;  
  int PhyL_PR_i;
  int PhyL_Chl_i;
  int MPB_PR_i;
  int MPB_Chl_i;
  int PhyD_PR_i;
  int PhyD_Chl_i;
  int Tricho_PR_i;
  int Tricho_Chl_i;
  int ZooL_sv_i;
  int PAR_i;


  /*
   * common cell variables
   */
  int Tfactor_i;
  int vis_i;
  int umax_i;
  int phi_MB_i;
  int phi_PL_i;
  int phi_DF_i;
  int phi_Tricho_i;
} workspace;

void zooplankton_large_spectral_grow_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*
   * tracers
   */
  ws->ZooL_N_i = e->find_index(tracers, "ZooL_N", e);
  ws->PhyL_N_i = e->find_index(tracers, "PhyL_N", e);
  ws->PhyL_NR_i = e->find_index(tracers, "PhyL_NR", e);
  ws->PhyL_I_i = e->find_index(tracers, "PhyL_I", e);
  ws->PhyL_PR_i = e->find_index(tracers, "PhyL_PR", e);
  ws->PhyL_Chl_i = e->find_index(tracers, "PhyL_Chl", e);

  // non essential 
  ws->PhyD_N_i = e->try_index(tracers, "PhyD_N", e);
  ws->PhyD_NR_i = e->try_index(tracers, "PhyD_NR", e);
  ws->PhyD_I_i = e->try_index(tracers, "PhyD_I", e);
  ws->PhyD_PR_i = e->try_index(tracers, "PhyD_PR", e);
  ws->PhyD_Chl_i = e->try_index(tracers, "PhyD_Chl", e);

  ws->MPB_N_i = e->try_index(tracers, "MPB_N", e);
  ws->MPB_NR_i = e->try_index(tracers, "MPB_NR", e);
  ws->MPB_I_i = e->try_index(tracers, "MPB_I", e);
  ws->MPB_PR_i = e->try_index(tracers, "MPB_PR", e);
  ws->MPB_Chl_i = e->try_index(tracers, "MPB_Chl", e);

  ws->Tricho_N_i = e->try_index(tracers, "Tricho_N", e);
  ws->Tricho_NR_i = e->try_index(tracers, "Tricho_NR", e);
  ws->Tricho_I_i = e->try_index(tracers, "Tricho_I", e);
  ws->Tricho_PR_i = e->try_index(tracers, "Tricho_PR", e);
  ws->Tricho_Chl_i = e->try_index(tracers, "Tricho_Chl", e);

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

  ws->ZooL_sv_i = e->try_index(tracers, "ZooL_sv", e);
  ws->PAR_i = e->try_index(tracers, "PAR", e);
  
  /*non essential diagnostic tracer */

  ws->ZooL_N_gr_i = e->try_index(tracers, "ZooL_N_gr", e);
  ws->ZooL_N_rm_i = e->try_index(tracers, "ZooL_N_rm", e);
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
  ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
  
  if (ws->PhyD_N_i > -1)
    ws->with_df = 1;
  else
    ws->with_df = 0;

  if (ws->Tricho_N_i > -1)
    ws->with_Tricho = 1;
  else
    ws->with_Tricho = 0;
  
  if (ws->MPB_N_i > -1)
    ws->with_mpb = 1;
  else
    ws->with_mpb = 0;

  ws->ZLdvmrate = try_parameter_value(e, "ZLdvmrate");
  if (isnan(ws->ZLdvmrate))
    ws->ZLdvmrate = 0.0;

  ws->ZLpar = try_parameter_value(e, "ZLpar");
  if (isnan(ws->ZLpar))
    ws->ZLpar = 1e-10;

  /*
   * parameters [KWA - note: if phyto_tracer exists Phyto_rad is essential for grazing to proceed!! ]
   */
  ws->umax_t0 = get_parameter_value(e, "ZLumax");
  ws->rad = get_parameter_value(e, "ZLrad");
  ws->meth = get_parameter_stringvalue(e, "ZLmeth");
  ws->swim_t0 = get_parameter_value(e, "ZLswim");
  ws->TKEeps = get_parameter_value(e, "TKEeps");
  ws->PLrad = get_parameter_value(e, "PLrad");
  
  if (ws->with_mpb)
    ws->MBrad = get_parameter_value(e, "MBrad");
  
  if (ws->with_df)
    ws->DFrad = get_parameter_value(e, "DFrad");

  if (ws->with_Tricho)
      ws->Tricho_rad = get_parameter_value(e, "Tricho_rad");

  ws->m = try_parameter_value(e, "ZLm");  /* mol P cell-1 */
  if (isnan(ws->m))
    ws->m = ZooCellMass(ws->rad);
  ws->E = get_parameter_value(e, "ZL_E");
  ws->FDG = get_parameter_value(e, "ZL_FDG");
  ws->KO_aer = get_parameter_value(e, "KO_aer");
  
  /*
   * common cell variables
   */
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->vis_i = find_index(e->cv_cell, "viscosity", e);
  ws->umax_i = find_index_or_add(e->cv_cell, "ZLumax", e);
  ws->phi_PL_i = find_index_or_add(e->cv_cell, "phi_PL", e);
  //if (ws->with_mpb)
  ws->phi_MB_i = find_index_or_add(e->cv_cell, "phi_MB", e);
  //if (ws->with_df)
  ws->phi_DF_i = find_index_or_add(e->cv_cell, "phi_DF", e);
  //if (ws->with_Tricho)
  ws->phi_Tricho_i = find_index_or_add(e->cv_cell, "phi_Tricho", e);
}

void zooplankton_large_spectral_grow_wc_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
  
  ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
  
  if(ws->with_df)
    ws->with_df = process_present(e,PT_WC,"dinoflagellate_spectral_grow_wc");
  emstag(LINFO,"eco:zooplankton_large_spectral_grow_wc:postinit","%sCalculating consumption of Dinoflagellates",(ws->with_df?"":"NOT "));
  
  if(ws->with_mpb)
    ws->with_mpb = process_present(e,PT_WC,"microphytobenthos_spectral_grow_wc");
  emstag(LINFO,"eco:zooplankton_large_spectral_grow_wc:postinit","%sCalculating consumption of MPB",(ws->with_mpb?"":"NOT "));
  
  if(ws->with_Tricho)
    ws->with_Tricho = process_present(e,PT_WC,"trichodesmium_grow_wc");
  emstag(LINFO,"eco:zooplankton_large_spectral_grow_wc:postinit","%sCalculating consumption of Tricho",(ws->with_Tricho?"":"NOT "));
}

void zooplankton_large_spectral_grow_wc_destroy(eprocess* p)
{
  free(p->workspace);
}

void zooplankton_large_spectral_grow_wc_precalc(eprocess* p, void* pp)
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
  
  cv[ws->phi_PL_i] = phi(ws->meth, ws->PLrad, 0.004 * pow(ws->PLrad, 0.26), 0.0, ws->rad, swim, 0.0, ws->TKEeps, vis, 1000.0, temp);
  cv[ws->phi_MB_i] = (ws->with_mpb) ? phi(ws->meth, ws->MBrad, 0.004 * pow(ws->MBrad, 0.26), 0.0, ws->rad, swim, 0.0, ws->TKEeps, vis, 1000.0, temp) : 0.0;
  cv[ws->phi_DF_i] = (ws->with_df) ? phi(ws->meth, ws->DFrad, 0.004 * pow(ws->DFrad, 0.26), 0.0, ws->rad, swim, 0.0, ws->TKEeps, vis, 1000.0, temp) : 0.0;
  cv[ws->phi_Tricho_i] = (ws->with_Tricho) ? phi(ws->meth, ws->Tricho_rad, 0.004 * pow(ws->Tricho_rad, 0.26), 0.0, ws->rad, swim, 0.0, ws->TKEeps, vis, 1000.0, temp) : 0.0;

  if (ws->do_mb) {
    double ZooL_N = y[ws->ZooL_N_i];
    
    y[ws->TN_i] += ZooL_N;
    y[ws->TP_i] += ZooL_N * red_W_P;
    y[ws->TC_i] += ZooL_N * red_W_C;
    y[ws->BOD_i] += ZooL_N * red_W_O;
  }
}

void zooplankton_large_spectral_grow_wc_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;

  /* Rearrange declarations to minimise if statements - 13/08/14 */

  double ZooL_N = y[ws->ZooL_N_i];

  double PhyL_N = y[ws->PhyL_N_i];
  double PhyL_NR = y[ws->PhyL_NR_i];
  double PhyL_I = y[ws->PhyL_I_i];
  double PhyL_PR = y[ws->PhyL_PR_i];
  double PhyL_Chl = y[ws->PhyL_Chl_i];
  double phi_PL = cv[ws->phi_PL_i];
  
  double PhyD_N,PhyD_NR,PhyD_I,PhyD_PR,PhyD_Chl,phi_DF;

  if (ws->with_df){
    PhyD_N = y[ws->PhyD_N_i];
    PhyD_NR = y[ws->PhyD_NR_i];
    PhyD_I = y[ws->PhyD_I_i];
    PhyD_PR = y[ws->PhyD_PR_i];
    PhyD_Chl = y[ws->PhyD_Chl_i];
    phi_DF = cv[ws->phi_DF_i];
  }else{
    PhyD_N = 0.0;
    PhyD_NR = 0.0;
    PhyD_I = 0.0;
    PhyD_PR = 0.0;
    PhyD_Chl = 0.0;
    phi_DF = 0.0;
  }

  double Tricho_N,Tricho_NR,Tricho_I,Tricho_PR,Tricho_Chl,phi_Tricho;

  if (ws->with_Tricho){
    Tricho_N = y[ws->Tricho_N_i];
    Tricho_NR = y[ws->Tricho_NR_i];
    Tricho_I = y[ws->Tricho_I_i];
    Tricho_PR = y[ws->Tricho_PR_i];
    Tricho_Chl = y[ws->Tricho_Chl_i];
    phi_Tricho = cv[ws->phi_Tricho_i];
  }else{
    Tricho_N = 0.0;
    Tricho_NR = 0.0;
    Tricho_I =  0.0;
    Tricho_PR = 0.0;
    Tricho_Chl = 0.0;
    phi_Tricho = 0.0;
  }

  double MPB_N,MPB_NR,MPB_I,MPB_PR,MPB_Chl,phi_MB;

  if (ws->with_mpb){
    MPB_N = y[ws->MPB_N_i];
    MPB_NR = y[ws->MPB_NR_i];
    MPB_I = y[ws->MPB_I_i];
    MPB_PR = y[ws->MPB_PR_i];
    MPB_Chl = y[ws->MPB_Chl_i];
    phi_MB = cv[ws->phi_MB_i];
  }else{
    MPB_N = 0.0;
    MPB_NR = 0.0;
    MPB_I = 0.0;
    MPB_PR = 0.0;
    MPB_Chl = 0.0;
    phi_MB = 0.0;
  };

  double Oxygen = y[ws->Oxygen_i];

  /* no zooplankton cells [cell m-3] = (mg N m-3) * (mol N mg N-1) / (mol P cell-1) / (mol N mol P-1) */ 

  double cells = ZooL_N * mgN2molN / ws->m / red_A_N;

  /* max_enc = per zooplankton cell number rate of encounter */

  double max_enc = PhyL_N * phi_PL + MPB_N * phi_MB + PhyD_N * phi_DF + Tricho_N * phi_Tricho;
  double umax = cv[ws->umax_i];

  /* Technical description is more obvious - max_ing is about taking umax and turning it into a clearance 
     rate per zooplankton cell so it can be compared with the max_enc */

  double max_ing = umax * ws->m * red_A_N / mgN2molN / ws->E;

  /* graze [mg N m-3 s-1] */

  double graze = cells * e_min(max_ing, max_enc);

  /*    if (graze == 0.0)
        max_enc = 1.0;
	
	double NH4_release = graze * (1.0 - ws->E) * (1.0 - ws->FDG);
	double growth = graze * ws->E;
	double DF_graze = (ws->with_df) ? graze * PhyD_N * phi_DF / max_enc : 0.0;
	double DIC_release = (ws->with_df) ? NH4_release * red_W_C + DF_graze * (PhyD_C / e_max(PhyD_N) - red_W_C) : NH4_release * red_W_C;
	double Oxy_pr = -DIC_release * red_W_O / red_W_C * Oxygen / (ws->KO_aer + e_max(Oxygen));
  */

  /*UR Changed order, strict compiler wants declarations first
   */
  double NH4_release ;
  double growth ;
  double DIC_release ;
  double Oxy_pr;
  
  if (graze == 0.0)
    max_enc = 1.0;
  
  NH4_release = graze * (1.0 - ws->E) * (1.0 - ws->FDG);
  growth = graze * ws->E;
  /*   DF_graze = (ws->with_df) ? graze * PhyD_N * phi_DF / max_enc : 0.0; */
  DIC_release = NH4_release * red_W_C ;
  
  /* need to do this later now that we have reserves */

  //Oxy_pr = -DIC_release * red_W_O / red_W_C * Oxygen / (ws->KO_aer + e_max(Oxygen));
  
  double grL = graze * phi_PL / max_enc ;
  double grMPB = (ws->with_mpb) ? graze * phi_MB / max_enc : 0.0;
  double grD = (ws->with_df) ? graze * phi_DF / max_enc : 0.0;
  double grTricho = (ws->with_Tricho) ? graze * phi_Tricho / max_enc : 0.0;

  /* all grazed nutrient reserve must go into dissolved pools as it is not at Redfield */
  
  y1[ws->ZooL_N_i] += growth;

  y1[ws->PhyL_N_i] -= grL * PhyL_N ;
  y1[ws->PhyL_NR_i] -= grL * PhyL_NR ;
  y1[ws->PhyL_I_i] -= grL * PhyL_I ;
  y1[ws->PhyL_PR_i] -= grL * PhyL_PR ;
  y1[ws->PhyL_Chl_i] -= grL * PhyL_Chl ;

  if (ws->with_df) {
    y1[ws->PhyD_N_i] -= grD * PhyD_N;
    y1[ws->PhyD_NR_i] -= grD * PhyD_NR;
    y1[ws->PhyD_I_i] -= grD * PhyD_I;
    y1[ws->PhyD_PR_i] -= grD * PhyD_PR ;
    y1[ws->PhyD_Chl_i] -= grD * PhyD_Chl ;
  }
  if (ws->with_Tricho) {
    y1[ws->Tricho_N_i] -= grTricho * Tricho_N;
    y1[ws->Tricho_NR_i] -= grTricho * Tricho_NR;
    y1[ws->Tricho_I_i] -= grTricho * Tricho_I;
    y1[ws->Tricho_PR_i] -= grTricho * Tricho_PR ;
    y1[ws->Tricho_Chl_i] -= grTricho * Tricho_Chl ;
  }
  if (ws->with_mpb) {
    y1[ws->MPB_N_i] -= grMPB * MPB_N ;
    y1[ws->MPB_NR_i] -= grMPB * MPB_NR ;
    y1[ws->MPB_I_i] -= grMPB * MPB_I ;
    y1[ws->MPB_PR_i] -= grMPB * MPB_PR ;
    y1[ws->MPB_Chl_i] -= grMPB * MPB_Chl ;
  }

  y1[ws->NH4_i] += NH4_release + (grL * PhyL_NR) + (grMPB * MPB_NR) + (grD * PhyD_NR) + (grTricho * Tricho_NR);
  y1[ws->DIP_i] += NH4_release * red_W_P + grL * PhyL_PR + grMPB * MPB_PR + grD * PhyD_PR + grTricho * Tricho_PR;

  DIC_release += (grL * PhyL_I + grMPB * MPB_I + grD * PhyD_I + grTricho * Tricho_I)* 106.0/1060.0 * 12.01;

  y1[ws->DIC_i] += DIC_release;

  Oxy_pr = -DIC_release * red_W_O / red_W_C * Oxygen / (ws->KO_aer + e_max(Oxygen));

  y1[ws->Oxygen_i] += Oxy_pr;

  y1[ws->DetPL_N_i] += graze * (1.0 - ws->E) * ws->FDG;
  
 if (ws-> NH4_pr_i > -1)
  y1[ws->NH4_pr_i] += (NH4_release + (grL * PhyL_NR) + (grMPB * MPB_NR) + (grD * PhyD_NR) + (grTricho * Tricho_NR)) * SEC_PER_DAY;
 
 if (ws->ZooL_N_rm_i  > -1)
   //  y1[ws->ZooL_N_rm_i] += graze * SEC_PER_DAY;
   y1[ws->ZooL_N_rm_i] += (grL * PhyL_N + grD * PhyD_N + grTricho * Tricho_N + grMPB * MPB_N) * red_W_C * SEC_PER_DAY;
 
  if (ws->ZooL_N_gr_i  > -1)
  y1[ws->ZooL_N_gr_i] += growth * SEC_PER_DAY;

  if (ws->Oxy_pr_i  > -1)
  y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY;


  /*UR 29/3/2005 changed to provide absolute growth rate
   * as advised by JP
   y1[ws->ZooL_N_gr_i] += growth / umax;
  */
}

void zooplankton_large_spectral_grow_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  double z_centre,z_bot;
  int wcbotk;

  double ZooL_N = y[ws->ZooL_N_i];
  
  y[ws->TN_i] += ZooL_N;
  y[ws->TP_i] += ZooL_N * red_W_P;
  y[ws->TC_i] += ZooL_N * red_W_C;
  y[ws->BOD_i] += ZooL_N * red_W_O;
  
  /* Vertical migration is light and total depth dependent. */

  /* only migrate if total depth greater than 100 m depth */

  if ((ws->ZooL_sv_i > -1) && (ws->PAR_i > -1)){

    /* Depths are -ve below mean sea level. */
    
    wcbotk = einterface_getwcbotk(c->col->model, c->b);
    z_bot = einterface_getcellz(c->col->model,c->b,wcbotk);
    
    if (z_bot > -100.0){
      y[ws->ZooL_sv_i] = 0.0;
      return;
    }
    
    z_centre = einterface_getcellz(c->col->model,c->b,c->k_wc);
    
    y[ws->ZooL_sv_i] = - ws->ZLdvmrate;
    
    if (y[ws->PAR_i]< ws->ZLpar){
      y[ws->ZooL_sv_i] = + ws->ZLdvmrate;
      if (z_centre > -10.0){  // allow diffusion to mix to surface.
	y[ws->ZooL_sv_i] = - 0.0;
      }
    }
  }
}
