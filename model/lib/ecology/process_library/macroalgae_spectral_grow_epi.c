/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_spectral_grow_epi.c
 *  
 *  Description: Macroalgae growth.
 *  
 *  Mass transfer limited nutrient uptake of nitrogen (NH4 preferentially), phosphorus.
 *  Respiration considered in mortality.  
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: macroalgae_spectral_grow_epi.c 5946 2018-09-14 00:21:27Z bai155 $
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
#include "macroalgae_spectral_grow_epi.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;                  /* flag */

  /*
   * parameters
   */
  double umax_t0;
  double Benth_resp;
  double MAleafden;

  /*
   * epis
   */
  int MA_N_i;
  int MA_N_pr_i;
  int MA_N_gr_i;
  int EpiTN_i;
  int EpiTP_i;
  int EpiTC_i;
  int EpiBOD_i;
  
  int ustrcw_skin_i;

  int KI_MA_i;
  
  /*
   * tracers
   */
  
  int DIC_wc_i;
  int NH4_wc_i;
  int NO3_wc_i;
  int DIP_wc_i;
  int Oxygen_wc_i;
  int Oxy_pr_wc_i;

  /*
   * common cell variables
   */
  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int umax_i;

} workspace;

void macroalgae_spectral_grow_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "MAumax");
    ws->Benth_resp = get_parameter_value(e, "Benth_resp") * atk_A_I;
    ws->MAleafden = get_parameter_value(e, "MAleafden");

    /*
     * epis
     */
    ws->MA_N_i = e->find_index(epis, "MA_N", e) + OFFSET_EPI;
   

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;

    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*
     * tracers
     */
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);
    ws->DIP_wc_i = e->find_index(tracers, "DIP", e);


    /*non essential diagnostic tracers*/

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

    ws->MA_N_pr_i = e->try_index(epis, "MA_N_pr", e);
    if (ws->MA_N_pr_i > -1) 
	 ws->MA_N_pr_i  += OFFSET_EPI;

    ws->MA_N_gr_i = e->try_index(epis, "MA_N_gr", e);
    if (ws->MA_N_gr_i > -1) 
	 ws->MA_N_gr_i  += OFFSET_EPI;

    /*
     * common variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "MAumax", e);

}

void macroalgae_spectral_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;
    ws->KI_MA_i = find_index_or_add(e->cv_cell, "KI_MA", e);
}

void macroalgae_spectral_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void macroalgae_spectral_grow_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MA_N = y[ws->MA_N_i] * 1000.0 ;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;
    
    if (ws->do_mb) {
        y[ws->EpiTN_i] += MA_N;
        y[ws->EpiTP_i] += MA_N * atk_W_P;
        y[ws->EpiTC_i] += MA_N * atk_W_C;
	y[ws->EpiBOD_i] += MA_N * atk_W_O;
    }
}

void macroalgae_spectral_grow_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double dz_wc = c->dz_wc;

    double MA_N = y[ws->MA_N_i];

    if (MA_N < 1e-12)
      return;

    double NO3_wc = y[ws->NO3_wc_i];
    double NH4_wc = y[ws->NH4_wc_i];
    double din_wc = NO3_wc + NH4_wc;
    double DIN_wc = (din_wc > 0.0) ? din_wc : EPS_DIN;
    double dip_wc = y[ws->DIP_wc_i];
    double DIP_wc = (dip_wc > 0.0) ? dip_wc : EPS_DIP;
    double Sc[2];
    
    /* Calculate maximum nutrient flux (Zhang 2011 Ecol. Mod. 222:1456-1470). */
    
    Sc[0] = 1.05e-6 / cv[ws->DNO3_i];
    Sc[1] = 1.05e-6 / cv[ws->DPO4_i];
    
    /* ws->tau in N m-2 */

    /* Density-specific shear stress from friction velocity. Density-specific 
       avoids multiplying by density before then dividing by it below */
    
    double tau = y[ws->ustrcw_skin_i]*y[ws->ustrcw_skin_i];
    
    double S_DIN = 2850.0*pow((2.0*tau),0.38)*pow(Sc[0],-0.6)/86400.0; /* m s-1 */
    double S_DIP = 2850.0*pow((2.0*tau),0.38)*pow(Sc[1],-0.6)/86400.0; /* m s-1 */

    double umax = cv[ws->umax_i];  /* s-1 */
    
    double kI = c->cv[ws->KI_MA_i] * ( atk_A_N / atk_A_I );  /* mol eqN m-2 s-1 */

    /* available surface area */
    
    double SA = 1.0-exp(-MA_N * ws->MAleafden);   /* dimensionless */
    
    double kN_mass = S_DIN * SA * DIN_wc * mgN2molN;  /* mol N m-2 s-1 */
    double kP_mass = S_DIP * SA * DIP_wc * mgP2molP * (atk_A_N / atk_A_P ); /* mol eqN m-2 s-1 */
    
    /* MA_N now in g N m-2, so need to multiply  mgN2molN by 1000 */
    
    double growthrate = min(kI,min(kN_mass, kP_mass)) / (MA_N * mgN2molN * 1000.0); /* s-1 */
    
    growthrate = min(umax,growthrate);
    
    double growth = MA_N * growthrate;

    /* preferential ammonia uptake - watch units between growth and uptake */

    double k_NH4_mass = S_DIN * SA * NH4_wc;    /* DNO3 only 4% different to DNH4 */
    double NH4uptake = min(k_NH4_mass,growth*1000.0);
    double NO3uptake = growth*1000.0 - NH4uptake;

    double Oxy_pr = (1000.0 * growth * atk_W_O + NO3uptake * 48.0/14.01) / dz_wc;
    
    y1[ws->MA_N_i] += growth;
    y1[ws->NH4_wc_i] -= NH4uptake / dz_wc;
    y1[ws->NO3_wc_i] -= NO3uptake / dz_wc;
    y1[ws->DIP_wc_i] -= growth * 1000.0 * atk_W_P / dz_wc;
    y1[ws->DIC_wc_i] -= growth * 1000.0 * atk_W_C / dz_wc;
    y1[ws->Oxygen_wc_i] += Oxy_pr;
    
    if (ws->MA_N_gr_i> -1)
      y1[ws->MA_N_gr_i] += growthrate;
    if (ws->MA_N_pr_i> -1)
      y1[ws->MA_N_pr_i] += growth * SEC_PER_DAY;
    if (ws->Oxy_pr_wc_i> -1)      
      y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY;     

}
void macroalgae_spectral_grow_epi_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double MA_N = y[ws->MA_N_i] * 1000.0;
  
  y[ws->EpiTN_i] += MA_N;
  y[ws->EpiTP_i] += MA_N * atk_W_P;
  y[ws->EpiTC_i] += MA_N * atk_W_C;
  y[ws->EpiBOD_i] += MA_N * atk_W_O;
}
