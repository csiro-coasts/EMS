/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_spectral_grow_epi.c
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
 *  $Id: macroalgae_spectral_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "macroalgae_spectral_grow_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;                  /* flag */

  /*
   * parameters
   */
  double umax_t0;
  double MAleafden;

  /*
   * wc 
   */
  int MA_N_i;
  int MA_N_pr_i;
  int MA_N_gr_i;
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;
  int tau_i;


  int KI_MA_i;
  
  /*
   * tracers
   */
  
  int DIC_i;
  int NH4_i;
  int NO3_i;
  int DIP_i;
  int Oxygen_i;
  int Oxy_pr_i;

  /*
   * common cell variables
   */
  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int umax_i;

} workspace;

void macroalgae_spectral_grow_wc_init(eprocess* p)
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
    ws->umax_t0 = get_parameter_value(e, "MAumax_wc");
    ws->MAleafden = get_parameter_value(e, "MAleafden_wc");

    /*
     * e
     */
    ws->MA_N_i= e->find_index(tracers, "MA_N_wc", e);
    ws->TN_i = e->find_index(tracers, "TN", e) ;
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);


    /*
     * tracers
     */
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);

    /*non essential diagnostic tracers*/

    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);

    ws->MA_N_pr_i = e->try_index(tracers, "MA_N_pr_wc", e);

    ws->MA_N_gr_i = e->try_index(tracers, "MA_N_gr_wc", e);

      ws->tau_i = e->try_index(tracers, "tau_wc", e);




    /*
     * common variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "MAumax_wc", e);

}

void macroalgae_spectral_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance", e) >= 0) ? 1 : 0;
    ws->KI_MA_i = find_index_or_add(e->cv_cell, "KI_MA_wc", e);
}

void macroalgae_spectral_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void macroalgae_spectral_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
  ecology* e = p->ecology;

    void* model = e->model;



    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MA_N = y[ws->MA_N_i] * 1000.0 ;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;
    
    if (ws->do_mb) {
        y[ws->TN_i] += MA_N;
        y[ws->TP_i] += MA_N * atk_W_P;
        y[ws->TC_i] += MA_N * atk_W_C;
	y[ws->BOD_i] += MA_N * atk_W_O;
    }
}

void macroalgae_spectral_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double dz_wc = c->dz_wc;
  ecology* e = p->ecology;

    void* model = e->model;

    double MA_N = y[ws->MA_N_i];

    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Sc[2];
    
    /* Calculate maximum nutrient flux (Zhang 2011 Ecol. Mod. 222:1456-1470). */
    
    Sc[0] = 1.05e-6 / cv[ws->DNO3_i];
    Sc[1] = 1.05e-6 / cv[ws->DPO4_i];
    
    /* ws->tau in N m-2 */

    /* Density-specific shear stress from friction velocity. Density-specific 
       avoids multiplying by density before then dividing by it below */
    
    /*    double tau = 0.01;*/
    double tau = 0.01;/*/einterface_getustrcw(model, c->b);*/

    double S_DIN = 2850.0*pow((2.0*tau),0.38)*pow(Sc[0],-0.6)/86400.0; /* m s-1 */
    double S_DIP = 2850.0*pow((2.0*tau),0.38)*pow(Sc[1],-0.6)/86400.0; /* m s-1 */

    double umax = cv[ws->umax_i];  /* s-1 */
    
    double kI = c->cv[ws->KI_MA_i] * ( atk_A_N / atk_A_I );  /* mol eqN m-2 s-1 */
    
    /* available surface area */
    
    double SA = 1.0-exp(-MA_N * ws->MAleafden);   /* dimensionless */
    
    double kN_mass = S_DIN * SA * DIN * mgN2molN;  /* mol N m-2 s-1 */
    double kP_mass = S_DIP * SA * DIP * mgP2molP * (atk_A_N / atk_A_P ); /* mol eqN m-2 s-1 */
    
    /* MA_N now in g N m-2, so need to multiply  mgN2molN by 1000 */
    
    double growthrate = min(kI,min(kN_mass, kP_mass)) / (MA_N * mgN2molN * 1000.0); /* s-1 */
    
    growthrate = min(umax,growthrate);
    
    //  printf("umax %e, kI %e, kN %e, kP %e \n", umax,kI/(MA_N * mgN2molN * 1000.0),kN_mass/(MA_N * mgN2molN * 1000.0),kP_mass/(MA_N * mgN2molN * 1000.0));
    
    double growth = MA_N * growthrate;

    /* preferential ammonia uptake - watch units between growth and uptake */

    double k_NH4_mass = S_DIN * SA * NH4;    /* DNO3 only 4% different to DNH4 */
    double NH4uptake = min(k_NH4_mass,growth*1000.0);
    double NO3uptake = growth*1000.0 - NH4uptake;

    double Oxy_pr = 1000.0 * growth * atk_W_O ;
    
    // printf("%d, %d, %d, %d, %d, %d, %d, %d, %d \n",ws->MA_N_i,ws->MA_N_gr_i,ws->MA_N_pr_i,ws->NH4_wc_i,ws->NO3_wc_i,ws->DIP_wc_i,ws->DIC_wc_i,ws->Oxygen_wc_i,ws->Oxy_pr_wc_i);
    
    y1[ws->MA_N_i] += growth;
    // y1[ws->NH4_wc_i] -= growth * 1000.0 * NH4_wc / DIN_wc / dz_wc;
    // y1[ws->NO3_wc_i] -= growth * 1000.0 * NO3_wc / DIN_wc / dz_wc;
    y1[ws->NH4_i] -= NH4uptake ;
    y1[ws->NO3_i] -= NO3uptake ;
    y1[ws->DIP_i] -= growth * 1000.0 * atk_W_P ;
    y1[ws->DIC_i] -= growth * 1000.0 * atk_W_C ;
    y1[ws->Oxygen_i] += Oxy_pr;
    
    if (ws->MA_N_gr_i> -1)
      y1[ws->MA_N_gr_i] += growthrate / umax;
    if (ws->MA_N_pr_i> -1)
      y1[ws->MA_N_pr_i] += growth * SEC_PER_DAY;
    if (ws->Oxy_pr_i> -1)      
      y1[ws->Oxy_pr_i] += Oxy_pr * SEC_PER_DAY;     

    if (ws->tau_i> -1)      
      y1[ws->tau_i] +=tau;     



}
void macroalgae_spectral_grow_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double MA_N = y[ws->MA_N_i] * 1000.0;
  
  y[ws->TN_i] += MA_N;
  y[ws->TP_i] += MA_N * atk_W_P;
  y[ws->TC_i] += MA_N * atk_W_C;
  y[ws->BOD_i] += MA_N * atk_W_O;
}
