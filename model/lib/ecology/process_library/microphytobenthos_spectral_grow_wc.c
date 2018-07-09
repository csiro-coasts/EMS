/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/microphytobenthos_spectral_grow_wc.c
 *  
 *  Description: Microphytobenthos growth model.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: microphytobenthos_spectral_grow_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "microphytobenthos_spectral_grow_wc.h"
#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;                  /* flag */  

  /* parameters */
  
  double umax_t0;
  double psi;
  double m;
  double rad;
  double vol;
  double Chlmax; 
  double Plank_resp;
  double C2Chlmin;

  /*  tracers   */

  int MPB_N_i;
  int MPB_NR_i;
  int MPB_I_i;
  int MPB_PR_i;
  int NO3_i;
  int NH4_i;
  int DIP_i;
  int DIC_i;
  int Oxygen_i;
  int MPB_N_pr_i;
  int MPB_N_gr_i;
  int MPB_Chl_i;
  int Oxy_pr_i;
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;

  int KI_MPB_i;
  int yCfac_MPB_i;

  /* common cell variables */

  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int umax_i;
} workspace;

void microphytobenthos_spectral_grow_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*  parameters */
  
  ws->umax_t0 = get_parameter_value(e, "MBumax");

  ws->rad = try_parameter_value(e, "MBrad");
  if(isnan(ws->rad))
    e->quitfn("'MBrad' needed in microphytobenthos_spectral_grow_wc_init");

  ws->psi = try_parameter_value(e, "MBpsi");
  if (isnan(ws->psi))
    ws->psi = psi(ws->rad);

  ws->m = try_parameter_value(e, "MBm");

  if (isnan(ws->m))
    ws->m = PhyCellMass(ws->rad);

  ws->Plank_resp = get_parameter_value(e, "Plank_resp");

  ws->vol = 4.0 / 3.0 * M_PI * ws->rad * ws->rad * ws->rad;

  ws->Chlmax = PhyCellChl(ws->rad);

  ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
  if (isnan(ws->C2Chlmin))
    ws->C2Chlmin = 20.0;
  
  /* tracers */

  ws->MPB_N_i = e->find_index(tracers, "MPB_N", e);
  ws->MPB_NR_i = e->find_index(tracers, "MPB_NR", e);
  ws->MPB_PR_i = e->find_index(tracers, "MPB_PR", e);
  ws->MPB_I_i = e->find_index(tracers, "MPB_I", e);
  ws->NO3_i = e->find_index(tracers, "NO3", e);
  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->DIP_i = e->find_index(tracers, "DIP", e);
  ws->DIC_i = e->find_index(tracers, "DIC", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
  ws->MPB_Chl_i = e->find_index(tracers, "MPB_Chl", e);
  ws->TN_i = e->find_index(tracers, "TN", e);
  ws->TP_i = e->find_index(tracers, "TP", e);
  ws->TC_i = e->find_index(tracers, "TC", e);
  ws->BOD_i = e->find_index(tracers, "BOD", e);

  /* non-essential diagnostic tracers */
  
  ws->MPB_N_pr_i = e->try_index(tracers, "MPB_N_pr", e);
  ws->MPB_N_gr_i = e->try_index(tracers, "MPB_N_gr", e);
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);

  /*  common cell variables */

  ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
  ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->umax_i = find_index_or_add(e->cv_cell, "MBumax", e);

}

void microphytobenthos_spectral_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    double v1,v2,v3,v4,v5;

  /*
   * Not valid during a pre_build (RECOM)
   */
  if (!e->pre_build) {
    v1 = einterface_gettracersvel(e->model,"MPB_N");
    v2 = einterface_gettracersvel(e->model,"MPB_NR");
    v3 = einterface_gettracersvel(e->model,"MPB_PR");
    v4 = einterface_gettracersvel(e->model,"MPB_Chl");
    v5 = einterface_gettracersvel(e->model,"MPB_I");
    
    if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)){
      e->quitfn("Sinking rates of MPB :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
    }
  }
    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

    ws->KI_MPB_i = -1;
    ws->yCfac_MPB_i = -1;

  /*
   * Define KI's here
   */
 
    if (!process_present(e,PT_WC, "light_spectral_wc"))
      emstag(LPANIC, "eco:microphytobenthos_spectral_grow_wc_init",
	  "MPB grow is specified to be spectral but light_spectral_wc process not found!");
    ws->KI_MPB_i = find_index_or_add(e->cv_cell, "KI_MPB", e);
    ws->yCfac_MPB_i = find_index_or_add(e->cv_cell, "yCfac_MPB", e);

  /*
   * Not valid during a pre_build (RECOM)
   */
  if (!e->pre_build) {
    v1 = einterface_gettracersvel(e->model,"MPB_N");
    v2 = einterface_gettracersvel(e->model,"MPB_NR");
    v3 = einterface_gettracersvel(e->model,"MPB_PR");
    v4 = einterface_gettracersvel(e->model,"MPB_Chl");
    v5 = einterface_gettracersvel(e->model,"MPB_I");
    
    if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)){
      printf("Mass conservation violation due to microalgae reserves \n");
      printf("and structural material sinking at different rates \n");
      printf("Editing of .prm file required. \n");
      printf("Sinking of MPB :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
      exit(-1);
    } 
  }
}

void microphytobenthos_spectral_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void microphytobenthos_spectral_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MPB_N = y[ws->MPB_N_i];
    double MPB_NR = y[ws->MPB_NR_i];
    double MPB_PR = y[ws->MPB_PR_i];
    double MPB_I = y[ws->MPB_I_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    // y[ws->Kd_i] += ws->aA * MPB_N * mgN2molN / ws->m / red_A_N;

    if (ws->do_mb) {
        y[ws->TN_i] += MPB_N + MPB_NR;
        y[ws->TP_i] += MPB_N * red_W_P + MPB_PR;
        y[ws->TC_i] += MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01;
	y[ws->BOD_i] += (MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01)*C_O_W;
    }
}

void microphytobenthos_spectral_grow_wc_calc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double PI_max = ws->m * red_A_I * 1000.0 ; /* mmol Q . cell -1 */
  double PN_max = ws->m * red_A_N * 1000.0 * MW_Nitr; /* mg N cell-1 */
  double PP_max = ws->m * red_A_P * 1000.0 * MW_Phos; /* mg P cell-1 */
  double MPB_N = y[ws->MPB_N_i];
  double MPB_NR = y[ws->MPB_NR_i];
  double MPB_PR = y[ws->MPB_PR_i];
  double MPB_I = y[ws->MPB_I_i];
  double MPB_Chl = y[ws->MPB_Chl_i];
  double NO3 = y[ws->NO3_i];
  double NH4 = y[ws->NH4_i];
  double din = NO3 + NH4;
  double DIN = (din > 0.0) ? din : EPS_DIN;
  double dip = y[ws->DIP_i];
  double DIP = (dip > 0.0) ? dip : EPS_DIP;
  double Iuptake;
  double Nuptake;
  double Puptake;
  double umax = cv[ws->umax_i];	
  double Iquota;
  double Nquota;
  double Pquota;
  double cellChl;
  
/* obtain cell properties */

  if (MPB_N > 1e-9){  /* 1e-9 is to avoid unstable small divisions */

    //double cellnum = MPB_N / (ws->m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellnum = MPB_N / PN_max; /* cell m-3 */
  /* calculate N and I quotas, and cellular chl concentration from 
     total concentration and cell number */

    Nquota = max(0.0,min(1.0,(MPB_NR / cellnum) / PN_max));   /* dimensionless */
    Iquota = max(0.0,min(1.0,(MPB_I / cellnum) / PI_max));    /* dimensionless */
    Pquota = max(0.0,min(1.0,(MPB_PR / cellnum) / PP_max));    /* dimensionless */

    cellChl = MPB_Chl / (ws->vol * cellnum);    /* mg Chl m-3  */

    /* LIGHT ABSORPTION */

    double KI;            // mol photon cell-1 s-1
    double Chlsynfactor;  // d aA / d Chl a 

    KI = c->cv[ws->KI_MPB_i];
    Chlsynfactor = c->cv[ws->yCfac_MPB_i];

    double KN = (ws->psi * cv[ws->DNO3_i] * DIN );   /* mg N cell-1 s-1 */
    double KP = (ws->psi * cv[ws->DPO4_i] * DIP );   /* mg P cell-1 s-1 */

    Iuptake = KI*(1.0-Iquota) * cellnum * 1000.0; /* mmol photon m-3 s-1 */
    Nuptake = KN*(1.0-Nquota) * cellnum; /* mg N m-3 s-1 */
    Puptake = KP*(1.0-Pquota) * cellnum; /* mg P m-3 s-1 */

    /* preferential ammonia uptake */

    double KNH4 = (ws->psi * cv[ws->DNO3_i] * NH4); /* DNO3 only 4% different to DNH4 */
    double NH4uptake = min(KNH4 * cellnum, Nuptake);
    double NO3uptake = Nuptake - NH4uptake;

    /* PHYTOPLANKTON GROWTH */

    double growthrate = umax * Nquota * Iquota * Pquota; /* s-1 */

    double growth = MPB_N * growthrate;         /* mg N m-3 s-1 */

    /* growth is transfer of nutrient reserve NR to structual cell nutrient N
     * with corresponding reduction in energy reserve */
     
     /* PHYTOPLANKTON RESPIRATION mmol Q m-3 s-1
     basal respiration as a fraction of umax for fixed carbon, with sqrt(Iquota) 
     to stop respiration sending MPB_I negative.  */
    
    // double Iresp = umax * ws->Plank_resp * MPB_N * red_A_I / (red_A_N * MW_Nitr) * sqrt(Iquota);

    double Iresp = umax * ws->Plank_resp * MPB_I;

    /* Addition of variable Chl content equations 30 June 2011 MB 

       Chlorophyll content shared between daughter cells during division. 

       Chlorophyll synthesised based on:

       1. a maximum synthesis of the average chl in cells of that vol, 
          2.06e7*pow(1e18*vol,-0.320) in a cell doubling period, umax
       2. incremental benefit of adding Chl to absorption, Chlsynfactor
       3. the benefit this passes onto to growth, (1.0 - Iquota)             */

    double dChldt_dilute = growthrate * cellChl ;

    double tmmp = max(Chlsynfactor * (1.0 - Iquota),0.0);

    double dChldt_syn = ws->Chlmax * umax * min(tmmp, 1.33) * Nquota * Pquota;

    if (PN_max*5.6786 < ws->C2Chlmin * (MPB_Chl / cellnum)){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
    }

    /*	UPDATE RATES of STATE VARIABLES */

    y1[ws->MPB_I_i] += Iuptake - growth * red_A_I / (red_A_N * MW_Nitr) - Iresp ;
    y1[ws->MPB_NR_i] += Nuptake - growth;
    y1[ws->MPB_PR_i] += Puptake - growth * red_W_P;
    y1[ws->MPB_N_i] += growth;
    // y1[ws->NH4_i] -= Nuptake * NH4 / DIN;
    // y1[ws->NO3_i] -= Nuptake * NO3 / DIN;
    y1[ws->NO3_i] -= NO3uptake ;
    y1[ws->NH4_i] -= NH4uptake ;
    y1[ws->DIP_i] -= Puptake;
    // y1[ws->DIC_i] -= growth * red_W_C;
    // y1[ws->Oxygen_i] += growth * red_W_O;
    y1[ws->DIC_i] += (Iresp - Iuptake) * 106.0/1060.0*12.01;
    y1[ws->Oxygen_i] -= (Iresp - Iuptake) * 138.0/1060.0*32.00;

    /* update water column chlorophyll concentration, units mg Chl m-3 s-1 */

    y1[ws->MPB_Chl_i] += (dChldt_syn - dChldt_dilute) * ws->vol * cellnum ;

    /* UPDATE DIAGNOSTICS */

    if (ws->MPB_N_pr_i > -1)
      y1[ws->MPB_N_pr_i] += growth * SEC_PER_DAY * red_W_C;
    if (ws->MPB_N_gr_i > -1)
      y1[ws->MPB_N_gr_i] = growthrate * SEC_PER_DAY;
    if (ws->Oxy_pr_i > -1)
      y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
  }
}

void microphytobenthos_spectral_grow_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double MPB_N = y[ws->MPB_N_i];
  double MPB_NR = y[ws->MPB_NR_i];
  double MPB_PR = y[ws->MPB_PR_i];
  double MPB_I = y[ws->MPB_I_i];
  
  y[ws->TN_i] += MPB_N + MPB_NR;
  y[ws->TP_i] += MPB_N * red_W_P + MPB_PR;
  y[ws->TC_i] += MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01;  
  y[ws->BOD_i] += (MPB_N * red_W_C + MPB_I*106.0/1060.0*12.01)*C_O_W;

} 
