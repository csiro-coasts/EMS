/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/phytoplankton_spectral_grow_wc.c
 *  
 *  Description: Phytoplankton growth model.
 *
 *  Options: phytoplankton_spectral_grow_wc(small|large)
 *
 *  Small - State variables PhyS_*, parameters PS*
 *  Large - State variables PhyL_*, parameters PL*
 *
 *  Growth processes include: nutrient uptake into reserves, growth from reserves into structural material, 
 *  pigment synthesis, photosynthesis / respiration.
 *
 *  Model description: See zooxanthallae equations in:
 * 
 *  Baird, M. E., M. Mongin, F. Rizwi, L. K. Bay, N. E. Cantin, M. Soja-Wozniak and J. Skerratt (2018) 
 *  A mechanistic model of coral bleaching due to temperature-mediated light-driven reactive oxygen 
 *  build-up in zooxanthellae. Ecol. Model 386: 20-37.
 *
 *  Note: The published model description describes fixed carbon reserves (mg C cell-1). The code is written as 
 *         energy reverse (mmol photon cell-1) - the published description is clearer but the difference is just 
 *         a scaling factor. 
 *   
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: phytoplankton_spectral_grow_wc.c 7221 2022-09-25 23:06:34Z bai155 $
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
// #include "einterface.h"
#include "phytoplankton_spectral_grow_wc.h"

double ginterface_get_svel(void* model, char *name);

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;                  /* flag */
  /*
   * flag: 1 for large phytoplankton, 0 for small phytoplankton
   */
  int large;

  /*
   * parameters
   */
  double umax_t0;
  double psi;
  double m;
  double rad;
  double vol;
  double Chlmax;
  double Plank_resp;
  double C2Chlmin;

  /*
   * tracers
   */
  int Phy_N_i;
  int Phy_I_i;
  int Phy_NR_i;
  int Phy_PR_i;
  int Phy_N_pr_i;
  int Phy_N_gr_i;
  int NO3_i;
  int NH4_i;
  int DIP_i;
  int DIC_i;
  int Oxygen_i;
  int Oxy_pr_i;
  int Phy_Chl_i;
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;
  /*
   * common cell variables
   */

  int KI_i;
  int yCfac_i;

  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int umax_i;
} workspace;

void phytoplankton_spectral_grow_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  int nprms = p->prms->n;
  char* prm0 = p->prms->se[0]->s;
  char* prm1 = NULL;

  int large;

  p->workspace = ws;

  if (toupper(prm0[0]) == 'L')
    ws->large = 1;
  else if (toupper(prm0[0]) == 'S')
    ws->large = 0;
  else
    e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected first parameter: expected \"small\" or \"large\"\n", e->processfname, p->name, prm0);
  large = ws->large;

  /*
   * parameters
   */
  ws->umax_t0 = get_parameter_value(e, (large) ? "PLumax" : "PSumax");
  ws->m = PhyCellMass(get_parameter_value(e, (large) ? "PLrad" : "PSrad"));
  ws->rad = get_parameter_value(e, (large) ? "PLrad" : "PSrad");
  ws->psi = try_parameter_value(e, (large) ? "PLpsi" : "PSpsi");
  if (isnan(ws->psi)){
    ws->psi = psi(ws->rad);
    eco_write_setup(e,"Code default of psi (large = %d), %e \n",large,ws->psi);
  }
  ws->m = try_parameter_value(e, (large) ? "PLm" : "PSm");
  if (isnan(ws->m)){
    ws->m = PhyCellMass(ws->rad);
    eco_write_setup(e,"Code default of m (large = %d), %e \n",large,ws->m);
  }
  ws->Chlmax = PhyCellChl(get_parameter_value(e, (large) ? "PLrad" : "PSrad"));
  ws->Plank_resp = get_parameter_value(e, "Plank_resp");
  ws->vol = 4.0 / 3.0 * M_PI * ws->rad * ws->rad * ws->rad;
  ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
  if (isnan(ws->C2Chlmin))
    ws->C2Chlmin = 20.0;
  /*
   * tracers
   */
  ws->Phy_N_i = e->find_index(tracers, (large) ? "PhyL_N" : "PhyS_N", e);
  ws->Phy_I_i = e->find_index(tracers, (large) ? "PhyL_I" : "PhyS_I", e);
  ws->Phy_NR_i = e->find_index(tracers, (large) ? "PhyL_NR" : "PhyS_NR", e);
  ws->Phy_PR_i = e->find_index(tracers, (large) ? "PhyL_PR" : "PhyS_PR", e);
  ws->NO3_i = e->find_index(tracers, "NO3", e);
  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->DIP_i = e->find_index(tracers, "DIP", e);
  ws->DIC_i = e->find_index(tracers, "DIC", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);


  /*non essential diagnostic tracers */
  ws->Phy_N_pr_i = e->try_index(tracers, (large) ? "PhyL_N_pr" : "PhyS_N_pr", e);
  ws->Phy_N_gr_i = e->try_index(tracers, (large) ? "PhyL_N_gr" : "PhyS_N_gr", e);
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);

 /*common tracers */
  ws->Phy_Chl_i = e->find_index(tracers, (large) ? "PhyL_Chl" : "PhyS_Chl", e);
  ws->TN_i = e->find_index(tracers, "TN", e);
  ws->TP_i = e->find_index(tracers, "TP", e);
  ws->TC_i = e->find_index(tracers, "TC", e);
  ws->BOD_i = e->find_index(tracers, "BOD", e);

  /*
   * common cell variables
   */
  ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
  ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->umax_i = find_index_or_add(e->cv_cell, (large) ? "PLumax" : "PSumax", e);
}

void phytoplankton_spectral_grow_wc_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;

  ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

  ws->KI_i = -1;
  ws->yCfac_i = -1;

  /*
   * Define KI's here
   */
    // Must have the associated light model
    if (!process_present(e,PT_WC, "light_spectral_wc") && !process_present(e,PT_COL, "light_spectral_col"))
      emstag(LPANIC, "eco:phytoplankton_spectral_grow_wc_init",
	     "Phytoplankton grow is specified to be spectral but light_spectral_wc or light_spectral_col process not found!");
    if (ws->large){
      ws->KI_i = find_index_or_add(e->cv_cell, "KI_l", e);
      ws->yCfac_i = find_index_or_add(e->cv_cell, "yCfac_l", e);
    }else{
      ws->KI_i = find_index_or_add(e->cv_cell, "KI_s", e);
      ws->yCfac_i = find_index_or_add(e->cv_cell, "yCfac_s", e);
    }

 /*
  * Not valid during a pre_build (RECOM)
  */
    if (!e->pre_build) {

      /* test for equal sinking rates of structural material and reserves */
      
      double v1,v2,v3,v4,v5;
      
      if (ws->large){
	
	v1 = ginterface_get_svel(e->model,"PhyL_N");
	v2 = ginterface_get_svel(e->model,"PhyL_NR");
	v3 = ginterface_get_svel(e->model,"PhyL_PR");
	v4 = ginterface_get_svel(e->model,"PhyL_Chl");
	v5 = ginterface_get_svel(e->model,"PhyL_I");
	
	if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)){
	  e->quitfn("Sinking rates of PhyL :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
	}
      } else {
	
	v1 = ginterface_get_svel(e->model,"PhyS_N");
	v2 = ginterface_get_svel(e->model,"PhyS_NR");
	v3 = ginterface_get_svel(e->model,"PhyS_PR");
	v4 = ginterface_get_svel(e->model,"PhyS_Chl");
	v5 = ginterface_get_svel(e->model,"PhyS_I");
	
	if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)){
	  e->quitfn("Sinking rates of PhyS :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
	}
      }
    }
}

void phytoplankton_spectral_grow_wc_destroy(eprocess* p)
{
  workspace* ws = (workspace*) p->workspace;
  free(ws);
}

void phytoplankton_spectral_grow_wc_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  double* y = c->y;

  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
  double Phy_N = y[ws->Phy_N_i];
  double Phy_NR = y[ws->Phy_NR_i];
  double Phy_PR = y[ws->Phy_PR_i];
  double Phy_I = y[ws->Phy_I_i];
 
  cv[ws->umax_i] = ws->umax_t0 * Tfactor;

  if (ws->do_mb) {
    y[ws->TN_i] += Phy_N + Phy_NR;
    y[ws->TP_i] += Phy_N * red_W_P + Phy_PR;
    y[ws->TC_i] += Phy_N * red_W_C + Phy_I*106.0/1060.0*12.01;
    y[ws->BOD_i] += (Phy_N * red_W_C + Phy_I*106.0/1060.0*12.01)*C_O_W;
  }
}

void phytoplankton_spectral_grow_wc_calc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;

  double PI_max = ws->m * red_A_I * 1000.0;           /* mmol Q cell-1 */
  double PN_max = ws->m * red_A_N * 1000.0 * MW_Nitr; /* mg N cell-1 */
  double PP_max = ws->m * red_A_P * 1000.0 * MW_Phos; /* mg P cell-1 */
  double Phy_N = y[ws->Phy_N_i];                      /* mg N m-3 */
  double Phy_I = y[ws->Phy_I_i];                      /* mmol Q m-3 */
  double Phy_NR = y[ws->Phy_NR_i];
  double Phy_PR = y[ws->Phy_PR_i];
  double Phy_Chl = y[ws->Phy_Chl_i];
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

  if (Phy_N > 1e-9){  /* 1e-9 is to avoid unstable small divisions */

    double cellnum = Phy_N / PN_max ; /* cell m-3 */

  /* calculate N and I quotas, and cellular chl concentration from 
     total concentration and cell number */

    Nquota = max(0.0,min(1.0,(Phy_NR / cellnum) / PN_max));   /* dimensionless */
    Iquota = max(0.0,min(1.0,(Phy_I / cellnum) / PI_max));    /* dimensionless */
    Pquota = max(0.0,min(1.0,(Phy_PR / cellnum) / PP_max));    /* dimensionless */

    cellChl = Phy_Chl / (ws->vol * cellnum);    /* mg Chl m-3  */

    /* LIGHT ABSORPTION */

    double KI = c->cv[ws->KI_i];               // mol photon cell-1 s-1
    double Chlsynfactor = c->cv[ws->yCfac_i];  // d aA / d Chl a 
    
    double KN = (ws->psi * cv[ws->DNO3_i] * DIN);   /* mg N cell-1 s-1 */
    double KP = (ws->psi * cv[ws->DPO4_i] * DIP);   /* mg P cell-1 s-1 */

    Iuptake = KI*(1.0-Iquota) * cellnum * 1000.0; /* mmol photon m-3 s-1 */
    Nuptake = KN*(1.0-Nquota) * cellnum; /* mg N m-3 s-1 */
    Puptake = KP*(1.0-Pquota) * cellnum; /* mg P m-3 s-1 */
    
    /* preferential ammonia uptake */

    double KNH4 = (ws->psi * cv[ws->DNO3_i] * NH4); /* DNO3 only 4% different to DNH4 */
    double NH4uptake = min(KNH4 * cellnum, Nuptake);
    double NO3uptake = Nuptake - NH4uptake;

    /* PHYTOPLANKTON GROWTH */
    
    double growthrate = umax * Nquota * Iquota * Pquota; /* s-1 */
    
    double growth = Phy_N * growthrate;         /* mg N m-3 s-1 */
    
    /* growth is transfer of nutrient reserve NR to structual cell nutrient N
     * with corresponding reduction in energy reserve */
    
    /* PHYTOPLANKTON RESPIRATION mmol Q m-3 s-1
       basal respiration as a fraction of umax for fixed carbon   */
    
    double Iresp = umax * ws->Plank_resp * Phy_I;
    
    /* Chlorophyll content shared between daughter cells during division. 
       
       Chlorophyll synthesised based on:
       
       1. a maximum synthesis of the average chl in cells of that vol, 
          2.06e7*pow(1e18*vol,-0.320) in a cell doubling period, umax
       2. incremental benefit of adding Chl to absorption, Chlsynfactor
       3. the benefit this passes onto to growth, (1.0 - Iquota)             */

    double dChldt_dilute = growthrate * cellChl ;
    
    double tmmp = max(Chlsynfactor * (1.0 - Iquota),0.0);
    double dChldt_syn = ws->Chlmax * umax * min(tmmp, 1.33) * Nquota * Pquota;
   
    if (PN_max*5.6786 < ws->C2Chlmin * (Phy_Chl / cellnum)){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
    }
    
    /*	UPDATE RATES of STATE VARIABLES */
    
    y1[ws->Phy_I_i] += Iuptake - growth * red_A_I / (red_A_N * MW_Nitr) - Iresp ;
    y1[ws->Phy_NR_i] += Nuptake - growth;
    y1[ws->Phy_PR_i] += Puptake - growth * red_W_P;
    y1[ws->Phy_N_i] += growth;
    y1[ws->NO3_i] -= NO3uptake ;
    y1[ws->NH4_i] -= NH4uptake ;
    y1[ws->DIP_i] -= Puptake;
    y1[ws->DIC_i] += (Iresp - Iuptake) * 106.0/1060.0*12.01;
    y1[ws->Oxygen_i] += - (Iresp - Iuptake) * 106.0/1060.0*32.00 + NO3uptake * 48.0/14.01; 
    
    /* update water column chlorophyll concentration, units mg Chl m-3 s-1 */

    y1[ws->Phy_Chl_i] += (dChldt_syn - dChldt_dilute) * ws->vol * cellnum ;

    /* UPDATE DIAGNOSTICS */
     if (ws-> Phy_N_pr_i > -1)
       y1[ws->Phy_N_pr_i] += growth * SEC_PER_DAY * red_W_C;
     if (ws-> Phy_N_gr_i > -1)
       y1[ws->Phy_N_gr_i] = growthrate * SEC_PER_DAY;
     if (ws-> Oxy_pr_i > -1)
       y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
  }
}

void phytoplankton_spectral_grow_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

  double Phy_N = y[ws->Phy_N_i];
  double Phy_NR = y[ws->Phy_NR_i];
  double Phy_Chl = y[ws->Phy_Chl_i];
  double Phy_PR = y[ws->Phy_PR_i];
  double Phy_I = y[ws->Phy_I_i];

  y[ws->TN_i] += Phy_N + Phy_NR;
  y[ws->TP_i] += Phy_N * red_W_P + Phy_PR;
  y[ws->TC_i] += Phy_N * red_W_C + Phy_I*106.0/1060.0*12.01;
  y[ws->BOD_i] += (Phy_N * red_W_C + Phy_I*106.0/1060.0*12.01)*C_O_W;
}
