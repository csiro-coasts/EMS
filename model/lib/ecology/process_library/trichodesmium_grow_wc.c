/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/trichodesmium_grow_wc.c
 *  
 *  Description: Trichodesmium growth model.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: trichodesmium_grow_wc.c 5908 2018-08-29 04:27:09Z bai155 $
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
#include "trichodesmium_grow_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;                  /* flag */    

  /* parameters */
  
  double umax_t0;
  double aA;
  double psi;
  double DINcrit;
  double N2;
  double Sh;
  double m;
  double rad;
  double p_max;
  double p_min;
  double vol;
  double Chlmax; 
  double Plank_resp;
  double C2Chlmin;

  /*  tracers   */

  int Tricho_N_i;
  int Tricho_NR_i;
  int Tricho_I_i;
  int Tricho_PR_i;
  int Tricho_sv_i;
  int NO3_i;
  int NH4_i;
  int DIP_i;
  int DIC_i;
  int Oxygen_i;
  int Tricho_N_pr_i;
  int Tricho_N_gr_i;
  int Nfix_i;
  int Tricho_Chl_i;
  int Oxy_pr_i;
  int dens_i;
  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;
  int KI_Tricho_i;
  int yCfac_Tricho_i;

  /* common cell variables */

  int DNO3_i;
  int DNH4_i;
  int DPO4_i;
  int DN2_i;
  int Tfactor_i;
  int vis_i;
  int umax_i;

} workspace;

void trichodesmium_grow_wc_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*  parameters */
  
  ws->p_min = get_parameter_value(e, "p_min");
  ws->p_max = get_parameter_value(e, "p_max");
  ws->umax_t0 = get_parameter_value(e, "Tricho_umax");
  ws->rad = get_parameter_value(e, "Tricho_rad");
  ws->N2 = get_parameter_value(e, "N2");
  ws->DINcrit = get_parameter_value(e, "DINcrit");
  ws->psi = try_parameter_value(e, "Tricho_psi");

  if (isnan(ws->psi))
    ws->psi = psi(ws->rad);

  ws->Sh = get_parameter_value(e, "Tricho_Sh");

  ws->m = try_parameter_value(e, "Tricho_m");

  if (isnan(ws->m))
    ws->m = PhyCellMass(ws->rad);

  ws->Plank_resp = get_parameter_value(e, "Plank_resp");

  ws->vol = 4.0 / 3.0 * M_PI * ws->rad * ws->rad * ws->rad;

  ws->Chlmax = PhyCellChl(ws->rad);

  ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
  if (isnan(ws->C2Chlmin))
    ws->C2Chlmin = 20.0;
  
  /* absorption done the old way until Kd removed from all code 
  ws->aA = try_parameter_value(e, "Tricho_aA");
  if (isnan(ws->aA)) {
    double absorb = get_parameter_value(e, "Tricho_absorb");
    ws->aA = aa(ws->rad, absorb);
  }
  */

  /* tracers */

  ws->Tricho_N_i = e->find_index(tracers, "Tricho_N", e);
  ws->Tricho_NR_i = e->find_index(tracers, "Tricho_NR", e);
  ws->Tricho_PR_i = e->find_index(tracers, "Tricho_PR", e);
  ws->Tricho_I_i = e->find_index(tracers, "Tricho_I", e);
  ws->Tricho_Chl_i = e->find_index(tracers, "Tricho_Chl", e);
  ws->NO3_i = e->find_index(tracers, "NO3", e);
  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->DIP_i = e->find_index(tracers, "DIP", e);
  ws->DIC_i = e->find_index(tracers, "DIC", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
  ws->Tricho_sv_i = e->find_index(tracers, "Tricho_sv", e);

  /* diagnostic tracers */

  ws->Tricho_N_pr_i = e->try_index(tracers, "Tricho_N_pr", e);
  ws->Tricho_N_gr_i = e->try_index(tracers, "Tricho_N_gr", e);
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
  ws->Nfix_i = e->find_index(tracers, "Nfix", e);
  
  /*
   * The density tracer, dens, will only be defined if CALCDENS YES is
   * specified in the prm file
   */
  ws->dens_i = e->find_index(e->tracers, "dens", e);
  ws->TN_i = e->find_index(tracers, "TN", e);
  ws->TP_i = e->find_index(tracers, "TP", e);
  ws->TC_i = e->find_index(tracers, "TC", e);
  ws->BOD_i = e->find_index(tracers, "BOD", e);

  /*  common cell variables */

  ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
  ws->DNH4_i = find_index(e->cv_cell, "DNH4", e);
  ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
  ws->DN2_i = find_index(e->cv_cell, "DN2", e);
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->vis_i = find_index(e->cv_cell, "viscosity", e);
  ws->umax_i = find_index_or_add(e->cv_cell, "Tricho_umax", e);
}

void trichodesmium_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

    ws->KI_Tricho_i = -1;
    ws->yCfac_Tricho_i = -1;

  /*
   * Define KI's here
   */
 
    if (!process_present(e,PT_WC, "light_spectral_wc"))
      emstag(LPANIC, "eco:trichodesmium_grow_wc_init",
	  "trichodesmium_grow_wc requires light_spectral_wc process!");
    ws->KI_Tricho_i = find_index_or_add(e->cv_cell, "KI_Tricho", e);
    ws->yCfac_Tricho_i = find_index_or_add(e->cv_cell, "yCfac_Tricho", e);

  /*
   * Not valid during a pre_build (RECOM)
   */
  if (!e->pre_build) {
    double v1,v2,v3,v4,v5;

    v1 = einterface_gettracersvel(e->model,"Tricho_N");
    v2 = einterface_gettracersvel(e->model,"Tricho_NR");
    v3 = einterface_gettracersvel(e->model,"Tricho_PR");
    v4 = einterface_gettracersvel(e->model,"Tricho_Chl");
    v5 = einterface_gettracersvel(e->model,"Tricho_I");
    
    if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)){
      e->quitfn("Sinking rates of Tricho :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
    } 
  }
}

void trichodesmium_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void trichodesmium_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double Tricho_N = y[ws->Tricho_N_i];
    double Tricho_NR = y[ws->Tricho_NR_i];
    double Tricho_PR = y[ws->Tricho_PR_i];
    double Tricho_I = y[ws->Tricho_I_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    if (ws->do_mb) {
        y[ws->TN_i] += Tricho_N + Tricho_NR;
        y[ws->TP_i] += Tricho_N * red_W_P + Tricho_PR;
        y[ws->TC_i] += Tricho_N * red_W_C + Tricho_I*106.0/1060.0*12.01;
	y[ws->BOD_i] += (Tricho_N * red_W_C + Tricho_I*106.0/1060.0*12.01)*C_O_W;
    }
}

void trichodesmium_grow_wc_calc(eprocess* p, void* pp)
{

  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double PI_max = ws->m * red_A_I * 1000.0 ; /* mmol Q . cell -1 */
  double PN_max = ws->m * red_A_N * 1000.0 * MW_Nitr; /* mg N cell-1 */
  double PP_max = ws->m * red_A_P * 1000.0 * MW_Phos; /* mg P cell-1 */
  double Tricho_N = y[ws->Tricho_N_i];
  double Tricho_NR = y[ws->Tricho_NR_i];
  double Tricho_PR = y[ws->Tricho_PR_i];
  double Tricho_I = y[ws->Tricho_I_i];
  double Tricho_Chl = y[ws->Tricho_Chl_i];
  double NO3 = y[ws->NO3_i];
  double NH4 = y[ws->NH4_i];
  double din = NO3 + NH4;
  double DIN = (din > 0.0) ? din : EPS_DIN;
  double dip = y[ws->DIP_i];
  double DIP = (dip > 0.0) ? dip : EPS_DIP;
  double Iuptake;
  //double NO3uptake;
  //double NH4uptake;
  //double N2uptake;
  double Puptake;
  double Nuptake;
  double umax = cv[ws->umax_i];	
  double Nquota;
  double Pquota;
  double Iquota;
  double cellChl;
  
/* obtain cell properties */

  if (Tricho_N > 1e-9){  /* to avoid unstable small divisions */

    // double cellnum = Tricho_N / (ws->m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellnum = Tricho_N / PN_max; /* cell m-3 */
  /* calculate N and I quotas, and cellular chl concentration from 
     total concentration and cell number */

    Nquota = max(0.0,min(1.0,(Tricho_NR / cellnum) / PN_max));   /* dimensionless */
    Iquota = max(0.0,min(1.0,(Tricho_I / cellnum) / PI_max));    /* dimensionless */
    Pquota = max(0.0,min(1.0,(Tricho_PR / cellnum) / PP_max));    /* dimensionless */

    // cellChl = min(PN_max*5.6786/12.0, Tricho_Chl / (ws->vol * cellnum));    /* mg Chl m-3  */

    cellChl = Tricho_Chl / (ws->vol * cellnum);

    /* LIGHT ABSORPTION */

    double KI;            // mol photon cell-1 s-1
    double Chlsynfactor;  // d aA / d Chl a 

    KI = c->cv[ws->KI_Tricho_i];
    Chlsynfactor = c->cv[ws->yCfac_Tricho_i];

    double KNH4 = (ws->psi * cv[ws->DNH4_i] * NH4 * ws->Sh);   /* mg N cell-1 s-1 */
    double KNO3 = (ws->psi * cv[ws->DNO3_i] * DIN * ws->Sh);   /* mg N cell-1 s-1 */
    double KP = (ws->psi * cv[ws->DPO4_i] * DIP * ws->Sh);   /* mg P cell-1 s-1 */
    double KN2 = (ws->psi * cv[ws->DNO3_i] * ws->DINcrit * ws->Sh);   /* mg N cell-1 s-1 */
    double Nitrogenase_prop = (DIN <= ws->DINcrit) ? 0.07 : 0.0;

    /* Uptake of dissolved nutrients and light/ATP */
    Iuptake = KI*(1.0-Iquota) * (1. - Nitrogenase_prop) * cellnum * 1000.0; /* mmol photon m-3 s-1 */
    Nuptake = KNO3 * (1.0-Nquota) * cellnum; /* mg N m-3 s-1 */
    Puptake = KP   * (1.0-Pquota) * cellnum; /* mg P m-3 s-1 */

    /* preferential ammonia uptake */

    double NH4uptake = min(KNH4 * cellnum, Nuptake);
    double NO3uptake = Nuptake - NH4uptake;
    //NH4uptake = KNH4*(1.0-Nquota) * cellnum; /* mg N m-3 s-1 */
    //NO3uptake = max(0.0,KNO3*(1.0-Nquota) * cellnum - NH4uptake); /* mg N m-3 s-1 */

    /* Nitrogen fixation */
    //    double N2uptake = max(0.0, (Iquota>0.0 & Pquota>0.0 & DIN <= ws->DINcrit) ? KN2 * (1.0-Nquota) * cellnum - NH4uptake - NO3uptake : 0.0);

    double N2uptake = 0.0;

    if (DIN <ws->DINcrit){
      N2uptake = max(0.0, KN2 * (1.0-Nquota) * Iquota * Pquota * cellnum - NH4uptake - NO3uptake);
    }
    
    /* TRICHODESMIUM GROWTH */

    double growthrate = umax * Nquota * Iquota * Pquota; /* s-1 */
    double growth = Tricho_N * growthrate;         /* mg N m-3 s-1 */

    /* Ifix is the energy loss associated with nitrogenase production */
    double Ifix = (DIN <= ws->DINcrit) ? Iuptake/3.0 : 0.0;

    /* growth is transfer of nutrient reserve NR to structual cell nutrient N
     * with corresponding reduction in energy reserve */
     
     /* TRICHODESMIUM RESPIRATION mmol Q m-3 s-1
     basal respiration as a fraction of umax for fixed carbon, with sqrt(Iquota) 
     to stop respiration sending Tricho_I negative.  */
    
    // double Iresp = umax * ws->Plank_resp * Tricho_N * red_A_I / (red_A_N * MW_Nitr) * sqrt(Iquota);

    /* Iresp is respiration in terms of energy usage. resp is used in conversion to O2 and DIC,
       and is analogous to growth */

    double Iresp = umax * ws->Plank_resp * Tricho_I;


    /* Chlorophyll content shared between daughter cells during division. 

       Chlorophyll synthesised based on:

       1. a maximum synthesis of the average chl in cells of that vol, 
          2.06e7*pow(1e18*vol,-0.320) in a cell doubling period, umax
       2. incremental benefit of adding Chl to absorption, Chlsynfactor
       3. the benefit this passes onto to growth, (1.0 - Iquota)             */

    double dChldt_dilute = growthrate * cellChl ;

    double tmmp = max(Chlsynfactor * (1.0 - Iquota),0.0);
    double dChldt_syn = ws->Chlmax * umax * min(tmmp, 1.33) * Nquota * Pquota;

    // double dChldt_syn = ws->Chlmax * umax * Chlsynfactor * (1.0 - ws->Iquota) * Nquota * Pquota;

    if (PN_max*5.6786 < ws->C2Chlmin * (Tricho_Chl / cellnum)){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
    }

    /*	UPDATE RATES of STATE VARIABLES */

    y1[ws->Tricho_I_i] += Iuptake - growth * red_A_I / (red_A_N * MW_Nitr) - Iresp - Ifix;
    y1[ws->Tricho_NR_i] += NH4uptake + NO3uptake + N2uptake - growth;
    y1[ws->Tricho_PR_i] += Puptake - growth * red_W_P;
    y1[ws->Tricho_N_i] += growth;

    y1[ws->NH4_i] -= NH4uptake;
    y1[ws->NO3_i] -= NO3uptake;
    y1[ws->DIP_i] -= Puptake;
    //    y1[ws->DIC_i] -= (growth - resp) * red_W_C;
    //    y1[ws->Oxygen_i] += (growth - resp) * red_W_O;
    y1[ws->DIC_i] += (Ifix + Iresp - Iuptake) * 106.0/1060.0*12.01;
    y1[ws->Oxygen_i] += - (Ifix + Iresp - Iuptake) * 106.0/1060.0*32.00 + NO3uptake * 48.0/14.01;

    /* update water column chlorophyll concentration, units mg Chl m-3 s-1 */

    y1[ws->Tricho_Chl_i] += (dChldt_syn - dChldt_dilute) * ws->vol * cellnum ;

    /* UPDATE DIAGNOSTICS */
    if (ws->Tricho_N_pr_i > -1)
      y1[ws->Tricho_N_pr_i] += growth * SEC_PER_DAY * red_W_C;
    if (ws->Tricho_N_gr_i > -1)
      y1[ws->Tricho_N_gr_i] = growthrate * SEC_PER_DAY;
    if (ws->Oxy_pr_i > -1)
      //   y1[ws->Oxy_pr_i] += (growth - resp) * red_W_O * SEC_PER_DAY;
      y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
    y1[ws->Nfix_i] += N2uptake;
  }
}

void trichodesmium_grow_wc_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;

  double* y = c->y;

  double g = 9.81;              /* Acceleration due to gravity */
//  double g = get_parameter_value(e, "g");
  double PI_max = ws->m * red_A_I * 1000.0 ; /* mmol Q . cell -1 */
  
  double Tricho_N = y[ws->Tricho_N_i];
  double Tricho_NR = y[ws->Tricho_NR_i];
  double Tricho_PR = y[ws->Tricho_PR_i];
  double Tricho_I = y[ws->Tricho_I_i];

  // double density = 1020.0; //y[ws->dens_i];
  double density = y[ws->dens_i];



  //  double Iquota = y[ws->Tricho_I_i]/PI_max;

  double cellnum = Tricho_N / (ws->m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */

  /* calculate N and I quotas, and cellular chl concentration from 
     total concentration and cell number */

  double Iquota = max(0.0,min(1.0,(Tricho_I / cellnum) / PI_max));  

  double Tricho_density = ws->p_min + Iquota * (ws->p_max - ws->p_min);

  y[ws->Tricho_sv_i] = -2.0/9.0 * (Tricho_density - density) * g * ws->rad * ws->rad / 1e-3; // dynamic_viscosity_w; 
//  changed units of sinking rate to m/s like all the other sinking rates that Nugzar's code handles 
//  also limited sinking to 3m/d (same as PhyL_N) as code was returning >1500m/d!  should still allow floating.

  // printf("Iquota %e, dens %e sink %e, dynamic_viscosity_w %e \n", Iquota,Tricho_density,y[ws->Tricho_sv_i],dynamic_viscosity_w);


// KWA additional non-conservative N2uptake accounted for in mass balance Nfix
  
  y[ws->TN_i] += Tricho_N + Tricho_NR;
  y[ws->TP_i] += Tricho_N * red_W_P + Tricho_PR;
  y[ws->TC_i] += Tricho_N * red_W_C + Tricho_I*106.0/1060.0*12.01;
  y[ws->BOD_i] += (Tricho_N * red_W_C + Tricho_I*106.0/1060.0*12.01)*C_O_W;
}
