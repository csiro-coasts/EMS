/*
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/symbiodinium_spectral_free.c
 *  
 *  Description:
 * 
 *  This process contains:
 *    1. Growth and photophysiology of symbiodinium with physiology, xanthophyll cycle, reaction centre dynamics and 
 *       reative oxygen build-up.
 *
 *  To be considered: - Grazing loss in sediments.
 *                    - Photosynthesis in top sediment layer. 
 *
 *  To be used for free-living symbiodinium in both WC and SED.
 *
 *  References:
 * 
 *  Quigley KM, Bay LK and Willis BL (2017) Temperature and Water Quality-Related Patterns in
 *  Sediment-Associated Symbiodinium Communities Impact Symbiont Uptake and Fitness of Juveniles in the
 *  Genus Acropora. Front. Mar. Sci. 4:401. doi: 10.3389/fmars.2017.00401.
 *
 *  Baird, M. E., M. Mongin, F. Rizwi, L. K. Bay, N. E. Cantin, M. Soja-Wozniak and J. Skerratt (2018) 
 *  A mechanistic model of coral bleaching due to temperature-mediated light-driven reactive oxygen 
 *  build-up in zooxanthellae. Ecol. Model 386: 20-37.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
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
#include "symbiodinium_spectral_free.h"

double ginterface_get_svel(void* model, char *name);

#define unitch 1000.0

typedef struct {
  int do_mb;                  /* flag */

  /*
   * parameters
   */
  double CSumax_t0;
  double CSm;
  double CSrad;
  double CSvol;
  double Chlmax; 
  double Plank_resp;
  double CSmort_t0;
  double C2Chlmin;
  double KO_aer;

  double psi;

  /* new parameters for xanthophyll / bleaching model */

  double Xanth_tau;
  double chla2rcii;
  double photon2rcii;
  double xan2chl_CS;
  double ROSmult;
  double RubiscoOffTemp;
  double photon2ros;

  /*
   * tracers
   */

  int CS_N_i;    /* Coral symbiont - microalgae */
  int CS_Chl_i;

  int CS_NR_i;
  int CS_PR_i;
  int CS_I_i;

  int CS_Xh_i;
  int CS_Xp_i;

  int CS_Qred_i;
  int CS_Qox_i;
  int CS_Qi_i;
  int CS_RO_i;

  int DetPL_N_i;

  int CS_tempfunc_i;  
  int temp_clim_wc_i;

  int TN_i;
  int TP_i;
  int TC_i;
  int BOD_i;

  /*
   * tracers
   */
  
  int temp_wc_i;
  int PAR_i;
  int DIC_i;
  int NH4_i;
  int NO3_i;
  int DIP_i;
  int Oxygen_i;
  int COD_i;
  int Oxy_pr_wc_i;

  int CS_N_pr_i;
  int CS_bleach_i;

  int KI_CS_free_i;
  int yCfac_CS_free_i;

  int do_light_spectral_col;
  
  /*
   * common cell variables
   */
  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int CSumax_i;
  int CSmort_i;

} workspace;

void symbiodinium_spectral_free_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;

    p->workspace = ws;

    /*
     * parameters
     */

    ws->CSumax_t0 = try_parameter_value(e, "CS_free_umax");
    if (isnan(ws->CSumax_t0)){
      ws->CSumax_t0 = 1.0/86400.0;
      eco_write_setup(e,"Code default of CS_free_umax = %e \n",ws->CSumax_t0);
    }

    ws->CSrad = try_parameter_value(e, "CSrad");
    if (isnan(ws->CSrad)){
      ws->CSrad = 5.0e-6;
      eco_write_setup(e,"Code default of CS_red = %e \n",ws->CSrad);
    }

    ws->Chlmax = PhyCellChl(ws->CSrad);
    ws->CSvol = (4.0*M_PI/3.0) * ws->CSrad * ws->CSrad * ws->CSrad;
    ws->CSm = PhyCellMass(ws->CSrad);

    eco_write_setup(e,"Calculated value of CSm = %e \n",ws->CSm);

    ws->psi = try_parameter_value(e, "CSpsi");
    if (isnan(ws->psi)){
      ws->psi = psi(ws->CSrad);
      eco_write_setup(e,"Code default of CSpsi = %e \n",ws->psi);
    }

    ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
    if (isnan(ws->C2Chlmin)){
      ws->C2Chlmin = 20.0;
      eco_write_setup(e,"Code default of C2Chlmin = %e \n",ws->C2Chlmin);
    }
    ws->CSmort_t0 = try_parameter_value(e, "CSmort");
    if (isnan(ws->CSmort_t0)){
      ws->CSmort_t0 = 0.1/86400.0;
      eco_write_setup(e,"Code default of CS_free_mort = %e \n",ws->CSmort_t0);
    }

    ws->Xanth_tau = try_parameter_value(e, "Xanth_tau");
    if (isnan(ws->Xanth_tau)){
      ws->Xanth_tau = 1.0/1200.0;   // 20 mins: Gustafsson et al., 2014
      eco_write_setup(e,"Code default of Xanth_tau = %e \n",ws->Xanth_tau);
    }

    ws->chla2rcii = try_parameter_value(e, "chla2rcii");
    if (isnan(ws->chla2rcii)){
      ws->chla2rcii = 0.002/893.49; // mol RCII g Chl-1 [Suggett et al., 2009], MW Chla 893.49.
      eco_write_setup(e,"Code default of chla2rcii = %e \n",ws->chla2rcii);
    }

    ws->photon2rcii = try_parameter_value(e, "photon2rcii");
    if (isnan(ws->photon2rcii)){
      ws->photon2rcii = 0.1e-6; // 3.93e-8; // mol RCII mol photon-1
      eco_write_setup(e,"Code default of photon2rcii = %e \n",ws->photon2rcii);
    }

    ws->ROSmult = try_parameter_value(e, "ROSmult");
    if (isnan(ws->ROSmult)){
      ws->ROSmult = 1.0; 
        // empirical fit to field observations
      eco_write_setup(e,"Code default of ROSmult = %e \n",ws->ROSmult);
    }

    ws->photon2ros = try_parameter_value(e, "CS_photon2ros");
    if (isnan(ws->photon2ros)){
      ws->photon2ros = 7000.0;  // empirical fit to field observations
      eco_write_setup(e,"Code default of photon2ros = %e \n",ws->photon2ros);
    }

    // Calculated assuming absorption cross section is pi r^2 m cell-1, and Suggett 2008 measured absorption cross section of PSII per photon.

    ws->xan2chl_CS = try_parameter_value(e, "xan2chl_CS");
    if (isnan(ws->xan2chl_CS)){
      ws->xan2chl_CS = 0.2448;   // from phytoplankton ???
      eco_write_setup(e,"Code default of xan2chl_CS = %e \n",ws->xan2chl_CS);
    }

    ws->Plank_resp = try_parameter_value(e, "Plank_resp");
    if (isnan(ws->Plank_resp)){
      ws->Plank_resp = 0.1;   // fraction of maximum growth rate.
      eco_write_setup(e,"Code default of Plank_resp = %e \n",ws->Plank_resp);
    }

    ws->RubiscoOffTemp = try_parameter_value(e, "RubiscoOffTemp");
    if (isnan(ws->RubiscoOffTemp)){
      ws->RubiscoOffTemp = 2.0;  // temperature above climatology at which bleaching rate maximum 
      eco_write_setup(e,"Code default of RubiscoOffTemp = %e \n",ws->RubiscoOffTemp);
    }

    ws->KO_aer = get_parameter_value(e, "KO_aer");

    /*
     * tracers
     */
 
    ws->CS_N_i = e->find_index(tracers, "CS_N_free", e);
    ws->CS_Chl_i = e->find_index(tracers, "CS_Chl_free", e);

    ws->CS_NR_i = e->find_index(tracers, "CS_NR_free", e);
    ws->CS_PR_i = e->find_index(tracers, "CS_PR_free", e);
    ws->CS_I_i = e->find_index(tracers, "CS_I_free", e);

    ws->CS_Xh_i = e->find_index(tracers, "CS_Xh_free", e);
    ws->CS_Xp_i = e->find_index(tracers, "CS_Xp_free", e);

    ws->CS_Qred_i = e->find_index(tracers, "CS_Qred_free", e);
    ws->CS_Qox_i = e->find_index(tracers, "CS_Qox_free", e);
    ws->CS_Qi_i = e->find_index(tracers, "CS_Qi_free", e);
    ws->CS_RO_i = e->find_index(tracers, "CS_RO_free", e);

    ws->TN_i = e->find_index(tracers, "TN", e); 
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);

    ws->temp_wc_i = e->find_index(tracers, "temp", e);
    ws->temp_clim_wc_i = e->try_index(tracers, "temp_clim", e);
    ws->PAR_i = e->find_index(tracers, "PAR", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->COD_i = e->find_index(tracers, "COD", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);

    /*
     * common variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);

    ws->CSumax_i = find_index_or_add(e->cv_cell, "CSumax_free", e);
    ws->CSmort_i = find_index_or_add(e->cv_cell, "CSmort_free", e);
    ws->CS_tempfunc_i = find_index_or_add(e->cv_cell, "CStempfunc", e);
   
    /*non essential diagnositc tracer*/

    ws->CS_N_pr_i = e->try_index(tracers, "CS_N_free_pr", e);
    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

}

void symbiodinium_spectral_free_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11;

    ws->KI_CS_free_i = find_index_or_add(e->cv_cell, "KI_CS_free", e);
    ws->yCfac_CS_free_i = find_index_or_add(e->cv_cell, "yCfac_CS_free", e);

    /*
     * Not valid during a pre_build (RECOM)
     */
    if (!e->pre_build) {
      v1 = ginterface_get_svel(e->model,"CS_N_free");
      v2 = ginterface_get_svel(e->model,"CS_NR_free");
      v3 = ginterface_get_svel(e->model,"CS_PR_free");
      v4 = ginterface_get_svel(e->model,"CS_Chl_free");
      v5 = ginterface_get_svel(e->model,"CS_I_free");
      v6 = ginterface_get_svel(e->model,"CS_Xh_free");
      v7 = ginterface_get_svel(e->model,"CS_Xp_free");
      v8 = ginterface_get_svel(e->model,"CS_RO_free");
      v9 = ginterface_get_svel(e->model,"CS_Qox_free");
      v10 = ginterface_get_svel(e->model,"CS_Qred_free");
      v11 = ginterface_get_svel(e->model,"CS_Qi_free");
      
      if ((v1!=v2)||(v1!=v3)||(v1!=v4)||(v1!=v5)||(v1!=v6)||(v1!=v7)||(v1!=v8)||(v1!=v9)||(v1!=v10)||(v1!=v11)){
	printf("Mass conservation violation due to Symbiodinium cell contents \n");
	printf("and structural material sinking at different rates \n");
	printf("Editing of .prm file required. \n");
	printf("Sinking of Symbiodinium :N %e, NR %e, PR %e, Chl %e, I %e are not equal",v1,v2,v3,v4,v5);
	printf("Sinking of Symbiodinium :Xh %e, Xp %e, RO %e, Qox %e, Qred %e, Qi %e are not equal",v6,v7,v8,v9,v10,v11);
	exit(-1);
      } 
    }
    ws->do_light_spectral_col = 0;
    if (process_present(e,PT_COL,"light_spectral_col")){
      ws->do_light_spectral_col = 1;
    }
}

void symbiodinium_spectral_free_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;
    free(ws);
}

void symbiodinium_spectral_free_precalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Coral_N =  y[ws->CS_N_i];

    /* Coupling with physical processes (temperature, shear stress, light) is calculated only 
       at the beginning of the ecological time increment */

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->CSumax_i] = ws->CSumax_t0 * Tfactor;
    cv[ws->CSmort_i] = ws->CSmort_t0 * Tfactor;
    cv[ws->CS_tempfunc_i] = 1.0;

    if (ws->temp_clim_wc_i > -1){
      if (y[ws->temp_clim_wc_i] > 10.0){  // sometimes zeros come thru - if so, leave as last time step.
	double deltemp = min(max(0.0,(y[ws->temp_wc_i] - max(y[ws->temp_clim_wc_i],26.0))),ws->RubiscoOffTemp);
	cv[ws->CS_tempfunc_i] = (1.0 - exp(-(ws->RubiscoOffTemp-deltemp)))/(1.0-exp(-(ws->RubiscoOffTemp)));
      }
    }

    /* Mass balance must include all state variables introduced in this routine that hold N */

    y[ws->TN_i] += Coral_N + y[ws->CS_NR_i];
    y[ws->TP_i] += Coral_N * red_W_P + y[ws->CS_PR_i];
    y[ws->TC_i] += Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01;
    y[ws->BOD_i] += (Coral_N * red_W_C + y[ws->CS_I_i] *106.0/1060.0*12.01)*C_O_W;

}

void symbiodinium_spectral_free_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;

    double CS_N = y[ws->CS_N_i];

    if (CS_N < 1e-9){
      return;
    }

    double porosity = c->porosity;
    double dz = 0.0; // layer thickness x porosity.

    double kI, Chlsynfactor;

    double* cv = c->cv;

    double CSumax = cv[ws->CSumax_i];  /* s-1 */
    double CSmort = cv[ws->CSmort_i];  /* s-1 */
  
    if (ws->do_light_spectral_col) {
      if (p->type == PT_WC){
	dz = c->dz_wc;
	kI = c->cv[ws->KI_CS_free_i]; // mol photon cell-1 s-1 */
	Chlsynfactor = max(0.0,min(c->cv[ws->yCfac_CS_free_i],1.33)); /* used below for chl synthesis calculation */
      }else{ 
	dz = c->dz_sed * porosity;
	kI = c->cv[ws->KI_CS_free_i]; // mol photon cell-1 s-1 */
	Chlsynfactor = max(0.0,min(c->cv[ws->yCfac_CS_free_i],1.33)); /* used below for chl synthesis calculation */
      }
    }else{
      if (p->type == PT_WC){
	dz = c->dz_wc;
	kI = 0.0; // mol photon cell-1 s-1 */
	Chlsynfactor = 0.0; /* used below for chl synthesis calculation */
      }else{ 
	dz = c->dz_sed * porosity;
	kI = 0.0; // mol photon cell-1 s-1 */
	Chlsynfactor = 0.0; /* used below for chl synthesis calculation */	
      }
    }

    double CS_Chl = y[ws->CS_Chl_i];

    double CS_NR = y[ws->CS_NR_i];
    double CS_PR = y[ws->CS_PR_i];
    double CS_I = y[ws->CS_I_i];
    double CS_Xp = y[ws->CS_Xp_i];
    double CS_Xh = y[ws->CS_Xh_i];

    double CS_Qi = y[ws->CS_Qi_i];
    double CS_Qred = y[ws->CS_Qred_i];
    double CS_Qox = y[ws->CS_Qox_i];
    double CS_RO = y[ws->CS_RO_i];

    double CS_Qt = CS_Qi + CS_Qred + CS_Qox; /* total reaction centres.*/
    double CS_Qa = CS_Qred + CS_Qox; /* active reaction centres.*/

    double PI_max = ws->CSm * red_A_I * 1000.0 ; /* mmol photon cell -1 */
    double PN_max = ws->CSm * red_A_N * 1000.0 * MW_Nitr; /* mg N cell-1 */
    double PP_max = ws->CSm * red_A_P * 1000.0 * MW_Phos; /* mg P cell-1 */

    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double DIN = NO3 + NH4;
    double DIP = y[ws->DIP_i];
    double DIC =  y[ws->DIC_i];

     /* Zoothanthellae - cells quantified per m2 of host tissue. */

    double cellnum = CS_N / (ws->CSm * red_A_N * 1000.0 * MW_Nitr); /* cell m-2 */
    double cellChl = CS_Chl / (ws->CSvol * cellnum);                /* mg Chl m-3  */
    double cellXp = CS_Xp / (ws->CSvol * cellnum);                /* mg Chl m-3  */
    double cellXh = CS_Xh / (ws->CSvol * cellnum);                /* mg Chl m-3  */

    /* now get normalised reserves. */

    double Nquota = max(0.0,min(1.0,(CS_NR / cellnum) / PN_max));   /* dimensionless */
    double Iquota = max(0.0,min(1.0,(CS_I / cellnum) / PI_max));    /* dimensionless */
    double Pquota = max(0.0,min(1.0,(CS_PR / cellnum) / PP_max));    /* dimensionless */

    /* Information from light_spectral_epi.c */

    double KN = (ws->psi * cv[ws->DNO3_i] * DIN);   /* mg N cell-1 s-1 */
    double KP = (ws->psi * cv[ws->DPO4_i] * DIP);   /* mg P cell-1 s-1 */

    double Nuptake = KN*(1.0-Nquota) * cellnum; /* mg N m-3 s-1 */
    double Puptake = KP*(1.0-Pquota) * cellnum; /* mg P m-3 s-1 */
    
    /* preferential ammonia uptake */

    double KNH4 = (ws->psi * cv[ws->DNO3_i] * NH4); /* DNO3 only 4% different to DNH4 */
    double NH4uptake = min(KNH4 * cellnum, Nuptake);
    double NO3uptake = Nuptake - NH4uptake;

    double Iuptake = 0.0;
    
    if (CS_Qt > 1.0e-9){
     Iuptake = kI * CS_Qox/CS_Qt * cv[ws->CS_tempfunc_i] * (1.0 - Iquota) * cellnum * 1000.0;
    }

    /* respiration based on internal reserves */

    double Iresp = CSumax * ws->Plank_resp * CS_I;  /* mmol photon m-2 s-1 */

    double symb_grow = CSumax*Nquota*Iquota*Pquota; /* s-1 */

    /* Death of polyps reduces CS_N */

    y1[ws->CS_N_i] += (CS_N * symb_grow);

    /* need to add lost reserves to the water column, including that from translocated cells and those lost to mucus */

    y1[ws->NH4_i] += - NH4uptake / porosity;
    y1[ws->NO3_i] += - NO3uptake / porosity;
    y1[ws->DIP_i] += - Puptake / porosity;

    /* respiration, carbon fixation and nitrate uptake */

    y1[ws->DIC_i] += (Iresp - Iuptake) * 106.0/1060.0*12.01 / porosity;
    y1[ws->Oxygen_i] += (-(Iresp - Iuptake) * 106.0/1060.0*32.00 + NO3uptake * 48.0/14.01) / porosity;

    /* Uptake from water column - growth of cells  */

    y1[ws->CS_I_i] += Iuptake - symb_grow * red_A_I / (red_A_N * MW_Nitr) * CS_N - Iresp ;
    y1[ws->CS_NR_i] += Nuptake - symb_grow * CS_N;
    y1[ws->CS_PR_i] += Puptake - symb_grow * CS_N * red_W_P;

    y1[ws->CS_N_i] -= CSmort * CS_N;
    y1[ws->CS_NR_i] -= CSmort * CS_NR;
    y1[ws->CS_PR_i] -= CSmort * CS_PR;
    y1[ws->CS_I_i] -= CSmort * CS_I;

    y1[ws->CS_Chl_i] -= CSmort * CS_Chl;
    y1[ws->CS_Xh_i] -= CSmort * CS_Xh;
    y1[ws->CS_Xp_i] -= CSmort * CS_Xp;
    y1[ws->CS_Qi_i] -= CSmort * CS_Qi;
    y1[ws->CS_Qox_i] -= CSmort * CS_Qox;
    y1[ws->CS_Qred_i] -= CSmort * CS_Qred;
    y1[ws->CS_RO_i] -= CSmort * CS_RO;
    
    y1[ws->DetPL_N_i] += CS_N * CSmort;

    y1[ws->NH4_i] += CS_NR * CSmort / porosity ;
    y1[ws->DIP_i] +=  CS_PR * CSmort / porosity;
    y1[ws->DIC_i] +=  CS_I * CSmort * 106.0/1060.0*12.01 / porosity;

    //    if (ws->NH4_pr_i > -1)
    //  y1[ws->NH4_pr_i] += CS_NR * CSmort * SEC_PER_DAY * dz * porosity;

    if (ws->CS_N_pr_i > -1)
      y1[ws->CS_N_pr_i] += CS_N * symb_grow * SEC_PER_DAY ;

    /* Okay, mortality not part of oxygen balance - need to accumulate oxy_pr through processes, then apply to oxygen state variable */
    
    double Oxygen2 = y[ws->Oxygen_i] * y[ws->Oxygen_i]; 
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );
    
    y1[ws->Oxygen_i] -= CS_I * CSmort * 106.0/1060.0*12.01 * C_O_W * sigmoid / porosity;
    
    if (ws->COD_i > -1){
      y1[ws->COD_i] += CS_I * CSmort * 106.0/1060.0*12.01 * C_O_W * (1.0-sigmoid)/ porosity;
    }
    
    if (ws->Oxy_pr_wc_i > -1)
      y1[ws->Oxy_pr_wc_i] = CS_N * symb_grow * SEC_PER_DAY * red_W_O;
  
    /* Pigment synthesis -----------------------------------------------------------------------------------*/

    /* only synthesis if there are low numbers of inhibited RCs. */ 

    double xans = CS_Xp+CS_Xh;

    double dChldt_syn = ws->Chlmax * CSumax * Chlsynfactor * (1.0 - Iquota) * (1.0 - (CS_Qi/CS_Qt));
    double dXpdt_syn = 0.0;
    double dXhdt_syn = 0.0;

    dXpdt_syn = dChldt_syn * ws->xan2chl_CS;
    
    if (CS_N * 5.6786 < ws->C2Chlmin * CS_Chl ){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
      dXpdt_syn  = 0.0;
      dXhdt_syn  = 0.0;
    }

    // This is dilution of intracellular concentration - not total chl.

    double dChldt_dilute = symb_grow * cellChl;
    double dXhdt_dilute = symb_grow * cellXh;
    double dXpdt_dilute = symb_grow * cellXp;

    y1[ws->CS_Chl_i] += (dChldt_syn - dChldt_dilute) * ws->CSvol * cellnum;
    y1[ws->CS_Xp_i] +=  (dXpdt_syn - dXpdt_dilute) * ws->CSvol * cellnum;
    y1[ws->CS_Xh_i] +=  (dXhdt_syn - dXhdt_dilute) * ws->CSvol * cellnum;

    /* Do xanthophyll cycle adjustments. - need to asymptote to zero change at each extreme */

    if (xans > 0.0){
      double metric = (CS_Qi / CS_Qt);
      double parabole = 1.0;
      if (metric > 0.5){
	if (CS_Xh > CS_Xp)
	  parabole = 1.0-4.0*pow(CS_Xp/xans-0.5,2.0);
      }else{ 
	if (CS_Xp > CS_Xh)
	  parabole = 1.0-4.0*pow(CS_Xh/xans-0.5,2.0);
      }     

      // factor of 8 = 2^3 so that metric-0.5 ranges between -1 and 1.

      double swwitch = 8.0*pow(metric-0.5,3.0) * ws->Xanth_tau * parabole * xans;
      
      y1[ws->CS_Xp_i] -= swwitch;
      y1[ws->CS_Xh_i] += swwitch;
    }

    /*********************************************************************************************************/

    /* units: mmol, mg, m-2, s-1 .*/

    /* Rate of repair from from Qi to Qox */


    if (CS_Qt > 1.0e-9){

      double Qi2Qox = 268.0 * ws->photon2rcii * CS_Qi; // this repairs 10 mol ph m-2 d-1 hitting 1 mmol of RCII.
      
      double ARO = CSumax * Nquota * Iquota * Pquota * CS_RO;  // rate of detoxification.

      /* We are assuming that the oxygen in this cascade is not part of the oxygen balance (i.e. came from H20) */

      double absorb = kI * cellnum * ws->photon2rcii * 1000.0; //   units now mmol rcii m-2 s-1 
      
      double C_fix = absorb * (CS_Qox/CS_Qt) * y[ws->CS_tempfunc_i] * (1.0 - Iquota);
      double C_notfixed = absorb * (CS_Qox/CS_Qt) - C_fix;
      
      y1[ws->CS_Qox_i]  += dChldt_syn *  ws->CSvol * cellnum * ws->chla2rcii; // - CSmort * CS_Qox;
      y1[ws->CS_Qox_i]  += - C_notfixed + Qi2Qox;  
      y1[ws->CS_Qred_i] += C_notfixed - absorb * (CS_Qred / CS_Qt);
      y1[ws->CS_Qi_i]   += - Qi2Qox + absorb * (CS_Qred / CS_Qt);
      y1[ws->CS_RO_i]   += - ARO + 0.5 * (CS_Qi / CS_Qt) * absorb / ws->photon2rcii / ws->photon2ros;
    }
}

void symbiodinium_spectral_free_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double Coral_N =  y[ws->CS_N_i];    
  
  /* Mass balance must include all state variables introduced in this routine that hold N */

  y[ws->TN_i] += Coral_N + y[ws->CS_NR_i];
  y[ws->TP_i] += Coral_N * red_W_P + y[ws->CS_PR_i];
  y[ws->TC_i] += Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01;
  y[ws->BOD_i] += (Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01)*C_O_W;

}
