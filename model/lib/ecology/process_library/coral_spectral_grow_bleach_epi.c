/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/coral_spectral_grow_bleach_epi.c
 *  
 *  Description:
 * 
 *  This process contains a coral model with autotrophic and heterotrophic growth, zooxanthellae physiology, 
 *  xanthophyll cycle, reaction centre dynamics and reative oxygen build-up.
 *
 *  References:
 * 
 *  Gustafsson et al. (2013) The interchangeability of autotrophic 
 *  and hetertrophic nitrogen sources in Scleractinian coral symbiotic relationships: a 
 *  numerical study Ecol. Model.250:183-194.
 *
 *  Gustafsson et al. (2014) Modelling photoinhibition-driven bleaching in Scleractinian 
 *  corals as a function of light, temperature and heterotrophy. LO 59:603-622.
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
 *  $Id: coral_spectral_grow_bleach_epi.c 6681 2021-01-08 02:36:02Z bai155 $
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
#include "coral_spectral_grow_bleach_epi.h"

#define unitch 1000.0

typedef struct {
  int do_mb;                  /* flag */

  /*
   * parameters
   */
  double CSumax_t0;
  double CHumax_t0;
  double CHpolypden;
  double MAleafden;

  double CSm;
  double CSrad;
  double CSvol;
  double Chlmax; 
  double Plank_resp;
  double CSmort_t0;
  double CHmort_t0;
  double C2Chlmin;
  double KO_aer;

  double CStoCHfrac;
  double CHremin;
  double Splank;

  double CHarea;

  /* new parameters for xanthophyll / bleaching model */

  double Xanth_tau;
  double chla2rcii;
  double photon2rcii;
  double xan2chl_CS;
  double ROSthreshold;
  double ROSmult;
  double RubiscoOffTemp;
  double CSmaxbleachrate;
  double photon2ros;

  /*
   * epis
   */

  int CH_N_i;    /* Coral host - animal */
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

  int CS_tempfunc_i;  

  int temp_clim_wc_i;

  int MA_N_i;
  int ustrcw_skin_i;

  int EpiTN_i;
  int EpiTP_i;
  int EpiTC_i;
  int EpiBOD_i;

  /*
   * tracers
   */
  
  int temp_wc_i;
  int PAR_i;
  int DIC_wc_i;
  int NH4_wc_i;
  int NO3_wc_i;
  int DIP_wc_i;
  int Oxygen_wc_i;
  int COD_wc_i;
  int Oxy_pr_wc_i;

  int ZooL_N_wc_i;
  int ZooS_N_wc_i;
  int PhyS_N_wc_i;
  int PhyL_N_wc_i;
  int DetPL_N_wc_i;

  int PhyS_NR_wc_i;
  int PhyS_PR_wc_i;
  int PhyS_I_wc_i;
  int PhyS_Chl_wc_i; 
  int PhyL_NR_wc_i;
  int PhyL_PR_wc_i;
  int PhyL_I_wc_i;
  int PhyL_Chl_wc_i;

  int IN_up_i;
  int ON_up_i;
  int CS_N_pr_i;
  int CH_N_pr_i;
  int mucus_i;
  int CS_bleach_i;

  int KI_CS_i;
  int yCfac_CS_i;
  int omega_ar_i;

  
  /*
   * common cell variables
   */
  int DNO3_i;
  int DPO4_i;
  int Tfactor_i;
  int CHumax_i;
  int CSumax_i;
  int CHmort_i;
  int CSmort_i;
  int CHarea_cv_i;

  int recom;   // flag
} workspace;

void coral_spectral_grow_bleach_epi_init(eprocess* p)
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
    ws->CHumax_t0 = get_parameter_value(e, "CHumax");
    ws->CSumax_t0 = get_parameter_value(e, "CSumax");

    ws->CHpolypden = get_parameter_value(e, "CHpolypden");
    
    ws->CSrad = get_parameter_value(e, "CSrad");
    ws->Chlmax = PhyCellChl(ws->CSrad);
    ws->CSvol = (4.0*M_PI/3.0) * ws->CSrad * ws->CSrad * ws->CSrad;
    ws->CSm = PhyCellMass(ws->CSrad);

    eco_write_setup(e,"Calculated value of CSm = %e \n",ws->CSm);

    ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
    if (isnan(ws->C2Chlmin)){
      ws->C2Chlmin = 20.0;
      eco_write_setup(e,"Code default of C2Chlmin = %e \n",ws->C2Chlmin);
    }
    ws->CSmort_t0 = get_parameter_value(e, "CSmort");
    ws->CHmort_t0 = get_parameter_value(e, "CHmort");

    ws->CHremin = get_parameter_value(e, "CHremin");
    ws->Splank = get_parameter_value(e, "Splank");
    
    ws->CHarea = try_parameter_value(e, "CHarea");
    if (isnan(ws->CHarea)){
      ws->CHarea = 1.0;
      eco_write_setup(e,"Code default of CHarea = %e \n",ws->CHarea);
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

    ws->ROSthreshold = try_parameter_value(e, "ROSthreshold");
    if (isnan(ws->ROSthreshold)){
      ws->ROSthreshold = 1.418e-14;  // empirical fit to field observations
      eco_write_setup(e,"Code default of ROSthreshold = %e \n",ws->ROSthreshold);
    }

    ws->CSmaxbleachrate = try_parameter_value(e, "CSmaxbleachrate");
    if (isnan(ws->CSmaxbleachrate)){
      ws->CSmaxbleachrate = 1.0;        // 1 d-1 sounds reasonable.
      eco_write_setup(e,"Code default of CSmaxbleachrate = %e \n",ws->CSmaxbleachrate);
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
     * epis
     */
    ws->CH_N_i = e->find_index(epis, "CH_N", e) + OFFSET_EPI;
    ws->CS_N_i = e->find_index(epis, "CS_N", e) + OFFSET_EPI;
    ws->CS_Chl_i = e->find_index(epis, "CS_Chl", e) + OFFSET_EPI;

    ws->CS_NR_i = e->find_index(epis, "CS_NR", e) + OFFSET_EPI;
    ws->CS_PR_i = e->find_index(epis, "CS_PR", e) + OFFSET_EPI;
    ws->CS_I_i = e->find_index(epis, "CS_I", e) + OFFSET_EPI;

    ws->CS_Xh_i = e->find_index(epis, "CS_Xh", e) + OFFSET_EPI;
    ws->CS_Xp_i = e->find_index(epis, "CS_Xp", e) + OFFSET_EPI;

    ws->CS_Qred_i = e->find_index(epis, "CS_Qred", e) + OFFSET_EPI;
    ws->CS_Qox_i = e->find_index(epis, "CS_Qox", e) + OFFSET_EPI;
    ws->CS_Qi_i = e->find_index(epis, "CS_Qi", e) + OFFSET_EPI;
    ws->CS_RO_i = e->find_index(epis, "CS_RO", e) + OFFSET_EPI;

    ws->CS_tempfunc_i = e->find_index(epis, "CS_tempfunc", e) + OFFSET_EPI;

    ws->MA_N_i = e->try_index(epis, "MA_N", e);

    if (ws->MA_N_i > -1){
	 ws->MA_N_i += OFFSET_EPI;
	 ws->MAleafden = get_parameter_value(e, "MAleafden");
    }

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;

    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*
     * tracers
     */
    ws->temp_wc_i = e->find_index(tracers, "temp", e);
    ws->temp_clim_wc_i = e->find_index(tracers, "temp_clim", e);
    ws->PAR_i = e->find_index(tracers, "PAR", e);
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->COD_wc_i = e->find_index(tracers, "COD", e);
    ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);
    ws->DIP_wc_i = e->find_index(tracers, "DIP", e);

    ws->ZooL_N_wc_i = e->find_index(tracers, "ZooL_N", e);
    ws->ZooS_N_wc_i = e->find_index(tracers, "ZooS_N", e);
    ws->PhyL_N_wc_i = e->find_index(tracers, "PhyL_N", e);
    ws->PhyS_N_wc_i = e->find_index(tracers, "PhyS_N", e);
    ws->DetPL_N_wc_i = e->find_index(tracers, "DetPL_N", e);
    
    ws->PhyS_I_wc_i = e->find_index(tracers, "PhyS_I", e);
    ws->PhyS_NR_wc_i = e->find_index(tracers,"PhyS_NR", e);
    ws->PhyS_PR_wc_i = e->find_index(tracers,"PhyS_PR", e);
    ws->PhyS_Chl_wc_i = e->find_index(tracers,"PhyS_Chl", e);

    ws->PhyL_I_wc_i = e->find_index(tracers, "PhyL_I", e);
    ws->PhyL_NR_wc_i = e->find_index(tracers,"PhyL_NR", e);
    ws->PhyL_PR_wc_i = e->find_index(tracers,"PhyL_PR", e);
    ws->PhyL_Chl_wc_i = e->find_index(tracers,"PhyL_Chl", e);
    /*
     * common variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);

    ws->CHumax_i = find_index_or_add(e->cv_cell, "CHumax", e);
    ws->CSumax_i = find_index_or_add(e->cv_cell, "CSumax", e);

    ws->CHmort_i = find_index_or_add(e->cv_cell, "CHmort", e);
    ws->CSmort_i = find_index_or_add(e->cv_cell, "CSmort", e);

    ws->CHarea_cv_i = find_index_or_add(e->cv_cell, "CHarea_cv", e);

    /* Calcification */

    ws->omega_ar_i = e->try_index(tracers,"omega_ar", e);
    
    /*non essential diagnositc tracer*/

    ws->IN_up_i = e->try_index(epis, "Coral_IN_up", e);
    if (ws->IN_up_i > -1) 
	 ws->IN_up_i += OFFSET_EPI;

    ws->ON_up_i = e->try_index(epis, "Coral_ON_up", e);
    if (ws->ON_up_i > -1) 
	 ws->ON_up_i += OFFSET_EPI;

    ws->mucus_i = e->try_index(epis, "mucus", e);
    if (ws->mucus_i > -1) 
	 ws->mucus_i += OFFSET_EPI;
    
    ws->CS_bleach_i = e->try_index(epis, "CS_bleach", e);
    if (ws->CS_bleach_i > -1) 
	 ws->CS_bleach_i += OFFSET_EPI;

    ws->CS_N_pr_i = e->try_index(epis, "CS_N_pr", e);
    if (ws->CS_N_pr_i > -1) 
	 ws->CS_N_pr_i += OFFSET_EPI;

    ws->CH_N_pr_i = e->try_index(epis, "CH_N_pr", e);
    if (ws->CH_N_pr_i > -1) 
	 ws->CH_N_pr_i += OFFSET_EPI;

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

}

void coral_spectral_grow_bleach_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;
    ws->KI_CS_i = find_index_or_add(e->cv_cell, "KI_CS", e);
    ws->yCfac_CS_i = find_index_or_add(e->cv_cell, "yCfac_CS", e);

    ws->recom = 0;
    if (process_present(e,PT_WC,"recom_extras"))
      ws->recom = 1;

}

void coral_spectral_grow_bleach_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;
    free(ws);
}

void coral_spectral_grow_bleach_epi_precalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Coral_N =  y[ws->CS_N_i] + y[ws->CH_N_i]*unitch;

    double RR;
    double area;

    if (e->nstep == 0){
      if (process_present(e,PT_WC,"recom_extras")){
	if (y[ws->CS_RO_i] == 0.0){
	  y[ws->CS_Xp_i] = y[ws->CS_Chl_i]/2.0 * ws->xan2chl_CS;
	  y[ws->CS_Xh_i] = y[ws->CS_Chl_i]/2.0 * ws->xan2chl_CS;
	  y[ws->CS_Qred_i] = ws->chla2rcii * y[ws->CS_Chl_i] / 3.0; // mmol m-2
	  y[ws->CS_Qox_i] =  ws->chla2rcii * y[ws->CS_Chl_i] / 3.0; // mmol m-2
	  y[ws->CS_Qi_i] =   ws->chla2rcii * y[ws->CS_Chl_i] / 3.0; // mmol m-2
	  y[ws->CS_RO_i] = 0.0;
	  y[ws->CS_NR_i] =  y[ws->CS_N_i]/2.0;
	  y[ws->CS_PR_i] =  y[ws->CS_N_i]/2.0/16.0*32.0/14.0;
	  y[ws->CS_I_i] =  y[ws->CS_N_i]/2.0/14.0*1060.0/16.0;
	}
      }
    }
    cv[ws->CHarea_cv_i] = ws->CHarea;
    
    if (ws->recom){
      
      /* do calculation based on area on present cell dimension */
      
      area = einterface_cellarea(c->col->model,c->b); 
      cv[ws->CHarea_cv_i] = 1.0;
      RR = sqrt(area/PI);
      if (RR > 200.0)
	cv[ws->CHarea_cv_i] = 1.0-(RR-200.0)*(RR-200.0)/(RR*RR);
    }

    /* Coupling with physical processes (temperature, shear stress, light) is calculated only 
       at the beginning of the ecological time increment */

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->CHumax_i] = ws->CHumax_t0 * Tfactor;
    cv[ws->CSumax_i] = ws->CSumax_t0 * Tfactor;
    cv[ws->CHmort_i] = ws->CHmort_t0 * Tfactor;
    cv[ws->CSmort_i] = ws->CSmort_t0 * Tfactor;

    if (ws->CS_tempfunc_i > -1){
      if (y[ws->temp_clim_wc_i] > 10.0){  // sometimes zeros come thru - if so, leave as last time step.
	double deltemp = min(max(0.0,(y[ws->temp_wc_i] - max(y[ws->temp_clim_wc_i],26.0))),ws->RubiscoOffTemp);
	y[ws->CS_tempfunc_i] = (1.0 - exp(-(ws->RubiscoOffTemp-deltemp)))/(1.0-exp(-(ws->RubiscoOffTemp)));
      }
    }

    /* Mass balance must include all state variables introduced in this routine that hold N */
    if (ws->do_mb) {
      y[ws->EpiTN_i] += Coral_N + y[ws->CS_NR_i];
      y[ws->EpiTP_i] += Coral_N * red_W_P + y[ws->CS_PR_i];
      y[ws->EpiTC_i] += Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01;
      y[ws->EpiBOD_i] += (Coral_N * red_W_C + y[ws->CS_I_i] *106.0/1060.0*12.01)*C_O_W;
    }
}

void coral_spectral_grow_bleach_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;

    double CH_N = y[ws->CH_N_i];
    double CS_N = y[ws->CS_N_i];

    if ((CS_N < 1e-12)||(CH_N < 1e-9)){
      return;
    }

    double* cv = c->cv;
    double dz_wc = c->dz_wc;

    // double temp_wc = y[ws->temp_wc_i]; - now used in precalc
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

    double Nquota;
    double Pquota;
    double Iquota;

    double NO3_wc = y[ws->NO3_wc_i];
    double NH4_wc = y[ws->NH4_wc_i];
    double DIN_wc = NO3_wc + NH4_wc;
    double DIP_wc = y[ws->DIP_wc_i];
    double Sc[2];
    double DIC_wc =  y[ws->DIC_wc_i];

    double CHumax = cv[ws->CHumax_i];  /* s-1 */
    double CSumax = cv[ws->CSumax_i];  /* s-1 */
    double CHmort = cv[ws->CHmort_i];  /* s-1 */
    double CSmort = cv[ws->CSmort_i];  /* s-1 */

    double CHarea = cv[ws->CHarea_cv_i];

    if (ws->omega_ar_i > -1 && y[ws->omega_ar_i] < 2.0){  // growth halted due to low calcification //
      CHumax = 0.0;
      CSumax = 0.0;
    }

    double MA_N = 0.0;

    if (ws->MA_N_i > -1){
      MA_N = y[ws->MA_N_i];
    }

    /* Zoothanthellae - cells quantified per m2 of host tissue. */

    double cellnum = CS_N / (ws->CSm * red_A_N * 1000.0 * MW_Nitr); /* cell m-2 */
    double cellChl = CS_Chl / (ws->CSvol * cellnum);                /* mg Chl m-3  */

    /* now get normalised reserves. */

    Nquota = max(0.0,min(1.0,(CS_NR / cellnum) / PN_max));   /* dimensionless */
    Iquota = max(0.0,min(1.0,(CS_I / cellnum) / PI_max));    /* dimensionless */
    Pquota = max(0.0,min(1.0,(CS_PR / cellnum) / PP_max));    /* dimensionless */

    /* Information from light_spectral_epi.c */
    
    double kI = c->cv[ws->KI_CS_i]; // mol photon cell-1 s-1 */

    double Chlsynfactor = c->cv[ws->yCfac_CS_i]; /* used below for chl synthesis calculation */

    /*  Calculate maximum nutrient flux (Zhang 2011 Ecol. Mod. 222:1456-1470). */
    /*  ---------------------------------------------------------------------- */

    /* Schmidt number  = diffusivity momentum / diffusivity of nutrient ions */

    Sc[0] = 1.05e-6 / cv[ws->DNO3_i];
    Sc[1] = 1.05e-6 / cv[ws->DPO4_i];

    /* Nutrient uptake per m2 a function of shear stress (tau - N m-2) and Sc */

    /* Density-specific shear stress from friction velocity. Density-specific 
       avoids multiplying by density before then dividing by it below */

    double tau = y[ws->ustrcw_skin_i]*y[ws->ustrcw_skin_i];

    double S_DIN = 2850.0*pow((2.0*tau),0.38)*pow(Sc[0],-0.6)/86400.0; /* m s-1 */
    double S_DIP = 2850.0*pow((2.0*tau),0.38)*pow(Sc[1],-0.6)/86400.0; /* m s-1 */

    /* available surface area for nutrient uptake by coral - reduced by macroalgae and 
       asymptote to 1 m2 m-2 as polyp density increases  */

    /*  SHOULD WE HAVE A DIFFERENT SA FOR NUTRIENTS AND LIGHT DUE TO MACROALAGE */

    double SA = CHarea * exp(-MA_N * ws->MAleafden) * (1.0-exp(-CH_N * ws->CHpolypden / CHarea)); /* dimensionless */

    double kN_mass = S_DIN * SA * DIN_wc;  /* mg N m-2 s-1 */
    double kP_mass = S_DIP * SA * DIP_wc;  /* mg P m-2 s-1 */
    double k_NH4_mass = S_DIN * SA * NH4_wc;    /* DNO3 only 4% different to DNH4 */

    double Iuptake = 0.0;

    if (CS_Qt > 1e-9){
     Iuptake = kI * CS_Qox/CS_Qt * y[ws->CS_tempfunc_i] * (1.0 - Iquota) * cellnum * 1000.0;
    }
    double Nuptake = kN_mass*(1.0-Nquota) ; /* mg N m-2 s-1 */
    double Puptake = kP_mass*(1.0-Pquota) ; /* mg P m-2 s-1 */

    /* uptake ammonia preferientially to the mass transfer limit */

    double NH4uptake = min(k_NH4_mass,Nuptake);
    double NO3uptake = Nuptake - NH4uptake;      /* mg N m-2 s-1 */

    /* respiration based on internal reserves */

    double Iresp = CSumax * ws->Plank_resp * CS_I;  /* mmol photon m-2 s-1 */

    /* Coral animal - polyps quantified per m2. */

    double PhyS = y[ws->PhyS_N_wc_i];
    double PhyL = y[ws->PhyL_N_wc_i];
    double ZooS = y[ws->ZooS_N_wc_i];
    double ZooL = y[ws->ZooL_N_wc_i];
    double DetPL = y[ws->DetPL_N_wc_i];

    double PhyS_NR = y[ws->PhyS_NR_wc_i];
    double PhyS_PR = y[ws->PhyS_PR_wc_i];
    double PhyS_I = y[ws->PhyS_I_wc_i];
    double PhyS_Chl = y[ws->PhyS_Chl_wc_i];

    double PhyL_NR = y[ws->PhyL_NR_wc_i];
    double PhyL_PR = y[ws->PhyL_PR_wc_i];
    double PhyL_I = y[ws->PhyL_I_wc_i];
    double PhyL_Chl = y[ws->PhyL_Chl_wc_i];

    /* Uptake of animals is only weakly depended on velocity and is probably not mass transfer limited. 
       Instead specificy a constant rate transfer coefficient */

    double Splank = ws->Splank; /* Ribes and Atkinson (2007) Coral Reefs 26:413-421. */

    /* Heterotrophic feeding on water column organic N during dark only (Heidelberg et al 2003) */
    
    if (c->cv[ws->KI_CS_i] > 1e-16){  /* to avoid using light level, use absorption */
      Splank = 0.0;
    }
    
    /* solve like nutrient uptake - SA a fraction of 1 m2 that the polyps covered */

    double pot_grazing = Splank * SA * (DetPL + PhyS + PhyL + ZooS + ZooL);

    /* should maximum grazing be CHumax be multiplied by CH, but no inefficiency term as the uptake 
       rate is what they saw removed from the water column */

    double total_grazing = min(CHumax*CH_N, pot_grazing / unitch); // g N m-2 s-1

    /* Corals have two layers of zoothanthellae held with host cells. The polyp surface 
       area, CH_N * ws->CHpolypden, can exceed SA.  If the cell projected area exceeds the polyp 
       surface area * 2, then all zoothanthellae growth is delivered to the host */

    /* reduced translocation if symbiont population is small -  need to account for potential area covered by polyps */
    
    double CStoCHfrac = (cellnum * ws->CSrad * ws->CSrad * M_PI) / (CH_N * ws->CHpolypden * 2.0 * CHarea);

    /* translocation of organic N from symbiont to host - in mg N */

    double CHremin = ws->CHremin;  /* fraction of host death that is passes to symbiont  */

    /* Mortality needs to be multipled by respective biomass to get rate
       Also, biomass concentrated so I need to divide by ws->CHarea */

    double polypmort = CHmort * CH_N / CHarea;  /* s-1 */

    double symb_grow_IN = CSumax*Nquota*Iquota*Pquota; /* s-1 */

    /* enhanced growth of symbiont due to CH remineralisation - units are s-1 */

    double symb_grow = min(CSumax,symb_grow_IN + CHremin * polypmort * CH_N * unitch / CS_N);

    /* Redfield ratio mucus production due to unwanted host death - at the moment remineralised N doesn't go through reserves - negative means produce.*/

    double mucus = max(0.0,symb_grow_IN * CS_N + CHremin * polypmort * CH_N * unitch - CSumax * CS_N);

    /* Component of mucus that has reserves attached to it */

    double mucus_attNR = max(0.0,mucus - CHremin * polypmort * CH_N * unitch);

    /* translocate (+ve symbiont to host), including growth translocates, 
       mortality of symbionts and mortality of host (-ve). */

    double translocate = symb_grow * CS_N * CStoCHfrac + CSmort * CS_N;

    double translocate_attNR = symb_grow * CS_N * CStoCHfrac;

    double host_growth = total_grazing + translocate / unitch;

    /* Mucus production in g N */ 
 
    double mucus2 = max(0.0,host_growth - CHumax*CH_N);

    host_growth = host_growth - mucus2;

    mucus = mucus/unitch + mucus2;

    /* Death of polyps reduces CS_N */

    y1[ws->CS_N_i] += (CS_N * symb_grow) - translocate - polypmort * CS_N;
    y1[ws->CH_N_i] += host_growth - polypmort * CH_N;

    /* need to add lost reserves to the water column, including that from translocated cells and those lost to mucus */

    y1[ws->NH4_wc_i] += (-NH4uptake + (polypmort * CS_N * Nquota + translocate * Nquota)) / dz_wc; // + symb_grow * CS_N * (1.0 - CStoCHfrac)
    y1[ws->NO3_wc_i] += - NO3uptake / dz_wc;
    y1[ws->DIP_wc_i] += (-Puptake + (polypmort * CS_N * red_W_P * Pquota + translocate * red_W_P * Pquota))/ dz_wc;
    y1[ws->DIC_wc_i] += (polypmort * CS_N  + translocate) * red_W_C * Iquota / dz_wc;

    // NEED TO INCLUDE ANEROBIC DEMAND ON RELEASED DIC - STARTED WITH EXPULSION CASE BELOW.

    double Oxygen2 = y[ws->Oxygen_wc_i] * y[ws->Oxygen_wc_i];
    double sigmoid = Oxygen2 / (ws->KO_aer * ws->KO_aer + Oxygen2 );

    y1[ws->Oxygen_wc_i] -= (32.00 / 12.01) * (polypmort * CS_N  + translocate) * red_W_C * Iquota / dz_wc * sigmoid;

    if (ws->COD_wc_i > -1){
       y1[ws->COD_wc_i] += (32.00 / 12.01) * (polypmort * CS_N  + translocate) * red_W_C * Iquota / dz_wc * (1.0 - sigmoid);
    }
    /* respiration, carbon fixation and nitrate uptake */

    y1[ws->DIC_wc_i] += (Iresp - Iuptake) * 106.0/1060.0*12.01 / dz_wc;
    y1[ws->Oxygen_wc_i] += (-(Iresp - Iuptake) * 106.0/1060.0*32.00 + NO3uptake * 48.0/14.01) / dz_wc;

    /* Uptake from water column - growth of cells  */

    y1[ws->CS_I_i] += Iuptake - symb_grow_IN * red_A_I / (red_A_N * MW_Nitr) * CS_N - translocate_attNR * (red_A_I / (red_A_N * MW_Nitr)) * Iquota - Iresp ;
    y1[ws->CS_NR_i] += Nuptake - symb_grow_IN * CS_N - translocate_attNR * Nquota;
    y1[ws->CS_PR_i] += Puptake - symb_grow_IN * CS_N * red_W_P - translocate_attNR * red_W_P * Pquota;

    y1[ws->CS_NR_i] -= (polypmort+CSmort) * CS_NR;
    y1[ws->CS_PR_i] -= (polypmort+CSmort) * CS_PR;
    y1[ws->CS_I_i] -= (polypmort+CSmort) * CS_I;

    total_grazing  = total_grazing * unitch; /* unit change for changing water column properties */    

    if (Splank > 0.0){

      double loss = total_grazing * Splank * SA / pot_grazing / dz_wc;

      y1[ws->DetPL_N_wc_i] -= loss * DetPL ;
      y1[ws->PhyL_N_wc_i] -=  loss * PhyL ;
      y1[ws->PhyS_N_wc_i] -=  loss * PhyS ;
      y1[ws->ZooS_N_wc_i] -=  loss * ZooS ;
      y1[ws->ZooL_N_wc_i] -=  loss * ZooL ;

      //  added bits for B3p0

      y1[ws->PhyS_NR_wc_i] -= loss * PhyS_NR ;
      y1[ws->PhyS_I_wc_i] -= loss * PhyS_I;
      y1[ws->PhyS_PR_wc_i] -= loss * PhyS_PR ;
      y1[ws->PhyS_Chl_wc_i] -= loss * PhyS_Chl;
      
      y1[ws->PhyL_NR_wc_i] -= loss * PhyL_NR ;
      y1[ws->PhyL_I_wc_i] -= loss * PhyL_I;
      y1[ws->PhyL_PR_wc_i] -= loss * PhyL_PR ;
      y1[ws->PhyL_Chl_wc_i] -= loss * PhyL_Chl;

      y1[ws->DIP_wc_i] += loss * (PhyS_PR + PhyL_PR) ;
      y1[ws->NH4_wc_i] += loss * (PhyS_NR + PhyL_NR); 

      y1[ws->DIC_wc_i] += loss * (PhyS_I + PhyL_I) * 106.0/1060.0 * 12.01;
      y1[ws->Oxygen_wc_i] -= 32.00 * loss * (PhyS_I + PhyL_I) * 106.0/1060.0;

    }

    y1[ws->DetPL_N_wc_i] += (((1.0 - CHremin) * polypmort * CH_N + mucus) * unitch + polypmort * CS_N) / dz_wc;

    /* Pigment synthesis -----------------------------------------------------------------------------------*/

    /* only synthesis if there are low numbers of inhibited RCs. */ 

    double tmmp = max(Chlsynfactor * (1.0 - Iquota) * (1.0 - (CS_Qi/CS_Qt)) ,0.0); /* zero may be unnecessary */

    tmmp = min(tmmp, 1.33);

    double xans = CS_Xp+CS_Xh;

    double dChldt_syn = ws->Chlmax * CSumax * tmmp;
    double dXpdt_syn = 0.0;
    double dXhdt_syn = 0.0;

    dXpdt_syn = dChldt_syn * ws->xan2chl_CS;
    
    if (CS_N * 5.6786 < ws->C2Chlmin * CS_Chl ){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
      dXpdt_syn  = 0.0;
      dXhdt_syn  = 0.0;
    }

    // Note: no dilution of pigment due to division as abundance is quantified per m2.
      
    y1[ws->CS_Chl_i] += dChldt_syn * ws->CSvol * cellnum - (polypmort+CSmort) * CS_Chl;
    y1[ws->CS_Xp_i] +=  dXpdt_syn * ws->CSvol * cellnum - (polypmort+CSmort) * CS_Xp;
    y1[ws->CS_Xh_i] +=  dXhdt_syn * ws->CSvol * cellnum - (polypmort+CSmort) * CS_Xh;

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

    double Qi2Qox = 268.0 * ws->photon2rcii * CS_Qi; // this repairs 10 mol ph m-2 d-1 hitting 1 mmol of RCII.

    double ARO = CSumax * Nquota * Iquota * Pquota * CS_RO;  // rate of detoxification.

    /* We are assuming that the oxygen in this cascade is not part of the oxygen balance (i.e. came from H20) */

    double absorb = kI * cellnum * ws->photon2rcii * 1000.0; //   units now mmol rcii m-2 s-1 

    double C_fix = absorb * (CS_Qox/CS_Qt) * y[ws->CS_tempfunc_i] * (1.0 - Iquota);
    double C_notfixed = absorb * (CS_Qox/CS_Qt) - C_fix;

    y1[ws->CS_Qox_i]  += dChldt_syn *  ws->CSvol * cellnum * ws->chla2rcii - (polypmort+CSmort) * CS_Qox;
    y1[ws->CS_Qred_i] += - (polypmort+CSmort) * CS_Qred;
    y1[ws->CS_Qi_i]   += - (polypmort+CSmort) * CS_Qi;
    y1[ws->CS_RO_i]   += - (polypmort+CSmort) * CS_RO;

    if (CS_Qt > 1.0e-9){

      y1[ws->CS_Qox_i]  += - C_notfixed + Qi2Qox;  
      y1[ws->CS_Qred_i] += C_notfixed - absorb * (CS_Qred / CS_Qt);
      y1[ws->CS_Qi_i]   += - Qi2Qox + absorb * (CS_Qred / CS_Qt);
      y1[ws->CS_RO_i]   += - ARO + 0.5 * (CS_Qi / CS_Qt) * absorb / ws->photon2rcii / ws->photon2ros;
    }

    /* symbiont expulsion - assuming a 1:1 relationship between ROS and C destruction, and only affects 
       pigments and cell properties since they are per m2 */ 

    /* (s-1) (mol ROS cell-1) (mol C mol ROS-1) (cell m-2) */

    /* don't expel to low ROS concentrations */

    // double ROSpercell = CS_RO * (14.0 / 32.0) * (16.0 / 138.0) / (PN_max * cellnum);

    // double expulsion = max(0.0,min(1.0, ws->ROSmult * (ROSpercell - ws->ROSthreshold))) / 86400.0;

    /* redo above line in a more robust way */

    /* 1.418e-14 of the order what Suggett has (0.3 pmol cell-1 per day H2O2 production) */

    double ROSpercell = CS_RO / cellnum;

    // double expulsion = max(0.0,min(1.0,(ROSpercell - 1.418e-14)/1.418e-14)) / 86400.0;

    double expulsion = max(0.0,min(ws->CSmaxbleachrate,ws->ROSmult * (ROSpercell - ws->ROSthreshold)/ws->ROSthreshold/86400.0));


    y1[ws->CS_N_i] -= expulsion * CS_N;
    y1[ws->CS_Chl_i] -= expulsion * CS_Chl;
    y1[ws->CS_Xp_i] -= expulsion * CS_Xp;
    y1[ws->CS_Xh_i] -= expulsion * CS_Xh;
    y1[ws->CS_I_i] -= expulsion * CS_I;
    y1[ws->CS_NR_i] -= expulsion * CS_NR;
    y1[ws->CS_PR_i] -= expulsion * CS_PR;
    
    y1[ws->CS_Qox_i] -= expulsion * CS_Qox;
    y1[ws->CS_Qred_i] -= expulsion * CS_Qred;
    y1[ws->CS_Qi_i] -= expulsion * CS_Qi;

    y1[ws->CS_RO_i] -= expulsion * CS_RO;

    /* put expelled symbionts back in water column to conserve mass */

    y1[ws->DetPL_N_wc_i] += expulsion * CS_N / dz_wc;

    y1[ws->NH4_wc_i] +=  (expulsion * CS_N * Nquota) / dz_wc;
    y1[ws->DIP_wc_i] +=  (expulsion * CS_N * red_W_P * Pquota)/ dz_wc;
    y1[ws->DIC_wc_i] +=  (expulsion * CS_N * red_W_C * Iquota) / dz_wc;

    /* need to careful when oxygen approaches zero */

    y1[ws->Oxygen_wc_i] -= (expulsion * CS_N * red_W_C * Iquota * C_O_W) / dz_wc * sigmoid;

    if (ws->COD_wc_i > -1){
      y1[ws->COD_wc_i] += (expulsion * CS_N * red_W_C * Iquota) / dz_wc * C_O_W * (1.0-sigmoid);
    }
    
    if (ws->CS_bleach_i > -1)
      y1[ws->CS_bleach_i] = expulsion * 86400.0;  /* d-1 */
    
    /* Update diagnostics */

    if (ws->IN_up_i > -1)
      y1[ws->IN_up_i] = S_DIN * SA * DIN_wc;
    if (ws->ON_up_i > -1)
      y1[ws->ON_up_i] = Splank * SA * (ZooS + ZooL + PhyL + PhyS + DetPL);
    if (ws->CS_N_pr_i > -1)
      y1[ws->CS_N_pr_i] += CS_N * symb_grow * SEC_PER_DAY ;
    if (ws->CH_N_pr_i > -1)
      y1[ws->CH_N_pr_i] += host_growth * SEC_PER_DAY ;
    if (ws->mucus_i > -1)
      y1[ws->mucus_i] = (((1.0 - CHremin) * polypmort * CH_N + mucus) * unitch + polypmort * CS_N);

    /* Okay, mortality not part of oxygen balance - need to accumulate oxy_pr through processes, then apply to oxygen state variable */

    if (ws->Oxy_pr_wc_i > -1)
      y1[ws->Oxy_pr_wc_i] = CS_N * symb_grow * SEC_PER_DAY * red_W_O / dz_wc;
}

void coral_spectral_grow_bleach_epi_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double Coral_N =  y[ws->CS_N_i] + y[ws->CH_N_i]*unitch;    
  
  /* Mass balance must include all state variables introduced in this routine that hold N */
  if (ws->do_mb) {
    y[ws->EpiTN_i] += Coral_N + y[ws->CS_NR_i];
    y[ws->EpiTP_i] += Coral_N * red_W_P + y[ws->CS_PR_i];
    y[ws->EpiTC_i] += Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01;
    y[ws->EpiBOD_i] += (Coral_N * red_W_C + y[ws->CS_I_i] * 106.0/1060.0*12.01)*C_O_W;
  }
}
