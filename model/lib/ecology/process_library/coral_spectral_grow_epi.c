/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/coral_spectral_grow_epi.c
 *  
 *  Description:
 *  Implementation of simplified version of the host-symbiont model developed at UTS 
 *  by Malin Gustafsson, Mark Baird and Peter Ralph.
 *
 *  Gustafsson et al. (2013) The interchangeability of autotrophic 
 *  and hetertrophic nitrogen sources in Scleractinian coral symbiotic relationships: a 
 *  numerical study Ecol. Model.250:183-194.
 *
 *  Coral host consumes water column zooplankton, phytoplankton, organic N and 
 *  the symbiont. 
 *
 *  Symbiont grows like water column microalgae.
 *
 *  Note that host is in g N, symbiont in mg N.
 *
 *  Coral symbiont's and their hosts are both at the Redfield ratio, and release DetPL. 
 *  For the moment, lets have no DetBL grazing.
 *
 *  - Nutrient uptake by corals is reduced by the presence of macroalgae
 *  - Corals exist below macroalgae and seagrass leaves.
 *  
 *  17/09/2012 Mark Baird (UTS,WfO), Mathieu Mongin (CSIRO), Malin Gustafsson (UTS).
 *  25/01/2013 MEB (CSIRO) - Add coral calcification (implemented only
 *  if omega_ar is a tracer).
 *  23/04/2013 MEB Add oxygen terms.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: coral_spectral_grow_epi.c 7216 2022-09-18 01:20:11Z bai155 $
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
#include "coral_spectral_grow_epi.h"

double ginterface_cellarea(void* hmodel, int b);

#define unitch 1000.0

typedef struct {
  int do_mb;                  /* flag */

  /*
   * parameters
   */
  double CSumax_t0;
  double CHumax_t0;
  double Benth_resp;
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

  double CStoCHfrac;
  double CHremin;
  double Splank;

  double CHarea;

  /*
   * epis
   */

  int CH_N_i;    /* Coral host - animal */
  int CS_N_i;    /* Coral symbiont - microalgae */
  int CS_Chl_i;

  int MA_N_i;
  int ustrcw_skin_i;

  int EpiTN_i;
  int EpiTP_i;
  int EpiTC_i;
  int EpiBOD_i;

  /*
   * tracers
   */
  
  int DIC_wc_i;
  int NH4_wc_i;
  int NO3_wc_i;
  int DIP_wc_i;
  int Oxygen_wc_i;
  int Oxy_pr_wc_i;

  int ZooL_N_wc_i;
  int ZooS_N_wc_i;
  int PhyS_N_wc_i;
  int PhyL_N_wc_i;
  int DetPL_N_wc_i;

  int IN_up_i;
  int ON_up_i;
  int CS_N_pr_i;
  int CH_N_pr_i;
  int mucus_i;

  int KI_CS_i;
  int yCfac_CS_i;

 /* additional tracers if doing calcification */

  int omega_ar_i;
  int ALK_wc_i;
  int Gnet_i;

  double k_day_coral;
  double k_night_coral;

  double dissCaCO3_reef;
  double dissCaCO3_shelf;
  
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
} workspace;

void coral_spectral_grow_epi_init(eprocess* p)
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

    ws->C2Chlmin = try_parameter_value(e,"C2Chlmin");
    if (isnan(ws->C2Chlmin))
      ws->C2Chlmin = 20.0;

    ws->CSmort_t0 = get_parameter_value(e, "CSmort");
    ws->CHmort_t0 = get_parameter_value(e, "CHmort");

    ws->CHremin = get_parameter_value(e, "CHremin");
    ws->Splank = get_parameter_value(e, "Splank");

    ws->CHarea = try_parameter_value(e, "CHarea");

    /*
     * epis
     */
    ws->CH_N_i = e->find_index(epis, "CH_N", e) + OFFSET_EPI;
    ws->CS_N_i = e->find_index(epis, "CS_N", e) + OFFSET_EPI;
    ws->CS_Chl_i = e->find_index(epis, "CS_Chl", e) + OFFSET_EPI;

    ws->MA_N_i = e->try_index(epis, "MA_N", e);

    if (ws->MA_N_i > -1){
	 ws->MA_N_i += OFFSET_EPI;
	 ws->MAleafden = get_parameter_value(e, "MAleafden");
    }

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;

    ws->Gnet_i = e->find_index(epis, "Gnet", e) + OFFSET_EPI;

    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*
     * tracers
     */
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
    ws->NO3_wc_i = e->find_index(tracers, "NO3", e);
    ws->DIP_wc_i = e->find_index(tracers, "DIP", e);

    ws->ZooL_N_wc_i = e->find_index(tracers, "ZooL_N", e);
    ws->ZooS_N_wc_i = e->find_index(tracers, "ZooS_N", e);
    ws->PhyL_N_wc_i = e->find_index(tracers, "PhyL_N", e);
    ws->PhyS_N_wc_i = e->find_index(tracers, "PhyS_N", e);
    ws->DetPL_N_wc_i = e->find_index(tracers, "DetPL_N", e);

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

    /* Calcification */

    ws->omega_ar_i = e->try_index(tracers,"omega_ar", e);
    ws->ALK_wc_i = e->try_index(tracers, "alk", e);

    if (ws->omega_ar_i > -1){
       ws->k_day_coral   = get_parameter_value(e, "k_day_coral");
       ws->k_night_coral = get_parameter_value(e, "k_night_coral");
       ws->dissCaCO3_reef = get_parameter_value(e, "dissCaCO3_sed");

       ws->dissCaCO3_shelf = try_parameter_value(e, "dissCaCO3_shelf");
       if (isnan(ws->dissCaCO3_shelf)){
	   ws->dissCaCO3_shelf = 0.0001;
	 }
    }
    
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

    ws->CS_N_pr_i = e->try_index(epis, "CS_N_pr", e);
    if (ws->CS_N_pr_i > -1) 
	 ws->CS_N_pr_i += OFFSET_EPI;

    ws->CH_N_pr_i = e->try_index(epis, "CH_N_pr", e);
    if (ws->CH_N_pr_i > -1) 
	 ws->CH_N_pr_i += OFFSET_EPI;

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);
}

void coral_spectral_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;
    ws->KI_CS_i = find_index_or_add(e->cv_cell, "KI_CS", e);
    ws->yCfac_CS_i = find_index_or_add(e->cv_cell, "yCfac_CS", e);

}

void coral_spectral_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;
    free(ws);
}

void coral_spectral_grow_epi_precalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double RR;
    double area;

    if (isnan(ws->CHarea))
      ws->CHarea = 1.0;

    if (process_present(e,PT_WC,"recom_extras")){
      
      /* do calculation based on area on present cell dimension */
      
      area = ginterface_cellarea(c->col->model,c->b); 
      ws->CHarea = 1.0;
      RR = sqrt(area/PI);
      if (RR > 200.0)
	ws->CHarea = 1.0-(RR-200.0)*(RR-200.0)/(RR*RR);
    }

    /* Coupling with physical processes (temperature, shear stress, light) is calculated only 
       at the beginning of the ecological time increment */

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->CHumax_i] = ws->CHumax_t0 * Tfactor;
    cv[ws->CSumax_i] = ws->CSumax_t0 * Tfactor;
    cv[ws->CHmort_i] = ws->CHmort_t0 * Tfactor;
    cv[ws->CSmort_i] = ws->CSmort_t0 * Tfactor;

    double Coral_N = y[ws->CS_N_i] + y[ws->CH_N_i]*unitch;

    /* Mass balance must include all state variables introduced in this routine that hold N */
    if (ws->do_mb) {
      y[ws->EpiTN_i] += Coral_N;
      y[ws->EpiTP_i] += Coral_N * red_W_P;
      y[ws->EpiTC_i] += Coral_N * red_W_C;
      y[ws->EpiBOD_i] += Coral_N * red_W_O;
    }
}

void coral_spectral_grow_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;

    double* cv = c->cv;
    double dz_wc = c->dz_wc;

    double CH_N = y[ws->CH_N_i];
    double CS_N = y[ws->CS_N_i];
    double CS_Chl = y[ws->CS_Chl_i];

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

    // int wcbotk;
    // double z_bot;

    if (ws->omega_ar_i > -1){

      if ((CS_N < 1e-12)||(CH_N < 1e-9)){
	
	y1[ws->DIC_wc_i] += 12.01 *  ws->dissCaCO3_shelf / dz_wc ;
	y1[ws->ALK_wc_i] += 2.0 *  ws->dissCaCO3_shelf / dz_wc ;
	if (ws->Gnet_i> -1)
	  y1[ws->Gnet_i] -= 12.01 * ws->dissCaCO3_shelf ;  
	return;
      }
      
      if (y[ws->omega_ar_i] < 2.0){  // growth halted due to low calcification //
	CHumax = 0.0;
	CSumax = 0.0;
      }
    }

    double MA_N = 0.0;

    if (ws->MA_N_i > -1){
      MA_N = y[ws->MA_N_i];
    }
    /* Zoothanthellae - cells quantified per m2 of host tissue. */

    double cellnum = CS_N / (ws->CSm * red_A_N * 1000.0 * MW_Nitr); /* cell m-2 */
    double cellChl = CS_Chl / (ws->CSvol * cellnum);    /* mg Chl m-3  */

    /* Light absorption of cells follows phytoplankton_spectral_grow_wc.c  */

    double kI;            // mol photon cell-1 s-1
    double Chlsynfactor;  // d aA / d Chl a

    /* Convert flux per cell of photons to flux of m2 N equivalent per m2*/
    
    kI = c->cv[ws->KI_CS_i] * ( red_A_N / red_A_I ) * cellnum ;  /* mol eqN m-2 s-1 */

    Chlsynfactor = c->cv[ws->yCfac_CS_i]; /* used below for chl synthesis calculation */

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

    double SA = ws->CHarea * exp(-MA_N * ws->MAleafden) * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)); /* dimensionless */
    double kN_mass = S_DIN * SA * DIN_wc * mgN2molN;  /* mol N m-2 s-1 */
    double kP_mass = S_DIP * SA * DIP_wc * mgP2molP * (red_A_N / red_A_P); /* mol eqN m-2 s-1 */
   
    /* Growth calculated and a turnover rate (i.e. production / biomass) - use a law of the minimum that 
       includes as one of the limiting factors the maximum growth rate (see Everett et al., 2007) */

    double symb_growth = min(kI,min(kN_mass, kP_mass)) / (CS_N * mgN2molN); /* s-1 */;

    /* now how fast they are actually growing */

    double symb_growth_IN = min(CSumax,symb_growth);

    /* Coral animal - polyps quantified per m2. */

    double PhyS = y[ws->PhyS_N_wc_i];
    double PhyL = y[ws->PhyL_N_wc_i];
    double ZooS = y[ws->ZooS_N_wc_i];
    double ZooL = y[ws->ZooL_N_wc_i];
    double DetPL = y[ws->DetPL_N_wc_i];

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
    
    double CStoCHfrac = (cellnum * ws->CSrad * ws->CSrad * M_PI) / (CH_N * ws->CHpolypden * 2.0 * ws->CHarea);

    /* translocation of organic N from symbiont to host - in mg N */

    double CHremin = ws->CHremin;  /* fraction of host death that is passes to symbiont  */

    /* Mortality needs to be multipled by respective biomass to get rate
       Also, biomass concentrated so I need to divide by ws->CHarea */

    double polypmort = CHmort * CH_N / ws->CHarea;  /* s-1 */

    /* enhanced growth of symbiont due to CH remineralisation - units are s-1 */

    symb_growth = min(CSumax,symb_growth_IN + CHremin * polypmort * CH_N * unitch / CS_N);

    /* mucus production due to unwanted host death */

    double mucus = max(0.0,symb_growth_IN * CS_N + CHremin * polypmort * CH_N * unitch - CSumax * CS_N);

    /* translocate (+ve symbiont to host), including growth translocates, 
       mortality of symbionts and mortality of host (-ve). */

    double translocate = symb_growth * CS_N * CStoCHfrac + CSmort * CS_N;

    double host_growth = total_grazing + translocate / unitch;

    /* Mucus production in g N */ 
 
    double mucus2 = max(0.0,host_growth - CHumax*CH_N);

    host_growth = host_growth - mucus2;

    mucus = mucus/unitch + mucus2;

    /* Death of polyps reduces CS_N */

    y1[ws->CS_N_i] += (CS_N * symb_growth) - translocate - polypmort * CS_N;

    y1[ws->CH_N_i] += host_growth - polypmort * CH_N;

    y1[ws->NH4_wc_i] -= CS_N * symb_growth_IN * NH4_wc / DIN_wc / dz_wc;
    y1[ws->NO3_wc_i] -= CS_N * symb_growth_IN * NO3_wc / DIN_wc / dz_wc;
    y1[ws->DIP_wc_i] -= CS_N * symb_growth_IN * red_W_P / dz_wc;
    y1[ws->DIC_wc_i] -= CS_N * symb_growth_IN * red_W_C / dz_wc;

    total_grazing  = total_grazing * unitch;    

    if (Splank > 0.0){
      y1[ws->DetPL_N_wc_i] -= total_grazing * Splank * SA * DetPL / pot_grazing / dz_wc;
      y1[ws->PhyL_N_wc_i] -= total_grazing * Splank * SA * PhyL / pot_grazing / dz_wc;
      y1[ws->PhyS_N_wc_i] -= total_grazing * Splank * SA * PhyS / pot_grazing / dz_wc;
      y1[ws->ZooS_N_wc_i] -= total_grazing * Splank * SA * ZooS / pot_grazing / dz_wc;
      y1[ws->ZooL_N_wc_i] -= total_grazing * Splank * SA * ZooL / pot_grazing / dz_wc;
    }
    y1[ws->DetPL_N_wc_i] += (((1.0 - CHremin) * polypmort * CH_N + mucus) * unitch + polypmort * CS_N) / dz_wc;

    /* Oxygen production due to symbiont growth - all respiration is presently in breakdown of labile detritus */

    y1[ws->Oxygen_wc_i] += (CS_N * symb_growth_IN * red_W_O) / dz_wc;

    /* Now that we know dividing rate of CS including space limitation, do chlorophyll calc */

    double dChldt_dilute = (1.0 - CStoCHfrac) * symb_growth * cellChl;

    /* Use an abridged version of the dependence of chlorophyll synthesis on light history - in 
       fact since there is no state variable for energy reserves, Iquota based on present state 
       Need to compare light capture to nutrient uptakes (N,P) and growth rate */

    double Iquota = min(1.0,kI/min(kN_mass, kP_mass));

    Iquota = min(Iquota,kI/(CS_N * mgN2molN)/CSumax);

    Iquota  = max(Iquota,0.0);

    double tmmp = max(Chlsynfactor * (1.0 - Iquota),0.0); /* zero may be unnecessary */
    
    double dChldt_syn = ws->Chlmax * CSumax * min(tmmp, 1.33) ;

    if (CS_N * 5.6786 < ws->C2Chlmin * CS_Chl ){  /* don't synthesise if can't fit in any more. */
      dChldt_syn = 0.0;
    }

    y1[ws->CS_Chl_i] += (dChldt_syn - dChldt_dilute) * ws->CSvol * cellnum - polypmort * CS_Chl;

    y1[ws->CS_Chl_i] -= CSmort* CS_Chl; // this is in translocate for CS_N 

    /* Now do calcification and dissolution of sand among corals */

    if (ws->omega_ar_i > -1 && y[ws->omega_ar_i] > 1.0 & dz_wc > 0.2){

      /* Net community calcification - i.e. calc - diss of community */
      
      double g_night_coral = ws->k_night_coral*(y[ws->omega_ar_i]-1.0);
      double g_day_coral   = ws->k_day_coral*(y[ws->omega_ar_i]-1.0);
      
      if (c->cv[ws->KI_CS_i] > 1e-16){
	g_night_coral = 0.0;
      }else{
	g_day_coral = 0.0;
      }
      
      /* Similar to Anthony 2011, except use Iquota rather than SWR */
      
      double g0_coral = (g_day_coral * Iquota  + g_night_coral) * SA;

      /* dissolution 7 mmol m-2 d-1 Cryonak LO, convert to seconds */

      double dissCaCO3_reef = ws->dissCaCO3_reef;
      
      /* Now do the effect of calcification and dissolution on water column carbon chemistry */
      /* Use SA as a proxy of percent coral coverage */
      
      /* Alkalinity is in mmol m-3, DIC is in mg C m-3 */
      
      y1[ws->DIC_wc_i] -= 12.01 * (g0_coral - dissCaCO3_reef ) / dz_wc ;
      y1[ws->ALK_wc_i] -= 2.0 * (g0_coral - dissCaCO3_reef  )/ dz_wc ;
     
      /* output calcification rate - mg C m-2 s-1 - so it can be used in mass balance */
      if (ws->Gnet_i> -1)
	y1[ws->Gnet_i] += 12.01 * (g0_coral - dissCaCO3_reef) ;

      // printf("g0_coral %e, dissCaCO3_reef %e, y1[ws->Gnet_i] %e \n",g0_coral,dissCaCO3_reef,y1[ws->Gnet_i]);
  
    }

    // printf("SA = %e old = %e new = %e \n",SA,1.0-exp(-CH_N * ws->CHpolypden),ws->CHarea * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)));
    
    /* Update diagnostics */

    if (ws->IN_up_i > -1)
      y1[ws->IN_up_i] = S_DIN * SA * DIN_wc;
    if (ws->ON_up_i > -1)
      y1[ws->ON_up_i] = Splank * SA * (ZooS + ZooL + PhyL + PhyS + DetPL);
    if (ws->CS_N_pr_i > -1)
      y1[ws->CS_N_pr_i] += CS_N * symb_growth_IN * SEC_PER_DAY ;
    if (ws->CH_N_pr_i > -1)
      y1[ws->CH_N_pr_i] += host_growth * SEC_PER_DAY ;
    if (ws->mucus_i > -1)
      y1[ws->mucus_i] = (((1.0 - CHremin) * polypmort * CH_N + mucus) * unitch + polypmort * CS_N);
    if (ws->Oxy_pr_wc_i > -1)
      y1[ws->Oxy_pr_wc_i] = CS_N * symb_growth_IN * SEC_PER_DAY * red_W_O / dz_wc;
}

void coral_spectral_grow_epi_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  double Coral_N = y[ws->CS_N_i] + y[ws->CH_N_i]*unitch;
 
  y[ws->EpiTN_i] += Coral_N;
  y[ws->EpiTP_i] += Coral_N * red_W_P;
  y[ws->EpiTC_i] += Coral_N * red_W_C;
  y[ws->EpiBOD_i] += Coral_N * red_W_O;
}
