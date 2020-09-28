/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/filter_feeder_epi.c
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
 *  $Id: filter_feeder_epi.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "stringtable.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "filter_feeder_epi.h"

  
typedef struct {
  int do_mb;                  /* flag */
  int with_df;                /* flag */
  int with_mpb;               /* flag */

/* parameters*/
  double FF_umax_t0;
  double FF_den;
  double FF_mort_t0;
  double FF_area;

  double MAleafden;
  double KO_aer;


/* epis */
    int nFF_i;
    int FF_N_i;
    int FF_N_pr_i;
    int FF_N_rm_i;

    int MA_N_i;

    int EpiTN_i;
    int EpiTP_i;
    int EpiTC_i;
    int EpiBOD_i;

/* Tracers */

  int ZooL_N_wc_i;
  int ZooS_N_wc_i;

  int PhyL_N_wc_i;
  int PhyL_NR_wc_i;
  int PhyL_I_wc_i;
  int PhyL_PR_wc_i;
  int PhyL_Chl_wc_i;

  int PhyS_N_wc_i;
  int PhyS_NR_wc_i;
  int PhyS_I_wc_i;
  int PhyS_PR_wc_i;
  int PhyS_Chl_wc_i;

  int PhyD_N_wc_i;
  int PhyD_NR_wc_i;
  int PhyD_I_wc_i;
  int PhyD_PR_wc_i;
  int PhyD_Chl_wc_i;

  int MPB_N_wc_i;
  int MPB_NR_wc_i;
  int MPB_I_wc_i;
  int MPB_PR_wc_i;
  int MPB_Chl_wc_i;

  int DetPL_N_wc_i;
  int Mud_wc_i;
  int FineSed_wc_i;

  int NH4_wc_i;
  int DIC_wc_i;
  int DIP_wc_i;
  int Oxygen_wc_i;

  /*  int temp_wc_i;
      int salt_wc_i; */
  int NH4_pr_wc_i;
  int Oxy_pr_wc_i;

   /*
   * common cell variables
   */
  int Tfactor_i;
  int FF_umax_i;
  int FF_mort_i;

} workspace;


/* This is only called once during the lifetime of a process/ecology. */
/*******************************************************************************/
void filter_feeder_epi_init(eprocess* p)
/*******************************************************************************/

{
	ecology* e = p->ecology;
	stringtable* tracers = e->tracers;
	workspace* ws = malloc(sizeof(workspace));
	 stringtable* epis = e->epis;
        int OFFSET_EPI = tracers->n * 2;

	 /* not sure if need this?
       int OFFSET_SED = tracers->n; */

	p->workspace = ws;

/*parameters*/
  ws->FF_umax_t0 = get_parameter_value(e, "FFumax");
  ws->FF_mort_t0 = get_parameter_value(e, "FFmort");

  ws->FF_den = get_parameter_value(e, "FFden");
  ws->FF_area = get_parameter_value(e, "FFarea");
  
  ws->KO_aer = get_parameter_value(e, "KO_aer");
  
/* epis */
    ws->nFF_i = e->try_index(epis, "nFF", e) + OFFSET_EPI;
    ws->FF_N_i = e->find_index(epis, "FF_N", e) + OFFSET_EPI;
    ws->FF_N_pr_i = e->find_index(epis, "FF_N_pr", e) + OFFSET_EPI;
    ws->FF_N_rm_i = e->find_index(epis, "FF_N_rm", e) + OFFSET_EPI;

    ws->MA_N_i = e->try_index(epis, "MA_N", e);

    if (ws->MA_N_i > -1){
	 ws->MA_N_i += OFFSET_EPI;
	 ws->MAleafden = get_parameter_value(e, "MAleafden");
    }

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;

/* TRACERS */
  ws->ZooS_N_wc_i = e->find_index(tracers, "ZooL_N", e);
  ws->ZooL_N_wc_i = e->find_index(tracers, "ZooL_N", e);

  ws->PhyL_N_wc_i = e->find_index(tracers, "PhyL_N", e);
  ws->PhyL_NR_wc_i = e->find_index(tracers, "PhyL_NR", e);
  ws->PhyL_I_wc_i = e->find_index(tracers, "PhyL_I", e);
  ws->PhyL_PR_wc_i = e->find_index(tracers, "PhyL_PR", e);
  ws->PhyL_Chl_wc_i = e->find_index(tracers, "PhyL_Chl", e);

  ws->PhyS_N_wc_i = e->find_index(tracers, "PhyS_N", e);
  ws->PhyS_NR_wc_i = e->find_index(tracers, "PhyS_NR", e);
  ws->PhyS_I_wc_i = e->find_index(tracers, "PhyS_I", e);
  ws->PhyS_PR_wc_i = e->find_index(tracers, "PhyS_PR", e);
  ws->PhyS_Chl_wc_i = e->find_index(tracers, "PhyS_Chl", e);

  ws->DetPL_N_wc_i = e->find_index(tracers, "DetPL_N", e);
  ws->Mud_wc_i = e->find_index(tracers, "Mud", e);
  ws->FineSed_wc_i = e->find_index(tracers, "FineSed", e);

  ws->NH4_wc_i = e->find_index(tracers, "NH4", e);
  ws->DIC_wc_i = e->find_index(tracers, "DIC", e);
  ws->DIP_wc_i = e->find_index(tracers, "DIP", e);
  ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);

  ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);
  ws->NH4_pr_wc_i = e->try_index(tracers, "NH4_pr", e);
 
 // non essential 
  ws->PhyD_N_wc_i = e->try_index(tracers, "PhyD_N", e);
  ws->PhyD_NR_wc_i = e->try_index(tracers, "PhyD_NR", e);
  ws->PhyD_I_wc_i = e->try_index(tracers, "PhyD_I", e);
  ws->PhyD_PR_wc_i = e->try_index(tracers, "PhyD_PR", e);
  ws->PhyD_Chl_wc_i = e->try_index(tracers, "PhyD_Chl", e);

  ws->MPB_N_wc_i = e->try_index(tracers, "MPB_N", e);
  ws->MPB_NR_wc_i = e->try_index(tracers, "MPB_NR", e);
  ws->MPB_I_wc_i = e->try_index(tracers, "MPB_I", e);
  ws->MPB_PR_wc_i = e->try_index(tracers, "MPB_PR", e);
  ws->MPB_Chl_wc_i = e->try_index(tracers, "MPB_Chl", e);

 if (ws->PhyD_N_wc_i > -1)
    ws->with_df = 1;
  else
    ws->with_df = 0;

  if (ws->MPB_N_wc_i > -1)
    ws->with_mpb = 1;
  else
    ws->with_mpb = 0;

 /*
   * common cell variables
   */
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->FF_umax_i = find_index_or_add(e->cv_cell, "FFumax", e);
  ws->FF_mort_i = find_index_or_add(e->cv_cell, "FFmort", e);

}
/*******************************************************************************/
void filter_feeder_epi_postinit(eprocess* p)
/*******************************************************************************/

{
	ecology* e = p->ecology;
	workspace* ws = p->workspace;
	
	/* KWA don't think we need this so long as salmon_waste process called before remineralization*/

        ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;

  if(ws->with_df)
    ws->with_df = process_present(e,PT_WC,"dinoflagellate_spectral_grow_wc");
  emstag(LINFO,"eco:zooplankton_large_spectral_grow_wc:postinit","%sCalculating filter feeder consumption of Dinoflagellates",(ws->with_df?"":"NOT "));
  
  if(ws->with_mpb)
    ws->with_mpb = process_present(e,PT_WC,"microphytobenthos_spectral_grow_wc");
  emstag(LINFO,"eco:zooplankton_large_spectral_grow_wc:postinit","%sCalculating filter feeder consumption of MPB",(ws->with_mpb?"":"NOT "));
  
}
/*******************************************************************************/
void filter_feeder_epi_destroy(eprocess* p)
/*******************************************************************************/
{
	free(p->workspace);
        
}
/*******************************************************************************/
void filter_feeder_epi_precalc(eprocess* p, void* pp)
/*******************************************************************************/
{
	workspace* ws = p->workspace;
	cell* c = (cell*) pp;
	double* y = c->y;
	double* cv = c->cv;

        double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

        cv[ws->FF_umax_i] = ws->FF_umax_t0 * Tfactor;
        cv[ws->FF_mort_i] = ws->FF_mort_t0 * Tfactor;
 
    // Bail out if deeper than 20m - does this also have to be in calc??
    double z_centre = einterface_getcellz(c->col->model,c->b,c->k_wc);
    if (z_centre < -24.0)
      return;

  if (ws->do_mb) {
    double FF_N = y[ws->FF_N_i];
    
      y[ws->EpiTN_i] += FF_N;
      y[ws->EpiTP_i] += FF_N * red_W_P;
      y[ws->EpiTC_i] += FF_N * red_W_C;
      y[ws->EpiBOD_i] += FF_N * red_W_O;
 }
}

/*******************************************************************************/
void filter_feeder_epi_calc(eprocess* p, void* pp)
/*******************************************************************************/
{ 
workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;
    double dz_wc = c->dz_wc;

    // Bail out if deeper than 20m - does this also have to be in calc??
    double z_centre = einterface_getcellz(c->col->model,c->b,c->k_wc);
    if (z_centre < -24.0)
      return;

/*LOCAL DECLARATION*/
    
  double FF_N = y[ws->FF_N_i];      /* Filter feeders assumed at Redfield Ratio; quantified mgN per m2. */
  double nFF = y[ws->nFF_i];        /* number of FF per m2 */

  double MA_N = 0.0;
    if (ws->MA_N_i > -1){
      MA_N = y[ws->MA_N_i];
    }

  double Mud = y[ws->Mud_wc_i];
  double FineSed = y[ws->FineSed_wc_i];

  double DetPL_N = y[ws->DetPL_N_wc_i];

  double ZooL_N = y[ws->ZooL_N_wc_i];
  double ZooS_N = y[ws->ZooL_N_wc_i];

  double PhyL_N = y[ws->PhyL_N_wc_i];
  double PhyL_NR = y[ws->PhyL_NR_wc_i];
  double PhyL_I = y[ws->PhyL_I_wc_i];
  double PhyL_PR = y[ws->PhyL_PR_wc_i];
  double PhyL_Chl = y[ws->PhyL_Chl_wc_i];

  double PhyS_N = y[ws->PhyS_N_wc_i];
  double PhyS_NR = y[ws->PhyS_NR_wc_i];
  double PhyS_I = y[ws->PhyS_I_wc_i];
  double PhyS_PR = y[ws->PhyS_PR_wc_i];
  double PhyS_Chl = y[ws->PhyS_Chl_wc_i];
  
  double PhyD_N,PhyD_NR,PhyD_I,PhyD_PR,PhyD_Chl;

  if (ws->with_df){
    PhyD_N = y[ws->PhyD_N_wc_i];
    PhyD_NR = y[ws->PhyD_NR_wc_i];
    PhyD_I = y[ws->PhyD_I_wc_i];
    PhyD_PR = y[ws->PhyD_PR_wc_i];
    PhyD_Chl = y[ws->PhyD_Chl_wc_i];
  }else{
    PhyD_N = 0.0;
    PhyD_NR = 0.0;
    PhyD_I = 0.0;
    PhyD_PR = 0.0;
    PhyD_Chl = 0.0;
  }

  double MPB_N,MPB_NR,MPB_I,MPB_PR,MPB_Chl;

  if (ws->with_mpb){
    MPB_N = y[ws->MPB_N_wc_i];
    MPB_NR = y[ws->MPB_NR_wc_i];
    MPB_I = y[ws->MPB_I_wc_i];
    MPB_PR = y[ws->MPB_PR_wc_i];
    MPB_Chl = y[ws->MPB_Chl_wc_i];
  }else{
    MPB_N = 0.0;
    MPB_NR = 0.0;
    MPB_I = 0.0;
    MPB_PR = 0.0;
    MPB_Chl = 0.0;
  };

  double Oxygen = y[ws->Oxygen_wc_i];




 /* like coral grazing  */

    double FFumax = cv[ws->FF_umax_i];  /* s-1 */

    /*  SA a fraction of 1 m2 that the FF cover */

    double SA = ws->FF_area * exp(-MA_N * ws->MAleafden) * (1.0-exp((-FF_N/1000.) * ws->FF_den / ws->FF_area)); /* dimensionless, but MA_N in g, so FF_N should be too! typically 1-4%*/

    /* mussel soft tissue dry weight in g from length in cm Navaro & Winter 1982 */

    double FF_len = 6.;
    double FF_dw = 0.0089 * pow(FF_len,2.82) ;  /* g dw */

    /* plankton DW ~ 2 x cell C (Lund 1964, Verduin et al, 1976) */

    double POM_N = (DetPL_N + PhyS_N + PhyL_N + PhyD_N + MPB_N + ZooS_N + ZooL_N);
    double POM_dw = (POM_N * red_W_C * 2.0) / 1000.; /* mg /L */
    double TPM_dw = POM_dw + ((Mud + FineSed) * 1000.0 * 1000.0) / 1000. ; /* mg /L */

    /* clearance rate depends on composition of sppm and is regulated up/down as a fn of organic matter content Velasco 2002; assumed FF_dw in g */
    /* CR_TPM capped at 4 l h-1 e.g. at spm < 5 mg/l */

    double CR_TPM = (min(pow(10,(2.94 - (0.61 * log10(TPM_dw)) - (1.33 * log10(POM_dw)))),4.0)/(1000*60*60)) ;  /* m3/s/g dw  */
     
    double pot_grazing = CR_TPM * POM_N * FF_dw * nFF * SA; /* mgN m-2 s-1 */
    
    double total_grazing = min(FFumax * FF_N, pot_grazing); /* mg N m-2 s-1 */


   /* loss terms = pseudofaeces, mortality, (assume grazed nutrient reserves & TSS immediately released)*/
   /* pseudofaeces form when TPM > 2.7 mg/l Velasco 2002; */

    double RR_TPM = 0.0;

    if (TPM_dw > 2.7)
       RR_TPM = (pow((14.71 - 7.08 * log10(TPM_dw) * log10(POM_dw)),2)/(1000*60*60)) ;  /* m3/s/g dw */

    double pot_RR_TPM = RR_TPM * POM_N * FF_dw * nFF * SA; /* mgN/m2/s */
    double total_RR_TPM = min(pot_RR_TPM,total_grazing*0.7); /* pseudofaeces formation < or = 90% of total grazing */

   /* mortality is quadratic */

   double FFmort = cv[ws->FF_mort_i];  /* s-1 (mgN/m2)-1 */
   double total_mort = FFmort * FF_N * FF_N; /* mg N m-2 s-1 */

   /* add partcles to FF_N */
 
    y1[ws->FF_N_i] += total_grazing - total_RR_TPM - total_mort ;

    /* increment number of individuals of size FF_len corresponding to biomass; assume 9% of soft tissue DW = organic N */
    y1[ws->nFF_i] = FF_N / (FF_dw * 0.09 * 1000) ;  /* number of ind m-2 - this suggests 125.3 mg N per 6 cm ind */

   /* remove particles from WC */
   /* adjust all fluxes to be distributed over depth of wc cell */

    y1[ws->PhyL_N_wc_i] -= total_grazing * PhyL_N/POM_N / dz_wc ;
    y1[ws->PhyL_Chl_wc_i] -= total_grazing * PhyL_N/POM_N * PhyL_Chl/PhyL_N / dz_wc;
    y1[ws->PhyL_NR_wc_i] -= total_grazing * PhyL_N/POM_N * PhyL_NR/PhyL_N / dz_wc;
    y1[ws->PhyL_PR_wc_i] -= total_grazing * PhyL_N/POM_N * PhyL_PR/PhyL_N / dz_wc;
    y1[ws->PhyL_I_wc_i] -= total_grazing * PhyL_N/POM_N * PhyL_I/PhyL_N / dz_wc ;

    y1[ws->PhyS_N_wc_i] -= total_grazing * PhyS_N/POM_N / dz_wc ;
    y1[ws->PhyS_Chl_wc_i] -= total_grazing * PhyS_N/POM_N * PhyS_Chl/PhyS_N / dz_wc;
    y1[ws->PhyS_NR_wc_i] -= total_grazing * PhyS_N/POM_N * PhyS_NR/PhyS_N / dz_wc;
    y1[ws->PhyS_PR_wc_i] -= total_grazing * PhyS_N/POM_N * PhyS_PR/PhyS_N / dz_wc;
    y1[ws->PhyS_I_wc_i] -= total_grazing * PhyS_N/POM_N * PhyS_I/PhyS_N / dz_wc ;

    y1[ws->PhyD_N_wc_i] -= total_grazing * PhyD_N/POM_N / dz_wc ;
    y1[ws->PhyD_Chl_wc_i] -= total_grazing * PhyD_N/POM_N * PhyD_Chl/PhyD_N / dz_wc;
    y1[ws->PhyD_NR_wc_i] -= total_grazing * PhyD_N/POM_N * PhyD_NR/PhyD_N / dz_wc;
    y1[ws->PhyD_PR_wc_i] -= total_grazing * PhyD_N/POM_N * PhyD_PR/PhyD_N / dz_wc;
    y1[ws->PhyD_I_wc_i] -= total_grazing * PhyD_N/POM_N * PhyD_I/PhyD_N / dz_wc ;

    y1[ws->MPB_N_wc_i] -= total_grazing * MPB_N/POM_N / dz_wc ;
    y1[ws->MPB_Chl_wc_i] -= total_grazing * MPB_N/POM_N * MPB_Chl/MPB_N / dz_wc;
    y1[ws->MPB_NR_wc_i] -= total_grazing * MPB_N/POM_N * MPB_NR/MPB_N / dz_wc;
    y1[ws->MPB_PR_wc_i] -= total_grazing * MPB_N/POM_N * MPB_PR/MPB_N / dz_wc;
    y1[ws->MPB_I_wc_i] -= total_grazing * MPB_N/POM_N * MPB_I/MPB_N / dz_wc ;

    y1[ws->ZooL_N_wc_i] -= total_grazing * ZooL_N/POM_N / dz_wc ;
    y1[ws->ZooS_N_wc_i] -= total_grazing * ZooS_N/POM_N / dz_wc ;

    y1[ws->DetPL_N_wc_i] += (- total_grazing * DetPL_N/POM_N  + total_RR_TPM + total_mort)/ dz_wc ;

   /* reserves released to NH4 & DIP & DIC */
 
    double GRF_N = (PhyL_N/POM_N * PhyL_NR/PhyL_N + PhyS_N/POM_N * PhyS_NR/PhyS_N + PhyD_N/POM_N * PhyD_NR/PhyD_N + MPB_N/POM_N * MPB_NR/MPB_N);
    double NH4_prod = total_grazing * GRF_N ;
    y1[ws->NH4_wc_i] += NH4_prod / dz_wc ;

    double GRF_P = (PhyL_N/POM_N * PhyL_PR/PhyL_N + PhyS_N/POM_N * PhyS_PR/PhyS_N + PhyD_N/POM_N * PhyD_PR/PhyD_N + MPB_N/POM_N * MPB_PR/MPB_N);
    y1[ws->DIP_wc_i] += total_grazing * GRF_P / dz_wc ;

    double GRF_C = (PhyL_N/POM_N * PhyL_I/PhyL_N + PhyS_N/POM_N * PhyS_I/PhyS_N + PhyD_N/POM_N * PhyD_I/PhyD_N + MPB_N/POM_N * MPB_I/MPB_N);
    double DIC_release = total_grazing * GRF_C * (106.0/1060.0) * 12.01 ;
    y1[ws->DIC_wc_i] += DIC_release / dz_wc;

    double Oxy_pr = -DIC_release * red_W_O / red_W_C * Oxygen / (ws->KO_aer + e_max(Oxygen));
    y1[ws->Oxygen_wc_i] += Oxy_pr / dz_wc;

    /* Update diagnostics */

    if (ws->FF_N_pr_i > -1)
      y1[ws->FF_N_pr_i] += (total_grazing - total_RR_TPM - total_mort) * red_W_C * SEC_PER_DAY ; /* mgC m-2 d-1 */
    if (ws->FF_N_rm_i > -1)
      y1[ws->FF_N_rm_i] += total_grazing * red_W_C* SEC_PER_DAY ; /* mgC m-2 d-1 */
    if (ws->Oxy_pr_wc_i  > -1)
      y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY / dz_wc;
    if (ws->NH4_pr_wc_i  > -1)
      y1[ws->NH4_pr_wc_i] += NH4_prod * SEC_PER_DAY / dz_wc;

}

void filter_feeder_epi_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;
    
    double FF_N = y[ws->FF_N_i];

    y[ws->EpiTN_i] += FF_N;
    y[ws->EpiTP_i] += FF_N * red_W_P;
    y[ws->EpiTC_i] += FF_N * red_W_C;
    y[ws->EpiBOD_i] += FF_N * red_W_O;

}
