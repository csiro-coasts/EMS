/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/filter_feeder_wc.c
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
 *  $Id: filter_feeder_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
//#include "einterface.h"
#include "filter_feeder_wc.h"

double ginterface_getcellz(void *model, int b, int k);
int ginterface_getwcbotk(void *model, int b);
double ginterface_cellarea(void* hmodel, int b);
  
typedef struct {
  int do_mb;                  /* flag */
  int with_df;                /* flag */
  int with_mpb;               /* flag */

/* parameters*/
  double FF_umax_t0;
  double FF_den;
  double FF_mort_t0;
  double FF_area;

  double KO_aer;


/* epis -> wcs */
    int n_Mussel_i;
    int Mussel_N_pr_i;
    int Mussel_N_rm_i;

    int Mussel_WW_i;

    int TN_i;
    int TP_i;
    int TC_i;
    int BOD_i;

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

  int Mussel_N_gr_i;

} workspace;


/* This is only called once during the lifetime of a process/ecology. */
/*******************************************************************************/
void filter_feeder_wc_init(eprocess* p)
/*******************************************************************************/

{
	ecology* e = p->ecology;
	stringtable* tracers = e->tracers;
	workspace* ws = malloc(sizeof(workspace));
        char* prm = p->prms->se[0]->s;

	p->workspace = ws;

/*parameters*/
  ws->FF_umax_t0 = get_parameter_value(e, "FFumax");
  ws->FF_mort_t0 = get_parameter_value(e, "FFmort");

  ws->FF_den = get_parameter_value(e, "FFden");
  ws->FF_area = get_parameter_value(e, "FFarea");
  
  ws->KO_aer = get_parameter_value(e, "KO_aer");
  
/* epis -> wcs */
    ws->n_Mussel_i = e->try_index(tracers, "n_Mussel", e);
    ws->Mussel_N_pr_i = e->find_index(tracers, "Mussel_N_pr", e);
    ws->Mussel_N_rm_i = e->find_index(tracers, "Mussel_N_rm", e);

    ws->Mussel_WW_i = e->find_index(tracers, "Mussel_WW", e);

    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);

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

  ws->Mussel_N_gr_i = find_index_or_add(e->cv_cell, "Mussel_N_gr", e);

}
/*******************************************************************************/
void filter_feeder_wc_postinit(eprocess* p)
/*******************************************************************************/

{
	ecology* e = p->ecology;
	workspace* ws = p->workspace;
	
	/* KWA don't think we need this so long as salmon_waste process called before remineralization*/

        ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

  if(ws->with_df)
    ws->with_df = process_present(e,PT_WC,"dinoflagellate_spectral_grow_wc");
  emstag(LINFO,"eco:filter_feeder_wc:postinit","%sCalculating filter feeder consumption of WC Dinoflagellates",(ws->with_df?"":"NOT "));
  
  if(ws->with_mpb)
    ws->with_mpb = process_present(e,PT_WC,"microphytobenthos_spectral_grow_wc");
  emstag(LINFO,"eco:filter_feeder_wc:postinit","%sCalculating filter feeder consumption of WC MPB",(ws->with_mpb?"":"NOT "));
  
}
/*******************************************************************************/
void filter_feeder_wc_destroy(eprocess* p)
/*******************************************************************************/
{
	free(p->workspace);
        
}
/*******************************************************************************/
void filter_feeder_wc_precalc(eprocess* p, void* pp)
/*******************************************************************************/
{
        ecology* e = p->ecology;
 	workspace* ws = p->workspace;
	cell* c = (cell*) pp;
	double* y = c->y;
	double* cv = c->cv;

        double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

        cv[ws->FF_umax_i] = ws->FF_umax_t0 * Tfactor;
        cv[ws->FF_mort_i] = ws->FF_mort_t0 * Tfactor;
 
    // Bail out if deeper than 24m - does this also have to be in calc??
    double z_centre = ginterface_getcellz(c->col->model,c->b,c->k_wc);
    if (z_centre < -24.0)
      return;

      if (ws->do_mb) {
      double Mussel_N_gr = cv[ws->Mussel_N_gr_i];
      y[ws->TN_i] += Mussel_N_gr;
      y[ws->TP_i] += Mussel_N_gr * red_W_P;
      y[ws->TC_i] += Mussel_N_gr * red_W_C;
      y[ws->BOD_i] += Mussel_N_gr * red_W_O;
      }

/*LOCAL DECLARATION*/
      double Mussel_WW = y[ws->Mussel_WW_i];  //  kg / farm site
 
    // Bail out if Mussel = NaN
    if (isnan(Mussel_WW))
      return;
    // Bail out if Mussel <= 0
    if ((Mussel_WW <= 0))
      return;

    /* if Mussel_WW > 1000 t /farm site assume Mussel_WW in g */
    if ((Mussel_WW/1000) > 1000)
		 Mussel_WW = Mussel_WW/1000;
    Mussel_WW = min(Mussel_WW,1000000);
	       
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
  double NH4 = y[ws->NH4_wc_i];
  double DIP = y[ws->DIP_wc_i];
  double DIC = y[ws->DIC_wc_i];

    int wcbotk; 
    wcbotk = ginterface_getwcbotk(c->col->model, c->b);
    double z_bot = ginterface_getcellz(c->col->model,c->b,wcbotk);
    double farm_volume = ginterface_cellarea(c->col->e->model, c->b) * -(max(-20,z_bot));

 // Bail out if dry cell
    if (farm_volume <= 0)
      return;

    /* mussel soft tissue dry weight in g from length in cm Navaro & Winter 1982 */

    double FF_len = 4.;
    double FF_dw = 0.0089 * pow(FF_len,2.82) ;  /* g dw */

    /* mussel shell weight in g from length in cm Navaro & Winter 1982 */

    double FF_sw = 0.550 * pow(FF_len,2.54) ;   /* g sw */

    /* assume mussel 60% water Marambio et al 2012 */

    double Mussel_N = Mussel_WW * 0.4 * (FF_dw / (FF_dw + FF_sw)) * 0.09 * 1000 * 1000 / farm_volume ; // mg N m-3

    /* number of individuals of size FF_len corresponding to biomass; assume 9% of soft tissue DW = organic N */
    y[ws->n_Mussel_i] += Mussel_N / (FF_dw * 0.09 * 1000) ;  /* number of ind m-3 - this suggests 125.3 mg N per 6 cm ind */
          
    double FFumax = cv[ws->FF_umax_i];  /* s-1 */
    
/*  SA a fraction of 1 m2 that the FF cover - not sure how to generalise to wc ropes..... 1st assume density not limiting*/
/*    double SA = ws->FF_area * (1.0-exp((-FF_WC_N/1000.) * ws->FF_den / ws->FF_area)); */
/* dimensionless, but MA_N in g, so FF_N should be too! typically 1-4%*/


    /* plankton DW ~ 2 x cell C (Lund 1964, Verduin et al, 1976) */

    double POM_N = (DetPL_N + PhyS_N + PhyL_N + PhyD_N + MPB_N + ZooS_N + ZooL_N);
    double POM_dw = (POM_N * red_W_C * 2.0) / 1000.; /* mg /L */
    double TPM_dw = POM_dw + ((Mud + FineSed) * 1000.0 * 1000.0) / 1000. ; /* mg /L */

    /* clearance rate depends on composition of sppm and is regulated up/down as a fn of organic matter content Velasco 2002; assumed FF_dw in g */
    /* CR_TPM capped at 4 l h-1 e.g. at spm < 5 mg/l */

    double CR_TPM = (min(pow(10,(2.94 - (0.61 * log10(TPM_dw)) - (1.33 * log10(POM_dw)))),4.0)/(1000*60*60)) ;  /* m3/s/g dw  */
     
    double pot_grazing = CR_TPM * POM_N * FF_dw * y[ws->n_Mussel_i] ; /* mgN m-3 s-1 */
    
    double total_grazing = min(FFumax * Mussel_N, pot_grazing); /* mg N m-3 s-1 */


   /* loss terms = pseudofaeces, mortality, (assume grazed nutrient reserves & TSS immediately released)*/
   /* pseudofaeces form when TPM > 2.7 mg/l Velasco 2002; */

    double RR_TPM = 0.0;

    if (TPM_dw > 2.7)
       RR_TPM = (pow((14.71 - 7.08 * log10(TPM_dw) * log10(POM_dw)),2)/(1000*60*60)) ;  /* m3/s/g dw */

    double pot_RR_TPM = RR_TPM * POM_N * FF_dw * y[ws->n_Mussel_i] ; /* mgN/m2/s */
    double total_RR_TPM = min(pot_RR_TPM,total_grazing); /* pseudofaeces formation < or = total grazing */

   /* assume mortality is small due to farm care, say 1% per day */

   double FFmort = 0.01/86400;  /* s-1 */
   double total_mort = FFmort * Mussel_N; /* mg N m-3 s-1 */

   /* for conservation of find mass N added to Mussels; but structural N & number of individuals remains determined by prescribed Mussel_WW */
 
   //y1[ws->FF_WC_N_i] += total_grazing - total_RR_TPM - total_mort ;
   cv[ws->Mussel_N_gr_i] += (total_grazing - total_RR_TPM - total_mort) * e->dt ;
 
   /* remove particles from WC */
   /* unit conversion to mg/m3/dt (calculation applied outside of integrator) */

    y[ws->PhyL_N_wc_i] -= total_grazing * (PhyL_N/POM_N) * e->dt ;
    y[ws->PhyL_Chl_wc_i] -= total_grazing * (PhyL_N/POM_N * PhyL_Chl/PhyL_N) * e->dt ;
    y[ws->PhyL_NR_wc_i] -= total_grazing * (PhyL_N/POM_N * PhyL_NR/PhyL_N) * e->dt ;
    y[ws->PhyL_PR_wc_i] -= total_grazing * (PhyL_N/POM_N * PhyL_PR/PhyL_N) * e->dt ;
    y[ws->PhyL_I_wc_i] -= total_grazing * (PhyL_N/POM_N * PhyL_I/PhyL_N) * e->dt  ;

    y[ws->PhyS_N_wc_i] -= total_grazing * (PhyS_N/POM_N) * e->dt  ;
    y[ws->PhyS_Chl_wc_i] -= total_grazing * (PhyS_N/POM_N * PhyS_Chl/PhyS_N) * e->dt ;
    y[ws->PhyS_NR_wc_i] -= total_grazing * (PhyS_N/POM_N * PhyS_NR/PhyS_N) * e->dt ;
    y[ws->PhyS_PR_wc_i] -= total_grazing * (PhyS_N/POM_N * PhyS_PR/PhyS_N) * e->dt ;
    y[ws->PhyS_I_wc_i] -= total_grazing * (PhyS_N/POM_N * PhyS_I/PhyS_N) * e->dt ;

    y[ws->PhyD_N_wc_i] -= total_grazing * (PhyD_N/POM_N) * e->dt ;
    y[ws->PhyD_Chl_wc_i] -= total_grazing * (PhyD_N/POM_N * PhyD_Chl/PhyD_N) * e->dt ;
    y[ws->PhyD_NR_wc_i] -= total_grazing * (PhyD_N/POM_N * PhyD_NR/PhyD_N) * e->dt ;
    y[ws->PhyD_PR_wc_i] -= total_grazing * (PhyD_N/POM_N * PhyD_PR/PhyD_N) * e->dt ;
    y[ws->PhyD_I_wc_i] -= total_grazing * (PhyD_N/POM_N * PhyD_I/PhyD_N) * e->dt ;

    y[ws->MPB_N_wc_i] -= total_grazing * (MPB_N/POM_N) * e->dt ;
    y[ws->MPB_Chl_wc_i] -= total_grazing * (MPB_N/POM_N * MPB_Chl/MPB_N) * e->dt ;
    y[ws->MPB_NR_wc_i] -= total_grazing * (MPB_N/POM_N * MPB_NR/MPB_N) * e->dt ;
    y[ws->MPB_PR_wc_i] -= total_grazing * (MPB_N/POM_N * MPB_PR/MPB_N) * e->dt ;
    y[ws->MPB_I_wc_i] -= total_grazing * (MPB_N/POM_N * MPB_I/MPB_N) * e->dt ;

    y[ws->ZooL_N_wc_i] -= total_grazing * (ZooL_N/POM_N) * e->dt ;
    y[ws->ZooS_N_wc_i] -= total_grazing * (ZooS_N/POM_N) * e->dt ;

    y[ws->DetPL_N_wc_i] += (- total_grazing * DetPL_N/POM_N  + total_RR_TPM + total_mort) * e->dt ;

   /* reserves released to NH4 & DIP & DIC */
 
    double GRF_N = (PhyL_N/POM_N * PhyL_NR/PhyL_N + PhyS_N/POM_N * PhyS_NR/PhyS_N + PhyD_N/POM_N * PhyD_NR/PhyD_N + MPB_N/POM_N * MPB_NR/MPB_N);
    double NH4_prod = total_grazing * GRF_N ;
    y[ws->NH4_wc_i] += NH4_prod * e->dt ;

    double GRF_P = (PhyL_N/POM_N * PhyL_PR/PhyL_N + PhyS_N/POM_N * PhyS_PR/PhyS_N + PhyD_N/POM_N * PhyD_PR/PhyD_N + MPB_N/POM_N * MPB_PR/MPB_N);
    y[ws->DIP_wc_i] += total_grazing * GRF_P * e->dt ;

    double GRF_C = (PhyL_N/POM_N * PhyL_I/PhyL_N + PhyS_N/POM_N * PhyS_I/PhyS_N + PhyD_N/POM_N * PhyD_I/PhyD_N + MPB_N/POM_N * MPB_I/MPB_N);
    double DIC_release = total_grazing * GRF_C * (106.0/1060.0) * 12.01 ;
    y[ws->DIC_wc_i] += DIC_release * e->dt ;

    double Oxy_pr = -DIC_release * red_W_O / red_W_C * Oxygen / (ws->KO_aer + e_max(Oxygen));
    y[ws->Oxygen_wc_i] += Oxy_pr * e->dt ;

    /* Update diagnostics */

    if (ws->Mussel_N_pr_i > -1)
      y[ws->Mussel_N_pr_i] += (total_grazing - total_RR_TPM - total_mort) * red_W_C * SEC_PER_DAY ; /* mgC m-3 d-1 */
    if (ws->Mussel_N_rm_i > -1)
      y[ws->Mussel_N_rm_i] += total_grazing * red_W_C* SEC_PER_DAY ; /* mgC m-3 d-1 */
    if (ws->Oxy_pr_wc_i  > -1)
      y[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY ;
    if (ws->NH4_pr_wc_i  > -1)
      y[ws->NH4_pr_wc_i] += NH4_prod * SEC_PER_DAY ;

}

/*******************************************************************************/
void filter_feeder_wc_calc(eprocess* p, void* pp)
/*******************************************************************************/
{ 
workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;
    double dz_wc = c->dz_wc;

}

void filter_feeder_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;
    double* cv = c->cv;
    
    double Mussel_N_gr = cv[ws->Mussel_N_gr_i];

    y[ws->TN_i] += Mussel_N_gr;
    y[ws->TP_i] += Mussel_N_gr * red_W_P;
    y[ws->TC_i] += Mussel_N_gr * red_W_C;
    y[ws->BOD_i] += Mussel_N_gr * red_W_O;

}
