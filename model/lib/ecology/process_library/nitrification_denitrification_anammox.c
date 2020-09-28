/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nitrification_denitrification_anammox.c
 *  
 *  Description:  Bacterial-mediated nitrogen / oxygen pool changes in both sediment and water column due to:
 *                Nitrification 
 *                Denitrification 
 *                Anammox (anaerobic ammonium oxidation)
 * 
 *  Note that dentirification and anammox remove N from the system. However, the oxygen liberated 
 *  during anammox is released into the water coulmn as dissolved oxygen while the N becomes N2 gas. 
 *                
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: nitrification_denitrification_anammox.c 6190 2019-03-28 23:27:59Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "nitrification_denitrification_anammox.h"

/* routines in old exchange files. Should be moved to ultils.h */

double sat_O2_percent(double DO, double salt, double temp);

typedef struct {
    /*
     * parameters
     */
  double r_nit_t0;
  double r_den_t0;
  double r_ana_t0;
  double KO_nit;
  double KO_den;
  double KO_ana;
  
  /*
   * tracers
   */
  int NH4_i;
  int NH4_pr_i;
  int NO3_i;
  int Oxygen_i;
  int COD_i;
  int Oxy_pr_i;
  int Den_fl_i;
  int Amm_fl_i;
  
  int salt_i; 
  int temp_i;	   /* temperature in degree c*/
  
  int oxy_sat_i;

  /*
   * common variables
   */
  int Tfactor_i;
  int r_nit_i;
  int r_den_i;
  int r_ana_i;
} workspace;

void nitrification_denitrification_anammox_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*
   * parameters
   */

  if (p->type == PT_SED){

    eco_write_setup(e,"nitrification_denitrification_anammox called in sediments \n");

    ws->r_nit_t0 = get_parameter_value(e, "r_nit_sed");
    ws->r_den_t0 = get_parameter_value(e, "r_den_sed");
    ws->KO_nit = get_parameter_value(e, "KO_nit_sed");
    ws->KO_den = get_parameter_value(e, "KO_den_sed");
    
    ws->KO_ana = try_parameter_value(e,"KO_ana_sed");
    if (isnan(ws->KO_ana)){
      ws->KO_ana = 10000.0;
      eco_write_setup(e,"Code default of ws->KO_ana_sed = %e \n",ws->KO_ana);
    }
    ws->r_ana_t0 = try_parameter_value(e,"r_ana_sed");
    if (isnan(ws->r_ana_t0)){
      ws->r_ana_t0 = 0.0;
     eco_write_setup(e,"Code default of  ws->r_ana_sed = %e \n", ws->r_ana_t0);
    }
  }
  
  if (p->type == PT_WC){

    eco_write_setup(e,"nitrification_denitrification_anammox called in water column \n");

    ws->r_nit_t0 = get_parameter_value(e, "r_nit_wc");
    ws->r_den_t0 = get_parameter_value(e, "r_den_wc");
    ws->KO_nit = get_parameter_value(e, "KO_nit_wc");
    ws->KO_den = get_parameter_value(e, "KO_den_wc");
   
    ws->KO_ana = try_parameter_value(e,"KO_ana_wc");
    if (isnan(ws->KO_ana)){
      ws->KO_ana = 10000.0;
      eco_write_setup(e,"Code default of ws->KO_ana_wc = %e \n",ws->KO_ana);
    }
    ws->r_ana_t0 = try_parameter_value(e,"r_ana_wc");
    if (isnan(ws->r_ana_t0)){
      ws->r_ana_t0 = 0.0;
     eco_write_setup(e,"Code default of  ws->r_ana_wc = %e \n", ws->r_ana_t0);
    }
    
  }
    
  /*
   * tracers
   */
  ws->NH4_i = e->find_index(tracers, "NH4", e);
  ws->NO3_i = e->find_index(tracers, "NO3", e);
  ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
  ws->COD_i = e->try_index(tracers, "COD", e);
  
  ws->salt_i = e->find_index(tracers, "salt", e);
  ws->temp_i = e->find_index(tracers, "temp", e);
  
  /*non essential diagnostic tracers*/
  ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
  ws->Den_fl_i = e->try_index(tracers, "Den_fl", e);
  ws->Amm_fl_i = e->try_index(tracers, "Amm_fl", e);
  ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
  ws->oxy_sat_i = e->try_index(tracers, "Oxy_sat", e);
  /*
   * common variables
   */
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->r_nit_i = find_index_or_add(e->cv_cell, "r_nit", e);
  ws->r_ana_i = find_index_or_add(e->cv_cell, "r_ana", e);
  ws->r_den_i = find_index_or_add(e->cv_cell, "r_den", e);
}

void nitrification_denitrification_anammox_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  if (p->type == PT_SED){
    if ( process_present(e,PT_SED,"nitrification_denitrification_fraction_sed")|| 
	 process_present(e,PT_SED,"nitrification_denitrification_sed"))    {
      emstag(LPANIC,"eco:nitrification_denitrification_anammox:postinit","nitrification_denitrification_fraction_sed enabled, these processes are mutually exclusive!");
    }
  }
  if (p->type == PT_WC){
    if ( process_present(e,PT_WC,"nitrification_wc")){
      emstag(LPANIC,"eco:nitrification_denitrification_anammox:postinit","nitrification_wc enabled, these processes are mutually exclusive!");
    }
  }
}

void nitrification_denitrification_anammox_destroy(eprocess* p)
{
  free(p->workspace);
}

void nitrification_denitrification_anammox_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  
  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
  
  cv[ws->r_nit_i] = ws->r_nit_t0 * Tfactor;
  cv[ws->r_den_i] = ws->r_den_t0 * Tfactor;
  cv[ws->r_ana_i] = ws->r_ana_t0 * Tfactor;
}

void nitrification_denitrification_anammox_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double NH4 = y[ws->NH4_i];
  double NO3 = y[ws->NO3_i];
  double Oxygen = e_max(y[ws->Oxygen_i]);

  double Nitrification = 0.0;
  if (NH4 > 0.0001) {
    Nitrification = cv[ws->r_nit_i] * NH4 * Oxygen * Oxygen / (ws->KO_nit * ws->KO_nit + Oxygen * Oxygen);}

  double Denitrification = 0.0;
  if (NO3 > 0.0001) {
    Denitrification =  cv[ws->r_den_i] * NO3 * ws->KO_den / (ws->KO_den + Oxygen);}

  double Anammox = 0.0;
  /* if ((NH4 > 0.0001) && (NO3 > 0.0001)) {
     Anammox = cv[ws->r_ana_i] * sqrt(max(0.0,NH4 * NO3)) * exp(- Oxygen / ws->KO_ana);} */

  y1[ws->NH4_i] += - Nitrification - Anammox/2.0;
  y1[ws->NO3_i] += Nitrification - Denitrification - Anammox/2.0;

  double dz = 0.0; // layer thickness x porosity.

  if (p->type == PT_WC){
    dz = c->dz_wc;
  }else{
    dz = c->dz_sed * c->porosity;
  }
    
  if (ws-> NH4_pr_i > -1)
    y1[ws->NH4_pr_i] -= Nitrification * SEC_PER_DAY * dz;
  
  if (ws-> Den_fl_i > -1)
    y1[ws->Den_fl_i] += Denitrification * SEC_PER_DAY * dz;

  if (ws-> Amm_fl_i > -1)
    y1[ws->Amm_fl_i] += Anammox * SEC_PER_DAY * dz;
  
  y1[ws->Oxygen_i] -= (Nitrification - Denitrification - Anammox/2.0) * 48.0/14.01 ;
  
  if (ws-> Oxy_pr_i > -1)
    y1[ws->Oxy_pr_i] -= (Nitrification - Denitrification - Anammox/2.0) * 48.0/14.01 * dz * SEC_PER_DAY;
}
void nitrification_denitrification_anammox_postcalc(eprocess* p, void* pp)
{


  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

  if (ws->oxy_sat_i > -1)
    y[ws->oxy_sat_i] = sat_O2_percent(y[ws->Oxygen_i],y[ws->salt_i],y[ws->temp_i]);  
}
