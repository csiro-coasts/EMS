/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nitrification_denitrification_sed.c
 *  
 *  Description:  Sediment processes include:
 *                Nitrification 
 *                Denitrification 
 *                Anammox (anaerobic ammonium oxidation)
 *  
 *                
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: nitrification_denitrification_sed.c 6546 2020-05-06 06:58:13Z bai155 $
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
#include "nitrification_denitrification_sed.h"

/* routines in old exchange files. Should be moved to ultils.h */

double sat_O2_percent(double DO, double salt, double temp);

typedef struct {
    /*
     * parameters
     */
  double r_nit_sed_t0;
  double r_den_t0;
  double KO_Nit;
  double KO_Den;
  double KO_Amm;
  double r_amm_t0;
  
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
  int r_nit_sed_i;
  int r_den_i;
  int r_amm_i;
} workspace;

void nitrification_denitrification_sed_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  
  p->workspace = ws;
  
  /*
   * parameters
   */
  ws->r_nit_sed_t0 = get_parameter_value(e, "r_nit_sed");
  ws->r_den_t0 = get_parameter_value(e, "r_den");
  ws->KO_Nit = get_parameter_value(e, "KO_Nit");
  ws->KO_Den = get_parameter_value(e, "KO_Den");
  
  ws->KO_Amm = try_parameter_value(e,"KO_Amm");
  if (isnan(ws->KO_Amm)){
      ws->KO_Amm = 10000.0;
      eco_write_setup(e,"Anammox not implemented: add KO_Amm to parameter file to implement. \n");
  }
  ws->r_amm_t0 = try_parameter_value(e,"r_amm");
  if (isnan(ws->r_amm_t0)){
      ws->r_amm_t0 = 0.0;
      eco_write_setup(e,"Anammox not implemented: add KO_Amm to parameter file to implement. \n");
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
  ws->r_nit_sed_i = find_index_or_add(e->cv_cell, "r_nit_sed", e);
  ws->r_amm_i = find_index_or_add(e->cv_cell, "r_amm", e);
  ws->r_den_i = find_index_or_add(e->cv_cell, "r_den", e);
}

void nitrification_denitrification_sed_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  
  if( process_present(e,PT_SED,"nitrification_denitrification_fraction_sed"))
    {
      emstag(LPANIC,"eco:nitrification_denitrification_sed:postinit","nitrification_denitrification_fraction_sed seems to be enabled, these processes are mutually exclusive!");
    }
}


void nitrification_denitrification_sed_destroy(eprocess* p)
{
  free(p->workspace);
}

void nitrification_denitrification_sed_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  
  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
  
  cv[ws->r_nit_sed_i] = ws->r_nit_sed_t0 * Tfactor;
  cv[ws->r_den_i] = ws->r_den_t0 * Tfactor;
  cv[ws->r_amm_i] = ws->r_amm_t0 * Tfactor;
}

void nitrification_denitrification_sed_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  double porosity = c->porosity;
  
  double NH4 = y[ws->NH4_i];
  double NO3 = y[ws->NO3_i];
  double Oxygen = e_max(y[ws->Oxygen_i]);
  double Nitrification = cv[ws->r_nit_sed_i] * NH4 * Oxygen * Oxygen / (ws->KO_Nit * ws->KO_Nit + Oxygen * Oxygen);
  double Denitrification = cv[ws->r_den_i] * NO3 * ws->KO_Den / (ws->KO_Den + Oxygen);


  // B3p0

  // double Anammox = cv[ws->r_amm_i] * NH4 * exp(- Oxygen / ws->KO_Amm);

  // y1[ws->NH4_i] += - Nitrification - Anammox;
  // y1[ws->NO3_i] += Nitrification - Denitrification;

  double Anammox = cv[ws->r_amm_i] * sqrt(max(0.0,NH4 * NO3)) * exp(- Oxygen / ws->KO_Amm);

  y1[ws->NH4_i] += - Nitrification - Anammox/2.0;
  y1[ws->NO3_i] += Nitrification - Denitrification - Anammox/2.0;
    
  if (ws-> NH4_pr_i > -1)
    y1[ws->NH4_pr_i] -= Nitrification * SEC_PER_DAY * c->dz_sed * porosity;
  
  if (ws-> Den_fl_i > -1)
    y1[ws->Den_fl_i] += Denitrification * SEC_PER_DAY * c->dz_sed * porosity;

  if (ws-> Amm_fl_i > -1)
    y1[ws->Amm_fl_i] += Anammox * SEC_PER_DAY * c->dz_sed * porosity;
  
  /*
   * KWA added : 2 moles of DO lost for every mole of N == 4.57 g (NIT_N_0)
   *             DO per g of N nitrified and not denitrified
   * MEB added : This does not include O2 going into water. NIT_N_O has to be changed when NO3 included in O2 budget.
   */
  
  y1[ws->Oxygen_i] -= (Nitrification - Denitrification) * 48.0/14.01 ;
  
  if (ws-> Oxy_pr_i > -1)
    y1[ws->Oxy_pr_i] -= (Nitrification - Denitrification) * 48.0/14.01 *  porosity * c->dz_sed * SEC_PER_DAY;
}
void nitrification_denitrification_sed_postcalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;  
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

  if (ws->oxy_sat_i > -1)
    y[ws->oxy_sat_i] = sat_O2_percent(y[ws->Oxygen_i],y[ws->salt_i],y[ws->temp_i]);  
}
