/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/gas_exchange_wc.c
 *  
 *  Description:
 *
 *  Wind speed dependent gas exchange calculation
 *
 *  Options: gas_exchange_wc(o|c,o|c|d) 
 *
 *  o - Oxygen
 *  c - Carbon dioxide
 *  d - dummy to fill 2nd argument.
 *
 *
 *  For CO2, sea-air fluxes determined in precalc, applied in calc to
 *      avoid recalculating carbon chemistry every ecological sub-step.
 *
 *  For O2, sea-air flux calculated every ecological sub-step.
 *
 *  Physical factors affecting sea-air flux (T, wind speed) used to
 *  calculated transfer coefficient *n precalc.
 *
 *  Units  Conc      Fluxes
 *
 *  DIC    mg m-3    mg m-2 s-1
 *
 *  O2     mg m-3    mg m-2 s-1
 *
 *  Calculated as a sea - air flux - negative is into the water. This
 *  is the climate community convention.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: gas_exchange_wc.c 5929 2018-09-10 06:22:56Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "eprocess.h"
#include "utils.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"

/* routines in old exchange files. Should be moved to ultils.h */

double sat_O2(double salt, double temp); // 
double sat_O2_percent(double DO, double salt, double temp);

typedef struct {

  int gas,gas1;
  
  int do_mb;
  int TC_i;
  
  int salt_i; 
  int temp_i;	   /* temperature in degree c*/

  int DIC_i;	   /* INPUT total inorganic carbon (mol/m^3) */
  int dco2star_i;  /* dco2star   = Delta CO2* (mol/m^3)*/

  int oxygen_i;   
  int oxy_sat_i;

  int O2_flux_i;   
  int CO2_flux_i;

  int layer_flag_i; /* flag to make sure gas flux is into the topmost layer thicker than 10 cm */
                  /* 0 - still to do calc this column */
                  /* 1 - do the calc in this layer */
                  /* 2 - calc already done in this column */

  /*
   * common cell variables
   */
  
  int O2_coeff_i;    /* rate coefficient determined in precalc and applied in calc */
  int oxy_conc_sat_i;   /* oxygen sat at last time point. */
  
} workspace;

void gas_exchange_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    char* prm = p->prms->se[0]->s;
    char* prm1 = p->prms->se[1]->s;

    p->workspace = ws;

    ws->gas = prm[0];
    ws->gas1 = prm1[0];

    ws->O2_flux_i = ws->oxy_sat_i = ws->CO2_flux_i = -1;
    
    if (ws->gas == 'o' || ws->gas1 == 'o' ){   /* Oxygen */   

      eco_write_setup(e,"Wind-speed dependent Oxygen fluxes resolved \n");

      ws->salt_i = e->find_index(tracers, "salt", e);
      ws->temp_i = e->find_index(tracers, "temp", e);
      ws->oxygen_i = e->find_index(tracers, "Oxygen", e);

      ws->oxy_sat_i = e->try_index(tracers, "Oxy_sat", e);
      ws->O2_flux_i = e->try_index(tracers, "O2_flux", e);

      ws->oxy_conc_sat_i = find_index_or_add(e->cv_cell, "oxy_conc_sat", e);
      ws->O2_coeff_i = find_index_or_add(e->cv_cell, "O2_coeff", e);

    }
    if (ws->gas == 'c' || ws->gas1 == 'c'){   /* Carbon dioxide */   

      eco_write_setup(e,"Wind-speed dependent Carbon fluxes resolved \n");

      ws->DIC_i = e->find_index(tracers, "DIC", e);
      ws->dco2star_i = e->find_index(tracers,"dco2star", e);
      ws->CO2_flux_i = e->find_index(tracers, "CO2_flux", e);
      
    }
    if (ws->gas != 'c' && ws->gas != 'o' && ws->gas != 'd'){
      e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected gas type: expected \"oxygen\", \"carbon\" or \"dummy\"\n", e->processfname, p->name, prm);
    }
    if (ws->gas1 != 'c' && ws->gas1 != 'o' && ws->gas1 != 'd'){
      e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected gas type: expected \"oxygen\", \"carbon\" or \"dummy\"\n", e->processfname, p->name, prm);
    }
    
    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

    if (ws->do_mb){
      ws->TC_i = e->find_index(tracers, "TC", e);
    }
    /* adds as a double - perhaps this could be changed to a int */

    ws->layer_flag_i = find_index_or_add(e->cv_column, "layer_flag", e);
}

void gas_exchange_wc_destroy(eprocess* p)
{
}

void gas_exchange_wc_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  column* col = c->col;
  double* y = c->y;

  double *cv_layer_flag = col->cv[ws->layer_flag_i];

  /* return if gas flux has already been put in the water column */
  
  if (isnan(cv_layer_flag[0]))
    cv_layer_flag[0] = 0.0;
      
  if (cv_layer_flag[0] != 0.0){
    if (ws->CO2_flux_i > -1)
      y[ws->CO2_flux_i] = 0.0;
    if (ws->O2_flux_i > -1)
      y[ws->O2_flux_i] = 0.0;
    return;
  }
  /* return if this is layer is less than 20 cm thick */
  
  if ((c->dz_wc < 0.2)){ // || (c->k_wc != einterface_getwcbotk(model, c->b)))
    if (ws->CO2_flux_i > -1)
      y[ws->CO2_flux_i] = 0.0;
    if (ws->O2_flux_i > -1)
      y[ws->O2_flux_i] = 0.0;
    return;
  }

  double temp= y[ws->temp_i];

  double U = einterface_get_windspeed(model, c->b);

  if (isnan(U))
    U = 5.0; // = gentle breeze;
  
  if (U > 20.0)
    U = 20.0;

  if (U < 0.0)
    U = 0.0;

  double Sc,salt;
    
  /* *********************************************************************
     Computes the Schmidt number of CO2 in seawater using the 
     formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
     7373-7382).  Input is temperature in deg C.
  *********************************************************************/
  
  Sc = 2073.1 - 125.62*temp + 3.6276 * pow(temp,2) - 0.043219 * pow(temp,3);
  
  /* *********************************************************************
     Compute the transfer velocity for CO2 in m/s [4]
     t       : Temperature in degrees Centigrade
     Sal         : Salinity (PSU)
     U           : Windspeed (m.s-1)
     Includes chemcical enhancement, windspeed and Schmidt number.
     Wanninkhof R. (1992). Relationship between wind speed and gas 
     exchange over the ocean. Journal of Geophysical Research
     97(C5), 7373-7382.

     Actually cubic relationship from Wanninkhof and McGillis (1999), 
     Schmidt numbers from Wanninhkof (1992).
  *********************************************************************/
  
  
  // Sea - air flux - negative is into the water.
  
  double tmp1;

  if (ws->CO2_flux_i > -1 && ws->DIC_i > -1) {
    tmp1 = ( 0.0283 / (3600.0*100.0)) * U * U * U * sqrt((660.0/Sc)); /*dco2star in mol m-3 so *1e3 to get 
                                                                                     dco2star in mmol m-3*/

    y[ws->CO2_flux_i]= - tmp1 * y[ws->dco2star_i]* 1.0e3 * 12.01;    
  }
  
  salt = y[ws->salt_i];

  c->cv[ws->oxy_conc_sat_i] = sat_O2(salt, temp);
  
  Sc = 1953.4 + 128.00 * temp + 3.9918 * temp * temp - 0.050091 * temp * temp * temp;
  
  // Sea - air flux - negative is into the water.
  
  c->cv[ws->O2_coeff_i] = ( 0.0283 / (3600.0 * 100.0) ) * U * U * U * sqrt((660.0/Sc));

  cv_layer_flag[0] = 1.0;  // so this becomes the calculation //

}

void gas_exchange_wc_calc(eprocess* p, void* pp)
{


  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  column* col = c->col;
  
  double *cv_layer_flag = col->cv[ws->layer_flag_i];

  /* only continue if this is calc layer */

  if (cv_layer_flag[0] != 1.0)
    return;

  double* y = ia->y;
  double* y1 = ia->y1;

  if (ws->CO2_flux_i > -1)
    y1[ws->DIC_i] -= y[ws->CO2_flux_i]/c->dz_wc;
  
  y1[ws->oxygen_i] -= c->cv[ws->O2_coeff_i] * (y[ws->oxygen_i] - c->cv[ws->oxy_conc_sat_i])/c->dz_wc;

  if (ws->O2_flux_i > -1)
    y1[ws->O2_flux_i] +=  c->cv[ws->O2_coeff_i] * (y[ws->oxygen_i] - c->cv[ws->oxy_conc_sat_i]);
}

void gas_exchange_wc_postcalc(eprocess* p, void* pp)
{


  cell* c = ((cell*) pp);
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;

  double *cv_layer_flag = col->cv[ws->layer_flag_i];

  if (ws->oxy_sat_i > -1){
    y[ws->oxy_sat_i] = sat_O2_percent(y[ws->oxygen_i],y[ws->salt_i],y[ws->temp_i]);
  }
  
  if (cv_layer_flag[0] == 1.0){
    cv_layer_flag[0] = 2.0;  // calc now done, so change flag 
  }
}
