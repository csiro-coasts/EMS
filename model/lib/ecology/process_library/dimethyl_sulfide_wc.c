/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dimethyl_sulfide_wc.c
 *  
 *  Description:
 *
 *  DMS model including Wind speed dependent gas exchange and generation.
 *
 *  Only works if gas_exchange is already setting up flux calculations.
 *
 *  Units  Conc      Fluxes
 *
 *  DMS    mg ??  m-3    mg ?? m-2 s-1
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
 *  $Id: dimethyl_sulfide_wc.c 6545 2020-05-06 06:57:00Z bai155 $
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

double ginterface_get_windspeed(void *model, int b);

/* routines in old exchange files. Should be moved to ultils.h */

typedef struct {
   
  // int salt_i; 
  int temp_i;	   /* temperature in degree c*/

  int dms_i;
  int dms_flux_i;

  int layer_flag_i; /* flag to make sure gas flux is into the topmost layer thicker than 10 cm */
                  /* 0 - still to do calc this column */
                  /* 1 - do the calc in this layer */
                  /* 2 - calc already done in this column */

  double PLmgN2molDMS;
  double BLmgN2molDMS;

  double r_DetBL_t0;
  double r_DetPL_t0;
  int DetPL_N_i;
  int DetBL_N_i;
  int r_DetBL_i;
  int r_DetPL_i;
  int Tfactor_i;

  /*
   * common cell variables
   */
  
  int dms_coeff_i;    /* rate coefficient determined in precalc and applied in calc */
  
} workspace;

void dimethyl_sulfide_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    ws->r_DetPL_t0 = try_parameter_value(e,"r_DetPL");
    if (isnan(ws->r_DetPL_t0)){
      ws->r_DetPL_t0 = 4.000000e-02/86400.0;
      eco_write_setup(e,"Code default of r_DetPL = %e \n",ws->r_DetPL_t0);
    }

    ws->r_DetBL_t0 = try_parameter_value(e,"r_DetBL");
    if (isnan(ws->r_DetBL_t0)){
      ws->r_DetBL_t0 = 1.0e-3/86400.0;
      eco_write_setup(e,"Code default of r_DetBL = %e \n",ws->r_DetBL_t0);
    }

    ws->PLmgN2molDMS = try_parameter_value(e,"PLmgN2molDMS");
    if (isnan(ws->PLmgN2molDMS)){
      ws->PLmgN2molDMS = 1.0;
      eco_write_setup(e,"Code default of PLmgN2molDMS = %e \n",ws->PLmgN2molDMS);
    }

    ws->BLmgN2molDMS = try_parameter_value(e,"BLmgN2molDMS");
    if (isnan(ws->BLmgN2molDMS)){
      ws->BLmgN2molDMS = 1.0;
      eco_write_setup(e,"Code default of BLmgN2molDMS = %e \n",ws->BLmgN2molDMS);
    }
    
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->DetBL_N_i = e->find_index(tracers, "DetBL_N", e);

    ws->r_DetPL_i = find_index_or_add(e->cv_cell, "r_DetPL_dms", e);
    ws->r_DetBL_i = find_index_or_add(e->cv_cell, "r_DetBL_dms", e);
      
    ws->dms_flux_i = -1;
      
    eco_write_setup(e,"Wind-speed dependent DMS fluxes resolved \n");
    
    // ws->salt_i = e->find_index(tracers, "salt", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
    ws->dms_i = e->find_index(tracers, "DMS", e);
    ws->dms_flux_i = e->try_index(tracers, "DMS_flux", e);
    ws->dms_coeff_i = find_index_or_add(e->cv_cell, "DMS_coeff", e);
    
    /* adds as a double - perhaps this could be changed to a int */

    ws->layer_flag_i = find_index_or_add(e->cv_column, "layer_flag", e);
    eco_write_setup(e,"dimethyl_sulfide_wc: layer flag %d \n",ws->layer_flag_i);
}

void dimethyl_sulfide_wc_destroy(eprocess* p)
{
}

void dimethyl_sulfide_wc_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  column* col = c->col;
  double* y = c->y;
  double* cv = c->cv;
  
  double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

  cv[ws->r_DetPL_i] = Tfactor * ws->r_DetPL_t0;
  cv[ws->r_DetBL_i] = Tfactor * ws->r_DetBL_t0;
  
  double *cv_layer_flag = col->cv[ws->layer_flag_i];

  /* return if gas flux has already been put in the water column */
  
  // if (isnan(cv_layer_flag[0]))
    // cv_layer_flag[0] = 0.0;
      
  //if (cv_layer_flag[0] != 0.0){
  //  c->cv[ws->dms_coeff_i] = 0.0;
  //  y[ws->dms_flux_i] = 0.0;
  //  return;
  //}
  /* return if this is layer is less than 20 cm thick */
  
  // if ((c->dz_wc < 0.20)){ // || (c->k_wc != ginterface_getwcbotk(model, c->b)))
  //  if (ws->dms_flux_i > -1)
  //    y[ws->dms_flux_i] = 0.0;
  //  return;
  // }

  double temp= y[ws->temp_i];

  double U = ginterface_get_windspeed(model, c->b);

  if (isnan(U))
    U = 5.0; // = gentle breeze;
  
  if (U > 20.0)
    U = 20.0;

  if (U < 5.0)
    U = 5.0;
    
  /* *********************************************************************
     Computes the Schmidt number of CO2 in seawater using the 
     formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
     7373-7382).  Input is temperature in deg C.
  *********************************************************************/
  
  // double Sc = 2073.1 - 125.62*temp + 3.6276 * pow(temp,2) - 0.043219 * pow(temp,3);
  
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
   
  double Sc = 2674.0 + 147.12 * temp + 3.726 * temp * temp - 0.038 * temp * temp * temp;
  
  // Sea - air flux - negative is into the water [m/s]
  
  c->cv[ws->dms_coeff_i] = max((2.1 * U - 2.8 ),0.0) / (3600.0 * 100.0) * sqrt((600.0/Sc));

  // printf("DMS precalc: ws->dms_coeff_i %d, temp %e, Sc %e, U %e c->cv[ws->dms_coeff_i] %e \n",ws->dms_coeff_i,temp,Sc,U,c->cv[ws->dms_coeff_i]);

  // cv_layer_flag[0] = 1.0;  // so this becomes the calculation //

}

void dimethyl_sulfide_wc_calc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  column* col = c->col;
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;
  
  /* generation of DMS */

  y1[ws->dms_i] += cv[ws->r_DetPL_i] * y[ws->DetPL_N_i] * ws->PLmgN2molDMS + cv[ws->r_DetBL_i] * y[ws->DetBL_N_i] * ws->BLmgN2molDMS;
  
  /* sea-air flux */
  
  double *cv_layer_flag = col->cv[ws->layer_flag_i];
  
  /* only continue if this is calc layer */

  if (cv_layer_flag[0] != 1.0)   
    return;

  // printf("dms: index %d, flag %e, dz %e, coeff %e, dms %e \n",ws->dms_flux_i,cv_layer_flag[0], c->dz_wc, c->cv[ws->dms_coeff_i],y[ws->dms_i]);
  
  y1[ws->dms_i] -= c->cv[ws->dms_coeff_i] * y[ws->dms_i] /c->dz_wc;

  if (ws->dms_flux_i > -1)
    y1[ws->dms_flux_i] += c->cv[ws->dms_coeff_i] * y[ws->dms_i];
}

void dimethyl_sulfide_wc_postcalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;  
  cell* c = ((cell*) pp);
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  
  // test for applying flux done in gas_exchange
  
  //double *cv_layer_flag = col->cv[ws->layer_flag_i];
 
  // if (cv_layer_flag[0] == 1.0){
  //  cv_layer_flag[0] = 2.0;  // calc now done, so change flag 
  //}
}
