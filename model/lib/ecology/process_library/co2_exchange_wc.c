/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/co2_exchange_wc.c
 *  
 *  Description:
 *  Process co2_exchange
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: co2_exchange_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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

typedef struct {

  //int  XkW_i;    /*  piston velocity*/

  int co2ex_i;      /* air to sea co2 fluxes in mol m-2 s-1 */
  int TEMP_i;	  /* temperature in degree c*/
  int DIC_i;	  /* INPUT total inorganic carbon (mol/m^3)*/
  int dco2star_i;   /* dco2star   = Delta CO2* (mol/m^3)*/
  
  double co2exx;    /* rate of exchange determined in precalc and applied in calc */
  
} workspace;

void co2_exchange_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    // ws->XkW_i =e->find_index(tracers, "xKw", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->TEMP_i = e->find_index(e->tracers, "temp", e);

    ws->co2ex_i = e->try_index(tracers,"co2ex", e);
    ws->dco2star_i = e->try_index(tracers,"dco2star", e);
}

void co2_exchange_wc_destroy(eprocess* p)
{
}


void co2_exchange_wc_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* y = c->y;

  /* only continue if in the top watercolumn layer */

  if (c->k_wc != einterface_getwctopk(model, c->b))
    return;

  double U = einterface_get_windspeed(model, c->b);
    
  /* define for co2 air sea fluxes*/

  double scco2;      /*internal parameter*/ 

  double temp= y[ws->TEMP_i];
  double dco2star=y[ws->dco2star_i];
   
  /*********************************************************************/	
  /*Calculate air to sea fluxes */
  /******************************************************************** 
  
  *********************************************************************
     Computes the Schmidt number of CO2 in seawater using the 
     formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
     7373-7382).  Input is temperature in deg C.
     t       : Temperature in degrees Centigrade
  *********************************************************************/
  
  scco2 = 2073.1 - 125.62*temp + 3.6276*pow(temp,2) - 0.043219*pow(temp,3);
  
  /* 
*********************************************************************
Compute the transfer velocity for CO2 in m/s [4]
t       : Temperature in degrees Centigrade
Sal         : Salinity (PSU)
U           : Windspeed (m.s-1)
Includes chemcical enhancement, windspeed and Schmidt number.
Wanninkhof R. (1992). Relationship between wind speed and gas 
exchange over the ocean. Journal of Geophysical Research
97(C5), 7373-7382.
  *********************************************************************/
  
  
  /*the xKW is the piston velocity ~
    [Xconv*a*(u2+v)] in ms-1
    u2 is the instantaneous ssmi wind speed v is the variance a is coefficient this follow the OCMIP guideline adpated from wanninkhof 1992.
    is read from a file */

  ws->co2exx = 0.0283 / (100 * 3600 ) * U * U * U * sqrt((660/scco2))*dco2star*1e3*12.01; /*dco2star in mol ,-2 so *1e3 to get in mmol m-3*/
   
  // printf("in precalc, U = %e, gas exchange = %e \n",U,ws->co2exx);

}

void co2_exchange_wc_calc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);

  /* only continue if in the top watercolumn layer */

  if (c->k_wc != einterface_getwctopk(model, c->b))
    return;

  double* y1 = ia->y1;
  // double* y = ia->y;

  // printf("in calc, gas exchange = %e, \n",ws->co2exx);

  y1[ws->DIC_i] += ws->co2exx/c->dz_wc; 
  //  y[ws->co2ex_i] = ws->co2exx/c->dz_wc;
}

void co2_exchange_wc_postcalc(eprocess* p, void* pp)
{ 
}
