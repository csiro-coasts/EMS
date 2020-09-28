/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/gas_exchange_epi.c
 *  
 *  Description:
 *  Bubble formation due to epibenthic plants and MPB.
 * 
 *  Calculated as a sea - air flux - negative is into the water. 
 *  This is the climate community convention.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: gas_exchange_epi.c 5846 2018-06-29 04:14:26Z riz008 $
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

  int salt_i; 
  int temp_i;	   /* temperature in degree c*/

  int oxygen_i;   
  int oxy_sat_i;
  int O2_flux_i;

  int oxy_conc_sat_i;
  
} workspace;

void gas_exchange_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    ws->salt_i = e->find_index(tracers, "salt", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
    ws->oxygen_i = e->find_index(tracers, "Oxygen", e);

    ws->oxy_sat_i = e->try_index(tracers, "Oxy_sat", e);
    ws->O2_flux_i = e->try_index(tracers, "O2_flux", e);

    ws->oxy_conc_sat_i = find_index_or_add(e->cv_cell, "oxy_conc_sat", e);

}

void gas_exchange_epi_destroy(eprocess* p)
{
}

void gas_exchange_epi_precalc(eprocess* p, void* pp)
{


  workspace* ws = p->workspace;
  cell* c = (cell*) pp;

  double* y = c->y;

  /* assume that the oxygen saturation is the same for plants and MPB */

  c->cv[ws->oxy_conc_sat_i] = sat_O2(y[ws->salt_i],y[ws->temp_i]);
}

void gas_exchange_epi_calc(eprocess* p, void* pp)
{


  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);

 
  double* y = ia->y;
  double* y1 = ia->y1;

  /* Bubbles keep oxygen below 120 % saturation over 1 hour */

  double deltaO2 = max(0.0,(y[ws->oxygen_i] - c->cv[ws->oxy_conc_sat_i] * 1.2)/3600.0);

  y1[ws->oxygen_i] -= deltaO2;
  y1[ws->O2_flux_i] += deltaO2*c->dz_wc;
}

void gas_exchange_epi_postcalc(eprocess* p, void* pp)
{






}
