/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/oxygen_exchange_wc.c
 *  
 *  Description:
 *  This procedure is supposed to replace gasExchange()
 *  from gas.c, old BM, by Stephen Walker. Here is the
 *  description of the old procedure:
 *  
 *  Gas exchange at the water surface.
 *  
 *  It is  assumed that the gas is freely available in the
 *  atmosphere as an 'infinite' source. The surface flux F is
 *  given by
 *  
 *  (G-S)
 *  F = -K ----- ,
 *  dz
 *  where K is the molecular diffusion coefficient, G is the gas
 *  concentration in the water, S is the saturation
 *  concentration, and dz is a small stagnant layer thickness.
 *  Both K and S will in general vary with temperature. dz
 *  varies with wind speed, oil slicks etc.
 *  
 *  The rate of change of concentration in the surface layer is
 *  then
 *  
 *  dG       (G-S)  1
 *  --  = -K -----  - ,
 *  dt         dz   D
 *  
 *  where D is the thickness of the surface layer. This equation
 *  implies that G exponentially approaches S (it is analogous
 *  to a decay equation in (G-S), assuming that S remains
 *  constant over a time step).
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: oxygen_exchange_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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

/*
 * must be a coefficient of molecular diffusion of oxygen in water
 */
#define DIFF_O2 2.0e-9
/*
 * stagnant layer thickness of Oxygen in water
 */
#define DZ_O2 2.0e-5

typedef struct {
    int salt_i;
    int temp_i;
    int Oxygen_i;
    int Oxy_sat_i;
} workspace;

static double sat_O2(double salt, double temp)
{
#if defined(O2_SAT_OLD)
    return 8000.0;
#else
    double t = (temp + 273.15) / 100.0;

    return 1426.73 * exp(-173.4292 + 249.6339 / t + 143.3483 * log(t) - 21.8492 * t + salt * (-0.033096 + 0.014259 * t - 0.0017 * t * t));
#endif
}

static double sat_O2_percent(double DO, double salt, double temp)
{
    double t = (temp + 273.15) / 100.0;
    return (DO/32.)*2.2414/(exp(-173.4292 + 249.6339 / t + 143.3483 * log(t) - 21.8492 * t + salt * (-0.033096 + 0.014259 * t - 0.0017 * t * t)));
    // returns DO in % (~60-120%)//
}


void oxygen_exchange_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    ws->salt_i = e->find_index(tracers, "salt", e);
    ws->temp_i = e->find_index(tracers, "temp", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Oxy_sat_i = e->try_index(tracers, "Oxy_sat", e);

}

void oxygen_exchange_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void oxygen_exchange_sed_calc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    void* model = e->model;
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    int b = c->b;
    int topk = einterface_getsedtopk(model, b);

    double salt, temp, sat, Oxygen;
    double* y;
    double* y1;


    if (c->k_sed != topk)
        return;

    y = ia->y;
    y1 = ia->y1;

    salt = y[ws->salt_i];
    temp = y[ws->temp_i];
    Oxygen = y[ws->Oxygen_i];

    sat = sat_O2(salt, temp);
    y1[ws->Oxygen_i] += (sat - Oxygen) * DIFF_O2 / DZ_O2 / c->dz_sed;

}


void oxygen_exchange_sed_postcalc(eprocess* p, void* pp)
{

  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  double salt, temp, Oxygen;

  if(ws->Oxy_sat_i <= -1)
    return;


  salt = y[ws->salt_i];
  temp = y[ws->temp_i];
  Oxygen = y[ws->Oxygen_i];
  y[ws->Oxy_sat_i] = sat_O2_percent(Oxygen,salt,temp);

}
