/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/age_wc.c
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
 *  $Id: age_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "age_wc.h"

typedef struct {
  /*tracer */
  int Age_wc_i;
  int source_i;
  int Age_rate_wc_i ;
  /*
   * parameters
   */
  double ageing;
  double anti_ageing;
  
  //  int do_age;
  
} workspace;

void age_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    // stringtable* epis = e->epis;
    
    /*get the tracer list*/
    p->workspace = ws;
    
    ws->Age_wc_i = e->find_index(tracers, "Age", e);
    ws->source_i = e->find_index(tracers, "source", e);
    
    /*non essential diagnostic tracer*/
    ws->Age_rate_wc_i = e->try_index(tracers, "Age_rate", e);
    
    /*get the parameters i.e age decay and anti_age decay*/
    ws->ageing = get_parameter_value(e, "ageing_decay");
    ws->anti_ageing = get_parameter_value(e, "anti_ageing_decay");
}

void age_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void age_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
    ecology* e = p->ecology;

    if (process_present(e,PT_WC,"recom_extras")){

      double z_bot;
      int wcbotk;

      wcbotk = einterface_getwcbotk(c->col->model, c->b);
      z_bot = einterface_getcellz(c->col->model,c->b,wcbotk);

      if (z_bot > -10.0){
	y[ws->source_i] = 1.0; 
      }
    }

    // int ij[2];

    /*  einterface_get_ij(p->ecology->model, c->b, ij);*/
}
  
void age_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* y = ia->y;
    double* y1 = ia->y1;
     double* cv = c->cv;

     // int ij[2];

    /* Age 1 day per day in a box */

    /*  einterface_get_ij(p->ecology->model, c->b, ij);*/

    /* Anti-ageing (decay of ageing effect) */

    y1[ws->Age_wc_i] -= y[ws->Age_wc_i] * (ws->anti_ageing) * (1.0-y[ws->source_i]);
   
    /* Ageing in the target area */

    y1[ws->Age_wc_i] += (ws->ageing) * y[ws->source_i];

    if (ws-> Age_rate_wc_i > -1)
      y1[ws->Age_rate_wc_i] += (ws->ageing)*y[ws->source_i] - y[ws->Age_wc_i] * (ws->anti_ageing)*(1.0-y[ws->source_i]);  
}
