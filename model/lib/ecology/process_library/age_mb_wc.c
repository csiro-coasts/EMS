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
 *  $Id: age_mb_wc.c 5846 2018-06-29 04:14:26Z riz008 $
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
 
  int Age_wc_i;
  int do_age;

} workspace;

void age_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;
    ws->Age_wc_i = e->find_index(tracers, "Age", e);

}

void age_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void age_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    
    int ij[2];

    einterface_get_ij(p->ecology->model, c->b, ij);

    ws->do_age = 0;

      if (ij[0] > 26){
	if (ij[0] < 32){
	  if (ij[1] > 14){
	    if(ij[1] < 18){
	      ws->do_age = 1;
	    }
	  }
	}
      }
}

void age_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* y = ia->y;
    double* y1 = ia->y1;

    int ij[2];

    /* Age 1 day per day in a box */

    einterface_get_ij(p->ecology->model, c->b, ij);

    /* Anti-ageing (decay of ageing effect) */

    y1[ws->Age_wc_i] -= y[ws->Age_wc_i]/86400.0/10.0;

    /* Ageing in the target area */
    
    if (ws->do_age){
      y1[ws->Age_wc_i] += 1.0/86400.0;
    }
}
