/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/carbon_leak_sed.c
 *  
 *  Description:
 *  Process implementation template
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: carbon_leak_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "einterface.h"
#include "carbon_leak_sed.h"


typedef struct {
    double KO;
    int DIC_i;
    int temp_i;
    int loc_i;
} workspace;



void carbon_leak_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    p->workspace = ws;
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->loc_i = e->find_index(tracers, "DIC_leak", e);

    double KOO = get_parameter_value(e, "DIC_leak_rate");
    ws->KO=KOO;
}

void carbon_leak_sed_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

}

void carbon_leak_sed_calc(eprocess* p, void* pp)
{


  ecology* e = p->ecology;
  void* model = e->model;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* y1 = ia->y1;
  double* y = c->y;
  double porosity = c->porosity;
  column *col = c->col;
  double ba;
  int topk_sed = einterface_getsedtopk(e->model, 0);
     if (c->k_sed>1)  // if not in the bottom layer  return
//    if (c->k_sed<topk_sed)  // if not in top layer return
	   return;
  y1[ws->DIC_i] +=(ws->KO  /   porosity/ c->dz_sed)*y[ws->loc_i];

}

void carbon_leak_sed_postcalc(eprocess* p, void* pp)
{

  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

}
