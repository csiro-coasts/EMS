/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/variable_parameter.c
 *  
 *  Description: rotuine to modify mBGC model parameter when used as diagnostics
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: variable_parameter.c 7198 2022-09-14 06:02:52Z bai155 $
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
//#include "einterface.h"
#include "variable_parameter.h"

int ginterface_getwcbotk(void *model, int b);
double ginterface_getcellz(void *model, int b, int k);

typedef struct {

 
  int PhyL_sv_i;


} workspace;

void variable_parameter_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));

  p->workspace = ws;

  ws->PhyL_sv_i = e->find_index(tracers, "PhyL_sv", e);


}




void variable_parameter_destroy(eprocess* p)
{
  workspace* ws = (workspace*) p->workspace;
  free(ws);
}

void variable_parameter_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* cv = c->cv;
  double* y = c->y;

}

void variable_parameter_calc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  cell* c = ((cell*) ia->media);
  double* cv = c->cv;
  double* y = ia->y;
  double* y1 = ia->y1;

}
void variable_parameter_postcalc(eprocess* p, void* pp)
{



  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

 
  double z_centre,z_bot;
  int wcbotk;
  if ((ws->PhyL_sv_i > -1) )
{

  wcbotk = ginterface_getwcbotk(c->col->model, c->b);
    z_bot = ginterface_getcellz(c->col->model,c->b,wcbotk);
    
  if (z_bot > -200.0){
    y[ws->PhyL_sv_i] = 0.0;}
      else{
	y[ws->PhyL_sv_i] =  -0.0001157;
    }
 } 
}

