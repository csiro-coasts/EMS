/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/mixed_layer_age_col.c
 *  
 *  Description:  Mixed_Layer_Age tracer implentation.
 *
 *  This routine implements the derivative of the age tracer (which is potentially advected and 
 *  diffused by the hydrodynamic model). Within the source region the derivative is 1 d d-1, 
 *  and the remaining area decays at a specificed rate.
 *
 *  For the case of RECOM, source area is botz < 10 m.
 *  
 *  Mongin, M., M. E. Baird. A. Lenton and S. Hadley (2016) Optimising reef-scale CO2 removal 
 *                    by seaweed to buffer ocean acidification. Environ. Res. Lett 11, 034023.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: mixed_layer_age_col.c 7212 2022-09-16 21:24:17Z bai155 $
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
#include "mixed_layer_age_col.h"

void ginterface_get_ij(void* model, int col, int *ij);
int ginterface_getwcbotk(void *model, int b);
double ginterface_getcellz(void *model, int b, int k);

typedef struct {
  /*tracer */
  int mixed_layer_i;
  int cv_mixed_layer_depth_i;
  int Mixed_Layer_Age_col_i;
  int Mixed_Layer_Age_rate_col_i;
  /*
   * parameters
   */
  double ageing;
  double anti_ageing;
  
} workspace;

void mixed_layer_age_col_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* epis = e->epis;
    
    /*get the tracer list*/
    p->workspace = ws;
    
    ws->mixed_layer_i = e->find_index(epis, "ecology_mixed_layer", e);
    
    // ws->cv_mixed_layer_depth_i = find_index_or_add(e->cv_cell, "cv_mixed_layer_depth", e);

    ws->Mixed_Layer_Age_col_i = -1;
    ws->Mixed_Layer_Age_col_i = e->find_index(tracers, "Mixed_Layer_Age", e);

    if (ws->Mixed_Layer_Age_col_i == -1)
      e_quit("No point in running process mixed_layer_age_col without tracer mixed_layer_age");
      
    /*non essential diagnostic tracer*/

    ws->Mixed_Layer_Age_rate_col_i = e->try_index(tracers, "Mixed_Layer_Age_rate", e);
    
    /*get the parameters i.e age decay and anti_age decay*/

    ws->ageing = try_parameter_value(e, "ageing_decay");
    if (isnan(ws->ageing)){
      ws->ageing = 1.0;
      eco_write_setup(e,"Code default of ageing = %e \n",ws->ageing);
    }

    ws->anti_ageing = try_parameter_value(e, "anti_ageing_decay");
    if (isnan(ws->anti_ageing)){
      ws->anti_ageing = 0.1;
      eco_write_setup(e,"Code default of anti-ageing = %e \n",ws->anti_ageing);
    }
}

void mixed_layer_age_col_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = (workspace *)p->workspace;

}

void mixed_layer_age_col_destroy(eprocess* p)
{
    free(p->workspace);
}

void mixed_layer_age_col_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  column *col = (column *) pp;
  cell* c = (cell*) pp;
  workspace* ws = p->workspace;
  double **y = col->y;
  double *zc = col->zc;
  double *dz = col->dz;
  double *y_epi = col->y_epi;
    
  double mld,zz;
  int n;
  
  // mld = c->cv[ws->cv_mixed_layer_depth_i];
  mld = y_epi[ws->mixed_layer_i];
  for (n = 0; n<col->n_wc; n++) { // 0 is top
    zz = -col->zc[n];
    if (zz < mld){
      y[ws->Mixed_Layer_Age_col_i][n] += e->dt/86400.0;    
    }
    // printf("ind %d, n %d, zz %e, mld %e age %e ecol mld %e \n",ws->cv_mixed_layer_depth_i,n,zz,mld,y[ws->Mixed_Layer_Age_col_i],y_epi[ws->mixed_layer_i]); 
  }
}
  
void mixed_layer_age_col_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* y = ia->y;
    double* y1 = ia->y1;
    double* cv = c->cv;
       
}
