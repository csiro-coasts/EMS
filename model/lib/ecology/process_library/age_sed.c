/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/age_sed.c
 *  
 *  Description:  Age tracer implentation in sediment porewaters.
 *
 *  This routine implements the derivative of the age tracer (which is potentially advected and 
 *  diffused by the hydrodynamic model). Within the source region (all sediments) the derivative is 1 d d-1.
 *
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
 *  $Id: age_sed.c 6056 2019-01-30 02:45:57Z riz008 $
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
#include "age_sed.h"

typedef struct {
  /*tracer */
  int Age_sed_i;
  /*
   * parameters
   */
  double ageing;
  
} workspace;

void age_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    
    /*get the tracer list*/
    p->workspace = ws;
    
    ws->Age_sed_i = e->find_index(tracers, "Age", e);
    
    /*get the parameters i.e age decay */

    ws->ageing = try_parameter_value(e, "ageing_sed");
    if (isnan(ws->ageing)){
      ws->ageing = 1.0/86400.0;
      eco_write_setup(e,"Code default of ageing in sediments = %e \n",ws->ageing);
    }
}

void age_sed_postinit(eprocess* p)
{


  
}

void age_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void age_sed_precalc(eprocess* p, void* pp)
{

    




}  
void age_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;


    double* y1 = ia->y1;

   
    /* Ageing in the target area */

    y1[ws->Age_sed_i] += (ws->ageing);
 
}
