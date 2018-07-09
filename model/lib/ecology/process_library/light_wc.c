/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_wc.c
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
 *  $Id: light_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ems.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "utils.h"
#include "eprocess.h"
#include "column.h"
#include "cell.h"
#include "einterface.h"
#include "light_wc.h"
#include "constants.h"

#define EPS_KD 1.0e-10


typedef struct {
    /*
     * parameters
     */
    double k_swr_par;
    
    /*
     * tracers
     */
    int Light_i;
    int Kd_i;

    /*
     * common variables
     */
    int lighttop_i;
    
} workspace;

void light_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->k_swr_par = try_parameter_value(e, "k_swr_par");
    if (isnan(ws->k_swr_par))
        ws->k_swr_par = SWR2PAR;
        
    /*
     * tracers
     */
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Kd_i = e->find_index(tracers, "Kd", e);
    
    
    /*
     * common column variables
     */
    ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);
    
}

void light_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;

    if( process_present(e,PT_WC,"light_wc_gradient"))
    {
      emstag(LPANIC,"eco:light_wc:postinit","light_wc_gradient seems to be enabled, these processes are mutually exclusive!");
    }
    
}

void light_wc_destroy(eprocess* p)
{
    free(p->workspace);
}

void light_wc_precalc(eprocess* p, void* pp)
{
    cell* c = (cell*) pp;
    column* col = c->col;
    workspace* ws = p->workspace;
    double* y = c->y;
    double dz = c->dz_wc;
    double lighttop, lightbot;

    double Kd = (y[ws->Kd_i] > 0.0) ? y[ws->Kd_i] : EPS_KD;
    double *cv_lighttop = col->cv[ws->lighttop_i];
		
    if (isnan(cv_lighttop[0]))
        cv_lighttop[0] = einterface_getlighttop(col->model, col->b) * ws->k_swr_par;
    lighttop = cv_lighttop[0];
    lightbot = lighttop * exp(-Kd * dz);
    y[ws->Light_i] = (lighttop - lightbot) / Kd / dz;

    cv_lighttop[0]= lightbot;
    
}


