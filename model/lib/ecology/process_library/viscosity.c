/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/viscosity.c
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
 *  $Id: viscosity.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "viscosity.h"

typedef struct {
    /*
     * common variables
     */
    int vis_i;

    /*
     * tracers
     */
    int temp_i;
    int salt_i;
} workspace;

void viscosity_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * common variables
     */
    ws->vis_i = find_index_or_add(e->cv_cell, "viscosity", e);

    /*
     * tracers
     */
    ws->temp_i = e->find_index(e->tracers, "temp", e);
    ws->salt_i = e->find_index(e->tracers, "salt", e);
}

void viscosity_destroy(eprocess* p)
{
    free(p->workspace);
}

void viscosity_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    cv[ws->vis_i] = Viscosity(y[ws->temp_i], y[ws->salt_i]);
}
