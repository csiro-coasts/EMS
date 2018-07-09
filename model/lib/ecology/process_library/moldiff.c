/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/moldiff.c
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
 *  $Id: moldiff.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "moldiff.h"

typedef struct {
    /*
     * common variables
     */
    int DNO3_i;
    int DNH4_i;
    int DPO4_i;
    int DN2_i;

    /*
     * tracers
     */
    int temp_i;
    int salt_i;
} workspace;

void moldiff_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * common variables
     */
    ws->DNO3_i = find_index_or_add(e->cv_cell, "DNO3", e);
    ws->DNH4_i = find_index_or_add(e->cv_cell, "DNH4", e);
    ws->DPO4_i = find_index_or_add(e->cv_cell, "DPO4", e);
    ws->DN2_i = find_index_or_add(e->cv_cell, "DN2", e);

    /*
     * tracers
     */
    ws->temp_i = e->find_index(e->tracers, "temp", e);
    ws->salt_i = e->find_index(e->tracers, "salt", e);
}

void moldiff_destroy(eprocess* p)
{
    free(p->workspace);
}

void moldiff_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;
    double diff_mult = MolDiff(1.0, y[ws->temp_i], y[ws->salt_i]);

    cv[ws->DNO3_i] = diff_mult * DNO3_25_0;
    cv[ws->DNH4_i] = diff_mult * DNH4_25_0;
    cv[ws->DPO4_i] = diff_mult * DPO4_25_0;
    cv[ws->DN2_i] = diff_mult * DN2_25_0;
}
