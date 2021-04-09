/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/values_common_epi.c
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
 *  $Id: values_common_epi.c 6692 2021-03-24 01:04:45Z wil00y $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "values_common_epi.h"

typedef struct {
    /*
     * common variables
     */

  int ustrcw_skin_i;
} workspace;

void values_common_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    int OFFSET_SED = tracers->n;
    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    /*
     * tracers
     */

    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*non essentail diagnostic tracer*/
}

void values_common_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void values_common_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* y = c->y;

    if (isnan(y[ws->ustrcw_skin_i]))
      y[ws->ustrcw_skin_i] = 0.0;
    /*
     * (already multipled by dz_sed)
     */
}

void values_common_epi_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* y = c->y;

}
