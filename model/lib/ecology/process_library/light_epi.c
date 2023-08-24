/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_epi.c
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
 *  $Id: light_epi.c 7211 2022-09-16 21:12:57Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "ecofunct.h"
#include "utils.h"
#include "eprocess.h"
#include "column.h"
#include "cell.h"
#include "light_epi.h"
#include "constants.h"

double ginterface_getlighttop(void *model, int b);

#define EPS_KD 1.0e-10

typedef struct {
    /*
     * parameters
     */
    double k_swr_par;

    /*
     * epibenthic variables
     */
    int Epilight_i;
    int Epilightatt_i;

    /*
     * common variables
     */
    int lighttop_i;
} workspace;

void light_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* epis = e->epis;
    workspace* ws = malloc(sizeof(workspace));
    int OFFSET_EPI = e->ntr * 2;

    p->workspace = ws;

    /*
     * parameters
     */
    ws->k_swr_par = try_parameter_value(e, "k_swr_par");
    if (isnan(ws->k_swr_par))
        ws->k_swr_par = SWR2PAR;

    /*
     * epis
     */
    ws->Epilight_i = e->try_index(epis, "Epilight", e);
    if (ws->Epilight_i >= 0)
      ws->Epilight_i += OFFSET_EPI;
    
    ws->Epilightatt_i = e->try_index(epis, "Epilightatt", e);
    if (ws->Epilightatt_i >= 0)
      ws->Epilightatt_i += OFFSET_EPI;

    /*
     * common column variables
     */
    ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);
}

void light_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void light_epi_precalc(eprocess* p, void* pp)
{
    cell* c = (cell*) pp;
    column* col = c->col;
    workspace* ws = p->workspace;
    double* y = c->y;
    double lighttop;
    double *cv_lighttop = col->cv[ws->lighttop_i];

    if (isnan(cv_lighttop[0]))
        cv_lighttop[0] = ginterface_getlighttop(col->model, col->b) * ws->k_swr_par;
    lighttop = cv_lighttop[0];

    /*
     * save the light intensity at the top of the epibenthos
     */
    if (ws->Epilight_i >= 0) {
        y[ws->Epilight_i] = lighttop;

        /*
         * calculate light attenuation by Macroalgae and Seagrass
         */
        if (ws->Epilightatt_i >= 0)
            lighttop = lighttop * exp(-y[ws->Epilightatt_i]);
    }

    cv_lighttop[0] = lighttop;
}
