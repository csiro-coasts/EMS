/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/diffusion_epi.c
 *  
 *  Description:
 *  This procedure is supposed to account for diffusion of
 *  dissolved tracers between bottom water layer and top
 *  sediment layer. Alternatively, this process can be
 *  modelled outside ecology; however, if the top sediment
 *  layer is quite thin (as is the case with the Nugzar's
 *  models), some stuff there may be almost completely depleted
 *  during the ecology step and the model becomes quite stiff.
 *  
 *  The process is quite basic int its formulation and applies
 *  uniform settings to all tracers. It is based on the Flick's
 *  law:
 *  
 *  area flux = - D grad(conc).
 *  
 *  If we assume that the gradient is equal to difference
 *  in tracer concentrations divided by thickness of the
 *  layer in which the process takes place, it comes to
 *  
 *  area flux = - D (c1 - c2) / dz.
 *  
 *  kg m-2 s-1 = (m2 s-1) (kg m-3) m-1
 *  
 *  The process looks for two parameters, EpiDiffCoeff and
 *  EpiDiffDz, with typical values of something 2e-9 and 2e-5.
 *  Strictly speaking, D and dz do depend on temperature,
 *  porosity and, perhaps, salinity.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: diffusion_epi.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <emslogger.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "constants.h"
#include "einterface.h"
#include "diffusion_epi.h"

#define NINC 5

typedef struct {
    /*
     * coeff = EpiDiffCoeff / EpiDiffDz
     */
    double coeff;

    /*
     * tracers
     */
    int ntr;
    int* tracer_wc_i;
    int* tracer_sed_i;
} workspace;

void diffusion_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    int i;

    int OFFSET_SED = tracers->n;

    if ((e->find_index(tracers, "porosity", e) == -1)){
      quit("Porosity tracer not found");
    }

    /*
     * parameters
     */
    double EpiDiffCoeff = get_parameter_value(e, "EpiDiffCoeff");
    double EpiDiffDz = get_parameter_value(e, "EpiDiffDz");

    ws->coeff = EpiDiffCoeff / EpiDiffDz;

    /*
     * tracers
     */
    ws->ntr = 0;
    ws->tracer_wc_i = NULL;
    ws->tracer_sed_i = NULL;
    for (i = 0; i < e->ntr; ++i) {
      if (einterface_gettracerparticflag(e->model, e->tracernames[i]) ==  0 &&
          einterface_gettracerdiagnflag(e->model, e->tracernames[i]) == 0) {
        if (ws->ntr % NINC == 0) {
          ws->tracer_wc_i = realloc(ws->tracer_wc_i, 
				              (ws->ntr + NINC) * sizeof(int));
          ws->tracer_sed_i = realloc(ws->tracer_sed_i, 
				              (ws->ntr + NINC) * sizeof(int));
        }
	// Ignore salt and temp
	if ( strcmp(e->tracernames[i], "salt") != 0 &&
	     strcmp(e->tracernames[i], "temp") != 0 ) {
	  ws->tracer_wc_i[ws->ntr] = i;
	  ws->tracer_sed_i[ws->ntr] = i + OFFSET_SED;
	  ws->ntr++;
	}
      }
    }

    for (i = 0; i < ws->ntr; ++i)
      emslog(LTRACE, "  %s: %s(%d)\n", e->tracers->name, e->tracernames[ws->tracer_wc_i[i]], ws->tracer_wc_i[i]);

    p->workspace = ws;
}

void diffusion_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->tracer_wc_i);
    free(ws->tracer_sed_i);
    free(ws);
}

void diffusion_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double porosity = c->porosity;
    int i;

    for (i = 0; i < ws->ntr; ++i) {
        int iwc = ws->tracer_wc_i[i];
        int ised = ws->tracer_sed_i[i];
        double flux = -(y[iwc] - y[ised]) * ws->coeff;

        y1[iwc] += flux / c->dz_wc;
        y1[ised] -= flux / c->dz_sed / porosity;
    }
}
