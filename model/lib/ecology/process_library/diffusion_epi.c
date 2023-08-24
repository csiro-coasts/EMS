/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/diffusion_epi.c
 *  
 *  Description:
 *  This procedure accounts for diffusion of
 *  dissolved tracers between bottom water layer and the top
 *  sediment layer. Alternatively, this process can be
 *  modelled outside ecology; however, if the top sediment
 *  layer is quite thin, the tracer in the sediment becomes completely depleted
 *  during the ecology step and the model becomes quite stiff.
 *  
 *  The process is quite basic in its formulation and applies
 *  uniform settings to all tracers. It is based on the Fick's
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
 *  If the tracer list contains a 2D variable <tracername>_sedflux, then the 
 *  flux of that tracer due to this processes is outputted. 
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: diffusion_epi.c 7199 2022-09-14 06:36:12Z bai155 $
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
//#include "einterface.h"
#include "diffusion_epi.h"

int ginterface_gettracerdiagnflag(void* model, char* name);
int ginterface_gettracerparticflag(void* model, char* name);

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

  int fluxout_oxygen_i;
  int flux_oxygen_wc_i;
  int flux_oxygen_sed_i;

  int fluxout_nitrate_i;
  int flux_nitrate_wc_i;
  int flux_nitrate_sed_i;

  int fluxout_dip_i;
  int flux_dip_wc_i;
  int flux_dip_sed_i;

  int fluxout_amm_i;
  int flux_amm_wc_i;
  int flux_amm_sed_i;

  int fluxout_dic_i;
  int flux_dic_wc_i;
  int flux_dic_sed_i;

  int fluxout_alk_i;
  int flux_alk_wc_i;
  int flux_alk_sed_i;

} workspace;

void diffusion_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;
    workspace* ws = malloc(sizeof(workspace));
    int i;

    int OFFSET_SED = tracers->n;
    int OFFSET_EPI = tracers->n * 2;

    if ((e->find_index(tracers, "porosity", e) == -1)){
      quit("Porosity tracer not found");
    }

    /*
     * parameters
     */

    double EpiDiffCoeff = try_parameter_value(e,"EpiDiffCoeff");
    if (isnan(EpiDiffCoeff)){
      EpiDiffCoeff = 3.000000e-07;
      eco_write_setup(e,"Code default of EpiDiffCoeff = %e \n",EpiDiffCoeff);
    }

    double EpiDiffDz = try_parameter_value(e,"EpiDiffDz");
    if (isnan(EpiDiffDz)){
      EpiDiffDz = 6.500000e-03;
      eco_write_setup(e,"Code default of EpiDiffDz = %e \n",EpiDiffDz);
    }
    
    // double EpiDiffCoeff = get_parameter_value(e, "EpiDiffCoeff");
    // double EpiDiffDz = get_parameter_value(e, "EpiDiffDz");

    ws->coeff = EpiDiffCoeff / EpiDiffDz;

    /*
     * tracers
     */
    ws->ntr = 0;
    ws->tracer_wc_i = NULL;
    ws->tracer_sed_i = NULL;
    for (i = 0; i < e->ntr; ++i) {
      if (ginterface_gettracerparticflag(e->model, e->tracernames[i]) ==  0 &&
          ginterface_gettracerdiagnflag(e->model, e->tracernames[i]) == 0) {
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

    ws->fluxout_oxygen_i = e->try_index(epis, "Oxygen_sedflux", e); 

    if (ws->fluxout_oxygen_i > -1){
      ws->fluxout_oxygen_i += OFFSET_EPI;
      ws->flux_oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
      ws->flux_oxygen_sed_i = e->find_index(tracers, "Oxygen", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water oxygen flux"); 
    }

    ws->fluxout_nitrate_i = e->try_index(epis, "NO3_sedflux", e);

    if (ws->fluxout_nitrate_i > -1){
      ws->fluxout_nitrate_i += OFFSET_EPI;
      ws->flux_nitrate_wc_i = e->find_index(tracers, "NO3", e);
      ws->flux_nitrate_sed_i = e->find_index(tracers, "NO3", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water nitrate flux"); 
    }

    ws->fluxout_dip_i = e->try_index(epis, "DIP_sedflux", e);

    if (ws->fluxout_dip_i > -1){
      ws->fluxout_dip_i += OFFSET_EPI;
      ws->flux_dip_wc_i = e->find_index(tracers, "DIP", e);
      ws->flux_dip_sed_i = e->find_index(tracers, "DIP", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water DIP flux"); 
    }

    ws->fluxout_amm_i = e->try_index(epis, "NH4_sedflux", e);

    if (ws->fluxout_amm_i > -1){
      ws->fluxout_amm_i += OFFSET_EPI;
      ws->flux_amm_wc_i = e->find_index(tracers, "NH4", e);
      ws->flux_amm_sed_i = e->find_index(tracers, "NH4", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water NH4 flux"); 
    }

    ws->fluxout_dic_i = e->try_index(epis, "DIC_sedflux", e);

    if (ws->fluxout_dic_i > -1){
      ws->fluxout_dic_i += OFFSET_EPI;
      ws->flux_dic_wc_i = e->find_index(tracers, "DIC", e);
      ws->flux_dic_sed_i = e->find_index(tracers, "DIC", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water DIC flux");
    }

    ws->fluxout_alk_i = e->try_index(epis, "alk_sedflux", e);

    if (ws->fluxout_alk_i > -1){
      ws->fluxout_alk_i += OFFSET_EPI;
      ws->flux_alk_wc_i = e->find_index(tracers, "alk", e);
      ws->flux_alk_sed_i = e->find_index(tracers, "alk", e) + OFFSET_SED;
      eco_write_setup(e,"\nOutputting sediment-water alk flux");
    }

    eco_write_setup(e,"\n");

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

void diffusion_epi_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;
  
  if (ws->fluxout_oxygen_i > -1)
    y[ws->fluxout_oxygen_i] = -(y[ws->flux_oxygen_wc_i] - y[ws->flux_oxygen_sed_i]) * ws->coeff;
  
  if (ws->fluxout_nitrate_i > -1)
    y[ws->fluxout_nitrate_i] = -(y[ws->flux_nitrate_wc_i] - y[ws->flux_nitrate_sed_i]) * ws->coeff;
  
  if (ws->fluxout_dip_i > -1)
    y[ws->fluxout_dip_i] = -(y[ws->flux_dip_wc_i] - y[ws->flux_dip_sed_i]) * ws->coeff;

  if (ws->fluxout_dic_i > -1)
    y[ws->fluxout_dic_i] = -(y[ws->flux_dic_wc_i] - y[ws->flux_dic_sed_i]) * ws->coeff;

  if (ws->fluxout_alk_i > -1)
    y[ws->fluxout_alk_i] = -(y[ws->flux_alk_wc_i] - y[ws->flux_alk_sed_i]) * ws->coeff;
 
  if (ws->fluxout_amm_i > -1)
    y[ws->fluxout_amm_i] = -(y[ws->flux_amm_wc_i] - y[ws->flux_amm_sed_i]) * ws->coeff;
}
