/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/diffusion_heat_epi.c
 *  
 *  Description:
 *  This procedure accounts for the diffusion heat and salt
 *  from the bottom water layer to the sediment.
 *
 *  Heat conduction is not presently considered, and we assume the specific heat of 
 *  water and particles are equal.
 *
 *  Heat and salinity do not move from the sediment to the water column, so the whole system 
 *  is not energy conserving - it is just a means to keep sediment temperature and salt close to 
 *  the hydrodynamic prediction of bottom conditions.
 *
 *  For more information see diffusion_epi.c which applies sediment - water fluxes for 
 *  all other dissolved tracers.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: diffusion_heat_epi.c 6551 2020-05-10 23:15:52Z bai155 $
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
#include "diffusion_heat_epi.h"

typedef struct {
  /*
   * coeff = EpiDiffCoeff / EpiDiffDz
   */
  double coeff;
  
  /*
   * tracers
   */

  int temp_wc_i;
  int temp_sed_i;
  int salt_wc_i;
  int salt_sed_i;

} workspace;

void diffusion_heat_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    int i;

    int OFFSET_SED = tracers->n;

    ws->temp_wc_i = e->find_index(tracers, "temp", e);
    ws->temp_sed_i = e->find_index(tracers, "temp", e) + OFFSET_SED;

    ws->salt_wc_i = e->find_index(tracers, "salt", e);
    ws->salt_sed_i = e->find_index(tracers, "salt", e) + OFFSET_SED;

    if ((e->find_index(tracers, "porosity", e) == -1)){
      quit("Porosity tracer not found");
    }

    /*
     * parameters
     */

    eco_write_setup(e,"\nHeat transfer changing only sediment temperature \n");

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

    ws->coeff = EpiDiffCoeff / EpiDiffDz;

    p->workspace = ws;
}

void diffusion_heat_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void diffusion_heat_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double porosity = c->porosity;

    double heatflux = -(y[ws->temp_wc_i] - y[ws->temp_sed_i]) * ws->coeff;

    double saltflux = -(y[ws->salt_wc_i] - y[ws->salt_sed_i]) * ws->coeff;

    // only update sediment

    // don't divide by porosity because we are heating the water and particles together?
    
    y1[ws->temp_sed_i] -= heatflux / c->dz_sed; // porosity;

    y1[ws->salt_sed_i] -= saltflux / c->dz_sed / porosity;
}

void diffusion_heat_epi_postcalc(eprocess* p, void* pp)
{
}
