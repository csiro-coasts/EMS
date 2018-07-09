/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_root_mortality_epi.c
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
 *  $Id: seagrass_root_mortality_epi.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "seagrass_root_mortality_epi.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  int do_mb;
  char species;

  /*
   * parameters
   */

  double mL_t0;
  double mL_ROOT_t0;
  
  /*
   * epis
   */

  int SGROOT_N_i;
  int SG_N_i;
  
  /*
   * tracers
   */
  int DetBL_N_sed_i;
  
  /*
   * common cell variables
   */
  int Tfactor_i;
  int mL_i;
  int mL_ROOT_i;

} workspace;

void seagrass_root_mortality_epi_init(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = malloc(sizeof(workspace));
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;

  char* prm = p->prms->se[0]->s;
  
  int OFFSET_SED = tracers->n;
  int OFFSET_EPI = tracers->n * 2;
  
  p->workspace = ws;

  ws->species = prm[0];
  
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->DetBL_N_sed_i = e->find_index(tracers, "DetBL_N", e) + OFFSET_SED;
  
  if (ws->species == 'Z'){   /* Zostera - standard terms */
    ws->mL_ROOT_t0 = get_parameter_value(e, "SGROOT_mL");
    ws->mL_t0 = get_parameter_value(e, "SG_mL");
    ws->SG_N_i = e->find_index(epis, "SG_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGROOT_N", e) + OFFSET_EPI;
    ws->mL_i = find_index_or_add(e->cv_cell, "SG_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGROOT_mL", e);

  } else if (ws->species == 'H'){   /* Halophila */
    ws->mL_ROOT_t0 = get_parameter_value(e, "SGHROOT_mL");  
    ws->mL_t0 = get_parameter_value(e, "SGH_mL");
    ws->SG_N_i = e->find_index(epis, "SGH_N", e) + OFFSET_EPI;
    ws->SGROOT_N_i = e->find_index(epis, "SGHROOT_N", e) + OFFSET_EPI;
    ws->mL_i = find_index_or_add(e->cv_cell, "SGH_mL", e);
    ws->mL_ROOT_i = find_index_or_add(e->cv_cell, "SGROOT_mL", e);

  } else if (ws->species == 'P'){   /* Posidonia */
    
  } else{
    e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected seagrass type: expected \"Zostera\", \"Halophila\" or \"Posidonia\"\n", e->processfname, p->name, prm);
  }
}

void seagrass_root_mortality_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void seagrass_root_mortality_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;
    cv[ws->mL_ROOT_i] = ws->mL_ROOT_t0 * Tfactor;

}

void seagrass_root_mortality_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;

    double SG_N = y[ws->SG_N_i];
    double SGROOT_N = y[ws->SGROOT_N_i];
    double mortality = cv[ws->mL_i] * SG_N;
    double mort_root = cv[ws->mL_ROOT_i] * SGROOT_N;

    y1[ws->SGROOT_N_i] -= mort_root;
    y1[ws->SG_N_i] -= mortality;
    y1[ws->DetBL_N_sed_i] += (mortality + mort_root) / c->dz_sed;
}

void seagrass_root_mortality_epi_postcalc(eprocess* p, void* pp)
{
}
