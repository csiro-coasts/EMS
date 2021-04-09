/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_mortality_epi.c
 *  
 *  Description: Macroalgae mortality (can be used for spectral / non-spectral model).
 *  
 *  Linear mortality rate.
 *
 *  25/06/2020 MEB - Added multiple seaweed types and shear stress mortality.  
 *  
 *  Options: macroalgae_spectral_grow_epi(b|g|r)
 *
 *  First argument:   b - brown algae (default).
 *                    g - green algae.
 *                    r - red algae.
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: macroalgae_mortality_epi.c 6580 2020-07-29 03:58:50Z bai155 $
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
#include "macroalgae_mortality_epi.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
  
    char species;
    /*
     * parameters
     */
    double mL_t0;
    double MA_tau_critical;
    double MA_tau_efold;

    /*
     * epis
     */
    int MA_N_i;
    int ustrcw_skin_i;

    /*
     * tracers
     */
    int DetBL_N_wc_i;

    /*
     * common cell variables
     */
    int Tfactor_i;
    int mL_i;

  /*
   * Spectral model uses different units for N
   */
  int unitch;
} workspace;

void macroalgae_mortality_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    char* prm = NULL;

    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    if (p->prms->n){
      prm = p->prms->se[0]->s;
      ws->species = prm[0];
    }else{
       ws->species = 'b';
    }

    switch (ws->species){

    case 'r':
      eco_write_setup(e,"Choosen red macroalgae \n");

      ws->mL_t0 = try_parameter_value(e, "MAR_mL");
      if (isnan(ws->mL_t0)){
	ws->mL_t0 = 0.01/86400;
	eco_write_setup(e,"Code default: MAR_mL  = %e \n",ws->mL_t0);
      }

      ws->MA_N_i = e->find_index(epis, "MAR_N", e) + OFFSET_EPI;

      ws->MA_tau_critical = try_parameter_value(e, "MAR_tau_critical");
      if (isnan(ws->MA_tau_critical)){
	ws->MA_tau_critical = 10000000.0;  // no loss
	eco_write_setup(e,"Code default:  MAR_tau_critical = %e \n",ws->MA_tau_critical);
      }

      ws->MA_tau_efold = try_parameter_value(e, "MAR_tau_efold");
      if (isnan(ws->MA_tau_efold)){
	ws->MA_tau_efold = 0.5 * 24.0 * 3600.0;
	eco_write_setup(e,"Code default: MAR_tau_efold = %e \n",ws->MA_tau_efold);
      }
      
      break;

    case 'g':

      eco_write_setup(e,"Choosen green macroalgae \n");

      ws->mL_t0 = try_parameter_value(e, "MAG_mL");
      if (isnan(ws->mL_t0)){
	ws->mL_t0 = 0.01/86400;
	eco_write_setup(e,"Code default: MAG_mL  = %e \n",ws->mL_t0);
      }
     
      ws->MA_N_i = e->find_index(epis, "MAG_N", e) + OFFSET_EPI;

      ws->MA_tau_critical = try_parameter_value(e, "MAG_tau_critical");
      if (isnan(ws->MA_tau_critical)){
	ws->MA_tau_critical = 10000000.0;  // no loss
	eco_write_setup(e,"Code default:  MAG_tau_critical = %e \n",ws->MA_tau_critical);
      }

      ws->MA_tau_efold = try_parameter_value(e, "MAG_tau_efold");
      if (isnan(ws->MA_tau_efold)){
	ws->MA_tau_efold = 0.5 * 24.0 * 3600.0;
	eco_write_setup(e,"Code default: MAG_tau_efold = %e \n",ws->MA_tau_efold);
      }

      break;

    case 'b':

      eco_write_setup(e,"Choosen (or default) brown macroalgae \n");
      ws->mL_t0 = get_parameter_value(e, "MA_mL");
      ws->MA_N_i = e->find_index(epis, "MA_N", e) + OFFSET_EPI;

      ws->MA_tau_critical = try_parameter_value(e, "MA_tau_critical");
      if (isnan(ws->MA_tau_critical)){
	ws->MA_tau_critical = 10000000.0;  // no loss;
	eco_write_setup(e,"Code default:  MA_tau_critical = %e \n",ws->MA_tau_critical);
      }

      ws->MA_tau_efold = try_parameter_value(e, "MA_tau_efold");
      if (isnan(ws->MA_tau_efold)){
	ws->MA_tau_efold = 0.5 * 24.0 * 3600.0;
	eco_write_setup(e,"Code default: MA_tau_efold = %e \n",ws->MA_tau_efold);
      }
      
      break;
    }
    
    /*
     * tracers
     */
    ws->DetBL_N_wc_i = e->find_index(tracers, "DetBL_N", e);
    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->mL_i = find_index_or_add(e->cv_cell, "MA_mL", e);
}

void macroalgae_mortality_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    /*
     * Key off the light spectral process
     */
    ws->unitch = 1;
    if (process_present(e, PT_EPI, "macroalgae_spectral_grow_epi")){
      ws->unitch = 1000.0;
      eco_write_setup(e,"Macroalgae mortality: spectral model so units of macroalgae grams \n");
    }else{
      eco_write_setup(e,"Macroalgae mortality: non-spectral model so units of macroalgae mg \n");
    }
}

void macroalgae_mortality_epi_destroy(eprocess* p)
{
    free(p->workspace);
}

void macroalgae_mortality_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->mL_i] = ws->mL_t0 * Tfactor;
}

void macroalgae_mortality_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;

    double MA_N = y[ws->MA_N_i];
    double mortality = c->cv[ws->mL_i] * MA_N;

    /* add shear stress mortality */

    double tau = y[ws->ustrcw_skin_i]*y[ws->ustrcw_skin_i];
    mortality += min(2.0/86400,max((tau - ws->MA_tau_critical)/ws->MA_tau_critical, 0.0) * (1.0 / ws->MA_tau_efold)) * MA_N ;
    
    y1[ws->MA_N_i] -= mortality;
    y1[ws->DetBL_N_wc_i] += ws->unitch * mortality / c->dz_wc;
}

void macroalgae_mortality_epi_postcalc(eprocess* p, void* pp)
{
}
