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
 *  $Id: values_common_epi.c 7343 2023-04-11 22:49:32Z bai155 $
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
  int Oxy_pr_wc_i;
  int Oxy_pr_sed_i;
  int EpiOxy_pr_i;
  int ustrcw_skin_i;

  int MA_N_i;
  int MAR_N_i; 
  int MAG_N_i;
  
  int SG_N_i;
  int SGH_N_i;
  int SGD_N_i;
  int SGP_N_i;

  int CS_N_i;

  int MA_kI_i;
  int MAG_kI_i;
  int MAR_kI_i;
  
  int SG_kI_i;
  int SGH_kI_i;
  int SGD_kI_i;
  int SGP_kI_i;

  int cv_MA_kI_i;
  int cv_MAG_kI_i;
  int cv_MAR_kI_i;

  int cv_SG_kI_i;
  int cv_SGH_kI_i;
  int cv_SGD_kI_i;
  int cv_SGP_kI_i;

  int cv_CS_kI_i;
  int cv_CS_yCfac_i;
  int CS_kI_i;
  int CS_yCfac_i;

  int do_light_spectral_col;
  
  int do_mixed_layer_age_col;
  int cv_mixed_layer_depth_i;
  int mixed_layer_depth_i;
  int ecol_mld_i;

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
    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

    ws->Oxy_pr_sed_i = e->try_index(tracers, "Oxy_pr", e);
    if (ws->Oxy_pr_sed_i > -1)
      ws->Oxy_pr_sed_i += OFFSET_SED;
    
    ws->EpiOxy_pr_i = e->try_index(epis, "EpiOxy_pr", e);
    if (ws->EpiOxy_pr_i > -1)
      ws->EpiOxy_pr_i += OFFSET_EPI;
}

void values_common_epi_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace *ws = p->workspace;
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;

  int OFFSET_SED = tracers->n;
  int OFFSET_EPI = tracers->n * 2;

    ws->do_light_spectral_col = 0;
    ws->do_mixed_layer_age_col = 0;

    ws->MA_N_i = -1;
    ws->MAR_N_i = -1; 
    ws->MAG_N_i = -1;
    
    ws->SG_N_i= -1;
    ws->SGH_N_i = -1;
    ws->SGD_N_i = -1;
    ws->SGP_N_i = -1;

    ws->MA_N_i = e->try_index(epis, "MA_N", e);
    if (ws->MA_N_i > -1)
      ws->MA_N_i += OFFSET_EPI;

    ws->MAG_N_i = e->try_index(epis, "MAG_N", e);
    if (ws->MAG_N_i > -1)
      ws->MAG_N_i += OFFSET_EPI;

    ws->MAR_N_i = e->try_index(epis, "MAR_N", e);
    if (ws->MAR_N_i > -1)
      ws->MAR_N_i += OFFSET_EPI;

    ws->SG_N_i = e->try_index(epis, "SG_N", e);
    if (ws->SG_N_i > -1)
      ws->SG_N_i += OFFSET_EPI;

    ws->SGH_N_i = e->try_index(epis, "SGH_N", e);
    if (ws->SGH_N_i > -1)
      ws->SGH_N_i += OFFSET_EPI;

    ws->SGD_N_i = e->try_index(epis, "SGD_N", e);
    if (ws->SGD_N_i > -1)
      ws->SGD_N_i += OFFSET_EPI;
    
    ws->SGP_N_i = e->try_index(epis, "SGP_N", e);
    if (ws->SGP_N_i > -1)
      ws->SGP_N_i += OFFSET_EPI;

    ws->CS_N_i = e->try_index(epis, "CS_N", e);
    if (ws->CS_N_i > -1)
      ws->CS_N_i += OFFSET_EPI;

    ws->MA_kI_i = -1;
    ws->MAR_kI_i = -1; 
    ws->MAG_kI_i = -1;
    
    ws->SG_kI_i= -1;
    ws->SGH_kI_i = -1;
    ws->SGD_kI_i = -1;
    ws->SGP_kI_i = -1;

    ws->CS_kI_i = -1;
    ws->CS_yCfac_i = -1;
    
    ws->cv_MA_kI_i = -1;
    ws->cv_MAR_kI_i = -1;
    ws->cv_MAG_kI_i = -1;

    ws->cv_SG_kI_i = -1;
    ws->cv_SGH_kI_i = -1;
    ws->cv_SGD_kI_i = -1;
    ws->cv_SGP_kI_i = -1;

    ws->cv_CS_kI_i = -1;
    ws->cv_CS_yCfac_i = -1;

    if (process_present(e,PT_COL,"light_spectral_col")){
      ws->do_light_spectral_col = 1;
      if (ws->SG_N_i > -1){
	ws->cv_SG_kI_i = find_index_or_add(e->cv_cell, "KI_SG", e);
	ws->SG_kI_i = e->find_index(epis, "SG_kI", e) + OFFSET_EPI;
      }
      if (ws->SGH_N_i > -1){
	ws->cv_SGH_kI_i = find_index_or_add(e->cv_cell, "KI_SGH", e);
	ws->SGH_kI_i = e->find_index(epis, "SGH_kI", e) + OFFSET_EPI;
      }
      if (ws->SGD_N_i > -1){
	ws->cv_SGD_kI_i = find_index_or_add(e->cv_cell, "KI_SGD", e);
	ws->SGD_kI_i = e->find_index(epis, "SGD_kI", e) + OFFSET_EPI;
      }
      if (ws->SGP_N_i > -1){
	ws->cv_SGP_kI_i = find_index_or_add(e->cv_cell, "KI_SGP", e);
	ws->SGP_kI_i = e->find_index(epis, "SGP_kI", e) + OFFSET_EPI;
      }
      if (ws->MA_N_i > -1){
	ws->cv_MA_kI_i = find_index_or_add(e->cv_cell, "KI_MA", e);
	ws->MA_kI_i = e->find_index(epis, "MA_kI", e) + OFFSET_EPI;
      }
      if (ws->MAR_N_i > -1){
	ws->cv_MAR_kI_i = find_index_or_add(e->cv_cell, "KI_MAR", e);
	ws->MAR_kI_i = e->find_index(epis, "MAR_kI", e) + OFFSET_EPI;
      }
      if (ws->MAG_N_i > -1){
	ws->cv_MAG_kI_i = find_index_or_add(e->cv_cell, "KI_MAG", e);
	ws->MAG_kI_i = e->find_index(epis, "MAG_kI", e) + OFFSET_EPI;
      }
      if (ws->CS_N_i > -1){
	ws->cv_CS_kI_i = find_index_or_add(e->cv_cell, "KI_CS", e);
	ws->cv_CS_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_CS", e);
	ws->CS_yCfac_i = e->find_index(epis, "CS_yCfac", e) + OFFSET_EPI;
	ws->CS_kI_i = e->find_index(epis, "CS_kI", e) + OFFSET_EPI;
      }
    }
    ws->mixed_layer_depth_i = -1;
    ws->ecol_mld_i = -1;
    
    if (process_present(e,PT_COL,"mixed_layer_age_col")){
      ws->mixed_layer_depth_i = e->try_index(epis, "mixed_layer", e);
      if (ws->mixed_layer_depth_i > -1){
	ws->mixed_layer_depth_i += OFFSET_EPI;
	ws->do_mixed_layer_age_col = 1;
	// ws->cv_mixed_layer_depth_i = find_index_or_add(e->cv_cell, "cv_mixed_layer_depth", e);
      }
      ws->ecol_mld_i = e->try_index(epis, "ecology_mixed_layer", e);
      if (ws->ecol_mld_i > -1)
	ws->ecol_mld_i += OFFSET_EPI;
    }
    // printf("ws->do_mixed_layer_age_col %d ws->mixed_layer_depth_i %d ws->cv_mixed_layer_depth_i %d \n",ws->do_mixed_layer_age_col,ws->mixed_layer_depth_i,ws->cv_mixed_layer_depth_i);
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
 
    if (ws->do_light_spectral_col == 1){
      if (ws->MA_kI_i > -1)
	c->cv[ws->cv_MA_kI_i] = y[ws->MA_kI_i];
      if (ws->MAR_kI_i > -1)
	c->cv[ws->cv_MAR_kI_i] = y[ws->MAR_kI_i];
      if (ws->MAG_kI_i > -1)
	c->cv[ws->cv_MAG_kI_i] = y[ws->MAG_kI_i];
      if (ws->SG_kI_i > -1)
	c->cv[ws->cv_SG_kI_i] =  y[ws->SG_kI_i];
      if (ws->SGH_kI_i > -1)
	c->cv[ws->cv_SGH_kI_i] = y[ws->SGH_kI_i];
      if (ws->SGD_kI_i > -1)
	c->cv[ws->cv_SGD_kI_i] = y[ws->SGD_kI_i];
      if (ws->SGP_kI_i > -1)
	c->cv[ws->cv_SGP_kI_i] = y[ws->SGP_kI_i];
      if (ws->CS_kI_i > -1){
	c->cv[ws->cv_CS_kI_i] = y[ws->CS_kI_i] * 8.359335857479461e-29;	
	c->cv[ws->cv_CS_yCfac_i] = y[ws->CS_yCfac_i];
      }
    }
    if (ws->do_mixed_layer_age_col == 1){
      if (ws->mixed_layer_depth_i > -1){
	// c->cv[ws->cv_mixed_layer_depth_i] = y[ws->mixed_layer_depth_i];
	y[ws->ecol_mld_i] = y[ws->mixed_layer_depth_i];
      }
    }
    // printf("in values_common_epi mld %e, c->cv[ws->cv_mixed_layer_depth_i] %e \n",y[ws->mixed_layer_depth_i],c->cv[ws->cv_mixed_layer_depth_i]);
}

void values_common_epi_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* y = c->y;

    if ( (ws->EpiOxy_pr_i > -1) && 
	 (ws->Oxy_pr_wc_i > -1) &&
	 (ws->Oxy_pr_sed_i > -1) )
      y[ws->EpiOxy_pr_i] = y[ws->Oxy_pr_wc_i] * c->dz_wc + y[ws->Oxy_pr_sed_i];
    /*
     * (already multipled by dz_sed)
     */
}
