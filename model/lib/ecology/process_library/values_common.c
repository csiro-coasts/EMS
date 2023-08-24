/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/values_common.c
 *  
 *  Description:
 *  
 *  Process that generates diagnostic variables (DIN,EFI, EPO) that are used in 
 *  ecology calculations, and recalculated in postcalc for output.  
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: values_common.c 7351 2023-04-30 03:42:14Z bai155 $
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
#include "values_common.h"

typedef struct {
  int with_df;                /* flag */
  int with_mpb;               /* flag */
  int with_Tricho;            /* flag */

    /*
     * tracers
     */
  int DIN_i;
  int NO3_i;
  int NH4_i;
  int Kd_i;
  int salt_i;
  int EFI_i;
  int Mud_i;
  int FineSed_i;
  int Dust_i;
  
  int Mud_mineral_i;
  int Mud_carbonate_i;

  int PhyS_N_i;
  int PhyL_N_i;
  int PhyD_N_i;
  int Tricho_N_i;
  int MPB_N_i;
  int ZooS_N_i;
  int ZooL_N_i;
  int DetPL_N_i;
  int DetBL_N_i;
  int DetR_C_i;
  int CS_N_free_i;
  int SMA_N_i;

  int EPO_i;

  int do_light_spectral_col;

  int PhyL_yCfac_i;
  int PhyL_kI_i;
  int PhyS_yCfac_i;
  int PhyS_kI_i;
  int PhyD_yCfac_i;
  int PhyD_kI_i;
  int Tricho_yCfac_i;
  int Tricho_kI_i;
  int MPB_yCfac_i;
  int MPB_kI_i;
  int SMA_kI_i;

  int cv_PhyL_yCfac_i;
  int cv_PhyL_kI_i;
  int cv_PhyS_yCfac_i;
  int cv_PhyS_kI_i;
  int cv_PhyD_yCfac_i;
  int cv_PhyD_kI_i;
  int cv_Tricho_yCfac_i;
  int cv_Tricho_kI_i;
  int cv_MPB_yCfac_i;
  int cv_MPB_kI_i;
  int cv_SMA_kI_i;

  int CS_free_yCfac_i;
  int CS_free_kI_i;
  int cv_CS_free_yCfac_i;
  int cv_CS_free_kI_i;

  /*
   * parameters
   */
  double k_w;
  double k_C_fw;
  double k_TSS;
  
} workspace;

void values_common_init(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = malloc(sizeof(workspace));
  stringtable* tracers = e->tracers;
  
  memset(ws, 0, sizeof(workspace));
  
  p->workspace = ws;
  
  /*
   * tracers
   */

  ws->DIN_i = -1;
  ws->NO3_i = -1;
  ws->NH4_i = -1;
  ws->salt_i = -1;
  ws->Mud_i  = -1;
  ws->FineSed_i = -1;
  ws->EFI_i = -1;
  ws->Dust_i = -1;
  ws->Mud_mineral_i = -1;
  ws->Mud_carbonate_i = -1;

  ws->DIN_i = e->find_index(e->tracers, "DIN", e);
  ws->NO3_i = e->find_index(e->tracers, "NO3", e);
  ws->NH4_i = e->find_index(e->tracers, "NH4", e);
  ws->salt_i = e->find_index(e->tracers, "salt", e);
  ws->Mud_i  = e->try_index(e->tracers, "Mud", e);
  ws->FineSed_i = e->find_index(e->tracers, "FineSed", e);
  ws->EFI_i = e->find_index(e->tracers, "EFI", e);
  ws->Dust_i = e->try_index(e->tracers, "Dust", e);

  ws->Mud_mineral_i = e->try_index(e->tracers, "Mud-mineral", e);
  ws->Mud_carbonate_i = e->try_index(e->tracers, "Mud-carbonate", e);

  ws->Kd_i = -1;
  ws->PhyS_N_i = -1;
  ws->PhyL_N_i = -1;
  ws->PhyD_N_i = -1;
  ws->Tricho_N_i = -1;
  ws->MPB_N_i = -1;
  ws->SMA_N_i = -1;
  ws->EPO_i = -1;

  ws->PhyS_N_i = e->find_index(e->tracers, "PhyS_N", e);
  ws->PhyL_N_i = e->find_index(e->tracers, "PhyL_N", e);
  ws->Tricho_N_i = e->try_index(e->tracers, "Tricho_N", e);
  ws->MPB_N_i = e->find_index(e->tracers, "MPB_N", e);
  ws->PhyD_N_i = e->try_index(e->tracers, "PhyD_N", e);

  ws->SMA_N_i = e->try_index(e->tracers, "SMA_N", e);

  ws->EPO_i = e->try_index(e->tracers, "EPO", e);

  ws->CS_N_free_i = -1;
  ws->CS_N_free_i = e->try_index(e->tracers, "CS_N_free", e);

  if (ws->EPO_i > -1){

    ws->ZooS_N_i = e->find_index(e->tracers, "ZooS_N", e);
    ws->ZooL_N_i = e->find_index(e->tracers, "ZooL_N", e);
    
    ws->DetPL_N_i = e->find_index(e->tracers, "DetPL_N", e);
    ws->DetBL_N_i = e->find_index(e->tracers, "DetBL_N", e);
    ws->DetR_C_i = e->find_index(e->tracers, "DetR_C", e);

 /*non essential diagnostic tracer */

    if (ws->PhyD_N_i > -1)
      ws->with_df = 1;
    else
      ws->with_df = 0;
    
    if (ws->Tricho_N_i > -1)
      ws->with_Tricho = 1;
    else
      ws->with_Tricho = 0;  
  }
  
  /* Output particles used in EFI */

  eco_write_setup(e,"\nParticles used in EFI: ");

  if (ws->FineSed_i > -1)
    eco_write_setup(e,"FineSed ");
  
  if (ws->Mud_i > -1)
    eco_write_setup(e,"Mud ");

  if (ws->Mud_carbonate_i > -1)
    eco_write_setup(e,"Mud_carbonate ");

  if (ws->Mud_mineral_i > -1)
    eco_write_setup(e,"Mud_mineral ");
  
  if (ws->Dust_i > -1)
    eco_write_setup(e,"Dust ");

  eco_write_setup(e,"\n");

  /* Output particles used in EFO */

  if (ws->EPO_i > -1){

    eco_write_setup(e,"\nOrganic particles used in EFO: PhyS_N PhyL_N MPB_N ZooS_N ZooL_N DetPL_N DetBL_N DetR_C ");

    if (ws->PhyD_N_i > -1)
      eco_write_setup(e,"PhyD_N ");

    if (ws->Tricho_N_i > -1)
      eco_write_setup(e,"Tricho_N  ");

    eco_write_setup(e,"\n");
  }
}

void values_common_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace *ws = p->workspace;
  stringtable* tracers = e->tracers;

  /*
   * parameters
   */
  if (process_present(e,PT_WC,"light_wc")) {
    ws->k_w = get_parameter_value(e, "k_w");
    ws->k_C_fw = get_parameter_value(e, "k_C_fw");
    ws->k_TSS = get_parameter_value(e, "k_TSS");
    ws->Kd_i = e->find_index(e->tracers, "Kd", e);
    eco_write_setup(e,"\nCalculating Kd for old non-spectrally-resolving model.");
  }

  ws->do_light_spectral_col = 0;

  ws->cv_PhyL_yCfac_i = -1;
  ws->cv_PhyL_kI_i = -1;
  ws->cv_PhyS_yCfac_i = -1;
  ws->cv_PhyS_kI_i = -1;
  ws->cv_PhyD_yCfac_i = -1;
  ws->cv_PhyD_kI_i = -1;
  ws->cv_Tricho_yCfac_i = -1;
  ws->cv_Tricho_kI_i = -1;
  ws->cv_MPB_yCfac_i = -1;
  ws->cv_MPB_kI_i = -1;
  ws->cv_SMA_kI_i = -1;
  ws->PhyL_yCfac_i = -1;
  ws->PhyL_kI_i = -1;
  ws->PhyS_yCfac_i = -1;
  ws->PhyS_kI_i = -1;
  ws->PhyD_yCfac_i = -1;
  ws->PhyD_kI_i = -1;
  ws->Tricho_yCfac_i = -1;
  ws->Tricho_kI_i = -1;
  ws->MPB_yCfac_i = -1;
  ws->MPB_kI_i = -1;

  ws->cv_CS_free_yCfac_i = -1;
  ws->cv_CS_free_kI_i = -1;

  ws->CS_free_yCfac_i = -1;
  ws->CS_free_kI_i = -1;

  if (process_present(e,PT_COL,"light_spectral_col")){
      ws->do_light_spectral_col = 1;
      if (ws->PhyL_N_i > -1){
	ws->cv_PhyL_kI_i = find_index_or_add(e->cv_cell, "KI_l", e);
	ws->cv_PhyL_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_l", e);
	ws->PhyL_kI_i = e->find_index(tracers, "PhyL_kI", e);
	ws->PhyL_yCfac_i = e->find_index(tracers, "PhyL_yCfac", e);
      }
      if (ws->PhyS_N_i >-1){
	ws->cv_PhyS_kI_i = find_index_or_add(e->cv_cell, "KI_s", e);
	ws->cv_PhyS_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_s", e);
	ws->PhyS_kI_i = e->find_index(tracers, "PhyS_kI", e);
	ws->PhyS_yCfac_i = e->find_index(tracers, "PhyS_yCfac", e);
      }
      if (ws->MPB_N_i >-1){
	ws->cv_MPB_kI_i = find_index_or_add(e->cv_cell, "KI_MPB", e);
	ws->cv_MPB_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_MPB", e);
	ws->MPB_kI_i = e->find_index(tracers, "MPB_kI", e);
	ws->MPB_yCfac_i = e->find_index(tracers, "MPB_yCfac", e);
      }
      if (ws->Tricho_N_i >-1){
	ws->cv_Tricho_kI_i = find_index_or_add(e->cv_cell, "KI_Tricho", e);
	ws->cv_Tricho_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_Tricho", e);
	ws->Tricho_kI_i = e->find_index(tracers, "Tricho_kI", e);
	ws->Tricho_yCfac_i = e->find_index(tracers, "Tricho_yCfac", e);
      }
      if (ws->PhyD_N_i >-1){
	ws->cv_PhyD_kI_i = find_index_or_add(e->cv_cell, "KI_PhyD", e);
	ws->cv_PhyD_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_PhyD", e);
	ws->PhyD_kI_i = e->find_index(tracers, "PhyD_kI", e);
	ws->PhyD_yCfac_i = e->find_index(tracers, "PhyD_yCfac", e);
      }
      if (ws->CS_N_free_i > -1){
	ws->cv_CS_free_kI_i = find_index_or_add(e->cv_cell, "KI_CS_free", e);
	ws->cv_CS_free_yCfac_i = find_index_or_add(e->cv_cell, "yCfac_CS_free", e);
	ws->CS_free_kI_i = e->find_index(tracers, "CS_free_kI", e);
	ws->CS_free_yCfac_i = e->find_index(tracers, "CS_free_yCfac", e);
      }
      if (ws->SMA_N_i > -1){
	ws->cv_SMA_kI_i = find_index_or_add(e->cv_cell, "KI_SMA", e);
	ws->SMA_kI_i = e->find_index(tracers, "SMA_kI", e);
      }
    }  
}

void values_common_destroy(eprocess* p)
{
    free(p->workspace);
}

void values_common_precalc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  cell* c = (cell*) pp;
  double* y = c->y;
  double salt = (y[ws->salt_i] > 35.0) ? 35.0 : y[ws->salt_i];
  double EFI = 0.0;

  // NOTE: kI from light_col multiplied by 1e20 to avoid reset for small numbers by sediment model, need to
  //       divide by that amount here,
  
  if (ws->do_light_spectral_col == 1){

    if (c->type == CT_WC){ // Water column

      if (ws->PhyL_kI_i > -1){
	c->cv[ws->cv_PhyL_kI_i] = y[ws->PhyL_kI_i] * 8.359335857479461e-29 ;
	c->cv[ws->cv_PhyL_yCfac_i] = y[ws->PhyL_yCfac_i];
      }
      if (ws->PhyS_kI_i > -1){
	c->cv[ws->cv_PhyS_kI_i] = y[ws->PhyS_kI_i] * 8.359335857479461e-29 ;
	c->cv[ws->cv_PhyS_yCfac_i] = y[ws->PhyS_yCfac_i];
      }
      if (ws->MPB_kI_i > -1){
	c->cv[ws->cv_MPB_kI_i] = y[ws->MPB_kI_i] * 8.359335857479461e-29 ;
	c->cv[ws->cv_MPB_yCfac_i] = y[ws->MPB_yCfac_i];
      }
      
      if (ws->Tricho_kI_i > -1){
	c->cv[ws->cv_Tricho_kI_i] = y[ws->Tricho_kI_i] * 8.359335857479461e-29 ;
	c->cv[ws->cv_Tricho_yCfac_i] = y[ws->Tricho_yCfac_i];
      }
      if (ws->PhyD_kI_i > -1){
	c->cv[ws->cv_PhyD_kI_i] = y[ws->PhyD_kI_i] * 8.359335857479461e-29 ;
	c->cv[ws->cv_PhyD_yCfac_i] = y[ws->PhyD_yCfac_i];
      }
      if (ws->CS_N_free_i > -1){
	c->cv[ws->cv_CS_free_kI_i] = y[ws->CS_free_kI_i] * 8.359335857479461e-29;
	c->cv[ws->cv_CS_free_yCfac_i] = y[ws->CS_free_yCfac_i];
      }
      if (ws->SMA_N_i > -1){
	c->cv[ws->cv_SMA_kI_i] = y[ws->SMA_kI_i];
      }
    }

    if (c->type == CT_SED){ // Sediments - use same index as WC because different y array passed.

      if (ws->cv_MPB_kI_i > -1){
	
	c->cv[ws->cv_MPB_kI_i] = 0.0;
	c->cv[ws->cv_MPB_yCfac_i] = 0.0;
	
	if (c->k_sed == c->col->topk_sed){  // top sediment layer
	  c->cv[ws->cv_MPB_kI_i] = y[ws->MPB_kI_i] * 8.359335857479461e-29 ;
	  c->cv[ws->cv_MPB_yCfac_i] = y[ws->MPB_yCfac_i];
	} 
      }

      if (ws->CS_N_free_i > -1){
	
	c->cv[ws->cv_CS_free_kI_i] = 0.0;
	c->cv[ws->cv_CS_free_yCfac_i] = 0.0;
	
	if (c->k_sed == c->col->topk_sed){
	  c->cv[ws->cv_CS_free_kI_i] = y[ws->CS_free_kI_i] * 8.359335857479461e-29 ;
	  c->cv[ws->cv_CS_free_yCfac_i] = y[ws->CS_free_yCfac_i];
	}
      }
    }
  }

  if (ws->FineSed_i > -1){
    EFI += y[ws->FineSed_i];
  }
  
  if (ws->Mud_i > -1){
    EFI += y[ws->Mud_i];
  }

  if (ws->Mud_carbonate_i > -1){
    EFI += y[ws->Mud_carbonate_i];
  }

  if (ws->Mud_mineral_i > -1){
    EFI += y[ws->Mud_mineral_i];
  }
  
  if (ws->Dust_i > -1){
    EFI += y[ws->Dust_i];
  }

  if (ws->Kd_i  > -1){
    if (c->type == CT_WC){
      y[ws->Kd_i] += ws->k_TSS *EFI+ ws->k_w + ws->k_C_fw * (35.0 - salt) / 35.0;
    }
  }
  y[ws->EFI_i] = EFI;
}

void values_common_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    double* y = ((cell*) pp)->y;

    if (ws->EPO_i > -1){
      y[ws->EPO_i] = ((y[ws->PhyS_N_i] + y[ws->PhyL_N_i] + y[ws->MPB_N_i] + y[ws->ZooS_N_i] + y[ws->ZooL_N_i] + y[ws->DetPL_N_i]) * red_W_C + y[ws->DetBL_N_i] * atk_W_C + y[ws->DetR_C_i]) / 1.0e6;  // kg C m-3

      if (ws->PhyD_N_i > -1){
	y[ws->EPO_i] +=  y[ws->PhyD_N_i] * (106.0*12.01/16.0/14.01)/1.0e6;
      }
      if (ws->Tricho_N_i > -1){
	y[ws->EPO_i] += y[ws->Tricho_N_i] * (106.0*12.01/16.0/14.01)/1.0e6;
      }
    }
    
    y[ws->DIN_i] = y[ws->NO3_i] + y[ws->NH4_i];

    y[ws->EFI_i] = y[ws->FineSed_i];

    if (ws->Mud_i > -1){
      y[ws->EFI_i] += y[ws->Mud_i];
    }
    
    if (ws->Mud_carbonate_i > -1){
      y[ws->EFI_i] += y[ws->Mud_carbonate_i];
    }
    
    if (ws->Mud_mineral_i > -1){
      y[ws->EFI_i] += y[ws->Mud_mineral_i];
    }
    
    if (ws->Dust_i > -1){
       y[ws->EFI_i] += y[ws->Dust_i];
  }
}
