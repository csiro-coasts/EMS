/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/values_common.c
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
 *  $Id: values_common.c 5846 2018-06-29 04:14:26Z riz008 $
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

  int EPO_i;
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

  ws->EPO_i = e->try_index(e->tracers, "EPO", e);

  if (ws->EPO_i > -1){

    ws->PhyS_N_i = e->find_index(e->tracers, "PhyS_N", e);
    ws->PhyL_N_i = e->find_index(e->tracers, "PhyL_N", e);
    ws->PhyD_N_i = e->try_index(e->tracers, "PhyD_N", e);
    ws->Tricho_N_i = e->try_index(e->tracers, "Tricho_N", e);
    ws->MPB_N_i = e->find_index(e->tracers, "MPB_N", e);
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

}

void values_common_post_init(eprocess* p)
{
  ecology* e = p->ecology;
  workspace *ws = p->workspace;

  /*
   * parameters
   */
  if (process_present(e,PT_WC,"light_wc")) {
    ws->k_w = get_parameter_value(e, "k_w");
    ws->k_C_fw = get_parameter_value(e, "k_C_fw");
    ws->k_TSS = get_parameter_value(e, "k_TSS");
    ws->Kd_i = e->find_index(e->tracers, "Kd", e);   
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
      /*y[ws->EPO_i] = ((y[ws->PhyS_N_i] + y[ws->PhyL_N_i] + y[ws->MPB_N_i] + y[ws->Tricho_N_i] + y[ws->ZooS_N_i] + y[ws->ZooL_N_i] + y[ws->DetPL_N_i]) * (106.0*12.01/16.0/14.01) + y[ws->DetBL_N_i] * (550.0*12.01/30.0/14.01) + y[ws->DetR_C_i])/1.0e6;*/
      y[ws->EPO_i] = ((y[ws->PhyS_N_i] + y[ws->PhyL_N_i] + y[ws->MPB_N_i] + y[ws->ZooS_N_i] + y[ws->ZooL_N_i] + y[ws->DetPL_N_i]) * (106.0*12.01/16.0/14.01) + y[ws->DetBL_N_i] * (550.0*12.01/30.0/14.01) + y[ws->DetR_C_i])/1.0e6;

      if (ws->PhyD_N_i > -1){
      y[ws->EPO_i] +=  y[ws->PhyD_N_i] * (106.0*12.01/16.0/14.01);
      }
      if (ws->Tricho_N_i > -1){
      y[ws->EPO_i] += y[ws->Tricho_N_i] * (106.0*12.01/16.0/14.01);
      }}

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
