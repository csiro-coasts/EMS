/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_sed.c
 *  
 *  Description:
 *  
 *  Sediment spectrally-resolved optical model. Model calculates:
 *  
 *  - scalar irradiance and absorption for MPB (output time - ECOLOGY_DT).
 *
 *  MPB assumed to sit as single layer above the sediments and below the water column.
 *
 *  WARNING: variables output in precalc are calculated at output time - ECOLOGY_DT, but are stored in
 *           output files at output time.
 *
 *  Options: light_spectral_wc(G|H)
 *  
 *  G - Gaussian approx. of chl-a and non chl-a specific absorption.
 *  H - HPLC-determined of pigment-specific absorption.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: light_spectral_sed.c 6238 2019-05-29 06:34:36Z bai155 $
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
#include "einterface.h"
#include "light_spectral_sed.h"
#include "constants.h"

typedef struct {
  /*
   * parameters
   */

  char pig;
  
  double m_MPB;
  double rad_MPB;
  double xan2chl_MPB;

  /*
   * tracers
   */
  int PAR_i;
  int PAR_z_i;
  int wPAR_bot;
  int wPAR_top;
  int MPB_N_i;
  int MPB_Chl_i;
  /*
   * common variables
   */
  int cv_lighttop_s_i;
  int KI_MPB_i;
  int yCfac_MPB_i;

  int num_waves;
  double *wave;
  
  double *yC_MPB;

} workspace;

void light_spectral_sed_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  workspace* ws = malloc(sizeof(workspace));
  


  char* prm = NULL;

  if (p->prms->n){
    prm = p->prms->se[0]->s;
    ws->pig = prm[0];
  }else{
    ws->pig = 'G';   // default is Guassian.
  }

  p->workspace = ws;
  int w;
  
  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");

  ws->yC_MPB = d_alloc_1d(ws->num_waves);
  
  for (w=0; w<ws->num_waves; w++){
    if (ws->wave[w] == 410.0)
      ws->wPAR_bot = w;
    if (ws->wave[w] == 690.0)
      ws->wPAR_top = w;
  }
  
  ws->MPB_N_i   = e->try_index(tracers, "MPB_N", e);
  ws->MPB_Chl_i = e->try_index(tracers, "MPB_Chl", e);

  /*
   * parameters
   */
  if (ws->MPB_N_i > -1){
    ws->rad_MPB = get_parameter_value(e, "MBrad");
    ws->m_MPB = PhyCellMass(ws->rad_MPB);
    ws->xan2chl_MPB = try_parameter_value(e, "MBxan2chl");
    if (isnan(ws->xan2chl_MPB))
      ws->xan2chl_MPB = 0.0;
  }
  /*
   * tracers
   */
  // ws->Kd_490_i = e->find_index(tracers, "Kd_490", e);
  
  ws->PAR_i = e->try_index(tracers, "PAR", e);

  ws->PAR_z_i = -1;
  ws->PAR_z_i = e->try_index(tracers, "Light_down", e);
  
  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");
  
  /*
   * common column variables
   */
  ws->cv_lighttop_s_i = find_index_or_add_vector(e->cv_column, "lighttop_s", e, ws->num_waves);
  ws->yCfac_MPB_i = find_index_or_add(e->cv_cell, "yCfac_MPB", e);
  ws->KI_MPB_i = find_index_or_add(e->cv_cell, "KI_MPB", e);

  if (ws->pig == 'H'){  // HPLC determined absorption coefficients

    // do in postinit

  }else{   // Guassian curve determined absorption coefficients 
    
    eco_write_setup(e,"\n Guassian curve determined absorption coefficients in sediments \n");
    
    double *yC_chl;
    double *yC_xanth;
    
    yC_chl = d_alloc_1d(ws->num_waves);
    yC_xanth = d_alloc_1d(ws->num_waves);

    absorbwave(yC_chl,1.0, ws->num_waves, ws->wave);
    absorbxanth(yC_xanth,1.0, ws->num_waves, ws->wave);
    
    for (w=0; w<ws->num_waves; w++){
      ws->yC_MPB[w] = yC_chl[w] + ws->xan2chl_MPB * yC_xanth[w];
    }

    d_free_1d(yC_chl);
    d_free_1d(yC_xanth);
  }
}

void light_spectral_sed_postinit(eprocess* p)
{
  ecology* e = p->ecology;
  workspace* ws = (workspace *)p->workspace;
  
  if (e->pre_build) return;
  
  if (ws->pig == 'H'){  // HPLC determined absorption coefficients

    eco_write_setup(e,"\nMass-specific absorption coefficients used for microalgae in sediments \n");
    
    if (e->bio_opt==NULL){
      ecology_find_rsr_waves(e);
      e->bio_opt = bio_opt_init(e);
    }

    bio_opt_prop *bio = e->bio_opt;
    
    int w;
    for (w=0; w<ws->num_waves; w++){
      ws->yC_MPB[w] = bio->yC_microplankton[w];
    }
  }
}


void light_spectral_sed_destroy(eprocess* p)
{
  free(p->workspace);
}

void light_spectral_sed_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  double lighttop;
  bio_opt_prop *bio = e->bio_opt;
  
  int w;
  double energy = 0.0;
  
  double* lighttop_s = col->cv[ws->cv_lighttop_s_i];
  double Kd_s[ws->num_waves];
  double light_s[ws->num_waves];
  
  for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
    energy += lighttop_s[w];
  }
  
  if (isnan(lighttop_s[0]) || energy < 0.01){
    
    for (w=0; w<ws->num_waves; w++) {
      lighttop_s[w] =  0.0;
    }
    if (ws->KI_MPB_i>-1){
      c->cv[ws->KI_MPB_i] = 0.0;
      c->cv[ws->yCfac_MPB_i] = 0.0;
    }
    //if (ws->Kd_490_i > -1){ // i.e. outputting
    //  y[ws->Kd_490_i] = 0.0;
    //}
    if (ws->PAR_i > -1){ // i.e. outputting
      y[ws->PAR_i] = 0.0;
    }
    if (ws->PAR_z_i > -1){ // i.e. outputting
      y[ws->PAR_z_i] = 0.0;
    }
    return;
  }

  double *aA_s_MPB;   
  aA_s_MPB = d_alloc_1d(ws->num_waves);
  
  if (ws->KI_MPB_i>-1){
    
    c->cv[ws->KI_MPB_i] = 0.0;
    c->cv[ws->yCfac_MPB_i] = 0.0;
    
    /* now do it for microphytobenthos */
    
    double rad = ws->rad_MPB;
    double vol = 4.0 / 3.0 * M_PI * rad * rad * rad;
    double m   = ws->m_MPB;    
    double sumlight,yCfac;
    
    double MPB_N = y[ws->MPB_N_i];
    double MPB_Chl = y[ws->MPB_Chl_i];
    
    double cellnum =  MPB_N / (m * red_A_N * 1000.0 * MW_Nitr); /* cell m-3 */
    double cellChl = MPB_Chl / (vol * cellnum); /* cellular Chl conc */
    
    // intermediate var
    
    double *absorbance;
        
    absorbance = d_alloc_1d(ws->num_waves);
    
    /* if it is dark, assume equal weighting for all wavelengths. */
    
    sumlight = 0.0;
    yCfac = 0.0;
    
    for (w=0; w<ws->num_waves; w++){
      absorbance[w] = ws->yC_MPB[w] * cellChl;
      yCfac = yCfac + diffaa(absorbance[w],rad)*max(lighttop_s[w]*ws->wave[w],0.0001);
      sumlight = sumlight + max(lighttop_s[w]*ws->wave[w],0.0001);
    }
    c->cv[ws->yCfac_MPB_i] = yCfac/sumlight;
    
    aawave(rad, absorbance, aA_s_MPB, ws->num_waves,ws->wave);
    d_free_1d(absorbance);
   
  }else{
    for (w=0; w<ws->num_waves; w++){
      aA_s_MPB[w] = 0.0;
    }
  }
  double dz = c->dz_sed;
  energy = 0.0;
  if (ws->PAR_z_i > -1){ /* i.e. outputting */
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      energy += lighttop_s[w] * 8.359335857479461e-09 * ws->wave[w];
    }	
    y[ws->PAR_z_i] = energy;
  }

  for (w=0; w<ws->num_waves; w++){
    lighttop = lighttop_s[w]; 
    Kd_s[w] = (bio->kw_s[w] + y[ws->MPB_N_i] / (ws->m_MPB*red_A_N * 1000.0 * MW_Nitr) * aA_s_MPB[w]) * dz;   /* vertical attenuation */
    lighttop_s[w] = lighttop * exp(-Kd_s[w]);
    light_s[w] = (lighttop - lighttop_s[w] ) / Kd_s[w];
    c->cv[ws->KI_MPB_i] += light_s[w] * 8.359335857479461e-09 * ws->wave[w] * aA_s_MPB[w];
  }

  energy = 0.0;
  
  if (ws->PAR_i > -1){ /* i.e. outputting */
    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      energy += light_s[w] * 8.359335857479461e-09 * ws->wave[w];
    }	
    y[ws->PAR_i] = energy;
  }

  for (w=0; w<ws->num_waves; w++){
    lighttop_s[w] = 0.0;  // only top layer has light.
  }

  // free intermediate var
    if (aA_s_MPB != NULL)
      d_free_1d(aA_s_MPB);

}
