/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_uq_epi.c
 *  
 *  Description:
 *
 *  Epi-benthic spectrally-resolved optical model. This is the standard routine for B3p0.
 *  
 *  Model calculates:
 *  
 *  - absorption by each model autotroph in precalc (output time minus ECOLOGY_DT).
 *  - the component of light returned to the surface from the bottom (at the output time).
 *  - since light_spectral_sed does not return light to the surface, all of the simulated 
 *    satellite products are caluclated here (output time).
 *  
 *  WARNING: variables last calculated in precalc (PAR etc.) are output at the output time but are the value from 
 *           output time - ECOLOGY_DT. 
 *           Variables calculated in postcalc are at the correct time.
 *  
 *  To output remote-sensing relectance at wavelenth XXX requires a 2D tracer named: R_XXX 
 *
 *  Optical model described in:
 *  
 *  Baird, M. E., N. Cherukuru, E. Jones, N. Margvelashvili, M. Mongin, K. Oubelkheir, P. J. Ralph, F. Rizwi, 
 *  B. J. Robson, T. Schroeder, J. Skerratt, A. D. L. Steven and K. A. Wild-Allen (2016) Remote-sensing 
 *  reflectance and true colour produced by a coupled hydrodynamic, optical, sediment, biogeochemical 
 *  model of the Great Barrier Reef, Australia: comparison with satellite data. Env. Model. Software 78: 79-96.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: light_spectral_uq_epi.c 6375 2019-11-06 02:15:11Z bai155 $
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
#include "light_spectral_uq_epi.h"
#include "constants.h"

typedef struct {
  /*
   * parameters
   */
  double CSm;
  double CSrad;
  double CSvol;

  double MPBm;
  double MPBrad;
  double MPBvol;

  /* epibenthic variables */

  int Epilight_s_i;
  int Epilight_i;

  int num_waves;
  double *wave;

  /* tracers */

  int MA_N_i;
  int SG_N_i;
  int SGH_N_i;
  int SGD_N_i;
  int CH_N_i;
  int CS_N_i;
  int CS_Chl_i;

  int CS_NR_i;
  int CS_Xp_i;
  int CS_Xh_i;

  int Mud_sed_i;  
  int Finesed_sed_i;
  int Sand_sed_i;
  int Dust_sed_i;
  int CarbSand_sed_i;
  
  int Mud_mineral_sed_i;
  int Mud_carbonate_sed_i;
  int Sand_carbonate_sed_i;
  int Sand_mineral_sed_i;

  int MPB_N_sed_i;
  int MPB_Chl_sed_i;

  /* common variables */

  int cv_lighttop_s_i;
  int lighttop_i;
  int KI_SG_i;
  int KI_SGH_i;
  int KI_SGD_i;
  int KI_MA_i;
  int KI_CS_i;
  int yCfac_CS_i;

  int wPAR_bot;
  int wPAR_top;

  int w443;  
  int w488;
  int w547;
  int w510;
  int w560;
  int w645;
  int w486;
  int w551;
  int w678;
  int w667;
  int w748;
  int w709;
  int w665;
  int w490;

  int OC3M_i;
  int TSSM_i;
  int KD490M_i;
  int OC3V_i;
  int OC4Me_i;
  int nFLH_i;
  int Hue_i;

  int *RRS_i;

  int w_bot_i;
  int u_surf_i;
  int SWR_bot_abs_i;
  int num_rrs_waves;
  double *rrs_wave;

  int z_secchi_i;
  
  double bbcell1;
  double bbcell2;
  double xan2chl_MPB;

  double g0;
  double g1;

  /* Local variables */

  double MAleafden;
  double SGleafden;
  double SGHleafden;
  double SGDleafden;
  double SGorient;
  double SGHorient;
  double SGDorient;
  double CHpolypden;
  double CHarea;

  /*
   * output variables
   */

  int EpiPAR_i;
  int EpiPAR_sg_i;
  int Secchi_i;
  int Zenith_i;

  int cv_moonlight_i;
  int cv_lunar_zenith_i;
  int cv_lunar_phase_i;
  int cv_moon_fulldisk_i;
  int Moonlight_i;
  int Lunar_zenith_i;
  int Lunar_phase_i;
  int Moon_fulldisk_i;
  

} workspace;

static void finalise_reflectances(workspace *ws, cell *c);

void light_spectral_uq_epi_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;
  workspace* ws = malloc(sizeof(workspace));
  int w, w2, count_Li;
  
  int OFFSET_SED = tracers->n;
  int OFFSET_EPI = tracers->n * 2;

  bio_opt_prop *bio = e->bio_opt;
  
  p->workspace = ws;

  /* try index in init for the benefit of RECOM */

  ws->SG_N_i = e->try_index(epis, "SG_N", e);
  if (ws->SG_N_i >=0)
    ws->SG_N_i += OFFSET_EPI;
  
  ws->MA_N_i = e->try_index(epis, "MA_N", e); 
  if (ws->MA_N_i >=0)
    ws->MA_N_i += OFFSET_EPI;
  
  ws->SGH_N_i = e->try_index(epis, "SGH_N", e);
  if (ws->SGH_N_i >=0)
    ws->SGH_N_i += OFFSET_EPI;

  ws->SGD_N_i = e->try_index(epis, "SGD_N", e);
  if (ws->SGD_N_i >=0)
    ws->SGD_N_i += OFFSET_EPI;

  ws->CH_N_i = e->try_index(epis, "CH_N", e);
  if (ws->CH_N_i >=0)
    ws->CH_N_i += OFFSET_EPI;

  ws->CS_N_i = e->try_index(epis, "CS_N", e);
  if (ws->CS_N_i >=0)
    ws->CS_N_i += OFFSET_EPI;

  ws->CS_Chl_i = e->try_index(epis, "CS_Chl", e);
  if (ws->CS_Chl_i >=0)
    ws->CS_Chl_i += OFFSET_EPI;

  ws->EpiPAR_i = e->try_index(epis, "EpiPAR", e);
  if (ws->EpiPAR_i >=0)
    ws->EpiPAR_i += OFFSET_EPI;
  
  ws->EpiPAR_sg_i = e->try_index(epis, "EpiPAR_sg", e);
  if (ws->EpiPAR_sg_i >=0)
    ws->EpiPAR_sg_i += OFFSET_EPI;

  ws->Secchi_i = e->try_index(epis, "Secchi", e);
  if (ws->Secchi_i >=0){
    ws->Secchi_i += OFFSET_EPI;
    ws->z_secchi_i = find_index_or_add(e->cv_column, "z_secchi", e);
  }
  
  ws->Zenith_i = e->try_index(epis, "Zenith2D", e);
  if (ws->Zenith_i >=0)
    ws->Zenith_i += OFFSET_EPI;

  ws->Moonlight_i = e->try_index(epis, "Moonlight", e);
  if (ws->Moonlight_i >=0){
    ws->Moonlight_i += OFFSET_EPI;
    ws->cv_moonlight_i = find_index_or_add(e->cv_column, "moonlight", e);
  }

  if (ws->Moonlight_i >=0){

    ws->Lunar_zenith_i = e->try_index(epis, "Lunar_zenith", e);
    if (ws->Lunar_zenith_i >=0){
      ws->Lunar_zenith_i += OFFSET_EPI;
      ws->cv_lunar_zenith_i = find_index_or_add(e->cv_column, "lunar_zenith", e);
    }
    
    ws->Lunar_phase_i = e->try_index(epis, "Lunar_phase", e);
    if (ws->Lunar_phase_i >=0){
      ws->Lunar_phase_i += OFFSET_EPI;
      ws->cv_lunar_phase_i = find_index_or_add(e->cv_column, "lunar_phase", e);
    }
    
    ws->Moon_fulldisk_i = e->try_index(epis, "Moon_fulldisk", e);
    if (ws->Moon_fulldisk_i >=0){
      ws->Moon_fulldisk_i += OFFSET_EPI;
      ws->cv_moon_fulldisk_i = find_index_or_add(e->cv_column, "moon_fulldisk", e);
    }
  }
  
  ws->SWR_bot_abs_i = e->try_index(epis, "SWR_bot_abs", e);
  if (ws->SWR_bot_abs_i >= 0)
    ws->SWR_bot_abs_i += OFFSET_EPI;
  
  ws->Epilight_i = e->try_index(epis, "Epilight", e);
  if (ws->Epilight_i >= 0)
    ws->Epilight_i += OFFSET_EPI;
  
  ws->OC3M_i = e->try_index(epis, "OC3M", e);
  if (ws->OC3M_i >= 0)
    ws->OC3M_i += OFFSET_EPI;

  ws->OC3V_i = e->try_index(epis, "OC3V", e);
  if (ws->OC3V_i >= 0)
    ws->OC3V_i += OFFSET_EPI;

  ws->OC4Me_i = e->try_index(epis, "OC4Me", e);
  if (ws->OC4Me_i >= 0)
    ws->OC4Me_i += OFFSET_EPI;

  ws->nFLH_i = e->try_index(epis, "nFLH", e);
  if (ws->nFLH_i >= 0)
    ws->nFLH_i += OFFSET_EPI;

  ws->Hue_i = e->try_index(epis, "Hue", e);
  if (ws->Hue_i >= 0)
    ws->Hue_i += OFFSET_EPI;

  ws->TSSM_i = e->try_index(epis, "TSSM", e);
  if (ws->TSSM_i >= 0)
    ws->TSSM_i += OFFSET_EPI;

  ws->KD490M_i = e->try_index(epis, "KD490M", e);
  if (ws->KD490M_i >= 0)
    ws->KD490M_i += OFFSET_EPI;
  
  ws->Mud_sed_i = e->try_index(tracers, "Mud", e);
  if (ws->Mud_sed_i >= 0)
    ws->Mud_sed_i += OFFSET_SED;
  
  ws->Sand_sed_i = e->try_index(tracers, "Sand", e);
  if (ws->Sand_sed_i >= 0)
    ws->Sand_sed_i += OFFSET_SED;
  
  ws->Finesed_sed_i = e->try_index(tracers, "Finesed", e);
  if (ws->Finesed_sed_i >= 0)
    ws->Finesed_sed_i += OFFSET_SED;
  
  ws->Dust_sed_i = e->try_index(tracers, "Dust", e);
  if (ws->Dust_sed_i >= 0)
    ws->Dust_sed_i += OFFSET_SED;

  ws->CarbSand_sed_i = e->try_index(tracers, "CarbSand", e);
  if (ws->CarbSand_sed_i >= 0)
    ws->CarbSand_sed_i += OFFSET_SED;

  ws->Mud_carbonate_sed_i = e->try_index(tracers, "Mud-carbonate", e);
  if (ws->Mud_carbonate_sed_i >= 0)
    ws->Mud_carbonate_sed_i += OFFSET_SED;
  
  ws->Sand_carbonate_sed_i = e->try_index(tracers, "Sand-carbonate", e);
  if (ws->Sand_carbonate_sed_i >= 0)
    ws->Sand_carbonate_sed_i += OFFSET_SED;

  ws->Mud_mineral_sed_i = e->try_index(tracers, "Mud-mineral", e);
  if (ws->Mud_mineral_sed_i >= 0)
    ws->Mud_mineral_sed_i += OFFSET_SED;
  
  ws->Sand_mineral_sed_i = e->try_index(tracers, "Sand-mineral", e);
  if (ws->Sand_mineral_sed_i >= 0)
    ws->Sand_mineral_sed_i += OFFSET_SED;

  ws->MPB_N_sed_i = e->try_index(tracers, "MPB_N", e);
  if (ws->MPB_N_sed_i >= 0)
    ws->MPB_N_sed_i += OFFSET_SED;
  
  ws->MPB_Chl_sed_i = e->try_index(tracers, "MPB_Chl", e);
  if (ws->MPB_Chl_sed_i >= 0)
    ws->MPB_Chl_sed_i += OFFSET_SED;
  

/* Get wavelengths for light */

  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");
  
  /*  common column variables */

  ws->cv_lighttop_s_i = find_index_or_add_vector(e->cv_column, "lighttop_s", e, ws->num_waves);
  ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);
  
}

void light_spectral_uq_epi_postinit(eprocess* p)
{
  double *absorbance;
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;
  workspace* ws = (workspace *)p->workspace;
  double NtoCHL = get_parameter_value(e, "NtoCHL");
  char buf[MAXSTRLEN];
  int OFFSET_EPI = tracers->n * 2;

  int w,w2;

/* for surface reflectance weighting */

  ws->num_rrs_waves = e->num_rsr_waves;
  ws->rrs_wave      = e->rsr_waves;

  ws->w_bot_i  = find_index_or_add_vector(e->cv_column, "w_bot", e, ws->num_rrs_waves);
  ws->u_surf_i = find_index_or_add_vector(e->cv_column, "u_surf", e, ws->num_rrs_waves);

  /*
   * Cache indices and wavelengths of the variables to output
   */
  ws->RRS_i = NULL;
  if (ws->num_rrs_waves) 
    ws->RRS_i = i_alloc_1d(ws->num_rrs_waves);
  for (w2=0; w2<ws->num_rrs_waves; w2++) {
    int i_tmp;
    sprintf(buf, "R_%d", (int)ws->rrs_wave[w2]);
    i_tmp = e->try_index(epis, buf, e);
    if (i_tmp > -1)
      ws->RRS_i[w2] = i_tmp + OFFSET_EPI;
    else
      e->quitfn("light_spectral_uq_epi: Could not find %s\n", buf);
  }

/* set light indexes for PAR and OC3M so that code isn't looking for them each time step. */

  for (w=0; w<ws->num_waves; w++){
    if (ws->wave[w] == 410.0)
      ws->wPAR_bot = w;
    if (ws->wave[w] == 690.0)
      ws->wPAR_top = w;
  }
  for (w2=0; w2<ws->num_rrs_waves; w2++){
    if (ws->rrs_wave[w2] == 443.0)
      ws->w443 = w2;
    if (ws->rrs_wave[w2] == 488.0)
      ws->w488 = w2;
    if (ws->rrs_wave[w2] == 510.0)
      ws->w510 = w2;
    if (ws->rrs_wave[w2] == 547.0)
      ws->w547 = w2;
    if (ws->rrs_wave[w2] == 560.0)
      ws->w560 = w2;
    if (ws->rrs_wave[w2] == 645.0)
      ws->w645 = w2;
    if (ws->rrs_wave[w2] == 486.0)
      ws->w486 = w2;
    if (ws->rrs_wave[w2] == 551.0)
      ws->w551 = w2;
    if (ws->rrs_wave[w2] == 667.0)
      ws->w667 = w2;
    if (ws->rrs_wave[w2] == 678.0)
      ws->w678 = w2;
    if (ws->rrs_wave[w2] == 748.0)
      ws->w748 = w2;
    if (ws->rrs_wave[w2] == 709.0)
      ws->w709 = w2;
    if (ws->rrs_wave[w2] == 665.0)
      ws->w665 = w2;
    if (ws->rrs_wave[w2] == 490.0)
      ws->w490 = w2;
  }

  ws->CS_Xp_i = e->try_index(epis, "CS_Xp", e);
  if (ws->CS_Xp_i >=0)
    ws->CS_Xp_i += OFFSET_EPI;
  
  ws->CS_NR_i = e->try_index(epis, "CS_NR", e);
  if (ws->CS_NR_i >=0)
    ws->CS_NR_i += OFFSET_EPI;
  
  ws->CS_Xh_i = e->try_index(epis, "CS_Xh", e);
  if (ws->CS_Xh_i >=0)
    ws->CS_Xh_i += OFFSET_EPI;
  
  ws->MAleafden = 0.0;
  if (ws->MA_N_i > -1){
    
    /* Put in macroalgae initialisations if there is a macroalgae tracer */

    ws->MAleafden = get_parameter_value(e, "MAleafden");
    ws->KI_MA_i = find_index_or_add(e->cv_cell, "KI_MA", e);

  }

  ws->SGleafden = 0.0;
  ws->SGorient  = 0.0;
  if (ws->SG_N_i > -1){

    /* Put in Zostera initialisations if there is a Zostera tracer */

    ws->SGleafden = get_parameter_value(e, "SGleafden");
    ws->SGorient = try_parameter_value(e, "SGorient");
    if (isnan(ws->SGorient))
      ws->SGorient = 1.0;
    ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SG", e);

  }

  ws->SGHleafden = 0.0;
  ws->SGHorient  = 0.0;
  if (ws->SGH_N_i > -1){

    /* Put in Halophila initialisations if there is a Halophila tracer */

    ws->SGHleafden = get_parameter_value(e, "SGHleafden");
    ws->SGHorient = try_parameter_value(e, "SGHorient");
    if (isnan(ws->SGHorient))
      ws->SGHorient = 1.0;
    ws->KI_SGH_i = find_index_or_add(e->cv_cell, "KI_SGH", e);
  }

  ws->SGDleafden = 0.0;
  ws->SGDorient  = 0.0;
  if (ws->SGD_N_i > -1){

    /* Put in Deep seagrass initialisations if there is a tracer */

    ws->SGDleafden = get_parameter_value(e, "SGDleafden");
    ws->SGDorient = try_parameter_value(e, "SGDorient");
    if (isnan(ws->SGDorient))
      ws->SGDorient = 1.0;
    ws->KI_SGD_i = find_index_or_add(e->cv_cell, "KI_SGD", e);
  }

  if (ws->MPB_N_sed_i > -1){ 
    
    ws->MPBrad = get_parameter_value(e, "MBrad");
    ws->MPBm = PhyCellMass(ws->MPBrad);
    ws->MPBvol = 4.0 / 3.0 * M_PI * ws->MPBrad * ws->MPBrad * ws->MPBrad;
    
  }

  if (ws->CS_N_i > -1){

    /* Only in coral initialisations if there is a coral tracer */

    ws->CHarea = try_parameter_value(e, "CHarea");
    if (isnan(ws->CHarea))
      ws->CHarea = 1.0;
    
    ws->CHpolypden = get_parameter_value(e, "CHpolypden");
    ws->CSrad = get_parameter_value(e, "CSrad");
    ws->CSm = PhyCellMass(ws->CSrad);
    ws->CSvol = 4.0 / 3.0 * M_PI * ws->CSrad * ws->CSrad * ws->CSrad;
    
    ws->KI_CS_i = find_index_or_add(e->cv_cell, "KI_CS", e);
    ws->yCfac_CS_i = find_index_or_add(e->cv_cell, "yCfac_CS", e);
  }

    /* ********************************************************************************/
    /* Initialise spectrally-resolved surface reflectance values not already obtained */

    /* air-water and solid angle reflectance parameters */

    ws->g0 = try_parameter_value(e, "g0");
    if (isnan(ws->g0)){
      ws->g0 = 0.0895;   // Brando 2012
      eco_write_setup(e,"Code default of g0 = %e \n",ws->g0);
    }
    ws->g1 = try_parameter_value(e, "g1");
    if (isnan(ws->g1)){
      ws->g1 = 0.1247;   // Brando 2012
      eco_write_setup(e,"Code default of g1 = %e \n",ws->g1);
    }
    ws->bbcell1 = try_parameter_value(e, "bbcell1");
    if (isnan(ws->bbcell1)){
      ws->bbcell1 = 4.26e-14;   // Whitmore 2012
      eco_write_setup(e,"Code default of bbcell1  = %e \n",ws->bbcell1);
    }
    ws->bbcell2 = try_parameter_value(e, "bbcell2");
    if (isnan(ws->bbcell2)){
      ws->bbcell2 = 2.028;   // Whitmore 2012
      eco_write_setup(e,"Code default of bbcell2  = %e \n",ws->bbcell2);
    }

    ws->xan2chl_MPB = get_parameter_value(e, "MBxan2chl");

    /**********************************************************************************/ 
}

void light_spectral_uq_epi_destroy(eprocess* p)
{
  workspace* ws = p->workspace;
  if (ws->RRS_i)
    i_free_1d(ws->RRS_i);
  free(ws);
}

void light_spectral_uq_epi_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  double* lighttop_s = col->cv[ws->cv_lighttop_s_i];
  
  int w,i;

  int num_waves     = ws->num_waves;
  int num_rrs_waves = ws->num_rrs_waves;

  double *wave  = ws->wave;
  double *rrs_wave  = ws->rrs_wave;

  double *aA_s_CS = NULL;
  double *aA_s_CS_cph = NULL;

  bio_opt_prop *bio = e->bio_opt;

  double energy = 0.0;

  /* Need this here for the refletance bit */

  double SA = 0.0;

  for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
    energy += lighttop_s[w];
  }

  if (energy < 0.01){
    if (ws->MA_N_i > -1){
      c->cv[ws->KI_MA_i] = 0.0;
    }
    if (ws->EpiPAR_sg_i > -1){
      y[ws->EpiPAR_sg_i] = 0.0;
    }
    if (ws->SG_N_i > -1){
      c->cv[ws->KI_SG_i] = 0.0;
    }
    if (ws->SGH_N_i > -1){
      c->cv[ws->KI_SGH_i] = 0.0;
    }
    if (ws->SGD_N_i > -1){
      c->cv[ws->KI_SGD_i] = 0.0;
    }
    if (ws->CS_N_i > -1){
      c->cv[ws->KI_CS_i] = 0.0;
      c->cv[ws->yCfac_CS_i] = 0.0;
    }

    double *cv_lighttop = col->cv[ws->lighttop_i];
    cv_lighttop[0] = 0.0;
   
    return;
  }

    /* Model state */

  double photons = 0.0;

  /* Sum up energy, convert to photons, and store flux for use in growth routines 
     
     Absorption by macroalgae before seagrass */

      /* photons = leaf absorbance * leaf area fraction

         % N dry weight = 1.92 +/- 0.05, Duarte (1990), MEPS 67: 201-207.

         leaf dry  weight ~ 20 mg Hasnsen (2000), MEPS 199, 83-96.

         leaf dimension = 40 cm x 4 mm Kemp (1987), MEPS 41, 79-86

         no.of leaves n = SG_N / (20 * 0.0192) 

         leaf area (m2) = n * 0.4 * 0.004
  
                        =  SG_N * 0.0042                          

         where 0.0042 has units m2 (mg N m-2)-1, and is in parameter file, 
	 
	 now changed to g N m-2 in parameter file */

    /* 8.359335857479461e-9 = (1e9 h c)-1 / Av 
       where h = Plank const, c = speed of light, Av = Avogadro const, 1e9 nm m-1  */

  if (ws->MA_N_i > -1){
    if (ws->KI_MA_i > -1)
      c->cv[ws->KI_MA_i] = 0.0;
    
    double MA_N = y[ws->MA_N_i];
    
    if (MA_N > 0.0){  /* macroalgae */
      
      photons = 0.0;
      
      for (w=0; w<num_waves; w++){
	photons += lighttop_s[w] * (1.0 - exp(-MA_N * ws->MAleafden * bio->MA_aAwave[w])) 
	  * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-MA_N * ws->MAleafden * bio->MA_aAwave[w]);
      }
      c->cv[ws->KI_MA_i] = photons;

    }
  }
  
  if (ws->EpiPAR_sg_i > -1) {

  double light_above_sg = 0.0;  // quantum-weighted and in photon m-2 d-1

    for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
      light_above_sg = light_above_sg + lighttop_s[w] * wave[w];
    }
    y[ws->EpiPAR_sg_i] = light_above_sg * 86400.0 * 8.359335857479461e-09;
  }
  
  /* Zostera first */

  if (ws->SG_N_i > -1){

    double SG_N = y[ws->SG_N_i];

    if (ws->KI_SG_i > -1)
      c->cv[ws->KI_SG_i] = 0.0;

    if (SG_N > 0.0){  /* seagrass */
      
      photons = 0.0;
      
      for (w=0; w<num_waves; w++){
	photons += lighttop_s[w] * (1.0 - exp(-SG_N * ws->SGorient * ws->SGleafden * bio->SG_aAwave[w])) 
	  * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SG_N * ws->SGorient * ws->SGleafden * bio->SG_aAwave[w]);  
      }
      c->cv[ws->KI_SG_i] = photons;
    }
  }
  
  /* Then Halophila if present */

  if (ws->SGH_N_i > -1){

    double SGH_N = y[ws->SGH_N_i];

    if (ws->KI_SGH_i > -1)
      c->cv[ws->KI_SGH_i] = 0.0;

    if (SGH_N > 0.0){  /* seagrass */

      photons = 0.0;

      for (w=0; w<num_waves; w++){

	photons += lighttop_s[w] * (1.0 - exp(-SGH_N * ws->SGHorient * ws->SGHleafden * bio->SGH_aAwave[w])) 
                                            * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SGH_N * ws->SGHorient * ws->SGHleafden * bio->SGH_aAwave[w]);
      }
      c->cv[ws->KI_SGH_i] = photons;

    }
  }

  /* Then Deep seagrass if present */

  if (ws->SGD_N_i > -1){

    double SGD_N = y[ws->SGD_N_i];

    if (ws->KI_SGD_i > -1)
      c->cv[ws->KI_SGD_i] = 0.0;

    if (SGD_N > 0.0){  /* seagrass */

      photons = 0.0;

      for (w=0; w<num_waves; w++){

	photons += lighttop_s[w] * (1.0 - exp(-SGD_N * ws->SGDorient * ws->SGDleafden * bio->SGH_aAwave[w])) 
                                            * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SGD_N * ws->SGDorient * ws->SGDleafden * bio->SGH_aAwave[w]);
      }
      c->cv[ws->KI_SGD_i] = photons;

    }
  }

  
  if (ws->CS_N_i > -1){

    double CH_N = y[ws->CH_N_i];
    double CS_N = y[ws->CS_N_i];
    double CS_Chl = y[ws->CS_Chl_i];
    if (ws->KI_CS_i > -1){
      c->cv[ws->KI_CS_i] = 0.0;
      c->cv[ws->yCfac_CS_i] = 0.0;
    }
    if (CS_N > 0.0){  /* Zoothanthellae in corals - represented like cells not leaves */
                      /* so need to find an average light in the layer like in wc */
      
      /* Need intracellular chlorophyll concentration - since cells aren't advected, 
	 I could have intra. chl conc as the state variable. But leave as is, incase we 
	 want to have the ability of zoothanthellae being expelled into water column */
      
      double cellnum = CS_N / (ws->CSm * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
      double cellChl = CS_Chl / (ws->CSvol * cellnum); /* cellular Chl conc */
      
      // intermediate var
  
      double *absorbance;
      
      absorbance = d_alloc_1d(num_waves);     
      aA_s_CS = d_alloc_1d(num_waves);

      if (ws->CS_Xh_i > -1){

	double cellXp = y[ws->CS_Xp_i] / (ws->CSvol * cellnum); /* cellular Xp conc */

	/* Only include photosynthetic carts in KI calulations */
	
	for (w=0; w<num_waves; w++){
	  absorbance[w] = cellXp * bio->yC_diadinoxanthin[w] + cellChl * bio->yC_symbiodinium[w];
	}
        /* absorption cross section for all pigments */ 

	aA_s_CS_cph = d_alloc_1d(num_waves);
      }else{
	absorbwave(absorbance,cellChl, num_waves, wave);  // absorbance just with Chl 
      }
      
      /* weight self-shading coefficient based on photons incident on the top of 
	 the layer - multiply by wave[w] to convert energy to photons, but ignore constant 
	 (1e9 h c)-1 / Av is constant across all wavelengths */
      
      double sumlight = 0.0;
      double yCfac = 0.0;
      double kd;
      double Itop;
      double meanlight;

      for (w=0; w<num_waves; w++){
	yCfac = yCfac + diffaa(absorbance[w],ws->CSrad)*max(lighttop_s[w]*wave[w],0.00001);
	sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
      }
      c->cv[ws->yCfac_CS_i] = yCfac/sumlight;
      aawave(ws->CSrad, absorbance, aA_s_CS, num_waves,wave);

      /* now do again including heat dissappating pigments for light attenuation  
	 - write over version of yCfac only if there are heat disspiating pigments */

      if (ws->CS_Xh_i > -1){
	sumlight = 0.0;
	yCfac = 0.0;
	double cellXh = y[ws->CS_Xh_i] / (ws->CSvol * cellnum); /* cellular Xh conc */
	for (w=0; w<num_waves; w++){
	  absorbance[w] += cellXh * bio->yC_diatoxanthin[w];
	  yCfac = yCfac + diffaa(absorbance[w],ws->CSrad)*max(lighttop_s[w]*wave[w],0.00001);
	  sumlight = sumlight + max(lighttop_s[w]*wave[w],0.00001);
	}
	c->cv[ws->yCfac_CS_i] = yCfac/sumlight;
	aawave(ws->CSrad, absorbance, aA_s_CS_cph, num_waves,wave);
      }

      // free intermediate var
      d_free_1d(absorbance);

      /* Integrate photons absorbed across all wavelengths */

      photons = 0.0;

      /* establish surface area of coral polyps - only attenuate within this area */ 

      SA = ws->CHarea * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)); /* dimensionless */

      for (w=0; w<num_waves; w++){

	/* 1e-12 - just very low light to avoid dividing by zero. */

	if (ws->CS_Xh_i > -1){
	  kd = max(aA_s_CS_cph[w] * cellnum / SA,1e-12);
	}else{
	  kd = max(aA_s_CS[w] * cellnum / SA,1e-12);
	}

        Itop = lighttop_s[w];

        /* use CHarea to establish what misses coral communities - also above needed to divide cellnum to account for squeezing. */ 

	lighttop_s[w] = (1.0 - SA) * lighttop_s[w] + SA * lighttop_s[w] * exp(-kd);
	
	// meanlight = (Itop - lighttop_s[w]) / kd;

        meanlight = Itop * (1.0 - exp(-kd)) / kd;

	photons += meanlight * 8.359335857479461e-09 * wave[w] * aA_s_CS[w];
      }    
      c->cv[ws->KI_CS_i] = photons;
    }
  }
  
  if (aA_s_CS != NULL){
    d_free_1d(aA_s_CS);
    if (ws->CS_Xp_i > -1)
      d_free_1d(aA_s_CS_cph);
  }
  
  
  /* sum up light in the PAR range (400-700 nm), assuming weighting by photons 
     (not energy) - remove when sediment is spectrally-resolved.  */
  
  double lighttop = 0.0;
  double wavesum = 0.0;
  
  double *cv_lighttop = col->cv[ws->lighttop_i];

  for (w=ws->wPAR_bot; w <= ws->wPAR_top; w++){
    lighttop = lighttop + lighttop_s[w]*wave[w];
    wavesum = wavesum + wave[w];
  }

  // mol photon m-2 d-1 

  if (ws->EpiPAR_i > -1)
    y[ws->EpiPAR_i] = lighttop * 86400.0 * 8.359335857479461e-09;
  
  lighttop = lighttop/wavesum; 
  cv_lighttop[0] = lighttop;
}

void light_spectral_uq_epi_postcalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  
  int w,i;

  int num_waves     = ws->num_waves;
  int num_rrs_waves = ws->num_rrs_waves;

  double *wave  = ws->wave;
  double *rrs_wave  = ws->rrs_wave;

  double *aA_s_CS = NULL;

  bio_opt_prop *bio = e->bio_opt;

  /* Moonlight */
  
  if (ws->Moonlight_i > -1){
    
    double* cv_moonlight = col->cv[ws->cv_moonlight_i];
    y[ws->Moonlight_i] = cv_moonlight[0];
   
    if (ws->Lunar_zenith_i > -1){
      double* cv_lunar_zenith = col->cv[ws->cv_lunar_zenith_i];
      y[ws->Lunar_zenith_i] = cv_lunar_zenith[0];
    }

    if (ws->Lunar_phase_i > -1){
      double* cv_lunar_phase = col->cv[ws->cv_lunar_phase_i];
      y[ws->Lunar_phase_i] = cv_lunar_phase[0];
    }

    if (ws->Moon_fulldisk_i > -1){
      double* cv_moon_fulldisk = col->cv[ws->cv_moon_fulldisk_i];
      y[ws->Moon_fulldisk_i] = cv_moon_fulldisk[0];
    }
  }

  double zenith = einterface_calc_zenith(e->model, e->t + e->dt, c->b);
  
  if (ws->Zenith_i > -1)
    y[ws->Zenith_i] = zenith;

  zenith = (fabs(zenith)<1.0e-15)? 1.0e-15:zenith;
  zenith = ((zenith-M_PI/2.0)-fabs(zenith-M_PI/2.0))/2.0+M_PI/2.0;

  if (e->num_rsr_waves) {
    double* u_surf = col->cv[ws->u_surf_i];
    int w2;
    
    if (zenith > (M_PI/2.0 - 0.01)){
      for (w2=0; w2<num_rrs_waves; w2++){
	y[ws->RRS_i[w2]] = 0.0;
      }
      if (ws->OC3M_i > -1)
	y[ws->OC3M_i]= 0.0;
      if (ws->OC4Me_i > -1)
	y[ws->OC4Me_i]= 0.0;
      if (ws->nFLH_i > -1)
	y[ws->nFLH_i]= 0.0;
      if (ws->Hue_i > -1)
	y[ws->Hue_i]= 0.0;
      if (ws->OC3V_i > -1)
	y[ws->OC3V_i]= 0.0;
      if (ws->TSSM_i > -1)
	y[ws->TSSM_i]= 0.0;
      if (ws->KD490M_i > -1)
	y[ws->KD490M_i]= 0.0;
      if (ws->Secchi_i > -1)
	y[ws->Secchi_i] =  0.0;
      return;
    }
    
    /*
     * Section to calculate reflectances
     */
    
    double* w_bot  = col->cv[ws->w_bot_i];
    double* u_bot = d_alloc_1d(num_rrs_waves);
    
    for (w2=0; w2<num_rrs_waves; w2++){
      u_bot[w2] = 0.0;
    }    
    double all_sed = 0.0;
    
    double *rhoskel    = bio->coral_skeleton_srs;
    double *rhosand    = bio->sand_srs;
    double *rhomud     = bio->mud_srs;
    double *rhofinesed = bio->finesed_srs;

    double *rhoseagrass = bio->seagrass_srs;
    double *rhomacroalgae = bio->macroalgae_srs;

    double Mud, Sand, Finesed, Dust, CarbSand;

    double f_MA = 0.0;
    double f_SG = 0.0;
    double f_SGH = 0.0;
    double f_SGD = 0.0;
    double f_polyp = 0.0;
    double f_mpb = 0.0;
    double f_zoo = 0.0;
    double f_skel = 0.0;
    double frac = 0.0;

    /*
     * Here we progress from the tallest species down to the
     * sediments, calculating the fraction of light that is still
     * penetrating down after each reflecting thing
     *
     * Note: variable dependent things like ws->MAleafden should've
     *       been set to zero above in the case where there is no MA
     *       and so on...
     * frac - what is available at that depth.
     * u_bot - addition to the bottom.
     */

    /* 
     * Macroalgae 
     */
    if (ws->MA_N_i > -1) {
      double MA_N = y[ws->MA_N_i];
      
      /* f_MA = (1 - exp(omega_MA * MA)) */
      f_MA = (1.0 - exp(-MA_N * ws->MAleafden));
      for (w2 = 0; w2 < ws->num_rrs_waves; w2++)
	u_bot[w2] += f_MA * rhomacroalgae[w2];
    }

    /*
     * Seagrass 
     */
    if (ws->SG_N_i > -1) {
      double SG_N = y[ws->SG_N_i];
      
      /* f_SG = (1 - f_MA) (1 - exp(omega_SG * SG)) */
      f_SG = (1.0 - f_MA) * (1.0 - exp(-SG_N * ws->SGorient * ws->SGleafden));
      for (w2 = 0; w2 < ws->num_rrs_waves; w2++)
	u_bot[w2] += f_SG * rhoseagrass[w2];
    }
    
    /*
     * Halophila
     */
    if (ws->SGH_N_i > -1) {
      double SGH_N = y[ws->SGH_N_i];
    
      /* f_SGH = (1 - f_MA - f_SG) (1 - exp(omega_SGH * SGH)) */
      f_SGH = (1.0 - f_MA - f_SG) * (1.0 - exp(-SGH_N * ws->SGHorient * ws->SGHleafden));
      for (w2 = 0; w2 < ws->num_rrs_waves; w2++)
	u_bot[w2] += f_SGH * rhoseagrass[w2];
    }
    /*
     * Deep seagrass
     */
    if (ws->SGD_N_i > -1) {
      double SGD_N = y[ws->SGD_N_i];
    
      /* f_SGD = (1 - f_MA - f_SG -f_SGH) (1 - exp(omega_SGD * SGD)) */
      f_SGD = (1.0 - f_MA - f_SG - f_SGH) * (1.0 - exp(-SGD_N * ws->SGDorient * ws->SGDleafden));
      for (w2 = 0; w2 < ws->num_rrs_waves; w2++)
	u_bot[w2] += f_SGD * rhoseagrass[w2];
    }

    /*
     * Corals : Polyps, zoos then skelton
     */       
    if (ws->CS_N_i > -1){

      double CH_N = y[ws->CH_N_i];
      double CS_N = y[ws->CS_N_i];
      double CS_Chl = y[ws->CS_Chl_i];
      
      if (CS_N > 0.0){
	
	double cellnum = CS_N / (ws->CSm * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
	double cellChl = CS_Chl / (ws->CSvol * cellnum); /* cellular Chl conc */
	
	// intermediate var
	
	double *absorbance;
	
	absorbance = d_alloc_1d(num_rrs_waves);     

	if (ws->CS_Xp_i > -1){

	  double cellXp = y[ws->CS_Xp_i] / (ws->CSvol * cellnum); /* cellular Xp conc */
	  double cellXh = y[ws->CS_Xh_i] / (ws->CSvol * cellnum); /* cellular Xh conc */
	  
	/* Add photosynthetic cart to absorption */
	
	  for (w2=0; w2<num_rrs_waves; w2++){
	    absorbance[w2] = cellXh * bio->yC_diatoxanthin_rsr[w2] + cellXp * bio->yC_diadinoxanthin_rsr[w2] + cellChl * bio->yC_symbiodinium_rsr[w2];
	  }
	}else{
	  absorbwave(absorbance,cellChl, num_rrs_waves, rrs_wave); // old way.
	}
	
	aA_s_CS = d_alloc_1d(num_rrs_waves);
	
	aawave(ws->CSrad, absorbance, aA_s_CS, num_rrs_waves, rrs_wave);
	
	// free intermediate var
	d_free_1d(absorbance);

	double bb_zoo,u_zoo;

	// frac = (1.0 - f_MA - f_SG - f_SGH - f_SGD) * (1 - exp(-omega_CH * CH);

	double SA = ws->CHarea * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)); /* dimensionless */

	f_polyp = (1.0 - f_MA - f_SG - f_SGH - f_SGD) * SA;  // polyps
	
	f_zoo = min(f_polyp, M_PI*M_PI/(2.0*sqrt(3.0)) * cellnum * ws->CSrad * ws->CSrad);

	f_skel = f_polyp - f_zoo;
      
	/* zoos need to consider bb/(a+bb) */

	for (w2 = 0; w2 < ws->num_rrs_waves; w2++) {
	  bb_zoo = bio->scatfrac[w2] * ws->bbcell1 * pow(ws->CSrad*2e6, ws->bbcell2);
	  u_zoo = bb_zoo/(bb_zoo + aA_s_CS[w2]);
	  u_bot[w2] += f_zoo * u_zoo + f_skel * rhoskel[w2];
	}
	
	if (aA_s_CS != NULL){
	  d_free_1d(aA_s_CS);
	}
      }
    }
    
    /*
     * Microphytobenthos
     */
    if (ws->MPB_N_sed_i > -1) {
      
      double *absorbance;
      double *absorbance_xanth;

      double xan2chl = ws->xan2chl_MPB;

      double cellnum = y[ws->MPB_N_sed_i] / (ws->MPBm * 1000.0 * red_A_N * MW_Nitr); /* cell m-2 */
      double cellChl = y[ws->MPB_Chl_sed_i] / (ws->MPBvol * cellnum); /* cellular Chl conc */

      double u_mpb,bb_mpb;
      double *aA_s_MPB = d_alloc_1d(num_rrs_waves);
      
      absorbance       = d_alloc_1d(num_rrs_waves);
      absorbance_xanth = d_alloc_1d(num_rrs_waves);
      
      absorbwave(absorbance,cellChl, num_rrs_waves, rrs_wave);
      absorbxanth(absorbance_xanth, cellChl * xan2chl, num_rrs_waves, rrs_wave);
      
      for (w2=0; w2<num_rrs_waves; w2++){
	absorbance[w2] += absorbance_xanth[w2];
      }
      aawave(ws->MPBrad, absorbance, aA_s_MPB, num_rrs_waves,rrs_wave);
    
      /* f_MPB = min((1 - f_MA - f_SG - f_SGH - -f_SGD - f_polyps), pi*n*pi*r2_mpb/(2*sqrt(3))) */
      f_mpb = min(1.0 - f_MA - f_SG - f_SGH - f_SGD - f_polyp, M_PI*M_PI/(2.0*sqrt(3.0)) * cellnum * ws->MPBrad * ws->MPBrad * c->dz_sed);
      for (w2 = 0; w2 < ws->num_rrs_waves; w2++) {
	bb_mpb = bio->scatfrac[w2] * ws->bbcell1 * pow(ws->MPBrad*2e6, ws->bbcell2);
	u_mpb = bb_mpb/(bb_mpb + aA_s_MPB[w2]);
	u_bot[w2] += f_mpb * u_mpb;
      }

      // free intermediate vars
      d_free_1d(absorbance);
      d_free_1d(absorbance_xanth);
      d_free_1d(aA_s_MPB);
    }

    frac = 1.0 - f_MA - f_SG - f_SGH - f_SGD - f_polyp - f_mpb;

    /* 
     * Sediment particulates 
     */
    Sand = 0.0;
    if (ws->Sand_sed_i > -1)
      Sand = y[ws->Sand_sed_i];
    
    Mud = 0.0;
    if (ws->Mud_sed_i > -1)
      Mud = y[ws->Mud_sed_i];

    Finesed = 0.0;
    if (ws->Finesed_sed_i > -1)
      Finesed = y[ws->Finesed_sed_i];

    Dust = 0.0;
    if (ws->Dust_sed_i > -1)
      Dust = y[ws->Dust_sed_i];

    CarbSand = 0.0;
    if (ws->CarbSand_sed_i > -1)
      CarbSand = y[ws->CarbSand_sed_i];

    /* add the mineral / carbonate fractions to the Mud and CarbSand */


    if (ws->Sand_carbonate_sed_i > -1)
      CarbSand += y[ws->Sand_carbonate_sed_i];

    if (ws->Mud_carbonate_sed_i > -1)
      CarbSand += y[ws->Mud_carbonate_sed_i];

    if (ws->Sand_mineral_sed_i > -1)
      Mud += y[ws->Sand_mineral_sed_i];

    if (ws->Mud_carbonate_sed_i > -1)
      Mud += y[ws->Mud_mineral_sed_i];

    /* Sum them all up */
    all_sed = Sand + CarbSand + Mud + Finesed + Dust;

    /*
     * Particulate contributions
     */
    for (w2=0; w2<num_rrs_waves; w2++) {
      u_bot[w2] += frac * ((Sand+CarbSand)    * rhosand[w2]    / all_sed);
      u_bot[w2] += frac * (Mud     * rhomud[w2]     / all_sed);
      u_bot[w2] += frac * ((Finesed + Dust) * rhofinesed[w2] / all_sed);
    }
    
    for (w2=0; w2<num_rrs_waves; w2++) {
      u_surf[w2] = u_surf[w2] + u_bot[w2] * w_bot[w2];
    }

    /* Save mean bottom reflectance - problem with this calculation is it is operating on w2, not w wavelenghts. */
    
    if (ws->SWR_bot_abs_i > -1){
      y[ws->SWR_bot_abs_i] = 0.0;
      for (w2=0; w2<ws->num_rrs_waves; w2++){
	y[ws->SWR_bot_abs_i] += u_bot[w2];
      }
      y[ws->SWR_bot_abs_i] = 1.0 - (y[ws->SWR_bot_abs_i] / ws->num_rrs_waves);
    }

    /* Secchi depth */
  
    if (ws->Secchi_i > -1){
      double* z_secchi_cv = col->cv[ws->z_secchi_i];
      y[ws->Secchi_i] = z_secchi_cv[0];
    }
   
    /* Call helper routine to finish up */
    finalise_reflectances(ws, c);
    d_free_1d(u_bot);

}   // if reflectances
}

static void finalise_reflectances(workspace *ws, cell *c)
{
  double rs;
  column* col = c->col;
  double* u_surf = col->cv[ws->u_surf_i];
  double* y = c->y;
  int w2;

  for (w2=0; w2<ws->num_rrs_waves; w2++) {
     rs = u_surf[w2] * (ws->g0 + ws->g1 * u_surf[w2]);
    //rs = u_surf[w2] / 4.0 / 3.145926535;
    y[ws->RRS_i[w2]] = 0.52 * rs / (1.0 - 1.7 * rs);
  }

  /* OC3M algorithm */
  
  if (ws->OC3M_i > -1){
    
    double ratio;
    double a_oc3[5] = {0.283, -2.753, 1.457, 0.659, -1.403};
    
    if (y[ws->RRS_i[ws->w443]] > y[ws->RRS_i[ws->w488]]) 
      ratio = log10(y[ws->RRS_i[ws->w443]]/y[ws->RRS_i[ws->w547]]);
    else
      ratio = log10(y[ws->RRS_i[ws->w488]]/y[ws->RRS_i[ws->w547]]);
    
    y[ws->OC3M_i] = pow(10.0,a_oc3[0] + ratio *(a_oc3[1] + ratio *(a_oc3[2] + ratio *(a_oc3[3] + ratio * a_oc3[4]))));
  }

  /* OC4Me algorithm used for both MERIS and Sentinal-3 */
  
  if (ws->OC4Me_i > -1){
    
    double ratio;
    double a_oc4[5] = {0.3255, -2.7677, 2.4409, -1.1288, -0.4990};

    if (y[ws->RRS_i[ws->w443]] > y[ws->RRS_i[ws->w488]] && y[ws->RRS_i[ws->w443]] > y[ws->RRS_i[ws->w510]]){
      ratio = log10(y[ws->RRS_i[ws->w443]]/y[ws->RRS_i[ws->w560]]);
    }else{
      if (y[ws->RRS_i[ws->w488]] > y[ws->RRS_i[ws->w510]])
	ratio = log10(y[ws->RRS_i[ws->w488]]/y[ws->RRS_i[ws->w560]]);
      else
	ratio = log10(y[ws->RRS_i[ws->w510]]/y[ws->RRS_i[ws->w560]]);
    }
    y[ws->OC4Me_i] = pow(10.0,a_oc4[0] + ratio *(a_oc4[1] + ratio *(a_oc4[2] + ratio *(a_oc4[3] + ratio * a_oc4[4]))));
  }

  if (ws->Hue_i > -1){     // Hue angle algorithm used for Sentinal-3.
   
    double x_tris_MSI = 11.756*y[ws->RRS_i[ws->w443]] + 6.423*y[ws->RRS_i[ws->w490]] + 53.696*y[ws->RRS_i[ws->w560]] + 32.028*y[ws->RRS_i[ws->w665]] + 0.529*y[ws->RRS_i[ws->w709]];
    double y_tris_MSI = 1.744*y[ws->RRS_i[ws->w443]] + 22.289*y[ws->RRS_i[ws->w490]] + 65.702*y[ws->RRS_i[ws->w560]] + 16.808*y[ws->RRS_i[ws->w665]] + 0.192*y[ws->RRS_i[ws->w709]];
    double z_tris_MSI = 62.696*y[ws->RRS_i[ws->w443]] + 31.101*y[ws->RRS_i[ws->w490]] + 1.778*y[ws->RRS_i[ws->w560]] + 0.015*y[ws->RRS_i[ws->w665]] + 0.000*y[ws->RRS_i[ws->w709]];

    double cx_MSI = x_tris_MSI/(x_tris_MSI+y_tris_MSI+z_tris_MSI+1.0e-6);
    double cy_MSI = y_tris_MSI/(x_tris_MSI+y_tris_MSI+z_tris_MSI+1.0e-6);
    
    double xw_MSI = cx_MSI - 0.333333333;
    double yw_MSI = cy_MSI - 0.333333333;

    double atanterm = atan(yw_MSI/(xw_MSI+1.0e-6));
    
    y[ws->Hue_i] = atanterm*180.0/M_PI;

    if (xw_MSI<0.0){
      if (yw_MSI<0.0){
	y[ws->Hue_i] = 180.0 + atanterm*180.0/M_PI;
      }else{
	y[ws->Hue_i] = 180.0 - atanterm*180.0/M_PI;
      }
    }
  }
  /* MODIS nFLH */
  
  /* nFLH uses nLw, which is normalised (to zero zenith) water leaving irradiance 
     for MODIS 678, and has units mW / cm2 / um / sr-1. This includes a factor fo = 148.097 mW / cm2 / um
     pi - sr-1 */
  
  if (ws->nFLH_i > -1){
    y[ws->nFLH_i] = (148.097 / 3.145926535) * (y[ws->RRS_i[ws->w678]] - (70.0/81.0) * y[ws->RRS_i[ws->w667]] - (11.0/81.0) * y[ws->RRS_i[ws->w748]]);
  }

  /* TSS algorithm - use local coastal relationship from Petus et al., 2014 */
  
  if (ws->TSSM_i > -1){
    double r645 = y[ws->RRS_i[ws->w645]];
    y[ws->TSSM_i] = (12450.0*r645*r645 + 666.0*r645 + 0.48)/1000.0;
  }
  /* KD490 algorithm - use NASA  */
  
  if (ws->KD490M_i > -1){

    double ratio1 = y[ws->RRS_i[ws->w488]]/y[ws->RRS_i[ws->w547]];

    double X = log10(ratio1);

    double a_kd2m[5] = {-0.8813, -2.0584, 2.5878, -3.4885, -1.5061};
 
    y[ws->KD490M_i] = pow(10.0,a_kd2m[0] + X * (a_kd2m[1] + X * (a_kd2m[2] + X * (a_kd2m[3] + X * a_kd2m[4])))) + 0.0166;

  }
  if (ws->OC3V_i > -1){
    
    double ratio2;
    double a_oc3[5] = {0.2228,-2.4683,1.5867,-0.4275,-0.7768};
    
    if (y[ws->RRS_i[ws->w443]] > y[ws->RRS_i[ws->w486]]) 
      ratio2 = log10(y[ws->RRS_i[ws->w443]]/y[ws->RRS_i[ws->w551]]);
    else
      ratio2 = log10(y[ws->RRS_i[ws->w486]]/y[ws->RRS_i[ws->w551]]);
    
    y[ws->OC3V_i] = pow(10.0,a_oc3[0] + ratio2 *(a_oc3[1] + ratio2 *(a_oc3[2] + ratio2 *(a_oc3[3] + ratio2 * a_oc3[4]))));
  }

}
