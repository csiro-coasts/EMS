/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_epi.c
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
 *  $Id: light_spectral_epi.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "light_spectral_epi.h"
#include "constants.h"

typedef struct {
  /*
   * parameters
   */
  double CSm;
  double CSrad;
  double CSvol;

  /* epibenthic variables */

  int Epilight_s_i;
  int Epilight_i;

  int num_waves;
  double *wave;

  /* tracers */

  int MA_N_i;
  int SG_N_i;
  int SGH_N_i;
  int CH_N_i;
  int CS_N_i;
  int CS_Chl_i;

  /* common variables */

  int cv_lighttop_s_i;
  int lighttop_i;
  int KI_SG_i;
  int KI_SGH_i;
  int KI_MA_i;
  int KI_CS_i;
  int yCfac_CS_i;

  int wPAR_bot;
  int wPAR_top;

  /* Local variables */

  double MAleafden;
  double SGleafden;
  double SGHleafden;
  double CHpolypden;
  double CHarea;
  double *aA_s_CS;

  /*
   * output variables
   */

  int EpiPAR_i;

  // int Epilight_novar;
  // int *Epilight_ovar_s_i;
  // int *Epilight_ovar_i;

} workspace;

void light_spectral_epi_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* tracers = e->tracers;
  stringtable* epis = e->epis;
  workspace* ws = malloc(sizeof(workspace));


  
  int OFFSET_EPI = tracers->n * 2;
  
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

  ws->Epilight_i = e->try_index(epis, "Epilight", e);
   if (ws->Epilight_i >= 0)
     ws->Epilight_i += OFFSET_EPI;

/* Get wavelengths for light */

  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");
  
  /*  common column variables */

  ws->cv_lighttop_s_i = find_index_or_add_vector(e->cv_column, "lighttop_s", e, ws->num_waves);
  ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);

  /*
   * Cache indices and wavelengths of the variables to output
   */
  
 /*  ws->Epilight_novar = 0; */
  
/*   for (w=0; w<ws->num_waves; w++) { */
/*     int index; */
/*     // Epilight_s */
/*     sprintf(buf, "Epilight_%d", (int)ws->wave[w]); */
/*     index = e->try_index(epis, buf, e); */
/*     if (index >= 0) { */
/*       ws->Epilight_novar++; */
/*     } */
/*   } */
  
/*   if (ws->Epilight_novar > 0) { */
/*     ws->Epilight_ovar_s_i = i_alloc_1d(ws->Epilight_novar); */
/*     ws->Epilight_ovar_i   = i_alloc_1d(ws->Epilight_novar); */
/*   } */
  
/*   count_Li = 0; */
/*   for (w=0; w<ws->num_waves; w++) { */
/*     // Light */
/*     int index; */
/*     sprintf(buf, "Epilight_%d", (int)ws->wave[w]); */
/*     index = e->try_index(epis, buf, e); */
    
/*     if (index >= 0) { */
/*       // Found it! */
/*       ws->Epilight_ovar_s_i[count_Li] = w; */
/*       ws->Epilight_ovar_i[count_Li]   = index + OFFSET_EPI; */
/*       count_Li++; */
/*     } */
/*   } */
}

void light_spectral_epi_postinit(eprocess* p)
{
  double *absorbance;
  ecology* e = p->ecology;


  workspace* ws = (workspace *)p->workspace;
  double NtoCHL = get_parameter_value(e, "NtoCHL");



  int w;

/* set light indexes so that code isn't looking for them each time step. */

  for (w=0; w<ws->num_waves; w++){
    if (ws->wave[w] == 410.0)
      ws->wPAR_bot = w;
    if (ws->wave[w] == 690.0)
      ws->wPAR_top = w;
  }
  
  if (ws->MA_N_i > -1){
    
    /* Put in macroalgae initialisations if there is a macroalgae tracer */

    ws->MAleafden = get_parameter_value(e, "MAleafden");
    ws->KI_MA_i = find_index_or_add(e->cv_cell, "KI_MA", e);

  }

  if (ws->SG_N_i > -1){

    /* Put in Zostera initialisations if there is a Zostera tracer */

    ws->SGleafden = get_parameter_value(e, "SGleafden");
    ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SG", e);

  }

  if (ws->SGH_N_i > -1){

    /* Put in Halophila initialisations if there is a Halophila tracer */

    ws->SGHleafden = get_parameter_value(e, "SGHleafden");
    ws->KI_SGH_i = find_index_or_add(e->cv_cell, "KI_SGH", e);

  }

  if (ws->CS_N_i > -1){

    /* Only in coral initialisations if there is a coral tracer */

    ws->CHarea = try_parameter_value(e, "CHarea");
    
    ws->CHpolypden = get_parameter_value(e, "CHpolypden");
    ws->CSrad = get_parameter_value(e, "CSrad");
    ws->CSm = PhyCellMass(ws->CSrad);
    ws->CSvol = 4.0 / 3.0 * M_PI * ws->CSrad * ws->CSrad * ws->CSrad;
    
    ws->KI_CS_i = find_index_or_add(e->cv_cell, "KI_CS", e);
    ws->yCfac_CS_i = find_index_or_add(e->cv_cell, "yCfac_CS", e);
    
    ws->aA_s_CS = NULL;
    
    if (ws->KI_CS_i > -1)
      ws->aA_s_CS = d_alloc_1d(ws->num_waves);
    
    // intermediate var
    absorbance = d_alloc_1d(ws->num_waves);
    
    if (ws->aA_s_CS != NULL) {
      absorbwave(absorbance, ws->CSm*6.625*14.0/(ws->CSvol*NtoCHL), ws->num_waves, ws->wave);
      aawave(ws->CSrad, absorbance, ws->aA_s_CS, ws->num_waves, ws->wave);
    }
    
    // free intermediate var
    d_free_1d(absorbance);
  }
  
}

void light_spectral_epi_destroy(eprocess* p)
{
  workspace* ws = p->workspace;
  
  if (ws->aA_s_CS != NULL)
    d_free_1d(ws->aA_s_CS);
  free(ws);
}

void light_spectral_epi_precalc(eprocess* p, void* pp)
{
  ecology* e = p->ecology;
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  double* lighttop_s = col->cv[ws->cv_lighttop_s_i];

  int num_waves = ws->num_waves;
  int w;

  double *wave  = ws->wave;
  bio_opt_prop *bio = e->bio_opt;

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
    
    double MA_N = y[ws->MA_N_i];
    c->cv[ws->KI_MA_i] = 0.0;
    
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

  /* Zostera first */

  if (ws->SG_N_i > -1){

    double SG_N = y[ws->SG_N_i];
    c->cv[ws->KI_SG_i] = 0.0;

    if (SG_N > 0.0){  /* seagrass */
      
      photons = 0.0;
      
      for (w=0; w<num_waves; w++){
	photons += lighttop_s[w] * (1.0 - exp(-SG_N * ws->SGleafden * bio->SG_aAwave[w])) 
	  * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SG_N * ws->SGleafden * bio->SG_aAwave[w]);  
      }
      c->cv[ws->KI_SG_i] = photons;
    }
  }

  /* Then Halophila if present */

  if (ws->SGH_N_i > -1){

    double SGH_N = y[ws->SGH_N_i];
    c->cv[ws->KI_SGH_i] = 0.0;

    if (SGH_N > 0.0){  /* seagrass */

      photons = 0.0;

      for (w=0; w<num_waves; w++){

	photons += lighttop_s[w] * (1.0 - exp(-SGH_N * ws->SGHleafden * bio->SGH_aAwave[w])) 
                                            * 8.359335857479461e-09 * wave[w];
	lighttop_s[w] = lighttop_s[w] * exp(-SGH_N * ws->SGHleafden * bio->SGH_aAwave[w]);
      }
      c->cv[ws->KI_SGH_i] = photons;
    }
  }
  
  if (ws->CS_N_i > -1){

    double CH_N = y[ws->CH_N_i];
    double CS_N = y[ws->CS_N_i];
    double CS_Chl = y[ws->CS_Chl_i];

    c->cv[ws->KI_CS_i] = 0.0;
    c->cv[ws->yCfac_CS_i] = 0.0;

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
      absorbwave(absorbance,cellChl, num_waves, wave);	
      
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
      
      aawave(ws->CSrad, absorbance, ws->aA_s_CS, num_waves,wave);
      
      // free intermediate var
      d_free_1d(absorbance);

      /* Integrate photons absorbed across all wavelengths */

      photons = 0.0;
      
      if (isnan(ws->CHarea))
	ws->CHarea = 1.0;

      double RR;
      double area;

      if (process_present(e,PT_WC,"recom_extras")){
	
	/* do calculation based on area on present cell dimension */
	
	area = einterface_cellarea(c->col->model,c->b); 
	
	RR = sqrt(area/PI);
	ws->CHarea = 1.0;
	if (RR > 200.0)
	  ws->CHarea = 1.0-(RR-200.0)*(RR-200.0)/(RR*RR);
      }
      

      /* establish surface area of coral polyps - only attenuate within this area */ 

      double SA = ws->CHarea * (1.0-exp(-CH_N * ws->CHpolypden / ws->CHarea)); /* dimensionless */

      for (w=0; w<num_waves; w++){

	/* 1e-12 - just very low light to avoid dividing by zero. */ 

	kd = max(ws->aA_s_CS[w] * cellnum / SA,1e-12);

        Itop = lighttop_s[w];

        /* use CHarea to establish what misses coral communities - also above needed to divide cellnum to account for squeezing. */ 

	lighttop_s[w] = (1.0 - SA) * lighttop_s[w] + SA * lighttop_s[w] * exp(-kd);
	
	// meanlight = (Itop - lighttop_s[w]) / kd;

        meanlight = Itop * (1.0 - exp(-kd)) / kd;

	photons += meanlight * 8.359335857479461e-09 * wave[w] * ws->aA_s_CS[w];

      }    
      c->cv[ws->KI_CS_i] = photons;
    }
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
  
  lighttop = lighttop/wavesum; 
  cv_lighttop[0] = lighttop;
  
  // Assign light above epi-pelagic to output. This is simply bottom of water column, but this 
  // is not stored in the water column (the average light in layer is).
  
  //  for (i=0; i<ws->Epilight_novar; i++) {
  //  y[ws->Epilight_ovar_i[i]] = lighttop_s[ws->Epilight_ovar_s_i[i]];
  

  // ** KWA - add following to save non-spectral Epilight (I think this is bottom of epi layer light

  if (ws->EpiPAR_i > -1) {
    y[ws->EpiPAR_i] = lighttop;
  }
}
