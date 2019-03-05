/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral2par_epi.c
 *  
 *  OBSOLETE: not necessary if using light_spectral_sed.c
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: light_spectral2par_epi.c 6046 2018-12-12 00:53:28Z bai155 $
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
#include "light_spectral2par_epi.h"
#include "constants.h"

#define EPS_KD 1.0e-10

typedef struct {

  /* epibenthic variables */

  int Epilight_i;
  int Epilightatt_i;

  int num_waves;
  double *wave;

  /* common variables */

  int lighttop_i;
  int lighttop_s_i;
  int cv_lighttop_s_i;

} workspace;

void light_spectral2par_epi_init(eprocess* p)
{
  ecology* e = p->ecology;
  stringtable* epis = e->epis;
  workspace* ws = malloc(sizeof(workspace));
  int OFFSET_EPI = e->ntr * 2;
  
  p->workspace = ws;

  /* spectrally-resolved parameters */
  
  ws->num_waves = get_parameter_num_values(e, "Light_lambda");
  ws->wave      = get_parameter_value_ptr(e, "Light_lambda");
  
  /* epis */
  
  ws->Epilight_i = e->try_index(epis, "Epilight", e) + OFFSET_EPI;
  if (ws->Epilight_i >= 0)
    ws->Epilightatt_i = e->try_index(epis, "Epilightatt", e) + OFFSET_EPI;
  else
    ws->Epilightatt_i = -1;
  
  /*  common column variables */

  ws->cv_lighttop_s_i = find_index_or_add_vector(e->cv_column, "lighttop_s", e,
					                       ws->num_waves);
  ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);
}

void light_spectral2par_epi_destroy(eprocess* p)
{
  free(p->workspace);
}

void light_spectral2par_epi_precalc(eprocess* p, void* pp)
{
  cell* c = (cell*) pp;
  column* col = c->col;
  workspace* ws = p->workspace;
  double* y = c->y;
  double lighttop;
  double *cv_lighttop = col->cv[ws->lighttop_i];
  double *lighttop_s = col->cv[ws->cv_lighttop_s_i];

  /* sum up light in the PAR range (400-700 nm) */

  int w;
  double wavesum;
  double *wave = ws->wave;

  lighttop = 0.0;
  wavesum = 0.0;

 /* sum up light in the PAR range (400-700 nm), 
     assuming weighting by photons (not energy)  */

  for (w=0; w<ws->num_waves; w++){
    if (wave[w] >= 400.0 & wave[w]<= 700.0){
      lighttop = lighttop + lighttop_s[w]*wave[w];
      wavesum = wavesum + wave[w];
    }
  }
  lighttop = lighttop/wavesum; 

  /* save the light intensity at the top of the epibenthos */

  if (ws->Epilight_i >= 0) {
    y[ws->Epilight_i] = lighttop;
    
    /* calculate light attenuation by Macroalgae and Seagrass */
    
    /*UR-what if we don't have Macroalgea and Seagrass */
    if (ws->Epilightatt_i >= 0)
      lighttop = lighttop * exp(-y[ws->Epilightatt_i]);
  }
  
  cv_lighttop[0] = lighttop;
}
