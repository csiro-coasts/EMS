/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/coral_spectral_carb_epi.c
 *  
 *  Description:
 *  Just coral calcification and dissolution terms.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: coral_spectral_carb_epi.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "constants.h"
#include "coral_spectral_carb_epi.h"

typedef struct {
  /*
   * parameters
   */

  double CHpolypden;
  double MAleafden;
  double CSrad;
  double CSm;
  double CHarea;

  int recom;

  /*
   * epis
   */

  int CH_N_i;    /* Coral host - animal */
  int CS_N_i;    /* Coral symbiont - microalgae */
  int CS_I_i;

  int MA_N_i;
  int ustrcw_skin_i;

  /*
   * tracers
   */
  
  int temp_wc_i;
  int DIC_wc_i;

  int KI_CS_i;

  int omega_ar_i;
  int ALK_wc_i;
  int Gnet_i;

  double k_day_coral;
  double k_night_coral;

  double dissCaCO3_reef;
  double dissCaCO3_shelf;
  
  /*
   * common cell variables
   */

  int Tfactor_i;
  int CHarea_cv_i;
 
} workspace;

void coral_spectral_carb_epi_init(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = malloc(sizeof(workspace));
    stringtable* tracers = e->tracers;
    stringtable* epis = e->epis;

    int OFFSET_EPI = tracers->n * 2;

    p->workspace = ws;

    /*
     * parameters
     */

    ws->CHpolypden = get_parameter_value(e, "CHpolypden");
    ws->CHarea = try_parameter_value(e, "CHarea"); 
    ws->CSrad = get_parameter_value(e, "CSrad");
    ws->CSm = PhyCellMass(ws->CSrad);

    // Calculated assuming absorption cross section is pi r^2 m cell-1, and Suggett 2008 measured absorption cross section of PSII per photon.

    /*
     * epis
     */

    ws->CH_N_i = e->find_index(epis, "CH_N", e) + OFFSET_EPI;
    ws->CS_N_i = e->find_index(epis, "CS_N", e) + OFFSET_EPI;
    ws->CS_I_i = e->find_index(epis, "CS_I", e) + OFFSET_EPI;

    ws->MA_N_i = e->try_index(epis, "MA_N", e);

    if (ws->MA_N_i > -1){
	 ws->MA_N_i += OFFSET_EPI;
	 ws->MAleafden = get_parameter_value(e, "MAleafden");
    }

    ws->Gnet_i = e->find_index(epis, "Gnet", e) + OFFSET_EPI;
    ws->ustrcw_skin_i = e->find_index(epis, "ustrcw_skin", e) + OFFSET_EPI;

    /*
     * tracers
     */

    ws->temp_wc_i = e->find_index(tracers, "temp", e);
    ws->DIC_wc_i = e->find_index(tracers, "DIC", e);

    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->CHarea_cv_i = find_index_or_add(e->cv_cell, "CHarea_cv", e);

    /* Calcification */

    ws->omega_ar_i = e->find_index(tracers,"omega_ar", e);
    ws->ALK_wc_i = e->find_index(tracers, "alk", e);

    if (ws->omega_ar_i > -1){
       ws->k_day_coral   = get_parameter_value(e, "k_day_coral");
       ws->k_night_coral = get_parameter_value(e, "k_night_coral");
       ws->dissCaCO3_reef = get_parameter_value(e, "dissCaCO3_sed");

       ws->dissCaCO3_shelf = try_parameter_value(e, "dissCaCO3_shelf");
       if (isnan(ws->dissCaCO3_shelf)){
	   ws->dissCaCO3_shelf = 0.0001;
	 }
    }
}

void coral_spectral_carb_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->KI_CS_i = find_index_or_add(e->cv_cell, "KI_CS", e);
    
    ws->recom = 0;
    if (process_present(e,PT_WC,"recom_extras"))
	ws->recom = 1;
}

void coral_spectral_carb_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;
    free(ws);
}

void coral_spectral_carb_epi_precalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double area,RR;
    
    cv[ws->CHarea_cv_i] = ws->CHarea; // get from bio.prm file.

    if (ws->recom){
      
      /* do calculation based on area on present cell dimension: This should be a column variable? */
      
      area = einterface_cellarea(c->col->model,c->b); 
      cv[ws->CHarea_cv_i] = 1.0;
      RR = sqrt(area/PI);
      if (RR > 200.0)
	cv[ws->CHarea_cv_i] = 1.0-(RR-200.0)*(RR-200.0)/(RR*RR);
    }

    /* Coupling with physical processes (temperature, shear stress, light) is calculated only 
       at the beginning of the ecological time increment */

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
}

void coral_spectral_carb_epi_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;

    double* cv = c->cv;
    double dz_wc = c->dz_wc;

    double CH_N = y[ws->CH_N_i];
    double CS_N = y[ws->CS_N_i];

    //double dissCaCO3_shelf = (-11.51 * y[ws->omega_ar_i] + 33.683)/86400.0;
    double  dissCaCO3_shelf=ws->dissCaCO3_shelf;
    // printf("ws->omega_ar_i %d, CH_N %e, y[ws->omega_ar_i] %e \n",ws->omega_ar_i,CH_N, y[ws->omega_ar_i]);


    if ((CS_N < 1e-12)||(CH_N < 1e-9)){

      /* Still need to dissolve carbonate sands if on shelf */
      
      y1[ws->DIC_wc_i] += 12.01 * dissCaCO3_shelf / dz_wc ;
      y1[ws->ALK_wc_i] += 2.0 *  dissCaCO3_shelf / dz_wc ;
      if (ws->Gnet_i> -1)
	y1[ws->Gnet_i] -= 12.01 * dissCaCO3_shelf ;  
      return;
    }
    
    double temp_wc = y[ws->temp_wc_i];
    double CS_I = y[ws->CS_I_i];
    double PI_max = ws->CSm * red_A_I * 1000.0 ; /* mmol photon cell -1 */
    double Iquota;
    double MA_N = 0.0;

    double CHarea = cv[ws->CHarea_cv_i];

    if (ws->MA_N_i > -1){
      MA_N = y[ws->MA_N_i];
    }

    /* Zoothanthellae - cells quantified per m2 of host tissue. */

    double cellnum = CS_N / (ws->CSm * red_A_N * 1000.0 * MW_Nitr); /* cell m-2 */
 
    /* now get normalised reserves. */

    // if (ws->omega_ar_i > -1){
      Iquota = max(0.0,min(1.0,(CS_I / cellnum) / PI_max));    /* bleach form */
      // }else{
      // Iquota = min(1.0,kI/min(kN_mass, kP_mass));
      // Iquota = min(Iquota,kI/(CS_N * mgN2molN)/CSumax);
      // Iquota  = max(Iquota,0.0); 
      // }
    /* Information from light_spectral_epi.c */
    
    double kI = c->cv[ws->KI_CS_i]; // mol photon cell-1 s-1 */

    double SA = CHarea * exp(-MA_N * ws->MAleafden) * (1.0-exp(-CH_N * ws->CHpolypden / CHarea)); /* dimensionless */

    // printf("CHarea %e, MA_N %e, ws->MAleafden %e, CH_N %e, ws->CHpolypden %e,\n",CHarea,MA_N,ws->MAleafden,CH_N,ws->CHpolypden);

    /* Now do calcification and dissolution of sand among corals */

    if (ws->omega_ar_i > -1 && y[ws->omega_ar_i] > 1.0){

      /* Net community calcification - i.e. calc - diss of community */;
      
      double g_night_coral = ws->k_night_coral*(y[ws->omega_ar_i]-1.0);
      double g_day_coral   = ws->k_day_coral*(y[ws->omega_ar_i]-1.0);
      
      if (c->cv[ws->KI_CS_i] > 1e-16){
	g_night_coral = 0.0;
      }else{
	g_day_coral = 0.0;
      }
      
      /* Similar to Anthony 2011, except use Iquota rather than SWR */
      
      //      double g0_coral = (g_day_coral * Iquota  + g_night_coral) * SA;

      double g0_coral = (g_day_coral * Iquota * Iquota  + g_night_coral) * SA;

      /* dissolution 7 mmol m-2 d-1 Cryonak LO, convert to seconds */

      double dissCaCO3_reef = ws->dissCaCO3_reef;
      
      /* Now do the effect of calcification and dissolution on water column carbon chemistry */
      /* Use SA as a proxy of percent coral coverage */
      
      /* Alkalinity is in mmol m-3, DIC is in mg C m-3 */
      
      y1[ws->DIC_wc_i] -= 12.01 * (g0_coral - dissCaCO3_reef ) / dz_wc ;
      y1[ws->ALK_wc_i] -= 2.0 * (g0_coral - dissCaCO3_reef  )/ dz_wc ;
     
      /* output calcification rate - mg C m-2 s-1 - so it can be used in mass balance */

      if (ws->Gnet_i> -1)
	y1[ws->Gnet_i] += 12.01 * (g0_coral - dissCaCO3_reef) ;
    }
}

void coral_spectral_carb_epi_postcalc(eprocess* p, void* pp)
{
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  double* y = c->y;

}
