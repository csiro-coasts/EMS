/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_spectral_grow_epi.c
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
 *  $Id: seagrass_spectral_grow_epi.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "seagrass_spectral_grow_epi.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20
#define unitch 1000.0

typedef struct {
  int do_mb;                  /* flags */
  char species;
  
  /*
   * parameters
   */
  double umax_t0;
  double omega;
  double KN;
  double KP;
  
  double SGfrac;
  double SGtransrate;
  
  /*
   * epis
   */
  int SG_N_i;
  int SG_N_pr_i;
  int SG_N_gr_i;
  int EpiTN_i;
  int EpiTP_i;
  int EpiTC_i;
  int EpiBOD_i;
  
  int KI_SG_i;
  int SGROOT_N_i;
  
  /*
   * tracers
   */
  int NH4_sed_i;
  int NO3_sed_i;
  int DIP_sed_i;
  int DIC_sed_i;
  int Oxygen_wc_i;
  int Oxy_pr_wc_i;
  
  /*
   * common cell variables
   */
  int Tfactor_i;
  int umax_i;

} workspace;

void seagrass_spectral_grow_epi_init(eprocess* p)
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

    if (ws->species == 'Z'){   /* Zostera - standard terms */

      printf("Seagrass species Zostera chosen \n");

      ws->umax_t0 = get_parameter_value(e, "SGumax");
      ws->omega = get_parameter_value(e, "SGleafden");
      ws->KN = get_parameter_value(e, "SG_KN");
      ws->KP = get_parameter_value(e, "SG_KP");
      ws->SGROOT_N_i = e->find_index(epis, "SGROOT_N", e) + OFFSET_EPI;
      ws->SGfrac = get_parameter_value(e, "SGfrac");
      ws->SGtransrate = get_parameter_value(e, "SGtransrate");

      ws->SG_N_i = e->find_index(epis, "SG_N", e) + OFFSET_EPI;

      ws->umax_i = find_index_or_add(e->cv_cell, "SGumax", e);

      /*non essential diagnostic tracers*/

      ws->SG_N_pr_i = e->try_index(epis, "SG_N_pr", e);
      if (ws->SG_N_pr_i > -1) 
	ws->SG_N_pr_i += OFFSET_EPI;

      ws->SG_N_gr_i = e->try_index(epis, "SG_N_gr", e);
      if (ws->SG_N_gr_i > -1) 
	ws->SG_N_gr_i += OFFSET_EPI;

    } else if (ws->species == 'H'){   /* Halophila */

      /* to avoid changing code, SGH terms are assigned to SG */

      printf("Seagrass species Halophila chosen \n");

      ws->umax_t0 = get_parameter_value(e, "SGHumax");
      ws->omega = get_parameter_value(e, "SGHleafden");
      ws->KN = get_parameter_value(e, "SGH_KN");
      ws->KP = get_parameter_value(e, "SGH_KP");
      ws->SGROOT_N_i = e->find_index(epis, "SGHROOT_N", e) + OFFSET_EPI;
      ws->SGfrac = get_parameter_value(e, "SGHfrac");
      ws->SGtransrate = get_parameter_value(e, "SGHtransrate");

      ws->SG_N_i = e->find_index(epis, "SGH_N", e) + OFFSET_EPI;

      ws->umax_i = find_index_or_add(e->cv_cell, "SGHumax", e);

      /*non essential diagnostic tracers*/

      ws->SG_N_pr_i = e->try_index(epis, "SGH_N_pr", e);
      if (ws->SG_N_pr_i > -1) 
	ws->SG_N_pr_i += OFFSET_EPI;

      ws->SG_N_gr_i = e->try_index(epis, "SGH_N_gr", e);
      if (ws->SG_N_gr_i > -1) 
	ws->SG_N_gr_i += OFFSET_EPI;
      
    } else if (ws->species == 'P'){   /* Posidonia */
 
      /* still to write posidonia code */


    } else{
      e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected seagrass type: expected \"Zostera\", \"Halophila\" or \"Posidonia\"\n", e->processfname, p->name, prm);

    }

    ws->EpiTN_i = e->find_index(epis, "EpiTN", e) + OFFSET_EPI;
    ws->EpiTP_i = e->find_index(epis, "EpiTP", e) + OFFSET_EPI;
    ws->EpiTC_i = e->find_index(epis, "EpiTC", e) + OFFSET_EPI;
    ws->EpiBOD_i = e->find_index(epis, "EpiBOD", e) + OFFSET_EPI;
    /*
     * tracers
     */
    ws->Oxygen_wc_i = e->find_index(tracers, "Oxygen", e);
    ws->NH4_sed_i = e->find_index(tracers, "NH4", e) + OFFSET_SED;
    ws->NO3_sed_i = e->find_index(tracers, "NO3", e) + OFFSET_SED;
    ws->DIP_sed_i = e->find_index(tracers, "DIP", e) + OFFSET_SED;
    ws->DIC_sed_i = e->find_index(tracers, "DIC", e) + OFFSET_SED;
   
    /*non essential diagnostic tracer*/

    ws->Oxy_pr_wc_i = e->try_index(tracers, "Oxy_pr", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
}

void seagrass_spectral_grow_epi_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_epi", e) >= 0) ? 1 : 0;

    /* In this routine KI_SG is generic, but species-specific column variables 
       must be added for light_spectral_epi.c */

    if (ws->species == 'Z'){   /* Zostera - standard terms */

      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SG", e);

    } else if (ws->species == 'H'){   /* Halophila */

      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SGH", e);

    } else if (ws->species == 'P'){   /* Posidonia */

      ws->KI_SG_i = find_index_or_add(e->cv_cell, "KI_SGP", e);

    }
}

void seagrass_spectral_grow_epi_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws);
}

void seagrass_spectral_grow_epi_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;



    double SG_N = y[ws->SG_N_i] * unitch;
    double SGROOT_N = y[ws->SGROOT_N_i] * unitch;

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    // printf("In SG: ws->umax_i %d, ws->umax_t0 %e, cv[ws->umax_i] %e, Tfactor %e \n",ws->umax_i,ws->umax_t0,cv[ws->umax_i],Tfactor);

    if (ws->do_mb) {
      y[ws->EpiTN_i] += SG_N + SGROOT_N;
      y[ws->EpiTP_i] += (SG_N + SGROOT_N)* atk_W_P;
      y[ws->EpiTC_i] += (SG_N + SGROOT_N)* atk_W_C;
      y[ws->EpiBOD_i] += (SG_N + SGROOT_N)* atk_W_O;
    }
}

void seagrass_spectral_grow_epi_calc(eprocess* p, void* pp)
{
  workspace* ws = p->workspace;
  intargs* ia = (intargs*) pp;
  double* y = ia->y;
  double* y1 = ia->y1;
  cell* c = (cell*) ia->media;
  double* cv = c->cv;
  double dz_sed = c->dz_sed;
  double dz_wc = c->dz_wc;
  double porosity = c->porosity;
  
  double SG_N = y[ws->SG_N_i];
  double SGROOT_N = y[ws->SGROOT_N_i];
  
  if (SG_N > 1e-12){

    double NO3_sed = y[ws->NO3_sed_i];
    double NH4_sed = y[ws->NH4_sed_i];
    double din_sed = NO3_sed + NH4_sed;
    double DIN_sed = (din_sed > 0.0) ? din_sed : EPS_DIN;
    double dip_sed = y[ws->DIP_sed_i];
    double DIP_sed = (dip_sed > 0.0) ? dip_sed : EPS_DIP;

    double umax = cv[ws->umax_i]; /* s-1 */

    double SGfrac = ws->SGfrac;
    double SGtransrate = ws->SGtransrate;

    /* kI is calculated in light_spectral_epi.c in mol photon m-2 s-1 and converted to turnover time */

    /* SG_N now in g N m-2, so need to multiply mgN2molN by 1000 */

    double kI = c->cv[ws->KI_SG_i] * atk_A_N / ( atk_A_I * SG_N * mgN2molN * 1000.0);

    /* Respiration up to 4 mol photon m-2 d-1, 0.7 is the wavelength average absorbance */

//    double resp = (30.0/5500.0)*14.0*(4.0*0.7/86400.0)*ws->omega;

//    kI = max(0.0,kI-resp);

    /* Uptake from sediment is done using a Michealis-Menton form */

    /* Modify half-sat constant to account for root reach */

    double kN = umax * DIN_sed / ( ws->KN / ( SGROOT_N / (SG_N + SGROOT_N)) + DIN_sed);
    double kP = umax * DIP_sed / ( ws->KP / ( SGROOT_N / (SG_N + SGROOT_N)) + DIP_sed);

    // Law of the minimum, but including maximum growth rate.

    double growthrate = min(kI,min(kN, kP));

    double growth = SG_N * growthrate;

    double kNH4 = umax * NH4_sed / ( ws->KN + NH4_sed);

    double NH4uptake = min(kNH4 * SG_N, growth);
    double NO3uptake = growth - NH4uptake;

    double Oxy_pr = growth * unitch * atk_W_O / dz_wc;
    
    y1[ws->SG_N_i] += growth;
    //y1[ws->NH4_sed_i] -= growth * unitch * NH4_sed / DIN_sed / dz_sed / porosity;
    //y1[ws->NO3_sed_i] -= growth * unitch * NO3_sed / DIN_sed / dz_sed / porosity;
    y1[ws->NH4_sed_i] -= NH4uptake * unitch / dz_sed / porosity;
    y1[ws->NO3_sed_i] -= NO3uptake * unitch / dz_sed / porosity;
    y1[ws->DIP_sed_i] -= growth * unitch * atk_W_P / dz_sed / porosity;
    y1[ws->DIC_sed_i] -= growth * unitch * atk_W_C / dz_sed / porosity;
    y1[ws->Oxygen_wc_i] += Oxy_pr;

    if (ws-> SG_N_gr_i > -1)
      	y1[ws->SG_N_gr_i] += growthrate;	
    if (ws-> SG_N_pr_i > -1)
      y1[ws->SG_N_pr_i] += growth * SEC_PER_DAY;
    if (ws->Oxy_pr_wc_i > -1)
      y1[ws->Oxy_pr_wc_i] += Oxy_pr * SEC_PER_DAY;
    
    /* Translocation of N between root and leaf */ 
    
    double translocate_to_roots = (SGfrac-SGROOT_N/(SG_N+SGROOT_N))*SGtransrate*(SG_N+SGROOT_N);

    y1[ws->SGROOT_N_i] += translocate_to_roots;
    y1[ws->SG_N_i] -= translocate_to_roots;
  }
}

void seagrass_spectral_grow_epi_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double SG_N = y[ws->SG_N_i] * unitch;
    double SGROOT_N = y[ws->SGROOT_N_i] * unitch;

    y[ws->EpiTN_i] += SG_N + SGROOT_N;
    y[ws->EpiTP_i] += (SG_N + SGROOT_N)* atk_W_P;
    y[ws->EpiTC_i] += (SG_N + SGROOT_N)* atk_W_C;
    y[ws->EpiBOD_i] += (SG_N + SGROOT_N)* atk_W_O;
}
