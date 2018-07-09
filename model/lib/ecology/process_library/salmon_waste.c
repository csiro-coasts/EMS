/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/salmon_waste.c
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
 *  $Id: salmon_waste.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "salmon_waste.h"
#include "einterface.h"

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double FCR_t0;
    double C_in_feed;
    double N_in_feed;
    double P_in_feed;
    double feed_FPOC;
    double feed_FPON;
    double feed_FDIN;
    double feed_FPOP;
    double feed_FDIP;
  /*double fish_resp_t0;*/

    /*
     * tracers
     */
    int DetPL_N_i;
    int DetR_C_i;
    int DetR_P_i;
    int NH4_i;
    int DIP_i;
  /*int Oxygen_i;
    int Oxy_pr_i;
    int NH4_pr_i;
    int TN_i;
    int TC_i;
    int TP_i;
    int BOD_i;
    int COD_i;*/
    int Fish_size_i;
    int Fish_num_i;
    int Fish_N_pr_i;
    int Fish_P_pr_i;

    /*
     * common variables
     */
    int Tfactor_i;
  /*int fish_resp_i;*/
    int FCR_i;
} workspace;

void salmon_waste_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->FCR_t0 = get_parameter_value(e, "FCR");
    /*ws->fish_resp_t0 = get_parameter_value(e, "fish_resp");*/
    ws->C_in_feed = get_parameter_value(e, "C_in_feed");
    ws->N_in_feed = get_parameter_value(e, "N_in_feed");
    ws->P_in_feed = get_parameter_value(e, "P_in_feed");
    ws->feed_FPOC = get_parameter_value(e, "feed_FPOC");
    ws->feed_FPON = get_parameter_value(e, "feed_FPON");
    ws->feed_FDIN = get_parameter_value(e, "feed_FDIN");
    ws->feed_FPOP = get_parameter_value(e, "feed_FPOP");
    ws->feed_FDIP = get_parameter_value(e, "feed_FDIP");
 
    /*
     * tracers
     */
    ws->DetPL_N_i = e->find_index(tracers, "DetPL_N", e);
    ws->DetR_C_i = e->find_index(tracers, "DetR_C", e);
    ws->DetR_P_i = e->find_index(tracers, "DetR_P", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
 /* ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Oxy_pr_i = e->try_index(tracers, "Oxy_pr", e);
    ws->NH4_pr_i = e->try_index(tracers, "NH4_pr", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TC_i = e->find_index(tracers, "TC", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->BOD_i = e->find_index(tracers, "BOD", e);
    ws->COD_i = e->try_index(tracers, "COD", e);*/
    ws->Fish_size_i = e->try_index(tracers, "Fish_size", e);
    ws->Fish_num_i = e->try_index(tracers, "Fish_num", e);
   /*ws->botz_i = e->find_index(tracers, "botz", e);*/

    /*non essential diagnostic tracer*/
    
    ws->Fish_N_pr_i = e->try_index(tracers, "Fish_N_pr", e);
    ws->Fish_P_pr_i = e->try_index(tracers, "Fish_P_pr", e);
 
    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->FCR_i = find_index_or_add(e->cv_cell, "FCR", e);
    /*ws->fish_resp_i = find_index_or_add(e->cv_cell, "fish_resp", e); */

    eco_write_setup(e,"Putting in salmon waste \n");

}

/* KWA don't think we need this so long as salmon_waste process called before remineralization

void salmon_waste_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, (p->type == PT_WC) ? "massbalance_wc" : "massbalance_sed", e) >= 0) ? 1 : 0;
}
*/
void salmon_waste_destroy(eprocess* p)
{
    free(p->workspace);
}

void salmon_waste_precalc(eprocess* p, void* pp)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    /*double layer_thickness = c->dz_wc;
      double botz = y[ws->botz_i]; */

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->FCR_i] = Tfactor * ws->FCR_t0;
    /*    cv[ws->fish_resp_i] = Tfactor * ws->fish_resp_t0;*/

    // Bail out if deeper than 20m
    double z_centre = einterface_getcellz(c->col->model,c->b,c->k_wc);
    if (z_centre < -20.0)
      return;

    double Fish_size = y[ws->Fish_size_i];
    double Fish_num = y[ws->Fish_num_i];

    // Bail out if Fish = NaN
    if (isnan(Fish_size) || isnan(Fish_num))
      return;

    // Bail out if Fish <= 0
    if ((Fish_size <= 0) || (Fish_num <=0))
      return;

    /* QC cap fish size at 8000g*/ 
    double Fish_size_qc = min(Fish_size, 8000.);

    /*    double DetPL_N = y[ws->DetPL_N_i];
    double DetR_C = y[ws->DetR_C_i];
    double DetR_P = y[ws->DetR_P_i];*/

    int wcbotk; 
    wcbotk = einterface_getwcbotk(c->col->model, c->b);
    double z_bot = einterface_getcellz(c->col->model,c->b,wcbotk);
    double pen_volume = einterface_cellarea(c->col->e->model, c->b) * -(max(-20,z_bot));
    /*  double pen_volume = einterface_cellarea(c->col->model, c->b) * -(max(-20,z_bot));*/
 // Bail out if dry cell
    if (pen_volume <= 0)
      return;

    /* assuming t:units = "days since 2010-01-01 00:00:00 +0" find week in year 0<week<=52 */
    double time = einterface_getmodeltime(e->model) / SEC_PER_DAY;
    double week = (time - (floor((time) / 365.25) * 365.25)) / 7.0;

    /* salmon growth model coefficients */
    double ASa0 = 0.119;
    double ASa1 = -0.149;
    double ASb0 = 0.705;
    double ASb1 = 0.713;
    double ASd = 0.092;
    double ASo = 7.54;

    /*Salmon growth in g/fish for a week by Shane Richards (his model requires kg) */
    double AS_growth = Fish_size_qc + ((ASa0/ASa1)*(1+ASd*cos((2.0*M_PI/52.0)*(week-ASo)))*(1.0-exp(-ASa1*(Fish_size_qc/1000.0)))*exp(-ASb0*pow((Fish_size_qc/1000.0),ASb1)))*1000.0;

    /* feed required to sustain growth for all fish in pen kg dry weight /week */
    double feed = (AS_growth - Fish_size_qc) * Fish_num * cv[ws->FCR_i] / 1000.0;

    /* Salmon waste kg/week for whole pen */
    /* ?? need to divide by volume of pen max(botz,-20); how do point source loads in mg/s work?? */

    double SW_POC = feed * ws->C_in_feed * ws->feed_FPOC;
    double SW_PON = feed * ws->N_in_feed * ws->feed_FPON;
    double SW_DIN = feed * ws->N_in_feed * ws->feed_FDIN;
    double SW_POP = feed * ws->P_in_feed * ws->feed_FPOP;
    double SW_DIP = feed * ws->P_in_feed * ws->feed_FDIP;

    /* unit conversion to mg/m3/dt (calculation applied outside of integrator) */
    y[ws->DetPL_N_i] += (SW_PON * ((1000.0*1000.0)/(7.0*SEC_PER_DAY)) * e->dt)/pen_volume;
    y[ws->DetR_C_i] += ((SW_POC - (SW_PON * red_W_C)) * ((1000.0*1000.0)/(7.0*SEC_PER_DAY)) * e->dt)/pen_volume;
    y[ws->DetR_P_i] += ((SW_POP - (SW_PON * red_W_P)) * ((1000.0*1000.0)/(7.0*SEC_PER_DAY)) * e->dt)/pen_volume;
    y[ws->NH4_i] += (SW_DIN * ((1000.0*1000.0)/(7.0*SEC_PER_DAY)) * e->dt)/pen_volume;
    y[ws->DIP_i] += (SW_DIP * ((1000.0*1000.0)/(7.0*SEC_PER_DAY)) * e->dt)/pen_volume;

    /* Salmon waste TN & TP in mg/m3/d */
    /*if (ws->Fish_N_pr_i > -1)*/
    y[ws->Fish_N_pr_i] = (SW_PON + SW_DIN)*(1000.*1000.)/(7.* pen_volume);
    /*if (ws->Fish_P_pr_i > -1)*/
    y[ws->Fish_P_pr_i] = (SW_POP - (SW_PON * red_W_P) + SW_DIP)*(1000.*1000.)/(7.* pen_volume);

}

void salmon_waste_postcalc(eprocess* p, void* pp)
{
   /* KWA don't think we need this so long as salmon_waste process called before remineralization
     workspace* ws = p->workspace;
    cell* c = ((cell*) pp);
    double* y = c->y;

    double porosity = (c->type == CT_SED) ? c->porosity : 1.0;

    double NH4 = y[ws->NH4_i];
    double NO3 = y[ws->NO3_i];
    double DetPL_N = y[ws->DetPL_N_i];
    double DetBL_N = y[ws->DetBL_N_i];
    double DetR_N = y[ws->DetR_N_i];
    double DIC = y[ws->DIC_i];
    double DIP = y[ws->DIP_i];
    double DetR_C = y[ws->DetR_C_i];
    double DetR_P = y[ws->DetR_P_i];
    double DOR_N = y[ws->DOR_N_i];
    double DOR_C = y[ws->DOR_C_i];
    double DOR_P = y[ws->DOR_P_i];

    y[ws->TN_i] += (NH4 + NO3 + DOR_N) * porosity + DetPL_N + DetBL_N + DetR_N;
    y[ws->TC_i] += (DIC + DOR_C) * porosity + DetPL_N * red_W_C + DetBL_N * atk_W_C + DetR_C;
    y[ws->TP_i] += (DIP + DOR_P) * porosity + DetPL_N * red_W_P + DetBL_N * atk_W_P + DetR_P;
    y[ws->BOD_i] += (DOR_C * C_O_W * porosity + DetPL_N * red_W_O + DetBL_N * atk_W_O + DetR_C * C_O_W);

   */
}
