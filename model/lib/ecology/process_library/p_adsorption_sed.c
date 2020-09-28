/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/p_adsorption_sed.c
 *  
 *  Description: Carries out P-adsorption-desorption for inorganic sediments in bed.
 *  Process implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: p_adsorption_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "p_adsorption_sed.h"

#define EFI_EPS 1.0e-4

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    
    double Pads_r_t0;
    double Pads_Ksed;
    double Pads_KO;
    double Pads_exp;
    double r_immob_PIP_t0;

    /*
     * tracers
     */
    int PIP_i;
    int PIPI_i;
    int DIP_i;
    int Oxygen_i;
    int TP_i;
    int EFI_i;

    /*
     * common variables
     */
    int Tfactor_i;
    int Pads_r_i;
    int r_immob_PIP_i;

} workspace;

void p_adsorption_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    
    ws->Pads_r_t0 = get_parameter_value(e, "Pads_r");
    ws->Pads_Ksed = get_parameter_value(e, "Pads_Ksed");
    ws->Pads_KO = get_parameter_value(e, "Pads_KO");
    ws->Pads_exp = get_parameter_value(e, "Pads_exp");
    ws->r_immob_PIP_t0 = get_parameter_value(e, "r_immob_PIP");

    /*
     * tracers
     */
    ws->EFI_i = e->find_index(e->tracers, "EFI", e);
    ws->PIP_i = e->find_index(tracers, "PIP", e);
    ws->PIPI_i = e->find_index(tracers, "PIPI", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->TP_i = e->find_index(tracers, "TP", e);

    /*
     * common variables
     */
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->Pads_r_i = find_index_or_add(e->cv_cell, "Pads_r", e);
    ws->r_immob_PIP_i = find_index_or_add(e->cv_cell, "r_immob_PIP", e);
}

void p_adsorption_sed_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
}

void p_adsorption_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void p_adsorption_sed_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;

    cv[ws->Pads_r_i] = ws->Pads_r_t0 * Tfactor;
    cv[ws->r_immob_PIP_i] = ws->r_immob_PIP_t0 * Tfactor;

    if (ws->do_mb) {
        double* y = c->y;
        double PIP = y[ws->PIP_i];
        double PIPI = y[ws->PIPI_i];

        y[ws->TP_i] += PIP + PIPI;
    }
}

void p_adsorption_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    double* y = ia->y;
    double* y1 = ia->y1;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double porosity = c->porosity;

    double Pads_r = cv[ws->Pads_r_i];
    double r_immob_PIP = cv[ws->r_immob_PIP_i];

    double Oxygen = y[ws->Oxygen_i];
    double DIP = y[ws->DIP_i];
    double EFI = y[ws->EFI_i];
    double PIP = y[ws->PIP_i];

    
    /*    double adsorption = Pads_r * EFI * ws->Pads_Ksed * DIP * (Oxygen / (ws->Pads_KO + e_max(Oxygen)));

	  double net_desorp = Pads_r * EFI * ws->Pads_Ksed * (pow(e_max(PIP) / (e_max(EFI + EFI_EPS) * ws->Pads_Ksed), npow)) - adsorption; */

    double net_desorp  = Pads_r * (PIP / (e_max(EFI) * ws->Pads_Ksed) -  DIP * (Oxygen / (ws->Pads_KO + e_max(Oxygen))));
       
    double PIP_immob = r_immob_PIP * PIP;
    
    y1[ws->DIP_i] += (net_desorp) / porosity;
    y1[ws->PIP_i] += -net_desorp - PIP_immob;
    y1[ws->PIPI_i] += PIP_immob;
}

void p_adsorption_sed_postcalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = ((cell*) pp);
    double* y = c->y;

    double PIP = y[ws->PIP_i];
    double PIPI = y[ws->PIPI_i];
    y[ws->TP_i] += PIP + PIPI;
}
