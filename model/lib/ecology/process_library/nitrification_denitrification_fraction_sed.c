/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nitrification_denitrification_fraction_sed.c
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
 *  $Id: nitrification_denitrification_fraction_sed.c 5846 2018-06-29 04:14:26Z riz008 $
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
#include "nitrification_denitrification_fraction_sed.h"

typedef struct {
    /*
     * parameters
     */
    double Fmax_Nit_sed;
    double KO_Nit;
    double KO_Den;

    /*
     * tracers
     */
    int NH4_i;
    int NH4_pr_i;
    int NO3_i;
    int Oxygen_i;
    int Oxy_pr_i;
    int Den_fl_i;

    /*
     * common variables
     */
    int NH4_remin_pr_i;
} workspace;

void nitrification_denitrification_fraction_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->Fmax_Nit_sed = get_parameter_value(e, "Fmax_Nit_sed");
    ws->KO_Nit = get_parameter_value(e, "KO_Nit");
    ws->KO_Den = get_parameter_value(e, "KO_Den");

    /*
     * tracers
     */
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->NH4_pr_i = e->find_index(tracers, "NH4_pr", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Oxy_pr_i = e->find_index(tracers, "Oxy_pr", e);
    ws->Den_fl_i = e->find_index(tracers, "Den_fl", e);

    /*
     * common variables
     */
    /*
     * (NH4_remin_pr should already be created by remineralization.)
     */
    ws->NH4_remin_pr_i = find_index(e->cv_cell, "NH4_remin_pr", e);
}


void nitrification_denitrification_fraction_sed_postinit(eprocess* p)
{
    ecology* e = p->ecology;

    if( process_present(e,PT_SED,"nitrification_denitrification_sed"))
    {
      emstag(LPANIC,"eco:nitrification_denitrification_fraction_sed:postinit","nitrification_denitrification_sed seems to be enabled, these processes are mutually exclusive!");
    }
}

void nitrification_denitrification_fraction_sed_destroy(eprocess* p)
{
    free(p->workspace);
}

void nitrification_denitrification_fraction_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = (cell*) ia->media;
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double porosity = c->porosity;

    double Oxygen = y[ws->Oxygen_i];

    double F_Nit = ws->Fmax_Nit_sed * Oxygen * Oxygen / (ws->KO_Nit * ws->KO_Nit + e_max(Oxygen) * e_max(Oxygen));
    double F_Den = ws->KO_Den / (ws->KO_Den + e_max(Oxygen));

    double Nitrification = F_Nit * cv[ws->NH4_remin_pr_i] / porosity;
    double Denitrification = F_Den * Nitrification;

    y1[ws->NH4_i] -= Nitrification;
    y1[ws->NH4_pr_i] -= Nitrification * SEC_PER_DAY * c->dz_sed * porosity;
    y1[ws->NO3_i] += Nitrification - Denitrification;
    y1[ws->Den_fl_i] += Denitrification * SEC_PER_DAY * c->dz_sed * porosity;

    /*
     * KWA added : 2 moles of DO lost for every mole of N == 4.57 g (NIT_N_0)
     *             DO per g of N nitrified and not denitrified 
     */
    y1[ws->Oxygen_i] -= (Nitrification - Denitrification) * NIT_N_0 ;
    y1[ws->Oxy_pr_i] -= (Nitrification - Denitrification) * NIT_N_0 
                                            * porosity * c->dz_sed * SEC_PER_DAY;
}
