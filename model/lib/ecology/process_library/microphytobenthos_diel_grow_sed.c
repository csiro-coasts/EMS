/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/microphytobenthos_diel_grow_sed.c
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
 *  $Id: microphytobenthos_diel_grow_sed.c 6694 2021-03-24 01:08:59Z wil00y $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "microphytobenthos_diel_grow_sed.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20
#define UROUND 1.11e-16

typedef struct {
    int do_mb;                  /* flag */

    /*
     * parameters
     */
    double umax_t0;
    double aA;
    double psi;
    double Sh;
    double m;
    int n;
    double Plank_resp;
    int len;
    double* grow;
    double* rate;
    double NtoCHL;

    /*
     * tracers
     */
    int MPB_N_i;
    int MPB_NR_i;
    int MPB_I_i;
    int MPB_N_pr_i;
    int MPB_N_gr_i;
    int NO3_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Oxy_pr_i;
    int Light_i;
    int Chl_a_i;
    int TN_i;
    int TP_i;
    int TC_i;

    /*
     * common cell variables
     */
    int DNO3_i;
    int DPO4_i;
    int Tfactor_i;
    int umax_i;
} workspace;

void microphytobenthos_diel_grow_sed_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));

    p->workspace = ws;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, "MBumax");
    ws->aA = try_parameter_value(e, "MBaA");
    if (isnan(ws->aA)) {
        /*UR CHANGED to exit with some meaningful error message */
        double absorb = try_parameter_value(e, "MBabsorb");
        double rad = try_parameter_value(e, "MBrad");
        if(isnan(rad) || isnan(absorb))
          e->quitfn("If tracer 'MBaA' is not provided 'MBrad' and 'MBabsorb' are required to calculate microphytobenthos_diel_grow_sed!");
        ws->aA = aa(rad, absorb);
    }
    ws->psi = try_parameter_value(e, "MBpsi");
    if (isnan(ws->psi))
    {
      if(isnan(try_parameter_value(e, "MBrad")))
        e->quitfn("Either tracer 'MBpsi' or 'MBrad' are required to calculate microphytobenthos_diel_grow_sed!");
      ws->psi = psi(get_parameter_value(e, "MBrad"));
    }
    ws->Sh = get_parameter_value(e, "MBSh");
    ws->m = try_parameter_value(e, "MBm");
    if (isnan(ws->m))
    {
        /*UR CHANGED to exit with some meaningful error message */
        if(isnan(try_parameter_value(e, "MBrad")))
          e->quitfn("Either tracer 'MBm' or 'MBrad' are required to calculate microphytobenthos_diel_grow_sed!");
        ws->m = PhyCellMass(try_parameter_value(e, "MBrad"));

    }
    ws->n = rint(get_parameter_value(e, "MBn"));
    ws->Plank_resp = get_parameter_value(e, "Plank_resp");
    ws->len = extract_grow(ws->n, get_parameter_stringvalue(e, "MBtable"), &ws->grow, &ws->rate);
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");

    /*
     * tracers
     */
    ws->MPB_N_i = e->find_index(tracers, "MPB_N", e);
    ws->MPB_NR_i = e->find_index(tracers, "MPB_NR", e);
    ws->MPB_I_i = e->find_index(tracers, "MPB_I", e);
    ws->MPB_N_pr_i = e->find_index(tracers, "MPB_N_pr", e);
    ws->MPB_N_gr_i = e->find_index(tracers, "MPB_N_gr", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Oxy_pr_i = e->find_index(tracers, "Oxy_pr", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    /*
     * common cell variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, "DFumax", e);
}

void microphytobenthos_diel_grow_sed_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;

    ws->do_mb = (try_index(e->cv_model, "massbalance_sed", e) >= 0) ? 1 : 0;
}

void microphytobenthos_diel_grow_sed_destroy(eprocess* p)
{
    free(((workspace*) (p->workspace))->grow);
    free(((workspace*) (p->workspace))->rate);
    free(p->workspace);
}

void microphytobenthos_diel_grow_sed_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double MPB_N = y[ws->MPB_N_i];
    double MPB_NR = y[ws->MPB_NR_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    if (ws->do_mb) {
        y[ws->TN_i] += MPB_N + MPB_NR;
        y[ws->TP_i] += (MPB_N + MPB_NR) * red_W_P;
        y[ws->TC_i] += MPB_N * red_W_C;
    }
}

void microphytobenthos_diel_grow_sed_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double porosity = c->porosity;
    double tort2 = 1.0 - log(porosity * porosity);
	double MI_max = ws->m * red_A_I; /* mol Q . cell -1 */
	double MN_max = ws->m * red_A_N; /* mol N . cell -1 */
    double MPB_N = y[ws->MPB_N_i];
    double MPB_NR = y[ws->MPB_NR_i];
    double MPB_I = y[ws->MPB_I_i];
    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Light = y[ws->Light_i];
    double cells = MPB_N * mgN2molN / ws->m / red_A_N;
    double Itop = 2.77e18 * Light / AV;
    double DNO3 = cv[ws->DNO3_i] * porosity / tort2;
    double DPO4 = cv[ws->DPO4_i] * porosity / tort2;
    double diff[2];
    double conc[2];
	double Iuptake;
	double Nuptake;
    double umax = cv[ws->umax_i];
	double Iquota;
	double Nquota;

    diff[0] = DNO3;
    diff[1] = DPO4;
    conc[0] = DIN * mgN2molN;
    conc[1] = DIP * mgP2molP;

    {
	double tmp = cells * ws->aA * c->dz_sed;
	double I = (tmp > UROUND) ? Itop * (1.0 - exp(-tmp)) / tmp : Itop;

	/* LIGHT ABSORPTION */
		double KI = (ws->aA * I);   /*molQ cell-1 s-1*/
	/*	MPB_I = mol Q/m3
		MPB_I / (MPB_N /(ws->m * redA_N * MW_Nitr)) = mol Q/cell */

		if (MPB_N > 0.0)
			{
			if ((MPB_I * ws->m * red_A_N * MW_Nitr / MPB_N) > MI_max)/* only uptake if not already full */
			{Iuptake = 0.0;}
			else
			{Iuptake = KI*(MI_max - (MPB_I * ws->m * red_A_N * MW_Nitr / MPB_N))/MI_max;} /*molQ cell s-1*/

			Iuptake = Iuptake * MPB_N / (ws->m * red_A_N * MW_Nitr); /*molQ m-3 s-1*/

			Iquota = (MPB_I/ (MPB_N /(ws->m * red_A_N * MW_Nitr)))/MI_max;
			}
		else
			{Iuptake = 0.0;
			Iquota = 0.0;}

	/*	printf("I quota  ");
		printf("%f \n", Iquota);
		printf("Iuptake /d ");
		printf("%f \n", Iuptake*86400.);*/

	/* NUTRIENT UPTAKE,
	 * internal reserves primarily controlled by N uptake with P uptake at Redfield
	 * MPB_NR = mg N/m3
	 * MPB_NR / (MPB_N/(ws->m * red_A_N)) = mol N/cell
	 */

		double KN = (ws->psi * diff[0] * conc[0] * ws->Sh);   /*molN cell-1 s-1*/

		if (MPB_N > 0.0)
			{
			if ((MPB_NR * ws->m * red_A_N / MPB_N) > MN_max)
			{Nuptake = 0.0;}
			else
			{Nuptake = KN*(MN_max - (MPB_NR * ws->m * red_A_N / MPB_N))/MN_max;}  /* mol N. cell-1 s-1 */

			Nuptake = Nuptake * MPB_N / (ws->m * red_A_N);  /* mg N m-3 s-1 */

	 /* when insufficient P then internal reserves controlled by P uptake with N uptake at Redfield */
	    	if ((Nuptake * red_W_P * 7200.) > DIP)
			{Nuptake = DIP / (red_W_P * 7200.);}

	/* check N uptake does not exceed resources during timestep (dt =7200s) and reduce if necessary */
			if ((Nuptake * 7200.) > DIN)
			{Nuptake = DIN / 7200.;}

			Nquota = (MPB_NR / (MPB_N/(ws->m * red_A_N))) / MN_max;
			}
		else
			{Nuptake = 0.0;
			Nquota = 0.0;}
	/*	printf("N quota  ");
		printf("%f \n", Nquota);
		printf("Nuptake /d ");
		printf("%f \n", Nuptake*86400. );*/

	 /* PHYTOPLANKTON GROWTH s-1
	  * controlled by internal reserve of nutrient (at Redfield) cf maximum internal nutrient
	  * and internal reserve of light energy cf maximum internal light energy reserve
	  * limit quota to max of 1 to avoid excessive growth
	  */
		if (Iquota > 1.0) {Iquota = 1.0;}
		if (Nquota > 1.0) {Nquota = 1.0;}

		double growthrate = umax * Nquota * Iquota;
	 /*	printf("New growthrate /d   ");
		printf("%f \n", growthrate*86400.); */

     /* old CR growth: */
     /* double growthrate = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, I, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
     */

        double growth = MPB_N * growthrate;

	 /* PHYTOPLANKTON RESPIRATION molQ m-3 s-1
	 respiration = fraction due to growth (+ basal respiration : assume this included in mortality term)
	 */
        double Iresp = growth * ws->Plank_resp * red_A_I / (red_A_N * MW_Nitr);

	/*	UPDATE STATE VARIABLES */
		y1[ws->MPB_I_i] += Iuptake - (growth * red_A_I / (red_A_N * MW_Nitr)) - Iresp ;
 		y1[ws->MPB_NR_i] += Nuptake - growth;
        y1[ws->MPB_N_i] += growth;
        y1[ws->NH4_i] -= Nuptake * NH4 / DIN / porosity;
        y1[ws->NO3_i] -= Nuptake * NO3 / DIN / porosity;
        y1[ws->DIP_i] -= Nuptake * red_W_P / porosity;
        y1[ws->DIC_i] -= growth * red_W_C / porosity;

        y1[ws->Oxygen_i] += growth * red_W_O / porosity;
        y1[ws->MPB_N_pr_i] += growth * SEC_PER_DAY * red_W_C * c->dz_sed;
        y1[ws->MPB_N_gr_i] = growthrate * SEC_PER_DAY ;
        y1[ws->Oxy_pr_i] += growth * red_W_O / porosity * SEC_PER_DAY * c->dz_sed;

/*        y1[ws->MPB_N_i] += growth;
        y1[ws->MPB_N_pr_i] += growth * SEC_PER_DAY * c->dz_sed;
        y1[ws->MPB_N_gr_i] = growthrate / umax;
        y1[ws->NH4_i] -= growth * NH4 / DIN / porosity;
        y1[ws->NO3_i] -= growth * NO3 / DIN / porosity;
        y1[ws->DIP_i] -= growth * red_W_P / porosity;
        y1[ws->DIC_i] -= growth * red_W_C / porosity;
        y1[ws->Oxygen_i] += growth * red_W_O / porosity;
        y1[ws->Oxy_pr_i] += growth * red_W_O / porosity * SEC_PER_DAY * c->dz_sed; */
    }
}

void microphytobenthos_diel_grow_sed_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double MPB_N = y[ws->MPB_N_i];
    double MPB_NR = y[ws->MPB_NR_i];

    y[ws->Chl_a_i] += MPB_N / ws->NtoCHL;
    y[ws->TN_i] += MPB_N + MPB_NR;
    y[ws->TP_i] += (MPB_N + MPB_NR) * red_W_P;
    y[ws->TC_i] += MPB_N * red_W_C;
}
