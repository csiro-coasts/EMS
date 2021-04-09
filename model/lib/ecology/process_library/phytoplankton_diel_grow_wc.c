/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/phytoplankton_diel_grow_wc.c
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
 *  $Id: phytoplankton_diel_grow_wc.c 6698 2021-03-24 01:11:43Z wil00y $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "stringtable.h"
#include "utils.h"
#include "ecofunct.h"
#include "constants.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "phytoplankton_diel_grow_wc.h"

#define EPS_DIN 1.0e-20
#define EPS_DIP 1.0e-20

typedef struct {
    int do_mb;                  /* flag */
    /*
     * flag: 1 for large phytoplankton, 0 for small phytoplankton
     */
    int large;

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
    int Phy_N_i;
    int Phy_I_i;
    int Phy_NR_i;
    int Phy_N_pr_i;
    int Phy_N_gr_i;
    int NO3_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;
    int Oxygen_i;
    int Oxy_pr_i;
    int Light_i;
    int Chl_a_i;
    int Kd_i;
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

void phytoplankton_diel_grow_wc_init(eprocess* p)
{
    ecology* e = p->ecology;
    stringtable* tracers = e->tracers;
    workspace* ws = malloc(sizeof(workspace));
    char* prm = p->prms->se[0]->s;
    int large;

    p->workspace = ws;

    if (toupper(prm[0]) == 'L')
        ws->large = 1;
    else if (toupper(prm[0]) == 'S')
        ws->large = 0;
    else
        e->quitfn("ecology: error: \"%s\": \"%s(%s)\": unexpected parameter: expected \"small\" or \"large\"\n", e->processfname, p->name, prm);
    large = ws->large;

    /*
     * parameters
     */
    ws->umax_t0 = get_parameter_value(e, (large) ? "PLumax" : "PSumax");
    ws->aA = try_parameter_value(e, (large) ? "PLaA" : "PSaA");
    if (isnan(ws->aA)) {
        double absorb = get_parameter_value(e, (large) ? "PLabsorb" : "PSabsorb");
        double rad = get_parameter_value(e, (large) ? "PLrad" : "PSrad");

        ws->aA = aa(rad, absorb);
    }
    ws->psi = try_parameter_value(e, (large) ? "PLpsi" : "PSpsi");
    if (isnan(ws->psi))
        ws->psi = psi(get_parameter_value(e, (large) ? "PLrad" : "PSrad"));
    ws->Sh = get_parameter_value(e, (large) ? "PLSh" : "PSSh");
    ws->m = try_parameter_value(e, (large) ? "PLm" : "PSs");
    if (isnan(ws->m))
        ws->m = PhyCellMass(get_parameter_value(e, (large) ? "PLrad" : "PSrad"));
    ws->n = rint(get_parameter_value(e, (large) ? "PLn" : "PSn"));
    ws->Plank_resp = get_parameter_value(e, "Plank_resp");
    ws->len = extract_grow(ws->n, get_parameter_stringvalue(e, (large) ? "PLtable" : "PStable"), &ws->grow, &ws->rate);
    ws->NtoCHL = get_parameter_value(e, "NtoCHL");

    /*
     * tracers
     */
    ws->Phy_N_i = e->find_index(tracers, (large) ? "PhyL_N" : "PhyS_N", e);
    ws->Phy_I_i = e->find_index(tracers, (large) ? "PhyL_I" : "PhyS_I", e);
    ws->Phy_NR_i = e->find_index(tracers, (large) ? "PhyL_NR" : "PhyS_NR", e);
    ws->NO3_i = e->find_index(tracers, "NO3", e);
    ws->NH4_i = e->find_index(tracers, "NH4", e);
    ws->DIP_i = e->find_index(tracers, "DIP", e);
    ws->DIC_i = e->find_index(tracers, "DIC", e);
    ws->Oxygen_i = e->find_index(tracers, "Oxygen", e);
    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Phy_N_pr_i = e->find_index(tracers, (large) ? "PhyL_N_pr" : "PhyS_N_pr", e);
    ws->Phy_N_gr_i = e->find_index(tracers, (large) ? "PhyL_N_gr" : "PhyS_N_gr", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->Kd_i = e->find_index(tracers, "Kd", e);
    ws->Oxy_pr_i = e->find_index(tracers, "Oxy_pr", e);
    ws->TN_i = e->find_index(tracers, "TN", e);
    ws->TP_i = e->find_index(tracers, "TP", e);
    ws->TC_i = e->find_index(tracers, "TC", e);

    /*
     * common cell variables
     */
    ws->DNO3_i = find_index(e->cv_cell, "DNO3", e);
    ws->DPO4_i = find_index(e->cv_cell, "DPO4", e);
    ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
    ws->umax_i = find_index_or_add(e->cv_cell, (large) ? "PLumax" : "PSumax", e);

}

void phytoplankton_diel_grow_wc_postinit(eprocess* p)
{
    ecology* e = p->ecology;
    workspace* ws = p->workspace;
    double v1,v2,v5;
    
    ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;
    
   /*
    * Not valid during a pre_build (RECOM)
    */
   if (!e->pre_build) {
    /* test for equal sinking rates of structural material and reserves */
    if (ws->large) {
      v1 = einterface_gettracersvel(e->model,"PhyL_N");
      v2 = einterface_gettracersvel(e->model,"PhyL_NR");
      v5 = einterface_gettracersvel(e->model,"PhyL_I");
      
      if ((v1!=v2)||(v1!=v5)){
	printf("Mass conservation violation due to microalgae reserves \n");
	printf("and structural material sinking at different rates \n");
	printf("Editing of .prm file required. \n");
	printf("Sinking of PhyL :N %e, NR %e, I %e are not equal  \n",v1,v2,v5);
	exit(-1);
      }
    } else {
      v1 = einterface_gettracersvel(e->model,"PhyS_N");
      v2 = einterface_gettracersvel(e->model,"PhyS_NR");
      v5 = einterface_gettracersvel(e->model,"PhyS_I");
      
      if ((v1!=v2)||(v1!=v5)){
	printf("Mass conservation violation due to microalgae reserves \n");
	printf("and structural material sinking at different rates \n");
	printf("Editing of .prm file required. \n");
	printf("Sinking of PhyS :N %e, NR %e, I %e are not equal \n ",v1,v2,v5);
	exit(-1);
      }
    }
   }
}

void phytoplankton_diel_grow_wc_destroy(eprocess* p)
{
    workspace* ws = (workspace*) p->workspace;

    free(ws->grow);
    free(ws->rate);
    free(ws);
}

void phytoplankton_diel_grow_wc_precalc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    cell* c = (cell*) pp;
    double* cv = c->cv;
    double* y = c->y;

    double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
    double Phy_N = y[ws->Phy_N_i];
    double Phy_NR = y[ws->Phy_NR_i];

    cv[ws->umax_i] = ws->umax_t0 * Tfactor;

    y[ws->Kd_i] += ws->aA * Phy_N * mgN2molN / ws->m / red_A_N;

    if (ws->do_mb) {
    	y[ws->TN_i] += Phy_N + Phy_NR;
    	y[ws->TP_i] += Phy_N * red_W_P + Phy_NR * red_W_P ;
    	y[ws->TC_i] += Phy_N * red_W_C ;
/*      y[ws->TN_i] += Phy_N;
        y[ws->TP_i] += Phy_N * red_W_P;
        y[ws->TC_i] += Phy_N * red_W_C;*/
    }
}

void phytoplankton_diel_grow_wc_calc(eprocess* p, void* pp)
{
    workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);
    double* cv = c->cv;
    double* y = ia->y;
    double* y1 = ia->y1;

    double PI_max = ws->m * red_A_I; /* mol Q . cell -1 */
    double PN_max = ws->m * red_A_N; /* mol N . cell -1 */
    double Phy_N = y[ws->Phy_N_i];
    double Phy_I = y[ws->Phy_I_i];
    double Phy_NR = y[ws->Phy_NR_i];
    double NO3 = y[ws->NO3_i];
    double NH4 = y[ws->NH4_i];
    double din = NO3 + NH4;
    double DIN = (din > 0.0) ? din : EPS_DIN;
    double dip = y[ws->DIP_i];
    double DIP = (dip > 0.0) ? dip : EPS_DIP;
    double Light = y[ws->Light_i];
    double I = 2.77e18 * Light / AV;
    double diff[2];
    double conc[2];
    double Iuptake;
    double Nuptake;
    double umax = cv[ws->umax_i];
    double Iquota;
    double Nquota;
    double e_dt = c->col->e->dt;

    diff[0] = cv[ws->DNO3_i];
    diff[1] = cv[ws->DPO4_i];
    conc[0] = DIN * mgN2molN;
    conc[1] = DIP * mgP2molP;

    {

	/* LIGHT ABSORPTION */
		double KI = (ws->aA * I);   /*molQ cell-1 s-1*/
	/*	Phy_I = mol Q/m3
		Phy_I / (Phy_N /(ws->m * redA_N * MW_Nitr)) = mol Q/cell */

		if (Phy_N > 0.0)
			{
			if ((Phy_I * ws->m * red_A_N * MW_Nitr / Phy_N) > PI_max)/* only uptake if not already full */
			{Iuptake = 0.0;}
			else
			{Iuptake = KI*(PI_max - (Phy_I * ws->m * red_A_N * MW_Nitr / Phy_N))/PI_max;} /*molQ cell s-1*/

			Iuptake = Iuptake * Phy_N / (ws->m * red_A_N * MW_Nitr); /*molQ m-3 s-1*/

			Iquota = (Phy_I/ (Phy_N /(ws->m * red_A_N * MW_Nitr)))/PI_max;
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
	 * Phy_NR = mg N/m3
	 * Phy_NR / (Phy_N/(ws->m * red_A_N)) = mol N/cell
	 */

		double KN = (ws->psi * diff[0] * conc[0] * ws->Sh);   /*molN cell-1 s-1*/

		if (Phy_N > 0.0)
			{
			if ((Phy_NR * ws->m * red_A_N / Phy_N) > PN_max)
			{Nuptake = 0.0;}
			else
			{Nuptake = KN*(PN_max - (Phy_NR * ws->m * red_A_N / Phy_N))/PN_max;}  /* mol N. cell-1 s-1 */

			Nuptake = Nuptake * Phy_N / (ws->m * red_A_N);  /* mg N m-3 s-1 */

	 /* when insufficient P then internal reserves controlled by P uptake with N uptake at Redfield */
	    	if ((Nuptake * red_W_P * e_dt) > DIP)
			{Nuptake = DIP / (red_W_P * e_dt);}

	/* check N uptake does not exceed resources during timestep (dt = ecology timestep) and reduce if necessary */
			if ((Nuptake * e_dt) > DIN)
			{Nuptake = DIN / e_dt;}

			Nquota = (Phy_NR / (Phy_N/(ws->m * red_A_N))) / PN_max;
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
     /*   growthrate = plankgrow(umax, ws->aA, ws->psi, ws->Sh, diff, ws->m, I, conc, ws->len, ws->n, ws->grow, ws->rate, ws->Plank_resp);
		printf("growthrate /d   ");
		printf("%f \n", growthrate*86400.);*/

        double growth = Phy_N * growthrate;
	/* growth is transfer of nutrient reserve NR to structual cell nutrient N
	 * with corresponding reduction in energy reserve
	 * C uptake is sufficient to keep Phy_N at Redfield */

	 /* PHYTOPLANKTON RESPIRATION molQ m-3 s-1
	 respiration = fraction due to growth (+ basal respiration : assume this included in mortality term)
	 */
        double Iresp = growth * ws->Plank_resp * red_A_I / (red_A_N * MW_Nitr);

	/*	UPDATE STATE VARIABLES */
		y1[ws->Phy_I_i] += Iuptake - (growth * red_A_I / (red_A_N * MW_Nitr)) - Iresp ;
 		y1[ws->Phy_NR_i] += Nuptake - growth;
        y1[ws->Phy_N_i] += growth;
        y1[ws->NH4_i] -= Nuptake * NH4 / DIN;
        y1[ws->NO3_i] -= Nuptake * NO3 / DIN;
        y1[ws->DIP_i] -= Nuptake * red_W_P;
        y1[ws->DIC_i] -= growth * red_W_C;

 /*       y1[ws->NH4_i] -= growth * NH4 / DIN;
        y1[ws->NO3_i] -= growth * NO3 / DIN;
        y1[ws->DIP_i] -= growth * red_W_P;
        y1[ws->DIC_i] -= growth * red_W_C;*/

        y1[ws->Oxygen_i] += growth * red_W_O;

        y1[ws->Phy_N_pr_i] += growth * SEC_PER_DAY * red_W_C;
        y1[ws->Phy_N_gr_i] = growthrate * SEC_PER_DAY;
        y1[ws->Oxy_pr_i] += growth * red_W_O * SEC_PER_DAY;
    }
}

void phytoplankton_diel_grow_wc_postcalc(eprocess* p, void* pp)
{
    cell* c = ((cell*) pp);
    workspace* ws = p->workspace;
    double* y = c->y;

    double Phy_N = y[ws->Phy_N_i];
    double Phy_NR = y[ws->Phy_NR_i];

    y[ws->Chl_a_i] += Phy_N / ws->NtoCHL;
    y[ws->TN_i] += Phy_N + Phy_NR;
    y[ws->TP_i] += Phy_N * red_W_P + Phy_NR * red_W_P ;
    y[ws->TC_i] += Phy_N * red_W_C ;
}
