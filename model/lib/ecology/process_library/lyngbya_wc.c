/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/lyngbya_growth_wc.c
 *  
 *  Description:
 *  Process implementation template
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: lyngbya_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

/**
 * Process Template
 * -  Make sure that all functions used are well documented, functionaly and
 *    scientificly!
 * -  Be aware of processes your process might depend on, even though there are
 *    means to check by name if a   process is active - that might not be a good
 *    longterm solution, names may change.
 * -  Keep I/O (reading files) to a minimum or out of here if possible! If you must,
 *    make sure it is nly reaDOR_n_id once
 *
 *
 *
 */
 /* Any key word alled <something>MARKER is used by the wizard and should be
  * removed if you create a process manually.
  */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ecology_internal.h"
#include "utils.h"
#include "ecofunct.h"                                  
#include "constants.h"
#include "eprocess.h"
#include "stringtable.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "lyngbya_wc.h"

#define EPS 1.0e-20

typedef struct {
	int do_mb;                  /* flag */


/* Parameters */
    double NtoCHL;
    double Phy_Kdin;
    double Phy_KP;
    double Phy_pmax;
    double Phy_alpha;
    double Plank_resp;
    double Phy_gamma;
    double Plank_mort;
    
/* Tracers */
    int PhyS_N_i;
    int NO3_i;
    int Oxygen_i;
    int Chl_a_i;
    int Light_i;
    int TN_i;
    int TP_i;
    int TC_i;
    int NH4_i;
    int DIP_i;
    int DIC_i;

    int Kd_i;
    int temp_i;
    int salt_i; 
/*diagnostic variable */
    int Phy_N_pr_i;
    int Phy_N_gr_i;

/* Common Cell Variables */
    int Tfactor_i;
    int Phy_pmax_i;

} workspace;


/* This is only called once during the lifetime of a process/ecology. */

void lyngbya_wc_init(eprocess* p)
{
	ecology* e = p->ecology;

	stringtable* tracers = e->tracers;
	workspace* ws = malloc(sizeof(workspace));

	p->workspace = ws;

/* PARAMETERS */
  ws->NtoCHL = get_parameter_value(e,"NtoCHL_L_N2");
  ws->Phy_Kdin = get_parameter_value(e,"Phy_L_N2_Kdin");
  ws->Phy_KP = get_parameter_value(e,"Phy_L_N2_KP");
  ws->Phy_pmax = get_parameter_value(e,"Phy_L_N2_pmax");
  ws->Phy_alpha = get_parameter_value(e,"Phy_L_N2_alpha");
  ws->Plank_resp = get_parameter_value(e, "Plank_L_N2_resp");
  ws->Plank_mort = get_parameter_value(e, "Plank_L_N2_mort");


/* TRACERS */
  ws->PhyS_N_i = e->find_index(tracers,"Phy_L_N2", e);
  ws->NO3_i = e->find_index(tracers,"NO3", e);
  ws->DIC_i = e->find_index(tracers,"DIC", e);
  ws->Oxygen_i = e->find_index(tracers,"Oxygen", e);
  ws->NH4_i = e->find_index(tracers,"NH4", e);
  ws->DIP_i = e->find_index(tracers,"DIP", e);
  ws->Light_i = e->find_index(tracers,"Light", e);

  ws->TN_i = e->find_index(tracers, "TN", e);
  ws->TP_i = e->find_index(tracers, "TP", e);
  ws->TC_i = e->find_index(tracers, "TC", e);

  ws->temp_i = e->find_index(tracers,"temp", e);
  ws->salt_i = e->find_index(tracers,"salt", e);
  ws->Chl_a_i = e->find_index(tracers,"Chl_a", e);


/* diagnostic fluxes TRACERS */ 

  ws->Phy_N_pr_i = e->find_index(tracers,"Phy_L_N2_fix", e);



/* COMMON CELL VARIABLES */
  ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
  ws->Phy_pmax_i = find_index_or_add(e->cv_cell, "Phy_pmax", e);
}
void lyngbya_wc_postinit(eprocess* p)
{
	ecology* e = p->ecology;
	workspace* ws = p->workspace;
	ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

}

void lyngbya_wc_destroy(eprocess* p)
{
	free(p->workspace);

}


void lyngbya_wc_precalc(eprocess* p, void* pp)
{

  /*mass balance and updade the light attenuation coefficent*/

	workspace* ws = p->workspace;
	cell* c = (cell*) pp;
	double* y = c->y;
	double* cv = c->cv;
        

	double Tfactor = (ws->Tfactor_i >= 0) ? cv[ws->Tfactor_i] : 1.0;
	cv[ws->Phy_pmax_i] = ws->Phy_pmax ;

	double PhyS_N = y[ws->PhyS_N_i];


	/*mass balance */
	if (ws->do_mb) {

		y[ws->TN_i] += PhyS_N ;
                y[ws->TP_i] += PhyS_N * atk_W_P;
  		y[ws->TC_i] += PhyS_N * atk_W_C;


	}

}




/**
 * Documentation:
 * Add context and scientific documentation here!
 *
 */
void lyngbya_wc_calc(eprocess* p, void* pp)
{
	workspace* ws = p->workspace;
	intargs* ia = (intargs*) pp;
	cell* c = ((cell*) ia->media);
	double* cv = c->cv;
	double* y = ia->y;
	double* y1 = ia->y1;

	/*Assign initial tracer values which is integrated i.e.*/
    double NH4 = y[ws->NH4_i]; 
    double NO3 = y[ws->NO3_i];
    double DIP = y[ws->DIP_i];
    double Phy_N = y[ws->PhyS_N_i];
    double Light = y[ws->Light_i];
    double Pmax = cv[ws->Phy_pmax_i] ;
    double temp = y[ws->temp_i]; 
    double salt = y[ws->salt_i]; 

    
/*local declaration */
   
double TSfact;
double limitE;
double utot;
double limitDIP;
double DIN=NH4+NO3;
double  limitDIN;  
double limitPE;
double limitnfix;
double limitEfix;
double u_a;
double uNfix;
double growth;
double resp;
double din_uptake;
double nh4_uptake;
double no3_uptake;
double dip_uptake;
double dic_uptake;






	/*temperature and salinity effect on Lym growth */

    TSfact=   -0.7+(0.005*pow(temp,2))-0.000157*pow(temp,3)+0.001206*pow(salt,2)- 0.00001*pow(salt,3)-0.00008053*temp*pow(salt,2)+(0.00009659*pow(temp,2))*salt;
    if (temp<=7.0) 
{
        TSfact=0.0;}
    
    if (salt<=9.0)
{
        TSfact=0.0;}
    
/**/

utot=TSfact*ws->Phy_pmax_i ;
/* Light*/

limitE= 1-exp(-ws->Phy_alpha*Light/utot); /*% former equation P= +uref*alpha*E/(sqrt(uref^2+alpha*E^2))*/

/*Nutrients*/

limitDIP =  DIP/(ws->Phy_KP+DIP); /* P limitation factor*/

limitDIN =  DIN/(ws->Phy_Kdin+DIN); /* DIN limitation factor*/

limitPE  =  e_min(limitE,limitDIP);/**/

limitPE=max(limitPE,0);

/*following law of the miniumu the combined light and Phos limitation factor  */
limitnfix  =  (limitPE-limitDIN);/* fractiona N fixation requried*/

limitnfix=max(limitnfix,0);

limitEfix =  limitE - Phy_N*limitnfix; /* n fixation should ahve an energy cost that will lower P where gamma is the fraction cost of fixing N*/

u_a   =  utot*e_min(limitEfix,limitDIP); /*%actual specific growth rate*/

uNfix =  (u_a - (utot*limitDIN));  /*% actual rate of N fixtation*/



/* growth fluxes and respiration fluxes  */
growth=u_a *Phy_N;
resp= growth* ws->Plank_resp;


/* nutrient uptake*/
din_uptake=-(u_a-uNfix)*Phy_N;
nh4_uptake=-growth*NH4/DIN;
no3_uptake= -growth*NO3/DIN;

dip_uptake= din_uptake * atk_W_P;
dic_uptake= din_uptake * atk_W_C;


    /* updating the change in the state variables*/
	
    y1[ws->PhyS_N_i] += growth;
    y1[ws->PhyS_N_i] -= resp;
    y1[ws->PhyS_N_i] -=  ws->Plank_mort*Phy_N*Phy_N;

 
    y1[ws->NH4_i]  -= nh4_uptake;
    y1[ws->NO3_i]  -= no3_uptake;
    y1[ws->DIP_i]  -= dip_uptake; 
    y1[ws->DIC_i]  -= (dic_uptake);
    y1[ws->DIC_i]  += resp * atk_W_C;
    y1[ws->Oxygen_i] += (growth - resp) * atk_W_O;



    /* upade the fluxes variables */
    y1[ws->Phy_N_pr_i] +=uNfix * Phy_N* SEC_PER_DAY ;
   
 }


/**
 * Documentation:
 * Add context and scientific documentation here!
 *
 */
void lyngbya_wc_postcalc(eprocess* p, void* pp)
{
	cell* c = ((cell*) pp);
	workspace* ws = p->workspace;
	double* y = c->y;

  	double Phy_N = y[ws->PhyS_N_i];
	/*calculate the chlorophyll*/
	y[ws->Chl_a_i]  =+(((Phy_N/ws->NtoCHL)));
}

