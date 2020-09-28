/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/coral_growth.c
 *  
 *  Description: Growth terms for Mathieu Mongin's Heron study. Not for use with spectral model.
 *
 *  Reference:
 *
 *  Mongin, M. and M. E. Baird (2014) The interacting effects of photosynthesis, calcification and 
 *   water circulation on carbon chemistry variability on a coral reef flat: a modelling study. Ecol. Mod. 284, 19-34.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: coral_growth.c 5948 2018-09-14 00:30:49Z bai155 $
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
#include "stringtable.h"
#include "cell.h"
#include "column.h"
#include "einterface.h"
#include "coral_growth.h"

  
typedef struct {


/* parameters*/
double p_max_coral;
double p_dark_coral;
double swr_p_k_coral;

double p_max_bma;
double p_dark_bma;
double swr_p_k_bma;

double p_max_alro;
double p_dark_alro;
double swr_p_k_alro;

double swr_g_k_coral;
double swr_g_k_sand;

double g_max_calc;
double k_day_coral;
double k_day_sand;

double k_night_coral;
double k_night_sand;
double g_dark_sand;

double alro_cover;
double bma_cover;
double bommies_cover;
double lcrcslop_cover;
double coral_cover;

double g_max_cal;


/* Tracers */

/*input*/
int CO32_i;   	  /* carbonate */
int HCO3_i;   	  /*bi carbonate*/
int PH_i;      	  /* PH*/
int TEMP_i;	  /*temperature in degree c*/
int SALT_i;	  /* salinity PSU*/
int omega_ca_i;   /* calcite saturation state*/ 
int omega_ar_i;   /* aragonite saturation state*/
 
int bommies_i;    /**/ 
int bma_i;        /*benthic macro algae*/
int lcrcslop_i;    /*corla on the slop */
int coral_i;
int alro_i;        /*rubble sand sediment */

int swr_i;        /* solar radiation ~ long wave*/
int Epilight_i;
int Epilightatt_i;

int p_bommies_i;     /*bommies photo*/ 
int p_bma_i;         /*benthic macro algae photo*/
int p_lcrcslop_i;    /*coral on the slope photo */
int p_alro_i;        /*rubble sand sediment photo */
int p_coral_i;       /*coral  photo */

int g_bommies_i;     /* bommies calci*/ 
int g_lcrcslop_i;    /*corla on the slope calci */
int g_alro_i;        /*rubble sand sediment calci */
int g_coral_i;       /*  coral calci */

/*output */
int gnet_i;
int pnet_i;
int DIC_i;	      /*    INPUT total inorganic carbon (mol/m^3)*/
int ALK_i;	      /*  total alkalinity (eq/m^3) */  




} workspace;


/* This is only called once during the lifetime of a process/ecology. */
/*******************************************************************************/
void coral_growth_init(eprocess* p)
/*******************************************************************************/

{
	ecology* e = p->ecology;
	stringtable* tracers = e->tracers;
	workspace* ws = malloc(sizeof(workspace));
	 stringtable* epis = e->epis;

        int OFFSET_EPI = tracers->n * 2;

	p->workspace = ws;

/*parameters*/
    ws->p_max_coral   = get_parameter_value(e, "p_max_coral");
    ws->p_dark_coral  = get_parameter_value(e, "p_dark_coral");
    ws->swr_p_k_coral = get_parameter_value(e, "swr_p_k_coral");
    
    ws->p_max_bma   = get_parameter_value(e, "p_max_bma");
    ws->p_dark_bma  = get_parameter_value(e, "p_dark_bma");
    ws->swr_p_k_bma = get_parameter_value(e, "swr_p_k_bma");
    
    ws->p_max_alro   = get_parameter_value(e, "p_max_alro");
    ws->p_dark_alro = get_parameter_value(e, "p_dark_alro");
    ws->swr_p_k_alro = get_parameter_value(e, "swr_p_k_alro");
    
    ws->swr_g_k_coral = get_parameter_value(e, "swr_g_k_coral");
    ws->swr_g_k_sand = get_parameter_value(e, "swr_g_k_sand");
    
    ws->k_day_coral   = get_parameter_value(e, "k_day_coral");
    ws->k_day_sand   = get_parameter_value(e, "k_day_sand");
    
    ws->k_night_coral = get_parameter_value(e, "k_night_coral");
    ws->k_night_sand = get_parameter_value(e, "k_night_sand");
    ws->g_dark_sand = get_parameter_value(e, "g_dark_sand");
   

 ws->lcrcslop_cover = get_parameter_value(e, "lcrcslop_cover");
 ws->bma_cover = get_parameter_value(e, "bma_cover");
 ws->bommies_cover= get_parameter_value(e, "bommies_cover");
 ws->alro_cover = get_parameter_value(e, "alro_cover");
 ws->coral_cover = get_parameter_value(e, "coral_cover");

/* TRACERS */
  ws->CO32_i     = e->find_index(tracers,"CO32", e);
  ws->HCO3_i     = e->find_index(tracers,"HCO3", e);
  ws->PH_i       = e->find_index(tracers,"PH", e);
  ws->TEMP_i     = e->find_index(tracers, "temp", e);
  ws->SALT_i     = e->find_index(tracers, "salt", e);
  ws->DIC_i      = e->find_index(tracers,"DIC", e);
  ws->ALK_i      = e->find_index(tracers,"alk", e);
  ws->omega_ca_i = e->find_index(tracers,"omega_ca", e);
  ws->omega_ar_i = e->find_index(tracers,"omega_ar", e);
  
  ws->bommies_i    = e->find_index(epis,"bommies", e)+ OFFSET_EPI;
  ws->bma_i        = e->find_index(epis,"bma", e)+ OFFSET_EPI;
  ws->alro_i       = e->find_index(epis,"alro", e)+ OFFSET_EPI;
  ws->lcrcslop_i   = e->find_index(epis,"lcrcslop", e)+ OFFSET_EPI;
  ws->coral_i      = e->find_index(epis,"coral", e)+ OFFSET_EPI;

  ws->p_bommies_i   = e->find_index(epis,"p_bommies", e)+ OFFSET_EPI;
  ws->p_bma_i       = e->find_index(epis,"p_bma", e)+ OFFSET_EPI;
  ws->p_alro_i      = e->find_index(epis,"p_alro", e)+ OFFSET_EPI;
  ws->p_lcrcslop_i  = e->find_index(epis,"p_lcrcslop", e)+ OFFSET_EPI;
  ws->p_coral_i     = e->find_index(epis,"p_coral", e)+ OFFSET_EPI;

  

  ws->g_bommies_i    = e->find_index(epis,"g_bommies", e)+ OFFSET_EPI;
  ws->g_alro_i       = e->find_index(epis,"g_alro", e)+ OFFSET_EPI;
  ws->g_lcrcslop_i   = e->find_index(epis,"g_lcrcslop", e)+ OFFSET_EPI;
  ws->g_coral_i      = e->find_index(epis,"g_coral", e)+ OFFSET_EPI;

  ws->gnet_i         = e->find_index(epis,"gnet", e)+ OFFSET_EPI;
  ws->pnet_i         = e->find_index(epis,"pnet", e)+ OFFSET_EPI;

  /* ws->swr_i          = e->find_index(tracers,"Light", e);*/

  /*ws->Epilight_i = e->find_index(epis, "Epilight", e) + OFFSET_EPI;*/
  ws->Epilight_i = e->find_index(tracers, "Light", e) ;

 ws->Epilightatt_i = e->find_index(epis, "Epilightatt", e) + OFFSET_EPI;

}
/*******************************************************************************/
void coral_growth_postinit(eprocess* p)
/*******************************************************************************/

{


	

	
}
/*******************************************************************************/
void coral_growth_destroy(eprocess* p)
/*******************************************************************************/
{
	free(p->workspace);

}
/*******************************************************************************/
void coral_growth_precalc(eprocess* p, void* pp)
/*******************************************************************************/
{
	workspace* ws = p->workspace;
	cell* c = (cell*) pp;
	double* y = c->y;

	
	 y[ws->Epilightatt_i] += 0;


}

/*******************************************************************************/
void coral_growth_calc(eprocess* p, void* pp)
/*******************************************************************************/
{ 
workspace* ws = p->workspace;
    intargs* ia = (intargs*) pp;
    cell* c = ((cell*) ia->media);

    double* y = ia->y;
    double* y1 = ia->y1;
    double dz_wc = c->dz_wc;

/*LOCAL DECLARATION*/
    






    double bommies =y[ws->bommies_i];
    double bma     =y[ws->bma_i];
    double alro    =y[ws->alro_i];
    double lcrcslop=y[ws->lcrcslop_i];
    double coral   =y[ws->coral_i];

    double Epilight = y[ws->Epilight_i];
    /* double swr = y[ws->Epilight_i];*/

    double omega_ar=y[ws->omega_ar_i];		
      if (omega_ar < 1)
	omega_ar = 1;

    /* mutliply the coral by its percentage cover    */

  alro=alro*ws->alro_cover;

  bommies=bommies*ws->bommies_cover;

  bma=bma*ws->bma_cover;

  lcrcslop=lcrcslop*ws->lcrcslop_cover;

  coral=coral*ws->coral_cover;



    /* define internal parameters */    
    
    /*local declaration */
    double p0_coral;
    double p0_bma;
    double p0_alro;

    double g0_coral;
    double g0_sand;
    
    double g_dark_coral;
  /*  double g_dark_sand;*/
    
    double g_max_coral;
    double g_max_sand;
    
    /*declare parameterds*/
    
    /*   photosynthesis  organic carbon production*/
    
    /*for coral i.e bommies and lcrcslop anbd coral */
    p0_coral=ws->p_max_coral*(1-exp(-2*Epilight/ws->swr_p_k_coral))-ws->p_dark_coral;
    
    /* for bma macro algae*/
    p0_bma=ws->p_max_bma*(1-exp(-2*Epilight/ws->swr_p_k_bma))-ws->p_dark_bma;

    
    /*for sand and sediment i.e alro */
    p0_alro=ws->p_max_alro*(1-exp(-2*Epilight/ws->swr_p_k_alro))-ws->p_dark_alro;

    
    
    /* calcification or inorganic carbon production*/
    
    /*for coral community bomies lcrcslop coral */
	g_dark_coral=ws->k_night_coral*pow((omega_ar-1),1);
	g_max_coral= ws->k_day_coral*pow((omega_ar-1),1);
	g0_coral=g_max_coral*(1-exp(-2*Epilight/ ws->swr_g_k_coral))+g_dark_coral;
  
    /*for sand and sediment  */

	g_max_sand= ws->k_day_sand*pow((omega_ar-1),1);
	g0_sand=g_max_sand*(1-exp(-2*Epilight/ ws->swr_g_k_sand))+ws->g_dark_sand;

/*********************************************************************/	
/*update the variables  */
/********************************************************************/

/*diagnostic variables*/
/*y1[ws-> gnet_i]         +=86400*(g0_coral*(bommies+lcrcslop+coral)+g0_sand*alro);*/

/*y1[ws-> pnet_i]         +=86400*(p0_coral*(bommies))+86400*(p0_coral*(lcrcslop))+86400*p0_bma*bma+86400*p0_alro*alro+86400*(p0_coral*(coral));**/


/*net rate of photosynthesis for each habitat and total*/
y1[ws-> p_bommies_i]    +=86400*(p0_coral*(bommies));
y1[ws-> p_lcrcslop_i]   +=86400*(p0_coral*(lcrcslop));
y1[ws-> p_bma_i]        +=86400*p0_bma*bma;
y1[ws-> p_alro_i]       +=86400*p0_alro*alro;
y1[ws-> p_coral_i]      +=86400*(p0_coral*(coral));

y1[ws-> pnet_i]  +=(86400*(p0_coral*(bommies)))+(86400*(p0_coral*(lcrcslop)))+(86400*p0_bma*bma)+(86400*p0_alro*alro)+(86400*(p0_coral*(coral)));



/*net rate of calcification for each habitat and total*/
y1[ws-> g_bommies_i]    +=86400*g0_coral*(bommies);
y1[ws-> g_lcrcslop_i]   +=86400*g0_coral*(lcrcslop);
y1[ws-> g_alro_i]       +=86400*g0_sand*alro;
y1[ws-> g_coral_i]    +=86400*g0_coral*(coral);
y1[ws-> gnet_i]       +=(86400*g0_coral*(bommies))+(86400*g0_coral*(lcrcslop))+(86400*g0_sand*alro)+(86400*g0_coral*(coral));



 /**/



/*state variables*/
/* if (dz_wc > 0.1)*/
/*{*/
  y1[ws-> DIC_i] -=(  (p0_coral*bommies)+(p0_coral*lcrcslop)+(p0_coral*coral)+(p0_bma*bma)+(p0_alro*alro)+(g0_coral*bommies)+ (g0_coral*coral)+(g0_coral*lcrcslop) +(g0_sand*alro))/ dz_wc;   


  y1[ws-> ALK_i]  -=2* (   (g0_coral*bommies)+(g0_coral*coral) +(g0_coral*lcrcslop) +(g0_sand*alro)  )/ dz_wc;  
  /* }
 else 
{

   y1[ws-> DIC_i]  -=0;
   y1[ws-> ALK_i]  -=0;
   }*/
 }


void coral_growth_postcalc(eprocess* p, void* pp)
{


    


    
  
    
    
    
    
}


