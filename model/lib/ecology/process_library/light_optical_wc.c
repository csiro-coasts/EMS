/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_optical_wc.c
 *  
 *  Description:
 *  Calculate light (not spectrally resolved) based on optical properties
 *  of Keppel Bay water
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: light_optical_wc.c 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

/**
 * Process Template
 * -  Make sure that all functions used are well documented, functionaly and
 *    scientificly!
 * -  Be aware of processes your process might depend on, even though there are
 *    means to check by name if a process is active - that might not be a good
 *    longterm solution, names may change.
 * -  Keep I/O (reading files) to a minimum or out of here if possible! If you must,
 *    make sure it is only read once
 *
 *
 *
 */
 /* Any key word alled <something>MARKER is used by the wizard and should be
  * removed if you create a process manually.
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
/* REPLACE filename below with the header file of you process */
#include "light_optical_wc.h"

#define EPS_KD 1.0e-10
#define AW_490 0.015
#define MNAP 0.0254046
#define C_NAP 0.0130861
#define BBW_490 0.00157132
#define DTOR 0.01745329

typedef struct {
    /*
     * parameters 
     */
    double k_swr_par;

    /*
     * tracers 
     */
    int Light_i;
    int Kd_i;
    int Chl_a_i;
    int tssf_i;
    int tssuf_i;
    int DOR_C_i;

    /*
     * common variables 
     */
    int ztop_i;
    int lighttop_i;
    int bbp490_i;
} workspace;


/* This is only called once during the lifetime of a process/ecology. */

void light_optical_wc_init(eprocess* p)
{
	ecology* e = p->ecology;
	stringtable* tracers = e->tracers;
	workspace* ws = malloc(sizeof(workspace));

	p->workspace = ws;

	/*
	 * tracers
	 * The Ecology struct will have the correct and save functions assigned
	 * to find_index or try_index, depending on the mode of the host model
	 * these can use a 'lazy' assignment (on a needs basis)
	 *
	 * Assign the tracer index, if you are sure that the tracer exists
	 * or the tracer is mandatory use;
	 * e->find_index(tracers, "<TRACER_NAME>", e);
	 *
	 * if you are not sure and don't want to fail the process and the
	 * program use;
	 * e->try_index(tracers, "<TRACER_NAME>", e);
	 * For NH4 that would look like this;
	 * ws->NH4_i = e->find_index(tracers, "NH4", e);
	 */


    ws->Light_i = e->find_index(tracers, "Light", e);
    ws->Chl_a_i = e->find_index(tracers, "Chl_a", e);
    ws->DOR_C_i = e->find_index(tracers, "DOR_C", e);
    ws->tssf_i = e->find_index(tracers, "TSSF", e);
    ws->tssuf_i = e->find_index(tracers, "TSSUF", e);
    ws->Kd_i = e->find_index(tracers, "Kd_optical", e);
    ws->bbp490_i = e->find_index(tracers, "bbp490", e);


	    /*
	 * Assign the parameter value as defined in 'biofname', if you are sure that the parameter exists
	 * or the variable is mandatory use;
	 * get_parameter_value(e, "<PARAMETER_NAME>", e);
	 *
	 * if you are not sure and don't want to fail the process and the
	 * program use;
	 * try_parameter_value(e, "<PARAMETER_NAME>", e);
	 * For DFrad that would look like this;
	 *
	 * ws->DFrad = get_parameter_value(e, "DFrad", e);
	 *
	 */

    ws->k_swr_par = try_parameter_value(e, "k_swr_par");
    if (isnan(ws->k_swr_par))
        ws->k_swr_par = SWR2PAR;
        

	/*Assign the tracer index for common tracer, common tracers variables are stored separately
	   * and are not integrated
	   * if you are sure that the tracer exists
	   * or the tracer is mandatory use;
	   * find_index(e->cv_cell, "<TRACER_NAME>", e);
	   *
	   * if you are not sure that it exists but it is mandetory don't want to fail the process and the
	   * program use;
	   * find_index_or_add(e->cv_cell, "<TRACER_NAME>", e);
	   *
	   * if you are not sure and don't want to fail the process and the
	   * program use;
	   * try_index(e->cv_cell, "<TRACER_NAME>", e);
	   *
	   * i.e.
	   * ws->Tfactor_i = try_index(e->cv_cell, "Tfactor", e);
	   *
	   * this might fail if the host model uses a 'lazy' creation of tracers.
	   *
	   * To cater for the above scenario you can use
	   * try_index_or_add(e->cv_cell, "<TRACER_NAME>", e);
	   *
	   * i.e.
	   * ws->Tfactor_i = try_index_or_add(e->cv_cell, "Tfactor", e);
	   *
	   * this will test the tracer first and add it if it is provided by the host model.
	   * (this is the savest)
	   */
    ws->lighttop_i = find_index_or_add(e->cv_column, "lighttop", e);
    ws->ztop_i = find_index_or_add(e->cv_column, "ztop", e);
    
    

	/* funtion parameters
	   * the option exists to pass parameters to the function in the 'processes.prm
     * file in the way of <process name >(argument-list)
     * The number of arguments as comma separated list in 'argument-list' must be given
     * as integer in the process entry in the file allprocesses.
     * The parameter can be retrieved lik
     * char* prm = p->prms->se[0]->s
     * for the first on
     * where p->prms returns a 'stringtable' (see stringtable.h
     * and p->prms->se[0] the first entry as 'stringentry' in the stringtable
     * and p->prms->se[0]->s the content of 'stringentry'
     * */


}

void light_optical_wc_postinit(eprocess* p)
{
	ecology* e = p->ecology;

	/*
	 * Add code here for calculations or assignments which should be executed after init
	 * but before any calculation. This is only called once during the lifetime of a process/ecology.
	 */

    if( process_present(e,PT_WC,"light_wc_gradient()")) {
      emstag(LPANIC,"eco:light_optical_properties:postinit","light_wc_gradient seems to be enabled, these processes are mutually exclusive!");
    }
    if( process_present(e,PT_WC,"light_wc()")) {
      emstag(LPANIC,"eco:light_optical_properties:postinit","light_wc seems to be enabled, these processes are mutually exclusive!");
    }


}

void light_optical_wc_destroy(eprocess* p)
{
	free(p->workspace);
	    /*
	 * Add code here for releasing any additional system resources you have requested.
	 * It is unlikely for standart processes that you would and it is expected
	 * that you know what you are doing if so.
	 *
	 * This is only called once during the lifetime of a process/ecology.
	 */

}

/* Add here the functions you like to be executed before the integration step.
 * This function is only executed once per time step.
 * ws   - pointer to your workspace defined at the begining of this file
 * cell - pointer to the current cell
 * y    - cell based tracer values
 * cv   - common tracer array as exchange between processes
 *        These are the tracers which might need to be adjusted but are
 *        not integrated
 *
 */
/**
 * Documentation:
 * Add context and scientific documentation here!
 *
 */
void light_optical_wc_precalc(eprocess* p, void* pp)
{
	workspace* ws = p->workspace;
	cell* c = (cell*) pp;
	column* col = c->col;
	ecology* e = p->ecology;
	double* y = c->y;
	double dz = c->dz_wc;
	double lighttop, lightbot, ztop, zbot;
	double ay440, bb_490, k1, k2, at_490, Kd;
	double TSS = y[ws->tssf_i] + y[ws->tssuf_i];
//BJR 2008/1/21
//	double DOR_C = y[ws->DOR_C_i];
	double DOR_C = y[ws->DOR_C_i]/1000;
//END BJR 2008/1/21
	double Chl_a = y[ws->Chl_a_i];
	
	double zenith = e->zenith;
	//    double zenith = 90.0;
	
	double ay_490 = 0.0;
	double bbp_490;
	
	double aph_490 = 0.0331*(pow(Chl_a,0.7574));  /* relationship from FEKB 2003 spectrophotometer data aph490_CHL.XLS */
	double Cnap = TSS*1000.0; /* TSS here does not include Chl_a and is in kg/m3.  We could possibly add detritus to TSS here. */

	double anap_490 = (Cnap > 0.7) ? MNAP*Cnap + Cnap : (Cnap > 0) ? MNAP*Cnap : 0.0;
	
	/*
	Conversion of DOC to ay440
	Kadija Oubelkheir et.al., (2004) Bio-optical characterisation of
	Australian ocean,coastal and estuarine waters, Ocean Optics conference Freemantle Australia
	*/
	if (DOR_C>0.56) 
	{
	  ay440 = (DOR_C-0.56)/4.38;
	  
	   /* Derived from 2003 Keppel Bay data */
	  ay_490 = 0.4173*ay440 + 0.0027;
	}
	if (isnan((col->cv[ws->ztop_i])[0]))
	    (col->cv[ws->ztop_i])[0] = 0;
	ztop = (col->cv[ws->ztop_i])[0];
	zbot = ztop + dz;
	
	at_490 = AW_490 + aph_490 + anap_490 + ay_490;
        bbp_490 = (TSS > 0.0) ? 36.5379 * TSS + 0.0310525 : 0.0; /* from FEKB 2003 ac9 data */

	bb_490 = BBW_490 + bbp_490; /* total backscattering */
	
	/* from Zong-Ping Lee JGR v110 c09019 */
	k1=(-0.057+0.482*(pow(at_490,0.5))+(4.221*bb_490))*(1+0.090*sin(zenith*DTOR));    
	k2=(0.183+(0.702*at_490)+(-2.567*bb_490))*(1.465+(-0.667)*cos(zenith*DTOR)); 
	Kd=k1+(k2/(pow((1+(ztop+zbot)/2),0.5)));              
	
	if (Kd<=0.0)
	   Kd = EPS_KD;
	
	if (isnan((col->cv[ws->lighttop_i])[0]))
	    (col->cv[ws->lighttop_i])[0] = einterface_getlighttop(col->model, col->b) * ws->k_swr_par;
	lighttop = (col->cv[ws->lighttop_i])[0];
	
	lightbot = lighttop * exp(-Kd * dz);
	y[ws->Light_i] = (lighttop - lightbot) / Kd / dz;
	
	(col->cv[ws->lighttop_i])[0] = lightbot;
	(col->cv[ws->ztop_i])[0] = zbot;
	y[ws->bbp490_i] = bbp_490;
	y[ws->Kd_i] = Kd;
	
}


/* Add here the functions you like to be executed during the integration step.
 * This function may be executed several times during a time step. Do not be
 * concerned with absolute values, use values from y to initialise and y1 as
 * target for your calculations.
 *
 * ws   - pointer to your workspace defined at the begining of this file
 * cell - pointer to the current cell
 * y    - working array which lifes ONLY over the calculation
 *        These reflect the initialisation values for calculations to
 *        be executed in this call of the function.
 * y1   - working array which lifes ONLY over the calculation
 *        Apply the result of the calculation to this array of the tracers!
 * cv   - common tracer array as exchange between processes
 *        These are the tracers which might need to be adjusted but are
 *        not integrated
 *
 */

/**
 * Documentation:
 * Add context and scientific documentation here!
 *
 */
void light_optical_wc_calc(eprocess* p, void* pp)
{

	/*Assign initial tracer values which is integrated i.e.
	double NH4 = y[ws->NH4_i]; */

	/*Assign a common tracer value which is not integrated
	double Tfactor = cv[ws->Tfactor_i];*/

	/*Assign a working variable
	double DFrad = ws->DFrad; */


	/* to calculate some value updating the change
	 * i.e.
	 * y1[ws->NH4_i] += NH4 * Tfactor / DFrad;
	 */

}


/* Add here the functions you like to be executed after the integration step.
 * This function is only executed once per time step.
 * ws   - pointer to your workspace defined at the begining of this file
 * cell - pointer to the current cell
 * y    - cell based tracer values
 * cv   - common tracer array as exchange between processes
 *        These are the tracers which might need to be adjusted but are
 *        not integrated
 *
 */
/**
 * Documentation:
 * Add context and scientific documentation here!
 *
 */
void light_optical_wc_postcalc(eprocess* p, void* pp)
{

}

