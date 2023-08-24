/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/chlorophyll_normalized_wc.c
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
 *  $Id: chlorophyll_normalized_wc.c 7203 2022-09-15 04:53:37Z bai155 $
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
 *    make sure it is nly read once
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
/* REPLACE filename below with the header file of you process */
#include "chlorophyll_normalized_wc.h"


typedef struct {
	int do_mb;                  /* flag */

	/*
	 * parameters
	 *
	 * This is the working set of state variables, it lives for the duration
	 * between process init and destroy. It does not survive a time step.
	 * Define the list of variables you require for this process i.e

	double DFrad;
	 */


	/*
	 * Define the index of your tracers in the host array here.
	 * The purpose of using indices is to save access time, otherwise
	 * the name has to be searched for every time you intend to access
	 * the variable. Using indices and save them shortens that time considerably.
	 * I.e. - index for the NH4 tracer
	 *
	 int NH4_i;
	 int Tfactor_i;
	 */

/* Parameters */

/* Tracers */
    int Chl_a_i;
    int NormChl_a_i;
    int OD_WeightNorm_i;

/* Common Column tracers */
    int OD_Weight5sum_col_i;


} workspace;


/* This is only called once during the lifetime of a process/ecology. */

void chlorophyll_normalized_wc_init(eprocess* p)
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



/* TRACERS */
  ws->Chl_a_i = e->find_index(tracers,"Chl_a", e);
  ws->NormChl_a_i = e->try_index(tracers,"NormChl_a", e);
  ws->OD_WeightNorm_i = e->try_index(tracers,"OD_WeightNorm", e);
  if(ws->NormChl_a_i < 0)
  	emstag(LERROR,"ecology:chlorophyll_normalized_wc","Tracer 'NormChl_a' not found, disabeling process!");
  
/* Column Tracers */
  ws->OD_Weight5sum_col_i = find_index_or_add(e->cv_column,"OD_Weight5sum", e);


}

void chlorophyll_normalized_wc_postinit(eprocess* p)
{
	ecology* e = p->ecology;
	workspace* ws = p->workspace;

	ws->do_mb = (try_index(e->cv_model, "massbalance_wc", e) >= 0) ? 1 : 0;

	/*
	 * Add code here for calculations or assignments which should be executed after init
	 * but before any calculation. This is only called once during the lifetime of a process/ecology.
	 */


}

void chlorophyll_normalized_wc_destroy(eprocess* p)
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
void chlorophyll_normalized_wc_precalc(eprocess* p, void* pp)
{
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
void chlorophyll_normalized_wc_calc(eprocess* p, void* pp)
{
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
 * Calculate the integrated colour ratio for chlorophyll based on the normalized weight calculated in the process light_wc, 
 * see: Karen Wild-Allen - Observations of diffuse upwelling  irradiance and chloropyll in case 1 waters near the Canary Islands (Spain)
 *  in  "Optics and Laser Technology; Vol. 29;  No. 1; pp 3-8; 1997"
 *
 */
void chlorophyll_normalized_wc_postcalc(eprocess* p, void* pp)
{
  int k=0;
  cell* c = ((cell*) pp);
  workspace* ws = p->workspace;
  column *col = c->col;
  double *cv_od = col->cv[ws->OD_Weight5sum_col_i];

  if(ws->Chl_a_i >= 0 && !isnan(cv_od[0]) && ws->OD_WeightNorm_i > 0 && ws->NormChl_a_i > 0)
  {
  	/* It is expected that if OD_Weight5sum is defined that the OD_WeightNorm is setup across 
  	 * the entire colum.
  	 */
  	
  	col->cells[k]->y[ws->NormChl_a_i]  =  col->cells[k]->y[ws->Chl_a_i]  * col->cells[k]->y[ws->OD_WeightNorm_i] ;
  	k++;
  	for(;k < col->ncells ;k++)
  		col->cells[k]->y[ws->NormChl_a_i]  =  col->cells[k]->y[ws->Chl_a_i]  * col->cells[k]->y[ws->OD_WeightNorm_i]  + col->cells[k-1]->y[ws->Chl_a_i];
  }

}

