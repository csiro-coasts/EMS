/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/parameter_info.h
 *  
 *  Description:
 *  Parameter info structure -- header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: parameter_info.h 6303 2019-09-11 00:52:58Z riz008 $
 *
 */

#if !defined(_PARAMETER_INFO_H)

#include "externallibs.h"
#include "declarations.h"

#define MAX_PARAMS (232)

struct parameter_info {
  int index;
  
  /*
   * The parameter name 
   */
  char* name;
  
  /*
   * A short description of the parameter 
   */
  char* desc;

  char* sym;
  
  char* units;
  
  /*
   * Parameter string converted to double 
   */
  double *value;
  int num_values; /* size of this parameter */
  
  double stderr;
  
  /*
   * Parameter string 
   */
  char* stringvalue;
  
  /*
   * Metadata
   */
  char* ref;
};

/** Reads parameter file into parameter_info array.
 * @param biofname bio parameter file
 * @param fp Parameter file
 * @param prefix Prefix, like "GRID0", or NULL
 * @param quitfn Quit error function
 * @param warnfn Warn error function
 * @param emptyfn Empty error function
 * @param nprm Pointer to number of parameters (output)
 * @param parameters Pointer to array of parameter_info (output)
 */
void read_parameter_info(char* fname, FILE *fp, char* prefix, error_fn quitfn, error_fn warnfn, error_fn emptyfn, int* npr, parameter_info* parameters[]);

/** Clears parameter_info array.
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */
void clear_parameter_info(int* nprm, parameter_info* parameters[]);

/** Finds parameter index in a parameter_info array.
 * @param name Parameter name
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 * @return Index id successful; -1 otherwise
 */
int find_parameter_index(const char* name, int nprm, parameter_info pars[]);

/** Transfers parameters defined in d-1 to s-1.
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */
void days2seconds(int nprm, parameter_info parameters[]);

void eco_params_std(parameter_info **parameters, int *nprm);
void eco_params_est(parameter_info **parameters, int *nprm);
void eco_params_gbr4(parameter_info **parameters, int *nprm);
void eco_params_bgc2p0(parameter_info **parameters, int *nprm);
void eco_params_bgc3p1(parameter_info **parameters, int *nprm);

#define _PARAMETER_INFO_H
#endif
