/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/variable_parameter.h
 *  
 *  Description:
 *  Process header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: variable_parameter.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_VARIABLE_PARAMETER_H)
void  variable_parameter_init(eprocess* p);
void  variable_parameter_destroy(eprocess* p);
void  variable_parameter_precalc(eprocess* p, void* pp);
void  variable_parameter_calc(eprocess* p, void* pp);
void  variable_parameter_postcalc(eprocess* p, void* pp);

#define _VARIABLE_PARAMETER_H
#endif
