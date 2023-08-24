/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/mixed_layer_age_col.h
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
 *  $Id: mixed_layer_age_col.h 6056 2019-01-30 02:45:57Z riz008 $
 *
 */

#if !defined(_MIXED_LAYER_AGE_COL_H)

void mixed_layer_age_col_init(eprocess* p);
void mixed_layer_age_col_postinit(eprocess* p);
void mixed_layer_age_col_destroy(eprocess* p);
void mixed_layer_age_col_precalc(eprocess* p, void* pp);
void mixed_layer_age_col_postcalc(eprocess* p, void* pp);
#define _MIXED_LAYER_AGE_COL_H
#endif
