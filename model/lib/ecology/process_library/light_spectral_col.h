/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_col.h
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
 *  $Id: light_spectral_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_LIGHT_SPECTRAL_COL_H)

void light_spectral_col_init(eprocess* p);
void light_spectral_col_destroy(eprocess* p);
void light_spectral_col_precalc(eprocess* p, void* pp);
void light_spectral_col_postcalc(eprocess* p, void* pp);


#define _LIGHT_SPECTRAL_COL_H
#endif

