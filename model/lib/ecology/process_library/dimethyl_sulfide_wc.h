/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dimethyl_sulfide_wc.h
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
 *  $Id: dimethyl_sulfide_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(DIMETHYL_SULFIDE_H)

void dimethyl_sulfide_wc_init(eprocess* p);
void dimethyl_sulfide_wc_destroy(eprocess* p);
void dimethyl_sulfide_wc_precalc(eprocess* p, void* pp);
void dimethyl_sulfide_wc_calc(eprocess* p, void* pp);
void dimethyl_sulfide_wc_postcalc(eprocess* p, void* pp);



#define DIMETHYL_SULFIDE_H
#endif
