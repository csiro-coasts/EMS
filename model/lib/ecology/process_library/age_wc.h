/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/age_wc.h
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
 *  $Id: age_wc.h 6056 2019-01-30 02:45:57Z riz008 $
 *
 */

#if !defined(_AGE_WC_H)

void age_wc_init(eprocess* p);
void age_wc_postinit(eprocess* p);
void age_wc_destroy(eprocess* p);
void age_wc_precalc(eprocess* p, void* pp);
void age_wc_calc(eprocess* p, void* pp);

#define _AGE_WC_H
#endif
