/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/filter_feeder_wc.h
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
 *  $Id: filter_feeder_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_FILTER_FEEDER_WC_H)

void filter_feeder_wc_init(eprocess* p);
void filter_feeder_wc_postinit(eprocess* p);
void filter_feeder_wc_destroy(eprocess* p);
void filter_feeder_wc_precalc(eprocess* p, void* pp);
void filter_feeder_wc_calc(eprocess* p, void* pp);
void filter_feeder_wc_postcalc(eprocess* p, void* pp);

#define _FILTER_FEEDER_WC_H
#endif
