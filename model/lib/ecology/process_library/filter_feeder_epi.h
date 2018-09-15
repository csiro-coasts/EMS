/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/filter_feeder_epi.h
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
 *  $Id: filter_feeder_epi.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_FILTER_FEEDER_EPI_H)

void filter_feeder_epi_init(eprocess* p);
void filter_feeder_epi_postinit(eprocess* p);
void filter_feeder_epi_destroy(eprocess* p);
void filter_feeder_epi_precalc(eprocess* p, void* pp);
void filter_feeder_epi_calc(eprocess* p, void* pp);
void filter_feeder_epi_postcalc(eprocess* p, void* pp);

#define _FILTER_FEEDER_EPI_H
#endif
