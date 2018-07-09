/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/values_common_epi.h
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
 *  $Id: values_common_epi.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_VALUES_COMMON_EPI_H)

void values_common_epi_init(eprocess* p);
void values_common_epi_destroy(eprocess* p);
void values_common_epi_precalc(eprocess* p, void* pp);
void values_common_epi_postcalc(eprocess* p, void* pp);

#define _VALUES_COMMON_EPI_H
#endif
