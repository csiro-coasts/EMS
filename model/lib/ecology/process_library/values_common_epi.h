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
 *  $Id: values_common_epi.h 6905 2021-09-07 07:03:11Z bai155 $
 *
 */

#if !defined(_VALUES_COMMON_EPI_H)

void values_common_epi_init(eprocess* p);
void values_common_epi_postinit(eprocess* p);
void values_common_epi_destroy(eprocess* p);
void values_common_epi_precalc(eprocess* p, void* pp);
void values_common_epi_postcalc(eprocess* p, void* pp);

#define _VALUES_COMMON_EPI_H
#endif
