/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/massbalance_epi.h
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
 *  $Id: massbalance_epi.h 5972 2018-09-25 05:50:32Z riz008 $
 *
 */


#include "eco_constants.h"
#if !defined(_MASSBALANCE_EPI_H)


void massbalance_epi_init(eprocess* p);
void massbalance_epi_postinit(eprocess* p);
void massbalance_epi_destroy(eprocess* p);
void massbalance_epi_precalc(eprocess* p, void* pp);
void massbalance_epi_postcalc(eprocess* p, void* pp);

#define _MASSBALANCE_EPI_H
#endif
