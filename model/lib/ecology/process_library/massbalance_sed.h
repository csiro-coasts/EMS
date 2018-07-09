/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/massbalance_sed.h
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
 *  $Id: massbalance_sed.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */


#if !defined(_MASSBALANCE_SED_H)


void massbalance_sed_init(eprocess* p);
void massbalance_sed_destroy(eprocess* p);
void massbalance_sed_precalc(eprocess* p, void* pp);
void massbalance_sed_postcalc(eprocess* p, void* pp);

#define _MASSBALANCE_SED_H
#endif
