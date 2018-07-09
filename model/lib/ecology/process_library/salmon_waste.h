/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/salmon_waste.h
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
 *  $Id: salmon_waste.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_SALMON_WASTE_H)

void salmon_waste_init(eprocess* p);
void salmon_waste_destroy(eprocess* p);
void salmon_waste_precalc(eprocess* p, void* pp);
void salmon_waste_postcalc(eprocess* p, void* pp);

#define _SALMON_WASTE_H
#endif
