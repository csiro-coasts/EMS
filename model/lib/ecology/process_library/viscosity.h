/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/viscosity.h
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
 *  $Id: viscosity.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_VISCOSITY_H)

void viscosity_init(eprocess* p);
void viscosity_destroy(eprocess* p);
void viscosity_precalc(eprocess* p, void* pp);

#define _VISCOSITY_H
#endif
