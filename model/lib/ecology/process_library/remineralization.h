/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/remineralization.h
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
 *  $Id: remineralization.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_REMINERALIZATION_H)

void remineralization_init(eprocess* p);
void remineralization_postinit(eprocess* p);
void remineralization_destroy(eprocess* p);
void remineralization_precalc(eprocess* p, void* pp);
void remineralization_calc(eprocess* p, void* pp);
void remineralization_postcalc(eprocess* p, void* pp);

#define _REMINERALIZATION_H
#endif
