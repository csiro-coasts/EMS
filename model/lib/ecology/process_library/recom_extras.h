/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/recom_extras.h
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
 *  $Id: recom_extras.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_RECOM_EXTRAS_H)

void recom_extras_init(eprocess* p);
void recom_extras_postinit(eprocess* p);
void recom_extras_destroy(eprocess* p);
void recom_extras_postcalc(eprocess* p, void* pp);

#define _RECOM_EXTRAS_H
#endif
