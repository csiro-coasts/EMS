/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/moldiff.h
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
 *  $Id: moldiff.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_MOLDIFF_H)

void moldiff_init(eprocess* p);
void moldiff_destroy(eprocess* p);
void moldiff_precalc(eprocess* p, void* pp);

#define _MOLDIFF_H
#endif
