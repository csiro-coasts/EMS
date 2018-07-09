/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/anm_epi.h
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
 *  $Id: anm_epi.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_ANM_EPI_H)

void anm_epi_init(eprocess* p);
void anm_epi_destroy(eprocess* p);
void anm_epi_calc(eprocess* p, void* pp);

#define _ANM_EPI_H
#endif
