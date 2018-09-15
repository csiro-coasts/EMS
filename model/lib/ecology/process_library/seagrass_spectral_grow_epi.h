/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/seagrass_spectral_grow_epi.h
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
 *  $Id: seagrass_spectral_grow_epi.h 5905 2018-08-29 01:14:21Z bai155 $
 *
 */

#if !defined(_SEAGRASS_SPECTRAL_GROW_EPI_H)

void seagrass_spectral_grow_epi_init(eprocess* p);
void seagrass_spectral_grow_epi_postinit(eprocess* p);
void seagrass_spectral_grow_epi_destroy(eprocess* p);
void seagrass_spectral_grow_epi_precalc(eprocess* p, void* pp);
void seagrass_spectral_grow_epi_calc(eprocess* p, void* pp);
void seagrass_spectral_grow_epi_postcalc(eprocess* p, void* pp);

#define _SEAGRASS_SPECTRAL_GROW_EPI_H
#endif
