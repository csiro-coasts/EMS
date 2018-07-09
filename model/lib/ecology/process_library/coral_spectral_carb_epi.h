/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/coral_spectral_carb_epi.h
 *  
 *  Description:
 *  Process heade
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: coral_spectral_carb_epi.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_CORAL_SPECTRAL_CARB_EPI_H)

void coral_spectral_carb_epi_init(eprocess* p);
void coral_spectral_carb_epi_postinit(eprocess* p);
void coral_spectral_carb_epi_destroy(eprocess* p);
void coral_spectral_carb_epi_precalc(eprocess* p, void* pp);
void coral_spectral_carb_epi_calc(eprocess* p, void* pp);
void coral_spectral_carb_epi_postcalc(eprocess* p, void* pp);

#define _CORAL_SPECTRAL_CARB_EPI_H
#endif
