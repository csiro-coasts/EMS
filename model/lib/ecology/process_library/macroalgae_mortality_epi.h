/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/macroalgae_mortality_epi.h
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
 *  $Id: macroalgae_mortality_epi.h 5972 2018-09-25 05:50:32Z riz008 $
 *
 */

#if !defined(_MACROALGAE_MORTALITY_EPI_H)

void macroalgae_mortality_epi_init(eprocess* p);
void macroalgae_mortality_epi_postinit(eprocess* p);
void macroalgae_mortality_epi_destroy(eprocess* p);
void macroalgae_mortality_epi_precalc(eprocess* p, void* pp);
void macroalgae_mortality_epi_calc(eprocess* p, void* pp);
void macroalgae_mortality_epi_postcalc(eprocess* p, void* pp);

#define _MACROALGAE_MORTALITY_EPI_H
#endif
