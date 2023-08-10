/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/symbiodinium_spectral_free.h
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
 *  $Id: symbiodinium_spectral_free.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_SYMBIODINIUM_SPECTRAL_FREE_H)

void symbiodinium_spectral_free_init(eprocess* p);
void symbiodinium_spectral_free_postinit(eprocess* p);
void symbiodinium_spectral_free_destroy(eprocess* p);
void symbiodinium_spectral_free_precalc(eprocess* p, void* pp);
void symbiodinium_spectral_free_calc(eprocess* p, void* pp);
void symbiodinium_spectral_free_postcalc(eprocess* p, void* pp);

#define _SYMBIODINIUM_SPECTRAL_FREE_H
#endif
