/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/dinoflagellate_spectral_mortality_wc.h
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
 *  $Id: dinoflagellate_spectral_mortality_wc.h 6150 2019-03-05 02:00:25Z riz008 $
 *
 */

#if !defined(_DINOFLAGELLATE_SPECTRAL_MORTALITY_WC_H)

void dinoflagellate_spectral_mortality_wc_init(eprocess* p);
void dinoflagellate_spectral_mortality_wc_postinit(eprocess* p);
void dinoflagellate_spectral_mortality_wc_destroy(eprocess* p);
void dinoflagellate_spectral_mortality_wc_precalc(eprocess* p, void* pp);
void dinoflagellate_spectral_mortality_wc_calc(eprocess* p, void* pp);
void dinoflagellate_spectral_mortality_wc_postcalc(eprocess* p, void* pp);

#define _DINOFLAGELLATE_SPECTRAL_MORTALITY_WC_H
#endif
