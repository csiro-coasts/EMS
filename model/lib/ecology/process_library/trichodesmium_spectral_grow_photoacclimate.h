/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/trichodesmium_spectral_grow_photoacclimate.h
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
 *  $Id: $
 *
 */

#if !defined(_TRICHODESMIUM_SPECTRAL_GROW_PHOTOACCLIMATE_H)

void trichodesmium_spectral_grow_photoacclimate_init(eprocess* p);
void trichodesmium_spectral_grow_photoacclimate_postinit(eprocess* p);
void trichodesmium_spectral_grow_photoacclimate_destroy(eprocess* p);
void trichodesmium_spectral_grow_photoacclimate_precalc(eprocess* p, void* pp);
void trichodesmium_spectral_grow_photoacclimate_calc(eprocess* p, void* pp);
void trichodesmium_spectral_grow_photoacclimate_postcalc(eprocess* p, void* pp);

#define _TRICHODESMIUM_SPECTRAL_GROW_PHOTOACCLIMATE_H
#endif
