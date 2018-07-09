/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/phytoplankton_diel_grow_wc.h
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
 *  $Id: phytoplankton_diel_grow_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_PHYTOPLANKTON_DIEL_GROW_WC_H)

void phytoplankton_diel_grow_wc_init(eprocess* p);
void phytoplankton_diel_grow_wc_postinit(eprocess* p);
void phytoplankton_diel_grow_wc_destroy(eprocess* p);
void phytoplankton_diel_grow_wc_precalc(eprocess* p, void* pp);
void phytoplankton_diel_grow_wc_calc(eprocess* p, void* pp);
void phytoplankton_diel_grow_wc_postcalc(eprocess* p, void* pp);

#define _PHYTOPLANKTON_DIEL_GROW_WC_H
#endif
