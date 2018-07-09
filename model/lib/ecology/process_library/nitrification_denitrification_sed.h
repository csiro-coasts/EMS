/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/nitrification_denitrification_sed.h
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
 *  $Id: nitrification_denitrification_sed.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_NITRIFICATION_DENITRIFICATION_SED_H)

void nitrification_denitrification_sed_init(eprocess* p);
void nitrification_denitrification_sed_postinit(eprocess* p);
void nitrification_denitrification_sed_destroy(eprocess* p);
void nitrification_denitrification_sed_precalc(eprocess* p, void* pp);
void nitrification_denitrification_sed_calc(eprocess* p, void* pp);
void nitrification_denitrification_sed_postcalc(eprocess* p, void* pp);

#define _NITRIFICATION_DENITRIFICATION_SED_H
#endif
