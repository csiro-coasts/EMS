/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/diffusion_heat_epi.h
 *  
 *  Description:
 *  Process implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: diffusion_heat_epi.h 6007 2018-10-30 00:08:59Z bai155 $
 *
 */

#if !defined(_DIFFUSION_HEAT_EPI_H)

void diffusion_heat_epi_init(eprocess* p);
void diffusion_heat_epi_destroy(eprocess* p);
void diffusion_heat_epi_calc(eprocess* p, void* pp);
void diffusion_heat_epi_postcalc(eprocess* p, void* pp);

#define _DIFFUSION_HEAT_EPI_H
#endif
