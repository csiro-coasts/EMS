/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/light_spectral_proto_epi.h
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
 *  $Id: light_spectral_proto_epi.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_LIGHT_SPECTRAL_PROTO_EPI_H)

void light_spectral_proto_epi_init(eprocess* p);
void light_spectral_proto_epi_postinit(eprocess* p);
void light_spectral_proto_epi_destroy(eprocess* p);
void light_spectral_proto_epi_precalc(eprocess* p, void* pp);

#define _LIGHT_SPECTRAL_PROTO_EPI_H
#endif
