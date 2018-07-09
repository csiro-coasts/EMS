/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/gas_exchange_wc.h
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
 *  $Id: gas_exchange_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(GAS_EXCHANGE_H)

void gas_exchange_wc_init(eprocess* p);
void gas_exchange_wc_destroy(eprocess* p);
void gas_exchange_wc_precalc(eprocess* p, void* pp);
void gas_exchange_wc_calc(eprocess* p, void* pp);
void gas_exchange_wc_postcalc(eprocess* p, void* pp);



#define GAS_EXCHANGE_H
#endif
