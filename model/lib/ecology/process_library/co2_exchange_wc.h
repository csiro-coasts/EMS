/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/ecology/process_library/co2_exchange_wc.h
 *  
 *  Description:
 *  Header file for process co2_exchange
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: co2_exchange_wc.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(CO2_EXCHANGE_H)

void co2_exchange_wc_init(eprocess* p);
void co2_exchange_wc_destroy(eprocess* p);
void co2_exchange_wc_precalc(eprocess* p, void* pp);
void co2_exchange_wc_calc(eprocess* p, void* pp);
void co2_exchange_wc_postcalc(eprocess* p, void* pp);



#define CO2_EXCHANGE_H
#endif
