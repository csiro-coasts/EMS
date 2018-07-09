/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/nan.c
 *
 *  \brief NaN definition
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: nan.c 5833 2018-06-27 00:21:35Z riz008 $
 */

/* UR
 * alternative together with nan.h 
 * the origenal version caused failure on itanium2 for ecology
 * tested against gcc on suse 7.x - 9.2, icc on Redhat E-64, Suse E-64 Altix 
 */

#include "nan.h"

double setNaN()
{
  return NaN; 
}

