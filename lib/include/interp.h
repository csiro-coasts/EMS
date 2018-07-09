/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/interp.h
 *
 *  \brief Header for 1d lookup table interpolation
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: interp.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#if !defined(_INTERP_H)
#define _INTERP_H

void interp1d(double *X, double *Y, int NX, double *X1, double *Y1, int NX1);

#endif
