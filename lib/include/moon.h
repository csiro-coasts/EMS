/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/moon.h
 *
 *  \brief Prototpyes for moon.c
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: moon.h 6356 2019-10-03 03:46:58Z riz008 $
 */

#if !defined(_MOON_H)
#define _MOON_H

#define r2r(x)  (2.0 * PI * (x - floor(x)))
void moonvars(double jdate, double *rasc, double *decl,
	      double *plon, double *plat, double *dist);

#endif
