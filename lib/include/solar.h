/**
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/solar.h
 *
 *  \brief Prototpyes for solar elevation calculation
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: solar.h 5900 2018-08-28 02:09:26Z riz008 $
 */

#if !defined(_SOLAR_H)
#define _SOLAR_H
void dtime(char *ounit, char *iunit, double time, int *year, double *day, double *lon);
double calc_solar_elevation(char *ounit, char *tunit, double time, double lat,
			    double *out_dec, double *lon);
#endif
