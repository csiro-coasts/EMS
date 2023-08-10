/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/autotracer.h
 *  
 *  Description:
 *  Primary sparse structures.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sparse.h 6886 2021-07-29 00:46:39Z her127 $
 *
 */

#if !defined(_AUTOTRACER_H)
#define _AUTOTRACER_H

struct autotracer {
  char *identifier;
  char *name;
  char *long_name;
  char *units;
  double valid_range[2];
  double fill_value;
  int type;
  int inwc;
  int dissol;
  int advect;
  int diffuse;
  int diagn;
};

tracer_info_t autotracer[] = {
  { .identifier = "TAG",
    .name = "salt",
    .long_name = "Salinity",
    .units = "PSU",
    .valid_range_0 = 0,
    .valid_range_1 = 40,
    .fill_value = 35,
    .type = WATER|HYDRO|PROGNOSTIC,
    .inwc = 1,
    .dissol = 1,
    .advect = 1,
    .diffuse = 1,
    .diagn = 0 },

  { .identifier = "TAG",
    .name = "temp",
    .long_name = "Temperature",
    .units = "degrees C",
    .valid_range_0 = -4,
    .valid_range_1 = 40,
    .fill_value = 20,
    .type = WATER|HYDRO|PROGNOSTIC,
    .inwc = 1,
    .dissol = 1,
    .advect = 1,
    .diffuse = 1,
    .diagn = 0 }

};

#endif
