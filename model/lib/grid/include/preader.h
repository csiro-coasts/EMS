/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/preader.h
 *  
 *  Description: A header file with preader.c
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: preader.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_PREADER_H)
#define _PREADER_H

#if !defined(_PREADER_STRUCT)
#define _PREADER_STRUCT
typedef struct preader preader;
#endif

preader* preader_create(double xmin, double xmax, double ymin, double ymax, int nx, int ny);
preader* preader_create_from_file(char* fname);
point* preader_getpoint(preader * pr);
void preader_destroy(preader * pr);

#endif
