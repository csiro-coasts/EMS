/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/gridaverager.h
 *  
 *  Description: Calculate average value of a 2D field a cell of a grid
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridaverager.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_GRIDAVERAGER_H)
#define _GRIDAVERAGER_H

struct gridaverager;
typedef struct gridaverager gridaverager;

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
    double x;
    double y;
    double z;
} point;
#endif

#include "gridmap.h"

gridaverager* ga_create(gridmap* gm);
void ga_destroy(gridaverager* ga);
void ga_addpoints(gridaverager* ga, int n, point points[]);
void ga_getvalue(gridaverager* ga, point * p);

#endif
