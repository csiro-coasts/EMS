/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/minell.h
 *  
 *  Description: A header for the minimal ellipse stuff
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: minell.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_MINELL_H)
#define _MINELL_H

#if !defined(_MINELL_STRUCT)
#define _MINELL_STRUCT
typedef struct minell minell;
#endif

#include "grid_utils.h"

/* Note that minell_build() shuffles the input point array */
minell* minell_build(int n, point p[]);
void minell_destroy(minell* me);
void minell_scalepoints(minell* me, int n, point p[]);
void minell_rescalepoints(minell* me, int n, point p[]);

#endif
