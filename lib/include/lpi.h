/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/lpi.h
 *  
 *  Description: A header file for lpi's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lpi.h 6595 2020-09-03 03:36:52Z riz008 $
 *
 */

#if !defined(_LPI_H)
#define _LPI_H

#if !defined(_LPI_STRUCT)
#define _LPI_STRUCT
typedef struct lpi lpi;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"

/* Builds linear interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
lpi* lpi_build(delaunay* d);

/* Destroys linear interpolator.
 *
 * @param l Structure to be destroyed
 */
void lpi_destroy(lpi* l);

/* Finds linearly interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void lpi_interpolate_point(lpi* l, point* p);
/* 3D linear interpolation */
void lpi_interpolate_point3d(lpi* l, point* p);

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void lpi_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void lpi_rebuild(lpi* l, point* p);
/* Rebuilds weights for an individual triangle */
void lpi_rebuild_tri(lpi* l, point* p, int tid, int k);



#endif
