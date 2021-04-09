/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/bilinear.h
 *  
 *  Description: A header file for bilinear's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: bilinear.h 6597 2020-09-03 05:29:26Z riz008 $
 *
 */

#if !defined(_BL_H)
#define _BL_H

#if !defined(_BL_STRUCT)
#define _BL_STRUCT
typedef struct bl bl;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"

/* Builds bi-linear interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
bl* bl_build(delaunay* d);

/* Destroys linear interpolator.
 *
 * @param l Structure to be destroyed
 */
void bl_destroy(bl* l);

/* Finds linearly interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void bl_interpolate_point(bl* l, point* p);
void bl_interpolate_point2(bl* l1, bl* l2, point* p);

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void bl_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void bl_rebuild(bl* l, point* p);
void bl_rebuild2(bl* l1, bl* l2, point* p);


#endif
