/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/baycentric.h
 *  
 *  Description: A header file for baycentric public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: baycentric.h 5859 2018-07-02 04:04:44Z her127 $
 *
 */

#if !defined(_BAL_H)
#define _BAL_H

#if !defined(_BAL_STRUCT)
#define _BAL_STRUCT
typedef struct bal bal;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"
#include "ems.h"

/* Builds baycentric interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
bal* bal_build(delaunay* d);

/* Destroys baycentric interpolator.
 *
 * @param l Structure to be destroyed
 */
void bal_destroy(bal* l);

/* Finds baycentric interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void bal_interpolate_point(bal* l, point* p);
void bal_interpolate_point2(bal* l1, bal* l2, point* p);

/* Baycentric interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void bal_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void bal_rebuild(bal* l, point* p);
void bal_rebuild2(bal* l1, bal* l2, point* p);

#endif
