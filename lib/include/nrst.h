/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/nearest.h
 *  
 *  Description: A header file for nearest's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lpi.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_NRST_H)
#define _NRST_H

#if !defined(_NRST_STRUCT)
#define _NRST_STRUCT
typedef struct nrst nrst;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"

/* Builds nearest neighbour interpolator.
 *
 * @param d Delaunay triangulation
 * @return nearest neighbour interpolator
 */
nrst* nrst_build(delaunay* d);

/* Destroys nearest neighbour interpolator.
 *
 * @param l Structure to be destroyed
 */
void nrst_destroy(nrst* n);

/* Finds nearest neighbour interpolated value in a point.
 *
 * @param l nearest neighbour interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void nrst_interpolate_point(nrst* n, point* p);

/* Nearest neighbour interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void nrst_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void nrst_rebuild(nrst* n, point* p);



#endif
