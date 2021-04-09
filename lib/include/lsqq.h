/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/lsqq.h
 *  
 *  Description: A header file for lsqq's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lsqq.h 6753 2021-04-07 00:56:29Z her127 $
 *
 */
#if !defined(_LSQQ_H)
#define _LSQQ_H

#if !defined(_LSQQ_STRUCT)
#define _LSQQ_STRUCT
typedef struct lsqq lsqq;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"
#include "svd.h"

#define nop 6 /* Number of polynomial coefficients             */

typedef struct {
    double **B;
    double w[nop];
    int ncells;
    int *cells;
    double tmx;
    double tmn;
    int orf;
    double xr;
    double yr;
    double beta;
    double alpha;
} qweights;

/* Builds linear interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
lsqq* lsqq_build(delaunay* d);

/* Destroys linear interpolator.
 *
 * @param l Structure to be destroyed
 */
void lsqq_destroy(lsqq* l);

/* Finds linearly interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void lsqq_interpolate_point(lsqq* l, point* p);

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void lsqq_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void lsqq_rebuild(lsqq* l, point* p);

/* Limiter */
void lsqq_limit(lsqq* l, int npoint, point* p);

void lsqq_minmax(lsqq* l, int id, double *mn, double *mx);

#endif
