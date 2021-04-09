/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/lsql.h
 *  
 *  Description: A header file for lsql's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lsql.h 6597 2020-09-03 05:29:26Z riz008 $
 *
 */
#if !defined(_LSQL_H)
#define _LSQL_H

#if !defined(_LSQL_STRUCT)
#define _LSQL_STRUCT
typedef struct lsql lsql;
#endif

#include "nan.h"
#include "delaunay.h"
#include "grid_utils.h"
#include "svd.h"

#define nol 3 /* Number of polynomial coefficients             */

typedef struct {
    double **B;
    double w[nol];
    int ncells;
    int *cells;
    double tmx;
    double tmn;
    int orf;
    double xr;
    double yr;
    double beta;
} lweights;

/* Builds linear interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
lsql* lsql_build(delaunay* d);

/* Destroys linear interpolator.
 *
 * @param l Structure to be destroyed
 */
void lsql_destroy(lsql* l);

/* Finds linearly interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void lsql_interpolate_point(lsql* l, point* p);

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void lsql_interpolate_points(int nin, point pin[], int nout, point pout[]);

/* Rebuild the wieights given new data */
void lsql_rebuild(lsql* l, point* p);

/* Limiter */
void lsql_limit(lsql* l, int npoint, point* p);

void lsql_minmax(lsql* l, int id, double *mn, double *mx);

#endif
