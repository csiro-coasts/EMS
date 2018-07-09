/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/nnpi.h
 *  
 *  Description: A header file for nnpi and nnhpi
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: nnpi.h 5859 2018-07-02 04:04:44Z her127 $
 *
 */

#if !defined(_NNPI_H)
#define _NNPI_H

#if !defined(_NNPI_STRUCT)
#define _NNPI_STRUCT
typedef struct nnpi nnpi;
#endif

#include "delaunay.h"

typedef enum { SIBSON, NON_SIBSONIAN } NN_RULE;

/********/
/* NNPI */
/********/

/* Creates Natural Neighbours point interpolator.
 *
 * @param d Delaunay triangulation
 * @return Natural Neighbours interpolation
 */
nnpi* nnpi_create(delaunay* d);

/* Destroys Natural Neighbours point interpolator.
 *
 * @param nn Structure to be destroyed
 */
void nnpi_destroy(nnpi* nn);

void nnpi_reset(nnpi* nn);

/* Performs Natural Neighbours interpolation in a point.
 *
 * @param nn NN interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void nnpi_interpolate_point(nnpi* nn, point* p);

/* Performs Natural Neighbours interpolation for an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param wmin Minimal allowed weight
 * @param nout Number of output points
 * @param pout Array of output points [nout]
 */
void nnpi_interpolate_points(int nin, point pin[], double wmin, int nout, point pout[]);

/* Sets minimal allowed weight for Natural Neighbours interpolation.
 *
 * For Sibson interpolation, setting wmin = 0 is equivalent to interpolating
 * inside convex hall of the data only (returning NaNs otherwise).
 *
 * @param nn Natural Neighbours point interpolator
 * @param wmin Minimal allowed weight
 */
void nnpi_setwmin(nnpi* nn, double wmin);

/* Gets number of data points involved in current interpolation. For use by
 * `nnai'.
 *
 * @return Number of data points involved in current interpolation
 */
int nnpi_get_nvertices(nnpi* nn);

/* Gets indices of data points involved in current interpolation. For use by
 * `nnai'.
 *
 * @return indices of data points involved in current interpolation
 */
int* nnpi_get_vertices(nnpi* nn);

/* Gets weights of data points involved in current interpolation. For use by
 * `nnai'.
 * @return weights of data points involved in current interpolation
 */
double* nnpi_get_weights(nnpi* nn);

/*
 * Sets the rule
 */
void nnpi_set_rule(nnpi *nn, NN_RULE rule);

/*
 * Sets test_vertice
 */
void nnpi_set_test_vertice(nnpi *nn, int val);

void nnpi_calculate_weights(nnpi* nn, point* p);

/*********/
/* NNHPI */
/*********/

#if !defined(_NNHPI_STRUCT)
#define _NNHPI_STRUCT
typedef struct nnhpi nnhpi;
#endif

/* Creates Natural Neighbours hashing point interpolator.
 *
 * @param d Delaunay triangulation
 * @param size Hash table size (should be of order of number of output points)
 * @return Natural Neighbours interpolation
 */
nnhpi* nnhpi_create(delaunay* d, int size);

/* Destroys Natural Neighbours hashing point interpolation.
 *
 * @param nn Structure to be destroyed
 */
void nnhpi_destroy(nnhpi* nn);
void nnhpi_destroy_weights(nnhpi* nn);

/* Finds Natural Neighbours-interpolated value in a point.
 *
 * @param nnhpi NN point hashing interpolator
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void nnhpi_interpolate(nnhpi* nnhpi, point* p);

/* Modifies interpolated data.
 *
 * Finds point* pd in the underlying Delaunay triangulation such that
 * pd->x = p->x and pd->y = p->y, and copies p->z to pd->z. Exits with error
 * if the point is not found.
 *
 * @param nnhpi Natural Neighbours hashing point interpolator
 * @param p New data
 */
void nnhpi_modify_data(nnhpi* nnhpi, point* p);

/* Sets minimal allowed weight for Natural Neighbours interpolation.
 *
 * For Sibson interpolation, setting wmin = 0 is equivalent to interpolating
 * inside convex hall of the data only (returning NaNs otherwise).
 *
 * @param nn Natural Neighbours point hashing interpolator
 * @param wmin Minimal allowed weight
 */
void nnhpi_setwmin(nnhpi* nn, double wmin);

/*
 * Sets the rule
 */
void nnhpi_set_rule(nnhpi *nn, NN_RULE rule);

#endif
