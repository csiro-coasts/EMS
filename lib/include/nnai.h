/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/nnai.h
 *  
 *  Description: A header file for nnai's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: nnai.h 6595 2020-09-03 03:36:52Z riz008 $
 *
 */

#if !defined(_NNAI_H)
#define _NNAI_H

#if !defined(_NNAI_STRUCT)
#define _NNAI_STRUCT
typedef struct nnai nnai;
#endif

#include "delaunay.h"
#include "nnpi.h"

/** Builds Natural Neighbours array interpolator.
 *
 * This includes calculation of weights used in nnai_interpolate().
 *
 * @param d Delaunay triangulation
 * @return Natural Neighbours interpolation
 */
nnai* nnai_build(delaunay* d, int n, double* x, double* y, NN_RULE nn_rule);

/* Destroys Natural Neighbours array interpolator.
 *
 * @param nn Structure to be destroyed
 */
void nnai_destroy(nnai* nn);

/* Conducts NN interpolation in a fixed array of output points using 
 * data specified in a fixed array of input points. Uses pre-calculated
 * weights.
 *
 * @param nn NN array interpolator
 * @param zin input data [nn->d->npoints]
 * @param zout output data [nn->n]. Must be pre-allocated!
 */
void nnai_interpolate(nnai* nn, double* zin, double* zout);

/* Sets minimal allowed weight for Natural Neighbours interpolation.
 *
 * For Sibson interpolation, setting wmin = 0 is equivalent to interpolating
 * inside convex hall of the data only (returning NaNs otherwise).
 *
 * @param nn Natural Neighbours array interpolator
 * @param wmin Minimal allowed weight
 */
void nnai_setwmin(nnai* nn, double wmin);

#endif
