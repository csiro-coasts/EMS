/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/interp.c
 *
 *  \brief Basic interpolation routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: interp.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>

/** interp1
 *
 * 1D linear interpolation. Finds Y1 at X1 given X & Y
 *
 * The output vector Y1 *must* be appropriately be pre-allocated
 *
 * @param X  Input X vector
 * @param Y  Input Y vector
 * @param NX Length of X and Y
 * @param X1 Output X vector to interpolate onto
 * @param Y1 Output Y interpolated values
 * @param NX1 Length of output vectors
 */
void interp1d(double *X, double *Y, int NX, double *X1, double *Y1, int NX1)
{
  int k, kk;

  /*
   * Solve for Y1
   */
  for (k=0; k<NX1; k++){
    // Zero order hold for x-values below the range
    if (X1[k] <= X[0])
      Y1[k] = Y[0];
    // Zero order hold for x-values above the range
    else if (X1[k] >= X[NX-1])
      Y1[k] = Y[NX-1];
    // Do the interpolation
    else {
      for (kk=0; kk<NX-1; kk++) {
	// Find the two-point bracket
	if (X1[k] >= X[kk] && X1[k] < X[kk+1]) {
	  Y1[k] = (Y[kk] * (X1[k] - X[kk+1]) + Y[kk+1] * (X[kk] - X1[k]))
	                                / (X[kk] - X[kk+1]);
	}
      }
    }
  }
}


