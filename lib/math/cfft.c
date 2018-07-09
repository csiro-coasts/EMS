/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/cfft.c
 *
 *  \brief FFT routine
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: cfft.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "ems.h"

/** FFT routine similar to four1() in "Numerical Recipes
  * in C" by Press et al.
  *
  * Here more error checking is done and arrays are zero based.
  *
  * @param data an array of ndata complex numbers, stored with
  * real and imaginary parts in adjacent locations.
  * @param ndata Number of data points - must be a power of 2
  * @param dirn If 1, the forward transform is done, if -1
  * the reverse transform is done (division by
  * ndata is not done by the routine in this case)
  * 
  * @author Stephen Walker.
  */
void cfft(double *data, int ndata, int dirn)
{
  int mmax;
  int i, j, m, n;
  int istep;
  double wtemp;
  double wr;
  double wpr;
  double wpi;
  double wi;
  double theta;
  double tempr;
  double tempi;

  /* NOTE - NR seems to have the transform direction backwards! */
  dirn = -dirn;

  /* Check that ndata is an integer power of 2 */
  n = ndata;
  i = 1;
  while (n >>= 1)
    i <<= 1;
  if (i != ndata)
    quit("cfft: Number of points must be a power of 2\n");

  n = ndata << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      double tmp;
      tmp = data[i - 1];
      data[i - 1] = data[j - 1];
      data[j - 1] = tmp;
      tmp = data[i];
      data[i] = data[j];
      data[j] = tmp;
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax = 2;
  while (n > mmax) {
    istep = 2 * mmax;
    theta = 2.0 * M_PI / (dirn * mmax);
    wtemp = sin(theta / 2.0);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j - 1] - wi * data[j];
        tempi = wr * data[j] + wi * data[j - 1];
        data[j - 1] = data[i - 1] - tempr;
        data[j] = data[i] - tempi;
        data[i - 1] += tempr;
        data[i] += tempi;
      }
      wtemp = wr;
      wr = wr * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}
