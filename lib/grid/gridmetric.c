/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/gridmetric.c
 *
 *  \brief Calculate grid metrics numerically
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridmetric.c 5861 2018-07-02 04:08:04Z her127 $
 */

#include <math.h>
#include <stdio.h>
#include "ems.h"

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (geod_inv_geod_fwd_sodanos(DEG2RAD(x1), DEG2RAD(y1),\
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))

/** Calculate grid metrics numerically.
  *
  * @param x grid x coordinates
  * @param y grid y coordinates
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param h1 storage for h1 values
  * @param h2 storage for h2 values
  */
void grid_get_metrics(double **x, double **y, long int nce1, long int nce2,
                      double **h1, double **h2)
{
  long i, j;
  double dx;
  double dy;

  if (h1 != NULL) {
    /* interior h1 values */
    for (j = 0; j < nce2 + 1; j++)
      for (i = 1; i < nce1; i++) {
        if (isnan(x[j][i])) {
          h1[j][i] = NaN;
        } else if (isnan(x[j][i + 1])) {
          dx = x[j][i] - x[j][i - 1];
          dy = y[j][i] - y[j][i - 1];
          h1[j][i] = hypot(dx, dy);
        } else if (isnan(x[j][i - 1])) {
          dx = x[j][i + 1] - x[j][i];
          dy = y[j][i + 1] - y[j][i];
          h1[j][i] = hypot(dx, dy);
        } else {
          dx = x[j][i + 1] - x[j][i - 1];
          dy = y[j][i + 1] - y[j][i - 1];
          h1[j][i] = hypot(dx, dy) / 2.0;
        }
      }
    /* boundary h1 values */
    for (j = 0; j < nce2 + 1; j++) {
      /* h1 value when i == 0 */
      dx = x[j][1] - x[j][0];
      dy = y[j][1] - y[j][0];
      h1[j][0] = hypot(dx, dy);
      /* h1 value when i == nce1 */
      dx = x[j][nce1] - x[j][nce1 - 1];
      dy = y[j][nce1] - y[j][nce1 - 1];
      h1[j][nce1] = hypot(dx, dy);
    }
  }

  if (h2 != NULL) {
    /* interior h2 values */
    for (j = 1; j < nce2; j++)
      for (i = 0; i < nce1 + 1; i++) {
        if (isnan(x[j][i])) {
          h2[j][i] = NaN;
        } else if (isnan(x[j + 1][i])) {
          dx = x[j][i] - x[j - 1][i];
          dy = y[j][i] - y[j - 1][i];
          h2[j][i] = hypot(dx, dy);
        } else if (isnan(x[j - 1][i])) {
          dx = x[j + 1][i] - x[j][i];
          dy = y[j + 1][i] - y[j][i];
          h2[j][i] = hypot(dx, dy);
        } else {
          dx = x[j + 1][i] - x[j - 1][i];
          dy = y[j + 1][i] - y[j - 1][i];
          h2[j][i] = hypot(dx, dy) / 2.0;
        }
      }

    /* boundary h2 values */
    for (i = 0; i < nce1 + 1; i++) {
      /* h2 value when j == 0 */
      dx = x[1][i] - x[0][i];
      dy = y[1][i] - y[0][i];
      h2[0][i] = hypot(dx, dy);
      /* h2 value when j == nce2 */
      dx = x[nce2][i] - x[nce2 - 1][i];
      dy = y[nce2][i] - y[nce2 - 1][i];
      h2[nce2][i] = hypot(dx, dy);
    }
  }
}


/** Calculate unstructured grid metrics numerically.
  *
  * @param npe number of grid nodes
  * @param x grid x coordinates
  * @param y grid y coordinates
  * @param ns2 number of cells
  * @param h1 storage for edge lengths
  * @param h2 storage for interior lengths
  */
void grid_get_metrics_us(int *npe, double **x, double **y, long int ns2,
                      double **h1, double **h2)
{
  long i, j, ii;
  double dx;
  double dy;

  if (h1 != NULL) {
    /* interior h1 values */
    for (j = 1; j <= ns2; j++)
      for (i = 1; i <= npe[j]; i++) {
        if (isnan(x[j][i])) {
          h1[j][i - 1] = NaN;
        } else {
	  ii = (i + 1 > npe[j]) ? 1 : i + 1;
          dx = x[j][ii] - x[j][i];
          dy = y[j][ii] - y[j][i];
          h1[j][i - 1] = hypot(dx, dy);
        }
      }
  }

  if (h2 != NULL) {
    /* interior h2 values */
    for (j = 1; j <= ns2; j++)
      for (i = 1; i <= npe[j]; i++) {
        if (isnan(x[j][i])) {
          h2[j][i - 1] = NaN;
        } else {
	  ii = (i + 1 > npe[j]) ? 1 : i + 1;
          dx = 0.5 * (x[j][i] + x[j][ii]) - x[j][0];
          dy = 0.5 * (y[j][i] + y[j][ii]) - y[j][0];
          h2[j][i - 1] = hypot(dx, dy);
        }
      }
  }
}


/** Calculate unstructured grid metrics numerically for a geographic grid.
  *
  * @param x grid longitude coordinates (degrees)
  * @param y grid latitude coordinates (degrees)
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param h1 storage for h1 values
  * @param h2 storage for h2 values
  */
void grid_get_geog_metrics(double **x, double **y, int nce1, int nce2,
                           double **h1, double **h2)
{
  int i, j;

  if (h1 != NULL) {
    /* interior h1 values */
    for (j = 0; j < nce2 + 1; j++)
      for (i = 1; i < nce1; i++)
        if (isnan(x[j][i])) {
          h1[j][i] = NaN;
        } else if (isnan(x[j][i + 1])) {
          h1[j][i] = GEODESIC(x[j][i], y[j][i], x[j][i - 1], y[j][i - 1]);
        } else if (isnan(x[j][i - 1])) {
          h1[j][i] = GEODESIC(x[j][i + 1], y[j][i + 1], x[j][i], y[j][i]);
        } else {
          h1[j][i] =
            GEODESIC(x[j][i + 1], y[j][i + 1], x[j][i - 1],
                     y[j][i - 1]) / 2.0;
        }

    /* boundary h1 values */
    for (j = 0; j < nce2 + 1; j++) {
      /* h1 value when i == 0 */
      h1[j][0] = GEODESIC(x[j][1], y[j][1], x[j][0], y[j][0]);
      /* h1 value when i == nce1 */
      h1[j][nce1] = GEODESIC(x[j][nce1], y[j][nce1],
                             x[j][nce1 - 1], y[j][nce1 - 1]);
    }
  }

  if (h2 != NULL) {
    /* interior h2 values */
    for (j = 1; j < nce2; j++)
      for (i = 0; i < nce1 + 1; i++)
        if (isnan(x[j][i])) {
          h2[j][i] = NaN;
        } else if (isnan(x[j + 1][i])) {
          h2[j][i] = GEODESIC(x[j][i], y[j][i], x[j - 1][i], y[j - 1][i]);
        } else if (isnan(x[j - 1][i])) {
          h2[j][i] = GEODESIC(x[j + 1][i], y[j + 1][i], x[j][i], y[j][i]);
        } else {
          h2[j][i] =
            GEODESIC(x[j + 1][i], y[j + 1][i], x[j - 1][i],
                     y[j - 1][i]) / 2.0;
        }

    /* boundary h2 values */
    for (i = 0; i < nce1 + 1; i++) {
      /* h2 value when j == 0 */
      h2[0][i] = GEODESIC(x[1][i], y[1][i], x[0][i], y[0][i]);
      /* h2 value when j == nce2 */
      h2[nce2][i] = GEODESIC(x[nce2][i], y[nce2][i],
                             x[nce2 - 1][i], y[nce2 - 1][i]);
    }
  }
}


/** Calculate unstructured grid metrics numerically.
  *
  * @param npe number of grid nodes
  * @param x grid longitude coordinates (degrees)
  * @param y grid latitude coordinates (degrees)
  * @param ns2 number of cells
  * @param h1 storage for edge lengths
  * @param h2 storage for interior lengths
  */
void grid_get_geog_metrics_us(int *npe, double **x, double **y, long int ns2,
			      double **h1, double **h2)
{
  long i, j, ii;
  double xm;
  double ym;

  if (h1 != NULL) {
    /* interior h1 values */
    for (j = 1; j <= ns2; j++)
      for (i = 1; i <= npe[j]; i++) {
        if (isnan(x[j][i])) {
          h1[j][i - 1] = NaN;
        } else {
	  ii = (i + 1 > npe[j]) ? 1 : i + 1;
          h1[j][i - 1] =
            GEODESIC(x[j][ii], y[j][ii], x[j][i], y[j][i]);
        }
      }
  }

  if (h2 != NULL) {
    /* interior h2 values */
    for (j = 1; j <= ns2; j++)
      for (i = 1; i <= npe[j]; i++) {
        if (isnan(x[j][i])) {
          h2[j][i - 1] = NaN;
        } else {
	  h2[j][i - 1] = 0.0;
	  ii = (i + 1 > npe[j]) ? 1 : i + 1;
          xm = 0.5 * (x[j][i] + x[j][ii]);
          ym = 0.5 * (y[j][i] + y[j][ii]);
	  /* Boundary edges where the centre may lie on the edge midpoint */
	  if (xm != x[j][0] || ym != y[j][0])
	    h2[j][i - 1] =
	      GEODESIC(xm, ym, x[j][0], y[j][0]);
        }
      }
  }
}
