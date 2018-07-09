/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/gridmisc.c
 *
 *  \brief Miscellaneous grid routines
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridmisc.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "grid.h"

/** Create a new grid of corner positions from the cell centers.
  */
void grid_centre_to_corner(int nce1, int nce2, double **cx, double **cy,
                           double **gx, double **gy)
{
  int i, j;

  /* Generate interior grid coordinates by interpolation */
  /* Interior points */
  for (j = 1; j < nce2; j++) {
    for (i = 1; i < nce1; i++) {
      gx[j][i] =
        (cx[j][i] + cx[j][i - 1] + cx[j - 1][i] + cx[j - 1][i - 1]) / 4;
      gy[j][i] =
        (cy[j][i] + cy[j][i - 1] + cy[j - 1][i] + cy[j - 1][i - 1]) / 4;
    }
  }
  /* Left and right sides of grid */
  for (j = 1; j < nce2; j++) {
    gx[j][0] = cx[j][0] + cx[j - 1][0] - gx[j][1];
    gy[j][0] = cy[j][0] + cy[j - 1][0] - gy[j][1];
    gx[j][nce1] = cx[j][nce1 - 1] + cx[j - 1][nce1 - 1] - gx[j][nce1 - 1];
    gy[j][nce1] = cy[j][nce1 - 1] + cy[j - 1][nce1 - 1] - gy[j][nce1 - 1];
  }

  /* Bottom and top sides of grid */
  for (i = 1; i < nce1; i++) {
    gx[0][i] = cx[0][i] + cx[0][i - 1] - gx[1][i];
    gy[0][i] = cy[0][i] + cy[0][i - 1] - gy[1][i];
    gx[nce2][i] = cx[nce2 - 1][i] + cx[nce2 - 1][i - 1] - gx[nce2 - 1][i];
    gy[nce2][i] = cy[nce2 - 1][i] + cy[nce2 - 1][i - 1] - gy[nce2 - 1][i];
  }

  /* Corners */
  gx[0][0] = gx[1][0] + gx[0][1] - gx[1][1];
  gy[0][0] = gy[1][0] + gy[0][1] - gy[1][1];
  gx[nce2][0] = gx[nce2][1] + gx[nce2 - 1][0] - gx[nce2 - 1][1];
  gy[nce2][0] = gy[nce2][1] + gy[nce2 - 1][0] - gy[nce2 - 1][1];
  gx[0][nce1] = gx[0][nce1 - 1] + gx[1][nce1] - gx[1][nce1 - 1];
  gy[0][nce1] = gy[0][nce1 - 1] + gy[1][nce1] - gy[1][nce1 - 1];
  gx[nce2][nce1] =
    gx[nce2][nce1 - 1] + gx[nce2 - 1][nce1] - gx[nce2 - 1][nce1 - 1];
  gy[nce2][nce1] =
    gy[nce2][nce1 - 1] + gy[nce2 - 1][nce1] - gy[nce2 - 1][nce1 - 1];
}


/** Expands the grid by one cell all around.
  */
void grid_expand(int nce1, int nce2, double **cx, double **cy, double **gx,
                 double **gy)
{
  int i, j;
  double m, b;

  /* Transfer cell centered points to inner region of the grid. */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      gx[j + 1][i + 1] = cx[j][i];
      gy[j + 1][i + 1] = cy[j][i];
    }
  }

  /* Project out along left and right sides of grid. */
  for (j = 0; j < nce2; j++) {
    m = (cy[j][1] - cy[j][0]) / (cx[j][1] - cx[j][0]);
    b = cy[j][0] - m * cx[j][0];
    gx[j + 1][0] = cx[j][0] - (cx[j][1] - cx[j][0]);
    gy[j + 1][0] = m * gx[j + 1][0] + b;

    m =
      (cy[j][nce1 - 1] - cy[j][nce1 - 2]) / (cx[j][nce1 - 1] -
                                             cx[j][nce1 - 2]);
    b = cy[j][nce1 - 2] - m * cx[j][nce1 - 2];
    gx[j + 1][nce1 + 1] =
      cx[j][nce1 - 1] + (cx[j][nce1 - 1] - cx[j][nce1 - 2]);
    gy[j + 1][nce1 + 1] = m * gx[j + 1][nce1 + 1] + b;
  }

  /* Project out along bottom and top sides of grid. */
  for (i = 0; i < nce1; i++) {
    double dx = cx[1][i] - cx[0][i];
    if (fabs(dx) > 1e-6) {
      m = (cy[1][i] - cy[0][i]) / dx;
      b = cy[0][i] - m * cx[0][i];
      gx[0][i + 1] = cx[0][i] - dx;
      gy[0][i + 1] = m * gx[0][i + 1] + b;
    } else {
      gx[0][i + 1] = cx[0][i];
      gy[0][i + 1] = 2 * cy[0][i] - cy[1][i];
    }

    dx = cx[nce2 - 1][i] - cx[nce2 - 2][i];
    if (fabs(dx) > 1e-6) {
      m = (cy[nce2 - 1][i] - cy[nce2 - 2][i]) / dx;
      b = cy[nce2 - 2][i] - m * cx[nce2 - 2][i];
      gx[nce2 + 1][i + 1] = cx[nce2 - 1][i] + dx;
      gy[nce2 + 1][i + 1] = m * gx[nce2 + 1][i + 1] + b;
    } else {
      gx[nce2 + 1][i + 1] = cx[nce2 - 1][i];
      gy[nce2 + 1][i + 1] = 2 * cy[nce2 - 1][i] - cy[nce2 - 2][i];
    }
  }

  /* Corners */
  gx[0][0] = gx[1][0] + gx[0][1] - gx[1][1];
  gy[0][0] = gy[1][0] + gy[0][1] - gy[1][1];
  gx[nce2 + 1][0] = gx[nce2 + 1][1] + gx[nce2][0] - gx[nce2][1];
  gy[nce2 + 1][0] = gy[nce2 + 1][1] + gy[nce2][0] - gy[nce2][1];
  gx[0][nce1 + 1] = gx[0][nce1] + gx[1][nce1 + 1] - gx[1][nce1];
  gy[0][nce1 + 1] = gy[0][nce1] + gy[1][nce1 + 1] - gy[1][nce1];
  gx[nce2 + 1][nce1 + 1] =
    gx[nce2 + 1][nce1] + gx[nce2][nce1 + 1] - gx[nce2][nce1];
  gy[nce2 + 1][nce1 + 1] =
    gy[nce2 + 1][nce1] + gy[nce2][nce1 + 1] - gy[nce2][nce1];
}
