/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/rect_coord.c
 *
 *  \brief Grid coordinates for rectangular geometries
 *
 *  Routines to calculate grid coordinates for model grids with
 *  rectangular geometries
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: rect_coord.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <math.h>
#include <stdio.h>
#include "grid.h"

#define DEG2RAD		(M_PI/180.0)

/** Calculate coordinates for rectangular grid.
  *
  * @param x where to store grid x values
  * @param y where to store grid y values
  * @param h1 where to store h1 metric values
  * @param h2 where to store h2 metric values
  * @param a1 where to store a1 angle values
  * @param a2 where to store a2 angle values
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param x00 x origin offset
  * @param y00 y origin offset
  * @param rotn angle (degrees) between East and e1 axis
  * @param xinc cell size in x direction
  * @param yinc cell size in y direction
  */
void grid_gen_rect_coord(double **x,  /* where to store grid x values */
                         double **y,  /* where to store grid y values */
                         double **h1, /* where to store h1 metric values */
                         double **h2, /* where to store h2 metric values */
                         double **a1, /* where to store a1 angle values */
                         double **a2, /* where to store a2 angle values */
                         long int nce1, /* number of cells in e1 direction 
                                         */
                         long int nce2, /* number of cells in e2 direction 
                                         */
                         double x00,  /* x origin offset */
                         double y00,  /* y origin offset */
                         double rotn, /* angle (degrees) between East and
                                         e1 axis */
                         double xinc, /* cell size in x direction */
                         double yinc  /* cell size in y direction */
  )
{
  long i, j;
  double xval, yval;
  double sinth;
  double costh;

  sinth = sin(rotn * DEG2RAD);
  costh = cos(rotn * DEG2RAD);
  for (j = 0; j < nce2 + 1; j++) {
    yval = j * yinc;
    for (i = 0; i < nce1 + 1; i++) {
      xval = i * xinc;
      x[j][i] = x00 + xval * costh - yval * sinth;
      y[j][i] = y00 + xval * sinth + yval * costh;
      /***** Analytic calculation deleted
	    h1[j][i] = fabs(xinc);
	    h2[j][i] = fabs(yinc);
	    *******/
    }
  }

  /* Calculate h1 and h2 numerically */
  grid_get_metrics(x, y, nce1, nce2, h1, h2);
  /* calculate a1, a2 numerically */
  grid_get_angle(x, y, nce1, nce2, a1, a2);
}
