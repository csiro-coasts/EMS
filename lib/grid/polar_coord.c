/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/polar_coord.c
 *
 *  \brief Grid coordinates for polar geometries
 *
 *  Routines to calculate grid coordinates for model grids with
 *  polar geometries
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: polar_coord.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <math.h>
#include <stdio.h>
#include "grid.h"

#define DEG2RAD		(M_PI/180.0)

/** Calculate grid coordinates for model grids with polar
  * coordinates.
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
  * @param rotn angle (degrees) between West and e1=0 radial
  * @param arc angular size of grid from origin
  * @param rmin  minimum grid radius
  */
void grid_gen_polar_coord(double **x, /* where to store grid x values */
                          double **y, /* where to store grid y values */
                          double **h1,  /* where to store h1 metric values 
                                         */
                          double **h2,  /* where to store h2 metric values 
                                         */
                          double **a1,  /* where to store a1 angle values */
                          double **a2,  /* where to store a2 angle values */
                          long int nce1,  /* number of cells in e1
                                             direction */
                          long int nce2,  /* number of cells in e2
                                             direction */
                          double x00, /* x origin offset */
                          double y00, /* y origin offset */
                          double rotn,  /* angle (degrees) between East
                                           and e1=0 radial */
                          double arc, /* angular size of grid from origin */
                          double rmin /* minimum grid radius */
  )
{
  int i, j;
  double fac;
  double r;
  double th0;
  double th;
  double dth;

  /* convert rotation to radians from x axis here rotn == 0.0 means the
     arc starts on -x axis */
  th0 = rotn * DEG2RAD;
  if (th0 > 2.0 * M_PI)
    th0 -= 2.0 * M_PI;
  if (th0 < 0.0)
    th0 += 2.0 * M_PI;
  dth = (arc / nce1) * DEG2RAD;

  /* initialise r and r scaling */
  r = rmin;
  fac = exp(dth);
  for (j = 0; j < nce2 + 1; j++) {
    th = th0;
    for (i = 0; i < nce1 + 1; i++) {
      x[j][i] = x00 + r * cos(th);
      y[j][i] = y00 + r * sin(th);
      /****** Analytic calculation deleted
	    h1[j][i] = fabs(r*dth);
	    h2[j][i] = fabs(r*dth);
	    *********/
      th -= dth;
      if (th > 2.0 * M_PI)
        th -= 2.0 * M_PI;
      if (th < 0.0)
        th += 2.0 * M_PI;
    }
    r *= fac;
  }

  /* Calculate h1 and h2 numerically */
  grid_get_metrics(x, y, nce1, nce2, h1, h2);
  /* calculate a1, a2 numerically */
  grid_get_angle(x, y, nce1, nce2, a1, a2);
}
