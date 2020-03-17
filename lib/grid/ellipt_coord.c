/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/ellipt_coord.c
 * 
 *  \brief Elliptic geometries
 *
 *  Routines to calculate grid coordinates for
 *  m3d model grids with elliptic geometries
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ellipt_coord.c 6437 2019-11-27 23:58:45Z riz008 $
 */


#include <math.h>
#include <stdio.h>
#include "emsmath.h"
#include "grid.h"

#define	DEG2RAD		(M_PI/180.0)

/** Calculates grid coordinates for model grids with elliptic
  * geometries.
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
  * @param rotn angle (degrees) between East and sigma=1 axis
  * @param ella distance from origin to either focus point
  * @param taumax max tau value ( 1 >= taumax > taumin >= -1 )
  * @param taumin min tau value ( 1 >= taumax > taumin >= -1 )
  * @param nsimm number of symmetric cells in e2 reflected about sigma = 1
  */
void grid_gen_elliptic_coord(double **x,  /* where to store grid x values */
                             double **y,  /* where to store grid y values */
                             double **h1, /* where to store h1 metric
                                             values */
                             double **h2, /* where to store h2 metric
                                             values */
                             double **a1, /* where to store a1 angle
                                             values */
                             double **a2, /* where to store a2 angle
                                             values */
                             long int nce1, /* number of cells in e1
                                               direction */
                             long int nce2, /* number of cells in e2
                                               direction */
                             double x00,  /* x origin offset */
                             double y00,  /* y origin offset */
                             double rotn, /* angle (degrees) between East
                                             and sigma=1 axis */
                             double ella, /* distance from origin to
                                             either focus point */
                             double taumax, /* max tau value ( 1 >= taumax 
                                               > taumin >= -1 ) */
                             double taumin, /* max tau value ( 1 >= taumax 
                                               > taumin >= -1 ) */
                             long int nsimm /* number of symmetric cells
                                               in e2 reflected about sigma 
                                               = 1 */
  )
{
  int i, j;
  double dtau;
  double sinth;
  double costh;
  double th;
  double sigma;
  double tau;
  double xval;
  double yval;
  double v;

  th = rotn * DEG2RAD;
  sinth = sin(th);
  costh = cos(th);

  /* tau spacing */
  dtau = (taumax - taumin) / nce1;
  sigma = 1.0;
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      tau = taumin + i * dtau;
      /* (x,y) before rotation and offset */
      x[j][i] = ella * sigma * tau;
      y[j][i] = ella * sqrt((sigma * sigma - 1) * (1 - tau * tau));
    }
    /* calculate new sigma to make cells roughly square on y axis */
    v = sigma * dtau + sqrt(sigma * sigma - 1);
    sigma = sqrt(v * v + 1);
  }

  /* shift and reflect as necessary for symmetric cells */
  for (j = nce2; j >= nsimm; j--)
    for (i = 0; i < nce1 + 1; i++) {
      x[j][i] = x[j - nsimm][i];
      y[j][i] = y[j - nsimm][i];
    }
  for (j = 0; j < nsimm; j++)
    for (i = 0; i < nce1 + 1; i++) {
      x[j][i] = x[2 * nsimm - j][i];
      y[j][i] = -y[2 * nsimm - j][i];
    }

  /* do rotation and offset */
  for (j = 0; j < nce2 + 1; j++)
    for (i = 0; i < nce1 + 1; i++) {
      xval = x[j][i];
      yval = y[j][i];
      x[j][i] = x00 + xval * costh - yval * sinth;
      y[j][i] = y00 + xval * sinth + yval * costh;
    }

  /* calculate h1 and h2 values numerically */
  grid_get_metrics(x, y, nce1, nce2, h1, h2);
  /* calculate a1 and a2 angle values numerically */
  grid_get_angle(x, y, nce1, nce2, a1, a2);
}
