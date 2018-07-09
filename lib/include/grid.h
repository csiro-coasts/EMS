/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/grid.h
 *
 *  \brief Protoyypes for grid coordinates calculation routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: grid.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#ifndef _GRID_H
#define _GRID_H

#include "xytoij.h"

void grid_gen_rect_coord(double **x, double **y, double **h1, double **h2,
                         double **a1, double **a2, long nce1, long nce2,
                         double x00, double y00, double rotn, double xinc,
                         double yinc);

void grid_gen_polar_coord(double **x, double **y, double **h1, double **h2,
                          double **a1, double **a2, long nce1, long nce2,
                          double x00, double y00, double rotn, double arc,
                          double rmin);

void grid_gen_elliptic_coord(double **x, double **y, double **h1,
                             double **h2, double **a1, double **a2,
                             long nce1, long nce2, double x00, double y00,
                             double rotn, double ella, double taumax,
                             double taumin, long nsimm);

void grid_get_metrics(double **x, double **y, long nce1, long nce2,
                      double **h1, double **h2);

void grid_get_geog_metrics(double **x, double **y, int nce1, int nce2,
                           double **h1, double **h2);

void grid_get_angle(double **x, double **y, long nce1, long nce2,
                    double **a1, double **a2);

void grid_get_geog_angle(double **x, double **y, int nce1, int nce2,
                         double **a1, double **a2);


void grid_centre_to_corner(int nce1, int nce2, double **cx, double **cy,
                           double **gx, double **gy);

void grid_expand(int nce1, int nce2, double **cx, double **cy,
                 double **gx, double **gy);

/* Unstructured versions */
void grid_get_metrics_us(int *npe, double **x, double **y, long int ns2,
			 double **h1, double **h2);
void grid_get_geog_metrics_us(int *npe, double **x, double **y, long int ns2,
			      double **h1, double **h2);
void grid_get_angle_us(int *npe, double **x, double **y, long int ns2,
		       double **a1, double **a2);
void grid_get_geog_angle_us(int *npe, double **x, double **y, long int ns2,
			    double **a1, double **a2);

#endif
