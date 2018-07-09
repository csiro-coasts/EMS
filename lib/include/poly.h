/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/poly.h
 *
 *  \brief Library routines for polylines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: poly.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _POLY_H
#define _POLY_H

#include <stdio.h>

typedef struct {
  double x;
  double y;
} point_t;

typedef struct {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
} extent_t;

typedef struct {
  int n;                        /* number of points */
  int nallocated;               /* number of allocated points */
  extent_t e;                   /* bounding rectangle */
  double *x;                    /* array of x coordinates [n] */
  double *y;                    /* array of y coordinates [n] */
} poly_t;

poly_t *poly_create();
void poly_destroy(poly_t *pl);

void poly_add_point(poly_t *pl, double x, double y);
void poly_add_points(poly_t *pl, int n, double x[], double y[]);
void poly_add_point_at(poly_t *pl, int index, double x, double y);
void poly_append(poly_t *pl1, poly_t *pl2);
double poly_area(poly_t *pl);
void poly_clear(poly_t *pl);
void poly_close(poly_t *pl);
int poly_contains_point(poly_t *pl, double x, double y);
poly_t *poly_copy(poly_t *pl);
void poly_delete_point(poly_t *pl, int index);
void poly_despike(poly_t *pl, double maxdist);
int poly_find_index(poly_t *pl, double x, double y);
int poly_is_closed(poly_t *pl, double eps);
int poly_read(poly_t *pl, FILE * fp);
void poly_resample(poly_t *pl, double eps);
void poly_reverse(poly_t *pl);
poly_t *poly_smooth(poly_t *pl, int ns);
void poly_write(poly_t *pl, FILE * fp);

poly_t* poly_formbound(int nce1, int nce2, double** x, double** y);
poly_t* poly_formboundij(int nce1, int nce2, double** x);
void poly_compact(poly_t* pl, double eps);

#endif
