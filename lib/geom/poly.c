/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/geom/poly.c
 *
 *  \brief Library routines for polylines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: poly.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <values.h>
#include "poly.h"

#define PL_NSTART 4

#if !defined(MAXLINELEN)
#define MAXLINELEN 2048
#endif

#if !defined(max)
#define max(x,y) ((x)>(y) ? (x) : (y) )
#endif

#if !defined(min)
#define min(x,y) ((x)<(y) ? (x) : (y) )
#endif

void quit(char *format, ...);
int prm_next_line(char *line, long n, FILE * fp);
void spline(double *x, double *y, long n, int derivspec,
            double start_deriv, double end_deriv, double *ydd);
void spline_interp(double *xa, double *ya, double *ydd, long n, double x,
                   double *y);

/* Static functions */

/** Clears extent_t.
 * @param e Extent
 */
static void extent_clear(extent_t *e)
{
  e->xmin = MAXDOUBLE;
  e->ymin = MAXDOUBLE;
  e->xmax = -MAXDOUBLE;
  e->ymax = -MAXDOUBLE;
}

/** Checks whether a point belongs to extent_t.
 * @param e Extent.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @return 1 for yes, 0 for no.
 */
static int extent_contains_point(extent_t *e, double x, double y)
{
  if (x >= e->xmin && x <= e->xmax && y >= e->ymin && y <= e->ymax)
    return 1;

  return 0;
}

/** Updates extent_t to include a specified point.
 * @param e Extent.
 * @param x X coordinate.
 * @param y Y coordinate.
 */
static void extent_update(extent_t *e, double x, double y)
{
  if (x < e->xmin)
    e->xmin = x;
  if (y < e->ymin)
    e->ymin = y;
  if (x > e->xmax)
    e->xmax = x;
  if (y > e->ymax)
    e->ymax = y;
}


static int ononeline(double* x, double* y, int i1, int i2, int i3, double eps)
{
    return fabs((x[i1] - x[i2]) * (y[i3] - y[i2]) - 
		            (x[i3] - x[i2]) * (y[i1] - y[i2])) <= eps;
}


/** Tests whether two points coincide. A threshold distance is used to
 * define the tolerance within which two points may be considered the same.
 *
 * @param p1 first point
 * @param p2 second point
 * @param eps threshold distance.
 * @return non-zero if successful.
 */
static int point_equal(double x1, double y1, double x2, double y2,
                       double eps)
{
  return ((fabs(x1 - x2) <= eps) && (fabs(y1 - y2) <= eps));
}

/** Re-calculates extent_t of a polyline.
 * @param pl Polyline.
 * @param x X coordinate.
 * @param y Y coordinate.
 */
static void recalc_extent(poly_t *pl)
{
  double *xs = pl->x;
  double *ys = pl->y;
  extent_t *e = &pl->e;
  int n = pl->n;
  int i;

  extent_clear(e);
  for (i = 0; i < n; ++i)
    extent_update(e, xs[i], ys[i]);
}

/** Appends a point to the tail of a polyline.
 * @param pl Polyline.
 * @param x X coordinate.
 * @param y Y coordinate.
 */
void poly_add_point(poly_t *pl, double x, double y)
{
  if (isnan(x) || isnan(y))
    quit("poly_add_point(): NaN detected\n");

  if (pl->n == pl->nallocated) {
    pl->x = realloc(pl->x, pl->nallocated * sizeof(double) * 2);
    pl->y = realloc(pl->y, pl->nallocated * sizeof(double) * 2);
    pl->nallocated *= 2;
  }

  pl->x[pl->n] = x;
  pl->y[pl->n] = y;
  pl->n++;
  extent_update(&pl->e, x, y);
}

/** Appends a point to the tail of a polyline.
 * @param pl Polyline.
 * @param n Number of points
 * @param x Array [n] of X coordinates
 * @param y Array [n] of Y coordinates
 */
void poly_add_points(poly_t *pl, int n, double x[], double y[])
{
  int nfin = pl->n + n;
  int i;

  if (nfin > pl->nallocated) {
    pl->x = realloc(pl->x, nfin * sizeof(double));
    pl->y = realloc(pl->y, nfin * sizeof(double));
    pl->nallocated = nfin;
  }

  for (i = 0; i < n; ++i) {
    double xx = x[i];
    double yy = y[i];

    if (isnan(xx) || isnan(yy))
      quit("poly_add_point(): NaN detected\n");

    pl->x[pl->n] = xx;
    pl->y[pl->n] = yy;
    pl->n++;

    extent_update(&pl->e, xx, yy);
  }

}

/** Adds a point to a polyline at a given position. No action for invalid 
 * index.
 * @param pl Polyline.
 * @param index Index.
 * @param x X coordinate.
 * @param y Y coordinate.
 */
void poly_add_point_at(poly_t *pl, int index, double x, double y)
{
  if (index > pl->n - 1) {
    poly_add_point(pl, x, y);
    return;
  }

  if (pl->n == pl->nallocated) {
    pl->x = realloc(pl->x, pl->nallocated * sizeof(double) * 2);
    pl->y = realloc(pl->y, pl->nallocated * sizeof(double) * 2);
    pl->nallocated *= 2;
  }

  memmove(&pl->x[index + 1], &pl->x[index],
          (pl->n - index) * sizeof(double));
  memmove(&pl->y[index + 1], &pl->x[index],
          (pl->n - index) * sizeof(double));
  pl->x[index] = x;
  pl->y[index] = y;
  pl->n++;
  extent_update(&pl->e, x, y);
}

/** Appends one polyline to another, like strcat().
 * @param pl1 Destination polyline.
 * @param pl2 Source polyline.
 */
void poly_append(poly_t *pl1, poly_t *pl2)
{
  int n = pl1->n + pl2->n;
  int sizechanged = 0;

  while (n < pl1->nallocated) {
    pl1->nallocated *= 2;
    sizechanged = 1;
  }

  if (sizechanged) {
    pl1->x = realloc(pl1->x, pl1->nallocated * sizeof(double));
    pl1->y = realloc(pl1->y, pl1->nallocated * sizeof(double));
  }

  memcpy(&pl1->x[pl1->n], pl2->x, pl2->n * sizeof(double));
  memcpy(&pl1->y[pl1->n], pl2->y, pl2->n * sizeof(double));
}

/** Computes the area of a polygon. For an open polyline, effectively closes
 * it. The area is positive for counterclockwise polygon, negative otherwise.
 * @param pl Polyline.
 */
double poly_area(poly_t *pl)
{
  double area = 0.0;
  double *x = pl->x;
  double *y = pl->y;
  int n = pl->n;
  int i;

  for (i = 0; i < n; ++i) {
    int i1 = (i + 1) % n;
    area += (x[i1] - x[i]) * (y[i1] + y[i]);
  }

  return area / 2.0;
}

/** Clears a polyline. Does not deallocate memory; use poly_destroy() for that.
 * @param pl Polyline.
 */
void poly_clear(poly_t *pl)
{
  pl->n = 0;
  extent_clear(&pl->e);
}

/** Makes a deep copy of a polyline.
 * @param pl Polyline.
 * @return Polyline.
 */
poly_t *poly_copy(poly_t *pl)
{
  poly_t *pl1 = malloc(sizeof(poly_t));

  pl1->n = pl->n;
  pl1->nallocated = pl->nallocated;
  pl1->e.xmin = pl->e.xmin;
  pl1->e.xmax = pl->e.xmax;
  pl1->e.ymin = pl->e.ymin;
  pl1->e.ymax = pl->e.ymax;

  pl1->x = malloc(pl1->nallocated * sizeof(double));
  pl1->y = malloc(pl1->nallocated * sizeof(double));
  memcpy(pl1->x, pl->x, pl->n * sizeof(double));
  memcpy(pl1->y, pl->y, pl->n * sizeof(double));

  return pl1;
}

/** Closes a polyline by adding the first point to the tail if necessary.
 * @param pl Polyline.
 * @return Polyline.
 */
void poly_close(poly_t *pl)
{
  if (!poly_is_closed(pl, 0.0))
    poly_add_point(pl, pl->x[0], pl->y[0]);
}

/** Tests whether a point is inside a polygon.
 * The polyline is assumed to be closed: an extra line segment from the end 
 * point back to the start point is assumed if necessary.
 * @param pl Polyline.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @return 1 for yes, 0 for no.
 */
int poly_contains_point(poly_t *pl, double x, double y)
{
  double *xs = pl->x;
  double *ys = pl->y;
  int n = pl->n;
  int hits;
  int i;

  if (n <= 1)
    return 0;
  if (!extent_contains_point(&pl->e, x, y))
    return 0;

  for (i = 0, hits = 0; i < n; ++i) {
    int i1 = (i + 1) % n;
    double x1 = xs[i] - x;
    double y1 = ys[i] - y;
    double x2 = xs[i1] - x;
    double y2 = ys[i1] - y;
    if (y1 == 0.0 && y2 == 0.0) {
      if (x1 * x2 <= 0)
        return 1;
    } else if (y1 == 0.0) {
      if (x1 == 0.0)
        return 1;
      if (x1 > 0.0)
        hits += (y2 > 0.0) ? 1 : -1;
    } else if (y2 == 0.0) {
      if (x2 == 0.0)
        return 1;
      if (x2 > 0.0)
        hits += (y1 < 0.0) ? 1 : -1;
    } else if (y1 * y2 < 0.0) {
      if (x1 > 0.0 && x2 > 0.0)
        hits += 2;
      else if (x1 * x2 <= 0.0) {
        double xx = x1 - (x2 - x1) * y1 / (y2 - y1);
        if (xx == 0)
          return 1;
        if (xx > 0.0)
          hits += 2;
      }
    }
  }

  if ((hits / 2) % 2)
    return 1;

  return 0;
}

/** Constructor.
 * @return Polyline.
 */
poly_t *poly_create(void)
{
  poly_t *pl = malloc(sizeof(poly_t));
  pl->x = malloc(PL_NSTART * sizeof(double));
  pl->y = malloc(PL_NSTART * sizeof(double));
  pl->n = 0;
  pl->nallocated = PL_NSTART;
  extent_clear(&pl->e);

  return pl;
}

/** Removes a point at a given position from a polyline.
 * @param pl Polyline.
 * @param index Position.
 */
void poly_delete_point(poly_t *pl, int index)
{
  double xx = pl->x[index];
  double yy = pl->y[index];
  extent_t *e = &pl->e;

  memmove(&pl->x[index], &pl->x[index + 1],
          (pl->n - index - 1) * sizeof(double));
  memmove(&pl->y[index], &pl->y[index + 1],
          (pl->n - index - 1) * sizeof(double));
  pl->n--;

  if (e->xmin == xx || e->xmax == xx || e->ymin == yy || e->ymax == yy) {
    double *x = pl->x;
    double *y = pl->y;
    int i;
    extent_clear(e);
    for (i = 0; i < pl->n; ++i)
      extent_update(e, x[i], y[i]);
  }
}

/** Destructor.
 * @param pl Polyline.
 */
void poly_destroy(poly_t *pl)
{
  if (pl == NULL)
    return;
 
  if (pl->x != NULL )
    free(pl->x);
  if (pl->y != NULL )
    free(pl->y);
 
  free(pl);
  pl = NULL;
  
}

/** Finds index of a point within polyline.
 * @param pl Polyline.
 * @param x X coordinate.
 * @param y Y coordinate.
 * @return Index if found; -1 otherwise.
 */
int poly_find_index(poly_t *pl, double x, double y)
{
  double *xs = pl->x;
  double *ys = pl->y;
  int n = pl->n;
  int i;

  if(isnan(x) || isnan(y))
    return -1;

  for (i = 0; i < n; ++i) {
    if (x == xs[i] && y == ys[i])
      return i;
  }

  return -1;
}

/** Checks whether the polyline is closed.
 * @param pl Polyline.
 * @param eps Distance tolerance.
 * @return 1 for yes, 0 for no.
 */
int poly_is_closed(poly_t *pl, double eps)
{
  return pl->n > 1 &&
    point_equal(pl->x[0], pl->y[0], pl->x[pl->n - 1], pl->y[pl->n - 1],
                eps);
}

/** Reads points from a file stream and appends them to a polyline.
 * @param pl Polyline.
 * @param fp File handle.
 * @return The number of points in the polyline.
 */
int poly_read(poly_t *pl, FILE * fp)
{
  double x, y;
  char buf[MAXLINELEN];

  /* Skip comments and blank lines */
  if (prm_next_line(buf, MAXLINELEN, fp) == 0)
    return 0;

  while (sscanf(buf, "%lf %lf", &x, &y) == 2) {
    poly_add_point(pl, x, y);
    if (fgets(buf, MAXLINELEN, fp) == NULL)
      break;
  }

  return (pl->n);
}

/** Resamples a polyline including only points more than the threshold
 * distance apart.
 * @param pl Polyline.
 * @param eps Threshold distance.
 */
void poly_resample(poly_t *pl, double eps)
{
  double *xs = pl->x;
  double *ys = pl->y;
  int n = pl->n;
  int i, now;

  if (pl->n <= 1)
    return;

  for (i = 1, now = 0; i < n; ++i) {
    if (!point_equal(xs[now], ys[now], xs[i], ys[i], eps)) {
      now++;
      xs[now] = xs[i];
      ys[now] = ys[i];
    }
  }

  if (pl->n != now + 1) {
    pl->n = now + 1;
    recalc_extent(pl);
  }
}

/** Reverse points in a polyline.
 * @param pl Polyline.
 */
void poly_reverse(poly_t *pl)
{
  double *x = pl->x;
  double *y = pl->y;
  int n = pl->n;
  int n2 = n / 2;
  int i, j;

  if (pl->n <= 1)
    return;

  for (i = 0, j = pl->n - 1; i < n2; ++i, --j) {
    double tmp = x[i];
    x[i] = x[j];
    x[j] = tmp;
    tmp = y[i];
    y[i] = y[j];
    y[j] = tmp;
  }
}

/** Smoothes a polyline. Fits a spline to each of x and y as functions of
 * distance along the line.
 * @param pl Polyline.
 * @param ns Number of points in smoothed polyline.
 * @return Smoothed polyline.
 */
poly_t *poly_smooth(poly_t *pl, int ns)
{
  int ns1 = ns - 1;
  int n1 = pl->n - 1;
  double *x = pl->x;
  double *y = pl->y;

  poly_t *pl1 = poly_create();
  double *x1;
  double *y1;

  double *d;                    /* distance from start values */
  double *ydd;                  /* spline computed second derivatives */
  int i;

  if (pl->n < 3)
    return pl1;

  if (ns < 2)
    quit("poly_smooth(): ns < 2\n");

  x1 = calloc(ns, sizeof(double));
  y1 = calloc(ns, sizeof(double));
  d = calloc(pl->n, sizeof(double));
  ydd = calloc(pl->n, sizeof(double));

  d[0] = 0.0;
  for (i = 1; i < pl->n; ++i)
    d[i] = d[i - 1] + hypot(x[i] - x[i - 1], y[i] - y[i - 1]);

  /* fit spline to x */
  spline(d, x, pl->n, 0, 1.0e+31, 1.0e+31, ydd);

  /* loop to evaluate interpolated x values */
  for (i = 0; i < ns; i++) {
    double dist = (i < ns1) ? i * d[n1] / ns1 : d[n1];
    spline_interp(d, x, ydd, pl->n, dist, &x1[i]);
  }

  /* fit spline to y */
  spline(d, y, pl->n, 0, 1.0e+31, 1.0e+31, ydd);

  /* loop to evaluate interpolated y values */
  for (i = 0; i < ns; i++) {
    double dist = (i < ns1) ? i * d[n1] / ns1 : d[n1];
    spline_interp(d, y, ydd, pl->n, dist, &y1[i]);
  }

  for (i = 0; i < ns; i++)
    poly_add_point(pl1, x1[i], y1[i]);

  free(x1);
  free(y1);
  free(d);
  free(ydd);

  return pl1;
}

/** Remove spikes from a polyline.
 * @param pl Polyline.
 * @param maxdist Maximal allowed distance between two adjacent points.
 */
void poly_despike(poly_t *pl, double maxdist)
{
  double *xs = pl->x;
  double *ys = pl->y;
  int n = pl->n;
  int i, now;

  if (pl->n <= 1)
    return;

  for (i = 1, now = 0; i < n; ++i) {
    if (hypot(xs[now] - xs[i], ys[now] - ys[i]) < maxdist) {
      now++;
      xs[now] = xs[i];
      ys[now] = ys[i];
    }
  }

  if (pl->n != now + 1) {
    pl->n = now + 1;
    recalc_extent(pl);
  }
}

/** Writes polyline to a file stream.
 * @param pl Polyline.
 * @param fp File handle.
 */
void poly_write(poly_t *pl, FILE * fp)
{
  int i;

  fprintf(fp, "# %d\n", pl->n);
  for (i = 0; i < pl->n; ++i)
    fprintf(fp, "%.4f %.4f\n", pl->x[i], pl->y[i]);
}

/* Forms polyline boundary around a grid.
 * Note: supposed to handle a grid of corner nodes only.
 * @param nce1 Number of cells in X direction
 * @param nce2 Number of cells in Y direction
 * @param x X coordinates of grid nodes
 * @param y Y coordinates of grid nodes
 * @return Boundary polyline
 */
poly_t* poly_formbound(int nce1, int nce2, double** x, double** y)
{
    poly_t* pl = poly_create();
    int direction = 0;          /* 0 - down, 1 - right, 2 - up, 3 -left */
    int iinc[] = { 0, 1, 0, -1 };
    int jinc[] = { 1, 0, -1, 0 };
    int i, j, istart, jstart;

    for (j = 0; j <= nce2; ++j)
        for (i = 0; i <= nce1; ++i)
            if (!isnan(x[j][i]))
                goto ok;

  ok:
    if (j > nce2)
        return pl;

    istart = i;
    jstart = j;
    poly_add_point(pl, x[j][i], y[j][i]);

    do {
        int direction_stop = (direction + 2) % 4;
        int inext, jnext;

        direction = (direction + 3) % 4;
        inext = i + iinc[direction];
        jnext = j + jinc[direction];

        while ((inext < 0 || jnext < 0 || inext > nce1 || jnext > nce2 || isnan(x[jnext][inext])) && direction != direction_stop) {
            direction = (direction + 1) % 4;
            inext = i + iinc[direction];
            jnext = j + jinc[direction];
        }

        if (direction == direction_stop)
            break;

        i = inext;
        j = jnext;

        poly_add_point(pl, x[j][i], y[j][i]);

    } while (j != jstart || i != istart);

    poly_close(pl);

    return (pl);
}

/* Forms polyline boundary around a grid in index space.
 * @param nce1 Number of cells in X direction
 * @param nce2 Number of cells in Y direction
 * @param x X coordinates of grid nodes
 * @return Boundary polyline
 */
poly_t* poly_formboundij(int nce1, int nce2, double** x)
{
    poly_t* pl = poly_create();
    int direction = 0;          /* 0 - down, 1 - right, 2 - up, 3 -left */
    int iinc[] = { 0, 1, 0, -1 };
    int jinc[] = { 1, 0, -1, 0 };
    int i, j, istart, jstart;

    for (j = 0; j <= nce2; ++j)
        for (i = 0; i <= nce1; ++i)
            if (!isnan(x[j][i]))
                goto ok;

  ok:
    if (j > nce2)
        return pl;

    istart = i;
    jstart = j;
    poly_add_point(pl, i, j);

    do {
        int direction_stop = (direction + 2) % 4;
        int inext, jnext;

        direction = (direction + 3) % 4;
        inext = i + iinc[direction];
        jnext = j + jinc[direction];

        while ((inext < 0 || jnext < 0 || inext > nce1 || jnext > nce2 || isnan(x[jnext][inext])) && direction != direction_stop) {
            direction = (direction + 1) % 4;
            inext = i + iinc[direction];
            jnext = j + jinc[direction];
        }

        if (direction == direction_stop)
            break;

        i = inext;
        j = jnext;

        poly_add_point(pl, i, j);

    } while (j != jstart || i != istart);

    poly_close(pl);

    return (pl);
}

/* Deletes redundant nodes.
 * @param pl Polyline
 * @param eps A small number used in tests on two points being the same or
 *            three points belonging to one line.
 */
void poly_compact(poly_t* pl, double eps)
{
    int n = pl->n;
    double* x = pl->x;
    double* y = pl->y;
    int* ids = NULL;
    int nnew;
    int ileft, imiddle, iright;
    int i;

    if (n <= 4)
        return;                 /* do not bother */

    ids = malloc(pl->n * sizeof(int));
    nnew = 0;

    for (i = 0, imiddle = 0, ileft = n - 1; i < n - 1; ++i)
        if (!point_equal(x[ileft], y[ileft], x[imiddle], y[imiddle], eps))
            break;
        else
            ileft--;

    for (i = 0, iright = 1; i < n; ++i) {
        if (!point_equal(x[iright], y[iright], x[imiddle], y[imiddle], eps)) {
            if (!ononeline(x, y, ileft, imiddle, iright, eps)) {
                ids[nnew++] = imiddle;
                ileft = imiddle;
            }
            imiddle = iright;
        }
        iright = (iright + 1) % n;
    }

    if (nnew != n) {
        for (i = 0; i < nnew; ++i) {
            x[i] = x[ids[i]];
            y[i] = y[ids[i]];
        }
        pl->n = nnew;
    }

    free(ids);
    pl->x = realloc(pl->x, sizeof(double) * pl->n);
    pl->y = realloc(pl->y, sizeof(double) * pl->n);
    pl->nallocated = pl->n;
}

// EOF
