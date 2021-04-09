/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/lsq_lin/lsql.c
 *  
 *  Description: 2D least squares linear interpolation
 *
 *  `lsql' -- "Least squares linear" -- is a structure for
 *  conducting least squares linear interpolation on a given data
 *  on a "point-to-point" basis.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lsql.c 6597 2020-09-03 05:29:26Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lsql.h"
#include "emsalloc.h"


struct lsql {
    delaunay* d;
    lweights* weights;
};

//
int lsql_verbose = 0;

/*-------------------------------------------------------------------*/
/* Gets the metrics of cells surrounding edge e using least squares  */
/* fitting, where the solution to npe equations of a polynomial of   */
/* degree nol is obtained via singular value decomposition.          */
/* The least squares fit uses cell centred values of cells           */
/* surrounding a given edge, e.                                      */
/* The polynomial used is:                                           */
/* t = co + c1.x + c2.y                                              */
/* Used (lat,lon) based coordinates.                                 */
/* @param d Delaunay triangulation                                   */
/* @return linear interpolator                                       */
/*-------------------------------------------------------------------*/
lsql* lsql_build(delaunay* d)
{
    int i, j;
    int npem;
    double **p, **b;
    double *f, *s, *w, *std;
    double z0, **coeff;
    lsql* l = malloc(sizeof(lsql));
    short *mask;
    int *npe;
    int trif = 0;
    int edef = 0;

    l->d = d;
    l->weights = malloc(d->npoints * sizeof(lweights));
    npe = i_alloc_1d(d->npoints);
    mask = s_alloc_1d(d->npoints);

    if (trif) {
      for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];
	/*if(t->vids[0]==784||t->vids[1]==784||t->vids[2]==784) {*/
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("%f %f\n",d->points[t->vids[1]].x, d->points[t->vids[1]].y);
	printf("%f %f\n",d->points[t->vids[2]].x, d->points[t->vids[2]].y);
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("NaN NaN\n");

      }
    }
    if (edef) {
      j = 0;
      for (i = 0; i < d->nedges; ++i) {
	printf("%f %f\n",d->points[d->edges[j]].x, d->points[d->edges[j]].y);
	printf("%f %f\n",d->points[d->edges[j+1]].x, d->points[d->edges[j+1]].y);
	printf("NaN NaN\n");
	j += 2;
      }
    }

    /*---------------------------------------------------------------*/
    /* Count the number of points associated with trianges having  i */
    /* as a vertex, and get the maximum number of points.            */
    npem = 0;
    memset(npe, 0, d->npoints * sizeof(int));
    /* If points used in the least squares have been set outside     */
    /* this routine (stored in d->point_trianges, and set d->ptf=1)  */
    /* then use these.                                               */
    if (d->ptf) {
      for (i = 0; i < d->npoints; ++i) {
	npe[i] = d->n_point_triangles[i];
	npem = max(npem, npe[i]);
      }
    } else {
      for (i = 0; i < d->npoints; ++i) {
	int nt = d->n_point_triangles[i];
	memset(mask, 0, d->npoints * sizeof(short));
	for (j = 0; j < nt; ++j) {
	  int n, ni, tid = d->point_triangles[i][j];
	  triangle* t = &d->triangles[tid];
	  for (n = 0; n < 3 ; n++) {
	    ni = t->vids[n];
	    if (!mask[ni]) {
	    npe[i]++;
	    mask[ni] = 1;
	    }
	  }
	}
	npem = max(npem, npe[i]);
	if(npe[i]==0)printf("zero points at id=%d [%f %f]: #tri=%d\n",i,d->points[i].x,d->points[i].y, nt);
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate                                                      */
    p = d_alloc_2d(nol, npem);
    b = d_alloc_2d(nol, npem);
    f = d_alloc_1d(nol);
    w = d_alloc_1d(nol);
    s = d_alloc_1d(npem);
    std = d_alloc_1d(npem);

    /*---------------------------------------------------------------*/
    /* Loop through all the points                                   */
    for (i = 0; i < d->npoints; ++i) {
      int nt = d->n_point_triangles[i];
      int jj, n;
      double *x, *y;

      lweights* lw = &l->weights[i];
      lw->ncells = npe[i];
      lw->orf = 3;
      lw->xr = d->points[i].x;
      lw->yr = d->points[i].y;
      lw->beta = 1.0;
      if (nt) {
	lw->B = d_alloc_2d(nol, npe[i]);
	lw->cells = i_alloc_1d(npe[i]);
	x = d_alloc_1d(npe[i]);
	y = d_alloc_1d(npe[i]);
      } else
	continue;

      /* Get the locations of surrounding points                     */
      jj = 0;
      memset(mask, 0, d->npoints * sizeof(short));
      for (j = 0; j < nt; ++j) {
	if (d->ptf) {
	  int ni = d->point_triangles[i][j];
	  if (ni >=0 && ni < d->npoints) {
	    x[jj] = d->points[ni].x - lw->xr;
	    y[jj] = d->points[ni].y - lw->yr;
	    lw->cells[jj++] = ni;
	  }
	} else {
	  int ni, tid = d->point_triangles[i][j];
	  triangle* t = &d->triangles[tid];
	  for (n = 0; n < 3 ; n++) {
	    ni = t->vids[n];
	    if (!mask[ni]) {
	      x[jj] = d->points[ni].x - lw->xr;
	      y[jj] = d->points[ni].y - lw->yr;
	      lw->cells[jj++] = ni;
	      mask[ni] = 1;
	    }
	  }
	}
      }

      /* Set up the weights matrix */
      for (j = 0; j < npe[i]; j++) {
	p[j][0] = 1.0;
	p[j][1] = x[j];
	p[j][2] = y[j];

	s[j] = 1.0;
	std[j] = 0.0;
      }

      /* Least squares fitting via singular value decomposition,     */
      /* where p = U * W * V^T, then f = B *s where b = V*W^-1*U^T.  */
      /* Note: these weights are computed using a scalar field (s)   */
      /* of 1.                                                       */
      svd_lsq_B(p, lw->orf, npe[i], s, NULL, b, f);

      /* Save the coefficients                                       */
      for (n = 0; n < lw->orf; n++) {
	for (j = 0; j < npe[i]; j++) {
	  lw->B[j][n] = b[j][n];
	}
      }

      /* Compute the coefficients                                    */
      memset(lw->w, 0, nol * sizeof(double));
      for (n = 0; n < lw->ncells; n++) { 
	j = lw->cells[n];
	z0 = d->points[j].z;
	for (jj = 0; jj < lw->orf; jj++) {
	  lw->w[jj] += z0 * lw->B[n][jj];
	}
      }
      d_free_1d(x);
      d_free_1d(y);
    }
    d_free_2d(b);
    d_free_2d(p);
    d_free_1d(f);
    d_free_1d(w);
    d_free_1d(s);
    d_free_1d(std);
    i_free_1d(npe);
    s_free_1d(mask);

    return l;
}

/*-------------------------------------------------------------------*/
/* Limit the weights to make a monotone solution                     */
/*-------------------------------------------------------------------*/
void lsql_limit(lsql* l, int npoint, point* p)
{
  int i, id;
  double z0, xp, yp, d;
  double beta = 1.0;

  /* Recalculate the weights using data in p.z */
  id = (int)p[0].z;
  lweights* lw = &l->weights[id];

  lw->beta = 1.0;

  for (i = 0; i < npoint; ++i) {
    lsql_interpolate_point(l, &p[i]);
    xp = p[i].x - lw->xr;
    yp = p[i].y - lw->yr;
    d = lw->w[1] * xp + lw->w[2] * yp;
    z0 = p[i].z;
    if (d != 0.0 && z0 > p[i].v[1]) {
      beta = max(min((p[i].v[1] - lw->w[0]) / d, beta), 0.0);
    }
    if (d != 0.0 && z0 < p[i].v[0]) {
      beta = max(min((p[i].v[0] - lw->w[0]) / d, beta), 0.0);
    }
  }
  lw->beta = min(max(beta, 0.0), 1.0);

}

/*-------------------------------------------------------------------*/
/* Remakes the weights from cached weights                           */
/*-------------------------------------------------------------------*/
void lsql_rebuild(lsql* l, point* p)
{
  int i, j, jj, n;
  double z0;
  delaunay *d = l->d;

  /* Recalculate the weights using data in p.z */
  for (i = 0; i < d->npoints; ++i) {
    lweights* lw = &l->weights[i];
    memset(lw->w, 0, nol * sizeof(double));
    lw->tmx = -1e10;
    lw->tmn = 1e10;
    for (n = 0; n < lw->ncells; n++) { 
      j = lw->cells[n];
      z0 = p[j].z;
      lw->tmx = max(lw->tmx, z0);
      lw->tmn = min(lw->tmn, z0);
      for (jj = 0; jj < lw->orf; jj++) {
	lw->w[jj] += z0 * lw->B[n][jj];
      }
    }
    /* Set the leading term to the cell mean so that the integral of */
    /* the lsq function over a cell equalls the cell mean value.     */
    lw->w[0] = p[i].z;
  }
}


/* Destroys quadratic interpolator.
 *
 * @param l Structure to be destroyed
 */
void lsql_destroy(lsql* l)
{
    d_free_2d(l->weights->B);
    i_free_1d(l->weights->cells);
    free(l->weights);
    free(l);
}

/* Finds lsq quadratic interpolated value in a point.
 *
 * @param l Quadratic interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void lsql_interpolate_point(lsql* l, point* p)
{
  int i;
  delaunay* d = l->d;
  lweights* lw;
  double v;
  int id = (int)p->z;
  double xp, yp;
  int jj;
  lw = &l->weights[id];
  xp = p->x - lw->xr;
  yp = p->y - lw->yr;

  v = lw->w[0] + lw->beta * (lw->w[1] * xp + lw->w[2] * yp);

  if (lw->orf == 6) {
    v += (lw->w[3] * xp * xp +
	  lw->w[4] * xp * yp +
	  lw->w[5] * yp * yp);
  }

  p->z = v;
}

void lsql_minmax(lsql* l, int id, double *mn, double *mx)
{
  lweights* lw = &l->weights[id];
  *mn = lw->tmn;
  *mx = lw->tmx;
}

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void lsql_interpolate_points(int nin, point pin[], int nout, point pout[])
{
    delaunay* d = delaunay_build(nin, pin, 0, NULL, 0, NULL);
    lsql* l = lsql_build(d);
    int seed = 0;
    int i;

    if (lsql_verbose) {
        fprintf(stderr, "xytoi:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];

            fprintf(stderr, "(%.7g,%.7g) -> %d\n", p->x, p->y, delaunay_xytoi(d, p, seed));
        }
    }

    for (i = 0; i < nout; ++i)
        lsql_interpolate_point(l, &pout[i]);

    if (lsql_verbose) {
        fprintf(stderr, "output:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];;
            fprintf(stderr, "  %d:%15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    lsql_destroy(l);
    delaunay_destroy(d);
}

// EOF
