/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/linear/baycentric.c
 *  
 *  Description: 2D baycentric linear interpolation
 *
 *  `baycentric' -- "Baycentric interpolation" -- is a structure for
 *  conducting linear baycentric interpolation on a given data on a
 *  "point-to-point" basis. This should be used with an underlying
 *  Voronoi mesh generated from COMPAS, where the Delaunay trianges
 *  are the dual. We find which triangle a point (x,y) falls within,
 *  then use a linear Baycentric interpolation using the three cell
 *  centres comprising that triangle.
 *
 *  Reference:
 *  https://codeplea.com/triangular-interpolation
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: baycentric.c 5859 2018-07-02 04:04:44Z her127 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "baycentric.h"

#define nop 3 /* Number of polynomial coefficients                   */

typedef struct {
  double w[nop];  /* Weights                                         */
  int ncells;     /* Number of cells surrounding a point             */
  int *cells;     /* Central cell and surrounding cell indices       */
  double *x;      /* Central cell and surrounding cell x location    */
  double *y;      /* Central cell and surrounding cell y location    */
  double tmx;     /* Maximum tracer value                            */
  double tmn;     /* Minimum tracer value                            */
  double xr, yr;  /* Central cell (x,y) location                     */
  int **tri;      /* Indices of the triangles                        */
  int idt;        /* Triangle that point (x,y) falls into            */
  int isghost;    /* = 1 for ghost cells                             */
  double origin;  /* Origin angle to find triangle to use            */
  double *angle;
} lweights;

struct bal {
    delaunay* d;
    lweights* weights;
};

//
int bal_verbose = 0;

/*-------------------------------------------------------------------*/
/* Gets the locations of cells surrounding cell c and the dimensions */
/* of cell c (dx and dy).                                            */
/* Used (lat,lon) based coordinates.                                 */
/* @param d Delaunay triangulation                                   */
/* @return baycentric interpolator                                   */
/*-------------------------------------------------------------------*/
bal* bal_build(delaunay* d)
  {
    int i, j, jp;
    int npem;
    bal* l = malloc(sizeof(bal));

    l->weights = malloc(d->npoints * sizeof(lweights));
    l->d = d;

    /*---------------------------------------------------------------*/
    /* Count the number of points associated with trianges having i  */
    /* as a vertex, and get the maximum number of points.            */
    npem = 0;
    for (i = 0; i < d->npoints; ++i) {
      npem = max(npem, d->n_point_triangles[i]);
    }

    /*---------------------------------------------------------------*/
    /* Loop through all the points                                   */
    for (i = 0; i < d->npoints; ++i) {
      int ncells = d->n_point_triangles[i];

      double x, y;

      int ntri = ncells - 1;

      lweights *lw = &l->weights[i];

      lw->ncells = ncells;
      lw->isghost = 0;
      if (ncells == 0) {
	lw->isghost = 1;
	continue;
      }
      if (d->point_triangles[i][0] == -1) {
	lw->isghost = 1;
	continue;
      }

      lw->cells = i_alloc_1d(ncells);
      lw->angle = d_alloc_1d(ncells);
      lw->x = d_alloc_1d(ncells);
      lw->y = d_alloc_1d(ncells);
      lw->xr = d->points[i].x;
      lw->yr = d->points[i].y;

      /* Get the locations of surrounding points                     */
      for (j = 0; j < ncells; ++j) {
	lw->cells[j] = d->point_triangles[i][j];
      }

      /* Initialize the coefficients                                 */
      memset(lw->w, nop, nop * sizeof(double));

      lw->tri = i_alloc_2d(3, ntri);

      /* Get the indices of the dual and coordinates of the          */
      /* circumcentre.                                               */
      lw->x[0] = d->points[lw->cells[0]].x;
      lw->y[0] = d->points[lw->cells[0]].y;
      for (j = 1; j < ncells; j++) {
	lw->x[j] = d->points[lw->cells[j]].x;
	lw->y[j] = d->points[lw->cells[j]].y;
	jp = (j == ncells-1) ? 1 : j + 1;
	lw->tri[j-1][0] = 0;
	lw->tri[j-1][1] = j;
	lw->tri[j-1][2] = jp;
	x = d->points[lw->cells[j]].x - lw->xr;
	y = d->points[lw->cells[j]].y - lw->yr;
	lw->angle[j] = atan2(y, x);
      }
      x = d->points[lw->cells[1]].x - lw->xr;
      y = d->points[lw->cells[1]].y - lw->yr;
      lw->origin = PI - atan2(y, x);
    }
    return l;
}

/*-------------------------------------------------------------------*/
/* Remakes the weights from cached weights                           */
/*-------------------------------------------------------------------*/
void bal_rebuild(bal* l, point* p)
{
  int i, j;
  int id, idt, ntri;
  double xp, yp;
  double d1;
  double x1, y1, x2, y2, x3, y3;
  delaunay *d = l->d;
  lweights *lw, *lws;

  for (i = 0; i < d->npoints; ++i) {
    id = (int)p[i].z;
    lw = &l->weights[i];
    lws = &l->weights[id];

    if(lw->isghost) continue;
    if(lws->isghost) printf("Ghost source = %d from %d\n",id,i);

    /* Get the location of the point to interpolate                  */
    xp = p[i].x;
    yp = p[i].y;
    ntri = lws->ncells - 1;

    /* Get the triangle the location (xp, yp) falls into. Use the    */
    /* angle relative to the central point, and find the triangle    */
    /* which brackets this angle.                                    */
    x1 = xp - lws->xr;
    y1 = yp - lws->yr;

    x1 = atan2(y1, x1);
    idt = 0;

    /*
    inc = 2.0 * PI / (double)ntri;
    d1 = lws->origin;
    for (j = 0; j < ntri; j++) {
      if (x1 >= d1 && x1 < d1 + inc) {
	idt = j;
	break;
      }
      d1 += inc;
    }
    */
    for (j = 1; j <= ntri; j++) {
      int jp =  j + 1;
      double sc = 1.0;
      if (j == ntri) {
	jp = 1;
	sc = -1.0;
      }
      if (x1 <= lws->angle[j] && x1 > sc * lws->angle[jp]) {
	idt = j - 1;
	break;
      }
    }

    /* Get the Baycentric weights                                    */
    lw->idt = idt;
    x1 = lws->x[lws->tri[idt][0]];
    y1 = lws->y[lws->tri[idt][0]];
    x2 = lws->x[lws->tri[idt][1]];
    y2 = lws->y[lws->tri[idt][1]];
    x3 = lws->x[lws->tri[idt][2]];
    y3 = lws->y[lws->tri[idt][2]];

    d1 = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3) ;

    lw->w[0] = ((y2 - y3) *(xp - x3) + (x3 - x2) * (yp - y3)) / d1;
    lw->w[1] = ((y3 - y1) *(xp - x3) + (x1 - x3) * (yp - y3)) / d1;
    lw->w[2] = 1.0 - lw->w[0] - lw->w[1];
  }
}

void bal_rebuild2(bal* l1, bal* l2, point* p)
{
  int i, j;
  int id, idt, ntri;
  double xp, yp;
  double d1;
  double x1, y1, x2, y2, x3, y3;
  delaunay *d = l1->d;
  lweights *lw, *lws;
  int vid = d->vid;

  i = (int)p->v[vid];
  id = (int)p->z;

  lw = &l1->weights[i];
  lws = &l2->weights[id];

  if(lw->isghost) return;

  /* Get the location of the point to interpolate                    */
  xp = p->x;
  yp = p->y;
  ntri = lws->ncells - 1;

  /* Get the triangle the location (xp, yp) falls into. Use the      */
  /* angle relative to the central point, and find the triangle      */
  /* which brackets this angle.                                      */
  x1 = xp - lws->xr;
  y1 = yp - lws->yr;

  x1 = atan2(y1, x1);
  idt = 0;

  for (j = 1; j <= ntri; j++) {
    int jp =  j + 1;
    double sc = 1.0;
    if (j == ntri) {
      jp = 1;
      sc = -1.0;
    }
    if (x1 <= lws->angle[j] && x1 > sc * lws->angle[jp]) {
      idt = j - 1;
      break;
    }
  }

  /* Get the Baycentric weights                                      */
  lw->idt = idt;
  x1 = lws->x[lws->tri[idt][0]];
  y1 = lws->y[lws->tri[idt][0]];
  x2 = lws->x[lws->tri[idt][1]];
  y2 = lws->y[lws->tri[idt][1]];
  x3 = lws->x[lws->tri[idt][2]];
  y3 = lws->y[lws->tri[idt][2]];
  
  d1 = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3) ;
  
  lw->w[0] = ((y2 - y3) *(xp - x3) + (x3 - x2) * (yp - y3)) / d1;
  lw->w[1] = ((y3 - y1) *(xp - x3) + (x1 - x3) * (yp - y3)) / d1;
  lw->w[2] = 1.0 - lw->w[0] - lw->w[1];
}


/* Destroys baycentric interpolator.
 *
 * @param l Structure to be destroyed
 */
void bal_destroy(bal* l)
{
    i_free_1d(l->weights->cells);
    i_free_2d(l->weights->tri);
    d_free_1d(l->weights->x);
    d_free_1d(l->weights->y);
    free(l->weights);
    free(l);
}

/* Finds linear baycentric interpolated value in a point.
 *
 * @param l Baycentric interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void bal_interpolate_point(bal* l, point* p)
{
  int i, j;
  delaunay* d = l->d;
  lweights *lw, *lws;
  double v, z0;
  int id = (int)p->x;
  int ids = (int)p->y;
  int vid = d->vid;
  int idt;

  lw = &l->weights[id];
  lws = &l->weights[ids];

  v = 0.0;
  lw->tmx = -1e10;
  lw->tmn = 1e10;
  for (i = 0; i < 3; i++) {
    idt = lw->idt;
    j = lws->cells[lws->tri[idt][i]];
    z0 = (j == -1) ?  d->points[ids].v[vid] : d->points[j].v[vid];
    v += lw->w[i] * z0;
    /*if(id==623&&vid==3)printf("interp %d %d i=%d j=%d w=%f z0=%f %f\n",ids,idt,i,j,lw->w[i],z0,v);*/
    lw->tmx = max(lw->tmx, z0);
    lw->tmn = min(lw->tmn, z0);
  }
  v = min(lw->tmx, max(lw->tmn, v));
  p->z = v;
}


void bal_interpolate_point2(bal* l1, bal* l2, point* p)
{
  int i, j;
  delaunay* d = l2->d;
  lweights *lw, *lws;
  double v, z0;
  int id = (int)p->x;
  int ids = (int)p->y;
  int vid = d->vid;
  int idt;

  lw = &l1->weights[id];
  lws = &l2->weights[ids];

  v = 0.0;
  lw->tmx = -1e10;
  lw->tmn = 1e10;
  for (i = 0; i < 3; i++) {
    idt = lw->idt;
    j = lws->cells[lws->tri[idt][i]];
    z0 = (j == -1) ?  d->points[ids].v[vid] : d->points[j].v[vid];
    v += lw->w[i] * z0;
    /*if(id==623&&vid==3)printf("interp %d %d i=%d j=%d w=%f z0=%f %f\n",ids,idt,i,j,lw->w[i],z0,v);*/
    lw->tmx = max(lw->tmx, z0);
    lw->tmn = min(lw->tmn, z0);
  }
  v = min(lw->tmx, max(lw->tmn, v));
  p->z = v;
}

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void bal_interpolate_points(int nin, point pin[], int nout, point pout[])
{
    delaunay* d = delaunay_build(nin, pin, 0, NULL, 0, NULL);
    bal* l = bal_build(d);
    int seed = 0;
    int i;

    if (bal_verbose) {
        fprintf(stderr, "xytoi:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];

            fprintf(stderr, "(%.7g,%.7g) -> %d\n", p->x, p->y, delaunay_xytoi(d, p, seed));
        }
    }

    for (i = 0; i < nout; ++i)
        bal_interpolate_point(l, &pout[i]);

    if (bal_verbose) {
        fprintf(stderr, "output:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];;
            fprintf(stderr, "  %d:%15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    bal_destroy(l);
    delaunay_destroy(d);
}

// EOF
