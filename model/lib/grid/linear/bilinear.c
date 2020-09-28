/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/linear/bilinear.c
 *  
 *  Description: 2D bilinear interpolation
 *
 *  `bilinear' -- "Bilinear interpolation" -- is
 *                 a structure for conducting bilinear interpolation on a given
 *                 data on a "point-to-point" basis.
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: bilinear.c 5859 2018-07-02 04:04:44Z her127 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bilinear.h"

#define nop 4 /* Number of polynomial coefficients             */

int qpoints[4][4] = {
  {0, 1, 5, 2},
  {0, 3, 6, 2},
  {0, 3, 7, 4},
  {0, 1, 8, 4}
};

typedef struct {
    double w[nop];
    int ncells;
    int *cells;
    double tmx;
    double tmn;
    double xr;
    double yr;
    double dx[nop];
    double dy[nop];
    int qf, *qa;
    int isghost;
} lweights;

struct bl {
    delaunay* d;
    lweights* weights;
};

//
int bl_verbose = 0;

/*-------------------------------------------------------------------*/
/* Gets the locations of cells surrounding cell c and the dimensions */
/* of cell c (dx and dy).                                            */
/* Used (lat,lon) based coordinates.                                 */
/* @param d Delaunay triangulation                                   */
/* @return bilinear interpolator                                     */
/*-------------------------------------------------------------------*/
bl* bl_build(delaunay* d)
{
    int i, j;
    int npem;
    bl* l = malloc(sizeof(bl));

    l->weights = malloc(d->npoints * sizeof(lweights));
    l->d = d;

    /*---------------------------------------------------------------*/
    /* Count the number of points associated with trianges having  i */
    /* as a vertex, and get the maximum number of points.            */
    npem = 0;
    for (i = 0; i < d->npoints; ++i) {
      npem = max(npem, d->n_point_triangles[i]);
    }

    /*---------------------------------------------------------------*/
    /* Loop through all the points                                   */
    for (i = 0; i < d->npoints; ++i) {
      int nt = d->n_point_triangles[i];

      lweights *lw = &l->weights[i];
      lw->ncells = nt;

      lw->isghost = 0;
      if (nt == 0) {
	lw->isghost = 1;
	continue;
      }
      if (d->point_triangles[i][0] == -1) {
	lw->isghost = 1;
	continue;
      }

      lw->cells = i_alloc_1d(nt);
      lw->xr = d->points[i].x;
      lw->yr = d->points[i].y;

      /* Get the locations of surrounding points                     */
      for (j = 0; j < nt; ++j) {
	lw->cells[j] = d->point_triangles[i][j];
      }

      /* Initialize the coefficients                                 */
      memset(lw->w, nop, nop * sizeof(double));

      /* Get the quadrant edge distances                             */
      /* Quadrant 0 and 3 x distance                                 */
      j = d->point_triangles[i][1];
      lw->dx[0] = lw->dx[3] = fabs(d->points[i].x - d->points[j].x);
      /* Quadrant 0 and 1 y distance                                 */
      j = d->point_triangles[i][2];
      lw->dy[0] = lw->dy[1] = fabs(d->points[i].y - d->points[j].y);
      /* Quadrant 1 and 2 x distance                                 */
      j = d->point_triangles[i][3];
      lw->dx[1] = lw->dx[2] = fabs(d->points[i].x - d->points[j].x);
      /* Quadrant 2 and 3 y distance                                 */
      j = d->point_triangles[i][4];
      lw->dy[2] = lw->dy[3] = fabs(d->points[i].y - d->points[j].y);
    }
    return l;
}


/* Remakes the weights from cached weights. The indexing of the      */
/* stencil in the Delaunay point_triangles list is as follows:       */
/* index 0 = central cell                                            */
/* index 1 = west cell                                               */
/* index 2 = north cell                                              */
/* index 3 = east cell                                               */
/* index 4 = south cell                                              */
/* index 5 = NW cell                                                 */
/* index 6 = NE cell                                                 */
/* index 7 = SE cell                                                 */
/* index 8 = SW cell                                                 */
/* Given a streamline falls within cell c, the quadrant and          */
/* indices to use are:                                               */
/* dx = x - cellx[c], dy = y - celly[c]                              */
/* dx < 0, dy > 0: quadrant 0, indices [0 1 5 2]                     */
/* dx > 0, dy > 0: quadrant 1, indices [0 2 6 3]                     */
/* dx > 0, dy < 0: quadrant 2, indices [0 3 7 4]                     */
/* dx < 0, dy < 0: quadrant 3, indices [0 1 8 4]                     */
/* Any index with value -1 gets the central index value.             */
void bl_rebuild(bl* l, point* p)
{
  int i;
  int id;
  double xp, yp;
  double q, r;
  delaunay *d = l->d;
  lweights *lw, *lws;

  for (i = 0; i < d->npoints; ++i) {
    id = (int)p[i].z;
    lw = &l->weights[i];
    lws = &l->weights[id];

    if(lw->isghost) continue;

    xp = p[i].x - lws->xr;
    yp = p[i].y - lws->yr;

    if (xp == 0) {
      if (yp >=0) 
	lw->qf = 1;
      else
	lw->qf = 2;
    } else if (yp == 0) {
      if (xp >=0) 
	lw->qf = 1;
      else
	lw->qf = 0;
    } else {
      if (xp < 0 && yp > 0) {
	lw->qf = 0;
      } else if (xp > 0 && yp > 0) {
	lw->qf = 1;
      } else if (xp > 0 && yp < 0) {
	lw->qf = 2;
      } else {
	lw->qf = 3;
      }
    }

    lw->qa = qpoints[lw->qf];
    q = fabs(xp) / lws->dx[lw->qf];
    r = fabs(yp) / lws->dy[lw->qf];
    q = max(min(q, 1.0), 0.0);
    r = max(min(r, 1.0), 0.0);
    /*if(i==623)printf("xp=%f yp=%f qf=%d q=%e r=%e %f %f\n",xp,yp,lw->qf,q,r,lws->dx[lw->qf],lws->dy[lw->qf]);*/
    lw->w[0] = (1.0 - q) * (1.0 - r);
    lw->w[1] = q * (1.0 - r);
    lw->w[2] = q * r;
    lw->w[3] = (1.0 - q) * r;
    /*if (i==623)printf("weights %f %f %f %f\n",lw->w[0],lw->w[1],lw->w[2],lw->w[3]);*/
  }
}

/* Accounts for interpolation on a streamline where the source and   */
/* destination locations may lie in different layers, and hence use  */
/* different grid_spec structures. Use l2 to get the coordinates and */
/* dimensions of the source cell, and l1 to place the weights for    */
/* the destination cell.                                             */
void bl_rebuild2(bl* l1, bl* l2, point* p)
{
  int i;
  int id;
  double xp, yp;
  double q, r;
  delaunay *d = l1->d;
  lweights *lw, *lws;

  int vid = d->vid;

  i = (int)p->v[vid];
  id = (int)p->z;

  lw = &l1->weights[i];
  lws = &l2->weights[id];

  if(lw->isghost) return;

  xp = p->x - lws->xr;
  yp = p->y - lws->yr;

  if (xp == 0) {
    if (yp >=0) 
      lw->qf = 1;
    else
      lw->qf = 2;
  } else if (yp == 0) {
    if (xp >=0) 
      lw->qf = 1;
    else
      lw->qf = 0;
  } else {
    if (xp < 0 && yp > 0) {
      lw->qf = 0;
    } else if (xp > 0 && yp > 0) {
      lw->qf = 1;
    } else if (xp > 0 && yp < 0) {
      lw->qf = 2;
    } else {
      lw->qf = 3;
    }
  }
  
  lw->qa = qpoints[lw->qf];
  q = fabs(xp) / lws->dx[lw->qf];
  r = fabs(yp) / lws->dy[lw->qf];
  q = max(min(q, 1.0), 0.0);
  r = max(min(r, 1.0), 0.0);
  /*if(i==623)printf("xp=%f yp=%f qf=%d q=%e r=%e %f %f\n",xp,yp,lw->qf,q,r,lws->dx[lw->qf],lws->dy[lw->qf]);*/
  lw->w[0] = (1.0 - q) * (1.0 - r);
  lw->w[1] = q * (1.0 - r);
  lw->w[2] = q * r;
  lw->w[3] = (1.0 - q) * r;
  /*if (i==623)printf("weights %f %f %f %f\n",lw->w[0],lw->w[1],lw->w[2],lw->w[3]);*/
}


/* Destroys quadratic interpolator.
 *
 * @param l Structure to be destroyed
 */
void bl_destroy(bl* l)
{
    i_free_1d(l->weights->cells);
    free(l->weights);
    free(l);
}

/* Finds lsq quadratic interpolated value in a point.
 *
 * @param l Quadratic interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void bl_interpolate_point(bl* l, point* p)
{
  int i, j;
  delaunay* d = l->d;
  lweights *lw, *lws;
  double v, z0;
  int id = (int)p->x;
  int ids = (int)p->y;
  int vid = d->vid;

  lw = &l->weights[id];
  lws = &l->weights[ids];
  if (id==643)printf("lin3 %x %x\n",lw->cells,lws->cells);
  v = 0.0;
  lw->tmx = -1e10;
  lw->tmn = 1e10;
  for (i = 0; i < 4; i++) {
    j = lws->cells[lw->qa[i]];
    z0 = (j == -1) ?  d->points[ids].v[vid] : d->points[j].v[vid];
    v += lw->w[i] * z0;
    /*if(vid==3&&id==623)printf("interp i=%d qf=%d qa=%d cells=%d z0=%f w=%e v=%f\n",i,lw->qf,lw->qa[i],j,z0,lw->w[i],v);*/
    lw->tmx = max(lw->tmx, z0);
    lw->tmn = min(lw->tmn, z0);
  }
  v = min(lw->tmx, max(lw->tmn, v));
  p->z = v;
}

void bl_interpolate_point2(bl* l1, bl* l2, point* p)
{
  int i, j;
  delaunay* d = l2->d;
  lweights *lw, *lws;
  double v, z0;
  int id = (int)p->x;
  int ids = (int)p->y;
  int vid = d->vid;

  lw = &l1->weights[id];
  lws = &l2->weights[ids];
  /*printf("a %x\n",lw->qa);*/
  /*if (lw->qa==NULL)return(d->points[ids].v[vid]);*/
  v = 0.0;
  lw->tmx = -1e10;
  lw->tmn = 1e10;
  for (i = 0; i < 4; i++) {
    /*printf(" a %d %d\n",i,lw->qa[i]);*/
    j = lws->cells[lw->qa[i]];
    /*printf(" b %d %d\n",i,j);*/
    z0 = (j == -1) ?  d->points[ids].v[vid] : d->points[j].v[vid];
    v += lw->w[i] * z0;
    /*if(vid==3&&id==623)printf("interp i=%d qf=%d qa=%d cells=%d z0=%f w=%e v=%f\n",i,lw->qf,lw->qa[i],j,z0,lw->w[i],v);*/
    lw->tmx = max(lw->tmx, z0);
    lw->tmn = min(lw->tmn, z0);
  }
  v = min(lw->tmx, max(lw->tmn, v));
  p->z = v;
  /*printf("ok\n");*/
}

/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void bl_interpolate_points(int nin, point pin[], int nout, point pout[])
{
    delaunay* d = delaunay_build(nin, pin, 0, NULL, 0, NULL);
    bl* l = bl_build(d);
    int seed = 0;
    int i;

    if (bl_verbose) {
        fprintf(stderr, "xytoi:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];

            fprintf(stderr, "(%.7g,%.7g) -> %d\n", p->x, p->y, delaunay_xytoi(d, p, seed));
        }
    }

    for (i = 0; i < nout; ++i)
        bl_interpolate_point(l, &pout[i]);

    if (bl_verbose) {
        fprintf(stderr, "output:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];;
            fprintf(stderr, "  %d:%15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    bl_destroy(l);
    delaunay_destroy(d);
}

// EOF
