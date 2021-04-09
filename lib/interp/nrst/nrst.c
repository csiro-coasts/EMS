/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/nrst/nrst.c
 *  
 *  Description: 2D nearest neighbour interpolation
 *
 * `nrst' -- "Nearest Neighbour Interpolator" -- is a structure for
 * conducting nearest neighbour interpolation on a given data on a
 * "point-to-point" basis. It finds a triangle a point resides in, 
 * resulted from the Delaunay triangluation of input
 * data and then finds the vertex closest to the point.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: nrst.c 6670 2020-09-23 06:30:22Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nrst.h"

typedef struct {
    double x[3];
    double y[3];
    double z[3];
    int n;
    double xp, yp;
} nweights;

struct nrst {
    delaunay* d;
    nweights* weights;
};

//
int nrst_verbose = 0;

/* Builds nearest neighbour interpolator.
 *
 * @param d Delaunay triangulation
 * @return Nearest Neighbour interpolator
 */
nrst* nrst_build(delaunay* d)
{
    int i;
    nrst* n = malloc(sizeof(nrst));

    n->d = d;
    n->weights = malloc(d->ntriangles * sizeof(nweights));

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];
        nweights* nw = &n->weights[i];
	nw->x[0] = d->points[t->vids[0]].x;
	nw->x[1] = d->points[t->vids[1]].x;
	nw->x[2] = d->points[t->vids[2]].x;
	nw->y[0] = d->points[t->vids[0]].y;
	nw->y[1] = d->points[t->vids[1]].y;
	nw->y[2] = d->points[t->vids[2]].y;
        nw->z[0] = d->points[t->vids[0]].z;
        nw->z[1] = d->points[t->vids[1]].z;
        nw->z[2] = d->points[t->vids[2]].z;
	nw->xp = NaN;
	nw->yp = NaN;
    }
    return n;
}


/* Remakes the weights from cached weights */
void nrst_rebuild(nrst* n, point* p)
{
    int i;
    delaunay *d = n->d;

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];
        nweights* nw = &n->weights[i];
        nw->z[0] = d->points[t->vids[0]].z;
        nw->z[1] = d->points[t->vids[1]].z;
        nw->z[2] = d->points[t->vids[2]].z;
    }
}


/* Destroys nearest neighbour interpolator.
 *
 * @param l Structure to be destroyed
 */
void nrst_destroy(nrst* n)
{
    free(n->weights);
    free(n);
}

/* Finds nearest neighbour interpolated value in a point.
 *
 * @param l Nearest Neighbour interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void nrst_interpolate_point(nrst* n, point* p)
{
    delaunay* d = n->d;
    int i, tid = delaunay_xytoi_ng(d, p, d->first_id);
    double x, y, dist, dm;

    if (tid >= 0) {
        nweights* nw = &n->weights[tid];
	if (p->x == nw->xp && p->y == nw->yp) {
	  p->z = nw->z[nw->n];
	  return;
	}
	dm = 1e30;
	nw->xp = p->x;
	nw->yp = p->y;
	for (i = 0; i < 3; i++) {
	  x = nw->x[i] - p->x;
	  y = nw->y[i] - p->y;
	  dist = sqrt(x * x + y * y);
	  if (dist < dm) {
	    dm = dist;
	    p->z = nw->z[i];
	    nw->n = i;
	  }
	}
    } else {
        p->z = NaN;
    }
}


/* Nearest Neighbour interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void nrst_interpolate_points(int nin, point pin[], int nout, point pout[])
{
    delaunay* d = delaunay_build(nin, pin, 0, NULL, 0, NULL);
    nrst* n = nrst_build(d);
    int seed = 0;
    int i;

    if (nrst_verbose) {
        fprintf(stderr, "xytoi:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];

            fprintf(stderr, "(%.7g,%.7g) -> %d\n", p->x, p->y, delaunay_xytoi(d, p, seed));
        }
    }

    for (i = 0; i < nout; ++i)
        nrst_interpolate_point(n, &pout[i]);

    if (nrst_verbose) {
        fprintf(stderr, "output:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];;
            fprintf(stderr, "  %d:%15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    nrst_destroy(n);
    delaunay_destroy(d);
}

// EOF
