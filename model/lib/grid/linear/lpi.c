/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/linear/lpi.c
 *  
 *  Description: 2D linear interpolation
 *
 * `lpi' -- "Linear Point Interpolator" -- is a structure for
 * conducting linear interpolation on a given data on a
 * "point-to-point" basis. It interpolates linearly within each
 * triangle resulted from the Delaunay triangluation of input
 * data. `lpi' is much faster than all Natural Neighbours
 * interpolators in `nn' library.
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lpi.c 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "lpi.h"

typedef struct {
    double w[3];
} lweights;

typedef struct {
  int flag;
  double w[6];
} lcweights;

struct lpi {
    delaunay* d;
    lweights* weights;
    lcweights* cweights;
};

//
int lpi_verbose = 0;

/* Builds linear interpolator.
 *
 * @param d Delaunay triangulation
 * @return Linear interpolator
 */
lpi* lpi_build(delaunay* d)
{
    int i;
    lpi* l = malloc(sizeof(lpi));

    l->d = d;
    l->weights = malloc(d->ntriangles * sizeof(lweights));
    l->cweights = malloc(d->ntriangles * sizeof(lcweights));

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];
        lweights* lw = &l->weights[i];
        lcweights* cw = &l->cweights[i];
        double x0 = d->points[t->vids[0]].x;
        double y0 = d->points[t->vids[0]].y;
        double z0 = d->points[t->vids[0]].z;
        double x1 = d->points[t->vids[1]].x;
        double y1 = d->points[t->vids[1]].y;
        double z1 = d->points[t->vids[1]].z;
        double x2 = d->points[t->vids[2]].x;
        double y2 = d->points[t->vids[2]].y;
        double z2 = d->points[t->vids[2]].z;
        double x02 = x0 - x2;
        double y02 = y0 - y2;
        double z02 = z0 - z2;
        double x12 = x1 - x2;
        double y12 = y1 - y2;
        double z12 = z1 - z2;

        if (y12 != 0.0) {
            double y0212 = y02 / y12;

            lw->w[0] = (z02 - z12 * y0212) / (x02 - x12 * y0212);
            lw->w[1] = (z12 - lw->w[0] * x12) / y12;
            lw->w[2] = (z2 - lw->w[0] * x2 - lw->w[1] * y2);

	    cw->w[0] = 1.0 / (x02 - x12 * y0212);
	    cw->w[1] = y0212;
	    cw->w[2] = 1.0 / y12;
	    cw->w[3] = x12;
	    cw->w[4] = x2;
	    cw->w[5] = y2;
	    cw->flag = 1;
        } else {
            double x0212 = x02 / x12;

            lw->w[1] = (z02 - z12 * x0212) / (y02 - y12 * x0212);
            lw->w[0] = (z12 - lw->w[1] * y12) / x12;
            lw->w[2] = (z2 - lw->w[0] * x2 - lw->w[1] * y2);

	    cw->w[0] = 1.0 / (y02 - y12 * x0212);
	    cw->w[1] = x0212;
	    cw->w[2] = 1.0 / x12;
	    cw->w[3] = y12;
	    cw->w[4] = x2;
	    cw->w[5] = y2;
	    cw->flag = 0;
        }
    }

    return l;
}


/* Remakes the weights from cached weights */
void lpi_rebuild(lpi* l, point* p)
{
    int i;
    delaunay *d = l->d;

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];
        lweights* lw = &l->weights[i];
        lcweights* cw = &l->cweights[i];
        double z0 = p[t->vids[0]].z;
        double z1 = p[t->vids[1]].z;
        double z2 = p[t->vids[2]].z;
        double z02 = z0 - z2;
        double z12 = z1 - z2;

	if (cw->flag) {
	  lw->w[0] = (z02 - z12 * cw->w[1]) * cw->w[0];
	  lw->w[1] = (z12 - lw->w[0] * cw->w[3]) * cw->w[2];
	  lw->w[2] = (z2 - lw->w[0] * cw->w[4] - lw->w[1] * cw->w[5]);
	} else {
	  lw->w[1] = (z02 - z12 * cw->w[1]) * cw->w[0];
	  lw->w[0] = (z12 - lw->w[1] * cw->w[3]) * cw->w[2];
	  lw->w[2] = (z2 - lw->w[0] * cw->w[4] - lw->w[1] * cw->w[5]);
	}
    }
}

void lpi_rebuild_tri(lpi* l, point* p, int tid, int k)
{

  delaunay *d = l->d;
  triangle* t = &d->triangles[tid];
  lweights* lw = &l->weights[tid];
  lcweights* cw = &l->cweights[tid];

  double z0 = p[t->vids[0]].v[k];
  double z1 = p[t->vids[1]].v[k];
  double z2 = p[t->vids[2]].v[k];
  double z02 = z0 - z2;
  double z12 = z1 - z2;
  printf("rebuild %f %f %f %d %d %d\n",z0,z1,z2,tid,k,t->vids[0]);
  if (cw->flag) {
    lw->w[0] = (z02 - z12 * cw->w[1]) * cw->w[0];
    lw->w[1] = (z12 - lw->w[0] * cw->w[3]) * cw->w[2];
    lw->w[2] = (z2 - lw->w[0] * cw->w[4] - lw->w[1] * cw->w[5]);
  } else {
    lw->w[1] = (z02 - z12 * cw->w[1]) * cw->w[0];
    lw->w[0] = (z12 - lw->w[1] * cw->w[3]) * cw->w[2];
    lw->w[2] = (z2 - lw->w[0] * cw->w[4] - lw->w[1] * cw->w[5]);
  }
}

/* Destroys linear interpolator.
 *
 * @param l Structure to be destroyed
 */
void lpi_destroy(lpi* l)
{
    free(l->weights);
    free(l);
}

/* Finds linearly interpolated value in a point.
 *
 * @param l Linear interpolation
 * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
 */
void lpi_interpolate_point(lpi* l, point* p)
{
    delaunay* d = l->d;
    int tid = delaunay_xytoi_ng(d, p, d->first_id);

    if (tid >= 0) {
        lweights* lw = &l->weights[tid];

        d->first_id = tid;
        p->z = p->x * lw->w[0] + p->y * lw->w[1] + lw->w[2];
    } else {
        p->z = NaN;
    }
}


/* Linearly interpolates data in an array of points.
 *
 * @param nin Number of input points
 * @param pin Array of input points [pin]
 * @param nout Number of ouput points
 * @param pout Array of output points [nout]
 */
void lpi_interpolate_points(int nin, point pin[], int nout, point pout[])
{
    delaunay* d = delaunay_build(nin, pin, 0, NULL, 0, NULL);
    lpi* l = lpi_build(d);
    int seed = 0;
    int i;

    if (lpi_verbose) {
        fprintf(stderr, "xytoi:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];

            fprintf(stderr, "(%.7g,%.7g) -> %d\n", p->x, p->y, delaunay_xytoi(d, p, seed));
        }
    }

    for (i = 0; i < nout; ++i)
        lpi_interpolate_point(l, &pout[i]);

    if (lpi_verbose) {
        fprintf(stderr, "output:\n");
        for (i = 0; i < nout; ++i) {
            point* p = &pout[i];;
            fprintf(stderr, "  %d:%15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    lpi_destroy(l);
    delaunay_destroy(d);
}

// EOF
