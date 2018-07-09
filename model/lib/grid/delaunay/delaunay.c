/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/delaunay/delaunay.c
 *  
 *  Description:
 *  Delaunay triangulation - a wrapper to triangulate()
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: delaunay.c 5859 2018-07-02 04:04:44Z her127 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "triangle.h"
#include "istack.h"
#include "nan.h"
#include "delaunay.h"


//
int verbose = 0;

int intri(delaunay *d, int id, point *p);

/*
 * This parameter is used in search of tricircles containing a given point:
 *   if there are no more triangles than N_SEARCH_TURNON
 *     do linear search
 *   else
 *     do more complicated stuff
 */
#define N_SEARCH_TURNON 20
#define N_FLAGS_TURNON 1000
#define N_FLAGS_INC 100

static void tio_init(struct triangulateio* tio)
{
    tio->pointlist = NULL;
    tio->pointattributelist = NULL;
    tio->pointmarkerlist = NULL;
    tio->numberofpoints = 0;
    tio->numberofpointattributes = 0;
    tio->trianglelist = NULL;
    tio->triangleattributelist = NULL;
    tio->trianglearealist = NULL;
    tio->neighborlist = NULL;
    tio->numberoftriangles = 0;
    tio->numberofcorners = 0;
    tio->numberoftriangleattributes = 0;
    tio->segmentlist = 0;
    tio->segmentmarkerlist = NULL;
    tio->numberofsegments = 0;
    tio->holelist = NULL;
    tio->numberofholes = 0;
    tio->regionlist = NULL;
    tio->numberofregions = 0;
    tio->edgelist = NULL;
    tio->edgemarkerlist = NULL;
    tio->normlist = NULL;
    tio->numberofedges = 0;
}

static void tio_destroy(struct triangulateio* tio)
{
    if (tio->pointlist != NULL)
        free(tio->pointlist);
    if (tio->pointattributelist != NULL)
        free(tio->pointattributelist);
    if (tio->pointmarkerlist != NULL)
        free(tio->pointmarkerlist);
    if (tio->trianglelist != NULL)
        free(tio->trianglelist);
    if (tio->triangleattributelist != NULL)
        free(tio->triangleattributelist);
    if (tio->trianglearealist != NULL)
        free(tio->trianglearealist);
    if (tio->neighborlist != NULL)
        free(tio->neighborlist);
    if (tio->segmentlist != NULL)
        free(tio->segmentlist);
    if (tio->segmentmarkerlist != NULL)
        free(tio->segmentmarkerlist);
    if (tio->holelist != NULL)
        free(tio->holelist);
    if (tio->regionlist != NULL)
        free(tio->regionlist);
    if (tio->edgelist != NULL)
        free(tio->edgelist);
    if (tio->edgemarkerlist != NULL)
        free(tio->edgemarkerlist);
    if (tio->normlist != NULL)
        free(tio->normlist);
}

delaunay* delaunay_create()
{
    delaunay* d = malloc(sizeof(delaunay));

    d->npoints = 0;
    d->points = NULL;
    d->xmin = DBL_MAX;
    d->xmax = -DBL_MAX;
    d->ymin = DBL_MAX;
    d->ymax = -DBL_MAX;
    d->ntriangles = 0;
    d->triangles = NULL;
    d->circles = NULL;
    d->neighbours = NULL;
    d->n_point_triangles = NULL;
    d->point_triangles = NULL;
    d->nedges = 0;
    d->edges = NULL;
    d->flags = NULL;
    d->first_id = -1;
    d->t_in = NULL;
    d->t_out = NULL;
    d->nflags = 0;
    d->nflagsallocated = 0;
    d->flagids = NULL;

    return d;
}

static void tio2delaunay(struct triangulateio* tio_out, delaunay* d)
{
    int i, j;

    /*
     * I assume that all input points appear in tio_out in the same order as 
     * they were written to tio_in. I have seen no exceptions so far, even
     * if duplicate points were presented. Just in case, let us make a couple
     * of checks. 
     */
    assert(tio_out->numberofpoints == d->npoints);
    assert(tio_out->pointlist[2 * d->npoints - 2] == d->points[d->npoints - 1].x && tio_out->pointlist[2 * d->npoints - 1] == d->points[d->npoints - 1].y);

    for (i = 0, j = 0; i < d->npoints; ++i) {
        point* p = &d->points[i];

        if (p->x < d->xmin)
            d->xmin = p->x;
        if (p->x > d->xmax)
            d->xmax = p->x;
        if (p->y < d->ymin)
            d->ymin = p->y;
        if (p->y > d->ymax)
            d->ymax = p->y;
    }
    if (verbose) {
        fprintf(stderr, "input:\n");
        for (i = 0, j = 0; i < d->npoints; ++i) {
            point* p = &d->points[i];

            fprintf(stderr, "  %d: %15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    d->ntriangles = tio_out->numberoftriangles;
    if (d->ntriangles > 0) {
        d->triangles = malloc(d->ntriangles * sizeof(triangle));
        d->neighbours = malloc(d->ntriangles * sizeof(triangle_neighbours));
        d->circles = malloc(d->ntriangles * sizeof(circle));
        d->n_point_triangles = calloc(d->npoints, sizeof(int));
        d->point_triangles = malloc(d->npoints * sizeof(int*));
        d->flags = calloc(d->ntriangles, sizeof(int));
    }

    if (verbose)
        fprintf(stderr, "triangles:\n");
    for (i = 0; i < d->ntriangles; ++i) {
        int offset = i * 3;
        triangle* t = &d->triangles[i];
        triangle_neighbours* n = &d->neighbours[i];
        circle* c = &d->circles[i];
        int status;

        t->vids[0] = tio_out->trianglelist[offset];
        t->vids[1] = tio_out->trianglelist[offset + 1];
        t->vids[2] = tio_out->trianglelist[offset + 2];
        n->tids[0] = tio_out->neighborlist[offset];
        n->tids[1] = tio_out->neighborlist[offset + 1];
        n->tids[2] = tio_out->neighborlist[offset + 2];

        status = circle_build(c, &d->points[t->vids[0]], &d->points[t->vids[1]], &d->points[t->vids[2]], 0);
        assert(status);

        if (verbose)
            fprintf(stderr, "  %d: (%d,%d,%d)\n", i, t->vids[0], t->vids[1], t->vids[2]);
    }

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j)
            d->n_point_triangles[t->vids[j]]++;
    }
    if (d->ntriangles > 0) {
        for (i = 0; i < d->npoints; ++i) {
            if (d->n_point_triangles[i] > 0)
                d->point_triangles[i] = malloc(d->n_point_triangles[i] * sizeof(int));
            else
                d->point_triangles[i] = NULL;
            d->n_point_triangles[i] = 0;
        }
    }
    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j) {
            int vid = t->vids[j];

            d->point_triangles[vid][d->n_point_triangles[vid]] = i;
            d->n_point_triangles[vid]++;
        }
    }

    if (tio_out->edgelist != NULL) {
        d->nedges = tio_out->numberofedges;
        d->edges = malloc(d->nedges * 2 * sizeof(int));
        memcpy(d->edges, tio_out->edgelist, d->nedges * 2 * sizeof(int));
    }
}

static void tio2delaunay_v(struct triangulateio* tio_out, struct triangulateio* vio_out, delaunay* d)
{
    int i, j;

    /*
     * I assume that all input points appear in tio_out in the same order as 
     * they were written to tio_in. I have seen no exceptions so far, even
     * if duplicate points were presented. Just in case, let us make a couple
     * of checks. 
     */
    /*
    assert(tio_out->numberofpoints == d->npoints);
    assert(tio_out->pointlist[2 * d->npoints - 2] == d->points[d->npoints - 1].x && tio_out->pointlist[2 * d->npoints - 1] == d->points[d->npoints - 1].y);

    free((point *)d->points);
    */
    d->npoints = tio_out->numberofpoints;
    d->points = malloc(d->npoints * sizeof(point));
    for (i = 0, j = 0; i < d->npoints; ++i) {
      point* p = &d->points[i];
      int offset = i * 2;
      p->x = tio_out->pointlist[offset];
      p->y = tio_out->pointlist[offset + 1];
    }

    if (vio_out != NULL) {
      d->nvoints = vio_out->numberofpoints;
      d->voints = malloc(d->nvoints * sizeof(point));
      for (i = 0, j = 0; i < d->nvoints; ++i) {
	point* p = &d->voints[i];
	int offset = i * 2;
	p->x = vio_out->pointlist[offset];
	p->y = vio_out->pointlist[offset + 1];
      }
      if (vio_out->edgelist != NULL) {
        d->nvdges = vio_out->numberofedges;
        d->vdges = malloc(d->nvdges * 2 * sizeof(int));
        memcpy(d->vdges, vio_out->edgelist, d->nvdges * 2 * sizeof(int));
      }
    }

    for (i = 0, j = 0; i < d->npoints; ++i) {
        point* p = &d->points[i];

        if (p->x < d->xmin)
            d->xmin = p->x;
        if (p->x > d->xmax)
            d->xmax = p->x;
        if (p->y < d->ymin)
            d->ymin = p->y;
        if (p->y > d->ymax)
            d->ymax = p->y;
    }
    if (verbose) {
        fprintf(stderr, "input:\n");
        for (i = 0, j = 0; i < d->npoints; ++i) {
            point* p = &d->points[i];

            fprintf(stderr, "  %d: %15.7g %15.7g %15.7g\n", i, p->x, p->y, p->z);
        }
    }

    d->ntriangles = tio_out->numberoftriangles;
    if (d->ntriangles > 0) {
        d->triangles = malloc(d->ntriangles * sizeof(triangle));
        d->neighbours = malloc(d->ntriangles * sizeof(triangle_neighbours));
        d->circles = malloc(d->ntriangles * sizeof(circle));
        d->n_point_triangles = calloc(d->npoints, sizeof(int));
        d->point_triangles = malloc(d->npoints * sizeof(int*));
        d->flags = calloc(d->ntriangles, sizeof(int));
    }

    if (verbose)
        fprintf(stderr, "triangles:\n");
    for (i = 0; i < d->ntriangles; ++i) {
        int offset = i * 3;
        triangle* t = &d->triangles[i];
        triangle_neighbours* n = &d->neighbours[i];
        circle* c = &d->circles[i];
        int status;

        t->vids[0] = tio_out->trianglelist[offset];
        t->vids[1] = tio_out->trianglelist[offset + 1];
        t->vids[2] = tio_out->trianglelist[offset + 2];

        n->tids[0] = tio_out->neighborlist[offset];
        n->tids[1] = tio_out->neighborlist[offset + 1];
        n->tids[2] = tio_out->neighborlist[offset + 2];

        status = circle_build(c, &d->points[t->vids[0]], &d->points[t->vids[1]], &d->points[t->vids[2]], 0);
        assert(status);

        if (verbose)
            fprintf(stderr, "  %d: (%d,%d,%d)\n", i, t->vids[0], t->vids[1], t->vids[2]);
    }

    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j)
            d->n_point_triangles[t->vids[j]]++;
    }
    if (d->ntriangles > 0) {
        for (i = 0; i < d->npoints; ++i) {
            if (d->n_point_triangles[i] > 0)
                d->point_triangles[i] = malloc(d->n_point_triangles[i] * sizeof(int));
            else
                d->point_triangles[i] = NULL;
            d->n_point_triangles[i] = 0;
        }
    }
    for (i = 0; i < d->ntriangles; ++i) {
        triangle* t = &d->triangles[i];

        for (j = 0; j < 3; ++j) {
            int vid = t->vids[j];

            d->point_triangles[vid][d->n_point_triangles[vid]] = i;
            d->n_point_triangles[vid]++;
        }
    }

    if (tio_out->edgelist != NULL) {
        d->nedges = tio_out->numberofedges;
        d->edges = malloc(d->nedges * 2 * sizeof(int));
        memcpy(d->edges, tio_out->edgelist, d->nedges * 2 * sizeof(int));
    }
}

/* Builds Delaunay triangulation of the given array of points.
 *
 * @param np Number of points
 * @param points Array of points [np] (input)
 * @param ns Number of forced segments
 * @param segments Array of (forced) segment endpoint indices [2*ns]
 * @param nh Number of holes
 * @param holes Array of hole (x,y) coordinates [2*nh]
 * @return Delaunay triangulation structure with triangulation results
 */
delaunay* delaunay_build(int np, point points[], int ns, int segments[], int nh, double holes[])
{
    delaunay* d = delaunay_create();
    struct triangulateio tio_in;
    struct triangulateio tio_out;
    char cmd[64] = "eznC";
    int i, j;

    assert(sizeof(REAL) == sizeof(double));

    tio_init(&tio_in);

    if (np == 0) {
        free(d);
        return NULL;
    }

    tio_in.pointlist = malloc(np * 2 * sizeof(double));
    tio_in.numberofpoints = np;
    for (i = 0, j = 0; i < np; ++i) {
        tio_in.pointlist[j++] = points[i].x;
        tio_in.pointlist[j++] = points[i].y;
    }

    if (ns > 0) {
        tio_in.segmentlist = malloc(ns * 2 * sizeof(int));
        tio_in.numberofsegments = ns;
        memcpy(tio_in.segmentlist, segments, ns * 2 * sizeof(int));
    }

    if (nh > 0) {
        tio_in.holelist = malloc(nh * 2 * sizeof(double));
        tio_in.numberofholes = nh;
        memcpy(tio_in.holelist, holes, nh * 2 * sizeof(double));
    }

    tio_init(&tio_out);

    if (!verbose)
        strcat(cmd, "Q");
    else if (verbose > 1)
        strcat(cmd, "VV");
    if (ns != 0)
        strcat(cmd, "p");

    if (verbose)
        fflush(stderr);

    /*
     * climax 
     */
    triangulate(cmd, &tio_in, &tio_out, NULL);

    if (verbose)
        fflush(stderr);

    d->npoints = np;
    d->points = points;

    tio2delaunay(&tio_out, d);

    tio_destroy(&tio_in);
    tio_destroy(&tio_out);

    return d;
}

/* Destroys Delaunay triangulation.
 *
 * @param d Structure to be destroyed
 */
void delaunay_destroy(delaunay* d)
{
    if (d == NULL)
        return;

    if (d->point_triangles != NULL) {
        int i;

        for (i = 0; i < d->npoints; ++i)
            if (d->point_triangles[i] != NULL)
                free(d->point_triangles[i]);
        free(d->point_triangles);
    }
    if (d->nedges > 0)
        free(d->edges);
    if (d->n_point_triangles != NULL)
        free(d->n_point_triangles);
    if (d->flags != NULL)
        free(d->flags);
    if (d->circles != NULL)
        free(d->circles);
    if (d->neighbours != NULL)
        free(d->neighbours);
    if (d->triangles != NULL)
        free(d->triangles);
    if (d->t_in != NULL)
        istack_destroy(d->t_in);
    if (d->t_out != NULL)
        istack_destroy(d->t_out);
    if (d->flagids != NULL)
        free(d->flagids);
    free(d);
}

/* Returns whether the point p is on the right side of the vector (p0, p1).
 */
static int onrightside(point* p, point* p0, point* p1)
{
    return (p1->x - p->x) * (p0->y - p->y) > (p0->x - p->x) * (p1->y - p->y);
}

/* Finds triangle specified point belongs to (if any).
 *
 * @param d Delaunay triangulation
 * @param p Point to be mapped
 * @param seed Triangle index to start with
 * @return Triangle id if successful, -1 otherwhile
 */
int delaunay_xytoi(delaunay* d, point* p, int id)
{
    triangle* t;
    int i;

    if (p->x < d->xmin || p->x > d->xmax || p->y < d->ymin || p->y > d->ymax)
        return -1;

    if (id < 0 || id > d->ntriangles)
        id = 0;
    t = &d->triangles[id];
    do {
        for (i = 0; i < 3; ++i) {
            int i1 = (i + 1) % 3;
            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
                id = d->neighbours[id].tids[(i + 2) % 3];
                if (id < 0)
                    return id;
                t = &d->triangles[id];
                break;
            }
        }
    } while (i < 3);

    return id;
}

/* Finds triangle specified point belongs to.The triangle is found  */
/* by walking through the triangulation. If the walk goes outside   */
/* the triangulation, then the previous triangle in the walk is     */
/* returned.                                                        */
int delaunay_xytoi_ng(delaunay* d, point* p, int id)
{
    triangle* t;
    int i;
    int idp = id;

    if (p->x < d->xmin || p->x > d->xmax || p->y < d->ymin || p->y > d->ymax)
        return id;

    if (id < 0 || id > d->ntriangles)
        id = 0;
    t = &d->triangles[id];
    do {
        for (i = 0; i < 3; ++i) {
            int i1 = (i + 1) % 3;

            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
	        idp = id;
                id = d->neighbours[id].tids[(i + 2) % 3];

                if (id < 0)
		  return idp;

                t = &d->triangles[id];
                break;
            }
        }
    } while (i < 3);

    return id;
}

/* Finds triangle specified point belongs to. This is intended for   */
/* use with the Lagrange advection scheme. The triangle is found by  */
/* walking through the triangulation, and is only guaranteed to      */
/* terminate in the target triangle for Delaunay triangulations; for */
/* arbitrary triangulations it may end up in an infinite loop. If    */
/* this happens the loop is broken by selecting a random edge for    */
/* mapping to a neighbouring triangle. An infinite loop may also     */
/* result if the point lies on a vertex (i.e. the source cell in the */
/* Lagrange streamline tracing is the same as the streamline origin  */
/* (the destination cell)). The baycentric method is used to check   */
/* if the point les within each triangle in the walk.                */
/* For Lagrange streamline tracking, it helps to seed the first_id   */
/* with the triangle corresponding to the destination cell to        */
/* improve efficiency (i.e. reduce the length of the walk).          */
/* If the walk goes outside the triangulation, then an edge that     */
/* sends the walk back into the triangulation is found. This may     */
/* make the walk inefficiently long.                                 */
/* Note this routine may end in an infinite loop for triangulations  */
/* built from quad meshes (non Delaunay).                           */
int delaunay_xytoi_lag(delaunay* d, point* p, int id)
{
    triangle* t;
    int i;
    int ids, idp=id;

    if (p->x < d->xmin || p->x > d->xmax || p->y < d->ymin || p->y > d->ymax)
        return id;

    if (id < 0 || id > d->ntriangles)
        id = 0;
    t = &d->triangles[id];

    if(fabs(p->x - d->points[t->vids[0]].x) < 1e-5 &&
       fabs(p->y - d->points[t->vids[0]].y) < 1e-5) return(id);
    do {
        for (i = 0; i < 3; ++i) {
	  int i1 = (i + 1) % 3, i2;
	    int ran = 0;

	    if (intri(d, id, p)) return(id);
            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
                ids = d->neighbours[id].tids[(i + 2) % 3];

		if (ids == idp) {
		  ids = d->neighbours[id].tids[ran3(&ran)%3];
		}

		/*
		  printf("%d id=%d idp=%d ids=%d\n",i,id,idp,ids);
		  printf(" %f %f\n",d->points[t->vids[0]].x,d->points[t->vids[0]].y);
		  printf(" %f %f\n",d->points[t->vids[1]].x,d->points[t->vids[1]].y);
		  printf(" %f %f\n",d->points[t->vids[2]].x,d->points[t->vids[2]].y);
		*/

                if (ids < 0) {
		  int i2 = 0;
		  while((ids = d->neighbours[id].tids[i2]) < 0)
		    i2 = i2 % 3;
		}
	        idp = id;
		id = ids;

                t = &d->triangles[id];
		/*printf("%d %f %f\n",id,d->points[t->vids[0]].x,d->points[t->vids[0]].y);*/

                break;
            }
        }
    } while (i < 3);

    return id;
}


int delaunay_xytoi_lago(delaunay* d, point* p, int id)
{
    triangle* t;
    int i;
    int ids, idp=id;

    if (p->x < d->xmin || p->x > d->xmax || p->y < d->ymin || p->y > d->ymax)
        return id;

    if (id < 0 || id > d->ntriangles)
        id = 0;
    t = &d->triangles[id];

    if(fabs(p->x - d->points[t->vids[0]].x) < 1e-5 &&
       fabs(p->y - d->points[t->vids[0]].y) < 1e-5) return(id);
    do {
        for (i = 0; i < 3; ++i) {
            int i1 = (i + 1) % 3;
	    int ran = 0;

	    if (intri(d, id, p)) return(id);
            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
                ids = d->neighbours[id].tids[(i + 2) % 3];

		if (ids == idp) {
		  ids = d->neighbours[id].tids[ran3(&ran)%3];
		}
		/*printf("%d %d %d %d %d\n",i,id,idp,ids,ran3(&ran)%3);*/

                if (ids < 0) {
		  triangle *t1 = &d->triangles[id];
		  printf(" %f %f\n",d->points[t1->vids[0]].x,d->points[t1->vids[0]].y);
		  printf(" %f %f\n",d->points[t1->vids[1]].x,d->points[t1->vids[1]].y);
		  printf(" %f %f\n",d->points[t1->vids[2]].x,d->points[t1->vids[2]].y);
		  while((ids = d->neighbours[id].tids[ran3(&ran)%3]) < 0) {
		    printf("idp %d %d\n",idp,ids);
		  }
		  printf("end %d %d\n",idp,ids);
		}
	        idp = id;
		id = ids;

                t = &d->triangles[id];
		printf("%f %f\n",d->points[t->vids[0]].x,d->points[t->vids[0]].y);

                break;
            }
        }
    } while (i < 3);

    return id;
}

/* Returns 1 if a point lies within a triangle. Uses Baycentric     */
/* coordinates to find if this is so. The area is > 0 for vertices  */
/* ordered in an anti-clockwise sense, and negative for a clockwise */
/* sense.                                                           */
int intri(delaunay *d, int id, point *p) {
  triangle *tr = &d->triangles[id];
  double p0x = d->points[tr->vids[0]].x;
  double p0y = d->points[tr->vids[0]].y;
  double p1x = d->points[tr->vids[1]].x;
  double p1y = d->points[tr->vids[1]].y;
  double p2x = d->points[tr->vids[2]].x;
  double p2y = d->points[tr->vids[2]].y;
  double area = 0.5 * (-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
  double ap = 2.0*fabs(area);
  double s = (p0y*p2x - p0x*p2y + (p2y - p0y)*p->x + (p0x - p2x)*p->y)/(2.0*area);
  double t = (p0x*p1y - p0y*p1x + (p0y - p1y)*p->x + (p1x - p0x)*p->y)/(2.0*area);
  double sp = s * ap;
  double tp = t * ap;
  int ret = (sp > 0.0 & tp > 0.0 & sp + tp < ap) ? 1 : 0;
  return(ret);
}



static void delaunay_addflag(delaunay* d, int i)
{
    if (d->nflags == d->nflagsallocated) {
        d->nflagsallocated += N_FLAGS_INC;
        d->flagids = realloc(d->flagids, d->nflagsallocated * sizeof(int));
    }
    d->flagids[d->nflags] = i;
    d->nflags++;
}

static void delaunay_resetflags(delaunay* d)
{
    int i;

    for (i = 0; i < d->nflags; ++i)
        d->flags[d->flagids[i]] = 0;
    d->nflags = 0;
}

/* Finds all tricircles specified point belongs to.
 *
 * @param d Delaunay triangulation
 * @param p Point to be mapped
 * @param n Pointer to the number of tricircles within `d' containing `p'
 *          (output)
 * @param out Pointer to an array of indices of the corresponding triangles 
 *            [n] (output)
 *
 * There is a standard search procedure involving search through triangle
 * neighbours (not through vertex neighbours). It must be a bit faster due to
 * the smaller number of triangle neighbours (3 per triangle) but may fail
 * for a point outside convex hall.
 *
 * We may wish to modify this procedure in future: first check if the point
 * is inside the convex hall, and depending on that use one of the two
 * search algorithms. It not 100% clear though whether this will lead to a
 * substantial speed gains because of the check on convex hall involved.
 */
void delaunay_circles_find(delaunay* d, point* p, int* n, int** out)
{
    /*
     * This flag was introduced as a hack to handle some degenerate cases. It 
     * is set to 1 only if the triangle associated with the first circle is
     * already known to contain the point. In this case the circle is assumed 
     * to contain the point without a check. In my practice this turned
     * useful in some cases when point p coincided with one of the vertices
     * of a thin triangle. 
     */
    int contains = 0;
    int i;

    if (d->t_in == NULL) {
        d->t_in = istack_create();
        d->t_out = istack_create();
    }

    /*
     * if there are only a few data points, do linear search
     */
    if (d->ntriangles <= N_SEARCH_TURNON) {
        istack_reset(d->t_out);

        for (i = 0; i < d->ntriangles; ++i) {
            if (circle_contains(&d->circles[i], p)) {
                istack_push(d->t_out, i);
            }
        }

        *n = d->t_out->n;
        *out = d->t_out->v;

        return;
    }
    /*
     * otherwise, do a more complicated stuff
     */

    /*
     * It is important to have a reasonable seed here. If the last search
     * was successful -- start with the last found tricircle, otherwhile (i) 
     * try to find a triangle containing p; if fails then (ii) check
     * tricircles from the last search; if fails then (iii) make linear
     * search through all tricircles 
     */
    if (d->first_id < 0 || !circle_contains(&d->circles[d->first_id], p)) {
        /*
         * if any triangle contains p -- start with this triangle 
         */
        d->first_id = delaunay_xytoi(d, p, d->first_id);
        contains = (d->first_id >= 0);

        /*
         * if no triangle contains p, there still is a chance that it is
         * inside some of circumcircles 
         */
        if (d->first_id < 0) {
            int nn = d->t_out->n;
            int tid = -1;

            /*
             * first check results of the last search 
             */
            for (i = 0; i < nn; ++i) {
                tid = d->t_out->v[i];
                if (circle_contains(&d->circles[tid], p))
                    break;
            }
            /*
             * if unsuccessful, search through all circles 
             */
            if (tid < 0 || i == nn) {
                double nt = d->ntriangles;

                for (tid = 0; tid < nt; ++tid) {
                    if (circle_contains(&d->circles[tid], p))
                        break;
                }
                if (tid == nt) {
                    istack_reset(d->t_out);
                    *n = 0;
                    *out = NULL;
                    return;     /* failed */
                }
            }
            d->first_id = tid;
        }
    }

    istack_reset(d->t_in);
    istack_reset(d->t_out);

    istack_push(d->t_in, d->first_id);
    d->flags[d->first_id] = 1;
    delaunay_addflag(d, d->first_id);

    /*
     * main cycle 
     */
    while (d->t_in->n > 0) {
        int tid = istack_pop(d->t_in);
        triangle* t = &d->triangles[tid];

        if (contains || circle_contains(&d->circles[tid], p)) {
            istack_push(d->t_out, tid);
            for (i = 0; i < 3; ++i) {
                int vid = t->vids[i];
                int nt = d->n_point_triangles[vid];
                int j;

                for (j = 0; j < nt; ++j) {
                    int ntid = d->point_triangles[vid][j];

                    if (d->flags[ntid] == 0) {
                        istack_push(d->t_in, ntid);
                        d->flags[ntid] = 1;
                        delaunay_addflag(d, ntid);
                    }
                }
            }
        }
        contains = 0;
    }

    *n = d->t_out->n;
    *out = d->t_out->v;
    delaunay_resetflags(d);
}


/* Builds Delaunay triangulation of the given array of points.
 *
 * @param np Number of points
 * @param points Array of points [np] (input)
 * @param ns Number of forced segments
 * @param segments Array of (forced) segment endpoint indices [2*ns]
 * @param nh Number of holes
 * @param holes Array of hole (x,y) coordinates [2*nh]
 * @return Delaunay triangulation structure with triangulation results
 */
delaunay* delaunay_voronoi_build(int np, point points[], int ns, int segments[], int nh, double holes[], double a, char *code)
{
    delaunay* d = delaunay_create();
    struct triangulateio tio_in;
    struct triangulateio tio_out;
    struct triangulateio vio_out;
    char cmd[64];
    int i, j;
    int dovor = 1;

    assert(sizeof(REAL) == sizeof(double));

    if (a)
      sprintf(cmd, "eznCDa%.8f", a);
    else if(code != NULL)
      strcpy(cmd, code);
    else
      strcpy(cmd, "eznCDq");

    tio_init(&tio_in);

    if (np > 0) {
      tio_in.pointlist = malloc(np * 2 * sizeof(double));
      tio_in.numberofpoints = np;
      for (i = 0, j = 0; i < np; ++i) {
        tio_in.pointlist[j++] = points[i].x;
        tio_in.pointlist[j++] = points[i].y;
      }
    }

    if (ns > 0) {
        tio_in.segmentlist = malloc(ns * 2 * sizeof(int));
        tio_in.numberofsegments = ns;
        memcpy(tio_in.segmentlist, segments, ns * 2 * sizeof(int));
    }

    if (nh > 0) {
        tio_in.holelist = malloc(nh * 2 * sizeof(double));
        tio_in.numberofholes = nh;
        memcpy(tio_in.holelist, holes, nh * 2 * sizeof(double));
    }

    tio_init(&tio_out);
    tio_init(&vio_out);

    if (!verbose)
        strcat(cmd, "Q");
    else if (verbose > 1)
        strcat(cmd, "VV");
    if (ns != 0)
        strcat(cmd, "p");
    if (dovor)
      strcat(cmd, "v");
    if (verbose)
        fflush(stderr);

    /*
     * climax 
     */
    triangulate(cmd, &tio_in, &tio_out, &vio_out);

    if (verbose)
        fflush(stderr);
    /*
    d->npoints = np;
    d->points = points;
    */
    tio2delaunay_v(&tio_out, &vio_out, d);

    tio_destroy(&tio_in);
    tio_destroy(&tio_out);
    tio_destroy(&vio_out);

    return d;
}

// EOF
