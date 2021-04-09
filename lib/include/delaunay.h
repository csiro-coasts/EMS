/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/delaunay.h
 *  
 *  Description: Header for delaunay triangulation wrapper
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: delaunay.h 6595 2020-09-03 03:36:52Z riz008 $
 *
 */

#if !defined(_DELAUNAY_H)
#define _DELAUNAY_H

#include "grid_utils.h"
#include "istack.h"

#if !defined(_DELAUNAY_STRUCT)
#define _DELAUNAY_STRUCT
typedef struct {
    int vids[3];
} triangle;

typedef struct {
    int tids[3];
} triangle_neighbours;

/** Structure to perform the Delaunay triangulation of a given array of points.
 *
 * Contains a deep copy of the input array of points.
 * Contains triangles, circles and edges resulted from the triangulation.
 * Contains neighbour triangles for each triangle.
 * Contains point to triangle map.
 */
typedef struct {
    int npoints;
    point* points;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int nvoints;
    point* voints;

    int ntriangles;
    triangle* triangles;
    circle* circles;
    triangle_neighbours* neighbours;    /* for delaunay_xytoi() */

    int ptf;                    /* Overwrite point_trianlges */
    int* n_point_triangles;     /* n_point_triangles[i] is number of
                                 * triangles i-th point belongs to */
    int** point_triangles;      /* point_triangles[i][j] is index of j-th
                                 * triangle i-th point belongs to */

    int nedges;
    int* edges;                 /* n-th edge is formed by points[edges[n*2]]
                                 * and points[edges[n*2+1]] */
    int nvdges;
    int* vdges;                 /* n-th edge is formed by points[edges[n*2]]
                                 * and points[edges[n*2+1]] */

    /*
     * Work data for delaunay_circles_find(). Placed here for efficiency
     * reasons. Should be moved to the procedure if parallelizable code
     * needed. 
     */
    int* flags;
    int first_id;               /* last search result, used in start up of a
                                 * new search */
    istack* t_in;
    istack* t_out;

    /*
     * to keep track of flags set to 1 in the case of very large data sets
     */
    int nflags;
    int nflagsallocated;
    int* flagids;
  int vid;    /* Id for interpolations */
} delaunay;
#endif

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
delaunay* delaunay_build(int np, point points[], int ns, int segments[], int nh, double holes[]);
delaunay* delaunay_voronoi_build(int np, point points[], int ns, int segments[], int nh, double holes[], double a, char *code);
delaunay* delaunay_create();

/* Destroys Delaunay triangulation.
 *
 * @param d Structure to be destroyed
 */
void delaunay_destroy(delaunay* d);


/* Finds triangle specified point belongs to (if any).
 *
 * @param d Delaunay triangulation
 * @param p Point to be mapped
 * @param seed Triangle index to start with
 * @return Triangle id if successful, -1 otherwhile
 */
int delaunay_xytoi(delaunay* d, point* p, int id);
/* Same as above, but attempts a no-gradient result if 
 * searches beyond the grid perimeter.
 */
int delaunay_xytoi_ng(delaunay* d, point* p, int id);
/* Same as above, designed to use with the semi-Lagrange
 * advection.
 */
int delaunay_xytoi_lag(delaunay* d, point* p, int id);

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
void delaunay_circles_find(delaunay* d, point* p, int* n, int** out);


#endif
