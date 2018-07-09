/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/utils.h
 *  
 *  Description: Header file for grid_utils
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: grid_utils.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_GRID_UTILS_H)
#define _GRID_UTILS_H

#define WMIN_DEF (-DBL_MAX)

/*********/
/* POINT */
/*********/

/* 
 * "point" is a basic data structure in this package.
 */
#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct {
  double x;
  double y;
  double z;
  double *v;
} point;
#endif

/*
 * Generic quit function for the grid library
 */
void grid_quit(char* format, ...);

/* Smoothes the input point array by averaging the input x, y and z values
 * for each cell within virtual rectangular nx by ny grid. The corners of the
 * grid are created from min and max values of the input array. It also frees
 * the original array and returns results and new dimension via original
 * data and size pointers. 
 *
 * @param pn Pointer to number of points (input/output)
 * @param ppoints Pointer to array of points (input/output) [*pn]
 * @param nx Number of x nodes in decimation
 * @param ny Number of y nodes in decimation
 */
void points_thingrid(int* pn, point** ppoints, int nx, int ny);

/* Smoothes the input point array by averaging the input data (X,Y and Z
 * values) until the sum of the distances between points does not exceed the
 * specified maximum value. It also frees the original array and returns
 * results and new dimension via original data and size pointers. 
 *
 * @param pn Pointer to number of points (input/output)
 * @param ppoints Pointer to array of points (input/output) [*pn]
 * @param rmax Maximum allowed accumulated distance
 */
void points_thinlin(int* nin, point** pin, double rmax);

/* Calculates X and/or Y ranges of the input array of points. If necessary,
 * adjusts the range according to the zoom value.
 *
 * @param n Number of points
 * @param points Array of points
 * @param xmin Min X value if *xmin = NaN on input, not changed otherwise
 * @param xmax Max X value if *xmax = NaN on input, not changed otherwise
 * @param ymin Min Y value if *ymin = NaN on input, not changed otherwise
 * @param ymax Max Y value if *ymax = NaN on input, not changed otherwise
 */
void points_getrange(int n, point points[], double zoom, double* xmin, double* xmax, double* ymin, double* ymax);

/* Generates rectangular grid nx by ny using specified min and max x and y 
 * values. Allocates space for the output point array, be sure to free it
 * when necessary!
 *
 * @param xmin Min x value
 * @param xmax Max x value
 * @param ymin Min y value
 * @param ymax Max y value
 * @param nx Number of x nodes
 * @param ny Number of y nodes
 * @param nout Pointer to number of output points
 * @param pout Pointer to array of output points [*nout]
 */
void points_generate(double xmin, double xmax, double ymin, double ymax, int nx, int ny, int* nout, point** pout);

/* Reads array of points from a columnar file.
 *
 * @param fname File name (can be "stdin" or "-" for standard input)
 * @param dim Number of dimensions (must be 2 or 3)
 * @param n Pointer to number of points (output)
 * @param points Pointer to array of points [*n] (output) (to be freed)
 * @param std Pointer to array of data std (to be freed)
 */
void points_read(char* fname, int dim, int* n, point** points, double** std);

/** Writes array of points to std out
 *
 * @param n Pointer to number of points
 * @param points Pointer to array of points [*n]
 */
void points_write(FILE *f, int n, point* points);

/** Scales Y coordinate so that the resulting set fits into square:
 ** xmax - xmin = ymax - ymin
 *
 * @param n Number of points
 * @param points The points to scale
 * @return Y axis compression coefficient
 */
double points_scaletosquare(int n, point* points);

/** Compresses Y domain by a given multiple.
 *
 * @param n Number of points
 * @param points The points to scale
 * @param Y axis compression coefficient as returned by points_scaletosquare()
 */
void points_scale(int n, point* points, double k);

/* 
 * Calculates whether the point p is on the left side of vector [p0, p1].
 * Returns 1 if on the left side, -1 if on the right side, 0 if on the line.
 */
int points_onleftside(point* p, point* p0, point* p1);

/* 
 * Moves points within in an array (of up to 5 points) so that they become
 * aranged in counterclockwise direction.
 */
void points_makeccw(int n, point* points[]);

/* 
 * Moves the last point to the front of the point pointer array.
 */
void points_movelasttofront(point** points, int n);

point** points_shuffle(int n, point* p, unsigned int seed);

double points_distance_xy(point* p1, point* p2);

point *points_create_and_fill(double *x, double *y, double *z, int n);

/**********/
/* CIRCLE */
/**********/
#if !defined(_CIRCLE_STRUCT)
#define _CIRCLE_STRUCT
typedef struct {
  double x;
  double y;
  double r;
} circle;
#endif

int circle_contains(circle* c, point* p);
int circle_build(circle* c, point* p1, point* p2, point* p3, int flag);


#endif /* _UTILS_H */
