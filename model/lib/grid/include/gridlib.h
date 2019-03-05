/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/gridlib.h
 *  
 *  Description: Public header file for the EMS grid library
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridlib.h 6102 2019-02-08 05:15:14Z her127 $
 *
 */

#if !defined(_GRIDLIB_H)
#define _GRIDLIB_H

#include "ems.h"/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File:
 *  
 *  Description:
 *  
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridlib.h 6102 2019-02-08 05:15:14Z her127 $
 *
 */


typedef enum { 
  GRID_CSA             = 0, 
  GRID_NN_SIBSON       = 1, 
  GRID_NN_NONSIBSONIAN = 2,
  GRID_LINEAR          = 4,
  GRID_AVERAGE         = 8,
  GRID_LSQQ            = 16,
  GRID_BL              = 32,
  GRID_BAL             = 64,
  GRID_LSQL            = 128
} INTERP_RULE;

#include <float.h>
#include <math.h>
#include "config.h"
#include "ems.h"
#include "csa.h"
#include "delaunay.h"
#include "grid_utils.h"
#include "gridaverager.h"
#include "gridmap.h"
#include "gridnodes.h"
#include "minell.h"
#include "nnpi.h"
#include "nnai.h"
#include "lpi.h"
#include "bilinear.h"
#include "baycentric.h"
#include "lsqq.h"
#include "lsql.h"
#include "preader.h"
#include "svd.h"
#include "triangle.h"
#include "gridlib_version.h"

// interpolator function declaration
typedef void (*interp_func) (void*, point *);
typedef void (*interp_func2) (void*, void*, point *);
typedef void (*rebuild_func) (void*, point *);
typedef void (*rebuild_func2) (void*, void*, point *);

typedef struct {
  int    nbathy;      /* number of input points */
  point *pbathy;      /* input points */

  int    npout;      /* number of output points */
  point *pout;       /* output points */

  int destroy_pbathy;  // whether or not to destroy the pbathy array

  /* Pointer to the grid node struct */
  gridnodes *gn; 

  /* Type of interpolation */
  INTERP_RULE type;

  /*
   * Optional arguments. i.e. ones that have a default
   */

  /* Any mask to apply */
  int **mask;

  /* Double density, center or corner. see gridnodes.h */
  NODETYPE node_type;

  /* Points per edge */
  int ppe;

  /* Depth range */
  double zmin;
  double zmax;

  /* Estimate depth for this cell only */
  int cell_i;
  int cell_j;

  /* Whether to do the interpolation in index space */
  int index_space;

  /* Cache of the interpolator and the interpolate function */
  interp_func interpolate_point;
  interp_func2 interpolate_point2;
  rebuild_func rebuild;
  rebuild_func2 rebuild2;
  void       *interpolator;
  delaunay   *d;
  int nz;      /* Number of layers */
  double *z;   /* Depth levels */
  int id;      /* Optional Delaunay id for interpolation */

  /* This is the output vector */
  FILE *output_fileId;

} GRID_SPECS;

#define PPE_DEF 3
#define PPE_MAX 10
#define ZMIN_DEF (-DBL_MAX)
#define ZMAX_DEF DBL_MAX

/* Main entry function */
int grid_interp(GRID_SPECS *gs);

/* Create and destroy functions for the spec structure */
GRID_SPECS *grid_spec_create(void);
void grid_specs_destroy(GRID_SPECS *gs);
void grid_spec_init(GRID_SPECS *gs);

/* Convenient wrapper functions */
GRID_SPECS *grid_interp_init(double *x, double *y, double *z, int npoints,
			     char *rule);
void grid_interp_init_t(GRID_SPECS *gs, delaunay *d, char *rule, int var);
double grid_interp_on_point(GRID_SPECS *gs, double xcoord, double ycoord);
double grid_interp_on_point2(GRID_SPECS **gs, int k1, int k2, double xcoord, double ycoord);
double grid_interp_on_point3d(GRID_SPECS **gs, double xcoord, double ycoord, double depth, double bot);
void grid_interp_reinit(GRID_SPECS *gs, double *z, int npoints);

#endif

// EOF
