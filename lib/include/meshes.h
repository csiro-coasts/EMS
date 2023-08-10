/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/meshes.h
 *  
 *  Description: Header for unstructured mesh construction
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: meshes.h 6595 2021-08-13 03:36:52Z her127 $
 *
 */

#if !defined(_MESHES_H)
#define _MESHES_H

typedef struct {
  int ns;                  /* # coordinates */
  int ns2;                 /* # cells */
  int mnpe;                /* Maximum vertices / cell */
  int *npe;                /* Vertices / cell array */
  int si;                  /* Start index */
  int type;                /* HEX or QUAD */
  double *xloc;            /* x coordinates */
  double *yloc;            /* y coordinates */
  int ***eloc;             /* Index to coordinate map */
  int nobc;                /* # Boundaries */
  int *npts;               /* # OBC cells */
  int **loc;               /* OBC index # */
  int ***obc;              /* OBC locations */  
  int nce1, nce2;          /* Cartesian grid size (optional) */
  int *iloc;               /* i location (optional) */
  int *jloc;               /* j location (optional) */
  int **neic;              /* Cell neighbour map */
  int **neij;              /* Cell neighbor index map */
  int *map;                /* Dummy mapping function */
  delaunay *d;             /* Triangulation for xytoc */
  int nd;                  /* Number of points in d */
  int *tri2c;              /* Map from d to c */
  poly_t *ppl;             /* Perimeter polygon */
  int np;                  /* Number of perimeter indices */     
  int *npedge;             /* Number of perimeter points */
  point **pedge;           /* Perimeter points */
} meshs_t;

void create_meshs(meshs_t *m, double **x, double **y);
void meshs_free(meshs_t *m);
void neighbour_finder_meshs(meshs_t *m);

#endif
