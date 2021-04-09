/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/gridmap.h
 *  
 *  Description: Header file for gridmap.c
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: gridmap.h 6597 2020-09-03 05:29:26Z riz008 $
 *
 */

#if !defined(_GRIDMAP_H)
#define _GRIDMAP_H

struct gridmap;
typedef struct gridmap gridmap;

// Make consistent with ems
typedef poly_t poly;
typedef extent_t extent;

gridmap* gridmap_build(int nce1, int nce2, double** gx, double** gy);
void gridmap_destroy(gridmap* gm);
int gridmap_fij2xy(gridmap* gm, double fi, double fj, double* x, double* y);
int gridmap_xy2ij(gridmap* gm, double x, double y, int* i, int* j);
int gridmap_xy2fij(gridmap* gm, double x, double y, double* fi, double* fj);
int gridmap_getnce1(gridmap* gm);
int gridmap_getnce2(gridmap* gm);
void gridmap_getextent(gridmap* gm, double* xmin, double* xmax, double* ymin, double* ymax);

#endif
