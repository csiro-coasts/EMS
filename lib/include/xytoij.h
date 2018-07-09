/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/xytoij.h
 *
 *  \brief xy2ij prototypes
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: xytoij.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _XYTOIJ_H
#define _XYTOIJ_H

#include "poly.h"
#include "hash.h"

typedef struct xytoij_leaf xytoij_leaf_t;

struct xytoij_leaf {
  poly_t *boundary;
  long mini;
  long maxi;
  long minj;
  long maxj;
  xytoij_leaf_t *half1;
  xytoij_leaf_t *half2;
};

typedef struct {
  long leaves;
  poly_t *outline;
  xytoij_leaf_t *trunk;

  double **gridx;               /* Array of X coords [nce2+1][nce2+1] */
  double **gridy;               /* Array of Y coords [nce2+1][nce2+1] */
  int nce1;                     /* Number of cells in e1 direction. */
  int nce2;                     /* Number of cells in e2 direction. */

/**** PRIVATE - used for housekeeping and efficiency reasons.
 ****/
  hash_table_t *ht;

} xytoij_tree_t;

/* Proto-types */
xytoij_tree_t *grid_xytoij_init(double **gx, double **gy, int nce1,
                                int nce2);
xytoij_tree_t *grid_xytoij_init_hash(double **gx, double **gy, int nce1,
                                     int nce2, int hsize);
int grid_xytoij(xytoij_tree_t *partition, double x, double y, int *ival,
                int *jval);
int grid_ijtoxy(xytoij_tree_t *partition, int ival, int jval, double *x,
                double *y);
int grid_xytofij(xytoij_tree_t *partition, double x, double y,
                 double *ival, double *jval);
int grid_fgrid_ijtoxy(xytoij_tree_t *partition, double ival, double jval,
                      double *x, double *y);
                      
/*UR-ADDED destroy funtions for tree and leaf */

void tree_destroy(xytoij_tree_t* tree);
void leaf_destroy(xytoij_leaf_t* leaf);

#endif
