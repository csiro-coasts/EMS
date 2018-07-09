/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/grid/xytoij.c
 *
 *  \brief Tree based (x,y) to (i,j) conversions
 *
 *  Calculates the indices (i,j) of a topologically rectangular grid 
 *  cell containing the point (x,y)
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: xytoij.c 5833 2018-06-27 00:21:35Z riz008 $
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>
#include "ems.h"


#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5
#define iszero(x) (((x) > -EPS_ZERO) && ((x) < EPS_ZERO))

void xytoij_hash_destroy(void* k, void* d);

typedef struct {
  int i;
  int j;
} xyij_hash_data_t;

int hash_opt = 1;

xytoij_leaf_t *create_leaf(poly_t *pl, long int i1, long int i2,
                           long int j1, long int j2, double **gx,
                           double **gy)
{
  int n = pl->n;
  xytoij_leaf_t *petal;
  double x, y;
  int i = 0, j = 0, ii = 0;

  if ((petal = (xytoij_leaf_t *)malloc(sizeof(xytoij_leaf_t))) == NULL)
    quit("create_leaf: memory allocation failure\n");
  memset(petal, 0, sizeof(xytoij_leaf_t));

  petal->boundary = pl;
  petal->mini = INT_MAX;
  petal->maxi = INT_MIN;
  petal->minj = INT_MAX;
  petal->maxj = INT_MIN;
  petal->half1 = NULL;
  petal->half2 = NULL;

  if (n == 0)
    return petal;

  x = pl->x[0];
  y = pl->y[0];

  for (j = j1; j <= j2; ++j) {
    for (i = i1; i <= i2; ++i) {
      if (x == gx[j][i] && y == gy[j][i])
        goto aftersearch;
    }
  }

aftersearch:

  if (j > j2)
    quit("create_leaf(): boundary vertex not in the grid, j = %i, j2 = %i\n", j, j2);

  petal->mini = i;
  petal->maxi = i;
  petal->minj = j;
  petal->maxj = j;

  for (ii = 1; ii < n; ++ii) {
    x = pl->x[ii];
    y = pl->y[ii];

    if (i > i1 && x == gx[j][i - 1] && y == gy[j][i - 1]) {
      emslog(LDEBUG, "xytoij.c aftersearch: ii = %i, decrementing i, x = %6.3f, gx[j][i+1] = %6.3f, y = %6.3f, gy[j][i-1] = %6.3f\n", ii, x, gx[j][i-1], y, gy[j][i-1]);
      /* Decrement after we've logged */
      i--;
    } else if (i < i2 && x == gx[j][i + 1] && y == gy[j][i + 1])
      i++;
    else if (j > j1 && x == gx[j - 1][i] && y == gy[j - 1][i])
      j--;
    else if (j < j2 && x == gx[j + 1][i] && y == gy[j + 1][i])
      j++;
    else if (x == gx[j][i] && y == gy[j][i])
      continue;
    else
      quit("create_leaf(): boundary vertex not in the grid, i = %i, i1 = %i, i2 = %i, j = %i, j1 = %i, j2 = %i, n = %i, x = %6.3f, y = %6.3f, gx = %6.3f, gy = %6.2f\n", i, i1, i2, j, j1, j2, n, x, y, gx[j][i], gy[j][i]);

    if (petal->mini > i)
      petal->mini = i;
    if (petal->maxi < i)
      petal->maxi = i;
    if (petal->minj > j)
      petal->minj = j;
    if (petal->maxj < j)
      petal->maxj = j;
  }

  return petal;
}

xytoij_tree_t *create_tree(xytoij_leaf_t *startleaf, poly_t *outline)
{
  xytoij_tree_t *conifer;

  if ((conifer = (xytoij_tree_t *)malloc(sizeof(xytoij_tree_t))) == NULL)
    quit("create_tree: memory allocation failure\n");
  memset(conifer, 0, sizeof(xytoij_tree_t));
  conifer->leaves = 1;
  conifer->outline = outline;
  conifer->trunk = startleaf;

  return (conifer);
}

void add_to_leaf(xytoij_leaf_t *currleaf, xytoij_leaf_t *leaf1,
                 xytoij_leaf_t *leaf2)
{
  currleaf->half1 = leaf1;
  currleaf->half2 = leaf2;
}

void add_to_tree(xytoij_tree_t *conifer, xytoij_leaf_t *currleaf,
                 xytoij_leaf_t *leaf1, xytoij_leaf_t *leaf2)
{
  add_to_leaf(currleaf, leaf1, leaf2);
  if (leaf1 != NULL)
    ++(conifer->leaves);
  if (leaf2 != NULL)
    ++(conifer->leaves);
}

int find_average(int n1, int n2)
{
  return (n1 + n2) / 2;
}

poly_t *form_boundary(int nce1, int nce2, double **x, double **y)
{
  poly_t *pl = poly_create();
  int i, j, istart, jstart;
  enum { DOWN = 0, RIGHT = 1, UP = 2, LEFT = 3 } direction;
  int iinc[] = { 0, 1, 0, -1 };
  int jinc[] = { 1, 0, -1, 0 };

  for (j = 0; j <= nce2; ++j)
    for (i = 0; i <= nce1; ++i)
      if (!isnan(x[j][i]))
        goto ok;

ok:
  if (j > nce2)
    return pl;

  istart = i;
  jstart = j;
  direction = DOWN;
  poly_add_point(pl, x[j][i], y[j][i]);

  do {
    int direction_stop = (direction + 2) % 4;
    int inext, jnext;

    direction = (direction + 3) % 4;
    inext = i + iinc[direction];
    jnext = j + jinc[direction];

    while ((inext < 0 || jnext < 0 || inext > nce1 || jnext > nce2 ||
            isnan(x[jnext][inext])) &&
           ((int)direction) != direction_stop) {
      direction = (direction + 1) % 4;
      inext = i + iinc[direction];
      jnext = j + jinc[direction];
    }

    if (((int)direction) == direction_stop)
      break;

    i = inext;
    j = jnext;

    poly_add_point(pl, x[j][i], y[j][i]);

  } while (j != jstart || i != istart);

  poly_close(pl);

  return (pl);
}

/** Cuts polygon in two by polyline.
 * The polyline goes either horizontally ([fixed][changes]) or vertically
 * ([changes][fixed]) in index space; the physical nodes are given by
 * input double arrays; first two intersections of the cutting polyline
 * with the polygon are used to form the new polygons.
 * @param pl Original polygon.
 * @param x Array of x node coordinates.
 * @param y Array of y node coordinates.
 * @param horiz flag: 1 for horizontal cut; 0 otherwise.
 * @param index Value of "fixed" index.
 * @param start Start value of "variable" index.
 * @param end End value of "variable" index.
 * @param pl1 Output polygon 1.
 * @param pl2 Output polygon 2.
 */
void cut_boundary(poly_t *pl, double **x, double **y, int horiz, int index,
                  int start, int end, poly_t **pl1, poly_t **pl2)
{
  int n = pl->n;
  int i = -1;
  int i1 = -1;                  /* array index of the first intersection */
  int i2 = -1;                  /* array index of the second intersection */
  int ii1 = -1;                 /* polygon index of the first intersection
                                 */
  int ii2 = -1;                 /* polygon index of the second
                                   intersection */
  int tmp = -1;

  /* if the polygon has been explicitely closed, ignore the last point */
  if (poly_is_closed(pl, 1.0e-8))
    n--;

  if (horiz) {
    /* find first intersection */
    for (i = start; i < end; ++i) {
      ii1 = poly_find_index(pl, x[index][i], y[index][i]);
      if (ii1 < 0)
        /* this node does not belong to the boundary */
        continue;
      /* this node belongs to the boundary polygon To accept it as a valid
         intersection, we must ensure that the very next node is either
         inside the polygon or belongs to the boundary but is not the next
         boundary node. */
      tmp = poly_find_index(pl, x[index][i + 1], y[index][i + 1]);
      if ((tmp < 0 &&
           poly_contains_point(pl, x[index][i + 1], y[index][i + 1])) ||
          (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
        break;
    }

    if (i < end)
      /* how the for-cycle ended */
      /* ok */
      i1 = i;
    else
      /* no intersection found */
      return;

    /* find second intersection start from the node next to the first
       intersection */
    for (i = i1 + 1; i <= end; ++i) {
      ii2 = poly_find_index(pl, x[index][i], y[index][i]);
      if (ii2 >= 0)
        /* this node must be inside the boundary polygon -- skip */
        break;
    }

    if (ii2 < 0)
      /* no intersection found */
      return;
    else
      /* ok */
      i2 = i;

    /* we found all necessary details, now form the new polygons */
    *pl1 = poly_create();
    *pl2 = poly_create();

    /* add the portion of perimeter */
    for (i = ii1; i != ii2; i = (i + 1) % n)
      poly_add_point(*pl1, pl->x[i], pl->y[i]);
    /* add the cutting section */
    for (i = i2; i > i1; --i)
      poly_add_point(*pl1, x[index][i], y[index][i]);

    /* add the portion of perimeter */
    for (i = ii2; i != ii1; i = (i + 1) % n)
      poly_add_point(*pl2, pl->x[i], pl->y[i]);
    /* add the cutting section */
    for (i = i1; i < i2; ++i)
      poly_add_point(*pl2, x[index][i], y[index][i]);

  } else {                      /* vertical cut */
    for (i = start; i < end; ++i) {
      ii1 = poly_find_index(pl, x[i][index], y[i][index]);
      if (ii1 < 0)
        continue;
      tmp = poly_find_index(pl, x[i + 1][index], y[i + 1][index]);
      if ((tmp < 0 &&
           poly_contains_point(pl, x[i + 1][index], y[i + 1][index])) ||
          (tmp >= 0 && abs(tmp - ii1) > 1 && abs(tmp - ii1) < n - 1))
        break;
    }

    if (i < end)
      i1 = i;
    else
      return;

    for (i = i1 + 1; i <= end; ++i) {
      ii2 = poly_find_index(pl, x[i][index], y[i][index]);
      if (ii2 >= 0)
        break;
    }

    if (ii2 < 0)
      return;
    else
      i2 = i;

    *pl1 = poly_create();
    *pl2 = poly_create();

    for (i = ii1; i != ii2; i = (i + 1) % n)
      poly_add_point(*pl1, pl->x[i], pl->y[i]);
    for (i = i2; i > i1; --i)
      poly_add_point(*pl1, x[i][index], y[i][index]);

    for (i = ii2; i != ii1; i = (i + 1) % n)
      poly_add_point(*pl2, pl->x[i], pl->y[i]);
    for (i = i1; i < i2; ++i)
      poly_add_point(*pl2, x[i][index], y[i][index]);

    // close the polygons
    poly_close(*pl1);
    poly_close(*pl2);
  }
}

void divide_leaf(xytoij_leaf_t *currleaf, xytoij_leaf_t **leaf1,
                 xytoij_leaf_t **leaf2, double **gx, double **gy)
{
  poly_t *pl1 = NULL;
  poly_t *pl2 = NULL;
  int index;

  if ((currleaf->maxi <= currleaf->mini + 1) &&
      (currleaf->maxj <= currleaf->minj + 1)) {
    *leaf1 = *leaf2 = NULL;
    return;
  }

  if (currleaf->maxi - currleaf->mini > currleaf->maxj - currleaf->minj) {
    /* divide "vertically" */
    index = find_average(currleaf->mini, currleaf->maxi);
    cut_boundary(currleaf->boundary, gx, gy, 0, index, currleaf->minj,
                 currleaf->maxj, &pl1, &pl2);
  } else {
    index = find_average(currleaf->minj, currleaf->maxj);
    cut_boundary(currleaf->boundary, gx, gy, 1, index, currleaf->mini,
                 currleaf->maxi, &pl1, &pl2);
  }

  if (pl1 == NULL || pl2 == NULL)
    quit("divide_leaf(): could not cut the boundary\n");

  *leaf1 =
    create_leaf(pl1, currleaf->mini, currleaf->maxi, currleaf->minj,
                currleaf->maxj, gx, gy);
  *leaf2 =
    create_leaf(pl2, currleaf->mini, currleaf->maxi, currleaf->minj,
                currleaf->maxj, gx, gy);
}

void sub_divide(xytoij_tree_t *partition, xytoij_leaf_t *currleaf,
                double **gx, double **gy)
{
  xytoij_leaf_t *leaf1;
  xytoij_leaf_t *leaf2;

  divide_leaf(currleaf, &leaf1, &leaf2, gx, gy);
  add_to_tree(partition, currleaf, leaf1, leaf2);
  if (leaf1 != NULL)
    sub_divide(partition, leaf1, gx, gy);
  if (leaf2 != NULL)
    sub_divide(partition, leaf2, gx, gy);
}

int xyij_compare(void *o1, void *o2)
{
  /*
  point_t *k1 = (point_t *)o1;
  point_t *k2 = (point_t *)o2;

  return (fabs(k1->x - k2->x) < EPS) && (fabs(k1->y - k2->y) < EPS);
  */
  double* k1 = (double*)o1;
  double* k2 = (double*)o2;

  return (fabs(k1[0] - k2[0]) < EPS) && (fabs(k1[1] - k2[1]) < EPS);
}

size_t xyij_hash(void *o)
{
/*UR-FIX
 * below will cause serious problems on 64-bit machines
 * consequence of using a long will be that lpy * 3 will be
 * placed outside the available space ad thus produce rubbish
 * the hashtable becomes a memory leak...!
 *   point_t *p = (point_t *)o;
  unsigned long *lpx = (unsigned long *)&p->x;
  unsigned long *lpy = (unsigned long *)&p->y;
  unsigned long bits = 0L;
  */
  double* k = (double*)o;
  unsigned int *lpx = (unsigned int *)&k[0];
  unsigned int *lpy = (unsigned int *)&k[1];
  unsigned int bits = 0;
  bits = *lpx++;
  bits ^= *lpy++ * 3;
  bits ^= *lpx * 31;
  bits ^= *lpy * 7;
  return (size_t)((bits) ^ (bits >> 31));
}

/** Initialises a xytoij_tree_t structure, given a pair of 2 dimensional
  * regular coordinates, to facilitate conversion from coordinate
  * to inidice space.
  *
  * @param gx array of X coordinates (of size (nce1+1)*(nce2+1)).
  * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1)).
  * @param nce1 number of cells in e1 direction.
  * @param nce2 number of cells in e2 direction.
  * @param hsize hash table size.
  * @return a partition xytoij_tree_t to be use by xytoij.
  */
xytoij_tree_t *grid_xytoij_init_hash(double **gx, double **gy, int nce1,
                                     int nce2, int hsize)
{
  xytoij_tree_t *partition = grid_xytoij_init(gx, gy, nce1, nce2);
  emstag(LTRACE,"lib:xytoij:grid_xytoij_init_hash","Creating tree with partition of size %d ",hsize);

  if (hsize > 0)
    partition->ht = ht_create_complex(hsize, xyij_hash, xyij_compare,
				      NULL, xytoij_hash_destroy);

  return partition;
}


/** Initialises a xytoij_tree_t structure, given a pair of 2 dimensional
  * regular coordinates, to facilitate conversion from coordinate
  * to inidice space.
  *
  * @param gx array of X coordinates (of size (nce1+1)*(nce2+1)).
  * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1)).
  * @param nce1 number of cells in e1 direction.
  * @param nce2 number of cells in e2 direction.
  * @return a partition xytoij_tree_t to be use by xytoij.
  */
xytoij_tree_t *grid_xytoij_init(double **gx, double **gy, int nce1,
                                int nce2)
{
  xytoij_tree_t *partition;

  xytoij_leaf_t *firstleaf;
  poly_t *firstpl;

  firstpl = form_boundary(nce1, nce2, gx, gy);
  firstleaf = create_leaf(firstpl, 0, nce1, 0, nce2, gx, gy);
  partition = create_tree(firstleaf, firstpl);
  sub_divide(partition, firstleaf, gx, gy);

  /* The values are not copies only referenced */
  partition->gridx = gx;
  partition->gridy = gy;
  partition->nce1 = nce1;
  partition->nce2 = nce2;
  partition->ht = NULL;         /* No hash by default */

  return partition;
}


/** calculates the indices (i,j) of a topologically rectangular
  * grid cell containing the point (x,y).
  *
  * @param partition a xytoij_tree_t structure returned from xytoij_init.
  * @param x X coordinate.
  * @param y Y coordinate.
  * @param ival pointer to returned I indice value.
  * @param jval pointer to returned J indice value.
  * @return non-zero if successful.
  */
int grid_xytoij(xytoij_tree_t *partition, double x, double y, int *ival,
                int *jval)
{
  xytoij_leaf_t *currleaf;

  /* Check to see that initialisation has been done */
  if (partition == NULL)
    quit("grid_xytoij: initialisation of partition not done\n");

  if (ival == NULL || jval == NULL)
    return (0);

  /* check the hash table to see if this is already in the hash table */
  if (partition->ht != NULL) {
      double xy[2] = { x, y };
    xyij_hash_data_t *hd = (xyij_hash_data_t *)ht_find(partition->ht, xy);
    if (hd != NULL) {
      *ival = hd->i;
      *jval = hd->j;
      return 1;
    }
  }

  /* check if point is in grid outline */
  currleaf = partition->trunk;
  if (!poly_contains_point(currleaf->boundary, x, y))
    return 0;

  /* do the full search */
  while (currleaf->half1 != NULL) {
    if (poly_contains_point(currleaf->half1->boundary, x, y))
      currleaf = currleaf->half1;
    else
      currleaf = currleaf->half2;
  }

  *ival = currleaf->mini;
  *jval = currleaf->minj;

  /* UR-ADDED 16/3/2005
   * disable the test after it failed once, otherwise the log  causes
   * a substantial performance hit
   */
  if( !hash_opt)
    return (1);

  /* Add to hash table if not too full (80 %) */
  if (partition->ht != NULL) {
    if (((partition->ht->noccupied * 100) / partition->ht->nelems) < 80) {
      xyij_hash_data_t *hd;

      double *xy = malloc(2 * sizeof(double));
      xy[0] = x;
      xy[1] = y;
      hd = (xyij_hash_data_t *)malloc(sizeof(xyij_hash_data_t));
      hd->i = *ival;
      hd->j = *jval;

      ht_add(partition->ht, xy, hd);

      /* if ((partition->ht->noccupied % 10000) == 0)
         ht_print_stats(partition->ht); */
    } else
    {
      warn
        ("grid_xytoij: hash table optimisation has been disabled - table full.");
      hash_opt =0;
    }
  }

  return (1);
}

/** Calculates the XY coordinate for the specified fractional (i,j)
  * within a topologically rectangular grid.
  *
  * The tranformation used to compute the coords is a forward
  * tetragonal bilinear texture mapping.
  *
  * @param partition a xytoij_tree_t structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int grid_fgrid_ijtoxy(xytoij_tree_t *partition, double ival, double jval,
                      double *x, double *y)
{
  int status = 1;
  int i, j;
  double u, v;
  double **gx = partition->gridx;
  double **gy = partition->gridy;
  double a, b, c, d, e, f, g, h;

  /* Trim I to range 0 to nce1 */
  if (ival < 0) {
    ival = 0;
    status = 0;
  }

  if (ival > partition->nce1) {
    ival = partition->nce1 - EPS;
    status = 0;
  }

  /* Trim J to range 0 to nce2 */
  if (jval < 0) {
    jval = 0;
    status = 0;
  }

  if (jval > partition->nce2) {
    jval = partition->nce2 - EPS;
    status = 0;
  }

  i = (int)ival;
  j = (int)jval;
  u = ival - i;
  v = jval - j;

  a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
  b = gx[j][i + 1] - gx[j][i];
  c = gx[j + 1][i] - gx[j][i];
  d = gx[j][i];
  e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
  f = gy[j][i + 1] - gy[j][i];
  g = gy[j + 1][i] - gy[j][i];
  h = gy[j][i];

  *x = a * u * v + b * u + c * v + d;
  *y = e * u * v + f * u + g * v + h;

  return status;
}


/** Calculates the XY coordinate for the specified (i,j) within a
  * topologically rectangular grid.
  *
  * @param partition a xytoij_tree_t structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int grid_ijtoxy(xytoij_tree_t *partition, int ival, int jval, double *x,
                double *y)
{
  int status = 1;

  /* Trim I to range 0 to nce1 */
  if (ival < 0) {
    ival = 0;
    status = 0;
  }

  if (ival > partition->nce1) {
    ival = partition->nce1;
    status = 0;
  }

  /* Trim J to range 0 to nce2 */
  if (jval < 0) {
    jval = 0;
    status = 0;
  }

  if (jval > partition->nce2) {
    jval = partition->nce2;
    status = 0;
  }

  *x = partition->gridx[jval][ival];
  *y = partition->gridy[jval][ival];

  return status;
}

/** Calculates the branch of sqrt() to be taken in grid_xytofij(). Has to be
   * called only once for a grid.
   *
   * @param partition a xytoij_tree_t structure returned from xytoij_init.
   * @param x Pointer to returned X coordinate.
   * @param y Pointer to returned Y coordinate.
   * @return 1 or -1 if successful; 0 otherwhile.
   */
static int calc_branch(xytoij_tree_t *partition, double x, double y)
{
  int i, j;
  double **gx = partition->gridx;
  double **gy = partition->gridy;
  int sign = 1;
  double error[2];

  /* normally one tries grid_xytoij() before calling calc_branch() */
  if (grid_xytoij(partition, x, y, &i, &j) == 0)
    return 0;                   /* failed */

  {
    double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
    double b = gx[j][i + 1] - gx[j][i];
    double c = gx[j + 1][i] - gx[j][i];
    double d = gx[j][i];
    double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
    double f = gy[j][i + 1] - gy[j][i];
    double g = gy[j + 1][i] - gy[j][i];
    double h = gy[j][i];

    double A = a * f - b * e;

    double B, C;
    int k;
    double u, v, v_denom;

    /* normally one checks A before calling calc_branch() */
    if (iszero(A))
      return 0;                 /* failed */

    B = e * x - a * y + a * h - d * e + c * f - b * g;
    C = g * x - c * y + c * h - d * g;

    for (k = 0; k < 2; ++k) {
      error[k] = 0.0;

      u = (-B + sign * sqrt(B * B - 4 * A * C)) / (2 * A);
      v_denom = a * u + c;
      v =
        (iszero(v_denom)) ? (y - f * u - h) / (e * u + g) : (x - b * u -
                                                             d) / v_denom;
      if (u < 0.0)
        error[k] -= u;
      else if (u > 1.0)
        error[k] += (u - 1.0);
      if (v < 0.0)
        error[k] -= v;
      else if (v > 1.0)
        error[k] += (v - 1.0);

      sign = -1;
    }
  }

  if (error[0] < error[1])
    return 1;
  return -1;
}

/** Calculates the XY coordinate for the specified fractional (i,j)
  * within a topologically rectangular grid.
  *
  * The tranformation used to compute the indices is an inverse
  * tetragonal bilinear texture mapping.
  *
  * At the moment we assume that either there is only one grid in the model or
  * all the grids use a uniform branch in xitofij() convertion.
  *
  * @param partition a xytoij_tree_t structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int grid_xytofij(xytoij_tree_t *partition, double x, double y,
                 double *ival, double *jval)
{
  int i, j;
  double **gx = partition->gridx;
  double **gy = partition->gridy;
  static int sign = 0;      /* SHOULD THIS BE STATIC ? JRW */

  if (grid_xytoij(partition, x, y, &i, &j) == 0)
    return 0;                   /* failed */

  {
    double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
    double b = gx[j][i + 1] - gx[j][i];
    double c = gx[j + 1][i] - gx[j][i];
    double d = gx[j][i];
    double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
    double f = gy[j][i + 1] - gy[j][i];
    double g = gy[j + 1][i] - gy[j][i];
    double h = gy[j][i];

    double A = a * f - b * e;
    double B = e * x - a * y + a * h - d * e + c * f - b * g;
    double C = g * x - c * y + c * h - d * g;

    double u, v, d1, d2;
    if (iszero(A))
      u = -C / B;
    else {
      if (sign == 0) {
        sign = calc_branch(partition, x, y);
        if (sign == 0)
          return 0;             /* failed */
      }
      u = (-B + sign * sqrt(B * B - 4 * A * C)) / (2 * A);
    }
    d1 = a * u + c;
    d2 = e * u + g;
    v =
      (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

    if (u < 0.0)
      u = 0.0;
    else if (u >= 1.0)
      u = 1.0 - EPS;
    if (v < 0.0)
      v = 0.0;
    else if (v >= 1.0)
      v = 1.0 - EPS;

    *ival = i + u;
    *jval = j + v;
  }

  return 1;
}


/*UR-ADDED adding destroy functions
 */
void leaf_destroy(xytoij_leaf_t* leave)
{
  if (leave->half1 != NULL)
    leaf_destroy(leave->half1);
  if (leave->half2 != NULL)
    leaf_destroy(leave->half2);
  poly_destroy(leave->boundary);
  free(leave);
}


void tree_destroy(xytoij_tree_t* tree)
{
  leaf_destroy(tree->trunk);

  if(tree->ht != NULL) {
    emstag(LTRACE,"lib:xytoij:grid_xytoij_init_hash","Destroy tree with partition of size %d occupied %d",tree->ht->nelems,tree->ht->noccupied);
    ht_destroy(tree->ht);
  }

  /*UR-TODO fix this - there is memory assigned here which is hard to track
   * - dmalloc fails * /
  if(tree->outline != NULL)
    poly_destroy(tree->outline);
  fprintf(stderr,"finished tree_destroy:outline \n");

  */
  /*UR-TODO are these destroyed anywhere ?
  if(tree->gridx != NULL)
    free_2d(tree->gridx);
  fprintf(stderr,"finished tree_destroy:gridx \n");

  if(tree->gridy != NULL)
    free_2d(tree->gridy);
  fprintf(stderr,"finished tree_destroy:gridy \n");
    */
  free(tree);
}


/*UR-ADDED ffree the key and data in an appropriate way */
void xytoij_hash_destroy(void* k, void* d)
{
  free(k);
  free(d);
}


