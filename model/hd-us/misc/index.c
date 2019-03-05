/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/misc/index.c
 *  
 *  Description:
 *  Routines which return grid or variable
 *  indices for meco model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: index.c 5996 2018-10-17 04:57:10Z her127 $
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hd.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5
#define iszero(x) (((x) > -EPS_ZERO) && ((x) < EPS_ZERO))
#define EPS_ZEROP 1.0e-10
#define iszerop(x) (((x) > -EPS_ZEROP) && ((x) < EPS_ZEROP))

static int calc_branch_p(xytoij_tree_t *partition, double x, double y);

/*   x3,y3------x2,y2
 *     |          |
 *     |          |
 *     |          |
 *     |          |
 *     |          |
 *   x0,y0------x1,y1 
 *
 * The tranformation used to compute the indices is an inverse
 * tetragonal bilinear texture mapping. */

static int inside_cell(double x, double y,
                       double x0, double y0, double x1, double y1,
                       double x2, double y2, double x3, double y3)
{
  double a = x0 - x1 - x3 + x2;
  double b = x1 - x0;
  double c = x3 - x0;
  double d = x0;
  double e = y0 - y1 - y3 + y2;
  double f = y1 - y0;
  double g = y3 - y0;
  double h = y0;

  double A = a * f - b * e;
  double B = e * x - a * y + a * h - d * e + c * f - b * g;
  double C = g * x - c * y + c * h - d * g;

  double u, v, d1, d2;
  if (iszero(A))
    u = -C / B;
  else
    u = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);

  d1 = a * u + c;
  d2 = e * u + g;
  v = (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

  if ((u < 0.0) || (u >= 1.0))
    return 0;
  if ((v < 0.0) || (v >= 1.0))
    return 0;

  return 1;
}


static int inside_cell_p(double x, double y, double *u, double *v,
                       double x0, double y0, double x1, double y1,
                       double x2, double y2, double x3, double y3)
{
  double a = x0 - x1 - x3 + x2;
  double b = x1 - x0;
  double c = x3 - x0;
  double d = x0;
  double e = y0 - y1 - y3 + y2;
  double f = y1 - y0;
  double g = y3 - y0;
  double h = y0;

  double A = a * f - b * e;
  double B = e * x - a * y + a * h - d * e + c * f - b * g;
  double C = g * x - c * y + c * h - d * g;

  double d1, d2;
  if (iszerop(A))
    *u = -C / B;
  else
    *u = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);

  d1 = a * (*u) + c;
  d2 = e * (*u) + g;
  *v = (fabs(d2) > fabs(d1)) ? (y - f * (*u) - h) / d2 : (x - b * (*u) - d) / d1;

  if ((*u < 0.0) || (*u >= 1.0))
    return 0;
  if ((*v < 0.0) || (*v >= 1.0))
    return 0;
  return 1;
}


/* Find the wet cell that that this maps into.
 */
int hd_xyztoindex_w(void *data, double x, double y, double z,
                    int *cr, int *cs, int *cb)
{
  geometry_t *window = (geometry_t *)data;

  *cr = *cs = *cb = NOTVALID;
  *cr = hd_grid_xyztoc_w(window, x, y, z);
  if (*cr <= 0 || *cr >= window->szc) return(-1);

  *cs = window->m2d[*cr];
  *cb = window->bot_t[window->c2cc[*cs]];
  return(*cr);
}

int hd_xyztoindex_m(void *data, double x, double y, double z,
                    int *cr, int *cs, int *cb)
{
  master_t *master = (master_t *)data;
  geometry_t *geom = master->geom;

  *cr = *cs = *cb = NOTVALID;
  *cr = hd_grid_xyztoc_m(master, x, y, z);
  if (*cr <= 0 || *cr >= geom->szc) return(-1);

  *cs = geom->m2d[*cr];
  *cb = geom->bot_t[geom->c2cc[*cs]];
  return(*cr);
}

int hd_get_tracer_index_m(void *data, char *name)
{
  master_t *master = (master_t *)data;
  int i = 0;

  for (i = 0; i < master->ntr; ++i) {
    if (strcasecmp(name, master->trname[i]) == 0)
      return i;
  }

  return -1;
}

int hd_get_tracer_index_w(void *data, char *name)
{
  geometry_t *window = (geometry_t *)data;
  win_priv_t *wincon = window->wincon;
  int i = 0;

  for (i = 0; i < wincon->ntr; ++i) {
    if (strcasecmp(name, wincon->trname[i]) == 0)
      return i;
  }

  return -1;
}

long ztoc_m_old(master_t *master, int c, double z)
{
  if (c > 0 && c <= geom->sgnum) {
    double dc = 0.0;
    while (c != geom->zm1[c]) {
      dc += master->dz[c];
      c = geom->zm1[c];
      if (dc >= fabs(z))
        break;
    }
  }
  return (geom->zp1[c]);
}

long ztoc_m(master_t *master, int c, double z)
{
  if (c > 0 && c <= geom->sgnum) {
    while (geom->gridz[c] > z && c != geom->zm1[c]) {
      c = geom->zm1[c];
    }
  }
  return (c);
}

long ztoc_w(geometry_t *window, int c, double z)
{
  win_priv_t *wincon = window->wincon;
  if (c > 0 && c <= window->enon) {
    double dc = 0.0;
    while (c != window->zm1[c]) {
      dc += (wincon->dz[c] * wincon->Ds[window->m2d[c]]);
      c = window->zm1[c];
      if (dc >= fabs(z))
        break;
    }
  }
  return (window->zp1[c]);
}

long ztoc(geometry_t *geom, int c, double z)
{
  if (c > 0 && c <= geom->enon) {
    while (geom->gridz[c] > z && c != geom->zm1[c]) {
      c = geom->zm1[c];
    }
    return (c);
  }
  return(0);
}

long ztoc_dz(geometry_t *window, int c, double z, double *dz, double *fz)
{
  if (c > 0 && c <= window->enon) {
    double dc = 0.0;
    while (dc < fabs(z)) {
      dc += dz[c];
      c = window->zm1[c];
    }
    *fz = dc + z;
    return (window->zp1[c]);
  }
  return(0);
}

long ztoc_z(geometry_t *window, int c, double z, double *cz)
{
  if (c > 0 && c <= window->enon) {
    while (cz[c] > z && c != window->zm1[c]) {
      c = window->zm1[c];
    }
    return(c);
  }
  return(0);
}

int xyztoindex_m(master_t *master, double x, double y, double z,
                 int *cr, int *cs, int *cb)
{
  geometry_t *geom = master->geom;
  double i, j;

  if(grid_xytofij(master->xyij_tree, x, y, &i, &j)) {
    *cs = geom->map[geom->nz-1][(int)j][(int)i];
    *cb = geom->bot_t[geom->c2cc[*cs]];

    /*if(ANY(c, geom->wsa, geom->a2_t)) {*/
    if(*cs > 0 && *cs < geom->sgnumS) {
      *cr = ztoc_m(master, *cs, z);
      return(*cr);
    }
    else
      return(-1);
  }
  else
    return(-1);
}


int xyztoindex(geometry_t *geom, xytoij_tree_t *xyij_tree, double *x, double *y, double z,
	       int *cs, int *cb)
{
  double d1, d2;
  int cr;

  if(grid_xytofij(xyij_tree, *x, *y, &d1, &d2)) {
    *cs = geom->map[geom->nz-1][(int)d1][(int)d2];
    *cb = geom->bot_t[geom->c2cc[*cs]];
    *x = d1;
    *y = d2;
    if(*cs > 0 && *cs < geom->sgnumS) {
      cr = *cs;
      while (geom->gridz[cr] > z && cr != geom->zm1[cr]) {
	cr = geom->zm1[cr];
      }
      return(cr);
    }
    else
      return(0);
  }
  else
    return(0);
}


/* Precision computation of fractional (i,j) given (x,y) */
/* (i.e. uses iszsero() to a precision of 1e-15).        */
/* Uses an inverse tetragonal bilinear texture mapping;  */
/* e.g. http://www.cescg.org/CESCG97/olearnik/txmap.htm  */
int xytoc(geometry_t *geom, xytoij_tree_t *xyij_tree, double x, double y, double *w1, double *w2)
{
  int i, j, c;
  double **gx = xyij_tree->gridx;
  double **gy = xyij_tree->gridy;
  static int sign = 0;

  if (grid_xytoij(xyij_tree, x, y, &i, &j) == 0)
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
    if (iszerop(A))
      u = -C / B;
    else {
      if (sign == 0) {
        sign = calc_branch_p(xyij_tree, x, y);
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
    else if (u >= 1.0 - EPS_ZEROP) {
      u = EPS_ZEROP;
      i += 1;
    }
    if (v < 0.0)
      v = 0.0;
    else if (v >= 1.0 - EPS_ZEROP) {
      v = EPS_ZEROP;
      j += 1;
    }
    *w1 = i + u;
    *w2 = j + v;
  }
  c = geom->map[geom->nz-1][j][i];
  if (i < 0 || i > geom->nce1) return(0);
  if (j < 0 || j > geom->nce2) return(0);
  if (c > 0 && c <= geom->sgnumS)
    return (c);
  else
    return(0);
}


/** Calculates the branch of sqrt() to be taken in grid_xytofij(). Has to be
   * called only once for a grid.
   *
   * @param partition a xytoij_tree_t structure returned from xytoij_init.
   * @param x Pointer to returned X coordinate.
   * @param y Pointer to returned Y coordinate.
   * @return 1 or -1 if successful; 0 otherwhile.
   */
static int calc_branch_p(xytoij_tree_t *partition, double x, double y)
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
    if (iszerop(A))
      return 0;                 /* failed */

    B = e * x - a * y + a * h - d * e + c * f - b * g;
    C = g * x - c * y + c * h - d * g;

    for (k = 0; k < 2; ++k) {
      error[k] = 0.0;

      u = (-B + sign * sqrt(B * B - 4 * A * C)) / (2 * A);
      v_denom = a * u + c;
      v =
        (iszerop(v_denom)) ? (y - f * u - h) / (e * u + g) : (x - b * u -
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

/* Wrapper to create a GRID_SPEC interpolation structure */
void hd_grid_interp_init(GRID_SPECS *gs, double *tr, char *method)
{
  double *x = geom->cellx;
  double *y = geom->celly;

  gs = grid_interp_init(x, y, tr, geom->enonS, method);
}


int hd_grid_interp_3d(geometry_t *geom, double *ret)
{
  int c, cc, cs;
  double x, y, z;
  if (geom->gs == NULL) return(0);

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    x = geom->cellx[cs];
    y = geom->celly[cs];
    z = geom->cellz[c] * master->Ds[cs];
    ret[c] = grid_interp_on_point(geom->gs[0], x, y);
  }
  return(1);
}

int hd_grid_interp_2d(GRID_SPECS *gs, double *ret)
{
  int c, cc;
  double x, y, z;
  if (geom->gs == NULL) return(0);

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    x = geom->cellx[c];
    y = geom->celly[c];
    ret[c] = grid_interp_on_point(gs, x, y);
  }
  return(1);
}

/*-------------------------------------------------------------------*/
/* Wrapper to find an (i,j) given (x,y) for structured grids, or the */
/* unstructured index for unstructured grids. Always returns the     */
/* unstructured (sparse) coordinate. Also returns i and j for        */
/* structured grids.                                                 */
/* Returns zero if unsuccessful.                                     */
/*-------------------------------------------------------------------*/
int hd_grid_xytoij(master_t *master, double x, double y, int *i, int *j) 
{
  geometry_t *geom = master->geom;
  int c, ret = 0;
  point p;

  if (geom->us_type & US_QUAD && geom->us_type & US_IJ) {
    if (grid_xytoij(master->xyij_tree, x, y, i, j))
      ret = geom->map[geom->nz - 1][*j][*i];
  } else if (geom->us_type & (US_TRI|US_HEX|US_MIX)) {
    delaunay *d = geom->d;
    p.x = x;
    p.y = y;
    c = delaunay_xytoi_lag(d, &p, d->first_id);
    if (c != -1)
      ret = *i = geom->tri2c[c];
  }
  return(ret);
}

/* Same as the above, but returns -1 if location can't be found */
int hd_grid_xytoij_r(master_t *master, double x, double y, int *i, int *j) 
{
  geometry_t *geom = master->geom;
  int c, ret = 0;
  point p;

  if (geom->us_type & US_QUAD && geom->us_type & US_IJ) {
    if (grid_xytoij(master->xyij_tree, x, y, i, j))
      ret = geom->map[geom->nz - 1][*j][*i];
  } else if (geom->us_type & (US_TRI|US_HEX|US_MIX)) {
    delaunay *d = geom->d;
    p.x = x;
    p.y = y;
    c = delaunay_xytoi(d, &p, d->first_id);
    if (c != -1)
      ret = *i = geom->tri2c[c];
  }
  return(ret);
}


/*-------------------------------------------------------------------*/
/* Wrapper to find the grid index given (x,y,z).                     */
/* Returns zero if unsuccessful.                                     */
/* Assumes set_dzz() has been called so wincon->cellz is populated   */
/* the the vertical cellz spacing (accounting for the free surface   */
/* and partial bottom).                                              */
/*-------------------------------------------------------------------*/
int hd_grid_xyztoc(geometry_t *window, double x, double y, double z, int ci)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  double *cellz = wincon->cellz;
  int c, cs;

  /* Find the horizontal position in the mesh                        */
  if (window->us_type & US_QUAD && window->us_type & US_IJ) {
    int i = -1, j = -1;
    if (grid_xytoij(window->xyij_tree, x, y, &i, &j)) {
      cs = window->map[geom->nz - 1][j][i];
    } else
      cs = window->m2d[ci];
    if (!cs) cs = window->m2d[ci];
  } else if (window->us_type & (US_TRI|US_HEX|US_MIX)) {
    delaunay *d = window->d;
    point p;
    p.x = x;
    p.y = y;
    cs = delaunay_xytoi_lag(d, &p, d->first_id);
    if (cs != -1)
      cs = geom->tri2c[cs];
    else
      cs = window->m2d[ci];
  }

  /* Find the vertical position in the mesh                          */
  c = cs;
  /* First get the layer of the free surface                         */
  while (windat->eta[cs] < window->gridz[c] && c != window->zm1[c]) {
    c = window->zm1[c];
  }
  /* Get the first cellz layer less than z                           */
  while (z < cellz[c] && c != window->zm1[c]) {
    c = window->zm1[c];
  }
  return(c);
}

int hd_grid_xyztoc_w(geometry_t *window, double x, double y, double z)
{
  window_t *windat = window->windat;
  int cc, c, cg, cs;

  /* Find the horizontal position in the mesh                        */
  if (window->us_type & US_QUAD && window->us_type & US_IJ) {
    int i = -1, j = -1;
    if (grid_xytoij(window->xyij_tree, x, y, &i, &j)) {
      int found = 0;
      cg = window->map[window->nz - 1][j][i]; /* Global coordinate */
      for (cc = 1; cc <= window->b2_t; cc++) {
	cs = window->w2_t[cc];
	c = window->wsa[cs];
	if (c == cg) {
	  found = 1;
	  break;
	}
      }
      if (!found) return(-1); /* Not in the window */
    } else
      return(-1);             /* Not in the grid */
  } else if (window->us_type & (US_TRI|US_HEX|US_MIX)) {
    delaunay *d = window->d;
    point p;
    int found = 0;
    p.x = x;
    p.y = y;

    cg = delaunay_xytoi_lag(d, &p, d->first_id);
    if (cg != -1) {
      cg = window->tri2c[cg];

      for (cc = 1; cc <= window->b2_t; cc++) {
	cs = window->w2_t[cc];
	c = window->wsa[cs];
	if (c == cg) {
	  found = 1;
	  break;
	}
      }
      if (!found) return(-1);
    } else
      return(-1);
  }

  /* Find the vertical position in the mesh                          */
  c = cs;
  /* First get the layer of the free surface                         */
  while (z < window->gridz[c] && c != window->zm1[c]) {
    c = window->zm1[c];
  }
  return(c);
}

int hd_grid_xyztoc_m(master_t *master, double x, double y, double z)
{
  geometry_t *geom = master->geom;
  double *cellz = geom->cellz;        /* cellz depth                 */
  int c, cs;

  /* Find the horizontal position in the mesh                        */
  if (geom->us_type & US_QUAD && geom->us_type & US_IJ) {
    int i = -1, j = -1;
    if (grid_xytoij(master->xyij_tree, x, y, &i, &j)) {
      cs = geom->map[geom->nz - 1][j][i];
    } else
      return(-1);
  } else if (geom->us_type & (US_TRI|US_HEX|US_MIX)) {
    delaunay *d = geom->d;
    point p;
    p.x = x;
    p.y = y;
    cs = delaunay_xytoi_lag(d, &p, d->first_id);
    if (cs != -1)
      cs = geom->tri2c[cs];
    else
      return(-1);
  }

  /* Find the vertical position in the mesh                          */
  c = cs;
  /* First get the layer of the free surface                         */
  while (z < geom->gridz[c] && c != geom->zm1[c]) {
    c = geom->zm1[c];
  }
  return(c);
}
