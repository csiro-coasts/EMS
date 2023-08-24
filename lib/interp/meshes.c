/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/interp/meshes.c
 *
 *  \brief Unstructured mesh construction routines
 *
 *  Routines which set up traingulationa nd corresponding dual indexing.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lagrange.c 6876 2021-07-29 00:40:40Z her127 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "ems.h"

int edgesort_double (void const *x1, void const *x2);
int edgesort_compare (void const *x1, void const *x2);
void centre_finder(meshs_t *m, delaunay *d, int mode);
point *concave_hull(delaunay *d, point *v, int *count, int *mask, int **map);
static int compare_int(const void* v1, const void* v2);

/*-------------------------------------------------------------------*/
/* Initialize a mesh structure.                                      */
/*-------------------------------------------------------------------*/
meshs_t *meshs_init(ugrid_t *u)
{
  int cc, n;

  meshs_t *m = (meshs_t *)malloc(sizeof(meshs_t));
  memset(m, 0, sizeof(meshs_t));

  m->ns2 = u->ns2;
  m->si = u->si;
  m->type = u->type;

  /*-----------------------------------------------------------------*/
  /* Vertices per cell                                               */
  m->npe = i_alloc_1d(m->ns2+1);
  m->mnpe = 0;
  for (cc = 1; cc <= u->ns2; cc++) {
    m->npe[cc] = u->npe2[cc];
    if (m->npe[cc] > m->mnpe)
      m->mnpe = m->npe[cc];
  }

  /*-----------------------------------------------------------------*/
  /* Mesh coordinates                                                */
  m->xloc = d_alloc_1d((m->ns2 * (m->mnpe+1)) + 1);
  m->yloc = d_alloc_1d((m->ns2 * (m->mnpe+1)) + 1);
  m->eloc = i_alloc_3d(m->mnpe+1, m->ns2+1, 2);

  return m;
}

/* END meshs_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates the mesh element connectivity including a list of all     */
/* geographic locations with corresponding index, and the indices of */
/* all vertices and centres of the mesh elements.                    */
/* Removes duplicate locations in a list. O(nlog(n)) operations.     */
/*-------------------------------------------------------------------*/
void create_meshs(meshs_t *m,   /* mesh structure                     */
		  double **x,
		  double **y
		  )
{
  int cc, j, jj, nn, n, next, nedge;
  double **edge;
  int **cc2n;
  int checkf = 0;

  /* Map the locations into a continuous vector                      */
  /* Get the vector size                                             */
  /* Allocate                                                        */
  nedge = m->ns2;
  for (cc = 1; cc <= m->ns2; cc++) {
    nedge += m->npe[cc];
  }
  edge = d_alloc_2d(4, nedge);

  if (m->xloc) d_free_1d(m->xloc);
  m->xloc = d_alloc_1d(nedge + m->ns2 + 1);
  if (m->yloc) d_free_1d(m->yloc);
  m->yloc = d_alloc_1d(nedge + m->ns2 + 1);
  cc2n = i_alloc_2d(m->mnpe+1, m->ns2+1);

  /* Make the vector                                                 */
  next = 0;
  for (cc = 1; cc <= m->ns2; cc++) {
    for (j = 1; j <= m->npe[cc]; j++) {
      edge[next][0] = x[cc][j];
      edge[next][1] = y[cc][j];
      edge[next][2] = (double)cc;
      edge[next][3] = (double)j;
      if (cc == checkf) printf("old%d = [%f %f]\n",j, x[cc][j], y[cc][j]);
      next++;
    }
  }
  for (cc = 1; cc <= m->ns2; cc++) {
      edge[next][0] = x[cc][0];
      edge[next][1] = y[cc][0];
      edge[next][2] = (double)cc;
      edge[next][3] = 0.0;
      if (cc == checkf) printf("old%d = [%f %f]\n",j, x[cc][0], y[cc][0]);
      next++;
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(double)*4, edgesort_double);

  if (checkf) {
    for (n = 0; n < nedge; n++) {
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      if (cc == checkf) printf("new%d = %d[%f %f]\n", j, n, edge[n][0], edge[n][1]);
    }
  }

  /* Remove the duplicates                                           */
  m->ns = 1;
  for (n = 0; n < nedge; n++) {
    nn = (n == nedge - 1) ? n - 1 : n + 1;
    if ((edge[n][0] == edge[nn][0]) &&
	(edge[n][1] == edge[nn][1])) {
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      cc2n[cc][j] = m->ns;
      /* If the last entries are duplicates, make sure these are     */
      /* also included in the mesh array.                            */
      if (n == nedge - 1) cc2n[cc][j]--;
      if (nn == nedge - 1) {
	m->xloc[m->ns] = edge[n][0];
	m->yloc[m->ns] = edge[n][1];
	m->ns++;
      }
    } else {
      m->xloc[m->ns] = edge[n][0];
      m->yloc[m->ns] = edge[n][1];
      cc = (int)edge[n][2];
      j = (int)edge[n][3];
      cc2n[cc][j] = m->ns;
      m->ns++;
    }
  }
  m->ns--;

  /* Make the mapping from indices to coordinates                    */
  for (cc = 1; cc <= m->ns2; cc++) {
     for (j = 1; j <= m->npe[cc]; j++) {
       jj = (j == m->npe[cc]) ? 1 : j + 1;
       m->eloc[0][cc][j] = cc2n[cc][j];
       m->eloc[1][cc][j] = cc2n[cc][jj];
     }
     m->eloc[0][cc][0] = m->eloc[1][cc][0] = cc2n[cc][0];
  }

  if (checkf) {
    for (cc = 1; cc <= m->ns2; cc++) {
      for (j = 0; j <= m->npe[cc]; j++) {
	double xc = x[cc][j];
	double yc = y[cc][j];
	int found = 0;
	for (n = 1; n <= m->ns; n++) {
	  if (xc == m->xloc[n] && yc == m->yloc[n]) {
	    found = 1;
	    break;
	  }
	}
	if (!found) printf("Can't find coordinate cc=%d, j=%d [%f %f]\n", cc, j, xc, yc);
	if (m->eloc[0][cc][j] < 1 || m->eloc[0][cc][j] > m->ns) 
	  printf("Invalid0 index cc=%d, j=%d [%f %f] : %d\n", cc, j, xc, yc, m->eloc[0][cc][j]);
	if (m->eloc[1][cc][j] < 1 || m->eloc[1][cc][j] > m->ns) 
	  printf("Invalid1 index cc=%d, j=%d [%f %f] : %d\n", cc, j, xc, yc, m->eloc[1][cc][j]);
	if (cc == checkf) {
	  printf("check: cc=%d j=%d %f %f\n",cc, j, m->xloc[m->eloc[0][cc][j]], m->yloc[m->eloc[0][cc][j]]);
	}

      }
    }
  }

  d_free_2d(edge);
  i_free_2d(cc2n);
}

/* END create_mesh()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees the meshs structure                                         */
/*-------------------------------------------------------------------*/
void meshs_free(meshs_t *m)
{
  i_free_1d(m->npe);
  d_free_1d(m->xloc);
  d_free_1d(m->yloc);
  i_free_3d(m->eloc);
  i_free_2d(m->neic);
  i_free_2d(m->neij);
  free((meshs_t *)m);
}

/* END meshs_free()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Edge sort function for dual row / index sorting                   */
/* Written by Darren Engwirda, 9/2017                                */
/* return -1 if "edge-1" < "edge-2"                                  */
/* return +1 if "edge-1" > "edge-2"                                  */
/* return +0 if "edge-1" = "edge-2"                                  */
/*-------------------------------------------------------------------*/
int edgesort_compare (void const *x1, void const *x2)
{
  int  const *e1 = (int *)x1 ;
  int  const *e2 = (int *)x2 ;
    
  if (e1[0]  < e2[0]) {         /* less on 1st col. */
    return -1 ;
  } else {
    if (e1[0]  > e2[0]) {       /* more on 1st col. */
      return +1 ;
    } else {
      if (e1[0] == e2[0]) {     /* test on 2nd col. */
	if (e1[1]  < e2[1]) {   /* less on 2nd col. */
	  return -1 ;
	} else {
	  if (e1[1]  > e2[1])   /* more on 2nd col. */
	    return +1 ;
	}
      }
    }
  }
  return +0 ;                   /* have exact match */
}

int edgesort_double (void const *x1, void const *x2)
{
  double  const *e1 = (double *)x1 ;
  double  const *e2 = (double *)x2 ;
    
  if (e1[0]  < e2[0]) {         /* less on 1st col. */
    return -1 ;
  } else {
    if (e1[0]  > e2[0]) {       /* more on 1st col. */
      return +1 ;
    } else {
      if (e1[0] == e2[0]) {     /* test on 2nd col. */
	if (e1[1]  < e2[1]) {   /* less on 2nd col. */
	  return -1 ;
	} else {
	  if (e1[1]  > e2[1])   /* more on 2nd col. */
	    return +1 ;
	}
      }
    }
  }
  return +0 ;                   /* have exact match */
}

/* END edge_sort_compare()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the neighbors in the mesh structure.              */
/* Uses O(nlog(n)) operations.                                       */
/* Originally written by Darren Engwirda.                            */
/*-------------------------------------------------------------------*/
void neighbour_finder_meshs(meshs_t *m)
{
  int cc, c1, c2, j, j1, j2, nedge, next;
  int **edge;

  /* Allocate                                                        */
  nedge = 0;
  for (cc = 1; cc <= m->ns2; cc++) {
    nedge += m->npe[cc];
  }
  edge = i_alloc_2d(4, nedge);

  /* Populate the edge structure                                     */
  next = 0;
  for (cc = 1; cc <= m->ns2; cc++) {
    for (j = 1; j <= m->npe[cc]; j++) {
      j1 = m->eloc[0][cc][j];
      j2 = m->eloc[1][cc][j];
      if (j1 < j2) {
	edge[next][0] = j1;
	edge[next][1] = j2;
      } else {
	edge[next][0] = j2;
	edge[next][1] = j1;
      }
      edge[next][2] = cc;
      edge[next][3] = j;
      next++;
    }
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(int)*4, edgesort_compare);

  /* Create the mappings                                             */
  m->neic = i_alloc_2d(m->ns2+1, m->mnpe+1);
  m->neij = i_alloc_2d(m->ns2+1, m->mnpe+1);
  for (cc = 1; cc <= m->ns2; cc++) {
    for (j = 1; j <= m->npe[cc]; j++) {
      m->neic[j][cc] = 0;
      m->neij[j][cc] = j;
    }
  }

  /* Create the mappings                                             */
  for (cc = 0; cc < nedge-1; cc++) {
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      c1 = edge[cc][2];
      c2 = edge[cc+1][2];
      j1 = edge[cc][3];
      j2 = edge[cc+1][3];
      m->neic[j1][c1] = c2;
      m->neij[j1][c1] = j2;
      m->neic[j2][c2] = c1;
      m->neij[j2][c2] = j1;
    }
  }
  i_free_2d(edge);
}    

/* END neighbour_finder_meshs()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds common points of a mesh centre locations and Delaunay       */
/* vertices.                                                         */
/*-------------------------------------------------------------------*/
void centre_finder(meshs_t *m, delaunay *d, int mode)
{
  int cc, c, i, n, next, nedge;
  double **edge;
  triangle *t;

  /* Map the locations into a continuous vector                      */
  /* Get the vector size                                             */
  /* Allocate                                                        */
  nedge = m->ns2 + (d->ntriangles * 3);
  edge = d_alloc_2d(4, nedge);

  /* Make the vector                                                 */
  next = 0;
  for (cc = 1; cc <= m->ns2; cc++) {
    edge[next][0] = m->xloc[m->eloc[0][cc][0]];
    edge[next][1] = m->yloc[m->eloc[0][cc][0]];
    edge[next][2] = (double)cc;
    edge[next][3] = -1.0;
    next++;
  }
  for (cc = 0; cc < d->ntriangles; cc++) {
    t = &d->triangles[cc];
    for (i = 0; i < 3; i++) {
      edge[next][0] = d->points[t->vids[i]].x;
      edge[next][1] = d->points[t->vids[i]].y;
      edge[next][2] = -1.0;
      edge[next][3] = (double)cc;
      next++;
    }
  }

  /* Sort edges, duplicates are consecutive!                         */
  qsort(edge[0], nedge, sizeof(double)*4, edgesort_double);

  /* Find common centres and make the map                            */
  /* Search forward for common Delaunay vertices                     */
  for (cc = 0; cc < nedge-1; cc++) {
    if ((edge[cc][0] == edge[cc+1][0]) &&
	(edge[cc][1] == edge[cc+1][1])) {
      if (edge[cc][2] != -1.0 && edge[cc][3] == -1.0 &&
	  edge[cc+1][2] == -1.0 && edge[cc+1][3] != -1.0) {
	c = (int)edge[cc][2];
	i = cc + 1;
	while (i < nedge && edge[i][3] != -1.0 && 
	       edge[cc][0] == edge[i][0] && edge[cc][1] == edge[i][1]) {
	  n = (int)edge[i][3];
	  m->tri2c[n] = c;
	  edge[i][3] = -1.0;
	  i++;
	}
      }
    }
  }
  /* Search backward for common Delaunay vertices                    */
  for (cc = 1; cc < nedge; cc++) {
    if ((edge[cc][0] == edge[cc-1][0]) &&
	(edge[cc][1] == edge[cc-1][1])) {
      if (edge[cc][2] != -1.0 && edge[cc][3] == -1.0 &&
	  edge[cc-1][2] == -1.0 && edge[cc-1][3] != -1.0) {
	c = (int)edge[cc][2];
	i = cc - 1;
	while (i >= 0 && edge[i][3] != -1.0 && 
	       edge[cc][0] == edge[i][0] && edge[cc][1] == edge[i][1]) {
	  n = (int)edge[i][3];
	  m->tri2c[n] = c;
	  i--;
	}
      }
    }
  }
  d_free_2d(edge);

  /* Triangulations for quads use cell centre, vertex and/or edge    */
  /* points, and these triangles dont always include the cell        */
  /* centre (e.g. a triangle may use a vertex and 2 edge points).    */
  /* These triangles will not be mapped to a centre via quicksort.   */
  /* If this is the case, then look for mapped trianges with 2 edge  */
  /* points (vid >= m->ns), then the other triangle within the cell  */
  /* sharing those edges and the cell vertex will not be mapped.     */
  /* First associate the mapped triangle edge points with the cell   */
  /* centre (there will only be two possible centres for each edge   */
  /* point), then look for unmapped trianges and map those trianges  */
  /* to the centre which has a common association from two different */
  /* edge points.                                                    */
  if (mode & I_QUAD) {
    int *mask, **vertex, ca[4], pid;
    vertex = i_alloc_2d(2, m->nd);
    memset(vertex[0], 0, m->nd * sizeof(int));
    memset(vertex[1], 0, m->nd * sizeof(int));
    mask = i_alloc_1d(m->ns2);
    memset(mask, 0, m->ns2 * sizeof(int));

    /* Associate edges in the triangulation with the mesh centre     */
    /* index for valid mapped triangles.                             */
    for (n = 0; n < d->ntriangles; n++) {
      if ((c = m->tri2c[n])) {    /* Triangle is mapped              */
	t = &d->triangles[n];
	for (i = 0; i < 3; i++) { /* Loop through triangle corners   */
	  pid = t->vids[i];       /* Index in the triangulation      */
	  if (pid >= m->ns) {     /* Edge point                      */
	    if (!vertex[pid][0]) 
	      vertex[pid][0] = c; /* Save the first mesh centre      */
	    else if(vertex[pid][0] && c != vertex[pid][0])
	      vertex[pid][1] = c; /* Save the second mesh centre     */
	  }
	}
      }
    }

    /* Loop through all trianges and set unmapped triangles          */
    for (n = 0; n < d->ntriangles; n++) {
      t = &d->triangles[n];
      if (!m->tri2c[n]) {     /* Triangle is unmapped               */
	cc = 0;
	/* Count the number of times this edge is associated with a  */
	/* centre. Keep track of the centres encountered (there will */
	/* be a maximum of four).                                    */
	memset(ca, 0, 4 * sizeof(int));
	for (i = 0; i < 3; i++) {
	  pid = t->vids[i];     
	  if ((c = vertex[pid][0])) { /* First associated centre     */
	    mask[c]++;
	    ca[cc++] = c;
	  }
	  if ((c = vertex[pid][1])) { /* Second associated centre    */
	    mask[c]++;
	    ca[cc++] = c;
	  }
	}
	/* Look for centres that are assocaited with 2 edges. The    */
	/* triangle then maps to this centre.                        */
	for (i = 0; i < 4; i++) {
	  if (mask[ca[i]] == 2) m->tri2c[n] = ca[i];
	  mask[ca[i]] = 0;
	}
      }
    }
    i_free_1d(mask);
    i_free_2d(vertex);
  }
}

/* END centre_finder()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set up a delaunay triangulation for xytoi interpolation using     */
/* delaunay_xytoi(). The triangulation is based on grid centres and  */
/* grid vertices.                                                    */
/*-------------------------------------------------------------------*/
void create_d_xytoc(ugrid_t *u)
{
  meshs_t *m = u->meshs;
  delaunay *d = NULL;
  int np = m->ns;
  int n, nn, c, i, nj;
  point *pin;
  point *pout;
  triangle *t;
  int **deic;
  int ntriangles;
  int verbose = 0;

  /* Include edge locations for quad grids                           */
  if (m->type & I_QUAD)
    np += u->nMesh2_edge;

  /* Make the points of every centre and vertex in the mesh          */
  nn = 0;
  pin = malloc(np * sizeof(point));
  for (n = 1; n <= m->ns; n++) {
    pin[nn].x = m->xloc[n];
    pin[nn++].y = m->yloc[n];
  }
  if (m->type & I_QUAD) {
    for (n = u->si; n < u->nMesh2_edge; n++) {
      pin[nn].x = u->u1x[n];
      pin[nn++].y = u->u1y[n];
    }
  }

  /* Build the triangulation                                         */
  m->nd = np;
  m->d = delaunay_build(np, pin, 0, NULL, 0, NULL);
  d = m->d;

  /* Make the mapping from triangle index to mesh index              */
  m->tri2c = i_alloc_1d(d->ntriangles);
  memset(m->tri2c, 0, d->ntriangles * sizeof(int));
  centre_finder(m, d, u->type);
  /* Old code.                                                       */
  /* Triangulations for quads dont always include the cell centre;   */
  /* use a (slower) explicit mapping method where x and y points     */
  /* are compared independently (this may not always work).          
  for (c = 1; c <= m->ns2; c++) {
    double x, y;
    x = m->xloc[m->eloc[0][c][0]];
    y = m->yloc[m->eloc[0][c][0]];
    for (n = 0; n < d->ntriangles; n++) {
      t = &d->triangles[n];
      if ((x == d->points[t->vids[0]].x || x == d->points[t->vids[1]].x ||
	   x == d->points[t->vids[2]].x) && 
	  (y == d->points[t->vids[0]].y || y == d->points[t->vids[1]].y ||
	   y == d->points[t->vids[2]].y)) {
	m->tri2c[n] = c;
      }
    }
  }
  */
  if (verbose) {
    for (n = 0; n < d->ntriangles; n++) {
      t = &d->triangles[n];
      if (verbose == 1 && m->tri2c[n]) {
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("%f %f\n",d->points[t->vids[1]].x, d->points[t->vids[1]].y);
	printf("%f %f\n",d->points[t->vids[2]].x, d->points[t->vids[2]].y);
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("NaN NaN\n");
      }
      if (verbose == 2) {
	double xm = (d->points[t->vids[0]].x + d->points[t->vids[1]].x + d->points[t->vids[2]].x) / 3.0;
	double ym = (d->points[t->vids[0]].y + d->points[t->vids[1]].y + d->points[t->vids[2]].y) / 3.0;
	printf("%f %f %d\n",xm, ym, n);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Make a polygon of the domain convex hull                        */
  /*  pout = concave_hull(pin, &np);*/
  deic = i_alloc_2d(d->ntriangles, 3);
  pout = concave_hull(m->d, NULL, &np, m->tri2c, deic);
  m->np = np;
  m->ppl = poly_create();
  for (c = 0; c < np; c++) {
    poly_add_point(m->ppl, pout[c].x, pout[c].y);
  }
  /* Print if required                                               */
  if (verbose == 3) {
    poly_t *pl = m->ppl;
    double *xs = pl->x;
    double *ys = pl->y;
    n = pl->n;
    for (i = 0; i < n; i++) {
      printf("%f %f\n", xs[i], ys[i]);
    }
  }
  free((point *)pout);

  /*-----------------------------------------------------------------*/
  /* Get a list of perimeter edges (mesh indexed).                   */
  /* Count the number of perimeter edges in each cell                */
  m->npedge = i_alloc_1d(m->ns2+1);
  memset(m->npedge, 0, (m->ns2+1) * sizeof(int));
  for (n = 0; n < d->ntriangles; ++n) {
    triangle* t = &d->triangles[n];
    c = m->tri2c[n];
    if (!c) continue;
    /* Only count triangles once, but count all perimeter edges      */
    for (nn = 0; nn < 3; nn++) {
      int nj = deic[nn][n];
      if (nj == -1 || (nj >= 0 && !m->tri2c[nj])) {
	m->npedge[c] += 1;
      }
    }
  }
  /* Get the maximum and allocate                                    */
  nn = 0;
  for (c = 1; c <= m->ns2; c++) {
    if (m->npedge[c] > nn) {
      nn = m->npedge[c];
    }
  }
  m->pedge = (point **)alloc_2d(2*nn, m->ns2+1, sizeof(point));
  /* Save the locations of the perimeter edges                       */
  memset(m->npedge, 0, (m->ns2+1) * sizeof(int));
  for (n = 0; n < d->ntriangles; ++n) {
    double x, y;
    triangle* t = &d->triangles[n];
    c = m->tri2c[n];
    if (!c) continue;
    /* Only count triangles once, but count all perimeter edges      */
    for (nn = 0; nn < 3; nn++) {
      nj = deic[nn][n];
      if (nj == -1 || (nj >= 0 && !m->tri2c[nj])) {
	point *p = m->pedge[c];
	i = (nn == 2) ? 0 : nn + 1;
	/* Loop over triange vertices, and save the edge points      */
	p[m->npedge[c]].x = d->points[t->vids[nn]].x;
	p[m->npedge[c]].y = d->points[t->vids[nn]].y;
	m->npedge[c] += 1;
	p[m->npedge[c]].x = d->points[t->vids[i]].x;
	p[m->npedge[c]].y = d->points[t->vids[i]].y;
	m->npedge[c] += 1;
      }
    }
  }
  /* Print if required                                               */
  if (verbose == 4) {
    for (c = 1; c <= m->ns2; c++) {
      if (m->npedge[c]) {
	for (n = 0; n < m->npedge[c]; n+=2) {
	  printf("%f %f\n",m->pedge[c][n].x, m->pedge[c][n].y);
	  printf("%f %f\n",m->pedge[c][n+1].x, m->pedge[c][n+1].y);
	  printf("NaN NaN\n");
	}
      }
    }
  }

  if(deic)i_free_2d(deic);
}

/* END create_d_xytoc()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds concave hull of a set of n points by finding the perimeter  */
/* of a Delaunay triangulation of those points.                      */
/* If a triangulation is supplied (v=NULL), then use this to find    */
/* the perimeter.                                                    */
/* If a set of points is supplied (din=NULL), then build the         */
/* triangulation with these.                                         */
/* msk is an optional array of triangles where msk > 1 indicates a   */
/* valid triangle is to be included in the perimeter calculation,    */
/* and all trianges with msk = 0 are excluded.                       */
/* An optional connectivity mapping may be returned (**map).         */
/*-------------------------------------------------------------------*/
point *concave_hull(delaunay *din, point *v, int *count, int *msk, int **map)
{
  int n = *count;
  delaunay *d;
  int **nei, *pe;
  int i, j, ii, nj, ni, npe = 0;
  int *index, *mask;
  point *stack;
  int dostack = 0;
  int verbose = 0;

  /* Set the triangulation                                           */
  if (din == NULL)
    d = delaunay_build(n, v, 0, NULL, 0, NULL);
  else
    d = din;
  delaunay_neighbour_finder(d, &nei);

  /* Count the perimeter edges. These will always be less than the   */
  /* number of trianges.                                             */
  npe = ii = 0;
  pe = i_alloc_1d(d->ntriangles);
  index = i_alloc_1d(d->ntriangles);
  memset(index, 0, d->ntriangles * sizeof(int));
  for (i = 0; i < d->ntriangles; ++i) {
    triangle* t = &d->triangles[i];
    if (msk && !msk[i]) continue;
    /* Only count triangles once, but count all perimeter edges      */
    for (n = 0; n < 3; n++) {
      nj = nei[n][i];
      if (nj == -1 || (msk != NULL && nj >= 0 && !msk[nj])) {
	if (!index[i]) {
	  pe[npe++] = i;
	  index[i] = 1;
	}
	ii++;
      }
    }
  }
  i_free_1d(index);

  /* Get the indices of the perimeter edge points. A mapping == -1   */
  /* indicates that triange lies on the perimeter of the             */
  /* triangulation (an interior triangle should have 3 valid maps).  */
  /* Alternitively, a mapping to a msk value of 0 also indicates a   */
  /* perimiter triangle.                                             */
  index = i_alloc_1d(2 * ii);
  ni = 0;
  for (i = 0; i < npe; ++i) {
    triangle* t;
    j = pe[i];
    t = &d->triangles[j];
    if (msk && !msk[j]) continue;
    for (n = 0; n < 3; n++) {
      nj = nei[n][j];
      if (nj == -1 || (msk != NULL && nj >= 0 && !msk[nj])) {
	int nn = (n == 2) ? 0 : n + 1;
	/* Note, this assumes triangles are ordered anti-clockwise,  */
	/* with the first index of the perimeter edge being in the   */
	/* mapping direction.                                        */   
	index[ni] = t->vids[n];
	index[ni+1] = t->vids[nn];
	if (verbose == 1) {
	  /* Points
	  printf("%f %f p%d\n",d->points[t->vids[n]].x, d->points[t->vids[n]].y, ni);
	  printf("%f %f p%d\n",d->points[t->vids[nn]].x, d->points[t->vids[nn]].y, ni);
	  */
	  /* Polygons                                                */
	  printf("%f %f\n",d->points[t->vids[n]].x, d->points[t->vids[n]].y);
	  printf("%f %f\n",d->points[t->vids[nn]].x, d->points[t->vids[nn]].y);
	  printf("NaN NaN\n");
	}
	ni+=2;
      }
    }
  }

  /* Make a continuous polygon                                       */
  i_free_1d(pe);
  pe = i_alloc_1d(d->npoints);
  memset(pe, 0, d->npoints * sizeof(int));
  mask = i_alloc_1d(2 * ni);
  memset(mask, 0, 2 * ni * sizeof(int));
  n = 0;
  pe[index[n]] = 1;
  mask[n++] = index[0];
  mask[n++] = index[1];
  for (j = 0; j < ni; j+=2) {
    for (i = 0; i < ni; i+=2) {
      ii = index[i];
      if (!pe[ii]) {
	if (mask[n-1] == index[i]) {
	  mask[n++] = index[i+1];
	  pe[index[i]] = 1;
	  break;
	} else if (mask[n-1] == index[i+1]) {
	  mask[n++] = index[i];
	  pe[index[i+1]] = 1;
	  break;
	}
      }
    }
  }

  /* Save the point locations                                        */
  stack = (point *)malloc((n + 1)* sizeof(point));
  for (i = 0; i < n; i++) {
    stack[i].x = d->points[mask[i]].x;
    stack[i].y = d->points[mask[i]].y;
    if (verbose == 2) printf("%f %f\n", stack[i].x, stack[i].y);
  }
  stack[n].x = stack[0].x;    /* Close the polygon                   */
  stack[n].y = stack[0].y;
  if (verbose == 2) printf("%f %f\n", stack[n].x, stack[n].y);
  *count = n;

  /* Save the neighbour maps                                         */
  if (map) {
    for (i = 0; i < d->ntriangles; ++i) {
      for (n = 0; n < 3; n++) map[n][i] = nei[n][i];
    }
  }

  i_free_1d(pe);
  i_free_1d(index);
  i_free_1d(mask);
  i_free_2d(nei);
  if(v != NULL) delaunay_destroy(d);
  return stack;
}

/* END concave_hull()                                                */
/*-------------------------------------------------------------------*/



