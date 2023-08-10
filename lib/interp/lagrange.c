/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/interp/lagrange.c
 *
 *  \brief Semi-Lagrange interpolation routines
 *
 *  Routines which deal setting up semi-Lagrange streamline tracing
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lagrange.c 6876 2021-07-29 00:40:40Z her127 $
 */

/* Interpolation is set up by providing a UGRID netCDF file.         */
/* Source location is the end point of a streamline.                 */
/* Destination location is the start point of a streamline.          */
/* Streamline source points may be computed (e.g. particle tracking) */
/* or tracer values may be interpolated at the source.               */
/* Note: the points used in velocity interpolation correspond to     */
/* every mesh point, so that the streamline can be well resolved     */
/* when traced. The locations used in interpolation are set up in    */
/* this routine, and streamline source points are computed in        */
/* semi_lagrange_source() (at every mesh location) or                */
/* semi_lagrange_atxy() (at a given (x,y,z) location).               */
/* Tracers may be interpolated using values on a subset of the full  */
/* mesh. The locations used in interpolation are set up in this      */
/* routine. The actual locations interpolated may be a different set */
/* of points again, and tracers are interpolated at these points in  */
/* semi_lagrange_tr() at a given mesh index and layer.               */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "ems.h"

/* Lagrangian method */
#define L_LINEAR      0x001
#define L_SIB         0x002
#define L_NONSIB      0x004
#define L_CUBIC       0x008
#define L_LSQUAD      0x010
#define L_BILIN       0x020
#define L_BAYLIN      0x040
#define L_FG          0x080
#define L_LSLIN       0x100
#define L_NRST        0x200

#define SMALL 1e-10             /* Minimum value for velocity bounds */
#define RADIUS 6370997.0
#define ECC 0.0
#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

ugrid_t *ugrid_alloc(void);
lagrange_t *lagrange_alloc(void);
void ugrid_maps(ugrid_t *u, int fid);
void set_counts(int dimf, size_t *count, int d1, int d2);
void delaunay_sed(lagrange_t *l, int vid, double **v);
void grid_spec_init_delaunay(lagrange_t *l, int vid);
void lagrange_free(lagrange_t *l);
static npts_t* npt_create(int n);
static void reset_npt(lagrange_t *l, delaunay **d, npts_t **npt) ;
static void reset_interp_points(meshs_t *mesh, lagrange_t *l, npts_t **npt);
static void rebuild_pt(meshs_t *mesh, lagrange_t *l, delaunay **d, npts_t **npt);
static double lag_interp(lagrange_t *l, GRID_SPECS **gs, double x, double y, double z, int c, int k, int co, int vid);
int lag_pos(meshs_t *m, lagrange_t *l, int *ci, int *ki, double u,
	    double v, double w, double *cx, double *cy, double *cz,
	    double dt, int is_geog);
void lag_index(meshs_t *m, lagrange_t *l, double x, double y, double z, 
	       int *c, int *k, int first);
int find_intersect(point s, point d, point e1, point e2, double *xi, double *yi);
int nearest_p(meshs_t *m, double xin, double yin, int *vec, int nvec);
int ugrid_read_3d(ugrid_t *u, int id, char *name, double **p, int dump);
int ugrid_read_2d(ugrid_t *u, int id, char *name, double *p, int dump);

/*-------------------------------------------------------------------*/
/* Reads in grid information for setting up semi-Lagrange streamline */
/* tracking from an unstructured UGRID input file.                   */
/*-------------------------------------------------------------------*/
ugrid_t *ugrid_open(char *name, int options)
{
  int fid;
  int ncerr;
  ugrid_t *u;
  char buf[MAXSTRLEN];

  u = ugrid_alloc();

  /* Clear the string arrays in case data in file isn't zero terminated */
  sprintf(u->vers, "%c", '\0');
  sprintf(u->chead, "%c", '\0');
  sprintf(u->phead, "%c", '\0');
  sprintf(u->rcode, "%c", '\0');
  sprintf(u->timeunit, "%c", '\0');
  sprintf(u->lenunit, "%c", '\0');
  sprintf(u->projection, "%c", '\0');
  u->type = I_HEX;
  u->options = options;

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(name, NC_NOWRITE, &fid)) != NC_NOERR) {
    warn("Can't find input file %s\n", name);
    quit((char *)nc_strerror(ncerr));
  }

  /* Get dimensions                                                  */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "Mesh2_layerfaces"), &u->kgridsize);
  if ((int)u->kgridsize == 0) quit("kgridsize = 0: is %s a UGRID file?", name);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "Mesh2_layers"), &u->kcentresize);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_node"), &u->nMesh2_node);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_edge"), &u->nMesh2_edge);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMesh2_face"), &u->nMesh2_face);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nMaxMesh2_face_nodes"), &u->nMaxMesh2_face_nodes);
  nc_get_att_int(fid, NC_GLOBAL, "start_index", &u->oset);
  u->oset = (u->oset) ? 0 : 1;
  u->si = (u->oset) ? 0 : 1;

  /* Get global attributes */
  if (options & I_GLOB) {
    nc_get_att_text(fid, NC_GLOBAL, "title", u->chead);
    nc_get_att_text(fid, NC_GLOBAL, "paramhead", u->phead);
    nc_get_att_text(fid, NC_GLOBAL, "paramfile", u->prmname);
    nc_get_att_text(fid, NC_GLOBAL, "version", u->vers);
    nc_get_att_text(fid, NC_GLOBAL, "Run_code", u->rcode);
  }
  if (nc_get_att_int(fid, NC_GLOBAL, "NCE1", &ncerr) >= 0)
    u->type = I_QUAD;
  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", u->timeunit);
  nc_get_att_text(fid, ncw_var_id(fid, "Mesh2_layerfaces"), "units",
		  u->lenunit);

  /* Check if geographic. */
  nc_get_att_text(fid, ncw_var_id(fid, "Mesh2_node_x"), "projection", u->projection);
  if (u->projection != NULL) {
    ts_set_default_proj_type(u->projection);
  }

  /* time independent variables */
  /* Vertical layer structure */
  if (u->cellz == NULL) {
    size_t start[4];
    size_t count[4];
    u->cellz = d_alloc_1d((int)u->kcentresize);
    count[0] = u->kcentresize;
    start[0] = 0;
    count[1] = start[1] = 0;
    count[2] = start[2] = 0;
    count[3] = start[3] = 0;
    nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_layers"), start, count,
		       u->cellz);
    u->gridz = d_alloc_1d(u->kgridsize);
    count[0] = u->kgridsize;
    nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_layerfaces"), start, count,
		       u->gridz);
  }

  /* Set the maps                                                    */
  ugrid_maps(u, fid);
  nc_close(fid);
  return (u);
}

/* END ugrid_open()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates ugrid memory                                            */
/*-------------------------------------------------------------------*/
ugrid_t *ugrid_alloc(void)
{
  ugrid_t *ug = (ugrid_t *)malloc(sizeof(ugrid_t));
  memset(ug, 0, sizeof(ugrid_t));
  return ug;
}

/* END ugrid_alloc()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the ugrid geometry from the information in the input      */
/* netCDF file.                                                      */
/*-------------------------------------------------------------------*/
void ugrid_maps(ugrid_t *u, int fid)
{
  int c, cc, v, n, nn;       /* Counters                             */
  size_t start[4];           /* netCDF start vector for reads        */
  size_t count[4];           /* netCDF count vector for reads        */
  int **ic2v;                /* Centre to vertex mapping             */
  int face_dim = 0;          /* Start index                          */
  int fill;                  /* Fill value for c2v                   */
  int i1;
  double *bathy;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;

  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*-----------------------------------------------------------------*/
  /* Get dimensions and allocate                                     */
  nc_get_att_int(fid, NC_GLOBAL, "face_dim", &face_dim);

  u->c2v = i_alloc_2d(u->nMesh2_face, u->nMaxMesh2_face_nodes);
  if (face_dim)
    ic2v = i_alloc_2d(u->nMesh2_face, u->nMaxMesh2_face_nodes);
  else
    ic2v = i_alloc_2d(u->nMaxMesh2_face_nodes, u->nMesh2_face);

  u->cellx = d_alloc_1d(u->nMesh2_face);
  u->celly = d_alloc_1d(u->nMesh2_face);
  u->gridx = d_alloc_1d(u->nMesh2_node);
  u->gridy = d_alloc_1d(u->nMesh2_node);
  u->u1x = d_alloc_1d(u->nMesh2_edge);
  u->u1y = d_alloc_1d(u->nMesh2_edge);
  u->index = i_alloc_1d(u->nMesh2_face);
  bathy = d_alloc_1d(u->nMesh2_face);

  u->ns2 = (int)u->nMesh2_face - u->si;
  u->npem = (int)u->nMaxMesh2_face_nodes - u->si;
  u->npe2 = i_alloc_1d(u->ns2 + 1);
  u->x = d_alloc_2d(u->npem + 1, u->ns2 + 1);
  u->y = d_alloc_2d(u->npem + 1, u->ns2 + 1);
  u->bathy = d_alloc_1d(u->ns2 + 1);

  /*-----------------------------------------------------------------*/
  /* Get the locations of the cell centre                            */
  count[0] = u->nMesh2_face;
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_face_x"), start, count,
                     u->cellx);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_face_y"), start, count,
                     u->celly);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_depth"), start, count,
                     bathy);
  nc_get_att_double(fid, ncw_var_id(fid, "Mesh2_depth"), "_FillValue", &u->bfill);

  /*
  if ((i1 = ncw_var_id(fid, "Mesh2_index")) >= 0)
    nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_index"), start, count, u->index);  
  else {
  */
    for (cc = 0; cc < u->nMesh2_face; cc++)
      u->index[cc] = cc + 1;

  /*-----------------------------------------------------------------*/
  /* Get the locations of the cell vertices                          */
  count[0] = u->nMesh2_node;
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_node_x"), start, count,
                     u->gridx);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_node_y"), start, count,
                     u->gridy);
  /* Get the locations of the cell edges                             */
  count[0] = u->nMesh2_edge;
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_edge_x"), start, count,
                     u->u1x);
  nc_get_vara_double(fid, ncw_var_id(fid, "Mesh2_edge_y"), start, count,
                     u->u1y);

  /*-----------------------------------------------------------------*/
  /* Get the mapping from centre to vertices and populate the        */
  /* parameters structure.                                           */
  /* Note that the input of the mesh topology to create the input    */
  /* file may list the mesh elements in arbitary order, and this     */
  /* will determine which faces are neighbours in the cell arrays,   */
  /* and which edges are associated with which cells. This in turn   */
  /* will determine on which cells or faces data read in from the    */
  /* input netCDF file is placed. Therefore, make sure to put the    */
  /* grid locations and bathymetry in the correct order via index[]  */
  /* so that the same mappings are created as was used to generate   */
  /* the input netCDF file.                                          */
  set_counts(face_dim, count, u->nMaxMesh2_face_nodes, u->nMesh2_face);
  nc_get_vara_int(fid, ncw_var_id(fid, "Mesh2_face_nodes"), start, count,
		   ic2v[0]);
  nc_get_att_int(fid, ncw_var_id(fid, "Mesh2_face_nodes"), "_FillValue", &fill);

  for (cc = u->si; cc < u->nMesh2_face; cc++)
    for (n = u->si; n < u->nMaxMesh2_face_nodes; n++) {
      if (face_dim)
	u->c2v[n][cc] = ic2v[n][cc] + u->oset;
      else
	u->c2v[n][cc] = ic2v[cc][n] + u->oset;
    }

  /* Copy the mesh information into x and y, with indices from 1:ns2 */
  for (cc = u->si; cc < u->nMesh2_face; cc++) {
    c = u->index[cc];
    u->npe2[c] = 0;
    u->x[c][0] = u->cellx[cc];
    u->y[c][0] = u->celly[cc];
    u->bathy[c] = bathy[cc];
    /*if(isnan(u->bathy[c])) u->bathy[c] = NaN;*/
    for (n = u->si; n < u->nMaxMesh2_face_nodes; n++) {
      nn = n + u->oset;
      if (u->c2v[n][cc] > 0 && u->c2v[n][cc] != fill + u->oset) {
	u->npe2[c]++;
	v = u->c2v[n][cc] - u->oset;
	u->x[c][nn] = u->gridx[v];
	u->y[c][nn] = u->gridy[v];
      }
    }
  }

  /* Create the mesh structure                                       */
  u->meshs = meshs_init(u);
  create_meshs(u->meshs, u->x, u->y);
  neighbour_finder_meshs(u->meshs);
  /* Create the Delaunay structure for mapping (x,y) to c            */
  create_d_xytoc(u);

  d_free_1d(bathy);
  i_free_2d(ic2v);
  d_free_1d(u->u1x);
  d_free_1d(u->u1y);
}

/* END ugrid_maps()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees the UGRID structure                                         */
/*-------------------------------------------------------------------*/
void ugrid_free(ugrid_t *u)
{
  int k;

  if (u->meshs) meshs_free(u->meshs);
  d_free_1d(u->cellz);
  d_free_1d(u->gridz);
  d_free_1d(u->cellx);
  d_free_1d(u->celly);
  d_free_1d(u->gridx);
  d_free_1d(u->gridy);
  i_free_1d(u->npe2);
  d_free_2d(u->x);
  d_free_2d(u->y);
  d_free_1d(u->bathy);
  i_free_1d(u->index);
  free((ugrid_t *)u);
}

/* END ugrid_free()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets count vector for netCDF input                                */
/*-------------------------------------------------------------------*/
void set_counts(int dimf, size_t *count, int d1, int d2)
{
  if (dimf) {
    count[0] = d1;
    count[1] = d2;
  } else {
    count[0] = d2;
    count[1] = d1;
  }
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up a grid_specs structure for Lagrangian tracking and        */
/* interpolation. The streamline is traced on the supplied UGRID     */
/* file geometry, tracers may be interpolated on a subset of these   */
/* if nvec != 0. In this case vec is a vector of the UGRID indices   */
/* that will be interpolted.                                         */
/*-------------------------------------------------------------------*/
void lagrange_init(ugrid_t *u,      /* UGRID structure               */
		   int nivec,       /* Number of tracer points       */
		   int *ivec        /* Tracer point vector           */
		   )
{
  int cc, c, k, n, m, nd;
  int *nk;
  point **p, **pb;
  double v;
  delaunay *d;
  int dot = (ivec != NULL && nivec != 0) ? 1 : 0;
  lagrange_t *l;
  meshs_t *mesh = u->meshs;
  int ns2  = mesh->ns2;
  int nvec, *vec;
  int nbm, *nb, **bvec;
  int *mask;
  int verbose = 0;

  /*-----------------------------------------------------------------*/  
  /* Allocate                                                        */
  l = lagrange_alloc();
  l->nz = u->kcentresize + 1;             /* Water column + sediment */
  l->bot = i_alloc_1d(ns2 + 1);
  l->d2u = i_alloc_1d(ns2 + 1);
  l->options = u->options;
  /* Note; l->layers are horizontal surfaces, not layers with an     */
  /* associated thickness. Each layer resides in the middle of the   */
  /* model layer, with depth u->cellz.                               */
  l->layers = d_alloc_1d(l->nz);
  for (k = 1; k < l->nz; k++)
    l->layers[k] = u->cellz[k-1];
  l->dol = i_alloc_1d(l->nz);
  memset(l->dol, 0, l->nz * sizeof(int));
  /* Set the momentum interpolation scheme                           */
  strcpy(l->momsr, "linear");
  l->mosl = L_LINEAR;
  /* Set the interpolation schemes                                   */
  strcpy(l->trasr, "linear");
  l->osl = L_LINEAR;

  /*-----------------------------------------------------------------*/
  /* If interpolation schemes are different, build a separate        */
  /* Delaunay structure for tracers.                                 */
  if (dot) {
    nvec = nivec;
    vec = ivec;
    strcpy(l->trasr, l->momsr);
    l->osl = l->mosl;
  } else if (l->osl != l->mosl) {
    /* Rebuild the tracer point trianges. Momentum and tracers may   */
    /* have different interpolation schemes, so tracers need its own */
    /* point triangles. Only for the full mesh (not subsets).        */
    dot = 2;
    nvec = ns2;
    vec = i_alloc_1d(nvec);
    for (c = 1; c <= ns2; c++) vec[c-1] = c;
  }

  /*-----------------------------------------------------------------*/
  /* Set up the triangulation for each layer (velocity).             */
  /* Allocate known sized arrays.                                    */
  nk = i_alloc_1d(l->nz);
  memset(nk, 0, l->nz * sizeof(int));
  nb = i_alloc_1d(l->nz);
  memset(nb, 0, l->nz * sizeof(int));
  for (c = 1; c <= ns2; c++) l->bot[c] = -1;

  /* Count the maximum number of points and perimeter segment        */
  /* vertices. This occurs in the top layer.                         */
  nd = nbm = 0;
  k = l->nz-1;
  for (c = 1; c <= ns2; c++) {
    if (u->bathy[c] == u->bfill || isnan(u->bathy[c])) continue;
    if (u->bathy[c] <= u->cellz[k-1]) {
      /* Wet cells                                                   */
      nd++;
      /* Check for perimeter cells (cn = 0)                          */
      for (n = 1; n <= mesh->npe[c]; n++) {
	int cn = mesh->neic[n][c];
	/* Perimeter cell below the surface                          */
	if (u->bathy[cn] > u->cellz[k-1]) cn = 0;
	if (cn == 0) {
	  nd += 3;     /* Number of points in the triangulation      */
	  nbm += 3;    /* Number of boundary points                  */
	}
      }
    }
  }

  /* Get the points and perimeter segments for each layer            */
  l->u2d = i_alloc_2d(nd + 1, l->nz);
  l->b2d = i_alloc_2d(nd + 1, l->nz);
  p = (point **)alloc_2d(nd, l->nz, sizeof(point));
  pb = (point **)alloc_2d(nbm, l->nz, sizeof(point));
  bvec = i_alloc_2d(nbm, l->nz);
  m = ns2 * (mesh->mnpe + 1);
  mask = i_alloc_1d(m);
  for (k = 1; k < l->nz; k++) {
    for (c = 1; c <= ns2; c++) {
      l->u2d[k][cc] = -1;
      l->b2d[k][cc] = -1;
    }
    memset(mask, 0, m * sizeof(int));
    if (l->options & I_2D && k != l->nz-1) continue;
    for (c = 1; c <= ns2; c++) {
      if (u->bathy[c] == u->bfill || isnan(u->bathy[c])) continue;
      if (u->bathy[c] <= u->cellz[k-1]) {
	/* Wet cells                                                 */	
	p[k][nk[k]].x = mesh->xloc[mesh->eloc[0][c][0]];
	p[k][nk[k]].y = mesh->yloc[mesh->eloc[0][c][0]];
	l->u2d[k][nk[k]] = c;	
	if (k == l->nz-1) l->d2u[c] = nk[k];
	nk[k]++;
	if (l->bot[c] < 0) l->bot[c] = k;

	/* Ghost cells, in 3 dimensions.                             */
	/* The triangulation used for interpolation does not have    */
	/* unique ghost cells mapping across land boundaries.        */
	/* The geographic location of ghost cells is set to the cell */
	/* edge so as to account for multiple ghost cells in the     */
	/* interpolation triangulation. This means the boundary      */
	/* condition for velocity is v=0 at ghost cells rather than  */
	/* v=-vi. For tracers we use v=vi. Ghost cells are given a   */
	/* mapping u2d=-c to distinguish them, where c is the        */
	/* interior wet cell.                                        */
	/* Check for perimeter cells (cn = 0)                        */
	for (n = 1; n <= mesh->npe[c]; n++) {
	  int cn = mesh->neic[n][c];

	  /* Perimeter cell below the surface                        */
	  if (cn && u->bathy[cn] > u->cellz[k-1]) cn = 0;

	  if (cn == 0) {
	    /* Set the Delaunay point to the edge centre             */
	    p[k][nk[k]].x = 0.5 * (mesh->xloc[mesh->eloc[0][c][n]] + 
				   mesh->xloc[mesh->eloc[1][c][n]]);
	    p[k][nk[k]].y = 0.5 * (mesh->yloc[mesh->eloc[0][c][n]] + 
				   mesh->yloc[mesh->eloc[1][c][n]]);
	    l->u2d[k][nk[k]] = -c;
	    /* Perimeter locations only                              */
	    pb[k][nb[k]].x = p[k][nk[k]].x;
	    pb[k][nb[k]].y = p[k][nk[k]].y;
	    bvec[k][nb[k]] = c;
	    nb[k]++;
	    nk[k]++;

	    /* Save the mesh vertex locations                        */
	    if (!mask[mesh->eloc[0][c][n]]) {
	      p[k][nk[k]].x = mesh->xloc[mesh->eloc[0][c][n]];
	      p[k][nk[k]].y = mesh->yloc[mesh->eloc[0][c][n]];
	      /* Save the centres assocaited with this vertex.       */
	      l->u2d[k][nk[k]] = -c;
	      /* Save the perimeter points (used only if dot != 0)   */
	      pb[k][nb[k]].x = p[k][nk[k]].x;
	      pb[k][nb[k]].y = p[k][nk[k]].y;
	      /* Save the centre assocaied with this point           */
	      bvec[k][nb[k]] = c;
	      nb[k]++;
	      /* Increment                                           */
	      mask[mesh->eloc[0][c][n]] = nk[k];
	      nk[k]++;
	    } else {
	      /* Vertices are associated with multiple mesh centres; */
	      /* save the first centre found                         */
	      l->b2d[k][mask[mesh->eloc[0][c][n]]] = -c;
	    }

	    if (!mask[mesh->eloc[1][c][n]]) {
	      p[k][nk[k]].x = mesh->xloc[mesh->eloc[1][c][n]];
	      p[k][nk[k]].y = mesh->yloc[mesh->eloc[1][c][n]];
	      if (l->u2d[k][nk[k]]) l->b2d[k][nk[k]] = l->u2d[k][nk[k]];
	      l->u2d[k][nk[k]] = -c;
	      pb[k][nb[k]].x = p[k][nk[k]].x;
	      pb[k][nb[k]].y = p[k][nk[k]].y;
	      bvec[k][nb[k]] = c;
	      nb[k]++;
	      mask[mesh->eloc[0][c][n]] = nk[k];
	      nk[k]++;
	    } else {
	      l->b2d[k][mask[mesh->eloc[1][c][n]]] = -c;
	    }
	  }
	}
      }
    }
  }

  /* Sediments                                                       */
  k = l->nz-1;
  nk[0] = nk[k];
  for (cc = 0; cc < nk[k]; cc++) {
    p[0][cc].x = p[k][cc].x;
    p[0][cc].y = p[k][cc].y;
  }

  /* Fill the points array and create the triangulation for every    */
  /* layer.                                                          */
  l->d = (delaunay **)calloc(l->nz, sizeof(delaunay *));
  for (k = 0; k < l->nz; k++) {
    l->d[k] = NULL;

    if (nk[k] == 0) continue;
    l->d[k] = delaunay_build(nk[k], p[k], 0, NULL, 0, NULL);
    l->dol[k] = nk[k];
    for (c = 0; c < nk[k]; c++) {
      l->d[k]->points[c].v = d_alloc_1d(4);
    }
  }
  for (k = 0; k < l->nz; k++) {
    if (l->dol[k]) {
      delaunay *d = l->d[k];
      d->xmin = l->d[l->nz-1]->xmin;
      d->xmax = l->d[l->nz-1]->xmax;
      d->ymin = l->d[l->nz-1]->ymin;
      d->ymax = l->d[l->nz-1]->ymax;
      d->ptf = 0;
    }
  }

  /* Rebuild the momentum point trianges                             */
  rebuild_pt(mesh, l, l->d, l->nptm);
  free((point **)p);
  p = NULL;

  /*-----------------------------------------------------------------*/
  /* Set up the triangulation for each layer (tracers).              */
  if (dot) {
    l->nvec = nvec;
    memset(nk, 0, l->nz * sizeof(int));

    p = (point **)alloc_2d(nvec + nbm, l->nz, sizeof(point));
    l->t2d = i_alloc_2d(ns2 + 1, l->nz);
    
    for (k = 1; k < l->nz; k++) {
      if (l->options & I_2D && k != l->nz-1) continue;

      /* Perimeter points                                            */
      for (cc = 0; cc < nb[k]; cc++) {
	p[k][nk[k]].x = pb[k][cc].x;
	p[k][nk[k]].y = pb[k][cc].y;
	if (l->options & I_NRST)
	  l->t2d[k][nk[k]] = -nearest_p(mesh, p[k][nk[k]].x, p[k][nk[k]].y, vec, nvec);
	else {
	  l->t2d[k][nk[k]] = -bvec[k][cc];
	}
	nk[k]++;
      }

      /* Designated points                                           */
      for (cc = 0; cc < nvec; cc++) {
	c = vec[cc];
	if (u->bathy[c] == u->bfill || isnan(u->bathy[c])) continue;
	if (u->bathy[c] <= u->cellz[k-1]) {
	  /*if (u->bathy[c] >= u->cellz[k-1] && u->bathy[c] < u->cellz[k]) {*/
	  /* Wet cells only for a tracer subset. The interpolation   */
	  /* scheme must account for extrapolation to points outside */
	  /* the triangulation via the perimeter points.             */
	  p[k][nk[k]].x = mesh->xloc[mesh->eloc[0][c][0]];
	  p[k][nk[k]].y = mesh->yloc[mesh->eloc[0][c][0]];
	  l->t2d[k][nk[k]] = c;
	  nk[k]++;
	}
      }
    }

    /* Sediments                                                     */
    k = l->nz-1;
    nk[0] = nk[k];
    for (cc = 0; cc < nk[k]; cc++) {
      p[0][cc].x = p[k][cc].x;
      p[0][cc].y = p[k][cc].y;
    }

    /* Fill the points array and create the triangulation for every  */
    /* layer.                                                        */
    l->td = (delaunay **)calloc(l->nz, sizeof(delaunay *));
    for (k = 0; k < l->nz; k++) {
      l->td[k] = NULL;
      if (nk[k] == 0) continue;
      l->td[k] = delaunay_build(nk[k], p[k], 0, NULL, 0, NULL);
      for (c = 0; c < nk[k]; c++) {
	l->td[k]->points[c].v = d_alloc_1d(4);
      }
    }
    for (k = 0; k < l->nz; k++) {
      if (l->td[k] != NULL) {
	delaunay *d = l->td[k];
	d->xmin = l->td[l->nz-1]->xmin;
	d->xmax = l->td[l->nz-1]->xmax;
	d->ymin = l->td[l->nz-1]->ymin;
	d->ymax = l->td[l->nz-1]->ymax;
	d->ptf = 0;
      }
    }
    if (dot == 2) rebuild_pt(mesh, l, l->td, l->nptt);

    if (verbose == 2) {
      k = l->nz-1;
      delaunay *d = l->td[k];
      for (n = 0; n < d->ntriangles; n++) {
	triangle *t = &d->triangles[n];
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("%f %f\n",d->points[t->vids[1]].x, d->points[t->vids[1]].y);
	printf("%f %f\n",d->points[t->vids[2]].x, d->points[t->vids[2]].y);
	printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
	printf("NaN NaN\n");
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the data values, for u, v, w and tracers                    */
  for (c = 0; c <= 3; c++)
    delaunay_fill(l, c, NULL, D_ZERO);

  /*-----------------------------------------------------------------*/
  /* Set up the GRID_SPEC structure for each layer and each          */
  /* variable, where all structures use the same triangulation, d.   */
  for (c = 0; c <= 3; c++)
    grid_spec_init_delaunay(l, c);

  /*-----------------------------------------------------------------*/
  /* Set the layer thickness                                         */
  l->dzz = d_alloc_2d(ns2 + 1, l->nz);
  for (k = 1; k < l->nz; k++) {
    for (c = 1; c <= ns2; c++) {
      if (k != l->bot[c]) {
	l->dzz[k][c] = l->layers[k] - l->layers[k-1];
      }
    }
  }
  /* Bottom thickness                                                */
  for (c = 1; c <= ns2; c++) {
    k = l->bot[c];
    l->dzz[k][c] = l->layers[k] - u->bathy[c];
  }

  l->cx = d_alloc_2d(ns2 + 1, l->nz);
  l->cy = d_alloc_2d(ns2 + 1, l->nz);
  l->cz = d_alloc_2d(ns2 + 1, l->nz);
  l->cl = i_alloc_2d(ns2 + 1, l->nz);
  l->ck = i_alloc_2d(ns2 + 1, l->nz);
  for (k = 1; k < l->nz; k++) {
    if (!(l->dol[k])) continue;
    memset(l->cx[k], 0, ns2 * sizeof(double));
    memset(l->cy[k], 0, ns2 * sizeof(double));
    memset(l->cz[k], 0, ns2 * sizeof(double));
    memset(l->cl[k], 0, ns2 * sizeof(int));
    memset(l->ck[k], 0, ns2 * sizeof(int));
  }

  if (verbose == 1) {
    k = l->nz - 1;
    delaunay *d = l->d[k];
    for (n = 0; n < d->ntriangles; n++) {
      triangle *t = &d->triangles[n];
      printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      printf("%f %f\n",d->points[t->vids[1]].x, d->points[t->vids[1]].y);
      printf("%f %f\n",d->points[t->vids[2]].x, d->points[t->vids[2]].y);
      printf("%f %f\n",d->points[t->vids[0]].x, d->points[t->vids[0]].y);
      printf("NaN NaN\n");
    }
  }

  if (p) free((point **)p);
  u->l = l;
  i_free_1d(nk);
  i_free_1d(nb);
  i_free_1d(mask);
  i_free_2d(bvec);
}

/* END lagrange_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the near point in a vector of locations to a given location */
/*-------------------------------------------------------------------*/
int nearest_p(meshs_t *m, double xin, double yin, int *vec, int nvec)
{
  int cc, c, cm;
  double dm = HUGE;

  for (cc = 0; cc < nvec; cc++) {
    c = vec[cc];
    double x = xin - m->xloc[m->eloc[0][c][0]];
    double y = xin - m->xloc[m->eloc[0][c][0]];
    double dist = sqrt(x * x + y * y);
    if (dist < dm) {
      dm = dist;
      cm = cc;
    }
  }
  return(cm);
}

/* END nearest_p()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates lagrange memory                                         */
/*-------------------------------------------------------------------*/
lagrange_t *lagrange_alloc(void)
{
  lagrange_t *l = (lagrange_t *)malloc(sizeof(lagrange_t));
  memset(l, 0, sizeof(lagrange_t));
  return l;
}

/* END lagrange_alloc()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Frees the Lagrange structure                                      */
/*-------------------------------------------------------------------*/
void lagrange_free(lagrange_t *l)
{
  int k;

  /* GRID_SPEC structures                                            */
  for (k = 0; k < l->nz; k++) {
    continue;  /* FR to check */
    if (l->dol[k]) {
      grid_specs_destroy(l->gsx[k]);
      grid_specs_destroy(l->gsy[k]);
      grid_specs_destroy(l->gsz[k]);
      grid_specs_destroy(l->gst[k]);
    }
  }
  /* Delaunay structures                                             */
  /* FR to check
  for (k = 0; k < l->nz; k++) delaunay_destroy(l->d[k]);
  if (l->nvec)
    for (k = 0; k < l->nz; k++) delaunay_destroy(l->td[k]);
  */
  i_free_1d(l->bot);
  i_free_2d(l->u2d);
  i_free_1d(l->d2u);
  if (l->t2d) i_free_2d(l->t2d);
  if (l->b2d) i_free_2d(l->b2d);
  d_free_1d(l->layers);
  i_free_1d(l->dol);
  d_free_2d(l->cx);
  d_free_2d(l->cy);
  d_free_2d(l->cz);
  d_free_2d(l->dzz);
  i_free_2d(l->cl);
  i_free_2d(l->ck);
  free((lagrange_t *)l);
}

/* END lagrange_free()                                               */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Resets the weights for interpolation                              */
/*-------------------------------------------------------------------*/
void delaunay_fill(lagrange_t *l, int vid, double **v, int mode)
{
  int c, cb, cc, k, n, m;
  GRID_SPECS **gs;
  int ef = 0;
  int osl;
  int **map = l->u2d;
  int **bap = l->b2d;
  delaunay **d = l->d;

  /*-----------------------------------------------------------------*/
  /* Set maps                                                        */
  if (vid == 3 && l->nvec) {
    map = l->t2d;
    bap = NULL;
    d = l->td;
  }

  /*-----------------------------------------------------------------*/
  /* Zero the Delaunay data if required                              */
  if (mode & D_ZERO) {
    for (k = 1; k < l->nz; k++) {
      if (l->dol[k]) {
	for (cc = 0; cc < d[k]->npoints; cc++) {
	  d[k]->points[cc].v[vid] = 0.0;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Fill the points array in the delaunay structure.                */
  if (mode & D_FILL) {
    for (k = 1; k < l->nz; k++) {
      if (!(l->dol[k])) continue;
      for (cc = 0; cc < d[k]->npoints; cc++) {
	c = map[k][cc];
	/* Wet cells                                                 */
	if (c > 0) {
	  d[k]->points[cc].z = v[k][c];
	} else {
	  /* Ghost cells (points on the mesh perimeter)              */
	  c = abs(c);
	  if (vid == 3) {
	    d[k]->points[cc].z = v[k][c];
	    if (bap && (cb = abs(bap[k][cc]))) {
	      d[k]->points[cc].z = 0.5 * (d[k]->points[cc].z + v[k][cb]);
	    }
	  } else 
	    d[k]->points[cc].z = 0.0;
	}
      }
    }
    delaunay_sed(l, vid, v);
  }

  /*-----------------------------------------------------------------*/
  /* Rebuild the weights                                             */
  if (mode & D_RESET) {
    osl = (vid == 3) ? l->osl : l->mosl;
    if (vid == 0)
      gs = l->gsx;
    else if (vid == 1)
      gs = l->gsy;
    else if (vid == 2)
      gs = l->gsz;
    else if (vid == 3) {
      gs = l->gst;
    }
    for (k = 0; k < l->nz; k++) {
      if (!(l->dol[k])) continue;
      if (osl & (L_LINEAR|L_CUBIC|L_LSQUAD|L_LSLIN)) {
	gs[k]->rebuild(gs[k]->interpolator, d[k]->points);
      } else {
	/* Not properly function: uses ht_delete() in nnpi_interpolate()
	   if (osl & (L_SIB|L_NONSIB))
	   nnhpi_destroy_weights(gs[k]->interpolator);
	*/
	for (cc = 0; cc < d[k]->npoints; cc++)
	  gs[k]->rebuild(gs[k]->interpolator, &d[k]->points[cc]);
      }
    }
  }
}

/* END delaunay_fill()                                               */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Sets the boundary condition for the sediment layer                */
/*-------------------------------------------------------------------*/
void delaunay_sed(lagrange_t *l, int vid, double **v)
{
  int cc, c, cb, kb;
  int nz = l->nz-1;
  double val;
  double s = -1.0;
  delaunay *d = l->d[nz];
  int **map = l->u2d;

  if (vid == 3) {
    s = 1.0;
    if (l->nvec) {
      d = l->td[nz];
      map = l->t2d;
    }
  }

  /* Fill the Delaunay values. d is the surface Delaunay structure.  */
  for (cc = 0; cc < d->npoints; cc++) {
    c = abs(map[nz][cc]);
    kb = l->bot[c];
    cb = abs(map[kb][cc]);
    val = s * v[kb][cb];
    l->d[0]->points[cc].v[vid] = val;
    l->d[0]->points[cc].z = val;
  }
}

/* END delaunay_sed()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initializes a transport grid_spec structure                       */
/*-------------------------------------------------------------------*/
void grid_spec_init_delaunay(lagrange_t *l, int vid)
{
  int k, kk, nz = l->nz;
  GRID_SPECS **gs;

  /*-----------------------------------------------------------------*/
  /* Allocate and initialize                                         */
  gs = (GRID_SPECS **)calloc(nz+1, sizeof(GRID_SPECS *));

  for (k = 0; k < nz; k++) {

    gs[k] = grid_spec_create();

    if (!(l->dol[k])) continue;

    if (vid == 3) {
      if (l->nvec)
	grid_interp_init_t(gs[k], l->td[k], l->trasr, vid);
      else
	grid_interp_init_t(gs[k], l->d[k], l->trasr, vid);
    } else {
      grid_interp_init_t(gs[k], l->d[k], l->momsr, vid);
    }

    /* Set the depth levels */
    gs[k]->nz = nz;
    gs[k]->z = d_alloc_1d(nz);
    for (kk = 1; kk < nz; kk++) {
      gs[k]->z[kk] = l->layers[kk];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the correct grid_spec structure                             */
  if (vid == 0)
    l->gsx = gs;
  else if (vid == 1)
    l->gsy = gs;
  else if (vid == 2)
    l->gsz = gs;
  else if (vid == 3)
    l->gst = gs;
}

/* END grid_spec_init_delaunay()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Rebuild the Delaunay point trianges.                              */
/* The Delaunay structure contains the points in adjacent trianges   */
/* emanating from a point. This can be used as the neighbourhood for */
/* least squares interpolation. Some meshes (e.g. COMPAS) have a     */
/* dual representing a high quality Delaunay triangulation, and      */
/* rather than have triangle() create the point trianges we over-    */
/* ride them using the neighbours of the mesh.                       */
/*-------------------------------------------------------------------*/
void rebuild_pt(meshs_t *mesh, lagrange_t *l, delaunay **d, npts_t **npt)
{
  int cc, k;
  int nz = l->nz;

  /* Allocate                                                        */
  npt = (npts_t **)calloc(nz, sizeof(npts_t *));
  memset(npt, 0, nz * sizeof(npts_t));

  /* Free the point_triangles array if necessary                     */
  for (k = 0; k < nz; k++) {
    if (l->d[k] != NULL) {
      npt[k] = npt_create(d[k]->npoints);
      if (d[k]->point_triangles != NULL) {
	for (cc = 0; cc < d[k]->npoints; ++cc)
	  if (d[k]->point_triangles[cc] != NULL)
	    free(d[k]->point_triangles[cc]);
	free(d[k]->point_triangles);
      }
      if (d[k]->n_point_triangles != NULL)
	free(d[k]->n_point_triangles);
    }
  }

  /* Tracer point triangles                                          */
  reset_interp_points(mesh, l, npt);
}

/* END rebuild_pt()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets the points in the Delaunay point trianges.                 */
/*-------------------------------------------------------------------*/
void reset_interp_points(meshs_t *mesh, 
			 lagrange_t *l,
			 npts_t **npt
			 )
{
  int cc, c, cn, k, i, j, n;
  int nz = l->nz;
  int **map = l->u2d;
  int **rmap;
  delaunay **d = l->d;

  /*-----------------------------------------------------------------*/
  /* Set maps                                                        */
  if (l->nvec) {
    map = l->t2d;
    rmap = i_alloc_2d(l->nvec + 1, nz + 1);
    d = l->td;
  } else
    rmap = i_alloc_2d(mesh->ns2 + 1, nz + 1);

  /*-----------------------------------------------------------------*/
  /* Make the reverse map; Delaunay points to ugrid coordinates      */
  for (k = 0; k < nz; k++) {
    if (!(l->dol[k])) continue;
    for(cc = 0; cc < d[k]->npoints; cc++) {
      c = abs(map[k][cc]);
      rmap[k][c] = cc;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for (k = 0; k < nz; k++) {
    if (!(l->dol[k])) continue;
    for(cc = 0; cc < d[k]->npoints; cc++) {
      c = abs(map[k][cc]);

      /* Set the flag in the Delaunay structure                      */
      d[k]->ptf = 1;

      /* Cell centre                                                 */
      npt[k]->npt[cc] = 1;

      /* Neighbour cells to include                                  */
      for (n = 1; n <= mesh->npe[c]; n++) {
	cn = mesh->neic[n][c];
	if (cn > 0) {
	  npt[k]->npt[cc]++;
	}
      }

      /* Allocate                                                    */
      npt[k]->pt[cc] = malloc(npt[k]->npt[cc] * sizeof(int));

      /* Fill the locations of Delaunay indices                      */
      /* Wet cells to include                                        */
      j = 0;
      npt[k]->pt[cc][j++] = c;
      for (i = 1; i <= mesh->npe[c]; i++) {
	cn = mesh->neic[i][c];
	if (cn > 0) {
	  npt[k]->pt[cc][j++] = rmap[k][cn];
	}
      }
    }
  }

  /* Set the sediments                                               */
  k = nz-1;
  d[0]->ptf = 1;
  npt[0]->npt = npt[k]->npt;
  npt[0]->pt = npt[k]->pt;

  i_free_2d(rmap);
}

/* END reset_interp_points()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets the point triangles in the delaunay structure with tracer  */
/* or momentum point trianges. The point trianges are lists of       */
/* points in the delaunay triangulation used in interpolations       */
/* around any given triangulation point, and are dependent on the    */
/* type of interpolation used.                                       */
/*-------------------------------------------------------------------*/
static void reset_npt(lagrange_t *l, delaunay **d, npts_t **npt) 
{
  int k;

  for (k = 0; k <= l->nz; k++) {
    if (d[k] != NULL) {
      d[k]->n_point_triangles = npt[k]->npt;
      d[k]->point_triangles = npt[k]->pt;
    }
  }
}

/* END reset_npt()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a point triangle structure                                */
/*-------------------------------------------------------------------*/
static npts_t* npt_create(int n)
{
  npts_t* npt = malloc(sizeof(npts_t));

  npt->npt = calloc(n, sizeof(int));
  memset(npt->npt, 0, n * sizeof(int));
  npt->pt = malloc(n * sizeof(int*));

  return npt;
}

/* END npt_create()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the          */
/* streamline source locations.                                      */
/* Assume nz = u->kcentresize+1 as specified in the ugrid file.      */
/* A vector kvec[0:nk-1] of size nk is supplied specifying the       */
/* layers for streamline tracing. The values of this vecor may range */
/* from 1 to nz-1 (note that layer 0 is the sediment and should not  */
/* be included for streamline tracing - only used a bottom BC).      */
/* A vector vec[nk][nvec] for each of those layers is supplied       */
/* specifying the model data index (1:mesh->ns2) of the destination  */
/* points for streamline tracing.                                    */
/* Note; if 2D tracing is required, then nk = 1, nk[0] = nz-1 and    */
/* w = 0.                                                            */
/* This algorithm assumes all variables have a no-gradient condition */
/* imposed above the free surface within the input data.             */
/* If mode = I_FIRST, then the velocity data is used to populate     */
/* the grid_spec structures and weights are rebuilt.                 */
/* If nu = nv = nw = NULL, then it is assumed the grid_spec has been */
/* previously populated, and velocities are interpolated at the      */
/* first iteration rather than using the velocity (u,v,w).           */
/* If eta = NULL, then the streamline is not trunctaed to the free   */
/* surface.                                                          */
/* The streamline is always truncated to the bottom depth.           */
/*-------------------------------------------------------------------*/
void semi_lagrange_source(ugrid_t *ug,  /* UGRID information         */
			  double **u,   /* u velocity for each layer */
			  double **v,   /* v velocity for each layer */
			  double **w,   /* w velocity for each layer */
			  double *eta,  /* Surface elevation         */
			  double dtin,  /* Time-step                 */
			  int nk,
			  int *kvec,
			  int *nvec,
			  int **vec,
			  int mode
			  )
{
  lagrange_t *l = ug->l;        /* Lagrange information              */
  meshs_t *m = ug->meshs;       /* Mesh information                  */
  delaunay **d = l->d;          /* Delaunay information              */
  double dt;                    /* Sub-time step                     */
  double time_left;             /* Number of sub-timesteps per dt    */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  int *cl;                      /* Cell index of source location     */
  int *ck;                      /* Cell layer of source location     */
  double *cx, *cy, *cz;         /* Source location                   */
  double un, vn, wn;            /* Interpolated velocity             */
  double p, q, r;               /* Sub-steps                         */
  double d1, d2;                /* Dummies                           */
  int nz = l->nz;               /* Number of layers                  */
  int cc, c, cn, cs, n, i;      /* Counters                          */
  int k, kk, ks;                /* Vertical counters                 */
  int first;                    /* First iteration flag              */
  double m2deg = 1.0 / (60.0 * 1852.0);  /* Degrees to m conversion  */
  int has_proj = (strlen(ug->projection) > 0); /* Projection flag    */
  int is_geog = has_proj && (strcasecmp(ug->projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int onestep = 0;  /* Find streamline in one step                   */
  if (!is_geog) m2deg = 1.0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  for (kk = 0; kk < nk; kk++) {
    k = kvec[kk];
    if (l->dol[k]) {
      /* Set the initial streamline locations                        */
      for (cc = 0; cc < nvec[k]; cc++) {
	c = vec[k][cc];
	l->cl[k][c] = c;
	l->ck[k][c] = k;
	l->cx[k][c] = m->xloc[m->eloc[0][c][0]];
	l->cy[k][c] = m->yloc[m->eloc[0][c][0]];
	l->cz[k][c] = l->layers[k];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Initialize the interpolation structure                          */
  if (mode & I_FIRST) {
    delaunay_fill(l, 0, u, D_FILL|D_RESET);
    delaunay_fill(l, 1, v, D_FILL|D_RESET);
    if (!(l->options & I_2D))
      delaunay_fill(l, 2, w, D_FILL|D_RESET);
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin                         */
  for (kk = 0; kk < nk; kk++) {
    k = kvec[kk];
    cl = l->cl[k];
    ck = l->ck[k];
    cx = l->cx[k];
    cy = l->cy[k];
    cz = l->cz[k];
    if (!(l->dol[k])) continue;

    for (cc = 0; cc < nvec[k]; cc++) {
      first = 0;
      c = vec[k][cc];
      time_left = dtin;
      if (u == NULL && v == NULL && w == NULL)
	first = 1;
      else {
	un = u[k][c];
	vn = v[k][c];
	wn = w[k][c];
      }
      n = 0;

      while (time_left > 0) {

	double un = max(fabs(un), SMALL);
	double vn = max(fabs(vn), SMALL);
	double wn = max(fabs(wn), SMALL);
	cs = cl[c];
	ks = abs(ck[c]);

	/* Get the velocties at the origin                           */
	if (!onestep && first) {
	  un = lag_interp(l, l->gsx, cx[c], cy[c], cz[c], cl[c], ck[c], c, 0);
	  vn = lag_interp(l, l->gsy, cx[c], cy[c], cz[c], cl[c], ck[c], c, 1);
	  if (!(l->options & I_2D))
	    wn = lag_interp(l, l->gsz, cx[c], cy[c], cz[c], cl[c], ck[c], c, 2);
	}
	first = 1;

	/* Get the sub-timestep for the next integration. Streamlines  */
	/* cannot cross more than one cell in one sub-timestep.        */
	d1 = d2 = 0.0;
	for (i = 1; i <= m->npe[cs]; i++) {
	  if ((cn = m->neic[i][cs])) {
	    p = m->xloc[m->eloc[0][cs][0]] - m->xloc[m->eloc[0][cn][0]];
	    q = m->yloc[m->eloc[0][cs][0]] - m->yloc[m->eloc[0][cn][0]];
	    d1 = sqrt(p * p + q * q) / m2deg;
	    d2 += 1.0;
	  }
	  d1 /= d2;
	}
	p = (d1) ? fabs(d1 / sqrt(un * un + vn * vn)) : dt;
	if (l->options & I_2D) {
	  wn = 0.0;
	  r = p;
	} else {
	  q = l->dzz[ks][cs];
	  r = (q) ? fabs(q / wn) : dt;
	}
	dt = min(SCALE * p, SCALE * r);
	dt = (onestep) ? time_left : min(dt, time_left);

	/* Get the new location of the streamline origin             */
	cl[c] = lag_pos(m, l, &cl[c], &ck[c], un, vn, wn, &cx[c], &cy[c], &cz[c], dt, is_geog);

	/* Truncate                                                  */
	if (cz[c] < ug->bathy[cl[c]]) cz[c] = ug->bathy[cl[c]];
	if (eta != NULL && cz[c] > eta[cl[c]]) cz[c] = eta[cl[c]];

	/* Get the number of sub-timesteps used                      */
	time_left = time_left - dt;
	n++;
      }
    }
  }
}

/* END semi_lagrange_source()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This finds the source of a streamline if the destiunation is      */
/* suppled; (xin, yin, zin).                                         */
/* Note; if 2D tracing is required, then zin = 0 and w = 0.          */
/* This algorithm assumes all variables have a no-gradient condition */
/* imposed above the free surface within the input data.             */
/* If mode = I_FIRST, then the velocity data is used to populate     */
/* the grid_spec structures and weights are rebuilt.                 */
/* If nu = nv = nw = NULL, then it is assumed the grid_spec has been */
/* previously populated, and velocities are interpolated at the      */
/* first iteration rather than using the velocity (u,v,w).           */
/* If eta = NULL, then the streamline is not trunctaed to the free   */
/* surface.                                                          */
/* The streamline is always truncated to the bottom depth.           */
/*-------------------------------------------------------------------*/
int semi_lagrange_atxy(ugrid_t *ug,     /* UGRID information         */
		       double **u,      /* u velocity for each layer */
		       double **v,      /* v velocity for each layer */
		       double **w,      /* w velocity for each layer */
		       double *eta,     /* Surface elevation         */
		       double dtin,     /* Time-step                 */
		       double *xin,     /* x location                */
		       double *yin,     /* y location                */
		       double *zin,     /* z location                */
		       int mode
		       )
{
  lagrange_t *l = ug->l;        /* Lagrange information              */
  meshs_t *m = ug->meshs;       /* Mesh information                  */
  delaunay **d = l->d;          /* Delaunay information              */
  double dt;                    /* Sub-time step                     */
  double time_left;             /* Number of sub-timesteps per dt    */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  int *cl;                      /* Cell index of source location     */
  int *ck;                      /* Cell layer of source location     */
  double *cx, *cy, *cz;         /* Source location                   */
  double un, vn, wn;            /* Interpolated velocity             */
  double p, q, r;               /* Sub-steps                         */
  double d1, d2;                /* Dummies                           */
  int nz = l->nz;               /* Number of layers                  */
  int cc, c, cn, cs, n;         /* Counters                          */
  int first = 0;                /* First iteration flag              */
  int i, j, k, kk, ks;
  double m2deg = 1.0 / (60.0 * 1852.0);
  int has_proj = (strlen(ug->projection) > 0);
  int is_geog = has_proj && (strcasecmp(ug->projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int onestep = 0;  /* Find streamline in one step                   */
  if (!is_geog) m2deg = 1.0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  cc = 0;
  lag_index(m, l, *xin, *yin, *zin, &c, &k, -1);
  l->cx[k][c] = *xin;
  l->cy[k][c] = *yin;
  l->cz[k][c] = *zin;
  l->cl[k][c] = c;
  l->ck[k][c] = k;

  /*-----------------------------------------------------------------*/
  /* Initialize the interpolation structure                          */
  if (mode & I_FIRST) {
    delaunay_fill(l, 0, u, D_FILL|D_RESET);
    delaunay_fill(l, 1, v, D_FILL|D_RESET);
    if (!(l->options & I_2D))
      delaunay_fill(l, 2, w, D_FILL|D_RESET);
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin                         */
  cl = l->cl[k];
  ck = l->ck[k];
  cx = l->cx[k];
  cy = l->cy[k];
  cz = l->cz[k];
  time_left = dtin;
  if (u == NULL && v == NULL && w == NULL)
    first = 1;
  else {
    un = u[k][c];
    vn = v[k][c];
    wn = w[k][c];
  }
  n = 0;

  while (time_left > 0) {

    double ua = max(fabs(un), SMALL);
    double va = max(fabs(vn), SMALL);
    double wa = max(fabs(wn), SMALL);
    cs = cl[c];
    ks = abs(ck[c]);

    /* Get the velocties at the origin                               */
    if (!onestep && first) {
      un = lag_interp(l, l->gsx, cx[c], cy[c], cz[c], cl[c], ck[c], c, 0);
      vn = lag_interp(l, l->gsy, cx[c], cy[c], cz[c], cl[c], ck[c], c, 1);
      if (!(l->options & I_2D))
	wn = lag_interp(l, l->gsz, cx[c], cy[c], cz[c], cl[c], ck[c], c, 2);
    }
    first = 1;

    /* Get the sub-timestep for the next integration. Streamlines    */
    /* cannot cross more than one cell in one sub-timestep.          */
    d1 = d2 = 0.0;
    for (i = 1; i <= m->npe[cs]; i++) {
      if ((cn = m->neic[i][cs])) {
	p = m->xloc[m->eloc[0][cs][0]] - m->xloc[m->eloc[0][cn][0]];
	q = m->yloc[m->eloc[0][cs][0]] - m->yloc[m->eloc[0][cn][0]];
	d1 = sqrt(p * p + q * q) / m2deg;
	d2 += 1.0;
      }
      d1 /= d2;
    }
    p = (d1) ? fabs(d1 / sqrt(ua * ua + va * va)) : dt;
    if (l->options & I_2D) {
      wn = 0.0;
      r = p;
    } else {
      q = l->dzz[ks][cs];
      r = (q) ? fabs(q / wa) : dt;
    }
    dt = min(SCALE * p, SCALE * r);
    dt = (onestep) ? time_left : min(dt, time_left);

    /* Get the new location of the streamline origin                 */
    cl[c] = lag_pos(m, l, &cl[c], &ck[c], un, vn, wn, &cx[c], &cy[c], &cz[c], dt, is_geog);

    /* Truncate                                                      */
    d1 = l->layers[l->bot[cl[c]]];
    if (cz[c] < d1) cz[c] = d1;
    if (eta != NULL && cz[c] > eta[cl[c]]) cz[c] = eta[cl[c]];

    /* Get the number of sub-timesteps used                          */
    time_left = time_left - dt;
    n++;
  }
  *xin = cx[c];
  *yin = cy[c];
  *zin = cz[c];
  return(cl[c]);
}

/* END semi_lagrange_atxy()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates tracer values at the streamline source location.     */
/* If mode = I_FIRST, then the tracer data is used to populate the   */
/* grid_spec structures and weights are rebuilt.                     */
/* If a valid layer and location of the streamline origin (i.e. the  */
/* destinationin point) in the mesh is supplied (ki >= 0 and ci >=0, */
/* then tracer is interpolated only at this location.                */
/* If a vector of ponts was submitted to lagrange_init(), then the   */
/* tracers are interpolated at the source points of streamlines      */
/* originating from those points. Otherwise the interpolation is     */
/* performed at the source location of all streamlines in the mesh.  */
/*-------------------------------------------------------------------*/
double semi_lagrange_tr(ugrid_t *ug, /* UGRID information            */
			double **tr, /* Tracer array                 */
			int ki,      /* Layer to interpolate         */
			int ci,      /* Cell index to interpolte     */
			int mode     /* Options                      */
			)
{
  lagrange_t *l = ug->l;        /* Lagrange information              */
  meshs_t *m = ug->meshs;       /* Mesh information                  */
  delaunay **d = l->d;          /* Delaunay information              */
  int **map = l->u2d;           /* Delaunay to mesh index map        */
  int cc, c, k;                 /* Counters                          */
  double val;                   /* Return value                      */

  /* Set pointers                                                    */
  if (l->nvec) {
    d = l->td;
    map = l->t2d;
  }

  /* Update the tracer values at cell centres in the interpolation   */
  /* structure.                                                      */
  if (mode & I_FIRST) {
    delaunay_fill(l, 3, tr, D_FILL|D_RESET);
  }

  /* Return the interpolated value at layer ki, cell ci              */
  if (ki >= 0 && ci >= 0) {
    if (!(l->dol[ki])) return(0.0);
    k = ki;
    c = ci;
    val = lag_interp(l, l->gst, l->cx[k][c], l->cy[k][c], l->cz[k][c], l->cl[k][c], l->ck[k][c], c, 3);
    return(val);
  }

  /* Do the interpolation for all centres                            */
  for (k = 1; k < l->nz; k++) {
    if (!(l->dol[ki])) continue;
    for (cc = 0; cc < d[k]->npoints; cc++) {
      c = abs(map[k][cc]);
      tr[k][c] = lag_interp(l, l->gst, l->cx[k][c], l->cy[k][c], l->cz[k][c], l->cl[k][c], l->ck[k][c], c, 3);
    }
  }
}

/* END semi_lagrange_tr()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates in 3D on an unstructured grid given an (x,y,z)       */
/* location. The cell, c, where (x,y,z) resides must also be         */
/* supplied.                                                         */
/*-------------------------------------------------------------------*/
static double lag_interp(lagrange_t *l, GRID_SPECS **gs, double x, double y, double z, int c, int k, int co, int vid)
{
  int **map = (l->nvec) ? l->t2d : l->u2d;
  int osl = (vid == 3) ? l->osl : l->mosl;
  int k2 = abs(k) + 1;
  int k1 = (k < 0) ? 0 : k;
  double v, v1, v2, d;
  double dzz;

  if (k2 > l->nz-1) k2 = l->nz-1; /* No-gradient above surface       */
  if (l->options & I_2D) k2 = l->nz - 1;
  dzz = l->dzz[k2][c];            /* cellz thickness                 */

  /* For least squares, the interpolation is centered on the cell    */
  /* centre (i.e. co). This is mapped to the Delaunay index and      */
  /* passed to the interpolation scheme directly in gs->id (i.e.     */
  /* there is no need to find the triangle (x,y) resides in within   */
  /* the interpolation scheme via delaunay_xytoi_ng().               */
  if (osl & (L_LSQUAD|L_LSLIN)) {
    gs[k2]->id = abs(map[k2][c]);
  }
  v2 = grid_interp_on_point(gs[k2], x, y);

  if(isnan(v2)){
    quit("lag_interp: NaN c=%d cl=%d(k%d) %d %f %f %f\n",co, c, k2, vid, x, y, v2);
  }
  if (l->options & I_2D) return(v2);

  /* Evaluate in the layer with a cellz level less than cz. If this  */
  /* is the bottom layer, then the ghost layer below the bottom will */
  /* take care of the vertical interpolation.                        */
  if (k1 == 0) {
    v1 = v2;
    if (vid < 3) {
      v1 *= -1.0;
    }
  } else {
    if (osl & (L_LSQUAD|L_LSLIN)) gs[k1]->id = abs(map[k1][c]);
    v1 = grid_interp_on_point(gs[k1], x, y);
  }

  /* Interpolate vertically                                          */
  d = z - l->layers[k];
  v = d * (v2 - v1) / dzz + v1;

  return(v);

}

/* END lag_interp()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the new location of the streamline origin. If this lands */
/* in a ghost cell, then shift the streamline in the wet interior    */
/* and re-compute a new location.                                    */
/* Destination cell is the origin of the streamline.                 */
/* Source cell is the end point of the streamline track.             */
/* Returns a new source cell after moving a distance increment.      */
/* The new geographic location is returned in cx, cy and cz.         */
/*-------------------------------------------------------------------*/
int lag_pos(meshs_t *m,         /* Mesh information                  */
	    lagrange_t *l,      /* Lagrange information              */
	    int *ci,            /* Current source cell               */
	    int *ki,            /* Current source layer              */
	    double u,           /* East velocity                     */
	    double v,           /* North velocity                    */
	    double w,           /* Vertical velocity                 */
            double *cx,         /* x location of streamline          */
	    double *cy,         /* y location of streamline          */
	    double *cz,         /* z location of streamline          */
	    double dt,          /* Time step                         */
	    int is_geog         /* Geographic flag                   */
	    )
{
  double xin = *cx;
  double yin = *cy;
  double zin = *cz;
  double dist, dir;
  double slat, slon;
  double sinth, costh;
  double speed;
  double eps = 1e-8;
  double m2deg = 1.0 / (60.0 * 1852.0);
  double fact = 0.99;
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int cc, c, j;
  point p;
  delaunay *d = l->d[l->nz-1];

  /* Get the new horizontal location                                 */
  if (!is_geog) m2deg = 1.0;
  if (sodanos) {
    speed = sqrt(u * u + v * v);
    dist = speed * dt;
    dir =  PI / 2.0 - atan2(v, u);  /* deg T                         */
    geod_fwd_sodanos(DEG2RAD(xin), DEG2RAD(yin), dir, dist, RADIUS, ECC, &slon, &slat);
    slon = RAD2DEG(slon);
    slat = RAD2DEG(slat);
  } else {
    /* Get the new distances                                         */
    speed = sqrt(u * u + v * v);
    dist = speed * dt * m2deg;

    dir =  atan2(v, u);
    /* Update the origin location                                    */
    sinth = sin(dir);
    costh = cos(dir);
    slon = xin - dist * costh;
    slat = yin - dist * sinth;
  }

  /* Update the streamline source location                           */
  *cx = slon;
  *cy = slat;
  *cz -= w * dt;

  /* Get the cell the new horizontal location resides in, i.e. at    */
  /* the same depth as that of the input location, c.                */
  c = *ci;
  lag_index(m, l, *cx, *cy, *cz, &c, ki, -1);

  /* If the source lies outside the mesh, and an index can't be      */
  /* found, then compute the distance to the perimeter and set the   */
  /* source to lie just inside the domain.                           */
  /*if (!poly_contains_point(m->ppl, slon, slat)) {*/
  if (c == 0) {
    /*printf("ghost %d (%f %f) : (%f %f) %f\n",*ci,xin,yin,slon,slat,dist);*/
    for (j = 0; j < m->npedge[*ci]; j+=2) {
      point s, d;
      double xi, yi, ndt;
      s.x = slon; s.y = slat;
      d.x = xin; d.y = yin;
      if (find_intersect(s, d, m->pedge[*ci][j], m->pedge[*ci][j+1], &xi, &yi)) {
	double dx = xin - xi;
	double dy = yin - yi;
	dist = fact * sqrt(dx * dx + dy * dy);
	ndt = (speed) ? dist / speed : 0.0;
	c = lag_pos(m, l, ci, ki, u, v, w, &xin, &yin, &zin, ndt, is_geog);
	*cx = xin; *cy = yin; *cz = zin;
	break;
      }
    }
  }

  /* If the source index still can't be found, then set the source   */
  /* to the last valid cell.                                         */
  if (c == 0) c = *ci;

  return(c);
}

/* END lag_pos()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the index of the polygon that (x,y,z) resides in            */
/*-------------------------------------------------------------------*/
void lag_index(meshs_t *m,
	       lagrange_t *l,
	       double x, 
	       double y, 
	       double z, 
	       int *c, 
	       int *k,
	       int first
	       )
{
  int cc;
  point p;
  delaunay *d = m->d;

  /* Get the cell the new horizontal location resides in, i.e. at    */
  /* the same depth as that of the input location, c.                */
  p.x = x;
  p.y = y;
  if (first >= 0) d->first_id = first;
  cc = delaunay_xytoi(m->d, &p, d->first_id);
  *c = m->tri2c[cc];

  /**c = abs(l->u2d[l->nz-1][cc]);*/

  /* Evaluate the layer above the bottom that cz resides in, with a  */
  /* layer level less than cz.                                       */
  *k = l->nz-1;
  if (z > l->layers[l->nz-1])    /* Set to the surface layer         */
    *k = l->nz-1;
  else
    while (*k >= 0 && z < l->layers[*k]) *k -= 1;
  if (*k < l->bot[*c]) *k = -*k;  /* Set to the sediment layer       */
}

/* END lag_index()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 if a line from point (xn, yn) to (xs, ys) intersects    */
/* edge j. If so, the location of the intersection are returned in   */
/* (xi, yi).                                                         */
/* This version is mre formally correct.                             */
/*-------------------------------------------------------------------*/
int find_intersect(point s,          /* source point                 */
		   point d,          /* destination point            */
		   point e1,         /* edge point 1                 */
		   point e2,         /* edge point 2                 */
		   double *xi,       /* x location of intersection   */
		   double *yi        /* y location of intersection   */
		   )
{
  double x1, y1, x2, y2, s1, s2;
  int found = 0;
  double v1x = e1.x;
  double v1y = e1.y;
  double v2x = e2.x;
  double v2y = e2.y;
  double dx = fabs(v1x - v2x);
  double dy = fabs(v1y - v2y);

  /* Get the location where edge and streamline intersect. Get        */
  /* equations for line between source to destination, and between    */
  /* the verticies, make them equal and solve for x and y.            */
  /* Souce-destination eqn: y=(x-d.x)s2+d.y                           */
  /* Vertices eqn: y = (x-v1x)s1 + v1y                                */
  /* Special treatment when v1x=v1y (dx=0) and s.x=d.x                */
  if (dx) {
    s1 = (v2y - v1y) / (v2x - v1x);
    if (s.x == d.x) 
      *xi = d.x;
    else {
      s2 = (s.y - d.y) / (s.x - d.x);
      *xi = (s2 * d.x - s1 * v1x + v1y - d.y) / (s2 - s1);
    }
    *yi = s1 * (*xi - v1x) + v1y;
  } else {
    s1 = (v2x - v1x) / (v2y - v1y);
    if(s.y == d.y)
      *yi = d.y;
    else {
      s2 = (s.x - d.x) / (s.y - d.y);
      *yi = (s1 * v1y - s2 * d.y + d.x - v1x) / (s1 - s2);
    }
    *xi = s1 * (*yi - v1y) + v1x;
  }
  /* If the intersection lies between the bounds of the vertices then */
  /*the streamline crosses this edge.                                 */
  if (*xi >= min(v1x, v2x) &&
      *xi <= max(v1x, v2x) &&
      *yi >= min(v1y, v2y) &&
      *yi <= max(v1y, v2y)) {
    if (*xi >= min(d.x, s.x) &&
	*xi <= max(d.x, s.x) &&
	*yi >= min(d.y, s.y) &&
	*yi <= max(d.y, s.y)) {
      return(1);
    }
  }
  return(found);
}

/* END find_intersect()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the ugrid and Lagrange structures, and fills the          */
/* triangulation with global data.                                   */
/*-------------------------------------------------------------------*/
void test_lagrange(char *infile, char *trname)
{
  ugrid_t *ug = NULL;
  lagrange_t *l;
  int fid, vid;
  int ncerr;
  int ns2;
  int dump = 0;
  int c, cc, c2, cb, k, kk;
  int nz;
  int nk, *kvec, *nvec, **vec;
  double **nu, **nv, **nw, **ntr, *eta;
  double dt = 3600.0;
  double xin, yin, zin;
  double wtop, wbot;
  size_t start[4];
  size_t count[4];
  int verbose = 1;

  /* Set up the ugrid structure based on infile                      */
  if (ug == NULL) {
    ug = ugrid_open(infile, I_GLOB);
    /* Set up the Lagrange structure                                 */
    lagrange_init(ug, 0, NULL);
  }
  l = ug->l;
  nz = l->nz;
  ns2 = ug->ns2;
  if (verbose) printf("UGRID and Lagrange structures initialized OK\n");

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(infile, NC_NOWRITE, &fid)) != NC_NOERR) {
    warn("Can't find input file %s\n", infile);
    quit((char *)nc_strerror(ncerr));
  }

  /* Allocate                                                        */
  nu = d_alloc_2d(ns2+1, nz);
  nv = d_alloc_2d(ns2+1, nz);
  nw = d_alloc_2d(ns2+1, nz);
  nw = d_alloc_2d(ns2+1, nz);
  ntr = d_alloc_2d(ns2+1, nz);
  eta = d_alloc_1d(ns2+1);

  /* Read the data from file. Note the first index of data in        */
  /* in UGRID may be eith 0 or 1, specified by start_index in the    */
  /* global attributes (usually 0). The data in the Lagrange         */
  /* structure is expected to index from 1 to ns2. This re-packing   */
  /* is handled in the ugrid_read routines.                          */
  ugrid_read_3d(ug, fid, "u", nu, dump);
  ugrid_read_3d(ug, fid, "v", nu, dump);
  ugrid_read_3d(ug, fid, "w", nw, dump);
  ugrid_read_3d(ug, fid, trname, ntr, dump);
  ugrid_read_2d(ug, fid, "eta", eta, dump);

  /* Center the vertical velocity                                    */
  for (c = 1; c <= ns2; c++) {
    k = l->bot[c];
    wbot = nw[l->bot[c]][c];
    for (k = l->bot[c]+1; k < nz; k++) {
      wtop = (k == l->nz-1) ? wbot : nw[k][c];
      nw[k][c] = 0.5 * (wbot + wtop);
      wbot = wtop;
    }
  }
  if (verbose) printf("Data read OK from file %s\n", infile);

  /* Fill the Delaunay triangulation points                          */
  delaunay_fill(l, 0, nu, D_FILL|D_RESET);
  delaunay_fill(l, 1, nv, D_FILL|D_RESET);
  delaunay_fill(l, 2, nw, D_FILL|D_RESET);
  delaunay_fill(l, 3, ntr, D_FILL|D_RESET);

  /* Allocate the points to process vectors                          */
  kvec = i_alloc_1d(nz);
  nvec = i_alloc_1d(nz);
  vec = i_alloc_2d(ns2, nz);
  /* Set up the points to process                                    */
  memset(nvec, 0, nz * sizeof(int));
  nk = 0;
  for (c = 1; c <= ns2; c++) {
    k = l->bot[c];
    for (k = l->bot[c]; k < nz; k++) {
      vec[k][nvec[k]++] = c;
    }
  }
  /* Set the layers to process                                       */
  nk = 0;
  for (k = 0; k < nz; k++) {
    if (nvec[k]) {
      kvec[nk++] = k;
    }
  }

  /* Track the stream to its source. Streamline source is stored in  */
  /* (u->l->cx, u->l->cy, u->l->cz).                                 */
  semi_lagrange_source(ug, NULL, NULL, NULL, eta, dt, nk, kvec, nvec, vec, 0);
  if (verbose) printf("Streamlines tracked OK\n");

  /* Update the tracers                                              */
  for (c = 1; c <= ns2; c++) {
    for (k = l->bot[c]; k < nz; k++) {
      ntr[k][c] = semi_lagrange_tr(ug, NULL, k, c, 0);
      if (verbose == 2) printf("k=%d c=%d cx=%f cy=%f cz=%f : %f\n",k, c, l->cx[k][c],
			       l->cy[k][c], l->cz[k][c], ntr[k][c]);
    }
  }
  if (verbose) printf("Tracer %s updated OK\n", trname);

  /* Streamline tracking from a single (x,y,z) location              */
  /*
  xin = u->cellx[1];
  yin = u->celly[1];
  zin = 0.0;
  c = semi_lagrange_atxy(ug, nu, nv, nw, eta, dt, &xin, &yin, &zin, 0);
  if (verbose) printf("Loaction (%f %f %f) tracked to (%f %f %f)\n", 
		      u->cellx[1], u->cellx[1], 0.0, xin, yin, zin);
  */

  d_free_2d(nu);
  d_free_2d(nv);
  d_free_2d(nw);
  d_free_2d(ntr);
  d_free_1d(eta);
  lagrange_free(l);
  ugrid_free(ug);
  if (verbose) printf("Largrange tracking successful\n");
}

/* END test_lagrange()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads a 3D variable from file and packs into Lagrange indexing    */
/*-------------------------------------------------------------------*/
int ugrid_read_3d(ugrid_t *u, int id, char *name, double **p, int dump)
{
  size_t start[4];
  size_t count[4];
  int ns2 = u->nMesh2_face;
  int nz = u->kcentresize;
  int vid;
  int k, cc, c;
  double **in = d_alloc_2d(ns2, nz);

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = u->kcentresize;
  count[2] = u->nMesh2_face;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
    for (k = 0; k < nz; k++)
      memset((void *)p[k], 0, sizeof(double) * ns2);
    return(1);
  }

  nc_get_vara_double(id, vid, start, count, in[0]);
  for (k = 0; k < nz; k++) {
    for (cc = u->si; cc < u->nMesh2_face; cc++) {
      c = u->index[cc];
      p[k+1][c] = in[k][cc];
    }
  }
  d_free_2d(in);
  return(0);
}

/* END ugrid_read_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads a 2D variable from file and packs into Lagrange indexing    */
/*-------------------------------------------------------------------*/
int ugrid_read_2d(ugrid_t *u, int id, char *name, double *p, int dump)
{
  size_t start[4];
  size_t count[4];
  int ns2 = u->nMesh2_face;
  int vid;
  int k, cc, c;
  double *in = d_alloc_1d(ns2);

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = u->nMesh2_face;
  count[2] = 0;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
      memset((void *)p, 0, sizeof(double) * ns2);
    return(1);
  }

  nc_get_vara_double(id, vid, start, count, in);
  for (cc = u->si; cc < u->nMesh2_face; cc++) {
    c = u->index[cc];
    p[c] = in[cc];
  }
  d_free_1d(in);
  return(0);
}

/* END ugrid_read_2d()                                               */
/*-------------------------------------------------------------------*/
