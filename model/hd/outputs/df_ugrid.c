/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/df_ugrid.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  Output is written in cf compliant unstructured coordinates.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: df_ugrid.c 6807 2021-06-23 01:37:45Z her127 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define Filli 999999999
#define Filld 9.9692099683868690E36;


// from dumpfile.c
int get_nc_mode(dump_file_t *df);
void pack_ugrid1(int *hmap, int mapsize, double *var, double *pack, int oset);
void pack_ugrid2(geometry_t *geom, int *hmap, int mapsize, double *var, double **pack, int oset);
void pack_ugrids(int *map, int sednz, int mapsize, double **var, double **pack, int oset);
void pack_ugrid_i2(int *map, int size, int np, int **var, int **pack, int oset);
void pack_ugrid_i2a(int *map, int size, int np, int **var, int **pack, int oset);
void pack_ugrid_ri2(int *map, int size, int np, int **var, int **pack, int oset);
void pack_ugrid_ri2a(int *map, int size, int np, int **var, int **pack, int oset);
void vel_tan3D(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset);
void vel_tan2D(dump_data_t *dumpdata, double *u1, double *u2, double *p1, double *p2,  int oset);
void vel_cen3D(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset);
void vel_cen2D(dump_data_t *dumpdata, double *u1, double *u2, double *p1, double *p2,  int oset);
void vel_scalar(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset);

#define UGRID_ALL_VARS "u1av u2av uav vav wind1 wtop eta patm u1 u2 u v w dens dens_0 Kz Vz Cd u1kh topz "

void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
			      nc_type fptype, char *modulo);
static void df_ugrid_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_ugrid_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_ugrid_get_dimids(int cdfid, df_ugrid_var_t var, int *ns, int *zid);
static void df_ugrid_get_dimsizes(dump_file_t *df, df_ugrid_var_t var, size_t *sz, size_t *ns, size_t *nz);
static void write_mean_atts(dump_data_t *dumpdata, int fid);
int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref);
void set_dims(int dimf, int *dims, int d1, int d2);
void set_count(int dimf, size_t *count, int d1, int d2);

/*-------------------------------------------------------------------*/
/* Creates the mappings required for UGRID output                    */
/*-------------------------------------------------------------------*/
void create_ugrid_maps(dump_data_t *dumpdata)
{
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int c, cc, cn, e, ee, v, vv;
  int szm, npem;
  int i, j, *c2i, *c2j;

  if (dumpdata->npe != 0 && dumpdata->u1e != NULL) return;

  /*-----------------------------------------------------------------*/
  /* Set dimensions                                                  */
  dumpdata->npe = 4;
  npem = dumpdata->npe + 1;
  dumpdata->start_index = 0;
  dumpdata->face_dim = 0;
  dumpdata->edge_dim = 0;
  dumpdata->nface2 = geom->b2_t;
  dumpdata->nface3 = geom->b3_t;

  /*-----------------------------------------------------------------*/
  /* Get the sparse to Cartesian mappings                            */
  dumpdata->c2i = i_alloc_1d(geom->sgsiz);
  dumpdata->c2j = i_alloc_1d(geom->sgsiz);
  for (cc = 1; cc < geom->sgsiz; cc++) {
    dumpdata->c2i[cc] = -1;
    dumpdata->c2j[cc] = -1;
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    dumpdata->c2i[c] = i;
    dumpdata->c2j[c] = j;

    cn = geom->xp1[c];
    dumpdata->c2i[cn]=i+1;
    dumpdata->c2j[cn]=j;

    cn = geom->yp1[c];
    dumpdata->c2i[cn]=i;
    dumpdata->c2j[cn]=j+1;

    cn = geom->yp1[geom->xp1[c]];
    dumpdata->c2i[cn]=i+1;
    dumpdata->c2j[cn]=j+1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the edge mappings                                           */
  dumpdata->nedge2 = geom->n2_e1 + geom->n2_e2;
  dumpdata->nedge3 = geom->n3_e1 + geom->n3_e2;
  dumpdata->w2_e1 = i_alloc_1d(dumpdata->nedge3+1);
  dumpdata->e2c = i_alloc_2d(2, dumpdata->nedge3+1);
  dumpdata->c2e = i_alloc_2d(geom->sgsiz, npem);
  dumpdata->u1e = i_alloc_1d(geom->sgsiz);
  dumpdata->u2e = i_alloc_1d(geom->sgsiz);
  dumpdata->u1v = i_alloc_1d(geom->sgsiz);
  memset(dumpdata->u1e, 0, (geom->sgsiz) * sizeof(int));
  memset(dumpdata->u2e, 0, (geom->sgsiz) * sizeof(int));

  /* Map all the edges left (u1) of the centre                       */
  e = dumpdata->nu1 = 1;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    /* Grid vertex is bottom left corner                             */
    dumpdata->w2_e1[dumpdata->nu1++] = e;
    dumpdata->u1e[c] = e;
    dumpdata->c2e[1][c] = e;
    dumpdata->e2c[e][0] = geom->xm1[c];
    dumpdata->e2c[e][1] = c;
    e++;
    dumpdata->w2_e1[dumpdata->nu1++] = -e;
    dumpdata->u2e[c] = e;
    dumpdata->c2e[4][c] = e;
    dumpdata->e2c[e][0] = geom->ym1[c];
    dumpdata->e2c[e][1] = c;
    e++;
    /* Check for ghost vertices                                      */
    cn = geom->xp1[c];
    if (geom->wgst[cn]) {
      dumpdata->w2_e1[dumpdata->nu1++] = e;
      dumpdata->u1e[cn] = e;
      dumpdata->c2e[3][c] = e;
      dumpdata->e2c[e][0] = c;
      dumpdata->e2c[e][1] = cn;
      e++;
    }
    cn = geom->yp1[c];
    if (geom->wgst[cn]) {
      dumpdata->w2_e1[dumpdata->nu1++] = -e;
      dumpdata->u2e[cn] = e;
      dumpdata->c2e[2][c] = e;
      dumpdata->e2c[e][0] = c;
      dumpdata->e2c[e][1] = cn;
      e++;
    }
  }
  /* R_EDGE and F_EDGE open boundaries                               */
  for (i = 0; i < geom->nobc; i++) {
    open_bdrys_t *open = geom->open[i];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no2_e1; cc++) {
        c = open->obc_e1[cc];
	cn = open->obc_t[cc];
	if(!dumpdata->u1e[c]) {
	  dumpdata->w2_e1[dumpdata->nu1++] = e;
	  dumpdata->u1e[c] = e;
	  dumpdata->c2e[3][cn] = e;
	  dumpdata->e2c[e][0] = cn;
	  dumpdata->e2c[e][1] = c;
	  e++;
	}
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no2_e2; cc++) {
        c = open->obc_e2[cc];
	cn = open->obc_t[cc];
	if(!dumpdata->u2e[c]) {
	  dumpdata->w2_e1[dumpdata->nu1++] = -e;
	  dumpdata->u2e[c] = e;
	  dumpdata->c2e[2][cn] = e;
	  dumpdata->e2c[e][0] = cn;
	  dumpdata->e2c[e][1] = c;
	  e++;
	}
      }
    }
  }
  dumpdata->nu1S = dumpdata->nu1-1;
  dumpdata->nedge2 = e-1;

  /*  Get the 3D edges                                               */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    if (c == geom->m2d[c]) continue;
    /* Grid vertex is bottom left corner                             */
    if (!dumpdata->u1e[c]) {
      dumpdata->w2_e1[dumpdata->nu1++] = e;
      dumpdata->u1e[c] = e;
      dumpdata->c2e[1][c] = e;
      dumpdata->e2c[e][0] = geom->xm1[c];
      dumpdata->e2c[e][1] = c;
      e++;
    }

    if (!dumpdata->u2e[c]) {
      dumpdata->w2_e1[dumpdata->nu1++] = -e;
      dumpdata->u2e[c] = e;
      dumpdata->c2e[4][c] = e;
      dumpdata->e2c[e][0] = geom->ym1[c];
      dumpdata->e2c[e][1] = c;
      e++;
    }
    /* Check for ghost vertices                                      */
    cn = geom->xp1[c];
    if (geom->wgst[cn] && !dumpdata->u1e[cn]) {
      dumpdata->w2_e1[dumpdata->nu1++] = e;
      dumpdata->u1e[cn] = e;
      dumpdata->c2e[3][c] = e;
      dumpdata->e2c[e][0] = c;
      dumpdata->e2c[e][1] = cn;
      e++;
    }
    cn = geom->yp1[c];
    if (geom->wgst[cn] && !dumpdata->u2e[cn]) {
      dumpdata->w2_e1[dumpdata->nu1++] = -e;
      dumpdata->u2e[cn] = e;
      dumpdata->c2e[2][c] = e;
      dumpdata->e2c[e][0] = c;
      dumpdata->e2c[e][1] = cn;
      e++;
    }
  }
  for (i = 0; i < geom->nobc; i++) {
    open_bdrys_t *open = geom->open[i];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	cn = open->obc_t[cc];
	if (c == geom->m2d[c]) continue;
	if(!dumpdata->u1e[c]) {
	  dumpdata->w2_e1[dumpdata->nu1++] = e;
	  dumpdata->u1e[c] = e;
	  dumpdata->c2e[3][cn] = e;
	  dumpdata->e2c[e][0] = cn;
	  dumpdata->e2c[e][1] = c;
	  e++;
	}
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	cn = open->obc_t[cc];
	if (c == geom->m2d[c]) continue;
	if(!dumpdata->u2e[c]) {
	  dumpdata->w2_e1[dumpdata->nu1++] = -e;
	  dumpdata->u2e[c] = e;
	  dumpdata->c2e[2][cn] = e;
	  dumpdata->e2c[e][0] = cn;
	  dumpdata->e2c[e][1] = c;
	  e++;
	}
      }
    }
  }
  dumpdata->nedge3 = e-1;
  dumpdata->nu1--;

  /*-----------------------------------------------------------------*/
  /* Get the vertex mappings                                         */
  /* Count the vertices                                              */
  v = vv = 1;
  dumpdata->v2c = i_alloc_1d(geom->sgsizS);
  memset(dumpdata->v2c, 0, geom->sgsizS*sizeof(int));
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    /* Grid vertex is bottom left corner                             */
    dumpdata->u1v[c] = v;
    dumpdata->v2c[v] = c;
    v++;
    /* Check for ghost vertices                                      */
    cn = geom->xp1[c];
    if (geom->wgst[cn]) {
      dumpdata->u1v[cn] = v;
      dumpdata->v2c[v] = cn;
      v++;
    }
    cn = geom->yp1[c];
    if (geom->wgst[cn]) {
      dumpdata->u1v[cn] = v;
      dumpdata->v2c[v] = cn;
      v++;
    }
    cn = geom->yp1[geom->xp1[c]];
    if (geom->wgst[cn]) {
      dumpdata->u1v[cn] = v;
      dumpdata->v2c[v] = cn;
      v++;
    }
  }
  /* R_EDGE and F_EDGE open boundaries                               */
  for (i = 0; i < geom->nobc; i++) {
    open_bdrys_t *open = geom->open[i];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no2_e1; cc++) {
        c = open->obc_e1[cc];
	dumpdata->u1v[c] = v;
	dumpdata->v2c[v] = c;
	v++;
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no2_e2; cc++) {
        c = open->obc_e2[cc];
	dumpdata->u1v[c] = v;
	dumpdata->v2c[v] = c;
	v++;
      }
    }
  }
  dumpdata->nvertex2 = v-1;

  dumpdata->e2v = i_alloc_2d(2, dumpdata->nedge2+1);
  dumpdata->c2v = i_alloc_2d(geom->sgsizS, npem);

  /* Get the centre to vertex map                                    */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = cn = geom->w2_t[cc];
    v = dumpdata->u1v[c];
    dumpdata->c2v[1][c] = v;

    cn = geom->yp1[c];
    v = dumpdata->u1v[cn];
    dumpdata->c2v[2][c] = v;

    cn = geom->yp1[geom->xp1[c]];
    v = dumpdata->u1v[cn];
    dumpdata->c2v[3][c] = v;

    cn = geom->xp1[c];
    v = dumpdata->u1v[cn];
    dumpdata->c2v[4][c] = v;
  }

  for (e = 1; e <= dumpdata->nu1S; e++) {
    c = dumpdata->e2c[e][1];
    v = dumpdata->u1v[c];
    dumpdata->e2v[e][0] = v;
    cn = geom->yp1[c];
    v = dumpdata->u1v[cn];
    dumpdata->e2v[e][1] = v;
  }
  /* Get the edge to vertex map on u2 edges                          */
  for (e = 1; e <= dumpdata->nu1S; e++) {
    c = dumpdata->e2c[e][1];
    v = dumpdata->u1v[c];
    dumpdata->e2v[e][0] = v;
    cn = geom->xp1[c];
    v = dumpdata->u1v[cn];
    dumpdata->e2v[e][1] = v;
  }

  /*-----------------------------------------------------------------*/
  /* Extra sparse array used to pack arrays in df_sparse             */
  dumpdata->szm = szm = max(dumpdata->nface3+1, dumpdata->nedge3+1);
  dumpdata->w1s = d_alloc_1d(szm);
  dumpdata->w2s = d_alloc_1d(szm);
  if (dumpdata->face_dim == 1)
    dumpdata->i2s = i_alloc_2d(szm, npem + dumpdata->start_index);
  else
    dumpdata->i2s = i_alloc_2d(npem + dumpdata->start_index, szm);
  dumpdata->wc1 = d_alloc_2d(dumpdata->nface2+1, geom->nz);
  dumpdata->wc2 = d_alloc_2d(dumpdata->nface2+1, geom->nz);
  dumpdata->we1 = d_alloc_2d(dumpdata->nedge2+1, geom->nz);
  dumpdata->we2 = d_alloc_2d(dumpdata->nedge2+1, geom->nz);
  dumpdata->wv = d_alloc_1d(dumpdata->nvertex2+1);
}

/* END create_ugrid_maps()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the UGRID netCDF attributes                                */
/*-------------------------------------------------------------------*/
void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
				 nc_type fptype, char *modulo)
{
  int vid;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double fill;
  int ifill;
  int nface2 = dumpdata->nface2;
  int nedge2 = dumpdata->nedge2;
  int nvertex2 = dumpdata->nvertex2;
  int nface3 = dumpdata->nface3;
  int nedge3 = dumpdata->nedge3;
  int nvertex3 = dumpdata->nvertex3;
  int two = 2;
  char start_index[MAXSTRLEN];

  sprintf(start_index, "%d", dumpdata->start_index);

  vid = ncw_var_id(cdfid, "Mesh2");
  write_text_att(cdfid, vid, "cf_role", "mesh_topology");
  write_text_att(cdfid, vid, "long_name", "Topology data of 2D unstructured mesh");
  nc_put_att_int(cdfid, vid, "topology_dimension", NC_INT, 1, &two);
  write_text_att(cdfid, vid, "node_coordinates", "Mesh2_node_x Mesh2_node_y") ;
  write_text_att(cdfid, vid, "face_node_connectivity", "Mesh2_face_nodes");
  write_text_att(cdfid, vid, "face_dimension", "nMesh2_face");
  write_text_att(cdfid, vid, "edge_node_connectivity", "Mesh2_edge_nodes");
  write_text_att(cdfid, vid, "edge_dimension", "nMesh2_edge");
  write_text_att(cdfid, vid, "edge_coordinates", "Mesh2_edge_x Mesh2_edge_y");
  write_text_att(cdfid, vid, "face_coordinates", "Mesh2_face_x Mesh2_face_y");
  write_text_att(cdfid, vid, "face_edge_connectivity", "Mesh2_face_edges");
  write_text_att(cdfid, vid, "face_face_connectivity", "Mesh2_face_links");

  /* time independent variables */
  vid = ncw_var_id(cdfid, "Mesh2_layerfaces");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer faces");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "Mesh2_layers");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer centre");
  write_text_att(cdfid, vid, "coordinate_type", "Z");
  write_text_att(cdfid, vid, "positive", "up");

  vid = ncw_var_id(cdfid, "Mesh2_node_x");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "longitude");
    write_text_att(cdfid, vid, "long_name", "Longitude of 2D mesh nodes.");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate of 2D mesh nodes.");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_node_y");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "latitude");
    write_text_att(cdfid, vid, "long_name", "Latitude of 2D mesh nodes.");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate of 2D mesh nodes.");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_face_x");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "longitude");
    write_text_att(cdfid, vid, "long_name", "Characteristics longitude of 2D mesh face.");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "bounds", "Mesh2_face_xbnds");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate of 2D mesh face.");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "bounds", "Mesh2_face_xbnds");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_face_y");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "latitude");
    write_text_att(cdfid, vid, "long_name", "Longitude bounds of 2D mesh face (i.e. corner coordinates).");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "bounds", "Mesh2_face_ybnds");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate of 2D mesh face.");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "bounds", "Mesh2_face_ybnds");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_face_xbnds");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "longitude");
    write_text_att(cdfid, vid, "long_name", "Longitude bounds of 2D mesh face (i.e. corner coordinates).");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
    fill = Filld;
    nc_put_att_double(cdfid, vid, "_FillValue", NC_DOUBLE, 1, &fill);
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate of bounds of 2D mesh face (i.e. corner coordinates).");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    fill = Filld;
    nc_put_att_double(cdfid, vid, "_FillValue", NC_DOUBLE, 1, &fill);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_face_ybnds");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "latitude");
    write_text_att(cdfid, vid, "long_name", "Latitude bounds of 2D mesh face (i.e. corner coordinates).");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
    fill = Filld;
    nc_put_att_double(cdfid, vid, "_FillValue", NC_DOUBLE, 1, &fill);
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate of bounds of 2D mesh face (i.e. corner coordinates).");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    fill = Filld;
    nc_put_att_double(cdfid, vid, "_FillValue", NC_DOUBLE, 1, &fill);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_edge_x");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "longitude");
    write_text_att(cdfid, vid, "long_name", "Characteristic longitude of 2D mesh edge (e.g. midpoint of the edge).");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate of 2D mesh edge (e.g. midpoint of the edge).");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_edge_y");
  if (is_geog) {
    write_text_att(cdfid, vid, "standard_name", "latitude");
    write_text_att(cdfid, vid, "long_name", "Characteristic latitude of 2D mesh edge (e.g. midpoint of the edge).");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate of 2D mesh edge (e.g. midpoint of the edge).");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  vid = ncw_var_id(cdfid, "Mesh2_face_nodes");
  write_text_att(cdfid, vid, "cf_role", "face_node_connectivity");
  write_text_att(cdfid, vid, "long_name", "Maps every face to its corner nodes.");
  ifill = Filli;
  nc_put_att_int(cdfid, vid, "_FillValue", NC_INT, 1, &ifill);
  write_text_att(cdfid, vid, "start_index", start_index);

  vid = ncw_var_id(cdfid, "Mesh2_edge_nodes");
  write_text_att(cdfid, vid, "cf_role", "edge_node_connectivity");
  write_text_att(cdfid, vid, "long_name", "Maps every edge to the two nodes that it connects.");
  write_text_att(cdfid, vid, "start_index", start_index);

  vid = ncw_var_id(cdfid, "Mesh2_face_edges");
  write_text_att(cdfid, vid, "cf_role", "face_edge_connectivity");
  write_text_att(cdfid, vid, "long_name", "Maps every face to its edges.");
  ifill = Filli;
  nc_put_att_int(cdfid, vid, "_FillValue", NC_INT, 1, &ifill);
  write_text_att(cdfid, vid, "start_index", start_index);

  vid = ncw_var_id(cdfid, "Mesh2_face_links");
  write_text_att(cdfid, vid, "cf_role", "face_face_connectivity");
  write_text_att(cdfid, vid, "long_name", "Indicates which faces are neighbors.");
  write_text_att(cdfid, vid, "start_index", start_index);

  vid = ncw_var_id(cdfid, "Mesh2_depth");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at sea-bed at cell centre");

  vid = ncw_var_id(cdfid, "h1au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Edge length between nodes");
  write_text_att(cdfid, vid, "coordinates", "Mesh2_edge_x, Mesh2_edge_y");

  vid = ncw_var_id(cdfid, "h2au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Length between face centres");
  write_text_att(cdfid, vid, "coordinates", "Mesh2_edge_x, Mesh2_edge_y");

  vid = ncw_var_id(cdfid, "thetau1");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation of normal vector");
  write_text_att(cdfid, vid, "coordinates", "Mesh2_edge_x, Mesh2_edge_y");

  vid = ncw_var_id(cdfid, "coriolis");
  write_text_att(cdfid, vid, "units", "s-1");
  write_text_att(cdfid, vid, "long_name", "Coriolis parameter");
  write_text_att(cdfid, vid, "coordinates", "Mesh2_face_x, Mesh2_face_y");

  /* time dependent variables */

  vid = ncw_var_id(cdfid, "t");
  write_text_att(cdfid, vid, "units", dumpdata->output_tunit);
  write_text_att(cdfid, vid, "long_name", "Time");
  write_text_att(cdfid, vid, "coordinate_type", "time");
  if (strlen(modulo))
    write_text_att(cdfid, vid, "modulo", modulo);

  if (geom->sednz > 0) {
    vid = ncw_var_id(cdfid, "Mesh2_layerfaces_sed");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate at sediment layer faces");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
    
    vid = ncw_var_id(cdfid, "Mesh2_layers_sed");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate at sediment layer centres");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  if ((vid = ncw_var_id(cdfid, "u1av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of depth averaged current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "u2av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "tangential component of depth averaged current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "uav")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "east component of depth averaged current at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "vav")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "north component of depth averaged current at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "wtop")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "Vertical velocity at surface");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "eta")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Surface Elevation");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "topz")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate for surface cell");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "wind1")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of wind stress at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "windx")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "east component of wind stress at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "windy")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "north component of wind stress at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "patm")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "u1")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "u2")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "tangential component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "umean")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "mean normal component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y, Mesh2_layers");
  }
  if ((vid = ncw_var_id(cdfid, "u1vmean")) >= 0) {
    write_text_att(cdfid, vid, "units", "m3s-1");
    write_text_att(cdfid, vid, "long_name",
                   "mean normal component of volume flux at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_edge_x, Mesh2_edge_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "u")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "east component of current at face centre");
    write_text_att(cdfid, vid, "standard_name","eastward_sea_water_velocity");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "v")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "north component of current at face centre");
    write_text_att(cdfid, vid, "standard_name","northward_sea_water_velocity");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "w")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "K component of current at cell centre and Z grid");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "dens")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Density");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "dens_0")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Potential density");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
    fill = 1025.0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Kz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Kz");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Vz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Vz");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "ptconc")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Particle concentration");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "Cd")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Bottom drag coefficient");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "u1kh")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre2 second-1");
    write_text_att(cdfid, vid, "long_name",
                   "horizontal diffusivity at centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "ustrcw")) >= 0) {
    write_text_att(cdfid, vid, "units", "ms-1");
    write_text_att(cdfid, vid, "long_name", "Bottom friction velocity");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "dz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m");
    write_text_att(cdfid, vid, "long_name", "Cell centered thickness");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }

  if ((vid = ncw_var_id(cdfid, "dzu1")) >= 0) {
    write_text_att(cdfid, vid, "units", "m");
    write_text_att(cdfid, vid, "long_name", "Face centered thickness");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
  }


  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "t, Mesh2_face_x, Mesh2_face_y, Mesh2_layers"}
    };
    tracer_write_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "t, Mesh2_face_x, Mesh2_face_y"}
    };
    tracer_write_nc(cdfid, dumpdata->ntrS, dumpdata->trinfo_2d, 2, attr);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", codeheader);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "UGRID-1.0");
  nc_put_att_int(cdfid, NC_GLOBAL, "nface2", NC_INT, 1, &nface2);
  nc_put_att_int(cdfid, NC_GLOBAL, "nedge2", NC_INT, 1, &nedge2);
  nc_put_att_int(cdfid, NC_GLOBAL, "nvertex2", NC_INT, 1, &nvertex2);
  nc_put_att_int(cdfid, NC_GLOBAL, "nface3", NC_INT, 1, &nface3);
  nc_put_att_int(cdfid, NC_GLOBAL, "nedge3", NC_INT, 1, &nedge3);
  nc_put_att_int(cdfid, NC_GLOBAL, "nvertex3", NC_INT, 1, &nvertex3);
  if (dumpdata->nce1 > 0 && dumpdata->nce2 > 0) {
    nc_put_att_int(cdfid, NC_GLOBAL, "NCE1", NC_INT, 1, &dumpdata->nce1);
    nc_put_att_int(cdfid, NC_GLOBAL, "NCE2", NC_INT, 1, &dumpdata->nce2);
  }
  nc_put_att_int(cdfid, NC_GLOBAL, "start_index", NC_INT, 1, &dumpdata->start_index);
  nc_put_att_int(cdfid, NC_GLOBAL, "face_dim", NC_INT, 1, &dumpdata->face_dim);
  nc_put_att_int(cdfid, NC_GLOBAL, "edge_dim", NC_INT, 1, &dumpdata->edge_dim);
  if (dumpdata->runno >= 0)
    write_text_att(cdfid, NC_GLOBAL, "Run_ID", dumpdata->runnoc);
  /*nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);*/
  if (strlen(dumpdata->runcode))
    write_text_att(cdfid, NC_GLOBAL, "Run_code", dumpdata->runcode);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);
  if (strlen(dumpdata->trl))
    write_text_att(cdfid, NC_GLOBAL, "Technology_Readiness_Level", dumpdata->trl);
  if (strlen(dumpdata->reference))
    write_text_att(cdfid, NC_GLOBAL, "Output_Reference", dumpdata->reference);
  write_mean_atts(dumpdata, cdfid);
}

/* END write_dump_attributes_ugrid()                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* UGRID reset function                                              */
/*-------------------------------------------------------------------*/
void df_ugrid_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_ugrid_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t < (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "t", df->tout);
}

/* END df_ugrid_reset()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates the UGRID file                                            */
/*-------------------------------------------------------------------*/
void *df_ugrid_create(dump_data_t *dumpdata, dump_file_t *df)
{
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int n;
  int cdfid;
  int dims[10];
  int vid;                      /* df_variable_t dimension */
  /* dimension ids */
  int recdimid;                 /* record dimension id */
  int nMesh2_nodeid;            /* Node dimension */
  int nMesh2_edgeid;            /* Edge dimension */
  int nMesh2_faceid;            /* Face dimension */
  int nMesh2_face_linksid;      /* Number of face pairs */
  int nMaxMesh2_face_nodesid;   /* Max # nodes per face */
  int Mesh2_layersid;           /* Number of layers */ 
  int Mesh2_layerfacesid;       /* Number of layer faces */ 
  int Twoid;                    /* Id for 2 */
  int nface2;
  int nedge2;
  int nvertex2;
  int npem;
  int kcentreid_sed;            /* K dimension id at grid centre for
                                   sediments */
  int kgridid_sed;              /* K dimension id at grid corner for
                                   sediments */
  FILE *fp;
  int nc_mode;

  df_ugrid_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);

  nface2 = dumpdata->nface2 - dumpdata->start_index;
  nedge2 = dumpdata->nedge2 - dumpdata->start_index;
  nvertex2 = dumpdata->nvertex2 - dumpdata->start_index;
  npem = dumpdata->npe - dumpdata->start_index;

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_ugrid_init_data(dumpdata, df, cdfid);
      data = (df_ugrid_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "t", dumpdata->t);
      return data;
    }
  }

  // Get nc_mode
  nc_mode = get_nc_mode(df);

  /* create the netCDF file */
  if (ncw_create(df->name, nc_mode, &cdfid) != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "record", NC_UNLIMITED, &recdimid);

  nc_def_dim(cdfid, "Mesh2_layers", df->nz, &Mesh2_layersid);
  nc_def_dim(cdfid, "Mesh2_layerfaces", df->nz + 1, &Mesh2_layerfacesid);
  if (df->nz_sed > 0) {
    nc_def_dim(cdfid, "Mesh2_layerfaces_sed", df->nz_sed + 1, &kgridid_sed);
    nc_def_dim(cdfid, "Mesh2_layers_sed", df->nz_sed, &kcentreid_sed);
  }
  nc_def_dim(cdfid, "nMesh2_node", nvertex2, &nMesh2_nodeid);
  nc_def_dim(cdfid, "nMesh2_edge", nedge2, &nMesh2_edgeid);
  nc_def_dim(cdfid, "nMesh2_face", nface2, &nMesh2_faceid);

  nc_def_dim(cdfid, "nMaxMesh2_face_nodes", npem, &nMaxMesh2_face_nodesid);
  nc_def_dim(cdfid, "nMesh2_face_links", nedge2, &nMesh2_face_linksid);
  nc_def_dim(cdfid, "Two", 2, &Twoid);

  /* Grid */
  nc_def_var(cdfid, "Mesh2", NC_INT, 0, NULL, &vid);

  /* Connectivity */
  set_dims(dumpdata->face_dim, dims, nMaxMesh2_face_nodesid, nMesh2_faceid);
  nc_def_var(cdfid, "Mesh2_face_nodes", NC_INT, 2, dims, &vid);
  set_dims(dumpdata->edge_dim, dims, Twoid, nMesh2_edgeid);
  nc_def_var(cdfid, "Mesh2_edge_nodes", NC_INT, 2, dims, &vid);
  set_dims(dumpdata->face_dim, dims, nMaxMesh2_face_nodesid, nMesh2_faceid);
  nc_def_var(cdfid, "Mesh2_face_edges", NC_INT, 2, dims, &vid);
  set_dims(dumpdata->edge_dim, dims, Twoid, nMesh2_face_linksid);
  nc_def_var(cdfid, "Mesh2_face_links", NC_INT, 2, dims, &vid);

  /* Mesh node coordinates */
  dims[0] = nMesh2_nodeid;
  nc_def_var(cdfid, "Mesh2_node_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh2_node_y", NC_DOUBLE, 1, dims, &vid);

  /* Mesh face and edge coordinate variables */
  dims[0] = nMesh2_faceid;
  nc_def_var(cdfid, "Mesh2_face_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh2_face_y", NC_DOUBLE, 1, dims, &vid);
  set_dims(dumpdata->face_dim, dims, nMaxMesh2_face_nodesid, nMesh2_faceid);
  nc_def_var(cdfid, "Mesh2_face_xbnds", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "Mesh2_face_ybnds", NC_DOUBLE, 2, dims, &vid);
  dims[0] = nMesh2_edgeid;
  nc_def_var(cdfid, "Mesh2_edge_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh2_edge_y", NC_DOUBLE, 1, dims, &vid);

  /* Vertical structure */
  dims[0] = Mesh2_layersid;
  nc_def_var(cdfid, "Mesh2_layers", NC_DOUBLE, 1, dims, &vid);
  dims[0] = Mesh2_layerfacesid;
  nc_def_var(cdfid, "Mesh2_layerfaces", NC_DOUBLE, 1, dims, &vid);
  dims[0] = nMesh2_faceid;
  nc_def_var(cdfid, "Mesh2_depth", NC_DOUBLE, 1, dims, &vid);

  /* Cell lengths */
  dims[0] = nMesh2_edgeid;
  nc_def_var(cdfid, "h2au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h1au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "thetau1", NC_DOUBLE, 1, dims, &vid);
  dims[0] = nMesh2_faceid;
  nc_def_var(cdfid, "coriolis", fptype, 1, dims, &vid);
  
  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid, "t", NC_DOUBLE, 1, dims, &vid);

  if (dumpdata->sednz > 0) {
    dims[0] = recdimid;
    dims[1] = kgridid_sed;
    dims[2] = nMesh2_faceid;
    nc_def_var(cdfid, "Mesh2_layerfaces_sed", NC_DOUBLE, 3, dims, &vid);

    dims[1] = kcentreid_sed;
    nc_def_var(cdfid, "Mesh2_layers_sed", NC_DOUBLE, 3, dims, &vid);
  }

  df_ugrid_init_data(dumpdata, df, cdfid);
  data = (df_ugrid_data_t *)df->private_data;

  /*
   * Set chunksizes to be a single horizontal layer
   */
  for (n = 0; n < df->nvars; n++) {
    size_t chunksize[] = {1,1,1};
    if (data->vars[n].ndims == 1) {
      df_ugrid_get_dimids(cdfid, data->vars[n], &dims[1], NULL);
      ncw_def_var2(df->name,cdfid, df->vars[n], data->vars[n].type, 2, dims, &vid, df->compress);
#ifdef NC_NETCDF4
      if (nc_mode & NC_NETCDF4) {
	chunksize[0] = 1;
	(void)nc_inq_dimlen(cdfid, dims[1], &chunksize[1]);
	ncw_def_var_chunking(df->name, cdfid, vid, chunksize);
      }
#endif
    } else if (data->vars[n].ndims == 2) {
      df_ugrid_get_dimids(cdfid, data->vars[n], &dims[2], &dims[1]);
      ncw_def_var2(df->name,cdfid, df->vars[n], data->vars[n].type, 3, dims, &vid, df->compress);
#ifdef NC_NETCDF4
      if (nc_mode & NC_NETCDF4) {
	chunksize[0] = 1;
	chunksize[1] = 1;
	(void)nc_inq_dimlen(cdfid, dims[2], &chunksize[2]);
	ncw_def_var_chunking(df->name, cdfid, vid, chunksize);
      }
#endif
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }
  }

  write_dump_attributes_ugrid(dumpdata, cdfid, fptype, df->modulo);

  nc_enddef(cdfid);

  df_ugrid_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}

/* END df_ugrid_create()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes data to UGRID                                              */
/*-------------------------------------------------------------------*/
void df_ugrid_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double st, newt = t;
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;
  int fid = data->fid;
  geometry_t *geom = dumpdata->master->geom;
  int c,cc;
  int oset = (dumpdata->start_index) ? 0 : 1;

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, dumpdata->output_tunit, 
		       &newt, 1);
  if (strlen(df->modulo)) {
    /*st = schedule->start_time;*/
    st = df->tstart;
    tm_change_time_units(dumpdata->timeunit, dumpdata->output_tunit, 
			 &st, 1);
    newt -= st;
  }
  nc_put_vara_double(fid, ncw_var_id(fid, "t"), start, count, &newt);

  if (dumpdata->sednz > 0) {
    count[1] = dumpdata->sednz;
    count[2] = dumpdata->nface2;
    pack_ugrids(geom->w2_t, dumpdata->sednz - 1, count[2], geom->cellz_sed, dumpdata->wc1, oset);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_layers_sed"), start, count, dumpdata->wc1);

    count[1] = dumpdata->sednz + 1;
    pack_ugrids(geom->w2_t, dumpdata->sednz, count[2], geom->gridz_sed, dumpdata->wc1, oset);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_layerfaces_sed"), start, count, dumpdata->wc1);
  }

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_ugrid_var_t vn = data->vars[n];
    void *p = vn.v;
    start[1] = df->ilower;
    if (vn.ndims == 1) {
      double *w1 = dumpdata->w1s;
      double *w2 = NULL;
      memset(w1,  0, dumpdata->nface2*sizeof(double));
      df_ugrid_get_dimsizes(df, vn, NULL, &count[1], NULL);
      if (vn.vector_mode != VM_NONE) {
	void *p1 = vn.v1;
	void *p2 = vn.v2;
	w2 = dumpdata->w2s;
	if (vn.func2 != NULL) {
	  vn.func2(dumpdata, (*(double **)p1), (*(double **)p2), w1, w2, oset);
	}
      } else
	pack_ugrid1(vn.hmap, count[1], (*(double **)p), w1, oset);
      if (vn.vector_mode & VM_EAST)
	nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count, w1);
      else if (vn.vector_mode & VM_NORTH)
	nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count, w2);
      else
	nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count, w1);

    }  else if (vn.ndims == 2) {
      size_t sz;
      int cc,c,k, kb;
      double **w1 = dumpdata->wc1;
      double **w2 = NULL;

      /* Initialize land cells */
      for (k = dumpdata->nz-1; k >= 0; k--)
	for (cc = 0; cc < count[2]; cc++)
	  w1[k][cc] = NaN;

      if (!vn.sediment) {
	if (vn.vector_mode != VM_NONE) {
	  void *p1 = vn.v1;
	  void *p2 = vn.v2;
	  w1 = dumpdata->we1;
	  w2 = dumpdata->we2;
	  if (vn.func3 != NULL) {
	    vn.func3(dumpdata, (*(double **)p1), (*(double **)p2), w1, w2, oset);
	  }
	} else 
	  pack_ugrid2(geom, vn.hmap, sz, (*(double **)p), w1, oset);
	start[1] = df->klower;
	df_ugrid_get_dimsizes(df, vn, &sz, &count[2], &count[1]);
	if (vn.vector_mode & VM_EAST)
	  nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, w1);
	else if (vn.vector_mode & VM_NORTH)
	  nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, w2);
	else
	  nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, w1);
      } else {
	start[1] = 0;
	df_ugrid_get_dimsizes(df, vn, NULL, &count[2], &count[1]);
	/* Initialize land cells */
	for (k = dumpdata->sednz-1; k >= 0; k--)
	  for (cc = 0; cc < count[2]; cc++)
	    w1[k][cc] = NaN;
	pack_ugrids(geom->w2_t, dumpdata->sednz, count[2], (*(double ***)p), w1, oset);
	nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, w1);
      }
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

/* END df_ugrid_write()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Closes the UGRID file                                             */
/*-------------------------------------------------------------------*/
void df_ugrid_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}

/* END df_ugrid_close()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes model geometry to UGRID                                    */
/*-------------------------------------------------------------------*/
static void df_ugrid_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int fid = data->fid;
  int c, cc, e, ee, vv, v;
  double **xbnds, **ybnds;
  int oset = (dumpdata->start_index) ? 0 : 1;
  void (*i2) (int *map, int size, int np, int **var, int **pack, int oset) = NULL;
  void (*i2a) (int *map, int size, int np, int **var, int **pack, int oset) = NULL;
  int n, i, j;
  double *we;
  int npem = dumpdata->npe;

  if (dumpdata->face_dim)
    i2 = pack_ugrid_i2;
  else
    i2 = pack_ugrid_ri2;
  if (dumpdata->edge_dim)
    i2a = pack_ugrid_i2a;
  else
    i2a = pack_ugrid_ri2a;

  set_longitude(dumpdata, df, 1);
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /* time independent variables */
  count[0] = dumpdata->nz + 1;
  nc_put_vara_double(fid, ncw_var_id(fid, "Mesh2_layerfaces"), start, count,
                     dumpdata->gridz);
  count[0] = dumpdata->nz;
  nc_put_vara_double(fid, ncw_var_id(fid, "Mesh2_layers"), start, count,
                     dumpdata->cellz);

  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  i2(geom->w2_t, dumpdata->nface2, npem,  dumpdata->c2v, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_nodes"), start, count, dumpdata->i2s);

  set_count(dumpdata->edge_dim, count, 2, dumpdata->nedge2);
  i2a(NULL, dumpdata->nedge2, 2, dumpdata->e2v, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_edge_nodes"), start, count, dumpdata->i2s);

  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  i2(geom->w2_t, dumpdata->nface2, npem, dumpdata->c2e, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_edges"), start, count, dumpdata->i2s);

  set_count(dumpdata->edge_dim, count, 2, dumpdata->nedge2+dumpdata->start_index);
  i2a(geom->c2cc, dumpdata->nedge2, 2, dumpdata->e2c, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_links"), start, count, dumpdata->i2s);

  count[0] = dumpdata->nvertex2;
  for (v = 1; v <= dumpdata->nvertex2; v++) {
    c = dumpdata->v2c[v];
    i = dumpdata->c2i[c];
    j = dumpdata->c2j[c];
    dumpdata->wv[v] = dumpdata->gridx[j][i];
  }
  pack_ugrid1(NULL, dumpdata->nvertex2, dumpdata->wv, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_node_x"), start, count, dumpdata->w1s);
  for (v = 1; v <= dumpdata->nvertex2; v++) {
    c = dumpdata->v2c[v];
    i = dumpdata->c2i[c];
    j = dumpdata->c2j[c];
    dumpdata->wv[v] = dumpdata->gridy[j][i];
  }
  pack_ugrid1(NULL, dumpdata->nvertex2, dumpdata->wv, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_node_y"), start, count, dumpdata->w1s);

  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->cellx, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_face_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->celly, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_face_y"), start, count, dumpdata->w1s);

  we = d_alloc_1d(dumpdata->nedge2+1);
  count[0] = dumpdata->nedge2;
  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      i = dumpdata->c2i[c];
      j = dumpdata->c2j[c];
      if (i >= 0 && j >= 0) we[e] = dumpdata->u1x[j][i];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      i = dumpdata->c2i[c];
      j = dumpdata->c2j[c];
      if (i >= 0 && j >= 0) we[e] = dumpdata->u2x[j][i];
    }
  }
  pack_ugrid1(NULL, dumpdata->nedge2, we, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_edge_x"), start, count, dumpdata->w1s);
  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      i = dumpdata->c2i[c];
      j = dumpdata->c2j[c];
      if (i >= 0 && j >= 0) we[e] = dumpdata->u1y[j][i];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      i = dumpdata->c2i[c];
      j = dumpdata->c2j[c];
      if (i >= 0 && j >= 0) we[e] = dumpdata->u2y[j][i];
    }
  }
  pack_ugrid1(NULL, dumpdata->nedge2, we, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_edge_y"), start, count, dumpdata->w1s);

  if (dumpdata->face_dim) {
    xbnds = d_alloc_2d(geom->sgsizS, npem+1);
    ybnds = d_alloc_2d(geom->sgsizS, npem+1);
  } else {
    xbnds = d_alloc_2d(npem+1, geom->sgsizS);
    ybnds = d_alloc_2d(npem+1, geom->sgsizS);
  }

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (vv = 1; vv <= npem; vv++) {
      v = dumpdata->c2v[vv][c];
      if (dumpdata->face_dim) {
	xbnds[vv-oset][cc-oset] = (vv <= npem) ? geom->gridx[v] : Filld;
	ybnds[vv-oset][cc-oset] = (vv <= npem) ? geom->gridy[v] : Filld;
      } else {
	xbnds[cc-oset][vv-oset] = (vv <= npem) ? geom->gridx[v] : Filld;
	ybnds[cc-oset][vv-oset] = (vv <= npem) ? geom->gridy[v] : Filld;
      }
    }
  }

  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_xbnds"), start, count, xbnds);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_ybnds"), start, count, ybnds);
  d_free_2d(xbnds);
  d_free_2d(ybnds);

  count[0] = dumpdata->nedge2;
  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      we[e] = geom->h2au1[c];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      we[e] = geom->h1au2[c];
    }
  }
  pack_ugrid1(NULL, dumpdata->nedge2, we, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au1"), start, count, dumpdata->w1s);

  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      we[e] = geom->h1au1[c];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      we[e] = geom->h2au2[c];
    }
  }
  pack_ugrid1(NULL, dumpdata->nedge2, we, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au1"), start, count, dumpdata->w1s);

  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      we[e] = geom->thetau1[c];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      we[e] = geom->thetau2[c] + PI / 2.0;
    }
  }
  pack_ugrid1(NULL, dumpdata->nedge2, we, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau1"), start, count, dumpdata->w1s);
  d_free_1d(we);

  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, master->coriolis, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "coriolis"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->botz, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_depth"), start, count, dumpdata->w1s);
  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}

/* END df_ugrid_writegeom()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Supplies variable information for UGRID                           */
/*-------------------------------------------------------------------*/
int df_ugrid_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
			 char *name, df_ugrid_var_t *var)
{
  int found = 1;
  int n = 0;
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;

  var->v = NULL;
  var->ndims = 1; // This is the only case for sparse
  var->type = nc_get_default_type(df->bpv);
  var->xylocation = CL_NONE;
  var->zlocation = CL_NONE;
  var->vector_mode = VM_NONE;
  var->sediment = 0;
  var->hmap = geom->w2_t;
  var->vmap = geom->s2k;
  var->m2d = geom->m2d;
  var->func2 = NULL;
  var->func3 = NULL;

  if (strcmp(name, "u1av") == 0) {
    var->v1 = (void *)&master->u1av;
    var->v2 = (void *)&master->u2av;
    var->xylocation = CL_SP2|CL_EDGE;
    var->vector_mode = VM_EAST;
    var->func2 = vel_tan2D;
  }

  else if (strcmp(name, "u2av") == 0) {
    var->v1 = (void *)&master->u1av;
    var->v2 = (void *)&master->u2av;
    var->xylocation = CL_SP2|CL_EDGE;
    var->vector_mode = VM_NORTH;
    var->func2 = vel_tan2D;
  }

  else if (strcmp(name, "uav") == 0) {
    var->v1 = (void *)&master->u1av;
    var->v2 = (void *)&master->u2av;
    var->xylocation = CL_SP2|CL_FACE;
    var->vector_mode = VM_EAST;
    var->func2 = vel_cen2D;
  }

  else if (strcmp(name, "vav") == 0) {
    var->v1 = (void *)&master->u1av;
    var->v2 = (void *)&master->u2av;
    var->xylocation = CL_SP2|CL_FACE;
    var->vector_mode = VM_NORTH;
    var->func2 = vel_cen2D;
  }

  else if (strcmp(name, "wtop") == 0) {
    var->v = (void *)&master->wtop;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "topz") == 0) {
    var->v = (void *)&master->topz;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "eta") == 0) {
    var->v = (void *)&master->eta;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "wind1") == 0) {
    var->v1 = (void *)&master->wind1;
    var->v2 = (void *)&master->wind2;
    var->xylocation = CL_SP2|CL_EDGE;
    var->vector_mode = VM_EAST;
    var->func2 = vel_tan2D;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void *)&master->patm;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "u") == 0) {
    var->v1 = (void *)&master->u1;
    var->v2 = (void *)&master->u2;
    var->xylocation = CL_SP3|CL_FACE;
    var->vector_mode = VM_EAST;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->func3 = vel_cen3D;
  }

  else if (strcmp(name, "v") == 0) {
    var->v1 = (void *)&master->u1;
    var->v2 = (void *)&master->u2;
    var->xylocation = CL_SP3|CL_FACE;
    var->vector_mode = VM_NORTH;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->func3 = vel_cen3D;
  }

  else if (strcmp(name, "u1") == 0) {
    var->v1 = (void *)&master->u1;
    var->v2 = (void *)&master->u2;
    var->xylocation = CL_SP3|CL_EDGE;
    var->vector_mode = VM_EAST;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->func3 = vel_tan3D;
  }

  else if (strcmp(name, "u2") == 0) {
    var->v1 = (void *)&master->u1;
    var->v2 = (void *)&master->u2;
    var->xylocation = CL_SP3|CL_EDGE;
    var->vector_mode = VM_NORTH;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->func3 = vel_tan3D;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void *)&master->w;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "u1mean") == 0) {
    var->v1 = (void *)&master->u1m;
    var->v2 = (void *)&master->u2m;
    var->xylocation = CL_SP3|CL_FACE;
    var->vector_mode = VM_EAST;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->func3 = vel_cen3D;
  }

  else if (strcmp(name, "u2mean") == 0) {
    var->v1 = (void *)&master->u1m;
    var->v2 = (void *)&master->u2m;
    var->xylocation = CL_SP3|CL_FACE;
    var->vector_mode = VM_NORTH;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->func3 = vel_cen3D;
  }

  else if (strcmp(name, "umean") == 0) {
    var->v1 = (void *)&master->u1m;
    var->v2 = (void *)&master->u2m;
    var->xylocation = CL_SP3|CL_EDGE;
    var->vector_mode = VM_MAG;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->func3 = vel_tan3D;
  }

  else if (strcmp(name, "u1vmean") == 0) {
    var->v1 = (void *)&master->u1vm;
    var->v2 = (void *)&master->u2vm;
    var->xylocation = CL_SP3|CL_EDGE;
    var->vector_mode = VM_MAG;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->func3 = vel_tan3D;
  }

  else if (strcmp(name, "dens") == 0) {
    var->v = (void *)&master->dens;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "dens_0") == 0) {
    var->v = (void *)&master->dens_0;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "Kz") == 0) {
    var->v = (void *)&master->Kz;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "Vz") == 0) {
    var->v = (void *)&master->Vz;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "dz") == 0) {
    var->v = (void *)&master->dz;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "dzu1") == 0) {
    var->v1 = (void *)&master->dzu1;
    var->v2 = (void *)&master->dzu2;
    var->vector_mode = VM_MAG;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->func3 = vel_scalar;
  }

  else if (strcmp(name, "Cd") == 0) {
    var->v = (void *)&master->Cd;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "u1kh") == 0) {
    var->v = (void *)&master->u1kh;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
    var->vmap = geom->s2k;
    var->m2d = geom->m2d;
  }

  else
    found = 0;

  if (!found) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) {
        var->v = (void **)&master->tr_wc[n];
        var->xylocation = CL_SP3|CL_FACE;
	var->ndims = 2;
	var->zlocation = CL_CENTRE;
	var->hmap = geom->w3_t;
        var->sediment = 0;
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, name) == 0) {
        var->v = (void **)&master->tr_wcS[n];
        var->xylocation = CL_SP2|CL_FACE;
        var->sediment = 0;
        found = 1;
        break;
      }
    }
  }

  if (!found) {
    /* Sediment tracers */
    for (n = 0; n < dumpdata->nsed; n++) {
      char name1[MAXSTRLEN];
      strcpy(name1, dumpdata->trinfo_sed[n].name);
      strcat(name1, "_sed");
      if (strcmp(name, name1) == 0) {
        var->ndims = 2;
        var->v = (void **)&master->tr_sed[n];
        var->xylocation = CL_SP2|CL_FACE;
        var->zlocation = CL_CENTRE;
        var->sediment = 1;
        found = 1;
        break;
      }
    }
  }

  return found;
}

/* END df_ugrid_get_varinfo()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* UGRID variable data initialisation                                */
/*-------------------------------------------------------------------*/
static void df_ugrid_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid)
{
  int i, n;
  df_ugrid_data_t *data = NULL;

  df_parse_vars(dumpdata,df,NULL,UGRID_ALL_VARS);

  // Clean up before allocating more memory
  if (df->private_data != NULL) {
    if (((df_ugrid_data_t *)df->private_data)->vars != NULL)
      free(((df_ugrid_data_t *)df->private_data)->vars);
    free(df->private_data);
  }
  
  data = (df_ugrid_data_t *)malloc(sizeof(df_ugrid_data_t));
  memset(data, 0, sizeof(df_ugrid_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_ugrid_var_t *)malloc(df->nvars * sizeof(df_ugrid_var_t))) == NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_ugrid_var_t *var = &data->vars[i];
    if (df_ugrid_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_ugrid_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

  if (df->diagn && DEBUG("dump"))
    dlog("dump",
         "Diagnostic tracers initialisation has been synchronised with \"%s\" output; period = %f s.\n",
         df->name, df->tinc);

}

/* END df_ugrid_init_data()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Supplies variable dimension ids                                   */
/*-------------------------------------------------------------------*/
static void df_ugrid_get_dimids(int cdfid, df_ugrid_var_t var, int *ns, int *zid)
{
  if (var.xylocation & CL_FACE) 
    *ns = ncw_dim_id(cdfid, "nMesh2_face");
  else if (var.xylocation & CL_EDGE) 
    *ns = ncw_dim_id(cdfid, "nMesh2_edge");
  else if (var.xylocation & CL_VERTEX) 
    *ns = ncw_dim_id(cdfid, "nMesh2_node");
  else
    hd_quit("df_ugrid_get_dimids: can't find dimensions.\n");

  if (!var.sediment) {
    switch (var.zlocation) {
    case CL_GRID:
      *zid = ncw_dim_id(cdfid, "Mesh2_layerfaces");
      break;

    case CL_CENTRE:
      *zid = ncw_dim_id(cdfid, "Mesh2_layers");
      break;

    default:
      break;
    }
  } else {
    switch (var.zlocation) {
    case CL_GRID:
      *zid = ncw_dim_id(cdfid, "Mesh2_layerfaces_sed");
      break;

    case CL_CENTRE:
    default:
      *zid = ncw_dim_id(cdfid, "Mesh2_layers_sed");
      break;
    }
  }

}

/* END df_ugrid_get_dimids()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Supplies variable dimension sizes                                 */
/*-------------------------------------------------------------------*/

static void df_ugrid_get_dimsizes(dump_file_t *df, 
				  df_ugrid_var_t var, 
				  size_t *sz,   /* 3D vector size    */
				  size_t * ns,  /* 2D vector size    */
				  size_t *nz)   /* Number of layers  */
{
  if (var.xylocation == (CL_SP2|CL_FACE))
    *ns = df->nface2;
  else if (var.xylocation == (CL_SP2|CL_EDGE))
    *ns = df->nedge2;
  else if (var.xylocation == (CL_SP2|CL_VERTEX))
    *ns = df->nvertex2;
  else if (var.xylocation == (CL_SP3|CL_FACE)) {
    *ns = df->nface2;
    *sz = df->nface3;
  } else if (var.xylocation == (CL_SP3|CL_EDGE)) {
    *ns = df->nedge2;
    *sz = df->nedge3;
  } else if (var.xylocation == (CL_SP3|CL_VERTEX)) {
    *ns = df->nvertex2;
    *sz = df->nvertex3;
  } else
    hd_quit("df_ugrid:df_sp_get_dimsizes error in xylocation\n");

  if (!var.sediment) {
    switch (var.zlocation) {
    case CL_GRID:
      *nz = df->nz + 1;
      break;

    case CL_CENTRE:
      *nz = df->nz;
      break;

    default:
      break;
    }
  }

  else {
    switch (var.zlocation) {
    case CL_GRID:
      *nz = df->nz_sed + 1;
      break;

    case CL_CENTRE:
      *nz = df->nz_sed;
      break;

    default:
      break;
    }
  }
}

/* END df_ugrid_get_dimsizes()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the sparse array into consecutive wet cells. If 'map' is    */
/* geom->w2_t then cell centres are mapped to a consecutive array.   */
/* Note; unlike COMPAS, the first b2_t centres are not consecutive,  */
/* since R_EDGE and F_EDGE OBC cells may be included, so the w2_t    */
/* mapping must be used to get wet cell centres.                     */
/* Edge and vertex cells are consecutive, and don't need a mapping.  */
/* Note the the first null cell in the sparse array is removed so    */
/* that the dumped variable indies go from 0:ns-1 whereas the sparse */
/* indicies go from 1:ns with value[0]=0.                            */
/*-------------------------------------------------------------------*/
void pack_ugrid1(int *map, int mapsize, double *var, double *pack, int oset)
{
  int cc, c;

  for (cc = 1; cc <= mapsize; cc++) {
    c = (map == NULL) ? cc : map[cc];
    pack[cc-oset] = var[c];
  }
}

/* END pack_ugrid1()                                                 */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Packs the sparse array into consecutive wet cells.                */
/* Note the the first null cell in the sparse array is removed so    */
/* that the dumped variable indies go from 0:ns-1 whereas the sparse */
/* indicies go from 1:ns with value[0]=0.                            */
/*-------------------------------------------------------------------*/
void pack_ugrid2(geometry_t *geom, int *hmap, int mapsize, double *var, double **pack, int oset)
{
  int cc, c, c2, k, cc2;

  for (cc = 1; cc <= mapsize; cc++) {
    c = hmap[cc];
    c2 = geom->m2d[c];
    cc2 = geom->c2cc[c2];
    k = geom->s2k[c];
    pack[k][cc2-oset] = var[c];
  }
}

/* END pack_ugrid2()                                                 */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Packs sediment sparse array into consecutive wet cells for each   */
/* layer.                                                            */
/*-------------------------------------------------------------------*/
void pack_ugrids(int *map, int sednz, int mapsize, double **var, double **pack, int oset)
{
  int c, cc, k;

  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    for (k = 0; k < sednz; k++) {
      pack[k][cc-oset] = var[k][c];
    }
  }
}

/* END pack_ugrids()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the integer sparse array into consecutive wet cells.        */
/* Note: the index values (var) are decremented by oset in order to  */
/* correctly visualize ugrid files. These indices are currently not  */
/* re-read from file: if they are then the offset must be accounted  */
/* for.                                                              */
/*-------------------------------------------------------------------*/
void pack_ugrid_i2(int *map, int size, int np, int **var, int **pack, int oset)
{
  int c, cc, i;
  int os = (oset) ? 0 : 1;

  for (cc = 1; cc <= size; cc++) {
    c = (map == NULL) ? cc : map[cc];
    for (i = 1; i <= np; i++) {
      pack[i-oset][cc-oset] = (var[i][c] > 0) ? var[i][c] - oset : Filli;
    }
  }
}

void pack_ugrid_ri2(int *map, int size, int np, int **var, int **pack, int oset)
{
  int c, cc, i;
  int os = (oset) ? 0 : 1;

  for (cc = 1; cc <= size; cc++) {
    c = (map == NULL) ? cc : map[cc];
    for (i = 1; i <= np; i++) {
      pack[cc-oset][i-oset] = (var[i][c] > 0) ? var[i][c] - oset : Filli;
    }
  }
}


/* END pack_ugrid_i2()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the integer sparse array into consecutive wet cells.        */
/*-------------------------------------------------------------------*/
void pack_ugrid_i2a(int *map, int size, int np, int **var, int **pack, int oset)
{
  int c, cc, i, val;

  for (cc = 1; cc <= size; cc++) {
    for (i = 0; i < np; i++) {
      val = var[cc][i];
      if (map != NULL) val = map[val];
      pack[i][cc-oset] = (val > 0) ? val - oset : Filli;
    }
  }
}

void pack_ugrid_ri2a(int *map, int size, int np, int **var, int **pack, int oset)
{
  int c, cc, i, val;

  for (cc = 1; cc <= size; cc++) {
    for (i = 0; i < np; i++) {
      val = var[cc][i];
      if (map != NULL) val = map[val];
      pack[cc-oset][i] = (val > 0) ? val - oset : Filli;
    }
  }
}

/* END pack_ugrid_i2b()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes 3D scalars on edges.                                     */
/*-------------------------------------------------------------------*/
void vel_scalar(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset) 
{
  master_t *master =  dumpdata->master;
  geometry_t *geom = master->geom;
  double vel;
  int c, cc, cm, ci, e, ee, c2, cc2, e2, k;

  for (ee = 1; ee <= dumpdata->nu1; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      ci = geom->wgst[c];
      c2 = geom->m2d[c];
      e2 = dumpdata->u1e[c2];
      k = (ci) ? geom->s2k[ci] : geom->s2k[c];
      p1[k][e2-oset] = u1[c];
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      ci = geom->wgst[c];
      c2 = geom->m2d[c];
      e2 = dumpdata->u2e[c2];
      k = (ci) ? geom->s2k[ci] : geom->s2k[c];
      p1[k][e2-oset] = u2[c];
    }
  }
}

/* END vel_scaler()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes 3D edge vector at edges and  tangential to edges.        */
/*-------------------------------------------------------------------*/
void vel_tan3D(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset) 
{
  master_t *master =  dumpdata->master;
  geometry_t *geom = master->geom;
  double vel;
  int c, cc, cm, ci, e, ee, e2, c2, k;

  for (ee = 1; ee <= dumpdata->nu1; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      ci = geom->wgst[c];
      c2 = geom->m2d[c];
      e2 = dumpdata->u1e[c2];
      k = (ci) ? geom->s2k[ci] : geom->s2k[c];
      cm = dumpdata->e2c[e][0];
      vel = 0.25 * (u2[c] + u2[geom->yp1[c]] +
		    u2[cm] + u2[geom->yp1[cm]]);

      p1[k][e2-oset] = u1[c];
      p2[k][e2-oset] = vel;
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      ci = geom->wgst[c];
      c2 = geom->m2d[c];
      e2 = dumpdata->u2e[c2];
      k = (ci) ? geom->s2k[ci] : geom->s2k[c];
      cm = dumpdata->e2c[e][0];
      vel = 0.25 * (u1[c] + u1[geom->xp1[c]] +
		    u1[cm] + u1[geom->xp1[cm]]);
      p1[k][e2-oset] = u2[c];
      p2[k][e2-oset] = vel;
    }
  }
}

/* END vel_tan3D()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes 2D edge vector at edges and  tangential to edges.        */
/*-------------------------------------------------------------------*/
void vel_tan2D(dump_data_t *dumpdata, double *u1, double *u2, double *p1, double *p2,  int oset) 
{
  master_t *master =  dumpdata->master;
  geometry_t *geom = master->geom;
  int c, cc, cm, e, ee;

  for (ee = 1; ee <= dumpdata->nu1S; ee++) {
    e = dumpdata->w2_e1[ee];
    if (e > 0) {
      c = dumpdata->e2c[e][1];
      cm = dumpdata->e2c[e][0];
      p1[e-oset] = u1[c];
      p2[e-oset] = 0.25 * (u2[c] + u2[geom->yp1[c]] +
			   u2[cm] + u2[geom->yp1[cm]]);
    } else {
      e = abs(e);
      c = dumpdata->e2c[e][1];
      cm = dumpdata->e2c[e][0];
      p1[e-oset] = u2[c];
      p2[e-oset] = 0.25 * (u1[c] + u1[geom->xp1[c]] +
			   u1[cm] + u1[geom->xp1[cm]]);
    }
  }
}

/* END vel_tan2D()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes 3D cell centered velocity.                               */
/*-------------------------------------------------------------------*/
void vel_cen3D(dump_data_t *dumpdata, double *u1, double *u2, double **p1, double **p2,  int oset) 
{
  master_t *master =  dumpdata->master;
  geometry_t *geom = master->geom;
  int cc;

  for (cc = 1; cc <= geom->n3_t; cc++) {
    int c = geom->w3_t[cc];
    int ci = geom->wgst[c];
    int xp = geom->xp1[c];
    int yp = geom->yp1[c];
    int c2 = geom->m2d[c];
    int cc2 = geom->c2cc[c2];
    int k = (ci) ? geom->s2k[ci] : geom->s2k[c];
    int xp2 = geom->m2d[xp];
    int yp2 = geom->m2d[yp];
    double u1val = (u1[c] + u1[xp]) / 2.0;
    double u2val = (u2[c] + u2[yp]) / 2.0;
    double sinth =
      (sin(geom->thetau1[c2]) + sin(geom->thetau1[xp2]) +
       sin(geom->thetau2[c2]) + sin(geom->thetau2[yp2])) / 4.0;
    double costh =
      (cos(geom->thetau1[c2]) + cos(geom->thetau1[xp2]) +
       cos(geom->thetau2[c2]) + cos(geom->thetau2[yp2])) / 4.0;
    p1[k][cc2-oset] = u1val * costh - u2val * sinth;
    p2[k][cc2-oset] = u1val * sinth + u2val * costh;
  }
}

/* END vel_cen3D()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes 2D cell centered velocity.                               */
/*-------------------------------------------------------------------*/
void vel_cen2D(dump_data_t *dumpdata, double *u1, double *u2, double *p1, double *p2,  int oset) 
{
  master_t *master =  dumpdata->master;
  geometry_t *geom = master->geom;
  int cc;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    int c = geom->w2_t[cc];
    int xp = geom->xp1[c];
    int yp = geom->yp1[c];
    int c2 = geom->m2d[c];
    int xp2 = geom->m2d[xp];
    int yp2 = geom->m2d[yp];
    double u1val = (u1[c] + u1[xp]) / 2.0;
    double u2val = (u2[c] + u2[yp]) / 2.0;
    double sinth =
      (sin(geom->thetau1[c2]) + sin(geom->thetau1[xp2]) +
       sin(geom->thetau2[c2]) + sin(geom->thetau2[yp2])) / 4.0;
    double costh =
      (cos(geom->thetau1[c2]) + cos(geom->thetau1[xp2]) +
       cos(geom->thetau2[c2]) + cos(geom->thetau2[yp2])) / 4.0;
    p1[cc-oset] = u1val * costh - u2val * sinth;
    p2[cc-oset] = u1val * sinth + u2val * costh;
  }
}

/* END vel_cen2D()                                                   */
/*-------------------------------------------------------------------*/

static void write_mean_atts(dump_data_t *dumpdata, int fid)
{
  master_t *master =  dumpdata->master;
  char buf[MAXSTRLEN], mcs[MAXSTRLEN];
  int i;

  if (!(master->means & NONE)) {
    nc_put_att_double(fid, NC_GLOBAL, "mean_c",
		      NC_DOUBLE, 1, &master->meanc[1]);
    nc_put_att_double(fid, NC_GLOBAL, "mean_next",
		      NC_DOUBLE, 1, &master->means_next);
    if (master->means_dt == SEASONAL || master->means_dt == MONTHLY || master->means_dt == DAILY) {
      sprintf(mcs, "%f ", master->meancs[1]);
      if (master->means_dt == SEASONAL || master->means_dt == MONTHLY) {
	for (i = 2; i <= 12; i++) {
	  sprintf(buf, "%f ", master->meancs[i]);
	  strcat(mcs, buf);
	}
      }
      if (master->means_dt == DAILY) {
	for (i = 2; i <= 365; i++) {
	  sprintf(buf, "%f ", master->meancs[i]);
	  strcat(mcs, buf);
	}
      }
      nc_put_att_text(fid, NC_GLOBAL, "mean_mc", strlen(mcs), mcs);
    }
  }
}


void set_dims(int dimf, int *dims, int d1, int d2)
{
  if (dimf) {
    dims[0] = d1;
    dims[1] = d2;
  } else {
    dims[0] = d2;
    dims[1] = d1;
  }
}

void set_count(int dimf, size_t *count, int d1, int d2)
{
  if (dimf) {
    count[0] = d1;
    count[1] = d2;
  } else {
    count[0] = d2;
    count[1] = d1;
  }
}
