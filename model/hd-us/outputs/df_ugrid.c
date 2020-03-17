/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/df_ugrid.c
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
 *  $Id: df_ugrid.c 6403 2019-11-21 23:01:37Z her127 $
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
void pack_ugrid2(int *hmap, int *vmap, int *m2d, int mapsize, double *var, double **pack, int oset);
void pack_ugrid3(int *map, int mapsize, double *var, double *pack, int oset);
void pack_ugrids(int sednz, int mapsize, double **var, double **pack, int oset);
void pack_ugridi(int size, int *var, int *pack, int oset);
void pack_ugrid_i2(int size, int np, int **var, int **pack, int oset);
void pack_ugrid_i2a(int size, int np, int **var, int **pack, int oset);
void pack_ugrid_ri2(int size, int np, int **var, int **pack, int oset);
void pack_ugrid_ri2a(int size, int np, int **var, int **pack, int oset);

void check_window_map_us(geometry_t **window, char *name);

#define UGRID_ALL_VARS "u1av u2av uav vav wind1 wtop eta patm u1 u2 u v w dens dens_0 Kz Vz Cd u1kh topz "

double percentiles1[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
			0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
			0.90, 0.95, 1.00};
#define Ns1 (sizeof(percentiles1) / sizeof(double))

void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
			      nc_type fptype, char *modulo);
static void df_ugrid_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_ugrid3_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_ugrid_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_ugrid_get_dimids(int cdfid, df_ugrid_var_t var, int *ns, int *zid);
static void df_ugrid_get_dimsizes(dump_file_t *df, df_ugrid_var_t var, size_t *sz, size_t *ns, size_t *nz);
static void df_ugrid3_get_dimids(int cdfid, df_ugrid_var_t var, int *ns);
static void df_ugrid3_get_dimsizes(dump_file_t *df, df_ugrid_var_t var, size_t *ns);
static void write_mean_atts(dump_data_t *dumpdata, int fid);
int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref);
void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
				 nc_type fptype, char *modulo)
{
  /* dimension ids */
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

  vid = ncw_var_id(cdfid, "Mesh2_index");
  write_text_att(cdfid, vid, "long_name", "Connectivity to unstructured index");
  write_text_att(cdfid, vid, "coordinates", "Mesh2_face_x, Mesh2_face_y");
  if (geom->us_type & US_IJ) {
    vid = ncw_var_id(cdfid, "Mesh2_iindex");
    write_text_att(cdfid, vid, "long_name", "Connectivity to i index");
    write_text_att(cdfid, vid, "coordinates", "Mesh2_face_x, Mesh2_face_y");

    vid = ncw_var_id(cdfid, "Mesh2_jindex");
    write_text_att(cdfid, vid, "long_name", "Connectivity to j index");
    write_text_att(cdfid, vid, "coordinates", "Mesh2_face_x, Mesh2_face_y");
  }

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

  if (dumpdata->nvars) {
    int i, j, dims[10];
    char buf[MAXSTRLEN];
    for (i = 0; i < dumpdata->nvars; ++i) {
      for (j = 0; j < Ns1; ++j) {
	sprintf(buf, "%s_%3.3d", dumpdata->vars[i], (int)(percentiles1[j]*100));
	dims[0] = ncw_dim_id(cdfid, "Mesh2_layers");
	dims[1] = ncw_dim_id(cdfid, "nMesh2_face");
	/*
	  dims[0] = geom->nz;
	  dims[1] = geom->szcS;
	*/
	dims[2] = 0;
	nc_def_var(cdfid, buf, NC_DOUBLE, 2, dims, &vid);
	write_text_att(cdfid, vid, "units", "");
	sprintf(buf, "%s %g percentile", dumpdata->vars[i], percentiles1[j]);
	write_text_att(cdfid, vid,"long_name", buf);
	write_text_att(cdfid, vid, "coordinates", "Mesh2_face_x, Mesh2_face_y, Mesh2_layers");
      }
    }
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
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);
  write_mean_atts(dumpdata, cdfid);
}

void write_dump_attributes_ugrid3(dump_data_t *dumpdata, int cdfid,
				 nc_type fptype, char *modulo)
{
  /* dimension ids */
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

  vid = ncw_var_id(cdfid, "Mesh3");
  write_text_att(cdfid, vid, "cf_role", "mesh_topology");
  write_text_att(cdfid, vid, "long_name", "Topology data of 2D unstructured mesh");
  nc_put_att_int(cdfid, vid, "topology_dimension", NC_INT, 1, &two);
  write_text_att(cdfid, vid, "node_coordinates", "Mesh3_node_x Mesh3_node_y") ;
  write_text_att(cdfid, vid, "face_dimension", "nMesh3_face");
  write_text_att(cdfid, vid, "edge_dimension", "nMesh3_edge");
  write_text_att(cdfid, vid, "volume_dimension", "nMesh3_vol");
  write_text_att(cdfid, vid, "edge_coordinates", "Mesh3_edge_x Mesh3_edge_y");
  write_text_att(cdfid, vid, "face_coordinates", "Mesh3_face_x Mesh3_face_y");

  /* time independent variables */
  vid = ncw_var_id(cdfid, "Mesh3_layerfaces");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer faces");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "Mesh3_layers");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer centre");
  write_text_att(cdfid, vid, "coordinate_type", "Z");
  write_text_att(cdfid, vid, "positive", "up");

  vid = ncw_var_id(cdfid, "Mesh3_node_x");
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

  vid = ncw_var_id(cdfid, "Mesh3_node_y");
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

  vid = ncw_var_id(cdfid, "Mesh3_face_x");
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

  vid = ncw_var_id(cdfid, "Mesh3_face_y");
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

  vid = ncw_var_id(cdfid, "Mesh3_face_xbnds");
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

  vid = ncw_var_id(cdfid, "Mesh3_face_ybnds");
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

  vid = ncw_var_id(cdfid, "Mesh3_edge_x");
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

  vid = ncw_var_id(cdfid, "Mesh3_edge_y");
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

  vid = ncw_var_id(cdfid, "Mesh3_depth");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at sea-bed at cell centre");

  vid = ncw_var_id(cdfid, "h1au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Edge length between nodes");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_edge_x, Mesh3_edge_y");

  vid = ncw_var_id(cdfid, "h2au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Length between face centres");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_edge_x, Mesh3_edge_y");

  vid = ncw_var_id(cdfid, "thetau1");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation of normal vector");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_edge_x, Mesh3_edge_y");

  vid = ncw_var_id(cdfid, "coriolis");
  write_text_att(cdfid, vid, "units", "s-1");
  write_text_att(cdfid, vid, "long_name", "Coriolis parameter");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_face_x, Mesh3_face_y");

  vid = ncw_var_id(cdfid, "Mesh3_index");
  write_text_att(cdfid, vid, "long_name", "Connectivity to unstructured index");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_face_x, Mesh3_face_y");
  if (geom->us_type & US_IJ) {
    vid = ncw_var_id(cdfid, "Mesh3_iindex");
    write_text_att(cdfid, vid, "long_name", "Connectivity to i index");
    write_text_att(cdfid, vid, "coordinates", "Mesh3_face_x, Mesh3_face_y");

    vid = ncw_var_id(cdfid, "Mesh3_jindex");
    write_text_att(cdfid, vid, "long_name", "Connectivity to j index");
    write_text_att(cdfid, vid, "coordinates", "Mesh3_face_x, Mesh3_face_y");
  }
  vid = ncw_var_id(cdfid, "Mesh3_f2k");
  write_text_att(cdfid, vid, "long_name", "Face to layer map for 3D variables");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_face_x, Mesh3_face_y");
  vid = ncw_var_id(cdfid, "Mesh3_e2k");
  write_text_att(cdfid, vid, "long_name", "Edge to layer map for 3D variables");
  write_text_att(cdfid, vid, "coordinates", "Mesh3_edge_x, Mesh3_edge_y");


  /* time dependent variables */

  vid = ncw_var_id(cdfid, "t");
  write_text_att(cdfid, vid, "units", dumpdata->output_tunit);
  write_text_att(cdfid, vid, "long_name", "Time");
  write_text_att(cdfid, vid, "coordinate_type", "time");
  if (strlen(modulo))
    write_text_att(cdfid, vid, "modulo", modulo);

  if ((vid = ncw_var_id(cdfid, "u1av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of depth averaged current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "u2av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "tangential component of depth averaged current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "uav")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "east component of depth averaged current at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "vav")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "north component of depth averaged current at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "wtop")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "Vertical velocity at surface");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "eta")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Surface Elevation");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "topz")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate for surface cell");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "wind1")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of wind stress at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y");
  }

  if ((vid = ncw_var_id(cdfid, "windx")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "east component of wind stress at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "windy")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "north component of wind stress at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "patm")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "u1")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "normal component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "u2")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "tangential component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "umean")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "mean normal component of current at edge");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_edge_x, Mesh3_edge_y, Mesh3_layers");
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
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "v")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "north component of current at face centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "w")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "K component of current at cell centre and Z grid");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "dens")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Density");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "dens_0")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Potential density");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
    fill = 1025.0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Kz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Kz");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Vz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Vz");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "ptconc")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Particle concentration");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "Cd")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Bottom drag coefficient");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "u1kh")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre2 second-1");
    write_text_att(cdfid, vid, "long_name",
                   "horizontal diffusivity at centre");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers");
  }

  if ((vid = ncw_var_id(cdfid, "ustrcw")) >= 0) {
    write_text_att(cdfid, vid, "units", "ms-1");
    write_text_att(cdfid, vid, "long_name", "Bottom friction velocity");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh3_face_x, Mesh3_face_y");
  }


  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "t, Mesh3_face_x, Mesh3_face_y, Mesh3_layers"}
    };
    tracer_write_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "t, Mesh3_face_x, Mesh3_face_y"}
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
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);

}


void df_ugrid_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_ugrid_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t <= (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "t", df->tout);
}

void *df_ugrid_create(dump_data_t *dumpdata, dump_file_t *df)
{
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
  int nface2 = dumpdata->nface2 - dumpdata->start_index;
  int nedge2 = dumpdata->nedge2 - dumpdata->start_index;
  int nvertex2 = dumpdata->nvertex2 - dumpdata->start_index;
  int npem = dumpdata->npe - dumpdata->start_index;
  int kcentreid_sed;            /* K dimension id at grid centre for
                                   sediments */
  int kgridid_sed;              /* K dimension id at grid corner for
                                   sediments */
  FILE *fp;
  int nc_mode;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;

  df_ugrid_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);

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
  nc_def_dim(cdfid, "nMesh2_face_links", geom->v2_e1, &nMesh2_face_linksid);
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

  /* Grid index (optional) */
  dims[0] = nMesh2_faceid;
  nc_def_var(cdfid, "Mesh2_index", NC_INT, 1, dims, &vid);
  if (geom->us_type & US_IJ) {
    nc_def_var(cdfid, "Mesh2_iindex", NC_INT, 1, dims, &vid);
    nc_def_var(cdfid, "Mesh2_jindex", NC_INT, 1, dims, &vid);
  }
  
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


void *df_ugrid3_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int n;
  int cdfid;
  int dims[10];
  int vid;                      /* df_variable_t dimension */
  /* dimension ids */
  int recdimid;                 /* record dimension id */

  int nMesh3_nodeid;            /* Node dimension */
  int nMesh3_edgeid;            /* Edge dimension */
  int nMesh3_faceid;            /* Face dimension */
  int nMesh3_volid;             /* Volume dimension */
  int nMesh3_vol_nodeid;        /* 3D node dimension */
  int nMesh3_vol_edgeid;        /* 3D edge dimension */
  int nMaxMesh3_face_nodesid;   /* Max # nodes per face */
  int Mesh3_layersid;           /* Number of layers */ 
  int Mesh3_layerfacesid;       /* Number of layer faces */ 
  int Twoid;                    /* Id for 2 */
  int nface2 = dumpdata->nface2 - dumpdata->start_index;
  int nedge2 = dumpdata->nedge2 - dumpdata->start_index;
  int nvertex2 = dumpdata->nvertex2 - dumpdata->start_index;
  int nface3 = dumpdata->nface3 - dumpdata->start_index;
  int nedge3 = dumpdata->nedge3 - dumpdata->start_index;
  int nvertex3 = dumpdata->nvertex3 - dumpdata->start_index;
  int npem = dumpdata->npe - dumpdata->start_index;

  FILE *fp;
  int nc_mode;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;

  df_ugrid_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);

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

  nc_def_dim(cdfid, "Mesh3_layers", df->nz, &Mesh3_layersid);
  nc_def_dim(cdfid, "Mesh3_layerfaces", df->nz + 1, &Mesh3_layerfacesid);
  nc_def_dim(cdfid, "nMesh3_node", nvertex2, &nMesh3_nodeid);
  nc_def_dim(cdfid, "nMesh3_edge", nedge2, &nMesh3_edgeid);
  nc_def_dim(cdfid, "nMesh3_face", nface2, &nMesh3_faceid);
  nc_def_dim(cdfid, "nMesh3_vol", nface3, &nMesh3_volid);
  nc_def_dim(cdfid, "nMesh3_vol_node", nvertex3, &nMesh3_vol_nodeid);
  nc_def_dim(cdfid, "nMesh3_vol_edge", nedge3, &nMesh3_vol_edgeid);
  nc_def_dim(cdfid, "nMaxMesh3_face_nodes", npem, &nMaxMesh3_face_nodesid);
  nc_def_dim(cdfid, "Two", 2, &Twoid);

  /* Grid */
  nc_def_var(cdfid, "Mesh3", NC_INT, 0, NULL, &vid);

  /* Mesh node coordinates */
  dims[0] = nMesh3_nodeid;
  nc_def_var(cdfid, "Mesh3_node_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh3_node_y", NC_DOUBLE, 1, dims, &vid);

  /* Mesh face and edge coordinate variables */
  dims[0] = nMesh3_faceid;
  nc_def_var(cdfid, "Mesh3_face_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh3_face_y", NC_DOUBLE, 1, dims, &vid);
  set_dims(dumpdata->face_dim, dims, nMaxMesh3_face_nodesid, nMesh3_faceid);
  nc_def_var(cdfid, "Mesh3_face_xbnds", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "Mesh3_face_ybnds", NC_DOUBLE, 2, dims, &vid);
  dims[0] = nMesh3_edgeid;
  nc_def_var(cdfid, "Mesh3_edge_x", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "Mesh3_edge_y", NC_DOUBLE, 1, dims, &vid);

  /* Vertical structure */
  dims[0] = Mesh3_layersid;
  nc_def_var(cdfid, "Mesh3_layers", NC_DOUBLE, 1, dims, &vid);
  dims[0] = Mesh3_layerfacesid;
  nc_def_var(cdfid, "Mesh3_layerfaces", NC_DOUBLE, 1, dims, &vid);
  dims[0] = nMesh3_faceid;
  nc_def_var(cdfid, "Mesh3_depth", NC_DOUBLE, 1, dims, &vid);

  /* Cell lengths */
  dims[0] = nMesh3_edgeid;
  nc_def_var(cdfid, "h2au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h1au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "thetau1", NC_DOUBLE, 1, dims, &vid);
  dims[0] = nMesh3_faceid;
  nc_def_var(cdfid, "coriolis", fptype, 1, dims, &vid);

  /* Grid index (optional) */
  dims[0] = nMesh3_faceid;
  nc_def_var(cdfid, "Mesh3_index", NC_INT, 1, dims, &vid);
  dims[0] = nMesh3_volid;
  nc_def_var(cdfid, "Mesh3_f2k", NC_DOUBLE, 1, dims, &vid);
  dims[0] = nMesh3_vol_edgeid;
  nc_def_var(cdfid, "Mesh3_e2k", NC_DOUBLE, 1, dims, &vid);
  if (geom->us_type & US_IJ) {
    nc_def_var(cdfid, "Mesh3_iindex", NC_INT, 1, dims, &vid);
    nc_def_var(cdfid, "Mesh3_jindex", NC_INT, 1, dims, &vid);
  }

  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid, "t", NC_DOUBLE, 1, dims, &vid);

  df_ugrid_init_data(dumpdata, df, cdfid);
  data = (df_ugrid_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    df_ugrid3_get_dimids(cdfid, data->vars[n], &dims[1]);
    ncw_def_var2(df->name,cdfid, df->vars[n], data->vars[n].type, 2, dims, &vid, df->compress);
  }

  write_dump_attributes_ugrid3(dumpdata, cdfid, fptype, df->modulo);

  nc_enddef(cdfid);

  df_ugrid3_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
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
    count[2] = df->nface2;
    pack_ugrids(dumpdata->sednz - 1, count[2], geom->cellz_sed, dumpdata->wc, oset);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_layers_sed"), start, count, dumpdata->wc);

    count[1] = dumpdata->sednz + 1;
    pack_ugrids(dumpdata->sednz, count[2], geom->gridz_sed, dumpdata->wc, oset);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_layerfaces_sed"), start, count, dumpdata->wc);
  }

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_ugrid_var_t vn = data->vars[n];
    void *p = vn.v;
    start[1] = df->ilower;

    if (vn.ndims == 1) {
      memset(dumpdata->w1s,  0, geom->szcS*sizeof(double));
      df_ugrid_get_dimsizes(df, vn, NULL, &count[1], NULL);

      pack_ugrid1(vn.hmap, count[1], (*(double **)p), dumpdata->w1s, oset);
      nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count,
		       dumpdata->w1s);
      /*
      nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count,
		       (*(double **)p));
      */
    }  else if (vn.ndims == 2) {
      size_t sz;
      int cc,c,k, kb;
      double **w = dumpdata->wc;

      if (!vn.sediment) {
	if (vn.xylocation & CL_EDGE) w = dumpdata->we;
	start[1] = df->klower;
	df_ugrid_get_dimsizes(df, vn, &sz, &count[2], &count[1]);
	
	/* Initialize land cells */
	for (k = dumpdata->nz-1; k >= 0; k--)
	  for (cc = 0; cc < count[2]; cc++)
	    w[k][cc] = NaN;
	pack_ugrid2(vn.hmap, vn.vmap, vn.m2d, sz, (*(double **)p), w, oset);
	nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count,
			 w);
      } else {
	start[1] = 0;
	df_ugrid_get_dimsizes(df, vn, NULL, &count[2], &count[1]);
	/* Initialize land cells */
	for (k = dumpdata->sednz-1; k >= 0; k--)
	  for (cc = 0; cc < count[2]; cc++)
	    w[k][cc] = NaN;
	pack_ugrids(dumpdata->sednz, count[2], (*(double ***)p), w, oset);
	nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, w);
      }
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}


void df_ugrid3_write(dump_data_t *dumpdata, dump_file_t *df, double t)
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

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_ugrid_var_t vn = data->vars[n];
    void *p = vn.v;
    start[1] = df->ilower;

    memset(dumpdata->w1,  0, geom->szm*sizeof(double));
    df_ugrid3_get_dimsizes(df, vn, &count[1]);
    pack_ugrid3(vn.hmap, count[1], (*(double **)p), dumpdata->w1, oset);
    nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count,
		       dumpdata->w1);
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}


void df_ugrid_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}


static void df_ugrid_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int fid = data->fid;
  int c, cc, vv, v;
  double **xbnds, **ybnds;
  int oset = (dumpdata->start_index) ? 0 : 1;
  void (*i2) (int size, int np, int **var, int **pack, int oset) = NULL;
  void (*i2a) (int size, int np, int **var, int **pack, int oset) = NULL;
  int n;

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
  i2(geom->b2_t, geom->npem,  geom->c2v, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_nodes"), start, count, dumpdata->i2s);

  set_count(dumpdata->edge_dim, count, 2, dumpdata->nedge2);
  i2a(geom->n2_e1, 2, geom->e2v, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_edge_nodes"), start, count, dumpdata->i2s);

  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  i2(geom->b2_t, geom->npem, geom->c2e, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_edges"), start, count, dumpdata->i2s);

  set_count(dumpdata->edge_dim, count, 2, geom->v2_e1+dumpdata->start_index);
  i2a(geom->v2_e1, 2, geom->e2c, dumpdata->i2s, oset);
  nc_i_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_links"), start, count, dumpdata->i2s);

  count[0] = dumpdata->nvertex2;
  pack_ugrid1(geom->w2_e2, geom->n2_e2, geom->gridx, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_node_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e2, geom->n2_e2, geom->gridy, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_node_y"), start, count, dumpdata->w1s);

  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->cellx, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_face_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->celly, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_face_y"), start, count, dumpdata->w1s);

  count[0] = dumpdata->nedge2;
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->u1x, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_edge_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->u1y, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_edge_y"), start, count, dumpdata->w1s);

  if (dumpdata->face_dim) {
    xbnds = d_alloc_2d(geom->szcS, geom->npem+1);
    ybnds = d_alloc_2d(geom->szcS, geom->npem+1);
  } else {
    xbnds = d_alloc_2d(geom->npem+1, geom->szcS);
    ybnds = d_alloc_2d(geom->npem+1, geom->szcS);
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (vv = 1; vv <= geom->npem; vv++) {
      v = geom->c2v[vv][c];
      if (dumpdata->face_dim) {
	xbnds[vv-oset][cc-oset] = (vv <= geom->npe[c]) ? geom->gridx[v] : Filld;
	ybnds[vv-oset][cc-oset] = (vv <= geom->npe[c]) ? geom->gridy[v] : Filld;
      } else {
	xbnds[cc-oset][vv-oset] = (vv <= geom->npe[c]) ? geom->gridx[v] : Filld;
	ybnds[cc-oset][vv-oset] = (vv <= geom->npe[c]) ? geom->gridy[v] : Filld;
      }
    }
  }
  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_xbnds"), start, count, xbnds);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh2_face_ybnds"), start, count, ybnds);
  d_free_2d(xbnds);
  d_free_2d(ybnds);

  count[0] = dumpdata->nedge2;
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->h1au1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au1"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->h2au1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au1"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->thetau1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau1"), start, count, dumpdata->w1s);


  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, master->coriolis, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "coriolis"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->botz, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh2_depth"), start, count, dumpdata->w1s);
  count[0] = dumpdata->nface2;
  pack_ugridi(geom->b2_t, geom->c2cc, dumpdata->i1s, oset);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh2_index"), start, count, dumpdata->i1s);
  if (geom->us_type & US_IJ) {
  pack_ugridi(geom->b2_t, geom->s2i, dumpdata->i1s, oset);
    nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh2_iindex"), start, count, dumpdata->i1s);
  pack_ugridi(geom->b2_t, geom->s2j, dumpdata->i1s, oset);
    nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh2_jindex"), start, count, dumpdata->i1s);
  }
  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}


static void df_ugrid3_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int fid = data->fid;
  int c, cc, vv, v;
  double **xbnds, **ybnds;
  int oset = (dumpdata->start_index) ? 0 : 1;
  void (*i2) (int size, int np, int **var, int **pack, int oset) = NULL;
  void (*i2a) (int size, int np, int **var, int **pack, int oset) = NULL;
  int n;

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
  nc_put_vara_double(fid, ncw_var_id(fid, "Mesh3_layerfaces"), start, count,
                     dumpdata->gridz);
  count[0] = dumpdata->nz;
  nc_put_vara_double(fid, ncw_var_id(fid, "Mesh3_layers"), start, count,
                     dumpdata->cellz);


  count[0] = dumpdata->nvertex2;
  pack_ugrid1(geom->w2_e2, geom->n2_e2, geom->gridx, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_node_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e2, geom->n2_e2, geom->gridy, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_node_y"), start, count, dumpdata->w1s);

  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->cellx, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_face_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->celly, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_face_y"), start, count, dumpdata->w1s);

  count[0] = dumpdata->nedge2;
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->u1x, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_edge_x"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->u1y, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_edge_y"), start, count, dumpdata->w1s);

  if (dumpdata->face_dim) {
    xbnds = d_alloc_2d(geom->szcS, geom->npem+1);
    ybnds = d_alloc_2d(geom->szcS, geom->npem+1);
  } else {
    xbnds = d_alloc_2d(geom->npem+1, geom->szcS);
    ybnds = d_alloc_2d(geom->npem+1, geom->szcS);
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (vv = 1; vv <= geom->npem; vv++) {
      v = geom->c2v[vv][c];
      if (dumpdata->face_dim) {
	xbnds[vv-oset][cc-oset] = (vv <= geom->npe[c]) ? geom->gridx[v] : Filld;
	ybnds[vv-oset][cc-oset] = (vv <= geom->npe[c]) ? geom->gridy[v] : Filld;
      } else {
	xbnds[cc-oset][vv-oset] = (vv <= geom->npe[c]) ? geom->gridx[v] : Filld;
	ybnds[cc-oset][vv-oset] = (vv <= geom->npe[c]) ? geom->gridy[v] : Filld;
      }
    }
  }
  set_count(dumpdata->face_dim, count, dumpdata->npe, dumpdata->nface2);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh3_face_xbnds"), start, count, xbnds);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "Mesh3_face_ybnds"), start, count, ybnds);
  d_free_2d(xbnds);
  d_free_2d(ybnds);

  count[0] = dumpdata->nedge2;
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->h1au1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au1"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->h2au1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au1"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_e1, geom->n2_e1, geom->thetau1, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau1"), start, count, dumpdata->w1s);


  count[0] = dumpdata->nface2;
  pack_ugrid1(geom->w2_t, geom->b2_t, master->coriolis, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "coriolis"), start, count, dumpdata->w1s);
  pack_ugrid1(geom->w2_t, geom->b2_t, geom->botz, dumpdata->w1s, oset);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_depth"), start, count, dumpdata->w1s);
  count[0] = dumpdata->nface2;
  pack_ugridi(geom->b2_t, geom->c2cc, dumpdata->i1s, oset);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh3_index"), start, count, dumpdata->i1s);
  if (geom->us_type & US_IJ) {
  pack_ugridi(geom->b2_t, geom->s2i, dumpdata->i1s, oset);
    nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh3_iindex"), start, count, dumpdata->i1s);
  pack_ugridi(geom->b2_t, geom->s2j, dumpdata->i1s, oset);
    nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh3_jindex"), start, count, dumpdata->i1s);
  }
  count[0] = dumpdata->nface3;
  pack_ugridi(geom->b3_t, geom->s2k, dumpdata->i1s, oset);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh3_f2k"), start, count, dumpdata->i1s);
  count[0] = dumpdata->nedge3;
  pack_ugridi(geom->b3_e1, geom->e2k, dumpdata->i1s, oset);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "Mesh3_e2k"), start, count, dumpdata->i1s);

  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}


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
  var->sediment = 0;
  var->hmap = geom->w2_t;
  var->vmap = geom->s2k;
  var->m2d = geom->m2d;

  if (strcmp(name, "u1av") == 0) {
    var->v = (void **)&master->u1av;
    var->xylocation = CL_SP2|CL_EDGE;
    var->hmap = geom->w2_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "u2av") == 0) {
    var->v = (void **)&master->u2av;
    var->xylocation = CL_SP2|CL_EDGE;
    var->hmap = geom->w2_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "uav") == 0) {
    var->v = (void *)&master->uav;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "vav") == 0) {
    var->v = (void *)&master->vav;
    var->xylocation = CL_SP2|CL_FACE;
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
    var->v = (void *)&master->wind1;
    var->xylocation = CL_SP2|CL_EDGE;
    var->hmap = geom->w2_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void *)&master->patm;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "u") == 0) {
    var->v = (void *)&master->u;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "v") == 0) {
    var->v = (void *)&master->v;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "u1") == 0) {
    var->v = (void *)&master->u1;
    var->xylocation = CL_SP3|CL_EDGE;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->hmap = geom->w3_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "u2") == 0) {
    var->v = (void *)&master->u2;
    var->xylocation = CL_SP3|CL_EDGE;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->hmap = geom->w3_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void *)&master->w;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "u1mean") == 0) {
    var->v = (void *)&master->u1m;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "u2mean") == 0) {
    var->v = (void *)&master->u2m;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "umean") == 0) {
    var->v = (void *)&master->ume;
    var->xylocation = CL_SP3|CL_EDGE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
  }

  else if (strcmp(name, "u1vmean") == 0) {
    var->v = (void *)&master->u1vm;
    var->xylocation = CL_SP3|CL_EDGE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
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

static void df_ugrid3_get_dimids(int cdfid, df_ugrid_var_t var, int *ns)
{
  if (var.xylocation & CL_FACE) {
    if (var.zlocation & (CL_CENTRE|CL_GRID)) 
      *ns = ncw_dim_id(cdfid, "nMesh3_vol");
    else
      *ns = ncw_dim_id(cdfid, "nMesh3_face");
  } else if (var.xylocation & CL_EDGE) {
    if (var.zlocation & (CL_CENTRE|CL_GRID)) 
      *ns = ncw_dim_id(cdfid, "nMesh3_vol_edge");
    else
      *ns = ncw_dim_id(cdfid, "nMesh3_edge");
  } else if (var.xylocation & CL_VERTEX) {
    if (var.zlocation & (CL_CENTRE|CL_GRID)) 
      *ns = ncw_dim_id(cdfid, "nMesh3_vol_node");
    else
      *ns = ncw_dim_id(cdfid, "nMesh2_node");
  } else
    hd_quit("df_ugrid3_get_dimids: can't find dimensions.\n");
}

static void df_ugrid_get_dimsizes(dump_file_t *df, 
				  df_ugrid_var_t var, 
				  size_t *sz,         /* 3D vector size           */
				  size_t * ns,        /* Surface (2D) vector size */
				  size_t *nz)         /* Number of layers         */
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


static void df_ugrid3_get_dimsizes(dump_file_t *df, 
				   df_ugrid_var_t var, 
				   size_t * ns)       /* Surface (2D) vector size */
{
  if (var.xylocation == (CL_SP2|CL_FACE))
    *ns = df->nface2;
  else if (var.xylocation == (CL_SP2|CL_EDGE))
    *ns = df->nedge2;
  else if (var.xylocation == (CL_SP2|CL_VERTEX))
    *ns = df->nvertex2;
  else if (var.xylocation == (CL_SP3|CL_FACE))
    *ns = df->nface3;
  else if (var.xylocation == (CL_SP3|CL_EDGE))
    *ns = df->nedge3;
  else if (var.xylocation == (CL_SP3|CL_VERTEX))
    *ns = df->nvertex3;
  else
    hd_quit("df_ugrid:df_sp_get_dimsizes error in xylocation\n");
}


/*-------------------------------------------------------------------*/
/* Packs the sparse array into consecutive wet cells. If the 'map'   */
/* is geom->wsa then the entire wet grid, boundaries included, is    */
/* packed.                                                           */
/* Note the the first null cell in the sparse array is removed so    */
/* that the dumped variable indies go from 0:ns-1 whereas the sparse */
/* indicies go from 1:ns with value[0]=0.                            */
/* Cells are written to file in the order of the sparse vector; this */
/* is wet cells surrouded by wet cells, wet cells adacent to land,   */
/* open boundary cells. Auxiliary and ghost cells not included.      */
/*-------------------------------------------------------------------*/
void pack_ugrid1(int *map, int mapsize, double *var, double *pack, int oset)
{
  int cc, c;

  for (cc = 1; cc <= mapsize; cc++) {
    /*c = map[cc];*/
    pack[cc-oset] = var[cc];
  }
}

/* END pack_ugrid1()                                                 */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Packs the sparse array into consecutive wet cells. If the 'map'   */
/* is geom->wsa then the entire wet grid, boundaries included, is    */
/* packed.                                                           */
/* Note the the first null cell in the sparse array is removed so    */
/* that the dumped variable indies go from 0:ns-1 whereas the sparse */
/* indicies go from 1:ns with value[0]=0.                            */
/*-------------------------------------------------------------------*/
void pack_ugrid2(int *hmap, int *vmap, int *m2d, int mapsize, double *var, double **pack, int oset)
{
  int cc, c, c2, k;

  for (cc = 1; cc <= mapsize; cc++) {
    c = hmap[cc];
    c2 = m2d[c];
    k = vmap[c];
    pack[k][c2-oset] = var[c];
  }
}

/* END pack_ugrid2()                                                 */
/*-------------------------------------------------------------------*/

void pack_ugrids(int sednz, int mapsize, double **var, double **pack, int oset)
{
  int cc, k;

  for (cc = 1; cc <= mapsize; cc++) {
    for (k = 0; k < sednz; k++) {
      pack[k][cc-oset] = var[k][cc];
    }
  }
}

void pack_ugrid3(int *map, int mapsize, double *var, double *pack, int oset)
{
  int cc, c;

  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    pack[cc-oset] = var[c];
  }
}

/*-------------------------------------------------------------------*/
/* Packs the integer sparse array into consecutive wet cells.        */
/*-------------------------------------------------------------------*/
void pack_ugridi(int size, int *var, int *pack, int oset)
{
  int c;

  for (c = 1; c <= size; c++) {
    pack[c-oset] = var[c];
  }
}

/* END pack_ugridi()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the integer sparse array into consecutive wet cells.        */
/* Note: the index values (var) are decremented by oset in order to  */
/* correctly visualize ugrid files. These indices are currently not  */
/* re-read from file: if they are then the offset must be accounted  */
/* for.                                                              */
/*-------------------------------------------------------------------*/
void pack_ugrid_i2(int size, int np, int **var, int **pack, int oset)
{
  int c, i;
  int os = (oset) ? 0 : 1;

  for (c = 1; c <= size; c++) {
    for (i = 1; i <= np; i++) {
      pack[i-oset][c-oset] = (var[i][c] > 0) ? var[i][c] - oset : Filli;
    }
  }
}

void pack_ugrid_ri2(int size, int np, int **var, int **pack, int oset)
{
  int c, i;
  int os = (oset) ? 0 : 1;

  for (c = 1; c <= size; c++) {
    for (i = 1; i <= np; i++) {
      pack[c-oset][i-oset] = (var[i][c] > 0) ? var[i][c] - oset : Filli;
    }
  }
}


/* END pack_ugrid_i2()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the integer sparse array into consecutive wet cells.        */
/*-------------------------------------------------------------------*/
void pack_ugrid_i2a(int size, int np, int **var, int **pack, int oset)
{
  int c, i;

  for (c = 1; c <= size; c++) {
    for (i = 0; i < np; i++) {
      pack[i][c-oset] = (var[c][i] > 0) ? var[c][i] - oset : Filli;
    }
  }
}

void pack_ugrid_ri2a(int size, int np, int **var, int **pack, int oset)
{
  int c, i;

  for (c = 1; c <= size; c++) {
    for (i = 0; i < np; i++) {
      pack[c-oset][i] = (var[c][i] > 0) ? var[c][i] - oset : Filli;
    }
  }
}

/* END pack_ugrid_i2b()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Dumps window geometry to file                                     */
/*-------------------------------------------------------------------*/
void dump_windows_us(master_t *master, geometry_t **window, char *name, char *iname)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int cdfid, n, c;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  int dims[4];
  int vid;
  long t;
  /* dimension ids */
  int szcid;            /* 3D sparse dimension           */
  int szcSid;           /* 2D sparse dimension           */
  int szeid;            /* 3D sparse dimension           */
  int szeSid;           /* 2D sparse dimension           */
  int szvid;            /* 3D sparse dimension           */
  int szvSid;           /* 2D sparse dimension           */
  int kgridid;            /* K dimension id at grid corner */
  int oid, tid;
  int nwinsid;
  int nwins = geom->nwindows;
  char key[MAXSTRLEN];
  int *d1;
  int ct;
  int ncmode = overwrite_output() ? NC_CLOBBER:NC_NOCLOBBER;
  int npem, nvem, nvcm;
  int szcw, szcSw;
  int szew, szeSw;
  int szvw, szvSw;

  /*-----------------------------------------------------------------*/
  /* Create the netCDF file                                          */
#ifdef NC_NETCDF4
  ncmode |= NC_NETCDF4;
#else
  ncmode |= NC_64BIT_OFFSET;
#endif
  if (ncw_create(name, ncmode, &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", name);

  /*--------------------------*/
  /* Define global dimensions */
  /*--------------------------*/
  ncw_def_dim(name, cdfid, "szc",    geom->szc,  &szcid);
  ncw_def_dim(name, cdfid, "szcS",   geom->szcS, &szcSid);
  ncw_def_dim(name, cdfid, "sze",    geom->sze,  &szeid);
  ncw_def_dim(name, cdfid, "szeS",   geom->szeS, &szeSid);
  ncw_def_dim(name, cdfid, "szv",    geom->szv,  &szvid);
  ncw_def_dim(name, cdfid, "szvS",   geom->szvS, &szvSid);
  ncw_def_dim(name, cdfid, "k_grid",   geom->nz + 1, &kgridid);
  ncw_def_dim(name, cdfid, "nwp1",     nwins+1,      &nwinsid);

  /* Variables of length 1 */
  ncw_def_dim(name, cdfid, "one", 1, &oid);
  ncw_def_dim(name, cdfid, "two", 2, &tid);
  
  /* Global variables */
  dims[0] = kgridid;
  ncw_def_var(name, cdfid, "z_grid", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate at grid layer faces");  
  write_text_att(cdfid, vid, "coordinate_type", "Z");
  dims[0] = szcid;
  ncw_def_var(name, cdfid, "wn", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Global-local window map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "sc", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Global-local wet coordinate map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "ac", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Global-local auxiliary coordinate map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  dims[0] = oid;
  ncw_def_var(name, cdfid, "nz", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Number of layers");
  ncw_def_var(name, cdfid, "sednz", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Number of Sediment layers");

  /*
   * Window based sizes as variables, 1 per window
   */
  dims[0] = nwinsid;
  ncw_def_var(name, cdfid, "enon_w",  NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "enonS_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "ewet_w",  NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "ewetS_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "snon_w",  NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "snonS_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "szc_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "szcS_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "sze_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "szeS_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "szv_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "szvS_w", NC_INT, 1, dims, &vid);

  ncw_def_var(name, cdfid, "npem_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nvem_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nvcm_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "neem_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "nbpt_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte1S_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "n2_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "n3_t_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "n2_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "n3_e1_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "n2_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "n3_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "nm2s_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "ns2m_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "nbptS_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbe1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbe1S_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbe2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbe2S_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v2_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b2_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a2_t_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v3_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b3_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a3_t_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v2_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b2_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a2_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "x2_e1_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v2_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b2_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a2_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "x3_e1_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v3_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b3_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a3_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "ns2mS_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nm2sS_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "ns2me1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nm2se1_w", NC_INT, 1, dims, &vid);

  ncw_def_var(name, cdfid, "ns2me1S_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nm2se1S_w", NC_INT, 1, dims, &vid);
  
  /*
   * 2D - One 1D array per window
   *      Figure out a max size/window for each variable
   *
   * Define the master dimensions as we go along
   */
  dims[0] = nwinsid; // always

  /* CENTRES                                                         */ 
  /* window->npe[szcS]                                               */
  /* window->bot_t[szcS]                                             */
  /* window->sur_t[szcS]                                             */
  /* window->nsur_t[szcS]                                            */
  /* window->w2_t[szcS]                                              */
  /* window->zm1[szc]                                                */
  /* window->zp1[szc]                                                */
  /* window->m2d[szc]                                                */
  /* window->wsa[szc]                                                */
  /* window->wgst[szc]                                               */
  /* window->w3_t[szc]                                               */
  /* window->c2c[npem][szc]                                          */
  /* window->c2e[npem][szc]                                          */
  /* window->c2v[npem][szc]                                          */
  /* window->vIc[npem][szcS]                                         */
  /* window->eSc[npem][szcS]                                         */
  /* window->bpt[nbpt]                                               */
  /* window->bin[nbpt]                                               */
  /* window->bin2[nbpt]                                              */
  /* window->cellarea[szcS]                                          */
  
  /* Max 3D centres */
  ct = window[1]->szc;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szc > ct ? window[n]->szc : ct);
  ncw_def_dim(name, cdfid, "szc_max", ct, &szcw);

  /* Max 2D centres */
  ct = window[1]->szcS;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szcS > ct ? window[n]->szcS : ct);
  ncw_def_dim(name, cdfid, "szcS_max", ct, &szcSw);

  /* Max npe */
  ct = window[1]->npem + 1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->npem+1 > ct ? window[n]->npem+1 : ct);
  ncw_def_dim(name, cdfid, "npe_max", ct, &npem);

  /*
   * Surface sparse arrays
   */
  dims[1] = szcSw;
  ncw_def_var2(name, cdfid, "npe",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centres surrounding a cell");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  /* cell area */
  ncw_def_var2(name, cdfid, "cellarea", NC_DOUBLE, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell area");
  write_text_att(cdfid, vid, "units", "m2");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  
  /* Max a2_t */
  ct = window[1]->a2_t+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->a2_t+1 > ct ? window[n]->a2_t+1 : ct);
  ncw_def_dim(name, cdfid, "a2_t_max", ct, &dims[1]);
  // dims = a2_t
  ncw_def_var2(name, cdfid, "bot_t",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Bottom centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "sur_t",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Surface centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  
  ncw_def_var2(name, cdfid, "nsur_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Updated surface centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  /*
   * 3D sparse arrays
   */
  dims[1] = szcw;
  ncw_def_var2(name, cdfid, "wsa",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "m2d", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "zp1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centre map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "zm1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centre map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "wgst", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "lateral boundary ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  /* Maps */
  dims[1] = npem;
  dims[2] = szcw;
  ncw_def_var2(name, cdfid, "c2c", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "c2e", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "c2v", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  dims[1] = npem;
  dims[2] = szcSw;
  ncw_def_var2(name, cdfid, "vIc", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex index");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "eSc", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge sign");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  /* Max nbt */
  ct = window[1]->nbpt+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpt+1 > ct ? window[n]->nbpt+1 : ct);
  ncw_def_dim(name, cdfid, "nbpt_max", ct, &dims[1]);

  ncw_def_var2(name, cdfid, "bpt",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "bin",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "bin2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Two interior map to ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "dbpt", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Direction of the ghost cell");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  
  /* EDGES                                                           */ 
  /* window->nee[sze]                                                */
  /* window->wse[sze]                                                */
  /* window->e2c[sze][2]                                             */
  /* window->e2v[sze][2]                                             */
  /* window->e2e[sze][2]                                             */
  /* window->eSe[neem][sze]                                          */
  /* window->wAe[neem][sze]                                          */
  /* window->ep[sze]                                                 */
  /* window->em[sze]                                                 */
  /* window->zm1e[sze]                                               */
  /* window->zp1e[sze]                                               */
  /* window->e2k[sze]                                                */
  /* window->m2de[sze]                                               */ 
  /* window->bot_e1[szeS]                                            */
  /* window->sur_e1[szeS]                                            */
  /* window->bpte1[nbpte1]                                           */
  /* window->bine1[nbpte1]                                           */
  /* window->w2_e1[szeS]                                             */
  /* window->w3_e1[sze]                                              */

  ct = window[1]->sze;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->sze > ct ? window[n]->sze : ct);
  ncw_def_dim(name, cdfid, "sze_max", ct, &szew);
  dims[1] = szew;

  ncw_def_var2(name, cdfid, "wse",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "nee",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edges in tangential velocity map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "m2de", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2k", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertical edge number");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2ijk", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to (i,j,k) map for structured grids");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "ep", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map (positive)");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "em", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map (negative)");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "zp1e", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "zm1e", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  dims[1] = szew;
  dims[2] = tid;
  ncw_def_var2(name, cdfid, "e2c", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2v", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2e", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = window[1]->neem+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->neem+1 > ct ? window[n]->neem+1 : ct);
  ncw_def_dim(name, cdfid, "neem_max", ct, &dims[1]);

  dims[2] = szew;
  ncw_def_var2(name, cdfid, "eSe", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Tangential velocity edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "wAe", NC_DOUBLE, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Tangential velocity weights");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = window[1]->szeS;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szeS > ct ? window[n]->szeS : ct);
  ncw_def_dim(name, cdfid, "szeS_max", ct, &szeSw);
  dims[1] = szeSw;
  ncw_def_var2(name, cdfid, "bot_e1",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Bottom edge coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");
  ncw_def_var2(name, cdfid, "sur_e1",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Surface edge coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ct = window[1]->nbpte1+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte1+1 > ct ? window[n]->nbpte1+1 : ct);
  ncw_def_dim(name, cdfid, "nbpte1_max",  ct,  &dims[1]);
  ncw_def_var2(name, cdfid, "bpte1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D ghost edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "bine1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to 3D ghost edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = window[1]->nbpte1S+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte1S+1 > ct ? window[n]->nbpte1S+1 : ct);
  ncw_def_dim(name, cdfid, "nbpte1S_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bpte1S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D ghost edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ncw_def_var2(name, cdfid, "bine1S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to 2D ghost edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  /* VERTICES                                                        */ 
  /* window->nve[szvS]                                               */
  /* window->nvc[szvS]                                               */
  /* window->wsv[szv]                                                */
  /* window->v2c[szv][nvcm]                                          */
  /* window->v2e[szv][nvem]                                          */
  /* window->eSv[nvem][szv]                                          */
  /* window->zm1v[szv]                                               */
  /* window->zp1v[szv]                                               */
  /* window->dualarea[szvS]                                          */
  /* window->dualareap[szvS][nvcm]                                   */
  /* window->m2dv[szv]                                               */
  /* window->w2_e2[szvS]                                             */
  /* window->w3_e2[szv]                                              */
  ct = window[1]->szv;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szv > ct ? window[n]->szv : ct);
  ncw_def_dim(name, cdfid, "szv_max", ct, &szvw);
  dims[1] = szvw;

  ncw_def_var2(name, cdfid, "wsv",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "m2dv", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "zp1v", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "zm1v", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "v2ijk", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex to (i,j,k) map for structured grids");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ct = window[1]->nvem + 1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nvem+1 > ct ? window[n]->nvem+1 : ct);
  ncw_def_dim(name, cdfid, "nve_max", ct, &nvem);
  dims[1] = szvw;
  dims[2] = nvem;
  ncw_def_var2(name, cdfid, "v2e", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ct = window[1]->nvcm + 1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nvcm+1 > ct ? window[n]->nvcm+1 : ct);
  ncw_def_dim(name, cdfid, "nvc_max", ct, &nvcm);
  dims[1] = szvw;
  dims[2] = nvcm;
  ncw_def_var2(name, cdfid, "v2c", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex to centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  dims[1] = nvem;
  dims[2] = szvw;
  ncw_def_var2(name, cdfid, "eSv", NC_INT, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex sign");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ct = window[1]->szvS;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szvS > ct ? window[n]->szvS : ct);
  ncw_def_dim(name, cdfid, "szvS_max", ct, &szvSw);
  dims[1] = szvSw;

  ncw_def_var2(name, cdfid, "nve",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number edges from a vertex");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  ncw_def_var2(name, cdfid, "nvc",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number centres around a vertex");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  ncw_def_var2(name, cdfid, "dualarea",  NC_DOUBLE, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Area of the dual");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  dims[1] = szvSw;
  dims[2] = nvcm;
  ncw_def_var2(name, cdfid, "dualareap",  NC_DOUBLE, 3, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Partial area of the dual");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");
  
  ct = window[1]->n2_t+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_t+1 > ct ? window[n]->n2_t+1 : ct);
  ncw_def_dim(name, cdfid, "n2_t_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  
  ct = window[1]->n2_e1+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_e1+1 > ct ? window[n]->n2_e1+1 : ct);
  ncw_def_dim(name, cdfid, "n2_e1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_e1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ct = window[1]->n2_e2+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_e2+1 > ct ? window[n]->n2_e2+1 : ct);
  ncw_def_dim(name, cdfid, "n2_e2_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_e2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for vertices");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");
  
  ct = window[1]->n3_t+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_t+1 > ct ? window[n]->n3_t+1 : ct);
  ncw_def_dim(name, cdfid, "n3_t_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ct = window[1]->n3_e1+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_e1+1 > ct ? window[n]->n3_e1+1 : ct);
  ncw_def_dim(name, cdfid, "n3_e1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_e1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = window[1]->n3_e2+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_e2+1 > ct ? window[n]->n3_e2+1 : ct);
  ncw_def_dim(name, cdfid, "n3_e2_max", geom->n3_e2, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_e2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for vertices");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  dims[1] = szcid;
  ncw_def_var2(name, cdfid, "s2i", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian i map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  ncw_def_var2(name, cdfid, "s2j", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian j map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  ncw_def_var2(name, cdfid, "s2k", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian k map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ct = window[1]->nm2s+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nm2s+1 > ct ? window[n]->nm2s+1 : ct);
  ncw_def_dim(name, cdfid, "nm2s_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "m2s",   NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Master to slave map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ct = window[1]->nm2se1+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nm2se1+1 > ct ? window[n]->nm2se1+1 : ct);
  ncw_def_dim(name, cdfid, "nm2se1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "m2se1",   NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Master to slave map for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");
  
  ct = window[1]->ns2m+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->ns2m+1 > ct ? window[n]->ns2m+1 : ct);
  ncw_def_dim(name, cdfid, "ns2m_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "s2m",   NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Slave to master map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ct = window[1]->ns2me1+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->ns2me1+1 > ct ? window[n]->ns2me1+1 : ct);
  ncw_def_dim(name, cdfid, "ns2me1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "s2me1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Slave to master map for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  dims[1] = szcSid;
  ncw_def_var(name, cdfid, "botz", NC_DOUBLE, 2, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Bottom depth");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  // } nwindows

  /* Global attributes */
  write_text_att(cdfid, NC_GLOBAL, "Parameter_filename", iname);
  nc_put_att_int(cdfid, NC_GLOBAL, "Windows", NC_INT, 1, &window[1]->nwindows);
  write_text_att(cdfid, NC_GLOBAL, "EMS Version", version);
  write_text_att(cdfid, NC_GLOBAL, "Executable", executable);
  getcwd(key, MAXSTRLEN);
  write_text_att(cdfid, NC_GLOBAL, "Directory", key);
  write_text_att(cdfid, NC_GLOBAL, "Parameterheader", params->parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "Grid_description", params->grid_desc);
  time(&t);
  write_text_att(cdfid, NC_GLOBAL, "Date_created", ctime(&t));

  ncw_enddef(name, cdfid);

  /* 
   * Write the data 
   */
  d1 = i_alloc_1d(geom->sgsiz);
  count[0] = dumpdata->nz + 1;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "z_grid"), start, count,
                     dumpdata->gridz);
  count[0] = geom->szc;
  for (c = 1; c < geom->szc; c++) d1[c] = geom->fm[c].wn;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wn"), start, count, d1);
  for (c = 1; c < geom->szc; c++) d1[c] = geom->fm[c].sc;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sc"), start, count, d1);
  for (c = 1; c < geom->szc; c++) d1[c] = geom->fm[c].ac;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ac"), start, count, d1);
  i_free_1d(d1);

  nc_put_var_int(cdfid, ncw_var_id(cdfid, "nz"), &geom->nz);
  nc_put_var_int(cdfid, ncw_var_id(cdfid, "sednz"), &geom->sednz);

  for (n = 1; n <= window[1]->nwindows; n++) {
    /* 1D arrays */
    start[0] = n;
    count[0] = 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "enon_w"), start, count, &window[n]->enon);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "enonS_w"), start, count, &window[n]->enonS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ewet_w"), start, count, &window[n]->ewet);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ewetS_w"), start, count, &window[n]->ewetS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "snon_w"), start, count, &window[n]->snon);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "snonS_w"), start, count, &window[n]->snonS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "szc_w"), start, count, &window[n]->szc);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "szcS_w"), start, count, &window[n]->szcS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sze_w"), start, count, &window[n]->sze);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "szeS_w"), start, count, &window[n]->szeS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "szv_w"), start, count, &window[n]->szv);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "szvS_w"), start, count, &window[n]->szvS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "npem_w"), start, count, &window[n]->npem);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nvem_w"), start, count, &window[n]->nvem);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nvcm_w"), start, count, &window[n]->nvcm);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "neem_w"), start, count, &window[n]->neem);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpt_w"), start, count, &window[n]->nbpt);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte1_w"), start, count, &window[n]->nbpte1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte1S_w"), start, count, &window[n]->nbpte1S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_t_w"), start, count, &window[n]->n2_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_t_w"), start, count, &window[n]->n3_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_e1_w"), start, count, &window[n]->n2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_e1_w"), start, count, &window[n]->n3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_e2_w"), start, count, &window[n]->n2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_e2_w"), start, count, &window[n]->n3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2s_w"), start, count, &window[n]->nm2s);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2m_w"), start, count, &window[n]->ns2m);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2se1_w"), start, count, &window[n]->nm2se1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2me1_w"), start, count, &window[n]->ns2me1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbptS_w"), start, count, &window[n]->nbptS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbe1_w"), start, count, &window[n]->nbe1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbe1S_w"), start, count, &window[n]->nbe1S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbe2_w"), start, count, &window[n]->nbe2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbe2S_w"), start, count, &window[n]->nbe2S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2_t_w"), start, count, &window[n]->v2_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b2_t_w"), start, count, &window[n]->b2_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a2_t_w"), start, count, &window[n]->a2_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v3_t_w"), start, count, &window[n]->v3_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b3_t_w"), start, count, &window[n]->b3_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a3_t_w"), start, count, &window[n]->a3_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2_e1_w"), start, count, &window[n]->v2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b2_e1_w"), start, count, &window[n]->b2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a2_e1_w"), start, count, &window[n]->a2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "x2_e1_w"), start, count, &window[n]->x2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2_e2_w"), start, count, &window[n]->v2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b2_e2_w"), start, count, &window[n]->b2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a2_e2_w"), start, count, &window[n]->a2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v3_e1_w"), start, count, &window[n]->v3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b3_e1_w"), start, count, &window[n]->b3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a3_e1_w"), start, count, &window[n]->a3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "x3_e1_w"), start, count, &window[n]->x3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v3_e2_w"), start, count, &window[n]->v3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b3_e2_w"), start, count, &window[n]->b3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a3_e2_w"), start, count, &window[n]->a3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2mS_w"), start, count, &window[n]->ns2mS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2sS_w"), start, count, &window[n]->nm2sS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2se1S_w"), start, count, &window[n]->nm2se1S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2me1S_w"), start, count, &window[n]->ns2me1S);

    /* 2D arrays */
    count[0] = 1; /* one window at a time */
    count[1] = window[n]->a2_t + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_t"), start, count, window[n]->bot_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_t"), start, count, window[n]->sur_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nsur_t"), start, count, window[n]->nsur_t);

    count[1] = window[n]->szeS;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_e1"), start, count, window[n]->bot_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_e1"), start, count, window[n]->sur_e1);

    count[1] = window[n]->szc;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2d"), start, count, window[n]->m2d);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wsa"), start, count, window[n]->wsa);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1"), start, count, window[n]->zp1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1"), start, count, window[n]->zm1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wgst"), start, count, window[n]->wgst);

    count[1] = window[n]->szcS;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "npe"), start, count, window[n]->npe);
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "cellarea"), start, count, window[n]->cellarea);

    count[1] = window[n]->npem+1;
    count[2] = window[n]->szc;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2c"), start, count, window[n]->c2c[0]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2e"), start, count, window[n]->c2e[0]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2v"), start, count, window[n]->c2v[0]);

    count[1] = window[n]->sze;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2de"), start, count, window[n]->m2de);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wse"), start, count, window[n]->wse);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1e"), start, count, window[n]->zp1e);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1e"), start, count, window[n]->zm1e);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ep"), start, count, window[n]->ep);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "em"), start, count, window[n]->em);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2k"), start, count, window[n]->e2k);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2ijk"), start, count, window[n]->e2ijk);

    count[1] = window[n]->sze;
    count[2] = 2;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2c"), start, count, window[n]->e2c[0]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2e"), start, count, window[n]->e2e[0]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2v"), start, count, window[n]->e2v[0]);

    count[1] = window[n]->neem+1;
    count[2] = window[n]->sze;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSe"), start, count, window[n]->eSe[0]);
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "wAe"), start, count, window[n]->wAe[0]);

    count[1] = window[n]->szeS;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nee"), start, count, window[n]->nee);

    count[1] = window[n]->szv;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2dv"), start, count, window[n]->m2dv);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wsv"), start, count, window[n]->wsv);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1v"), start, count, window[n]->zp1v);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1v"), start, count, window[n]->zm1v);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2ijk"), start, count, window[n]->v2ijk);

    count[1] = window[n]->szv;
    count[2] = window[n]->nvcm+1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2c"), start, count, window[n]->v2c[0]);

    count[1] = window[n]->szv;
    count[2] = window[n]->nvem + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2e"), start, count, window[n]->v2e[0]);

    count[1] = window[n]->npem + 1;
    count[2] = window[n]->szcS;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "vIc"), start, count, window[n]->vIc[0]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSc"), start, count, window[n]->eSc[0]);

    count[1] = window[n]->szvS;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nve"), start, count, window[n]->nve);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nvc"), start, count, window[n]->nvc);

    count[1] = window[n]->nvem+1;
    count[2] = window[n]->szv;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSv"), start, count, window[n]->eSv[0]);

    count[1] = window[n]->szvS;
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "dualarea"), start, count, window[n]->dualarea);

    count[1] = window[n]->szvS;
    count[2] = window[n]->nvcm+1;
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "dualareap"), start, count,
			window[n]->dualareap[0]);
    
    count[1] = window[n]->nbpt + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpt"), start, count, window[n]->bpt);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bin"), start, count, window[n]->bin);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bin2"), start, count, window[n]->bin2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "dbpt"), start, count, window[n]->dbpt);

    count[1] = window[n]->nbpte1 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1"), start, count, window[n]->bpte1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1"), start, count, window[n]->bine1);

    count[1] = window[n]->nbpte1S + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1S"), start, count, window[n]->bpte1S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1S"), start, count, window[n]->bine1S);

    count[1] = window[n]->n2_t + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_t"), start, count, window[n]->w2_t);
    count[1] = window[n]->n3_t + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_t"), start, count, window[n]->w3_t);
    count[1] = window[n]->n2_e1 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e1"), start, count, window[n]->w2_e1);
    count[1] = window[n]->n3_e1 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e1"), start, count, window[n]->w3_e1);
    count[1] = window[n]->n2_e2 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e2"), start, count, window[n]->w2_e2);
    count[1] = window[n]->n3_e2 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e2"), start, count, window[n]->w3_e2);

    count[1] = window[n]->szc;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2i"), start, count, window[n]->s2i);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2j"), start, count, window[n]->s2j);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2k"), start, count, window[n]->s2k);

    count[1] = window[n]->nm2s + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2s"), start, count, window[n]->m2s);
    count[1] = window[n]->nm2se1 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2se1"), start, count, window[n]->m2se1);

    count[1] = window[n]->ns2m + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2m"), start, count, window[n]->s2m);
    count[1] = window[n]->ns2me1 + 1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2me1"), start, count, window[n]->s2me1);
    count[1] = window[n]->szcS;
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "botz"), start, count, window[n]->botz);
  }

  ncw_close(name, cdfid);
  /* check_window_map_us(window, name); */
}

/* END dump_windows_us()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads window geometry from file                                   */
/* Modified Windows based format which avoids excessively large      */
/* large numbers of dimension for number of windows > 38 and hence   */
/* does not run into the NC_MAX_DIMS limit of 1024.                  */
/*-------------------------------------------------------------------*/
void read_windows_us(geometry_t *geom, geometry_t **window, char *name)
{
  int fid, n;
  int c, cc, cg, *d1;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  size_t d;
  char key[MAXSTRLEN];
  int dims[10];
  int vid, oid, tid;
  long t;
  double *layers;
  int nz, sednz;
  int nwin;

  /* Open the dump file for reading                                  */
  ncw_open(name, NC_NOWRITE, &fid);

  d1 = i_alloc_1d(geom->szc);
  count[0] = geom->szc;
  sprintf(key, "wn");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c < geom->szc; c++) geom->fm[c].wn = d1[c];
  sprintf(key, "sc");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c < geom->szc; c++) geom->fm[c].sc = d1[c];
  sprintf(key, "ac");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c < geom->szc; c++) geom->fm[c].ac = d1[c];
  i_free_1d(d1);

  /* Get layers info and check                                       */
  layers = d_alloc_1d(geom->nz + 1);
  count[0] = geom->nz + 1;;
  sprintf(key, "z_grid");
  nc_get_vara_double(fid, ncw_var_id(fid, key), start, count, layers);
  for (cc = 0; cc < geom->nz; cc++) {
    if (layers[cc] != geom->layers[cc]) {
      hd_warn("Layer %d : layer depths = %f (.prm) %f (file)\n", cc,
	     geom->layers[cc], layers[cc]);
      hd_quit("Incompatible layer depths in window map file %s with parameter file.\n", name);
    }
  }
  d_free_1d(layers);

  /* Get global variables, i.e. ones that are the same across all    */
  /* windows.                                                        */
  ncw_inq_dimlen(name, fid, ncw_dim_id(fid, "nwp1"), &d);
  nwin = d-1; // nwp1 is numnber of windows plus 1
  if (nwin != geom->nwindows)
    hd_quit("Incompatible number of windows in file %s %d vs %d\n", name, nwin, geom->nwindows);

  oid = ncw_dim_id(fid, "one");
  tid = ncw_dim_id(fid, "two");
  nc_get_var_int(fid, ncw_var_id(fid, "nz"), &nz);
  nc_get_var_int(fid, ncw_var_id(fid, "sednz"), &sednz);


  /* Get dimensions */
  for (n = 1; n <= geom->nwindows; n++) {
    double *botz;

    emstag(LINFO,"read_windows_wb","Reading window %d map\n", n);

    window[n] = window_alloc();
    window[n]->wn = n;

    window[n]->nwindows = nwin;
    window[n]->nz    = nz;
    window[n]->sednz = sednz;

    start[0] = n;
    count[0] = 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "enon_w"), start, count, &window[n]->enon);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ewet_w"), start, count, &window[n]->ewet);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "snon_w"), start, count, &window[n]->snon);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "enonS_w"), start, count, &window[n]->enonS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ewetS_w"), start, count, &window[n]->ewetS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "snonS_w"), start, count, &window[n]->snonS);

    ncw_get_vara_int(name, fid, ncw_var_id(fid, "szc_w"), start, count, &window[n]->szc);
    window[n]->sgnum = window[n]->szc - 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "szcS_w"), start, count, &window[n]->szcS);
    window[n]->sgnumS = window[n]->szcS - 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "sze_w"), start, count, &window[n]->sze);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "szeS_w"), start, count, &window[n]->szeS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "szv_w"), start, count, &window[n]->szv);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "szvS_w"), start, count, &window[n]->szvS);

    ncw_get_vara_int(name, fid, ncw_var_id(fid, "npem_w"), start, count, &window[n]->npem);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nvem_w"), start, count, &window[n]->nvem);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nvcm_w"), start, count, &window[n]->nvcm);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "neem_w"), start, count, &window[n]->neem);

    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_t_w"), start, count, &window[n]->n2_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_t_w"), start, count, &window[n]->n3_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a2_t_w"), start, count, &window[n]->a2_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a3_t_w"), start, count, &window[n]->a3_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b2_t_w"), start, count, &window[n]->b2_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b3_t_w"), start, count, &window[n]->b3_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2_t_w"), start, count, &window[n]->v2_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v3_t_w"), start, count, &window[n]->v3_t);

    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_e1_w"), start, count, &window[n]->n2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_e1_w"), start, count, &window[n]->n3_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "x2_e1_w"), start, count, &window[n]->x2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "x3_e1_w"), start, count, &window[n]->x3_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a2_e1_w"), start, count, &window[n]->a2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a3_e1_w"), start, count, &window[n]->a3_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b2_e1_w"), start, count, &window[n]->b2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b3_e1_w"), start, count, &window[n]->b3_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2_e1_w"), start, count, &window[n]->v2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v3_e1_w"), start, count, &window[n]->v3_e1);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_e2_w"), start, count, &window[n]->n2_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_e2_w"), start, count, &window[n]->n3_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a2_e2_w"), start, count, &window[n]->a2_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "a3_e2_w"), start, count, &window[n]->a3_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b2_e2_w"), start, count, &window[n]->b2_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "b3_e2_w"), start, count, &window[n]->b3_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2_e2_w"), start, count, &window[n]->v2_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v3_e2_w"), start, count, &window[n]->v3_e2);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpt_w"), start, count, &window[n]->nbpt);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbptS_w"), start, count, &window[n]->nbptS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte1_w"), start, count, &window[n]->nbpte1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte1S_w"), start, count, &window[n]->nbpte1S);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbe1_w"), start, count, &window[n]->nbe1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbe1S_w"), start, count, &window[n]->nbe1S);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbe2_w"), start, count, &window[n]->nbe2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbe2S_w"), start, count, &window[n]->nbe2S);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nm2s_w"), start, count, &window[n]->nm2s);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ns2m_w"), start, count, &window[n]->ns2m);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nm2se1_w"), start, count, &window[n]->nm2se1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ns2me1_w"), start, count, &window[n]->ns2me1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nm2sS_w"), start, count, &window[n]->nm2sS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ns2mS_w"), start, count, &window[n]->ns2mS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nm2se1S_w"), start, count, &window[n]->nm2se1S);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ns2me1S_w"), start, count, &window[n]->ns2me1S);

    /* CENTRES                                                       */
    /* window->npe[szcS]                                             */
    /* window->c2c[npem][szc]                                        */
    /* window->c2e[npem][szc]                                        */
    /* window->c2v[npem][szc]                                        */
    /* window->vIc[npem][szcS]                                       */
    /* window->eSc[npem][szcS]                                       */
    /* window->zm1[szc]                                              */
    /* window->zp1[szc]                                              */
    /* window->m2d[szc]                                              */
    /* window->wsa[szc]                                              */
    /* window->wgst[szc]                                             */
    /* window->bot_t[szcS]                                           */
    /* window->sur_t[szcS]                                           */
    /* window->nsur_t[szcS]                                          */
    /* window->bpt[nbpt]                                             */
    /* window->bin[nbpt]                                             */
    /* window->bin2[nbpt]                                            */
    /* window->w2_t[szcS]                                            */
    /* window->w3_t[szc]                                             */
    alloc_geom_us(window[n], CENTRE_A);
    window[n]->w2_t = i_alloc_1d(window[n]->szcS);
    window[n]->w3_t = i_alloc_1d(window[n]->szc);
    window[n]->sur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->nsur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->bot_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->s2i = i_alloc_1d(window[n]->szc);
    window[n]->s2j = i_alloc_1d(window[n]->szc);
    window[n]->s2k = i_alloc_1d(window[n]->szc);
    window[n]->bpt = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin2 = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->dbpt = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->wgst = i_alloc_1d(window[n]->enon + 1);
    window[n]->m2d = i_alloc_1d(window[n]->szc);
    window[n]->c2e = i_alloc_2d(window[n]->szc, window[n]->npem+1);
    window[n]->c2v = i_alloc_2d(window[n]->szc, window[n]->npem+1);
    window[n]->vIc = i_alloc_2d(window[n]->szcS, window[n]->npem + 1);
    window[n]->eSc = i_alloc_2d(window[n]->szcS, window[n]->npem + 1);
    window[n]->m2s = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->s2m = i_alloc_1d(window[n]->ns2m + 1);

    count[1] = window[n]->a2_t + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bot_t"), start, count, window[n]->bot_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "sur_t"), start, count, window[n]->sur_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nsur_t"), start, count, window[n]->nsur_t);

    count[1] = window[n]->szc;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2d"), start, count, window[n]->m2d);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "wsa"), start, count, window[n]->wsa);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1"), start, count, window[n]->zp1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1"), start, count, window[n]->zm1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "wgst"), start, count, window[n]->wgst);

    count[1] = window[n]->szcS;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "npe"), start, count, window[n]->npe);

    count[1] = window[n]->npem + 1;
    count[2] = window[n]->szc;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2c"), start, count, window[n]->c2c[0]);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2e"), start, count, window[n]->c2e[0]);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2v"), start, count, window[n]->c2v[0]);

    count[1] = window[n]->nbpt + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpt"), start, count, window[n]->bpt);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bin"), start, count, window[n]->bin);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bin2"), start, count, window[n]->bin2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "dbpt"), start, count, window[n]->dbpt);

    count[1] = window[n]->n2_t + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_t"), start, count, window[n]->w2_t);
    count[1] = window[n]->n3_t + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_t"), start, count, window[n]->w3_t);
 
    count[1] = window[n]->szc;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2i"), start, count, window[n]->s2i);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2j"), start, count, window[n]->s2j);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2k"), start, count, window[n]->s2k);

    count[1] = window[n]->npem + 1;
    count[2] = window[n]->szcS;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "vIc"), start, count, window[n]->vIc[0]);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSc"), start, count, window[n]->eSc[0]);

    count[1] = window[n]->nm2s + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2s"), start, count, window[n]->m2s);

    count[1] = window[n]->ns2m + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2m"), start, count, window[n]->s2m);

    /* EDGES                                                         */ 
    /* window->nee[sze]                                              */
    /* window->wse[sze]                                              */
    /* window->e2c[sze][2]                                           */
    /* window->e2v[sze][2]                                           */
    /* window->e2e[sze][2]                                           */
    /* window->eSe[neem][sze]                                        */
    /* window->wAe[neem][sze]                                        */
    /* window->ep[sze]                                               */
    /* window->em[sze]                                               */
    /* window->zm1e[sze]                                             */
    /* window->zp1e[sze]                                             */
    /* window->e2k[sze]                                              */
    /* window->m2de[sze]                                             */ 
    /* window->bot_e1[szeS]                                          */
    /* window->sur_e1[szeS]                                          */
    /* window->bpte1[nbpte1]                                         */
    /* window->bine1[nbpte1]                                         */
    /* window->w2_e1[szeS]                                           */
    /* window->w3_e1[sze]                                            */
    alloc_geom_us(window[n], EDGE_A);
    window[n]->e2c = i_alloc_2d(2, window[n]->sze);
    window[n]->e2v = i_alloc_2d(2, window[n]->sze);
    window[n]->e2e = i_alloc_2d(2, window[n]->sze);
    window[n]->eSe = i_alloc_2d(window[n]->sze, window[n]->neem + 1);
    window[n]->wAe = d_alloc_2d(window[n]->sze, window[n]->neem + 1);
    window[n]->ep = i_alloc_1d(window[n]->sze);
    window[n]->em = i_alloc_1d(window[n]->sze);
    window[n]->zm1e = i_alloc_1d(window[n]->sze);
    window[n]->zp1e = i_alloc_1d(window[n]->sze);
    window[n]->e2k = i_alloc_1d(window[n]->sze);
    window[n]->e2ijk = i_alloc_1d(window[n]->sze);
    window[n]->m2de = i_alloc_1d(window[n]->sze);
    window[n]->sur_e1 = i_alloc_1d(window[n]->szeS);
    window[n]->bot_e1 = i_alloc_1d(window[n]->szeS);
    window[n]->bpte1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bine1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bpte1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->bine1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->w2_e1 = i_alloc_1d(window[n]->szeS);
    window[n]->w3_e1 = i_alloc_1d(window[n]->sze);
    window[n]->m2se1 = i_alloc_1d(window[n]->nm2se1 + 1);
    window[n]->s2me1 = i_alloc_1d(window[n]->ns2me1 + 1);

    count[1] = window[n]->szeS;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bot_e1"), start, count, window[n]->bot_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "sur_e1"), start, count, window[n]->sur_e1);
   
    count[1] = window[n]->sze;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2de"), start, count, window[n]->m2de);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "wse"), start, count, window[n]->wse);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1e"), start, count, window[n]->zp1e);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1e"), start, count, window[n]->zm1e);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ep"), start, count, window[n]->ep);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "em"), start, count, window[n]->em);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2k"), start, count, window[n]->e2k);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2ijk"), start, count, window[n]->e2ijk);

    count[1] = window[n]->sze;
    count[2] = 2;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2c"), start, count, window[n]->e2c[0]);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2e"), start, count, window[n]->e2e[0]);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2v"), start, count, window[n]->e2v[0]);

    count[1] = window[n]->neem + 1;
    count[2] = window[n]->sze;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSe"), start, count, window[n]->eSe[0]);
    ncw_get_vara_double(name, fid, ncw_var_id(fid, "wAe"), start, count, window[n]->wAe[0]);

    count[1] = window[n]->szeS;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nee"), start, count, window[n]->nee);

    count[1] = window[n]->nbpte1 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpte1"), start, count, window[n]->bpte1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bine1"), start, count, window[n]->bine1);

    count[1] = window[n]->nbpte1S + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpte1S"), start, count, window[n]->bpte1S);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "bine1S"), start, count, window[n]->bine1S);

    count[1] = window[n]->n2_e1 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_e1"), start, count, window[n]->w2_e1);
    count[1] = window[n]->n3_e1 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_e1"), start, count, window[n]->w3_e1);

    count[1] = window[n]->nm2se1 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2se1"), start, count, window[n]->m2se1);
    count[1] = window[n]->ns2me1 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2me1"), start, count, window[n]->s2me1);

    /* VERTICES                                                       */ 
    /* window->nve[szvS]                                              */
    /* window->nvc[szvS]                                              */
    /* window->wsv[szv]                                               */
    /* window->v2c[szv][nvcm]                                         */
    /* window->v2e[szv][nvem]                                         */
    /* window->eSv[nvem][szv]                                         */
    /* window->zm1v[szv]                                              */
    /* window->zp1v[szv]                                              */
    /* window->dualarea[szvS]                                         */
    /* window->dualareap[szvS][nvcm]                                  */
    /* window->m2dv[szv]                                              */
    /* window->bot_e2[szvS]                                           */
    /* window->sur_e2[szvS]                                           */
    /* window->w2_e2[szvS]                                            */
    /* window->w3_e2[szv]                                             */
    alloc_geom_us(window[n], VERTEX_A);
    window[n]->w2_e2 = i_alloc_1d(window[n]->szvS);
    window[n]->w3_e2 = i_alloc_1d(window[n]->szv);
    window[n]->m2dv = i_alloc_1d(window[n]->szv);
    window[n]->zm1v = i_alloc_1d(window[n]->szv);
    window[n]->zp1v = i_alloc_1d(window[n]->szv);
    window[n]->v2ijk = i_alloc_1d(window[n]->szv);
    window[n]->v2c = i_alloc_2d(window[n]->nvcm+1, window[n]->szv);
    window[n]->v2e = i_alloc_2d(window[n]->nvem+1, window[n]->szv);
    window[n]->eSv = i_alloc_2d(window[n]->szv, window[n]->nvem + 1);

    count[1] = window[n]->szv;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2dv"), start, count, window[n]->m2dv);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "wsv"), start, count, window[n]->wsv);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1v"), start, count, window[n]->zp1v);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1v"), start, count, window[n]->zm1v);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2ijk"), start, count, window[n]->v2ijk);

    count[1] = window[n]->szv;
    count[2] = window[n]->nvcm + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2c"), start, count, window[n]->v2c[0]);
    count[2] = window[n]->nvem + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2e"), start, count, window[n]->v2e[0]);

    count[1] = window[n]->szvS;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nve"), start, count, window[n]->nve);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nvc"), start, count, window[n]->nvc);

    count[1] = window[n]->nvem + 1;
    count[2] = window[n]->szv;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSv"), start, count, window[n]->eSv[0]);

    count[1] = window[n]->szvS;
    ncw_get_vara_double(name, fid, ncw_var_id(fid, "dualarea"), start, count, window[n]->dualarea);

    count[1] = window[n]->szvS;
    count[2] = window[n]->nvcm + 1;
    ncw_get_vara_double(name, fid, ncw_var_id(fid, "dualareap"), start, count,
			    window[n]->dualareap[0]);

    count[1] = window[n]->n2_e2 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_e2"), start, count, window[n]->w2_e2);
    count[1] = window[n]->n3_e2 + 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_e2"), start, count, window[n]->w3_e2);

    botz = d_alloc_1d(window[n]->szcS);
    count[1] = window[n]->szcS;
    sprintf(key, "botz");
    nc_get_vara_double(fid, ncw_var_id(fid, key), start, count, botz);
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      cg = window[n]->wsa[c];
      if (botz[c] != geom->botz[cg]) {
	hd_warn("Window%d (%d %d)=%5.2f : Global (%d %d)=%5.2f\n", n,
		window[n]->s2i[c], window[n]->s2j[c],  botz[c], 
		geom->s2i[cg], geom->s2j[cg], geom->botz[cg]);
	hd_quit("Incompatible bottom depths in window map file %s with parameter file.\n", name);
      }
    }
    d_free_1d(botz);
    
    /*-----------------------------------------------------------------*/
    /* Make the global to local map for this window. This differs from */
    /* the global to local map in geom in that it is defined only over */
    /* all global cells in the window, including ghost cells which are */
    /* associated with zero window and local coordinate in geom->fm.  */
    window[n]->fm =
      (global_map_t *)malloc(sizeof(global_map_t) * (geom->enon + 1));
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      window[n]->fm[c].wn = window[n]->wn;
      window[n]->fm[c].sc = cc;
      if (geom->fm[c].wn == 0)
	window[n]->fm[c].ac = 0;
      else if (geom->fm[c].wn == window[n]->wn)
	window[n]->fm[c].ac = 1;
      else
	window[n]->fm[c].ac = 2;
    }

    // Define max sizes
    window[n]->szm = max(max(window[n]->szc, window[n]->sze), window[n]->szv);
    window[n]->szmS = max(max(window[n]->szcS, window[n]->szeS), window[n]->szvS);

    emstag(LINFO,"read_windows_wb","Window %d map OK\n", n);
  }
  ncw_close(name, fid);
}

/* END read_windows_us()                                             */
/*-------------------------------------------------------------------*/

/*************************/
/* Check window map      */
/* NOT fully implemented */
/*************************/
void check_window_map_us(geometry_t **window, char *name)
{
  geometry_t **win;
  int e, ee, n, nn, j, c, c2, cc, v;

  win = (geometry_t **)p_alloc_1d(window[1]->nwindows);
  read_windows_us(geom, win, name);

  for (n = 1; n <= geom->nwindows; n++) {
    printf("Checking window %d\n", n);
    if (window[n]->b2_e1 != win[n]->b2_e1)
      printf("b2_e1 %d:%d\n", window[n]->b2_e1, win[n]->b2_e1);
    for (ee = 1; ee <= window[n]->b2_e1; ee++) {
      if (window[n]->w2_e1[ee] != win[n]->w2_e1[ee])
	printf("w2_e1 %d %d:%d\n", ee, window[n]->w2_e1[ee], win[n]->w2_e1[ee]);
    }
    if (window[n]->b2_t != win[n]->b2_t)
      printf("b2_t %d:%d\n", window[n]->b2_t, win[n]->b2_t);
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      if (window[n]->w2_t[cc] != win[n]->w2_t[cc])
	printf("w2_t %d %d:%d\n", cc, window[n]->w2_t[cc], win[n]->w2_t[cc]);
      c = window[n]->w2_t[cc];
      if (window[n]->npe[c] != win[n]->npe[c]) printf ("npe error\n");
      for (nn = 1; nn <= window[n]->npem; nn++) {
	if (window[n]->c2c[nn][c] != win[n]->c2c[nn][c]) printf("c2c error\n");
	if (window[n]->c2e[nn][c] != win[n]->c2e[nn][c]) printf("c2e error\n");
	if (window[n]->c2v[nn][c] != win[n]->c2v[nn][c]) printf("c2v error\n");
	if (window[n]->eSc[nn][c] != win[n]->eSc[nn][c]) printf("eSc error\n");
	if (window[n]->vIc[nn][c] != win[n]->vIc[nn][c]) printf("vIc error\n");
      }
    }
    for (e = 1; e < window[n]->sze; e++) {
      if (window[n]->e2c[e][0] != win[n]->e2c[e][0]) printf("e2c[0] error\n");
      if (window[n]->e2c[e][1] != win[n]->e2c[e][1]) printf("e2c[1] error\n");
      if (window[n]->e2e[e][0] != win[n]->e2e[e][0]) printf("e2e[0] error\n");
      if (window[n]->e2e[e][1] != win[n]->e2e[e][1]) printf("e2e[1] error\n");
      if (window[n]->e2v[e][0] != win[n]->e2v[e][0]) printf("e2v[0] error\n");
      if (window[n]->e2v[e][1] != win[n]->e2v[e][1]) printf("e2v[1] error\n");
      if (window[n]->ep[e] != win[n]->ep[e]) printf("ep error\n");
      if (window[n]->zm1e[e] != win[n]->zm1e[e]) printf("zm1e error\n");
      if (window[n]->zp1e[e] != win[n]->zp1e[e]) printf("zp1e error\n");
      if (window[n]->m2de[e] != win[n]->m2de[e]) printf("m2de error\n");
      for (nn = 1; nn <= window[n]->neem; nn++) {
	if (window[n]->eSe[nn][e] != win[n]->eSe[nn][e]) printf("eSe error\n");
	if (window[n]->wAe[nn][e] != win[n]->wAe[nn][e]) printf("wAe error\n");
      }
    }

    /* 3D */
    if (window[n]->b3_t != win[n]->b3_t)
      printf("b3_t %d:%d\n", window[n]->b3_t, win[n]->b3_t);
    for (cc = 1; cc <= window[n]->b3_t; cc++) {
      if (window[n]->w3_t[cc] != win[n]->w3_t[cc])
	printf("w3_t %d %d:%d\n", cc, window[n]->w3_t[cc], win[n]->w3_t[cc]);
    }

    /* Vertices */
    for (v = 1; v < window[n]->szv; v++) {
      if (window[n]->m2dv[v] != win[n]->m2dv[v]) printf("m2dv error\n");
      if (window[n]->wsv[v]  != win[n]->wsv[v])  printf("wsv error\n");
      if (window[n]->zp1v[v] != win[n]->zp1v[v]) printf("zp1v error\n");
      if (window[n]->zm1v[v] != win[n]->zm1v[v]) printf("zm1v error\n");
    }
  }
  /* unlink(name); */
  hd_quit("done check_window\n");
}

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


void read_mean_atts(master_t *master, int fid)
{
  geometry_t *geom = master->geom;
  parameters_t *params = master->params;
  double d1;
  int cc;

  if (!(params->means & NONE)) {
    nc_get_att_double(fid, NC_GLOBAL, "mean_c", &d1);
    for (cc = 1; cc < geom->szcS; cc++)
      master->meanc[cc] = d1;
    nc_get_att_double(fid, NC_GLOBAL, "mean_next", &master->means_next);
    if (nc_get_att_text(fid, NC_GLOBAL, "mean_mc", params->means_mc) != NC_NOERR)
      sprintf(params->means_mc, "%c", '\0');
  }
}
