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
 *  $Id: df_ugrid.c 7336 2023-04-11 02:35:01Z her127 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define Filli 999999
#define Filld 9.9692099683868690E36;
#define MAXCONSTIT 30

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
void pack_obc_t(int mapsize, double *var, double *pack);
int pack_obc3_t(open_bdrys_t *open, double *var, double *pack);
int pack_obc2_e1(open_bdrys_t *open, double *var, double *pack);
int pack_obc3_e1(open_bdrys_t *open, double *v1, double *v2, double *pack);

#define UGRID_ALL_VARS "u1av u2av uav vav wind1 wtop eta patm u1 u2 u v w dens dens_0 Kz Vz Cd u1kh topz "

double percentiles1[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
			0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
			0.90, 0.95, 1.00};
#define Ns1 (sizeof(percentiles1) / sizeof(double))

void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
				 nc_type fptype, char *modulo);
void write_dump_attributes_ugrid3(dump_data_t *dumpdata, dump_file_t *df, int cdfid,
				  nc_type fptype, char *modulo);
static void df_ugrid_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_ugrid3_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_obc_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_ugrid_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_ugrid_get_dimids(int cdfid, df_ugrid_var_t var, int *ns, int *zid);
static void df_ugrid_get_dimsizes(dump_file_t *df, df_ugrid_var_t var, size_t *sz, size_t *ns, size_t *nz);
static void df_ugrid3_get_dimids(int cdfid, df_ugrid_var_t var, int *ns);
static void df_ugrid3_get_dimsizes(dump_file_t *df, df_ugrid_var_t var, size_t *ns);
static void write_mean_atts(dump_data_t *dumpdata, int fid);
int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref);
void check_print(geometry_t *window, char *tag, int cc, int c, int v1, int v2, int mode);
void check_print_obc(char *tag, int nb, int cc, int v1, int v2);

void write_dump_attributes_ugrid(dump_data_t *dumpdata, int cdfid,
				 nc_type fptype, char *modulo)
{
  master_t *master = dumpdata->master;
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

  if ((vid = ncw_var_id(cdfid, "air_temp")) >= 0) {
    write_text_att(cdfid, vid, "units", "degrees C");
    write_text_att(cdfid, vid, "long_name", "Air temperature");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "cloud")) >= 0) {
    write_text_att(cdfid, vid, "units", "oktas");
    write_text_att(cdfid, vid, "long_name", "Cloud amount");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "precipitation")) >= 0) {
    write_text_att(cdfid, vid, "units", "mm day-1");
    write_text_att(cdfid, vid, "long_name", "Precipitation");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "wet_bulb")) >= 0) {
    write_text_att(cdfid, vid, "units", "degrees C");
    write_text_att(cdfid, vid, "long_name", "Wet bulb temperature");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "dew_point")) >= 0) {
    write_text_att(cdfid, vid, "units", "degrees C");
    write_text_att(cdfid, vid, "long_name", "Dew point temperature");
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
  write_text_att(cdfid, NC_GLOBAL, "UGRID", "3D layered");
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

void write_dump_attributes_ugrid3(dump_data_t *dumpdata, dump_file_t *df, int cdfid,
				 nc_type fptype, char *modulo)
{
  master_t *master = dumpdata->master;
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

  if (df->flag & DF_OBC) {
    nface2 = df->nface2;
    nface3 = df->nface3;
    nvertex2 = 0;
    nedge2 = df->nedge2;
    nedge3 = df->nedge3;
    nvertex3 = 0;
  }

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

  if (!(df->flag & DF_OBC)) {
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
  }

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

  if (!(df->flag & DF_OBC)) {
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
  }

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

  if (!(df->flag & DF_OBC)) {
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
  }

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

  if ((vid = ncw_var_id(cdfid, "air_temp")) >= 0) {
    write_text_att(cdfid, vid, "units", "degrees C");
    write_text_att(cdfid, vid, "long_name", "Air temperature");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "cloud")) >= 0) {
    write_text_att(cdfid, vid, "units", "oktas");
    write_text_att(cdfid, vid, "long_name", "Cloud amount");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if ((vid = ncw_var_id(cdfid, "precipitation")) >= 0) {
    write_text_att(cdfid, vid, "units", "mm day-1");
    write_text_att(cdfid, vid, "long_name", "Precipitation");
    write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
  }

  if (master->sh_f & WETBULB) {
    if ((vid = ncw_var_id(cdfid, "wet_bulb")) >= 0) {
      write_text_att(cdfid, vid, "units", "degrees C");
      write_text_att(cdfid, vid, "long_name", "Wet bulb temperature");
      write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
    }
  } else if (master->sh_f & DEWPOINT) {
    if ((vid = ncw_var_id(cdfid, "dew_point")) >= 0) {
      write_text_att(cdfid, vid, "units", "degrees C");
      write_text_att(cdfid, vid, "long_name", "Dew point temperature");
      write_text_att(cdfid, vid, "coordinates", "t, Mesh2_face_x, Mesh2_face_y");
    }
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
  write_text_att(cdfid, NC_GLOBAL, "UGRID", "3D");
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

}


void df_ugrid_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_ugrid_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t < (df->tout - df->tinc))
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

  if (df->flag & DF_OBC) {
    nvertex2 = 0;
    nvertex3 = 0;
    nface2 = df->nface2;
    nface3 = df->nface3;
    nedge2 = df->nedge2;
    nedge3 = df->nedge3;
  }

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
  if (!(df->flag & DF_OBC)) {
    set_dims(dumpdata->face_dim, dims, nMaxMesh3_face_nodesid, nMesh3_faceid);
    nc_def_var(cdfid, "Mesh3_face_xbnds", NC_DOUBLE, 2, dims, &vid);
    nc_def_var(cdfid, "Mesh3_face_ybnds", NC_DOUBLE, 2, dims, &vid);
  }
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
  if (!(df->flag & DF_OBC)) {
    dims[0] = nMesh3_faceid;
    nc_def_var(cdfid, "coriolis", fptype, 1, dims, &vid);
  }

  /* Grid index (optional) */
  if (!(df->flag & DF_OBC)) {
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

  write_dump_attributes_ugrid3(dumpdata, df, cdfid, fptype, df->modulo);

  nc_enddef(cdfid);

  if (df->flag & DF_OBC)
    df_obc_writegeom(dumpdata, df);
  else
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
  if (df->flag & DF_OBC) newt -= df->tinc;
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

    if (df->flag & DF_OBC) {
      open_bdrys_t *open = geom->open[df->obcid];
      if (p == NULL) continue;
      df_ugrid3_get_dimsizes(df, vn, &count[1]);
      if (vn.xylocation & CL_EDGE) {
	void *v1 = vn.v;
	void *v2 = vn.vt;
	pack_obc3_e1(open, (*(double **)v1), (*(double **)v2), dumpdata->w1);
	nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count, dumpdata->w1);
      } else {
	pack_obc_t(count[1], (*(double **)p), dumpdata->w1);
	nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count, dumpdata->w1);
      }
    } else {
      memset(dumpdata->w1,  0, geom->szm*sizeof(double));
      df_ugrid3_get_dimsizes(df, vn, &count[1]);
      pack_ugrid3(vn.hmap, count[1], (*(double **)p), dumpdata->w1, oset);
      nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count,
		       dumpdata->w1);
    }
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


static void df_obc_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_ugrid_data_t *data = (df_ugrid_data_t *)df->private_data;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  open_bdrys_t *open = geom->open[df->obcid];
  int fid = data->fid;
  int c, cc, ee, e;
  int m, n;

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

  count[0] = df->nface2;
  pack_ugrid3(open->obc_t, open->no2_t, geom->cellx, dumpdata->w1s, 1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_face_x"), start, count, dumpdata->w1s);
  pack_ugrid3(open->obc_t, open->no2_t, geom->celly, dumpdata->w1s, 1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_face_y"), start, count, dumpdata->w1s);

  count[0] = pack_obc2_e1(open, geom->u1x, dumpdata->w1s);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_edge_x"), start, count, dumpdata->w1s);
  count[0] = pack_obc2_e1(open, geom->u1y, dumpdata->w1s);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_edge_y"), start, count, dumpdata->w1s);

  count[0] = pack_obc2_e1(open, geom->h1au1, dumpdata->w1s);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au1"), start, count, dumpdata->w1s);
  count[0] = pack_obc2_e1(open, geom->h2au1, dumpdata->w1s);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au1"), start, count, dumpdata->w1s);
  count[0] = pack_obc2_e1(open, geom->thetau1, dumpdata->w1s);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau1"), start, count, dumpdata->w1s);

  count[0] = df->nface2;
  pack_ugrid3(open->obc_t, open->no2_t, geom->botz, dumpdata->w1s, 1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "Mesh3_depth"), start, count, dumpdata->w1s);

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
    if (df->flag & DF_OBC) {
      open_bdrys_t *open = geom->open[df->obcid];
      var->v = (void *)&open->transfer_eta;
    }
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

  else if (strcmp(name, "air_temp") == 0) {
    var->v = (void *)&master->airtemp;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "cloud") == 0) {
    var->v = (void *)&master->cloud;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "dew_point") == 0) {
    var->v = (void *)&master->wetb;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "wet_bulb") == 0) {
    var->v = (void *)&master->wetb;
    var->xylocation = CL_SP2|CL_FACE;
  }

  else if (strcmp(name, "precipitation") == 0) {
    var->v = (void *)&master->precip;
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
    if (df->flag & DF_OBC) {
      open_bdrys_t *open = geom->open[df->obcid];
      var->v = (void *)&open->transfer_u1;
      var->vt = (void *)&open->transfer_u2;
    }
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

  else if (strcmp(name, "dz") == 0) {
    var->v = (void *)&master->dz;
    var->xylocation = CL_SP3|CL_FACE;
    var->ndims = 2;
    var->zlocation = CL_CENTRE;
    var->hmap = geom->w3_t;
  }

  else if (strcmp(name, "dzu1") == 0) {
    var->v = (void *)&master->dzu1;
    var->xylocation = CL_SP3|CL_EDGE;
    var->zlocation = CL_CENTRE;
    var->ndims = 2;
    var->hmap = geom->w3_e1;
    var->vmap = geom->e2k;
    var->m2d = geom->m2de;
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
	if (df->flag & DF_OBC) {
	  open_bdrys_t *open = geom->open[df->obcid];
	  int tm = open->trm[n];
	  var->v = (void *)&open->t_transfer[tm];
	}
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
    hd_quit("df_ugrid:df_ugrid_get_dimsizes error in xylocation\n");

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
    hd_quit("df_ugrid:df_ugrid3_get_dimsizes error in xylocation\n");
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

int pack_obc3_t(open_bdrys_t *open, double *var, double *pack)
{
  int n, m, cc, c, ee, e;

  n = 0;
  for (cc = 1; cc <= open->no2_t; cc++) {
    c = open->obc_t[cc];
    pack[n++] = var[c];
  }
  for (ee = 1; ee <= open->no2_e1; ee++) {
    c = open->obc_e2[ee];
    e = open->obc_e1[ee];
    for(m = 1; m <= open->bgz; m++) {
      c = open->omape[ee][c];
      pack[n++] = var[c];
    }
  }
  return(n);
}

void pack_obc_t(int mapsize, double *var, double *pack)
{
  int cc;

  for (cc = 1; cc <= mapsize; cc++) {
    pack[cc-1] = var[cc];
  }
}

/*-------------------------------------------------------------------*/
/* Packs 2D normal and tangential edges into a contiguous array      */
int pack_obc2_e1(open_bdrys_t *open, double *var, double *pack)
{
  int n, m, cc, c, ee, e;

  n = 0;
  for (ee = 1; ee <= open->no2_e1; ee++) {
    e = open->obc_e1[ee];
    pack[n++] = var[e];
  }
  for (ee = open->no3_e1 + 1; ee <= open->to2_e1; ee++) {
    e = open->obc_e1[ee];
    pack[n++] = var[e];
  }
  return(n);
}

/*-------------------------------------------------------------------*/
/* Packs the 3D u1 and u2 transfer vectors into a contiguous array   */
int pack_obc3_e1(open_bdrys_t *open, double *v1, double *v2, double *pack)
{
  int n, ee;

  n = 0;
  for (ee = 1; ee <= open->no3_e1; ee++) {
    pack[n++] = v1[ee];
  }
  for (ee = open->no3_e1 + 1; ee <= open->to3_e1; ee++) {
    pack[n++] = v2[ee];
  }
  return(n);
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
  parameters_t *params = master->params;
  int cdfid, n, nb, c, cc, vc, i;
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
  int kgridid;          /* K dimension id at grid corner */
  int oid, tid;
  int nwinsid;
  int nwins = geom->nwindows;
  char key[MAXSTRLEN];
  char tname[MAXSTRLEN], tnameu[MAXSTRLEN], tnamev[MAXSTRLEN];
  int *d1;
  int ct;
  int ncmode = overwrite_output() ? NC_CLOBBER:NC_NOCLOBBER;
  int npem, nvem, nvcm;
  int szcw, szcSw;
  int szew, szeSw;
  int szvw, szvSw;

  int nobc;
  int no2_t, no2_e1;
  int tf = 0;
  int tmask[MAXCONSTIT+1];

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
  
  ncw_def_var(name, cdfid, "nobc_w", NC_INT, 1, dims, &vid);

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
  ct = window[1]->szc+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szc > ct ? window[n]->szc+1 : ct);
  ncw_def_dim(name, cdfid, "szc_max", ct, &szcw);

  /* Max 2D centres */
  ct = window[1]->szcS+1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->szcS > ct ? window[n]->szcS+1 : ct);
  ncw_def_dim(name, cdfid, "szcS_max", ct, &szcSw);

  /* Max npe */
  ct = window[1]->npem + 1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->npem+1 > ct ? window[n]->npem+1 : ct);
  ncw_def_dim(name, cdfid, "npe_max", ct, &npem);

  /* Max obcs */
  ct = window[1]->nobc;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nobc > ct ? window[n]->nobc : ct);
  ncw_def_dim(name, cdfid, "nobc_max", ct, &nobc);

  ct = 0;
  for (n=1; n<=nwins; n++)
    for (nb = 0; nb < window[n]->nobc; nb++)
      ct = max(ct, window[n]->open[nb]->no2_t);
  ncw_def_dim(name, cdfid, "no2_t_max", ct, &no2_t);

  ct = 0;
  for (n=1; n<=nwins; n++)
    for (nb = 0; nb < window[n]->nobc; nb++) {
      i =  window[n]->open[nb]->no2_e1 + window[n]->open[nb]->to2_e1 - 
	window[n]->open[nb]->no3_e1;
      ct = max(ct, i);
    }
  ncw_def_dim(name, cdfid, "no2_e1_max", ct, &no2_e1);

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
  /* window->nee[szeS]                                                */
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

  /* OBC */
  dims[1] = nobc;
  ncw_def_var2(name, cdfid, "no2_t_w",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number of OBC cell centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "no2_e1_w",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number of OBC normal and tangential  edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  dims[1] = nobc;
  memset(tmask, 0, (MAXCONSTIT+1) * sizeof(int));
  for (n=2; n<=nwins; n++) {
    for (nb = 0; nb < window[n]->nobc; nb++) {
      open_bdrys_t *open = window[n]->open[nb];
      tidal_consts_t *tc = &open->tc;
      if (tc == NULL) continue;

      for (i = 1; i <= tc->nt; i++) {
	char units[MAXSTRLEN], lname[MAXSTRLEN];
	int sz = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	if (tmask[i]) continue;
	strncpy(tname, tc->tname[i], sz);
	tname[sz] = '\0';

	/* Amplitude and phase */
	dims[2] = no2_t;
	sprintf(key, "%s_amp",tname);
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	sprintf(key, "%s amplitude",tname);
	write_text_att(cdfid, vid, "long_name", key);
	write_text_att(cdfid, vid, "units", "m");
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");

	sprintf(key, "%s_phase",tname);
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	sprintf(key, "%s phase",tname);
	write_text_att(cdfid, vid, "long_name", key);
	write_text_att(cdfid, vid, "units", "degrees");
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");
	tmask[i] = 1;

	/* East velocity amplitude and phase */
	dims[2] = no2_e1;
	if (master->tidef & TD_TRAN) {
	  sprintf(key, "%s_amp_U",tname);
	  sprintf(lname, "%s eastward transport amplitude", tname);
	  strcpy(units, "m2s-1");
	  strcpy(tnameu, "U");
	} else {
	  sprintf(key, "%s_amp_u",tname);
	  sprintf(lname, "%s eastward velocity amplitude", tname);
	  strcpy(units, "ms-1");
	  strcpy(tnameu, "u");
	}
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	write_text_att(cdfid, vid, "long_name", lname);
	write_text_att(cdfid, vid, "units", units);
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");

	if (master->tidef & TD_TRAN) {
	  sprintf(key, "%s_phase_U",tname);
	  sprintf(lname, "%s eastward transport phase", tname);
	} else {
	  sprintf(key, "%s_phase_u",tname);
	  sprintf(lname, "%s eastward velocity phase", tname);
	}
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	write_text_att(cdfid, vid, "long_name", lname);
	write_text_att(cdfid, vid, "units", "degrees");
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");

	/* North velocity amplitude and phase */
	if (master->tidef & TD_TRAN) {
	  sprintf(key, "%s_amp_V",tname);
	  sprintf(lname, "%s northward transport amplitude", tname);
	  strcpy(units, "m2s-1");
	  strcpy(tnamev, "V");
	} else {
	  sprintf(key, "%s_amp_v",tname);
	  sprintf(lname, "%s northward velocity amplitude", tname);
	  strcpy(units, "ms-1");
	  strcpy(tnamev, "v");
	}
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	write_text_att(cdfid, vid, "long_name", lname);
	write_text_att(cdfid, vid, "units", units);
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");

	if (master->tidef & TD_TRAN) {
	  sprintf(key, "%s_phase_V",tname);
	  sprintf(lname, "%s northward transport phase", tname);
	} else {
	  sprintf(key, "%s_phase_v",tname);
	  sprintf(lname, "%s northward velocity phase", tname);
	}
	ncw_def_var2(name, cdfid, key,  NC_DOUBLE, 3, dims, &vid, 1);
	write_text_att(cdfid, vid, "long_name", lname);
	write_text_att(cdfid, vid, "units", "degrees");
	write_text_att(cdfid, vid, "cartesian_axis", "C2D");
      }
    }
  }

  /* Global attributes */
  write_text_att(cdfid, NC_GLOBAL, "Parameter_filename", iname);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "WINMAP");
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

    /* OBC */
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nobc_w"), start, count, &window[n]->nobc);

    count[1] = window[n]->nobc;
    for (nb = 0; nb < window[n]->nobc; nb++) {
     window[n]->wincon->i7[nb] = window[n]->open[nb]->no2_t;
    }
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "no2_t_w"), start, count, window[n]->wincon->i7);

    for (nb = 0; nb < window[n]->nobc; nb++) {
      i = window[n]->open[nb]->no2_e1 + window[n]->open[nb]->to2_e1 -
	window[n]->open[nb]->no3_e1;
      window[n]->wincon->i7[nb] = i;
    }
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "no2_e1_w"), start, count, window[n]->wincon->i7);

    for (nb = 0; nb < window[n]->nobc; nb++) {
      open_bdrys_t *open = window[n]->open[nb];
      tidal_consts_t *tc = &open->tc;

      if (tc == NULL) continue;
      start[1] = nb;
      count[1] = 1L;
      memset(window[n]->wincon->d1, 0, window[n]->szm * sizeof(double));
      memset(window[n]->wincon->d1, 0, window[n]->szm * sizeof(double));

      for (i = 1; i <= tc->nt; i++) {
	char tname[5];
	int sz = strlen(tc->tname[i]) - strlen(strpbrk(tc->tname[i], " "));
	strncpy(tname, tc->tname[i], sz);
	tname[sz] = '\0';

	/* Elevation amplitude and phase */
	tc = &open->tc;
	count[2] = window[n]->open[nb]->no2_t;
	if (tc != NULL) {
	  for (cc = 1; cc <= open->no2_t; cc++) {
	    int ci = (tc->map != NULL) ? tc->map[cc] : cc;
	    window[n]->wincon->d1[cc-1] = tc->amp[ci][i];
	    window[n]->wincon->d2[cc-1] = tc->pha[ci][i];
	  }
	  sprintf(key, "%s_amp",tname);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, count,
			      window[n]->wincon->d1);
	  sprintf(key, "%s_phase",tname);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, count,
			      window[n]->wincon->d2);
	}

	/* Eastward normal and tangential velocity amplitude and     */
	/* phase.                                                    */
	vc = 0;
	tc = &open->tun;
	count[2] = window[n]->open[nb]->no2_e1;
	if (tc != NULL) {
	  for (cc = 1; cc <= open->no2_e1; cc++) {
	    int ci = (tc->map != NULL) ? tc->map[cc] : cc;
	    window[n]->wincon->d1[vc] = tc->amp[ci][i];
	    window[n]->wincon->d2[vc] = tc->pha[ci][i];
	    vc++;
	  }
	}
	tc = &open->tut;
	count[2] += (window[n]->open[nb]->to2_e1 - 
		     window[n]->open[nb]->no3_e1);
	if (tc != NULL) {
	  for (cc = open->no3_e1+1; cc <= open->to2_e1; cc++) {
	    int ci = (tc->map != NULL) ? tc->map[cc] : cc;
	    window[n]->wincon->d1[vc] = tc->amp[ci][i];
	    window[n]->wincon->d2[vc] = tc->pha[ci][i];
	    vc++;
	  }
	  sprintf(key, "%s_amp_%s",tname, tnameu);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, 
			      count, window[n]->wincon->d1);
	  sprintf(key, "%s_phase_%s",tname, tnameu);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, 
			      count, window[n]->wincon->d2);
	}

	/* Northward normal and tangential velocity amplitude and    */
	/* phase.                                                    */
	vc = 0;
	tc = &open->tvn;
	count[2] = window[n]->open[nb]->no2_e1;
	if (tc != NULL) {
	  for (cc = 1; cc <= open->no2_e1; cc++) {
	    int ci = (tc->map != NULL) ? tc->map[cc] : cc;
	    window[n]->wincon->d1[vc] = tc->amp[ci][i];
	    window[n]->wincon->d2[vc] = tc->pha[ci][i];
	    vc++;
	  }
	}
	tc = &open->tvt;
	count[2] += (window[n]->open[nb]->to2_e1 - 
		     window[n]->open[nb]->no3_e1);
	if (tc != NULL) {
	  for (cc = open->no3_e1+1; cc <= open->to2_e1; cc++) {
	    int ci = (tc->map != NULL) ? tc->map[cc] : cc;
	    window[n]->wincon->d1[vc] = tc->amp[ci][i];
	    window[n]->wincon->d2[vc] = tc->pha[ci][i];
	    vc++;
	  }
	  sprintf(key, "%s_amp_%s",tname, tnamev);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, 
			      count, window[n]->wincon->d1);
	  sprintf(key, "%s_phase_%s",tname, tnamev);
	  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, key), start, 
			      count, window[n]->wincon->d2);
	}
      }
    }
  }

  ncw_close(name, cdfid);
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
  int c, cc, cg, ee, *d1;
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

  /* Set the global map geom->fm for OBC ghost cells. This is        */
  /* reset in OBC_build after obc_e2 is found.                       */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      geom->fm[c].wn = 0;
      geom->fm[c].sc = 0;
    }
  }

  /* Get layers info and check                                       */
  layers = d_alloc_1d(geom->nz + 1);
  count[0] = geom->nz + 1;;
  sprintf(key, "z_grid");
  nc_get_vara_double(fid, ncw_var_id(fid, key), start, count, layers);
  for (cc = 0; cc < geom->nz; cc++) {
    if (layers[cc] != geom->layers[cc]) {
      hd_warn("Layer %d : layer depths = %f (.prm) %f (file)\n", cc,
	     geom->layers[cc], layers[cc]);
      hd_quit("read_windows_us: Incompatible layer depths in window map file %s with parameter file.\n", name);
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
    /* window->nee[szeS]                                              */
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
	hd_warn("Window%d (%d)=%5.2f : Global (%d)=%5.2f\n", n,
		c, botz[c], 
		cg, geom->botz[cg]);
	hd_quit("read_windows_us: Incompatible bottom depths in window map file %s with parameter file: check runlog.\n", name);
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


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/


#define C2D  1
#define E2D  2
#define V2D  4
#define C3D  8
#define E3D  16
#define V3D  32

/*-------------------------------------------------------------------*/
/* Checks a window map against the maps generated inline.            */
/*-------------------------------------------------------------------*/
void check_window_map_us(geometry_t *geom, geometry_t **window, char *name)
{
  geometry_t **win;
  int e, ee, n, nn, j, c, c1, c2, cc, vv, v, wn;

  win = (geometry_t **)p_alloc_1d(window[1]->nwindows);
  read_windows_us(geom, win, name);

  for (n = 1; n <= geom->nwindows; n++) {
    printf("Checking window %d\n", n);
    if (window[n]->npem != win[n]->npem)
      hd_quit("npem error %d:%d\n", window[n]->npem, win[n]->npem);
    if (window[n]->neem != win[n]->neem)
      hd_quit("neem error %d:%d\n", window[n]->neem, win[n]->neem);
    if (window[n]->nvem != win[n]->nvem)
      hd_quit("nvem error %d:%d\n", window[n]->nvem, win[n]->nvem);
    if (window[n]->nvcm != win[n]->nvcm)
      hd_quit("nvcm error %d:%d\n", window[n]->nvcm, win[n]->nvcm);

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
    if (window[n]->v2_t != win[n]->v2_t)
      hd_quit("v2_t %d:%d\n", window[n]->v2_t, win[n]->v2_t);
    if (window[n]->b2_t != win[n]->b2_t)
      hd_quit("b2_t %d:%d\n", window[n]->b2_t, win[n]->b2_t);
    if (window[n]->a2_t != win[n]->a2_t)
      hd_quit("a2_t %d:%d\n", window[n]->a2_t, win[n]->a2_t);
    if (window[n]->n2_t != win[n]->n2_t)
      hd_quit("n2_t %d:%d\n", window[n]->n2_t, win[n]->n2_t);
    for (cc = 1; cc <= window[n]->n2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (window[n]->w2_t[cc] != win[n]->w2_t[cc]) 
	check_print(window[n], "w2_t", cc, c, window[n]->w2_t[cc], win[n]->w2_t[cc], C2D);
      if (window[n]->npe[c] != win[n]->npe[c]) 
	check_print(window[n], "npe", cc, c, window[n]->npe[c], win[n]->npe[c], C2D);
      if (window[n]->cellarea[c] != win[n]->cellarea[c]) 
	check_print(window[n], "cellarea", cc, c, window[n]->cellarea[c], win[n]->cellarea[c], C2D);
      for (nn = 1; nn <= window[n]->npem; nn++) {
	if (window[n]->vIc[nn][c] != win[n]->vIc[nn][c]) 
	  check_print(window[n], "vIc", cc, c, window[n]->vIc[nn][c], win[n]->vIc[nn][c], C2D);
	if (window[n]->eSc[nn][c] != win[n]->eSc[nn][c]) 
	  check_print(window[n], "eSc", cc, c, window[n]->eSc[nn][c], win[n]->eSc[nn][c], C2D);
      }
    }
    for (cc = 1; cc <= window[n]->a2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (window[n]->bot_t[cc] != win[n]->bot_t[cc]) 
	check_print(window[n], "bot_t", cc, c, window[n]->bot_t[cc], win[n]->bot_t[cc], C2D);
      /* Only nsur_t is set at initialisation
      if (window[n]->sur_t[cc] != win[n]->sur_t[cc]) hd_quit ("sur_t error @ c = %d (%d & %d)\n", c, window[n]->sur_t[cc], win[n]->sur_t[cc]);
      */
      if (window[n]->nsur_t[cc] != win[n]->nsur_t[cc]) 
	check_print(window[n], "sur_t", cc, c, window[n]->nsur_t[cc], win[n]->nsur_t[cc], C2D);
    }

    if (window[n]->v3_t != win[n]->v3_t)
      hd_quit("v3_t %d:%d\n", window[n]->v3_t, win[n]->v3_t);
    if (window[n]->b3_t != win[n]->b3_t)
      hd_quit("b3_t %d:%d\n", window[n]->b3_t, win[n]->b3_t);
    if (window[n]->a3_t != win[n]->a3_t)
      hd_quit("a3_t %d:%d\n", window[n]->a3_t, win[n]->a3_t);
    if (window[n]->n3_t != win[n]->n3_t)
      hd_quit("n3_t %d:%d\n", window[n]->n3_t, win[n]->n3_t);
    for (cc = 1; cc <= window[n]->n3_t; cc++) {
      c = window[n]->w3_t[cc];
      if (window[n]->w3_t[cc] != win[n]->w3_t[cc])
	check_print(window[n], "w3_t", cc, c, window[n]->w3_t[cc], win[n]->w3_t[cc], C3D);
      if (window[n]->zm1[c] != win[n]->zm1[c])
	check_print(window[n], "zm1", cc, c, window[n]->zm1[c], win[n]->zm1[c], C3D);
      if (window[n]->zp1[c] != win[n]->zp1[c])
	check_print(window[n], "zp1", cc, c, window[n]->zp1[c], win[n]->zp1[c], C3D);
      if (window[n]->m2d[c] != win[n]->m2d[c]) 
	check_print(window[n], "m2d", cc, c, window[n]->m2d[c], win[n]->m2d[c], C3D);
      if (window[n]->wsa[c] != win[n]->wsa[c]) 
	check_print(window[n], "wsa", cc, c, window[n]->wsa[c], win[n]->wsa[c], C3D);
      if (window[n]->wgst[c] != win[n]->wgst[c]) 
	check_print(window[n], "wgst", cc, c, window[n]->wgst[c], win[n]->wgst[c], C3D);
      if (window[n]->s2k[c] != win[n]->s2k[c]) 
	check_print(window[n], "s2k", cc, c, window[n]->s2k[c], win[n]->s2k[c], C3D);
      for (nn = 1; nn <= window[n]->npem; nn++) {
	if (window[n]->c2c[nn][c] != win[n]->c2c[nn][c]) hd_quit ("c2c error @ c = %d\n", c);
	if (window[n]->c2e[nn][c] != win[n]->c2e[nn][c]) hd_quit ("c2e error @ c = %d\n", c);
	if (window[n]->c2v[nn][c] != win[n]->c2v[nn][c]) hd_quit ("c2v error @ c = %d\n", c);
      }
    }
    if (window[n]->nbpt != win[n]->nbpt)
      hd_quit("nbpt %d:%d\n", window[n]->nbpt, win[n]->nbpt);
    for (cc = 1; cc <= window[n]->nbpt; cc++) {
      if (window[n]->bpt[cc] != win[n]->bpt[cc]) hd_quit ("bpt error @ c = %d\n", cc);
      if (window[n]->bin[cc] != win[n]->bin[cc]) hd_quit ("bin error @ c = %d\n", cc);
      if (window[n]->bin2[cc] != win[n]->bin2[cc]) hd_quit ("bin2 error @ c = %d\n", cc);
      if (window[n]->dbpt[cc] != win[n]->dbpt[cc]) hd_quit ("dbpt error @ c = %d\n", cc);
    }
    if (window[n]->nm2s != win[n]->nm2s)
      hd_quit("nm2s %d:%d\n", window[n]->nm2s, win[n]->nm2s);
    for (cc = 1; cc <= window[n]->nm2s; cc++)
      if (window[n]->m2s[cc] != win[n]->m2s[cc]) hd_quit ("m2s error @ c = %d\n", cc);
    if (window[n]->ns2m != win[n]->ns2m)
      hd_quit("ns2m %d:%d\n", window[n]->ns2m, win[n]->ns2m);
    for (cc = 1; cc <= window[n]->ns2m; cc++)
      if (window[n]->s2m[cc] != win[n]->s2m[cc]) hd_quit ("s2m error @ c = %d\n", cc);

    /* EDGES                                                         */ 
    /* window->nee[szeS]                                              */
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
    if (window[n]->v2_e1 != win[n]->v2_e1)
      hd_quit("v2_e1 %d:%d\n", window[n]->v2_e1, win[n]->v2_e1);
    if (window[n]->b2_e1 != win[n]->b2_e1)
      hd_quit("b2_e1 %d:%d\n", window[n]->b2_e1, win[n]->b2_e1);
    if (window[n]->a2_e1 != win[n]->a2_e1)
      hd_quit("a2_e1 %d:%d\n", window[n]->a2_e1, win[n]->a2_e1);
    if (window[n]->n2_e1 != win[n]->n2_e1)
      hd_quit("n2_e1 %d:%d\n", window[n]->n2_e1, win[n]->n2_e1);
    for (ee = 1; ee <= window[n]->n2_e1; ee++) {
      e = window[n]->w2_e1[ee];
      if (window[n]->w2_e1[ee] != win[n]->w2_e1[ee]) 
	check_print(window[n], "w2_e1", ee, e, window[n]->w2_e1[ee], win[n]->w2_e1[ee], E2D);
      if (window[n]->bot_e1[ee] != win[n]->bot_e1[ee])
	check_print(window[n], "bot_e1", ee, e, window[n]->bot_e1[ee], win[n]->bot_e1[ee], E2D);
      if (window[n]->sur_e1[ee] != win[n]->sur_e1[ee])
	check_print(window[n], "sur_e1", ee, e, window[n]->sur_e1[ee], win[n]->sur_e1[ee], E2D);
      if (window[n]->nee[e] != win[n]->nee[e]) 
	check_print(window[n], "nee", ee, e, window[n]->nee[e], win[n]->nee[e], E2D);
    }
    if (window[n]->v3_e1 != win[n]->v3_e1)
      hd_quit("v3_e1 %d:%d\n", window[n]->v3_e1, win[n]->v3_e1);
    if (window[n]->b3_e1 != win[n]->b3_e1)
      hd_quit("b3_e1 %d:%d\n", window[n]->b3_e1, win[n]->b3_e1);
    if (window[n]->a3_e1 != win[n]->a3_e1)
      hd_quit("a3_e1 %d:%d\n", window[n]->a3_e1, win[n]->a3_e1);
    if (window[n]->n3_e1 != win[n]->n3_e1)
      hd_quit("n3_e1 %d:%d\n", window[n]->n3_e1, win[n]->n3_e1);
    for (ee = 1; ee <= window[n]->n3_e1; ee++) {
      e = window[n]->w3_e1[ee];
      if (window[n]->w3_e1[ee] != win[n]->w3_e1[ee]) 
	check_print(window[n], "w3_e1", ee, e, window[n]->w3_e1[ee], win[n]->w3_e1[ee], E3D);
      if (window[n]->wse[e] != win[n]->wse[e]) hd_quit("wse error @ %d\n", e);
      if (window[n]->ep[e] != win[n]->ep[e]) hd_quit("ep error @ %d\n", e);
      if (window[n]->em[e] != win[n]->em[e]) hd_quit("em error @ %d\n", e);
      if (window[n]->zm1e[e] != win[n]->zm1e[e]) hd_quit("zm1e error @ %d\n", e);
      if (window[n]->zp1e[e] != win[n]->zp1e[e]) hd_quit("zp1e error @ %d\n", e);
      if (window[n]->m2de[e] != win[n]->m2de[e]) hd_quit("m2de error @ %d\n", e);
      if (window[n]->e2k[e] != win[n]->e2k[e]) hd_quit("e2k error @ %d\n", e);
      if (window[n]->e2c[e][0] != win[n]->e2c[e][0]) hd_quit("e2c[0] error @ %d\n", e);
      if (window[n]->e2c[e][1] != win[n]->e2c[e][1]) hd_quit("e2c[1] error @ %d\n", e);
      if (window[n]->e2e[e][0] != win[n]->e2e[e][0]) hd_quit("e2e[0] error @ %d\n", e);
      if (window[n]->e2e[e][1] != win[n]->e2e[e][1]) hd_quit("e2e[1] error @ %d\n", e);
      if (window[n]->e2v[e][0] != win[n]->e2v[e][0]) hd_quit("e2v[0] error @ %d\n", e);
      if (window[n]->e2v[e][1] != win[n]->e2v[e][1]) hd_quit("e2v[1] error @ %d\n", e);
      for (nn = 1; nn <= window[n]->neem; nn++) {
	if (window[n]->eSe[nn][e] != win[n]->eSe[nn][e]) hd_quit ("eSe error @ e = %d\n", e);
	if (window[n]->wAe[nn][e] != win[n]->wAe[nn][e]) hd_quit ("wAe error @ e = %d\n", e);
      }
    }
    if (window[n]->nbpte1 != win[n]->nbpte1)
      hd_quit("nbpte1 %d:%d\n", window[n]->nbpte1, win[n]->nbpte1);
    for (ee = 1; ee <= window[n]->nbpt; ee++) {
      if (window[n]->bpte1[ee] != win[n]->bpte1[ee]) hd_quit ("bpte1 error @ e = %d\n", ee);
      if (window[n]->bine1[ee] != win[n]->bine1[ee]) hd_quit ("bine1 error @ e = %d\n", ee);
    }
    if (window[n]->nm2se1 != win[n]->nm2se1)
      hd_quit("nm2se1 %d:%d\n", window[n]->nm2se1, win[n]->nm2se1);
    for (ee = 1; ee <= window[n]->nm2se1; ee++)
      if (window[n]->m2se1[ee] != win[n]->m2se1[ee]) hd_quit ("m2se1 error @ e = %d\n", ee);
    if (window[n]->ns2me1 != win[n]->ns2me1)
      hd_quit("ns2me1 %d:%d\n", window[n]->ns2me1, win[n]->ns2me1);
    for (ee = 1; ee <= window[n]->ns2me1; ee++)
      if (window[n]->s2me1[ee] != win[n]->s2me1[ee]) hd_quit ("s2me1 error @ e = %d\n", ee);

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
    if (window[n]->v2_e2 != win[n]->v2_e2)
      hd_quit("v2_e2 %d:%d\n", window[n]->v2_e2, win[n]->v2_e2);
    if (window[n]->b2_e2 != win[n]->b2_e2)
      hd_quit("b2_e2 %d:%d\n", window[n]->b2_e2, win[n]->b2_e2);
    if (window[n]->a2_e2 != win[n]->a2_e2)
      hd_quit("a2_e2 %d:%d\n", window[n]->a2_e2, win[n]->a2_e2);
    if (window[n]->n2_e2 != win[n]->n2_e2)
      hd_quit("n2_e2 %d:%d\n", window[n]->n2_e2, win[n]->n2_e2);
    for (vv = 1; vv <= window[n]->n2_e2; vv++) {
      vv = window[n]->w2_e2[vv];
      if (window[n]->w2_e2[vv] != win[n]->w2_e2[vv]) 
	check_print(window[n], "w2_e2", vv, v, window[n]->w2_e2[vv], win[n]->w2_e2[vv], V2D);
      /*
      if (window[n]->bot_e2[vv] != win[n]->bot_e2[vv]) 
	check_print(window[n], "bot_e2", vv, v, window[n]->bot_e2[vv], win[n]->bot_e2[vv], V2D;
      if (window[n]->sur_e2[vv] != win[n]->sur_e2[vv]) 
	check_print(window[n], "sur_e2", vv, v, window[n]->sur_e2[vv], win[n]->sur_e2[vv], V2D;
      */
      if (window[n]->nvc[v] != win[n]->nvc[v]) hd_quit ("nvc error @ v = %d\n", v);
      if (window[n]->nve[v] != win[n]->nve[v]) hd_quit ("nve error @ v = %d\n", v);
      if (window[n]->dualarea[v] != win[n]->dualarea[v]) hd_quit ("dualarea error @ v = %d\n", v);
    }
    if (window[n]->v3_e2 != win[n]->v3_e2)
      hd_quit("v3_e2 %d:%d\n", window[n]->v3_e2, win[n]->v3_e2);
    if (window[n]->b3_e2 != win[n]->b3_e2)
      hd_quit("b3_e2 %d:%d\n", window[n]->b3_e2, win[n]->b3_e2);
    if (window[n]->a3_e2 != win[n]->a3_e2)
      hd_quit("a3_e2 %d:%d\n", window[n]->a3_e2, win[n]->a3_e2);
    if (window[n]->n3_e2 != win[n]->n3_e2)
      hd_quit("n3_e2 %d:%d\n", window[n]->n3_e2, win[n]->n3_e2);
    for (vv = 1; vv <= window[n]->n3_e2; vv++) {
      vv = window[n]->w3_e2[vv];
      if (window[n]->w3_e2[vv] != win[n]->w3_e2[vv]) 
	check_print(window[n], "w3_e2", vv, v, window[n]->w3_e2[vv], win[n]->w3_e2[vv], V3D);
      if (window[n]->wsv[v] != win[n]->wsv[v]) hd_quit ("wsv error @ v = %d\n", v);
      if (window[n]->zm1v[v] != win[n]->zm1v[v]) hd_quit ("zm1v error @ v = %d\n", v);
      if (window[n]->zp1v[v] != win[n]->zp1v[v]) hd_quit ("zp1v error @ v = %d\n", v);
      if (window[n]->m2dv[v] != win[n]->m2dv[v]) hd_quit ("m2dv error @ v = %d\n", v);
      for (nn = 1; nn <= window[n]->nvcm; nn++) {
	if (window[n]->v2c[v][nn] != win[n]->v2c[v][nn]) hd_quit ("v2c error @ v = %d\n", v);
	if (window[n]->dualareap[v][nn] != win[n]->dualareap[v][nn]) hd_quit ("dualareap error @ v = %d\n", v);
      }
      for (nn = 1; nn <= window[n]->nvem; nn++) {
	if (window[n]->v2e[v][nn] != win[n]->v2e[v][nn]) hd_quit ("v2e error @ v = %d\n", v);
	if (window[n]->eSv[nn][v] != win[n]->eSv[nn][v]) hd_quit ("eSv error @ v = %d\n", v);
      }
    }

  }

  /* Open boundaries                                                  */ 
  OBC_build(geom->open, geom, win, geom->nwindows);

  for (wn = 1; wn <= geom->nwindows; wn++) {
    printf("Checking boundaries for window %d\n", wn);
    if (window[wn]->nobc != win[wn]->nobc)
      check_print_obc("nobc", n, 0, window[n]->nobc, win[n]->nobc);
    for (n = 0; n < window[wn]->nobc; n++) {
      int cs;
      open_bdrys_t *open = window[wn]->open[n];
      open_bdrys_t *openw = win[wn]->open[n];
      if (open->bgz != openw->bgz)
	check_print_obc("bgz", n, 0, open->bgz, openw->bgz);
      if(open->no2_t != openw->no2_t)
	check_print_obc("no2_t", 0, n, open->no2_t, openw->no2_t);
      if(open->no2_e1 != openw->no2_e1)
	check_print_obc("no2_e1", n, 0, open->no2_e1, openw->no2_e1);
      if(open->no2_e2 != openw->no2_e2)
	check_print_obc("no2_e2", n, 0, open->no2_e2, openw->no2_e2);
      if(open->no3_t != openw->no3_t)
	check_print_obc("no3_t", n, 0, open->no3_t, openw->no3_t);
      for (cc = 1; cc <= open->no3_t; cc++) {
	if (open->obc_t[cc] != openw->obc_t[cc]) 
	  check_print_obc("obc_t", n, cc, open->obc_t[cc], openw->obc_t[cc]);
	if (open->oi1_t[cc] != openw->oi1_t[cc]) 
	  check_print_obc("oi1_t", n, cc, open->oi1_t[cc], openw->oi1_t[cc]);
	if (open->oi2_t[cc] != openw->oi2_t[cc])
	  check_print_obc("oi2_t", n, cc, open->oi2_t[cc], openw->oi2_t[cc]);
	if (open->nepc[cc] != openw->nepc[cc])
	  check_print_obc("nepc", n, cc, open->nepc[cc], openw->nepc[cc]);
	if (open->inc[cc] != openw->inc[cc])
	  check_print_obc("inc", n, cc, open->inc[cc], openw->inc[cc]);
	if (open->olap[cc] != openw->olap[cc])
	  check_print_obc("olap", n, cc, open->olap[cc], openw->olap[cc]);
	if (cc <= open->no2_t && open->bot_t[cc] != openw->bot_t[cc])
	  check_print_obc("bot_t", n, cc, open->bot_t[cc], openw->bot_t[cc]);
	c = open->obc_t[cc];
	cs = window[wn]->m2d[c];
	for (j = 1; j <= window[wn]->npe[cs]; j++) {
	  if (open->bec[j][cc] != openw->bec[j][cc])
	    check_print_obc("bec", n, cc, open->bec[j][cc], openw->bec[j][cc]);
	  if (open->bcc[j][cc] != openw->bcc[j][cc])
	    check_print_obc("bcc", n, cc, open->bcc[j][cc], openw->bcc[j][cc]);
	}
      }
      if(open->no2_e1 != openw->no2_e1)
	check_print_obc("no2_e1", n, 0, open->no2_e1, openw->no2_e1);
      if(open->to2_e1 != openw->to2_e1)
	check_print_obc("to2_e1", n, 0, open->to2_e1, openw->to2_e1);
      if(open->no3_e1 != openw->no3_e1)
	check_print_obc("no3_e1", n, 0, open->no3_e1, openw->no3_e1);
      if(open->to3_e1 != openw->to3_e1)
	check_print_obc("to3_e1", n, 0, open->to3_e1, openw->to3_e1);
      for (ee = 1; ee <= open->to3_e1; ee++) {
	if (open->obc_e1[ee] != openw->obc_e1[ee]) 
	  check_print_obc("obc_e1", n, ee, open->obc_e1[ee], openw->obc_e1[ee]);
	if (open->obc_e2[ee] != openw->obc_e2[ee]) 
	  check_print_obc("obc_e2", n, ee, open->obc_e2[ee], openw->obc_e2[ee]);
	if (open->oi1_e1[ee] != openw->oi1_e1[ee]) 
	  check_print_obc("oi1_e1", n, ee, open->oi1_e1[ee], openw->oi1_e1[ee]);
	if (open->oi2_e1[ee] != openw->oi2_e1[ee]) 
	  check_print_obc("oi2_e1", n, ee, open->oi2_e1[ee], openw->oi2_e1[ee]);
	if (ee <= open->no3_e1 && open->ogc_t[ee] != openw->ogc_t[ee]) 
	  check_print_obc("ogc_t", n, ee, open->ogc_t[ee], openw->ogc_t[ee]);
	if (ee <= open->no3_e1 && open->ceni[ee] != openw->ceni[ee]) 
	  check_print_obc("ceni", n, ee, open->ceni[ee], openw->ceni[ee]);
	if (ee <= open->no3_e1 && open->ini[ee] != openw->ini[ee]) 
	  check_print_obc("ini", n, ee, open->ini[ee], openw->ini[ee]);
	if (ee <= open->no3_e1 && open->outi[ee] != openw->outi[ee]) 
	  check_print_obc("outi", n, ee, open->outi[ee], openw->outi[ee]);
	if (ee <= open->no3_e1 && open->dir[ee] != openw->dir[ee]) 
	  check_print_obc("dir", n, ee, open->dir[ee], openw->dir[ee]);

      }
    }
  }

  /* unlink(name); */
  /*
  hd_quit("done check_window\n");
  */
  printf("done check_window\n");
}

/* END check_window_map_us()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Dumps window geometry to file                                     */
/* Note: some of the dimensions declared are of the array size, e.g. */
/* n2_t dimension is written as n2_t+1, this must be decremented if  */
/* read back in.                                                     */
/*-------------------------------------------------------------------*/
void dump_geom_us(master_t *master, char *iname)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  parameters_t *params = master->params;
  mesh_t *m = params->mesh;
  int cdfid, n, nb, c, cc, vc, i;
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
  int kgridid;          /* K dimension id at grid corner */
  int oid, tid;
  int nwins = geom->nwindows;
  char key[MAXSTRLEN];
  char name[MAXSTRLEN];
  int *d1;
  int ct;
  int ncmode = overwrite_output() ? NC_CLOBBER:NC_NOCLOBBER;
  int npem, nvem, nvcm;
  int nobc;
  int no2_t, no2_e1;
  int ns2 = m->ns2;

  /*-----------------------------------------------------------------*/
  /* Create the netCDF file                                          */
#ifdef NC_NETCDF4
  ncmode |= NC_NETCDF4;
#else
  ncmode |= NC_64BIT_OFFSET;
#endif
  for (i = 0; i < strlen(iname) - 4; i++)
    key[i] = iname[i];
  sprintf(name, "%s.map", key);

  if (ncw_create(name, ncmode, &cdfid) != NC_NOERR)
    hd_quit("Couldn't create output dump file %s\n", name);

  /* Define global dimensions                                        */
  ncw_def_dim(name, cdfid, "szc",    geom->szc,  &szcid);
  ncw_def_dim(name, cdfid, "szcS",   geom->szcS, &szcSid);
  ncw_def_dim(name, cdfid, "sze",    geom->sze,  &szeid);
  ncw_def_dim(name, cdfid, "szeS",   geom->szeS, &szeSid);
  ncw_def_dim(name, cdfid, "szv",    geom->szv,  &szvid);
  ncw_def_dim(name, cdfid, "szvS",   geom->szvS, &szvSid);
  ncw_def_dim(name, cdfid, "k_grid",   geom->nz + 1, &kgridid);

  /* Variables of length 1                                           */
  ncw_def_dim(name, cdfid, "one", 1, &oid);
  ncw_def_dim(name, cdfid, "two", 2, &tid);
  
  /* Global variables                                                */
  dims[0] = kgridid;
  ncw_def_var(name, cdfid, "z_grid", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate at grid layer faces");  
  write_text_att(cdfid, vid, "coordinate_type", "Z");
  dims[0] = szcid;
  ncw_def_var(name, cdfid, "cellz", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate at cell centre");  
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  ncw_def_var(name, cdfid, "gridz", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate at grid layer faces");  
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  dims[0] = oid;
  ncw_def_var(name, cdfid, "sednz", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Number of Sediment layers");

  /* Geometry variables                                              */
  ncw_def_dim(name, cdfid, "enon", geom->enon, &vid);
  ncw_def_dim(name, cdfid, "enonS", geom->enonS, &vid);
  ncw_def_dim(name, cdfid, "ewet", geom->enon, &vid);
  ncw_def_dim(name, cdfid, "ewetS", geom->enon, &vid);
  ncw_def_dim(name, cdfid, "snon", geom->enon, &vid);
  ncw_def_dim(name, cdfid, "snonS", geom->enon, &vid);
  ncw_def_dim(name, cdfid, "v2_t", geom->v2_t, &vid);
  ncw_def_dim(name, cdfid, "b2_t", geom->b2_t, &vid);
  ncw_def_dim(name, cdfid, "v3_t", geom->v3_t, &vid);
  ncw_def_dim(name, cdfid, "b3_t", geom->b3_t, &vid);
  ncw_def_dim(name, cdfid, "v2_e1", geom->v2_e1, &vid);
  ncw_def_dim(name, cdfid, "b2_e1", geom->b2_e1, &vid);
  ncw_def_dim(name, cdfid, "x2_e1", geom->x2_e1, &vid);
  ncw_def_dim(name, cdfid, "v3_e1", geom->v3_e1, &vid);
  ncw_def_dim(name, cdfid, "b3_e1", geom->b3_e1, &vid);
  ncw_def_dim(name, cdfid, "x3_e1", geom->x3_e1, &vid);
  ncw_def_dim(name, cdfid, "v2_e2", geom->v2_e2, &vid);
  ncw_def_dim(name, cdfid, "b2_e2", geom->b2_e2, &vid);
  ncw_def_dim(name, cdfid, "v3_e2", geom->v3_e2, &vid);
  ncw_def_dim(name, cdfid, "b3_e2", geom->b3_e2, &vid);
  ncw_def_dim(name, cdfid, "nobc", geom->nobc, &vid);
  ncw_def_dim(name, cdfid, "bdry", geom->bdry, &vid);

  /* CENTRES                                                         */ 
  /* geom->npe[szcS]                                                 */
  /* geom->bot_t[szcS]                                               */
  /* geom->sur_t[szcS]                                               */
  /* geom->w2_t[szcS]                                                */
  /* geom->zm1[szc]                                                  */
  /* geom->zp1[szc]                                                  */
  /* geom->m2d[szc]                                                  */
  /* geom->wsa[szc]                                                  */
  /* geom->wgst[szc]                                                 */
  /* geom->w3_t[szc]                                                 */
  /* geom->c2c[npem][szc]                                            */
  /* geom->c2e[npem][szc]                                            */
  /* geom->c2v[npem][szc]                                            */
  /* geom->vIc[npem][szcS]                                           */
  /* geom->eSc[npem][szcS]                                           */
  /* geom->bpt[nbpt]                                                 */
  /* geom->bin[nbpt]                                                 */
  /* geom->cellarea[szcS]                                            */
  
  /* npe                                                             */
  ct = geom->npem + 1;
  ncw_def_dim(name, cdfid, "npe", ct, &npem);

  /* obcs                                                            */
  ct = 0;
  for (nb = 0; nb < geom->nobc; nb++)
    ct = max(ct, geom->open[nb]->no2_t);
  ncw_def_dim(name, cdfid, "no2_t", ct, &no2_t);

  ct = 0;
  for (nb = 0; nb < geom->nobc; nb++) {
    i =  geom->open[nb]->no2_e1 + geom->open[nb]->to2_e1 - 
      geom->open[nb]->no3_e1;
    ct = max(ct, i);
  }
  ncw_def_dim(name, cdfid, "no2_e1", ct, &no2_e1);

  /* Surface sparse arrays                                           */
  dims[0] = szcSid;
  ncw_def_var2(name, cdfid, "npe",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centres surrounding a cell");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "cc2s", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "cell index to sparse map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  /* cell area                                                       */
  ncw_def_var2(name, cdfid, "cellarea", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell area");
  write_text_att(cdfid, vid, "units", "m2");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "cellx", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell centre x location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "celly", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell centre y location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  dims[0] = szeSid;
  ncw_def_var2(name, cdfid, "u1x", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge centre x location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "u1y", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge centre y location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "h1au1", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge length between nodes");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "h2au1", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Length between face centres");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "h1acell", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Upstrean edge length between nodes");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "edgearea", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge area");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "thetau1", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Angle between u1 vector and x axis");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "thetau2", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Angle between u1 vector and y axis");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  dims[0] = szvSid;
  ncw_def_var2(name, cdfid, "gridx", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex x location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  ncw_def_var2(name, cdfid, "gridy", NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex y location");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  dims[0] = szcSid;
  ncw_def_var2(name, cdfid, "bot_t",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Bottom centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "sur_t",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Surface centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  /*
  ncw_def_var2(name, cdfid, "nsur_t", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Updated surface centre coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  */
  /* 3D sparse arrays                                                */
  dims[0] = szcid;
  ncw_def_var2(name, cdfid, "wsa",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "m2d", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "zp1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centre map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "zm1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Centre map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "wgst", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "lateral boundary ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "c2cc", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "sparse to cell index map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  /* Maps                                                            */
  dims[0] = npem;
  dims[1] = szcid;
  ncw_def_var2(name, cdfid, "c2c", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "c2e", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "c2v", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Cell to vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  dims[0] = npem;
  dims[1] = szcSid;
  ncw_def_var2(name, cdfid, "vIc", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex index");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "eSc", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge sign");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");


  /* nbt                                                             */
  ct = geom->nbptS+1;
  ncw_def_dim(name, cdfid, "nbptS", ct, &dims[0]);
  ct = geom->nbpt+1;
  ncw_def_dim(name, cdfid, "nbpt", ct, &dims[0]);

  ncw_def_var2(name, cdfid, "bpt",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "bin",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "dbpt", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Direction of the ghost cell");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ncw_def_var2(name, cdfid, "dbin", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Direction of interior");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  /* EDGES                                                             */ 
  /* geom->nee[szeS]                                                   */
  /* geom->wse[sze]                                                    */
  /* geom->e2c[sze][2]                                                 */
  /* geom->e2v[sze][2]                                                 */
  /* geom->e2e[sze][2]                                                 */
  /* geom->eSe[neem][sze]                                              */
  /* geom->wAe[neem][sze]                                              */
  /* geom->ep[sze]                                                     */
  /* geom->em[sze]                                                     */
  /* geom->zm1e[sze]                                                   */
  /* geom->zp1e[sze]                                                   */
  /* geom->e2k[sze]                                                    */
  /* geom->m2de[sze]                                                   */
  /* geom->bot_e1[szeS]                                                */
  /* geom->sur_e1[szeS]                                                */
  /* geom->bpte1[nbpte1]                                               */
  /* geom->bine1[nbpte1]                                               */
  /* geom->bine2[nbpte1]                                               */
  /* geom->w2_e1[szeS]                                                 */
  /* geom->w3_e1[sze]                                                  */

  dims[0] = szeid;

  ncw_def_var2(name, cdfid, "wse",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "nee",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edges in tangential velocity map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "m2de", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2k", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertical edge number");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2ijk", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to (i,j,k) map for structured grids");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "ep", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map (positive)");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "em", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map (negative)");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "zp1e", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "zm1e", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  if (geom->nobc) {
    ncw_def_var2(name, cdfid, "maske", NC_INT, 1, dims, &vid, 1);
    write_text_att(cdfid, vid, "long_name", "Open boundary edge mask");
    write_text_att(cdfid, vid, "units", "");
    write_text_att(cdfid, vid, "cartesian_axis", "E3D");
  }

  dims[0] = szeid;
  dims[1] = tid;
  ncw_def_var2(name, cdfid, "e2c", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2v", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "e2e", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Edge to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = geom->neem+1;
  ncw_def_dim(name, cdfid, "neem", ct, &dims[0]);

  dims[1] = szeid;
  ncw_def_var2(name, cdfid, "eSe", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Tangential velocity edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "wAe", NC_DOUBLE, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Tangential velocity weights");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  dims[0] = szeSid;
  ncw_def_var2(name, cdfid, "bot_e1",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Bottom edge coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");
  ncw_def_var2(name, cdfid, "sur_e1",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Surface edge coordinate");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ct = geom->nbpte1+1;
  ncw_def_dim(name, cdfid, "nbpte1",  ct,  &dims[0]);
  ncw_def_var2(name, cdfid, "bpte1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D ghost edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "bine1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to 3D ghost edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ncw_def_var2(name, cdfid, "bine2", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to 3D ghost centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = geom->nbpte1S+1;
  ncw_def_dim(name, cdfid, "nbpte1S", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "bpte1S", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D ghost edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ncw_def_var2(name, cdfid, "bine1S", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to 2D ghost edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ct = geom->nbe1;
  ncw_def_dim(name, cdfid, "nbe1",  ct,  &dims[0]);
  ct = geom->nbe1S;
  ncw_def_dim(name, cdfid, "nbe1S",  ct,  &dims[0]);

  /* VERTICES                                                          */ 
  /* geom->nve[szvS]                                                   */
  /* geom->nvc[szvS]                                                   */
  /* geom->wsv[szv]                                                    */
  /* geom->v2c[szv][nvcm]                                              */
  /* geom->v2e[szv][nvem]                                              */
  /* geom->eSv[nvem][szv]                                              */
  /* geom->zm1v[szv]                                                   */
  /* geom->zp1v[szv]                                                   */
  /* geom->dualarea[szvS]                                              */
  /* geom->dualareap[szvS][nvcm]                                       */
  /* geom->m2dv[szv]                                                   */
  /* geom->w2_e2[szvS]                                                 */
  /* geom->w3_e2[szv]                                                  */

  dims[0] = szvid;
  ncw_def_var2(name, cdfid, "wsv",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "m2dv", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D vertex map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "zp1v", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ncw_def_var2(name, cdfid, "zm1v", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  if (geom->v2ijk) {
    ncw_def_var2(name, cdfid, "v2ijk", NC_INT, 1, dims, &vid, 1);
    write_text_att(cdfid, vid, "long_name", "Vertex to (i,j,k) map for structured grids");
    write_text_att(cdfid, vid, "units", "");
    write_text_att(cdfid, vid, "cartesian_axis", "V3D");
  }

  ct = geom->nvem + 1;
  ncw_def_dim(name, cdfid, "nve", ct, &nvem);
  dims[0] = szvid;
  dims[1] = nvem;
  ncw_def_var2(name, cdfid, "v2e", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex to edge map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  ct = geom->nvcm + 1;
  ncw_def_dim(name, cdfid, "nvc", ct, &nvcm);
  dims[0] = szvid;
  dims[1] = nvcm;
  ncw_def_var2(name, cdfid, "v2c", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex to centre map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  dims[0] = nvem;
  dims[1] = szvid;
  ncw_def_var2(name, cdfid, "eSv", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Vertex sign");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  dims[0] = szvSid;
  ncw_def_var2(name, cdfid, "nve",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number edges from a vertex");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  ncw_def_var2(name, cdfid, "nvc",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number centres around a vertex");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  ncw_def_var2(name, cdfid, "dualarea",  NC_DOUBLE, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Area of the dual");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");

  dims[0] = szvSid;
  dims[1] = nvcm;
  ncw_def_var2(name, cdfid, "dualareap",  NC_DOUBLE, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Partial area of the dual");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");
  
  ct = geom->n2_t+1;
  ncw_def_dim(name, cdfid, "n2_t", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w2_t", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  dims[1] = dims[0];
  dims[0] = npem;
  ncw_def_var2(name, cdfid, "hacell", NC_DOUBLE, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Half distance between faces");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");
  
  ct = geom->n2_e1+1;
  ncw_def_dim(name, cdfid, "n2_e1", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w2_e1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E2D");

  ct = geom->n2_e2+1;
  ncw_def_dim(name, cdfid, "n2_e2", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w2_e2", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "2D work array for vertices");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V2D");
  
  ct = geom->n3_t+1;
  ncw_def_dim(name, cdfid, "n3_t", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w3_t", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  ct = geom->n3_e1+1;
  ncw_def_dim(name, cdfid, "n3_e1", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w3_e1", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "E3D");

  ct = geom->n3_e2+1;
  ncw_def_dim(name, cdfid, "n3_e2", ct, &dims[0]);
  ncw_def_var2(name, cdfid, "w3_e2", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D work array for vertices");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "V3D");

  dims[0] = szcid;
  ncw_def_var2(name, cdfid, "s2i", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian i map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  ncw_def_var2(name, cdfid, "s2j", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian j map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");
  ncw_def_var2(name, cdfid, "s2k", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Unstructured to Cartesian k map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C3D");

  dims[0] = szcSid;
  ncw_def_var(name, cdfid, "botz", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Bottom depth");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  /* OBC                                                               */
  dims[0] = nobc;
  ncw_def_var2(name, cdfid, "no2_t",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number of OBC cell centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  ncw_def_var2(name, cdfid, "no2_e1",  NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Number of OBC normal and tangential  edges");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "C2D");

  /* Global attributes                                                 */
  write_text_att(cdfid, NC_GLOBAL, "Parameter_filename", iname);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "WINMAP");
  nc_put_att_int(cdfid, NC_GLOBAL, "Windows", NC_INT, 1, &geom->nwindows);
  write_text_att(cdfid, NC_GLOBAL, "EMS Version", version);
  write_text_att(cdfid, NC_GLOBAL, "Executable", executable);
  getcwd(key, MAXSTRLEN);
  write_text_att(cdfid, NC_GLOBAL, "Directory", key);
  write_text_att(cdfid, NC_GLOBAL, "Parameterheader", params->parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "Grid_description", params->grid_desc);
  time(&t);
  write_text_att(cdfid, NC_GLOBAL, "Date_created", ctime(&t));

  nc_enddef(cdfid);
  /*ncw_enddef(name, cdfid);*/

  /*-----------------------------------------------------------------*/
  /* Write the data                                                  */

  count[0] = geom->nz+1;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "z_grid"), start, count, geom->layers);

  count[0] = geom->szc;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "cellz"), start, count,
                     geom->cellz);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "gridz"), start, count,
                     geom->gridz);

  /* Constrants                                                      */
  start[0] = 0;
  count[0] = 1;
  nc_put_var_int(cdfid, ncw_var_id(cdfid, "sednz"), &geom->sednz);

  /* 1D arrays                                                         */
  count[0] = geom->b2_t + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_t"), start, count, geom->bot_t);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_t"), start, count, geom->sur_t);
  /*ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nsur_t"), start, count, geom->nsur_t);*/

  count[0] = geom->szeS;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_e1"), start, count, geom->bot_e1);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_e1"), start, count, geom->sur_e1);

  count[0] = geom->szc;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2d"), start, count, geom->m2d);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wsa"), start, count, geom->wsa);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1"), start, count, geom->zp1);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1"), start, count, geom->zm1);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wgst"), start, count, geom->wgst);
  
  count[0] = geom->szcS;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "npe"), start, count, geom->npe);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "cellarea"), start, count, geom->cellarea);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "cellx"), start, count, geom->cellx);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "celly"), start, count, geom->celly);
  count[0] = geom->szeS;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "u1x"), start, count, geom->u1x);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "u1y"), start, count, geom->u1y);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "h1au1"), start, count, geom->h1au1);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "h2au1"), start, count, geom->h2au1);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "h1acell"), start, count, geom->h1acell);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "edgearea"), start, count, geom->edgearea);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "thetau1"), start, count, geom->thetau1);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "thetau2"), start, count, geom->thetau2);
  count[0] = geom->szvS;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "gridx"), start, count, geom->gridx);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "gridy"), start, count, geom->gridy);
  count[0] = ns2+1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "cc2s"), start, count, geom->cc2s);

  count[0] = geom->npem+1;
  count[1] = geom->szc;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2c"), start, count, geom->c2c[0]);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2e"), start, count, geom->c2e[0]);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2v"), start, count, geom->c2v[0]);
  count[1] = geom->n2_t+1;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "hacell"), start, count, geom->hacell[0]);

  count[0] = geom->sze;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2de"), start, count, geom->m2de);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wse"), start, count, geom->wse);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1e"), start, count, geom->zp1e);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1e"), start, count, geom->zm1e);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ep"), start, count, geom->ep);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "em"), start, count, geom->em);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2k"), start, count, geom->e2k);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2ijk"), start, count, geom->e2ijk);
  if (geom->nobc)
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "maske"), start, count, geom->maske);

  count[0] = geom->sze;
  count[1] = 2;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2c"), start, count, geom->e2c[0]);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2e"), start, count, geom->e2e[0]);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "e2v"), start, count, geom->e2v[0]);

  count[0] = geom->neem+1;
  count[1] = geom->sze;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSe"), start, count, geom->eSe[0]);
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "wAe"), start, count, geom->wAe[0]);

  count[0] = geom->szeS;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nee"), start, count, geom->nee);

  count[0] = geom->szv;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2dv"), start, count, geom->m2dv);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wsv"), start, count, geom->wsv);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1v"), start, count, geom->zp1v);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1v"), start, count, geom->zm1v);
  if (geom->v2ijk)
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2ijk"), start, count, geom->v2ijk);

  count[0] = geom->szv;
  count[1] = geom->nvcm+1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2c"), start, count, geom->v2c[0]);

  count[0] = geom->szv;
  count[1] = geom->nvem + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v2e"), start, count, geom->v2e[0]);

  count[0] = geom->npem + 1;
  count[1] = geom->szcS;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "vIc"), start, count, geom->vIc[0]);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSc"), start, count, geom->eSc[0]);

  count[0] = geom->szvS;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nve"), start, count, geom->nve);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nvc"), start, count, geom->nvc);

  count[0] = geom->nvem+1;
  count[1] = geom->szv;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "eSv"), start, count, geom->eSv[0]);

  count[0] = geom->szvS;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "dualarea"), start, count, geom->dualarea);

  count[0] = geom->szvS;
  count[1] = geom->nvcm+1;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "dualareap"), start, count,
			geom->dualareap[0]);
    
  count[0] = geom->nbpt + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpt"), start, count, geom->bpt);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bin"), start, count, geom->bin);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "dbpt"), start, count, geom->dbpt);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "dbin"), start, count, geom->dbin);

  count[0] = geom->nbpte1 + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1"), start, count, geom->bpte1);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1"), start, count, geom->bine1);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine2"), start, count, geom->bine2);

  count[0] = geom->nbpte1S + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1S"), start, count, geom->bpte1S);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1S"), start, count, geom->bine1S);

  count[0] = geom->n2_t + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_t"), start, count, geom->w2_t);
  count[0] = geom->n3_t + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_t"), start, count, geom->w3_t);
  count[0] = geom->n2_e1 + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e1"), start, count, geom->w2_e1);
  count[0] = geom->n3_e1 + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e1"), start, count, geom->w3_e1);
  count[0] = geom->n2_e2 + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e2"), start, count, geom->w2_e2);
  count[0] = geom->n3_e2 + 1;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e2"), start, count, geom->w3_e2);

  count[0] = geom->szc;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2i"), start, count, geom->s2i);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2j"), start, count, geom->s2j);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2k"), start, count, geom->s2k);
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "c2cc"), start, count, geom->c2cc);

  count[0] = geom->szcS;
  ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "botz"), start, count, geom->botz);

  /* OBC                                                               */
  /*ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nobc"), start, count, &geom->nobc);*/

  if (geom->nobc) {
    d1 = i_alloc_1d(geom->nobc);
    count[0] = geom->nobc;
    for (nb = 0; nb < geom->nobc; nb++) {
      d1[nb] = geom->open[nb]->no2_t;
    }
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "no2_t"), start, count, d1);

    for (nb = 0; nb < geom->nobc; nb++) {
      i = geom->open[nb]->no2_e1 + geom->open[nb]->to2_e1 -
	geom->open[nb]->no3_e1;
      d1[nb] = i;
    }
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "no2_e1"), start, count, d1);
    i_free_1d(d1);
  }

  ncw_close(name, cdfid);
}

/* END dump_geom_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the global geometry from file.                              */
/* Note: many of the dimensions declared are of the array size, e.g. */
/* n2_t dimension is written as n2_t+1, so decrement these values    */
/* after arrays are read in to reflext their actual value.           */
/*-------------------------------------------------------------------*/
void read_geom_us(parameters_t *params, geometry_t *geom, char *name)
{
  int fid, n;
  int c, cc, cs, cg, cb, c1, c2, tn, *d1, ns2;
  int ee, e, i, j, k;
  int *mask, *maske, *imape, *imapc, *omapc;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  size_t d;
  size_t szcid;            /* 3D sparse dimension           */
  size_t szcSid;           /* 2D sparse dimension           */
  size_t szeid;            /* 3D sparse dimension           */
  size_t szeSid;           /* 2D sparse dimension           */
  size_t szvid;            /* 3D sparse dimension           */
  size_t szvSid;           /* 2D sparse dimension           */
  size_t kgridid;          /* K dimension id at grid corner */
  char key[MAXSTRLEN];
  int dims[10];
  int vid, oid, tid;
  long t;
  int nz, sednz;
  mesh_t *m;         /* Input mesh                                   */
  double bmax;
  unsigned long ***flag = NULL; /* Flag for Cartesian grid           */
  double **topo = NULL;         /* Cartesian bathymetry              */
  int laus = 2;      /* Ghost cells adjacent to OBCs                 */

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  m = params->mesh;
  ns2 = m->ns2;
  geom->us_type = params->us_type;
  geom->nwindows = params->nwindows;
  if (params->win_size) {
    geom->win_size = d_alloc_1d(geom->nwindows + 1);
    for (n = 1; n <= geom->nwindows; n++) {
      geom->win_size[n] = params->win_size[n - 1];
    }
  }
  geom->compatible = params->compatible;

  /* Open the dump file for reading                                  */
  ncw_open(name, NC_NOWRITE, &fid);

  /* Get the dimensions                                              */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_grid"), &d);
  geom->nz = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "szc"), &d);
  geom->szc = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "szcS"), &d);
  geom->szcS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "sze"), &d);
  geom->sze = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "szeS"), &d);
  geom->szeS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "szv"), &d);
  geom->szv = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "szvS"), &d);
  geom->szvS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "enon"), &d);
  geom->enon = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "enonS"), &d);
  geom->enonS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "ewet"), &d);
  geom->ewet = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "ewetS"), &d);
  geom->ewetS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "snon"), &d);
  geom->snon = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "snonS"), &d);
  geom->snonS = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v2_t"), &d);
  geom->v2_t = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b2_t"), &d);
  geom->b2_t = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n2_t"), &d);
  geom->n2_t = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v3_t"), &d);
  geom->v3_t = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b3_t"), &d);
  geom->b3_t = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n3_t"), &d);
  geom->n3_t = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v2_e1"), &d);
  geom->v2_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b2_e1"), &d);
  geom->b2_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "x2_e1"), &d);
  geom->x2_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n2_e1"), &d);
  geom->n2_e1 = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v3_e1"), &d);
  geom->v3_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b3_e1"), &d);
  geom->b3_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "x3_e1"), &d);
  geom->x3_e1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n3_e1"), &d);
  geom->n3_e1 = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v2_e2"), &d);
  geom->v2_e2 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b2_e2"), &d);
  geom->b2_e2 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n2_e2"), &d);
  geom->n2_e2 = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "v3_e2"), &d);
  geom->v3_e2 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "b3_e2"), &d);
  geom->b3_e2 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "n3_e2"), &d);
  geom->n3_e2 = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nobc"), &d);
  geom->nobc = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "bdry"), &d);
  geom->bdry = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "npe"), &d);
  geom->npem = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbptS"), &d);
  geom->nbptS = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbpt"), &d);
  geom->nbpt = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbpte1S"), &d);
  geom->nbpte1S = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbpte1"), &d);
  geom->nbpte1 = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbe1"), &d);
  geom->nbe1 = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nbe1S"), &d);
  geom->nbe1S = (long)d;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "neem"), &d);
  geom->neem = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nve"), &d);
  geom->nvem = (long)d-1;
  nc_inq_dimlen(fid, ncw_dim_id(fid, "nvc"), &d);
  geom->nvcm = (long)d-1;
  geom->sgnum = geom->szc - 1;
  geom->sgnumS = geom->szcS - 1;
  geom->sgsiz = geom->szc;
  geom->sgsizS = geom->szcS;

  /* Check                                                           */
  if (geom->npem != m->mnpe)
    hd_quit("Maximum number of edges differs : %d vs %d\n", geom->npem, m->mnpe);
  if (geom->nobc != m->nobc)
    hd_quit("Number of OBCs differs : %d vs %d\n", geom->nobc, m->nobc);

  oid = ncw_dim_id(fid, "one");
  tid = ncw_dim_id(fid, "two");
  nc_get_var_int(fid, ncw_var_id(fid, "sednz"), &sednz);
  geom->sednz = sednz;

  /* Get dimensions                                                  */
  emstag(LINFO,"read_windows","Reading geometry map\n");

  /* CENTRES                                                         */
  /* geom->npe[szcS]                                                 */
  /* geom->c2c[npem][szc]                                            */
  /* geom->c2e[npem][szc]                                            */
  /* geom->c2v[npem][szc]                                            */
  /* geom->vIc[npem][szcS]                                           */
  /* geom->eSc[npem][szcS]                                           */
  /* geom->zm1[szc]                                                  */
  /* geom->zp1[szc]                                                  */
  /* geom->m2d[szc]                                                  */
  /* geom->wsa[szc]                                                  */
  /* geom->wgst[szc]                                                 */
  /* geom->bot_t[szcS]                                               */
  /* geom->sur_t[szcS]                                               */
  /* geom->nsur_t[szcS]                                              */
  /* geom->bpt[nbpt]                                                 */
  /* geom->bin[nbpt]                                                 */
  /* geom->w2_t[szcS]                                                */
  /* geom->w3_t[szc]                                                 */
  alloc_geom_us(geom, CENTRE_A);
  geom->w2_t = i_alloc_1d(geom->szcS);
  geom->w3_t = i_alloc_1d(geom->szc);
  geom->sur_t = i_alloc_1d(geom->b2_t + 1);
  geom->bot_t = i_alloc_1d(geom->b2_t + 1);
  geom->s2i = i_alloc_1d(geom->szc);
  geom->s2j = i_alloc_1d(geom->szc);
  geom->s2k = i_alloc_1d(geom->szc);
  geom->c2cc = i_alloc_1d(geom->szc);
  geom->cc2s = i_alloc_1d(ns2+1);
  geom->bpt = i_alloc_1d(geom->nbpt + 1);
  geom->bin = i_alloc_1d(geom->nbpt + 1);
  geom->dbin = i_alloc_1d(geom->nbpt + 1);
  geom->dbpt = i_alloc_1d(geom->nbpt + 1);
  geom->wgst = i_alloc_1d(geom->enon + 1);
  geom->m2d = i_alloc_1d(geom->szc);
  geom->c2e = i_alloc_2d(geom->szc, geom->npem+1);
  geom->c2v = i_alloc_2d(geom->szc, geom->npem+1);
  geom->vIc = i_alloc_2d(geom->szcS, geom->npem + 1);
  geom->eSc = i_alloc_2d(geom->szcS, geom->npem + 1);
  geom->layers = d_alloc_1d(geom->nz + 1);

  count[0] = geom->nz+1;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "z_grid"), start, count, geom->layers);
  for (cc = 0; cc <= geom->nz; cc++) {
    if (params->layers[cc] != geom->layers[cc]) {
      hd_warn("Layer %d : layer depths = %f (.prm) %f (file)\n", cc,
	     params->layers[cc], geom->layers[cc]);
      hd_quit("Incompatible layer depths in window map file %s with parameter file.\n", name);
    }
  }

  count[0] = geom->szc;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "cellz"), start, count, geom->cellz);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "gridz"), start, count, geom->gridz);

  count[0] = geom->b2_t + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bot_t"), start, count, geom->bot_t);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "sur_t"), start, count, geom->sur_t);

  count[0] = geom->szc;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2d"), start, count, geom->m2d);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "wsa"), start, count, geom->wsa);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1"), start, count, geom->zp1);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1"), start, count, geom->zm1);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "wgst"), start, count, geom->wgst);
  
  count[0] = geom->szcS;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "npe"), start, count, geom->npe);
  count[0] = ns2+1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "cc2s"), start, count, geom->cc2s);
  count[0] = geom->npem + 1;
  count[1] = geom->szc;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2c"), start, count, geom->c2c[0]);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2e"), start, count, geom->c2e[0]);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2v"), start, count, geom->c2v[0]);

  count[0] = geom->nbpt + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpt"), start, count, geom->bpt);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bin"), start, count, geom->bin);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "dbpt"), start, count, geom->dbpt);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "dbin"), start, count, geom->dbin);

  count[0] = geom->n2_t + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_t"), start, count, geom->w2_t);
  count[0] = geom->n3_t + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_t"), start, count, geom->w3_t);
 
  count[0] = geom->szc;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2i"), start, count, geom->s2i);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2j"), start, count, geom->s2j);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "s2k"), start, count, geom->s2k);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "c2cc"), start, count, geom->c2cc);

  count[0] = geom->npem + 1;
  count[1] = geom->szcS;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "vIc"), start, count, geom->vIc[0]);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSc"), start, count, geom->eSc[0]);

  count[0] = geom->szcS;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "cellx"), start, count, geom->cellx);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "celly"), start, count, geom->celly);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "cellarea"), start, count, geom->cellarea);

  /* EDGES                                                           */ 
  /* geom->nee[szeS]                                                 */
  /* geom->wse[sze]                                                  */
  /* geom->e2c[sze][2]                                               */
  /* geom->e2v[sze][2]                                               */
  /* geom->e2e[sze][2]                                               */
  /* geom->eSe[neem][sze]                                            */
  /* geom->wAe[neem][sze]                                            */
  /* geom->ep[sze]                                                   */
  /* geom->em[sze]                                                   */
  /* geom->zm1e[sze]                                                 */
  /* geom->zp1e[sze]                                                 */
  /* geom->e2k[sze]                                                  */
  /* geom->m2de[sze]                                                 */ 
  /* geom->bot_e1[szeS]                                              */
  /* geom->sur_e1[szeS]                                              */
  /* geom->bpte1[nbpte1]                                             */
  /* geom->bine1[nbpte1]                                             */
  /* geom->bine2[nbpte1]                                             */
  /* geom->w2_e1[szeS]                                               */
  /* geom->w3_e1[sze]                                                */
  alloc_geom_us(geom, EDGE_A);
  geom->e2c = i_alloc_2d(2, geom->sze);
  geom->e2v = i_alloc_2d(2, geom->sze);
  geom->e2e = i_alloc_2d(2, geom->sze);
  geom->eSe = i_alloc_2d(geom->sze, geom->neem + 1);
  geom->wAe = d_alloc_2d(geom->sze, geom->neem + 1);
  geom->ep = i_alloc_1d(geom->sze);
  geom->em = i_alloc_1d(geom->sze);
  geom->zm1e = i_alloc_1d(geom->sze);
  geom->zp1e = i_alloc_1d(geom->sze);
  geom->e2k = i_alloc_1d(geom->sze);
  geom->e2ijk = i_alloc_1d(geom->sze);
  geom->m2de = i_alloc_1d(geom->sze);
  geom->sur_e1 = i_alloc_1d(geom->szeS);
  geom->bot_e1 = i_alloc_1d(geom->szeS);
  geom->bpte1 = i_alloc_1d(geom->nbpte1 + 1);
  geom->bine1 = i_alloc_1d(geom->nbpte1 + 1);
  geom->bine2 = i_alloc_1d(geom->nbpte1 + 1);
  geom->bpte1S = i_alloc_1d(geom->nbpte1S + 1);
  geom->bine1S = i_alloc_1d(geom->nbpte1S + 1);
  geom->w2_e1 = i_alloc_1d(geom->szeS);
  geom->w3_e1 = i_alloc_1d(geom->sze);

  count[0] = geom->szeS;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bot_e1"), start, count, geom->bot_e1);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "sur_e1"), start, count, geom->sur_e1);
   
  count[0] = geom->sze;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2de"), start, count, geom->m2de);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "wse"), start, count, geom->wse);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1e"), start, count, geom->zp1e);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1e"), start, count, geom->zm1e);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "ep"), start, count, geom->ep);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "em"), start, count, geom->em);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2k"), start, count, geom->e2k);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2ijk"), start, count, geom->e2ijk);
  if (geom->nobc) {
    maske = i_alloc_1d(geom->sze);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "maske"), start, count, maske);
  }

  count[0] = geom->sze;
  count[1] = 2;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2c"), start, count, geom->e2c[0]);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2e"), start, count, geom->e2e[0]);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "e2v"), start, count, geom->e2v[0]);
  
  count[0] = geom->neem + 1;
  count[1] = geom->sze;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSe"), start, count, geom->eSe[0]);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "wAe"), start, count, geom->wAe[0]);

  count[0] = geom->szeS;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "nee"), start, count, geom->nee);

  count[0] = geom->nbpte1 + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpte1"), start, count, geom->bpte1);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bine1"), start, count, geom->bine1);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bine2"), start, count, geom->bine2);

  count[0] = geom->nbpte1S + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bpte1S"), start, count, geom->bpte1S);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "bine1S"), start, count, geom->bine1S);

  count[0] = geom->n2_e1 + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_e1"), start, count, geom->w2_e1);
  count[0] = geom->n3_e1 + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_e1"), start, count, geom->w3_e1);

  count[0] = geom->szeS;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "u1x"), start, count, geom->u1x);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "u1y"), start, count, geom->u1y);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "h1au1"), start, count, geom->h1au1);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "h2au1"), start, count, geom->h2au1);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "h1acell"), start, count, geom->h1acell);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "edgearea"), start, count, geom->edgearea);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "thetau1"), start, count, geom->thetau1);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "thetau2"), start, count, geom->thetau2);
  count[0] = geom->npem + 1;
  count[1] = geom->n2_t + 1;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "hacell"), start, count, geom->hacell[0]);

  /* VERTICES                                                        */ 
  /* geom->nve[szvS]                                                 */
  /* geom->nvc[szvS]                                                 */
  /* geom->wsv[szv]                                                  */
  /* geom->v2c[szv][nvcm]                                            */
  /* geom->v2e[szv][nvem]                                            */
  /* geom->eSv[nvem][szv]                                            */
  /* geom->zm1v[szv]                                                 */
  /* geom->zp1v[szv]                                                 */
  /* geom->dualarea[szvS]                                            */
  /* geom->dualareap[szvS][nvcm]                                     */
  /* geom->m2dv[szv]                                                 */
  /* geom->bot_e2[szvS]                                              */
  /* geom->sur_e2[szvS]                                              */
  /* geom->w2_e2[szvS]                                               */
  /* geom->w3_e2[szv]                                                */
  alloc_geom_us(geom, VERTEX_A);
  geom->w2_e2 = i_alloc_1d(geom->szvS);
  geom->w3_e2 = i_alloc_1d(geom->szv);
  geom->m2dv = i_alloc_1d(geom->szv);
  geom->zm1v = i_alloc_1d(geom->szv);
  geom->zp1v = i_alloc_1d(geom->szv);
  if (params->us_type & US_IJ) geom->v2ijk = i_alloc_1d(geom->szv);
  geom->v2c = i_alloc_2d(geom->nvcm+1, geom->szv);
  geom->v2e = i_alloc_2d(geom->nvem+1, geom->szv);
  geom->eSv = i_alloc_2d(geom->szv, geom->nvem + 1);

  count[0] = geom->szv;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "m2dv"), start, count, geom->m2dv);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "wsv"), start, count, geom->wsv);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zp1v"), start, count, geom->zp1v);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "zm1v"), start, count, geom->zm1v);

  if (params->us_type & US_IJ)
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2ijk"), start, count, geom->v2ijk);

  count[0] = geom->szv;
  count[1] = geom->nvcm + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2c"), start, count, geom->v2c[0]);
  count[1] = geom->nvem + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "v2e"), start, count, geom->v2e[0]);

  count[0] = geom->szvS;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "nve"), start, count, geom->nve);
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "nvc"), start, count, geom->nvc);

  count[0] = geom->nvem + 1;
  count[1] = geom->szv;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "eSv"), start, count, geom->eSv[0]);

  count[0] = geom->szvS;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "dualarea"), start, count, geom->dualarea);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "gridx"), start, count, geom->gridx);
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "gridy"), start, count, geom->gridy);

  count[0] = geom->szvS;
  count[1] = geom->nvcm + 1;
  ncw_get_vara_double(name, fid, ncw_var_id(fid, "dualareap"), start, count,
		      geom->dualareap[0]);

  count[0] = geom->n2_e2 + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w2_e2"), start, count, geom->w2_e2);
  count[0] = geom->n3_e2 + 1;
  ncw_get_vara_int(name, fid, ncw_var_id(fid, "w3_e2"), start, count, geom->w3_e2);

  geom->botz = d_alloc_1d(geom->szcS);
  count[0] = geom->szcS;
  nc_get_vara_double(fid, ncw_var_id(fid, "botz"), start, count, geom->botz);

  /* Max sizes                                                       */
  geom->szm = max(max(geom->szc, geom->sze), geom->szv);
  geom->szmS = max(max(geom->szcS, geom->szeS), geom->szvS);

  emstag(LINFO,"read_geom_us","Global geometry map OK\n", n);

  ncw_close(name, fid);

  if (params->us_type & US_IJ) {
    geom->nce1 = m->nce1;
    geom->nce2 = m->nce2;
    geom->nfe1 = m->nce1 + 1;
    geom->nfe2 = m->nce2 + 1;
    geom->map = (unsigned long ***)l_alloc_3d(m->nce1+1, m->nce2+1, nz);
    flag = (unsigned long ***)l_alloc_3d(geom->nce1 + 1, geom->nce2 + 1, nz);
    topo = d_alloc_2d(geom->nce1, geom->nce2);

    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      geom->map[geom->s2k[c]][geom->s2j[c]][geom->s2i[c]] = c;
    }
    /* Initialize                                                      */
    for (j = 0; j < geom->nce2; j++) {
      for (i = 0; i < geom->nce1; i++) {
	if (params->topo)
	  topo[j][i] = params->topo[j][i];
	else
	  topo[j][i] = NaN;
	for (k = geom->nz - 1; k >= 0; k--) {
	  if (params->flag)
	    flag[k][j][i] = params->flag[j][i];
	  else
	    flag[k][j][i] = SOLID;
	}
      }
    }
    for (cc = 1; cc <= ns2; cc++) {
      c = geom->cc2s[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      topo[j][i] = params->bathy[cc];
    }
    /* Ghost cells */
    for (cc = 1; cc <= geom->n2_t; cc++) {
      int c1;
      c = geom->w2_t[cc];
      c1 = geom->wgst[c];
      if (c1) {
	geom->map[geom->s2k[c]][geom->s2j[c]][geom->s2i[c]] = c;
      }
    }
    /* Vertex maps to (i,j)                                           */
    if (geom->npem == 4)
      vertex_map_4(geom);
  }

  /*-----------------------------------------------------------------*/
  /* OBCs                                                            */
  get_mesh_obc(params, params->mesh->neic);
  make_geom_obc(params, geom, laus, maske);

  if (params->map_type & GEOM_CHECK) return;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      e = open->obc_e1[ee];
      /* OBC lengths currently hold the distance from the centre to  */
      /* the edge. In master_setghosts() we assume that the ghost    */
      /* centre is the same distance from the edge as this wet       */
      /* centre; i.e. we multiply by 2. Using the -g option, which   */
      /* is used to dump a geom file, master_setghosts() is called   */
      /* twice, so here we divide by 4 so that the original length   */
      /* is sent to window_build() below.                            */
      /* This length is only used for a default flux adjustment.     */
      geom->h2au1[e] /= 4.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the sponge zones                                            */
  set_sponge_cells(geom);

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the geometry vectors                        */
  alloc_geom_us(geom, (MAP_A | GRID_A | MASTER_A));

  /*-----------------------------------------------------------------*/
  /* Find the deepest coordinte in the grid                          */
  bmax = 1e10;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->botz[c] < bmax) {
      bmax = geom->botz[c];
      geom->cdeep = geom->bot_t[cc];
    }
  }
  if (params->sigma) {
    for (n = 0; n < params->nz; n++)
      params->layers[n] /= bmax;
  }

  /* Set a no-gradient condition over the sediment for gridz */
  geom->topgrid = geom->layers[geom->nz];

  if (DEBUG("init_m"))
    dlog("init_m", "\nMaster geometry created OK\n");

  /*-----------------------------------------------------------------*/
  /* Set up the windows                                              */
  TIMING_SET;
  window_build(geom, params);
  TIMING_DUMP(1, "  window_build");
  
  /*-----------------------------------------------------------------*/
  /* Get the centre and edge angles                                  */
  write_us_map(geom, params);

  /*-----------------------------------------------------------------*/
  /* Make and initialise the master data structure                   */
  master = master_build(params, geom);
  master->geom = geom;
  master->sgrid = geom;
  get_filter(geom);
  if (DEBUG("init_m"))
    dlog("init_m", "\nMaster data created OK\n");

  /*-----------------------------------------------------------------*/
  /* Set up the dump data structure                                  */
  dumpdata = dumpdata_build(params, geom, master, flag, params->layers, topo);
  if (DEBUG("init_m"))
    dlog("init_m", "\nDumpdata created OK\n");

  /*-----------------------------------------------------------------*/
  /* Make a map from Delaunay triangulation to unstructured          */
  /* coordinate. Used for mapping (x,y) to i.                        */
  if (!(params->us_type & US_IJ))
    create_delaunay_cell(geom, params);

  /* Set geographical flag */
  geom->is_geog = (strlen(params->projection) > 0) &&
    (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);

  if (flag) l_free_3d((long ***)flag);
  if (topo) d_free_2d(topo);
  if(maske) i_free_1d(maske);
}

/* END read_geom_us()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks a window map against the maps generated inline.            */
/*-------------------------------------------------------------------*/
void check_geom_map_us(parameters_t *params, geometry_t *geom, char *name)
{
  geometry_t *sgrid;
  int e, ee, n, nn, j, c, c1, c2, cc, vv, v, wn;
  int laus = 2;      /* Ghost cells adjacent to OBCs                 */

  sgrid = (geometry_t *)malloc(sizeof(geometry_t));
  memset(sgrid, 0, sizeof(geometry_t));
  read_geom_us(params, sgrid, name);

  printf("Check geom for map file %s\n", name);
  if (geom->npem != sgrid->npem)
    hd_quit("npem error %d:%d\n", geom->npem, sgrid->npem);
  if (geom->neem != sgrid->neem)
    hd_quit("neem error %d:%d\n", geom->neem, sgrid->neem);
  if (geom->nvem != sgrid->nvem)
    hd_quit("nvem error %d:%d\n", geom->nvem, sgrid->nvem);
  if (geom->nvcm != sgrid->nvcm)
    hd_quit("nvcm error %d:%d\n", geom->nvcm, sgrid->nvcm);
  
  /* CENTRES                                                         */
  /* window->npe[szcS]                                               */
  /* window->c2c[npem][szc]                                          */
  /* window->c2e[npem][szc]                                          */
  /* window->c2v[npem][szc]                                          */
  /* window->vIc[npem][szcS]                                         */
  /* window->eSc[npem][szcS]                                         */
  /* window->zm1[szc]                                                */
  /* window->zp1[szc]                                                */
  /* window->m2d[szc]                                                */
  /* window->wsa[szc]                                                */
  /* window->wgst[szc]                                               */
  /* window->bot_t[szcS]                                             */
  /* window->sur_t[szcS]                                             */
  /* window->bpt[nbpt]                                               */
  /* window->bin[nbpt]                                               */
  /* window->dbpt[nbpt]                                              */
  /* window->dbin[nbpt]                                              */
  /* window->w2_t[szcS]                                              */
  /* window->w3_t[szc]                                               */
  if (geom->v2_t != sgrid->v2_t)
    hd_quit("v2_t %d:%d\n", geom->v2_t, sgrid->v2_t);
  if (geom->b2_t != sgrid->b2_t)
    hd_quit("b2_t %d:%d\n", geom->b2_t, sgrid->b2_t);
  if (geom->n2_t != sgrid->n2_t)
    hd_quit("n2_t %d:%d\n", geom->n2_t, sgrid->n2_t);
  for (cc = 1; cc <= geom->n2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->w2_t[cc] != sgrid->w2_t[cc]) 
      check_print(geom, "w2_t", cc, c, geom->w2_t[cc], sgrid->w2_t[cc], C2D);
    if (geom->npe[c] != sgrid->npe[c]) 
      check_print(geom, "npe", cc, c, geom->npe[c], sgrid->npe[c], C2D);
    if (geom->cellarea[c] != sgrid->cellarea[c]) 
      check_print(geom, "cellarea", cc, c, geom->cellarea[c], sgrid->cellarea[c], C2D);
    for (nn = 1; nn <= geom->npem; nn++) {
      if (geom->vIc[nn][c] != sgrid->vIc[nn][c]) 
	check_print(geom, "vIc", cc, c, geom->vIc[nn][c], sgrid->vIc[nn][c], C2D);
      if (geom->eSc[nn][c] != sgrid->eSc[nn][c]) 
	check_print(geom, "eSc", cc, c, geom->eSc[nn][c], sgrid->eSc[nn][c], C2D);
    }
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->bot_t[cc] != sgrid->bot_t[cc]) 
      check_print(geom, "bot_t", cc, c, geom->bot_t[cc], sgrid->bot_t[cc], C2D);
    if (geom->sur_t[cc] != sgrid->sur_t[cc]) 
      check_print(geom, "sur_t", cc, c, geom->sur_t[cc], sgrid->sur_t[cc], C2D);
  }
  
  if (geom->v3_t != sgrid->v3_t)
    hd_quit("v3_t %d:%d\n", geom->v3_t, sgrid->v3_t);
  if (geom->b3_t != sgrid->b3_t)
    hd_quit("b3_t %d:%d\n", geom->b3_t, sgrid->b3_t);
  if (geom->n3_t != sgrid->n3_t)
    hd_quit("n3_t %d:%d\n", geom->n3_t, sgrid->n3_t);
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    if (geom->w3_t[cc] != sgrid->w3_t[cc])
      check_print(geom, "w3_t", cc, c, geom->w3_t[cc], sgrid->w3_t[cc], C3D);
    if (geom->zm1[c] != sgrid->zm1[c])
      check_print(geom, "zm1", cc, c, geom->zm1[c], sgrid->zm1[c], C3D);
    if (geom->zp1[c] != sgrid->zp1[c])
      check_print(geom, "zp1", cc, c, geom->zp1[c], sgrid->zp1[c], C3D);
    if (geom->m2d[c] != sgrid->m2d[c]) 
      check_print(geom, "m2d", cc, c, geom->m2d[c], sgrid->m2d[c], C3D);
    if (geom->wsa[c] != sgrid->wsa[c]) 
      check_print(geom, "wsa", cc, c, geom->wsa[c], sgrid->wsa[c], C3D);
    if (geom->wgst[c] != sgrid->wgst[c]) 
      check_print(geom, "wgst", cc, c, geom->wgst[c], sgrid->wgst[c], C3D);
    if (geom->s2k[c] != sgrid->s2k[c]) 
      check_print(geom, "s2k", cc, c, geom->s2k[c], sgrid->s2k[c], C3D);
    for (nn = 1; nn <= geom->npem; nn++) {
      if (geom->c2c[nn][c] != sgrid->c2c[nn][c]) hd_quit ("c2c error @ c = %d\n", c);
      if (geom->c2e[nn][c] != sgrid->c2e[nn][c]) hd_quit ("c2e error @ c = %d\n", c);
      if (geom->c2v[nn][c] != sgrid->c2v[nn][c]) hd_quit ("c2v error @ c = %d\n", c);
    }
  }
  if (geom->nbpt != sgrid->nbpt)
    hd_quit("nbpt %d:%d\n", geom->nbpt, sgrid->nbpt);
  for (cc = 1; cc <= geom->nbpt; cc++) {
    if (geom->bpt[cc] != sgrid->bpt[cc]) hd_quit ("bpt error @ c = %d\n", cc);
    if (geom->bin[cc] != sgrid->bin[cc]) hd_quit ("bin error @ c = %d\n", cc);
    if (geom->dbpt[cc] != sgrid->dbpt[cc]) hd_quit ("dbpt error @ c = %d\n", cc);
    if (geom->dbin[cc] != sgrid->dbin[cc]) hd_quit ("dbin error @ c = %d\n", cc);
  }

  /* EDGES                                                           */ 
  /* window->nee[szeS]                                               */
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
  if (geom->v2_e1 != sgrid->v2_e1)
    hd_quit("v2_e1 %d:%d\n", geom->v2_e1, sgrid->v2_e1);
  if (geom->b2_e1 != sgrid->b2_e1)
    hd_quit("b2_e1 %d:%d\n", geom->b2_e1, sgrid->b2_e1);
  if (geom->n2_e1 != sgrid->n2_e1)
    hd_quit("n2_e1 %d:%d\n", geom->n2_e1, sgrid->n2_e1);
  for (ee = 1; ee <= geom->n2_e1; ee++) {
    e = geom->w2_e1[ee];
    if (geom->w2_e1[ee] != sgrid->w2_e1[ee]) 
      check_print(geom, "w2_e1", ee, e, geom->w2_e1[ee], sgrid->w2_e1[ee], E2D);
    if (geom->bot_e1[ee] != sgrid->bot_e1[ee])
      check_print(geom, "bot_e1", ee, e, geom->bot_e1[ee], sgrid->bot_e1[ee], E2D);
    if (geom->sur_e1[ee] != sgrid->sur_e1[ee])
      check_print(geom, "sur_e1", ee, e, geom->sur_e1[ee], sgrid->sur_e1[ee], E2D);
    if (geom->nee[e] != sgrid->nee[e]) 
      check_print(geom, "nee", ee, e, geom->nee[e], sgrid->nee[e], E2D);
  }
  if (geom->v3_e1 != sgrid->v3_e1)
    hd_quit("v3_e1 %d:%d\n", geom->v3_e1, sgrid->v3_e1);
  if (geom->b3_e1 != sgrid->b3_e1)
    hd_quit("b3_e1 %d:%d\n", geom->b3_e1, sgrid->b3_e1);
  if (geom->n3_e1 != sgrid->n3_e1)
    hd_quit("n3_e1 %d:%d\n", geom->n3_e1, sgrid->n3_e1);
  for (ee = 1; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    if (geom->w3_e1[ee] != sgrid->w3_e1[ee]) 
      check_print(geom, "w3_e1", ee, e, geom->w3_e1[ee], sgrid->w3_e1[ee], E3D);
    if (geom->wse[e] != sgrid->wse[e]) hd_quit("wse error @ %d\n", e);
    if (geom->ep[e] != sgrid->ep[e]) hd_quit("ep error @ %d\n", e);
    if (geom->em[e] != sgrid->em[e]) hd_quit("em error @ %d\n", e);
    if (geom->zm1e[e] != sgrid->zm1e[e]) hd_quit("zm1e error @ %d\n", e);
    if (geom->zp1e[e] != sgrid->zp1e[e]) hd_quit("zp1e error @ %d\n", e);
    if (geom->m2de[e] != sgrid->m2de[e]) hd_quit("m2de error @ %d\n", e);
    if (geom->e2k[e] != sgrid->e2k[e]) hd_quit("e2k error @ %d\n", e);
    if (geom->e2c[e][0] != sgrid->e2c[e][0]) hd_quit("e2c[0] error @ %d\n", e);
    if (geom->e2c[e][1] != sgrid->e2c[e][1]) hd_quit("e2c[1] error @ %d\n", e);
    if (geom->e2e[e][0] != sgrid->e2e[e][0]) hd_quit("e2e[0] error @ %d\n", e);
    if (geom->e2e[e][1] != sgrid->e2e[e][1]) hd_quit("e2e[1] error @ %d\n", e);
    if (geom->e2v[e][0] != sgrid->e2v[e][0]) hd_quit("e2v[0] error @ %d\n", e);
    if (geom->e2v[e][1] != sgrid->e2v[e][1]) hd_quit("e2v[1] error @ %d\n", e);
    for (nn = 1; nn <= geom->neem; nn++) {
      if (geom->eSe[nn][e] != sgrid->eSe[nn][e]) hd_quit ("eSe error @ e = %d\n", e);
      if (geom->wAe[nn][e] != sgrid->wAe[nn][e]) hd_quit ("wAe error @ e = %d\n", e);
    }
  }
  if (geom->nbpte1 != sgrid->nbpte1)
    hd_quit("nbpte1 %d:%d\n", geom->nbpte1, sgrid->nbpte1);
  for (ee = 1; ee <= geom->nbpt; ee++) {
    if (geom->bpte1[ee] != sgrid->bpte1[ee]) hd_quit ("bpte1 error @ e = %d\n", ee);
    if (geom->bine1[ee] != sgrid->bine1[ee]) hd_quit ("bine1 error @ e = %d\n", ee);
  }

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
  /* window->bot_e2[szvS]                                            */
  /* window->sur_e2[szvS]                                            */
  /* window->w2_e2[szvS]                                             */
  /* window->w3_e2[szv]                                              */
  if (geom->v2_e2 != sgrid->v2_e2)
    hd_quit("v2_e2 %d:%d\n", geom->v2_e2, sgrid->v2_e2);
  if (geom->b2_e2 != sgrid->b2_e2)
    hd_quit("b2_e2 %d:%d\n", geom->b2_e2, sgrid->b2_e2);
  if (geom->n2_e2 != sgrid->n2_e2)
    hd_quit("n2_e2 %d:%d\n", geom->n2_e2, sgrid->n2_e2);
  for (vv = 1; vv <= geom->n2_e2; vv++) {
    vv = geom->w2_e2[vv];
    if (geom->w2_e2[vv] != sgrid->w2_e2[vv]) 
      check_print(geom, "w2_e2", vv, v, geom->w2_e2[vv], sgrid->w2_e2[vv], V2D);
    /*
      if (geom->bot_e2[vv] != sgrid->bot_e2[vv]) 
      check_print(geom, "bot_e2", vv, v, geom->bot_e2[vv], sgrid->bot_e2[vv], V2D;
      if (geom->sur_e2[vv] != sgrid->sur_e2[vv]) 
      check_print(geom, "sur_e2", vv, v, geom->sur_e2[vv], sgrid->sur_e2[vv], V2D;
    */
    if (geom->nvc[v] != sgrid->nvc[v]) hd_quit ("nvc error @ v = %d\n", v);
    if (geom->nve[v] != sgrid->nve[v]) hd_quit ("nve error @ v = %d\n", v);
    if (geom->dualarea[v] != sgrid->dualarea[v]) hd_quit ("dualarea error @ v = %d\n", v);
  }
  if (geom->v3_e2 != sgrid->v3_e2)
    hd_quit("v3_e2 %d:%d\n", geom->v3_e2, sgrid->v3_e2);
  if (geom->b3_e2 != sgrid->b3_e2)
    hd_quit("b3_e2 %d:%d\n", geom->b3_e2, sgrid->b3_e2);
  if (geom->n3_e2 != sgrid->n3_e2)
    hd_quit("n3_e2 %d:%d\n", geom->n3_e2, sgrid->n3_e2);
  for (vv = 1; vv <= geom->n3_e2; vv++) {
    vv = geom->w3_e2[vv];
    if (geom->w3_e2[vv] != sgrid->w3_e2[vv]) 
      check_print(geom, "w3_e2", vv, v, geom->w3_e2[vv], sgrid->w3_e2[vv], V3D);
    if (geom->wsv[v] != sgrid->wsv[v]) hd_quit ("wsv error @ v = %d\n", v);
    if (geom->zm1v[v] != sgrid->zm1v[v]) hd_quit ("zm1v error @ v = %d\n", v);
    if (geom->zp1v[v] != sgrid->zp1v[v]) hd_quit ("zp1v error @ v = %d\n", v);
    if (geom->m2dv[v] != sgrid->m2dv[v]) hd_quit ("m2dv error @ v = %d\n", v);
    for (nn = 1; nn <= geom->nvcm; nn++) {
      if (geom->v2c[v][nn] != sgrid->v2c[v][nn]) hd_quit ("v2c error @ v = %d\n", v);
      if (geom->dualareap[v][nn] != sgrid->dualareap[v][nn]) hd_quit ("dualareap error @ v = %d\n", v);
    }
    for (nn = 1; nn <= geom->nvem; nn++) {
      if (geom->v2e[v][nn] != sgrid->v2e[v][nn]) hd_quit ("v2e error @ v = %d\n", v);
      if (geom->eSv[nn][v] != sgrid->eSv[nn][v]) hd_quit ("eSv error @ v = %d\n", v);
    }
  }

  printf("Checking boundaries\n");
  if (geom->nobc != sgrid->nobc)
    check_print_obc("nobc", n, 0, geom->nobc, sgrid->nobc);
  for (n = 0; n < geom->nobc; n++) {
    int cs;
    open_bdrys_t *open = geom->open[n];
    open_bdrys_t *openw = sgrid->open[n];
    if (open->bgz != openw->bgz)
      check_print_obc("bgz", n, 0, open->bgz, openw->bgz);
    if(open->no2_t != openw->no2_t)
      check_print_obc("no2_t", 0, n, open->no2_t, openw->no2_t);
    if(open->no2_e1 != openw->no2_e1)
      check_print_obc("no2_e1", n, 0, open->no2_e1, openw->no2_e1);
    if(open->no2_e2 != openw->no2_e2)
      check_print_obc("no2_e2", n, 0, open->no2_e2, openw->no2_e2);
    if(open->no3_t != openw->no3_t)
      check_print_obc("no3_t", n, 0, open->no3_t, openw->no3_t);
    for (cc = 1; cc <= open->no3_t; cc++) {
      if (open->obc_t[cc] != openw->obc_t[cc]) 
	check_print_obc("obc_t", n, cc, open->obc_t[cc], openw->obc_t[cc]);
      if (open->oi1_t[cc] != openw->oi1_t[cc]) 
	check_print_obc("oi1_t", n, cc, open->oi1_t[cc], openw->oi1_t[cc]);
      if (open->oi2_t[cc] != openw->oi2_t[cc])
	check_print_obc("oi2_t", n, cc, open->oi2_t[cc], openw->oi2_t[cc]);
      if (open->nepc[cc] != openw->nepc[cc])
	check_print_obc("nepc", n, cc, open->nepc[cc], openw->nepc[cc]);
      if (open->inc[cc] != openw->inc[cc])
	check_print_obc("inc", n, cc, open->inc[cc], openw->inc[cc]);
      if (open->olap[cc] != openw->olap[cc])
	check_print_obc("olap", n, cc, open->olap[cc], openw->olap[cc]);
      if (cc <= open->no2_t && open->bot_t[cc] != openw->bot_t[cc])
	check_print_obc("bot_t", n, cc, open->bot_t[cc], openw->bot_t[cc]);
      c = open->obc_t[cc];
      cs = geom->m2d[c];
      for (j = 1; j <= geom->npe[cs]; j++) {
	if (open->bec[j][cc] != openw->bec[j][cc])
	  check_print_obc("bec", n, cc, open->bec[j][cc], openw->bec[j][cc]);
	if (open->bcc[j][cc] != openw->bcc[j][cc])
	  check_print_obc("bcc", n, cc, open->bcc[j][cc], openw->bcc[j][cc]);
      }
    }
    if(open->no2_e1 != openw->no2_e1)
      check_print_obc("no2_e1", n, 0, open->no2_e1, openw->no2_e1);
    if(open->to2_e1 != openw->to2_e1)
      check_print_obc("to2_e1", n, 0, open->to2_e1, openw->to2_e1);
    if(open->no3_e1 != openw->no3_e1)
      check_print_obc("no3_e1", n, 0, open->no3_e1, openw->no3_e1);
    if(open->to3_e1 != openw->to3_e1)
      check_print_obc("to3_e1", n, 0, open->to3_e1, openw->to3_e1);
    for (ee = 1; ee <= open->to3_e1; ee++) {
      if (open->obc_e1[ee] != openw->obc_e1[ee]) 
	check_print_obc("obc_e1", n, ee, open->obc_e1[ee], openw->obc_e1[ee]);
      if (open->obc_e2[ee] != openw->obc_e2[ee]) 
	check_print_obc("obc_e2", n, ee, open->obc_e2[ee], openw->obc_e2[ee]);
      if (open->oi1_e1[ee] != openw->oi1_e1[ee]) 
	check_print_obc("oi1_e1", n, ee, open->oi1_e1[ee], openw->oi1_e1[ee]);
      if (open->oi2_e1[ee] != openw->oi2_e1[ee]) 
	check_print_obc("oi2_e1", n, ee, open->oi2_e1[ee], openw->oi2_e1[ee]);
      if (ee <= open->no3_e1 && open->ogc_t[ee] != openw->ogc_t[ee]) 
	check_print_obc("ogc_t", n, ee, open->ogc_t[ee], openw->ogc_t[ee]);
      if (ee <= open->no3_e1 && open->ceni[ee] != openw->ceni[ee]) 
	check_print_obc("ceni", n, ee, open->ceni[ee], openw->ceni[ee]);
      if (ee <= open->no3_e1 && open->ini[ee] != openw->ini[ee]) 
	check_print_obc("ini", n, ee, open->ini[ee], openw->ini[ee]);
      if (ee <= open->no3_e1 && open->outi[ee] != openw->outi[ee]) 
	check_print_obc("outi", n, ee, open->outi[ee], openw->outi[ee]);
      if (ee <= open->no3_e1 && open->dir[ee] != openw->dir[ee]) 
	check_print_obc("dir", n, ee, open->dir[ee], openw->dir[ee]);
      }
  }
  printf("Done check_geom\n");
}

/* END check_geom_map_us()                                           */
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

void check_print(geometry_t *window, char *tag, int cc, int c, int v1, int v2, int mode)
{
  int vn, bn, an, nn, dn;
  char index[MAXSTRLEN], pos[MAXSTRLEN];

  if(mode & (C2D|C3D)) {
    strcpy(index, "cc");
    strcpy(pos, "c");
  }
  if(mode & (E2D|E3D)) {
    strcpy(index, "ee");
    strcpy(pos, "e");
  }
  if(mode & (V2D|V3D)) {
    strcpy(index, "vv");
    strcpy(pos, "v");
  }
  dn = 2;
  vn = window->v2_t;
  bn = window->b2_t;
  an = window->a2_t;
  nn = window->n2_t;
  if (mode & (C3D|E3D|V3D)) {
    dn = 3;
    vn = window->v3_t;
    bn = window->b3_t;
    an = window->a3_t;
    nn = window->n3_t;
  }
  printf(" v%d_t = %d\n", dn, vn);
  printf(" b%d_t = %d\n", dn, bn);
  printf(" a%d_t = %d\n", dn, an);
  printf(" n%d_t = %d\n", dn, nn);
  hd_quit ("%s error @ %s=%d %s=%d (%d & %d)\n", tag, index, cc, pos, c, v1, v2);
}

void check_print_obc(char *tag, int nb, int cc, int v1, int v2)
{
  if (cc !=0)
    hd_quit("%s error bdry%d @ %d : %d & %d\n", tag, nb, cc, v1, v2);
  else
    hd_quit("%s error bdry%d : %d & %d\n", tag, nb, v1, v2);
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
