/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/df_sparse.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  Output is written in sparse coordinates.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: df_sparse.c 6242 2019-06-12 06:24:07Z riz008 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

// from dumpfile.c
int get_nc_mode(dump_file_t *df);

#define SP_ALL_VARS "u1av u2av wtop topz eta wind1 wind2 patm u1 u2 w dens dens_0 Kz Vz Cd u1vh u2vh "

/* Structure to describe each dump file time dep variable */
typedef struct {
  void **v;                     /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  nc_type type;                 /* netCDF type of this variable */
  int xylocation;
} df_sp_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_sp_var_t *vars;           /* List of dump file variables */
} df_sp_data_t;


void write_dump_attributes_sp(dump_data_t *dumpdata, int cdfid,
			      nc_type fptype, int ilower, int nce1, 
			      int nfe1, int jlower, int nce2, int nfe2, 
			      int ewt, int ns2, char *modulo);
static void df_sp_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static int df_sp_get_varinfo(dump_data_t *dumpdata, dump_file_t *df, 
			     char *name, df_sp_var_t *var);
static void df_sp_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_sp_get_dimids(int cdfid, df_sp_var_t var, int *sz);
static void df_sp_get_dimsizes(dump_file_t *df, df_sp_var_t var, size_t *ne1);
int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref);
void check_window_map(geometry_t **window, char *name);


void write_dump_attributes_sp(dump_data_t *dumpdata, int cdfid,
			      nc_type fptype, int ilower, int nce1, 
			      int nfe1, int jlower, int nce2, int nfe2, 
			      int ns2, int ns3, char *modulo)
{
  /* dimension ids */
  int vid;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double fill;

  /* time independent variables */
  vid = ncw_var_id(cdfid, "z_grid");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer faces");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "z_centre");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at grid layer centre");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  vid = ncw_var_id(cdfid, "x_grid");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Longitude at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "X coordinate at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1,
                     jlower, nce2, CL_GRID);

  vid = ncw_var_id(cdfid, "y_grid");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Latitude at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name",
                   "Y coordinate at grid corners");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_GRID);

  vid = ncw_var_id(cdfid, "x_centre");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Longitude at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
    write_text_att(cdfid, vid, "units", "degrees_east");
  } else {
    write_text_att(cdfid, vid, "long_name", "X coordinate at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "X");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_CENTRE);

  vid = ncw_var_id(cdfid, "y_centre");
  if (is_geog) {
    write_text_att(cdfid, vid, "long_name", "Latitude at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
    write_text_att(cdfid, vid, "units", "degrees_north");
  } else {
    write_text_att(cdfid, vid, "long_name", "Y coordinate at cell centre");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, ilower, nce1, jlower, nce2,
                     CL_CENTRE);

  vid = ncw_var_id(cdfid, "botz");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Z coordinate at sea-bed at cell centre");

  vid = ncw_var_id(cdfid, "h1au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell width at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "h2au1");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell height at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "h1au2");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell width at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "h2au2");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name",
                 "Cell height at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "h1acell");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell width at cell centre");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "h2acell");
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Cell height at cell centre");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "thetau1");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation at centre of left face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "thetau2");
  write_text_att(cdfid, vid, "units", "radian");
  write_text_att(cdfid, vid, "long_name",
                 "Cell rotation at centre of back face");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "s2i");
  write_text_att(cdfid, vid, "units", " ");
  write_text_att(cdfid, vid, "long_name", "Sparse to i coordinate map");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "s2j");
  write_text_att(cdfid, vid, "units", " ");
  write_text_att(cdfid, vid, "long_name", "Sparse to j coordinate map");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "s2k");
  write_text_att(cdfid, vid, "units", " ");
  write_text_att(cdfid, vid, "long_name", "Sparse to k coordinate map");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

  vid = ncw_var_id(cdfid, "coriolis");
  write_text_att(cdfid, vid, "units", " ");
  write_text_att(cdfid, vid, "long_name", "Coriolis parameter");
  write_text_att(cdfid, vid, "coordinates", "sparseS");

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
                   "I component of depth averaged current at left face");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "u2av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "J component of depth averaged current at back face");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "wtop")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "Vertical velocity at surface");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "topz")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name",
                   "Z coordinate for surface cell");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "eta")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Surface Elevation");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "wind1")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "I component of wind stress at left face");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "wind2")) >= 0) {
    write_text_att(cdfid, vid, "units", "Nm-2");
    write_text_att(cdfid, vid, "long_name",
                   "J component of wind stress at back face");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "patm")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "u1")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "I component of current at left face");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "u2")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "J component of current at back face");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "w")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "K component of current at cell centre and Z grid");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "dens")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Density");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "dens_0")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Potential density");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
    fill = 1025.0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Kz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Kz");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "Vz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Vz");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
    fill = 0;
    nc_put_att_double(cdfid, vid, "_FillValue", fptype, 1, &fill);
  }

  if ((vid = ncw_var_id(cdfid, "ptconc")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Particle concentration");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "Cd")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Bottom drag coefficient");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "ustrcw")) >= 0) {
    write_text_att(cdfid, vid, "units", "ms-1");
    write_text_att(cdfid, vid, "long_name", "Bottom friction velocity");
    write_text_att(cdfid, vid, "coordinates", "t, sparseS");
  }

  if ((vid = ncw_var_id(cdfid, "origin")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Streamline origin");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
    if (dumpdata->togn & TOPRIGHT)
      write_text_att(cdfid, vid, "datum", "top right");
    else
      write_text_att(cdfid, vid, "datum", "bottom left");
  }

   if ((vid = ncw_var_id(cdfid, "p")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "x offset in origin cell");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "q")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "y offset in origin cell");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  if ((vid = ncw_var_id(cdfid, "r")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "z offset in origin cell");
    write_text_att(cdfid, vid, "coordinates", "t, sparse");
  }

  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "t, sparse"}
    };
    tracer_write_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "t, sparseS"}
    };
    tracer_write_nc(cdfid, dumpdata->ntrS, dumpdata->trinfo_2d, 2, attr);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", dumpdata->grid_name);
  write_text_att(cdfid, NC_GLOBAL, "description", dumpdata->grid_desc);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "codehead", codeheader);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_date_created(cdfid);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CMR/Timeseries/SHOC");
  if (dumpdata->runno >= 0)
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);
  nc_put_att_int(cdfid, NC_GLOBAL, "nce1", NC_INT, 1, &nce1);
  nc_put_att_int(cdfid, NC_GLOBAL, "nce2", NC_INT, 1, &nce2);
  nc_put_att_int(cdfid, NC_GLOBAL, "nfe1", NC_INT, 1, &nfe1);
  nc_put_att_int(cdfid, NC_GLOBAL, "nfe2", NC_INT, 1, &nfe2);
  nc_put_att_int(cdfid, NC_GLOBAL, "ns3", NC_INT, 1, &ns3);
  nc_put_att_int(cdfid, NC_GLOBAL, "ns2", NC_INT, 1, &ns2);
  if (dumpdata->tmode & FFSL)
    nc_put_att_double(cdfid, NC_GLOBAL, "dt", NC_INT, 1, &master->grid_dt);

}

void df_sp_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_sp_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t <= (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "t", df->tout);
}

void *df_sp_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int n;
  int cdfid;
  int dims[10];
  int vid;                      /* df_variable_t dimension */
  /* dimension ids */
  int recdimid;                 /* record dimension id */
  int igridid;                  /* I dimension id at grid corner */
  int jgridid;                  /* J dimension id at grid corner */
  int kgridid;                  /* K dimension id at grid corner */
  int icentreid;                /* I dimension id at grid centre */
  int jcentreid;                /* J dimension id at grid centre */
  int kcentreid;                /* K dimension id at grid centre */
  int ileftid;                  /* I dimension id at left face */
  int jleftid;                  /* J dimension id at left face */
  int ibackid;                  /* I dimension id at back face */
  int jbackid;                  /* J dimension id at back face */
  int ns3;                      /* 3D sparse size */
  int ns2;                      /* 2D sparse size */
  FILE *fp;
  int nc_mode;

  df_sp_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);

  if (dumpdata->master->mpi_rank) return(NULL);

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_sp_init_data(dumpdata, df, cdfid);
      data = (df_sp_data_t *)df->private_data;
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

  nc_def_dim(cdfid, "k_grid", df->nz + 1, &kgridid);
  nc_def_dim(cdfid, "j_grid", df->nce2 + 1, &jgridid);
  nc_def_dim(cdfid, "i_grid", df->nce1 + 1, &igridid);

  nc_def_dim(cdfid, "k_centre", df->nz, &kcentreid);
  nc_def_dim(cdfid, "j_centre", df->nce2, &jcentreid);
  nc_def_dim(cdfid, "i_centre", df->nce1, &icentreid);

  nc_def_dim(cdfid, "j_left", df->nce2, &jleftid);
  nc_def_dim(cdfid, "i_left", df->nce1 + 1, &ileftid);

  nc_def_dim(cdfid, "j_back", df->nce2 + 1, &jbackid);
  nc_def_dim(cdfid, "i_back", df->nce1, &ibackid);

  nc_def_dim(cdfid, "ns3", df->ns3, &ns3);
  nc_def_dim(cdfid, "ns2", df->ns2, &ns2);

  dims[0] = kgridid;
  nc_def_var(cdfid, "z_grid", NC_DOUBLE, 1, dims, &vid);
  dims[0] = kcentreid;
  nc_def_var(cdfid, "z_centre", NC_DOUBLE, 1, dims, &vid);
  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "x_grid", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_grid", NC_DOUBLE, 2, dims, &vid);
  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "x_centre", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_centre", NC_DOUBLE, 2, dims, &vid);
  dims[0] = ns2;
  nc_def_var(cdfid, "botz", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h1au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h1au2", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h1acell", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h2au1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h2au2", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "h2acell", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "thetau1", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "thetau2", NC_DOUBLE, 1, dims, &vid);
  nc_def_var(cdfid, "coriolis", fptype, 1, dims, &vid);

  dims[0] = ns3;
  nc_def_var(cdfid, "s2i", NC_INT, 1, dims, &vid);
  nc_def_var(cdfid, "s2j", NC_INT, 1, dims, &vid);
  nc_def_var(cdfid, "s2k", NC_INT, 1, dims, &vid);

  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid, "t", NC_DOUBLE, 1, dims, &vid);

  df_sp_init_data(dumpdata, df, cdfid);
  data = (df_sp_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    if (data->vars[n].ndims == 1) {
      df_sp_get_dimids(cdfid, data->vars[n], &dims[1]);
      ncw_def_var2(df->name,cdfid, df->vars[n], data->vars[n].type, 2, dims, &vid, df->compress);
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }
  }

  write_dump_attributes_sp(dumpdata, cdfid, fptype, df->ilower, df->nce1,
                   df->nfe1, df->jlower, df->nce2, df->nfe2, df->ns2, df->ns3, df->modulo);

  nc_enddef(cdfid);

  df_sp_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}



void df_sp_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double st, newt = t;
  df_sp_data_t *data = (df_sp_data_t *)df->private_data;
  int fid = data->fid;
  geometry_t *geom = dumpdata->master->geom;

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
    df_sp_var_t vn = data->vars[n];
    void *p = vn.v;

    if (vn.ndims == 1) {
      int nx = geom->nce1;
      int ny = geom->nce2;

      /* Get size of the sparse array, i.e. ns2 or ns3 */
      df_sp_get_dimsizes(df, vn, &count[1]);
      
      /* Get any face lengths */
      if (vn.xylocation & CL_LEFT)
	nx = geom->nfe1;
      if (vn.xylocation & CL_BACK)
	ny = geom->nfe2;

      memset(dumpdata->w1s, 0, geom->sgsiz*sizeof(double));
      memset(dumpdata->w1,  0, geom->sgsiz*sizeof(double));

      /* Cart to sparse */
      if (vn.xylocation & CL_SP2)
	c2s_2d(geom, dumpdata->w1s, (*(double ***)p), nx, ny);
      else if (vn.xylocation & CL_SP3)
	c2s_3d(geom, dumpdata->w1s, (*(double ****)p), nx, ny, geom->nz);
      else
	hd_quit("df_sparse:df_sp_write error in xylocation\n");
      
      /* Pack into dummy array and dump to file */
      pack_sparse(dumpdata->i1, count[1], dumpdata->w1s, dumpdata->w1);
      start[1] = df->ilower;
      nc_d_writesub_1d(fid, ncw_var_id(fid, df->vars[n]), start, count,
                       dumpdata->w1);
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

void df_sp_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_sp_data_t *data = (df_sp_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}


static void df_sp_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_sp_data_t *data = (df_sp_data_t *)df->private_data;
  master_t *master= dumpdata->master;
  geometry_t *geom = master->geom;
  int fid = data->fid;
  int c, cc, *s1, *s2, *s3;

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
  nc_put_vara_double(fid, ncw_var_id(fid, "z_grid"), start, count,
                     dumpdata->gridz);
  count[0] = dumpdata->nz;
  nc_put_vara_double(fid, ncw_var_id(fid, "z_centre"), start, count,
                     dumpdata->cellz);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2 + 1;
  count[1] = df->nce1 + 1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_grid"), start, count,
                   dumpdata->gridx);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_grid"), start, count,
                   dumpdata->gridy);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_centre"), start, count,
                   dumpdata->cellx);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_centre"), start, count,
                   dumpdata->celly);

  start[0] = 0;
  start[1] = 0;
  count[0] = df->ns2;
  count[1] = 0;
  pack_sparse(dumpdata->i1, count[0], geom->h1acell, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1acell"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->h2acell, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2acell"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], master->coriolis, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "coriolis"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->botz, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "botz"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->h1au1, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au1"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->h2au1, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au1"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->thetau1, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau1"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->h1au2, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h1au2"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->h2au2, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "h2au2"), start, count,
                   dumpdata->w1);
  pack_sparse(dumpdata->i1, count[0], geom->thetau2, dumpdata->w1);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "thetau2"), start, count,
                   dumpdata->w1);

  count[0] = df->ns3;
  s1 = i_alloc_1d(geom->sgsiz);
  s2 = i_alloc_1d(geom->sgsiz);
  s3 = i_alloc_1d(geom->sgsiz);
  for(cc = 1; cc <= count[0]; cc++) {
    c = geom->wsa[cc];
    s1[cc - 1] = geom->s2i[c];
    s2[cc - 1] = geom->s2j[c];
    s3[cc - 1] = geom->s2k[c];
  }
  nc_i_writesub_1d(fid, ncw_var_id(fid, "s2i"), start, count, s1);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "s2j"), start, count, s2);
  nc_i_writesub_1d(fid, ncw_var_id(fid, "s2k"), start, count, s3);
  i_free_1d(s1);
  i_free_1d(s2);
  i_free_1d(s3);
               
  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}

static int df_sp_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
                              char *name, df_sp_var_t *var)
{
  int found = 1;
  int n = 0;

  var->v = NULL;
  var->ndims = 1; // This is the only case for sparse
  var->type = nc_get_default_type(df->bpv);
  var->xylocation = CL_NONE;

  if (strcmp(name, "u1av") == 0) {
    var->v = (void **)&dumpdata->u1av;
    var->xylocation = CL_SP2|CL_LEFT;
  }

  else if (strcmp(name, "u2av") == 0) {
    var->v = (void **)&dumpdata->u2av;
    var->xylocation = CL_SP2|CL_BACK;
  }

  else if (strcmp(name, "wtop") == 0) {
    var->v = (void **)&dumpdata->wtop;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "topz") == 0) {
    var->v = (void **)&dumpdata->topz;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "eta") == 0) {
    var->v = (void **)&dumpdata->eta;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "wind1") == 0) {
    var->v = (void **)&dumpdata->wind1;
    var->xylocation = CL_SP2|CL_LEFT;
  }

  else if (strcmp(name, "wind2") == 0) {
    var->v = (void **)&dumpdata->wind2;
    var->xylocation = CL_SP2|CL_BACK;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void **)&dumpdata->patm;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "u1") == 0) {
    var->v = (void **)&dumpdata->u1;
    var->xylocation = CL_SP3|CL_LEFT;
  }

  else if (strcmp(name, "u2") == 0) {
    var->v = (void **)&dumpdata->u2;
    var->xylocation = CL_SP3|CL_BACK;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void **)&dumpdata->w;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "u1mean") == 0) {
    var->v = (void **)&dumpdata->u1m;
    var->xylocation = CL_SP3|CL_LEFT;
  }

  else if (strcmp(name, "u2mean") == 0) {
    var->v = (void **)&dumpdata->u2m;
    var->xylocation = CL_SP3|CL_BACK;
  }

  else if (strcmp(name, "u1vmean") == 0) {
    var->v = (void **)&dumpdata->u1vm;
    var->xylocation = CL_SP3|CL_LEFT;
  }

  else if (strcmp(name, "u2vmean") == 0) {
    var->v = (void **)&dumpdata->u2vm;
    var->xylocation = CL_SP3|CL_BACK;
  }

  else if (strcmp(name, "dens") == 0) {
    var->v = (void **)&dumpdata->dens;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "dens_0") == 0) {
    var->v = (void **)&dumpdata->dens_0;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "Kz") == 0) {
    var->v = (void **)&dumpdata->Kz;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "Vz") == 0) {
    var->v = (void **)&dumpdata->Vz;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "Cd") == 0) {
    var->v = (void **)&dumpdata->Cd;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "u1vh") == 0) {
    var->v = (void **)&dumpdata->u1vh;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "u2vh") == 0) {
    var->v = (void **)&dumpdata->u2vh;
    var->xylocation = CL_SP2|CL_CENTRE;
  }

  else if (strcmp(name, "origin") == 0) {
    var->v = (void **)&dumpdata->origin;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "p") == 0) {
    var->v = (void **)&dumpdata->pc;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "q") == 0) {
    var->v = (void **)&dumpdata->qc;
    var->xylocation = CL_SP3|CL_CENTRE;
  }

  else if (strcmp(name, "r") == 0) {
    var->v = (void **)&dumpdata->rc;
    var->xylocation = CL_SP3|CL_CENTRE;
  }
  else
    found = 0;

  if (!found) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) {
        var->v = (void **)&dumpdata->tr_wc[n];
        var->xylocation = CL_SP3|CL_CENTRE;
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, name) == 0) {
        var->v = (void **)&dumpdata->tr_wcS[n];
        var->xylocation = CL_SP2|CL_CENTRE;
        found = 1;
        break;
      }
    }
  }

  return found;
}


static void df_sp_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid)
{
  int i;
  df_sp_data_t *data = NULL;

  df_parse_vars(dumpdata,df,NULL,SP_ALL_VARS);
  
  // Clean up before allocating more memory
  if (df->private_data != NULL) {
    if (((df_sp_data_t *)df->private_data)->vars != NULL)
      free(((df_sp_data_t *)df->private_data)->vars);
    free(df->private_data);
  }
  
  data = (df_sp_data_t *)malloc(sizeof(df_sp_data_t));
  memset(data, 0, sizeof(df_sp_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_sp_var_t *)malloc(df->nvars * sizeof(df_sp_var_t))) == NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_sp_var_t *var = &data->vars[i];
    if (df_sp_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_sp_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

  if (df->diagn && DEBUG("dump"))
    dlog("dump",
         "Diagnostic tracers initialisation has been synchronised with \"%s\" output; period = %f s.\n",
         df->name, df->tinc);

}


static void df_sp_get_dimids(int cdfid, df_sp_var_t var, int *ns)
{
  if (var.xylocation & CL_SP2)
    *ns = ncw_dim_id(cdfid, "ns2");
  else if (var.xylocation & CL_SP3)
    *ns = ncw_dim_id(cdfid, "ns3");
  else
    hd_quit("df_sparse:df_sp_get_dimids: xylocation has to be either ns2 or ns3 - got %d\n", var.xylocation);
}

static void df_sp_get_dimsizes(dump_file_t *df, df_sp_var_t var, size_t *ns)
{
  if (var.xylocation & CL_SP2)
    *ns = df->ns2;
  else if (var.xylocation & CL_SP3) {
    *ns = df->ns3;
    if (df->klower == df->nz - 1)
      *ns = df->ns2;
  } else
    hd_quit("df_sparse:df_sp_get_dimsizes error in xylocation\n");
}

/*-------------------------------------------------------------------*/
/* Dumps window geometry to file                                     */
/*-------------------------------------------------------------------*/
void dump_windows(master_t *master, geometry_t **window, char *name, char *iname)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int cdfid, n, c;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  int dims[2];
  int vid;
  long t;
  /* dimension ids */
  int sgsizid;            /* 3D sparse dimension           */
  int sgsizSid;           /* 2D sparse dimension           */
  int kgridid;            /* K dimension id at grid corner */
  int oid;
  int nwinsid;
  int nwins = geom->nwindows;
  char key[MAXSTRLEN];
  int *d1;
  int ct;
  int ncmode = overwrite_output() ? NC_CLOBBER:NC_NOCLOBBER;

  if (dumpdata->master->mpi_rank) return;

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
  ncw_def_dim(name, cdfid, "sgsiz",    geom->sgsiz,  &sgsizid);
  ncw_def_dim(name, cdfid, "sgsizS",   geom->sgsizS, &sgsizSid);
  ncw_def_dim(name, cdfid, "k_grid",   geom->nz + 1, &kgridid);
  ncw_def_dim(name, cdfid, "nwp1",     nwins+1,      &nwinsid);
  /* Variables of lenght 1 */
  ncw_def_dim(name, cdfid, "one", 1, &oid);
  
  /* Global variables */
  dims[0] = kgridid;
  ncw_def_var(name, cdfid, "z_grid", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate at grid layer faces");  
  write_text_att(cdfid, vid, "coordinate_type", "Z");
  dims[0] = sgsizid;
  ncw_def_var(name, cdfid, "wn", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Global-local window map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  ncw_def_var2(name, cdfid, "sc", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Global-local wet coordinate map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  ncw_def_var2(name, cdfid, "ac", NC_INT, 1, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Global-local auxiliary coordinate map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  dims[0] = oid;
  ncw_def_var(name, cdfid, "nz", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Number of layers");
  ncw_def_var(name, cdfid, "sednz", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Number of Sediment layers");
  ncw_def_var(name, cdfid, "zoomf", NC_INT, 1, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Zoom factor");

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
  
  ncw_def_var(name, cdfid, "sb_t_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "sb_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "sb_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "sgsiz_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "sgsizS_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "nbpt_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte1S_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nbpte2S_w", NC_INT, 1, dims, &vid);
  
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
  ncw_def_var(name, cdfid, "x2_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a3_e1_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "x3_e1_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "v3_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "b3_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "a3_e2_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "x3_e2_w", NC_INT, 1, dims, &vid);
  
  ncw_def_var(name, cdfid, "ns2mS_w", NC_INT, 1, dims, &vid);
  ncw_def_var(name, cdfid, "nm2sS_w", NC_INT, 1, dims, &vid);
  
  /*
   * 2D - One 1D array per window
   *      Figure out a max size/window for each variable
   *
   * Define the master dimensions as we go along
   */
  dims[0] = nwinsid; // always
  ct = window[1]->a2_t;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->a2_t > ct ? window[n]->a2_t : ct);
  ncw_def_dim(name, cdfid, "a2_t_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bot_t",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of cell centre bottom");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  ncw_def_var2(name, cdfid, "sur_t",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of cell centre surface");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");
  
  ncw_def_var2(name, cdfid, "nsur_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of cell centre updated surface");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  ct = window[1]->x2_e1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->x2_e1 > ct ? window[n]->x2_e1 : ct);
  ncw_def_dim(name, cdfid, "x2_e1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bot_e1",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of e1 bottom");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");
  ncw_def_var2(name, cdfid, "sur_e1",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of e1 surface");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  ct = window[1]->x2_e2;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->x2_e2 > ct ? window[n]->x2_e2 : ct);
  ncw_def_dim(name, cdfid, "x2_e2_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bot_e2",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of e2 bottom");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");
  ncw_def_var2(name, cdfid, "sur_e2",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse coordinate of e2 surface");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparseS");

  /*
   * sgsiz can be very big for large domains, so figure out a max
   */
  ct = window[1]->sgsiz;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->sgsiz > ct ? window[n]->sgsiz : ct);
  ncw_def_dim(name, cdfid, "sgsiz_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "wsa",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Local - global sparse map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "m2d", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "3D to 2D sparse map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "xp1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in +e1 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "xm1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in -e1 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "yp1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in +e2 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "ym1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in -e2 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "zp1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in +z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "zm1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in -z direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "xmyp1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in -e1+e2 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "xpym1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse map in +e1-e2 direction");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "wgst", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "lateral boundary ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->nbpt;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpt > ct ? window[n]->nbpt : ct);
  ncw_def_dim(name, cdfid, "nbpt_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bpt",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "bin",  NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "bin2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Two interior map to ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->nbpte1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte1 > ct ? window[n]->nbpte1 : ct);
  ncw_def_dim(name, cdfid, "nbpte1_max",  ct,  &dims[1]);
  ncw_def_var2(name, cdfid, "bpte1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "e1 3D ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "bine1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to e1 3D ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->nbpte1S;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte1S > ct ? window[n]->nbpte1S : ct);
  ncw_def_dim(name, cdfid, "nbpte1S_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bpte1S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "e1 2D ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "bine1S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to e1 2D ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->nbpte2;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte2 > ct ? window[n]->nbpte2 : ct);
  ncw_def_dim(name, cdfid, "nbpte2_max",  ct,  &dims[1]);
  ncw_def_var2(name, cdfid, "bpte2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "e2 3D ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ncw_def_var2(name, cdfid, "bine2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to e2 3D ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->nbpte2S;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nbpte2S > ct ? window[n]->nbpte2S : ct);
  ncw_def_dim(name, cdfid, "nbpte2S_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "bpte2S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "e2 2D ghost cell map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "bine2S", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Interior map to e2 2D ghost cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  
  ct = window[1]->n2_t;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_t > ct ? window[n]->n2_t : ct);
  ncw_def_dim(name, cdfid, "n2_t_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 2D work array for cell centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  
  ct = window[1]->n2_e1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_e1 > ct ? window[n]->n2_e1 : ct);
  ncw_def_dim(name, cdfid, "n2_e1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_e1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 2D work array for e1 cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->n2_e2;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n2_e2 > ct ? window[n]->n2_e2 : ct);
  ncw_def_dim(name, cdfid, "n2_e2_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w2_e2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 2D work array for e2 cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  
  ct = window[1]->n3_t;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_t > ct ? window[n]->n3_t : ct);
  ncw_def_dim(name, cdfid, "n3_t_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_t", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 3D work array for cell centres");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->n3_e1;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_e1 > ct ? window[n]->n3_e1 : ct);
  ncw_def_dim(name, cdfid, "n3_e1_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_e1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 3D work array for e1 cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  ct = window[1]->n3_e2;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->n3_e2 > ct ? window[n]->n3_e2 : ct);
  ncw_def_dim(name, cdfid, "n3_e2_max", geom->n3_e2, &dims[1]);
  ncw_def_var2(name, cdfid, "w3_e2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse 3D work array for e2 cells");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  dims[1] = sgsizid;
  ncw_def_var2(name, cdfid, "s2i", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse to Cartesian i map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "s2j", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse to Cartesian j map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "s2k", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Sparse to Cartesian k map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  
  ct = window[1]->nm2s;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->nm2s > ct ? window[n]->nm2s : ct);
  ncw_def_dim(name, cdfid, "nm2s_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "m2s",   NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Master to slave map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "m2se1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Master to slave map for e1 faces");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "m2se2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Master to slave map for e2 faces");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  
  ct = window[1]->ns2m;
  for (n=2; n<=nwins; n++)
    ct = (window[n]->ns2m > ct ? window[n]->ns2m : ct);
  ncw_def_dim(name, cdfid, "ns2m_max", ct, &dims[1]);
  ncw_def_var2(name, cdfid, "s2m",   NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Slave to master map");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "s2me1", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Slave to master map for e1 faces");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");
  ncw_def_var2(name, cdfid, "s2me2", NC_INT, 2, dims, &vid, 1);
  write_text_att(cdfid, vid, "long_name", "Slave to master map for e2 faces");
  write_text_att(cdfid, vid, "units", "");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  dims[1] = sgsizSid;
  ncw_def_var(name, cdfid, "botz", NC_DOUBLE, 2, dims, &vid);
  write_text_att(cdfid, vid, "long_name", "Bottom depth");
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "cartesian_axis", "sparse");

  // } nwindows

  /* Global attributes */
  write_text_att(cdfid, NC_GLOBAL, "Parameter_filename", iname);
  nc_put_att_int(cdfid, NC_GLOBAL, "Windows", NC_INT, 1, &window[1]->nwindows);
  write_text_att(cdfid, NC_GLOBAL, "Version", version);
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
  count[0] = geom->sgsiz;
  for (c = 1; c <= geom->sgnum; c++) d1[c] = geom->fm[c].wn;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wn"), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) d1[c] = geom->fm[c].sc;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sc"), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) d1[c] = geom->fm[c].ac;
  ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ac"), start, count, d1);
  i_free_1d(d1);

  nc_put_var_int(cdfid, ncw_var_id(cdfid, "nz"), &geom->nz);
  nc_put_var_int(cdfid, ncw_var_id(cdfid, "sednz"), &geom->sednz);
  nc_put_var_int(cdfid, ncw_var_id(cdfid, "zoomf"), &window[1]->zoomf);

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
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sgsiz_w"), start, count, &window[n]->sgsiz);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sgsizS_w"), start, count, &window[n]->sgsizS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpt_w"), start, count, &window[n]->nbpt);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte1_w"), start, count, &window[n]->nbpte1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte1S_w"), start, count, &window[n]->nbpte1S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte2_w"), start, count, &window[n]->nbpte2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nbpte2S_w"), start, count, &window[n]->nbpte2S);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_t_w"), start, count, &window[n]->n2_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_t_w"), start, count, &window[n]->n3_t);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_e1_w"), start, count, &window[n]->n2_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_e1_w"), start, count, &window[n]->n3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n2_e2_w"), start, count, &window[n]->n2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "n3_e2_w"), start, count, &window[n]->n3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2s_w"), start, count, &window[n]->nm2s);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2m_w"), start, count, &window[n]->ns2m);
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
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "x2_e2_w"), start, count, &window[n]->x2_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v3_e1_w"), start, count, &window[n]->v3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b3_e1_w"), start, count, &window[n]->b3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a3_e1_w"), start, count, &window[n]->a3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "x3_e1_w"), start, count, &window[n]->x3_e1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "v3_e2_w"), start, count, &window[n]->v3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "b3_e2_w"), start, count, &window[n]->b3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "a3_e2_w"), start, count, &window[n]->a3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "x3_e2_w"), start, count, &window[n]->x3_e2);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ns2mS_w"), start, count, &window[n]->ns2mS);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nm2sS_w"), start, count, &window[n]->nm2sS);

    /* 2D arrays */
    count[0] = 1;
    start[1] = 0;
    count[1] = window[n]->a2_t;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_t"), start, count, &window[n]->bot_t[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_t"), start, count, &window[n]->sur_t[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "nsur_t"), start, count, &window[n]->nsur_t[1]);

    count[1] = window[n]->x2_e1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_e1"), start, count, &window[n]->bot_e1[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_e1"), start, count, &window[n]->sur_e1[1]);

    count[1] = window[n]->x2_e2;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bot_e2"), start, count, &window[n]->bot_e2[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "sur_e2"), start, count, &window[n]->sur_e2[1]);

    count[1] = window[n]->sgsiz;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2d"), start, count, window[n]->m2d);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wsa"), start, count, window[n]->wsa);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "xp1"), start, count, window[n]->xp1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "xm1"), start, count, window[n]->xm1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "yp1"), start, count, window[n]->yp1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "ym1"), start, count, window[n]->ym1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zp1"), start, count, window[n]->zp1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "zm1"), start, count, window[n]->zm1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "xmyp1"), start, count, window[n]->xmyp1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "xpym1"), start, count, window[n]->xpym1);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "wgst"), start, count, window[n]->wgst);

    count[1] = window[n]->nbpt;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpt"), start, count, &window[n]->bpt[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bin"), start, count, &window[n]->bin[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bin2"), start, count, &window[n]->bin2[1]);

    count[1] = window[n]->nbpte1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1"), start, count, &window[n]->bpte1[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1"), start, count, &window[n]->bine1[1]);

    count[1] = window[n]->nbpte1S;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte1S"), start, count, &window[n]->bpte1S[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine1S"), start, count, &window[n]->bine1S[1]);

    count[1] = window[n]->nbpte2;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte2"), start, count, &window[n]->bpte2[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine2"), start, count, &window[n]->bine2[1]);

    count[1] = window[n]->nbpte2S;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bpte2S"), start, count, &window[n]->bpte2S[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "bine2S"), start, count, &window[n]->bine2S[1]);

    count[1] = window[n]->n2_t;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_t"), start, count, &window[n]->w2_t[1]);
    count[1] = window[n]->n3_t;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_t"), start, count, &window[n]->w3_t[1]);
    count[1] = window[n]->n2_e1;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e1"), start, count, &window[n]->w2_e1[1]);
    count[1] = window[n]->n3_e1;
    sprintf(key, "w3_e1_%d", n);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e1"), start, count, &window[n]->w3_e1[1]);

    count[1] = window[n]->n2_e2;
    sprintf(key, "w2_e2_%d", n);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w2_e2"), start, count, &window[n]->w2_e2[1]);
    count[1] = window[n]->n3_e2;
    sprintf(key, "w3_e2_%d", n);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "w3_e2"), start, count, &window[n]->w3_e2[1]);

    count[1] = window[n]->sgsiz;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2i"), start, count, window[n]->s2i);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2j"), start, count, window[n]->s2j);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2k"), start, count, window[n]->s2k);

    count[1] = window[n]->nm2s;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2s"), start, count, &window[n]->m2s[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2se1"), start, count, &window[n]->m2se1[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "m2se2"), start, count, &window[n]->m2se2[1]);

    count[1] = window[n]->ns2m;
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2m"), start, count, &window[n]->s2m[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2me1"), start, count, &window[n]->s2me1[1]);
    ncw_put_vara_int(name, cdfid, ncw_var_id(cdfid, "s2me2"), start, count, &window[n]->s2me2[1]);

    count[1] = window[n]->sgsizS;
    sprintf(key, "botz_%d", n);
    ncw_put_vara_double(name, cdfid, ncw_var_id(cdfid, "botz"), start, count, window[n]->botz);
  }

  ncw_close(name, cdfid);
  /* check_window_map(window, name); */
}

/* END dump_windows()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads window geometry from file                                   */
/* This is the original flat file format                             */
/*-------------------------------------------------------------------*/
static void read_windows_flat(geometry_t *geom, geometry_t **window, char *name)
{
  int fid, n;
  int c, cc, cg, *d1;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  size_t d;
  char key[MAXSTRLEN];



  double *layers;
  int ncerr;

  /* Open the dump file for reading */
  if ((ncerr = nc_open(name, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("Can't find window map file %s\n", name);
    hd_quit((char *)nc_strerror(ncerr));
  }

  geom->zoomf = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  d1 = i_alloc_1d(geom->sgsiz);
  count[0] = geom->sgsiz;
  sprintf(key, "wn");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].wn = d1[c];
  sprintf(key, "sc");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].sc = d1[c];
  sprintf(key, "ac");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].ac = d1[c];
  i_free_1d(d1);

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

  /* Get dimensions */
  for (n = 1; n <= geom->nwindows; n++) {
    double *botz;

    emstag(LINFO,"read_windows_flat","Reading window %d map\n", n);

    window[n] = window_alloc();
    window[n]->wn = n;

    sprintf(key, "nwindows_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nwindows = (int)d;
    if (window[n]->nwindows != geom->nwindows)
      hd_quit("Incompatible number of windows in file %s\n", name);

    sprintf(key, "enon_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->enon = (int)d;
    sprintf(key, "ewet_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->ewet = (int)d;
    sprintf(key, "snon_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->snon = (int)d;
    sprintf(key, "enonS_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->enonS = (int)d;
    sprintf(key, "ewetS_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->ewetS = (int)d;
    sprintf(key, "snonS_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->snonS = (int)d;
    sprintf(key, "nz_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nz = (int)d;
    sprintf(key, "sednz_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->sednz = (int)d;
    sprintf(key, "zoomf_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->zoomf = (int)d;

    sprintf(key, "sgsiz_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->sgsiz = (int)d;
    window[n]->sgnum = window[n]->sgsiz - 1;
    sprintf(key, "sgsizS_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->sgsizS = (int)d;
    window[n]->sgnumS = window[n]->sgsizS - 1;

    sprintf(key, "n2_t_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n2_t = (int)d - 1;
    sprintf(key, "n3_t_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n3_t = (int)d - 1;

    sprintf(key, "n2_e1_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n2_e1 = (int)d - 1;
    sprintf(key, "n3_e1_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n3_e1 = (int)d - 1;

    sprintf(key, "n2_e2_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n2_e2 = (int)d - 1;
    sprintf(key, "n3_e2_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->n3_e2 = (int)d - 1;

    sprintf(key, "nbpt_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nbpt = (int)d - 1;

    sprintf(key, "nbpte1_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nbpte1 = (int)d - 1;
    sprintf(key, "nbpte1S_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nbpte1S = (int)d - 1;


    sprintf(key, "nbpte2_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nbpte2 = (int)d - 1;
    sprintf(key, "nbpte2S_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nbpte2S = (int)d - 1;

    sprintf(key, "nm2s_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->nm2s = (int)d - 1;
    sprintf(key, "ns2m_%d", n);
    nc_inq_dimlen(fid, ncw_dim_id(fid, key), &d);
    window[n]->ns2m = (int)d - 1;

    count[0] = 1;
    sprintf(key, "nbptS_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbptS);
    sprintf(key, "nbe1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe1);
    sprintf(key, "nbe1S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe1S);
    sprintf(key, "nbe2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe2);
    sprintf(key, "nbe2S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe2S);
    sprintf(key, "v2_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_t);
    sprintf(key, "b2_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_t);
    sprintf(key, "a2_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_t);
    sprintf(key, "v3_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_t);
    sprintf(key, "b3_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_t);
    sprintf(key, "a3_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_t);

    sprintf(key, "v2_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_e1);
    sprintf(key, "b2_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_e1);
    sprintf(key, "a2_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_e1);
    sprintf(key, "x2_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x2_e1);
    sprintf(key, "v3_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_e1);
    sprintf(key, "b3_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_e1);
    sprintf(key, "a3_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_e1);
    sprintf(key, "x3_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x3_e1);

    sprintf(key, "v2_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_e2);
    sprintf(key, "b2_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_e2);
    sprintf(key, "a2_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_e2);
    sprintf(key, "x2_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x2_e2);
    sprintf(key, "v3_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_e2);
    sprintf(key, "b3_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_e2);
    sprintf(key, "a3_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_e2);
    sprintf(key, "x3_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x3_e2);

    sprintf(key, "nm2sS_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nm2sS);
    sprintf(key, "ns2mS_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->ns2mS);

    alloc_geom(window[n], (MAP_A | GRID_A | WINDOW_A));
    window[n]->w2_t = i_alloc_1d(window[n]->n2_t + 1);
    window[n]->w3_t = i_alloc_1d(window[n]->n3_t + 1);
    window[n]->w2_e1 = i_alloc_1d(window[n]->n2_e1 + 1);
    window[n]->w3_e1 = i_alloc_1d(window[n]->n3_e1 + 1);
    window[n]->w2_e2 = i_alloc_1d(window[n]->n2_e2 + 1);
    window[n]->w3_e2 = i_alloc_1d(window[n]->n3_e2 + 1);
    window[n]->sur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->nsur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->bot_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->sur_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
    window[n]->bot_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
    window[n]->sur_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
    window[n]->bot_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
    window[n]->m2d = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2i = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2j = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2k = i_alloc_1d(window[n]->enon + 1);
    window[n]->m2s = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->m2se1 = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->m2se2 = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->s2m = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->s2me1 = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->s2me2 = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->bpt = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin2 = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bpte1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bine1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bpte1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->bine1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->bpte2 = i_alloc_1d(window[n]->nbpte2 + 1);
    window[n]->bine2 = i_alloc_1d(window[n]->nbpte2 + 1);
    window[n]->bpte2S = i_alloc_1d(window[n]->nbpte2S + 1);
    window[n]->bine2S = i_alloc_1d(window[n]->nbpte2S + 1);
    window[n]->wgst = i_alloc_1d(window[n]->enon + 1);

    count[0] = window[n]->a2_t+1;
    sprintf(key, "bot_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bot_t);
    sprintf(key, "sur_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->sur_t);
    sprintf(key, "nsur_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->nsur_t);

    count[0] = window[n]->x2_e1+1;
    sprintf(key, "bot_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bot_e1);
    sprintf(key, "sur_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->sur_e1);

    count[0] = window[n]->x2_e2+1;
    sprintf(key, "bot_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bot_e2);
    sprintf(key, "sur_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->sur_e2);

    count[0] = window[n]->sgsiz;
    sprintf(key, "m2d_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->m2d);
    sprintf(key, "wsa_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->wsa);
    sprintf(key, "xp1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xp1);
    sprintf(key, "xm1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xm1);
    sprintf(key, "yp1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->yp1);
    sprintf(key, "ym1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->ym1);
    sprintf(key, "zp1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->zp1);
    sprintf(key, "zm1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->zm1);
    sprintf(key, "xmyp1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xmyp1);
    sprintf(key, "xpym1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xpym1);
    sprintf(key, "wgst_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->wgst);

    count[0] = window[n]->nbpt+1;
    sprintf(key, "bpt_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bpt);
    sprintf(key, "bin_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bin);
    sprintf(key, "bin2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bin2);

    count[0] = window[n]->nbpte1+1;
    sprintf(key, "bpte1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bpte1);
    sprintf(key, "bine1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bine1);
    count[0] = window[n]->nbpte1S+1;
    sprintf(key, "bpte1S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bpte1S);
    sprintf(key, "bine1S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bine1S);

    count[0] = window[n]->nbpte2+1;
    sprintf(key, "bpte2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bpte2);
    sprintf(key, "bine2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bine2);
    count[0] = window[n]->nbpte2S+1;
    sprintf(key, "bpte2S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bpte2S);
    sprintf(key, "bine2S_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->bine2S);

    count[0] = window[n]->n2_t+1;
    sprintf(key, "w2_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w2_t);
    count[0] = window[n]->n3_t+1;
    sprintf(key, "w3_t_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w3_t);

    count[0] = window[n]->n2_e1+1;
    sprintf(key, "w2_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w2_e1);
    count[0] = window[n]->n3_e1+1;
    sprintf(key, "w3_e1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w3_e1);

    count[0] = window[n]->n2_e2+1;
    sprintf(key, "w2_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w2_e2);
    count[0] = window[n]->n3_e2+1;
    sprintf(key, "w3_e2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->w3_e2);

    count[0] = window[n]->sgsiz;
    sprintf(key, "s2i_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2i);
    count[0] = window[n]->sgsiz;
    sprintf(key, "s2j_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2j);
    count[0] = window[n]->sgsiz;
    sprintf(key, "s2k_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2k);

    count[0] = window[n]->nm2s+1;
    sprintf(key, "m2s_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->m2s);
    sprintf(key, "m2se1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->m2se1);
    sprintf(key, "m2se2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->m2se2);
    count[0] = window[n]->ns2m+1;
    sprintf(key, "s2m_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2m);
    sprintf(key, "s2me1_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2me1);
    sprintf(key, "s2me2_%d", n);
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2me2);

    botz = d_alloc_1d(window[n]->sgsizS);
    count[0] = window[n]->sgsizS;
    sprintf(key, "botz_%d", n);
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
    emstag(LINFO,"read_windows_flat","Window %d map OK\n", n);
  }
  ncw_close(name, fid);
}

/* END read_windows_flat()                                           */
/*-------------------------------------------------------------------*/

/* 
 * Reads window geometry from file
 *
 * This is the modified Windows based format. This avoids excessively
 * large numbers of dimension for number of windows > 38 and hence
 * does not run into the NC_MAX_DIMS limit of 1024.
 */
static void read_windows_wb(geometry_t *geom, geometry_t **window, char *name)
{
  int fid, n;
  int c, cc, cg, *d1;
  size_t start[4] = {0, 0, 0, 0};
  size_t count[4] = {0, 0, 0, 0};
  size_t d;
  char key[MAXSTRLEN];

  int oid;

  double *layers;
  int nz, sednz, zoomf;
  int nwin;

  /* Open the dump file for reading */
  ncw_open(name, NC_NOWRITE, &fid);

  geom->zoomf = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  d1 = i_alloc_1d(geom->sgsiz);
  count[0] = geom->sgsiz;
  sprintf(key, "wn");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].wn = d1[c];
  sprintf(key, "sc");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].sc = d1[c];
  sprintf(key, "ac");
  nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, d1);
  for (c = 1; c <= geom->sgnum; c++) geom->fm[c].ac = d1[c];
  i_free_1d(d1);

  /*
   * Get layers info and check
   */
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

  /*
   * Get global variables, i.e. ones that are the same across all windows
   */
  ncw_inq_dimlen(name, fid, ncw_dim_id(fid, "nwp1"), &d);
  nwin = d-1; // nwp1 is numnber of windows plus 1
  if (nwin != geom->nwindows)
    hd_quit("Incompatible number of windows in file %s %d vs %d\n", name, nwin, geom->nwindows);

  oid = ncw_dim_id(fid, "one");
  nc_get_var_int(fid, ncw_var_id(fid, "nz"), &nz);
  nc_get_var_int(fid, ncw_var_id(fid, "sednz"), &sednz);
  nc_get_var_int(fid, ncw_var_id(fid, "zoomf"), &zoomf);

  /* Get dimensions */
  for (n = 1; n <= geom->nwindows; n++) {
    double *botz;

    emstag(LINFO,"read_windows_wb","Reading window %d map\n", n);

    window[n] = window_alloc();
    window[n]->wn = n;

    window[n]->nwindows = nwin;
    window[n]->nz    = nz;
    window[n]->sednz = sednz;
    window[n]->zoomf = zoomf;

    start[0] = n;
    count[0] = 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "enon_w"), start, count, &window[n]->enon);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ewet_w"), start, count, &window[n]->ewet);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "snon_w"), start, count, &window[n]->snon);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "enonS_w"), start, count, &window[n]->enonS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ewetS_w"), start, count, &window[n]->ewetS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "snonS_w"), start, count, &window[n]->snonS);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "sgsiz_w"), start, count, &window[n]->sgsiz);
    window[n]->sgnum = window[n]->sgsiz - 1;
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "sgsizS_w"), start, count, &window[n]->sgsizS);
    window[n]->sgnumS = window[n]->sgsizS - 1;

    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_t_w"), start, count, &window[n]->n2_t);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_t_w"), start, count, &window[n]->n3_t);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_e1_w"), start, count, &window[n]->n2_e1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_e1_w"), start, count, &window[n]->n3_e1);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n2_e2_w"), start, count, &window[n]->n2_e2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "n3_e2_w"), start, count, &window[n]->n3_e2);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpt_w"), start, count, &window[n]->nbpt);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte1_w"), start, count, &window[n]->nbpte1);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte1S_w"), start, count, &window[n]->nbpte1S);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte2_w"), start, count, &window[n]->nbpte2);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nbpte2S_w"), start, count, &window[n]->nbpte2S);
    
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "nm2s_w"), start, count, &window[n]->nm2s);
    ncw_get_vara_int(name, fid, ncw_var_id(fid, "ns2m_w"), start, count, &window[n]->ns2m);
    
    count[0] = 1;
    sprintf(key, "nbptS_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbptS);
    sprintf(key, "nbe1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe1);
    sprintf(key, "nbe1S_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe1S);
    sprintf(key, "nbe2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe2);
    sprintf(key, "nbe2S_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nbe2S);
    sprintf(key, "v2_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_t);
    sprintf(key, "b2_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_t);
    sprintf(key, "a2_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_t);
    sprintf(key, "v3_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_t);
    sprintf(key, "b3_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_t);
    sprintf(key, "a3_t_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_t);

    sprintf(key, "v2_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_e1);
    sprintf(key, "b2_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_e1);
    sprintf(key, "a2_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_e1);
    sprintf(key, "x2_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x2_e1);
    sprintf(key, "v3_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_e1);
    sprintf(key, "b3_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_e1);
    sprintf(key, "a3_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_e1);
    sprintf(key, "x3_e1_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x3_e1);

    sprintf(key, "v2_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v2_e2);
    sprintf(key, "b2_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b2_e2);
    sprintf(key, "a2_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a2_e2);
    sprintf(key, "x2_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x2_e2);
    sprintf(key, "v3_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->v3_e2);
    sprintf(key, "b3_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->b3_e2);
    sprintf(key, "a3_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->a3_e2);
    sprintf(key, "x3_e2_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->x3_e2);

    sprintf(key, "nm2sS_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nm2sS);
    sprintf(key, "ns2mS_w");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->ns2mS);

    alloc_geom(window[n], (MAP_A | GRID_A | WINDOW_A));
    window[n]->w2_t = i_alloc_1d(window[n]->n2_t + 1);
    window[n]->w3_t = i_alloc_1d(window[n]->n3_t + 1);
    window[n]->w2_e1 = i_alloc_1d(window[n]->n2_e1 + 1);
    window[n]->w3_e1 = i_alloc_1d(window[n]->n3_e1 + 1);
    window[n]->w2_e2 = i_alloc_1d(window[n]->n2_e2 + 1);
    window[n]->w3_e2 = i_alloc_1d(window[n]->n3_e2 + 1);
    window[n]->sur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->nsur_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->bot_t = i_alloc_1d(window[n]->a2_t + 1);
    window[n]->sur_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
    window[n]->bot_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
    window[n]->sur_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
    window[n]->bot_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
    window[n]->m2d = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2i = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2j = i_alloc_1d(window[n]->enon + 1);
    window[n]->s2k = i_alloc_1d(window[n]->enon + 1);
    window[n]->m2s = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->m2se1 = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->m2se2 = i_alloc_1d(window[n]->nm2s + 1);
    window[n]->s2m = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->s2me1 = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->s2me2 = i_alloc_1d(window[n]->ns2m + 1);
    window[n]->bpt = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bin2 = i_alloc_1d(window[n]->nbpt + 1);
    window[n]->bpte1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bine1 = i_alloc_1d(window[n]->nbpte1 + 1);
    window[n]->bpte1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->bine1S = i_alloc_1d(window[n]->nbpte1S + 1);
    window[n]->bpte2 = i_alloc_1d(window[n]->nbpte2 + 1);
    window[n]->bine2 = i_alloc_1d(window[n]->nbpte2 + 1);
    window[n]->bpte2S = i_alloc_1d(window[n]->nbpte2S + 1);
    window[n]->bine2S = i_alloc_1d(window[n]->nbpte2S + 1);
    window[n]->wgst = i_alloc_1d(window[n]->enon + 1);

    start[1] = 0;
    count[1] = window[n]->a2_t;
    sprintf(key, "bot_t");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bot_t[1]);
    sprintf(key, "sur_t");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->sur_t[1]);
    sprintf(key, "nsur_t");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->nsur_t[1]);

    count[1] = window[n]->x2_e1;
    sprintf(key, "bot_e1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bot_e1[1]);
    sprintf(key, "sur_e1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->sur_e1[1]);

    count[1] = window[n]->x2_e2;
    sprintf(key, "bot_e2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bot_e2[1]);
    sprintf(key, "sur_e2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->sur_e2[1]);

    count[1] = window[n]->sgsiz;
    sprintf(key, "m2d");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->m2d);
    sprintf(key, "wsa");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->wsa);
    sprintf(key, "xp1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xp1);
    sprintf(key, "xm1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xm1);
    sprintf(key, "yp1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->yp1);
    sprintf(key, "ym1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->ym1);
    sprintf(key, "zp1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->zp1);
    sprintf(key, "zm1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->zm1);
    sprintf(key, "xmyp1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xmyp1);
    sprintf(key, "xpym1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->xpym1);
    sprintf(key, "wgst");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->wgst);

    count[1] = window[n]->nbpt;
    sprintf(key, "bpt");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bpt[1]);
    sprintf(key, "bin");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bin[1]);
    sprintf(key, "bin2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bin2[1]);
    
    count[1] = window[n]->nbpte1;
    sprintf(key, "bpte1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bpte1[1]);
    sprintf(key, "bine1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bine1[1]);

    count[1] = window[n]->nbpte1S;
    sprintf(key, "bpte1S");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bpte1S[1]);
    sprintf(key, "bine1S");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bine1S[1]);

    count[1] = window[n]->nbpte2;
    sprintf(key, "bpte2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bpte2[1]);
    sprintf(key, "bine2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bine2[1]);

    count[1] = window[n]->nbpte2S;
    sprintf(key, "bpte2S");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bpte2S[1]);
    sprintf(key, "bine2S");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->bine2S[1]);

    count[1] = window[n]->n2_t;
    sprintf(key, "w2_t");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w2_t[1]);

    count[1] = window[n]->n3_t;
    sprintf(key, "w3_t");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w3_t[1]);

    count[1] = window[n]->n2_e1;
    sprintf(key, "w2_e1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w2_e1[1]);

    count[1] = window[n]->n3_e1;
    sprintf(key, "w3_e1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w3_e1[1]);

    count[1] = window[n]->n2_e2;
    sprintf(key, "w2_e2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w2_e2[1]);
    count[1] = window[n]->n3_e2;
    sprintf(key, "w3_e2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->w3_e2[1]);

    count[1] = window[n]->sgsiz;
    sprintf(key, "s2i");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2i);
    count[1] = window[n]->sgsiz;
    sprintf(key, "s2j");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2j);
    count[1] = window[n]->sgsiz;
    sprintf(key, "s2k");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, window[n]->s2k);

    count[1] = window[n]->nm2s;
    sprintf(key, "m2s");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->m2s[1]);
    sprintf(key, "m2se1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->m2se1[1]);
    sprintf(key, "m2se2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->m2se2[1]);
    count[1] = window[n]->ns2m;
    sprintf(key, "s2m");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->s2m[1]);
    sprintf(key, "s2me1");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->s2me1[1]);
    sprintf(key, "s2me2");
    nc_get_vara_int(fid, ncw_var_id(fid, key), start, count, &window[n]->s2me2[1]);

    botz = d_alloc_1d(window[n]->sgsizS);
    count[1] = window[n]->sgsizS;
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
    emstag(LINFO,"read_windows_wb","Window %d map OK\n", n);
  }
  ncw_close(name, fid);
}

/* END read_windows_wb()                                             */
/*-------------------------------------------------------------------*/

/*
 * Gateway function to figure out which format the window map file is
 * in
 */
void read_windows(geometry_t *geom, geometry_t **window, char *name)
{
  int ncid, wb=0;
  
  ncw_open(name, NC_NOWRITE, &ncid);
  wb = ncw_var_exists(ncid, "enon_w");
  ncw_close(name, ncid);

  if (wb)
    read_windows_wb(geom, window, name);
  else
    read_windows_flat(geom, window, name);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void check_window_map(geometry_t **window, char *name)
{
  geometry_t **win;
  int n, cc;

  win = (geometry_t **)p_alloc_1d(window[1]->nwindows);
  read_windows(geom, win, name);

  for (n = 1; n <= geom->nwindows; n++) {
    printf("Checking window %d\n", n);
    if(window[n]->enon != win[n]->enon)printf("enon_%d %d\n", n, win[n]->enon);
    if(window[n]->ewet != win[n]->ewet)printf("ewet_%d %d\n", n, win[n]->ewet);
    if(window[n]->snon != win[n]->snon)printf("snon_%d %d\n", n, win[n]->snon);
    if(window[n]->enonS != win[n]->enonS)printf("enonS_%d %d\n", n, win[n]->enonS);
    if(window[n]->ewetS != win[n]->ewetS)printf("ewetS_%d %d\n", n, win[n]->ewetS);
    if(window[n]->snonS != win[n]->snonS)printf("snonS_%d %d\n", n, win[n]->snonS);
    if(window[n]->nz != win[n]->nz)printf("nz_%d %d\n", n, win[n]->nz);
    if(window[n]->sednz != win[n]->sednz)printf("sednz_%d %d\n", n, win[n]->sednz);
    if(window[n]->zoomf != win[n]->zoomf)printf("zoomf_%d %d\n", n, win[n]->zoomf);
    if(window[n]->sgsiz != win[n]->sgsiz)printf("sgsiz_%d %d\n", n, win[n]->sgsiz);
    if(window[n]->sgnum != win[n]->sgnum)printf("sgnum_%d %d\n", n, win[n]->sgnum);
    if(window[n]->sgsizS != win[n]->sgsizS)printf("sgsizS_%d %d\n", n, win[n]->sgsizS);
    if(window[n]->sgnumS != win[n]->sgnumS)printf("sgnumS_%d %d\n", n, win[n]->sgnumS);
    if(window[n]->n2_t != win[n]->n2_t)printf("n2_t_%d %d\n", n, win[n]->n2_t);
    if(window[n]->n3_t != win[n]->n3_t)printf("n3_t_%d %d\n", n, win[n]->n3_t);
    if(window[n]->n2_e1 != win[n]->n2_e1)printf("n2_e1_%d %d\n", n, win[n]->n2_e1);
    if(window[n]->n3_e1 != win[n]->n3_e1)printf("n3_e1_%d %d\n", n, win[n]->n3_e1);
    if(window[n]->n2_e2 != win[n]->n2_e2)printf("n2_e2_%d %d\n", n, win[n]->n2_e2);
    if(window[n]->n3_e2 != win[n]->n3_e2)printf("n3_e2_%d %d\n", n, win[n]->n3_e2);
    if(window[n]->nbpt != win[n]->nbpt)printf("nbpt_%d %d\n", n, win[n]->nbpt);
    if(window[n]->nbpte1 != win[n]->nbpte1)printf("nbpte1_%d %d\n", n, win[n]->nbpte1);
    if(window[n]->nbpte1S != win[n]->nbpte1S)printf("nbpte1S_%d %d\n", n, win[n]->nbpte1S);
    if(window[n]->nbpte2 != win[n]->nbpte2)printf("nbpte2_%d %d\n", n, win[n]->nbpte2);
    if(window[n]->nbpte2S != win[n]->nbpte2S)printf("nbpte2S_%d %d\n", n, win[n]->nbpte2S);
    if(window[n]->nm2s != win[n]->nm2s)printf("nm2s_%d %d\n", n, win[n]->nm2s);
    if(window[n]->ns2m != win[n]->ns2m)printf("ns2m_%d %d\n", n, win[n]->ns2m);    

    if(window[n]->nbptS != win[n]->nbptS)printf("nbptS_%d %d\n", n, win[n]->nbptS);
    if(window[n]->nbe1 != win[n]->nbe1)printf("nbe1_%d %d\n", n, win[n]->nbe1);
    if(window[n]->nbe1S != win[n]->nbe1S)printf("nbe1S_%d %d\n", n, win[n]->nbe1S);
    if(window[n]->nbe2 != win[n]->nbe2)printf("nbe2_%d %d\n", n, win[n]->nbe2);
    if(window[n]->nbe2S != win[n]->nbe2S)printf("nbe2S_%d %d\n", n, win[n]->nbe2S);
    if(window[n]->v2_t != win[n]->v2_t)printf("v2_t_%d %d\n", n, win[n]->v2_t);
    if(window[n]->b2_t != win[n]->b2_t)printf("b2_t_%d %d\n", n, win[n]->b2_t);
    if(window[n]->a2_t != win[n]->a2_t)printf("a2_t_%d %d\n", n, win[n]->a2_t);
    if(window[n]->v3_t != win[n]->v3_t)printf("v3_t_%d %d\n", n, win[n]->v3_t);
    if(window[n]->b3_t != win[n]->b3_t)printf("b3_t_%d %d\n", n, win[n]->b3_t);
    if(window[n]->a3_t != win[n]->a3_t)printf("a3_t_%d %d\n", n, win[n]->a3_t);

    if(window[n]->v2_e1 != win[n]->v2_e1)printf("v2_e1_%d %d\n", n, win[n]->v2_e1);
    if(window[n]->b2_e1 != win[n]->b2_e1)printf("b2_e1_%d %d\n", n, win[n]->b2_e1);
    if(window[n]->a2_e1 != win[n]->a2_e1)printf("a2_e1_%d %d\n", n, win[n]->a2_e1);
    if(window[n]->x2_e1 != win[n]->x2_e1)printf("x2_e1_%d %d\n", n, win[n]->x2_e1);
    if(window[n]->v3_e1 != win[n]->v3_e1)printf("v3_e1_%d %d\n", n, win[n]->v3_e1);
    if(window[n]->b3_e1 != win[n]->b3_e1)printf("b3_e1_%d %d\n", n, win[n]->b3_e1);
    if(window[n]->a3_e1 != win[n]->a3_e1)printf("a3_e1_%d %d\n", n, win[n]->a3_e1);
    if(window[n]->x3_e1 != win[n]->x3_e1)printf("x3_e1_%d %d\n", n, win[n]->x3_e1);

    if(window[n]->v2_e2 != win[n]->v2_e2)printf("v2_e2_%d %d\n", n, win[n]->v2_e2);
    if(window[n]->b2_e2 != win[n]->b2_e2)printf("b2_e2_%d %d\n", n, win[n]->b2_e2);
    if(window[n]->a2_e2 != win[n]->a2_e2)printf("a2_e2_%d %d\n", n, win[n]->a2_e2);
    if(window[n]->x2_e2 != win[n]->x2_e2)printf("x2_e2_%d %d\n", n, win[n]->x2_e2);
    if(window[n]->v3_e2 != win[n]->v3_e2)printf("v3_e2_%d %d\n", n, win[n]->v3_e2);
    if(window[n]->b3_e2 != win[n]->b3_e2)printf("b3_e2_%d %d\n", n, win[n]->b3_e2);
    if(window[n]->a3_e2 != win[n]->a3_e2)printf("a3_e2_%d %d\n", n, win[n]->a3_e2);
    if(window[n]->x3_e2 != win[n]->x3_e2)printf("x3_e2_%d %d\n", n, win[n]->x3_e2);

    if(window[n]->nm2sS != win[n]->nm2sS)printf("nm2sS_%d %d\n", n, win[n]->nm2sS);
    if(window[n]->ns2mS != win[n]->ns2mS)printf("ns2mS_%d %d\n", n, win[n]->ns2mS);

    for (cc = 1; cc <= window[n]->a2_t; cc++) {
      if(window[n]->bot_t[cc] != win[n]->bot_t[cc])
	printf("bot_t %d : %d %d\n",cc, window[n]->bot_t[cc], win[n]->bot_t[cc]);
      if(window[n]->sur_t[cc] != win[n]->sur_t[cc])
	printf("sur_t %d : %d %d\n",cc, window[n]->sur_t[cc], win[n]->sur_t[cc]);
      if(window[n]->nsur_t[cc] != win[n]->nsur_t[cc])
	printf("nsur_t %d : %d %d\n",cc, window[n]->nsur_t[cc], win[n]->nsur_t[cc]);
    }
    for (cc = 1; cc <= window[n]->x2_e1; cc++) {
      if(window[n]->bot_e1[cc] != win[n]->bot_e1[cc])
	printf("bot_e1 %d : %d %d\n",cc, window[n]->bot_e1[cc], win[n]->bot_e1[cc]);
      if(window[n]->sur_e1[cc] != win[n]->sur_e1[cc])
	printf("sur_e1 %d : %d %d\n",cc, window[n]->sur_e1[cc], win[n]->sur_e1[cc]);
    }
    for (cc = 1; cc <= window[n]->x2_e2; cc++) {
      if(window[n]->bot_e2[cc] != win[n]->bot_e2[cc])
	printf("bot_e2 %d : %d %d\n",cc, window[n]->bot_e2[cc], win[n]->bot_e2[cc]);
      if(window[n]->sur_e2[cc] != win[n]->sur_e2[cc])
	printf("sur_e2 %d : %d %d\n",cc, window[n]->sur_e2[cc], win[n]->sur_e2[cc]);
    }

    for (cc = 1; cc <= window[n]->sgnum; cc++) {
      if(window[n]->m2d[cc] != win[n]->m2d[cc])
	printf("m2d %d : %d %d\n",cc, window[n]->m2d[cc], win[n]->m2d[cc]);
      if(window[n]->wsa[cc] != win[n]->wsa[cc])
	printf("wsa %d : %d %d\n",cc, window[n]->wsa[cc], win[n]->wsa[cc]);
      if(window[n]->xp1[cc] != win[n]->xp1[cc])
	printf("xp1 %d : %d %d\n",cc, window[n]->xp1[cc], win[n]->xp1[cc]);
      if(window[n]->xm1[cc] != win[n]->xm1[cc])
	printf("xm1 %d : %d %d\n",cc, window[n]->xm1[cc], win[n]->xm1[cc]);
      if(window[n]->yp1[cc] != win[n]->yp1[cc])
	printf("yp1 %d : %d %d\n",cc, window[n]->yp1[cc], win[n]->yp1[cc]);
      if(window[n]->ym1[cc] != win[n]->ym1[cc])
	printf("ym1 %d : %d %d\n",cc, window[n]->ym1[cc], win[n]->ym1[cc]);
      if(window[n]->zp1[cc] != win[n]->zp1[cc])
	printf("zp1 %d : %d %d\n",cc, window[n]->zp1[cc], win[n]->zp1[cc]);
      if(window[n]->xmyp1[cc] != win[n]->xmyp1[cc])
	printf("xmyp1 %d : %d %d\n",cc, window[n]->xmyp1[cc], win[n]->xmyp1[cc]);
      if(window[n]->xpym1[cc] != win[n]->xpym1[cc])
	printf("xpym1 %d : %d %d\n",cc, window[n]->xpym1[cc], win[n]->xpym1[cc]);
      if(window[n]->wgst[cc] != win[n]->wgst[cc])
	printf("wgst %d : %d %d\n",cc, window[n]->wgst[cc], win[n]->wgst[cc]);
      if(window[n]->s2i[cc] != win[n]->s2i[cc])
	printf("s2i %d : %d %d\n",cc, window[n]->s2i[cc], win[n]->s2i[cc]);
      if(window[n]->s2j[cc] != win[n]->s2j[cc])
	printf("s2j %d : %d %d\n",cc, window[n]->s2j[cc], win[n]->s2j[cc]);
      if(window[n]->s2k[cc] != win[n]->s2k[cc])
	printf("s2k %d : %d %d\n",cc, window[n]->s2k[cc], win[n]->s2k[cc]);
    }
    for (cc = 1; cc <= window[n]->nbpt; cc++) {
      if(window[n]->bpt[cc] != win[n]->bpt[cc])
	printf("bpt %d : %d %d\n",cc, window[n]->bpt[cc], win[n]->bpt[cc]);
      if(window[n]->bin[cc] != win[n]->bin[cc])
	printf("bin %d : %d %d\n",cc, window[n]->bin[cc], win[n]->bin[cc]);
      if(window[n]->bin2[cc] != win[n]->bin2[cc])
	printf("bin2 %d : %d %d\n",cc, window[n]->bin2[cc], win[n]->bin2[cc]);
    }
    for (cc = 1; cc <= window[n]->nbpte1; cc++) {
      if(window[n]->bpte1[cc] != win[n]->bpte1[cc])
	printf("bpte1 %d : %d %d\n",cc, window[n]->bpte1[cc], win[n]->bpte1[cc]);
      if(window[n]->bine1[cc] != win[n]->bine1[cc])
	printf("bine1 %d : %d %d\n",cc, window[n]->bine1[cc], win[n]->bine1[cc]);
    }
    for (cc = 1; cc <= window[n]->nbpte1S; cc++) {
      if(window[n]->bpte1S[cc] != win[n]->bpte1S[cc])
	printf("bpte1S %d : %d %d\n",cc, window[n]->bpte1S[cc], win[n]->bpte1S[cc]);
      if(window[n]->bine1S[cc] != win[n]->bine1S[cc])
	printf("bine1S %d : %d %d\n",cc, window[n]->bine1S[cc], win[n]->bine1S[cc]);
    }
    for (cc = 1; cc <= window[n]->nbpte2; cc++) {
      if(window[n]->bpte2[cc] != win[n]->bpte2[cc])
	printf("bpte2 %d : %d %d\n",cc, window[n]->bpte2[cc], win[n]->bpte2[cc]);
      if(window[n]->bine2[cc] != win[n]->bine2[cc])
	printf("bine2 %d : %d %d\n",cc, window[n]->bine2[cc], win[n]->bine2[cc]);
    }
    for (cc = 1; cc <= window[n]->nbpte2S; cc++) {
      if(window[n]->bpte2S[cc] != win[n]->bpte2S[cc])
	printf("bpte2S %d : %d %d\n",cc, window[n]->bpte2S[cc], win[n]->bpte2S[cc]);
      if(window[n]->bine2S[cc] != win[n]->bine2S[cc])
	printf("bine2S %d : %d %d\n",cc, window[n]->bine2S[cc], win[n]->bine2S[cc]);
    }
    for (cc = 1; cc <= window[n]->n2_t; cc++)
      if(window[n]->w2_t[cc] != win[n]->w2_t[cc])
	printf("w2_t %d : %d %d\n",cc, window[n]->w2_t[cc], win[n]->w2_t[cc]);
    for (cc = 1; cc <= window[n]->n3_t; cc++)
      if(window[n]->w3_t[cc] != win[n]->w3_t[cc])
	printf("w3_t %d : %d %d\n",cc, window[n]->w3_t[cc], win[n]->w3_t[cc]);
    for (cc = 1; cc <= window[n]->n2_e1; cc++)
      if(window[n]->w2_e1[cc] != win[n]->w2_e1[cc])
	printf("w2_e1 %d : %d %d\n",cc, window[n]->w2_e1[cc], win[n]->w2_e1[cc]);
    for (cc = 1; cc <= window[n]->n3_e1; cc++)
      if(window[n]->w3_e1[cc] != win[n]->w3_e1[cc])
	printf("w3_e1 %d : %d %d\n",cc, window[n]->w3_e1[cc], win[n]->w3_e1[cc]);
    for (cc = 1; cc <= window[n]->n2_e2; cc++)
      if(window[n]->w2_e2[cc] != win[n]->w2_e2[cc])
	printf("w2_e2 %d : %d %d\n",cc, window[n]->w2_e2[cc], win[n]->w2_e2[cc]);
    for (cc = 1; cc <= window[n]->n3_e2; cc++)
      if(window[n]->w3_e2[cc] != win[n]->w3_e2[cc])
	printf("w3_e2 %d : %d %d\n",cc, window[n]->w3_e2[cc], win[n]->w3_e2[cc]);

    for (cc = 1; cc <= window[n]->nm2s; cc++) {
      if(window[n]->m2s[cc] != win[n]->m2s[cc])
	printf("m2s %d : %d %d\n",cc, window[n]->m2s[cc], win[n]->m2s[cc]);
      if(window[n]->m2se1[cc] != win[n]->m2se1[cc])
	printf("m2se1 %d : %d %d\n",cc, window[n]->m2se1[cc], win[n]->m2se1[cc]);
      if(window[n]->m2se2[cc] != win[n]->m2se2[cc])
	printf("m2se2 %d : %d %d\n",cc, window[n]->m2se2[cc], win[n]->m2se2[cc]);
    }
    for (cc = 1; cc <= window[n]->ns2m; cc++) {
      if(window[n]->s2m[cc] != win[n]->s2m[cc])
	printf("s2m %d : %d %d\n",cc, window[n]->s2m[cc], win[n]->s2m[cc]);
      if(window[n]->s2me1[cc] != win[n]->s2me1[cc])
	printf("s2me1 %d : %d %d\n",cc, window[n]->s2me1[cc], win[n]->s2me1[cc]);
      if(window[n]->s2me2[cc] != win[n]->s2me2[cc])
	printf("s2me2 %d : %d %d\n",cc, window[n]->s2me2[cc], win[n]->s2me2[cc]);
    }
    printf("Window %d OK\n", n);
  }
}

/*-------------------------------------------------------------------*/
