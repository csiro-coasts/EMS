/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/outputs/df_roms.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  Output is written in ROMS compatible netCDF.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: df_roms.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#define ROMS_ALL_VARS "u ubar v vbar salt temp zeta pair sustr svstr "

/* Hard-coded vertical coordinate parameters */
#define ROMS_VTRANSFORM 2
#define ROMS_VSTRETCHING 4
#define ROMS_HC 20
/* 8 and 4 degrees in radians */
#define ROMS_THETA_S 8
#define ROMS_THETA_B 4
/* ROMS is stuck in the past */
#define ROMS_TUNIT "days since 1968-05-23"
#define ROMS_IC_VARS "zeta ubar vbar salt temp u v"
#define ROMS_BDRY_VARS "v_east v_west v_north v_south u_east u_west u_north u_south temp_east temp_west temp_north temp_south salt_east salt_west salt_north salt_south zeta_east zeta_west zeta_north zeta_south ubar_east ubar_west ubar_north ubar_south vbar_east vbar_west vbar_north vbar_south"
#define ROMS_BDRY_ORIENT {"", "east", "west", "north", "south"}

#define ROMS_NONE  0x0000
#define ROMS_EAST  0x1000
#define ROMS_WEST  0x0100
#define ROMS_NORTH 0x0010
#define ROMS_SOUTH 0x0001

/* Structure to describe each dump file time dep variable */
typedef struct {
  void **v;                     /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  nc_type type;                 /* netCDF type of this variable */
  int xylocation;
  int zlocation;
  double scale;
  int bor; /* boundary orientation */
} df_roms_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_roms_var_t *vars;           /* List of dump file variables */
} df_roms_data_t;


static void write_dump_attributes_roms(dump_data_t *dumpdata, int cdfid,
				       char *fname, int write_coord);
static void df_roms_writegeom(dump_data_t *dumpdata, dump_file_t *df, int si);
static int df_roms_get_varinfo(dump_data_t *dumpdata, dump_file_t *df, 
			     char *name, df_roms_var_t *var);
static void df_roms_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_roms_get_dimids(int cdfid, df_roms_var_t var, 
			       int *e1id, int *e2id, int *zid);
static void df_roms_get_dimids_bdry(int cdfid, df_roms_var_t var, 
				    int *e1id, int *e2id, int *zid);
static void df_roms_get_dimsizes(dump_file_t *df, df_roms_var_t var, 
				romsgrid_t *romsgrid, 
				size_t *ne1, size_t *ne2, size_t *nz);
static void df_roms_get_dimsizes_bdry(dump_file_t *df, df_roms_var_t var, 
				      romsgrid_t *romsgrid, 
				      size_t start[], size_t count[], size_t *nz);
static void tracer_write_roms_nc(int fid, int ntr, tracer_info_t tracers[], 
				 int nattr, tracer_att_t attr[], 
				 int write_coord);
static void roms_compute_depths(int N,
				float *sc, float *Cs,
				double hc, double h, double zeta,
				double *depths);
static void roms_nc_d_writesub_sigma_3d(int fid, int varid, size_t *start,
					size_t *count, double ***values,
					double *cellz, double *gridz,
					double **eta, double **h,
					float *sc, float *Cs, int nz);
static void roms_write_text_att_with_bdry(int cdfid,
					  char *str,   /* variable name */
					  char *units,
					  char *long_name,
					  char *field,
					  int write_coord,
					  char *coords);
static void roms_nc_d_writesub_sigma_bdry_1d(int fid, int varid, size_t *start,
					     size_t *count, double **values);
static void roms_nc_d_writesub_sigma_bdry_2d(int fid, int varid, size_t *start,
					     size_t *count, size_t nz, size_t nsigma,
					     double ***values,
					     double *cellz, double *gridz,
					     double **eta, double **h,
					     float *sc, float *Cs);
static int roms_get_orient(char *name);
void write_roms_bathy(int fid, int nce2, int nce1, double bmin, double **h);

static void write_dump_attributes_roms(dump_data_t *dumpdata, int cdfid,
				       char *fname, int write_coord)
{
  /* dimension ids */
  int vid;


  /* time independent variables */
  vid = ncw_var_id(cdfid, "Vtransform");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, "long_name", 
		   "vertical terrain-following transformation equation");
    vid = ncw_var_id(cdfid, "Vstretching");
    write_text_att(cdfid, vid, "long_name", 
		   "vertical terrain-following stretching function");
  }

  vid = ncw_var_id(cdfid, "theta_b");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, 
		   "long_name", "S-coordinate bottom control parameter");
    write_text_att(cdfid, vid, "units", "nondimensional");
  }

  vid = ncw_var_id(cdfid, "theta_s");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, 
		   "long_name", "S-coordinate surface control parameter");
    write_text_att(cdfid, vid, "units", "nondimensional");
  }

  vid = ncw_var_id(cdfid, "hc");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, 
		   "long_name", "S-coordinate parameter, critical depth");
    write_text_att(cdfid, vid, "units", "meter");
  }

  /* Horizontal grid */
  vid = ncw_var_id(cdfid, "h");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, "long_name", "Final bathymetry at RHO-points");
    write_text_att(cdfid, vid, "units", "meter");
    write_text_att(cdfid, vid, "field", "bath, scalar");
  }
  vid = ncw_var_id(cdfid, "lon_rho");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "longitude of RHO-points");
    write_text_att(cdfid, vid, "units", "degree_east");
    write_text_att(cdfid, vid, "standard_name", "longitude");
  }
  vid = ncw_var_id(cdfid, "lat_rho");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "latitude of RHO-points");
    write_text_att(cdfid, vid, "units", "degree_north");
    write_text_att(cdfid, vid, "standard_name", "latitude");
  }

  vid = ncw_var_id(cdfid, "lon_u");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "longitude of U-points");
    write_text_att(cdfid, vid, "units", "degree_east");
    write_text_att(cdfid, vid, "standard_name", "longitude");
  }
  vid = ncw_var_id(cdfid, "lat_u");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "latitude of U-points");
    write_text_att(cdfid, vid, "units", "degree_north");
    write_text_att(cdfid, vid, "standard_name", "latitude");
  }

  vid = ncw_var_id(cdfid, "lon_v");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "longitude of V-points");
    write_text_att(cdfid, vid, "units", "degree_east");
    write_text_att(cdfid, vid, "standard_name", "longitude");
  }
  vid = ncw_var_id(cdfid, "lat_v");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "latitude of V-points");
    write_text_att(cdfid, vid, "units", "degree_north");
    write_text_att(cdfid, vid, "standard_name", "latitude");
  }

  /* Vertical grid */
  vid = ncw_var_id(cdfid, "s_rho");
  if (vid > -1) {
    write_text_att(cdfid, vid, "long_name", "S-coordinate at RHO-points");
    write_text_att(cdfid, vid, "valid_min", "-1.");
    write_text_att(cdfid, vid, "valid_max", "0.");
    write_text_att(cdfid, vid, "positive", "up");
    if (ROMS_VTRANSFORM == 1) {
      write_text_att(cdfid, vid, "standard_name", "ocean_s_coordinate_g1");
    } else if (ROMS_VTRANSFORM == 2) {
      write_text_att(cdfid, vid, "standard_name", "ocean_s_coordinate_g2");
    } else
      hd_quit("df_roms:write_dump_attributes_roms: Illegal VTRANSFORMATION of '%d' specified\n", ROMS_VTRANSFORM);
    write_text_att(cdfid, vid, "formula_terms", 
		   "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc");
  }
  vid = ncw_var_id(cdfid, "s_w");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, "long_name", "S-coordinate at W-points");
    write_text_att(cdfid, vid, "valid_min", "-1.");
    write_text_att(cdfid, vid, "valid_max", "0.");
    write_text_att(cdfid, vid, "positive", "up");
    if (ROMS_VTRANSFORM == 1) {
      write_text_att(cdfid, vid, "standard_name", "ocean_s_coordinate_g1");
    } else if (ROMS_VTRANSFORM == 2) {
      write_text_att(cdfid, vid, "standard_name", "ocean_s_coordinate_g2");
    } else
      hd_quit("df_roms:write_dump_attributes_roms: Illegal VTRANSFORMATION of '%d' specified\n", ROMS_VTRANSFORM);
    write_text_att(cdfid, vid, "formula_terms", 
		   "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc");
  }
  
  vid = ncw_var_id(cdfid, "Cs_r");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, 
		   "long_name", "S-coordinate stretching curves at RHO-points");
    write_text_att(cdfid, vid, "units", "nondimensional");
    write_text_att(cdfid, vid, "valid_min", "-1.");
    write_text_att(cdfid, vid, "valid_max", "0.");
    write_text_att(cdfid, vid, "field", "Cs_r, scalar");
  }

  vid = ncw_var_id(cdfid, "Cs_w");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, 
		   "long_name", "S-coordinate stretching curves at W-points");
    write_text_att(cdfid, vid, "units", "nondimensional");
    write_text_att(cdfid, vid, "valid_min", "-1.");
    write_text_att(cdfid, vid, "valid_max", "0.");
    write_text_att(cdfid, vid, "field", "Cs_w, scalar");
  }
  /* time dependent variables */
  vid = ncw_var_id(cdfid, "ocean_time"); 
  if (vid > -1 ) {
    write_text_att(cdfid, vid, "long_name", "Modified Julian days");
    write_text_att(cdfid, vid, "units", ROMS_TUNIT);
    write_text_att(cdfid, vid, "field", "ocean_time, scalar, series");
  }

  /* UR added */
  vid = ncw_var_id(cdfid, "bry_time");
  if (vid > -1 ) {
    write_text_att(cdfid, vid, "long_name", "Modified Julian days");
    write_text_att(cdfid, vid, "units", ROMS_TUNIT);
    /* write_text_att(cdfid, vid, "field", "ocean_time, scalar, series"); */
  }

  roms_write_text_att_with_bdry(cdfid, "ubar", "meter second-1",
				"vertically integrated u-momentum component",
				"ubar-velocity, scalar, series",
				write_coord, "lon_u lat_u ocean_time");

  roms_write_text_att_with_bdry(cdfid, "vbar", "meter second-1",
				"vertically integrated v-momentum component",
				"vbar-velocity, scalar, series",
				write_coord, "lon_v lat_v ocean_time");
  
  roms_write_text_att_with_bdry(cdfid, "zeta", dumpdata->lenunit,
				"free surface", 
				"free-surface, scalar, series", 
				write_coord, "lon_rho lat_rho ocean_time");

  if ((vid = ncw_var_id(cdfid, "sustr")) >= 0) {
    write_text_att(cdfid, vid, "units", "N m-2");
    write_text_att(cdfid, vid, "long_name",
		   "Kinematic wind stress, u-component (m2 s-2)");
    if (write_coord)
      write_text_att(cdfid, vid, "coordinates", "lon_u lat_u ocean_time");
  }

  if ((vid = ncw_var_id(cdfid, "svstr")) >= 0) {
    write_text_att(cdfid, vid, "units", "N m-2");
    write_text_att(cdfid, vid, "long_name",
		   "Kinematic wind stress, v-component (m2 s-2)");
    if (write_coord)
      write_text_att(cdfid, vid, "coordinates", "lon_v lat_v ocean_time");
  }

  /* Units for this ?? */
  if ((vid = ncw_var_id(cdfid, "pair")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    if (write_coord)
      write_text_att(cdfid, vid, "coordinates", "lon_rho lat_rho ocean_time");
  }

  roms_write_text_att_with_bdry(cdfid, "u", "meter second-1",
				"u-momentum component",
				"u-velocity, scalar, series",
				write_coord, "lon_u lat_u ocean_time");

  roms_write_text_att_with_bdry(cdfid, "v", "meter second-1",
				"v-momentum component",
				"v-velocity, scalar, series",
				write_coord, "lon_v lat_v ocean_time");

  /* Write tracer attributes */
  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "lon_rho lat_rho s_rho ocean_time"}
    };
    tracer_write_roms_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr,
			 write_coord);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "lon_rho lat_rho ocean_time"}
    };
    tracer_write_roms_nc(cdfid, dumpdata->ntrS, dumpdata->trinfo_2d, 2, attr,
			 write_coord);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "filename", fname);
  write_text_att(cdfid, NC_GLOBAL, "type", "ROMS");
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CF-1.2");

}

void df_roms_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_roms_data_t *data = df->private_data;
  df->tout = t;
  data->nextrec = find_next_restart_record(df, data->fid, "ocean_time", t);
}

void *df_roms_create(dump_data_t *dumpdata, dump_file_t *df)
{
  romsgrid_t *romsgrid = dumpdata->romsgrid;
  int cdfid;
  int dims[4];
  int vid,n;
  /* dimension ids */
  int recdimid;
  int xi_rhoid;                 /* I dimension id at grid centre */
  int eta_rhoid;                /* J dimension id at grid centre */
  int xi_uid;                   /* I dimension id at u edges     */
  int eta_uid;                  /* J dimension id at u edges     */
  int xi_vid;                   /* I dimension id at v edges     */
  int eta_vid;                  /* J dimension id at v edges     */
  int s_rhoid;                  /* K dimension id at grid centre */
  int s_wid;                    /* K dimension id at grid layer  */
  int nfe1, nfe2;
  FILE *fp;

  int si;
  /* This flag needs to be consistent with roms_grid */
  int inface = 1;

  df_roms_data_t *data = NULL;


  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_roms_init_data(dumpdata, df, cdfid);
      data = (df_roms_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "t", dumpdata->t);
      return data;
    }
  }

  df->nz = dumpdata->nz;

  /* Set the ROMS subsection defaults */
  if (df->nce1 == dumpdata->nce1) df->nce1 = romsgrid->xi_rho;
  if (df->nce2 == dumpdata->nce2) df->nce2 = romsgrid->eta_rho;
  si = (inface) ? 1 : 0;
  nfe1 = (inface) ? romsgrid->xi_rho : romsgrid->xi_rho + 1;
  if (df->nfe1 == dumpdata->nfe1) 
    df->nfe1 = nfe1;
  else
    df->nfe1 = min(df->nfe1 - 1, nfe1);

  nfe2 = (inface) ? romsgrid->eta_rho : romsgrid->eta_rho + 1;
  if (df->nfe2 == dumpdata->nfe2) 
    df->nfe2 = nfe2;
  else
    df->nfe2 = min(df->nfe2 - 1, nfe2);

  if (df->nce1 != romsgrid->nce1 && df->ilower + df->nce1 > romsgrid->nce1)
    hd_quit("df_create_roms : Invalid i_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->ilower, df->ilower + df->nce1 - 1, romsgrid->nce1 - 1);
  if (df->nce2 != romsgrid->nce2 && df->jlower + df->nce2 > romsgrid->nce2)
    hd_quit("df_create_roms : Invalid j_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->jlower, df->jlower + df->nce2 - 1, romsgrid->nce2 - 1);

  /* create the netCDF file */
  if (ncw_create(df->name, (overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("df_roms_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "xi_rho",  romsgrid->xi_rho,  &xi_rhoid);
  nc_def_dim(cdfid, "eta_rho", romsgrid->eta_rho, &eta_rhoid);
  nc_def_dim(cdfid, "xi_u",    romsgrid->xi_u,    &xi_uid);
  nc_def_dim(cdfid, "eta_u",   romsgrid->eta_u,   &eta_uid);
  nc_def_dim(cdfid, "xi_v",    romsgrid->xi_v,    &xi_vid);
  nc_def_dim(cdfid, "eta_v",   romsgrid->eta_v,   &eta_vid);

  /* sigma depths */
  nc_def_dim(cdfid, "s_rho", romsgrid->nsigma,   &s_rhoid);
  nc_def_dim(cdfid, "s_w",   romsgrid->nsigma+1, &s_wid);
  
  /* Unlimited time dimension */
  nc_def_dim(cdfid, "ocean_time", NC_UNLIMITED, &recdimid);

  /*
   * Define variables 
   */
  // scalars
  nc_def_var(cdfid, "Vtransform", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "Vstretching", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_b", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_s", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "hc",      NC_FLOAT, 0, dims, &vid);
  
  // grid
  dims[0] = eta_rhoid;
  dims[1] = xi_rhoid;
  nc_def_var(cdfid, "h", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "lon_rho", NC_FLOAT, 2, dims, &vid);
  nc_def_var(cdfid, "lat_rho", NC_FLOAT, 2, dims, &vid);
  dims[0] = eta_uid;
  dims[1] = xi_uid;
  nc_def_var(cdfid, "lon_u", NC_FLOAT, 2, dims, &vid);
  nc_def_var(cdfid, "lat_u", NC_FLOAT, 2, dims, &vid);
  dims[0] = eta_vid;
  dims[1] = xi_vid;
  nc_def_var(cdfid, "lon_v", NC_FLOAT, 2, dims, &vid);
  nc_def_var(cdfid, "lat_v", NC_FLOAT, 2, dims, &vid);

  // vectors on rho grid
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "s_rho", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "s_w", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "Cs_r", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "Cs_w", NC_FLOAT, 1, dims, &vid);

  /* time */
  dims[0] = recdimid;
  nc_def_var(cdfid, "ocean_time", NC_DOUBLE, 1, dims, &vid);

  /* time dependent variables */
  df_roms_init_data(dumpdata, df, cdfid);
  data = (df_roms_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    if (data->vars[n].ndims == 2) {
      df_roms_get_dimids(cdfid, data->vars[n], &dims[2], &dims[1], NULL);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 3, dims, &vid);
    } else if (data->vars[n].ndims == 3) {
      df_roms_get_dimids(cdfid, data->vars[n],&dims[3],&dims[2],&dims[1]);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 4, dims, &vid);
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }
  }

  write_dump_attributes_roms(dumpdata, cdfid, df->name, 1);

  nc_enddef(cdfid);

  df_roms_writegeom(dumpdata, df, si);

  return data;
}


/*
 * Write out the actual data
 */
void df_roms_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = t;
  df_roms_data_t *data = (df_roms_data_t *)df->private_data;
  int fid = data->fid;
  romsgrid_t *romsgrid = dumpdata->romsgrid;


  /* Always cascade search */
  cs_fill_land(dumpdata);

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, ROMS_TUNIT, &newt, 1);
  nc_put_vara_double(fid, ncw_var_id(fid, "ocean_time"), start, count, &newt);

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_roms_var_t vn = data->vars[n];
    void *p = vn.v;
    double **pdn, **pd = *((double ***)p);
    if (vn.ndims == 2) {
      int j,i;
      start[1] = df->jlower;
      start[2] = df->ilower;
      df_roms_get_dimsizes(df, vn, romsgrid, &count[2], &count[1], NULL);
      pdn = d_alloc_2d(count[2]-start[2], count[1]-start[1]);
      /* Apply any scaling */
      for (j=start[1]; j<count[1]; j++)
	for (i=start[2]; i<count[2]; i++)
	  pdn[j][i] = pd[j][i] * vn.scale;
      /* Write out */
      nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, pdn);
      d_free_2d(pdn);

    } else if (vn.ndims == 3) {
      start[1] = df->klower;
      start[2] = df->jlower;
      start[3] = df->ilower;
      df_roms_get_dimsizes(df, vn, romsgrid, &count[3], &count[2], &count[1]);
      roms_nc_d_writesub_sigma_3d(fid, ncw_var_id(fid, df->vars[n]), start, count,
				  *((double ****)p), dumpdata->cellz, 
				  dumpdata->gridz, dumpdata->eta, romsgrid->h,
				  romsgrid->sc_rho, romsgrid->Cs_rho,
				  df->nz);
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

void df_roms_reset_bdry(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_roms_data_t *data = df->private_data;
  df->tout = t;
  data->nextrec = find_next_restart_record(df, data->fid, "bdry_time", t);
}

/*
 * Write out ROMS boundary data
 */
void df_roms_write_bdry(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = t;
  df_roms_data_t *data = (df_roms_data_t *)df->private_data;
  int fid = data->fid;
  romsgrid_t *romsgrid = dumpdata->romsgrid;


  /* Always cascade search */
  cs_fill_land(dumpdata);

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, ROMS_TUNIT, &newt, 1);
  nc_put_vara_double(fid, ncw_var_id(fid, "bry_time"), start, count, &newt);

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_roms_var_t vn = data->vars[n];
    void *p = vn.v;
    if (vn.ndims == 2) {
      df_roms_get_dimsizes_bdry(df, vn, romsgrid, start, count, NULL);
      roms_nc_d_writesub_sigma_bdry_1d(fid, ncw_var_id(fid, df->vars[n]), 
				       start, count, *((double ***)p));
    } else if (vn.ndims == 3) {
      size_t nsigma;
      df_roms_get_dimsizes_bdry(df, vn, romsgrid, start, count, &nsigma);
      roms_nc_d_writesub_sigma_bdry_2d(fid, ncw_var_id(fid, df->vars[n]), 
				       start, count, df->nz, nsigma,
				       *((double ****)p), dumpdata->cellz, 
				       dumpdata->gridz, dumpdata->eta, romsgrid->h,
				       romsgrid->sc_rho, romsgrid->Cs_rho);
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

/*
 * Simplified version of the ROMS output file for ocean initial conditions
 */
void df_roms_write_ic(dump_data_t *dumpdata, char *name)
{
  romsgrid_t *romsgrid = dumpdata->romsgrid;
  int cdfid;
  int dims[4];
  size_t start[4];
  size_t count[4];
  int vid,n;
  /* dimension ids */
  int recdimid;
  int xi_rhoid;                 /* I dimension id at grid centre */
  int eta_rhoid;                /* J dimension id at grid centre */
  int xi_uid;                   /* I dimension id at u edges     */
  int eta_uid;                  /* J dimension id at u edges     */
  int xi_vid;                   /* I dimension id at v edges     */
  int eta_vid;                  /* J dimension id at v edges     */
  int s_rhoid;                  /* K dimension id at grid centre */
  int s_wid;                  /* K dimension id at grid centre */
  int nfe1, nfe2;


  
  /* Dummy struct */
  dump_file_t *df = calloc(1, sizeof(dump_file_t));

  sprintf(df->name, "%s", name);

  int si;
  /* This flag needs to be consistent with roms_grid */
  int inface = 1;

  df_roms_data_t *data = NULL;


  /* Set the ROMS subsection defaults */
  if (df->nce1 == dumpdata->nce1) df->nce1 = romsgrid->xi_rho;
  if (df->nce2 == dumpdata->nce2) df->nce2 = romsgrid->eta_rho;
  si = (inface) ? 1 : 0;
  nfe1 = (inface) ? romsgrid->xi_rho : romsgrid->xi_rho + 1;
  if (df->nfe1 == dumpdata->nfe1) 
    df->nfe1 = nfe1;
  else
    df->nfe1 = min(df->nfe1 - 1, nfe1);

  nfe2 = (inface) ? romsgrid->eta_rho : romsgrid->eta_rho + 1;
  if (df->nfe2 == dumpdata->nfe2) 
    df->nfe2 = nfe2;
  else
    df->nfe2 = min(df->nfe2 - 1, nfe2);

  if (df->nce1 != romsgrid->nce1 && df->ilower + df->nce1 > romsgrid->nce1)
    hd_quit("df_create_roms : Invalid i_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->ilower, df->ilower + df->nce1 - 1, romsgrid->nce1 - 1);
  if (df->nce2 != romsgrid->nce2 && df->jlower + df->nce2 > romsgrid->nce2)
    hd_quit("df_create_roms : Invalid j_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->jlower, df->jlower + df->nce2 - 1, romsgrid->nce2 - 1);

  /* create the netCDF file */
  if (ncw_create(df->name, (overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("df_roms_create_ic: Couldn't create output dump file %s\n", 
	    df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "xi_rho",  romsgrid->xi_rho,  &xi_rhoid);
  nc_def_dim(cdfid, "eta_rho", romsgrid->eta_rho, &eta_rhoid);
  nc_def_dim(cdfid, "xi_u",    romsgrid->xi_u,    &xi_uid);
  nc_def_dim(cdfid, "eta_u",   romsgrid->eta_u,   &eta_uid);
  nc_def_dim(cdfid, "xi_v",    romsgrid->xi_v,    &xi_vid);
  nc_def_dim(cdfid, "eta_v",   romsgrid->eta_v,   &eta_vid);
  nc_def_dim(cdfid, "s_rho",   romsgrid->nsigma,  &s_rhoid);
  nc_def_dim(cdfid, "s_w",     romsgrid->nsigma+1,&s_wid);

  /* Unlimited time dimension */
  nc_def_dim(cdfid, "ocean_time", NC_UNLIMITED, &recdimid);

  /*
   * Define variables 
   */
  // UR added
  // scalars
  nc_def_var(cdfid, "Vtransform", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "Vstretching", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_b", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_s", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "hc",      NC_FLOAT, 0, dims, &vid);

  // vectors on rho grid
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "s_rho", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "s_w", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "Cs_r", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "Cs_w", NC_FLOAT, 1, dims, &vid);

  // Bottom depths
  dims[0] = eta_rhoid;
  dims[1] = xi_rhoid;
  nc_def_var(cdfid, "h", NC_DOUBLE, 2, dims, &vid);

  /* time */
  dims[0] = recdimid;
  nc_def_var(cdfid, "ocean_time", NC_DOUBLE, 1, dims, &vid);

  /* Set up variables */
  df->nvars = parseline(strdup(ROMS_IC_VARS), df->vars, MAXNUMVARS);

  /* time dependent variables */
  df_roms_init_data(dumpdata, df, cdfid);
  data = (df_roms_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    if (data->vars[n].ndims == 2) {
      df_roms_get_dimids(cdfid, data->vars[n], &dims[2], &dims[1], NULL);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 3, dims, &vid);
    } else if (data->vars[n].ndims == 3) {
      df_roms_get_dimids(cdfid, data->vars[n], &dims[3],&dims[2],&dims[1]);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 4, dims, &vid);
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'", data->vars[n].ndims);
    }
  }

  write_dump_attributes_roms(dumpdata, cdfid, df->name, 0);
  nc_enddef(cdfid);

  /* Write one time point */
  df->nz = dumpdata->nz;
  df_roms_write(dumpdata, df, df->tstart);

  /* Consts */
  {
    int ival;
    float fval;
    ival = ROMS_VTRANSFORM;
    nc_put_var_int(cdfid, ncw_var_id(cdfid, "Vtransform"), &ival);
    ival = ROMS_VSTRETCHING;
    nc_put_var_int(cdfid, ncw_var_id(cdfid, "Vstretching"), &ival);
    fval = ROMS_THETA_B;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "theta_b"), &fval);
    fval = ROMS_THETA_S;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "theta_s"), &fval);
    fval = ROMS_HC;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "hc"), &fval);
  }

  /* Vertical geometry */
  start[0] = df->klower;
  count[0] = romsgrid->nsigma;
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "s_rho"), start, count,
                    romsgrid->sc_rho);
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "Cs_r"),  start, count,
                    romsgrid->Cs_rho);

  count[0] = romsgrid->nsigma+1;
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "s_w"),  start, count,
                    romsgrid->sc_w);
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "Cs_w"), start, count,
                    romsgrid->Cs_w);

  /* Add bathymetry */
  write_roms_bathy(cdfid, romsgrid->eta_rho, 
		   romsgrid->xi_rho, 1.0, romsgrid->h);
  nc_close(cdfid);

  free(df);
}

/*
 * ROMS boundary data file
 */
void *df_roms_create_bdry(dump_data_t *dumpdata, dump_file_t *df)
{
  romsgrid_t *romsgrid = dumpdata->romsgrid;
  int cdfid;
  int dims[4];
  size_t start[4];
  size_t count[4];
  int vid,n;
  /* dimension ids */
  int recdimid;
  int xi_rhoid;                 /* I dimension id at grid centre */
  int eta_rhoid;                /* J dimension id at grid centre */
  int xi_uid;                   /* I dimension id at u edges     */
  int eta_uid;                  /* J dimension id at u edges     */
  int xi_vid;                   /* I dimension id at v edges     */
  int eta_vid;                  /* J dimension id at v edges     */
  int s_rhoid;                  /* K dimension id at grid centre */
  int s_wid;                    /* K dimension id at grid corner */
  int nfe1, nfe2;


  
  int si;
  /* This flag needs to be consistent with roms_grid */
  int inface = 1;

  df_roms_data_t *data = NULL;


  /* Set the ROMS subsection defaults */
  if (df->nce1 == dumpdata->nce1) df->nce1 = romsgrid->xi_rho;
  if (df->nce2 == dumpdata->nce2) df->nce2 = romsgrid->eta_rho;
  si = (inface) ? 1 : 0;
  nfe1 = (inface) ? romsgrid->xi_rho : romsgrid->xi_rho + 1;
  if (df->nfe1 == dumpdata->nfe1) 
    df->nfe1 = nfe1;
  else
    df->nfe1 = min(df->nfe1 - 1, nfe1);

  nfe2 = (inface) ? romsgrid->eta_rho : romsgrid->eta_rho + 1;
  if (df->nfe2 == dumpdata->nfe2) 
    df->nfe2 = nfe2;
  else
    df->nfe2 = min(df->nfe2 - 1, nfe2);

  if (df->nce1 != romsgrid->nce1 && df->ilower + df->nce1 > romsgrid->nce1)
    hd_quit("df_create_roms : Invalid i_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->ilower, df->ilower + df->nce1 - 1, romsgrid->nce1 - 1);
  if (df->nce2 != romsgrid->nce2 && df->jlower + df->nce2 > romsgrid->nce2)
    hd_quit("df_create_roms : Invalid j_range for file %s (%d-%d), valid range is 0 to %d\n", df->name, df->jlower, df->jlower + df->nce2 - 1, romsgrid->nce2 - 1);

  /* create the netCDF file */
  if (ncw_create(df->name, (overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("df_roms_create_ic: Couldn't create output dump file %s\n", 
	    df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "xi_rho",  romsgrid->xi_rho,  &xi_rhoid);
  nc_def_dim(cdfid, "eta_rho", romsgrid->eta_rho, &eta_rhoid);
  nc_def_dim(cdfid, "xi_u",    romsgrid->xi_u,    &xi_uid);
  nc_def_dim(cdfid, "eta_u",   romsgrid->eta_u,   &eta_uid);
  nc_def_dim(cdfid, "xi_v",    romsgrid->xi_v,    &xi_vid);
  nc_def_dim(cdfid, "eta_v",   romsgrid->eta_v,   &eta_vid);
  nc_def_dim(cdfid, "s_rho",   romsgrid->nsigma,  &s_rhoid);
  nc_def_dim(cdfid, "s_w",   romsgrid->nsigma+1,  &s_wid);
  
  /* Apr 20, 2014, UR; changed time dimension and ocean_time variable to bry_time as defined by Rutgers */ 
  /* Unlimited time dimension */
  nc_def_dim(cdfid, "bry_time", NC_UNLIMITED, &recdimid);

  /*
   * Define variables 
   */

  // UR added 
  // scalars
  nc_def_var(cdfid, "Vtransform", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "Vstretching", NC_INT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_b", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "theta_s", NC_FLOAT, 0, dims, &vid);
  nc_def_var(cdfid, "hc",      NC_FLOAT, 0, dims, &vid);

  // vectors on rho grid
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "s_rho", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "s_w", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_rhoid;
  nc_def_var(cdfid, "Cs_r", NC_FLOAT, 1, dims, &vid);
  dims[0] = s_wid;
  nc_def_var(cdfid, "Cs_w", NC_FLOAT, 1, dims, &vid);

  /* time */
  dims[0] = recdimid;
  nc_def_var(cdfid, "bry_time", NC_DOUBLE, 1, dims, &vid);

  /* Set up variables */
  df->nvars = parseline(strdup(ROMS_BDRY_VARS), df->vars, MAXNUMVARS);

  /* time dependent variables */
  df_roms_init_data(dumpdata, df, cdfid);
  data = (df_roms_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    // Reduced dimensions
    if (data->vars[n].ndims == 2) {
      df_roms_get_dimids_bdry(cdfid, data->vars[n], &dims[1], NULL, NULL);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 2, dims, &vid);
    } else if (data->vars[n].ndims == 3) {
      df_roms_get_dimids_bdry(cdfid, data->vars[n],&dims[2],NULL, &dims[1]);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 3, dims, &vid);

      /* Special case for T&S */
      roms_write_text_att_with_bdry(cdfid, "temp", "degrees C",
				    "Temperture", "", 0, "");
      roms_write_text_att_with_bdry(cdfid, "salt", "PSU",
				    "Salinity", "", 0, "");
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'", data->vars[n].ndims);
    }
  }

  write_dump_attributes_roms(dumpdata, cdfid, df->name, 0);
  nc_enddef(cdfid);

  /* Consts */
  {
    int ival;
    float fval;
    ival = ROMS_VTRANSFORM;
    nc_put_var_int(cdfid, ncw_var_id(cdfid, "Vtransform"), &ival);
    ival = ROMS_VSTRETCHING;
    nc_put_var_int(cdfid, ncw_var_id(cdfid, "Vstretching"), &ival);
    fval = ROMS_THETA_B;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "theta_b"), &fval);
    fval = ROMS_THETA_S;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "theta_s"), &fval);
    fval = ROMS_HC;
    nc_put_var_float(cdfid, ncw_var_id(cdfid, "hc"), &fval);
  }


  /* Vertical geometry */
  df->nz = dumpdata->nz;

  start[0] = df->klower;
  count[0] = romsgrid->nsigma;
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "s_rho"), start, count,
                    romsgrid->sc_rho);
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "Cs_r"),  start, count,
                    romsgrid->Cs_rho);

  count[0] = romsgrid->nsigma + 1;
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "s_w"),  start, count,
                    romsgrid->sc_w);
  nc_put_vara_float(cdfid, ncw_var_id(cdfid, "Cs_w"), start, count,
                    romsgrid->Cs_w);

  return(data);
}


void df_roms_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_roms_data_t *data = (df_roms_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}

static void df_roms_writegeom(dump_data_t *dumpdata, dump_file_t *df, int si)
{
  size_t start[4];
  size_t count[4];
  df_roms_data_t *data = (df_roms_data_t *)df->private_data;
  romsgrid_t *romsgrid = dumpdata->romsgrid;
  int fid = data->fid;
  int i,j;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /* Consts */
  {
    int ival;
    float fval;
    ival = ROMS_VTRANSFORM;
    nc_put_var_int(fid, ncw_var_id(fid, "Vtransform"), &ival);
    ival = ROMS_VSTRETCHING;
    nc_put_var_int(fid, ncw_var_id(fid, "Vstretching"), &ival);
    fval = ROMS_THETA_B;
    nc_put_var_float(fid, ncw_var_id(fid, "theta_b"), &fval);
    fval = ROMS_THETA_S;
    nc_put_var_float(fid, ncw_var_id(fid, "theta_s"), &fval);
    fval = ROMS_HC;
    nc_put_var_float(fid, ncw_var_id(fid, "hc"), &fval);
  }

  /* Horizontal geometry */
  if (df->long0_360) {
    // rho points
    for (j = 0; j < df->nce2; j++)
      for (i = 0; i < df->nce1; i++)
	romsgrid->lon_rho[j][i] += 180.0;
    // u-points
    for (j = 0; j < df->nce2; j++)
      for (i = si; i < df->nfe1; i++)
	romsgrid->lon_u[j][i] += 180.0;
    // v-points
    for (j = si; j < df->nfe2; j++)
      for (i = 0; i < df->nce1; i++)
	romsgrid->lon_v[j][i] += 180.0;
  }

  /* Add bathymetry */
  write_roms_bathy(fid, romsgrid->eta_rho, romsgrid->xi_rho, 
		   1.0, romsgrid->h);
  
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lon_rho"), start, count,
                   romsgrid->lon_rho);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lat_rho"), start, count,
                   romsgrid->lat_rho);

  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = romsgrid->eta_u;
  count[1] = romsgrid->xi_u;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lon_u"), start, count,
                   romsgrid->lon_u);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lat_u"), start, count,
                   romsgrid->lat_u);

  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = romsgrid->eta_v;
  count[1] = romsgrid->xi_v;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lon_v"), start, count,
                   romsgrid->lon_v);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "lat_v"), start, count,
                   romsgrid->lat_v);

  /* Vertical geometry */
  start[0] = df->klower;
  count[0] = romsgrid->nsigma;
  nc_put_vara_float(fid, ncw_var_id(fid, "s_rho"), start, count, 
		    romsgrid->sc_rho);
  nc_put_vara_float(fid, ncw_var_id(fid, "Cs_r"),  start, count,
		    romsgrid->Cs_rho);

  count[0] = romsgrid->nsigma + 1;
  nc_put_vara_float(fid, ncw_var_id(fid, "s_w"),  start, count,
		    romsgrid->sc_w);
  nc_put_vara_float(fid, ncw_var_id(fid, "Cs_w"), start, count,
		    romsgrid->Cs_w);
  
  /* Write to file */
  ncw_sync(df->name,fid);
  
  /* Undo longitude shift */
  if (df->long0_360) {
    // rho points
    for (j = 0; j < df->nce2; j++)
      for (i = 0; i < df->nce1; i++)
	romsgrid->lon_rho[j][i] -= 180.0;
    // u-points
    for (j = 0; j < df->nce2; j++)
      for (i = si; i < df->nfe1; i++)
	romsgrid->lon_u[j][i] -= 180.0;
    // v-points
    for (j = si; j < df->nfe2; j++)
      for (i = 0; i < df->nce1; i++)
	romsgrid->lon_v[j][i] -= 180.0;
  }
}



static int df_roms_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
			       char *name, df_roms_var_t *var)
{
  int found = 1;
  int n = 0;
  double fact = 4.0e3*1025.0;  /* Conversion Wm-2 to ms-1K */


  var->v = NULL;
  var->ndims = 2;
  var->type = nc_get_default_type(df->bpv);
  var->xylocation = CL_CENTRE;
  var->zlocation = CL_NONE;
  var->scale = 1.0;

  if (strncmp(name, "ubar", 4) == 0) {
    var->v = (void **)&dumpdata->u1av;
    var->xylocation = CL_GRID|CL_LEFT;
    var->bor = roms_get_orient(name);
  }

  else if (strncmp(name, "vbar", 4) == 0) {
    var->v = (void **)&dumpdata->u2av;
    var->xylocation = CL_GRID|CL_BACK;
    var->bor = roms_get_orient(name);
  }

  else if (strncmp(name, "zeta", 4) == 0) {
    var->v = (void **)&dumpdata->eta;
    var->bor = roms_get_orient(name);
  }

  else if (strcmp(name, "sustr") == 0) {
    var->v = (void **)&dumpdata->wind1;
    var->xylocation = CL_GRID|CL_LEFT;
  }

  else if (strcmp(name, "svstr") == 0) {
    var->v = (void **)&dumpdata->wind2;
    var->xylocation = CL_GRID|CL_BACK;
  }

  else if (strcmp(name, "pair") == 0) {
    var->v = (void **)&dumpdata->patm;
  }

  else if (strcmp(name, "u") == 0 ||
	   strncmp(name, "u_", 2) == 0) {
    var->v = (void **)&dumpdata->u1;
    var->ndims = 3;
    var->xylocation = CL_GRID|CL_LEFT;
    var->zlocation = CL_CENTRE;
    var->bor = roms_get_orient(name);
  }

  else if (strcmp(name, "v") == 0 ||
	   strncmp(name, "v_", 2) == 0) {
    var->v = (void **)&dumpdata->u2;
    var->ndims = 3;
    var->xylocation = CL_GRID|CL_BACK;
    var->zlocation = CL_CENTRE;
    var->bor = roms_get_orient(name);
  }

  else
    found = 0;

  if (!found) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if ( (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) ||
	   (strcmp(dumpdata->trinfo_3d[n].name, "salt") == 0 &&
	    strncmp(name, "salt_", 5) == 0) ||
	   (strcmp(dumpdata->trinfo_3d[n].name, "temp") == 0 &&
	    strncmp(name, "temp_", 5) == 0) ) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_wc[n];
        var->xylocation = CL_CENTRE;
	var->zlocation = CL_CENTRE;
	var->bor = roms_get_orient(name);
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, name) == 0) {
        var->ndims = 2;
	var->v = (void **)&dumpdata->tr_wcS[n];
        var->xylocation = CL_CENTRE;
	/* Apply any ROMS scaling */
	if (strcmp(name, "lwr") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "shf") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "lhf") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "precip") == 0) var->scale = 1000.0 * 86400.0;
	if (strcmp(name, "lprec") == 0) var->scale = 1000.0 * 86400.0;
	if (strcmp(name, "evap") == 0) var->scale = 1000.0 * 86400.0;
        found = 1;
        break;
      }
    }
  }

  return found;
}


static void df_roms_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid)
{
  int i, n;
  df_roms_data_t *data = NULL;
  char line[MAXNUMVARS * 20] = "";

  for (i = 0; i < df->nvars; ++i) {
    if (contains_token(df->vars[i], "ALL")) {
      sprintf(line, ROMS_ALL_VARS);
      /*
       * On the one hand this is good as all tracers are added
       * automatically but on the other hand extraneous unwanted ones
       * can also appear so maybe its better to maintain a list here?
       */
      for (n = 0; n < dumpdata->ntr; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);

      for (n = 0; n < dumpdata->ntrS; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_2d[n].name);
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    } else if (contains_token(df->vars[i], "NONTRACERS"))
      sprintf(line, ROMS_ALL_VARS);
    else if (contains_token(df->vars[i], "TRACERS_WC")) {
      for (n = 0; n < dumpdata->ntr; n++)
        if (!dumpdata->trinfo_3d[n].diagn)
          sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);
    } else if (contains_token(df->vars[i], "TRACERS_DIAGN_WC")) {
      for (n = 0; n < dumpdata->ntr; n++)
        if (dumpdata->trinfo_3d[n].diagn)
          sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    } else
      sprintf(line + strlen(line), "%s ", df->vars[i]);
  }

  /* Clean duplicates */
  df->nvars = parseline(strdup(line), df->vars, MAXNUMVARS);
  strcpy(line, "");
  for (i = 0; i < df->nvars; ++i)
    if (!contains_token(line, df->vars[i]))
      sprintf(line + strlen(line), "%s ", df->vars[i]);
  free(df->vars[0]);

  /* Re-parse the variable list, and check number */
  df->nvars = parseline(strdup(line), df->vars, MAXNUMVARS);

  data = (df_roms_data_t *)malloc(sizeof(df_roms_data_t));
  memset(data, 0, sizeof(df_roms_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_roms_var_t *)malloc(df->nvars * sizeof(df_roms_var_t))) == NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_roms_var_t *var = &data->vars[i];
    if (df_roms_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_roms_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

}

static void df_roms_get_dimids(int cdfid, df_roms_var_t var, 
			      int *e1id, int *e2id, int *zid)
{
  /* Horizontal grid */
  if (var.xylocation & CL_CENTRE) {
    *e1id = ncw_dim_id(cdfid, "xi_rho");
    *e2id = ncw_dim_id(cdfid, "eta_rho");
  } else if (var.xylocation & CL_GRID) {
    if (var.xylocation & CL_LEFT) {
      /* These are U-points */
      *e1id = ncw_dim_id(cdfid, "xi_u");
      *e2id = ncw_dim_id(cdfid, "eta_u");
    } else {
      /* These are V-points */
      *e1id = ncw_dim_id(cdfid, "xi_v");
      *e2id = ncw_dim_id(cdfid, "eta_v");
    }
  }

  /* We only have centres */
  if (zid) {
    *zid = ncw_dim_id(cdfid, "s_rho");
    if (!(var.zlocation & CL_CENTRE))
      hd_quit("dumpfile:df_roms_get_dimids: Only RHO-points are supported in the vertical\n");
  }
}


static void df_roms_get_dimids_bdry(int cdfid, df_roms_var_t var, 
				    int *e1id, int *e2id, int *zid)
{
  /* Horizontal grid */
  if (var.xylocation & CL_CENTRE) {
    if (var.bor & (ROMS_EAST|ROMS_WEST))
      *e1id = ncw_dim_id(cdfid, "eta_rho");
    else
      *e1id = ncw_dim_id(cdfid, "xi_rho");
  } else if (var.xylocation & CL_GRID) {
    if (var.xylocation & CL_LEFT) {
      /* These are U-points */
      if (var.bor & (ROMS_EAST|ROMS_WEST))
	*e1id = ncw_dim_id(cdfid, "eta_u");
      else
	*e1id = ncw_dim_id(cdfid, "xi_u");
    } else {
      /* These are V-points */
      if (var.bor & (ROMS_EAST|ROMS_WEST))
	*e1id = ncw_dim_id(cdfid, "eta_v");
      else
	*e1id = ncw_dim_id(cdfid, "xi_v");
    }
  }

  /* We only have centres */
  if (zid) {
    *zid = ncw_dim_id(cdfid, "s_rho");
    if (!(var.zlocation & CL_CENTRE))
      hd_quit("dumpfile:df_roms_get_dimids: Only RHO-points are supported in the vertical\n");
  }

}


static void df_roms_get_dimsizes(dump_file_t *df, df_roms_var_t var,
				 romsgrid_t *romsgrid,
				 size_t *ne1, size_t *ne2, size_t *nz)
{
  if (var.xylocation & CL_CENTRE) {
    *ne1 = romsgrid->xi_rho;
    *ne2 = romsgrid->eta_rho;
  } else if (var.xylocation & CL_GRID) {
    if (var.xylocation & CL_LEFT) {
      /* These are U-points */
      *ne1 = romsgrid->xi_u;
      *ne2 = romsgrid->eta_u;
    } else {
      /* These are V-points */
      *ne1 = romsgrid->xi_v;
      *ne2 = romsgrid->eta_v;
    }
  }

  /* We only have centres */
  if (nz) {
    *nz = romsgrid->nsigma;
    if (!(var.zlocation & CL_CENTRE))
      hd_quit("dumpfile:df_roms_get_dimsizes: Only RHO-points are supported in the vertical\n");
  }
}

static void df_roms_get_dimsizes_bdry(dump_file_t *df, df_roms_var_t var,
				      romsgrid_t *romsgrid,
				      size_t start[], size_t count[], size_t *nz)
{
  /* 
   * Fill in start/count J,I (ie. eta,xi) 
   *
   * Also note that first index is always time, if this assumption
   * changes we can make it an offset, similar to other places
   */
  if (var.xylocation & CL_CENTRE) {
    if (var.bor & ROMS_EAST) {
      start[1] = 0;
      count[1] = romsgrid->eta_rho;
      start[2] = romsgrid->xi_rho - 1;
      count[2] = 1;
    }
    else if (var.bor & ROMS_WEST) {
      start[1] = 0;
      count[1] = romsgrid->eta_rho;
      start[2] = 0;
      count[2] = 1;
    }
    else if (var.bor & ROMS_NORTH) {
      start[1] = romsgrid->eta_rho - 1;
      count[1] = 1;
      start[2] = 0;
      count[2] = romsgrid->xi_rho;
    } else {
      /* South */
      start[1] = 0;
      count[1] = 1;
      start[2] = 0;
      count[2] = romsgrid->xi_rho;
    }
  } else if (var.xylocation & CL_GRID) {
    if (var.xylocation & CL_LEFT) {
      /* These are U-points */
      if (var.bor & ROMS_EAST) {
	start[1] = 0;
	count[1] = romsgrid->eta_u;
	start[2] = romsgrid->xi_u - 1;
	count[2] = 1;
      }
      else if (var.bor & ROMS_WEST) {
	start[1] = 0;
	count[1] = romsgrid->eta_u;
	start[2] = 0;
	count[2] = 1;
      }
      else if (var.bor & ROMS_NORTH) {
	start[1] = romsgrid->eta_u - 1;
	count[1] = 1;
	start[2] = 0;
	count[2] = romsgrid->xi_u;
      } else {
	/* South */
	start[1] = 0;
	count[1] = 1;
	start[2] = 0;
	count[2] = romsgrid->xi_u;
      }
    } else {
      /* These are V-points */
      if (var.bor & ROMS_EAST) {
	start[1] = 0;
	count[1] = romsgrid->eta_v;
	start[2] = romsgrid->xi_v - 1;
	count[2] = 1;
      }
      else if (var.bor & ROMS_WEST) {
	start[1] = 0;
	count[1] = romsgrid->eta_v;
	start[2] = 0;
	count[2] = 1;
      }
      else if (var.bor & ROMS_NORTH) {
	start[1] = romsgrid->eta_v - 1;
	count[1] = 1;
	start[2] = 0;
	count[2] = romsgrid->xi_v;
      } else {
	/* South */
	start[1] = 0;
	count[1] = 1;
	start[2] = 0;
	count[2] = romsgrid->xi_v;
      }
    }
  }

  /* We only have centres */
  if (nz) {
    *nz = romsgrid->nsigma;
    if (!(var.zlocation & CL_CENTRE))
      hd_quit("dumpfile:df_roms_get_dimsizes: Only RHO-points are supported in the vertical\n");
  }
}

/*
 * See https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
 * for more details. There are several different stretching functions
 * but I have only implemented 1 for now. SHOC also has a function
 * so it might be worth consolidating these so that we have both the
 * ability to dump the native sigma coordinates when SHOC runs in
 * sigma and for SHOC to use different vertical stretching depending
 * on these new stretching schemes
 */
void roms_compute_vert(int N, float *sc, float *Cs, int isW)
{
  int k=0;
  double cs1;

  /* Bail out if not A. Shchepetkin (2010) UCLA-ROMS */
  if (ROMS_VSTRETCHING != 4)
    hd_quit("roms_compute_vert: Only stretching 4 is currently supported\n");

  if (isW) {
    /* W-points */
    for (k=0; k<N; k++) {
      sc[k] = (k-N-0.)/N; // force floating point calculation
      cs1 = (1 - cosh(ROMS_THETA_S * sc[k])) / (cosh(ROMS_THETA_S) - 1);
      Cs[k] = (exp(ROMS_THETA_B * cs1) - 1) / (1 - exp(-ROMS_THETA_B));
    }
  } else {
    /* RHO-points */
    for (k=1; k<=N; k++) {
      sc[k-1] = (k-N-0.5)/N;
      cs1 = (1 - cosh(ROMS_THETA_S * sc[k-1])) / (cosh(ROMS_THETA_S) - 1);
      Cs[k-1] = (exp(ROMS_THETA_B * cs1) - 1) / (1 - exp(-ROMS_THETA_B));
    }
  }
}


/*
 *
 * See https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
 * for more details. There are two transforms referenced there but
 * there is only one version (2) that's implemented here
 *
 * Given the vector of S computes the corresponding depths
 */
static void roms_compute_depths(int N,
				float *sc, float *Cs,
				double hc, double h, double zeta,
				double *depths) // out
{
  double S;
  int i;

  /* Bail out if not using Transform 2 */
  if (ROMS_VTRANSFORM != 2)
    hd_quit("roms_compute_depths: Only Transform 2 is currently supported\n");
  
  for (i=0; i<N; i++) {
    S = (hc*sc[i] + h*Cs[i]) / (hc + h);
    depths[i] = zeta + (zeta + h) * S;
  }
}

/** Writes tracer attributes to a NetCDF file.
 * `attr' and `attrval' are common attributes added to all tracers.
 * @param fid NetCDF file
 * @param ntr number of tracers
 * @param tracers array of `tracer_desc'
 * @param attr array of additional attributes to be added to all tracers
 */
static void tracer_write_roms_nc(int fid, int ntr, tracer_info_t tracers[],
				 int nattr, tracer_att_t attr[], 
				 int write_coord)
{
  int n;
  char name[MAXSTRLEN];
  char lname[MAXSTRLEN];
  char units[MAXSTRLEN];

  for (n = 0; n < ntr; n++) {
    int vid,ai;
    tracer_info_t *tr = &tracers[n];
    strcpy(name,tr->name);
    if (strcmp(attr[0].attr,"tracer_sed") == 0)
      strcat(name, "_sed");
    strcpy(lname, tr->long_name);
    strcpy(units, tr->units);
    /*
      if (strcmp(name, "lhf") == 0) {
      strcpy(lname, "Evaporation rate");
      strcpy(units, "kg/m2/s");
      }
    */
    if (strcmp(name, "evap") == 0) strcpy(units, "mm day-1");
    if (strcmp(name, "precip") == 0) strcpy(units, "mm day-1");
    if (strcmp(name, "lprec") == 0) strcpy(units, "mm day-1");
    vid = ncw_var_id(fid, name);
    if (vid >= 0) {
      nc_put_att_text(fid, vid, "long_name", strlen(lname), lname);
      nc_put_att_text(fid, vid, "units", strlen(units), units);
      nc_put_att_double(fid, vid, "_FillValueWC", NC_DOUBLE, 1,
                        &tr->fill_value_wc);
      nc_put_att_double(fid, vid, "valid_range", NC_DOUBLE, 2,
                        tr->valid_range_wc);
      for (ai=0; ai<nattr; ai++) {
	if ((strcmp(attr[ai].attr, "coordinates") == 0) && !write_coord)
	  continue;
	nc_put_att_text(fid, vid, attr[ai].attr, 
			strlen(attr[ai].attrval), attr[ai].attrval);
      }
    }
  }
}


/* 
 * Write out a double 3d subsection in sigma coordinates
 * 
 * The tracers are defined on cell centres but the sigma levels can
 * end up outside these points at the top and bottom. Not sure exactly
 * what the best solution is but this functions simply holds the
 * values at the top and bottom
 */
static void roms_nc_d_writesub_sigma_3d(int fid, int varid, size_t *start,
					size_t *count, double ***values,
					double *cellz, double *gridz, 
					double **eta, double **h,
					float *sc, float *Cs, int nz)
{
  int nd = 0;
  unsigned int offset = 0;  /* Has a record dimension = 1, else 0 */
  double ***nvals = NULL;
  double *z_depths = NULL, *s_depths = NULL;
  double *z_vals = NULL, *s_vals = NULL;
  size_t nstart[4];
  unsigned int i, j, k, I, J, K;
  int status;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 3);

  I = count[offset + 2];  /* number of cells in the i-direction */
  J = count[offset + 1];  /* number of cells in the j-direction */
  K = count[offset];      /* number of cells in the k-direction */

  /* Final 3D array */
  nvals = d_alloc_3d(I, J, K);

  /* column based work arrays */
  z_depths  = d_alloc_1d(nz);
  s_depths  = d_alloc_1d(K);
  z_vals    = d_alloc_1d(nz);
  s_vals    = d_alloc_1d(K);
  /*
   * Loop over each column
   */
  for (j = 0; j < J; ++j)
    for (i = 0; i < I; ++i) {
      int ncells = 0;
      int jidx = start[offset + 1] + j;
      int iidx = start[offset + 2] + i;

      /* Get sigma depths */
      roms_compute_depths(K, sc, Cs, ROMS_HC, h[j][i], eta[j][i], s_depths);

      /* Fill in all cell centres in the water column */
      memset(z_vals,   0, nz*sizeof(double));
      memset(z_depths, 0, nz*sizeof(double));
      for (k = 0; k < nz; k++) {
	/* 
	 * Commented out so that we interpolate on land as well
	 */
	//if (-h[jidx][iidx] <= cellz[k]) {
	z_vals[ncells]   = values[start[offset] + k][jidx][iidx];
	z_depths[ncells] = cellz[k];
	ncells++;
	//}
      }
      
      if (ncells) {
	/*
	 * Interpolate onto sigma depths
	 * Note: This function automatically applies a ZOH above and below
	 */
	interp1d(z_depths, z_vals, ncells, s_depths, s_vals, K);

      } else {
	/* Land or outside */
	memcpy(s_vals, z_vals, K*sizeof(double));
      }
      
      /* Fill this column */
      for (k = 0; k < K; ++k)
	nvals[k][j][i] = s_vals[k];
    }
  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;
  
  status = nc_put_vara_double(fid, varid, nstart, count, nvals[0][0]);
  if (status != NC_NOERR)
    hd_quit("df_roms:roms_nc_d_writesub_sigma_3d: %s\n", nc_strerror(status));
  
  d_free_3d(nvals);
  d_free_1d(z_depths);
  d_free_1d(s_depths);
  d_free_1d(z_vals);
  d_free_1d(s_vals);
}

/*
 * Extracts and write out a single column/row vector out of the given 2D variable
 */
static void roms_nc_d_writesub_sigma_bdry_1d(int fid, int varid, size_t *start,
					     size_t *count, double **values)
{
  double *nvals = NULL;
  size_t nstart[4], ncount[4];
  int status, offset;
  int i, I, j, J, l, len;

  offset = 1; /* Always TIME varying */

  I = count[offset + 1];  /* number of cells in the i-direction */
  J = count[offset];      /* number of cells in the j-direction */

  /* Length of this array */
  len = max(I, J);

  /* Final 1D array */
  nvals = d_alloc_1d(len);

  /*
   * Loop over 2D array and pluck out the appropriate column/row
   */
  l = 0;
  for (j = start[offset]; j < start[offset] + J; ++j)
    for (i = start[offset+1]; i < start[offset+1] + I; ++i)
      /* Fill this vector */
      nvals[l++] = values[j][i];
  
  /* Time */
  nstart[0] = start[0]; ncount[0] = 1L;
  /* Vector */
  nstart[1] = 0; ncount[1] = len;
 
  /* Write data */
  status = nc_put_vara_double(fid, varid, nstart, ncount, nvals);
  if (status != NC_NOERR)
    hd_quit("df_roms:roms_nc_d_writesub_sigma_bdry_1d: %s\n", nc_strerror(status));
  
  d_free_1d(nvals);
}


/* 
 * Extract and write out a single row/column vertical plane out of the
 * given 3D variable
 */
static void roms_nc_d_writesub_sigma_bdry_2d(int fid, int varid, size_t *start,
					     size_t *count, size_t nz,
					     size_t nsigma,
					     double ***values,
					     double *cellz, double *gridz,
					     double **eta, double **h,
					     float *sc, float *Cs)
{
  double **nvals = NULL;
  size_t nstart[4], ncount[4];
  int status, offset;
  int i, I, j, J, l, len, k, K;
  double *z_depths = NULL, *s_depths = NULL;
  double *z_vals = NULL, *s_vals = NULL;

  offset = 1; /* Always TIME varying */

  I = count[offset + 1];  /* number of cells in the i-direction */
  J = count[offset];      /* number of cells in the j-direction */
  K = nsigma;             /* number of cells in the k-direction */

  /* Length of this array */
  len = max(I, J);

  /* Final 2D array */
  nvals = d_alloc_2d(len,K);

  /* column based work arrays */
  z_depths = d_alloc_1d(nz);
  s_depths = d_alloc_1d(K);
  z_vals   = d_alloc_1d(nz);
  s_vals   = d_alloc_1d(K);

  /*
   * Loop over 2D array and pluck out the appropriate column/row
   */
  l = 0;
  for (j = start[offset]; j < start[offset] + J; ++j) {
    for (i = start[offset+1]; i < start[offset+1] + I; ++i) {
      int ncells = 0;

      /* Get sigma depths */
      roms_compute_depths(K, sc, Cs, ROMS_HC, h[j][i], eta[j][i], s_depths);
      
      /* Fill in all cell centres in the water column */
      memset(z_vals,   0, nz*sizeof(double));
      memset(z_depths, 0, nz*sizeof(double));
      for (k = 0; k < nz; k++) {
	/* 
	 * Commented out so that we interpolate on land as well
	 */
	// if (-h[j][i] <= cellz[k]) {
	z_vals[ncells]   = values[k][j][i];
	z_depths[ncells] = cellz[k];
	ncells++;
	//}
      }

      if (ncells) {
	/*
	 * Interpolate onto sigma depths
	 * Note: This function automatically applies a ZOH above and below
	 */
	interp1d(z_depths, z_vals, ncells, s_depths, s_vals, K);

      } else {
	/* Land or outside */
	memcpy(s_vals, z_vals, K*sizeof(double));
      }
      
      /* Fill this column */
      for (k = 0; k < K; ++k)
	nvals[k][l] = s_vals[k];
      l++;
    }
  }
  
  /* Time */
  nstart[0] = start[0]; ncount[0] = 1L;
  /* Depth */
  nstart[1] = 0; ncount[1] = K;
  /* Row or column length */
  nstart[2] = 0; ncount[2] = len;

  /* Write data */
  status = nc_put_vara_double(fid, varid, nstart, ncount, nvals[0]);
  if (status != NC_NOERR)
    hd_quit("df_roms:roms_nc_d_writesub_sigma_bdry_2d: %s\n", nc_strerror(status));
  
  d_free_2d(nvals);
  d_free_1d(z_depths);
  d_free_1d(s_depths);
  d_free_1d(z_vals);
  d_free_1d(s_vals);
}


static void roms_write_text_att_with_bdry(int cdfid,
					  char *str,   /* variable name */
					  char *units,
					  char *long_name,
					  char *field,
					  int write_coord,
					  char *coords)
{
  char *orient[] = ROMS_BDRY_ORIENT;
  int i, vid;
  char buf[MAXSTRLEN];
  
  // Base
  for (i=0; i<5; i++) {
    if (strlen(orient[i]))
      sprintf(buf, "%s_%s", str, orient[i]);
    else
      sprintf(buf, "%s", str);
    if ( (vid = ncw_var_id(cdfid, buf)) >= 0 ) {
      if (strlen(units))
	write_text_att(cdfid, vid, "units",     units);
      if (strlen(orient[i])) {
	sprintf(buf, "%s %s boundary", long_name, orient[i]);
	write_text_att(cdfid, vid, "long_name", buf);
      }
      else
	write_text_att(cdfid, vid, "long_name", long_name);
      if (strlen(field))
	write_text_att(cdfid, vid, "field",     field);
      if (write_coord)
	write_text_att(cdfid, vid, "coordinates", coords);
      /* Point all time to the same time variable for boundaries only */
      if (i)
	write_text_att(cdfid, vid, "time", "bry_time");
    }
  }
}

static int roms_get_orient(char *name)
{
  if (strstr(name, "_east"))
    return(ROMS_EAST);
  else if (strstr(name, "_west"))
    return(ROMS_WEST);
  else if (strstr(name, "_north"))
    return(ROMS_NORTH);
  else if (strstr(name, "_south"))
    return(ROMS_SOUTH);
  else
    return(ROMS_NONE);
}

/* Write bathymetry, with land = bmin */
void write_roms_bathy(int fid, int nce2, int nce1, double bmin, double **h)
{
  size_t start[2] = {0, 0}, count[2];
  int j,i;
  count[0] = nce2;
  count[1] = nce1;
  double **hnew = d_alloc_2d(nce1, nce2);
  for (j=0; j<nce2; j++) {
    for (i=0; i<nce1; i++) {
      if (isnan(h[j][i]))
	hnew[j][i] = bmin;
      else
	hnew[j][i] = h[j][i];
    }
  }
  nc_put_vara_double(fid, ncw_var_id(fid, "h"), start, count, hnew[0]);

  d_free_2d(hnew);
}
