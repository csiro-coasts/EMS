/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/df_mom.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  Output is written in MOM compatible netCDF.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: df_mom.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define MOMID   1
#define EMSID   2
#define NOID    4

#define MOM_ALL_VARS "u1av u2av wtop eta_t taux tauy patm u1 u2 w dens dens_0 Kz Vz Cd u1vh u2vh "

/* Structure to describe each dump file time dep variable */
typedef struct {
  void **v;                     /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  nc_type type;                 /* netCDF type of this variable */
  int xylocation;
  int zlocation;
  double scale;
} df_mom_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_mom_var_t *vars;           /* List of dump file variables */
} df_mom_data_t;


void write_dump_attributes_mom(dump_data_t *dumpdata, int cdfid,
			       nc_type fptype, int nce1, int nce2, int nz, 
			       char *fname);
static void df_mom_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static int df_mom_get_varinfo(dump_data_t *dumpdata, dump_file_t *df, 
			     char *name, df_mom_var_t *var);
static void df_mom_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid);
static void df_mom_get_dimids(int cdfid, df_mom_var_t var, 
			      int *e1id, int *e2id, int *zid);
static void df_mom_get_dimsizes(dump_file_t *df, df_mom_var_t var, 
				momgrid_t *momgrid, 
				size_t *ne1, size_t *ne2, size_t *nz);
void pack_mom2d(int gid, dump_file_t *df, int nx, int ny, int nz,
		double **in, double **out, double sc);
void pack_mom3d(int gid, dump_file_t *df, int nx, int ny, int nz, 
		double ***in, double ***out, double sc);
void tracer_write_mom_nc(int fid, int ntr, tracer_info_t tracers[], int nattr,
			 tracer_att_t attr[]);


void write_dump_attributes_mom(dump_data_t *dumpdata, int cdfid,
			       nc_type fptype, int nce1, int nce2, int nz, 
			       char *fname)
			      
			      
{
  /* dimension ids */
  int vid, n;
  int coid = MOMID; /* coid = MOMID : MOM type coordinate definition */
                    /* coid = EMSID : EMS type coordinate definition */
                    /* coid = NOID : No coordinate definition */
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];

  /* time independent variables */
  vid = ncw_var_id(cdfid, "grid_x_T");
  write_text_att(cdfid, vid, "long_name", "Nominal Longitude of T-cell center");
  write_text_att(cdfid, vid, "units", "degree_east");
  write_text_att(cdfid, vid, "cartesian_axis", "X");
  vid = ncw_var_id(cdfid, "grid_y_T");
  write_text_att(cdfid, vid, "long_name", "Nominal Latitude of T-cell center");
  write_text_att(cdfid, vid, "units", "degree_north");
  write_text_att(cdfid, vid, "cartesian_axis", "Y");
  vid = ncw_var_id(cdfid, "grid_x_C");
  write_text_att(cdfid, vid, "long_name", "Nominal Longitude of C-cell center");
  write_text_att(cdfid, vid, "units", "degree_east");
  write_text_att(cdfid, vid, "cartesian_axis", "X");
  vid = ncw_var_id(cdfid, "grid_y_C");
  write_text_att(cdfid, vid, "long_name", "Nominal Latitude of C-cell center");
  write_text_att(cdfid, vid, "units", "degree_north");
  write_text_att(cdfid, vid, "cartesian_axis", "Y");

  /* Vertical grid */
  vid = ncw_var_id(cdfid, "zt");
  write_text_att(cdfid, vid, "long_name", "zt");
  write_text_att(cdfid, vid, "units", "meters");
  write_text_att(cdfid, vid, "cartesian_axis", "z");
  write_text_att(cdfid, vid, "positive", "down");
  vid = ncw_var_id(cdfid, "zb");
  write_text_att(cdfid, vid, "long_name", "zb");
  write_text_att(cdfid, vid, "units", "meters");
  write_text_att(cdfid, vid, "cartesian_axis", "z");
  write_text_att(cdfid, vid, "positive", "down");

  /* T grid cells */
  vid = ncw_var_id(cdfid, "x_T");
  write_text_att(cdfid, vid, "long_name", "Geographic longitude of T_cell centers");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "y_T");
  write_text_att(cdfid, vid, "long_name", "Geographic latitude of T_cell centers");
  write_text_att(cdfid, vid, "units", "degree_north");

  /* C grid cells */
  vid = ncw_var_id(cdfid, "x_C");
  write_text_att(cdfid, vid, "long_name", "Geographic longitude of C_cell centers");
  write_text_att(cdfid, vid, "units", "degree_east");
  vid = ncw_var_id(cdfid, "y_C");
  write_text_att(cdfid, vid, "long_name", "Geographic latitude of C_cell centers");
  write_text_att(cdfid, vid, "units", "degree_north");


  /* time dependent variables */
  vid = ncw_var_id(cdfid, "TIME");
  write_text_att(cdfid, vid, "long_name", "TIME");
  /* Strip the time zone */
  strcpy(buf,dumpdata->output_tunit);
  n = parseline(buf, fields, MAXNUMARGS);
  sprintf(buf,"%s %s %s %s",fields[0],fields[1],fields[2],fields[3]);
  write_text_att(cdfid, vid, "units", buf);
  write_text_att(cdfid, vid, "cartesian_axis", "T");
  write_text_att(cdfid, vid, "calendar", "julian");

  if ((vid = ncw_var_id(cdfid, "u1av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "I component of depth averaged current at left face");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C");
  }

  if ((vid = ncw_var_id(cdfid, "u2av")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "J component of depth averaged current at back face");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C");
  }

  if ((vid = ncw_var_id(cdfid, "wtop")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "Vertical velocity at surface");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_T, y_T");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_T, grid_y_T");
  }

  if ((vid = ncw_var_id(cdfid, "eta_t")) >= 0) {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "surface height on T cells");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_T, y_T");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_T, grid_y_T");
  }

  if ((vid = ncw_var_id(cdfid, "taux")) >= 0) {
    write_text_att(cdfid, vid, "units", "N/m^2");
    write_text_att(cdfid, vid, "long_name", "ZONAL WIND STRESS");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C");
  }

  if ((vid = ncw_var_id(cdfid, "tauy")) >= 0) {
    write_text_att(cdfid, vid, "units", "N/m^2");
    write_text_att(cdfid, vid, "long_name", "MERIDIONAL WIND STRESS");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C");
  }

  if ((vid = ncw_var_id(cdfid, "patm")) >= 0) {
    write_text_att(cdfid, vid, "units", "Pa");
    write_text_att(cdfid, vid, "long_name", "Atmospheric pressure");
    write_text_att(cdfid, vid, "coordinates", "t, grid_x_T, grid_y_T");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_T, y_T");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_T, grid_y_T");
  }

  if ((vid = ncw_var_id(cdfid, "u1")) >= 0) {
    write_text_att(cdfid, vid, "units", "m/sec");
    write_text_att(cdfid, vid, "long_name", "zonal current");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
  }

  if ((vid = ncw_var_id(cdfid, "u2")) >= 0) {
    write_text_att(cdfid, vid, "units", "m/sec");
    write_text_att(cdfid, vid, "long_name", "meridional current");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
  }

  if ((vid = ncw_var_id(cdfid, "w")) >= 0) {
    write_text_att(cdfid, vid, "units", "metre second-1");
    write_text_att(cdfid, vid, "long_name",
                   "K component of current at cell centre and Z grid");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zb");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zb");
  }

  if ((vid = ncw_var_id(cdfid, "dens")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Density");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
  }

  if ((vid = ncw_var_id(cdfid, "dens_0")) >= 0) {
    write_text_att(cdfid, vid, "units", "kg metre-3");
    write_text_att(cdfid, vid, "long_name", "Potential density");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
  }

  if ((vid = ncw_var_id(cdfid, "Kz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Kz");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
  }

  if ((vid = ncw_var_id(cdfid, "Vz")) >= 0) {
    write_text_att(cdfid, vid, "units", "m2 s-1");
    write_text_att(cdfid, vid, "long_name", "Vz");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C, zt");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C, zt");
 }

  if ((vid = ncw_var_id(cdfid, "Cd")) >= 0) {
    write_text_att(cdfid, vid, "units", " ");
    write_text_att(cdfid, vid, "long_name", "Bottom drag coefficient");
    if (coid == EMSID)
      write_text_att(cdfid, vid, "coordinates", "t, x_C, y_C");
    else if (coid == MOMID)
      write_text_att(cdfid, vid, "coordinates", "t, grid_x_C, grid_y_C");
  }

  {
    tracer_att_t attr[] = { {"tracer", "true"},
    {"coordinates", "t, grid_x_C, grid_y_C, zt"}
    };
    tracer_write_mom_nc(cdfid, dumpdata->ntr, dumpdata->trinfo_3d, 2, attr);
  }

  {
    tracer_att_t attr[] = { {"tracer2D", "true"},
    {"coordinates", "t, grid_x_C, grid_y_C"}
    };
    tracer_write_mom_nc(cdfid, dumpdata->ntrS, dumpdata->trinfo_2d, 2, attr);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "filename", fname);

}

void df_mom_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_mom_data_t *data = df->private_data;
  df->tout = t;
  data->nextrec = find_next_restart_record(df, data->fid, "t", t);
}

void *df_mom_create(dump_data_t *dumpdata, dump_file_t *df)
{
  momgrid_t *momgrid = dumpdata->momgrid;
  int n;
  int cdfid;
  int dims[10];
  int vid;                      /* df_variable_t dimension */
  /* dimension ids */
  int recdimid;                 /* record dimension id */
  int grid_x_C;                 /* I dimension id at grid corner */
  int grid_y_C;                 /* J dimension id at grid corner */
  int grid_x_T;                 /* I dimension id at grid centre */
  int grid_y_T;                 /* J dimension id at grid centre */
  int zt;                       /* K dimension id at grid centre */
  int zb;                       /* K dimension id at grid face   */
  FILE *fp;

  df_mom_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_mom_init_data(dumpdata, df, cdfid);
      data = (df_mom_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "t", dumpdata->t);
      return data;
    }
  }

  /* Set the MOM subsection defaults */
  if (df->nce1 == dumpdata->nce1) df->nce1 = momgrid->grid_x_T;
  if (df->nce2 == dumpdata->nce2) df->nce2 = momgrid->grid_y_T;
  if (df->nfe1 == dumpdata->nfe1) 
    df->nfe1 = momgrid->grid_x_C;
  else
    df->nfe1 = min(df->nfe1 - 1, momgrid->grid_x_C);

  if (df->nfe2 == dumpdata->nfe2)
    df->nfe2 = momgrid->grid_y_C;
  else
    df->nfe2 = min(df->nfe2 - 1, momgrid->grid_y_C);
  /*
  if (momgrid->cytype & U1BDRY) {
    if (df->nce1 == momgrid->grid_x_T) df->nce1 -= 1;
    if (df->nfe1 == momgrid->grid_x_C) df->nfe1 -= 1;
  }
  if (momgrid->cytype & U2BDRY) {
    if (df->nce2 == momgrid->grid_y_T) df->nce2 -= 1;
    if (df->nfe2 == momgrid->grid_y_C) df->nfe2 -= 1;
  }
  */

  if (df->nce1 != momgrid->nce1 && df->ilower + df->nce1 > momgrid->nce1)
    hd_quit("df_create_mom : Invalid i_range for file %s (%d-%d), valid range is 0 to %d\n", 
	    df->name, df->ilower, df->ilower + df->nce1 - 1, momgrid->grid_x_T - 1);
  if (df->nce2 != momgrid->nce2 && df->jlower + df->nce2 > momgrid->nce2)
    hd_quit("df_create_mom : Invalid j_range for file %s (%d-%d), valid range is 0 to %d\n", 
	    df->name, df->jlower, df->jlower + df->nce2 - 1, momgrid->grid_y_T - 1);

  /* create the netCDF file */
  if (ncw_create(df->name, (overwrite_output()?NC_CLOBBER:NC_NOCLOBBER), &cdfid) != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "grid_x_T", df->nce1, &grid_x_T);
  nc_def_dim(cdfid, "grid_y_T", df->nce2, &grid_y_T);
  nc_def_dim(cdfid, "grid_x_C", df->nfe1, &grid_x_C);
  nc_def_dim(cdfid, "grid_y_C", df->nfe2, &grid_y_C);
  nc_def_dim(cdfid, "zt", df->nz, &zt);
  nc_def_dim(cdfid, "zb", df->nz, &zb);
  nc_def_dim(cdfid, "TIME", NC_UNLIMITED, &recdimid);
  dims[0] = grid_x_T;
  nc_def_var(cdfid, "grid_x_T", NC_FLOAT, 1, dims, &vid);
  dims[0] = grid_y_T;
  nc_def_var(cdfid, "grid_y_T", NC_FLOAT, 1, dims, &vid);
  dims[0] = grid_x_C;
  nc_def_var(cdfid, "grid_x_C", NC_FLOAT, 1, dims, &vid);
  dims[0] = grid_y_C;
  nc_def_var(cdfid, "grid_y_C", NC_FLOAT, 1, dims, &vid);
  dims[0] = zt;
  nc_def_var(cdfid, "zt", NC_FLOAT, 1, dims, &vid);
  dims[0] = zb;
  nc_def_var(cdfid, "zb", NC_FLOAT, 1, dims, &vid);

  dims[0] = grid_y_T;
  dims[1] = grid_x_T;
  nc_def_var(cdfid, "x_T", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_T", NC_DOUBLE, 2, dims, &vid);
  dims[0] = grid_y_C;
  dims[1] = grid_x_C;
  nc_def_var(cdfid, "x_C", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_C", NC_DOUBLE, 2, dims, &vid);

  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid, "TIME", NC_DOUBLE, 1, dims, &vid);

  df_mom_init_data(dumpdata, df, cdfid);
  data = (df_mom_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    if (data->vars[n].ndims == 2) {
      df_mom_get_dimids(cdfid, data->vars[n], &dims[2], &dims[1], NULL);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 3, dims, &vid);
    } else if (data->vars[n].ndims == 3) {
      df_mom_get_dimids(cdfid, data->vars[n], &dims[3], &dims[2], &dims[1]);
      nc_def_var(cdfid, df->vars[n], data->vars[n].type, 4, dims, &vid);
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }
  }

  write_dump_attributes_mom(dumpdata, cdfid, fptype, df->nce1,
			    df->nce2, df->nz, df->name);

  nc_enddef(cdfid);

  df_mom_writegeom(dumpdata, df);

  return data;
}



void df_mom_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = t;
  df_mom_data_t *data = (df_mom_data_t *)df->private_data;
  int fid = data->fid;
  momgrid_t *momgrid = dumpdata->momgrid;
  double **d2, ***d3;

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, dumpdata->output_tunit, 
		       &newt, 1);
  nc_put_vara_double(fid, ncw_var_id(fid, "TIME"), start, count, &newt);

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_mom_var_t vn = data->vars[n];
    void *p = vn.v;
    int gid = vn.xylocation;
    int ne1 = 0, ne2 = 0, nz = momgrid->nz;
    if (vn.ndims == 2) {
      start[1] = df->jlower;
      start[2] = df->ilower;
      if (gid & CL_CENTRE) {
	d2 = momgrid->T2;
	ne1 = momgrid->grid_x_T;
	ne2 = momgrid->grid_y_T;
      } else {
	d2 = momgrid->C2;
	ne1 = momgrid->grid_x_C;
	ne2 = momgrid->grid_y_C;
      }
      df_mom_get_dimsizes(df, vn, momgrid, &count[2], &count[1], NULL);
      pack_mom2d(gid, df, ne1, ne2, nz, (*((double ***)p)), d2, vn.scale);
      nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count, d2);
    } else if (vn.ndims == 3) {
      start[1] = df->klower;
      start[2] = df->jlower;
      start[3] = df->ilower;
      if (gid & CL_CENTRE)
	d3 = momgrid->T3;
      else
	d3 = momgrid->C3;
      df_mom_get_dimsizes(df, vn, momgrid, &count[3], &count[2], &count[1]);
      pack_mom3d(gid, df, ne1, ne2, nz, (*((double ****)p)), d3, vn.scale);
      nc_d_writesub_3d(fid, ncw_var_id(fid, df->vars[n]), start, count, d3);
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

void df_mom_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_mom_data_t *data = (df_mom_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}

void pack_mom2d(int gid, dump_file_t *df, int nx, int ny, int nz,
		double **in, double **out, double sc)
{
  int i, j;
  int im, jm;
  int is = df->ilower;
  int js = df->jlower;
  int ie = df->ilower + df->nce1;
  int je = df->jlower + df->nce2;

  if (gid & CL_CENTRE) {
    for (j = js; j < je; j++) {
      for (i = is; i < ie; i++) {
	out[j][i] = in[j][i] * sc;
      }
    }
  } else if (gid & CL_LEFT) {
    for (j = js; j < je; j++) {
      for (i = is; i < ie; i++) {
	jm = (j < ny - 1) ? j + 1 : j;
	out[j][i] = 0.5 * (in[j][i+1] + in[jm][i+1]) * sc;
      }
    }
  } else if (gid & CL_BACK) {
    for (j = js; j < je; j++) {
      for (i = is; i < ie; i++) {
	im = (i < nx - 1) ? i : i + 1;
	out[j][i] = 0.5 * (in[j+1][i] + in[j+1][im]) * sc;
      }
    }
  }
  /*
  for (j = js; j < je; j++) {
    for (i = is; i < ie; i++) {
      if (dumpdata->flag[nz-1][j][i] & (SOLID | OUTSIDE) || isnan(out[j][i]))
        out[j][i] = 0.0;
    }
  }
  */
}


void pack_mom3d(int gid, dump_file_t *df, int nx, int ny, int nz, 
		double ***in, double ***out, double sc)
{
  int i, j, k, kk;
  int im, jm;
  int is = df->ilower;
  int js = df->jlower;
  int ks = df->klower;
  int ie = df->ilower + df->nce1;
  int je = df->jlower + df->nce2;
  int ke = df->klower + df->nz;

  if (gid & CL_CENTRE) {
    for (k = ks; k < ke; k++) {
      kk = nz - 1 - k;
      for (j = js; j < je; j++) {
	for (i = is; i < ie; i++) {
	  out[kk][j][i] = in[k][j][i] * sc;
	}
      }
    }
  } else if (gid & CL_LEFT) {
    for (k = ks; k < ke; k++) {
      kk = nz - 1 - k;
      for (j = js; j < je; j++) {
	for (i = is; i < ie; i++) {
	  jm = (j < ny - 1) ? j : j + 1;
	  out[kk][j][i] = 0.5 * (in[k][j][i+1] + in[k][jm][i+1]) * sc;
	}
      }
    }
  } else if (gid & CL_BACK) {
    for (k = ks; k < ke; k++) {
      kk = nz - 1 - k;
      for (j = js; j < je; j++) {
	for (i = is; i < ie; i++) {
	  im = (i < nx - 1) ? i : i + 1;
	  out[kk][j][i] = 0.5 * (in[k][j+1][i] + in[k][j+1][im]) * sc;
	}
      }
    }
  }
  /*
  for (k = ks; k < ke; k++) {
    kk = nz - 1 - k;
    for (j = js; j < je; j++) {
      for (i = is; i < ie; i++) {
	if (dumpdata->flag[k][j][i] & (SOLID | OUTSIDE) || isnan(out[kk][j][i]))
	  out[kk][j][i] = 0.0;
      }
    }
  }
  */
}


static void df_mom_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_mom_data_t *data = (df_mom_data_t *)df->private_data;
  momgrid_t *momgrid = dumpdata->momgrid;
  int fid = data->fid;
  int i;
  float x[df->nfe1], y[df->nfe2];

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /* time independent variables */

  count[0] = df->nce1;
  for(i = df->ilower; i < df->ilower + df->nce1; i++)
    x[i - df->ilower] = momgrid->grid_x_Ta[i];
  nc_put_vara_float(fid, ncw_var_id(fid, "grid_x_T"), start, count, x);

  count[0] = df->nce2;
  for(i = df->jlower; i < df->jlower + df->nce2; i++)
    y[i - df->jlower] = momgrid->grid_y_Ta[i];
  nc_put_vara_float(fid, ncw_var_id(fid, "grid_y_T"), start, count, y);

  count[0] = df->nfe1;
  for(i = df->ilower; i < df->ilower + df->nfe1; i++)
    x[i - df->ilower] = momgrid->grid_x_Ca[i];
  nc_put_vara_float(fid, ncw_var_id(fid, "grid_x_C"), start, count, x);

  count[0] = df->nfe2;
  for(i = df->jlower; i < df->jlower + df->nfe2; i++)
    y[i - df->jlower] = momgrid->grid_y_Ca[i];
  nc_put_vara_float(fid, ncw_var_id(fid, "grid_y_C"), start, count, y);

  start[0] = df->klower;
  count[0] = df->nz;
  nc_put_vara_float(fid, ncw_var_id(fid, "zt"), start, count,
                     momgrid->zt);
  nc_put_vara_float(fid, ncw_var_id(fid, "zb"), start, count,
                     momgrid->zb);

  if (df->long0_360) {
    int j;
    for (j = 0; j < df->nce2; j++)
      for (i = 0; i < df->nce1; i++)
	momgrid->x_T[j][i] += 180.0;
    for (j = 0; j < df->nfe2; j++)
      for (i = 0; i < df->nfe1; i++)
	momgrid->x_C[j][i] += 180.0;
  }
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_T"), start, count,
                   momgrid->x_T);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_T"), start, count,
                   momgrid->y_T);
  count[0] = df->nfe2;
  count[1] = df->nfe1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_C"), start, count,
                   momgrid->x_C);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_C"), start, count,
                   momgrid->y_C);
               
  ncw_sync(df->name,fid);
  if (df->long0_360) {
    int j;
    for (j = 0; j < df->nce2; j++)
      for (i = 0; i < df->nce1; i++)
	momgrid->x_T[j][i] -= 180.0;
    for (j = 0; j < df->nfe2; j++)
      for (i = 0; i < df->nfe1; i++)
	momgrid->x_C[j][i] -= 180.0;
  }
}


static int df_mom_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
                              char *name, df_mom_var_t *var)
{
  int found = 1;
  int n = 0;
  double fact = 4.0e3*1025.0;  /* Conversion Wm-2 to ms-1K */
  double hlv = 2.500e6;        /* MOM latent heat of evaporation (J/kg) */

  var->v = NULL;
  var->ndims = 2;
  var->type = nc_get_default_type(df->bpv);
  var->xylocation = CL_CENTRE;
  var->zlocation = CL_NONE;
  var->scale = 1.0;

  if (strcmp(name, "u1av") == 0) {
    var->v = (void **)&dumpdata->u1av;
    var->xylocation = CL_GRID|CL_LEFT;
  }

  else if (strcmp(name, "u2av") == 0) {
    var->v = (void **)&dumpdata->u2av;
    var->xylocation = CL_GRID|CL_BACK;
  }

  else if (strcmp(name, "wtop") == 0) {
    var->v = (void **)&dumpdata->wtop;
  }

  else if (strcmp(name, "eta_t") == 0) {
    var->v = (void **)&dumpdata->eta;
  }

  else if (strcmp(name, "taux") == 0) {
    var->v = (void **)&dumpdata->wind1;
    var->xylocation = CL_GRID|CL_LEFT;
  }

  else if (strcmp(name, "tauy") == 0) {
    var->v = (void **)&dumpdata->wind2;
    var->xylocation = CL_GRID|CL_BACK;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void **)&dumpdata->patm;
  }

  else if (strcmp(name, "u1") == 0) {
    var->v = (void **)&dumpdata->u1;
    var->ndims = 3;
    var->xylocation = CL_GRID|CL_LEFT;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "u2") == 0) {
    var->v = (void **)&dumpdata->u2;
    var->ndims = 3;
    var->xylocation = CL_GRID|CL_BACK;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void **)&dumpdata->w;
    var->ndims = 3;
    var->xylocation = CL_CENTRE;
    var->zlocation = CL_GRID;
  }

  else if (strcmp(name, "dens") == 0) {
    var->v = (void **)&dumpdata->dens;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "dens_0") == 0) {
    var->v = (void **)&dumpdata->dens_0;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "Kz") == 0) {
    var->v = (void **)&dumpdata->Kz;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "Vz") == 0) {
    var->v = (void **)&dumpdata->Vz;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "Cd") == 0) {
    var->v = (void **)&dumpdata->Cd;
  }

  else if (strcmp(name, "u1vh") == 0) {
    var->v = (void **)&dumpdata->u1vh;
  }

  else if (strcmp(name, "u2vh") == 0) {
    var->v = (void **)&dumpdata->u2vh;
  }
  else
    found = 0;

  if (!found) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_wc[n];
        var->xylocation = CL_CENTRE;
	var->zlocation = CL_CENTRE;
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
	/* Apply any MOM scaling */
	if (strcmp(name, "lhf") == 0) var->scale = -1.0 / hlv;
	if (strcmp(name, "shf") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "swr") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "lwr") == 0) var->scale = -1.0 * fact;
	if (strcmp(name, "precip") == 0) var->scale = 1025.0;
	if (strcmp(name, "lprec") == 0) var->scale = 1025.0;
	if (strcmp(name, "evap") == 0) var->scale = 1025.0;
        found = 1;
        break;
      }
    }
  }

  return found;
}


static void df_mom_init_data(dump_data_t *dumpdata, dump_file_t *df, int fid)
{
  int i, n;
  df_mom_data_t *data = NULL;
  char line[MAXNUMVARS * 20] = "";

  for (i = 0; i < df->nvars; ++i) {
    if (contains_token(df->vars[i], "ALL")) {
      sprintf(line, MOM_ALL_VARS);
      for (n = 0; n < dumpdata->ntr; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);

      for (n = 0; n < dumpdata->ntrS; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_2d[n].name);
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    } else if (contains_token(df->vars[i], "NONTRACERS"))
      sprintf(line, MOM_ALL_VARS);
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

  data = (df_mom_data_t *)malloc(sizeof(df_mom_data_t));
  memset(data, 0, sizeof(df_mom_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_mom_var_t *)malloc(df->nvars * sizeof(df_mom_var_t))) == NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_mom_var_t *var = &data->vars[i];
    if (df_mom_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_mom_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

}


static void df_mom_get_dimids(int cdfid, df_mom_var_t var, 
			      int *e1id, int *e2id, int *zid)
{
  if (var.xylocation & CL_CENTRE) {
    *e1id = ncw_dim_id(cdfid, "grid_x_T");
    *e2id = ncw_dim_id(cdfid, "grid_y_T");
  } else if (var.xylocation & CL_GRID) {
    *e1id = ncw_dim_id(cdfid, "grid_x_C");
    *e2id = ncw_dim_id(cdfid, "grid_y_C");
  }

  switch (var.zlocation) {
  case CL_GRID:
    *zid = ncw_dim_id(cdfid, "zb");
    break;
  case CL_CENTRE:
    *zid = ncw_dim_id(cdfid, "zt");
    break;
    
  default:
    break;
  }
}

static void df_mom_get_dimsizes(dump_file_t *df, df_mom_var_t var,
				momgrid_t *momgrid,
				size_t *ne1, size_t *ne2, size_t *nz)

{
  if (var.xylocation & CL_CENTRE) {
    *ne1 = df->nce1;
    *ne2 = df->nce2;
  } else if (var.xylocation & CL_GRID) {
    *ne1 = df->nfe1;
    *ne2 = df->nfe2;
  }
  switch (var.zlocation) {
  case CL_GRID:
    *nz = momgrid->nz;
    break;

  case CL_CENTRE:
    *nz = momgrid->nz;
    break;

  default:
    break;
  }
}


/** Writes tracer attributes to a NetCDF file.
 * `attr' and `attrval' are common attributes added to all tracers.
 * @param fid NetCDF file
 * @param ntr number of tracers
 * @param tracers array of `tracer_desc'
 * @param attr array of additional attributes to be added to all tracers
 */
void tracer_write_mom_nc(int fid, int ntr, tracer_info_t tracers[], int nattr,
			 tracer_att_t attr[])
{
  int n;
  char name[MAXSTRLEN];
  char lname[MAXSTRLEN];
  char units[MAXSTRLEN];

  for (n = 0; n < ntr; n++) {
    int vid;
    tracer_info_t *tr = &tracers[n];
    strcpy(name,tr->name);
    if (strcmp(attr[0].attr,"tracer_sed") == 0)
      strcat(name, "_sed");
    strcpy(lname, tr->long_name);
    strcpy(units, tr->units);
    if (strcmp(name, "lhf") == 0) {
      strcpy(lname, "Evaporation rate");
      strcpy(units, "kg/m2/s");
    }
    if (strcmp(name, "evap") == 0) strcpy(units, "(kg/m3)(m/s)");
    if (strcmp(name, "precip") == 0) strcpy(units, "(kg/m3)(m/s)");
    if (strcmp(name, "lprec") == 0) strcpy(units, "(kg/m3)(m/s)");
    vid = ncw_var_id(fid, name);
    if (vid >= 0) {
      nc_put_att_text(fid, vid, "long_name", strlen(lname), lname);
      nc_put_att_text(fid, vid, "units", strlen(units), units);
      nc_put_att_double(fid, vid, "_FillValueWC", NC_DOUBLE, 1,
                        &tr->fill_value_wc);
      nc_put_att_double(fid, vid, "valid_range", NC_DOUBLE, 2,
                        tr->valid_range_wc);
    }
  }
}


