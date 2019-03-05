/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/inputs/readdump.c
 *  
 *  Description:
 *  Routine to read model netCDF dump file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: readdump.c 6147 2019-03-05 01:58:38Z her127 $
 *
 */

#include <string.h>
#include <netcdf.h>
#include "hd.h"
#include "tracer.h"


void read1d(geometry_t *geom, int id, char *name, double *p,
            int ni, int nj, int nk);
int read2d(geometry_t *geom, int id, char *name, double *p, int dump,
	   int nj, int ni);
int dumpdata_read_2d(dump_data_t *dumpdata, int id, char *name,
		     double **p, int dump, int nj, int ni);
int dumpdata_read_3d(dump_data_t *dumpdata, int id, char *name,
		     double ***p, int dump, int nk, int nj, int ni);
int dumpdata_read_3d_tr(dump_data_t *dumpdata, int id, char *name,
			double ***p, int dump, int nk, int nj, int ni);
void zero2d(dump_data_t *dumpdata, double **dst, int n1, int n2);
void zero3d(dump_data_t *dumpdata, double ***dst, int n1, int n2, int n3);
double smoothbathy(double **a, int c, int p, int m, int pp, int mm, int cc,
                   int mode, double bmin);
void intp_undef_s(geometry_t *geom, double *a);
void mean_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy);
void set_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy, double val);
void adjust_outside_obc(parameters_t * params, double **bathy);
int iswet(int bathy);
void read_sed_layers(geometry_t *geom, parameters_t *params, master_t *master, dump_data_t *dumpdata);
void remove_channel(parameters_t *params, unsigned long ***flag, double **bathy);
void data_infill_3d(master_t *master, double ***a);
void data_infill_2d(master_t *master, double **a);
double smooth_bathy(unsigned long ***flag, int nz, int i, int j, int **xm, int **xp, 
		    int **ym, int **yp, double **bathy, double bmin, double bmax);

/*-----------------------------------------------------------------------*/
/* Reads the input data from netCDF file                                 */
/*-----------------------------------------------------------------------*/
int dump_read(geometry_t *geom, /* Sparse global geometry structure */
              parameters_t *params, /* Input parameter data structure */
              master_t *master, /* Master data structure */
              int cdfid, int ti)
{
  int c, cc, i, j, k, n;
  size_t ndumps;
  size_t start[4];
  size_t count[4];

  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    hd_warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  if (DEBUG("dump"))
    dlog("dump", "Reading dump %d.\n", ti);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = geom->nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*---------------------------------------------------------------------*/
  /* Initialise : this sets ghost cells to zero which are set later      */
  memset(geom->h1au1, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2au2, 0, geom->sgsizS * sizeof(double));
  memset(geom->h1acell, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2acell, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2au1, 0, geom->sgsizS * sizeof(double));
  memset(geom->h1au2, 0, geom->sgsizS * sizeof(double));
  memset(geom->botz, 0, geom->sgsizS * sizeof(double));

  /*---------------------------------------------------------------------*/
  /* Read the 1 dimensional input data into the geometric sparse array   */
  read1d(geom, cdfid, "z_grid", geom->gridz, geom->nce1, geom->nce2,
         geom->nz + 1);
  read1d(geom, cdfid, "z_centre", geom->cellz, geom->nce1, geom->nce2,
         geom->nz);

  /*---------------------------------------------------------------------*/
  /* Read the 2 dimensional input data into the geometric sparse array   */
  read2d(geom, cdfid, "x_grid", geom->gridx, -1, geom->nfe2, geom->nfe1);
  read2d(geom, cdfid, "y_grid", geom->gridy, -1, geom->nfe2, geom->nfe1);
  /* 
  read2d(geom,cdfid,"h1agrid",geom->h1agrid,0,geom->nfe2,geom->nfe1);
  read2d(geom,cdfid,"h2agrid",geom->h2agrid,0,geom->nfe2,geom->nfe1);
  */
  read2d(geom, cdfid, "x_centre", geom->cellx, -1, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "y_centre", geom->celly, -1, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "h1acell", geom->h1acell, -1, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "h2acell", geom->h2acell, -1, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "botz", geom->botz, -1, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "x_left", geom->u1x, -1, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "y_left", geom->u1y, -1, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "h1au1", geom->h1au1, -1, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "h2au1", geom->h2au1, -1, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "thetau1", geom->thetau1, -1, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "x_back", geom->u2x, -1, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "y_back", geom->u2y, -1, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "h1au2", geom->h1au2, -1, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "h2au2", geom->h2au2, -1, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "thetau2", geom->thetau2, -1, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "coriolis", master->coriolis, -1, geom->nce2, geom->nce1);

  /*---------------------------------------------------------------------*/
  /* Time dependent variables                                            */
  start[0] = ti;
  count[0] = 1L;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "t"), start, count,
                     &params->t);

  if (params->timeunit) {
    char timeunits[MAXSTRLEN];
    memset(timeunits, 0, MAXSTRLEN);
    nc_get_att_text(cdfid, ncw_var_id(cdfid, "t"), "units", timeunits);
    tm_change_time_units(timeunits, params->timeunit, &params->t, 1);
  }

  /*
  nc_get_vara_long(cdfid, ncw_var_id(cdfid, "flag"), start, count,
		   (long *)dumpdata->flag[0][0]);
  if (params->sigma) {
    int nz = dumpdata->nz - 1;
    for (i = 0; i < dumpdata->nfe1; i++)
      for (j = 0; j < dumpdata->nfe2; j++)
	if (!(dumpdata->flag[nz][j][i] & (SOLID | OUTSIDE))) {
	  for (k = 0; k < nz; k++)
	    dumpdata->flag[k][j][i] = dumpdata->flag[nz][j][i];
	}
  }
  */

  /* Read the 2 dimensional input data into the master sparse array */
  read2d(geom, cdfid, "u1av", master->u1av, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "u1bot", master->u1bot, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "wind1", master->wind1, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "u2av", master->u2av, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "u2bot", master->u2bot, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "wind2", master->wind2, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "wtop", master->wtop, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "topz", master->topz, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "eta", master->eta, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "patm", master->patm, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "Cd", master->Cd, ti, geom->nce2, geom->nce1);

  /* Read the 3 dimensional input data into the master sparse array */
  read3d(geom, cdfid, "u1", master->u1, ti, geom->nz, geom->nce2,
         geom->nfe1);
  read3d(geom, cdfid, "u2", master->u2, ti, geom->nz, geom->nfe2,
         geom->nce1);
  read3d(geom, cdfid, "w", master->w, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "dens", master->dens, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "dens_0", master->dens_0, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "Kz", master->Kz, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "Vz", master->Vz, ti, geom->nz, geom->nce2,
         geom->nce1);

  /* Check for NaNs                                                  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    int ci, cj;
    c = geom->w2_t[cc];
    if (isnan(master->eta[c])) {
      hd_warn("dump_read: Found NaN at (%d %d).\n", geom->s2i[c], geom->s2j[c]);
    }
  }

  /* Read the 3D tracer data into the master sparse array */
  for (n = 0; n < master->ntr; n++) {
    /* If a 3D tracer is included in the input file then initialise  */
    /* with this data. If not use the specification in the parameter */
    /* file tracer list; TRACER#.data                                */
    if (read3d(geom, cdfid, master->trinfo_3d[n].name, master->tr_wc[n], ti,
	       geom->nz, geom->nce2, geom->nce1))
      load_wc_tracer_name(master, params->prmfd, master->trinfo_3d[n].name, WATER);

    /* Ignore DA variables */
    if (strncasecmp(master->trinfo_3d[n].tag, "DA_OBS", 6) == 0 ) continue;

    /* Check for NaN's in the initial condition */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      if (isnan(master->tr_wc[n][c]))
	hd_warn
	  ("dump_read: NaN found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
    }
  }
  /* Read the 2D tracer data into the master sparse array */
  for (n = 0; n < master->ntrS; n++) {
    if (read2d(geom, cdfid, master->trinfo_2d[n].name, master->tr_wcS[n], ti,
	       geom->nce2, geom->nce1))
      load_wc_tracer_name(master, params->prmfd, master->trinfo_2d[n].name, INTER);
  }

  /* Sediment thickness */
  if (params->sednz > 0) {
    double ***d1;
    d1 = d_alloc_3d(geom->nce1, geom->nce2, geom->sednz + 1);
    dumpdata_read_3d(dumpdata, cdfid, "z_grid_sed", d1, ti,
                     dumpdata->sednz + 1, dumpdata->nce2, dumpdata->nce1);
    for (k = 0; k < params->sednz + 1; k++)
      c2s_2d(geom, geom->gridz_sed[k], d1[k], geom->nce1, geom->nce2);
    d_free_3d(d1);

    d1 = d_alloc_3d(geom->nce1, geom->nce2, geom->sednz);
    dumpdata_read_3d(dumpdata, cdfid, "z_centre_sed", d1, ti,
                     dumpdata->sednz, dumpdata->nce2, dumpdata->nce1);
    for (k = 0; k < params->sednz; k++)
      c2s_2d(geom, geom->cellz_sed[k], d1[k], geom->nce1, geom->nce2);

    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      if (dumpdata_read_3d(dumpdata, cdfid, name, d1, ti,
			   dumpdata->sednz, dumpdata->nce2, dumpdata->nce1)) {
	load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_sed[n].name, SEDIM);
      } else {
	for (k = 0; k < params->sednz; k++)
	  c2s_2d(geom, master->tr_sed[n][k], d1[k], geom->nce1, geom->nce2);
      }
    }
    d_free_3d(d1);
  }
  read_mean_atts(master, cdfid);
  return (1);
}

/* END dump_read()                                                   */
/*-------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/* Reads the time dependent input data from netCDF file                  */
/*-----------------------------------------------------------------------*/
int dump_re_read(master_t *master, /* Master data structure */
		 int cdfid, int ti)
{
  geometry_t *geom = master->geom;
  int c, cc, i, j, k, n;
  size_t ndumps;
  size_t start[4];
  size_t count[4];

  int time_vid;
  char ftimeunits[MAXSTRLEN];
  memset(ftimeunits,0,MAXSTRLEN);
  
  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    hd_warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  if (DEBUG("dump"))
    dlog("dump", "Reading dump %d.\n", ti);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* Time independent variables */
  count[0] = geom->nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  /*---------------------------------------------------------------------*/
  /* Time dependent variables                                            */
  start[0] = ti;
  count[0] = 1L;
  time_vid = ncw_var_id(cdfid, "t");
  nc_get_vara_double(cdfid, time_vid, start, count, &master->t);

  /* Make sure time value is consistent with the time units */
  nc_get_att_text(cdfid, time_vid, "units", ftimeunits);
  if (strcmp(ftimeunits, master->timeunit) != 0)
    tm_change_time_units(ftimeunits, master->timeunit, &master->t, 1);
  
  /* Read the 2 dimensional input data into the master sparse array */
  read2d(geom, cdfid, "u1av", master->u1av, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "u1bot", master->u1bot, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "wind1", master->wind1, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "u2av", master->u2av, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "u2bot", master->u2bot, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "wind2", master->wind2, ti, geom->nfe2, geom->nce1);
  read2d(geom, cdfid, "wtop", master->wtop, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "topz", master->topz, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "eta", master->eta, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "patm", master->patm, ti, geom->nce2, geom->nce1);
  read2d(geom, cdfid, "Cd", master->Cd, ti, geom->nce2, geom->nce1);
  /* Note: use the dumpdata->u1av array to read the backward vel */
  read2d(geom, cdfid, "u1av", master->u1avb, ti, geom->nce2, geom->nfe1);
  read2d(geom, cdfid, "u2av", master->u2avb, ti, geom->nfe2, geom->nce1);

  /* Read the 3 dimensional input data into the master sparse array */
  read3d(geom, cdfid, "u1", master->u1, ti, geom->nz, geom->nce2,
         geom->nfe1);
  read3d(geom, cdfid, "u2", master->u2, ti, geom->nz, geom->nfe2,
         geom->nce1);
  read3d(geom, cdfid, "w", master->w, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "dens", master->dens, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "dens_0", master->dens_0, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "Kz", master->Kz, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "Vz", master->Vz, ti, geom->nz, geom->nce2,
         geom->nce1);
  read3d(geom, cdfid, "u1", master->u1b, ti, geom->nz, geom->nce2,
         geom->nfe1);
  read3d(geom, cdfid, "u2", master->u2b, ti, geom->nz, geom->nfe2,
         geom->nce1);

  /* Check for NaNs                                                  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    int ci, cj;
    c = geom->w2_t[cc];
    if (isnan(master->eta[c])) {
      hd_warn("dump_read: Found NaN at (%d %d).\n", geom->s2i[c], geom->s2j[c]);
    }
  }

  /* Read the 3D tracer data into the master sparse array */
  for (n = 0; n < master->ntr; n++) {
    /* Ignore DA variables */
    if (strncmp(master->trinfo_3d[n].tag, "DA_", 3) == 0 ) continue;

    read3d(geom, cdfid, master->trinfo_3d[n].name, master->tr_wc[n], ti,
	   geom->nz, geom->nce2, geom->nce1);
    /* Check for NaN's in the initial condition */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      if (isnan(master->tr_wc[n][c]))
	hd_warn
	  ("dump_re_read: NaN found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
    }
  }
  /* Read the 2D tracer data into the master sparse array */
  for (n = 0; n < master->ntrS; n++) {
    read2d(geom, cdfid, master->trinfo_2d[n].name, master->tr_wcS[n], ti,
           geom->nce2, geom->nce1);
  }

  /* Sediment thickness */
  if (geom->sednz > 0) {
    double ***d1 = d_alloc_3d(dumpdata->nce1, dumpdata->nce2, dumpdata->sednz);
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      dumpdata_read_3d(dumpdata, cdfid, name, d1, ti,
                       dumpdata->sednz, dumpdata->nce2, dumpdata->nce1);
      for (k = 0; k < geom->sednz; k++)
        c2s_2d(geom, master->tr_sed[n][k], d1[k], geom->nce1, geom->nce2);
    }
    d_free_3d(d1);
  }

  return (1);
}

/* END dump_re_read()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 2d cartesian array from id, zeroing it if an error occurs  */
/* and copying it into the sparse array.                             */
/*-------------------------------------------------------------------*/
void read1d(geometry_t *geom, int id, char *name, double *p,
            int ni, int nj, int nk)
{
  double *array = d_alloc_1d(nk); /* Dummy cartesian array */
  int c;                        /* Sparse couters */
  int i, j, k;                  /* Cartesian couters */

  size_t start[4];
  size_t count[4];

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = nk;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  if (nc_get_vara_double(id, ncw_var_id(id, name), start, count, array) <
      0) {
    memset(p, 0, geom->sgsiz * sizeof(double));
    return;
  }

  /* Copy into the sparse array */
  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < ni && j < nj)
      p[c] = array[k];
  }
  d_free_1d(array);
}

/* END read1d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 2d cartesian array from id, zeroing it if an error occurs  */
/* and copying it into the sparse array.                             */
/*-------------------------------------------------------------------*/
int read2d(geometry_t *geom, int id, char *name, double *p, int dump,
	   int nj, int ni)
{
  double **array = d_alloc_2d(ni, nj);  /* Dummy cartesian array */
  int c;                        /* Sparse couters */
  int i, j, k;                  /* Cartesian couters */

  size_t start[4];
  size_t count[4];

  if (dump >= 0) {
    start[0] = dump;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = 1L;
    count[1] = nj;
    count[2] = ni;
    count[3] = 0;
  } else {
    start[0] = 0L;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = nj;
    count[1] = ni;
    count[2] = 0;
    count[3] = 0;
  }
  if (nc_get_vara_double(id, ncw_var_id(id, name), start, count, array[0])
      < 0) {
    memset(p, 0, geom->sgsizS * sizeof(double));
    return(1);
  }

  /* Copy into the sparse array */
  for (c = 1; c <= geom->sgnumS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < ni && j < nj)
      p[c] = array[j][i];
  }
  d_free_2d(array);
  return(0);
}

/* END read2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 3d cartesian array from id, zeroing it if an error occurs  */
/* and copying it into the sparse array.                             */
/*-------------------------------------------------------------------*/
int read3d(geometry_t *geom, int id, char *name, double *p, int dump,
            int nk, int nj, int ni)
{
  double ***array = d_alloc_3d(ni, nj, nk); /* Dummy cartesian array */
  int c;                        /* Sparse couters */
  int i, j, k;                  /* Cartesian couters */

  size_t start[4];
  size_t count[4];

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nk;
  count[2] = nj;
  count[3] = ni;
  if (nc_get_vara_double
      (id, ncw_var_id(id, name), start, count, array[0][0]) < 0) {
    memset(p, 0, geom->sgsiz * sizeof(double));
    return(1);
  }

  /* Copy into the sparse array */
  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < ni && j < nj)
      p[c] = array[k][j][i];
  }
  d_free_3d(array);
  return(0);
}

/* END read3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set variables at ghost cells                           */
/*-------------------------------------------------------------------*/
void master_setghosts(geometry_t *geom, /* Sparse global geometry    */
                      master_t *master, /* Master data               */
		      double **cellx,   /* Cell centre x location    */
		      double **celly,   /* Cell centre y location    */
		      double **gridx,   /* Cell corner x location    */
		      double **gridy,   /* Cell corner y location    */
		      int *s2i,
		      int *s2j
  )
{
  int cb, ci, c, cc;            /* Sparse couters */
  int n;
  double mx = 0;
  double *eta;

  /* Set a no-gradient condition over lateral boundaries */
  for (cc = 1; cc <= geom->nbptS; cc++) {
    cb = geom->bpt[cc];
    ci = geom->bin[cc];
    geom->h1au1[cb] = geom->h1au1[ci];
    geom->h2au2[cb] = geom->h2au2[ci];
    geom->h1acell[cb] = geom->h1acell[ci];
    geom->h2acell[cb] = geom->h2acell[ci];
    geom->h2au1[cb] = geom->h2au1[ci];
    geom->h1au2[cb] = geom->h1au2[ci];
    geom->botz[cb] = geom->botz[ci];
    if (master != NULL)
      master->coriolis[cb] = master->coriolis[ci];
  }

  /* Set the bottom depth on e1 and e2 open boundaries (no gradient) */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no2_e1; cc++) {
        c = open->obc_e1[cc];
        ci = open->oi1_e1[cc];
        geom->botz[c] = geom->botz[ci];
	geom->h1au1[c] = geom->h1au1[ci];
	geom->h1au2[c] = geom->h1au2[ci];
	geom->h2au2[c] = geom->h2au2[ci];
	geom->h2au1[c] = geom->h2au1[ci];
	geom->thetau1[c] = geom->thetau1[ci];
	geom->thetau2[c] = geom->thetau2[ci];
	geom->h1acell[c] = geom->h1acell[ci];
	geom->h2acell[c] = geom->h2acell[ci];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no2_e2; cc++) {
        c = open->obc_e2[cc];
        ci = open->oi1_e2[cc];
        geom->botz[c] = geom->botz[ci];
	geom->h1au1[c] = geom->h1au1[ci];
	geom->h1au2[c] = geom->h1au2[ci];
	geom->h2au2[c] = geom->h2au2[ci];
	geom->h2au1[c] = geom->h2au1[ci];
	geom->thetau1[c] = geom->thetau1[ci];
	geom->thetau2[c] = geom->thetau2[ci];
	geom->h1acell[c] = geom->h1acell[ci];
	geom->h2acell[c] = geom->h2acell[ci];
      }
    }
    /* Boundary length */
    open->length = 0.0;
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      if (open->type & U1BDRY)
	open->length += geom->h2au1[c];
      else
	open->length += geom->h1au2[c];
    }

    /* Boundary area */
    open->area = 0.0;
    for (cc = 1; cc <= open->no3_t; cc++) {
      int cs, zp1;
      double dz;
      c = open->obc_t[cc];
      cs = geom->m2d[c];
      zp1 = geom->zp1[c];
      dz = min(master->eta[cs], geom->gridz[zp1]) -
	max(geom->botz[cs], geom->gridz[c]);
      if (open->type & U1BDRY)
	open->area += dz * geom->h2au1[cs];
      else
	open->area += dz * geom->h1au2[cs];
    }

    /* Boundary depths */
    open->meandep = open->maxdep = 0.0;
    open->mindep = -1e10;
    for (cc = 1; cc <= open->no2_t; cc++) {
      int cs;
      c = open->obc_t[cc];
      open->maxdep = min(open->maxdep, geom->botz[c]);
      open->mindep = max(open->mindep, geom->botz[c]);
      open->meandep += geom->botz[c];
    }
    open->meandep /= (double)open->no2_t;

    /* Cell location and bottom depth */
    for (cc = 1; cc <= open->no2_t; cc++) {
      int i, j, bi, bj, m, m1, m2;
      c = ci = open->obc_t[cc];
      i = bi = geom->s2i[c];
      j = bj = geom->s2j[c];
      for (m = 0; m < open->bgz; m++) {
	m1 = c;
	m2 = open->nmap[c];
	c = open->omap[c];
	if (open->type & U1BDRY && open->ocodex & R_EDGE) {
	  i = (i == geom->nce1) ? i : i+1;
	  if(i >= bi || (isnan(cellx[j][i]))) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m2], geom->cellx[m1], 0, 1, 2);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m2], geom->celly[m1], 0, 1, 2);
	  } else {
	    geom->cellx[c] = cellx[j][i];
	    geom->celly[c] = celly[j][i];
	  }
	}
	if (open->type & U1BDRY && open->ocodex & L_EDGE) {
	  i = (i == 0) ? i : i-1;
	  if(i <= bi || (isnan(cellx[j][i]))) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m1], geom->cellx[m2], 1, 2, 0);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m1], geom->celly[m2], 1, 2, 0);
	    else {
	      geom->cellx[c] = cellx[j][i];
	      geom->celly[c] = celly[j][i];
	    }
	  }
	}
	if (open->type & U2BDRY && open->ocodey & F_EDGE) {
	  j = (j == geom->nce2) ? j : j+1;
	  if(j >= bj || (isnan(celly[j][i]))) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m2], geom->cellx[m1], 0, 1, 2);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m2], geom->celly[m1], 0, 1, 2);
	  } else {
	    geom->cellx[c] = cellx[j][i];
	    geom->celly[c] = celly[j][i];
	  }
	}
	if (open->type & U2BDRY && open->ocodey & B_EDGE) {
	  j = (j == 0) ? j : j-1;
	  if(j <= bj || (isnan(celly[j][i]))) {
	    if(!isnan(geom->cellx[m1]) && !isnan(geom->cellx[m2]))
	      geom->cellx[c] = intp(geom->cellx[m1], geom->cellx[m2], 1, 2, 0);
	    if(!isnan(geom->celly[m1]) && !isnan(geom->celly[m2]))
	      geom->celly[c] = intp(geom->celly[m1], geom->celly[m2], 1, 2, 0);
	    else {
	      geom->cellx[c] = cellx[j][i];
	      geom->celly[c] = celly[j][i];
	    }
	  }
	}
	/* Uniform depth and metrics over the OBC ghosts */
	geom->botz[c] = geom->botz[ci];
	geom->h1au1[c] = geom->h1au1[ci];
	geom->h1au2[c] = geom->h1au2[ci];
	geom->h2au2[c] = geom->h2au2[ci];
	geom->h2au1[c] = geom->h2au1[ci];
	geom->thetau1[c] = geom->thetau1[ci];
	geom->thetau2[c] = geom->thetau2[ci];
	geom->h1acell[c] = geom->h1acell[ci];
	geom->h2acell[c] = geom->h2acell[ci];
      }
    }
  }
  /* Median filter eta if required */
  eta = d_alloc_1d(geom->sgsizS);
  memcpy(eta, master->eta, geom->sgsizS * sizeof(double));
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      eta[open->omap[c]] = eta[open->nmap[c]];
    }
  }
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->options & OP_ETAFIL) {
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	master->eta[c] = con_median(geom, eta, 10.0, c);
      }
    }
  }
  d_free_1d(eta);

  /* Set the lateral boundaries again to pick up cells defined on */
  /* R_EDGE or F_EDGE OBCs. */
  for (cc = 1; cc <= geom->nbptS; cc++) {
    cb = geom->bpt[cc];
    ci = geom->bin[cc];
    geom->h1au1[cb] = geom->h1au1[ci];
    geom->h2au2[cb] = geom->h2au2[ci];
    geom->h1acell[cb] = geom->h1acell[ci];
    geom->h2acell[cb] = geom->h2acell[ci];
    geom->h2au1[cb] = geom->h2au1[ci];
    geom->h1au2[cb] = geom->h1au2[ci];
    geom->botz[cb] = geom->botz[ci];
    if (master != NULL)
      master->coriolis[cb] = master->coriolis[ci];
  }
  c = geom->map[geom->nz-1][geom->nce2][geom->nce1];
  ci = geom->xm1[geom->ym1[c]];
  if (geom->h1au1[c] == 0.0) geom->h1au1[c] = geom->h1au1[ci];
  if (geom->h2au2[c] == 0.0) geom->h2au2[c] = geom->h2au2[ci];

  /* Cell location */
  for (cc = 1; cc <= geom->nbptS; cc++) {
    int i, j, ig, jg;
    c = geom->bpt[cc];
    ci = geom->bin[cc];
    i = ig = geom->s2i[ci];
    j = jg = geom->s2j[ci];

    if (geom->xp1[ci] == c) ig = i + 1;
    else if (geom->xm1[ci] == c) ig = i - 1;
    else if (geom->yp1[ci] == c) jg = j + 1;
    else if (geom->ym1[ci] == c) jg = j - 1;
    else if (geom->yp1[geom->xp1[ci]] == c || geom->xp1[geom->yp1[ci]] == c) {
      ig = i + 1;
      jg = j + 1;
    } else if (geom->yp1[geom->xm1[ci]] == c || geom->xm1[geom->yp1[ci]] == c) {
      ig = i - 1;
      jg = j + 1;
    } else if (geom->ym1[geom->xp1[ci]] == c || geom->xp1[geom->ym1[ci]] == c) {
      ig = i + 1;
      jg = j - 1;
    } else if (geom->ym1[geom->xm1[ci]] == c || geom->xm1[geom->ym1[ci]] == c) {
      ig = i - 1;
      jg = j - 1;
    }

    ig = min(max(ig, 0), geom->nce1-1);
    jg = min(max(jg, 0), geom->nce2-1);

    geom->cellx[c] = cellx[jg][ig];
    geom->celly[c] = celly[jg][ig];
    geom->gridx[c] = gridx[jg][ig];
    geom->gridy[c] = gridy[jg][ig];
    if (s2i != NULL) {
      s2i[c] = ig;
      cellx[jg][ig] = cellx[j][i];
      gridx[jg][ig] = gridx[j][i];
    }
    if (s2j != NULL) {
      s2j[c] = jg;
      celly[jg][ig] = celly[j][i];
      gridy[jg][ig] = gridy[j][i];
    }
  }

  intp_undef_s(geom, geom->cellx);
  intp_undef_s(geom, geom->celly);
  intp_undef_s(geom, geom->gridx);
  intp_undef_s(geom, geom->gridy);

  /*
  intp_undef_s(geom, geom->thetau1);
  intp_undef_s(geom, geom->thetau2);
  */

  /* Get the surface coordinate corresponding to the maximum bottom */
  /* depth.  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->botz[c] < mx) {
      geom->cmx = cc;
      mx = geom->botz[c];
    }
  }
}

/* END master_setghosts()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate grid location at undefined cells           */
/*-------------------------------------------------------------------*/
void intp_undef_s(geometry_t *geom, double *a)
{
  int c;
  double *aa = d_alloc_1d(geom->sgnumS);
  memcpy(aa, a, geom->sgnumS * sizeof(double));

  for (c = 1; c <= geom->sgnumS; c++) {
    int xm1 = geom->xm1[c];
    int xp1 = geom->xp1[c];
    int ym1 = geom->ym1[c];
    int yp1 = geom->yp1[c];
    int xm2 = geom->xm1[xm1];
    int xp2 = geom->xp1[xp1];
    int ym2 = geom->ym1[ym1];
    int yp2 = geom->yp1[yp1];
    if (isnan(a[c])) {
      if (!isnan(aa[xm1]) && !isnan(aa[xm2]))
	a[c] = intp(a[xm2], a[xm1], 0, 1, 2);
      else if (!isnan(aa[xp1]) && !isnan(aa[xp2]))
	a[c] = intp(aa[xp1], aa[xp2], 1, 2, 0);
      else if (!isnan(aa[ym1]) && !isnan(aa[ym2]))
	a[c] = intp(aa[ym2], aa[ym1], 0, 1, 2);
      else if (!isnan(aa[yp1]) && !isnan(aa[yp2]))
	a[c] = intp(aa[yp1], aa[yp2], 1, 2, 0);
    }
  }
  d_free_1d(aa);
}

/* END intp_undef_s()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate grid location at undefined cells           */
/*-------------------------------------------------------------------*/
void intp_undef(int nce1, int nce2, double **a)
{
  int i, j;

  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (isnan(a[j][i])) {
	if (i > 1 && !isnan(a[j][i-1]) && !isnan(a[j][i-2]))
	  a[j][i] = intp(a[j][i-2], a[j][i-1], i-2, i-1, i);
	else if (i < nce1-2 && !isnan(a[j][i+1]) && !isnan(a[j][i+2]))
	  a[j][i] = intp(a[j][i+1], a[j][i+2], i+1, i+2, i);
	else if (j > 1 && !isnan(a[j-1][i]) && !isnan(a[j-2][i]))
	  a[j][i] = intp(a[j-2][i], a[j-1][i], j-2, j-1, j);
	else if (j < nce2-2 && !isnan(a[j+1][i]) && !isnan(a[j+2][i]))
	  a[j][i] = intp(a[j+1][i], a[j+2][i], j+1, j+2, j);
      }
    }
  }
}

/* END intp_undef()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set R_EDGE and F_EDGE interior boundaries to OUTSIDE   */
/*-------------------------------------------------------------------*/
void set_bdry_flags(geometry_t *geom, dump_data_t *dumpdata)
{
  int n;
  int c, cc, xp1;
  int i, j, k;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
        c = open->obc_e1[cc];
        i = geom->s2i[c];
        j = geom->s2j[c];
        k = geom->s2k[c];
	xp1 = geom->xp1[c];
	if (open->stagger & INFACE) i++;
        /*if (i < geom->nce1 && xp1 == geom->xp1[xp1]) {*/
        if (i < geom->nce1) {
          dumpdata->flag[k][j][i] |= OUTSIDE;
          dumpdata->botz[j][i] = NaN;
        }
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
        c = open->obc_e2[cc];
        i = geom->s2i[c];
        j = geom->s2j[c];
        k = geom->s2k[c];
	if (open->stagger & INFACE) j++;
        /*if (j < geom->nce2 && geom->yp1[geom->yp1[c]] == c) {*/
        if (j < geom->nce2) {
          dumpdata->flag[k][j][i] |= OUTSIDE;
          dumpdata->botz[j][i] = NaN;
        }
      }
    }
  }
}

/* END set_bdry_flags()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the input data into the dump data structure and   */
/* copy to the sparse array.                                         */
/*-------------------------------------------------------------------*/
int dumpdata_read(geometry_t *geom, /* Sparse global geometry structure */
                  parameters_t *params, /* Input parameter data structure */
                  master_t *master, /* Master data structure */
                  dump_data_t *dumpdata,  /* Dump data structure */
                  int cdfid, int ti)
{
  int c, i, j, k, n;
  size_t ndumps;
  size_t start[4];
  size_t count[4];

  master->dumpdata = dumpdata;
  master->data_infill = params->data_infill;
  master->i1 = NULL;
  master->i2 = NULL;
  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* time independent variables */
  count[0] = dumpdata->nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "z_grid"), start, count,
                     dumpdata->gridz);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "z_centre"), start, count,
                     dumpdata->cellz);
  /* Copy into the sparse array */
  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    geom->gridz[c] = dumpdata->gridz[k];
    geom->cellz[c] = dumpdata->cellz[k];
  }

  count[0] = dumpdata->nce2 + 1;
  count[1] = dumpdata->nce1 + 1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "x_grid"), start, count,
                     dumpdata->gridx[0]);
  c2s_2d(geom, geom->gridx, dumpdata->gridx, geom->nfe1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "y_grid"), start, count,
                     dumpdata->gridy[0]);
  c2s_2d(geom, geom->gridy, dumpdata->gridy, geom->nfe1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1agrid"), start, count,
                     dumpdata->h1agrid[0]);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2agrid"), start, count,
                     dumpdata->h2agrid[0]);
  count[0] = dumpdata->nce2;
  count[1] = dumpdata->nce1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "x_centre"), start, count,
                     dumpdata->cellx[0]);
  c2s_2d(geom, geom->cellx, dumpdata->cellx, geom->nce1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "y_centre"), start, count,
                     dumpdata->celly[0]);
  c2s_2d(geom, geom->celly, dumpdata->celly, geom->nce1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1acell"), start, count,
                     dumpdata->h1acell[0]);
  c2s_2d(geom, geom->h1acell, dumpdata->h1acell, geom->nce1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2acell"), start, count,
                     dumpdata->h2acell[0]);
  c2s_2d(geom, geom->h2acell, dumpdata->h2acell, geom->nce1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "coriolis"), start, count,
                     dumpdata->coriolis[0]);
  c2s_2d(geom, master->coriolis, dumpdata->coriolis, geom->nce1,
         geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "botz"), start, count,
                     dumpdata->botz[0]);
  if (!params->reset_bathyf)
    c2s_2d(geom, geom->botz, dumpdata->botz, geom->nce1, geom->nce2);
  count[0] = dumpdata->nce2;
  count[1] = dumpdata->nce1 + 1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "x_left"), start, count,
                     dumpdata->u1x[0]);
  c2s_2d(geom, geom->u1x, dumpdata->u1x, geom->nfe1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "y_left"), start, count,
                     dumpdata->u1y[0]);
  c2s_2d(geom, geom->u1y, dumpdata->u1y, geom->nfe1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1au1"), start, count,
                     dumpdata->h1au1[0]);
  c2s_2d(geom, geom->h1au1, dumpdata->h1au1, geom->nfe1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2au1"), start, count,
                     dumpdata->h2au1[0]);
  c2s_2d(geom, geom->h2au1, dumpdata->h2au1, geom->nfe1, geom->nce2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "thetau1"), start, count,
                     dumpdata->thetau1[0]);
  c2s_2d(geom, geom->thetau1, dumpdata->thetau1, geom->nfe1, geom->nce2);
  count[0] = dumpdata->nce2 + 1;
  count[1] = dumpdata->nce1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "x_back"), start, count,
                     dumpdata->u2x[0]);
  c2s_2d(geom, geom->u2x, dumpdata->u2x, geom->nce1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "y_back"), start, count,
                     dumpdata->u2y[0]);
  c2s_2d(geom, geom->u2y, dumpdata->u2y, geom->nce1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1au2"), start, count,
                     dumpdata->h1au2[0]);
  c2s_2d(geom, geom->h1au2, dumpdata->h1au2, geom->nce1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2au2"), start, count,
                     dumpdata->h2au2[0]);
  c2s_2d(geom, geom->h2au2, dumpdata->h2au2, geom->nce1, geom->nfe2);
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "thetau2"), start, count,
                     dumpdata->thetau2[0]);
  c2s_2d(geom, geom->thetau2, dumpdata->thetau2, geom->nce1, geom->nfe2);
  count[0] = dumpdata->nce1;
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "crci"), start, count,
                    dumpdata->crci);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "clci"), start, count,
                    dumpdata->clci);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "frci"), start, count,
                    dumpdata->frci);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "flci"), start, count,
                    dumpdata->flci);
  count[0] = dumpdata->nfe1;
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "crfi"), start, count,
                    dumpdata->crfi);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "clfi"), start, count,
                    dumpdata->clfi);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "frfi"), start, count,
                    dumpdata->frfi);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "flfi"), start, count,
                    dumpdata->flfi);
  count[0] = dumpdata->nce2;
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "cfcj"), start, count,
                    dumpdata->cfcj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "cbcj"), start, count,
                    dumpdata->cbcj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "ffcj"), start, count,
                    dumpdata->ffcj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "fbcj"), start, count,
                    dumpdata->fbcj);
  count[0] = dumpdata->nfe2;
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "cffj"), start, count,
                    dumpdata->cffj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "cbfj"), start, count,
                    dumpdata->cbfj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "fffj"), start, count,
                    dumpdata->fffj);
  nc_get_vara_short(cdfid, ncw_var_id(cdfid, "fbfj"), start, count,
                    dumpdata->fbfj);

  {
    /* Read land fill map, if available */
    int fmi = ncw_var_id(cdfid, "fmap_i");
    int fmj = ncw_var_id(cdfid, "fmap_j");
    if (fmi > -1 && fmj > -1) {
      dumpdata->fmap_i = i_alloc_3d(dumpdata->nce1,dumpdata->nce2,dumpdata->nz);
      dumpdata->fmap_j = i_alloc_3d(dumpdata->nce1,dumpdata->nce2,dumpdata->nz);

      count[0] = dumpdata->nz;
      count[1] = dumpdata->nce2;
      count[2] = dumpdata->nce1;
      nc_get_vara_int(cdfid, fmi, start, count, dumpdata->fmap_i[0][0]);
      nc_get_vara_int(cdfid, fmj, start, count, dumpdata->fmap_j[0][0]);
    }
  }

  /* time dependent variables */
  start[0] = ti;
  count[0] = 1L;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "t"), start, count,
                     &dumpdata->t);

  if (params->timeunit) {
    char timeunits[MAXSTRLEN];
    memset(timeunits, 0, MAXSTRLEN);
    nc_get_att_text(cdfid, ncw_var_id(cdfid, "t"), "units", timeunits);
    tm_change_time_units(timeunits, params->timeunit, &dumpdata->t, 1);
  }
  params->t = dumpdata->t;

  count[1] = dumpdata->nz;
  count[2] = dumpdata->nfe2;
  count[3] = dumpdata->nfe1;

  nc_get_vara_long(cdfid, ncw_var_id(cdfid, "flag"), start, count,
		   (long *)dumpdata->flag[0][0]);
  /* Adjust the flag array for sigma coordinates */
  if (params->sigma) {
    int nz = dumpdata->nz - 1;
    for (i = 0; i < dumpdata->nfe1; i++)
      for (j = 0; j < dumpdata->nfe2; j++)
	if (!(dumpdata->flag[nz][j][i] & (SOLID | OUTSIDE))) {
	  for (k = 0; k < nz; k++)
	    dumpdata->flag[k][j][i] = dumpdata->flag[nz][j][i];
	}
  }

  /* Note: use the dumpdata->u1av array to read the backward vel */
  dumpdata_read_2d(dumpdata, cdfid, "u1avb", dumpdata->u1av, ti,
                   dumpdata->nce2, dumpdata->nfe1);
  c2s_2d(geom, master->u1avb, dumpdata->u1av, geom->nfe1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "u1av", dumpdata->u1av, ti,
                   dumpdata->nce2, dumpdata->nfe1);
  data_infill_2d(master, dumpdata->u1av);
  c2s_2d(geom, master->u1av, dumpdata->u1av, geom->nfe1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "u1bot", dumpdata->u1bot, ti,
                   dumpdata->nce2, dumpdata->nfe1);
  data_infill_2d(master, dumpdata->u1bot);
  c2s_2d(geom, master->u1bot, dumpdata->u1bot, geom->nfe1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "wind1", dumpdata->wind1, ti,
                   dumpdata->nce2, dumpdata->nfe1);
  data_infill_2d(master, dumpdata->wind1);
  c2s_2d(geom, master->wind1, dumpdata->wind1, geom->nfe1, geom->nce2);
  /* Note: use the dumpdata->u1av array to read the backward vel */
  dumpdata_read_2d(dumpdata, cdfid, "u2avb", dumpdata->u2av, ti,
                   dumpdata->nfe2, dumpdata->nce1);
  c2s_2d(geom, master->u2avb, dumpdata->u2av, geom->nce1, geom->nfe2);
  dumpdata_read_2d(dumpdata, cdfid, "u2av", dumpdata->u2av, ti,
                   dumpdata->nfe2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->u2av);
  c2s_2d(geom, master->u2av, dumpdata->u2av, geom->nce1, geom->nfe2);
  dumpdata_read_2d(dumpdata, cdfid, "u2bot", dumpdata->u2bot, ti,
                   dumpdata->nfe2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->u2bot);
  c2s_2d(geom, master->u2bot, dumpdata->u2bot, geom->nce1, geom->nfe2);
  dumpdata_read_2d(dumpdata, cdfid, "wind2", dumpdata->wind2, ti,
                   dumpdata->nfe2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->wind2);
  c2s_2d(geom, master->wind2, dumpdata->wind2, geom->nce1, geom->nfe2);
  dumpdata_read_2d(dumpdata, cdfid, "wtop", dumpdata->wtop, ti,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->wtop);
  c2s_2d(geom, master->wtop, dumpdata->wtop, geom->nce1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "topz", dumpdata->topz, ti,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->topz);
  c2s_2d(geom, master->topz, dumpdata->topz, geom->nce1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "eta", dumpdata->eta, ti,
		   dumpdata->nce2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->eta);
  c2s_2d(geom, master->eta, dumpdata->eta, geom->nce1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "patm", dumpdata->patm, ti,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_2d(master, dumpdata->patm);
  c2s_2d(geom, master->patm, dumpdata->patm, geom->nce1, geom->nce2);
  dumpdata_read_2d(dumpdata, cdfid, "Cd", dumpdata->Cd, ti, dumpdata->nce2,
                   dumpdata->nce1);
  data_infill_2d(master, dumpdata->Cd);
  c2s_2d(geom, master->Cd, dumpdata->Cd, geom->nce1, geom->nce2);
  /* Note: use the dumpdata->u1 array to read the backward vel */
  dumpdata_read_3d(dumpdata, cdfid, "u1b", dumpdata->u1, ti, dumpdata->nz,
                   dumpdata->nce2, dumpdata->nfe1);
  c2s_3d(geom, master->u1b, dumpdata->u1, geom->nfe1, geom->nce2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "u2b", dumpdata->u2, ti, dumpdata->nz,
                   dumpdata->nfe2, dumpdata->nce1);
  c2s_3d(geom, master->u2b, dumpdata->u2, geom->nce1, geom->nfe2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "u1", dumpdata->u1, ti, dumpdata->nz,
                   dumpdata->nce2, dumpdata->nfe1);
  data_infill_3d(master, dumpdata->u1);
  c2s_3d(geom, master->u1, dumpdata->u1, geom->nfe1, geom->nce2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "u2", dumpdata->u2, ti, dumpdata->nz,
                   dumpdata->nfe2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->u2);
  c2s_3d(geom, master->u2, dumpdata->u2, geom->nce1, geom->nfe2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "w", dumpdata->w, ti, dumpdata->nz,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->w);
  c2s_3d(geom, master->w, dumpdata->w, geom->nce1, geom->nce2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "dens", dumpdata->dens, ti,
                   dumpdata->nz, dumpdata->nce2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->dens);
  c2s_3d(geom, master->dens, dumpdata->dens, geom->nce1, geom->nce2,
         geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "dens_0", dumpdata->dens_0, ti,
                   dumpdata->nz, dumpdata->nce2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->dens_0);
  c2s_3d(geom, master->dens_0, dumpdata->dens_0, geom->nce1, geom->nce2,
         geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "Kz", dumpdata->Kz, ti, dumpdata->nz,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->Kz);
  c2s_3d(geom, master->Kz, dumpdata->Kz, geom->nce1, geom->nce2, geom->nz);
  dumpdata_read_3d(dumpdata, cdfid, "Vz", dumpdata->Vz, ti, dumpdata->nz,
                   dumpdata->nce2, dumpdata->nce1);
  data_infill_3d(master, dumpdata->Vz);
  c2s_3d(geom, master->Vz, dumpdata->Vz, geom->nce1, geom->nce2, geom->nz);

  /* Check for NaNs                                                  */
  for (i = 1; i <= geom->b2_t; i++) {
    int ci, cj;
    c = geom->w2_t[i];
    if (isnan(master->eta[c])) {
      find_closest_nonnan(dumpdata, dumpdata->eta, geom->s2i[c], geom->s2j[c],
			  geom->nz-1, &ci, &cj);
      master->eta[c] = dumpdata->eta[cj][ci];
      hd_warn("eta_init: Found NaN at (%d %d); replaced with (%d %d).\n", 
	      geom->s2i[c], geom->s2j[c], ci, cj);
    }
  }

  for (n = 0; n < dumpdata->ntr; n++) {
    /* If a 3D tracer is included in the input file then initialise  */
    /* with this data. If not use the specification in the parameter */
    /* file tracer list; TRACER#.data                                */
    if (dumpdata_read_3d_tr(dumpdata, cdfid, dumpdata->trinfo_3d[n].name,
			    dumpdata->tr_wc[n], ti, dumpdata->nz, 
			    dumpdata->nce2, dumpdata->nce1)) {
      load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_3d[n].name, WATER);
    } else {
      data_infill_3d(master, dumpdata->tr_wc[n]);
      c2s_3d(geom, master->tr_wc[n], dumpdata->tr_wc[n], geom->nce1,
      geom->nce2, geom->nz);
    }
    /* Ignore DA variables */
    if (strncasecmp(master->trinfo_3d[n].tag, "DA_OBS", 6) == 0 ) continue;

    /* Check for NaN's in the initial condition */
    for (i = 1; i <= geom->b3_t; i++) {
      c = geom->w3_t[i];
      /*
      if (master->trinfo_3d[n].valid_range_wc[0] >= 0.0 &&
	  master->tr_wc[n][c] == 0.0) {
	hd_warn
	  ("dumpdata_read: Zero found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
	master->tr_wc[n][c] = master->trinfo_3d[n].fill_value_wc;
      }
      */
      if (isnan(master->tr_wc[n][c]))
	hd_warn
	  ("dumpdata_read: NaN found in tracer '%s' at (%d %d %d).\n",
	   master->trinfo_3d[n].name, geom->s2i[c], geom->s2j[c], geom->s2k[c]);
      master->tr_wc[n][c] = min(master->tr_wc[n][c], master->trinfo_3d[n].valid_range_wc[1]);
    }
  }

  for (n = 0; n < dumpdata->ntrS; n++) {
    if (dumpdata_read_2d(dumpdata, cdfid, dumpdata->trinfo_2d[n].name,
                     dumpdata->tr_wcS[n], ti, dumpdata->nce2,
			 dumpdata->nce1)) {
      load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_2d[n].name, INTER);
    } else {
      data_infill_2d(master, dumpdata->tr_wcS[n]);
      c2s_2d(geom, master->tr_wcS[n], dumpdata->tr_wcS[n], geom->nce1,
	     geom->nce2);
    }
  }

  if (dumpdata->sednz > 0) {
    int flg = 0;
    if (dumpdata_read_3d(dumpdata, cdfid, "z_grid_sed", dumpdata->gridz_sed,
			 ti, dumpdata->sednz + 1, dumpdata->nce2,
			 dumpdata->nce1)) flg = 1;
    for (k = 0; k < params->sednz + 1; k++) {
      c2s_2d(geom, geom->gridz_sed[k], dumpdata->gridz_sed[k], geom->nce1,
             geom->nce2);
    }
    if (dumpdata_read_3d(dumpdata, cdfid, "z_centre_sed", dumpdata->cellz_sed,
			 ti, dumpdata->sednz, dumpdata->nce2, dumpdata->nce1))
      flg = 1;
    for (k = 0; k < params->sednz; k++) {
      c2s_2d(geom, geom->cellz_sed[k], dumpdata->cellz_sed[k], geom->nce1,
             geom->nce2);
    }
    /* If sediment structure does not exist in the input file, then  */
    /* use the structure defined in the parameter file.              */
    if (flg) read_sed_layers(geom, params, master, dumpdata);
		     
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      if (dumpdata_read_3d(dumpdata, cdfid, name, dumpdata->tr_sed[n], ti,
			   dumpdata->sednz, dumpdata->nce2, dumpdata->nce1)) {
	load_wc_tracer_name(master, params->prmfd, dumpdata->trinfo_sed[n].name, SEDIM);
      } else {
	for (k = 0; k < params->sednz; k++)
	  c2s_2d(geom, master->tr_sed[n][k], dumpdata->tr_sed[n][k],
		 geom->nce1, geom->nce2);
      }
    }
  }

  /* If bathymetry is changed with reset_bathyf (smoothed, masked)  */
  /* when using -p mode, then the input data should be cascade      */
  /* filled onto the updated grid using data_infill. dumpdata->botz */
  /* contains the original bathymetry from the INPUT_FILE at this   */
  /* stage, hence discrepencies in bathymetry can be found by       */
  /* comparing dumdata->botz and params->topo (the latter contains  */
  /* the modified bathymetry). Reset the dumpdata bathymetry after  */
  /* arrays have been filled.                                       */
  if(params->data_infill) {
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	dumpdata->botz[j][i] = params->topo[j][i];
    i_free_3d(master->i1);
    i_free_3d(master->i2);
  }
  read_mean_atts(master, cdfid);
  return (1);
}

/* END dumpdata_read()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the input data into the dump data structure and   */
/* copy to the sparse array.                                         */
/*-------------------------------------------------------------------*/
double **bathy_read(parameters_t *params, /* Input parameter data structure */
	       int cdfid, int ti)
{
  double **bathy;               /* Bathymetry array */
  int i,j;

  size_t ndumps;
  size_t start[4];
  size_t count[4];

  /* Check that the requested dump exists */
  nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "record"), &ndumps);
  if (ti > (int)(ndumps - 1)) {
    warn("dump_read: dump %d does not exist in file!\n", ti);
    return (0);
  }

  /* Allocate memory */
  bathy = d_alloc_2d(params->nce1, params->nce2);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  /* time independent variables */
  count[2] = 0;
  count[3] = 0;

  if (DEBUG("init_m"))
    dlog("init_m", "\nReading bathymetry from input file\n\n");

  /* Read the bathymetry */
  count[0] = params->nce2;
  count[1] = params->nce1;
  nc_get_vara_double(cdfid, ncw_var_id(cdfid, "botz"), start, count,
                     bathy[0]);



  for (j = 0; j < params->nce2; j++) {
    for (i = 0; i < params->nce1; i++) {
      if(isnan(bathy[j][i]))
	bathy[j][i] = NOTVALID;
    }
  }

  /* Check for uniform bathymetry over zoomed grids */
  if(params->dozoom) {
    int n, ii, jj;
    double d1;
    for (n = 1; n <= params->nzoom; n++) {
      i = params->zci[n] - params->zmfe1 / 2;
      j = params->zcj[n] - params->zmfe2 / 2;
      d1 = bathy[params->zcj[n]][params->zci[n]];
      for (ii = i; ii < i + params->zmfe1; ii++)
	for (jj = j; jj < j + params->zmfe2; jj++)
	  if (bathy[jj][ii] != d1) {
	    hd_warn("Zoomed cells must have uniform bathymetry. Remake the input file.\n");
	    hd_quit("Invalid bathymetry in zoomed cell %4.1f@(%d %d) should be %4.1f.\n", bathy[jj][ii], ii, jj, d1);
	  }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set boundary OUTSIDE OBC cells if required                      */
  adjust_outside_obc(params, bathy);
  return (bathy);
}

/* END bathy_read()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 2d array from id                                           */
/*-------------------------------------------------------------------*/
int dumpdata_read_2d(dump_data_t *dumpdata, int id, char *name,
		     double **p, int dump, int nj, int ni)
{
  size_t start[4];
  size_t count[4];
  int vid;

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nj;
  count[2] = ni;
  count[3] = 0;
  
  vid = ncw_var_id(id, name);
  if (vid < 0) {
    zero2d(dumpdata, p, ni, nj);
    return(1);
  }
  nc_get_vara_double(id, vid, start, count, p[0]);
  return(0);
}

/* END dumpdata_read_2d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 3d array from id                                           */
/*-------------------------------------------------------------------*/
int dumpdata_read_3d(dump_data_t *dumpdata, int id, char *name,
		     double ***p, int dump, int nk, int nj, int ni)
{
  size_t start[4];
  size_t count[4];
  int vid;

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nk;
  count[2] = nj;
  count[3] = ni;

  vid = ncw_var_id(id, name);
  if (vid < 0) {
    int k;
    for (k = 0; k < nk; k++)
      zero2d(dumpdata, p[k], ni, nj);
    return(1);
  }
  nc_get_vara_double(id, vid, start, count, p[0][0]);
  return(0);
}

/* END dumpdata_read_3d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read a 3d array from id                                           */
/*-------------------------------------------------------------------*/
int dumpdata_read_3d_tr(dump_data_t *dumpdata, int id, char *name,
			 double ***p, int dump, int nk, int nj, int ni)
{
  size_t start[4];
  size_t count[4];
  int vid;

  start[0] = dump;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nk;
  count[2] = nj;
  count[3] = ni;

  vid = ncw_var_id(id, name);
  if (vid < 0)
    return(1);

  nc_get_vara_double(id, vid, start, count, p[0][0]);
  return(0);
}

/* END dumpdata_read_3d_tr()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Zeros a 2d array                                                  */
/*-------------------------------------------------------------------*/
void zero2d(dump_data_t *dumpdata, double **dst, int n1, int n2)
{
  /* The fast way */
  memset((void *)dst[0], 0, sizeof(double) * n1 * n2);

  if (dst[0][0] != 0.0 || dst[n2 - 1][n1 - 1] != 0.0) {
    int i, j;
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++)
        dst[j][i] = 0.0;
  }
}

/* END zero2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Zeros a 2d array                                                  */
/*-------------------------------------------------------------------*/
void zero3d(dump_data_t *dumpdata, double ***dst, int n1, int n2, int n3)
{
  int k;

  for (k = 0; k < n3; k++)
    zero2d(dumpdata, dst[k], n1, n2);
}

/* END zero3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
double **set_bathy(parameters_t *params,  /* Input parameter data
                                             structure */
                   int cdfid    /* Input netCDF file handle */
  )
{
  int i, j, k, n;               /* Counters */
  int is, ie, js, je;           /* Limits of grid */
  int xsize, ysize;             /* Size of original grid */
  int bathylimit = 0;           /* Flag to set min and max bathymetry */
  int percent = 0;              /* Flag to set maximum elevation */
  double min_cell_thickness = 0.0;  /* Minimum cell thickness */
  unsigned long topflag = 0;    /* Surface flag value */
  double bmin, bmax;            /* Minimum and maximum bathymetry */
  double etamax;                /* Maximum surface elevation */
  double *mincelldz;            /* Minimum layer thickness */
  unsigned long ***flag;        /* Flag for Cartesian grid */
  double *layers;               /* Layer spacing array */
  double maxgrad;               /* Maximum allowable bottom gradient */
  double mdh;                   /* Maximum bottom gradient found in file
                                   topo */
  double d1, d2;                /* Dummies */
  double dist = 0.0;            /* Grid spacing over maximum bottom
                                   gradient */
  int ic = 0, jc = 0;           /* Grid position at maximum bottom
                                   gradient */
  int nps = 0;                  /* Smoothing passes to reduce bottom
                                   gradient */
  int **xp, **yp, **xm, **ym;
  int ip, jp, im, jm;
  int ipp, jpp, imm, jmm;
  double **h1au1, **h1au2;
  double **h2au2, **h2au1;
  double **bathy;               /* Bathymetry array */
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int nz = params->nz;
  int sf = 1;
  int maskf = 0;
  char dir = '\0';

  /* Allocate memory */
  bathy = d_alloc_2d(nce1, nce2);
  layers = params->layers;
  flag = (unsigned long ***)l_alloc_3d(nce1 + 1, nce2 + 1, nz);
  h1au1 = d_alloc_2d(nce1 + 1, nce2);
  h2au1 = d_alloc_2d(nce1 + 1, nce2);
  h2au2 = d_alloc_2d(nce1, nce2 + 1);
  h1au2 = d_alloc_2d(nce1, nce2 + 1);
  maxgrad = params->maxgrad;

  /*-----------------------------------------------------------------*/
  /* Read parameters required for adjusting the bathymetry */
  is = js = 0;
  ie = xsize = nce1;
  je = ysize = nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = nce1 - params->edgef;
    je = nce2 - params->edgef;
    xsize = nce1 - 2 * params->edgef;
    ysize = nce2 - 2 * params->edgef;
  }

  etamax = params->etamax;
  min_cell_thickness = atof(params->mct);
  if (params->mct[strlen(params->mct) - 1] == '%')
    percent = 1;
  mincelldz = d_alloc_1d(nz);
  for (k = 0; k < nz; k++) {
    if (percent)
      mincelldz[k] =
        (layers[k + 1] - layers[k]) * min_cell_thickness / 100.0;
    else
      mincelldz[k] = min_cell_thickness;
  }

  /*-----------------------------------------------------------------*/
  /* Read the bathymetry array */
  bmin = bmax = 0.0;
  if (params->bmin)
    bmin = params->bmin;
  if (params->bmax)
    bmax = params->bmax;
  if (bmin && bmax) {
    if (bmin > 0 && bmax > 0 && bmin > bmax)
      hd_quit("ERROR: set_bathy: BATHYMIN (%4.2f) > BATHYMAX (%4.2f)\n", bmin, bmax);
    else
      bathylimit = 1;
  }

  if (params->reset_bathyf) {
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
	bathy[j][i] = params->topo[j][i];
  } else {
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
	bathy[j][i] = LANDCELL;

    if (params->nvals == 0) {
      /* Use default bathymetry */
      for (j = js; j < je; j++)
	for (i = is; i < ie; i++)
	  bathy[j][i] = layers[0];
    } else if (params->nvals == (xsize) * (ysize)) {
      /* Use bathymetry as read */
      n = 0;
      for (j = js; j < je; j++)
	for (i = is; i < ie; i++)
	  bathy[j][i] = -params->bathy[n++];
    } else
      hd_quit
	("set_bathy: bad bathymetry data (neither 0 nor number of cells)\n");

    /*---------------------------------------------------------------*/
    /* Set boundary OUTSIDE OBC cells if required                    */
    adjust_outside_obc(params, bathy);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise the flag array */
  for (j = 0; j <= nce2; j++)
    for (i = 0; i <= nce1; i++)
      for (k = 0; k < nz; k++) {
        flag[k][j][i] = 0;
      }

  /*-----------------------------------------------------------------*/
  /* Set the flags to solid on the nce1 and nce2 edges */
  for (k = 0; k < nz; k++) {
    i = nce1;
    for (j = 0; j < nce2; j++) {
      if (bathy[j][i - 1] >= layers[k + 1] || bathy[j][i - 1] >= etamax)
        flag[k][j][i] |= SOLID;
      else
        flag[k][j][i] |= U1OUTSIDE;
    }
    j = nce2;
    for (i = 0; i < nce1; i++) {
      if (bathy[j - 1][i] >= layers[k + 1] || bathy[j - 1][i] >= etamax)
        flag[k][j][i] |= SOLID;
      else
        flag[k][j][i] |= U2OUTSIDE;
    }
  }

  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {

      /* Adjust botz to min and max bathymetry, checking for outside */
      /* cells (botz < gridz[0]) and land cells (botz > etamax).  */
      if (bathylimit && bathy[j][i] < etamax && bathy[j][i] > -bmin)
        bathy[j][i] = -bmin;
      if (bathylimit && bathy[j][i] >= layers[0] && bathy[j][i] < -bmax)
        bathy[j][i] = -bmax;

      /* Mark outside cells */
      if (bathy[j][i] < layers[0]) {
        flag[nz - 1][j][i] |= OUTSIDE;
      }

      /* Load the cell values */
      topflag = flag[nz - 1][j][i];
      for (k = 0; k < nz; k++) {
        if (topflag & OUTSIDE) {
          flag[k][j][i] |= OUTSIDE;
        } else {
          /* Don't allow cells to be less than MIN_CELL_THICKNESS */
          if (bathy[j][i] < layers[k + 1] &&
              bathy[j][i] > layers[k + 1] - mincelldz[k]) {
            double newbotz;
            if ((layers[k + 1] - bathy[j][i]) >
                (layers[k + 1] - layers[k]) / 2)
              newbotz = layers[k];
            else
              newbotz = layers[k + 1];
            bathy[j][i] = newbotz;
          }
          if (bathy[j][i] >= layers[k + 1] || bathy[j][i] >= etamax) {
            /* A solid cell */
            flag[k][j][i] |= SOLID;
          }
        }
      }
    }
  }

  /* Get the horizontal grid spacings */
  if (params->runmode & (AUTO | DUMP)) {
    for (j = 0; j < params->nce2 + 1; j++) {
      for (i = 0; i < params->nce1; i++) {
        h2au2[j][i] = 2.0 * params->h2[j * 2][i * 2 + 1];
        h1au2[j][i] = 2.0 * params->h1[j * 2][i * 2 + 1];
      }
    }
    for (j = 0; j < params->nce2; j++) {
      for (i = 0; i < params->nce1 + 1; i++) {
        h1au1[j][i] = 2.0 * params->h1[j * 2 + 1][i * 2];
        h2au1[j][i] = 2.0 * params->h2[j * 2 + 1][i * 2];
      }
    }
  } else {
    size_t start[4];
    size_t count[4];

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    count[0] = nce2;
    count[1] = nce1 + 1;
    count[2] = 0;
    count[3] = 0;
    nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h1au1"), start, count,
                       h1au1[0]);
    count[0] = nce2 + 1;
    count[1] = nce1;
    nc_get_vara_double(cdfid, ncw_var_id(cdfid, "h2au2"), start, count,
                       h2au2[0]);
  }

  /* Clean up isolated cells */
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	  ip = (i < nce1 - 1) ? i+1 : i;
	  im = (i > 0) ? i-1 : i;
	  jp = (j < nce2 - 1) ? j+1 : j;
	  jm = (j > 0) ? j-1 : j;
	  if (flag[nz - 1][j][im] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][j][ip] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][jm][i] & (SOLID | OUTSIDE) &&
	      flag[nz - 1][jp][i] & (SOLID | OUTSIDE)) {
	    bathy[j][i] = LANDCELL;
	    for (k = 0; k < nz; k++) {
	      flag[k][j][i] |= SOLID;
	    }
	    hd_warn("Filling in isolated wet cell at (%d %d)\n",i, j);
	  }
	}
      }
    }

    /*remove_channel(params, flag, bathy);*/

  /*-----------------------------------------------------------------*/
  /* Set any bathymetry masking                                      */
  d1 = bmin;
  is = ie = js = je = -1;
  prm_set_errfn(hd_silent_warn);
  if (prm_read_int(params->prmfd, "BATHY_MASK_IS", &is)) {
    prm_set_errfn(hd_quit);
    prm_read_int(params->prmfd, "BATHY_MASK_IE", &ie);
    prm_read_double(params->prmfd, "BATHY_MASK_DS", &d1);
    prm_read_double(params->prmfd, "BATHY_MASK_DE", &d2);
    d1 *= -1.0; d2 *= -1.0;
    prm_set_errfn(hd_silent_warn);
    maskf = 1;
  } else if (prm_read_int(params->prmfd, "BATHY_MASK_JS", &js)) {
    prm_set_errfn(hd_quit);
    prm_read_int(params->prmfd, "BATHY_MASK_JE", &je);
    prm_read_double(params->prmfd, "BATHY_MASK_DS", &d1);
    prm_read_double(params->prmfd, "BATHY_MASK_DE", &d2);
    d1 *= -1.0; d2 *= -1.0;
    prm_set_errfn(hd_silent_warn);
    maskf = 2;
  }

  /* Single block - gradient or single value */
  mdh = 0.0;
  if (prm_read_double(params->prmfd, "BATHY_MASK_VAL", &mdh)) {
    maskf = 3;
    mdh *= -1.0;
  }
  if (maskf) {
    int *bmi, *bmj;
    read_blocks(params->prmfd, "BATHY_MASK", &n, &bmi, &bmj);
    for (k = 1; k <= n; k++) {
      i = bmi[k];
      j = bmj[k];
      if (maskf == 1)
	bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
      else if (maskf == 2)
	bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
      else 
	bathy[j][i] = mdh;
    }
  }

  /* Multiple blocks - single value only */
  if (prm_read_int(params->prmfd, "BATHY_NBLOCKS", &im)) {
    char buf[MAXSTRLEN];
    int *bmi, *bmj;
    for (ip = 0; ip < im; ip++) {
      mdh = 0.0;
      sprintf(buf, "BATHY_MASK_VAL%d", ip);
      if (prm_read_double(params->prmfd, buf, &mdh))
	mdh *= -1.0;
      sprintf(buf, "BATHY_MASK%d", ip);
      read_blocks(params->prmfd, buf, &n, &bmi, &bmj);
      for (k = 1; k <= n; k++) {
	i = bmi[k];
	j = bmj[k];
	bathy[j][i] = mdh;
      }
    }
  }

  /*
  if (prm_read_int(params->prmfd, "BATHY_MASK", &n)) {
    int i1, i2, j1, j2;
    prm_flush_line(params->prmfd);
    if (fscanf(params->prmfd, "(%d,%d)-(%d,%d)", &i1, &i2, &j1, &j2) == 4) {
      prm_flush_line(params->prmfd);
      for (j = 0; j < params->nce2; j++) {
	for (i = 0; i < params->nce1; i++) {
	  if (i >= i1 && i <= i2 && j >= j1 && j <= j2) {
	    bathy[j][i] = mdh;
	    if (i >= is && i <= ie)
	      bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
	    if (j >= js && j <= je)
	      bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
	  }
	}
      }
    }
  }

  if (prm_read_int(params->prmfd, "BATHY_MASK", &n)) {
    prm_flush_line(params->prmfd);
    for (k = 0; k < n; ++k) {
      if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	hd_quit("set_bathy: Can't read i j in BATHY_MASK points list.\n");
      prm_flush_line(params->prmfd);
      bathy[j][i] = mdh;
      if (i >= is && i <= ie)
	bathy[j][i] = (d2 - d1) * (double)(i - is) / (double)(ie - is) + d1;
      if (j >= js && j <= je)
	bathy[j][i] = (d2 - d1) * (double)(j - js) / (double)(je - js) + d1;
    }
  }
  */
  /*-----------------------------------------------------------------*/
  /* If the zoom points were exclusive (zoomf < 0) then reset the    */
  /* zoom points as all wet cells not listed.                        */
  /* Not implemented yet.
  if (params->zoomf < 0) {
    int nzp, *zci, *zcj;
    nzp = params->nzoom;
    zci = i_alloc_1d(nzp + 1);
    zcj = i_alloc_1d(nzp + 1);
    memcpy(zci, params->zci, (nzp + 1) * sizeof(int));
    memcpy(zcj, params->zcj, (nzp + 1) * sizeof(int));
    i_free_1d(params->zci);
    i_free_1d(params->zcj);
    params->zoomf = -params->zoomf;
    params->nzoom = 0;
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
        if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) params->nzoom++;
    params->nzoom -= nzp;
    params->zci = i_alloc_1d(params->nzoom + 1);
    params->zcj = i_alloc_1d(params->nzoom + 1);
    n = 0;
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
        if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	  if (!ANY(i, zci, nzp) && !ANY(j, zcj, nzp)) {
	    n++;
	    params->zci[n] = i;
	    params->zcj[n] = j;
	  }
	}
  }
  */

  /* Set up mapping arrays */
  xp = i_alloc_2d(nce1, nce2);
  xm = i_alloc_2d(nce1, nce2);
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (i < nce1 - 1)
	xp[j][i] =
	  ((flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE) ||
	    i == nce1 - 1) ? i : i + 1);
      else
	xp[j][i] = i;
      if (i > 0)
	xm[j][i] =
	  ((flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE) ||
	    i == 0) ? i : i - 1);
      else
	xm[j][i] = i;
    }
  }
  for (i = 0; i < params->nemx; i++) {
    xp[params->emjsx[i]][params->emisx[i]] = params->emidx[i];
    xm[params->emjdx[i]][params->emidx[i]] = params->emisx[i];
  }

  yp = i_alloc_2d(nce1, nce2);
  ym = i_alloc_2d(nce1, nce2);
  for (i = 0; i < nce1; i++) {
    for (j = 0; j < nce2; j++) {
      if (j < nce2 - 1)
	yp[j][i] =
	  ((flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE) ||
	    j == nce2 - 1) ? j : j + 1);
      else
	yp[j][i] = j;
      if (j > 0)
	ym[j][i] =
	  ((flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE) ||
	    j == 0) ? j : j - 1);
      else
	ym[j][i] = j;
    }
  }
  for (j = 0; j < params->nemy; j++) {
    yp[params->emjsy[j]][params->emisy[j]] = params->emidy[j];
    ym[params->emjdy[j]][params->emidy[j]] = params->emisy[j];
  }

  /*-----------------------------------------------------------------*/
  /* Smooth the bathymetry with a convolution filter if required.    */
  /* Smoothing over a selected area.                                 */
  if(params->smooth) {
    int m;
    if (prm_read_int(params->prmfd, "SMOOTH_MASK", &m)) {
      double **nbthy = d_alloc_2d(nce1, nce2);
      double kk[4][4];
      /* Set the filter weights */
      kk[1][3] = 0.0625;
      kk[2][3] = 0.125;
      kk[3][3] = 0.0625;
      kk[1][2] = 0.125;
      kk[2][2] = 0.25;
      kk[3][2] = 0.125;
      kk[1][1] = 0.0625;
      kk[2][1] = 0.125;
      kk[3][1] = 0.0625;
      n = params->smooth; 
      while(n) {
	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	    hd_quit("set_bathy: Can't read i j in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	    im = xm[j][i];
	    ip = xp[j][i];
	    jm = ym[j][i];
	    jp = yp[j][i];
	    /* Perform the convolution */
	    nbthy[j][i] = kk[1][3] * min(bathy[jp][im], 0.0) + 
	                  kk[2][3] * min(bathy[jp][i], 0.0) + 
	                  kk[3][3] * min(bathy[jp][ip], 0.0) +
	                  kk[1][2] * min(bathy[j][im], 0.0) + 
	                  kk[2][2] * min(bathy[j][i], 0.0) + 
	                  kk[3][2] * min(bathy[j][ip], 0.0) +
                          kk[1][1] * min(bathy[jm][im], 0.0) + 
	                  kk[2][1] * min(bathy[jm][i], 0.0) + 
	                  kk[3][1] * min(bathy[jm][ip], 0.0);
	    nbthy[j][i] = min(nbthy[j][i], -bmin);
	    nbthy[j][i] = max(nbthy[j][i], -bmax);
	  }
	}
	prm_flush_line(params->prmfd);
	prm_read_int(params->prmfd, "SMOOTH_MASK", &m);
	for (k = 0; k < m; ++k) {
	  if (fscanf(params->prmfd, "%d %d", &i, &j) != 2)
	    hd_quit("set_bathy: Can't read i j in SMOOTH_MASK points list.\n");
	  prm_flush_line(params->prmfd);
	  bathy[j][i] = nbthy[j][i];      
	}
	n--;
      }
      d_free_2d(nbthy);
    } else {
      /* Smoothing over a the whole domain.                          */
      double **nbthy = d_alloc_2d(nce1, nce2);
      double kk[4][4];
      n = params->smooth;

      /* Set the filter weights */
      kk[1][3] = 0.0625;
      kk[2][3] = 0.125;
      kk[3][3] = 0.0625;
      kk[1][2] = 0.125;
      kk[2][2] = 0.25;
      kk[3][2] = 0.125;
      kk[1][1] = 0.0625;
      kk[2][1] = 0.125;
      kk[3][1] = 0.0625;

      while(n) {
	for (j = 0; j < nce2; j++) {
	  for (i = 0; i < nce1; i++) {
	    if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	      double k0 = kk[2][2];
	      double k1 = kk[1][2], k2 = kk[3][2];
	      double k3 = kk[1][1], k4 = kk[2][1], k5 = kk[3][1];
	      double k6 = kk[1][3], k7 = kk[2][3], k8 = kk[3][3];
	      if (flag[nz - 1][j][i-1] & (SOLID | OUTSIDE)) {
		k0 += k1;
		k1 = 0.0;
	      }
	      if (flag[nz - 1][j][i+1] & (SOLID | OUTSIDE)) {
		k0 += k2;
		k2 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i-1] & (SOLID | OUTSIDE)) {
		k0 += k3;
		k3 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i] & (SOLID | OUTSIDE)) {
		k0 += k4;
		k4 = 0.0;
	      }
	      if (flag[nz - 1][j-1][i+1] & (SOLID | OUTSIDE)) {
		k0 += k5;
		k5 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i-1] & (SOLID | OUTSIDE)) {
		k0 += k6;
		k6 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i] & (SOLID | OUTSIDE)) {
		k0 += k7;
		k7 = 0.0;
	      }
	      if (flag[nz - 1][j+1][i+1] & (SOLID | OUTSIDE)) {
		k0 += k8;
		k8 = 0.0;
	      }
	      im = xm[j][i];
	      ip = xp[j][i];
	      jm = ym[j][i];
	      jp = yp[j][i];
	      /* Perform the convolution */
	      nbthy[j][i] = k6 * bathy[jp][im] + 
	                    k7 * bathy[jp][i] + 
	                    k8 * bathy[jp][ip] +
	                    k1 * bathy[j][im] + 
	                    k0 * bathy[j][i] + 
	                    k2 * bathy[j][ip] +
                            k3 * bathy[jm][im] + 
	                    k4 * bathy[jm][i] + 
	                    k5 * bathy[jm][ip];
	      nbthy[j][i] = min(nbthy[j][i], -bmin);
	      nbthy[j][i] = max(nbthy[j][i], -bmax);
	    }
	  }
	}
	for (j = 0; j < nce2; j++) {
	  for (i = 0; i < nce1; i++) {
	    if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
	      bathy[j][i] = nbthy[j][i];
	    }
	  }
	}
	n--;
      }
      d_free_2d(nbthy);
    }
  }

  /* Check that the bathymetry gradient is OK */
  while (maxgrad > 0 && sf) {
    mdh = 0.0;
    sf = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
        if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
          im = xm[j][i];
          ip = xp[j][i];
          jm = ym[j][i];
          jp = yp[j][i];
          /* Backward x gradient */
          d1 = fabs(bathy[j][i] - bathy[j][im]) / h1au1[j][i];
          d2 = fabs(bathy[j][ip] - bathy[j][i]) / h1au1[j][ip];
          if (d1 > maxgrad || d2 > maxgrad) {
            ipp = xp[j][ip];
            imm = xm[j][im];
            if (params->runmode & (AUTO | DUMP)) {
              bathy[j][i] =
                smoothbathy(bathy, i, ip, im, ipp, imm, j, 1, bmin);
              if (!(flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE)))
                bathy[j][ip] =
                  smoothbathy(bathy, ip, ipp, i, xp[j][ipp], im, j, 1,
                              bmin);
              if (!(flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE)))
                bathy[j][im] =
                  smoothbathy(bathy, im, i, imm, ip, xm[j][imm], j, 1,
                              bmin);
            }
            if (d1 > mdh) {
              mdh = d1;
              ic = i;
              jc = j;
              dist = h1au1[j][i];
              dir = 'x';
            }
          }
          /* Backward y gradient */
          d2 = fabs(bathy[j][i] - bathy[jm][i]) / h2au2[j][i];
          d1 = fabs(bathy[jp][i] - bathy[j][i]) / h2au2[jp][i];
          if (d1 > maxgrad || d2 > maxgrad) {
            jpp = yp[jp][i];
            jmm = ym[jm][i];
            if (params->runmode & (AUTO | DUMP)) {
              bathy[j][i] =
                smoothbathy(bathy, j, jp, jm, jpp, jmm, i, 0, bmin);
              if (!(flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE)))
                bathy[jp][i] =
                  smoothbathy(bathy, jp, jpp, j, yp[jpp][i], jm, i, 0,
                              bmin);
              if (!(flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE)))
                bathy[jm][i] =
                  smoothbathy(bathy, jm, j, jmm, jp, ym[jmm][i], i, 0,
                              bmin);
            }
            sf = 1;
            if (d2 > mdh) {
              mdh = d2;
              ic = i;
              jc = j;
              dist = h2au2[j][i];
              dir = 'y';
            }
          }
        }
      }
    }
    nps++;
    if (mdh != 0.0) {
      printf
        ("Warning : Pass %d %c maximum bottom topography gradient too large : %4.3f\n",
         nps, dir, mdh);
      printf("%6.2fm in %8.2fm at point (%d,%d)\n", mdh * dist, dist, ic,
             jc);
    }
    if (nps > 10)
      break;
  }

  /* Check that uniform flow over topography will not violate the */
  /* vertical velocity Courant number.  */
  if (params->bvel > 0.0) {
    double ***dz, w, w1, w2, ***u1flux, ***u2flux, top, bot;
    double fcbotx, fctopx, fcboty, fctopy;
    short **kb;
    double vel = params->bvel;
    dz = d_alloc_3d(nce1, nce2, nz);
    u1flux = d_alloc_3d(nce1 + 1, nce2, nz);
    u2flux = d_alloc_3d(nce1, nce2 + 1, nz);
    kb = s_alloc_2d(nce1, nce2);
    sf = 1;
    nps = 0;

    while (sf) {
      sf = 0;
      /* Get the cell thickness */
      for (j = 0; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          kb[j][i] = nz;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            kb[j][i] = 1;
            while (kb[j][i] < nz && layers[kb[j][i]] < bathy[j][i])
              kb[j][i]++;
            kb[j][i]--;
            for (k = kb[j][i]; k < nz; k++) {
              top = (k == nz - 1) ? 0 : layers[k + 1];
              bot = (k == kb[j][i]) ? bathy[j][i] : layers[k];
              dz[k][j][i] = top - bot;
            }
          }
        }
      }

      /* Get the fluxes */
      for (j = 0; j < nce2; j++) {
        for (i = 1; i < nce1; i++) {
          for (k = 0; k < nz; k++)
            u1flux[k][j][i] = 0.0;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            for (k = kb[j][i]; k < nz; k++) {
              if (!(flag[k][j][i] & (SOLID | OUTSIDE)) &&
                  !(flag[k][j][i - 1] & (SOLID | OUTSIDE))) {
                u1flux[k][j][i] = vel * dz[k][j][i] * h2au1[j][i];
              }
            }
          }
        }
      }
      for (j = 1; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          for (k = 0; k < nz; k++)
            u2flux[k][j][i] = 0.0;
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            for (k = kb[j][i]; k < nz; k++) {
              if (!(flag[k][j][i] & (SOLID | OUTSIDE)) &&
                  !(flag[k][j - 1][i] & (SOLID | OUTSIDE))) {
                u2flux[k][j][i] = vel * dz[k][j][i] * h1au2[j][i];
              }
            }
          }
        }
      }

      /* Get the vertical velocity */
      for (j = 0; j < nce2; j++) {
        for (i = 0; i < nce1; i++) {
          if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
            /* x velocity */
            w1 = w2 = fcbotx = fcboty = 0.0;
            for (k = kb[j][i]; k < nz - 1; k++) {

              fctopx = fcbotx + u1flux[k][j][i] - u1flux[k][j][i + 1];
              w1 = fabs(fctopx / (h2au1[j][i] * h1au2[j][i]));
              fcbotx = fctopx;
              fctopy = fcboty + u2flux[k][j][i] - u2flux[k][j + 1][i];
              w2 = fabs(fctopy / (h2au1[j][i] * h1au2[j][i]));
              fcboty = fctopy;
              w = max(w1, w2);
              d1 = w * params->grid_dt / dz[k][j][i];
              if (d1 > 1.0) {
                sf = 1;
                if (params->runmode == AUTO) {

                  if (!(flag[nz - 1][j][i + 1] & (SOLID | OUTSIDE)) &&
                      !(flag[nz - 1][j][i - 1] & (SOLID | OUTSIDE))) {
                    printf
                      ("Pass %d e1 Courant violation at (%d %d %d) %f %f %f\n",
                       nps, i, j, k, d1, w1, bathy[j][i]);
                    im = xm[j][i];
                    ip = xp[j][i];
                    ipp = xp[j][ip];
                    imm = xm[j][im];
                    bathy[j][i] =
                      smoothbathy(bathy, i, ip, im, ipp, imm, j, 1, bmin);
                    bathy[j][ip] =
                      smoothbathy(bathy, ip, ipp, i, xp[j][ipp], im, j, 1,
                                  bmin);
                    bathy[j][im] =
                      smoothbathy(bathy, im, i, imm, ip, xm[j][imm], j, 1,
                                  bmin);
                  }
                  if (!(flag[nz - 1][j + 1][i] & (SOLID | OUTSIDE)) &&
                      !(flag[nz - 1][j - 1][i] & (SOLID | OUTSIDE))) {
                    printf
                      ("Pass %d e2 Courant violation at (%d %d %d) %f %f %f\n",
                       nps, i, j, k, d1, w2, bathy[j][i]);
                    jm = ym[j][i];
                    jp = yp[j][i];
                    jpp = yp[jp][i];
                    jmm = ym[jm][i];
                    bathy[j][i] =
                      smoothbathy(bathy, j, jp, jm, jpp, jmm, i, 0, bmin);
                    bathy[jp][i] =
                      smoothbathy(bathy, jp, jpp, j, yp[jpp][i], jm, i, 0,
                                  bmin);
                    bathy[jm][i] =
                      smoothbathy(bathy, jm, j, jmm, jp, ym[jmm][i], i, 0,
                                  bmin);
                  }
                }
              }
            }
          }
        }
      }
      nps++;
      if (nps > 10)
        break;
    }

    d_free_3d(dz);
    d_free_3d(u1flux);
    d_free_3d(u2flux);
    s_free_2d(kb);
  }

  /* Set uniform bathymetry within bathycon cells from the boundary  */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->bathycon) {
      int m;
      /* Loop into the last non-ghost cell bathycon cells into the   */
      /* interior and get the bathymetry at this cell.               */
      for (m = 0; m < open->npts; m++) {	  
	i = open->iloc[m];
	j = open->jloc[m];
	if (open->type & U2BDRY) {
	  jc = j;
	  if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			   (flag[nz - 1][j][i] & OUTSIDE))) {
	    jc--;
	    for (k = 0; k < open->bathycon; k++) {
	      jc = ym[jc][i];
	    }
	    d1 = bathy[jc][i];
	    jc = j - 1;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[jc][i] = d1;
	      jc = ym[jc][i];
	    }
	  } else {
	    for (k = 0; k < open->bathycon; k++)
	      jc = yp[jc][i];
	    d1 = bathy[jc][i];
	    jc = j;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[jc][i] = d1;
	      jc = yp[jc][i];
	    }
	  }
	}
	if (open->type & U1BDRY) {
	  ic = i;
	  if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			   (flag[nz - 1][j][i] & OUTSIDE))) {
	    ic--;
	    for (k = 0; k < open->bathycon; k++)
	      ic = xm[j][ic];
	    d1 = bathy[j][ic];
	    ic = i - 1;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[j][ic] = d1;
	      ic = xm[j][ic];
	    }
	  } else {
	    for (k = 0; k < open->bathycon; k++)
	      ic = xp[j][ic];
	    d1 = bathy[j][ic];
	    ic = i;
	    for (k = 0; k < open->bathycon; k++) {
	      bathy[j][ic] = d1;
	      ic = xp[j][ic];
	    }
	  }
	}
      }
    }
  }

  /* Smooth bathymetry within smooth_z cells from the boundary  */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    if (open->smooth_z && open->smooth_n) {
      double **nbthy = d_alloc_2d(nce1, nce2);
      int m, np;
      /* Loop into the last non-ghost cell bathycon cells into the   */
      /* interior and get the bathymetry at this cell.               */
      for (np = 0; np < open->smooth_n; np++) {
	/* Smooth */
	for (m = 0; m < open->npts; m++) {	  
	  i = open->iloc[m];
	  j = open->jloc[m];
	  if (open->type & U2BDRY) {
	    jc = j;
	    if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      jc--;
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[jc][i] = smooth_bathy(flag, nz, i, jc, xm, xp, ym, yp, bathy, bmin, bmax);
		jc = ym[jc][i];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[jc][i] = smooth_bathy(flag, nz, i, jc, xm, xp, ym, yp, bathy, bmin, bmax);
		jc = yp[jc][i];
	      }
	    }
	  }
	  if (open->type & U1BDRY) {
	    ic = i;
	    if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      ic--;
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[j][ic] = smooth_bathy(flag, nz, ic, j, xm, xp, ym, yp, bathy, bmin, bmax);
		ic = xm[j][ic];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		nbthy[j][ic] = smooth_bathy(flag, nz, ic, j, xm, xp, ym, yp, bathy, bmin, bmax);
		ic = xp[j][ic];
	      }
	    }
	  }
	}
	/* Copy */
	for (m = 0; m < open->npts; m++) {	  
	  i = open->iloc[m];
	  j = open->jloc[m];
	  if (open->type & U2BDRY) {
	    jc = j;
	    if (j == nce2 || (j != nce2 && j == yp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      jc--;
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[jc][i] = nbthy[jc][i];
		jc = ym[jc][i];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[jc][i] = nbthy[jc][i];
		jc = yp[jc][i];
	      }
	    }
	  }
	  if (open->type & U1BDRY) {
	    ic = i;
	    if (i == nce1 || (i != nce1 && i == xp[j][i] &&
			      (flag[nz - 1][j][i] & OUTSIDE))) {
	      ic--;
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[j][ic] = nbthy[j][ic];
		ic = xm[j][ic];
	      }
	    } else {
	      for (k = 0; k < open->smooth_z; k++) {
		bathy[j][ic] = nbthy[j][ic];
		ic = xp[j][ic];
	      }
	    }
	  }
	}
      }
      d_free_2d(nbthy);
    }
  }

  /* Set uniform bathymetry over a zoomed cells if required */
  if(params->dozoom) {
    int ii, jj;
    double d1, d2;
    int **mask;

    mask = i_alloc_2d(nce1, nce2);
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++) {
	mask[j][i] = 0;
	if (bathy[j][i] == LANDCELL || bathy[j][i] == NOTVALID) {
	  mask[j][i] = 1;
	}
      }
    for (n = 1; n <= params->nzoom; n++) {
      i = params->zci[n] - params->zmfe1 / 2;
      j = params->zcj[n] - params->zmfe2 / 2;
      d1 = d2 = 0;
      for (ii = i; ii < i + params->zmfe1; ii++)
	for (jj = j; jj < j + params->zmfe2; jj++)
	  if (bathy[jj][ii] != LANDCELL && bathy[jj][ii] != NOTVALID) {
	    d1 += bathy[jj][ii];
	    d2 += 1.0;
	  }
      d1 /= d2;
      /*d1 = bathy[params->zcj[n]][params->zci[n]];*/
      for (ii = i; ii < i + params->zmfe1; ii++)
	for (jj = j; jj < j + params->zmfe2; jj++) {
	  if (!mask[jj][ii]) {
	    bathy[jj][ii] = d1;
	    mask[jj][ii] = 1;
	  }
	}
      mask[params->zcj[n]][params->zci[n]] = 2;
    }
    if (!(params->dozoom & PRECOND)) {
      /* Set bathymetry at coarse-fine transition boundaries */
      /* Left boundaries */
      for (n = 1; n <= params->nzoom; n++) {
	i = params->zci[n] - params->zmfe1;
	j = params->zcj[n];
	if (i >= 0 && !mask[j][i])
	  mean_z_bathy(i, j, params->zmfe1, params->zmfe2, bathy);
      }
      /* Back boundaries */
      for (n = 1; n <= params->nzoom; n++) {
	i = params->zci[n];
	j = params->zcj[n] - params->zmfe2;
	if (j >= 0 && !mask[j][i])
	  mean_z_bathy(i, j, params->zmfe1, params->zmfe2, bathy);
      }
      /* Right boundaries */
      for (n = 1; n <= params->nzoom; n++) {
	i = params->zci[n] + params->zmfe1;
	j = params->zcj[n];
	if (i < nce1 && !mask[j][i])
	  mean_z_bathy(i, j, params->zmfe1, params->zmfe2, bathy);
      }
      /* Front boundaries */
      for (n = 1; n <= params->nzoom; n++) {
	i = params->zci[n];
	j = params->zcj[n] + params->zmfe2;
	if (j < nce2 && !mask[j][i])
	  mean_z_bathy(i, j, params->zmfe1, params->zmfe2, bathy);
      }
    } else {
      if (params->zmfe1 > 1) {
	ii = 0;
	jj = 1;
      }
      if (params->zmfe2 > 1) {
	ii = 1;
	jj = 0;
      }
      i = params->zci[1];
      j = params->zcj[1];
      ip = i - ii;
      jp = j - jj;
      set_z_bathy(ip, jp, params->zmfe1, params->zmfe2, bathy, bathy[j][i]);
      /*mean_z_bathy(ip, jp, params->zmfe1, params->zmfe2, bathy);*/
      i = params->zci[params->nzoom];
      j = params->zcj[params->nzoom];
      ip = i + ii;
      jp = j + jj;
      set_z_bathy(ip, jp, params->zmfe1, params->zmfe2, bathy, bathy[j][i]);
      /*mean_z_bathy(ip, jp, params->zmfe1, params->zmfe2, bathy);*/
    }
    i_free_2d(mask);
  }

  i_free_2d(xp);
  i_free_2d(xm);
  i_free_2d(yp);
  i_free_2d(ym);
  d_free_2d(h1au1);
  d_free_2d(h2au2);
  d_free_2d(h2au1);
  d_free_2d(h1au2);
  l_free_3d((long ***)flag);

  return (bathy);
}

void mean_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy)
{
  int ii, jj;
  double d1, d2;
  i -= zmfe1 / 2;
  j -= zmfe2 / 2;
  d1 = d2 = 0;
  for (ii = i; ii < i + zmfe1; ii++)
    for (jj = j; jj < j + zmfe2; jj++) {
      if (bathy[jj][ii] != LANDCELL && bathy[jj][ii] != NOTVALID) {
	d1 += bathy[jj][ii];
	d2 += 1.0;
      }
    }
  d1 /= d2;
  for (ii = i; ii < i + zmfe1; ii++)
    for (jj = j; jj < j + zmfe2; jj++)
      if (bathy[jj][ii] != LANDCELL && bathy[jj][ii] != NOTVALID)
	bathy[jj][ii] = d1;
}

void set_z_bathy(int i, int j, int zmfe1, int zmfe2, double **bathy, double val)
{
  int ii, jj;
  i -= zmfe1 / 2;
  j -= zmfe2 / 2;
  for (ii = i; ii < i + zmfe1; ii++)
    for (jj = j; jj < j + zmfe2; jj++) {
      bathy[jj][ii] = val;
    }
}


/* END set_bathy()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths the bathymetry at a point                                 */
/*-------------------------------------------------------------------*/
double smooth_bathy(unsigned long ***flag, int nz, int i, int j, int **xm, int **xp, 
		    int **ym, int **yp, double **bathy, double bmin, double bmax)
{
  double nbthy, kk[4][4];
  int im, ip, jm, jp;

  /* Set the filter weights */
  kk[1][3] = 0.0625;
  kk[2][3] = 0.125;
  kk[3][3] = 0.0625;
  kk[1][2] = 0.125;
  kk[2][2] = 0.25;
  kk[3][2] = 0.125;
  kk[1][1] = 0.0625;
  kk[2][1] = 0.125;
  kk[3][1] = 0.0625;


  if (!(flag[nz - 1][j][i] & (SOLID | OUTSIDE))) {
    double k0 = kk[2][2];
    double k1 = kk[1][2], k2 = kk[3][2];
    double k3 = kk[1][1], k4 = kk[2][1], k5 = kk[3][1];
    double k6 = kk[1][3], k7 = kk[2][3], k8 = kk[3][3];
    if (flag[nz - 1][j][i-1] & (SOLID | OUTSIDE)) {
      k0 += k1;
      k1 = 0.0;
    }
    if (flag[nz - 1][j][i+1] & (SOLID | OUTSIDE)) {
      k0 += k2;
      k2 = 0.0;
    }
    if (flag[nz - 1][j-1][i-1] & (SOLID | OUTSIDE)) {
      k0 += k3;
      k3 = 0.0;
    }
    if (flag[nz - 1][j-1][i] & (SOLID | OUTSIDE)) {
      k0 += k4;
      k4 = 0.0;
    }
    if (flag[nz - 1][j-1][i+1] & (SOLID | OUTSIDE)) {
      k0 += k5;
      k5 = 0.0;
    }
    if (flag[nz - 1][j+1][i-1] & (SOLID | OUTSIDE)) {
      k0 += k6;
      k6 = 0.0;
    }
    if (flag[nz - 1][j+1][i] & (SOLID | OUTSIDE)) {
      k0 += k7;
      k7 = 0.0;
    }
    if (flag[nz - 1][j+1][i+1] & (SOLID | OUTSIDE)) {
      k0 += k8;
      k8 = 0.0;
    }
    im = xm[j][i];
    ip = xp[j][i];
    jm = ym[j][i];
    jp = yp[j][i];
    /* Perform the convolution */
    nbthy = k6 * bathy[jp][im] + 
      k7 * bathy[jp][i] + 
      k8 * bathy[jp][ip] +
      k1 * bathy[j][im] + 
      k0 * bathy[j][i] + 
      k2 * bathy[j][ip] +
      k3 * bathy[jm][im] + 
      k4 * bathy[jm][i] + 
      k5 * bathy[jm][ip];
    nbthy = min(nbthy, -bmin);
    nbthy = max(nbthy, -bmax);
  }
  return (nbthy);
}

/* END smooth_bathy()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to apply a convolution smoothing filter to an element of  */
/* a 2-D array.                                                      */
/*-------------------------------------------------------------------*/
double
smoothbathy(double **a, int c, int p, int m, int pp, int mm, int cc,
            int mode, double bmin)
{
  double kk[5];
  double fi;

  /* Set the filter weights */
  kk[0] = 0.1;
  kk[1] = 0.2;
  kk[2] = 0.4;
  kk[3] = 0.2;
  kk[4] = 0.1;

  /* Perform the convolution */
  if (mode)
    fi = kk[0] * a[cc][mm] + kk[1] * a[cc][m] + kk[2] * a[cc][c] +
      kk[3] * a[cc][p] + kk[4] * a[cc][pp];
  else
    fi = kk[0] * a[mm][cc] + kk[1] * a[m][cc] + kk[2] * a[c][cc] +
      kk[3] * a[p][cc] + kk[4] * a[pp][cc];
  if (fi > -bmin)
    fi = -bmin;
  return (fi);
}

/* END smoothbathy()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Moves an open boundary into the interior                          */
/*-------------------------------------------------------------------*/
void adjust_outside_obc(parameters_t * params, double **bathy)
{
  int n, m, i, j, k, ef, m1, m2;
  int *iloc, *jloc, nb;

  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];

    if (open->bout < 0) continue;

    nb = 0;
    iloc = i_alloc_1d(open->npts);
    jloc = i_alloc_1d(open->npts);

    for (m = 0; m < open->npts; m++) {	  
      i = open->iloc[m];
      j = open->jloc[m];
      
      /* U1 boundaries                                               */
      if (open->type & U1BDRY) {
	ef = 0;
	if (i < params->nce1 && i > 0) {
	  if (!iswet(bathy[j][i])) ef = R_EDGE;
	  if (!iswet(bathy[j][i-1])) ef = L_EDGE;
	}
	if (i == params->nce1) ef = R_EDGE;
	if (i == 0) ef = L_EDGE;
	/* R_EDGE                                                    */
	if (ef == R_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    i = (i > 0) ? i-1 : i;
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;		  
	    open->iloc[m] = i;
	  }
	  m1 = (i > 0) ? i-1 : i;
	  if (iswet(bathy[j][m1])) {
	    /* Use below to enforce no OBCs 1 cell wide              */
	    m2 = (m1 > 0) ? m1-1 : m1;
	    if (!iswet(bathy[j][m2]))
	      bathy[j][m1] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
	/* L_EDGE                                                    */
	if (ef == L_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;
	    i = (i < params->nce1) ? i+1 : i;
	    open->iloc[m] = i;
	  }
	  if (iswet(bathy[j][i])) {
	    /* Use below to enforce no OBCs 1 cell wide              */
	    m1 = (i < params->nce1) ? i+1 : i;
	    if (!iswet(bathy[j][m1]))
	      bathy[j][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i-1] = NOTVALID;
	      nb++;}
	  }
	}
      }
      /* U2 boundaries                                               */
      if (open->type & U2BDRY) {
	ef = 0;
	if (j < params->nce2 && j > 0) {
	  if (!iswet(bathy[j][i])) ef = F_EDGE;
	  if (!iswet(bathy[j-1][i])) ef = B_EDGE;
	}
	if (j == params->nce2) ef = F_EDGE;
	if (j == 0) ef = B_EDGE;
	/* F_EDGE                                                    */
	if (ef == F_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    j = (j > 0) ? j-1 : j;
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;		  
	    open->jloc[m] = j;
	  }
	  m1 = (j > 0) ? j-1 : j;
	  if (iswet(bathy[m1][i])) {
	    m2 = (m1 > 0) ? m1-1 : m1;
	    if (!iswet(bathy[m2][i]))
	      bathy[m1][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
	/* B_EDGE                                                    */
	if (ef == B_EDGE) {
	  for (k = 1; k <= open->bout; k++) {
	    if (iswet(bathy[j][i])) bathy[j][i] = NOTVALID;
	    j = (j < params->nce2) ? j+1 : j;
	    open->jloc[m] = j;
	  }
	  if (iswet(bathy[j][i])) {
	    m1 = (j < params->nce2) ? j+1 : j;
	    if (!iswet(bathy[m1][i]))
	      bathy[j][i] = NOTVALID;
	    else {
	      iloc[nb] = open->iloc[m];
	      jloc[nb] = open->jloc[m];
	      bathy[j-1][i] = NOTVALID;
	      nb++;
	    }
	  }
	}
      }
    }
    for (m = 0; m < nb; m++) {
      open->iloc[m] = iloc[m];
      open->jloc[m] = jloc[m];
      open->npts = nb;
    }
    i_free_1d(iloc);
    i_free_1d(jloc);
  }
  /* Omit the OBC points if the OBC cell is defined as NOTVALID above */
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    int m, dof, nb = 0;
    if (open->bout < 0) continue;
    for (m = 0; m < open->npts; m++) {	  
      i = open->iloc[m];
      j = open->jloc[m];
      dof = 1;
      if (open->type & U1BDRY) {
	if (i > 0 && i < params->nce1 && bathy[j][i] == NOTVALID && 
	    bathy[j][i-1] == NOTVALID) 
	  dof = 0;
	else if (i == 0 && bathy[j][i] == NOTVALID)
	  dof = 0;
	else if (i == params->nce1 && bathy[j][i-1] == NOTVALID)
	  dof = 0;
      }
      if (open->type & U2BDRY) {
	if (j > 0 && j < params->nce2 && bathy[j][i] == NOTVALID && 
	    bathy[j-1][i] == NOTVALID) 
	  dof = 0;
	else if (j == 0 && bathy[j][i] == NOTVALID)
	  dof = 0;
	else if (j == params->nce2 && bathy[j-1][i] == NOTVALID)
	  dof = 0;
      }
      if (dof) {
	open->iloc[nb] = i;
	open->jloc[nb] = j;
	nb++;
      }
    }
    open->npts = nb;
    if (open->bout >= 0)
      hd_warn("Boundary%d (%s) new RANGE (%d,%d)-(%d,%d)\n",
	      n,open->name,open->iloc[0],open->jloc[0],
	      open->iloc[nb-1],open->jloc[nb-1]);
  }
}

/* END adjust_outside_obc()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the sediment layer configuration                            */
/*-------------------------------------------------------------------*/
void read_sed_layers(geometry_t *geom,      /* Global geometry       */
		     parameters_t *params,  /* Input parameters      */
		     master_t *master,      /* Master data           */
		     dump_data_t *dumpdata  /* Dump data             */
		     )
{
  FILE *fp = params->prmfd;
  int i, j, k, m;

  /*-----------------------------------------------------------------*/
  /* Read the sediment layer structure from file                     */
  /* Already read in in read_params()
  if (params->gridz_sed) d_free_1d(params->gridz_sed);
  if (prm_read_darray(fp, "LAYERFACES_SED", &params->gridz_sed,
                      &params->sednz)) {
    if (--params->sednz < 1)
      hd_quit("Number of sediment layers must be 2 or more\n");
  } else {
    double *dz_sed = NULL;
    if (prm_read_int(fp, "NSEDLAYERS", &params->sednz)) {
      dz_sed = d_alloc_1d(params->sednz);
      if (prm_read_darray(fp, "NSEDLAYERS", &dz_sed, &params->sednz)) {
        if (params->sednz) {
          params->gridz_sed = d_alloc_1d(params->sednz + 1);
          params->gridz_sed[params->sednz] = 0.0;
          for (m = params->sednz - 1; m >= 0; m--) {
            params->gridz_sed[m] = params->gridz_sed[m + 1] -
              dz_sed[params->sednz - m - 1];
	  }
        }
      }
      if(dz_sed)
        d_free_1d(dz_sed);
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Transfer to the dumpdata variables                              */
  for (j = 0; j < geom->nce2; j++)
    for (i = 0; i < geom->nce1; i++) {
        for (k = 0; k < params->sednz; k++) {
	  dumpdata->gridz_sed[k][j][i] = params->gridz_sed[k];
	  dumpdata->cellz_sed[k][j][i] = 0.5 * (params->gridz_sed[k] +
						params->gridz_sed[k + 1]);
	}
	dumpdata->gridz_sed[params->sednz][j][i] = params->gridz_sed[params->sednz];
    }

  /*-----------------------------------------------------------------*/
  /* Transfer to the window geometry variables                       */
  for (k = 0; k < params->sednz + 1; k++) {
    c2s_2d(geom, geom->gridz_sed[k], dumpdata->gridz_sed[k], geom->nce1,
	   geom->nce2);
  }
  for (k = 0; k < params->sednz; k++) {
    c2s_2d(geom, geom->cellz_sed[k], dumpdata->cellz_sed[k], geom->nce1,
	   geom->nce2);
  }
}

/* END read_sed_layers()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Removes channels one cell wide from the grid                      */
/*-------------------------------------------------------------------*/
void remove_channel(parameters_t *params, unsigned long ***flag, double **bathy)
{
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int nz = params->nz;
  int i, j, k;
  int ip, im, jp, jm;

  /* Identify one cell channels or canyons */
  for (k = 0; k < nz - 1; k++) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (!(flag[k][j][i] & (SOLID | OUTSIDE))) {
	  ip = (i < nce1 - 1) ? i+1 : i;
	  im = (i > 0) ? i-1 : i;
	  jp = (j < nce2 - 1) ? j+1 : j;
	  jm = (j > 0) ? j-1 : j;

          if (flag[k][j][im] & (SOLID | OUTSIDE) &&
              flag[k][j][ip] & (SOLID | OUTSIDE) &&
              !(flag[k][jm][i] & (SOLID | OUTSIDE)) &&
              !(flag[k][jp][i] & (SOLID | OUTSIDE)))
	    hd_warn("remove_channel: one-grid channel TYPEjj, %d %d %d \n",i, j, k);
          else if (!(flag[k][j][im] & (SOLID | OUTSIDE)) &&
                   !(flag[k][j][ip] & (SOLID | OUTSIDE)) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEii, %d %d %d\n",i, j, k);  
          else if (flag[k][j][im] & (SOLID | OUTSIDE) &&
                   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   !(flag[k][jm][i] & (SOLID | OUTSIDE)) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEjm, %d %d %d\n",i,j, k);
          else if (flag[k][j][im] & (SOLID | OUTSIDE) &&
                   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   !(flag[k][jp][i] & (SOLID | OUTSIDE))) 
	    hd_warn("remove_channel: one-grid channel TYPEjp, %d %d %d\n",i, j, k);
          else if (!(flag[k][j][im] & (SOLID | OUTSIDE)) &&
		   flag[k][j][ip] & (SOLID | OUTSIDE) &&
                   flag[k][jm][i] & (SOLID | OUTSIDE) &&
                   flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEim, %d %d %d\n",i, j, k); 
          else if ( flag[k][j][im] & (SOLID | OUTSIDE) &&
		    !(flag[k][j][ip] & (SOLID | OUTSIDE)) &&
		    flag[k][jm][i] & (SOLID | OUTSIDE) &&
		    flag[k][jp][i] & (SOLID | OUTSIDE))
	    hd_warn("remove_channel: one-grid channel TYPEip, %d %d %d\n",i, j, k);
	}
      }
    }
  }

  /* Fill in channels or canyons one cell wide at the surface */
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (!(flag[nz-1][j][i] & (SOLID | OUTSIDE))) {
	ip = (i < nce1 - 1) ? i+1 : i;
	im = (i > 0) ? i-1 : i;
	jp = (j < nce2 - 1) ? j+1 : j;
	jm = (j > 0) ? j-1 : j;
	
	if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
	    flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
	    !(flag[nz-1][jm][i] & (SOLID | OUTSIDE)) &&
	    !(flag[nz-1][jp][i] & (SOLID | OUTSIDE))) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (!(flag[nz-1][j][im] & (SOLID | OUTSIDE)) &&
		 !(flag[nz-1][j][ip] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 !(flag[nz-1][jm][i] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 !(flag[nz-1][jp][i] & (SOLID | OUTSIDE))) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if (!(flag[nz-1][j][im] & (SOLID | OUTSIDE)) &&
		 flag[nz-1][j][ip] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		 flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
	else if ( flag[nz-1][j][im] & (SOLID | OUTSIDE) &&
		  !(flag[nz-1][j][ip] & (SOLID | OUTSIDE)) &&
		  flag[nz-1][jm][i] & (SOLID | OUTSIDE) &&
		  flag[nz-1][jp][i] & (SOLID | OUTSIDE)) {
	  bathy[j][i] = LANDCELL;
	  for (k = 0; k < nz; k++) {
	    flag[k][j][i] |= SOLID;
	  }
	}
      }
    }
  }
}

/* END remove_channel()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Cascade fills a 2D array over new bathymetry                      */
/*-------------------------------------------------------------------*/
void data_infill_2d(master_t *master, double **a)
{
  int i, j, ci, cj;
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  if(master->data_infill) {
    parameters_t *params = master->params;
    double **bathy = params->topo;
    for (j = 0; j < geom->nce2; j++)
      for (i = 0; i < geom->nce1; i++) {
	if (bathy[j][i] < dumpdata->gridz[0] || bathy[j][i] >= master->etamax)
	  if (find_closest_nonnan(dumpdata, a, i, j, geom->nz-1, &ci, &cj))
	    a[j][i] = a[cj][ci];
      }
  }
}

/* END data_infill_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Cascade fills a 2D array over new bathymetry                      */
/*-------------------------------------------------------------------*/
void data_infill_3d(master_t *master, double ***a)
{
  int i, j, k, ci, cj;
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  /*double **bathy = dumpdata->botz;*/

  if(master->data_infill) {
    parameters_t *params = master->params;
    double **bathy = params->topo;
    if (master->i1 == NULL && master->i2 == NULL) {
      master->i1 = i_alloc_3d(geom->nce1, geom->nce2, geom->nz+1);
      master->i2 = i_alloc_3d(geom->nce1, geom->nce2, geom->nz+1);
      for (k = 0; k < geom->nz; k++)
	for (j = 0; j < geom->nce2; j++)
	  for (i = 0; i < geom->nce1; i++) {
	    master->i1[k][j][i] = master->i2[k][j][i] = -1;
	    if (bathy[j][i] < dumpdata->gridz[0] || bathy[j][i] >= master->etamax ||
		bathy[j][i] >= dumpdata->gridz[k + 1])
	      if (find_closest_nonnan(dumpdata, a[k], i, j, k, &ci, &cj)) {
		master->i1[k][j][i] = ci;
		master->i2[k][j][i] = cj;
		a[k][j][i] = a[k][cj][ci];
	      }
	  }
    } else {
      for (k = 0; k < geom->nz; k++)
	for (j = 0; j < geom->nce2; j++)
	  for (i = 0; i < geom->nce1; i++) {
	    if (bathy[j][i] < dumpdata->gridz[0] || bathy[j][i] >= master->etamax ||
		bathy[j][i] >= dumpdata->gridz[k + 1])
	      ci = master->i1[k][j][i];
	      cj = master->i2[k][j][i];
	      if (ci >= 0 && cj >= 0)
		a[k][j][i] = a[k][cj][ci];
	  }
    }
  }
}

/* END data_infill_3d()                                              */
/*-------------------------------------------------------------------*/


int iswet(int bathy)
{
  if (bathy == NOTVALID || bathy == LANDCELL)
    return 0;
  return 1;
}
