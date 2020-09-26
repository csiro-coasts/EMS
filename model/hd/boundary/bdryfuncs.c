/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/boundary/bdryfuncs.c
 *  
 *  Description:
 *  Standard boundary functions.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: bdryfuncs.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include "eqn_parser.h"
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Custom data associated with the uav custom boundaries             */
typedef struct {
  int *nmap;
  int *tmap;
} uav_data_t;

typedef struct {
  int ntsfiles;
  timeseries_t **tsfiles;
  cstring *filenames;
} tsfiles_t;

void getuv_from_hdstd(double t, double x, double y, double z,
                      int nts, timeseries_t **ts, cstring * filenames,
                      double *u, double *v);
void getuv(double t, double x, double y, double z, int nts,
           timeseries_t **ts, cstring * filenames, double *u, double *v,
	   char *uname, char * vname);
tsfiles_t *tsfiles_alloc(master_t *master, bdry_details_t *data);
void read_bdry_zone_std(master_t *master, open_bdrys_t *open, int cc, int mode);

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Free data for the custom routines                                */
/*------------------------------------------------------------------*/
/* Multifile time series                                            */
/*------------------------------------------------------------------*/
void bf_ts_free(master_t *master, bdry_details_t *data)
{
  int i;
  tsfiles_t *tsf;

  tsf = data->custdata;

  if (tsf->ntsfiles) {
    free(tsf->filenames);
    for (i = 0; i < tsf->ntsfiles; ++i)
      hd_ts_free(master, tsf->tsfiles[i]);    
  }
}

/* Maps for uv_to_u1av, uv_to_u2av routines */
void bf_c2cc_free(geometry_t *window, bdry_details_t *data)
{
  uav_data_t *c2cc = data->custdata;
  if (c2cc->nmap)
    i_free_1d(c2cc->nmap);
  if (c2cc->tmap)
    i_free_1d(c2cc->tmap);
  free((uav_data_t *)c2cc);
}

/* Null routine */
void bf_void_free(geometry_t *window, bdry_details_t *data)
{
}

/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert a HD standard dump file velocity vector into U1          */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u1_init_m(master_t *master, open_bdrys_t *open,
                           bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;
  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames,
                          "u1", schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames,
                          "u2", schedule->start_time, schedule->stop_time);
  }
}

void bf_hdstd_to_u1_init_w(geometry_t *window, open_bdrys_t *open,
                           bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_hdstd_to_u1_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  int zone = open->relax_zone_nor;
  int mode = U1BDRY|U1GEN;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
    zone = open->relax_zone_tan;
    mode = U2BDRY|U1GEN;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t, geom->u1x[cs], geom->u1y[cs],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u1[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
  /* If a block of data is read, then read the rest */
  if (zone) {
    for (cc = sb; cc <= eb; cc++) {
      c = open->obc_e1[cc];
      read_bdry_zone_std(master, open, cc, mode);
    }
  }
}

double bf_hdstd_to_u1_w(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, open_bdrys_t *open, double t,
                        int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u1[cc]);
}

void bf_hdstd_to_u1_t(master_t *master, open_bdrys_t *open_w,
                      bdry_details_t *data, geometry_t *window,
                      window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1; 
  int eb = open_w->no3_e1;
  int ebm = open_m->no3_e1;
  int zone = open_w->relax_zone_nor;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no3_e1 + 1;
    eb = open_w->to3_e1;
    ebm = open_m->to3_e1;
    zone = open_w->relax_zone_nor;
  }
  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u1[cc];
    open_w->transfer_u1[cc] = open_m->transfer_u1[c];
  }
  if (zone) {
    int zn;
    for (cc = sb; cc <= eb; cc++) {
      c = open_w->tmap_u1[cc];
      for (zn = 1; zn < zone; zn++) {
	open_w->transfer_u1[cc+eb] = open_m->transfer_u1[c+ebm];
      }
    }
  }
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert a HD standard dump file velocity vector into U2          */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u2_init_m(master_t *master, open_bdrys_t *open,
                           bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;
  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames,
                          "u1", schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames,
                          "u2", schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_hdstd_to_u2_init_w(geometry_t *window, open_bdrys_t *open,
                           bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_hdstd_to_u2_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  int zone = open->relax_zone_nor;
  int mode = U2BDRY|U2GEN;

  if(data->type & TAN) {
    sb = open->no3_e2 + 1;
    eb = open->to3_e2;
    zone = open->relax_zone_tan;
    mode = U1BDRY|U2GEN;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    cs = geom->m2d[c];
    thetau2 = geom->thetau2[cs];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t, geom->u2x[cs], geom->u2y[cs],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u2[cc] = cos(thetau2) * v - sin(thetau2) * u;

    if (zone)
      read_bdry_zone_std(master, open, cc, mode);
  }
}

double bf_hdstd_to_u2_w(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, open_bdrys_t *open, double t,
                        int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u2[cc]);
}

void bf_hdstd_to_u2_t(master_t *master, open_bdrys_t *open_w,
                      bdry_details_t *data, geometry_t *window,
                      window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1; 
  int eb = open_w->no3_e2;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no3_e2 + 1;
    eb = open_w->to3_e2;
  }
  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u2[cc];
    open_w->transfer_u2[cc] = open_m->transfer_u2[c];
  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert a HD standard dump file U1/U2 velocity vector into U1AV  */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u1av_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t3d, geom->u1x[cs], geom->u1y[cs],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u1av[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert a HD standard dump file U1/U2 velocity vector into U2AV  */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u2av_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e2 + 1;
    eb = open->to3_e2;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    cs = geom->m2d[c];
    thetau2 = geom->thetau2[cs];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t3d, geom->u2x[cs], geom->u2y[cs],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u2av[cc] = cos(thetau2) * v - sin(thetau2) * u;
  }
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert U/V components of velocity into a U1                     */
/*------------------------------------------------------------------*/
void bf_uv_to_u1_init_m(master_t *master, open_bdrys_t *open,
                        bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;

  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "u",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "v",
                          schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_uv_to_u1_init_w(geometry_t *window, open_bdrys_t *open,
                        bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_uv_to_u1_m(geometry_t *geom, master_t *master,
                   open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;
  double z, eps = 1e-10;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0, z = geom->cellz[c] * master->Ds[cs];
    z = (z) ? z : eps; /* If z = 0.0, then the weighted interpolation returns strange values */
    getuv(master->t, geom->u1x[cs], geom->u1y[cs], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

/*
 * Performs the same function as bf_uv_to_u1_m for MPI
 */
void bf_uv_to_u1_mw(geometry_t *window, master_t *master, 
		    open_bdrys_t *open_w, bdry_details_t *data)
{
  geometry_t   *geom   = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  win_priv_t *wincon = window->wincon;
  window_t   *windat = window->windat;
  int c, cc, cs;
  int sb = 1; 
  int eb = open_w->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;
  double z, eps = 1e-10;

  if(data->type & TAN) {
    sb = open_w->no3_e1 + 1;
    eb = open_w->to3_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open_w->obc_e1[cc];
    cs = window->m2d[c];
    thetau1 = window->thetau1[cs];
    u = 0.0, v = 0.0, z = window->cellz[c] * wincon->Ds[cs];
    z = (z) ? z : eps; /* If z = 0.0, then the weighted interpolation returns strange values */
    getuv(windat->t, window->u1x[cs], window->u1y[cs], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open_m->transfer_u1[open_w->tmap_u1[cc]] = cos(thetau1) * u + sin(thetau1) * v;
  }
}


double bf_uv_to_u1_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u1[cc]);
}

void bf_uv_to_u1_t(master_t *master, open_bdrys_t *open_w,
                   bdry_details_t *data, geometry_t *window,
                   window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1; 
  int eb = open_w->no3_e1;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no3_e1 + 1;
    eb = open_w->to3_e1;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u1[cc];
    open_w->transfer_u1[cc] = open_m->transfer_u1[c];
  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert U/V components of velocity into a U2                     */
/*------------------------------------------------------------------*/
void bf_uv_to_u2_init_m(master_t *master, open_bdrys_t *open,
                        bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;
  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "u",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "v",
                          schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_uv_to_u2_init_w(geometry_t *window, open_bdrys_t *open,
                        bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_uv_to_u2_m(geometry_t *geom, master_t *master,
                   open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1;
  int eb = open->no3_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  double z, eps = 1e-10;

  if(data->type & TAN) {
    sb = open->no3_e2 + 1;
    eb = open->to3_e2;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    cs = geom->m2d[c];
    thetau2 = geom->thetau2[cs];
    u = 0.0, v = 0.0, z = geom->cellz[c] * master->Ds[cs];
    z = (z) ? z : eps;
    getuv(master->t, geom->u2x[cs], geom->u2y[cs], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u2[cc] = cos(thetau2) * v - sin(thetau2) * u;
  }
}

/*
 * Performs the same function as bf_uv_to_u1_m for MPI
 */
void bf_uv_to_u2_mw(geometry_t *window, master_t *master, 
		    open_bdrys_t *open_w, bdry_details_t *data)
{
  geometry_t   *geom   = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  win_priv_t *wincon = window->wincon;
  window_t   *windat = window->windat;
  int c, cc, cs;
  int sb = 1;
  int eb = open_w->no3_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  double z, eps = 1e-10;

  if(data->type & TAN) {
    sb = open_w->no3_e2 + 1;
    eb = open_w->to3_e2;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open_w->obc_e2[cc];
    cs = window->m2d[c];
    thetau2 = window->thetau2[cs];
    u = 0.0, v = 0.0, z = window->cellz[c] * wincon->Ds[cs];
    z = (z) ? z : eps;
    getuv(windat->t, geom->u2x[cs], geom->u2y[cs], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open_m->transfer_u2[open_w->tmap_u2[cc]] = cos(thetau2) * v - sin(thetau2) * u;
  }
}

double bf_uv_to_u2_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u2[cc]);
}

void bf_uv_to_u2_t(master_t *master, open_bdrys_t *open_w,
                   bdry_details_t *data, geometry_t *window,
                   window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1;
  int eb = open_w->no3_e2;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no3_e2 + 1;
    eb = open_w->to3_e2;
  }
  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u2[cc];
    open_w->transfer_u2[cc] = open_m->transfer_u2[c];
  }
}


void getuv_from_hdstd(double t, double x, double y, double z,
                      int nts, timeseries_t **ts, cstring * filenames,
                      double *u, double *v)
{

  /* Read the u1 and u2 from the 'existing' model grid. */
  double u1val = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames,
                                          "u1", t, x, y, z);
  double u2val = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames,
                                          "u2", t, x, y, z);
  double thetau1 = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames,
                                            "thetau1", t, x, y, z);
  double thetau2 = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames,
                                            "thetau2", t, x, y, z);
  double sinth = (sin(thetau1) + sin(thetau2)) / 2;
  double costh = (cos(thetau1) + cos(thetau2)) / 2;
  *u = u1val * costh - u2val * sinth;
  *v = u1val * sinth + u2val * costh;
}


void getuv(double t, double x, double y, double z,
           int nts, timeseries_t **ts, cstring * filenames, double *u,
           double *v, char *uname, char * vname)
{
  *u = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames, uname, 
					t, x, y, z);
  *v = hd_ts_multifile_eval_xyz_by_name(nts, ts, filenames, vname, 
					t, x, y, z);
}


tsfiles_t *tsfiles_alloc(master_t *master, bdry_details_t *data)
{
  tsfiles_t *tsf = (tsfiles_t *)malloc(sizeof(tsfiles_t));
  int i;
  char buf[MAXSTRLEN];

  for (i = 0; i < data->nargs; ++i) {
    /* don't double append on boundary copy */
    char *bp = NULL;
    if ( (bp = strstr(data->args[i], master->bdrypath)) == NULL)
      if (strlen(master->bdrypath)) {
	sprintf(buf, "%s%s", master->bdrypath, data->args[i]);
	strcpy(data->args[i], buf);
      }
  }

  tsf->ntsfiles = data->nargs;
  tsf->tsfiles = hd_ts_multifile_read(master, data->nargs, data->args);
  tsf->filenames = (cstring *) malloc(sizeof(cstring) * tsf->ntsfiles);
  for (i = 0; i < tsf->ntsfiles; ++i) {
    /* strcpy(tsf->filenames[i], ((char**)data->args)[i]); */
    strcpy(tsf->filenames[i], data->args[i]);
  }
  return tsf;
}

/*------------------------------------------------------------------*/
/* Custom routines to convert (u,v) data into u1av or u2av.         */
/* The eastward (u) and northward (v) 3D currents are rotated onto  */
/* the grid and depth averaged.                                     */
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/* Window initialisation. Make a c2cc map and store in the custom   */
/* pointer in the custom data structure; (int *)custdata.           */
/* This initialisation can be used for U1 and U2 boundaries.        */
/*------------------------------------------------------------------*/
void bf_uv_to_uav_init_w(geometry_t *window, open_bdrys_t *open,
			 bdry_details_t *data, bdry_details_t *data_in)
{
  int c, cc;
  uav_data_t *c2cc = (uav_data_t *)malloc(sizeof(uav_data_t));

  /* Initialise */
  memset(c2cc, 0, sizeof(uav_data_t));
  c2cc->nmap = i_alloc_1d(window->enon + 1);
  c2cc->tmap = i_alloc_1d(window->enon + 1);

  if (open->type & U1BDRY) {
    for (cc = 1; cc <= open->no3_e1; cc++) {
      c = open->obc_e1[cc];
      c2cc->nmap[c] = cc;
    }
    for (cc = open->no3_e2 + 1; cc <= open->to3_e2; cc++) {
      c = open->obc_e2[cc];
      c2cc->tmap[c] = cc;
    }
  } else if (open->type & U2BDRY) {
    for (cc = 1; cc <= open->no3_e2; cc++) {
      c = open->obc_e2[cc];
      c2cc->nmap[c] = cc;
    }
    for (cc = open->no3_e1 + 1; cc <= open->to3_e1; cc++) {
      c = open->obc_e1[cc];
      c2cc->tmap[c] = cc;
    }
  }
  data->custdata = (void *)c2cc;
}

/*------------------------------------------------------------------*/
/* Read in rotated (u,v) 3D currents to transfer_u1av               */
/*------------------------------------------------------------------*/
void bf_uv_to_u1av_m(geometry_t *geom, master_t *master,
		     open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[cs], geom->u1y[cs], geom->cellz[c],
	  tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1av[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

/*------------------------------------------------------------------*/
/* Vertically integrate rotated (u,v) 3D currents to a normal u1av  */
/*------------------------------------------------------------------*/
double bf_uv_to_u1av_w(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, open_bdrys_t *open, 
		       double t, int c, int cc, bdry_details_t *data)
{
  int zm1, c1;                         /* Sparse location / counter */
  int sb = 1;
  int eb = open->no3_e1;
  double depth = 0.0;
  double uav = 0.0;
  int *c2cc;
  uav_data_t *udata = (uav_data_t *)data->custdata;
  /*double *eta = (windat->eta_rlx) ? windat->eta_rlx : windat->eta;*/

  c2cc = udata->nmap;
  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
    c2cc = udata->tmap;
  }

  /* Vertically integrate the velocity */
  zm1 = window->zm1[c];
  while(c != zm1) {
    c1 = c2cc[c];
    uav += open->transfer_u1av[c1] * windat->dzu1[c];
    depth += windat->dzu1[c];
    c = zm1;
    zm1 = window->zm1[zm1];
  }
  /* Divide by total depth and reset the transfer vector */
  open->transfer_u1av[cc] = uav / depth;

  return (open->transfer_u1av[cc]);
}

/*------------------------------------------------------------------*/
/* Read in rotated (u,v) 3D currents to transfer_u2av               */
/*------------------------------------------------------------------*/
void bf_uv_to_u2av_m(geometry_t *geom, master_t *master,
		     open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1;
  int eb = open->no3_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e2 + 1;
    eb = open->to3_e2;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    cs = geom->m2d[c];
    thetau2 = geom->thetau2[cs];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u2x[cs], geom->u2y[cs], geom->cellz[c],
	  tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u2av[cc] = cos(thetau2) * v - sin(thetau2) * u;
  }
}

/*------------------------------------------------------------------*/
/* Vertically integrate rotated (u,v) 3D currents to a normal u2av  */
/*------------------------------------------------------------------*/
double bf_uv_to_u2av_w(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, open_bdrys_t *open, 
		       double t, int c, int cc, bdry_details_t *data)
{
  int zm1, c1;                         /* Sparse location / counter */
  int sb = 1;
  int eb = open->no3_e2;
  double depth = 0.0;
  double vav = 0.0;
  int *c2cc;
  uav_data_t *udata = (uav_data_t *)data->custdata;

  c2cc = udata->nmap;
  if(data->type & TAN) {
    sb = open->no3_e2 + 1;
    eb = open->to3_e2;
    c2cc = udata->tmap;
  }

  /* Vertically integrate the velocity */
  zm1 = window->zm1[c];
  while(c != zm1) {
    c1 = c2cc[c];
    vav += open->transfer_u2av[c1] * windat->dzu2[c];
    depth += windat->dzu2[c];
    c = zm1;
    zm1 = window->zm1[zm1];
  }

  /* Divide by total depth and reset the transfer vector */
  open->transfer_u2av[cc] = vav / depth;

  return (open->transfer_u2av[cc]);
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert UAV/VAV components of velocity into a U1AV               */
/*------------------------------------------------------------------*/
void bf_uvav_to_u1av_init_m(master_t *master, open_bdrys_t *open,
			    bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;

  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "uav",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "vav",
                          schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_uvav_to_u1av_init_w(geometry_t *window, open_bdrys_t *open,
			    bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_uvav_to_u1av_m(geometry_t *geom, master_t *master,
		       open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no2_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;

  if(data->type & TAN) {
    sb = open->no2_e1 + 1;
    eb = open->to2_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[cs], geom->u1y[cs], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u1av[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}


double bf_uvav_to_u1av_w(geometry_t *window, window_t *windat,
			 win_priv_t *wincon, open_bdrys_t *open, double t,
			 int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u1av[cc]);
}

void bf_uvav_to_u1av_t(master_t *master, open_bdrys_t *open_w,
		       bdry_details_t *data, geometry_t *window,
		       window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1; 
  int eb = open_w->no2_e1;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no2_e1 + 1;
    eb = open_w->to2_e1;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u1[cc];
    open_w->transfer_u1av[cc] = open_m->transfer_u1av[c];
  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert UAV/VAV components of velocity into a U2AV               */
/*------------------------------------------------------------------*/
void bf_uvav_to_u2av_init_m(master_t *master, open_bdrys_t *open,
			    bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;
  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "uav",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "vav",
                          schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_uvav_to_u2av_init_w(geometry_t *window, open_bdrys_t *open,
			    bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_uvav_to_u2av_m(geometry_t *geom, master_t *master,
		       open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1;
  int eb = open->no2_e2;
  double thetau2, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no2_e2 + 1;
    eb = open->to2_e2;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    cs = geom->m2d[c];
    thetau2 = geom->thetau2[cs];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u2x[cs], geom->u2y[cs], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u2av[cc] = cos(thetau2) * v - sin(thetau2) * u;
  }
}

double bf_uvav_to_u2av_w(geometry_t *window, window_t *windat,
			 win_priv_t *wincon, open_bdrys_t *open, double t,
			 int c, int cc, bdry_details_t *data)
{
  return (open->transfer_u2av[cc]);
}

void bf_uvav_to_u2av_t(master_t *master, open_bdrys_t *open_w,
		       bdry_details_t *data, geometry_t *window,
		       window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1;
  int eb = open_w->no2_e2;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no2_e2 + 1;
    eb = open_w->to2_e2;
  }
  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u2[cc];
    open_w->transfer_u2av[cc] = open_m->transfer_u2av[c];
  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert U/V components of velocity into a U1 and adjust with     */
/* data read from a UAV/VAV to U1AV computation.                    */
/*------------------------------------------------------------------*/
void bf_uv_adj_u2_init_m(master_t *master, open_bdrys_t *open,
			 bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;

  if (schedule != NULL) { 
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "u",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "v",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "uav",
                          schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "vav",
                          schedule->start_time, schedule->stop_time);
  }
}

/* Convert a HD standard dump file velocity vector into U1 */
void bf_uv_adj_u1_init_w(geometry_t *window, open_bdrys_t *open,
			 bdry_details_t *data, bdry_details_t *data_in)
{
}

void bf_uv_adj_u1_m(geometry_t *geom, master_t *master,
		    open_bdrys_t *open, bdry_details_t *data)
{
  int c, cc, cs;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;

  /* Read the 3D velocities */
  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv(master->t, geom->u1x[cs], geom->u1y[cs], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
  /* Read the 2D velocities */
  if(data->type & TAN) {
    sb = open->no2_e1 + 1;
    eb = open->to2_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv(master->t, geom->u1x[cs], geom->u1y[cs], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u1av[cc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

double bf_uv_adj_u1_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data)
{

  /* Vertically integrate the 3D velocity */
  bf_uv_to_u2av_w(window, windat, wincon, open, t, c, cc, data);

  /* Adjust the 3D velocity to 2D velocity */

  return (open->transfer_u1[cc]);
}

void bf_uv_adj_u1_t(master_t *master, open_bdrys_t *open_w,
		    bdry_details_t *data, geometry_t *window,
		    window_t *windat)
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  int c, cc;                    /* Sparse location / counter */
  int sb = 1; 
  int eb = open_w->no3_e1;

  /* Master and slave transfer vectors point to each other for 1 win */
  if (master->nwindows == 1)
    return;

  /* Transfer data read from file from the master to the slave */
  if(data->type & TAN) {
    sb = open_w->no3_e1 + 1;
    eb = open_w->to3_e1;
  }

  for (cc = sb; cc <= eb; cc++) {
    c = open_w->tmap_u1[cc];
    open_w->transfer_u1[cc] = open_m->transfer_u1[c];
  }
}

/*
 * Local helper structs
 */
typedef struct {
  int cc;
  open_bdrys_t *open;
} bf_use_eqn_helper;


/*
 * Helper function for the use_eqn to find the tracer names
 */
static double *bf_use_eqn_find_tracer(const char *str, void *data)
{
  bf_use_eqn_helper *h = (bf_use_eqn_helper *)data;

  int tm;
  double *valPtr = NULL;

  /*
   * Call out to the library function
   */
  tm = tracer_find_index(str, master->ntr, master->trinfo_3d);
  
  /* See if we found anything */
  if (tm > -1) {
    int c  = h->open->obc_t[h->cc];
    valPtr = &master->tr_wc[tm][c];
  } else {
    tm = tracer_find_index(str, master->ntrS, master->trinfo_2d);
    if (tm > -1) {
      int c  = h->open->obc_t[h->cc];
      c = master->geom->m2d[c];
      valPtr = &master->tr_wcS[tm][c];
    }
    if (strcmp(str,"eta") == 0) {
      int c  = h->open->obc_t[h->cc];
      c = master->geom->m2d[c];
      valPtr = &master->eta[c];
    }
  }

  return(valPtr);
}

/*
 * Functions that use an equation to calculate the boundary value
 *
 * Init function is called once
 */
void bf_use_eqn_init_m(master_t       *master, 
		       open_bdrys_t   *open, 
		       bdry_details_t *data)
{
  int cc;
  int num = open->no3_t;
  void **pdata = p_alloc_1d(num); // xxx gotta free this somewhere
  bf_use_eqn_helper h;
  char err[MAXSTRLEN];

  sprintf(err, "boundary %d, tracer %s", open->id, data->name);

  h.open = open;
  /* loop through and initialise the parser for each boundary cell */
  for (cc = 1; cc <= num; cc++) {
    h.cc = cc;
    pdata[cc-1] = EqnCreateParser(data->args[0], bf_use_eqn_find_tracer, 
				  (void *)&h, err);
  }

  /* Set custom data */
  data->custdata = pdata;
}

/*
 * The actual evaluate function is called every boundary dt
 */
void bf_use_eqn_m(geometry_t     *geom,
		  master_t       *master,
		  open_bdrys_t   *open,
		  bdry_details_t *data)
{
  int cc, tm, tn; /* counters */
  /* 
   * These are the parser object pointers 
   */
  void **ptr = (void **)data->custdata;

  /*
   * Find this index
   * Note: we could cache this
   */
  tn = tracer_find_index(data->name, master->ntr, master->trinfo_3d);

  /* Fill the boundary */
  for (cc=1; cc<open->no3_t; cc++) {
    tm = open->trm[tn];
    open->t_transfer[tm][cc] = EqnGetValue(ptr[cc-1]);
  }
}

/*
 * We're using the master directly so there is no need for any
 * transfers, however we still need this function as is called
 */
void bf_use_eqn_trans(master_t *master, open_bdrys_t *open_w,
		      bdry_details_t *data, geometry_t *window,
		      window_t *windat)
{
  /* no op */
}


double bf_use_eqn_w(geometry_t *window, window_t *windat,
		    win_priv_t *wincon, open_bdrys_t *open, double t,
		    int c, int cc, bdry_details_t *data)
{
  /*
   * Find this index
   * Note: we could cache this
   */
  int tn = tracer_find_index(data->name, wincon->ntr, wincon->trinfo_3d);
  int tm = open->trm[tn];
  int oc = 1; // xxx bad value
  for (oc=1; oc <= open->no3_t; oc++) {
    if (c == open->obc_t[oc])
      break;
  }

  return(open->t_transfer[tm][oc]);
}

void bf_eqn_free(master_t *master, bdry_details_t *data)
{
  free((void **)data->custdata);
}


/*-------------------------------------------------------------------*/
/* Reads in 3D blocks of data for flow relaxation and nudging        */
/*-------------------------------------------------------------------*/
void read_bdry_zone_std(master_t *master, open_bdrys_t *open, int cc, int mode)
{
  geometry_t *geom = master->geom;
  bdry_details_t *data;
  tsfiles_t *tsf;
  double *tvec, thetau1, u, v;
  int *obc;
  int tinc;
  int c, cs, i;
  int zone;

  if (mode & U1BDRY && mode & U1GEN) {
    data = &open->datau1;
    tvec = open->transfer_u1;
    tinc = open->no3_e1;
    zone = open->relax_zone_nor;    
    obc = open->obc_e1;
  } else if (mode & U1BDRY && mode & U2GEN) {
    data = &open->datau2;
    tvec = open->transfer_u2;
    tinc = open->to3_e2;
    obc = open->obc_e2;
    zone = open->relax_zone_tan;
  } if (mode & U2BDRY && mode & U2GEN) {
    data = &open->datau2;
    tvec = open->transfer_u2;
    tinc = open->no3_e2;
    zone = open->relax_zone_nor;
    obc = open->obc_e2;
  } else if (mode & U2BDRY && mode & U1GEN) {
    data = &open->datau1;
    tvec = open->transfer_u1;
    tinc = open->to3_e1;
    zone = open->relax_zone_tan;
    obc = open->obc_e1;
  }
  tsf = (tsfiles_t *)data->custdata;

  c = obc[cc];
  for (i = 1; i < zone; i++) {
    c = open->nmap[c];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[cs];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t, geom->u1x[cs], geom->u1y[cs],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    tvec[cc + i * tinc] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

/* END read_bdry_zone_std()                                          */
/*-------------------------------------------------------------------*/

/* EOF */
