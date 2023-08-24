/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/boundary/bdryfuncs.c
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
 *  $Id: bdryfuncs.c 7151 2022-07-07 02:30:05Z her127 $
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
void build_gauss_ref(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     double dist, double gscale, double tstart);

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
  int ee, e, es, c;
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
    mode = U1BDRY|U1GEN;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t, geom->u1x[es], geom->u1y[es],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u1[ee] = cos(thetau1) * u + sin(thetau1) * v;
  }
  /* If a block of data is read, then read the rest */
  if (zone) {
    for (ee = sb; ee <= eb; ee++) {
      e = open->obc_e1[ee];
      read_bdry_zone_std(master, open, ee, mode);
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
  int ee, e, es, c;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  int zone = open->relax_zone_nor;
  int mode = U2BDRY|U2GEN;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
    zone = open->relax_zone_tan;
    mode = U1BDRY|U2GEN;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t, geom->u1x[es], geom->u1y[es],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u2[ee] = cos(thetau1) * u + sin(thetau1) * v;

    if (zone)
      read_bdry_zone_std(master, open, ee, mode);
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
  int ee, e, es, c;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t3d, geom->u1x[es], geom->u1y[es],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u1av[ee] = cos(thetau1) * u + sin(thetau1) * v;
  }
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert a HD standard dump file U1/U2 velocity vector into U2AV  */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u2av_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data)
{
  int ee, e, es, c;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv_from_hdstd(master->t3d, geom->u1x[es], geom->u1y[es],
                     geom->cellz[c], tsf->ntsfiles, tsf->tsfiles,
                     tsf->filenames, &u, &v);
    open->transfer_u2av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int ee, e, es, c, cs;
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
  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0, z = geom->cellz[c] * master->Ds[cs];
    z = (z) ? z : eps; /* If z = 0.0, then the weighted interpolation returns strange values */
    getuv(master->t, geom->u1x[es], geom->u1y[es], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int e, ee, es, c, cs;
  int sb = 1;
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  double z, eps = 1e-10;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    cs = geom->m2d[c];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0, z = geom->cellz[c] * master->Ds[cs];
    z = (z) ? z : eps;
    getuv(master->t, geom->u1x[es], geom->u1y[es], z,
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u2[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
    c = open_w->tmap_u2[cc];
    open_w->transfer_u2[cc] = open_m->transfer_u2[c];
  }
}


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Convert UGRID OBC velocity data into a U1                        */
/*------------------------------------------------------------------*/
void bf_ug_to_u1_init_m(master_t *master, open_bdrys_t *open,
                        bdry_details_t *data)
{
  tsfiles_t *tsf = tsfiles_alloc(master, data);

  data->custdata = tsf;

  if (schedule != NULL) {
    hd_ts_multifile_check(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "u1",
                          schedule->start_time, schedule->stop_time);
  }
}

/*------------------------------------------------------------------*/
/* UGRID to normal velocity. Normal transfer vectors are filled     */
/* first, so read the normal and tangential velocities as a block   */
/* into the u1 transfer vector. The tangential velocities are       */
/* copied to the transfer_u2 vector if bf_ug_to_u2_m() is called.   */
/*------------------------------------------------------------------*/
void bf_ug_to_u1_m(geometry_t *geom, master_t *master,
                   open_bdrys_t *open, bdry_details_t *data)
{
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;
  int ee;

  hd_ts_multifile_eval_isparse(tsf->ntsfiles, tsf->tsfiles, tsf->filenames, "u1",
			      master->d3, master->t, open->to3_e1, master->thIO);
  for (ee = 1; ee <= open->to3_e1; ee++)
    open->transfer_u1[ee] = master->d3[ee - 1];
}

/*------------------------------------------------------------------*/
/* UGRID to tangential velocity. Copy to the transfer_u2 vector     */
/* from the transfer_u1 vector (previously read as a data block).   */
/*------------------------------------------------------------------*/
void bf_ug_to_u2_m(geometry_t *geom, master_t *master,
                   open_bdrys_t *open, bdry_details_t *data)
{
  int ee;

  for (ee = open->no3_e1 + 1; ee <= open->to3_e1; ee++) {
    open->transfer_u2[ee] = open->transfer_u1[ee];
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
    if (strlen(master->bdrypath)) {
      sprintf(buf, "%s%s", master->bdrypath, data->args[i]);
      strcpy(data->args[i], buf);
    }
  }
  if (endswith(data->args[0], ".mpk"))
    hd_warn("Attempting to read file %s for OBC function %s variable %s\n", data->args[0], data->custom_tag, data->name);

  tsf->ntsfiles = data->nargs;
  tsf->tsfiles = hd_ts_multifile_read_us(master, data->nargs, data->args,
					 data->i_rule);
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

  for (cc = 1; cc <= open->no3_e1; cc++) {
    c = open->obc_e1[cc];
    c2cc->nmap[c] = cc;
  }
  for (cc = open->no3_e1 + 1; cc <= open->to3_e1; cc++) {
    c = open->obc_e1[cc];
    c2cc->tmap[c] = cc;
  }
  data->custdata = (void *)c2cc;
}

/*------------------------------------------------------------------*/
/* Read in rotated (u,v) 3D currents to transfer_u1av               */
/*------------------------------------------------------------------*/
void bf_uv_to_u1av_m(geometry_t *geom, master_t *master,
		     open_bdrys_t *open, bdry_details_t *data)
{
  int ee, e, es, c;
  int sb = 1; 
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[es], geom->u1y[es], geom->cellz[c],
	  tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int ee, e, es, c;
  int sb = 1;
  int eb = open->no3_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no3_e1 + 1;
    eb = open->to3_e1;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[es], geom->u1y[es], geom->cellz[c],
	  tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u2av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int eb = open->no3_e1;
  double depth = 0.0;
  double vav = 0.0;
  int *c2cc;
  uav_data_t *udata = (uav_data_t *)data->custdata;

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
  int ee, e, es, c;
  int sb = 1; 
  int eb = open->no2_e1;
  double thetau1, u, v;
  tsfiles_t *tsf;

  if(data->type & TAN) {
    sb = open->no2_e1 + 1;
    eb = open->to2_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[es], geom->u1y[es], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u1av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int ee, e, es, c;
  int sb = 1;
  int eb = open->no2_e1;
  double thetau1, u, v;
  tsfiles_t *tsf = (tsfiles_t *)data->custdata;

  if(data->type & TAN) {
    sb = open->no2_e1 + 1;
    eb = open->to2_e1;
  }

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t3d, geom->u1x[es], geom->u1y[es], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u2av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  int ee, e, es, c;
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
  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t, geom->u1x[es], geom->u1y[es], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "u", "v");
    open->transfer_u1[ee] = cos(thetau1) * u + sin(thetau1) * v;
  }
  /* Read the 2D velocities */
  if(data->type & TAN) {
    sb = open->no2_e1 + 1;
    eb = open->to2_e1;
  }
  tsf = (tsfiles_t *)data->custdata;
  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = geom->m2de[e];
    c = open->obc_e2[ee];
    thetau1 = geom->thetau1[es];
    u = 0.0, v = 0.0;
    getuv(master->t, geom->u1x[es], geom->u1y[es], geom->cellz[c],
          tsf->ntsfiles, tsf->tsfiles, tsf->filenames, &u, &v, "uav", "vav");
    open->transfer_u1av[ee] = cos(thetau1) * u + sin(thetau1) * v;
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
  double ramp = (master->rampf & FILEIN) ? master->rampval : 1.0;
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
    tinc = open->to3_e1;
    obc = open->obc_e1;
    zone = open->relax_zone_tan;
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


/*-------------------------------------------------------------------*/
/* Master initialisation routine for u1 river forcing                */
/*-------------------------------------------------------------------*/
void bf_gauss_init_m(master_t *master, open_bdrys_t *open,
		     bdry_details_t *d)
{
  tra_data_t *data = (tra_data_t *)malloc(sizeof(tra_data_t));
  geometry_t *geom = master->geom;
  int cc, c, c2;
  double d1;

  /* Initialise */
  memset(data, 0, sizeof(tra_data_t));
  d->custdata = (void *)data;

  /* Mid point of the boundary */
  /* R5000 
  data->tstart = 10.0;
  data->mid = open->no2_t / 2 + 1;
  */
  /* R200 
  data->tstart = 10.0 - 18.0/3600.0;
  data->mid = open->no2_t / 2;
  */

  data->tstart = 10.0 - 1.0/24.0;
  data->mid = open->no2_t / 2;

  data->x = d_alloc_1d(open->no2_t + 1);
  data->cells = i_alloc_1d(open->no2_t + 1);
  data->ncells = open->no2_t;
  strcpy(data->name, "passive");
  
  /* Get the distances along the boundary                            */
  d1 = 0.0;
  for (cc = data->mid; cc <= open->no2_t; cc++) {
    c = open->obc_t[cc];
    data->x[cc] = d1;
    data->cells[cc] = c;
    d1 += sqrt(geom->cellarea[c]);
  }
  data->dist = d1;

  cc = data->mid;
  c = open->obc_t[cc];
  d1 = sqrt(geom->cellarea[c]);
  for (cc = data->mid-1; cc >= 1; cc--) {
    c = open->obc_t[cc];
    data->x[cc] = d1;
    data->cells[cc] = c;
    d1 += sqrt(geom->cellarea[c]);
  }
}

/* END bf_gauss_init_m()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window initialisation routine for u1 river forcing                */
/*-------------------------------------------------------------------*/
void bf_gauss_init_w(geometry_t *window,  /* Window geometry structure */
                       open_bdrys_t *open,  /* OBC data */
                       bdry_details_t *d, /* Custom boundary data for
                                             window */
                       bdry_details_t *din  /* Custom boundary data from
                                               master */
  )
{
  tra_data_t *data = (tra_data_t *)malloc(sizeof(tra_data_t));
  tra_data_t *data_in = din->custdata;
  int cc, c, cn, i;
  double d1, d2;

  /* Initialise */
  memset(data, 0, sizeof(tra_data_t));
  d->custdata = (void *)data;
  data->dist = data_in->dist;
  data->mid = data_in->mid;
  data->x = d_alloc_1d(open->no2_t + 1);
  data->cells = i_alloc_1d(open->no2_t + 1);
  data->tstart = data_in->tstart;
  data->flag = 0x01;
  strcpy(data->name, data_in->name);
  for (cc = 1; cc <= open->no2_t; cc++) {
    c = open->obc_t[cc];
    if((cn = ANY(c, data_in->cells, data_in->ncells))) {
      data->cells[cc] = cn;
      data->x[cc] = data_in->x[cn];
    }
  }
}

/* END bf_gauss_init_w()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calling routine from the master for river forcing. This routine   */
/* reads the flow rate from file and stores it in the custom data    */
/* structure. The window and master share the same bdry_details_t    */
/* structure, hence all windows may also access this flow rate.      */
/*-------------------------------------------------------------------*/
void bf_gauss_m(geometry_t *geom, master_t *master,
		open_bdrys_t *open, bdry_details_t *d)
{
  return;
}

/* END bf_gauss_m()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calling routine from the window for u1 river forcing. This        */
/* calculates the river profile and returns a velocity for the u1    */
/* boundary variable. This routine uses the flow rate read from file */
/* on the master and stored in data->flow.                           */
/*-------------------------------------------------------------------*/
double bf_gauss_w(geometry_t *window, window_t *windat,
		  win_priv_t *wincon, open_bdrys_t *open, double t,
		  int c, int c2, bdry_details_t *d)
{
  tra_data_t *data = (tra_data_t *)d->custdata;
  int cc, cc2, tn, tm, j, e;
  double x, y, d1, vel, dist, gscale = 0.25 * data->dist;

  if (data->flag & 0x01) {
    if (windat->dum1) build_gauss_ref(window, windat, wincon, 
				      data->dist, gscale, data->tstart);
    data->flag = 0;
  }
  /* Find the tracer */
  /*
  tn = tracer_find_index(data->name, wincon->ntr, wincon->trinfo_3d);
  tm = open->trm[tn];
  */

  /* Find the index in the boundary array */
  for (cc = 1; cc <= open->no3_t; cc++) {
    if (c == open->obc_t[cc])
      break; 
  }
  for (cc2 = 1; cc2 <= open->no2_t; cc2++) {
    if (c2 == open->obc_t[cc2])
      break; 
  }

  /* Get the average inward boundary velocity */
  vel = d1 = 0.0;
  for (j = 1; j <= window->npe[c2]; j++) {
    if ((e = open->bec[j][cc]) && open->bcc[j][cc] > 0) {
      vel += windat->u1[e];
      d1 += 1.0;
    }
  }
  vel /= d1;

  /* Reconstruct the Gaussian function */
  d1 = 0.0;
  if (t >= data->tstart * 86400.0) {

    x = data->x[cc2];
    dist = vel * (t - data->tstart * 86400);
    if (dist < data->dist) 
      y = data->dist - dist;
    else if (dist > 2.0 * data->dist)
      y = 2.0 * data->dist;
    else
      y = dist - data->dist;

    d1 = wgt_gaussian_2d(x, y, gscale);    
  }
  return(d1);
}

void bf_gauss_t(master_t *master,  /* Master data */
               open_bdrys_t *open_w,  /* OBC structure in the window */
               bdry_details_t *data,  /* Custom data for this OBC */
               geometry_t *window,  /* Window geometry */
               window_t *windat /* Window data */
  )
{
  tra_data_t *d = (tra_data_t *)data->custdata;
  d->flag = 0x01;
  return;
}

void bf_gauss_free(master_t *master, bdry_details_t *d)
{
  tra_data_t *data = (tra_data_t *)d->custdata;
  d_free_1d(data->x);
  i_free_1d(data->cells);
  free((void **)d->custdata);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void build_gauss_ref(geometry_t *window,
		     window_t *windat,
		     win_priv_t *wincon,
		     double dist,
		     double gscale,
		     double tstart
		     )
{
  int c, cs, cc;
  int cm;
  int tn;
  double d1 = 0.0, d2, d3, d4;
  double *tr;
  double x0, y0;
  double d, x, y;
  double L2, Linf;
  double u, v;
  double ts = 10.5, te = 12.5;
  double tps = 11.5;
  double tpe = 12.5;
  int type = 5;
  int res = 500.0; /* Only for type = 4 */
  int as = 0; /* 0=ffsl, 1=quick, 2=vl */

  if (as == 0)
    tn = tracer_find_index("passivet", wincon->ntr, wincon->trinfo_3d);
  else
    tn = tracer_find_index("passive", wincon->ntr, wincon->trinfo_3d);
  tr =  windat->tr_wc[tn];

  if (type == 0) {
    /*---------------------------------------------------------------*/
    /* Use the location of maximum concentration as the location of  */
    /* the center for the analytic solution.                         */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      if (tr[c] > d1) {
	d1 = tr[c];
	cm = c;
      }
    }
    x0 = window->cellx[cm];
    y0 = window->celly[cm];
    /*
    tps = 11.125;
    tpe = 11.5;
    tps=11.709;
    tpe=11.75;
    */
  } else if (type == 4) {
    /*---------------------------------------------------------------*/
    /* Track the position of the quickest solution for each          */
    /* resolution from when the center enters the domain to 12.5     */
    /* days.                                                         */
    double ts, xs, ys, xe, ye, tinc;
    if (res == 5000) {
      ts = 10.75;
      xs = 2500.0;
      ys = 12500.0;
      xe = 37500.0;
      ye = 12500.0;
    }
    if (res == 2500) {
      /* ffsl: c = 52, ts = 10.5, te = 12.5 */
      /* q : c=52, ts = 10.5, te = 12.5 */
      /* vl: c = 52, ts = 10.375 te = 12 */
      ts = (as == 0) ? 10.5 : 10.375;
      te = (as == 0) ? 12.5 : 12.0;
      xs = 500.0;
      ys = 11500.0;
      xe = 38750.0;
      ye = 11250.0;
    }
    if (res == 1000) {
      /* ffsl: c = 468, ts = 10.5, te = 12.5 */
      /* q : c=469, ts = 10.5, te = 12.5 */
      /* vl: c = 469, ts = 10.375 te = 12.5 */
      ts = 10.5;
      te = 12.5;
      xs = 500.0;
      ys = 11500.0;
      xe = (as == 0) ? 36500.0 : 37500.0;
      ye = 11500.0;
    }
    if (res == 500) {
      ts = 10.5;
      xs = 250.0;
      ys = 12250.0;
      xe = (as == 0) ? 34250 : 35250.0;
      ye = 12750.0;
    }
    if (res == 200) {
      /*
      ts = 10.5;
      xs = 100.0;
      ys = 12300.0;
      */
      ts = 11.5;
      xs = 16900.0;
      ys = 13100.0;
      xe = 35900.0;
      ye = 13100.0;
    }

    tinc = 86400.0 * (te - ts);
    u = (xe - xs) / tinc;
    v = (ye - ys) / tinc;
    if (windat->days < ts) {
      x0 = xs;
      y0 = ys;
    } else {
      x0 = xs + u * 86400.0 * (windat->days - ts);
      y0 = ys + v * 86400.0 * (windat->days - ts);
      d1 = 1.0;
    }
  } else if (type == 5) {
    /*---------------------------------------------------------------*/
    /* Track the position of the 200m quickest solution from 11.5 to */
    /* 12.5 days.                                                    */
    double ts, xs, ys, xe, ye, tinc;
    ts = 11.5;
    te = 12.5;
    xs = 16900.0;
    ys = 13100.0;
    xe = 35900.0;
    ye = 13100.0;

    tinc = 86400.0 * (te - ts);
    u = (xe - xs) / tinc;
    v = (ye - ys) / tinc;
    if (windat->days < ts) {
      x0 = xs;
      y0 = ys;
    } else {
      x0 = xs + u * 86400.0 * (windat->days - ts);
      y0 = ys + v * 86400.0 * (windat->days - ts);
      d1 = 1.0;
    }
  } else {
    /*---------------------------------------------------------------*/
    /* Lagrangian tracking of the center using the velocity mean     */
    /* over the whole domain.                                        */
    if (type == 1) {
      u = v = 0.0;
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	u += windat->u[c];
	v += windat->v[c];
      }
      u /= (double)window->b3_t;
      v /= (double)window->b3_t;
    }
    /*---------------------------------------------------------------*/
    /* Lagrangian tracking of the center using the mean of depth     */
    /* average velocity.                                             */
    if (type == 2) { /* Mean of depth average */
      u = v = d1 = d2 = d3 = 0.0;
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	while (c != window->zm1[c]) {
	  u += (windat->u[c] * wincon->dz[c]);
	  v += (windat->v[c] * wincon->dz[c]);
	  d1 += wincon->dz[c];
	  c = window->zm1[c];
	}
	d2 += (u / d1);
	d3 += (v / d1);
      }
      u = d2 / (double)window->b2_t;
      v = d3 / (double)window->b2_t;
    }
    /*---------------------------------------------------------------*/
    /* Lagrangian tracking of the center using the 3D current        */
    if (type == 3) {
      u = windat->u[c];
      v = windat->v[c];
    }
    if (windat->days < tstart) {
      wincon->a1 = 1125;
      wincon->a2 = 500.0;
      wincon->a3 = 11500.0;
    } else {
      int i = -1, j = -1;
      c = cs = (int)wincon->a1;
      d1 = u * windat->dttr;
      d2 = v * windat->dttr;
      wincon->a2 += d1;
      wincon->a3 += d2;
      grid_xytoij(window->xyij_tree, wincon->a2, wincon->a3, &i, &j);
      wincon->a1 = (double)window->map[geom->nz - 1][j][i];
    }

    /* Get the cell with maximum concentration                       */
    x0 = wincon->a2;
    y0 = wincon->a3;
  }

  /*printf("%f %d : %f %f\n", windat->days,(int)wincon->a1,x0,y0);*/
  if (d1 == 0.0) return;

  /* Make a gaussian concentration centered on cm                    */
  d1 = d2 = d3 = d4 = 0.0;
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    windat->dum1[c] = windat->dum2[c] = 0.0;
    x = fabs(window->cellx[c] - x0);
    y = fabs(window->celly[c] - y0);
    d = sqrt(x * x + y * y);
    if (d <= dist) {
      windat->dum1[c] = wgt_gaussian_2d(x, y, gscale);
    }
    windat->dum2[c] = fabs(tr[c] - windat->dum1[c]);
    d1 += windat->dum2[c] * windat->dum2[c];
    d2 += windat->dum1[c] * windat->dum1[c];
    if (windat->dum2[c] > d3) d3 = windat->dum2[c];
    if (windat->dum1[c] > d4) d4 = windat->dum1[c];
  }
  /*
  L2 = sqrt(d1) / sqrt(d2);
  Linf = d3 / d4;
  */
  L2 = L2_norm(window, tr, windat->dum1);
  Linf = Linf_norm(window, tr, windat->dum1);

  if (windat->days > tps && windat->days <= tpe) 
    printf("%f %f %f\n",windat->days, L2, Linf);

}

double L2_norm(geometry_t *window, double *mod, double *ref)
{
  int c, cc;
  double L2, dif;
  double d1, d2, sa;

  d1 = d2 = sa = 0.0;
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    dif = fabs(mod[c] - ref[c]);
    sa += window->cellarea[c];
    d1 += dif * dif  * window->cellarea[c];
    d2 += ref[c] * ref[c] * window->cellarea[c];
  }
  L2 = sqrt(d1/sa) / sqrt(d2/sa);
  return L2;
}

double Linf_norm(geometry_t *window, double *mod, double *ref)
{
  int c, cc;
  double Linf, dif;
  double d1, d2, sa;

  d1 = d2 = sa = 0.0;
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    dif = fabs(mod[c] - ref[c]);
    if (dif > d1) d1 = dif;
    if (ref[c] > d2) d2 = ref[c];
  }
  Linf = d1 / d2;
  return Linf;
}

/* END build_gauss_ref()                                             */
/*-------------------------------------------------------------------*/

/* EOF */
