/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/forcings.c
 *  
 *  Description:
 *  Read external forcings variable parameters and
 *  apply forcings at stipulated interval.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: forcings.c 7156 2022-07-07 02:31:40Z her127 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"


void forcings_init(master_t *master)
{
  sched_register(schedule, "forcings:wind", wind_init,
                 wind_event, wind_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:patm", patm_init,
                 patm_event, patm_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:precip", precip_init,
                 precip_event, precip_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:evap", evap_init,
                 evap_event, evap_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:airtemp", airtemp_init,
                 airtemp_event, airtemp_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:rh", rh_init,
                 rh_event, rh_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:wetbulb", wetbulb_init,
                 wetbulb_event, wetbulb_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:cloud", cloud_init,
                 cloud_event, cloud_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:swr", swr_init,
                 swr_event, swr_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:light", light_init,
                 light_event, light_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:webf", webf_init,
                 webf_event, webf_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:storm", storm_init,
                 storm_event, storm_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:eta_relax", eta_relax_init,
                 eta_relax_event, eta_relax_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:vel_relax", vel_relax_init,
                 vel_relax_event, vel_relax_cleanup, master, NULL, NULL);
  sched_register(schedule, "forcings:regulate", regulate_init,
                 regulate_event, regulate_cleanup, master, NULL, NULL);

  if (master->heatflux & NET_HEAT) {
    sched_register(schedule, "forcings:heatf", heatf_init,
                   heatf_event, heatf_cleanup, master, NULL, NULL);
  } else if (master->heatflux & COMP_HEAT) {
    sched_register(schedule, "forcings:longout", longout_init,
                   longout_event, longout_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:longin", longin_init,
                   longin_event, longin_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:sheat", sheat_init,
                   sheat_event, sheat_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:lheat", lheat_init,
                   lheat_event, lheat_cleanup, master, NULL, NULL);
  } else if (master->heatflux & COMP_HEAT_MOM) {
    sched_register(schedule, "forcings:swr_mean", swr_mean_init,
                   swr_mean_event, swr_mean_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:longout", longout_init,
                   longout_event, longout_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:sheat", sheat_init,
                   sheat_event, sheat_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:lheat", lheat_init,
                   lheat_event, lheat_cleanup, master, NULL, NULL);
  } else if (master->heatflux & COMP_HEAT_NONE) {
    sched_register(schedule, "forcings:swr_mean", swr_mean_init,
                   swr_mean_event, swr_mean_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:longout", longout_init,
                   longout_event, longout_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:sheat", sheat_init,
                   sheat_event, sheat_cleanup, master, NULL, NULL);
    sched_register(schedule, "forcings:lheat", lheat_init,
                   lheat_event, lheat_cleanup, master, NULL, NULL);
  } else if (master->heatflux & (SURF_RELAX|INVERSE)) {
    sched_register(schedule, "forcings:hftemp", hftemp_init,
                   hftemp_event, hftemp_cleanup, master, NULL, NULL);
  }
  if (master->heatflux & ADVANCED) {
    if (master->airtemp == NULL)
      hd_quit("Heat flux requires an air temperature file");
    if (master->patm == NULL)
      hd_quit("Heat flux requires an air pressure file");
    if (master->rh == NULL && master->wetb == NULL)
      hd_quit
        ("Heat flux requires an humidity (wet bulb or relative) file");
  }
}

void forcings_end(void)
{
  sched_deregister(schedule, "forcings:da");
  sched_deregister(schedule, "forcings:webf");
  sched_deregister(schedule, "forcings:swr");
  sched_deregister(schedule, "forcings:light");
  sched_deregister(schedule, "forcings:cloud");
  sched_deregister(schedule, "forcings:rh");
  sched_deregister(schedule, "forcings:airtemp");
  sched_deregister(schedule, "forcings:wetb");
  sched_deregister(schedule, "forcings:evap");
  sched_deregister(schedule, "forcings:precip");
  sched_deregister(schedule, "forcings:patm");
  sched_deregister(schedule, "forcings:wind");
  sched_deregister(schedule, "forcings:storm");
  sched_deregister(schedule, "forcings:eta_relax");
  sched_deregister(schedule, "forcings:heatf");
  sched_deregister(schedule, "forcings:hftemp");
  sched_deregister(schedule, "forcings:swan");
}


timeseries_t *frc_read_cell_ts_o(master_t *master, char *key,
				 char *varname, char *varunit, double *dt,
				 int *id, double **p)
{
  timeseries_t *ts;
  char buf[MAXLINELEN];
  char buf2[MAXSTRLEN];
  char key1[MAXLINELEN];
  geometry_t *geom = master->geom;

  /* Make key for input time step */
  sprintf(key1, "%s_INPUT_DT", key);

  /* Find time series name in parameter file */
  prm_set_errfn(hd_silent_warn);
  if (prm_read_char(master->prmfd, key, buf) <= 0)
    return NULL;

  prm_set_errfn(hd_quit);
  prm_get_time_in_secs(master->prmfd, key1, dt);
  ts = hd_ts_read(master, buf, 1);
  if ((*id = ts_get_index(ts, fv_get_varname(buf, varname, buf2))) < 0)
    hd_quit("frc_read_cell_ts: Can't find %s variable in %s\n", varname,
            ts->name);
  if (strcasecmp(ts->varunit[*id], varunit) != 0)
    hd_quit("frc_read_cell_ts: %s units must be %s\n", varname, varunit);
  if (p)
    *p = d_alloc_1d(geom->sgsizS);

  return ts;
}

timeseries_t *frc_read_cell_ts(master_t *master, char *fname, double i_dt,
			       char *varname, char *varunit, double *dt,
			       int *id, double **p)
{
  timeseries_t *ts;
  char buf[MAXLINELEN];
  geometry_t *geom = master->geom;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  if (strlen(fname) == 0)
    return NULL;

  *dt = i_dt;
  ts = hd_ts_read(master, fname, 1);
  if ((*id = ts_get_index(ts, fv_get_varname(fname, varname, buf))) < 0)
    hd_quit("frc_read_cell_ts: Can't find %s variable in %s\n", varname,
            ts->name);
  if (strcasecmp(ts->varunit[*id], varunit) != 0)
    hd_quit("frc_read_cell_ts: %s units must be %s\n", varname, varunit);
  if (p)
    *p = d_alloc_1d(geom->szcS);

  return ts;
}


timeseries_t **frc_read_cell_ts_mult(master_t *master, char *fname, double i_dt,
				     char *varname, char *varunit, double *dt,
				     int *id, int *ntsfiles, double **p, int quitmode)
{
  timeseries_t **ts;
  char buf[MAXLINELEN];
  geometry_t *geom = master->geom;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  if (strlen(fname) == 0)
    return NULL;

  *dt = i_dt;

  *ntsfiles = parseline(fname, (char **)files, MAXNUMTSFILES);
  filenames = (cstring *)malloc(sizeof(cstring) * *ntsfiles);
  for (t = 0; t < *ntsfiles; ++t)
     strcpy(filenames[t], ((char **)files)[t]);

  ts = hd_ts_multifile_read(master, *ntsfiles, filenames);
  if (ts == NULL) {
    free(filenames);
    free(ts);
    return NULL;
  }
  if (hd_ts_multifile_get_index(*ntsfiles, ts, filenames, varname, id) > 0) {
    for(t = 0; t < *ntsfiles; t++) {
      if (id[t] >= 0) {
	if (strcasecmp(ts[t]->varunit[id[t]], varunit) != 0) {
	  if (quitmode)
	    hd_quit("forcings: file %s %s units must be %s\n", filenames[t],
		    varname, varunit);
	  else {
	    hd_warn("forcings: file %s %s units must be %s\n", filenames[t],
		    varname, varunit);
	    return NULL;
	  }
	}
      } else
	hd_warn("forcings: Can't find variable '%s' in timeseries file %s.\n", 
		varname, filenames[t]);
    }
    hd_ts_multifile_check(*ntsfiles, ts, filenames, varname,
			  schedule->start_time, schedule->stop_time);
    if (p)
      *p = d_alloc_1d(geom->szcS);

  } else {
    if (quitmode)
      hd_quit("forcings: Can't find variable '%s' in any timeseries file %s.\n", 
	      varname, fname);
    else {
      hd_warn("forcings: Can't find variable '%s' in any timeseries file %s.\n", 
	      varname, fname);
      return NULL;
    }
  }
  return ts;
}


timeseries_t **frc_read_cell_ts_mult_us(master_t *master, char *fname, double i_dt,
					char *i_rule, char *varname, char *varunit, 
					double *dt, int *id, int *ntsfiles, double **p,
					int quitmode)
{
  timeseries_t **ts;
  char buf[MAXLINELEN];
  geometry_t *geom = master->geom;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  if (strlen(fname) == 0)
    return NULL;

  *dt = i_dt;

  *ntsfiles = parseline(fname, (char **)files, MAXNUMTSFILES);
  filenames = (cstring *)malloc(sizeof(cstring) * *ntsfiles);
  for (t = 0; t < *ntsfiles; ++t)
     strcpy(filenames[t], ((char **)files)[t]);

  if (strlen(i_rule))
    ts = hd_ts_multifile_read_us(master, *ntsfiles, filenames, i_rule);
  else
    ts = hd_ts_multifile_read(master, *ntsfiles, filenames);
  if (ts == NULL) {
    free(filenames);
    free(ts);
    return NULL;
  }
  if (hd_ts_multifile_get_index(*ntsfiles, ts, filenames, varname, id) > 0) {
    for(t = 0; t < *ntsfiles; t++) {
      if (id[t] >= 0) {
	if (strcasecmp(ts[t]->varunit[id[t]], varunit) != 0) {
	  if (quitmode)
	    hd_quit("forcings: file %s %s units must be %s\n", filenames[t],
		    varname, varunit);
	  else {
	    hd_warn("forcings: file %s %s units must be %s\n", filenames[t],
		    varname, varunit);
	    return NULL;
	  }
	}
      } else
	hd_warn("forcings: Can't find variable '%s' in timeseries file %s.\n", 
		varname, filenames[t]);
    }
    hd_ts_multifile_check(*ntsfiles, ts, filenames, varname,
			  schedule->start_time, schedule->stop_time);
    if (p)
      *p = d_alloc_1d(geom->szcS);

  } else {
    if (quitmode)
      hd_quit("forcings: Can't find variable '%s' in any timeseries file %s.\n", 
	      varname, fname);
    else {
      hd_warn("forcings: Can't find variable '%s' in any timeseries file %s.\n", 
	      varname, fname);
      return NULL;
    }
  }
  return ts;
}


timeseries_t *frcw_read_cell_ts(master_t *master, char *fname, double i_dt,
				char *varname, char *varunit, double *dt,
				int *id, double **p)
{
  timeseries_t *ts;
  char buf[MAXLINELEN];
  geometry_t *geom = master->geom;

  if (strlen(fname) == 0)
    return NULL;

  *dt = i_dt;
  ts = hd_ts_read(master, fname, 1);
  if ((*id = ts_get_index(ts, fv_get_varname(fname, varname, buf))) < 0) {
    hd_warn("frc_read_cell_ts: Can't find %s variable in %s\n", varname,
            ts->name);
    return NULL;
  }
  if (strcasecmp(ts->varunit[*id], varunit) != 0)
    hd_warn("frcw_read_cell_ts: %s units must be %s\n", varname, varunit);
  if (p)
    *p = d_alloc_1d(geom->szcS);

  return ts;
}

void frc_ts_eval_grid(master_t *master, double t, timeseries_t *ts, int id,
                      double *p, double conv)
{
  int cc, c;
  geometry_t *geom = master->geom;

  /* Sanity checks */
  if (ts == NULL)
    return;

  if (ts->nt <= 0)
    hd_quit("frc_ts_eval_grid: Bad time series\n");

  /* Evaluate values at cell centres */
  for (cc = 1; cc <= geom->b2_t; ++cc) {
    c = geom->w2_t[cc];
    p[c] = conv * ts_eval_xyz(ts, id, t, geom->cellx[c], geom->celly[c], 0.0);
  }
}


void frc_ts_eval_grid_mult(master_t *master, double t, timeseries_t **ts, int *id,
			   int ntsfiles, double *p, double conv)
{
  int c, cc;
  geometry_t *geom = master->geom;

  /* Sanity checks */
  if (ts == NULL)
    return;

  if (ts[0]->nt <= 0)
    hd_quit("frc_ts_eval_grid: Bad time series\n");

  /* If the file is ugrid, read it straight in                       */
  if (df_is_ugrid(ts[0]->df)) {
    frc_multifile_eval_ugrid2D(master, ntsfiles, ts, id,
			       p, t, geom->b2_t, master->thIO);
    return;
  }
  /*
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "eta", master->eta, tsin,
				vc2, nc2, (mode|VEL2D));

void hd_trans_multifile_eval(master_t *master, 
			     int ntsfiles, timeseries_t **tsfiles,
			     cstring * names, char *var, double *v, double t,
			     int *vec, int nvec, int mode)
    vc3 = geom->w2_t;
    nc3 = geom->b2_t;
    ve3 = geom->w2_e1;
    ne3 = geom->n2_e1;
  */

  /* Evaluate values at cell centres */
  for (cc = 1; cc <= geom->b2_t; ++cc) {
    c = geom->w2_t[cc];
    p[c] = conv * hd_ts_multifile_eval_xy(ntsfiles, ts, id,
				      t, geom->cellx[c], geom->celly[c]);
  }
}


/*-------------------------------------------------------------------*/
/* Routine to calculate a time series input time-step that is an     */
/* integral multiple of the model time-step. Input and output in s.  */
/*-------------------------------------------------------------------*/
double frc_get_input_dt(double dt,  /* Model time step */
                        double wdt, /* Desired time series time-step */
                        char *name  /* Name of the time series */
  )
{
  int n, n1, nt;
  int tol = 20;
  int m, m1, m2;
  char unt[MAXSTRLEN];
  double t;

  if (wdt < dt)
    return (dt);
  n = (int)(wdt / dt);
  n1 = tol * n / 100;
  m1 = m2 = n;
  for (m = 1; m <= n1; m++) {
    if (m1 * (int)dt % 60 == 0) {
      nt = (int)m1 *(int)dt;
      if (nt < 60) {
        strcpy(unt, "s");
        t = (double)nt;
      } else if (nt >= 60 && nt < 3600) {
        strcpy(unt, "min");
        t = (double)nt / 60;
      } else if (nt >= 3600 && nt < 86400) {
        strcpy(unt, "hr");
        t = (double)nt / 3600;
      } else {
        strcpy(unt, "day");
        t = (double)nt / 86400;
      }
      if ((double)m1 * dt != wdt)
        hd_warn("Changed input DT of time series %s to %6.3f %s\n", name,
                t, unt);
      return (m1 * (int)dt);
    }
    if (m2 * (int)dt % 60 == 0) {
      nt = (int)m2 *(int)dt;
      if (nt < 60) {
        strcpy(unt, "s");
        t = (double)nt;
      } else if (nt >= 60 && nt < 3600) {
        strcpy(unt, "min");
        t = (double)nt / 60;
      } else if (nt >= 3600 && nt < 86400) {
        strcpy(unt, "hr");
        t = (double)nt / 3600;
      } else {
        strcpy(unt, "day");
        t = (double)nt / 86400;
      }
      if ((double)m1 * dt != wdt)
        hd_warn("Changed input DT of time series %s to %6.3f %s\n", name,
                t, unt);
      return (m2 * (int)dt);
    }
    m1++;
    m2--;
  }
  hd_quit("Can't find input DT of time series %s divisable by %d\n", name,
          (int)dt);
  return (0);
}



/*-------------------------------------------------------------------*/
/* Truncates the precision to n bits                                 */
/*-------------------------------------------------------------------*/
void trunc_data(geometry_t *window,      /* Window geometry            */
	       window_t *windat,       /* Window data                */
	       win_priv_t *wincon,     /* Window constants           */
	      int n)
{
  int c, cc;
  double mask;

  mask = pow(10.0,n);
  for(cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    windat->eta[c] = (double)((int)(windat->eta[c] * mask) / mask);
    windat->u1av[c] = (double)((int)(windat->u1av[c] * mask) / mask);
    windat->u2av[c] = (double)((int)(windat->u2av[c] * mask) / mask);
  }
  for(cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    windat->sal[c] = (double)((int)(windat->sal[c] * mask) / mask);
    windat->temp[c] = (double)((int)(windat->temp[c] * mask) / mask);
    windat->u1[c] = (double)((int)(windat->u1[c] * mask) / mask);
    windat->u2[c] = (double)((int)(windat->u2[c] * mask) / mask);
  }
}

/* END truncate()                                                    */
/*-------------------------------------------------------------------*/
