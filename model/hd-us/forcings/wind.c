/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/wind.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric wind.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: wind.c 7362 2023-06-09 03:23:05Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

#define RHO_0 (1025.0) /* Reference density */

/* Prototypes */
typedef struct wind_data wind_data_t;

double windstress_orig(master_t *master, wind_data_t *data, 
		     double *wx, double *wy, int c);
double windstress_bunker(master_t *master, wind_data_t *data, 
		       double *wx, double *wy, int c);
double windstress_largepond(master_t *master, wind_data_t *data, 
			  double *wx, double *wy, int c);
double windstress_kitiag(master_t *master, wind_data_t *data, 
		       double *wx, double *wy, int c);
double windstress_kondo(master_t *master, wind_data_t *data, 
		      double *wx, double *wy, int c);
void wind_center(master_t *master, double *wx, double *wy);

struct wind_data{
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids_u[MAXNUMTSFILES];  /* Variable id's for u/e1 */
  int varids_v[MAXNUMTSFILES];  /* Variable id's for v/e1 */
  int u_id;                     /* TS id for u variable */
  int v_id;                     /* TS id for v variable */
  double wind_scale;            /* Scale factor for wind speeds */
  double dlv0;                  /* Drag law v0 */
  double dlv1;                  /* Drag law v1 */
  double dlc0;                  /* Drag law c0 */
  double dlc1;                  /* Drag law c1 */
  int type;                     /* Wind speed or stress */
  int wcs;                      /* Wind coordinate system */
  int neutral;                  /* Drag under neutral conditions */
  double (*windstress) (master_t *, wind_data_t *, double *, double *, int);
};


/* Functions for reading the schedule the wind forcings. */
int wind_init(sched_event_t *event)
{
  /*UR-FIX not used char buf[MAXSTRLEN]; */
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  wind_data_t *data = NULL;

  /* Read parameters */
  if (strlen(params->wind) == 0)
    return 0;

  data = (wind_data_t *)malloc(sizeof(wind_data_t));
  schedSetPrivateData(event, data);
  data->dt = params->wind_dt;

  data->ntsfiles = parseline(params->wind, (char **)files, MAXNUMTSFILES);
  filenames = (cstring *)malloc(sizeof(cstring) * data->ntsfiles);
  for (t = 0; t < data->ntsfiles; ++t)
     strcpy(filenames[t], ((char **)files)[t]);

  if (strlen(params->wind_interp))
    data->tsfiles = hd_ts_multifile_read_us(master, data->ntsfiles, filenames,
					    params->wind_interp);
  else
    data->tsfiles = hd_ts_multifile_read(master, data->ntsfiles, filenames);
  if (data->tsfiles == NULL) {
    free(filenames);
    free(data->tsfiles);
    free(data);
    return 0;
  }

  data->wcs = 0;
  /* Note: the name used in the input file for wind_u and wind_v */
  /* is defined in misc/filevars.c.                              */
  if ((hd_ts_multifile_get_index(data->ntsfiles, data->tsfiles,
                  filenames, "wind_u", data->varids_u) > 0) &&
      (hd_ts_multifile_get_index(data->ntsfiles, data->tsfiles,
                  filenames, "wind_v", data->varids_v) > 0)) {
      hd_ts_multifile_check(data->ntsfiles, data->tsfiles,
                            filenames, "wind_u",
                            schedule->start_time, schedule->stop_time);
      hd_ts_multifile_check(data->ntsfiles, data->tsfiles,
                            filenames, "wind_v",
                            schedule->start_time, schedule->stop_time);

  } else if ((hd_ts_multifile_get_index(data->ntsfiles, data->tsfiles,
                  filenames, "wind1", data->varids_u) > 0)) {
      hd_ts_multifile_check(data->ntsfiles, data->tsfiles,
                            filenames, "wind1",
                            schedule->start_time, schedule->stop_time);
      data->wcs = 1;
  } else {
    hd_quit("wind_init: Can't find wind variables ('wind_u','wind_v') or ('wind_e1', 'wind_e2') in timeseries file.\n");
  }

  data->wind_scale = params->wind_scale;
  data->dlv0 = params->dlv0;
  data->dlv1 = params->dlv1;
  data->dlc0 = params->dlc0;
  data->dlc1 = params->dlc1;
  data->type = params->wind_type;
  data->windstress = NULL;
  data->neutral = 0;
  if (params->stress_fn & BUNKER)
    data->windstress = windstress_bunker;
  else if (params->stress_fn & LARGEPOND)
    data->windstress = windstress_largepond;
  else if (params->stress_fn & KITIAG)
    data->windstress = windstress_kitiag;
  else if (params->stress_fn & KONDO) {
    data->windstress = windstress_kondo;
    if (strlen(params->airtemp) == 0)
      hd_warn("wind_init: AIRTEMP is required for WIND_STRESS_FCTN = KONDO.\n");
    if (strlen(params->wetb) == 0) {
      if (strlen(params->rh) == 0)
	hd_warn("wind_init: HUMIDITY is required for WIND_STRESS_FCTN = KONDO.\n");
      else
	hd_warn("wind_init: WETBULB is required for WIND_STRESS_FCTN = KONDO.\n");
    }
    data->neutral = params->neutral;
  }
  else
    data->windstress = windstress_orig;
  return 1;
}

double wind_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wind_data_t *data = (wind_data_t *)schedGetPrivateData(event);
  geometry_t *geom = master->sgrid;

  if (t >= (event->next_event - SEPS)) {
    int ii, i, cc, c, ee, e;
    double windx, windy;
    double cd;

    /*---------------------------------------------------------------*/
    /* Read the wind speed and convert to wind stress                */
    /*---------------------------------------------------------------*/
    if(data->type & SPEED) {
      if (!data->wcs) {
	/*-----------------------------------------------------------*/
	/* Read the u and v wind speed components and rotate onto    */
	/* the grid.                                                 */
	/* Wind e1 component at u1 points                            */
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  i = geom->w2_e1[ee];
	  windx = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_u,
				      t, geom->u1x[i], geom->u1y[i]);
	  windy = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_v,
				      t, geom->u1x[i], geom->u1y[i]);
	  cd = data->windstress(master, data, &windx, &windy, i);
	  master->wind1[i] = windx * geom->costhu1[i]
	    + windy * geom->sinthu1[i];
	  if (master->numbers & WIND_CD) master->wind_Cd[i] = cd;
	}
      } else {
	/*-----------------------------------------------------------*/
	/* Read the normal wind speed component and directly apply   */
	/* on the grid.                                              */
	/* Wind e1 component at u1 points                            */
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  i = geom->w2_e1[ee];
	  c = geom->e2c[i][0];
	  windx = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_u,
				      t, geom->u1x[i], geom->u1y[i]);


	  windy = 0.0;
	  cd = data->windstress(master, data, &windx, &windy, c);
	  master->wind1[i] = windx;
	  if (master->numbers & WIND_CD) master->wind_Cd[i] = cd;
	}
      }
      /*-------------------------------------------------------------*/
      /* Wind speed at cell edges                                    */
      for (ee = 1; ee <= geom->b2_e1; ++ee) {
	i = geom->w2_e1[ee];
	windx = data->wind_scale
          * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				    data->varids_u,
				    t, geom->u1x[i], geom->u1y[i]);
	windy = data->wind_scale
          * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				    data->varids_v,
				    t, geom->u1x[i], geom->u1y[i]);
	master->windspeed[i] = sqrt(windx * windx + windy * windy);
	master->winddir[i] = 0.0;
	if (master->windspeed[i] > 0.0) {
	  master->winddir[i] = acos(windy / master->windspeed[i]);
	  if (windx < 0)
	    master->winddir[i] *= -1.0;
	}
      }
      if (master->storm_dt) {
	memcpy(master->swind1, master->wind1, geom->sgsizS * sizeof(double));
      }
    } else {
      /*-------------------------------------------------------------*/
      /* Read the wind stress components directly                    */
      /*-------------------------------------------------------------*/
      if (!data->wcs) {
	/*-----------------------------------------------------------*/
	/* Read the u and v wind stress components and rotate onto   */
	/* the grid.                                                 */
	/* Wind stress e1 component at u1 points                     */
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  i = geom->w2_e1[ee];
	  windx = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_u,
				      t, geom->u1x[i], geom->u1y[i]);
	  windy = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_v,
				      t, geom->u1x[i], geom->u1y[i]);
	  master->wind1[i] = windx * geom->costhu1[i]
	    + windy * geom->sinthu1[i];

	  /* Wind speed at cell centres from wind stress             */
	  stresswind(&windx, &windy, data->dlv0, data->dlv1, data->dlc0,
		     data->dlc1);
	  master->windspeed[i] = sqrt(windx * windx + windy * windy);
	  master->winddir[i] = 0.0;
	  if (master->windspeed[i] > 0.0) {
	    master->winddir[i] = acos(windy / master->windspeed[i]);
	    if (windx < 0)
	      master->winddir[i] *= -1.0;
	  }
	}
      } else {
	/*-----------------------------------------------------------*/
	/* Read the normal wind stress components and directly apply */
	/* on the grid.                                              */
	/* Wind stress e1 component at u1 points                     */
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  i = geom->w2_e1[ee];
	  windx = data->wind_scale
	    * hd_ts_multifile_eval_xy(data->ntsfiles, data->tsfiles,
				      data->varids_u,
				      t, geom->u1x[i], geom->u1y[i]);
	  master->wind1[i] = windx;
	}
	/* Wind speed at cell centres from wind stress             */
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  int n, eoe;
	  double fs = 1.0;
	  i = geom->w2_e1[ee];
	  windx = master->wind1[i];
	  windy = 0.0;
	  for (n = 1; n <= geom->nee[i]; n++) {
	    eoe = geom->eSe[n][i];
	    fs = geom->h1au1[geom->m2de[eoe]] / geom->h2au1[i];
	    if (!eoe) continue;
	    windy += fs * geom->wAe[n][i] * master->wind1[eoe];
	  }
	  stresswind(&windx, &windy, data->dlv0, data->dlv1, data->dlc0,
		     data->dlc1);
	  master->windspeed[i] = sqrt(windx * windx + windy * windy);
	  master->winddir[i] = 0.0;
	  if (master->windspeed[i] > 0.0) {
	    master->winddir[i] = acos(windy / master->windspeed[i]);
	    if (windx < 0)
	      master->winddir[i] *= -1.0;
	  }
	}
      }
    }
    event->next_event += data->dt;

    /* Remove the contribution that goes into the waves if required. */
    /* Wave form drag (tau_w) and breaking disapation (tau_diss) are */
    /* input directly (e.g. from WWIII).                             */
    /* Note: tau_diss should be negative (overall we add the         */
    /* dissapation to stress); see webf.c.                           */
    if (master->waves & STOKES_DRIFT) {
      if (master->tau_w1 && master->tau_diss1) {
	double tau_w, tau_diss;
	for (ee = 1; ee <= geom->b2_e1; ++ee) {
	  i = geom->w2_e1[ee];
	  tau_w = vel_c2e(geom, master->tau_w1, master->tau_w2, i);
	  tau_diss = vel_c2e(geom, master->tau_diss1, master->tau_diss2, i);
	  master->wind1[i] -= (tau_w + tau_diss);
	}
      }
    }
    /* Wave form drag and breaking (whitecapping, depth induced wave */
    /* breaking and rollers) are computed by SWAN and passed to the  */
    /* wave variables. Wave form drag and wave surface streaming are */
    /* is also computed by SWAN and subtracted from wind stress.     */
    if (master->waves & NEARSHORE) {
      double tau_w, tau_diss;
      for (ee = 1; ee <= geom->b2_e1; ++ee) {
	i = geom->w2_e1[ee];
	tau_w = vel_c2e(geom, master->wave_wfdx, master->wave_wfdy, i);
	tau_diss = fabs(vel_c2e(geom, master->wave_fwcapx, master->wave_fwcapy, i));
	tau_diss += fabs(vel_c2e(geom, master->wave_fbrex, master->wave_fbrey, i));
	tau_diss += fabs(vel_c2e(geom, master->wave_fsurx, master->wave_fsury, i));
	/*master->wind1[i] -= ((tau_w + tau_diss) / RHO_0);*/
	master->wind1[i] -= (tau_w + tau_diss);
      }
    }

    /* Get the tangential component of the wind                      */
    for (ee = 1; ee <= geom->b2_e1; ee++) {
      double fs;
      e = geom->w2_e1[ee];
      master->wind2[e] = 0.0;
      for (i = 1; i <= geom->nee[e]; i++) {
	ii = geom->eSe[i][e];
	fs = geom->h1au1[ii] / geom->h2au1[e];
	if (!ii) continue;
	master->wind2[e] += fs * geom->wAe[i][e] * master->wind1[ii];
      }
    }

    /* Get the cell centered components of the wind                  */
#if defined(HAVE_WAVE_MODULE)
    if (master->do_wave & W_SWAN) {
      wind_center(master, master->swind1, master->swind2);
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	stresswind(&master->swind1[c], &master->swind2[c], 
		   data->dlv0, data->dlv1, data->dlc0, data->dlc1);
      }
    }
#endif
  }
  return event->next_event;
}

void wind_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wind_data_t *data = (wind_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}

/* Routine to calculate wind stress */
void windstress(double *wx, double *wy,
                double v0, double v1, double c0, double c1)
{
  double cd;
  double v = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  if (v <= v0)
    cd = c0;
  else if (v >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (v - v0) / (v1 - v0);

  *wx = cd * v * air_dens * (*wx);
  *wy = cd * v * air_dens * (*wy);
}

double windstress_orig(master_t *master, wind_data_t *data, 
		     double *wx, double *wy, int c)
{
  double v0 = data->dlv0;
  double v1 = data->dlv1;
  double c0 = data->dlc0;
  double c1 = data->dlc1;
  double cd;
  double v = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  if (v <= v0)
    cd = c0;
  else if (v >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (v - v0) / (v1 - v0);

  *wx = cd * v * air_dens * (*wx);
  *wy = cd * v * air_dens * (*wy);
  return(cd);
}

double windstress_bunker(master_t *master, wind_data_t *data, 
		       double *wx, double *wy, int c)
{
  double cd, c_h, c_e;
  double v = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  double ta = (master->airtemp) ? master->airtemp[c] : master->temp[c];
  cd = bunker(master->temp[c], ta, v, data->dlv0, &c_e, &c_h);

  *wx = cd * v * air_dens * (*wx);
  *wy = cd * v * air_dens * (*wy);
  return(cd);
}

double windstress_largepond(master_t *master, wind_data_t *data, 
			 double *wx, double *wy, int c)
{
  double cd, c_h, c_e;
  double v = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  double ta = (master->airtemp) ? master->airtemp[c] : master->temp[c];
  cd = large_pond(master->temp[c], ta, v, data->dlv0, &c_e, &c_h);

  *wx = cd * v * air_dens * (*wx);
  *wy = cd * v * air_dens * (*wy);
  return(cd);
}

double windstress_kitiag(master_t *master, wind_data_t *data, 
		       double *wx, double *wy, int c)
{
  double cd, c_h, c_e;
  double v = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  cd = kitiagorodskii(v, data->dlv0, &c_e, &c_h);

  *wx = cd * v * air_dens * (*wx);
  *wy = cd * v * air_dens * (*wy);
  return(cd);
}

double windstress_kondo(master_t *master, wind_data_t *data, 
		      double *wx, double *wy, int c)
{
  double cd, c_h, c_e;
  double sst = master->temp[c];
  double sal = master->sal[c];
  double ta = (master->airtemp) ? master->airtemp[c] : master->temp[c];
  double pres = (master->patm) ? master->patm[c] : master->ambpress;
  double wspd = sqrt((*wx) * (*wx) + (*wy) * (*wy));
  double dtw;            /* Wet bulb temperature                     */
  double qs = NOTVALID;  /* Specific humidity at surface (kg/kg)     */
  double q = NOTVALID;   /* Specific humidity at 10m (kg/kg)         */
  double es = 0.;        /* Vapour pressure at the air temp (HPa)    */
  double esat;     /* Saturation vapour pressure (HPa)               */
  double ew;       /* Vapour pressure at water temp. (Hpa)           */
  double rh;       /* Relative humidity (%)                          */

  /* Convert pressure to HPa                                         */
  pres /= 100.0;

  hd_quit("windstress kondo not supported. Please consider using a different scheme. eg. L&P\n");
  
  if (!data->neutral && master->sh_f !=NONE) {
    /* Saturation vapour pressure over the water. Since the wet bulb */
    /* temp. over water is not available, ew represents the maximum  */
    /* amount of water the air can accomodate ie. assume that wet    */
    /* temp. = SST and rel humitity = 100%                           */
    ew = sst * (sst * (sst * (sst * (sst * (sst * 6.136821e-11 + 2.034081e-8) +
				     3.03124e-6) + 2.650648e-4) +
		       1.428946e-2) + 0.4436519) + 6.1078;
    ew *= (1.0 - 0.000537 * sal); /* ew at salinity sal */

    /* Get the specific humidities (kg/kg)                           */
    qs = 0.622 * ew / (pres - (0.378 * ew));
    q = 0.622 * es / (pres - (0.378 * es));
  }
  cd = kondo(sst, ta, wspd, qs, q, data->dlv0, &c_e, &c_h);

  *wx = cd * wspd * air_dens * (*wx);
  *wy = cd * wspd * air_dens * (*wy);
  return(cd);
}

/* Routine to calculate wind stress */
void stresswind(double *wx, double *wy,
                double v0, double v1, double c0, double c1)
{
  double cd = 1.4e-3, w0;
  double sgn;

  sgn = (*wx > 0.0) ? 1.0 : -1.0;
  w0 = sqrt((fabs(*wx)) / (cd * air_dens));
  if (w0 <= v0)
    cd = c0;
  else if (w0 >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (w0 - v0) / (v1 - v0);
  *wx = sgn * sqrt((fabs(*wx)) / (cd * air_dens));

  sgn = (*wy > 0.0) ? 1.0 : -1.0;
  w0 = sqrt((fabs(*wy)) / (cd * air_dens));
  if (w0 <= v0)
    cd = c0;
  else if (w0 >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (w0 - v0) / (v1 - v0);
  *wy = sgn * sqrt((fabs(*wy)) / (cd * air_dens));
}

void wind_center(master_t *master, double *wx, double *wy)
{
  geometry_t *geom = master->geom;
  int cc, c, ee, e;
  double nu, nv;
  double a;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];

    wx[c] = 0.0;
    wy[c] = 0.0;
    nu = nv = 0.0;

    for (ee = 1; ee <= geom->npe[c]; ee++) {
      e = geom->c2e[ee][c];
      a = 0.5 * geom->h1au1[e] * geom->h2au1[e];
      /* Get the cell centered east and north velocity               */
      wx[c] += a * (master->wind1[e] * geom->costhu1[e] + master->wind2[e] * geom->costhu2[e]);
      nu += a;
      wy[c] += a * (master->wind1[e] * geom->sinthu1[e] + master->wind2[e] * geom->sinthu2[e]);
      nv += a;
    }
    wx[c] = (nu) ? wx[c] / nu : 0.0;
    wy[c] = (nv) ? wy[c] / nv : 0.0;
  }
}
