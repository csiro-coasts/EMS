/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/webf.c
 *  
 *  Description:
 *  Apply wave enhanced bottom friction.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: webf.c 7361 2023-06-09 03:22:40Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include <time.h>
#include "hd.h"

#define RHO_0 (1025.0) /* Reference density */

int amp_id, per_id, dir_id, ub_id, fx_id, fy_id, ste1_id, ste2_id;
int aw1_id, aw2_id, wo1_id, wo2_id, ss1_id, ss2_id;

/* Multi-file wave input data */
typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **ts;            /* Timeseries file */
  int ub_id[MAXNUMTSFILES];     /* TS id for orbital velocity variable */
  int period_id[MAXNUMTSFILES]; /* TS id for orbital period variable */
  int dir_id[MAXNUMTSFILES];    /* TS id for orbital direction variable */
  int amp_id[MAXNUMTSFILES];    /* TS id for wave amplitude variable */
  int fx_id[MAXNUMTSFILES];     /* TS id for wave forcing variables */
  int fy_id[MAXNUMTSFILES];     /* TS id for wave forcing variables */
  int ste1_id[MAXNUMTSFILES];   /* TS id for stokes forcing variables */
  int ste2_id[MAXNUMTSFILES];   /* TS id for stokes forcing variables */
  int aw1_id[MAXNUMTSFILES];    /* TS id for wave stress forcing variables */
  int aw2_id[MAXNUMTSFILES];    /* TS id for wave stress forcing variables */
  int wo1_id[MAXNUMTSFILES];    /* TS id for wave breaking forcing variables */
  int wo2_id[MAXNUMTSFILES];    /* TS id for wave breaking forcing variables */
  int **ss1_id;                 /* Ids for spectral Sokes drift */
  int **ss2_id;                 /* Ids for spectral Sokes direction */
  int wcs;                      /* Wave force coordinate system */
  int stcs;                     /* Stokes coordinate system */
  char interp_type[MAXSTRLEN];  /* Grid library interp type to use */
} webf_data_t;

int webf_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  webf_data_t *data = NULL;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  char buf[MAXSTRLEN];
  int t, n;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);
  amp_id = per_id = dir_id = ub_id = fx_id = fy_id = ste1_id = ste2_id = -1;
  aw1_id = aw2_id = wo1_id = wo2_id = -1;

  /* Read parameters */
  if (strlen(params->webf) == 0 || params->webf_dt == 0.0)
    return 0;

  data = (webf_data_t *)malloc(sizeof(webf_data_t));
  schedSetPrivateData(event, data);
  data->wcs = 0;
  data->stcs = 0;
  sprintf(data->interp_type, "%c", '\0');
  if (strlen(params->webf_interp))
    sprintf(data->interp_type, "%s", params->webf_interp);

  /* Open multifile timeseries */
  data->dt = params->webf_dt;
  data->ntsfiles = parseline(params->webf, (char **)files, MAXNUMTSFILES);
  filenames = (cstring *)malloc(sizeof(cstring) * data->ntsfiles);
  for (t = 0; t < data->ntsfiles; ++t)
     strcpy(filenames[t], ((char **)files)[t]);

  /* Initialise all ids */
  for (t=0; t<MAXNUMTSFILES; t++) {
    data->ub_id[t] = -1; 
    data->period_id[t] = -1;
    data->dir_id[t] = -1; 
    data->amp_id[t] = -1;
    data->fx_id[t] = -1;
    data->fy_id[t] = -1;
    data->ste1_id[t] = -1;
    data->ste2_id[t] = -1;
    data->aw1_id[t] = -1;
    data->aw2_id[t] = -1;
    data->wo1_id[t] = -1;
    data->wo2_id[t] = -1;
  }

  data->ts = hd_ts_multifile_read(master, data->ntsfiles, filenames);
  if (data->ts) { 
    for (t = 0; t < data->ntsfiles; t++) {
      if ((data->ub_id[t] = ts_get_index(data->ts[t],
					 fv_get_varname(filenames[t],
							"ub",buf))) < 0)
	hd_warn("webf_init: Can't find ub variable in %s\n",
		data->ts[t]->name);

      if ((data->amp_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "amplitude",buf))) < 0)
	hd_warn("webf_init: Can't find amplitude variable in %s\n",
		data->ts[t]->name);

      if ((data->period_id[t] = ts_get_index(data->ts[t],
					     fv_get_varname(filenames[t],
							    "period",buf))) < 0)

      hd_warn("webf_init: Can't find period variable in %s\n",
	      data->ts[t]->name);

      if ((data->dir_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "direction",buf))) < 0)
        hd_warn("webf_init: Can't find direction variable in %s\n",
		data->ts[t]->name);

      if (((data->fx_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "force_e1",buf))) >= 0) &&
	  ((data->fy_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "force_e2",buf))) >= 0))
	 data->wcs = 1;
      else if (((data->fx_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "force_x",buf))) < 0) ||
               ((data->fy_id[t] = ts_get_index(data->ts[t],
					  fv_get_varname(filenames[t],
							 "force_y",buf))) < 0))
           hd_warn("webf_init: Can't find wave forcing variables in %s\n",
		data->ts[t]->name);

      if (((data->ste1_id[t] = ts_get_index(data->ts[t],
					    fv_get_varname(filenames[t],
							   "stokes_e1",buf))) >= 0) &&
	  ((data->ste2_id[t] = ts_get_index(data->ts[t],
					    fv_get_varname(filenames[t],
							   "stokes_e2",buf))) >= 0))
	data->stcs = 1;
      else if (((data->ste1_id[t] = ts_get_index(data->ts[t],
						 fv_get_varname(filenames[t],
								"stokes_x",buf))) >= 0) &&
               ((data->ste2_id[t] = ts_get_index(data->ts[t],
						 fv_get_varname(filenames[t],
								"stokes_y",buf))) >= 0))
	data->stcs = 0;
      else
	hd_warn("webf_init: Can't find stokes velocity variables in %s\n",
		data->ts[t]->name);

      if (((data->aw1_id[t] = ts_get_index(data->ts[t],
					   fv_get_varname(filenames[t],
							  "utaw",buf))) < 0) ||
	  ((data->aw2_id[t] = ts_get_index(data->ts[t],
					    fv_get_varname(filenames[t],
							   "vtaw",buf))) < 0))
	hd_warn("webf_init: Can't find wave supported wind stress variables in %s\n",
		data->ts[t]->name);

      if (((data->wo1_id[t] = ts_get_index(data->ts[t],
					   fv_get_varname(filenames[t],
							  "utwo",buf))) < 0) ||
	  ((data->wo2_id[t] = ts_get_index(data->ts[t],
					    fv_get_varname(filenames[t],
							   "vtwo",buf))) < 0))
	hd_warn("webf_init: Can't find wave to ocean stress variables in %s\n",
		data->ts[t]->name);

      if (master->waves & SPECTRAL) {
	size_t nsfr = 0;
	size_t start[4];
	size_t count[4];
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;
	count[0] = nsfr;
	count[1] = 0;
	count[2] = 0;
	count[3] = 0;
	master->nsfr = (int)nsfr;
	master->freq = d_alloc_1d(master->nsfr);
	/*
	nc_inq_dimlen(cdfid, ncw_dim_id(cdfid, "nfreq"), &nsfr);
	nc_get_vara_double(cdfid, ncw_var_id(cdfid, "frequency"), start, count,
			   master->freq);
	*/
	hd_quit("webf_init: SPECTRAL not supported\n");
	
	data->ss1_id = i_alloc_2d(MAXNUMTSFILES, master->nsfr);
	data->ss2_id = i_alloc_2d(MAXNUMTSFILES, master->nsfr);
	for (n = 0; n < master->nsfr; n++) {
	    char sname[MAXSTRLEN], dname[MAXSTRLEN];
	  sprintf(sname, "stokes_x%d", n);
	  sprintf(dname, "stokes_y%d", n);
	  if (((data->ss1_id[n][t] = ts_get_index(data->ts[t],
					       fv_get_varname(filenames[t],
							      sname, buf))) < 0) ||
	    ((data->ss2_id[n][t] = ts_get_index(data->ts[t],
					     fv_get_varname(filenames[t],
							    dname, buf))) < 0))
	    hd_warn("webf_init: Can't find spectral Stokes drift variables in %s\n",
		    data->ts[t]->name);
	}
      }
    }
    for (t = 0; t < data->ntsfiles; t++) {
      if (data->amp_id[t] >= 0 && amp_id == -1)
	amp_id = data->amp_id[t];
      if (data->period_id[t] >= 0 && per_id == -1)
	per_id = data->period_id[t];
      if (data->dir_id[t] >= 0 && dir_id == -1)
	dir_id = data->dir_id[t];
      if (data->ub_id[t] >= 0 && ub_id == -1)
	ub_id = data->ub_id[t];
      if (data->fx_id[t] >= 0 && fx_id == -1)
	fx_id = data->fx_id[t];
      if (data->fy_id[t] >= 0 && fy_id == -1)
	fy_id = data->fy_id[t];
      if (data->ste1_id[t] >= 0 && ste1_id == -1)
	ste1_id = data->ste1_id[t];
      if (data->ste2_id[t] >= 0 && ste2_id == -1)
	ste2_id = data->ste2_id[t];
      if (data->aw1_id[t] >= 0 && aw1_id == -1)
	aw1_id = data->aw1_id[t];
      if (data->aw2_id[t] >= 0 && aw2_id == -1)
	aw2_id = data->aw2_id[t];
      if (data->wo1_id[t] >= 0 && wo1_id == -1)
	wo1_id = data->wo1_id[t];
      if (data->wo2_id[t] >= 0 && wo2_id == -1)
	wo2_id = data->wo2_id[t];
      for (n = 0; n < master->nsfr; n++) {
	if (data->ss1_id[n][t] >= 0 && ss1_id == -1)
	  ss1_id = data->ss1_id[n][t];
	if (data->ss1_id[n][t] >= 0 && ss1_id == -1)
	  ss2_id = data->ss2_id[n][t];
      }
    }
  } else {
    /* data-ts == NULL */
    free(data);
    return 0;
  }
  return 1;
}


double webf_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  webf_data_t *data = (webf_data_t *)schedGetPrivateData(event);
  geometry_t *geom = master->geom;

  if (t >= (event->next_event - SEPS)) {
    int ii, i;
    double wavex, wavey;

    if (strlen(data->interp_type)) {
      /* Use the grid library as per the interpolation type          */
      /* specified.                                                  */
      if (ub_id >= 0)
	hd_ts_grid_interp_multifile(master, data->ts, data->ntsfiles, data->ub_id,
				    master->wave_ub, t, geom->w2_t, geom->b2_t, 
				    data->interp_type);
      if (amp_id >= 0)
	hd_ts_grid_interp_multifile(master, data->ts, data->ntsfiles, data->amp_id,
				    master->wave_amp, t, geom->w2_t, geom->b2_t, 
				    data->interp_type);
      if (per_id >= 0)
	hd_ts_grid_interp_multifile(master, data->ts, data->ntsfiles, data->period_id,
				    master->wave_period, t, geom->w2_t, geom->b2_t, 
				    data->interp_type);
      if (dir_id >= 0)
	hd_ts_grid_interp_multifile(master, data->ts, data->ntsfiles, data->dir_id,
				    master->wave_dir, t, geom->w2_t, geom->b2_t, 
				    data->interp_type);
    } else {
      /* The standard EMS library functions                          */
      if (ub_id >= 0)
	frc_ts_eval_grid_mult(master, t, data->ts, data->ub_id,data->ntsfiles, 
			      master->wave_ub, 1.0);
      if (amp_id >= 0)
	frc_ts_eval_grid_mult(master, t, data->ts, data->amp_id,data->ntsfiles, 
			      master->wave_amp, 1.0);
      if (per_id >= 0)
	frc_ts_eval_grid_mult(master, t, data->ts, data->period_id, 
			      data->ntsfiles, master->wave_period, 1.0);
      if (dir_id >= 0)
	frc_ts_eval_grid_mult(master, t, data->ts, data->dir_id,
			      data->ntsfiles, master->wave_dir, 1.0);
    }

    /* Read the wave-induced forcing data, and if necessary rotate   */
    /* onto the grid                                                 */
    if (fx_id >= 0 && fy_id >= 0) {
      if (!data->wcs) {
	/* Wave parameters are stored as east and north cell         */
	/* centered components, and rotated onto the edges as        */
	/* required. Save the east and north components directly.    */
        for (ii = 1; ii < geom->b2_t; ++ii) {
          i = geom->w2_t[ii];
          wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->fx_id,
                                          t, geom->cellx[i], geom->celly[i]);
          wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->fy_id,
                                          t, geom->cellx[i], geom->celly[i]);
          master->wave_Fx[i];
          master->wave_Fy[i];
	}
      } else {
	/* Rotate the components to east and north, using the cell   */
	/* centered angle.                                           */
        for (ii = 1; ii < geom->b2_t; ++ii) {
          i = geom->w2_t[ii];
          wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->fx_id,
                                          t, geom->cellx[i], geom->celly[i]);
          wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->fy_id,
                                          t, geom->cellx[i], geom->celly[i]);
          master->wave_Fx[i] = wavex * geom->costhcell[i] - wavey * geom->sinthcell[i];
          master->wave_Fy[i] = wavex * geom->sinthcell[i] + wavey * geom->costhcell[i];
	}
      }
    }

    /* Read the stokes velocity forcing data, and if necessary       */
    /* rotate onto the grid.                                         */
    if (master->waves & STOKES && ste1_id >= 0 && ste2_id >= 0) {
      if (!data->stcs) {
        for (ii = 1; ii < geom->b2_t; ++ii) {
          i = geom->w2_t[ii];
          wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->ste1_id,
                                          t, geom->cellx[i], geom->celly[i]);
          wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->ste2_id,
                                          t, geom->cellx[i], geom->celly[i]);
          master->wave_ste1[i] = wavex;
          master->wave_ste2[i] = wavey;
	}
      } else {
        for (ii = 1; ii < geom->b2_t; ++ii) {
          i = geom->w2_t[ii];
          wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->ste1_id,
                                          t, geom->cellx[i], geom->celly[i]);
          wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->ste2_id,
                                          t, geom->cellx[i], geom->celly[i]);
          master->wave_ste1[i] = wavex * geom->costhcell[i] - wavey * geom->sinthcell[i];
          master->wave_ste2[i] = wavex * geom->sinthcell[i] + wavey * geom->costhcell[i];
	}
      }
    }

    /* Read the wave supported wind friction velocities and convert  */
    /* to a stress, rotating onto the grid.                          */
    if (master->waves & STOKES_WIND && aw1_id >= 0 && aw2_id >= 0) {
      /* Only east and noth vector component import currently        */
      for (ii = 1; ii < geom->b2_t; ++ii) {
	i = geom->w2_t[ii];
	wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->aw1_id,
					t, geom->cellx[i], geom->celly[i]);
	wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->aw2_id,
					t, geom->cellx[i], geom->celly[i]);
	master->tau_w1[i] = RHO_0 * wavex;
	master->tau_w2[i] = RHO_0 * wavey;
      }
      /* Rotate to east and north
        for (ii = 1; ii < geom->b2_t; ++ii) {
	  i = geom->w2_t[ii];
	  wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->aw1_id,
					  t, geom->cellx[i], geom->celly[i]);
	  wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->aw2_id,
					  t, geom->cellx[i], geom->celly[i]);
          master->tau_w1[i] = RHO_0 * (wavex * geom->costhcell[i] - wavey * geom->sinthcell[i]);
          master->tau_w2[i] = RHO_0 * (wavex * geom->sinthcell[i] + wavey * geom->costhcell[i]);
	}
      */
    }

    /* Read the wave to ocean friction velocities and convert to a   */
    /* stress, rotating onto the grid.                               */
    /* Note: stress = dens * (friction velocity)**2                  */
    /* tau_w & tau_diss2 have units of m2s-2 (friction velocity**2)  */
    /* These stresses should be negative.                            */
    if (master->waves & STOKES_WIND && wo1_id >= 0 && wo2_id >= 0) {
      /* Only east and noth vector component import currently        */
      for (ii = 1; ii < geom->b2_t; ++ii) {
	i = geom->w2_t[ii];
	wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->wo1_id,
					t, geom->cellx[i], geom->celly[i]);
	wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->wo2_id,
					t, geom->cellx[i], geom->celly[i]);
	master->tau_diss1[i] = -RHO_0 * wavex;
	master->tau_diss2[i] = -RHO_0 * wavey;
      }
      /* Rotate to east and north
        for (ii = 1; ii < geom->b2_t; ++ii) {
	  i = geom->w2_t[ii];
	  wavex = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->wo1_id,
					  t, geom->cellx[i], geom->celly[i]);
	  wavey = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, data->wo2_id,
					  t, geom->cellx[i], geom->celly[i]);
          master->tau_diss1[i] = RHO_0 * (wavex * geom->costhcell[i] - wavey * geom->sinthcell[i]);
          master->tau_diss2[i] = RHO_0 * (wavex * geom->sinthcell[i] + wavey * geom->costhcell[i]);
	}
      */
    }

    /* Read the Stokes spectrum and Stokes spectrum dimension        */
    if (master->waves & STOKES && master->waves & SPECTRAL && ss1_id >= 0 && ss2_id >= 0) {
      memset(master->wave_ste1, 0, geom->sgsizS * sizeof(double));
      memset(master->wave_ste2, 0, geom->sgsizS * sizeof(double));
      memset(master->wave_stke1, 0, geom->sgsiz * sizeof(double));
      memset(master->wave_stke2, 0, geom->sgsiz * sizeof(double));
      for (ii = 1; ii < geom->b2_t; ++ii) {
	int c;
	double drift, dir;
	i = c = geom->w2_t[ii];
	for (ii = 0; ii < master->nsfr; ii++) {
	  int *s_id = data->ss1_id[ii];
	  int *d_id = data->ss2_id[ii];
	  drift = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, s_id,
					  t, geom->cellx[i], geom->celly[i]);
	  dir = hd_ts_multifile_eval_xy(data->ntsfiles, data->ts, d_id,
					t, geom->cellx[i], geom->celly[i]);
	  master->wave_ste1[i] += drift * sin(dir);
	  master->wave_ste2[i] += drift * cos(dir);
	  while (c != geom->zm1[c]) {
	    double depth = min(0.0, geom->cellz[c]);
	    double dz = geom->gridz[c] - geom->gridz[geom->zm1[c]];
	    double attn = exp(master->freq[ii] * depth);
	    if (attn * dz < 100)
	      attn *= sinh(0.5 * attn / dz) / (0.5 * attn / dz);
	    master->wave_stke1[c] += attn * drift * sin(dir);
	    master->wave_stke2[c] += attn * drift * cos(dir);
	  }
	}
      }
    }
    event->next_event += data->dt;
  }
  return event->next_event;
}


void webf_cleanup(sched_event_t *event, double t)
{
  webf_data_t *data = (webf_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->ts[i]);
    free(data);
  }
}

/* Called by the waves library */
void w_check_wave_data(int *wa, int *wp, int *wd, int *wu, int *wif, int *stv)
{
  *wa = *wp = *wd = *wu = *wif = *stv = 0;
  if (per_id >= 0) *wp = 1;
  if (amp_id >= 0) *wa = 1;
  if (dir_id >= 0) *wd = 1;
  if (ub_id >= 0) *wu = 1;
  if (fx_id >= 0 && fy_id >= 0) *wif = 1;
  if (ste1_id >= 0 && ste2_id >= 0) *stv = 1;
  /*
  if (aw1_id >= 0 && aw2_id >= 0) *aw = 1;
  if (wo1_id >= 0 && wo2_id >= 0) *wo = 1;
  */
}

/* Called by run_setup to write the setup.txt file */
void webf_write_setup(FILE *fp, sched_event_t * event)
{
  char buf[MAXSTRLEN];
  char files[MAXNUMTSFILES][MAXSTRLEN];
  int t;
  int count = 0;
  webf_data_t *data = (webf_data_t *)schedGetPrivateData(event);
  
  fprintf(fp,"Waves variables read from %d file(s):\n", data->ntsfiles);
  for (t=0; t<data->ntsfiles; t++) {
    timeseries_t *ts = data->ts[t];
    fprintf(fp, "\t%d : %s", t+1, ts->name);
    /*
     * Print out all variables that are read from file and their substitutions, if any
     */
    if (ub_id >= 0) {
      fprintf(fp, "(%s)", ts->varname[data->ub_id[t]]);
      count++;
    }
    
    if (per_id >= 0) {
      fprintf(fp, "(%s)", ts->varname[data->period_id[t]]);
      count++;
    }
    
    if (amp_id >= 0) {
      fprintf(fp, "(%s)", ts->varname[data->amp_id[t]]);
      count++;
    }
    
    if (dir_id >= 0) {
      fprintf(fp, "(%s)", ts->varname[data->dir_id[t]]);
      count++;
    }
    fprintf(fp,"\n");
  }

  if (fx_id >= 0 && fy_id >= 0) {
    fprintf(fp, "\tRadiation stresses read\n");
    count++;
  }

  if (ste1_id >= 0 && ste2_id >= 0) {
    fprintf(fp, "\tStokes velocity read\n");
    count++;
  }

  if (aw1_id >= 0 && aw2_id >= 0) {
    fprintf(fp, "\tWave supported wind stress read\n");
    count++;
  }
  if (wo1_id >= 0 && wo2_id >= 0) {
    fprintf(fp, "\tWave to ocean stress read\n");
    count++;
  }
  if (ste1_id >= 0 && ste2_id >= 0 && aw1_id < 0 && aw2_id < 0 &&
      wo1_id < 0 && wo2_id < 0)
      fprintf(fp, "\tNote: Atmospheric stress is NOT compensated by the contribution to the waves.\n");

  if (!count)
    fprintf(fp, "  NO wave variables read from file\n");

  if (strlen(data->interp_type))
    fprintf(fp,"\nWaves variables interpolation scheme : %s\n", 
	    data->interp_type);
}



