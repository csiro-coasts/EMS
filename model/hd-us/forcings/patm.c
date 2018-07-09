/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/patm.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric patm.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: patm.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t *ts;             /* Timeseries file */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
  int id;                       /* TS id for the variable */
} patm_data_t;

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
} patm_mdata_t;

/* Functions for reading the schedule the patm forcings. */
int patm_init_single(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  patm_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (patm_data_t *)malloc(sizeof(patm_data_t));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->patm, params->patm_dt,
			      "pressure", "Pa", &data->dt, &data->id,
			      NULL);

  if (data->ts == NULL) {
    geometry_t *geom = master->geom;
    int c, cc;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->patm[c] = master->ambpress;
    }
    free(data);
    return 0;
  }
  /* 
     data->dt=frc_get_input_dt(master->grid_dt,data->dt,"air pressure"); */
  return 1;
}


double patm_event_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  patm_data_t *data = (patm_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid(master, t, data->ts, data->id, master->patm, 1.0);
    event->next_event += data->dt;
  }

  return event->next_event;
}


void patm_cleanup_single(sched_event_t *event, double t)
{
  patm_data_t *data = (patm_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    hd_ts_free(master, data->ts);
    free(data);
  }
}


int patm_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  patm_mdata_t *data = NULL;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (patm_mdata_t *)malloc(sizeof(patm_mdata_t));
  schedSetPrivateData(event, data);
  data->tsfiles = frc_read_cell_ts_mult(master, params->patm, params->patm_dt,
					"pressure", "Pa", &data->dt, data->varids, 
					&data->ntsfiles, NULL, 1);

  if (data->tsfiles == NULL) {
    geometry_t *geom = master->geom;
    int c, cc;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->patm[c] = master->ambpress;
    }
    free(data);
    return 0;
  }
  /* 
     data->dt=frc_get_input_dt(master->grid_dt,data->dt,"air pressure"); */
  return 1;
}


double patm_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  patm_mdata_t *data = (patm_mdata_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid_mult(master, t, data->tsfiles, data->varids, data->ntsfiles, 
			  master->patm, 1.0);
    event->next_event += data->dt;
  }
  return event->next_event;
}


void patm_cleanup(sched_event_t *event, double t)
{
  patm_mdata_t *data = (patm_mdata_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}
