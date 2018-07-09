/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/cloud.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric cloud.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: cloud.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} cloud_data_t;

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
} cloud_mdata_t;

/* Functions for reading the schedule the cloud forcings. */
int cloud_init_single(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  cloud_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (cloud_data_t *)malloc(sizeof(cloud_data_t));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->cloud, params->cloud_dt,
                              "cloud", "oktas", &data->dt, &data->id,
                              &master->cloud);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double cloud_event_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  cloud_data_t *data = (cloud_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    /* Cloud cover, convert to fraction (0-1) */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->cloud, 0.125);
    event->next_event += data->dt;
  }

  return event->next_event;
}

void cloud_cleanup_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  cloud_data_t *data = (cloud_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->cloud);
    hd_ts_free(master, data->ts);
    free(data);
  }
}


int cloud_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  cloud_mdata_t *data = NULL;
  char files[MAXNUMTSFILES][MAXSTRLEN];
  cstring *filenames;
  int t;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (cloud_mdata_t *)malloc(sizeof(cloud_mdata_t));
  schedSetPrivateData(event, data);
  data->tsfiles = frc_read_cell_ts_mult(master, params->cloud, params->cloud_dt,
					"cloud", "oktas", &data->dt, data->varids, 
					&data->ntsfiles, &master->cloud, 1);

  if (data->tsfiles == NULL) {
    free(data);
    return 0;
  }
  return 1;
}


double cloud_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  cloud_mdata_t *data = (cloud_mdata_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid_mult(master, t, data->tsfiles, data->varids, data->ntsfiles, 
			  master->cloud, 0.125);
    event->next_event += data->dt;
  }
  return event->next_event;
}


void cloud_cleanup(sched_event_t *event, double t)
{
  cloud_mdata_t *data = (cloud_mdata_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}
