/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/airtemp.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric airtemp.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: airtemp.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} airtemp_data_t;

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
} airtemp_mdata_t;

/* Functions for reading the schedule the airtemp forcings. */
int airtemp_init_single(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  airtemp_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (airtemp_data_t *)malloc(sizeof(airtemp_data_t));
  schedSetPrivateData(event, data);
  if ((data->ts = frcw_read_cell_ts(master, params->airtemp, 
				   params->airtemp_dt,
				   "air_temp", "degrees C", &data->dt,
				    &data->id, &master->airtemp)) == NULL)
    data->ts = frc_read_cell_ts(master, params->airtemp, 
				params->airtemp_dt,
				"air_temp", "degrees_C", &data->dt,
				&data->id, &master->airtemp);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double airtemp_event_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  airtemp_data_t *data = (airtemp_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid(master, t, data->ts, data->id, master->airtemp, 1.0);
    event->next_event += data->dt;
  }

  return event->next_event;
}

void airtemp_cleanup_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  airtemp_data_t *data = (airtemp_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->airtemp);
    hd_ts_free(master, data->ts);
    free(data);
  }
}


int airtemp_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  airtemp_mdata_t *data = NULL;




  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (airtemp_mdata_t *)malloc(sizeof(airtemp_mdata_t));
  schedSetPrivateData(event, data);
  if((data->tsfiles = frc_read_cell_ts_mult(master, params->airtemp, params->airtemp_dt,
					   "air_temp", "degrees C", &data->dt, data->varids, 
					    &data->ntsfiles, &master->airtemp, 0)) == NULL)
    data->tsfiles = frc_read_cell_ts_mult(master, params->airtemp, params->airtemp_dt,
					  "air_temp", "degrees_C", &data->dt, data->varids, 
					  &data->ntsfiles, &master->airtemp, 1);
  if (data->tsfiles == NULL) {
    free(data);
    return 0;
  }
  return 1;
}


double airtemp_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  airtemp_mdata_t *data = (airtemp_mdata_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid_mult(master, t, data->tsfiles, data->varids, data->ntsfiles, 
			  master->airtemp, 1.0);
    event->next_event += data->dt;
  }
  return event->next_event;
}


void airtemp_cleanup(sched_event_t *event, double t)
{
  airtemp_mdata_t *data = (airtemp_mdata_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}
