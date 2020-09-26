/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/precip.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric precip.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: precip.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} precip_data_t;

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
} precip_mdata_t;

/* Functions for reading the schedule the precip forcings. */
int precip_init_single(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  precip_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (precip_data_t *)malloc(sizeof(precip_data_t));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->precip, params->precip_dt,
                              "precipitation", "mm day-1", &data->dt,
                              &data->id, &master->precip);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double precip_event_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  precip_data_t *data = (precip_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    /* Precipitation, convert to m/s */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->precip,
                     1.0 / (1000.0 * 86400.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

void precip_cleanup_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  precip_data_t *data = (precip_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->precip);
    hd_ts_free(master, data->ts);
    free(data);
  }
}


int precip_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  precip_mdata_t *data = NULL;


  

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (precip_mdata_t *)malloc(sizeof(precip_mdata_t));
  schedSetPrivateData(event, data);
  data->tsfiles = frc_read_cell_ts_mult(master, params->precip, params->precip_dt,
					"precipitation", "mm day-1", &data->dt, data->varids, 
					&data->ntsfiles, &master->precip, 1);

  if (data->tsfiles == NULL) {
    free(data);
    return 0;
  }
  return 1;
}


double precip_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  precip_mdata_t *data = (precip_mdata_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid_mult(master, t, data->tsfiles, data->varids, data->ntsfiles, 
			  master->precip, 1.0 / (1000.0 * 86400.0));
    event->next_event += data->dt;
  }
  return event->next_event;
}


void precip_cleanup(sched_event_t *event, double t)
{
  precip_mdata_t *data = (precip_mdata_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}
