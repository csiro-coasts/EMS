/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/wetbulb.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric wetb.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: wetbulb.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} wetbulb_data_t;

typedef struct {
  double dt;                    /* Time step */
  int ntsfiles;                 /* Number of timeseries files */
  timeseries_t **tsfiles;       /* Timeseries files */
  int varids[MAXNUMTSFILES];    /* Variable id's */
} wetbulb_mdata_t;

/* Functions for reading the schedule the rh forcings. */
int wetbulb_init_single(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  wetbulb_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (wetbulb_data_t *)malloc(sizeof(wetbulb_data_t));
  schedSetPrivateData(event, data);
  if (data->ts = frcw_read_cell_ts(master, params->wetb, params->wetb_dt,
                              "wet_bulb", "Degrees C", &data->dt,
                              &data->id, &master->wetb))
    master->sh_f = WETBULB;
  else if (data->ts = frc_read_cell_ts(master, params->wetb, params->wetb_dt,
                              "dew_point", "degrees C", &data->dt,
                              &data->id, &master->wetb))
    master->sh_f = DEWPOINT;
  else
    master->sh_f = NONE;

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double wetbulb_event_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wetbulb_data_t *data = (wetbulb_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid(master, t, data->ts, data->id, master->wetb, 1.0);
    event->next_event += data->dt;
  }

  return event->next_event;
}

void wetbulb_cleanup_single(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wetbulb_data_t *data = (wetbulb_data_t *)schedGetPrivateData(event);
  if (data != NULL) {
    free(master->wetb);
    hd_ts_free(master, data->ts);
    free(data);
  }

}


int wetbulb_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  wetbulb_mdata_t *data = NULL;

  char wetb[MAXSTRLEN];



  /* Read parameters */
  prm_set_errfn(hd_silent_warn);
  data = (wetbulb_mdata_t *)malloc(sizeof(wetbulb_mdata_t));
  schedSetPrivateData(event, data);

  strcpy(wetb, params->wetb);
  if (data->tsfiles = frc_read_cell_ts_mult(master, wetb, params->wetb_dt,
					    "wet_bulb", "Degrees C", &data->dt, data->varids, 
					    &data->ntsfiles, &master->wetb, 0))
    master->sh_f = WETBULB;
  else if (data->tsfiles = frc_read_cell_ts_mult(master, params->wetb, params->wetb_dt,
						 "dew_point", "degrees C", &data->dt, data->varids, 
						 &data->ntsfiles, &master->wetb, 1))
    master->sh_f = DEWPOINT;
  else
    master->sh_f = NONE;

  if (data->tsfiles == NULL) {
    free(data);
    return 0;
  }
  return 1;
}


double wetbulb_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wetbulb_mdata_t *data = (wetbulb_mdata_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid_mult(master, t, data->tsfiles, data->varids, data->ntsfiles, 
			  master->wetb, 1.0);
    event->next_event += data->dt;
  }
  return event->next_event;
}


void wetbulb_cleanup(sched_event_t *event, double t)
{
  wetbulb_mdata_t *data = (wetbulb_mdata_t *)schedGetPrivateData(event);

  if (data != NULL) {
    int i;
    for (i=0; i<data->ntsfiles; ++i)
       hd_ts_free(master, data->tsfiles[i]);
    free(data);
  }
}
