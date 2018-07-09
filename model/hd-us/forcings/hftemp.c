/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/hftemp.c
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
 *  $Id: hftemp.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct hftemp_data {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} hftemp_data;

/* Functions for reading the schedule the hftemp forcings. */
int hftemp_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  hftemp_data *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (hftemp_data *) malloc(sizeof(hftemp_data));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->hftemp, params->hftemp_dt,
                              "heatflux_temp", "Degrees C", &data->dt,
                              &data->id, &master->hftemp);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double hftemp_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  hftemp_data *data = (hftemp_data *) schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid(master, t, data->ts, data->id, master->hftemp, 1.0);
    event->next_event += data->dt;
  }

  return event->next_event;
}

void hftemp_cleanup(sched_event_t *event, double t)
{
  hftemp_data *data = (hftemp_data *) schedGetPrivateData(event);

  if (data != NULL) {
    free(master->hftemp);
    hd_ts_free(master, data->ts);
    free(data);
  }
}
