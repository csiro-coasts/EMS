/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/evap.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric evap.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: evap.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} evap_data_t;

/* Functions for reading the schedule the evap forcings. */
int evap_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  evap_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (evap_data_t *)malloc(sizeof(evap_data_t));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->evap, params->evap_dt,
                              "evaporation", "mm day-1", &data->dt,
                              &data->id, &master->evap);

  if (data->ts == NULL) {
    free(data);
    return 0;
  }

  return 1;
}

double evap_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  evap_data_t *data = (evap_data_t *)schedGetPrivateData(event);

  if (t >= (event->next_event - SEPS)) {
    /* Evaporation, convert to m/s */
    frc_ts_eval_grid(master, t, data->ts, data->id, master->evap,
                     1.0 / (1000.0 * 86400.0));
    event->next_event += data->dt;
  }

  return event->next_event;
}

void evap_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  evap_data_t *data = (evap_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->evap);
    hd_ts_free(master, data->ts);
    free(data);
  }
}
