/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/swr.c
 *  
 *  Description:
 *  Event scheduler routines for Atmospheric swr.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: swr.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

/* The swr timeseries is used to calculate the surface heat flux */
/* at every time-step */
typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} swr_data_t;

/* The light timeseries is assumed to be a daily mean swr used in */
/* the ecology module. */
typedef struct {
  double dt;                    /* Time step */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
} light_data_t;

/* Functions for reading the schedule the swr forcings. */
int swr_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  swr_data_t *data = NULL;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  data = (swr_data_t *)malloc(sizeof(swr_data_t));
  schedSetPrivateData(event, data);
  data->ts = frc_read_cell_ts(master, params->swr, params->swr_dt,
                              "swr", "W m-2", &data->dt, &data->id,
                              &master->swr);

  if (data->ts == NULL) {
    free(data);
    if (master->heatflux & NET_HEAT && master->swr_attn)
      hd_quit("Attenuation of swr requires RADIATION file.\n");
    return 0;
  }

  if (master->albedo < -1.0 || master->albedo > 1.0)
    hd_quit("RADIATION requires ALBEDO parameter.\n");
  /*prm_read_double(master->prmfd, "SWR_ATTENUATION", &master->swr_attn);*/

  return 1;
}

double swr_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  swr_data_t *data = (swr_data_t *)schedGetPrivateData(event);
  geometry_t *geom = master->geom;
  int c, cc;

  if (t >= (event->next_event - SEPS)) {
    frc_ts_eval_grid(master, t, data->ts, data->id, master->swr, 1.0);
    event->next_event += data->dt;
  }

  /* Adjust for albedo */
  if (master->albedo >= 0.0) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->swr[c] *= (1.0 - master->albedo);
    }
  }

  return event->next_event;
}

void swr_cleanup(sched_event_t *event, double t)
{
  swr_data_t *data = (swr_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    if (master->swr)
      free(master->swr);
    hd_ts_free(master, data->ts);
    free(data);
  }
}

/* Functions for reading the schedule the light forcings. */
int light_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  light_data_t *data = NULL;
  FILE *fp = NULL;
  int index = -1;

  /* Read parameters */
  prm_set_errfn(hd_silent_warn);

  if ( (params->light == NULL) || (strlen(params->light) == 0) )
    return 0;

  /* Light is read from 'swr' in a transport file in trans_reset_event(); */
  /* no scheduled light function is required.                             */
  if (strcmp(params->light, "file") == 0) {
    master->light = d_alloc_1d(master->geom->sgsizS);
    params->albedo_l = master->albedo_l = 0.0;
    return 0;
  }

  data = (light_data_t *)malloc(sizeof(light_data_t));
  schedSetPrivateData(event, data);
  
  /* First check to see if its a file name */
  if ( (fp = fopen(params->light, "r")) != NULL ) {
    fclose(fp);
    data->ts = frc_read_cell_ts(master, params->light, params->light_dt,
				"swr", "W m-2", &data->dt, &data->id,
				&master->light);

  } else if ( (index = tracer_find_index(params->light, master->ntrS,
					 master->trinfo_2d)) > -1 ) {
    /* See if its a tracer */
    data->id = index;
    data->ts = NULL;
    data->dt = params->light_dt;
    master->light = d_alloc_1d(master->geom->sgsizS);

  } else {
    /* nothing found */
    hd_warn("LIGHT parameter '%s' is neither a file nor a 2D tracer in the model\n");
    free(data);
    return 0;
  }

  if (master->albedo_l < 0.0)
    hd_quit("LIGHT requires ALBEDO_LIGHT parameter.\n");

  return 1;
}

double light_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  light_data_t *data = (light_data_t *)schedGetPrivateData(event);
  int c;
  double conv = 1.0;

  if (t >= (event->next_event - SEPS)) {
    if (data->ts != NULL)
      /* Read from file */
      frc_ts_eval_grid(master, t, data->ts, data->id, master->light, conv);
    else {
      /* wire up tracer */
      for (c = 1; c <= geom->sgnumS; ++c) {
	master->light[c] = conv * master->tr_wcS[data->id][c];
      }
    }
    event->next_event += data->dt;
  }

  return event->next_event;
}

void light_cleanup(sched_event_t *event, double t)
{
  light_data_t *data = (light_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    free(master->light);
    if (data->ts != NULL)
      hd_ts_free(master, data->ts);
    free(data);
  }
}
