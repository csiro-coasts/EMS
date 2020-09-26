/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/dhw.c
 *  
 *  Description:
 *  Routines relating to Degree Heating Weeks calculations
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dhw.c 6423 2019-11-22 00:23:33Z her127 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/* Hotspot increment: Hotspot = MMM + thresh                         */
double thresh = 1.0;

/* Local functions */
static int dhw_init(sched_event_t *event);
double dhw_event(sched_event_t *event, double t);
static void dhw_cleanup(sched_event_t *event, double t);

typedef struct {
  master_t *master;             /* Grid associated with */
  int ntsfiles;                 /* Number of time-series files */
  timeseries_t **tsfiles;       /* Array of time-series files */
  cstring *tsnames;             /* Array of time-series files */
  double dt;                    /* Relaxation time step */
  char dhdo[MAXSTRLEN];         /* Offset dhd file */
  double offset;                /* Offset dhd increment (sec) */
  int mtype;                    /* Type of daily mean */
  int checked;
  int first;
} tr_dhw_data_t;


/*-------------------------------------------------------------------*/
/* Opens files for tracer resets                                     */
/*-------------------------------------------------------------------*/
static void open_dhd_tsfiles(master_t *master, char *fnames,
                                  tr_dhw_data_t *dhw)
{
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int i;
  int nf = parseline(fnames, fields, MAXNUMARGS);

  dhw->ntsfiles = nf;

  if (nf > 0) {
    dhw->tsfiles = (timeseries_t **)malloc(sizeof(timeseries_t *) * nf);
    dhw->tsnames = (cstring *) malloc(sizeof(cstring) * nf);
    memset(dhw->tsfiles, 0, sizeof(timeseries_t *) * nf);
    memset(dhw->tsnames, 0, sizeof(cstring) * nf);
    for (i = 0; i < nf; ++i) {
      strcpy(dhw->tsnames[i], fields[i]);
      dhw->tsfiles[i] = hd_ts_read(master, dhw->tsnames[i], 0);
    }
  } else
    hd_quit("Must specify at least one dhw time-series file.\n");
}

/* END open_dhd_tsfiles()                                       */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void tracer_dhw_init(master_t *master)
{
  parameters_t *params = master->params;


  char *oset = "12 week";
  double dt;
 
  tr_dhw_data_t *dhw = NULL;
  

  if (!(master->dhwf & DHW_NOAA)) return;

  /* Allocate memory for dhw structure, populate and           */
  /* register the scheduler events.                            */
  dhw = (tr_dhw_data_t *)malloc(sizeof(tr_dhw_data_t));
  memset(dhw, 0, sizeof(tr_dhw_data_t));

  dhw->master = master;
  dhw->dt = params->dhw_dt;
  strcpy(dhw->dhdo, params->dhdf);
  tm_scale_to_secs(oset, &dt);
  dhw->offset = dt;
  dhw->checked = 0;
  dhw->first = 1;
  if (params->dhwf & DHW_INT)
    dhw->mtype = DHW_INT;
  else if (params->dhwf & DHW_MEAN)
    dhw->mtype = DHW_MEAN;
  else
    dhw->mtype = DHW_SNAP;

  /* Register the scheduled function */
  sched_register(schedule, "dhw", dhw_init, dhw_event, dhw_cleanup,
		 dhw, NULL, NULL);
}

/* END tracer_dhw_init()                                             */
/*-------------------------------------------------------------------*/

static int dhw_init(sched_event_t *event)
{
  return 1;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
double dhw_event(sched_event_t *event, double t)
{
  tr_dhw_data_t *dhw = (tr_dhw_data_t *)schedGetPublicData(event);
  master_t *master = dhw->master;
  int c, cs, cc, n;
  double dhdo, tin;
  double fact = 1.0 / 7.0;
  int varid;
  int found = 1;
  int nanf = 0;

  /* Only increment the dhw after the first day of iteration         */
  if(dhw->first) {
    dhw->first = 0;
    event->next_event += dhw->dt;
    return event->next_event;
  }

  if (t >= (event->next_event - SEPS)) {

    /* Open the dhd offset file                                      */
    open_dhd_tsfiles(master, dhw->dhdo, dhw);

    /* Check that the offset dhd data exists in file */
    for (n = 0; n < dhw->ntsfiles; ++n) {
      char buf[MAXSTRLEN];
      timeseries_t *ts = dhw->tsfiles[n];
      varid = ts_get_index(ts, fv_get_varname(dhw->tsnames[n], "dhd", buf));
      if (varid < 0) found = 0;
    }
    if (found && !dhw->checked) {
      prm_set_errfn(hd_silent_warn);
      hd_ts_multifile_check(dhw->ntsfiles, dhw->tsfiles, dhw->tsnames, "dhd",
					  schedule->start_time, schedule->stop_time);
      dhw->checked = 1;
    }

    /* Check if the offset dhd value is in the file                  */
    tin = t - dhw->offset;
    for (n = 0; n < dhw->ntsfiles; n++) {
      if (ts_has_time(dhw->tsfiles[n], tin) == 0)
	found = 0;
    }

    /* Reset the dhw value. This is the sum of dhd values in the dhw */
    /* time window (usually 12 weeks). Therefore, subtract the       */
    /* offset dhd value (e.g. value 12 weeks ago) and add the        */
    /* current value.                                                */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];

      /* Get the offset dhd value from file                          */      
      dhdo = 0.0;
      if (found)
	dhdo = hd_ts_multifile_eval_xyz_by_name(dhw->ntsfiles, dhw->tsfiles,
						dhw->tsnames, "dhd", tin,
						geom->cellx[cs], 
						geom->celly[cs],
						geom->cellz[c] * master->Ds[cs]);

      /* Update the dhw value                                        */
      if(isnan(dhdo)) {
	dhdo = 0.0;
	nanf = 1;
      }

      /* DHD is the integral of temp-MMM over the day                */
      if (dhw->mtype & (DHW_INT|DHW_SNAP))
	master->dhw[c] += fact * (master->dhd[c] - dhdo);
      /* DHD is the daily mean temp - MMM                            */
      if (dhw->mtype & DHW_MEAN) {
	if (master->dhd[c] > master->dhwc[c] + thresh)
	  master->dhw[c] += fact * ((master->dhd[c] - master->dhwc[c]) - dhdo);
      }

      /* Reinitialize the dhd value                                  */
      master->dhd[c] = 0.0;
    }

    event->next_event += dhw->dt;
  }
  if(nanf)
    hd_warn("dhw_event: NaN found in file %s at %f days\n", dhw->dhdo, master->days);
    
  /* Close the dhd file */
  for (cc = 0; cc < dhw->ntsfiles; ++cc)
    hd_ts_free(master, dhw->tsfiles[cc]);
  free((cstring *)dhw->tsnames);

  return event->next_event;
}

/* END dhw_event()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void dhw_cleanup(sched_event_t *event, double t)
{
  tr_dhw_data_t *dhw = (tr_dhw_data_t *)schedGetPrivateData(event);

  if (dhw != NULL) {
    int i;
    for (i=0; i<dhw->ntsfiles; ++i)
       hd_ts_free(master, dhw->tsfiles[i]);
    free(dhw);
  }
}

/* END dhw_cleanup()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the dhd value.                                         */
/*-------------------------------------------------------------------*/
void calc_dhd(geometry_t *window,       /* Window geometry       */
	      window_t *windat,         /* Window data           */
	      win_priv_t *wincon        /* Window constants      */
	      )
{
  int c, cc;
  double fact = 1.0 / 86400.0;

  if (wincon->dhwf & DHW_NOAA) {
    if (windat->dhd && windat->dhwc) {
      /* DHD is the integral of temp-MMM over the day                */
      if (wincon->dhwf & DHW_INT) {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (windat->temp[c] > windat->dhwc[c] + thresh) {
	    windat->dhd[c] += fact * (windat->temp[c] - windat->dhwc[c]) * windat->dt;
	  }
	}
      }
      /* DHD is the daily mean temp                                  */
      if (wincon->dhwf & DHW_MEAN) {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  windat->dhd[c] += fact * windat->temp[c] * windat->dt;
	}
      }
      /* DHD is the temp at dhwh hours                               */
      if (wincon->dhwf & DHW_SNAP) {
	double hrs = (double)((int)windat->days + fabs(wincon->dhwh));
	if (windat->days < hrs && wincon->dhwf & DHW_SET)
	  wincon->dhwf &= ~DHW_SET;
	if (windat->days >= hrs && !(wincon->dhwf & DHW_SET)) {
	  wincon->dhwf |= DHW_SET;
	  for (cc = 1; cc <= window->b3_t; cc++) {
	    c = window->w3_t[cc];
	    if (windat->temp[c] > windat->dhwc[c] + thresh) {
	      windat->dhd[c] = windat->temp[c] - windat->dhwc[c];
	    }
	  }
	}
      }
    }
  }
}

/* END calc_dhd()                                                    */
/*-------------------------------------------------------------------*/
