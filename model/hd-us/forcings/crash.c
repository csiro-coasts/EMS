/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/da.c
 *  
 *  Description:
 *  Event scheduler routines for automated crash recovery.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: crash.c 6264 2019-08-08 04:23:54Z her127 $
 *
 */

#include <math.h>
#include <string.h>
#include <libgen.h>
#include "hd.h"

#define MAXCRASH 5

/* Local functions */
static int crash_init(sched_event_t *event);
static void crash_cleanup(sched_event_t *event, double t);

typedef struct {
  hd_data_t *hd_data;           /* Grid associated with             */
  double dt;                    /* Restart time step                */
  double stop;                  /* Stop time                        */
  double tnext;                 /* Next event time                  */
  char rsfname[MAXSTRLEN];      /* Restart file name                */
  int flag;                     /* Control flag                     */
  int cm;                       /* Timestep reduction method        */
  double odt;                   /* Original timestep                */
  double ou1vh;                 /* Original u1vh                    */ 
  double ou2vh;                 /* Original u2vh                    */ 
  double fact;                  /* Timestep reduction factor        */
  int ncr;                      /* Crash recovery iteration         */
} crash_data_t;


/*------------------------------------------------------------------*/
/* Automated crash recovery initialisation                          */
/*------------------------------------------------------------------*/
void crash_recovery_init(hd_data_t *hd_data)
{
  master_t *master = hd_data->master;
  parameters_t *params = hd_data->params;
  crash_data_t *crash = NULL;

  /*----------------------------------------------------------------*/
  /* Allocate memory                                                */
  crash = (crash_data_t *)malloc(sizeof(crash_data_t));
  memset(crash, 0, sizeof(crash_data_t));
  crash->hd_data = hd_data;

  /*----------------------------------------------------------------*/
  /* Check restarts are set.                                        */
  if (strlen(params->restart_name) == 0 || master->restart_dt == 0.0)
    hd_quit("crash_init: Restart files must be written for automatic crash restarts.\n");
  else {
    strcpy(crash->rsfname, master->restart_name);
    crash->dt = master->restart_dt;
  }
  if (master->u1vh0 > 0.0) crash->ou1vh = master->u1vh0;
  if (master->u1vh0 > 0.0) crash->ou2vh = master->u2vh0;
  crash->cm = RS_PREV;
  master->crf = crash->flag = NONE;
  crash->odt = master->grid_dt;
  crash->fact = 5.0;
  if (crash->cm & RS_PREV) crash->fact = 2.0;
  crash->ncr = 0;
  tm_scale_to_secs(params->stop_time, &crash->stop);
  crash->tnext = crash->stop;

  sched_register(schedule, "crash",
		 crash_init, crash_event, crash_cleanup,
		 crash, NULL, NULL);
}

/* END crash_revovery_init()                                        */
/*------------------------------------------------------------------*/

static int crash_init(sched_event_t *event)
{
  return 1;
}


/*------------------------------------------------------------------*/
/* Automated crash recovery                                         */
/* crash_event(sched_get_even_by_name(schedule, "crash"),master->t) */
/*------------------------------------------------------------------*/
double crash_event(sched_event_t *event, double t)
{
  crash_data_t *crash = (crash_data_t *)schedGetPublicData(event);
  hd_data_t *hd_data = crash->hd_data;
  master_t *master = hd_data->master;

  dump_data_t *dumpdata = master->dumpdata;
  char restart_fname[MAXSTRLEN];
  int fid = 0;
  double newt;
  int i;

  /* Crash has occurred : read the restart and reset master constants */
  if (master->crf == RS_RESTART) {

    strcpy(restart_fname, crash->rsfname);
    /* Open the file                                                  */
    if (nc_open(restart_fname, NC_NOWRITE, &fid) != NC_NOERR)
      hd_quit("Can't find crash restart file %s\n", restart_fname);
        
    /* Reset the master state                                         */
    dump_re_read(master, fid, 0);

    /* Get the new start time                                         */
    schedule->t = newt = get_restart_time(restart_fname, schedule->units);
    /* Reset the dumpfile dump times, except the restart file         */
    for (i = 0; i < dumpdata->ndf - 1; ++i) {
      dumpdata->dumplist[i].reset(dumpdata, &dumpdata->dumplist[i], newt);
    }
    /* Reset timeseries files (not required)                          */
    for (i = 0; i < nts; ++i) {
      /* Rewind to the first next time after newt */
      while (newt <= (tslist[i].tsout - tslist[i].tsdt))
	tslist[i].tsout -= tslist[i].tsdt;
      fflush(tslist[i].fp);
      fclose(tslist[i].fp);
    }
    forced_restart = 1;
    master->t = newt;
    ts_init(sched_get_even_by_name(schedule, "timeseries"));
    forced_restart = 0;
    tm_change_time_units(master->timeunit, master->output_tunit, &newt, 1);

    /* Reset the timestep and friction                                */
    master->grid_dt /= crash->fact;
    if (crash->ncr > MAXCRASH) {
      crash_restart = 0;
      hd_quit("Can't recover from instability: %d attempts at %5.1f days\n", crash->ncr, master->days);
    }
    reset_hor_diff(master, 0.0, master->diff_scale);
    reset_obc_adjust(master->geom, master->grid_dt / master->iratio);
    crash->ncr++;

    /* Set the time when normal simulation commences                  */
    master->crf = crash->flag = (RS_RESUME|RS_WINSET);
    /*event->next_event = schedule->t + crash->dt;*/
    crash->tnext = schedule->t + crash->dt;

    hd_warn("CRASHED: Recovery #%d; restarting at %5.1f days, timestep = %5.2f, resume = %5.1f days\n", 
	    crash->ncr, newt, master->grid_dt, crash->tnext / 86400.0);

    /* Reset the windows                                              */
    for (i = 1; i <=master->nwindows; i++)
      window_reset(master, hd_data->window[i], hd_data->windat[i], hd_data->wincon[i], RS_ALL);

  } else if (master->crf == RS_RESUME && t >= (crash->tnext - SEPS)) {
    newt = t;
    tm_change_time_units(master->timeunit, master->output_tunit, &newt, 1);

    /* Reset the timestep and friction                                */
    if (crash->cm == RS_ORIG)
      master->grid_dt = crash->odt;
    else {
      master->grid_dt = min(crash->odt,  master->grid_dt * crash->fact);
      crash->tnext = schedule->t + crash->dt;
      master->crf = crash->flag = (RS_RESUME|RS_WINSET);
    }
    reset_hor_diff(master, master->u1vh0, master->diff_scale);
    reset_obc_adjust(master->geom, master->grid_dt / master->iratio);
    if (master->grid_dt == crash->odt) {
      master->crf = crash->flag = (NONE|RS_RESET);
      crash->ncr = 0;
      crash->tnext = crash->stop;
    }
    for (i = 1; i <=master->nwindows; i++)
      window_reset(master, hd_data->window[i], hd_data->windat[i], hd_data->wincon[i], RS_VH);

    if (crash->cm == RS_ORIG)
      hd_warn("RECOVERY: Crash recovery %d complete : commencing at %5.1f days\n", crash->ncr, newt);
    else
      hd_warn("RECOVERY: Crash recovery %d complete : commencing at %5.1f days with dt=%f\n", crash->ncr, newt, master->grid_dt);
  } else {
    /*hd_warn("crash: no action @ %f, next = %f\n",master->days,crash->tnext/86400.0);*/
  }
  return crash->tnext;
}

/* END crash_event()                                                */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Clean up the crash recovery                                      */
/*------------------------------------------------------------------*/
void crash_cleanup(sched_event_t *event, double t)
{
  crash_data_t *crash = (crash_data_t *)schedGetPublicData(event);

  if (crash != NULL) {
    free(crash);
  }
}

/* END crash_cleanup()                                              */
/*------------------------------------------------------------------*/
