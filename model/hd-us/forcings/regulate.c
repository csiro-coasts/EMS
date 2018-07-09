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
 *  $Id: regulate.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include <libgen.h>
#include "hd.h"

typedef struct {
  master_t *master;             /* Grid associated with             */
  double dt;                    /* Time increment to read file      */
  FILE *fp;                     /* Regulate file pointer            */
  char filename[MAXSTRLEN];     /* Regulate file name               */
  char cmnd[MAXSTRLEN];         /* Previous command                 */
  char prmname[MAXSTRLEN];      /* Parameter filename               */
} regulate_t;

void reset_dt(parameters_t *params, master_t *master, FILE *fp);

/*------------------------------------------------------------------*/
/* User run regulation initialisation                               */
/*------------------------------------------------------------------*/
int regulate_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  regulate_t *regulate = NULL;

  if (params->regulate_dt == 0)
    return 0;

  /*----------------------------------------------------------------*/
  /* Allocate memory                                                */
  regulate = (regulate_t *)malloc(sizeof(regulate_t));
  schedSetPrivateData(event, regulate);
  memset(regulate, 0, sizeof(regulate_t));
  regulate->master = master;
  regulate->dt = params->regulate_dt;
  strcpy(regulate->filename, params->regulate);
  sprintf(regulate->cmnd, "%c", '\0');
  strcpy(regulate->prmname, params->prmname);
  return 1;
}

/* END regulate_init()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* User run regulation events.                                      */
/* The regulate syntax is:                                          */
/* REGULATE CODE @ x units                                          */
/* x units is a time, e.g. '10 days'. Regulate action is only taken */
/* after the model time >= x units.                                 */
/* CODE is the regulation code:                                     */
/* 'STOP' terminates the run                                        */
/* 'PAUSE' suspends the run until 'RESUME' is read                  */
/* 'RESUME' continues a paused run                                  */ 
/* 'REINIT_OBC <obc_name> COND' reinitialises an OBC to a standard  */
/*     OBC of type COND if it exists, if not switch between the     */
/*     standard conditions 0 and 1.                                 */
/* 'DT_REINIT' re-initialises the timestep.                         */
/* 'HVISC_REINIT' re-initialises the horizontal viscosity.          */
/* 'DUMP_REINIT' re-initialises the output netCDF files.            */
/* 'TS_REINIT' re-initialises the ascii time-series files.          */
/* 'WIN_REINIT' re-initialises the window partitioning.             */
/* 'PSS_REINIT' re-initialises the source/sinks.                    */
/*------------------------------------------------------------------*/
double regulate_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  regulate_t *regulate = (regulate_t *)schedGetPrivateData(event);
  geometry_t *geom = master->geom;
  FILE *fp;
  char reg[MAXSTRLEN], buf[MAXSTRLEN], cmnd[MAXSTRLEN];
  char *files[MAXSTRLEN * MAXNUMARGS];
  double time;
  int m, n, i;

  if (t >= (event->next_event - SEPS)) {
    master->regf = 0;
    if (regulate->fp = fopen(regulate->filename, "r+")) {
      prm_set_errfn(errfn_quiet);
      if (prm_read_char(regulate->fp, "REGULATE", reg)) {

	strcpy(cmnd, reg);
	m = parseline(reg, files, MAXNUMARGS);

	/*----------------------------------------------------------*/
	/* Return if the regulation time is > than model time       */
	time = master->t;
	for (n = 1; n < m; n++) {
	  if (strcmp(files[n], "@") == 0 && n+2 < m) {
	    sprintf(buf, "%s %s", files[n+1], files[n+2]);
	    tm_scale_to_secs(buf, &time);
	    if (time > master->t) {
	      fclose(regulate->fp);
	      event->next_event += regulate->dt;
	      return event->next_event;
	    }
	  }
	}

	if ((fp = fopen(regulate->prmname, "r")) == NULL) {
	  hd_warn("regulate: Can't open file %s\n", regulate->prmname);
	  event->next_event += regulate->dt;
	  return event->next_event;
	}

	/*----------------------------------------------------------*/
	/* Stop the run                                             */
	if (strcmp(files[0], "STOP") == 0) {
	  hd_quit_and_dump("regulate: User stopped run at %f days\n",
			   master->days);
	/*----------------------------------------------------------*/
	/* Pause and resume                                         */
	} else if (strcmp(files[0], "PAUSE") == 0) {
	  hd_warn("regulate: User paused run at %5.2f days\n", master->days);
	  while(strcmp(files[0], "RESUME") != 0) {
	    fclose(regulate->fp);
	    if (regulate->fp = fopen(regulate->filename, "r+")) {
	      if (prm_read_char(regulate->fp, "REGULATE", buf))
		m = parseline(buf, files, MAXNUMARGS);
	    }
	  }
	  hd_warn("regulate: User resumed run at %5.2f days\n", master->days);
	/*----------------------------------------------------------*/
	/* Reconfigure an OBC                                       */
	} else if (strcmp(files[0], "OBC_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  i = 1;
	  for (n = 0; n < geom->nobc; n++) {
	    open_bdrys_t *open = geom->open[n];
	    if (m >= 2 && strcmp(files[1], open->name) == 0) {
	      i = 0;
	      break;
	    }
	  }
	  if (i) {
	    hd_warn("regulate: Can't find OBC %s to reconfigure\n", files[1]); 
	    event->next_event += regulate->dt;
	    return event->next_event;
	  }
	  master->regf |= RS_OBCSET;
	  sprintf(master->runerror, "OBC %d", n);
	  if (m >= 3)
	    sprintf(master->runerror, "%s %s", master->runerror, files[2]);
	  hd_warn("regulate: User requested OBC %s reconfigure\n",
		  files[1]);
	/*----------------------------------------------------------*/
	/* Reconfigure dumpfiles                                    */
	} else if (strcmp(files[0], "DUMP_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  hd_warn("regulate: User requested dumpfile reconfiguration at %.2f days\n", 
		  master->days);
	  dumpfile_resetup(master);
	/*----------------------------------------------------------*/
	/* Reconfigure time-series files                            */
	} else if (strcmp(files[0], "TS_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  hd_warn("regulate: User requested timeseries reconfiguration at %.2f days\n", 
		  master->days);
	  ts_resetup(master, fp);
	  master->regf |= RS_TSSET;
	/*----------------------------------------------------------*/
	/* Reconfigure point source/sinks                           */
	} else if (strcmp(files[0], "PSS_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  hd_warn("regulate: User requested point source/sink reconfiguration at %.2f days\n", 
		  master->days);
	  master->regf |= RS_PSSSET;
	/*----------------------------------------------------------*/
	/* Reconfigure time-step                                    */
	} else if (strcmp(files[0], "DT_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  reset_dt(params, master, fp);
	  reset_hor_diff(master, 0.0, 0.0, master->diff_scale);
	  master->regf |= RS_RESET;
	  hd_warn("regulate: User requested timestep = %.2f at %.2f days\n", 
		  master->grid_dt, master->days);
	/*----------------------------------------------------------*/
	/* Reconfigure horizontal mixing                            */
	} else if (strcmp(files[0], "HVISC_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  read_hdiff(params, fp, 0);
	  reset_hor_diff(master, params->u1vh, params->u2vh, params->diff_scale);
	  master->regf |= RS_RESET;
	  hd_warn("regulate: User requested horizontal mixing reconfiguration at %.2f days\n", 
		  master->days);
	/*----------------------------------------------------------*/
	/* Reconfigure windows                                      */
	} else if (strcmp(files[0], "WIN_REINIT") == 0 &&
		   strcmp(cmnd, regulate->cmnd) != 0) {
	  parameters_t *params = master->params;
	  if (params->nwindows > 1) {
	    read_window_info(params, fp);
	    if (geom->win_size) d_free_1d(geom->win_size);
	    if (params->win_size) {
	      geom->win_size = d_alloc_1d(geom->nwindows + 1);
	      for (n = 1; n <= geom->nwindows; n++) {
		geom->win_size[n] = params->win_size[n - 1];
	      }
	    }
	    if (params->nwn) {
	      geom->nwn = params->nwn;
	      geom->wnx = params->wnx;  
	      geom->wny = params->wny;
	    }
	    if (m > 4)
	      master->win_reset = atoi(files[1]);
	    else
	      master->win_reset = -1;
	    hd_warn("regulate: User requested %d windows at %.2f days\n", params->nwindows,
		    master->days);
	  } else {
	    hd_warn("regulate: WINDOWS must be initially > 1 for window resets.\n");
	  }
	}
	fclose(fp);
      }
      strcpy(regulate->cmnd, cmnd);
      fclose(regulate->fp);
    }
    event->next_event += regulate->dt;
  }
  return event->next_event;
}

/* END regulate_event()                                             */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Clean up the regulate recovery                                   */
/*------------------------------------------------------------------*/
void regulate_cleanup(sched_event_t *event, double t)
{
  regulate_t *regulate = (regulate_t *)schedGetPrivateData(event);
  master_t *master = regulate->master;

  if (regulate != NULL) {
    free(regulate);
  }
}

/* END regulate_cleanup()                                           */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Re-read the time-step                                            */
/*------------------------------------------------------------------*/
void reset_dt(parameters_t *params, master_t *master, FILE *fp)
{
  prm_get_time_in_secs(fp, "DT", &params->grid_dt);
  prm_read_int(fp, "IRATIO", &params->iratio);
  master->grid_dt = params->grid_dt;
  master->iratio = params->iratio;
}

/* END reset_dt()                                                   */
/*------------------------------------------------------------------*/
