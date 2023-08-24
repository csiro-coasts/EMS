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
 *  $Id: crash.c 6982 2022-02-27 23:39:50Z her127 $
 *
 */

#include <math.h>
#include <string.h>
#include <libgen.h>
#include <sys/time.h>
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
  int type;                     /* Recovery type                    */
  double odt;                   /* Original timestep                */
  double ou1vh;                 /* Original u1vh                    */ 
  double ou2vh;                 /* Original u2vh                    */ 
  double fact;                  /* Timestep reduction factor        */
  double vhi;                   /* Viscosity increment              */
  double vhl;                   /* Viscosity upper limit            */
  double *afact;                /* Array of timestep factors        */
  int nfact;                    /* Size of afact                    */
  int nc;                       /* afact counter                    */
  int ncr;                      /* Crash recovery iteration         */
  char cname[MAXSTRLEN];        /* Crash log name                   */
} crash_data_t;

int reset_vhreg(master_t *master, geometry_t **window, 
		window_t **windat, crash_data_t *crash, int mode, double *vh);
void crash_log(parameters_t *params, crash_data_t *crash);
int get_afact(crash_data_t *crash, char *key);

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

  /* Set the recovery mode                                          */
  crash->type = RS_REC;
  if (strlen(params->crf)) {
    int i,n;
    char *fields[MAXSTRLEN * MAXNUMARGS];
    char key[MAXSTRLEN];
    n = parseline(params->crf, fields, MAXNUMARGS);
    /* Set defaults                                                 */
    for (i = 0; i < n; i++) {
      if (strcmp(fields[i], "OPTIMIZE") == 0) {
	crash->type = RS_OPT;
	crash->flag = RS_VHSET;
	crash->fact = 1.25;
	crash->vhi = 10.0;
	crash->vhl = 100.0;
      }
      if (strcmp(fields[i], "RECOVER") == 0) {
	crash->type = RS_REC;
	crash->fact = 2.0;
      }
    }
    /* Set options                                                  */
    crash->afact = NULL;
    crash->nfact = crash->nc = 0;
    for (i = 0; i < n; i++) {
      if (decode_tag(fields[i], "DT", key)) {
	if (!get_afact(crash, key))
	  crash->fact = 1.0 / atof(key);
      }
      if (decode_tag(fields[i], "VH", key))
	crash->vhi = atof(key);
      if (decode_tag(fields[i], "VHL", key))
	crash->vhl = atof(key);
    }
  }

  crash_log(params, crash);

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
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  char restart_fname[MAXSTRLEN], buf[MAXSTRLEN];
  int errf, fid = 0;
  double newt, vhn;
  int c, cc, c2, i, ren;
  long tm;
  FILE *fp, *sp, *dp;

  /* Crash has occurred : read the restart and reset master constants */
  if (master->crf == RS_RESTART) {
    fp = fopen(crash->cname, "a");
    strcpy(restart_fname, crash->rsfname);
    /* Open the file                                                  */
    if ((errf = nc_open(restart_fname, NC_NOWRITE, &fid)) != NC_NOERR)
      hd_quit("Can't find crash restart file %s (errf=%d)\n", restart_fname, errf);
        
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

    /* Reset the timestep and friction for recovery                   */
    if (crash->type & RS_REC) {
      master->grid_dt /= crash->fact;
      if (crash->ncr > MAXCRASH) {
	crash_restart = 0;
	fprintf(fp, "Can't recover from instability: %d attempts at %5.1f days\n", crash->ncr, master->days);
	hd_quit("Can't recover from instability: %d attempts at %5.1f days\n", crash->ncr, master->days);
      }
      set_hor_diff(master);
      /*reset_hor_diff(master, 0.0, master->diff_scale);*/
      reset_obc_adjust(master->geom, master->grid_dt / master->iratio);
      crash->ncr++;

      /* Set the time when normal simulation commences                */
      master->crf = crash->flag = (RS_RESUME|RS_WINSET);
      /*event->next_event = schedule->t + crash->dt;*/
      crash->tnext = schedule->t + crash->dt;

      fprintf(fp, "CRASHED: Recovery #%d; restarting at %5.1f days, timestep = %5.2f, resume = %5.1f days\n", 
	    crash->ncr, newt, master->grid_dt, crash->tnext / 86400.0);
      hd_warn("CRASHED: Recovery #%d; restarting at %5.1f days, timestep = %5.2f, resume = %5.1f days\n", 
	    crash->ncr, newt, master->grid_dt, crash->tnext / 86400.0);

    }

    /* Optimise the timestep and friction                             */
    if (crash->type & RS_OPT) {
      /* Increase the friction (only use with diff_scale=AUTO)        */
      if (crash->flag & RS_VHSET) {
	if (!(master->diff_scale & VH_REG) && master->u1vh0 >= crash->vhl) {
	  crash->flag &= ~RS_VHSET;
	} else {
	  if (master->diff_scale & VH_REG)
	    ren = reset_vhreg(master, hd_data->window, hd_data->windat, crash, 0, &vhn);
	  else
	    master->u1vh0 = min(crash->vhl, master->u1vh0 + crash->vhi);
	  set_hor_diff(master);
	}
      } else {
	/* Decrease the timestep                                      */
	if (crash->nfact) {
	  crash->fact = crash->afact[crash->nc];
	  crash->nc++;
	  if (crash->nc == crash->nfact) crash->nc--;
	}
	master->grid_dt /= crash->fact;
	if (master->grid_dt < 1.0) {
	  fprintf(fp, "CRASH RECOVERY: dt < 1\n");
	  hd_quit("CRASH RECOVERY: dt < 1\n");
	}
	if (master->diff_scale & VH_REG)
	  ren = reset_vhreg(master, hd_data->window, hd_data->windat, crash, 1, &vhn);
	else
	  master->u1vh0 = crash->vhi;
	set_hor_diff(master);
	crash->flag |= RS_VHSET;
      }
      if (master->diff_scale & VH_REG) {
	hd_warn("CRASHED: Optimisation #%d; restarting at %5.2f days, dt=%5.1f U1VH=%5.1f%%; region%d\n", 
		crash->ncr, newt, master->grid_dt, vhn, ren);
	fprintf(fp, "CRASHED: Optimisation #%d; restarting at %5.2f days, dt=%5.1f\n",
		crash->ncr, newt, master->grid_dt          );
	fprintf(fp, "         region%d U1VH=%5.1f%%\n",  ren, vhn);
	fprintf(fp, "         U1VH %s\n", master->u1vhci);
      } else {
	hd_warn("CRASHED: Optimisation #%d; restarting at %5.2f days, dt=%5.1f U1VH=%5.1f%%\n", 
		crash->ncr, newt, master->grid_dt, master->u1vh0);
	fprintf(fp, "CRASHED: Optimisation #%d; restarting at %5.2f days, dt=%5.1f U1VH=%5.1f%%\n", 
		crash->ncr, newt, master->grid_dt, master->u1vh0);
      }
      schedule->exec_start_time = time(NULL);
      master->crf = RS_WINSET;
      crash->ncr++;
    }

    /* Reset the windows                                              */
    sp = fopen("cr.site", "a");
    for (i = 1; i <= master->nwindows; i++) {
      geometry_t *window = hd_data->window[i];
      window_t *windat = hd_data->windat[i];
      window_reset(master, hd_data->window[i], hd_data->windat[i], hd_data->wincon[i], RS_ALL);
      if ((c = windat->cc2)) {
	int cg = window->wsa[c];
	fprintf(sp, "%f %f c%d\n", geom->cellx[cg], geom->celly[cg], crash->ncr);
	fprintf(fp, "         window%d %d[%f %f], %5.2f days\n", i, c, 
		geom->cellx[cg], geom->celly[cg], windat->days);
      }
    }
    if ((dp = fopen("diag.txt", "r")) != NULL) {
      if (prm_read_char(dp, "Simulation time  =", buf)) {
	fprintf(fp, "         Crash time  = %s\n", buf);
      }
      if (prm_read_char(dp, "Total time ratio =", buf)) {
	fprintf(fp, "         Run time ratio = %s\n", buf);
      }
      fclose(dp);
    }
    time(&tm);
    fprintf(fp, "         Restart time: %s\n", ctime(&tm));
    fclose(fp);
    fclose(sp);

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
    if (crash->type & RS_REC) {
      reset_hor_diff(master, master->u1vh0, master->diff_scale);
      reset_obc_adjust(master->geom, master->grid_dt / master->iratio);
    }

    if (master->grid_dt == crash->odt) {
      master->crf = crash->flag = (NONE|RS_RESET);
      crash->ncr = 0;
      crash->tnext = crash->stop;
    }
    for (i = 1; i <=master->nwindows; i++)
      window_reset(master, hd_data->window[i], hd_data->windat[i], hd_data->wincon[i], RS_VH);

    fp = fopen(crash->cname, "a");
    if (crash->cm == RS_ORIG) {
      fprintf(fp, "RECOVERY: Crash recovery %d complete : commencing at %5.1f days\n", crash->ncr, newt);
      hd_warn("RECOVERY: Crash recovery %d complete : commencing at %5.1f days\n", crash->ncr, newt);
    } else {
      fprintf(fp, "RECOVERY: Crash recovery %d complete : commencing at %5.1f days with dt=%f\n", crash->ncr, newt, master->grid_dt);
      hd_warn("RECOVERY: Crash recovery %d complete : commencing at %5.1f days with dt=%f\n", crash->ncr, newt, master->grid_dt);
    }
    fclose(fp);
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


/*------------------------------------------------------------------*/
/* Routine to increase viscosity in a regional representation in    */
/* the region where instability occurred.                           */
/*------------------------------------------------------------------*/
int reset_vhreg(master_t *master, 
		geometry_t **window, 
		window_t **windat, 
		crash_data_t *crash,
		int mode,
		double *vhn
		)
{
  geometry_t *geom = master->geom;
  char buf[MAXSTRLEN];
  int nwindows = master->nwindows;
  int nw, n, nr, nf, c, cg, ret;
  int nreg, *rgn;
  double *rgv, sgn, smg;
  int *mask;
  char *fields[MAXSTRLEN * MAXNUMARGS];

  strcpy(buf, master->u1vhci);
  nf = parseline(buf, fields, MAXNUMARGS);
  sprintf(master->u1vhci, "%s %s", fields[0], fields[1]);
  nreg = nf - 2;

  /* Save the region viscosities                                    */
  rgn = i_alloc_1d(nreg);
  rgv = d_alloc_1d(nreg);
  mask = i_alloc_1d(nreg);
  for (n = 0; n < nreg; n++)
    sscanf(fields[n+2], "%d:%lf", &rgn[n], &rgv[n]);

  /* Loop through windows looking for instabilities                 */
  memset(mask, 0, nreg * sizeof(int));
  for (nw = 1; nw <= nwindows; nw++) {
    if ((c = windat[nw]->cc2)) {
      cg = window[nw]->wsa[c];
      /* Get the region where instability occurred                  */
      nr = master->vhreg[geom->m2d[cg]];
      /* Alter the unstable region viscosity                        */
      if (!mask[nr] && nr >= 0 && nr < nreg) {
	sgn = 1.0;
	if (rgv[nr] < 0.0) {
	  sgn = -1.0;
	  rgv[nr] = fabs(rgv[nr]);
	}
	smg = rgv[nr] - floor(rgv[nr]);
	rgv[nr] = floor(rgv[nr]);
	if (mode == 1)
	  rgv[nr] = crash->vhi;
	else {
	  rgv[nr] = min(crash->vhl, rgv[nr] + crash->vhi);
	  if (rgv[nr] >= crash->vhl) crash->flag &= ~RS_VHSET;
	}
	rgv[nr] = sgn * (rgv[nr] + smg);
	mask[nr] = 1;
	*vhn = rgv[nr];
	ret = nr;
      }
    }
  }
  /* Rebuild the regional input                                     */
  for (n = 0; n < nreg; n++) {
    sprintf(buf, " %d:%3.1f", rgn[n], rgv[n]);
    strcat(master->u1vhci, buf);
  }
  i_free_1d(mask);
  i_free_1d(rgn);
  d_free_1d(rgv);
  return(ret);
}

/* END reset_vhreg()                                                */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Opens a crash recovery logfile and writes a header.              */
/*------------------------------------------------------------------*/
void crash_log(parameters_t *params, crash_data_t *crash)
{
  FILE *fp;
  char buf[MAXSTRLEN], pname[MAXSTRLEN];
  int n, i, j;
  long t;

  /* Get the name of the parameter file                             */
  if (endswith(params->prmname, ".tran")) {
    j = 5;
    strcpy(buf, params->prmname);
  } else if (endswith(params->prmname, ".prm")) {
    j = 4;
    strcpy(buf, params->prmname);
  } else
    return;
  n = strlen(buf);
  for (i = 0; i < n-j; i++)
    pname[i] = buf[i];
  pname[i] = '\0';

  sprintf(buf, "%s.cr", pname);
  strcpy(crash->cname, buf);
  strcpy(params->crashname, buf);
  fp = fopen(buf, "w");

  fprintf(fp, "\nCRASH LOG for parameter file %s\n", params->prmname);
  fprintf(fp, "\n---------------------------\n");
  if (crash->type & RS_REC)
    fprintf(fp, "Operating in RECOVERY mode.\n");
  if (crash->type & RS_OPT) {
    fprintf(fp, "Operating in OPTIMIZE mode.\n");
    fprintf(fp, "Viscosity increment = %4.1f\n", crash->vhi);
    fprintf(fp, "Viscosity upper limit (%) = %4.1f\n", crash->vhl);
  }
  if (crash->nfact) {
    fprintf(fp, "Time-step reduction increment = ", 1.0 / crash->fact);
    for (i = 0; i < crash->nfact; i++)
      fprintf(fp, "%2.1f ", 1.0 / crash->afact[i]);
    fprintf(fp, "\n");
  } else
    fprintf(fp, "Time-step reduction increment = %2.1f\n", 1.0 / crash->fact);
  fprintf(fp, "---------------------------\n");
  time(&t);
  fprintf(fp, "Run:               %s", ctime(&t));
  fprintf(fp, "EMS Version:       %s\n", version);
  fprintf(fp, "Executable file:   %s\n",executable );
  getcwd(buf, MAXSTRLEN);
  fprintf(fp, "Working directory: %s\n",buf);
  fprintf(fp, "ID_CODE:           %s\n", params->runcode);
  fprintf(fp, "Input file:        %s\n", params->idumpname);
  if (strlen(params->opath)) fprintf(fp, "Output path:       %s\n", params->opath);
  fprintf(fp, "Parameter header:  %s\n", params->parameterheader);
  fprintf(fp, "Start time:        %s\n", params->start_time);
  fprintf(fp, "Stop time:         %s\n", params->stop_time);
  fprintf(fp, "\n");
  fflush(fp);
  fclose(fp);
  fp = fopen("cr.site", "w");
  /*fprintf(fp, "\n");*/
  fclose(fp);
}

/* END crash_log()                                                  */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Decodes a string, e.g. DTV:0.5,0.6,0.9 into an array             */
/*------------------------------------------------------------------*/
int get_afact(crash_data_t *crash, char *key)
{
  char buf[MAXSTRLEN];
  char *tok;
  int n = 0;

  strcpy(buf, key);
  tok = strtok(buf, ",");
  while (tok != NULL) {
    crash->nfact++;
    tok = strtok(NULL, ",");
  }
  if (crash->nfact) {
    crash->afact = d_alloc_1d(crash->nfact);
    strcpy(buf, key);
    tok = strtok(buf, ",");
    n = 0;
    while (tok != NULL) {
      crash->afact[n++] = 1.0 / atof(tok);
      tok = strtok(NULL, ",");
    }
    return(1);
  }
  return(0);
}

/* END get_afact()                                                  */
/*------------------------------------------------------------------*/
