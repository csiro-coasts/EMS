/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/tracers/relaxation.c
 *  
 *  Description:
 *  Tracer relaxation.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: relaxation.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/* Local functions */
static int tr_relax_init(sched_event_t *event);
static double tr_relax_event(sched_event_t *event, double t);
static void tr_relax_cleanup(sched_event_t *event, double t);

typedef struct {
  master_t *master;             /* Grid associated with */
  int tid;                      /* Tracer id */
  int ntsfiles;                 /* Number of time-series files */
  timeseries_t **tsfiles;       /* Array of time-series files */
  cstring *tsnames;             /* Array of time-series files */
  int varids[MAXNUMTSFILES];    /* Variable ids. */
  double dt;                    /* Relaxation time step */
  timeseries_t *tc;             /* Relaxation rate timeseries */
  int tcid;                     /* Relaxation rate id */
  double relax_rate;            /* Relaxation rate used */
  double dt0, dt1, rlx0, rlx1;  /* Adaptive relaxation endpoints */
  double slope;                 /* Adaptive relaxation slope */
  int tctype;                   /* Adaptive relaxation type */
  int da_started;               /* Start of a DA cycle */
} tr_relax_data_t;

static double get_relax_rate(master_t *master, tr_relax_data_t *relax, double dtr);

static void open_tr_relax_tsfiles(master_t *master, char *fnames,
                                  tr_relax_data_t *relax)
{
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int i;
  int nf = parseline(fnames, fields, MAXNUMARGS);

  relax->ntsfiles = nf;

  if (nf > 0) {
    relax->tsfiles = (timeseries_t **)malloc(sizeof(timeseries_t *) * nf);
    relax->tsnames = (cstring *) malloc(sizeof(cstring) * nf);
    memset(relax->tsfiles, 0, sizeof(timeseries_t *) * nf);
    memset(relax->tsnames, 0, sizeof(cstring) * nf);
    for (i = 0; i < nf; ++i) {
      strcpy(relax->tsnames[i], fields[i]);
      relax->tsfiles[i] = hd_ts_read(master, relax->tsnames[i], 0);
    }
  } else
    hd_quit("Must specify at least one relaxation time-series file.\n");
}


/* Check whether the tr_relax_data_t are to be included, if so, then
 * register them with the scheduler */
void tracer_relax_init(master_t *master)
{
  char buf[MAXSTRLEN];
  char tu0[MAXSTRLEN], tu1[MAXSTRLEN];
  char fname[MAXSTRLEN];
  double dt, tconst = 0.0;
  int t;
  tr_relax_data_t *relax = NULL;
  int rlx[master->ntr];

  master->nrlx = 0;
  prm_set_errfn(hd_silent_warn);
  for (t = master->atr; t < master->ntr; t++) {

    /* Read the 'tag' after the 'data', if this is a non-numeric,    */
    /* assume it is a file that can be relaxed.                      */
    strcpy(fname, master->trinfo_3d[t].relax_file);

    if (strlen(fname)) {
      rlx[master->nrlx] = t;
      master->nrlx++;

      /* If the data was read from a file, then check whether a      */
      /* relaxation is required.                                     */
      if (strlen(master->trinfo_3d[t].relax_dt)) {
        tm_scale_to_secs(master->trinfo_3d[t].relax_dt, &dt);
	if (strlen(master->trinfo_3d[t].r_rate)) {
	  double r0, r1, d0, d1;

	  /* Allocate memory for relaxation structure                */
	  relax = (tr_relax_data_t *)malloc(sizeof(tr_relax_data_t));
	  memset(relax, 0, sizeof(tr_relax_data_t));

	  /* Assign the time constant                                */
	  relax->tc = NULL;
	  if (sscanf(master->trinfo_3d[t].r_rate, "%lf %s", &tconst, buf) 
	      == 2) {
	    /* Relaxation constant is in a number with time units */
	    tm_scale_to_secs(master->trinfo_3d[t].r_rate, &tconst);
	    master->trinfo_3d[t].relax_rate = 1.0 / tconst;
	    hd_warn("Relaxation constant = %s for tracer %s\n",
		    master->trinfo_3d[t].r_rate, master->trinfo_3d[t].name);
	  } else if (sscanf(master->trinfo_3d[t].r_rate, "%s %lf %lf %s %lf %lf %s", 
			  buf, &d0, &r0, tu0, &d1, &r1, tu1) == 7) {
	    if(strcmp(buf, "linear") == 0) {
	      /* Adaptive relaxation */
	      hd_warn("Linear adaptive relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      relax->dt0 = d0;
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      sprintf(buf, "%f %s", r1, tu1);
	      tm_scale_to_secs(buf, &r1);
	      relax->rlx0 = r0;
	      relax->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	      relax->tctype = (RLX_ADPT|RLX_LINR);
	    } else if(strcmp(buf, "temporal") == 0) {
	      /* Temporally linear relaxation */
	      hd_warn("Temporally linear relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      relax->dt0 = d0;
	      relax->dt1 = d1;
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      sprintf(buf, "%f %s", r1, tu1);
	      tm_scale_to_secs(buf, &r1);
	      relax->rlx0 = r0;
	      relax->rlx1 = r1;
	      relax->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	      relax->tctype = (RLX_TIME);
	    } else if(strcmp(buf, "depth") == 0) {
	      /* Depth scaled linear relaxation */
	      hd_warn("Depth scaled linear relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      relax->dt0 = d0;
	      relax->dt1 = d1;
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      sprintf(buf, "%f %s", r1, tu1);
	      tm_scale_to_secs(buf, &r1);
	      relax->rlx0 = r0;
	      relax->rlx1 = r1;
	      relax->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	      relax->tctype = (RLX_ADPT|RLX_DEP);
	    } else if(strcmp(buf, "exp_depth") == 0) {
	      /* Depth scaled exponential relaxation                    */
	      /* rate = (r0 - r1)exp(-depth/d0) + (r1 - r0.exp(-d1/d0)) */
	      /* d0 = scaling factor (e.g. 500)                         */
	      /* r0 = rate at depth equal 0                             */
	      /* d1 = maximum depth                                     */
              /* r1 = rate at depth d1                                  */
	      /* e.g. exp_depth 500 100 days 4000 7 days                */
	      /* d0 ~ 40 for d1 = 200                                   */
	      /* d0 ~ 500 for d1 = 4000                                 */
	      hd_warn("Depth scaled exponential relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      relax->dt0 = fabs(d0);
	      relax->dt1 = fabs(d1);
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      sprintf(buf, "%f %s", r1, tu1);
	      tm_scale_to_secs(buf, &r1);
	      relax->rlx0 = r0;
	      relax->rlx1 = r1;
	      relax->slope = r0 * exp(-relax->dt1 / relax->dt0);
	      relax->tctype = (RLX_ADPT|RLX_EDEP);
	    } else if(strcmp(buf, "cos_depth") == 0) {
	      double d, r;
	      /* Depth scaled cosine relaxation                    */
	      /* rate = 0.5*((r0-r1)*cos(d*pi/(d1-d0)-d0*PI/(d1-d0))+(t0+t1)); */
	      /* This is designed for d0 < d1 */
	      hd_warn("Depth scaled cosine relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      d0 = fabs(d0);
	      d1 = fabs(d1);
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      sprintf(buf, "%f %s", r1, tu1);
	      tm_scale_to_secs(buf, &r1);
	      if (d0 > d1) {
		d = d0;
		d0 = d1;
		d1 = d;
		r = r0;
		r0 = r1;
		r1 = r;
	      }
	      relax->dt0 = fabs(d0);
	      relax->dt1 = fabs(d1);
	      relax->rlx0 = r0;
	      relax->rlx1 = r1;
	      relax->slope = d0 * PI / (d1 - d0);
	      relax->tctype = (RLX_ADPT|RLX_CDEP);
	    } else
	      hd_quit("Linear adaptive relaxation format, e.g; 'linear 10 1 day 0.1 10 days\n");
	  } else if (sscanf(master->trinfo_3d[t].r_rate, "%s %lf %lf %s", buf, &d0, &r0, tu0) == 4) {
	    if (strcmp(buf, "exponential") == 0) {
	      hd_warn("Exponential adaptive relaxation performed for tracer %s\n",
		      master->trinfo_3d[t].name);
	      tconst = 1.0;
	      sprintf(buf, "%f %s", r0, tu0);
	      tm_scale_to_secs(buf, &r0);
	      relax->dt0 = d0 * log(r0);
	      relax->rlx0 = r0;
	      relax->tctype = (RLX_ADPT|RLX_EXP);
	    } else
	      hd_quit("Exponential adaptive relaxation format, e.g; 'exponential 10 1 day'\n");
	  } else {
	    /* Relaxation constant is in a file */
	    relax->tc = hd_ts_read(master, master->trinfo_3d[t].r_rate, 0);
	    relax->tcid = ts_get_index(relax->tc, 
			     fv_get_varname(master->trinfo_3d[t].r_rate, 
					    "relaxation_time_constant", buf));
	    if (relax->tcid < 0)
	      hd_quit("%s does not contain 'relaxation_time_constant'.\n", 
		      master->trinfo_3d[t].r_rate);
	    tconst = ts_eval(relax->tc, relax->tcid, master->t);
	    sprintf(buf, "%f %s\n",tconst, relax->tc->varunit[relax->tcid]);
	    tm_scale_to_secs(buf, &tconst);
	    hd_warn("Relaxation constant read from file %s for tracer %s\n",
		    master->trinfo_3d[t].r_rate, master->trinfo_3d[t].name);
	  }
	  if (tconst > 0) {
            char schedName[MAXSTRLEN];
	    int tn, tnf = 1;

	    /* Populate and register the scheduler events.           */
            relax->master = master;
            relax->tid = t;
            relax->dt = dt;
	    /* Is the relaxation file is the name of a tracer in the list? */
	    for (tn = master->atr; tn < master->ntr; tn++) {
	      if (strcmp(fname, master->trinfo_3d[tn].name) == 0) {
		if (relax->tsnames == NULL)
		  relax->tsnames = (cstring *)malloc(sizeof(cstring));
		strcpy(relax->tsnames[0], fname);
		relax->varids[0] = tn;
		tnf = 0;
		break;
	      }
	    }
	    /* The relaxation file is indeed a file */
	    if (tnf) {
	      open_tr_relax_tsfiles(master, fname, relax);
	      hd_ts_multifile_check(relax->ntsfiles, relax->tsfiles,
				    relax->tsnames, master->trinfo_3d[t].name,
				    schedule->start_time,
				    schedule->stop_time);
	      hd_ts_multifile_get_index(relax->ntsfiles, relax->tsfiles,
					relax->tsnames, master->trinfo_3d[t].name,
					relax->varids);
	    }
            strcpy(schedName, "tracer_relax:");
            strcat(schedName, master->trinfo_3d[t].name);

            sched_register(schedule, schedName,
                           tr_relax_init, tr_relax_event, tr_relax_cleanup,
                           relax, NULL, NULL);
          } else
            hd_quit
              ("tr_relax_init: Invalid relaxation time constant specified for tracer '%s'.\n",
               master->trinfo_3d[t].name);

	  master->trinfo_3d[t].tctype = relax->tctype;
        } else
          hd_quit
            ("tr_relax_init: Relaxation time constant was not specified for tracer '%s'.\n",
             master->trinfo_3d[t].name);
      } else
          hd_quit
            ("tr_relax_init: Relaxation input time was not specified for tracer '%s'.\n",
             master->trinfo_3d[t].name);
    }
  }
  if (master->nrlx) {
    master->relax = i_alloc_1d(master->nrlx);
    for (t = 0; t < master->nrlx; t++)
      master->relax[t] = rlx[t];
  }
}


/* Free all remaining relaxation events */
void tracer_relax_end(master_t *master)
{
  int t;

  for (t = 0; t < master->ntr; ++t) {
    char schedName[MAXSTRLEN];
    strcpy(schedName, "tracer_relax:");
    strcat(schedName, master->trinfo_3d[t].name);
    sched_deregister(schedule, schedName);
  }
}


static int tr_relax_init(sched_event_t *event)
{
  return 1;
}

/* Read the relaxation data in. */
static double tr_relax_event(sched_event_t *event, double t)
{
  tr_relax_data_t *relax = (tr_relax_data_t *)schedGetPublicData(event);
  master_t *master = relax->master;
  geometry_t *geom = master->geom;
  double tsout = schedule->stop_time;
  int cc, c, cs;
  int tid = relax->tid;
  double ntval, dtr;
  double sf = 1.0 / 86400.0;

  /*
   * Only relax on reanalysis
   * Note1: We actually end up relaxing on the very first timestep
   *        regardless of DA type because master->da = 0. This may be bad for DA
   *        So, we have 2 choices, either ignore first OR
   *        Rejig the input file so that analysis field is the same as 
   *        the state variables.
   * Note2: The second time around tsout will cause it not to called
   *        until da rejigs all the stop times
   *
   */
  if (master->da && !(master->da & (NONE|DO_DA))) {
    relax->da_started = 0;
    return(tsout);
  }
  
  /* Check whether the current time exceeds the next update time */
  if ((t + relax->dt / 10.0) >= (event->next_event - SEPS)) {
    double *trv = master->tr_wc[tid];
    double relax_rate;
    relax->relax_rate = 0.0;
    if (relax->tctype & RLX_TIME) master->trinfo_3d[tid].relax_rate = 0.0;

    if (master->rtemp && tid == master->tno) {
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	if (relax->ntsfiles == 0 && relax->varids[0] > -1)
	  master->rtemp[c] = master->tr_wc[relax->varids[0]][c];
	else 
	  master->rtemp[c] = hd_ts_multifile_eval_xyz(relax->ntsfiles,
			   relax->tsfiles, relax->varids, t,
			   geom->cellx[cs], geom->celly[cs],
			   geom->cellz[c] * master->Ds[cs]);
	if (relax->tctype & (RLX_DEP|RLX_EDEP|RLX_CDEP))
	  dtr = geom->botz[geom->m2d[c]];
	else
	  dtr = trv[c] - master->rtemp[c];
	relax_rate = get_relax_rate(master, relax, dtr);
	/* Set for do_ts_relax */
	if (relax->tctype & (RLX_LINR|RLX_EXP|RLX_DEP|RLX_EDEP|RLX_CDEP)) {
	  master->trinfo_3d[tid].relax_rate = 0.0;
	  master->temp_tc[c] = sf / relax_rate;
	} else
	  master->trinfo_3d[tid].relax_rate = relax_rate;
      }
    } else if (master->rsalt && tid == master->sno) {
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	if (relax->ntsfiles == 0 && relax->varids[0] > -1)
	  master->rsalt[c] = master->tr_wc[relax->varids[0]][c];
	else
	  master->rsalt[c] = hd_ts_multifile_eval_xyz(relax->ntsfiles,
			   relax->tsfiles, relax->varids, t,
			   geom->cellx[cs], geom->celly[cs],
			   geom->cellz[c] * master->Ds[cs]);
	if (relax->tctype & (RLX_DEP|RLX_EDEP|RLX_CDEP))
	  dtr = geom->botz[geom->m2d[c]];
	else
	  dtr = trv[c] - master->rsalt[c];
	relax_rate = get_relax_rate(master, relax, dtr);
	/* Set for do_ts_relax */
	if (relax->tctype & (RLX_LINR|RLX_EXP|RLX_DEP|RLX_EDEP|RLX_CDEP)) {
	  master->trinfo_3d[tid].relax_rate = 0.0;
	  master->salt_tc[c] = sf / relax_rate;
	} else
	  master->trinfo_3d[tid].relax_rate = relax_rate;
      }
    } else {
      int trn = relax->varids[0];
      if (trn < master->ntr && 
	  strcmp(relax->tsnames[0], master->trinfo_3d[trn].name) == 0) {
	for (cc = 1; cc <= geom->b3_t; cc++) {
	  c = geom->w3_t[cc];
	  cs = geom->m2d[c];
	  ntval = master->tr_wc[trn][c];
	  if (relax->tctype & (RLX_DEP|RLX_EDEP|RLX_CDEP))
	    dtr = geom->botz[geom->m2d[c]];
	  else
	    dtr = trv[c] - ntval;
	  relax_rate = get_relax_rate(master, relax, dtr);
	  trv[c] -= relax->dt * relax_rate * (trv[c] - ntval);
	}
      } else {
	for (cc = 1; cc <= geom->b3_t; cc++) {
	  c = geom->w3_t[cc];
	  cs = geom->m2d[c];
	  ntval = hd_ts_multifile_eval_xyz(relax->ntsfiles,
					   relax->tsfiles,
					   relax->varids, t,
					   geom->cellx[cs], geom->celly[cs],
					   geom->cellz[c] * master->Ds[cs]);
	  if (relax->tctype & (RLX_DEP|RLX_EDEP|RLX_CDEP))
	    dtr = geom->botz[geom->m2d[c]];
	  else
	    dtr = trv[c] - ntval;
	  relax_rate = get_relax_rate(master, relax, dtr);
	  trv[c] -= relax->dt * relax_rate * (trv[c] - ntval);
	}
      }
    }

    /* Set the next update time */
    event->next_event += relax->dt;
    tsout = event->next_event;
  }
  return tsout;
}


/* Cleanup routines */
static void tr_relax_cleanup(sched_event_t *event, double t)
{
  tr_relax_data_t *relax = (tr_relax_data_t *)schedGetPublicData(event);
  master_t *master = relax->master;
  int i;

  if (relax == NULL)
    return;

  for (i = 0; i < relax->ntsfiles; ++i)
    hd_ts_free(master, relax->tsfiles[i]);
  free(relax->tsfiles);
  free(relax->tsnames);
  free(relax);
}

static double get_relax_rate(master_t *master, tr_relax_data_t *relax, double dtr)
{
  int tid = relax->tid;
  double rr;

  if (relax->relax_rate) {
    /* Relaxation rate has been set for this timestep */
    return(relax->relax_rate);
  } else if (relax->tc != NULL) {
    /* Relaxation rate read from file. Temporally variable only. */
    char buf[MAXSTRLEN];
    double tt = master->t;
    if (master->da && !relax->da_started && ts_is_modulo(relax->tc)) {
      relax->da_started = 1;
      tt += 10.; // forces wrap around of modulo file
    }
    rr = ts_eval(relax->tc, relax->tcid, tt);
    sprintf(buf, "%f %s\n",rr, relax->tc->varunit[relax->tcid]);
    tm_scale_to_secs(buf, &rr);
    relax->relax_rate = 1.0 / rr;
    //printf("(%d) Relax rate at %f(%f) is %f\n", (int)master->da, master->t,tt, rr);
    //fflush(NULL);
    return(relax->relax_rate);
  } else if (master->trinfo_3d[tid].relax_rate) {
    /* Constant relaxation rate */
    relax->relax_rate = master->trinfo_3d[tid].relax_rate;
    return(relax->relax_rate);
  } else {
    /* Adaptive relaxation */
    if (relax->tctype & RLX_LINR) {
      double dmin = min(relax->dt0, relax->dt1);
      double dmax = max(relax->dt0, relax->dt1);
      dtr = min(max(fabs(dtr), dmin), dmax);
      rr = (fabs(dtr) - relax->dt0) * relax->slope + relax->rlx0;
      rr = 1.0 / rr;
    } else if (relax->tctype & RLX_TIME) {
      rr = min((master->days - relax->dt0) * relax->slope + relax->rlx0, relax->rlx1);
      rr = relax->relax_rate = 1.0 / rr;
    } else if (relax->tctype & RLX_DEP) {
      double dmin = min(relax->dt0, relax->dt1);
      double dmax = max(relax->dt0, relax->dt1);
      dtr = min(max(dtr, dmin), dmax);
      rr = (dtr - relax->dt0) * relax->slope + relax->rlx0;
      rr = 1.0 / rr;
    } else if (relax->tctype & RLX_EDEP) {
      double rmin = min(relax->rlx0, relax->rlx1);
      double rmax = max(relax->rlx0, relax->rlx1);
      rr = max(rmin,
	       (rmax - rmin) * exp(-fabs(dtr) / relax->dt0) + (rmin - relax->slope));
      rr = 1.0 / rr;
    } else if (relax->tctype & RLX_CDEP) {
      double r0 = relax->rlx0;
      double r1 = relax->rlx1;
      double d = fabs(dtr);
      double dd = relax->dt1 - relax->dt0;
      if (d < relax->dt0)
	rr = r0;
      else if (d > relax->dt1)
	rr = r1;
      else
	rr = 0.5 * ((r0 - r1) * cos(d * PI / dd - relax->slope) + (r0 + r1));
      rr = 1.0 / rr;
    } else if (relax->tctype & RLX_EXP) {
      rr = exp(relax->dt0 / fabs(dtr));
      rr = min(max(rr, relax->rlx0 / 100.0), relax->rlx0 * 100.0);
      rr = 1.0 / rr;
    }
    return(rr);
  }
}


/*------------------------------------------------------------------*/
/* Routine to relax tracers at a given cell.                        */
/* Usage:                                                           */
/* tr_relax_event_c(                                                */
/* sched_get_even_by_name(schedule, "tracer_relax:trname"),         */
/* windat->t, c);                                                   */
/* trname = temp, salt etc.                                         */
/*------------------------------------------------------------------*/
void tr_relax_event_c(sched_event_t *event, double relax_rate,
		      double dt, double t, int c)
{
  tr_relax_data_t *relax = (tr_relax_data_t *)schedGetPublicData(event);
  master_t *master = relax->master;
  geometry_t *geom = master->geom;
  int cs = geom->m2d[c];
  int tid = relax->tid;
  double *trv = master->tr_wc[tid];
  double ntval;

  ntval = hd_ts_multifile_eval_xyz(relax->ntsfiles,
				   relax->tsfiles,
				   relax->varids, t,
				   geom->cellx[cs], geom->celly[cs],
				   geom->cellz[c] * master->Ds[cs]);
  trv[c] -= dt * (trv[c] - ntval) / relax_rate;
}

/* END tr_relax_event_c()                                           */
/*------------------------------------------------------------------*/
