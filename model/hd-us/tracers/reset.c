/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/tracers/reset.c
 *  
 *  Description:
 *  reset.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: reset.c 7247 2022-10-25 00:28:46Z her127 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"


/* Local functions */
static int tr_reset_init(sched_event_t *event);
double tr_reset_event(sched_event_t *event, double t);
double tr_reset2d_event(sched_event_t *event, double t);
static void tr_reset_cleanup(sched_event_t *event, double t);
static double trans_reset_event(sched_event_t *event, double t);

int get_datum(int ntsfiles, timeseries_t **tsfiles, cstring * names,
	      char *varname,  char *attname);

typedef struct {
  master_t *master;             /* Grid associated with */
  int tid;                      /* Tracer id */
  int ntsfiles;                 /* Number of time-series files */
  timeseries_t **tsfiles;       /* Array of time-series files */
  cstring *tsnames;             /* Array of time-series files */
  double dt;                    /* Relaxation time step */
  double tnext;                 /* Next time step */
  int in_type;                  /* Type of reset file to read */
  double fact;                  /* Scaling factor */
  char interp_type[MAXSTRLEN];  /* Interpolation type (optional) */
  int start_index;              /* Start index for UGRID */
  poly_t *pl;                   /* Inclusion /exclusion polygon */
} tr_reset_data_t;


/*-------------------------------------------------------------------*/
/* Opens files for tracer resets                                     */
/*-------------------------------------------------------------------*/
static void open_tr_reset_tsfiles(master_t *master, char *fnames,
                                  tr_reset_data_t *reset)
{
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int i;
  int nf = parseline(fnames, fields, MAXNUMARGS);

  reset->ntsfiles = nf;

  if (nf > 0) {
    reset->tsfiles = (timeseries_t **)malloc(sizeof(timeseries_t *) * nf);
    reset->tsnames = (cstring *) malloc(sizeof(cstring) * nf);
    memset(reset->tsfiles, 0, sizeof(timeseries_t *) * nf);
    memset(reset->tsnames, 0, sizeof(cstring) * nf);
    if (is_set_variable(fnames)) {
      reset->in_type |= SP_SET;
      for (i = 0; i < nf; ++i)
	strcpy(reset->tsnames[i], fields[i]);
    } else {
      for (i = 0; i < nf; ++i) {
	strcpy(reset->tsnames[i], fields[i]);
	reset->tsfiles[i] = hd_ts_read_us(master, reset->tsnames[i], 0, reset->interp_type);
      }
    }
  } else
    hd_quit("Must specify at least one reset time-series file.\n");
}

/* END open_tr_reset_tsfiles()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Check whether the tr_reset_data_t are to be included, if so, then */
/* register them with the scheduler.                                 */
/*-------------------------------------------------------------------*/
void tracer_reset_init(master_t *master)
{
  char fname[MAXSTRLEN];
  char buf[MAXSTRLEN];
  double dt;
  int t;
  tr_reset_data_t *reset = NULL;
  int rst[master->ntr];

  master->nres = 0;
  prm_set_errfn(hd_silent_warn);

  for (t = 0; t < master->ntr; t++) {

    /* Read the 'tag' after the 'data', if this is a non-numeric,    */
    /* assume it is a file that can be reset.                        */
    strcpy(fname, master->trinfo_3d[t].reset_file);
    if (strlen(fname)) {
      rst[master->nres] = t;
      master->nres++;

      /* If the data was read from a file, then check whether a      */
      /* resetation is required.                                     */
      if (strlen(master->trinfo_3d[t].reset_dt)) {
	char schedName[MAXSTRLEN];
        tm_scale_to_secs(master->trinfo_3d[t].reset_dt, &dt);

	if (master->save_force &&
	    (strcmp(master->trinfo_3d[t].name, "otemp") == 0 ||
	     strcmp(master->trinfo_3d[t].name, "osalt") == 0 ||
	     strcmp(master->trinfo_3d[t].name, "ovelu") == 0 ||
	     strcmp(master->trinfo_3d[t].name, "ovelv") == 0))
	  dt =86400.0;

	/* Allocate memory for resetation structure, populate and    */
	/* register the scheduler events.                            */
	reset = (tr_reset_data_t *)malloc(sizeof(tr_reset_data_t));
	memset(reset, 0, sizeof(tr_reset_data_t));

	reset->master = master;
	reset->tid = t;
	reset->dt = dt;

	/* Interpolation rule */
	strcpy(reset->interp_type, master->trinfo_3d[t].reset_interp);

	/* Set the initial time one timestep previous. See           */
	/* tr_reset_event() for the rationale for this.              */
	reset->tnext = schedule->start_time - master->grid_dt;
	open_tr_reset_tsfiles(master, fname, reset);
	hd_ts_multifile_check(reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, master->trinfo_3d[t].name,
			      schedule->start_time,
			      schedule->stop_time);
	strcpy(schedName, "tracer_reset:");
	strcat(schedName, master->trinfo_3d[t].name);

	/* Check if the input file is sparse format                  */
	reset->in_type = XYZ_TINT;
	if(check_sparse_dumpfile(geom, reset->ntsfiles, 
				 reset->tsfiles, reset->tsnames) == 0)
	  if (master->tmode & SP_STRUCT) 
	    reset->in_type = SP_TINT;

	/* Check if there is a polygon zone                          */
	reset->pl = NULL;
	if (master->trinfo_3d[t].flag & (P_IN|P_EX)) {
	  FILE *fp;
	  reset->pl = poly_create();
	  if ((fp = fopen(master->trinfo_3d[t].tag, "r")) != NULL) {
	    poly_read(reset->pl, fp);
	    reset->in_type |= master->trinfo_3d[t].flag;
	    fclose(fp);
	  }
	}

	sched_register(schedule, schedName,
		       tr_reset_init, tr_reset_event, tr_reset_cleanup,
		       reset, NULL, NULL);

	if (strcmp(master->trinfo_3d[t].tag, "SWR_INVERSE") == 0)
	  master->trinfo_3d[t].flag |= SWR_INVERSE;
      } else
          hd_quit
            ("tr_reset_init: Resetation input time was not specified for tracer '%s'.\n",
             master->trinfo_3d[t].name);
    }
  }

  if (master->nres) {
    master->reset = i_alloc_1d(master->nres);
    for (t = 0; t < master->nres; t++)
      master->reset[t] = rst[t];
  }
}

/* END tracer_reset_init()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Free all remaining reset events                                   */
/*-------------------------------------------------------------------*/
void tracer_reset_end(master_t *master)
{
  int t;

  for (t = 0; t < master->ntr; ++t) {
    char schedName[MAXSTRLEN];
    strcpy(schedName, "tracer_reset:");
    strcat(schedName, master->trinfo_3d[t].name);
    sched_deregister(schedule, schedName);
  }
}

/* END tracer_reset_end()                                            */
/*-------------------------------------------------------------------*/


static int tr_reset_init(sched_event_t *event)
{
  return 1;
  /*  return tr_reset_event(event, schedule->start_time);*/
}


/*-------------------------------------------------------------------*/
/* Read the reset data in                                            */
/*-------------------------------------------------------------------*/
double tr_reset_event(sched_event_t *event, double t)
{
  tr_reset_data_t *reset = (tr_reset_data_t *)schedGetPublicData(event);
  master_t *master = reset->master;
  geometry_t *geom = master->geom;
  double tsout = schedule->stop_time;
  double tsin;
  int tid = reset->tid;

  /* Check whether the current time exceeds the next update time     */
  if ((t + reset->dt / 10.0) >= (reset->tnext - master->dt)) {
    double *trv = master->tr_wc[tid];

    /* Note : there is a sequence where reset updates, tracer        */
    /* transport, tracer statistics and output dumps are performed.  */
    /* If a reset is used in conjunction with tracerstats to get a   */
    /* difference between a reset and defined tracer at a time t,    */
    /* then the order of the sequence must be :                      */
    /* (a) Reset the tracer at time t                                */
    /* (b) Update the defined tracer via transport to time t         */
    /* (c) Perform the difference in tracerstats at time t           */
    /* (d) Output the dumpfile at time t                             */
    /* The only way to syncronize this sequence using the scheduler  */
    /* is to read the reset tracer in for time t at the simulation   */
    /* time t-dt. Therefore, in tracer_reset_init() the tnext time   */
    /* is initialised to (start_time - grid_dt), and when the reset  */
    /* event is triggered at this time, use (tnext + grid_dt) as the */
    /* input time for the eval() function.                           */
    tsin = reset->tnext + master->grid_dt;
    hd_ts_multifile_evalp(master, reset->ntsfiles, reset->tsfiles,
			  reset->tsnames, master->trinfo_3d[tid].name,
			  trv, tsin, geom->wsa, geom->a3_t, 
			  reset->in_type, reset->pl);

    /* Set the next update time                                      */
    reset->tnext += reset->dt;
    tsout = reset->tnext;

    /* Scale if required                                             */
    scale_tracer(master->trinfo_3d, geom, master->tr_wc, tid);

    /* Set the flag for increment updates                            */
    if (master->trinfo_3d[tid].increment) master->trinc[tid] = 1;

    if (master->trinfo_3d[tid].flag & SWR_INVERSE) {
      int n = tracer_find_index("temp", master->ntr, master->trinfo_3d);
      master->trinfo_3d[n].flag = DO_SWR_INVERSE;
    }

    if (strcmp(master->trinfo_3d[tid].long_name, "u1") == 0)
      memcpy(master->u1, trv, geom->sgsiz * sizeof(double));
    if (strcmp(master->trinfo_3d[tid].long_name, "u2") == 0)
      memcpy(master->u2, trv, geom->sgsiz * sizeof(double));
  }

  return tsout;
}

/* END tr_reset_event()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Cleanup routines                                                  */
/*-------------------------------------------------------------------*/
static void tr_reset_cleanup(sched_event_t *event, double t)
{
  tr_reset_data_t *reset = (tr_reset_data_t *)schedGetPublicData(event);
  master_t *master = reset->master;
  int i;

  if (reset == NULL)
    return;

  if(reset->pl) poly_destroy(reset->pl); 
  for (i = 0; i < reset->ntsfiles; ++i)
    hd_ts_free(master, reset->tsfiles[i]);
  free(reset->tsfiles);
  free(reset->tsnames);
  free(reset);
}

/* END tr_reset_cleanup()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Check whether the tr_reset_data_t are to be included, if so, then */
/* register them with the scheduler.                                 */
/*-------------------------------------------------------------------*/
void tracer_reset2d_init(master_t *master)
{
  char fname[MAXSTRLEN];
  double dt;
  int t;
  tr_reset_data_t *reset = NULL;
  int rst[master->ntrS];

  master->nres2d = 0;
  prm_set_errfn(hd_silent_warn);
  for (t = 0; t < master->ntrS; t++) {

    /* Read the 'tag' after the 'data', if this is a non-numeric,    */
    /* assume it is a file that can be reset.                        */
    strcpy(fname, master->trinfo_2d[t].reset_file);
    if (strlen(fname)) {
      rst[master->nres2d] = t;
      master->nres2d++;

      /* If the data was read from a file, then check whether a      */
      /* resetation is required.                                     */
      if (strlen(master->trinfo_2d[t].reset_dt)) {
	char schedName[MAXSTRLEN];
        tm_scale_to_secs(master->trinfo_2d[t].reset_dt, &dt);

	/* Allocate memory for resetation structure, populate and    */
	/* register the scheduler events.                            */
	reset = (tr_reset_data_t *)malloc(sizeof(tr_reset_data_t));
	memset(reset, 0, sizeof(tr_reset_data_t));

	reset->master = master;
	reset->tid = t;
	reset->dt = dt;
	reset->in_type = (XYZ_TINT|VEL2D);

	/* Interpolation rule */
	strcpy(reset->interp_type, master->trinfo_2d[t].reset_interp);

	/* Set the initial time one timestep previous. See           */
	/* tr_reset_event() for the rationale for this.              */
	reset->tnext = schedule->start_time - master->grid_dt;
	open_tr_reset_tsfiles(master, fname, reset);
	if (reset->in_type & SP_SET) {
	  if (!strlen(reset->interp_type))
	    hd_quit("tracer_reset2d(): Requite 'reset_interp' for %s\n", fname);
	} else
	  hd_ts_multifile_check(reset->ntsfiles, reset->tsfiles,
				reset->tsnames, master->trinfo_2d[t].name,
				schedule->start_time,
				schedule->stop_time);

	strcpy(schedName, "tracer_reset2d:");
	strcat(schedName, master->trinfo_2d[t].name);

	if (strcmp(master->trinfo_2d[t].name, "ghrsst") == 0)
	  reset->in_type |= GHRSST;

	if(!(reset->in_type & SP_SET) && check_sparse_dumpfile(geom, reset->ntsfiles, 
				 reset->tsfiles, reset->tsnames) == 0)
	  reset->in_type = SP_TINT;

	/* Check if there is a polygon zone                          */
	reset->pl = NULL;
	if (master->trinfo_2d[t].flag & (P_IN|P_EX)) {
	  FILE *fp;
	  reset->pl = poly_create();
	  if ((fp = fopen(master->trinfo_2d[t].tag, "r")) != NULL) {
	    poly_read(reset->pl, fp);
	    reset->in_type |= master->trinfo_2d[t].flag;
	    fclose(fp);
	  }
	}

        sched_register(schedule, schedName,
                       tr_reset_init, tr_reset2d_event, tr_reset_cleanup,
                       reset, NULL, NULL);
      } else
          hd_quit
            ("tr_reset2d_init: Resetation input time was not specified for tracer '%s'.\n",
             master->trinfo_2d[t].name);
    }
  }
  if (master->nres2d) {
    master->reset2d = i_alloc_1d(master->nres2d);
    for (t = 0; t < master->nres2d; t++)
      master->reset2d[t] = rst[t];
  }
}

/* END tracer_reset2d_init()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Free all remaining reset events                                   */
/*-------------------------------------------------------------------*/
void tracer_reset2d_end(master_t *master)
{
  int t;

  for (t = 0; t < master->ntr; ++t) {
    char schedName[MAXSTRLEN];
    strcpy(schedName, "tracer_reset2d:");
    strcat(schedName, master->trinfo_2d[t].name);
    sched_deregister(schedule, schedName);
  }
}

/* END tracer_reset2d_end()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read the reset data in                                            */
/*-------------------------------------------------------------------*/
double tr_reset2d_event(sched_event_t *event, double t)
{
  tr_reset_data_t *reset = (tr_reset_data_t *)schedGetPublicData(event);
  master_t *master = reset->master;
  geometry_t *geom = master->geom;
  double tsout = schedule->stop_time;
  double tsin;
  int tid = reset->tid;

  /* Check whether the current time exceeds the next update time     */
  if ((t + reset->dt / 10.0) >= reset->tnext) {
    double *trv = master->tr_wcS[tid];

    tsin = reset->tnext + master->grid_dt;

    /* Switch on interp type, if specified */
    if (strlen(reset->interp_type)) {
      char buf[MAXLINELEN];

      /*-----------------------------------------------------------------*/
      /* Check for regionalized formats                                  */
      if (reset->in_type & SP_SET) {
	sprintf(buf, "[var=%s] [i_rule=%s] %s", master->trinfo_2d[tid].name, reset->interp_type, 
		master->trinfo_2d[tid].reset_file);
	set_variable(master, buf, trv, &tsin);
      } else {
	/* Use the grid interp library */
	/* NOTE: Only reads the first file!! */
	/*
	char *vname = fv_get_varname(reset->tsnames[0], master->trinfo_2d[tid].name, buf);
	hd_ts_grid_interp(master, reset->tsfiles[0], vname,
			  trv, tsin, geom->wsa, geom->a2_t, reset->interp_type);
	*/
	hd_ts_multifile_evalp(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, master->trinfo_2d[tid].name,
			      trv, tsin, geom->wsa, geom->a2_t, 
			      reset->in_type, reset->pl);

      }
			
    } else {
      /* Standard emslib interpolation */
      hd_ts_multifile_evalp(master, reset->ntsfiles, reset->tsfiles,
			    reset->tsnames, master->trinfo_2d[tid].name,
			    trv, tsin, geom->wsa, geom->a2_t, 
			    reset->in_type, reset->pl);
    }
    /* Scale if required                                             */
    scale_tracer(master->trinfo_2d, geom, master->tr_wcS, tid);

    /* If we read state variable quantities into the model (e.g. in  */
    /* transport mode = NONE, then copy them into relevant state     */
    /* variable arrays.                                              */
    /* Note: for 2D (u,v) variables copy to the surface layer of the */
    /* 3D state variables (u,v); for SWAN coupling.                  */
    if (strcmp(master->trinfo_2d[tid].long_name, "eta") == 0)
      memcpy(master->eta, trv, geom->szcS * sizeof(double));
    if (strcmp(master->trinfo_2d[tid].long_name, "u") == 0)
      memcpy(master->u, trv, geom->szcS * sizeof(double));
    if (strcmp(master->trinfo_2d[tid].long_name, "v") == 0)
      memcpy(master->v, trv, geom->szcS * sizeof(double));

    /* Set the next update time                                      */
    reset->tnext += reset->dt;
    tsout = reset->tnext;

    /* Set the flag for increment updates                            */
    if (master->trinfo_2d[tid].increment) master->trincS[tid] = 1;

    /*
    if (strcmp(master->trinfo_2d[tid].name, "eta") == 0) 
      memcpy(master->eta, trv, geom->szcS * sizeof(double));
    if (strcmp(master->trinfo_2d[tid].name, "uav") == 0) 
      memcpy(master->uav, trv, geom->szcS * sizeof(double));
    if (strcmp(master->trinfo_2d[tid].name, "vav") == 0) 
      memcpy(master->vav, trv, geom->szcS * sizeof(double));
    */
  }

  return tsout;
}

/* END tr_reset2d_event()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Check whether the transport variables are required to be reset    */
/* and register them with the scheduler.                             */
/*-------------------------------------------------------------------*/
void trans_reset_init(parameters_t *params, master_t *master)
{
  char fname[MAXSTRLEN];
  int varids[MAXNUMTSFILES];
  double dt;
  int tn;
  tr_reset_data_t *reset_trans = NULL;

  prm_set_errfn(hd_silent_warn);
  strcpy(fname, params->trans_data);
  dt = master->grid_dt;
  if (master->tmode & (NONE|SP_DUMP)) return;

  /* Allocate memory for resetation structure, populate and    */
  /* register the scheduler events.                            */
  reset_trans = (tr_reset_data_t *)malloc(sizeof(tr_reset_data_t));
  memset(reset_trans, 0, sizeof(tr_reset_data_t));
  reset_trans->master = master;
  reset_trans->dt = dt;
  reset_trans->tnext = schedule->start_time;
  if (master->tmode & SP_SIMPLE) strcpy(reset_trans->interp_type, "nearest_eps");
  if (master->tmode & SP_SIMPLEU) {
    strcpy(reset_trans->interp_type, "nearest");
    master->tmode |= SP_SIMPLE;
  }
  open_tr_reset_tsfiles(master, fname, reset_trans);

  if (master->tmode & SP_ORIGIN) {
    if (hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			      reset_trans->tsnames, "origin",
			      schedule->start_time, schedule->stop_time)) {
      /* Read the datum for the streamline origin from the attributes */
      master->togn = get_datum(reset_trans->ntsfiles, reset_trans->tsfiles,
			       reset_trans->tsnames, "origin", "datum");
      if (master->togn == -1)
	hd_warn("trans_reset_init: can't find 'datum' attribute in variable 'origin'.\n");
      
    } else
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "p",
			  schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "q",
			  schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "r",
			  schedule->start_time, schedule->stop_time);
    /*if (master->fillf & MONOTONIC && !(master->fillf & SET_BDRY)) {*/
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u1",
			  schedule->start_time, schedule->stop_time);
    /*
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u2",
			  schedule->start_time, schedule->stop_time);
    */
  } else if (master->tmode & (GLOBAL|SP_SIMPLE)) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u",
			  schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "v",
			  schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "eta",
			  schedule->start_time, schedule->stop_time);
    if (master->tmode & SP_SIMPLE)
      hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			    reset_trans->tsnames, "w",
			    schedule->start_time, schedule->stop_time);
  } else {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u1",
			  schedule->start_time, schedule->stop_time);
    /*
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u2",
			  schedule->start_time, schedule->stop_time);
    */
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "w",
			  schedule->start_time, schedule->stop_time);
  }
  if (master->tmode & SP_FFSL) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u1vmean",
			  schedule->start_time, schedule->stop_time);
    /*
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "u2vmean",
			  schedule->start_time, schedule->stop_time);
    */
  }
  if (!(master->tmode & (GLOBAL|SP_SIMPLE))) {
    if (!master->do_closure)
      hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			    reset_trans->tsnames, "Kz",
			    schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "eta",
			  schedule->start_time, schedule->stop_time);
  }
  for (tn = 0; tn < master->ntrvars; tn++) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, master->trvars[tn],
			  schedule->start_time, schedule->stop_time);    
  }
  for (tn = 0; tn < master->ntrvarsS; tn++) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, master->trvarsS[tn],
			  schedule->start_time, schedule->stop_time);    
  }
  if (!master->trinfo_3d[master->tno].advect && 
      !master->trinfo_3d[master->tno].diffuse)
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "temp",
			  schedule->start_time, schedule->stop_time);
  if (!master->trinfo_3d[master->sno].advect && 
      !master->trinfo_3d[master->sno].diffuse)
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "salt",
			  schedule->start_time, schedule->stop_time);

  if (master->light && hd_ts_multifile_get_index(reset_trans->ntsfiles, reset_trans->tsfiles,
						 reset_trans->tsnames, "swr", varids) != 0) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "swr",
			  schedule->start_time, schedule->stop_time);
    master->tmode |= DO_SWR;
  }

  if (hd_ts_multifile_get_index(reset_trans->ntsfiles, reset_trans->tsfiles,
				reset_trans->tsnames, "dz", varids) != 0 &&
      hd_ts_multifile_get_index(reset_trans->ntsfiles, reset_trans->tsfiles,
				reset_trans->tsnames, "dzu1", varids) != 0) {
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "dz",
			  schedule->start_time, schedule->stop_time);
    hd_ts_multifile_check(reset_trans->ntsfiles, reset_trans->tsfiles,
			  reset_trans->tsnames, "dzu1",
			  schedule->start_time, schedule->stop_time);
    master->tmode |= DO_DZ;
  }

  if (master->tmode & SP_UGRID) {
    reset_trans->start_index = 1;
    if (strcmp(hd_ts_multifile_get_text_att(reset_trans->ntsfiles, reset_trans->tsfiles,
					    reset_trans->tsnames, "Mesh2_face_nodes",
					    "start_index"), "0") == 0)
      reset_trans->start_index = 0;
    if (reset_trans->start_index == 1) master->tmode |= SP_START1;
  }
  sched_register(schedule, "transport_reset",
		 tr_reset_init, trans_reset_event, tr_reset_cleanup,
		 reset_trans, NULL, NULL);

  if (master->tmode & (SP_EXACT|SP_TINT)) {
    if(check_sparse_dumpfile(master->geom->trans->tpg, 
			     reset_trans->ntsfiles, 
			     reset_trans->tsfiles,
			     reset_trans->tsnames))
      hd_quit("trans_reset_init: Input TRANS_DATA file incompatible with grid.\n");
  }

}

/* END trans_reset_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Free all remaining reset events                                   */
/*-------------------------------------------------------------------*/
void trans_reset_end(master_t *master)
{

  sched_deregister(schedule, "transport_reset");
}

/* END trans_reset_end()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Read the reset data in                                            */
/*-------------------------------------------------------------------*/
static double trans_reset_event(sched_event_t *event, double t)
{
  tr_reset_data_t *reset = (tr_reset_data_t *)schedGetPublicData(event);
  master_t *master = reset->master;
  geometry_t *geom = master->geom;
  double tsout = schedule->stop_time;
  int mode = master->tmode, n, tn;
  int *vc3 = geom->w3_t;
  int *vc2 = geom->w2_t;
  int nc2 = geom->b2_t;
  int nc3 = geom->b3_t;
  int *ve3 = geom->w3_e1;
  int *ve2 = geom->w2_e1;
  int ne2 = geom->n2_e1;
  int ne3 = geom->n3_e1;
  int c,cs,cc;
  double tsin = t;
  double tsync = 0.0;

  /* Check whether the current time exceeds the next update time     */

  /* The dump, hd_step and reset scheduling are synced in that order */
  /* (see comments in tr_reset_event). This means that dumps or      */
  /* resets are one time-step out of sync and the (p,q,r) read at    */
  /* next tiume-step should be used for the current hd_step().       */
  if (master->tmode & SP_ORIGIN) {
    tsync = master->dt;
    tsin = reset->tnext + master->grid_dt;
  }
  if (master->tmode & SP_STRUCT) {
    ve3 = geom->s2c;
    ve2 = geom->s2c;
    ne2 = geom->ns2;
    ne3 = geom->ns3;
    vc3 = geom->s2c;
    vc2 = geom->s2c;
    nc2 = geom->ns2;
    nc3 = geom->ns3;
  }
  if (master->tmode & SP_UGRID) {
    vc3 = geom->w2_t;
    nc3 = geom->b2_t;
    ve3 = geom->w2_e1;
    ne3 = geom->n2_e1;
  }

  if ((t + reset->dt / 10.0) >= (reset->tnext - tsync)) {
    if (master->tmode & SP_ORIGIN) {
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "origin", master->origin, tsin,
			      vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "p", master->pc, tsin,
			      vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "q", master->qc, tsin,
			      vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "r", master->rc, tsin,
			      vc3, nc3, mode);
      /*if (master->fillf & MONOTONIC && !(master->fillf & SET_BDRY)) {*/
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "u1", master->u1, tsin,
			      ve3, ne3, (mode|U1GEN));
    }
    if (master->tmode & (GLOBAL|SP_SIMPLE)) {
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "u", master->u, tsin, vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "v", master->v, tsin, vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "eta", master->eta, tsin,
			      vc2, nc2, (mode|VEL2D));
      if (master->tmode & SP_SIMPLE)
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "w", master->w, tsin, vc3, nc3, 
				mode);
    }
    if (master->tmode & SP_FFSL) {
      int nmode = (mode|U1GEN);
      if (master->tmode & SP_STRUCT) {
	nmode |= U2STRUCT;
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "u2vmean", master->u1vm, tsin,
				ve3, ne3, nmode);
	nmode &= ~U2STRUCT;
	nmode |= U1STRUCT;
      }
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "u1vmean", master->u1vm, tsin,
			      ve3, ne3, nmode);
      if (master->tmode & SP_STRUCT) nmode &= ~U1STRUCT;
    }
    /* Read in mandatory variables                                   */
    if (!(master->tmode & (SP_ORIGIN|GLOBAL|SP_SIMPLE))) {
      int nmode = (mode|U1GEN);
      if (master->tmode & SP_STRUCT) {
	nmode |= U2STRUCT;
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "u2", master->u1, tsin,
				ve3, ne3, nmode);
	nmode &= ~U2STRUCT;
	nmode |= U1STRUCT;
      }
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "u1", master->u1, tsin,
			      ve3, ne3, nmode);
      if (master->tmode & SP_STRUCT) nmode &= ~U1STRUCT;
      if (!(master->conserve & CONS_W) || master->conserve & (CONS_PSS|CONS_WS)) {
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "w", master->w, tsin,
				vc3, nc3, mode);
	master->w[0] = 1.0;
      }
    }
    if (!(master->tmode & (GLOBAL|SP_SIMPLE))) {
      if (!master->do_closure)
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "Kz", master->Kz, tsin,
				vc3, nc3, mode);

      if (!(master->conserve & NOINIT)) {
	memcpy(master->oldeta, master->d2, geom->szcS * sizeof(double));
	hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
				reset->tsnames, "eta", master->eta, tsin,
				vc2, nc2, (mode|VEL2D));
	memcpy(master->d2, master->eta, geom->szcS * sizeof(double));
      }
    }
    /* Read in additional specified variabes if required             */
    for (tn = 0; tn < master->ntrvars; tn++) {
      int trn = master->trvm[tn];
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, master->trvars[tn], 
			      master->tr_wc[trn], tsin, vc3, nc3, mode);
    }
    for (tn = 0; tn < master->ntrvarsS; tn++) {
      int trn = master->trvm[tn];
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, master->trvarsS[tn], 
			      master->tr_wcS[trn], tsin, vc2, nc2, (mode|VEL2D));
    }
    if (master->tmode & DO_SWR) {
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "swr", master->light, tsin,
			      vc2, nc2, (mode|VEL2D));
    }
    if (master->tmode & DO_DZ) {
      int nmode = (mode|U1GEN);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "dz", master->dz, tsin,
			      vc3, nc3, mode);
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "dzu1", master->dzu1, tsin,
			      ve3, ne3, nmode);
    }
    /* Read in temp and salt if these variables aren't advected      */
    if (!master->trinfo_3d[master->tno].advect && 
	!master->trinfo_3d[master->tno].diffuse)
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "temp", master->temp, tsin,
			      vc3, nc3, mode);
    if (!master->trinfo_3d[master->sno].advect && 
	!master->trinfo_3d[master->sno].diffuse)
      hd_trans_multifile_eval(master, reset->ntsfiles, reset->tsfiles,
			      reset->tsnames, "salt", master->sal, tsin,
			      vc3, nc3, mode);

    /* Set the next update time                                      */
    reset->tnext += reset->dt;
    tsout = reset->tnext;
  }
  return tsout;
}

/* END trans_reset_event()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the origin datum read from the file attributes            */
/*-------------------------------------------------------------------*/
int get_datum(int ntsfiles, timeseries_t **tsfiles,
	      cstring * names, char *varname, 
	      char *attname)
{
  int varids[MAXNUMTSFILES];
  char *text;
  if (hd_ts_multifile_get_index(ntsfiles, tsfiles, 
				names, varname, varids) == 0)
    hd_quit("hd_ts_multifile_get_att: Unable to evaluate the variable '%s'.\n", varname);

  if (ntsfiles > 0) {
    int i;
    int index = -1;
    timeseries_t *ts;
    datafile_t *df;
    df_variable_t *v;
    df_attribute_t *a;

    for (i = 0; i < ntsfiles; ++i) {
      ts = tsfiles[i];
      if (varids[i] < 0) continue;
      index = i;
    }
    ts = tsfiles[index];
    v = df_get_variable(ts->df, varids[index]);
    df = ts->df;
    for (i = 0; i < v->na; ++i) {
      df_attribute_t *a = &v->attributes[i];
      if (strcasecmp(a->name, attname) == 0) {
	text = ATT_TEXT(a);
	if (strcmp(text, "top right") == 0) 
	  return(TOPRIGHT);
	else if (strcmp(text, "bottom left") == 0) 
	  return(BOTMLEFT);
	else {
	  hd_warn("trans_reset_init: can't infer streamline origin datum.\n");
	  return(-1);
	}
      }
    }
  }
  return(-1);
}

/* END get_datum()                                                   */
/*-------------------------------------------------------------------*/
