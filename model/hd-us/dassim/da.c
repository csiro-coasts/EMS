/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/da.c
 *  
 *  Description:
 *  Event scheduler routines for data assimilation.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

/*
 * Some notes:
 *  o) The da event is initialised first before anything (eg. dumps)
 *     that may want to make use of it
 *  o) The da event fires off last, after anything that uses it. This
 *     way dumps are consistent whether it is a renanlysis cycle or
 *     not
 *  o) Additional logic has been added to the following in order to
 *     handle DA:
 *     - dumps (dump_init, dump_event)
 *     - timeseries output (ts_init, ts_event)
 *     - relaxtion event (tr_relax_event)
 *     - the function that calls do_ts_relax (auxiliary_routines)
 */

#include <math.h>
#include <string.h>
#include <libgen.h>
#include "hd.h"

#ifdef HAVE_DA
#include "da_interface.h"

typedef struct {
  char name[MAXSTRLEN];
  int tr_id;          /* tracer id */
  int tr_id_anl;      /* tracer id for the analysis */
  int tr_id_obs_val;  /* tracer id for the observation values */
  int tr_id_obs_err;  /* tracer id for the observation errors */
  int is2D;           /* 0 == 3D, 1 == 2D  */
  int offset;         /* Row offset in DA vectors */
  int size;           /* Length of this state */
} da_state;


typedef struct {
  double dt;                    /* Analysis time step               */
  double fcst_dt;               /* Forecast time step               */
  char rsfname[MAXSTRLEN];      /* Restart file name                */
  int rfn;                      /* Dump number of restart file      */
  int first;                    /* Flag for first iteration         */
  double prev;                  /* Previous time of event           */

  /* Pointer to interface with the da library */
  void *daobj;

  /* Array of da_states in the Anomaly matrix */
  int       nstates;
  da_state *states;
  
} da_data_t;

/* Local function declarations */
static void read_previous_state(master_t *master, da_data_t *data);
static void read_anom_file(da_data_t *data, master_t *master, char *fname);
static void init_obs(da_data_t *data,master_t *master,parameters_t *params);
static void init_maps(da_data_t *data, master_t *master);


/*------------------------------------------------------------------*/
/* Functions for reading the schedule the data assimilation.        */
/*------------------------------------------------------------------*/
int da_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  da_data_t *data = NULL;
  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  char *tags[MAXSTRLEN * MAXNUMARGS];
  int t, n, nf;
  int nrows;
  int offset;
  int nda = 0;

  /*
   * Do nothing if not doing da
   */
  if (!params->da) {
    master->da = NONE;
    return 0;
  }

  /*----------------------------------------------------------------*/
  /* Allocate memory                                                */
  data = (da_data_t *)malloc(sizeof(da_data_t));

  /*----------------------------------------------------------------*/
  /* Check restarts are set.                                        */
  if (strlen(params->restart_name) == 0 || master->restart_dt == 0.0)
    hd_quit("Restart files must be written for data assimilation.\n");
  else {
    strcpy(data->rsfname, dirname(strdup(master->restart_name)));
    strcat(data->rsfname, "/da_restart.nc");
  }

  /*----------------------------------------------------------------*/
  /* Read parameters                                                */
  prm_set_errfn(hd_silent_warn);

  /* Set the initial state of DA */
  master->da = 0;

  master->da_dt = data->dt = params->da_dt;

  if (fmod(master->da_dt, master->restart_dt) != 0)
    hd_quit("Data assimilation dt (%f) must be divisible by restart_dt (%f).\n",
	    master->da_dt, master->restart_dt);
  master->da_fcst_dt = data->fcst_dt = params->da_fcst_dt;
  if (master->da_dt > master->da_fcst_dt) {
    hd_warn("DA: Forecast cycle must be > reanalysis cycle. No forecast attempted.\n");
    master->da_fcst_dt = data->fcst_dt = 0.0;
  }
  data->first = 1;

  /* Dumpfile number for restarts */
  data->rfn = params->ndf;

  /* Create the da object */
  data->daobj = da_create_object(params->da_ls);

  /*
   * Parse the anomaly fields and create & fill in the da_states
   */
  nf = parseline(params->da_anom_states, fields, MAXNUMARGS);
  if (nf < 1)
    hd_quit("DA: Must have atleast one anomaly state\n");
  
  data->nstates = nf;
  data->states  = (da_state *)calloc(nf, sizeof(da_state));
  
  offset = 0;
  for (n = 0; n < nf; n++) {
    int ti, sz, is2D;
    char     *tr  = fields[n];
    da_state *das = data->states;
    int t;
    int tr_id_obs_val, tr_id_obs_err, tr_id_anl;

    /* Check if its 3D or 2D */
    is2D = 0;
    if ((ti = tracer_find_index(tr, master->ntr, master->trinfo_3d)) > -1) {
      sz = master->geom->b3_t;
    } else 
      if((ti = tracer_find_index(tr, master->ntrS, master->trinfo_2d)) > -1) {
	sz   = master->geom->b2_t;
	is2D = 1;
      } else if (strcmp(tr, "eta") == 0) {
	/* Special cases */
	sz   = master->geom->b2_t;
	is2D = 1;
      } else {
	hd_quit("DA: No tracer with the name '%s' found in the model\n", tr);
      }
    
    /*
     * Now loop through all the tracers and find its corresponding
     * analysis and observation ids
     */
    tr_id_anl = -1;
    tr_id_obs_val = -1;
    tr_id_obs_err = -1;
    for (t = 0; t < master->ntr; t++) {
      int ntags;
      if (strlen(master->trinfo_3d[t].tag)) {
	strcpy(buf, master->trinfo_3d[t].tag);
	ntags = parseline(buf, tags, MAXNUMARGS);
	/* 
	 * If tags is referring to this state 
	 * Note: There can be no other tags
	 */
	if (strcmp(tags[1], tr) == 0 || strcmp(tags[1], "eta") == 0) {
	  /* then see whether its an analysis or obs value/error */
	  if (strcmp(tags[0], "DA_ANL") == 0) {
	    tr_id_anl = t;
	    master->trinfo_3d[t].flag |= DA_ANL;
	    nda++;
	  } else if (strcmp(tags[0], "DA_OBS_VAL") == 0) {
	    tr_id_obs_val = t;
	    master->trinfo_3d[t].flag |= DA_OBS_VAL;
	    nda++;
	  } else if (strcmp(tags[0], "DA_OBS_ERR") == 0) {
	    tr_id_obs_err = t;
	    master->trinfo_3d[t].flag |= DA_OBS_ERR;
	    nda++;
	  }
	}
      }
    }

    /*
     * Fill in the da_states struct
     */
    sprintf(das[n].name, "%s", tr);
    das[n].tr_id = ti;
    das[n].tr_id_anl     = tr_id_anl;
    das[n].tr_id_obs_val = tr_id_obs_val;
    das[n].tr_id_obs_err = tr_id_obs_err;
    das[n].is2D   = is2D;
    das[n].size   = sz;
    das[n].offset = offset;

    /* Add this state to da obj */
    da_add_state(data->daobj, das[n].name, das[n].offset);

    /* Increment offset */
    offset += sz;

  } /* end for (n < nf) */

  /* We're now ready to Read from file and create the A vector */
  if (params->da_anom_file != NULL && strlen(params->da_anom_file))
    read_anom_file(data, master, params->da_anom_file);  
  else
    hd_quit("DA: Anomaly file not specified");

  /* Set up the maps */
  init_maps(data, master);

  /* Now initialise the observations */
  init_obs(data, master, params);

  schedSetPrivateData(event, data);
  event->next_event += data->dt;

  return 1;
}

/* END dassim_init()                                                */
/*------------------------------------------------------------------*/

/*
 * The main analysis function
 */
static void do_analysis(da_data_t *data, master_t *master)
{
  int n, offset = 0;
  int nobs;
  int cc, c, c2;
  int index = 0;
    
  /* Set the background fields */
  for (n=0; n<data->nstates; n++) {
    da_state *das = &data->states[n];
    for (cc = 1; cc <= das->size; cc++) {
      double val;
      if (das->is2D) {
	c  = geom->w2_t[cc];
	c2 = c;
	if (strcmp(das->name, "eta") == 0)
	  val = master->eta[c];
	else
	  val = master->tr_wcS[das->tr_id][c];
      } else {
	c  = geom->w3_t[cc];
	c2 = geom->m2d[c];
	val = master->tr_wc[das->tr_id][c];
      }

      /* Set the background vector */
      da_set_wb(data->daobj, index, val);

      /* Set the location vectors */
      da_set_xy(data->daobj, index, geom->cellx[c2], geom->celly[c2]);

      index++;
    }
  }
  
  /* Read all observations for this time */
  nobs = da_read_all_obs(data->daobj, master->t);
  hd_warn("DA: Found %d observations at %.4f days\n", nobs, master->t/86400.);

  /*
   * Call the main work horse function, if there are any observations
   * Otherwise, set the analysis the same as the background - no da
   */
  if (nobs > 0) {
    da_do_analysis(data->daobj);
  } else {
    da_copy_wb_wa(data->daobj);
  }
  
  /* Now copy stuff back */
  for (n=0; n<data->nstates; n++) {
    da_state *das = &data->states[n];
    /* 
     * Copy analysis fields back, if asked 
     */
    if (das->tr_id_anl > -1) {
      index = das->offset;
      for (cc = 1; cc <= das->size; cc++) {
	double *val_ptr;
	if (das->is2D) {
	  c = geom->w2_t[cc];
	  if (strcmp(das->name, "eta") == 0)
	    val_ptr = &master->eta[c];
	  else
	    val_ptr = &master->tr_wcS[das->tr_id_anl][c];
	} else {
	  c = geom->w3_t[cc];
	  val_ptr = &master->tr_wc[das->tr_id_anl][c];
	}
	/* Get the analysis value */
	*val_ptr = da_get_wa(data->daobj, index++);
      }
    }
    /* 
     * Copy observation value back, if asked 
     */
    if (das->tr_id_obs_val > -1) {
      index = das->offset;
      for (cc = 1; cc <= das->size; cc++) {
	double *val_ptr;
	if (das->is2D) {
	  c = geom->w2_t[cc];
	  if (strcmp(das->name, "eta") == 0)
	    val_ptr = &master->eta[c];
	  else
	    val_ptr = &master->tr_wcS[das->tr_id_obs_val][c];
	} else {
	  c = geom->w3_t[cc];
	  val_ptr = &master->tr_wc[das->tr_id_obs_val][c];
	}
	/* Get the obs value */
	*val_ptr = da_get_wo_val(data->daobj, index++);
      }
    }
    /* 
     * Copy observation errors back, if asked 
     */
    if (das->tr_id_obs_err > -1) {
      index = das->offset;
      for (cc = 1; cc <= das->size; cc++) {
	double *val_ptr;
	if (das->is2D) {
	  c = geom->w2_t[cc];
	  if (strcmp(das->name, "eta") == 0)
	    val_ptr = &master->eta[c];
	  else
	    val_ptr = &master->tr_wcS[das->tr_id_obs_err][c];
	} else {
	  c = geom->w3_t[cc];
	  val_ptr = &master->tr_wc[das->tr_id_obs_err][c];
	}
	/* Get the obs value */
	*val_ptr = da_get_wo_err(data->daobj, index++);
      }
    }
  }
}


/*------------------------------------------------------------------*/
/* Data assimilation event.                                         */
/* During the forecast cycle master->da = NO_DA                     */
/* During the reanalysis cycle master->da = DO_DA                   */
/* Analyses are made during the forecast cycle at da_dt intervals.  */
/* Restart dumps are made during the reanalysis cycle at da_dt.     */
/* Model state is reset to the restart file at fcst_dt intervals.   */
/*------------------------------------------------------------------*/
double da_event(sched_event_t *event, double t)
{
  master_t   *master = (master_t *)schedGetPublicData(event);
  da_data_t  *data   = (da_data_t *)schedGetPrivateData(event);
  geometry_t *geom   = master->geom;
  char restart_fname[MAXSTRLEN];
  struct timeval tm1, tm2;

  if (t >= (event->next_event - SEPS)) {
    
    /*--------------------------------------------------------------*/
    if (master->da & NO_DA) {  /* Forecast cycle                    */
      if (master->da & FCST_DA) {
	/* Reload the state from the previous DA restart            */
	read_previous_state(master, data);

	/* Set the next event code (for the reanalysis)             */
	event->next_event = data->prev;
	master->da = DO_DA;
      } else {
	/* 
	 * Perform the analysis
	 */
	gettimeofday(&tm1, NULL);
	do_analysis(data, master);
	gettimeofday(&tm2, NULL);
	hd_warn("DA analysis took %.5f secs\n", 
		((tm2.tv_sec  - tm1.tv_sec)+(tm2.tv_usec - tm1.tv_usec)*1e-6));
	
	if (data->fcst_dt) {
	  /* Setup any extended forecast */
	  master->da |= FCST_DA;
	  data->prev = event->next_event;
	  event->next_event += data->fcst_dt;
	} else {
	  /* Rewind and update da flag */
	  read_previous_state(master, data);
	  master->da = DO_DA;
	}
      }
    /*--------------------------------------------------------------*/
    } else {  /* Reanalysis cycle                                   */
      FILE *fp;
      char buf[MAXSTRLEN];
      
      /* Dump output file */
      // xxx should've already dumped
      /*
	dump_f_snapshot(sched_get_even_by_name(schedule, "dumps"),
	data->rfn, master->t);
      */

      strcpy(restart_fname, master->restart_name);
      fp = fopen(data->rsfname, "r");
      if (fp != NULL) {
	fclose(fp);
	unlink(data->rsfname);
      }
      // xxx copy file instead of renaming to facilitate crash recovery
      //rename(restart_fname, data->rsfname);
      sprintf(buf, "cp %s %s", restart_fname, data->rsfname);
      if (system(buf))
	hd_quit("DA error copying restart file %s to %s\n",
		restart_fname, data->rsfname);
      
      /* Set the next event (for the forecast cycle) */
      event->next_event += data->dt;
      master->da = NO_DA;
    }
  }

  return event->next_event;
}

/* END da_event()                                                   */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Clean up the data assimilation                                   */
/*------------------------------------------------------------------*/
void da_cleanup(sched_event_t *event, double t)
{
  da_data_t *data = (da_data_t *)schedGetPrivateData(event);

  if (data != NULL) {
    if (data->states != NULL)
      free(data->states);

    /* Call out to the da library for the objects cleanup */
    da_destroy_object(data->daobj);

    free(data);
  }
}

/* END da_cleanup()                                                 */
/*------------------------------------------------------------------*/

/*******************/
/* LOCAL Functions */
/*******************/

/*------------------------------------------------------------------*/
/* Reload the state from the previous DA restart                    */
/*------------------------------------------------------------------*/
static void read_previous_state(master_t *master, da_data_t *data) {
  int dafid = 0;                  /* input cdf file id              */
  char restart_fname[MAXSTRLEN];  /* Restart filename               */
  char buf[MAXSTRLEN];
  int i;
  dump_data_t *dumpdata = master->dumpdata;
  FILE *fp = NULL;

  /* Use the DA restart to initialize the reanalysis                */
  strcpy(restart_fname, data->rsfname);
  /* The first reanalysis is re-initialized with the input file     */
  if (data->first) {
    strcpy(restart_fname, master->idumpname);
    data->first = 0;
  }
  
  /* Open the file                                                  */
  if (nc_open(restart_fname, NC_NOWRITE, &dafid) != NC_NOERR)
    hd_quit("Can't find DA restart file %s\n", restart_fname);

  /* Reset the master state                                         */
  dump_re_read(master, dafid, 0);

  /* We're done with the file */
  nc_close(dafid);

  /* Reset the windows                                              */
  for (i=1; i<=master->nwindows; i++)
    window_reset(master, hd_data->window[i], hd_data->windat[i],
		 hd_data->wincon[i], RS_ALL);

  /* Get the new start time                                         */
  schedule->t = get_restart_time(restart_fname, schedule->units);

  /* Reset time for restart */
  dumpdata->dumplist[dumpdata->ndf-1].tout = schedule->t;
}

/* END read_previous_state()                                        */
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/* Creates the mapping objects                                      */
/*------------------------------------------------------------------*/
static void init_maps(da_data_t *data, master_t *master)
{
  geometry_t  *geom     = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int c,cc;
  int i,j,k;
  
  /* Initialise the maps object */
  da_create_maps(data->daobj, dumpdata->gridx, dumpdata->gridy, geom->layers,
		 geom->nce1, geom->nce2, geom->nz);

  /* 
   * Now fill in the cart-sparse maps 
   */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    /* 2D */
    c = geom->w2_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    da_fill_map(data->daobj, cc-1, i, j, -1);
  }
  /* 3D */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    da_fill_map(data->daobj, cc-1, i, j, k);
  }
}

/* END init_maps()                                                  */
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/* Read in the anomaly file and fill in the A matrix                */
/*------------------------------------------------------------------*/
static void read_anom_file(da_data_t *data, master_t *master, char *fname)
{
  geometry_t *geom = master->geom;
  datafile_t *df   = df_alloc();
  int n, i, col, nrows, row;
  int index;
  int row_offset;

  /* Initialise the datafile struct */
  df_read(fname, df);

  /* Let the library know which is the record variable */
  if ( (index = df_get_index(df, "nmember") ) > -1 ) {
    df_set_record(df, index);
  } else {
    hd_quit("DA: 'nmember' field not found in Anomaly file '%s'\n", fname);
  }

  /* Count the number or rows */
  nrows = 0;
  for (n=0; n<data->nstates; n++)
    nrows += data->states[n].size;
  
  /* Allocate memory for the A matrix */
  if (da_allocA(data->daobj, nrows, df->nrecords))
    hd_quit("DA: Error in A matrix allocation");

  /* Fill in the A matrix */
  row_offset = 0;
  for (n=0; n<data->nstates; n++) {
    da_state *das    = &data->states[n];
    int      varid   = df_get_index(df, das->name);
    df_variable_t *v = &df->variables[varid];
    int c,cc;
    double val;
    int coords[3];

    /* Read the entire file in one hit */
    df_read_all_records(df, v);

    /*
     * Note: We could use the x,y coordinates instead so that
     * anomaly files on a different grid can be interpolated onto the
     * model grid (eg. Setas's onto Storm) but the libraries need
     * updating for this to happen. Project for another day.
     */

    /* Over every column, i.e. member */
    for (col = 0; col<df->nrecords; col++) {
      row = row_offset;
      for (cc = 1; cc <= das->size; cc++) {
	if (das->is2D) {
	  c = geom->w2_t[cc];
	  coords[0] = geom->s2j[c];
	  coords[1] = geom->s2i[c];
	} else {
	  c  = geom->w3_t[cc];
	  coords[0] = geom->s2k[c];
	  coords[1] = geom->s2j[c];
	  coords[2] = geom->s2i[c];
	}
	/* Get the value */
	val = df_get_data_value(df, v, col, &coords[0]);
	
	/* Fill in A */
	da_fill_A(data->daobj, row++, col, val);
      }
    }
    row_offset = row;
  }

  /* 
   * useful debugging
   * da_writeA(data->daobj, "anm.txt");
   */
  
  /* 
   * Sanity check to make sure we've initialised every single element
   */
  if (!da_check_A(data->daobj))
    hd_quit("DA: Error in initialising A\n");

  /* All done we can close the file and finish up */
  df_free(df);
}

/* END read_anom_file */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Read in the observations' info and initialise the da_obs objects */
/*------------------------------------------------------------------*/
static void init_obs(da_data_t *data, master_t *master, parameters_t *params)
{
  int i;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  char buf[MAXSTRLEN];
    
  /* s/c code */
  if (params->da_nobs <= 0)
    return;

  /*
   * Now loop over and create the observation objects
   */
  for(i=0; i<params->da_nobs; i++) {
    char *type = params->da_obs_types[i];
    char *name = params->da_obs_names[i];;
    char *file = params->da_obs_files[i];
    char *loc  = params->da_obs_locs[i];
    double dt  = params->da_obs_dt[i];
    double err = params->da_obs_errs[i];

    void *obs   = NULL;
    char *fname = fv_get_filename(file, buf);

    int nf,n;

    /*
     * Check name
     */
    if (name == NULL || strlen(name) == 0) {
      hd_warn("DA: No name supplied for DA_OBS%1.1d\n", i);
    }

    /* Switch on observation type */
    if (strcasecmp(type, "fixed") == 0) {
      double pos[3];

      /* Decode location */
      nf = parseline(loc, fields, MAXNUMARGS);
      if (nf != 3)
	hd_quit("DA: Incorrect location specified for '%s'\n", name);

      /*
       * Create the fixed observation object
       */
      obs = da_create_fixed_obs(name, fname, master->timeunit, err, dt);

      if (!str2double(fields[0], &pos[0]))
	hd_quit("DA: Incorrect longitude '%s' specified for '%s'\n", fields[0], name);
      if (!str2double(fields[1], &pos[1]))
	hd_quit("DA: Incorrect latitude '%s' specified for '%s'\n", fields[1], name);
      if (!str2double(fields[2], &pos[2]))
	hd_quit("DA: Incorrect depth '%s' specified for '%s'\n", fields[2], name);

      /* Set location */
      da_fixed_obs_set_location(obs, data->daobj, pos[0], pos[1], pos[2]);

    } else if (strcasecmp(type, "glider") == 0) {
      /*
       * Create the glider observation object
       */
      obs = da_create_glider_obs(name, fname, master->timeunit, err, dt);
      
      /* Decode location */
      nf = parseline(loc, fields, MAXNUMARGS);
      if (nf != 3)
	hd_quit("DA: Incorrect location specified for '%s'\n", name);
      
      /* Pass in the strings */
      da_glider_obs_set_location(obs, fields[0], fields[1], fields[2]);

    } else if (strcasecmp(type, "2d-field") == 0) {
      /*
       * Create the 2D field observation object
       */
      obs = da_create_2d_field_obs(name, fname, master->timeunit, err);

    } else {
      hd_quit("DA: Unknown observation type '%s' specified for DA_OBS%1.1d\n", 
	                                                      type, i);
    }

    /* Loop over to find all the variables in this file */
    for (n=0; n<data->nstates; n++) {
      char *name = data->states[n].name;
      char *var  = fv_get_varname(file, name, buf);
      int  is3D  = strcmp(name, "eta") == 0 ? 0 : 1;
      /* Attach this state to name given in the prm file */
      da_obs_add_var(data->daobj, obs, var, name, is3D);
    }

    /* Add to the main object */
    da_add_obs(data->daobj, obs);
  }
}


/***********************/
/* END local functions */
/***********************/

/*
 * Data assimilation main exported functions
 */
void dassim_init(master_t *master)
{
  sched_register(schedule, "dassim", da_init,
                 da_event, da_cleanup, master, NULL, NULL);
}

/*
 * This is needed to initialise DA to its proper state. If DA is off,
 * then this flag will already be non-zero, specifically NONE
 */
void dassim_post_init(master_t *master)
{
  if (master->da == 0)
    master->da = NO_DA;
}

void dassim_end(void)
{
  sched_deregister(schedule, "dassim");
}

void dassim_run_setup(FILE *fp, sched_event_t *event)
{
  master_t     *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  da_data_t    *data   = (da_data_t *)schedGetPrivateData(event);
  int n,r, nobs;

  fprintf(fp,"  Analysis cycle DT = %6.3f hours\n", data->dt/3600);
  if (data->fcst_dt)
    fprintf(fp,"  Forecast cycle DT = %6.3f hours\n", data->fcst_dt/3600);

  fprintf(fp, "\n  Number of analysis states = %d\n", data->nstates);
  for (n=0; n<data->nstates; n++) {
    da_state *das = &data->states[n];
    tracer_info_t *trinfo = (das->is2D) ? master->trinfo_2d : master->trinfo_3d;
    char *anl     = (das->tr_id_anl < 0)     ? "none" : trinfo[das->tr_id_anl].name;
    char *obs_val = (das->tr_id_obs_val < 0) ? "none" : trinfo[das->tr_id_obs_val].name;
    char *obs_err = (das->tr_id_obs_err < 0) ? "none" : trinfo[das->tr_id_obs_err].name;

    int isrlx = 0;
    /* Loop to find which tracer this is being relaxed to */
    for (r=0; r<master->nrlx; r++) {
      if (master->relax[r] == das->tr_id) {
	/* 
	 * Now compare the relax filename with the name of the
	 * analysis tracer. Note: this is not the actual relax_info
	 * from the scheduler
	 */
	if (strcmp(master->trinfo_3d[das->tr_id].relax_file, anl) == 0) {
	  isrlx = 1;
	  break;
	}
      }
    }
    fprintf(fp, "   State #%2d (%s) :\n", n, das->name);
    fprintf(fp, "    Tracer for analysis    values = %s\n", anl);
    fprintf(fp, "    Tracer for observation values = %s\n", obs_val);
    fprintf(fp, "    Tracer for observation errors = %s\n", obs_err);
    if (!isrlx)
      fprintf(fp,"    WARNING : STATE NOT BEING RELAXED\n"); 
  }
  
  /* Name of anomaly file */
  fprintf(fp, "\n   Name of Anomaly file = %s\n", params->da_anom_file);
  
  /* Localisation factor */
  if (params->da_ls > 0.0)
    fprintf(fp, "\n   Using localisation factor of %.2f degrees\n", params->da_ls);
  else
    fprintf(fp, "\n   Localisation not invoked\n");

  /* Names of all obervations */
  fprintf(fp, "\n");
  nobs = da_get_num_obs(data->daobj);
  fprintf(fp, "   Number of Observation files = %d\n", nobs);
  for (n=0; n<nobs;n++) {
    void *obs = da_get_obs_ptr(data->daobj, n);
    double dt = da_obs_get_dt(obs);
    fprintf(fp, "    Obs #%2d :\n", n);
    fprintf(fp, "     Name = %s\n", da_obs_get_name(obs));
    fprintf(fp, "     File = %s\n", da_obs_get_fname(obs));
    fprintf(fp, "     Vars = %s\n", da_obs_get_vars(obs));
    fprintf(fp, "     Type = %s\n", da_obs_get_type(obs));
    fprintf(fp, "     Error = %g\n", da_obs_get_error(obs));
    if (dt > 0)
      fprintf(fp, "     DT    = %g secs\n", dt);
    else
      fprintf(fp, "     DT    = not time averaged\n");
    fprintf(fp, "     %s\n", da_obs_get_location(obs));
  }
  fprintf(fp, "\n");
}
#endif /* HAVE_DA */

/* EOF */
