/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/outputs/dumpfile.c
 *  
 *  Description:
 *  Routines for meco model to allow handling
 *  multiple netCDF output files, with different
 *  sets of time dependent variables and different
 *  output schedules.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dumpfile.c 6341 2019-09-16 23:57:49Z her127 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include <libgen.h>
#include "hd.h"
#include "tracer.h"

#define MAXDUMPFILES (300)
#define TYPE_FROM_FILE (-9999)
#define SIMPLE_MISSING_VALUE 1e35
#define SIMPLE_SHORT_MISSING_VALUE -32767
#define STD_ALL_VARS "u1av u2av wtop topz eta wind1 wind2 patm u1 u2 w u v dens dens_0 Kz Vz u1bot u2bot Cd flag u1vh u2vh "
#define SIMPLE_ALL_VARS "uav vav avg_speed avg_dir u v current_speed current_dir w eta wind_u wind_v wind_mag wind_dir patm dens dens_0 Kz Vz bottom_u bottom_v bottom_speed bottom_dir Cd "

#define PARRAY_MISSING_VALUE 1e35
#define PARRAY_SHORT_MISSING_VALUE -32767
#define PARRAY_ALL_VARS "u1av u2av wtop topz eta wind1 wind2 patm u1 u2 w dens dens_0 Kz Vz u1bot u2bot Cd u v uav vav "
#define PARRAY_EXCLUDE_VARS "botz"
#define DF_PARRAY    1
#define DF_MEMORY    2

#define WILDCARD '*'

/* Standard SHOC dumpfile format. */
/* Structure to describe each dump file time dep variable */
typedef struct {
  void **v;                     /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  nc_type type;                 /* netCDF type of this variable */
  int xylocation;
  int zlocation;
  int sediment;                 /* Flag; = 1 if variable is in sediment */
} df_std_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_std_var_t *vars;           /* List of dump file variables */
} df_std_data_t;

/*Structure for data when dispatching event */

typedef struct {
  double t;                      /* Netcdf file id */
  void* data ;                  /* Netcdf record number */
} df_dispatch_data_t;


/* Prototypes for routines below */
static int read_subsection_range(FILE * fp, char *key, int nc, int nf,
         int *sslower, int *ssnc, int *ssnf);
         
static void read_dumpfile_points(FILE * fp, char *key, dump_file_t* file);
         
static void vector_component(dump_data_t *dumpdata, double **vi,
                             double **vj, int ilower, int jlower, int nce1,
                             int nce2, int mode, double **nv);
void write_text_att(int cdfid, int varid, const char *name,
                    const char *text);
int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref);
void df_null_reset(dump_data_t *dumpdata, dump_file_t *df, double t);

/* History output file prototypes */
static void *df_restart_create(dump_data_t *dumpdata, dump_file_t *df);
static void df_restart_write(dump_data_t *dumpdata, dump_file_t *df, double t);
static void df_restart_close(dump_data_t *dumpdata, dump_file_t *df);

/* Standard output file prototypes */
static void df_std_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_std_init_data(dump_data_t *dumpdata, dump_file_t *df,
                             int fid);
static void df_std_get_dimids(int cdfid, df_std_var_t var, int *e1id,
                              int *e2id, int *zid);
static void df_std_get_dimsizes(dump_file_t *df, df_std_var_t var,
                                size_t * ne1, size_t * ne2, size_t * nz);

static void *df_simple_create(dump_data_t *dumpdata, dump_file_t *df);
static void *df_simple_cf_create(dump_data_t *dumpdata, dump_file_t *df);
static void df_simple_write(dump_data_t *dumpdata, dump_file_t *df,
                            double t);
static void df_simple_close(dump_data_t *dumpdata, dump_file_t *df);
static void df_simple_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_simple_init_data(dump_data_t *dumpdata, dump_file_t *df,
                                int fid);
static void df_simple_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, size_t * start, size_t * count,
                                  double **values);
static void df_simple_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, size_t * start, size_t * count,
                                  double ***values, int klower);

static int decode_args(char *oprm, char *args[], int maxn);
static int contains_string(char *p, char *tag);
int get_transfiles(dump_data_t *dumpdata, FILE * fp);
int set_transfiles(int fn, char *key, dump_data_t *dumpdata, dump_file_t *list, 
		   double t, FILE *fp);
void read_dumpfiles(FILE *fp, dump_file_t *list, dump_data_t *dumpdata,
		    double t, int f);
void read_filter(dump_file_t *df, FILE *fp, char *key);
double df_filter_ij(dump_file_t *df, double **a, unsigned long **flag, 
		    int ne1, int ne2, int i, int j);
double df_filter(dump_file_t *df, double *a, int c);
void df_filter_1d(dump_file_t *df, dump_data_t *dumpdata, double ***a, int i, int j);
void df_filter_2d(dump_file_t *df, dump_data_t *dumpdata, double **a, 
		  int is, int ie, int js, int je);
void df_filter_3d(dump_file_t *df, dump_data_t *dumpdata, double ***a, 
		  int is, int ie, int js, int je, int ks, int ke);

void parray_grid_init(dump_data_t *dumpdata, dump_file_t *df);
void parray_gride_init(dump_data_t *dumpdata, dump_file_t *df);

/* Memory output prototypes */
void *df_memory_create(dump_data_t *dumpdata, dump_file_t *df);
void df_memory_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_memory_close(dump_data_t *dumpdata, dump_file_t *df);

/* ROMS output as defined in df_roms.c */
void *df_roms_create(dump_data_t *dumpdata, dump_file_t *df);
void df_roms_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_roms_close(dump_data_t *dumpdata, dump_file_t *df);
void df_roms_reset(dump_data_t *dumpdata, dump_file_t *df, double t);

/* ROMS boundary output as defined in df_roms.c */
void *df_roms_create_bdry(dump_data_t *dumpdata, dump_file_t *df);
void df_roms_write_bdry(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_roms_reset_bdry(dump_data_t *dumpdata, dump_file_t *df, double t);

static struct {
  char *tag;
  void *(*create) (dump_data_t *dumpdata, dump_file_t *);
  void (*write) (dump_data_t *dumpdata, dump_file_t *df, double t);
  void (*close) (dump_data_t *dumpdata, dump_file_t *df);
  void (*reset) (dump_data_t *dumpdata, dump_file_t *df, double t);
  int xyGridOutput;
} df_func_map[] = {
  {
    "standard", df_std_create, df_std_write, df_std_close, df_std_reset, 1}, {
    "simple", df_simple_create, df_simple_write, df_simple_close, df_simple_reset, 1}, {
    "simple_cf", df_simple_cf_create, df_simple_write, df_simple_close, df_simple_reset, 1}, {
    "parray", df_parray_create, df_parray_write, df_parray_close, df_parray_reset, 0}, {
    "memory", df_memory_create, df_memory_write, df_memory_close, df_null_reset, 0}, {
    "sparse", df_ugrid3_create, df_ugrid3_write, df_ugrid_close, df_ugrid_reset, 1}, {
    "ugrid", df_ugrid_create, df_ugrid_write, df_ugrid_close, df_ugrid_reset, 1}, {
    "mom", df_mom_create, df_mom_write, df_mom_close, df_mom_reset, 1}, {
    "roms", df_roms_create, df_roms_write, df_roms_close, df_roms_reset, 1}, {
    "roms-bdry", df_roms_create_bdry, df_roms_write_bdry, df_roms_close, df_roms_reset_bdry, 1}, {
    "restart", df_restart_create, df_restart_write, df_restart_close, df_ugrid_reset, 1}, {
    NULL, NULL, NULL, NULL, NULL}
};

int barof = 0;  /* Global flag for standard OBCs */

void df_null_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  return;
}

/*------------------------------------------------------------------*/
/* Time scheduler_t functions for outputing the model results       */
/*------------------------------------------------------------------*/
double dump_event(sched_event_t *event, double t)
{
  int i;
  double tout = schedule->stop_time;
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  df_dispatch_data_t* dispatch_data = malloc(sizeof(df_dispatch_data_t));
  int debug = 0;

  /* Output dump and test point values if required */
  master_fill(master, window, windat, wincon);

  /*
   * Wait for dumps to finish before filling buffer
   */
#ifdef HAVE_PTHREADS
  if (event->sm_dump) {
    if (debug) {
      while (sem_trywait(event->sm_dump)) {
	if (errno == EAGAIN)
	  hd_warn("dump_event: waiting at %f\n", t);
	else
	  hd_quit("dump_event: bad errno\n");
	sleep(1);
      }
    } else {
      /* Block until buffer is filled */
      sem_wait(event->sm_dump);
    }
  }
#endif

  dumpdata_fill(geom, master, dumpdata);
  dispatch_data->t = t;
  schedSetPrivateData(event,dispatch_data);

  /* Run through the list of dumpfiles and find the the one that        */
  /* next fires off - after the current one
   * the function dumpfiles is called in the dispatch part of the event,
   * that means after the next event time is established below          */
  for (i = 0; i < dumpdata->ndf; ++i) {
    int mon;
    /*
     * Setup the increments
     */
    if (dumpdata->dumplist[i].incf == YEARLY)
      dumpdata->dumplist[i].tinc = next_year(t, master->timeunit);
    if (dumpdata->dumplist[i].incf == SEASONAL)
      dumpdata->dumplist[i].tinc = next_season(t, master->timeunit, &mon);
    if (dumpdata->dumplist[i].incf == MONTHLY)
      dumpdata->dumplist[i].tinc = next_month(t, master->timeunit, &mon);
    if (dumpdata->dumplist[i].incf == DAILY)
      dumpdata->dumplist[i].tinc = next_day(t, master->timeunit, &mon);

    /* Skip dump if the DA flag doesn't match for this dumpfile */
    if (master->da && !(dumpdata->dumplist[i].da_cycle & master->da)) continue;

    /*
     * The first if statement is needed as at this stage we may not
     * have updated dumpdata.tout as that happens in dumpfiles via the
     * dispatcher.
     */
    if(abs(dumpdata->dumplist[i].tout - t) < DT_EPS) {
      if (dumpdata->dumplist[i].tout + dumpdata->dumplist[i].tinc < tout &&
          t <= (dumpdata->dumplist[i].tstop + DT_EPS))
	tout = dumpdata->dumplist[i].tout+dumpdata->dumplist[i].tinc;
    } else {
      /* Here we're just looking for the next one that fires */
      if (dumpdata->dumplist[i].tout < tout &&
          t <= (dumpdata->dumplist[i].tstop + DT_EPS))
	tout = dumpdata->dumplist[i].tout;
    }
  }

#ifdef HAVE_PTHREADS
  if (event->sm_fill) sem_post(event->sm_fill);
#endif

  emstag(LDEBUG,"hd:output:dumpfile:dump_event","Next dump event at %f days.\n", tout / 86400.0);

  event->next_event = tout;
  return event->next_event;
}

/* END dump_event()                                                 */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Time scheduler_t functions for outputing the model results       */
/* Usage:                                                           */
/*dump_snapshot(sched_get_even_by_name(schedule, "dumps"),windat->t)*/
/*------------------------------------------------------------------*/
void dump_snapshot(sched_event_t *event, double t)
{
  int f;
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  /* Output dump and test point values if required */
  master_fill(master, window, windat, wincon);
  dumpdata_fill(geom, master, dumpdata);
  master->df_diagn_set = dumpfiles(dumpdata, t, master->prmfd);

  /* Loop through dump file dumplist, writing all files */
  for (f = 0; f < dumpdata->ndf; f++) {

    if (!(dumpdata->dumplist[f].finished)) {

      emstag(LTRACE,"hd:output:dumpfile:dump_snapshot","Writing dump %f days to '%s'\n", t / 86400.0,
             dumpdata->dumplist[f].name);

      if (dumpdata->us_type & US_IJ)
	dumpdata->dumplist[f].landfill(dumpdata);
      dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
      dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], t);
    }
  }
}

void dump_mp_snapshot(master_t *master)
{
  int f;
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  /* Output dump and test point values if required */
  master_fill(master, window, windat, wincon);
  dumpdata_fill(geom, master, dumpdata);

  /* Loop through dump file dumplist, writing all files */
  for (f = 0; f < dumpdata->ndf; f++) {
    if (strcmp(dumpdata->dumplist[f].type, "memory") != 0) continue;
    if (!(dumpdata->dumplist[f].finished)) {
      if (dumpdata->us_type & US_IJ)
	dumpdata->dumplist[f].landfill(dumpdata);
      dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
      dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], 
				  master->t);
    }
  }
}

/* END dump_snapshot()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Time scheduler_t functions for outputing the model results to    */
/* dumpfile number f.                                               */
/* Usage:                                                           */
/* dump_snapshot(sched_get_even_by_name(schedule, "dumps"), 1,      */
/*   windat->t)                                                     */
/*------------------------------------------------------------------*/
void dump_f_snapshot(sched_event_t *event, int f, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  /* Output dump and test point values if required */
  master_fill(master, window, windat, wincon);
  dumpdata_fill(geom, master, dumpdata);

  /* Loop through dump file dumplist, writing all files */
  if (!(dumpdata->dumplist[f].finished)) {

    emstag(LTRACE,"hd:output:dumpfile:dump_snapshot","Writing dump %f days to '%s'\n", t / 86400.0,
	   dumpdata->dumplist[f].name);

    if (dumpdata->us_type & US_IJ)
      dumpdata->dumplist[f].landfill(dumpdata);
    dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
    dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], t);
  }
}

/* END dump_f_snapshot()                                            */
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/* Output the model results to for dumpfiles of type 'memory'.      */
/*------------------------------------------------------------------*/
void dump_m_snapshot(master_t *master)
{
  int f;
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  /* Output dump and test point values if required */
  /*
  master_fill(master, window, windat, wincon);
  dumpdata_fill(geom, master, dumpdata);
  */
  /* Loop through dump file dumplist, writing all files */
  for (f = 0; f < dumpdata->ndf; f++) {
    if (strcmp(dumpdata->dumplist[f].type, "memory") != 0) continue;
    if (!(dumpdata->dumplist[f].finished)) {
      if (dumpdata->us_type & US_IJ)
	dumpdata->dumplist[f].landfill(dumpdata);
      dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
      dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], 
				  master->t);
    }
  }
}

/* END dump_m_snapshot()                                            */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Output the model results to for dumpfiles of type 'memory'.      */
/*------------------------------------------------------------------*/
void dump_eta_snapshot(master_t *master,    /* Master data          */
		       geometry_t **window, /* Slave geometry       */
		       window_t **windat,
		       win_priv_t **wincon
		       )
{
  int n, m, bn, cc, c, gc, i, f;
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  dump_file_t *df;
  int doeta = 1;
  int dov2d = 1;

  /* Loop through dump file dumplist, writing all files */
  for (f = 0; f < dumpdata->ndf; f++) {
    df = &dumpdata->dumplist[f];

    if (df->flag & DF_MPK && df->flag & (DF_ETA|DF_V2D)) {
      /*if (master->t < df->tout /10) continue;*/
      /* Transfer to master */
      if (master->nwindows > 1 && (doeta || dov2d)) {
	for (n = 1; n <= master->nwindows; n++) {
	  for (bn = 0; bn < window[n]->nobc; bn++) {
	    open_bdrys_t *open = window[n]->open[bn];
	    if (df->flag & DF_ETA) {
	      for (cc = 1; cc <= open->no2_t; cc++) {
		c = open->obc_t[cc];
		gc = window[n]->wsa[c];
		master->eta[gc] = windat[n]->eta[c];
	      }
	      for (cc = 1; cc <= open->no3_t; cc++) {
		c = open->obc_t[cc];
		gc = window[n]->wsa[c];
		master->temp[gc] = windat[n]->temp[c];
	      }
	    }
	    if (df->flag & DF_V2D) {
	      for (cc = 1; cc <= open->no2_e1; cc++) {
		c = open->obc_e1[cc];
		gc = window[n]->wsa[c];
		master->u1av[gc] = windat[n]->u1av[c];
	      }
	      for (cc = open->no3_e1 + 1; cc <= open->to2_e1; cc++) {
		c = open->obc_e1[cc];
		gc = window[n]->wsa[c];
		master->u1av[gc] = windat[n]->u1av[c];
	      }
	      for (cc = 1; cc <= open->to3_e1; cc++) {
		c = open->obc_e1[cc];
		gc = window[n]->wsa[c];
		master->u1[gc] = windat[n]->u1[c];
	      }
	    }
	  }
	}
      }
      /* Transfer to dumpdata */
      if (doeta && df->flag & DF_ETA) {
	s2c_2d(geom, master->eta, dumpdata->eta, geom->nce1, geom->nce2);
	s2c_3d(geom, master->temp, dumpdata->tr_wc[master->tno], geom->nce1, geom->nce2, geom->nz);
	doeta = 0;
      }
      if (dov2d && df->flag & DF_V2D) {
	s2c_2d(geom, master->u1av, dumpdata->u1av, geom->nfe1, geom->nce2);
	s2c_3d(geom, master->u1, dumpdata->u1, geom->nfe1, geom->nce2, geom->nz);
	dov2d = 0;
      }
      /* Dump */
      if (!(df->finished)) {
	if (dumpdata->us_type & US_IJ)
	  df->landfill(dumpdata);
	dump_bathy_mask(dumpdata, df->bathymask);
	df->write(dumpdata, df, master->t);
      }
    }
  }
}

/* END dump_eta_snapshot()                                            */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Set up the output dumps                                          */
/*------------------------------------------------------------------*/
void dumpfile_setup(master_t *master)
{
  dump_data_t *dumpdata = master->dumpdata;

  /* Get dump file parameters, create files and write model geometry */
  dumpfile_init(dumpdata, dumpdata->t, master->prmfd, &dumpdata->ndf,
                &dumpdata->dumplist);
}

/* END dumpfile_setup()                                             */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialise the output dumps                                      */
/*------------------------------------------------------------------*/
int dump_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  dump_data_t *dumpdata = master->dumpdata;
  char buf[MAXSTRLEN];

#ifdef HAVE_PTHREADS
  if (prm_read_char(master->prmfd, "SCHED_MODE", buf))
    if (strcasecmp(buf, "pthreads") == 0) {
      /* Allocate and initialise semaphores */
      event->sm_fill = malloc(sizeof(sem_t));
      event->sm_dump = malloc(sizeof(sem_t));

      if (sem_init(event->sm_fill, 0, 0))
	hd_quit("dump_init: error initialising semaphore fill\n");
      if (sem_init(event->sm_dump, 0, 1))
	hd_quit("dump_init: error initialising semaphore dump\n");
    }
#endif

  return 1;
}

/* END dump_init()                                                  */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Reads the dump tinc for sparse files and copies to grid_dt. Only */
/* resets grid_dt if the filetype = sparse, filename = TRANS_DATA   */
/* filename and a valid tinc is defined. Also sets dumpfiles = 0.   */
/*------------------------------------------------------------------*/
void trans_check_dump(master_t *master, dump_data_t *dumpdata, char *trdata)
{
  FILE *fp = master->prmfd;
  int i, ndf;
  double tinc;
  char key[MAXSTRLEN], buf[MAXSTRLEN], fname[MAXSTRLEN];

  dumpdata->ndf = 0;
  if (prm_read_int(fp, "OutputFiles", &ndf)) {
    for (i = 0; i < ndf; i++) {
      sprintf(key, "file%d.tinc", i);
      prm_skip_to_end_of_key(fp, key);
      if (fgets(buf, MAXSTRLEN, fp) != NULL && tm_scale_to_secs(buf, &tinc)) {
	sprintf(key, "file%d.filetype", i);
	prm_read_char(fp, key, buf);
	if (strcmp(buf, "sparse") == 0) {
	  sprintf(key, "file%d.name", i);
	  prm_read_char(fp, key, buf);
	  if (strcmp(buf, trdata) == 0)
	    master->grid_dt = tinc;
	}
      }
    }
  }
}

/*------------------------------------------------------------------*/



/*------------------------------------------------------------------*/
/* Initialise the output dumps                                      */
/*------------------------------------------------------------------*/
void* dump_dispatch(void* data)
{
  sched_event_t* event = (sched_event_t*) data;
  master_t *master = (master_t *)schedGetPublicData(event);
  df_dispatch_data_t* df_data = (df_dispatch_data_t*) schedGetPrivateData(event);
  int debug = 0;

#ifdef HAVE_PTHREADS
  if (event->sm_fill) {
    if (debug) {
      while (sem_trywait(event->sm_fill)) {
	if (errno == EAGAIN)
	  hd_warn("dump_dispatch: waiting at %f\n", df_data->t);
	else
	  hd_quit("dump_dispatch: bad errno\n");
	sleep(1);
      }
    } else {
      /* Block until buffer is filled */
      sem_wait(event->sm_fill);
    }
  }
#endif
  
  master->df_diagn_set = dumpfiles(master->dumpdata, df_data->t, master->prmfd);

#ifdef HAVE_PTHREADS
  /* Done, signal dumped */
  if (event->sm_dump) sem_post(event->sm_dump);
#endif

  return 0;
}


/*------------------------------------------------------------------*/
/* Initialise the output dumps                                      */
/*------------------------------------------------------------------*/
int dump_progress(sched_event_t *event)
{
  /*
   * Synchronisation is now achieved via a pari of semaphores so this
   * function is no longer used. It's possible we'll need this in the
   * future if we move to multi dispatch threads
   */
  hd_quit("dumpfile:dump_progress called when it should not have!\n");
  return 0;
}


/*------------------------------------------------------------------*/
/* Close the output dumps                                           */
/*------------------------------------------------------------------*/
void dump_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;

  /* Fill the dumpdata structure with master data */
  master_fill(master, window, windat, wincon);
  dumpdata_fill(geom, master, dumpdata);

  dumpfile_close(dumpdata, dumpdata->t, master->prmfd);
}

/* END dump_cleanup()                                               */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Routine to fill the master with dump data.                       */
/*------------------------------------------------------------------*/
void master_fill_dump(master_t *master, double t) 
{
  int i, j, k, c, n, f;
  int is, ie, js, je, ks, ke;
  dump_file_t *df;
  dump_data_t *dumpdata = master->dumpdata;
  geometry_t *geom = master->geom;

  /* Loop through dump file dumplist, writing any file which needs it */
  for (f = 0; f < dumpdata->ndf; f++) {
    /* Write to file if time is between start and stop time, and if file
       is open */
    /*
    if (!(dumpdata->dumplist[f].da_cycle & dumpdata->master->da)) continue;

    if (t >= (dumpdata->dumplist[f].tout - DT_EPS) &&
        t <= dumpdata->dumplist[f].tstop &&
        !(dumpdata->dumplist[f].finished)) {
      df = &dumpdata->dumplist[f];
      for (n = 0; n < df->nvars; n++) {
	df_std_var_t vn = df->vars[n];
	ks = 0; ke = geom->nz - 1;
	if (vn.ndims == 2) {
	  start[1] = df->jlower;
	  start[2] = df->ilower;
	  df_std_get_dimsizes(df, vn, &count[2], &count[1], NULL);
	  nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count,
			   (*((double ***)p)));
	} else if (vn.ndims == 3) {
	  start[1] = (!vn.sediment) ? df->klower : 0;
	  start[2] = df->jlower;
	  start[3] = df->ilower;
	  df_std_get_dimsizes(df, vn, &count[3], &count[2], &count[1]);
	}

      }
    }
    */
  }
}

/* END master_fill_dump()                                           */
/*------------------------------------------------------------------*/


int dumpfiles(dump_data_t *dumpdata, double t, FILE * fp)
{
  int f;
  int i, j, k, n, wn;
  int init_diagn = 0;
  int ymd[3];
  double newt;

  /* Loop through dump file dumplist, writing any file which needs it */
  for (f = 0; f < dumpdata->ndf; f++) {
    /* Skip dump if the DA flag doesn't match for this dumpfile */
    if (dumpdata->master->da && 
	!(dumpdata->dumplist[f].da_cycle & dumpdata->master->da)) continue;

    /* 
     * Write to file if time is between start and stop time, and if file
     * is open 
     */
    if (t >= (dumpdata->dumplist[f].tout - DT_EPS) &&
        t <= dumpdata->dumplist[f].tstop &&
        !(dumpdata->dumplist[f].finished)) {

      emstag(LDEBUG,"hd:output:dumpfile:dumpfiles", "Writing dump %f days to '%s' of type '%s'", t / 86400.0,
             dumpdata->dumplist[f].name,dumpdata->dumplist[f].type);

      wn = 0;
      /*
         for(wn=1; wn<=window[1]->nwindows; wn++) { if( (f == 0) &&
         VERBOSE(window[wn], "mass"))
         print_tracer_mass(window[wn],window[wn]->windat,window[wn]->wincon);
         } */

      /*
       * See if we need to start the next chunk 
       * Key off the timeunit for dumping rather than master timeunit
       */
      newt = t;
      tm_change_time_units(dumpdata->timeunit, dumpdata->dumplist[f].tunit, &newt, 1);
      tm_time_to_ymd(newt, dumpdata->dumplist[f].tunit, &ymd[0], &ymd[1], &ymd[2]);

      if (dumpdata->dumplist[f].chunk != NONE) {
	int chunk_next = 0;
	if (dumpdata->dumplist[f].chunk == YEARLY) {
	  /* Check Year only */
	  if (dumpdata->dumplist[f].curr_ymd[0] != ymd[0])
	    chunk_next = 1;
	} else if (dumpdata->dumplist[f].chunk == MONTHLY) {
	  /* Check Year and Month */
	  if (dumpdata->dumplist[f].curr_ymd[0] != ymd[0] ||
	      dumpdata->dumplist[f].curr_ymd[1] != ymd[1])
	    chunk_next = 1;
	} else if (dumpdata->dumplist[f].chunk == DAILY) {
	  /* Check Year, Month and Day */
	  if (dumpdata->dumplist[f].curr_ymd[0] != ymd[0] ||
	      dumpdata->dumplist[f].curr_ymd[1] != ymd[1] ||
	      dumpdata->dumplist[f].curr_ymd[2] != ymd[2])
	    chunk_next = 1;
	}
	/*
	 * Close and re-create this dumpfile, if next chunk
	 */
	if (chunk_next) {
	  dumpdata->dumplist[f].close(dumpdata, &dumpdata->dumplist[f]);
	  /*
	   * Need to do this here as dumpdata->t is not quite valid in
	   * the create functions for the transport mode
	   */
	  set_chunk_name(dumpdata, &dumpdata->dumplist[f], t);

	  dumpdata->dumplist[f].private_data = 
	    dumpdata->dumplist[f].create(dumpdata, &dumpdata->dumplist[f]);
	}
      }
      if (dumpdata->us_type & US_IJ)
	dumpdata->dumplist[f].landfill(dumpdata);
      dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
      dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], t);
      dumpdata->dumplist[f].tout += dumpdata->dumplist[f].tinc;
      if (dumpdata->dumplist[f].diagn)
        init_diagn = 1;

    }

    /* Close file if time is past stop time */
    if (t > (dumpdata->dumplist[f].tstop - DT_EPS)
        && !(dumpdata->dumplist[f].finished))
      dumpdata->dumplist[f].close(dumpdata, &dumpdata->dumplist[f]);
  }
  /* Re-initialise diagnostic tracers if necessary */
  if (init_diagn) {
    for (n = 0; n < dumpdata->ntr; ++n) {
      tracer_info_t *tr = &dumpdata->trinfo_3d[n];
      if (tr->diagn) {
        for (k = 0; k < dumpdata->nz; ++k)
          for (j = 0; j < dumpdata->nce2; ++j)
            for (i = 0; i < dumpdata->nce1; ++i)
              dumpdata->tr_wc[n][k][j][i] = 0.0;
      }
    }

    for (n = 0; n < dumpdata->nsed; ++n) {
      tracer_info_t *tr = &dumpdata->trinfo_sed[n];
      if (tr->diagn) {
        for (k = 0; k < dumpdata->sednz; ++k)
          for (j = 0; j < dumpdata->nce2; ++j)
            for (i = 0; i < dumpdata->nce1; ++i)
              dumpdata->tr_sed[n][k][j][i] = 0.0;
      }
    }
  }
  return (init_diagn);
}


/**** Run through the list of dump files closing them down.
 ****/
void dumpfile_close(dump_data_t *dumpdata, double t, FILE * fp)
{
  int f;

  /* Loop through dump file dumplist, writing any file which needs it */
  for (f = 0; f < dumpdata->ndf; f++) {
    /* Write to file if time is between start and stop time, and if file
       is open */
/* TODO UR - check if this spoils thread dispatch 
--was
*/
if (!(dumpdata->dumplist[f].finished)) {

/*    if ((t >= (dumpdata->dumplist[f].tout - DT_EPS)) &&
        (t <= dumpdata->dumplist[f].tstop) &&
        !(dumpdata->dumplist[f].finished)) {
*/
      if (dumpdata->us_type & US_IJ)
	dumpdata->dumplist[f].landfill(dumpdata);
      dump_bathy_mask(dumpdata, dumpdata->dumplist[f].bathymask);
      dumpdata->dumplist[f].write(dumpdata, &dumpdata->dumplist[f], t);
      dumpdata->dumplist[f].close(dumpdata, &dumpdata->dumplist[f]);
    }
  }
}


/*------------------------------------------------------------------*/
/* Counts the number of dumpfiles                                   */
/*------------------------------------------------------------------*/
int count_dumpfiles(FILE *fp, int *type)
{
  FILE *ip, *op;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  int nf, n, m, i;
  double d1;

  if (prm_read_char(fp, "OutputFiles", buf)) {
    /* Is it a file? */
    if ((ip = fopen(buf, "r")) != NULL) {
      if (prm_read_double(ip, "multi-dumpfile-version", &d1)) {
	nf = 0;
	if (prm_read_int(ip, "nfiles", &n)) {
	  for (m = 0; m < n; m++) {
	    sprintf(key, "file%d", m);
	    if (prm_read_char(ip, key, buf)) {
	      if ((op = fopen(buf, "r")) != NULL) {
		if (prm_read_int(op, "OutputFiles", &i)) {
		  nf += i;
		}
		fclose (op);
	      }
	    }
	  }
	} else {
	  hd_warn("No multi-dumpfiles specified.\n");
	  fclose(ip);
	  *type = 0;
	  return(0);
	}
      } else {
	hd_warn("Incorrect format for multi-dumpfile.\n");
	fclose(ip);
	*type = 0;
	return(0);
      }
      *type = 1;
      fclose(ip);
      return(nf);
    } else {
      *type = 0;
      return(atoi(buf));
    }
  } else
    hd_warn("No output files specified\n");
  return(0);
}

/* END count_dumpfiles()                                            */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Reads the dumpfile attributes from file                          */
/*------------------------------------------------------------------*/
void read_dumpfiles(FILE *fp, dump_file_t *list, dump_data_t *dumpdata, 
		    double t, int f)
{
  char buf[MAXSTRLEN], key[MAXSTRLEN], buf2[MAXLINELEN];
  int i, mapIndex = 0;

  list->flag = 0;
  list->append = forced_restart;

  /* Get file name */
  sprintf(key, "file%d.name", f);
  prm_skip_to_end_of_key(fp, key);
  if (fscanf(fp, "%s", list->name) != 1) {
    hd_quit("dumpfile_init: Can't read output file %d name\n", f);
  }

  if (strlen(dumpdata->opath)) {
    sprintf(buf, "%s%s", dumpdata->opath, list->name);
    strcpy(list->name, buf);
  }

  if(!endswith(list->name,".nc"))
    strcat(list->name,".nc");

  /* Cache away the original name, omit .nc */
  strncpy(list->origname, list->name, strlen(list->name)-3);

  /* Set dump limits for near real-time operation */
  if (nrt_restart) {
    list->tout = t;
    list->tstop = schedule->stop_time;
  } else {
    /* Output start time */
    sprintf(key, "file%d.tstart", f);
    prm_skip_to_end_of_key(fp, key);
    if (fgets(buf, MAXSTRLEN, fp) == NULL ||
	!tm_scale_to_secs(buf, &list->tout))
      list->tout = dumpdata->start_time;
    /*hd_quit("dumpfile_init: Can't read file %d start time\n", f);*/

    /* Output stop time */
    sprintf(key, "file%d.tstop", f);
    prm_skip_to_end_of_key(fp, key);
    if (fgets(buf, MAXSTRLEN, fp) == NULL ||
	!tm_scale_to_secs(buf, &list->tstop))
      list->tstop = dumpdata->stop_time;
    /*hd_quit("dumpfile_init: Can't read file %d stop time\n", f);*/
  }
  list->tstart = list->tout;

  /* Output time increment */
  sprintf(key, "file%d.tinc", f);
  /*
    prm_skip_to_end_of_key(fp, key);
    if (fgets(buf, MAXSTRLEN, fp) == NULL ||
    !tm_scale_to_secs(buf, &list->tinc)) {
  */
  list->incf = NONE;
  if (prm_read_char(fp, key, buf)) {
    if (!tm_scale_to_secs(buf, &list->tinc)) {
      if (contains_token(buf, "YEARLY")) {
	list->incf = YEARLY;
	list->tinc = 365.0;
      } else if (contains_token(buf, "SEASONAL")) {
	list->incf = SEASONAL;
	list->tinc = 91.0;
      } else if (contains_token(buf, "MONTHLY")) {
	list->incf = MONTHLY;
	list->tinc = 30.0;
      } else if (contains_token(buf, "DAILY")) {
	list->incf = DAILY;
	list->tinc = 1.0;
      } else
	hd_quit("dumpfile_init: Can't read file %d time increment\n", f);
    }
  }

  /* Fit the start and stop times to the actual simulation period */
  if (t > list->tout) {
    list->tout = ceil((t - list->tout)
		     / list->tinc)*list->tinc + list->tout;
  }

  if (schedule->stop_time <= list->tstop) {
    list->tstop = floor((schedule->stop_time - list->tstop)
		       / list->tinc)*list->tinc + list->tstop;
  }

  /* Bytes per value */
  sprintf(key, "file%d.bytespervalue", f);
  prm_skip_to_end_of_key(fp, key);
  if (fscanf(fp, "%d", &list->bpv) != 1)
    hd_quit("dumpfile_init: Can't read file %d bytes per value\n", f);

  /* Time dependent variables */
  sprintf(key, "file%d.vars", f);
  prm_skip_to_end_of_key(fp, key);
  if (fgets(buf2, MAXLINELEN, fp) == NULL)
    hd_quit("dumpfile_init: Can't read file %d variables\n", f);
  else
    list->nvars = parseline(strdup(buf2), list->vars, MAXNUMVARS);

  /* Set the flags for barotropic variables. If eta, T and S exist  */
  /* in the same file, then sinc = tinc, and synchronisation will   */
  /* occur at t = tinc + dt2d. If eta and T & S are in seperate     */
  /* files, then it's possible to set sinc = tinc - SEPS for T & S, */
  /* then synchronisation will occur at t = tinc.                   */
  /* Set flag & DF_BARO for barotropic variables if the latter is   */
  /* the case.                                                      */
  for (i = 0; i < list->nvars; ++i) {
    if (strcmp(list->vars[i], "eta") == 0) {
      int ii, bf = 0;
      list->flag |= DF_ETA;
      for (ii = 0; ii < list->nvars; ++ii) {
	if (strcmp(list->vars[ii], "temp") == 0 ||
	    strcmp(list->vars[ii], "salt") == 0) bf = 1;
      }
      if (!bf) {
	barof = 1;
	list->flag |= DF_BARO;
      }
    }
    if (strcmp(list->vars[i], "uav") == 0 ||
	strcmp(list->vars[i], "vav") == 0) {
      int ii, bf = 0;
      list->flag |= DF_V2D;
      for (ii = 0; ii < list->nvars; ++ii) {
	if (strcmp(list->vars[ii], "u") == 0 ||
	    strcmp(list->vars[ii], "v") == 0) bf = 1;
      }
      if (!bf && !(list->flag & DF_BARO) && !barof) {
	barof = 1;
	list->flag |= DF_BARO;
      }
    }
  }

  prm_set_errfn(hd_silent_warn);
  /* Read the time modulo */
  memset(list->modulo, 0, strlen(list->modulo));
  sprintf(key, "file%d.modulo", f);
  prm_read_char(fp, key, list->modulo);

  /* Read the output file type. */
  sprintf(key, "file%d.filetype", f);
  if (prm_read_char(fp, key, buf) == 0)
    strcpy(buf, df_func_map[0].tag);
  if (!(dumpdata->us_type & US_IJ) && (strcmp(buf, "ugrid") != 0 && 
				       strcmp(buf, "parray") != 0 &&
				       strcmp(buf, "memory") != 0))
    hd_quit("read_dumpfiles: Unstructured grid can only dump filetype 'ugrid' or 'parray'.\n");

  mapIndex = 0;
  while (df_func_map[mapIndex].tag != NULL) {
    if (strcasecmp(buf, df_func_map[mapIndex].tag) == 0) {
      list->create = df_func_map[mapIndex].create;
      list->write = df_func_map[mapIndex].write;
      list->reset = df_func_map[mapIndex].reset;
      list->close = df_func_map[mapIndex].close;
      strcpy(list->type,buf);
      break;
    }
    ++mapIndex;
  }

  if ( (df_func_map[mapIndex].tag == NULL) && !(params->runmode & PRE_MARVL) )
    hd_quit("dumpfile_init: '%s' is unknown output type.\n", buf);
  if (strcmp(list->type, "memory") == 0) list->flag |= DF_MPK;

  /* Read the time 2-way nest sync time increment */
  sprintf(key, "file%d.sync_dt", f);
  if (prm_read_char(fp, key, buf)) {
    if (!tm_scale_to_secs(buf, &list->sinc))
      list->sinc = list->tinc;
  } else {
    if (list->flag & DF_MPK && !(list->flag & DF_BARO)) {
      /*list->sinc = list->tinc-dumpdata->master->grid_dt/dumpdata->master->iratio;*/
      list->sinc = list->tinc;
    } else
      list->sinc = list->tinc;
  }
  if (list->flag & DF_MPK && barof) {
    if (!(list->flag & (DF_ETA|DF_V2D)) && !(list->flag & DF_BARO))
      list->sinc -= SEPS;
  }

  sprintf(key, "file%d.filter", f);
  read_filter(list, fp, key);

  list->long0_360 = 0;
  sprintf(key, "file%d.long0_360", f);
  if (prm_read_char(fp, key, buf)) {
    if (strcmp(buf, "0:360") == 0)
      list->long0_360 = O360;
    else if (strcmp(buf, "-180:180") == 0)
      list->long0_360 = O180;
    else if (strcmp(buf, "yes") == 0 || strcmp(buf, "no") == 0) {
      if (is_true(buf))
	list->long0_360 = O360;
    }
  }

  list->compress = 0;
  sprintf(key, "file%d.compress", f);
  if (prm_read_char(fp, key, buf))
    list->compress = is_true(buf);
  list->da_cycle = NONE;
  if (params->da) {
    /* Make forecast mode the default */
    list->da_cycle = NO_DA;
      
    sprintf(key, "file%d.da", f);
    if (prm_read_char(fp, key, buf)) {
      /*
       * DA_ flags
       *
       * Please make sure any changes are in sync with ts_init
       */
      if (contains_token(buf, "forecast"))
	list->da_cycle = NO_DA;
      else if (contains_token(buf, "reanalysis"))
	list->da_cycle = DO_DA;
      else
	hd_quit("dumpfile_init: unknown da_cycle '%s' specified\n", buf);
    }
  }
  sprintf(key, "file%d.fill_rule", f);
  if (prm_read_char(fp, key, buf))
    list->landfill = locate_landfill_function(buf);
  else
    list->landfill = locate_landfill_function("default");
  sprintf(key, "file%d.i_rule", f);
  if (!(prm_read_char(fp, key, list->irule)))
    strcpy(list->irule, "linear");
  if(strcmp(list->irule, "linear") == 0)
    list->osl = L_LINEAR;
  else if(strcmp(list->irule, "bilinear") == 0)
    list->osl = L_BILIN;
  else if(strcmp(list->irule, "baylinear") == 0)
    list->osl = L_BAYLIN;
  else if(strcmp(list->irule, "nn_sibson") == 0)
    list->osl = L_SIB;
  else if(strcmp(list->irule, "nn_non_sibson") == 0)
    list->osl = L_NONSIB;
  else if(strcmp(list->irule, "cubic") == 0)
    list->osl = L_CUBIC;
  else if(strcmp(list->irule, "quadratic") == 0)
    list->osl = L_LSQUAD;
  else if(strcmp(list->irule, "linearlsq") == 0)
    list->osl = L_LSLIN;
  else if(strcmp(list->irule, "nearest") == 0)
    list->osl = L_NRST;

  sprintf(key, "file%d.bathy_mask", f);
  if (!(prm_read_double(fp, key, &list->bathymask)))
    list->bathymask = 9999.0;

  list->chunk = NONE; /* default to none */
  sprintf(key, "file%d.chunk", f);
  if (prm_read_char(fp, key, buf)) {
    if (contains_token(buf, "DAILY"))
      list->chunk = DAILY;
    else if (contains_token(buf, "MONTHLY"))
      list->chunk = MONTHLY;
    else if (contains_token(buf, "YEARLY"))
      list->chunk = YEARLY;
    else
      hd_quit("dumpfile_init: unknown outfile chunk '%s' specified\n", buf);
  }
  
  /* Read the subsection */
  list->ilower = 0;
  list->jlower = 0;
  list->klower = 0;
  list->nce1 = dumpdata->nce1;
  list->nfe1 = dumpdata->nfe1;
  list->nce2 = dumpdata->nce2;
  list->nfe2 = dumpdata->nfe2;
  list->nz = dumpdata->nz;
  list->nz_sed = dumpdata->sednz;
  list->ns2 = dumpdata->ns2;
  list->ns3 = dumpdata->ns3;
  list->nface2 = dumpdata->nface2;
  list->nedge2 = dumpdata->nedge2;
  list->nvertex2 = dumpdata->nvertex2;
  list->nface3 = dumpdata->nface3;
  list->nedge3 = dumpdata->nedge3;
  list->nvertex3 = dumpdata->nvertex3;
  if (df_func_map[mapIndex].xyGridOutput) {
    sprintf(key, "file%d.i_range", f);
    read_subsection_range(fp, key, dumpdata->nce1, dumpdata->nfe1,
			  &list->ilower, &list->nce1, &list->nfe1);
    
    sprintf(key, "file%d.j_range", f);
    read_subsection_range(fp, key, dumpdata->nce2, dumpdata->nfe2,
			  &list->jlower, &list->nce2, &list->nfe2);
  }
  sprintf(key, "file%d.k_range", f);
  if (read_subsection_range(fp, key, dumpdata->nz, dumpdata->nz + 1,
			    &list->klower, &list->nz, NULL)) {
    strcpy(buf, list->name);
    if (contains_string(buf, "surf")) {
      list->klower = list->nz - 1;
      list->nz = 1;
    }
  }

  /* Read the list of points to output. Only used by parray */
  if (!df_func_map[mapIndex].xyGridOutput) {
    sprintf(key, "file%d.points", f);
    if(prm_read_char(fp,key,buf))
      read_dumpfile_points(fp,key,list);
    else {
      sprintf(key, "file%d.pointfile", f);
      if(prm_read_char(fp,key,buf)) {
	FILE* fpp = fopen(buf,"r");
	if(fpp == NULL)
	  hd_quit("Points file doesnot exist for: %s %s",key,buf);
	
	sprintf(key, "file%d.points", f);
	read_dumpfile_points(fpp,key,list);
	fclose(fpp);
      }
    }
  }

  /* Read output time units, if specified */
  sprintf(key, "file%d.timeunit", f);
  if (prm_read_char(fp, key, buf))
    strcpy(list->tunit, buf);
  else
    strcpy(list->tunit, dumpdata->output_tunit);

  /* Overwrite file name for chunking */
  set_chunk_name(dumpdata, list, t);

  if (!(params->runmode & PRE_MARVL))
    list->private_data = list->create(dumpdata, list);
  list->finished = 0;
}

/* END read_dumpfiles()                                             */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Sets up and populates a data filtering structure                 */
/*------------------------------------------------------------------*/
void read_filter(dump_file_t *df, FILE *fp, char *key)
{
  int n;
  char buf[MAXSTRLEN];
  df_filter_t *f;
  df->filter = NULL;
  int *st = NULL, sz;

  if (prm_read_char(fp, key, buf)) {
    if (strcmp(buf, "none") == 0 || strcmp(buf, "copy") == 0) return;
    df->filter = (df_filter_t *)malloc(sizeof(df_filter_t));
    f = df->filter;
    if (strcmp(buf, "average3") == 0) {
      strcpy(f->name, "average 3 point");
      f->type = M3PT;
      f->size = 9;
      f->f = 0.0;
    }
    else if (strcmp(buf, "average5") == 0) {
      strcpy(f->name, "average 5 point");
      f->type = M5PT;
      f->size = 25;
      f->f = 0.0;
    }
    else if (strcmp(buf, "weighted3") == 0 || strcmp(buf, "shapiro3") == 0) {
      strcpy(f->name, "shapiro 3 point");
      f->type = W3PT;
      f->size = 9;
      f->f = 0.0;
    }
    else if (strcmp(buf, "weighted5") == 0) {
      strcpy(f->name, "weighted 5 point");
      f->type = W5PT;
      f->size = 25;
      f->f = 0.0;
      f->v = 1.0 / 9.0;
    }
    else if (strcmp(buf, "highpass3") == 0) {
      strcpy(f->name, "hi-pass 3 point");
      f->type = HP3PT;
      f->size = 9;
      f->f = 0.25;
    }
    else if (strcmp(buf, "shuman3") == 0) {
      strcpy(f->name, "shuman 3 point");
      f->type = HP3PT;
      f->size = 9;
      f->f = 1.0;
      f->v = 0.25;
    } else {
      hd_warn("read_filter: Can't recognise filter %s for dumpfile %s\n", key, df->name);
      free ((df_filter_t *)df->filter);
      df->filter = NULL;
      return;
    }

    sz = sqrt(f->size) / 2;
    sz = 1 + sz * geom->npem;

    f->k = d_alloc_1d(sz);
    for (n = 0; n < sz; n++)
      f->k[n] = 1.0 / (double)f->size;
    f->map = i_alloc_1d(sz);
    
    if (f->type & W3PT) {
      f->k[0] = f->k[2] = f->k[6] = f->k[8] = 0.0625;
      f->k[1] = f->k[3] = f->k[5] = f->k[7] = 0.125;
      f->k[4] = 0.25;
    }
    if (f->type & HP3PT) {
      f->k[0] = f->k[2] = f->k[6] = f->k[8] = 0.0;
      f->k[1] = f->k[3] = f->k[5] = f->k[7] = 1.0;
      f->k[4] = -4.0;
    }
    if (f->type & SHU3PT) {
      f->k[0] = f->k[2] = f->k[6] = f->k[8] = 0.25 * f->v * f->v;
      f->k[1] = f->k[3] = f->k[5] = f->k[7] = 0.5 * f->v * (1.0 - f->v);
      f->k[4] = -2.0 * f->v * (1.0 - f->v) - f->v * f->v;
    }
    if (f->type & W5PT) {
      f->k[0] = f->k[4] = f->k[20] = f->k[24] = 1.0 * f->v * f->v;
      f->k[1] = f->k[3] = f->k[21] = f->k[23] = 2.0 * f->v * f->v;
      f->k[5] = f->k[9] = f->k[15] = f->k[19] = 2.0 * f->v * f->v;
      f->k[2] = f->k[10] = f->k[14] = f->k[22] = 3.0 * f->v * f->v;
      f->k[6] = f->k[8] = f->k[16] = f->k[18] = 4.0 * f->v * f->v;
      f->k[7] = f->k[11] = f->k[13] = f->k[17] = 6.0 * f->v * f->v;
      f->k[12] = 9.0 * f->v * f->v;
    }
  }
}

/* END read_filter()                                                */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialise the dumpfile structure and populate                   */
/*------------------------------------------------------------------*/
void dumpfile_init(dump_data_t *dumpdata, double t, FILE * fp, int *n,
		   dump_file_t **dumplist)
{
  char key[50];
  char buf[MAXSTRLEN];
  char buf2[MAXLINELEN];
  char path[MAXLINELEN];
  dump_file_t *list = NULL;
  int f;
 /* int i = 0; */
  int nfiles = 0;
  int transfiles = 0;
  double restart_dt = 0.0;
  int type;

  prm_set_errfn(hd_silent_warn);

  if (prm_read_char(fp, "restart_dt", buf)) {
     tm_scale_to_secs(buf, &restart_dt);
  }

  /* Get number of global output files */
  nfiles = count_dumpfiles(fp, &type);

  prm_get_time_in_secs(fp, "START_TIME", &dumpdata->start_time);
  prm_get_time_in_secs(fp, "STOP_TIME", &dumpdata->stop_time);

  if (strlen(dumpdata->trkey)) {
    transfiles = get_transfiles(dumpdata, fp);
    nfiles += transfiles;
  }

  prm_set_errfn(hd_quit);
  *n = nfiles + (restart_dt > 0.0);
  if (*n <= 0) {
    hd_warn("dumpfile_init: No netCDF output files have been specified.\n");
    return;
  }
  if (*n > MAXDUMPFILES) {
    hd_warn("dumpfile_init: Maximum of %d output files, rest ignored\n",
            MAXDUMPFILES);
    *n = MAXDUMPFILES;
  }

/* Allocate memory for dump_file_t array */
  list = (dump_file_t *)malloc(sizeof(dump_file_t) * (*n));
  memset(list, 0, sizeof(dump_file_t) * (*n));

  prm_set_errfn(hd_quit);

  /* Set up each file list entry */
  strcpy(dumpdata->rev, params->rev);
  sprintf(dumpdata->grid_name, "%s", params->grid_name);
  if (params->runno > 0) dumpdata->runno = params->runno;
  if (type == 1) {
    FILE *ip, *op;
    int nf, nn, m, i;
    double d1;
    prm_read_char(fp, "OutputFiles", buf);
    if ((ip = fopen(buf, "r")) != NULL) {
      if (prm_read_double(ip, "multi-dumpfile-version", &d1)) {
	f = 0;
	if (prm_read_int(ip, "nfiles", &nf)) {
	  for (m = 0; m < nf; m++) {
	    sprintf(key, "file%d", m);
	    if (prm_read_char(ip, key, buf)) {
	      if ((op = fopen(buf, "r")) != NULL) {
		if (prm_read_int(op, "OutputFiles", &i)) {
		  for (nn = 0; nn < i; nn++) {
		    read_dumpfiles(op, &list[f], dumpdata, t, nn);
		    f++;
		  }
		}
		fclose(op);
	      }
	    }
	  }
	}
      }
      fclose(ip);
    }
  } else if (type == 0) {
    for (f = 0; f < nfiles - transfiles; f++) {
      read_dumpfiles(fp, &list[f], dumpdata, t, f);
    }
  }
  if (transfiles) {
    f = set_transfiles(nfiles - transfiles, dumpdata->trkey, dumpdata, list, t, fp);
    if (f != nfiles)
      hd_quit("Cannot infer correct number of output transport files.\n");
  }

  /* Add the restart file to the end */
  if (restart_dt > 0.0) {
    if (prm_read_char(fp, "restart_name", buf))
      strcpy(list[nfiles].name, buf);
    else
      strcpy(list[nfiles].name, "restart.nc");
    list[nfiles].tout = t;
    list[nfiles].tstop = schedule->stop_time;
    list[nfiles].tinc = restart_dt;
    list[nfiles].bpv = 8;
    list[nfiles].nvars = parseline(strdup("ALL"), list[nfiles].vars,
				   MAXNUMVARS);

    list[nfiles].create = df_restart_create;
    list[nfiles].write = df_restart_write;
    list[nfiles].close = df_restart_close;
    list[nfiles].landfill = locate_landfill_function("default");
    
    list[nfiles].ilower = 0;
    list[nfiles].jlower = 0;
    list[nfiles].klower = 0;
    list[nfiles].nce1 = dumpdata->nce1;
    list[nfiles].nfe1 = dumpdata->nfe1;
    list[nfiles].nce2 = dumpdata->nce2;
    list[nfiles].nfe2 = dumpdata->nfe2;
    list[nfiles].nz = dumpdata->nz;
    list[nfiles].nz_sed = dumpdata->sednz;
    list[nfiles].ns2 = dumpdata->ns2;
    list[nfiles].ns3 = dumpdata->ns3;
    list[nfiles].nface2 = dumpdata->nface2;
    list[nfiles].nedge2 = dumpdata->nedge2;
    list[nfiles].nvertex2 = dumpdata->nvertex2;
    list[nfiles].nface3 = dumpdata->nface3;
    list[nfiles].nedge3 = dumpdata->nedge3;
    list[nfiles].nvertex3 = dumpdata->nvertex3;
    list[nfiles].da_cycle = (NO_DA|DO_DA|NONE); // write in all cases
    list[nfiles].private_data = list[nfiles].create(dumpdata, &list[nfiles]);
    list[nfiles].finished = 0;
    list[nfiles].compress = 0;
    list[nfiles].bathymask = 9999.0;
    if (prm_read_char(fp, "restart_compress", buf))
      list[nfiles].compress = is_true(buf);
    strcpy(list[nfiles].tunit, dumpdata->output_tunit);
  }

  /* If re-initialisation of diagnostic tracers has not been synced with
     one of the output files, sync it with the very first one. */
  if (!dumpdata->df_diagn_set && *n > 0) {
    list[0].diagn = 1;
    dumpdata->df_diagn_set = 1;
  }
#if defined(HAVE_ECOLOGY_MODULE)
  /* Check if diagnostic tracers re-initialisation is synced with ecology.
   */
  if (dumpdata->do_eco) {
    double eco_out;
    for (f = 0; f < *n; ++f)
      if (list[f].diagn)
        break;
    eco_out =
      (dumpdata->t +
       (int)((list[f].tout + DT_EPS -
              dumpdata->t) / dumpdata->ecodt) * dumpdata->ecodt);
    if ((list[f].tinc / dumpdata->ecodt) !=
        (int)(list[f].tinc / dumpdata->ecodt))
      hd_warn
        ("Diagnostic tracers initialisation can not be synchronized with ecology step: ecology step = %.1fs, output step = %.1fs.\n",
         dumpdata->ecodt, list[f].tinc);
    else if ((list[f].tout - eco_out) > DT_EPS) {
      hd_warn
        ("Changed output start for \"%s\" from %.1fs to %.1fs to syncronise initialisation of diagnostic tracers with ecology stepping.\n",
         list[f].name, list[f].tout, eco_out);
      list[f].tout = eco_out;
    }
  }
#endif

  *dumplist = list;

}

/* END dumpfile_init()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Close dumpfiles, re-initialise and populate the dumpfile         */
/* structure.                                                       */
/*------------------------------------------------------------------*/
void dumpfile_resetup(master_t *master)
{
  geometry_t *geom = master->geom;
  dump_data_t *dumpdata = master->dumpdata;
  int n, f, ndf, fr = forced_restart;
  char **dname;
  int *next;

  ndf = dumpdata->ndf;
  dname = (char **)malloc(ndf * sizeof(char *));
  next = i_alloc_1d(ndf);

  /*----------------------------------------------------------------*/
  /* Close the dumpfiles and free the dumplist                      */
  for (f = 0; f < ndf; f++) {
    dname[f] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    strcpy(dname[f], dumpdata->dumplist[f].name);
    next[f] = dumpdata->dumplist[f].tout;
    if (!(dumpdata->dumplist[f].finished)) {
      dumpdata->dumplist[f].close(dumpdata, &dumpdata->dumplist[f]);
    }
  }
  free((dump_file_t *) dumpdata->dumplist);

  /*----------------------------------------------------------------*/
  /* Reread and initialise the new dumplist                         */
  forced_restart = 1;
  dumpfile_init(dumpdata, dumpdata->t-DT_EPS, master->prmfd, &dumpdata->ndf,
                &dumpdata->dumplist);

  /*----------------------------------------------------------------*/
  /* Reset the next dump events                                     */
  for (f = 0; f < dumpdata->ndf; f++) {
    dump_file_t *df = &dumpdata->dumplist[f];

    /* Skip dump if the DA flag doesn't match for this dumpfile     */
    if (dumpdata->master->da && !(df->da_cycle & dumpdata->master->da)) 
      continue;
	
    for (n = 0; n < ndf; n++) {
      if (strcmp(dname[n], df->name) == 0) {
	df->tout = next[n];
      }
    }
    /* Set the record number for the next dump                      */
    df->reset(dumpdata, df, dumpdata->t+df->tinc);
  }

  forced_restart = fr;
  free((char **)dname);
  i_free_1d(next);
}

/* END dumpfile_resetup()                                           */
/*------------------------------------------------------------------*/


static void read_dumpfile_points(FILE * fp, char *key, dump_file_t* file)
{
  int i;
  char dummy[MAXSTRLEN];
  long fpos;
  prm_read_int(fp, key, &file->npoints);
  prm_flush_line(fp);
  if (file->npoints > 0) {
    file->x = d_alloc_1d(file->npoints);
    file->y = d_alloc_1d(file->npoints);
    for (i = 0; i < file->npoints; ++i)
      {
	fpos = ftell(fp);
	/* Read the x,y points first */
	if (fscanf(fp, "%lf %lf", &file->x[i], &file->y[i]) != 2) {
	  /* Reset file position */
	  fseek(fp, fpos, SEEK_SET);
	  /* Now try with labels */
	  if (fscanf(fp, "%s %lf %lf", dummy, &file->x[i], &file->y[i]) != 3)
	    hd_quit
	      ("dumpfile_init: Output point array incorrect formed for %s ",key);
	  else {
	    if (file->label == NULL)
	      file->label = c_alloc_2d(MAXSTRLEN, file->npoints);
	    sprintf(file->label[i], "%s", dummy);
	  }
	}
      }
  }
}



/**
 * parse the assigned variables output and rebuild the 
 * 
 */
void df_parse_vars(dump_data_t *dumpdata, dump_file_t *df, char* excludevars, char* all_vars)
{
	
	int i, n;
  char subset[MAXNUMVARS][MAXSTRLEN];
  char buffer[MAXSTRLEN]; 
  char line[MAXNUMVARS * 20] = "";

  for (i = 0; i < df->nvars; ++i) {
    if (contains_token(df->vars[i], "ALL")) {
      /* sprintf(line, STD_ALL_VARS); */
      sprintf(line, "%s", all_vars);
      for (n = 0; n < dumpdata->ntr; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);

      for (n = 0; n < dumpdata->nsed; n++) {
			  char name[MAXSTRLEN];
			  strcpy(name, dumpdata->trinfo_sed[n].name);
			  strcat(name, "_sed");
			  sprintf(line + strlen(line), "%s ", name);
      }

      for (n = 0; n < dumpdata->ntrS; n++)
        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_2d[n].name);
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    } else if (contains_token(df->vars[i], "NONTRACERS"))
    {
      /* sprintf(line, STD_ALL_VARS); */
      sprintf(line, "%s", all_vars);
    }
    else if (contains_token(df->vars[i], "TRACERS_WC")) {
      for (n = 0; n < dumpdata->ntr; n++)
        if (!dumpdata->trinfo_3d[n].diagn)
          sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);
    } else if (contains_token(df->vars[i], "TRACERS_DIAGN_WC")) {
      for (n = 0; n < dumpdata->ntr; n++)
        if (dumpdata->trinfo_3d[n].diagn)
          sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    } else if (contains_token(df->vars[i], "TRACERS_SED")) {
      for (n = 0; n < dumpdata->nsed; n++)
        if (!dumpdata->trinfo_sed[n].diagn) {
          char name[MAXSTRLEN];
          strcpy(name, dumpdata->trinfo_sed[n].name);
          strcat(name, "_sed");
          sprintf(line + strlen(line), "%s ", name);
        }
    } else if (contains_token(df->vars[i], "TRACERS_DIAGN_SED")) {
      for (n = 0; n < dumpdata->nsed; n++) 
	    {
  		  if (dumpdata->trinfo_sed[n].diagn) 
        {
	        char name[MAXSTRLEN];
	        strcpy(name, dumpdata->trinfo_sed[n].name);
	        strcat(name, "_sed");
	        sprintf(line + strlen(line), "%s ", name);
	      }
      }
      if (!dumpdata->df_diagn_set) {
        df->diagn = 1;
        dumpdata->df_diagn_set = 1;
      }
    }else if (contains_char( WILDCARD,df->vars[i]))
    {
    	/* Scan all names for the pattern given as name 
    	 * and add them
    	 */
    	int n,nel;
    	char name[MAXSTRLEN];
    	strcpy(buffer,STD_ALL_VARS);
    	nel = parseline(buffer,(char **)subset,MAXNUMVARS);
    	

    	for (n = 0; n < nel; n++)
    	{
    		if(contains_pattern(((char **)subset)[n],df->vars[i]))
    			sprintf(line + strlen(line), "%s ",((char **)subset)[n]);
    	}
    	
    	for (n = 0; n < dumpdata->ntr; n++)
    	{
    		if(contains_pattern(dumpdata->trinfo_3d[n].name,df->vars[i]))
        	sprintf(line + strlen(line), "%s ", dumpdata->trinfo_3d[n].name);
    	}
    	
    	for (n = 0; n < dumpdata->nsed; n++) 
    	{
    		if(contains_pattern(dumpdata->trinfo_sed[n].name,df->vars[i]))
    		{
    			strcpy(name, dumpdata->trinfo_sed[n].name);
			  	strcat(name, "_sed");
        	sprintf(line + strlen(line), "%s ", name);
    		}
    	}
    	
    	for (n = 0; n < dumpdata->ntrS; n++)
    	{
    		if(contains_pattern(dumpdata->trinfo_2d[n].name,df->vars[i]))
	        sprintf(line + strlen(line), "%s ", dumpdata->trinfo_2d[n].name);
    	} 
    
    }else
      sprintf(line + strlen(line), "%s ", df->vars[i]);
  }

  /* Start Cleaning duplicates */
  df->nvars = parseline(strdup(line), df->vars, MAXNUMVARS);
	
	strcpy(line, "");
	if(excludevars == NULL)
	{
		for (i = 0; i < df->nvars; ++i)
	    if (!contains_token(line, df->vars[i]))
	      sprintf(line + strlen(line), "%s ", df->vars[i]);
	}else{
	  for (i = 0; i < df->nvars; ++i)
	    if (!contains_token(line, df->vars[i]) && !contains_token(excludevars, df->vars[i]))
	      sprintf(line + strlen(line), "%s ", df->vars[i]);
	}
	
	/*
	 * (FR) I think this is here in an attempt to free the strdup
	 *      in the call to parseline above. In which case the free
	 *      is wrong as the pointer has effectively been stolen by
	 *      df->vars[0]. Parseline divvy's up 'line' into tokens
	 *      without allocation so, for now, we must live with
	 *      sending in duplicates of the line. This also means
	 *      that programs such as valgrind get very confused
	 */
	// free(df->vars[0]);

  /* Re-parse the variable list, and check number */
  df->nvars = parseline(strdup(line), df->vars, MAXNUMVARS);
  
}



/* History output file format. in large part this is a front to the standard
 * file format, except that the file is over-written each time.
 */
static void *df_restart_create(dump_data_t *dumpdata, dump_file_t *df) {

  /* Force a FULL dump */
  if ((df->nvars != 1) && (strcmp(df->vars[0], "ALL") != 0)) {
     hd_quit("All variables must be specified for restart file '%s'.\n",
           df->name);
  }

  if ( (df->ilower != 0) || (df->jlower != 0) || (df->klower != 0)
    || (df->nce1 != dumpdata->nce1) || (df->nce2 != dumpdata->nce2)
    || (df->nz != dumpdata->nz)) {
    hd_quit("Full domain must be specified for restart file '%s'.\n", df->name);
  }

  if (df->bpv != 8)
    hd_quit("History file '%s' must use 8 bytes per value'.\n", df->name);

  df->append = 0;       /* Ensure that file is not appended to */

  return NULL;
}

static void df_restart_write(dump_data_t *dumpdata, dump_file_t *df, double t) {
  FILE *fp = NULL;
  char restart_fname[MAXSTRLEN];
  char new_restart_fname[MAXSTRLEN];

  /* Check that the dumptime is one of the valid dump times for
   * a restart file (i.e. not an error dump). */
  if (t < (df->tout - DT_EPS)) return;

  /* Write the restart file. */
  strcpy(restart_fname, df->name);
  strcpy(new_restart_fname, dirname(strdup(df->name)));
  strcat(new_restart_fname, "/new_restart.nc");
  strcpy(df->name, new_restart_fname);
  df_ugrid_create(dumpdata, df);
  df_ugrid_write(dumpdata, df, t);
  df_ugrid_close(dumpdata, df);
  df->finished = 0;
  strcpy(df->name, restart_fname);
  /* Check if old restart file exists. If so remove it,
   * and replace with the new restart file. */
  fp = fopen(restart_fname, "r");
  if (fp != NULL) {
    fclose(fp);
    unlink(restart_fname);
  }
  rename(new_restart_fname, restart_fname);
}

static void df_restart_close(dump_data_t *dumpdata, dump_file_t *df) {
  /* Don't ever finish/close */
  // df->finished = 1;
}

void df_std_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_std_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t <= (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "t", df->tout);
}

/* Standard output file format */
void *df_std_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int n;
  int cdfid;
  int dims[10];
  int vid;                      /* df_variable_t dimension */
  /* dimension ids */
  int recdimid;                 /* record dimension id */
  int igridid;                  /* I dimension id at grid corner */
  int jgridid;                  /* J dimension id at grid corner */
  int kgridid;                  /* K dimension id at grid corner */
  int icentreid;                /* I dimension id at grid centre */
  int jcentreid;                /* J dimension id at grid centre */
  int kcentreid;                /* K dimension id at grid centre */
  int ileftid;                  /* I dimension id at left face */
  int jleftid;                  /* J dimension id at left face */
  int ibackid;                  /* I dimension id at back face */
  int jbackid;                  /* J dimension id at back face */
  int kcentreid_sed;            /* K dimension id at grid centre for
                                   sediments */
  int kgridid_sed;              /* K dimension id at grid corner for
                                   sediments */
  int status;
  int nc_mode;

  df_std_data_t *data = NULL;
  nc_type fptype = nc_get_default_type(df->bpv);
  FILE *fp;

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_std_init_data(dumpdata, df, cdfid);
      data = (df_std_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "t", dumpdata->t);
      df->finished = 0;
      return data;
    }
  }

  // Get the nc_mode flag
  nc_mode = get_nc_mode(df);
  
  /* Create the netCDF file using netCDF4 format */
  status = ncw_create(df->name, nc_mode, &cdfid);
  if (status != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "record", NC_UNLIMITED, &recdimid);

  nc_def_dim(cdfid, "k_grid", df->nz + 1, &kgridid);
  nc_def_dim(cdfid, "j_grid", df->nce2 + 1, &jgridid);
  nc_def_dim(cdfid, "i_grid", df->nce1 + 1, &igridid);

  nc_def_dim(cdfid, "k_centre", df->nz, &kcentreid);
  nc_def_dim(cdfid, "j_centre", df->nce2, &jcentreid);
  nc_def_dim(cdfid, "i_centre", df->nce1, &icentreid);

  nc_def_dim(cdfid, "j_left", df->nce2, &jleftid);
  nc_def_dim(cdfid, "i_left", df->nce1 + 1, &ileftid);

  nc_def_dim(cdfid, "j_back", df->nce2 + 1, &jbackid);
  nc_def_dim(cdfid, "i_back", df->nce1, &ibackid);

  if (df->nz_sed > 0) {
    nc_def_dim(cdfid, "k_grid_sed", df->nz_sed + 1, &kgridid_sed);
    nc_def_dim(cdfid, "k_centre_sed", df->nz_sed, &kcentreid_sed);
  }

  /* time independent variables */
  dims[0] = kgridid;
  nc_def_var(cdfid, "z_grid", NC_DOUBLE, 1, dims, &vid);

  dims[0] = kcentreid;
  nc_def_var(cdfid, "z_centre", NC_DOUBLE, 1, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "x_grid", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_grid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "x_centre", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_centre", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "x_left", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_left", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "x_back", NC_DOUBLE, 2, dims, &vid);
  nc_def_var(cdfid, "y_back", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "botz", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "h1au1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "h1au2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "h1acell", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "h1agrid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "h2au1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "h2au2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "h2acell", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jgridid;
  dims[1] = igridid;
  nc_def_var(cdfid, "h2agrid", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jleftid;
  dims[1] = ileftid;
  nc_def_var(cdfid, "thetau1", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jbackid;
  dims[1] = ibackid;
  nc_def_var(cdfid, "thetau2", NC_DOUBLE, 2, dims, &vid);

  dims[0] = jcentreid;
  dims[1] = icentreid;
  nc_def_var(cdfid, "coriolis", fptype, 2, dims, &vid);

  dims[0] = icentreid;
  nc_def_var(cdfid, "crci", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "clci", NC_SHORT, 1, dims, &vid);

  dims[0] = igridid;
  nc_def_var(cdfid, "crfi", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "clfi", NC_SHORT, 1, dims, &vid);

  dims[0] = icentreid;
  nc_def_var(cdfid, "frci", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "flci", NC_SHORT, 1, dims, &vid);

  dims[0] = igridid;
  nc_def_var(cdfid, "frfi", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "flfi", NC_SHORT, 1, dims, &vid);

  dims[0] = jcentreid;
  nc_def_var(cdfid, "cfcj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "cbcj", NC_SHORT, 1, dims, &vid);

  dims[0] = jgridid;
  nc_def_var(cdfid, "cffj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "cbfj", NC_SHORT, 1, dims, &vid);

  dims[0] = jcentreid;
  nc_def_var(cdfid, "ffcj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "fbcj", NC_SHORT, 1, dims, &vid);

  dims[0] = jgridid;
  nc_def_var(cdfid, "fffj", NC_SHORT, 1, dims, &vid);
  nc_def_var(cdfid, "fbfj", NC_SHORT, 1, dims, &vid);

  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid, "t", NC_DOUBLE, 1, dims, &vid);

  if (dumpdata->sednz > 0) {
    dims[0] = recdimid;
    dims[1] = kgridid_sed;
    dims[2] = jcentreid;
    dims[3] = icentreid;
    nc_def_var(cdfid, "z_grid_sed", NC_DOUBLE, 4, dims, &vid);

    dims[1] = kcentreid_sed;
    nc_def_var(cdfid, "z_centre_sed", NC_DOUBLE, 4, dims, &vid);
  }

  df_std_init_data(dumpdata, df, cdfid);
  data = (df_std_data_t *)df->private_data;
  for (n = 0; n < df->nvars; n++) {
    if (data->vars[n].ndims == 2) {
      df_std_get_dimids(cdfid, data->vars[n], &dims[2], &dims[1], NULL);
      ncw_def_var2(df->name, cdfid, df->vars[n], data->vars[n].type, 
			    3, dims, &vid, df->compress);
    } else if (data->vars[n].ndims == 3) {
      df_std_get_dimids(cdfid, data->vars[n],
                        &dims[3], &dims[2], &dims[1]);
      ncw_def_var2(df->name, cdfid, df->vars[n], data->vars[n].type, 
			    4, dims, &vid, df->compress);
    } else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }
  }

  write_dump_attributes(dumpdata, cdfid, fptype, df->ilower, df->nce1,
                        df->nfe1, df->jlower, df->nce2, df->nfe2, df->modulo,
			df->tunit);

  nc_enddef(cdfid);

  df_std_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}


double df_filter_ij(dump_file_t *df, 
		    double **a, 
		    unsigned long **flag,
		    int ne1, int ne2,
		    int i, int j)
{
  int n, s;
  int ii, jj, ip, jp;
  double b;
  double aa, ks;
  df_filter_t *f = df->filter;

  s = sqrt(f->size) / 2;
  if (flag[j][i] & (SOLID|OUTSIDE)) {
    b = a[j][i];
    return(b);
  }
  b = (f->f) ? f->f * a[j][i] : 0.0;
  n = 0;
  ks = 0.0;
  for (jj = j-s; jj <= j+s; jj++) {
    for (ii = i-s; ii <= i+s; ii++) {
      ip = max(ii, 0);
      ip = min(ip, ne1-1);
      jp = max(jj, 0);
      jp = min(jp, ne2-1);
      if (!(flag[jp][ip] & (SOLID|OUTSIDE))) {
	b += f->k[n] * a[jp][ip];
	ks += f->k[n];
	n++;
      }
    }
  }
  if (ks) b *= 1.0 / ks;
  return(b);
}
 
double df_filter(dump_file_t *df, 
		 double *a, 
		 int c)
{
  int n, s, cn;
  int ii, jj, ip, jp;
  double b;
  double aa, ks;
  df_filter_t *f = df->filter;
  int *st = NULL, sz;

  sz = sqrt(f->size);
  st = stencil(geom, c, &sz, ST_SIZED, 0);

  if (geom->wgst[c]) {
    b = a[c];
    return(b);
  }
  b = (f->f) ? f->f * a[c] : 0.0;
  ks = 0.0;

  for (n = 0; n < sz; n++) {
    cn = st[n];
    b += f->k[n] * a[c];
    ks += f->k[n];
  }
  if (ks) b *= 1.0 / ks;
  return(b);
}

void df_filter_1d(dump_file_t *df, 
		  dump_data_t *dumpdata, 
		  double ***a, 
		  int i, int j)
{
  int k;

  for (k = 0; k < dumpdata->nz; k++) {
    dumpdata->w3[k][j][i] = df_filter_ij(df, a[k], dumpdata->flag[k], 
					 dumpdata->nce1, dumpdata->nce2, i, j);
  }
}


void df_filter_2d(dump_file_t *df, 
		  dump_data_t *dumpdata, 
		  double **a, 
		  int is, int ie, int js, int je)
{
  int i, j;
  int ii, jj, ip, jp;
  double **b = dumpdata->w2;
  double aa, ks;
  df_filter_t *f = df->filter;

  for (j = js; j < je; j++)
    for (i = is; i < ie; i++)
      dumpdata->w2[j][i] = df_filter_ij(df, a, dumpdata->flag[dumpdata->nz-1], 
					dumpdata->nce1, dumpdata->nce2, i, j);
}


void df_filter_3d(dump_file_t *df, 
		  dump_data_t *dumpdata, 
		  double ***a, 
		  int is, int ie, int js, int je, int ks, int ke)
{
  int i, j, k;

  for (k = ks; k < ke; k++)
    for (j = js; j < je; j++)
      for (i = is; i < ie; i++)
	dumpdata->w3[k][j][i] = df_filter_ij(df, a[k], dumpdata->flag[k], 
					     dumpdata->nce1, dumpdata->nce2, i, j);
}


void df_std_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double st, newt = t;
  df_std_data_t *data = (df_std_data_t *)df->private_data;
  int fid = data->fid;

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, df->tunit, &newt,
                         1);
  if (strlen(df->modulo)) {
    st = df->tstart;
    tm_change_time_units(dumpdata->timeunit, df->tunit, 
			 &st, 1);
    newt -= st;
  }
  nc_put_vara_double(fid, ncw_var_id(fid, "t"), start, count, &newt);

  if (dumpdata->sednz > 0) {
    count[1] = dumpdata->sednz;
    count[2] = dumpdata->nce2;
    count[3] = dumpdata->nce1;
    nc_d_writesub_3d(fid, ncw_var_id(fid, "z_centre_sed"), start, count,
                     dumpdata->cellz_sed);

    count[1] = dumpdata->sednz + 1;
    nc_d_writesub_3d(fid, ncw_var_id(fid, "z_grid_sed"), start, count,
                     dumpdata->gridz_sed);
  }

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_std_var_t vn = data->vars[n];
    void *p = vn.v;

    if (strcmp("flag", df->vars[n]) == 0) {
      start[1] = df->klower;
      start[2] = df->jlower;
      start[3] = df->ilower;
      df_std_get_dimsizes(df, vn, &count[3], &count[2], &count[1]);
      nc_l_writesub_3d(fid, ncw_var_id(fid, df->vars[n]), start, count,
                       dumpdata->flag);
    } else if (vn.ndims == 2) {
      start[1] = df->jlower;
      start[2] = df->ilower;
      df_std_get_dimsizes(df, vn, &count[2], &count[1], NULL);
      if (df->filter) {
	df_filter_2d(df, dumpdata, *((double ***)p), start[2], count[2], 
		     start[1], count[1]);
	p = &dumpdata->w2;
      }
      nc_d_writesub_2d(fid, ncw_var_id(fid, df->vars[n]), start, count,
                       (*((double ***)p)));
    } else if (vn.ndims == 3) {
      start[1] = (!vn.sediment) ? df->klower : 0;
      start[2] = df->jlower;
      start[3] = df->ilower;
      df_std_get_dimsizes(df, vn, &count[3], &count[2], &count[1]);
      if (df->filter) {
	df_filter_3d(df, dumpdata, *((double ****)p), start[3], count[3], 
		     start[2], count[2], start[1], count[1]);
	p = &dumpdata->w3;
      }
      nc_d_writesub_3d(fid, ncw_var_id(fid, df->vars[n]), start, count,
                       (*((double ****)p)));
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

void df_std_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_std_data_t *data = (df_std_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}


static void df_std_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[4];
  size_t count[4];
  df_std_data_t *data = (df_std_data_t *)df->private_data;
  int fid = data->fid;

  set_longitude(dumpdata, df, 1);
  start[0] = df->klower;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = df->nz + 1;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  nc_d_writesub_1d(fid, ncw_var_id(fid, "z_grid"), start, count,
                   dumpdata->gridz);
  start[0] = df->klower;
  count[0] = df->nz;
  nc_d_writesub_1d(fid, ncw_var_id(fid, "z_centre"), start, count,
                   dumpdata->cellz);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2 + 1;
  count[1] = df->nce1 + 1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_grid"), start, count,
                   dumpdata->gridx);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_grid"), start, count,
                   dumpdata->gridy);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h1agrid"), start, count,
                   dumpdata->h1agrid);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h2agrid"), start, count,
                   dumpdata->h2agrid);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_centre"), start, count,
                   dumpdata->cellx);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_centre"), start, count,
                   dumpdata->celly);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h1acell"), start, count,
                   dumpdata->h1acell);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h2acell"), start, count,
                   dumpdata->h2acell);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "coriolis"), start, count,
                   dumpdata->coriolis);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "botz"), start, count,
                   dumpdata->botz);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1 + 1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_left"), start, count,
                   dumpdata->u1x);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_left"), start, count,
                   dumpdata->u1y);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h1au1"), start, count,
                   dumpdata->h1au1);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h2au1"), start, count,
                   dumpdata->h2au1);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "thetau1"), start, count,
                   dumpdata->thetau1);
  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2 + 1;
  count[1] = df->nce1;
  nc_d_writesub_2d(fid, ncw_var_id(fid, "x_back"), start, count,
                   dumpdata->u2x);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "y_back"), start, count,
                   dumpdata->u2y);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h1au2"), start, count,
                   dumpdata->h1au2);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "h2au2"), start, count,
                   dumpdata->h2au2);
  nc_d_writesub_2d(fid, ncw_var_id(fid, "thetau2"), start, count,
                   dumpdata->thetau2);

  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}

static int df_std_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
                              char *name, df_std_var_t *var)
{
  int found = 1;
  int n = 0;

  var->v = NULL;
  var->ndims = 2;
  var->type = nc_get_default_type(df->bpv);
  var->xylocation = CL_CENTRE;
  var->zlocation = CL_NONE;
  var->sediment = 0;

  if (strcmp(name, "u1av") == 0) {
    var->v = (void **)&dumpdata->u1av;
    var->xylocation = CL_LEFT;
  }

  else if (strcmp(name, "u2av") == 0) {
    var->v = (void **)&dumpdata->u2av;
    var->xylocation = CL_BACK;
  }

  else if (strcmp(name, "wtop") == 0) {
    var->v = (void **)&dumpdata->wtop;
  }

  else if (strcmp(name, "topz") == 0) {
    var->v = (void **)&dumpdata->topz;
  }

  else if (strcmp(name, "eta") == 0) {
    var->v = (void **)&dumpdata->eta;
  }

  else if (strcmp(name, "wind1") == 0) {
    var->v = (void **)&dumpdata->wind1;
    var->xylocation = CL_LEFT;
  }

  else if (strcmp(name, "wind2") == 0) {
    var->v = (void **)&dumpdata->wind2;
    var->xylocation = CL_BACK;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void **)&dumpdata->patm;
  }

  else if (strcmp(name, "u1") == 0) {
    var->v = (void **)&dumpdata->u1;
    var->ndims = 3;
    var->xylocation = CL_LEFT;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "u2") == 0) {
    var->v = (void **)&dumpdata->u2;
    var->ndims = 3;
    var->xylocation = CL_BACK;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void **)&dumpdata->w;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "u") == 0) {
    var->v = (void **)&dumpdata->u;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "v") == 0) {
    var->v = (void **)&dumpdata->v;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "dens") == 0) {
    var->v = (void **)&dumpdata->dens;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "dens_0") == 0) {
    var->v = (void **)&dumpdata->dens_0;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "Kz") == 0) {
    var->v = (void **)&dumpdata->Kz;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "Vz") == 0) {
    var->v = (void **)&dumpdata->Vz;
    var->ndims = 3;
    var->zlocation = CL_CENTRE;
  }

  else if (strcmp(name, "u1bot") == 0) {
    var->v = (void **)&dumpdata->u1bot;
    var->xylocation = CL_LEFT;
  }

  else if (strcmp(name, "u2bot") == 0) {
    var->v = (void **)&dumpdata->u2bot;
    var->xylocation = CL_BACK;
  }

  else if (strcmp(name, "Cd") == 0) {
    var->v = (void **)&dumpdata->Cd;
  }

  else if (strcmp(name, "flag") == 0) {
    var->v = (void **)&dumpdata->flag;
    var->ndims = 3;
    var->xylocation = CL_GRID;
    var->zlocation = CL_CENTRE;
    var->type = NC_INT;
  }

  else if (strcmp(name, "u1vh") == 0) {
    var->v = (void **)&dumpdata->u1vh;
    /* var->xylocation = CL_LEFT; */
  }

  else if (strcmp(name, "u2vh") == 0) {
    var->v = (void **)&dumpdata->u2vh;
    /* var->xylocation = CL_BACK; */
  }

  else
    found = 0;

  if (!found || dumpdata->tmode & NONE) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_wc[n];
        var->xylocation = CL_CENTRE;
        var->zlocation = CL_CENTRE;
        var->sediment = 0;
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, name) == 0) {
        var->ndims = 2;
        var->v = (void **)&dumpdata->tr_wcS[n];
        var->xylocation = CL_CENTRE;
        var->sediment = 0;
        found = 1;
        break;
      }
    }
  }

  if (!found) {
    /* Sediment tracers */
    for (n = 0; n < dumpdata->nsed; n++) {
      char name1[MAXSTRLEN];
      strcpy(name1, dumpdata->trinfo_sed[n].name);
      strcat(name1, "_sed");
      if (strcmp(name, name1) == 0) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_sed[n];
        var->xylocation = CL_CENTRE;
        var->zlocation = CL_CENTRE;
        var->sediment = 1;
        found = 1;
        break;
      }
    }
  }

  return found;
}


static void df_std_init_data(dump_data_t *dumpdata, dump_file_t *df,
                             int fid)
{
  int i;
  df_std_data_t *data = NULL;
  
  df_parse_vars(dumpdata,df,NULL,STD_ALL_VARS);
  
  // Clean up before allocating more memory
  if (df->private_data != NULL) {
    if (((df_std_data_t *)df->private_data)->vars != NULL)
      free(((df_std_data_t *)df->private_data)->vars);
    free(df->private_data);
  }
  
  data = (df_std_data_t *)malloc(sizeof(df_std_data_t));
  memset(data, 0, sizeof(df_std_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_std_var_t *)malloc(df->nvars * sizeof(df_std_var_t))) == NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_std_var_t *var = &data->vars[i];
    if (df_std_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_std_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

  if (df->diagn && DEBUG("dump"))
    dlog("dump",
         "Diagnostic tracers initialisation has been synchronised with \"%s\" output; period = %f s.\n",
         df->name, df->tinc);

}


static void df_std_get_dimids(int cdfid, df_std_var_t var,
                              int *e1id, int *e2id, int *zid)
{
  switch (var.xylocation) {
  case CL_GRID:
    *e2id = ncw_dim_id(cdfid, "j_grid");
    *e1id = ncw_dim_id(cdfid, "i_grid");
    break;

  case CL_CENTRE:
    *e2id = ncw_dim_id(cdfid, "j_centre");
    *e1id = ncw_dim_id(cdfid, "i_centre");
    break;

  case CL_LEFT:
    *e2id = ncw_dim_id(cdfid, "j_left");
    *e1id = ncw_dim_id(cdfid, "i_left");
    break;

  case CL_BACK:
    *e2id = ncw_dim_id(cdfid, "j_back");
    *e1id = ncw_dim_id(cdfid, "i_back");
    break;

  default:
    break;
  }
  if (!var.sediment) {
    switch (var.zlocation) {
    case CL_GRID:
      *zid = ncw_dim_id(cdfid, "k_grid");
      break;

    case CL_CENTRE:
      *zid = ncw_dim_id(cdfid, "k_centre");
      break;

    default:
      break;
    }
  } else {
    switch (var.zlocation) {
    case CL_GRID:
      *zid = ncw_dim_id(cdfid, "k_grid_sed");
      break;

    case CL_CENTRE:
    default:
      *zid = ncw_dim_id(cdfid, "k_centre_sed");
      break;
    }
  }

}


static void df_std_get_dimsizes(dump_file_t *df, df_std_var_t var,
                                size_t * ne1, size_t * ne2, size_t * nz)
{
  switch (var.xylocation) {
  case CL_GRID:
    *ne2 = df->nfe2;
    *ne1 = df->nfe1;
    break;

  case CL_CENTRE:
    *ne2 = df->nce2;
    *ne1 = df->nce1;
    break;

  case CL_LEFT:
    *ne2 = df->nce2;
    *ne1 = df->nfe1;
    break;

  case CL_BACK:
    *ne2 = df->nfe2;
    *ne1 = df->nce1;
    break;

  default:
    break;
  }

  if (!var.sediment) {
    switch (var.zlocation) {
    case CL_GRID:
      *nz = df->nz + 1;
      break;

    case CL_CENTRE:
      *nz = df->nz;
      break;

    default:
      break;
    }
  }

  else {
    switch (var.zlocation) {
    case CL_GRID:
      *nz = df->nz_sed + 1;
      break;

    case CL_CENTRE:
    default:
      *nz = df->nz_sed;
      break;
    }
  }
}


/* Simple grid centred file format */

#define VM_NONE   0x00
#define VM_EAST   0x01
#define VM_NORTH  0x02
#define VM_MAG    0x04
#define VM_DIRN   0x08
#define VM_NOR    0x10
#define VM_TAN    0x20

/* Structure to describe each scalar time dep variable */
typedef struct {
  int ndims;                    /* Number of spatial dimensions */
  char *name;                   /* df_variable_t name */
  char *units;                  /* df_variable_t units */
  char *long_name;              /* Descriptive name */
  void **v;                     /* Pointer to values */
  void **vi;                    /* Pointer to values associated with i
                                   axis */
  void **vj;                    /* Pointer to values associated with j
                                   axis */
  int vector_mode;              /* NONE, east, north, mag, dirn */
  nc_type fptype;               /* Data type. */
  double valid_range[2];        /* df_variable_t valid range */
  double scale;                 /* Scale factor */
  double offset;                /* Scale factor */
  char *std_name;               /* The CF standard_name. */
  char *vector_name;            /* Vector name if appropriate. */
  char *vector_components;      /* Vector component if appropriate. */
  int sediment;                 /* Flag; = 1 if variable is tracer */
} df_simple_var_t;

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_simple_var_t *vars;        /* List of simple dump file variables */
} df_simple_data_t;

static void df_simple_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid);
static void df_simple_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid);

void df_simple_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_simple_data_t *data = df->private_data;
  /* Rewind to the first next time after t */
  while (t <= (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "time", df->tout);
}

static void *df_simple_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int nid, kid, jid, iid, vid, tid, n;
  int cdfid;
  int dims[10];
  df_simple_data_t *data = NULL;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double missing_value = SIMPLE_MISSING_VALUE;
  double sb_missing_value = -99;
  double vr[2];
  int kid_sed;
  FILE *fp;
  int nc_mode;

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_simple_init_data(dumpdata, df, cdfid);
      data = (df_simple_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "time", dumpdata->t);
      df->finished = 0;
      return data;
    }
  }

  // Get the nc_mode flag
  nc_mode = get_nc_mode(df);
  
  /* create the netCDF file */
  if (ncw_create(df->name, nc_mode, &cdfid) != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "time", NC_UNLIMITED, &nid);

  nc_def_dim(cdfid, "i", df->nce1, &iid);
  nc_def_dim(cdfid, "j", df->nce2, &jid);
  nc_def_dim(cdfid, "k", df->nz, &kid);

  if (df->nz_sed > 0)
    nc_def_dim(cdfid, "k_sed", df->nz_sed, &kid_sed);

  /* time independent variables */
  dims[0] = kid;
  nc_def_var(cdfid, "zc", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  dims[0] = jid;
  dims[1] = iid;
  nc_def_var(cdfid, "xc", NC_DOUBLE, 2, dims, &vid);
  if (is_geog) {
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
  } else {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "X");
    write_text_att(cdfid, vid, "coordinate_type", "X");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, df->ilower, df->nce1,
                     df->jlower, df->nce2, CL_CENTRE);


  nc_def_var(cdfid, "yc", NC_DOUBLE, 2, dims, &vid);
  if (is_geog) {
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
  } else {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Y");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, df->ilower, df->nce1,
                     df->jlower, df->nce2, CL_CENTRE);

  if (has_proj) {
    nc_def_var(cdfid, "longitude", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");

    nc_def_var(cdfid, "latitude", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
  }

  nc_def_var(cdfid, "botz", NC_DOUBLE, 2, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Depth of sea-bed");
  nc_put_att_double(cdfid, vid, "missing_value", NC_DOUBLE, 1,
                    &sb_missing_value);
  write_text_att(cdfid, vid, "positive", "down");
  write_text_att(cdfid, vid, "coordinates", "latitude, longitude");

  /* time dependent variables */
  dims[0] = nid;
  nc_def_var(cdfid, "time", NC_DOUBLE, 1, dims, &tid);
  write_text_att(cdfid, tid, "units", df->tunit);
  write_text_att(cdfid, tid, "long_name", "Time");
  write_text_att(cdfid, tid, "coordinate_type", "time");

  if (df->nz_sed > 0) {
    dims[0] = nid;
    dims[1] = kid_sed;
    dims[2] = jid;
    dims[3] = iid;
    nc_def_var(cdfid, "zc_sed", NC_DOUBLE, 4, dims, &vid);
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate, sediment");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  df_simple_init_data(dumpdata, df, cdfid);
  data = (df_simple_data_t *)df->private_data;

  for (n = 0; n < df->nvars; n++) {
    if (!data->vars[n].sediment && data->vars[n].ndims == 2) {
      dims[1] = jid;
      dims[2] = iid;

      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 3, dims, &vid, df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, latitude, longitude, yc, xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time, yc, xc");
    } else if (!data->vars[n].sediment && data->vars[n].ndims == 3) {
      dims[1] = kid;
      dims[2] = jid;
      dims[3] = iid;
      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 4, dims, &vid, df->compress);

      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, zc, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, zc, latitude, longitude, yc, xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time, zc, yc, xc");
    }
    else if ((df->nz_sed > 0) && data->vars[n].sediment &&
             data->vars[n].ndims == 3) {
      dims[1] = kid_sed;
      dims[2] = jid;
      dims[3] = iid;
      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 4, dims,&vid,df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, zc_sed, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, zc_sed, latitude, longitude, yc, xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time, zc_sed, yc, xc");
    }

    else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }

    if (data->vars[n].fptype == NC_SHORT) {
      nc_put_att_double(cdfid, vid, "scale_factor", NC_DOUBLE,
                        1, &data->vars[n].scale);
      nc_put_att_double(cdfid, vid, "add_offset", NC_DOUBLE,
                        1, &data->vars[n].offset);
    }

    write_text_att(cdfid, vid, "units", data->vars[n].units);
    if (!data->vars[n].sediment)
      write_text_att(cdfid, vid, "long_name", data->vars[n].long_name);
    else {
      char long_name[MAXSTRLEN];
      strcpy(long_name, data->vars[n].long_name);
      strcat(long_name, " in Sediment");
      write_text_att(cdfid, vid, "long_name", long_name);
    }
    if (data->vars[n].vector_name != NULL) {
      write_text_att(cdfid, vid, "vector_name", data->vars[n].vector_name);
      write_text_att(cdfid, vid, "vector_components",
                     data->vars[n].vector_components);
    }

    vr[0] =
      (data->vars[n].valid_range[0] -
       data->vars[n].offset) / data->vars[n].scale;
    vr[1] =
      (data->vars[n].valid_range[1] -
       data->vars[n].offset) / data->vars[n].scale;
    nc_put_att_double(cdfid, vid, "valid_range", data->vars[n].fptype, 2,
                      vr);
    if (data->vars[n].fptype == NC_SHORT)
      missing_value = SIMPLE_SHORT_MISSING_VALUE;
    else
      missing_value = SIMPLE_MISSING_VALUE;
    nc_put_att_double(cdfid, vid, "missing_value", data->vars[n].fptype,
                      1, &missing_value);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", dumpdata->grid_name);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CMR/Timeseries");
  if (dumpdata->runno >= 0)
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);

  nc_enddef(cdfid);

  df_simple_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}

static void *df_simple_cf_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int nid, kid, jid, iid, vid, tid, n;
  int cdfid;
  int dims[10];
  df_simple_data_t *data = NULL;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double missing_value = SIMPLE_MISSING_VALUE;
  double sb_missing_value = -99;
  double vr[2];
  int kid_sed;
  FILE *fp;
  int nc_mode;
                                
  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_simple_init_data(dumpdata, df, cdfid);
      data = (df_simple_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "time", dumpdata->t);
      df->finished = 0;
      return data;
    }
  }

  // Get the nc_mode flag
  nc_mode = get_nc_mode(df);

  /* create the netCDF file */
  if (ncw_create(df->name, nc_mode, &cdfid) != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "time", NC_UNLIMITED, &nid);

  nc_def_dim(cdfid, "i", df->nce1, &iid);
  nc_def_dim(cdfid, "j", df->nce2, &jid);
  nc_def_dim(cdfid, "k", df->nz, &kid);

  if (df->nz_sed > 0)
    nc_def_dim(cdfid, "k_sed", df->nz_sed, &kid_sed);

  /* time independent variables */
  dims[0] = kid;
  nc_def_var(cdfid, "zc", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", "m");
  write_text_att(cdfid, vid, "positive", "up");
  write_text_att(cdfid, vid, "long_name", "Z coordinate");
  write_text_att(cdfid, vid, "axis", "Z");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  dims[0] = jid;
  dims[1] = iid;
  /*UR 2/2010 if this is geographic stick with 
   nc_def_var(cdfid, "xc", NC_DOUBLE, 2, dims, &vid);
  */
  if (is_geog) {
    nc_def_var(cdfid, "longitude", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "standard_name", "longitude");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
  } else {
  	nc_def_var(cdfid, "xc", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "X");
    write_text_att(cdfid, vid, "axis", "X");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, df->ilower, df->nce1,
                     df->jlower, df->nce2, CL_CENTRE);


  nc_def_var(cdfid, "latitude", NC_DOUBLE, 2, dims, &vid);
  if (is_geog) {
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "standard_name", "latitude");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
  } else {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Y");
    write_text_att(cdfid, vid, "axis", "Y");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);
  write_analytic_att(dumpdata, cdfid, vid, df->ilower, df->nce1,
                     df->jlower, df->nce2, CL_CENTRE);

  if (has_proj && !is_geog) {
    nc_def_var(cdfid, "longitude", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "standard_name", "longitude");

    nc_def_var(cdfid, "latitude", NC_DOUBLE, 2, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "standard_name", "latitude");
  }

  nc_def_var(cdfid, "botz", NC_DOUBLE, 2, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Depth of sea-bed");
  write_text_att(cdfid, vid, "standard_name", "depth");
  nc_put_att_double(cdfid, vid, "missing_value", NC_DOUBLE, 1,
                    &sb_missing_value);
  write_text_att(cdfid, vid, "positive", "down");
  write_text_att(cdfid, vid, "coordinates", "latitude longitude");

  /*UR outside value for open water boundaries, TODO - replace with constant definition */
  write_text_att(cdfid, vid, "outside", "9999");

  /* time dependent variables */
  dims[0] = nid;
  nc_def_var(cdfid, "time", NC_DOUBLE, 1, dims, &tid);
  write_text_att(cdfid, tid, "units", df->tunit);
  write_text_att(cdfid, tid, "long_name", "Time");
  write_text_att(cdfid, tid, "standard_name", "time");
  write_text_att(cdfid, tid, "coordinate_type", "time");

  if (df->nz_sed > 0) {
    dims[0] = nid;
    dims[1] = kid_sed;
    dims[2] = jid;
    dims[3] = iid;
    nc_def_var(cdfid, "zc_sedim", NC_DOUBLE, 4, dims, &vid);
    write_text_att(cdfid, vid, "units", "m");
    write_text_att(cdfid, vid, "positive", "up");
    write_text_att(cdfid, vid, "long_name", "Z coordinate of sediment centre");
    // write_text_att(cdfid, vid, "axis", "Z");
    /* Add 1d var so simple tools can handle sediment variables */
    dims[0] = kid_sed;
    nc_def_var(cdfid, "zcsed", NC_DOUBLE, 1, dims, &vid);
    write_text_att(cdfid, vid, "units", "m");
    write_text_att(cdfid, vid, "positive", "up");
    write_text_att(cdfid, vid, "long_name", "Z coordinate, sediment");
    write_text_att(cdfid, vid, "axis", "Z");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  df_simple_init_data(dumpdata, df, cdfid);
  data = (df_simple_data_t *)df->private_data;

  dims[0] = nid;
  for (n = 0; n < df->nvars; n++) {
    if (!data->vars[n].sediment && data->vars[n].ndims == 2) {
      dims[1] = jid;
      dims[2] = iid;

      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 3, dims,&vid,df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time latitude longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time yc xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time yc xc");
    } else if (!data->vars[n].sediment && data->vars[n].ndims == 3) {
      dims[1] = kid;
      dims[2] = jid;
      dims[3] = iid;
      ncw_def_var2(df->name,cdfid, data->vars[n].name, data->vars[n].fptype, 4, dims,&vid,df->compress);

      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time zc latitude longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time zc yc xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time zc yc xc");
    }
    else if ((df->nz_sed > 0) && data->vars[n].sediment &&
             data->vars[n].ndims == 3) {
      dims[1] = kid_sed;
      dims[2] = jid;
      dims[3] = iid;
      ncw_def_var2(df->name,cdfid, data->vars[n].name, data->vars[n].fptype, 4, dims, &vid,df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time zcsed latitude longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time zcsed yc xc");
      else
        write_text_att(cdfid, vid, "coordinates", "time zcsed yc xc");
    }

    else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }

    if (data->vars[n].fptype == NC_SHORT) {
      nc_put_att_double(cdfid, vid, "scale_factor", NC_DOUBLE,
                        1, &data->vars[n].scale);
      nc_put_att_double(cdfid, vid, "add_offset", NC_DOUBLE,
                        1, &data->vars[n].offset);
    }

    write_text_att(cdfid, vid, "units", data->vars[n].units);
    if (!data->vars[n].sediment)
      write_text_att(cdfid, vid, "long_name", data->vars[n].long_name);
    else {
      char long_name[MAXSTRLEN];
      strcpy(long_name, data->vars[n].long_name);
      strcat(long_name, " in Sediment");
      write_text_att(cdfid, vid, "long_name", long_name);
    }

    if (data->vars[n].std_name != NULL && strlen(data->vars[n].std_name) > 0) {
       write_text_att(cdfid, vid, "standard_name", data->vars[n].std_name);
    }

    if (data->vars[n].vector_name != NULL) {
      write_text_att(cdfid, vid, "vector_name", data->vars[n].vector_name);
      write_text_att(cdfid, vid, "vector_components",
                     data->vars[n].vector_components);
    }

    vr[0] =
      (data->vars[n].valid_range[0] -
       data->vars[n].offset) / data->vars[n].scale;
    vr[1] =
      (data->vars[n].valid_range[1] -
       data->vars[n].offset) / data->vars[n].scale;
    nc_put_att_double(cdfid, vid, "valid_range", data->vars[n].fptype, 2,
                      vr);
    if (data->vars[n].fptype == NC_SHORT)
      missing_value = SIMPLE_SHORT_MISSING_VALUE;
    else
      missing_value = SIMPLE_MISSING_VALUE;
    nc_put_att_double(cdfid, vid, "missing_value", data->vars[n].fptype,
                      1, &missing_value);

    /*UR June 2010 - not the right place, but there is no other right place either */
    if (strcmp(data->vars[n].name, "eta") == 0)
    {
    	write_text_att(cdfid, vid, "positive", "up");
    }
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", dumpdata->grid_name);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CF-1.0");
  if (dumpdata->runno >= 0)
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);

  nc_enddef(cdfid);

  df_simple_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}

static void df_simple_write(dump_data_t *dumpdata, dump_file_t *df,
                            double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = t;
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;
  int fid = data->fid;

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;

  tm_change_time_units(dumpdata->timeunit, df->tunit, &newt,
                         1);
  nc_put_vara_double(fid, ncw_var_id(fid, "time"), start, count, &newt);
  count[0] = 1L;

  if (df->nz_sed != 0) {
    count[1] = df->nz_sed;
    count[2] = df->nce2;
    count[3] = df->nce1;
    nc_put_vara_double(fid, ncw_var_id(fid, "zc_sedim"), start, count,
                       dumpdata->cellz_sed[0][0]);
  }

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_simple_var_t *var = &data->vars[n];
    void *p = var->v;
    if (var->ndims == 2) {
      start[1] = df->jlower;
      start[2] = df->ilower;
      count[1] = df->nce2;
      count[2] = df->nce1;
      if (var->vector_mode == VM_NONE)
        df_simple_writesub_2d(dumpdata, df, n,
                              start, count, (*((double ***)p)));
      else
        df_simple_writevec_2d(dumpdata, df, n);
    } else if (var->ndims == 3) {
      
      if (var->sediment) {
        start[1] = 0;
        count[1] = df->nz_sed;
      } else {
        start[1] = df->klower;
        count[1] = df->nz;
      }
      start[2] = df->jlower;
      start[3] = df->ilower;
      count[2] = df->nce2;
      count[3] = df->nce1;
      if (var->vector_mode == VM_NONE) {
        df_simple_writesub_3d(dumpdata, df, n,
                              start, count, (*((double ****)p)), start[1]);
      }
      else
        df_simple_writevec_3d(dumpdata, df, n);
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

static void df_simple_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}


static void df_simple_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[2];
  size_t count[2];
  df_std_data_t *data = (df_std_data_t *)df->private_data;
  int fid = data->fid;
  int i, j;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double **botz = d_alloc_2d(dumpdata->nce1, dumpdata->nce2);

  set_longitude(dumpdata, df, 1);
  start[0] = df->klower;
  start[1] = 0;
  count[0] = df->nz;
  count[1] = 0;
  nc_d_writesub_1d(fid, ncw_var_id(fid, "zc"), start, count,
                   dumpdata->cellz);
  if (df->nz_sed) {
    start[0] = 0;
    count[0] = df->nz_sed;
    nc_d_writesub_1d(fid, ncw_var_id(fid, "zcsed"), start, count,
		     dumpdata->cellzcsed);
  }

  start[0] = df->jlower;
  start[1] = df->ilower;
  count[0] = df->nce2;
  count[1] = df->nce1;
  if (!is_geog) {
    nc_d_writesub_2d(fid, ncw_var_id(fid, "xc"), start, count,
		     dumpdata->cellx);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "yc"), start, count,
		     dumpdata->celly);
  }
  for (j = 0; j < dumpdata->nce2; j++)
    for (i = 0; i < dumpdata->nce1; i++)
      botz[j][i] = -dumpdata->botz[j][i];
  nc_d_writesub_2d(fid, ncw_var_id(fid, "botz"), start, count, botz);
  d_free_2d(botz);

  /* Write lat/long if native, or convert from projection */
  if (is_geog) {
    nc_d_writesub_2d(fid, ncw_var_id(fid, "longitude"), start, count,
                     dumpdata->cellx);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "latitude"), start, count,
                     dumpdata->celly);

  } else if (has_proj) {
    char *args[256];
    int nargs = parseline(strdup(projection), args, 256);
    map_proj_t *mp = mp_init(nargs, args);
    double **lon = d_alloc_2d(dumpdata->nce1, dumpdata->nce2);
    double **lat = d_alloc_2d(dumpdata->nce1, dumpdata->nce2);

    for (j = 0; j < dumpdata->nce2; j++)
      for (i = 0; i < dumpdata->nce1; i++)
        mp_inverse(mp, dumpdata->cellx[j][i], dumpdata->celly[j][i],
                   &lat[j][i], &lon[j][i]);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "longitude"), start, count, lon);
    nc_d_writesub_2d(fid, ncw_var_id(fid, "latitude"), start, count, lat);
    d_free_2d(lat);
    d_free_2d(lon);
  }

  ncw_sync(df->name,fid);
  set_longitude(dumpdata, df, 0);
}

static int df_simple_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
                                 char *name, df_simple_var_t *var)
{
  master_t *master =  dumpdata->master;
  int found = 1;
  int n = 0;
  int nargs;
  char *args[256];

  var->ndims = 2;
  var->units = NULL;
  var->long_name = NULL;
  var->std_name = NULL;
  var->v = NULL;
  var->vi = NULL;
  var->vj = NULL;
  var->vector_mode = VM_NONE;
  var->fptype = nc_get_default_type(df->bpv);
  var->valid_range[0] = 0.0;
  var->valid_range[1] = 0.0;
  var->scale = 1.0;
  var->offset = 0.0;
  var->vector_name = NULL;
  var->vector_components = NULL;
  var->sediment = 0;

  /* If the variable name has an alternative name, then use this. */
  nargs = decode_args(name, args, 256);
  if (nargs == 1)
    var->name = args[0];
  else if (nargs == 2)
    var->name = args[1];
  else
    hd_quit("df_simple_get_varinfo: Invalid output variable '%s'.\n",name);

  if (strcmp(args[0], "eta") == 0) {
    var->v = (void **)&dumpdata->eta;
    var->units = dumpdata->lenunit;
    var->long_name = "Surface elevation";
    var->std_name = "sea_surface_height_above_sea_level";
    var->valid_range[0] = -master->etamax;
    var->valid_range[1] = master->etamax;
  }

  else if (strcmp(args[0], "patm") == 0) {
    var->v = (void **)&dumpdata->patm;
    var->units = "Pa";
    var->long_name = "Atmospheric pressure";
    var->std_name = "surface_air_pressure";
    var->valid_range[0] = 91000.0;
    var->valid_range[1] = 111000.0;
  }

  else if (strcmp(args[0], "Cd") == 0) {
    var->v = (void **)&dumpdata->Cd;
    var->units = "";
    var->long_name = "Bottom drag coefficient";
    var->std_name = "bottom_drag_coefficient";  
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 1.0;
  }

  else if (strcmp(args[0], "Kz") == 0) {
    var->v = (void **)&dumpdata->Kz;
    var->ndims = 3;
    var->units = "m2s-1";
    var->long_name = "Kz";
    var->std_name = "sea_water_vertical_diffusivity"; 
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 0.1;
  }

  else if (strcmp(args[0], "Vz") == 0) {
    var->v = (void **)&dumpdata->Vz;
    var->ndims = 3;
    var->units = "m2s-1";
    var->long_name = "Vz";
    var->std_name = "sea_water_vertical_eddy_viscocity"; 
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 0.1;
  }

  else if (strcmp(args[0], "dens") == 0) {
    var->v = (void **)&dumpdata->dens;
    var->ndims = 3;
    var->units = "kg m-3";
    var->long_name = "Density";
    var->std_name = "sea_water_density";
    var->valid_range[0] = 980.0;
    var->valid_range[1] = 1070.0;
  }

  else if (strcmp(args[0], "dens_0") == 0) {
    var->v = (void **)&dumpdata->dens_0;
    var->ndims = 3;
    var->units = "kg m-3";
    var->long_name = "Potential density";
    var->std_name = "sea_water_potential_density";
    var->valid_range[0] = 980.0;
    var->valid_range[1] = 1070.0;
  }
  else if (strcmp(args[0], "w") == 0) {
    var->v = (void **)&dumpdata->w;
    var->ndims = 3;
    var->units = "ms-1";
    var->long_name = "Vertical current";
    var->std_name = "upward_sea_water_velocity";
    var->valid_range[0] = -0.1;
    var->valid_range[1] = 0.1;
  }

  else if (strcmp(args[0], "uav") == 0) {
    var->vector_mode = VM_EAST;
    var->vi = (void **)&dumpdata->u1av;
    var->vj = (void **)&dumpdata->u2av;
    var->units = "ms-1";
    var->long_name = "Eastward depth averaged current";
    var->std_name = "eastward_sea_water_velocity";
    if (master->fatal & VEL2D) {
      var->valid_range[0] = -master->velmax2d;
      var->valid_range[1] = master->velmax2d;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Depth average currents";
    var->vector_components = "uav vav";
  }

  else if (strcmp(args[0], "vav") == 0) {
    var->vector_mode = VM_NORTH;
    var->vi = (void **)&dumpdata->u1av;
    var->vj = (void **)&dumpdata->u2av;
    var->units = "ms-1";
    var->long_name = "Northward depth averaged current";
    var->std_name = "northward_sea_water_velocity";
    if (master->fatal & VEL2D) {
      var->valid_range[0] = -master->velmax2d;
      var->valid_range[1] = master->velmax2d;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Depth average currents";
    var->vector_components = "uav vav";
  }

  else if (strcmp(args[0], "avg_speed") == 0) {
    var->vector_mode = VM_MAG;
    var->vi = (void **)&dumpdata->u1av;
    var->vj = (void **)&dumpdata->u2av;
    var->units = "ms-1";
    var->long_name = "Depth averaged current speed";
/*    var->std_name = "barotropic_sea_water_speed"; */
    if (master->fatal & VEL2D) {
      var->valid_range[0] = -master->velmax2d;
      var->valid_range[1] = master->velmax2d;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
  }

  else if (strcmp(args[0], "avg_dir") == 0) {
    var->vector_mode = VM_DIRN;
    var->vi = (void **)&dumpdata->u1av;
    var->vj = (void **)&dumpdata->u2av;
    var->units = "degrees";
    var->long_name = "Depth averaged current direction";
/*    var->std_name = "barotropic_direction_of_sea_water";  */
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 360.0;
  }

  else if (strcmp(args[0], "u") == 0) {
    var->vector_mode = VM_EAST;
    var->vi = (void **)&dumpdata->u1;
    var->vj = (void **)&dumpdata->u2;
    var->units = "ms-1";
    var->long_name = "Eastward current";
    var->std_name = "eastward_sea_water_velocity";
    var->ndims = 3;
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Currents";
    var->vector_components = "u v";
  }

  else if (strcmp(args[0], "v") == 0) {
    var->vector_mode = VM_NORTH;
    var->vi = (void **)&dumpdata->u1;
    var->vj = (void **)&dumpdata->u2;
    var->units = "ms-1";
    var->long_name = "Northward current";
    var->std_name = "northward_sea_water_velocity";
    var->ndims = 3;
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Currents";
    var->vector_components = "u v";
  }

  else if (strcmp(args[0], "current_speed") == 0) {
    var->vector_mode = VM_MAG;
    var->vi = (void **)&dumpdata->u1;
    var->vj = (void **)&dumpdata->u2;
    var->units = "ms-1";
    var->long_name = "Current speed";
/*    var->std_name = "baroclinic_sea_water_speed"; */
    var->ndims = 3;
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
  }

  else if (strcmp(args[0], "current_dir") == 0) {
    var->vector_mode = VM_DIRN;
    var->vi = (void **)&dumpdata->u1;
    var->vj = (void **)&dumpdata->u2;
    var->units = "degrees";
    var->long_name = "Current direction";
/*    var->std_name = "baroclinic_direction_of_sea_water_speed"; */
    var->ndims = 3;
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 360.0;
  }

  else if (strcmp(args[0], "wind_u") == 0) {
    var->vector_mode = VM_EAST;
    var->vi = (void **)&dumpdata->wind1;
    var->vj = (void **)&dumpdata->wind2;
    var->units = "Nm-2";
    var->long_name = "Eastward wind stress";
    var->std_name = "surface_eastward_wind_stress"; 
    var->valid_range[0] = -4.5;
    var->valid_range[1] = 4.5;
    var->vector_name = "Wind stress";
    var->vector_components = "wind_u wind_v";
  }

  else if (strcmp(args[0], "wind_v") == 0) {
    var->vector_mode = VM_NORTH;
    var->vi = (void **)&dumpdata->wind1;
    var->vj = (void **)&dumpdata->wind2;
    var->units = "Nm-2";
    var->long_name = "Northward wind stress";
    var->std_name = "surface_northward_wind_stress"; 
    var->valid_range[0] = -4.5;
    var->valid_range[1] = 4.5;
    var->vector_name = "Wind stress";
    var->vector_components = "wind_u wind_v";
  }

  else if (strcmp(args[0], "wind_mag") == 0) {
    var->vector_mode = VM_MAG;
    var->vi = (void **)&dumpdata->wind1;
    var->vj = (void **)&dumpdata->wind2;
    var->units = "Nm-2";
    var->long_name = "Wind stress magnitude";
    var->std_name = "surface_wind_stress_speed";
    var->valid_range[0] = -50.0;
    var->valid_range[1] = 50.0;
  }

  else if (strcmp(args[0], "wind_dir") == 0) {
    var->vector_mode = VM_DIRN;
    var->vi = (void **)&dumpdata->wind1;
    var->vj = (void **)&dumpdata->wind2;
    var->units = "degrees";
    var->long_name = "Wind stress direction";
    var->std_name = "surface_direction_of_wind_stress"; 
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 360.0;
  }

  else if (strcmp(args[0], "bottom_u") == 0) {
    var->vector_mode = VM_EAST;
    var->vi = (void **)&dumpdata->u1bot;
    var->vj = (void **)&dumpdata->u2bot;
    var->units = "ms-1";
    var->long_name = "Eastward bottom current deviation";
/*    var->std_name = "bottom_eastward_sea_water_velocity"; */
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Bottom current";
    var->vector_components = "bottom_u bottom_v";
  }

  else if (strcmp(args[0], "bottom_v") == 0) {
    var->vector_mode = VM_NORTH;
    var->vi = (void **)&dumpdata->u1bot;
    var->vj = (void **)&dumpdata->u2bot;
    var->units = "ms-1";
    var->long_name = "Northward bottom current deviation";
/*    var->std_name = "bottom_northward_sea_water_velocity"; */
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
    var->vector_name = "Bottom current";
    var->vector_components = "bottom_u bottom_v";
  }

  else if (strcmp(args[0], "bottom_speed") == 0) {
    var->vector_mode = VM_MAG;
    var->vi = (void **)&dumpdata->u1bot;
    var->vj = (void **)&dumpdata->u2bot;
    var->units = "ms-1";
    var->long_name = "Bottom current deviation speed";
/*    var->std_name = "bottom_sea_water_speed"; */
    if (master->fatal & VEL3D) {
      var->valid_range[0] = -master->velmax;
      var->valid_range[1] = master->velmax;
    } else {
      var->valid_range[0] = -100.0;
      var->valid_range[1] = 100.0;
    }
  }

  else if (strcmp(args[0], "bottom_dir") == 0) {
    var->vector_mode = VM_DIRN;
    var->vi = (void **)&dumpdata->u1bot;
    var->vj = (void **)&dumpdata->u2bot;
    var->units = "degrees";
    var->long_name = "Bottom current deviation direction";
/*    var->std_name = "bottom_direction_of_sea_water"; */
    var->valid_range[0] = 0.0;
    var->valid_range[1] = 360.0;
  }
  else
    found = 0;

  if (!found || dumpdata->tmode & NONE) {
    // var->vector_mode = VM_NONE;
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, args[0]) == 0) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_wc[n];
        var->units = dumpdata->trinfo_3d[n].units;
        var->long_name = dumpdata->trinfo_3d[n].long_name;
        if (strlen(dumpdata->trinfo_3d[n].std_name) != 0) {
           var->std_name = dumpdata->trinfo_3d[n].std_name;
        }
        if (strlen(dumpdata->trinfo_3d[n].vector_name) != 0) {
	  var->vector_name = dumpdata->trinfo_3d[n].vector_name;
        }
	/* If this tracer is a vector component, then set attributes for rotation */
        if (strlen(dumpdata->trinfo_3d[n].vector_components) != 0) {
	  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
	  int m;
	  var->vector_components = dumpdata->trinfo_3d[n].vector_components;
	  strcpy(buf, dumpdata->trinfo_3d[n].vector_components);
	  m = parseline(buf, fields, MAXNUMARGS);
	  if (strcmp(fields[0], args[0]) == 0) {
	    if ((m = tracer_find_index(fields[1], dumpdata->ntr, dumpdata->trinfo_3d)) >= 0) {
	      var->vector_mode = VM_EAST;
	      var->vi = (void **)&dumpdata->tr_wc[n];
	      var->vj = (void **)&dumpdata->tr_wc[m];
	    }
	  } else {
	    if ((m = tracer_find_index(fields[0], dumpdata->ntr, dumpdata->trinfo_3d)) >= 0) {
	      var->vector_mode = VM_NORTH;
	      var->vi = (void **)&dumpdata->tr_wc[m];
	      var->vj = (void **)&dumpdata->tr_wc[n];
	    }
	  }
        }
        var->valid_range[0] = dumpdata->trinfo_3d[n].valid_range_wc[0];
        var->valid_range[1] = dumpdata->trinfo_3d[n].valid_range_wc[1];

        var->sediment = 0;
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, args[0]) == 0) {
        var->ndims = 2;
        var->v = (void **)&dumpdata->tr_wcS[n];
        var->units = dumpdata->trinfo_2d[n].units;
        var->long_name = dumpdata->trinfo_2d[n].long_name;
        if (strlen(dumpdata->trinfo_2d[n].std_name) != 0) {
           var->std_name = dumpdata->trinfo_2d[n].std_name;
        }
        if (strlen(dumpdata->trinfo_2d[n].vector_name) != 0) {
	  var->vector_name = dumpdata->trinfo_2d[n].vector_name;
        }
	/* If this tracer is a vector component, then set attributes for rotation */
        if (strlen(dumpdata->trinfo_2d[n].vector_components) != 0) {
	  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
	  int m;
	  var->vector_components = dumpdata->trinfo_2d[n].vector_components;
	  strcpy(buf, dumpdata->trinfo_2d[n].vector_components);
	  m = parseline(buf, fields, MAXNUMARGS);
	  if (strcmp(fields[0], args[0]) == 0) {
	    if ((m = tracer_find_index(fields[1], dumpdata->ntrS, dumpdata->trinfo_2d)) >= 0) {
	      var->vector_mode = VM_EAST;
	      var->vi = (void **)&dumpdata->tr_wcS[n];
	      var->vj = (void **)&dumpdata->tr_wcS[m];
	    }
	  } else {
	    if ((m = tracer_find_index(fields[0], dumpdata->ntrS, dumpdata->trinfo_2d)) >= 0) {
	      var->vector_mode = VM_NORTH;
	      var->vi = (void **)&dumpdata->tr_wcS[m];
	      var->vj = (void **)&dumpdata->tr_wcS[n];
	    }
	  }
        }
        var->valid_range[0] = dumpdata->trinfo_2d[n].valid_range_wc[0];
        var->valid_range[1] = dumpdata->trinfo_2d[n].valid_range_wc[1];
        var->sediment = 0;
        found = 1;
        break;
      }
    }
  }

  if (!found) {
    /* Sediment tracers */
    for (n = 0; n < dumpdata->nsed; n++) {
      char name[MAXSTRLEN];
      char long_name[MAXSTRLEN];
      strcpy(name, dumpdata->trinfo_sed[n].name);
      strcat(name, "_sed");
      if (strcmp(name, args[0]) == 0) {
        var->ndims = 3;
        var->v = (void **)&dumpdata->tr_sed[n];
        var->units = dumpdata->trinfo_sed[n].units;
        strcpy(long_name, dumpdata->trinfo_sed[n].long_name);
        var->long_name = dumpdata->trinfo_sed[n].long_name;
        if (strlen(dumpdata->trinfo_sed[n].std_name) != 0) {
           var->std_name = dumpdata->trinfo_sed[n].std_name;
        }
        var->valid_range[0] = dumpdata->trinfo_sed[n].valid_range_sed[0];
        var->valid_range[1] = dumpdata->trinfo_sed[n].valid_range_sed[1];
        var->sediment = 1;
        found = 1;
        break;
      }
    }
  }

  if (found && var->fptype == NC_SHORT) {
    var->scale = ((var->valid_range[1] - var->valid_range[0])) / 64000.0;
    var->offset = ((var->valid_range[1] - var->valid_range[0])) / 2.0
      + var->valid_range[0];
  }

  return found;
}

static void df_simple_init_data(dump_data_t *dumpdata, dump_file_t *df,
                                int fid)
{
  int i;
  df_simple_data_t *data = NULL;

  df_parse_vars(dumpdata,df,NULL,SIMPLE_ALL_VARS);
  
  // Clean up before allocating more memory
  if (df->private_data != NULL) {
    if (((df_simple_data_t *)df->private_data)->vars != NULL)
      free(((df_simple_data_t *)df->private_data)->vars);
    free(df->private_data);
  }
  
  data = (df_simple_data_t *)malloc(sizeof(df_simple_data_t));
  memset(data, 0, sizeof(df_simple_data_t));
  df->private_data = data;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_simple_var_t *)malloc(df->nvars * sizeof(df_simple_var_t))) ==
      NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_simple_var_t *var = &data->vars[i];
    if (df_simple_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_simple_init_data: Unknown variable '%s'.", df->vars[i]);
  }
  data->fid = fid;

  /* Set next record to zero */
  data->nextrec = 0;

  if (df->diagn && DEBUG("dump"))
    dlog("dump",
         "Diagnostic tracers initialisation has been synchronised with \"%s\" output; period = %f s.\n",
         df->name, df->tinc);

}

static void df_simple_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid)
{
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;
  df_simple_var_t *var = &data->vars[vid];
  double **c = d_alloc_2d(dumpdata->nce1, dumpdata->nce2);
  double **vi = *((double ***)var->vi);
  double **vj = *((double ***)var->vj);
  size_t start[4];
  size_t count[4];

  start[0] = data->nextrec;
  start[1] = df->jlower;
  start[2] = df->ilower;
  count[0] = 1;
  count[1] = df->nce2;
  count[2] = df->nce1;

  vector_component(dumpdata, vi, vj, df->ilower, df->jlower,
                   df->nce1, df->nce2, var->vector_mode, c);
  df_simple_writesub_2d(dumpdata, df, vid, start, count, c);
  d_free_2d(c);
}

static void df_simple_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid)
{
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;
  df_simple_var_t *var = &data->vars[vid];
  double ***c = d_alloc_3d(dumpdata->nce1, dumpdata->nce2, dumpdata->nz);
  double ***vi = *((double ****)var->vi);
  double ***vj = *((double ****)var->vj);
  int klower = df->klower;
  int kupper = klower + df->nz - 1;
  int k, nk;
  size_t start[4];
  size_t count[4];

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = df->jlower;
  start[3] = df->ilower;
  count[0] = 1;
  count[1] = df->nz;
  count[2] = df->nce2;
  count[3] = df->nce1;
  for (k = klower, nk = 0; k <= kupper; ++k, ++nk)
    vector_component(dumpdata, vi[k], vj[k], df->ilower, df->jlower,
                     df->nce1, df->nce2, var->vector_mode, c[nk]);

  df_simple_writesub_3d(dumpdata, df, vid, start, count, c, df->klower);
  d_free_3d(c);
}

static void df_simple_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, size_t * start, size_t * count,
                                  double **values)
{
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;
  df_simple_var_t *var = &data->vars[vid];
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  double **nvals = NULL;
  size_t nstart[4];
  unsigned int i, j;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  unsigned int nzm1 = dumpdata->nz - 1;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 2);
  nvals = d_alloc_2d(count[offset + 1], count[offset]);
  for (j = 0; j < count[offset]; ++j) {
    for (i = 0; i < count[offset + 1]; ++i) {
      unsigned int oi = start[offset + 1] + i;
      unsigned int oj = start[offset] + j;
      if (var->fptype == NC_SHORT) {
        if (dumpdata->flag[nzm1][oj][oi] & (SOLID | OUTSIDE) &&
	    df->landfill == locate_landfill_function("default"))
          nvals[j][i] = SIMPLE_SHORT_MISSING_VALUE;
        else
          nvals[j][i] = (values[oj][oi] - var->offset) / var->scale;

      } else {
        if (dumpdata->flag[nzm1][oj][oi] & (SOLID | OUTSIDE) &&
	    df->landfill == locate_landfill_function("default"))
          nvals[j][i] = SIMPLE_MISSING_VALUE;
        else {
          nvals[j][i] = values[oj][oi];
	}
      }
    }
  }

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;
  nc_put_vara_double(fid, varid, nstart, count, nvals[0]);

  d_free_2d(nvals);
}


/* Write out a double 3d subsection.
 */
static void df_simple_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, size_t * start, size_t * count,
                                  double ***values, int klower)
{
  df_simple_data_t *data = (df_simple_data_t *)df->private_data;
  df_simple_var_t *var = &data->vars[vid];
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  double ***nvals = NULL;
  size_t nstart[4];
  int fid = data->fid;
  unsigned int varid = ncw_var_id(fid, df->vars[vid]);
  unsigned int i, j, k;


  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 3);
  nvals = d_alloc_3d(count[offset + 2], count[offset + 1], count[offset]);
  if (!var->sediment) {         /* watercolumn */
    for (k = 0; k < count[offset]; ++k) {
      for (j = 0; j < count[offset + 1]; ++j) {
        for (i = 0; i < count[offset + 2]; ++i) {
          int oi = start[offset + 2] + i;
          int oj = start[offset + 1] + j;
          int ok = start[offset] + k;
          if (var->fptype == NC_SHORT) {
            if (dumpdata->flag[klower + k][oj][oi] & (SOLID | OUTSIDE) &&
		df->landfill == locate_landfill_function("default"))
              nvals[k][j][i] = SIMPLE_SHORT_MISSING_VALUE;
            else
              nvals[k][j][i] =
                (values[ok][oj][oi] - var->offset) / var->scale;
          } else {
            if (dumpdata->flag[klower + k][oj][oi] & (SOLID | OUTSIDE) &&
		df->landfill == locate_landfill_function("default"))
              nvals[k][j][i] = SIMPLE_MISSING_VALUE;
            else
              nvals[k][j][i] = values[ok][oj][oi];
          }
        }
      }
    }
  } else {                      /* sediment */
    for (k = 0; k < count[offset]; ++k) {
      for (j = 0; j < count[offset + 1]; ++j) {
        for (i = 0; i < count[offset + 2]; ++i) {
          int oi = start[offset + 2] + i;
          int oj = start[offset + 1] + j;
          int ok = start[offset] + k;
          if (var->fptype == NC_SHORT) {
            if (dumpdata->
                flag[dumpdata->nz - 1][oj][oi] & (SOLID | OUTSIDE) &&
		df->landfill == locate_landfill_function("default"))
              nvals[k][j][i] = SIMPLE_SHORT_MISSING_VALUE;
            else
              nvals[k][j][i] =
                (values[ok][oj][oi] - var->offset) / var->scale;
          } else {
            if (dumpdata->
                flag[dumpdata->nz - 1][oj][oi] & (SOLID | OUTSIDE) &&
		df->landfill == locate_landfill_function("default"))
              nvals[k][j][i] = SIMPLE_MISSING_VALUE;
            else
              nvals[k][j][i] = values[ok][oj][oi];
          }
        }
      }
    }
  }
  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  nc_put_vara_double(fid, varid, nstart, count, nvals[0][0]);

  d_free_3d(nvals);
}



static int read_subsection_range(FILE * fp, char *key, int nc, int nf,
                                  int *sslower, int *ssnc, int *ssnf)
{
  int lower, upper;

  prm_set_errfn(hd_silent_warn);
  if (prm_skip_to_end_of_key(fp, key) > 0) {
    if (fscanf(fp, "%d %d", &lower, &upper) == 2) {
      if ((lower >= 0) && (lower < nc) && (upper >= 0) && (upper < nc)) {
        if (sslower != NULL)
          *sslower = lower;
        if (ssnc != NULL)
          *ssnc = (upper - lower) + 1;
        if (ssnf != NULL)
          *ssnf = *ssnc + (nf - nc);
      } else {
        hd_quit
          ("read_subsection_range: Invalid subsection range specified for '%s', valid range is %d to %d\n",
           key, 0, nc - 1);
      }
      return(0);
    } else
      hd_quit
        ("read_subsection_range: Invalid subsection range specifed for '%s', valid range is %d to %d\n",
         key, 0, nc - 1);
  }
  return(1);
}



/* figure out the floating point type we are using */
nc_type nc_get_default_type(int bpv)
{

  switch (bpv) {
  case 8:
    return NC_DOUBLE;

  case 4:
    return NC_FLOAT;

  case 2:
    return NC_SHORT;
  }

  return NC_DOUBLE;
}

/* Write out a short 1d subsection.
 */
#if 0
void nc_s_writesub_1d(int fid, int varid, size_t * start,
                             size_t * count, short *values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  size_t nstart[4];
  unsigned int i;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 1);

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  nc_put_vara_short(fid, varid, nstart, count, &values[start[offset]]);
}
#endif

/* Write out an int 1d subsection. */
void nc_i_writesub_1d(int fid, int varid, size_t * start,
                             size_t * count, int *values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  size_t nstart[4];
  unsigned int i;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 1);

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  nc_put_vara_int(fid, varid, nstart, count, &values[start[offset]]);
}


/* Write out a long 3d subsection.
 */
void nc_l_writesub_3d(int fid, int varid, size_t * start,
                             size_t * count, unsigned long ***values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  long ***nvals = NULL;
  size_t nstart[4];
  unsigned int i, j, k;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 3);
  nvals =
    (long ***)l_alloc_3d(count[offset + 2], count[offset + 1],
                                  count[offset]);
  for (k = 0; k < count[offset]; ++k)
    for (j = 0; j < count[offset + 1]; ++j)
      for (i = 0; i < count[offset + 2]; ++i)
        nvals[k][j][i] =
          values[start[offset] + k][start[offset + 1] +
                                    j][start[offset + 2] + i];

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  nc_put_vara_long(fid, varid, nstart, count, nvals[0][0]);

  l_free_3d((long ***)nvals);
}

/* Write out a double 1d subsection.
 */
void nc_d_writesub_1d(int fid, int varid, size_t * start,
                             size_t * count, double *values)
{
  /*UR-FIX for icc */
  int nd = 0;
  int offset = 0;               /* Has a record dimension = 1, else 0 */
  size_t nstart[4];
  int i;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 1);
  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  nc_put_vara_double(fid, varid, nstart, count, &values[start[offset]]);
}


/* Write out a int 2d subsection.
 */
void nc_i_writesub_2d(int fid, int varid, size_t * start,
                             size_t * count, int **values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  int **nvals = NULL;
  size_t nstart[4];
  unsigned int i, j;
  int status;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 2);
  nvals = i_alloc_2d(count[offset + 1], count[offset]);
  for (j = 0; j < count[offset]; ++j)
    for (i = 0; i < count[offset + 1]; ++i)
      nvals[j][i] = values[start[offset] + j][start[offset + 1] + i];

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;
  /* 
   * Comment out the following call and uncomment the lines below for
   * debugging purposes. Presently checking the return status is bogus
   * as in the case for nan's
   */
  nc_put_vara_int(fid, varid, nstart, count, nvals[0]);
  /*
    status = nc_put_vara_double(fid, varid, nstart, count, nvals[0]);
    if (status != NC_NOERR)
    hd_quit("dumpfile:nc_d_writesub_2d: %s\n", nc_strerror(status));
  */
  i_free_2d(nvals);
}

/* Write out a double 2d subsection.
 */
void nc_d_writesub_2d(int fid, int varid, size_t * start,
                             size_t * count, double **values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  double **nvals = NULL;
  size_t nstart[4];
  unsigned int i, j;
  int status;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 2);
  nvals = d_alloc_2d(count[offset + 1], count[offset]);
  for (j = 0; j < count[offset]; ++j)
    for (i = 0; i < count[offset + 1]; ++i)
      nvals[j][i] = values[start[offset] + j][start[offset + 1] + i];

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;
  /* 
   * Comment out the following call and uncomment the lines below for
   * debugging purposes. Presently checking the return status is bogus
   * as in the case for nan's
   */
  nc_put_vara_double(fid, varid, nstart, count, nvals[0]);
  /*
    status = nc_put_vara_double(fid, varid, nstart, count, nvals[0]);
    if (status != NC_NOERR)
    hd_quit("dumpfile:nc_d_writesub_2d: %s\n", nc_strerror(status));
  */
  d_free_2d(nvals);
}


/* Write out a double 3d subsection.
 */
void nc_d_writesub_3d(int fid, int varid, size_t * start,
                             size_t * count, double ***values)
{
  /*UR-FIX for icc */
  int nd = 0;
  unsigned int offset = 0;      /* Has a record dimension = 1, else 0 */
  double ***nvals = NULL;
  size_t nstart[4];
  unsigned int i, j, k;
  int status;

  nc_inq_varndims(fid, varid, &nd);
  offset = (nd > 3);
  nvals = d_alloc_3d(count[offset + 2], count[offset + 1], count[offset]);
  for (k = 0; k < count[offset]; ++k)
    for (j = 0; j < count[offset + 1]; ++j)
      for (i = 0; i < count[offset + 2]; ++i)
        nvals[k][j][i] =
          values[start[offset] + k][start[offset + 1] +
                                    j][start[offset + 2] + i];

  for (nstart[0] = start[0], i = offset; i < nd; ++i)
    nstart[i] = 0;

  /* 
   * Comment out the following call and uncomment the lines below for
   * debugging purposes. Presently checking the return status is bogus
   * as in the case for nan's
   */
  nc_put_vara_double(fid, varid, nstart, count, nvals[0][0]);
  /*
    status = nc_put_vara_double(fid, varid, nstart, count, nvals[0][0]);
    if (status != NC_NOERR)
    hd_quit("dumpfile:nc_d_writesub_3d: %s\n", nc_strerror(status));
  */
  d_free_3d(nvals);
}


static void vector_component(dump_data_t *dumpdata, double **vi,
                             double **vj, int ilower, int jlower, int nce1,
                             int nce2, int mode, double **nv)
{
  int i;
  int j;
  int iupper = ilower + nce1 - 1;
  int jupper = jlower + nce2 - 1;

  for (j = jlower; j <= jupper; j++) {
    for (i = ilower; i <= iupper; i++) {
      int fl = i;
      int fr = i <= dumpdata->nce1 ? i + 1 : i;
      int fb = j;
      int ff = j <= dumpdata->nce2 ? j + 1 : j;
      double u1val = (vi[j][fl] + vi[j][fr]) / 2.0;
      double u2val = (vj[fb][i] + vj[ff][i]) / 2.0;
      double sinth =
        (sin(dumpdata->thetau1[j][fl]) + sin(dumpdata->thetau1[j][fr]) +
         sin(dumpdata->thetau2[fb][i]) +
         sin(dumpdata->thetau2[ff][i])) / 4.0;
      double costh =
        (cos(dumpdata->thetau1[j][fl]) + cos(dumpdata->thetau1[j][fr]) +
         cos(dumpdata->thetau2[fb][i]) +
         cos(dumpdata->thetau2[ff][i])) / 4.0;
      double x = u1val * costh - u2val * sinth;
      double y = u1val * sinth + u2val * costh;

      switch (mode) {
      case VM_EAST:
        nv[j][i] = x;
        break;

      case VM_NORTH:
        nv[j][i] = y;
        break;

      case VM_MAG:
        nv[j][i] = sqrt(x * x + y * y);
        break;

      case VM_DIRN:
        nv[j][i] = fmod(atan2(x, y) * 180.0 / M_PI + 180.0, 360.0);
        break;
      }
    }
  }
}

static int decode_args(char *oprm, char *args[], int maxn)
{
  char *p = strdup(oprm);
  int len = strlen(p);
  int i;

  for (i = 0; i < len; ++i) {
    if (p[i] == '(' || p[i] == ')' || p[i] == ',')
      p[i] = ' ';
  }

  return parseline(p, args, maxn);
}

static int contains_string(char *p, char *tag)
{
  char *args[256];
  int len = strlen(p);
  int i, n;

  for (i = 0; i < len; ++i) {
    if (p[i] == '_' || p[i] == '.')
      p[i] = ' ';
  }

  i = parseline(p, args, 256);
  for (n = 0; n < i; n++)
    if (strcmp(args[n], tag) == 0)
      return 1;
  return 0;
}


/* Prototypes for routines below */
static void df_parray_writegeom(dump_data_t *dumpdata, dump_file_t *df);
static void df_parray_init_data(dump_data_t *dumpdata, dump_file_t *df,
                                int fid);
static void df_parray_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int n, double *values);
static void df_parray_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int n, double *values);
static void df_parray_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz);
static void df_parray_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz);
static void df_parray_writevec_cen_2d(dump_data_t *dumpdata, dump_file_t *df,
				      int n, double *values_i,
                                  double *values_j);
static void df_parray_writevec_cen_3d(dump_data_t *dumpdata, dump_file_t *df,
				      int vid, double *values_i,
				      double *values_j, int klower, int nz);
static double get_vector_component(dump_data_t *dumpdata, double vi,
                                   double vj, int c, int mode);


/* Structure to describe each scalar time dep variable */
typedef struct {
  int ndims;                    /* Number of spatial dimensions */
  char name[MAXSTRLEN];         /* df_variable_t name */
  char units[MAXSTRLEN];        /* df_variable_t units */
  char long_name[MAXSTRLEN];    /* Descriptive name */
  int vector_mode;              /* NONE, east, north, mag, dirn */
  void **data[2];               /* Pointer to values (scalar or vector) */
  nc_type fptype;               /* Data type. */
  double valid_range[2];        /* df_variable_t valid range */
  double scale;                 /* Scale factor */
  double offset;                /* Scale factor */
  int xylocation[2];
  int zlocation[2];
  int sediment;                 /* Flag; = 1 if variable is tracer */
} df_parray_var_old_t;

typedef struct {
  void **v;                     /* Pointer to values */
  int ndims;                    /* Number of spatial dimensions */
  char name[MAXSTRLEN];         /* df_variable_t name */
  char units[MAXSTRLEN];        /* df_variable_t units */
  char long_name[MAXSTRLEN];    /* Descriptive name */
  int vector_mode;              /* NONE, east, north, mag, dirn */
  nc_type fptype;               /* Data type. */
  double valid_range[2];        /* df_variable_t valid range */
  double scale;                 /* Scale factor */
  double offset;                /* Scale factor */
  int xylocation[2];
  int zlocation[2];
  int sediment;                 /* Sediment flag */
  void **data[2];               /* Pointer to values (scalar or vector) */

} df_parray_var_t;

/* Structure to describe each dump file. Also used for memory output */
typedef struct {
  int type;                     /* DF_PARRAY or DF_MEMORY */
  int fid;                      /* Netcdf file id */
  int nextrec;                  /* Netcdf record number */
  df_parray_var_t *vars;        /* List of point array dump file variables */
  df_mempack_t *data;           /* Output data structure for memory output */

  GRID_SPECS **gs;              /* Grid spec for parray interpolation */
  delaunay **d;                 /* Delaunay data structure for centre interpolation */
  int *cells;                   /* Cell locations used in the interpolation */
  int *ids;                     /* Delaunay indices used in the interpolation */
  int *c2k;                     /* Interpolation cells to layer map */
  int ncells;                   /* Number of cells used in the interpolation */
  int nz;                       /* Surface layer number */
  int ncs;                      /* Cell centres in surface layer */
  int nc;                       /* Total cell centres */

  GRID_SPECS **ge;              /* Grid spec for parray interpolation */
  delaunay **de;                /* Delaunay data structure for edge interpolation */
  int *eells;                   /* Edge locations used in the interpolation */
  int *eds;                     /* Delaunay indices used in the interpolation */
  int *e2k;                     /* Interpolation edge to layer map */
  int neells;                   /* Number of edgess used in the interpolation */
  int nes;                      /* Edges in surface layer */
  int ec;                       /* Total edge centres */

} df_parray_data_t;

/* Cell corners */
static double xcorner[4] = { 0, 1, 1, 0 };
static double ycorner[4] = { 0, 0, 1, 1 };

/*UR declare here to use in write_geom */
static double get_var_value_2d(dump_file_t *df, df_parray_var_t *var, int id, int k);
static double get_var_value(dump_file_t *df, df_parray_var_t *var, int id, int i, int k);

static double get_missing(df_parray_var_t *var);

INLINE void set_data(df_parray_var_t *var, void **data)
{
  var->vector_mode = VM_NONE;
  var->data[0] = data;
  var->data[1] = data;
}

INLINE void set_data_vector(df_parray_var_t *var, int vm, void *data_i,
                            void *data_j)
{
  var->vector_mode = vm;
  var->data[0] = data_i;
  var->data[1] = data_j;
}

INLINE double *get_data(df_parray_var_t *var)
{
  void *p = var->v;
  return (*(double **)p);
}

INLINE double **get_data_2d(df_parray_var_t *var)
{
  return (*(double ***)var->data[0]);
}

INLINE double *get_data_vector(df_parray_var_t *var, int comp)
{
  void *p = var->data[comp];
  return (*(double **)p);
}

INLINE double **get_data_vector_2d(df_parray_var_t *var, int comp)
{
  return (*(double ***)var->data[comp]);
}

INLINE double ***get_data_3d(df_parray_var_t *var)
{
  return (*(double ****)var->data[0]);
}

INLINE double ***get_data_vector_3d(df_parray_var_t *var, int comp)
{
  return (*(double ****)var->data[comp]);
}

INLINE void set_xyloc(df_parray_var_t *var, int xyloc)
{
  var->xylocation[0] = xyloc;
  var->xylocation[1] = xyloc;
}

INLINE void set_xyloc_vector(df_parray_var_t *var, int xyloc_i,
                             int xyloc_j)
{
  var->xylocation[0] = xyloc_i;
  var->xylocation[1] = xyloc_j;
}

INLINE int get_xyloc(df_parray_var_t *var)
{
  return var->xylocation[0];
}

INLINE int get_xyloc_vector(df_parray_var_t *var, int comp)
{
  return var->xylocation[comp];
}

INLINE void set_zloc(df_parray_var_t *var, int zloc)
{
  var->zlocation[0] = zloc;
  var->zlocation[1] = zloc;
}

INLINE void set_zloc_vector(df_parray_var_t *var, int zloc_i, int zloc_j)
{
  var->zlocation[0] = zloc_i;
  var->zlocation[1] = zloc_j;
}

INLINE int get_zloc(df_parray_var_t *var)
{
  return var->zlocation[0];
}

INLINE int get_zloc_vector(df_parray_var_t *var, int comp)
{
  return var->zlocation[comp];
}


void df_parray_reset(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  df_parray_data_t *data = df->private_data;
  /* Rewind to the first next time after time */
  while (t <= (df->tout - df->tinc))
    df->tout -= df->tinc;
  data->nextrec = find_next_restart_record(df, data->fid, "time", df->tout);
}


void *df_parray_create(dump_data_t *dumpdata, dump_file_t *df)
{
  int nid, kid, vid, npid, tid, n;
  int cdfid;
  int dims[10];
  df_parray_data_t *data = NULL;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double missing_value = PARRAY_MISSING_VALUE;
  double vr[2];
  int kid_sed;
  FILE *fp;
  int nc_mode;

  if (df->npoints <= 0)
    hd_quit
      ("df_parray_create: No valid points have been specified for output.");

  /* If can output can be appended to, then re-open the file, and resync
   * the to next record. */
  fp = fopen(df->name, "rb");
  if (fp != NULL) {
    fclose(fp);
    if (df->append) {
      ncw_open(df->name, NC_WRITE, &cdfid);
      df_parray_init_data(dumpdata, df, cdfid);
      data = (df_parray_data_t *)df->private_data;
      data->nextrec = find_next_restart_record(df, cdfid, "time", dumpdata->t);
      df->finished = 0;
      return data;
    }
  }

  // Get the nc_mode flag
  nc_mode = get_nc_mode(df);

  /* create the netCDF file */
  if (ncw_create(df->name, nc_mode, &cdfid) != NC_NOERR)
    hd_quit("dumpfile_create: Couldn't create output dump file %s\n",
            df->name);

  /* define dimensions */
  nc_def_dim(cdfid, "nrecords", NC_UNLIMITED, &nid);

  nc_def_dim(cdfid, "npoints", df->npoints, &npid);
  nc_def_dim(cdfid, "k", df->nz, &kid);

  if (df->nz_sed > 0)
    nc_def_dim(cdfid, "k_sed", df->nz_sed, &kid_sed);

  /* time independent variables */
  dims[0] = kid;
  nc_def_var(cdfid, "z", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", dumpdata->lenunit);
  write_text_att(cdfid, vid, "long_name", "Z coordinate");
  write_text_att(cdfid, vid, "coordinate_type", "Z");

  dims[0] = npid;
  nc_def_var(cdfid, "x", NC_DOUBLE, 1, dims, &vid);
  if (is_geog) {
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");
  } else {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "X");
    write_text_att(cdfid, vid, "coordinate_type", "X");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);


  nc_def_var(cdfid, "y", NC_DOUBLE, 1, dims, &vid);
  if (is_geog) {
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
  } else {
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Y");
    write_text_att(cdfid, vid, "coordinate_type", "Y");
  }
  if (has_proj)
    write_text_att(cdfid, vid, "projection", projection);

  if (has_proj) {
    nc_def_var(cdfid, "longitude", NC_DOUBLE, 1, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_east");
    write_text_att(cdfid, vid, "long_name", "Longitude");
    write_text_att(cdfid, vid, "coordinate_type", "longitude");

    nc_def_var(cdfid, "latitude", NC_DOUBLE, 1, dims, &vid);
    write_text_att(cdfid, vid, "units", "degrees_north");
    write_text_att(cdfid, vid, "long_name", "Latitude");
    write_text_att(cdfid, vid, "coordinate_type", "latitude");
  }
  /*UR-ADDED botz */
  nc_def_var(cdfid, "botz", NC_DOUBLE, 1, dims, &vid);
  write_text_att(cdfid, vid, "units", "metres");
  write_text_att(cdfid, vid, "long_name", "Depth");
  write_text_att(cdfid, vid, "positive", "down");

  /* Define label variables */
  if (df->label) {
    /* This is the max length we'll use in this file */
    int maxstrlen_id;
    nc_def_dim(cdfid, "MAXSTRLEN", 100, &maxstrlen_id);
    dims[0] = npid;
    dims[1] = maxstrlen_id;
    nc_def_var(cdfid, "label", NC_CHAR, 2, dims, &vid);
    write_text_att(cdfid, vid, "long_name", "Point label");
  }

  /* time dependent variables */
  dims[0] = nid;
  nc_def_var(cdfid, "time", NC_DOUBLE, 1, dims, &tid);
  write_text_att(cdfid, tid, "units", df->tunit);
  write_text_att(cdfid, tid, "long_name", "Time");
  write_text_att(cdfid, tid, "coordinate_type", "time");

  if (df->nz_sed > 0) {
    dims[0] = nid;
    dims[1] = kid_sed;
    dims[2] = npid;
    nc_def_var(cdfid, "z_sed", NC_DOUBLE, 3, dims, &vid);
    write_text_att(cdfid, vid, "units", dumpdata->lenunit);
    write_text_att(cdfid, vid, "long_name", "Z coordinate, sediment");
    write_text_att(cdfid, vid, "coordinate_type", "Z");
  }

  df_parray_init_data(dumpdata, df, cdfid);
  data = (df_parray_data_t *)df->private_data;

  for (n = 0; n < df->nvars; n++) {
    if (!data->vars[n].sediment && data->vars[n].ndims == 2) {
      dims[1] = npid;

      ncw_def_var2(df->name,cdfid, data->vars[n].name, data->vars[n].fptype, 2, dims,&vid, df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, latitude, longitude, y, x");
      else
        write_text_att(cdfid, vid, "coordinates", "time, y, x");
    } else if (!data->vars[n].sediment && data->vars[n].ndims == 3) {
      dims[1] = kid;
      dims[2] = npid;
      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 3, dims, &vid, df->compress);

      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, z, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, z, latitude, longitude, y, x");
      else
        write_text_att(cdfid, vid, "coordinates", "time, z, y, x");
    }

    else if ((df->nz_sed > 0) && data->vars[n].sediment &&
             data->vars[n].ndims == 3) {
      dims[1] = kid_sed;
      dims[2] = npid;
      ncw_def_var2(df->name, cdfid, data->vars[n].name, data->vars[n].fptype, 3, dims, &vid, df->compress);
      if (is_geog)
        write_text_att(cdfid, vid, "coordinates",
                       "time, z_sed, latitude, longitude");
      else if (has_proj)
        write_text_att(cdfid, vid, "coordinates",
                       "time, z_sed, latitude, longitude, y, x");
      else
        write_text_att(cdfid, vid, "coordinates", "time, z_sed, y, x");
    }

    else {
      hd_quit
        ("dumpfile_create: Unable to create variable, incorrect number of dimensions '%d'",
         data->vars[n].ndims);
    }

    if (data->vars[n].fptype == NC_SHORT) {
      nc_put_att_double(cdfid, vid, "scale_factor", NC_DOUBLE,
                        1, &data->vars[n].scale);
      nc_put_att_double(cdfid, vid, "add_offset", NC_DOUBLE,
                        1, &data->vars[n].offset);
    }

    write_text_att(cdfid, vid, "units", data->vars[n].units);
    if (!data->vars[n].sediment)
      write_text_att(cdfid, vid, "long_name", data->vars[n].long_name);
    else {
      char long_name[MAXSTRLEN];
      strcpy(long_name, data->vars[n].long_name);
      strcat(long_name, " in Sediment");
      write_text_att(cdfid, vid, "long_name", long_name);
    }

    vr[0] =
      (data->vars[n].valid_range[0] -
       data->vars[n].offset) / data->vars[n].scale;
    vr[1] =
      (data->vars[n].valid_range[1] -
       data->vars[n].offset) / data->vars[n].scale;
    nc_put_att_double(cdfid, vid, "valid_range", data->vars[n].fptype, 2,
                      vr);
    if (data->vars[n].fptype == NC_SHORT)
      missing_value = PARRAY_SHORT_MISSING_VALUE;
    else
      missing_value = PARRAY_MISSING_VALUE;
    nc_put_att_double(cdfid, vid, "missing_value", data->vars[n].fptype,
                      1, &missing_value);
  }

  /* global attributes */
  write_text_att(cdfid, NC_GLOBAL, "title", codeheader);
  write_text_att(cdfid, NC_GLOBAL, "paramhead", parameterheader);
  write_text_att(cdfid, NC_GLOBAL, "paramfile", dumpdata->prmname);
  write_text_att(cdfid, NC_GLOBAL, "ems_version", version);
  write_text_att(cdfid, NC_GLOBAL, "Conventions", "CMR/Timeseries");
  if (dumpdata->runno >= 0)
    nc_put_att_double(cdfid, NC_GLOBAL, "Run_ID", NC_DOUBLE, 1, &dumpdata->runno);
  if (strlen(dumpdata->rev))
    write_text_att(cdfid, NC_GLOBAL, "Parameter_File_Revision", dumpdata->rev);
 
  nc_enddef(cdfid);

  df_parray_writegeom(dumpdata, df);

  df->finished = 0;

  return data;
}

void df_parray_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  size_t start[4];
  size_t count[4];
  double newt = t;
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  int fid = data->fid;

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = 1L;
  tm_change_time_units(dumpdata->timeunit, df->tunit, &newt,
                         1);
  nc_put_vara_double(fid, ncw_var_id(fid, "time"), start, count, &newt);
  count[0] = 1L;

  if (df->nz_sed != 0) {
    count[1] = df->nz_sed;
    count[2] = df->npoints;

/*  FIX
  nc_put_vara_double(fid,ncw_var_id(fid,"z_sed"),start,count,dumpdata->cellz_sed[0][0]);
*/
  }

  /* Loop over each variable */
  for (n = 0; n < df->nvars; n++) {
    df_parray_var_t *var = &data->vars[n];

    if (var->ndims == 2) {
      if (var->vector_mode == VM_NONE)
        df_parray_writesub_2d(dumpdata, df, n, get_data(var));
      else {
        df_parray_writevec_2d(dumpdata, df, n, get_data(var));
	/*df_parray_writevec_2d(dumpdata, df, n,
                              get_data_vector(var, 0),
                              get_data_vector(var, 1));*/
      }

    } else if (var->ndims == 3) {
      int klower = 0;
      int nz = 0;

      if (var->sediment) {
        nz = df->nz_sed;
        quit("Sed variable output not support for parray file format.");
      } else

      {
        klower = df->klower;
        nz = df->nz;
      }
      if (var->vector_mode & VM_NONE)
        df_parray_writesub_3d(dumpdata, df, n, get_data(var),
                              klower, nz);
      else {
        df_parray_writevec_3d(dumpdata, df, n, get_data(var),
                              klower, nz);
        /*df_parray_writevec_3d(dumpdata, df, n,
                              get_data_vector(var, 0),
                              get_data_vector(var, 1), klower, nz);*/
      }
    }
  }

  data->nextrec++;

  ncw_sync(df->name,fid);
}

void df_parray_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;

  ncw_close(df->name,data->fid);
  data->fid = -1;
  df->finished = 1;
}


static void df_parray_writegeom(dump_data_t *dumpdata, dump_file_t *df)
{
  size_t start[2];
  size_t count[2];
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  int fid = data->fid;
  int i, c, id, k;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);

  double *v = d_alloc_1d(df->npoints);
  double fi = 0.0;
  double fj = 0.0;

  start[0] = df->klower;
  start[1] = 0;
  count[0] = df->nz;
  count[1] = 0;
  nc_d_writesub_1d(fid, ncw_var_id(fid, "z"), start, count,
                   dumpdata->cellz);
  start[0] = 0;
  count[0] = df->npoints;
  nc_d_writesub_1d(fid, ncw_var_id(fid, "x"), start, count, df->x);
  nc_d_writesub_1d(fid, ncw_var_id(fid, "y"), start, count, df->y);

  /* Write lat/long if native, or convert from projection */
  if (is_geog) {
    nc_d_writesub_1d(fid, ncw_var_id(fid, "longitude"), start, count,
                     df->x);
    nc_d_writesub_1d(fid, ncw_var_id(fid, "latitude"), start, count,
                     df->y);

  } else if (has_proj) {
    char *args[256];
    int nargs = parseline(strdup(projection), args, 256);
    map_proj_t *mp = mp_init(nargs, args);
    double *lon = d_alloc_1d(df->npoints);
    double *lat = d_alloc_1d(df->npoints);

    for (i = 0; i < df->npoints; i++)
      mp_inverse(mp, df->x[i], df->y[i], &lat[i], &lon[i]);
    nc_d_writesub_1d(fid, ncw_var_id(fid, "longitude"), start, count, lon);
    nc_d_writesub_1d(fid, ncw_var_id(fid, "latitude"), start, count, lat);
    d_free_1d(lat);
    d_free_1d(lon);
  }

  /* Write labels */
  if (df->label)
    for (i = 0; i < df->npoints; ++i) {
      start[0] = i; /* point i */
      start[1] = 0; /* beginning of the string */
      count[0] = 1; /* write one point */
      count[1] = strlen(df->label[i]) + 1;
      nc_put_vara_text(fid, ncw_var_id(fid, "label"), start, count, df->label[i]);
    }


  /*UR-ADDED botz - write this here and not as a var */
  if(params->parray_inter_botz) {
    delaunay **d = data->d;
    /* Fill the Delaunay structure with values */
    for (i = 0; i < data->ncells; i++) {
      c = data->cells[i];
      if (c == geom->m2d[c]) {
	id = data->ids[i];
	k = data->c2k[i];
	d[k]->points[id].v[0] = geom->botz[c];
      }
    }

    /*UR interpolate in the same fashion as other variables */
    for (i = 0; i < df->npoints; ++i) {
      id = data->ids[i];
      v[i] = get_var_value(df, NULL, id, i, geom->nz-1);
    }
  } else {
    /*UR don't interpolate, just output the value of the column the location falls into*/
    for (i = 0; i < df->npoints; ++i) {
      int c, ic, jc;
      /*if (grid_xytofij(master->xyij_tree, df->x[i], df->y[i], &fi, &fj) >= 0) {
	v[i] = dumpdata->botz[(int)fj][(int)fi];*/
      if ((c = hd_grid_xytoij(master, df->x[i], df->y[i], &ic, &jc))) {
	v[i] = master->geom->botz[c];
      } else
        v[i] = 0.0;
    }
    start[0] = 0;
    start[1] = 0;
    count[0] = df->npoints;
    count[1] = 0;
  }
  nc_put_vara_double(fid, ncw_var_id(fid, "botz"), start, count, v);
  d_free_1d(v);
  /*UR end */

  ncw_sync(df->name,fid);
}


int df_parray_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
                          char *name, df_parray_var_t *var)
{
  int found = 1;
  int n = 0;
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;

  var->v = NULL;
  var->ndims = 2;
  strcpy(var->name, name);
  strcpy(var->units, "");
  strcpy(var->long_name, "");
  var->fptype = nc_get_default_type(df->bpv);
  var->valid_range[0] = 0.0;
  var->valid_range[1] = 0.0;
  var->scale = 1.0;
  var->offset = 0.0;
  var->sediment = 0;
  set_data(var, NULL);
  set_xyloc(var, CL_CENTRE);
  set_zloc(var, CL_NONE);

  if (strcmp(name, "u1av") == 0) {
    set_data_vector(var, VM_NOR, (void *)&master->u1av,
                    (void *)&master->u2av);
    var->v = (void *)&master->u1av;
    set_xyloc_vector(var, CL_SP2|CL_EDGE, CL_SP2|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
    var->vector_mode = VM_NOR;
  }

  else if (strcmp(name, "u2av") == 0) {
    set_data_vector(var, VM_TAN, (void *)&master->u1av,
                    (void *)&master->u2av);
    var->v = (void *)&master->u2av;
    set_xyloc_vector(var, CL_SP2|CL_EDGE, CL_SP2|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
    var->vector_mode = VM_TAN;
  }

  else if (strcmp(name, "uav") == 0) {
    var->v = (void *)&master->uav;
    var->ndims = 2;
    /* Cell centered east velocity is saved in uav, hence does not   */
    /* need to be reconstructed from edge components.                */
    var->vector_mode = VM_NONE;
  }

  else if (strcmp(name, "vav") == 0) {
    var->v = (void *)&master->vav;
    var->ndims = 2;
    /* Cell centered north velocity is saved in vav, hence does not  */
    /* need to be reconstructed from edge components.                */
    var->vector_mode = VM_NONE;
  }

  else if (strcmp(name, "wtop") == 0) {
    var->v = (void *)&master->wtop;
    set_xyloc(var, CL_SP2|CL_FACE);
  }

  else if (strcmp(name, "topz") == 0) {
    var->v = (void *)&master->topz;
    set_xyloc(var, CL_SP2|CL_FACE);
  }

  else if (strcmp(name, "eta") == 0) {
    var->v = (void *)&master->eta;
    set_xyloc(var, CL_SP2|CL_FACE);
    strcpy(var->units, dumpdata->lenunit);
  }

  else if (strcmp(name, "wind1") == 0) {
    hd_quit("Variable 'wind1' cannot be dumped to point array files. Remove from list %s\n", df->name);
    set_data_vector(var, VM_NOR, (void *)&master->wind1,
                    (void *)&master->wind2);
    set_xyloc_vector(var, CL_SP2|CL_EDGE, CL_SP2|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
    var->vector_mode = VM_NOR;
  }

  else if (strcmp(name, "patm") == 0) {
    var->v = (void *)&master->patm;
    set_xyloc(var, CL_SP2|CL_FACE);
  }

  else if (strcmp(name, "u") == 0) {
    var->v = (void *)&master->u;
    set_xyloc(var, CL_SP3|CL_EDGE);
    var->ndims = 3;
    /* Cell centered east velocity is saved in u, hence does not     */
    /* need to be reconstructed from edge components.                */
    var->vector_mode = VM_NONE;
    strcpy(var->units,"ms-1");
  }

  else if (strcmp(name, "v") == 0) {
    var->v = (void *)&master->v;
    set_xyloc(var, CL_SP3|CL_EDGE);
    var->ndims = 3;
    /* Cell centered north velocity is saved in v, hence does not    */
    /* need to be reconstructed from edge components.                */
    var->vector_mode = VM_NONE;
    strcpy(var->units,"ms-1");
  }

  else if (strcmp(name, "u1") == 0) {
    var->ndims = 3;
    var->v = (void *)&master->u1;
    set_data_vector(var, VM_EAST, (void *)&master->u1,
                    (void *)&master->u2);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
    var->vector_mode = VM_NOR;
  }

  else if (strcmp(name, "u2") == 0) {
    var->ndims = 3;
    var->v = (void *)&master->u2;
    set_data_vector(var, VM_NORTH, (void *)&master->u1,
                    (void *)&master->u2);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
    var->vector_mode = VM_TAN;
  }

  else if (strcmp(name, "w") == 0) {
    var->v = (void *)&master->w;
    set_xyloc(var, CL_SP3|CL_FACE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else if (strcmp(name, "u1mean") == 0) {
    var->ndims = 3;
    set_data_vector(var, VM_EAST, (void *)&master->u1m,
                    (void *)&master->u2m);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
  }

  else if (strcmp(name, "u2mean") == 0) {
    var->ndims = 3;
    set_data_vector(var, VM_NORTH, (void *)&master->u1m,
                    (void *)&master->u2m);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
  }

  else if (strcmp(name, "u1vmean") == 0) {
    var->ndims = 3;
    set_data_vector(var, VM_EAST, (void *)&master->u1vm,
                    (void *)&master->u2vm);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
  }

  else if (strcmp(name, "u2vmean") == 0) {
    var->ndims = 3;
    set_data_vector(var, VM_NORTH, (void *)&master->u1vm,
                    (void *)&master->u2vm);
    set_xyloc_vector(var, CL_SP3|CL_EDGE, CL_SP3|CL_EDGE);
    set_zloc_vector(var, CL_CENTRE, CL_CENTRE);
  }

  else if (strcmp(name, "dens") == 0) {
    var->v = (void *)&master->dens;
    set_xyloc(var, CL_SP3|CL_FACE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else if (strcmp(name, "dens_0") == 0) {
    var->v = (void *)&master->dens_0;
    set_xyloc(var, CL_SP3|CL_FACE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else if (strcmp(name, "Kz") == 0) {
    var->v = (void *)&master->Kz;
    set_xyloc(var, CL_SP3|CL_FACE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else if (strcmp(name, "Vz") == 0) {
    var->v = (void *)&master->Vz;
    set_xyloc(var, CL_SP3|CL_FACE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else if (strcmp(name, "Cd") == 0) {
    var->v = (void *)&master->Cd;
    set_xyloc(var, CL_SP2|CL_FACE);
  }

  else if (strcmp(name, "u1vh") == 0) {
    hd_quit("Variable 'u1vh' cannot be dumped to point array files. Remove from list %s\n", df->name);
    var->v = (void *)&master->u1vh;
    set_xyloc(var, CL_SP3|CL_EDGE);
    set_zloc(var, CL_CENTRE);
    var->ndims = 3;
  }

  else
    found = 0;

  if (!found) {
    /* 3D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntr; n++) {
      if (strcmp(dumpdata->trinfo_3d[n].name, name) == 0) {
        var->v = (void **)&master->tr_wc[n];
	set_xyloc(var, CL_SP3|CL_FACE);
	set_zloc(var, CL_CENTRE);
	var->ndims = 3;
        strcpy(var->units, dumpdata->trinfo_3d[n].units);
        found = 1;
        break;
      }
    }
    /* 2D Watercolumn tracers */
    for (n = 0; n < dumpdata->ntrS; n++) {
      if (strcmp(dumpdata->trinfo_2d[n].name, name) == 0) {
        var->v = (void **)&master->tr_wcS[n];
	set_xyloc(var, CL_SP2|CL_FACE);
        found = 1;
        break;
      }
    }
  }

  return found;
}

static void df_parray_init_data(dump_data_t *dumpdata, dump_file_t *df,
                                int fid)
{
  int i;
  df_parray_data_t *data = NULL;

  df_parse_vars(dumpdata,df,PARRAY_EXCLUDE_VARS,PARRAY_ALL_VARS);
  
  // Clean up before allocating more memory
  if (df->private_data != NULL) {
    if (((df_parray_data_t *)df->private_data)->vars != NULL)
      free(((df_parray_data_t *)df->private_data)->vars);
    free(df->private_data);
  }
  
  data = (df_parray_data_t *)malloc(sizeof(df_parray_data_t));
  memset(data, 0, sizeof(df_parray_data_t));
  df->private_data = data;
  data->type = DF_PARRAY;

  /* Allocate memory for variable information */
  if ((data->vars =
       (df_parray_var_t *)malloc(df->nvars * sizeof(df_parray_var_t))) ==
      NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_parray_var_t *var = &data->vars[i];
    if (df_parray_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_parray_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  data->fid = fid;

  /* Allocate and initialize the grid_spec structures for            */
  /* interpolation. Note: two of these exist for cell centered       */
  /* variables (which can be iunterpolated onto a geographic         */
  /* location using the i_rule) or edge variables (which are dumped  */
  /* as the nearest neighbour edge to the geographic location).      */
  /* The var->vector_mode distinguished between centred and edge     */
  /* variables.                                                      */
  parray_grid_init(dumpdata, df);
  parray_gride_init(dumpdata, df);

  /* Set next record to zero */
  data->nextrec = 0;

  if (df->diagn && DEBUG("dump"))
    dlog("dump",
         "Diagnostic tracers initialisation has been synchronised with \"%s\" output; period = %f s.\n",
         df->name, df->tinc);
}

static double get_var_value(dump_file_t *df, df_parray_var_t *var, int id, int i, int k)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  GRID_SPECS **gs = (var->vector_mode & (VM_NOR|VM_TAN)) ? data->ge : data->gs;
  double x, y;
  double v = 0.0;

  if (df->osl & (L_BAYLIN|L_BILIN)) {
    x = (int)id;
    y = (int)id;
  } else {
    x = df->x[i];
    y = df->y[i];
  }
  v = grid_interp_on_point(gs[k], x, y);

  if (df->osl & (L_BAYLIN|L_BILIN)) {
    x = df->x[i];
    y = df->y[i];
  }
  if (var != NULL && var->fptype == NC_SHORT)
    return (v - var->offset) / var->scale;

  return v;
}

static double get_var_value_2d(dump_file_t *df, df_parray_var_t *var, int id, int k)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  GRID_SPECS **gs;
  double x, y;
  double v = 0.0;

  if (var == NULL)
    gs = data->gs;
  else
    gs = (var->vector_mode & (VM_NOR|VM_TAN)) ? data->ge : data->gs;
  
  if (df->osl & (L_BAYLIN|L_BILIN)) {
    x = (int)id;
    y = (int)id;
  } else {
    x = data->d[k]->points[id].x;
    y = data->d[k]->points[id].y;
  }
  v = grid_interp_on_point(gs[k], x, y);

  if (df->osl & (L_BAYLIN|L_BILIN)) {
    x = data->d[k]->points[id].x;
    y = data->d[k]->points[id].y;
  }

  if(strcmp(df->name,"/home/her127/work/meco/compas/est/tile0/bdry0-1_uv_nor.mpk.nc")==0) {
    int i;
    for (i=0; i<data->d[k]->npoints;i++)
      printf("%f %f\n",data->d[k]->points[i].x,data->d[k]->points[i].y);
  }
    
  if (var != NULL && var->fptype == NC_SHORT)
    return (v - var->offset) / var->scale;

  return v;
}

static double get_var_value_3d(dump_file_t *df, df_parray_var_t *var, int id, int k)
{
  double v;

  v = get_var_value_2d(df, var, id, k);
  return(v);

}

static double get_missing(df_parray_var_t *var)
{
  if (var->fptype == NC_SHORT)
    return PARRAY_SHORT_MISSING_VALUE;
  return NaN;
}

static void df_parray_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d =  data->d;
  GRID_SPECS **gs = data->gs;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, k, c, id;
  int fi, fj;
  double *v = d_alloc_1d(df->npoints);
  size_t start[2];
  size_t count[2];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->ncells; i++) {
    c = data->cells[i];
    if (c == geom->m2d[c]) {
      id = data->ids[i];
      k = data->c2k[i];
      d[k]->vid = 0;
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, c);
      else
	d[k]->points[id].v[0] = values[c];
      d[k]->points[id].z = d[k]->points[id].v[0];
      if (df->osl & (L_BAYLIN|L_BILIN)) {
	d[k]->points[id].z = (double)id;
	d[k]->points[id].v[0] = (double)id;
      }
    }
  }

  /* Rebuild the weights.                                            */
  k = geom->nz-1;
  if (df->osl & (L_SIB|L_NONSIB))
    for (i = 0; i < d[k]->npoints; i++)
      gs[k]->rebuild(gs[k]->interpolator, &d[k]->points[i]);
  else
    gs[k]->rebuild(gs[k]->interpolator, d[k]->points);
  if(df->osl & (L_BILIN|L_BAYLIN)) {
    for (i = 0; i < data->ncells; i++) {
      c = data->cells[i];
      if (c == geom->m2d[c]) {
	id = data->ids[i];
	k = data->c2k[i];
	if (df->filter)
	  d[k]->points[id].v[0] = df_filter(df, values, c);
	else
	  d[k]->points[id].v[0] = values[c];
      }
    }
  }

  for(i = 0; i < df->npoints; i++) {
    v[i] = get_missing(var);
  }

  for (i = 0; i < df->npoints; ++i) {
    id = data->ids[i];
    v[i] = get_var_value(df, var, id, i, k);
  }

  start[0] = data->nextrec;
  start[1] = 0;
  count[0] = 1;
  count[1] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v);

  d_free_1d(v);

}

static void df_parray_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = data->d;
  GRID_SPECS **gs = data->gs;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, k, c, id;
  int fi, fj;
  double **v = d_alloc_2d(df->npoints, nz);
  size_t start[3];
  size_t count[3];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->ncells; i++) {
    c = data->cells[i];
    id = data->ids[i];
    k = data->c2k[i];
    d[k]->vid = 0;
    if (df->filter)
      d[k]->points[id].v[0] = df_filter(df, values, c);
    else
      d[k]->points[id].v[0] = values[c];
    d[k]->points[id].z = d[k]->points[id].v[0];
    if (df->osl & (L_BAYLIN|L_BILIN)) {
      d[k]->points[id].z = (double)id;
      d[k]->points[id].v[0] = (double)id;
    }
  }

  /* Rebuild the weights.                                            */
  for (k = 0; k < geom->nz; k++) {
    if (d[k] == NULL) continue;
    if (df->osl & (L_SIB|L_NONSIB))
      for (i = 0; i < d[k]->npoints; i++)
	data->gs[k]->rebuild(data->gs[k]->interpolator, &d[k]->points[i]);
    else
      data->gs[k]->rebuild(data->gs[k]->interpolator, d[k]->points);
  }
  if(df->osl & (L_BILIN|L_BAYLIN)) {
    for (i = 0; i < data->ncells; i++) {
      c = data->cells[i];
      id = data->ids[i];
      k = data->c2k[i];
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, c);
      else
	d[k]->points[id].v[0] = values[c];
    }
  }

  for(i = 0; i < df->npoints; i++)
    for (k = 0; k < nz; ++k)
      v[k][i] = 1.0;
  /*v[k][i] = get_missing(var);*/
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < nz; ++k){
      if (d[k] == NULL) continue;
      id = data->ids[i];
      v[k][i] = get_var_value(df, var, id, i, k);
    }
  }

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = nz;
  count[2] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v[0]);

  d_free_2d(v);
}

static void df_parray_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = data->de;
  GRID_SPECS **gs = data->ge;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, k, e, id;
  int fi, fj;
  double **v = d_alloc_2d(df->npoints, nz);
  size_t start[3];
  size_t count[3];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->neells; i++) {
    e = data->eells[i];
    id = data->eds[i];
    k = data->e2k[i];
    d[k]->vid = 0;
    if (df->filter)
      d[k]->points[id].v[0] = df_filter(df, values, e);
    else
      d[k]->points[id].v[0] = values[e];
    d[k]->points[id].z = d[k]->points[id].v[0];
  }

  /* Rebuild the weights.                                            */
  for (k = 0; k < geom->nz; k++) {
    if (d[k] == NULL) continue;
    data->gs[k]->rebuild(data->gs[k]->interpolator, d[k]->points);
  }

  for(i = 0; i < df->npoints; i++)
    for (k = 0; k < nz; ++k)
      v[k][i] = 1.0;
  /*v[k][i] = get_missing(var);*/
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < nz; ++k){
      if (d[k] == NULL) continue;
      id = data->eds[i];
      v[k][i] = get_var_value(df, var, id, i, k);
    }
  }

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = nz;
  count[2] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v[0]);

  d_free_2d(v);
}


static void df_parray_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = data->de;
  GRID_SPECS **gs = data->ge;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, k, e, id;
  int fi, fj;
  double *v = d_alloc_1d(df->npoints);
  size_t start[2];
  size_t count[2];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->neells; i++) {
    e = data->eells[i];
    if (e == geom->m2de[e]) {
      id = data->eds[i];
      k = data->e2k[i];
      d[k]->vid = 0;
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, e);
      else
	d[k]->points[id].v[0] = values[e];
      d[k]->points[id].z = d[k]->points[id].v[0];
    }
  }

  /* Rebuild the weights.                                            */
  k = geom->nz-1;
  gs[k]->rebuild(gs[k]->interpolator, d[k]->points);

  for(i = 0; i < df->npoints; i++) {
    v[i] = get_missing(var);
  }
  for (i = 0; i < df->npoints; ++i) {
    id = data->eds[i];
    v[i] = get_var_value(df, var, id, i, k);
  }

  start[0] = data->nextrec;
  start[1] = 0;
  count[0] = 1;
  count[1] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v);

  d_free_1d(v);
}


static void df_parray_writevec_cen_2d(dump_data_t *dumpdata, dump_file_t *df,
				      int vid, double *values_i,
				      double *values_j)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = data->d;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, k, c, id;
  int ee, e, es;
  int fi, fj;
  double nu, nv, a;
  double *v = d_alloc_1d(df->npoints);
  size_t start[2];
  size_t count[2];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->ncells; i++) {
    c = data->cells[i];
    if (c == geom->m2d[c]) {
      id = data->ids[i];
      k = data->c2k[i];
      d[k]->points[id].v[0] = 0.0;
      d[k]->points[id].v[1] = 0.0;
      nu = nv = 0.0;
      for (ee = 1; ee <= geom->npe[c]; ee++) {
	e = geom->c2e[ee][c];
	/* Get the cell centered east and north velocity               */
	a = 0.5 * geom->h1au1[es] * geom->h2au1[es];
	d[k]->points[id].v[0] += a * (values_i[e] * geom->costhu1[es] + values_j[e] * geom->costhu2[es]);
	nu += a;
	d[k]->points[id].v[1] += a * (values_i[e] * geom->sinthu1[es] + values_j[e] * geom->sinthu2[es]);
	nv += a;
      }
      d[k]->points[id].v[0] = (nu) ? d[k]->points[id].v[0] / nu : 0.0;
      d[k]->points[id].v[1] = (nv) ? d[k]->points[id].v[1] / nv : 0.0;
    }
  }

  k = geom->nz-1;
  for(i = 0; i < df->npoints; i++)
    v[i] = get_missing(var);

  for (i = 0; i < df->npoints; ++i) {
    double vi, vj;
    id = data->ids[i];
    c = data->cells[i];
    if (c == geom->m2d[c]) {
      d[k]->vid = 0;
      vi = get_var_value(df, var, id, i, k);
      d[k]->vid = 1;
      vj = get_var_value(df, var, id, i, k);
      v[i] =
        get_vector_component(dumpdata, vi, vj, c,
                             var->vector_mode);
    }
  }

  start[0] = data->nextrec;
  start[1] = 0;
  count[0] = 1;
  count[1] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v);

  d_free_1d(v);
}

static void df_parray_writevec_cen_3d(dump_data_t *dumpdata, dump_file_t *df,
				      int vid, double *values_i,
				      double *values_j, int klower, int nz)
{
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &data->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = data->d;
  int fid = data->fid;
  int varid = ncw_var_id(fid, df->vars[vid]);
  int i, j, k, c, c2, id;
  int ee, e, es;
  int fi, fj;
  double nu, nv, a;
  double **v = d_alloc_2d(df->npoints, nz);
  size_t start[3];
  size_t count[3];

  /* Fill the Delaunay structure with values */
  for (i = 0; i < data->ncells; i++) {
    c = data->cells[i];
    c2 = geom->m2d[c];
    id = data->ids[i];
    k = data->c2k[i];
    d[k]->points[id].v[0] = 0.0;
    d[k]->points[id].v[1] = 0.0;
    nu = nv = 0.0;
    for (ee = 1; ee <= geom->npe[c2]; ee++) {
      e = geom->c2e[ee][c];
      es = geom->m2de[e];
      /* Get the cell centered east and north velocity               */
      a = 0.5 * geom->h1au1[es] * geom->h2au1[es];
      d[k]->points[id].v[0] += a * (values_i[e] * geom->costhu1[es] + values_j[e] * geom->costhu2[es]);
      nu += a;
      d[k]->points[id].v[1] += a * (values_i[e] * geom->sinthu1[es] + values_j[e] * geom->sinthu2[es]);
      nv += a;
    }
    d[k]->points[id].v[0] = (nu) ? d[k]->points[id].v[0] / nu : 0.0;
    d[k]->points[id].v[1] = (nv) ? d[k]->points[id].v[1] / nv : 0.0;
  }

  for(i = 0; i < df->npoints; i++)
    for (k = 0; k < nz; ++k)
      v[k][i] = get_missing(var);
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < nz; ++k){
      double vi, vj;
      if (d[k] == NULL) continue;
      id = data->ids[i];
      c = data->cells[i];
      d[k]->vid = 0;
      vi = get_var_value(df, var, id, i, klower + k);
      d[k]->vid = 1;
      vj = get_var_value(df, var, id, i, klower + k);
      v[k][i] = get_vector_component(dumpdata, vi, vj, c,
				     var->vector_mode);
    }
  }

  start[0] = data->nextrec;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = nz;
  count[2] = df->npoints;

  nc_put_vara_double(fid, varid, start, count, v[0]);

  d_free_2d(v);
}

static double get_vector_component(dump_data_t *dumpdata, double vi,
                                   double vj, int c, int mode)
{
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  int e, j;
  double sinth = 0.0, costh = 0.0;
  int npe = geom->npe[geom->m2d[c]];

  for (j = 1; j <= npe; j++) {
    e = geom->m2de[geom->c2e[j][c]];
    sinth += geom->sinthu1[e];
    costh += geom->costhu1[e];
  }

  sinth /= (double)npe;
  costh /= (double)npe;
  double x = vi * costh + vj * sinth;
  double y = vj * costh - vi * sinth;
  double nv = 0.0;

  switch (mode) {
  case VM_NOR:
    nv = x;
    break;

  case VM_TAN:
    nv = y;
    break;

  case VM_EAST:
    nv = vi;
    break;

  case VM_NORTH:
    nv = vj;
    break;

  case VM_MAG:
    nv = sqrt(x * x + y * y);
    break;

  case VM_DIRN:
    nv = fmod(atan2(x, y) * 180.0 / M_PI + 180.0, 360.0);
    break;
  }

  return nv;
}

int find_next_restart_record(dump_file_t *df, int cdfid,
                                    char *timevar, double tref) {
  int ndims = 0;
  int nextrec = 0;
  int vid = ncw_var_id(cdfid, timevar);
  if ((nc_inq_varndims(cdfid, vid, &ndims) == NC_NOERR) && (ndims == 1)) {
     int dimid = 0;
     size_t start = 0;
     size_t nt = 0;
     nc_inq_vardimid(cdfid, vid, &dimid);
     nc_inq_dimlen(cdfid, dimid, &nt);

     if (nt > 0) {
        double *t = d_alloc_1d(nt);
        int i;
        nc_get_vara_double(cdfid,vid,&start,&nt,t);
        if (master->timeunit) {
           char timeunits[MAXSTRLEN];
           memset(timeunits,0,MAXSTRLEN);
           nc_get_att_text(cdfid,vid,"units",timeunits);
           tm_change_time_units(timeunits, master->timeunit, t, nt);    
        }

        for (i=0; i<nt; ++i, ++nextrec) {
          if (tref <= t[i])
             break;
        }

        d_free_1d(t);
     }
  } else {
    hd_quit("Dumpfile(%s): Expected variable '%s' to have single dimension.\n",
             df->name, timevar);
  }

  return nextrec;
}


/*-------------------------------------------------------------------*/
/* Creates a GRID_SPEC structure and Delaunay triangulation for      */
/* triangular linear interpolation of output.                        */
/*-------------------------------------------------------------------*/
void parray_grid_init(dump_data_t *dumpdata, dump_file_t *df)
{
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  char *uvrule = df->irule;
  int cc, c, c2, cn, ci, k, kk;
  int i, j, id, fi, fj;
  int *nk, nz = geom->nz;
  int np;
  int *cells;
  point **p;
  double v;
  delaunay **d, *dp;
  int n, ee, e;
  int *mask, *n2i;
  GRID_SPECS **gs;
  int isalloc;

  isalloc = (data->gs && data->d) ? 1 : 0;

  if (!isalloc) {
    /* Get the number of points for the triangulation                */
    np = 0;
    data->ncells = 0;
    data->nz = geom->nz - 1;
    mask = i_alloc_1d(geom->szc);
    memset(mask, 0, geom->szc * sizeof(int));
    for (i = 0; i < df->npoints; ++i) {
      c = hd_grid_xytoij(master, df->x[i], df->y[i], &fi, &fj);
      if (c > 0 && c < geom->szcS) {
	c2 = geom->m2d[c];
	while (c != geom->zm1[c]) {
	  if (!mask[c]) {
	    if (c == geom->m2d[c]) np++;
	    data->ncells ++;
	    mask[c] = 1;
	  }

	  for (j = 1; j <= geom->npe[c2]; j++) {
	    cn = geom->c2c[j][c];
	    if (geom->wgst[cn]) continue;
	    if (!mask[cn]) {
	      if (cn == geom->m2d[cn]) np++;
	      data->ncells ++;
	      mask[cn] = 1;
	    }
	  }

	  data->nz = min(data->nz, geom->s2k[c]);
	  c = geom->zm1[c];
	}
      }
    }
    if (!data->ncells) hd_quit("parray_grid_init: Can't find any cells for file %s Delaunay interpolation.\n", df->name);

    /* Allocate                                                      */
    data->cells = i_alloc_1d(data->ncells);
    data->ids = i_alloc_1d(data->ncells);
    data->c2k = i_alloc_1d(data->ncells);
    n2i = i_alloc_1d(data->ncells);
    p = (point **)alloc_2d(np, nz, sizeof(point));
    nk = i_alloc_1d(nz);
    memset(nk, 0, nz * sizeof(int));

    /* Set the surface cell centres for the triangulation            */
    n = 0;
    memset(mask, 0, geom->szc * sizeof(int));
    for (i = 0; i < df->npoints; ++i) {
      c = hd_grid_xytoij(master, df->x[i], df->y[i], &fi, &fj);
      if (c > 0 && c < geom->szcS) {
	c2 = geom->m2d[c];
	k = geom->s2k[c];
	if (!mask[c]) {
	  n2i[n] = i;
	  data->ids[n] = nk[k];
	  data->c2k[n] = k;
	  data->cells[n++] = c;
	  mask[c] = nk[k] + 1;
	  p[k][nk[k]].x = geom->cellx[c2];
	  p[k][nk[k]++].y = geom->celly[c2];
	  /*
	  if(strcmp(df->name,"/home/her127/work/meco/compas/est/tile0/bdry0-1_uv_nor.mpk.nc")==0)
	    if(c==c2)printf("%f %f\n",p[k][nk[k]-1].x,p[k][nk[k]-1].y);
	  */
	}
      }
    }
    data->ncs = n;
    /* Set the sub-surface cell centres for the triangulation        */
    for (i = 0; i < data->ncs; i++) {
      c = data->cells[i];
      c2 = geom->m2d[c];
      c = geom->zm1[c];
      id = data->ids[i];
      while (c != geom->zm1[c]) {
	k = geom->s2k[c];
	if (!mask[c]) {
	  n2i[n] = i;
	  data->ids[n] = nk[k];
	  data->c2k[n] = k;
	  data->cells[n++] = c;
	  mask[c] = nk[k] + 1;
	  p[k][nk[k]].x = geom->cellx[c2];
	  p[k][nk[k]++].y = geom->celly[c2];
	}
	c = geom->zm1[c];	
      }
    }
    data->nc = n;
    /* Set the cells surrounding the cell centres                    */
    for (i = 0; i < data->ncs; i++) {
      c = data->cells[i];
      c2 = geom->m2d[c];
      while (c != geom->zm1[c]) {
	k = geom->s2k[c];
	for (j = 1; j <= geom->npe[c2]; j++) {
	  cn = geom->c2c[j][c];
	  if (geom->wgst[cn]) continue;
	  if (!mask[cn]) {
	    data->ids[n] = nk[k];
	    data->c2k[n] = k;
	    data->cells[n++] = cn;
	    mask[cn] = nk[k] + 1;
	    /* The geographic location of ghost cells is set to the  */
	    /* cell edge so as to account for multiple ghost cells   */
	    /* in the triangulation.                                 */
	    if (geom->wgst[cn]) {
	      e = geom->m2de[geom->c2e[j][c]];
	      p[k][nk[k]].x = geom->u1x[e];
	      p[k][nk[k]].y = geom->u1y[e];
	    } else {
	      p[k][nk[k]].x = geom->cellx[geom->m2d[cn]];
	      p[k][nk[k]].y = geom->celly[geom->m2d[cn]];
	    }
	    /*
	  if(strcmp(df->name,"/home/her127/work/meco/compas/est/tile0/bdry0-1_uv_nor.mpk.nc")==0)
	    if(c==c2)printf("%f %f\n",p[k][nk[k]].x,p[k][nk[k]].y);
	    */
	    nk[k]++;
	  }
	}
	c = geom->zm1[c];
      }
    }

    /*---------------------------------------------------------------*/
    /* Fill the points array and create the triangulation for every  */
    /* layer. Note this is designed to be used with COMPAS arrays,   */
    /* which start at index 1, hence use [c+1].                      */
    data->d = (delaunay **)calloc(nz, sizeof(delaunay *));
    for (k = 0; k < nz; k++) {
      delaunay *dk;
      data->d[k] = NULL;
      if (nk[k] == 0) continue;
      /*
      if(strcmp(df->name,"/home/her127/work/meco/compas/est/tile0/bdry0-1_uv_nor.mpk.nc")==0)
	if(k==nz-1)for(i=0;i<nk[k];i++)printf("%f %f\n",p[k][i].x,p[k][i].y);
      */
      data->d[k] = delaunay_build(nk[k], p[k], 0, NULL, 0, NULL);
      dk = data->d[k];
      dk->vid = 0;
      dk->ptf = 1;

      /* Free the point triangle arrays                              */
      if (dk->point_triangles != NULL) {
	for (i = 0; i < dk->npoints; ++i)
	  if (dk->point_triangles[i] != NULL)
	    free(dk->point_triangles[i]);
      }
      memset(dk->n_point_triangles, 0, dk->npoints * sizeof(int));

      for (c = 0; c < nk[k]; c++) {
	dk->points[c].v = d_alloc_1d(2);
      }
    }

    /* Set the new Delaunay point triangles for the interpolation    */
    for (i = 0; i < data->ncells; i++) {
      delaunay *dk;
      c = data->cells[i];
      c2 = geom->m2d[c];
      k = geom->s2k[c];
      id = data->ids[i];
      dk = data->d[k];

      /* Cell centres in the triangulation                           */
      if (i < data->nc) {
	dk->n_point_triangles[id] = geom->npe[c2] + 1;
	dk->point_triangles[id] = malloc(dk->n_point_triangles[id] * sizeof(int));
	dk->point_triangles[id][0] = id;
	for (j = 1; j <= geom->npe[c2]; j++) {
	  cn = geom->c2c[j][c];
	  dk->point_triangles[id][j] = mask[cn] - 1;
	}
      } else {
	/* Cells surrounding the cell centres                        */
	dk->n_point_triangles[id] = 1;
	dk->point_triangles[id] = malloc(dk->n_point_triangles[id] * sizeof(int));
	dk->point_triangles[id][0] = -1;
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate and initialize                                       */
    gs = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    for (k = 0; k < nz; k++) {
      
      gs[k] = grid_spec_create();
      if (data->d[k] == NULL) continue;
      
      grid_interp_init_t(gs[k], data->d[k], uvrule, 0);
    }

    /*---------------------------------------------------------------*/
    /* Set the grid_spec structure                                   */
    data->gs = gs;

  }

  /*-----------------------------------------------------------------*/
  /* Fill the points array in the delaunay structure.                */
  for (i = 0; i < data->nc; i++) {
    c = data->cells[i];
    k = geom->s2k[c];
    id = data->ids[i];
    /*
    data->d[k]->points[id].x = df->x[n2i[i]];
    data->d[k]->points[id].y = df->y[n2i[i]];
    */
    data->d[k]->points[id].z = (double)id;
  }

  /*-----------------------------------------------------------------*/
  /* Rebuild the weights. Note: this is cumulative (i.e. weights are */
  /* added to new cells) for each parray point.                      */
  for (k = 0; k < nz; k++) {
    if (data->d[k] == NULL) continue;
    data->gs[k]->rebuild(data->gs[k]->interpolator, data->d[k]->points);
  }
  i_free_1d(mask);
  i_free_1d(n2i);
}

/* END parray_grid_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a GRID_SPEC structure and Delaunay triangulation for      */
/* interpolation of edge output. Edges included in the triangulation */
/* are the nearest edge to a dump location, and all edges            */
/* surrounding vertices at either end of that edge.                  */
/*-------------------------------------------------------------------*/
void parray_gride_init(dump_data_t *dumpdata, dump_file_t *df)
{
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  df_parray_data_t *data = (df_parray_data_t *)df->private_data;
  char *uvrule = "nearest";
  int c, c2, ee, e, e2, en, ei, k, kk;
  int i, j, id, fi, fj, n;
  int *nk, nz = geom->nz;
  int np;
  int *eells, *eloc;
  point **p;
  int vv, v, vs;
  double val;
  delaunay **d;
  int *mask, *n2i;
  GRID_SPECS **ge;
  int isalloc;
  double x, y, dist, dm;

  isalloc = (data->gs && data->de) ? 1 : 0;
  /* Check for edge vectors requiring output                         */
  if (!isalloc) {
    isalloc = 1;
    for (n = 0; n < df->nvars; n++) {
      if (data->vars[n].vector_mode & (VM_NOR|VM_TAN))
	isalloc = 0;
    }
  }
  if (isalloc) return;

    /* Get the number of points for the triangulation                */
    np = 0;
    data->neells = 0;
    data->nz = geom->nz - 1;
    mask = i_alloc_1d(geom->sze);
    memset(mask, 0, geom->sze * sizeof(int));
    eloc = i_alloc_1d(df->npoints);
    for (i = 0; i < df->npoints; ++i) {
      c = hd_grid_xytoij(master, df->x[i], df->y[i], &fi, &fj);
      if (c > 0 && c < geom->szcS) {
	/* Find the edge closest to the dump location                */
	c2 = geom->m2d[c];
	dm = HUGE;
	for (j = 1; j <= geom->npe[c2]; j++) {
	  x = df->x[i] - geom->u1x[geom->c2e[j][c2]];
	  x = df->y[i] - geom->u1y[geom->c2e[j][c2]];
	  dist =sqrt(x * x + y * y);
	  if (dist < dm) {
	    dm = dist;
	    e = geom->c2e[j][c];
	  }
	}
	eloc[i] = e;
	e2 = geom->m2de[e];
	while (e != geom->zm1e[e]) {
	  if (!mask[e]) {
	    if (e == geom->m2de[e]) np++;
	    data->neells ++;
	    mask[e] = 1;
	  }
	  for (j = 0; j < 2; j++) {
	    v = geom->e2v[e][j];
	    vs = geom->m2dv[v];
	    for (vv = 1; vv <= geom->nve[vs]; vv++) {
	      en = geom->v2e[v][vv];
	      if (!mask[en]) {
		if (en == geom->m2de[en]) np++;
		data->neells ++;
		mask[en] = 1;
	      }
	    }
	  }
	  data->nz = min(data->nz, geom->e2k[e]);
	  e = geom->zm1e[e];
	}
      }
    }

    /* Allocate                                                      */
    data->eells = i_alloc_1d(data->neells);
    data->eds = i_alloc_1d(data->neells);
    data->e2k = i_alloc_1d(data->neells);
    n2i = i_alloc_1d(data->neells);
    p = (point **)alloc_2d(np, nz, sizeof(point));
    nk = i_alloc_1d(nz);
    memset(nk, 0, nz * sizeof(int));

    /* Set the surface cell centres for the triangulation            */
    n = 0;
    memset(mask, 0, geom->sze * sizeof(int));
    for (i = 0; i < df->npoints; ++i) {
      e = eloc[i];
      e2 = geom->m2de[e];
      k = geom->e2k[e];
      if (!mask[e]) {
	n2i[n] = i;
	data->eds[n] = nk[k];
	data->e2k[n] = k;
	data->eells[n++] = e;
	mask[e] = nk[k] + 1;
	p[k][nk[k]].x = geom->u1x[e2];
	p[k][nk[k]++].y = geom->u1y[e2];
      }
    }
    data->nes = n;
    /* Set the sub-surface cell centres for the triangulation        */
    for (i = 0; i < data->nes; i++) {
      e = data->eells[i];
      e2 = geom->m2de[e];
      e = geom->zm1e[e];
      while (e != geom->zm1e[e]) {
	k = geom->e2k[e];
	if (!mask[e]) {
	  n2i[n] = i;
	  data->eds[n] = nk[k];
	  data->e2k[n] = k;
	  data->eells[n++] = e;
	  mask[e] = nk[k] + 1;
	  p[k][nk[k]].x = geom->u1x[e2];
	  p[k][nk[k]++].y = geom->u1y[e2];
	}
	e = geom->zm1e[e];
      }
    }
    data->ec = n;
    /* Set the edges surrounding the edge vertices                   */
    for (i = 0; i < data->nes; i++) {
      e = data->eells[i];
      e2 = geom->m2de[e];
      while (e != geom->zm1e[e]) {
	k = geom->e2k[e];
	for (j = 0; j < 2; j++) {
	  v = geom->e2v[e][j];
	  vs = geom->m2dv[v];
	  for (vv = 1; vv <= geom->nve[vs]; vv++) {
	    en = geom->v2e[v][vv];
	    if (!mask[en]) {
	      e2 = geom->m2de[en];
	      data->eds[n] = nk[k];
	      data->e2k[n] = k;
	      data->eells[n++] = en;
	      mask[en] = nk[k] + 1;
	      p[k][nk[k]].x = geom->u1x[e2];
	      p[k][nk[k]].y = geom->u1y[e2];
	      nk[k]++;
	    }
	  }
	}
	e = geom->zm1e[e];
      }
    }

    /*---------------------------------------------------------------*/
    /* Fill the points array and create the triangulation for every  */
    /* layer. Note this is designed to be used with COMPAS arrays,   */
    /* which start at index 1, hence use [c+1].                      */
    data->de = (delaunay **)calloc(nz, sizeof(delaunay *));
    for (k = 0; k < nz; k++) {
      delaunay *dk;
      data->de[k] = NULL;
      if (nk[k] == 0) continue;

      data->de[k] = delaunay_build(nk[k], p[k], 0, NULL, 0, NULL);
      dk = data->de[k];

      dk->vid = 0;
      dk->ptf = 1;

      /* Free the point triangle arrays                              */
      if (dk->point_triangles != NULL) {
	for (i = 0; i < dk->npoints; ++i)
	  if (dk->point_triangles[i] != NULL)
	    free(dk->point_triangles[i]);
      }
      memset(dk->n_point_triangles, 0, dk->npoints * sizeof(int));

      for (e = 0; e < nk[k]; e++) {
	dk->points[e].v = d_alloc_1d(2);
      }
    }

    /* Set the new Delaunay point triangles for the interpolation    */
    for (i = 0; i < data->neells; i++) {
      delaunay *dk;
      e = data->eells[i];
      e2 = geom->m2de[e];
      k = geom->e2k[e];
      id = data->eds[i];
      dk = data->de[k];

      /* Edges in the triangulation; these point_trianges are the    */
      /* indices corresponding to all vertices in the triangulation  */
      /* (i.e. edge locations) whose triangles share a common edge   */
      /* index.                                                      */
      if (i < data->ec) {
	n = 0;
	dk->n_point_triangles[id] = 4 + 1;
	dk->point_triangles[id] = malloc(dk->n_point_triangles[id] * sizeof(int));
	dk->point_triangles[id][n++] = id;

	for (j = 0; j < 2; j++) {
	  v = geom->e2v[e][j]; /* Vertices of edge e                 */
	  vs = geom->m2dv[v];
	  for (vv = 1; vv <= geom->nve[vs]; vv++) {
	    int cc, jj, ec;
	    en = geom->v2e[v][vv]; /* Edges surrounding v            */
	    if (en == e) continue;
	    for (jj = 0; jj < 2; jj++) {
	      c = geom->e2c[e][jj]; /* Centres either side of e      */
	      c2 = geom->m2d[c];
	      for (cc = 1; cc <= geom->npe[c2]; cc++) {
		ec = geom->c2e[cc][c]; /* Edges surrounding c        */
		if (en == ec) {
		  dk->point_triangles[id][n++] = mask[en] - 1;
		  continue;
		}
	      }
	    }
	  }
	}
      } else {
	/* Cells surrounding the cell centres                        */
	dk->n_point_triangles[id] = 1;
	dk->point_triangles[id] = malloc(dk->n_point_triangles[id] * sizeof(int));
	dk->point_triangles[id][0] = -1;
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate and initialize                                       */
    ge = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    for (k = 0; k < nz; k++) {
      ge[k] = grid_spec_create();
      if (data->de[k] == NULL) continue;
      
      grid_interp_init_t(ge[k], data->de[k], uvrule, 0);
    }

    /*---------------------------------------------------------------*/
    /* Set the grid_spec structure                                   */
    data->ge = ge;



  /*-----------------------------------------------------------------*/
  /* Fill the points array in the delaunay structure.                */
  for (i = 0; i < data->ec; i++) {
    e = data->eells[i];
    k = geom->e2k[e];
    id = data->eds[i];
    /*
    data->de[k]->points[id].x = df->x[n2i[i]];
    data->de[k]->points[id].y = df->y[n2i[i]];
    */
    data->de[k]->points[id].z = (double)id;
  }

  /*-----------------------------------------------------------------*/
  /* Rebuild the weights. Note: this is cumulative (i.e. weights are */
  /* added to new cells) for each parray point.                      */
  for (k = 0; k < nz; k++) {
    if (data->de[k] == NULL) continue;
    data->ge[k]->rebuild(data->ge[k]->interpolator, data->de[k]->points);
  }

  i_free_1d(mask);
  i_free_1d(n2i);
}

/* END parray_gride_init()                                           */
/*-------------------------------------------------------------------*/


/* Simple grid centred file format */

static void df_memory_init_data(dump_data_t *dumpdata, dump_file_t *df);
static void df_memory_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int n, double *values, double *v2d);
static void df_memory_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz, double *v2d);
static void df_memory_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int n, double *values, double *v2d);
static void df_memory_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz, double *v2d);
static void df_memory_print(dump_file_t *df);
static void df_memory_dump(dump_file_t *df);

/* Structure to describe each dump file */
typedef struct {
  int fid;                      /* Netcdf file id */
  df_parray_var_t *vars;        /* List of point array dump file variables
                                 */
  df_mempack_t *data;           /* Output data structure */
  GRID_SPECS **gs;              /* Grid spec for parray interpolation */
  delaunay **d;                 /* Delaunay data structure for interpolation */
  int *cells;                   /* Cell locations used in the interpolation */
  int *ids;                     /* Delaunay indices used in the interpolation */
  int *c2k;                     /* Interpolation cells to layer map */
  int ncells;                   /* Number of cells used in the interpolation */
  int nz;                       /* Surface layer number */
  int ncs;                      /* Cell centres in surface layer */
  int nc;                       /* Total cell centres */

  GRID_SPECS **ge;              /* Grid spec for parray interpolation */
  delaunay **de;                /* Delaunay data structure for edge interpolation */
  int *eells;                   /* Edge locations used in the interpolation */
  int *eds;                     /* Delaunay indices used in the interpolation */
  int *e2k;                     /* Interpolation edge to layer map */
  int neells;                   /* Number of edgess used in the interpolation */
  int nes;                      /* Edges in surface layer */
  int ec;                       /* Total edge centres */

} df_memory_data_t;  

static void df_memory_init_data(dump_data_t *dumpdata, dump_file_t *df)
{
  int i;
  df_parray_data_t *mem = NULL;

  df_parse_vars(dumpdata,df,PARRAY_EXCLUDE_VARS,PARRAY_ALL_VARS);
  
  mem = (df_parray_data_t *)malloc(sizeof(df_parray_data_t));
  memset(mem, 0, sizeof(df_parray_data_t));
  df->private_data = mem;
  mem->type = DF_MEMORY;

  /* Allocate memory for variable information */
  if ((mem->vars =
       (df_parray_var_t *)malloc(df->nvars * sizeof(df_parray_var_t))) ==
      NULL)
    hd_quit("dumpfile_vars: Can't allocate memory for variable info.\n");

  /* Locate and store info for each variable */
  for (i = 0; i < df->nvars; i++) {
    df_parray_var_t *var = &mem->vars[i];
    if (df_parray_get_varinfo(dumpdata, df, df->vars[i], var) == 0)
      hd_quit("dumpfile:df_memory_init_data: Unknown variable '%s'.", df->vars[i]);
  }

  /* Allocate and initialize the grid_spec structures for            */
  /* interpolation.                                                  */
  parray_grid_init(dumpdata, df);
  parray_gride_init(dumpdata, df);

}

void *df_memory_create(dump_data_t *dumpdata, dump_file_t *df)
{
  char buf[MAXSTRLEN];
  int n2d, n3d;
  int i, n;
  df_parray_data_t *mem = NULL;
  df_mempack_t *data = NULL;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);

  if (df->npoints <= 0)
    hd_quit
      ("df_memory_create: No valid points have been specified for output.");


  df_memory_init_data(dumpdata, df);
  mem = (df_parray_data_t *)df->private_data;
  mem->data = (df_mempack_t *)malloc(sizeof(df_mempack_t));
  data = mem->data;
  data->runcode =  (dumpdata->runcode == RS_FAIL) ? DF_FAIL : DF_RUN;

  /* Get the number of variables and populate the headers */
  n2d = n3d = 0;    /* x and y */
  memset(data->vars, 0, sizeof(char) * MAXSTRLEN);
  memset(data->units, 0, sizeof(char) * MAXSTRLEN);
  memset(data->tunits, 0, sizeof(char) * MAXSTRLEN);
  memset(data->xyzunits, 0, sizeof(char) * MAXSTRLEN);

  for (n = 0; n < df->nvars; n++) {
    df_parray_var_t *var = &mem->vars[n];
    if (var->ndims == 2)
      n2d++;
    else if (var->ndims == 3)
      n3d++;
    strcpy(buf, var->units);
    strip(buf, " ");
    strcat(data->vars, var->name);
    strcat(data->vars, " ");
    strcat(data->units, buf);
    strcat(data->units, " ");
  }

  data->npoints = df->npoints;
  data->nz = df->nz;
  data->n2d= n2d;
  data->n3d= n3d;
  data->dt = df->sinc;
  strcpy(data->tunits, dumpdata->timeunit);
  data->xyz = d_alloc_1d(2 * df->npoints + df->nz);
  if (n2d) data->v2d = d_alloc_1d(n2d * df->npoints);
  if (n3d) data->v3d = d_alloc_1d(n3d * df->npoints * df->nz);

  /* Write the position data */
  if (is_geog)
    sprintf(data->xyzunits,"degrees_east degrees_north m");
  else
    sprintf(data->xyzunits,"%s %s m",dumpdata->lenunit, dumpdata->lenunit);

  if (is_geog) {
    n = 0;
    for (i = 0; i < df->npoints; i++) {
      data->xyz[n] = df->x[i];
      n++;
    }
    for (i = 0; i < df->npoints; i++) {
      data->xyz[n] = df->y[i];
      n++;
    }
    for (i = 0; i < df->nz; i++) {
      data->xyz[n] = dumpdata->cellz[i];
      n++;
    }
  } else if (has_proj) {
    char *args[256];
    int nargs = parseline(strdup(projection), args, 256);
    map_proj_t *mp = mp_init(nargs, args);
    double *lon = d_alloc_1d(df->npoints);
    double *lat = d_alloc_1d(df->npoints);

    for (i = 0; i < df->npoints; i++)
      mp_inverse(mp, df->x[i], df->y[i], &lat[i], &lon[i]);
    n = 0;
    for (i = 0; i < df->npoints; i++) {
      data->xyz[n] = lon[i];
      data->xyz[n] = df->x[i];
      n++;
    }
    for (i = 0; i < df->npoints; i++) {
      data->xyz[n] = lat[i];
      data->xyz[n] = df->y[i];
      n++;
    }
    for (i = 0; i < df->nz; i++) {
      data->xyz[n] = dumpdata->cellz[i];
      n++;
    }
    d_free_1d(lat);
    d_free_1d(lon);
  }
  df->finished = 0;
  return mem;
}


void df_memory_write(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int n;
  double newt = t;
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_mempack_t *data = mem->data;
  double *v2d, *v3d;

  /*
  tm_change_time_units(dumpdata->timeunit, df->tunit, &newt, 1);
  */
  data->time = newt;
  data->runcode =  (dumpdata->runcode == RS_FAIL) ? DF_FAIL : DF_RUN;

  /* Loop over each variable */
  v2d = data->v2d;
  v3d = data->v3d;
  for (n = 0; n < df->nvars; n++) {
    df_parray_var_t *var = &mem->vars[n];
    if (var->ndims == 2) {
      if (var->vector_mode == VM_NONE)
        df_memory_writesub_2d(dumpdata, df, n, get_data(var), v2d);
      else
        df_memory_writevec_2d(dumpdata, df, n, get_data(var), v2d);
      v2d += data->npoints;
    } else if (var->ndims == 3) {
      int klower = 0;
      int nz = 0;
      if (var->sediment) {
        nz = df->nz_sed;
        quit("Sed variable output not support for parray file format.");
      } else {
        klower = df->klower;
        nz = df->nz;
      }
      if (var->vector_mode == VM_NONE) {
        df_memory_writesub_3d(dumpdata, df, n, get_data(var),
                              klower, nz, v3d);
      } else {
        df_memory_writevec_3d(dumpdata, df, n, get_data(var),
                              klower, nz, v3d);
      }
      v3d += (data->npoints * data->nz);
    }
  }
  /* Print to ascii file if required */
  df_memory_dump(df);
}

void df_memory_close(dump_data_t *dumpdata, dump_file_t *df)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_mempack_t *data = mem->data;

  d_free_1d(data->xyz);
  if (data->v2d)
    d_free_1d(data->v2d);
  if (data->v3d)
    d_free_1d(data->v3d);
  df->finished = 1;
}


static void df_memory_writesub_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, double *v2d)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &mem->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = mem->d;
  GRID_SPECS **gs = mem->gs;
  int i, k, c, id;
  int fi, fj;

  /* Fill the Delaunay structure with values */
  for (i = 0; i < mem->ncells; i++) {
    c = mem->cells[i];
    if (c == geom->m2d[c]) {
      id = mem->ids[i];
      k = mem->c2k[i];
      d[k]->vid = 0;
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, c);
      else
	d[k]->points[id].v[0] = values[c];
      d[k]->points[id].z = d[k]->points[id].v[0];
      if (df->osl & (L_BAYLIN|L_BILIN)) {
	d[k]->points[id].z = (double)id;
	d[k]->points[id].v[0] = (double)id;
      }
    }
  }

  /* Rebuild the weights.                                            */
  k = geom->nz-1;
  if (df->osl & (L_SIB|L_NONSIB))
    for (i = 0; i < d[k]->npoints; i++)
      gs[k]->rebuild(gs[k]->interpolator, &d[k]->points[i]);
  else
    gs[k]->rebuild(gs[k]->interpolator, d[k]->points);
  if(df->osl & (L_BILIN|L_BAYLIN)) {
    for (i = 0; i < mem->ncells; i++) {
      c = mem->cells[i];
      if (c == geom->m2d[c]) {
	id = mem->ids[i];
	k = mem->c2k[i];
	if (df->filter)
	  d[k]->points[id].v[0] = df_filter(df, values, c);
	else
	  d[k]->points[id].v[0] = values[c];
      }
    }
  }

  for (i = 0; i < df->npoints; ++i) {
    id = mem->ids[i];
    v2d[i] = get_var_value(df, var, id, i, k);
  }
}



static void df_memory_writesub_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz, double *v3d)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &mem->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = mem->d;
  GRID_SPECS **gs = mem->gs;
  int i, k, c, id;
  int fi, fj;

  /* Fill the Delaunay structure with values */
  for (i = 0; i < mem->ncells; i++) {
    c = mem->cells[i];
    id = mem->ids[i];
    k = mem->c2k[i];
    d[k]->vid = 0;
    if (df->filter)
      d[k]->points[id].v[0] = df_filter(df, values, c);
    else
      d[k]->points[id].v[0] = values[c];
    d[k]->points[id].z = d[k]->points[id].v[0];
    if (df->osl & (L_BAYLIN|L_BILIN)) {
      d[k]->points[id].z = (double)id;
      d[k]->points[id].v[0] = (double)id;
    }
  }

  /* Rebuild the weights.                                            */
  for (k = 0; k < geom->nz; k++) {
    if (d[k] == NULL) continue;
    if (df->osl & (L_SIB|L_NONSIB))
      for (i = 0; i < d[k]->npoints; i++)
	mem->gs[k]->rebuild(mem->gs[k]->interpolator, &d[k]->points[i]);
    else
      mem->gs[k]->rebuild(mem->gs[k]->interpolator, d[k]->points);
  }
  if(df->osl & (L_BILIN|L_BAYLIN)) {
    for (i = 0; i < mem->ncells; i++) {
      c = mem->cells[i];
      id = mem->ids[i];
      k = mem->c2k[i];
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, c);
      else
	d[k]->points[id].v[0] = values[c];
    }
  }
  /*
  c = 0;
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < geom->nz; ++k){
      v3d[c++] = 0.0;
      v3d[c++] = get_missing(var);
    }
  }
  */
  c = 0;
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < geom->nz; ++k){
      v3d[c] = 0.0;
      if (d[k] != NULL) {
	id = mem->ids[i];
	v3d[c] = get_var_value(df, var, id, i, k);
      }
      c++;
    }
  }
}


static void df_memory_writevec_2d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, double *v2d)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &mem->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = mem->de;
  GRID_SPECS **gs = mem->ge;
  int i, k, e, id;

  /* Fill the Delaunay structure with values */
  for (i = 0; i < mem->neells; i++) {
    e = mem->eells[i];
    if (e == geom->m2de[e]) {
      id = mem->eds[i];
      k = mem->e2k[i];
      d[k]->vid = 0;
      if (df->filter)
	d[k]->points[id].v[0] = df_filter(df, values, e);
      else
	d[k]->points[id].v[0] = values[e];
      d[k]->points[id].z = d[k]->points[id].v[0];
    }
  }

  /* Rebuild the weights.                                            */
  k = geom->nz-1;
  gs[k]->rebuild(gs[k]->interpolator, d[k]->points);

  /*
  for(i = 0; i < df->npoints; i++) {
    v2d[i] = get_missing(var);
  }
  */
  for (i = 0; i < df->npoints; ++i) {
    id = mem->eds[i];
    v2d[i] = get_var_value(df, var, id, i, k);
  }
}

static void df_memory_writevec_3d(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz, double *v3d)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &mem->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = mem->de;
  GRID_SPECS **gs = mem->ge;
  int i, k, e, id;

  /* Fill the Delaunay structure with values */
  for (i = 0; i < mem->neells; i++) {
    e = mem->eells[i];
    id = mem->eds[i];
    k = mem->e2k[i];
    d[k]->vid = 0;
    d[k]->points[id].z = d[k]->points[id].v[0];
  }
}

static void df_memory_writevec_3do(dump_data_t *dumpdata, dump_file_t *df,
                                  int vid, double *values, int klower,
                                  int nz, double *v3d)
{
  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_parray_var_t *var = &mem->vars[vid];
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;
  delaunay **d = mem->de;
  GRID_SPECS **gs = mem->ge;
  int i, k, e, id;
  printf("a\n");
  /* Fill the Delaunay structure with values */
  for (i = 0; i < mem->neells; i++) {
    e = mem->eells[i];
    id = mem->eds[i];
    k = mem->e2k[i];
    d[k]->vid = 0;
    if (df->filter)
      d[k]->points[id].v[0] = df_filter(df, values, e);
    else
      d[k]->points[id].v[0] = values[e];
    d[k]->points[id].z = d[k]->points[id].v[0];
  }
  printf("b\n");
  /* Rebuild the weights.                                            */
  for (k = 0; k < geom->nz; k++) {
    if (d[k] == NULL) continue;
    mem->gs[k]->rebuild(mem->gs[k]->interpolator, d[k]->points);
  }
  printf("b\n");
  e = 0;
  for (i = 0; i < df->npoints; ++i) {
    for (k = 0; k < nz; ++k){
      v3d[e] = get_missing(var);
      if (d[k] != NULL) {
	id = mem->eds[i];
	v3d[e] = get_var_value(df, var, id, i, k);
      }
      e++;
    }
  }
  printf("c\n");
}


static void df_memory_print(dump_file_t *df)
{
  char buf[MAXSTRLEN];
  int i, k, m, n;
  FILE *fp;
  char *vars[256], *units[256];

  if (endswith(df->name, ".nc")) {
    n = strlen(df->name);
    for (i = 0; i < n-3; i++)
      buf[i] = df->name[i];
    buf[i] = '\0';
  } else
    strcpy(buf, df->name);
  n = strlen(buf);
  for (i = n - 1; i >= 0; i--)
    if (buf[i] == '.')
      break;
  buf[i] = '\0';
  if(!endswith(buf,".mpk"))
    strcat(buf,".mpk");
  if ((fp = fopen(buf, "w")) == NULL)
    hd_quit("Can't open file %s\n", df->name);

  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_mempack_t *data = mem->data;

  fprintf(fp, "Runcode = %d\n", data->runcode);
  fprintf(fp, "Number of points = %d\n", data->npoints);
  fprintf(fp, "Number of layers = %d\n", data->nz);
  fprintf(fp, "Number of 2d variables = %d\n", data->n2d);
  fprintf(fp, "Number of 3d variables = %d\n", data->n3d);
  fprintf(fp, "Time of dump = %f %s\n\n", data->time, data->tunits);

  n = parseline(strdup(data->xyzunits), units, 256);
  if (n != 3) hd_quit("Inconsistent xyz information in memory output data.\n");
  m = 0;
  fprintf(fp, "x (%s)\n", units[0]);
  for(n = 0; n < data->npoints; n++) {
    fprintf(fp, "%f ", data->xyz[m]);
    m++;
  }
  fprintf(fp, "\n\ny (%s)\n", units[1]);
  for(n = 0; n < data->npoints; n++) {
    fprintf(fp, "%f ", data->xyz[m]);
    m++;
  }
  fprintf(fp, "\n\nz (%s)\n", units[2]);
  for(n = 0; n < data->nz; n++) {
    fprintf(fp, "%f ", data->xyz[m]);
    m++;
  }
  n = parseline(strdup(data->vars), vars, 256);
  if (n != data->n2d + data->n3d) 
    hd_quit("Inconsistent variables numbers in memory output data (%d vs %d).\n", n, data->n2d + data->n3d);

  n = parseline(strdup(data->units), units, 256);
  if (n != data->n2d + data->n3d) 
    hd_quit("Inconsistent units numbers in memory output data (%d vs %d).\n", n, data->n2d + data->n3d);
  m = 0;
  for (i = 0; i < data->n2d; i++) {  
    fprintf(fp, "\n\n%s (%s)\n", vars[i], units[i]);
    for(n = 0; n < data->npoints; n++) {
      fprintf(fp, "%f ", data->v2d[m]);
      m++;
    }
  }
  m = 0;
  for (i = data->n2d; i < data->n2d + data->n3d; i++) {  
    fprintf(fp, "\n\n%s (%s)\n", vars[i], units[i]);
    for(n = 0; n < data->npoints; n++) {
      for(k = 0; k < data->nz; k++) {
	fprintf(fp, "%f ", data->v3d[m]);
	m++;
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
}


static void df_memory_dump(dump_file_t *df)
{
  char buf[MAXSTRLEN];
  int i, k, m, n;
  FILE *fp;
  char *vars[256], *units[256];

  if (endswith(df->name, ".nc")) {
    n = strlen(df->name);
    for (i = 0; i < n-3; i++)
      buf[i] = df->name[i];
    buf[i] = '\0';
  } else
    strcpy(buf, df->name);
  n = strlen(buf);
  for (i = n - 1; i >= 0; i--)
    if (buf[i] == '.')
      break;
  buf[i] = '\0';
  if(!endswith(buf,".mpk"))
    strcat(buf,".mpk");
  if ((fp = fopen(buf, "w")) == NULL)
    hd_quit("Can't open file %s\n", df->name);

  df_parray_data_t *mem = (df_parray_data_t *)df->private_data;
  df_mempack_t *data = mem->data;
  fprintf(fp, "mempack-version 1.0\n");
  fprintf(fp, "runcode %d\n", data->runcode);
  fprintf(fp, "npoints %d\n", data->npoints);
  fprintf(fp, "nz %d\n", data->nz);
  fprintf(fp, "n2d %d\n", data->n2d);
  fprintf(fp, "n3d %d\n", data->n3d);
  fprintf(fp, "time %f\n", data->time);
  fprintf(fp, "dt %f\n", data->dt);
  fprintf(fp, "tunits %s\n", data->tunits);
  fprintf(fp, "xyzunits %s\n", data->xyzunits);
  fprintf(fp, "vars %s\n", data->vars);
  fprintf(fp, "units %s\n", data->units);
  m = 2 * data->npoints + data->nz;
  fprintf(fp, "xyzdata %d\n", m);
  for(n = 0; n < m; n++)
    fprintf(fp, "%f\n", data->xyz[n]);
  m = data->n2d * data->npoints;
  fprintf(fp, "var2d %d\n", m);
  for(n = 0; n < m; n++)
    fprintf(fp, "%f\n", data->v2d[n]);
  m =  data->n3d * data->npoints * data->nz;
  fprintf(fp, "var3d %d\n", m);
  for(n = 0; n < m; n++)
    fprintf(fp, "%f\n", data->v3d[n]);
  fprintf(fp, "end\n");
  fflush(fp);
  fclose(fp);
}


/* Shift the longitude range from -180 - 180 deg to 0 - 360 deg. Note: this */
/* option should not be generally used since visualisation will not work    */
/* using the 0-360 range (coastline is not aligned).                        */
void set_longitude(dump_data_t *dumpdata, dump_file_t *df, int mode)
{
  int i, j, cc, c;
  double offset;
  master_t *master = dumpdata->master;
  geometry_t *geom = master->geom;

  if (df->long0_360 & (O360|O180)) {
    if (df->long0_360 & O360)
      offset = (mode) ? 180.0 : -180.0;
    if (df->long0_360 & O180)
      offset = (mode) ? -360.0 : 360.0;
    if (geom->us_type & US_IJ) {
      for (j = 0; j < dumpdata->nce2; j++) {
	for (i = 0; i < dumpdata->nce1; i++)
	  dumpdata->cellx[j][i] += offset;
	for (i = 0; i < dumpdata->nfe1; i++)
	  dumpdata->u1x[j][i] += offset;
      }
      for (j = 0; j < dumpdata->nfe2; j++) {
	for (i = 0; i < dumpdata->nfe1; i++)
	  dumpdata->gridx[j][i] += offset;
	for (i = 0; i < dumpdata->nce1; i++)
	  dumpdata->u2x[j][i] += offset;
      }
    } else {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	geom->cellx[c] += offset;
      }
      for (cc = 1; cc <= geom->n2_e1; cc++) {
	c = geom->w2_e1[cc];
	geom->u1x[c] += offset;
      }
      for (cc = 1; cc <= geom->n2_e2; cc++) {
	c = geom->w2_e2[cc];
	geom->gridx[c] += offset;
      }
    }
  }
}

/*------------------------------------------------------------------*/
/* Finds the time corresponding to the next month                   */
/*------------------------------------------------------------------*/
double next_year(double time, char *unit)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  sprintf(buf, "%d-01-01 00:00:00", yr+1);
  next = tm_datestr_to_julsecs(buf, unit);
  return(next-time);
}

/*------------------------------------------------------------------*/
/* Finds the time corresponding to the next season                  */
/*------------------------------------------------------------------*/
double next_season(double time, char *unit, int *smon)
{
  int season;
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  if (mon >= 3 && mon <= 5) {
    sprintf(buf, "%d-06-01 00:00:00", yr);
    season = AUTUMN;
  } else if (mon >= 6 && mon <= 8) {
    sprintf(buf, "%d-09-01 00:00:00", yr);
    season = WINTER;
  } else if (mon >= 9 && mon <= 11) {
    sprintf(buf, "%d-12-01 00:00:00", yr);
    season = SPRING;
  } else {
    yr = (mon == 12) ? yr+1 : yr;
    sprintf(buf, "%d-03-01 00:00:00", yr);
    season = SUMMER;
  }
  next = tm_datestr_to_julsecs(buf, unit);
  *smon = mon;
  return(next-time);
}

/*------------------------------------------------------------------*/
/* Finds the time corresponding to the end of the equivalent season */
/* in the previous year.                                            */
/*------------------------------------------------------------------*/
double prev_season(double time, char *unit, int *season, int *smon)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double prev;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  if (mon >= 3 && mon <= 5) {
    sprintf(buf, "%d-06-01 00:00:00", yr-1);
    *season = AUTUMN;
    *smon = 6;
  } else if (mon >= 6 && mon <= 8) {
    sprintf(buf, "%d-09-01 00:00:00", yr-1);
    *season = WINTER;
    *smon = 9;
  } else if (mon >= 9 && mon <= 11) {
    sprintf(buf, "%d-12-01 00:00:00", yr-1);
    *season = SPRING;
    *smon = 12;
  } else {
    yr = (mon == 12) ? yr+1 : yr;
    sprintf(buf, "%d-03-01 00:00:00", yr-1);
    *season = SUMMER;
    *smon = 3;
  }
  prev = tm_datestr_to_julsecs(buf, unit);
  return(prev);
}


/*------------------------------------------------------------------*/
/* Finds the time corresponding to the next month                   */
/*------------------------------------------------------------------*/
double next_month(double time, char *unit, int *smon)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  *smon = mon;
  if (mon == 12) {
    mon = 1;
    yr++;
  } else
    mon++;
  sprintf(buf, "%d-%d-01 00:00:00", yr, mon);
  next = tm_datestr_to_julsecs(buf, unit);
  return(next-time);
}


/*------------------------------------------------------------------*/
/* Finds the time corresponding to the end of the equivalent month  */
/* in the previous year.                                            */
/*------------------------------------------------------------------*/
double prev_month(double time, char *unit, int *smon)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double prev;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  if (mon == 12)
    mon = 1;
  else {
    mon++;
    yr--;
  }
  *smon = mon;
  sprintf(buf, "%d-%d-01 00:00:00", yr, mon);
  prev = tm_datestr_to_julsecs(buf, unit);
  return(prev);
}


/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Finds the time corresponding to the next day                     */
/*------------------------------------------------------------------*/
double next_day(double time, char *unit, int *sday)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  *sday = yrday(yr, mon, day);
  if (mon == 12 && day == 31) {
    day = 1;
    mon = 1;
    yr++;
  } else
    day++;
  sprintf(buf, "%d-%d-%d 00:00:00", yr, mon, day);
  next = tm_datestr_to_julsecs(buf, unit);
  return(next-time);
}

/*------------------------------------------------------------------*/
/* Finds the time corresponding to the end of the equivalent day    */
/* in the previous year.                                            */
/*------------------------------------------------------------------*/
double prev_day(double time, char *unit, int *sday)
{
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double prev;
  int yr, mon, day;

  date = tm_time_to_datestr(time, unit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
  if (mon == 12 && day == 31) {
    mon = 1;
    day = 1;
  } else {
    day++;
    yr--;
  }
  *sday = yrday(yr, mon, day);
  sprintf(buf, "%d-%d-%d 00:00:00", yr, mon, day);
  prev = tm_datestr_to_julsecs(buf, unit);
  return(prev);
}


/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Appends to the dumplist to output sparse transport files on a    */
/* monthly basis with a 1 hour dump increment. Only transport       */
/* variables are dumped.                                            */
/*------------------------------------------------------------------*/
int get_transfiles(dump_data_t *dumpdata, FILE * fp) {
  double time;
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int n, yr, mon, pmon, day;
  int st, et, nf;

  /* Get the month of the start time */
  date = tm_time_to_datestr(dumpdata->start_time, dumpdata->timeunit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &pmon, &day, buf);
  /* Get start and stop time as integers in seconds */
  st = (int)(dumpdata->start_time);
  et = (int)(dumpdata->stop_time);
  nf = 1;
  /* Daily loop from start to stop and count the months */
  for (n = st; n <= et; n += 86400) {
    time = (double)n;
    date = tm_time_to_datestr(time, dumpdata->timeunit);
    sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
    if (mon != pmon) nf++;
    pmon = mon;
  }
  return(nf);
}

int set_transfiles(int fn, char *key, dump_data_t *dumpdata, dump_file_t *list,
		   double t, FILE *fp) {
  FILE *tp;
  double time;
  char datestr[MAXSTRLEN], buf[MAXSTRLEN];
  char *date;
  double next;
  int n, m, yr, mon, pmon, nmon, day, tr;
  int st, et, nf;
  char *mons[13] = {"dum","jan","feb","mar","apr","may","jun","jul","aug",
		    "sep","oct","nov","dec"};

  /* Get the month of the start time */
  date = tm_time_to_datestr(dumpdata->start_time - 86400.0, dumpdata->timeunit);
  sscanf(date, "%d-%d-%d %s\n",&yr, &pmon, &day, buf);
  /* Get start and stop time as integers in seconds */
  st = (int)(dumpdata->start_time);
  et = (int)(dumpdata->stop_time);
  /* Set the list file number */
  nf = fn;
  /* Daily loop from start to stop and count the months */
  for (n = st; n <= et; n += 86400) {
    time = (double)n;
    date = tm_time_to_datestr(time, dumpdata->timeunit);
    sscanf(date, "%d-%d-%d %s\n",&yr, &mon, &day, buf);
    if (mon != pmon || n == st) {
      /* Set the filename */
      sprintf(list[nf].name,"%s_trans_%4d-%02d.nc", key, yr, mon);
      if (strlen(dumpdata->opath)) {
	sprintf(buf, "%s%s", dumpdata->opath, list[nf].name);
	strcpy(list[nf].name, buf);
      }
      /* Set the start time for the dump */
      list[nf].tout = (double)n;
      /* Set the stop time for the dump */
      for (m = n; m <= et; m += 86400) {
	time = (double)m;
	date = tm_time_to_datestr(time, dumpdata->timeunit);
	sscanf(date, "%d-%d-%d %s\n",&yr, &nmon, &day, buf);
	if (nmon != mon) break;
      }
      /*list[nf].tstop = (double)(m - 86400);*/
      list[nf].tstop = (double)(m);
      list[nf].append = forced_restart;
      /* Dump increment of 1 hour */
      if (prm_read_char(fp, "TRANS_DT", buf)) {
	tm_scale_to_secs(buf, &list[nf].tinc);
      } else
	list[nf].tinc = 3600.0;
      /* Fit the start and stop times to the actual simulation period */
      if (t > list[nf].tout) {
	list[nf].tout = ceil((t - list[nf].tout)
			    / list[nf].tinc)*list[nf].tinc + list[nf].tout;
      }      
      if (schedule->stop_time <= list[nf].tstop) {
	list[nf].tstop = floor((schedule->stop_time - list[nf].tstop)
			      / list[nf].tinc)*list[nf].tinc + list[nf].tstop;
      }
      list[nf].bpv = 8;
      /* Dump variables. Include swr if present (for bgc) */
      m = 0;
      for (tr = 0; tr < dumpdata->ntrS; tr++) {
	if (strcmp(dumpdata->trinfo_2d[tr].name, "swr") == 0) {
	  m = 1;
	}
      }
      if (m)
	strcpy(buf, "eta umean wmean temp salt Kzmean swr");
      else
	strcpy(buf, "eta umean wmean temp salt Kzmean");
      if (dumpdata->tmode & INEXACT) {
	if (m)
	  strcpy(buf, "eta u1 temp salt Kz origin p q r swr");
	else
	  strcpy(buf, "eta u1 temp salt Kz origin p q r");
      }
      if (dumpdata->tmode & SP_FFSL) {
	if (m)
	  strcpy(buf, "eta umean wmean temp salt Kzmean u1vmean swr");
	else
	  strcpy(buf, "eta umean wmean temp salt Kzmean u1vmean");
      }
      list[nf].nvars = parseline(strdup(buf), list[nf].vars, MAXNUMVARS);
      strcpy(list[nf].type, "sparse");
      list[nf].create = df_ugrid3_create;
      list[nf].write = df_ugrid3_write;
      list[nf].close = df_ugrid_close;
      list[nf].landfill = locate_landfill_function("default");
      
      list[nf].ilower = 0;
      list[nf].jlower = 0;
      list[nf].klower = 0;
      list[nf].nce1 = dumpdata->nce1;
      list[nf].nfe1 = dumpdata->nfe1;
      list[nf].nce2 = dumpdata->nce2;
      list[nf].nfe2 = dumpdata->nfe2;
      list[nf].nface2 = dumpdata->nface2;
      list[nf].nedge2 = dumpdata->nedge2;
      list[nf].nvertex2 = dumpdata->nvertex2;
      list[nf].nface3 = dumpdata->nface3;
      list[nf].nedge3 = dumpdata->nedge3;
      list[nf].nvertex3 = dumpdata->nvertex3;
      list[nf].nz = dumpdata->nz;
      list[nf].nz_sed = dumpdata->sednz;
      list[nf].ns2 = dumpdata->ns2;
      list[nf].ns3 = dumpdata->ns3;
      list[nf].da_cycle = NONE;
      list[nf].bathymask = 9999.0;
      list[nf].private_data = list[nf].create(dumpdata, &list[nf]);
      list[nf].finished = 0;
      strcpy(list[nf].tunit, dumpdata->output_tunit);
      nf++;
    }
    pmon = mon;
  }

  if ((tp = fopen("trans.mnc", "w")) == NULL)
    hd_warn("Can't create tranport .mnc file 'trans.txt'\n");
  else {
    fprintf(tp, "multi-netcdf-version 1.0\n");
    fprintf(tp, "nfiles %d\n", nf - fn);
    for (n = fn; n < nf; n++)
      fprintf(tp, "file%1.1d.filename  %s\n", n-fn, list[n].name);
    fclose(tp);
  }
  return(nf);
}

/* END get_transfiles()                                             */
/*------------------------------------------------------------------*/

/*
 * Gets the modified name of the output file based on the chunk type
 */
void set_chunk_name(dump_data_t *dumpdata, dump_file_t *df, double t)
{
  int yr, mon, day;
  double newt = t;

  /* Early exit, if there is no chunking */
  if (df->chunk == NONE)
    return;
  
  /* Key off the timeunit for dumping rather than master timeunit */
  tm_change_time_units(dumpdata->timeunit, df->tunit, &newt, 1);

  /* Figure out the chunk file name */
  if (!tm_time_to_ymd(newt, df->tunit, &yr, &mon, &day))
    hd_quit("dumpfile:set_chunk_name: Cannont determine y-m-d\n");

  if (df->chunk == YEARLY) {
    // eg. gbr_all-2011.nc
    sprintf(df->name, "%s_%4d.nc", df->origname, yr);
  } else if (df->chunk == MONTHLY) {
    // eg. gbr_all_2011-04.nc
    sprintf(df->name, "%s_%4d-%02d.nc", df->origname, yr, mon);
  } else if (df->chunk == DAILY) {
    // eg. gbr_all_2011-04-24.nc
    sprintf(df->name, "%s_%4d-%02d-%02d.nc", df->origname, yr, mon, day);
  }
  // else, as is
  
  /* Cache date stamp */
  df->curr_ymd[0] = yr;
  df->curr_ymd[1] = mon;
  df->curr_ymd[2] = day;
}

/*
 * Return the nc_mode for argument into ncw_create
 */
int get_nc_mode(dump_file_t *df)
{
  int nc_mode;

  /* Check the clobber flag */
  nc_mode = (overwrite_output() ? NC_CLOBBER : NC_NOCLOBBER);

  /* If we're compressing, need to specify netcdf4 file format */
  if (df->compress) {
#ifdef NC_NETCDF4
    /* set netCDF4 format */
    nc_mode |= NC_NETCDF4;
#else
    hd_quit("dumpfile_create: Compression flag has been specified for file %s but this netcdf library (version = %s) does not support it - must be atleast 4.x\n", df->name, nc_inq_libvers());
#endif
  } else {
    /* Make 64bit classical the default */
    if (!(master->compatible & V4201))
      nc_mode |= NC_64BIT_OFFSET;
  }
  return(nc_mode);
}

