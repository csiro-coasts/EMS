/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/master/hd_init.c
 *  
 *  Description:
 *  This file contains the main skeleton
 *  of the 3-d non-linear hydrodynamic
 *  model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd_init.c 6468 2020-02-18 23:47:02Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "hd.h"

#define NW 1
#define NE 2
#define SE 4
#define SW 8
#define EPS 1e-10

/*
#define nm3d 22
#define nm2d 5
char *mn_3d[nm3d] = {
  "u1mean", "u2mean", "wmean", "Kzmean", "temp_mean", "salt_mean",
  "flux_e1", "flux_e2", "flux_w", "flux_kz",
  "u1_adv", "u1_hdif", "u1_vdif", "u1_cor", "u1_btp", "u1_bcp",
  "u2_adv", "u2_hdif", "u2_vdif", "u2_cor", "u2_btp", "u2_bcp"
};
int mf_3d[nm3d] = {
  VEL3D, VEL3D, VEL3D, KZ_M, TS, TS,
  FLUX, FLUX, FLUX, FLUX,
  TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, 
  TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY, TENDENCY
};
char *mn_2d[nm2d] = {
  "eta_mean",
  "u1av_mean", "u2av_mean",
  "w1mean", "w2mean"
};
int mf_2d[nm2d] = {
  ETA_M,
  VEL2D, VEL2D,
  WIND, WIND
};
*/

void set_geo2index(master_t *master, geometry_t **window);
xytoij_tree_t *xytoij_init_sparse(geometry_t *geom);
void sigma_tracer_step_3d(master_t *master);
void inner_mean(double *var, int cs, int cd);
void init_fetch(master_t *master, parameters_t *params, unsigned long **flag);
int sign(double x);
void trans_sourcegrid_init(parameters_t *params, master_t *master, geometry_t *geom,
			   geometry_t **window, window_t **windat, win_priv_t **wincon,
			   dump_data_t *dumpdata);
void set_tp(geometry_t *geom, master_t *master, geometry_t **window);
int find_nearest(geometry_t *geom, double x, double y);
int oedge(int npe, int n);

/*------------------------------------------------------------------*/
/* Initialise the hydrodynamic structure                            */
/*------------------------------------------------------------------*/
hd_data_t *hd_init(FILE * prmfd)
{
  int icdfid = 0;               /* input cdf file id                */
  int timeindex = 0;            /* record index in dump file        */
  geometry_t *sgrid;


  /*----------------------------------------------------------------*/
  /* Read the input parameters for this run                         */
  if (autof == 0) {
    params = params_read(prmfd);
  } else if (autof == 8)
    params = params_read_t(prmfd);
  else {
    params = auto_params(prmfd, autof);
  }

  /*----------------------------------------------------------------*/
  /* Read in the initial conditions for -p, or reset the bathymetry */
  /* for -g.                                                        */
  if (params->runmode & (MANUAL | RE_ROAM | TRANS)) {
    /* Open input dump, get grid size and check compatibilities     */
    icdfid = dump_open_us(params, params->idumpname, 1);

    /* Using the model start time, confirm that we have a dump in   */
    /* the input file.                                              */
    timeindex = dump_choose_by_time(params, icdfid, schedule->start_time);

  } else {
    /* Set and check the bathymetry array                           */
    if (params->us_type & US_IJ)
      params->topo = set_bathy(params, icdfid);
  }

  /*----------------------------------------------------------------*/
  /* Make the mesh structure for -g. This is created from data in   */
  /* the input file for -p.                                         */
  if (params->runmode & (AUTO | DUMP)) {
    if (params->us_type & US_RS) {
      build_sparse_map(params);
      meshstruct_s(params, geom);
    } else
      meshstruct_us(params);
  }

  /*----------------------------------------------------------------*/   
  /* Create the unstructured mappings                               */
  sgrid = (geometry_t *)malloc(sizeof(geometry_t));
  memset(sgrid, 0, sizeof(geometry_t));
  TIMING_SET;
  build_sparse_grid_us(params, sgrid, params->nz, params->bathy,
		       params->layers);
  TIMING_DUMP(1, "  build_sparse_grid_us");

  /*----------------------------------------------------------------*/
  /* Set the automated parameters for -g                            */
  if (params->runmode & (AUTO | DUMP)) {
    /* Set the bottom depth on open boundaries                      */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);

    /* Set parameters to default values                             */
    autoset(params, master, sgrid);

    /* Reset the post-tracer ROAM parameterisation if required      */
    if (params->runmode & ROAM) autoset_roam(params, master, window);

    /* Set the ghost cells with valid data                          */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);
  }

  /*----------------------------------------------------------------*/
  /* Read the input dump data into the geometry and master          */
  if (params->runmode & (MANUAL | RE_ROAM | TRANS)) {
    dumpdata_read_us(geom, params, master, dumpdata, icdfid, timeindex);
    dump_close(icdfid);

    /* Initialise the velocity if required                          */
    vel_init(geom, params, master);

    /* Set the ghost and OBC cells with valid data                  */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);
  }

  if (params->runmode & RE_NOTS)
    temp_salt_init(params, master);

  /* Initialize the geometry and master structures */
  TIMING_SET;
  compute_constants(params, sgrid, master);
  TIMING_DUMP(1, " compute_constants");

#ifdef HAVE_WAVE_MODULE
  /* Initialise the fetch if required                                */
  if (params->do_wave)
    init_fetch(master, params, dumpdata->flag[geom->nz-1]);
#endif

  /* Initialize the mixing scheme */
  closure_init(params, master);

  /* Initialise the x,y to index conversion routine */
  set_geo2index(master, window);

  /* Set up the dumpfile structures */
  if (!(params->runmode & DUMP))
    dumpfile_setup(master);
  /* Dump and 'memory' output files */
  dump_m_snapshot(master);

  /* Initialize the diagostic time counters */
  monitor(master, NULL, -1);

  /* Set initial model time and output time */
  if (master->t >= schedule->stop_time)
    hd_quit_and_dump("Model start time is past stop time\n");

  /* Calculate the 2d time step */
  if (master->iratio < 0) {
    hd_warn("Illegal IRATIO value, using 1 instead\n");
    master->iratio = 1;
  }
  master->dt2d = master->grid_dt / master->iratio;
  master->dtb = master->dtf = master->grid_dt;
  master->dttc = 0.0;

  /* Initialze the boundary conditions */
  bdry_custom_m(params, sgrid, master);
  bdry_custom_w(sgrid, window);

  /* Initialise particle tracking */
  pt_params_init(master, master->prmfd);
  /* 
   * Initialise data assimilation - make sure this comes BEFORE
   * anything else that wants to use. eg. tracer_relax, timeseries and dumps
   */

#ifdef HAVE_DA
  dassim_init(master);
#endif
  /* Initialise the tracer relaxation.  */
  tracer_relax_init(master);

  /* Initialise the tracer resetting.  */
  tracer_reset_init(master);
  tracer_reset2d_init(master);

  /* Initialise the DHW diagnostic.  */
  tracer_dhw_init(master);

  /* Initialise the geometric arrays in the windows with master */
  /* geometry data.  */
  window_init(sgrid, window);

  /* Initialise the tracers */
  tracer_step_init(master);

  /* Initialise any external forcings */
  forcings_init(master);

  /* Initialise any tidal forcing */
  csr_tide_init(master, window);
  custom_tide_init(master, window, TD_ETA);
  custom_tide_init(master, window, TD_NU);
  custom_tide_init(master, window, TD_NV);
  custom_tide_init(master, window, TD_TU);
  custom_tide_init(master, window, TD_TV);

  /* Make and initialise the window data structures */
  windat = win_data_build(master, window);

  /* Make and initialise the window private data structures */
  wincon = win_consts_init(master, window);
  
  /* Calculate required initial conditions */
  TIMING_SET;
  pre_run_setup(master, window, windat, wincon);
  TIMING_DUMP(1, " pre_run_setup");

  /* Initialise the source/sink variables */
  sourcesink_init(params, master, window, windat, wincon);

  /* Initialise the distributed processing */
  dp_init(master, window, windat, wincon, geom->nwindows);

  /* Initialise the alerts log */
  alerts_init(master, window);

  /* Initialise the regions */
  init_regions_g(master, params);
  init_regions_w(master, window);

  /* Initialise the source grid for transport mode */
  if (params->runmode & TRANS)
    trans_sourcegrid_init(params, master, geom, window, windat, wincon, dumpdata);

  /* Initialise the timeseries events */
  timeseries_init(prmfd, master, geom, window, dumpdata);

  /* Write out the win_mp file */
  write_window_map(window, params);

  /* Set up 2 way nesting flags */
  init_2way(master, window, windat, wincon);

#ifdef HAVE_MPI
  if (mpi_check_multi_windows_sparse_arrays(geom))
    hd_quit("Single and multli-windows sparse arrays mismatch\n");
#endif
  

  hd_data = (hd_data_t *)malloc(sizeof(hd_data_t));
  hd_data->master = master;
  hd_data->geom = geom;
  hd_data->params = params;
  hd_data->dumpdata = dumpdata;
  hd_data->window = window;
  hd_data->windat = windat;
  hd_data->wincon = wincon;

  /* Initialise the crash recovery if required */
  if (crash_restart) crash_recovery_init(hd_data);
  bathy_compare(master);

  /* Write the output dump if required */
  trans_write(hd_data);
  if (params->runmode & (AUTO|DUMP)) {
    char tag[MAXSTRLEN], sag[MAXSTRLEN];

    if (params->runmode & AUTO)
      sprintf(tag, "%s.nc", params->oname);
    else
      strcpy(tag, params->oname);

    create_df(dumpdata, tag, US_WUS);
    if (dumpdata->us_type & US_WS) {
      sprintf(sag, "s_%s", tag);
      create_df(dumpdata, sag, US_WS);
    }

    if (params->runmode & AUTO)
      params_write(params, dumpdata);
    if (strlen(params->cookiecut)) cookie_cut(master, params);
    if (DEBUG("init_m"))
      dlog("init_m", "\nDumpfile %s created OK\n\n", tag);
  }

  write_mom_grid(dumpdata);
  write_roms_grid(dumpdata);
  if (strlen(params->wind_file))
    dump_windows_us(master, window, params->wind_file, params->prmname);

  if (params->runmode & EXIT) {
    /*windows_clear(hd_data);*/
    geom_free_us(master, sgrid, UNUSED);
    timeseries_end();
    master_end(master);
    sched_end(schedule);
    fclose(prmfd);
    master_free(master);
    debug_end();
    exit(0);
  }

  custom_init(hd_data);

  /* Register an interest with the scheduler for dumping the */
  /* output to file.  */
  if (!(params->runmode & DUMP)) {
    if (params->tmode & SP_CHECK)
      trans_check_dump(master, dumpdata, params->trans_data);
    else
      sched_register(schedule, "dumps", dump_init,
		     dump_event, dump_cleanup, master,
		     dump_dispatch, dump_progress);
  }
  /* Print the runtime diagnostics */
  write_run_setup(hd_data);
  if (strlen(master->opath)) {
    /* Write also to the output path */
    char buf[MAXSTRLEN];
    sprintf(buf, "%s%s", master->opath, setup_logfile);
    if (params->runno != 0.0)
      sprintf(buf, "%ssetup%2.1f.txt", master->opath, params->runno);
    strcpy(setup_logfile, buf);
    write_run_setup(hd_data);
  }

  /* Setup flag */
#ifdef HAVE_DA
  dassim_post_init(master);
#endif

  /* Register an interest with the scheduler for resetting   */
  /* the transport variables. This should always occur before */
  /* output dumps occur (hence register an interest after dumps). */
  if (params->runmode & TRANS) trans_reset_init(params, master);     

  geom_free_us(master, geom, UNUSED);
  if (DEBUG("init_m"))
    dlog("init_m", "\nInitialisation complete\n\n");

  return hd_data;
}

/* END hd_init()                                                    */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Close down the grid, freeing and allocate memory, etc.           */
/*------------------------------------------------------------------*/
void master_end(master_t *master)
{

  forcings_end();
  /*
  dp_cleanup();
  */
  sched_deregister(schedule, "dumps");

  tracer_reset_end(master);

  tracer_relax_end(master);

  ptrack_end();

}

/* END master_end()                                                  */
/*-------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Added, the scheduler was using the master after it was already   */
/*  dealocated (ie. forcings/wetbulb.c:wetbulb_cleanup)             */
/* Close down the grid, freeing and allocate memory, etc.           */
/* cleanup the late destroy stack                                   */
/*------------------------------------------------------------------*/
void master_free(master_t *master)
{
  free_stack_t* fs ;
  free_stack_t* fs2;
  if(master->free_stack != NULL)
  {
    emstag(LTRACE,"hd:hd_init:master_free","starting calling free stack ");
    fs = master->free_stack;
    do {
      fs2 = fs;
      fs = fs->next;
      if(fs2->free != NULL)
      {
        fs2->free(fs2->object);
        free(fs2);
      }else
      {
         free(fs2->object);
         free(fs2);
      }
    } while(fs != NULL);
    emstag(LDEBUG,"hd:hd_init:master_free","finished calling free stack ");
  }

/*UR-TODO - dmalloc fails on the master with error 22
 * meaning there is memory which was allocated by other means than malloc,
 * calloc, reallloc
 * */
  emstag(LDEBUG,"hd:hd_init:master_free","destroying master");
  if (master != NULL)
    free(master);
  emstag(LDEBUG,"hd:hd_init:master_free","destroyed master");
}

/* END master_free()                                                 */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Function to add a late reference for cleanup                      */
void free_stack_add(void* obj, void (*free) (void*) )
{
  free_stack_t* fs = (free_stack_t*) malloc(sizeof(free_stack_t));
  fs->last= NULL;
  fs->next= NULL;
  fs->object = obj;
  fs->free = free;
  if(master->free_stack == NULL)
  {
    master->free_stack = fs;
    master->free_stack->last = fs;
  }else
  {
    master->free_stack->last->next = fs;
    master->free_stack->last = fs;
  }
}
/* END free_stack_add()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up conversion trees, Delaunay structures and mapping         */
/* functions to convert geographic locations to index space.         */
/*-------------------------------------------------------------------*/
void set_geo2index(master_t *master, geometry_t **window)
{
  geometry_t *geom = master->geom;
  int n;

  if (geom->us_type & US_IJ) {
    master->xyij_tree = grid_xytoij_init(dumpdata->gridx, dumpdata->gridy,
					 geom->nce1, geom->nce2);
    geom->xyij_tree = master->xyij_tree;
  }
  if (geom->nwindows == 1) {
    window[1]->d = geom->d;
    window[1]->tri2c = geom->tri2c;
    window[1]->c2tri = geom->c2tri;
    window[1]->xyij_tree = master->xyij_tree;
  } else {
    for (n = 1; n <= geom->nwindows; n++) {
      if (window[n]->us_type & US_IJ) {
	window[n]->xyij_tree = master->xyij_tree;
	window[n]->map = geom->map;
      } else {
	/* Using xyztoc() for windows can be unsuccessful if the     */
	/* windows are discontinuous, as the triangulation will span */
	/* the 'gaps' of the discontinuity and may return a triange  */
	/* index that doesn't have a valid map to the mesh. Better   */
	/* to use the global triangulation then find which window    */
	/* the global cells resides in. See hd_grid_xyztoc_w in      */
	/* misc/index.c to see how this is handled.                  */
	window[n]->d = geom->d;
	window[n]->tri2c = geom->tri2c;
	window[n]->c2tri = geom->c2tri;
	/* This creates triangulations for each window
	create_delaunay_cell_w(window[n]);
	*/
      }
    }
  }
}

/* END set_geo2index()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the constants in the master geometry and data file */
/* from information read in through the input dumpfile and parameter */
/* file.                                                             */
/*-------------------------------------------------------------------*/
void compute_constants(parameters_t *params, /* Parameter structure  */
                       geometry_t *geom,     /* Global geometry      */
                       master_t *master      /* Master data          */
  )
{
  int c, cc, m;                 /* Centre coordinates                */
  int ee, e, es;                /* Edge coordinate                   */
  int vv, v;                    /* Vertex coordinate                 */
  int c1, c2, c3;               /* 2D and 3D sparse coordinates      */
  int cb;                       /* Bottom sparse coordinates         */
  double dh;                    /* Thickness of bottom cell          */
  double d1;                    /* Dummy                             */
  double top, bot;              /* Surface and bottom depths         */
  int n = 0;                    /* Mean grid spacing counter         */
  int tn;                       /* Tracer counter                    */
  int ns;                       /* Number of smoothing passes        */
  double *bfc;                  /* Bottom friction coefficient       */


  /* Set the master parameters from the parameter data structure     */
  /* Time stepping constants and variables                           */
  master->t = master->tstart = params->t;
  master->grid_dt = master->ogrid_dt = master->dt = params->grid_dt;
  master->iratio = params->iratio;
  master->tratio = params->tratio;
  master->trsplit = params->trsplit;
  master->nstep = master->nstep2d = 0;
  strcpy(master->timeunit, params->timeunit);
  strcpy(master->output_tunit, params->output_tunit);
  strcpy(master->prmname, params->prmname);
  strcpy(master->idumpname, params->idumpname);
  strcpy(master->opath, params->opath);
  master->mrtr = 0.0;

  /* Flags                                                           */
  master->trasc = params->trasc;
  master->momsc = params->momsc;
  master->momsc2d = master->momsc;
  master->ultimate = params->ultimate;
  master->osl = params->osl;
  master->smagorinsky = params->smagorinsky;
  master->sue1 = params->sue1;
  master->kue1 = params->kue1;
  master->bsue1 = params->bsue1;
  master->bkue1 = params->bkue1;
  master->smag_smooth = params->smag_smooth;
  master->diff_scale = params->diff_scale;
  master->visc_method = params->visc_method;
  master->stab = params->stab;
  master->thin_merge = params->thin_merge;
  master->sigma = params->sigma;
  master->nonlinear = params->nonlinear;
  master->calc_dens = params->calc_dens;
  master->roammode = params->roammode;
  if (strlen(params->densname)) {
    tn = tracer_find_index(params->densname, master->ntr, master->trinfo_3d);
    if (tn >= 0 && tn < master->ntr)
      master->calc_dens = -tn;
    else
      hd_warn("Can't find substitute density tracer %s.\n", params->densname);
  }
  master->slip = params->slipprm;
  master->mode2d = params->mode2d;
  master->tidef = params->tidef;
  master->tidep = params->tidep;
  master->nwindows = params->nwindows;
  master->dbc = 0;
  master->dbgtime = 0.0;
  m = 0;
  if(params->dbk >= 0 && params->dbj >= 0 && params->dbi >= 0) {
    if(params->dbk < params->nz && params->dbj < params->nce2 && 
       params->dbi < params->nce1) {
      master->dbc = geom->map[params->dbk][params->dbj][params->dbi];
      master->dbj = 1;
    } else
      hd_warn("Can't assign debug location (%d %d %d) to sparse location.\n",
	      params->dbi, params->dbj, params->dbk);
    m = 1;
  } else if (params->dbi >= 0 && params->dbi < geom->szc) {
    master->dbc = params->dbi;
    master->dbj = params->dbj;
    if (master->dbj > geom->npe[geom->m2d[master->dbc]]) master->dbj = 1;
    m = 1;
  }
  if (m) {
    if (params->dbgtime)
      master->dbgtime = params->t + params->dbgtime;
    master->dbgf = params->dbgf;
  }

  master->cfl = params->cfl;
  if (!(master->cfl & NONE)) {
    tm_scale_to_secs(params->cfl_dt, &master->cfl_dt);
  }

  if (master->nwindows > 1)
    master->win_reset = params->win_reset;
  else
    master->win_reset = 0;

  master->decorr = params->decorr;
  master->decf = params->decf;
  if (params->decf & DEC_ETA)
    master->decv = master->eta;
  else if (params->decf & DEC_U1)
    master->decv = master->u1;
  else {
    tn = tracer_find_index(params->decv, master->ntr, master->trinfo_3d);
    if (tn >= 0) {
      master->decv = master->tr_wc[tn];
      master->decn = tn;
    } else {
      master->decf = NONE;
    }
  }
  master->decs = 1.0;
  if (strcmp(params->decs, "km") == 0) master->decs = 1e3;

  master->monomn = params->monomn;
  master->monomx = params->monomx;
  master->monon = tracer_find_index(params->monotr, master->ntr, master->trinfo_3d);
  master->mixlayer = params->mixlayer;
  master->show_layers = params->show_layers;
  master->lnm = fabs(params->lnm);
  master->vorticity = params->vorticity;
  master->numbers = params->numbers;
  master->numbers1 = params->numbers1;
  master->tendf = params->tendf;
  master->robust = params->robust;
  master->fatal = params->fatal;
  master->trout = params->trout;
  master->tmode = params->tmode;
  master->fillf = params->fillf;
  master->pssinput = params->pssinput;
  master->conserve = params->conserve;
  master->do_closure = params->do_closure;
  master->lyear = params->lyear;
  master->trflsh = params->trflsh;
  strcpy(master->trage, params->trage);
  master->totals = params->totals;
  master->totals_dt = params->totals_dt;
  master->compatible = params->compatible;
  master->filter = params->filter;
  master->porusplate = params->porusplate;
  strcpy(master->reef_frac, params->reef_frac);
  master->gint_errfcn = params->gint_errfcn;
  master->dhwf = params->dhwf;
  master->dhwh = params->dhwh;
  master->swr_type = params->swr_type;
  master->togn = TOPRIGHT;
  master->crf = NONE;
  master->regf = NONE;
  strcpy(master->bathystats, params->bathystats);
  if (params->ntot) {
    master->ntot = params->ntot;
    master->totname = (char **)malloc(master->ntot * sizeof(char *));
    for (tn = 0; tn < master->ntot; tn++) {
      master->totname[tn] = (char *)malloc(sizeof(char)*MAXSTRLEN);	
      strcpy(master->totname[tn], params->totname[tn]);
    }
  }
  master->nprofn = -1;
  if (strlen(params->nprof))
    master->nprofn = tracer_find_index(params->nprof, master->ntr, master->trinfo_3d);
  master->nprofn2d = -1;
  if (strlen(params->nprof2d))
    master->nprofn2d = tracer_find_index(params->nprof2d, master->ntrS, master->trinfo_2d);

  master->trtend = -1;
  if (strlen(params->trtend)) {
    tn = tracer_find_index(params->trtend, master->ntr, master->trinfo_3d);
    if (tn >= 0 && tn < master->ntr) {
      master->trtend = tn;
      for (cc = 0; cc < master->ntr; cc++) {
	if (strcmp(master->trinfo_3d[cc].name, "tra_adv") == 0) {
	  sprintf(master->trinfo_3d[cc].name, "%s_adv", master->trinfo_3d[tn].name);
	  sprintf(master->trinfo_3d[cc].long_name, "%s advective tendency", master->trinfo_3d[tn].long_name);
	  sprintf(master->trinfo_3d[cc].units, "%s", master->trinfo_3d[tn].units);
	}
	if (strcmp(master->trinfo_3d[cc].name, "tra_hdif") == 0) {
	  sprintf(master->trinfo_3d[cc].name, "%s_hdif", master->trinfo_3d[tn].name);
	  sprintf(master->trinfo_3d[cc].long_name, "%s horizontal diffusive tendency", master->trinfo_3d[tn].long_name);
	  sprintf(master->trinfo_3d[cc].units, "%s", master->trinfo_3d[tn].units);
	}
	if (strcmp(master->trinfo_3d[cc].name, "tra_vdif") == 0) {
	  sprintf(master->trinfo_3d[cc].name, "%s_vdif", master->trinfo_3d[tn].name);
	  sprintf(master->trinfo_3d[cc].long_name, "%s vertical diffusive tendency", master->trinfo_3d[tn].long_name);
	  sprintf(master->trinfo_3d[cc].units, "%s", master->trinfo_3d[tn].units);
	}
	if (strcmp(master->trinfo_3d[cc].name, "tra_ncon") == 0) {
	  sprintf(master->trinfo_3d[cc].name, "%s_ncon", master->trinfo_3d[tn].name);
	  sprintf(master->trinfo_3d[cc].long_name, "%s non-conservative tendency", master->trinfo_3d[tn].long_name);
	  sprintf(master->trinfo_3d[cc].units, "%s", master->trinfo_3d[tn].units);
	}
      }
    } else
      hd_warn("Can't find tracer %s for tendencies.\n", params->trtend);
  }

  /* Means                                                           */
  init_means(master, params);

  if (master->mono && master->monomn == master->monomx) {
    master->monomn = HUGE;
    master->monomx = -HUGE;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      master->monomn = min(master->monomn, master->tr_wc[master->monon][c]);
      master->monomx = max(master->monomx, master->tr_wc[master->monon][c]);
    }
  }

  /* Particles                                                       */
  master->do_pt = params->do_pt;
  if (master->do_pt)
    strcpy(master->ptinname, params->ptinname);

  master->u1_f = params->u1_f;
  master->u1av_f = params->u1av_f;
  master->orbital = params->orbital;
  master->waves = params->waves;
  master->save_force = params->save_force;
  master->exmapf = params->exmapf;
  geom->nwindows = params->nwindows;
  master->show_win = params->show_win;
  strcpy(master->orthoweights, params->orthoweights);
  strcpy(master->nodal_dir, params->nodal_dir);
  strcpy(master->tide_con_file, params->tide_con_file);
  strcpy(master->bdrypath, params->bdrypath);
  master->restart_dt = params->restart_dt;
  strcpy(master->restart_name, params->restart_name);
  master->runmode = VEL3D;
  master->sh_f = NONE;
  if (master->mode2d)
    master->runmode = VEL2D;
  if(params->runmode & TRANS)
    master->runmode = TRANS;
  if(params->runmode & DUMP)
    master->runmode |= DUMP;
  master->vinit = NONE;

  if (strlen(params->vel_init)) {
    master->vinit = VINIT_FLOW;
    if (strcmp(params->vel_init, "GEOSTROPHIC") == 0) master->vinit = VINIT_GEO;
  }

  /* Constants                                                       */
  master->g = 9.81;
  master->ambpress = params->ambpress;
  master->hmin = params->hmin;
  master->uf = params->uf;
  /*master->quad_bfc = params->quad_bfc;*/
  master->etamax = params->etamax;
  master->velmax = params->velmax;
  master->velmax2d = params->velmax2d;
  master->etadiff = params->etadiff;
  master->wmax = params->wmax;
  master->rampstart = params->rampstart;
  master->rampend = params->rampend;
  master->rampf = params->rampf;
  master->eqt_alpha = params->eqt_alpha;
  master->eqt_beta = params->eqt_beta;

  /* Tracer constants and variables                                  */
  master->trperc = -1;
  master->trflux = -1;
  master->vf = NULL;
  for (tn = 0; tn < master->ntr; tn++) {
    tracer_info_t *tracer = &master->trinfo_3d[tn];
    master->advect[tn] = tracer->advect;
    master->diffuse[tn] = tracer->diffuse;
    master->fill_value_wc[tn] = tracer->fill_value_wc;
    master->mintr[tn] = tracer->valid_range_wc[0];
    master->maxtr[tn] = tracer->valid_range_wc[1];
    strcpy(master->trname[tn], tracer->name);
    if (strcmp(params->trflux, tracer->name) == 0) {
      master->trflux = tn;
      master->trfd1 = params->trfd1;
      master->trfd2 = params->trfd2;
    }
    if (strcmp(params->trperc, tracer->name) == 0)
      master->trperc = tn;
  }
  strcpy(master->trpercr, params->trpercr);
  if ((strcmp(params->trperc, "NONE") != 0) && master->trperc < 0)
    hd_quit("Can't find tracer %s for percentiles.\n", params->trperc);
  if ((strcmp(params->trflux, "NONE") != 0) && master->trflux < 0)
    hd_quit("Can't find tracer %s for flux diagnostic.\n", params->trflux);

  for (tn = 0; tn < master->ntdif_h; tn++)
    master->tdif_h[tn] = params->tdif_h[tn];
  for (tn = 0; tn < master->ntdif_v; tn++)
    master->tdif_v[tn] = params->tdif_v[tn];

  /* Transport mode : get additional variables to be reset           */
  master->ntrvars = master->ntrvarsS = 0;
  if (strlen(params->trvars)) {
    int nfiles, tt, *tmap, *tmapS;
    char files[MAXNUMTSFILES][MAXSTRLEN], trfile[MAXSTRLEN];
    nfiles = parseline(params->trvars, (char **)files, MAXNUMTSFILES);
    if(nfiles) {
      tmap = i_alloc_1d(nfiles);
      tmapS = i_alloc_1d(nfiles);
      /* Count the number of valid tracers */      
      for (tt = 0; tt < nfiles; ++tt) {
	strcpy(trfile, ((char **)files)[tt]);
	tmap[tt] = -1;
	tmapS[tt] = -1;
	tn = tracer_find_index(trfile, master->ntr, master->trinfo_3d);
	if (tn >= 0) {
	  if ((strcmp(trfile, "temp") == 0) &&
	      (strcmp(trfile, "salt") == 0)) {
	    master->trinfo_3d[tn].advect = 0;
	    master->trinfo_3d[tn].diffuse = 0;
	  } else {
	    master->ntrvars++;
	    tmap[tt] = tn;
	  }
	} else if ((tn = tracer_find_index(trfile, master->ntrS, master->trinfo_2d)) >= 0) {
	  master->ntrvarsS++;
	  tmapS[tt] = tn;
	} else {
	  hd_warn("Reset transport variable %s is not included in tracer list.\n",
		  trfile);
	}
      }
      /* Fill the arrays with valid tracers                          */
      /* 3D variables                                                */
      if (master->ntrvars) {
	master->trvars =
	  (cstring *) malloc(sizeof(cstring) * master->ntrvars);
	master->trvm = i_alloc_1d(master->ntrvars);
	for (tt = 0; tt < master->ntrvars; ++tt)
	  master->trvm[tt] = -1;
	master->ntrvars = 0;
	for (tt = 0; tt < nfiles; ++tt) {
	  tn = tmap[tt];
	  if (tn >= 0) {
	    strcpy(master->trvars[master->ntrvars], ((char **)files)[tt]);
	    master->trvm[master->ntrvars] = tn;
	    master->trinfo_3d[tn].advect = 0;
	    master->trinfo_3d[tn].diffuse = 0;
	    master->ntrvars++;
	  }
	}
      }
      i_free_1d(tmap);
      /* 2D variables                                                */
      if (master->ntrvarsS) {
	master->trvarsS =
	  (cstring *) malloc(sizeof(cstring) * master->ntrvarsS);
	master->trvmS = i_alloc_1d(master->ntrvarsS);
	for (tt = 0; tt < master->ntrvarsS; ++tt)
	  master->trvmS[tt] = -1;
	master->ntrvarsS = 0;
	for (tt = 0; tt < nfiles; ++tt) {
	  tn = tmapS[tt];
	  if (tn >= 0) {
	    strcpy(master->trvarsS[master->ntrvarsS], ((char **)files)[tt]);
	    master->trvmS[master->ntrvarsS] = tn;
	    master->ntrvarsS++;
	  }
	}
      }
      i_free_1d(tmapS);
    }
  }

  /* Mixing constants and variables                                  */
  master->vz0 = params->vz0;
  master->kz0 = params->kz0;
  master->zs = params->zs;
  master->min_tke = params->min_tke;
  master->min_diss = params->min_diss;
  master->Lmin = params->Lmin;
  master->smooth_VzKz = params->smooth_VzKz;
  master->eparam = params->eparam;
  master->kz_alpha = params->kz_alpha;
  master->vz_alpha = params->vz_alpha;
  master->wave_alpha = params->wave_alpha;
  master->wave_b1 = params->wave_b1;
  master->wave_hf = params->wave_hf;

  /* Heatflux variables                                              */
  master->heatflux = params->heatflux;
  master->albedo = params->albedo;
  master->albedo_l = params->albedo_l;
  master->zref = params->zref;
  master->bulkf = params->bulkf;
  master->hfadf = params->hfadf;
  master->hftc = params->hftc;
  master->hf_ramp = params->hf_ramp;
  /* Override input file swr parameters if they are numbers          */
  if (master->swr_attn) {
    double d1;
    int cc, c;
    if (sscanf(params->swr_attn, "%lf", &d1) == 1) {
      if (params->swr_type & SWR_2D) {
        for (cc = 1; cc <= geom->b2_t; cc++) {
	  c = geom->w2_t[cc];
	  master->swr_attn[c] = d1;
	}
      }
      if (params->swr_type & SWR_3D) {
        for (cc = 1; cc <= geom->b3_t; cc++) {
	  c = geom->w3_t[cc];
	  master->swr_attn[c] = d1;
	}
      }
    }
    if (master->swr_attn1 && sscanf(params->swr_attn1, "%lf", &d1) == 1) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	master->swr_attn1[c] = d1;
      }
    }
    if (master->swr_tran && sscanf(params->swr_tran, "%lf", &d1) == 1) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	master->swr_tran[c] = d1;
      }
    }
    if (master->swr_babs && sscanf(params->swr_babs, "%lf", &d1) == 1) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	master->swr_babs[c] = d1;
      }
    }
  }

  /* Saltflux variables                                              */
  master->saltflux = params->saltflux;

  /* Wind variables                                                  */
  master->wind_dt = params->wind_dt;
  master->storm_dt = params->storm_dt;
  master->nstorm = params->nstorm;

  /* Elevation relaxation variables */
  master->etarlx = params->etarlx;
  if (params->etarlx & (RELAX|ALERT|BOUNDARY|ETA_TPXO)) {
    master->eta_rlx = relax_info_init(params->etarlxn, params->etarlxtcs, 
				      params->etarlxdt, 0, 0);
    if (params->etarlx & (RELAX|ALERT|BOUNDARY)) {
      if ((tn = tracer_find_index("oeta", master->ntrS, master->trinfo_2d)) >= 0)
	master->eta_rlx->val1 = master->tr_wcS[tn];
    }
    if (params->etarlx & ETA_TPXO) {
      if ((tn = tracer_find_index("tpxotide", master->ntrS, master->trinfo_2d)) >= 0) {
	master->eta_rlx->val1 = master->tr_wcS[tn];
      } else
	hd_quit("Can't find TPXO eta relaxation tracer: use NUMBERS TPXO.\n");
    }
    master->eta_rlx->rlx = params->etarlx;
  }

  /* Velocity relaxation                                             */
  master->velrlx = params->velrlx;
  if (params->velrlx & RELAX) {
    master->vel_rlx = relax_info_init(params->velrlxn, params->velrlxtcs, 
				      params->velrlxdt, geom->sgsiz, geom->sgsiz);
    master->vel_rlx->rlx = params->velrlx;
  }
  master->tide_r = params->tide_r;

  /* Alerts                                                          */
  master->alertf = NONE;
  master->alert_dt = 0.0;
  master->ef = NULL;
  if (strlen(params->alert)) {
    char *fields[MAXSTRLEN * 2];
    int nf = parseline(params->alert, fields, 2);
    if (strcmp(fields[0], "ACTIVE") == 0)
      master->alertf = ACTIVE;
    if (strcmp(fields[0], "PASSIVE") == 0)
      master->alertf = PASSIVE;

    if(nf > 1)
      sprintf(master->alertname, "%s", fields[1]);
    else
      strcpy(master->alertname, "alert");
    if (strlen(params->alert_dt))
      tm_scale_to_secs(params->alert_dt, &master->alert_dt);

    master->eta_f = params->eta_f;
    master->vel2d_f = params->vel2d_f;
    master->vel3d_f = params->vel3d_f;
    master->wvel_f = params->wvel_f;
    master->tend_f = params->tend_f;
    master->div2d_f = params->div2d_f;
    master->div3d_f = params->div3d_f;
    master->cfl_f = params->cfl_f;
    master->ts_f = params->ts_f;
    master->shear_f = params->shear_f;
    master->hdiff_f = params->hdiff_f;
    master->amax = params->amax;
    master->hmax = params->hmax;
    master->vmax = params->vmax;
    master->btmax = params->btmax;
    master->bcmax = params->bcmax;
    master->cmax = params->cmax;
    master->detamax = params->detamax;
    master->dwmax = params->dwmax;
    master->dtmax = params->dtmax;
    master->dsmax = params->dsmax;
    master->smax = params->smax;
    if (master->alertf & ACTIVE) {
      if ((master->vel2d_f || master->vel3d_f || master->wvel_f ||
	   master->div2d_f || master->div3d_f) &&
	  master->smagcode & (U1_A|U2_A))
	hd_warn("Active diffusion ALERT requires Smagorinsky diffusion\n");
      if ((master->eta_f || master->div2d_f) && !(master->etarlx & ALERT))
	hd_warn("Active relaxation ALERT requires eta_relaxation_file\n");
      if ((master->eta_f || master->div2d_f) && master->tide_r & MEAN_R &&
	  !(master->means & ETA_M))
	hd_warn("Active relaxation ALERT requires mean elevation specified\n");
    }
  }

#if defined(HAVE_SEDIMENT_MODULE)
  /* Sediments                                                       */
  master->do_sed = params->do_sed;
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  /* Ecology                                                         */
  master->do_eco = params->do_eco;
  master->ecodt = params->ecodt;
#endif

#if defined(HAVE_WAVE_MODULE)
  /* Waves                                                           */
  master->do_wave = params->do_wave;
  master->wavedt = params->wavedt;
#endif

  /*-----------------------------------------------------------------*/
  /* Initialize the sigma variables                                  */
  for (c = 1; c <= geom->enonS; c++) {  /* Sigma variables           */
    /*
       master->Ds[c]=1.0; master->Hs[c]=0.0; master->Hn1[c]=1.0;
       master->Hn2[c]=1.0; */
  }

  /*-----------------------------------------------------------------*/
  /* Smooth the topography if required                               */
  if (params->sigma) {
    /* testtopo(params,master); */
    for (cc = 1; cc <= geom->a2_t; cc++) {
      c = geom->w2_t[cc];
      master->Ds[c] = master->eta[c] - geom->botz[c];
      master->Hs[c] = geom->botz[c];
      master->topz[c] = 0.0;
      geom->botz[c] = -1.0;
    }
    memcpy(master->dz, geom->cellz, geom->sgsiz * sizeof(double));
    for (cc = 1; cc <= geom->a3_t; cc++) {
      c = geom->w3_t[cc];
      c2 = geom->s2k[c];
      geom->gridz[c] = params->layers[c2];
      geom->cellz[c] = 0.5 * (params->layers[c2] + params->layers[c2 + 1]);
    }
    sigma_tracer_step_3d(master);
  }
  if (params->show_layers)
    save_sigma_layers(params, master);

  /*-----------------------------------------------------------------*/
  /* Set the cell thickness for the master                           */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = c3 = geom->w2_t[cc];
    c2 = geom->m2d[c3];
    cb = geom->bot_t[cc];
    top = master->eta[c2];

    while (c3 != cb) {
      bot = geom->gridz[c3] * master->Ds[c2];
      master->dz[c3] = top - bot;
      top = bot;
      c3 = geom->zm1[c3];
    }

    /* Set the cell thickness at the bottom */
    master->dz[cb] = top - geom->botz[c2] * master->Ds[c2];
  }

  /*-----------------------------------------------------------------*/
  /* Set the surface roughness                                       */
  if (params->z0s) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->z0[c] = params->z0s[geom->c2cc[c]];
    }
  } else {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      c2 = geom->m2d[c];
      master->z0[c2] = params->z0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the drag coefficient                                        */
  bfc = d_alloc_1d(geom->szcS);
  d1 = 0.003;
  if (sscanf(params->quad_bfc, "%lf", &dh) == 1)
    d1 = atof(params->quad_bfc);
  value_init_2d(master, bfc, params->prmfd, params->quad_bfc,
		  "Cd", "QBFC", d1, "linear"); 

  master->quad_bfc = 0.0;
  for (cc = 1; cc < geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    master->Cd[c] = bfc[c];
    master->quad_bfc += bfc[c];
  }
  master->quad_bfc /= (double)geom->b2_t;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    c2 = geom->m2d[c];
    cb = geom->bot_t[cc];
    dh = max(0.5 * master->dz[cb] * master->Ds[c2], master->hmin);
    if (master->mode2d) dh = 0.5 * (master->eta[c2] - geom->botz[c2]);
    d1 = log((dh + master->z0[c2]) / master->z0[c2]) / VON_KAR;
    /*if (master->quad_bfc < 0.0) {*/
    if (bfc[c] < 0.0) {
      if (master->mode2d) {
	/* Assume quad_bfc is the Manning coefficient, n[s.m^(-1/3)] */
	/* , where Cd = g.n^2.D^(-1/3), e.g. see:                    */
	/* http://www.marinespecies.org/introduced/wiki/Bed_roughness_and_friction_factors_in_estuaries */
	/* This makes the bottom drag depth dependent.               */
	dh = master->eta[c2] - geom->botz[c2];
	d1 = -1.0 / 3.0;
	master->Cd[c2] = master->g * bfc[c] * bfc[c] * pow(dh, d1);
      } else
	master->Cd[c2] = fabs(bfc[c]);
    } else
      master->Cd[c2] = max(bfc[c2], 1.0 / (d1 * d1));
  }

  ns = get_smoothing(params->smooth_v, "CD");
  for (n = 0; n < ns; n++)
    smooth3(master, master->Cd, geom->w2_t, geom->b2_t, geom->szcS, 0);
  /*smooth(master, master->Cd, geom->w2_t, geom->b2_t);*/
  if ((d1 = get_scaling(params->scale_v, "CD")) != 1.0) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->Cd[c] *= d1;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for cell centers                        */
  master->hmean1 = master->hmean2 = master->amean = master->edmean = master->maxres = 1e-10;
  master->minres = 1e10;
  geom->totarea = 0.0;
  for (cc = 1; cc <= geom->n2_t; cc++) {
    c = geom->w2_t[cc];

    if (cc <= geom->b2_t) 
      geom->totarea += geom->cellarea[c];

    geom->sinthcell[c] = 0.0;
    geom->costhcell[c] = 0.0;
    if (master->cellres)
      master->cellres[c] = 0.0;
    for (m = 1; m <= geom->npe[c]; m++) {
      e = geom->c2e[m][c];
      geom->sinthcell[c] += sin(geom->thetau1[e]);
      geom->costhcell[c] += cos(geom->thetau1[e]);
      if (master->cellres)
	master->cellres[c] += (2.0 * geom->hacell[m][c]); 
    }
    geom->sinthcell[c] /= (double)geom->npe[c];
    geom->costhcell[c] /= (double)geom->npe[c];
    if (master->cellres)
      master->cellres[c] /= (double)geom->npe[c];
    geom->dHde1[e] = 0.5 * (geom->botz[geom->e2c[e][0]] - geom->botz[geom->e2c[e][1]]);
  }
  master->amean = geom->totarea / (double)geom->b2_t;
  dh = 0.0;
  for (ee = 1; ee <= geom->b2_e1; ee++) {
    e = geom->w2_e1[ee];
    dh += geom->edgearea[e];
  }
  master->edmean = dh / (double)geom->b2_e1;

  /*-----------------------------------------------------------------*/
  /* Scale horizontal mixing coefficients by cell size.              */ 
  /* Cell centered horizontal diffusivity.                           */
  memset(master->u1kh, 0, geom->szc * sizeof(double));
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    if (params->u1kh > 0.0) {
      if (params->diff_scale & NONE)
	master->u1kh[c] = fabs(params->u1kh);
      if (params->diff_scale & LINEAR) {
	master->u1kh[c] = fabs(params->u1kh * sqrt(geom->cellarea[c2]) /
			       sqrt(master->amean));
      }
      if (params->diff_scale & (NONLIN|CUBIC|AREAL))
	master->u1kh[c] = fabs(params->u1kh * 
			       geom->cellarea[c2] / master->amean);
    }
  }

  master->u1vh0 = params->u1vh;
  master->u2vh0 = params->u2vh;
  master->u1kh0 = params->u1kh;
  master->u2kh0 = params->u2kh;
  memset(geom->dHde1, 0, geom->szeS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for u1 momentum equations               */
  for (ee = 1; ee <= geom->v2_e1; ee++) {
    e = geom->w2_e1[ee];
    master->hmean1 += geom->h1au1[e];
    master->hmean2 += geom->h2au1[e];
    if (geom->h2au1[e] < master->minres) {
      master->minres = geom->h2au1[e];
      master->minrese = e;
    }
    if (geom->h2au1[e] > master->maxres) {
      master->maxres = geom->h2au1[e];
      master->maxrese = e;
    }
  }
  master->hmean1 /= (double)geom->v2_e1;
  master->hmean2 /= (double)geom->v2_e1;

  for (ee = 1; ee <= geom->n2_e1; ee++) {
    int em, ep;
    double sgn;
    double eps = 1e-15;
    e = geom->w2_e1[ee];
    c1 = geom->e2c[e][0];
    c2 = geom->e2c[e][1];
    sgn = (geom->e2e[e][0] == 1) ? 1.0 : -1.0;

    master->u1c1[e] =
      -1.0 / (geom->h2au1[e] * geom->h2au1[e] * geom->h1au1[e]);

    em = e2e(geom, e, geom->e2e[e][0]);
    ep = e2e(geom, e, geom->e2e[e][1]);
    master->u1c3[e] =
      1.0 * (geom->h1acell[e] -
             geom->h1acell[em]) / (geom->h2au1[e] * geom->h2au1[e]);
    master->u1c4[e] =
      1.0 * (geom->h1au1[ep] -
             geom->h1au1[em]) / (geom->h1au1[e] * geom->h2au1[e]);
    master->u1c5[e] =
      sgn * (master->coriolis[c1] + master->coriolis[c2]) / 2;
    master->u1c6[e] = -1.0 / geom->h2au1[e];
    geom->sinthu1[e] = sin(geom->thetau1[e]);
    geom->costhu1[e] = cos(geom->thetau1[e]);
    geom->sinthu2[e] = sin(geom->thetau2[e]);
    geom->costhu2[e] = cos(geom->thetau2[e]);
    if (fabs(geom->sinthu1[e]) < eps) geom->sinthu1[e] = 0.0;
    if (fabs(geom->costhu1[e]) < eps) geom->costhu1[e] = 0.0;
    if (fabs(geom->sinthu2[e]) < eps) geom->sinthu2[e] = 0.0;
    if (fabs(geom->costhu2[e]) < eps) geom->costhu2[e] = 0.0;
    /*
    if (ee > geom->v2_e1 && ee <= geom->b2_e1)
      geom->botzu1[e] = min(geom->botz[c1], geom->botz[c2]);
    else
    */
      geom->botzu1[e] = max(geom->botz[c1], geom->botz[c2]);
  }
  /* Set the normal boundary cells to the cell centre depth          */
  for (m = 0; m < geom->nobc; m++) {
    open_bdrys_t *open = geom->open[m];
    /* Set the flux adjustment as a ratio of 2D timestep if required */
    if (open->afr) {
      open->adjust_flux = open->afr * master->grid_dt / master->iratio;
      hd_warn("OBC %s uses flux adjustment = %f s\n", open->name, open->adjust_flux);
    }
    for (ee = 1; ee <= open->no2_e1; ee++) {
      e = open->obc_e1[ee];
      c = open->obc_e2[ee];
      geom->botzu1[e] = geom->botz[c];
    }
  }

  /* Horizontal viscosity                                            */
  if (params->u1vh > 0.0 || params->bsue1 > 0.0) {
    double vh = (params->u1vh > 0.0) ? params->u1vh : params->bsue1;
    int etype = 1;   /* 0 : use h1au1 as the edge length             */
                     /* 1 : use sqrt(edgearea) as the edge length    */
                     /* 2 : use sqrt(cellarea) of common cells       */
    memset(master->u1vh, 0, geom->sze * sizeof(double));
    memset(master->u2kh, 0, geom->sze * sizeof(double));
    for (ee = 1; ee <= geom->n3_e1; ee++) {
      e = geom->w3_e1[ee];
      es = geom->m2de[e];
      c1 = geom->e2c[e][0];
      c2 = geom->e2c[e][1];
      if (etype == 0) {
	d1 = geom->h1au1[es];
	dh = master->hmean1;
      } else if (etype == 1) {
	d1 = sqrt(geom->edgearea[es]);
	dh = sqrt(master->edmean);
      } else if (etype == 2) {
	d1 = 0.5 * (sqrt(geom->cellarea[geom->m2d[c1]]) + 
		    sqrt(geom->cellarea[geom->m2d[c2]]));
	dh = sqrt(master->amean);
      }

      if (params->diff_scale & NONE)
	master->u1vh[e] = fabs(vh);
      /* Note : horizontal diffusion coeffients are scaled to the grid */
      if (params->diff_scale & LINEAR)
	master->u1vh[e] = fabs(vh * d1 / dh);
      if (params->diff_scale & NONLIN)
	master->u1vh[e] = fabs(vh * d1 * d1 / (dh * dh));
      if (params->diff_scale & CUBIC)
	master->u1vh[e] = fabs(vh * d1 * d1 * d1 / (dh * dh * dh));
      if (params->diff_scale & AREAL) {
	/*
	d1 = 0.5 * (geom->cellarea[geom->m2d[c1]] + geom->cellarea[geom->m2d[c2]]);
	master->u1vh[e] = fabs(vh * d1 / master->amean);
	*/
	master->u1vh[e] = fabs(vh * geom->edgearea[es] / master->edmean);
      }

      /* Scale the Laplacian viscosity to a biharmonic value, see    */
      /* Griffies and Hallberg (2000) Mon. Wea. Rev. 128. Section 2a */
      /* where we use edge area for del^2.                           */
      if (params->diff_scale & SCALEBI)
	master->u1vh[e] *= (0.125 * geom->edgearea[es]);
    }

    /* Cell centered horizontal viscosity diagnostic                 */
    if (master->u1vhc) {
      for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      c2 = geom->m2d[c];
      master->u1vhc[c] = 0.0;
      d1 = 0.0;
      for (ee = 1; ee <= geom->npe[c2]; ee++) {
	e = geom->c2e[ee][c];
	es = geom->m2de[e];
	dh = geom->edgearea[es];
	if (params->diff_scale & LINEAR) dh = sqrt(dh);
	if (params->diff_scale & CUBIC) dh = dh * sqrt(dh);
	master->u1vhc[c] += (dh * master->u1vh[e]);
	d1 += dh;
      }
      master->u1vhc[c] /= d1;
      if (params->diff_scale & SCALEBI)
	master->u1vhc[c] /= (0.125 * geom->cellarea[c2]);
      }
    }

    /* Smoothing                                                     */
    if (params->diff_scale & CUBIC) {
      int cs, j, eoe;
      double *u1vh = d_alloc_1d(geom->sze);


      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	d1 = sqrt(geom->cellarea[cs]);
	dh = sqrt(master->amean);
	u1vh[c] = fabs(vh * d1 * d1 / (dh * dh));
	if (params->diff_scale & SCALEBI)
	  u1vh[c] *= (0.125 * d1 * d1);
      }
      smooth3(master, u1vh, geom->w3_t, geom->b3_t, geom->sze, -1);

      for (ee = 1; ee <= geom->n3_e1; ee++) {
	e = geom->w3_e1[ee];
	c1 = geom->wgst[geom->e2c[e][0]] ? geom->e2c[e][1] : geom->e2c[e][0];
	c2 = geom->wgst[geom->e2c[e][1]] ? geom->e2c[e][0] : geom->e2c[e][1];
	master->u1vh[e] = 0.5 * (u1vh[c1] + u1vh[c2]);
      }

      for (m = 0; m < 6; m++) {
	memcpy(u1vh, master->u1vh, geom->sze * sizeof(double));
	memset(master->u1vh, 0, geom->sze * sizeof(double));
	for (ee = 1; ee <= geom->b3_e1; ee++) {
	  e = geom->w3_e1[ee];
	  es = geom->m2de[e];
	  d1 = 0.0;
	  for (j = 1; j <= geom->nee[es]; j++) {
	    eoe = geom->eSe[j][e];
	    if (eoe) {
	      master->u1vh[e] += u1vh[eoe];
	      d1 += 1.0;
	    }
	  }
	  master->u1vh[e] /= d1;
	}
      }
      d_free_1d(u1vh);
    }
  }

  /* Set the edge centered horizontal diffusivity                    */
  if (params->u1kh > 0.0 || params->bkue1 > 0.0) {
    double vh = (params->u1kh > 0.0) ? params->u1kh : params->bkue1;
    for (ee = 1; ee <= geom->n3_e1; ee++) {
      int etype = 1;
      e = geom->w3_e1[ee];
      es = geom->m2de[e];
      c1 = geom->e2c[e][0];
      c2 = geom->e2c[e][1];
      if (etype == 0) {
	d1 = geom->h1au1[es];
	dh = master->hmean1;
      } else if (etype == 1) {
	d1 = sqrt(geom->edgearea[es]);
	dh = sqrt(master->edmean);
      } else if (etype == 2) {
	d1 = 0.5 * (sqrt(geom->cellarea[geom->m2d[c1]]) + 
		    sqrt(geom->cellarea[geom->m2d[c2]]));
	dh = sqrt(master->amean);
      }
      if (params->diff_scale & NONE)
	master->u2kh[e] = fabs(vh);
      /* Note : horizontal diffusion coeffients are scaled to the grid */
      if (params->diff_scale & LINEAR)
	master->u2kh[e] = fabs(vh * d1 / dh);
      if (params->diff_scale & (NONLIN|CUBIC))
	master->u2kh[e] = fabs(vh * d1 * d1 / (dh * dh));
      if (params->diff_scale & AREAL) {
	master->u2kh[e] = fabs(vh * geom->edgearea[es] / master->edmean);
      }
    }
  }

  /* Set the ghost cells                                             */
  for (cc = 1; cc <= geom->nbpte1S; cc++) {
    c = geom->bpte1[cc];
    c2 = geom->bine1[cc];
    geom->botzu1[c] = geom->botzu1[c2];
    geom->botzu1[c] = HUGE_VAL;
  }

  /*-----------------------------------------------------------------*/
  /* Optimised horizontal mixing                                     */
  if(params->diff_scale & AUTO) {
    double hf = 0.05;             /* Factor for horizontal diffusion */
    double step = 1;              /* Integral step of diffusion > 1  */
    double hmax = 1e10;
    double d1;
    int u1khf = 0, u1vhf = 0, u2khf = 0, i1, cs;
    if (params->u1kh >= 0.0)
      u1khf = u2khf = 1;
    if (params->u1vh >= 0.0)
      u1vhf = 1;

    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      /* Note: Stability criterion for a quad is:                    */
      /* 1/[(1/h1*h1 + 1/h2*h2) * 4 * dt]. For arbitary polgons      */
      /* assume h1 = h2 = sqrt(area), then the stability criterion   */
      /* is 1/[2/area] * 4 * dt.                                     */
      hmax = 1.0 / ((2.0 / geom->cellarea[cs]) * 4.0 * params->grid_dt);
      i1 = (int)hmax / (int)step;
      hmax = step * (double)i1;
      d1 = 0.01 * 4.0 * geom->cellarea[cs] / params->grid_dt;
      i1 = (int)d1 / (int)step;
      if (u1khf)
	master->u1kh[c] = step * (double)i1;
      if (master->u1kh[c] == 0.0) master->u1kh[c] = d1;
      /* Set limits                                                  */
      if (u1khf) {
	if (master->u1kh[c] > hmax)
	  master->u1kh[c] = hmax;
	if (master->u1kh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u1kh[c] = step * (double)i1;
	}
      }
    }

    for (ee = 1; ee <= geom->n3_e2; ee++) {
      e = geom->w3_e1[ee];
      es = geom->m2de[e];
      /* Use the same stability approach as for centres above. Could */
      /* also assume h1 = h1au1 and h2 = h2au1 for edges and use the */
      /* criterion for quads.                                        */
      hmax = 1.0 / ((2.0 / geom->edgearea[es]) * 4.0 * params->grid_dt);
      i1 = (int)hmax / (int)step;
      hmax = step * (double)i1;
      d1 = 0.01 * geom->edgearea[es] / params->grid_dt;
      i1 = (int)d1 / (int)step;
      if (u1vhf)
	master->u1vh[e] = master->u2kh[e] = step * (double)i1;
      if (master->u1vh[e] == 0.0) master->u1vh[e] = d1;
      if (master->u2kh[e] == 0.0) master->u2kh[e] = d1;
      /* Set limits */
      if (u1vhf) {
	if (master->u1vh[e] > hmax)
	  master->u1vh[e] = hmax;
	if (master->u1vh[e] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u1vh[e] = step * (double)i1;
	}
	if (params->diff_scale & SCALEBI || params->visc_method & US_BIHARMONIC)
	  master->u1vh[e] *= (0.125 * geom->edgearea[es]);
      }
      if (u2khf) {
	if (master->u2kh[e] > hmax)
	  master->u2kh[e] = hmax;
	if (master->u2kh[e] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u2kh[e] = step * (double)i1;
	}
      }
    }
  }

  /* Scale horizontal mixing                                         */
  if ((d1 = get_scaling(params->scale_v, "U1VH")) != 1.0) {
    for (ee = 1; ee <= geom->n3_e1; ee++) {
      e = geom->w3_e1[ee];
      master->u1vh[e] *= d1;
    }
  }
  if ((d1 = get_scaling(params->scale_v, "U1KH")) != 1.0) {
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      master->u1kh[c] *= d1;
    }
    for (ee = 1; ee <= geom->n3_e1; ee++) {
      e = geom->w3_e1[ee];
      master->u2kh[e] *= d1;
    }
  }

  for (ee = 1; ee <= geom->n2_e1; ee++) {
    e = geom->w2_e1[ee];
    if (master->u1vhin)
      master->u1vhin[e] = master->u1vh[e];
  }

  /* Set the sponge zones if required                                */
  /*
  w1 = d_alloc_1d(geom->sze);
  set_sponge(geom, master->u1vh, params->grid_dt, w1);
  d_free_1d(w1);
  */

  /* Smooth horizontal mixing                                        */
  ns = get_smoothing(params->smooth_v, "U1VH");
  for (n = 0; n < ns; n++) {
    /*smooth3(master, master->u1vh, geom->w3_t, geom->b3_t, geom->sze, -1);*/
    smooth3e(master, master->u1vh, geom->w3_e1, geom->b3_e1, geom->sze);
  }
  ns = get_smoothing(params->smooth_v, "U1KH");
  for (n = 0; n < ns; n++) {
    smooth3(master, master->u1kh, geom->w3_t, geom->n3_t, geom->szc, 0);
    smooth3(master, master->u2kh, geom->w3_t, geom->b3_t, geom->sze, -1);
  }

  /* Set the ghost cells                                             */
  for (ee = 1; ee <= geom->nbpte1; ee++) {
    e = geom->bpte1[ee];
    c2 = geom->bine1[ee];
    master->u1vh[e] = master->u1vh[c2];
    master->u2kh[e] = master->u2kh[c2];
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for grid corners                        */
  for (vv = 1; vv <= geom->n2_e2; vv++) {
    double maxbz = -1e10;
    double sum = 0.0;
    double d1 = 0.0;
    v = geom->w2_e2[vv];

    for (cc = 1; cc <= geom->nvc[v]; cc++) {
      c = geom->v2c[v][cc];
      if (!c) continue;
      maxbz = max(maxbz, geom->botz[c]);
      if (!isnan(geom->botz[c])) {
	sum += geom->botz[c];
	d1 += 1.0;
      }
    }
    geom->botzgrid[v] = maxbz;       /* Minimum depth                */
    if (d1)
      geom->botzgrid[v] = sum / d1;  /* Average depth                */
  }

  /*-----------------------------------------------------------------*/
  /* Initialize the tidal energy extraction                          */
  if (master->nturb = params->nturb) {
    master->turb = i_alloc_1d(master->nturb);
    master->cturb = d_alloc_1d(master->nturb);
    for (m = 0; m < params->nturb; m++) {
      master->turb[m] = hd_grid_xyztoc_m(master, params->turbv[0][m], 
					 params->turbv[1][m], params->turbv[2][m]);
      master->cturb[m] = params->turbv[3][m];
      if (master->turb[m] <= 0) hd_warn("Can't find #%d tidal turbine location: [%f %f]\n",
					m, params->turbv[0][m], params->turbv[1][m]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Initialize arrays                                               */
  memset(master->waterss2d, 0, geom->szcS * sizeof(double));
  memset(master->wbot, 0, geom->szcS * sizeof(double));
  init_wvel_bounds(geom, master);        /* Initialize wbot & detadt */

  /* Get the indicator of mesh uniformity if required.               */
  /* Target edge and tangential velocities are chosen for each edge, */
  /* and these are rotated east and north. For every edge nee edges  */
  /* of eSe a normal vector is reconstructed from these east and     */
  /* north components, and the tangential velocity using wAe is      */
  /* calculated. On a uniform mesh this tangential velocity should   */
  /* equal the target tangential velocity; the uniformity is the %   */
  /* deviation from this value averaged onto the cell centre.        */
  if (params->numbers1 & MESHUN) {
    memset(master->meshun, 0, geom->szcS * sizeof(double));
    for (cc = 1; cc <= geom->b2_t; cc++) {
      int eoe;
      double ue, ve, vt, sum = 0.0;
      double u1 = 1.0;     /* Target normal velocity                 */
      double u2 = 0.0;     /* Target tangential velocity             */

      c = geom->w2_t[cc];
      for (n = 1; n <= geom->npe[c]; n++) {
	e = geom->c2e[n][c];
	vt = 0.0;
	double ux = u1 * geom->costhu1[e] - u2 * geom->sinthu1[e];
	double uy = u1 * geom->sinthu1[e] + u2 * geom->costhu1[e];
	for (ee = 1; ee <= geom->nee[e]; ee++) {
	  eoe = geom->eSe[ee][e];	  
	  ue = ux * geom->costhu1[eoe] + uy * geom->sinthu1[eoe];
	  ve = -ux * geom->sinthu2[eoe] + uy * geom->costhu2[eoe];
	  vt += ue * geom->wAe[ee][e];
	}
	master->meshun[c] += 100.0 * (vt - u2) * geom->edgearea[e];
	sum += geom->edgearea[e];
      }
      master->meshun[c] /= sum;
    }
  }
}

/* END compute_constants()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the master data structure          */
/*-------------------------------------------------------------------*/
master_t *master_build(parameters_t *params, geometry_t *geom)
{
  master_t *master;             /* Master data structure             */
  int tn,tt;

  /* Allocate memory for the master                                  */
  master = master_alloc();

  /* Transfer input data to the master                               */
  master->geom = geom;
  master->params = params;
  master->prmfd = params->prmfd;
  master->ntdif_h = params->ntdif_h;
  master->ntdif_v = params->ntdif_v;
  master->ntr = params->ntr;
  master->atr = params->atr;
  master->ntrS = params->ntrS;
  master->atrS = params->atrS;
  master->nsed = params->nsed;
  master->ntre = params->ntre;
  master->ntreS = params->ntreS;
  strcpy(master->tracerdata, params->tracerdata);
  strcpy(master->timeunit, params->timeunit);
  master->tsfile_caching = params->tsfile_caching;

  /* Constants                                                       */
  master->Cd = d_alloc_1d(geom->szcS);
  master->u1c1 = d_alloc_1d(geom->szeS);
  master->u1c3 = d_alloc_1d(geom->szeS);
  master->u1c4 = d_alloc_1d(geom->szeS);
  master->u1c5 = d_alloc_1d(geom->szeS);
  master->u1c6 = d_alloc_1d(geom->szeS);
  master->coriolis = d_alloc_1d(geom->szcS);
  master->fv = d_alloc_1d(geom->szvS);
  master->one = d_alloc_1d(geom->szm);
  for (tn = 1; tn < geom->szm; tn++)
    master->one[tn] = 1.0;

  /* Vertical grid geometry                                          */
  master->dz = d_alloc_1d(geom->szc);
  master->dzu1 = d_alloc_1d(geom->sze);
  master->depth_e1 = d_alloc_1d(geom->szeS);
  master->topz = d_alloc_1d(geom->szcS);
  master->topk = i_alloc_1d(geom->szcS);
  master->botk = i_alloc_1d(geom->szcS);

  /* SIGMA variables                                                 */
  if (params->sigma) {
    master->Ds = d_alloc_1d(geom->szcS);
    master->Hs = d_alloc_1d(geom->szcS);
    master->Hn1 = d_alloc_1d(geom->szeS);
  } else {
    master->Ds = master->one;
    master->Hs = master->one;
    master->Hn1 = master->one;
    master->Hn2 = master->one;
  }

  /* Velocity variables                                              */
  master->u1 = d_alloc_1d(geom->sze);
  master->u2 = d_alloc_1d(geom->sze);
  master->nu1 = d_alloc_1d(geom->sze);
  master->u1bot = d_alloc_1d(geom->szeS);
  master->u1flux3d = d_alloc_1d(geom->sze);
  master->u1adv = d_alloc_1d(geom->szeS);
  master->u1inter = d_alloc_1d(geom->szeS);
  master->u1av = d_alloc_1d(geom->szeS);
  master->u2av = d_alloc_1d(geom->szeS);
  master->nu1av = d_alloc_1d(geom->szeS);
  master->u1flux = d_alloc_1d(geom->szeS);
  master->u = d_alloc_1d(geom->szc);
  master->v = d_alloc_1d(geom->szc);
  master->uav = d_alloc_1d(geom->szeS);
  master->vav = d_alloc_1d(geom->szeS);

  master->d2 = d_alloc_1d(geom->szmS);
  master->d3 = d_alloc_1d(geom->szm);

  /* Vertical velocity variables */
  master->w = d_alloc_1d(geom->szc);
  master->wtop = d_alloc_1d(geom->szcS);
  master->wbot = d_alloc_1d(geom->szcS);

  /* Unstructured                                                    */
  master->circ = d_alloc_1d(geom->szv);
  master->kec = d_alloc_1d(geom->szc);
  master->div = d_alloc_1d(geom->szc);
  master->rvor = d_alloc_1d(geom->szv);
  master->nrvor = d_alloc_1d(geom->szv);
  master->npvor = d_alloc_1d(geom->szv);
  master->rvor = d_alloc_1d(geom->szv);
  master->nrvore = d_alloc_1d(geom->sze);
  master->npvore = d_alloc_1d(geom->sze);
  master->nrvorc = d_alloc_1d(geom->sgsiz);

  master->u1b = d_alloc_1d(geom->sze);
  master->u2b = d_alloc_1d(geom->sze);
  master->u1avb = d_alloc_1d(geom->szeS);
  master->etab = d_alloc_1d(geom->szcS);

  /* 3D Tracer constants and variables                               */
  if (params->ntr) {
    master->tr_wc = d_alloc_2d(geom->szc, master->ntr);
    master->trinfo_3d =
      (tracer_info_t *)malloc(sizeof(tracer_info_t) * master->ntr);
    memset(master->trinfo_3d, 0, sizeof(tracer_info_t) * master->ntr);
    for (tn = 0; tn < master->ntr; tn++) {
      tracer_copy(&master->trinfo_3d[tn], &params->trinfo_3d[tn]);
      strcpy(master->trinfo_3d[tn].name, params->trinfo_3d[tn].name);
    }

    master->mintr = d_alloc_1d(master->ntr);
    master->maxtr = d_alloc_1d(master->ntr);
    master->fill_value_wc = d_alloc_1d(master->ntr);
    master->advect = i_alloc_1d(master->ntr);
    master->diffuse = i_alloc_1d(master->ntr);

    master->trinc = i_alloc_1d(params->ntr);
    memset(master->trinc, 0 , params->ntr * sizeof(int));
  }
  if (params->ntrS) {
    master->trincS = i_alloc_1d(params->ntrS);
    memset(master->trincS, 0 , params->ntrS * sizeof(int));
  }
  if (master->ntdif_h)
    master->tdif_h = i_alloc_1d(master->ntdif_h);
  if (master->ntdif_v)
    master->tdif_v = i_alloc_1d(master->ntdif_v);
  master->trname = (char **)malloc(master->ntr * sizeof(char *));
  for (tn = 0; tn < master->ntr; tn++) {
     master->trname[tn] = (char *)malloc(sizeof(char)*MAXSTRLEN);
     memset(master->trname[tn], 0, sizeof(char)*MAXSTRLEN);
  }
  master->df_diagn_set = 0;
  master->trfilter = params->trfilter;

  /* Set the 3D tracer pointers                                      */
  if (master->ntr) {
    init_tracer_3d(params, master);
    /* Set the mean tracer to 3D auto tracers if required            */
    if (params->means & MTRA3D) {
      for (tn = 0; tn < master->ntr; tn++) {
	if (contains_token(params->means_tra, master->trinfo_3d[tn].name) != NULL) {
	  strcpy(params->means_tra, master->trinfo_3d[tn].name);
	}
      }
    }
  }

  /* 2D Tracer constants and variables                               */
  if (master->ntrS) {
    master->tr_wcS = d_alloc_2d(geom->szmS, master->ntrS);
    master->trinfo_2d =
      (tracer_info_t *)malloc(sizeof(tracer_info_t) * master->ntrS);
    memset(master->trinfo_2d, 0, sizeof(tracer_info_t) * master->ntrS);
    for (tn = 0; tn < master->ntrS; tn++) {
      tracer_copy(&master->trinfo_2d[tn], &params->trinfo_2d[tn]);
      strcpy(master->trinfo_2d[tn].name, params->trinfo_2d[tn].name);
    }
    init_tracer_2d(params, master);
    /* Set the mean tracer to 2D auto tracers if required            */
    if (params->means & MTRA2D) {
      for (tn = 0; tn < master->ntrS; tn++) {
	if (contains_token(params->means_tra, master->trinfo_2d[tn].name) != NULL) {
	  strcpy(params->means_tra, master->trinfo_2d[tn].name);
	}
      }
    }
  }

  /* Sediment tracer constants and variables                         */
  if (master->nsed) {
    master->tr_sed = d_alloc_3d(geom->szcS, geom->sednz, master->nsed);
    master->trinfo_sed =
      (tracer_info_t *)malloc(sizeof(tracer_info_t) * master->nsed);
    memset(master->trinfo_sed, 0, sizeof(tracer_info_t) * master->nsed);
    for (tn = 0; tn < master->nsed; tn++) {
      tracer_copy(&master->trinfo_sed[tn], &params->trinfo_sed[tn]);
      strcpy(master->trinfo_sed[tn].name, params->trinfo_sed[tn].name);
    }
    init_tracer_sed(params, master);
  }

  /* Do any custom initialisation of the sediments tracer            */
#if defined(HAVE_SEDIMENT_MODULE)
  sed_tracer_custom(master);
#endif

  /* Mixing constants and variables                                  */
  master->Kz = d_alloc_1d(geom->szc);
  master->Vz = d_alloc_1d(geom->szc);
  master->z0 = d_alloc_1d(geom->szcS);

  /* Surface elevation variables                                     */
  master->eta = d_alloc_1d(geom->szcS);
  master->oldeta = d_alloc_1d(geom->szcS);
  master->waterss = d_alloc_1d(geom->szc);
  master->waterss2d = d_alloc_1d(geom->szcS);
  master->detadt = d_alloc_1d(geom->szcS);
  master->wdiff2d = d_alloc_1d(geom->szcS);

  /* Density variables                                               */
  master->dens = d_alloc_1d(geom->szc);
  master->dens_0 = d_alloc_1d(geom->szc);
  master->topdensu1 = d_alloc_1d(geom->szeS);
  master->densavu1 = d_alloc_1d(geom->szeS);

  /* Atmospheric variables                                           */
  master->wind1 = d_alloc_1d(geom->szeS);
  master->windspeed = d_alloc_1d(geom->szeS);
  master->winddir = d_alloc_1d(geom->szeS);
  if (params->wind_dt && params->storm_dt) {
    master->swind1 = d_alloc_1d(geom->szeS);
  }
  master->patm = d_alloc_1d(geom->szcS);

  /* Heatflux variables                                              */
  if (params->heatflux & (NET_HEAT | ADVANCED | INVERSE | COMP_HEAT | COMP_HEAT_MOM)) {
    master->heatf = d_alloc_1d(geom->szcS);
    if (!master->swr)
      master->swr = d_alloc_1d(geom->szcS);
  }

  /* Horizontal diffusion variables                                  */
  /* smagcode = U1_A  : No Smagorinsky; u1vh allocated               */
  /* smagcode = U1_AK : No Smagorinsky; u1kh, u2kh allocated         */
  /* smagcode = U1_SP : Smagorinsky; u1vh = sdc, no sponges          */
  /* smagcode = U1_SA : Smagorinsky; u1vh allocated, uses sponges    */
  /* smagcode = U1_SPK: Smagorinsky; u1kh = sdc, no sponges          */
  /* smagcode = U1_SAK: Smagorinsky; u1kh allocated, uses sponges    */
  master->smagcode = 0;
  /* Uncomment this code if separate sponges apply in the e1 and e2  */
  /* directions.                                                     */

  if (params->smagorinsky > 0.0 && params->u1vh < 0.0) {
    int n;
    if (params->smagorinsky == 1.0 && params->sue1 != 1.0) 
      master->smagcode |= U1_SA;
    for (n = 0; n < params->nobc; ++n)
      if (params->open[n]->sponge_zone_h)
	master->smagcode |= U1_SA;
    if (master->smagcode & U1_SA) {
      master->u1vh = d_alloc_1d(geom->sze);
      master->sde = d_alloc_1d(geom->sze);
    } else {
      /*master->u1vh = d_alloc_1d(geom->sze);*/
      master->sde = d_alloc_1d(geom->sze);
      master->u1vh = master->sde;
      master->smagcode |= U1_SP;
    }
  } else {
    master->u1vh = d_alloc_1d(geom->sze);
    master->smagcode |= U1_A;
  }
  if (params->smagorinsky > 0.0 && params->u1kh < 0.0) {
    int n;
    tn = tracer_find_index("smagorinsky", master->ntr, master->trinfo_3d);
    if (params->smagorinsky == 1.0 && params->kue1 != 1.0) 
      master->smagcode |= U1_SAK;
    for (n = 0; n < params->nobc; ++n) {
      if (params->open[n]->sponge_zone_h)
	master->smagcode |= U1_SAK;
      if (params->open[n]->bcond_tra[tn] == FILEIN)
      	master->smagcode |= U1_SBC;
    }
    if (master->smagcode & U1_SAK) {
      master->u1kh = d_alloc_1d(geom->szc);
      if (!master->sde) master->sde = d_alloc_1d(geom->sze);
      master->u2kh = d_alloc_1d(geom->sze);
    } else {
      master->u1kh = master->sdc;
      if (!master->sde) master->sde = d_alloc_1d(geom->sze);
      master->u2kh = master->sde;
      master->smagcode |= U1_SPK;
    }
  } else {
    master->u1kh = d_alloc_1d(geom->szc);
    master->u2kh = d_alloc_1d(geom->sze);
    master->smagcode |= U1_AK;
    if (params->sigma)
      hd_warn("** Smagorinsky e1 diffusion works best with sigma. **\n");
  }
  master->t11 = d_alloc_1d(geom->szm);
  master->t22 = d_alloc_1d(geom->szm);
  master->t12 = d_alloc_1d(geom->szm);

  /* Miscillaneous */
  if (params->thin_merge) {
    master->kth_e1 = i_alloc_1d(geom->szeS);
  }
  if (!(params->means & NONE)) {
    master->meanc = d_alloc_1d(geom->szcS);
    memset(master->meanc, 0, geom->szcS * sizeof(double));
    if (params->means & TIDAL) {
      master->odeta = d_alloc_1d(geom->szcS);
      memset(master->odeta, 0, geom->szcS * sizeof(double));
    }
  }
  if (params->runmode & TRANS) {
    /* Use the viscosity tensors to store streamline Courant nos.    */
    master->origin = master->Vz;
    master->pc = master->t11;
    master->qc = master->t12;
    master->rc = master->t22;
  }

  master->wclk = d_alloc_1d(geom->nwindows + 1);
  memset(master->wclk, 0, geom->nwindows * sizeof(double));

  /* Initialize                                                      */
  memset(master->wtop, 0, geom->szcS * sizeof(double));
  memset(master->wbot, 0, geom->szcS * sizeof(double));
  memset(master->waterss, 0, geom->szc * sizeof(double));
  memset(master->waterss2d, 0, geom->szcS * sizeof(double));
  memset(master->eta, 0, geom->szcS * sizeof(double));
  memset(master->detadt, 0, geom->szcS * sizeof(double));
  memset(master->u1av, 0, geom->szeS * sizeof(double));
  memset(master->w, 0, geom->szc * sizeof(double));
  memset(master->u1, 0, geom->sze * sizeof(double));

  /* Initialise s2m_3d tracer maps                                   */
  if (master->ntr) {
    master->trmap_s2m_3d = i_alloc_1d(master->ntr);
    tt = 0;
    for (tn = 0; tn < master->ntr; tn++) {
      if (strncmp(master->trinfo_3d[tn].tag, "DA_", 3)) {
	master->trmap_s2m_3d[tt] = tn;
	tt++;
      }
    }
    master->ntrmap_s2m_3d = tt;
  }

#ifdef HAVE_OMP
  master->trans_num_omp = params->trans_num_omp;
#endif

  if (DEBUG("init_w"))
    dlog("init_w", "Master data initialised OK\n");
  return (master);
}

/* END master_build()                                                */
/*-------------------------------------------------------------------*/


/*UR-CHANGED renamed */
/*-------------------------------------------------------------------*/
/* Frees memory from arrays which reside in the master               */
/*-------------------------------------------------------------------*/
void master_free_nwin(master_t *master)
{
  d_free_1d(master->dz);
  d_free_1d(master->q);
  d_free_1d(master->t11);
  d_free_1d(master->t12);
  d_free_1d(master->t22);
  d_free_1d(master->oldeta);
  d_free_1d(master->topdensu1);
  d_free_1d(master->topdensu2);
  d_free_1d(master->densavu1);
  d_free_1d(master->densavu2);
  d_free_1d(master->u1adv);
  d_free_1d(master->u2adv);
  d_free_1d(master->u1inter);
  d_free_1d(master->u2inter);
  i_free_1d(master->botk);
  i_free_1d(master->topk);
  d_free_1d(master->topz);
  d_free_1d(master->Cd);
  d_free_1d(master->z0);
  d_free_1d(master->patm);
  d_free_1d(master->wind1);
  d_free_1d(master->wind2);
  d_free_1d(master->windspeed);
  d_free_1d(master->winddir);
  if (master->wind_dt && master->storm_dt) {
    d_free_1d(master->swind1);
    d_free_1d(master->swind2);
  }
  d_free_1d(master->waterss);
  d_free_1d(master->waterss2d);
  d_free_1d(master->u1c1);
  d_free_1d(master->u1c3);
  d_free_1d(master->u1c4);
  d_free_1d(master->u1c5);
  d_free_1d(master->u1c6);
  d_free_1d(master->u1vh);
  d_free_1d(master->u2c1);
  d_free_1d(master->u2c3);
  d_free_1d(master->u2c4);
  d_free_1d(master->u2c5);
  d_free_1d(master->u2c6);
  d_free_1d(master->u2vh);
  if (master->sigma) {
    d_free_1d(master->Ds);
    d_free_1d(master->Hs);
    d_free_1d(master->Hn1);
    d_free_1d(master->Hn2);
  }
  d_free_1d(master->one);
  if (!(master->means & NONE)) {
    d_free_1d(master->meanc);
    d_free_1d(master->meancs);
    if (master->means & TIDAL)
      d_free_1d(master->odeta);
  }
  if (master->ef)
    d_free_1d(master->ef);
  if (master->d2) d_free_1d(master->d2);
  if (master->d3) d_free_1d(master->d3);
  if (master->trtot) d_free_1d(master->trtot);
  if (master->trinc) i_free_1d(master->trinc);
  if (master->trincS) i_free_1d(master->trincS);
  if (master->vf)
    d_free_1d(master->vf);
  if (master->ntrmap_s2m_3d)
    i_free_1d(master->trmap_s2m_3d);
}

/* END master_free_nwin()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the master data structure          */
/*-------------------------------------------------------------------*/
master_t *master_alloc(void)
{
  master_t *master = (master_t *)malloc(sizeof(master_t));
  memset(master, 0, sizeof(master_t));
  return master;
}

/* END master_alloc()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the parameter data structure       */
/*-------------------------------------------------------------------*/
parameters_t *params_alloc(void)
{
  parameters_t *params = (parameters_t *)malloc(sizeof(parameters_t));
  memset(params, 0, sizeof(parameters_t));
  return params;
}

/* END params_alloc()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the OBC data structure             */
/*-------------------------------------------------------------------*/
open_bdrys_t *OBC_alloc(void)
{
  open_bdrys_t *OBC = (open_bdrys_t *)malloc(sizeof(open_bdrys_t));
  memset(OBC, 0, sizeof(open_bdrys_t));
  return OBC;
}

/* END OBC_alloc()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate tracers from Cartesian to sigma            */
/*-------------------------------------------------------------------*/
void sigma_tracer_step_3d(master_t *master)
{
  geometry_t *geom = master->geom;
  int c, cc;
  int cs, c2, ci;
  int k, kk, n;
  int nz = geom->nz - 1;
  double tr[geom->nz + 1];
  double zs;

  for (n = 0; n < master->ntr; n++) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = c2 = geom->w2_t[cc];
      cs = geom->m2d[c];

      /* Copy the tracer values into a buffer                        */
      kk = nz;
      while (c != geom->zm1[c]) {
        tr[kk] = master->tr_wc[n][c];
        /* Set values below the 'z' bottom to the last value         */
        if ((!tr[kk] || isnan(tr[kk])) && kk < nz)
          tr[kk] = tr[kk + 1];
        kk--;
        c = geom->zm1[c];
      }

      c = c2;
      while (c != geom->zm1[c]) {
        zs = geom->cellz[c] * master->Ds[cs];

        /* Find the z layers bracketing the sigma layer.             */
        /* Note : master->dz contains the original cellz depths as   */
        /* read from the input file.                                 */
        ci = c2;
        kk = nz;
        while (master->dz[ci] > zs && ci != geom->zm1[ci]) {
          ci = geom->zm1[ci];
          kk--;
        }

        /* Do the interpolation                                      */
        if (ci == geom->zp1[ci])
          master->tr_wc[n][c] = tr[nz];
        else if (ci == geom->zm1[ci])
          master->tr_wc[n][c] = tr[0];
        else {
          /* Linear interpolation                                    */
          master->tr_wc[n][c] =
            (tr[kk] - tr[kk + 1]) * (zs -
                                     master->dz[ci]) / (master->dz[ci] -
                                                        master->dz[geom->
                                                                   zp1
                                                                   [ci]]) +
            tr[kk];
          /* 2nd order extrapolation */
          /*
             kp1=geom->zp1[ci]; kp2=geom->zp1[kp1]; kp3=geom->zp1[kp2];
             km1=geom->zm1[kk]; km2=geom->zm1[km1];
             master->tr_wc[n][c]=0.5*(bc_polint(tr[kp3],tr[kp2],tr[kp1])+
             bc_polint(tr[km2],tr[km1],tr[kk])); */
        }
        c = geom->zm1[c];
      }

      /* Apply a smoothing filter to the vertical profile            */
      master->tr_wc[n][c] = master->tr_wc[n][geom->zp1[c]];
      c = c2;
      while (c != geom->zm1[c]) {
        k = geom->s2k[c];
        tr[k] = 0.25 * (master->tr_wc[n][geom->zp1[c]] +
                        master->tr_wc[n][geom->zm1[c]]) +
          0.5 * master->tr_wc[n][c];
        c = geom->zm1[c];
      }

      c = c2;
      while (c != geom->zm1[c]) {
        k = geom->s2k[c];
        master->tr_wc[n][c] = tr[k];
        c = geom->zm1[c];
      }
    }
  }
}

/* END sigma_tracer_step_3d()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the surface and bottom vertical velocity     */
/*-------------------------------------------------------------------*/
void init_wvel_bounds(geometry_t *geom, /* Global geometry           */
                      master_t *master  /* Master data               */
  )
{
  int c, cc, c1, c2, n;         /* Sparse coodinate / counter        */
  int cs;                       /* 2D sparse coordinate              */
  int cb;                       /* Bottom sparse coordinate          */
  double eta_l;                 /* Surface elevation at e2c[0]       */
  double eta_r;                 /* Surface elevation at e2c[1]       */

  int e, ee, es, eb;
  double *u1bot;                /* Bottom e1 velocity                */


  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  u1bot = d_alloc_1d(geom->szeS);
  memset(u1bot, 0, geom->szeS * sizeof(double));
  for (ee = 1; ee <= geom->b2_e1; ee++) {
    eb = geom->bot_e1[ee];    /* 3D bottom coordinate                */
    es = geom->m2de[eb];
    u1bot[es] = master->u1[eb];
  }
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                     geom->bin2, 0);

  /*-----------------------------------------------------------------*/
  /* SIGMA : zero velocity at the sigma boundaries                   */
  if (master->sigma)
    return;

  /* In the linear case wtop and detadt are the same, whereas in the */
  /* non-linear case, they are related by terms involving the        */
  /* surface slope and surface horizontal velocities.                */
  memcpy(master->wtop, master->detadt, geom->enonS * sizeof(double));
  if (!master->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.                                   */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    cb = geom->bot_t[cc];
    cs = geom->m2d[c];

    master->detadt[cs] = 0.0;
    master->wbot[cs] = 0.0;
    for (n = 1; n <= geom->npe[c] / 2; n++) {
      e = geom->c2e[n][c];
      es = geom->m2de[e];
      c1 = geom->e2c[e][0];
      c2 = geom->e2c[e][1];
      /* Calculate the surface vertical velocity                     */
      eta_l = master->eta[c1];
      eta_r = master->eta[c2];
      master->detadt[cs] += (master->u1[es] * (eta_r - eta_l) / 
			     (geom->h2au1[es] * (double)geom->npe[cs]));

      /* Calculate the bottom vertical velocity                      */
      master->wbot[cs] -= (master->u1bot[es] * geom->dHde1[es] / 
			   (geom->h2au1[es]* (double)geom->npe[cs]));
    }
  }
  d_free_1d(u1bot);
}

/* END init_wvel_bounds()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to re-initialize the window partitioning                  */
/*-------------------------------------------------------------------*/
void reset_windows(hd_data_t *hd_data)
{
  geometry_t *geom = hd_data->geom;
  master_t *master = hd_data->master;
  parameters_t *params = hd_data->params;
  int n;
  double d1 = 0.0, d2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the new window sizes                                        */
  if (geom->nwindows > 1 && geom->win_size) {
    for (n = 1; n <= geom->nwindows; n++)
      d1 += master->wclk[n];
    for (n = 1; n <= geom->nwindows; n++) {
      geom->win_size[n] *= (d1 / (master->wclk[n] * geom->nwindows));
      d2 += geom->win_size[n];
    }
    for (n = 1; n <= geom->nwindows; n++) {
      geom->win_size[n] /= d2;
      master->wclk[n] = 0.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Close existing windows                                          */
  windows_clear(hd_data);

  /*-----------------------------------------------------------------*/
  /* Set up the windows                                              */
  master->nwindows = geom->nwindows = params->nwindows;
  window_build(geom, params);

  /* Initialze the boundary conditions                               */
  bdry_custom_w(geom, window);

  /* Initialise the geometric arrays in the windows with master      */
  /* geometry data.                                                  */
  window_init(geom, window);

  /* Initialise any tidal forcing                                    */
  csr_tide_init(master, window);

  /* Make and initialise the window data structures                  */
  windat = win_data_build(master, window);

  /* Make and initialise the window private data structures          */
  wincon = win_consts_init(master, window);

  /* Calculate required initial conditions                           */
  pre_run_setup(master, window, windat, wincon);

  /* Initialise the source/sink variables                            */
  sourcesink_init(params, master, window, windat, wincon);

  /* Initialise the distributed processing                           */
  dp_init(master, window, windat, wincon, geom->nwindows);
  hd_data->window = window;
  hd_data->windat = windat;
  hd_data->wincon = wincon;
}

/* END reset_windows()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the fetch at 8 compass points                */
/*                                                                   */
/*                              fetch[0]                             */
/*                                 |                                 */
/*                              7  |  1                              */
/*                                 |                                 */
/*                        6 ------ o ------ 2                        */
/*                                 |                                 */
/*                              5  |  3                              */
/*                                 |                                 */
/*                                 4                                 */
/*                                                                   */
/*-------------------------------------------------------------------*/
void init_fetch(master_t *master, parameters_t *params, unsigned long **flag)
{
  geometry_t *geom = master->geom;
  int nn, maxnn = 5000;
  int ns = 4;
  int db = 0;
  int c, cc, c1, n;
  int ii, jj, x, y;
  int nz = geom->nz - 1;
  /* 0 1 2 3 4 5 6 7 : Mask to set fetch = 0 if mask=0 */
  int mask[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  double **fdata = NULL;
  char buf[MAXSTRLEN];
  char *fname[8] = {"fetch_n","fetch_ne","fetch_e","fetch_se",
		    "fetch_s","fetch_sw","fetch_w","fetch_nw"};

  not_included("wave fetch");
  return;

  master->fetch = d_alloc_2d(8, geom->sgsizS);

  /* Read the fetch data from file if required, including OBCs       */
  if (strlen(params->fetch_init)) {
    timeseries_t *ts;
    int id;

    /* hd_ts_read() is used to read eta from file and this requires  */
    /* time variables to be set on the master. These are only set in */
    /* compute_constants(), called after this routine, so preset     */
    /* them here.                                                    */
    master->t = params->t;
    strcpy(master->timeunit, params->timeunit);
    strcpy(master->output_tunit, params->output_tunit);

    fdata = d_alloc_2d(geom->sgsizS, 8);
    ts = hd_ts_read(master, params->fetch_init, 0);
    for (ii = 0; ii < 8; ii++) {
      memset(fdata[ii], 0, geom->sgsizS * sizeof(double));
      id = ts_get_index(ts, fv_get_varname(params->fetch_init, fname[ii], buf));
      if (id < 0) {
	hd_warn("init_fetch: The file '%s' does not contain the tracer '%s'.\n", params->fetch_init, fname[ii]);
	continue;
      }

      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	fdata[ii][c] = ts_eval_xy(ts, id, master->t, geom->cellx[c],
			      geom->celly[c]) * 1e3;
      }
    }
    hd_ts_free(master, ts);
  }

  /* Loop over the water cells                                       */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (n = 0; n < 8; n++)
      master->fetch[c][n] = 0.0;
    if(db)printf("start %d\n",c);
    /* North                                                         */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[0])
        master->fetch[c][0] += geom->h2acell[c1];
      ii += sign(geom->sinthcell[c1]);
      jj += sign(geom->costhcell[c1]);
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[0] && fdata)
      master->fetch[c][0] += fdata[0][geom->map[nz][jj][ii]];
    if(db)printf("n OK\n");
    /* North-east */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[1])
        master->fetch[c][1] += sqrt(geom->h1acell[c1] *
                                    geom->h1acell[c1] +
                                    geom->h2acell[c1] * geom->h2acell[c1]);
      x = (sign(geom->sinthcell[c1]) + sign(geom->costhcell[c1]));
      y = (sign(geom->costhcell[c1]) + sign(-geom->sinthcell[c1]));
      if (x > 1)
        x = 1;
      if (x < -1)
        x = -1;
      if (y > 1)
        y = 1;
      if (y < -1)
        y = -1;
      ii += x;
      jj += y;
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[1] && fdata)
      master->fetch[c][1] += fdata[1][geom->map[nz][jj][ii]];
    if(db)printf("ne OK\n");
    /* East */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[2])
        master->fetch[c][2] += geom->h1acell[c1];
      ii += sign(geom->costhcell[c1]);
      jj += sign(-geom->sinthcell[c1]);
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[2] && fdata)
      master->fetch[c][2] += fdata[2][geom->map[nz][jj][ii]];
    if(db)printf("e OK\n");
    /* South-east */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[3])
        master->fetch[c][3] += sqrt(geom->h1acell[c1] *
                                    geom->h1acell[c1] +
                                    geom->h2acell[c1] * geom->h2acell[c1]);
      x = (sign(geom->sinthcell[c1]) - sign(geom->costhcell[c1]));
      y = (sign(geom->costhcell[c1]) - sign(-geom->sinthcell[c1]));

      if (x > 1)
        x = 1;
      if (x < -1)
        x = -1;
      if (y > 1)
        y = 1;
      if (y < -1)
        y = -1;
      ii += x;
      jj += y;
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[3] && fdata)
      master->fetch[c][3] += fdata[3][geom->map[nz][jj][ii]];
    if(db)printf("se OK\n");
    /* South */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    c1 = geom->map[nz][jj][ii];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      if (mask[4])
        master->fetch[c][4] += geom->h2acell[c1];
      ii -= sign(geom->sinthcell[c1]);
      jj -= sign(geom->costhcell[c1]);
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[4] && fdata)
      master->fetch[c][4] += fdata[4][geom->map[nz][jj][ii]];
    if(db)printf("s OK\n");
    /* South-west */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn++;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[5])
        master->fetch[c][5] += sqrt(geom->h1acell[c1] *
                                    geom->h1acell[c1] +
                                    geom->h2acell[c1] * geom->h2acell[c1]);
      x = (-sign(geom->sinthcell[c1]) - sign(geom->costhcell[c1]));
      y = (-sign(geom->costhcell[c1]) - sign(-geom->sinthcell[c1]));
      if (x > 1)
        x = 1;
      if (x < -1)
        x = -1;
      if (y > 1)
        y = 1;
      if (y < -1)
        y = -1;
      ii += x;
      jj += y;
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[5] && fdata)
      master->fetch[c][5] += fdata[5][geom->map[nz][jj][ii]];
    if(db)printf("sw OK\n");
    /* West */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[6])
        master->fetch[c][6] += geom->h1acell[c1];
      ii -= sign(geom->costhcell[c1]);
      jj -= sign(-geom->sinthcell[c1]);
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[6] && fdata)
      master->fetch[c][6] += fdata[6][geom->map[nz][jj][ii]];
    if(db)printf("w OK\n");
    /* North-west */
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    nn = 0;
    while (!(flag[jj][ii] & (SOLID | OUTSIDE | U1BDRY | U2BDRY)) &&
           ii < geom->nce1 - 1 && jj < geom->nce2 - 1 && nn < maxnn) {
      nn++;
      c1 = geom->map[nz][jj][ii];
      if (mask[7])
        master->fetch[c][7] += sqrt(geom->h1acell[c1] *
                                    geom->h1acell[c1] +
                                    geom->h2acell[c1] * geom->h2acell[c1]);
      x = (sign(geom->sinthcell[c1]) - sign(geom->costhcell[c1]));
      y = (sign(geom->costhcell[c1]) - sign(-geom->sinthcell[c1]));
      if (x > 1)
        x = 1;
      if (x < -1)
        x = -1;
      if (y > 1)
        y = 1;
      if (y < -1)
        y = -1;
      ii += x;
      jj += y;
    }
    if (flag[jj][ii] & (U1BDRY | U2BDRY) && mask[7] && fdata)
      master->fetch[c][7] += fdata[7][geom->map[nz][jj][ii]];
    if(db)printf("nw OK\n");
  }
  /* Smooth the fetch                                                */
  if (fdata == NULL) fdata = d_alloc_2d(geom->sgsizS, 8);
  for (n = 0; n < 8; n++) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      fdata[n][c] = master->fetch[c][n];
    }
    for (ii = 0; ii < ns; ii++)
      smooth(master, fdata[n], geom->w2_t, geom->b2_t);
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->fetch[c][n] = fdata[n][c];
    }
  }
  if (fdata) d_free_2d(fdata);    
}

/* END init_fetch()                                                  */
/*-------------------------------------------------------------------*/

int sign(double x)
{
  return ((x < 0.5 && x > -0.5) ? 0 : ((x > 0.5) ? 1 : -1));
}


/*-------------------------------------------------------------------*/
/* Routine to return the fetch given a wind direction                */
/*-------------------------------------------------------------------*/
double get_fetch(double *fetch, double dir)
{

  dir =
    (dir < 0.0) ? (dir * 180 / 3.1415926) + 360.0 : dir * 180 / 3.1415926;

  if (dir < 30 || dir >= 330)
    return (fetch[4]);
  else if (dir >= 30 && dir < 60)
    return (fetch[5]);
  else if (dir >= 60 && dir < 120)
    return (fetch[6]);
  else if (dir >= 120 && dir < 150)
    return (fetch[7]);
  else if (dir >= 150 && dir < 210)
    return (fetch[0]);
  else if (dir >= 210 && dir < 240)
    return (fetch[1]);
  else if (dir >= 240 && dir < 300)
    return (fetch[2]);
  else if (dir >= 300 && dir < 330)
    return (fetch[3]);
  else
    return (0.0);
}

/* END get_fetch()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialise the variables required to read velocity data in the    */
/* transport mode. Note: in the structured model these data may be   */
/* read from a different grid to that which the transport model is   */
/* run on. While this is unlikely to be necessary on unstructured    */
/* we retain the structures necessary to do so if this is required   */
/* in the future.                                                    */
/*-------------------------------------------------------------------*/
void trans_sourcegrid_init(parameters_t *params, 
			   master_t *master, 
			   geometry_t *geom,
			   geometry_t **window, 
			   window_t **windat, 
			   win_priv_t **wincon,
			   dump_data_t *dumpdata)
{
  geometry_t *tpg;
  window_t *tpd;
  win_priv_t *tpc;
  transport_t *tp;
  int m, n;

  if (!(master->tmode & SP_ORIGIN) && geom->nwindows > 1)
    hd_quit("Transport mode using different source and target grids must have WINDOWS = 1.\n");

  geom->trans = (transport_t *)malloc(sizeof(transport_t));
  memset(geom->trans, 0, sizeof(transport_t));
  tp = geom->trans;

  tp->tpg = (geometry_t *)malloc(sizeof(geometry_t));
  memset(tp->tpg, 0, sizeof(geometry_t));
  tpg = tp->tpg;

  tp->tpd = (window_t *)malloc(sizeof(window_t));
  memset(tp->tpd, 0, sizeof(window_t));
  tpd = tp->tpd;

  tp->tpc = (win_priv_t *)malloc(sizeof(win_priv_t));
  memset(tp->tpc, 0, sizeof(win_priv_t));
  tpc = tp->tpc;

  /* Set defaults                                                    */
  tpg->windat = tpd;
  tpg->wincon = tpc;
  tpc->nonlinear = master->nonlinear;
  tpc->velmax = master->velmax;
  tpc->etamax = master->etamax;
  tpc->sigma = master->sigma;
  tpc->slip = master->slip;

  /*-----------------------------------------------------------------*/
  /* Check if open boundaries are specified. Used for MONOTONIC      */
  /* global fills.                                                   */
  if (master->fillf & MONOTONIC) {
    m = 0;
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (!(open->bcond_nor & NOTHIN)) m = 1;
      if (!(open->bcond_ele & NOTHIN)) m = 2;
      if (!(open->bcond_nor & (FILEIN|CUSTOM|NOTHIN)))
	hd_warn("MONOTONIC global fill: FILEIN or CUSTOM OBCs required for OBC %s - errors possible.\n", open->name); 
    }
    if (m == 1) {
      master->fillf |= SET_BDRY;
      for (n = 1; n <= geom->nwindows; n++)
	wincon[n]->fillf |= SET_BDRY;
	hd_warn("MONOTONIC global fill: boundary fluxes derived via OBCs.\n");
    }
    if (m == 2) {
      master->fillf |= SET_BDRY_ETA;
      for (n = 1; n <= geom->nwindows; n++)
	wincon[n]->fillf |= SET_BDRY_ETA;
	hd_warn("Surface elevation boundaries derived via OBCs.\n");
    }
  }

  /*-----------------------------------------------------------------*/
  /* If source and target grids are the same, point the transport    */
  /* variables at the relevant hydro variables.                      */
  if (master->tmode & SP_ORIGIN && (strcmp(params->idumpname, params->sourcefile) != 0)) {
    hd_warn("Source and target grid must be the same for STREAMLINE mode. Setting source = target.\n");
    strcpy(params->sourcefile, params->idumpname);
  }
  set_tp(geom, master, window);
  return;
}

/* END trans_sourcegrid_init()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the transport source grid                          */
/*-------------------------------------------------------------------*/
void set_tp(geometry_t *geom, master_t *master, geometry_t **window)
{
  int n, m;
  transport_t *tp = geom->trans;
  geometry_t *tpg = tp->tpg;
  window_t *tpd = tp->tpd;
  win_priv_t *tpc = tp->tpc;
  int do_multi = 0;

  tpg->nce1 = geom->nce1;
  tpg->nce2 = geom->nce2;
  tpg->nz = geom->nz;
  tpg->a2_t = geom->a2_t;
  tpg->a3_t = geom->a3_t;
  tpg->cellx = geom->cellx;
  tpg->celly = geom->celly;
  tpg->cellz = geom->cellz;
  tpg->gridz = geom->gridz;
  tpg->h1au1 = geom->h1au1;
  tpg->h2au2 = geom->h2au2;
  tpg->h1acell = geom->h1acell;
  tpg->h2acell = geom->h2acell;
  tpg->zp1 = geom->zp1;
  tpg->zm1 = geom->zm1;
  tpg->m2d = geom->m2d;
  tpg->nbpt = geom->nbpt;
  tpg->bpt = geom->bpt;
  tpg->bin = geom->bin;
  tpg->wsa = geom->wsa;
  tp->ns2 = tpg->ns2 = geom->ns2;
  tp->ns3 = tpg->ns3 = geom->ns3;
  tp->w1 = dumpdata->w1;
  if (do_multi && master->tmode & (NONE|SP_ORIGIN)) {
    tpg->sgsiz = geom->sgsiz;
    tpg->sgsizS = geom->sgsizS;
    tpg->b2_t = geom->b2_t;
    tpg->b3_t = geom->b3_t;
    tpg->n2_t = geom->n2_t;
    tpg->n3_t = geom->n3_t;
    tpg->v2_t = geom->v2_t;
    tpg->v3_t = geom->v3_t;
    tpg->w2_t = geom->w2_t;
    tpg->w3_t = geom->w3_t;
    tpg->s2i = geom->s2i;
    tpg->s2j = geom->s2j;
    tpg->s2k = geom->s2k;
    tpg->sur_t = geom->sur_t;
    tpg->nsur_t = i_alloc_1d(geom->b2_t + 1);
    tpg->bot_t = geom->bot_t;
    tpg->dHde1 = geom->dHde1;
    tpg->dHde2 = geom->dHde2;
    tpg->botz = geom->botz;
    
    tpd->u1 = master->u1;
    tpd->u2 = master->u2;
    tpd->w = master->w;
    tpd->wtop = master->wtop;
    tpd->wbot = master->wbot;
    tpd->Kz = master->Kz;
    tpd->eta = master->eta;
    tpd->etab = master->etab;
    tpd->light = master->light;
    tpd->detadt = master->detadt;
    tpd->temp = master->temp;
    tpd->sal = master->sal;
    tpd->tr_wc = (double **)malloc(sizeof(double *) * (2+master->ntrvars));
    m = tracer_find_index("salt", master->ntr, master->trinfo_3d);
    tpd->tr_wc[0] = master->tr_wc[m];
    m = tracer_find_index("temp", master->ntr, master->trinfo_3d);
    tpd->tr_wc[1] = master->tr_wc[m];
    for (m = 0; m < master->ntrvars; m++)
      tpd->tr_wc[2+m] = master->tr_wc[master->trvm[m]];
    if (master->ntrvarsS) {
      tpd->tr_wcS = (double **)malloc(sizeof(double *) * (master->ntrvarsS));
      for (m = 0; m < master->ntrvarsS; m++)
	tpd->tr_wcS[m] = master->tr_wcS[master->trvmS[m]];
    }
    for (n = 1; n <= geom->nwindows; n++) {
      win_priv_t *wincon = window[n]->wincon;
      window[n]->trans = geom->trans;
      tpc->d1 = wincon->d1;
      tpc->d2 = wincon->d2;
      tpc->d3 = wincon->d3;
      tpc->d4 = wincon->d4;
      tpc->w4 = wincon->w4;
      tpc->w5 = wincon->w5;
      tpc->w6 = wincon->w6;
      tpc->w7 = wincon->w7;
      tpc->w8 = wincon->w8;
      tpc->w9 = wincon->w9;
      tpc->w10 = wincon->w10;
      tpc->s1 = wincon->s1;
      tpc->s2 = wincon->s2;
      tpc->s3 = wincon->s3;
      tpc->i1 = wincon->i1;
      tpc->i2 = wincon->i2;
      tpc->i3 = wincon->i3;
      tpc->dz = wincon->dz;
      tpc->m2d = wincon->m2d;
      tpc->oldeta = wincon->oldeta;
      tpc->wgt = wincon->wgt;
      tpc->cdry_e1 = wincon->cdry_e1;
      tpc->tgrid = wincon->tgrid = NONE;
    }
  } else {
    tpd->tr_wc = (double **)malloc(sizeof(double *) * (2+master->ntrvars));
    if (master->ntrvarsS) {
      tpd->tr_wcS = (double **)malloc(sizeof(double *) * (master->ntrvarsS));
      for (m = 0; m < master->ntrvarsS; m++)
	tpd->tr_wcS[m] = master->tr_wcS[master->trvmS[m]];
    }
    for (n = 1; n <= geom->nwindows; n++) {
      window_t *windat = window[n]->windat;
      win_priv_t *wincon = window[n]->wincon;
      tpg->sgsiz = window[n]->sgsiz;
      tpg->sgsizS = window[n]->sgsizS;
      tpg->a2_t = window[n]->a2_t;
      tpg->a3_t = window[n]->a3_t;
      tpg->b2_t = window[n]->b2_t;
      tpg->b3_t = window[n]->b3_t;
      tpg->n2_t = window[n]->n2_t;
      tpg->n3_t = window[n]->n3_t;
      tpg->v2_t = window[n]->v2_t;
      tpg->v3_t = window[n]->v3_t;
      tpg->w2_t = window[n]->w2_t;
      tpg->w3_t = window[n]->w3_t;
      tpg->s2i = window[n]->s2i;
      tpg->s2j = window[n]->s2j;
      tpg->s2k = window[n]->s2k;
      tpg->sur_t = window[n]->sur_t;
      tpg->nsur_t = window[n]->nsur_t;
      tpg->bot_t = window[n]->bot_t;
      tpg->dHde1 = window[n]->dHde1;
      tpg->dHde2 = window[n]->dHde2;
      tpg->botz = window[n]->botz;
      tpg->costhcell = window[n]->costhcell;
      tpg->sinthcell = window[n]->sinthcell;

      tpd->u1 = windat->u1;
      tpd->u2 = windat->u2;
      tpd->w = windat->w;
      tpd->u1vm = windat->u1vm = master->u1vm;
      tpd->u2vm = windat->u2vm = master->u2vm;
      tpd->wtop = windat->wtop;
      tpd->wbot = windat->wbot;
      tpd->Kz = windat->Kz;
      tpd->eta = windat->eta;
      tpd->etab = windat->etab;
      tpd->detadt = windat->detadt;
      tpd->light = windat->light;
      tpd->temp = windat->temp;
      tpd->sal = windat->sal;
      m = tracer_find_index("salt", master->ntr, master->trinfo_3d);
      tpd->tr_wc[0] = windat->tr_wc[m];
      m = tracer_find_index("temp", master->ntr, master->trinfo_3d);
      tpd->tr_wc[1] = windat->tr_wc[m];
      for (m = 0; m < master->ntrvars; m++)
	tpd->tr_wc[2+m] = windat->tr_wc[master->trvm[m]];

      tpc->d1 = wincon->d1;
      tpc->d2 = wincon->d2;
      tpc->d3 = wincon->d3;
      tpc->d4 = wincon->d4;
      tpc->d7 = wincon->d7;
      tpc->w4 = wincon->w4;
      tpc->w5 = wincon->w5;
      tpc->w6 = wincon->w6;
      tpc->w7 = wincon->w7;
      tpc->w8 = wincon->w8;
      tpc->w9 = wincon->w9;
      tpc->w10 = wincon->w10;
      tpc->s1 = wincon->s1;
      tpc->s2 = wincon->s2;
      tpc->s3 = wincon->s3;
      tpc->i1 = wincon->i1;
      tpc->i2 = wincon->i2;
      tpc->i3 = wincon->i3;
      tpc->dz = wincon->dz;
      tpc->m2d = wincon->m2d;
      tpc->oldeta = wincon->oldeta;
      tpc->wgt = wincon->wgt;
      tpc->cdry_e1 = wincon->cdry_e1;
      tpc->tgrid = wincon->tgrid = NONE;
      window[n]->trans = geom->trans;

      tpc->nu = wincon->nu;
      tpc->nv = wincon->nv;
      tpc->nw = wincon->nw;
      tpc->crfxc = wincon->crfxc;
      tpc->crfyc = wincon->crfyc;
      tpc->crfzc = wincon->crfzc;
      tpc->crfxf = wincon->crfxf;
      tpc->crfyf = wincon->crfyf;
      tpc->crfzf = wincon->crfzf;
      tpc->clxc = wincon->clxc;
      tpc->clyc = wincon->clyc;
      tpc->clzc = wincon->clzc;
      tpc->clxf = wincon->clxf;
      tpc->clyf = wincon->clyf;
      tpc->clzf = wincon->clzf;
      tpc->tr_mod = wincon->tr_mod;
      tpc->tr_mod_x = wincon->tr_mod_x;
      tpc->tr_mod_y = wincon->tr_mod_y;
      tpc->tr_mod_z = wincon->tr_mod_z;
    }
  }
  set_dz(tpg, tpd, tpc);
  master->togn = TOPRIGHT;

}

/* END set_tp()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialize a relaxation structure.                                */
/*-------------------------------------------------------------------*/
relax_info_t *relax_info_init(char *rname,  /* Relaxation filename   */
			      char *tcname, /* Time constant string  */
			      double dt,    /* Input dt              */
			      int sz1,      /* Size of array 1       */
			      int sz2       /* Size of array 2       */
			      )
{
  relax_info_t *relax = (relax_info_t *)malloc(sizeof(relax_info_t));
  memset(relax, 0, sizeof(relax_info_t));
  relax->val1 = NULL;
  relax->val2 = NULL;
  strcpy(relax->rlxn, rname);
  strcpy(relax->rlxtc, tcname); 
  relax->rlxdt = dt;
  relax->size = 0;
  if (sz1) {
    relax->val1 = d_alloc_1d(sz1);
    memset(relax->val1, 0, sz1 * sizeof(double));
    relax->size = sz1 - 1;
  }
  if (sz2) {
    relax->val2 = d_alloc_1d(sz2);
    memset(relax->val2, 0, sz2 * sizeof(double));
  }
  return relax;
}

/* END relax_info_init()                                             */
/*-------------------------------------------------------------------*/


/* initialise the custom exchanges */
void custom_init(hd_data_t* hdata)
{
  custom_function_t* custom_fnc;
  if(hdata->master->custom_fnstack != NULL && hdata->master->custom_fnstack->nfunctions)
  {
      emstag(LTRACE,"hd:hd_init:custom_init","starting custom init for n %d",hdata->master->custom_fnstack->nfunctions);
      custom_fnc = hdata->master->custom_fnstack->functions;
      while(custom_fnc != NULL)
      {
        if(custom_fnc->init != NULL)
          custom_fnc->init(hdata,custom_fnc);
        custom_fnc = custom_fnc->next;
      }
  }else
    emstag(LTRACE,"hd:hd_ini:custom_init","No functions assigned!");
}
