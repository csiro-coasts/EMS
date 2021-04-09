/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/master/hd_init.c
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
 *  $Id: hd_init.c 6712 2021-03-29 00:55:24Z her127 $
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

#define nm3d 24
#define nm2d 5
char *mn_3d[nm3d] = {
  "u1mean", "u2mean", "wmean", "Kzmean", "temp_mean", "salt_mean",
  "u1vmean", "u2vmean",
  "flux_e1", "flux_e2", "flux_w", "flux_kz",
  "u1_adv", "u1_hdif", "u1_vdif", "u1_cor", "u1_btp", "u1_bcp",
  "u2_adv", "u2_hdif", "u2_vdif", "u2_cor", "u2_btp", "u2_bcp"
};
int mf_3d[nm3d] = {
  VEL3D, VEL3D, VEL3D, KZ_M, TS, TS,
  VOLFLUX, VOLFLUX,
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

xytoij_tree_t *xytoij_init_sparse(geometry_t *geom);
void sigma_tracer_step_3d(master_t *master);
void reset_inner_consts(master_t *master, geometry_t *geom);
void inner_mean(double *var, int cs, int cd);
void init_fetch(master_t *master, parameters_t *params, unsigned long **flag);
int sign(double x);
void trans_sourcegrid_init(parameters_t *params, master_t *master, geometry_t *geom,
			   geometry_t **window, window_t **windat, win_priv_t **wincon,
			   dump_data_t *dumpdata);
int get_st_map(master_t *master, geometry_t *wd, geometry_t *ws, double **sx, double **sy,
	       double **gridx, double **gridy, int *s2t, double *xoset, double *yoset);
void get_ts_map(geometry_t *geom, geometry_t **window, win_priv_t **wincon,
		geometry_t *ws, win_priv_t *wc, double **gridx, double **gridy, int of);
void set_trans_grid(int nce1, int nce2, double **sx, double **sy, 
		    double **cellx, double **celly, 
		    double **gridx, double **gridy,
		    double **u1x, double **u1y, double **u2x, double **u2y, geometry_t *geom);
void s2ij_setghosts(geometry_t *geom);
int check_quad(geometry_t *wt, xytoij_tree_t *t_xyij_tree, double **sx, double **sy,
	       int i, int j, double *x, double *y, int quad);	       
int s2t_weightsi(geometry_t *ws, double *dz, double x, double y, double z, 
		 double **wgt, int ci);
void set_tp(geometry_t *geom, master_t *master, geometry_t **window);
int find_nearest(geometry_t *geom, double x, double y);

/*------------------------------------------------------------------*/
/* Initialise the hydrodynamic structure                            */
/*------------------------------------------------------------------*/
hd_data_t *hd_init(FILE * prmfd)
{
  int dumpf = 0;
  int icdfid = 0;               /* input cdf file id */
  int timeindex = 0;            /* record index in dump file */

  /* Read the input parameters for this run */
  if (autof == 0) {
    params = params_read(prmfd);
  } else if (autof == 8)
    params = params_read_t(prmfd);
  else {
    params = auto_params(prmfd, autof);
  }

  if (params->runmode & (MANUAL | RE_ROAM | TRANS)) {
    /* Open input dump, get grid size and check compatibilities */
    icdfid = dump_open(params, params->idumpname, 1);
    /* Using the model start time, confirm that we have a dump in */
    /* the input file.  */
    timeindex = dump_choose_by_time(params, icdfid, schedule->start_time);
    /* Read bathymetry */
    params->topo = bathy_read(params, icdfid, timeindex);
    if (params->reset_bathyf) params->topo = set_bathy(params, icdfid);
  } else {
    /* Set and check the bathymetry array */
    params->topo = set_bathy(params, icdfid);
  }
  /* Create the id code                                             */
  create_code(params);

  /* Make the sparse map */
  build_sparse_map(params);

  /* Read the input dump data into the geometry and master */
  if (params->runmode & (AUTO | DUMP)) {
    /* Set the bottom depth on open boundaries */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);
    /* Set parameters to default values */
    autoset(params, master);

    /* Reset the post-tracer ROAM parameterisation if required */
    if (params->runmode & ROAM) autoset_roam(params, master, window);

    /* Set the ghost cells with valid data */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);

  }
  if (params->runmode & (MANUAL | RE_ROAM | TRANS)) {
    if (dumpf)
      dump_read(geom, params, master, icdfid, timeindex);
    else
      dumpdata_read(geom, params, master, dumpdata, icdfid, timeindex);
    dump_close(icdfid);

    /* Initialise the velocity if required */
    vel_init(geom, params, master);
    /* Set the ghost and OBC cells with valid data */
    master_setghosts(geom, master, dumpdata->cellx, dumpdata->celly,
		     dumpdata->gridx, dumpdata->gridy, NULL, NULL);
  }
  if (params->runmode & RE_NOTS)
    temp_salt_init(params, master);

  if (params->runmode & (AUTO | DUMP) && params->dozoom & PRECOND)
    hgrid_zoom_m(geom, window, dumpdata);

  /* Initialize the geometry and master structures */
  compute_constants(params, geom, master);

#ifdef HAVE_WAVE_MODULE
  /* Initialise the fetch if required                                */
  if (params->do_wave)
    init_fetch(master, params, dumpdata->flag[geom->nz-1]);
#endif

  /* Initialize the dump data structure */
  if (dumpf && mpi_rank == 0)
    dumpdata_init(dumpdata, geom, master);

  /* Initialize the mixing scheme */
  closure_init(params, master);

  /* Initialise the x,y to i,j conversion routine */
  /* xytoij_init in clibs uses 2D array as arguments */
  /* master->xyij_tree = xytoij_init_sparse(geom); */
  master->xyij_tree = grid_xytoij_init(dumpdata->gridx, dumpdata->gridy,
                                       geom->nce1, geom->nce2);

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

  /* Initialze the boundary conditions */
  bdry_custom_m(params, geom, master);
  bdry_custom_w(geom, window);

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
  /*tracer_dhw_init(master);*/

  /* Initialise the geometric arrays in the windows with master */
  /* geometry data.  */
  window_init(geom, window);

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
  pre_run_setup(master, window, windat, wincon);

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

  /* Initialise the window data transfer buffers' function pointers */
  win_data_init_transfer_buf_funcs(master, window, windat, wincon, geom->nwindows);

  write_window_map(window, params);

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
    char tag[MAXSTRLEN];
    int ocdfid;
    if (params->runmode & AUTO)
      sprintf(tag, "%s.nc", params->oname);
    else
      strcpy(tag, params->oname);
    ocdfid = dump_create(dumpdata, tag);
    dumpdata_fill(geom, master, dumpdata);
    set_bdry_flags(geom, dumpdata);
    if (params->runmode & DUMP)
      write_bdry(params, geom);
    dump_write(dumpdata, ocdfid, 0L);
    dump_close(ocdfid);
    if (params->runmode & AUTO)
      params_write(params, dumpdata);
    if (DEBUG("init_m"))
      dlog("init_m", "\nDumpfile %s created OK\n\n", tag);
  }
  write_mom_grid(dumpdata);
  write_roms_grid(dumpdata);
  if (strlen(params->wind_file))
    dump_windows(master, window, params->wind_file, params->prmname);

  if (params->runmode & EXIT) {
    windows_clear(hd_data);
    geom_free(master, UNUSED);
    timeseries_end();
    master_end(master);
    sched_end(schedule);
    fclose(prmfd);
    /*UR-ADDED */
    master_free(master);
    debug_end();
    exit(0);
  }

  custom_init(hd_data);

  /* Register an interest with the scheduler for dumping the */
  /* output to file.  */
  if (!(params->runmode & DUMP)) {
    /*
    sched_register(schedule, "dumps", dump_init,
                   dump_event, dump_cleanup, master, NULL, NULL);
    */
    if (params->tmode & SP_CHECK)
      trans_check_dump(master, dumpdata, params->trans_data);
    else
      sched_register(schedule, "dumps", dump_init,
		     dump_event, dump_cleanup, master,
		     dump_dispatch, dump_progress);
  }

  /* Initialise the DHW diagnostic. This must occur after the dumps */
  /* have been scheduled so that the scheduler performs it prior to */
  /* dumping. This is important for daily NRT operation.            */
  tracer_dhw_init(master);

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

  geom_free(master, UNUSED);
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

  sched_deregister(schedule, "dumps");

  tracer_reset_end(master);

  tracer_relax_end(master);

  ptrack_end();

}

/* END master_end()                                                  */
/*-------------------------------------------------------------------*/


/*UR----------------------------------------------------------------*/
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

/* END master_free()                                                  */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/*UR-ADDED function to add a late reference for cleanup              */
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
/* Routine to set up the x,y to i,j conversion xytoij_tree_t. Grid   */
/* arrays must be stored in 2D temporary arrays for the conversion.  */
/*-------------------------------------------------------------------*/
xytoij_tree_t *xytoij_init_sparse(geometry_t *geom  /* Global geometry
                                                       structure */
  )
{
  double **gridx = d_alloc_2d(geom->nfe1, geom->nfe2);
  double **gridy = d_alloc_2d(geom->nfe1, geom->nfe2);
  xytoij_tree_t *xytoij;
  s2c_2d(geom, geom->gridx, gridx, geom->nfe1, geom->nfe2);
  s2c_2d(geom, geom->gridy, gridy, geom->nfe1, geom->nfe2);
  xytoij = grid_xytoij_init(gridx, gridy, geom->nce1, geom->nce2);

  /*UR-TODO these are used in xytoij but released here ?? */
  d_free_2d(gridx);
  d_free_2d(gridy);
  return (xytoij);
}

/* END xytoij_init_sparse()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the constants in the master geometry and data file */
/* from information read in through the input dumpfile and parameter */
/* file.                                                             */
/*-------------------------------------------------------------------*/
void compute_constants(parameters_t *params,  /* Input parameter data
                                                 structure */
                       geometry_t *geom,  /* Sparse global geometry
                                             structure */
                       master_t *master /* Master data structure */
  )
{
  int c, cc;                    /* Cell coordinates */
  int c2, c3;                   /* 2D and 3D sparse coordinates */
  int cb;                       /* Bottom sparse coordinates */
  int xm1, xp1;                 /* Sparse locations at i+1 and i-1 */
  int ym1, yp1;                 /* Sparse locations at j+1 and j-1 */
  double dh;                    /* Thickness of bottom cell */
  double v;                     /* Dummy */
  double top, bot;              /* Surface and bottom depths */
  int n = 0;                    /* Mean grid spacing counter */
  int tn;                       /* Tracer counter */
  int ns;                       /* Number of smoothing passes */
  double *w1;                   /* Dummy */

  /* Set the master parameters from the parameter data structure */
  /* Time stepping constants and variables */
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

  /* Flags */
  master->trasc = params->trasc;
  master->momsc = params->momsc;
  if (master->momsc & ANGULAR3D)
    master->momsc2d = ORDER2;
  else
    master->momsc2d = master->momsc;
  master->ultimate = params->ultimate;
  master->osl = params->osl;
  master->smagorinsky = params->smagorinsky;
  master->sue1 = params->sue1;
  master->sue2 = params->sue2;
  master->kue1 = params->kue1;
  master->kue2 = params->kue2;
  master->bsue1 = params->bsue1;
  master->bsue2 = params->bsue2;
  master->bkue1 = params->bkue1;
  master->bkue2 = params->bkue2;
  master->smag_smooth = params->smag_smooth;
  master->diff_scale = params->diff_scale;
  master->visc_method = params->visc_method;
  if (params->roammode & (A_ROAM_R1|A_ROAM_R2|A_ROAM_R3|A_ROAM_R4) && params->nwindows > 1)
    master->visc_method |= ROAM_WIN;
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
  master->nwindows = params->nwindows;
  master->dbc = 0;
  master->dbgtime = 0.0;
  if(params->dbk >= 0 && params->dbj >= 0 && params->dbi >= 0) {
    if(params->dbk < params->nz && params->dbj < params->nce2 && 
       params->dbi < params->nce1)
      master->dbc = geom->map[params->dbk][params->dbj][params->dbi];
    else
      hd_warn("Can't assign debug location (%d %d %d) to sparse location.\n",
	      params->dbi, params->dbj, params->dbk);
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
  else if (params->decf & DEC_U2)
    master->decv = master->u2;
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

  master->mixlayer = params->mixlayer;
  master->show_layers = params->show_layers;
  master->lnm = fabs(params->lnm);
  master->vorticity = params->vorticity;
  master->numbers = params->numbers;
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
  master->swr_type = params->swr_type;
  master->ndhw = params->ndhw;
  if (master->ndhw) {
    int m;
    master->dhwf = i_alloc_1d(master->ndhw);
    master->dhwh = d_alloc_1d(master->ndhw);
    for (m = 0; m < params->ndhw; m++) {
      master->dhwf[m] = params->dhwf[m];
      master->dhwh[m] = params->dhwh[m];
    }
  }
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

  /* Means */
  master->means = params->means;
  if (strlen(params->means_dt)) {
    if (tm_scale_to_secs(params->means_dt, &master->means_dt))
      master->means_next = master->t + master->means_dt;
    else if (contains_token(params->means_dt, "YEARLY")) {
      master->means_dt = YEARLY;
      master->means_next = master->t + next_year(master->t, master->timeunit);
    } else if (contains_token(params->means_dt, "SEASONAL")) {
      master->means_dt = SEASONAL;
      master->meancs = d_alloc_1d(13);
      memset(master->meancs, 0, 13 * sizeof(double));
      master->means_next = master->t + next_season(master->t, 
						   master->timeunit, &c);
    } else if (contains_token(params->means_dt, "MONTHLY")) {
      master->means_dt = MONTHLY;
      master->meancs = d_alloc_1d(13);
      memset(master->meancs, 0, 13 * sizeof(double));
      master->means_next = master->t + next_month(master->t, 
						  master->timeunit, &c);
    } else if (contains_token(params->means_dt, "DAILY")) {
      master->means_dt = DAILY;
      master->meancs = d_alloc_1d(366);
      memset(master->meancs, 0, 366 * sizeof(double));
      master->means_next = master->t + next_day(master->t, 
						master->timeunit, &c);
    } else
      master->means_dt = 0.0;
    if (strlen(params->means_mc) && master->means_dt == SEASONAL || 
	       master->means_dt == MONTHLY || master->means_dt == DAILY) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      cc = parseline(params->means_mc, fields, MAXNUMARGS);
      for (c = 1; c <= cc; c++) {
	master->meancs[c] = atof(fields[c-1]);
      }
    }
    if (strlen(params->means_os)) {
      tm_scale_to_secs(params->stop_time, &master->means_os);
      for (cc = 1; cc <= geom->enonS; cc++) {
        master->meanc[cc] = master->means_os;
      }
    } else
      master->means_os = 0.0;
  } else
    master->means_dt = 0.0;

  if (master->means & TRANSPORT) {
    if (master->means_dt == 0.0)
      hd_warn("constants : Must set MEAN_DT using TRANSPORT with MEAN\n");
    else {
      if (fmod(master->means_dt, master->dt) != 0.0) {
        hd_warn
          ("constants : MEAN_DT is not integrally divisable by DT : resetting.\n");
        c = (int)master->means_dt / master->dt;
        master->means_dt = (double)c *master->dt;
      }
    }
    master->tratio = (int)master->means_dt / master->dt;
  }
  /* Store all 3d mean tracers in tm_3d */
  master->ntm_3d = 0;
  if (master->means & VEL3D) master->ntm_3d += 3;
  if (master->means & TS) master->ntm_3d += 2;
  if (master->means & KZ_M) master->ntm_3d += 1;
  if (master->means & FLUX) master->ntm_3d += 4;
  if (master->means & TENDENCY) master->ntm_3d += 12;
  if (master->means & VOLFLUX) master->ntm_3d += 2;
  if (master->means & MTRA3D) {
    master->ntm_3d += 1;
    if ((master->means_tra = tracer_find_index(params->means_tra, master->ntr, master->trinfo_3d)) < 0) {
      hd_warn("compute_constants: Can't find 3D tracer %s for MEAN tracer.\n", params->means_tra);
      master->means &= ~MTRA3D;
    }
  }
  if (master->ntm_3d) {
    master->tm_3d = i_alloc_1d(master->ntm_3d);
    ns = 0;
    for (tn = 0; tn < nm3d; tn++) {
      c = tracer_find_index(mn_3d[tn], master->ntr, master->trinfo_3d);
      if (master->means & mf_3d[tn] && c >= 0) {
	master->tm_3d[ns] = c;
	ns++;
      }
    }
  }

  /* Store all 2d mean tracers in tm_2d */
  master->ntm_2d = 0;
  if (master->means & VEL2D) master->ntm_2d += 2;
  if (master->means & ETA_M) master->ntm_2d += 1;
  if (master->means & WIND) master->ntm_2d += 2;
  if (master->means & MTRA2D) {
    master->ntm_2d += 1;
    if ((master->means_tra = tracer_find_index(params->means_tra, master->ntrS, master->trinfo_2d)) < 0) {
      hd_warn("compute_constants: Can't find 2D tracer %s for MEAN tracer.\n", params->means_tra);
      master->means &= ~MTRA2D;
    }
  }
  if (master->ntm_2d) {
    master->tm_2d = i_alloc_1d(master->ntm_2d);
    ns = 0;
    for (tn = 0; tn < nm2d; tn++) {
      c = tracer_find_index(mn_2d[tn], master->ntrS, master->trinfo_2d);
      if (master->means & mf_2d[tn] && c >= 0) {
	master->tm_2d[ns] = c;
	ns++;
      }
    }
  }

  /* Particles */
  master->do_pt = params->do_pt;
  if (master->do_pt)
    strcpy(master->ptinname, params->ptinname);

  master->u1_f = params->u1_f;
  master->u1av_f = params->u1av_f;
  master->u2_f = params->u2_f;
  master->u2av_f = params->u2av_f;
  master->orbital = params->orbital;
  master->waves = params->waves;
  master->save_force = params->save_force;
  master->dozoom = params->dozoom;
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

  master->e1bs = master->e1be = master->e2bs = master->e2be = NULL;
  if (params->nbl1) {
    master->nbl1 = params->nbl1;
    master->e1bs = i_alloc_1d(master->nbl1);
    master->e1be = i_alloc_1d(master->nbl1);
    master->e2bs = i_alloc_1d(master->nbl1);
    master->e2be = i_alloc_1d(master->nbl1);
    for (tn = 0; tn < master->nbl1; tn++) {
      if (sscanf(params->e1_blend[tn], "%d %d %d %d", &master->e1bs[tn], 
		 &master->e1be[tn], &master->e2bs[tn], &master->e2be[tn]) != 4) {
	hd_warn("Blending specification %s not recognized: no e1 blending performed.\n",
		params->e1_blend[tn]);
	master->nbl1 = 0;
      }
    }
  }
  if (params->nbl2) {
    master->nbl2 = params->nbl2;
    master->e2bs = i_alloc_1d(master->nbl2);
    master->e2be = i_alloc_1d(master->nbl2);
    master->e1bs = i_alloc_1d(master->nbl2);
    master->e1be = i_alloc_1d(master->nbl2);
    for (tn = 0; tn < master->nbl2; tn++) {
      if (sscanf(params->e2_blend[tn], "%d %d %d %d", &master->e2bs[tn], 
		 &master->e2be[tn], &master->e1bs[tn], &master->e1be[tn]) != 4) {
	hd_warn("Blending specification %s not recognized: no e2 blending performed.\n",
		params->e2_blend[tn]);
	master->nbl2 = 0;
      }
    }
  }

  /* Constants */
  master->g = 9.81;
  master->ambpress = params->ambpress;
  master->hmin = params->hmin;
  master->uf = params->uf;
  master->quad_bfc = params->quad_bfc;
  master->etamax = params->etamax;
  master->velmax = params->velmax;
  master->velmax2d = params->velmax2d;
  master->etadiff = params->etadiff;
  master->wmax = params->wmax;
  master->rampstart = params->rampstart;
  master->rampend = params->rampend;
  master->rampf = params->rampf;

  /* Tracer constants and variables */
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
    if (strcmp(params->trflux, tracer->name) == 0)
      master->trflux = tn;
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

  /* Transport mode : get additional variables to be reset */
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
      /* Fill the arrays with valid tracers */
      /* 3D variables */
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
      /* 2D variables */
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

  /* Mixing constants and variables */
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

  /* Heatflux variables */
  master->heatflux = params->heatflux;
  master->albedo = params->albedo;
  master->albedo_l = params->albedo_l;
  master->zref = params->zref;
  master->bulkf = params->bulkf;
  master->hfadf = params->hfadf;
  master->hftc = params->hftc;
  master->hf_ramp = params->hf_ramp;
  /* Override input file swr parameters if they are numbers */
  if (master->swr_attn) {
    double d1;
    int cc, c;
    master->attn_tr = -1;
    if ((tn = tracer_find_index(params->swr_attn, master->ntrS, master->trinfo_2d)) >= 0) {
      memcpy(master->swr_attn, master->tr_wcS[tn], geom->sgsizS * sizeof(double));
      master->attn_tr = tn;
    }
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
    master->tran_tr = -1;
    if (master->swr_tran && (tn = tracer_find_index(params->swr_tran, master->ntrS, master->trinfo_2d)) >= 0) {
      memcpy(master->swr_tran, master->tr_wcS[tn], geom->sgsizS * sizeof(double));
      master->tran_tr = tn;
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

  /* Saltflux variables */
  master->saltflux = params->saltflux;

  /* Wind variables */
  master->wind_dt = params->wind_dt;
  master->storm_dt = params->storm_dt;
  master->nstorm = params->nstorm;

  /* Elevation relaxation variables */
  master->etarlx = params->etarlx;
  if (params->etarlx & (RELAX|ALERT|BOUNDARY)) {
    master->eta_rlx = relax_info_init(params->etarlxn, params->etarlxtcs, 
				      params->etarlxdt, 0, 0);
    if ((tn = tracer_find_index("oeta", master->ntrS, master->trinfo_2d)) >= 0)
      master->eta_rlx->val1 = master->tr_wcS[tn];
    master->eta_rlx->rlx = params->etarlx;
  }
  /*
  strcpy(master->etarlxn, params->etarlxn);
  master->etarlxdt = params->etarlxdt;
  master->etarlxtc = params->etarlxtc;
  strcpy(master->etarlxtcs, params->etarlxtcs);
  */

  /* Velocity relaxation */
  master->velrlx = params->velrlx;
  if (params->velrlx & RELAX) {
    master->vel_rlx = relax_info_init(params->velrlxn, params->velrlxtcs, 
				      params->velrlxdt, geom->sgsiz, geom->sgsiz);
    master->vel_rlx->rlx = params->velrlx;
  }
  master->tide_r = params->tide_r;

  /* Alerts */
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
  /* Sediments */
  master->do_sed = params->do_sed;
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  /* Ecology */
  master->do_eco = params->do_eco;
  master->ecodt = params->ecodt;
#endif

#if defined(HAVE_WAVE_MODULE)
  /* Waves */
  master->do_wave = params->do_wave;
  master->wavedt = params->wavedt;
#endif

  /*-----------------------------------------------------------------*/
  /* Initialize the sigma variables */
  for (c = 1; c <= geom->enonS; c++) {  /* Sigma variables */
    /*
       master->Ds[c]=1.0; master->Hs[c]=0.0; master->Hn1[c]=1.0;
       master->Hn2[c]=1.0; */
  }

  /*-----------------------------------------------------------------*/
  /* Smooth the topography if required */
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
  /* Set the cell thickness for the master */
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
    int i, j;
    double **z0 = d_alloc_2d(params->nce1, params->nce2);
    c = 0;
    for (j = 0; j < geom->nce2; j++)
      for (i = 0; i < geom->nce1; i++)
	z0[j][i] = params->z0s[c++];
    c2s_2d(geom, master->z0, z0, geom->nce1, geom->nce2);
    d_free_2d(z0);
  } else {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      c2 = geom->m2d[c];
      master->z0[c2] = params->z0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the drag coefficient */
  for (c = 1; c < geom->enonS; c++)
    master->Cd[c] = master->quad_bfc;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    c2 = geom->m2d[c];
    cb = geom->bot_t[cc];
    dh = max(0.5 * master->dz[cb] * master->Ds[c2], master->hmin);
    v = log((dh + master->z0[c2]) / master->z0[c2]) / VON_KAR;
    if (master->quad_bfc < 0.0)
      master->Cd[c2] = -master->quad_bfc;
    else
      master->Cd[c2] = max(master->quad_bfc, 1.0 / (v * v));
  }
  ns = get_smoothing(params->smooth_v, "CD");
  for (n = 0; n < ns; n++)
    smooth(master, master->Cd, geom->w2_t, geom->b2_t);
  if ((v = get_scaling(params->scale_v, "CD")) != 1.0) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->Cd[c] *= v;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set variables at multiple ghost cells                           */
  for (cc = 1; cc <= geom->enonS; cc++) {
    c = cc;
    c3 = geom->mgc[c];
    while (c3) {
      geom->thetau1[geom->mgc[c]] = geom->thetau1[c];
      geom->thetau2[geom->mgc[c]] = geom->thetau2[c];
      c = c3;
      c3 = geom->mgc[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for cell centers */
  master->hmean1 = master->hmean2 = 1e-10;
  geom->totarea = 0.0;
  for (cc = 1; cc <= geom->n2_t; cc++) {
    c = geom->w2_t[cc];
    xp1 = geom->xp1[c];
    xm1 = geom->xm1[c];
    yp1 = geom->yp1[c];
    ym1 = geom->ym1[c];

    if(isnan(geom->h1acell[c])) {
      c3 = find_non_nan_c(geom, geom->h1acell, c);
      geom->h1acell[c] = geom->h1acell[c3];
      hd_warn("compute_constants: replacing h1acell NaN @ (%d %d) with %5.2f from (%d %d)", geom->s2i[c], geom->s2j[c], geom->h1acell[c], geom->s2i[c3], geom->s2j[c3]);
    }
    if(isnan(geom->h2acell[c])) {
      c3 = find_non_nan_c(geom, geom->h2acell, c);
      geom->h2acell[c] = geom->h2acell[c3];
      hd_warn("compute_constants: replacing h2acell NaN @ (%d %d) with %5.2f from (%d %d)", geom->s2i[c], geom->s2j[c], geom->h2acell[c], geom->s2i[c3], geom->s2j[c3]);
    }
    geom->cellarea[c] = geom->h1acell[c] * geom->h2acell[c];
    if (cc <= geom->b2_t)
      geom->totarea += geom->cellarea[c];

    geom->sinthcell[c] = (sin(geom->thetau1[c]) +
                          sin(geom->thetau2[c])) / 2.0;
    if (!isnan(geom->thetau1[xp1]) && geom->thetau1[xp1] != 0.)
      geom->sinthcell[c] = (geom->sinthcell[c] + sin(geom->thetau1[xp1])) / 2.0;
    if (!isnan(geom->thetau2[yp1]) && geom->thetau2[yp1] != 0.)
      geom->sinthcell[c] = (geom->sinthcell[c] + sin(geom->thetau2[yp1])) / 2.0;
    geom->costhcell[c] = (cos(geom->thetau1[c]) +
                          cos(geom->thetau2[c])) / 2.0;
    if (!isnan(geom->thetau1[xp1]) && geom->thetau1[xp1] != 0.)
      geom->costhcell[c] = (geom->costhcell[c] + cos(geom->thetau1[xp1])) / 2.0;
    if (!isnan(geom->thetau2[yp1]) && geom->thetau2[yp1] != 0.)
      geom->costhcell[c] = (geom->costhcell[c] + cos(geom->thetau2[yp1])) / 2.0;

    /*
    geom->sinthcell[c] = (sin(geom->thetau1[c]) +
                          sin(geom->thetau1[xp1]) +
                          sin(geom->thetau2[c]) +
                          sin(geom->thetau2[yp1])) / 4.0;
    geom->costhcell[c] = (cos(geom->thetau1[c]) +
                          cos(geom->thetau1[xp1]) +
                          cos(geom->thetau2[c]) +
                          cos(geom->thetau2[yp1])) / 4.0;
    */
    geom->dHde1[c] = 0.5 * (geom->botz[xp1] - geom->botz[xm1]);
    geom->dHde2[c] = 0.5 * (geom->botz[yp1] - geom->botz[ym1]);

    master->hmean1 += geom->h1acell[c];

    if(!isnan(geom->h2acell[c]))
      master->hmean2 += geom->h2acell[c];
    else
      hd_warn("compute_constants: h2acell NaN found (%d %d)", geom->s2i[c], geom->s2j[c]);
    n += 1;
  }

  /*-----------------------------------------------------------------*/
  /* Scale horizontal mixing coefficients by cell size */
  
  if (n) {
    master->hmean1 /= (double)n;
    master->hmean2 /= (double)n;
  }
  memset(master->u1kh, 0, geom->sgsiz * sizeof(double));
  memset(master->u2kh, 0, geom->sgsiz * sizeof(double));
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    c2 = geom->m2d[c];
    if (params->u1kh > 0.0) {
      if (params->diff_scale == LINEAR)
	master->u1kh[c] = fabs(params->u1kh * geom->h1acell[c2] /
			       master->hmean1);
      if (params->diff_scale == NONLIN)
	master->u1kh[c] = fabs(params->u1kh * 
			       geom->h1acell[c2] * geom->h1acell[c2] /
			       (master->hmean1 * master->hmean1));
    }
    if (params->u2kh > 0.0) {
      if (params->diff_scale == LINEAR)
	master->u2kh[c] = fabs(params->u2kh * geom->h2acell[c2] /
			       master->hmean2);
      if (params->diff_scale == NONLIN)
	master->u2kh[c] = fabs(params->u2kh * 
			       geom->h2acell[c2] * geom->h2acell[c2] /
			       (master->hmean2 * master->hmean2));
    }
  }

  master->u1vh0 = params->u1vh;
  master->u2vh0 = params->u2vh;
  master->u1kh0 = params->u1kh;
  master->u2kh0 = params->u2kh;
  memset(geom->dHde1, 0, geom->sgsizS * sizeof(double));
  memset(geom->dHde2, 0, geom->sgsizS * sizeof(double));
  for (c = 1; c < geom->enonS; c++) {
    if(geom->h1au1[c] == 0)geom->h1au1[c] = master->hmean1;
    if(geom->h2au2[c] == 0)geom->h2au2[c] = master->hmean2;
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for u1 momentum equations */
  memset(master->u1vh, 0, geom->sgsiz * sizeof(double));
  for (cc = 1; cc <= geom->b2_e1; cc++) {
    c = geom->w2_e1[cc];
    xp1 = geom->xp1[c];
    xm1 = geom->xm1[c];
    master->u1c1[c] =
      -1.0 / (geom->h1au1[c] * geom->h1au1[c] * geom->h2au1[c]);
    master->u1c3[c] =
      1.0 * (geom->h1acell[c] -
             geom->h1acell[xm1]) / (geom->h1au1[c] * geom->h1au1[c]);
    master->u1c4[c] =
      1.0 * (geom->h2acell[c] -
             geom->h2acell[xm1]) / (geom->h1au1[c] * geom->h2au1[c]);
    master->u1c5[c] =
      1.0 * (master->coriolis[xm1] + master->coriolis[c]) / 2;
    master->u1c6[c] = -1.0 / geom->h1au1[c];
    geom->sinthu1[c] = sin(geom->thetau1[c]);
    geom->costhu1[c] = cos(geom->thetau1[c]);
    geom->botzu1[c] = max(geom->botz[xm1], geom->botz[c]);
  }
  if (params->u1vh > 0.0) {
    for (cc = 1; cc <= geom->n3_e1; cc++) {
      c = geom->w3_e1[cc];
      c2 = geom->m2d[c];
      /* Note : horizontal diffusion coeffients are scaled to the grid */
      if (params->diff_scale == LINEAR)
	master->u1vh[c] = fabs(params->u1vh * geom->h1au1[c2] / 
			       master->hmean1);
      if (params->diff_scale == NONLIN)
	master->u1vh[c] = fabs(params->u1vh * 
			       geom->h1au1[c2] * geom->h1au1[c2] / 
			       (master->hmean1 * master->hmean1));
    }
  } else {
    if (params->bsue1)
      scale_hdiff_m(master, master->u1vh, params->sue1, params->bsue1, 1);
  }

  /* Set the ghost cells */
  for (cc = 1; cc <= geom->nbpte1S; cc++) {
    c = geom->bpte1S[cc];
    c2 = geom->bine1S[cc];
    geom->botzu1[c] = geom->botzu1[c2];
    geom->botzu1[c] = HUGE_VAL;
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for u2 momentum equations */
  memset(master->u2vh, 0, geom->sgsiz * sizeof(double));
  for (cc = 1; cc <= geom->b2_e2; cc++) {
    c = geom->w2_e2[cc];
    yp1 = geom->yp1[c];
    ym1 = geom->ym1[c];
    master->u2c1[c] =
      -1.0 / (geom->h1au2[c] * geom->h2au2[c] * geom->h2au2[c]);
    master->u2c3[c] =
      1.0 * (geom->h1acell[c] -
             geom->h1acell[ym1]) / (geom->h1au2[c] * geom->h2au2[c]);
    master->u2c4[c] =
      1.0 * (geom->h2acell[c] -
             geom->h2acell[ym1]) / (geom->h2au2[c] * geom->h2au2[c]);
    master->u2c5[c] =
      -1.0 * (master->coriolis[c] + master->coriolis[ym1]) / 2;
    master->u2c6[c] = -1.0 / geom->h2au2[c];
    geom->sinthu2[c] = sin(geom->thetau2[c]);
    geom->costhu2[c] = cos(geom->thetau2[c]);
    geom->botzu2[c] = max(geom->botz[ym1], geom->botz[c]);
  }

  if (params->u2vh > 0.0) {
    for (cc = 1; cc <= geom->n3_e2; cc++) {
      c = geom->w3_e2[cc];
      c2 = geom->m2d[c];
      if (params->diff_scale == LINEAR)
	master->u2vh[c] = fabs(params->u2vh * geom->h2au2[c2] / 
			       master->hmean2);
      if (params->diff_scale == NONLIN)
	master->u2vh[c] = fabs(params->u2vh * 
			       geom->h2au2[c2] * geom->h2au2[c2] / 
			       (master->hmean2 * master->hmean2));
    }
  } else {
    if (params->bsue2)
      scale_hdiff_m(master, master->u2vh, params->sue2, params->bsue2, 2);
  }

  /* Set the ghost cells */
  for (cc = 1; cc <= geom->nbpte2S; cc++) {
    c = geom->bpte2S[cc];
    c2 = geom->bine2S[cc];
    geom->botzu2[c] = geom->botzu2[c2];
    geom->botzu2[c] = HUGE_VAL;
  }

  /*-----------------------------------------------------------------*/
  /* Reset constants over inner explicit maps                        */
  if (master->exmapf)
    reset_inner_consts(master, geom);

  /*-----------------------------------------------------------------*/
  /* Optimised horizontal mixing */
  if(params->diff_scale == AUTO) {
    double hf = 0.05;             /* Factor for horizontal diffusion */
    double step = 1;              /* Integral step of diffusion > 1  */
    double hmax = 1e10;
    double d1, d2;
    int u1khf = 0, u1vhf = 0, u2khf = 0, u2vhf = 0, i1, cs;
    if (params->u1kh >= 0.0)
      u1khf = 1;
    if (params->u1vh >= 0.0)
      u1vhf = 1;
    if (params->u2kh >= 0.0)
      u2khf = 1;
    if (params->u2vh >= 0.0)
      u2vhf = 1;
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      hmax = 1.0 / ((1.0 / (geom->h1acell[cs] * geom->h1acell[cs]) +
		   1.0 / (geom->h2acell[cs] * geom->h2acell[cs])) * 4.0 *
		  params->grid_dt);
      i1 = (int)hmax / (int)step;
      hmax = step * (double)i1;
      d1 = 0.01 * geom->h1acell[cs] * geom->h1acell[cs] / params->grid_dt;
      d2 = 0.01 * geom->h2acell[cs] * geom->h2acell[cs] / params->grid_dt;
      i1 = (int)d1 / (int)step;
      if (u1khf)
	master->u1kh[c] = step * (double)i1;
      if (u1vhf)
	master->u1vh[c] = step * (double)i1;
      if (master->u1kh[c] == 0.0) master->u1kh[c] = d1;
      if (master->u1vh[c] == 0.0) master->u1vh[c] = d1;
      i1 = (int)d2 / (int)step;
      if (u2khf)
	master->u2kh[c] = step * (double)i1;
      if (u2vhf)
	master->u2vh[c] = step * (double)i1;
      if (master->u2kh[c] == 0.0) master->u2kh[c] = d2;
      if (master->u2vh[c] == 0.0) master->u2vh[c] = d2;
      /* Set limits */
      if (u1khf) {
	if (master->u1kh[c] > hmax)
	  master->u1kh[c] = hmax;
	if (master->u1kh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u1kh[c] = step * (double)i1;
	}
      }
      if (u1vhf) {
	if (master->u1vh[c] > hmax)
	  master->u1vh[c] = hmax;
	if (master->u1vh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
        master->u1vh[c] = step * (double)i1;
	}
      }
      if (u2khf) {
	if (master->u2kh[c] > hmax)
	  master->u2kh[c] = hmax;
	if (master->u2kh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u2kh[c] = step * (double)i1;
	}
      }
      if (u2vhf) {
	if (master->u2vh[c] > hmax)
	  master->u2vh[c] = hmax;
	if (master->u2vh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u2vh[c] = step * (double)i1;
	}
      }
    }
  }

  /* Scale horizontal mixing */
  if ((v = get_scaling(params->scale_v, "U1VH")) != 1.0) {
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      master->u1vh[c] *= v;
    }
  }
  if ((v = get_scaling(params->scale_v, "U2VH")) != 1.0) {
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      master->u2vh[c] *= v;
    }
  }
  if ((v = get_scaling(params->scale_v, "U1KH")) != 1.0) {
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      master->u1kh[c] *= v;
    }
  }
  if ((v = get_scaling(params->scale_v, "U2KH")) != 1.0) {
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      master->u2kh[c] *= v;
    }
  }

  for (c = 1; c < geom->enonS; c++) {
    if (master->u1vhin)
      master->u1vhin[c] = master->u1vh[c];
    if (master->u2vhin)
      master->u2vhin[c] = master->u2vh[c];
  }

  /* Set the sponge zones if required */
  w1 = d_alloc_1d(geom->sgsiz);
  set_sponge(geom, master->u1vh, master->u2vh, params->grid_dt, w1);
  d_free_1d(w1);

  /* Smooth horizontal mixing                                        */
  ns = get_smoothing(params->smooth_v, "U1VH");
  for (n = 0; n < ns; n++)
    smooth3(master, master->u1vh, geom->w3_e1, geom->n3_e1);
  ns = get_smoothing(params->smooth_v, "U2VH");
  for (n = 0; n < ns; n++)
    smooth3(master, master->u2vh, geom->w3_e2, geom->n3_e2);
  ns = get_smoothing(params->smooth_v, "U1KH");
  for (n = 0; n < ns; n++)
    smooth3(master, master->u1kh, geom->w3_t, geom->n3_t);
  ns = get_smoothing(params->smooth_v, "U2KH");
  for (n = 0; n < ns; n++)
    smooth3(master, master->u2kh, geom->w3_t, geom->n3_t);

  /* Set the ghost cells */
  for (cc = 1; cc <= geom->nbpte1; cc++) {
    c = geom->bpte1[cc];
    c2 = geom->bine1[cc];
    master->u1vh[c] = master->u1vh[c2];
  }
  for (cc = 1; cc <= geom->nbpte2; cc++) {
    c = geom->bpte2[cc];
    c2 = geom->bine2[cc];
    master->u2vh[c] = master->u2vh[c2];
  }

  /*-----------------------------------------------------------------*/
  /* Precalculated constants for grid corners */
  for (cc = 1; cc <= geom->n2_t; cc++) {
    double maxbz = -1e10;
    double sum = 0.0;
    double d1 = 0.0;
    c = geom->w2_t[cc];
    xm1 = geom->xm1[c];
    ym1 = geom->ym1[c];
    /* Cell location (i-1,j-1) */
    maxbz = max(maxbz, geom->botz[geom->ym1[xm1]]);
    if (!isnan(geom->botz[geom->ym1[xm1]])) {
      sum += geom->botz[geom->ym1[xm1]];
      d1 += 1.0;
    }
    /* Cell location (i-1,j) */
    maxbz = max(maxbz, geom->botz[xm1]);
    if (!isnan(geom->botz[xm1])) {
      sum += geom->botz[xm1];
      d1 += 1.0;
    }
    /* Cell location (i,j-1) */
    maxbz = max(maxbz, geom->botz[ym1]);
    if (!isnan(geom->botz[ym1])) {
      sum += geom->botz[ym1];
      d1 += 1.0;
    }
    /* Cell location (i,j) */
    maxbz = max(maxbz, geom->botz[c]);
    if (!isnan(geom->botz[c])) {
      sum += geom->botz[c];
      d1 += 1.0;
    }
    geom->botzgrid[c] = maxbz;  /* Minimum depth */
    if (d1)
      geom->botzgrid[c] = sum / d1;  /* Average depth */
  }

  /*-----------------------------------------------------------------*/
  /* Initialize arrays */
  memset(master->waterss2d, 0, (geom->enonS + 1) * sizeof(double));
  memset(master->wbot, 0, (geom->enonS + 1) * sizeof(double));
  init_wvel_bounds(geom, master); /* Initialize wbot & detadt */

  /* Scale tracers if required (hardwire option)                     */
  /*
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    bot = min(fabs(geom->gridz[c]), 23.0);
    master->sal[c] += (1.15 - 0.52 * bot / 23.0);
    master->temp[c] += (1.5 - 2.0 * bot / 23.0);
  }
  */
}

/* END compute_constants()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the master data structure          */
/*-------------------------------------------------------------------*/
master_t *master_build(parameters_t *params, geometry_t *geom)
{
  master_t *master;             /* Master data structure */
  int tn,tt;

  /* Allocate memory for the master */
  master = master_alloc();

  /* Transfer input data to the master */
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
  strcpy(master->tracerdata, params->tracerdata);
  strcpy(master->timeunit, params->timeunit);
  master->tsfile_caching = params->tsfile_caching;
  master->zoomf = geom->zoomf;
  master->zmfe1 = geom->zmfe1;
  master->zmfe2 = geom->zmfe2;

  /* Constants */
  master->Cd = d_alloc_1d(geom->sgsizS);
  master->u1c1 = d_alloc_1d(geom->sgsizS);
  master->u1c3 = d_alloc_1d(geom->sgsizS);
  master->u1c4 = d_alloc_1d(geom->sgsizS);
  master->u1c5 = d_alloc_1d(geom->sgsizS);
  master->u1c6 = d_alloc_1d(geom->sgsizS);
  master->u2c1 = d_alloc_1d(geom->sgsizS);
  master->u2c3 = d_alloc_1d(geom->sgsizS);
  master->u2c4 = d_alloc_1d(geom->sgsizS);
  master->u2c5 = d_alloc_1d(geom->sgsizS);
  master->u2c6 = d_alloc_1d(geom->sgsizS);
  master->coriolis = d_alloc_1d(geom->sgsizS);
  master->one = d_alloc_1d(geom->sgsizS);
  for (tn = 1; tn <= geom->enonS; tn++)
    master->one[tn] = 1.0;

  /* Vertical grid geometry */
  master->dz = d_alloc_1d(geom->sgsiz);
  master->dzu1 = d_alloc_1d(geom->sgsiz);
  master->dzu2 = d_alloc_1d(geom->sgsiz);
  master->depth_e1 = d_alloc_1d(geom->sgsizS);
  master->depth_e2 = d_alloc_1d(geom->sgsizS);
  master->topz = d_alloc_1d(geom->sgsizS);
  master->topk = i_alloc_1d(geom->sgsizS);
  master->botk = i_alloc_1d(geom->sgsizS);

  /* SIGMA variables */
  if (params->sigma) {
    master->Ds = d_alloc_1d(geom->sgsizS);
    master->Hs = d_alloc_1d(geom->sgsizS);
    master->Hn1 = d_alloc_1d(geom->sgsizS);
    master->Hn2 = d_alloc_1d(geom->sgsizS);
  } else {
    master->Ds = master->one;
    master->Hs = master->one;
    master->Hn1 = master->one;
    master->Hn2 = master->one;
  }

  /* e1 velocity variables */
  master->u1 = d_alloc_1d(geom->sgsiz);
  master->nu1 = d_alloc_1d(geom->sgsiz);
  master->u1bot = d_alloc_1d(geom->sgsizS);
  master->u1flux3d = d_alloc_1d(geom->sgsiz);
  master->u1adv = d_alloc_1d(geom->sgsizS);
  master->u1inter = d_alloc_1d(geom->sgsizS);
  master->u1av = d_alloc_1d(geom->sgsizS);
  master->nu1av = d_alloc_1d(geom->sgsizS);
  master->u1flux = d_alloc_1d(geom->sgsizS);

  /* e2 velocity variables */
  master->u2 = d_alloc_1d(geom->sgsiz);
  master->nu2 = d_alloc_1d(geom->sgsiz);
  master->u2bot = d_alloc_1d(geom->sgsizS);
  master->u2flux3d = d_alloc_1d(geom->sgsiz);
  master->u2adv = d_alloc_1d(geom->sgsizS);
  master->u2inter = d_alloc_1d(geom->sgsizS);
  master->u2av = d_alloc_1d(geom->sgsizS);
  master->nu2av = d_alloc_1d(geom->sgsizS);
  master->u2flux = d_alloc_1d(geom->sgsizS);

  /* Grid refinement variables */
  if (params->dozoom) {
    master->u1flux_z = d_alloc_1d(geom->sgsizS);
    master->u2flux_z = d_alloc_1d(geom->sgsizS);
    master->d1 = d_alloc_1d(geom->sgsizS);
    master->d2 = d_alloc_1d(geom->sgsizS);
    master->d3 = d_alloc_1d(geom->sgsizS);
    master->d4 = d_alloc_1d(geom->sgsizS);
    master->d5 = d_alloc_1d(geom->sgsizS);
    master->d6 = d_alloc_1d(geom->sgsizS);
  }

  /* Vertical velocity variables */
  master->w = d_alloc_1d(geom->sgsiz);
  master->wtop = d_alloc_1d(geom->sgsizS);
  master->wbot = d_alloc_1d(geom->sgsizS);

  master->u1b = d_alloc_1d(geom->sgsiz);
  master->u2b = d_alloc_1d(geom->sgsiz);
  master->u1avb = d_alloc_1d(geom->sgsizS);
  master->u2avb = d_alloc_1d(geom->sgsizS);
  master->etab = d_alloc_1d(geom->sgsizS);

  /* 3D Tracer constants and variables */
  master->tr_wc = d_alloc_2d(geom->sgsiz, master->ntr);
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
  if (params->ntr) {
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

  /* Set the 3D tracer pointers */
  init_tracer_3d(params, master);
  /* Set the mean tracer to 3D auto tracers if required              */
  if (params->means & MTRA3D) {
    for (tn = 0; tn < master->ntr; tn++) {
      if (contains_token(params->means_tra, master->trinfo_3d[tn].name) != NULL) {
	strcpy(params->means_tra, master->trinfo_3d[tn].name);
      }
    }
  }

  /* 2D Tracer constants and variables */
  if (master->ntrS) {
    master->tr_wcS = d_alloc_2d(geom->sgsizS, master->ntrS);
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

  /* Sediment tracer constants and variables */
  if (master->nsed) {
    master->tr_sed = d_alloc_3d(geom->sgsizS, geom->sednz, master->nsed);
    master->trinfo_sed =
      (tracer_info_t *)malloc(sizeof(tracer_info_t) * master->nsed);
    memset(master->trinfo_sed, 0, sizeof(tracer_info_t) * master->nsed);
    for (tn = 0; tn < master->nsed; tn++) {
      tracer_copy(&master->trinfo_sed[tn], &params->trinfo_sed[tn]);
      strcpy(master->trinfo_sed[tn].name, params->trinfo_sed[tn].name);
    }
    init_tracer_sed(params, master);
  }

  /* Do any custom initialisation of the sediments tracer */
#if defined(HAVE_SEDIMENT_MODULE)
  sed_tracer_custom(master);
#endif

  /* Mixing constants and variables */
  master->Kz = d_alloc_1d(geom->sgsiz);
  master->Vz = d_alloc_1d(geom->sgsiz);
  master->z0 = d_alloc_1d(geom->sgsizS);

  /* Surface elevation variables */
  master->eta = d_alloc_1d(geom->sgsizS);
  master->oldeta = d_alloc_1d(geom->sgsizS);
  master->waterss = d_alloc_1d(geom->sgsiz);
  master->waterss2d = d_alloc_1d(geom->sgsizS);
  master->detadt = d_alloc_1d(geom->sgsizS);
  master->wdiff2d = d_alloc_1d(geom->sgsizS);

  /* Density variables */
  master->dens = d_alloc_1d(geom->sgsiz);
  master->dens_0 = d_alloc_1d(geom->sgsiz);
  master->topdensu1 = d_alloc_1d(geom->sgsizS);
  master->topdensu2 = d_alloc_1d(geom->sgsizS);
  master->densavu1 = d_alloc_1d(geom->sgsizS);
  master->densavu2 = d_alloc_1d(geom->sgsizS);

  /* Atmospheric variables */
  master->wind1 = d_alloc_1d(geom->sgsizS);
  master->wind2 = d_alloc_1d(geom->sgsizS);
  master->windspeed = d_alloc_1d(geom->sgsizS);
  master->winddir = d_alloc_1d(geom->sgsizS);
  if (params->wind_dt && params->storm_dt) {
    master->swind1 = d_alloc_1d(geom->sgsizS);
    master->swind2 = d_alloc_1d(geom->sgsizS);
  }
  master->patm = d_alloc_1d(geom->sgsizS);

  /* Heatflux variables */
  if (params->heatflux & (NET_HEAT | ADVANCED | INVERSE | COMP_HEAT | COMP_HEAT_MOM)) {
    master->heatf = d_alloc_1d(geom->sgsizS);
    if (!master->swr)
      master->swr = d_alloc_1d(geom->sgsizS);
  }

  /* Horizontal diffusion variables                                  */
  /* smagcode = U1_A  : No Smagorinsky; u1vh allocated               */
  /* smagcode = U2_A  : No Smagorinsky; u2vh allocated               */
  /* smagcode = U1_AK : No Smagorinsky; u1kh allocated               */
  /* smagcode = U2_AK : No Smagorinsky; u2kh allocated               */
  /* smagcode = U1_SP : Smagorinsky; u1vh = sdc, no sponges          */
  /* smagcode = U1_SA : Smagorinsky; u1vh allocated, uses sponges    */
  /* smagcode = U2_SP : Smagorinsky; u2vh = sdc, no sponges          */
  /* smagcode = U2_SA : Smagorinsky; u2vh allocated, uses sponges    */
  /* smagcode = U1_SPK: Smagorinsky; u1kh = sdc, no sponges          */
  /* smagcode = U1_SAK: Smagorinsky; u1kh allocated, uses sponges    */
  /* smagcode = U2_SPK: Smagorinsky; u2kh = sdc, no sponges          */
  /* smagcode = U2_SAK: Smagorinsky; u2kh allocated, uses sponges    */
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
    if (master->smagcode & U1_SA)
      master->u1vh = d_alloc_1d(geom->sgsiz);
    else {
      master->u1vh = master->sdc;
      master->smagcode |= U1_SP;
    }
  } else {
    master->u1vh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U1_A;
  }
  if (params->smagorinsky > 0.0 && params->u2vh < 0.0) {
    int n;
    if (params->smagorinsky == 1.0 && params->sue2 != 1.0) 
      master->smagcode |= U2_SA;
    for (n = 0; n < params->nobc; ++n)
      if (params->open[n]->sponge_zone_h)
	master->smagcode |= U2_SA;
    if (master->smagcode & U2_SA)
      master->u2vh = d_alloc_1d(geom->sgsiz);
    else {
      master->u2vh = master->sdc;
      master->smagcode |= U2_SP;
    }
  } else {
    master->u2vh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U2_A;
  }
  /*
  if (params->smagorinsky > 0.0 && params->u1vh < 0.0) {
    master->u1vh = master->sdc;
    master->smagcode |= U1_SP;
  } else {
    master->u1vh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U1_A;
  }
  if (params->smagorinsky > 0.0 && params->u2vh < 0.0) {
    master->u2vh = master->sdc;
    master->smagcode |= U2_F1;
  } else {
    master->u2vh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U2_A;
  }
  */
  if (params->smagorinsky > 0.0 && params->u1kh < 0.0) {
    int n;
    if (params->smagorinsky == 1.0 && params->kue1 != 1.0) 
      master->smagcode |= U1_SAK;
    for (n = 0; n < params->nobc; ++n)
      if (params->open[n]->sponge_zone_h)
	master->smagcode |= U1_SAK;
    if (master->smagcode & U1_SAK)
      master->u1kh = d_alloc_1d(geom->sgsiz);
    else {
      master->u1kh = master->sdc;
      master->smagcode |= U1_SPK;
    }
  } else {
    master->u1kh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U1_AK;
    if (params->sigma)
      hd_warn("** Smagorinsky e1 diffusion works best with sigma. **\n");
  }
  if (params->smagorinsky > 0.0 && params->u2kh < 0.0) {
    int n;
    if (params->smagorinsky == 1.0 && params->kue2 != 1.0) 
      master->smagcode |= U2_SAK;
    for (n = 0; n < params->nobc; ++n)
      if (params->open[n]->sponge_zone_h)
	master->smagcode |= U2_SAK;
    if (master->smagcode & U2_SAK)
      master->u2kh = d_alloc_1d(geom->sgsiz);
    else {
      master->u2kh = master->sdc;
      master->smagcode |= U2_SPK;
    }
  } else {
    master->u2kh = d_alloc_1d(geom->sgsiz);
    master->smagcode |= U2_AK;
    if (params->sigma)
      hd_warn("** Smagorinsky e2 diffusion works best with sigma. **\n");
  }
  master->t11 = d_alloc_1d(geom->sgsiz);
  master->t22 = d_alloc_1d(geom->sgsiz);
  master->t12 = d_alloc_1d(geom->sgsiz);

  /* Miscillaneous */
  if (params->thin_merge) {
    master->kth_e1 = i_alloc_1d(geom->sgsizS);
    master->kth_e2 = i_alloc_1d(geom->sgsizS);
  }
  if (!(params->means & NONE)) {
    master->meanc = d_alloc_1d(geom->sgsizS);
    memset(master->meanc, 0, geom->sgsizS * sizeof(double));
    if (params->means & TIDAL) {
      master->odeta = d_alloc_1d(geom->sgsizS);
      memset(master->odeta, 0, geom->sgsizS * sizeof(double));
    }
  }
  if (params->runmode & TRANS) {
    /* Use the viscosity tensors to store streamline Courant nos.  */
    master->origin = master->Vz;
    master->pc = master->t11;
    master->qc = master->t12;
    master->rc = master->t22;
  }

  /*UR-FIX out of bounds
   * org
   *  master->wclk = d_alloc_1d(geom->nwindows);
   *
   *  since master->wclk is referenced as master->wclk[nwindows]
   */
  master->wclk = d_alloc_1d(geom->nwindows + 1);
  memset(master->wclk, 0, geom->nwindows * sizeof(double));

  /* Initialize */
  memset(master->wtop, 0, geom->sgsizS * sizeof(double));
  memset(master->wbot, 0, geom->sgsizS * sizeof(double));
  memset(master->waterss, 0, geom->sgsiz * sizeof(double));
  memset(master->waterss2d, 0, geom->sgsizS * sizeof(double));
  memset(master->eta, 0, geom->sgsizS * sizeof(double));
  memset(master->detadt, 0, geom->sgsizS * sizeof(double));
  memset(master->u1av, 0, geom->sgsizS * sizeof(double));
  memset(master->u2av, 0, geom->sgsizS * sizeof(double));
  memset(master->w, 0, geom->sgsiz * sizeof(double));
  memset(master->u1, 0, geom->sgsiz * sizeof(double));
  memset(master->u2, 0, geom->sgsiz * sizeof(double));
  /*
  if (master->eta_rlx)
    memset(master->eta_rlx, 0, geom->sgsizS * sizeof(double));
  */

  /*
   * Initialise s2m_3d tracer maps
   */
  master->trmap_s2m_3d = i_alloc_1d(master->ntr);
  tt = 0;
  for (tn = 0; tn < master->ntr; tn++) {
    if (strncmp(master->trinfo_3d[tn].tag, "DA_", 3)) {
      master->trmap_s2m_3d[tt] = tn;
      tt++;
    }
  }
  master->ntrmap_s2m_3d = tt;

  /* mpi_rank is a global set in main.c */
  master->mpi_rank = mpi_rank;

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
  if (master->d1) d_free_1d(master->d1);
  if (master->d2) d_free_1d(master->d2);
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

      /* Copy the tracer values into a buffer */
      kk = nz;
      while (c != geom->zm1[c]) {
        tr[kk] = master->tr_wc[n][c];
        /* Set values below the 'z' bottom to the last value */
        if ((!tr[kk] || isnan(tr[kk])) && kk < nz)
          tr[kk] = tr[kk + 1];
        kk--;
        c = geom->zm1[c];
      }

      c = c2;
      while (c != geom->zm1[c]) {
        zs = geom->cellz[c] * master->Ds[cs];

        /* Find the z layers bracketing the sigma layer.  */
        /* Note : master->dz contains the original cellz depths as */
        /* read from the input file.  */
        ci = c2;
        kk = nz;
        while (master->dz[ci] > zs && ci != geom->zm1[ci]) {
          ci = geom->zm1[ci];
          kk--;
        }

        /* Do the interpolation */
        if (ci == geom->zp1[ci])
          master->tr_wc[n][c] = tr[nz];
        else if (ci == geom->zm1[ci])
          master->tr_wc[n][c] = tr[0];
        else {
          /* Linear interpolation */
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

      /* Apply a smoothing filter to the vertical profile */
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
void init_wvel_bounds(geometry_t *geom, /* Sparse global geometry
                                           structure */
                      master_t *master  /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coodinate / counter */
  int cs;                       /* 2D sparse coordinate */
  int cb;                       /* Bottom sparse coordinate */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1 */
  double eta_l;                 /* Surface elevation at i-1 */
  double eta_r;                 /* Surface elevation at i+1 */
  double eta_b;                 /* Surface elevation at j-1 */
  double eta_f;                 /* Surface elevation at j+1 */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise */
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                     geom->bin2, 0);

  /*-----------------------------------------------------------------*/
  /* SIGMA : zero velocity at the sigma boundaries */
  if (master->sigma)
    return;

  /* In the linear case wtop and detadt are the same, whereas in the */
  /* non-linear case, they are related by terms involving the */
  /* surface slope and surface horizontal velocities.  */
  memcpy(master->wtop, master->detadt, geom->enonS * sizeof(double));
  if (!master->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    cb = geom->bot_t[cc];
    cs = geom->m2d[c];
    xp1 = geom->xp1[c];
    xm1 = geom->xm1[cs];
    yp1 = geom->yp1[c];
    ym1 = geom->ym1[cs];

    /* Calculate the surface vertical velocity */
    eta_l = master->eta[xm1];
    eta_r = master->eta[geom->xp1[cs]];
    eta_b = master->eta[ym1];
    eta_f = master->eta[geom->yp1[cs]];

    master->detadt[cs] -=
      (master->u1[c] + master->u1[xp1]) * (eta_r -
                                           eta_l) / (4.0 *
                                                     geom->h1acell[cs]) +
      (master->u2[c] + master->u2[yp1]) * (eta_f -
                                           eta_b) / (4.0 *
                                                     geom->h2acell[cs]);

    /* Calculate the bottom vertical velocity */
    master->wbot[cs] = -(master->u1[cb] + master->u1[geom->xp1[cb]]) *
      geom->dHde1[cs] / (2.0 * geom->h1acell[cs]) -
      (master->u2[cb] + master->u2[geom->yp1[cb]]) *
      geom->dHde2[cs] / (2.0 * geom->h2acell[cs]);
  }
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
/* Reset the geometric parameters on common faces across inner       */
/* explicit maps.                                                    */
/*-------------------------------------------------------------------*/
void reset_inner_consts(master_t *master,   /* Master data structure */
			geometry_t *geom    /* Global geometry       */
			)
{
  int cc, cs, cd;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    cs = geom->w2_t[cc];
    if ((cd = geom->emmask[cs])) {
      if (geom->xm1[cs] == cd) {
	inner_mean(geom->u1x, cs, cd);
	inner_mean(geom->u1y, cs, cd);
	inner_mean(geom->h1au1, cs, cd);
	inner_mean(geom->thetau1, cs, cd);
	inner_mean(geom->sinthu1, cs, cd);
	inner_mean(geom->costhu1, cs, cd);
	inner_mean(master->u1c1, cs, cd);
	inner_mean(master->u1c3, cs, cd);
	inner_mean(master->u1c4, cs, cd);
	inner_mean(master->u1c5, cs, cd);
	inner_mean(master->u1c6, cs, cd);
      } else {
	inner_mean(geom->u2x, cs, cd);
	inner_mean(geom->u2y, cs, cd);
	inner_mean(geom->h2au2, cs, cd);
	inner_mean(geom->thetau2, cs, cd);
	inner_mean(geom->sinthu2, cs, cd);
	inner_mean(geom->costhu2, cs, cd);
	inner_mean(master->u2c1, cs, cd);
	inner_mean(master->u2c3, cs, cd);
	inner_mean(master->u2c4, cs, cd);
	inner_mean(master->u2c5, cs, cd);
	inner_mean(master->u2c6, cs, cd);
      }
    }
  }
}
/*
void inner_mean(double *var, int cs, int cd)
{
  double vs = var[cs];
  double vd = var[cd];

  if ((vs > 0 && vd > 0) || (vs < 0 && vd < 0))
    var[cs] = var[cd] = 0.5 * (vs + vd);
  else
    inner_mean_neg(var, cs, cd);
}
*/
void inner_mean(double *var, int cs, int cd)
{
  double val;
  val = 0.5 * (fabs(var[cs]) + fabs(var[cd]));
  var[cs] = (var[cs] > 0) ? val : -val;
  var[cd] = (var[cd] > 0) ? val : -val;
}

/* END reset_inner_consts()                                          */
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
  int xp1, xm1, yp1, ym1;
  double **fdata = NULL;
  char buf[MAXSTRLEN];
  char *fname[8] = {"fetch_n","fetch_ne","fetch_e","fetch_se",
		    "fetch_s","fetch_sw","fetch_w","fetch_nw"};

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

  /* Loop over the water cells */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (n = 0; n < 8; n++)
      master->fetch[c][n] = 0.0;
    if(db)printf("start %d\n",c);
    /* North */
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
/* transport mode. Note: these data may be read from a different     */
/* grid to that which the transport model is run on.                 */
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
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  int c, i, j, k, m, n, cb, ci, cc;
  int zp1, zm1;
  size_t start[4];
  size_t count[4];
  int fid;
  int ncerr;
  size_t nz, kz;
  size_t nce1, nce2;
  size_t nfe1, nfe2;
  double **bathy;
  double *layers;
  double *cellz;
  double **d1, **d2;
  double **u1x, **u1y;
  double **u2x, **u2y;
  double **cellx, **celly;
  double **gridx, **gridy;
  double **sx, **sy;
  int nobc;       /* No. boundaries in source grid */
  int *npts;      /* No. cells in each boundary */
  int **iloc;     /* Boundary x location list   */
  int **jloc;     /* Boundary y location list   */
  int *type;      /* Boundary orientation       */
  double w1, w2;
  short **mask;
  int xm1, xp1, yp1, ym1;

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
  if (strcmp(params->idumpname, params->sourcefile) == 0) {
    set_tp(geom, master, window);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Read the dimensions of the velocity grid                        */
  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(params->sourcefile, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("Can't find transport source grid file %s\n", params->sourcefile);
    hd_quit((char *)nc_strerror(ncerr));
  }

  /* Get dimensions                                                  */
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_grid"), &kz);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "k_centre"), &nz);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_centre"), &nce2);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_centre"), &nce1);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "i_left"), &nfe1);
  nc_inq_dimlen(fid, ncw_dim_id(fid, "j_back"), &nfe2);

  /* Read the vertical layer structure */
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = kz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  layers =  d_alloc_1d(kz);
  nc_get_vara_double(fid, ncw_var_id(fid, "z_grid"), start, count, layers);
  count[0] = nz;
  cellz =  d_alloc_1d(nz);
  nc_get_vara_double(fid, ncw_var_id(fid, "z_centre"), start, count, cellz);

  /* Read the bathymetry of the source grid */
  count[0] = nce2;
  count[1] = nce1;
  bathy =  d_alloc_2d(nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "botz"), start, count, bathy[0]);

  /*-----------------------------------------------------------------*/
  /* Read open boundary information for the source grid              */
  nobc = m = 0;
  prm_read_int(params->prmfd, "SBOUNDARIES", &nobc);
  if (nobc > 0) {
    open_bdrys_t **open;
    open = (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * nobc);      
    npts = i_alloc_1d(nobc);
    type = i_alloc_1d(nobc);
    for (n = 0; n < nobc; n++) {
      open[n] = OBC_alloc();
      sprintf(keyword, "SBOUNDARY%1d.TYPE", n);
      prm_read_char(params->prmfd, keyword, buf);
      if (strcmp(buf, "u1") == 0)
	open[n]->type = U1BDRY;
      if (strcmp(buf, "u2") == 0)
	open[n]->type = U2BDRY;
      get_obc_list(open[n], params->prmfd, n, "SBOUNDARY");
      type[n] = open[n]->type;
      npts[n] = open[n]->npts;
      if (npts[n] > m) m = npts[n];
    }
    iloc = i_alloc_2d(m, nobc);
    jloc = i_alloc_2d(m, nobc);
    for (n = 0; n < nobc; n++) {
      for (m = 0; m < npts[n]; m++) {
	iloc[n][m] = open[n]->iloc[m];
	jloc[n][m] = open[n]->jloc[m];
      }
      free((open_bdrys_t *)open[n]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Make sparse maps for the velocity grid */
  build_sparse_grid(params, tpg, nce1, nce2, nz, bathy, layers, nobc,
		    npts, iloc, jloc, type);
  if (nobc) {
    i_free_1d(npts);
    i_free_1d(type);
    i_free_2d(iloc);
    i_free_2d(jloc);
  }
  hd_warn("Number of wet cells in source grid: 2D = %d, 3D = %d\n", 
	  tpg->b2_t, tpg->b3_t);
  hd_warn("Number of wet cells in target grid: 2D = %d, 3D = %d\n", 
	  geom->b2_t, geom->b3_t);

  /* Allocate memory for the transport variables */
  tp->ns2 = tpg->ns2;
  tp->ns3 = tpg->ns3;
  /* Source grid 3D */
  n = tpg->enon + 1;
  tpd->u1 = d_alloc_1d(n);
  tpd->u2 = d_alloc_1d(n);
  tpd->w = d_alloc_1d(n);
  if (master->tmode & SP_FFSL) {
    tpd->u1vm = d_alloc_1d(n);
    tpd->u2vm = d_alloc_1d(n);
  }
    tpc->nu = d_alloc_1d(n);
    tpc->nv = d_alloc_1d(n);
    tpc->nw = d_alloc_1d(n);
    tpc->clxf = i_alloc_1d(n);
    tpc->clyf = i_alloc_1d(n);
    tpc->clzf = i_alloc_1d(n);
    tpc->clxc = i_alloc_1d(n);
    tpc->clyc = i_alloc_1d(n);
    tpc->clzc = i_alloc_1d(n);
    tpc->crfxf = d_alloc_1d(n);
    tpc->crfyf = d_alloc_1d(n);
    tpc->crfzf = d_alloc_1d(n);
    tpc->crfxc = d_alloc_1d(n);
    tpc->crfyc = d_alloc_1d(n);
    tpc->crfzc = d_alloc_1d(n);
    tpc->tr_mod = d_alloc_1d(n);
    tpc->tr_mod_x = d_alloc_1d(n);
    tpc->tr_mod_y = d_alloc_1d(n);
    tpc->tr_mod_z = d_alloc_1d(n);
  tpd->Kz = d_alloc_1d(n);
  tpg->cellz = d_alloc_1d(n);
  tpg->gridz = d_alloc_1d(n);
  tpg->mgm = i_alloc_1d(n);
  tpc->dz = d_alloc_1d(n);
  tpc->w4 = d_alloc_1d(n);
  tpc->w5 = d_alloc_1d(n);
  tpc->w6 = d_alloc_1d(n);
  tpc->w7 = d_alloc_1d(n);
  tpc->w8 = d_alloc_1d(n);
  tpc->w9 = d_alloc_1d(n);
  tpc->w10 = d_alloc_1d(n);
  tpc->s1 = i_alloc_1d(n);
  tpc->s2 = i_alloc_1d(n);
  tpc->s3 = i_alloc_1d(n);
  tpc->m2d = i_alloc_1d(n);
  tpc->c1 = c_alloc_1d(n);

  /* Source grid tracers */
  tpd->ntr = 2 + master->ntrvars; /* temp & salt + additional tracers */
  tpd->tr_wc = d_alloc_2d(n, tpd->ntr);
  tpd->ntrS = master->ntrvarsS;
  if (tpd->ntrS)
    tpd->tr_wcS = d_alloc_2d(n, tpd->ntrS);
  /* Tracer 0 is salt */
  tpd->sal = tpd->tr_wc[0];
  tpd->sno = 0;
  /* Tracer 0 is temp */
  tpd->temp = tpd->tr_wc[1];
  tpd->tno = 1;

  /* Source grid 2D */
  n = tpg->enonS + 1;
  tpg->botz = d_alloc_1d(n);
  tpg->cellx = d_alloc_1d(n);
  tpg->celly = d_alloc_1d(n);
  tpg->gridx = d_alloc_1d(n);
  tpg->gridy = d_alloc_1d(n);
  tpg->h1au1 = d_alloc_1d(n);
  tpg->h2au2 = d_alloc_1d(n);
  tpg->h1au2 = d_alloc_1d(n);
  tpg->h2au1 = d_alloc_1d(n);
  tpg->h1acell = d_alloc_1d(n);
  tpg->h2acell = d_alloc_1d(n);
  tpg->thetau1 = d_alloc_1d(n);
  tpg->thetau2 = d_alloc_1d(n);
  tpg->sinthcell = d_alloc_1d(n);
  tpg->costhcell = d_alloc_1d(n);
  tpg->nsur_t = i_alloc_1d(n);
  tpg->dHde1 = d_alloc_1d(n);
  tpg->dHde2 = d_alloc_1d(n);
  tpd->eta = d_alloc_1d(n);
  tpd->etab = d_alloc_1d(n);
  tpd->detadt = d_alloc_1d(n);
  tpd->wtop = d_alloc_1d(n);
  tpd->wbot = d_alloc_1d(n);
  tpd->light = d_alloc_1d(n);
  tpc->d1 = d_alloc_1d(n);
  tpc->d2 = d_alloc_1d(n);
  tpc->d3 = d_alloc_1d(n);
  tpc->d4 = d_alloc_1d(n);
  tpc->d7 = d_alloc_1d(n);
  tpc->d8 = d_alloc_1d(n);
  tpc->oldeta = d_alloc_1d(n);
  tpc->i1 = i_alloc_1d(n);
  tpc->i2 = i_alloc_1d(n);
  tpc->i3 = i_alloc_1d(n);
  tpc->cdry_e1 = i_alloc_1d(n);

  /* Target grid 3D */
  n = geom->sgsizS;  
  tpc->wgt = d_alloc_2d(8, n); /* Not required */

  tpc->velmax = master->velmax;
  tpc->etamax = master->etamax;
  tp->w1 = d_alloc_1d(tpg->sgsiz);
  memcpy(tpg->nsur_t, tpg->sur_t, (tpg->b2_t + 1) * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Read the cell information off the source grid                   */
  for (c = 1; c <= tpg->sgnum; c++) {
    i = tpg->s2i[c];
    j = tpg->s2j[c];
    k = tpg->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    tpg->cellz[c] = cellz[k];
    tpg->gridz[c] = layers[k];
  }
  c2s_2d(tpg, tpg->botz, bathy, nce1, nce2);

  /* Get the grid spacings at cell faces                             */
  d1 = d_alloc_2d(nfe1, nce2);
  count[0] = nce2;
  count[1] = nfe1;
  nc_get_vara_double(fid, ncw_var_id(fid, "h1au1"), start, count, d1[0]);
  c2s_2d(tpg, tpg->h1au1, d1, nfe1, nce2);                     
  d_free_2d(d1);
  d1 = d_alloc_2d(nce1, nfe2);
  count[0] = nfe2;
  count[1] = nce1;
  nc_get_vara_double(fid, ncw_var_id(fid, "h2au2"), start, count, d1[0]);
  c2s_2d(tpg, tpg->h2au2, d1, nce1, nfe2);                     
  d_free_2d(d1);

  /*-----------------------------------------------------------------*/
  /* Set up the source grid mesh (with origin in the top-right cell  */
  /* corner) for the xyij_tree.                                      */
  sx = d_alloc_2d(nfe1+1, nfe2+1);
  sy = d_alloc_2d(nfe1+1, nfe2+1);
  /* Use (u2x,u2y) for the top and bottom of the grid, i = 1:nce1    */
  u2x = d_alloc_2d(nce1, nfe2);
  u2y = d_alloc_2d(nce1, nfe2);
  count[0] = nfe2;
  count[1] = nce1;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_back"), start, count, u2x[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_back"), start, count, u2y[0]);
  /*
  d_free_2d(d1);
  d_free_2d(d2);
  */
  /* Use (u1x,u1y) for the left and right of the grid, j = 1:nce2    */
  u1x = d_alloc_2d(nfe1, nce2);
  u1y = d_alloc_2d(nfe1, nce2);
  count[0] = nce2;
  count[1] = nfe1;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_left"), start, count, u1x[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_left"), start, count, u1y[0]);
  /*
  d_free_2d(d1);
  d_free_2d(d2);
  */
  /* Use (gridx,gridy) for the corners                              */
  gridx = d_alloc_2d(nfe1, nfe2);
  gridy = d_alloc_2d(nfe1, nfe2);
  count[0] = nfe2;
  count[1] = nfe1;
  nc_get_vara_double(fid, ncw_var_id(fid, "x_grid"), start, count, gridx[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_grid"), start, count, gridy[0]);

  /* Use (cellx,celly) for all other cells                           */
  d1 = d_alloc_2d(nce1, nce2);
  d2 = d_alloc_2d(nce1, nce2);
  cellx = d_alloc_2d(nce1, nce2);
  celly = d_alloc_2d(nce1, nce2);
  count[0] = nce2;
  count[1] = nce1;
  nc_get_vara_double(fid, ncw_var_id(fid, "h1acell"), start, count, d1[0]);
  c2s_2d(tpg, tpg->h1acell, d1, nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "h2acell"), start, count, d2[0]);
  c2s_2d(tpg, tpg->h2acell, d2, nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "x_centre"), start, count, cellx[0]);
  c2s_2d(tpg, tpg->cellx, cellx, nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "y_centre"), start, count, celly[0]);
  c2s_2d(tpg, tpg->celly, celly, nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "thetau1"), start, count, d1[0]);
  c2s_2d(tpg, tpg->thetau1, d1, nce1, nce2);
  nc_get_vara_double(fid, ncw_var_id(fid, "thetau2"), start, count, d2[0]);
  c2s_2d(tpg, tpg->thetau2, d2, nce1, nce2);

  /* Generate the source mesh and allocate the undefined cells */
  set_trans_grid(nce1, nce2, sx, sy, cellx, celly, gridx, gridy, u1x, u1y, u2x, u2y, tpg);

  /* Generate the source mesh tree */
  if(!(master->tmode & INEXACT))
    tp->xyij_tree = grid_xytoij_init(sx, sy, nfe1, nfe2);
  /*tpg->xyij_tree = tp->xyij_tree;*/
  tpg->xyij_tree = grid_xytoij_init(gridx, gridy, tpg->nce1, tpg->nce2);

  /*
  c = geom->map[nz-1][27][17];
  printf("a %d\n",c);
  c = xytoc(geom, master->xyij_tree, geom->u2x[c], geom->u2y[c], &w1, &w2);
  printf("a %d (%d %d) : %f %f (%d %d)\n",c,geom->s2i[c],geom->s2j[c],w1,w2,(int)w1,(int)w2);
  grid_xytofij(master->xyij_tree, geom->u2x[c], geom->u2y[c], &w1, &w2);
  printf("a %f %f (%d %d) : (%d %d)\n",w1,w2,(int)w1,(int)w2,i,j);
  grid_xytofij(tp->xyij_tree, sx[28][19], sy[28][19], &w1, &w2);
  c = tpg->map[nz-1][(int)w2][(int)w1];
  printf("%d (%f %f) (%d %d)\n",c,w1,w2,tpg->s2i[c],tpg->s2j[c]);
  printf("%f %f\n",cellx[27][18],celly[27][18]);
  */

  /*-----------------------------------------------------------------*/
  /* Reset cellz to account for the bottom                           */
  for (n = 1; n <= geom->nwindows; n++) {
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      cb = window[n]->bot_t[cc];
      zp1 = window[n]->zp1[cb];
      zm1 = window[n]->zm1[cb];
      window[n]->cellz[cb] = window[n]->botz[c] + 
	0.5 * (window[n]->gridz[zp1] - window[n]->botz[c]);
      window[n]->cellz[zm1] = window[n]->gridz[cb]; /* Sediment cell */
    }
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    cb = geom->bot_t[cc];
    zp1 = geom->zp1[cb];
    zm1 = geom->zm1[cb];
    geom->cellz[cb] = geom->botz[c] + 
      0.5 * (geom->gridz[zp1] - geom->botz[c]);
    geom->cellz[zm1] = geom->gridz[cb];
  }
  for (cc = 1; cc <= tpg->b2_t; cc++) {
    c = tpg->w2_t[cc];
    cb = tpg->bot_t[cc];
    zp1 = tpg->zp1[cb];
    zm1 = tpg->zm1[cb];
    tpg->cellz[cb] = tpg->botz[c] + 
      0.5 * (tpg->gridz[zp1] - tpg->botz[c]);
    tpg->cellz[zm1] = tpg->gridz[cb];
  }

  /*-----------------------------------------------------------------*/
  /* Set a no-gradient condition over lateral boundaries             */
  s2ij_setghosts(tpg);
  for (cc = 1; cc <= geom->nbpt; cc++) {
    cb = geom->bpt[cc];
    ci = geom->bin[cc];
    geom->cellz[cb] = geom->cellz[ci];
  }
  for (cc = 1; cc <= tpg->nbpt; cc++) {
    cb = tpg->bpt[cc];
    ci = tpg->bin[cc];
    tpg->cellz[cb] = tpg->cellz[ci];
    tpg->gridz[cb] = tpg->gridz[ci];
    if (cc <= tpg->nbptS) {
      tpg->h1au1[cb] = tpg->h1au1[ci];
      tpg->h2au2[cb] = tpg->h2au2[ci];
      tpg->botz[cb] = tpg->botz[ci];
    }
  }

  /* Set the bottom depth on e1 and e2 open boundaries (no gradient) */
  for (n = 0; n < tpg->nobc; n++) {
    open_bdrys_t *open = tpg->open[n];
    open->stagger = OUTFACE;
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no2_e1; cc++) {
        c = open->obc_e1[cc];
        ci = open->oi1_e1[cc];
        tpg->botz[c] = tpg->botz[ci];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no2_e2; cc++) {
        c = open->obc_e2[cc];
        ci = open->oi1_e2[cc];
        tpg->botz[c] = tpg->botz[ci];
      }
    }
  }

  /* Get the source - target map */
  set_map_l(tpg);
  tpc->tgrid = get_st_map(master, geom, tpg, sx, sy, gridx, gridy, tpg->mgm, tpc->d7, tpc->d8);
  for (n = 1; n <= geom->nwindows; n++) {
    wincon[n]->tgrid = tpc->tgrid;
    wincon[n]->d8 = d_alloc_1d(window[n]->sgsizS);
  }

  /* Get the target - source map */
  get_ts_map(geom, window, wincon, tpg, tpc, gridx, gridy, master->togn);

  /* Get the sin and cos at cell centres */
  for (cc = 1; cc <= tpg->b2_t; cc++) {
    c = tpg->w2_t[cc];
    xp1 = tpg->xp1[c];
    xm1 = tpg->xm1[c];
    yp1 = tpg->yp1[c];
    ym1 = tpg->ym1[c];

    tpg->sinthcell[c] = (sin(tpg->thetau1[c]) +
			 sin(tpg->thetau2[c])) / 2.0;
    if (!isnan(tpg->thetau1[xp1]) && tpg->thetau1[xp1] != 0.)
      tpg->sinthcell[c] = (tpg->sinthcell[c] + sin(tpg->thetau1[xp1])) / 2.0;
    if (!isnan(tpg->thetau2[yp1]) && tpg->thetau2[yp1] != 0.)
      tpg->sinthcell[c] = (tpg->sinthcell[c] + sin(tpg->thetau2[yp1])) / 2.0;
    tpg->costhcell[c] = (cos(tpg->thetau1[c]) +
			 cos(tpg->thetau2[c])) / 2.0;
    if (!isnan(tpg->thetau1[xp1]) && tpg->thetau1[xp1] != 0.)
      tpg->costhcell[c] = (tpg->costhcell[c] + cos(tpg->thetau1[xp1])) / 2.0;
    if (!isnan(tpg->thetau2[yp1]) && tpg->thetau2[yp1] != 0.)
      tpg->costhcell[c] = (tpg->costhcell[c] + cos(tpg->thetau2[yp1])) / 2.0;
  }

  /* Sanity checks */
  if (tpc->tgrid == EXACT) {
    for (n = 1; n <= geom->nwindows; n++) {
      double *xinit = wincon[n]->d7;
      double *yinit = wincon[n]->d8;
      if (fabs(xinit[c]) < EPS) xinit[c] = 0.0;
      if (fabs(yinit[c]) < EPS) yinit[c] = 0.0;
      for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (xinit[c] != 0.0 || yinit[c] != 0.0)
	hd_quit("Can't define source start trajectory (%d %d) = (%e %e)\n",
		window[n]->s2i[c], window[n]->s2j[c], xinit[c], yinit[c]);
      if (c != window[n]->mgm[c])
	hd_quit("Target to source maps differ t(%d %d) => s(%d)\n",
		window[n]->s2i[c], window[n]->s2j[c], window[n]->mgm[c]);
      }
    }
    for (cc = 1; cc <= tpg->b2_t; cc++) {
      c = tpg->w2_t[cc];
      i = tpg->s2i[c];
      j = tpg->s2j[c];
      if (tpc->s2[c]) continue;
      if (c != tpg->mgm[c])
	hd_quit("Source to target maps differ s(%d %d) => t(%d)\n", i, j, tpg->mgm[c]);
      if (tpc->d7[c] != 0.0 || tpc->d8[c] != 0.0)
	hd_quit("Can't define source cell offset in corresponding target cell (%d %d) = (%e %e)\n",
		i, j, tpc->d7[c], tpc->d8[c]);
    }
  }

  d_free_2d(d1);
  d_free_2d(d2);
  /*
  d_free_2d(sx);
  d_free_2d(sy);
  */
  d_free_2d(u1x);
  d_free_2d(u1y);
  d_free_2d(u2x);
  d_free_2d(u2y);
  /*
  d_free_2d(gridx);
  d_free_2d(gridy);
  */
  d_free_2d(cellx);
  d_free_2d(celly);

  set_dz(tpg, tpd, tpc);
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->trans = geom->trans;

  d_free_1d(layers);
  d_free_1d(cellz);
  d_free_2d(bathy);
}

/* END trans_vel_init()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find the cells to process and offsets to apply to the  */
/* source grid.                                                      */
/*-------------------------------------------------------------------*/
void get_ts_map(geometry_t *geom, geometry_t **window, win_priv_t **wincon,
		geometry_t *ws, win_priv_t *wc, double **gridx, double **gridy, 
		int of)
{
  xytoij_tree_t *s_xyij_tree;
  int c, cs, c2, c1, cc, n, i, j, clast;
  double x, y, z, d1, d2, top, bot;
  double *dzz = wc->w9;
  int bl_flag;
  double eps = 1.e-9;

  /* Get the tree for the source grid                                */
  /*s_xyij_tree = grid_xytoij_init(gridx, gridy, ws->nce1, ws->nce2);*/
  s_xyij_tree = ws->xyij_tree; 
  bl_flag = (of & BOTMLEFT) ? 1 : 0;

  for (n = 1; n <= geom->nwindows; n++) {
    double *xinit = wincon[n]->d7;
    double *yinit = wincon[n]->d8;
    window[n]->mgm = i_alloc_1d(window[n]->sgsiz);
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      /* Get the wet points in the source grid                       */
      c = c2 = window[n]->w2_t[cc];
      /* The source grid tree is supplied to xytoc() and the cell    */
      /* centre, is returned. The quadrant that (cellx,celly) falls  */
      /* into (relative to the source grid cell centre) must be      */
      /* located.                                                    */
      cs = xytoc(ws, s_xyij_tree, window[n]->cellx[c], window[n]->celly[c], &x, &y);
      i = (int)x; j = (int)y;
      /* Find the closest wet cell if c maps onto land               */
      if (cs == 0) {
	int in, jn, k = ws->nz - 1;
	in = (int)floor(x);
	jn = (int)floor(y);
	if (in < 0 || in > ws->nce1 || jn < 0 || jn > ws->nce2)
	  hd_quit("Can't map target-source for target %d(%d %d)\n",
		  c, window[n]->s2i[c], window[n]->s2j[c]);
	/*cs = find_closest_w(ws, in, jn, &i, &j);*/
	cs = find_nearest(ws, window[n]->cellx[c2], window[n]->celly[c2]);
	if (cs == 0)
	  hd_quit("Can't map target@%d(%d %d) to source@(%d %d)\n",
		  c,window[n]->s2i[c], window[n]->s2j[c], i, j);
	else
	  hd_warn("get_ts_map: target (%d %d) maps to source land: using source (%d %d)\n",
		  window[n]->s2i[c], window[n]->s2j[c], ws->s2i[cs], ws->s2j[cs]);
	x = y = 0.5;
	i = j = 0;
      }
      /* Locate the quadrant and set the t2s map to the source mesh  */
      /* origin. Note xinit and yinit must be > 0 and < 1.           */
      if (bl_flag) {
	d1 = (x - (double)i) - 0.5;
	d2 = (y - (double)j) - 0.5;
      } else {
	d1 = 0.5 - (x - (double)i);
	d2 = 0.5 - (y - (double)j);
      }
      d1 = (fabs(d1) < eps) ? 0.0 : d1;   /* Precision errors in xytoc */
      d2 = (fabs(d2) < eps) ? 0.0 : d2;
      if (fabs(d1) < 2.0*EPS && fabs(d2) < 2.0*EPS) { /* Centre      */
	xinit[c] = 0.0;
	yinit[c] = 0.0;
      } else if (d1 >= 0.0 && d2 >= 0.0) {       /* SW quadrant (tr) */
	xinit[c] = d1;                           /* NE quadrant (bl) */
	yinit[c] = d2;
      } else if (d1 >= 0.0 && d2 < 0.0) {        /* NW quadrant (tr) */
	if (bl_flag)                             /* SE quadrant (bl) */
	  cs = ws->ym1[cs];
	else
	  cs = ws->yp1[cs];
	xinit[c] = d1;
	yinit[c] = 1.0 + d2;
      } else if (d1 < 0.0 && d2 < 0.0) {        /* NE quadrant (tr) */
	if (bl_flag)                            /* SW quadrant (bl) */
	  cs = ws->xm1[ws->ym1[cs]];
	else
	  cs = ws->xp1[ws->yp1[cs]];
	xinit[c] = 1.0 + d1;
	yinit[c] = 1.0 + d2;
      } else if (d1 < 0.0 && d2 >= 0.0) {       /* SE quadrant (tr) */
	if (bl_flag)                            /* SW quadrant (bl) */
	  cs = ws->xm1[cs];
	else
	  cs = ws->xp1[cs];
	xinit[c] = 1.0 + d1;
	yinit[c] = d2;
      }
      window[n]->mgm[c] = cs;
      clast = cs;

      /* Get the t2s map for all surface and sub-surface cells       */
      top = 0.0;
      while (c != window[n]->zm1[c]) {
	bot = max(window[n]->gridz[c], window[n]->botz[c2]);
	z = 0.5 * (top + bot);
	c1 = ztoc(ws, cs, z);
	/*if (z < ws->cellz[c1]) c1 = ws->zm1[c1];*/
	if (z < ws->botz[ws->m2d[c1]]) {
	  hd_warn("get_ts_map: target cellz %f %d(%d %d %d) is below source bottom %f at %d(%d %d %d): mapped to source %d\n",
		  z,c,window[n]->s2i[c],window[n]->s2j[c],window[n]->s2k[c],
		  ws->botz[ws->m2d[c1]],c1,ws->s2i[c1],ws->s2j[c1],ws->s2k[c1], clast);
	  c1 = clast;
	}
	window[n]->mgm[c] = c1;
	c = window[n]->zm1[c];
	top = bot;
	clast = c1;
      }
    }
  }
  /*free((xytoij_tree_t *)s_xyij_tree);*/
}

/* END get_ts_map()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate grid location at undefined cells           */
/*-------------------------------------------------------------------*/
void set_trans_grid(int nce1, int nce2, double **sx, double **sy, 
		    double **cellx, double **celly, 
		    double **gridx, double **gridy,
		    double **u1x, double **u1y, double **u2x, double **u2y,
		    geometry_t *geom)
		
{
  int i, j, ii, jj, c, cc;
  int nfe1 = nce1 + 1;
  int nfe2 = nce2 + 1;
  short **mask;

  mask = s_alloc_2d(nfe1+1,nfe2+1);
  for (j = 0; j <= nfe2; j++)
    for (i = 0; i <= nfe1; i++) {
      mask[j][i] = 0;
    }

  /* Use (u2x,u2y) for the top and bottom of the grid, i = 1:nce1    */
  for (i = 1; i < nfe1; i++) {
    sx[0][i] = u2x[0][i-1];
    sy[0][i] = u2y[0][i-1];
    sx[nfe2][i] = u2x[nce2][i-1];
    sy[nfe2][i] = u2y[nce2][i-1];
    if (!isnan(u2x[0][i-1]) && !isnan(u2y[0][i-1]))mask[0][i] = 1;
    if (!isnan(u2x[nce2][i-1]) && !isnan(u2y[nce2][i-1]))mask[nfe2][i] = 1;
  }
  /* Use (u1x,u1y) for the left and right of the grid, j = 1:nce2    */
  for (j = 1; j < nfe2; j++) {
    sx[j][0] = u1x[j-1][0];
    sy[j][0] = u1y[j-1][0];
    sx[j][nfe1] = u1x[j-1][nce1];
    sy[j][nfe1] = u1y[j-1][nce1];
    if (!isnan(u1x[j-1][0]) && !isnan(u1y[j-1][0]))mask[j][0] = 1;
    if (!isnan(u1x[j-1][nce1]) && !isnan(u1y[j-1][nce1]))mask[j][nfe1] = 1;
  }
  /* Use (gridx,gridy) for the corners                              */
  sx[0][0] = gridx[0][0]; sy[0][0] = gridy[0][0];
  sx[0][nfe1] = gridx[0][nce1]; sy[0][nfe1] = gridy[0][nce1];
  sx[nfe2][0] = gridx[nce2][0]; sy[nfe2][0] = gridy[nce2][0];
  sx[nfe2][nfe1] = gridx[nce2][nce1]; sy[nfe2][nfe1] = gridy[nce2][nce1];
  if (!isnan(gridx[0][0]) && !isnan(gridy[0][0]))mask[0][0] = 1;
  if (!isnan(gridx[0][nce1]) && !isnan(gridy[0][nce1]))mask[0][nfe1] = 1;
  if (!isnan(gridx[nce2][0]) && !isnan(gridy[nce2][0]))mask[nfe2][0] = 1;
  if (!isnan(gridx[nce2][nce1]) && !isnan(gridy[nce2][nce1]))mask[nfe2][nfe1] = 1;
  /* Use (cellx,celly) for all other cells                           */
  for (j = 1; j < nfe2; j++)
    for (i = 1; i < nfe1; i++) {
      sx[j][i] = cellx[j-1][i-1];
      sy[j][i] = celly[j-1][i-1];
      if(!isnan(cellx[j-1][i-1]) && !isnan(celly[j-1][i-1])) mask[j][i] = 1;
    }
  /* Interpolate over undefined cells                                */
  for (j = 1; j < nfe2; j++) {
    for (i = 1; i < nfe1; i++) {
      ii = i - 1;
      jj = j - 1;
      if (isnan(sx[j][i])) {
	if (mask[j][i-1])
	  sx[j][i] = intpf(cellx[jj][ii-1], u1x[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j][i+1])
	  sx[j][i] = intpf(cellx[jj][ii+1], u1x[jj][ii+1], 0.0, -0.5, -1.0);
	else if (mask[j-1][i])
	  sx[j][i] = intpf(cellx[jj-1][ii], u2x[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j+1][i])
	  sx[j][i] = intpf(cellx[jj+1][ii], u2x[jj+1][ii], 0.0, -0.5, -1.0);
	else if (mask[j-1][i-1])
	  sx[j][i] = intpf(cellx[jj-1][ii-1], gridx[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j+1][i+1])
	  sx[j][i] = intpf(cellx[jj+1][ii+1], gridx[jj+1][ii+1], 0.0, -0.5, -1.0);
	else if (mask[j+1][i-1]) {
	  sx[j][i] = intpf(cellx[jj+1][ii-1], gridx[jj+1][ii], 0.0, 0.5, 1.0);
	}
	else if (mask[j-1][i+1]) {
	  sx[j][i] = intpf(cellx[jj-1][ii+1], gridx[jj][ii+1], 0.0, -0.5, -1.0);
	}
      }
    }
  }
  for (j = 1; j < nfe2; j++) {
    for (i = 1; i < nfe1; i++) {
      ii = i - 1;
      jj = j - 1;
      if (isnan(sy[j][i])) {
	if (mask[j][i-1])
	  sy[j][i] = intpf(celly[jj][ii-1], u1y[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j][i+1])
	  sy[j][i] = intpf(celly[jj][ii+1], u1y[jj][ii+1], 0.0, -0.5, -1.0);
	else if (mask[j-1][i])
	  sy[j][i] = intpf(celly[jj-1][ii], u2y[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j+1][i])
	  sy[j][i] = intpf(celly[jj+1][ii], u2y[jj+1][ii], 0.0, -0.5, -1.0);
	else if (mask[j-1][i-1])
	  sy[j][i] = intpf(celly[jj-1][ii-1], gridy[jj][ii], 0.0, 0.5, 1.0);
	else if (mask[j+1][i+1])
	  sy[j][i] = intpf(celly[jj+1][ii+1], gridy[jj+1][ii+1], 0.0, -0.5, -1.0);
	else if (mask[j+1][i-1])
	  sy[j][i] = intpf(celly[jj+1][ii-1], gridy[jj+1][ii], 0.0, 0.5, 1.0);
	else if (mask[j-1][i+1])
	  sy[j][i] = intpf(celly[jj-1][ii+1], gridy[jj][ii+1], 0.0, -0.5, -1.0);
      }
    }
  }

  for (j = 0; j <= nfe2; j++) {
    for (i = 0; i <= nfe1; i++) {
      int ci, cj;
      if(isnan(sx[j][i])) {
	find_non_nan(sx, nfe1, nfe2, i, j, &ci, &cj);
	sx[j][i] = sx[cj][ci];
      }
      if(isnan(sy[j][i])) {
	find_non_nan(sy, nfe1, nfe2, i, j, &ci, &cj);
	sy[j][i] = sy[cj][ci];
      }
    }
  }

  /* Check for nans                                                  */
  for (cc = 1; cc <= geom->v2_t; cc++) {
    c = geom->w2_t[cc];
    ii = geom->s2i[c];
    jj = geom->s2j[c];
    i = ii+1;
    j = jj+1;
    /* Check the south-west quadrant */
    if(isnan(sx[j][i]) || isnan(sy[j][i])) {
      printf("SW NaN found at 1(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j][i-1]) || isnan(sy[j][i-1])) {
      printf("SW NaN found at 2(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j-1][i]) || isnan(sy[j-1][i])) {
      printf("SW NaN found at 3(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j-1][i-1]) || isnan(sy[j-1][i-1])) {
      printf("SW NaN found at 4(%d %d)\n",ii,jj);
    }
    /* Check the north-west quadrant */
    if(isnan(sx[j][i-1]) || isnan(sy[j][i-1])) {
      printf("NW NaN found at 2(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j+1][i]) || isnan(sy[j+1][i])) {
      printf("NW NaN found at 3(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j+1][i-1]) || isnan(sy[j+1][i-1])) {
      printf("NW NaN found at 4(%d %d)\n",ii,jj);    
    }
    /* Check the north-east quadrant */
    if(isnan(sx[j][i+1]) || isnan(sy[j][i+1])) {
      printf("NE NaN found at 2(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j+1][i]) || isnan(sy[j+1][i])) {
      printf("NE NaN found at 3(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j+1][i+1]) || isnan(sy[j+1][i+1])) {
      printf("NE NaN found at 4(%d %d) %f %f\n",ii,jj,sx[j+1][i+1],sy[j+1][i+1]);
    }
    /* Check the south-east quadrant */
    if(isnan(sx[j][i+1]) || isnan(sy[j][i+1])) {
      printf("SE NaN found at 2(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j-1][i]) || isnan(sy[j-1][i])) {
      printf("SE NaN found at 3(%d %d)\n",ii,jj);    
    }
    if(isnan(sx[j-1][i+1]) || isnan(sy[j-1][i+1])) {
      printf("SE NaN found at 4(%d %d)\n",ii,jj);    
    }
  }
  s_free_2d(mask);
}

/* END set_trans_grid()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the mapping between source and target grids if  */
/* the source grid is a subset of the target grid.                   */
/*-------------------------------------------------------------------*/
int get_st_map(master_t *master, geometry_t *wt, geometry_t *ws, double **sx, double **sy,
	       double **gridx, double **gridy, int *s2t, double *xoset, double *yoset)
{
  double **d1, **d2;
  double **tx, **ty;
  xytoij_tree_t *t_xyij_tree;
  int nce1, nce2, nfe1, nfe2;
  int i, j, k, c, cc, ii, jj, cs, ct, c1, c2, clast;
  double w1, w2;
  double z, top, bot;
  double x, y;
  double xo, xe, xn, xne;
  double yo, ye, yn, yne;
  int co, ce, cn, cne;
  int io, ie, in, ine;
  int jo, je, jn, jne;
  int org, tgrid;
  int *mask;

  nce1 = wt->nce1;
  nce2 = wt->nce2;
  nfe1 = wt->nfe1;
  nfe2 = wt->nfe2;

  /* The preprocessor includes open boundary cells on F_EDGE and     */
  /* R_EDGE in the interior of the domain (i.e. OUTSIDE cells in the */
  /* domain interior) in the arrays b2_t and b3_t. These must be     */
  /* excluded from the s2t map; make a mask of these locations.      */
  /* Use the dummy array in ws->wincon->s2 for the mask.             */
  mask = ws->wincon->s2;
  memset(mask, 0, ws->sgsizS * sizeof(int));
  for (i = 0; i < ws->nobc; i++) {
    open_bdrys_t *open = ws->open[i];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no2_e1; cc++) {
        c = open->obc_e1[cc];
	mask[c] = 1;
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no2_e2; cc++) {
        c = open->obc_e2[cc];
	mask[c] = 1;
      }
    }
  }

  /* Find the origin */
  i = j = 0;
  while (i < nce1 && j < nce2) {
    if (!isnan(dumpdata->cellx[j][i]) && !isnan(dumpdata->celly[j][i])) {
      xo = dumpdata->cellx[j][i];
      yo = dumpdata->celly[j][i];
      break;
    }
    i += 1;
    j += 1;
  }
  i = nce1 - 1; j = nce2 - 1;
  while (i >= 0 && j >= 0) {
    if (!isnan(dumpdata->cellx[j][i]) && !isnan(dumpdata->celly[j][i])) {
      xe = dumpdata->cellx[j][i];
      ye = dumpdata->celly[j][i];
      break;
    }
    i -= 1;
    j -= 1;
  }
  org = 0;
  if (xo < xe)
    org |= L_EDGE;
  else
    org |= R_EDGE;
  if (yo < ye)
    org |= B_EDGE;
  else
    org |= F_EDGE;
  if (org == (L_EDGE|B_EDGE))
    hd_warn("origin located in SW corner of target grid\n");
  else if (org == (L_EDGE|F_EDGE))
    hd_warn("origin located in NW corner of target grid\n");
  else if (org == (R_EDGE|B_EDGE))
    hd_warn("origin located in SE corner of target grid\n");
  else if (org == (R_EDGE|F_EDGE))
    hd_warn("origin located in NE corner of target grid\n");

  /*-----------------------------------------------------------------*/
  /* Set up the target grid mesh (with origin in the top-right cell  */
  /* corner) for the xyij_tree.                                      */
  /* Get the grid spacings at cell faces                             */
  if(!(master->tmode & INEXACT)) {
    tx = d_alloc_2d(nfe1+1, nfe2+1);
    ty = d_alloc_2d(nfe1+1, nfe2+1);
    set_trans_grid(nce1, nce2, tx, ty, dumpdata->cellx, dumpdata->celly, 
		   dumpdata->gridx, dumpdata->gridy, dumpdata->u1x, dumpdata->u1y, 
		   dumpdata->u2x, dumpdata->u2y, geom);
    t_xyij_tree = grid_xytoij_init(tx, ty, nfe1, nfe2);
    k = wt->nz-1;
    
    /*---------------------------------------------------------------*/
    /* Determine if the source grid is equal to, a subset or superset*/
    /* of the target grid.                                           */
    tgrid = 0;
    for (cc = 1; cc <= ws->b2_t; cc++) {
      c = cs = ws->w2_t[cc];
      i = ws->s2i[c];
      j = ws->s2j[c];
      /* Check the grid origin (SE corner) */
      if (mask[c]) continue;
      co = xytoc(wt, master->xyij_tree, gridx[j][i], gridy[j][i], &xo, &yo);
      io = (int)xo; jo = (int)yo;
      /* Check the eastern grid corner */
      ce = xytoc(wt, master->xyij_tree, gridx[j][i+1], gridy[j][i+1], &xe, &ye);
      ie = (int)xe; je = (int)ye;
      xe -= floor(xe); ye -= floor(ye);
      if (xe <= 2*EPS && ie-1 == io) ie--;
      ce = wt->map[wt->nz-1][je][ie];
      /* Check the northern grid corner */
      cn = xytoc(wt, master->xyij_tree, gridx[j+1][i], gridy[j+1][i], &xn, &yn);
      in = (int)xn; jn = (int)yn;
      xn -= floor(xn); yn -= floor(yn);
      if (yn <= 2*EPS && jn-1 == jo) jn--;
      cn = wt->map[wt->nz-1][jn][in];
      /* Check the north-eastern grid corner */
      cne = xytoc(wt, master->xyij_tree, gridx[j+1][i+1], gridy[j+1][i+1], &xne, &yne);
      ine = (int)xne; jne = (int)yne;
      xne -= floor(xne); yne -= floor(yne);
      if (xne <= 2*EPS && ine-1 == io) ine--;
      if (yne <= 2*EPS && jne-1 == jo) jne--;
      cne = wt->map[wt->nz-1][jne][ine];
      if (co == ce && co == cn && co == cne) {
	if (xe <= 2*EPS && xn <= 2*EPS && xne <= 2*EPS &&
	    ye <= 2*EPS && yn <= 2*EPS && yne <= 2*EPS) {
	  if (!(tgrid & (INEXACT|SUBSET))) 
	    tgrid = EXACT;
	  else if (!(tgrid & INEXACT)) 
	    tgrid = SUBSET;
	}
      } else {
	tgrid = INEXACT;
      }
    }
  } else
    tgrid = INEXACT;

  if (tgrid == EXACT) {
    master->togn = TOPRIGHT;
    hd_warn("Multi grid transport : exact source and target grids\n");
  } else if(tgrid == SUBSET) {
    master->togn = TOPRIGHT;
    hd_warn("Multi grid transport : source grid is a subset of targer grid\n");
  } else {
    hd_warn("Multi grid transport : inexact, non-subset, source and target grids\n");
    master->togn = BOTMLEFT;
    if (!(master->tmode & INEXACT)) {
      d_free_2d(tx);
      d_free_2d(ty);
      free((xytoij_tree_t *)t_xyij_tree);
    }
    /* Get the s2t map                                               */
    memset(s2t, 0, ws->sgsiz * sizeof(int));
    for (cc = 1; cc <= ws->b2_t; cc++) {
      c = cs = ws->w2_t[cc];
      ct = xytoc(wt, master->xyij_tree, ws->cellx[c], ws->celly[c], &x, &y);
      if (ct == 0) {
	ct = find_nearest(wt, ws->cellx[cs], ws->celly[cs]);
	if (ct == 0) {     /* Not in the target grid                 */
	  hd_warn("Can't map source@%d(%d %d) to target\n",
		  cs,ws->s2i[c], ws->s2j[c]);
	  continue;
	}
	else {
	  hd_warn("get_ts_map: source (%d %d) maps to source land: using target (%d %d)\n",
		  ws->s2i[cs], ws->s2j[cs], wt->s2i[ct], wt->s2j[ct]);
	  x = (double)wt->s2i[ct]+0.5;
	  y = (double)wt->s2j[ct]+0.5;
	}
      }
      /* Transform to the transport grid                             */
      x -= floor(x);
      y -= floor(y);
      r2tij(wt, ct, x, y, &xo, &yo);
      s2t[c] = ct;
      /* Get the t2s map for all surface and sub-surface cells       */
      clast = ct;
      c2 = wt->m2d[ct];
      top = 0.0;
      while (c != ws->zm1[c]) {
	bot = max(ws->gridz[c], ws->botz[cs]);
	z = 0.5 * (top + bot);
	w1 = wt->botz[c2];
	ct = ztoc(wt, c2, z);
	/*if (z < wt->cellz[ct]) ct = ws->zm1[ct];*/
	if (z < w1) {
	  hd_warn("get_st_map: source cellz %f %d(%d %d %d) is below target bottom %f at %d(%d %d %d): mapped to target %d\n",
		  z,c,ws->s2i[c],ws->s2j[c],ws->s2k[c],
		  w1,ct,wt->s2i[ct],wt->s2j[ct],wt->s2k[ct],clast);
	  ct = clast;
	}
	s2t[c] = ct;
	c = ws->zm1[c];
	top = bot;
	clast = ct;
      }
    }
    return(tgrid);
    wt->gs = NULL;
    wt->gs = (GRID_SPECS **)malloc(sizeof(GRID_SPECS *) * (wt->nz + 1));
    for (i = 0; i < wt->nz; i++)
      hd_grid_interp_init(wt->gs[i], master->temp, "GRID_LINEAR");

    return(tgrid);
  }

  /*-----------------------------------------------------------------*/
  /* Get the map from source to target grids                         */
  for (cc = 1; cc <= ws->b2_t; cc++) {
    /* Get the wet points in the source grid                         */
    c = ws->w2_t[cc];
    ii = ws->s2i[c];
    jj = ws->s2j[c];
    if (mask[c]) continue;
    /* Get the corresponding (i,j) in the source mesh                */
    i = ii + 1;
    j = jj + 1;

    /* Check the south-west quadrant                                 */
    check_quad(wt, t_xyij_tree, sx, sy, i, j, &x, &y, SW);
    cs = c;
    ct = wt->map[wt->nz-1][(int)y][(int)x];
    xoset[cs] = (x - (int)x);
    yoset[cs] = (y - (int)y);
    /* Get the s2t map for all surface and sub-surface cells       */
    while (cs != ws->zm1[cs]) {
      ct = ztoc(wt, wt->m2d[ct], ws->cellz[cs]);
      if (ws->cellz[cs] < wt->cellz[ct])
	ct = wt->zm1[ct];
      s2t[cs] = ct;
      cs = wt->zm1[cs];
    }
    /* s2t map for the sediment cell, interpolating the cell centre  */
    /* into the sediment.                                            */
    w1 = 2.0 * ws->gridz[ws->zp1[cs]] - ws->cellz[ws->zp1[cs]];
    ct = ztoc(wt, wt->m2d[ct], w1);
    if (ws->cellz[cs] < wt->cellz[ct])
      ct = wt->zm1[ct];
    s2t[cs] = ct;

    /* Check the north-east quadrant                                 */
    check_quad(wt, t_xyij_tree, sx, sy, i, j, &x, &y, NE);
    cs = c;
    ct = wt->map[wt->nz-1][(int)y][(int)x];
    xoset[ws->xp1[ws->yp1[cs]]] = (x - (int)x);
    yoset[ws->xp1[ws->yp1[cs]]] = (y - (int)y);
    /* Get the s2t map for all surface and sub-surface cells       */
    while (cs != ws->zm1[cs]) {
      ct = ztoc(wt, wt->m2d[ct], ws->cellz[cs]);
      if (ws->cellz[cs] < wt->cellz[ct])
	ct = wt->zm1[ct];
      s2t[ws->xp1[ws->yp1[cs]]] = wt->xp1[wt->yp1[ct]];
      cs = wt->zm1[cs];
    }
    /* s2t map for the sediment cell                                 */
    w1 = 2.0 * ws->gridz[ws->zp1[cs]] - ws->cellz[ws->zp1[cs]];
    ct = ztoc(wt, wt->m2d[ct], w1);
    if (ws->cellz[cs] < wt->cellz[ct])
      ct = wt->zm1[ct];
    s2t[ws->xp1[ws->yp1[cs]]] = wt->xp1[wt->yp1[ct]];

    /* Check the north-west quadrant                                 */
    check_quad(wt, t_xyij_tree, sx, sy, i, j, &x, &y, NW);
    cs = c;
    ct = wt->map[wt->nz-1][(int)y][(int)x];
    xoset[ws->yp1[cs]] = (x - (int)x);
    yoset[ws->yp1[cs]] = (y - (int)y);
    while (cs != ws->zm1[cs]) {
      ct = ztoc(wt, wt->m2d[ct], ws->cellz[cs]);
      if (ws->cellz[cs] < wt->cellz[ct])
	ct = wt->zm1[ct];
      s2t[ws->yp1[cs]] = wt->yp1[ct];
      cs = wt->zm1[cs];
    }
    /* s2t map for the sediment cell                                 */
    w1 = 2.0 * ws->gridz[ws->zp1[cs]] - ws->cellz[ws->zp1[cs]];
    ct = ztoc(wt, wt->m2d[ct], w1);
    if (ws->cellz[cs] < wt->cellz[ct])
      ct = wt->zm1[ct];
    s2t[ws->yp1[cs]] = wt->yp1[ct];

    /* Check the south-east quadrant                                 */
    check_quad(wt, t_xyij_tree, sx, sy, i, j, &x, &y, SE);
    cs = c;
    ct = wt->map[wt->nz-1][(int)y][(int)x];
    xoset[ws->xp1[cs]] = (x - (int)x);
    yoset[ws->xp1[cs]] = (y - (int)y);
    while (cs != ws->zm1[cs]) {
      ct = ztoc(wt, wt->m2d[ct], ws->cellz[cs]);
      if (ws->cellz[cs] < wt->cellz[ct])
	ct = wt->zm1[ct];
      s2t[ws->xp1[cs]] = wt->xp1[ct];
      cs = wt->zm1[cs];
    }
    /* s2t map for the sediment cell                                 */
    w1 = 2.0 * ws->gridz[ws->zp1[cs]] - ws->cellz[ws->zp1[cs]];
    ct = ztoc(wt, wt->m2d[ct], w1);
    if (ws->cellz[cs] < wt->cellz[ct])
      ct = wt->zm1[ct];
    s2t[ws->xp1[cs]] = wt->xp1[ct];
  }

  /* Capture multiple ghost cells approached from alternative directions */
  for (cc = 1; cc <= ws->nbpt; cc++) {
    int ci, cg;
    cg = ws->bpt[cc];
    ci = ws->bin[cc];
    ct = s2t[ci];

    if (ct != 0 && s2t[cg] == 0) {
      int *igm, *gim;
      int c1, c2, c3, c4;
      /* Get the direction the ghost is approached from                  */
      if (cg == ws->xm1[ci]) { igm = wt->xm1; gim = wt->xp1; }
      if (cg == ws->xp1[ci]) { igm = wt->xp1; gim = wt->xm1; }
      if (cg == ws->ym1[ci]) { igm = wt->ym1; gim = wt->yp1; }
      if (cg == ws->yp1[ci]) { igm = wt->yp1; gim = wt->ym1; }

      /* Check that this cell is a multiple ghost cell where one of the  */
      /* multiple ghosts is associated with a valid s2t map.             */
      c1 = ws->yp1[ws->ym1[cg]];
      c2 = ws->ym1[ws->yp1[cg]];
      c3 = ws->xp1[ws->xm1[cg]];
      c4 = ws->xm1[ws->xp1[cg]];
      if (c1 != cg && s2t[c1] != 0) s2t[cg] = igm[gim[s2t[c1]]];
      else if (c2 != cg && s2t[c2] != 0) s2t[cg] = igm[gim[s2t[c2]]];
      else if (c3 != cg && s2t[c3] != 0) s2t[cg] = igm[gim[s2t[c3]]];
      else if (c4 != cg && s2t[c4] != 0) s2t[cg] = igm[gim[s2t[c4]]];
      /*
      c = ws->yp1[ws->ym1[cg]];
      if (c != cg && s2t[c] != 0) s2t[cg] = igm[gim[s2t[c]]];
      c = ws->ym1[ws->yp1[cg]];
      if (c != cg && s2t[c] != 0) s2t[cg] = igm[gim[s2t[c]]];
      c = ws->xp1[ws->xm1[cg]];
      if (c != cg && s2t[c] != 0) s2t[cg] = igm[gim[s2t[c]]];
      c = ws->xm1[ws->xp1[cg]];
      if (c != cg && s2t[c] != 0) s2t[cg] = igm[gim[s2t[c]]];
      */
    }
  }

  /*-----------------------------------------------------------------*/
  /* Check every source cell maps to a target cell                   */
  for (cc = 1; cc <= ws->b3_t; cc++) {
    cs = ws->w3_t[cc];            /* Source cell location            */
    ct = s2t[cs];                 /* Target cell location            */
    if ((ct <= 0 || ct >= wt->sgnum) && DEBUG("init_m"))
      dlog("init_m", "Undefined source-target map : s %d (%d %d %d) -> t %d (%d %d %d)\n",
	     cs,ws->s2i[cs],ws->s2j[cs],ws->s2k[cs],ct,wt->s2i[ct],wt->s2j[ct],wt->s2k[ct]);
    /* NW quadrant                                                    */
    ct = s2t[ws->yp1[cs]];
    if ((ct <= 0 || ct >= wt->sgnum) && DEBUG("init_m"))
      dlog("init_m", "Undefined source-target map at NW quadrant : s %d (%d %d %d) -> t %d (%d %d %d)\n",
	     ws->yp1[cs],ws->s2i[cs],ws->s2j[cs]+1,ws->s2k[cs],ct,wt->s2i[ct],wt->s2j[ct],wt->s2k[ct]);
    /* SE quadrant                                                    */
    ct = s2t[ws->xp1[cs]];
    if ((ct <= 0 || ct >= wt->sgnum) && DEBUG("init_m"))
      dlog("init_m", "Undefined source-target map at SE quadrant : s %d (%d %d %d) -> t %d (%d %d %d)\n",
	     ws->xp1[cs],ws->s2i[cs]+1,ws->s2j[cs],ws->s2k[cs],ct,wt->s2i[ct],wt->s2j[ct],wt->s2k[ct]);
    /* NE quadrant                                                    */
    ct = s2t[ws->xp1[ws->yp1[cs]]];
    if ((ct <= 0 || ct >= wt->sgnum) && DEBUG("init_m"))
      dlog("init_m", "Undefined source-target map at NE quadrant : s %d (%d %d %d) -> t %d (%d %d %d)\n",
	     ws->xp1[ws->yp1[cs]],ws->s2i[cs]+1,ws->s2j[cs]+1,ws->s2k[cs],ct,wt->s2i[ct],wt->s2j[ct],wt->s2k[ct]);
  }

  d_free_2d(tx);
  d_free_2d(ty);
  free((xytoij_tree_t *)t_xyij_tree);
  return(tgrid);
}

/* END get_st_map()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check that the corners of a source mesh for a given    */
/* quadrant with origin in the bottom-left corner lies within a      */
/* target mesh.                                                      */
/* (x,y) contain the locations of the original wet cell from which   */
/* the quadrants were derived, relative to the original (non-mesh)   */
/* grid; i.e.                                                        */
/* quad = SW => (x,y) = NE corner                                    */
/* quad = NE => (x,y) = SW corner                                    */
/* quad = NW => (x,y) = SE corner                                    */
/* quad = SE => (x,y) = NW corner                                    */
/*                                                                   */
/*        |----|----|                                                */
/*        | NW | NE |                                                */
/*        |--(i,j)--|                                                */
/*        | SW | SE |                                                */
/*        |----|----|                                                */
/*                                                                   */
/*-------------------------------------------------------------------*/
int check_quad(geometry_t *wt,             /* Target geometry        */
	       xytoij_tree_t *t_xyij_tree, /* Target mesh tree       */
	       double **sx,                /* Source x mesh coords   */
	       double **sy,                /* Source y mesh coords   */
	       int i, int j,               /* Quadrant origin        */
	       double *x, double *y,       /* Fractions for ne       */
	       int quad)
{
  double xo, xe, xn, xne;
  double yo, ye, yn, yne;
  int ret, c, ir, jr;
  int io, ie, in, ine;
  int jo, je, jn, jne;

  if (quad == SW) {
    i -= 1; j -= 1;
  }
  if (quad == NW) {
    i -= 1;
  }
  if (quad == SE) {
    j -= 1;
  }

  /* Set the quadrant origin                                         */
  c = xytoc(wt, t_xyij_tree, sx[j][i], sy[j][i], &xo, &yo);
  io = (int)xo; jo = (int)yo;
  xo -= (double)io; yo -= (double)jo;
  if (quad == NE) {
    ir = io; jr = jo;
  }
  /* Quadrant north-east location                                    */
  c = xytoc(wt, t_xyij_tree, sx[j+1][i+1], sy[j+1][i+1], &xne, &yne);
  ine = (int)xne; jne = (int)yne;
  xne -= (double)ine; yne -= (double)jne;
  *x = (xne <= 2*EPS) ? 0.0 : 1.0 - xne;
  *y = (yne <= 2*EPS) ? 0.0 : 1.0 - yne;
  if (quad == SW) {
    ir = ine; jr = jne;
  }
  if (ine - 1 == io && xne <= 2*EPS) ine--;
  if (jne - 1 == jo && yne <= 2*EPS) jne--;
  /* Quadrant east location                                          */
  c = xytoc(wt, t_xyij_tree, sx[j][i+1], sy[j][i+1], &xe, &ye);
  ie = (int)xe; je = (int)ye;
  xe -= (double)ie; ye -= (double)je;
  if (quad == NW) {
    ir = ie; jr = je;
  }
  if (ie - 1 == io && xe <= 2*EPS)ie--;
  /* Quadrant north location                                     */
  c = xytoc(wt, t_xyij_tree, sx[j+1][i], sy[j+1][i], &xn, &yn);
  in = (int)xn; jn = (int)yn;
  xn -= (double)in; yn -= (double)jn;
  if (quad == SE) {
    ir = in; jr = jn;
  }
  if (jn - 1 == jo && yn <= 2*EPS)jn--;
  /* Add the grid coordinate (grid coord = mesh coord - 1)       */
  *x += (double)(ir - 1);
  *y += (double)(jr - 1);
  /* Check if all corners of the source mesh lie in the same     */
  /* target mesh.                                                */
  if (ie == io && in == io && ine == io && je == jo && jn == jo && jne == jo) {
    if (xe <= 2*EPS && xn <= 2*EPS && xne <= 2*EPS &&
	ye <= 2*EPS && yn <= 2*EPS && yne <= 2*EPS) {
      return(EXACT);
    } else {
      return(SUBSET);
    }
  } else {
    return(INEXACT);
  }
}

/* END check_quad()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets dimensional and (i,j) vectors at ghost cells                 */
/*-------------------------------------------------------------------*/
void s2ij_setghosts(geometry_t *geom)
{
  int cc, c, ci;
  int i, j, ig, jg;

  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];
    ci = geom->bin[cc];
    i = ig = geom->s2i[ci];
    j = jg = geom->s2j[ci];

    if (geom->xp1[ci] == c) ig = i + 1;
    else if (geom->xm1[ci] == c) ig = i - 1;
    else if (geom->yp1[ci] == c) jg = j + 1;
    else if (geom->ym1[ci] == c) jg = j - 1;
    else if (geom->yp1[geom->xp1[ci]] == c || geom->xp1[geom->yp1[ci]] == c) {
      ig = i + 1;
      jg = j + 1;
    } else if (geom->yp1[geom->xm1[ci]] == c || geom->xm1[geom->yp1[ci]] == c) {
      ig = i - 1;
      jg = j + 1;
    } else if (geom->ym1[geom->xp1[ci]] == c || geom->xp1[geom->ym1[ci]] == c) {
      ig = i + 1;
      jg = j - 1;
    } else if (geom->ym1[geom->xm1[ci]] == c || geom->xm1[geom->ym1[ci]] == c) {
      ig = i - 1;
      jg = j - 1;
    }

    ig = min(max(ig, 0), geom->nce1-1);
    jg = min(max(jg, 0), geom->nce2-1);

    geom->s2i[c] = ig;
    geom->s2j[c] = jg;
  }
}

/* END s2ij_setghosts()                                              */
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
  tpg->xp1 = geom->xp1;
  tpg->xm1 = geom->xm1;
  tpg->yp1 = geom->yp1;
  tpg->ym1 = geom->ym1;
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

/*-------------------------------------------------------------------*/
/* Returns the nearest wet cell to geographic location (x,y)         */
/*-------------------------------------------------------------------*/
int find_nearest(geometry_t *geom, double x, double y)
{
  int c, cc, c1, cs;
  double dx, dy, dist, mindist = HUGE;

  for (cc = 1; cc < geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    dx = x - geom->cellx[c];
    dy = y - geom->celly[c];
    dist = dx * dx + dy * dy;
    if (dist < mindist) {
      cs = c;
      mindist = dist;
    }
  }

  if (cs >0 && cs < geom->enonS)
    return(cs);
  else
    return(0);
}
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
