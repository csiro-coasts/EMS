/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/particles/pt.c
 *  
 *  Description:
 *  particle_t tracking code for meco
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: pt.c 5916 2018-09-05 03:51:46Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <netcdf.h>
#include "hd.h"

#define RADIUS 6370997.0
#define ECC 0.0

#define MAX_COL    127
#define HIST_SCALE 5
#define CONSTANT 2 /*settling velocity tpes*/
#define PSTOKES  4
#define DIURNAL  8

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

/* Prototypes */
void particles_to_conc(master_t *master, long np, particle_t *p);
void ptgrid_xytoij(master_t *master, long np, particle_t *p);
void ptgrid_ijtoxy(master_t *master, long np, particle_t *p);
void pt_fit_to_grid(master_t *master, long np, particle_t *p);
void pt_write_at_t(master_t *master, double t);
void pt_move(master_t *master, dump_data_t *dumpdata, particle_t *p,
             double dt, double maxdh);
void pt_set_age(master_t *master, particle_t *p, double dt);
void pt_set_size(master_t *master, particle_t *p, double dt, int n);
void pt_set_svel(master_t *master, particle_t *p, double dt, int n);
void pt_new(master_t *master, long np, particle_t *p);
int hd_pt_create(master_t *master, char *name);
int get_pos_m(master_t *master, particle_t *p, int c, int ci, double u, double v, double w,
	      double *cx, double *cy, double *cz, double dt);

/*-------------------------------------------------------------------*/
/* Event called by the time scheduler. This routine writes out the   */
/* particle positions, and possibly resets the particle positions,   */
/* as required.                                                      */
/*-------------------------------------------------------------------*/
double ptrack_event(sched_event_t *event, double t)
{
  double tout = schedule->stop_time;
  master_t *master = (master_t *)schedGetPublicData(event);
    int n;
  /* If no particle output file is open, then assume nothing needs */
  /* to be done.  */
  if (master->ptfid < 0)
    return schedule->stop_time;

  /* If before the particle tracking start time, then nothing to do */
  if (t + master->dt / 10.0 < master->ptstart)
    return master->ptstart;

  /* Write next output if required */
  if (t + master->dt / 10.0 >= master->ptout_t)
    pt_write_at_t(master, t);

  /* Reset particles if required */
  if (t >= master->ptreset_t) {

    if (DEBUG("particles"))
      dlog("particles", "Reseting particles at t=%f\n", master->t);

    pt_read(master->ptinname, master->ptinrec, &master->ptn,
            &master->ptp, NULL, NULL, NULL);

    if (master->pt_nsource <= 0) {
      ptgrid_xytoij(master, master->ptn, master->ptp);
      pt_fit_to_grid(master, master->ptn, master->ptp);
    }
    master->ptreset_t += master->ptreset;
  }

  /* Determine when next output of particles should be */
  if (master->ptout_t < tout)
    tout = master->ptout_t;

  /* Determine when next reset of particles should be */
  if (master->ptreset_t < tout)
    tout = master->ptreset_t;

  return tout;
}

/* END ptrack_event()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* scheduler_t initialisation routine for particles                  */
/*-------------------------------------------------------------------*/
int ptrack_init(sched_event_t *event)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  char buf[MAXSTRLEN];

  /* Set output and reset times */
  master->ptout_t = master->ptstart;
  master->ptreset_t = master->ptstart + master->ptreset;
  master->ptnext_t = master->ptstart + master->ptstep / 2;

  /* Create the output file */
  strcpy(buf, master->ptoutname);
  if(!endswith(master->ptoutname,".nc"))
    strcat(buf,".nc");

  /*
  master->ptfid = pt_create(buf, master->ptn, master->output_tunit,
			    master->pt_dumpf);
  */
  if (!(master->runmode & DUMP))
    master->ptfid = hd_pt_create(master, buf);

  return ptrack_event(event, schedule->t);
}

/* END ptrack_init()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* scheduler_t termination routine for particles                     */
/*-------------------------------------------------------------------*/
void ptrack_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);

  if (master->ptfid >= 0) {
    if (t + master->dt / 10 >= master->ptout_t)
      pt_write_at_t(master, t);
    if (master->ptmsk) s_free_1d(master->ptmsk);
    nc_close(master->ptfid);
  }
  master->ptfid = -1;
}

/* END ptrack_cleanup()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* End particle tracking                                             */
/*-------------------------------------------------------------------*/
void ptrack_end(void)
{
  sched_deregister(schedule, "particles");
}

/* END ptrack_end()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* particle_t tracking initialisation routine                        */
/*-------------------------------------------------------------------*/
void pt_params_init(master_t *master, FILE * fp)
{
  char buf[MAXSTRLEN];
  char key[MAXSTRLEN];
  int i;
  geometry_t *geom = master->geom;
  int n;
  int restart = 0;

  master->mage = master->magec = 0.0;
  master->shist = 0;
  master->phist = NULL;
  master->pt_svel = NULL;

  /* Read particle tracking file parameters                          */
  prm_set_errfn(hd_silent_warn);
  if (master->do_pt == 0) {
    /* No particle tracking input file - so nothing further to do    */
    return;
  }

  prm_set_errfn(hd_quit);
  prm_read_int(fp, "PT_InputRecord", &master->ptinrec);
  if (prm_read_char(fp, "PT_Restart", buf))
    restart = is_true(buf);
  prm_read_char(fp, "PT_OutputFile", master->ptoutname);
  prm_get_time_in_secs(fp, "PT_StartTime", &master->ptstart);
  prm_get_time_in_secs(fp, "PT_StopTime", &master->ptend);
  prm_get_time_in_secs(fp, "PT_OutputStep", &master->ptoutinc);
  prm_get_time_in_secs(fp, "PT_TimeStep", &master->ptstep);
  prm_get_time_in_secs(fp, "PT_ResetStep", &master->ptreset);

  if (master->ptstart < master->t)
    master->ptstart = master->t;

  /* Read horizontal diffusion coefficient                           */
  prm_read_double(fp, "PT_KH", &master->pt_kh);

  /* Read vertical diffusion coefficient multiplier                  */
  prm_read_double(fp, "PT_KZ_MULT", &master->pt_kz_mult);

  /* Read mass of each new particle in kg                            */
  prm_read_double(fp, "PT_Mass", &master->pt_mass);

  /* Read flag which determines behaviour on reaching an open        */
  /* boundary.  */
  prm_read_char(fp, "PT_StickyBoundary", buf);
  master->pt_stickybdry = is_true(buf);

  prm_set_errfn(hd_silent_warn);

  master->pt_dumpf = 0;
  /* Read the age colour stretch                                     */
  if (prm_get_time_in_secs(fp, "PT_AgeLimit", &master->pt_agelim)) {
    master->pt_dumpf |= PT_AGE;
    master->shist = HIST_SCALE * master->pt_agelim / 86400;
  }

  /* Read the size colour stretch                                    */
  if (prm_read_double(fp, "PT_SizeLimit", &master->pt_sizelim)) {
    master->pt_dumpf |= PT_SIZE;
  }

  /* Set up a histogram array                                        */
  sprintf(key, "PT_Histogram_DT");
  if (prm_get_time_in_secs(fp, key, &tsphist.tsdt) && master->shist) {
    master->phist = i_alloc_1d(master->shist);
    memset(master->phist, 0,  master->shist * sizeof(int));
    if (tsphist.fp == NULL) {
      strcpy(tsphist.pname, "part_hist.ts");
      tsphist.fp = fopen("part_hist.ts", "w");
      tsphist.tsdt = 3600.0;
      tsphist.master = master;
      tsphist.i = 0;
      tsphist.j = 0;
      tsphist.k = 0;
      tsphist.tsout = master->t;
      fprintf(tsphist.fp, "## COLUMNS %d\n", master->shist + 1);
      fprintf(tsphist.fp, "##\n");
      fprintf(tsphist.fp, "## COLUMN1.name  Time\n");
      fprintf(tsphist.fp, "## COLUMN1.long_name  Time\n");
      fprintf(tsphist.fp,
              "## COLUMN1.units  days since 1990-01-01 00:00:00 +10\n");
      for (i = 0; i < master->shist; i++)
	fprintf(tsphist.fp, "## COLUMN%1.1d.name  day%1.1dto%1.1d\n",
		i + 2, i, i+1);
      fprintf(tsphist.fp, "##\n");
    }
  }

  /* Read the mean age region */
  sprintf(key, "PT_age_region");
  if (prm_read_char(fp, key, buf)) {
    char *files[MAXSTRLEN * MAXNUMARGS];
    double d1, d2, top, bot;
    int cc, c, m, nf;

    master->ptmsk = s_alloc_1d(geom->sgsiz);
    memset(master->ptmsk, 0, geom->sgsiz * sizeof(short));
    nf = parseline(buf, files, MAXNUMARGS);
    /* Initialize the flushing tracer */
    /* Get the range */
    if (prm_read_char(fp, "PT_age_range", buf)) {
      if (sscanf(buf, "%lf %lf", &d1, &d2) == 2) {
	top = d1; bot = d2;
	if (bot > top) {
	  d1 = bot;
	  bot = top;
	  top = d1;
	}
      }
    } else {
      top = HUGE; 
      bot = -HUGE;
    }
    if (nf == 1) {
      int *iloc, *jloc;
      read_blocks(fp, key, &n, &iloc, &jloc);
      for (m = 1; m <= n; m++) {
	c = geom->cc2s[iloc[m]];
	while (c != geom->zm1[c] && geom->cellz[c] <= top && geom->cellz[c] >= bot) {
	  master->ptmsk[c] = 1;
	  c = geom->zm1[c];
	}
      }
      i_free_1d(iloc);
    } else {
      int nr;
      double *regionid;
      regionid = d_alloc_1d(geom->sgsiz);
      nr = read_regioni(master, files[0], regionid);
      if (nr) {
	int rgn;
	for (m = 1; m < nf; m++) {
	  sscanf(files[m], "%d", &rgn);
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c]) {
	      while (c != geom->zm1[c] && geom->cellz[c] <= top && geom->cellz[c] >= bot) {
		master->ptmsk[c] = 1;
		c = geom->zm1[c];
	      }
	    }
	  }
	}
      } else
	hd_warn("Particle age: No points specified for age region.\n");
    }
  }

  /* Read the vertical swimming velocity file                        */
  master->wvel_i = NULL;
  sprintf(key, "PT_w_file");
  if(prm_read_char(fp, key, buf)){
    master->wvel_i = (pt_ts_t *)malloc(sizeof(pt_ts_t));
    pt_ts_t *wdata = master->wvel_i;
    wdata->val = NULL;
    wdata->ts = hd_ts_read(master, buf, 0);
    wdata->id = ts_get_index(wdata->ts, fv_get_varname(buf, "wpt", buf));
    if (wdata->id < 0)
      hd_quit("pt_params_init: The file '%s' does not contain the settling variable 'wpt' .\n", buf);
    wdata->val = d_alloc_1d(geom->sgsiz);
    memset(wdata->val, 0, geom->sgsiz * sizeof(double)); 
  }

  /* Read the horizontal swimming velocity file                      */
  master->uvel_i = NULL;
  master->vvel_i = NULL;
  sprintf(key, "PT_uv_file");
  if(prm_read_char(fp, key, buf)){
    master->uvel_i = (pt_ts_t *)malloc(sizeof(pt_ts_t));
    master->vvel_i = (pt_ts_t *)malloc(sizeof(pt_ts_t));
    pt_ts_t *udata = master->uvel_i;
    pt_ts_t *vdata = master->vvel_i;
    udata->val = NULL;
    vdata->val = NULL;
    udata->ts = hd_ts_read(master, buf, 0);
    udata->id = ts_get_index(udata->ts, fv_get_varname(buf, "upt", buf));
    if (udata->id < 0)
      hd_quit("pt_params_init: The file '%s' does not contain the settling variable 'upt' .\n", buf);
    udata->val = d_alloc_1d(geom->sgsiz);
    memset(udata->val, 0, geom->sgsiz * sizeof(double)); 
    vdata->ts = hd_ts_read(master, buf, 0);
    vdata->id = ts_get_index(vdata->ts, fv_get_varname(buf, "vpt", buf));
    if (vdata->id < 0)
      hd_quit("pt_params_init: The file '%s' does not contain the settling variable 'vpt' .\n", buf);
    vdata->val = d_alloc_1d(geom->sgsiz);
    memset(vdata->val, 0, geom->sgsiz * sizeof(double)); 
  }

  /* Reads the mortality file                                        */
  master->mort_i = NULL;
  sprintf(key,"PT_mortality_file");
  if(prm_read_char(fp, key, buf)){
    master->mort_i = (pt_ts_t *)malloc(sizeof(pt_ts_t));
    pt_ts_t *wdata = master->mort_i;
    wdata->ts = hd_ts_read(master, buf, 0);
    wdata->id = ts_get_index(wdata->ts, fv_get_varname(buf, "mpt", buf));
    if(wdata->id < 0)
      hd_quit("pt_params_init: The file '%s' does not contain the mortality variable 'mpt' .\n", buf);
  }
  
  /* If the release rate is specified, then a start and end          */
  /* coordinates for the continous line must be specified.           */
  master->pt_nsource = 0;
  prm_set_errfn(hd_silent_warn);
  prm_read_int(fp, "PT_NSOURCE", &master->pt_nsource);
  prm_set_errfn(hd_quit);

  if (master->pt_nsource > 0) {
    master->pt_rate = d_alloc_1d(master->pt_nsource);
    master->pt_colour =
      (short *)malloc(sizeof(short) * master->pt_nsource);
    master->pt_x1 = d_alloc_1d(master->pt_nsource);
    master->pt_y1 = d_alloc_1d(master->pt_nsource);
    master->pt_z1 = d_alloc_1d(master->pt_nsource);
    master->pt_x2 = d_alloc_1d(master->pt_nsource);
    master->pt_y2 = d_alloc_1d(master->pt_nsource);
    master->pt_z2 = d_alloc_1d(master->pt_nsource);
    master->pt_accum = d_alloc_1d(master->pt_nsource);
    master->pt_size = d_alloc_1d(master->pt_nsource);
    master->pt_decay = d_alloc_1d(master->pt_nsource);
    master->pt_svel = d_alloc_1d(master->pt_nsource);
    master->pt_dens = d_alloc_1d(master->pt_nsource);
    master->svel_type = i_alloc_1d(master->pt_nsource);
    master->pt_sper = d_alloc_1d(master->pt_nsource);
    memset(master->pt_svel, 0, master->pt_nsource*sizeof(double)); 

    for (i = 0; i < master->pt_nsource; ++i) {
      int colourBit = 0;
      master->pt_accum[i] = master->pt_size[i] = 0.0;

      /* Read rate                                                   */
      sprintf(key, "PT_Source%d.Rate", i);
      prm_read_double(fp, key, &master->pt_rate[i]);

      /* Read the colour bit value                                   */
      sprintf(key, "PT_Source%d.ColourBit", i);
      prm_read_int(fp, key, &colourBit);
      master->pt_colour[i] = 1 << colourBit;
      master->pt_colour[i] &= ~(PT_ACTIVE | PT_LOST);

      /* Read the initial size (diameter (m))                        */
      prm_set_errfn(hd_silent_warn);
      if (master->pt_sizelim) {
	master->pt_size[i] = master->pt_decay[i] = 0.0;
	sprintf(key, "PT_Source%d.Size", i);
	if (!prm_read_double(fp, key, &master->pt_size[i]))
	  hd_warn("pt_params_init: Particle size zero -> particle is always lost.\n");
	sprintf(key, "PT_Source%d.Growth", i);
	if (!prm_get_time_in_secs(fp, key, &master->pt_decay[i])) {
	  sprintf(key, "PT_Source%d.Decay", i);
	  if (prm_get_time_in_secs(fp, key, &master->pt_decay[i]))
	    master->pt_decay[i] *= -1.0;
	}
      }

      /* Read the settling velocity type (ms-1)                      */
      /* Set the default value                                       */
      master->svel_type[i] = NONE;
      /* Goes through the options                                    */
      sprintf(key, "PT_Source%d.Stype",i);
      prm_read_char(fp, key, buf);
      if (strcmp(buf, "CONSTANT") == 0) {
	sprintf(key, "PT_Source%d.Svel", i);        
	if (prm_read_double(fp, key, &master->pt_svel[i])) {
	  master->svel_type[i]= CONSTANT;
	} else {
	  hd_quit("pt_params_init: Constant settling velocity type requires settling velocity attribute (PT_Source.Svel)#. \n");
	}
      } else if (strcmp(buf, "STOKES") == 0) {
	sprintf(key, "PT_Source%d.Dens", i);
	if (prm_read_double(fp, key, &master->pt_dens[i])) {
	  if (master->pt_size[i] == 0.0) 
	    hd_quit("pt_params_init: Stokes settling velocity type requires Size attribute(PT_Source#.Size).\n"); 
	  master->svel_type[i] = PSTOKES;
	} else { hd_quit("pt_params_init: Stokes settling velocity type requires Density attribute (PT_Source#.Dens).\n");
	}
      } else if (strcmp(buf, "DIURNAL") == 0) {
        sprintf(key, "PT_Source%d.Svel",i);
	if (prm_read_double(fp, key, &master->pt_svel[i])) {
	  sprintf(key, "PT_Source%d.Sper",i);
	  if (prm_read_char(fp, key, buf)){
	    tm_scale_to_secs(buf, &master->pt_sper[i]);
	    master->svel_type[i]= DIURNAL;
          } else {
	    hd_quit("pt_params_init: Diurnal settling velocity type requires Periodicity attribute (PT_Source#.Sper).\n ");
	  }
	} else {
	  hd_quit("pt_params_init: Diurnal settling velocity type requires Settling Velocity attribute (PT_Source#.Svel).\n ");
	}
      } else if(strcmp(buf, "FILEIN") == 0) {
	char filename[MAXSTRLEN];
	sprintf(key, "PT_Source%d.File" ,i);
	if (prm_read_char(fp, key, filename)){
	  if (i == 0)
	    master->wvel_s = (pt_ts_t **)malloc(sizeof(pt_ts_t *) * master->pt_nsource);
	  master->wvel_s[i] = (pt_ts_t *)malloc(sizeof(pt_ts_t *));
	  pt_ts_t *wdata = master->wvel_s[i];
	  wdata->ts = hd_ts_read(master, filename, 0);
	  wdata->id = ts_get_index(wdata->ts, fv_get_varname(filename, "wpt", buf));
	  if (wdata->id < 0)
	    hd_quit("pt_params_init: The file '%s' does not contain the settling variable 'wpt' .\n", filename);
	  master->svel_type[i] = FILEIN;
	} else {
	  hd_quit("pt_params_init: File input requires settling velocity time series file.(PT_Source#.File)\n");
	}
      }

      /* Read start coordinate for particle source */
      sprintf(key, "PT_Source%d.StartLocation", i);
      prm_read_char(fp, key, buf);
      if (sscanf(buf, "%lf %lf %lf", &master->pt_x1[i], &master->pt_y1[i],
                 &master->pt_z1[i]) != 3)
        hd_quit("pt_params_init: Can't read StartLocation values\n");

      /* Read end coordinate for particle source */
      sprintf(key, "PT_Source%d.EndLocation", i);
      prm_read_char(fp, key, buf);
      if (sscanf(buf, "%lf %lf %lf", &master->pt_x2[i], &master->pt_y2[i],
                 &master->pt_z2[i]) != 3)
        hd_quit("pt_params_init: Can't read EndLocation values\n");

    }

  } else {                      /* Default to old behaviour */

    /* Read input rate of new particles per second */
    double rate = 0.0;
    prm_set_errfn(hd_silent_warn);
    prm_read_double(fp, "PT_Rate", &rate);
    prm_set_errfn(hd_quit);
    if (rate > 0.0) {
      master->pt_nsource = 1;
      master->pt_rate = d_alloc_1d(master->pt_nsource);
      master->pt_colour =
        (short *)malloc(sizeof(short) * master->pt_nsource);
      master->pt_x1 = d_alloc_1d(master->pt_nsource);
      master->pt_y1 = d_alloc_1d(master->pt_nsource);
      master->pt_z1 = d_alloc_1d(master->pt_nsource);
      master->pt_x2 = d_alloc_1d(master->pt_nsource);
      master->pt_y2 = d_alloc_1d(master->pt_nsource);
      master->pt_z2 = d_alloc_1d(master->pt_nsource);
      master->pt_accum = d_alloc_1d(master->pt_nsource);
      master->pt_rate[0] = rate;
      master->pt_colour[0] = 0;
      master->pt_accum[0] = 0.0;
      master->pt_size[0] = 0.0;

      /* Read start coordinate for particle source */
      prm_read_char(fp, "PT_StartLocation", buf);
      if (sscanf(buf, "%lf %lf %lf", &master->pt_x1[0], &master->pt_y1[0],
                 &master->pt_z1[0]) != 3)
        hd_quit("pt_params_init: Can't read PT_StartLocation values\n");

      /* Read end coordinate for particle source */
      prm_read_char(fp, "PT_EndLocation", buf);
      if (sscanf(buf, "%lf %lf %lf", &master->pt_x2[0], &master->pt_y2[0],
                 &master->pt_z2[0]) != 3)
        hd_quit("pt_params_init: Can't read PT_EndLocation values\n");
    }
  }

  /* Allocate the memory for the particle concentrations, and */
  /* initialise to NaN.  */
  for (i = 1; i < geom->sgsiz; i++) {
    master->ptconc[i] = NaN;
  }

  /* Read particles and convert to grid coords */
  pt_read(master->ptinname,
          master->ptinrec, &master->ptn, &master->ptp, NULL, NULL, NULL);

  /* Allocate memory for the particle to source map */
  master->pt_sm = s_alloc_1d(master->ptn);

  if (master->pt_nsource > 0) {
    pt_new(master, master->ptn, master->ptp);
    for (i = 0; i < master->ptn; i++) {
      /*
      master->ptp[i].e1 = master->pt_x1[0];
      master->ptp[i].e2 = master->pt_y1[0];
      */
      master->ptp[i].e1 = 0.0;
      master->ptp[i].e2 = 0.0;
      master->ptp[i].e3 = 0.0;
    }
    if (restart) {
      ptgrid_xytoij(master, master->ptn, master->ptp);
      pt_fit_to_grid(master, master->ptn, master->ptp);
    }
  } else {
    master->pt_sizelim = 0.0;
    ptgrid_xytoij(master, master->ptn, master->ptp);
    pt_fit_to_grid(master, master->ptn, master->ptp);
  }

  /* Initialise the particle concentrations */
  particles_to_conc(master, master->ptn, master->ptp);


  /* Allocate array for max vert diffusion values */
  master->maxdiffw = f_alloc_1d(geom->sgsizS);

  /* Register an interest with the tscheduler for dumping the output */
  /* output particle file.  */
  sched_register(schedule, "particles", ptrack_init,
                 ptrack_event, ptrack_cleanup, master, NULL, NULL);

  /* Throw away some random numbers. This is an attempt to avoid rather
     bizzare behaviour which can occur under certain circustances. If the
     utility used to generate initial particle positions uses the same set
     of random numbers as those used on the first timestep for
     diffusivities (see ptmove below), the particle positions and
     diffusivities are correlated, which can lead to unexpected and
     entertaining results for the first time step! It probably isn't
     necessary to throw away as many as shown below, but it probably
     doesn't hurt, either. */
  for (i = 0; i < master->ptn; i++) {
    int raninit = 0;
    ran3(&raninit);
    ran3(&raninit);
    ran3(&raninit);
  }
}

/* END pt_params_init()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up global arrays required for velocity interpolation from    */
/* the windows.                                                      */
/*-------------------------------------------------------------------*/
void pt_setup(master_t *master,
	      geometry_t **window,
	      window_t **windat,
	      win_priv_t **wincon
	      )
{
  geometry_t *geom = master->geom;
  int wn, cc, c, k;

  if (master->do_pt) {
    /* Get the global ghost cell mask from the windows               */
    geom->gcm = i_alloc_2d(geom->szcS, geom->nz+1);
    for (wn = 1; wn <= geom->nwindows; wn++) {
      for (cc = 1; cc < window[wn]->szcS; cc++) {
	c = window[wn]->wsa[cc];
	for (k = 0; k < window[wn]->nz; k++) {
	  geom->gcm[k][c] = window[wn]->gcm[k][cc];
	}
      }
    }
  }
}

/* END pt_setup()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Convert horizontal particle world coords to grid coords           */
/*-------------------------------------------------------------------*/
void ptgrid_xytoij(master_t *master, long np, particle_t *p)
{
  int n, i, j;
  for (n = 0; n < np; n++) {
    p[n].c = hd_grid_xytoij(master, p[n].e1, p[n].e2, &i, &j);
  }
}

/* END ptgrid_xytoij()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Convert horizontal particle grid coords to world coords           */
/*-------------------------------------------------------------------*/
void ptgrid_ijtoxy(master_t *master, long np, particle_t *p)
{
  int n;
  return;
}

/* END ptgrid_ijtoxy()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to activate more particles randomly along a line source   */
/* defined by the two points (pt_x1,pt_y1,pt_z1) and                 */
/* (pt_x2,pt_y2,pt_z2).                                              */
/*-------------------------------------------------------------------*/
void pt_new(master_t *master, long np, particle_t *p)
{
  int i = 0;
  int n = 0;
  int raninit = 0;

  for (i = 0; i < master->pt_nsource; ++i) {
    /* Number of particles introduced per particle step */
    double diffx = (master->pt_x2[i] - master->pt_x1[i]);
    double diffy = (master->pt_y2[i] - master->pt_y1[i]);
    double diffz = (master->pt_z2[i] - master->pt_z1[i]);
    double npart =
      master->pt_rate[i] * master->ptstep + master->pt_accum[i];
    int ptint = (int)npart;
    int count = 0;

    /* Accumulate the fractions of particles to be released */
    master->pt_accum[i] = npart - ptint;

    if(master->svel_type[i] & FILEIN) {
      pt_ts_t *w_s = master->wvel_s[i];
      /*Evaluates the time series of settling velocities*/
      master->pt_svel[i] = ts_eval_xyz(w_s->ts, w_s->id,
				       master->t, master->pt_x2[i],
				       master->pt_y2[i],
				       master->pt_z2[i]);
    }

    /* Process the new release particles */
    while (count < ptint) {
      /* Look for inactive or lost particles */

      if (!(p[n].flag & PT_ACTIVE) || (p[n].flag & PT_LOST)) {
        double r1 = ran3(&raninit);
        double r2 = ran3(&raninit);

        /* Having computed the frational hrizontal and vertical */
        /* distances, compute the XYZ location and convert the XY */
        /* into an appropraite IJ position.  */
        p[n].e1 = master->pt_x1[i] + r1 * diffx;
        p[n].e2 = master->pt_y1[i] + r1 * diffy;
        ptgrid_xytoij(master, 1, &master->ptp[n]);
        p[n].e3 = master->pt_z1[i] + r2 * diffz;

        /* Set flag to indicate that this particle is now active     */
        p[n].flag = PT_ACTIVE | master->pt_colour[i];
	p[n].dumpf = 0;

	/* Set the age equal to zero. Only the variable out_age is   */
	/* written to the output file and this is scaled to a short  */
	/* integer using the scaling agelim. This is done to reduce  */
	/* the size of the output file.                              */
	p[n].age = 0.0;
	p[n].out_age = 0;
	if (master->pt_agelim)
	  p[n].dumpf |= PT_AGE;

	/* Set the size equal to the initial size. Size is scaled as */
	/* per age to sizelim.                                       */
	p[n].size = master->pt_size[i];
	if(master->pt_decay[i] < 0.0)
	  p[n].out_size = MAX_COL;
	else
	  p[n].out_size = 0;
	if (master->pt_sizelim)
	  p[n].dumpf |= PT_SIZE;

	/* Initialise the settling velocity */
	p[n].svel = 0.0;

	/* Get the source map for this particle                      */
	master->pt_sm[n] = i;

        count++;
      }
      n++;
      if (n >= np) {
        hd_warn("pt_new: All particles already in use\n");
        return;
      }
    }
  }

}

/* END pt_new()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to activate more particles randomly along a line source   */
/* defined by the two points (pt_x1,pt_y1,pt_z1) and                 */
/* (pt_x2,pt_y2,pt_z2).                                              */
/*-------------------------------------------------------------------*/
void pt_split(master_t *master, long np, particle_t *p, 
	      double x, double y, double z)
{
  int i = 0;
  int n = 0;
  int raninit = 0;

  int count = 0;

  /* Process the new release particles */
  while (count < 1) {
    /* Look for inactive or lost particles */
    if (!(p[n].flag & PT_ACTIVE) || (p[n].flag & PT_LOST)) {
      double r1 = ran3(&raninit);
      double r2 = ran3(&raninit);

      /* Having computed the frational hrizontal and vertical */
      /* distances, compute the XYZ location and convert the XY */
      /* into an appropraite IJ position.  */
      p[n].e1 = x;
      p[n].e2 = y;
      ptgrid_xytoij(master, 1, &master->ptp[n]);
      p[n].e3 = z;

      /* Set flag to indicate that this particle is now active     */
      p[n].flag = PT_ACTIVE | master->pt_colour[i];
      p[n].dumpf = 0;

      /* Set the age equal to zero. Only the variable out_age is   */
      /* written to the output file and this is scaled to a short  */
      /* integer using the scaling agelim. This is done to reduce  */
      /* the size of the output file.                              */
      p[n].age = 0.0;
      p[n].out_age = 0;
      if (master->pt_agelim)
	p[n].dumpf |= PT_AGE;

      /* Set the size equal to the initial size. Size is scaled as */
      /* per age to sizelim.                                       */
      p[n].size = 1e-5;;
      p[n].out_size = 0;
      if (master->pt_sizelim)
	p[n].dumpf |= PT_SIZE;

      /* Initialise the settling velocity */
      p[n].svel = 0.0;

      count++;
    }
    n++;
    if (n >= np) {
      hd_warn("pt_split: All particles already in use\n");
      return;
    }
  }
}

/* END pt_split()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate concentrations from particle positions.      */
/*-------------------------------------------------------------------*/
void particles_to_conc(master_t *master, long np, particle_t *p)
{
  int c, n;
  int cc;
  double mass = master->pt_mass;
  double *dz = master->dz;
  geometry_t *geom = master->geom;

  /* Clear the particle concentration array */
  memset(master->ptconc, 0, geom->sgsiz * sizeof(double));

  /* Count the number of particles per cell.  */
  for (n = 0; n < np; n++) {
    unsigned long f = p[n].flag;
    c = p[n].c;
    if ((f & PT_ACTIVE) && (!(f & PT_LOST))) {
      master->ptconc[c] += mass;
    }
  }

  /* Calculate cell concentrations */
  for (cc = 1; cc < geom->b3_t; cc++) {
    int c = geom->w3_t[cc];
    int cs = geom->m2d[c];
    double a = geom->cellarea[cs];
    master->ptconc[c] /= a * dz[c];
  }
}
/* END particles_to_conc()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to mark only those particles in the water as active       */
/*-------------------------------------------------------------------*/
void pt_fit_to_grid(master_t *master, long np, particle_t *p)
{
  int n;

  for (n = 0; n < np; n++) {
    int c = p[n].c;
    int cs = geom->m2d[c];
    if (!c || p[n].e3 < geom->botz[cs] || p[n].e3 > master->topz[cs])
      p[n].flag &= ~PT_ACTIVE;
    else
      p[n].flag |= PT_ACTIVE;

    /* Set the age equal to zero. Only the variable out_age is       */
    /* written to the output file and this is scaled to a short      */
    /* integer using the scaling agelim. This is done to reduce      */
    /* the size of the output file.                                  */
    p[n].age = 0.0;
    p[n].out_age = 0;
    if (master->pt_agelim)
      p[n].dumpf |= PT_AGE;

    /* Set the flag whether a particle is within the age region      */
    if (master->ptmsk) {
      if(master->ptmsk[c]) 
	p[n].flag |= PT_IN;
      else
	p[n].flag |= PT_OUT;
    }
  }
}

/* END pt_fit_to_grid()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Copy particles and convert to real-world coordinates.             */
/*-------------------------------------------------------------------*/
void pt_write_at_t(master_t *master, double t)
{
  particle_t *pcopy;
  double newt = t;
  int n;
  if (DEBUG("particles"))
    dlog("particles", "Writing particles to file at t=%f\n", master->t);

  /* Write particles and free memory */
  tm_change_time_units(master->timeunit, master->output_tunit, &newt, 1);
  pt_write(master->ptfid, master->ptrec, newt, master->ptn, master->ptp);
  master->ptrec++;
  master->ptout_t += master->ptoutinc;
  if (master->ptout_t > master->ptend) {
    nc_close(master->ptfid);
    master->ptfid = -1;
  }
}

/* END pt_write_at_t()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to track all the particles for a single timestep          */
/*-------------------------------------------------------------------*/
void pt_update(master_t *master,
	       geometry_t **window,
	       window_t **windat,
	       win_priv_t **wincon
	       )
{
  geometry_t *geom = master->geom;

/*-------------------------------------------------------------------*/
  /* If no particles, or before start time, or output file not open, */
  /* nothing to do.                                                  */
  if (!master->ptn || master->t < master->ptstart || master->ptfid == -1)
    return;

  /*-----------------------------------------------------------------*/
  /* Move particles if required                                      */
  if (master->t + master->dt / 2 >= master->ptnext_t) {
    int n;
    int cc, c, c2;
    double x, y, z;
    double maxdh;
    pt_ts_t *u_i = master->uvel_i;
    pt_ts_t *v_i = master->vvel_i;
    pt_ts_t *w_i = master->wvel_i;
    pt_ts_t *m_i = master->mort_i;

    if (DEBUG("particles"))
      dlog("particles", "Moving particles, t=%.0f\n", master->t);

    /* Find min Kz values in each column, and scale by PT_KZ_MULT */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      int c = geom->w3_t[cc];
      int cs = geom->m2d[c];
      double min = 1;
      if (master->Kz[c] < min)
        min = master->Kz[c];
      master->maxdiffw[cs] =
        sqrt(6 * min * master->pt_kz_mult / master->ptstep);
    }

    /* Maximum possible horizontal diffusive velocity */
    maxdh = sqrt(6 * master->pt_kh / master->ptstep);

    /*---------------------------------------------------------------*/
    /* Set the settling velocity of each particle from an initial    */
    /* release if required.                                          */
    if(u_i && v_i) {
      double u, v;
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	c2 = geom->m2d[c];
	u = ts_eval_xyz(u_i->ts, u_i->id, master->t, 
			geom->cellx[c2], geom->celly[c2], 
			geom->cellz[c]);
	v = ts_eval_xyz(v_i->ts, v_i->id, master->t, 
			geom->cellx[c2], geom->celly[c2], 
			geom->cellz[c]);
	u_i->val[c] = u * geom->costhcell[c2] + v * geom->sinthcell[c2];
	v_i->val[c] = -u * geom->sinthcell[c] + v * geom->costhcell[c2];
      }
    }
    if(w_i) {
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	c2 = geom->m2d[c];
	w_i->val[c] = ts_eval_xyz(w_i->ts, w_i->id, master->t, 
				  geom->cellx[c2], geom->celly[c2], 
				  geom->cellz[c]);
      }
    }

    if (master->pt_nsource > 0)
      pt_new(master, master->ptn, master->ptp);

    /*---------------------------------------------------------------*/
    /* Set the settling velocity of each particle from a source if   */
    /* required.                                                     */
    if(master->pt_svel) {
      for (n = 0; n < master->ptn; n++) {
      	pt_set_svel(master, &master->ptp[n], master->ptstep, n);
      }
    }

    /*---------------------------------------------------------------*/
    /* Reads the mortality percentage from file if required          */
    if (m_i) {
      int part_numb = 0;
      int part_numb_act = 0;
      int part_numb_lost = 0;
      m_i->data = ts_eval(m_i->ts, m_i->id, master->t);
      for (n = 0; n < master->ptn; n++) {
	part_numb++;
	if (((master->ptp[n].flag) & PT_ACTIVE) &&
	    !((master->ptp[n].flag) & PT_LOST)) {
	  part_numb_act++;
	} else if(((master->ptp[n].flag) & PT_LOST)) {
	  part_numb_lost++;
	}
      }
      /* Make some of the particles lost due to mortality            */
      double dp = round((m_i->data /100)* part_numb_act);
      if(dp > 0){
	int interval = round(part_numb_act/dp);
	for (n = 0; n < master->ptn; n += interval) {
	  master->ptp[n].flag |= PT_LOST;
	}
      }
    }

    /* Initialize cell centered velocities                           */
    for (n = 1; n <= geom->nwindows; n++) {
      ffsl_init(window[n], windat[n], wincon[n]);
    }

    /*---------------------------------------------------------------*/
    /* Loop to move each particle                                    */
    for (n = 0; n < master->ptn; n++) {
      if (((master->ptp[n].flag) & PT_ACTIVE) &&
	  !((master->ptp[n].flag) & PT_LOST)) {
        pt_move(master, dumpdata, &master->ptp[n], master->ptstep, maxdh);
      }
    }

    /*---------------------------------------------------------------*/
    /* Set the age of each particle if required                      */
    if (master->pt_agelim) {
      for (n = 0; n < master->ptn; n++) {
	pt_set_age(master, &master->ptp[n], master->ptstep);
      }
    }

    /*---------------------------------------------------------------*/
    /* Set the size of each particle if required                     */
    if (master->pt_sizelim) {
      for (n = 0; n < master->ptn; n++) {
	pt_set_size(master, &master->ptp[n], master->ptstep, n);
      }
    }

    master->ptnext_t += master->ptstep;

    /*---------------------------------------------------------------*/
    /* Calculate new particle concentrations in model cells          */
    particles_to_conc(master, master->ptn, master->ptp);
  }

  /*-----------------------------------------------------------------*/
  /* Update the particle histogram if required                       */
  if (master->phist && master->t >= tsphist.tsout - DT_EPS) {
    int n;
    fprintf(tsphist.fp, "%f ", master->t / 86400);
    for (n = 0; n < master->shist; n++)
      fprintf(tsphist.fp, "%d ", master->phist[n]);
    fprintf(tsphist.fp, "\n");
    tsphist.tsout += tsphist.tsdt;
  }
}

/* END pt_update()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* This routine sets the age of all particles                        */
/*-------------------------------------------------------------------*/
void pt_set_age(master_t *master, /*  Master data                    */
		particle_t *p,    /* Pointer to particle to be moved */
		double dt         /* Time interval for which to move */
  )
{

  /* Set the age                                                     */
  if ((p->flag & PT_ACTIVE) && !(p->flag & PT_LOST)) {
    p->age += dt;
    p->out_age = min((short)(MAX_COL * p->age / master->pt_agelim),
		     MAX_COL);
  } else {
    /* Particle is lost - reset age to zero                          */
    p->age = 0.0;
    p->out_age = 0;
  }
}

/* END pt_set_age()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* This routine sets the size of all particles                       */
/*-------------------------------------------------------------------*/
void pt_set_size(master_t *master, /* Master data                    */
		 particle_t *p,    /* Pointer to particle            */
		 double dt,        /* Time interval for moving       */
		 int n             /* Particle number                */
  )
{
  double tol = 1e-10;
  double rate;

  /* Set the size                                                    */
  if ((p->flag & PT_ACTIVE) && !(p->flag & PT_LOST)) {
    /* If the size decreases below the minimum set the particle as   */
    /* lost.                                                         */
    rate = master->pt_decay[master->pt_sm[n]];
    if(p->size < tol) {
      p->flag |= PT_LOST;
    } else {
      /* Increase / decrease the size                                */
      p->size += p->size * dt / rate;
      if(rate < 0)
	p->out_size = max((short)(MAX_COL * p->size / master->pt_sizelim), 0);
      else {
	double size = master->pt_size[master->pt_sm[n]];
	p->out_size = max((short)(MAX_COL * (p->size - size) /
				  (master->pt_sizelim - size)), 0);
      }
    }
  } else {
    /* Particle is lost - reset size to zero                         */
    p->out_size = 0;
  }
}

/* END pt_set_size()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*This routine sets the settling velocity of all particles           */
/*-------------------------------------------------------------------*/
void pt_set_svel(master_t *master, /* Master data                    */
                 particle_t *p,    /* Pointer to particle            */
                 double dt,        /* Time interval for moving)      */
                 int n             /* Particle number                */
)

{
  int sm = master->pt_sm[n];         /* Particle source map          */
  int stype = master->svel_type[sm]; /* Settling velocity type       */
  double vel = master->pt_svel[sm];  /* Fixed settling velocity      */
  double per = master->pt_sper[sm];  /* Fixed settling period        */
  double dens = master->pt_dens[sm]; /* Water density in cell        */
  double mu = 1.4e-3; /* Molecular viscocity of water at 20C (kg/m/s)*/

  /* Set the settling type velocity                                  */
  if ((p->flag & PT_ACTIVE) && !(p->flag & PT_LOST)) {
    if (stype == NONE) {
      /* No settling                                                 */
      p->svel = 0.0;
    } else if (stype == CONSTANT) {
      /* Constant settling                                           */
      p->svel = vel;
    } else if (stype == PSTOKES) {
      /* Stokes settling                                             */
      int c = ztoc_m(master, floor(p->e1), p->e3);
      p->svel = -master->g * (dens - master->dens[c]) *
	p->size * p->size / (18.0 * mu);
    } else if (stype == DIURNAL) {
      /* Diurnal settling                                            */
      double frac = (master->t/(per))-(floor(master->t/(per)));
      p->svel = (-vel * cos(frac * 2*PI));
    } else if (stype == FILEIN) {
      /* Reading svel from a ts file                                 */
      p->svel = vel;
    }
  }
}

/* END pt_set_vel()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Values indicating which face a particle goes through              */
#define F_NONE     0
#define F_LEFT     1
#define F_RIGHT    2
#define F_BACK     3
#define F_FORWARD  4
#define F_BOTTOM   5
#define F_TOP      6

#define EPS         (1e-8)
#define TINY_INSIDE (1e-8)
/*-------------------------------------------------------------------*/
/* This routine moves a single particle for a specified time         */
/* The position of the particle is :                                 */
/* p->c = mesh index of the particle                                 */
/* p->e1 = x location of the particle                                */
/* p->e2 = y location of the particle                                */
/* p->e3 = depth of particle                                         */
/*-------------------------------------------------------------------*/
void pt_move(master_t *master,    /* Pointer to model grid structure */
             dump_data_t *dumpdata, /* Dump data structure           */
             particle_t *p, /* Pointer to the particle to be moved   */
             double dt,     /* The amount of time for which to move  */
             double maxdh   /* Maximum possible horizontal diffusion */
  )
{
  geometry_t *geom = master->geom;
  double diffu1;
  double diffu2;
  double diffw;
  double tleft = dt;
  double u, v, w;
  double upt, vpt, wpt;
  double cx, cy, cz;
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  int raninit = 0;
  int nloop = 0;
  int c, c2, co, wn, cl;
  pt_ts_t *u_i = master->uvel_i;
  pt_ts_t *v_i = master->vvel_i;
  pt_ts_t *w_i = master->wvel_i;
  int *c2cc = i_alloc_1d(geom->szcS);

  for (c2 = 1; c2 <= geom->n2_t; c2++) {
    c = geom->w2_t[c2];
    c2cc[c] = c2;
  }

  /*-----------------------------------------------------------------*/
  /* Find the horizontal (i,j) indices for the model cell that the   */
  /* particle is currently in, and check that it is sensible.        */
  c = co = p->c;
  c2 = geom->m2d[c];
  cx = p->e1;
  cy = p->e2;
  cz = p->e3;

  /*-----------------------------------------------------------------*/
  /* Convert the p.e1 sparse coordinate to Cartesian                 */
  if (c < 0 || c >= geom->szc) {
    /* Something wrong - print a message, flag this particle as      */
    /* lost, and give up. Note -9999 are the (i,j) coordinates for   */
    /* ghost cells.                                                  */
    hd_warn
      ("ptmove start: Particle not in water horizontally, %d %.12f %.12f\n",
       c, p->e1, p->e2);
    p->flag |= PT_LOST;
    return;
  }

  /* Surface may have moved since we last did a particle step, so    */
  /* truncate particles to surface.                                  */
  if (p->e3 > master->eta[c2])
    p->e3 = master->eta[c2];

  if (p->e3 < geom->botz[c2]) {
    /* Something wrong - print a message, flag this particle as lost */
    /* and give up.                                                  */
    if (DEBUG("particles"))
      dlog("particles", "Particle not in water vertically, c=%d k=%d x=%.12f y=%.12f z=%.12f\n",
	   c, geom->s2k[c], p->e1, p->e2, p->e3);
    p->flag |= PT_LOST;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate random diffusion velocities for this step             */
  diffu1 = maxdh * (2 * ran3(&raninit) - 1);
  diffu2 = maxdh * (2 * ran3(&raninit) - 1);
  diffw = master->maxdiffw[c2] * (2 * ran3(&raninit) - 1);

  /*-----------------------------------------------------------------*/
  /* Set the settling/swimming velocity for each particle from an    */
  /* initial release if required.                                    */
  if (u_i && v_i) {
    upt = u_i->val[c];
    vpt = v_i->val[c];
  } else {
    upt = 0.0;
    vpt = 0.0;
  }
  wpt = (w_i) ? w_i->val[c] : 0.0;

  /*-----------------------------------------------------------------*/
  /* Loop to move particle through grid. Each pass through this loop */
  /* (usually) moves the particle across one model cell.             */
  while (tleft > 0) {
    double tmin;
    double q, r, ua;
    int cl, wn;

    /* Get the local window amd coordinate                           */
    wn = geom->fm[c].wn;     /* Window for this particle             */
    cl = geom->fm[c].sc;     /* Local coordinate for this particle   */

    /*---------------------------------------------------------------*/
    /* Get the velocities at the particle location                   */
    u = hd_trans_interp(window[wn], windat[wn]->gsx, cx, cy, cz, cl, 0, 0);
    v = hd_trans_interp(window[wn], windat[wn]->gsy, cx, cy, cz, cl, 0, 1);
    w = hd_trans_interp(window[wn], windat[wn]->gsz, cx, cy, cz, cl, 0, 2);

    /*---------------------------------------------------------------*/
    /* Add diffusion velocities                                      */
    u = u + diffu1 + upt;
    v = v + diffu2 + vpt;
    w = w + diffw + wpt;

    /*---------------------------------------------------------------*/
    /* Get the sub time-step for the move                            */
    q = sqrt(geom->cellarea[c2]);
    ua = sqrt(u * u + v * v);
    r = (w < 0.0) ? wincon[wn]->dzz[c] : wincon[wn]->dzz[window[wn]->zm1[c]];
    tmin = SCALE * fabs(q / ua);
    tmin = (r) ? min(tmin, SCALE * fabs(r / w)) : tmin;

    /*---------------------------------------------------------------*/
    /* Get the new location and cell the particle resides in         */
    c = get_pos_m(master, p, co, c, u, v, w, &cx, &cy, &cz, tmin);
    c2 = geom->m2d[c];

    /*---------------------------------------------------------------*/
    /* Check if it landed in unknown territory                       */
    if (c <= 0 || c >= geom->szc) {
      hd_warn
	("ptmove end: Particle not in water horizontally, c=%d x=%.12f y=%.12f z=%.12f\n",c, p->e1, p->e2, p->e3);
      p->flag |= PT_LOST;
      return;
    }

    /* Check the particle didn't go through the bottom or surface    */
    if (cz <= geom->botz[c2]) {
      cz = geom->botz[c2] + TINY_INSIDE;
      c = geom->bot_t[c2cc[c2]];
    }
    if (cz >= master->eta[c2]) {
      cz = master->eta[c2] - TINY_INSIDE;
      c = geom->sur_t[c2cc[c2]];
    }

    /*---------------------------------------------------------------*/
    /* If the particle lands in a ghost cell, flag it as lost        */
    if (geom->wgst[c]) {
      p->flag |= PT_LOST;
      tmin = tleft;
      if (master->ptmsk == NULL) {
	master->mage += p->age;
	master->magec += 1.0;
	if(master->phist)
	  master->phist[(int)min(master->shist,
				 floor(p->age / 86400))]++;
      }
    }

    /*---------------------------------------------------------------*/
    /* Update the particle position                                  */
    p->e1 = cx;
    p->e2 = cy;
    p->e3 = cz;
    p->c = c;

    /*---------------------------------------------------------------*/
    /* Sanity checks                                                 */
    if (++nloop > 100) {
      if (DEBUG("particles")) {
	dlog("particles", "nloop = %d, tmin = %g, %d %.12f %.12f %.12f\n",nloop, tmin, c, p->e1, p->e2, p->e3);
      }
      if (nloop > 110) {
	if (DEBUG("particles"))
	  dlog("particles", "Abandoning this particle\n");
        p->flag |= PT_LOST;
        tmin = tleft;
      }
    }

    /*---------------------------------------------------------------*/
    /* Check whether it has left the age region                      */
    if (master->ptmsk) {
      if (p->flag & PT_IN && !master->ptmsk[c]) {
	master->mage += p->age;
	master->magec += 1.0;
	if(master->phist)
	  master->phist[(int)min(master->shist,
				 floor(p->age) / 86400)]++;
	p->flag &= ~ PT_IN;
	p->flag |= PT_OUT;
      }
      if (p->flag & PT_OUT && master->ptmsk[c]) {
	p->flag &= ~ PT_OUT;
	p->flag |= PT_IN;
	p->age = 0;
      }
    }

    /* Decrease time left for this step                              */
    tleft -= tmin;
  
  }
  i_free_1d(c2cc);
}

/* END pt_move()                                                     */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Computes the new location of the streamline origin. If this lands */
/* in a ghost cell, then shift the streamline in the wet interior    */
/* and re-compute a new location.                                    */
/* Destination cell is the origin of the streamline.                 */
/* Source cell is the end point of the streamline track.             */
/* Returns a new source cell after moving a distance increment.      */
/* The new geographic location is returned in cx, cy and cz.         */
/*-------------------------------------------------------------------*/
int get_pos_m(master_t *master,   /* Window geometry                 */
	      particle_t *p,      /* Particle                        */
	      int c,              /* Location of destination         */
	      int ci,             /* Current source cell             */
	      double u,           /* East velocity                   */
	      double v,           /* North velocity                  */
	      double w,           /* Vertical velocity               */
	      double *cx,         /* x location of streamline        */
	      double *cy,         /* y location of streamline        */
	      double *cz,         /* z location of streamline        */
	      double dt           /* Time step                       */
	      )
{
  geometry_t *geom = master->geom;
  int cs = geom->m2d[c];
  int cis = geom->m2d[ci];
  int cn, cns, zm1;
  double xin = *cx;
  double yin = *cy;
  double zin = *cz;
  double dist, dir, s1, s2;
  double slat, slon;
  double sinth, costh;
  double dx, dy;
  double nx, ny;
  double m2deg = 1.0 / (60.0 * 1852.0);
  int isghost = 0;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int j, k = geom->s2k[ci];
  int found = 0;

  /* Get the new horizontal location                                 */
  if (!is_geog) m2deg = 1.0;
  if (sodanos) {
    dist = sqrt(u * u + v * v) * dt;
    dir =  PI / 2.0 - atan2(v, u);  /* deg T                         */
    geod_fwd_sodanos(DEG2RAD(xin), DEG2RAD(yin), dir, dist, RADIUS, ECC, &slon, &slat);
    slon = RAD2DEG(slon);
    slat = RAD2DEG(slat);
  } else {
    /* Get the new distances                                         */
    dist = sqrt(u * u + v * v) * dt * m2deg;
    dir =  atan2(v, u);
    /* Update the origin location                                    */
    sinth = sin(dir);
    costh = cos(dir);
    slon = xin - dist * costh;
    slat = yin - dist * sinth;
  }

  /* Get the cell the new horizontal location resides in, i.e. at    */
  /* the same depth as that of the input location, c.                */
  if (geom->us_type & US_IJ) {
    /* Structured grids: use the xytoi tree                          */
    int i = -1, j = -1;
    if (grid_xytoij(master->xyij_tree, slon, slat, &i, &j))
      cn = geom->map[geom->nz - 1][j][i];
    else {
      cn = geom->m2d[ci];
      isghost = 1;
    }
    cns = geom->m2d[cn];
  } else {
    /* Unstructured meshes: walk through the Voronoi mesh            */
    cn = found = find_cell(geom, c, slon, slat, &nx, &ny);
    cns = geom->m2d[cn];
  }

  /* Get the vertical layer of the source cell                       */
  if (cn == geom->zm1[cn]) cn = geom->zp1[cn];
  if (geom->wgst[cn]) cn = geom->wgst[cn];
  k = geom->s2k[cn];
  if (k > geom->nz - 1) 
    hd_quit("Can't find streamline position: destination = %d[%f %f]\n",
	    c, geom->cellx[cs], geom->celly[cs]);

  /* If the new location ixs a ghost cell, send it inside the mesh.   */
  if (geom->gcm[k][cns] & L_GHOST) isghost = 1;
  if (isghost) {
    if (master->pt_stickybdry) {
      if (geom->us_type & US_IJ) {
	int es, vs;
	for (j = 1; j <= geom->npe[cis]; j++) {
	  if ((found = intersect(geom, cis, j, xin, yin, slon, slat, &nx, &ny)))
	    break;
	}
	if (!found) {
	  nx = slon;
	  ny = slat;
	}
      }
      dx = nx - xin;
      dy = ny - yin;
      dist = sqrt(dx * dx + dy * dy) - TINY_INSIDE;
      slon = xin - dist * costh;
      slat = yin - dist * sinth;

      if (geom->us_type & US_IJ) {
	int i = -1, j = -1;
	if (grid_xytoij(master->xyij_tree, slon, slat, &i, &j))
	  cn = geom->map[geom->nz - 1][j][i];
	else
	  cn = geom->m2d[ci];
	cns = geom->m2d[cn];
      } else {
	cns = found = find_cell(geom, c, slon, slat, &nx, &ny);
	cns = geom->m2d[cns];
      }
    } else {
      p->flag |= PT_LOST;
      if (master->ptmsk == NULL) {
	master->mage += p->age;
	master->magec += 1.0;
	if(master->phist)
	  master->phist[(int)min(master->shist,
				 floor(p->age / 86400))]++;
      }
    }
  }
  if (isghost && !found) p->flag |= PT_LOST;
  *cx = slon;
  *cy = slat;
  *cz -= w * dt;

  /* Find the vertical position in the mesh                          */
  cn = cns;
  /* First get the layer of the free surface                         */
  while (master->eta[cns] < geom->gridz[cn] && c != geom->zm1[cn]) {
    cn = geom->zm1[cn];
  }
  /* Get the first cellz layer less than z                           */
  while (*cz <= max(geom->botz[cns], geom->gridz[cn]) && cn != geom->zm1[cn]) {
    cn = geom->zm1[cn];
  }
  return(cn);
}

/* END get_pos_m()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Create particle tracking output file.                             */
/*-------------------------------------------------------------------*/
int hd_pt_create(master_t *master,
		 char *name)               /* name particle filename */
{
  int ncid;                     /* netCDF id */
  /* dimension ids */
  int t_dim, n_dim, m_dim, l_dim;
  /* variable ids */
  int t_id, x_id, y_id, z_id, flag_id;
  int age_id, size_id, source_id, vid;
  /* variable shapes */
  int n, dims[2];
  int np = master->ptn;
  char *t_units = master->output_tunit;
  int dumpf = master->pt_dumpf;
  int age = master->pt_agelim;
  int size = master->pt_sizelim;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  size_t start[4];
  size_t count[4];

  /* enter define mode */
  if (nc_create(name, NC_NOCLOBBER, &ncid) != NC_NOERR)
    quit("pt_create: Unable to create the particle file '%s'\n", name);

  /* define dimensions */
  nc_def_dim(ncid, "t", NC_UNLIMITED, &t_dim);
  nc_def_dim(ncid, "n", np, &n_dim);
  nc_def_dim(ncid, "m", master->pt_nsource, &m_dim);
  nc_def_dim(ncid, "l", 2, &l_dim);

  /* define variables */
  dims[0] = t_dim;
  nc_def_var(ncid, "t", NC_DOUBLE, 1, dims, &t_id);
  dims[1] = n_dim;
  nc_def_var(ncid, "x", NC_DOUBLE, 2, dims, &x_id);
  vid = ncw_var_id(ncid, "x");
  if (is_geog) {
    write_text_att(ncid, vid, "long_name", "X coordinate of particle");
    write_text_att(ncid, vid, "coordinate_type", "longitude");
    write_text_att(ncid, vid, "units", "degrees_east");
  } else {
    write_text_att(ncid, vid, "long_name", "X coordinate of particle");
    write_text_att(ncid, vid, "coordinate_type", "X");
    write_text_att(ncid, vid, "units", "metre");
    }
  if (has_proj) write_text_att(ncid, vid, "projection", projection);
  nc_def_var(ncid, "y", NC_DOUBLE, 2, dims, &y_id);
  vid = ncw_var_id(ncid, "y");
  if (is_geog) {
    write_text_att(ncid, vid, "long_name", "Y coordinate of particle");
    write_text_att(ncid, vid, "coordinate_type", "latitude");
    write_text_att(ncid, vid, "units", "degrees_north");
  } else {
    write_text_att(ncid, vid, "long_name", "Y coordinate of particle");
    write_text_att(ncid, vid, "coordinate_type", "Y");
    write_text_att(ncid, vid, "units", "metre");
    }
  if (has_proj) write_text_att(ncid, vid, "projection", projection);
  nc_def_var(ncid, "z", NC_DOUBLE, 2, dims, &z_id);
  vid = ncw_var_id(ncid, "z");
  write_text_att(ncid, vid, "long_name", "Z coordinate of particle");
  write_text_att(ncid, vid, "units", "m");
  write_text_att(ncid, vid, "coordinate_type", "Z");
  write_text_att(ncid, vid, "positive", "up");

  nc_def_var(ncid, "flag", NC_SHORT, 2, dims, &flag_id);
  vid = ncw_var_id(ncid, "flag");
  write_text_att(ncid, vid, "long_name", "Particle status flag");

  if (dumpf & PT_AGE) {
    strcpy(key, "age");
    nc_def_var(ncid, key, NC_BYTE, 2, dims, &age_id);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "long_name", "Age of particle since release");
    write_text_att(ncid, vid, "units", "scaled byte");
    sprintf(buf, "%d seconds", age);
    write_text_att(ncid, vid, "Agelimit", buf);
  }
  if (dumpf & PT_SIZE) {
    strcpy(key, "size");
    nc_def_var(ncid, key, NC_BYTE, 2, dims, &size_id);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "long_name", "Size of particle");
    write_text_att(ncid, vid, "units", "scaled byte");
    sprintf(buf, "%d m", size);
    write_text_att(ncid, vid, "Sizelimit", buf);
  }
  /* Time variable units attribute */
  nc_put_att_text(ncid, t_id, "units", strlen(t_units), t_units);

  if (master->pt_nsource) {
    dims[0] = l_dim;
    dims[1] = m_dim;
    strcpy(key, "source_locationX");
    nc_def_var(ncid, key, NC_DOUBLE, 2, dims, &vid);
    vid = ncw_var_id(ncid, key);
    if (is_geog) {
      write_text_att(ncid, vid, "long_name", "Longitude of source");
      write_text_att(ncid, vid, "coordinate_type", "longitude");
      write_text_att(ncid, vid, "units", "degrees_east");
    } else {
      write_text_att(ncid, vid, "long_name",
		     "X coordinate of source");
      write_text_att(ncid, vid, "coordinate_type", "X");
      write_text_att(ncid, vid, "units", "metre");
    }
    if (has_proj)
      write_text_att(ncid, vid, "projection", projection);
    strcpy(key, "source_locationY");
    nc_def_var(ncid, key, NC_DOUBLE, 2, dims, &vid);
    vid = ncw_var_id(ncid, key);
    if (is_geog) {
      write_text_att(ncid, vid, "long_name", "Latitude of source");
      write_text_att(ncid, vid, "coordinate_type", "latitude");
      write_text_att(ncid, vid, "units", "degrees_north");
    } else {
      write_text_att(ncid, vid, "long_name",
		     "Y coordinate of source");
      write_text_att(ncid, vid, "coordinate_type", "Y");
      write_text_att(ncid, vid, "units", "metre");
    }
    if (has_proj)
      write_text_att(ncid, vid, "projection", projection);
    strcpy(key, "source_locationZ");
    nc_def_var(ncid, key, NC_DOUBLE, 2, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "long_name", "Depth of source");
    write_text_att(ncid, vid, "coordinate_type", "Z");
    write_text_att(ncid, vid, "units", "metre");
    dims[0] = m_dim;
    strcpy(key, "source_rate");
    nc_def_var(ncid, key, NC_DOUBLE, 1, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "units", "pt/hr");
    strcpy(key, "source_settling_velocity");
    nc_def_var(ncid, key, NC_DOUBLE, 1, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "units", "ms-1");
    write_text_att(ncid, vid, "positive", "up");
    strcpy(key, "source_density");
    nc_def_var(ncid, key, NC_DOUBLE, 1, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "units", "kgm-3");
    strcpy(key, "source_size");
    nc_def_var(ncid, key, NC_DOUBLE, 1, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "units", "m");
    strcpy(key, "source_decay");
    nc_def_var(ncid, key, NC_DOUBLE, 1, dims, &vid);
    vid = ncw_var_id(ncid, key);
    write_text_att(ncid, vid, "units", "s-1");
    strcpy(key, "source_colourbit");
    nc_def_var(ncid, key, NC_INT, 1, dims, &vid);
  }

  write_text_att(ncid, NC_GLOBAL, "title", codeheader);
  write_text_att(ncid, NC_GLOBAL, "paramhead", parameterheader);
  getcwd(buf, MAXSTRLEN);
  sprintf(buf, "%s/%s", buf, master->prmname);
  write_text_att(ncid, NC_GLOBAL, "paramfile", buf);
  write_text_att(ncid, NC_GLOBAL, "version", version);
  write_text_att(ncid, NC_GLOBAL, "Conventions", "CMR/Timeseries/SHOC");
  write_text_att(ncid, NC_GLOBAL, "InputFile", master->ptinname);
  sprintf(buf,"%6.2f days", master->ptstart/86400.0);
  write_text_att(ncid, NC_GLOBAL, "StartTime", buf);
  sprintf(buf,"%6.2f days", master->ptend/86400.0);
  write_text_att(ncid, NC_GLOBAL, "StopTime", buf);
  sprintf(buf,"%6.4f hours", master->ptoutinc/3600.0);
  write_text_att(ncid, NC_GLOBAL, "OutputStep", buf);
  sprintf(buf,"%6.4f hours", master->ptstep/3600.0);
  write_text_att(ncid, NC_GLOBAL, "TimeStep", buf);
  sprintf(buf,"%6.2f days", master->ptreset/86400.0);
  write_text_att(ncid, NC_GLOBAL, "ResetStep", buf);
  nc_put_att_double(ncid, NC_GLOBAL, "KH", NC_DOUBLE, 1, &master->pt_kh);
  nc_put_att_double(ncid, NC_GLOBAL, "KZ_MULT", NC_DOUBLE, 1, &master->pt_kz_mult);
  nc_put_att_double(ncid, NC_GLOBAL, "Mass", NC_DOUBLE, 1, &master->pt_mass);
  nc_put_att_int(ncid, NC_GLOBAL, "StickyBoundary", NC_INT, 1, &master->pt_stickybdry);
  write_text_att(ncid, NC_GLOBAL, "flag:PT_ACTIVE", "0x001, particle is active");
  write_text_att(ncid, NC_GLOBAL, "flag:PT_LOST", "0x002, particle is lost");
  write_text_att(ncid, NC_GLOBAL, "flag:PT_AGE", "0x004, particle has age attribute");
  write_text_att(ncid, NC_GLOBAL, "flag:PT_SIZE", "0x008, particle has size attribute");

  /* leave define mode */
  nc_enddef(ncid);

  if (master->pt_nsource) {
    int *cbit;
    double **xloc, **yloc, **zloc;
    xloc = d_alloc_2d(master->pt_nsource, 2);
    yloc = d_alloc_2d(master->pt_nsource, 2);
    zloc = d_alloc_2d(master->pt_nsource, 2);
    cbit = i_alloc_1d(master->pt_nsource);
    for (n = 0; n < master->pt_nsource; n++) {
      cbit[n] = (int)master->pt_colour[n];
      xloc[0][n] = master->pt_x1[n];
      xloc[1][n] = master->pt_x2[n];
      yloc[0][n] = master->pt_y1[n];
      yloc[1][n] = master->pt_y2[n];
      zloc[0][n] = master->pt_z1[n];
      zloc[1][n] = master->pt_z2[n];
    }
    start[0] = 0;
    start[1] = 0;
    count[0] = 2;
    count[1] = master->pt_nsource;
    nc_d_writesub_2d(ncid, ncw_var_id(ncid, "source_locationX"), start, count, xloc);
    nc_d_writesub_2d(ncid, ncw_var_id(ncid, "source_locationY"), start, count, yloc);
    nc_d_writesub_2d(ncid, ncw_var_id(ncid, "source_locationZ"), start, count, zloc);
    start[0] = 0;
    count[0] = master->pt_nsource;
    nc_d_writesub_1d(ncid, ncw_var_id(ncid, "source_rate"), start, count,
		     master->pt_rate);
    nc_d_writesub_1d(ncid, ncw_var_id(ncid, "source_settling_velocity"), start, count,
		     master->pt_svel);
    nc_d_writesub_1d(ncid, ncw_var_id(ncid, "source_density"), start, count,
		     master->pt_dens);
    nc_d_writesub_1d(ncid, ncw_var_id(ncid, "source_size"), start, count,
		     master->pt_size);
    nc_d_writesub_1d(ncid, ncw_var_id(ncid, "source_decay"), start, count,
		     master->pt_decay);
    nc_i_writesub_1d(ncid, ncw_var_id(ncid, "source_colourbit"), start, count,
		     cbit);
    d_free_2d(xloc);
    d_free_2d(yloc);
    i_free_1d(cbit);
  }
  return (ncid);
}

/* END hd_pt_create()                                                */
/*-------------------------------------------------------------------*/

