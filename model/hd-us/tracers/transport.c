/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/tracers/transport.c
 *  
 *  Description: Transport model routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: transport.c 6479 2020-02-18 23:54:44Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define SMALL 1e-10             /* Minimum value for velocity bounds */
#define RADIUS 6370997.0
#define ECC 0.0
#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

int sls[5] = {8, 8, 27, 64, 125};

void global_fill(geometry_t *window, window_t *windat, 
		 win_priv_t *wincon, double *dtracer, int n, int mode);
void monotonic_fill(geometry_t *window, window_t *windat, 
		    win_priv_t *wincon, int n, int mode);
void global_tfill(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, double *dtracer, int n);
double max_lag(geometry_t *window, double *tr, int c);
double min_lag(geometry_t *window, double *tr, int c);
void check_tracer_eta(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_tracer_eta(geometry_t *window, window_t *windat, 
		    win_priv_t *wincon);
double obc_tr_mass(geometry_t *window, window_t *windat, 
		   win_priv_t *wincon, int n, int ord, int vecf);
double obc_flux(geometry_t *window, window_t *windat, 
		win_priv_t *wincon, int vecf);
void obc_scale(geometry_t *window, window_t *windat, 
	       win_priv_t *wincon, double msf, int vecf);
double obc_msf(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       double obcvol, int *vec, int nvec, int vecf, double *msf);
void check_warn(geometry_t *geom, char *var, double val, double t, 
		int c, int dim);
void trans_data_check(master_t *master, geometry_t **window,
		      window_t **windat, win_priv_t **wincon);
void trans_data_nan(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_dzt(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_merged_dz(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_lateral_OBC_tr(geometry_t *window, int ntr, double **tr);
void calc_volume(geometry_t *window, window_t *windat, win_priv_t *wincon);
void check_transport(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_means_t(geometry_t *window, window_t *windat, win_priv_t *wincon);
void recalc_vel(geometry_t *window, window_t *windat, win_priv_t *wincon);
void flux_to_velocity(geometry_t *window, window_t *windat, win_priv_t *wincon);
void merge_volflux(geometry_t *window, window_t *windat, win_priv_t *wincon);
int cg;
void vel_center_w(geometry_t *window, window_t *windat, win_priv_t *wincon, double *nw);
void delaunay_reinit(geometry_t *window, delaunay **d, int vid, double *v);
void grid_spec_init_tran(geometry_t *window, window_t *windat, win_priv_t *wincon,
			 int vid);
void set_trans_sed(geometry_t *window, int vid, double *v);
void weights_v(geometry_t *window, int co, double rin,double *dzz, double *wgt, int osl);
void print_tri_k(delaunay *d, int k);
void print_source_k(geometry_t *window, int nvec, int *vec, double *cx, double *cy, 
		    int k, int *s2k);
void set_tr_nograd(geometry_t *window, double *tr);
void set_interp_points(geometry_t *window, window_t *windat, win_priv_t *wincon,
		       delaunay **d, int size);
void set_interp_square(geometry_t *window, window_t *windat, win_priv_t *wincon,
		       delaunay **d);
void set_interp_ring(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     delaunay **d);
void set_interp_points_n(geometry_t *window, window_t *windat, win_priv_t *wincon,
			 delaunay **d, npt_t **npt, int size);
void set_interp_square_n(geometry_t *window, window_t *windat, win_priv_t *wincon,
			 delaunay **d, npt_t **npt);
void set_interp_ring_n(geometry_t *window, window_t *windat, win_priv_t *wincon,
		       delaunay **d, npt_t **npt);
void reset_npt(geometry_t *window, delaunay **d, npt_t **npt);
npt_t* npt_create(int n);
int get_pos(geometry_t *window, int c, int ci, double u, double v, double w,
            double *cx, double *cy, double *cz, double dt, int mode);
int get_posc(geometry_t *window, int *ei, int ci, double u, double v, 
            double *cx, double *cy, double dt, int mode, double *odist);
int find_cell(geometry_t *window, int c, double x, double y, double *xi, double *yi);
double get_dist(double x1, double y1, double x2, double y2);
void ffsl_doc(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	      double *ntr, double dtu, int mode);
void check_monotone(geometry_t *window, window_t *windat, win_priv_t *wincon, double tr, int c);
void universal_limit(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     int cin, int co, int cs, int cb, double dt, double *tr, double *Fx, double *Fz, int *mask);

/*-------------------------------------------------------------------*/
/* Transport step                                                    */
/*-------------------------------------------------------------------*/
void transport_step(master_t *master, geometry_t **window,
		    window_t **windat, win_priv_t **wincon)
{
  geometry_t *geom = master->geom;
  int nwindows = master->nwindows;
  geometry_t *tpg;
  window_t *tpd;
  win_priv_t *tpc;
  int nn, n, s, nb;
  int cc, c;
  double clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Reset the means on the master if required.                      */
  /* This is now done within a scheduled function
     reset_means_m_o(master);*/

  /*-----------------------------------------------------------------*/
  /* Do the custom tracer routines on the master                     */
  TIMING_SET;
  bdry_eval_tr_m(geom, master);
  TIMING_DUMP(2, "  bdry3d_tr ");

  /* Read in boundary velocities for STREAMLINE mode; used to get    */
  /* fluxes for global fills.                                        */
  if (master->tmode & DO_OBC || master->fillf & SET_BDRY) {
    TIMING_SET;
    bdry_eval_u1_m(geom, master);
    TIMING_DUMP(2,"  bdry_eval_u1_u2");
  }
  if (master->fillf & SET_BDRY_ETA) {
    TIMING_SET;
    bdry_eval_eta_m(geom, master);
    TIMING_DUMP(2,"  bdry_eval_eta");
  }

  /*-----------------------------------------------------------------*/
  /* Solve the tracer equation in each window                        */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];

    /*---------------------------------------------------------------*/
    /* Set pointers                                                  */
    tpg = window[n];
    tpd = windat[n];
    tpc = wincon[n];
    wincon[n]->togn = master->togn;

    /*---------------------------------------------------------------*/
    /* Interpolate onto the target grid (master)                     */
    memcpy(master->eta, tpd->eta, tpg->szcS * sizeof(double));
    memcpy(master->Kz, tpd->Kz, tpg->sgsiz * sizeof(double));
    /* Velocities are required for MONOTONIC global fills          */
    if (master->fillf & MONOTONIC || master->smagorinsky > 0.0) {
      memcpy(master->u1, tpd->u1, tpg->sgsiz * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Fill 3D  velocities into the window data structures.          */
    win_data_fill_3d(master, window[n], windat[n], master->nwindows);

    /*---------------------------------------------------------------*/
    /* Calculate a heat and salt flux if required                    */
    calc_heatf(window[n], windat[n], wincon[n]);
    calc_saltf(window[n], windat[n], wincon[n]);
    
    init_sigma(window[n], windat[n], wincon[n]);
    trans_data_nan(window[n], windat[n], wincon[n]);

    /* Calculate the 3D velocity alert diagnostics                   */
    if (!(master->alertf & NONE)) {
      vel2D_lbc(windat[n]->u1, window[n]->nbpte1, window[n]->nbe1,
		window[n]->bpte1, window[n]->bine1, wincon[n]->slip);
      alerts_w(window[n], VEL3D);
    }

    /*---------------------------------------------------------------*/
    /* Read the boundary velocities from file for global fills if    */
    /* required.                                                     */
    if (master->tmode & DO_OBC) {
      for (nb = 0; nb < window[n]->nobc; nb++)
	reset_bdry_eta(window[n], windat[n], wincon[n], window[n]->open[nb], windat[n]->eta);
      windat[n]->nu1 = wincon[n]->w9;
      /* Transfer any custom data from the master to the slaves      */
      bdry_transfer_u1(master, window[n], windat[n]);
      memcpy(windat[n]->nu1, windat[n]->u1, window[n]->sgsiz * sizeof(double));
      bdry_u1_3d(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->u1, windat[n]->nu1, window[n]->sgsiz * sizeof(double));
    }
    if (master->fillf & SET_BDRY) {
      windat[n]->nu1 = wincon[n]->w9;
      /* Transfer any custom data from the master to the slaves      */
      bdry_transfer_u1(master, window[n], windat[n]);
      memset(windat[n]->nu1, 0, window[n]->sgsiz * sizeof(double));
      bdry_u1_3d(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->u1, windat[n]->nu1, window[n]->sgsiz * sizeof(double));
    }
    if (master->fillf & SET_BDRY_ETA) {
      wincon[n]->neweta = wincon[n]->d4;
      bdry_transfer_eta(master, window[n], windat[n]);
      memset(wincon[n]->neweta, 0, window[n]->sgsizS * sizeof(double));
      bdry_eta(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->eta, wincon[n]->neweta, window[n]->sgsizS * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Set the cell centered velocities. Note: tangential velocity   */
    /* is not read from file and must be computed.                   */
    vel_cen(window[n], windat[n], wincon[n], windat[n]->u1, NULL,
	    windat[n]->u, windat[n]->v, NULL, NULL, 0);

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.  */
#if !GLOB_BC
    set_lateral_OBC_tr(window[n], windat[n]->ntr, windat[n]->tr_wc);
    set_lateral_BC_tr(windat[n]->tr_wc, windat[n]->ntr, window[n]->nbpt,
                      window[n]->bpt, window[n]->bin);
    set_lateral_BC_vel(windat[n]->u1flux3d, window[n]->nbpt,
                       window[n]->bpte1, window[n]->bine1);
#endif

    /*---------------------------------------------------------------*/
    /* Set the vertical  mixing coefficients. This  is done  in this */
    /* window loop so that transferred velocities from other windows */
    /* can be used in the velocity shear term.                       */
    if (wincon[n]->do_closure) {
      density_w(window[n], windat[n], wincon[n]);
      wincon[n]->calc_closure(window[n], windat[n], wincon[n]);
      bdry_closure(window[n], windat[n], wincon[n]);
    }

    /*---------------------------------------------------------------*/
    /* Get the surface boundary condition for vertical velocity.     */
    if(wincon[n]->tgrid & (EXACT|SUBSET|INEXACT)) {
      tpd->t = master->t;
      tpd->dt = master->dt;
      tpd->dtf = master->dtf;
      tpd->dtb = master->dtb;
      tpd->dttr = master->dttr;
      tpd->rampval = master->rampval;
      tpd->nstep = master->nstep;
      /*tpd->etarlxtc = master->etarlxtc;*/
      tpd->df_diagn_set = master->df_diagn_set;
      if (!(master->tmode & SP_ORIGIN)) {
	for (cc = 1; cc <= tpg->b2_t; cc++) {
	  c = tpg->w2_t[cc];
	  tpd->detadt[c] = (tpd->eta[c] - tpd->etab[c]) / tpd->dttr;	
	}
	vel_w_bounds(tpg, tpd, tpc);
	memcpy(tpd->etab, tpd->eta, tpg->sgsizS * sizeof(double));
      }
    } else {
      if (!(master->tmode & SP_ORIGIN)) {
	for (cc = 1; cc <= window[n]->b2_t; cc++) {
	  c = window[n]->w2_t[cc];
	  windat[n]->detadt[c] = (windat[n]->eta[c] - windat[n]->etab[c]) / 
	    windat[n]->dttr;
	}
	vel_w_bounds(window[n], windat[n], wincon[n]);
	memcpy(windat[n]->etab, windat[n]->eta, 
	       window[n]->sgsizS * sizeof(double));
      }
    }

    /*---------------------------------------------------------------*/
    /* When using the flux-form semi-lagrange scheme, the velocities */
    /* must be set to zero on the first time step of the transport   */
    /* model run. On the first time step, the transport model uses   */
    /* the same value for the old elevation (oldeta) and the new     */
    /* elevation (eta), so in effect the sea surface elevation does  */
    /* not change on the first time step. The flux divergence        */
    /* therefore also needs to be zero for the first time step only, */
    /* which is achieved by setting u1, u2, w and volume fluxes = 0. */
    if ((wincon[n]->trasc == FFSL) && (windat[n]->nstep == 0)) {
      hd_warn("Setting initial velocity and volume fluxes to zero for FFSL advection scheme (t = %f)\n", windat[n]->t / 86400);
      memset(windat[n]->u1, 0, window[n]->sze * sizeof(double));
      memset(windat[n]->w, 0, window[n]->szc * sizeof(double));
      memset(windat[n]->u1vm, 0, window[n]->sze * sizeof(double));
      memset(windat[n]->u1flux3d, 0, window[n]->sze * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Evalulate the sources and sinks of water                      */
    /* If using the flux-form semi-lagrange advection scheme in the  */
    /* transport model, don't call ss_water. The volume flux due to  */
    /* the point source should have been stored (as a vertical       */
    /* velocity) in wmean and transferred above to waterss.          */
    if (wincon[n]->npss)
      ss_water(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Get the cfl time-steps if required                            */
    if (!(wincon[n]->cfl & NONE))
    calc_cfl(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Get the mean velocity and elevation if required               */
    calc_means_t(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Increment the mean counter if required                        */
    reset_means(window[n], windat[n], wincon[n], RESET);

    /*-----------------------------------------------------------------*/
    /* Get the initial total mass of tracer using cells to process and */
    /* dz corresponding to the previous time-step.                     */
    calc_volume(window[n], windat[n], wincon[n]);

    /*-----------------------------------------------------------------*/
    /* Fill options:                                                   */
    /* keyword = GLOBAL, flag = GLOBAL; non-monotonic global fill      */
    /* keyword = MONOTONIC, flag = MONOTONIC; monotonic global fill    */
    /* keyword = OBC_ADJUST, flag |= OBC_ADJUST; adjust OBC fluxes for */
    /*           global fills.                                         */
    /* keyword = WEIGHTED, flag = WEIGHTED; monotonic global fill      */
    /* keyword = DIAGNOSE, flag |= DIAGNOSE; print diagnostics to file */
    /*           trans.ts                                              */
    /* keyword = DIAGNOSE_BGC, flag |= DIAGNOSE_BGC; print diagnostics */
    /*           to file trans.ts                                      */
    if (!(wincon[n]->fillf & NONE))
      global_fill(window[n], windat[n], wincon[n], NULL, -1, 0);

    /*---------------------------------------------------------------*/
    /* Set up non semi-lagrange advection schemes                    */
    if (!(wincon[n]->trasc & LAGRANGE)) {
      if (wincon[n]->conserve & CONS_MRG)
	merge_volflux(window[n], windat[n], wincon[n]);

      /* Set the vertical grid spacings at e1 and e2 faces.          */
      set_dz_at_u1(window[n], windat[n], wincon[n]);
      set_flux_3d(window[n], windat[n], wincon[n], VEL3D);

      /* Compute eta based on depth averaged velocity divergence.    */
      /* This corresponds to the end of the time-step.               */     
      if (wincon[n]->conserve & CONS_ETA) {
	if (wincon[n]->trasc == FFSL)
	  check_tracer_eta(window[n], windat[n], wincon[n]);
	else
	  calc_tracer_eta(window[n], windat[n], wincon[n]);
      }
      /* Recompute the fluxes and velocity                           */
      if (wincon[n]->conserve & CONS_SUB)
	recalc_vel(window[n], windat[n], wincon[n]);

      /* Compute the vertical velocity based on velocity divergence. */
      /* This makes velocity and elevation dynamically consistent,   */
      /* hence advection is conservative.                            */  
      set_dz(window[n], windat[n], wincon[n]);
      if (wincon[n]->conserve & CONS_W) {
	if (wincon[n]->trasc == FFSL)
	  ff_sl_w_update(window[n], windat[n], wincon[n]);
	else
	  vel_w_update(window[n], windat[n], wincon[n]);
      }
    } else {
      set_dz_at_u1(window[n], windat[n], wincon[n]);
      /* Set up the cell centered surface vectors and dz arrays      */
      set_dz(window[n], windat[n], wincon[n]);
      /* Compute eta based on depth averaged velocity divergence.    */
      /* This corresponds to the end of the time-step.               */     
      if (wincon[n]->conserve & CONS_W) {
	set_flux_3d(window[n], windat[n], wincon[n], VEL3D);
	calc_tracer_eta(window[n], windat[n], wincon[n]);
	vel_w_update(window[n], windat[n], wincon[n]);
      }

      /* Set the dz arrays for eta at t+1 (nsur_t)                   */
      reset_dzt(window[n], windat[n], wincon[n]);
      memcpy(wincon[n]->s1, wincon[n]->s2, window[n]->sgsiz * sizeof(int));
      wincon[n]->vc = wincon[n]->vc2;
      wincon[n]->vcs = wincon[n]->vcs2;
      wincon[n]->vcs1 = wincon[n]->vca2;
      wincon[n]->vc1 = wincon[n]->vci2;
      if(wincon[n]->tgrid & (EXACT|SUBSET|INEXACT)) {
	set_dz(tpg, tpd, tpc);
	reset_dzt(tpg, tpd, tpc);
	memcpy(tpc->s1, tpc->s2, tpg->sgsiz * sizeof(int));
	tpc->vc = tpc->vc2;
	tpc->vcs = tpc->vcs2;
	tpc->vcs1 = tpc->vca2;
	tpc->vc1 = tpc->vci2;
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the Smagorinsky mixing if required                        */
    if (wincon[n]->smagorinsky != 0.0) {
      memcpy(windat[n]->u1b, windat[n]->u1, window[n]->sgsiz * sizeof(double));
      wincon[n]->hor_mix->pre(window[n], windat[n], wincon[n]);
      wincon[n]->hor_mix->setup(window[n], windat[n], wincon[n]);
    }

    /*---------------------------------------------------------------*/
    /* Compute the tangential component of velocity                  */
    vel_tan_3d(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Calculate the alert diagnostics if required                   */
    if (wincon[n]->numbers & SPEED_2D || wincon[n]->means & VEL2D)
      vint_3d(window[n], windat[n], wincon[n]);
    if (!(master->alertf & NONE)) {
      vint_3d(window[n], windat[n], wincon[n]);
      alerts_w(window[n], VEL3D);
      wincon[n]->neweta = windat[n]->eta;
      alerts_w(window[n], VEL2D);
      alerts_w(window[n], ETA_A);
      alerts_w(window[n], WVEL);
      alerts_w(window[n], TRACERS);
      master_alert_fill(master, window[n], windat[n]);
    }
    windat[n]->wclk = (dp_clock() - clock);
    debug_c(window[n], D_INIT, D_POST);
  }
}

/* END transport_step()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Post transport step                                               */
/*-------------------------------------------------------------------*/
void transport_post(master_t *master, geometry_t **window,
		    window_t **windat, win_priv_t **wincon)
{
  int nwindows = master->nwindows;
  int nn, n;
  double clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Solve the tracer equation in each window                        */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];

    /* Save the value of eta for this iteration                      */
    memcpy(wincon[n]->oldeta, windat[n]->eta, 
	   window[n]->sgsizS * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Adjust the total wc mass to conserve                          */
    /*
    global_tfill(window[n], windat[n], wincon[n], NULL, -1);
    */

    /* Transfers for regions                                         */
    if (geom->nregions)
      region_transfer(master, window[n]);

    /* Set the CFL diagnostics                                         */
    if (!(master->cfl & NONE)) {
      if (windat[n]->mcfl2d < master->mcfl2d)
	master->mcfl2d = windat[n]->mcfl2d;
      if (windat[n]->mcfl3d < master->mcfl3d) {
	master->mcfl3d = windat[n]->mcfl3d;
	master->cflc = windat[n]->cflc;
      }
    }
    
    windat[n]->wclk += (dp_clock() - clock);
  }
}

/* END transport_post()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Transport streamline origin initialisation                        */
/*-------------------------------------------------------------------*/
void transport_init(master_t *master, geometry_t *window,
		    window_t *windat, win_priv_t *wincon)
{
  geometry_t *geom = master->geom;
  geometry_t *tpg;
  window_t *tpd;
  win_priv_t *tpc;
  int cc, c;

  /*-----------------------------------------------------------------*/
  /* Get the streamline origin in each window. This is only required */
  /* to be performed so that non-zero origins exist if streamline    */
  /* infromation is dumped to file.                                  */
  tpg = window->trans->tpg;
  tpd = window->trans->tpd;
  tpc = window->trans->tpc;

  /*-----------------------------------------------------------------*/
  /* Fill 3D  velocities into the window data structures.            */
  /*
  win_data_fill_3d(master, window, windat, master->nwindows);
  
  init_sigma(window, windat, wincon);
  */
  /*-----------------------------------------------------------------*/
  /* Get the surface boundary condition for vertical velocity.       */
  if(wincon->tgrid & (EXACT|SUBSET)) {
    tpd->dttr = master->dttr;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      tpd->detadt[c] = (tpd->eta[c] - tpd->etab[c]) / windat->dttr;	
    }
    vel_w_bounds(tpg, tpd, tpc);
    memcpy(tpd->etab, tpd->eta, tpg->sgsizS * sizeof(double));
  } else {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      windat->detadt[c] = (windat->eta[c] - windat->etab[c]) / 
	windat->dttr;
    }
    vel_w_bounds(window, windat, wincon);
    memcpy(windat->etab, windat->eta, 
	   window->sgsizS * sizeof(double));
  }

  /*-----------------------------------------------------------------*/
  /* Set up the cell centered surface vectors and dz arrays          */
  /* Set the dz arrays for eta at t+1 (nsur_t)                       */
  reset_dzt(window, windat, wincon);
  memcpy(wincon->s1, wincon->s2, window->sgsiz * sizeof(int));
  wincon->vc = wincon->vc2;
  wincon->vcs = wincon->vcs2;
  if(wincon->tgrid & (EXACT|SUBSET)) {
    reset_dzt(tpg, tpd, tpc);
    memcpy(tpc->s1, tpc->s2, tpg->sgsiz * sizeof(int));
    tpc->vc = tpc->vc2;
    tpc->vcs = tpc->vcs2;
  }

  /*-----------------------------------------------------------------*/
  /* Get the streamline origin                                       */
  semi_lagrange_c(window, windat, wincon);

  /* Get the boundary fluxes                                         */
  obc_flux(window, windat, wincon, 1);
}

/* END transport_init()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion for the semi-Lagrange scheme.  */
/* This routine is optimized for speed.                              */
/*-------------------------------------------------------------------*/
void advect_diffuse_lag(geometry_t *window, /* Window geometry       */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon  /* Window constants      */
			)
{
  int nn;                 /* Tracer index                            */
  int c, cc;              /* Sparse indices                          */
  int c2;                 /* 2D cell corresponding to 3D cell        */
  int cb;                 /* Bottom sparse coordinate (cell centre)  */
  int vc;                 /* Tracer cells to provess counter         */
  int vcs;                /* Surface tracer cells to process counter */
  int vc1;                /* As for vc excluding OBC cells           */
  int vcs1;               /* As for vcs excluding OBC cells          */
  int zp1;                /* k index at k+1                          */
  double *Fx;             /* x tracer flux                           */
  double *Fy;             /* y tracer flux                           */
  double *Fz;             /* z tracer flux                           */
  double *dtracer;        /* Advective terms for x,y directions      */
  double *hflux;          /* Horizontal diffusion flux               */
  double dtu;             /* Sub-time step to use                    */
  int slf = 1;            /* Set to zero on the first sub-step       */
  int hdif = 0;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* w1 = dtracer                                                    */
  /* w2 = Fx                                                         */
  /* w3 = Fy                                                         */
  /* w4 = Fz                                                         */
  /* d3 = mean cell spacing in the x direction                       */
  /* d4 = mean cell spacing in the y direction                       */
  /* s1 = wet cells to process for tracers                           */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  dtracer = wincon->w1;
  Fx = wincon->w2;
  Fy = wincon->w3;
  Fz = wincon->w4;
  hflux = wincon->w9;
  vc = wincon->vc;
  vc1 = wincon->vc1;
  vcs = wincon->vcs;
  vcs1 = wincon->vcs1;
  if (!(wincon->alertf & NONE)) {
    memcpy(windat->salb, windat->sal, window->sgsiz * sizeof(double));
    memcpy(windat->tempb, windat->temp, window->sgsiz * sizeof(double));
  }
  if (wincon->u1kh0 > 0.0) hdif = 1;

  /* Find the origin of the streamline & grid Courant numbers        */
  semi_lagrange_c(window, windat, wincon);

  if (wincon->tmode & TR_CHECK)
    check_transport(window, windat, wincon);

  /* No ultimate limiting                                            */
  wincon->ultimate = 0;
  /* Time step                                                       */
  dtu = windat->dttr;

  /*-----------------------------------------------------------------*/
  /* Set a no-gradient OBC for vertical velocity; since u & v OBCs   */
  /* are applied independently, the divergence in w boundary cells   */
  /* does not obey continuity, hence w computed in boundary cells    */
  /* are erroneous. These are, however used in the lagrange scheme.  */
  /*bdry_w(window, windat, wincon);*/

  /*-----------------------------------------------------------------*/
  /* SIGMA : Get the water depth at the forward time                 */
  if (wincon->sigma)
    for (cc = 1; cc <= vcs; cc++) {
      c = wincon->s1[cc];       /* 3D surface wet cell to process    */
      c2 = window->m2d[c];      /* 2D cell corresponding to 3D cell  */
      wincon->Hn1[c2] = windat->eta[c2] - wincon->Hs[c2];
    }

  /*-----------------------------------------------------------------*/
  /* Get the total mass in regions                                   */
  region_volume_flux_trans(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Tracer loop                                                     */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values                 */
    double *trt = wincon->w10;      /* Tracer value at time t        */

    /*---------------------------------------------------------------*/
    /* Initialise the advective terms                                */
    memset(Fx, 0, window->sgsiz * sizeof(double));
    memset(Fy, 0, window->sgsiz * sizeof(double));
    memset(dtracer, 0, window->sgsiz * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Set the bottom boundary condition (no-gradient)               */
    for (cc = 1; cc <= window->a2_t; cc++) {
      cb = window->bot_t[cc];
      tr[window->zm1[cb]] = tr[cb];
    }
    /* Ghost cells for sediments                                     */
    for (cc = 1; cc <= window->ngsed; cc++) {
      cb = window->gsed_t[cc];
      tr[cb] = tr[window->zp1[cb]];
    }
    /* Set a no-gradient above the surface (note: the surface may    */
    /* have dropped during this time-step, and the no gradient set   */
    /* in tr_bounds() at the previous time-step may not be valid.    */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = c2 = wincon->s1[cc];
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	tr[c] = tr[c2];
      }
    }

    /*---------------------------------------------------------------*/
    /* Evaluate the sources and sinks before advection               */
    if (wincon->npss) {
      ss_tracer(window, windat, wincon, n, dtracer, dtu);
      if (!(wincon->fillf & NONE))
	global_fill(window, windat, wincon, dtracer, n, 1);

      if (wincon->pssinput & WIMPLICIT) {
	double *ovol = wincon->p1;
  	/* This pss input computes the concentration assuming the    */
	/* tracer mass input is diluted by the volume input over the */
	/* timestep.                                                 */
	for (cc = 1; cc <= vc1; cc++) {
	  c = wincon->s1[cc];
	  if (dtracer[c]) {
	    double Vtm1, Vpss, Vtot;
	    double *Vi = wincon->w8;
	    c2 = window->m2d[c];
	    /* SIGMA : Adjust tracer values for the depth            */
	    if (slf)
	      tr[c] *= wincon->Ds[c2];
	    else
	      tr[c] *= wincon->Hn1[c2];
	    /* Correct formulation implemented 12.04.11. Note: ovol  */
	    /* is the volume at the previous timestep + pss volume   */
	    /* input. The required Vtm1 is the volume at t-1, i.e.   */
	    /* without pss input, hence waterss is subtracted below. */
	    Vtm1 = ovol[c];   /* Vol at t-1 + pss input              */
	    Vpss = windat->waterss[c] * dtu; /* pss input            */
	    Vtot = (Vtm1 + Vpss);
	    /* New concentration = mass at t-1 + pss mass divided    */
	    /* by total volume (i.e. volume at t-1 + volume input    */
	    /* over the timestep, approximated here as waterss).     */
	    /* Note it is possible a cell was dry at the previous    */
	    /* timestep (i.e. newly wetted this timestep), and zero  */
	    /* volume is input from the sourcesink, hence must check */
	    /* for Vtm1+Vd > 0.                                      */
	    tr[c] = (Vtot) ? (Vtm1 * tr[c] - dtracer[c]) / Vtot : tr[c];
	    tr[c] /= wincon->Hn1[c2];
	  }
	}
      } else {
	/* This pss input injects the tracer mass directly into the  */
	/* cell.                                                     */
	for (cc = 1; cc <= vc1; cc++) {
	  c = wincon->s1[cc];
	  if (dtracer[c]) {
	    c2 = window->m2d[c];
	    /* SIGMA : Adjust tracer values for the depth            */
	    if (slf)
	      tr[c] *= wincon->Ds[c2];
	    else
	      tr[c] *= wincon->Hn1[c2];
	    tr[c] -= (dtracer[c] / (window->cellarea[c2]) * wincon->dz[c]);
	    tr[c] /= wincon->Hn1[c2];
	  }
	}
      }
    }
    memcpy(trt, tr, window->sgsiz * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Set the bottom boundary condition (no-gradient)               */
    for (cc = 1; cc <= window->a2_t; cc++) {
      cb = window->bot_t[cc];
      trt[window->zm1[cb]] = trt[cb];
    }
    /* Ghost cells for sediments                                     */
    for (cc = 1; cc <= window->ngsed; cc++) {
      cb = window->gsed_t[cc];
      trt[cb] = trt[window->zp1[cb]];
    }

    /*---------------------------------------------------------------*/
    /* Get the mass fluxes in the regions                            */
    region_flux_trans(window, windat, wincon, n, dtracer);

    /*---------------------------------------------------------------*/
    /* Get the updated concentration                                 */
    if (wincon->advect[n])
      semi_lagrange(window, windat, wincon, tr);

    /*---------------------------------------------------------------*/
    /* Save the advective fluxes if required                         */
    if (wincon->trflux == n)
      calc_flux(window, windat, wincon, Fx, Fz, dtu, 1, 4);

    /*---------------------------------------------------------------*/
    /* Get the horizontal diffusive fluxes                           */
    if (wincon->diffuse[n] && hdif) {
      hor_diffuse(window, windat, wincon, tr, Fx);
      memset(hflux, 0, window->sgsiz * sizeof(double));
      for (cc = 1; cc <= window->v3_t; cc++) {
	int e, j;
	c = window->w3_t[cc];
	c2 = window->m2d[c];
	for (j = 1; j <= window->npe[c2]; j++) {
	  e = window->c2e[j][c];
	  hflux[c] += (window->eSc[j][c2] * Fx[e] * dtu);
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Evaluate the horizontal diffusion                             */
    if (hdif) {
      for (cc = 1; cc <= vc1; cc++) {
	c = wincon->s1[cc];     /* Wet cell to process               */
	c2 = window->m2d[c];    /* 2D cell corresponding to 3D cell  */
	/* SIGMA : Adjust tracer values for the depth                */
	if (slf)
	  tr[c] *= wincon->Ds[c2];
	else
	  tr[c] *= wincon->Hn1[c2];
	tr[c] -= (hflux[c] / (window->cellarea[c2]) * wincon->dz[c]);
	tr[c] /= wincon->Hn1[c2];
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the total mass prior to filling                           */
    region_mass_tr(window, windat, wincon, -1.0, n, RG_GLOB);

    /*---------------------------------------------------------------*/
    /* Apply a multiplicative global fill to conserve mass           */
    if (!(wincon->fillf & NONE)) {
      int mode = (nn == 0) ? 2 : 3;
      mode = (nn == wincon->ntbdy - 1) ? 4 : mode;
      mode = (wincon->ntbdy == 1) ? 6 : mode;
      global_fill(window, windat, wincon, dtracer, n, mode);
    }
  }                         /* Tracer loop end                       */

  /* Clip the tracer concentration if required                       */
  tr_bounds(window, windat, wincon);
  check_bounds("transport", window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Reset the new surface coordinate at thin layers if required     */
  if (wincon->thin_merge) {
    for (vc = 1; vc <= wincon->nkth_e1; vc++) {
      cc = (int)wincon->kth_e1[vc];
      c = window->nsur_t[cc];
      window->nsur_t[cc] = window->zp1[c];
    }
  }
}

/* END advect_diffuse_lag()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the grid     */
/* Courant numbers.                                                  */
/*-------------------------------------------------------------------*/
void semi_lagrange_c(geometry_t *window,  /* Window geometry         */
                     window_t *windat,    /* Window data             */
                     win_priv_t *wincon   /* Window constants        */
  )
{
  double dt;                    /* Sub-time step                     */
  double time_left;             /* Number of sub-timesteps per dt    */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double *nu, *nv;              /* Cell centered horizontal velocity */
  double *nw = wincon->w10;     /* Cell centered vertical velocity   */
  double *cx = wincon->w6;      /* Streamline x location             */
  double *cy = wincon->w7;      /* Streamline y location             */
  double *cz = wincon->w8;      /* Streamline z location             */
  double *dzz = wincon->dzz;    /* Cell thickness                    */
  double *mask = wincon->w9;    /* Mask for sub-step scaling         */
  int *cl = wincon->s2;         /* Cell index of source location     */
  int *c2cc = wincon->s4;       /* Index to counter mapping          */
  double u, v, w;               /* Interpolated velocity             */
  double p, q, r;               /* Sub-steps                         */
  int cc, c, ci, c2, kb, zm1;
  int n;
  double slon, slat, dist, dir, bot;
  int i, j, k;
  double sinth, costh;
  double m2deg = 1.0 / (60.0 * 1852.0);
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int onestep = 0;  /* Find streamline in one step                   */
  if (!is_geog) m2deg = 1.0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(cx, 0, window->szc * sizeof(double));
  memset(cy, 0, window->szc * sizeof(double));
  memcpy(cl, wincon->s1, window->szc * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing                                   */
  set_dzz(window, dzz);

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face / center          */
  nu = windat->u;
  nv = windat->v;
  vel_center_w(window, windat, wincon, nw);

  /* Set the ghost cells                                             */
  for (cc = 1; cc <= window->nbpt; cc++) {
    int e, es;
    c = window->bpt[cc];
    ci = window->bin[cc];
    j = window->dbpt[cc];
    e = window->c2e[j][c];
    es = window->m2de[e];

    /* Rotate the edge velocity to east and north, and use this at   */
    /* at the ghost cell.                                            */
    nu[c] = windat->u1[e] * window->costhu1[es] + windat->u2[e] * window->costhu2[es];
    nv[c] = windat->u1[e] * window->sinthu1[es] + windat->u2[e] * window->sinthu2[es];
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value                                          */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    zm1 = window->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
  }
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    ci = window->zp1[c];
    nu[c] = nu[ci];
    nv[c] = nv[ci];
    nw[c] = nw[ci];
  }
  /* Set the open boundary ghost cells                               */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    int ee;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      c2 = open->obc_e2[ee];
      do {
	nu[c] = nu[c2];
	nv[c] = nv[c2];
	nw[c] = nw[c2];
	c = open->omape[ee][c];
      } while (c != open->omape[ee][c]);
    }
    /*OBC_bgz_nogradb(window, open, tr);*/
  }

  /* Set the cell width  mask                                        */
  for (c = 1; c <= window->enon; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bin[cc];
    mask[c] = 0.2;
  }

  /*-----------------------------------------------------------------*/
  /* Set the initial streamline location                             */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    cx[c] = window->cellx[c2];
    cy[c] = window->celly[c2];
  }
  memcpy(cz, wincon->cellz, window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Initialize the interpolation structure                          */
  delaunay_reinit(window, windat->d, 0, nu);
  delaunay_reinit(window, windat->d, 1, nv);
  delaunay_reinit(window, windat->d, 2, nw);

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin                         */
  for (cc = 1; cc <= wincon->vc; cc++) {
    int first = 0;
    c = cg = wincon->s1[cc];
    time_left = windat->dttr;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    n=0;
    /*
    if (window->s2i[c] == 22 && window->s2j[c] == 13 && window->s2k[c] == 22) {
      int ii =  window->s2i[cl[cc]];
      int jj =  window->s2j[cl[cc]];
      int kk =  window->s2k[cl[cc]];
      int id = window->c2p[window->s2k[c]][window->m2d[c]];
      printf("%d %d\n",c,id);
    }
    */
    while (time_left > 0) {

      int ins;
      double ua = max(fabs(u), SMALL);
      double va = max(fabs(v), SMALL);
      double wa = max(fabs(w), SMALL);

      /* Get the sub-timestep for the next integration. Streamlines  */
      /* cannot cross more than one cell in one sub-timestep.        */
      ci = cl[cc];
      ins = (ci == wincon->m2d[ci]) ? 1 : 0;
      c2 = (window->wgst[ci]) ? window->m2d[window->wgst[ci]] : window->m2d[ci];

      p = q = sqrt(window->cellarea[c2]) * mask[c];
      ua = sqrt(ua * ua + va * va);
      r = (w < 0.0) ? dzz[ci] : dzz[window->zm1[ci]];
      dt = min(SCALE * fabs(p / ua),
               SCALE * fabs(q / va));

      dt = (r) ? ((ins && nw[c] < 0.0) ? dt : min(dt, SCALE * fabs(r / wa))) : dt;
      dt = (onestep) ? time_left : min(dt, time_left);

      /* Get the velocties at the origin                             */
      if (!onestep && first) {
	u = hd_trans_interp(window, windat->gsx, cx[c], cy[c], cz[c], cl[cc],c,0);
	v = hd_trans_interp(window, windat->gsy, cx[c], cy[c], cz[c], cl[cc],c,1);
	w = hd_trans_interp(window, windat->gsz, cx[c], cy[c], cz[c], cl[cc],c,2);
      }
      first = 1;

      /* Get the new location of the streamline origin               */
      cl[cc] = get_pos(window, c, cl[cc], u, v, w, &cx[c], &cy[c], &cz[c], dt, 0);

      /* Get the number of sub-timesteps used                        */
      time_left = time_left - dt;
      n++;
    }

    /* Copy the streamline origin and offsets into t11, t12 and t22  */
    /* respectively, so that these quantities are available for      */
    /* dumping to sparse output files. Note that these variables may */
    /* only be dumped to sparse file since interpolation from        */
    /* (x,y,z) format is not possible for the sparse coordinate of   */
    /* the streamline origin.                                        */
    /*
    windat->pc[c] = cx[c];
    windat->qc[c] = cy[c];
    windat->rc[c] = cz[c];
    */
  }
  print_source_k(window, wincon->vc, wincon->s1, cx, cy, 22, window->s2k);
  debug_c(window, D_TS, D_STRML);
}

/* END semi_lagrange_c()                                             */
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
int get_pos(geometry_t *window, /* Window geometry                   */
	    int c,              /* Location of destination           */
	    int ci,             /* Current source cell               */
	    double u,           /* East velocity                     */
	    double v,           /* North velocity                    */
	    double w,           /* Vertical velocity                 */
            double *cx,         /* x location of streamline          */
	    double *cy,         /* y location of streamline          */
	    double *cz,         /* z location of streamline          */
	    double dt,          /* Time step                         */
	    int mode            /* mode=0, centres, 1=layer faces    */
	    )
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  delaunay *d = window->d;
  int cs = window->m2d[c];
  int cis = window->m2d[ci];
  int cn, cns, zm1;
  double xin = *cx;
  double yin = *cy;
  double zin = *cz;
  double dist, dir, s1, s2;
  double slat, slon;
  double sinth, costh;
  double dx, dy;
  double nx, ny;
  double eps = 1e-8;
  double m2deg = 1.0 / (60.0 * 1852.0);
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int j, k = window->s2k[ci];
  point p;

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
  if (window->us_type & US_IJ) {
    /* Structured grids: use the xytoi tree                          */
    cn = hd_grid_xyztoc(window, slon, slat, zin, ci);
    cns = window->m2d[cn];
  } else {
    /*
    p.x = slon;
    p.y = slat;
    d->first_id = window->c2tri[cis];
    cns = delaunay_xytoi_lag(d, &p, d->first_id);
    cns = geom->tri2c[cns];
    */
    /* Unstructured meshes: walk through the Voronoi mesh            */
    cn = find_cell(window, c, slon, slat, &nx, &ny);
    if (!cn) cn = c;
    cns = window->m2d[cn];
  }

  /* Get the vertical layer of the source cell                       */
  if (cn == window->zm1[cn]) cn = window->zp1[cn];
  if (window->wgst[cn]) cn = window->wgst[cn];
  k = window->s2k[cn];
  if (k > window->nz - 1) 
    hd_quit("Can't find streamline position: destination = %d[%f %f]\n",
	    c, window->cellx[cs], window->celly[cs]);

  /* If the new location is a ghost cell, send it inside the mesh.   */
  if (window->gcm[k][cns] & L_GHOST) {
    if (window->us_type & US_IJ) {
      int es, vs, found = 0;
      for (j = 1; j <= window->npe[cis]; j++) {
	if (intersect(window, cis, j, xin, yin, slon, slat, &nx, &ny))
	  break;
      }
    }

    dx = nx - xin;
    dy = ny - yin;
    dist = sqrt(dx * dx + dy * dy) - eps;
    slon = xin - dist * costh;
    slat = yin - dist * sinth;
    if (window->us_type & US_IJ) {
      cn = hd_grid_xyztoc(window, slon, slat, zin, ci);
      cns = window->m2d[cn];
    } else {
      /*
      p.x = slon;
      p.y = slat;
      d->first_id = window->c2tri[cis];
      cns = delaunay_xytoi_lag(d, &p, d->first_id);
      cns = geom->tri2c[cns];
      */
      cns = find_cell(window, c, slon, slat, &nx, &ny);
      if (!cns) cns = c;
      cns = window->m2d[cns];
    }
  }
  *cx = slon;
  *cy = slat;
  *cz -= w * dt;

  /* Find the vertical position in the mesh                          */
  cn = cns;
  /* First get the layer of the free surface                         */
  while (windat->eta[cns] < window->gridz[cn] && c != window->zm1[cn]) {
    cn = window->zm1[cn];
  }
  /* Get the first cellz layer less than z                           */
  if (mode) {
    while (*cz <= window->gridz[cn] && cn != window->zm1[cn]) {
      cn = window->zm1[cn];
    }
  } else {
    while (*cz <= wincon->cellz[cn] && cn != window->zm1[cn]) {
      cn = window->zm1[cn];
    }
  }
  return(cn);
}

/* END get_pos()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the new location of the streamline origin such that the  */
/* streamline traverses only one cell, i.e. the source location lies */
/* on an edge. Returns the cell centre of the source location.       */
/* If the streamline lands in a ghost cell, then shift into the wet  */
/* interior and re-compute a new location.                           */
/* Destination cell is the origin of the streamline.                 */
/* Source cell is the end point of the streamline track.             */
/* Returns a new source cell after moving a distance increment.      */
/* The new geographic location is returned in cx and cy.             */
/*-------------------------------------------------------------------*/
int get_posc(geometry_t *window, /* Window geometry                  */
	     int *ei,            /* Source edge                      */
	     int ci,             /* Current source cell              */
	     double u,           /* East velocity                    */
	     double v,           /* North velocity                   */
	     double *cx,         /* x location of streamline         */
	     double *cy,         /* y location of streamline         */
	     double dt,          /* Time step                        */
	     int mode,           /* mode=0, centres, 1=layer faces   */
	     double *odist       /* Distance of this segment         */
	    )
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int cis = window->m2d[ci];
  int cn, cns, zm1;
  double xin = *cx;
  double yin = *cy;
  double dist, dir, s1, s2;
  double slat, slon;
  double sinth, costh;
  double dx, dy;
  double nx, ny;
  double m2deg = 1.0 / (60.0 * 1852.0);
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int sodanos = 0;  /* Compute distances on the spheroid             */
  int j, k = window->s2k[ci], e;
  double xi, yi;
  int found;

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

  /* Return current cell if the streamline is at the end             */
  if(dist == 0.0) {
    *odist = 0.0;
    return(-ci);
  }

  /* Find the first edge the streamline crosses                      */
  cn = ci;
  for (j = 1; j <= window->npe[cis]; j++) {
    e = window->c2e[j][ci];
    if (e != *ei && (found = intersect(window, ci, j, *cx, *cy, slon, slat, &xi, &yi))) {
      cn = window->c2c[j][ci];
      break;
    }
  }
  *ei = e;

  /* Send back if the next cell is a ghost cell */
  k = window->s2k[cn];
  cns = window->m2d[cn];
  if (window->gcm[k][cns] & L_GHOST) {
    *odist = get_dist(*cx, *cy, xi, yi) / m2deg;
    *cx = xi;
    *cy = yi;
    return(-ci);
  }
  /* Find the first edge the streamline crosses                      */
  /* If the streamline ends in a cell (i.e. doesn't cross an edge)   */
  /* then return its source location, otherwise return the location  */
  /* of the edge it crosses.                                         */
  if (!found) {
    *odist = get_dist(*cx, *cy, slon, slat) / m2deg;
    *cx = slon;
    *cy = slat;
    return(-ci);
  } else {
    *odist = get_dist(*cx, *cy, xi, yi) / m2deg;
    *cx = xi;
    *cy = yi;
    return(cn);
  }
}

/* END get_posc()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed in the tracer loop to get the updated      */
/* concentration given the grid Courant numbers.                     */
/*-------------------------------------------------------------------*/
void semi_lagrange(geometry_t *window,  /* Processing window */
		   window_t *windat,  /* Window data structure */
		   win_priv_t *wincon,  /* Window geometry / constants */
		   double *tr   /* Tracer array */
		   )
{
  int cc, c, cs, c1;
  double *ntr = wincon->w5;
  int *cl = wincon->s2;
  int *c2cc = wincon->s4;
  double *cx = wincon->w6;
  double *cy = wincon->w7;
  double *cz = wincon->w8;
  double *botz = wincon->d6;

  /* Update the tracer values at cell centres in the interpolation   */
  /* structure.                                                      */
  delaunay_reinit(window, windat->d, 3, tr);

  memcpy(ntr, tr, window->szc * sizeof(double));

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    ntr[c] = hd_trans_interp(window, windat->gst, cx[c], cy[c], cz[c], cl[cc],c,3);
  }

  memcpy(tr, ntr, window->szc * sizeof(double));
}

/* END semi_lagrange()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the cell a point resides in for a Voronoi mesh. Uses the    */
/* property of the Voronoi diagram that all points in a Voronoi cell */
/* have a minimum distance to that cell centre. Here we walk through */
/* the mesh following the path between destination (origin) and      */
/* source cells as the crow flies, looking for the minimum distance  */
/* from the source location to the cell centre.   .                  */
/* If the cell centre of the walk arrives in a ghost cell and the    */
/* source-destination distance is less than the                      */
/* intersection-destination distance, then return the wet cell       */
/* interior to the ghost cell. If it's greater, then return the      */
/* ghost cell.                                                       */
/*-------------------------------------------------------------------*/
int find_cell(geometry_t *window, /* Window geometry                 */
	      int cd,             /* Destination (origin) cell       */
	      double xs,          /* x location of source cell       */
	      double ys,          /* y location of source cell       */
	      double *xi,    /* x location of intersection with edge */
	      double *yi     /* y location of intersection with edge */
	      )
{
  int cs = window->m2d[cd];      /* Surface cell                     */
  double xd = window->cellx[cs]; /* x location of destination        */
  double yd = window->celly[cs]; /* y location of destination        */
  double xn, yn;                 /* New cell centre location         */
  double dmx, dmn, d, di;
  int j, cn, cr, cg;

  cn = cr = cd; /* Start walk from the current destination cell      */
  xn = xd;   
  yn = yd;
  /* Note: as we're only interested in relative distances, we do not */
  /* need to convert distance to m.                                  */
  dmx = dmn = d = get_dist(xs, ys, xn, yn);
  while (d < dmn) {
    int found = 0;
    for (j = 1; j <= window->npe[cs]; j++) {
      if ((found = intersect(window, cn, j, xn, yn, xs, ys, xi, yi))) {
	cn = window->c2c[j][cn];
	break;
      }
    }
    if (!found) return(0);
    cs = window->m2d[cn];
    xn = window->cellx[cs];
    yn = window->celly[cs];
    di = get_dist(xd, yd, *xi, *yi);
    if ((d = get_dist(xs, ys, xn, yn)) < dmn) {
      cr = cn;
      dmn = d;
    }
    if ((cg = window->wgst[cn])) {
      if (dmx < di)
	return(cg);
      else
	return(cn);
    }
  }
  return(cr);
}

/* END find_cell()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the distance between two points                             */
/*-------------------------------------------------------------------*/
double get_dist(double x1, double y1, double x2, double y2) {
  double dx = x2 - x1;
  double dy = y2 - y1;
  return(sqrt(dx * dx + dy * dy));
}

/* END get_dist()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 if a line from cell centre c to (xs, ys) intersects     */
/* edge j. If so, the location of the intersection are returned in   */
/* (xi, yi).                                                         */
/*-------------------------------------------------------------------*/
int intersect(geometry_t *window,    /* Window geometry              */
	      int c,                 /* Cell containing destination  */
	      int j,                 /* Edge to check                */
	      double xn,             /* x location of destination    */
	      double yn,             /* y location of destination    */
	      double xs,             /* x location of source         */
	      double ys,             /* y location of source         */
	      double *xi,            /* x location of intersection   */
	      double *yi             /* y location of intersection   */
	      )
{
  int cs = window->m2d[c];
  int es = window->c2e[j][cs];   /* Edge of surface source cell      */
  double x1, y1, x2, y2, s1, s2;
  int found = 0;
  if (es) {
    int v1 = window->e2v[es][0];       /* Vertex 1 of edge           */
    int v2 = window->e2v[es][1];       /* Vertex 2 of edge           */
    if (v1 && v2) {
      double v1x = window->gridx[v1];
      double v1y = window->gridy[v1];
      double v2x = window->gridx[v2];
      double v2y = window->gridy[v2];
      double x1, x2, y1, y2;
      double dx = fabs(v1x - v2x);
      double dy = fabs(v1y - v2y);

      /* Get the location where edge and streamline intersect       */
      if (dx) {
	if (xs == xn) {
	  *xi = xn;
	  *yi = yn;
	  return(1);
	}
	s1 = (v2y - v1y) / (v2x - v1x);
	s2 = (ys - yn) / (xs - xn);
	*xi = (s2 * xn - s1 * v1x + v1y - yn) / (s2 - s1);
	*yi = s1 * (*xi - v1x) + v1y;
      } else {
	s2 = (ys - yn) / (xs - xn);
	*xi = v1x;
	*yi = s2 * (v1x - xn) + yn;
      }
      /* If the intersection lies between the bounds of the vertices */
      /* then the streamline crosses this edge.                      */
      if (*xi >= min(v1x, v2x) &&
	  *xi <= max(v1x, v2x) &&
	  *yi >= min(v1y, v2y) &&
	  *yi <= max(v1y, v2y)) {
	if (*xi >= min(xn, xs) &&
	    *xi <= max(xn, xs) &&
	    *yi >= min(yn, ys) &&
	    *yi <= max(yn, ys)) {
	  return(1);
	}
      }
      /*
      x1 = fabs(*xi - v1x);
      x2 = fabs(*xi - v2x);
      y1 = fabs(*yi - v1y);
      y2 = fabs(*yi - v2y);
      if (x1 <= dx && x2 <= dx && y1 <= dy && y2 <= dy) {
	return(1);
      }
      */
    }
  }
  return(found);
}

/* END intersect()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the cell a point resides in for a Voronoi mesh. Uses the    */
/* property of the Voronoi diagram that all points in a Voronoi      */
/* have a minimum distance to that cell centre. Here we search the   */
/* cells surrounding the destination (origin) cell  and find the     */
/* minimum distance to the source cell. This assumes the source cell */
/* is not more than one cell distant from the destination, which in  */
/* most cases will be true as semi_lagrange_c() sub-steps to ensure  */
/* this. However, if neighbour cell with is much less than the       */
/* destination cell, then it is possible this assumption does not    */
/* hold.                                                             */
/*-------------------------------------------------------------------*/
int find_cellq(geometry_t *window, /* Window geometry                */
	       int c,              /* Destination (origin) cell      */
	       double x,           /* x location of source cell      */
	       double y            /* y location of source cell      */
	       )
{
  int cr, cs = window->m2d[c];
  double dx = window->cellx[cs] - x;
  double dy = window->celly[cs] - y;
  double dist = sqrt(dx * dx + dy * dy);
  double r;
  int j, cc, cn, sz;

  cr = c;
  r = dist;
  for (j = 1; j <= window->npe[cs]; j++) {
    cn = window->c2c[j][cs];
    dx = window->cellx[cs] - x;
    dy = window->celly[cs] - y;
    dist = sqrt(dx * dx + dy * dy);
    if (dist < r) {
      cr = cn;
      r = dist;
    }
  }
  return(cr);
}

/* END find_cellq()                                                  */
/*-------------------------------------------------------------------*/

#define T_NONE    1
#define T_VERT    2
#define T_HORZ    4
#define T_HV      8

/*-------------------------------------------------------------------*/
/* Flux Form semi-Lagrange: tracks streamlines from cell edges and   */
/* piecewise linearly integrates tracer values to get a streamline   */
/* mean value, which it then applies to that edge. Also invokes a    */
/* limiter (using streamline mean cell centered values) to keep the  */
/* solution monotonic.                                               */
/*-------------------------------------------------------------------*/
void ffsl_don(geometry_t *window,           /* Window geometry        */
	     window_t *windat,             /* Window data            */
	     win_priv_t *wincon,           /* Window constants       */
	     double *tr,                   /* Tracer array           */
	     double *Fx,                   /* Horizontal flux        */
	     double *Fz,                   /* Vertical flux          */
	     double dtu                    /* Timestep               */
	  )
{
  double dt;                  /* Sub-time step for Euler velocity    */
  double time_left;           /* Number of sub-timesteps per dt      */
  double *u1;                 /* Horizontal edge velocity            */
  double *w;                  /* Vertical velocity                   */
  double *u1v;                /* Horizontal volume flux              */
  double *nu, *nv;            /* Cell centered horizontal velocity   */
  double *nw = wincon->nw;    /* Cell centered vertical velocity     */
  double *cx = wincon->tr_mod_x;  /* x position of streamline        */
  double *cy = wincon->tr_mod_y;  /* y position of streamline        */
  double *cz = wincon->tr_mod_z;  /* z poition of streamline         */
  double *tr_mod = wincon->tr_mod; /* z transverse tracer            */
  double *crfzf;   /* Fractional factors of face trajectories        */
  double *crfzc;   /* Fractional factors of centre trajectories      */
  int *clzf;       /* Cell counter at face source cell               */
  int *clzc;       /* Cell counter cell centre source                */
  double *ntr = wincon->w5;   /* New tracer value                    */
  double *dint = wincon->w3;  /* Streamline integrated tracer value  */
  double *smin = wincon->crfxc; /* Minimum streamline value          */
  double *smax = wincon->crfyc; /* Maximum streamline value          */
  double *sminz = wincon->crfxf; /* Minimum vertical streamline value*/
  double *smaxz = wincon->crfyf; /* Maximum vertical streamline value*/
  double *tr_moc = wincon->Fzh;
  int *cl = wincon->s2;       /* Cell location of streamline         */
  int *e2ee = wincon->s4;     /* Mapping from edge to index          */
  int *mask = wincon->s3;
  double u, v;                /* Velocities                          */
  double trs, tre;            /* Tracer at start and end of segments */
  double dz;                  /* Cell thickness                      */
  double dist, d1;            /* Lengths of segment                  */
  double sinth, costh;        /* sin and cos of edge angle           */
  double m2deg = 1.0 / (60.0 * 1852.0);  /* Meter to degree          */
  double px, py, pz;
  int cc, c, cs, ci, kb, zm1, zp1;
  int ee, e, e2, es, ei, pei;
  int c1, c2;
  int n;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  double tf = 0.55; /* Factor for implicitness of cross terms (0.5:1)*/
  int tmode = (wincon->means & TRANSPORT) ? 1 : 0;
  int sodanos = 0;            /* Compute distances on spheriod       */
  int tranf;                  /* Transverse term flag                */
  int doint = 1;              /* Integrate tracers along streamline  */
  int momo2 = 0;              /* Second order momentum estimation    */

  /* Set which transverse terms are used                             */
  tranf = (T_VERT|T_HORZ|T_HV);
  if (!is_geog) m2deg = 1.0;
  if (tmode) {
    u1 = windat->ume;
    w = windat->wm;
    u1v = windat->u1vm;
  } else {
    u1 = windat->u1;
    w = windat->w;
    u1v = windat->u1flux3d;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(cx, 0, window->sze * sizeof(double));
  memset(cy, 0, window->sze * sizeof(double));
  memset(cz, 0, window->sze * sizeof(double));
  memset(ntr, 0.0, window->sze * sizeof(double));
  memset(dint, 0.0, window->sze * sizeof(double));
  memset(tr_moc, 0.0, window->szc * sizeof(double));
  memset(mask, 0, window->szc * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Define a new tracer with z-direction transverse terms           */
  clzf = wincon->clzf;
  clzc = wincon->clzc;
  crfzf = wincon->crfzf;
  crfzc = wincon->crfzc;
  memcpy(tr_mod, tr, window->szc * sizeof(double));
  if (tranf & T_VERT) {
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      zp1 =  (nw[c] < 0.0) ? clzc[c] : window->zm1[clzc[c]];
      ci = (nw[c] < 0.0) ?  window->zm1[zp1] : clzc[c];
      tr_mod[c] = tf * tr[ci] + (1.0 - tf) * tr[c];
      /* Note: tf determines how much of the tracer at the forward   */
      /* timestep is used (i.e. the 'implicitness' of the solution). */
      /* This should lie between 0.5 (centered) or 1 (implicit). The */
      /* implicit solutions tend to be more diffuse.                 */
      if (crfzc[c] > 0.) {
	tr_mod[c] += tf * (crfzc[c] * (tr[zp1] - tr[ci]));
      }
    }
    set_tr_nograd(window, tr_mod);
    /* OBC ghost cells                                               */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int trn = (int)tr[0];
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = c2 = open->obc_e2[ee];
	zp1 =  (nw[c] < 0.0) ? clzc[c] : window->zm1[clzc[c]];
	ci = (nw[c] < 0.0) ?  window->zm1[zp1] : clzc[c];
	do {
	  c = open->omape[ee][c];
	  zp1 = open->omape[ee][zp1];
	  ci = open->omape[ee][ci];
	  tr_mod[c] = tf * tr[ci] + (1.0 - tf) * tr[c];
	  if (crfzc[c2] > 0.)
	    tr_mod[c] += tf * (crfzc[c2] * (tr[zp1] - tr[ci]));
	} while (c != open->omape[ee][c]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the initial streamline location                             */
  for (ee = 1; ee <= window->a3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    cx[e] = window->u1x[e2];
    cy[e] = window->u1y[e2];
    cz[e] = 0.5 * (wincon->cellz[c1] + wincon->cellz[c2]);
    e2ee[e] = ee;
  }
  /*
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    smin[e] = HUGE;
    smax[e] = -HUGE;
  }
  */
  /* Get the initial cell centre; this is the cell upstream of the   */
  /* edge.                                                           */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    for (n = 1; n <= window->npe[c2]; n++) {
      e = window->c2e[n][c];
      ee = e2ee[e];
      if (u1[e] > 0.0)
	cl[ee] = (window->eSc[n][c2] < 0) ? window->c2c[n][c] : c;
      else
	cl[ee] = (window->eSc[n][c2] < 0) ? c : window->c2c[n][c];
    }
  }
  /* Use interior cells for normal open boundary edges               */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    int trn = (int)tr[0];
    if (!(open->bcond_tra[trn] & TRCONC)) {
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
       cl[e2ee[e]] = open->obc_e2[ee];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Limit least squares interpolation functions                     */
  get_linear_limit(window, windat, wincon, tr_mod);

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin                         */
  for (ee = 1; ee <= window->a3_e1; ee++) {
    int first = 0;
    e = window->w3_e1[ee];
    if (window->eask[e] & (W_NOBC|W_TOBC)) continue;
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    time_left = dtu;
    u = windat->u1[e] * window->costhu1[es] + windat->u2[e] * window->costhu2[es];
    v = windat->u1[e] * window->sinthu1[es] + windat->u2[e] * window->sinthu2[es];
    n = 0;
    ei = e;
 
    if (wincon->fillf & CLIP) smin[e] = smax[e] = tr[cl[ee]];
    if (wincon->ultimate) smin[e] = smax[e] = tr_mod[cl[ee]];

    if (tmode) first = 1;

    while (cl[ee] > 0) {

      /* Set the current cell the streamline resides in. Note; these */
      /* locations are cell centre locations, stored for every edge. */
      c = ci = cl[ee];
      pz = wincon->cellz[ci];

      /* Get the velocties at the origin. For the first segment      */
      /* sub-step use (rotated u1) velocities at the edge location.  */
      if (first) {
	u = hd_trans_interp(window, windat->gsx, cx[e], cy[e], pz, c, 0, 0);
	v = hd_trans_interp(window, windat->gsy, cx[e], cy[e], pz, c, 0, 1);
      }
      first = 1;
      /* Assume only normal flow through normal OBC edges            */
      if (window->eask[e] & W_NOBC) {
	u = windat->u1[e] * window->costhu1[es];
	v = windat->u1[e] * window->sinthu1[es];
      }

      /* Get the new location of the streamline origin               */
      px = cx[e];
      py = cy[e];
      pei = ei;
      cl[ee] = get_posc(window, &ei, c, u, v, &cx[e], &cy[e], time_left, 1, &dist);

      /* Use a second order approximation for the velocity used to   */
      /* track the streamline.                                       */
      if (momo2) {
	double u2 = hd_trans_interp(window, windat->gsx, cx[e], cy[e], pz, c, 0, 0);
	double v2 = hd_trans_interp(window, windat->gsy, cx[e], cy[e], pz, c, 0, 1);
	u = 0.5 * (u + u2);
	v = 0.5 * (v + v2);
	cx[e] = px;
	cy[e] = py;
	ei = pei;
	cl[ee] = get_posc(window, &ei, c, u, v, &cx[e], &cy[e], time_left, 1, &dist);
      }

      /* If the source lies in a wet cell, then get the time it's    */
      /* taken to cross the cell.                                    */
      if (cl[ee] > 0) {
	d1 = sqrt(u * u + v * v);
	dt = dist / d1;
	time_left -= dt;
      }
      /*
      if (window->wgst[cl[ee]] && window->zp1[cl[ee]] != window->wgst[cl[ee]] &&
	  ee <= window->v3_e1) 
	printf("Streamline in ghost @ %f %d : %f %f\n",windat->days, e, window->cellx[window->m2d[cl[ee]]], 
				       window->celly[window->m2d[cl[ee]]]);

      if (cl[ee] <=0 || cl[ee] >= window->szc) 
	printf("Streamline out of bounds @ %f %d %f %f\n", windat->days, e, window->cellx[window->m2d[cl[ee]]], 
				       window->celly[window->m2d[cl[ee]]]);
      */
      trs = get_linear_value(window, windat, wincon, tr_mod, ci, px, py, pz);
      tre = get_linear_value(window, windat, wincon, tr_mod, ci, cx[e], cy[e], pz);

      /* Get the minimum and maximum values encountered in the       */
      /* stencil along the streamline.                               */
      if (wincon->ultimate) {
	get_local_bounds(window, windat, wincon, tr_mod[cl[cc]], NULL, e, abs(cl[ee]), smin, smax, 0);
	get_local_bounds(window, windat, wincon, tr_mod[cl[cc]], NULL, e, abs(cl[ee]), smin, smax, 0);
      } else if (wincon->fillf & CLIP) {
	get_local_bounds(window, windat, wincon, tre, tr, e, abs(cl[ee]), smin, smax, 2);
	get_local_bounds(window, windat, wincon, tre, tr, e, abs(cl[ee]), smin, smax, 2);
      }

      /* Integrate (linearly) the tracer along the streamline        */
      /* segment.                                                    */
      dint[e] += dist;

      if (doint) {
	ntr[e] += 0.5 * dist * (trs + tre);
      } else
	ntr[e] += (tr_mod[ci] * dist);

      /* Get the number of sub-timesteps used                        */
      n++;
      if(n > 100) {
	printf("stuck in loop e=%d [%f %f]\n",e,window->u1x[window->m2de[e]],window->u1y[window->m2de[e]]);
	exit(0);
      }
      /* Get the contributions to the horizontal transverse terms    */
      /* for the vertical advection, this is the tracer value at the */
      /* current source location.                                    */
      if (tranf & T_HV) {
	if (cl[ee] < 0) {
	  tre = get_linear_value(window, windat, wincon, tr, ci, cx[e], cy[e], pz);
	  tr_moc[c1] += tre;
	  tr_moc[c2] += tre;
	  mask[c1]++;
	  mask[c2]++;
	}
      }
    }
    /* Get the mean tracer value along the streamline                */
    Fx[e] = (dint[e]) ? ntr[e] / dint[e] : Fx[e];

    /* Universal flux limiter Eq 34 & 35 of Thuburn (1995)           */
    if (wincon->ultimate) {
      Fx[e] = min(Fx[e], smax[e]);
      Fx[e] = max(Fx[e], smin[e]);
    }
  }

  debug_c(window, D_TS, D_STRML);

  /* Multiply the tracer by the volume flux                          */
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    Fx[e] *= u1v[e];
  }

  /* Get the horizontal transverse terms by interpolating tr at the  */
  /* source location. Note: the weights and slope limiter for the    */
  /* least squares interpolation must be reset using tr (rather than */
  /* tr_mod).                                                        */
  if (!(tranf & T_HV)) {
    get_linear_limit(window, windat, wincon, tr);
    for (ee = 1; ee <= window->a3_e1; ee++) {
      e = window->w3_e1[ee];
      c = -cl[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      pz = wincon->cellz[c];
      tre = get_linear_value(window, windat, wincon, tr, c, cx[e], cy[e], pz);
      tr_moc[c1] += tre;
      tr_moc[c2] += tre;
      mask[c1]++;
      mask[c2]++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Define a new tracer with horizontal-direction transverse terms  */
  /* (note: we use the interpolated value at the streamline end      */
  /* directly in tr_moc).                                            */
  if (tranf & T_HORZ) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      tr_moc[c] /= (double)mask[c];
      tr_moc[c] = tf * tr_moc[c] + (1.0 - tf) * tr[c];
    }
    set_tr_nograd(window, tr_moc);
  } else
    memcpy(tr_moc, tr, window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes                                                 */
  /* Calculate the fractional part of the volume flux, using the     */
  /* Van Leer algorithm to get the tracer value on the face.         */
  ff_sl_van_leer(Fz, tr_moc, w, crfzf, clzf, 1, wincon->vc, wincon->s1, window->zp1, window->zm1);

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    ci = (w[c] > 0.0) ? c : window->zm1[c];
    c2 = (w[c] > 0.0) ? window->zm1[c] : c;
    if (wincon->fillf & CLIP) sminz[c] = smaxz[c] = tr[c2];
    if (wincon->ultimate) sminz[c] = smaxz[c] = tr_moc[c2];
    if (w[c] == 0.0) {
      Fz[c] = 0.0;
      continue;
    }
    ci = clzf[c];
    Fz[c] *= (crfzf[c] * wincon->dz[ci]);

    /* Integrate over the "integer" component of the trajectory      */
    ci = (w[c] > 0.0) ? window->zm1[c] : wincon->s1[cc];
    dz = wincon->dz[ci];
    dist = fabs(w[c] * dtu);
    while (dz < dist) {
      /* Get the minimum and maximum values encountered in the       */
      /* stencil along the streamline.                               */
      if (wincon->ultimate) {
	sminz[c] = min(sminz[c], tr_moc[ci]);
	smaxz[c] = max(smaxz[c], tr_moc[ci]);
      } else if (wincon->fillf & CLIP) {
	sminz[c] = min(sminz[c], tr[ci]);
	smaxz[c] = max(smaxz[c], tr[ci]);
      }
      Fz[c] += (dz * tr_moc[ci]);
      dist -= dz;
      ci = (w[c] < 0.0) ? window->zp1[ci] : window->zm1[ci];
      if (wincon->dz[ci] > 0)
	dz = wincon->dz[ci];
    }

    /* Ensure Fz has correct sign */
    if (w[c] < 0)
      Fz[c] *= -1.0;

    /* Universal flux limiter Eq 34 & 35 of Thuburn (1995)           */
    if (wincon->ultimate) {
      Fz[e] = min(Fz[c], smaxz[c] * w[c] * dtu);
      Fz[e] = max(Fz[c], sminz[c] * w[c] * dtu);
    }
  }

  /* Limit the mean tracer value if rquired                          */
  if (wincon->ultimate) {
    memset(mask, 0, window->sze * sizeof(int));
    /* Limit the horizontal fluxes in the order of the sparse        */
    /* vector. Some edges may be missed doing this if a cell is      */
    /* limited and the incoming flux is an outgoing flux for a       */
    /* subsequent cell.                                              */
    /*
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      int co, cs, cb;
      c = wincon->s1[cc];
      co = wincon->i3[cc];
      cs = wincon->i2[cc];
      cb = wincon->i1[cc];
      if(c==599)printf("a %f %d %d %d\n",windat->days,c,co,cs);
      while (c != window->zm1[c]) {
	universal_limit(window, windat, wincon, c, co, cs, cb, dtu, tr, Fx, Fz, mask);
	c = window->zm1[c];
      }
    }

    /* Limit the horizontal fluxes following outgoing flow from      */
    /* cells.                                                        */
    /*
    for (cc = 1; cc <= window->b2_t; cc++) {
      int cco, co, cs, cb;
      c = wincon->i6[cc];
      if (!(window->cask[c] & (W_WET|W_AUX))) continue;
      cco = window->c2cc[c];
      co = window->sur_t[cco];
      cs = window->nsur_t[cco];
      cb = window->bot_t[cco];
      while (c != window->zm1[c]) {
	universal_limit(window, windat, wincon, c, co, cs, cb, dtu, tr, Fx, Fz, mask);
	c = window->zm1[c];
      }
    }
    */
    for (cc = 1; cc <= wincon->s6[0]; cc++) {
      int cco, co, cs, cb, c2;
      c = wincon->s6[cc];
      c2 = window->m2d[c];
      if (!(window->cask[c] & (W_WET|W_AUX))) continue;
      cco = wincon->s7[c2];
      co = wincon->i3[cco];
      cs = max(co, wincon->i2[cco]);
      cb = wincon->i1[cco];
      /*
      co = window->sur_t[cco];
      cs = window->nsur_t[cco];
      cb = window->bot_t[cco];
      */
      /*if(c2==1&&windat->days>0.8)printf("aa %f %d %d %d %d\n",windat->days,cco,c,co,cs);*/
      universal_limit(window, windat, wincon, c, co, cs, cb, dtu, tr, Fx, Fz, mask);
    }
  }
}

/* END ffsl_don()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialises velocity field for FFSL                               */
/*-------------------------------------------------------------------*/
void ffsl_init(geometry_t *window,  /* Window geometry               */
	       window_t *windat,    /* Window data                   */
	       win_priv_t *wincon   /* Window constants              */
  )
{
  double *nu, *nv;
  double *dzz = wincon->dzz;
  double *nw = wincon->nw;
  double *mask = wincon->tr_mod;
  double *u1, *u2;
  int cc, c, ci, c2;
  int ee, e, es;
  int n, j, zm1;
  int tmode = (wincon->means & TRANSPORT) ? 1 : 0;

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing                                   */
  set_dzz(window, dzz);

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face / center          */
  if (tmode) {
    nu = windat->u1m;
    nv = windat->u2m;
    u1 = windat->ume;
    vel_w_trans(window, windat, wincon);
  } else {
    nu = windat->u;
    nv = windat->v;
    u1 = windat->u1;
  }
  vel_center_w(window, windat, wincon, nw);

  /* Set the ghost cells                                             */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    j = window->dbpt[cc];
    e = window->c2e[j][c];
    es = window->m2de[e];

    /* Rotate the edge velocity to east and north, and use this at   */
    /* at the ghost cell.                                            */
    nu[c] = u1[e] * window->costhu1[es] + windat->u2[e] * window->costhu2[es];
    nv[c] = u1[e] * window->sinthu1[es] + windat->u2[e] * window->sinthu2[es];
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value                                          */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    zm1 = window->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
  }
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    ci = window->zp1[c];
    nu[c] = nu[ci];
    nv[c] = nv[ci];
    nw[c] = nw[ci];
  }

  /* Set the open boundary ghost cells                               */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    int ee;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      c2 = open->obc_e2[ee];
      do {
	nu[c] = nu[c2];
	nv[c] = nv[c2];
	nw[c] = nw[c2];
	c = open->omape[ee][c];
      } while (c != open->omape[ee][c]);
    }
    /*OBC_bgz_nogradb(window, open, tr);*/
  }

  /* Set the cell width  mask                                        */
  for (c = 1; c <= window->enon; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bin[cc];
    mask[c] = 0.2;
  }

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    wincon->s2[cc] = c;
  }

  /*-----------------------------------------------------------------*/
  /* Initialize the interpolation structure                          */
  delaunay_reinit(window, windat->d, 0, nu);
  delaunay_reinit(window, windat->d, 1, nv);
  delaunay_reinit(window, windat->d, 2, nw);

  /* Make a contiguous mapping of cell locations                     */
  if (wincon->ultimate)
    mesh_expand_3d(window, u1);
}

/* END ffsl_init()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prepares arrays for FFSL with inline mean velocity transport      */
/*-------------------------------------------------------------------*/
void ffsl_trans_prep(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
		     )
{
  int tmode = (wincon->means & TRANSPORT) ? 1 : 0;

  if (!tmode) return;

  /* When using the flux-form semi-lagrange scheme, the velocities   */
  /* must be set to zero on the first time step of the transport     */
  /* model run. On the first time step, the transport model uses     */
  /* the same value for the old elevation (oldeta) and the new       */
  /* elevation (eta), so in effect the sea surface elevation does    */
  /* not change on the first time step. The flux divergence          */
  /* therefore also needs to be zero for the first time step only,   */
  /* which is achieved by setting u1, w and volume fluxes = 0.       */
  if (!wincon->suro[0]) {
    memset(windat->ume, 0, window->sze * sizeof(double));
    memset(windat->wm, 0, window->szc * sizeof(double));
    memset(windat->u1vm, 0, window->sze * sizeof(double));
    memcpy(wincon->etao, windat->eta, window->szcS * sizeof(double));
    memcpy(wincon->suro, window->sur_t, (window->a2_t + 1) * sizeof(int));
    wincon->suro[0] = 1;
  }

  /* Save the old surface from the 3D step in buffers                */
  memcpy(wincon->s5, window->sur_t, (window->a2_t + 1) * sizeof(int));
  memcpy(wincon->d4, wincon->oldeta, window->szcS * sizeof(double));
  memcpy(wincon->dzo, wincon->dz, window->szc * sizeof(double));

  /* Set the old elevation from the tracer step for dz calculation   */
  memcpy(window->nsur_t, wincon->suro, (window->a2_t + 1) * sizeof(int));
  memcpy(wincon->oldeta, wincon->etao, window->szcS * sizeof(double));

  /* Calculate layer thickness for the tracer step                   */
  set_dz(window, windat, wincon);

  /* Save the surface into tracer step buffers                       */
  memcpy(wincon->suro, window->nsur_t, (window->a2_t + 1) * sizeof(int));
  memcpy(wincon->etao, windat->eta, window->szcS * sizeof(double));
}

/* END ffsl_trans_prep()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reverts arrays for FFSL with inline mean velocity transport       */
/*-------------------------------------------------------------------*/
void ffsl_trans_post(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
		     )
{
  if (windat->dttr == 0.0) return;
  memcpy(window->sur_t, wincon->s5, (window->a2_t + 1) * sizeof(int));
  memcpy(wincon->oldeta, wincon->d4, window->szcS * sizeof(double));
  memcpy(wincon->dz, wincon->dzo, window->szc * sizeof(double));
}

/* END ffsl_trans_post()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Universal flux limiter, see Thuburn (1996), JCP.                  */
/* Assumes the upstream biased minimum and maximum values are stored */
/* smin and smax. Limits horizontal fluxes only.                     */
/*-------------------------------------------------------------------*/
void universal_limith(geometry_t *window,       /* Window geometry    */
		     window_t *windat,         /* Window data        */
		     win_priv_t *wincon,       /* Window constants   */
		     int c,
		     int co,
		     int cs,
		     int cb,
		     double dt,
		     double *tr,
		     double *Fx,
		     double *Fz,
		     int *mask
		     )
{
  int cv, e, es, ee;
  int c2 = window->m2d[c];
  int zp1 = window->zp1[c];
  double cout, cin;             /* Sum of out / inflow volume fluxes */
  double qn, qx;                /* Minimum / maximum tracer values   */
  double cn, cx;                /* Minimum / maximum fluxes          */
  double qmin, qmax;            /* Outflow bounding tracer values    */
  int ni, no;                   /* Number of inflow / outflow edges  */
  int sf[window->npem + 1];     /* Sign of flux (1=out, -1=in)       */
  double vol, volo;             /* Surface volume at time t and t-1  */
  double tro;                   /* Tracer mass at t-1                */
  double flux, dtracer, vflux;  /* Horizontal and vertical fluxes    */
  double *hflux;                /* Horizontal volume flux pointer    */
  double *smin = wincon->crfxc; /* Minimum streamline value          */
  double *smax = wincon->crfyc; /* Maximum streamline value          */
  double *sminz = wincon->crfxf; /* Minimum vertical streamline value*/
  double *smaxz = wincon->crfyf; /* Maximum vertical streamline value*/
  double *osubeta = wincon->d1;  /* Surface elevation at t-1         */
  double *subeta = wincon->d2;   /* Surface elevation at t           */
  double *tr_mod = wincon->tr_mod; /* z transverse tracer            */
  double *tr_moc = wincon->Fzh;    /* Horizontal transverse tracer   */
  double *w;                     /* Vertical velocity                */
  double d1;                     /* Dummy                            */
  int checkf = 0;
  int dolimit = 1;
  int cd = 0;

  if (checkf && window->s2i[c] == 47&& window->s2j[c] == 25 && window->s2k[c] == 22) cd = c;

  /* Set the pointers                                                */
  if (wincon->means & TRANSPORT) {
    hflux = windat->u1vm;
    w = windat->wm;
  } else {
    hflux = windat->u1flux3d;
    w = windat->w;
  }

  /* Get the cell volumes at time t and t-1                          */
  if (c == cs) {
    double watertop;
    double waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
    double dz = subeta[c2] - waterbot;
    vol = dz * window->cellarea[c2];
    /* First increment the old volume from the current layer to the  */
    /* layer below that containing the old surface. If elevation     */
    /* has risen this loop is skipped.                               */
    volo = tro = 0.0;
    cv = c;
    while (c != co) {
      zp1 = window->zp1[c];
      watertop = window->gridz[zp1];
      d1 = window->cellarea[c2] * (watertop - waterbot);
      volo += d1;
      tro += tr[c] * d1;
      c = zp1;
      waterbot = window->gridz[c];
    }

    /* Now increment the old volume for the layer containing the old */
    /* surface.                                                        */
    waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
    d1 = (osubeta[c2] - waterbot) * window->cellarea[c2];
    volo += d1;
    tro += tr[c] * d1;
  } else {
    cv = c;
    vol = volo = window->cellarea[c2] * wincon->dz[c];
    tro = tr[c] * volo;
  }

  /* Diagnostics                                                     */
  if(checkf && c == cd) {
    printf("-----------------------------\n");
    printf("start %s %f %d %d[%d %d %d](%f %f) cs=%d co=%d tr=%f dt=%f\n",windat->days,
	   wincon->trinfo_3d[(int)tr[0]].name, 
	   windat->nstep, c, window->s2i[c], window->s2j[c], window->s2k[c],
	   window->cellx[c], window->celly[c], cs, co, tr[c], dt);
    printf("eta(t)=%f eta(t-1)=%f volo=%f vol=%f\n", subeta[c2], osubeta[c2],
	   volo/window->cellarea[c2], vol/window->cellarea[c2]);
  }

  /* Get the sum of fluxes                                           */
  qn = qx = tr_mod[c];
  cn = cx = cout = cin = 0.0;
  ni = no = 0;

  /* Increment the fluxes from the current layer to the layer below  */
  /* that containing the old surface. If elevation has risen through */
  /* layers then this loop is skipped.                               */
  c = cv;
  if (cv == cs) {
    while (c != co) {
      zp1 = window->zp1[c];
      for (ee = 1; ee <= window->npe[c2]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	flux = window->eSc[ee][c2] * hflux[e] * dt;
	if(c==cd)printf("drop\n");
	if (flux > 0.0) {
	  cout += flux;
	  no++;
	} else if (flux < 0.0) {
	  mask[e] = 1;
	  qn = min(qn, smin[e]);
	  qx = max(qx, smax[e]);
	  cn -= smin[e] * flux;
	  cx -= smax[e] * flux;
	  cin -= flux;
	  ni++; 
	}
      }
      c = zp1;
    }
  }

  /* Increment the fluxes for the layer containing the old surface   */
  for (ee = 1; ee <= window->npe[c2]; ee++) {
    e = window->c2e[ee][c];
    es = window->m2de[e];
    /* Get the edge signs                                            */
    sf[ee] = 0;
    flux = window->eSc[ee][c2] * hflux[e] * dt;
    if (flux > 0.0) {              /* Outflow                        */
      sf[ee] = 1;
      cout += flux;
      if(c==cd)printf("aout %d e=%d min=%f Fx=%f\n",ee, e, 
		      smin[e] * flux, Fx[e]);
      no++;
    } else if (flux < 0.0) {
      sf[ee] = -1;
      mask[e] = 1;
      qn = min(qn, smin[e]);
      qx = max(qx, smax[e]);
      /* Thuburn (1996) Eq. 48 & 49 */
      cn -= min(window->eSc[ee][c2] * Fx[e]*dt, smin[e] * flux);
      cx -= max(window->eSc[ee][c2] * Fx[e]*dt, smax[e] * flux);
      cin -= flux;
      if(c==cd)printf("ain %d e=%d min=%f max=%f Fx=%f\n",ee,e,
		      cn, cx, Fx[e]*dt,smin[e]*flux);
      ni++; 
    }
  }

  /* Increment the fluxes from above the old surface elevation to    */
  /* the top of the vertical grid. If the surface has risen then     */
  /* these cells now contain water and horizontal divergence may     */
  /* alter the total amount of tracer. If elevation has dropped then */
  /* there is no water (and hence no divergence) in these cells.     */
  if (cv == cs && c != c2) {
    do {
      c = window->zp1[c];
      for (ee = 1; ee <= window->npe[c2]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	flux = window->eSc[ee][c2] * hflux[e] * dt;
	if(c==cd)printf("rise\n");
	if (flux > 0.0) {
	  cout += flux;
	  no++;
	} else if (flux < 0.0) {
	  mask[e] = 1;
	  qn = min(qn, smin[e]);
	  qx = max(qx, smax[e]);
	  cn -= smin[e] * flux;
	  cx -= smax[e] * flux;
	  cin -= flux;
	  ni++; 
	}
      }
    } while (c != c2);
  }
  c = cv;

  if (!cout) return;
  if (!ni) {
    qn = qx = tr[c];
  }

  /* Get the vertical fluxes for this cell                           */
  if (c == cs) 
    vflux = - Fz[c] * window->cellarea[c2] ;
  else
    vflux = (Fz[zp1] - Fz[c]) * window->cellarea[c2];

  /* Get the outflow bounding values                                 */
  qmin = (tro - vol * qx + cx - vflux) / (cout);
  qmax = (tro - vol * qn + cn - vflux) / (cout);

  /* Check                                                           */
  if (c == cd) {
    double d2, d3;
    double dtracern = 0.0;
    double dtracerx = 0.0;
    dtracer = 0.0;
    /* Get the continuity (volume) balance. This is independent of   */
    /* tracer value and checks volume fluxes and sea levels are      */
    /* correct. Should be ~1e-10.                                    */
    if (c == cs) {
      d2 = -w[c] * window->cellarea[c2] * dt;
    } else
      d2 = (w[zp1] - w[c]) * window->cellarea[c2] * dt;
    d3 = vflux / tr[c];
    d1 = window->cellarea[c2] * (osubeta[c2] - subeta[c2]) - ((cout-cin) + d2);
    d1 = (volo - vol) - ((cout-cin) + d2);
    printf("cnt=%e vflux=%f vflux=%f hflux=%f\n", d1, d2, d3, cout-cin);
    /* Print the vertical mass fluxes                                */
    printf("  vert flux wm=%f zp1=%f z=%f vflux=%f tr=%f\n",w[c], Fz[zp1], Fz[c], vflux, vflux / d2);
    /* Print the horizontal mass fluxes                              */
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      flux = window->eSc[ee][c2] * hflux[e] * dt;
      dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
      dtracerx += qmin * flux;
      dtracern += qmax * flux;
      if (sf[ee] == 1)printf("  out%d Fx=%f vol=%f min=%f max=%f tr=%f\n",ee, 
			     window->eSc[ee][c2] * Fx[e], hflux[e], hflux[e] * smin[e], 
			     hflux[e] * smax[e], Fx[e]/hflux[e]);
      if (sf[ee] == -1)printf("  in%d Fx=%f vol=%f min=%f max=%f tr=%f\n",ee, 
			      window->eSc[ee][c2] * Fx[e], hflux[e], hflux[e] * smin[e], 
			      hflux[e] * smax[e], Fx[e]/hflux[e]);
    }
    d1 = (volo * tr[c] - (dtracer + vflux)) / vol;
    d2 = (volo * tr[c] - (dtracern + vflux)) / vol;
    d3 = (volo * tr[c] - (dtracerx + vflux)) / vol;
    windat->dum1[c] = d1;
    printf("  cn=%f cx=%f min=%f max=%f tr=%f qmin=%f qmax=%f\n",cn, cx, qn, qx, d1, qmin, qmax);
    printf("  tr=%f trn=%f trx=%f\n",d1, d2, d3);
  }

  /* Limit the outflow edges                                         */
  if (dolimit) {
    dtracer = 0.0;

    /* Limit the fluxes from the current layer to the layer below    */
    /* that containing the old surface.                              */
    if (cv == cs) {
      while (c != co) {
	zp1 = window->zp1[c];
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][c];
	  if (!mask[e] && window->eSc[ee][c2] * hflux[e] > 0) {
	    mask[e] = 1;
	    d1 = window->eSc[ee][c2] * hflux[e];
	    Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	    Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	  }
	  dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
	}
	c = zp1;
      }
    }

    /* Limit the fluxes for the layer containing the old surface     */
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      if (!mask[e] && sf[ee] == 1) {    /* Outflow                   */
	mask[e] = 1;
	d1 = window->eSc[ee][c2] * hflux[e];
	flux =  Fx[e];
	if (c == cd) printf("limit %f : ",Fx[e]);
	Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	if (c == cd) printf("%f : ",Fx[e]);
	Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	if (c == cd) printf("%f\n",Fx[e]);
	/*if (fabs(flux - Fx[e]) > 1e-8) hd_quit("limit at c=%d d=%d %f %f\n",c, e, flux, Fx[e]);*/
      }
      dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
    }

    /* Limit the fluxes from above the old surface elevation to the  */
    /* top of the vertical grid.                                     */
    if (cv == cs && c != c2) {
      do {
	c = window->zp1[c];
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][c];
	  if (!mask[e] && window->eSc[ee][c2] * hflux[e] > 0) {
	    mask[e] = 1;
	    d1 = window->eSc[ee][c2] * hflux[e];
	    Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	    Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	  }
	  dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
	}
      } while (c != c2);
    }
    c = cv;
  }

  /* Check                                                           */
  if (checkf && c == cd) {
    printf("new dtracer=%f tr=%f\n", dtracer, d1);
  }
}

/* END universal_limit()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Universal flux limiter, see Thuburn (1996), JCP.                  */
/* Assumes the upstream biased minimum and maximum values are stored */
/* in smin and smax. Limits horizontal and vertical fluxes.          */
/*-------------------------------------------------------------------*/
void universal_limit(geometry_t *window,       /* Window geometry    */
		     window_t *windat,         /* Window data        */
		     win_priv_t *wincon,       /* Window constants   */
		     int ci,
		     int co,
		     int cs,
		     int cb,
		     double dt,
		     double *tr,
		     double *Fx,
		     double *Fz,
		     int *mask
		     )
{
  int c, e, es, ee;
  int c2 = window->m2d[ci];
  int zp1 = window->zp1[ci];
  double cout, cin;             /* Sum of out / inflow volume fluxes */
  double qn, qx;                /* Minimum / maximum tracer values   */
  double cn, cx;                /* Minimum / maximum fluxes          */
  double qmin, qmax;            /* Outflow bounding tracer values    */
  int ni, no;                   /* Number of inflow / outflow edges  */
  int sf[window->npem + 1];     /* Sign of flux (1=out, -1=in)       */
  double vol, volo;             /* Surface volume at time t and t-1  */
  double tro;                   /* Tracer mass at t-1                */
  double flux, dtracer, vflux;  /* Horizontal and vertical fluxes    */
  double *hflux;                /* Horizontal volume flux pointer    */
  double *smin = wincon->crfxc; /* Minimum streamline value          */
  double *smax = wincon->crfyc; /* Maximum streamline value          */
  double *sminz = wincon->crfxf; /* Minimum vertical streamline value*/
  double *smaxz = wincon->crfyf; /* Maximum vertical streamline value*/
  double *osubeta = wincon->d1;  /* Surface elevation at t-1         */
  double *subeta = wincon->d2;   /* Surface elevation at t           */
  double *tr_mod = wincon->tr_mod; /* z transverse tracer            */
  double *tr_moc = wincon->Fzh;    /* Horizontal transverse tracer   */
  double *w;                     /* Vertical velocity                */
  double d1;                     /* Dummy                            */
  int checkf = 0;
  int dolimit = 1;
  int cd = 0;

  c = ci;
  /*if (checkf && window->s2i[c] == 30 && window->s2j[c] == 10 && window->s2k[c] == 22) cd = c;*/
  if (checkf && c == 1) cd = c;

  /* Set the pointers                                                */
  if (wincon->means & TRANSPORT) {
    hflux = windat->u1vm;
    w = windat->wm;
  } else {
    hflux = windat->u1flux3d;
    w = windat->w;
  }

  /* Get the cell volumes at time t and t-1                          */
  if (c == cs) {
    double watertop;
    double waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
    double dz = subeta[c2] - waterbot;
    vol = dz * window->cellarea[c2];
    /* First increment the old volume from the current layer to the  */
    /* layer below that containing the old surface. If elevation     */
    /* has risen this loop is skipped.                               */
    volo = tro = 0.0;
    while (c != co) {
      zp1 = window->zp1[c];
      watertop = window->gridz[zp1];
      d1 = window->cellarea[c2] * (watertop - waterbot);
      volo += d1;
      tro += tr[c] * d1;
      c = zp1;
      waterbot = window->gridz[c];
    }
    if(checkf && c == cd) printf("a2 %d %d\n",c,co);
    /* Now increment the old volume for the layer containing the old */
    /* surface.                                                        */
    waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
    d1 = (osubeta[c2] - waterbot) * window->cellarea[c2];
    volo += d1;
    tro += tr[c] * d1;
  } else {
    vol = volo = window->cellarea[c2] * wincon->dz[c];
    tro = tr[c] * volo;
  }

  /* Diagnostics                                                     */
  if(checkf && c == cd) {
    printf("-----------------------------\n");
    printf("start %s %f %d %d[%d %d %d](%f %f) cs=%d co=%d tr=%f dt=%f\n",windat->days,
	   wincon->trinfo_3d[(int)tr[0]].name, 
	   windat->nstep, c, window->s2i[c], window->s2j[c], window->s2k[c],
	   window->cellx[c2], window->celly[c2], cs, co, tr[c], dt);
    printf("eta(t)=%f eta(t-1)=%f volo=%f vol=%f\n", subeta[c2], osubeta[c2],
	   volo/window->cellarea[c2], vol/window->cellarea[c2]);
  }

  /* Get the sum of fluxes                                           */
  qn = min(tr_mod[c], tr_moc[c]);
  qx = max(tr_mod[c], tr_moc[c]);
  cn = cx = cout = cin = 0.0;
  ni = no = 0;
  dtracer = 0.0;
  c = ci;

  /* Increment the fluxes from the current layer to the layer below  */
  /* that containing the old surface. If elevation has risen through */
  /* layers then this loop is skipped.                               */
  if (c == cs) {
    while (c != co) {
      zp1 = window->zp1[c];
      for (ee = 1; ee <= window->npe[c2]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	flux = window->eSc[ee][c2] * hflux[e] * dt;
	dtracer += window->eSc[ee][c2] * Fx[e] * dt;
	if(c==cd)printf("drop\n");
	if (flux > 0.0) {
	  cout += flux;
	  no++;
	} else if (flux < 0.0) {
	  mask[e] = 1;
	  qn = min(qn, smin[e]);
	  qx = max(qx, smax[e]);
	  cn -= smin[e] * flux;
	  cx -= smax[e] * flux;
	  cin -= flux;
	  ni++; 
	}
      }
      c = zp1;
    }
  }

  /* Increment the fluxes for the layer containing the old surface   */
  for (ee = 1; ee <= window->npe[c2]; ee++) {
    e = window->c2e[ee][c];
    es = window->m2de[e];
    /* Get the edge signs                                            */
    sf[ee] = 0;
    flux = window->eSc[ee][c2] * hflux[e] * dt;
    dtracer += window->eSc[ee][c2] * Fx[e] * dt;
    if (flux > 0.0) {              /* Outflow                        */
      sf[ee] = 1;
      cout += flux;
      if(c==cd)printf("aout %d e=%d min=%f Fx=%f\n",ee, e, 
		      smin[e] * flux, Fx[e]);
      no++;
    } else if (flux < 0.0) {
      sf[ee] = -1;
      mask[e] = 1;
      qn = min(qn, smin[e]);
      qx = max(qx, smax[e]);
      /* Thuburn (1996) Eq. 48 & 49 */
      cn -= min(window->eSc[ee][c2] * Fx[e]*dt, smin[e] * flux);
      cx -= max(window->eSc[ee][c2] * Fx[e]*dt, smax[e] * flux);
      cin -= flux;
      if(c==cd)printf("ain %d e=%d min=%f max=%f Fx=%f\n",ee,e,
		      cn, cx, Fx[e]*dt);
      ni++; 
    }
  }

  /* Increment the fluxes from above the old surface elevation to    */
  /* the top of the vertical grid. If the surface has risen then     */
  /* these cells now contain water and horizontal divergence may     */
  /* alter the total amount of tracer. If elevation has dropped then */
  /* there is no water (and hence no divergence) in these cells.     */
  if (ci == cs && c != c2) {
    do {
      c = window->zp1[c];
      for (ee = 1; ee <= window->npe[c2]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	flux = window->eSc[ee][c2] * hflux[e] * dt;
	dtracer += window->eSc[ee][c2] * Fx[e] * dt;
	if(c==cd)printf("rise\n");
	if (flux > 0.0) {
	  cout += flux;
	  no++;
	} else if (flux < 0.0) {
	  mask[e] = 1;
	  qn = min(qn, smin[e]);
	  qx = max(qx, smax[e]);
	  cn -= smin[e] * flux;
	  cx -= smax[e] * flux;
	  cin -= flux;
	  ni++; 
	}
      }
    } while (c != c2);
  }

  /* Vertical fluxes into the cell                                   */
  c = ci;
  flux = w[c] * window->cellarea[c2] * dt;
  if (flux > 0.0) {
    qn = min(qn, sminz[c]);
    qx = max(qx, smaxz[c]);
    /* Thuburn (1996) Eq. 48 & 49 */
    cn += min(Fz[c] * window->cellarea[c2], sminz[c] * flux);
    cx += max(Fz[c] * window->cellarea[c2], smaxz[c] * flux);
    cin += flux;
    ni++;
  } else if (flux < 0.0) {
    cout -= flux;
    no++;
  }
  if (c != cs) {
    zp1 = window->zp1[c];
    flux = w[zp1] * window->cellarea[c2] * dt;
    if (flux < 0.0) {
      qn = min(qn, sminz[zp1]);
      qx = max(qx, smaxz[zp1]);
      cn -= min(Fz[zp1] * window->cellarea[c2], sminz[zp1] * flux);
      cx -= max(Fz[zp1] * window->cellarea[c2], smaxz[zp1] * flux);
      cin -= flux;
      ni++;
    } else if (flux > 0.0) {
      cout += flux;
      no++;
    }
  }
  if (!cout) return;
  if (!ni) {
    qn = qx = tr[c];
  }

  /* Get the horizontal outflow bounding values                      */
  qmin = (tro - vol * qx + cx) / (cout);
  qmax = (tro - vol * qn + cn) / (cout);

  /* Check                                                           */
  if (c == cd) {
    double d2, d3;
    dtracer = 0.0;
    /* Get the continuity (volume) balance. This is independent of   */
    /* tracer value and checks volume fluxes and sea levels are      */
    /* correct. Should be ~1e-10.                                    */
    d1 = (volo - vol) - ((cout-cin));
    printf("cnt=%e flux=%f\n", d1, cout-cin);
    /* Print the vertical mass fluxes                                */

    vflux = - Fz[c] * window->cellarea[c2] ;
    flux = w[c] * dt;
    if (Fz[c] >= 0.0)
      printf("  vin Fz=%f vol=%f min=%f max=%f tr=%f\n", Fz[c], flux, 
	     flux * sminz[c], flux * smaxz[c], Fz[c]/flux);
    else
      printf("  vout Fz=%f vol=%f min=%f max=%f tr=%f\n", Fz[c], flux, 
	     flux * sminz[c], flux * smaxz[c], Fz[c]/flux);
    if (c != cs) {
      vflux += Fz[zp1] * window->cellarea[c2];
      flux = w[zp1] * window->cellarea[c2] * dt;

      if (Fz[zp1] < 0.0)
	printf("  vp1in Fzp1=%f vol=%f min=%f max=%f tr=%f\n", Fz[zp1], flux, 
	       flux * sminz[zp1], flux * smaxz[zp1], Fz[zp1]/flux);
      else
	printf("  vp1out Fzp1=%f vol=%f min=%f max=%f tr=%f\n", Fz[zp1], flux, 
	       flux * sminz[zp1], flux * smaxz[zp1], Fz[zp1]/flux);
    }
    /* Print the horizontal mass fluxes                              */
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      flux = window->eSc[ee][c2] * hflux[e] * dt;
      dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
      if (sf[ee] == 1)printf("  out%d Fx=%f vol=%f min=%f max=%f tr=%f\n",ee, 
			     window->eSc[ee][c2] * Fx[e], hflux[e], hflux[e] * smin[e], 
			     hflux[e] * smax[e], Fx[e]/hflux[e]);
      if (sf[ee] == -1)printf("  in%d Fx=%f vol=%f min=%f max=%f tr=%f\n",ee, 
			      window->eSc[ee][c2] * Fx[e], hflux[e], hflux[e] * smin[e], 
			      hflux[e] * smax[e], Fx[e]/hflux[e]);
    }

    d1 = (volo * tr[c] - (dtracer + vflux)) / vol;
    printf("  cn=%f cx=%f min=%f max=%f tr=%f qmin=%f qmax=%f\n",cn, cx, qn, qx, d1, qmin, qmax);
  }

  /* Limit the outflow edges                                         */
  if (dolimit) {
    dtracer = 0.0;
    /* Limit the fluxes from the current layer to the layer below    */
    /* that containing the old surface.                              */
    if (c == cs) {
      while (c != co) {
	zp1 = window->zp1[c];
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][c];
	  if (!mask[e] && window->eSc[ee][c2] * hflux[e] > 0) {
	    mask[e] = 1;
	    d1 = window->eSc[ee][c2] * hflux[e];
	    Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	    Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	  }
	  dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
	}
	c = zp1;
      }
    }

    /* Limit the fluxes for the layer containing the old surface     */
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      if (!mask[e] && sf[ee] == 1) {    /* Outflow                   */
	mask[e] = 1;
	d1 = window->eSc[ee][c2] * hflux[e];
	flux =  Fx[e];
	if (c == cd) printf("limit %f : ",Fx[e]);
	Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	if (c == cd) printf("%f : ",Fx[e]);
	Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	if (c == cd) printf("%f\n",Fx[e]);
	/*if (fabs(flux - Fx[e]) > 1e-8) hd_quit("limit at c=%d d=%d %f %f\n",c, e, flux, Fx[e]);*/
      }
      dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
    }

    /* Limit the fluxes from above the old surface elevation to the  */
    /* top of the vertical grid.                                     */
    if (ci == cs && c != c2) {
      do {
	c = window->zp1[c];
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][c];
	  if (!mask[e] && window->eSc[ee][c2] * hflux[e] > 0) {
	    mask[e] = 1;
	    d1 = window->eSc[ee][c2] * hflux[e];
	    Fx[e] = min(Fx[e] * window->eSc[ee][c2], qmax * d1) * window->eSc[ee][c2];
	    Fx[e] = max(Fx[e] * window->eSc[ee][c2], qmin * d1) * window->eSc[ee][c2];
	  }
	  dtracer += (window->eSc[ee][c2] * Fx[e] * dt);
	}
      } while (c != c2);
    }

    /* Limit the vertical outgoing fluxes                            */
    c = ci;
    if (Fz[c] < 0) {
      if (c == cd) printf("limit lower %f : ",Fz[c]);
      Fz[c] = max(Fz[c], qmax * w[c] * dt);
      if (c == cd) printf("%f : ",Fz[c]);
      Fz[c] = min(Fz[c], qmin * w[c] * dt);
      if (c == cd) printf("%f\n",Fz[c]);
    }
    zp1 = window->zp1[c];
    if (c != cs && Fz[zp1] > 0.0) {
      if (c == cd) printf("limit upper %f : ",Fz[zp1]);
      Fz[zp1] = min(Fz[zp1], qmax * w[zp1] * dt);
      if (c == cd) printf("%f : ",Fz[zp1]);
      Fz[zp1] = max(Fz[zp1], qmin * w[zp1] * dt);
      if (c == cd) printf("%f\n",Fz[zp1]);
    }
  }

  /* Check                                                           */
  if (checkf && c == cd) {
    d1 = (tro - (dtracer + vflux)) / vol;
    printf("new dtracer=%f tr=%f\n", dtracer, d1);
  }
}

/* END universal_limit()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* gets the local bounds for cell c, i.e. the minimum and maximum    */
/* tracer values at all stencil points used to interpolate a value.  */
/*-------------------------------------------------------------------*/
void get_local_bounds(geometry_t *window, 
		      window_t *windat, 
		      win_priv_t *wincon, 
		      double tr, 
		      double *trv,
		      int e,               /* Destination edge       */
		      int co,              /* Source centre          */
		      double *smin,
		      double *smax,
		      int mode /* 0 = stencil, 1 = interpolated value*/
		      )
{
  qweights *lw = &wincon->lw[co];


  if (mode == 0) {
    smin[e] = min(smin[e], tr);
    smax[e] = max(smax[e], tr);
  } else if (mode == 1) {
    smin[e] = min(smin[e], lw->tmn);
    smax[e] = max(smax[e], lw->tmx);
  } else if (mode == 2) {
    int j, c;
    for (j = 0; j < lw->ncells; j++) {
      c = lw->cells[j];
      smin[e] = min(smin[e], trv[c]);
      smax[e] = max(smax[e], trv[c]);
    }
  }
}

/* END get_local_bounds()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks that a tracer value with streamline origin c does not      */
/* exceed local minima or maxima.                                    */
/*-------------------------------------------------------------------*/
void check_monotone(geometry_t *window,        /* Window geometry    */
		    window_t *windat,          /* Window data        */
		    win_priv_t *wincon,        /* Window constants   */
		    double tr,                 /* Tracer value       */
		    int c
		    )
{
  qweights *lw = &wincon->lw[c];
  double mx = lw->tmx;
  double mn = lw->tmn;

  int tn = wincon->monon;

  if (tr > mx) {
    /*if (windat->mono) windat->mono[c] = max(windat->mono[c], fabs(tr - mx));*/
    /*
    hd_quit("Non monotonic max %f @ %d(%d %d %d) %12.10f : %f %f\n",windat->days,c, window->s2i[c], window->s2j[c], window->s2k[c],mx,mn);
    */
  }
  if (tr < mn) {
    /*if (windat->mono) windat->mono[c] = max(windat->mono[c], fabs(tr - mn));*/
    /*
    hd_quit("Non monotonic min %f @ %d(%d %d %d) %12.10f : %f %f\n",windat->days,c, window->s2i[c], window->s2j[c], window->s2k[c],mx,mn);
    */
  }
}

/* END check_monotone()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clips tracer values to the minimum and maximum values set in the  */
/* quadtratic interpolation function. If FILL_METHOD = MONOTONIC,    */
/* any residual mass created or destroyed during clipping is         */
/* iteratively filled throughout the whole domain.                   */
/*-------------------------------------------------------------------*/
void clip_ffsl(geometry_t *window,        /* Window geometry         */
	       window_t *windat,          /* Window data             */
	       win_priv_t *wincon,        /* Window constants        */
	       double *tr                 /* Tracer array            */
	       )
{
  int cc, c, cs, c2, k, id, m, j, cn, ee, e;
  int itmax = 20;
  int *cl = wincon->s2;
  double tn1 = HUGE, tx1 = -HUGE;
  double tn2 = HUGE, tx2 = -HUGE;
  double ptr;
  double nfilled = 0;
  double *mask = wincon->tr_mod;
  double *mn = wincon->w5;
  double *mx = wincon->w3;
  double omass, emass;
  double msf = 0.0;
  double eps = 1e-5;
  double msfmax = 2.0;
  double v2, v3, vol;
  double *mxc = wincon->w4;    /* Local tracer maximum at source     */
  double *mnc = wincon->w5;    /* Local tracer minimum at source     */
  double m1, m2;
  double *smin = wincon->crfxc; /* Minimum streamline value          */
  double *smax = wincon->crfyc; /* Maximum streamline value          */
  double *sminz = wincon->crfxf; /* Minimum vertical streamline value*/
  double *smaxz = wincon->crfyf; /* Maximum vertical streamline value*/

  /* Get the min/max along all streamlines                           */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];

    mn[c] = sminz[c];
    mx[c] = smaxz[c];
    for (ee = 1; ee <= window->npe[c2]; ee++) {
      e = window->c2e[ee][c];
      mn[c] = min(mn[c], smin[e]);
      mx[c] = max(mx[c], smax[e]);
    }
  }
  if (wincon->fillf & CLIP) {
    for (cc = 1; cc <= window->v3_t; cc++) {
      c = window->w3_t[cc];
      tr[c] =  max(min(mx[c], tr[c]), mn[c]);
    }
    return;
  }

  /* Set ghost cells                                                 */
  /*
  set_tr_nograd(window, tr);
  */
  /*-----------------------------------------------------------------*/
  /* Clip tracer values                                              */
  omass = emass = 0.0;
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = cl[cc];
    c2 = window->m2d[cs];

    mn[c] = sminz[c];
    mx[c] = smaxz[c];
    for (ee = 1; ee <= window->npe[window->m2d[c]]; ee++) {
      e = window->c2e[ee][c];
      mn[c] = min(mn[c], smin[e]);
      mx[c] = max(mx[c], smax[e]);
    }

    /* Get the minimum and maximum from the quadratic interpolator   */
    /*
    k = window->s2k[cs] + 1;
    id = window->c2p[k][c2];
    if (id >= 0) {
      lsqq_minmax(windat->gst[k]->interpolator, id, &tn1, &tx1);
    }
    */
    /* Minimum and maximum in the layer above                        */
    /*
    k = window->s2k[window->zp1[cs]] + 1;
    id = window->c2p[k][c2];
    if (id >= 0) {
      lsqq_minmax(windat->gst[k]->interpolator, id, &tn2, &tx2);
    }
    mn[c] = min(tn1, tn2);
    mx[c] = max(tx1, tx2);

    mn[c] = wincon->lw[cs].tmn;
    mx[c] = wincon->lw[cs].tmx;
    */
    /* Clip the tracer and get mass before and after clipping        */
    vol = window->cellarea[window->m2d[c]] * wincon->dz[c];
    m1 = tr[c] * vol;
    omass += m1;
    tr[c] =  max(min(mx[c], tr[c]), mn[c]);
    m2 = tr[c] * vol;
    emass += m2;
    windat->vol_cons[window->m2d[c]] = 0.5 * (windat->vol_cons[window->m2d[c]] + (m1 - m2) / m2);

  }

  /*-----------------------------------------------------------------*/
  /* Redistribute residual mass. Mass is to distributed to cells     */
  /* that are situated in uniform tracer distributions (mn == mx),   */
  /* are at the maximum value (tr == mx) or at the minimum value     */
  /* (tr == mn). In these cases, residual mass is iteratively        */
  /* distributed to cells that do not satisfy these conditions.      */
  if (wincon->fillf & MONOTONIC) {
    /* Get the mass scaling factor                                   */
    msf = omass / emass;
    if (msf == 1.0) return;
    
    m = 0;
    memset(mask, 0, window->szc * sizeof(double));
    while (msf != 1.0 && m < itmax) {

      /* Sanity checks                                               */
      if (msf > msfmax) {
	hd_warn("global_fill: Tracer global scaling %f considered too large at %5.2f days: no filling\n",
		msf, windat->t/86400);
	return;
      }
      /* Re-distribute the change in mass                            */
      emass = 0.0;
      for(cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	vol = window->cellarea[cs] * wincon->dz[c];
	v2 = tr[c] * msf;
	v3 = fabs(mx[c] - mn[c]);
	if (v3 && v3 < eps && !mask[c]) {
	  mask[c] = 1.0;
	} else if (v2 > mx[c] && !mask[c]) {
	  /* Check for adjusted mass > local maximum                 */
	  tr[c] = mx[c];
	  mask[c] = 1.0;
	  /*if(c==3540)printf("fill %e %e\n",v2,mx[c]);*/
	} else if (v2 < mn[c] && !mask[c]) {
	  /* Check for adjusted mass < local minimum                 */
	  tr[c] = mn[c];
	  mask[c] = 1.0;
	} else if (!mask[c]) {
	  /* Adjust the mass                                         */
	  tr[c] = v2;
	}
	emass += tr[c] * vol;
      }
      /* Calculate the change in mass over the time-step.            */
      msf = (emass) ? omass / emass : 1.0;
      m++;
    }
    if (m >= itmax) 
      hd_warn("global_fill: Tracer did not converge at %5.2f, residual mass = %5.3f\n",
	      windat->t/86400, omass - emass);
  }
  set_tr_nograd(window, tr);
}

/* END clip_ffsl()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up a grid_specs structure for variables using LAGRANGE that  */
/* require interpolation.                                            */
/*-------------------------------------------------------------------*/
void tran_grid_init(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int cc, c, c2, ci, k;
  int *nk, nz = window->nz;
  point **p;
  double v;
  delaunay *d;
  int n, m, ee, e;
  int ef = 0;
  int sf = 3;  /* 3 = one row, 5 = two rows                          */
  int *mask;
  int newcode = 1;
  /*
  if (window->nwindows > 1)
    hd_quit("Lagrangian tracking does not operate with more than one window: exiting.\n");
  */
  /* Allocate                                                        */
  n = window->szcS + window->szvS;
  p = (point **)alloc_2d(n, nz+1, sizeof(point));
  /* c2p is a map from the mesh index c to the delaunay points index */
  window->c2p = i_alloc_2d(window->szcS, nz+1);
  window->gcm = i_alloc_2d(window->szcS, nz+1);
  for (k = 0; k <= nz; k++)
    for (c = 0; c < window->szcS; c++) {
      window->c2p[k][c] = -1;
      window->gcm[k][c] = L_OUT;
    }

  /*-----------------------------------------------------------------*/
  /* Set the ghost cell mask                                         */
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c];
    window->gcm[k][c2] = L_WET;
  }
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    c2 = window->m2d[c];
    k = window->s2k[ci];
    window->gcm[k][c2] = L_GHOST;
  }
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    c2 = window->m2d[c];
    ci = window->zm1[c];
    k = window->s2k[c] - 1;
    if (k >= 0) window->gcm[k][c2] = L_SED;
  }

  /*-----------------------------------------------------------------*/
  /* Set up the triangulation for each layer.                        */
  /* First count the number of points used in the triangulation in   */
  /* each layer. Note that the number of ghost cells used in the     */
  /* triangulation become less with depth.                           */
  /* For quadratic least squares make a mapping from the centre of   */
  /* the Voronoi dual to the Delaunay index (c2p).                   */
  nk = i_alloc_1d(nz+1);
  memset(nk, 0, (nz+1) * sizeof(int));
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    p[k][nk[k]].x = window->cellx[c2];
    p[k][nk[k]].y = window->celly[c2];
    window->c2p[k][c2] = nk[k];
    nk[k]++;
  }

  /* Ghost cells, in 3 dimensions.                                   */
  /* Although the mesh has unique ghost cells mapping across all     */
  /* edges, the triangulation used for interpolation does not.       */
  /* The geographic location of ghost cells is set to the cell edge  */
  /* so as to account for multiple ghost cells in the interpolation  */
  /* triangulation. This means the boundary condition for velocity   */
  /* is v=0 at ghost cells rather than v=-vi.                        */
  for (cc = window->a3_t+1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    ci = window->wgst[c];
    c2 = window->m2d[ci];
    k = window->s2k[ci] + 1;
    window->s2k[c] = window->s2k[ci];
    for (n = 1; n <= window->npe[c2]; n++) {
      if (c == window->c2c[n][ci])
	break;
    }
    e = window->m2de[window->c2e[n][ci]];
    p[k][nk[k]].x = window->u1x[e];
    p[k][nk[k]].y = window->u1y[e];
    window->c2p[k][window->m2d[c]] = nk[k];
    nk[k]++;
  }

  /* Open boundary ghosts                                            */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    int nbgz = (open->bgz) ? open->bgz : 1;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      ci = open->obc_e2[ee];
      for (m = 0; m < nbgz; m++) {
	c2 = window->m2d[c];
	k = window->s2k[ci] + 1;
	window->s2k[c] = window->s2k[ci];
	if (window->c2p[k][c2] == -1) {
	  p[k][nk[k]].x = window->cellx[c2];
	  p[k][nk[k]].y = window->celly[c2];
	  window->c2p[k][c2] = nk[k];

	  nk[k]++;
	}
	c = open->omape[ee][c];
      }
    }
  }
  wincon->osl |= L_FG;

  /* Cell vertices                                                   */
  if (ef) {
    mask = i_alloc_1d(window->szv);
    memset(mask, 0, window->szv * sizeof(int));
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      k = window->s2k[c] + 1;
      for (n = 1; n <= window->npe[c2]; n++) {
	int v = window->c2v[n][c];
	int vs = window->m2dv[v];
	if (v && !mask[v]) {
	  p[k][nk[k]].x = window->gridx[vs];
	  p[k][nk[k]].y = window->gridy[vs];
	  nk[k]++;
	  mask[v] = 1;
	}
      }
    }
    i_free_1d(mask);
  }

  /* Sediments                                                       */
  k = 0;
  nk[k] = nk[nz];
  for (n = 0; n < nk[k]; n++) {
    p[k][n].x = p[nz][n].x;
    p[k][n].y = p[nz][n].y;
  }
  for (c = 0; c < window->szcS; c++)
    window->c2p[0][c] = window->c2p[nz][c];

  /* Fill the points array and create the triangulation for every    */
  /* layer. Note this is designed to be used with COMPAS arrays,     */
  /* which start at index 1, hence use [c+1].                        */
  windat->d = (delaunay **)calloc(nz+1, sizeof(delaunay *));

  for (k = 0; k <= nz; k++) {
    windat->d[k] = NULL;
    if (nk[k] == 0) continue;
    windat->d[k] = delaunay_build(nk[k], p[k], 0, NULL, 0, NULL);
    /*printf("%d %d %f %f\n",k,windat->d[k]->n_point_triangles[2],windat->d[k]->points[2].x,windat->d[k]->points[2].y);*/

    /*
    if(k==0) {
      for (c = 0; c < windat->d[k]->ntriangles; c++) {
        triangle* t = &windat->d[k]->triangles[c];
	printf("%f %f\n",windat->d[k]->points[t->vids[0]].x, windat->d[k]->points[t->vids[0]].y);
	printf("%f %f\n",windat->d[k]->points[t->vids[1]].x, windat->d[k]->points[t->vids[1]].y);
	printf("%f %f\n",windat->d[k]->points[t->vids[2]].x, windat->d[k]->points[t->vids[2]].y);
	printf("%f %f\n",windat->d[k]->points[t->vids[0]].x, windat->d[k]->points[t->vids[0]].y);
	printf("NaN NaN\n");
      }

      printf("%d\n",windat->d[k]->ntriangles);
      for (c = 0; c < windat->d[k]->npoints; c++)
	printf("%f %f\n",windat->d[k]->points[c].x, windat->d[k]->points[c].y);
      cc = 0;
      for (c = 0; c < windat->d[k]->nedges; c++) {
	printf("%f %f\n",windat->d[k]->points[windat->d[k]->edges[cc]].x, windat->d[k]->points[windat->d[k]->edges[cc]].y);
	printf("%f %f\n",windat->d[k]->points[windat->d[k]->edges[cc+1]].x, windat->d[k]->points[windat->d[k]->edges[cc+1]].y);
	cc+=2;
	printf("NaN NaN\n");
      } 

    }
    */
    for (c = 0; c < nk[k]; c++) {
      windat->d[k]->points[c].v = d_alloc_1d(4);
    }
  }

  for (k = 0; k <= nz; k++) {
    if (windat->d[k] != NULL) {
      delaunay *d = windat->d[k];
      d->xmin = windat->d[nz]->xmin;
      d->xmax = windat->d[nz]->xmax;
      d->ymin = windat->d[nz]->ymin;
      d->ymax = windat->d[nz]->ymax;
      d->ptf = 0;
    }
  }

  /* Set the momentum interpolation scheme                           */
  wincon->mosl = wincon->osl;
  strcpy(wincon->momsr, wincon->trasr);

  if (newcode) {

  wincon->osl = L_LINEAR;
  strcpy(wincon->trasr, "linear");
    /*
  wincon->mosl = L_LSQUAD;
  strcpy(wincon->momsr, "quadratic");

  wincon->mosl = L_LSLIN;
  strcpy(wincon->momsr, "linearlsq");

  wincon->mosl = L_LINEAR;
  strcpy(wincon->momsr, "linear");

  wincon->mosl = L_SIB;
  strcpy(wincon->momsr, "nn_sibson");
    */
  /* Set up the point triangles                                      */
  windat->nptt = (npt_t **)calloc(nz+1, sizeof(npt_t *));
  windat->nptm = (npt_t **)calloc(nz+1, sizeof(npt_t *));

  /* Free the point_triangles array if necessary                     */
  for (k = 0; k <= nz; k++) {
    if (windat->d[k] != NULL) {
      windat->nptt[k] = npt_create(windat->d[k]->npoints);
      windat->nptm[k] = npt_create(windat->d[k]->npoints);
      if (windat->d[k]->point_triangles != NULL) {
	for (cc = 0; cc < windat->d[k]->npoints; ++cc)
	  if (windat->d[k]->point_triangles[cc] != NULL)
	    free(windat->d[k]->point_triangles[cc]);
	free(windat->d[k]->point_triangles);
      }
      if (windat->d[k]->n_point_triangles != NULL)
	free(windat->d[k]->n_point_triangles);
    }
  }

  /* Tracer point triangles                                          */
  if (wincon->osl & L_BILIN) {
    set_interp_square_n(window, windat, wincon, windat->d, windat->nptt);
  } else if (wincon->osl & L_BAYLIN) {
    set_interp_ring_n(window, windat, wincon, windat->d, windat->nptt);
  } else if (!(wincon->osl & (L_BILIN|L_BAYLIN))) {
    if (sf > 0)
      set_interp_points_n(window, windat, wincon, windat->d, windat->nptt, sf);
  }
  /* Momentum point triangles                                        */
  if (wincon->mosl & L_BILIN) {
    set_interp_square_n(window, windat, wincon, windat->d, windat->nptm);
  } else if (wincon->mosl & L_BAYLIN) {
    set_interp_ring_n(window, windat, wincon, windat->d, windat->nptm);
  } else if (!(wincon->mosl & (L_BILIN|L_BAYLIN))) {
    if (sf > 0)
      set_interp_points_n(window, windat, wincon, windat->d, windat->nptm, sf);
  }
  } else {
    /* Original code                                                   */
  if (wincon->osl & L_BILIN) {
    set_interp_square(window, windat, wincon, windat->d);
  } else if (wincon->osl & L_BAYLIN) {
    set_interp_ring(window, windat, wincon, windat->d);
  } else {
    if (sf > 0)
      set_interp_points(window, windat, wincon, windat->d, sf);
  }
  }

  set_tr_nograd(window, windat->temp);

  /* Set the data values, for u, v, w and tracers                    */
  memset(nk, 0, (nz+1) * sizeof(int));
  for (cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    k = window->s2k[c];
    windat->d[k+1]->points[nk[k+1]].v[0] = windat->u[c];
    windat->d[k+1]->points[nk[k+1]].v[1] = windat->v[c];
    windat->d[k+1]->points[nk[k+1]].v[2] = windat->w[c];
    windat->d[k+1]->points[nk[k+1]].v[3] = windat->temp[c];
    nk[k+1]++;
  }

  if (ef) {
    mask = i_alloc_1d(window->szv);
    memset(mask, 0, window->szv * sizeof(int));
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      k = window->s2k[c];
      for (n = 1; n <= window->npe[c2]; n++) {
	int v = window->c2v[n][c];
	if (v && !mask[v]) {
	  windat->d[k+1]->points[nk[k+1]].v[0] = vertex_weighted_tr(window, windat, wincon, windat->u, v);
	  windat->d[k+1]->points[nk[k+1]].v[1] = vertex_weighted_tr(window, windat, wincon, windat->v, v);
	  windat->d[k+1]->points[nk[k+1]].v[2] = vertex_weighted_tr(window, windat, wincon, windat->w, v);
	  windat->d[k+1]->points[nk[k+1]].v[3] = vertex_weighted_tr(window, windat, wincon, windat->temp, v);
	  nk[k+1]++;
	  mask[v] = 1;
	}
      }
    }
    i_free_1d(mask);
  }

  /* Sediments                                                       */
  set_trans_sed(window, 0, windat->u);
  set_trans_sed(window, 1, windat->v);
  set_trans_sed(window, 2, windat->w);
  set_trans_sed(window, 3, windat->temp);

  /*-----------------------------------------------------------------*/
  /* Set up the GRID_SPEC structure for each layer and each          */
  /* variable, where all structures use the same triangulation, d.   */
  grid_spec_init_tran(window, windat, wincon, 0);
  grid_spec_init_tran(window, windat, wincon, 1);
  grid_spec_init_tran(window, windat, wincon, 2);
  grid_spec_init_tran(window, windat, wincon, 3);
}

/* END tran_grid_init()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initializes a transport grid_spec structure                       */
/*-------------------------------------------------------------------*/
void grid_spec_init_tran(geometry_t *window, 
			 window_t *windat, 
			 win_priv_t *wincon, 
			 int vid)
{
  char *uvrule = "linear";
  char rule[MAXSTRLEN];
  int k, kk, nz = window->nz;
  GRID_SPECS **gs;

  strcpy(rule, uvrule);
  /*if (wincon->osl & (L_BILIN|L_LSQUAD|L_LSLIN)) strcpy(rule, wincon->trasr);*/
  if (wincon->osl & (L_BILIN|L_BAYLIN)) strcpy(rule, wincon->trasr);

  /*
  if (vid == 3 && wincon->osl & L_BILIN) strcpy(rule, "baylinear");
  wincon->mosl = (wincon->osl & L_BILIN) ? 1 : 0;
  wincon->nosl = (strcmp(uvrule, "bilinear") == 0) ? 1 : 0;
  */

  /*-----------------------------------------------------------------*/
  /* Allocate and initialize                                         */
  gs = (GRID_SPECS **)calloc(nz+1, sizeof(GRID_SPECS *));
  if (vid == 3)
    reset_npt(window, windat->d, windat->nptt);
  else
    reset_npt(window, windat->d, windat->nptm);

  for (k = 0; k <= nz; k++) {

    gs[k] = grid_spec_create();

    if (windat->d[k] == NULL) continue;
    /*grid_spec_init(gs[k]);*/

    if (vid == 3) {
      grid_interp_init_t(gs[k], windat->d[k], wincon->trasr, vid);
    } else {
      grid_interp_init_t(gs[k], windat->d[k], wincon->momsr, vid);
    }

    /* Set the depth levels */
    gs[k]->nz = nz;
    gs[k]->z = d_alloc_1d(nz);
    for (kk = 0; kk < nz; kk++) {
      gs[k]->z[kk] = 0.5 * (window->layers[kk] + window->layers[kk+1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the correct grid_spec structure                             */
  if (vid == 0)
    windat->gsx = gs;
  else if (vid == 1)
    windat->gsy = gs;
  else if (vid == 2)
    windat->gsz = gs;
  else if (vid == 3)
    windat->gst = gs;
}

/* END grid_spec_init_tran()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets the weights for interpolation                              */
/*-------------------------------------------------------------------*/
void delaunay_reinit(geometry_t *window, delaunay **d, int vid, double *v)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c, cc, c2, k, n, m;
  GRID_SPECS **gs;
  int nk[window->nz+1];
  int ef = 0;
  int *mask;
  int id;
  int osl;

  /*-----------------------------------------------------------------*/
  /* Retrieve the correct grid_spec structure                        */
  osl = (vid == 3) ? wincon->osl : wincon->mosl;
  /*osl = wincon->osl;*/
  if (vid == 0)
    gs = windat->gsx;
  else if (vid == 1)
    gs = windat->gsy;
  else if (vid == 2)
    gs = windat->gsz;
  else if (vid == 3) {
    gs = windat->gst;
    set_tr_nograd(window, v);
  }

  /*
  if (vid == 3)
    reset_npt(window, d, windat->nptt);
  else
    reset_npt(window, d, windat->nptm);
  */
  /*if ((vid == 3 && wincon->osl) || (vid < 3 && wincon->mosl)) {*/
  if (osl & (L_BILIN|L_BAYLIN)) {
    /* Set the structures for quad mesh bilinear interpolations      */
    double *cx = wincon->w6;
    double *cy = wincon->w7;
    int *cl = wincon->s2;

    /* Initialise. Only wet points are filled below, but set all     */
    /* points (including ghosts).                                    */
    for (k = 0; k <= window->nz; k++) {
      if (d[k] != NULL) {
	for (cc = 0; cc < d[k]->npoints; cc++)
	  d[k]->points[cc].v[vid] = 0.0;
      }
    }
    /* Fill the points with source and destination information       */
    for (cc = 1; cc <= wincon->vc; cc++) {
      int ks, cs;
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      k = window->s2k[c] + 1;
      id = window->c2p[k][c2];
      d[k]->points[id].x = cx[c];
      d[k]->points[id].y = cy[c];

      cs = cl[cc];
      ks = window->s2k[cs] + 1;
      if (cs == window->zm1[cs]) ks = 0;
      c2 = window->m2d[cs];
      if (window->wgst[c2]) c2 = window->wgst[c2];
      d[k]->points[id].z = (double)window->c2p[ks][c2];
      d[k]->points[id].v[vid] = (double)ks;
      /*if (window->wgst[cs]) d[k]->points[id].v[vid] = (double)window->s2k[window->wgst[cs]] + 1;*/

      /*if (vid==3&&c==1443)printf("pos %d %d %d %f %f : %f %f %d\n",c,k,id,window->cellx[c2],window->celly[c2],d[k]->points[id].x,d[k]->points[id].y,id);*/
    }
  } else {

    /*---------------------------------------------------------------*/
    /* Fill the points array in the delaunay structure. Note this is */
    /* designed to be used with COMPAS arrays, which start at index  */
    /* 1, hence use c2-1.                                            */
    memset(nk, 0, (window->nz+1) * sizeof(int));
    for (cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      k = window->s2k[c] + 1;
      c2 = window->m2d[c];
      id = window->c2p[k][c2];
      d[k]->points[id].z = v[c];
    }
    /* OBC ghosts                                                    */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	c2 = open->obc_e2[ee];
	k = window->s2k[c2] + 1;
	for (m = 0; m < open->bgz; m++) {
	  c2 = window->m2d[c];
	  id = window->c2p[k][c2];
	  if (id != -1) d[k]->points[id].z = v[c];
	  c = open->omape[ee][c];
	}
      }
    }
  }
  set_trans_sed(window, vid, v);

  if (ef) {
    int n;
    mask = wincon->s4;
    memset(mask, 0, window->szv * sizeof(int));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      k = window->s2k[c] + 1;
      for (n = 1; n <= window->npe[c2]; n++) {
	int vi = window->c2v[n][c];
	if (vi && !mask[vi]) {
	  d[k]->points[nk[k]].z = vertex_weighted_tr(window, windat, wincon, v, vi);
	  /*if(v==windat->sal && c==534 && d[k]->points[nk[k]].z!=35.0)printf("%f %d %d %f : %f %d\n",windat->days,c,vi,d[k]->points[nk[k]].z,v[c],window->c2p[k][c2]);*/
	  nk[k]++;
	  mask[vi] = 1;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Rebuild the weights                                             */
  if (vid <= 3) {

    for (k = 0; k <= window->nz; k++) {
      if (d[k] == NULL) continue;
      if (osl & (L_LINEAR|L_CUBIC|L_LSQUAD|L_LSLIN)) {
	gs[k]->rebuild(gs[k]->interpolator, d[k]->points);
      } else if (osl & (L_BILIN|L_BAYLIN)) {
	windat->d[k]->vid = vid;
	if (k == 0) continue;
	for (cc = 0; cc < d[k]->npoints; cc++) {
	  int ks = (k == 0) ? 0 : (int)d[k]->points[cc].v[vid];
	  d[k]->points[cc].v[vid] = (double)cc;
	  gs[k]->rebuild2(gs[k]->interpolator, gs[ks]->interpolator, &d[k]->points[cc]);
	}
      } else {
	/* Not properly function: uses ht_delete() in nnpi_interpolate()
	if (osl & (L_SIB|L_NONSIB))
	  nnhpi_destroy_weights(gs[k]->interpolator);
	*/
	for (cc = 0; cc < d[k]->npoints; cc++)
	  gs[k]->rebuild(gs[k]->interpolator, &d[k]->points[cc]);
      }
    }
  } else {
    for (k = 0; k <= window->nz; k++) {
      if (d[k] == NULL) continue;
      gs[k]->rebuild(gs[k]->interpolator, d[k]->points);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reset values for quad mesh bilinear interpolations              */
  /*if ((vid == 3 && wincon->mosl) || (vid < 3 && wincon->nosl)) {*/
  if(osl & (L_BILIN|L_BAYLIN)) {
    for (cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      k = window->s2k[c] + 1;
      c2 = window->m2d[c];
      id = window->c2p[k][c2];
      d[k]->points[id].v[vid] = v[c];
    }
    /* OBC ghosts                                                    */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	c2 = open->obc_e2[ee];
	k = window->s2k[c2] + 1;
	for (m = 0; m < open->bgz; m++) {
	  c2 = window->m2d[c];
	  id = window->c2p[k][c2];
	  if (id != -1) d[k]->points[id].v[vid] = v[c];
	  c = open->omape[ee][c];
	}
      }
    }
  }
}

/* END delaunay_reinit()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets the point triangles in the delaunay structure with tracer  */
/* or momentum point trianges. The point trianges are lists of       */
/* points in the delaunay triangulation used in interpolations       */
/* around any given triangulation point, and are dependent on the    */
/* type of interpolation used.                                       */
/*-------------------------------------------------------------------*/
void reset_npt(geometry_t *window, delaunay **d, npt_t **npt) 
{
  int k;

  for (k = 0; k <= window->nz; k++) {
    if (d[k] != NULL) {
      d[k]->n_point_triangles = npt[k]->npt;
      /*memcpy(d[k]->n_point_triangles, npt[k]->npt, d[k]->npoints * sizeof(int));*/
      d[k]->point_triangles = npt[k]->pt;
      /*
      for (i = 0; i < d[k]->npoints; i++)
        memcpy(d[k]->point_triangles[i], npt[k]->pt[i], npt[k]->npt[i] * sizeof(int));
      */
    }
  }
}

/* END reset_npt()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Creates a point triangle structure                                */
/*-------------------------------------------------------------------*/
npt_t* npt_create(int n)
{
  npt_t* npt = malloc(sizeof(npt_t));

  npt->npt = calloc(n, sizeof(int));
  memset(npt->npt, 0, n * sizeof(int));
  npt->pt = malloc(n * sizeof(int*));

  return npt;
}

/* END npt_create()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the boundary condition for the sediment layer                */
/*-------------------------------------------------------------------*/
void set_trans_sed(geometry_t *window, int vid, double *v)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int ee, cc, c, c2, n, id;
  double s = -1.0;
  double *var = wincon->d6;
  int osl;

  if (vid == 3) s = 1.0;
  osl = (vid == 3) ? wincon->osl : wincon->mosl;
  /*osl = wincon->osl;*/

  /* Set the sediment cells beneath wet cells                        */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    c2 = window->m2d[c];
    var[c2] = s * v[c];
  }

  /* Set the ghost cells adjacent to wet cells                       */
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    c2 = window->bin[cc];
    var[c] = var[c2];
  }
  /* OBC ghosts                                                      */
  /*
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      c2 = open->obc_e2[ee];
      v[c] = v[c2];
    }
  }
  */
  /* Fill the Delaunay values                                        */
  for (cc = 1; cc <= window->n2_t; cc++) {
    c = window->w2_t[cc];
    id = window->c2p[0][c];
    windat->d[0]->points[id].v[vid] = var[c];
    windat->d[0]->points[id].z = var[c];
  }

  if (osl & (L_BILIN|L_BAYLIN)) {
    int *cl = wincon->s2;
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      id = window->c2p[0][c2];
      c2 = window->m2d[cl[cc]];
      if (window->wgst[c2]) c2 = window->wgst[c2];
      windat->d[0]->points[id].z = (double)window->c2p[0][c2];
    }
  }
}

/* END set_trans_sed()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a no-gradient condition over all ghost cells for a tracer    */
/*-------------------------------------------------------------------*/
void set_tr_nograd(geometry_t *window, double *tr)
{
  int n, ee, cc, c, c2;

  /* Lateral ghost cells                                             */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    c2 = window->bin[cc];
    tr[c] = tr[c2];
  }

  /* Sediment cells                                                  */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    tr[window->zm1[c]] = tr[c];
  }
  /* OBC ghosts                                                      */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      c2 = c = open->obc_e2[ee];
      do {
	c = open->omape[ee][c];
	tr[c] = tr[c2];
      } while (c != open->omape[ee][c]);
    }
    OBC_bgz_nogradb(window, open, tr);
  }
}

/* END set_tr_nograd()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates in 3D on an unstructured grid given an (x,y,z)       */
/* location. The cell where (x,y,z) resides must also be supplied,   */
/* computed with hd_grid_xyztoc() outside this routine since this    */
/* routine may be called many times for different variables.         */
/*-------------------------------------------------------------------*/
double hd_trans_interp(geometry_t *window, GRID_SPECS **gs, double x, double y, double z, int c, int co, int vid)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int osl = (vid == 3) ? wincon->osl : wincon->mosl;
  int k1, k2, kd;
  int cs = window->m2d[c];
  double v, v1, v2, d;
  int zp1 = window->zp1[c];
  double dzz = wincon->dzz[c];        /* cellz thickness             */
  double *cellz = wincon->cellz;      /* cellz depth                 */
  int is_sed = (c == window->zm1[c]) ? 1 : 0;
  double bfact = (osl & L_FG) ? 0.0 : -1.0;

  /*
  if (window->wgst[c]) {
    c = window->wgst[c];
    cs = window->m2d[c];
    zp1 = window->zp1[c];
    dzz = wincon->dzz[c];
  }
  */
  /* Evaluate in the layer with a cellz level greater than cz. If    */
  /* this is the surface layer then self mappings or a no-gradient   */
  /* above the surface will take care of the vertical interpolation. */
  k2 = window->s2k[zp1] + 1;

  /* For least squares, the interpolation is centered on the cell    */
  /* centre (i.e. co). This is mapped to the Delaunay index and      */
  /* passed to the interpolation scheme directly in gs->id (i.e.     */
  /* there is no need to find the triangle (x,y) resides in within   */
  /* the interpolation scheme via delaunay_xytoi_ng().               */
  if (osl & (L_LSQUAD|L_LSLIN)) {
    gs[k2]->id = window->c2p[k2][cs];
  }

  /* For bilinear and baylinear interpolation, supply the source and */
  /* destination indicies.                                           */
  if (osl & (L_BILIN|L_BAYLIN)) {
    kd = window->s2k[co] + 1;
    if (window->wgst[c]) {
      c = window->wgst[c];
      cs = window->m2d[c];
    }

    x = (double)window->c2p[kd][window->m2d[co]];
    y = (double)window->c2p[k2][cs];
    windat->d[k2]->vid = vid;
    v2 = grid_interp_on_point2(gs, kd, k2, x, y);
  } else {
    v2 = grid_interp_on_point(gs[k2], x, y);
    /*if(vid==3&&co==2031)printf("interp1 %d %f : %d %d %d %f %f\n",k2-1,v2,window->s2i[c],window->s2j[c],window->s2k[c],x,y);*/
  }
  /* If the streamline is in the sediment, use the layer above       */
  /*
  */

  if(isnan(v2)){
    v2 = gs[k2]->d->points[window->c2p[k2][cs]].v[vid];
    hd_quit("hd_tran_interp: NaN c=%d(k%d) cl=%d(k%d) k=%d %d %f %f %f %d\n",co,window->s2k[co],zp1,window->s2k[zp1],k2, vid,x,y,v2,window->wgst[c]);
    /*
    print_tri_k(gs[k2]->d, k2);
    */

  }

  /* Evaluate in the layer with a cellz level less than cz. If this  */
  /* is the bottom layer, then the ghost layer below the bottom will */
  /* take care of the vertical interpolation.                        */
  if (is_sed) {
    k1 = k2;
    v1 = v2;
    if (vid < 3) {
      v1 *= -1.0;
      /*dzz *= 2.0;*/
    }
  } else {
    k1 = window->s2k[c] + 1;
    if (osl & (L_LSQUAD|L_LSLIN)) gs[k1]->id = (double)window->c2p[k1][cs];
    if (osl & (L_BILIN|L_BAYLIN)) {
      x = (double)window->c2p[kd][window->m2d[co]];
      y = (double)window->c2p[k1][cs];
      windat->d[k1]->vid = vid;
      v1 = grid_interp_on_point2(gs, kd, k1, x, y);
    } else {
      v1 = grid_interp_on_point(gs[k1], x, y);
      /*if(vid==3&&c==3743)printf("interp2 %d %f %d\n",k1-1,v1,gs[k1]->id);*/
    }
    /*if(isnan(v1)) v1 = gs[k1]->d->points[window->c2p[k1][cs]].v[vid];*/
  }

  /* Interpolate vertically                                          */
  /* Linear scheme */

  d = z - cellz[c];
  v = d * (v2 - v1) / dzz + v1;

  /*
  if (v > 35.0 && window->s2i[c]!=1) {
    printf("%d %d %d %f\n",window->s2i[c],window->s2j[c],window->s2k[c], v);
    exit(0);
  }
  */
  /* Upwind scheme */
  /*
  if (z < window->gridz[c])
    return(v1);
  else
    return(v2);
  */
  /*if (vid==3&&c==575)printf("%f c=[%d %d] k=[%d %d] z=%f cellz=%f %f d=%f v1=%f v2=%f : %f %f %f\n",windat->days,c,zp1,k1,k2,z,cellz[c],cellz[window->zp1[c]],d,v1,v2,v,dzz,windat->eta[window->m2d[c]]);*/


  return(v);

}

/* END hd_trans_interp()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Same as hd_trans_interp() but uses a higher order interpolation   */
/* in the vertical.                                                  */
/*-------------------------------------------------------------------*/
double hd_trans_interpo(geometry_t *window, GRID_SPECS **gs, double x, double y, double z, int c, int co, int vid)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int k1, k2;
  double v, v1, v2, d;
  int zp1 = window->zp1[c];
  int zm1 = window->zm1[c];
  double dzz = wincon->dzz[c];        /* cellz thickness             */
  double *cellz = wincon->cellz;      /* cellz depth                 */
  int is_sed = (c == window->zm1[c]) ? 1 : 0;
  int osl = 1;
  int n;
  double wgt[5], va[5];
  double mn = 1e10, mx = -1e10;

  if(vid==3)osl = 0;

  /* Evaluate in the layer with a cellz level greater than cz. If    */
  /* this is the surface layer then self mappings or a no-gradient   */
  /* above the surface will take care of the vertical interpolation. */
  k2 = window->s2k[zp1] + 1;
  v2 = grid_interp_on_point(gs[k2], x, y);

  /* Evaluate in the layer with a cellz level less than cz. If this  */
  /* is the bottom layer, then the ghost layer below the bottom will */
  /* take care of the vertical interpolation.                        */
  if (is_sed) {
    k1 = k2;
    v1 = v2;
    if (vid < 3) {
      v1 *= -1.0;
      /*dzz *= 2.0;*/
    }
  } else {
    k1 = window->s2k[c] + 1;
    v1 = grid_interp_on_point(gs[k1], x, y);
  }
  /* Interpolate vertically                                          */
  d = z - cellz[c];
  weights_v(window, c, d, wincon->dzz, wgt, osl);

  if (osl == 1) {
    va[1] = v2;
    va[0] = v1;
  } else {
    va[2] = v2;
    va[1] = v1;
    if (is_sed)
      va[0] = v1;
    else {
      if (zm1 == window->zm1[zm1]) {
	va[0] = v1;
      } else {
	k2 = window->s2k[zm1] + 1;
	va[0] = grid_interp_on_point(gs[k2], x, y);
      }
    }
    zp1 = window->zp1[zp1];
    k2 = window->s2k[zp1] + 1;
    va[3] = grid_interp_on_point(gs[k2], x, y);
  }
  v = 0.0;
  for(n = 0; n <= osl; n++) {
    v += (wgt[n] * va[n]);
    mx = max(va[n], mx);
    mn = min(va[n], mn);
  }
  v = max(min(v, mx), mn);

  /*
  if (co==window->zm1[5])printf("%f c=[%d %d] k=[%d %d] z=%f cellz=%f %f d=%f v1=%f v2=%f : %f : %f %f\n",windat->days,c,zp1,k1,k2,z,cellz[c],cellz[window->zp1[c]],d,va[0],va[1],v,wgt[0],wgt[1]);

  if(co==window->zm1[5]) {
    printf("%f %f %f %f\n",wgt[0],wgt[1],wgt[2],wgt[3]);
    printf("%f %f %f %f\n",va[0],va[1],va[2],va[3]);
  }
  */
  /*
  if(co==7) {
    printf("%f %f\n",wgt[0],wgt[1]);
    printf("%f %f\n",va[0],va[1]);
  }
  */
  return(v);

}

/* END hd_trans_interpo()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets a set of points to use in the least squares fit, and saves   */
/* the Delaunay indicies of these points in the Delaunay structure;  */
/* d->point_triangles. The flag d->ptf is set to 1 to indicate this  */
/* has been done.                                                    */
/*-------------------------------------------------------------------*/
void set_interp_points(geometry_t *window, 
		       window_t *windat, 
		       win_priv_t *wincon,
		       delaunay **d,
		       int size
		       )
{
  int cc, c, c2, k, i, j, n, m, id;
  int cn, cns, cg, cgs;
  int *st = NULL, sz, np;
  int nz = window->nz;
  int *cells, *mask;

  /*-----------------------------------------------------------------*/
  /* Free the point_triangles array if necessary                     */
  for (k = 0; k <= nz; k++) {
    if (d[k] != NULL) {
      if (d[k]->point_triangles != NULL) {
	for (i = 0; i < d[k]->npoints; ++i)
	  if (d[k]->point_triangles[i] != NULL)
	    free(d[k]->point_triangles[i]);
      }
      memset(d[k]->n_point_triangles, 0, d[k]->npoints * sizeof(int));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Set the flag in the Delaunay structure                        */
    d[k]->ptf = 1;
    /* Get the stencil                                               */
    sz = size;
    st = stencil(window, c, &sz, ST_SIZED, 0);
    /* Wet cells to include                                          */
    d[k]->n_point_triangles[id] = sz;
    /* Ghost cells to include                                        */
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      for(n = 1; n <= window->npe[cns]; n++) {
	cg = window->c2c[n][cn];
	cgs = window->m2d[cg];
	if (window->wgst[cg] && window->c2p[k][cgs] >= 0)
	  d[k]->n_point_triangles[id]++;
      }
    }

    /* Allocate                                                      */
    d[k]->point_triangles[id] = malloc(d[k]->n_point_triangles[id] * sizeof(int));
    cells = i_alloc_1d(d[k]->n_point_triangles[id]);

    /* Fill the locations of Delaunay indices                        */
    /* Wet cells to include                                          */
    j = 0;
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      cells[j] = cn;
      d[k]->point_triangles[id][j++] = window->c2p[k][cns];
      if(window->c2p[k][cns] == -1) 
	hd_quit("Invalid wet c2p map at c=%d: s[%d]=%d\n",c, i, cn);
    }
    /* Ghost cells to include                                        */
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      for(n = 1; n <= window->npe[cns]; n++) {
	cg = window->c2c[n][cn];
	cgs = window->m2d[cg];
	if (window->wgst[cg] && window->c2p[k][cgs] >= 0) {
	  cells[j] = cg;
	  d[k]->point_triangles[id][j++] = window->c2p[k][cgs];
	  if(window->c2p[k][cns] == -1) 
	    hd_quit("Invalid ghost c2p map at c=%d: s[%d]=%d\n",c, i, cn);
	}
      }
    }
    i_free_1d(st);
  }

  /* Set the sediments                                               */
  d[0]->ptf = 1;
  d[0]->n_point_triangles = d[nz]->n_point_triangles;
  d[0]->point_triangles = d[nz]->point_triangles;
}

void set_interp_points_n(geometry_t *window, 
			 window_t *windat, 
			 win_priv_t *wincon,
			 delaunay **d,
			 npt_t **npt,
			 int size
			 )
{
  int cc, c, c2, k, i, j, n, m, id;
  int cn, cns, cg, cgs;
  int *st = NULL, sz, np;
  int nz = window->nz;
  int *cells, *mask;

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Set the flag in the Delaunay structure                        */
    d[k]->ptf = 1;
    /* Get the stencil                                               */
    sz = size;
    st = stencil(window, c, &sz, ST_SIZED, 0);

    /* Wet cells to include                                          */
    npt[k]->npt[id] = sz;
    /* Ghost cells to include                                        */
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      for(n = 1; n <= window->npe[cns]; n++) {
	cg = window->c2c[n][cn];
	cgs = window->m2d[cg];
	if (window->wgst[cg] && window->c2p[k][cgs] >= 0)
	  npt[k]->npt[id]++;
      }
    }

    /* Allocate                                                      */
    npt[k]->pt[id] = malloc(npt[k]->npt[id] * sizeof(int));
    cells = i_alloc_1d(npt[k]->npt[id]);

    /* Fill the locations of Delaunay indices                        */
    /* Wet cells to include                                          */
    j = 0;
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      cells[j] = cn;
      npt[k]->pt[id][j++] = window->c2p[k][cns];
      if(window->c2p[k][cns] == -1) 
	hd_quit("Invalid wet c2p map at c=%d k=%d: s[%d]=%d cns=%d\n",c, k, i, cn, cns);
    }
    /* Ghost cells to include                                        */
    for (i = 0; i < sz; i++) {
      cn = st[i];
      cns = window->m2d[cn];
      for(n = 1; n <= window->npe[cns]; n++) {
	cg = window->c2c[n][cn];
	cgs = window->m2d[cg];
	if (window->wgst[cg] && window->c2p[k][cgs] >= 0) {
	  cells[j] = cg;
	  npt[k]->pt[id][j++] = window->c2p[k][cgs];
	  if(window->c2p[k][cns] == -1) 
	    hd_quit("Invalid ghost c2p map at c=%d: s[%d]=%d\n",c, i, cn);
	}
      }
    }
    i_free_1d(st);
  }

  /* OBC ghosts                                                      */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    int ee;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      cn = open->obc_e2[ee];
      for (m = 0; m < open->bgz; m++) {
	c2 = window->m2d[c];

	k = window->s2k[cn] + 1;
	id = window->c2p[k][c2];

	sz = size;
	st = stencil(window, c, &sz, ST_SIZED, 0);
	/* Wet cells to include                                      */
	npt[k]->npt[id] = sz;

	/* Allocate                                                  */
	npt[k]->pt[id] = malloc(npt[k]->npt[id] * sizeof(int));
	cells = i_alloc_1d(npt[k]->npt[id]);

	/* Fill the locations of Delaunay indices                    */
	/* Wet cells to include                                      */
	j = 0;
	for (i = 0; i < sz; i++) {
	  cn = st[i];
	  cns = window->m2d[cn];
	  cells[j] = cn;
	  npt[k]->pt[id][j++] = window->c2p[k][cns];
	  if(window->c2p[k][cns] == -1) 
	    hd_quit("Invalid wet c2p map at c=%d: s[%d]=%d\n",c, i, cn);
	}
	i_free_1d(st);
	c = open->omape[ee][c];
      }
    }
  }

  /* Set the sediments                                               */
  d[0]->ptf = 1;
  npt[0]->npt = npt[nz]->npt;
  npt[0]->pt = npt[nz]->pt;
}

/* END set_interp_points()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to generate a stencil for quad grids. The indexing of the */
/* stencil in the Delaunay point_triangles list is as follows:       */
/* index 0 = central cell                                            */
/* index 1 = west cell                                               */
/* index 2 = north cell                                              */
/* index 3 = east cell                                               */
/* index 4 = south cell                                              */
/* index 5 = NW cell                                                 */
/* index 6 = NE cell                                                 */
/* index 7 = SE cell                                                 */
/* index 8 = SW cell                                                 */
/* Given a streamline falls within cell c, the quadrant and          */
/* indices to use are:                                               */
/* dx = x - cellx[c], dy = y - celly[c]                              */
/* dx < 0, dy > 0: quadrant 1, indices [0 1 5 2]                     */
/* dx > 0, dy > 0: quadrant 2, indices [0 2 6 3]                     */
/* dx > 0, dy < 0: quadrant 3, indices [0 3 8 4]                     */
/* dx < 0, dy < 0: quadrant 4, indices [0 4 2 1]                     */
/* Any index with value -1 gets the central index value.             */
/*-------------------------------------------------------------------*/
void set_interp_square(geometry_t *window, 
		       window_t *windat, 
		       win_priv_t *wincon,
		       delaunay **d
		       )
{
  int cc, c, c2, k, i, j, jj, n, id, idn;
  int cn, cns, cg, cgs;
  int *st = NULL, sz, np;
  int nz = window->nz;

  /*-----------------------------------------------------------------*/
  /* Free the point_triangles array if necessary                     */
  for (k = 0; k <= nz; k++) {
    if (d[k] != NULL) {
      if (d[k]->point_triangles != NULL) {
	for (i = 0; i < d[k]->npoints; ++i)
	  if (d[k]->point_triangles[i] != NULL)
	    free(d[k]->point_triangles[i]);
      }

      if (k > 0) {
	d[k]->ptf = 1;
	for (i = 0; i < d[k]->npoints; i++) {
	  d[k]->n_point_triangles[i] = 9;
	  d[k]->point_triangles[i] = malloc(d[k]->n_point_triangles[i] * sizeof(int));
	  for (j = 0; j < d[k]->n_point_triangles[i]; j++)
	    d[k]->point_triangles[i][j] = -1;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];
    /* Set the flag in the Delaunay structure                        */
    jj = 0;
    d[k]->point_triangles[id][jj++] = id;

    for (j = 1; j <= 4; j++) {
      cn = window->c2c[j][c];
      cns = window->m2d[cn];
      idn = window->c2p[k][cns];
      d[k]->point_triangles[id][jj++] = idn;
    }
    for (j = 1; j <= 4; j++) {
      int c1 = window->c2c[j][c];
      n = (j == 4) ? 1 : j + 1;
      c2 = window->c2c[n][c];
      if (!window->wgst[c1]) {
	cgs = window->c2c[n][window->m2d[c1]];
	idn = window->c2p[k][cgs];
	d[k]->point_triangles[id][jj++] = idn;
      } else if (!window->wgst[c2]) {
	cgs = window->c2c[j][window->m2d[c2]];
	idn = window->c2p[k][cgs];
	d[k]->point_triangles[id][jj++] = idn;
      } else
	d[k]->point_triangles[id][jj++] = -1;
    }
  }

  /* Ghost cells. Theoretically the streamline should never end up   */
  /* in a ghost cell, but we include a stencil for interpolation     */
  /* around ghost cells for safety.                                  */
  /*
  for(cc = window->b3_t + 1; cc <= window->n3_t; cc++) {
    int j1;
    cg = window->w3_t[cc];
    cgs = window->m2d[cg];
    c = window->wgst[cg];
    k = window->s2k[c] + 1;
    id = window->c2p[k][cgs];
    for (j = 1; j <= window->npe[window->m2d[c]]; j++) {
      if(window->c2c[j][c] == cg) break;
    }
    jj = 0;
    d[k]->point_triangles[id][jj++] = id;
    j1 = (j == 4) ? 1 : j + 1;
    c2 = window->c2c[j1][c];
    d[k]->point_triangles[id][jj++] = window->c2p[k][c2];
    if (!window->wgst[c2]) {
      c2 = window->c2c[j][c2];
      d[k]->point_triangles[id][jj++] = window->c2p[k][c2];
    }
    j1 = (j == 1) ? 4 : j - 1;
    c2 = window->c2c[j1][c];
    d[k]->point_triangles[id][jj++] = window->c2p[k][c2];
    if (!window->wgst[c2]) {
      c2 = window->c2c[j][c2];
      d[k]->point_triangles[id][jj++] = window->c2p[k][c2];
    }
    npt[k]->npt[i] = jj;
  }
  */
  /* Set the sediments                                               */
  d[0]->ptf = 1;
  d[0]->n_point_triangles = d[nz]->n_point_triangles;
  d[0]->point_triangles = d[nz]->point_triangles;
}


void set_interp_square_n(geometry_t *window, 
			 window_t *windat, 
			 win_priv_t *wincon,
			 delaunay **d,
			 npt_t **npt
			 )
{
  int cc, c, c2, k, i, j, jj, n, id, idn;
  int cn, cns, cg, cgs;
  int *st = NULL, sz, np;
  int nz = window->nz;

  /*-----------------------------------------------------------------*/
  /* Allocate the point_triangles array                              */
  for (k = 0; k <= nz; k++) {
    if (d[k] != NULL) {
      if (k > 0) {
	d[k]->ptf = 1;
	for (i = 0; i < d[k]->npoints; i++) {
	  npt[k]->npt[i] = 9;
	  npt[k]->pt[i] = malloc(npt[k]->npt[i] * sizeof(int));
	  for (j = 0; j < npt[k]->npt[i]; j++)
	    npt[k]->pt[i][j] = -1;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];
    /* Set the flag in the Delaunay structure                        */
    jj = 0;
    npt[k]->pt[id][jj++] = id;

    for (j = 1; j <= 4; j++) {
      cn = window->c2c[j][c];
      cns = window->m2d[cn];
      idn = window->c2p[k][cns];
      npt[k]->pt[id][jj++] = idn;
    }
    for (j = 1; j <= 4; j++) {
      int c1 = window->c2c[j][c];
      n = (j == 4) ? 1 : j + 1;
      c2 = window->c2c[n][c];
      if (!window->wgst[c1]) {
	cgs = window->c2c[n][window->m2d[c1]];
	idn = window->c2p[k][cgs];
	npt[k]->pt[id][jj++] = idn;
      } else if (!window->wgst[c2]) {
	cgs = window->c2c[j][window->m2d[c2]];
	idn = window->c2p[k][cgs];
	npt[k]->pt[id][jj++] = idn;
      } else {
	npt[k]->pt[id][jj++] = -1;
      }
    }
  }

  /* Ghost cells. Theoretically the streamline should never end up   */
  /* in a ghost cell, but we include a stencil for interpolation     */
  /* around ghost cells for safety.                                  */
  /*
  for(cc = window->b3_t + 1; cc <= window->n3_t; cc++) {
    int j1;
    cg = window->w3_t[cc];
    cgs = window->m2d[cg];
    c = window->wgst[cg];
    k = window->s2k[c] + 1;
    id = window->c2p[k][cgs];
    for (j = 1; j <= window->npe[window->m2d[c]]; j++) {
      if(window->c2c[j][c] == cg) break;
    }
    jj = 0;
    npt[k]->pt[id][jj++] = id;
    j1 = (j == 4) ? 1 : j + 1;
    c2 = window->c2c[j1][c];
    npt[k]->pt[id][jj++] = window->c2p[k][c2];
    if (!window->wgst[c2]) {
      c2 = window->c2c[j][c2];
      npt[k]->pt[id][jj++] = window->c2p[k][c2];
    }
    j1 = (j == 1) ? 4 : j - 1;
    c2 = window->c2c[j1][c];
    npt[k]->pt[id][jj++] = window->c2p[k][c2];
    if (!window->wgst[c2]) {
      c2 = window->c2c[j][c2];
      npt[k]->pt[id][jj++] = window->c2p[k][c2];
    }
    d[k]->n_point_triangles[i] = jj;
  }
  */
  /* Set the sediments                                               */
  d[0]->ptf = 1;
  npt[0]->npt = npt[nz]->npt;
  npt[0]->pt = npt[nz]->pt;
}

/* END set_interp_square()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets a set of points to use in baycentric linear interpolation,   */
/* and saves the Delaunay indicies of these points in the Delaunay   */
/* structure; d->point_trianges. The flag d->ptf is set to 1 to      */
/* indicate this has been done. Indicies are consecutive around the  */
/* centre.                                                           */
/*-------------------------------------------------------------------*/
void set_interp_ring(geometry_t *window, 
		     window_t *windat, 
		     win_priv_t *wincon,
		     delaunay **d
		     )
{
  int cc, c, c2, k, i, j, n, id;
  int cn, cns;
  int nz = window->nz;

  /*-----------------------------------------------------------------*/
  /* Free the point_triangles array if necessary                     */
  for (k = 0; k <= nz; k++) {
    if (d[k] != NULL) {
      if (d[k]->point_triangles != NULL) {
	for (i = 0; i < d[k]->npoints; ++i)
	  if (d[k]->point_triangles[i] != NULL)
	    free(d[k]->point_triangles[i]);
      }
      memset(d[k]->n_point_triangles, 0, d[k]->npoints * sizeof(int));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate and initialize                                         */
  for(cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Set the flag in the Delaunay structure                        */
    d[k]->ptf = 1;

    d[k]->n_point_triangles[id] = window->npe[c2] + 1;
    d[k]->point_triangles[id] = malloc(d[k]->n_point_triangles[id] * sizeof(int));
    d[k]->point_triangles[id][0] = -1;
  }
  wincon->nlmap = i_alloc_1d(window->szc);
  wincon->lmap = i_alloc_2d(window->npem+1, window->szc);

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Fill the locations of Delaunay indices. The indices are       */
    /* consecutive around the centre.                                */
    d[k]->point_triangles[id][0] = window->c2p[k][c2];
    wincon->nlmap[c] = window->npe[c2] + 1;
    wincon->lmap[c][0] = c;
    for (j = 1; j <= window->npe[c2]; j++) {
      cn = window->c2c[j][c];
      cns = window->m2d[cn];
      d[k]->point_triangles[id][j] = window->c2p[k][cns];
      wincon->lmap[c][j] = cn;
      if(window->c2p[k][cns] == -1) 
	hd_quit("Invalid wet c2p map at c=%d: j=%d cn=%d %d %d\n",c, j, cn, window->s2i[c], window->s2j[c]);
    }
  }

  /* Set the sediments                                               */
  d[0]->ptf = 1;
  d[0]->n_point_triangles = d[nz]->n_point_triangles;
  d[0]->point_triangles = d[nz]->point_triangles;
}

void set_interp_ring_n(geometry_t *window, 
		       window_t *windat, 
		       win_priv_t *wincon,
		       delaunay **d,
		       npt_t **npt
		       )
{
  int cc, c, c2, k, i, j, n, id;
  int cn, cns;
  int nz = window->nz;

  /*-----------------------------------------------------------------*/
  /* Allocate and initialize                                         */
  for(cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Set the flag in the Delaunay structure                        */
    d[k]->ptf = 1;

    npt[k]->npt[id] = window->npe[c2] + 1;
    npt[k]->pt[id] = malloc(npt[k]->npt[id] * sizeof(int));
    npt[k]->pt[id][0] = -1;
  }
  /* Open boundary ghosts                                            */
  for (n = 0; n < window->nobc; n++) {
    int ee, e, ci;
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      ci = open->obc_e2[ee];
      c2 = window->m2d[c];
      k = window->s2k[ci] + 1;
      id = window->c2p[k][c2];
      npt[k]->npt[id] = window->npe[c2] + 1;
      npt[k]->pt[id] = malloc(npt[k]->npt[id] * sizeof(int));
      npt[k]->pt[id][0] = -1;
    }
  }
  wincon->nlmap = i_alloc_1d(window->szc);
  wincon->lmap = i_alloc_2d(window->npem+1, window->szc);

  /*-----------------------------------------------------------------*/
  /* Set the new points                                              */
  for(cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    k = window->s2k[c] + 1;
    id = window->c2p[k][c2];

    /* Fill the locations of Delaunay indices. The indices are       */
    /* consecutive around the centre.                                */
    npt[k]->pt[id][0] = window->c2p[k][c2];
    wincon->nlmap[c] = window->npe[c2] + 1;
    wincon->lmap[c][0] = c;
    for (j = 1; j <= window->npe[c2]; j++) {
      cn = window->c2c[j][c];
      cns = window->m2d[cn];
      npt[k]->pt[id][j] = window->c2p[k][cns];
      wincon->lmap[c][j] = cn;
      if(window->c2p[k][cns] == -1) 
	hd_quit("Invalid wet c2p map at c=%d: j=%d cn=%d %d %d\n",c, j, cn, window->s2i[c], window->s2j[c]);
    }
  }

  /* Set the sediments                                               */
  d[0]->ptf = 1;
  npt[0]->npt = npt[nz]->npt;
  npt[0]->pt = npt[nz]->pt;
}

/* END set_interp_ring()                                             */
/*-------------------------------------------------------------------*/


void print_tri_k(delaunay *d, int k)
{
  int n;
  FILE *fp;
  char buf[MAXSTRLEN];

  sprintf(buf, "tri%d.txt", k);
  fp = fopen(buf, "w");

  for (n = 0; n < d->ntriangles; n++) {
    triangle* t = &d->triangles[n];
    fprintf(fp, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
    fprintf(fp, "%f %f\n", d->points[t->vids[1]].x, d->points[t->vids[1]].y);
    fprintf(fp, "%f %f\n", d->points[t->vids[2]].x, d->points[t->vids[2]].y);
    fprintf(fp, "%f %f\n", d->points[t->vids[0]].x, d->points[t->vids[0]].y);
    fprintf(fp, "NaN NaN\n");
  }
  fclose(fp);  
}

void print_source_k(geometry_t *window,
		    int nvec,
		    int *vec,
		    double *cx,
		    double *cy,
		    int k,
		    int *s2k
		    )
{
  
  int cc, c, cs;
  FILE *fp;
  char buf[MAXSTRLEN];

  sprintf(buf, "tri%d.txt", k);
  fp = fopen(buf, "w");

  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    if (s2k[c] == k)
      fprintf(fp, "%f %f\n", cx[c], cy[c]);
  }
  fclose(fp);
}


/*-------------------------------------------------------------------*/
/* Routine to centre a vertical velocity.                            */
/*-------------------------------------------------------------------*/
void vel_center_w(geometry_t *window, /* Window geometry             */
		  window_t *windat,   /* Window data                 */
		  win_priv_t *wincon, /* Window constants            */
		  double *nw          /* Centered z velocity         */
		  )
{
  int cc, c;
  int p1;
  double wtop;
  double *w = (wincon->means & TRANSPORT) ? windat->wm : windat->w;

  /* Calculate the w values at the cell centre                       */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    p1 = window->zp1[c];
    wtop =
      (cc <= wincon->vcs) ? windat->wtop[window->m2d[c]] : w[p1];
    nw[c] = 0.5 * (w[c] + wtop);
    if (fabs(nw[c]) < SMALL)
      nw[c] = SMALL;
  }
}

/* END vel_center_w()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the vertical cell spacing centred on the grid layers          */
/*-------------------------------------------------------------------*/
void set_dzz(geometry_t *window, double *dzz)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c, cc, c2, cb, zm1, nn;
  double *cellz = wincon->cellz;
  double *eta = wincon->oldeta;

  /* Initialize in case a streamline lands in a dry cell, a valid    */
  /* (i.e. non inf) increment is computed.                           */
  memcpy(dzz, wincon->dz, window->szc * sizeof(double));
  memcpy(cellz, window->cellz, window->szc * sizeof(double));
  /* Overwrite with valid wet cell thicknesses                       */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    zm1 = window->zm1[c];
    dzz[zm1] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
    cellz[zm1] = 0.5 * (window->gridz[c] + max(window->botz[c2], window->gridz[zm1]));
  }
  /* Set the surface layer                                           */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = c2 = wincon->s1[cc];
    dzz[c] = 0.5 * wincon->dz[c];
    cellz[c] = 0.5 * (eta[c2] + window->gridz[c]);
    /* Set the surface and bottom boundary conditions                */
    while (c2 != window->zp1[c2]) {
      c2 = window->zp1[c2];
      dzz[c2] = dzz[c];
      cellz[c2] = cellz[c];
    }
    cb = wincon->i1[cc];
    dzz[window->zm1[cb]] = 0.5 * wincon->dz[cb];
    dzz[window->zm1[cb]] = wincon->dz[cb];
    cellz[window->zm1[cb]] = window->botz[c2];
  }
  /* Surface and bottom conditions over dry cells */
  for (cc = 1; cc <= wincon->ncdry_e1; cc++) {
    c = c2 = wincon->cdry_e1[cc];
    dzz[window->zm1[c]] = dzz[c];
    while (c2 != window->zp1[c2]) {
      c2 = window->zp1[c2];
      dzz[c2] = dzz[c];
    }
  }
  /* Set ghost cells                                                 */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    c2 = window->bin[cc];
    dzz[c] = dzz[c2];
  }
  /* OBC ghosts                                                      */

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    OBC_bgz_nogradb(window, open, dzz);
  }
}

void set_dzz_o(geometry_t *window, double *dzz)
{
  win_priv_t *wincon = window->wincon;
  int c, cc, c2, zm1;
  double top, bot;

  /* Initialize in case a streamline lands in a dry cell, a valid    */
  /* (i.e. non inf) increment is computed.                           */
  memcpy(dzz, wincon->dz, window->sgsiz * sizeof(double));
  /* Overwrite with valid wet cell thicknesses                       */
  /* Since a no gradient condition for tracer is created at the      */
  /* surface and bottom (i.e. no vertical interpolation is performed */
  /* at surface and bottom), then it doesn't matter whether the      */
  /* streamline increments are relative to (eta-cellz) & (cellz -    */
  /* botz) or (cellz[zp1] - cellz) at the surface and bottom. The    */
  /* latter is more convenient, since it avoids having different dzz */
  /* values in each quadrant of the transport grid.                  */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    zm1 = window->zm1[c];
    top = wincon->etamax;
    while(c != zm1) {
      bot = window->cellz[c];
      dzz[c] = top - bot;
      top = bot;
      c = zm1;
      zm1 =  window->zm1[c];
    }
    dzz[c] = dzz[window->zp1[c]];
  }
}

/* END set_dzz()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness and cells to process to for the */
/* layer conifuration corresponding to the updated elevation         */
/* (windat->eta and window->nsur_t). These dz and cells to process   */
/* vectors should be used in any routines called after tracers are   */
/* updated, e.g. diagnostics and interface routines (ecology,        */
/* sediments).                                                       */
/*-------------------------------------------------------------------*/
void reset_dzt(geometry_t *window,  /* Processing window */
	       window_t *windat,  /* Window data structure */
	       win_priv_t *wincon /* Window geometry / constants */
	       )
{
  int cc;                       /* Counter */
  int c, cs, cn;                /* Sparse coordinate */
  int zp1, zm1;                 /* Sparse cell above / below cell c */
  double top;                   /* Top face of a cell */
  double bot;                   /* Bottom face of a cell */
  int vc;                       /* Tracer cells to provess counter */
  int vcs;                      /* Surface tracer cells to process counter */

  /*-----------------------------------------------------------------*/
  /* Set dz at the surface to account for the updated elevation. */
  for(cc = 1; cc <= window->b2_t; cc++) {
    c = max(window->sur_t[cc], window->nsur_t[cc]);
    cs = window->m2d[c];
    cn = window->nsur_t[cc];
    zp1 = window->zp1[c];
    bot = max(window->gridz[c], window->botz[cs]);
    while(c > cn) {
      top = window->gridz[zp1];
      wincon->dz[c] = top - bot;
      bot = top;
      c = zp1;
      zp1 =  window->zp1[c];
    }
    wincon->dz[c] = windat->eta[cs] - bot;

    /* Set to the layer below the surface if eta coincides with a    */
    /* layer surface.                                                */
    if (!wincon->dz[c]) {
      cn = (window->zm1[cn] > window->bot_t[cc]) ? cn : window->zm1[cn];
      window->nsur_t[cc] = cn;
    }
    /* Set all cell thickness above the old surface equal to the     */
    /* surface thickness.                                            */
    c = cn;
    while (c > cs) {
      c = window->zp1[c];
      wincon->dz[c] = wincon->dz[cn];
    }
  }

  wincon->ncdry_e1 = 0;
  memset(wincon->cdry_e1, 0, window->sgsizS * sizeof(int));

  if (!wincon->sigma) {
    /*---------------------------------------------------------------*/
    /* Set the cells to process vector for the updated elevation. */
    /* Cells are arranged thus : */
    /* 1 to wincon->vca2 = surface layer wet cells */
    /* wincon->vca2+1 to wincon->vcs2 = surface layer OBC cells */
    /* wincon->vcs2+1 to wincon->vci2 = sub-surface wet cells */
    /* wincon->vci2+1 to wincon->vc2 = sub-surface OBC cells */
    vc = 1;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->nsur_t[cc];
      cs = window->m2d[c];
      top = windat->eta[cs];
      bot = DRY_FRAC * wincon->hmin + window->botz[cs];
      if (top > bot && c != window->zm1[c]) {
        wincon->s2[vc] = c;
	wincon->m2d[c] = c;
	if (wincon->trasc & (LAGRANGE|FFSL))
	  wincon->i1[vc] = window->bot_t[cc];
        wincon->i2[vc] = window->nsur_t[cc];
	wincon->s4[c] = vc;
        if (cc == window->v2_t)
          wincon->vca2 = vc;
        vc++;
      } else {
	wincon->ncdry_e1++;
	wincon->cdry_e1[wincon->ncdry_e1] = c;
	wincon->dz[c] = DRY_FRAC * wincon->hmin;
      }
    }
    vcs = vc - 1;
    wincon->vcs2 = vcs;

    /* Loop from the layer below the surface to the bottom and get */
    /* the cells to process.  */
    for (cc = 1; cc <= vcs; cc++) {
      cs = wincon->s2[cc];
      c = window->zm1[cs];
      zm1 = window->zm1[c];
      while (c != zm1) {
        wincon->s2[vc] = c;
	wincon->m2d[c] = cs;
        vc++;
        c = zm1;
        zm1 = window->zm1[c];
      }
      if (cc == wincon->vca2)
        wincon->vci2 = vc - 1;
    }
    vc--;
    wincon->vc2 = vc;
  }

  if (wincon->osl == 2 || wincon->osl == 4) {
    /* Set the surface layer                                         */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = cn = wincon->s1[cc];
      /* Set the surface and bottom boundary conditions              */
      while (cn != window->zp1[cn]) {
	cn = window->zp1[cn];
	wincon->dz[cn] = wincon->dz[c];
      }
      cn = wincon->i1[cc];
      wincon->dz[window->zm1[cn]] = wincon->dz[cn];
    }
    /* Set ghost cells                                               */
    for (cc = 1; cc <= window->nbpt; cc++) {
      c = window->bpt[cc];
      cn = window->bin[cc];
      wincon->dz[c] = wincon->dz[cn];
    }
    for (cc = 1; cc <= window->ngsed; cc++) {
      c = window->gsed_t[cc];
      cn = window->zp1[c];
      wincon->dz[c] = wincon->dz[cn];
    }
  }

  /* Set the wet cell mask                                           */
  memset(wincon->c1, 0, window->sgsiz * sizeof(char));
  for (cc = 1; cc <= vc; cc++) {
    c = wincon->s2[cc];
    wincon->c1[c] = 1;
  }
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    cs = window->bin[cc];
    wincon->c1[c] = wincon->c1[cs];
  }
}

/* END reset_dzt()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* If newly wetted cells were merged, then reset to the original     */
/* configuration. This is required for so that the newly wetted cell */
/* is included in vertical diffusion.                                */
/*-------------------------------------------------------------------*/
void reset_merged_dz(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
	       )
{
  int c, cc, cs, zp1;
  double top, bot;

  for(cc = 1; cc <= window->b2_t; cc++) {
    c = max(window->sur_t[cc], window->nsur_t[cc]);
    cs = window->m2d[c];
    zp1 = window->zp1[c];
    bot = max(window->gridz[c], window->botz[cs]);
    if (wincon->c1[c]) {
      while(c > window->nsur_t[cc]) {
	top = window->gridz[zp1];
	wincon->dz[c] = top - bot;
	bot = top;
	c = zp1;
	zp1 =  window->zp1[c];
      }
      wincon->dz[c] = windat->eta[cs] - bot;
    }
  }
}

/* END reset_merged_dz()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the total mass (water column and sediment) at   */
/* a cell. Uses a global multiplicitive fill method (e.g. Rood       */
/* (1987) p83. Mass added through OBC's is accounted for via         */
/* Marchesiello et al (2001), eqn. 16.                               */
/* Global filling is required in the transport mode since:           */
/* (a) The semi-Lagrangian scheme is not conservative,               */
/* (b) Using a mean or snapshot does not guarentee continuity. e.g.  */
/*     consider the case of a symmetrical tide at maximum flood      */
/*     where velocities are zero at the current time-step, but       */
/*     deta/dt should be non-zero. Using mean velocities and snap-   */
/*     shots, consider a progressive wave (no phase difference       */
/*     between eta and v) where mean v is zero and change in eta is  */
/*     non-zero.                                                     */
/* i.e. integral(deta/dt).dt != -d/dx(integral(DU).dt) since the     */
/* integral() on the rhs must be solved by parts because D=H+eta.    */
/*-------------------------------------------------------------------*/
void global_fill(geometry_t *window,  /* Window geometry             */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon,  /* Window constants            */
		 double *dtracer,     /* Source / sinks              */
		 int n,               /* Tracer number               */
		 int mode             /* Fill mode                   */
		 )
{
  int itmax = 10;
  int c, cc, cs;
  int nn, tn, m;
  double tmass, tm, tvol, dtmass, dtvol, dobc;
  double msf = 0.0;
  double eps = 1e-5;
  double msfmax = 2.0;
  double vol, v1, v2, v3;
  double *mxc = wincon->w4;
  double *mnc = wincon->w5;
  double *mask = wincon->w10;
  double *frac = wincon->w8;
  double mintr = (n > -1) ? wincon->mintr[n] : 0.0;
  double maxtr = (n > -1) ? wincon->maxtr[n] : 0.0;
  int *vec = wincon->s3;
  int nvec, nvecS;
  int *c2cc = wincon->s4;
  int first = 0;
  int nv = windat->ntr + 1;  /* Total vol at start of the time-step  */
  int nb = windat->ntr + 2;  /* Total bdry vol fluxes for time-step  */
  int errf = 0;              /* Produce volume errors                */
  int ord = 2;               /* Order of OBC interpolations          */
  int fbf = 1;               /* fbf = 0 : Forward differentiation    */
                             /* fbf = 1 : Backward differentiation   */
  int vecf = 1;              /* vecf = 0 : fill all cells            */
                             /* vecf = 1 : fill non-OBC cells only   */

  /*-----------------------------------------------------------------*/
  /* Get the cells to process                                        */
  if (vecf) {
    nvec = 1;
    for(cc = 1; cc <= wincon->vcs1; cc++) {
      vec[nvec] = wincon->s1[cc];
      nvec++;
    }
    nvecS = nvec - 1;
    for(cc = wincon->vcs + 1; cc <= wincon->vc1; cc++) {
      vec[nvec] = wincon->s1[cc];
      nvec++;
    }
    nvec--;
  } else {
    nvec = wincon->vc;
    nvecS = wincon->vcs;
    vec = wincon->s1;
  }
  /* Set the flag to compute the conservation error                  */
  if (wincon->fillf & (WEIGHTED|MONOTONIC)) errf = 1;

  if (mode == 0) {

    /*---------------------------------------------------------------*/
    /* Save the total mass corresponding to the start of the time-   */
    /* step.                                                         */
    /* Initialise.                                                   */
    memset(wincon->tmass, 0, windat->ntr * sizeof(double));
    /* Get total mass for each tracer and total volume.              */
    for(cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      cs = window->m2d[c];
      vol = window->cellarea[cs] * wincon->dz[c];
      for (nn = 0; nn < wincon->ntbdy; nn++) {
	tn = wincon->tbdy[nn];
	wincon->tmass[tn] += (windat->tr_wc[tn][c] * vol);
      }
    }
    /* Add the boundary fluxes                                       */
    if (!fbf) {
      for (nn = 0; nn < wincon->ntbdy; nn++) {
	tn = wincon->tbdy[nn];
        wincon->tmass[tn] += obc_tr_mass(window, windat, wincon, tn, 
					 ord, vecf);
      }
    }
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      tn = wincon->tbdy[nn];
      if (isnan(wincon->tmass[tn])) {
	hd_warn("global_fill: Tracer %s total mass = NaN at %5.2f\n",
		wincon->trname[tn], windat->t/86400);
	wincon->tmass[tn] = 0.0;
      }
    }
    /* Get the total volume                                          */
    if (wincon->fillf & (OBC_ADJUST|DIAGNOSE|DIAGNOSE_BGC)) {
      wincon->tmass[nv] = 0.0;
      for(cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	cs = window->m2d[c];
	wincon->tmass[nv] += window->cellarea[cs] * wincon->dz[c];
      }
    } else
      wincon->tmass[nv] = 1.0;
    /* Get the local volume budget                                   */
    if (errf) {
      memset(windat->vol_cons, 0, window->sgsizS * sizeof(double));
      windat->vol_cons[0] = 1;
      for(cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	cs = window->m2d[c];
	windat->vol_cons[cs] += window->cellarea[cs] * wincon->dz[c];
      }
    }
    return;

  } else if (mode == 1) {
    /*---------------------------------------------------------------*/
    /* Add the contributions from sources and sinks to the total     */
    /* mass.                                                         */
    wincon->b1 = wincon->tmass[n];
    for(cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      cs = window->m2d[c];
      wincon->tmass[n] -= dtracer[c];
    }
    return;

  } else {
    /*---------------------------------------------------------------*/
    /* Save the total mass corresponding to the end of the time-step */
    /* (i.e. when elevation is read in).                             */

    /*---------------------------------------------------------------*/
    /* Calculate the total mass.                                     */
    tmass = 0.0;
    for(cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      cs = window->m2d[c];
      vol = window->cellarea[cs] * wincon->dz[c];
      tmass += (windat->tr_wc[n][c] * vol);

      /* Maximum concentration due to advection.                     */
      /* Note: tracer value at time t was stored in wincon->w10 in   */
      /* advect_diffuse_lag().                                       */
      v1 = mxc[c] = max_lag(window, wincon->w10, wincon->s2[c2cc[c]]);
      /* Maximum concentration due to source/sinks                   */
      mxc[c] = max(wincon->w10[c], v1);
      /* Minimum concentration due to advection                      */
      v1 = mnc[c] = min_lag(window, wincon->w10, wincon->s2[c2cc[c]]);
      /* Minimum concentration due to source/sinks                   */
      mnc[c] = min(wincon->w10[c], v1);
      v1 = fabs(mxc[c] - mnc[c]);
      if(v1 && v1 < eps) eps = v1;
    }
    eps *= 10.0;

    /*---------------------------------------------------------------*/
    /* Get the total volume. This is only done once, i.e. not for    */
    /* every tracer.                                                 */
    if (mode & 2) {
      first = 1;
      /* Get the OBC volume flux. For fbf = 0 this is used at the    */
      /* next time-step to get mass fluxes.                          */
      vol = wincon->tmass[nb];
      wincon->tmass[nb] = obc_flux(window, windat, wincon, vecf);
      vol = (fbf) ? wincon->tmass[nb] : vol;
      /* Get the OBC volume adjustment factor if required            */
      if (wincon->fillf & (OBC_ADJUST|WEIGHTED|DIAGNOSE|DIAGNOSE_BGC)) {
	tvol = obc_msf(window, windat, wincon, vol, vec, nvec, vecf, &msf);
	dtvol = tvol - wincon->tmass[nv] - msf * vol;
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the mass flux through open boundaries                     */
    if (fbf) dobc = obc_tr_mass(window, windat, wincon, n, ord, vecf);

    /*---------------------------------------------------------------*/
    /* Get the conservation error for each water column as a         */
    /* fraction of the total volume at time t.                       */
    /* Note : V(t) - V(t-1) + div(U) = 0                             */
    if (errf) {
      memcpy(wincon->d1, windat->vol_cons, window->sgsizS * sizeof(double));
      if (windat->vol_cons[0]) {

	double *u1 = windat->u1;
	double *dzu1 = windat->dzu1;
	/* Get the local volume budget                               */
	for(cc = 1; cc <= nvec; cc++) {
	  int j, e, es;
	  double div;
	  c = vec[cc];
	  cs = window->m2d[c];
	  div = 0.0;
	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][c];
	    es = window->m2de[e];
	    div += u1[e] * dzu1[e] * window->h1au1[es];
	  }
	  windat->vol_cons[cs] -= (window->cellarea[cs] * wincon->dz[c] +
				   div * windat->dttr);
	}
	/* Convert to % of water column volume                       */
	if (wincon->fillf & MONOTONIC) {
	  for(cc = 1; cc <= nvecS; cc++) {
	    c = vec[cc];
	    cs = window->m2d[c];
	    /* Remove additional mass from the water column          */
	    v1 = window->cellarea[cs] * wincon->dz[c];
	  }
	  windat->vol_cons[0] = 0;
	}
	/* Convert to % of total volume change                       */
	if (wincon->fillf & WEIGHTED) {
	  /* Sum over the whole non-boundary domain                  */
	  windat->vol_cons[0] = 0.0;
	  for(cc = 1; cc <= nvecS; cc++) {
	    c = vec[cc];
	    cs = window->m2d[c];
	    windat->vol_cons[0] += windat->vol_cons[cs];
	  }
	  v1 = 0.0; v3 = 1e10;
	  for(cc = 1; cc <= nvecS; cc++) {
	    c = vec[cc];
	    cs = window->m2d[c];
	    windat->vol_cons[cs] *= (100.0/windat->vol_cons[0]);
	    v3 = min(v1, 0.01 * windat->vol_cons[cs]);
	    v1 += windat->vol_cons[cs];
	  }
	  windat->vol_cons[0] = 0;
	  if (fabs(v1 - 100.0) > SMALL)
	    hd_warn("global_fill: volume fractions do not sum to 100%% : %f\n", v1);

	  /* Get the weights                                         */
	  v2 = 0.0;	
	  for(cc = 1; cc <= nvec; cc++) {
	    c = vec[cc];
	    cs = window->m2d[c];
	    vol = window->cellarea[cs] * wincon->dz[c];

	    /* Cell volume weighted                                  */
	    frac[c] = vol / tvol;

	    /* Cell mass weighted                                    */
	    /*frac[c] = windat->tr_wc[n][c] * vol / tmass;*/

	    /* Positive volume error weighted cell thickness         */
	    /* Normalize to positive fractions                       */
	    v1 = (0.01 * windat->vol_cons[cs] - v3) / (1.0 - (double)nvecS * v3);
	    /* Scale to water depth                                  */
	    /*
	    frac[c] = v1 * wincon->dz[c] / 
	      max(windat->eta[cs] - window->botz[cs], wincon->hmin);
	    */

	    /* Volume error weighted to cell thickness               */
	    /*
	    frac[c] = 0.01 * windat->vol_cons[cs] * wincon->dz[c] /
	      max(windat->eta[cs] - window->botz[cs], wincon->hmin);
	    */
	    v2 += frac[c];
	  }
	  if (fabs(v2 - 1.0) > SMALL)
	    printf("global_fill: fill weights do not sum to 1 : %f\n", v2);
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Calculate the change in mass over the time-step               */
    /* Note: M(t)=msf.M(t+1)-OBC                                     */
    /*    => msf = (M(t)+OBC)/M(t+1)                                 */
    /*wincon->tmass[n] += dobc;*/
    if (fbf) tmass -= dobc;
    if (tmass)
      msf = wincon->tmass[n] / tmass;

    /* Total mass to add/subtract - diagnostic only                  */
    dtmass = tmass - wincon->tmass[n];

    /*---------------------------------------------------------------*/
    /* Print diagnostics to file if required                         */
    if (wincon->fillf & (DIAGNOSE|DIAGNOSE_BGC)) {
      if (master->t >= tstrans.tsout - DT_EPS) {
	if (first)
	  fprintf(tstrans.fp, "%f %f %f ", windat->days, 100.0*dtvol/tvol, 
		  wincon->tmass[nv]);
	if (wincon->fillf & DIAGNOSE)
	  fprintf(tstrans.fp, "%f %f %f %f ", 100.0*dobc/tmass, 
		  100.0*dtmass/tmass, msf, wincon->tmass[n] - wincon->b1);
	wincon->b2 += (wincon->tmass[n] - wincon->b1);
	if (wincon->fillf & DIAGNOSE_BGC)
	  fprintf(tstrans.fp, "%f %f %f %f ", dobc, tmass, msf, 
		  wincon->b2);
      }
      if (mode & 4) {
	fprintf(tstrans.fp, "\n");
	fflush(tstrans.fp);
	tstrans.tsout += tstrans.tsdt;
      }
    }

    /*---------------------------------------------------------------*/
    /* Redistribute the mass, keeping concentrations monotonic. This */
    /* is done iteratively, where if the adjusted mass is greater    */
    /* (or less) than the maximum (minimum) local concentrations,    */
    /* concentration is capped to the maximum (minimum) and residual */
    /* mass is spread over all cells not subject to capping at the   */
    /* next iteration.                                               */
    if (wincon->fillf & MONOTONIC) {
      /*-------------------------------------------------------------*/
      /* The MONOTONIC option is a multiplicative fill.              */
      double vol;
      double nfilled = 0;
      m = 0;
      v1 = 1.0;
      memset(mask, 0, window->sgsiz * sizeof(double));
      while (msf != 1.0 && v1 != 0.0 && m < itmax) {
	/* Sanity checks                                             */
	if (isnan(msf) || msf <= 0.0) {
	  hd_warn("global_fill: Tracer %s couldn't compute global scaling at %5.2f days: no filling\n",
		  wincon->trname[n], windat->t/86400);
	  return;
	}
	if (msf > msfmax) {
	  hd_warn("global_fill: Tracer %s global scaling %f considered too large at %5.2f days: no filling\n",
		  wincon->trname[n], msf, windat->t/86400);
	  return;
	}
	/* Re-distribute the change in mass                          */
	v1 = tm = tmass = 0.0;
	for(cc = 1; cc <= nvec; cc++) {
	  c = vec[cc];
	  cs = window->m2d[c];
	  vol = window->cellarea[cs] * wincon->dz[c];
	  v2 = windat->tr_wc[n][c] * msf;
	  v3 = fabs(mxc[c] - mnc[c]);
	  if (v3 && v3 < eps && !mask[c]) {
	    v1 += (v2 - windat->tr_wc[n][c]) * vol;
	    mask[c] = 1.0;
	  } else if (v2 > mxc[c] && !mask[c]) {
	    /* Check for adjusted mass > local maximum               */
	    v1 += (v2 - mxc[c]) * vol;
	    windat->tr_wc[n][c] = mxc[c];
	    mask[c] = 1.0;
	  } else if (v2 < mnc[c] && !mask[c]) {
	    /* Check for adjusted mass < local minimum               */
	    v1 -= (mnc[c] - v2) * vol;
	    windat->tr_wc[n][c] = mnc[c];
	    mask[c] = 1.0;
	  } else if (!mask[c]) {
	    /* Adjust the mass                                       */
	    windat->tr_wc[n][c] = v2;
	    tm += (windat->tr_wc[n][c] * vol);
	  }
	  /* Get the new toal mass (diagnostic only)                 */
	  tmass += (windat->tr_wc[n][c] * vol);
	}
	/* Calculate the change in mass over the time-step.          */
	/* Note : wincon->tmass[n] = tmass + v1, and                 */
	/*        wincon->tmass[n] = tmass at end of the while loop. */ 
	/* Also : SUM(tr_wc[n][c]*(1.0-msf)*cellarea[cs]*dz[c]) =    */
	/*        tmass - wincon->tmass[n]                           */
	/*if(n==2)printf("%f %d %f %f %f %f %20.18lf\n",windat->days,m,wincon->tmass[n],tmass + v1,tmass,v1,msf);*/
	msf = 1.0;
	if (tmass && !isnan(tmass)) msf = (tm + v1) / tm;
	m++;
      }

      if (m >= itmax) 
	hd_warn("global_fill: Tracer %s did not converge at %5.2f, residual mass = %5.3f\n",
		wincon->trname[n], windat->t/86400, v1);

    } else if (wincon->fillf & WEIGHTED) {
      /*-------------------------------------------------------------*/
      /* The WEIGHTED option is a additive fill.                     */
      double dtr, vol;
      /* Re-distribute the change in mass                            */
      m = 0;
      v1 = 1.0;
      memset(mask, 0, window->sgsiz * sizeof(double));
      while (v1 != 0.0 && m < itmax) {
	v1 = v3 = tmass = 0.0;
	for(cc = 1; cc <= nvec; cc++) {
	  c = vec[cc];
	  cs = window->m2d[c];
	  vol = window->cellarea[cs] * wincon->dz[c];
	  mxc[c] = min(mxc[c], maxtr);
	  mnc[c] = max(mnc[c], mintr);
	  /* Get the concentration to remove from each cell          */
	  dtr = frac[c] * dtmass / vol;
	  /* Get the new concentration                               */
	  v2 = windat->tr_wc[n][c] - dtr;
	  v3 += dtr * vol;
	  if (fabs(mxc[c] - mnc[c]) && 
	      fabs(mxc[c] - mnc[c]) < eps && !mask[c]) {
	    v1 += (v2 - windat->tr_wc[n][c]) * vol;
	    mask[c] = 1.0;
	  } else if (v2 > mxc[c] && !mask[c]) {
	    /* Check for adjusted mass > local maximum               */
	    v1 += (v2 - mxc[c]) * vol;
	    windat->tr_wc[n][c] = mxc[c];
	    mask[c] = 1.0;
	  } else if (v2 < mnc[c] && !mask[c]) {
	    /* Check for adjusted mass < local minimum               */
	    v1 -= (mnc[c] - v2) * vol;
	    windat->tr_wc[n][c] = mnc[c];
	    mask[c] = 1.0;
	  } else if (!mask[c]) {
	    /* Adjust the mass                                       */
	    windat->tr_wc[n][c] = v2;
	  }
	  /* Get the new toal mass                                   */
	  tmass += (windat->tr_wc[n][c] * vol);
	}
	/* Note : SUM(dtr*cellarea[cs]*dz[c])=dtmass                 */
	/*if(n==4)printf("%f %d %f %f %f : %f %f\n",windat->days,m,dtmass,v1,v3,tmass,dobc);*/
	m++;
	/* Get the new Total mass to add/subtract - diagnostic only. */
	dtmass = tmass - wincon->tmass[n];
      }
      if (m >= itmax) hd_warn("global_fill: Tracer %s did not converge at %5.2f, residual mass = %5.3f\n",
			    wincon->trname[n], windat->t/86400, v1);

    } else if (wincon->fillf & GLOBAL) {
      /* Re-distribute the change in mass                            */
      for(cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	windat->tr_wc[n][c] *= msf;
      }
    }
    for (cc = 1; cc < window->b3_t; cc++) {
      c = window->w3_t[cc];
      if (isnan(windat->tr_wc[n][c]))
	hd_quit("Global fill: Tracer %d NaN at %d, %f days\n",n,c,windat->days);
    }
  }
}


/* END global_fill()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the tracer concentration monotonically given a  */
/* multiplicative factor msf. Note; msf = 1 implies no adjustment to */
/* the tracer.                                                       */
/*-------------------------------------------------------------------*/
void monotonic_fill(geometry_t *window,  /* Window geometry          */
		    window_t *windat,    /* Window data              */
		    win_priv_t *wincon,  /* Window constants         */
		    int n,               /* Tracer number            */
		    int mode
		    )
{
  int c, cs, c2, cc, m;
  double eps;
  double vol, v1, v2, v3, tme, tm, tmass;
  int *S = wincon->s2;
  double *mxc = wincon->w4;
  double *mnc = wincon->w5;
  double *done = wincon->w10;
  double *otr = wincon->w10;
  double *dtracer = wincon->w1;
  double *omass; 
  double *ovol = wincon->p1;
  int *mask = wincon->s3;
  int *c2cc = wincon->s4;
  double msf = 1.0;      /* Mass scaling factor                      */
  double merr = 0.0;     /* Computed mass error                      */
  double mer = 0.0;      /* Diagnosed mass error                     */
  double m_in = 0.0;     /* Total inflow mass flux                   */
  double m_out = 0.0;    /* Total outflow mass flux                  */
  double m0 = 0.0;       /* Total mass at time t-1                   */
  double m1 = 0.0;       /* Total mass at time t                     */
  double mpss = 0.0;     /* Pointsource mass                         */
  double msfmax = 0.5;
  int itmax = 10;
  int pssex = 0;

  /*-----------------------------------------------------------------*/
  /* Calculate the minimum and maximum values                        */
  omass = wincon->d7;
  memset(omass, 0, window->sgsizS * sizeof(double));
  for(cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    /* Maximum concentration due to advection.                       */
    /* Note: tracer value at time t was stored in wincon->w10 in     */
    /* advect_diffuse_lag().                                         */
    v1 = mxc[c] = max_lag(window, otr, cs);
    /* Maximum concentration due to source/sinks                     */
    /*mxc[c] = max(otr[c], v1);*/
    /* Minimum concentration due to advection                        */
    v1 = mnc[c] = min_lag(window, otr, c);
    /* Minimum concentration due to source/sinks                     */
    /*mnc[c] = min(otr[c], v1);*/
    v1 = fabs(mxc[c] - mnc[c]);
    if(v1 && v1 < eps) eps = v1;
  }
  eps *= 10.0;
  eps = SMALL;

  /*-----------------------------------------------------------------*/
  /* Get the mass diagnostics for this tracer                        */
  /* The mass balance is:                                            */
  /*             m0 + m_in - m_out + mer = m1                        */
  /* The correct mass at t+1 is m1 - merr. We scale mass at t+1 so   */
  /* as to achieve the correct mass, i.e.                            */
  /*                   msf.m1 = m1 - merr                            */
  /* or:              msf = (m1 - merr) / m1                         */
  /* If pss is excluded from the scaling, then:                      */
  /*           mpss + (m1 - mpss)msf = m1 - merr     or              */
  /*         msf = (m1 - merr - mpss) / (m1 - mpss)                  */
  /*-----------------------------------------------------------------*/
  /* Surface interior cells (don't include boundary cells)           */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    vol = wincon->dz[c] * window->cellarea[c2];
    m_out -= otr[c] * windat->Vi[c];
    m0 += otr[c] * ovol[c];
    m1 += windat->tr_wc[n][c] * vol;
    if(dtracer[c]) mpss += windat->tr_wc[n][c] * vol;
    omass[c2] += windat->tr_wc[n][c] * vol;
  }
  /* Sub-surface interior cells                                      */
  for (cc = wincon->vcs + 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    vol = wincon->dz[c] * window->cellarea[c2];
    m_out -= otr[c] * windat->Vi[c];
    m0 += otr[c] * ovol[c];
    m1 += windat->tr_wc[n][c] * vol;
    if(dtracer[c]) mpss += windat->tr_wc[n][c] * vol;
    omass[c2] += windat->tr_wc[n][c] * vol;
  }
  /* Boundary cells                                                  */
  for (cc = wincon->vcs1+1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    m_in += otr[c] * (windat->Vi[c] + ovol[c]);
  }
  for (cc = wincon->vc1+1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    m_in += otr[c] * (windat->Vi[c] + ovol[c]);
  }
  /* Diagnosed mass error                                            */
  mer = m1 - m0 - m_in + m_out;
  /*merr = mer;*/

  /*-----------------------------------------------------------------*/
  /* Multiply the tracer concentrations by the volume error to get   */
  /* the exact mass conserving distribution, then adjust to maintain */
  /* monotonicity, where if the adjusted mass is greater (or less)   */
  /* than the maximum (minimum) local concentrations, concentration  */
  /* is capped to the maximum (minimum).                             */
  /* First set the mask of cells to scale                            */
  memset(done, 0, window->sgsiz * sizeof(double));
  if (pssex) {
    for(cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      if (dtracer[c]) done[c] = 1;
    }
  } else
    mpss = 0.0;

  /* Get the scaling factor to achieve perfect conservation          */
  msf = (m1 - mpss) ? (m1 - merr - mpss) / (m1 - mpss) : 1.0;
  if (fabs(msf - 1.0) > msfmax) {
    hd_warn("monotonic_fill: Tracer %s global scaling %f considered too large at %5.2f days: no global filling\n",
	    wincon->trname[n], msf, windat->t/86400);
    msf = 1.0;
  }

  /*-----------------------------------------------------------------*/
  /* Re-distribute the change in mass.                               */
  tm = tme = 0.0;
  memset(windat->vol_cons, 0, window->sgsizS * sizeof(double));
  for(cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];

    /* Get the cell volume                                           */
    vol = window->cellarea[cs] * wincon->dz[c];
    /* Get the mass conserving concentration                         */
    v2 = windat->tr_wc[n][c] * msf;
    tme += (v2 * vol);
    /* get the difference between minimum and maximum concentrations */
    v3 = fabs(mxc[c] - mnc[c]);
    if (v3 && v3 < eps && !done[c]) {
      done[c] = 1.0;
      windat->tr_wc[n][c] = 0.5 * (mxc[c] + mnc[c]);
      windat->vol_cons[cs] += (v2 - windat->tr_wc[n][c]) * vol;
    } else if (v2 > mxc[c] && !done[c]) {
      /* Check for adjusted mass > local maximum                     */
      windat->vol_cons[cs] += (v2 - mxc[c]) * vol;
      windat->tr_wc[n][c] = mxc[c];
      done[c] = 1.0;
    } else if (v2 < mnc[c] && !done[c]) {
      /* Check for adjusted mass < local minimum                     */
      windat->vol_cons[cs] -= (mnc[c] - v2) * vol;
      windat->tr_wc[n][c] = mnc[c];
      done[c] = 1.0;
    } else if (!done[c]) {
      /* Adjust the mass                                             */
      windat->tr_wc[n][c] = v2;
    }
    tm += (windat->tr_wc[n][c] * vol);
  }

  /* Get the volume error */
  for(cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    if (omass[cs])
      windat->vol_cons[cs] =  100.0 * windat->vol_cons[cs] / omass[cs];
    else
      windat->vol_cons[cs] =  0.0;
  }
}

/* END monotonic_fill()                                             */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Computes the volume at time t-1. This uses the dz and cells to   */
/* process vector set in reset_dz(); i.e. wincon->dz and            */
/* wincon->s2.                                                      */
/*------------------------------------------------------------------*/
void calc_volume(geometry_t *window, /* Window geometry             */
		 window_t *windat,   /* Window data                 */
		 win_priv_t *wincon  /* Window constants            */
		)
{
  int c, cc, cs;
  int vc = (wincon->vc2) ? wincon->vc2 : wincon->vc;
  int *cells = (wincon->vc2) ? wincon->s2 : wincon->s1;

  /* Calculate the volume                                            */
  memset(wincon->p1, 0, window->sgsiz * sizeof(double));
  for(cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    wincon->p1[c] = window->cellarea[cs] * wincon->dz[c];
  }
  /* Newly dried cells */
  for (cc = 1; cc <= wincon->ncdry_e1; cc++) {
    c = wincon->cdry_e1[cc];
    cs = window->m2d[c];
    wincon->p1[c] = window->cellarea[cs] * wincon->dz[c];
  }
}

/* END calc_volume()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find a maximum value in cells surrounding the origin   */
/* of the streamline.                                                */
/*-------------------------------------------------------------------*/
double max_lag(geometry_t *window, double *tr, int c)
{
  win_priv_t *wincon= window->wincon;
  int c1, j;
  double v1 = -1e10;

  for (j = 0; j <= wincon->nlmap[c]; j++) {
    c1 = wincon->lmap[c][j];
    v1 = max(v1, tr[c1]);
  }
  return(v1);
}

/* END max_lag()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to find a minimum value in cells surrounding the origin   */
/* of the streamline.                                                */
/*-------------------------------------------------------------------*/
double min_lag(geometry_t *window, double *tr, int c)
{
  win_priv_t *wincon= window->wincon;
  int c1, j;
  double v1 = 1e10;

  for (j = 0; j <= wincon->nlmap[c]; j++) {
    c1 = wincon->lmap[c][j];
    v1 = min(v1, tr[c1]);
  }
  return(v1);
}

/* END min_lag()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute totals mass flux through open boundaries                  */
/*-------------------------------------------------------------------*/
double obc_tr_mass(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon, /* Window constants           */
		   int n,              /* Tracer number              */
		   int ord,            /* Order of interpolation     */
		   int vecf            /* Cells to process flag      */
		   )
{
  int c, cc, nn;
  double t;
  int ci, cp, cm;
  int *ct, *cp1, *cm1, cv;
  double *tr = windat->tr_wc[n];
  double *vel;
  double dobc = 0.0;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    if(!wincon->trinfo_3d[n].advect) continue; 

    if (vecf) {
      ct = open->oi1_t;
      cp1 = open->oi2_t;
      cm1 = open->obc_t;
    } else {
      ct = open->obc_t;
      cp1 = open->oi1_t;
      cm1 = open->obc_t;
    }
    if (ord == 1) {
      vel = windat->u1;
      cv = (vecf) ? open->ocodex : open->ocodey;
    }
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = ct[cc];
      cp = cp1[cc];
      cm = cm1[cc];
      if (ord == 1) {          /* 1st order                          */
	ci = window->c2e[cv][c];
	t = (window->eSc[cv][window->m2d[c]] * vel[ci] > 0.0) ? tr[cm] : tr[c];
      } else if (ord == 2) {   /* 2nd order                          */
	t = 0.5 *(tr[c] + tr[cm]);
      } else if (ord == 4) {   /* 4th order                          */
	t = (7.0 * (tr[cm] + tr[c]) - (tr[cp] + tr[cm])) / 12.0;
      } else if (ord == 0)
	t = 1.0;
      dobc += (open->dum[cc] * t);
    }
  }
  return(dobc);
}

/* END obc_tr_mass()                                                    */
/*-------------------------------------------------------------------*


/*-------------------------------------------------------------------*/
/* Compute totals mass flux through open boundaries                  */
/*-------------------------------------------------------------------*/
double obc_flux(geometry_t *window, /* Window geometry               */
		window_t *windat,   /* Window data                   */
		win_priv_t *wincon, /* Window constants              */
		int vecf            /* Cells to process flag         */
		)
{
  int c, cc, e, es, nn, j;
  double *vel;
  double *hat, *dz;
  int *cv;
  double dobc = 0.0;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    open->bflux_3d = 0.0;

    vel = windat->u1;
    hat = window->h1au1;
    dz = windat->dzu1;
    cv = (vecf) ? open->oi1_t : open->obc_t;
    j = (vecf) ? open->ocodex : open->ocodey;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = cv[cc];
      e = window->c2e[j][c];
      es = window->m2de[e];
      open->dum[cc] = (window->eSc[j][window->m2d[c]] * vel[e] * hat[es] * dz[e] * windat->dttr);
      open->bflux_3d += open->dum[cc];
    }
    dobc += open->bflux_3d;
  }
  return(dobc);
}

/* END obc_flux()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Scales the totals mass flux through open boundaries               */
/*-------------------------------------------------------------------*/
void obc_scale(geometry_t *window, /* Window geometry                */
	       window_t *windat,   /* Window data                    */
	       win_priv_t *wincon, /* Window constants               */
	       double msf,         /* Boundary scaling               */
	       int vecf            /* Cells to process flag          */
	       )
{
  int c, cc, nn;
  int *cv;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    cv = (vecf) ? open->oi1_t : open->obc_t;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = cv[cc];
      open->dum[cc] *= msf;
    }
  }
}

/* END obc_scale()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the scaling factor to be applied to OBC volume     */
/* fluxes so that domain volume is conserved.                        */
/* Ideally tvol - wincon->tmass[nv] = wincon->tmass[nb], but as a    */
/* result of point (b) above, error is introduced.                   */
/* Since volume fluxes are not associated with error from the        */
/* Lagrange scheme (point (a) above), the boundary volume error      */
/* (i.e. the domain volume divergence) can be computed and used      */
/* to adjust boundary tracer fluxes.                                 */
/*-------------------------------------------------------------------*/
double obc_msf(geometry_t *window, /* Window geometry                */
	       window_t *windat,   /* Window data                    */
	       win_priv_t *wincon, /* Window constants               */
	       double obcvol,      /* OBC volume flux                */
	       int *vec,           /* Cells to process flag          */
	       int nvec,           /* Number of cells to process     */
	       int vecf,           /* OBC location flag              */
	       double *msf         /* OBC volume scaling             */
	       )
{
  int c, cc, cs;
  int nv = windat->ntr + 1;  /* Total vol at start of the time-step  */
  int nb = windat->ntr + 2;  /* Total bdry vol fluxes for time-step  */
  double tvol = 0.0;

  *msf = 1.0;

  /*-----------------------------------------------------------------*/
  /* Get the total volume in the domain interior                     */
  for(cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cs = window->m2d[c];
    tvol += window->cellarea[cs] * wincon->dz[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get the adjustment factor for boundary fluxes.                  */
  /* If wincon->tmass[nb]>0 then there is a net input of volume.     */
  if (wincon->fillf & OBC_ADJUST) {
    *msf = (obcvol) ? (tvol - wincon->tmass[nv]) / obcvol : 1.0;
    if(*msf < 0)
      hd_warn("global_fill: boundary flux is opposite sign to volume change : %4.2f days\n", windat->days);
    obc_scale(window, windat, wincon, *msf, vecf);
  }
  return(tvol);
}

/* END obc_msf()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to force mass conservation over a timestep                */
/*-------------------------------------------------------------------*/
void global_tfill(geometry_t *window,  /* Window geometry            */
		  window_t *windat,    /* Window data                */
		  win_priv_t *wincon,  /* Window constants           */
		  double *dtracer,     /* Sourcesink & hdiff fluxes  */
		  int n                /* Tracer number              */
		  )
{
  int c, cc, e, ee, es;
  int nn, tn;
  double vol;
  double tmass[windat->ntr];
  double tsmass[windat->ntr];

  /* Add mass from sources/sinks to the total mass at the previous   */
  /* timestep.                                                       */
  if (dtracer && n >= 0) {
    for(cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      wincon->tmass[n] += dtracer[c];
    }
    return;
  }

  calc_tmass(window, windat, wincon, tmass, tsmass);

  /* Adjust the mass at the previous timestep for sediment fluxes    */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    tn = wincon->tbdy[nn];
    wincon->tmass[tn] -= (tsmass[tn] - wincon->tsmass[tn]);
  }
  /* Adjust the mass at the previous timestep for open boundary      */
  /* fluxes.                                                         */
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    double *vel, v;
    double *hat;
    int sc, *cv;

    vel = windat->u1;
    hat = window->h2au1;

    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      c = open->obc_e2[ee];
      v = vel[e];
      es = window->m2de[e];
      sc = -(double)open->dir[ee];
      vol = sc * v * hat[es] * wincon->dz[e] * windat->dttr;
      for (nn = 0; nn < wincon->ntbdy; nn++) {
	tn = wincon->tbdy[nn];
	wincon->tmass[tn] += vol * windat->tr_wc[tn][c];
      }
    }
  }

  /* Calculate the change in mass over the time-step                 */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    double msf = 1.0;
    tn = wincon->tbdy[nn];
    if (tmass[tn]) msf = wincon->tmass[tn] / tmass[tn];
    /* Re-distribute the change in mass                              */
    for(cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      windat->tr_wc[tn][c] *= msf;
    }
  }

  /* Calculate the new (adjusted) total mass                         */
  calc_tmass(window, windat, wincon, wincon->tmass, wincon->tsmass);
}

/* END global_tfill()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates total mass in the water column and sediments           */
/*-------------------------------------------------------------------*/
void calc_tmass(geometry_t *window, /* Window geometry               */
		window_t *windat,   /* Window data                   */
		win_priv_t *wincon, /* Window constants              */
		double *tmass,      /* Total WC mass                 */
		double *tsmass      /* Total sediment mass           */
		)
{
  int c, cc, cs, k;
  int nn, tn;
  double vol;

  memset(tmass, 0, windat->ntr * sizeof(double));
  memset(tsmass, 0, windat->ntr * sizeof(double));
  /* Save the total mass from this time-step                         */
  for(cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    vol = window->cellarea[cs] * wincon->dz[c];
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      tn = wincon->tbdy[nn];
      tmass[tn] += windat->tr_wc[tn][c] * vol;
    }
  }

  /* Save the total sediment mass                                    */
  for(cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    for (k = 0; k <= window->sednz-1; k++) {
      for (nn = 0; nn < wincon->ntbdy; nn++) {
	tn = wincon->tbdy[nn];
	if (wincon->trinfo_3d[tn].type & SEDIM) {
	  int tn_sed = tracer_find_index(wincon->trname[tn], wincon->nsed,
					 wincon->trinfo_sed);
	  tsmass[tn] += windat->tr_sed[tn_sed][k][cs] * window->cellarea[cs]*
	    (window->gridz_sed[k+1][cs] - window->gridz_sed[k][cs]);
	}
      }
    }
  }
}

/* END calc_tmass()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the elevation required for conservation on the tracer  */
/* timestep, given a 3D velocity field for the transport mode. Note  */
/* that elevation is still read from file every time-step; this eta  */
/* is assigned the elevation at the start of the time-step, and the  */
/* elevation at the end of the time-step is computed from the        */
/* divergence of the depth averaged 3D velocity fluxes. The vertical */
/* velocity may be similarly computed from the divergence of the 3D  */
/* velocity flux. This ensures dynamic consistency over the timestep */
/* and thus conservation of tracers. Without re-initializing eta to  */
/* the file value at the start of the timestep, the solution may     */
/* excessively drift over time. Note that conservation may not       */
/* be achieved in the long term in the elevation read from file is   */
/* not consistent with the velocity field, e.g. using means rather   */
/* than snapshots.                                                   */
/*-------------------------------------------------------------------*/
void calc_tracer_eta(geometry_t *window, /* Window geometry          */
		     window_t *windat,   /* Window data              */
		     win_priv_t *wincon  /* Window constants         */
		     )
{
  double *f1;                   /* Edge transports                   */
  double colflux;               /* Flux divergence                   */
  int c, cc, cs, ee, e, es;     /* Sparse coodinate / counter        */
  int bn;                       /* OBC counter                       */
  open_bdrys_t **open = window->open; /* Window OBC structure        */
  double *owaterss2d = wincon->d7;

  /*-----------------------------------------------------------------*/
  /* Set the elevation at the start of the time-step                 */
  memcpy(windat->etab, windat->eta, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Integrate the e1 3D fluxes through the water column.            */
  f1 = wincon->d1;
  memset(f1, 0, window->sgsizS * sizeof(double));
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    f1[es] += windat->u1flux3d[e] * windat->dttr;
  }

  /* Calculate the new elevation                                     */
  /*set_map_eta(window);*/
  for (cc = 1; cc <= window->b2_t; cc++) {
    int e, j;
    double eta;
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = 0.0;
    for (j = 1; j <= window->npe[c]; j++) {
      e = window->c2e[j][c];
      colflux += window->eSc[j][c] * f1[e];
    }

    /* Add source/sink flows to colflux. Note that colflux is the    */
    /* accumulated volume flowing out of this column, whereas        */
    /* waterss2d is the flow (m3 s-1) into the column. Hence the     */
    /* minus sign.  */
    /* colflux -= owaterss2d[c] * windat->dttr; */
    colflux -= windat->waterss2d[c] * windat->dttr;

    /* Calculate new etat value                                      */
    windat->eta[c] = max(wincon->oldeta[c] - colflux / window->cellarea[c],
			    window->botz[c]);
  }
  /* Use the pointsource volume at the end of the time-step          */
  memcpy(owaterss2d, windat->waterss2d, window->sgsizS * sizeof(double));
}

/* END calc_tracer_eta()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks that the provided sea level is equal to the vertical       */
/* integral of fluxes. If not, then the sea level is reset to that   */
/* consistent with the flux divergence.                              */
/*-------------------------------------------------------------------*/
void check_tracer_eta(geometry_t *window, /* Window geometry         */
		      window_t *windat,   /* Window data             */
		      win_priv_t *wincon  /* Window constants        */
		      )
{
  double *f1;                   /* Edge transports                   */
  double colflux;               /* Flux divergence                   */
  int c, cc, cs, ee, e, es;     /* Sparse coodinate / counter        */
  double eta;                   /* Computed sea level                */
  double eps = 1e-5;            /* Tolerance                         */

  /*-----------------------------------------------------------------*/
  /* Integrate the e1 3D fluxes through the water column.            */
  f1 = wincon->d1;
  memset(f1, 0, window->szeS * sizeof(double));
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    f1[es] += windat->u1flux3d[e] * windat->dttr;
  }

  /* Calculate the new elevation                                     */
  set_map_eta(window);
  for (cc = 1; cc <= window->b2_t; cc++) {
    int e, j;
    double eta;
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = 0.0;
    for (j = 1; j <= window->npe[c]; j++) {
      e = window->c2e[j][c];
      colflux += window->eSc[j][c] * f1[e];
    }

    /* Add source/sink flows to colflux. Note that colflux is the    */
    /* accumulated volume flowing out of this column, whereas        */
    /* waterss2d is the flow (m3 s-1) into the column. Hence the     */
    /* minus sign.  */
    colflux -= windat->waterss2d[c] * windat->dttr;

    /* Calculate new etat value                                      */
    eta = max(wincon->oldeta[c] - colflux / window->cellarea[c],
			    window->botz[c]);
    windat->eta[c] = eta;

    if (fabs(eta-windat->eta[c]) > eps) {
      emstag(LDEBUG,"check_tracer_eta","Sea level reset to flux divergence at %4.2f days: %f replaced with %f at (%d %d)\n", windat->days, windat->eta[c], eta, window->s2i[c], window->s2j[c]);
      windat->vol_cons[c] = (windat->vol_cons[c] == 1.0) ? 3.0 : 2.0;
      windat->eta[c] = eta;
    }
  }
}

/* END check_tracer_eta()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the totals output timeseries file.              */
/*-------------------------------------------------------------------*/
void init_trans(master_t *master,    /* Master data structure */
		geometry_t **window, /* Processing window */
		win_priv_t **wincon  /* Window geometry / constants */
  )
{
  int n, nn, m, trn;
  int ntr = wincon[1]->ntbdy;

  /*---------------------------------------------------------------*/
  /* Initialize the output file */
  if ((master->fillf & (DIAGNOSE|DIAGNOSE_BGC)) && tstrans.fp == NULL) {
    if (strlen(master->opath))
      sprintf(tstrans.pname, "%strans.ts", master->opath);
    else
      sprintf(tstrans.pname, "trans.ts");
    if ((tstrans.fp = fopen(tstrans.pname, "w")) == NULL)
      hd_quit("init_trans: Can't open file %s\n", tstrans.pname);
    tstrans.tsdt = 3600.0;
    tstrans.master = master;
    tstrans.i = 1;
    tstrans.j = 1;
    tstrans.k = 1;
    tstrans.tsout = master->t;
    fprintf(tstrans.fp, "## COLUMNS %d\n", 3 + 4 * ntr);
    fprintf(tstrans.fp, "##\n");
    fprintf(tstrans.fp, "## COLUMN1.name           Time\n");
    fprintf(tstrans.fp, "## COLUMN1.long_name      Time\n");
    fprintf(tstrans.fp, 
	    "## COLUMN1.units          days since 1990-01-01 00:00:00 +10\n");
    fprintf(tstrans.fp, "## COLUMN1.missing_value  -999\n");
    fprintf(tstrans.fp, "##\n");

    fprintf(tstrans.fp, "## COLUMN2.name           dvol\n");
    fprintf(tstrans.fp, "## COLUMN2.long_name      Volume error as %% of total volume\n");
    
    fprintf(tstrans.fp, "## COLUMN2.units          %%\n");
    fprintf(tstrans.fp, "## COLUMN2.missing_value  0.000000\n");
    fprintf(tstrans.fp, "##\n");
    m = 3;
    for (nn = 0; nn < ntr; nn++) {
      n = wincon[1]->tbdy[nn];
      if((trn = tracer_find_index(master->trname[n], master->ntr,
				  master->trinfo_3d)) >= 0) {
	if (master->fillf & DIAGNOSE) {
	  fprintf(tstrans.fp, "## COLUMN%d.name           %s_OBC%%\n", m, 
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.long_name      OBC mass flux for %s as %% of total mass\n", m,
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.units          %%\n", m);
	} else {
	  fprintf(tstrans.fp, "## COLUMN%d.name           %s_OBC\n", m, 
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.long_name      OBC mass flux for %s\n", m, master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.units          kg\n", m);
	}
	fprintf(tstrans.fp, "## COLUMN%d.missing_value  0.000000\n", m);
	fprintf(tstrans.fp, "##\n");
	m++;

	if (master->fillf & DIAGNOSE) {
	  fprintf(tstrans.fp, "## COLUMN%d.name           %s_error\n", m, 
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.long_name      Mass error for %s as %% of total mass\n", m,
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.units          %%\n", m);
	} else {
	  fprintf(tstrans.fp, "## COLUMN%d.name           %s_mass\n", m, 
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.long_name      Mass for %s \n", m,
		  master->trname[n]);
	  fprintf(tstrans.fp, "## COLUMN%d.units          kg\n", m);
	}
	/*fprintf(tstrans.fp, "## COLUMN%d.units          %s\n", m,
	  master->trinfo_3d[trn].units);*/
	fprintf(tstrans.fp, "## COLUMN%d.missing_value  0.000000\n", m);
	fprintf(tstrans.fp, "##\n");
	m++;
	  
	fprintf(tstrans.fp, "## COLUMN%d.name           %s_msf\n", m, 
		master->trname[n]);
	fprintf(tstrans.fp, "## COLUMN%d.long_name      Scaling for %s\n", m,
		master->trname[n]);
	fprintf(tstrans.fp, "## COLUMN%d.units            \n", m);
	fprintf(tstrans.fp, "## COLUMN%d.missing_value  0.000000\n", m);
	fprintf(tstrans.fp, "##\n");
	m++;
	
	fprintf(tstrans.fp, "## COLUMN%d.name           %s_input\n", m, 
		master->trname[n]);
	fprintf(tstrans.fp, "## COLUMN%d.long_name      Mass pss input for %s\n", m,
		master->trname[n]);
	fprintf(tstrans.fp, "## COLUMN%d.units          kg\n", m);
	    fprintf(tstrans.fp, "## COLUMN%d.missing_value  0.000000\n", m);
	    fprintf(tstrans.fp, "##\n");
	    m++;
      } else
	hd_quit("Can't recognize tracer %s for totals.\n", master->trname[n]);
    }
  }
}

/* END init_trans()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void check_transport(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int c, cc, cs, c1, j;
  int *S = wincon->s2;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    if (wincon->dz[c] <= 0)
      hd_quit("transport check: dz (%5.2f) <= 0 at %c(%d %d %d)\n",wincon->dz[c],
	      c, window->s2i[c], window->s2j[c], window->s2k[c]);
    if (cs <=0 || cs > window->enon)
      hd_quit("transport check: invalid source cell %d at %c(%d %d %d)\n", cs,
	      c, window->s2i[c], window->s2j[c], window->s2k[c]);
  }
}

/* END check_transport()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void trans_data_nan(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  geometry_t *tpg;
  window_t *tpd;
  win_priv_t *tpc;
  int cc, c, ee, e;

  tpg = window->trans->tpg;
  tpd = window->trans->tpd;
  tpc = window->trans->tpc;

  if (wincon->trasc == FFSL) {
    for (cc = 1; cc <= tpg->n2_t; cc++) {
      c = tpg->w2_t[cc];
      if (isnan(tpd->eta[c])) tpd->eta[c] = 0.0;
    }
    for (cc = 1; cc <= tpg->n3_t; cc++) {
      c = tpg->w3_t[cc];
      if (isnan(tpd->w[c])) tpd->w[c] = 0.0;
    }
    for (ee = 1; ee <= tpg->n3_e1; ee++) {
      e = tpg->w3_e1[ee];
      if (isnan(tpd->u1[e])) tpd->u1[e] = 0.0;
    }
  } else {
    for (cc = 1; cc <= tpg->n2_t; cc++) {
      c = tpg->w2_t[cc];
      if (isnan(tpd->eta[c]) || fabs(tpd->eta[c]) > tpc->etamax)
	tpd->eta[c] = 0.0;
    }
    for (cc = 1; cc <= tpg->n3_t; cc++) {
      c = tpg->w3_t[cc];
      if (isnan(tpd->w[c]) || fabs(tpd->w[c]) > tpc->velmax)
	tpd->w[c] = 0.0;
    }
    for (ee = 1; ee <= tpg->n3_e1; ee++) {
      c = tpg->w3_e1[ee];
      if (isnan(tpd->u1[e]) || fabs(tpd->u1[e]) > tpc->velmax)
	tpd->u1[e] = 0.0;
    }
  }
}


void trans_data_check(master_t *master, geometry_t **window,
		      window_t **windat, win_priv_t **wincon)
{
  geometry_t *geom = master->geom;
  int nwindows = master->nwindows;
  geometry_t *tpg;
  window_t *tpd;
  win_priv_t *tpc;
  int nn, n;
  int cc, c, cp, ee, e;

  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];
    tpg = window[n]->trans->tpg;
    tpd = window[n]->trans->tpd;
    tpc = window[n]->trans->tpc;

    for (cc = 1; cc <= tpg->n2_t; cc++) {
      c = tpg->w2_t[cc];
      if (isnan(tpd->eta[c]) || fabs(tpd->eta[c]) > tpc->etamax)
	check_warn(tpg, "sea level", tpd->eta[c], master->days, c, 2);
    }
    for (ee = 1; ee <= tpg->n3_e1; ee++) {
      c = tpg->w3_e1[ee];
      if (isnan(tpd->u1[e]) || fabs(tpd->u1[e]) > tpc->velmax)
	check_warn(tpg, "u1 velocity", tpd->u1[e], master->days, e, 3);
    }
    for (cc = 1; cc <= tpg->n3_t; cc++) {
      c = tpg->w3_t[cc];
      if (isnan(tpd->Kz[c]) || tpd->Kz[c] < 0.0)
	check_warn(tpg, "Kz", tpd->Kz[c], master->days, c, 3);
    }
  }
}

void check_warn(geometry_t *geom, 
		char *var, 
		double val, 
		double t, 
		int c, int dim)
{
  printf("Corrupt %s transport data (%f) supplied\n", var, val);
  hd_warn("Corrupt %s transport data (%f) supplied\n", var, val);
  if (dim == 2)
    hd_quit("Time = %f, location = %d(%d %d)\n", t, c, geom->s2i[c],
	    geom->s2j[c]);
  else
    hd_quit("Time = %f, location = %d(%d %d %d)\n", t, c, geom->s2i[c],
	    geom->s2j[c], geom->s2k[c]);

}



/*-------------------------------------------------------------------*/
/* Gets the weights for interpolations of order 1 & 3. Uses Eqn. 5   */ 
/* McDonald (1984) Mon. Wea. Rev. 112, 1267-1275. The distances are  */
/* re-mapped to local distances, where h[0] = 0. Note that (Eqn. 7)  */
/* hx[i-1] < p < h[i], hy[j-1] < q < hy[j] and hz[i] < r < hz[i+1].  */
/*-------------------------------------------------------------------*/
void weights_v(geometry_t *window, 
	       int co,             /* Source cell                    */
	       double rin,         /* z distance from dest origin    */
	       double *dzz,        /* Distance between layer centres */
	       double *wgt,
	       int osl
		)
{
  win_priv_t *wincon = window->wincon;
  int n = osl;
  int k, m;
  int ck;
  double r;
  double hz[5];

  /* Set up the normalized distances */
  ck = co;
  r = rin;
  if (n > 1) {
    ck = window->zm1[ck];
    r += dzz[ck];
  }
  hz[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hz[m] = hz[m-1] + dzz[ck];
    ck = window->zp1[ck];
  }

  /* Loop to get the weights */
  for (k = 0; k <= n; k++) {
    wgt[k] = 1.0;
    for (m = 0; m <= n; m++) {
      if (m != k)
	wgt[k] *= (r - hz[m]) / (hz[k] - hz[m]);
    }
  }
}

/* END weights_v()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set a no-gradient over open boundaries at RIGHT and    */
/* FRONT edges. Although the maps are made self mapping over these   */
/* boundaries, some situations may occur where a ghost cell adjacent */
/* to the OBC is the source cell for the tri-diagonal (whose maps    */
/* are not made self-mapping) and hence some OBC cells may be        */
/* included in the tri-diagonal. The altenative is to make these     */
/* ghost cells self-mapping.                                         */
/*-------------------------------------------------------------------*/
void set_lateral_OBC_tr(geometry_t *window, int ntr, double **tr)
{
  win_priv_t *wincon = window->wincon;
  int ee, c, c2, tn, nn, n;
 
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      c2 = open->obc_e2[ee];
      for (n = 0; n < wincon->ntbdy; n++) {
	tn = wincon->tbdy[n];
	if (!(open->bcond_tra[tn] & (TRCONC|TRCONF|TRFLUX)))
	  tr[tn][c] = tr[tn][c2];
      }
    }
  }
}

/* END set_lateral_OBC_tr()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the mean velocity and elevation from transport input     */
/*-------------------------------------------------------------------*/
void calc_means_t(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int c, cc, cs;
  double t = windat->dttr;

  if (wincon->means & VEL3D) {
    if (windat->u1m && windat->u2m) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	windat->u1m[c] = (windat->u1m[c] * windat->meanc[cs] + 
			  windat->u[c] * t) / (windat->meanc[cs] + t);
	windat->u2m[c] = (windat->u2m[c] * windat->meanc[cs] + 
			  windat->v[c] * t) / (windat->meanc[cs] + t);
      }
    }
    if (windat->wm) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	windat->wm[c] = (windat->wm[c] * windat->meanc[cs] + 
			 windat->w[c] * t)  / (windat->meanc[cs] + t);
      }
    }
  }
  if (wincon->means & ETA_M && windat->etam) {
    double to;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      to = fabs(windat->meanc[window->m2d[c]]);
      windat->etam[c] = (windat->etam[c] * windat->meanc[c] + 
			 windat->eta[c] * t) / (windat->meanc[c] + t);
    }
  }
  if (wincon->means & VEL2D) {
    if (windat->u1am) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	windat->u1am[c] = (windat->u1am[c] * windat->meanc[c] + 
			  windat->uav[c] * t) / (windat->meanc[c] + t);
	windat->u2am[c] = (windat->u2am[c] * windat->meanc[c] + 
			  windat->vav[c] * t) / (windat->meanc[c] + t);
      }
    }
  }
}

/* END calc_means_t()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-calculates fluxes and velocities for tratio < 1.               */
/*-------------------------------------------------------------------*/
void recalc_vel(geometry_t *window, /* Processing window      */
		window_t *windat,   /* Window data structure  */
		win_priv_t *wincon  /* Window constants       */
		)
{
  int zm1, ee, e, es, eb, e2;
  double sum, f, dz;

  /*-----------------------------------------------------------------*/
  /* Reset the e1 fluxes and velocity                                */
  memcpy(windat->u1flux3d, windat->u1vm, window->sgsiz * sizeof(double));
  /* Get the cell spacing and new surface coordinate for u1          */
  set_dz_at_u1(window, windat, wincon);
  /* Reset the flux at the surface cell, and zero above the surface  */
  for (ee = 1; ee < window->b2_e1; ee++) {
    e = es = window->sur_e1[ee];
    eb = window->bot_e1[ee];
    e2 = window->m2de[e];
    sum = 0.0;
    /* Reset the flux above the surface                              */
    while (e != window->zp1e[e]) {
      e = window->zp1e[e];
      sum += windat->u1flux3d[e];
      windat->u1flux3d[e] = 0.0;
    }
    /* Recalculate the velocity                                      */
    f = (double)(window->e2k[es] - window->e2k[eb] + 1);
    e = es;
    while (f && e != window->zm1e[eb]) {
      windat->u1flux3d[e] += sum / f;
      windat->u1[e] = (windat->dzu1[e]) ? windat->u1flux3d[e] / (windat->dzu1[e] * window->h1au1[e2]) : 0.0;
      e = window->zm1e[e];
    }
    e = es; 
    /* Adjust thin layers if required                                */
    if (wincon->thin_merge && windat->dzu1[e] < wincon->hmin && e != eb) {
      zm1 = window->zm1e[e];
      sum = windat->u1flux3d[e] + windat->u1flux3d[zm1];
      dz = windat->dzu1[e] + windat->dzu1[zm1];
      windat->u1flux3d[e] = windat->u1flux3d[zm1] = 0.5 * sum;
      windat->u1[e] = windat->u1[zm1] = (dz) ? sum / (dz * window->h1au1[e2]) : 0.0;
    }
    /* Set a no-gradient over the surface                           */
    while (e != window->zp1e[e]) {
      e = window->zp1e[e];
      windat->u1[e] = windat->u1[es];
    }
  }
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  set_lateral_BC_vel(windat->u1flux3d, window->nbpte1, window->bpte1, window->bine1);

}

/* END recalc_vel()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Merges volume fluxes of surface layer and layer beneath, and      */
/* resets a uniform velocity over these layers.                      */
/*-------------------------------------------------------------------*/
void merge_volflux(geometry_t *window, /* Processing window      */
		   window_t *windat,   /* Window data structure  */
		   win_priv_t *wincon  /* Window constants       */
		   )
{
  int ee, e, es, eb, zm1;
  double us;

  /* e1 fluxes                                                       */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->sur_e1[ee];
    eb = window->bot_e1[ee];
    es = window->m2d[e];
    zm1 = window->zm1e[e];
    if (e < eb && windat->dzu1[zm1]) {
      /*window->s2k[c] != window->s2k[wincon->i6[cc]]) {*/
      window->sur_e1[ee] = zm1;
      windat->u1vm[zm1] += windat->u1vm[e];
      windat->u1vm[e] = 0.0;
      merge_thin_layers(e, zm1, windat->u1, windat->dzu1);
      us = windat->u1[e];
    }
    /* Set a no-gradient to the surface                              */
    while (e != window->zp1e[e]) {
      windat->u1[e] = us;
      e = window->zp1e[e];
    }
  }
}

/* END merge_volflux()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the velocities given the volume fluxes                    */
/*-------------------------------------------------------------------*/
void flux_to_velocity(geometry_t *window,  /* Processing window      */
		      window_t *windat,    /* Window data structure  */
		      win_priv_t *wincon   /* Window constants       */
		      )
{
  int ee, e, es, eb, e1, e2, zm1;
  double sum, dz; 

  /*-----------------------------------------------------------------*/
  /* Pre-filter the e1 volume fluxes                                 */
  memcpy(windat->u1flux3d, windat->u1vm, window->sze * sizeof(double));
  vel2D_lbc(windat->u1vm, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->sur_e1[ee];
    zm1 = window->zm1e[e];
    while (e != zm1) {
      /*windat->u1flux3d[e] = cvol1w(window, windat->u1vm, e);*/
      /*windat->u1flux3d[e] = shuman(window, windat->u1vm, e);*/
      /*windat->u1flux3d[e] = shapiro_smooth(window, windat->u1vm, e, 0);*/
      /*windat->u1flux3d[e] = con_med(window, windat->u1vm, e);*/
      /*windat->u1av[e] = con_median(window, windat->u1flux3d, 1e10, 
	e, ST_SQ3|ST_EDGE);*/
      e = zm1;
      zm1 = window->zm1e[e];
    }
  }
  memcpy(windat->u1vm, windat->u1flux3d, window->sze * sizeof(double));

  /* Get the cell spacing and new surface coordinate for u1          */
  /*set_dz_at_u1(window, windat, wincon);*/

  /* Recompute the velocities                                        */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = es = e1 = window->sur_e1[ee];
    eb = window->bot_e1[ee];
    e2 = window->m2de[e];
    zm1 = window->zm1e[e];
    while (e != zm1) {
      if (wincon->thin_merge && windat->dzu1[e] < wincon->hmin && e == es && e != eb) {
	sum = windat->u1vm[e] + windat->u1vm[zm1];
	dz = windat->dzu1[e] + windat->dzu1[zm1];
	/*windat->u1vm[e] = windat->u1vm[zm1] = 0.5 * sum;*/
	windat->u1vm[zm1] = sum;
	windat->u1vm[e] = 0.0;
	windat->u1[e] = windat->u1[zm1] = (dz) ? sum / (dz * window->h1au1[e2]) : 0.0;
      } else {
	dz = windat->dzu1[e];
	windat->u1[e] = (dz) ? windat->u1vm[e] / (dz * window->h1au1[e2]) : 0.0;
      }
      e = zm1;
      zm1 = window->zm1e[e];
    }
    while (e1 != window->zp1e[e1]) {
      windat->u1[e1] = windat->u1[es];
      e1 = window->zp1e[e1];
    }
  }
}

/* END flux_to_velocity()                                            */
/*-------------------------------------------------------------------*/
