/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/tracers/transport.c
 *  
 *  Description: Transport model routines
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: transport.c 6720 2021-03-29 00:59:30Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define SMALL 1e-10             /* Minimum value for velocity bounds */
int sls[5] = {8, 8, 27, 64, 125};

int get_pos(geometry_t *window, int c, double nu, double nv, double nw,
            double *cx, double *cy, double *cz, double dt);
int get_pos_bl(geometry_t *window, int c, double nu, double nv, double nw,
	       double *cx, double *cy, double *cz, double dt);
int get_pos_diag(geometry_t *window, int c, double nu, double nv, double nw,
		 double *cx, double *cy, double *cz, double dt);
void global_fill(geometry_t *window, window_t *windat, 
		 win_priv_t *wincon, double *dtracer, int n, int mode);
void monotonic_fill(geometry_t *window, window_t *windat, 
		    win_priv_t *wincon, int n, int mode);
void global_tfill(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, double *dtracer, int n);
double max_lag(geometry_t *window, double *tr, int c);
double min_lag(geometry_t *window, double *tr, int c);
double max_lag_w(geometry_t *window, double *tr, int c, int cs);
double min_lag_w(geometry_t *window, double *tr, int c, int cs);
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
void weights_oo(geometry_t *window, int c, int co, double pin, double qin, 
		double rin, double *dzz);		
void weights_eo(geometry_t *window, int c, int co, double pin, double qin, 
		double rin, double *dz);	     
void set_eo(geometry_t *window, win_priv_t *wincon, double *dzz, int *c, 
	    double *cx, double *cy, double *cz);
int s2t_weights(geometry_t *ws, geometry_t *wt, double *dz, double z, 
		double *wgt, int ct);
int t2rij(geometry_t *window, int c, double p, double q, double *x, double *y);
int t2rc(geometry_t *window, int c, double p, double q);
void print_tridiagonal(geometry_t *window, double *tr, int c, int cs);
void check_tridiagonal(geometry_t *window, double *tr, char *text, int c, int cs);
double sumAij(geometry_t *window, int c);
int iswetAij(geometry_t *window, int cc);
int isnanAij(geometry_t *window, int c);
void local_fill(geometry_t *window, window_t *windat, win_priv_t *wincon,
		int mode);
void set_monotonic(int n, int co, double *ntr, double *tr, int *lmap);
void set_lateral_OBC_tr(geometry_t *window, int ntr, double **tr);
void calc_volume(geometry_t *window, window_t *windat, win_priv_t *wincon);
void check_transport(geometry_t *window, window_t *windat, win_priv_t *wincon);
double aij_1(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag);
double aij_2(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag);
double aij_3(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag);
double aij_4(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag);
void calc_means_t(geometry_t *window, window_t *windat, win_priv_t *wincon);
void recalc_vel(geometry_t *window, window_t *windat, win_priv_t *wincon);
void flux_to_velocity(geometry_t *window, window_t *windat, win_priv_t *wincon);
void merge_volflux(geometry_t *window, window_t *windat, win_priv_t *wincon);
int cg;

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
  reset_means_m(master);

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
    bdry_eval_u2_m(geom, master);
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
    tpg = window[n]->trans->tpg;
    tpd = window[n]->trans->tpd;
    tpc = window[n]->trans->tpc;
    wincon[n]->togn = master->togn;

    /*---------------------------------------------------------------*/
    /* Interpolate onto the target grid (master)                     */
    if(tpc->tgrid & (EXACT|SUBSET)) {
      memcpy(master->eta, tpd->eta, tpg->sgsizS * sizeof(double));
      memcpy(master->Kz, tpd->Kz, tpg->sgsiz * sizeof(double));
      /* Velocities are required for MONOTONIC global fills          */
      if (master->fillf & MONOTONIC || master->smagorinsky > 0.0) {
	memcpy(master->u1, tpd->u1, tpg->sgsiz * sizeof(double));
	memcpy(master->u2, tpd->u2, tpg->sgsiz * sizeof(double));
      }
    } else if (tpc->tgrid & INEXACT) {
      int c2, cs, cg, tn, m;
      double *dzz = tpc->w9;
      double *wgt = wincon[1]->wgt[1];
      double fi,fj;

      /* Set up arrays                                               */
      set_dz(tpg, tpd, tpc);
      set_dzz(tpg, dzz);
      /* Interpolate onto the target grid (master).                  */
      /* Set the lateral boundary conditions on the source grid.     */
      set_lateral_bc_eta(tpd->eta, tpg->nbptS, tpg->bpt, tpg->bin, tpg->bin2, 1);
      vel2D_lbc(tpd->u1, tpg->nbpte1, tpg->nbe1, tpg->bpte1, tpg->bine1, tpc->slip);
      vel2D_lbc(tpd->u2, tpg->nbpte2, tpg->nbe2, tpg->bpte2, tpg->bine2, tpc->slip);
      vel2D_lbc(tpd->Kz, 0, tpg->nbpt, tpg->bpt, tpg->bin, 1.0);
      set_lateral_OBC_tr(tpg, tpd->ntr, tpd->tr_wc);
      set_lateral_BC_tr(tpd->tr_wc, tpd->ntr, tpg->nbpt, tpg->bpt, tpg->bin);

      /* Set the bottom boundary condition (no-gradient)             */
      for (cc = 1; cc <= tpg->a2_t; cc++) {
	c = tpg->bot_t[cc];
	c2 = tpg->zm1[c];
	tpd->u1[c2] = -tpd->u1[c];
	tpd->u2[c2] = -tpd->u2[c];
	tpd->Kz[c2] = tpd->Kz[c];
	for (tn = 0; tn < tpd->ntr; tn++)
	  tpd->tr_wc[tn][c2] = tpd->tr_wc[tn][c];
      }

      for (cc = 1; cc <= tpg->ngsed; cc++) {
	c = tpg->gsed_t[cc];
	c2 = tpg->zp1[c];
	tpd->u1[c] = -tpd->u1[c2];
	tpd->u2[c] = -tpd->u2[c2];
	tpd->Kz[c] = tpd->Kz[c2];
	for (tn = 0; tn < tpd->ntr; tn++)
	  tpd->tr_wc[tn][c] = tpd->tr_wc[tn][c2];
      }

      /* Interpolate                                                 */
      for (cc = 1; cc <= window[n]->b2_t; cc++) {
	c = window[n]->w2_t[cc];
	cs = s2t_weights(tpg, window[n], dzz, 0.0, wgt, c);
	cg = window[n]->wsa[c];
	master->eta[cg] = int_val_bl(tpg, cs, wgt, tpd->eta);
	/* Interpolate additional tracers                            */
	for (tn = 0; tn < master->ntrvarsS; tn++) {
	  m = master->trvmS[tn];
	  master->tr_wcS[m][cg] = int_val_bl(tpg, cs, wgt, tpd->tr_wcS[tn]);
	}
      }
      /* Rotate source vectors to cell centered east/north */
      rotate_vel(tpg, tpc->w4, tpc->w5);
      /* Interpolate */
      for (cc = 1; cc <= window[n]->b3_t; cc++) {
	double z;
	c = window[n]->w3_t[cc];
	c2 = window[n]->m2d[c];
	z = max(window[n]->cellz[c], window[n]->botz[c2]);
	cs = s2t_weights(tpg, window[n], dzz, z, wgt, c);
	cg = window[n]->wsa[c];
	if (cs) {
	  master->Kz[cg] = int_val_bl(tpg, cs, wgt, tpd->Kz);
	  /*
	  master->u1[cg] = int_val_bl(tpg, cs, wgt, tpd->u1);
	  master->u2[cg] = int_val_bl(tpg, cs, wgt, tpd->u2);
	  */
	  /* Need to rorate these vectors? */
	  wincon[n]->w4[cg] = int_val_bl(tpg, cs, wgt, tpc->w4);
	  wincon[n]->w5[cg] = int_val_bl(tpg, cs, wgt, tpc->w5);
	  /* Interpolate T/S if required. Note these cannot be read  */
	  /* onto the master in reset_trans_event() since transport  */
	  /* data (in sparse format) is for the source grid only.    */
	  if (!master->trinfo_3d[master->tno].advect && 
	      !master->trinfo_3d[master->tno].diffuse)
	    master->temp[cg] = int_val_bl(tpg, cs, wgt, tpd->temp);
	  if (!master->trinfo_3d[master->sno].advect && 
	      !master->trinfo_3d[master->sno].diffuse)
	    master->sal[cg] = int_val_bl(tpg, cs, wgt, tpd->sal);
	  /* Interpolate additional tracers (tracers 0 & 1 are salt  */
	  /* & temp)                                                 */
	  for (tn = 0; tn < master->ntrvars; tn++) {
	    m = master->trvm[tn];
	    master->tr_wc[m][cg] = int_val_bl(tpg, cs, wgt, tpd->tr_wc[2+tn]);
	  }
	}
      }
      /* Rotate target vectors onto the grid */
      rerotate_vel(window[n], wincon[n]->w4, wincon[n]->w5);
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
      vel2D_lbc(windat[n]->u2, window[n]->nbpte2, window[n]->nbe2,
		window[n]->bpte2, window[n]->bine2, wincon[n]->slip);
      alerts_w(window[n], VEL3D);
    }

    /*---------------------------------------------------------------*/
    /* Read the boundary velocities from file for global fills if    */
    /* required.                                                     */
    if (master->tmode & DO_OBC) {
      for (nb = 0; nb < window[n]->nobc; nb++)
	reset_bdry_eta(window[n], windat[n], wincon[n], window[n]->open[nb], windat[n]->eta);
      windat[n]->nu1 = wincon[n]->w9;
      windat[n]->nu2 = wincon[n]->w10;
      /* Transfer any custom data from the master to the slaves      */
      bdry_transfer_u1(master, window[n], windat[n]);
      bdry_transfer_u2(master, window[n], windat[n]);
      memcpy(windat[n]->nu1, windat[n]->u1, window[n]->sgsiz * sizeof(double));
      memcpy(windat[n]->nu2, windat[n]->u2, window[n]->sgsiz * sizeof(double));
      bdry_u1_3d(window[n], windat[n], wincon[n]);
      bdry_u2_3d(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->u1, windat[n]->nu1, window[n]->sgsiz * sizeof(double));
      memcpy(windat[n]->u2, windat[n]->nu2, window[n]->sgsiz * sizeof(double));
    }
    if (master->fillf & SET_BDRY) {
      windat[n]->nu1 = wincon[n]->w9;
      windat[n]->nu2 = wincon[n]->w10;
      /* Transfer any custom data from the master to the slaves      */
      bdry_transfer_u1(master, window[n], windat[n]);
      bdry_transfer_u2(master, window[n], windat[n]);
      memset(windat[n]->nu1, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->nu2, 0, window[n]->sgsiz * sizeof(double));
      bdry_u1_3d(window[n], windat[n], wincon[n]);
      bdry_u2_3d(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->u1, windat[n]->nu1, window[n]->sgsiz * sizeof(double));
      memcpy(windat[n]->u2, windat[n]->nu2, window[n]->sgsiz * sizeof(double));
    }
    if (master->fillf & SET_BDRY_ETA) {
      wincon[n]->neweta = wincon[n]->d4;
      bdry_transfer_eta(master, window[n], windat[n]);
      memset(wincon[n]->neweta, 0, window[n]->sgsizS * sizeof(double));
      bdry_eta(window[n], windat[n], wincon[n]);
      memcpy(windat[n]->eta, wincon[n]->neweta, window[n]->sgsizS * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.  */
#if !GLOB_BC
    set_lateral_OBC_tr(window[n], windat[n]->ntr, windat[n]->tr_wc);
    set_lateral_BC_tr(windat[n]->tr_wc, windat[n]->ntr, window[n]->nbpt,
                      window[n]->bpt, window[n]->bin);
    set_lateral_BC_vel(windat[n]->u1flux3d, window[n]->nbpt,
                       window[n]->bpte1, window[n]->bine1);
    set_lateral_BC_vel(windat[n]->u2flux3d, window[n]->nbpt,
                       window[n]->bpte2, window[n]->bine2);
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
      memset(windat[n]->u1, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->u2, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->w, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->u1vm, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->u2vm, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->u1flux3d, 0, window[n]->sgsiz * sizeof(double));
      memset(windat[n]->u2flux3d, 0, window[n]->sgsiz * sizeof(double));
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
    /* Reset occurs @ TS at the end of tracer_step_3d()              */
    /*reset_means(window[n], windat[n], wincon[n], RESET);*/

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
    /* keyword = GLOB+VERR, flag = LOCALER; global fill using local    */
    /*           volume errors to compute global scaling factor.       */
    /* keyword = LOCAL, flag = LOCAL; local fill                       */
    /* keyword = MONOTONIC, flag |= MONOTONIC; local fill bails at     */
    /*           perfect mass balance, then less restrictive monotonic */
    /*           contraint is imposed.                                 */
    /* keyword = MONO+GLOB, flag |= MONGLOB|MONOTONIC; same as         */
    /*           monotonic, but residual mass is iteratively spread    */
    /*           globally.                                             */
    /* keyword = DIAGNOSE, flag |= DIAGNOSE; print diagnostics to file */
    /*           trans.ts                                              */
    /* keyword = DIAGNOSE_BGC, flag |= DIAGNOSE_BGC; print diagnostics */
    /*           to file trans.ts                                      */
    if (!(wincon[n]->fillf & (NONE|LOCAL)))
      global_fill(window[n], windat[n], wincon[n], NULL, -1, 0);
    if (wincon[n]->fillf & (LOCAL|LOCALER))
      local_fill(window[n], windat[n], wincon[n], 1);

    /*---------------------------------------------------------------*/
    /* Set up non semi-lagrange advection schemes                    */
    if (!(wincon[n]->trasc & LAGRANGE)) {
      if (wincon[n]->conserve & CONS_MRG)
	merge_volflux(window[n], windat[n], wincon[n]);

      /* Set the vertical grid spacings at e1 and e2 faces.          */
      set_dz_at_u1(window[n], windat[n], wincon[n]);
      set_dz_at_u2(window[n], windat[n], wincon[n]);
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
      set_dz_at_u2(window[n], windat[n], wincon[n]);
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
      memcpy(windat[n]->u2b, windat[n]->u2, window[n]->sgsiz * sizeof(double));
      wincon[n]->hor_mix->pre(window[n], windat[n], wincon[n]);
      wincon[n]->hor_mix->setup(window[n], windat[n], wincon[n]);
    }

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
  if(wincon->tgrid & NONE)
    semi_lagrange_c(window, windat, wincon);
  else
    semi_lagrange_t(window, wincon, tpg, tpd, tpc);

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
  int nn;                       /* Tracer index */
  int c, cc;                    /* Sparse indices */
  int c2;                       /* 2D cell corresponding to 3D cell */
  int cb;                       /* Bottom sparse coordinate (cell centre) */
  int vc;                       /* Tracer cells to provess counter */
  int vcs;                      /* Surface tracer cells to process counter */                                 
  int vc1;                      /* As for vc excluding OBC cells */
  int vcs1;                     /* As for vcs excluding OBC cells */
  int zp1;                      /* k index at k+1 */
  double *Fx;                   /* x tracer flux */
  double *Fy;                   /* y tracer flux */
  double *Fz;                   /* z tracer flux */
  double *dtracer;              /* Advective terms for x,y directions */
  double *hflux;                /* Horizontal diffusion flux */
  double dtu;                   /* Sub-time step to use */
  int slf = 1;                  /* Set to zero on the first sub-step */
  int hdif = 0;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* w1 = dtracer */
  /* w2 = Fx */
  /* w3 = Fy */
  /* w4 = Fz */
  /* d3 = mean cell spacing in the x direction */
  /* d4 = mean cell spacing in the y direction */
  /* s1 = wet cells to process for tracers */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
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
  if (wincon->u1kh0 > 0.0 && wincon->u2kh0 > 0.0) hdif = 1;

  /* Find the origin of the streamline & grid Courant numbers */
  if (wincon->tmode & SP_ORIGIN) {
    get_weights(window, windat, wincon);
  } else {
    if(wincon->tgrid & NONE)
      semi_lagrange_c(window, windat, wincon);
    else if(wincon->tgrid & (EXACT|SUBSET))
      semi_lagrange_t(window, wincon, window->trans->tpg, 
		      window->trans->tpd, window->trans->tpc);
    else
      semi_lagrange_tu(window, wincon, window->trans->tpg, 
		       window->trans->tpd, window->trans->tpc);
  }
  if (wincon->tmode & TR_CHECK)
    check_transport(window, windat, wincon);

  /* Reset the interpolation weights and get volume errors           */
  if (wincon->tmode & SET_AIJ) {
    reset_Aij(window, windat, wincon);
    if (wincon->fillf & LOCAL)
      local_fill(window, windat, wincon, 0);
    get_verr(window, windat, wincon);
  }

  /* No ultimate limiting */
  wincon->ultimate = 0;
  /* Time step */
  dtu = windat->dttr;

  /*-----------------------------------------------------------------*/
  /* Set a no-gradient OBC for vertical velocity; since u & v OBCs   */
  /* are applied independently, the divergence in w boundary cells   */
  /* does not obey continuity, hence w computed in boundary cells    */
  /* are erroneous. These are, however used in the lagrange scheme.  */
  bdry_w(window, windat, wincon);

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

    /*---------------------------------------------------------------*/
    /* Evaluate the sources and sinks before advection               */
    if (wincon->npss) {
      ss_tracer(window, windat, wincon, n, dtracer, dtu);
      if (!(wincon->fillf & (NONE|LOCAL)))
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
	    if (wincon->fillf & LOCAL) Vtm1 -= Vpss;   /* Vol at t-1 */
	    /* Total volume. For local fills the contribution from   */
	    /* all destination cells is used.                        */
	    /*Vtot = (wincon->fillf & LOCAL) ? Vi[c] : (Vtm1 + Vpss);*/
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
      calc_flux_old(window, windat, wincon, dtu);

    /*---------------------------------------------------------------*/
    /* Get the horizontal diffusive fluxes                           */
    if (wincon->diffuse[n] && hdif) {
      hor_diffuse(window, windat, wincon, tr, Fx, Fy);
      memset(hflux, 0, window->sgsiz * sizeof(double));
      for (cc = 1; cc <= window->v3_t; cc++) {
	c = window->w3_t[cc];
	hflux[c] = (Fx[window->xp1[c]] - Fx[c] +
		    Fy[window->yp1[c]] - Fy[c]) * dtu;
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
    if (!(wincon->fillf & (NONE|LOCAL|LOCALER))) {
      int mode = (nn == 0) ? 2 : 3;
      mode = (nn == wincon->ntbdy - 1) ? 4 : mode;
      mode = (wincon->ntbdy == 1) ? 6 : mode;
      global_fill(window, windat, wincon, dtracer, n, mode);
    }
    if ((wincon->fillf & LOCAL && wincon->fillf & MONOTONIC) || 
	wincon->fillf & LOCALER) {
      int mode = (nn == 0) ? 1 : 0;
      mode = (nn == wincon->ntbdy - 1) ? 2 : mode;
      monotonic_fill(window, windat, wincon, n, mode);
    }
  }                         /* Tracer loop end */

  /* Clip the tracer concentration if required */
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

  /*-----------------------------------------------------------------*/
  /* Reset layer thickness if newly wetted cells were merged         */
  if (wincon->fillf & (LOCAL|LOCALER) && wincon->tmode & SET_AIJ)
    reset_merged_dz(window, windat, wincon);
}

/* END advect_diffuse_lag()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the grid     */
/* Courant numbers. Used for single grids (source =  target).        */
/*-------------------------------------------------------------------*/
void semi_lagrange_c(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->w4;
  double *nv = wincon->w5;
  double *nw = wincon->w10;
  double *cx = wincon->w6;
  double *cy = wincon->w7;
  double *cz = wincon->w8;
  double *dzz = wincon->w9;
  double *mask = wincon->d3;
  double u, v, w;
  double p, q, r;
  int cc, c, ci, c2, zm1;
  int *cl = wincon->s2;        /* Streamline origin */
  int *c2cc = wincon->s4;      /* c to cc map */
  int n;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(cx, 0, window->sgsiz * sizeof(double));
  memset(cy, 0, window->sgsiz * sizeof(double));
  memset(cz, 0, window->sgsiz * sizeof(double));
  memset(nu, 0, window->sgsiz * sizeof(double));
  memset(nv, 0, window->sgsiz * sizeof(double));
  memset(nw, 0, window->sgsiz * sizeof(double));
  memset(cl, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing */
  set_dzz(window, dzz);

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face / center */
  vel_center(window, windat, wincon, nu, nv, nw);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    nu[c] = -nu[ci];
    nv[c] = -nv[ci];
    nw[c] = -nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
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
  /* Set the cell width  mask */
  for (c = 1; c <= window->enonS; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    mask[c] = 0.5;
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = cg = cl[cc] = wincon->s1[cc];
    cg = c;
    time_left = windat->dttr;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    n=0;

    while (time_left > 0) {

      int ins = (cl[cc] == wincon->m2d[cl[cc]]) ? 1 : 0;
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      c2 = window->m2d[cl[cc]];
      ci = (nu[cl[cc]] > 0.0) ? c2 : window->xp1[c2];
      p = window->h1au1[ci] * mask[ci];
      ci = (nv[cl[cc]] > 0.0) ? c2 : window->yp1[c2];
      q = window->h2au2[ci] * mask[ci];
      r = (nw[cl[cc]] < 0.0) ? dzz[cl[cc]] : dzz[window->zm1[cl[cc]]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = (r) ? ((ins && nw[cl[cc]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[cc]]))) : dt;
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position.  */
      cl[cc] =
        get_pos(window, cl[cc], u, v, w, &cx[c], &cy[c], &cz[c], dt);

      c2 = window->m2d[cl[cc]];
      p = cx[c] / window->h1au1[c2];
      q = cy[c] / window->h2au2[c2];
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;

      /* Since the velocity field may be non-conservative, it is     */
      /* possible that the streamline may arrive above the surface   */
      /* or below the sea bed, hence if the origin is in the top or  */
      /* bottom layer, the distance should be limited to eta-cellz   */
      /* or cellz-botz (i.e. dzz). Limiting 0<r<1 is equivalent to   */
      /* this.                                                       */
      /* This is currently controlled with a limit on *cz in         */
      /* get_pos().                                                  */
      /*r = min(max(0.0, r), 1.0);*/

      wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - r) * (1.0 - p) * q;
      wgt[2] = (1.0 - r) * p * (1.0 - q);
      wgt[3] = (1.0 - r) * p * q;
      wgt[4] = r * (1.0 - p) * (1.0 - q);
      wgt[5] = r * (1.0 - p) * q;
      wgt[6] = r * p * (1.0 - q);
      wgt[7] = r * p * q;

      /* Interpolate x,y,z velocities */
      u = int_val(window, cl[cc], wgt, nu);
      v = int_val(window, cl[cc], wgt, nv);
      w = int_val(window, cl[cc], wgt, nw);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      n++;

      /*
      if(p > 1.0 || p < 0.0)
	hd_warn("x cell distance outside of bounds at %d : %f\n", c, p);
      if(q > 1.0 || q < 0.0)
	hd_warn("y cell distance outside of bounds at %d : %f\n", c, q);
      if(r > 1.0 || r < 0.0)
	hd_warn("r cell distance outside of bounds at %d : %f\n", c, r);
      */
    }

    sl_check(window, &cl[cc], &cx[c], &cy[c], &cz[c]);
    if (wincon->osl == 0) {

      c2 = window->m2d[cl[cc]];
      p = cx[c] / window->h1au1[c2];
      q = cy[c] / window->h2au2[c2];
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;
      /*r = min(max(0.0, r), 1.0);*/

      /* Get the interpolation weights */
      wincon->wgt[c][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][1] = (1.0 - r) * (1.0 - p) * q;
      wincon->wgt[c][2] = (1.0 - r) * p * (1.0 - q);
      wincon->wgt[c][3] = (1.0 - r) * p * q;
      wincon->wgt[c][4] = r * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][5] = r * (1.0 - p) * q;
      wincon->wgt[c][6] = r * p * (1.0 - q);
      wincon->wgt[c][7] = r * p * q;
    } else {
      if (wincon->osl == 2 || wincon->osl == 4) {
	set_eo(window, wincon, dzz, &cl[cc], &cx[c], &cy[c], &cz[c]);
	weights_eo(window, c, cl[cc], cx[c], cy[c], cz[c], dzz);
      } else
	weights_oo(window, c, cl[cc], cx[c], cy[c], cz[c], dzz);
    }

    /* Copy the streamline origin and offsets into t11, t12 and t22  */
    /* respectively, so that these quantities are available for      */
    /* dumping to sparse output files. Note that these variables may */
    /* only be dumped to sparse file since interpolation from        */
    /* (x,y,z) format is not possible for the sparse coordinate of   */
    /* the streamline origin.                                        */
    if (windat->origin) {
      windat->origin[c] = (double)cl[cc];
      windat->pc[c] = p;
      windat->qc[c] = q;
      windat->rc[c] = r;
    }
  }
  debug_c(window, D_TS, D_STRML);
}

/* END semi_lagrange_c()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the grid     */
/* Courant numbers. Used for equal|subset source and target grids.   */
/* This uses the origin for interpolation as the top right corner    */
/* of the transport grid. This has the advantage that transport grid */
/*-------------------------------------------------------------------*/
void semi_lagrange_t(geometry_t *window,  /* Target geometry         */
                     win_priv_t *wincon,  /* Target constants        */
		     geometry_t *tpg,     /* Source geometry         */
                     window_t *tpd,       /* Source data             */
                     win_priv_t *tpc      /* Source constants        */
  )
{
  window_t *windat = window->windat;
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = tpc->w4;
  double *nv = tpc->w5;
  double *nw = tpc->w10;
  double *cx = wincon->w6;
  double *cy = wincon->w7;
  double *cz = wincon->w8;
  double *dzt = wincon->w9;
  double *dzz = tpc->w9;
  double *xinit = wincon->d7;
  double *yinit = wincon->d8;
  double *xoset = tpc->d7;
  double *yoset = tpc->d8;
  double *mask = tpc->d3;
  double u, v, w;
  double p, q, r;
  double fi, fj, x, y, z, zb;
  int cc, c, ci, c2, cs, ct, ctc, zm1;
  double *clt = wincon->w10;
  int *cl = wincon->s2;
  int *c2cc = wincon->s4;
  int *t2s = window->mgm;
  int *s2t = tpg->mgm;
  int ins;                     /* Is streamline in surface cell? */
  int n;

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across tracer OBC's             */
  set_map_t(tpg);

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(nu, 0, tpg->sgsiz * sizeof(double));
  memset(nv, 0, tpg->sgsiz * sizeof(double));
  memset(nw, 0, tpg->sgsiz * sizeof(double));
  memset(cl, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing, source grid                      */
  set_dzz(tpg, dzz);

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing, target grid                      */
  set_dzz(window, dzt);
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    dzt[c] = dzt[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face / center */
  vel_center(tpg, tpd, tpc, nu, nv, nw);

  /* Set the ghost cells */
  for (cc = 1; cc <= tpg->nbpt; cc++) {
    c = tpg->bpt[cc];
    ci = tpg->bin[cc];
    nu[c] = -nu[ci];
    nv[c] = -nv[ci];
    nw[c] = -nw[ci];
    dzz[c] = dzz[ci];
    tpc->m2d[c] = (ci == tpc->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
  for (cc = 1; cc <= tpg->b2_t; cc++) {
    c = tpg->bot_t[cc];
    zm1 = tpg->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get the increments of the target grid in the source grid        */ 
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];  /* Valid wet cell in target grid            */
    c2 = window->m2d[c]; /* Surface cell correspoding to c           */
    cs = t2s[c];         /* Correspoding cell in source grid         */
    cx[c] = xinit[c2];
    cy[c] = yinit[c2];
    cz[c] = (window->cellz[c] - tpg->cellz[cs]) / dzz[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the cell width  mask                                        */
  for (c = 1; c <= tpg->enonS; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= tpg->nbptS; cc++) {
    c = tpg->bpt[cc];
    mask[c] = 0.5;
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin. Loop over valid wet    */
  /* cells in the target grid (wincon->s1) and find the              */
  /* corresponding cell in the source grid (window->mgm(wincon->s1)).*/
  for (cc = 1; cc <= wincon->vc; cc++) {
    ct = wincon->s1[cc];         /* Target grid wet cell             */
    c2cc[ct] = cc;               /* Target grid index map            */
    c = cl[cc] = t2s[ct];        /* Source grid wet cell             */ 
    clt[cc] = (double)c;

    time_left = tpd->dttr;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    if(!c) continue;

    while (time_left > 0) {

      int ins = (cl[cc] == tpc->m2d[cl[cc]]) ? 1 : 0;
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      c2 = tpg->m2d[cl[cc]];
      ci = (nu[cl[cc]] > 0.0) ? c2 : tpg->xp1[c2];
      p = tpg->h1au1[ci] * mask[ci];
      ci = (nv[cl[cc]] > 0.0) ? c2 : tpg->yp1[c2];
      q = tpg->h2au2[ci] * mask[ci];
      r = (nw[cl[cc]] < 0.0) ? dzz[c] : dzz[tpg->zm1[c]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = (r) ? ((ins && nw[cl[cc]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[cc]]))) : dt;
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position in the source grid.  */
      cl[cc] =
        get_pos(tpg, cl[cc], u, v, w, &cx[ct], &cy[ct], &cz[ct], dt);

      /* Get the interpolation weights */
      c2 = tpg->m2d[cl[cc]];

      p = cx[ct] / tpg->h1au1[c2];
      q = cy[ct] / tpg->h2au2[c2];
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;

      wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - r) * (1.0 - p) * q;
      wgt[2] = (1.0 - r) * p * (1.0 - q);
      wgt[3] = (1.0 - r) * p * q;
      wgt[4] = r * (1.0 - p) * (1.0 - q);
      wgt[5] = r * (1.0 - p) * q;
      wgt[6] = r * p * (1.0 - q);
      wgt[7] = r * p * q;

      /* Interpolate x,y,z velocities */
      u = int_val(tpg, cl[cc], wgt, nu);
      v = int_val(tpg, cl[cc], wgt, nv);
      w = int_val(tpg, cl[cc], wgt, nw);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;

      /* Update the target streamline for valid source cells only */
      clt[cc] = (s2t[cl[cc]]) ? (double)cl[cc] : clt[cc];

      /*
      if(p > 1.0 || p < 0.0)
	hd_warn("Source x cell distance outside of bounds at %d : %f\n", c, p);
      if(q > 1.0 || q < 0.0)
	hd_warn("Source y cell distance outside of bounds at %d : %f\n", c, q);
      if(r > 1.0 || r < 0.0)
	hd_warn("Source r cell distance outside of bounds at %d : %f\n", c, r);
      */
    }
  }

  /* Get the weights in the target grid. Note: cl[cc] is the sparse  */
  /* location of the origin of the streamline in the source grid.    */
  /* tpg->mgm maps this to the corresponding target location.        */
  memset(windat->pc, 0, window->sgsiz * sizeof(double));
  memset(windat->qc, 0, window->sgsiz * sizeof(double));
  memset(windat->rc, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    ct = wincon->s1[cc];    /* Target wet cell                     */
    cs = (int)clt[cc];      /* Source cell streamline origin       */
    ci = tpg->m2d[cs];      /* 2D source cell location             */
    c = s2t[cs];            /* Target cell streamline origin       */
    c2 = window->m2d[c];    /* 2D target cell location             */
    
    p = (cx[ct] + xoset[ci]) / window->h1au1[c2];
    q = (cy[ct] + yoset[ci]) / window->h2au2[c2];
    r = cz[ct] + tpg->cellz[cs] - window->cellz[c];
    r = (r < dzt[c]) ? r / dzt[c] : (r - dzt[c]) / dzt[window->zp1[c]];

    wincon->wgt[ct][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
    wincon->wgt[ct][1] = (1.0 - r) * (1.0 - p) * q;
    wincon->wgt[ct][2] = (1.0 - r) * p * (1.0 - q);
    wincon->wgt[ct][3] = (1.0 - r) * p * q;
    wincon->wgt[ct][4] = r * (1.0 - p) * (1.0 - q);
    wincon->wgt[ct][5] = r * (1.0 - p) * q;
    wincon->wgt[ct][6] = r * p * (1.0 - q);
    wincon->wgt[ct][7] = r * p * q;
    cl[cc] = c;
    /* Copy the streamline origin and offsets into origin, p,q and */
    /* r respectively, so that these quantities are available to   */
    /* dumping to sparse output files. Note that these variables   */
    /* may only be dumped to sparse file since interpolation from  */
    /* (x,y,z) format is not possible for the sparse coordinate of */
    /* the streamline origin.                                      */
    windat->origin[ct] = (double)cl[cc];
    windat->pc[ct] = p;
    windat->qc[ct] = q;
    windat->rc[ct] = r;
  }
  debug_c(window, D_TS, D_STRML);
}

/* END semi_lagrange_t()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the grid     */
/* Courant numbers. Used for inexact source and target grids.        */
/* This uses the origin for interpolation as the bottom left corner  */
/* of the transport grid. This has the advantage that transport grid */
/* locations (i.e. using cell centres for the mesh) are the same as  */
/* the regular grid locations (which uses grid corners for the mesh) */
/* but the disadvantage that transport grid cell dimensions are      */
/* given at xp1[c] and yp1[c] rather than at c.                      */
/*-------------------------------------------------------------------*/
void semi_lagrange_tu(geometry_t *window,  /* Target geometry        */
		      win_priv_t *wincon,  /* Target constants       */
		      geometry_t *tpg,     /* Source geometry        */
		      window_t *tpd,       /* Source data            */
		      win_priv_t *tpc      /* Source constants       */
  )
{
  window_t *windat = window->windat;
  transport_t *tp = window->trans;
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = tpc->w4;
  double *nv = tpc->w5;
  double *nw = tpc->w10;
  double *cx = wincon->w6;
  double *cy = wincon->w7;
  double *cz = wincon->w8;
  double *dzt = wincon->w9;
  double *dzz = tpc->w9;
  double *xinit = wincon->d7;
  double *yinit = wincon->d8;
  double *xoset = tpc->d7;
  double *yoset = tpc->d8;
  double *mask = tpc->d3;
  double u, v, w;
  double p, q, r;
  double fi, fj, x, y, z, zb;
  int cc, c, ci, c2, cs, ct, ctc, zm1;
  double *clt = wincon->w10;
  int *cl = wincon->s2;
  int *c2cc = wincon->s4;
  int *t2s = window->mgm;
  int *s2t = tpg->mgm;
  int ins;                     /* Is streamline in surface cell? */
  int n;

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across tracer OBC's             */
  set_map_t(tpg);

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(nu, 0, tpg->sgsiz * sizeof(double));
  memset(nv, 0, tpg->sgsiz * sizeof(double));
  memset(nw, 0, tpg->sgsiz * sizeof(double));
  memset(cl, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing, source grid                      */
  set_dzz(tpg, dzz);

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing, target grid                      */
  set_dzz(window, dzt);
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    dzt[c] = dzt[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face / center */
  vel_center(tpg, tpd, tpc, nu, nv, nw);

  /* Set the ghost cells */
  for (cc = 1; cc <= tpg->nbpt; cc++) {
    c = tpg->bpt[cc];
    ci = tpg->bin[cc];
    nu[c] = -nu[ci];
    nv[c] = -nv[ci];
    nw[c] = -nw[ci];
    dzz[c] = dzz[ci];
    tpc->m2d[c] = (ci == tpc->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
  for (cc = 1; cc <= tpg->b2_t; cc++) {
    c = tpg->bot_t[cc];
    zm1 = tpg->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get the increments of the target grid in the source grid        */ 
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];  /* Valid wet cell in target grid            */
    c2 = window->m2d[c]; /* Surface cell correspoding to c           */
    cs = t2s[c];         /* Correspoding cell in source grid         */
    ci = tpg->m2d[cs];   /* 2D source cell                           */
    cx[c] = xinit[c2] * tpg->h1au1[tpg->xp1[ci]];
    cy[c] = yinit[c2] * tpg->h2au2[tpg->yp1[ci]];
    cz[c] = (window->cellz[c] - tpg->cellz[cs]) / dzz[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the cell width  mask                                        */
  for (c = 1; c <= tpg->enonS; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= tpg->nbptS; cc++) {
    c = tpg->bpt[cc];
    mask[c] = 0.5;
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin. Loop over valid wet    */
  /* cells in the target grid (wincon->s1) and find the              */
  /* corresponding cell in the source grid (window->mgm(wincon->s1)).*/
  for (cc = 1; cc <= wincon->vc; cc++) {
    ct = wincon->s1[cc];         /* Target grid wet cell             */
    c2cc[ct] = cc;               /* Target grid index map            */
    c = cl[cc] = t2s[ct];        /* Source grid wet cell             */ 
    clt[cc] = (double)c;

    time_left = tpd->dttr;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    if(!c) continue;

    while (time_left > 0) {

      int ins = (cl[cc] == tpc->m2d[cl[cc]]) ? 1 : 0;
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      c2 = tpg->m2d[cl[cc]];
      ci = (nu[cl[cc]] > 0.0) ? c2 : tpg->xp1[c2];
      p = tpg->h1au1[ci] * mask[ci];
      ci = (nv[cl[cc]] > 0.0) ? c2 : tpg->yp1[c2];
      q = tpg->h2au2[ci] * mask[ci];
      r = (nw[cl[cc]] < 0.0) ? dzz[c] : dzz[tpg->zm1[c]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = (r) ? ((ins && nw[cl[cc]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[cc]]))) : dt;
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position in the source grid.  */
      cl[cc] =
        get_pos_bl(tpg, cl[cc], u, v, w, &cx[ct], &cy[ct], &cz[ct], dt);

      /* Get the interpolation weights */
      c2 = tpg->m2d[cl[cc]];
      p = cx[ct] / tpg->h1au1[tpg->xp1[c2]];
      q = cy[ct] / tpg->h2au2[tpg->yp1[c2]];
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;

      wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - r) * (1.0 - p) * q;
      wgt[2] = (1.0 - r) * p * (1.0 - q);
      wgt[3] = (1.0 - r) * p * q;
      wgt[4] = r * (1.0 - p) * (1.0 - q);
      wgt[5] = r * (1.0 - p) * q;
      wgt[6] = r * p * (1.0 - q);
      wgt[7] = r * p * q;

      /* Interpolate x,y,z velocities */
      u = int_val_bl(tpg, cl[cc], wgt, nu);
      v = int_val_bl(tpg, cl[cc], wgt, nv);
      w = int_val_bl(tpg, cl[cc], wgt, nw);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      /* Convert to the regular grid so cl is always wet */
      ctc = t2rc(tpg, cl[cc], p, q);
      /* Update the target streamline for valid source cells only */
      clt[cc] = (s2t[ctc]) ? (double)cl[cc] : clt[cc];

      /*
      if(p > 1.0 || p < 0.0)
	hd_warn("Source x cell distance outside of bounds at %d : %f\n", c, p);
      if(q > 1.0 || q < 0.0)
	hd_warn("Source y cell distance outside of bounds at %d : %f\n", c, q);
      if(r > 1.0 || r < 0.0)
	hd_warn("Source r cell distance outside of bounds at %d : %f\n", c, r);
      */
    }
    cx[ct] = p;
    cy[ct] = q;
  }

  /*-----------------------------------------------------------------*/
  /* Convert cl[cc], cx, cy, cz to (lat,long) using the source       */
  /* tp->xyij_tree, then interpolate onto the target grid.           */
  memset(windat->pc, 0, window->sgsiz * sizeof(double));
  memset(windat->qc, 0, window->sgsiz * sizeof(double));
  memset(windat->rc, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    ct = wincon->s1[cc];    /* Target wet cell                       */
    cs = (int)clt[cc];      /* Source cell streamline origin         */
    cl[cc] = ct;

    memset(wincon->wgt[ct], 0, 8 * sizeof(double));
    wincon->wgt[ct][0] = 1.0;
    if(!cs) continue;

    /* Get the lat/long/depth of the source streamline               */
    c = t2rij(tpg, cs, cx[ct], cy[ct], &fi, &fj); /* c is always wet */
    c2 = tpg->m2d[c];
    /* Get the depth relative to mean sea level                      */
    z = cz[ct] + max(tpg->cellz[c], tpg->botz[c2]);
    /* Long / lat of the streamline origin                           */
    grid_fgrid_ijtoxy(tpg->xyij_tree, fi, fj, &x, &y);
    /* Get new value at ct, trn = interp(gs, tr, x, y, z) */
    /* Get the sparse location and increments in the target grid     */
    c2 = xytoc(geom, master->xyij_tree, x, y, &fi, &fj);
    c2 = (c2) ? c2 : window->m2d[s2t[cs]];
    /* Can't find location in target grid -> continue                */
    if(c2==0) {
      if (!(cs = find_nearest(window, x, y))) {
	hd_warn("Lost streamline at target coordinate %d\n",ct);
	continue;
      } else
	hd_warn("Can't map streamline from %d; using origin at %d\n",ct,c2);
    }
    /* Truncate depth to target cell depth                           */
    z = max(z, window->botz[c2]);
    /* Get the depth of the streamline in the target grid            */
    c = ztoc_z(window, c2, z, window->cellz);
    /* Get the vertical increment                                    */
    /* Note: if eta is not the same in corresponding source and      */
    /* target locations (e.g. due to interpolation errors) then dzz  */
    /* will differ in the surface and layer below the surface. This  */
    /* can cause r > 1.                                              */
    r = (z - max(window->cellz[c], window->botz[c2])) / dzt[c];
    /* Transform to the transport grid                               */
    fi -= floor(fi);
    fj -= floor(fj);
    c = r2tij(window, c, fi, fj, &p, &q);

    /*
    if(p > 1.0 || p < 0.0)
      hd_warn("x cell distance outside of bounds at %d(%d %d %d) : %f\n", c,
	      window->s2i[c],window->s2j[c],window->s2k[c],p);
    if(q > 1.0 || q < 0.0)
      hd_warn("y cell distance outside of bounds at %d(%d %d %d) : %f\n", c,
	      window->s2i[c],window->s2j[c],window->s2k[c],q);
    if(r > 1.0 || r < 0.0)
      hd_warn("r cell distance outside of bounds at %d(%d %d %d) : %e\n", ct,
	      window->s2i[c],window->s2j[c],window->s2k[c],r);
    */
    r = min(max(0.0, r), 1.0);

    /* Get the weights                                               */
    wincon->wgt[ct][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
    wincon->wgt[ct][1] = (1.0 - r) * (1.0 - p) * q;
    wincon->wgt[ct][2] = (1.0 - r) * p * (1.0 - q);
    wincon->wgt[ct][3] = (1.0 - r) * p * q;
    wincon->wgt[ct][4] = r * (1.0 - p) * (1.0 - q);
    wincon->wgt[ct][5] = r * (1.0 - p) * q;
    wincon->wgt[ct][6] = r * p * (1.0 - q);
    wincon->wgt[ct][7] = r * p * q;
    cl[cc] = c;

    /* Copy the streamline origin and offsets into origin, p,q and */
    /* r respectively, so that these quantities are available to   */
    /* dumping to sparse output files. Note that these variables   */
    /* may only be dumped to sparse file since interpolation from  */
    /* (x,y,z) format is not possible for the sparse coordinate of */
    /* the streamline origin.                                      */
    windat->origin[ct] = (double)cl[cc];
    windat->pc[ct] = p;
    windat->qc[ct] = q;
    windat->rc[ct] = r;
  }
  debug_c(window, D_TS, D_STRML);
}

/* END semi_lagrange_tu()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts sparse location and increments on the grid used for      */
/* semi-Lagrange (where grid corners lie on cell centers of the      */
/* regular grid) to a sparse location and increments on the regular  */
/* grid.                                                             */
/*-------------------------------------------------------------------*/
int t2rij(geometry_t *window, 
	  int c,         /* Sparse location in the transport grid    */
	  double p,      /* x increment in the transport grid        */
	  double q,      /* y increment in the transport grid        */
	  double *x,     /* x increment in the regular grid          */
	  double *y      /* y increment in the regular grid          */
	  )     
{
  int cr;
  double eps = 1e-10;

  cr = c;
  if (fabs(p) < 2.0*eps && fabs(q) < 2.0*eps) {   /* Centre          */
    *x = 0.5;
    *y = 0.5;
  } else if (p < 0.5 && q < 0.5) {                /* SW quadrant     */
    *x = p + 0.5;
    *y = q + 0.5;
  } else if (p < 0.5 && q >= 0.5) {               /* NW quadrant     */
    cr = window->yp1[c];
    *x = p + 0.5;
    *y = q - 0.5;
  } else if (p >= 0.5 && q >= 0.5) {              /* NE quadrant     */
    cr = window->xp1[window->yp1[c]];
    *x = p - 0.5;
    *y = q - 0.5;	
  } else if (p >= 0.5 && q < 0.5) {               /* SE quadrant     */
    cr = window->xp1[c];
    *x = p - 0.5;
    *y = q + 0.5;
  }
  /* cr is still a ghost cell, possible due to tri-diagonal error    */
  /* (e.g. on a non-orthogonal grid). Nugde the streamline into a    */
  /* wet cell.                                                       */
  if (window->wgst[cr]) {
    int cg = cr;
    hd_warn("Streamline from %d crossed coast.\n",c);
    cr = window->wgst[cr];
    if (window->xp1[cg] == cr || window->xm1[cg] == cr) *x = fabs(0.5 - p);
    if (window->yp1[cg] == cr || window->ym1[cg] == cr) *y = fabs(0.5 - q);
  }
  *x += window->s2i[cr];
  *y += window->s2j[cr];

  return(cr);
}

int t2rc(geometry_t *window, 
	 int c,          /* Sparse location in the transport grid    */
	 double p,       /* x increment in the transport grid        */
	 double q        /* y increment in the transport grid        */
	 )     
{
  int cr;
  double eps = 1e-10;

  cr = c;
  if (fabs(p) < 2.0*eps && fabs(q) < 2.0*eps) {   /* Centre          */
    return(cr);
  } else if (p < 0.5 && q < 0.5) {                /* SW quadrant     */
    return(cr);
  } else if (p < 0.5 && q >= 0.5) {               /* NW quadrant     */
    cr = window->yp1[c];
  } else if (p >= 0.5 && q >= 0.5) {              /* NE quadrant     */
    cr = window->xp1[window->yp1[c]];
  } else if (p >= 0.5 && q < 0.5) {               /* SE quadrant     */
    cr = window->xp1[c];
  }
  return(cr);
}

/* END t2rij()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts sparse location and increments on the regular grid to    */
/* the semi-Lagrange (where grid corners lie on cell centers of the  */
/* regular grid).                                                    */
/*-------------------------------------------------------------------*/
int r2tij(geometry_t *window, int c, double p, double q, double *x, double *y)
{
  int cr;
  double eps = 2.e-10;
  int pf = (p >= 0.5) ? 1 : 0;
  int qf = (q >= 0.5) ? 1 : 0;

  cr = c;
  if (fabs(p-0.5) < eps && fabs(q-0.5) < eps) {   /* Centre          */
    *x = 0.0;
    *y = 0.0;
  } else if (pf && qf) {                          /* NE quadrant     */
    *x = p - 0.5;
    *y = q - 0.5;
  } else if (pf && !qf) {                         /* SE quadrant     */
    cr = window->ym1[c];
    *x = p - 0.5;
    *y = q + 0.5;
  } else if (!pf && !qf) {                        /* SW quadrant     */
    cr = window->xm1[window->ym1[c]];
    *x = p + 0.5;
    *y = q + 0.5;	
  } else if (!pf && qf) {                         /* NW quadrant     */
    cr = window->xm1[c];
    *x = p + 0.5;
    *y = q - 0.5;
  }
  return(cr);
}

/* END r2tij()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the weights on the source grid given coordinates   */
/* (x,y) in the traget grid.                                         */
/* Use newval = int_val(ws, c, wgt, tr) to interpolate.              */
/*-------------------------------------------------------------------*/
int s2t_weights(geometry_t *ws, geometry_t *wt, double *dz, double z, 
		double *wgt, int ct)
{
  win_priv_t *wincon = wt->wincon;
  int c, cs, c2;
  double p, q, r;

  /* Get the sparse location in the source grid                      */
  cs = wt->mgm[ct] ;

  /* Get the horizontal increments                                   */
  c2 = wt->m2d[ct];
  p = wincon->d7[c2];
  q = wincon->d8[c2];

  /* Get the depth of the streamline in the target grid              */ 
  if (z) {
    /* Truncate depth to source cell depth                           */
    cs = ws->m2d[cs];
    z = max(z, ws->botz[cs]);
    c = ztoc_z(ws, cs, z, ws->cellz);
    r = (z - max(ws->cellz[c], ws->botz[cs])) / dz[c];
    /* If cs is a ghost cell, there is no sediment ghost and the     */
    /* potential exists for r<0, hence limit.                        */
    r = min(max(r, 0.0), 1.0);
  }
  else {
    c = cs;
    r = 0.0;
  }

  /* Get the weights                                                 */
  wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
  wgt[1] = (1.0 - r) * (1.0 - p) * q;
  wgt[2] = (1.0 - r) * p * (1.0 - q);
  wgt[3] = (1.0 - r) * p * q;
  wgt[4] = r * (1.0 - p) * (1.0 - q);
  wgt[5] = r * (1.0 - p) * q;
  wgt[6] = r * p * (1.0 - q);
  wgt[7] = r * p * q;
  return(c);
}

/* END s2t_weights()                                                 */
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
  int cc, c, cs;
  double *ntr = wincon->w5;
  int *cl = wincon->s2;
  int *c2cc = wincon->s4;

  /* Get the tracer concentrations at the forward timestep */
  memcpy(ntr, tr, window->sgsiz * sizeof(double));
  if (wincon->osl == 0) {
    if(wincon->togn & TOPRIGHT) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	c2cc[c] = cc;

	/*check_tridiagonal(window,tr,"before",c,cl[cc]);*/

	ntr[c] = int_val(window, cl[cc], wincon->wgt[c], tr);

	/*
	if(isnan(ntr[c])) {
	  print_tridiagonal(window,tr,c,cl[cc]);
	  s2ijk(window,c);
	  hd_quit("NaN at %f days\n",windat->days);
	}
	*/
      }
    } else {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	c2cc[c] = cc;
	ntr[c] = int_val_bl(window, cl[cc], wincon->wgt[c], tr);
      }
    }
  } else {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = cg = wincon->s1[cc];
      c2cc[c] = cc;
      ntr[c] = int_valo(window, wincon->wgt[c], wincon->lmap[cl[cc]], 
			tr, wincon->osl);
    }
    if (wincon->osl > 1) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = cg = wincon->s1[cc];
	set_monotonic(wincon->osl, c, &ntr[c], tr, wincon->lmap[cl[cc]]);
      }
    }
  }
  memcpy(tr, ntr, window->sgsiz * sizeof(double));
}

/* END semi_lagrange()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the weights from offsets read from file                  */
/*-------------------------------------------------------------------*/
void get_weights(geometry_t *window,     /* Window geometry          */
		 window_t *windat,       /* Window data              */
		 win_priv_t *wincon)     /* Window constants         */
{
  int c, cc;
  int *cl = wincon->s2;
  double p, q, r;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cl[cc] = (int)windat->origin[c];
    p = windat->pc[c];
    q = windat->qc[c];
    r = windat->rc[c];

    wincon->wgt[c][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
    wincon->wgt[c][1] = (1.0 - r) * (1.0 - p) * q;
    wincon->wgt[c][2] = (1.0 - r) * p * (1.0 - q);
    wincon->wgt[c][3] = (1.0 - r) * p * q;
    wincon->wgt[c][4] = r * (1.0 - p) * (1.0 - q);
    wincon->wgt[c][5] = r * (1.0 - p) * q;
    wincon->wgt[c][6] = r * p * (1.0 - q);
    wincon->wgt[c][7] = r * p * q;
  }
}

/* END get_weights()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates a value given the top right sparse grid              */
/*-------------------------------------------------------------------*/
double int_val(geometry_t *window, int c, double *wgt, double *in)
{
  double out;

  out = in[c] * wgt[0] +
    in[window->ym1[c]] * wgt[1] +
    in[window->xm1[c]] * wgt[2] +
    in[window->xm1[window->ym1[c]]] * wgt[3] +
    in[window->zp1[c]] * wgt[4] +
    in[window->ym1[window->zp1[c]]] * wgt[5] +
    in[window->xm1[window->zp1[c]]] * wgt[6] +
    in[window->ym1[window->xm1[window->zp1[c]]]] * wgt[7];

  return (out);
}

/* END int_val()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates a value given the bottom left sparse grid            */
/*-------------------------------------------------------------------*/
double int_val_bl(geometry_t *window, int c, double *wgt, double *in)
{
  double out;
  int **map = window->wincon->lmap;
  /*
  out = in[map[c][0]] * wgt[0] +
    in[map[c][1]] * wgt[1] +
    in[map[c][2]] * wgt[2] +
    in[map[c][3]] * wgt[3] +
    in[map[c][4]] * wgt[4] +
    in[map[c][5]] * wgt[5] +
    in[map[c][6]] * wgt[6] +
    in[map[c][7]] * wgt[7];
  */

  out = in[c] * wgt[0] +
    in[window->yp1[c]] * wgt[1] +
    in[window->xp1[c]] * wgt[2] +
    in[window->xp1[window->yp1[c]]] * wgt[3] +
    in[window->zp1[c]] * wgt[4] +
    in[window->yp1[window->zp1[c]]] * wgt[5] +
    in[window->xp1[window->zp1[c]]] * wgt[6] +
    in[window->yp1[window->xp1[window->zp1[c]]]] * wgt[7];

  return (out);
}


/* END int_val_bl()                                                     */
/*-------------------------------------------------------------------*/

#define TINY_INSIDE (1e-8)

/*-------------------------------------------------------------------*/
/* Locates the position of the streamline in sparse coordinates.     */
/* The origin for interpolations is in the top right horizontal and  */
/* bottom right vertical corner of the cell. Due to the orthogonal   */
/* grid (non-uniform grid spacing) the streamline distance must be   */
/* decremented for each individual cell.                             */
/*-------------------------------------------------------------------*/
int get_pos(geometry_t *window, int c, double nu, double nv, double nw,
            double *cx, double *cy, double *cz, double dt)
{
  win_priv_t *wincon = window->wincon;
  double *dzz = wincon->w9;
  int co, cm;
  double d;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;
  if (nw == SMALL)
    nw = 0.0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  d = nu * dt + (*cx);
  if (d < 0.0) {
    co = window->xp1[co];
    while (fabs(d) >= window->h1au1[window->m2d[co]]) {
      d += window->h1au1[window->m2d[co]];
      co = window->xp1[co];
    }
    d += window->h1au1[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= window->h1au1[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      cm = window->m2d[co];
      if (wincon->gmap[co] & L_EDGE && fabs(d) >= 0.5 * window->h1au1[cm]) {
	d = 0.5 * window->h1au1[cm] - TINY_INSIDE;
	break;
      }
      d -= window->h1au1[window->m2d[co]];
      co = window->xm1[co];
    }
  }
  *cx = fabs(d);

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = nv * dt + (*cy);
  if (d < 0.0) {
    co = window->yp1[co];
    while (fabs(d) >= window->h2au2[window->m2d[co]]) {
      d += window->h2au2[window->m2d[co]];
      co = window->yp1[co];
    }
    d += window->h2au2[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= window->h2au2[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      cm = window->m2d[co];
      if (wincon->gmap[co] & B_EDGE && fabs(d) >= 0.5 * window->h2au2[cm]) {
	d = 0.5 * window->h2au2[cm] - TINY_INSIDE;
	break;
      }
      d -= window->h2au2[window->m2d[co]];
      co = window->ym1[co];
    }
  }
  *cy = fabs(d);

  /*-----------------------------------------------------------------*/
  /* Adjust for converging curvilinear grids                         */
  while (*cx >= window->h1au1[window->m2d[co]] ||
	 *cy >= window->h2au2[window->m2d[co]]) {
    if (*cx >= window->h1au1[window->m2d[co]]) {
      cm = window->xm1[co];
      if (cm != window->xm1[cm]) {
	*cx -= window->h1au1[window->m2d[co]];
	co = cm;
      } else /* Streamline heading into land - keep it wet           */
	*cx -= 0.5 * window->h1au1[window->m2d[co]];
    }
    if (*cy >= window->h2au2[window->m2d[co]]) {
      cm = window->ym1[co];
      if (cm != window->ym1[cm]) {
	*cy -= window->h2au2[window->m2d[co]];
	co = window->ym1[co];
      } else
	*cy -= 0.5 * window->h2au2[window->m2d[co]];
    }
  }
  /* If the streamline moves horizontally into a wet cell above the  */
  /* free surface, then place it downwards into the first wet cell.  */
  if (!wincon->c1[co]) {
    while(!wincon->c1[co] && co != window->zm1[co])
      co = window->zm1[co];
  }

  /*-----------------------------------------------------------------*/
  /* z direction.                                                    */
  d = -nw * dt + (*cz);
  if (d < 0.0) {
    co = window->zm1[co];
    while (fabs(d) >= dzz[co]) {
      if (co == window->zm1[co]) {
	*cz = fabs(fmod(d, dzz[co]));
	d += dzz[co];
	return(co);
      }
      d += dzz[co];
      co = window->zm1[co];
    }
    d += dzz[co];
  }
  if (d > 0.0) {
    while (d >= dzz[co] && co != window->zp1[co]) {
      if (co == wincon->m2d[co] || !wincon->c1[window->zp1[co]]) {
	*cz = fmod(d, dzz[co]);
	return(co);
      }
      d -= dzz[co];
      co = window->zp1[co];
    }
  }
  *cz = fabs(d);
  /* Sanity check: *cz must be < dzz */
  if (dzz[co] && *cz > dzz[co]) *cz = fmod(*cz, dzz[co]);

  return (co);
}


/* Note: diagnostics are written to the logfile */
int get_pos_diag(geometry_t *window, int c, double nu, double nv, double nw,
		 double *cx, double *cy, double *cz, double dt)
{
  win_priv_t *wincon = window->wincon;
  double *dzz = wincon->w9;
  int co, cm, n;
  int maxiter = 100;
  double d, dd;
  double t = 0.0;
  int cd = 0;

  /* Initialise */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;
  if (nw == SMALL)
    nw = 0.0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  if (cg == cd && window->windat->days >= t) hd_warn("Starting cell %d @ %d(%d %d %d)\n",c, co, window->s2i[co],window->s2j[co],window->s2k[co]);
  d = dd = nu * dt + (*cx);
  if (d < 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting x direction, d < 0 : %f\n", d);
    n = 0;
    co = window->xp1[co];
    while (fabs(d) >= window->h1au1[window->m2d[co]]) {
      d += window->h1au1[window->m2d[co]];
      co = window->xp1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +x direction d=%f, c=%d(%d,%d,%d)\n",
		window->windat->days, dd, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
    d += window->h1au1[window->m2d[co]];
  }
  if (d > 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting x direction, d > 0 : %f\n", d);
    n = 0;
    while (fabs(d) >= window->h1au1[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      cm = window->m2d[co];
      if (wincon->gmap[co] & L_EDGE && fabs(d) >= 0.5 * window->h1au1[cm]) {
	d = 0.5 * window->h1au1[cm] - TINY_INSIDE;
	break;
      }
      d -= window->h1au1[window->m2d[co]];
      co = window->xm1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -x direction d=%f, c=%d(%d,%d,%d)\n",
		window->windat->days, dd, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
  }
  *cx = fabs(d);
  if (cg == cd && window->windat->days >= t) {
    hd_warn("x direction done, co = %d(%d %d %d), cx = %f\n", co, window->s2i[co],window->s2j[co],window->s2k[co], *cx);
    s2ijk(window, co);
  }

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = dd = nv * dt + (*cy);
  if (d < 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting y direction, d < 0 : %f\n", d);
    n = 0;
    co = window->yp1[co];
    while (fabs(d) >= window->h2au2[window->m2d[co]]) {
      d += window->h2au2[window->m2d[co]];
      co = window->yp1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +y direction d=%f, c=%d(%d,%d,%d)\n",
		window->windat->days, dd, c, window->s2i[c], window->s2j[c], window->s2k[c]);

    }
    if (cg == cd && window->windat->days >= t) hd_warn("End y direction : %d %f %f\n", co, d, window->h2au2[window->m2d[co]]);
    d += window->h2au2[window->m2d[co]];
  }
  if (d > 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting y direction, d > 0 : %f\n", d);
    n = 0;
    while (fabs(d) >= window->h2au2[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      cm = window->m2d[co];
      if (wincon->gmap[co] & B_EDGE && fabs(d) >= 0.5 * window->h2au2[cm]) {
	d = 0.5 * window->h2au2[cm] - TINY_INSIDE;
	break;
      }
      d -= window->h2au2[window->m2d[co]];
      co = window->ym1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -y direction d=%f, c=%d(%d,%d,%d)\n",
		window->windat->days, dd, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
  }
  *cy = fabs(d);
  if (cg == cd && window->windat->days >= t) {
    hd_warn("y direction done, co = %d(%d %d %d), cy = %f\n", co, window->s2i[co],window->s2j[co],window->s2k[co], *cy);
    s2ijk(window, co);
  }

  /*-----------------------------------------------------------------*/
  /* Adjust for converging curvilinear grids                         */
  if (cg == cd && window->windat->days >= t) hd_warn("Starting curvilinear adjustment\n");
  n = 0;
  while (*cx >= window->h1au1[window->m2d[co]] ||
	 *cy >= window->h2au2[window->m2d[co]]) {
    if (*cx >= window->h1au1[window->m2d[co]]) {
      cm = window->xm1[co];
      if (cm != window->xm1[cm]) {
	*cx -= window->h1au1[window->m2d[co]];
	co = cm;
      } else /* Streamline heading into land - keep it wet           */
	*cx -= 0.5 * window->h1au1[window->m2d[co]];
    }
    if (*cy >= window->h2au2[window->m2d[co]]) {
      cm = window->ym1[co];
      if (cm != window->ym1[cm]) {
	*cy -= window->h2au2[window->m2d[co]];
	co = window->ym1[co];
      } else
	*cy -= 0.5 * window->h2au2[window->m2d[co]];
    }
    n++;
    if (n >= maxiter) 
      hd_quit("Maximum iterations exceeded at %f days: converging grids, c=%d(%d,%d,%d)\n",
	      window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);
  }
  if (cg == cd && window->windat->days >= t) 
    hd_warn("Curvilinear adjustment OK, co=%d(%d %d %d) cx=%f cy=%f\n",co,window->s2i[co],window->s2j[co],window->s2k[co], *cx, *cy);
  /* If the streamline moves horizontally into a wet cell above the  */
  /* free surface, then place it downwards into the first wet cell.  */
  if (!wincon->c1[co]) {
    while(!wincon->c1[co] && co != window->zm1[co])
      co = window->zm1[co];
  }

  /*-----------------------------------------------------------------*/
  /* z direction.                                                    */
  d = dd = -nw * dt + (*cz);
  if (d < 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting z direction, d < 0 : %f %f %e\n", d, dzz[co], nw);
    n = 0;
    co = window->zm1[co];
    while (fabs(d) >= dzz[co]) {
      if (co == window->zm1[co]) {
	*cz = fabs(fmod(d, dzz[co]));
	d += dzz[co];
	return(co);
      }
      d += dzz[co];
      co = window->zm1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -z direction d=%f, c=%d(%d,%d,%d), co=%d(%d,%d,%d), w = %f, dz = %f\n",
		window->windat->days, dd,
		c, window->s2i[c], window->s2j[c], window->s2k[c], 
		co, window->s2i[co], window->s2j[co], window->s2k[co], 
		nw, dzz[co]);
    }
    d += dzz[co];
  }

  if (d > 0.0) {
    if (cg == cd && window->windat->days >= t) hd_warn("Starting z direction, d > 0 : %f %f %e\n",d, dzz[co], nw);
    n = 0;
    while (d >= dzz[co] && co != window->zp1[co]) {
      if (co == wincon->m2d[co] || !wincon->c1[window->zp1[co]]) {
	*cz = fmod(d, dzz[co]);
	return(co);
      }
      d -= dzz[co];
      co = window->zp1[co];
      n++;
      if (n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +z direction d=%f, c=%d(%d,%d,%d), co=%d(%d,%d,%d), w = %f, dz = %f\n",
		window->windat->days, dd,
		c, window->s2i[c], window->s2j[c], window->s2k[c], 
		co, window->s2i[co], window->s2j[co], window->s2k[co], 
		nw, dzz[co]);
    }
  }
  *cz = fabs(d);

  /* Sanity check: *cz must be < dzz */
  if (*cz > dzz[co]) *cz = fmod(*cz, dzz[co]);

  if (cg == cd && window->windat->days >= t) {
    hd_warn("z direction done, co = %d(%d %d %d), cz = %f\n", co, window->s2i[co],window->s2j[co],window->s2k[co], *cz);
    s2ijk(window, co);
  }
  return (co);
}

/* END get_pos()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Locates the position of the streamline in sparse coordinates.     */
/* The origin for interpolations is in the bottom left horizontal    */
/* and bottom left vertical corner of the cell. Due to the orthogonal*/
/* grid (non-uniform grid spacing) the streamline distance must be   */
/* decremented for each individual cell.                             */
/*-------------------------------------------------------------------*/
int get_pos_bl(geometry_t *window, int c, double nu, double nv, double nw,
	       double *cx, double *cy, double *cz, double dt)
{
  win_priv_t *wincon = window->wincon;
  double *dzz = wincon->w9;
  int co, cm, cp, n;
  double d;
  int maxiter = 100;
  int qf = 0;

  /* Initialise */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;
  if (nw == SMALL)
    nw = 0.0;

  /* x direction.  */
  d = nu * dt - (*cx);
  if (d > 0.0) {
    cp = window->m2d[co];
    co = window->xm1[co];
    n = 0;
    while (fabs(d) >= window->h1au1[cp]) {
      d -= window->h1au1[cp];
      cp = window->m2d[co];
      co = window->xm1[co];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -x direction c=%d(%d,%d,%d)\n",
		window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);

    }
    d -= window->h1au1[cp];
  }
  if (d < 0.0) {
    cp = window->m2d[window->xp1[co]];
    n = 0;
    while (fabs(d) >= window->h1au1[cp]) {
      d += window->h1au1[cp];
      co = window->xp1[co];
      cp = window->m2d[window->xp1[co]];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +x direction c=%d(%d,%d,%d)\n",
		window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
  }
  *cx = fabs(d);

  /* y direction.  */
  d = nv * dt - (*cy);
  if (d > 0.0) {
    cp = window->m2d[co];
    co = window->ym1[co];
    n = 0;
    while (fabs(d) >= window->h2au2[cp]) {
      d -= window->h2au2[cp];
      cp = window->m2d[co];
      co = window->ym1[co];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -y direction, c=%d(%d,%d,%d)\n",
		window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
    d -= window->h2au2[cp];
  }
  if (d < 0.0) {
    cp = window->m2d[window->yp1[co]];
    n = 0;
    while (fabs(d) >= window->h2au2[cp]) {
      d += window->h2au2[cp];
      co = window->yp1[co];
      cp = window->m2d[window->yp1[co]];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +y direction, c=%d(%d,%d,%d)\n",
		window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
  }
  *cy = fabs(d);

  /* Adjust for converging curvilinear grids */
  n = 0;
  while (*cx >= window->h1au1[window->xp1[window->m2d[co]]] ||
	 *cy >= window->h2au2[window->yp1[window->m2d[co]]]) {
    if (*cx >= window->h1au1[window->xp1[window->m2d[co]]]) {
      cm = window->xp1[co];
      if (cm != window->xp1[cm]) {
	*cx -= window->h1au1[window->xp1[window->m2d[co]]];
	co = cm;
      } else /* Streamline heading into land - keep it wet */
	*cx -= 0.5 * window->h1au1[window->xp1[window->m2d[co]]];
    }
    if (*cy >= window->h2au2[window->yp1[window->m2d[co]]]) {
      cm = window->yp1[co];
      if (cm != window->yp1[cm]) {
	*cy -= window->h2au2[window->yp1[window->m2d[co]]];
	co = window->yp1[co];
      } else
	*cy -= 0.5 * window->h2au2[window->yp1[window->m2d[co]]];
    }
    n++;
    if (qf && n >= maxiter) 
      hd_quit("Maximum iterations exceeded at %f days: converging grids, c=%d(%d,%d,%d)\n",
	      window->windat->days, c, window->s2i[c], window->s2j[c], window->s2k[c]);
  }
  /* If the streamline moves horizontally into a wet cell above the  */
  /* free surface, then place it downwards into the first wet cell.  */
  if (!wincon->c1[co]) {
    while(!wincon->c1[co] && co != window->zm1[co])
      co = window->zm1[co];
  }

  /* z direction.  */
  d = -nw * dt + (*cz);
  if (d < 0.0) {
    co = window->zm1[co];
    n = 0;
    while (fabs(d) >= dzz[co]) {
      if (co == window->zm1[co]) {
	*cz = fabs(fmod(d, dzz[co]));
	d += dzz[co];
	return(co);
      }
      d += dzz[co];
      co = window->zm1[co];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: -z direction c=%d(%d,%d,%d), co=%d(%d,%d,%d), w = %f, dz = %f\n",
		window->windat->days,
		c, window->s2i[c], window->s2j[c], window->s2k[c], 
		co, window->s2i[co], window->s2j[co], window->s2k[co], 
		nw, dzz[co]);
    }
    d += dzz[co];
  }
  if (d > 0.0) {
    n = 0;
    while (d >= dzz[co] && co != window->zp1[co]) {
      if (co == wincon->m2d[co] || !wincon->c1[window->zp1[co]]) {
	*cz = fmod(d, dzz[co]);
	return(co);
      }
      d -= dzz[co];
      co = window->zp1[co];
      n++;
      if (qf && n >= maxiter) 
	hd_quit("Maximum iterations exceeded at %f days: +z direction c=%d(%d,%d,%d), co=%d(%d,%d,%d), w = %f, dz = %f\n",
		window->windat->days,
		c, window->s2i[c], window->s2j[c], window->s2k[c], 
		co, window->s2i[co], window->s2j[co], window->s2k[co], 
		nw, dzz[co]);
    }
  }
  *cz = fabs(d);
  /* Sanity check: *cz must be < dzz */
  if (*cz > dzz[co]) *cz = fmod(*cz, dzz[co]);
  return (co);
}

/* END get_pos_bl()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to centre a velocity on the u1 or u2 point                */
/*-------------------------------------------------------------------*/
void vel_center(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon, /* Window geometry / constants */
                double *nu,     /* Centered u1 velocity */
                double *nv,     /* Centered u2 velocity */
                double *nw      /* Centered z velocity */
  )
{
  int cc, c;
  int p1;
  double wtop;
  double *u1, *u2;
  int n;

  if (wincon->means & TRANSPORT) {
    u1 = windat->u1m;
    u2 = windat->u2m;
  } else {
    u1 = windat->u1;
    u2 = windat->u2;
  }

  /* Calculate the u1 and u2 values at the cell centre */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    p1 = window->xp1[c];
    nu[c] = 0.5 * (u1[c] + u1[p1]);
    if (fabs(nu[c]) < SMALL)
      nu[c] = SMALL;
    p1 = window->yp1[c];
    nv[c] = 0.5 * (u2[c] + u2[p1]);
    if (fabs(nv[c]) < SMALL)
      nv[c] = SMALL;
  }

  /* Calculate the w values at the cell centre */
  if (wincon->means & TRANSPORT) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      p1 = window->zp1[c];
      nw[c] = 0.5 * (windat->wm[c] + windat->wm[p1]);
      if (fabs(nw[c]) < SMALL)
        nw[c] = SMALL;
    }
  } else {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      p1 = window->zp1[c];
      wtop =
        (cc <= wincon->vcs) ? windat->wtop[window->m2d[c]] : windat->w[p1];
      nw[c] = 0.5 * (windat->w[c] + wtop);
      if (fabs(nw[c]) < SMALL)
        nw[c] = SMALL;
    }
  }

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	nu[c] = (nu[c] < 0.0) ? 0.0 : nu[c];
      }
    }
    if (open->ocodex & L_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	nu[c] = (nu[c] > 0.0) ? 0.0 : nu[c];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	nv[c] = (nv[c] < 0.0) ? 0.0 : nv[c];
      }
    }
    if (open->ocodey & B_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	nv[c] = (nv[c] > 0.0) ? 0.0 : nv[c];
      }
    }
  }

  /* Check for undefined velocities */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      if (fabs(nu[c]) > wincon->velmax || isnan(nu[c])) {
	hd_warn("Undefined u1 velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nu[c]);
	nu[c] = SMALL;
      }
      if (fabs(nv[c]) > wincon->velmax || isnan(nv[c])) {
	hd_warn("Undefined u2 velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nv[c]);
	nv[c] = SMALL;
      }
      if (fabs(nw[c]) > wincon->velmax || isnan(nw[c])) {
	hd_warn("Undefined w velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nw[c]);
	nw[c] = SMALL;
      }
    }
}


void vel_center_w(geometry_t *window, /* Processing window */
		  window_t *windat, /* Window data structure */
		  win_priv_t *wincon, /* Window geometry / constants */
		  double *nw      /* Centered z velocity */
		  )
{
  int cc, c;
  int p1;
  double wtop;

  /* Calculate the w values at the cell centre */
  if (wincon->means & TRANSPORT) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      p1 = window->zp1[c];
      nw[c] = 0.5 * (windat->wm[c] + windat->wm[p1]);
      if (fabs(nw[c]) < SMALL)
        nw[c] = SMALL;
    }
  } else {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      p1 = window->zp1[c];
      wtop =
        (cc <= wincon->vcs) ? windat->wtop[window->m2d[c]] : windat->w[p1];
      nw[c] = 0.5 * (windat->w[c] + wtop);
      if (fabs(nw[c]) < SMALL)
        nw[c] = SMALL;
    }
  }
}

/* END vel_center()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the vertical cell spacing centred on the grid layers          */
/*-------------------------------------------------------------------*/
void set_dzz(geometry_t *window, double *dzz)
{
  win_priv_t *wincon = window->wincon;
  int c, cc, c2, zm1, nn;

  /* Initialize in case a streamline lands in a dry cell, a valid    */
  /* (i.e. non inf) increment is computed.                           */
  memcpy(dzz, wincon->dz, window->sgsiz * sizeof(double));
  /* Overwrite with valid wet cell thicknesses                       */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zm1 = window->zm1[c];
    dzz[zm1] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
  }
  /* Set the surface layer                                           */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = c2 = wincon->s1[cc];
    dzz[c] = 0.5 * wincon->dz[c];
    /* Set the surface and bottom boundary conditions                */
    while (c2 != window->zp1[c2]) {
      c2 = window->zp1[c2];
      dzz[c2] = dzz[c];
    }
    c2 = wincon->i1[cc];
    dzz[window->zm1[c2]] = 0.5 * wincon->dz[c2];
  }
  /* FRONT and RIGHT open boundaries                                 */
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	c2 = open->oi1_e1[cc];
	dzz[c] = dzz[c2];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	c2 = open->oi1_e2[cc];
	dzz[c] = dzz[c2];
      }
    }
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
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    c2 = window->zp1[c];
    dzz[c] = dzz[c2];
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
  int nwf = 0;                  /* Newly wetted cell flag */

  /*-----------------------------------------------------------------*/
  /* Merge newly wetted cells with the cell below if required.       */
  if (wincon->fillf & (LOCAL|LOCALER) && wincon->tmode & SET_AIJ) 
    nwf = 1;

  /*-----------------------------------------------------------------*/
  /* Set dz at the surface to account for the updated elevation. */
  for(cc = 1; cc <= window->b2_t; cc++) {
    c = max(window->sur_t[cc], window->nsur_t[cc]);
    cs = window->m2d[c];
    cn = (nwf) ? max(window->sur_t[cc], window->nsur_t[cc]) : window->nsur_t[cc];
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
      c = (nwf) ? max(window->sur_t[cc], window->nsur_t[cc]) : window->nsur_t[cc];
      cs = window->m2d[c];
      top = windat->eta[cs];
      bot = DRY_FRAC * wincon->hmin + window->botz[cs];
      if (top > bot && c != window->zm1[c]) {
        wincon->s2[vc] = c;
	wincon->m2d[c] = c;
	if (wincon->trasc & LAGRANGE)
	  wincon->i1[vc] = window->bot_t[cc];
        wincon->i2[vc] = window->nsur_t[cc];
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
	double *u2 = windat->u2;
	double *dzu1 = windat->dzu1;
	double *dzu2 = windat->dzu2;
	/* Get the local volume budget                               */
	for(cc = 1; cc <= nvec; cc++) {
	  int xp1, yp1, xp1s, yp1s;
	  double div;
	  c = vec[cc];
	  cs = window->m2d[c];

	  xp1 = window->xp1[c];
	  yp1 = window->yp1[c];
	  xp1s = window->m2d[xp1];
	  yp1s = window->m2d[yp1];

	  div = u1[xp1] * dzu1[xp1] * window->h2au1[xp1s] * wincon->mdx[xp1s] -
	    u1[c] * dzu1[c] * window->h2au1[cs] * wincon->mdx[cs] +
	    u2[yp1] * dzu2[yp1] * window->h1au2[yp1s] * wincon->mdy[yp1s] -
	    u2[c] * dzu2[c] * window->h1au2[cs] * wincon->mdy[cs];
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

#define LF_VOL   1
#define LF_AIJ   2
#define LF_FOF   4
#define LF_WETM  8
#define LF_AIJM  16
#define LF_OBCM  32
#define LF_OBCF  64
#define LF_DESF  128
#define LF_DESM  256
#define LF_DRY   -1.0

/*-------------------------------------------------------------------*/
/* Routine to implement local mass conservation                      */
/*-------------------------------------------------------------------*/
void local_fill(geometry_t *window,  /* Window geometry             */
		window_t *windat,    /* Window data                 */
		win_priv_t *wincon,  /* Window constants            */
		int mode
		)
{
  double (*aij) (geometry_t *, window_t *, win_priv_t *, int);
  int c, cs, ci, c1, c2, cc, co;
  int i, j, k, n = 0, ce;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double *Vd = wincon->w7;
  double *Vi = wincon->w8;
  double **A = wincon->wgt;
  double *ovol = wincon->p1;
  double aerr, verr, tverr, vrelerr;
  double cellvol, d1;
  int flag;
  double th = 0.001;
  int maxiter = 20;
  int forcev = 0;

  aij = aij_1;
  if (wincon->fillf & MONOTONIC) forcev = 1;

  /*-----------------------------------------------------------------*/
  /* Get the volume at the previous timestep                         */
  if (mode) {
    for(cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      ovol[c] += windat->waterss[c] * windat->dttr;	
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Do the iterations                                               */  
  th /= 100.0;
  do {
    memset(Vi, 0, window->sgsiz * sizeof(double));
    memset(Vd, 0, window->sgsiz * sizeof(double));
    if (wincon->fillf & DIAGNOSE) {
      memset(windat->vcorr, 0, window->sgsiz * sizeof(double));
      memset(windat->acorr, 0, window->sgsiz * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Eqn2: add sumi(aijVi(t)) to Vi                                */
    /* Vi is destination volume for each source cell.                */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      cs = S[cc];
      cellvol = window->cellarea[c2] * wincon->dz[c];
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	/* Only compute Vi for valid source cells, including OB      */
	/* source cells.                                             */
	if (mask[c1] & LF_WETM ) {
	  if (mask[c1] & LF_OBCM) {
	    /* Do not include OB destination cells in Vi             */
	    Vi[c1] += (A[c][j] * cellvol);
	    /* Vd is the total volume of valid interior destination  */
	    /* cells.                                                */
	    Vd[c1] += cellvol;
	  }
	  /* mask is non-zero if any destination cell is an obc       */
	  if (!(mask[c] & LF_OBCM)) mask[c1] |= LF_DESF;
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the fractional volume correction                          */
    /* set mask = LF_AIJM if cell meets criteria for vol corr        */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      if (mask[c] & LF_AIJM) mask[c] &= ~LF_AIJM;
      if (!(mask[c] & LF_DESF) && mask[c] & LF_OBCF && Vi[c]) 
	mask[c] |= LF_AIJM;
    }

    /*---------------------------------------------------------------*/    
    /* Diagnostics                                                   */
    /* tverr = total volume error                                    */
    /* vrelerr is relative cell volume error                         */
    /* verr is max absolute relative cell vol error at cell ce       */
    /* vcorr contains % cell volume errors                           */
    /* acorr contains cell volume errors                             */
    ce = 1;
    tverr = verr = 0.0;
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      if (mask[c] & LF_AIJM) {
	tverr += (ovol[c] - Vi[c]);
	vrelerr = ovol[c] / Vi[c] - 1.0;
	if (fabs(vrelerr) > verr) {
	  verr = fabs(vrelerr);
	  ce = c;
	}
	if (windat->vcorr)
	  windat->vcorr[c] = 100.0 * vrelerr;
	if (windat->acorr)
	  windat->acorr[c] = (ovol[c] - Vi[c]);
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the fractional coefficient correction and the new Aij     */
    flag = (forcev && n == maxiter - 1) ? LF_FOF : 0;
    aerr = aij(window, windat, wincon, flag);

    /*
    printf("%f %d %f %f %f %d(%d %d %d)\n",windat->days,n,aerr,verr,tverr,ce,window->s2i[ce],window->s2j[ce],window->s2k[ce]);
    */
    n++;
  } while (verr > th && n < maxiter);

  /* Set the final Vi for use later */
  memset(Vi, 0, window->sgsiz * sizeof(double));
  if (wincon->fillf & DIAGNOSE)
    memset(windat->vcorr, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    cellvol = window->cellarea[c2] * wincon->dz[c];
    d1 = 0.0;
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      d1 += A[c][j];
      if (mask[c1] & LF_WETM) {
	if (mask[c1] & LF_OBCM)
	  Vi[c1] += (A[c][j] * cellvol);
	if (!(mask[c] & LF_OBCM))
	  mask[c1] |= LF_DESF;
      }
      /*
      if (A[c][j] < 0.0 || A[c][j] > 1.0) hd_quit("Aij out of bounds : %d %d %f\n",
        c, j, A[c][j]);
      */
    }
    /*
    if (fabs(d1-1.0) > SMALL) {
      print_tridiagonal(window,windat->tr_wc[6],c,S[cc]);
      hd_quit("Aij does not sum to 1 : %d %d %f\n", c, cc, d1);
    }
    */
  }
  /* set mask = LF_AIJM if cell meets criteria for vol corr             */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    if (mask[c] & LF_AIJM) mask[c] &= ~LF_AIJM;
    if (!(mask[c] & LF_DESF) && mask[c] & LF_OBCF && Vi[c]) mask[c] |= LF_AIJM;
  }
  /* Set diagnostic arrays for use later */
  if (wincon->fillf & DIAGNOSE) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      if (mask[c] & LF_AIJM) {
	windat->vcorr[c] = 100.0 * (ovol[c] / Vi[c] - 1.0);	\
	windat->acorr[c] = (ovol[c] - Vi[c]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Print diagnostics to file if required                           */
  if (wincon->fillf & (DIAGNOSE|DIAGNOSE_BGC)) {
    if (windat->t >= tstrans.tsout - DT_EPS) {
      double ver = 0.0;      /* Computed volume error                */
      double v_in = 0.0;     /* Total inflow volume                  */
      double v_out = 0.0;    /* Total outflow volume                 */
      double v0 = 0.0;       /* Total volume at time t-1             */
      double v1 = 0.0;       /* Total volume at time t               */
      double vr = 0.0, vt = 0.0;

      /* Volume diagnostics for surface interior cells               */
      for (cc = 1; cc <= wincon->vcs1; cc++) {
	c = wincon->s1[cc];
	if (mask[c] & LF_AIJM)
	  ver += (Vi[c] - ovol[c]);
	else
	  v_out -= (Vi[c] - ovol[c]);
	v0 += ovol[c];
	v1 += window->cellarea[window->m2d[c]] * wincon->dz[c];
      }
      /* Sub-surface interior cells                                  */
      for (cc = wincon->vcs + 1; cc <= wincon->vc1; cc++) {
	c = wincon->s1[cc];
	if (mask[c] & LF_AIJM)
	  ver += (Vi[c] - ovol[c]);
	else
	  v_out -= (Vi[c] - ovol[c]);
	v0 += ovol[c];
	v1 += window->cellarea[window->m2d[c]] * wincon->dz[c];
      }
      /* Boundary cells                                              */
      for (cc = wincon->vcs1+1; cc <= wincon->vcs; cc++) {
	c = wincon->s1[cc];
	v_in += Vi[c];
	/*
	if(window->s2i[c]==1&&window->s2j[c]==13)vr += Vi[c];
	if(window->s2i[c]==48)vt += Vi[c];
	*/
      }
      for (cc = wincon->vc1+1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	v_in += Vi[c];
	/*
	if(window->s2i[c]==1&&window->s2j[c]==13)vr += Vi[c];
	if(window->s2i[c]==48)vt += Vi[c];
	*/
      }

      fprintf(tstrans.fp, "%f %f %f %f %f %f %f %f %f\n", windat->days, 
	      100.0 * verr, 
	      100.0 * aerr, (double)n,
	      v1 - v0 - v_in + v_out, ver,
	      v1 - v0, v0, v_out);
      fflush(tstrans.fp);
      tstrans.tsout += tstrans.tsdt;
    }
  }

  /* Compute local fill diagnostics */
  /* 1 to wincon->vcs1 = surface layer wet cells */
  /* wincon->vcs1+1 to wincon->vcs = surface layer OBC cells */
  /* wincon->vcs+1 to wincon->vc1 = sub-surface wet cells */
  /* wincon->vc1+1 to wincon->vc = sub-surface OBC cells */
}

/* Proportional */
double aij_1(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag) {
  int c, cc, c1, c2, cs, j;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double **A = wincon->wgt;
  double *ovol = wincon->p1;
  double *Vi = wincon->w8;
  double obcf, Acorr, Atot;
  double aerr = 0.0;

  /* Correct aij for volume errors */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    Atot = 0.0;

    if (!(mask[c] & LF_OBCM)) continue;
    if (!(mask[c] & LF_DESM)) continue;

    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (mask[c1] & LF_AIJM)
	A[c][j] *= ovol[c1] / Vi[c1];
      if(mask[c1] & LF_WETM)
	Atot += A[c][j];
    }

    Acorr = Atot - 1.0;
    aerr = max(aerr, fabs(Acorr));
    /*  if LF_FOF, stop at volume correction                          */
    if (flag & LF_FOF) continue;

    if (Atot > 1.0 || mask[c] & LF_OBCF) {
      /* distribute sum A error proportionally                        */
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM)
	  A[c][j] = A[c][j] / Atot;
      }
    } else {
      /* assign (+ve) sum A error to open bdy cell                      */
      /* NB Assigns all error to first OB source cell encountered.      */
      /* The order could matter if cell c has inputs from 2 boundaries. */
      /* We could put sum A error on any source cell with mask[]!=LF_AIJM*/
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM && !(mask[c1] & LF_OBCM)) {
	  A[c][j] -= Acorr;
	  continue;
	}
      }
    }
  }
  return(aerr);
}

/* Signed inverse */
double aij_2(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag) {
  int c, cc, c1, c2, cs, j;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double **A = wincon->wgt;
  double *Vd = wincon->w7;
  double *Vi = wincon->w8;
  double *ovol = wincon->p1;
  double Acorr, Atot, Nsource, obcf;
  double aerr = 0.0;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    
    if (!(mask[c] & LF_OBCM)) continue;
    if (!(mask[c] & LF_DESM)) continue;

    Acorr = Atot = Nsource = 0.0;
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (mask[c1] & LF_AIJM) {
	if (Vi[c1] > ovol[c1])
	  A[c][j] *= ovol[c1] / Vi[c1];
        else if (Vd[c1] > 0.0)
	  A[c][j] += (ovol[c1] - Vi[c1]) / Vd[c1];
        if (mask[c1] & LF_WETM) {
	  Atot += A[c][j];
	  Nsource += 1;
	}
      }
    }
    Acorr = Atot - 1.0;
    aerr = max(aerr, fabs(Acorr));
    
    /* if flag & LF_FOF dont do sum A correction */
    if (flag & LF_FOF) continue;

    if (Atot > 1.0) {
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM)
	  A[c][j] /= Atot;
      }
    } else if (Atot < 1.0 && mask[c] & LF_OBCF && Nsource > 0.0) {
      Acorr = Acorr / Nsource;
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM)
	  A[c][j] -= Acorr;
      }
    } else {
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (!(mask[c1] & LF_OBCM)) {           
	  A[c][j] += (1.0 - Atot);
	  continue;
	}
      }
    }
  }
  return(aerr);
}

/* equal split                                                      */
double  aij_3(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag) {
  int c, cc, c1, c2, cs, j;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double **A = wincon->wgt;
  double *ovol = wincon->p1;
  double *Vd = wincon->w7;
  double *Vi = wincon->w8;
  double Acorr, Atot, Nsource, obcf;
  double aerr = 0.0;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];

    if (!(mask[c] & LF_OBCM)) continue;
    if (!(mask[c] & LF_DESM)) continue;

    Acorr = Atot = Nsource = 0.0;
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      /* Spread volume correction over destination cells              */      
      if (mask[c1] & LF_AIJM && Vd[c1] > 0.0)
	A[c][j] += (ovol[c1] - Vi[c1]) / Vd[c1];
      if (mask[c1] & LF_WETM) {
	Atot += A[c][j];
	Nsource += 1;
      }
    }
    /* Error in sum aij                                              */
    Acorr = Atot - 1.0;
    aerr = max(aerr, fabs(Acorr));
    if (flag & LF_FOF || Nsource == 0) continue;
      
    /* if obcf, or Acorr > 0, spread Acorr over source cells         */
    if (mask[c] & LF_OBCF || Acorr > 0) {
      Acorr = Acorr / Nsource;
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM)
	  A[c][j] -= Acorr;
      }
    } else {
      /* put Acorr on OB cell */
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM && !(mask[c1] & LF_OBCM)) {
	  A[c][j] -= Acorr;
	  continue;
	}
      }
    }
  }
  return(aerr);
}

/* equal split in vol, sum aij weighted inversely with vol error   */
double  aij_4(geometry_t *window, window_t *windat, win_priv_t *wincon, int flag) {
  int c, cc, c1, c2, cs, j;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double **A = wincon->wgt;
  double *ovol = wincon->p1;
  double *Vd = wincon->w7;
  double *Vi = wincon->w8;
  double Acorr, Atot, Nsource, obcf, Wtot, TOL;
  double aerr = 0.0;
  double Wt[8];

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    Acorr = 0.0;

    if (!(mask[c] & LF_OBCM)) continue;
    if (!(mask[c] & LF_DESM)) continue;

    Acorr = Atot = Nsource = Wtot = 0.0;
    /* TOL sets limit to inverse weighting. Should be target rel error. */
    TOL = 0.001;
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      /*   Spread volume correction over destination cells              */	
      if (mask[c1] & LF_AIJM && Vd[c1] > 0.0)
	A[c][j] += (ovol[c1] - Vi[c1]) / Vd[c1];
      
      if (mask[c1] & LF_WETM) {
	Atot += A[c][j];
	Nsource += 1;

	/* Compute weights for sum A error distribution                   */
	Wt[j] = 1.0 / TOL;
	if (mask[c1] & LF_AIJM && Vi[c1] > 0.0)
	  Wt[j] = 1.0 / (TOL + fabs(ovol[c1] / Vi[c1] - 1.0));
	Wtot += Wt[j];
      }
    }

    /* Error in sum aij */
    Acorr = Atot - 1.0;
    aerr = max(aerr, fabs(Acorr));
    if (flag & LF_FOF || Nsource == 0) continue;
      
    /* if obcf or Acorr>0, spread Acorr over source cells with weights Wt */
    if (mask[c] & LF_OBCF || Acorr > 0) {
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM)
	  A[c][j] -= Acorr*Wt[j]/Wtot;
      }
    } else {
      /* put Acorr on OB cell */
      for (j = 0; j < wincon->nosl; j++) {
	c1 = wincon->lmap[cs][j];
	if (mask[c1] & LF_WETM && !(mask[c1] & LF_OBCM)) {
	  A[c][j] -= Acorr;
	  continue;
	}
      }
    }
  }
 return(aerr);
}


/* END local_fill()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reset the Aij weights at ghost cells                              */
/*-------------------------------------------------------------------*/
void reset_Aij(geometry_t *window,  /* Window geometry     */
	       window_t *windat,    /* Window data         */
	       win_priv_t *wincon   /* Window constants    */
	       )
{
  int c, cs, cc, j, k, c1, c2, ci, zp1, obcf;
  int *S = wincon->s2;
  int *mask = wincon->s3;
  double **A= wincon->wgt;
  double *msk = wincon->w6;
  double *ovol = wincon->p1;
  int osl = (wincon->osl) ? wincon->osl : 1;
  int nlay = wincon->nosl / (osl + 1);
  double Atot, Anot, Ack;
  int check = 0;

  /*-----------------------------------------------------------------*/
  /* Set the mask at non-wet cells. These are cells whose weights    */
  /* will be transferred to a wet cell. The mask contains the sparse */
  /* location of the wet cell. Ghost cells, sediment cells and       */
  /* surface cells below the top of the grid fall into this          */
  /* category.                                                       */
  memset(msk, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Initialise the mask for vaid wet cells. This mask is turned off */
  /* below for cells in the tridiagonal that are not wet cells.      */
  memset(mask, 0, window->sgsiz * sizeof(int));
  for (cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];
    mask[c] |= LF_WETM;
  }

  /*-----------------------------------------------------------------*/
  /* Redistribute the Aijs by adding weights associated with ghost   */
  /* cells to interior Aij.                                          */

  /* First check surface cells at the top of the grid.               */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    mask[c] |= LF_DESM;   /* Initialise all destination cells        */
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      /* Cells at the grid top                                       */
      if (j > nlay-1 && cs == window->zp1[cs]) {
	/* Put weights in cells above the surface into wet cells     */
	A[c][j-nlay] += A[c][j];
	A[c][j] = 0.0;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Next check for lateral ghost cells, i.e. if bin[] map lands on  */
  /* any of any of the Aijs.                                         */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    msk[c] = (double)window->bin[cc];
  }
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if ((ci = (int)msk[c1])) {
	for (k = 0; k < wincon->nosl; k++) {
	  if (ci == wincon->lmap[cs][k]) {
	    A[c][k] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    break;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Check the surface below the top of the grid. Note; the mask is  */
  /* reset for these cells.                                          */
  memset(msk, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = cs = wincon->s1[cc];
    zp1 = window->zp1[c];
    while (c != zp1) {
      c = zp1;
      zp1 = window->zp1[zp1];
      msk[c] = (double)cs;
    }
  }
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if ((ci = (int)msk[c1])) {
	for (k = 0; k < wincon->nosl; k++) {
	  if (ci == wincon->lmap[cs][k]) {
	    A[c][k] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    break;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Check for sediments and sediment ghost cells                    */
  memset(msk, 0, window->sgsiz * sizeof(double));
  /* Sediments                                                       */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->i1[cc];
    cs = window->zm1[c];
    msk[cs] = (double)c;
  }
  /* Sediment ghosts                                                 */
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    msk[c] = (double)window->ised_t[cc];
  }
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];

    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if ((ci = (int)msk[c1])) {
	int found = 1;
	for (k = 0; k < wincon->nosl; k++) {
	  if (ci == wincon->lmap[cs][k]) {
	    A[c][k] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    found = 0;
	    break;
	  }
	}
	/* Check the Aij above ci                                    */
	if (found && j < nlay) {
	  ci = window->zp1[c1];
	  if (!msk[ci]) {
	    A[c][j+nlay] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    found = 0;	    
	  }
	}
	/* Look for any wet cell Aij                                 */
	if (found) {
	  for (k = 0; k < wincon->nosl; k++) {
	    ci = wincon->lmap[cs][k];
	    if (!msk[ci]) {
	      A[c][k] += A[c][j];
	      A[c][j] = 0.0;
	      if (mask[c1] & LF_WETM)
		mask[c1] &= ~LF_WETM;
	      found = 0;
	      break;
	    }
	  }
	}
	/* Report and unassigned cells if required                   */
	if (found && check) {
	  hd_warn("reset_Aij: Cant find wet Aij at dest c=%d(%d %d %d), source cs=%d(%d %d %d), j=%d %f %f\n",c,
		  window->s2i[c],window->s2j[c],window->s2k[c],
		  cs,window->s2i[cs],window->s2j[cs],window->s2k[cs],j,windat->eta[window->m2d[cs]],wincon->dz[cs]);
	  printf("rest_Aij: Cant find wet Aij at dest c=%d(%d %d %d), source cs=%d(%d %d %d), j=%d %f %f\n",c,
		 window->s2i[c],window->s2j[c],window->s2k[c],
		 cs,window->s2i[cs],window->s2j[cs],window->s2k[cs],j,windat->eta[window->m2d[cs]],wincon->dz[cs]);
	  for (k = 0; k < wincon->nosl; k++) {
	    ci = wincon->lmap[cs][k];
	    printf(" k=%d ci=%d(%d %d %d) %x %f\n",k,ci,window->s2i[ci],window->s2j[ci],window->s2k[ci],(int)msk[ci],A[c][k]);
	  }
	  printf("source ");
	  s2ijk(window,cs);
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Check for dry c1: try to find the nearest wet cell              */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (!wincon->c1[c1]) {
	/* Search forwards from lmap[j]                              */
	for (k = j; k < wincon->nosl; k++) {
	  ci = wincon->lmap[cs][k];
	  if (wincon->c1[ci]) {
	    A[c][k] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    break;
	  }
	}
	/* Search backwards from lmap[j]                             */
	for (k = j; k >= 0; k--) {
	  ci = wincon->lmap[cs][k];
	  if (wincon->c1[ci]) {
	    A[c][k] += A[c][j];
	    A[c][j] = 0.0;
	    if (mask[c1] & LF_WETM)
	      mask[c1] &= ~LF_WETM;
	    break;
	  }
	}
      }
    }
  }

  /* Return at this stage for global fills (this routine may be      */
  /* called with global fills if REGIONS are invoked).               */
  if (!(wincon->fillf & (LOCAL|LOCALER))) return;

  /*-----------------------------------------------------------------*/
  /* Sanity check : sum(Aij) should equal 1 for LF_WETM cells. There */
  /* should be no non-valid cells (!LF_WETM) with nonzero Aij. If    */
  /* there are, then redistribute weights from non-valid cells to    */
  /* any valid wet cell (with mask & LF_WETM).                       */
  /* Note that Aij *= ovol/Vi in local_fill(), hence if ovol is zero */
  /* for all nonzero (valid) Aij then it is possible that Atot may   */
  /* still be zero in aij_*(), which results in division by zero. In */
  /* this case also redistribute. If no valid, nonzero Aij exit with */
  /* ovol != 0, then exclude the cell from the local fill.           */
  for (cc = 1; cc <= wincon->vc; cc++) {    
    c = wincon->s1[cc];
    cs = S[cc];
    Atot = Anot = Ack = 0.0;

    if(check) {
      if(isnanAij(window, c))
	hd_quit("reset_Aij: nan Aij at c=%d, %5.2f days\n", c, windat->days);
      if(!iswetAij(window, cc))
	hd_warn("reset_Aij: Pre no wet cells in Aij at c=%d, %5.2f days\n", c, windat->days);
      if(fabs(sumAij(window, c) - 1.0) > SMALL)
	hd_quit("reset_Aij: Pre Aij doesn't sum to 1 at c=%d, %5.2f days\n", c, windat->days);
    }

    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (mask[c1] & LF_WETM && ovol[c1])
	Atot += A[c][j];
    }
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      /* Sum the Aij for valid wet cells                             */
      if (!(mask[c1] & LF_WETM) || !ovol[c1]) {
	for (k = 0; k < wincon->nosl; k++) {
	  ci = wincon->lmap[cs][k];
	  if (mask[ci] & LF_WETM && ovol[ci]) {
	    A[c][k] += A[c][j];
	    Anot += A[c][j];
	    A[c][j] = 0.0;
	    break;
	  }
	}
      }
      Ack += A[c][j] * ovol[c1];
    }
    if (Anot && fabs(Atot+Anot-1.0) > SMALL)
      hd_warn("reset_Aij: Aij != 0 at non-valid cells, c=%d, Atot=%f Anot=%f at %5.2f days\n", c, Atot, Anot, windat->days);
    /* No valid cells found => exclude these cells                   */
    if (Ack == 0.0 || !iswetAij(window, cc)) {
      hd_warn("reset_Aij: can't find valid wet cell: excluding weights for c=%d at %5.2f days\n", c, windat->days);
      mask[c] &= ~LF_DESM;
    }
    if(check) {
      if(fabs(sumAij(window, c) - 1.0) > SMALL)
	hd_quit("reset_Aij: Post Aij doesn't sum to 1 at c=%d, %5.2f days\n", c, windat->days);
      if (!(mask[c] & LF_DESM) && !iswetAij(window, cc))
	hd_warn("reset_Aij: Post no wet cells in Aij at c=%d, %5.2f days\n", c, windat->days);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the masks                                                   */
  /* mask & LF_WETM : source cell is a valid wet cell                */
  /* mask & LF_OBCM : source cell is a non open boundary cell        */
  /* mask & LF_OBCF : all source cells are valid wet, non-obc cells  */
  /* mask & LF_AIJM : source cells with non-zero Aij, non-obc        */
  /* mask & LF_DESM : destination cell is volume corrected
  /* Valid wet cell mask, LF_WETM and LF_DESM is set in reset_Aij(). */
  /* Open boundary mask: set to LF_OBCM for non-boundary cells       */
  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    mask[c] |= LF_OBCM;
  }
  /* Wet source cells on the open boundary: set to LF_OBCF for valid */
  /* wet cells that are not on the boundary.                         */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    obcf = 1;
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (mask[c1] & LF_WETM && !(mask[c1] & LF_OBCM)) obcf = 0;
    }
    if (obcf) mask[c] |= LF_OBCF;
  }
}

/* END reset_Aij()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the volume error                                         */
/*-------------------------------------------------------------------*/
void get_verr(geometry_t *window,  /* Window geometry             */
	      window_t *windat,    /* Window data                 */
	      win_priv_t *wincon)  /* Window constants            */
{
  int c, cc, c1, c2, cs, j;
  double cellvol, d1;
  double *obc = wincon->w5;
  double **A = wincon->wgt;
  double *ovol = wincon->p1;
  int *mask = wincon->s3;
  int *S = wincon->s2;

  /* Open boundary mask                                              */
  /*
  memset(obc, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    if (mask[c] & LF_WETM) obc[c] = 1.0;
  }

  memset(windat->Vi, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    cellvol = window->cellarea[c2] * wincon->dz[c];
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      windat->Vi[c1] += (A[c][j] * cellvol * obc[c]);
    }
  }
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    windat->Vi[c] -= (ovol[c] * obc[c]);
  }
  */
  memset(windat->Vi, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cs = S[cc];
    cellvol = window->cellarea[c2] * wincon->dz[c];
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (mask[c1] & LF_WETM) {
	if (mask[c1] & LF_OBCM)
	  windat->Vi[c1] += (A[c][j] * cellvol);
	if (!(mask[c] & LF_OBCM))
	  mask[c1] |= LF_DESF;
      }
    }
  }
  /* set mask = LF_AIJM if cell meets criteria for vol corr             */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    if (mask[c] & LF_AIJM) mask[c] &= ~LF_AIJM;
    if (!(mask[c] & LF_DESF) && mask[c] & LF_OBCF && windat->Vi[c]) 
      mask[c] |= LF_AIJM;
  }
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    if (mask[c] & LF_AIJM)
      windat->Vi[c] -= ovol[c];
  }
}


/* END get_verr()                                                    */
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
    v1 = mxc[c] = max_lag_w(window, otr, c, cs);
    /* Maximum concentration due to source/sinks                     */
    /*mxc[c] = max(otr[c], v1);*/
    /* Minimum concentration due to advection                        */
    v1 = mnc[c] = min_lag_w(window, otr, c, cs);
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
    if (mask[c] & LF_AIJM)
      merr += otr[c] * windat->Vi[c];
    else
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
    if (mask[c] & LF_AIJM)
      merr += otr[c] * windat->Vi[c];
    else
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
  /* is capped to the maximum (minimum). Note: windat->Vi ~ 0 (hence */
  /* merr ~ 0) if the local fill terminated at exact volume          */
  /* conservation but sum(Aij) != 1.                                 */
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
  /* Diagnostics for global fill using local volume errors           */
  if (wincon->fillf & LOCALER && wincon->fillf & (DIAGNOSE|DIAGNOSE_BGC)) {
    if (mode == 1) {
      wincon->b1 = v1 = 0.0;
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	wincon->b1 += windat->Vi[c];
	v1 += window->cellarea[cs] * wincon->dz[c];
      }
      fprintf(tstrans.fp, "%f %f ", windat->days, 100.0 * wincon->b1 / v1);
    }
    if (wincon->fillf & (DIAGNOSE))
      fprintf(tstrans.fp, "%f %f ", 100.0 * merr / m1, msf);
    else
      fprintf(tstrans.fp, "%f %f ", merr, merr / m1);
    if (mode == 2) {
      fprintf(tstrans.fp, "\n");
      fflush(tstrans.fp);
    }
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

  /*-----------------------------------------------------------------*/
  /* Spread residual mass globally, where further residual mass is   */
  /* accounted for iteratively.                                      */
  if (wincon->fillf & (MONGLOB|LOCALER)) {
    m = 0;
    v1 = 1.0;
    msf = (tm) ? tme / tm : 1.0;
    /*memset(done, 0, window->sgsiz * sizeof(double));*/
    while (msf != 1.0 && v1 != 0.0 && m < itmax) {
      /* Sanity checks                                             */
      if (isnan(msf) || msf <= 0.0) {
	hd_warn("monotonic_fill: Tracer %s couldn't compute global scaling at %5.2f days (%f %f): no global filling\n",
		wincon->trname[n], windat->t/86400, tm, tme);
	return;
      }
      if (fabs(msf - 1.0) > msfmax) {
	hd_warn("monotonic_fill: Tracer %s global scaling %f considered too large at %5.2f days: no global filling\n",
		wincon->trname[n], msf, windat->t/86400);
	return;
      }
      /* Re-distribute the change in mass                          */
      v1 = tm = tmass = 0.0;
      for(cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	vol = window->cellarea[cs] * wincon->dz[c];
	v2 = windat->tr_wc[n][c] * msf;
	v3 = fabs(mxc[c] - mnc[c]);
	if (v3 && v3 < eps && !done[c]) {
	  v1 += (v2 - windat->tr_wc[n][c]) * vol;
	  done[c] = 1.0;
	} else if (v2 > mxc[c] && !done[c]) {
	  /* Check for adjusted mass > local maximum               */
	  v1 += (v2 - mxc[c]) * vol;
	  windat->tr_wc[n][c] = mxc[c];
	  done[c] = 1.0;
	} else if (v2 < mnc[c] && !done[c]) {
	  /* Check for adjusted mass < local minimum               */
	  v1 -= (mnc[c] - v2) * vol;
	  windat->tr_wc[n][c] = mnc[c];
	  done[c] = 1.0;
	} else if (!done[c]) {
	  /* Adjust the mass                                       */
	  windat->tr_wc[n][c] = v2;
	  tm += (windat->tr_wc[n][c] * vol);
	}
	/* Get the new toal mass (diagnostic only)                 */
	tmass += (windat->tr_wc[n][c] * vol);
      }
      /* Calculate the change in mass over the time-step.          */
      msf = 1.0;
      if (tm && !isnan(tm)) msf = (tm + v1) / tm;
      m++;
    }
    if (m >= itmax) 
      hd_warn("monotonic_fill: Tracer %s did not converge at %5.2f, residual mass = %5.3f\n",
	      wincon->trname[n], windat->t/86400, v1);
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
  int xm1 = window->xm1[c];
  int ym1 = window->ym1[c];
  int xmym1 = window->xm1[ym1];
  double v1, v2;
  v1 = max(tr[xmym1], max(tr[ym1], max(tr[c], tr[xm1])));
  c = window->zp1[c];
  xm1 = window->xm1[c];
  ym1 = window->ym1[c];
  xmym1 = window->xm1[ym1];
  v2 = max(tr[xmym1], max(tr[ym1], max(tr[c], tr[xm1])));
  return(max(v1, v2));
}
double max_lag_w(geometry_t *window, double *tr, int c, int cs)
{
  win_priv_t *wincon= window->wincon;
  int c1, j;
  double **A = wincon->wgt;
  double v1 = -1e10;

  for (j = 0; j < wincon->nosl; j++) {
    c1 = wincon->lmap[cs][j];
    if (A[c][j]) v1 = max(v1, tr[c1]);
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
  int xm1 = window->xm1[c];
  int ym1 = window->ym1[c];
  int xmym1 = window->xm1[ym1];
  double v1, v2;
  v1 = min(tr[xmym1], min(tr[ym1], min(tr[c], tr[xm1])));
  c = window->zm1[c];
  xm1 = window->xm1[c];
  ym1 = window->ym1[c];
  xmym1 = window->xm1[ym1];
  v2 = min(tr[xmym1], min(tr[ym1], min(tr[c], tr[xm1])));
  return(min(v1, v2));
}
double min_lag_w(geometry_t *window, double *tr, int c, int cs)
{
  win_priv_t *wincon= window->wincon;
  int c1, j;
  double **A = wincon->wgt;
  double v1 = 1e10;

  for (j = 0; j < wincon->nosl; j++) {
    c1 = wincon->lmap[cs][j];
    if (A[c][j]) v1 = min(v1, tr[c1]);
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
  double t, sc;
  int ci, cp, cm;
  int *ct, *cp1, *cm1, *cv;
  double *tr = windat->tr_wc[n];
  double *vel;
  double dobc = 0.0;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    if(!wincon->trinfo_3d[n].advect) continue; 
    sc = 1.0;
    if (vecf) {
      ct = open->oi1_t;
      cp1 = open->oi2_t;
      cm1 = open->obc_t;
    } else {
      ct = open->obc_t;
      cp1 = open->oi1_t;
      cm1 = open->obc_t;
    }
    if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) sc = -1.0;
    if (ord == 1) {
      vel = (open->ocodex & (L_EDGE | R_EDGE)) ? windat->u1 : windat->u2;
      cv = (vecf) ? open->oi1_e1 : open->obc_e1;
    }
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = ct[cc];
      cp = cp1[cc];
      cm = cm1[cc];
      if (ord == 1) {          /* 1st order                          */
	ci = cv[cc];
	t = (sc * vel[ci] > 0.0) ? tr[cm] : tr[c];
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
  int c, cc, cs, nn;
  double *vel;
  double *hat, *dz, sc;
  int *cv;
  double dobc = 0.0;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    open->bflux_3d = 0.0;
    sc = 1.0;
    if (open->ocodex & (L_EDGE | R_EDGE)) {
      vel = windat->u1;
      hat = window->h2au1;
      dz = windat->dzu1;
      cv = (vecf) ? open->oi1_e1 : open->obc_e1;
    } else {
      vel = windat->u2;
      hat = window->h1au2;
      dz = windat->dzu2;
      cv = (vecf) ? open->oi1_e2 : open->obc_e2;
    }
    if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) sc = -1.0;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = cv[cc];
      cs = window->m2d[c];
      open->dum[cc] = (sc * vel[c] * hat[cs] * dz[c] * windat->dttr);
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
    if (open->ocodex & (L_EDGE | R_EDGE)) {
      cv = (vecf) ? open->oi1_e1 : open->obc_e1;
    } else {
      cv = (vecf) ? open->oi1_e2 : open->obc_e2;
    }
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
  int c, cc, cs;
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
    int sc = 1.0, *cv;
    if (open->ocodex & (L_EDGE | R_EDGE)) {
      vel = windat->u1;
      hat = window->h2au1;
      cv = open->obc_e1;
    } else {
      vel = windat->u2;
      hat = window->h1au2;
      cv = open->obc_e2;
    }
    if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) sc = -1.0;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      v = vel[cv[cc]];
      cs = window->m2d[cv[cc]];
      vol = sc * v * hat[cs] * wincon->dz[c] * windat->dttr;
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
  double *f1, *f2;              /* e1 and e2 transports              */
  double colflux;               /* Flux divergence                   */
  int c, cc, cs;                /* Sparse coodinate / counter        */
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
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    f1[cs] += windat->u1flux3d[c] * windat->dttr;
  }
  /* Integrate the e2 3D fluxes through the water column.            */
  f2 = wincon->d2;
  memset(f2, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    f2[cs] += windat->u2flux3d[c] * windat->dttr;
  }

  /* Calculate the new elevation                                     */
  set_map_eta(window);
  for (cc = 1; cc <= window->b2_t; cc++) {
    double eta;
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = f1[window->xp1[c]] - f1[c] + f2[window->yp1[c]] - f2[c];

    /* Add source/sink flows to colflux. Note that colflux is the    */
    /* accumulated volume flowing out of this column, whereas        */
    /* waterss2d is the flow (m3 s-1) into the column. Hence the     */
    /* minus sign.  */
    /* colflux -= owaterss2d[c] * windat->dttr; */
    colflux -= windat->waterss2d[c] * windat->dttr;

    /* Calculate new etat value                                      */
    /* windat->eta[c] = max(windat->etab[c] - colflux / window->cellarea[c],
       window->botz[c]); */

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
  double *f1, *f2;              /* e1 and e2 transports              */
  double colflux;               /* Flux divergence                   */
  int c, cc, cs;                /* Sparse coodinate / counter        */
  double eta;                   /* Computed sea level                */
  double eps = 1e-5;            /* Tolerance                         */

  /*-----------------------------------------------------------------*/
  /* Integrate the e1 3D fluxes through the water column.            */
  f1 = wincon->d1;
  memset(f1, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    f1[cs] += windat->u1flux3d[c] * windat->dttr;
  }
  /* Integrate the e2 3D fluxes through the water column.            */
  f2 = wincon->d2;
  memset(f2, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    f2[cs] += windat->u2flux3d[c] * windat->dttr;
  }

  /* Calculate the new elevation                                     */
  set_map_eta(window);
  for (cc = 1; cc <= window->b2_t; cc++) {
    double eta;
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = f1[window->xp1[c]] - f1[c] + f2[window->yp1[c]] - f2[c];

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
    if (master->fillf & LOCAL)
      fprintf(tstrans.fp, "## COLUMNS 9\n");
    else if (master->fillf & LOCALER)
      fprintf(tstrans.fp, "## COLUMNS %d\n", 2 + 2 * ntr);
    else
      fprintf(tstrans.fp, "## COLUMNS %d\n", 3 + 4 * ntr);
    fprintf(tstrans.fp, "##\n");
    fprintf(tstrans.fp, "## COLUMN1.name           Time\n");
    fprintf(tstrans.fp, "## COLUMN1.long_name      Time\n");
    fprintf(tstrans.fp, 
	    "## COLUMN1.units          days since 1990-01-01 00:00:00 +10\n");
    fprintf(tstrans.fp, "## COLUMN1.missing_value  -999\n");
    fprintf(tstrans.fp, "##\n");

    if (master->fillf & LOCAL) {
      fprintf(tstrans.fp, "## COLUMN2.name           vcorr\n");
      fprintf(tstrans.fp, "## COLUMN2.long_name      Maxmim volume error\n");
      
      fprintf(tstrans.fp, "## COLUMN2.units          %%\n");
      fprintf(tstrans.fp, "## COLUMN2.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN3.name           acorr\n");
      fprintf(tstrans.fp, "## COLUMN3.long_name      Maximum Aij error\n");
      fprintf(tstrans.fp, "## COLUMN3.units          %%\n");
      fprintf(tstrans.fp, "## COLUMN3.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN4.name           niter\n");
      fprintf(tstrans.fp, "## COLUMN4.long_name      Number of iterations\n");
      fprintf(tstrans.fp, "## COLUMN4.units          \n");
      fprintf(tstrans.fp, "## COLUMN4.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN5.name           verr_d\n");
      fprintf(tstrans.fp, "## COLUMN5.long_name      Diagnosed volume error\n");
      fprintf(tstrans.fp, "## COLUMN5.units          m^3\n");
      fprintf(tstrans.fp, "## COLUMN5.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN6.name           verr_c\n");
      fprintf(tstrans.fp, "## COLUMN6.long_name      Computed volume error\n");
      fprintf(tstrans.fp, "## COLUMN6.units          m^3\n");
      fprintf(tstrans.fp, "## COLUMN6.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN7.name           vdiff\n");
      fprintf(tstrans.fp, "## COLUMN7.long_name      Volume difference over timestep\n");
      fprintf(tstrans.fp, "## COLUMN7.units          m^3\n");
      fprintf(tstrans.fp, "## COLUMN7.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN8.name           vt-1\n");
      fprintf(tstrans.fp, "## COLUMN8.long_name      Volume at time t-1\n");
      fprintf(tstrans.fp, "## COLUMN8.units          m^3\n");
      fprintf(tstrans.fp, "## COLUMN8.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      fprintf(tstrans.fp, "## COLUMN9.name           vout\n");
      fprintf(tstrans.fp, "## COLUMN9.long_name      Total outflow volume\n");
      fprintf(tstrans.fp, "## COLUMN9.units          m^3\n");
      fprintf(tstrans.fp, "## COLUMN9.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
    } else {
      fprintf(tstrans.fp, "## COLUMN2.name           dvol\n");
      fprintf(tstrans.fp, "## COLUMN2.long_name      Volume error as %% of total volume\n");
      
      fprintf(tstrans.fp, "## COLUMN2.units          %%\n");
      fprintf(tstrans.fp, "## COLUMN2.missing_value  0.000000\n");
      fprintf(tstrans.fp, "##\n");
      m = 3;
      if (!(master->fillf & LOCALER)) {
	fprintf(tstrans.fp, "## COLUMN3.name           msf_OBC\n");
	fprintf(tstrans.fp, "## COLUMN3.long_name      Open boundary scaling\n");
	fprintf(tstrans.fp, "## COLUMN3.units            \n");
	fprintf(tstrans.fp, "## COLUMN3.missing_value  0.000000\n");
	fprintf(tstrans.fp, "##\n");
	m++;
      }
      for (nn = 0; nn < ntr; nn++) {
	n = wincon[1]->tbdy[nn];
	if((trn = tracer_find_index(master->trname[n], master->ntr,
				    master->trinfo_3d)) >= 0) {
	  if (!(master->fillf & LOCALER)) {
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
	  }
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
	  
	  if (!(master->fillf & LOCALER)) {
	    fprintf(tstrans.fp, "## COLUMN%d.name           %s_input\n", m, 
		    master->trname[n]);
	    fprintf(tstrans.fp, "## COLUMN%d.long_name      Mass pss input for %s\n", m,
		    master->trname[n]);
	    fprintf(tstrans.fp, "## COLUMN%d.units          kg\n", m);
	    fprintf(tstrans.fp, "## COLUMN%d.missing_value  0.000000\n", m);
	    fprintf(tstrans.fp, "##\n");
	    m++;
	  }
	} else
	  hd_quit("Can't recognize tracer %s for totals.\n", master->trname[n]);
      }
    }
  }
}

/* END init_trans()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check that at least one of the cells in the tri-       */
/* diagonal matrix is wet. If not, this means the streamline origin  */
/* has landed in a dry cell, and in this case move it to the centre  */
/* of the nearest wet cell.                                          */
/*-------------------------------------------------------------------*/
void sl_check(geometry_t *window, int *c, double *cx, double *cy, double *cz)
{
  int **map = window->wincon->lmap;
  int c1, cs, n, fl = 1;
  int m = window->wincon->nosl;

  cs = *c;
  for (n = 0; n < m; n++) {
    c1 = map[cs][n];
    if (!window->wgst[c1]) fl = 0;
  }
  if (fl) {
    *c = window->wgst[cs];
    *cx = 0.0;
    *cy = 0.0;
    if (cs == window->zm1[cs]) *cz = 0.0;
  }
}

/* END sl_check()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to sum the weights                                        */
/*-------------------------------------------------------------------*/
double sumAij(geometry_t *window, int c)
{
  win_priv_t *wincon = window->wincon;
  double **A = wincon->wgt;
  double d1 = 0.0;
  int j;

  for (j = 0; j < wincon->nosl; j++)
    d1 += A[c][j];

  return(d1);
}

int iswetAij(geometry_t *window, int cc)
{
  win_priv_t *wincon = window->wincon;
  int *S = wincon->s2; 
  int ret = 0;
  int cs, c1, j;

  cs = S[cc];
  for (j = 0; j < wincon->nosl; j++) {
    c1 = wincon->lmap[cs][j];
    ret = max(ret, wincon->c1[c1]);
  }
  return(ret);
}


int isnanAij(geometry_t *window, int c)
{
  win_priv_t *wincon = window->wincon;
  double **A = wincon->wgt;
  int j;

  for (j = 0; j < wincon->nosl; j++) {
    if (isnan(A[c][j]))
      return(1);
  }
  return(0);
}

      
/* END sumAij()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void check_transport(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int c, cc, cs, c1, j;
  int *S = wincon->s2;
  double **A= wincon->wgt;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = S[cc];
    if (wincon->dz[c] <= 0)
      hd_quit("transport check: dz (%5.2f) <= 0 at %c(%d %d %d)\n",wincon->dz[c],
	      c, window->s2i[c], window->s2j[c], window->s2k[c]);
    if (cs <=0 || cs > window->enon)
      hd_quit("transport check: invalid source cell %d at %c(%d %d %d)\n", cs,
	      c, window->s2i[c], window->s2j[c], window->s2k[c]);
    for (j = 0; j < wincon->nosl; j++) {
      c1 = wincon->lmap[cs][j];
      if (c1 <=0 || c1 > window->enon)
	hd_quit("transport check: invalid tridiagonal cell %d, j=%d at destination cell %c(%d %d %d)\n", 
		c1, j, c, window->s2i[c], window->s2j[c], window->s2k[c]);
      if (A[c][j] < 0 || A[c][j] > 1)
	hd_warn("transport check: invalid tridiagonal weight %5.2, j=%d at destination cell %c(%d %d %d)\n",
		A[c][j], j, c, window->s2i[c], window->s2j[c], window->s2k[c]);
    }
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
  int cc, c;

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
      if (isnan(tpd->u1[c])) tpd->u1[c] = 0.0;
      if (isnan(tpd->u2[c])) tpd->u2[c] = 0.0;
      if (isnan(tpd->w[c])) tpd->w[c] = 0.0;
    }
  } else {
    for (cc = 1; cc <= tpg->n2_t; cc++) {
      c = tpg->w2_t[cc];
      if (isnan(tpd->eta[c]) || fabs(tpd->eta[c]) > tpc->etamax)
	tpd->eta[c] = 0.0;
    }
    for (cc = 1; cc <= tpg->n3_t; cc++) {
      c = tpg->w3_t[cc];
      if (isnan(tpd->u1[c]) || fabs(tpd->u1[c]) > tpc->velmax)
	tpd->u1[c] = 0.0;
      if (isnan(tpd->u2[c]) || fabs(tpd->u2[c]) > tpc->velmax)
	tpd->u2[c] = 0.0;
      if (isnan(tpd->w[c]) || fabs(tpd->w[c]) > tpc->velmax)
	tpd->w[c] = 0.0;
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
  int cc, c, cp;

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
    for (cc = 1; cc <= tpg->n3_t; cc++) {
      c = tpg->w3_t[cc];
      if (isnan(tpd->u1[c]) || fabs(tpd->u1[c]) > tpc->velmax)
	check_warn(tpg, "u1 velocity", tpd->u1[c], master->days, c, 3);
      if (isnan(tpd->u2[c]) || fabs(tpd->u2[c]) > tpc->velmax)
	check_warn(tpg, "u2 velocity", tpd->u2[c], master->days, c, 3);
      cp = tpg->xp1[c];
      if (isnan(tpd->u1[cp]) || fabs(tpd->u1[cp]) > tpc->velmax)
	check_warn(tpg, "u1 velocity", tpd->u1[cp], master->days, cp, 3);
      cp = tpg->yp1[c];
      if (isnan(tpd->u2[cp]) || fabs(tpd->u2[cp]) > tpc->velmax)
	check_warn(tpg, "u2 velocity", tpd->u2[cp], master->days, cp, 3);
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

void print_tridiagonal(geometry_t *window, double *tr, int c, int cs)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  printf("%4.2f : %f %f %f %f %f %f %f %f : %d(%d %d %d)\n", windat->days, 
	 wincon->wgt[c][0], wincon->wgt[c][1], wincon->wgt[c][2],
	 wincon->wgt[c][3], wincon->wgt[c][4], wincon->wgt[c][5],
	 wincon->wgt[c][6], wincon->wgt[c][7],
	 cs, window->s2i[cs], window->s2j[cs], window->s2k[cs]);
  if(wincon->togn & TOPRIGHT) {
    printf("          %f %f %f %f %f %f %f %f\n", tr[cs], 
	   tr[window->ym1[cs]], tr[window->xm1[cs]],
	   tr[window->xm1[window->ym1[cs]]], tr[window->zp1[cs]], 
	   tr[window->ym1[window->zp1[cs]]],
	   tr[window->xm1[window->zp1[cs]]], 
	   tr[window->ym1[window->xm1[window->zp1[cs]]]]);
    printf("          %d %d %d %d %d %d %d %d\n", cs, window->ym1[cs], 
	   window->xm1[cs], window->xm1[window->ym1[cs]], 
	   window->zp1[cs], window->ym1[window->zp1[cs]],
	   window->xm1[window->zp1[cs]], 
	   window->ym1[window->xm1[window->zp1[cs]]]);
  } else {
    printf("          %f %f %f %f %f %f %f %f\n", tr[cs], 
	   tr[window->yp1[cs]], tr[window->xp1[cs]],
	   tr[window->xp1[window->yp1[cs]]], tr[window->zp1[cs]], 
	   tr[window->yp1[window->zp1[cs]]],
	   tr[window->xp1[window->zp1[cs]]], 
	   tr[window->yp1[window->xp1[window->zp1[cs]]]]);
  }
}

void check_tridiagonal(geometry_t *window, double *tr, char *text, int c, int cs)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int in[8];
  int n, trn;
  double t;

  for (trn = 0; trn < windat->ntr; trn++)
    if (tr == windat->tr_wc[trn]) break;
  in[0] = cs;
  in[1] = window->ym1[cs];
  in[2] = window->xm1[cs];
  in[3] = window->xm1[window->ym1[cs]];
  in[4] = window->zp1[cs]; 
  in[5] = window->ym1[window->zp1[cs]];
  in[6] = window->xm1[window->zp1[cs]];
  in[7] = window->ym1[window->xm1[window->zp1[cs]]];

  for(n=0; n<8; n++) {
    t = tr[in[n]];
    if(isnan(wincon->wgt[c][n]))
      hd_quit("%s: wgt nan error at %d %d=%f\n",text,c,n,wincon->wgt[c][n]);
    if(wincon->wgt[c][n] < 0.0 || wincon->wgt[c][n] > 1.0)
      hd_warn("%s: wgt bounds error at %d %d=%f\n",text,c,n,wincon->wgt[c][n]);
    if (isnan(t))
      hd_quit("%s: Tracer nan error at %d %d=%f\n",text,in[n],n,t);
    if(t < wincon->mintr[trn] || t > wincon->maxtr[trn])
      hd_warn("%s: Tracer bounds error at %d %d=%f\n",text,in[n],n,t);
  }    
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the weights for interpolations of order 2 & 4. Uses Eqn. 5   */ 
/* McDonald (1984) Mon. Wea. Rev. 112, 1267-1275. The distances are  */
/* re-mapped to local distances, where h[0] = 0. Note that (Eqn. 8)  */
/* hx[i-1/2] < p < h[i+1/2], hy[j-1/2] < q < hy[j+1/2] and           */
/* hz[i-1/2] < r < hz[i+1/2].                                        */
/* Note: weights can be negative.                                    */
/*-------------------------------------------------------------------*/
void weights_eo(geometry_t *window, 
		int c,             /* Destination cell               */
		int co,            /* Source cell                    */
		double pin,        /* x distance from dest origin    */
		double qin,        /* y distance from dest origin    */
		double rin,        /* z distance from dest origin    */
		double *dzz)       /* Distance between layer faces   */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj, ck;
  double p, q, r, d;
  double hx[5], hy[5], hz[5];

  /* Set up the normalized distances */
  ck = co;
  ci = cj = window->m2d[co];
  p = window->h1au1[ci] - pin;
  q = window->h2au2[ci] - qin;
  r = dzz[ck] + rin;
  ci = window->xm1[ci];
  cj = window->ym1[cj];
  ck = window->zm1[ck];
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    ck = window->zm1[ck];
    p += window->h1au1[ci];
    q += window->h2au2[cj];
    r += dzz[ck];
  }
  hx[0] = hy[0] = hz[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + window->h1au1[ci];
    hy[m] = hy[m-1] + window->h2au2[cj];
    hz[m] = hz[m-1] + dzz[ck];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
    ck = window->zp1[ck];
  }

  /* Loop to get the weights */
  in = 0;
  d = 0.0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	wincon->wgt[c][in] = 1.0;
	for (m = 0; m <= n; m++) {
	  if (m != k)
	    wincon->wgt[c][in] *= (r - hz[m]) / (hz[k] - hz[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != j)
	    wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != i)
	    wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	}
	d += wincon->wgt[c][in];
	in++;
      }
  }
  /*
  if (fabs(d - 1.0) > SMALL) hd_warn("Weights do not sum to 1 at %d : %e\n",c, d);
  */
}

/* END weights_eo()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the weights for interpolations of order 1 & 3. Uses Eqn. 5   */ 
/* McDonald (1984) Mon. Wea. Rev. 112, 1267-1275. The distances are  */
/* re-mapped to local distances, where h[0] = 0. Note that (Eqn. 7)  */
/* hx[i-1] < p < h[i], hy[j-1] < q < hy[j] and hz[i] < r < hz[i+1].  */
/*-------------------------------------------------------------------*/
void weights_oo(geometry_t *window, 
		int c,             /* Destination cell               */
		int co,            /* Source cell                    */
		double pin,        /* x distance from dest origin    */
		double qin,        /* y distance from dest origin    */
		double rin,        /* z distance from dest origin    */
		double *dzz)       /* Distance between layer centres */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj, ck;
  double p, q, r;
  double hx[5], hy[5], hz[5];

  /* Set up the normalized distances */
  ck = co;
  ci = cj = window->m2d[co];
  p = window->h1au1[ci] - pin;
  q = window->h2au2[cj] - qin;
  r = rin;
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    ck = window->zm1[ck];
    p += window->h1au1[ci];
    q += window->h2au2[cj];
    r += dzz[ck];
  }
  hx[0] = hy[0] = hz[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + window->h1au1[ci];
    hy[m] = hy[m-1] + window->h2au2[cj];
    hz[m] = hz[m-1] + dzz[ck];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
    ck = window->zp1[ck];
  }

  /* Loop to get the weights */
  in = 0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	wincon->wgt[c][in] = 1.0;
	for (m = 0; m <= n; m++) {
	  if (m != k)
	    wincon->wgt[c][in] *= (r - hz[m]) / (hz[k] - hz[m]);
	}

	for (m = 0; m <= n; m++) {
	  if (m != j)
	    wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != i)
	    wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	}
	in++;
      }
  }
}

/* END weights_oo()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the interpolation on array in[] of order n, using the    */
/* weights wgt[] and mappings in lmap[].                             */
/*-------------------------------------------------------------------*/
double int_valo(geometry_t *window, double *wgt, int *lmap, double *in, int n)
{
  double out;
  int m;
  int nosl = window->wincon->nosl;

  out = 0.0;
  for (m = 0; m < nosl; m++) {
    out += wgt[m] * in[lmap[m]];
  }
  return (out);
}

/* END int_valo()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* The origin is in the front right corner horizontally, and bottom  */
/* right vertically. The indicies run from:                          */
/* O1. i&j :  c-1, c                | k : c, c+1                     */
/* O2.        c-1, c, c+1           |    c-1, c, c+1                 */
/* O3.        c-2, c-1, c, c+1      |    c-1, c, c+1, c+2            */
/* O4.        c-2, c-1, c, c+1, c+2 |    c-2, c-1, c, c+1, c+2       */
/* See McDonald (1984) Mon. Wea. Rev. 112, 1267-1275.                */
/* This routine only needs be done once at initialisation.           */
/*-------------------------------------------------------------------*/
void set_lmap(geometry_t *window, win_priv_t *wincon) {
  int c, co, i, j, k, m;
  int n = wincon->osl;
  int size;
  int ci, cj, ck;
  int ns, ne, nk;
  int mi, mj, mk;
  int kk, ckk;
  int klay[5] = {0, 0, 1, 1, 2};

  size = wincon->nosl = sls[n];
  wincon->mosl = klay[n];

  if (n == 0) {
    /* Original formulation                                          */
    wincon->wgt = d_alloc_2d(8, window->sgsiz);
    wincon->lmap = i_alloc_2d(8, window->sgsiz);

    for(c = 1; c <= window->enon; c++) {
      if (master->togn & TOPRIGHT) {
	wincon->lmap[c][0] = c;
	wincon->lmap[c][1] = window->ym1[c];
	wincon->lmap[c][2] = window->xm1[c];
	wincon->lmap[c][3] = window->xm1[window->ym1[c]];
	wincon->lmap[c][4] = window->zp1[c];
	wincon->lmap[c][5] = window->ym1[window->zp1[c]];
	wincon->lmap[c][6] = window->xm1[window->zp1[c]];
	wincon->lmap[c][7] = window->ym1[window->xm1[window->zp1[c]]];
      } else {
	wincon->lmap[c][0] = c;
	wincon->lmap[c][1] = window->yp1[c];
	wincon->lmap[c][2] = window->xp1[c];
	wincon->lmap[c][3] = window->xp1[window->yp1[c]];
	wincon->lmap[c][4] = window->zp1[c];
	wincon->lmap[c][5] = window->yp1[window->zp1[c]];
	wincon->lmap[c][6] = window->xp1[window->zp1[c]];
	wincon->lmap[c][7] = window->yp1[window->xp1[window->zp1[c]]];
      }
    }
    return;
  }

  /* Allocate memory                                                 */
  wincon->wgt = d_alloc_2d(size, window->sgsiz);
  wincon->lmap = i_alloc_2d(size, window->sgsiz);

  /* Get the number of cells to the start cell                       */
  ne = (n > 2) ? 2 : 1;
  nk = (n == 1 || n == 3) ? ne - 1 : ne;

  /* Get the maps for each origin                                    */
  /* Note: if a cell is self mapping backwards, then remember this   */
  /* so it does not map forward (mi,mj,mk) times.                    */
  for(c = 1; c <= window->enon; c++) {
    ci = co = c;
    m = mi = 0;
    for (ns = 0; ns < ne; ns++) {
      ci = window->xm1[ci];
      if (ci == co) mi++;
      co = ci;
    }
    for (i = 0; i <= n; i++) {
      cj = co = ci;
      mj = 0;
      for (ns = 0; ns < ne; ns++) {
	cj = window->ym1[cj];
	if (cj == co) mj++;
	co = cj;
      }
      for (j = 0; j <= n; j++) {
	ck = co = cj;
	mk = 0;
	for (ns = 0; ns < nk; ns++) {
	  ck = window->zm1[ck];
	  if (ck == co) mk++;
	  co = ck;
	}
	for (k = 0; k <= n; k++) {
	  ckk = ck;
	  /* Use the mirror image at boundaries                      */
	  if (wincon->momsc & LAGRANGE) {
	    /* Sediments                                             */
	    if (ck == window->zm1[ck]) {
	      for (kk = 0; kk < 2 - mk; kk++)
		ckk = window->zp1[ckk];
	    }
	    /* Surface                                               */
	    if (k > nk && ck == window->zp1[ck]) {  
	      /*for (kk = 0; kk <= n - k; kk++)*/
	      for (kk = 0; kk < k - 1; kk++)
		ckk = window->zm1[ckk];
	    }
	  }
	  wincon->lmap[c][m] = ckk;
	  m++;
	  if (k >= mk)
	    ck = window->zp1[ck];
	}
	if (j >= mj)
	  cj = window->yp1[cj];
      }
      if (i >= mi)
	ci = window->xp1[ci];
    }
  }
}

/* END set_lmap()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the origin  for tri-quadratic and tri-quartic    */
/* interpolations. See McDonald (1984) Mon. Wea. Rev. 112, 1267-1275 */
/* Eqn. 8.                                                           */ 
/*-------------------------------------------------------------------*/
void set_eo(geometry_t *window, win_priv_t *wincon, double *dzz,
	    int *c, double *cx, double *cy, double *cz)
{
  int co = *c;
  int cs = window->m2d[co];

  if (*cx > 0.5 * window->h1acell[cs]) {
    *cx -= window->h1au1[cs];
    co = window->xm1[co];
    cs = window->m2d[co];
  }
  if (*cy > 0.5 * window->h2acell[cs]) {
    *cy -= window->h2au2[cs];
    co = window->ym1[co];
  }
  if (*cz > 0.5 * wincon->dz[co]) {
    *cz -= dzz[co];
    co = window->zp1[co];
  }
  *c = co;
}

/* END set_eo()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to enforce monotonicity. The value from the interpolation */
/* cannot be more/less than the max/min of contributing values and   */
/* any sourcesink input.                                             */
/*-------------------------------------------------------------------*/
void set_monotonic(int n, int co, double *ntr, double *tr, int *lmap){

  int m, c;
  double mn = 1e10, mx = -1e10;

  for (m = 0; m < sls[n]; m++) {
    c = lmap[m];
    mx = max(tr[c], mx);
    mn = min(tr[c], mn);
  }
  mx = max(mx, tr[co]);
  mn = min(mn, tr[co]);
  *ntr = max(min(*ntr, mx), mn);
}

/* END set_monotonic()                                               */
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
  int cc, c, c2, tn, nn, n;
 
  /* FRONT and RIGHT open boundaries                               */
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	c2 = open->oi1_e1[cc];
	for (n = 0; n < wincon->ntbdy; n++) {
	  tn = wincon->tbdy[n];
	  if (!(open->bcond_tra[tn] & (TRCONC|TRCONF|TRFLUX)))
	    tr[tn][c] = tr[tn][c2];
	}
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	c2 = open->oi1_e2[cc];
	for (n = 0; n < wincon->ntbdy; n++) {
	  tn = wincon->tbdy[n];
	  if (!(open->bcond_tra[tn] & (TRCONC|TRCONF|TRFLUX)))
	    tr[tn][c] = tr[tn][c2];
	}
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
    if (windat->u1m) {
      for (cc = 1; cc <= window->b3_e1; cc++) {
	c = window->w3_e1[cc];
	cs = window->m2d[c];
	windat->u1m[c] = (windat->u1m[c] * windat->meanc[cs] + 
			  windat->u1[c] * t) / (windat->meanc[cs] + t);
      }
    }
    if (windat->u2m) {
      for (cc = 1; cc <= window->b3_e2; cc++) {
	c = window->w3_e2[cc];
	cs = window->m2d[c];
	windat->u2m[c] = (windat->u2m[c] * windat->meanc[cs] + 
			  windat->u2[c] * t) / (windat->meanc[cs] + t);
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
      for (cc = 1; cc <= window->b2_e1; cc++) {
	c = window->w2_e1[cc];
	windat->u1am[c] = (windat->u1am[c] * windat->meanc[c] + 
			  windat->u1av[c] * t) / (windat->meanc[c] + t);
      }
    }
    if (windat->u2am) {
      for (cc = 1; cc <= window->b2_e2; cc++) {
	c = window->w2_e2[cc];
	windat->u2am[c] = (windat->u2am[c] * windat->meanc[c] + 
			  windat->u2av[c] * t) / (windat->meanc[c] + t);
      }
    }
  }
}

/* END calc_means_t()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Rotates velocity components to east & north                       */
/*-------------------------------------------------------------------*/
void rotate_vel(geometry_t *window, double *u, double *v)
{
  window_t *windat = window->windat;
  int c, cc, cs, xp1, yp1;
  double u1val, u2val, sinth, costh;

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    u1val = 0.5 * (windat->u1[c] + windat->u1[xp1]);
    u2val = 0.5 * (windat->u2[c] + windat->u2[yp1]);
    sinth = window->sinthcell[cs];
    costh = window->costhcell[cs];
    u[c] = u1val * costh - u2val * sinth;
    v[c] = u1val * sinth + u2val * costh;
  }
}

/* END rotate_vel()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Rotates east & north velocity components onto the grid            */
/*-------------------------------------------------------------------*/
void rerotate_vel(geometry_t *window, double *u, double *v)
{
  window_t *windat = window->windat;
  int c, cc, cs, xm1, ym1;
  double uval, vval, sinth, costh;

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    sinth = window->sinthcell[cs];
    costh = window->costhcell[cs];
    windat->u1[c] = u[c] * costh + v[c] * sinth;
    windat->u2[c] = -u[c] * sinth + v[c] * costh;
    /* Face centered
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    uval = 0.5 * (u[c] + u[xm1]);
    vval = 0.5 * (v[c] + v[ym1]);
    windat->u1[c] = uval * costh + vval * sinth;
    windat->u2[c] = -uval * sinth + vval * costh;
    */
  }
}

/* END rerotate_vel()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-lagrangian advection scheme. See Casulli and Cheng (1992),   */
/* Int. J. Num. Meth. Fluids, 15, 629-648.                           */
/* This part is performed before the tracer loop to get the grid     */
/* Courant numbers. Used for single grids (source =  target).        */
/*-------------------------------------------------------------------*/
void streamline_atc(geometry_t *window,  /* Processing window     */
		    window_t *windat,    /* Window data structure */
		    win_priv_t *wincon,  /* Window constants      */
		    int c                /* Location to trace     */
  )
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->nu;
  double *nv = wincon->nv;
  double *nw = wincon->nw;
  double cx;
  double cy;
  double cz;
  double *dzz = wincon->w9;
  double *mask = wincon->d3;
  double u, v, w;
  double p, q, r;
  int ci, c2, zm1, cc;
  int *cl = wincon->s2;        /* Streamline origin */

  if (wincon->trasc != FFSL) return;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  cx = cy = cz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  cg = cl[c] = c;
  time_left = windat->dttr;
  u = nu[c];
  v = nv[c];
  w = nw[c];

  while (time_left > 0) {

    int ins = (cl[c] == wincon->m2d[cl[c]]) ? 1 : 0;
    /* Get the sub-timestep for the next integration. Streamlines    */
    /* cannot cross more than one cell in one sub-timestep.          */
    c2 = window->m2d[cl[c]];
    ci = (nu[cl[c]] > 0.0) ? c2 : window->xp1[c2];
    p = window->h1au1[ci] * mask[ci];
    ci = (nv[cl[c]] > 0.0) ? c2 : window->yp1[c2];
    q = window->h2au2[ci] * mask[ci];
    r = (nw[cl[c]] < 0.0) ? dzz[cl[c]] : dzz[window->zm1[cl[c]]];
    dt = min(SCALE * fabs(p / nu[cl[c]]),
	     SCALE * fabs(q / nv[cl[c]]));
    dt = (r) ? ((ins && nw[cl[c]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[c]]))) : dt;
    dt = min(dt, time_left);

    /* Get the bottom left grid coordinate defining the streamline   */
    /* position.                                                     */
    cl[c] = get_pos(window, cl[c], u, v, w, &cx, &cy, &cz, dt);

    c2 = window->m2d[cl[c]];
    p = cx / window->h1au1[c2];
    q = cy / window->h2au2[c2];
    r = (dzz[cl[c]]) ? cz / dzz[cl[c]] : 0.0;

    wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
    wgt[1] = (1.0 - r) * (1.0 - p) * q;
    wgt[2] = (1.0 - r) * p * (1.0 - q);
    wgt[3] = (1.0 - r) * p * q;
    wgt[4] = r * (1.0 - p) * (1.0 - q);
    wgt[5] = r * (1.0 - p) * q;
    wgt[6] = r * p * (1.0 - q);
    wgt[7] = r * p * q;

    /* Interpolate x,y,z velocities */
    u = int_val(window, cl[c], wgt, nu);
    v = int_val(window, cl[c], wgt, nv);
    w = int_val(window, cl[c], wgt, nw);

    /* Get the number of sub-timesteps used */
    time_left = time_left - dt;
  }

  /*if (wincon->osl == 0) {*/
    c2 = window->m2d[cl[c]];
    p = cx / window->h1au1[c2];
    q = cy / window->h2au2[c2];
    r = (dzz[cl[c]]) ? cz / dzz[cl[c]] : 0.0;
    
    /* Get the interpolation weights */
    wincon->wgt[c][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
    wincon->wgt[c][1] = (1.0 - r) * (1.0 - p) * q;
    wincon->wgt[c][2] = (1.0 - r) * p * (1.0 - q);
    wincon->wgt[c][3] = (1.0 - r) * p * q;
    wincon->wgt[c][4] = r * (1.0 - p) * (1.0 - q);
    wincon->wgt[c][5] = r * (1.0 - p) * q;
    wincon->wgt[c][6] = r * p * (1.0 - q);
    wincon->wgt[c][7] = r * p * q;
/*
  } else {
    
    if (wincon->osl == 2 || wincon->osl == 4) {
      set_eo(window, wincon, dzz, &cl[c], &cx[c], &cy[c], &cz[c]);
      weights_eo(window, c, cl[c], cx[c], cy[c], cz[c], dzz);
    } else
      weights_oo(window, c, cl[c], cx[c], cy[c], cz[c], dzz);
  }
    */
}

/* END streamline_atc()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays to use semi_lagrange_atc()                         */
/*-------------------------------------------------------------------*/
void prep_semi_lagrange_atc(geometry_t *window)
{
  win_priv_t *wincon = window->wincon;
  int c, cc, ci, zm1;
  double *dzz = wincon->w9;
  double *nu = wincon->nu;
  double *nv = wincon->nv;
  double *nw = wincon->nw;
  double *mask = wincon->d3;
  int *cl = wincon->s2;

  if (wincon->trasc != FFSL) return;

  reset_map_t_all(window);

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing */
  memset(cl, 0, window->sgsiz * sizeof(int));
  set_dzz(window, dzz);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    nu[c] = -nu[ci];
    nv[c] = -nv[ci];
    nw[c] = -nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
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
  /* Set the cell width  mask */
  for (c = 1; c <= window->enonS; c++) {
    mask[c] = 1.0;
  }
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    mask[c] = 0.5;
  }
  /* Set the wet cell mask                                           */
  for (c = 1; c <= window->enon; c++) {
    wincon->c1[c] = 1;
  }
}

/* END prep_semi_lagrange_atc()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Evaluates a tracer using semi-lagrange                            */
/*-------------------------------------------------------------------*/
void semi_lagrange_atc(geometry_t *window, /* Processing window      */
		       window_t *windat,   /* Window data structure  */
		       win_priv_t *wincon, /* Window constants       */
		       double *tr,         /* Tracer array           */
		       double *ntr,        /* Updated tracer         */
		       int c             
		       )
{
  int cs;
  int *cl = wincon->s2;

  /* Get the tracer concentrations at the forward timestep */
  ntr[c] = int_val(window, cl[c], wincon->wgt[c], tr);
}

double semi_lagrange_rtc(geometry_t *window, /* Processing window      */
		       window_t *windat,   /* Window data structure  */
		       win_priv_t *wincon, /* Window constants       */
		       double *tr,         /* Tracer array           */
		       int c             
		       )
{
  int cs;
  int *cl = wincon->s2;
  double val;
  /* Get the tracer concentrations at the forward timestep */
  val = int_val(window, cl[c], wincon->wgt[c], tr);
  /*
  val = int_valo(window, wincon->wgt[c], wincon->lmap[cl[c]], 
		 tr, wincon->osl);
  */
  return(val);
}

/* END semi_lagrange_atc()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-calculates fluxes and velocities for tratio < 1.               */
/*-------------------------------------------------------------------*/
void recalc_vel(geometry_t *window, /* Processing window      */
		window_t *windat,   /* Window data structure  */
		win_priv_t *wincon  /* Window constants       */
		)
{
  int cc, c, cs, cb, c2, zm1;
  double sum, f, dz;

  /*-----------------------------------------------------------------*/
  /* Reset the e1 fluxes and velocity                                */
  memcpy(windat->u1flux3d, windat->u1vm, window->sgsiz * sizeof(double));
  /* Get the cell spacing and new surface coordinate for u1          */
  set_dz_at_u1(window, windat, wincon);
  /* Reset the flux at the surface cell, and zero above the surface  */
  for (cc = 1; cc < window->b2_e1; cc++) {
    c = cs = window->sur_e1[cc];
    cb = window->bot_e1[cc];
    c2 = window->m2d[c];
    sum = 0.0;
    /* Reset the flux above the surface                              */
    while (c != window->zp1[c]) {
      c = window->zp1[c];
      sum += windat->u1flux3d[c];
      windat->u1flux3d[c] = 0.0;
    }
    /* Recalculate the velocity                                      */
    f = (double)(window->s2k[cs] - window->s2k[cb] + 1);
    c = cs;
    while (f && c != window->zm1[cb]) {
      windat->u1flux3d[c] += sum / f;
      windat->u1[c] = (windat->dzu1[c]) ? windat->u1flux3d[c] / (windat->dzu1[c] * window->h2au1[c2]) : 0.0;
      c = window->zm1[c];
    }
    c = cs; 
    /* Adjust thin layers if required                                */
    if (wincon->thin_merge && windat->dzu1[c] < wincon->hmin && c != cb) {
      zm1 = window->zm1[c];
      sum = windat->u1flux3d[c] + windat->u1flux3d[zm1];
      dz = windat->dzu1[c] + windat->dzu1[zm1];
      windat->u1flux3d[c] = windat->u1flux3d[zm1] = 0.5 * sum;
      windat->u1[c] = windat->u1[zm1] = (dz) ? sum / (dz * window->h2au1[c2]) : 0.0;
    }
    /* Set a no-gradient over the surface                           */
    while (c != window->zp1[c]) {
      c = window->zp1[c];
      windat->u1[c] = windat->u1[cs];
    }
  }
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  set_lateral_BC_vel(windat->u1flux3d, window->nbpte1, window->bpte1, window->bine1);

  /*-----------------------------------------------------------------*/
  /* Reset the e2 fluxes and velocity                                */
  memcpy(windat->u2flux3d, windat->u2vm, window->sgsiz * sizeof(double));
  /* Get the cell spacing and new surface coordinate for u2          */
  set_dz_at_u2(window, windat, wincon);
  /* Reset the flux at the surface cell, and zero above the surface  */
  for (cc = 1; cc < window->b2_e2; cc++) {
    c = cs = window->sur_e2[cc];
    cb = window->bot_e2[cc];
    c2 = window->m2d[c];
    sum = 0.0;
    while (c != window->zp1[c]) {
      c = window->zp1[c];
      sum += windat->u2flux3d[c];
      windat->u2flux3d[c] = 0.0;
    }
    f = (double)(window->s2k[cs] - window->s2k[cb] + 1);
    c = cs;
    while (f && c != window->zm1[cb]) {
      windat->u2flux3d[c] += sum / f;
      windat->u2[c] = (windat->dzu2[c]) ? windat->u2flux3d[c] / (windat->dzu2[c] * window->h1au2[c2]) : 0.0;
      c = window->zm1[c];
    }
    c = cs;
    if (wincon->thin_merge && windat->dzu2[c] < wincon->hmin && c != cb) {
      zm1 = window->zm1[c];
      sum = windat->u2flux3d[c] + windat->u2flux3d[zm1];
      dz = windat->dzu2[c] + windat->dzu2[zm1];
      windat->u2flux3d[c] = windat->u2flux3d[zm1] = 0.5 * sum;
      windat->u2[c] = windat->u2[zm1] = (dz) ? sum / (dz * window->h1au2[c2]) : 0.0;
    }
    while (c != window->zp1[c]) {
      c = window->zp1[c];
      windat->u2[c] = windat->u2[cs];
    }
  }
  vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
	    window->bpte2, window->bine2, wincon->slip);
  set_lateral_BC_vel(windat->u2flux3d, window->nbpte2, window->bpte2, window->bine2);
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
  int cc, c, cs, cb, zm1;
  double us;

  /* e1 fluxes                                                       */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->sur_e1[cc];
    cb = window->bot_e1[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    if (c < cb && windat->dzu1[zm1]) {
      /*window->s2k[c] != window->s2k[wincon->i6[cc]]) {*/
      window->sur_e1[cc] = zm1;
      windat->u1vm[zm1] += windat->u1vm[c];
      windat->u1vm[c] = 0.0;
      merge_thin_layers(c, zm1, windat->u1, windat->dzu1);
      us = windat->u1[c];
    }
    /* Set a no-gradient to the surface                              */
    while (c != window->zp1[c]) {
      windat->u1[c] = us;
      c = window->zp1[c];
    }
  }

  /* e2 fluxes                                                       */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->sur_e2[cc];
    cb = window->bot_e2[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    if (c < cb && windat->dzu2[zm1]) {
      /*window->s2k[c] != window->s2k[wincon->i7[cc]]) {*/
      window->sur_e2[cc] = zm1;
      windat->u2vm[zm1] += windat->u2vm[c];
      windat->u2vm[c] = 0.0;
      merge_thin_layers(c, zm1, windat->u2, windat->dzu2);
      us = windat->u2[c];
    }
    /* Set a no-gradient to the surface                              */
    while (c != window->zp1[c]) {
      windat->u2[c] = us;
      c = window->zp1[c];
    }
  }
}

/* END merge_volflux()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the velocities given the volume fluxes                    */
/*-------------------------------------------------------------------*/
void flux_to_velocity(geometry_t *window, /* Processing window      */
		      window_t *windat,   /* Window data structure  */
		      win_priv_t *wincon  /* Window constants       */
		      )
{
  int cc, c, cs, cb, c1, c2, zm1;
  double sum, dz; 

  /*-----------------------------------------------------------------*/
  /* Pre-filter the e1 volume fluxes                                 */
  memcpy(windat->u1flux3d, windat->u1vm, window->sgsiz * sizeof(double));
  vel2D_lbc(windat->u1vm, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->sur_e1[cc];
    zm1 = window->zm1[c];
    while (c != zm1) {
      /*windat->u1flux3d[c] = cvol1w(window, windat->u1vm, c);*/
      /*windat->u1flux3d[c] = shuman(window, windat->u1vm, c);*/
      /*windat->u1flux3d[c] = shapiro_smoothxr(window, windat->u1vm, c);*/
      /*windat->u1flux3d[c] = con_med(window, windat->u1vm, c);*/
      c = zm1;
      zm1 = window->zm1[c];
    }
  }
  memcpy(windat->u1vm, windat->u1flux3d, window->sgsiz * sizeof(double));

  /* Get the cell spacing and new surface coordinate for u1          */
  /*set_dz_at_u1(window, windat, wincon);*/

  /* Recompute the velocities                                        */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = cs = c1 = window->sur_e1[cc];
    cb = window->bot_e1[cc];
    c2 = window->m2d[c];
    zm1 = window->zm1[c];
    while (c != zm1) {
      if (wincon->thin_merge && windat->dzu1[c] < wincon->hmin && c == cs && c != cb) {
	sum = windat->u1vm[c] + windat->u1vm[zm1];
	dz = windat->dzu1[c] + windat->dzu1[zm1];
	/*windat->u1vm[c] = windat->u1vm[zm1] = 0.5 * sum;*/
	windat->u1vm[zm1] = sum;
	windat->u1vm[c] = 0.0;
	windat->u1[c] = windat->u1[zm1] = (dz) ? sum / (dz * window->h2au1[c2]) : 0.0;
      } else {
	dz = windat->dzu1[c];
	windat->u1[c] = (dz) ? windat->u1vm[c] / (dz * window->h2au1[c2]) : 0.0;
      }
      c = zm1;
      zm1 = window->zm1[c];
    }
    while (c1 != window->zp1[c1]) {
      windat->u1[c1] = windat->u1[cs];
      c1 = window->zp1[c1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Pre-filter the e2 volume fluxes                                 */
  memcpy(windat->u2flux3d, windat->u2vm, window->sgsiz * sizeof(double));
  vel2D_lbc(windat->u2vm, window->nbpte2, window->nbe2,
	    window->bpte2, window->bine2, wincon->slip);
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->sur_e2[cc];
    zm1 = window->zm1[c];
    while (c != zm1) {
      /*windat->u2flux3d[c] = cvol1w(window, windat->u2vm, c);*/
      /*windat->u2flux3d[c] = shuman(window, windat->u2vm, c);*/
      /*windat->u2flux3d[c] = shapiro_smoothyr(window, windat->u2vm ,c);*/
      /*windat->u2flux3d[c] = con_med(window, windat->u2vm, c);*/
      c = zm1;
      zm1 = window->zm1[c];
    }
  }
  memcpy(windat->u2vm, windat->u2flux3d, window->sgsiz * sizeof(double));

  /* Get the cell spacing and new surface coordinate for u2          */
  /*set_dz_at_u2(window, windat, wincon);*/

  /* Recompute the velocities                                        */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = cs = c1 = window->sur_e2[cc];
    cb = window->bot_e2[cc];
    c2 = window->m2d[c];
    zm1 = window->zm1[c];
    while (c != zm1) {
      if (wincon->thin_merge && windat->dzu2[c] < wincon->hmin && c == cs && c != cb) {
	sum = windat->u2vm[c] + windat->u2vm[zm1];
	dz = windat->dzu2[c] + windat->dzu2[zm1];
	/*windat->u2vm[c] = windat->u2vm[zm1] = 0.5 * sum;*/
	windat->u2vm[zm1] = sum;
	windat->u2vm[c] = 0.0;
	windat->u2[c] = windat->u2[zm1] = (dz) ? sum / (dz * window->h1au2[c2]) : 0.0;
      } else {
	dz = windat->dzu2[c];
	windat->u2[c] = (dz) ? windat->u2vm[c] / (dz * window->h1au2[c2]) : 0.0;
      }
      c = zm1;
      zm1 = window->zm1[c];
    }
    while (c1 != window->zp1[c1]) {
      windat->u2[c1] = windat->u2[cs];
      c1 = window->zp1[c1];
    }
  }
}

/* END flux_to_velocity()                                            */
/*-------------------------------------------------------------------*/
