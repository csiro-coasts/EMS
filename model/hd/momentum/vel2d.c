/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/vel2d.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vel2d.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* mode2d_step()      : Updates the 2D mode variables                */
/* vel_u1av_update()  : Solves the 2D u1 velocity equation           */
/* vel_u2av_update()  : Solves the 2D u2 velocity equation           */
/* advect_u1_2d()     : Solves the u1 advection equation             */
/* advect_u2_2d()     : Solves the u2 advection equation             */
/* bdry_u1_2d()       : Sets normal and tangential u1 OBC's          */
/* bdry_u2_2d()       : Sets normal and tangential u2 OBC's          */
/* set_map_e1()       : Sets maps self pointing across u1 boundaries */
/* set_map_e2()       : Sets maps self pointing across u1 boundaries */
/* extract_velocity_2d() : Gets u1&u2 velocity from u1&u2 transport  */
/* set_OBC()          : Implements various OBC conditions            */
/* store_old_bdry_values() : Populates the old values data structure */
/* eta_step()         : Calculates free surface height               */
/* bdry_eta()         : Sets the free surface boundary condition     */
/* vel2D_lbc()        : Sets 2D velocity lateral boundary conditions */
/* mdxs()             : Returns total depth at e1 faces at time t    */
/* mdxns()            : Returns total depth at e1 faces at time t+1  */
/* mdys()             : Returns total depth at e2 faces at time t    */
/* mdyns()            : Returns total depth at e2 faces at time t+1  */
/*-------------------------------------------------------------------*/

void precalc_u1_2d(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);
void precalc_u2_2d(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);
void streamfunction(geometry_t *window, window_t *windat,
                    win_priv_t *wincon);
void set_viscosity_2d(geometry_t *window);

/*-------------------------------------------------------------------*/
/* 2D window step part 1                                             */
/*-------------------------------------------------------------------*/
void mode2d_step_window_p1(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon)
{
  double clock;
  
  /*-------------------------------------------------------------*/
  /* Fill 2D variables into the window data structures.          */
  master->win_data_fill_2d(master, window, windat, master->nwindows);

  clock = dp_clock();

  /*-------------------------------------------------------------*/
  /* Set the lateral boundary conditions                         */
#if !GLOB_BC
  vel2D_lbc(windat->u1av, window->nbpte1S, window->nbe1S,
	    window->bpte1S, window->bine1S, wincon->slip);
  vel2D_lbc(windat->u1bot, window->nbpte1S, window->nbe1S,
	    window->bpte1S, window->bine1S, wincon->slip);
  vel2D_lbc(windat->u2av, window->nbpte2S, window->nbe2S,
	    window->bpte2S, window->bine2S, wincon->slip);
  vel2D_lbc(windat->u2bot, window->nbpte2S, window->nbe2S,
	    window->bpte2S, window->bine2S, wincon->slip);
  set_lateral_bc_eta(windat->eta, window->nbptS, window->bpt,
		     window->bin, window->bin2, 1);
  if (!(master->compatible & V1283)) {
    vel2D_lbc(windat->u1avb, window->nbpte1S, window->nbe1S,
	      window->bpte1S, window->bine1S, wincon->slip);
    vel2D_lbc(windat->u2avb, window->nbpte2S, window->nbe2S,
	      window->bpte2S, window->bine2S, wincon->slip);
  }
#endif
  
  /*-------------------------------------------------------------*/
  /* Set auxiliary cells that have been updated to the multi dt  */
  /* buffers.                                                    */
  fill_multidt_aux_2d(window, windat, master->nu1av,
		      windat->u1av_ae);
  fill_multidt_aux_2d(window, windat, master->nu2av,
		      windat->u2av_ae);
  
#if TR_CK
    check_transfers(geom, window, windat, wincon, nwindows, VEL2D);
#endif

  /*-----------------------------------------------------------------*/
  /* Update the 2D velocity                                          */
  vel_u1av_update(window, windat, wincon);
  vel_u2av_update(window, windat, wincon);

  /* Calculate fluxes on fine grid transfer cells for zoomed grids   */
  set_zoom_flux(window, windat, wincon);

  /* Set sources of momentum                                         */
  ss_momentum(window, windat, wincon, VEL2D);

  windat->wclk += (dp_clock() - clock);

  /*-----------------------------------------------------------------*/
  /* Refill the master with updated velocity from the window         */
  /* data structure. This is required at auxiliary cells for the     */
  /* time filtering of u1av & u2av in asselin(). These filtered      */
  /* velocities are then used to get fluxes to update eta, and are   */
  /* required to be defined at front and right edge auxiliary cells. */
  if (master->nwindows > 1)
    master->win_data_empty_2d(master, window, windat, NVELOCITY);
}

/* END mode2d_step_window_p1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2D window step part 2                                             */
/*-------------------------------------------------------------------*/
void mode2d_step_window_p2(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon)
{
  double clock;

  /*-----------------------------------------------------------------*/
  /* Fill the window with updated velocities.                        */
  if (master->nwindows > 1)
    master->win_data_refill_2d(master, window, windat, master->nwindows,
			       NVELOCITY);

  clock = dp_clock();
  
  /*-----------------------------------------------------------------*/
  /* Apply the Asselin filter to remove the computational mode       */
  asselin(window, windat, wincon);

  /* Add the velocity increments to 2D velocity if required          */
  do_vel_relax(window, windat, wincon, VEL2D);
  /*do_vel_increment_2d(window);*/

  /* Calculate the 2D velocity alert diagnostics                     */
  alerts_w(window, VEL2D);

  /* Check for fatal instabilities                                   */
  if (check_unstable(window, windat, wincon, VEL2D)) return;

  /*-----------------------------------------------------------------*/
  /* Get the surface elevation and set the elevation OBC.            */
  windat->t += windat->dtf2 / 2.0;
  eta_step(window, windat, wincon, wincon->ic);
  windat->t += windat->dtf2 / 2.0;
  bdry_eta(window, windat, wincon);
  /* Calculate the elevation alert diagnostics                       */
  alerts_w(window, ETA_A);

  /* Check for fatal instabilities                                   */
  if (check_unstable(window, windat, wincon, ETA_A)) return;

  /*-----------------------------------------------------------------*/
  /* Update the velocities                                           */
  leapfrog_update_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the total depths at the e1 and e2 faces                     */
  /*get_depths(window, windat, wincon);*/

  windat->wclk += (dp_clock() - clock);

  /*-----------------------------------------------------------------*/
  /* Refill the master with elevation and velocity data from the     */
  /* window data structure.                                          */
  if (master->nwindows > 1)
    master->win_data_empty_2d(master, window, windat, VELOCITY);

}

/* END mode2d_step_window_p2()                                       */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* 2D window step part 3                                             */
/*-------------------------------------------------------------------*/
void mode2d_step_window_p3(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon)
{
  /* Set the master time here so that boundary data is read in   */
  /* on the 2D time-step. This is over-written to remove see     */
  /* main mode2d_step loop                                       */
#pragma omp critical
  {
    /* The following two lines must execute atomically */
    master->t    = windat->t;
    master->days = master->t / 86400.0;
  }
  windat->days = master->days;

  if (master->nwindows > 1)
    master->win_data_refill_2d(master, window, windat, master->nwindows, ETA_A);
  {
    double clock = dp_clock();
    get_depths(window, windat, wincon);
    windat->wclk += (dp_clock() - clock);
  }
  if (master->nwindows > 1)
    master->win_data_empty_2d(master, window, windat, DEPTH);

}
/* END mode2d_step_window_p3()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to update the 2D mode                                     */
/*-------------------------------------------------------------------*/
void mode2d_step(geometry_t *geom,    /* Global geometry             */
                 master_t *master,    /* Master data structure       */
                 geometry_t **window, /* Slave geometry              */
                 window_t **windat,   /* Slave data                  */
                 win_priv_t **wincon, /* Slave constants             */
                 int nwindows         /* Num. of windows to process  */
  )
{
  int n;                              /* Counters                    */
  int ic;                             /* 2D mode counter             */
  double oldtime = master->t;
  int mrnk;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  mrnk = master->mpi_rank + 1;

  /* Loop over all windows */
  for( n = 1; n <= nwindows; n++){
    if( (master->dp_mode & DP_MPI) && (n != mrnk ) ){
	    /* If we are using MPI then only process the window
	       associated with our rank. */
	    continue;
    }
    memset(windat[n]->u1flux, 0, window[n]->sgsizS * sizeof(double));
    memset(windat[n]->u2flux, 0, window[n]->sgsizS * sizeof(double));
    memcpy(wincon[n]->oldeta, windat[n]->eta,
      window[n]->sgsizS * sizeof(double));
    wincon[n]->neweta = wincon[n]->d4;
    memset(wincon[n]->neweta, 0, window[n]->sgsizS * sizeof(double));
    if (!(wincon[n]->etarlx & NONE))
      memset(wincon[n]->eta_rlx3d, 0, window[n]->sgsizS * sizeof(double));
    set_viscosity_2d(window[n]);
  }
  /* Zero the 2D fluxes summed over the 3D step on the master if     */
  /* zoomed grids are used. These fluxes are used to adjust          */
  /* interpolated 3D velocities.                                     */
  if (master->dozoom != NOZOOM) {
    memset(master->u1flux, 0, geom->sgsizS * sizeof(double));
    memset(master->u2flux, 0, geom->sgsizS * sizeof(double));
  }

  /*-----------------------------------------------------------------*/
  /* Loop to step the 2D mode iratio times                           */
  for (ic = 0; ic < windat[1]->iratio; ic++) {
    /*---------------------------------------------------------------*/
    /* Do the custom elevation routines on the master                */
    master->t3d = oldtime + (double)(ic+1) * master->dt2d;
    // TIMING_SET;
    if (master->dp_mode & DP_MPI) {
      bdry_eval_eta_w(geom, master, window[mrnk]);
      bdry_eval_u1av_w(geom, master, window[mrnk]);
      bdry_eval_u2av_w(geom, master, window[mrnk]);
    } else {
      bdry_eval_eta_m(geom, master);
      bdry_eval_u1av_m(geom, master);
      bdry_eval_u2av_m(geom, master);
    }
    // TIMING_DUMP(2, "   bdry_eval_2d ");
    if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions using the master          */
#if GLOB_BC
    vel2D_lbc(master->u1av, geom->nbptS, geom->nbe1S,
              geom->bpte1S, geom->bine1S, master->slip);
    vel2D_lbc(master->u1bot, geom->nbptS, geom->nbe1S,
              geom->bpte1S, geom->bine1S, master->slip);
    vel2D_lbc(master->u2av, geom->nbptS, geom->nbe2S,
              geom->bpte2S, geom->bine2S, master->slip);
    vel2D_lbc(master->u2bot, geom->nbptS, geom->nbe2S,
              geom->bpte2S, geom->bine2S, master->slip);
    set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                       geom->bin2, 1);
    /* Backward velocities are required for stress tensors. Ghost    */
    /* cells are not set accurately since a lateral BC is not set on */
    /* nu1av/nu2av before the asselin() filtering; setting BCs here  */
    /* overcomes this.                                               */
    if (!(master->compatible & V1283)) {
      vel2D_lbc(master->u1avb, geom->nbptS, geom->nbe1S,
		geom->bpte1S, geom->bine1S, master->slip);
      vel2D_lbc(master->u2avb, geom->nbptS, geom->nbe2S,
		geom->bpte2S, geom->bine2S, master->slip);
    }
#endif

    /* 
     * Set the ic for each window 
     * Note: For MPI we need only set the ic for this process
     */
    for (n = 1; n <= nwindows; n++)
      wincon[n]->ic = ic;

    /* Do the 1st part of the 2D step : velocity                     */
    dp_vel2d_step_p1();

    /* Set the updated velocity on the master for zoomed grids       */
    if (master->dozoom != NOZOOM)
      global_interp_vel2d(master, geom);

    /* Do the 2nd part of the 2D step : elevation and fluxes         */
    dp_vel2d_step_p2();
    if (master->crf == RS_RESTART) return;

    /* Set eta and 2D fluxes on the master for zoomed grids          */
    if (master->dozoom != NOZOOM)
      global_interp_pre2d(master, geom, ic, windat[1]->dt2d);

#if TR_CK
    check_transfers(geom, window, windat, wincon, nwindows, VEL2D|FLUX);
#endif

    /* Do the 3rd part of the 2D step : depth at e1 and e2 faces     */
    dp_vel2d_step_p3();

    /* Set the face depth on the master for zoomed grids             */
    if (master->dozoom != NOZOOM)
      global_interp_post2d(master, geom);

    /* Couple at the barotropic level for 2-way nesting              */
    /*dump_eta_snapshot(master, window, windat, wincon);*/
  }

  /* This is to synchronise the master time onto 3D time-steps.      */
  /* If we don't do this then over time accumulation errors build up */
  /* and we get out of sync, causing, for one, the mean velocities   */
  /* to go out of whack.                                             */
  master->t = oldtime + master->dtf;
  master->days = master->t / 86400;

  for (n = 1; n <= master->nwindows; n++) {
    windat[n]->t    = master->t;
    windat[n]->days = master->days;
  }
}

/* END mode2d_step()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the lateral boundary conditions for velocity.      */
/*-------------------------------------------------------------------*/
void vel2D_lbc(double *vel,     /* Velocity array                    */
               int sgbpt,       /* Number of boundary cells          */
               int sgn,         /* Number of normal boundary cells   */
               int *bpt,        /* Boundary cell vector              */
               int *bin,        /* Interior cells to bpt             */
               double slip      /* Slip condition                    */
  )
{
  int c, cc, lc;                /* Counters                          */

  /* Set the normal velocity boundary conditions (zero flux)         */
  for (cc = 1; cc <= sgn; cc++) {
    c = bpt[cc];
    vel[c] = 0.0;
  }
  /* Set the tangential velocity boundary conditions (slip)          */
  for (; cc <= sgbpt; cc++) {
    c = bpt[cc];
    lc = bin[cc];
    vel[c] = slip * vel[lc];
  }
}

/* END vel2D_lbc()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the lateral boundary condition for elevation.      */
/*-------------------------------------------------------------------*/
void set_lateral_bc_eta(double *eta,  /* Elevation array             */
                        int sgbpt,    /* Number of boundary cells    */
                        int *bpt,     /* Boundary cell vector        */
                        int *bin,     /* Interior cells to bpt       */
                        int *bin2,    /* Two interior cells to bpt   */
                        int mode      /* mode=1 : zero               */
				      /* mode=0 : zero flux          */ 
  )
{
  int c1, c2, c3, cc;           /* Counters                          */

  if (mode == 1) {
    /* Set the elevation equal to zero at the land cells             */
    for (cc = 1; cc <= sgbpt; cc++) {
      c1 = bpt[cc];
      eta[c1] = 0.0;
    }
  } else if (mode == 0) {
    /* Set a zero flux condition across the land cells               */
    for (cc = 1; cc <= sgbpt; cc++) {
      c1 = bpt[cc];
      c2 = bin[cc];
      eta[c1] = eta[c2];
    }
  } else {
    /* Set a condition so that if the gradient is calculated across  */
    /* land cells than this is equivalent to a one sided gradient.   */
    for (cc = 1; cc <= sgbpt; cc++) {
      c1 = bpt[cc];
      c2 = bin[cc];
      c3 = bin2[cc];
      eta[c1] = 2.0 * eta[c2] - eta[c3];
    }
  }
}

/* END set_lateral_bc_eta()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void vel_u1av_update(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc;                    /* Sparse coodinate / counter */
  int xp1, xm1;                 /* Sparse locations at i+1, i-1 */
  int yp1, ym1;                 /* Sparse locations at j+1, j-1 */
  int xmyp1;                    /* Sparse location at (i-1,y+1) */
  double pgt = 0.0;             /* Pressure gradient term */
  double cot = 0.0;             /* Coriolis term */
  double bft;                   /* Bottom friction term */
  double rst = 0.0;             /* Radiation stress term */
  double *tzp;                  /* Surface height array */
  double *depth;                /* Depth of the water column */
  double midx;                  /* Water depth at the cell face */
  double val;                   /* Dummy */
  double u2au1;                 /* u2 value at the e1 face */
  double botu1;                 /* Bottom velocity at e1 face */
  double botu2;                 /* Bottom velocity at e2 face */
  double Cdu1;                  /* Bottom drag coefficient at the e1 face */
  double rho0 = 1024.0;         /* Reference density */
  int dbc;                      /* Debugging coordinate */

  dbc = (wincon->dbc) ? window->m2d[wincon->dbc] : 0;
  depth = windat->depth_e1;

  /* tzp is a pointer which points to either topz[][] or eta[][].  */
  /* topz[][] only gets updated once every 3d step. For the non- */
  /* linear case, we need to use eta instead, which gets updated */
  /* every 2d step. This value is only used to calculate transport */
  /* (nothing to do with the surface slope term, which always uses */
  /* eta).  */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation */
  precalc_u1_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across u2 OBC's */
  set_map_e1(window);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection */
  if (!(wincon->u1av_f & ADVECT)) {
    if(wincon->momsc2d & ANGULAR) {
      if(!(wincon->ic % 2))
	advect_u1_2d_ang_flux_b(window, windat, wincon);
      else
	advect_u1_2d_ang_flux_f(window, windat, wincon);
    } else if (wincon->momsc2d & LAGRANGE)
      semi_lagrange_u1av(window, windat, wincon);
    else
      advect_u1_2d(window, windat, wincon);
  }

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included) */
  wincon->hor_mix->setup_av(window, windat, wincon, tzp);
  wincon->hor_mix->u1av(window, windat, wincon, tzp);

  /*-----------------------------------------------------------------*/
  /* Add the body forces and step forward in time */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    xmyp1 = window->xmyp1[c];

    midx = wincon->mdx[c];

    /*---------------------------------------------------------------*/
    /* Pressure gradient term. Note : differences bracketed to */
    /* preserve accuracy.  */
    if (!(wincon->u1av_f & PRESS_BT))
      pgt = -((windat->eta[c] - windat->eta[xm1]) * 
	      wincon->topdensu1[c] +
	      (windat->patm[c] - windat->patm[xm1])) / 
	wincon->densavu1[c];

    /*---------------------------------------------------------------*/
    /* Coriolis */
    u2au1 = wincon->w7[c];      /* Set in precalc_u1_2d() */
    if (!(wincon->u1av_f & CORIOLIS))
      cot = wincon->u1c5[c] * u2au1;

    /*---------------------------------------------------------------*/
    /* Bottom friction */
    u2au1 = 0.25 * (windat->u2avb[c] + windat->u2avb[xm1] +
                    windat->u2avb[yp1] + windat->u2avb[xmyp1]);
    botu2 = 0.25 * (windat->u2bot[c] + windat->u2bot[xm1] +
                    windat->u2bot[yp1] + windat->u2bot[xmyp1]);
    botu1 = windat->u1avb[c] + windat->u1bot[c];
    botu2 = u2au1 + botu2;

    Cdu1 = 0.5 * (wincon->Cd[xm1] + wincon->Cd[c]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu1 * max(wincon->uf, val);
    /* Truncate to ensure stability */
    if (val > depth[c] * midx / windat->dt2d)
      val = depth[c] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case */
    bft = -val * botu1 / depth[c];
    if (c == dbc) {
      wincon->b1 = pgt;
      wincon->b2 = cot;
      wincon->b3 = bft;
    }

    /*---------------------------------------------------------------*/
    /* Calculate nu1av value */
    /* SIGMA : multiply new velocity by depth (midx) */
    windat->nu1av[c] +=
      windat->dt2d * (midx * (pgt + cot) + bft + wincon->u1inter[c]);
  }
  debug_c(window, D_UA, D_POST);

  /*-----------------------------------------------------------------*/
  /* Add the non-linear terms */
  /* SIGMA : No extra terms to include for the sigma case.  */
  if (wincon->nonlinear && !wincon->sigma) {
    for (cc = 1; cc <= window->v2_e1; cc++) {
      c = window->w2_e1[cc];
      xm1 = window->xm1[c];
      windat->nu1av[c] -= windat->dt2d * windat->u1av[c] *
        (windat->detadt[xm1] + windat->detadt[c]) /
        (2.0 * max(depth[c], wincon->hmin));
    }
  }
  
  /*-----------------------------------------------------------------*/
  /* Add the radiation stresses                                      */
  if (wincon->waves & WAVE_FOR) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s3[cc];
      xm1 = window->xm1[c];
      rst = windat->dt2d * windat->wave_Fx[c] /
	                 (rho0 * max(depth[c], wincon->hmin));
      if (windat->u1_rad) windat->u1_rad[c] = rst;
      windat->nu1av[c] += rst;
    }
  }
  else if (wincon->waves & TAN_RAD) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s3[cc];
      rst = windat->dt2d * (windat->wave_Sxy[window->yp1[c]] - 
	                    windat->wave_Sxy[window->ym1[c]]) / 
	    (window->h2au2[c] + window->h2au2[window->yp1[c]]);
      if (windat->u1_rad) windat->u1_rad[c] = rst;
      windat->nu1av[c] += rst;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Restrict flow if column nearly dry */
  for (cc = 1; cc <= window->v2_e1; cc++) {
    c = window->w2_e1[cc];
    xm1 = window->xm1[c];
    if (windat->nu1av[c] > 0.0 &&
        (val = (tzp[xm1] - window->botz[xm1])) < wincon->hmin)
      windat->nu1av[c] *= max(val, 0.0) / wincon->hmin;
    else if (windat->nu1av[c] < 0.0 &&
             (val = (tzp[c] - window->botz[c])) < wincon->hmin)
      windat->nu1av[c] *= max(val, 0.0) / wincon->hmin;
  }

  /*-----------------------------------------------------------------*/
  /* Blend velocities over a subregion if required                   */
  blend_vel(window, windat, wincon, U12D, windat->nu1av);

  /*-----------------------------------------------------------------*/
  /* Calculate u1av boundary values */
  bdry_u1_2d(window, windat, wincon);
  debug_c(window, D_UA, D_BDRY);
}

/* END vel_u1av_update()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/*-------------------------------------------------------------------*/
void advect_u1_2d(geometry_t *window, /* Processing window           */
                  window_t *windat,   /* Window data structure       */
                  win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc;                /* Sparse coordinate / counter           */
  int xp1, xm1;             /* Sparse coordinate at i+1, i-1         */
  int yp1, ym1;             /* Sparse coordinate at j+1, j-1         */
  int xmyp1;                /* Sparse coordinate at (i-1,j+1)        */
  int *ctp;                 /* Cells to process vector               */
  int vcs;                  /* Number of cells to process            */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double u1val;             /* Cell centered u1 value                */
  double u2val;             /* Cell cornered u2 value                */
  double h1;                /* Face centered grid spacing            */
  double h1av, h2av;        /* Mean cell spacings x and y directions */
  double m, d1;             /* Dummy variables                       */
  double *Du1;              /* Dummy for e1 advective fluxes         */
  double *vel;              /* Velocity to use in spatial gradients  */
  double *area;             /* Cell area                             */
  double *u1sh1h2_av;       /* u1 * u1 * h1 * h2                     */
  double *u1u2h1s_av;       /* u1 * u2 * h1 * h1                     */
  double *depth;            /* Depth of the water column             */
  double *dvel;             /* Velocity to interpolate               */
  double *vface;            /* Velocity centered on the grid face    */
  double *cn;               /* Courant number                        */
  double midx;              /* Depth at the e1 face                  */
  double *u2au1;            /* u2 value at the e1 face               */
  double lr;                /* Ratio of dtf : dt2d                   */

  if (!wincon->nonlinear)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  Du1 = wincon->d1;
  u1sh1h2_av = wincon->d2;
  u1u2h1s_av = wincon->d3;
  area = window->cellarea;
  depth = windat->depth_e1;
  vel = windat->u1av;
  u2au1 = wincon->w7;
  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  if(windat->dtu1 == windat->dtf)
    dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping  */
  else
    dtu *= 2.0;           /* 2D leapfrog time-step with sub-stepping */
  lr = windat->dtf2 / windat->dt2d;
  dtm = dtu * lr;
  trem = windat->dt2d;
  tremf = windat->dtf2;
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Initialise                                                  */
      /* memset(depth,0,window->sgsizS*sizeof(double)); */
      memset(Du1, 0, window->sgsizS * sizeof(double));
      memset(u1sh1h2_av, 0, window->sgsizS * sizeof(double));
      memset(u1u2h1s_av, 0, window->sgsizS * sizeof(double));

      /*-------------------------------------------------------------*/
      /* Get the tracer values in multi-dt auxiliary cells. This is  */
      /* a linear interpolation between the value at the start of    */
      /* the timestep and the end of longer timesteps.               */
      set_multidt_tr_2d(window, windat, trem, vel, windat->u1av_as,
                        windat->u1av_ae);

      /*-------------------------------------------------------------*/
      /* Loop over the u1 points to calculate depths. Note: Du1 at   */
      /* the coordinate ym1 is required for the cell corner fluxes,  */
      /* and this cell is not included in the e1 cells to process    */
      /* vector, hence Du1 must be calculated over auxiliary cells   */
      /* on all sides in order to accomodate this.                   */
      for (cc = 1; cc <= window->x2_e1; cc++) {
        c = window->w2_e1[cc];
        Du1[c] = depth[c] * vel[c];
      }
      if (wincon->momsc2d & VANLEER) {
	for (cc = window->b2_e1 + 1; cc <= window->n2_e1; cc++) {
	  c = window->w2_e1[cc];
	  xm1 = window->xm1[c];
	  Du1[xm1] = depth[xm1] * vel[xm1];
	}
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres. This must be        */
      /* calculated at ghost cells also since these locations will   */
      /* give a non-zero cell centered flux.                         */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xp1 = window->xp1[c];
          u1val = 0.5 * (vel[xp1] + vel[c]);

          /* SIGMA : Multiply by the cell centre depth               */
          if (u1val > 0)
            m = Du1[c] * wincon->mdx[c];
          else
            m = Du1[xp1] * wincon->mdx[xp1];
          u1sh1h2_av[c] = m * u1val * area[c];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xp1 = window->xp1[c];
          u1val = 0.5 * (vel[xp1] + vel[c]);

          /* SIGMA : Multiply by the cell centre depth               */
          u1sh1h2_av[c] =
            u1val * area[c] * 0.5 * (Du1[c] * wincon->mdx[c] +
                                     Du1[xp1] * wincon->mdx[xp1]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xp1 = window->xp1[c];
          vface[c] = 0.5 * (vel[xp1] + vel[c]);
          cn[c] = vface[c] * dtu / window->h1acell[c];
          /* The interpolated velocity must be shifted by xp1 when   */
          /* using van_leer_do() since this routine is set up to     */
          /* interpolate onto the cell face and we need to           */
          /* interpolate to the cell center.                         */
          dvel[c] = Du1[xp1] * wincon->mdx[xp1];
        }
	/* Ghost cells 2 cells in */
        for (cc = window->b2_e1 + 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xm1 = window->xm1[c];
          dvel[xm1] = Du1[c] * wincon->mdx[c];
	}
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array u2h1h2.                                       */
        van_leer_do(u1sh1h2_av, dvel, vface, cn, 1, window->n2_e1,
                    window->w2_e1, window->xp1, window->xm1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          u1sh1h2_av[c] *= vface[c] * area[c];
          /* 
             xp1=window->xp1[c]; u1val = 0.5*(vel[xp1]+vel[c]);
             u1sh1h2_av[c]=u1val*area[c]*0.5*(Du1[c]*wincon->mdx[c]+
             Du1[xp1]*wincon->mdx[xp1]); */
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners. This must be        */
      /* calculated at ghost cells also since these locations will   */
      /* give a non-zero flux at ouside corners.                     */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          h1 = 0.5 * (window->h1au2[xm1] + window->h1au2[c]);
          u2val = 0.5 * (windat->u2av[xm1] * wincon->mdy[xm1] +
                         windat->u2av[c] * wincon->mdy[c]);
          if (u2val > 0)
            m = Du1[ym1];
          else
            m = Du1[c];
          u1u2h1s_av[c] = m * h1 * h1 * u2val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          h1 = 0.5 * (window->h1au2[xm1] + window->h1au2[c]);
          u2val = 0.5 * (windat->u2av[xm1] * wincon->mdy[xm1] +
                         windat->u2av[c] * wincon->mdy[c]);
          h1av =
            window->h1au2[xm1] / (window->h1au2[xm1] + window->h1au2[c]);
          h2av =
            window->h2au1[ym1] / (window->h2au1[ym1] + window->h2au1[c]);
          /* SIGMA : Multiply by the cell corner depth               */
          u2val = h1av * (windat->u2av[c] * wincon->Ds[c] -
                          windat->u2av[xm1] * wincon->Ds[xm1]) +
            windat->u2av[xm1] * wincon->Ds[xm1];
          m = h2av * (Du1[c] - Du1[ym1]) + Du1[ym1];

          u1u2h1s_av[c] = m * h1 * h1 * u2val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          xm1 = window->xm1[c];
          vface[c] = 0.5 * (windat->u2av[xm1] + windat->u2av[c]);
          cn[c] =
            2.0 * vface[c] * dtu / (window->h2au2[c] + window->h2au2[xm1]);
          dvel[c] = Du1[c] * wincon->mdx[c];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array dvel.                                         */
        van_leer_do(u1u2h1s_av, dvel, vface, cn, 1, window->n2_e1,
                    window->w2_e1, window->yp1, window->ym1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n2_e1; cc++) {
          c = window->w2_e1[cc];
          u1u2h1s_av[c] *= vface[c] * area[c];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the u1 updated solution                                 */
      for (cc = 1; cc <= vcs; cc++) {
        c = ctp[cc];
        xm1 = window->xm1[c];
        yp1 = window->yp1[c];
        xmyp1 = window->xmyp1[c];
        midx = wincon->mdx[c];

        /* SIGMA : Multiply by cell depth                            */
        d1 = wincon->u1c1[c] *
          (u1sh1h2_av[c] - u1sh1h2_av[xm1] + u1u2h1s_av[yp1] -
           u1u2h1s_av[c]) / max(depth[c], wincon->hmin) +
          wincon->u1c3[c] * vel[c] * vel[c] * midx +
          wincon->u1c4[c] * u2au1[c] * u2au1[c] * midx;

        windat->nu1av[c] += d1 * dtu;

        /* Integrate the dispersion terms over the sub-timestep.     */
        /* u1adv is divided by dt to get the mean over the 2D mode   */
        /* in vel_u1_update().                                       */
        wincon->u1adv[c] -= (d1 * dtm);
      }
    }

    /* Get the timestep for the next iteration.                      */
    /* Decrement the time remaining. This must be computed for       */
    /* both the leapfrog step (trem) used to update velocity, amd    */
    /* the forward step (tremf) used to integrate u1inter[] over     */
    /* forward time-step windat->dtf.                                */
    trem -= dtu;
    tremf -= dtm;
    dtu = dtm / lr;
    if (trem < dtu)
      dtu = trem;
    if (tremf < dtm)
      dtm = tremf;
    if(!wincon->sigma)
      vel = windat->nu1av;

    if (wincon->nbl1) {
      vel2D_lbc(vel, window->nbpte1S, window->nbe1S,
		window->bpte1S, window->bine1S, wincon->slip);
      blend_vel(window, windat, wincon, U12D, vel);
    }

  }
  debug_c(window, D_UA, D_ADVECT);
}

/* END advect_u1_2d()                                          (slf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the correct map across open boundaries for u1      */
/*-------------------------------------------------------------------*/
void set_map_e1(geometry_t *window  /* Window data structure */
  )
{
  int n;
  int cc, cs, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries */
  /* (set in the tracer routine) to be self-mapping on the open */
  /* boundary.  */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodex & R_EDGE && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e1; cc++) {
        c = open[n]->obc_e1[cc];
        cs = window->xm1[c];
        window->xp1[cs] = c;
      }
    }

    /* Set maps to be self-mapping interior to northern open */
    /* boundaries.  */
    if (open[n]->ocodey & F_EDGE && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e2; cc++) {
        c = open[n]->obc_e2[cc];
        cs = window->ym1[c];
        window->yp1[cs] = cs;
      }
    }
  }
}

/* END set_map_e1()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the u1 velocities on all open boundaries in a  */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_u1_2d(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int n;                        /* Counters */
  int c, cc, c1;                /* Sparse coordinates / counters */
  double *depth;                /* Depth of the water column */
  double midx, midx1;           /* Depth at the e1 face */
  open_bdrys_t **open = window->open;

  depth = windat->depth_e1;

  /*-----------------------------------------------------------------*/
  /* Set the u1av velocity where u1av is normal to the boundary */
  for (n = 0; n < window->nobc; n++) {
    /* Get the u1av velocity derived from a no gradient condition */
    /* on u1flux.  */
    if (open[n]->bcond_nor2d & NOGRAD) {
      /*
      obc_nograd_2d(windat->nu1av[c], depth, wincon->mdx[c], 
		    window->h2au1[c1], wincon->hmin, 1, open[n]->no2_e1, 
		    open[n]->obc_e1, open[n]->oi1_e1[cc]);
      */
      for (cc = 1; cc <= open[n]->no2_e1; cc++) {
        c = open[n]->obc_e1[cc];
        c1 = open[n]->oi1_e1[cc];
        midx = wincon->mdx[c];
        midx1 = wincon->mdx[c1];
        /* Set the no gradient condition on velocity */
        windat->nu1av[c] =
          windat->nu1av[c1] * depth[c1] * midx1 * window->h2au1[c1] /
          (max(depth[c], wincon->hmin) * midx * window->h2au1[c]);
      }
    }
    /* Set the u1av normal velocity for all other BC's.  */
    else {
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e1,
              open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
              open[n]->cyc_e1, windat->nu1av, windat->u1av, 
	      windat->u1avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau1av, 
	      open[n]->transfer_u1av, open[n]->relax_zone_nor, U1BDRY);    
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the u1av velocity where u1av is tangential to the boundary */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	    open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	    open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	    windat->u1av, windat->u1avb, open[n]->bcond_tan2d, 
	    windat->dtf2, &open[n]->datau1av, open[n]->transfer_u1av, 
	    open[n]->relax_zone_tan, U1BDRY);	    
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      int *mi = (open[n]->ocodex & L_EDGE) ? window->xm1 : window->xp1;
      for (cc = 1; cc <= open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	c1 = open[n]->obc_t[cc];
	windat->nu1av[c1] = windat->nu1av[c];
      }
    }
  }

  /* Rescale for sigma active OBCs */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->type & U1BDRY && open->bcond_nor2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1av, mdxns,
			1, open->no2_e1, open->obc_e1);
      if (open->type & U2BDRY && open->bcond_tan2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1av, mdxns,
			open->no2_e1 + 1, open->to2_e1, open->obc_e1);
    }
  }
}

/* END bdry_u1_2d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the boundary conditions on the tracers                        */
/*-------------------------------------------------------------------*/
void set_OBC(geometry_t *window,   /* Processing window              */
             window_t *windat,     /* Window data structure          */
             win_priv_t *wincon,   /* Window geometry / constants    */
             open_bdrys_t *open,   /* OBC structure                  */
             int sb,               /* Start of OBC vector            */
             int eb,               /* End of OBC vector              */
             int *obc,             /* Open boundary vector           */
             int *oi1,             /* One cell interior to obc cells */
             int *oi2,             /* Two cell interior to obc cells */
             int *cyc,             /* Cyclic sparse coordinate       */
             double *vel,          /* Variable array at time t+1     */
             double *vel_t,        /* Variable array at time t       */
             double *vel_b,        /* Variable array at time t-1     */
             int bcond,            /* Open boundary condition        */
             double dt,            /* Time step                      */
             bdry_details_t *data, /* Custom data structure          */
	     double *transfer,     /* Transfer data                  */
	     int rlxn,             /* Flow relaxation zone           */
	     int code              /* General purpose code           */
  )
{
  int c, cc;                       /* Sparse coorindate / counter    */
  int c1, c2, c3, c4;              /* Sparse coordinates             */
  double *newval;                  /* New velocity values            */
  double *fval;                    /* Forced velocity values         */
  int *imap = NULL;                /* Interior cell map              */
  double *cs;                      /* Phase speed                    */
  double sf;                       /* Phase speed smoothing factor   */
  double rts;                      /* Relaxation time scale          */

  /* Return if no valid cells on this boundary                       */
  if (eb < sb) return;

  /* Initialise the pointers                                         */
  newval = wincon->w5;
  fval = wincon->w6;
  memset(fval, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Set the pointer to phase speed to allow phase speed smoothing   */
  /* for elevation if required.                                      */
  if (open->spf && vel_t == windat->eta) {
    cs = open->sphase;
    sf = open->spf;
  } else {
    cs = open->dum;
    memset(cs, 0, (open->no3_t + 1) * sizeof(double));
    sf = 0.0;
  }

  /*-----------------------------------------------------------------*/
  /* No action taken                                                 */
  /*if (bcond & (NOTHIN|LINEAR))*/
  if (bcond == NOTHIN || bcond == (NOTHIN|TIDALH) || bcond == LINEAR ||
      bcond == (NOTHIN|TIDALC))
    return;
  if (open->adjust_flux && code == ETASPEC) {
    if (bcond == (NOTHIN|FILEIN) || bcond == (NOTHIN|TIDALH|FILEIN) || 
	bcond == (NOTHIN|TIDALC|FILEIN) || bcond == (NOTHIN|TIDEBC|FILEIN))
      return;
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient condition                                     */
  else if (bcond & NOGRAD) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      newval[c] = vel[c1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Cyclic                                                          */
  else if (bcond & (CYCLIC | CYCLED)) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = cyc[cc];      
      newval[c] = vel[c1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Clamped                                                         */
  else if (bcond & CLAMPD) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      newval[c] = 0.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Gravity wave radiation                                          */
  else if (bcond & GRAVTY) {
    double *hat = NULL;
    if (open->ocodex & (L_EDGE | R_EDGE))
      hat = window->h1acell;
    else if (open->ocodey & (B_EDGE | F_EDGE))
      hat = window->h2acell;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = window->m2d[c1];
      newval[c] = bc_gravity(window, windat, wincon, hat[c2], dt,
			     c2, vel_t[c], vel[c1], &cs[cc], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set an Orlanski radiation boundary condition                    */
  else if (bcond & ORLANS) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      newval[c] = bc_orlanski(vel_t[c], vel_t[c1], vel_t[c2], vel[c1],
			      vel_b[c], vel_b[c1], &cs[cc], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Camerlengo & O'Brien radiation boundary condition         */
  else if (bcond & CAMOBR) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      newval[c] = bc_camerlengo(vel_t[c], vel_t[c1], vel_t[c2], vel[c1],
				vel_b[c], vel_b[c1], &cs[cc], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Miller and Thorpe radiation boundary condition            */
  else if (bcond & MILLER) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      newval[c] = bc_miller(vel_t[c], vel_t[c1], vel_t[c2], vel[c1],
                            vel_b[c], vel_b[c1], vel_b[c2],
			    &cs[cc], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Raymond and Kuo radiation condition                       */
  else if (bcond & RAYMND) {
    int tp, tm, tp1, tm1;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      tp = open->tmpp[c];
      tm = open->tmpm[c];
      tp1 = open->tmpp[c1];
      tm1 = open->tmpm[c1];
      newval[c] = bc_raymond(vel_t[c], vel_t[c1], vel_t[c2], vel[c1],
                             vel[c2], vel_t[tp], vel_t[tm],
                             vel_t[tp1], vel_t[tm1], &cs[cc], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Least squares linear interpolation                              */
  if (bcond & LINEXT) {
    imap = open->nmap;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      c3 = imap[c2];
      c4 = imap[c3];
      newval[c] = bc_leastsq(vel[c4], vel[c3], vel[c2], vel[c1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Polynomial extrapolation (2nd order)                            */
  if (bcond & POLEXT) {
    imap = open->nmap;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = oi2[cc];
      c3 = imap[c2];
      newval[c] = bc_polint(vel[c3], vel[c2], vel[c1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Vertically  averaged velocity  from the vertically  integrated  */
  /* velocity. Note, velocities  above the  free surface  are set to */
  /* zero.                                                           */
  if (bcond & VERTIN) {
    double *sum = wincon->d1;
    int zm1;
    double *vel = NULL;
    double *dz = NULL;
    double *bottom = NULL;
    double *depth = NULL;
    double *md = NULL;
    int *m1 = NULL;

    /*UR-TODO  - nu1 and nu2 seem to be un-initialised under pthread*/
    /* Set the pointers and initialise                               */
    memset(sum, 0, window->sgsizS * sizeof(double));
    if (code == U1BDRY) {
     /* vel = wincon->w9; */
      vel = windat->nu1; 
      dz = windat->dzu1;
      bottom = window->botzu1;
      depth = windat->depth_e1;
      md = wincon->mdx;
      m1 = window->xm1;
    } else if (code == U2BDRY) {
/*      vel = wincon->w10; */
      vel = windat->nu2; 
      dz = windat->dzu2; 
      bottom = window->botzu2;
      depth = windat->depth_e2;
      md = wincon->mdy;
      m1 = window->ym1;
    }

    /* Do the vertical integration                                   */
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = window->m2d[c];
      zm1 = window->zm1[c];
      while(c != zm1) {
	sum[c1] += vel[c] * dz[c] * md[c1];
	c = zm1;
	zm1 = window->zm1[c];
      }
    }

    /* Get the  vertically integrated  velocity. Note: this velocity */
    /* is invariant over the 2D time-step, hence the depth used must */
    /* be  the  3D  depth calculated  from  eta at the start  of the */
    /* time-step (i.e. oldeta  =  topz). The  2D  depth (depth_e1 or */
    /* depth_e2) must  also  be set to this depth so that this depth */
    /* is used in  eta_step()  to  calculate the flux, which is used */
    /* in velocity_adjust() to adjust the 3D velocities. */
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      depth[c] = (max(windat->topz[m1[c]], windat->topz[c]) - 
		  bottom[c]);/* * md[window->m2d[c]];*/
      newval[c] = sum[c] / max(depth[c] * md[window->m2d[c]], wincon->hmin);
    }
    /* Set the depth at OBC auxiliary cells */
    for (cc = 1; cc <= open->no2_a; cc++) {
      c = open->obc_a[cc];
      depth[c] = (max(windat->topz[m1[c]], windat->topz[c]) - 
		  bottom[c]);
    }
  }


  /*-----------------------------------------------------------------*/
  /* Local solution for tangential depth averaged velocity           */
  if (bcond & LOCALT) {
    if (open->type == U1BDRY)
      u2av_local(window, windat, wincon, open, newval, sb, eb, 1);
    if (open->type == U2BDRY) {
      u1av_local(window, windat, wincon, open, newval, sb, eb, 1);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Local solution for normal depth averaged velocity               */
  if (bcond & LOCALN) {
    if (open->type == U1BDRY)
      u1av_local(window, windat, wincon, open, fval, sb, eb, 0);
    if (open->type == U2BDRY)
      u2av_local(window, windat, wincon, open, fval, sb, eb, 0);
  }

  /*-----------------------------------------------------------------*/
  /* Local solution for elevation                                    */
  if (bcond & LOCALE) {
    eta_local(window, windat, wincon, open, fval, sb, eb);
  }

  /*-----------------------------------------------------------------*/
  /* Tidal synthesis from harmonics                                  */
  if (bcond & (TIDALH|TIDALC)) {
    double gmt = windat->t - wincon->tz;
    double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
    if (code == ETASPEC) {
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	/* Note; argument to csr_tide_eval() is GMT; correct for       */
	/* timezone.                                                   */
	fval[c] += (ramp * csr_tide_eval(&open->tc, cc, gmt));
      }
    } else if (code & U1BDRY) {
      double fvx, fvy;
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	c2 = window->m2d[c];
	fvx = (ramp * csr_tide_eval(&open->t1u, cc, gmt));
	fvy = (ramp * csr_tide_eval(&open->t1v, cc, gmt));
	fval[c] += (fvx * window->costhu1[c2] + fvy * window->sinthu1[c2]);
      }
    } else {
      double fvx, fvy;
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	c2 = window->m2d[c];
	fvx = (ramp * csr_tide_eval(&open->t2u, cc, gmt));
	fvy = (ramp * csr_tide_eval(&open->t2v, cc, gmt));
	fval[c] += (fvx * window->costhu2[c2] + fvy * window->sinthu2[c2]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Custom specification. Note CUSTOM is handled by the master      */
  if (bcond & CUSTOM) {
    double ramp = (wincon->rampf & CUSTOM) ? windat->rampval : 1.0;

    /* Blend geostrophic velocities with forcing over the ramp if    */
    /* required (i.e. if VELOCITY GEOSTROPHIC is ued with CUSTOM u1  */
    /* or u2 forcing).                                               */
    if (wincon->vinit & VINIT_GEO && ramp < 1.0 && (code & (U1GEN|U2GEN))) {
      double *bvel = (code & U1GEN) ? open->u1d : open->u2d;
      double fvel;
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	c2 = window->m2d[c];
	fvel = bdry_value_w(window, windat, wincon, open, data,
			    c, cc, windat->t);
	fval[c] = (1.0 - ramp) * bvel[cc] + ramp * fvel;
      }
    } else {
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	fval[c] += ramp * bdry_value_w(window, windat, 
				       wincon, open, data,
				       c, cc, windat->t);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* File input. Note FILEN is handled by the master, including ramp */
  if (bcond & FILEIN) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      fval[c] += transfer[cc];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Geostropic 3D current                                           */
  if (bcond & CUSTOM && open->options & (OP_GEOSTR|OP_YANKOVSKY)) {
    int m1, ci;
    double *dz, *md, *cf, *hat, *depth;
    double d1, d2, d3, d4, fa, fm, fb, d, sgn;
    double *sum = wincon->d1;
    sgn = 1.0;
    if (open->type == U1BDRY) {
      dz = windat->dzu1;
      md = wincon->mdx;
      cf = wincon->u1c6;
      hat = window->h2au1;
      depth = windat->depth_e1;
      if (open->ocodex & R_EDGE) sgn = -1.0;
    } else {
      dz = windat->dzu2;
      md = wincon->mdy;
      cf = wincon->u2c6;
      hat = window->h1au2;
      depth = windat->depth_e2;
      if (open->ocodey & F_EDGE) sgn = -1.0;
    }
    /* Integrate the flow over the face (barotropic component) */
    d1 = 0.0;
    d = 0.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      d1 += fval[c] * dz[c] * hat[window->m2d[c]];
      d += dz[c];
      vel[c] = fval[c];
      open->flow[cc] = fval[c];
    }
    /* Get the baroclinic component */
    memset(sum, 0, window->sgsizS * sizeof(double));
    d2 = d4 = 0.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      ci = open->obc_t[cc];
      c2 = window->m2d[ci];
      if (window->gridz[ci] > windat->eta[c2]) continue;
      m1 = (open->ocodex & L_EDGE || open->ocodey & B_EDGE) ? open->ogc_t[cc] : open->nmap[c];
      sum[c2] += wincon->g * md[c2] * (windat->dens[c] - windat->dens[m1]) * dz[ci];
      d3 = windat->dt * cf[c2] * sum[c2] / (0.5 * (windat->dens[c] + windat->dens[m1]));
      fval[c] += d3;
      newval[c] = d3;
      /* Yankovsky (2000) OBC */
      if (open->options & OP_YANKOVSKY)
	fval[c] = vel_t[oi1[cc]];
      d2 += fval[c] * dz[c] * hat[window->m2d[c]];
      d4 += d3 * dz[c] * hat[window->m2d[c]];
    }
    /* Adjust the profile so that the total flux = barotropic component */
    d3 = 0.0;
    c2 = window->m2d[obc[sb]];
    fa = (d1 - d2) / (depth[c2] * hat[c2]);
    fm = (d2) ? fabs(d1/d2) : 1.0;
    fb = (d1) ? (d1 - d4) / d1 : 0.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      ci = open->obc_t[cc];
      c2 = window->m2d[ci];
      if (window->gridz[ci] > windat->eta[c2]) continue;
      if (open->options & (OP_YANKOVSKY|OP_MULTF)) {
	/* Multiplicative scaling */
	fval[c] = (d2 * sgn > 0.0) ? fval[c] * fm : fval[c] + fa;
	/*fval[c] = (d1) ? fb * vel[c] + newval[c] : newval[c];*/
      } else {
	/* Additive scaling */
	fval[c] += fa;
      }
      d3 += fval[c] * dz[c] * hat[window->m2d[c]];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Flather radiation (elevation components).                       */
  /* NOTHIN for eta                                                  */
  if (bcond == (FLATHE|FILEIN|NOTHIN))
    return;
  /* FLATHR + Local Solution for eta forcing                         */
  if (bcond == (FLATHE|LOCALE)) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      vel[c] = fval[cc];
    }
    return;
  }
  /* Use FILEIN on eta                                               */
  if (bcond == (FLATHE|FILEIN)) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      vel[c] = transfer[cc];
    }
    return;
  }
  /* Use FILEIN + tide on eta                                        */
  if (bcond == (FLATHE|FILEIN|TIDALH)) {
    double gmt = windat->t - wincon->tz;
    double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      vel[c] = transfer[cc] + (ramp * csr_tide_eval(&open->tc, cc, gmt));
    }
    return;
  }
  /* The (FILEIN|CUSTOM) component is used for the 2D velocity OBC,  */
  /* the radiation component is used for the eta OBC.                */
  if (bcond & FLATHE) {
    if (!open->relax_time)open->relax_time = open->relax_timei = 1.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      fval[c] = vel_t[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Tidal constituent specification                                 */
  if (bcond & TIDEBC) {
    double depth, spd;
    int c1, c2;
    double ramp = (wincon->rampf & TIDEBC) ? windat->rampval : 1.0;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      fval[c] += (ramp * bc_tidal(window, windat, wincon, open->tideforce, c));
    }
  }

  /*-----------------------------------------------------------------*/
  /* LINEAR computation or NOTHIN.                                   */
  /* NOTHIN may be used with a radiation condition (NOTHIN|GRAVTY),  */
  /* where the value on the boundary is partially passive via        */
  /* radiation, or with data (NOTHIN|FILEIN), where the value on the */
  /* boundary is implicitly relaxed to data.                         */
  if (bcond & (LINEAR|NOTHIN)) {
    if (bcond & (FILEIN|CUSTOM|TIDALH|TIDALC|TIDEBC)) {
      /* Relax to data implicitly                                    */
      for (cc = sb; cc <= eb; cc++) {
        c = obc[cc];
	rts = dt / open->relax_time;
	vel[c] = fval[c] = (vel_t[c] + rts * fval[c]) / (1.0 + rts);
      }
    } else {
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	fval[c] += vel[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Flather radiation (velocity components)                         */
  if (bcond & FLATHR) {
    double spd;
    double depth;
    double sc = 1.0;
    double feta;
    double *hat = NULL;
    double *uval;
    int ci, ce, *cim = oi1;

    if (open->ocodex & (L_EDGE | R_EDGE))
      hat = window->h1acell;
    else if (open->ocodey & (B_EDGE | F_EDGE))
      hat = window->h2acell;
    if (open->ocodex & L_EDGE || open->ocodey & B_EDGE) {
      sc = -1.0;
      cim = obc;
    }
    if (bcond & (FILEIN|CUSTOM|LOCALN))
      uval = fval;
    else
      uval = newval;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = window->m2d[c1];
      ci = window->m2d[cim[cc]];
      depth = windat->eta[c2] - window->botz[c2] * wincon->Ds[c2];
      spd = (depth > 0.0) ? sqrt(wincon->g / depth) : 0.0;
      spd = (open->spf) ? open->sphase[cc] * hat[c2] / (depth * dt) : spd;
      /* Note eta is use at the interior cell (Palma & Matano (1998) */
      /* A.1 p 1340. These authors used GRAVTY on eta and tangential */
      /* component of velocity).                                     */
      feta = 0.0;
      if (open->bcond_ele & (FILEIN|CUSTOM))
	feta += open->transfer_eta[cc];
      if (open->bcond_ele & (TIDALH|TIDALC)) {
	double gmt = windat->t - wincon->tz;
	double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
	feta += ramp * csr_tide_eval(&open->tc, cc, gmt);
      }
      if (!(open->bcond_ele & (FILEIN|CUSTOM|TIDALH|TIDALC))) {
	ce = window->m2d[open->obc_t[cc]];
	feta += windat->eta[ce];
      }
      fval[c] = uval[c] + (sc * spd * (windat->eta[ci] - feta));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary value to the calculated condition.             */
  if (bcond == FILEIN || bcond == TIDEBC || 
      bcond == TIDALH || bcond == TIDALC ||
      bcond == CUSTOM || bcond == (FILEIN|TIDALH) || 
      bcond == (FILEIN|TIDALC) || bcond & FLATHR) {
    /* Active conditions - set directly to supplied data.            */
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      vel[c] = fval[c];
    }
  } else if (bcond & (FILEIN | CUSTOM | TIDEBC | TIDALH | TIDALC | 
		      LINEAR | NOTHIN)) {
    /* Partially passive conditions - relax with specified radiation */
    /* over the specified relaxation time.                           */
    /* Note : LINEAR for velocity and NOTHIN for elevation assume a  */
    /* radiation condition applied to a local solution.              */
    if (bcond & (MILLER | NOGRAD)) {
      /* Explicit fromulation.                                       */
      for (cc = sb; cc <= eb; cc++) {
        c = obc[cc];
        vel[c] =
          newval[c] + dt * (fval[c] - vel_t[c]) / open->relax_time;
      }
    } else if (bcond & (ORLANS | CAMOBR)) {
      /* Implicit formulations.                                      */
      for (cc = sb; cc <= eb; cc++) {
        c = obc[cc];
	rts = (cs[cc] > 0.0) ? open->relax_time : open->relax_timei;
        vel[c] = newval[c] + 2.0 * dt * (fval[c] - vel_b[c]) /
          (rts * (1.0 + cs[cc]));
        /*vel[c] = newval[c] + 2.0 * dt * (fval[c] - vel_t[c]) /
          (rts * (1.0 + cs[cc]));*/
      }
    } else if (bcond & GRAVTY) {
      /* Implicit formulation. Wave speed is always positive here    */
      /* (outgoing).                                                 */
      for (cc = sb; cc <= eb; cc++) {
        c = obc[cc];
        vel[c] = newval[c] + dt * (fval[c] - vel_t[c]) /
          (open->relax_time * (1.0 + cs[cc]));
      }
    } else if (bcond & RAYMND) {
      /* Use the adaptive relaxation (Marchesiello et al, 2001)      */
      /* where if cs>0 (outward propagation) the relax time is long  */
      /* (about 1 year) and if cs<0 no radiation is applied with     */
      /* strong nudging (relax of a few days). Note that cx = cy = 0 */
      /* for cx < 0 (Eqn 15), hence the RAYMND scheme returns the    */
      /* value vel_t[c] when cx < 0.                                 */
      for (cc = sb; cc <= eb; cc++) {
        c = obc[cc];
	rts = (cs[cc] > 0.0) ? open->relax_time : open->relax_timei;
        vel[c] = newval[c] + dt * (fval[c] - vel_t[c]) /
	                          (rts * (1.0 + cs[cc]));
      }
    }
  } else {
    /* Passive conditions - set to passive OBC.                      */
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      vel[c] = newval[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Impose a flow relaxtion condition if required                   */
  /* Relaxation to external solutions                                */
  if (rlxn) {
    double alpha;
    int bn;
    imap = open->nmap;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      bn = 2;
      while (c1 != imap[c1] && bn <= rlxn) {
        alpha = 1 - tanh(0.5 * (double)(bn - 1));
	sf = (code == 0) ? transfer[cc + (bn - 1) * eb] : vel[c];
	/*
	alpha = (double)(rlxn - bn + 1) / (double)rlxn;
	alpha *= alpha;
	vel[c1] = alpha * windat->rampval * vel[c] + (1.0 - alpha) * vel[c1];
	*/
	vel[c1] = alpha * windat->rampval * sf + (1.0 - alpha) * vel[c1];
        c1 = imap[c1];
        bn++;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Save the phase speed if required. Note: phase speed is bounded  */
  /* by the CFL condition (0 <= obc_phase <= 1). If waves are        */
  /* in-coming then obc_phase is negative, hence bounded to zero.    */
  /* Out-going waves have obc_phase > 0.                             */
  if (windat->obc_phase && vel_t == windat->eta && !open->adjust_flux) {
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      windat->obc_phase[c] = cs[cc];
    }
  }
}

/* END set_OBC()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays for the u1av calculation                           */
/*-------------------------------------------------------------------*/
void precalc_u1_2d(geometry_t *window,  /* Processing window */
                   window_t *windat,  /* Window data structure */
                   win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc;
  int xm1;                      /* Sparse coordinate at i-1 */
  int yp1;                      /* Sparse coordinate at j+1 */
  int xmyp1;                    /* Sparse coordinate at (i-1,j+1) */
  double *u2au1;
  int cw = 1;
  double top, bot;

  /* Eliminate dry cells from the cells to process list */
  for(cc = 1; cc <= wincon->vcs1; cc++) {
    c = window->m2d[wincon->s1[cc]];
    xm1 = window->xm1[c];
    top = max(windat->eta[c], windat->eta[xm1]);
    bot = DRY_FRAC * wincon->hmin + window->botzu1[c];
    if (top > bot) {
      wincon->s3[cw] = c;
      cw++;
    }
    else {
      windat->nu1av[c] = windat->u1av[c] = windat->u1avb[c] = 0.0;
    }
  }
  wincon->vcs = cw - 1;
  if(wincon->dolin_u1) {
    linear_bdry_cell(window, windat, wincon, wincon->s3, wincon->vcs,
		     wincon->linmask_u1, NULL);
    /*wincon->vcs = wincon->aclS;*/
  }

  /* Set pointers */
  u2au1 = wincon->w7;

  /* Copy the current velocities into the update array */
  memcpy(windat->nu1av, windat->u1avb, window->enonS * sizeof(double));

  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations */
  if (wincon->sigma) {
    for (cc = 1; cc <= window->v2_e1; cc++) {
      c = window->w2_e1[cc];
      windat->nu1av[c] *= mdxbs(window, windat, wincon, c);
    }
  }

  /* Get the u2 velocity at the e1 face */
  for (cc = 1; cc <= window->v2_e1; cc++) {
    c = window->w2_e1[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    u2au1[c] = 0.25 * (windat->u2av[c] + windat->u2av[xm1] +
                       windat->u2av[yp1] + windat->u2av[xmyp1]);
  }

  /* Re-calculate the u2 velocity at the e1 face for inner cyclic    */
  /* explicit maps. The coordinates systems are inverted for the     */
  /* ym1, xpym1 coordinates (i.e. the same magniture and direction   */
  /* velocity at c and ym1 is +ve at c and -ve at ym1, or vice       */
  /* versa) hence ym1 abd xpym1 must be subtracted from u1au2 (e.g.  */
  /* identical velocities at all four cells must add to that         */
  /* velocity rather than zero).                                     */
  for (cc = 1; cc <= window->neimS; cc++) {
    c = window->eims[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    u2au1[c] = 0.25 * (windat->u2av[c] + windat->u2av[yp1] -
		       windat->u2av[xm1] - windat->u2av[xmyp1]);
  }
}

/* END precalc_u1_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the values of a given array at previous time-steps */
/* on the boundaries and store in the old data structure.            */
/*-------------------------------------------------------------------*/
void store_old_bdry_values(double *F, /* df_variable_t array */
                           open_bdrys_t *open,  /* OBC structure */
                           bdry_old_values_t *old,  /* Old values data
                                                       structure */
                           int sb,  /* Start of OBC vector */
                           int eb,  /* End of OBC vector */
                           int *obc,  /* Open boundary vector */
                           int *oi1,  /* One cell interior to obc cells */
                           int *oi2 /* Two cell interior to obc cells */
  )
{
  int co, cc, c, c1, c2;        /* Sparse counters / coordinates */
  int tp, tm, tp1, tm1;
  for (cc = sb; cc <= eb; cc++) {
    co = cc - sb + 1;
    c = obc[cc];                /* Sparse coordinate on the boundary */
    c1 = oi1[cc];               /* One cell interior to c */
    c2 = oi2[cc];               /* Two cells interior to c */
    tp = open->tmpp[c];         /* Tangential coordinates (RAYMND) */
    tm = open->tmpm[c];
    tp1 = open->nmap[tp];
    tm1 = open->nmap[tm];
    old[co].t1i = old[co].ti;   /* F at time t-1 on the boundary */
    old[co].t1i1 = old[co].ti1; /* F at time t-1 one cell interior to c */
    old[co].t1i2 = old[co].ti2; /* F at t-1 two cells interior to c */
    old[co].ti = F[c];          /* F at time t on the boundary */
    old[co].ti1 = F[c1];        /* F at time t one cell interior to c */
    old[co].ti2 = F[c2];        /* F at time t two cells interior to c */
    old[co].tjp = F[tp];        /* F, time t tangential+1 on boundary */
    old[co].tjm = F[tm];        /* F, time t tangential-1 on boundary */
    old[co].tj1p = F[tp1];      /* F, time t tangential+1 on boundary-1 */
    old[co].tj1m = F[tm1];      /* F, time t tangential-1 on boundary-1 */
  }
}

/* END store_old_bdry_values()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e1 cell face           */
/*-------------------------------------------------------------------*/
void mdxs(geometry_t *window,   /* Processing window */
          win_priv_t *wincon    /* Window geometry / constants */
  )
{
  int c;                        /* Sparse coordinate/counter */
  int xm1;                      /* Sparse coordinate at i-1 */
  double h1av;

  for (c = 1; c <= window->enonS; c++) {
    xm1 = window->xm1[c];
    h1av =
      window->h1acell[xm1] / (window->h1acell[xm1] + window->h1acell[c]);
    wincon->mdx[c] =
      h1av * (wincon->Ds[c] - wincon->Ds[xm1]) + wincon->Ds[xm1];
  }
}

/* END mdxs()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e1 cell face at the    */
/* forward time step.                                                */
/*-------------------------------------------------------------------*/
double mdxns(geometry_t *window,  /* Processing window */
             window_t *windat,  /* Window data structure */
             win_priv_t *wincon,  /* Window geometry / constants */
             int c              /* Sparse coordiante */
  )
{
  int c2 = window->m2d[c];      /* 2D coordinate */
  int xm1 = window->xm1[c2];
  double h1av;

  if (!wincon->sigma)
    return (1.0);
  h1av =
    window->h1acell[xm1] / (window->h1acell[xm1] + window->h1acell[c2]);
  return (h1av *
          (windat->eta[c2] - wincon->Hs[c2] -
           (windat->eta[xm1] - wincon->Hs[xm1])) + windat->eta[xm1] -
          wincon->Hs[xm1]);
}

/* END mdxns()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e1 cell face at the    */
/* backward time step.                                               */
/*-------------------------------------------------------------------*/
double mdxbs(geometry_t *window,  /* Processing window */
             window_t *windat,  /* Window data structure */
             win_priv_t *wincon,  /* Window geometry / constants */
             int c              /* Sparse coordiante */
  )
{
  int c2 = window->m2d[c];      /* 2D coordinate */
  int xm1 = window->xm1[c2];
  double h1av;

  if (!wincon->sigma)
    return (1.0);
  h1av =
    window->h1acell[xm1] / (window->h1acell[xm1] + window->h1acell[c2]);
  return (h1av *
          (windat->etab[c2] - wincon->Hs[c2] -
           (windat->etab[xm1] - wincon->Hs[xm1])) + windat->etab[xm1] -
          wincon->Hs[xm1]);
}

/* END mdxbs()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e2 cell face           */
/*-------------------------------------------------------------------*/
void mdys(geometry_t *window,   /* Processing window */
          win_priv_t *wincon    /* Window geometry / constants */
  )
{
  int c;                        /* Sparse coordinate/counter */
  int ym1;                      /* Sparse coordinate at y-1 */
  double h2av;

  for (c = 1; c <= window->enonS; c++) {
    ym1 = window->ym1[c];
    h2av =
      window->h2acell[ym1] / (window->h2acell[ym1] + window->h2acell[c]);
    wincon->mdy[c] =
      h2av * (wincon->Ds[c] - wincon->Ds[ym1]) + wincon->Ds[ym1];
  }
}

/* END mdys()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e2 cell face at the    */
/* forward time step.                                                */
/*-------------------------------------------------------------------*/
double mdyns(geometry_t *window,  /* Processing window */
             window_t *windat,  /* Window data structure */
             win_priv_t *wincon,  /* Window geometry / constants */
             int c              /* Sparse coordiante */
  )
{
  int c2 = window->m2d[c];      /* 2D coordinate */
  int ym1 = window->ym1[c2];
  double h2av;

  if (!wincon->sigma)
    return (1.0);
  h2av =
    window->h2acell[ym1] / (window->h2acell[ym1] + window->h2acell[c2]);
  return (h2av *
          (windat->eta[c2] - wincon->Hs[c2] -
           (windat->eta[ym1] - wincon->Hs[ym1])) + windat->eta[ym1] -
          wincon->Hs[ym1]);
}

/* END mdyns()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an e2 cell face at the    */
/* backward time step.                                               */
/*-------------------------------------------------------------------*/
double mdybs(geometry_t *window,  /* Processing window */
             window_t *windat,  /* Window data structure */
             win_priv_t *wincon,  /* Window geometry / constants */
             int c              /* Sparse coordiante */
  )
{
  int c2 = window->m2d[c];      /* 2D coordinate */
  int ym1 = window->ym1[c2];
  double h2av;

  if (!wincon->sigma)
    return (1.0);
  h2av =
    window->h2acell[ym1] / (window->h2acell[ym1] + window->h2acell[c2]);
  return (h2av *
          (windat->etab[c2] - wincon->Hs[c2] -
           (windat->etab[ym1] - wincon->Hs[ym1])) + windat->etab[ym1] -
          wincon->Hs[ym1]);
}

/* END mdybs()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void vel_u2av_update(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc;                    /* Sparse coodinate / counter */
  int xp1, xm1;                 /* Sparse locations at i+1, i-1 */
  int yp1, ym1;                 /* Sparse locations at j+1, j-1 */
  int xpym1;                    /* Sparse location at (i-1,y+1) */
  double pgt = 0.0;             /* Pressure gradient term */
  double cot = 0.0;             /* Coriolis term */
  double bft;                   /* Bottom friction term */
  double rst = 0.0;             /* Radiation stress term */
  double *tzp;                  /* Surface height array */
  double *depth;                /* Depth of the water column */
  double midy;                  /* Water depth at the cell face */
  double val;                   /* Dummy */
  double u1au2;                 /* u1 value at the e2 face */
  double botu1;                 /* Bottom velocity at e1 face */
  double botu2;                 /* Bottom velocity at e2 face */
  double Cdu2;                  /* Bottom drag coefficient at the e2 face */
  double rho0 = 1024.0;         /* Reference density */
  int dbc;                      /* Debugging coordinate */

  dbc = (wincon->dbc) ? window->m2d[wincon->dbc] : 0;
  depth = windat->depth_e2;

  /* tzp is a pointer which points to either topz[][] or eta[][].  */
  /* topz[][] only gets updated once every 3d step. For the non- */
  /* linear case, we need to use eta instead, which gets updated */
  /* every 2d step. This value is only used to calculate transport */
  /* (nothing to do with the surface slope term, which always uses */
  /* eta).  */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u2 calculation */
  precalc_u2_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across u2 OBC's */
  set_map_e2(window);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection */
  if (!(wincon->u2av_f & ADVECT)) {
    if(wincon->momsc2d & ANGULAR) {
      if(!(wincon->ic % 2))
	advect_u2_2d_ang_flux_b(window, windat, wincon);
      else
	advect_u2_2d_ang_flux_f(window, windat, wincon);
    } else if (wincon->momsc2d & LAGRANGE)
      semi_lagrange_u2av(window, windat, wincon);
    else
      advect_u2_2d(window, windat, wincon);
  }

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included) */
  wincon->hor_mix->u2av(window, windat, wincon, tzp);

  /*-----------------------------------------------------------------*/
  /* Add the body forces and step forward in time */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    xpym1 = window->xpym1[c];

    midy = wincon->mdy[c];

    /*---------------------------------------------------------------*/
    /* Pressure gradient term. Note : differences bracketed to */
    /* preserve accuracy.  */
    if (!(wincon->u2av_f & PRESS_BT))
      pgt = -((windat->eta[c] - windat->eta[ym1]) * 
	      wincon->topdensu2[c] +
	      (windat->patm[c] - windat->patm[ym1])) / 
	wincon->densavu2[c];

    /*---------------------------------------------------------------*/
    /* Coriolis */
    u1au2 = wincon->w7[c];      /* Set in precalc_u2_2d() */
    if (!(wincon->u2av_f & CORIOLIS))
      cot = wincon->u2c5[c] * u1au2;

    /*---------------------------------------------------------------*/
    /* Bottom friction */
    u1au2 = 0.25 * (windat->u1avb[c] + windat->u1avb[ym1] +
                    windat->u1avb[xpym1] + windat->u1avb[xp1]);
    botu1 = 0.25 * (windat->u1bot[c] + windat->u1bot[ym1] +
                    windat->u1bot[xpym1] + windat->u1bot[xp1]);

    botu2 = windat->u2avb[c] + windat->u2bot[c];
    botu1 = u1au2 + botu1;

    Cdu2 = 0.5 * (wincon->Cd[ym1] + wincon->Cd[c]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu2 * max(wincon->uf, val);
    /* Truncate to ensure stability */
    if (val > depth[c] * midy / windat->dt2d)
      val = depth[c] * midy / windat->dt2d;
    /* Note: depth=1 in the sigma case */
    bft = -val * botu2 / depth[c];
    if (c == dbc) {
      wincon->b1 = pgt;
      wincon->b2 = cot;
      wincon->b3 = bft;
    }

    /*---------------------------------------------------------------*/
    /* Calculate nu1av value */
    /* SIGMA : multiply new velocity by depth (midy) */
    windat->nu2av[c] +=
      windat->dt2d * (midy * (pgt + cot) + bft + wincon->u2inter[c]);
  }
  debug_c(window, D_VA, D_POST);

  /*-----------------------------------------------------------------*/
  /* Add the non-linear terms */
  /* SIGMA : No extra terms to include for the sigma case.  */
  if (wincon->nonlinear && !wincon->sigma) {
    for (cc = 1; cc <= window->v2_e2; cc++) {
      c = window->w2_e2[cc];
      ym1 = window->ym1[c];
      windat->nu2av[c] -= windat->dt2d * windat->u2av[c] *
        (windat->detadt[ym1] + windat->detadt[c]) /
        (2.0 * max(depth[c], wincon->hmin));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add the radiation stresses                                      */
  if (wincon->waves & WAVE_FOR) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s3[cc];
      ym1 = window->ym1[c];
      rst = windat->dt2d * windat->wave_Fy[c] / 
	                   (rho0 * max(depth[c], wincon->hmin));
      if (windat->u2_rad) windat->u2_rad[c] = rst;
      windat->nu2av[c] += rst;
    }
  }
  else if (wincon->waves & TAN_RAD) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s3[cc];
      rst = windat->dt2d * (windat->wave_Syx[window->xp1[c]] - 
			    windat->wave_Syx[window->xm1[c]]) / 
	    (window->h1au1[c] + window->h1au1[window->xp1[c]]);
      if (windat->u2_rad) windat->u2_rad[c] = rst;
      windat->nu2av[c] += rst;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Restrict flow if column nearly dry */
  for (cc = 1; cc <= window->v2_e2; cc++) {
    c = window->w2_e2[cc];
    ym1 = window->ym1[c];
    if (windat->nu2av[c] > 0.0 &&
        (val = (tzp[ym1] - window->botz[ym1])) < wincon->hmin)
      windat->nu2av[c] *= max(val, 0.0) / wincon->hmin;
    else if (windat->nu2av[c] < 0.0 &&
             (val = (tzp[c] - window->botz[c])) < wincon->hmin)
      windat->nu2av[c] *= max(val, 0.0) / wincon->hmin;
  }

  /*-----------------------------------------------------------------*/
  /* Blend velocities over a subregion if required                   */
  blend_vel(window, windat, wincon, U22D, windat->nu2av);

  /*-----------------------------------------------------------------*/
  /* Calculate u1av boundary values */
  bdry_u2_2d(window, windat, wincon);
  debug_c(window, D_VA, D_BDRY);

  /*
  for (cc = 1; cc <= window->neimS; cc++) {
    int cs, cd;
    double vs, vd;
    val = 0.5 * (windat->nu2av[cs] - windat->nu2av[cd]);    
    cs = window->eims[cc];
    cd = window->eimd[cc];
    vs = windat->nu2av[cs];
    vd = windat->nu2av[cd] ;
    val = fabs(0.5 * (vs - vd));
    windat->nu2av[cs] = vs > 0 ? val : -val;
    windat->nu2av[cd] = vd > 0 ? val : -val;
  }
  */
}

/* END vel_u2av_update()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 2D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/*-------------------------------------------------------------------*/
void advect_u2_2d(geometry_t *window, /* Processing window           */
                  window_t *windat,   /* Window data structure       */
                  win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc;                /* Sparse coordinate / counter           */
  int xp1, xm1;             /* Sparse coordinate at i+1, i-1         */
  int yp1, ym1;             /* Sparse coordinate at j+1, j-1         */
  int xpym1;                /* Sparse coordinate at (i-1,j+1)        */
  int *ctp;                 /* Cells to process vector               */
  int vcs;                  /* Number of cells to process            */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double u1val;             /* Cell centered u1 value                */
  double u2val;             /* Cell cornered u2 value                */
  double h2;                /* Face centered grid spacing            */
  double h1av, h2av;        /* Mean cell spacings x and y directions */
  double m, d2;             /* Dummy variables                       */
  double *Du2;              /* Dummy for e1 advective fluxes         */
  double *vel;              /* Velocity to use in spatial gradients  */
  double *area;             /* Cell area                             */
  double *u2sh1h2_av;       /* u2 * u2 * h1 * h2                     */
  double *u1u2h2s_av;       /* u1 * u2 * h2 * h2                     */
  double *depth;            /* Depth of the water column             */
  double *dvel;             /* Velocity to interpolate               */
  double *vface;            /* Velocity centered on the grid face    */
  double *cn;               /* Courant number                        */
  double midy;              /* Depth at the e1 face                  */
  double *u1au2;            /* u1 value at the e2 face               */
  double lr;                /* Ratio of dtf : dt2d                   */

  if (!wincon->nonlinear)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  Du2 = wincon->d1;
  u2sh1h2_av = wincon->d2;
  u1u2h2s_av = wincon->d3;
  area = window->cellarea;
  depth = windat->depth_e2;
  vel = windat->u2av;
  u1au2 = wincon->w7;
  if(wincon->dolin_u2) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d2 = windat->dt / windat->dtf2;
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu2 / d2;       /* 2D forward time-step           */

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  if(windat->dtu1 == windat->dtf)
    dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping  */
  else
    dtu *= 2.0;           /* 2D leapfrog time-step with sub-stepping */
  lr = windat->dtf2 / windat->dt2d;
  dtm = dtu * lr;
  trem = windat->dt2d;
  tremf = windat->dtf2;
  while (trem > 0) {
    if (trem < dtu)
      dtu = trem;

    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Initialise                                                  */
      memset(Du2, 0, window->sgsizS * sizeof(double));
      memset(u2sh1h2_av, 0, window->sgsizS * sizeof(double));
      memset(u1u2h2s_av, 0, window->sgsizS * sizeof(double));

      /*-------------------------------------------------------------*/
      /* Get the tracer values in multi-dt auxiliary cells. This is  */
      /* a linear interpolation between the value at the start of    */
      /* the timestep and the end of longer timesteps.               */
      set_multidt_tr_2d(window, windat, trem, vel, windat->u2av_as,
                        windat->u2av_ae);

      /*-------------------------------------------------------------*/
      /* Loop over the u1 points to calculate depths                 */
      for (cc = 1; cc <= window->x2_e2; cc++) {
        c = window->w2_e2[cc];
        Du2[c] = depth[c] * vel[c];
      }
      if (wincon->momsc2d & VANLEER) {
	for (cc = window->b2_e2 + 1; cc <= window->n2_e2; cc++) {
	  c = window->w2_e2[cc];
	  ym1 = window->ym1[c];
	  Du2[ym1] = depth[ym1] * vel[ym1];
	}
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres                      */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          yp1 = window->yp1[c];
          u2val = 0.5 * (vel[yp1] + vel[c]);

          /* SIGMA : Multiply by the cell centre depth               */
          if (u2val > 0)
            m = Du2[c] * wincon->mdy[c];
          else
            m = Du2[yp1] * wincon->mdy[yp1];
          u2sh1h2_av[c] = m * u2val * area[c];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          yp1 = window->yp1[c];
          u2val = 0.5 * (vel[yp1] + vel[c]);

          /* SIGMA : Multiply by the cell centre depth               */
          u2sh1h2_av[c] =
            u2val * area[c] * 0.5 * (Du2[c] * wincon->mdy[c] +
                                     Du2[yp1] * wincon->mdy[yp1]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          yp1 = window->yp1[c];
          vface[c] = 0.5 * (vel[yp1] + vel[c]);
          cn[c] = vface[c] * dtu / window->h2acell[c];
          dvel[c] = Du2[yp1] * wincon->mdy[yp1];
        }
	/* Ghost cells 2 cells in */
        for (cc = window->b2_e2 + 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          ym1 = window->ym1[c];
          dvel[ym1] = Du2[c] * wincon->mdy[c];
	}
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array dvel.                                         */
        van_leer_do(u2sh1h2_av, dvel, vface, cn, 1, window->n2_e2,
                    window->w2_e2, window->yp1, window->ym1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          u2sh1h2_av[c] *= vface[c] * area[c];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners                      */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          h2 = 0.5 * (window->h2au1[ym1] + window->h2au1[c]);
          u1val = 0.5 * (windat->u1av[ym1] * wincon->mdx[ym1] +
                         windat->u1av[c] * wincon->mdx[c]);
          if (u1val > 0)
            m = Du2[xm1];
          else
            m = Du2[c];
          u1u2h2s_av[c] = m * h2 * h2 * u1val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];

          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          h2 = 0.5 * (window->h2au1[ym1] + window->h2au1[c]);
          h1av =
            window->h1au2[xm1] / (window->h1au2[xm1] + window->h1au2[c]);
          h2av =
            window->h2au1[ym1] / (window->h2au1[ym1] + window->h2au1[c]);
          /* SIGMA : Multiply by the cell corner depth               */
          u1val = h2av * (windat->u1av[c] * wincon->Ds[c] -
                          windat->u1av[ym1] * wincon->Ds[ym1]) +
            windat->u1av[ym1] * wincon->Ds[ym1];
          m = h1av * (Du2[c] - Du2[xm1]) + Du2[xm1];
          u1u2h2s_av[c] = m * h2 * h2 * u1val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          ym1 = window->ym1[c];
          vface[c] = 0.5 * (windat->u1av[ym1] + windat->u1av[c]);
          cn[c] =
            2.0 * vface[c] * dtu / (window->h1au1[c] + window->h1au1[ym1]);
          dvel[c] = Du2[c] * wincon->mdy[c];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array dvel.                                         */
        van_leer_do(u1u2h2s_av, dvel, vface, cn, 1, window->n2_e2,
                    window->w2_e2, window->xp1, window->xm1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n2_e2; cc++) {
          c = window->w2_e2[cc];
          u1u2h2s_av[c] *= vface[c] * area[c];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the u2 updated solution                                 */
      for (cc = 1; cc <= vcs; cc++) {
        c = ctp[cc];
        ym1 = window->ym1[c];
        xp1 = window->xp1[c];
        xpym1 = window->xpym1[c];
        midy = wincon->mdy[c];

        /* SIGMA : Multiply by cell depth                            */
        d2 = wincon->u2c1[c] *
          (u2sh1h2_av[c] - u2sh1h2_av[ym1] + u1u2h2s_av[xp1] -
           u1u2h2s_av[c]) / max(depth[c], wincon->hmin) +
          wincon->u2c3[c] * u1au2[c] * u1au2[c] * midy +
          wincon->u2c4[c] * vel[c] * vel[c] * midy;

        windat->nu2av[c] += d2 * dtu;

        /* Integrate the dispersion terms over the sub-timestep.     */
        /* u2adv is divided by dt to get the mean over the 2D mode   */
        /* in vel_u2_update().                                       */
        wincon->u2adv[c] -= (d2 * dtm);
      }
    }

    /* Decrement the time remaining                                  */
    trem -= dtu;
    tremf -= dtm;
    dtu = dtm / lr;
    if (trem < dtu)
      dtu = trem;
    if (tremf < dtm)
      dtm = tremf;
    if(!wincon->sigma)
      vel = windat->nu2av;

    if (wincon->nbl2) {
      vel2D_lbc(vel, window->nbpte2S, window->nbe2S,
		window->bpte2S, window->bine2S, wincon->slip);
      blend_vel(window, windat, wincon, U22D, vel);
    }

  }
  debug_c(window, D_VA, D_ADVECT);
}

/* END advect_u2_2d()                                          (slf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the correct map across open boundaries for u2      */
/*-------------------------------------------------------------------*/
void set_map_e2(geometry_t *window  /* Window data structure */
  )
{
  int n;
  int cc, cs, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to northern open boundaries */
  /* (set in the tracer routine) to be self-mapping on the open */
  /* boundary.  */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodey & F_EDGE && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e2; cc++) {
        c = open[n]->obc_e2[cc];
        cs = window->ym1[c];
        window->yp1[cs] = c;
      }
    }

    /* Set maps to be self-mapping interior to eastern open */
    /* boundaries.  */
    if (open[n]->ocodex & R_EDGE && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e1; cc++) {
        c = open[n]->obc_e1[cc];
        cs = window->xm1[c];
        window->xp1[cs] = cs;
      }
    }
  }
}

/* END set_map_e2()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the u2 velocities on all open boundaries in a  */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_u2_2d(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int n;                        /* Counters */
  int c, cc, c1;                /* Sparse coordinates / counters */
  double *depth;                /* Depth of the water column */
  double midy, midy1;           /* Depth at the e1 face */
  open_bdrys_t **open = window->open;

  depth = windat->depth_e2;

  /*-----------------------------------------------------------------*/
  /* Set the u2av velocity where u2av is normal to the boundary */
  for (n = 0; n < window->nobc; n++) {
    /* Get the u2av velocity derived from a no gradient condition */
    /* on u2flux.  */
    if (open[n]->bcond_nor2d & NOGRAD) {
      for (cc = 1; cc <= open[n]->no2_e2; cc++) {
        c = open[n]->obc_e2[cc];
        c1 = open[n]->oi1_e2[cc];
        midy = wincon->mdy[c];
        midy1 = wincon->mdy[c1];

        /* Set the no gradient condition on velocity */
        windat->nu2av[c] =
          windat->nu2av[c1] * depth[c1] * midy1 * window->h1au2[c1] /
          (max(depth[c], wincon->hmin) * midy * window->h1au2[c]);
      }
    }
    /* Set the u2av normal velocity for all other BC's.  */
    else {
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e2,
              open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
              open[n]->cyc_e2, windat->nu2av, windat->u2av, 
	      windat->u2avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau2av, 
	      open[n]->transfer_u2av, open[n]->relax_zone_nor, U2BDRY);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the u2av velocity where u2av is tangential to the boundary */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	    open[n]->to2_e2, open[n]->obc_e2, open[n]->oi1_e2,
	    open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2av,
	    windat->u2av, windat->u2avb, open[n]->bcond_tan2d, 
	    windat->dtf2, &open[n]->datau2av, open[n]->transfer_u2av, 
	    open[n]->relax_zone_tan, U2BDRY);	    
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      int *mi = (open[n]->ocodey & B_EDGE) ? window->ym1 : window->yp1;
      for (cc = 1; cc <= open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	c1 = open[n]->obc_t[cc];
	windat->nu2av[c1] = windat->nu2av[c];
      }
    }
  }

  /* Rescale for sigma active OBCs */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->type & U2BDRY && open->bcond_nor2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu2av, mdyns,
			1, open->no2_e2, open->obc_e2);
      if (open->type & U1BDRY && open->bcond_tan2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu2av, mdyns,
			open->no2_e2 + 1, open->to2_e2, open->obc_e2);
    }
  }
}

/* END bdry_u2_2d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays for the u2av calculation                           */
/*-------------------------------------------------------------------*/
void precalc_u2_2d(geometry_t *window,  /* Processing window */
                   window_t *windat,  /* Window data structure */
                   win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc;
  int xp1;                      /* Sparse coordinate at i+1 */
  int ym1;                      /* Sparse coordinate at j-1 */
  int xpym1;                    /* Sparse coordinate at (i+1,j-1) */
  double *u1au2;
  int cw = 1;
  double top, bot;

  /* Eliminate dry cells from the cells to process list */
  for(cc = 1; cc <= wincon->vcs2; cc++) {
    c = window->m2d[wincon->s2[cc]];
    ym1 = window->ym1[c];
    top = max(windat->eta[c], windat->eta[ym1]);
    bot = DRY_FRAC * wincon->hmin + window->botzu2[c];
    if (top > bot) {
      wincon->s3[cw] = c;
      cw++;
    }
    else {
      windat->nu2av[c] = windat->u2av[c] = windat->u2avb[c] = 0.0;
    }
  }
  wincon->vcs = cw - 1;
  if(wincon->dolin_u2) {
    linear_bdry_cell(window, windat, wincon, wincon->s3, wincon->vcs,
		     wincon->linmask_u2, NULL);
    /*wincon->vcs = wincon->aclS;*/
  }

  /* Set pointers */
  u1au2 = wincon->w7;

  /* Copy the current velocities into the update array */
  memcpy(windat->nu2av, windat->u2avb, window->sgsizS * sizeof(double));

  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations */
  if (wincon->sigma) {
    for (cc = 1; cc <= window->v2_e2; cc++) {
      c = window->w2_e2[cc];
      windat->nu2av[c] *= mdybs(window, windat, wincon, c);
    }
  }

  /* Get the u1 velocity at the e2 face */
  for (cc = 1; cc <= window->v2_e2; cc++) {
    c = window->w2_e2[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    u1au2[c] = 0.25 * (windat->u1av[c] + windat->u1av[ym1] +
                       windat->u1av[xpym1] + windat->u1av[xp1]);
  }

  /* Re-calculate the u1 velocity at the e2 face for inner cyclic    */
  /* explicit maps. The coordinates systems are inverted for the     */
  /* ym1, xpym1 coordinates (i.e. the same magniture and direction   */
  /* velocity at c and ym1 is +ve at c and -ve at ym1, or vice       */
  /* versa) hence ym1 abd xpym1 must be subtracted from u1au2 (e.g.  */
  /* identical velocities at all four cells must add to that         */
  /* velocity rather than zero).                                     */
  for (cc = 1; cc <= window->neimS; cc++) {
    c = window->eims[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    u1au2[c] = 0.25 * (windat->u1av[c] + windat->u1av[xp1] - 
		       windat->u1av[ym1] - windat->u1av[xpym1]);
  }
}

/* END precalc_u2_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to update the surface elevation                           */
/*-------------------------------------------------------------------*/
void eta_step(geometry_t *window,   /* Window geometry               */
              window_t *windat,     /* Window data                   */
              win_priv_t *wincon,   /* Window constants              */
              int ic                /* 2D mode step                  */
  )
{
  int c, cc;                    /* Sparse coordinate / counter       */
  double colflux;               /* Velocity transport divergence     */
  double *u1flux = wincon->d2;  /* Flux in e1 direction              */
  double *u2flux = wincon->d3;  /* Flux in e2 direction              */

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->sgsizS * sizeof(double));
  memset(u2flux, 0, window->sgsizS * sizeof(double));
  /* Calculate the flux at e1 wet and boundary cells                 */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    u1flux[c] = windat->u1av[c] * windat->depth_e1[c] * window->h2au1[c] *
      wincon->mdx[c] * windat->dt2d;
  }
  
  /* Calculate the flux at e2 wet and boundary cells                 */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    u2flux[c] = windat->u2av[c] * windat->depth_e2[c] * window->h1au2[c] *
      wincon->mdy[c] * windat->dt2d;
  }

  /* Calculate the flux at auxiliary cells using the tracer arrays.  */
  /* These arrays are used since some front and right cell faces     */
  /* required for cell center operations are not included in the e1  */
  /* and e2 cells to process vectors. e.g. e1/e2 cells on left/back  */
  /* faces adjacent to land are ghost cells, but wet cell centers.   */
  /* Auxiliary cells to the right/front of these cells are not       */
  /* included in the e1/e2 cells to process vectors, but are         */
  /* included in the cell center cells to process.                   */
  for (cc = window->b2_t + 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    u1flux[c] = windat->u1av[c] * windat->depth_e1[c] * window->h2au1[c] *
      wincon->mdx[c] * windat->dt2d;
    u2flux[c] = windat->u2av[c] * windat->depth_e2[c] * window->h1au2[c] *
	wincon->mdy[c] * windat->dt2d;
  }

  /* Use the fluxes summed over the fine grid faces for the coarse   */
  /* grid auxiliary fluxes if required.                              */
  if (window->zoomf > 1) {
    for (cc = window->b2_t + 1; cc <= window->a2_t; cc++) {
      c = window->w2_t[cc];
      u1flux[c] = windat->u1flux_z[c];
      u2flux[c] = windat->u2flux_z[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the fluxes for correction of the 3D mode. These fluxes must */
  /* be such that eta(t+1)=oldeta+divergence(fluxes) otherwise the   */
  /* divergence of vertical and horizontal 3D fluxes (using the 2D   */
  /* corrected velocities) in the surface layer will not result in   */
  /* the correct elevation change (eta(t+1)-oldeta) and the tracer   */
  /* equation will not be conservative at the surface. If every 2nd  */
  /* flux is saved (ic is odd) then fluxes are consistent - this can */
  /* only happen for odd iratios.                                    */
  if (((ic + 1) % 2) == 0) {
    for (cc = 1; cc <= window->b2_e1; cc++) {
      c = window->w2_e1[cc];
      windat->u1flux[c] += u1flux[c];
    }
    for (cc = 1; cc <= window->b2_e2; cc++) {
      c = window->w2_e2[cc];
      windat->u2flux[c] += u2flux[c];
    }
    /*
    for (cc = window->b2_t + 1; cc <= window->a2_t; cc++) {
      c = window->w2_t[cc];
      windat->u1flux[c] += u1flux[c];
      windat->u2flux[c] += u2flux[c];
    }
    */
  }

  /*-----------------------------------------------------------------*/
  /* Debugging info                                                  */
  if (wincon->dbc && window->wn == wincon->dbw) {
    c = wincon->dbc;
    wincon->b1 = u1flux[window->xp1[c]];
    wincon->b2 = u1flux[c];
    wincon->b3 = u2flux[window->yp1[c]];
    wincon->b4 = u2flux[c];
    debug_c(window, D_ETA, D_POST);
  }

  /*-----------------------------------------------------------------*/
  /* Update the elevation.                                           */
  /* Note : elevations are calculated on the boundaries here and     */
  /* overwritten in the boundary routine. A boundary condition of    */
  /* NOTHIN will use the elevations calculated here.                 */
  set_map_eta(window);
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = u1flux[window->xp1[c]] - u1flux[c] +
      u2flux[window->yp1[c]] - u2flux[c];

    /* Add source/sink flows to colflux. Note that colflux is the    */
    /* accumulated volume flowing out of this column, whereas        */
    /* waterss2d is the flow (m3 s-1) into the column. Hence the     */
    /* minus sign.  */
    colflux -= windat->waterss2d[c] * windat->dt2d;

    /* Calculate new etat value                            (etamark) */
    wincon->neweta[c] = max(windat->etab[c] - colflux / window->cellarea[c],
			    window->botz[c]);

    /* Rate of change of eta                                         */
    windat->detadt[c] =
      (wincon->neweta[c] - windat->etab[c]) / windat->dt2d;

    /* Check for large surface elevation value, or NaN               */
    /*
    if (fabs(wincon->neweta[c]) > wincon->etamax) {
      hd_quit_and_dump
        ("etastep: Surface exceeded ETAMAX (%f) at i=%ld j=%ld t=%f days\n",
         wincon->neweta[c], window->s2i[c], window->s2j[c],
         windat->t / 86400);
    }
    if (isnan(wincon->neweta[c])) {
      if (window->nwindows == 1)
	hd_quit_and_dump("etastep: NaN at i=%ld j=%ld t=%f days\n",
			 window->s2i[c], window->s2j[c], windat->t / 86400);
      else
	hd_quit_and_dump("etastep: NaN at i=%ld j=%ld t=%f days : window %d\n",
			 window->s2i[c], window->s2j[c], windat->t / 86400, 
			 window->wn);
    }
    */
  }

  if (wincon->dozoom != NOZOOM) {
    memcpy(windat->u1flux_z, u1flux, window->sgsizS * sizeof(double));
    memcpy(windat->u2flux_z, u2flux, window->sgsizS * sizeof(double));
  }

  /*-----------------------------------------------------------------*/
  /* Adjust the updated elevation due to eta relaxation.             */
  if (wincon->etarlx & (RELAX|ALERT)) {
    double rr = 0.0;
    /*if (wincon->etarlx & RELAX) rr = windat->etarlxtc;*/
    if (wincon->etarlx & RELAX) {
      relax_info_t *rlx = windat->eta_rlx;
      if (rlx->tctype & (RLX_CONS|RLX_FILE))
	rr = rlx->rate;
      else rr = 1.0;
    }

    /* This allows hard relaxation (starting at 2D timestep) over    */
    /* ramp period for ROAM startups.                                */
    if (wincon->rampf & ETA_RELAX && windat->rampval < 1.0)
      rr = windat->dt2d / (1.0 - windat->rampval);
    if (rr) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	/*do_eta_relax(window, c, windat->etarlxtc, 0);*/
	do_eta_relax(window, c, rr, wincon->tide_r);
      }
    }
  }
  /*-----------------------------------------------------------------*/
  /* Adjust the updated elevation due to eta increments.             */
  if (wincon->etarlx & INCREMENT)
    do_eta_increment(window);
}

/* END eta_step()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the elevations on all open boundaries in a     */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_eta(geometry_t *window, /* Processing window */
              window_t *windat, /* Window data structure */
              win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int n, m;                      /* Counters */
  int c, cc;

  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    scale_details_t *scale = open[n]->sdata_e;

    if (open[n]->bcond_ele != NOTHIN) {

      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_t,
              open[n]->obc_t, open[n]->oi1_t, open[n]->oi2_t,
              open[n]->cyc_t, wincon->neweta, windat->eta,
	      windat->etab, open[n]->bcond_ele,
              windat->dtf2, &open[n]->etadata, 
	      open[n]->transfer_eta, open[n]->relax_zone_ele, ETASPEC);

      /* Correct for atmospheric pressure */
      if (open[n]->inverse_barometer && !open[n]->adjust_flux) {
	double ramp = (wincon->rampf & INV_BARO) ? windat->rampval : 1.0;
        for (cc = 1; cc <= open[n]->no2_t; cc++) {
          c = open[n]->obc_t[cc];
          wincon->neweta[c] += ramp * (wincon->ambpress - windat->patm[c]) /
            (windat->dens[c] * wincon->g);
        }
      }

      /* Set eta to zero at boundary cells one cell deep */
      if (open[n]->stagger & INFACE) {
	int c1, c2;
        for (cc = 1; cc <= open[n]->no2_t; cc++) {
          c = open[n]->obc_t[cc];
          c1 = open[n]->oi1_t[cc];
          c2 = open[n]->oi2_t[cc];
	  if (c1 == c2) wincon->neweta[c] = 0.0;
	}
      }

      /* Multiply by the ramping factor */
      /*
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
        c = open[n]->obc_t[cc];
        wincon->neweta[c] *= windat->rampval;
      }
      */

      /* Scaling                                                     */
      if (!open[n]->adjust_flux) {
	if (scale->type == (TRSC_SUM | TRSC_NUM)) {
	  double fact = scale->fact;
	  for (cc = 1; cc <= open[n]->no2_t; cc++) {
	    c = open[n]->obc_t[cc];
	    wincon->neweta[c] += fact;
	  }
	} else if (scale->type == (TRSC_SUM | TRSC_TRA)) {
	  int trn = scale->ntr;
	  for (cc = 1; cc <= open[n]->no2_t; cc++) {
	    c = open[n]->obc_t[cc];
	    wincon->neweta[c] += windat->tr_wcS[trn][c];
	  }
	}
	if (scale->type == (TRSC_PCT | TRSC_NUM)) {
	  double fact = scale->fact;
	  for (cc = 1; cc <= open[n]->no2_t; cc++) {
	    c = open[n]->obc_t[cc];
	    wincon->neweta[c] *= fact;
	  }
	} else if (scale->type == (TRSC_PCT | TRSC_TRA)) {
	  int trn = scale->ntr;
	  for (cc = 1; cc <= open[n]->no2_t; cc++) {
	    c = open[n]->obc_t[cc];
	    wincon->neweta[c] *= windat->tr_wcS[trn][c];
	  }
	}
      }

      /* Average the sea level over intersecting OBCs if required   */
      if (open[n]->options & OP_OBCCM)
	average_OBC_corner(window, open[n], wincon->neweta, 0);

      /* Get the rate of change of elevation on the open boundaries */
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
        c = open[n]->obc_t[cc];
        windat->detadt[c] =
          (wincon->neweta[c] - windat->etab[c]) / windat->dt2d;
      }
    }

    /* Set the elevation at the normal velocity open boundary */
    /* locations equal to that at the elevation open boundary */
    /* location. For left and back boundaries these locations are */
    /* the same (i.e. the same value is overwritten - inefficient */
    /* but acceptable) and for right and front boundaries the */
    /* elevation boundary location is one cell into the interior */
    /* from the velocity boundary cell. This is required to set */
    /* dzu1[] and dzu2[] on the velocity open boundaries, which is */
    /* in turn used in velocity_adjust() to set vertical integrals */
    /* of 3D velocity equal to the 2D velocity.  */
    reset_bdry_eta(window, windat, wincon, open[n], wincon->neweta);
    /*
    if (open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	wincon->neweta[c] = wincon->neweta[open[n]->obc_t[cc]];
	for (m = 0; m < open[n]->bgz; m++) {
	  c = open[n]->omap[c];
	  wincon->neweta[c] = wincon->neweta[open[n]->obc_t[cc]];
	}
      }
      for (cc = 1; cc <= open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	wincon->neweta[c] = wincon->neweta[open[n]->obc_t[cc]];
	for (m = 1; m < open[n]->bgz; m++) {
	  c = open[n]->omap[c];
	  wincon->neweta[c] = wincon->neweta[open[n]->obc_t[cc]];
	}
      }
      for (cc = 1; cc <= open[n]->no2_a; cc++) {
	c = open[n]->obc_a[cc];
	for (m = 1; m < open[n]->bgz; m++) {
	  c = open[n]->omap[c];
	  wincon->neweta[c] = wincon->neweta[open[n]->obc_t[cc]];
	}
      }
    }
    */

    /* Set relaxation on the boundary                                */
    if (open[n]->relax_ele) {
      double rb = open[n]->rele_b;
      double ri = open[n]->rele_i;
      int nr = open[n]->relax_ele;
      int *imap = open[n]->nmap;
      double rc;
      int bn;
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
	c = open[n]->obc_t[cc];
	bn = 1;
	while (c != imap[c] && bn <= nr) {
	/* Linear */
	rc = (double)(bn - 1) * (ri - rb) / (double)(nr - 1) + rb;
	/* Hyperbolic tangent */
        /* rc = (1 - tanh(0.5 * (double)(bn - 1))) * rb; */
	do_eta_relax(window, c, rc, wincon->tide_r);
        c = imap[c];
        bn++;
	}
      }
    }
  }
}

/* END bdry_eta()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the elevation at boundary cells (no gradient)                 */
/*-------------------------------------------------------------------*/
void reset_bdry_eta(geometry_t *window, /* Window geometry           */
		    window_t *windat,   /* Window data               */
		    win_priv_t *wincon, /* Window constants          */
		    open_bdrys_t *open, /* Open boundary             */
		    double *eta         /* eta variable              */
		    )
{
  int n, m, c, cc;

  if (open->stagger & OUTFACE) {
    for (cc = 1; cc <= open->no2_e1; cc++) {
      c = open->obc_e1[cc];
      eta[c] = eta[open->obc_t[cc]];
      for (m = 0; m < open->bgz; m++) {
	c = open->omap[c];
	eta[c] = eta[open->obc_t[cc]];
      }
    }
    for (cc = 1; cc <= open->no2_e2; cc++) {
      c = open->obc_e2[cc];
      eta[c] = eta[open->obc_t[cc]];
      for (m = 1; m < open->bgz; m++) {
	c = open->omap[c];
	eta[c] = eta[open->obc_t[cc]];
      }
    }
    /* Set the elevation at auxiliary OBC cells */
    for (cc = 1; cc <= open->no2_a; cc++) {
      c = open->obc_a[cc];
      for (m = 1; m < open->bgz; m++) {
	c = open->omap[c];
	eta[c] = eta[open->obc_t[cc]];
      }
    }
  }
}

/* END reset_bdry_eta()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to reset the maps across all open boundaries              */
/*-------------------------------------------------------------------*/
void set_map_eta(geometry_t *window /* Window data structure */
  )
{
  int n;
  int cc, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries */
  /* to be self-mapping on the open boundary.  */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & OUTFACE) {
      /*if (open[n]->bcond_ele != NOTHIN && open[n]->stagger & OUTFACE) {*/
      if (open[n]->type & U1BDRY && open[n]->ocodex & R_EDGE) {
        for (cc = 1; cc <= open[n]->no2_e1; cc++) {
          c = open[n]->obc_e1[cc];
          window->xp1[window->xm1[c]] = c;
        }
      }
      if (open[n]->type & U2BDRY && open[n]->ocodey & F_EDGE) {
        for (cc = 1; cc <= open[n]->no2_e2; cc++) {
          c = open[n]->obc_e2[cc];
          window->yp1[window->ym1[c]] = c;
        }
      }
    }
  }
}

/* END set_map_eta()                                                  */
/*--------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to reset the maps across all interior open boundaries     */
/*-------------------------------------------------------------------*/
void set_map_inside(geometry_t *window /* Window data structure */
  )
{
  int n;
  int cc, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries */
  /* to be self-mapping on the open boundary.  */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      if (open[n]->type & U1BDRY && open[n]->ocodex & R_EDGE) {
        for (cc = 1; cc <= open[n]->no2_e1; cc++) {
          c = open[n]->obc_e1[cc];
          window->xp1[c] = c;
        }
      }
      if (open[n]->type & U2BDRY && open[n]->ocodey & F_EDGE) {
        for (cc = 1; cc <= open[n]->no2_e2; cc++) {
          c = open[n]->obc_e2[cc];
          window->yp1[c] = c;
        }
      }
    }
  }
}

/* END set_map_inside()                                               */
/*--------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/* Routine to apply the Asselin filter to remove the computational    */
/* mode. This is done for velocity only before the elevation          */
/* calculation so that elevation is computed with the filtered fluxes */
/* and thus the divergence of 3D adjusted velocities is consistent    */
/* with elevation change in the surface layer (conservation is        */
/* maintained for tracers).                                           */
/*--------------------------------------------------------------------*/
void asselin(geometry_t *window,  /* Processing window */
             window_t *windat,  /* Window data structure */
             win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, n;                 /* Sparse coordinate, counter */
  double aconst = 0.1;          /* Time filter constant */
  double mdx, mdy;              /* Total depth at time t */
  double mdxb, mdyb;            /* Total depth at time t-1 */
  double *u1av = wincon->d1;
  double *u2av = wincon->d2;
  int *mask = wincon->i7;
  double yr = 86400.0 * 365.0;

  memcpy(u1av, windat->u1av, window->sgsizS * sizeof(double));
  memcpy(u2av, windat->u2av, window->sgsizS * sizeof(double));
  memset(mask, 0, window->sgsizS * sizeof(int));

  /* For sigma the total depth at time t+1 is not available until */
  /* after eta_step(), hence the time filtering is done in terms of */
  /* velocity*depth instead of velocity. Note that nu1av=vel*depth at */
  /* this stage.  */
  if (wincon->sigma) {
    /* u1 velocity.  */
    for (c = 1; c <= window->enonS; c++) {
      mdx = wincon->mdx[c];
      mdxb = mdxbs(window, windat, wincon, c);
      windat->u1av[c] = (windat->u1av[c] * mdx + 0.5 * aconst *
                         (windat->u1avb[c] * mdxb -
                          2.0 * windat->u1av[c] * mdx +
                          windat->nu1av[c])) / mdx;
    }

    /* u2 velocity.  */
    for (c = 1; c <= window->enonS; c++) {
      mdy = wincon->mdy[c];
      mdyb = mdybs(window, windat, wincon, c);
      windat->u2av[c] = (windat->u2av[c] * mdy + 0.5 * aconst *
                         (windat->u2avb[c] * mdyb -
                          2.0 * windat->u2av[c] * mdy +
                          windat->nu2av[c])) / mdy;
    }
  } else {
    /* u1 velocity.  */
    for (c = 1; c <= window->enonS; c++) {
      windat->u1av[c] = windat->u1av[c] + 0.5 * aconst *
        (windat->u1avb[c] - 2.0 * windat->u1av[c] + windat->nu1av[c]);
    }
    /* u2 velocity.  */
    for (c = 1; c <= window->enonS; c++) {
      windat->u2av[c] = windat->u2av[c] + 0.5 * aconst *
        (windat->u2avb[c] - 2.0 * windat->u2av[c] + windat->nu2av[c]);
    }
  }

  /* If the open boundary condition for 2D velocity is VERTIN, then   */
  /* do not perform time filtering (restore boundary velocities to    */
  /* those before filtering is performed).                            */
  reset_map_t_all(window);
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    scale_details_t *scale = open->sdata_e;

    if (open->bcond_tan2d & (VERTIN|FILEIN|CUSTOM|TIDALC)) {
      for (cc = open->no3_e1+1; cc <= open->to2_e1; cc++) {
	c = open->obc_e1[cc];
	if (!mask[c])
	  windat->u1av[c] = u1av[c];
      }
      for (cc = open->no3_e2+1; cc <= open->to2_e2; cc++) {
	c = open->obc_e2[cc];
	if (!mask[c])
	  windat->u2av[c] = u2av[c];
      }
    }

    if (open->bcond_nor2d & (VERTIN|FILEIN|CUSTOM|TIDALC)) {
      int xp1, yp1, ci;
      double f1, fp1, f2, fp2, nvel, v1, v2;
      double rts, eta, tide, sgn;
      double rtsh = 1.2;

      /*--------------------------------------------------------------*/
      /* U1 boundaries                                                */
      if (open->type == U1BDRY) {
	if (!(open->bcond_nor2d & FLATHR)) {
	  for (cc = 1; cc <= open->no2_e1; cc++) {
	    c = open->obc_e1[cc];
	    ci = (open->ocodex & R_EDGE) ? open->nmap[c] : c;
	    if (!mask[open->nmap[ci]])
	      windat->u1av[c] = u1av[c];
	  }
	}
        /*------------------------------------------------------------*/
	/* Adjust the boundary velocity if required                   */
	if (open->adjust_flux) {
	  double adjust_flux;
	  if (wincon->rampf & FLUX_ADJUST && windat->rampval < 1.0) {
	    double dtr = wincon->rampend - wincon->rampstart;
	    adjust_flux = (windat->t - wincon->rampstart) * 
	      (open->adjust_flux - yr) / dtr + yr;
	  } else
	    adjust_flux = open->adjust_flux;
	  rts = windat->dtb2 / adjust_flux;

	  /* Loop through the U1 boundary cells                       */
	  for (cc = 1; cc <= open->no2_e1; cc++) {
	    c = open->obc_e1[cc];
	    xp1 = open->nmap[c];
	    /* Note : xp1 is the cell interior to c */
	    ci = (open->ocodex & R_EDGE) ? xp1 : c;
	    sgn = (open->ocodex & L_EDGE) ? 1.0 : -1.0;
	    yp1 = window->yp1[ci];
	    if (mask[ci]) continue;

	    /* Get the default flux adjustment */
	    if (open->adjust_flux <= 0) {
	      adjust_flux = window->h1acell[ci] / sqrt(wincon->g * 
						       windat->depth_e1[c] * wincon->Ds[c]);
	      rts = windat->dtb2 / adjust_flux;
	    }
	    /* Hard relaxation on corners with overlap                */
	    if ((open->options & OP_OBCCM) && open->olap[cc]) {
	      adjust_flux = rtsh * windat->dtb2;
	      rts = windat->dtb2 / adjust_flux;
	    }
	    /* Save the flux adjustment diagnostic if required        */
	    if (windat->obc_phase) {
	      if (open->adjust_flux_s)
		windat->obc_phase[ci] = windat->u1av[c];
	      else
		windat->obc_phase[ci] = adjust_flux;
	    }

	    /* Get the fluxes                                         */
	    fp1 = windat->u1av[xp1] * windat->depth_e1[xp1] * 
	      window->h2au1[xp1] * wincon->mdx[xp1] * windat->dt2d;
	    f2 = windat->u2av[ci] * windat->depth_e2[ci] * 
	      window->h1au2[ci] * wincon->mdy[ci] * windat->dt2d;
	    fp2 = windat->u2av[yp1] * windat->depth_e2[yp1] * 
	      window->h1au2[yp1] * wincon->mdy[yp1] * windat->dt2d;	

	    /* Get the boundary flux that would result in a flux      */
	    /* divergence sufficient to result in sea level equal to  */
	    /* eta_rlx[].                                             */
	    eta = 0.0;
	    if (open->bcond_ele & (FILEIN|CUSTOM))
	      eta += (wincon->compatible & V1670) ? windat->eta_rlx->val1[ci] :
		open->transfer_eta[cc];

	    /* Scaling                                                */
	    if (scale->type == (TRSC_SUM | TRSC_NUM)) {
	      eta += scale->fact;
	    } else if (scale->type == (TRSC_SUM | TRSC_TRA)) {
	      int trn = scale->ntr;
	      eta += windat->tr_wcS[trn][ci];
	    } else if (scale->type == (TRSC_PCT | TRSC_NUM)) {
	      eta *= scale->fact;
	    } else if (scale->type == (TRSC_PCT | TRSC_TRA)) {
	      int trn = scale->ntr;
	      eta *= windat->tr_wcS[trn][ci];
	    }

	    /* Tidal elevations                                       */
	    tide = 0.0;
	    if (open->bcond_ele & (TIDALH|TIDALC)) {
	      double gmt = windat->t - wincon->tz;
	      double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? 
		windat->rampval : 1.0;
	      tide = ramp * csr_tide_eval(&open->tc, cc, gmt);
	    }
	    if (open->bcond_ele & TIDEBC) {
	      double ramp = (wincon->rampf & TIDEBC) ? 
		windat->rampval : 1.0;
	      tide = ramp * bc_tidal(window, windat, wincon, open->tideforce, ci);
	    }
	    if (open->inverse_barometer) {
	      double ramp = (wincon->rampf & INV_BARO) ? windat->rampval : 1.0;
	      eta += ramp * (wincon->ambpress - windat->patm[ci]) /
		(windat->dens[ci] * wincon->g);
	    }
	    /* Adaptive flux adjustment
	    adjust_flux = get_flux_adjust(window, windat, wincon, open, 
					  windat->eta[c] - tide - eta, c, ci);
	    rts = windat->dtb2 / adjust_flux;
	    */

	    /* Dual relaxation: Relax hard to the tidal signal */
	    if (open->adjust_flux_s && open->bcond_ele & (TIDALH|TIDALC|TIDEBC)) {
	      double depth, etat, df, etadiff;
	      double rtst = (open->adjust_flux_s < 0.0) ? rts : windat->dtb2 / open->adjust_flux_s;
	      f2 *= 0.5; fp2 *= 0.5; fp1 *= 0.5;
	      df = fp1 + sgn * (fp2 - f2);

	      /* Get the tidal flux                                     */
	      f1 = sgn * ((tide - windat->etab[ci]) * 
			  window->cellarea[ci]) + df;
	      /* Convert this flux to a 2D velocity                     */
	      nvel = v1 = f1 / (windat->depth_e1[c] * window->h2au1[c] * 
				wincon->mdx[c] * windat->dt2d);
	      /* Relax the boundary velocity to this value              */
	      windat->u1av[c] -= rtst * (windat->u1av[c] - nvel);
	      /* Get the flux due to this velocity                      */
	      f1 = windat->u1av[c] * windat->depth_e1[c] * 
		window->h2au1[c] * wincon->mdx[c] * windat->dt2d;
	      /* Get the elevation due to this flux                     */
	      etat = windat->etab[ci]  + ((f1 - fp1)/sgn + (f2 - fp2)) / window->cellarea[ci];

	      /* Get the low frequency flux                             */
	      f1 = sgn * ((eta + tide - etat) * window->cellarea[ci]) + df;
	      depth = etat - window->botzu1[c];
	      /* Convert this flux to a 2D velocity                     */
	      nvel = v2 = f1 / (depth * window->h2au1[c] * 
				wincon->mdx[c] * windat->dt2d);

	      /* Get the velocity required for relaxation               */
	      nvel = (windat->u1av[c] + rts * nvel) / rts;
	    } else {  /* Single relaxation                              */
	      f1 = sgn * ((eta + tide - windat->etab[ci])* 
	      window->cellarea[ci] + fp2 - f2) + fp1;
	      /* Convert this flux to a 2D velocity                     */
	      nvel = f1 / (windat->depth_e1[c] * window->h2au1[c] * 
			   wincon->mdx[c] * windat->dt2d);
	    }

	    /* Relax the boundary velocity to this value                */
	    windat->u1av[c] -= rts * (windat->u1av[c] - nvel);

	    mask[ci] = 1;
	    /*windat->dum1[ci] = -rts * (windat->u1av[c] - nvel);*/

	    /* Save the flux adjustment diagnostic if required        */
	    if (windat->obc_phase && open->adjust_flux_s)
	      windat->obc_phase[ci] = -windat->dtb2 * (windat->obc_phase[ci] - v1 - v2) /
		(windat->u1av[c] - windat->obc_phase[ci]);
	  }
	}
      }

      /*--------------------------------------------------------------*/
      /* U2 boundaries                                                */
      if (open->type == U2BDRY) {
	if (!(open->bcond_nor2d & FLATHR)) {
	  for (cc = 1; cc <= open->no2_e2; cc++) {
	    c = open->obc_e2[cc];
	    ci = (open->ocodex & F_EDGE) ? open->nmap[c] : c;
	    if (!mask[open->nmap[ci]])
	      windat->u2av[c] = u2av[c];
	  }
	}
        /*------------------------------------------------------------*/
	/* Adjust the boundary velocity if required                   */
	if (open->adjust_flux) {
	  double adjust_flux;
	  if (wincon->rampf & FLUX_ADJUST && windat->rampval < 1.0) {
	    double dtr = wincon->rampend - wincon->rampstart;
	    adjust_flux = (windat->t - master->rampstart) * 
	      (open->adjust_flux - yr) / dtr + yr;
	  } else
	    adjust_flux = open->adjust_flux;
	  rts = windat->dtb2 / adjust_flux;

	  /* Loop through the U1 boundary cells                       */
	  for (cc = 1; cc <= open->no2_e2; cc++) {
	    c = open->obc_e2[cc];
	    yp1 = open->nmap[c];
	    ci = (open->ocodey & F_EDGE) ? yp1 : c;
	    sgn = (open->ocodey & B_EDGE) ? 1.0 : -1.0;
	    xp1 = window->xp1[ci];
	    if (mask[ci]) continue;
	    /* Get the default flux adjustment */
	    if (open->adjust_flux <= 0) {
	      adjust_flux = window->h2acell[ci] / sqrt(wincon->g * 
						       windat->depth_e2[c] * wincon->Ds[c]);
	      rts = windat->dtb2 / adjust_flux;
	    }
	    /* Hard relaxation on corners with overlap                */
	    if ((open->options & OP_OBCCM) && open->olap[cc]) {
	      adjust_flux = rtsh * windat->dtb2;
	      rts = windat->dtb2 / adjust_flux;
	    }
	    /* Save the flux adjustment diagnostic if required        */
	    if (windat->obc_phase) {
	      if (open->adjust_flux_s)
		windat->obc_phase[ci] = windat->u1av[c];
	      else
		windat->obc_phase[ci] = adjust_flux;
	    }

	    /* Get the fluxes                                         */
	    f1 = windat->u1av[ci] * windat->depth_e1[ci] * 
	      window->h2au1[ci] * wincon->mdx[ci] * windat->dt2d;
	    fp1 = windat->u1av[xp1] * windat->depth_e1[xp1] * 
	      window->h2au1[xp1] * wincon->mdx[xp1] * windat->dt2d;
	    fp2 = windat->u2av[yp1] * windat->depth_e2[yp1] * 
	      window->h1au2[yp1] * wincon->mdy[yp1] * windat->dt2d;	

	    /* Get the boundary flux that would result in a flux      */
	    /* divergence sufficient to result in sea level equal to  */
	    /* eta_rlx[].                                             */
	    eta = 0.0;
	    if (open->bcond_ele & (FILEIN|CUSTOM))
	      eta += (wincon->compatible & V1670) ? windat->eta_rlx->val1[ci] :
		open->transfer_eta[cc];

	    /* Scaling                                                */
	    if (scale->type == (TRSC_SUM | TRSC_NUM)) {
	      eta += scale->fact;
	    } else if (scale->type == (TRSC_SUM | TRSC_TRA)) {
	      int trn = scale->ntr;
	      eta += windat->tr_wcS[trn][ci];
	    } else if (scale->type == (TRSC_PCT | TRSC_NUM)) {
	      eta *= scale->fact;
	    } else if (scale->type == (TRSC_PCT | TRSC_TRA)) {
	      int trn = scale->ntr;
	      eta *= windat->tr_wcS[trn][ci];
	    }

	    tide = 0.;
	    /* Tidal elevations                                       */
	    tide = 0.0;
	    if (open->bcond_ele & (TIDALH|TIDALC)) {
	      double gmt = windat->t - wincon->tz;
	      double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? 
		windat->rampval : 1.0;
	      tide = ramp * csr_tide_eval(&open->tc, cc, gmt);
	    }
	    if (open->bcond_ele & TIDEBC) {
	      double ramp = (wincon->rampf & TIDEBC) ? 
		windat->rampval : 1.0;
	      tide = ramp * bc_tidal(window, windat, wincon, open->tideforce, ci);
	    }
	    if (open->inverse_barometer) {
	      double ramp = (wincon->rampf & INV_BARO) ? windat->rampval : 1.0;
	      eta += ramp * (wincon->ambpress - windat->patm[ci]) /
		(windat->dens[ci] * wincon->g);
	    }

	    /* Dual relaxation: Relax hard to the tidal signal */
	    if (open->adjust_flux_s && open->bcond_ele & (TIDALH|TIDALC|TIDEBC)) {
	      double depth, etat, df;
	      double rtst = (open->adjust_flux_s < 0.0) ? rts : windat->dtb2 / open->adjust_flux_s;
	      f1 *= 0.5; fp1 *= 0.5; fp2 *= 0.5;
	      df = fp2 + sgn * (fp1 - f1);

	      /* Get the tidal flux                                     */
	      f2 = sgn * ((tide - windat->etab[ci]) * 
			  window->cellarea[ci]) + df;
	      /* Convert this flux to a 2D velocity                     */
	      nvel = v1 = f2 / (windat->depth_e2[c] * window->h1au2[c] * 
				wincon->mdy[c] * windat->dt2d);
	      /* Relax the boundary velocity to this value              */
	      windat->u2av[c] -= rtst * (windat->u2av[c] - nvel);
	      /* Get the flux due to this velocity                      */
	      f2 = windat->u2av[c] * windat->depth_e2[c] * 
		window->h1au2[c] * wincon->mdy[c] * windat->dt2d;
	      /* Get the elevation due to this flux                     */
	      etat = windat->etab[ci]  + (f1 - fp1 + (f2 - fp2)/sgn) / window->cellarea[ci];

	      /* Get the low frequency flux                             */
	      f2 = sgn * ((eta + tide - etat) * window->cellarea[ci]) + df;
	      depth = etat - window->botzu2[c];
	      /* Convert this flux to a 2D velocity                     */
	      nvel = v2 = f2 / (depth * window->h1au2[c] * 
				wincon->mdy[c] * windat->dt2d);
	      /* Get the velocity required for relaxation               */
	      nvel = (windat->u2av[c] + rts * nvel) / rts;
	    } else {  /* Single relaxation                              */
	      f2 = sgn * ((eta + tide - windat->etab[ci])* 
			  window->cellarea[ci] + fp1 - f1) + fp2;
	      /* Convert this flux to a 2D velocity                     */
	      nvel = f2 / (windat->depth_e2[c] * window->h1au2[c] * 
			   wincon->mdy[c] * windat->dt2d);
	    }

	    /* Relax the boundary velocity to this value                */
	    windat->u2av[c] -= rts * (windat->u2av[c] - nvel);
	    mask[ci] = 1;
	    /*windat->dum1[ci] = -rts * (windat->u1av[c] - nvel);*/

	    /* Save the flux adjustment diagnostic if required        */
	    if (windat->obc_phase && open->adjust_flux_s)
	      windat->obc_phase[ci] = -windat->dtb2 * (windat->obc_phase[ci] - v1 - v2) /
		(windat->u2av[c] - windat->obc_phase[ci]);
	  }
	}
      }
    }
  }

  /* Elevation is not filtered since it is calculated using the */
  /* filtered velocities.  */
  /*
  if(((wincon->ic+1)%2) == 0) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc]; 
      windat->eta[c] = windat->eta[c] + 0.5 * aconst *
	(windat->etab[c] - 2.0 * windat->eta[c] + wincon->neweta[c]); 
    } 
  }
  */
}

/* END asselin()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the u1 velocity from the u1 velocity transport */
/*-------------------------------------------------------------------*/
void extract_velocity_2d(geometry_t *window,  /* Processing window */
                         window_t *windat,  /* Window data structure */
                         win_priv_t *wincon /* Window geometry / constants 
                                             */
  )
{
  int c, cc, n;                 /* Local sparse coordinate / counter */

  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    wincon->Hn1[c] = mdxns(window, windat, wincon, c);
    windat->u1av[c] /= wincon->Hn1[c];
  }
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    wincon->Hn2[c] = mdyns(window, windat, wincon, c);
    windat->u2av[c] /= wincon->Hn2[c];
  }
}

/* END extract_velocity_2d()                                         */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* LEAPFROG : Reset the velocities for the new time level.           */
/*-------------------------------------------------------------------*/
void leapfrog_update_2d(geometry_t *window, /* Processing window */
                        window_t *windat, /* Window data structure */
                        win_priv_t *wincon  /* Window geometry / constants 
                                             */
  )
{

  /*-----------------------------------------------------------------*/
  /* Set the velocities for the new time level */
  memcpy(windat->u1avb, windat->u1av, window->sgsizS * sizeof(double));
  memcpy(windat->u1av, windat->nu1av, window->sgsizS * sizeof(double));
  memcpy(windat->u2avb, windat->u2av, window->sgsizS * sizeof(double));
  memcpy(windat->u2av, windat->nu2av, window->sgsizS * sizeof(double));
  memcpy(windat->etab, windat->eta, window->sgsizS * sizeof(double));
  memcpy(windat->eta, wincon->neweta, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Extract the velocities from the updated solution for sigma */
  if (wincon->sigma)
    extract_velocity_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the mean u1av velocity if required */
  if (wincon->means & ETA_M) {
    int cc, c;
    double t = windat->dtf2, to;
    if (windat->etam) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	/* The mean for eta is calculated every time-step            */
        to = fabs(windat->meanc[window->m2d[c]]);
	windat->etam[c] = (windat->etam[c] * windat->meanc[c] + 
			   windat->eta[c] * t) / (windat->meanc[c] + t);
      }
    }
  }
  if (wincon->means & VEL2D) {
    int cc, c;
    double t = windat->dtf2;
    if (windat->u1am) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	windat->u1am[c] = (windat->u1am[c] * windat->meanc[c] + 
			   windat->u1av[c] * t) / (windat->meanc[c] + t);
      }
    }

    if (windat->u2am) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	windat->u2am[c] = (windat->u2am[c] * windat->meanc[c] + 
			   windat->u2av[c] * t) / (windat->meanc[c] + t);
      }
    }
  }
}

/* END leapfrog_update_2d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to calculate the total water depth at e1 and e2 faces.    */
/*-------------------------------------------------------------------*/
void get_depths(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc;                    /* Local sparse coordinate / counter */
  int xm1, ym1;                 /* Sparse location at i-1, j-1 */
  double top;                   /* Surface elevation at the cell face */
  double *tzp;                  /* Surface height array */
  double *depth;                /* Pointer to depth array */
  int n, c1;
  open_bdrys_t **open = window->open;

  /*-----------------------------------------------------------------*/
  /* tzp is a pointer which points to either topz[][] or eta[][].  */
  /* topz[][] only gets updated once every 3d step. For the non- */
  /* linear case, we need to use eta instead, which gets updated */
  /* every 2d step. This value is only used to calculate transport */
  /* (nothing to do with the surface slope term, which always uses */
  /* eta).  */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Total depth at e1 faces.  */
  /* Set the pointers.  */
  depth = windat->depth_e1;

  /*-----------------------------------------------------------------*/
  /* Loop over the u1 cells. Note the depth on the boundary is set */
  /* in bdry_u1_2d().  */
  for (cc = 1; cc <= window->v2_e1; cc++) {
    c = window->w2_e1[cc];
    xm1 = window->xm1[c]; 
    top = max(tzp[xm1], tzp[c]);
    depth[c] = top - window->botzu1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the total depth on the boundary */
  for (n = 0; n < window->nobc; n++) {
    for (cc = 1; cc <= open[n]->no2_e1; cc++) {
      c = open[n]->obc_e1[cc];
      c1 = open[n]->oi1_e1[cc];
      depth[c] = tzp[c] - window->botzu1[c];
    }
  }
  if (wincon->dozoom != NOZOOM) {
    for (n = 0; n < window->nobc; n++) {
      for (cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	xm1 = window->xm1[c]; 
	top = max(tzp[xm1], tzp[c]);
	depth[c] = top - window->botzu1[c];
      }
    }
  }

  /* Decrease the depth for sub-surface explicit maps if required    */
  if (window->sm_e1) {
    int cs;
    for (cc = 1; cc <= window->b2_e1; cc++) {
      cs = window->w2_e1[cc];
      c = window->sm_e1[cs];
      if (c)
	/*depth[cs] += window->gridz[window->zp1[c]];*/
      depth[cs] = window->gridz[window->zp1[c]] - window->botzu1[cs];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Total depth at e2 faces.  */
  /* Set the pointers.  */
  depth = windat->depth_e2;

  /*-----------------------------------------------------------------*/
  /* Loop over the u2 cells */
  for (cc = 1; cc <= window->v2_e2; cc++) {
    c = window->w2_e2[cc];
    ym1 = window->ym1[c];
    top = max(tzp[ym1], tzp[c]);
    depth[c] = top - window->botzu2[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the total depth on the boundary */
  for (n = 0; n < window->nobc; n++) {
    for (cc = 1; cc <= open[n]->no2_e2; cc++) {
      c = open[n]->obc_e2[cc];
      c1 = open[n]->oi1_e2[cc];
      depth[c] = tzp[c] - window->botzu2[c];
    }
  }
  if (wincon->dozoom != NOZOOM) {
    for (n = 0; n < window->nobc; n++) {
      for (cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	xm1 = window->xm1[c]; 
	top = max(tzp[xm1], tzp[c]);
	depth[c] = top - window->botzu2[c];
      }
    }
  }

  /* Decrease the depth for sub-surface explicit maps if required    */
  if (window->sm_e2) {
    int cs;
    for (cc = 1; cc <= window->b2_e2; cc++) {
      cs = window->w2_e2[cc];
      c = window->sm_e2[cs];
      if (c)
	/*depth[cs] += window->gridz[window->zp1[c]];*/
	depth[cs] = window->gridz[window->zp1[c]] - window->botzu2[cs];
    }
  }
}

/* END get_depths()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_tr_2d(geometry_t *window,  /* Window structure */
                       window_t *windat,  /* Window data structure */
                       double trem, /* Time remaining in the sub-step loop 
                                     */
                       double *vel, /* Velocity array */
                       double *as,  /* Multi-dt cells velocity time t
                                       values */
                       double *ae /* Multi-dt cells velocity time t+1
                                     values */
  )
{
  int cc, c;                    /* Sparse counter / coordinate */

  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    if (c < window->enonS)
      vel[c] =
        (windat->dt - trem) * (ae[cc] - as[cc]) / window->taux_t[cc] +
        as[cc];
  }
}

/* END set_multidt_tr_2d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the elevation for the 3D mode. Follows POM method. */
/*-------------------------------------------------------------------*/
void get_mode3d_eta(geometry_t *window, /* Processing window */
                    window_t *windat, /* Window data structure */
                    win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int ic = wincon->ic;          /* 2D step counter */
  int iratio = windat->iratio;  /* 3D/2D ratio */
  int c, cc;                    /* Sparse coordinate / counter */
  double aconst = 0.1;          /* Time filter constant */
  double eta3;                  /* 3D mode elevation */

  if (ic < iratio - 3)
    return;
  else if (ic == iratio - 3) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      eta3 = 0.25 * aconst * windat->eta[c];
    }
  } else if (ic == iratio - 2) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      eta3 += 0.5 * (1.0 - 0.5 * aconst) * windat->eta[c];
    }
  } else if (ic == iratio - 1) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      eta3 += 0.5 * windat->eta[c];
    }
  }
}

/* END get_mode3d_eta()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to vertically integrate 3D velocity                       */
/*-------------------------------------------------------------------*/
void vint_3d(geometry_t *window,   /* Window geometry          */
	    window_t *windat,     /* Window data              */
	    win_priv_t *wincon)   /* Window constants         */
{
  int cc, c, cs, zm1;
  double depth;

  memset(windat->u1av, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    depth = 0.0;
    while(c != zm1) {
      if (!isnan(windat->u1[c])) {
	windat->u1av[cs] += windat->u1[c] * windat->dzu1[c];
	depth += windat->dzu1[c];
      }
      c = zm1;
      zm1 = window->zm1[zm1];
    }
    windat->u1av[cs] /= depth;
  }
  memset(windat->u2av, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    depth = 0.0;
    while(c != zm1) {
      if (!isnan(windat->u2[c])) {
	windat->u2av[cs] += windat->u2[c] * windat->dzu2[c];
	depth += windat->dzu2[c];
      }
      c = zm1;
      zm1 = window->zm1[zm1];
    }
    windat->u2av[cs] /= depth;

  }
}

/* END vint_3d()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check if the model is going unstable and exit if so    */
/*-------------------------------------------------------------------*/
int check_unstable(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   int mode)             /* Instability mode         */
{
  double maxwind = 100.0;   /* Maximum wind allowed (ms-1)           */
  int c, cc;

  /* Check for sea level fatalities                                  */
  if (wincon->fatal & ETA_A && mode & ETA_A) {

    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];

      if (fabs(wincon->neweta[c]) > wincon->etamax) {
	hd_quit_and_dump
	  ("etastep: Surface exceeded ETAMAX (%f) at (%d %d) t=%8.3f days\n",
	   wincon->neweta[c], window->s2i[c], window->s2j[c],
	   windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(wincon->neweta[c])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("etastep: NaN at %d (%d %d) t=%8.3f days\n", c,
			     window->s2i[c], window->s2j[c], windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("etastep: NaN at %d (%ld %d) t=%8.3f days : window %d\n",
			     c, window->s2i[c], window->s2j[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
  }

  /* Check for 2D velocity fatalities                                */
  if (wincon->fatal & VEL2D && mode & VEL2D) {
    for (cc = 1; cc <= window->b2_e1; cc++) {
      c = window->w2_e1[cc];

      if (fabs(windat->u1av[c]) > wincon->velmax2d) {
	hd_quit_and_dump
	  ("vel2d: u1av velocity exceeded velmax (%f) at (%d %d) t=%8.3f days\n",
	   windat->u1av[c], window->s2i[c], window->s2j[c], windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u1av[c])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("vel2d: u1av NaN at (%d %d) t=%8.3f days\n",
			     window->s2i[c], window->s2j[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("vel2d: u1av NaN at (%d %d) t=%8.3f days : window %d\n",
			     window->s2i[c], window->s2j[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
    for (cc = 1; cc <= window->b2_e2; cc++) {
      c = window->w2_e2[cc];

      if (fabs(windat->u2av[c]) > wincon->velmax2d) {
	hd_quit_and_dump
	  ("vel2d: u2av velocity exceeded velmax (%f) at (%d %d) t=%8.3f days\n",
	   windat->u2av[c], window->s2i[c], window->s2j[c], windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u2av[c])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("vel2d: u2av NaN at (%d %d) t=%8.3f days\n",
			     window->s2i[c], window->s2j[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("vel2d: u2av NaN at (%d %d) t=%8.3f days : window %d\n",
			     window->s2i[c], window->s2j[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
  }

  /* Check for 3D velocity fatalities                                */
  if (wincon->fatal & VEL3D && mode & VEL3D) {
    for (cc = 1; cc <= window->b3_e1; cc++) {
      c = window->w3_e1[cc];

      if (fabs(windat->u1[c]) > wincon->velmax) {
	hd_quit_and_dump
	  ("vel3d: u1 velocity exceeded velmax (%f) at (%d %d %d) t=%8.3f days\n",
	   windat->u1[c], window->s2i[c], window->s2j[c], window->s2k[c], 
	   windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u1[c])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("vel3d: u1 NaN at (%d %d %d) t=%8.3f days\n",
			     window->s2i[c], window->s2j[c], window->s2k[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("vel3d: u1 NaN at (%d %d %d) t=%8.3f days : window %d\n",
			     window->s2i[c], window->s2j[c], window->s2k[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
    for (cc = 1; cc <= window->b3_e2; cc++) {
      c = window->w3_e2[cc];

      if (fabs(windat->u2[c]) > wincon->velmax) {
	hd_quit_and_dump
	  ("vel3d: u2 velocity exceeded velmax (%f) at (%d %d %d) t=%8.3f days\n",
	   windat->u2[c], window->s2i[c], window->s2j[c], window->s2k[c], 
	   windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u2[c])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("vel3d: u2 NaN at (%d %d %d) t=%8.3f days\n",
			     window->s2i[c], window->s2j[c], window->s2k[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("vel3d: u2 NaN at (%d %d %d) t=%8.3f days : window %d\n",
			     window->s2i[c], window->s2j[c], window->s2k[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
  }

  /* Check for wind stress fatalities                                */
  if (wincon->fatal & WIND && mode & WIND) {
    double maxstr = 0.00218 * maxwind * maxwind * air_dens;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];

      if (fabs(windat->wind1[c]) > maxstr) {
	hd_quit_and_dump
	  ("wind1: Wind exceeded 100ms-1 (%f) at (%d %d) t=%8.3f days\n",
	   windat->wind1[c], window->s2i[c], window->s2j[c],
	   windat->t / 86400);
	return(1);
      }
      if (fabs(windat->wind2[c]) > maxstr) {
	hd_quit_and_dump
	  ("wind2: Wind exceeded 100ms-1 (%f) at (%d %d) t=%8.3f days\n",
	   windat->wind2[c], window->s2i[c], window->s2j[c],
	   windat->t / 86400);
	return(1);
      }
    }
  }

  /* Check for temp & sal fatalities                                 */
  if (wincon->fatal & TS && mode & TS && wincon->fatal & NANF) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];

      if (isnan(windat->temp[c])) {
	hd_quit_and_dump
	  ("temp: NaN found  at (%d %d %d) t=%8.3f days : window %d\n",
	   window->s2i[c], window->s2j[c], window->s2k[c],
	   windat->t / 86400, window->wn);
	return(1);
      }
      if (isnan(windat->sal[c])) {
	hd_quit_and_dump
	  ("salt: NaN found  at (%d %d %d) t=%8.3f days : window %d\n",
	   window->s2i[c], window->s2j[c], window->s2k[c],
	   windat->t / 86400, window->wn);
	return(1);
      }
    }
  }
  return(0);    
}

/* END check_unstable()                                              */
/*-------------------------------------------------------------------*/


