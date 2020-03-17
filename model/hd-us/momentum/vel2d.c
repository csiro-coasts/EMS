/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/momentum/vel2d.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vel2d.c 6471 2020-02-18 23:50:34Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#define RADIUS 6370997.0

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* mode2d_step()      : Updates the 2D mode variables                */
/* vel_u1av_update()  : Solves the 2D u1 velocity equation           */
/* advect_u1_2d()     : Solves the u1 advection equation             */
/* bdry_u1_2d()       : Sets normal and tangential u1 OBC's          */
/* set_map_e1()       : Sets maps self pointing across u1 boundaries */
/* extract_velocity_2d() : Gets u1&u2 velocity from u1&u2 transport  */
/* eta_step()         : Calculates free surface height               */
/* bdry_eta()         : Sets the free surface boundary condition     */
/* vel2D_lbc()        : Sets 2D velocity lateral boundary conditions */
/* mdxs()             : Returns total depth at e1 faces at time t    */
/* mdxns()            : Returns total depth at e1 faces at time t+1  */
/* mdys()             : Returns total depth at e2 faces at time t    */
/* mdyns()            : Returns total depth at e2 faces at time t+1  */
/*-------------------------------------------------------------------*/


void set_viscosity_2d(geometry_t *window);
void precalc_u1_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void vel_tan_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void vel_components_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
double est_bot_stress(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      int e, double botu1, double val, double pgt, double cot, double uin);
double tpxo_error(geometry_t *window, window_t *windat, win_priv_t *wincon);
void pressure_u1av(geometry_t *window, window_t *windat, win_priv_t *wincon);
void coriolis_u1av(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bottom_u1av(geometry_t *window, window_t *windat, win_priv_t *wincon);
void nonlinear_u1av(geometry_t *window, window_t *windat, win_priv_t *wincon);

/*-------------------------------------------------------------------*/
/* 2D window step part 1                                             */
/*-------------------------------------------------------------------*/
void mode2d_step_window_p1(master_t *master,
			   geometry_t *window,
			   window_t *windat, win_priv_t *wincon)
{

  /*-----------------------------------------------------------------*/
  /* Get the cell centred east and north velocities                  */
  vel_cen(window, windat, wincon, windat->u1av, windat->u2av, 
	  windat->uav, windat->vav, NULL, NULL, 1);

  /*-----------------------------------------------------------------*/
  /* Update the 2D velocity                                          */
  if (wincon->compatible & V6257)
    vel_u1av_update_seq(window, windat, wincon);
  else
    vel_u1av_update(window, windat, wincon);

  /* Set sources of momentum                                         */
  ss_momentum(window, windat, wincon, VEL2D);

  /*-----------------------------------------------------------------*/
  /* Refill the master with updated velocity from the window         */
  /* data structure. This is required at auxiliary cells for the     */
  /* time filtering of u1av asselin(). This filtered velocity is     */
  /* then used to get fluxes to update eta.                          */
  if (master->nwindows > 1)
    win_data_empty_2d(master, window, windat, NVELOCITY);
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
  double clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Fill the window with updated velocities.                        */
  if (master->nwindows > 1)
    win_data_refill_2d(master, window, windat, master->nwindows,
                       NVELOCITY);
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

  /* Calculate the velocity components faces                         */
  /*vel_components_2d(window, windat, wincon);*/

  /* Get the tangential velocity to the edge                         */
  vel_tan_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Refill the master with elevation and velocity data from the     */
  /* window data structure.                                          */
  if (master->nwindows > 1)
    win_data_empty_2d(master, window, windat, VELOCITY);

  windat->wclk += (dp_clock() - clock);
}

/* END mode2d_step_window_p2()                                       */
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
  int n, nn;                          /* Counters                    */
  int ic;                             /* 2D mode counter             */
  double oldtime = master->t;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */

#if defined(HAVE_OMP)
#pragma omp parallel for private(n)
#endif

  for (n = 1; n <= nwindows; n++) {
    memset(windat[n]->u1flux, 0, window[n]->szeS * sizeof(double));
    memcpy(wincon[n]->oldeta, windat[n]->eta, window[n]->szcS * sizeof(double));
    wincon[n]->neweta = wincon[n]->d4;
    memset(wincon[n]->neweta, 0, window[n]->szcS * sizeof(double));
    if (!(wincon[n]->etarlx & NONE))
      memset(wincon[n]->eta_rlx3d, 0, window[n]->szcS * sizeof(double));
    set_viscosity_2d(window[n]);
  }

  /*-----------------------------------------------------------------*/
  /* Loop to step the 2D mode iratio times                           */
  for (ic = 0; ic < windat[1]->iratio; ic++) {

    /*---------------------------------------------------------------*/
    /* Do the custom elevation routines on the master                */
    master->t3d = oldtime + (double)(ic+1) * master->dt2d;
    bdry_eval_eta_m(geom, master);
    bdry_eval_u1av_m(geom, master);
    if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions using the master          */
#if GLOB_BC
    vel2D_lbc(master->u1av, geom->nbpte1S, geom->nbe1S,
              geom->bpte1S, geom->bine1S, master->slip);
    vel2D_lbc(master->u1bot, geom->nbpte1S, geom->nbe1S,
              geom->bpte1S, geom->bine1S, master->slip);
    set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                       geom->bin2, 1);
    /* Backward velocities are required for stress tensors. Ghost    */
    /* cells are not set accurately since a lateral BC is not set on */
    /* nu1av before the asselin() filtering; setting BCs here        */
    /* overcomes this.                                               */
    if (!(master->compatible & V1283)) {
      vel2D_lbc(master->u1avb, geom->nbpte1S, geom->nbe1S,
		geom->bpte1S, geom->bine1S, master->slip);
    }
#endif

#if defined(HAVE_OMP)
#pragma omp parallel for private(n,nn)
#endif

    for (nn = 1; nn <= nwindows; nn++) {
      n = wincon[1]->twin[nn];
      wincon[n]->ic = ic;

      /*-------------------------------------------------------------*/
      /* Fill 2D variables into the window data structures.          */
      win_data_fill_2d(master, window[n], windat[n], master->nwindows);

      /*-------------------------------------------------------------*/
      /* Set the lateral boundary conditions                         */
#if !GLOB_BC
      vel2D_lbc(windat[n]->u1av, window[n]->nbpte1S, window[n]->nbe1S,
                window[n]->bpte1S, window[n]->bine1S, wincon[n]->slip);
      vel2D_lbc(windat[n]->u1bot, window[n]->nbpte1S, window[n]->nbe1S,
                window[n]->bpte1S, window[n]->bine1S, wincon[n]->slip);
      set_lateral_bc_eta(windat[n]->eta, window[n]->nbptS, window[n]->bpt,
                         window[n]->bin, window[n]->bin2, 1);
      if (!(master->compatible & V1283)) {
	vel2D_lbc(windat[n]->u1avb, window[n]->nbpte1S, window[n]->nbe1S,
		  window[n]->bpte1S, window[n]->bine1S, wincon[n]->slip);
      }
#endif

      /*-------------------------------------------------------------*/
      /* Set auxiliary cells that have been updated to the multi dt  */
      /* buffers.                                                    */
      fill_multidt_aux_2d(window[n], windat[n], master->nu1av,
                          windat[n]->u1av_ae);

    }

#if TR_CK
    check_transfers(geom, window, windat, wincon, nwindows, VEL2D);
#endif

    /* Do the 1st part of the 2D step : velocity                     */
    dp_vel2d_step_p1();

    /* Do the 2nd part of the 2D step : elevation and fluxes         */
    dp_vel2d_step_p2();
    if (master->crf == RS_RESTART) return;

#if TR_CK
    check_transfers(geom, window, windat, wincon, nwindows, VEL2D|FLUX);
#endif

#if defined(HAVE_OMP)
#pragma omp parallel for private(n,nn)
#endif

    /* Do the 3rd part of the 2D step : depth at e1 and e2 faces     */
    for (nn = 1; nn <= nwindows; nn++) {
      n = wincon[1]->twin[nn];

      /* Set the master time here so that boundary data is read in   */
      /* on the 2D time-step. This is over-written below to remove   */
      /* precision issues.                                           */
#pragma omp critical
      {
	/* The following two lines must execute atomically */
	master->t    = windat[n]->t;
	master->days = master->t / 86400.0;
      }
      windat[n]->days = master->days;

      if (master->nwindows > 1)
	win_data_refill_2d(master, window[n], windat[n], master->nwindows,
			   ETA_A);

      get_depths(window[n], windat[n], wincon[n]);
      if (master->nwindows > 1)
	win_data_empty_2d(master, window[n], windat[n], DEPTH);
    }

    /* Couple at the barotropic level for 2-way nesting              */
    if (master->obcf & DF_BARO) {
      dump_eta_snapshot(master, window, windat, wincon);
    }
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
/* Update the 2D velocities                                          */
/*-------------------------------------------------------------------*/
void vel_u1av_update(geometry_t *window,  /* Window geometry         */
                     window_t *windat,    /* Window data             */
                     win_priv_t *wincon   /* Window constants        */
		     )
{
  int e, ee;                    /* Edge coordinate / counter         */
  int c1, c2, n;                /* Cell coordinate / counter         */
  double val;                   /* Dummy                             */
  double midx;                  /* Water depth at the edge           */
  double *tzp;                  /* Surface height array              */
  int dbc;                      /* Debugging coordinate              */
  int dbj;                      /* Debugging edge direction          */
  double *depth = windat->depth_e1;

  dbc = (wincon->dbc) ? window->m2d[wincon->dbc] : 0;
  dbj = wincon->dbj;

  /* tzp is a pointer which points to either topz[][] or eta[][].  */
  /* topz[][] only gets updated once every 3d step. For the non- */
  /* linear case, we need to use eta instead, which gets updated */
  /* every 2d step. This value is only used to calculate transport */
  /* (nothing to do with the surface slope term, which always uses */
  /* eta).  */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation          */
  precalc_u1_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection                                     */
  if (wincon->momsc2d & RINGLER) {
    if (nonlin_coriolis_2d(window, windat, wincon)) return;
  } else {
    if (!(wincon->u1av_f & ADVECT))
      advect_u1_2d(window, windat, wincon);
  }

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included)                 */
  wincon->hor_mix->setup_av(window, windat, wincon, tzp);
  wincon->hor_mix->u1av(window, windat, wincon, tzp);

  /*-----------------------------------------------------------------*/
  /* Get the horizontal pressure gradient                            */
  pressure_u1av(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the Coriolis force                                          */
  coriolis_u1av(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the bottom stress                                           */
  bottom_u1av(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the non-linear terms and radiation stresses                 */
  nonlinear_u1av(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Sum the tendencies                                              */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s3[ee];
    /* Note: T_ADV is not added, as the velocity is already updated  */
    /* with this tendedcy in momentum advection to account for       */
    /* sub-stepping.                                                 */
    for (n = 1; n < TEND2D; n++) {
      double dt = (n == T_HDF || n == T_ADV) ? 1.0 : windat->dt2d;
      midx = (n == T_BTP || n == T_COR) ? wincon->mdx[e] : 1.0;
      windat->nu1av[e] += dt * midx * wincon->tend2d[n][e];
    }
    windat->nu1av[e] += windat->dt2d * wincon->u1inter[e];
  }

  /*-----------------------------------------------------------------*/
  /* Restrict flow if column nearly dry                              */
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    if (windat->nu1av[e] > 0.0 &&
        (val = (tzp[c2] - window->botz[c2])) < wincon->hmin)
      windat->nu1av[e] *= max(val, 0.0) / wincon->hmin;
    else if (windat->nu1av[e] < 0.0 &&
             (val = (tzp[c1] - window->botz[c1])) < wincon->hmin)
      windat->nu1av[e] *= max(val, 0.0) / wincon->hmin;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate u1av boundary values                                  */
  bdry_u1_2d(window, windat, wincon);
  debug_c(window, D_UA, D_BDRY);

  /*-----------------------------------------------------------------*/
  /* Debugging                                                       */
  if (dbc) {
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      if (window->e2e[e][0] == dbj && window->e2c[e][dbj] == dbc) {
	wincon->b1 = wincon->tend2d[T_BTP][e];
	wincon->b2 = wincon->tend2d[T_COR][e];
	wincon->b3 = wincon->tend2d[T_BOT][e];
      }
    }
  }
  if (wincon->mode2d) debug_c(window, D_INIT, D_POST);
  debug_c(window, D_UA, D_POST);
  /* Get the bottom stress from the 2D mode if required              */
  if (wincon->mode2d && wincon->numbers & BOTSTRESS) {
    memcpy(wincon->w8, wincon->tend2d[T_BOT], window->szeS * sizeof(double));
    vel_cen(window, windat, wincon, wincon->w8, NULL, windat->tau_be1, windat->tau_be2,
	    windat->tau_bm, NULL, 1);
  }
  /*tpxo_error(window, windat, wincon);*/
}

/* END vel_u1av_update()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2D pressure tendency                                              */
/*-------------------------------------------------------------------*/
void pressure_u1av(geometry_t *window,    /* Window geometry         */
		   window_t *windat,      /* Window data             */
		   win_priv_t *wincon     /* Window constants        */
		   )
{
  int e, ee;
  int c1, c2;
  double alpha = 1.0;           /* Tide SAL constant                 */
  double beta = 0.0;            /* Body tide effect constant         */

  if (wincon->u1av_f & PRESS_BT) return;

  /* Set the tidal potential term if required                        */
  if (wincon->tidep) {
    equ_tide_eval(window, windat, wincon, windat->equitide);
    alpha = wincon->eqt_alpha;
    beta = wincon->eqt_beta;
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];

      /* Equilibrium tide body force                                 */
      wincon->tend2d[T_BTP][e] -= beta * (windat->equitide[c1] - windat->equitide[c2]) * 
	wincon->topdensu1[e] / wincon->densavu1[e];
    }
  }

  /* Pressure gradient term.                                         */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s3[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];

    wincon->tend2d[T_BTP][e] -= (alpha * (windat->eta[c1] - windat->eta[c2]) * 
				 wincon->topdensu1[e] +
				 (windat->patm[c1] - windat->patm[c2])) / 
      wincon->densavu1[e];
  }
}

/* END pressure_u1av()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2D coriolis tendency                                              */
/*-------------------------------------------------------------------*/
void coriolis_u1av(geometry_t *window,    /* Window geometry         */
		   window_t *windat,      /* Window data             */
		   win_priv_t *wincon     /* Window constants        */
		   )
{
  int e, ee;

  if (wincon->u1av_f & CORIOLIS || wincon->momsc & RINGLER) return;

  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s3[ee];
    wincon->tend2d[T_COR][e] = wincon->u1c5[e] * windat->u2av[e];
  }
}

/* END coriolis_u1av()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2D bottom friction tendency                                       */
/*-------------------------------------------------------------------*/
void bottom_u1av(geometry_t *window,      /* Window geometry         */
		 window_t *windat,        /* Window data             */
		 win_priv_t *wincon       /* Window constants        */
		 )
{
  int n, e, ee, eoe;
  int c1, c2;
  double val;                   /* Dummy                             */
  double *depth;                /* Depth of the water column         */
  double midx;                  /* Water depth at the edge           */
  double u2au1;                 /* Tangential velocity               */
  double botu1;                 /* Normal bottom velocity            */
  double botu2;                 /* Tangential Bottom velocity        */
  double Cdu1;                  /* Bottom drag coeff at the edge     */

  depth = windat->depth_e1;

  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s3[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    midx = wincon->mdx[e];

    u2au1 = botu2 = 0;
    for (n = 1; n <= window->nee[e]; n++) {
      eoe = window->eSe[n][e];
      if (!eoe) continue;
      u2au1 += window->wAe[n][e] * windat->u1avb[eoe];
      botu2 += window->wAe[n][e] * windat->u1bot[eoe];
    }
    botu1 = windat->u1avb[e] + windat->u1bot[e];
    botu2 = u2au1 + botu2;

    Cdu1 = 0.5 * (wincon->Cd[c1] + wincon->Cd[c2]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu1 * max(wincon->uf, val);
    /* Truncate to ensure stability                                  */
    if (val > depth[e] * midx / windat->dt2d)
      val = depth[e] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case                               */
    wincon->tend2d[T_BOT][e] = -val * botu1 / depth[e];
  }
}

/* END bottom_u1av()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* 2D non-linear tendency. Radiation stresses are also added if      */
/* required.                                                         */
/*-------------------------------------------------------------------*/
void nonlinear_u1av(geometry_t *window,   /* Window geometry         */
		    window_t *windat,     /* Window data             */
		    win_priv_t *wincon    /* Window constants        */
		    )
{
  int n, e, ee, eoe;
  int c1, c2;
  double *depth;                /* Depth of the water column         */
  double rst = 0.0;             /* Radiation stress term             */
  double rho0 = 1024.0;         /* Reference density                 */

  depth = windat->depth_e1;

  /*-----------------------------------------------------------------*/
  /* Add the non-linear terms                                        */
  /* SIGMA : No extra terms to include for the sigma case.           */
  if (wincon->nonlinear && !wincon->sigma) {
    for (ee = 1; ee <= window->v2_e1; ee++) {
      e = window->w2_e1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];

      wincon->tend2d[T_NLI][e] -= windat->u1av[e] * (windat->detadt[c1] + windat->detadt[c2]) /
        (2.0 * max(depth[e], wincon->hmin));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add the radiation stresses (radiation stresses are cell         */
  /* centered).                                                      */
  if (wincon->waves & WAVE_FOR) {
    double fe;
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      fe = vel_c2e(window, windat->wave_Fx, windat->wave_Fy, e);
      rst = windat->dt2d * fe / (rho0 * max(depth[e], wincon->hmin));
      if (windat->u1_rad) windat->u1_rad[e] = rst;
      wincon->tend2d[T_NLI][e] += rst;
    }
  }
  /* Not implemented for unstructured
  else if (wincon->waves & TAN_RAD) {
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      rst = windat->dt2d * (windat->wave_Sxy[c1] - windat->wave_Sxy[c2]) / 
	window->h2au1[e];
      if (windat->u1_rad) windat->u1_rad[e] = rst;
      wincon->tend2d[T_NLI][e] += rst;
      }
    }
  */
}

/* END nonlinear_u1av()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Update the 2D velocities seqentially (without using tendencies)   */
/*-------------------------------------------------------------------*/
void vel_u1av_update_seq(geometry_t *window,  /* Window geometry     */
			 window_t *windat,    /* Window data         */
			 win_priv_t *wincon   /* Window constants    */
			 )
{
  int e, ee, ep, em, eoe;       /* Edge coordinate / counter         */
  int c1, c2;                   /* Cell coordinate / counter         */
  int n, j, npe;                /* Counters                          */
  double pgt = 0.0;             /* Pressure gradient term            */
  double cot = 0.0;             /* Coriolis term                     */
  double bft;                   /* Bottom friction term              */
  double rst = 0.0;             /* Radiation stress term             */
  double *tzp;                  /* Surface height array              */
  double *depth;                /* Depth of the water column         */
  double midx;                  /* Water depth at the edge           */
  double val;                   /* Dummy                             */
  double u2au1;                 /* Tangential velocity               */
  double botu1;                 /* Normal bottom velocity            */
  double botu2;                 /* Tangential Bottom velocity        */
  double Cdu1;                  /* Bottom drag coeff at the edge     */
  double rho0 = 1024.0;         /* Reference density                 */
  int dbc;                      /* Debugging coordinate              */
  int dbj;                      /* Debugging edge direction          */
  double alpha = 1.0;           /* Tide SAL constant                 */
  double beta = 0.0;            /* Body tide effect constant         */

  dbc = (wincon->dbc) ? window->m2d[wincon->dbc] : 0;
  dbj = wincon->dbj;
  depth = windat->depth_e1;

  /* tzp is a pointer which points to either topz[][] or eta[][].  */
  /* topz[][] only gets updated once every 3d step. For the non- */
  /* linear case, we need to use eta instead, which gets updated */
  /* every 2d step. This value is only used to calculate transport */
  /* (nothing to do with the surface slope term, which always uses */
  /* eta).  */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation          */
  precalc_u1_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection                                     */
  if (wincon->momsc2d & RINGLER) {
    if (nonlin_coriolis_2d(window, windat, wincon)) return;
  } else {
    if (!(wincon->u1av_f & ADVECT))
      advect_u1_2d(window, windat, wincon);
  }

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included)                 */
  wincon->hor_mix->setup_av(window, windat, wincon, tzp);
  wincon->hor_mix->u1av(window, windat, wincon, tzp);

  /*-----------------------------------------------------------------*/
  /* Set the tidal potential term if required                        */
  if (wincon->tidep) {
    alpha = wincon->eqt_alpha;
    beta = wincon->eqt_beta;
    equ_tide_eval(window, windat, wincon, windat->equitide);
  }

  /*-----------------------------------------------------------------*/
  /* Add the body forces and step forward in time                    */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s3[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    midx = wincon->mdx[e];

    /*---------------------------------------------------------------*/
    /* Pressure gradient term. Note : differences bracketed to       */
    /* preserve accuracy.  */
    if (!(wincon->u1av_f & PRESS_BT)) {
      pgt = -(alpha * (windat->eta[c1] - windat->eta[c2]) * 
	      wincon->topdensu1[e] +
	      (windat->patm[c1] - windat->patm[c2])) / 
	wincon->densavu1[e];

      if (wincon->tidep) {
	pgt -= beta * (windat->equitide[c1] - windat->equitide[c2]) * 
	  wincon->topdensu1[e] / wincon->densavu1[e];
      }
    }

    /*---------------------------------------------------------------*/
    /* Coriolis                                                      */
    if (!(wincon->u1av_f & CORIOLIS) && !(wincon->momsc & RINGLER))
      cot = wincon->u1c5[e] * windat->u2av[e];

    /*---------------------------------------------------------------*/
    /* Bottom friction                                               */
    u2au1 = botu2 = 0;
    for (n = 1; n <= window->nee[e]; n++) {
      eoe = window->eSe[n][e];
      if (!eoe) continue;
      u2au1 += window->wAe[n][e] * windat->u1avb[eoe];
      botu2 += window->wAe[n][e] * windat->u1bot[eoe];
    }
    botu1 = windat->u1avb[e] + windat->u1bot[e];
    botu2 = u2au1 + botu2;

    Cdu1 = 0.5 * (wincon->Cd[c1] + wincon->Cd[c2]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu1 * max(wincon->uf, val);
    /* Truncate to ensure stability                                  */
    if (val > depth[e] * midx / windat->dt2d)
      val = depth[e] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case                               */
    bft = -val * botu1 / depth[e];
    if (window->e2e[e][0] == dbj && window->e2c[e][dbj] == dbc) {
      wincon->b1 = pgt;
      wincon->b2 = cot;
      wincon->b3 = bft;
    }
    if (wincon->mode2d && wincon->numbers & BOTSTRESS) windat->tau_be1[e] = bft;

    /*---------------------------------------------------------------*/
    /* Calculate nu1av value                                         */
    /* SIGMA : multiply new velocity by depth (midx)                 */
    windat->nu1av[e] +=
      windat->dt2d * (midx * (pgt + cot) + bft + wincon->u1inter[e]);
  }

  if (wincon->mode2d) debug_c(window, D_INIT, D_POST);
  debug_c(window, D_UA, D_POST);
  /* Get the bottom stress from the 2D mode if required              */
  if (wincon->mode2d && wincon->numbers & BOTSTRESS) {
    memcpy(wincon->w8, windat->tau_be1, window->szeS * sizeof(double));
    vel_cen(window, windat, wincon, wincon->w8, NULL, windat->tau_be1, windat->tau_be2,
	    windat->tau_bm, NULL, 1);
  }

  /*-----------------------------------------------------------------*/
  /* Add the non-linear terms                                        */
  /* SIGMA : No extra terms to include for the sigma case.           */
  if (wincon->nonlinear && !wincon->sigma) {
    for (ee = 1; ee <= window->v2_e1; ee++) {
      e = window->w2_e1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      windat->nu1av[e] -= windat->dt2d * windat->u1av[e] *
        (windat->detadt[c1] + windat->detadt[c2]) /
        (2.0 * max(depth[e], wincon->hmin));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add the radiation stresses (radiation stresses are cell         */
  /* centered).                                                      */
  if (wincon->waves & WAVE_FOR) {
    double fe;
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      fe = vel_c2e(window, windat->wave_Fx, windat->wave_Fy, e);
      rst = windat->dt2d * fe / (rho0 * max(depth[e], wincon->hmin));
      if (windat->u1_rad) windat->u1_rad[e] = rst;
      windat->nu1av[e] += rst;
    }
  }
  /* Not implemented for unstructured
  else if (wincon->waves & TAN_RAD) {
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s3[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      rst = windat->dt2d * (windat->wave_Sxy[c1] - windat->wave_Sxy[c2]) / 
	window->h2au1[e];
      if (windat->u1_rad) windat->u1_rad[e] = rst;
      windat->nu1av[e] += rst;
    }
  }
  */

  /*-----------------------------------------------------------------*/
  /* Restrict flow if column nearly dry                              */
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    if (windat->nu1av[e] > 0.0 &&
        (val = (tzp[c2] - window->botz[c2])) < wincon->hmin)
      windat->nu1av[e] *= max(val, 0.0) / wincon->hmin;
    else if (windat->nu1av[e] < 0.0 &&
             (val = (tzp[c1] - window->botz[c1])) < wincon->hmin)
      windat->nu1av[e] *= max(val, 0.0) / wincon->hmin;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate u1av boundary values                                  */
  bdry_u1_2d(window, windat, wincon);
  debug_c(window, D_UA, D_BDRY);
}

/* END vel_u1av_update_seq()                                         */
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
/*-------------------------------------------------------------------*/
/* Routine to set the correct map across open boundaries for u1      */
/*-------------------------------------------------------------------*/
void set_map_e1(geometry_t *window  /* Window geometry               */
  )
{
  int n;
  int ee, e, es;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries   */
  /* (set in the tracer routine) to be self-mapping on the open      */
  /* boundary.                                                       */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & OUTFACE) {
      for (ee = 1; ee <= open[n]->no3_e1; ee++) {
        e = open[n]->obc_e1[ee];
        es = open[n]->oi1_e1[ee];
	if (open[n]->ceni[ee])
	  window->ep[es] = e;
	else
	  window->em[es] = e;
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
void bdry_u1_2d(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int n;                        /* Counters                          */
  int e, ee, e1;                /* Sparse coordinates / counters     */
  double *depth;                /* Depth of the water column         */
  double midx, midx1;           /* Depth at the e1 face              */
  open_bdrys_t **open = window->open;

  depth = windat->depth_e1;

  /*-----------------------------------------------------------------*/
  /* Set the u1av velocity where u1av is normal to the boundary      */
  for (n = 0; n < window->nobc; n++) {
    /* Get the u1av velocity derived from a no gradient condition    */
    /* on u1flux.                                                    */
    if (open[n]->bcond_nor2d & NOGRAD) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
        e = open[n]->obc_e1[ee];
        e1 = open[n]->oi1_e1[ee];
        midx = wincon->mdx[e];
        midx1 = wincon->mdx[e1];
        /* Set the no gradient condition on velocity */
        windat->nu1av[e] =
          windat->nu1av[e1] * depth[e1] * midx1 * window->h1au1[e1] /
          (max(depth[e], wincon->hmin) * midx * window->h1au1[e]);
      }
    }
    /* Set the u1av normal velocity for all other BC's.              */
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
  /* Set the u1av velocity where u1av is tangential to the boundary  */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	    open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	    open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	    windat->u1av, windat->u1avb, open[n]->bcond_tan2d, 
	    windat->dtf2, &open[n]->datau2av, open[n]->transfer_u2av, 
	    open[n]->relax_zone_tan, U2BDRY);	    
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  /*
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	e1 = open[n]->omape[ee];
	windat->nu1av[e1] = windat->nu1av[e];
      }
    }
  }
  */
  /* Rescale for sigma active OBCs                                   */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->bcond_nor2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1av, mdxns,
			1, open->no2_e1, open->obc_e1);
      if (open->bcond_tan2d & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1av, mdxns,
			open->no2_e1 + 1, open->to2_e1, open->obc_e1);
    }
  }
}

/* END bdry_u1_2d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays for the u1av calculation                           */
/*-------------------------------------------------------------------*/
void precalc_u1_2d(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
		   )
{
  int e, ee;
  int n, eoe, c1, c2;
  int cw = 1;
  double top, bot;

  for(ee = 1; ee <= wincon->vcs1; ee++) {
    e = window->m2de[wincon->s1[ee]];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    top = max(windat->eta[c1], windat->eta[c2]);
    bot = DRY_FRAC * wincon->hmin + window->botzu1[e];
    if (top > bot) {
      wincon->s3[cw] = e;
      cw++;
    }
    else {
      windat->nu1av[e] = windat->u1av[e] = windat->u1avb[e] = 0.0;
    }
  }
  wincon->vcs = cw - 1;
  if(wincon->dolin_u1) {
    linear_bdry_cell(window, windat, wincon, wincon->s3, wincon->vcs,
		     wincon->linmask_u1, NULL);
  }

  /* Copy the current velocities into the update array               */
  memcpy(windat->nu1av, windat->u1avb, window->szeS * sizeof(double));
  for (n = 0; n < TEND2D; n++)
    memset(wincon->tend2d[n], 0, window->szeS * sizeof(double));

  /* Eliminate dry cells from the cells to process list              */
  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations                                              */
  if (wincon->sigma) {
    for (ee = 1; ee <= window->v2_e1; ee++) {
      e = window->w2_e1[ee];
      windat->nu1av[e] *= mdxbs(window, windat, wincon, e);
    }
  }
}

/* END precalc_u1_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute velocity components                            */
/*-------------------------------------------------------------------*/
void vel_components_2d(geometry_t *window,  /* Window geometry       */
		       window_t *windat,    /* Window data           */
		       win_priv_t *wincon   /* Window constants      */
		       )
{

  /* Get the tangential velocity to the edge                         */
  vel_tan_2d(window, windat, wincon);

  /* Get the cell centred east and north velocities                  */
  vel_cen(window, windat, wincon, windat->u1av, windat->u2av, 
	  windat->uav, windat->vav, NULL, NULL, 1);
}

/* END vel_components_2d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute tangential velocity at faces                   */
/*-------------------------------------------------------------------*/
void vel_tan_2d(geometry_t *window,    /* Window geometry            */
		window_t *windat,      /* Window data                */
		win_priv_t *wincon     /* Window constants           */
  )
{
  int e, ee, es, eoe, n;
  double fs = 1.0;

  /*-----------------------------------------------------------------*/
  /* Get the u2 velocity at the e1 face. Note : yp1 map must be able */
  /* to map to the u2 boundary using set_map_e2. Maps must be re-set */
  /* so as to be self mapping  over this boundary after this routine */
  /* is called.                                                      */
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    windat->u2av[e] = 0.0;
    for (n = 1; n <= window->nee[e]; n++) {
      eoe = window->eSe[n][e];
      if (!eoe) continue;
      fs = window->h1au1[eoe] / window->h2au1[e];
      windat->u2av[e] += fs * window->wAe[n][e] * windat->u1av[eoe];
    }
  }
}

/* END vel_tan_2d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an edge                   */
/*-------------------------------------------------------------------*/
void mdxs(geometry_t *window,   /* Window geometry                   */
          win_priv_t *wincon    /* Window constants                  */
	  )
{
  int e, ee, em;                    /* Sparse coordinate/counter     */
  int c1, c2;
  double h1av;

  for (ee = 1; ee <= window->n2_e1; ee++) {
    e = window->w2_e1[ee];
    em = window->em[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    h1av =
      window->h1acell[em] / (window->h1acell[em] + window->h1acell[e]);
    wincon->mdx[e] =
      h1av * (wincon->Ds[c1] - wincon->Ds[c2]) + wincon->Ds[c2];
  }
}

/* END mdxs()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an edge at the forward    */
/* time step.                                                        */
/*-------------------------------------------------------------------*/
double mdxns(geometry_t *window,  /* Window geometry                 */
             window_t *windat,    /* Window data                     */
             win_priv_t *wincon,  /* Window constants                */
             int e                /* Sparse coordiante               */
  )
{
  int e2 = window->m2de[e];      /* 2D coordinate                    */
  int em = window->em[e2];
  int c1 = window->e2c[e2][0];
  int c2 = window->e2c[e2][1];
  double h1av;

  if (!wincon->sigma)
    return (1.0);
  h1av =
    window->h1acell[em] / (window->h1acell[em] + window->h1acell[e2]);
  return (h1av *
          (windat->eta[c1] - wincon->Hs[c1] -
           (windat->eta[c2] - wincon->Hs[c2])) + windat->eta[c2] -
          wincon->Hs[c2]);
}

/* END mdxns()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the water depth at an edge at the backward   */
/* time step.                                                        */
/*-------------------------------------------------------------------*/
double mdxbs(geometry_t *window,  /* Window geometry                 */
             window_t *windat,    /* Window data                     */
             win_priv_t *wincon,  /* Window constants                */
             int e                /* Sparse coordiante */
  )
{
  int e2 = window->m2de[e];       /* 2D coordinate                   */
  int em = window->em[e2];
  int c1 = window->e2c[e2][0];
  int c2 = window->e2c[e2][1];
  double h1av;

  if (!wincon->sigma)
    return (1.0);
  h1av =
    window->h1acell[em] / (window->h1acell[em] + window->h1acell[e2]);
  return (h1av *
          (windat->etab[c1] - wincon->Hs[c1] -
           (windat->etab[c2] - wincon->Hs[c2])) + windat->etab[c2] -
          wincon->Hs[c2]);
}

/* END mdxbs()                                                       */
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
  int n, ee, e;                 /* Edge coordinates, counters        */
  double colflux;               /* Velocity transport divergence     */
  double *u1flux = wincon->d2;  /* Flux in e1 direction              */
  double d1;
  /*
  FILE *fp = fopen("aa.ts","a");
  FILE *op = fopen("bb.ts","w");
  */
  /*-----------------------------------------------------------------*/
  /* Get the fluxes at time t over the whole window                  */
  memset(u1flux, 0, window->szeS * sizeof(double));
  /* Calculate the flux at e1 wet and boundary cells                 */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    u1flux[e] = windat->u1av[e] * windat->depth_e1[e] * window->h1au1[e] *
      wincon->mdx[e] * windat->dt2d;
  }

  /* Calculate the flux at auxiliary cells.                          */
  for (ee = window->b2_e1 + 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    if (!u1flux[e])
      u1flux[e] = windat->u1av[e] * windat->depth_e1[e] * window->h1au1[e] *
	wincon->mdx[e] * windat->dt2d;
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
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      windat->u1flux[e] += u1flux[e];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Debugging info                                                  */
  if (wincon->dbc && window->wn == wincon->dbw) {
    c = wincon->dbc;
    for (n = 1; n <= window->npe[c]; n++) {
      wincon->ba[n] = u1flux[window->c2e[n][c]];
    }
    debug_c(window, D_ETA, D_POST);
  }

  /*-----------------------------------------------------------------*/
  /* Update the elevation.                                           */
  /* Note : elevations are calculated on the boundaries here and     */
  /* overwritten in the boundary routine. A boundary condition of    */
  /* NOTHIN will use the elevations calculated here.                 */
  /*set_map_eta(window);*/
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = 0.0;
    for (n = 1; n <= window->npe[c]; n++) {
      e = window->c2e[n][c];
      colflux += window->eSc[n][c] * u1flux[e];
    }

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

  }
  /*
  colflux=0.0;
  c=999;
  while(windat->eta[c]>0.05) {
    c=window->c2c[2][c];
    colflux+=200.0;
    fprintf(op,"%d %f\n",window->s2j[c]*200,windat->eta[c]);
  }
  d1 = windat->u1[window->c2e[4][999]];
  fprintf(fp,"%f %f %f\n",windat->days,colflux,wincon->g*windat->eta[999]*windat->eta[999]/(2.0*wincon->Cd[999]*d1*d1));
  fclose(fp);
  fclose(op);
  */

  /*-----------------------------------------------------------------*/
  /* Adjust the updated elevation due to eta relaxation.             */
  if (wincon->etarlx & (RELAX|ALERT|ETA_TPXO)) {
    double rr = 0.0;
    /*if (wincon->etarlx & RELAX) rr = windat->etarlxtc;*/
    if (wincon->etarlx & (RELAX|ETA_TPXO)) {
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
void bdry_eta(geometry_t *window,   /* Window geometry               */
              window_t *windat,     /* Window data                   */
              win_priv_t *wincon    /* Window constants              */
	      )
{
  int n, m;                         /* Counters                      */
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

      /* Correct for atmospheric pressure                            */
      if (open[n]->inverse_barometer && !open[n]->adjust_flux) {
	double ramp = (wincon->rampf & INV_BARO) ? windat->rampval : 1.0;
        for (cc = 1; cc <= open[n]->no2_t; cc++) {
          c = open[n]->obc_t[cc];
          wincon->neweta[c] += ramp * (wincon->ambpress - windat->patm[c]) /
            (windat->dens[c] * wincon->g);
        }
      }

      /* Set eta to zero at boundary cells one cell deep             */
      if (open[n]->stagger & INFACE) {
	int c1, c2;
        for (cc = 1; cc <= open[n]->no2_t; cc++) {
          c = open[n]->obc_t[cc];
          c1 = open[n]->oi1_t[cc];
          c2 = open[n]->oi2_t[cc];
	  if (c1 == c2) wincon->neweta[c] = 0.0;
	}
      }

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

      /* Get the rate of change of elevation on the open boundaries  */
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
        c = open[n]->obc_t[cc];
        windat->detadt[c] =
          (wincon->neweta[c] - windat->etab[c]) / windat->dt2d;
      }
    }


    /* Average the sea level over intersecting OBCs if required      */
    if (open[n]->options & OP_OBCCM)
      average_OBC_corner(window, open[n], wincon->neweta, 0);

    /* Set the elevation at the OBC ghost cells equal to that at the */
    /* elevation open boundary location. This is required to set     */
    /* dzu1[] and dzu2[] on the velocity open boundaries, which is   */
    /* in turn used in velocity_adjust() to set vertical integrals   */
    /* of 3D velocity equal to the 2D velocity.                      */
    reset_bdry_eta(window, windat, wincon, open[n], wincon->neweta);

    /* Set relaxation on the boundary                                */
    if (open[n]->relax_ele) {
      double rb = open[n]->rele_b;
      double ri = open[n]->rele_i;
      int nr = open[n]->relax_ele;
      int *imap;
      double rc;
      int e, ee, bn;
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	c = open[n]->obc_e2[ee];
	imap = open[n]->nmape[e];
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
  int n, m, j, c, c1, cc, e, ee;

  if (open->stagger & OUTFACE) {
    for (ee = 1; ee <= open->no2_e1; ee++) {
      c = open->ogc_t[ee];
      c1 = open->obc_e2[ee];
      eta[c] = eta[c1];
      for (m = 0; m < open->bgz; m++) {
	c = open->omape[ee][c];
	eta[c] = eta[c1];
      }
    }
    /* Set the elevation at auxiliary OBC cells                      */
    for (cc = 1; cc <= open->no2_a; cc++) {
      c = open->obc_a[cc];
      for (j = 1; j <= window->npe[c]; j++) {
	if ((e = open->bec[j][cc])) {
	  c1 = open->omape[e][c];
	  for (m = 0; m < open->bgz; m++) {
	    eta[c1] = eta[c];
	    c1 = open->omape[e][c];
	  }
	}
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
void set_map_eta(geometry_t *window /* Window geometry               */
  )
{
  int n, j;
  int ee, e, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries   */
  /* to be self-mapping on the open boundary.                        */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & OUTFACE) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	j = (open[n]->ceni[ee]) ? 0 : 1;
	c = window->m2d[window->e2c[e][j]];
	window->c2c[open[n]->outi[ee]][c] = c;
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
void set_map_inside(geometry_t *window /* Window geometry            */
  )
{
  int n, j;
  int ee, e, c;
  open_bdrys_t **open = window->open;

  /* Set any self-mapping maps interior to eastern open boundaries   */
  /* to be self-mapping on the open boundary.                        */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	j = (open[n]->ceni[ee]) ? 0 : 1;
	c = window->m2d[window->e2c[e][j]];
	window->c2c[open[n]->outi[ee]][c] = c;
      }
    }
  }
}

/* END set_map_inside()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to apply the Asselin filter to remove the computational   */
/* mode. This is done for velocity only before the elevation         */
/* calculation so that elevation is computed with the filtered       */
/* fluxes and thus the divergence of 3D adjusted velocities is       */
/* consistent with elevation change in the surface layer             */
/* (conservation is maintained for tracers).                         */
/*-------------------------------------------------------------------*/
void asselin(geometry_t *window,  /* Window geometry                 */
             window_t *windat,    /* Window data                     */
             win_priv_t *wincon   /* Window constants                */
  )
{
  int e, ee, e1;                /* Edge coordinate                   */
  int c, cc, ci;                /* Cell coordinate                   */
  int n, j;                     /* Counters                          */
  double aconst = 0.1;          /* Time filter constant              */
  double mdx, mdy;              /* Total depth at time t             */
  double mdxb, mdyb;            /* Total depth at time t-1           */
  double *u1av = wincon->d1;
  int *mask = wincon->i7;
  double yr = 86400.0 * 365.0;
  int checkf = 0;
  /*int checkf = 95784;*/
  /*int checkf = 49149;*/
  double mf=0.0, mfe[window->npem+1], af[window->npem+1];

  for (j = 1; j <= window->npem; j++) mfe[j] = 0.0;

  memcpy(u1av, windat->u1av, window->szeS * sizeof(double));
  memset(mask, 0, window->szcS * sizeof(int));
  
  /* For sigma the total depth at time t+1 is not available until    */
  /* after eta_step(), hence the time filtering is done in terms of  */
  /* velocity*depth instead of velocity. Note that nu1av=vel*depth   */
  /* at this stage.                                                  */
  if (wincon->sigma) {
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      mdx = wincon->mdx[e];
      mdxb = mdxbs(window, windat, wincon, e);
      windat->u1av[e] = (windat->u1av[e] * mdx + 0.5 * aconst *
                         (windat->u1avb[e] * mdxb -
                          2.0 * windat->u1av[e] * mdx +
                          windat->nu1av[e])) / mdx;
    }
  } else {
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      windat->u1av[e] = windat->u1av[e] + 0.5 * aconst *
        (windat->u1avb[e] - 2.0 * windat->u1av[e] + windat->nu1av[e]);
    }
  }

  /* If the open boundary condition for 2D velocity is VERTIN, then  */
  /* do not perform time filtering (restore boundary velocities to   */
  /* those before filtering is performed).                           */
  reset_map_t_all(window);
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    scale_details_t *scale = open->sdata_e;

    if (open->bcond_tan2d & (VERTIN|FILEIN|CUSTOM|TIDALC)) {
      for (ee = open->no3_e1+1; ee <= open->to2_e1; ee++) {
	e = open->obc_e1[ee];
	ci = open->obc_e2[ee];
	if (!mask[ci])
	  windat->u1av[e] = u1av[e];
      }
    }
    if (open->bcond_nor2d & (VERTIN|FILEIN|CUSTOM|TIDALC)) {
      double f1, f2, df; 
      double nvel[window->npem+1], v1[window->npem+1], v2[window->npem+1];
      double rts, eta, tide, sgn[window->npem+1];
      double rtsh = 1.2;

      /*-------------------------------------------------------------*/
      /* Edge boundaries                                             */
      if (!(open->bcond_nor2d & FLATHR)) {
	for (ee = 1; ee <= open->no2_e1; ee++) {
	  e = open->obc_e1[ee];
	  ci = open->obc_e2[ee];
	  if (!mask[ci])
	    windat->u1av[e] = u1av[e];
	}
      }
      /*-------------------------------------------------------------*/
      /* Adjust the boundary velocity if required                    */
      if (open->adjust_flux) {
	double adjust_flux;

	if (wincon->rampf & FLUX_ADJUST && windat->rampval < 1.0) {
	  double dtr = wincon->rampend - wincon->rampstart;
	  adjust_flux = (windat->t - wincon->rampstart) * 
	    (open->adjust_flux - yr) / dtr + yr;
	} else
	  adjust_flux = open->adjust_flux;
	rts = windat->dtb2 / adjust_flux;
	
	/* Loop through the U1 boundary cells                        */
	for (cc = 1; cc <= open->no2_t; cc++) {
	  c = open->obc_t[cc];

	  if (mask[c]) continue;
	  
	  /* Get the default flux adjustment                         */
	  if (open->adjust_flux <= 0) {
	    adjust_flux = 0.0;
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		/*
		adjust_flux += sqrt(window->cellarea[c]) /
		  sqrt(wincon->g * windat->depth_e1[e] * wincon->Ds[c]);
		*/
		adjust_flux += (window->h2au1[e] /
				sqrt(wincon->g * windat->depth_e1[e] * wincon->Ds[c]));
	      }
	    }
	    adjust_flux *= open->nepc[cc];
	    rts = windat->dtb2 / adjust_flux;
	    /*
	    if (wincon->rampf & FLUX_ADJUST && windat->rampval < 1.0) {
	      double dtr = wincon->rampend - wincon->rampstart;
	      adjust_flux = (windat->t - wincon->rampstart) * 
		(adjust_flux - 1.0 / rtsh) / dtr + 1.0 / rtsh;
	    }
	    */
	  }
	  /* Hard relaxation on corners with overlap                 */
	  if ((open->options & OP_OBCCM) && !open->olap[cc]) {
	    adjust_flux = rtsh * windat->dtb2;
	    rts = windat->dtb2 / adjust_flux;
	  }
	  /* Save the flux adjustment diagnostic if required         */
	  if (windat->obc_phase) {
	    if (open->adjust_flux_s)
	      windat->obc_phase[c] = windat->u1av[e];
	    else
	      windat->obc_phase[c] = adjust_flux;
	  }

	  /* Get the fluxes                                          */
	  f2 = 0.0;
	  for (j = 1; j <= window->npe[c]; j++) {
	    e = open->bec[j][cc];
	    e1 = window->c2e[j][c];
	    if (e == e1) {
	      sgn[j] = -1.0 * (double)window->eSc[j][c];
	      mfe[j] += window->eSc[j][c] * windat->u1av[e1] * windat->depth_e1[e1] * 
		window->h1au1[e1] * wincon->mdx[e1] * windat->dt2d;
	      mf += mfe[j];
	    } else {
	      f2 += window->eSc[j][c] * windat->u1av[e1] * windat->depth_e1[e1] * 
		window->h1au1[e1] * wincon->mdx[e1] * windat->dt2d;
	      sgn[j] = -1.0 * (double)window->eSc[j][c];
	    }
	  }

	  /* Get the boundary flux that would result in a flux       */
	  /* divergence sufficient to result in sea level equal to   */
	  /* eta_rlx[].                                              */
	  eta = 0.0;
	  if (open->bcond_ele & (CUSTOM|FILEIN))
	    eta += ((wincon->compatible & V1670) ? windat->eta_rlx->val1[c] :
		    open->transfer_eta[cc]);

	  /* Scaling                                                 */
	  if (scale->type == (TRSC_SUM | TRSC_NUM)) {
	    eta += scale->fact;
	  } else if (scale->type == (TRSC_SUM | TRSC_TRA)) {
	    int trn = scale->ntr;
	    eta += windat->tr_wcS[trn][c];
	  } else if (scale->type == (TRSC_PCT | TRSC_NUM)) {
	    eta *= scale->fact;
	  } else if (scale->type == (TRSC_PCT | TRSC_TRA)) {
	    int trn = scale->ntr;
	    eta *= windat->tr_wcS[trn][c];
	  }

	  /* Tidal elevations                                        */
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
	    tide = ramp * bc_tidal(window, windat, wincon, open->tideforce, c);
	  }
	  if (open->inverse_barometer) {
	    double ramp = (wincon->rampf & INV_BARO) ? windat->rampval : 1.0;
	    eta += ramp * (wincon->ambpress - windat->patm[c]) /
	      (windat->dens[c] * wincon->g);
	  }
	  /* Adaptive flux adjustment
	     adjust_flux = get_flux_adjust(window, windat, wincon, open, 
	     windat->eta[c] - tide - eta, e, c);
	     rts = windat->dtb2 / adjust_flux;
	  */

	  /* Dual relaxation: Relax hard to the tidal signal         */
	  if (open->adjust_flux_s && open->bcond_ele & (TIDALH|TIDALC|TIDEBC)) {
	    double depth, etat, df, etadiff;
	    double rtst = (open->adjust_flux_s < 0.0) ? rts : windat->dtb2 / open->adjust_flux_s;
	    f2 *= 0.5;
	    df = f2;
	    
	    /* Get the tidal flux                                    */
	    f1 = (tide - windat->etab[c]) * window->cellarea[c] + f2;
	    /* Convert this flux to a 2D velocity                    */
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		nvel[j] = open->nepc[cc] * sgn[j] * f1 / 
		  (windat->depth_e1[e] * window->h1au1[e] * 
		   wincon->mdx[e] * windat->dt2d);
		v1[j] = nvel[j];
	      }
	    }
	    /* Relax the boundary velocity to this value             */
	    f1 = 0.0;
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		windat->u1av[e] -= rtst * (windat->u1av[e] - nvel[j]);
		/* Get the flux due to this velocity                 */
		f1 = windat->u1av[e] * windat->depth_e1[e] * 
		  window->h1au1[e] * wincon->mdx[e] * windat->dt2d;
		/* Get the elevation due to this flux                */
		f2 -= sgn[j] * f1;
	      }
	    }
	    etat = windat->etab[c]  - f2 / window->cellarea[c];
	    
	    /* Get the low frequency flux                            */
	    f1 = (eta + tide - etat) * window->cellarea[c] + df;
	    /* Convert this flux to a 2D velocity                    */
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		depth = etat - window->botzu1[e];
		nvel[j] = open->nepc[cc] * sgn[j] * f1 / 
		  (depth * window->h1au1[e] * 
		   wincon->mdx[e] * windat->dt2d);
		v2[j] = nvel[j];
	      }
	    }
	    
	    /* Get the velocity required for relaxation              */
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		nvel[j] = (windat->u1av[e] + rts * v2[j]) / rts;
	      }
	    }
	  } else {  /* Single relaxation                             */
	    f1 = (eta + tide - windat->etab[c]) * window->cellarea[c];
	    /* Convert this flux to a 2D velocity                    */
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		nvel[j] = open->nepc[cc] * sgn[j] * (f1 + f2) / 
		  (windat->depth_e1[e] * window->h1au1[e] * 
		   wincon->mdx[e] * windat->dt2d);
	      }
	    }
	  }

	  /* Relax the boundary velocity to this value               */
	  for (j = 1; j <= window->npe[c]; j++) {
	    if ((e = open->bec[j][cc])) {
	      windat->u1av[e] -= rts * (windat->u1av[e] - nvel[j]);
	    }
	  }
	  mask[c] = 1;

	  /* Save the flux adjustment diagnostic if required         */
	  if (windat->obc_phase && open->adjust_flux_s) {
	    for (j = 1; j <= window->npe[c]; j++) {
	      if ((e = open->bec[j][cc])) {
		windat->obc_phase[c] = -windat->dtb2 * (windat->obc_phase[c] - 
				      open->nepc[cc] * (v1[j] + v2[j])) /
		  (windat->u1av[e] - windat->obc_phase[c]);
	      }
	    }
	  }
	  if (wincon->ic == 0 && window->wsa[c] == checkf) {
	    double vel;
	    f2 = 0.0;
	    for (j = 1; j <= window->npe[c]; j++) {
	      e = window->c2e[j][c];
	      vel = windat->u1av[e];
	      /*if (open->bec[j][cc]) vel = nvel[j];*/
	      f2 += window->eSc[j][c] * vel * windat->depth_e1[e] * 
		window->h1au1[e] * wincon->mdx[e] * windat->dt2d;
	    }
	    f1 = windat->etab[c] - f2 / window->cellarea[c];
	    printf("check: %f %d[%f %f] actual=%f eta=%f tide=%f new=%f\n",windat->days, c, 
		   window->cellx[c], window->celly[c], windat->eta[c], eta, tide, f1);
	  }
	}
      }
    }
  }
}

/* END asselin()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the u1 velocity from the u1 velocity transport */
/*-------------------------------------------------------------------*/
void extract_velocity_2d(geometry_t *window,  /* Window geometry     */
			 window_t *windat,    /* Window data         */
			 win_priv_t *wincon   /* Window constants    */
  )
{
  int e, ee, n;                 /* Local sparse coordinate / counter */

  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    wincon->Hn1[e] = mdxns(window, windat, wincon, e);
    windat->u1av[e] /= wincon->Hn1[e];
  }
}

/* END extract_velocity_2d()                                         */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* LEAPFROG : Reset the velocities for the new time level.           */
/*-------------------------------------------------------------------*/
void leapfrog_update_2d(geometry_t *window,   /* Window geometry     */
			window_t *windat,     /* Window data         */
			win_priv_t *wincon    /* Window constants    */
			)
{

  /*-----------------------------------------------------------------*/
  /* Set the velocities for the new time level                       */
  memcpy(windat->u1avb, windat->u1av, window->szeS * sizeof(double));
  memcpy(windat->u1av, windat->nu1av, window->szeS * sizeof(double));
  memcpy(windat->etab, windat->eta, window->szcS * sizeof(double));
  memcpy(windat->eta, wincon->neweta, window->szcS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Extract the velocities from the updated solution for sigma      */
  if (wincon->sigma)
    extract_velocity_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the mean u1av velocity if required                          */
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
    int ee, e, c, cc;
    double t = windat->dtf2;
    if (windat->u1am && windat->u2am) {
      for (cc = 1; cc <= window->b2_t; cc++) {
        c = window->w2_t[cc];
	windat->u1am[c] = (windat->u1am[c] * windat->meanc[c] + 
			   windat->uav[c] * t) / (windat->meanc[c] + t);
	windat->u2am[c] = (windat->u2am[c] * windat->meanc[c] + 
			   windat->vav[c] * t) / (windat->meanc[c] + t);
      }
    }
    if (windat->uame) {
      for (ee = 1; ee <= window->b2_e1; ee++) {
        e = window->w2_e1[ee];
	c = window->e2c[e][0];
	if (window->wgst[c]) c = window->e2c[e][1];
	windat->uame[e] = (windat->uame[e] * windat->meanc[c] + 
			   windat->u1av[e] * t) / (windat->meanc[c] + t);
      }
    }
  }
}

/* END leapfrog_update_2d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to calculate the total water depth at edges               */
/*-------------------------------------------------------------------*/
void get_depths(geometry_t *window,   /* Window geometry             */
		window_t *windat,     /* Window data                 */
		win_priv_t *wincon    /* Window constants            */
  )
{
  int c, cc, c1, c2;           /* Local sparse coordinate / counter  */
  int e, ee, ei;               /* Edge coordinate                    */
  double top;                  /* Surface elevation at the edge      */
  double *tzp;                 /* Surface height array               */
  double *depth;               /* Pointer to depth array             */
  int n;
  open_bdrys_t **open = window->open;

  /*-----------------------------------------------------------------*/
  /* tzp is a pointer which points to either topz[][] or eta[][].    */
  /* topz[][] only gets updated once every 3d step. For the non-     */
  /* linear case, we need to use eta instead, which gets updated     */
  /* every 2d step. This value is only used to calculate transport   */
  /* (nothing to do with the surface slope term, which always uses   */
  /* eta).                                                           */
  tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : windat->topz;

  /*-----------------------------------------------------------------*/
  /* Total depth at edges                                            */
  /* Set the pointers.                                               */
  depth = windat->depth_e1;

  /*-----------------------------------------------------------------*/
  /* Loop over the u1 edges. Note the depth on the boundary is set   */
  /* in bdry_u1_2d().                                                */
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    top = max(tzp[c1], tzp[c2]);
    depth[e] = top - window->botzu1[e];
  }

  /*-----------------------------------------------------------------*/
  /* Set the total depth on the boundary                             */
  for (n = 0; n < window->nobc; n++) {
    /* For TIDALC velocity OBCs with TPXO transport input, revert to */
    /* the transport provided by the custom tide file (i.e. do not   */
    /* account for elevation in the total depth).                    */
    if (open[n]->bcond_nor2d & TIDALC && wincon->tidef & TD_TRAN) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	depth[e] = -window->botzu1[e];
      }
    } else {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	c = open[n]->obc_e2[ee];
	depth[e] = tzp[c] - window->botzu1[e];
      }
    }
  }

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->bcond_tan2d & TIDALC && wincon->tidef & TD_TRAN) {
      for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	depth[e] = -window->botzu1[e];
      }
    } else {
      for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to2_e1; ee++) {
	e = open[n]->obc_e1[ee];
	c1 = window->e2c[e][0];
	c2 = window->e2c[e][1];
	top = max(tzp[c1], tzp[c2]);
	depth[e] = top - window->botzu1[e];
      }
    }
  }

  /* Set the depth at ghost cells                                    */
  for (ee = 1; ee <= window->nbpte1S; ee++) {
    e = window->bpte1S[ee];
    ei = window->bine1S[ee];
    depth[e] = depth[ei];
  }

  /* Decrease the depth for sub-surface explicit maps if required    */
  if (window->sm_e1) {
    int cs;
    for (cc = 1; cc <= window->b2_e1; cc++) {
      cs = window->w2_e1[cc];
      c = window->sm_e1[cs];
      if (c)
      depth[cs] = window->gridz[window->zp1[c]] - window->botzu1[cs];
    }
  }
}

/* END get_depths()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_tr_2d(geometry_t *window,  /* Window geometry       */
                       window_t *windat,    /* Window data           */
                       double trem, /* Time remain in sub-step       */
                       double *vel, /* Velocity array                */
                       double *as,  /* Multi-dt velocity at time t   */
                       double *ae   /* Multi-dt velocity at time t+1 */
		       )
{
  int cc, c;                        /* Sparse counter / coordinate   */

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
void get_mode3d_eta(geometry_t *window,       /* Window geometry     */
		    window_t *windat,         /* Window data         */
		    win_priv_t *wincon        /* Window constants    */
		    )
{
  int ic = wincon->ic;          /* 2D step counter                   */
  int iratio = windat->iratio;  /* 3D/2D ratio                       */
  int c, cc;                    /* Sparse coordinate / counter       */
  double aconst = 0.1;          /* Time filter constant              */
  double eta3;                  /* 3D mode elevation                 */

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
void vint_3d(geometry_t *window,       /* Window geometry            */
	     window_t *windat,         /* Window data                */
	     win_priv_t *wincon        /* Window constants           */
	     )
{
  int ee, e, es, zm1;
  double depth;

  memset(windat->u1av, 0, window->szeS * sizeof(double));
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    es = window->m2de[e];
    zm1 = window->zm1e[e];
    depth = 0.0;
    while(e != zm1) {
      if (!isnan(windat->u1[e])) {
	windat->u1av[es] += windat->u1[e] * windat->dzu1[e];
	depth += windat->dzu1[e];
      }
      e = zm1;
      zm1 = window->zm1e[zm1];
    }
    windat->u1av[es] /= depth;
  }
}

/* END vint_3d()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check if the model is going unstable and exit if so    */
/*-------------------------------------------------------------------*/
int check_unstable(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon, /* Window constants           */
		   int mode)           /* Instability mode           */
{
  double maxwind = 100.0;   /* Maximum wind allowed (ms-1)           */
  int c, cc, e, ee;

  /* Check for sea level fatalities                                  */
  if (wincon->fatal & ETA_A && mode & ETA_A) {

    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];

      if (fabs(wincon->neweta[c]) > wincon->etamax) {
	if (window->nwindows > 1) {
	  write_site(window, window->cellx[c], window->celly[c], "eta");
	  if (window->us_type & US_IJ)
	    hd_quit_and_dump
	      ("etastep: Surface exceeded ETAMAX (%5.2f) at c=%d(i=%d j=%d [%f %f]) cg=%d wn=%d t=%8.3f days\n",
	       wincon->neweta[c], c, window->s2i[c], window->s2j[c],
	       window->cellx[c], window->celly[c], window->wsa[c], window->wn, windat->t / 86400);
	  else
	    hd_quit_and_dump
	      ("etastep: Surface exceeded ETAMAX (%5.2f) at c=%d(c2cc=%d cg=%d [%f %f]) wn=%d t=%8.3f days\n",
	       wincon->neweta[c], c, window->s2i[c], window->wsa[c], 
	       window->cellx[c], window->celly[c], window->wn, windat->t / 86400);
	} else {
	  write_site(window, window->cellx[c], window->celly[c], "eta");
	  if (window->us_type & US_IJ)
	    hd_quit_and_dump
	      ("etastep: Surface exceeded ETAMAX (%5.2f) at c=%d(i=%d j=%d [%f %f]) t=%8.3f days\n",
	       wincon->neweta[c], c, window->s2i[c], window->s2j[c],
	       window->cellx[c], window->celly[c], windat->t / 86400);
	  else
	    hd_quit_and_dump
	      ("etastep: Surface exceeded ETAMAX (%5.2f) at c=%d(c2cc=%d cg=%d [%f %f]) t=%8.3f days\n",
	       wincon->neweta[c], c, window->s2i[c], window->wsa[c],
	       window->cellx[c], window->celly[c], windat->t / 86400);
	  return(1);
	}
      }
      if (wincon->fatal & NANF) {
	if(isnan(wincon->neweta[c])) {
	  if (window->nwindows == 1) {
	    write_site(window, window->cellx[c], window->celly[c], "eta NaN");
	    hd_quit_and_dump("etastep: NaN at %d (%d %d) t=%8.3f days\n", c,
			     window->s2i[c], window->s2j[c], windat->t / 86400);
	    return(1);
	  } else {
	    write_site(window, window->cellx[c], window->celly[c], "eta NaN");
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
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      c = window->e2c[e][0];
      if (fabs(windat->u1av[e]) > wincon->velmax2d) {
	write_site(window, window->u1x[e], window->u1y[e], "u1av");
	hd_quit_and_dump
	  ("vel2d: u1av velocity exceeded velmax (%f) at %d(%d %d) t=%8.3f days\n",
	   windat->u1av[e], e, window->s2i[c], window->s2j[c], windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u1av[e])) {
	  if (window->nwindows == 1) {
	    write_site(window, window->u1x[e], window->u1y[e], "u1av NaN");
	    hd_quit_and_dump("vel2d: u1av NaN at %d(%d %d) t=%8.3f days\n",
			     e, window->s2i[c], window->s2j[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    write_site(window, window->u1x[e], window->u1y[e], "u1av NaN");
	    hd_quit_and_dump("vel2d: u1av NaN at %d (%d %d) t=%8.3f days : window %d\n",
			     e, window->s2i[c], window->s2j[c], 
			     windat->t / 86400, window->wn);
	    return(1);
	  }
	}
      }
    }
  }

  /* Check for 3D velocity fatalities                                */
  if (wincon->fatal & VEL3D && mode & VEL3D) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      c = window->e2c[e][0];
      if (fabs(windat->u1[e]) > wincon->velmax) {
	write_site(window, window->u1x[window->m2de[e]], window->u1y[window->m2de[e]], "u1");
	hd_quit_and_dump
	  ("vel3d: u1 velocity exceeded velmax (%f) at %d(%d %d %d) t=%8.3f days\n",
	   windat->u1[e], e, window->s2i[c], window->s2j[c], window->s2k[c], 
	   windat->t / 86400);
	return(1);
      }
      if (wincon->fatal & NANF) {
	if(isnan(windat->u1[e])) {
	  if (window->nwindows == 1) {
	    hd_quit_and_dump("vel3d: u1 NaN at %d(%d %d %d) t=%8.3f days\n",
			     e, window->s2i[c], window->s2j[c], window->s2k[c], 
			     windat->t / 86400);
	    return(1);
	  } else {
	    hd_quit_and_dump("vel3d: u1 NaN at %d(%d %d %d) t=%8.3f days : window %d\n",
			     e, window->s2i[c], window->s2j[c], window->s2k[c], 
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
    for (ee = 1; ee <= window->b2_t; ee++) {
      e = window->w2_t[ee];
      c = window->e2c[e][0];
      if (fabs(windat->wind1[e]) > maxstr) {
	hd_quit_and_dump
	  ("wind1: Wind exceeded 100ms-1 (%f) at %d(%d %d) t=%8.3f days\n",
	   windat->wind1[e], e, window->s2i[c], window->s2j[c],
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
	write_site(window, window->cellx[window->m2d[c]], window->celly[window->m2d[c]], "temp");
	hd_quit_and_dump
	  ("temp: NaN found  at (%d %d %d) t=%8.3f days : window %d\n",
	   window->s2i[c], window->s2j[c], window->s2k[c],
	   windat->t / 86400, window->wn);
	return(1);
      }
      if (isnan(windat->sal[c])) {
	write_site(window, window->cellx[window->m2d[c]], window->celly[window->m2d[c]], "salt");
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


/*-------------------------------------------------------------------*/
/* Estimates the bottom stress for tidal flows, where an ensemble of */
/* bottom drag is created and velocities at the next time-step are   */
/* optimized against TPXO tidal velocity.                            */
/*-------------------------------------------------------------------*/
double est_bot_stress(geometry_t *window, /* Window geometry         */
		      window_t *windat,   /* Window data             */
		      win_priv_t *wincon, /* Window constants        */
		      int e,
		      double botu1,
		      double val,
		      double pgt,
		      double cot,
		      double uin
		      )
{
  int e2 = window->m2de[e];
  int c1 = window->e2c[e2][0];
  int c2 = window->e2c[e2][1];
  double ux, uy, ut, nut;
  double *depth = windat->depth_e1;
  double midx = wincon->mdx[e2];
  double bd, bft, u, v, nu;
  double bds = 0.0;      /* Ensemble start bottom drag               */
  double bde = 0.0001;   /* Bottom drag increment                    */
  double Cd, rmse, rm;
  double gmt = windat->t + windat->dt2d - wincon->tz;
  int nb = 50;
  int i;
  /* Get the TPXO edge tidal velocity at the next time-step          */
  /*
  ux = 0.5 * (windat->tpxotideu[c1] + windat->tpxotideu[c2]);
  uy = 0.5 * (windat->tpxotidev[c1] + windat->tpxotidev[c2]);
  */
  ux = csr_tide_eval(&wincon->tcu, e2, gmt);
  uy = csr_tide_eval(&wincon->tcv, e2, gmt);
  ut = (ux * window->costhu1[e2] + uy * window->sinthu1[e2]);

  /* Get the velocity (minus bottom friction at the next time-step   */
  nu = uin + windat->dt2d * (midx * (pgt + cot) + wincon->u1inter[e]);

  /* Create the ensemble and optimize                                */
  rm = HUGE;
  for (i = 1; i <= 50; i++) {
    bd = bds + (double)i * bde;
    v = bd * max(wincon->uf, val);
    /* Truncate to ensure stability                                  */
    if (v > depth[e] * midx / windat->dt2d)
      v = depth[e] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case                               */
    bft = -v * botu1 / depth[e];
    u = nu + windat->dt2d * bft;
    rmse = sqrt((u-ut) * (u-ut));
    if (rmse < rm) {
      rm = rmse;
      nut = u;
      Cd = bd;
    }
  }
  wincon->d2[e] = rm;
  wincon->d3[e] = Cd;
  return(Cd);
}

/* END est_bot_stress()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the error from TPXO velocity divergence relative to the   */
/* TPXO surface height.                                              */
/*-------------------------------------------------------------------*/
double tpxo_error(geometry_t *window, /* Window geometry             */
		  window_t *windat,   /* Window data                 */
		  win_priv_t *wincon  /* Window constants            */
		  )
{
  int cc, c, ee, e, n;
  double ux, uy, ut;
  double eta, oeta, neta, colflux, depth;
  double *u1flux = wincon->d2;
  double gmt0 = windat->t - wincon->tz;
  double gmt1 = windat->t + windat->dt2d - wincon->tz;

  /* Get the TPXO velocities                                         */
  memset(u1flux, 0, window->szeS * sizeof(double));
  for (ee = 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    depth = windat->depth_e1[e];
    ux = csr_tide_eval(&wincon->tcu, e, gmt0);
    uy = csr_tide_eval(&wincon->tcv, e, gmt0);
    ut = (ux * window->costhu1[e] + uy * window->sinthu1[e]);
    u1flux[e] = ut * depth * window->h1au1[e] * wincon->mdx[e] * windat->dt2d;
  }

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];

    /*---------------------------------------------------------------*/
    /* Get the divergence of velocity transport                      */
    colflux = 0.0;
    for (n = 1; n <= window->npe[c]; n++) {
      e = window->c2e[n][c];
      colflux += window->eSc[n][c] * u1flux[e];
    }

    /* Calculate new etat value                            (etamark) */
    oeta = csr_tide_eval(&wincon->tc, c, gmt0);
    neta = csr_tide_eval(&wincon->tc, c, gmt1);
    eta = max(oeta - colflux / window->cellarea[c], window->botz[c]);
    windat->dum2[c] = neta - eta;
  }
}

/* END tpxo_error()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a site file for the location of instabilities.             */
/*-------------------------------------------------------------------*/
void write_site(geometry_t *window, double x, double y, char *tag)
{
  FILE *fp;
  char buf[MAXSTRLEN];

  sprintf(buf, "%scrash.site", master->opath); 
  if ((fp = fopen(buf, "w")) == NULL) return;
  fprintf(fp, "%f %f %s\n", x, y, tag);
  fclose(fp);
}

/* END write_site()                                                  */
/*-------------------------------------------------------------------*/
