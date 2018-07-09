/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/vel3d.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vel3d.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* mode3d_step()      : Updates the 2D mode variables                */
/* vel_u1_update()    : Solves the 2D u1 velocity equation           */
/* vel_u2_update()    : Solves the 2D u2 velocity equation           */
/* advect_u1_3d()      : Solves the u1 advection equation            */
/* advect_u2_3d()      : Solves the u2 advection equation            */
/* bdry_u1_3d()        : Sets normal and tangential u1 OBC's         */
/* bdry_u2_3d()        : Sets normal and tangential u2 OBC's         */
/* extract_u1_3d()     : Gets u1 velocity from u1 transport          */
/* extract_u2_3d()     : Gets u2 velocity from u2 transport          */
/*-------------------------------------------------------------------*/

void precalc_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void precalc_u2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void merge_thin_layers(int c, int zm1, double *vel, double *dz);
void reset_thin_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_thin_u2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_thin_wtop(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     int *ctp, int vcs, double *wtop);
void set_dry_bdry(geometry_t *window, int nb, int *obc, double *vel);
void set_w_dry(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_wvel(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	      double *wvel, int mode);
void mom_balance(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_sdc_e1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_sdc_e2(geometry_t *window, window_t *windat, win_priv_t *wincon);

/*-------------------------------------------------------------------*/
/* Window step part 1                                                */
/*-------------------------------------------------------------------*/
void mode3d_step_window_p1(master_t *master,
                           geometry_t *window,
                           window_t *windat, win_priv_t *wincon)
{

  double clock;

  /*---------------------------------------------------------------*/
  /* Fill 3D  velocities into the window data structures. This can */
  /* be done inside this  window loop  (contrary to tracers) since */
  /* velocity  is  only updated  after the 2D  mode and  hence all */
  /* windows are filled with velocities at the current time-step.  */
  TIMING_SET;
  master->win_data_fill_3d(master, window, windat, master->nwindows);
  TIMING_DUMP(2, "   fill_3d ");

  clock = dp_clock();
  
  /*---------------------------------------------------------------*/
  /* Calculate a heat and salt flux if required                    */
  TIMING_SET;
  calc_heatf(window, windat, wincon);
  calc_saltf(window, windat, wincon);
  TIMING_DUMP(2, "   h&s fluxes ");

  /* Set self-mapping maps over the blend zone if required         */
  blend_vel(window, windat, wincon, BLM1, NULL);

  init_sigma(window, windat, wincon);

#if !GLOB_BC
  /* Set the lateral boundary conditions for velocity at t for     */
  /* multiple windows. These velocities are required for the       */
  /* advective terms, and require transfers before ghost cells     */
  /* adjacent to auxiliary cells can be set.                       */
  if (master->nwindows > 1) {
    vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	      window->bpte1, window->bine1, wincon->slip);
    vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
	      window->bpte2, window->bine2, wincon->slip);
  }
  /* Set the lateral boundary conditions for velocity at t-1.      */
  /* These velocities are required for stress tensors. Ghost cells */
  /* are not set accurately since a lateral BC is not set on       */
  /* nu1/nu2 before the asselin() filtering; setting BCs here      */
  /* overcomes this.                                               */
  if (!(master->compatible & V1283)) {
    vel2D_lbc(windat->u1b, window->nbpte1, window->nbe1,
	      window->bpte1, window->bine1, wincon->slip);
    vel2D_lbc(windat->u2b, window->nbpte2, window->nbe2,
	      window->bpte2, window->bine2, wincon->slip);
  }
#endif

  /*---------------------------------------------------------------*/
  /* Set the vertical  mixing coefficients. This  is done  in this */
  /* window loop so that transferred velocities from other windows */
  /* can be used in the velocity shear term.                       */
  wincon->calc_closure(window, windat, wincon);
  bdry_closure(window, windat, wincon);

  /*---------------------------------------------------------------*/
  /* Evalulate the sources and sinks of water                      */
  ss_water(window, windat, wincon);

  /*---------------------------------------------------------------*/
  /* Get the cfl time-steps if required                            */
  if (!(wincon->cfl & NONE))
    calc_cfl(window, windat, wincon);

  /*---------------------------------------------------------------*/
  /* Set the vertical grid spacings at e1 and e2 faces.            */
  set_dz_at_u1(window, windat, wincon);
  set_dz_at_u2(window, windat, wincon);

  debug_c(window, D_INIT, D_POST);

  /*---------------------------------------------------------------*/
  /* Set the stress tensors and Smagorinsky diffusion.             */
  wincon->hor_mix->pre(window, windat, wincon);

  /*---------------------------------------------------------------*/
  /* Check for fatal wind instabilities                            */
  if (check_unstable(window, windat, wincon, WIND)) return;

  windat->wclk = (dp_clock() - clock);

  /*---------------------------------------------------------------*/
  /* Transfer grid spacings and mixing coefficients to the master. */
  //TIMING_SET;
  if (master->nwindows > 1)
    master->win_data_empty_3d(master, window, windat, MIXING);
  //  TIMING_DUMP(2, "   empty_3d ");

}

/* END mode3d_step_window_p1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window step part 2                                                */
/*-------------------------------------------------------------------*/
void mode3d_step_window_p2(master_t *master,
                           geometry_t *window, window_t *windat,
                           win_priv_t *wincon)
{

  double clock;
  windat->dt = windat->dtf + windat->dtb;

  /*-----------------------------------------------------------------*/
  /* Fill the window with u1 data (dzu1, dzu2, Vz, Kz) from the      */
  /* master.                                                         */
  master->win_data_refill_3d(master, window, windat, master->nwindows, MIXING);

  clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Update the 3D velocity                                          */
  vel_u1_update(window, windat, wincon);
  if (master->crf == RS_RESTART) return;
  vel_u2_update(window, windat, wincon);
  if (master->crf == RS_RESTART) return;

  /* Set sources of momentum                                         */
  ss_momentum(window, windat, wincon, VEL3D);
  if (wincon->tendf && wincon->tendency)
    mom_balance(window, windat, wincon);

  windat->dt = windat->dtf;
  windat->wclk += (dp_clock() - clock);
}

/* END mode3d_step_window_p2()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window post step part 1                                           */
/* Adjust velocities in each window so that the vertical integral    */
/* of the 3D velocity equalls the 2D velocity and transfer the       */
/* adjusted 3D velocities back to the master.                        */
/*-------------------------------------------------------------------*/
void mode3d_post_window_p1(master_t *master,
                           geometry_t *window, window_t *windat,
                           win_priv_t *wincon)
{

  double clock = dp_clock();

  /* Extract the velocities from the updated solution                */
  extract_u1_3d(window, windat, wincon);
  extract_u2_3d(window, windat, wincon);

  /* Invoke time filtering and step forward                          */
  leapfrog_update_3d(window, windat, wincon);

  /* Add the velocity increments to 3D velocity if required          */
  do_vel_relax(window, windat, wincon, VEL3D);
  /*do_vel_increment_3d(window);*/

  /* Calculate the 3D velocity alert diagnostics                     */
  alerts_w(window, VEL3D);

  /* Check for fatal instabilities                                   */
  if (check_unstable(window, windat, wincon, VEL3D)) return;

  /* Adjust  the  velocities. 3D  fluxes calculated  below must  use */
  /* the adjusted velocities.                                        */
  velocity_adjust(window, windat, wincon);

  /* Calculate the fluxes.                                           */
  set_flux_3d(window, windat, wincon, VEL3D);

  /* Set the velocities at newly wetted cells                        */
  set_new_cells_u1(window, windat, wincon);
  set_new_cells_u2(window, windat, wincon);

  /* Set the lateral boundary conditions for velocity.               */
#if !GLOB_BC
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
            window->bpte1, window->bine1, wincon->slip);
  vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
            window->bpte2, window->bine2, wincon->slip);
#endif

  windat->wclk += (dp_clock() - clock);

  /* Transfer adjusted 3D velocities and 3D fluxes to the master     */
  master->win_data_empty_3d(master, window, windat, VELOCITY);
  master->win_data_empty_3d(master, window, windat, CFL);
}

/* END mode3d_post_window_p1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window post step part 2                                           */
/* Set the 3D fluxes through the cell faces, calculate the           */
/* vertical velocity and transfer this to the master.                */
/*-------------------------------------------------------------------*/
void mode3d_post_window_p2(master_t *master,
                           geometry_t *window, window_t *windat,
                           win_priv_t *wincon)
{

  double clock;

  /* Fill  the window  with adjusted  velocities and 3D  fluxes from */
  /* the master.  These are  required in  the vertical  velocity (3D */
  /* fluxes) routine  and tracer  (velocity and 3D fluxes)  routine. */
  master->win_data_refill_3d(master, window, windat, master->nwindows, VELOCITY);

  clock = dp_clock();
 
  /* Set up the cell centered surface vectors and dz arrays          */
  set_dz(window, windat, wincon);

  /* Calculate the vertical velocity                                 */
  vel_w_update(window, windat, wincon);

  /* Calculate the vertical velocity alert diagnostics               */
  alerts_w(window, WVEL);

  /* Transfer vertical velocity to the master                        */
  if (master->nwindows > 1)
    master->win_data_empty_3d(master, window, windat, WVEL);

  /* Transfer alert data if required                                 */
  if (!(master->alertf & NONE))
    master_alert_fill(master, window, windat);

  /* Reset the maps for blended zones if required                    */
  blend_vel(window, windat, wincon, BLR1, NULL);

  windat->wclk += (dp_clock() - clock);

}

/* END mode3d_post_window_p2()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to update the 3D mode in each window                      */
/*-------------------------------------------------------------------*/
void mode3d_step(geometry_t *geom,    /* Global geometry             */
                 master_t *master,    /* Master data                 */
                 geometry_t **window, /* Window geometry             */
                 window_t **windat,   /* Window data                 */
                 win_priv_t **wincon, /* Window constants            */
                 int nwindows         /* Number of windows           */
  )
{
  /*-----------------------------------------------------------------*/
  /* Reconfigure routines on the windows if required                 */
  /* Note: Not supported by MPI */
  if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);
  if (master->regf & RS_TSSET) timeseries_init_w(master, window);
  if (master->regf & RS_PSSSET) sourcesink_reinit(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Reset the means on the master if required.                      */
  reset_means_m(master);

  /*-----------------------------------------------------------------*/
  /* Do the custom tracer routines on the master                     */
  TIMING_SET;
  if (master->dp_mode & DP_MPI)
    bdry_eval_tr_w(geom, master, window[master->mpi_rank+1]);
  else
    bdry_eval_tr_m(geom, master);
  TIMING_DUMP(2, "  bdry3d_tr ");

  /*-----------------------------------------------------------------*/
  /* Invoke part 1 of the velocity step - mixing, vertical spacing   */
  TIMING_SET;
  dp_vel3d_step_p1();
  TIMING_DUMP(2, "  vel3d_step_1 ");

#if TR_CK
  check_transfers(geom, window, windat, wincon, nwindows, VEL3D);
#endif

  /*-----------------------------------------------------------------*/
  /* Interpolate mixing coefficients and vertical grid spacing on    */
  /* the master for zoomed grids.                                    */
  if (master->dozoom != NOZOOM) {
    global_interp(geom, master->Vz, geom->nzin);
    global_interp(geom, master->Kz, geom->nzin);
    interp_dz_at_u1(master, geom);
    interp_dz_at_u2(master, geom);
  }

  /*-----------------------------------------------------------------*/
  /* Do the custom velocity routines on the master                   */
  TIMING_SET;
  if (master->dp_mode & DP_MPI) {
    bdry_eval_u1_w(geom, master, window[master->mpi_rank+1]);
    bdry_eval_u2_w(geom, master, window[master->mpi_rank+1]);
  } else {
    bdry_eval_u1_m(geom, master);
    bdry_eval_u2_m(geom, master);
  }
  TIMING_DUMP(2,"  bdry3d_u1u2 ");

  if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);

  /*-----------------------------------------------------------------*/
  /* Invoke part 2 of the velocity step - update velocities          */
  TIMING_SET;
  dp_vel3d_step_p2();
  TIMING_DUMP(2,"  vel3d_step_2 ");
  if (master->crf == RS_RESTART) return;

#if TR_CK
  check_transfers(geom, window, windat, wincon, nwindows, VEL3D|MIXING);
#endif
}

/* END mode3d_step()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Do the post 3D mode calculations. This includes adjusting the     */
/* vertical average of 3D velocities to equal the 2D velocities, set */
/* the lateral boundary conditions on velocity, calculate the 3D     */
/* fluxes through the cell faces and calculate vertical velocity.    */
/* The velocities in the master must be updated and each window is   */
/* filled with tracer and velocity data from the master.             */
/*-------------------------------------------------------------------*/
void mode3d_post(geometry_t *geom,     /* Global geometry            */
                 master_t *master,     /* Master data structure      */
                 geometry_t **window,  /* Window geometry            */
                 window_t **windat,    /* Window data                */
                 win_priv_t **wincon,  /* Window constants           */
                 int nwindows          /* Number of windows          */
  )
{

  /* Step forward in time at interpolated cells for zoomed grids     */
  if (master->dozoom != NOZOOM)
    global_interp_pre(master, geom);

  /*-----------------------------------------------------------------*/
  /* Adjust velocities.                                              */
  dp_vel3d_post_p1();

  if (!master->mpi_rank && master->update_master)
    master->update_master(master, windat, CFL);

  /*-----------------------------------------------------------------*/
  /* Interpolate and adjust velocity on the master for zoomed grids  */
  if (master->dozoom != NOZOOM) {
    interp_dz_at_u1(master, geom);
    interp_dz_at_u2(master, geom);
    global_interp_e1(geom, master->u1, geom->nzine1);
    global_interp_e2(geom, master->u2, geom->nzine2);
    if (!(master->dozoom & PRECOND)) {
      interp_adjust_e1(master, geom);
      interp_adjust_e2(master, geom);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the  lateral boundary conditions for velocity. This is done */
  /* at this stage for velocity so that the velocities calculated in */
  /* velocity_adjust()  above  are used  for the tangential  lateral */
  /* boundary ghost cells.                                           */
  /* NOTE : this should be also done for velocity after any velocity */
  /* initialisation.                                                 */
#if GLOB_BC
  vel2D_lbc(master->u1, geom->nbpt, geom->nbe1,
            geom->bpte1, geom->bine1, master->slip);
  vel2D_lbc(master->u2, geom->nbpt, geom->nbe2,
            geom->bpte2, geom->bine2, master->slip);
#endif

  /* Set 3D fluxes                                                   */
  dp_vel3d_post_p2();

#if TR_CK
  check_transfers(geom, window, windat, wincon, nwindows, VEL3D|FLUX);
#endif

  /*-----------------------------------------------------------------*/
  /* Interpolate vertical velocity on the master for zoomed grids    */
  if (master->dozoom != NOZOOM)
    global_interp_post(master, geom);

}

/* END mode3d_post()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set up quantities required for the 2D mode             */
/*-------------------------------------------------------------------*/
void mode3d_prep(geometry_t *geom,      /* Global geometry           */
                 master_t *master,      /* Master data               */
                 geometry_t **window,   /* Window geometry           */
                 window_t **windat,     /* Window data               */
                 win_priv_t **wincon,   /* Window constants          */
                 int nwindows           /* Number of windows         */
  )
{
  int n, nn;                    /* Counters                          */
  int c, cc, cs;                /* Sparse coordinates / counters     */
  double ramp;                  /* Ramp value                        */

  /*-----------------------------------------------------------------*/
  /* Do the custom routines on the master                            */
  bdry_eval_tr_m(geom, master);
  bdry_eval_u1_m(geom, master);
  bdry_eval_u2_m(geom, master);

  /*-----------------------------------------------------------------*/
  /* Interpolate bottom velocity on the master for zoomed grids      */
  if (master->dozoom != NOZOOM) {
    global_interp_e1(geom, master->u1bot, geom->nzinS);
    global_interp_e2(geom, master->u2bot, geom->nzinS);
  }

  /*-----------------------------------------------------------------*/
  /* Loop over all windows to calculate the vertical mixing          */
  /* coefficients and set the vertical grid spacings.                */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];
    ramp = (wincon[n]->rampf & WIND) ? windat[n]->rampval : 1.0;

    /*---------------------------------------------------------------*/
    /* Set the vertical grid spacings at e1 and e2 faces.            */
    master->win_data_fill_3d(master, window[n], windat[n], nwindows);

    /*---------------------------------------------------------------*/
    /* Calculate a heatflux if required */
    calc_heatf(window[n], windat[n], wincon[n]);

    set_dz_at_u1(window[n], windat[n], wincon[n]);
    set_dz_at_u2(window[n], windat[n], wincon[n]);
    init_sigma(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Get the e1 vertical pressure integrals for the 2D mode.       */
    /* Set pointers and initialise.                                  */
    windat[n]->nu1 = wincon[n]->w9;
    windat[n]->nu2 = wincon[n]->w10;
    memset(wincon[n]->u1inter, 0, window[n]->sgsizS * sizeof(double));
    memset(windat[n]->u1bot, 0, window[n]->sgsizS * sizeof(double));
    /* Set the maps to be self-mapping across u2 OBC's.              */
    set_map_e1(window[n]);
    /* Get the cells to process at the e1 face                       */
    cells2process_e1(window[n], windat[n], wincon[n]);
    /* Transfer any custom data                                      */
    bdry_transfer_u1(master, window[n], windat[n]);
    /* Get the vertical pressure integrals.                          */
    pressure_u1(window[n], windat[n], wincon[n]);

    /* Get the momentum flux OUT of the water due to the wind        */
    for (cc = 1; cc <= wincon[n]->vcs; cc++) {
      c = wincon[n]->s1[cc];    /* 3D cell to process                */
      cs = window[n]->m2d[c];   /* 2D cell corresponding to c        */

      wincon[n]->u1inter[cs] += (ramp * windat[n]->wind1[cs] /        
				 wincon[n]->topdensu1[cs]);

      wincon[n]->u1inter[cs] /= windat[n]->depth_e1[cs];
      wincon[n]->topdensu1[cs] *= wincon[n]->g;
      wincon[n]->densavu1[cs] *= window[n]->h1au1[cs];
    }
    memset(windat[n]->u1, 0, window[n]->sgsiz * sizeof(double));
    /* Get the open boundary conditions                              */
    bdry_u1_3d(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Get the e2 vertical pressure integrals for the 2D mode.       */
    /* Set pointers and initialise.                                  */
    memset(wincon[n]->u2inter, 0, window[n]->sgsizS * sizeof(double));
    memset(windat[n]->u2bot, 0, window[n]->sgsizS * sizeof(double));
    /* Set the maps to be self-mapping across u2 OBC's.              */
    set_map_e2(window[n]);
    /* Get the cells to process at the e2 face                       */
    cells2process_e2(window[n], windat[n], wincon[n]);
    /* Transfer any custom data                                      */
    bdry_transfer_u2(master, window[n], windat[n]);
    /* Get the vertical pressure integrals.                          */
    pressure_u2(window[n], windat[n], wincon[n]);
    /* Get the momentum flux OUT of the water due to the wind        */
    for (cc = 1; cc <= wincon[n]->vcs; cc++) {
      c = wincon[n]->s2[cc];    /* 3D cell to process                */
      cs = window[n]->m2d[c];   /* 2D cell corresponding to c        */

      wincon[n]->u2inter[cs] += (ramp * windat[n]->wind2[cs] /
				 wincon[n]->topdensu2[cs]);

      wincon[n]->u2inter[cs] /= windat[n]->depth_e2[cs];
      wincon[n]->topdensu2[cs] *= wincon[n]->g;
      wincon[n]->densavu2[cs] *= window[n]->h2au2[cs];
    }
    memset(windat[n]->u2, 0, window[n]->sgsiz * sizeof(double));
    /* Get the open boundary conditions                              */
    bdry_u2_3d(window[n], windat[n], wincon[n]);
    /* Transfer alert data if required                               */
    if (!(master->alertf & NONE))
      master_alert_fill(master, window[n], windat[n]);
  }
}

/* END mode3d_prep()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void vel_u1_update(geometry_t *window,  /* Window geometry           */
                   window_t *windat,    /* Window data               */
                   win_priv_t *wincon   /* Window constants          */
  )
{
  int c, cc, cs;                /* Sparse coodinate / counter        */
  double *depth;                /* Depth of the water column         */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  depth = wincon->d5;

  /*-----------------------------------------------------------------*/
  /* Get the cells to process at the e1 face                         */
  cells2process_e1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation          */
  precalc_u1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across u2 OBC's                 */
  set_map_e1(window);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection                                     */
  if(wincon->momsc & (ANGULAR|ANGULAR3D)) {
    if(!(windat->nstep % 2))
      advect_u1_3d_ang_flux_b(window, windat, wincon);
    else
      advect_u1_3d_ang_flux_b(window, windat, wincon);
  }
  else {
    if (advect_u1_3d(window, windat, wincon)) return;
  }
  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s1, wincon->vc, windat->nu1,
             wincon->tendency, windat->u1_adv);
  debug_c(window, D_U, D_ADVECT);

  /*-----------------------------------------------------------------*/
  /* Add the dispersion terms to u1inter. Note: u1adv is multiplied  */
  /* by the 2D sub-timestep in the 2D mode in order to integrate     */
  /* over the 2D mode, hence u1adv needs to be divided by dt here to */
  /* complete the mean.                                              */
  if (wincon->nonlinear && !(wincon->u1_f & ADVECT)) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      /* Add the dispersion terms to u1inter                         */
      wincon->u1inter[cs] =
        (wincon->u1inter[cs] + wincon->u1adv[cs]) / windat->dtf;
    }
  }

  /* Initialise the dispersion terms                                 */
  memset(wincon->u1adv, 0, window->sgsizS * sizeof(double));

  /* Set the blended zones if required                               */
  blend_vel(window, windat, wincon, U1_VH, wincon->u1vh);

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included)                 */
  wincon->hor_mix->setup(window, windat, wincon);
  wincon->hor_mix->u1(window, windat, wincon);

  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s1, wincon->vc, windat->nu1,
             wincon->tendency, windat->u1_hdif);
  debug_c(window, D_U, D_HDIFF);

  /* Reset thin layers                                               */
  reset_thin_u1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Include the pressure gradient term                              */
  pressure_u1(window, windat, wincon);
  debug_c(window, D_U, D_PRESSURE);

  /*-----------------------------------------------------------------*/
  /* Include the Coriolis term                                       */
  if (!(wincon->u1_f & CORIOLIS))
    coriolis_u1(window, windat, wincon);
  if (wincon->tendf && wincon->tendency) {
    get_tend(window, wincon->s1, wincon->vc, windat->nu1,
             wincon->tendency, windat->u1_cor);
    if (wincon->waves & STOKES_DRIFT) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	windat->u1_cor[c] -= windat->u1_sto[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion                                       */
  if (vdiff_u1(window, windat, wincon)) return;
  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s1, wincon->vc, windat->nu1,
             wincon->tendency, windat->u1_vdif);
  debug_c(window, D_U, D_VZ);

  /*-----------------------------------------------------------------*/
  /* Divide the 3D dispersion term by the depth and 3D time-step to  */
  /* get the mean over the sub-timestep and complete the 3D          */
  /* advection vertical integral. The mean of the horizontal 2D      */
  /* advection over the 2D mode is subtracted from u1adv in          */
  /* advect_u1_2d(). This involves multiplying the 2D advective terms */
  /* by the 2D sub-timestep (if sub-timestepping is performed) or    */
  /* dt2d (if sub-timestepping is not performed) then dividing by dt */
  /* (this division is performed above) to get the mean over the 2D  */
  /* mode. In order to isolate the 3D part from this 2D time         */
  /* integration, multiply by dt here: this will cancel the division */
  /* required to get the time mean of the 3D part, so no operation   */
  /* by dt is required here.                                         */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];

    /* Divide internal terms for 2d part by depth                    */
    wincon->u1inter[cs] /= max(depth[cs], wincon->hmin);
    /* Other calculations to save time                               */
    wincon->topdensu1[cs] *= wincon->g;
    wincon->densavu1[cs] *= window->h1au1[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Blend velocities over a subregion if required                   */
  blend_vel(window, windat, wincon, U13D, windat->nu1);

  /*-----------------------------------------------------------------*/
  /* Calculate u1 boundary values                                    */
  bdry_u1_3d(window, windat, wincon);
  debug_c(window, D_U, D_POST);
}

/* END vel_u1_update()                                               */
/*-------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness at the e1 faces.                */
/*-------------------------------------------------------------------*/
void set_dz_at_u1(geometry_t *window,  /* Window geometry            */
                  window_t *windat,    /* Window data                */
                  win_priv_t *wincon   /* Window constants           */
  )
{
  int cc;                       /* Counter                           */
  int c, c3;                    /* Sparse coordinate                 */
  int xm1;                      /* Sparse cell at (i-1)              */
  int zm1;                      /* Sparse cell below cell c3         */
  int cs, cb;                   /* Sparse coordinate of the bottom   */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  double *maxeta;               /* Maximum elevation at the e1 face  */

  /*-----------------------------------------------------------------*/
  /* Set the pointers                                                */
  maxeta = wincon->d5;

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->x2_e1 + 1;
    memcpy(window->sur_e1, window->w2_e1, size * sizeof(int));

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom. This is performed over all wet + auxiliary cells      */
    /* since dzu1[] is required at the auxiliary cells in the        */
    /* advection routine advect_u1_3d().                             */
    memset(windat->dzu1, 0, window->sgsiz * sizeof(double));
    memset(maxeta, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->b2_e1; cc++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      c = c3 = window->sur_e1[cc];  /* 3D surface coordinate         */
      cs = window->m2d[c];          /* 2D cell corresponding to c    */
      cb = window->bot_e1[cc];      /* 3D bottom coordinate          */
      xm1 = window->xm1[cs];        /* Cell to the west of cs        */

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = windat->topz[cs];
      while (c != cb) {
        bot = window->gridz[c];
        windat->dzu1[c] = top - bot;
        top = bot;
        c = window->zm1[c];
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu1[cb] = top - window->botzu1[cs];
    }
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Find the 3D sparse coordinate corresponding to the free       */
    /* surface.                                                      */
    memset(windat->sur_e1, 0, window->sgsizS * sizeof(int));
    for (cc = 1; cc <= window->b2_e1; cc++) {

      /* Get the 2D and 3D sparse coordinates                        */
      cs = c = window->w2_e1[cc]; /* 2D coordinate                   */
      xm1 = window->xm1[cs];      /* Cell to the west of cs          */

      /* Get the new sparse location of the surface. The surface at  */
      /* the e1 face is defined as the higher of the cell centered   */
      /* elevations bordering the face. For elevations below the     */
      /* bottom this is set to the sediment coordinate.              */
      top = maxeta[cs] = max(windat->eta[cs], windat->eta[xm1]);
      top = maxeta[cs] = max(top, window->botzu1[cs]);
      zm1 = window->zm1[c];

      /* Loop down the water column                                  */
      while (c != zm1 && window->gridz[c] > top) {
        c = zm1;
        zm1 = window->zm1[c];
      }
      window->sur_e1[cc] = c;
      /* windat->sur_e1 is the k level of the surface and is         */
      /* tranferred across windows. Used to get the surface layer    */
      /* in auxiliary cells (note window->sur_e1 is only defined     */
      /* over cells b2_e1).                                          */
      windat->sur_e1[cs] = window->s2k[c];
    }

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom. This is performed over all wet + auxiliary cells      */
    /* since dzu1[] is required at the auxiliary cells in the        */
    /* advection routine advect_u1_3d().                             */
    memset(windat->dzu1, 0, window->sgsiz * sizeof(double));
    wincon->ncdry_e1 = 0;
    for (cc = 1; cc <= window->b2_e1; cc++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      c = c3 = window->sur_e1[cc];  /* 3D surface coordinate         */
      cs = window->m2d[c];          /* 2D cell corresponding to c    */
      cb = window->bot_e1[cc];      /* 3D bottom coordinate          */

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = maxeta[cs];

      while (c != cb) {
        bot = window->gridz[c];
        windat->dzu1[c] = top - bot;
        top = bot;
        c = window->zm1[c];
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu1[cb] = top - window->botzu1[cs];

      /* Save the locations of dry water columns                     */
      if (windat->dzu1[cb] == 0.0) {
	wincon->ncdry_e1++;
	wincon->cdry_e1[wincon->ncdry_e1] = cs;
      }

      /* Set the cell thickness above the old surface equal to zero. */
      /* (Note : a no-gradient condition may be required for higher  */
      /* order vertical advection schemes.                           */
      c = c3;
      top = windat->dzu1[c];
      while (c != cs) {
        c = window->zp1[c];
        windat->dzu1[c] = top;
        /*windat->dzu1[c] = 0.0;*/
      }
    }
  }
}

/* END set_dz_at_u1()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set a value in dry cells above the surface             */
/*-------------------------------------------------------------------*/
void set_surf_cond(geometry_t *window,  /* Window geometry            */
		   double *a,
		   double val,
		   int *vec,
		   int nvec)
{
  int *csur = window->wincon->i7;
  int c, cs;
  double sval;

  for (cs = 1; cs <= window->enonS; cs++) {
    if ((c = csur[cs])) {
      sval = (val < 0) ? a[c] : val;    
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	a[c] = sval;
      }
    }
  }
}


/* END set_surf_cond()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the surface coordinate over wet, auxiliary and     */
/* thin cells. Note csur = wincon->i7 uses the 2D sparse coordinate  */
/* as the index.                                                     */
/*-------------------------------------------------------------------*/
void set_surf_cells(geometry_t *window,  /* Window geometry          */
		    window_t *windat,    /* Window data              */
		    win_priv_t *wincon,  /* Window constants         */
		    int mode)
{
  int c, cc, cs, ks;
  int *csur = wincon->i7;
  int *cells, *acells, vc, *kth, nkth, *bpt, *bin, nbpt;

  /* Set pointers                                                    */
  if (mode) {
    vc = window->b2_e1;
    cells = window->sur_e1;
    acells = windat->sur_e1;
    nkth = wincon->nkth_e1;
    kth = wincon->kth_e1;
    bpt = window->bpte1S;
    bin = window->bine1S;
    nbpt = window->nbpte1S;
  } else {
    vc = window->b2_e2;
    cells = window->sur_e2;
    acells = windat->sur_e2;
    nkth = wincon->nkth_e2;
    kth = wincon->kth_e2;
    bpt = window->bpte1S;
    bin = window->bine1S;
    nbpt = window->nbpte1S;
  }

  /* Initialise and get wet window cells (OBC cells included)        */
  memset(csur, 0, window->sgsizS * sizeof(int));
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    csur[cs] = c;
  }

  /* Get auxiliary cells                                             */
  if (window->nwindows > 1) {
    for (cs = 1; cs <= window->enonS; cs++) {
      ks = acells[cs];
      if (ks) {
	c = get_local_sur(window, cs, ks);
	csur[cs] = c;
      }
    }
  }

  /* Get cells reset at thin layers                                  */
  if (wincon->thin_merge) {
    for (cc = 1; cc <= nkth; cc++) {
      c = window->zm1[kth[cc]];
      cs = window->m2d[c];
      csur[cs] = c;
    }
  }

  /* Get ghost cells                                                 */
  for (cc = 1; cc <= nbpt; cc++) {
    cs = bin[cc];
    c = bpt[cc];
    csur[c] = csur[cs];
  }
}

/* END set_surf_cells()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to exclude boundary cells from the cells to process so as */
/* to linearize the momentum equations on the boundary or within N   */
/* cells from the boundary.                                          */
/*-------------------------------------------------------------------*/
void linear_bdry_cell(
		 geometry_t *window,    /* Window geometry           */
		 window_t *windat,      /* Window data               */
		 win_priv_t *wincon,    /* Window constants          */
		 int *ctp,              /* Cells to process          */
		 int nctp,              /* Number of surface ctp     */
		 int *mask,             /* linear cell mask          */
		 int *bottom            /* Bottom cell vector        */
  )
{
  int c, cc, cb;
  int vc, vcs;
  int *cells = wincon->s4;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memcpy(cells, ctp, window->sgsizS* sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the surface cells to process vector                         */
  vcs = 1;
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    if(!mask[c]) {
      cells[vcs] = c;
      vcs++;
    }
  }
  vcs--;
  wincon->aclS = vcs;

  /*-----------------------------------------------------------------*/
  /* Return if 2d cells only are processed                           */
  if (bottom == NULL) {
    /*memcpy(wincon->s3, cells, window->sgsizS * sizeof(int));*/
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Set the sub-surface cells to process vector                     */
  vc = vcs + 1;
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    cb = bottom[cc];
    if(!mask[c]) {
      c = window->zm1[c];
      while (c <= cb) {
	cells[vc] = c;
	vc++;
	c = window->zm1[c];
      }
    }
  }

  vc = vcs + 1;
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    cb = bottom[cc];
    if(!mask[c]) {
      c = window->zm1[c];
      while (c <= cb) {
	cells[vc] = c;
	vc++;
	c = window->zm1[c];
      }
    }
  }
  vc--;
  /*
  memcpy(wincon->s3, cells, window->sgsiz * sizeof(int));
  */
  wincon->acl = vc;
}

/* END linear_bdry_cell()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the velocity over thin layers. The non-linear      */
/* terms are calculated using merged layers after which the thin     */
/* layer is set to the layer below. Pressure and Coriolis are        */
/* calculated normally, then the thin layer is again merged for the  */
/* vertical diffusion calculation.                                   */
/*-------------------------------------------------------------------*/
void set_thin_u1(geometry_t *window,    /* Window geometry           */
		 window_t *windat,      /* Window data               */
		 win_priv_t *wincon     /* Window constants          */
  )
{
  int cc, bn;                   /* Counter                           */
  int c, cs;                    /* Sparse coordinate                 */
  int zm1;                      /* Sparse cell below cell c          */
  int cb;                       /* Bottom sparse coordinate          */
  int *cells;                   /* Cells to process                  */
  int *bottom;                  /* Bottom cells                      */
  double *u1o = wincon->d7;     /* Original u1 at zm1                */
  open_bdrys_t **open = window->open; /* Window OBC structure        */
  int ks;

  /*-----------------------------------------------------------------*/
  /* Set the thin layer vector and reset the surface layer.          */
  /* Note : if the loop is done to a2_e1 then dzu1=0 at the first    */
  /* timestep for multiple windows in auxiliary cells since it has   */
  /* not yet been copied from the master.                            */
  /* Thin layers are not adjusted for SIGMA.                         */
  cells = wincon->s3;
  bottom = wincon->i5;
  wincon->nkth_e1 = 1;

  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cb = bottom[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    if (c < cb && windat->dzu1[zm1] && windat->dzu1[c] < wincon->hmin) {
      cells[cc] = zm1;
      /* Save the sparse location of the thin layer                  */
      wincon->kth_e1[wincon->nkth_e1] = c;
      wincon->nkth_e1++;
      u1o[window->m2d[c]] = windat->u1[zm1];	
      merge_thin_layers(c, zm1, windat->u1, windat->dzu1);
    } else
      cells[cc] = c;
    windat->sur_e1[cs] = 0;
  }

  /* Set thin layers for open boundary cells. These are used in the  */
  /* non-linear terms.                                               */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	zm1 = window->zm1[c];
	ks = windat->sur_e1[cs];
	if (windat->dzu1[zm1] && windat->dzu1[c] < wincon->hmin) {
	  merge_thin_layers(c, zm1, windat->u1, windat->dzu1);
	}
	windat->sur_e1[cs] = 0;
      }
    }
  }

  /* Set thin layers for auxiliary cells. These are used in the      */
  /* non-linear terms.                                               */
  if (window->nwindows > 1) {
    for (cs = 1; cs <= window->enonS; cs++) {
      ks = windat->sur_e1[cs];
      if (ks) {
	c = get_local_sur(window, cs, ks);
	zm1 = window->zm1[c];
	if (windat->dzu1[zm1] && windat->dzu1[c] < wincon->hmin) {
	  wincon->kth_e1[wincon->nkth_e1] = c;
	  wincon->nkth_e1++;
	  u1o[cs] = windat->u1[zm1];
	  merge_thin_layers(c, zm1, windat->u1, windat->dzu1);
	}
      }
    }
  }
  wincon->nkth_e1--;

  /* Set a no-gradient of the merged velocity above the surface      */
  for (cc = 1; cc <= wincon->nkth_e1; cc++) {
    c = wincon->kth_e1[cc];
    zm1 = window->zm1[c];
    for(cs = window->m2d[c]; cs != c; cs = window->zm1[cs]) {
      windat->u1[cs] = windat->u1[zm1];
    }
  }

  /* Fill the cells to process with sub-surface cells. If no thin    */
  /* cells were encountered make a direct copy of s1 to save time.   */
  if (wincon->nkth_e1 == 0) {
    memcpy(wincon->s3, wincon->s1, window->sgsiz * sizeof(int));
    wincon->ncl = wincon->vc;
  } else {
    wincon->ncl = wincon->vcs + 1;
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->zm1[cells[cc]];
      cb = bottom[cc];
      while (c <= cb) {
        cells[wincon->ncl] = c;
        wincon->ncl++;
        c = window->zm1[c];
      }
    }
    wincon->ncl--;
  }
}

/* END set_thin_u1()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset u1 velocity in thin layers                       */
/*-------------------------------------------------------------------*/
void reset_thin_u1(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
  )
{
  int cc;                       /* Counter                           */
  int c,zm1;                    /* Sparse coordinate                 */
  double *u1o = wincon->d7;     /* Original u1 at zm1                */

  /* Set a uniform velocity over thin layers if required */
  if (wincon->thin_merge) {
    for (cc = 1; cc <= wincon->nkth_e1; cc++) {
      c = wincon->kth_e1[cc];
      zm1 = window->zm1[c];
      windat->nu1[c] = windat->nu1[zm1];
      windat->dzu1[zm1] -= windat->dzu1[c];
      /* Reset u1 to the un-merged velocity                          */
      windat->u1[c] = (windat->dzu1[c]) ? windat->u1[zm1] + 
	               windat->dzu1[zm1] * 
 	              (windat->u1[zm1] - u1o[window->m2d[c]]) / 
 	               windat->dzu1[c] : 0.0;
      windat->u1[zm1] = u1o[window->m2d[c]];
    }
  }
}

/* END reset_thin_u1()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset wtop velocity in thin layers                     */
/*-------------------------------------------------------------------*/
void reset_thin_wtop(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon,  /* Window constants        */
		     int *ctp,            /* Cells to process vector */
		     int vcs,             /* Limit of ctp            */
		     double *wtop         /* Surface vertical vel    */
  )
{
  int cc;                       /* Counter                           */
  int c, cs;                    /* Sparse coordinates                */
  int xm1, xp1, ym1, yp1;       /* Spatial sparse coordinates        */
  double eta_l;                 /* Surface elevation at i-1          */
  double eta_r;                 /* Surface elevation at i+1          */
  double eta_b;                 /* Surface elevation at j-1          */
  double eta_f;                 /* Surface elevation at j+1          */

  memcpy(wtop, windat->wtop, window->sgsizS * sizeof(double));

  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[cs];
    xm1 = window->xm1[cs];
    yp1 = window->yp1[cs];
    ym1 = window->ym1[cs];

    /* Calculate the surface vertical velocity                       */
    /* Note : a no gradient condition exists above the surface for   */
    /* velocity, hence cells at cs will always contain a copy of the */
    /* surface velocity.                                             */
    eta_l = windat->eta[xm1];
    eta_r = windat->eta[xp1];
    eta_b = windat->eta[ym1];
    eta_f = windat->eta[yp1];

    wtop[cs] = windat->detadt[cs] +
      (windat->u1[cs] + windat->u1[xp1]) * (eta_r - eta_l) /
      (4.0 * window->h1acell[cs]) +
      (windat->u2[cs] + windat->u2[yp1]) * (eta_f - eta_b) /
      (4.0 * window->h2acell[cs]);

  }
}
/* END reset_thin_wtop()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the sparse coordinates of cells to process at the  */
/* e1 faces.                                                         */
/*-------------------------------------------------------------------*/
void cells2process_e1(geometry_t *window, /* Window geometry         */
                      window_t *windat,   /* Window data             */
                      win_priv_t *wincon  /* Window constants        */
  )
{
  int cc;                       /* Counter                           */
  int c, c2;                    /* Sparse coordinate                 */
  int xm1;                      /* Sparse cell at (i-1)              */
  int cs, cb;                   /* Sparse coordinate of the bottom   */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  double *maxeta;               /* Maximum elevation at the e1 face  */
  int vc, vcs;                  /* Counter                           */

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->x2_e1 + 1;
    wincon->vc = wincon->vc1 = window->v3_e1;
    wincon->vcs = wincon->vcs1 = window->v2_e1;
    wincon->vca = wincon->vca1 = 0;
    memcpy(wincon->s1, window->w3_e1, (window->a3_e1 + 1) * sizeof(int));
    memcpy(wincon->i1, window->sur_e1, size * sizeof(int));
    memcpy(wincon->i2, window->sur_e1, size * sizeof(int));
    memcpy(wincon->i3, window->sur_e1, size * sizeof(int));
    memcpy(wincon->i4, window->sur_e1, size * sizeof(int));
    memcpy(wincon->i5, window->bot_e1, size * sizeof(int));
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Set the pointers                                              */
    int *partial = wincon->i1;
    int *at_ele = wincon->i2;  /* 3D cell containing the elevation   */
    int *hi_ele = wincon->i3;  
    int *lo_ele = wincon->i4;
    int *bottom = wincon->i5;
    maxeta = wincon->d5;
    vcs = (wincon->dobdry_u1) ? window->b2_e1 : window->v2_e1;

    /* Set the sub-surface explicit map if required                  */
    if (window->sm_e1) {
      for (cc = 1; cc <= window->b2_e1; cc++) {
	cs = window->w2_e1[cc];
	c = window->sm_e1[cs];
	if (c) {
	  at_ele[cc] = window->sur_e1[cc];
	  window->sur_e1[cc] = c;
	  maxeta[cs] = min(maxeta[cs], window->gridz[window->zp1[c]]);
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the surface wet cells to process vector for u1 velocity   */
    /* and store in the buffer wincon->s1. Do not include cells      */
    /* where the free surface is below the bottom (i.e. correspond   */
    /* to sediment sparse locations).                                */
    vc = 1;
    wincon->ncbot_e1 = 0;
    for (cc = 1; cc <= vcs; cc++) {
      c = window->sur_e1[cc];   /* Surface coordinate                */
      cs = window->w2_e1[cc];   /* 2D coordinate                     */
      xm1 = window->xm1[cs];    /* Cell to the west of cs            */
      top = maxeta[cs];         /* Highest elevation at the face     */
      bot = DRY_FRAC * wincon->hmin + window->botzu1[cs];
      if (top > bot) {
        wincon->s1[vc] = c;
        cb = bottom[vc] = window->bot_e1[cc];
	/* Save cells which are one cell deep                        */
	if(c == cb) {
	  wincon->ncbot_e1++;
	  wincon->cbot_e1[wincon->ncbot_e1] = c;
	}
	/*
        lo_ele[vc] = c;
        if (top > windat->eta[xm1])
          lo_ele[vc] = window->xm1[c];
	*/
        vc++;
      }
    }
    wincon->vcs = vc - 1;

    /*---------------------------------------------------------------*/
    /* Loop from the layer below the surface to the bottom and get   */
    /* the cells to process. The are arranged so that cells 1 to     */
    /* wincon->vcs contain surface cells, cells wincon->vcs+1 to     */
    /* wincon->vca contain cells where one side of the cell face is  */
    /* dry and cells wincon->vca+1 to wincon->vc contain cells where */
    /* both sides of the cell face is wet.                           */
    /* First get the cells where one side of the face is dry. The    */
    /* last cell that falls into this category is temporarily saved  */
    /* in partial[] for use in the next loop.                        */
    for (cc = 1; cc <= wincon->vcs; cc++) {

      c = c2 = wincon->s1[cc];  /* 3D surface cell                   */
      cs = window->m2d[c];      /* 2D cell corresponding to c        */
      cb = bottom[cc];          /* 3D bottom cell                    */
      xm1 = window->xm1[cs];    /* Cell to the west of cs            */
      at_ele[cc] = c;

      top = min(windat->eta[cs], windat->eta[xm1]);

      /* Cells where one face is completely dry are added to s1.     */
      /* The cell containing the surface is stored in at_ele, i.e    */
      /* if eta[c] > eta[xm1] at_ele[cc] contains the higher cell    */
      /* location, s1[cc]. If eta[xm1] > eta[c] at_ele[cc] contains  */
      /* the lower cell location which will be the same as           */
      /* partial[cc].                                                */
      while (c < cb && (bot = max(window->gridz[c], 
				  window->botzu1[cs])) > top) {

        if (c != c2) {
          wincon->s1[vc] = c;
          vc++;
        }
        if (bot > windat->eta[cs])
          at_ele[cc] = window->zm1[at_ele[cc]];
        c = window->zm1[c];
      }

      /* Get the cell location where the surface resides. This is    */
      /* current value of c after the loop above (i.e when gridz     */
      /* becomes less than eta). Add this cell to s1 if not already  */
      /* accounted for. hi_ele = partial and lo_ele = xm1 if eta[c]  */
      /* > eta[xm1], and hi_ele = xm1 and lo_ele = partial if        */
      /* eta[xm1] > eta[c].                                          */
      if (c <= cb) {
        if (c != c2) {
          wincon->s1[vc] = c;
          vc++;
        }
        partial[cc] = hi_ele[cc] = c;
        lo_ele[cc] = window->xm1[c];
	/* Swap if eta[xm1] > eta[c]                                 */
        if (top < windat->eta[xm1]) {
          lo_ele[cc] = c;
          hi_ele[cc] = window->xm1[c];
        }
      }
    }
    wincon->vca = vc - 1;

    /* Next get the cells where both sides of the face is wet.       */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->zm1[partial[cc]];
      cb = bottom[cc];
      while (c <= cb) {
        wincon->s1[vc] = c;
        vc++;
        c = window->zm1[c];
      }
    }
    wincon->vc = vc - 1;
    /* Remember the number of u1 cells to process for the 2D mode    */
    wincon->vc1 = wincon->vc;
    wincon->vcs1 = wincon->vcs;
    wincon->vca1 = wincon->vca;
  }
}

/* END cells2process_e1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Vertical diffusion                                                */
/*-------------------------------------------------------------------*/
int vdiff_u1(geometry_t *window,  /* Window geometry                */
	     window_t *windat,    /* Window data                    */
	     win_priv_t *wincon   /* Window constants               */
  )
{
  int c, cc;               /* Sparse coordinates / counters          */
  int cs, cb;              /* Surface / bottom sparse coordinate     */
  int xm1;                 /* Cell location to the west of c         */
  int *at_ele;             /* Array of fully wet cell coordinates    */
  int *partial;            /* Coordinate of cell containing eta[c]   */
  int *bottom;             /* Bottom coordinate for cells to process */
  double *Vzav;            /* Viscosity at the e1 face               */
  double *f_top;           /* Flux out of the water at the surface   */
  double *f_bot;           /* Flux into the water at the bottom      */
  double *u2au1;           /* u2 velocity at e1 face                 */
  double Cdu1;             /* Bottom drag at the e1 face             */
  double botdz;            /* Thickness of the bottom layer          */
  double val;              /* Bottom current speed                   */
  double *w;               /* Vertical velocity                      */

  int yp1, xmyp1;
  double ramp = (wincon->rampf & WIND) ? windat->rampval : 1.0;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* s1 = wet cells to process for u1 velocity                       */
  partial = wincon->i1;
  bottom = wincon->i5;
  at_ele = wincon->i2;
  f_top = wincon->d1;
  f_bot = wincon->d2;
  Vzav = wincon->w1;
  u2au1 = wincon->w7;
  w = wincon->w8;

  /*-----------------------------------------------------------------*/
  /* Set the vertical viscosity at the e1 face. Where one side of    */
  /* the water column is higher, set the viscosity at the cell face  */
  /* equal to the viscosity of the higher water column. If the       */
  /* current water column, cs, is the higher then no action is taken */
  /* (c == cs), if the water column to the west of face cs is higher */
  /* then the viscosity at c is set to the viscosity at xm1. This is */
  /* only performed to the layer above that containing the free      */
  /* surface in column cs; below this the mean is taken.             */
  memcpy(Vzav, windat->Vz, window->sgsiz * sizeof(double));
  /* Set the viscosity on completely dry faces                       */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cs = at_ele[cc];
    while (c < cs) {
      xm1 = window->xm1[c];
      Vzav[c] = windat->Vz[xm1];
      c = window->zm1[c];
    }
  }

  /* Set the viscosity for cells where the face is partially dry     */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = partial[cc];
    xm1 = window->xm1[c];
    Vzav[c] = 0.5 * (windat->Vz[c] + windat->Vz[xm1]);
  }
  /* Set the viscosity for cells where the face completely wet       */
  for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    xm1 = window->xm1[c];
    Vzav[c] = 0.5 * (windat->Vz[c] + windat->Vz[xm1]);
  }

  /*-----------------------------------------------------------------*/
  /* Set the vertical velocity                                       */
  memset(w, 0, window->sgsiz * sizeof(double));
  if (wincon->momsc & WIMPLICIT) {
    int zm1;
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];         /* 3D cell to process              */
      xm1 = window->xm1[c];
      zm1 = window->zm1[c];
      w[zm1] = 0.5 * (windat->w[c] + windat->w[xm1]);
    }
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];       /* 3D cell to process                */
      cs = window->m2d[c];      /* 2D cell corresponding to c        */
      cb = bottom[cc];          /* 3D bottom coordinate              */
      zm1 = window->zm1[cb];
      xm1 = window->xm1[cs];
      /* Bottom vertical velocity                                    */
      w[zm1] = 0.5 * (windat->wbot[xm1] + windat->wbot[cs]);
      /* Surface vertical velocity                                   */
      w[c] = 0.5 * (windat->wtop[xm1] + windat->wtop[cs]);
    }
    set_wvel(window, windat, wincon, w, 0);
  }

  /*-----------------------------------------------------------------*/
  /* Set the momentum flux out of the water at the top (due to the   */
  /* wind) and flux into the water at the bottom (due to friction).  */
  memset(f_top, 0, window->sgsizS * sizeof(double));
  memset(f_bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];         /* 3D cell to process                */
    cs = window->m2d[c];        /* 2D cell corresponding to c        */
    cb = bottom[cc];            /* 3D bottom coordinate              */
    xm1 = window->xm1[cb];
    yp1 = window->yp1[cb];
    xmyp1 = window->xmyp1[cb];

    /* Surface flux                                                  */
    f_top[cs] = -ramp * windat->wind1[cs] / wincon->topdensu1[cs];      

    /* Add surface stress to 2d forcing term.                        */
    wincon->u1inter[cs] -= f_top[cs];

    /* Bottom flux */
    u2au1[cb] = 0.25 * (windat->u2b[xm1] + windat->u2b[xmyp1] +
                        windat->u2b[cb] + windat->u2b[yp1]);
    Cdu1 = 0.5 * (wincon->Cd[cs] + wincon->Cd[window->xm1[cs]]);
    val = sqrt(windat->u1b[cb] * windat->u1b[cb] + u2au1[cb] * u2au1[cb]);
    val = Cdu1 * max(wincon->uf, val);
    botdz = max(wincon->hmin, windat->dzu1[cb] * wincon->Ds[cs]);
    /* Quadratic bottom friction, truncated to ensure stability      */
    if (val > botdz / windat->dt)
      val = botdz / windat->dt;
    f_bot[cs] = -val * windat->u1b[cb];
  }
  if (wincon->vorticity & TENDENCY)
    memcpy(windat->rv_bsc, f_bot, window->sgsizS * sizeof(double));
  if (wincon->numbers & BOTSTRESS)
    memcpy(windat->tau_be1, f_bot, window->sgsizS * sizeof(double));
  if (wincon->numbers & EKPUMP)
    ekman_pump_e1(window, windat, wincon, windat->wind1, f_bot);

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion.                                      */
  /* SIGMA : multiply by depth to give D.u1 at t+1. Note: should be  */
  /* dividing by depth at t+1 but this is not possible since         */
  /* elevation is only solved after 3D velocity.                     */
  /* Note : the velocity at time t-1 is passed to the vertical       */
  /* diffusion routine. This appears to be more stable than passing  */
  /* the velocity at time t+1 (see Science manual).                  */
  if (!(wincon->u1_f & VDIFF)) {
    if (implicit_vdiff(window, windat, wincon, windat->u1b, windat->nu1, Vzav,
		       windat->dzu1, f_bot, f_top, wincon->s1, wincon->i5,
		       wincon->vcs, wincon->Ds, wincon->Ds, w)) return(1);

  }
  return(0);
}

/* END vdiff_u1()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to add the Coriolis term to the u1 velocity               */
/*-------------------------------------------------------------------*/
void coriolis_u1(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc, cs;                /* Sparse coordinates / counters     */
  double *u2au1;                /* u2 velocity at e1 face            */
  double *midx;                 /* Total depth at e1 face            */

  /* Set pointers */
  midx = wincon->d4;            /* Set in pressure_u1()              */
  u2au1 = wincon->w7;           /* Set in precalc_u1()               */

  /*-----------------------------------------------------------------*/
  /* Add the Coriolis term to the u1 velocity                        */
  /* SIGMA : multiply by total depth at u1                           */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    windat->nu1[c] += windat->dt * wincon->u1c5[cs] * u2au1[c] * midx[cs];
  }
  if(wincon->waves & STOKES_DRIFT) {
    double sdc, vf;
    double *tend;
    /* Stokes drift                                                  */
    tend = (wincon->tendf) ? windat->u1_sto : wincon->w2;
    memset(tend, 0.0, window->sgsiz * sizeof(double));
    get_sdc_e2(window, windat, wincon);
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      sdc = wincon->u1c5[cs] * wincon->w1[c] * midx[cs];
      windat->nu1[c] += windat->dt * sdc;
      wincon->u1inter[cs] += sdc * windat->dzu1[c];
      tend[c] += windat->dt * sdc;
    }
    /* Vortex force                                                  */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      sdc = wincon->w1[c] * midx[cs] * 
	((u2au1[window->xp1[c]] - u2au1[window->xm1[c]]) / (2.0 * window->h1au1[cs]) -
	 (windat->u1[window->yp1[c]] - windat->u1[window->ym1[c]]) / (2.0 * window->h2au1[cs]));
      windat->nu1[c] += windat->dt * sdc;
      wincon->u1inter[cs] += sdc * windat->dzu1[c];
      tend[c] += windat->dt * sdc;
    }
  }
}

/* END coriolis_u1()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U1 horizontal pressure gradient term. Here we can assume that     */
/* the u1 column is not U1OUTSIDE, U1SOLID, or U1BDRY, and is not    */
/* totally dry.                                                      */
/*                                                                   */
/* If eta is treated in linear fashion, then the term looks like     */
/*                                               _                   */
/*                  dP             deta         |  dro               */
/*                  --   =  g*ro0* ----  +   g* |  --- dz            */
/*                  de1             de1         |  de1               */
/*                                             -                     */
/* where ro0 is the surface density                                  */
/*             ((dens[nz-1][j][cl]+dens[nz-1][j][cr])/2)             */
/* and the integral is from z=0.0 down to the bottom of the layer    */
/* currently being calculated.                                       */
/*                                                                   */
/* If eta is nonlinear, the water surface may be in different        */
/* layers on either side of the u1 column. Here, the pressure        */
/* gradient term is:                                                 */
/*                                _                                  */
/*                  dP           |  dro                              */
/*                  --   =    g* |  --- dz                           */
/*                  de1          |  de1                              */
/*                              -                                    */
/* where the integral is from the higher water surface               */
/* (max(topz[j][cl],topz[j][cr])) to the bottom of the layer         */
/* currently being calculated.                                       */
/*-------------------------------------------------------------------*/
void pressure_u1(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc;              /* Local sparse coordinate / counter       */
  int cs;                 /* 2D local sparse coordinate              */
  int cl;                 /* Partially wet partial cell coordinate   */
  int ch;                 /* Fully wet partial cell coordinate       */
  int zm1;                /* Cell below cell c                       */
  int zp1, xmzp1;         /* Cell above and left-above cell c        */
  int xm1, xm1s;          /* Cell location to the west of c          */
  int *partial;           /* Coordinate of cell containing eta[c]    */
  int *at_ele;            /* Array of fully wet cell coordinates     */
  int *hi_ele;            /* Array of fully wet cell coordinates     */
  int *lo_ele;            /* Array of partially wet cell coordinates */
  double *dhdiv;          /* Vertically integrated depth (hlower)    */
  /* double *dhdiu; *//* Vertically integrated depth (hupper)        */
  double *dd;             /* Vertically integrated pressure gradient */
  double *dinter;         /* Pressure gradient for the 2D mode       */
  double *u1inter;        /* Vertical integrals for 2d mode          */
  double *hidens;         /* Density corresponding to coordinate ch  */
  double *dav;            /* Mean density at the cell face           */
  double *dzu1;           /* Cell theickness ( = windat->dzu1)       */
  double *hupper;         /* Partial cell upper thickness work array */
  double *hlower;         /* Partial cell lower thickness work array */
  double *dens;           /* Density ( = windat->dens)               */
  double *dzsum;          /* Sum of dz for dry / partial layers      */
  double *mask;           /* Mask out all fully wet or dry faces     */
  double *btp;            /* Barotropic pressure tendency            */
  double top;             /* Top cell thickness                      */
  double bot;             /* Bottom cell thickness                   */
  double *midx;           /* Total depth at the cell face            */
  double d1;              /* Dummy                                   */
  double tol = 1e-10;     /* Tolerence value                         */
  double *corr;           /* Baroclinic correction for sigma         */
  double *ddbc;           /* Baroclinic contribution to u1 velocity  */
  double *ddbt;           /* Barotropic contribution to u1 velocity  */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  hidens = wincon->w1;
  hupper = wincon->w2;
  hlower = wincon->w3;
  dav = wincon->w4;
  dzsum = wincon->w5;
  mask = wincon->w6;
  corr = wincon->w8;
  btp = wincon->w3;

  dd = wincon->d1;
  dhdiv = wincon->d2;
  u1inter = wincon->d3;
  midx = wincon->d4;
  /* dhdiu=wincon->d7; */

  dens = windat->dens;
  dzu1 = windat->dzu1;

  partial = wincon->i1;
  at_ele = wincon->i2;
  hi_ele = wincon->i3;
  lo_ele = wincon->i4;

  if (wincon->u1_f & PRESS_BC)
    ddbc = wincon->d7;
  else
    ddbc = dd;
  if (wincon->u1_f & PRESS_BT)
    ddbt = wincon->d7;
  else
    ddbt = dd;
  if (wincon->u1av_f & PRESS_BC)
    dinter = wincon->d7;
  else
    dinter = u1inter;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->topdensu1, 0, window->sgsizS * sizeof(double));
  memset(wincon->densavu1, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Linear eta case; calculate surface pressure gradient due to eta */
  /* slope and atmospheric pressure gradient.                        */
  if (!wincon->nonlinear) {
    double *surgrad = wincon->d1;

    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      xm1 = window->xm1[c];
      cs = window->m2d[c];
      xm1s = window->xm1[cs];

      midx[cs] = wincon->mdx[cs];
      wincon->topdensu1[cs] = 0.5 * (windat->dens[xm1] + windat->dens[c]);
      surgrad[cs] = wincon->g * wincon->topdensu1[cs] *
        (windat->eta[cs] - windat->eta[xm1s]) +
        (windat->patm[cs] - windat->patm[xm1s]);
      dd[cs] = 0.0;
    }
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      xm1 = window->xm1[c];
      cs = window->m2d[c];
      xm1s = window->xm1[cs];

      /* Density gradient integral                                   */
      dav[c] = 0.5 * (windat->dens[c] + windat->dens[xm1]);
      dd[cs] +=
        wincon->g * (windat->dens[c] -
                     windat->dens[xm1]) * windat->dzu1[c];

      /* Add term to new u1 value                                    */
      windat->nu1[c] +=
        windat->dt * wincon->u1c6[cs] * (surgrad[cs] + dd[cs]) / dav[c];

      /* Integrate internal term for 2d part                         */
      wincon->u1inter[cs] +=
        wincon->u1c6[cs] * dd[cs] * windat->dzu1[c] / dav[c];
      wincon->densavu1[cs] += dav[c] * windat->dzu1[c];
    }
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      xm1s = window->xm1[cs];
      top = max(windat->eta[xm1s], windat->eta[cs]);
      bot = window->botzu1[cs];
      wincon->densavu1[cs] /= (top - bot);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Non-linear eta case                                             */
  else {

    /*---------------------------------------------------------------*/
    /* Precondition work arrays for the pressure gradient integral.  */
    /* These are established for the surface cells (i.e. cells that  */
    /* have one face completely or partially dry) as follows:        */
    /*                                                               */
    /* work array   | completely dry face  | partially dry face      */
    /* -------------|----------------------|----------------------   */
    /* hidens       | rho at wet cell      | rho at full wet cell    */
    /* hupper       | cell thickness, dzu1 | upper thickness, hu     */
    /* hlower       | cell thickness, dzu1 | lower thickness, hl     */
    /* dav          | rho at wet cell      | 0.5(rho(c)+rho(xm1))    */
    /* dzsum        | rho(c) or -rho(xm1)  | rho(c) or -rho(xm1)     */
    /* mask         | zero                 | one                     */

    /* Initialise the work arrays                                    */
    memcpy(hidens, dens, window->sgsiz * sizeof(double));
    memcpy(hupper, dzu1, window->sgsiz * sizeof(double));
    memcpy(hlower, dzu1, window->sgsiz * sizeof(double));
    memcpy(dav, dens, window->sgsiz * sizeof(double));
    /*memset(dzsum, 0, window->sgsiz * sizeof(double));*/
    memcpy(dzsum, dens, window->sgsiz * sizeof(double));
    memset(mask, 0, window->sgsiz * sizeof(double));

    /* Set the work arrays for cells where the face is completely    */
    /* dry. Where one side of the water column is higher, set the    */
    /* density at the cell face equal to the density of the higher   */
    /* water column. If the current water column, cs, is the higher  */
    /* then no action is taken (c == cs), if the water column to the */
    /* west of face cs is higher then the density at c is set to the */
    /* density at xm1. This is only performed to the layer above     */
    /* that containing the free surface in column cs.                */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = at_ele[cc];
      midx[window->m2d[cs]] =
        0.5 * (wincon->Ds[window->m2d[cs]] + 
	       wincon->Ds[window->xm1[window->m2d[cs]]]);
      while (c < cs) {
        xm1 = window->xm1[c];
        hidens[c] = dav[c] = dens[xm1];
	dzsum[c] = -dens[xm1];
        c = window->zm1[c];
      }
    }
    if (!wincon->sigma)
      memcpy(midx, wincon->one, window->sgsizS * sizeof(double));

    /* Set the work arrays for cells where the face is partially dry */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = partial[cc];         /* 3D sparse of partially dry cell    */
      ch = hi_ele[cc];         /* 3D location of higher water column */
      cl = lo_ele[cc];         /* 3D location of lower water column  */
      cs = window->m2d[cl];    /* 2D loaction of lower water column  */
      xm1 = window->xm1[c];    /* Cell to the west of c              */

      bot = max(window->gridz[c], window->botzu1[window->m2d[c]]);
      hlower[c] = windat->eta[cs] - bot;
      hupper[c] = windat->dzu1[c] - hlower[c];
      hidens[c] = dens[ch];
      dav[c] = 0.5 * (dens[c] + dens[xm1]);
      mask[c] = 1.0;
      dzsum[c] = (at_ele[cc] == cl) ? -dens[xm1] : dens[c];
      /* Special treatment for hlower when the cell is completely    */
      /* dry (i.e. when eta <= botzu1).                              */
      if (windat->eta[cs] <= window->botzu1[window->m2d[c]]) { 
	hlower[c] = windat->dzu1[c];
	mask[c] = 0.0;
      }
    }

    /* Set the sub-surface explicit map if required                  */
    if (window->sm_e1) {
      for (cc = 1; cc <= wincon->vcs; cc++) {
	cs = window->m2d[wincon->s1[cc]];
	c = window->sm_e1[cs];
	if (c) {
	  hlower[c] = windat->dzu1[c];
	  hupper[c] = windat->eta[window->m2d[hi_ele[cc]]] - 
	    windat->eta[window->m2d[lo_ele[cc]]];
	  while (window->gridz[cs] > windat->eta[cs]) {
	    cs = window->zm1[cs];
	  }
	  dzsum[c] = (at_ele[cc] == lo_ele[cc]) ? -dens[window->xm1[cs]] : 
	    dens[cs];
	}
      }
    }

    /* Set the work arrays for cells where the face is wet. Where    */
    /* water exists on both sides of the cell face, set the density  */
    /* at the cell face to the mean of that at cell c and at cell    */
    /* xm1.                                                          */
    for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      xm1 = window->xm1[c];
      dav[c] = 0.5 * (dens[c] + dens[xm1]);
    }

    /*---------------------------------------------------------------*/
    /* SIGMA : Add the surface elevation gradient and apply the      */
    /* gradient correction. Note that sigma (gridz) takes on values  */
    /* from 0 to -1.                                                 */
    memset(corr, 0, window->sgsiz * sizeof(double));
    memset(dd, 0, window->sgsizS * sizeof(double));
    if (wincon->sigma) {
      int xmzm1;

      /* Subtract the mean density                                   */
      for (c = 1; c <= window->enon; c++)
        dens[c] -= wincon->rdens[c];

      for (cc = 1; cc <= wincon->vcs; cc++) {
        c = wincon->s1[cc];
        zm1 = window->zm1[c];
        cs = window->m2d[c];
        xm1s = window->xm1[cs];

        /* Surface pressure gradient                                 */
        ddbt[cs] = wincon->g * (windat->eta[cs] - windat->eta[xm1s]) +
          (windat->patm[cs] - windat->patm[xm1s]);

	/* Update velocity                                           */
	windat->nu1[c] +=
	  midx[cs] * windat->dt * wincon->u1c6[cs] * dd[cs];
        c = zm1;
        zm1 = window->zm1[c];

        while (c != zm1) {
          xm1 = window->xm1[c];
          zp1 = window->zp1[c];
          xmzp1 = window->zp1[xm1];
          xmzm1 = window->zm1[xm1];

          /* Density gradient correction                             */
          corr[c] = -wincon->g * midx[cs] * window->cellz[c] *
            (wincon->Ds[cs] - wincon->Ds[xm1s]) *
            0.5 * (dens[zp1] + dens[xmzp1] - dens[zm1] - dens[xmzm1]);

          /* Update velocity                                         */
          windat->nu1[c] +=
            midx[cs] * windat->dt * wincon->u1c6[cs] * dd[cs];

          c = zm1;
          zm1 = window->zm1[c];
        }
      }
      /* Get the barotropic pressure gradient tendency if required   */
      if (windat->u1_btp) {
        memset(windat->u1_btp, 0, window->sgsiz * sizeof(double));
        for (cc = 1; cc <= wincon->vc; cc++) {
          c = wincon->s1[cc];
          cs = window->m2d[c];
          windat->u1_btp[c] =
            midx[cs] * windat->dt * wincon->u1c6[cs] * dd[cs];
        }
      }
    }

    /*---------------------------------------------------------------*/
    /* Initialise vertical integral arrays                           */
    memset(dd, 0, window->sgsizS * sizeof(double));
    memset(u1inter, 0, window->sgsizS * sizeof(double));
    memset(dhdiv, 0, window->sgsizS * sizeof(double));
    /* memset(dhdiu,0,window->sgsizS*sizeof(double)); */

    /*---------------------------------------------------------------*/
    /* Reset the layer thickness over thin layers. The velocity      */
    /* change due to the pressure gradient in the thin layer is used */
    /* to update the velocity in the layer below the thin layer.     */
    /* This pressure gradient may not be the same as the layer below */
    /* since density may be different.                               */
    /*
    if (wincon->thin_merge) { 
      for (cc = 1; cc <= wincon->nkth_e1; cc++) {
	c = wincon->kth_e1[cc]; 
	zm1 = window->zm1[c];
	windat->dzu1[zm1] -= windat->dzu1[c]; 
	hupper[zm1] -= windat->dzu1[c]; 
      }
    }
    */

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the surface cells (completely or   */
    /* partially dry faces). The completely dry faces are the        */
    /* 'surface' part, and corresponds to the deta/de1 term in the   */
    /* 2D mode. For the sigma case this loop is always skipped.      */
    for (cc = 1; cc <= wincon->vca; cc++) {
      c = wincon->s1[cc];
      xm1 = window->xm1[c];
      cs = window->m2d[c];

      /* Get the pressure / surface gradients                        */
      d1 = midx[cs] * midx[cs] * (dens[c] - dens[xm1]) * mask[c];
      ddbc[cs] += wincon->g * d1 * hlower[c];
      dinter[cs] += wincon->g * mask[c] * d1 * hlower[c];
      ddbt[cs] += wincon->g * dzsum[c] * hupper[c];

      /* Update the u1 velocity                                      */
      windat->nu1[c] += windat->dt * wincon->u1c6[cs] * dd[cs] / dav[c];

      /* Integrate surface and average density for 2D mode           */
      wincon->topdensu1[cs] += hidens[c] * hupper[c];
      wincon->densavu1[cs] += dav[c] * hlower[c];
      wincon->u1inter[cs] +=
        wincon->u1c6[cs] * u1inter[cs] * hlower[c] / dav[c];
      dhdiv[cs] += hlower[c];
      /* dhdiu[cs] += hupper[c]; */
    }
    /* mark1 */

    /*---------------------------------------------------------------*/
    /* Get the barotropic pressure gradient tendency if required.    */
    if (windat->u1_btp && !wincon->sigma) {
      memset(wincon->d7, 0, window->sgsizS * sizeof(double));
      memset(btp, 0, window->sgsiz * sizeof(double));
      /* Surface                                                     */
      for (cc = 1; cc <= wincon->vca; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	wincon->d7[cs] += wincon->g * dzsum[c] * hupper[c];
	btp[c] = windat->dt * wincon->u1c6[cs] * wincon->d7[cs] / dav[c];
      }
      /* Sub-surface                                                 */
      for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	btp[c] = windat->dt * wincon->u1c6[cs] * wincon->d7[cs] / dav[c];
      }
      if (wincon->means & TENDENCY) {
	double t = windat->dtf;
	for (cc = 1; cc <= wincon->vc; cc++) {
	  c = wincon->s1[cc];
	  cs = window->m2d[c];
	  windat->u1_btp[c] = (windat->u1_btp[c] * windat->meanc[cs] + 
			       btp[c] * t) / (windat->meanc[cs] + t);
	}
      } else
	memcpy(windat->u1_btp, btp, window->sgsiz * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Add the thin layer density gradient contribution to velocity  */
    /* and set the thin layer thickness.                             */
    /*
    if (wincon->thin_merge) { 
      for (cc = 1; cc <= wincon->nkth_e1; cc++) {
       c = wincon->kth_e1[cc]; 
       xm1 = window->xm1[c]; 
       cs = window->m2d[c];
       xm1s = window->xm1[cs]; 
       zm1 = window->zm1[c];

       d1 = (windat->eta[cs]<windat->eta[xm1s]) ? -dens[xm1] : dens[c];
       dd[cs] += wincon->g*(d1*hupper[c]); 
       windat->nu1[zm1] += windat->dt*wincon->u1c6[cs]*dd[cs]/dav[c]; 
       wincon->topdensu1[cs] += hidens[c]*hupper[c]; 
       windat->dzu1[zm1]+=windat->dzu1[c]; 
      } 
    } 
    */

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the wet cell faces                 */
    for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {

      c = wincon->s1[cc];
      xm1 = window->xm1[c];
      cs = window->m2d[c];

      /* Get the pressure / surface gradients                        */
      d1 =
        wincon->g * midx[cs] * midx[cs] * (dens[c] - dens[xm1]) * dzu1[c];
      ddbc[cs] += (d1 + corr[c]);
      dinter[cs] += (d1 + corr[c]);

      /* Update the u1 velocity                                      */
      windat->nu1[c] += windat->dt * wincon->u1c6[cs] * dd[cs] / dav[c];

      /* Integrate surface and average density for 2D mode           */
      wincon->u1inter[cs] +=
        wincon->u1c6[cs] * u1inter[cs] * dzu1[c] / dav[c];
      wincon->densavu1[cs] += dav[c] * dzu1[c];
      dhdiv[cs] += dzu1[c];
    }

    /*---------------------------------------------------------------*/
    /* Divide mean densities by depth                                */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      xm1 = window->xm1[cs];
      wincon->densavu1[cs] /= dhdiv[cs];
      d1 = fabs(windat->topz[cs] - windat->topz[xm1]);
      /* d1=fabs(windat->eta[cs]-windat->eta[xm1]); */
      wincon->topdensu1[cs] =
        (d1 > tol) ? wincon->topdensu1[cs] / d1 : dav[c];
    }

    /* Reset the density by adding the mean density                  */
    if (wincon->sigma) {
      for (c = 1; c <= window->enon; c++)
        dens[c] += wincon->rdens[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the baroclinic pressure gradient tendency if required       */
  if (wincon->tendency && windat->u1_bcp) {
    int mf = wincon->means;        /* Save the means flag            */
    if (wincon->means & TENDENCY)  /* Turn mean TENDENCY off         */
      wincon->means &= ~TENDENCY;
    /* Get the baroclinic + barotropic tendency and store in w6      */
    get_tend(window, wincon->s1, wincon->vc, windat->nu1,
             wincon->tendency, wincon->w6);
    wincon->means = mf;            /* Reset the mean flag            */
    if (wincon->means & TENDENCY) {
      d1 = windat->dtf;
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	windat->u1_bcp[c] = (windat->u1_bcp[c] * windat->meanc[cs] + 
			     (wincon->w6[c] - btp[c]) * d1) / 
	  (windat->meanc[cs] + d1);
      }
    } else {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	windat->u1_bcp[c] = wincon->w6[c] - btp[c];
      }
    }
  }
}

/* END pressure_u1()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays for the u1 calculation                             */
/*-------------------------------------------------------------------*/
void precalc_u1(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int c, cc, cs;
  int xm1;                         /* Sparse coordinate at i-1       */
  int yp1;                         /* Sparse coordinate at j+1       */
  int xmyp1;                       /* Sparse coordinate at (i-1,j+1) */
  double *depth;
  double *u2au1;

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  windat->nu1 = wincon->w9;
  depth = wincon->d5;
  u2au1 = wincon->w7;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->u1inter, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Copy the current velocities into the update array               */
  memcpy(windat->nu1, windat->u1b, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations.                                             */
  if (wincon->sigma) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      windat->nu1[c] *= mdxbs(window, windat, wincon, c);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the depth at the e1 face.                                   */
  /* u1inter  is divided  by the depth  in vel_u1_update()  (for the */
  /* pressure and surface stress terms) but the dispersion terms are */
  /* not  to  be divided  by depth, so  multiply  by  depth here  to */
  /* compensate. Note: the work array wincon->d5 contains the height */
  /* of the surface at the e1 face at this stage, set in             */
  /* set_dz_at_u1().                                                 */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    /* Get the water depth at the e1 face                            */
    depth[cs] -= window->botzu1[cs];

    wincon->u1adv[cs] *= max(depth[cs], wincon->hmin);
  }

  /*-----------------------------------------------------------------*/
  /* Get the u2 velocity at the e1 face. Note : yp1 map must be able */
  /* to map to the u2 boundary using set_map_e2. Maps must be re-set */
  /* so as to be self mapping  over this boundary after this routine */
  /* is called.                                                      */
  set_map_e2(window);
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    u2au1[c] = 0.25 * (windat->u2[xm1] + windat->u2[xmyp1] +
                       windat->u2[c] + windat->u2[yp1]);
  }

  /* Initialise the tendency array of required                       */
  if (wincon->tendf && wincon->tendency)
    memcpy(wincon->tendency, windat->nu1, window->sgsiz * sizeof(double));
}

/* END precalc_u1()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the u1 velocities on all open boundaries in a  */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_u1_3d(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int n;                              /* Counters                    */
  open_bdrys_t **open = window->open; /* Open boundary structure     */

  /*-----------------------------------------------------------------*/
  /* Set the u1 velocity where u1 is normal to the boundary          */
  for (n = 0; n < window->nobc; n++) {

    set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
            open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
            open[n]->cyc_e1, windat->nu1, windat->u1, 
	    windat->u1b, open[n]->bcond_nor,
            windat->dtf, &open[n]->datau1, 
	    open[n]->transfer_u1, open[n]->relax_zone_nor, U1GEN);

    /* Set open boundary velocities above the free surface equal to  */
    /* zero.                                                         */
    set_dry_bdry(window, open[n]->no3_e1, open[n]->obc_e1, windat->nu1);
  }

  /*-----------------------------------------------------------------*/
  /* Set the u1 velocity where u1 is tangential to the boundary      */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
            open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
            open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
            windat->u1, windat->u1b, open[n]->bcond_tan, 
	    windat->dtf, &open[n]->datau1, open[n]->transfer_u1, 
	    open[n]->relax_zone_tan, U1GEN);            
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      int c, cc, c1;
      int *mi = (open[n]->ocodex & L_EDGE) ? window->xm1 : window->xp1;
      for (cc = 1; cc <= open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	c1 = open[n]->obc_t[cc];
	windat->nu1[c1] = windat->nu1[c];
      }
    }
  }

  /* Rescale for sigma active OBCs */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->type & U1BDRY && open->bcond_nor & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1, mdxns,
			1, open->no3_e1, open->obc_e1);
      if (open->type & U2BDRY && open->bcond_tan & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1, mdxns,
			open->no3_e1 + 1, open->to3_e1, open->obc_e1);
    }
  }
}

/* END bdry_u1_3d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the u1 velocity from the u1 velocity transport */
/*-------------------------------------------------------------------*/
void extract_u1_3d(geometry_t *window, /* Window geometry            */
                   window_t *windat,   /* Window data                */
                   win_priv_t *wincon  /* Window constants           */
  )
{
  int c, cc, cs, n;             /* Local sparse coordinate / counter */

  /* Extract the u1 velocity from velocity transport.                */
  /* Note : the water depths at the forward time steps have been     */
  /* calculated in the last 2D iteration in velocity_update2D()      */
  /* and can be used here.                                           */
  if (wincon->sigma) {
    for (cc = 1; cc <= window->b3_e1; cc++) {
      c = window->w3_e1[cc];
      cs = window->m2d[c];
      windat->nu1[c] /= wincon->Hn1[cs];
    }
  }
}

/* END extract_u1_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the u2 velocity from the u2 velocity transport */
/*-------------------------------------------------------------------*/
void extract_u2_3d(geometry_t *window,  /* Window geometry           */
                   window_t *windat,    /* Window data               */
                   win_priv_t *wincon   /* Window constants          */
  )
{
  int c, cc, cs, n;             /* Local sparse coordinate / counter */

  if (wincon->sigma) {
    for (cc = 1; cc <= window->b3_e2; cc++) {
      c = window->w3_e2[cc];
      cs = window->m2d[c];
      windat->nu2[c] /= wincon->Hn2[cs];
    }
  }
}

/* END extract_u2_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* LEAPFROG : Reset the velocities for the new time level, apply the */
/* Asselin time filter.                                              */
/*-------------------------------------------------------------------*/
void leapfrog_update_3d(geometry_t *window, /* Window geometry       */
                        window_t *windat,   /* Window data           */
                        win_priv_t *wincon  /* Window constants      */
  )
{
  int c, cc, n;                 /* Local sparse coordinate / counter */
  double aconst = 0.1;
  double *u1 = wincon->w1;
  double *u2 = wincon->w2;

  memcpy(u1, windat->u1, window->sgsiz * sizeof(double));
  memcpy(u2, windat->u2, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Apply the Asselin filter to remove the computational mode.      */
  /* u1 velocity.                                                    */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    windat->u1[c] = windat->u1[c] + 0.5 * aconst *
      (windat->u1b[c] - 2.0 * windat->u1[c] + windat->nu1[c]);
  }

  /* u2 velocity.                                                    */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    windat->u2[c] = windat->u2[c] + 0.5 * aconst *
      (windat->u2b[c] - 2.0 * windat->u2[c] + windat->nu2[c]);
  }

  /* If the open boundary condition for 2D velocity is FILEIN, then   */
  /* do not perform time filtering (restore boundary velocities to    */
  /* those before filtering is performed).                            */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->bcond_nor & (CUSTOM|FILEIN)) {
      if (open->type == U1BDRY) {
	for (cc = 1; cc <= open->no3_e1; cc++) {
	  c = open->obc_e1[cc];
	  windat->u1[c] = u1[c];
	}
      }
      if (open->type == U2BDRY) {
	for (cc = 1; cc <= open->no3_e2; cc++) {
	  c = open->obc_e2[cc];
	  windat->u2[c] = u2[c];
	}
      }
    }
    if (open->bcond_tan & (CUSTOM|FILEIN)) {
      if (open->type == U1BDRY) {
	for (cc = open->no3_e1 + 1; cc <= open->to3_e1; cc++) {
	  c = open->obc_e1[cc];
	  windat->u1[c] = u1[c];
	}
      }
      if (open->type == U2BDRY) {
	for (cc = open->no3_e2 + 1; cc <= open->to3_e2; cc++) {
	  c = open->obc_e2[cc];
	  windat->u2[c] = u2[c];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the velocities for the new time level                       */
  memcpy(windat->u1b, windat->u1, window->sgsiz * sizeof(double));
  memcpy(windat->u1, windat->nu1, window->sgsiz * sizeof(double));
  memcpy(windat->u2b, windat->u2, window->sgsiz * sizeof(double));
  memcpy(windat->u2, windat->nu2, window->sgsiz * sizeof(double));
}

/* END leapfrog_update_3d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset open boundary velocities                         */
/*-------------------------------------------------------------------*/
void reset_sigma_OBC(geometry_t *window, 
		     window_t *windat,
		     win_priv_t *wincon,
		     double *vel,
		     double (*Hn) (geometry_t *, window_t *, win_priv_t *, int),
		     int sb, int eb, int *obc)
{
  int c, cs, cc;

  for (cc = sb; cc <= eb; cc++) {
    c = obc[cc];
    cs = window->m2d[c];
    vel[c] *= Hn(window, windat, wincon, cs);
  }
}

/* END reset_sigma_OBC()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to pre-condition vertical velocity for the advection      */
/* computation.                                                      */
/*-------------------------------------------------------------------*/
void set_wvel(geometry_t *window,           /* Window geometry       */
	      window_t *windat,             /* Window data           */
	      win_priv_t *wincon,           /* Window constants      */
	      double *wvel,                 /* Conditioned w         */
	      int mode                      /* Mode=0:e1, mode=1:e2  */
	      )
{
  int c, cc, cs, cb, zm1, zp1;
  int *w2, *vec, *bot, *sur, a2, n2;
  int *map, m1;
  double top;

  /* Set pointers                                                    */
  if (mode == 0) {
    w2 = window->w2_e1;
    vec = wincon->s1;
    sur = window->sur_e1;
    bot = window->bot_e1;
    a2 = window->a2_e1;
    n2 = window->n2_e1;
    map = window->xm1;
  } else {
    w2 = window->w2_e2;
    vec = wincon->s2;
    sur = window->sur_e2;
    bot = window->bot_e2;
    a2 = window->a2_e2;
    n2 = window->n2_e2;
    map = window->ym1;
  }

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity. Note that w is cell centered, hence the centered      */
  /* bottom coordinate must be used. Auxiliary cells are not defined */
  /* for cell centers on western and southern boundaries, so the _t  */
  /* vectors cannot be used here and the cell center must be found   */
  /* by mapping downwards from bot_e1 until self-mapping is located. */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= a2; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = w2[cc];            /* 2D coordinate                     */
    cb = bot[cc];               /* Bottom coordinate                 */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs];
      c = window->zm1[c];
    }

    /* Set the bottom boundary condition. Find the cell centered     */
    /* bottom coordinate and set w to wbot.                          */
    c = cb;                       /* e1 face centered bottom cell    */
    while (c != window->zm1[c])
      c = window->zm1[c];         /* Cell centered sediment location */
    cb = window->zp1[c];
    wvel[cb] = windat->wbot[cs];
    wvel[c] = wvel[cb];
  }
  /* Set wvel at the surface for ghost cells. If window partitions   */
  /* lie next to ghost cells then these are not captured in the loop */
  /* above but are required for face centered vertical velocity.     */
  /* Note that wvel[cb] always = 0 and is not explicitly set.        */
  for (cc = a2 + 1; cc <= n2; cc++) {
    cs = c = w2[cc];
    m1 = map[cs];
    top = windat->eta[cs];
    while (c != window->zm1[c] && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs];
      c = window->zm1[c];
    }
  }
  /* Get the face average velocity (mean of cell c and xm1 or ym1   */
  /* cell). This is staggered so that the average at layer c is     */
  /* stored in the layer below. This allows values at every         */
  /* vertical face to be represented, with the bottom face stored   */
  /* in the sediment cell and sea surface face in the surface cell. */
  /* Loop from the bottom upwards and set the staggered average.    */
  for (cc = 1; cc <= a2; cc++) {
    c = bot[cc];               /* Bottom coordinate                 */
    cs = sur[cc];              /* Surface coordinate                */
    m1 = map[c];
    zm1 = window->zm1[c];
    wvel[zm1] = 0.5 * (wvel[c] + wvel[m1]);   /* Bottom             */
    while (c != cs && c != window->zp1[c]) {  /* Mid water          */
      zm1 = c;
      c = window->zp1[c];
      m1 = map[c];
      wvel[zm1] = 0.5 * (wvel[c] + wvel[m1]);
    }
    cs = window->m2d[c];
    wvel[c] = windat->wtop[cs];
  }
  /* Set the surface staggered average                               */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = vec[cc];
    m1 = map[c];
    wvel[c] = 0.5 * (wvel[c] + wvel[m1]);
  }
}

/* END set_wvel()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_tr_3d(geometry_t *window,  /* Window geometry       */
                       window_t *windat,    /* Window data           */
                       double trem, /* Time remaining in sub-step loop */
                       double *vel,  /* Velocity array               */
                       double *as,   /* Multi-dt cells velocity time t
                                       values */
                       double *ae    /* Multi-dt cells velocity time t+1
                                     values */
  )
{
  int cc, c;                    /* Sparse counter / coordinate       */

  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    if (c < window->enonS)
      vel[c] =
        (windat->dt - trem) * (ae[cc] - as[cc]) / window->taux_t[cc] +
        as[cc];
  }
}

/* END set_multidt_tr_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Update the 3D e2 velocity                                         */
/*-------------------------------------------------------------------*/
void vel_u2_update(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
  )
{
  int c, cc, cs;                /* Sparse coodinate / counter        */
  double *depth;                /* Depth of the water column         */
  

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  depth = wincon->d6;

  /*-----------------------------------------------------------------*/
  /* Get the cells to process at the e2 face                         */
  cells2process_e2(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation          */
  precalc_u2(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across u2 OBC's                 */
  set_map_e2(window);

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection                                     */
  if(wincon->momsc & (ANGULAR|ANGULAR3D)) {
    if(!(windat->nstep % 2))
      advect_u2_3d_ang_flux_b(window, windat, wincon);
    else
      advect_u2_3d_ang_flux_b(window, windat, wincon);
  }
  else {
    if (advect_u2_3d(window, windat, wincon)) return;
  }
  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s2, wincon->vc, windat->nu2,
             wincon->tendency, windat->u2_adv);
  debug_c(window, D_V, D_ADVECT);

  /*-----------------------------------------------------------------*/
  /* Add the dispersion terms to u2inter.                            */
  if (wincon->nonlinear && !(wincon->u2_f & ADVECT)) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      /* Add the dispersion terms to u2inter */
      wincon->u2inter[cs] =
        (wincon->u2inter[cs] + wincon->u2adv[cs]) / windat->dtf;
    }
  }
  /* Initialise the dispersion terms                                 */
  memset(wincon->u2adv, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included)                 */
  wincon->hor_mix->u2(window, windat, wincon);
  blend_vel(window, windat, wincon, U2_VH, wincon->u2vh);

  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s2, wincon->vc, windat->nu2,
             wincon->tendency, windat->u2_hdif);
  debug_c(window, D_V, D_HDIFF);

  /* Reset thin layers                                               */
  reset_thin_u2(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Include the pressure gradient term                              */
  pressure_u2(window, windat, wincon);
  debug_c(window, D_V, D_PRESSURE);

  /*-----------------------------------------------------------------*/
  /* Include the Coriolis term                                       */
  if (!(wincon->u2_f & CORIOLIS))
    coriolis_u2(window, windat, wincon);
  if (wincon->tendf && wincon->tendency) {
    get_tend(window, wincon->s2, wincon->vc, windat->nu2,
             wincon->tendency, windat->u2_cor);
    if (wincon->waves & STOKES_DRIFT) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s2[cc];
	windat->u2_cor[c] -= windat->u2_sto[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the vertical advection                                       */
  if (vdiff_u2(window, windat, wincon)) return;
  if (wincon->tendf && wincon->tendency)
    get_tend(window, wincon->s2, wincon->vc, windat->nu2,
             wincon->tendency, windat->u2_vdif);
  debug_c(window, D_V, D_VZ);

  /*-----------------------------------------------------------------*/
  /* Set up variables for the 2D mode                                */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];

    /* Divide internal terms for 2d part by depth                    */
    wincon->u2inter[cs] /= max(depth[cs], wincon->hmin);

    /* Other calculations to save time                               */
    wincon->topdensu2[cs] *= wincon->g;
    wincon->densavu2[cs] *= window->h2au2[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Blend velocities over a subregion if required                   */
  blend_vel(window, windat, wincon, U23D, windat->nu2);

  /*-----------------------------------------------------------------*/
  /* Calculate u2 boundary values                                    */
  bdry_u2_3d(window, windat, wincon);
  debug_c(window, D_V, D_POST);
}

/* END vel_u2_update()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness at the e1 faces.                */
/*-------------------------------------------------------------------*/
void set_dz_at_u2(geometry_t *window, /* Window geometry             */
                  window_t *windat,   /* Window data                 */
                  win_priv_t *wincon  /* Window constants            */
  )
{
  int cc;                       /* Counter                           */
  int c, c3;                    /* Sparse coordinate                 */
  int ym1;                      /* Sparse cell at (j-1)              */
  int zm1;                      /* Sparse cell below cell c3         */
  int cs, cb;                   /* Sparse coordinate of the bottom   */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  double *maxeta;               /* Maximum elevation at the e2 face  */

  /*-----------------------------------------------------------------*/
  /* Set the pointers                                                */
  maxeta = wincon->d6;

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->x2_e2 + 1;
    memcpy(window->sur_e2, window->w2_e2, size * sizeof(int));

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom.                                                       */
    memset(windat->dzu2, 0, window->sgsiz * sizeof(double));
    memset(maxeta, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->b2_e2; cc++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      c = c3 = window->sur_e2[cc];  /* 3D surface coordinate         */
      cs = window->m2d[c];          /* 2D cell corresponding to c    */
      ym1 = window->ym1[cs];        /* Cell to the south of cs       */
      cb = window->bot_e2[cc];      /* 3D bottom coordinate          */

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = windat->topz[cs];
      while (c != cb) {
        bot = window->gridz[c];
        windat->dzu2[c] = top - bot;
        top = bot;
        c = window->zm1[c];
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu2[cb] = top - window->botzu2[cs];

    }
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Find the 3D sparse coordinate corresponding to the free       */
    /* surface.                                                      */
    memset(windat->sur_e2, 0, window->sgsizS * sizeof(int));
    for (cc = 1; cc <= window->b2_e2; cc++) {

      /* Get the 2D and 3D sparse coordinates                        */
      cs = c = window->w2_e2[cc]; /* 2D coordinate                   */
      ym1 = window->ym1[cs];      /* Cell to the south of cs         */

      /* Get the new sparse location of the surface. The surface at  */
      /* the e2 face is defined as the higher of the cell centered   */
      /* elevations bordering the face. For elevations below the     */
      /* bottom this is set to the sediment coordinate.              */
      top = maxeta[cs] = max(windat->eta[cs], windat->eta[ym1]);
      top = maxeta[cs] = max(top, window->botzu2[cs]);
      zm1 = window->zm1[c];

      /* Loop down the water column                                  */
      while (c != zm1 && window->gridz[c] > top) {
        c = zm1;
        zm1 = window->zm1[c];
      }
      window->sur_e2[cc] = c;
      windat->sur_e2[cs] = window->s2k[c];
    }

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom.                                                       */
    memset(windat->dzu2, 0, window->sgsiz * sizeof(double));
    wincon->ncdry_e2 = 0;
    for (cc = 1; cc <= window->b2_e2; cc++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      c = c3 = window->sur_e2[cc];  /* 3D surface coordinate         */
      cs = window->m2d[c];          /* 2D cell corresponding to c    */
      ym1 = window->ym1[cs];        /* Cell to the south of cs       */
      cb = window->bot_e2[cc];      /* 3D bottom coordinate          */

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = maxeta[cs];
      while (c != cb) {
        bot = window->gridz[c];
        windat->dzu2[c] = top - bot;
        top = bot;
        c = window->zm1[c];
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu2[cb] = top - window->botzu2[cs];

      /* Save the locations of dry water columns                     */
      if (windat->dzu2[cb] == 0.0) {
	wincon->ncdry_e2++;
	wincon->cdry_e2[wincon->ncdry_e2] = cs;
      }

      /* Set the cell thickness above the old surface equal to zero. */
      /* (Note : a no-gradient condition may be required for higher  */
      /* order vertical advection schemes.                           */
      c = c3;
      top = windat->dzu2[c];
      while (c != cs) {
        c = window->zp1[c];
        windat->dzu2[c] = top;
        /*windat->dzu2[c] = 0.0;*/
      }
    }
  }
}

/* END set_dz_at_u2()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the velocity over thin layers                      */
/*-------------------------------------------------------------------*/
void set_thin_u2(geometry_t *window,    /* Window geometry           */
		 window_t *windat,      /* Window data               */
		 win_priv_t *wincon     /* Window constants          */
  )
{
  int cc, bn;                   /* Counter                           */
  int c, cs;                    /* Sparse coordinate                 */
  int cb;                       /* Bottom sparse coordinate          */
  int zm1;                      /* Sparse cell below cell c          */
  int *cells;                   /* Cells to process                  */
  int *bottom;                  /* Bottom cells                      */
  double *u2o = wincon->d7;     /* Original u2 at zm1                */
  open_bdrys_t **open = window->open; /* Window OBC structure        */

  /*-----------------------------------------------------------------*/
  /* Set the thin layer vector and reset the surface layer.          */
  /* Note : if the loop is done to a2_e2 then dzu2=0 at the first    */
  /* timestep for multiple windows in auxiliary cells since it has   */
  /* not yet been copied from the master.                            */
  cells = wincon->s3;
  bottom = wincon->i6;
  wincon->nkth_e2 = 1;

  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s2[cc];
    cb = bottom[cc];
    cs = window->m2d[c];
    zm1 = window->zm1[c];
    if (c < cb && windat->dzu2[zm1] && windat->dzu2[c] < wincon->hmin) {
      cells[cc] = zm1;
      /* Save the sparse location of the thin layer                  */
      wincon->kth_e2[wincon->nkth_e2] = c;
      wincon->nkth_e2++;
      u2o[window->m2d[c]] = windat->u2[zm1];
      merge_thin_layers(c, zm1, windat->u2, windat->dzu2);
    } else
      cells[cc] = c;
    windat->sur_e2[cs] = 0;
  }

  /* Set thin layers for open boundary cells. These are used in the  */
  /* non-linear terms.                                               */
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	zm1 = window->zm1[c];
	if (windat->dzu2[zm1] && windat->dzu2[c] < wincon->hmin)
	  merge_thin_layers(c, zm1, windat->u2, windat->dzu2);
	windat->sur_e2[cs] = 0;
      }
    }
  }

  /* Set thin layers for auxiliary cells. These are used in the      */
  /* non-linear terms.                                               */
  if (window->nwindows > 1) {
    for (cs = 1; cs <= window->enonS; cs++) {
      int ks = windat->sur_e2[cs];
      if (ks) {
	c = get_local_sur(window, cs, ks);
	zm1 = window->zm1[c];
	if (windat->dzu2[zm1] && windat->dzu2[c] < wincon->hmin) {
	  wincon->kth_e2[wincon->nkth_e2] = c;
	  wincon->nkth_e2++;
	  u2o[cs] = windat->u2[zm1];
	  merge_thin_layers(c, zm1, windat->u2, windat->dzu2);
	}
      }
    }
  }
  wincon->nkth_e2--;

  /* Set a no-gradient of the merged velocity above the surface      */
  for (cc = 1; cc <= wincon->nkth_e2; cc++) {
    c = wincon->kth_e2[cc];
    zm1 = window->zm1[c];
    for(cs = window->m2d[c]; cs != c; cs = window->zm1[cs]) {
      windat->u2[cs] = windat->u2[zm1];
    }
  }

  /* Fill the cells to process with sub-surface cells. If no thin    */
  /* cells were encountered make a direct copy of s1 to save time.   */
  if (wincon->nkth_e2 == 0) {
    memcpy(wincon->s3, wincon->s2, window->sgsiz * sizeof(int));
    wincon->ncl = wincon->vc;
  } else {
    wincon->ncl = wincon->vcs + 1;
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->zm1[cells[cc]];
      cb = bottom[cc];
      while (c <= cb) {
        cells[wincon->ncl] = c;
        wincon->ncl++;
        c = window->zm1[c];
      }
    }
    wincon->ncl--;
  }
}

/* END set_thin_u2()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset u2 velocity in thin layers                       */
/*-------------------------------------------------------------------*/
void reset_thin_u2(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
  )
{
  int cc;                       /* Counter                           */
  int c,zm1;                    /* Sparse coordinate                 */
  double *u2o = wincon->d7;     /* Original u2 at zm1                */

  /* Set a uniform velocity over thin layers if required             */
  if (wincon->thin_merge) {
    for (cc = 1; cc <= wincon->nkth_e2; cc++) {
      c = wincon->kth_e2[cc];
      zm1 = window->zm1[c];
      windat->nu2[c] = windat->nu2[zm1];
      windat->dzu2[zm1] -= windat->dzu2[c];
      /* Reset u2 to the un-merged velocity                          */
      windat->u2[c] = (windat->dzu2[c]) ? windat->u2[zm1] + 
	               windat->dzu2[zm1] * 
 	              (windat->u2[zm1] - u2o[window->m2d[c]]) / 
	               windat->dzu2[c] : 0.0;
      windat->u2[zm1] = u2o[window->m2d[c]];
    }
  }
}

/* END reset_thin_u2()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to perform thin layer merging                             */
/*-------------------------------------------------------------------*/
void merge_thin_layers(int c, int zm1,  /* Sparse coordinates        */
                       double *vel,     /* Velocity                  */
                       double *dz       /* Layer thickness           */
  )
{
  /* Set the velocity equal to the mean of the two layers            */
  vel[zm1] = (dz[c] * vel[c] + dz[zm1] * vel[zm1]) / (dz[zm1] + dz[c]);
  /* Merge the two layers                                            */
  vel[c] = vel[zm1];
  dz[zm1] += dz[c];
}

/* END merge_thin_layers()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the sparse coordinates of cells to process at the  */
/* e2 faces.                                                         */
/*-------------------------------------------------------------------*/
void cells2process_e2(geometry_t *window, /* Window geometry         */
                      window_t *windat,   /* Window data             */
                      win_priv_t *wincon  /* Window constants        */
  )
{
  int cc;                       /* Counter                           */
  int c, c2;                    /* Sparse coordinate                 */
  int ym1;                      /* Sparse cell at (j-1)              */
  int cs, cb;                   /* Sparse coordinate of the bottom   */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  double *maxeta;               /* Maximum elevation at the e2 face  */
  int vc, vcs;                  /* Counter                           */

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->x2_e2 + 1;
    wincon->vc = wincon->vc2 = window->v3_e2;
    wincon->vcs = wincon->vcs2 = window->v2_e2;
    wincon->vca = wincon->vca2 = 0;
    memcpy(wincon->s2, window->w3_e2, (window->a3_e2 + 1) * sizeof(int));
    memcpy(wincon->i1, window->sur_e2, size * sizeof(int));
    memcpy(wincon->i2, window->sur_e2, size * sizeof(int));
    memcpy(wincon->i3, window->sur_e2, size * sizeof(int));
    memcpy(wincon->i4, window->sur_e2, size * sizeof(int));
    memcpy(wincon->i6, window->bot_e2, size * sizeof(int));
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Set the pointers                                              */
    int *partial = wincon->i1;
    int *at_ele = wincon->i2;
    int *hi_ele = wincon->i3;
    int *lo_ele = wincon->i4;
    int *bottom = wincon->i6;
    maxeta = wincon->d6;
    vcs = (wincon->dobdry_u2) ? window->b2_e2 : window->v2_e2;

    /* Set the sub-surface explicit map if required                  */
    if (window->sm_e2) {
      for (cc = 1; cc <= window->b2_e2; cc++) {
	cs = window->w2_e2[cc];
	c = window->sm_e2[cs];
	if (c) {
	  window->sur_e2[cc] = c;
	  maxeta[cs] = min(maxeta[cs], window->gridz[window->zp1[c]]);
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the surface wet cells to process vector for u2 velocity   */
    /* and store in the buffer wincon->s2. Do not include cells      */
    /* where the free surface is below the bottom (i.e. correspond   */
    /* to sediment sparse locations).                                */
    vc = 1;
    wincon->ncbot_e2 = 0;
    for (cc = 1; cc <= vcs; cc++) {
      c = window->sur_e2[cc];   /* Surface coordinate                */
      cs = window->w2_e2[cc];   /* 2D coordinate                     */
      ym1 = window->ym1[cs];    /* Cell to the south of cs           */
      top = maxeta[cs];         /* Highest elevation at the face     */
      bot = DRY_FRAC * wincon->hmin + window->botzu2[cs];
      if (top > bot) {
        wincon->s2[vc] = c;
        cb = bottom[vc] = window->bot_e2[cc];
	/* Save cells which are one cell deep                        */
	if(c == cb) {
	  wincon->ncbot_e2++;
	  wincon->cbot_e2[wincon->ncbot_e2] = c;
	}
	/*
	lo_ele[vc] = (top > windat->eta[ym1]) ? window->ym1[c] : c;
	hi_ele[cc] = (c == lo_ele[cc]) ? window->ym1[c] : c;
        lo_ele[vc] = c;
        if (top > windat->eta[ym1])
          lo_ele[vc] = window->ym1[c];
	*/
        vc++;
      }
    }
    wincon->vcs = vc - 1;

    /*---------------------------------------------------------------*/
    /* Loop from the layer below the surface to the bottom and get   */
    /* the cells to process. The are arranged so that cells 1 to     */
    /* wincon->vcs contain surface cells, cells wincon->vcs+1 to     */
    /* wincon->vca contain cells where one side of the cell face is  */
    /* dry and cells wincon->vca+1 to wincon->vc contain cells where */
    /* both sides of the cell face is wet.                           */
    /* First get the cells where one side of the face is dry. The    */
    /* last cell that falls into this category is temporarily saved  */
    /* in partial[] for use in the next loop.                        */
    for (cc = 1; cc <= wincon->vcs; cc++) {

      c = c2 = wincon->s2[cc];  /* 3D surface cell                   */
      cs = window->m2d[c];      /* 2D cell corresponding to c        */
      cb = bottom[cc];          /* 3D bottom cell                    */
      ym1 = window->ym1[cs];    /* Cell to the south of cs           */
      at_ele[cc] = c;           /* Cell containing sea level         */

      top = min(windat->eta[cs], windat->eta[ym1]);

      while (c < cb && (bot = max(window->gridz[c], 
				  window->botzu2[cs])) > top) {
        /* Cells where one face is completely dry                    */

        if (c != c2) {
          wincon->s2[vc] = c;
          vc++;
        }
        if (bot > windat->eta[cs])
          at_ele[cc] = window->zm1[at_ele[cc]];
        c = window->zm1[c];
      }

      if (c <= cb) {
        /* Cells where one face is partially dry                     */
        if (c != c2) {
          wincon->s2[vc] = c;
          vc++;
        }
        partial[cc] = hi_ele[cc] = c;
        lo_ele[cc] = window->ym1[c];
        if (top < windat->eta[ym1]) {
          lo_ele[cc] = c;
          hi_ele[cc] = window->ym1[c];
        }
      }
    }
    wincon->vca = vc - 1;
    /* Next get the cells where both sides of the face is wet.       */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = window->zm1[partial[cc]];
      cb = bottom[cc];
      while (c <= cb) {
        wincon->s2[vc] = c;
        vc++;
        c = window->zm1[c];
      }
    }
    wincon->vc = vc - 1;
    /* Remember the number of u2 cells to process for the 2D mode    */
    wincon->vc2 = wincon->vc;
    wincon->vcs2 = wincon->vcs;
    wincon->vca2 = wincon->vca;
  }
}

/* END cells2process_e2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Vertical diffusion                                                */
/*-------------------------------------------------------------------*/
int vdiff_u2(geometry_t *window, /* Window geometry */
	     window_t *windat, /* Window data           */
	     win_priv_t *wincon  /* Window constants */
  )
{
  int c, cc;                    /* Sparse coordinates / counters */
  int cs, cb;                   /* Surface / bottom sparse coordinate */
  int ym1;                      /* Cell location to the west of c */
  int *at_ele;                  /* Array of fully wet cell coordinates */
  int *partial;                 /* Coordinate of cell containing eta[c] */
  int *bottom;                  /* Bottom coordinate for cells to process */
  double *Vzav;                 /* Viscosity at the e1 face */
  double *f_top;                /* Flux out of the water at the surface */
  double *f_bot;                /* Flux into the water at the bottom */
  double *u1au2;                /* u2 velocity at e1 face */
  double Cdu2;                  /* Bottom drag at the e1 face */
  double botdz;                 /* Thickness of the bottom layer */
  double val;                   /* Bottom current speed */
  double *w;                    /* Vertical velocity                 */
  int xp1;                      /* Sparse coordinate at i+1 */
  int xpym1;                    /* Sparse coordinate at (i+1,j-1) */
  double ramp = (wincon->rampf & WIND) ? windat->rampval : 1.0;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* s2 = wet cells to process for u2 velocity */
  partial = wincon->i1;
  bottom = wincon->i6;
  at_ele = wincon->i2;
  f_top = wincon->d1;
  f_bot = wincon->d2;
  Vzav = wincon->w1;
  u1au2 = wincon->w7;
  w = wincon->w8;

  /*-----------------------------------------------------------------*/
  /* Set the vertical viscosity at the e1 face. Where one side of */
  /* the water column is higher, set the viscosity at the cell face */
  /* equal to the viscosity of the higher water column. If the */
  /* current water column, cs, is the higher then no action is taken */
  /* (c == cs), if the water column to the west of face cs is higher */
  /* then the viscosity at c is set to the viscosity at ym1. This is */
  /* only performed to the layer above that containing the free */
  /* surface in column cs; below this the mean is taken.  */
  memcpy(Vzav, windat->Vz, window->sgsiz * sizeof(double));
  /* Set the viscosity on completely dry faces */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s2[cc];
    cs = at_ele[cc];
    while (c < cs) {
      ym1 = window->ym1[c];
      Vzav[c] = windat->Vz[ym1];
      c = window->zm1[c];
    }
  }
  /* Set the viscosity for cells where the face is partially dry */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = partial[cc];
    ym1 = window->ym1[c];
    Vzav[c] = 0.5 * (windat->Vz[c] + windat->Vz[ym1]);
  }
  /* Set the viscosity for cells where the face completely wet */
  for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
    c = wincon->s2[cc];
    ym1 = window->ym1[c];
    Vzav[c] = 0.5 * (windat->Vz[c] + windat->Vz[ym1]);
  }

  /*-----------------------------------------------------------------*/
  /* Set the vertical velocity                                       */
  memset(w, 0, window->sgsiz * sizeof(double));
  if (wincon->momsc & WIMPLICIT) {
    int zm1;
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];         /* 3D cell to process              */
      ym1 = window->ym1[c];
      zm1 = window->zm1[c];
      w[zm1] = 0.5 * (windat->w[c] + windat->w[ym1]);
    }
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];       /* 3D cell to process                */
      cs = window->m2d[c];      /* 2D cell corresponding to c        */
      cb = bottom[cc];          /* 3D bottom coordinate              */
      zm1 = window->zm1[cb];
      ym1 = window->ym1[cs];
      /* Bottom vertical velocity                                    */
      w[zm1] = 0.5 * (windat->wbot[ym1] + windat->wbot[cs]);
      /* Surface vertical velocity                                   */
      w[c] = 0.5 * (windat->wtop[ym1] + windat->wtop[cs]);
    }
    set_wvel(window, windat, wincon, w, 1);
  }

  /*-----------------------------------------------------------------*/
  /* Set the momentum flux out of the water at the top (due to the */ 
  /* wind) and flux into the water at the bottom (due to friction).  */
  memset(f_top, 0, window->sgsizS * sizeof(double));
  memset(f_bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s2[cc];         /* 3D cell to process */
    cs = window->m2d[c];        /* 2D cell corresponding to c */
    cb = bottom[cc];            /* 3D bottom coordinate */
    ym1 = window->ym1[cb];
    xp1 = window->xp1[cb];
    xpym1 = window->xpym1[cb];

    /* Surface flux */
    f_top[cs] = -ramp * windat->wind2[cs] / wincon->topdensu2[cs];

    /* Add surface stress to 2d forcing term.  */
    wincon->u2inter[cs] -= f_top[cs];

    /* Bottom flux */
    u1au2[cb] = 0.25 * (windat->u1b[ym1] + windat->u1b[xpym1] +
                        windat->u1b[cb] + windat->u1b[xp1]);
    Cdu2 = 0.5 * (wincon->Cd[cs] + wincon->Cd[window->ym1[cs]]);
    val = sqrt(windat->u2b[cb] * windat->u2b[cb] + u1au2[cb] * u1au2[cb]);
    val = Cdu2 * max(wincon->uf, val);
    botdz = max(wincon->hmin, windat->dzu2[cb] * wincon->Ds[cs]);
    /* Quadratic bottom friction, truncated to ensure stability */
    if (val > botdz / windat->dt)
      val = botdz / windat->dt;
    f_bot[cs] = -val * windat->u2b[cb];
  }

  if (wincon->vorticity & TENDENCY)
    memcpy(windat->rv_wsc, f_bot, window->sgsizS * sizeof(double));
  if (wincon->numbers & BOTSTRESS)
    memcpy(windat->tau_be2, f_bot, window->sgsizS * sizeof(double));
  if (wincon->numbers & EKPUMP)
    ekman_pump_e2(window, windat, wincon, windat->wind2, f_bot);

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion.  */
  /* SIGMA : multiply by depth to give D.u1 at t+1. Note: should be */
  /* dividing by depth at t+1 but this is not possible since */
  /* elevation is only solved after 3D velocity.  */
  if (!(wincon->u2_f & VDIFF)) {
    if (implicit_vdiff(window, windat, wincon, windat->u2b, windat->nu2, Vzav,
		       windat->dzu2, f_bot, f_top, wincon->s2, wincon->i6,
		       wincon->vcs, wincon->Ds, wincon->Ds, w)) return(1);

  }
  return(0);
}

/* END vdiff_u2()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to add the Coriolis term to the u2 velocity               */
/*-------------------------------------------------------------------*/
void coriolis_u2(geometry_t *window,  /* Window geometry */
                 window_t *windat,  /* Window data           */
                 win_priv_t *wincon /* Window constants */
  )
{
  int c, cc, cs;                /* Sparse coordinates / counters */
  double *u1au2;                /* u1 velocity at e2 face */
  double *midy;                 /* Total depth at e2 face */

  /* Set pointers */
  midy = wincon->d4;            /* Set in pressure_u2() */
  u1au2 = wincon->w7;           /* Set in precalc_u2() */

  /*-----------------------------------------------------------------*/
  /* Add the Coriolis term to the u1 velocity */
  /* SIGMA : multiply by total depth at u1 */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->nu2[c] += windat->dt * wincon->u2c5[cs] * u1au2[c] * midy[cs];
  }
  if(wincon->waves & STOKES_DRIFT) {
    double sdc;
    double *tend;
    /* Stokes drift                                                  */
    tend = (wincon->tendf) ? windat->u2_sto : wincon->w2;
    memset(tend, 0.0, window->sgsiz * sizeof(double));
    get_sdc_e1(window, windat, wincon);
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      sdc = wincon->u2c5[cs] * wincon->w1[c] * midy[cs];
      windat->nu2[c] += windat->dt * sdc;
      wincon->u2inter[cs] += sdc * windat->dzu2[c];
      tend[c] += windat->dt * sdc;
    }
    /* Vortex force                                                  */
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      sdc = wincon->w1[c] * midy[cs] * 
	((windat->u2[window->xp1[c]] - windat->u2[window->xm1[c]]) / (2.0 * window->h1au2[cs]) -
	 (u1au2[window->yp1[c]] - u1au2[window->ym1[c]]) / (2.0 * window->h2au2[cs]));
      windat->nu2[c] += windat->dt * sdc;
      wincon->u2inter[cs] += sdc * windat->dzu2[c];
      tend[c] += windat->dt * sdc;
    }
  }
}

/* END coriolis_u2()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U2 horizontal pressure gradient term.                             */
/*                                                                   */
/* If eta is treated in linear fashion, then the term looks like     */
/*                                               _                   */
/*                  dP             deta         |  dro               */
/*                  --   =  g*ro0* ----  +   g* |  --- dz            */
/*                  de2             de2         |  de2               */
/*                                             -                     */
/* where ro0 is the surface density                                  */
/*             ((dens[nz-1][j][cl]+dens[nz-1][j][cr])/2)             */
/* and the integral is from z=0.0 down to the bottom of the layer    */
/* currently being calculated.                                       */
/*                                                                   */
/* If eta is nonlinear, the water surface may be in different        */
/* layers on either side of the u1 column. Here, the pressure        */
/* gradient term is:                                                 */
/*                                _                                  */
/*                  dP           |  dro                              */
/*                  --   =    g* |  --- dz                           */
/*                  de2          |  de2                              */
/*                              -                                    */
/* where the integral is from the higher water surface               */
/* (max(topz[j][cl],topz[j][cr])) to the bottom of the layer         */
/* currently being calculated.                                       */
/*-------------------------------------------------------------------*/
void pressure_u2(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc;              /* Local sparse coordinate / counter       */
  int cs;                 /* 2D local sparse coordinate              */
  int cl;                 /* Partially wet partial cell coordinate   */
  int ch;                 /* Fully wet partial cell coordinate       */
  int ym1;                /* Cell location to the west of c          */
  int ym1s;               /* Surface cell location to the west of c  */
  int zm1;                /* Cell below cell c                       */
  int zp1, ymzp1;         /* Cell above and back-above cell c        */
  int *partial;           /* Coordinate of cell containing eta[c]    */
  int *at_ele;            /* Array of fully wet cell coordinates     */
  int *hi_ele;            /* Array of fully wet cell coordinates     */
  int *lo_ele;            /* Array of partially wet cell coordinates */
  double *dhdiv;          /* Vertically integrated depth             */
  double *dd;             /* Vertically integrated pressure gradient */
  double *dinter;         /* Pressure gradient for the 2D mode       */
  double *u2inter;        /* Vertical integrals for 2d mode          */
  double *hidens;         /* Density corresponding to coordinate ch  */
  double *dav;            /* Mean density at the cell face           */
  double *dzu2;           /* Cell theickness ( = windat->dzu2)       */
  double *hupper;         /* Partial cell upper thickness work array */
  double *hlower;         /* Partial cell lower thickness work array */
  double *dens;           /* Density ( = windat->dens)               */
  double *dzsum;          /* Sum of dz for dry / partial layers      */
  double *mask;           /* Mask out all fully wet or dry faces     */
  double *btp;            /* Barotropic pressure tendency            */
  double top;             /* Top cell thickness                      */
  double bot;             /* Bottom cell thickness                   */
  double *midy;           /* Total depth at the cell face            */
  double d1;              /* Dummy                                   */
  double tol = 1e-10;     /* Tolerence value                         */
  double *corr;           /* Baroclinic correction for sigma         */
  double *ddbc;           /* Baroclinic contribution to u2 velocity  */
  double *ddbt;           /* Barotropic contribution to u2 velocity  */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  hidens = wincon->w1;
  hupper = wincon->w2;
  hlower = wincon->w3;
  dav = wincon->w4;
  dzsum = wincon->w5;
  mask = wincon->w6;
  corr = wincon->w8;
  btp = wincon->w3;

  dd = wincon->d1;
  dhdiv = wincon->d2;
  u2inter = wincon->d3;
  midy = wincon->d4;

  dens = windat->dens;
  dzu2 = windat->dzu2;

  partial = wincon->i1;
  at_ele = wincon->i2;
  hi_ele = wincon->i3;
  lo_ele = wincon->i4;

  if (wincon->u2_f & PRESS_BC)
    ddbc = wincon->d7;
  else
    ddbc = dd;
  if (wincon->u2_f & PRESS_BT)
    ddbt = wincon->d7;
  else
    ddbt = dd;
  if (wincon->u2av_f & PRESS_BC)
    dinter = wincon->d7;
  else
    dinter = u2inter;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->topdensu2, 0, window->sgsizS * sizeof(double));
  memset(wincon->densavu2, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Linear eta case; calculate surface pressure gradient due to eta */
  /* slope and atmospheric pressure gradient.                        */
  if (!wincon->nonlinear) {
    double *surgrad = wincon->d1;

    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];
      ym1 = window->ym1[c];
      cs = window->m2d[c];
      ym1s = window->ym1[cs];

      midy[cs] = wincon->mdy[cs];
      wincon->topdensu2[cs] = 0.5 * (windat->dens[ym1] + windat->dens[c]);
      surgrad[cs] = wincon->g * wincon->topdensu2[cs] *
        (windat->eta[cs] - windat->eta[ym1s]) +
        (windat->patm[cs] - windat->patm[ym1s]);
      dd[cs] = 0.0;
    }
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      ym1 = window->ym1[c];
      cs = window->m2d[c];
      ym1s = window->ym1[cs];

      /* Density gradient integral                                   */
      dav[c] = 0.5 * (windat->dens[c] + windat->dens[ym1]);
      dd[cs] +=
        wincon->g * (windat->dens[c] - windat->dens[ym1]) * dzu2[c];

      /* Add term to new u2 value                                    */
      windat->nu2[c] +=
        windat->dt * wincon->u2c6[cs] * (surgrad[cs] + dd[cs]) / dav[c];

      /* Integrate internal term for 2d part                         */
      wincon->u2inter[cs] += wincon->u2c6[cs] * dd[cs] * dzu2[c] / dav[c];
      wincon->densavu2[cs] += dav[c] * dzu2[c];
    }
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      ym1s = window->ym1[cs];
      top = max(windat->eta[ym1s], windat->eta[cs]);
      bot = window->botzu2[cs];
      wincon->densavu2[cs] /= (top - bot);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Non-linear eta case                                             */
  else {

    /*---------------------------------------------------------------*/
    /* Precondition work arrays for the pressure gradient integral.  */
    /* These are established for the surface cells (i.e. cells that  */
    /* have one face completely or partially dry) as follows:        */
    /*                                                               */
    /* work array   | completely dry face  | partially dry face      */
    /* -------------|----------------------|----------------------   */
    /* hidens       | rho at wet cell      | rho at full wet cell    */
    /* hupper       | cell thickness, dzu2 | upper thickness, hu     */
    /* hlower       | cell thickness, dzu2 | lower thickness, hl     */
    /* dav          | rho at wet cell      | 0.5(rho(c)+rho(ym1))    */
    /* dzsum        | rho(c) or -rho(ym1)  | rho(c) or -rho(ym1)     */
    /* mask         | zero                 | one                     */

    /* Initialise the work arrays                                    */
    memcpy(hidens, dens, window->sgsiz * sizeof(double));
    memcpy(hupper, dzu2, window->sgsiz * sizeof(double));
    memcpy(hlower, dzu2, window->sgsiz * sizeof(double));
    memcpy(dav, dens, window->sgsiz * sizeof(double));
    /*memset(dzsum, 0, window->sgsiz * sizeof(double));*/
    memcpy(dzsum, dens, window->sgsiz * sizeof(double));
    memset(mask, 0, window->sgsiz * sizeof(double));

    /* Set the work arrays for cells where the face is completely    */
    /* dry. Where one side of the water column is higher, set the    */
    /* density at the cell face equal to the density of the higher   */
    /* water column. If the current water column, cs, is the higher  */
    /* then no action is taken (c == cs), if the water column to the */
    /* south of face cs is higher then the density at c is set to    */
    /* the density at ym1. This is only performed to the layer above */
    /* that containing the free surface in column cs.                */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];
      cs = at_ele[cc];
      midy[window->m2d[cs]] =
        0.5 * (wincon->Ds[window->m2d[cs]] + 
	       wincon->Ds[window->ym1[window->m2d[cs]]]);
      while (c < cs) {
        ym1 = window->ym1[c];
        hidens[c] = dav[c] = dens[ym1];
	dzsum[c] = -dens[ym1];
        c = window->zm1[c];
      }
    }
    if (!wincon->sigma)
      memcpy(midy, wincon->one, window->sgsizS * sizeof(double));

    /* Set the work arrays for cells where the face is partially dry */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = partial[cc];
      ch = hi_ele[cc];
      cl = lo_ele[cc];
      cs = window->m2d[cl];
      ym1 = window->ym1[c];
      bot = max(window->gridz[c], window->botzu2[window->m2d[c]]);
      hlower[c] = max(windat->eta[cs] - bot, 0.0);
      hupper[c] = windat->dzu2[c] - hlower[c];
      hidens[c] = dens[ch];
      mask[c] = 1.0;
      dav[c] = 0.5 * (dens[c] + dens[ym1]);
      dzsum[c] = (at_ele[cc] == cl) ? -dens[ym1] : dens[c];
      /* Special treatment for hlower when the cell is completely    */
      /* dry (i.e. when eta <= botzu2).                              */
      if(windat->eta[cs] <= window->botzu2[window->m2d[c]]) {
	hlower[c] = windat->dzu2[c];
	mask[c] = 0.0;
      }
    }

    /* Set the sub-surface explicit map if required                  */
    if (window->sm_e2) {
      for (cc = 1; cc <= wincon->vcs; cc++) {
	cs = window->m2d[wincon->s2[cc]];
	c = window->sm_e2[cs];
	if (c) {
	  hlower[c] = windat->dzu2[c];
	  hupper[c] = windat->eta[window->m2d[hi_ele[cc]]] - 
	    windat->eta[window->m2d[lo_ele[cc]]];
	  while (window->gridz[cs] > windat->eta[cs]) {
	    cs = window->zm1[cs];
	  }
	  dzsum[c] = (at_ele[cc] == lo_ele[cc]) ? -dens[window->ym1[cs]] : 
	    dens[cs];
	}
      }
    }

    /* Set the work arrays for cells where the face is wet. Where   */
    /* water exists on both sides of the cell face, set the density */
    /* at the cell face to the mean of that at cell c and at cell   */
    /* ym1.                                                         */
    for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      ym1 = window->ym1[c];
      dav[c] = 0.5 * (dens[c] + dens[ym1]);
    }

    /*---------------------------------------------------------------*/
    /* SIGMA : Add the surface elevation gradient and apply the      */
    /* gradient correction. Note that sigma (gridz) takes on values  */
    /* from 0 to -1.                                                 */
    memset(corr, 0, window->sgsiz * sizeof(double));
    memset(dd, 0, window->sgsizS * sizeof(double));

    if (wincon->sigma) {
      int ymzm1;

      /* Subtract the mean density                                   */
      for (c = 1; c <= window->enon; c++)
        dens[c] -= wincon->rdens[c];

      for (cc = 1; cc <= wincon->vcs; cc++) {
        c = ch = wincon->s2[cc];
        zm1 = window->zm1[c];
        cs = window->m2d[c];
        ym1s = window->ym1[cs];

        /* Surface pressure gradient                                 */
        ddbt[cs] = wincon->g * (windat->eta[cs] - windat->eta[ym1s]) +
          (windat->patm[cs] - windat->patm[ym1s]);

	/* Update velocity                                           */
	windat->nu2[c] +=
	  midy[cs] * windat->dt * wincon->u2c6[cs] * dd[cs];
        c = zm1;
        zm1 = window->zm1[c];

        while (c != zm1) {
          ym1 = window->ym1[c];
          zp1 = window->zp1[c];
          ymzp1 = window->zp1[ym1];
          ymzm1 = window->zm1[ym1];

          /* Density gradient correction                             */
          corr[c] = -wincon->g * midy[cs] * window->cellz[c] *
            (wincon->Ds[cs] - wincon->Ds[ym1s]) *
            0.5 * (dens[zp1] + dens[ymzp1] - dens[zm1] - dens[ymzm1]);

          /* Update velocity                                         */
          windat->nu2[c] +=
            midy[cs] * windat->dt * wincon->u2c6[cs] * dd[cs];

          c = zm1;
          zm1 = window->zm1[c];
        }
      }
      /* Get the barotropic pressure gradient tendency if required   */
      if (windat->u2_btp) {
        memset(windat->u2_btp, 0, window->sgsiz * sizeof(double));
        for (cc = 1; cc <= wincon->vc; cc++) {
          c = wincon->s2[cc];
          cs = window->m2d[c];
          windat->u2_btp[c] =
            midy[cs] * windat->dt * wincon->u2c6[cs] * dd[cs];
        }
      }
    }

    /*---------------------------------------------------------------*/
    /* Initialise vertical integral arrays                           */
    memset(dd, 0, window->sgsizS * sizeof(double));
    memset(u2inter, 0, window->sgsizS * sizeof(double));
    memset(dhdiv, 0, window->sgsizS * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the surface cells (completely or   */
    /* partially dry faces). The completely dry faces are the        */
    /* 'surface' part, and corresponds to the deta/de2 term in the   */
    /* 2D mode. For the sigma case this loop is always skipped.      */
    for (cc = 1; cc <= wincon->vca; cc++) {
      c = wincon->s2[cc];
      ym1 = window->ym1[c];
      cs = window->m2d[c];

      /* Get the pressure / surface gradients                        */
      d1 = midy[cs] * midy[cs] * (dens[c] - dens[ym1]) * mask[c];
      /* dd[cs] += wincon->g*(d1*dzu2[c]+dzsum[c]*hupper[c]); */
      ddbc[cs] += wincon->g * d1 * hlower[c];
      dinter[cs] += wincon->g * mask[c] * d1 * hlower[c];
      ddbt[cs] += wincon->g * dzsum[c] * hupper[c];

      /* Update the u2 velocity                                      */
      windat->nu2[c] += windat->dt * wincon->u2c6[cs] * dd[cs] / dav[c];

      /* Integrate surface and average density for 2D mode           */
      wincon->topdensu2[cs] += hidens[c] * hupper[c];
      wincon->densavu2[cs] += dav[c] * hlower[c];
      wincon->u2inter[cs] +=
        wincon->u2c6[cs] * u2inter[cs] * hlower[c] / dav[c];
      dhdiv[cs] += hlower[c];
    }
    /* mark2 */

    /*---------------------------------------------------------------*/
    /* Get the barotropic pressure gradient tendency if required.    */
    /* Note: the barotropic pressure gradient is depth independent.  */
    if (windat->u2_btp && !wincon->sigma) {
      memset(wincon->d7, 0, window->sgsizS * sizeof(double));
      memset(btp, 0, window->sgsiz * sizeof(double));
      /* Surface                                                     */
      for (cc = 1; cc <= wincon->vca; cc++) {
	c = wincon->s2[cc];
	cs = window->m2d[c];
	wincon->d7[cs] += wincon->g * dzsum[c] * hupper[c];
	btp[c] = windat->dt * wincon->u2c6[cs] * wincon->d7[cs] / dav[c];
      }
      /* Sub-surface                                                 */
      for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
	c = wincon->s2[cc];
	cs = window->m2d[c];
	btp[c] = windat->dt * wincon->u2c6[cs] * wincon->d7[cs] / dav[c];
      }
      if (wincon->means & TENDENCY) {
	double t = windat->dtf;
	for (cc = 1; cc <= wincon->vc; cc++) {
	  c = wincon->s2[cc];
	  cs = window->m2d[c];
	  windat->u2_btp[c] = (windat->u2_btp[c] * windat->meanc[cs] + 
			       btp[c] * t) / (windat->meanc[cs] + t);
	}
      } else
	memcpy(windat->u2_btp, btp, window->sgsiz * sizeof(double));
    }

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the wet cell faces                 */
    for (cc = wincon->vca + 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      ym1 = window->ym1[c];
      cs = window->m2d[c];

      /* Get the pressure / surface gradients                        */
      d1 =
        wincon->g * midy[cs] * midy[cs] * (dens[c] - dens[ym1]) * dzu2[c];
      ddbc[cs] += (d1 + corr[c]);
      dinter[cs] += (d1 + corr[c]);

      /* Update the u2 velocity */
      windat->nu2[c] += windat->dt * wincon->u2c6[cs] * dd[cs] / dav[c];

      /* Integrate surface and average density for 2D mode           */
      wincon->u2inter[cs] +=
        wincon->u2c6[cs] * u2inter[cs] * dzu2[c] / dav[c];
      wincon->densavu2[cs] += dav[c] * dzu2[c];
      dhdiv[cs] += dzu2[c];
    }

    /*---------------------------------------------------------------*/
    /* Divide means densities by depth                               */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      ym1 = window->ym1[cs];
      wincon->densavu2[cs] /= dhdiv[cs];
      d1 = fabs(windat->topz[cs] - windat->topz[ym1]);
      /* d1=fabs(windat->eta[cs]-windat->eta[ym1]); */
      wincon->topdensu2[cs] =
        (d1 > tol) ? wincon->topdensu2[cs] / d1 : dav[c];
    }

    /* Reset the density by adding the mean density                  */
    if (wincon->sigma) {
      for (c = 1; c <= window->enon; c++)
        dens[c] += wincon->rdens[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the baroclinic pressure gradient tendency if required       */
  if (wincon->tendency && windat->u1_bcp) {
    int mf = wincon->means;        /* Save the means flag            */
    if (wincon->means & TENDENCY)  /* Turn mean TENDENCY off         */
      wincon->means &= ~TENDENCY;
    /* Get the baroclinic + barotropic tendency and store in w6      */
    get_tend(window, wincon->s2, wincon->vc, windat->nu2,
	     wincon->tendency, wincon->w6);
    wincon->means = mf;            /* Reset the mean flag            */
    if (wincon->means & TENDENCY) {
      d1 = windat->dtf;
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s2[cc];
	cs = window->m2d[c];
	windat->u2_bcp[c] = (windat->u2_bcp[c] * windat->meanc[cs] + 
			     (wincon->w6[c] - btp[c]) * d1) / 
	  (windat->meanc[cs] + d1);
      } 
    } else {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s2[cc];
	windat->u2_bcp[c] = wincon->w6[c] - btp[c];
      }
    }
  }
}

/* END pressure_u2()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up arrays for the u2 calculation                             */
/*-------------------------------------------------------------------*/
void precalc_u2(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int c, cc, cs;
  int xp1;                      /* Sparse coordinate at i+1          */
  int ym1;                      /* Sparse coordinate at j-1          */
  int xpym1;                    /* Sparse coordinate at (i+1,j-1)    */
  double *depth;
  double *u1au2;

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  windat->nu2 = wincon->w10;
  depth = wincon->d6;
  u1au2 = wincon->w7;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->u2inter, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Copy the current velocities into the update array               */
  memcpy(windat->nu2, windat->u2b, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations.                                             */
  if (wincon->sigma) {
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s2[cc];
      windat->nu2[c] *= mdybs(window, windat, wincon, c);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the depth at the e1 face.                                   */
  /* Note: the work array wincon->d6 contains the height of the      */
  /* surface at the e2 face at this stage, set in set_dz_at_u2().    */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    /* Get the water depth at the e2 face                            */
    depth[cs] -= window->botzu2[cs];
    wincon->u2adv[cs] *= max(depth[cs], wincon->hmin);
  }

  /*-----------------------------------------------------------------*/
  /* Get the u1 velocity at the e2 face                              */
  set_map_e1(window);
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s2[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    u1au2[c] = 0.25 * (windat->u1[ym1] + windat->u1[xpym1] +
                       windat->u1[c] + windat->u1[xp1]);
  }

  /* Initialise the tendency array of required                       */
  if (wincon->tendf && wincon->tendency)
    memcpy(wincon->tendency, windat->nu2, window->sgsiz * sizeof(double));
}

/* END precalc_u2()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the u2 velocities on all open boundaries in a  */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_u2_3d(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int n;                            /* Counter                       */
  open_bdrys_t **open = window->open;

  /*-----------------------------------------------------------------*/
  /* Set the u2 velocity where u2 is normal to the boundary          */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e2,
            open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
            open[n]->cyc_e2, windat->nu2, windat->u2,
	    windat->u2b, open[n]->bcond_nor,
            windat->dtf, &open[n]->datau2, 
	    open[n]->transfer_u2, open[n]->relax_zone_nor, U2GEN);

    /* Set open boundary velocities above the free surface equal to  */
    /* zero.                                                         */
    set_dry_bdry(window, open[n]->no3_e2, open[n]->obc_e2, windat->nu2);
  }

  /*-----------------------------------------------------------------*/
  /* Set the u2 velocity where u2 is tangential to the boundary      */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
            open[n]->to3_e2, open[n]->obc_e2, open[n]->oi1_e2,
            open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2,
	    windat->u2, windat->u2b, open[n]->bcond_tan, windat->dtf,
            &open[n]->datau2, open[n]->transfer_u2, 
	    open[n]->relax_zone_tan, U2GEN);
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      int c, cc, c1;
      int *mi = (open[n]->ocodey & B_EDGE) ? window->ym1 : window->yp1;
      for (cc = 1; cc <= open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	c1 = open[n]->obc_t[cc];
	windat->nu2[c1] = windat->nu2[c];
      }
    }
  }

  /* Rescale for sigma active OBCs */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->type & U2BDRY && open->bcond_nor & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu2, mdyns,
			1, open->no3_e2, open->obc_e2);
      if (open->type & U1BDRY && open->bcond_tan & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu2, mdyns,
			open->no3_e2 + 1, open->to3_e2, open->obc_e2);
    }
  }
}

/* END bdry_u2_3d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set velocities above the free surface equal to zero on */
/* open boundaries.                                                  */
/*-------------------------------------------------------------------*/
void set_dry_bdry(geometry_t *window,             /* Window geometry */
	     int nb, int *obc, double *vel)
{
  int c, cc, cs;
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;

  if (wincon->sigma)
    return;

  for (cc = 1; cc <= nb; cc++) {
    c = obc[cc];
    cs = window->m2d[c];
    vel[c] = (windat->eta[cs] > window->gridz[c] * wincon->Ds[cs]) ? 
      vel[c] : 0.0;
  }
}

/* END set_dry_bdry()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the 3D fluxes through the cell faces. This is      */
/* performed over all wet cells only; auxiliary cell 3D fluxes are   */
/* set through a transfer from the master. This is preferrable to    */
/* calculating the 3D fluxes over all wet + auxiliary cells in the   */
/* window and dispensing with the need for transfers of 3D fluxes,   */
/* because complications can arise in the cells to process vectors   */
/* for bottom auxiliary cells where the bathymetry is stepped at the */
/* corresponding wet cell. In this case the wet cell may be a ghost  */
/* cell and thus not included in the cells to process vector v3_e1   */
/* or v3_e2, hence no corresponding auxiliary cell exists in a3_e1   */
/* or a3_e2 and 3D fluxes will not be set at the bottom, even though */
/* there may be a non-zero velocity associated with the bottom       */
/* auxiliary cell. This upsets the vertical velocity calculation at  */
/* the wet cell which subsequently results in incorrect surface      */
/* tracer calculations (continuity is violated).                     */
/*                                                                   */
/*              w  |  w  |  a  |    Side view.                       */
/*            -------------------   wet is a ghost cell for u1 &     */
/*                 |     |     |    aux is hence not defined.        */
/*             land| wet | aux |                                     */
/*                 |     |     |                                     */
/*                 --------------                                    */
/*                 |     |     |                                     */
/*                                                                   */
/*-------------------------------------------------------------------*/
void set_flux_3d(geometry_t *window,    /* Window geometry           */
                 window_t *windat,      /* Window data               */
                 win_priv_t *wincon,    /* Window constants          */
                 int mode               /* Velocity arrays to use    */
  )
{
  int c, cc, cs;                /* Sparse coodinate / counter        */
  double *u, *v;
  int n;

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  if (mode & TRANSPORT) {
    u = windat->u1m;
    v = windat->u2m;
  } else {
    u = windat->u1;
    v = windat->u2;
  }
  
  memset(windat->u1flux3d, 0, window->sgsiz * sizeof(double));
  memset(windat->u2flux3d, 0, window->sgsiz * sizeof(double));

  /* If using flux-form advection scheme in transport mode, transfer */
  /* u1vm to u1flux3d and u2vm to u2flux3d.                          */
  if (wincon->trasc == FFSL) {

    /* Set the velocities at auxillary cells to zero */
    for (cc = window->b3_e1 + 1; cc <= window->n3_e1; cc++) {
      c = window->w3_e1[cc];
      windat->u1[c] = 0.0;
      windat->u1vm[c] = 0.0;
    }
    for (cc = window->b3_e2 + 1; cc <= window->n3_e2; cc++) {
      c = window->w3_e2[cc];
      windat->u2[c] = 0.0;
      windat->u2vm[c] = 0.0;
    }

    for (cc = 1; cc <= window->b3_e1; cc++) {
      c = window->w3_e1[cc];
      windat->u1flux3d[c] = windat->u1vm[c];
    }
    for (cc = 1; cc <= window->b3_e2; cc++) {
      c = window->w3_e2[cc];
      windat->u2flux3d[c] = windat->u2vm[c];
    }

    return;
  }
  
  /*-----------------------------------------------------------------*/
  /* Fluxes through the e1 faces                                     */
    for (cc = 1; cc <= window->b3_e1; cc++) {
      c = window->w3_e1[cc];
      cs = window->m2d[c];
      windat->u1flux3d[c] = u[c] * windat->dzu1[c] * window->h2au1[cs] *
	wincon->mdx[cs];    
    }

    /* Set the fluxes for cells one cell deep                          */
    for (cc = 1; cc <= wincon->ncbot_e1; cc++) {
      c = wincon->cbot_e1[cc];
      cs = window->m2d[c];
      windat->u1flux3d[c] = windat->u1flux[cs]/windat->dt;
    }
    /* Zero fluxes above the surface                                   */
    for (cc = 1; cc <= window->b2_e1; cc++) {
      c = window->w2_e1[cc];
      cs = window->sur_e1[cc];
      while (c < cs) {
        windat->u1flux3d[c] = 0.0;
        c = window->zm1[c];
      }
    }

    /* Get the mean u1 velocity if required                            */
    if (wincon->means & VEL3D && windat->u1m) {
      double t = windat->dtf;
      for (cc = 1; cc <= window->b3_e1; cc++) {
	c = window->w3_e1[cc];
	cs = window->m2d[c];
	/*if (windat->u1flux3d[c])
	  Averaging when layers are wet can lead to time normanlization  */
	/* issues when free surface drops out of layers; use the         */
	/* no-gradient on u1 to get a surface mean in the top layer.     */
	  windat->u1m[c] = (windat->u1m[c] * windat->meanc[cs] + 
			    windat->u1[c] * t) / (windat->meanc[cs] + t);
      }
    }

    /* Get the mean u1 volume flux if required                         */
    if (wincon->means & VOLFLUX && windat->u1vm) {
      double t = windat->dtf;
      for (cc = 1; cc <= window->b3_e1; cc++) {
	c = window->w3_e1[cc];
	cs = window->m2d[c];
	windat->u1vm[c] = (windat->u1vm[c] * windat->meanc[cs] + 
			   windat->u1flux3d[c] * t) / (windat->meanc[cs] + t);
      }
    }

  /*-----------------------------------------------------------------*/
  /* Fluxes through the e2 faces                                     */
    for (cc = 1; cc <= window->b3_e2; cc++) {
      c = window->w3_e2[cc];
      cs = window->m2d[c];
      windat->u2flux3d[c] = v[c] * windat->dzu2[c] * window->h1au2[cs] *
        wincon->mdy[cs];
    }
    /* Set the fluxes for cells one cell deep                          */
    for (cc = 1; cc <= wincon->ncbot_e2; cc++) {
      c = wincon->cbot_e2[cc];
      cs = window->m2d[c];
      windat->u2flux3d[c] = windat->u2flux[cs]/windat->dt;
    }
    /* Zero fluxes above the surface                                   */
    for (cc = 1; cc <= window->b2_e2; cc++) {
      c = window->w2_e2[cc];
      cs = window->sur_e2[cc];
      while (c < cs) {
        windat->u2flux3d[c] = 0.0;
        c = window->zm1[c];
      }
    }
  
    /* Get the mean u2 velocity if required                            */
    if (wincon->means & VEL3D && windat->u2m) {
      double t = windat->dtf;
      for (cc = 1; cc <= window->b3_e2; cc++) {
	c = window->w3_e2[cc];
	cs = window->m2d[c];
	/*if (windat->u2flux3d[c])*/
	  windat->u2m[c] = (windat->u2m[c] * windat->meanc[cs] + 
			    windat->u2[c] * t) / (windat->meanc[cs] + t);
      }
    }

    /* Get the mean u2 volume flux if required                         */
    if (wincon->means & VOLFLUX && windat->u2vm) {
      double t = windat->dtf;
      for (cc = 1; cc <= window->b3_e2; cc++) {
	c = window->w3_e2[cc];
	cs = window->m2d[c];
	windat->u2vm[c] = (windat->u2vm[c] * windat->meanc[cs] + 
			   windat->u2flux3d[c] * t) / (windat->meanc[cs] + t);
      }
    }
}

/* END set_flux_3d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the 3D velocity so that the vertical integral   */
/* equals the 2D velocity. This routine assumes water depth and      */
/* grid spacings correspond to the start of the time-step, t.        */
/*-------------------------------------------------------------------*/
void velocity_adjust(geometry_t *window,    /* Window geometry       */
                     window_t *windat,      /* Window data           */
                     win_priv_t *wincon     /* Window constants      */
  )
{
  int c, cc, cs, cb;            /* Sparse coodinate / counter        */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *midx;                 /* SIGMA : depth at the e1 face      */
  double *midy;                 /* SIGMA : depth at the e1 face      */
  double *depth;                /* Water depth at the cell face      */
  double *adjust;               /* Velocity adjustment               */
  int bn;                       /* OBC counter                       */
  open_bdrys_t **open = window->open; /* Window OBC structure        */

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e1 faces.                                 */
  /* Set pointers and initialise                                     */
  sum = wincon->d1;
  midx = wincon->d2;
  depth = wincon->d5;           /* Set in precalc_u1()               */
  adjust = wincon->d3;

  memset(sum, 0, window->sgsizS * sizeof(double));
  /* The  number of surface  u1 cells to process (vcs1) may decrease */
  /* if a cell dries  during the 2D mode - hence adjust[] may not be */
  /* set  but u1 is still updated (with an old value residing in d3. */
  /* Therefore, zero adjust[] before use here.                       */
  memset(adjust, 0, window->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column. Note: normal */
  /* velocity open boundary sparse coordinates are not included in   */
  /* the cells to process vectors (whereas tangential boundary       */
  /* velocities are) and must be processed separately.               */
  /* Wet cells to process.                                           */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    midx[cs] = wincon->mdx[cs];
    sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
  }

  /* Open boundary normal velocities. Velocities above the surface   */
  /* are set to zero, so it is possible to integrate over the whole  */
  /* water column rather than just to the free surface.              */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	midx[cs] = wincon->mdx[cs];
	sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
      }
    }
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u1flux[cs] -
       sum[cs] * windat->dt * window->h2au1[cs]) / (depth[cs] *
                                                    windat->dt * midx[cs] *
                                                    window->h2au1[cs]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e1; cc++) {
	cs = open[bn]->obc_e1[cc];
	depth[cs] -= window->botzu1[cs];
	adjust[cs] =
	  (windat->u1flux[cs] -
	   sum[cs] * windat->dt * window->h2au1[cs]) / (depth[cs] *
							windat->dt *
							midx[cs] *
							window->h2au1[cs]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    windat->u1[c] += adjust[cs];
  }

  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	windat->u1[c] += adjust[cs];
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e1; cc++) {
    c = cs = wincon->cdry_e1[cc];
    windat->u1av[cs] = 0.0;
    windat->u1flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u1[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u1bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i5[cc];        /* 3D bottom coordinate              */
    windat->u1bot[cs] = windat->u1b[cb] - windat->u1avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e1; cc++) {
    c = wincon->cbot_e1[cc];
    cs = window->m2d[c];
    windat->u1[c] = windat->u1av[cs];
    windat->u1bot[cs] = 0.0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e2 faces.                                 */
  /* Set pointers and initialise                                     */
  midy = wincon->d2;
  depth = wincon->d6;           /* Set in precalc_u2()               */
  memset(sum, 0, window->sgsizS * sizeof(double));
  memset(adjust, 0, window->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column.              */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    midy[cs] = wincon->mdy[cs];
    sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	midy[cs] = wincon->mdy[cs];
	sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
      }
    }
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u2flux[cs] -
       sum[cs] * windat->dt * window->h1au2[cs]) / (depth[cs] *
                                                    windat->dt * midy[cs] *
                                                    window->h1au2[cs]);
  }

  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e2; cc++) {
	cs = open[bn]->obc_e2[cc];
	depth[cs] -= window->botzu2[cs];
	adjust[cs] =
	  (windat->u2flux[cs] -
	   sum[cs] * windat->dt * window->h1au2[cs]) / (depth[cs] *
							windat->dt *
							midy[cs] *
							window->h1au2[cs]);
      }

      /* Compute the mean integrated flux normal to the boundary     */
      /* Boundary adjusted volume constraint : Marchesiello et al    */
      /* (2001) Section 4.2.                                         */
      /*
      if (open[bn]->type & U2BDRY && open[bn]->adjust_flux) {
	double flux2d = 0.0;
	double flux3d = 0.0;
	double area = 0.0;
	for (cc = 1; cc <= open[bn]->no2_e2; cc++) {
	  cs = open[bn]->obc_e2[cc];
	  flux3d += sum[cs] * windat->dt * window->h1au2[cs];
	  flux2d += windat->u2flux[cs];
	  area += window->h1au2[cs] * depth[cs] * midy[cs];
	}
	open[bn]->bflux_2d = (open[bn]->bflux_2d * open[bn]->meanc +
			      flux2d * tm) / (open[bn]->meanc + tm);
	open[bn]->bflux_3d = (open[bn]->bflux_3d * open[bn]->meanc +
			      flux3d * tm) / (open[bn]->meanc + tm);
	if (windat->nstep) open[bn]->meanc += tm;
      }
      */
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->u2[c] += adjust[cs];
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	windat->u2[c] += adjust[cs];
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e2; cc++) {
    c = cs = wincon->cdry_e2[cc];
    windat->u2av[cs] = 0.0;
    windat->u2flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u2[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u2bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i6[cc];        /* 3D bottom coordinate              */
    windat->u2bot[cs] = windat->u2b[cb] - windat->u2avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e2; cc++) {
    c = wincon->cbot_e2[cc];
    cs = window->m2d[c];
    windat->u2[c] = windat->u2av[cs];
    windat->u2bot[cs] = 0.0;
  }
}

/* END velocity_adjust()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the u1 velocities at newly wetted cells.           */
/*-------------------------------------------------------------------*/
void set_new_cells_u1(geometry_t *window, /* Window geometry         */
                      window_t *windat,   /* Window data             */
                      win_priv_t *wincon  /* Window constants        */
  )
{
  int c, cc;                /* Sparse coodinate / counter            */
  int cn, co;               /* New and old surface sparse coodinates */
  int *oldsur;              /* Sparse coordinate of old surface      */

  /* Set pointers and initialise                                     */
  oldsur = wincon->i7;

  /*UR-FIX added/changed to prevent memory violation - read failure
   * orig
   * memcpy(oldsur, window->sur_e1, window->sgsizS * sizeof(int));
   * should be
   * memcpy(oldsur, window->sur_e1, (1+window->x2_e1) * sizeof(int));
   * this seems to be feasible, arrays only accessed at first level
   * window->x2_e1+1 - the length of 1d int window->sur_e1 is allocated
   */ 
  if(window->x2_e1 < window->b2_e1 ) {
    hd_warn("vel3d:set_new_cells_u1 - Assigned size is smaller than required copy %d < %d (grid %d )",
	    window->x2_e1,window->b2_e1,window->sgsizS);
  }
  memcpy(oldsur, window->sur_e1, (1+window->x2_e1) * sizeof(int));

  /* Get the location of the new surface. Note; by including thin    */
  /* layers in sur_e1 it appears as if elevation has risen if the    */
  /* layer was thin at the start of the time-step. The thin layer    */
  /* now will receive a copy of the (2D adjusted) velocity from the  */
  /* layer below.                                                    */
  set_dz_at_u1(window, windat, wincon);

  /* Set the velocities at newly wetted u1 cells                     */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];

    /* Set a no-gradient  velocity condition  above the new surface. */
    /* When sea  level rises  into new layers  it is possible that a */
    /* cell face at c is  wet but at xm1 (or ym1) is dry. This makes */
    /* the difference in u2h1h2 between  cell c and xm1  in the non- */
    /* linear advective  terms  large, which  can lead to velocities */
    /* that  are   unrealistically  large. A  no-gradient  condition */
    /* above  the  surface  minimises  this   effect  and  maintains */
    /* stability.                                                    */
    cn = window->sur_e1[cc];
    while (c < cn) {
      /*windat->u1[c]=0.0;*/
      windat->u1[c] = windat->u1[oldsur[cc]];
      c = window->zm1[c];
    }
    /* Set newly wetted cells equal to velocity at the old surface   */
    c = co = oldsur[cc];
    if (c > cn)
      c = window->zp1[c];
    while (c > cn) {
      windat->u1[c] = windat->u1[co];
      c = window->zp1[c];
    }
    windat->u1[c] = windat->u1[co];
  }

  /* Set the velocity at solid faces above explicit maps equal zero  */
  if (window->sm_e1) {
    int cb;
    for (cc = 1; cc <= window->b2_e1; cc++) {
      c = window->w2_e1[cc];
      cb = window->sm_e1[c];
      if (cb) {
	while (c < cb) {
	  windat->u1[c] = 0.0;
	  c = window->zm1[c];
	}
      }
    }
  }
}

/* END set_new_cells_u1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the u1 velocities at newly wetted cells.           */
/*-------------------------------------------------------------------*/
void set_new_cells_u2(geometry_t *window,   /* Window geometry       */
                      window_t *windat,     /* Window data           */
                      win_priv_t *wincon    /* Window constants      */
  )
{
  int c, cc;                /* Sparse coodinate / counter            */
  int cn, co;               /* New and old surface sparse coodinates */
  int *oldsur;              /* Sparse coordinate of old surface      */

  /* Set pointers and initialise                                     */
  oldsur = wincon->i7;
  
    /*UR added/changed to prevent memory violation - read failure
   * orig
   * memcpy(oldsur, window->sur_e2, window->sgsizS * sizeof(int));
   * should be 
   *  memcpy(oldsur, window->sur_e2, ( 1 + window->x2_e2) * sizeof(int));
   * this seems to be feasible, arrays only accessed at first level
   * window->x2_e2+1 - the length of 1d int window->sur_e2 is allocated
   */ 
  if(window->x2_e2 < window->b2_e2)
  {
  	hd_warn("vel3d:set_new_cells_u2 - Assigned size is smaller than required copy %d < %d (grid %d )",
  			window->x2_e2,window->b2_e2,window->sgsizS);
  }
  memcpy(oldsur, window->sur_e2, (1+window->x2_e2) * sizeof(int));

  /* Get the location of the new surface                             */
  set_dz_at_u2(window, windat, wincon);

  /* Set the velocities at newly wetted u2 cells                     */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    /* Set a no-gradient velocity condition above the new surface    */
    cn = window->sur_e2[cc];
    while (c < cn) {
      /*windat->u2[c]=0.0;*/
      windat->u2[c] = windat->u2[oldsur[cc]];
      c = window->zm1[c];
    }

    /* Set newly wetted cells equal to velocity at the old surface   */
    c = co = oldsur[cc];
    while (c > cn) {
      windat->u2[c] = windat->u2[co];
      c = window->zp1[c];
    }
    windat->u2[c] = windat->u2[co];
  }

  /* Set the velocity at solid faces above explicit maps equal zero  */
  if (window->sm_e2) {
    int cb;
    for (cc = 1; cc <= window->b2_e2; cc++) {
      c = window->w2_e2[cc];
      cb = window->sm_e2[c];
      if (cb) {
	while (c < cb) {
	  windat->u2[c] = 0.0;
	  c = window->zm1[c];
	}
      }
    }
  }
}

/* END set_new_cells_u2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the vertical velocity through continuity.              */
/*-------------------------------------------------------------------*/
void vel_w_update(geometry_t *window, /* Window geometry             */
                  window_t *windat,   /* Window data                 */
                  win_priv_t *wincon  /* Window constants            */
  )
{
  int c, cc, cs;      /* Sparse coodinate / counter                  */
  int *bottom;        /* Bottom sparse coordinate                    */
  int *sur;           /* Minimum of surface coordinate               */
  int *nsur;          /* Surface sparse coordinate after the 2D mode */
  double fctop;       /* Flux at the top face of the cell            */
  double fcbot;       /* Flux at the bottom face of the cell         */
  int xp1, yp1;       /* Sparse coordinates east and south of cell c */
  int zp1, zm1;       /* Sparse coordinate below cell c              */
  double hf, vf, cnt; /* Diagnostics for continuity                  */
  int dc = 0;         /* Diagnostic continuity surface coordinate    */
  int dw = 1;         /* Window for diagnostic continuity            */
  double detadt;      /* Rate of change of elevation                 */
  double *deta;       /* Elevation rate of change for SIGMA          */
  double d1;          /* Dummy                                       */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  memset(windat->w, 0, window->sgsiz * sizeof(double));
  bottom = wincon->i1;          /* Set in set_dz()                   */
  nsur = wincon->i2;            /* Set in set_dz()                   */
  sur = wincon->i3;             /* Set in set_dz()                   */
  deta = wincon->w1;

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across e1 and e2 OBC's          */
  reset_map_t(window);

  /*-----------------------------------------------------------------*/
  /* Calculate the vertical velocity at the surface and bottom       */
  /* boundaries.                                                     */
  if (wincon->momsc & WTOP_O4)
    vel_w_bounds_hiorder(window, windat, wincon);
  else
    vel_w_bounds(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Add the contribution from eta relaxation. Note that these       */
  /* fluxes are in/out of the surface cell and are not actually used */
  /* in the calculation of w, which loops from the bottom coordinate */
  /* to the cell below the surface. The surface condition on w       */
  /* implicitly includes these fluxes through the gradients of eta.  */
  /* Surface tracer mass fluxes due to eta relaxation must be        */
  /* explicitly added and this is handled in the routine ss_tracer() */
  /* in inputs/sourcesink.c. The above also applies for evaporation, */
  /* precipitation and source/sinks input into the surface.          */
  /* The total volume added to the cell in the 2D mode must be       */
  /* divided by the 3D time-step to get the flux over the 3D step.   */
  /* Note that waterss is the flow into the cell, wheras eta_rlx3d   */
  /* (m3) is the flow out of the cell, hence the minus sign.         */
  if (!(wincon->etarlx & NONE)) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      windat->waterss[c] -= wincon->eta_rlx3d[cs] / windat->dt;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add the contribution from precipitation and evaporation. The 2D */
  /* part was added in ss_water() before elevation was calculated in */
  /* the 2D mode. Note that the mass flux of salinity due to salt    */
  /* fluxes is handled as the surface boundary condition to vertical */
  /* diffusion.                                                      */
  if (wincon->saltflux & ORIGINAL) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      /* Convert from m s-1 to m3 s-1                                */
      windat->waterss[c] += (windat->nsfd[cs] * window->cellarea[cs]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the vertical velocity from the bottom to the surface. This  */
  /* is performed over all wet and open boundary tracer cells.       */
  /* Note: sur[] is lower of the surface coordinate before and after */
  /* the elevation is updated. Also, the cells to process vectors    */
  /* for tracers include the wet + open boundary cells (set in       */
  /* set_dz()).                                                      */
  /* Calculate velocity from the bottom upwards : this is consistent */
  /* with the original MECO formulation.                             */
  if (!wincon->sigma) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = bottom[cc];
      cs = window->m2d[c];
      windat->w[c] = windat->wbot[cs];
      fcbot = windat->w[c] * window->cellarea[cs];

      while (c > sur[cc]) {
        xp1 = window->xp1[c];
        yp1 = window->yp1[c];
        zp1 = window->zp1[c];
        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
        fctop = fcbot + windat->u1flux3d[c] - windat->u1flux3d[xp1] +
          windat->u2flux3d[c] - windat->u2flux3d[yp1] + windat->waterss[c];
        /* Velocity at cell top                                      */
        windat->w[zp1] = fctop / window->cellarea[cs];

	/*-----------------------------------------------------------*/
        /* Surface diagnostic checks.                                */
        if (window->wn == dw && zp1 == dc) {
          /* Check that the elevation change over the 3D timestep is */
          /* consistent with the divergence of the integrated 2D     */
          /* fluxes.                                                 */
          /* Change in elevation over the 3D timestep due to the     */
          /* divergence of 2D fluxes integrated over the 3D step.    */
	  /* Note : the fluxes here are the 2D fluxes integrated     */
	  /* over the 2D timestep.                                   */
          d1 = wincon->oldeta[dc] - (windat->u1flux[window->xp1[dc]] -
                                     windat->u1flux[dc] +
                                     windat->u2flux[window->yp1[dc]] -
                                     windat->u2flux[dc] -
				     windat->waterss[dc] * windat->dt) /
                                     window->cellarea[dc];

          /* Elevation difference shound be around 1e-17             */
          emstag(LDEBUG,"hd:vel3d:vel_w_update", "ele %f %e\n", windat->t / 86400,
                  windat->eta[dc] - d1);
          printf("ele %f %e\n", windat->t / 86400, windat->eta[cs] - d1);

          /* Check that the vertical integral of the 3D fluxes is    */
          /* equal to the 2D fluxes (dc must be a 2D coordinate).    */
          check_flux(window, windat, wincon, dw, dc);

          /* Calculate and print the continuity balance if required, */
	  /* i.e. the sum of horizontal divergence in the surface    */
	  /* layer plus the vertical flux into the surface layer     */
	  /* should equal the change in elevation over the timestep. */
          /* Output value should be around 1e-8.                     */
          hf =
            (windat->u1flux3d[window->xp1[zp1]] - windat->u1flux3d[zp1] +
             windat->u2flux3d[window->yp1[zp1]] -
             windat->u2flux3d[zp1]) * windat->dt;
          vf =
            (0.0 - windat->w[zp1]) * windat->dt * window->cellarea[cs] -
	    windat->waterss[zp1] * windat->dt;

          cnt =
            window->cellarea[cs] * (wincon->oldeta[cs] - windat->eta[cs]) -
            (hf + vf);
          emstag(LDEBUG,"hd:vel3d:vel_w_update", "cons %f : %e\n", windat->t / 86400, cnt);
          printf("cons %f : %e\n", windat->t / 86400, cnt);
        }
        /* Sub-surface diagnostic checks.                            */
        if (window->wn == dw && c == dc) {
          /* Calculate and print the continuity balance if required. */
          /* Output value should be around 1e-8.                     */
          double wtop = (c == cs) ? windat->wtop[cs] : windat->w[zp1];
          hf = windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
            windat->u2flux3d[window->yp1[c]] - windat->u2flux3d[c];
          vf = ((wtop - windat->w[c])) * window->cellarea[cs];
          emstag(LDEBUG,"hd:vel3d:vel_w_update", "cons sub %f %d : %e\n", windat->t / 86400,
                  window->s2k[c], hf + vf);
          printf("cons sub %f %d : %e\n", windat->t / 86400, window->s2k[c], hf + vf);
        }
	
        /* Transfer top flux to bottom for next cell                 */
        fcbot = fctop;
        c = zp1;
      }

      /* Calculations for those cells through which the surface      */
      /* moves; use value from wtop calculated above.                */
      while (c > nsur[cc]) {
        zp1 = window->zp1[c];
        windat->w[zp1] = windat->wtop[cs];
        c = zp1;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Calculate velocity from the surface downwards. Due to   */
  /* precision errors, if w is calculated upwards and wbot is set to */
  /* zero then wtop is not identically zero (usually around 1e-20).  */
  /* This can lead to conservation errors in the tracer calculation, */
  /* since the surface flux is assumed to equal zero (there is no    */
  /* ghost layer above the surface to save the surface flux in).     */
  /* To overcome this set wtop=0 and calculate w downards such that  */
  /* w(0)~1e-20. This bottom value can be used in the tracer         */
  /* vertical flux calculations, surface flux is equal to zero, thus */
  /* continutity is maintained.                                      */
  else {
    /* Calculate the rate of change of elevation to subtract from    */
    /* the vertical velocity.                                        */
    memset(deta, 0, window->sgsiz * sizeof(double));
    if (wincon->sigma) {
      for (cc = 1; cc <= wincon->vcs; cc++) {
        c = bottom[cc];
        cs = window->m2d[c];
        detadt = (windat->eta[cs] - wincon->oldeta[cs]) / windat->dt;
        while (c > sur[cc]) {
          zp1 = window->zp1[c];
          deta[c] = detadt * wincon->dz[c] * window->cellarea[cs];
          c = zp1;
        }
        deta[c] = detadt * wincon->dz[c] * window->cellarea[cs];
      }
    }
    /* Calculate the vertical velocity downwards                     */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = sur[cc];
      cs = window->m2d[c];
      fctop = 0.0;

      while (c <= bottom[cc]) {
        xp1 = window->xp1[c];
        yp1 = window->yp1[c];
        zm1 = window->zm1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
        fcbot = fctop + windat->u1flux3d[xp1] - windat->u1flux3d[c] +
          windat->u2flux3d[yp1] - windat->u2flux3d[c] + deta[c] +
          windat->waterss[c];

        /* Velocity at cell top                                      */
        windat->w[c] = fcbot / window->cellarea[cs];

        /* Diagnostic checks.                                        */
        if (window->wn == dw && cs == dc) {
          /* Calculate and print the continuity balance if required. */
          /* Output value should be around 1e-8.                     */
          double wtop =
            (c == cs) ? windat->wtop[cs] : windat->w[window->zp1[c]];
          hf =
            windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
            windat->u2flux3d[window->yp1[c]] - windat->u2flux3d[c] +
            deta[c];
          vf = ((wtop - windat->w[c])) * window->cellarea[cs];
          if (c == bottom[cc])
            emstag(LDEBUG,"hd:vel3d:vel_w_update", "cons sub %f %d : %e %e\n", windat->t / 86400,
                    window->s2k[c], hf + vf, windat->w[c]);
          else
            emstag(LDEBUG,"hd:vel3d:vel_w_update", "cons sub %f %d : %e\n", windat->t / 86400,
                    window->s2k[c], hf + vf);
        }

        /* Transfer top flux to bottom for next cell                 */
        fctop = fcbot;
        c = zm1;
      }
    }
  }

  /* Set the vertical velocity at dry faces                          */
  set_w_dry(window, windat, wincon);

  /* Set the open boundary condition if required                     */
  bdry_w(window, windat, wincon);

  /* Get the mean w velocity if required                             */
  /* If volume fluxes are being stored for the transport model using */
  /* the flux-form semi-lagrange (FFSL) advection scheme then the    */
  /* waterss array is stored in wm. This is to accomodate point      */
  /* sources. In the subsequent transport model run, with the FFSL   */
  /* scheme, the values from wmean are copied into waterss, and w is */
  /* re-calculated to ensure continuity.                             */
  /* This is only performed if point sources are included.           */
  /* In all other cases, mean values of w are calculated in wm.      */
  
  if (wincon->means & VEL3D && windat->wm) {
    double t = windat->dtf;
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      if (wincon->means & (VOLFLUX|PSSFLUX)) {
	windat->wm[c] = (windat->wm[c] * windat->meanc[cs] + 
	 	        (windat->waterss[c] * t / window->cellarea[cs])) /
                        (windat->meanc[cs] + t);
      } else {
	windat->wm[c] = (windat->wm[c] * windat->meanc[cs] + 
	 	         windat->w[c] * t)  / (windat->meanc[cs] + t);
      }
    }
  }
}

/* END vel_w_update()                                                */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Calculates the surface layer vertical velocity through continuity */
/* specifically for the flux-form semi-lagrange advection scheme in  */
/* the transport model.                                              */
/*-------------------------------------------------------------------*/
void ff_sl_w_update(geometry_t *window, /* Window geometry           */
		    window_t *windat,   /* Window data               */
		    win_priv_t *wincon  /* Window constants          */
  )
{
  int c, cc, cs;      /* Sparse coodinate / counter                  */
  int *bottom;        /* Bottom sparse coordinate                    */
  int *sur;           /* Minimum of surface coordinate               */
  int *nsur;          /* Surface sparse coordinate after the 2D mode */
  double fctop;       /* Flux at the top face of the cell            */

  double fcbot;       /* Flux at the bottom face of the cell         */
  int xp1, yp1;       /* Sparse coordinates east and south of cell c */
  int zp1, zm1;       /* Sparse coordinate below cell c              */
  int dw = 0;         /* Window for diagnostic continuity            */
  double detadt;      /* Rate of change of elevation                 */
  double *deta;       /* Elevation rate of change for SIGMA          */
  double *oldw;       /* w at start of step                          */
  double *Fx, *Fy;

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  nsur = wincon->i2;            /* Set in set_dz()                   */
  sur = wincon->i3;             /* Set in set_dz()                   */
  deta = wincon->w1;
  Fx = windat->u1flux3d;
  Fy = windat->u2flux3d;
  oldw = wincon->w2;
  memcpy(oldw, windat->w, window->sgsiz * sizeof(double));
  /* Uncomment if snapshots are required rather than cumulative
  if(windat->vol_cons) memset(windat->vol_cons, 0, window->sgsizS * sizeof(double)); */
  memset(windat->w, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across e1 and e2 OBC's          */
  reset_map_t(window);

  /*-----------------------------------------------------------------*/
  /* Calculate the vertical velocity at the surface and bottom       */
  /* boundaries.                                                     */
  if (wincon->momsc & WTOP_O4)
    vel_w_bounds_hiorder(window, windat, wincon);
  else
    vel_w_bounds(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the vertical velocity from the bottom to the surface. This  */
  /* is performed over all wet and open boundary tracer cells.       */
  /* Note: sur[] is lower of the surface coordinate before and after */
  /* the elevation is updated. Also, the cells to process vectors    */
  /* for tracers include the wet + open boundary cells (set in       */
  /* set_dz()).                                                      */
  /* Calculate velocity from the bottom upwards : this is consistent */
  /* with the original MECO formulation.                             */
  if (!wincon->sigma) {
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = bottom[cc];
      cs = window->m2d[c];
      windat->w[c] = windat->wbot[cs];
      fcbot = windat->w[c] * window->cellarea[cs];

      while (c > sur[cc]) {
        xp1 = window->xp1[c];
        yp1 = window->yp1[c];
        zp1 = window->zp1[c];
        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
        fctop = fcbot + Fx[c] - Fx[xp1] + Fy[c] - Fy[yp1] + 
          (windat->waterss[c] * window->cellarea[cs]);
        /* Velocity at cell top                                      */
        windat->w[zp1] = fctop / window->cellarea[cs];
        /* Transfer top flux to bottom for next cell                 */
        fcbot = fctop;
        c = zp1;
      }

      /* Calculations for those cells from the surface to the        */
      /* top of the grid.                                            */
      zp1 = window->zp1[c];
      while (c != zp1) {
	xp1 = window->xp1[c];
	yp1 = window->yp1[c];
        fctop = fcbot + Fx[c] - Fx[xp1] + Fy[c] - Fy[yp1] + 
          (windat->waterss[c] * window->cellarea[cs]);
        windat->w[zp1] = fctop / window->cellarea[cs];
        /* Transfer top flux to bottom for next cell                 */
        fcbot = fctop;
        c = zp1;
        zp1 = window->zp1[c];
      }
    }

    /* Only update the vertical velocity if the new velocity does    */
    /* not violate the Lipschitz condition. If it does, then report  */
    /* the location in vol_cons[c] = 1 and continue using the        */
    /* original velocity (conservation may be violated in the column */
    /* in this case).                                                */
    if (wincon->conserve & CONS_WS) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	double iratio = 200;
	double sf = 0.8;
	double minval = 1e-10;
	double ws, dto, dtn;
	c = wincon->s1[cc];
	zp1 = window->zp1[c];
	cs = window->m2d[c];
	ws = (cc <= wincon->vcs) ? minval : fabs(oldw[zp1] - oldw[c]);
	ws = max(minval, ws);
	dto = sf * fabs(wincon->dz[c] * wincon->Ds[cs] / ws);
	ws = (cc <= wincon->vcs) ? minval : fabs(windat->w[zp1] - windat->w[c]);
	ws = max(minval, ws);
	dtn = sf * fabs(wincon->dz[c] * wincon->Ds[cs] / ws);
	if (dto > windat->dt && dtn < windat->dt) {
	  if(windat->vol_cons) windat->vol_cons[cs] = 1;
	}
      }
      for (cc = 1; cc <= wincon->vcs; cc++) {
	cs = c = sur[cc];
	if (windat->vol_cons[cs] == 1.0) {
	  while (c != window->zm1[c]) {
	    windat->w[c] = oldw[c];
	    c = window->zm1[c];
	  }
	}
      }
    } else {
      for (cc = 1; cc <= wincon->vc; cc++) {
	double cons, minval = 1e-10;
	c = wincon->s1[cc];
	cs = window->m2d[c];
        xp1 = window->xp1[c];
        yp1 = window->yp1[c];
        zp1 = window->zp1[c];
        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
        cons = Fx[xp1] - Fx[c] + Fy[yp1] - Fy[c] + (windat->w[zp1] - 
		windat->w[c] - windat->waterss[c]) * window->cellarea[cs];
	if (fabs(cons) > minval && windat->vol_cons) windat->vol_cons[cs] = 1;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Calculate velocity from the surface downwards. Due to   */
  /* precision errors, if w is calculated upwards and wbot is set to */
  /* zero then wtop is not identically zero (usually around 1e-20).  */
  /* This can lead to conservation errors in the tracer calculation, */
  /* since the surface flux is assumed to equal zero (there is no    */
  /* ghost layer above the surface to save the surface flux in).     */
  /* To overcome this set wtop=0 and calculate w downards such that  */
  /* w(0)~1e-20. This bottom value can be used in the tracer         */
  /* vertical flux calculations, surface flux is equal to zero, thus */
  /* continutity is maintained.                                      */
  else {
    /* Calculate the rate of change of elevation to subtract from    */
    /* the vertical velocity.                                        */
    memset(deta, 0, window->sgsiz * sizeof(double));
    if (wincon->sigma) {
      for (cc = 1; cc <= wincon->vcs; cc++) {
        c = bottom[cc];
        cs = window->m2d[c];
        detadt = (windat->eta[cs] - wincon->oldeta[cs]) / windat->dt;
        while (c > sur[cc]) {
          zp1 = window->zp1[c];
          deta[c] = detadt * wincon->dz[c] * window->cellarea[cs];
          c = zp1;
        }
        deta[c] = detadt * wincon->dz[c] * window->cellarea[cs];
      }
    }
    /* Calculate the vertical velocity downwards                     */
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = sur[cc];
      cs = window->m2d[c];
      fctop = 0.0;

      while (c <= bottom[cc]) {
        xp1 = window->xp1[c];
        yp1 = window->yp1[c];
        zm1 = window->zm1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
        fcbot = fctop + windat->u1flux3d[xp1] - windat->u1flux3d[c] +
          windat->u2flux3d[yp1] - windat->u2flux3d[c] + deta[c] +
          windat->waterss[c];

        /* Velocity at cell top                                      */
        windat->w[c] = fcbot / window->cellarea[cs];

        /* Transfer top flux to bottom for next cell                 */
        fctop = fcbot;
        c = zm1;
      }
    }
  }

  /* Set the vertical velocity at dry faces                          */
  /* set_w_dry(window, windat, wincon); */

  /* Set the open boundary condition if required                     */
  /* bdry_w(window, windat, wincon); */

}

/* END ff_sl_w_update()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets vertical velocity above the surface to wtop. This aids       */
/* stability in the vertical momentum advection when cells have one  */
/* face dry.                                                         */
/*-------------------------------------------------------------------*/
void set_w_dry(geometry_t *window,    /* Window geometry             */
	       window_t *windat,      /* Window data                 */
	       win_priv_t *wincon     /* Window constants            */
	       )
{
  int c, cc, cs, zp1;

  if (!wincon->sigma) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->nsur_t[cc];
      cs = window->m2d[c];
      zp1 = window->zp1[c];
      while (c != zp1) {
        windat->w[zp1] = windat->wtop[cs];
        c = zp1;
	zp1 = window->zp1[c];
      }
    }
  }
}

/* END set_w_dry()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the surface and bottom vertical velocity     */
/*-------------------------------------------------------------------*/
void vel_w_bounds(geometry_t *window, /* Window geometry             */
                  window_t *windat,   /* Window data                 */
                  win_priv_t *wincon  /* Window constants            */
  )
{
  int c, cc;                    /* Sparse coodinate / counter        */
  int cs;                       /* 2D sparse coordinate              */
  int cb;                       /* Bottom sparse coordinate          */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1     */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1     */
  int *bottom;                  /* Bottom sparse coordinate          */
  double eta_l;                 /* Surface elevation at i-1          */
  double eta_r;                 /* Surface elevation at i+1          */
  double eta_b;                 /* Surface elevation at j-1          */
  double eta_f;                 /* Surface elevation at j+1          */
  double *u1top;                /* Surface e1 velocity               */
  double *u2top;                /* Surface e2 velocity               */
  double *u1bot;                /* Bottom e1 velocity                */
  double *u2bot;                /* Bottom e2 velocity                */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  u1top = wincon->d1;
  u1bot = wincon->d3;
  u2top = wincon->d2;
  u2bot = wincon->d4;
  memset(windat->wtop, 0, window->sgsizS * sizeof(double));
  memset(windat->wbot, 0, window->sgsizS * sizeof(double));
  memset(u1top, 0, window->sgsizS * sizeof(double));
  memset(u2top, 0, window->sgsizS * sizeof(double));
  memset(u1bot, 0, window->sgsizS * sizeof(double));
  memset(u2bot, 0, window->sgsizS * sizeof(double));

  /* SIGMA : zero velocity at the sigma boundaries                   */
  if (wincon->sigma)
    return;

  /* Set the surface velocity. Note: this cannot be retreived from   */
  /* the u1 array since neighbours (eg xp1) are required and with a  */
  /* cell drying the xp1 neighbour may not be the surface.           */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    cs = window->m2d[c];
    u1top[cs] = windat->u1[c];
    u2top[cs] = windat->u2[c];
  }
  for (cc = 1; cc <= window->b2_e1; cc++) {
    cb = window->bot_e1[cc];  /* 3D bottom coordinate                */
    cs = window->m2d[cb];
    u1bot[cs] = windat->u1[cb];
  }
  for (cc = 1; cc <= window->b2_e2; cc++) {
    cb = window->bot_e2[cc];  /* 3D bottom coordinate                */
    cs = window->m2d[cb];
    u2bot[cs] = windat->u2[cb];
  }

  /* Note : should relocate this call so that copies from the master */
  /* are made for multiple windows.                                  */
#if GLOB_BC
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                     geom->bin2, 0);
#endif

#if !GLOB_BC
  set_lateral_bc_eta(windat->eta, window->nbptS, window->bpt, window->bin,
                     window->bin2, 0);
#endif

  /* In the linear case wtop and detadt are the same, whereas in the */
  /* non-linear case, they are related by terms involving the        */
  /* surface slope and surface horizontal velocities.                */
  memcpy(windat->wtop, windat->detadt, window->sgsizS * sizeof(double));
  if (!wincon->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.                                   */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->i3[cc];
    cs = window->m2d[c]; 
    xp1 = window->xp1[cs];
    xm1 = window->xm1[cs];
    yp1 = window->yp1[cs];
    ym1 = window->ym1[cs];

    /* Calculate the surface vertical velocity                       */
    eta_l = windat->eta[xm1];
    eta_r = windat->eta[xp1];
    eta_b = windat->eta[ym1];
    eta_f = windat->eta[yp1];

    windat->wtop[cs] += (u1top[cs] + u1top[xp1]) * (eta_r - eta_l) /
      (4.0 * window->h1acell[cs]) +
      (u2top[cs] + u2top[yp1]) * (eta_f - eta_b) /
      (4.0 * window->h2acell[cs]);

    /* Calculate the bottom vertical velocity                        */
    windat->wbot[cs] = -(u1bot[cs] + u1bot[xp1]) *
      window->dHde1[cs] / (2.0 * window->h1acell[cs]) -
      (u2bot[cs] + u2bot[yp1]) *
      window->dHde2[cs] / (2.0 * window->h2acell[cs]);

    /* if (cs == 312) 
       printf("wtop1: t = %f, c = %d, wtop = %f\n",windat->t/86400, c, windat->wtop[cs]); */
  }
}

/* END vel_w_bounds()                                                */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to calculate the surface and bottom vertical velocity     */
/*-------------------------------------------------------------------*/
void vel_w_bounds_hiorder(geometry_t *window, /* Window geometry     */
			  window_t *windat,   /* Window data         */
			  win_priv_t *wincon  /* Window constants    */
			  )
{
  int c, cc;                    /* Sparse coodinate / counter        */
  int cs;                       /* 2D sparse coordinate              */
  int cb;                       /* Bottom sparse coordinate          */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1     */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1     */
  int xp2, xm2;                 /* Sparse coordinate at i+2, i-2     */
  int yp2, ym2;                 /* Sparse coordinate at j+2, j-2     */
  int *bottom;                  /* Bottom sparse coordinate          */
  double *u1bot;                /* Bottom e1 velocity                */
  double *u2bot;                /* Bottom e2 velocity                */
  double *detadt;               /* eta change over the 3D timestep   */
  double detade1;               /* e1 gradient of eta                */
  double detade2;               /* e2 gradient of eta                */
  double u1c, u2c;              /* Cell centered u1 and u2 velocity  */
  double dh1, dh2;              /* Grid spacings                     */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  detadt = wincon->d1;
  u1bot = wincon->d2;
  u2bot = wincon->d4;
  memset(windat->wtop, 0, window->sgsizS * sizeof(double));
  memset(windat->wbot, 0, window->sgsizS * sizeof(double));
  memset(u1bot, 0, window->sgsizS * sizeof(double));
  memset(u2bot, 0, window->sgsizS * sizeof(double));
  memset(detadt, 0, window->sgsizS * sizeof(double));

  /* SIGMA : zero velocity at the sigma boundaries                   */
  if (wincon->sigma)
    return;

  /* Set the change in elevation over the 3D timestep                */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    detadt[c] = (windat->eta[c] - wincon->oldeta[c]) / windat->dt;
  }
  /* Set the bottom velocity.                                        */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    cb = window->bot_e1[cc];  /* 3D bottom coordinate                */
    cs = window->m2d[cb];
    u1bot[cs] = windat->u1[cb];
  }
  for (cc = 1; cc <= window->b2_e2; cc++) {
    cb = window->bot_e2[cc];  /* 3D bottom coordinate                */
    cs = window->m2d[cb];
    u2bot[cs] = windat->u2[cb];
  }

  /* Note : should relocate this call so that copies from the master */
  /* are made for multiple windows.                                  */
#if GLOB_BC
  set_lateral_bc_eta(master->eta, geom->nbptS, geom->bpt, geom->bin,
                     geom->bin2, 0);
#endif

#if !GLOB_BC
  set_lateral_bc_eta(windat->eta, window->nbptS, window->bpt, window->bin,
                     window->bin2, 0);
#endif

  /* In the linear case wtop and detadt are the same, whereas in the */
  /* non-linear case, they are related by terms involving the        */
  /* surface slope and surface horizontal velocities.                */
  memcpy(windat->wtop, detadt, window->sgsizS * sizeof(double));
  if (!wincon->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.                                   */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->i3[cc];
    cs = window->m2d[c]; 
    xp1 = window->xp1[cs];
    xm1 = window->xm1[cs];
    xp2 = window->xp1[xp1];
    xm2 = window->xm1[xm1];
    yp1 = window->yp1[cs];
    ym1 = window->ym1[cs];
    yp2 = window->yp1[yp1];
    ym2 = window->ym1[yp2];

    /* Calculate the surface vertical velocity using sea level       */
    /* gradients and cell centered velocity computed with 4th order  */
    /* approximations.                                               */
    dh1 = window->h1au1[cs] + window->h1au1[xp1];
    dh2 = dh1 + window->h1au1[xm1] + window->h1au1[xp2];
    detade1 = 4.0 * (windat->eta[xp1] - windat->eta[xm1]) / (3.0 * dh1) -
      (windat->eta[xp2] - windat->eta[xm2]) / (3.0 * dh2);
    dh1 = window->h2au2[cs] + window->h2au2[yp1];
    dh2 = dh1 + window->h2au2[ym1] + window->h2au2[yp2];
    detade2 = 4.0 * (windat->eta[yp1] - windat->eta[ym1]) / (3.0 * dh1) -      
      (windat->eta[yp2] - windat->eta[ym2]) / (3.0 * dh2);
    u1c = (7.0 * (windat->u1[cs] + windat->u1[xp1]) - 
	   (windat->u1[xp2] + windat->u1[xm1])) / 12.0;
    u2c = (7.0 * (windat->u2[cs] + windat->u2[yp1]) - 
	   (windat->u2[yp2] + windat->u2[ym1])) / 12.0;
    windat->wtop[cs] += u1c * detade1 + u2c * detade2;

    /* Calculate the bottom vertical velocity                        */
    windat->wbot[cs] = -(u1bot[cs] + u1bot[xp1]) *
      window->dHde1[cs] / (2.0 * window->h1acell[cs]) -
      (u2bot[cs] + u2bot[yp1]) *
      window->dHde2[cs] / (2.0 * window->h2acell[cs]);
  }
}

/* END vel_w_bounds_hiorder()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the w velocities on all open boundaries in a   */
/* given window.                                                     */
/*-------------------------------------------------------------------*/
void bdry_w(geometry_t *window, /* Window geometry                   */
	    window_t *windat,   /* Window data                       */
	    win_priv_t *wincon  /* Window constants                  */
	    )
{
  int n;                              /* Counters                    */
  open_bdrys_t **open = window->open; /* Open boundary structure     */

  /*-----------------------------------------------------------------*/
  /* Set the u1 velocity where u1 is normal to the boundary          */
  for (n = 0; n < window->nobc; n++) {

    set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_t,
            open[n]->obc_t, open[n]->oi1_t, open[n]->oi2_t,
            open[n]->cyc_t, windat->w, windat->w, 
	    windat->w, open[n]->bcond_w, windat->dtf, NULL, NULL, 0.0, 0);

    /* Average the sea level over intersecting OBCs if required      */
    if (open[n]->options & OP_OBCCM)
      average_OBC_corner(window, open[n], windat->w, 1);
  }
}


/* END bdry_w()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implicit vertical diffusion scheme                                */
/*-------------------------------------------------------------------*/
int implicit_vdiff(geometry_t *window, /* Window geometry           */
                    window_t *windat,   /* Window data               */
                    win_priv_t *wincon, /* Window constants          */
                    double *var,  /* Array of variable values to diffuse */
                    double *nvar, /* Array of variable values to diffuse */
                    double *Kzin,  /* Mixing coefficient values      */
                    double *dzin,  /* Cell thicknesses               */
                    double *fb,    /* Flux out of bottom             */
                    double *ft,    /* Flux out of top                */
                    int *ctp,      /* Surface cells to process       */
                    int *cbm,      /* Bottom cells to process        */
                    int vcs,       /* Last index of surface cells    */
                    double *depth, /* SIGMA : Water depth            */
                    double *scale,  /* SIGMA : (depth) scaling for velocity */
		    double *w
  )
{
  int c, k;                  /* Sparse coordinate                    */
  int cs, ks, ks1;           /* Surface sparse coordinate            */
  int cb, kb;                /* Bottom sparse coordinate             */
  int zm1;                   /* Sparse cell below c                  */
  int zp1;                   /* Sparse cell above c                  */
  int cc;                    /* Sparse coordinate counter            */
  int c2;                    /* 2D cell corresponding to 3D location */
  double dt = windat->dt;    /* Time step for the window             */
  double dzdt;               /* dz / dt                              */
  double div;                /* Constant                             */
  double dzs;                /* Thickness of the surface layer       */
  double vel1, vel2;         /* Surface and zm1 velocities, var[]    */
  int thinf;                 /* Thin layer flag                      */

  /*-----------------------------------------------------------------*/
  /* Set pointers.                                                   */
  /* Note: the 3D  work  arrays  wincon->w#  could be  used  for the */
  /* dummy arrays below, but  execution speed is considerably faster */
  /* when  work array  access  is  sequential in  memory, hence  the */
  /* mapping to  a  contiguous  vertical  1D work array wincon->v#.  */
  double *Cm1 = wincon->v1;
  double *C = wincon->v2;
  double *Cp1 = wincon->v3;
  double *rhs = wincon->v4;
  double *sol = wincon->v5;
  double *ud = wincon->v6;
  double *dzface = wincon->v7;
  double *dvar = wincon->v8;
  double *dz = wincon->v9;
  double *Kz = wincon->v10;
  double *dzf = wincon->v11;

  /* Loop through the surface cells in this window                   */
  for (cc = 1; cc <= vcs; cc++) {

    cs = c = ctp[cc];     /* Set cs to the surface sparse coordinate */
    cb = cbm[cc];         /* Bottom sparse coordinate                */
    c2 = window->m2d[c];  /* 2D sparse location corresponding to c   */
    zm1 = window->zm1[c];
    dzs = dzin[c];
    ks = ks1 = window->s2k[cs];
    kb = window->s2k[cb];
    thinf = 0;

    /*---------------------------------------------------------------*/
    /* Merge thin  layers with layer  below if  required by  setting */
    /* velocity in the  layers cs and zm1 to a mass weighted mean of */
    /* these  two layers, then  assuming zm1  is the layer below the */
    /* thin layer.                                                   */
    if (wincon->thin_merge && dzin[cs] < wincon->hmin) {
      thinf = 1;
      cs = zm1;
      ks--;
      zm1 = window->zm1[zm1];
      vel1 = var[c];
      vel2 = var[cs];
      /* Thin layer lies in a single layer -> continue               */
      if (cs == zm1)
        continue;
      var[c] = var[cs] =
        (var[c] * dzs + var[cs] * dzin[cs]) / (dzin[cs] + dzs);
      dzin[cs] += dzs;
    }

    /*---------------------------------------------------------------*/
    /* Single layer  case  (i.e. surface lies  in the  bottom layer) */
    if (zm1 == window->zm1[zm1] || !dzin[zm1]) {
      nvar[cs] +=
        windat->dt * (fb[c2] - ft[c2]) / max(dzin[cs], wincon->hmin);
      if (c != cs) {
        var[c] = var[cs];
        dzin[cs] -= dzs;
      }
      continue;
    }

    /*---------------------------------------------------------------*/
    /* Calculate constants  required from  the surface  to  bottom.  */
    /* The sparse coordinate, c, lies in  the sediment at the end of */
    /* the loop.                                                     */
    c = cs;
    for (k = ks; k >= kb; k--) {
      zm1 = window->zm1[c];
      dzface[k] = 0.5 * (dzin[zm1] + dzin[c]);
      /* SIGMA : Precondition mixing coefficient and the cell for    */
      /* the sigma calculation.                                      */
      dz[k] = dzin[c] * depth[c2];
      Kz[k] = Kzin[c] / depth[c2];
      dzf[k] = dzin[zm1] / (dzin[zm1] + dzin[c]);
      c = zm1;
    }

    /*---------------------------------------------------------------*/
    /* Set up tri-diagonal set of equations.                         */
    /* Bottom layer.                                                 */
    zm1 = c;
    c = cb;
    zp1 = window->zp1[c];
    dzdt = dz[kb] / dt;

    Cm1[kb] = 0.0;
    Cp1[kb] = -Kz[kb + 1] / dzface[kb + 1] + w[c] * dzf[kb+1];
    C[kb] = dzdt - Cp1[kb] + w[c] - w[zm1];
    rhs[kb] = dzdt * var[cb] + fb[c2];

    /* Mid-water layers                                              */
    zm1 = c;
    c = zp1;
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / dt;
      Cm1[k] = -Kz[k] / dzface[k] + w[zm1] * dzf[k];
      Cp1[k] = -Kz[k + 1] / dzface[k + 1] + w[c] * dzf[k+1];
      C[k] = dzdt - Cm1[k] - Cp1[k] + w[c];
      Cm1[k] -= w[zm1];
      rhs[k] = dzdt * var[c];
      zm1 = c;
      c = window->zp1[c];
    }

    /* Surface layer                                                 */
    dzdt = dz[ks] / dt;
    Cm1[ks] = -Kz[ks] / dzface[ks] + w[zm1] * dzf[ks];
    Cp1[ks] = 0.0;
    C[ks] = dzdt - Cm1[ks] + w[c];
    rhs[ks] = dzdt * var[cs] - ft[c2];
    Cm1[ks] -= w[zm1];

    /*---------------------------------------------------------------*/
    /* Solve tridiagonal system                                      */
    div = C[kb];
    sol[kb] = rhs[kb] / div;
    for (k = kb + 1; k <= ks; k++) {
      ud[k] = Cp1[k - 1] / div;
      div = C[k] - Cm1[k] * ud[k];
      if (div == 0.0) {
        hd_quit_and_dump("Momentum diffusion : zero divisor\n");
	return(1);
      }
      sol[k] = (rhs[k] - Cm1[k] * sol[k - 1]) / div;
    }
    dvar[ks] = sol[ks] - var[cs];

    c = cs;
    for (k = ks - 1; k >= kb; k--) {
      c = window->zm1[c];
      sol[k] -= ud[k + 1] * sol[k + 1];
      dvar[k] = sol[k] - var[c];
    }

    /*---------------------------------------------------------------*/
    /* Reset the surface dz for thin layers                          */
    c = ctp[cc];
    if (c != cs) {
      dzin[cs] -= dzs;
      dvar[ks1] = dvar[ks];
    }

    /*---------------------------------------------------------------*/
    /* Update the variable                                           */
    c = ctp[cc];
    /* Restore the backward velocities to their un-merged values if  */
    /* merging occurred. If the bottom lies in the layer below the   */
    /* surface, the u1b values are used in vdiff_u2() to get the     */
    /* bottom stress, and the un-merged values are to be used.       */
    if (thinf) {
      var[c] = vel1;
      var[window->zm1[c]] = vel2;
    }
    for (k = ks1; k >= kb; k--) {
      nvar[c] += (dvar[k] * scale[c2]);
      c = window->zm1[c];
    }
  }
  return(0);
}

/* END implicit_vdiff()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set variables over blended zones                       */
/*-------------------------------------------------------------------*/
void blend_vel(geometry_t *window, 
	       window_t *windat, 
	       win_priv_t *wincon,
	       int mode,
	       double *vel)
{
  int n, bn = 0, c, cc, c3;
  int nb;
  int *p1, *m1, *index, i, i1, i2, c1, c2;
  double frac, v1, v2;
  double *dx, *dy;
  blend_t *blend;

  if (mode & (U12D|U13D|U1_VH|BLM1|BLR1))
    bn = wincon->nbl1;
  else if (mode & (U22D|U23D|U2_VH|BLM2|BLR2))
    bn = wincon->nbl2;
  else
    return;

  for (n = 0; n < bn; n++) {
    if (mode & (U12D|U13D|U1_VH|BLM1|BLR1)) {
      blend = wincon->ble1[n];
      if (mode & U12D)
	nb = blend->ne1bS;
      else
	nb = blend->ne1b;
      if (mode == U1_VH) {
	dx = window->h1au1;
	dy = window->h2au1;
      }
      p1 = window->xp1;
      m1 = window->xm1;
      index = window->s2i;
    } else if (mode & (U22D|U23D|U2_VH|BLM2|BLR2)) {
      blend = wincon->ble2[n];
      if (mode & U22D)
	nb = blend->ne2bS;
      else
	nb = blend->ne2b;
      if (mode == U2_VH) {
	dx = window->h1au2;
	dy = window->h2au2;
      }
      p1 = window->yp1;
      m1 = window->ym1;
      index = window->s2j;
    }

    if (mode & (U12D|U22D|U13D|U23D)) {
      for (cc = 1; cc <= nb; cc++) {
	c = blend->c[cc];
	c1 = blend->cs[cc];
	c2 = blend->ce[cc];
	i = index[c];
	i1 = index[c1] - 1;
	i2 = index[c2] + 1;
	v1 = vel[blend->os[cc]];
	v2 = vel[blend->oe[cc]];
	frac = (v2 - v1) / (double)(i2 - i1);
	vel[c] = frac * (double)(i - i1) + v1;
      }
    } else if (mode & (BLM1|BLM2)) {
      for (cc = 1; cc <= nb; cc++) {
	c = blend->c[cc];
	p1[c] = m1[c] = c;
      }
    } else if (mode & (BLR1|BLR2)) {
      for (cc = 1; cc <= nb; cc++) {
	c = blend->c[cc];
	p1[c] = blend->p[cc];
	m1[c] = blend->m[cc];;
      }
    } else {
      double vm, sfact = 0.9;
      int cs, im;
      int vmode = 2;   /* 0 = max, 1 = zero, 2 = original */
      for (cc = 1; cc <= nb; cc++) {
	c = blend->c[cc];
	cs = window->m2d[c];
	c1 = blend->cs[cc];
	c2 = blend->ce[cc];
	i = index[c];
	i1 = index[c1] - 1;
	i2 = index[c2] + 1;
	v1 = vel[blend->os[cc]];
	v2 = vel[blend->oe[cc]];
	im = (i2 + i1) / 2;
	vm = 1.0 / (dx[cs] * dx[cs]) + 1.0 / (dy[cs] * dy[cs]);
	vm = sfact / (4.0 * vm * windat->dtb);
	if (vmode == 0) {
	  if (i <= im) {
	    frac = (v1 - vm) / (double)(i1 - im);
	    vel[c] = frac * (double)(i - i1) + v1;
	  } else {
	    frac = (v2 - vm) / (double)(i2 - im);
	    vel[c] = frac * (double)(i - i2) + v2;
	  }
	} else if (vmode == 1)
	  vel[c] = 0.0;
      }
    }
  }
}

/* END blend_vel()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes geostrophically balanced velocity                        */
/*-------------------------------------------------------------------*/
void calc_geostrophic(geometry_t *window, 
		      window_t *windat, 
		      win_priv_t *wincon)
{
  int c, cs, cc, n;

  reset_map_t_all(window);

  cells2process_e1(window, windat, wincon);
  precalc_u1(window, windat, wincon);
  set_map_e1(window);
  pressure_u1(window, windat, wincon);
  extract_u1_3d(window, windat, wincon);

  cells2process_e2(window, windat, wincon);
  precalc_u2(window, windat, wincon);
  set_map_e2(window);
  pressure_u2(window, windat, wincon);
  extract_u2_3d(window, windat, wincon);

  for (cc = 1; cc < window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    windat->u1[c] = windat->nu2[c] / (windat->dt * wincon->u1c5[cs]);
  }
  for (cc = 1; cc < window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    windat->u2[c] = windat->nu1[c] / (windat->dt * wincon->u2c5[cs]);
  }

  /* Save the velocities for boundary ramping                        */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (cc = 1; cc <= open->no3_e1; cc++) {
      c = (open->ocodex & R_EDGE) ? open->obc_t[cc] : open->obc_e1[cc];
      open->u1d[cc] = windat->u1[c];
      windat->u1[open->obc_e1[cc]] = windat->u1[c];
    }
    for (cc = open->no3_e1 + 1; cc <=open->to3_e1; cc++) {
      c = open->obc_e1[cc];
      if (open->ocodey & B_EDGE) c = open->nmap[c];
      open->u1d[cc] = windat->u1[c];
    }
    for (cc = 1; cc <= open->no3_e2; cc++) {
      c = (open->ocodey & F_EDGE) ? open->obc_t[cc] : open->obc_e2[cc];
      open->u2d[cc] = windat->u2[c];
      windat->u2[open->obc_e2[cc]] = windat->u2[c];
    }
    for (cc = open->no3_e2 + 1; cc <=open->to3_e2; cc++) {
      c = open->obc_e2[cc];
      if (open->ocodex & L_EDGE) c = open->nmap[c];
      open->u2d[cc] = windat->u2[c];
    }
  }
  set_map_t(window);
}

/* END calc_geostrophic()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the geostropic current at a water column                  */
/*-------------------------------------------------------------------*/
void calc_geostrophic_c(geometry_t *window, 
			window_t *windat, 
			win_priv_t *wincon,
			int cs,
			double *vel,
			int mode
			)
{
  int c, ci, zm1;
  int c2 = cs;
  int *m1, *p1;
  double *dz, *md, dd = 0.0;

  if (mode) {
    m1 = window->ym1;
    p1 = window->yp1;
    dz = windat->dzu2;
    md = wincon->mdy;
  } else {
    m1 = window->xm1;
    p1 = window->xp1;
    dz = windat->dzu1;
    md = wincon->mdx;
  }

  if (window->wgst[cs]) c2 = m1[cs];
  ci = c2;
  c = cs;

  zm1 = window->zm1[ci];
  while (ci != zm1) {
    dd += wincon->g * md[c2] * md[c2] * (windat->dens[c] - windat->dens[m1[c]]) * dz[ci];
    vel[c] = windat->dtf * wincon->u2c6[c2] * dd / 
      (0.5 * (windat->dens[c] + windat->dens[m1[c]]));
    ci = zm1;
    zm1 = window->zm1[ci];
    c = (c2 == cs) ? ci : p1[ci];
  }
}

/* END calc_geostrophic_c()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the geostropic current at a water column                  */
/*-------------------------------------------------------------------*/
void calc_geostrophic_obc(geometry_t *window, 
			  window_t *windat, 
			  win_priv_t *wincon,
			  open_bdrys_t *open,
			  double *vel,
			  int mode
			  )
{
  int ce;
  int cc, c, cs, ci, c2;
  int m1, *obc;
  double *dz, *md;

  if (mode) {
    ce = open->no3_e2;
    obc = open->obc_e1;
    dz = windat->dzu2;
    md = wincon->mdy;
  } else {
    ce = open->no3_e1;
    obc = open->obc_e1;
    dz = windat->dzu1;
    md = wincon->mdx;
  }

  memset(open->dum, 0, (open->no3_t + 1) * sizeof(double));
  for (cc = 1; cc <= ce; cc++) {
    c = obc[cc];
    cs = window->m2d[c];
    ci = open->obc_t[cc];
    c2 = window->m2d[ci];
    m1 = (open->ocodex & L_EDGE || open->ocodey & B_EDGE) ? open->ogc_t[cc] : open->nmap[c];
    open->dum[cs] += wincon->g * md[c2] * md[c2] * (windat->dens[c] - windat->dens[m1]) * dz[ci];
    vel[c] = windat->dtf * wincon->u2c6[c2] * open->dum[cs] / 
      (0.5 * (windat->dens[c] + windat->dens[m1]));
  }
}

/* END calc_geostrophic_obc()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void mom_balance(geometry_t *window, 
		 window_t *windat, 
		 win_priv_t *wincon)
{
  int c, cc;
  double v, v1, v2, v3, v4, v5, v6, vm, vt;

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    v1 = sqrt(windat->u1_adv[c] * windat->u1_adv[c] +
	      windat->u2_adv[c] * windat->u2_adv[c]);
    v2 = sqrt(windat->u1_hdif[c] * windat->u1_hdif[c] +
	      windat->u2_hdif[c] * windat->u2_hdif[c]);
    v3 = sqrt(windat->u1_vdif[c] * windat->u1_vdif[c] +
	      windat->u2_vdif[c] * windat->u2_vdif[c]);
    v4 = sqrt(windat->u1_cor[c] * windat->u1_cor[c] +
	      windat->u2_cor[c] * windat->u2_cor[c]);
    v5 = sqrt(windat->u1_btp[c] * windat->u1_btp[c] +
	      windat->u2_btp[c] * windat->u2_btp[c]);
    v6 = sqrt(windat->u1_bcp[c] * windat->u1_bcp[c] +
	      windat->u2_bcp[c] * windat->u2_bcp[c]);
    vt = v1 + v2 + v3 + v4 + v5 + v6;
    v = 1.0;                    /* 1 = advection */
    vm = v1;
    v = (v2 > vm) ? 2.0 : v;    /* 2 = horizontal diffusion */
    vm = (v2 > vm) ? v2 : vm;
    v = (v3 > vm) ? 4.0 : v;    /* 4 = vertical diffusion */
    vm = (v3 > vm) ? v3 : vm;
    v = (v4 > vm) ? 8.0 : v;    /* 8 = Coriolis */
    vm = (v4 > vm) ? v4 : vm;
    v = (v5 > vm) ? 16.0 : v;   /* 16 = barotropic pressure */
    vm = (v5 > vm) ? v5 : vm;
    v = (v6 > vm) ? 32.0 : v;   /* 32 = baroclinic pressure */
    vm = (v6 > vm) ? v6 : vm;
    windat->mom_bal[c] = v;
  }
}

/* END mom_balance()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the stokes velocity profile at e1 faces                   */
/*-------------------------------------------------------------------*/
void get_sdc_e1(geometry_t *window, 
		   window_t *windat, 
		   win_priv_t *wincon)
{
  int c, cc, cs;
  double w;                                      /* Angular freq     */
  double k;                                      /* Wave number      */
  double uso, u1, u2, depth;
  double ramp = (wincon->rampf & STOKES) ? windat->rampval : 1.0;

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    cs = window->m2d[c];
    w = (windat->wave_period[cs]) ? 2.0 * PI / windat->wave_period[cs] : 0.0;
    k = w * w / wincon->g;
    /*
    u1 = 0.5 * (windat->wave_ste1[cs] * windat->wave_ste1[window->xm1[cs]]);
    u2 = 0.5 * (windat->wave_ste2[cs] * windat->wave_ste2[window->ym1[cs]]);
    uso = sqrt(u1 * u1 + u2 * u2);
    */
    uso = ramp * windat->wave_ste1[cs];
    while (c != window->zm1[c]) {
      depth = min(0.0, window->cellz[c]);
      wincon->w1[c] = uso * exp(2.0 * k * depth);
      c = window->zm1[c];
    }
  }
}

/* END get_sdc_e1()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the stokes velocity profile at e1 faces                   */
/*-------------------------------------------------------------------*/
void get_sdc_e2(geometry_t *window, 
		   window_t *windat, 
		   win_priv_t *wincon)
{
  int c, cc, cs;
  double w;                                      /* Angular freq     */
  double k;                                      /* Wave number      */
  double uso, u1, u2, depth;
  double ramp = (wincon->rampf & STOKES) ? windat->rampval : 1.0;

  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    cs = window->m2d[c];
    w = (windat->wave_period[cs]) ? 2.0 * PI / windat->wave_period[cs] : 0.0;
    k = w * w / wincon->g;
    /*
    u1 = 0.5 * (windat->wave_ste1[cs] * windat->wave_ste1[window->xm1[cs]]);
    u2 = 0.5 * (windat->wave_ste2[cs] * windat->wave_ste2[window->ym1[cs]]);
    uso = sqrt(u1 * u1 + u2 * u2);
    */
    uso = ramp * windat->wave_ste2[cs];
    while (c != window->zm1[c]) {
      depth = min(0.0, window->cellz[c]);
      wincon->w1[c] = uso * exp(2.0 * k * depth);
      c = window->zm1[c];
    }
  }
}

/* END get_sdc_e2()                                                  */
/*-------------------------------------------------------------------*/
