/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/momentum/vel3d.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vel3d.c 6741 2021-03-30 00:42:20Z her127 $
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
/* advect_u1_3d()      : Solves the u1 advection equation            */
/* bdry_u1_3d()        : Sets normal and tangential u1 OBC's         */
/* extract_u1_3d()     : Gets u1 velocity from u1 transport          */
/*-------------------------------------------------------------------*/

void precalc_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void merge_thin_layers(int c, int zm1, double *vel, double *dz);
void reset_thin_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_thin_wtop(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     int *ctp, int vcs, double *wtop);
void set_dry_bdry(geometry_t *window, int nb, int *obc, int *ceni, double *vel);
void set_w_dry(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_wvel(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	      double *wvel, int mode);
void mom_balance(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_sdc_e1(geometry_t *window, window_t *windat, win_priv_t *wincon);
double gridze(geometry_t *geom, int e);
void set_sur_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void turbines(geometry_t *window, window_t *windat, win_priv_t *wincon);

	
/*-------------------------------------------------------------------*/
/* Window step part 1                                                */
/*-------------------------------------------------------------------*/
void mode3d_step_window_p1(master_t *master,   /* Master data        */
                           geometry_t *window, /* Window geometry    */
                           window_t *windat,   /* Window data        */
			   win_priv_t *wincon  /* Window constants   */
			   )
{

  double clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Fill 3D  velocities into the window data structures. This can   */
  /* be done inside this  window loop  (contrary to tracers) since   */
  /* velocity  is  only updated  after the 2D  mode and  hence all   */
  /* windows are filled with velocities at the current time-step.    */
  win_data_fill_3d(master, window, windat, master->nwindows);

  /* Get the cell centred east and north velocities                  */
  vel_cen(window, windat, wincon, windat->u1, windat->u2, 
	  windat->u, windat->v, NULL, NULL, 0);

  /*-----------------------------------------------------------------*/
  /* Calculate a heat and salt flux if required                      */
  calc_heatf(window, windat, wincon);
  calc_saltf(window, windat, wincon);

  init_sigma(window, windat, wincon);

#if !GLOB_BC
  /* Set the lateral boundary conditions for velocity at t for       */
  /* multiple windows. These velocities are required for the         */
  /* advective terms, and require transfers before ghost cells       */
  /* adjacent to auxiliary cells can be set.                         */
  if (master->nwindows > 1) {
    vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	      window->bpte1, window->bine1, wincon->slip);
  }
  /* Set the lateral boundary conditions for velocity at t-1.        */
  /* These velocities are required for stress tensors. Ghost cells   */
  /* are not set accurately since a lateral BC is not set on         */
  /* nu1 before the asselin() filtering; setting BCs here overcomes  */
  /* this.                                                           */
  if (!(master->compatible & V1283)) {
    vel2D_lbc(windat->u1b, window->nbpte1, window->nbe1,
	      window->bpte1, window->bine1, wincon->slip);
  }
#endif

  /*-----------------------------------------------------------------*/
  /* Set the vertical  mixing coefficients. This  is done  in this   */
  /* window loop so that transferred velocities from other windows   */
  /* can be used in the velocity shear term.                         */
  wincon->calc_closure(window, windat, wincon);
  bdry_closure(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Evalulate the sources and sinks of water                        */
  ss_water(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the cfl time-steps if required                              */
  if (!(wincon->cfl & NONE))
    calc_cfl(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the vertical grid spacings at edges                         */
  set_dz_at_u1(window, windat, wincon);
  debug_c(window, D_INIT, D_POST);

  /*-----------------------------------------------------------------*/
  /* Set the stress tensors and Smagorinsky diffusion.               */
  wincon->hor_mix->pre(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Check for fatal wind instabilities                              */
  if (check_unstable(window, windat, wincon, WIND)) return;

  /*-----------------------------------------------------------------*/
  /* Transfer grid spacings and mixing coefficients to the master.   */
  if (master->nwindows > 1)
    win_data_empty_3d(master, window, windat, MIXING);

  windat->wclk = (dp_clock() - clock);
}

/* END mode3d_step_window_p1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window step part 2                                                */
/*-------------------------------------------------------------------*/
void mode3d_step_window_p2(master_t *master,   /* Master data        */
                           geometry_t *window, /* Window geometry    */
                           window_t *windat,   /* Window data        */
			   win_priv_t *wincon  /* Window constants   */
			   )
{

  double clock = dp_clock();
  windat->dt = windat->dtf + windat->dtb;

  /*-----------------------------------------------------------------*/
  /* Fill the window with u1 data (dzu1, dzu2, Vz, Kz) from the      */
  /* master.                                                         */
  win_data_refill_3d(master, window, windat, master->nwindows, MIXING);

  /*-----------------------------------------------------------------*/
  /* Update the 3D velocity                                          */
  vel_u1_update(window, windat, wincon);
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
void mode3d_post_window_p1(master_t *master,   /* Master data        */
                           geometry_t *window, /* Window geometry    */
                           window_t *windat,   /* Window data        */
			   win_priv_t *wincon  /* Window constants   */
			   )
{

  double clock = dp_clock();

  /* Extract the velocities from the updated solution                */
  extract_u1_3d(window, windat, wincon);

  /* Invoke time filtering and step forward                          */
  leapfrog_update_3d(window, windat, wincon);

  /* Add the velocity increments to 3D velocity if required          */
  do_vel_relax(window, windat, wincon, VEL3D);
  /*do_vel_increment_3d(window);*/

  /* Calculate the 3D velocity alert diagnostics                     */
  alerts_w(window, VEL3D);

  /* Check for fatal instabilties                                   */
  if (check_unstable(window, windat, wincon, VEL3D)) return;

  /* Adjust  the  velocities. 3D  fluxes calculated  below must  use */
  /* the adjusted velocities.                                        */
  velocity_adjust(window, windat, wincon);

  /* Calculate the fluxes.                                           */
  set_flux_3d(window, windat, wincon, VEL3D);

  /* Set the velocities at newly wetted cells                        */
  set_new_cells_u1(window, windat, wincon);

  /* Set the lateral boundary conditions for velocity.               */
#if !GLOB_BC
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
            window->bpte1, window->bine1, wincon->slip);
#endif

  /* Transfer adjusted 3D velocities and 3D fluxes to the master     */
  win_data_empty_3d(master, window, windat, VELOCITY|CFL);

  windat->wclk += (dp_clock() - clock);

}

/* END mode3d_post_window_p1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window post step part 2                                           */
/* Set the 3D fluxes through the cell edges, calculate the           */
/* vertical velocity and transfer this to the master.                */
/*-------------------------------------------------------------------*/
void mode3d_post_window_p2(master_t *master,   /* Master data        */
                           geometry_t *window, /* Window geometry    */
                           window_t *windat,   /* Window data        */
			   win_priv_t *wincon  /* Window constants   */
			   )
{

  double clock = dp_clock();

  /* Fill  the window  with adjusted  velocities and 3D  fluxes from */
  /* the master.  These are  required in  the vertical  velocity (3D */
  /* fluxes) routine  and tracer  (velocity and 3D fluxes)  routine. */
  win_data_refill_3d(master, window, windat, master->nwindows, VELOCITY);

  /* Set up the cell centered surface vectors and dz arrays          */
  set_dz(window, windat, wincon);

  /* Calculate the vertical velocity                                 */
  vel_w_update(window, windat, wincon);

  /* Calculate the velocity components edges                         */
  /*vel_components_3d(window, windat, wincon);*/

  /* Get the tangential velocity to the edge                         */
  vel_tan_3d(window, windat, wincon);

  /* Calculate the vertical velocity alert diagnostics               */
  alerts_w(window, WVEL);

  /* Transfer vertical velocity to the master                        */
  if (master->nwindows > 1)
    win_data_empty_3d(master, window, windat, WVEL);

  /* Transfer alert data if required                                 */
  if (!(master->alertf & NONE))
    master_alert_fill(master, window, windat);

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
  if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);
  if (master->regf & RS_TSSET) timeseries_init_w(master, window);
  if (master->regf & RS_PSSSET) sourcesink_reinit(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Reset the means on the master if required.                      */
  /* This is now done within a scheduled function
     reset_means_m_o(master);*/
  /* Reset means if FFSL is computed using mean velocities           */
  if (master->means & TRANSPORT) {
    if (master->nstep % (int)master->tratio == 0) {
      reset_means_m(master);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Read the elevation at the start of the timestep for tiled       */
  /* coupling.                                                       */
  if (master->obcf & DF_TILE) 
    bdry_tiled_eta(geom, master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the custom tracer routines on the master                     */
  TIMING_SET;
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
  /* Do the custom velocity routines on the master                   */
  TIMING_SET;
  bdry_eval_u1_m(geom, master);
  if (master->regf & RS_OBCSET) bdry_reconfigure(master, window);
  TIMING_DUMP(2,"  bdry3d_u1_u2 ");

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
/* fluxes through the cell edges and calculate vertical velocity.    */
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

  /*-----------------------------------------------------------------*/
  /* Adjust velocities.                                              */
  dp_vel3d_post_p1();

  /*-----------------------------------------------------------------*/
  /* Set the  lateral boundary conditions for velocity. This is done */
  /* at this stage for velocity so that the velocities calculated in */
  /* velocity_adjust()  above  are used  for the tangential  lateral */
  /* boundary ghost cells.                                           */
  /* NOTE : this should be also done for velocity after any velocity */
  /* initialisation.                                                 */
#if GLOB_BC
  vel2D_lbc(master->u1, geom->nbpte1, geom->nbe1,
            geom->bpte1, geom->bine1, master->slip);
#endif

  /* Set 3D fluxes                                                   */
  dp_vel3d_post_p2();

#if TR_CK
  check_transfers(geom, window, windat, wincon, nwindows, VEL3D|FLUX);
#endif

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
  int e, ee, es;                /* Sparse coordinates / counters     */
  double ramp;                  /* Ramp value                        */

  /*-----------------------------------------------------------------*/
  /* Do the custom routines on the master                            */
  bdry_eval_tr_m(geom, master);
  bdry_eval_u1_m(geom, master);

  /*-----------------------------------------------------------------*/
  /* Loop over all windows to calculate the vertical mixing          */
  /* coefficients and set the vertical grid spacings.                */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];
    ramp = (wincon[n]->rampf & WIND) ? windat[n]->rampval : 1.0;

    /*---------------------------------------------------------------*/
    /* Set the vertical grid spacings edges                          */
    win_data_fill_3d(master, window[n], windat[n], nwindows);

    /*---------------------------------------------------------------*/
    /* Calculate a heatflux if required                              */
    calc_heatf(window[n], windat[n], wincon[n]);
    set_dz_at_u1(window[n], windat[n], wincon[n]);
    init_sigma(window[n], windat[n], wincon[n]);
    wincon[n]->hor_mix->pre(window[n], windat[n], wincon[n]);
    wincon[n]->hor_mix->setup(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Get the e1 vertical pressure integrals for the 2D mode.       */
    /* Set pointers and initialise.                                  */
    windat[n]->nu1 = wincon[n]->w9;
    memset(wincon[n]->u1inter, 0, window[n]->szeS * sizeof(double));
    memset(windat[n]->u1bot, 0, window[n]->szeS * sizeof(double));
    /* Set the maps to be self-mapping                               */
    set_map_e1(window[n]);
    /* Get the cells to process at the edge                          */
    cells2process_e1(window[n], windat[n], wincon[n]);
    /* Transfer any custom data                                      */
    bdry_transfer_u1(master, window[n], windat[n]);
    /* Get the vertical pressure integrals.                          */
    pressure_u1(window[n], windat[n], wincon[n]);
    /* Allow for tidal energy extraction                             */
    if (wincon[n]->nturb > 0) turbines(window[n], windat[n], wincon[n]);

    /* Get the momentum flux OUT of the water due to the wind        */
    for (ee = 1; ee <= wincon[n]->vcs; ee++) {
      e = wincon[n]->s1[ee];    /* 3D cell to process                */
      es = window[n]->m2de[e];  /* 2D cell corresponding to e        */

      wincon[n]->u1inter[es] += (ramp * windat[n]->wind1[es] /        
				 wincon[n]->topdensu1[es]);

      wincon[n]->u1inter[es] /= windat[n]->depth_e1[es];
      wincon[n]->topdensu1[es] *= wincon[n]->g;
      wincon[n]->densavu1[es] *= window[n]->h2au1[es];
    }
    memset(windat[n]->u1, 0, window[n]->sze * sizeof(double));
    /* Get the open boundary conditions                              */
    bdry_u1_3d(window[n], windat[n], wincon[n]);

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
  int e, ee, es;                /* Sparse coodinate / counter        */
  double *depth;                /* Depth of the water column         */
  int dwn = 0;                  /* Debugging window                  */
  int dwe = 0;                  /* Global debugging edge             */
  int dwew = geom->fm[dwe].ec;  /* Local debugging edge              */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  depth = wincon->d5;

  /*-----------------------------------------------------------------*/
  /* Get the cells to process at the edges                           */
  cells2process_e1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Precalculate variables required for the u1 calculation          */
  precalc_u1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across OBC's                    */
  set_map_e1(window);

  /*-----------------------------------------------------------------*/
  /* Read the 3D state variables at the start of the timestep for    */
  /* tiled coupling.                                                 */
  if (master->obcf & DF_TILE) bdry_tiled_3d(window, windat, wincon);

  /* Debug output                                                    */
  if (dwe > 0) {
    if(window->nwindows == 1)
      printf("start %d %f %f\n",master->nstep,windat->days,windat->nu1[dwe]);
    if(window->nwindows > 1 && window->wn==dwn)
      printf("start %d %f %f\n",master->nstep,windat->days,windat->nu1[dwew]);
  }

  /*-----------------------------------------------------------------*/
  /* Do the horizontal advection                                     */
  if (wincon->momsc & RINGLER) {
    if (nonlin_coriolis_3d(window, windat, wincon)) return;
  } else {
    if (advect_u1_3d(window, windat, wincon)) return;
  }
  if (wincon->tendf && wincon->tendency)
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_ADV], windat->u1_adv, windat->u2_adv);
  debug_c(window, D_U, D_ADVECT);

  if (dwe > 0) {
    if(window->nwindows == 1)
      printf("advect %e\n",windat->nu1[dwe]);
    if(window->nwindows > 1 && window->wn==dwn)
      printf("advect %e\n",windat->nu1[dwew]);
  }

  /*-----------------------------------------------------------------*/
  /* Add the dispersion terms to u1inter. Note: u1adv is multiplied  */
  /* by the 2D sub-timestep in the 2D mode in order to integrate     */
  /* over the 2D mode, hence u1adv needs to be divided by dt here to */
  /* complete the mean.                                              */
  if (wincon->nonlinear && !(wincon->u1_f & ADVECT)) {
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[e];
      /* Add the dispersion terms to u1inter                         */
      wincon->u1inter[es] =
        (wincon->u1inter[es] + wincon->u1adv[es]) / windat->dtf;
    }
  }

  /* Initialise the dispersion terms                                 */
  memset(wincon->u1adv, 0, window->szeS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Do horizontal diffusion (with metrics included)                 */
  wincon->hor_mix->setup(window, windat, wincon);
  wincon->hor_mix->u1(window, windat, wincon);

  if (wincon->tendf && wincon->tendency)
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_HDF], windat->u1_hdif, windat->u2_hdif);
  debug_c(window, D_U, D_HDIFF);

  /* Reset thin layers                                               */
  reset_thin_u1(window, windat, wincon); 

  if (dwe > 0) {
    if(window->nwindows == 1)
      printf("hdiff %e\n",windat->nu1[dwe]);
    if(window->nwindows > 1 && window->wn==dwn)
      printf("hdiff %e\n",windat->nu1[dwew]);
  }

  /*-----------------------------------------------------------------*/
  /* Include the pressure gradient term                              */
  pressure_u1(window, windat, wincon);
  debug_c(window, D_U, D_PRESSURE);

  if (dwe > 0) {
    if(window->nwindows == 1)
      printf("press %e\n",windat->nu1[dwe]);
    if(window->nwindows > 1 && window->wn==dwn)
      printf("press %e\n",windat->nu1[dwew]);
  }

  /*-----------------------------------------------------------------*/
  /* Include the Coriolis term                                       */
  if (!(wincon->momsc & RINGLER)) {
    if (!(wincon->u1_f & CORIOLIS))
      coriolis_u1(window, windat, wincon);
  }

  /*-----------------------------------------------------------------*/
  /* Include Stokes forcing                                          */
  if (!(wincon->u1_f & CORIOLIS))
    stokes_u1(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion                                       */
  if (vdiff_u1(window, windat, wincon)) return;
  if (wincon->tendf && wincon->tendency)
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_VDF], windat->u1_vdif, windat->u2_vdif);
  debug_c(window, D_U, D_VZ);

  if (dwe > 0) {
    if(window->nwindows == 1)
      printf("vdiff %16.14e\n",windat->nu1[dwe]);
    if(window->nwindows > 1 && window->wn==dwn)
      printf("vdiff %16.14e\n",windat->nu1[dwew]);
  }

  /*-----------------------------------------------------------------*/
  /* Sum the tendencies                                              */
  if (!(wincon->compatible & V6257)) {
    for (ee = 1; ee <= wincon->vc; ee++) {
      int n;
      e = wincon->s1[ee];
      /* Note: T_ADV is not added, as this velocity already is       */
      /* updated with this tendedcy in momentum advection to account */
      /* for sub-stepping.                                           */
      for (n = 1; n < TEND3D; n++)
	windat->nu1[e] += wincon->tend3d[n][e];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allow for tidal energy extraction                               */
  if (wincon->nturb > 0) turbines(window, windat, wincon);

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
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    /* Divide internal terms for 2d part by depth                    */
    wincon->u1inter[es] /= max(depth[es], wincon->hmin);

    /* Other calculations to save time                               */
    wincon->topdensu1[es] *= wincon->g;
    wincon->densavu1[es] *= window->h2au1[es];
  }

  /*-----------------------------------------------------------------*/
  /* Calculate u1 boundary values                                    */
  bdry_u1_3d(window, windat, wincon);

  debug_c(window, D_U, D_POST);
}

/* END vel_u1_update()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness at the edges.                   */
/*-------------------------------------------------------------------*/
void set_dz_at_u1(geometry_t *window,  /* Window geometry            */
                  window_t *windat,    /* Window data                */
                  win_priv_t *wincon   /* Window constants           */
  )
{
  int ee;                       /* Counter                           */
  int e, e3;                    /* Sparse coordinate                 */
  int zm1;                      /* Sparse cell below cell c3         */
  int es, eb;                   /* Sparse coordinate of the bottom   */
  int c1, c2;                   /* Cell centres straddling e         */
  double top;                   /* Top edge of a cell                */
  double bot;                   /* Bottom edge of a cell             */
  double *maxeta;               /* Maximum elevation at the edge     */
  double gridz;
  int n, nc = window->nz + 1;

  /*-----------------------------------------------------------------*/
  /* Set the pointers                                                */
  maxeta = wincon->d5;

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->szeS;
    memcpy(window->sur_e1, window->w2_e1, size * sizeof(int));

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom. This is performed over all wet + auxiliary cells      */
    /* since dzu1[] is required at the auxiliary cells in the        */
    /* advection routine advect_u1_3d().                             */
    memset(windat->dzu1, 0, window->sze * sizeof(double));
    memset(maxeta, 0, window->szeS * sizeof(double));
    for (ee = 1; ee <= window->b2_e1; ee++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      e = e3 = window->sur_e1[ee];  /* 3D surface coordinate         */
      es = window->m2de[e];         /* 2D cell corresponding to e    */
      eb = window->bot_e1[ee];      /* 3D bottom coordinate          */
      c1 = window->e2c[es][0];
      c2 = window->e2c[es][1];

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = windat->topz[es];
      while (e != eb) {
        bot = gridze(window, e);
        windat->dzu1[e] = top - bot;
        top = bot;
        e = window->zm1e[e];
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu1[eb] = top - window->botzu1[es];
    }
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Find the 3D sparse coordinate corresponding to the free       */
    /* surface.                                                      */
    set_sur_u1(window, windat, wincon);

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells from the surface to the  */
    /* bottom. This is performed over all wet + auxiliary cells      */
    /* since dzu1[] is required at the auxiliary cells in the        */
    /* advection routine advect_u1_3d().                             */
    memset(windat->dzu1, 0, window->sze * sizeof(double));
    wincon->ncdry_e1 = 0;
    for (ee = 1; ee <= window->b2_e1; ee++) {

      /* Get the 2D and 3D sparse coordinates of the surface         */
      e = e3 = window->sur_e1[ee];  /* 3D surface coordinate         */
      es = window->m2de[e];         /* 2D cell corresponding to e    */
      eb = window->bot_e1[ee];      /* 3D bottom coordinate          */

      /* Set the cell thickness from the surface to the layer above  */
      /* the bottom.                                                 */
      top = maxeta[es];
      if (nc) n = 0;

      while (e != eb) {
        bot = gridze(window, e);
        windat->dzu1[e] = top - bot;
        top = bot;
        e = window->zm1e[e];

	if (nc) {
	  n++;
	  if (n > nc) {
	    printf("stuck in set_dz_at_u1: wn=%d ee=%d e=%d(%d) eb=%d c=(%d %d) %f %f\n",
		   window->wn, ee, es, window->wse[e], eb, 
		   window->s2i[window->e2ijk[e]], window->s2i[window->e2ijk[e]], 
		   window->u1x[es], window->u1y[es]);
	    exit(0);
	  }
	}
      }

      /* Set the cell thickness at the bottom                        */
      windat->dzu1[eb] = top - window->botzu1[es];

      /* Save the locations of dry water columns                     */
      if (windat->dzu1[eb] == 0.0) {
	wincon->ncdry_e1++;
	wincon->cdry_e1[wincon->ncdry_e1] = es;
      }

      /* Set the cell thickness above the old surface equal to zero. */
      /* (Note : a no-gradient condition may be required for higher  */
      /* order vertical advection schemes.                           */
      e = e3;
      top = windat->dzu1[e];
      while (e != es) {
        e = window->zp1e[e];
        windat->dzu1[e] = top;
      }
    }
  }
}

/* END set_dz_at_u1()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the location of the surface.                       */
/*-------------------------------------------------------------------*/
void set_sur_u1(geometry_t *window,  /* Window geometry              */
		window_t *windat,    /* Window data                  */
		win_priv_t *wincon   /* Window constants             */
  )
{
  int ee;                       /* Counter                           */
  int e, e3;                    /* Sparse coordinate                 */
  int es, eb;                   /* Sparse coordinate of the bottom   */
  int c1, c2;                   /* Cell centres straddling e         */
  int zm1;                      /* Sparse cell below cell e3         */
  double top;                   /* Top edge of a cell                */
  double *maxeta;               /* Maximum elevation at edge      e  */
  double g1,g2;
  /*-----------------------------------------------------------------*/
  /* Set the pointers                                                */
  maxeta = wincon->d5;

  memset(windat->sur_e1, 0, window->szeS * sizeof(int));
  for (ee = 1; ee <= window->b2_e1; ee++) {

    /* Get the 2D and 3D sparse coordinates                          */
    es = e = window->w2_e1[ee]; /* 2D coordinate                     */
    c1 = window->e2c[es][0];
    c2 = window->e2c[es][1];

    /* Get the new sparse location of the surface. The surface at    */
    /* the edge is defined as the higher of the cell centered        */
    /* elevations bordering the edge. For elevations below the       */
    /* bottom this is set to the sediment coordinate.                */
    top = maxeta[es] = max(windat->eta[c1], windat->eta[c2]);
    top = maxeta[es] = max(top, window->botzu1[es]);
    zm1 = window->zm1e[e];

    /* Loop down the water column                                    */
    while (e != zm1 && gridze(window, e) > top) {
      e = zm1;
      zm1 = window->zm1e[e];
    }
    window->sur_e1[ee] = e;

    /* windat->sur_e1 is the k level of the surface and is           */
    /* tranferred across windows. Used to get the surface layer      */
    /* in auxiliary cells (note window->sur_e1 is only defined       */
    /* over cells b2_e1).                                            */
    windat->sur_e1[es] = max(window->s2k[window->e2c[e][0]],
			     window->s2k[window->e2c[e][1]]);
  }
}

/* END set_surf_u1()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set a value in dry cells above the surface             */
/*-------------------------------------------------------------------*/
void set_surf_cond(geometry_t *window,  /* Window geometry           */
		   double *a,
		   double val,
		   int *vec,
		   int nvec)
{
  int *csur = window->wincon->i7;
  int ee, e, es;
  double sval;

  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    if ((e = csur[es])) {
      sval = (val < 0) ? a[e] : val;    
      while (e != window->zp1e[e]) {
	e = window->zp1e[e];
	a[e] = sval;
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
  int e, ee, es, ks;
  int *csur = wincon->i7;
  int *cells, *acells, vc, *kth, nkth, *bpt, *bin, nbpt;

  /* Set pointers                                                    */
  vc = window->b2_e1;
  cells = window->sur_e1;
  acells = windat->sur_e1;
  nkth = wincon->nkth_e1;
  kth = wincon->kth_e1;
  bpt = window->bpte1S;
  bin = window->bine1S;
  nbpt = window->nbpte1S;

  /* Initialise and get wet window cells (OBC cells included)        */
  memset(csur, 0, window->szeS * sizeof(int));
  for (ee = 1; ee <= vc; ee++) {
    e = cells[ee];
    es = window->m2de[e];
    csur[es] = e;
  }

  /* Get auxiliary cells                                             */
  if (window->nwindows > 1) {
    for (es = 1; es < window->szeS; es++) {
      ks = acells[es];
      if (ks) {
	e = get_local_sur(window, es, ks);
	csur[es] = e;
      }
    }
  }

  /* Get cells reset at thin layers                                  */
  if (wincon->thin_merge) {
    for (ee = 1; ee <= nkth; ee++) {
      e = window->zm1e[kth[ee]];
      es = window->m2de[e];
      csur[es] = e;
    }
  }

  /* Get ghost cells                                                 */
  for (ee = 1; ee <= nbpt; ee++) {
    es = bin[ee];
    e = bpt[ee];
    csur[e] = csur[es];
  }
}

/* END set_surf_cells()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to exclude boundary cells from the cells to process so as */
/* to linearize the momentum equations on the boundary or within N   */
/* cells from the boundary.                                          */
/*-------------------------------------------------------------------*/
void linear_bdry_cell(geometry_t *window,  /* Window geometry        */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int *ctp,            /* Cells to process       */
		      int nctp,            /* Number of surface ctp  */
		      int *mask,           /* linear cell mask       */
		      int *bottom          /* Bottom cell vector     */
		      )
{  int e, ee, eb;
  int vc, vcs;
  int *cells = wincon->s4;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memcpy(cells, ctp, window->szeS* sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the surface cells to process vector                         */
  vcs = 1;
  for (ee = 1; ee <= nctp; ee++) {
    e = ctp[ee];
    if(!mask[e]) {
      cells[vcs] = e;
      vcs++;
    }
  }
  vcs--;
  wincon->aclS = vcs;

  /*-----------------------------------------------------------------*/
  /* Return if 2d cells only are processed                           */
  if (bottom == NULL) {
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Set the sub-surface cells to process vector                     */
  vc = vcs + 1;
  for (ee = 1; ee <= nctp; ee++) {
    e = ctp[ee];
    eb = bottom[ee];
    if(!mask[e]) {
      e = window->zm1e[e];
      while (e <= eb) {
	cells[vc] = e;
	vc++;
	e = window->zm1e[e];
      }
    }
  }

  vc--;
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
  int ee, bn;                   /* Counter                           */
  int e, es;                    /* Sparse coordinate                 */
  int zm1;                      /* Sparse cell below cell c          */
  int eb;                       /* Bottom sparse coordinate          */
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

  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];
    eb = bottom[ee];
    es = window->m2de[e];
    zm1 = window->zm1e[e];
    if (e < eb && windat->dzu1[zm1] && windat->dzu1[e] < wincon->hmin) {
      cells[ee] = zm1;
      /* Save the sparse location of the thin layer                  */
      wincon->kth_e1[wincon->nkth_e1] = e;
      wincon->nkth_e1++;
      u1o[window->m2de[e]] = windat->u1[zm1];	
      merge_thin_layers(e, zm1, windat->u1, windat->dzu1);
    } else
      cells[ee] = e;
    windat->sur_e1[es] = 0;
  }

  /* Set thin layers for open boundary cells. These are used in the  */
  /* non-linear terms.                                               */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (ee = 1; ee <= open[bn]->no3_e1; ee++) {
	e = open[bn]->obc_e1[ee];
	es = window->m2de[e];
	zm1 = window->zm1e[e];
	ks = windat->sur_e1[es];
	if (windat->dzu1[zm1] && windat->dzu1[e] < wincon->hmin) {
	  merge_thin_layers(e, zm1, windat->u1, windat->dzu1);
	}
	windat->sur_e1[es] = 0;
      }
    }
  }

  /* Set thin layers for auxiliary cells. These are used in the      */
  /* non-linear terms.                                               */
  if (window->nwindows > 1) {
    for (es = 1; es < window->szeS; es++) {
      ks = windat->sur_e1[es];
      if (ks) {
	e = get_local_sur(window, es, ks);
	zm1 = window->zm1e[e];
	if (windat->dzu1[zm1] && windat->dzu1[e] < wincon->hmin) {
	  wincon->kth_e1[wincon->nkth_e1] = e;
	  wincon->nkth_e1++;
	  u1o[es] = windat->u1[zm1];
	  merge_thin_layers(e, zm1, windat->u1, windat->dzu1);
	}
      }
    }
  }
  wincon->nkth_e1--;

  /* Set a no-gradient of the merged velocity above the surface      */
  for (ee = 1; ee <= wincon->nkth_e1; ee++) {
    e = wincon->kth_e1[ee];
    zm1 = window->zm1e[e];
    for(es = window->m2de[e]; es != e; es = window->zm1e[es]) {
      windat->u1[es] = windat->u1[zm1];
    }
  }

  /* Fill the cells to process with sub-surface cells. If no thin    */
  /* cells were encountered make a direct copy of s1 to save time.   */
  if (wincon->nkth_e1 == 0) {
    memcpy(wincon->s3, wincon->s1, window->sze * sizeof(int));
    wincon->ncl = wincon->vc;
  } else {
    wincon->ncl = wincon->vcs + 1;
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = window->zm1e[cells[ee]];
      eb = bottom[ee];
      while (e <= eb) {
        cells[wincon->ncl] = e;
        wincon->ncl++;
        e = window->zm1e[e];
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
  int ee;                       /* Counter                           */
  int e, zm1;                   /* Sparse coordinate                 */
  double *u1o = wincon->d7;     /* Original u1 at zm1                */

  /* Set a uniform velocity over thin layers if required             */
  if (wincon->thin_merge) {
    for (ee = 1; ee <= wincon->nkth_e1; ee++) {
      e = wincon->kth_e1[ee];
      zm1 = window->zm1e[e];
      windat->nu1[e] = windat->nu1[zm1];
      windat->dzu1[zm1] -= windat->dzu1[e];
      /* Reset u1 to the un-merged velocity                          */
      windat->u1[e] = (windat->dzu1[e]) ? windat->u1[zm1] + 
	               windat->dzu1[zm1] * 
 	              (windat->u1[zm1] - u1o[window->m2de[e]]) / 
 	               windat->dzu1[e] : 0.0;
      windat->u1[zm1] = u1o[window->m2de[e]];
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
  int cc, ee, e, es;            /* Counter                           */
  int c, cs, c1, c2;            /* Sparse coordinates                */
  double eta_l;                 /* Surface elevation at i-1          */
  double eta_r;                 /* Surface elevation at i+1          */
  double eta_b;                 /* Surface elevation at j-1          */
  double eta_f;                 /* Surface elevation at j+1          */

  memcpy(wtop, windat->wtop, window->szcS * sizeof(double));

  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cs = window->m2d[c];
    wtop[cs] = windat->detadt[cs];
    for (ee = 1; ee <= window->npe[cs]; ee++) {
      e = window->c2e[ee][cs];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      /* Calculate the surface vertical velocity                     */
      eta_l = windat->eta[c1];
      eta_r = windat->eta[c2];
      windat->wtop[cs] += (windat->u1[es] * (eta_r - eta_l) / 
			   (window->h2au1[es] * (double)window->npe[cs]));
    }
  }
}
/* END reset_thin_wtop()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the sparse coordinates of cells to process at the  */
/* edges.                                                            */
/*-------------------------------------------------------------------*/
void cells2process_e1(geometry_t *window, /* Window geometry         */
                      window_t *windat,   /* Window data             */
                      win_priv_t *wincon  /* Window constants        */
		      )
{
  int ee;                       /* Counter                           */
  int e, e2;                    /* Sparse coordinate                 */
  int c1, c2, c1s, c2s;         /* Sparse cell at (i-1)              */
  int es, eb;                   /* Sparse coordinate of the bottom   */
  double top;                   /* Top edge of a cell                */
  double bot;                   /* Bottom edge of a cell             */
  double *maxeta;               /* Maximum elevation at the edge     */
  int vc, vcs;                  /* Counter                           */
  int c;
  int nc = 0;

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->szeS;
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
      for (ee = 1; ee <= window->b2_e1; ee++) {
	es = window->w2_e1[ee];
	e = window->sm_e1[es];
	c1 = window->e2c[e][0];
	if (e) {
	  at_ele[ee] = window->sur_e1[ee];
	  window->sur_e1[ee] = e;
	  maxeta[es] = min(maxeta[es], window->gridz[window->zp1[c1]]);
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
    for (ee = 1; ee <= vcs; ee++) {
      e = window->sur_e1[ee];   /* Surface coordinate                */
      es = window->w2_e1[ee];   /* 2D coordinate                     */
      top = maxeta[es];         /* Highest elevation at the edge     */
      bot = DRY_FRAC * wincon->hmin + window->botzu1[es];
      if (top > bot) {
        wincon->s1[vc] = e;
        eb = bottom[vc] = window->bot_e1[ee];
	/* Save cells which are one cell deep                        */
	if(e == eb) {
	  wincon->ncbot_e1++;
	  wincon->cbot_e1[wincon->ncbot_e1] = e;
	}
        vc++;
      }
    }
    wincon->vcs = vc - 1;

    /*---------------------------------------------------------------*/
    /* Loop from the layer below the surface to the bottom and get   */
    /* the cells to process. The are arranged so that cells 1 to     */
    /* wincon->vcs contain surface cells, cells wincon->vcs+1 to     */
    /* wincon->vca contain cells where one side of the cell edge is  */
    /* dry and cells wincon->vca+1 to wincon->vc contain cells where */
    /* both sides of the cell edge is wet.                           */
    /* First get the cells where one side of the edge is dry. The    */
    /* last cell that falls into this category is temporarily saved  */
    /* in partial[] for use in the next loop.                        */
    for (ee = 1; ee <= wincon->vcs; ee++) {

      e = e2 = wincon->s1[ee];  /* 3D surface cell                   */
      es = window->m2de[e];     /* 2D cell corresponding to c        */
      eb = bottom[ee];          /* 3D bottom cell                    */
      c1 = window->e2c[es][0];
      c2 = window->e2c[es][1];

      at_ele[ee] = e;
      top = min(windat->eta[c1], windat->eta[c2]);

      /* Cells where one edge is completely dry are added to s1.     */
      /* The cell containing the surface is stored in at_ele, i.e    */
      /* if eta[c1s] > eta[c2s] at_ele[ee] contains the higher cell  */
      /* location, s1[ee]. If eta[c2s] > eta[c1s] at_ele[ee] contains*/
      /* the lower cell location which will be the same as           */
      /* partial[ee].                                                */
      c1 = window->e2c[e][0];
      c1s = window->e2c[es][0];
      bot = max(window->gridz[c1], window->botzu1[es]);
      while (e < eb && (bot = max(window->gridz[c1], 
				  window->botzu1[es])) > top) {
        if (e != e2) {
          wincon->s1[vc] = e;
          vc++;
        }
        if (bot > windat->eta[c1s])
          at_ele[ee] = window->zm1e[at_ele[ee]];
        e = window->zm1e[e];
	c1 = window->e2c[e][0];
      }

      /* Get the cell location where the surface resides. This is    */
      /* current value of e after the loop above (i.e when gridz     */
      /* becomes less than eta). Add this cell to s1 if not already  */
      /* accounted for. hi_ele = partial and lo_ele = c2s if eta[c1s]*/
      /* > eta[c2s], and hi_ele = c2s and lo_ele = partial if        */
      /* eta[c2s] > eta[c1s].                                        */
      if (e <= eb) {
        if (e != e2) {
          wincon->s1[vc] = e;
          vc++;
        }
	c1 = window->e2c[e][0];
	c2 = window->e2c[e][1];
	c2s = window->m2d[c2];
        partial[ee] = e;
	hi_ele[ee] = c1;
        lo_ele[ee] = c2;
	/* Swap if eta[c2s] > eta[c1s]                                 */
        if (top < windat->eta[c2s]) {
          lo_ele[ee] = c1;
          hi_ele[ee] = c2;
        }
      }
    }
    wincon->vca = vc - 1;

    /* Next get the cells where both sides of the edge is wet.       */
    for (ee = 1; ee <= wincon->vcs; ee++) {
      int n = 0;
      e = window->zm1e[partial[ee]];
      eb = bottom[ee];
      while (e <= eb) {
        wincon->s1[vc] = e;
        vc++;
        e = window->zm1e[e];
	if (nc) {
	  n++;
	  if (n > nc) {
	    printf("stuck in cells2process_e1: wn=%d ee=%d e=%d(%d) eb=%d c=(%d %d)\n",
		   window->wn, ee, e, window->wse[e], eb, 
		   window->s2i[window->e2ijk[e]], window->s2i[window->e2ijk[e]]);
	    exit(0);
	  }
	}
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
int vdiff_u1(geometry_t *window,  /* Window geometry                 */
	     window_t *windat,    /* Window data                     */
	     win_priv_t *wincon   /* Window constants                */
  )
{
  int e, ee, eoe;          /* Sparse coordinates / counters          */
  int es, eb;              /* Surface / bottom sparse coordinate     */
  int c1, c2;              /* Cell location to the west of c         */
  int *at_ele;             /* Array of fully wet cell coordinates    */
  int *partial;            /* Coordinate of cell containing eta[c]   */
  int *bottom;             /* Bottom coordinate for cells to process */
  double *Vzav;            /* Viscosity at the edge                  */
  double *f_top;           /* Flux out of the water at the surface   */
  double *f_bot;           /* Flux into the water at the bottom      */
  double *u2au1;           /* Tangential velocity at edge            */
  double Cdu1;             /* Bottom drag at the edge                */
  double botdz;            /* Thickness of the bottom layer          */
  double val;              /* Bottom current speed                   */
  double *w;               /* Vertical velocity                      */
  int n, yp1, xmyp1;
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
  /* Set the vertical viscosity at the edge. Where one side of the   */
  /* water column is higher, set the viscosity at the cell edge      */
  /* equal to the viscosity of the higher water column. If the       */
  /* current water column, cs, is the higher then no action is taken */
  /* (c == cs), if the water column to the west of edge cs is higher */
  /* then the viscosity at c is set to the viscosity at m1. This is  */
  /* only performed to the layer above that containing the free      */
  /* surface in column cs; below this the mean is taken.             */
  /* Set the viscosity on completely dry edges                       */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];
    es = at_ele[ee];
    c1 = window->e2c[e][0];
    Vzav[e] = windat->Vz[c1];
    while (e < es) {
      c2 = window->e2c[e][1];
      Vzav[e] = windat->Vz[c2];
      e = window->zm1e[e];
    }
  }

  /* Set the viscosity for cells where the edge is partially dry     */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = partial[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    Vzav[e] = 0.5 * (windat->Vz[c1] + windat->Vz[c2]);
  }

  /* Set the viscosity for cells where the edge completely wet       */ 
  for (ee = wincon->vca + 1; ee <= wincon->vc; ee++) {
    e = wincon->s1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    Vzav[e] = 0.5 * (windat->Vz[c1] + windat->Vz[c2]);
  }

  /*-----------------------------------------------------------------*/
  /* Set the vertical velocity                                       */
  memset(w, 0, window->sze * sizeof(double));
  if (wincon->momsc & WIMPLICIT) {
    int zm1;
    for (ee = 1; ee <= wincon->vc; ee++) {
      e = wincon->s1[ee];         /* 3D cell to process              */
      zm1 = window->zm1e[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      w[zm1] = 0.5 * (windat->w[c1] + windat->w[c2]);
    }
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];       /* 3D cell to process                */
      es = window->m2de[e];     /* 2D cell corresponding to e        */
      eb = bottom[ee];          /* 3D bottom coordinate              */
      zm1 = window->zm1e[eb];
      c1 = window->e2c[es][0];
      c2 = window->e2c[es][1];
      /* Bottom vertical velocity                                    */
      w[zm1] = 0.5 * (windat->wbot[c2] + windat->wbot[c1]);
      /* Surface vertical velocity                                   */
      w[e] = 0.5 * (windat->wtop[c2] + windat->wtop[c1]);
    }
    set_wvel(window, windat, wincon, w, 0);
  }

  /*-----------------------------------------------------------------*/
  /* Set the momentum flux out of the water at the top (due to the   */
  /* wind) and flux into the water at the bottom (due to friction).  */
  memset(f_top, 0, window->szeS * sizeof(double));
  memset(f_bot, 0, window->szeS * sizeof(double));
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];         /* 3D cell to process                */
    es = window->m2de[e];       /* 2D cell corresponding to e        */
    eb = bottom[ee];            /* 3D bottom coordinate              */
    c1 = window->e2c[es][0];
    c2 = window->e2c[es][1];

    /* Surface flux                                                  */
    f_top[es] = -ramp * windat->wind1[es] / wincon->topdensu1[es];

    /* Add surface stress to 2d forcing term.                        */
    if (!(wincon->u1av_f & VDIFF))
      wincon->u1inter[es] -= f_top[es];

    /* Bottom flux */
    u2au1[eb] = 0.0;
    for (n = 1; n <= window->nee[es]; n++) {
      eoe = window->eSe[n][eb];
      if (!eoe) continue;
      u2au1[eb] += window->wAe[n][eb] * windat->u1b[eoe];
    }
    Cdu1 = 0.5 * (wincon->Cd[c1] + wincon->Cd[c2]);
    val = sqrt(windat->u1b[eb] * windat->u1b[eb] + u2au1[eb] * u2au1[eb]);
    val = Cdu1 * max(wincon->uf, val);
    botdz = max(wincon->hmin, windat->dzu1[eb] * wincon->mdx[es]);
    /* Quadratic bottom friction, truncated to ensure stability      */
    if (val > botdz / windat->dt)
      val = botdz / windat->dt;
    f_bot[es] = -val * windat->u1b[eb];
  }
  if (wincon->vorticity & TENDENCY)
    memcpy(windat->rv_bsc, f_bot, window->szeS * sizeof(double));
  if (wincon->numbers & BOTSTRESS)
    memcpy(windat->tau_be1, f_bot, window->szeS * sizeof(double));
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
    if (wincon->compatible & V6257) {
      if (implicit_vdiff(window, windat, wincon, windat->u1b, windat->nu1, Vzav,
			 windat->dzu1, f_bot, f_top, wincon->s1, wincon->i5,
		       wincon->vcs, wincon->Ds, wincon->Ds, w)) return(1);
    } else {
      if (implicit_vdiff(window, windat, wincon, windat->u1b, wincon->tend3d[T_VDF], Vzav,
			 windat->dzu1, f_bot, f_top, wincon->s1, wincon->i5,
			 wincon->vcs, wincon->Ds, wincon->Ds, w)) return(1);
    }
  }
  return(0);
}

/* END vdiff_u1()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to compute the Coriolis tendency                          */
/*-------------------------------------------------------------------*/
void coriolis_u1(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int e, ee, es;                /* Sparse coordinates / counters     */
  double *midx;                 /* Total depth at edge               */
  double *tend;

  /* Set pointers */
  midx = wincon->d4;            /* Set in pressure_u1()              */
  if (wincon->compatible & V6257)
    tend = windat->nu1;
  else
    tend = wincon->tend3d[T_COR];

  /*-----------------------------------------------------------------*/
  /* Add the Coriolis term to the u1 velocity                        */
  /* SIGMA : multiply by total depth at u1                           */
  for (ee = 1; ee <= wincon->vc; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    tend[e] += windat->dt * wincon->u1c5[es] * windat->u2[e] * midx[es];
  }

  if (wincon->tendf && wincon->tendency) {
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_COR], windat->u1_cor, windat->u2_cor);
  }
}

/* END coriolis_u1()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to compute the Stokes Coriolis and vortex forcing         */
/* tendency.                                                         */
/*-------------------------------------------------------------------*/
void stokes_u1(geometry_t *window,  /* Window geometry               */
	       window_t *windat,    /* Window data                   */
	       win_priv_t *wincon   /* Window constants              */
  )
{
  int e, ee, es;                /* Sparse coordinates / counters     */
  double *midx;                 /* Total depth at edge               */
  int v, vv, vs, v1, v2;
  double d2, iarea;
  double sdc, vf;
  double *tend;

  if(!(wincon->waves & STOKES_DRIFT)) return;

  /* Set pointers */
  midx = wincon->d4;            /* Set in pressure_u1()              */
  if (wincon->compatible & V6257)
    tend = windat->nu1;
  else
    tend = wincon->tend3d[T_STK];

  /*-----------------------------------------------------------------*/
  /* Stokes drift                                                    */
  get_sdc_e1(window, windat, wincon);
  for (ee = 1; ee <= wincon->vc; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    sdc = wincon->u1c5[es] * wincon->w1[e] * midx[es];
    tend[e] += windat->dt * sdc;
    wincon->u1inter[es] += sdc * windat->dzu1[e];
  }

  /*-----------------------------------------------------------------*/
  /* Relative vorticity at vertices                                  */
  memset(windat->rvor, 0, window->szv * sizeof(double));
  for (vv = 1; vv <= window->b3_e2; vv++) {
    v = window->w3_e2[vv];
    vs = window->m2dv[v];
    iarea = 1.0 / window->dualarea[vs];
    for (ee = 1; ee <= window->nve[vs]; ee++) {
      e = window->v2e[v][ee];
      d2 = window->h2au1[e] * windat->u1[e];
      windat->rvor[v] += window->eSv[ee][v] * iarea * d2;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Vortex force                                                    */
  for (ee = 1; ee <= wincon->vc; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    v1 = window->e2v[e][0];
    v2 = window->e2v[e][1];
    sdc = wincon->w1[e] * midx[es] * 0.5 * (windat->nrvor[v1] + windat->nrvor[v2]);
    tend[e] += windat->dt * sdc;
    wincon->u1inter[es] += sdc * windat->dzu1[e];
  }
  if (wincon->tendf)
    vel_cen(window, windat, wincon, wincon->tend3d[T_STK], NULL, windat->u1_sto, windat->u2_sto, 
	    NULL, NULL, 0);
}

/* END stokes_u1()                                                   */
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
  int e, ee;              /* Local sparse coordinate / counter       */
  int es;                 /* 2D local sparse coordinate              */
  int cl, cs;             /* Partially wet partial cell coordinate   */
  int ch;                 /* Fully wet partial cell coordinate       */
  int c, c1, c2, c1s, c2s; /* Cell centre counters                   */
  int zm1;                /* Cell below cell c                       */
  int zp1, xmzp1;         /* Cell above and left-above cell c        */
  int xm1, xm1s;          /* Cell location to the west of c          */
  int *partial;           /* Coordinate of cell containing eta[c]    */
  int *at_ele;            /* Array of fully wet cell coordinates     */
  int *hi_ele;            /* Array of fully wet cell coordinates     */
  int *lo_ele;            /* Array of partially wet cell coordinates */
  double *dhdiv;          /* Vertically integrated depth (hlower)    */
  /* double *dhdiu; *//* Vertically integrated depth (hupper)        */
  double *dd, *dc;        /* Vertically integrated pressure gradient */
  double *dinter;         /* Pressure gradient for the 2D mode       */
  double *u1inter;        /* Vertical integrals for 2d mode          */
  double *hidens;         /* Density corresponding to coordinate ch  */
  double *dav;            /* Mean density at the cell edge           */
  double *dzu1;           /* Cell theickness ( = windat->dzu1)       */
  double *hupper;         /* Partial cell upper thickness work array */
  double *hlower;         /* Partial cell lower thickness work array */
  double *dens;           /* Density ( = windat->dens)               */
  double *dzsum;          /* Sum of dz for dry / partial layers      */
  double *mask;           /* Mask out all fully wet or dry edges     */
  double top;             /* Top cell thickness                      */
  double bot;             /* Bottom cell thickness                   */
  double *midx;           /* Total depth at the cell edge            */
  double d1;              /* Dummy                                   */
  double tol = 1e-10;     /* Tolerence value                         */
  double *corr;           /* Baroclinic correction for sigma         */
  double *ddbc;           /* Baroclinic contribution to u1 velocity  */
  double *ddbt;           /* Barotropic contribution to u1 velocity  */
  double *tend1, *tend2;  /* Barotropic / baroclinic tendeny         */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  hidens = wincon->w1;
  hupper = wincon->w2;
  hlower = wincon->w3;
  dav = wincon->w4;
  dzsum = wincon->w5;
  mask = wincon->w6;
  corr = wincon->w8;

  dd = wincon->d1;
  dc = wincon->w7;
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
    ddbc = dc;
  if (wincon->u1_f & PRESS_BT)
    ddbt = wincon->d7;
  else
    ddbt = dd;
  if (wincon->u1av_f & PRESS_BC)
    dinter = wincon->d7;
  else
    dinter = u1inter;

  if (wincon->compatible & V6257) {
    tend1 = windat->nu1;
    tend2 = windat->nu1;
  } else {
    tend1 = wincon->tend3d[T_BTP];
    tend2 = wincon->tend3d[T_BCP];
  }

  /* For the 2D mode when baroclinic pressure is omitted, set        */
  /* surface and average density to the surface density, and return. */
  if (wincon->mode2d && wincon->u1av_f & PRESS_BC) {
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      wincon->densavu1[e] = wincon->topdensu1[e] = 0.5 * (windat->dens[c1] + windat->dens[c2]);
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->topdensu1, 0, window->szeS * sizeof(double));
  memset(wincon->densavu1, 0, window->szeS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Linear eta case; calculate surface pressure gradient due to eta */
  /* slope and atmospheric pressure gradient.                        */
  if (!wincon->nonlinear) {
    double *surgrad = wincon->d1;

    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];

      midx[es] = wincon->mdx[es];
      wincon->topdensu1[es] = 0.5 * (windat->dens[c1] + windat->dens[c2]);
      surgrad[es] = wincon->g * wincon->topdensu1[es] *
        (windat->eta[c1s] - windat->eta[c2s]) +
        (windat->patm[c1s] - windat->patm[c2s]);
      dc[es] = 0.0;
    }
    for (ee = 1; ee <= wincon->vc; ee++) {
      e = wincon->s1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];

      /* Density gradient integral                                   */
      dav[e] = 0.5 * (windat->dens[c1] + windat->dens[c2]);
      dc[es] +=
        wincon->g * (windat->dens[c1] -
                     windat->dens[c2]) * windat->dzu1[e];

      /* Add term to new u1 value                                    */
      /*
      windat->nu1[e] +=
        windat->dt * wincon->u1c6[es] * (surgrad[es] + dc[es]) / dav[e];
      */
      tend1[e] += windat->dt * wincon->u1c6[es] * surgrad[es] / dav[e];
      tend2[e] += windat->dt * wincon->u1c6[es] * dc[es] / dav[e];

      /* Integrate internal term for 2d part                         */
      wincon->u1inter[es] +=
        wincon->u1c6[es] * dc[es] * windat->dzu1[e] / dav[e];
      wincon->densavu1[es] += dav[e] * windat->dzu1[e];
    }
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];

      top = max(windat->eta[c1s], windat->eta[c2s]);
      bot = window->botzu1[es];
      wincon->densavu1[es] /= (top - bot);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Non-linear eta case                                             */
  else {

    /*---------------------------------------------------------------*/
    /* Precondition work arrays for the pressure gradient integral.  */
    /* These are established for the surface cells (i.e. cells that  */
    /* have one edge completely or partially dry) as follows:        */
    /*                                                               */
    /* work array   | completely dry edge  | partially dry edge      */
    /* -------------|----------------------|----------------------   */
    /* hidens       | rho at wet cell      | rho at full wet cell    */
    /* hupper       | cell thickness, dzu1 | upper thickness, hu     */
    /* hlower       | cell thickness, dzu1 | lower thickness, hl     */
    /* dav          | rho at wet cell      | 0.5(rho(c)+rho(xm1))    */
    /* dzsum        | rho(c) or -rho(xm1)  | rho(c) or -rho(xm1)     */
    /* mask         | zero                 | one                     */
    /* Note; c refers to the cell centre in the upwind direction     */
    /* from the edge, and xm1 refers to the cell centre in the       */
    /* downwind direction.                                           */

    /* Initialise the work arrays                                    */
    memcpy(hupper, dzu1, window->sze * sizeof(double));
    memcpy(hlower, dzu1, window->sze * sizeof(double));
    memset(mask, 0, window->szm * sizeof(double));
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      c1 = window->e2c[e][0];
      hidens[e] = dens[c1];
      dav[e] = dens[c1];
      dzsum[e] = dens[c1];
    }

    /* Set the work arrays for cells where the edge is completely    */
    /* dry. Where one side of the water column is higher, set the    */
    /* density at the cell edge equal to the density of the higher   */
    /* water column. If the current water column, cs, is the higher  */
    /* then no action is taken (c == cs), if the water column to the */
    /* west of edge cs is higher then the density at c is set to the */
    /* density at xm1. This is only performed to the layer above     */
    /* that containing the free surface in column cs.                */
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];
      es = at_ele[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      midx[window->m2de[es]] = 0.5 * (wincon->Ds[window->m2d[c1]] + 
				      wincon->Ds[window->m2d[c2]]);
      while (e < es) {
	c1 = window->e2c[e][0];
	c2 = window->e2c[e][1];
        hidens[e] = dav[e] = dens[c2];
	dzsum[e] = -dens[c2];
        e = window->zm1e[e];
      }
    }
    if (!wincon->sigma)
      memcpy(midx, wincon->one, window->szeS * sizeof(double));

    /* Set the work arrays for cells where the edge is partially dry */
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = partial[ee];         /* 3D sparse of partially dry cell    */
      ch = hi_ele[ee];         /* 3D location of higher water column */
      cl = lo_ele[ee];         /* 3D location of lower water column  */
      cs = window->m2d[cl];    /* 2D loaction of lower water column  */
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];

      bot = max(window->gridz[c1], window->botzu1[window->m2de[e]]);
      hlower[e] = windat->eta[cs] - bot;
      hupper[e] = windat->dzu1[e] - hlower[e];
      hidens[e] = dens[ch];
      dav[e] = 0.5 * (dens[c1] + dens[c2]);
      mask[e] = 1.0;
      dzsum[e] = (window->e2c[at_ele[ee]][0] == cl) ? -dens[c2] : dens[c1];
      /* Special treatment for hlower when the cell is completely    */
      /* dry (i.e. when eta <= botzu1).                              */
      if (windat->eta[c1s] <= window->botzu1[window->m2de[e]]) { 
	hlower[e] = windat->dzu1[e];
	mask[e] = 0.0;
      }
    }

    /* Set the sub-surface explicit map if required                  */
    if (window->sm_e1) {
      for (ee = 1; ee <= wincon->vcs; ee++) {
	es = window->m2de[wincon->s1[ee]];
	e = window->sm_e1[es];
	if (e) {
	  hlower[e] = windat->dzu1[e];
	  c1 = window->e2c[window->m2d[hi_ele[ee]]][0];
	  c2 = window->e2c[window->m2d[lo_ele[ee]]][0];
	  hupper[e] = windat->eta[c1] - windat->eta[c2];
	  while (window->gridz[c1] > windat->eta[c1]) {
	    c1 = window->zm1[c1];
	  }
	  c1 = window->e2c[es][0];
	  c2 = window->e2c[es][0];
	  dzsum[e] = (at_ele[ee] == lo_ele[ee]) ? -dens[c2] : 
	    dens[c1];
	}
      }
    }

    /* Set the work arrays for cells where the edge is wet. Where    */
    /* water exists on both sides of the cell edge, set the density  */
    /* at the cell edge to the mean of that at cell c and at cell    */
    /* xm1.                                                          */
    for (ee = wincon->vca + 1; ee <= wincon->vc; ee++) {
      e = wincon->s1[ee];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      dav[e] = 0.5 * (dens[c1] + dens[c2]);
    }

    /*---------------------------------------------------------------*/
    /* SIGMA : Add the surface elevation gradient and apply the      */
    /* gradient correction. Note that sigma (gridz) takes on values  */
    /* from 0 to -1.                                                 */
    memset(corr, 0, window->sze * sizeof(double));
    memset(dd, 0, window->szeS * sizeof(double));
    memset(dc, 0, window->szeS * sizeof(double));
    if (wincon->sigma) {

      /* Subtract the mean density                                   */
      for (c = 1; c <= window->enon; c++)
        dens[c] -= wincon->rdens[c];

      for (ee = 1; ee <= wincon->vcs; ee++) {
        e = wincon->s1[ee];
        zm1 = window->zm1e[e];
        es = window->m2de[e];
	c1s = window->e2c[es][0];
	c2s = window->e2c[es][1];

        /* Surface pressure gradient                                 */
        ddbt[es] = wincon->g * (windat->eta[c1s] - windat->eta[c2s]) +
          (windat->patm[c1s] - windat->patm[c2s]);

	/* Update velocity                                           */
	/*
	windat->nu1[e] +=
	  midx[es] * windat->dt * wincon->u1c6[es] * dd[es];
	*/
	tend1[e] += midx[es] * windat->dt * wincon->u1c6[es] * dd[es];
        e = zm1;
        zm1 = window->zm1e[e];

        while (e != zm1) {
	  int zp1, zp2, zm1, zm2;
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
          zp1 = window->zp1[c1];
          zp2 = window->zp1[c2];
          zm1 = window->zm1[c1];
          zm2 = window->zm1[c2];

          /* Density gradient correction                             */
          corr[e] = -wincon->g * midx[es] * window->cellz[c1] *
            (wincon->Ds[c1s] - wincon->Ds[c2s]) *
            0.5 * (dens[zp1] + dens[zp2] - dens[zm1] - dens[zm2]);

          /* Update velocity                                         */
	  /*
          windat->nu1[e] +=
            midx[es] * windat->dt * wincon->u1c6[es] * dd[es];
	  */
          tend1[e] += midx[es] * windat->dt * wincon->u1c6[es] * dd[es];

          e = zm1;
          zm1 = window->zm1e[e];
        }
      }
    }

    /*---------------------------------------------------------------*/
    /* Initialise vertical integral arrays                           */
    memset(dd, 0, window->szeS * sizeof(double));
    memset(dc, 0, window->szeS * sizeof(double));
    memset(u1inter, 0, window->szeS * sizeof(double));
    memset(dhdiv, 0, window->szeS * sizeof(double));

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the surface cells (completely or   */
    /* partially dry edges). The completely dry edges are the        */
    /* 'surface' part, and corresponds to the deta/de1 term in the   */
    /* 2D mode. For the sigma case this loop is always skipped.      */
    for (ee = 1; ee <= wincon->vca; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->e2c[es][0];
      c2s = window->e2c[es][1];

      /* Get the pressure / surface gradients                        */
      d1 = midx[es] * midx[es] * (dens[c1] - dens[c2]) * mask[e];
      ddbc[es] += wincon->g * d1 * hlower[e];
      dinter[es] += wincon->g * mask[e] * d1 * hlower[e];
      ddbt[es] += wincon->g * dzsum[e] * hupper[e];

      /* Upadte the u1 velocity                                      */
      /*
      windat->nu1[e] += windat->dt * wincon->u1c6[es] * (dd[es] + dc[es]) / dav[e];
      */
      tend1[e] += windat->dt * wincon->u1c6[es] * dd[es] / dav[e];
      tend2[e] += windat->dt * wincon->u1c6[es] * dc[es] / dav[e];

      /* Integrate surface and average density for 2D mode           */
      wincon->topdensu1[es] += hidens[e] * hupper[e];
      wincon->densavu1[es] += dav[e] * hlower[e];
      wincon->u1inter[es] +=
        wincon->u1c6[es] * u1inter[es] * hlower[e] / dav[e];

      dhdiv[es] += hlower[e];
    }
    /* mark1 */

    /*---------------------------------------------------------------*/
    /* Integrate the pressure for the wet cell edges                 */
    for (ee = wincon->vca + 1; ee <= wincon->vc; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];

      /* Get the pressure / surface gradients                        */
      d1 =
        wincon->g * midx[es] * midx[es] * (dens[c1] - dens[c2]) * dzu1[e];
      ddbc[es] += (d1 + corr[e]);
      dinter[es] += (d1 + corr[e]);

      /* Update the u1 velocity                                      */
      /*
      windat->nu1[e] += windat->dt * wincon->u1c6[es] * (dd[es] + dc[es]) / dav[e];
      */
      tend1[e] += windat->dt * wincon->u1c6[es] * dd[es] / dav[e];
      tend2[e] += windat->dt * wincon->u1c6[es] * dc[es] / dav[e];

      /* Integrate surface and average density for 2D mode           */
      wincon->u1inter[es] +=
        wincon->u1c6[es] * u1inter[es] * dzu1[e] / dav[e];

      wincon->densavu1[es] += dav[e] * dzu1[e];
      dhdiv[es] += dzu1[e];
    }

    /*---------------------------------------------------------------*/
    /* Divide mean densities by depth                                */
    for (ee = 1; ee <= wincon->vcs; ee++) {
      e = wincon->s1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->e2c[es][0];
      c2s = window->e2c[es][1];

      wincon->densavu1[es] /= dhdiv[es];
      d1 = fabs(windat->topz[c1s] - windat->topz[c2s]);
      wincon->topdensu1[es] =
        (d1 > tol) ? wincon->topdensu1[es] / d1 : dav[e];
    }

    /* Reset the density by adding the mean density                  */
    if (wincon->sigma) {
      for (c = 1; c <= window->enon; c++)
        dens[c] += wincon->rdens[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the pressure gradient tendencies if required.               */
  if (wincon->tendency && windat->u1_btp) {
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_BTP], windat->u1_btp, windat->u2_btp);
  }
  if (wincon->tendency && windat->u1_bcp) {
    get_tendv(window, wincon->s1, wincon->vc,
	      wincon->tend3d[T_BCP], windat->u1_bcp, windat->u2_bcp);
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
  int e, ee, es, eoe;
  int n;
  double *depth;

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  windat->nu1 = wincon->w9;
  depth = wincon->d5;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(wincon->u1inter, 0, window->szeS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Copy the current velocities into the update array               */
  memcpy(windat->nu1, windat->u1b, window->sze * sizeof(double));
  for (n = 0; n < TEND3D; n++)
    memset(wincon->tend3d[n], 0, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Multiply the velocity by the depth at the backward timestep for */
  /* SIGMA calculations.                                             */
  if (wincon->sigma) {
    for (ee = 1; ee <= wincon->vc; ee++) {
      e = wincon->s1[ee];
      windat->nu1[e] *= mdxbs(window, windat, wincon, e);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the depth at the edge.                                      */
  /* u1inter  is divided  by the depth  in vel_u1_update()  (for the */
  /* pressure and surface stress terms) but the dispersion terms are */
  /* not  to  be divided  by depth, so  multiply  by  depth here  to */
  /* compensate. Note: the work array wincon->d5 contains the height */
  /* of the surface at the edge at this stage, set in                */
  /* set_dz_at_u1().                                                 */
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    /* Get the water depth at the edge                               */
    depth[es] -= window->botzu1[es];

    wincon->u1adv[es] *= max(depth[es], wincon->hmin);
  }

  /* Initialise the tendency array of required                       */
  if (wincon->tendf && wincon->tendency)
    memcpy(wincon->tendency, windat->nu1, window->sze * sizeof(double));
}

/* END precalc_u1()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute velocity components                            */
/*-------------------------------------------------------------------*/
void vel_components_3d(geometry_t *window,  /* Window geometry       */
		       window_t *windat,    /* Window data           */
		       win_priv_t *wincon   /* Window constants      */
		       )
{

  /* Get the tangential velocity to the edge                         */
  vel_tan_3d(window, windat, wincon);

  /* Get the cell centred east and north velocities                  */
  vel_cen(window, windat, wincon, windat->u1, windat->u2, 
	  windat->u, windat->v, NULL, NULL, 0);
}

/* END vel_components_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute tangential velocity at edges                   */
/*-------------------------------------------------------------------*/
void vel_tan_3d(geometry_t *window,    /* Window geometry            */
		window_t *windat,      /* Window data                */
		win_priv_t *wincon     /* Window constants           */
  )
{
  int e, ee, es, eoe, n;
  double fs = 1.0;

  /*-----------------------------------------------------------------*/
  /* Get the u2 (tangential) velocity at the edge.                   */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    windat->u2[e] = 0.0;
    for (n = 1; n <= window->nee[es]; n++) {
      eoe = window->eSe[n][e];
      fs = window->h1au1[window->m2de[eoe]] / window->h2au1[es];
      if (!eoe) continue;
      windat->u2[e] += fs * window->wAe[n][e] * windat->u1[eoe];
    }
  }
}

/* END vel_tan_3d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute velocity components at the cell centre         */
/*-------------------------------------------------------------------*/
void vel_cen(geometry_t *window,  /* Window geometry                 */
	     window_t *windat,    /* Window data                     */
	     win_priv_t *wincon,  /* Window constants                */
	     double *u1,          /* Nor edge velocity to centre     */
	     double *u2,          /* Tan edge velocity to centre     */
	     double *u,           /* Cell centered u velocity        */
	     double *v,           /* Cell centered v velocity        */
	     double *mag,         /* Cell centred velocity magnitude */
	     double *dir,         /* Cell centred velocity direction */
	     int mode             /* 2D or 3D                        */
	     )
{
  int e, ee, es;
  int c, cs, cc;
  int *vec, nvec, sz;
  double a, nu, nv, *ut;
  geometry_t *geom=master->geom;

  /* Set pointers                                                    */
  if (mode) {
    vec = window->w2_t;
    nvec = window->a2_t;
    sz = window->szcS;
  } else {
    vec = window->w3_t;
    nvec = window->a3_t;
    sz = window->szc;
  }

  memset(u, 0, sz * sizeof(double));
  memset(v, 0, sz * sizeof(double));

  /* Compute the tangential component if not supplied                */
  if (u2 == NULL) {
    int *vee, nvee, nveg, n, eoe;
    ut = wincon->w2;
    memset(ut, 0, window->sze * sizeof(double));
    if (mode) {
      vee = window->w2_e1;
      nvee = window->a2_e1;
      nveg = window->n2_e1;
    } else {
      vee = window->w3_e1;
      nvee = window->a3_e1;
      nveg = window->n3_e1;
    }
    for (ee = 1; ee <= nvee; ee++) {
      e = vee[ee];
      es = window->m2de[e];
      ut[e] = 0.0;
      for (n = 1; n <= window->nee[es]; n++) {
	eoe = window->eSe[n][e];
	if (!eoe) continue;
	ut[e] += window->wAe[n][e] * u1[eoe];
      }
    }
    if (wincon->slip == 1) {
      for (ee = nvee + 1; ee <= nveg; ee++) {
	e = vee[ee];
	es = window->m2de[e];
	ut[e] = 0.0;
	for (n = 1; n <= window->nee[es]; n++) {
	  eoe = window->eSe[n][e];
	  if (!eoe) continue;
	  ut[e] += 2.0 * window->wAe[n][e] * u1[eoe];
	}
      }
    }
  } else {
    ut = u2;
  }

  /* Average u and v around all edges of the cell                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cs = window->m2d[c];
    u[c] = 0.0;
    v[c] = 0.0;
    nu = nv = 0.0;

    for (ee = 1; ee <= window->npe[cs]; ee++) {
      e = window->c2e[ee][c];
      es = window->m2de[e];
      a = 0.5 * window->h1au1[es] * window->h2au1[es];
      /* Get the cell centered east and north velocity               */
      u[c] += a * (u1[e] * window->costhu1[es] + ut[e] * window->costhu2[es]);
      nu += a;
      v[c] += a * (u1[e] * window->sinthu1[es] + ut[e] * window->sinthu2[es]);
      nv += a;
    }
    u[c] = (nu) ? u[c] / nu : 0.0;
    v[c] = (nv) ? v[c] / nv : 0.0;

    if (mag != NULL)
      mag[c] = sqrt(u[c] * u[c] + v[c] * v[c]);
    if (dir != NULL)
      dir[c] =  fmod(atan2(u[c], v[c]) * 180.0 / M_PI + 180.0, 360.0);
  }
}

/* END vel_cen()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute velocity gradients at the cell centre          */
/*-------------------------------------------------------------------*/
void vel_grad(geometry_t *window,  /* Window geometry                */
	      window_t *windat,    /* Window data                    */
	      win_priv_t *wincon,  /* Window constants               */
	      double *u1,          /* Nor edge velocity to centre    */
	      double *u2,          /* Tan edge velocity to centre    */
	      double *ung,         /* Normal velocity gradient       */
	      double *utg,         /* Tangential velocity gradient   */
	      int mode             /* 2D or 3D, nor or tan           */
	      )
{
  int e, ee, es;
  int i, si, ii;
  int ne;
  int *vec, nvec, *m2d;
  double nu, nv, *ut;
  double *uv = wincon->w2;
  double *vv = wincon->w3;

  /* Get the (u,v) compoinents on vertices                           */
  if (mode & GRAD_NOR) {
    ne = window->npem;
    if (mode & GRAD_3D) {
      vec = window->w3_t;
      nvec = window->n3_t;
    }
    if (mode & GRAD_2D) {
      vec = window->w2_t;
      nvec = window->n2_t;
    }
    m2d = window->m2d;
    memset(uv, 0, window->szc * sizeof(double));
    memset(vv, 0, window->szc * sizeof(double));
  }
  if (mode & GRAD_TAN) {
    ne = window->nvem;
    if (mode & GRAD_3D) {
      vec = window->w3_e2;
      nvec = window->n3_e2;
    }
    if (mode & GRAD_2D) {
      vec = window->w2_e2;
      nvec = window->n2_e2;
    }
    m2d = window->m2dv;
    memset(uv, 0, window->szv * sizeof(double));
    memset(vv, 0, window->szv * sizeof(double));
  }

  /* Compute the tangential component if not supplied                */
  if (u2 == NULL) {
    int *vee, nvee, n, eoe;
    ut = wincon->w2;
    if (mode) {
      vee = window->w2_e1;
      nvee = window->a2_e1;
    } else {
      vee = window->w3_e1;
      nvee = window->a3_e1;
    }
    for (ee = 1; ee <= nvee; ee++) {
      e = vee[ee];
      es = window->m2de[e];
      ut[e] = 0.0;
      for (n = 1; n <= window->nee[es]; n++) {
	eoe = window->eSe[n][e];
	if (!eoe) continue;
	ut[e] += window->wAe[n][e] * u1[eoe];
      }
    }
  } else {
    ut = u2;
  }

  /* Get the east and north components at centres or vertices        */
  for (ii = 1; ii <= nvec; ii++) {
    i = vec[ii];
    si = m2d[i];
    uv[i] = 0.0;
    vv[i] = 0.0;
    nu = nv = 0.0;
    ne = (mode & GRAD_NOR) ? window->npe[si] : window->nve[si];
    for (ee = 1; ee <= ne; ee++) {
      e = (mode & GRAD_NOR) ? window->c2e[ee][i] : window->v2e[i][ee];
      es = window->m2de[e];
      /* Get the cell centered east and north velocity               */
      if (e) {
	/* Get the cell centered east and north velocity             */
	uv[i] += (u1[e] * window->costhu1[es] + ut[e] * window->costhu2[es]);
	nu += 1.0;
	vv[i] += (u1[e] * window->sinthu1[es] + ut[e] * window->sinthu2[es]);
	nv += 1.0;
      }
    }
    uv[i] = (nu) ? uv[i] / nu : 0.0;
    vv[i] = (nv) ? vv[i] / nv : 0.0;
  }

  /* Get the normal velocity gradient dun/dn, dut/dn. Rotate the     */
  /* east and north vectors normal and tangential, then get the      */
  /* differences.                                                    */
  if (mode & GRAD_NOR) {
    int c1, c2;
    for (ee = 1; ee <= window->a3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      ung[e] = ((uv[c1] * window->costhu1[es] + vv[c1] * window->sinthu1[es]) -
		(uv[c2] * window->costhu1[es] + vv[c2] * window->sinthu1[es])) /
	window->h2au1[es];

      utg[e] = ((vv[c1] * window->costhu1[es] - uv[c1] * window->sinthu1[es]) - 
		(vv[c2] * window->costhu1[es] - uv[c2] * window->sinthu1[es])) /
	window->h2au1[es];
    }
  }

  /* Get the tangential velocity gradient dun/dt, dut/dt             */
  if (mode & GRAD_TAN) {
    int v1, v2;
    for (ee = 1; ee <= window->a3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      v1 = window->e2v[e][0];
      v2 = window->e2v[e][1];
      ung[e] = (uv[v1] * window->costhu1[es] + vv[v1] * window->sinthu1[es] -
		uv[v2] * window->costhu1[es] + vv[v2] * window->sinthu1[es]) /
	window->h1au1[es];
      utg[e] = (vv[v1] * window->costhu1[es] - uv[v1] * window->sinthu1[es] - 
		vv[v2] * window->costhu1[es] - uv[v2] * window->sinthu1[es]) /
	window->h1au1[es];
    }
  }
}

/* END vel_grad()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to compute a scalar gradient at the cell centre           */
/*-------------------------------------------------------------------*/
void tra_grad(geometry_t *window,  /* Window geometry                */
	      window_t *windat,    /* Window data                    */
	      win_priv_t *wincon,  /* Window constants               */
	      double *tr,          /* tracer to get gradients of     */
	      double *ung,         /* Normal velocity gradient       */
	      double *utg,         /* Tangential velocity gradient   */
	      int mode             /* 2D or 3D, nor or tan           */
	      )
{
  int cc, c, cs, cn;
  int e, ee, es;
  int *vec, nvec, *m2d, ne;
  double grad, nu, nv;

  /* Get the (u,v) compoinents on vertices                           */
  ne = window->npem;
  if (mode & GRAD_3D) {
    vec = window->w3_t;
    nvec = window->n3_t;
    memset(ung, 0, window->szc * sizeof(double));
    memset(utg, 0, window->szc * sizeof(double));
  }
  if (mode & GRAD_2D) {
    vec = window->w2_t;
    nvec = window->n2_t;
    memset(ung, 0, window->szcS * sizeof(double));
    memset(utg, 0, window->szcS * sizeof(double));
  }
  m2d = window->m2d;

  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cs = m2d[c];
    ung[c] = 0.0;
    utg[c] = 0.0;
    nu = nv = 0.0;
    ne = window->npe[cs];
    for (ee = 1; ee <= ne; ee++) {
      cn = window->c2c[ee][c];
      es = window->m2de[window->c2e[ee][c]];
      grad = window->eSc[ee][c] * (tr[cn] - tr[c]) / window->h2au1[es];

      /* Get the cell centered east and north velocity               */
      if (e) {
	ung[c] += grad * window->costhu1[es];
	nu += ceil(fabs(window->costhu1[es]));
	utg[c] += grad * window->sinthu1[es];
	nv += ceil(fabs(window->sinthu1[es]));
      }
    }
    ung[c] = (nu) ? ung[c] / nu : 0.0;
    utg[c] = (nv) ? utg[c] / nv : 0.0;
  }
}

/* END tra_grad()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Get an edge vector component from cell centered east and north    */
/* components.                                                       */
/*-------------------------------------------------------------------*/
double vel_c2e(geometry_t *window, double *u, double *v, int e)
{
  int es = window->m2de[e];
  int c1 = window->e2c[e][0];
  int c2 = window->e2c[e][1];
  double u1, u2;

  if (window->wgst[c1]) c1 = window->wgst[c1];
  if (window->wgst[c2]) c2 = window->wgst[c2];

  u1 = u[c1] * window->costhu1[es] + v[c1] * window->sinthu1[es];
  u2 = u[c2] * window->costhu1[es] + v[c2] * window->sinthu1[es];
  return(0.5 * (u1 + u2));
}

/* END vel_c2e()                                                     */
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
  /* Set the normal u1 velocities in the boundary cell               */
  for (n = 0; n < window->nobc; n++) {

    set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
            open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
            open[n]->cyc_e1, windat->nu1, windat->u1, 
	    windat->u1b, open[n]->bcond_nor,
            windat->dtf, &open[n]->datau1, 
	    open[n]->transfer_u1, open[n]->relax_zone_nor, (U1BDRY|U1GEN));

    /* Set open boundary velocities above the free surface equal to  */
    /* zero.                                                         */
    set_dry_bdry(window, open[n]->no3_e1, open[n]->obc_e1, open[n]->ceni, windat->nu1);
  }

  /*-----------------------------------------------------------------*/
  /* Set the tangential u1 velocities in the boundary cell           */
  for (n = 0; n < window->nobc; n++) {
    set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
            open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
            open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
            windat->u1, windat->u1b, open[n]->bcond_tan, 
	    windat->dtf, &open[n]->datau2, open[n]->transfer_u2, 
	    open[n]->relax_zone_tan, (U2BDRY|U1GEN));            
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient over the unused cell for interior staggers.   */
  /* This is for output appearance only.                             */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & INFACE) {
      int e, ee, e1;
      for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	e = open[n]->obc_e1[ee];
	e1 = e2e(window, e, open[n]->outi[ee]);
	windat->nu1[e1] = windat->nu1[e];
      }
    }
  }

  /* Rescale for sigma active OBCs                                   */
  if (wincon->sigma) {
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->type & U1BDRY && open->bcond_nor & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1, mdxns,
			1, open->no3_e1, open->obc_e1);
      if (open->type & U1BDRY && open->bcond_tan & (FILEIN|CUSTOM))
	reset_sigma_OBC(window, windat, wincon, windat->nu1, mdxns,
			open->no3_e1 + 1, open->to3_e1, open->obc_e1);
    }
  }
}

/* END bdry_u1_3d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Updates velocity in the open boundaries for tiled 2-way coupling  */
/* (at the 3D time-step). Here the velocity data from an alternative */
/* tile should be updated in u1 at the start of the time-step,       */
/* rather than updating nu1 at the end of the time-step.             */
/* This must be done because the scheduled dump for data exchange    */
/* occurs after the model hd_step for u1, so nu1 data from an        */
/* alternative tile is not available for an update prior to          */
/* leapfrog_update.                                                  */
/*-------------------------------------------------------------------*/
void bdry_tiled_3d(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon  /* Window constants           */
  )
{
  int cc, c, ee, e, n, nn, tn;         /* Counters                   */

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];

    if (open->options & OP_TILED) {
      /* Normal boundary velocity                                    */
      if (open->bcond_nor & FILEIN) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = open->obc_e1[ee];
	  windat->u1[e] = windat->nu1[e] = open->transfer_u1[ee];
	  /* Updates the backward velocity with the Asslin time      */
	  /* filter.                                                 */
	  windat->u1b[e] = windat->u1b[e] + 0.5 * 0.1 *
	    (open->d3[ee] - 2.0 * windat->u1b[e] + windat->u1[e]);
	}
      }
      /* Tangential boundary velocity                                */
      if (open->bcond_tan & FILEIN) {
	for (ee = open->no3_e1 + 1; ee <= open->to3_e1; ee++) {
	  e = open->obc_e1[ee];
	  windat->u1[e] = windat->nu1[e] = open->transfer_u2[ee];
	  windat->u1b[e] = windat->u1b[e] + 0.5 * 0.1 *
	    (open->d3[ee] - 2.0 * windat->u1b[e] + windat->u1[e]);
	}
      }

      /* Tracers                                                     */
      for (nn = 0; nn < wincon->ntbdy; nn++) {
	tn = wincon->tbdy[nn];
	if (open->bcond_tra[tn] & FILEIN) {
	  for (cc = 1; cc <= open->no3_t; cc++) {
	    c = open->obc_t[cc];
	    windat->tr_wc[tn][c] = open->t_transfer[open->trm[tn]][cc];
	  }
	}
      }
      /* Normal depth averaged velocity                              */
      if (open->bcond_nor2d & VERTIN) {
	double *sum = wincon->d1;
	int e1, c1, c2, zm1;
	memset(sum, 0, window->szeS * sizeof(double));	  
	for (ee = 1; ee <= open->no2_e1; ee++) {
	  e = open->obc_e1[ee];
	  e1 = window->m2de[e];
	  zm1 = window->zm1e[e];
	  while(e != zm1) {
	    sum[e1] += windat->u1[e] * windat->dzu1[e] * wincon->mdx[e1];
	    e = zm1;
	    zm1 = window->zm1e[e];
	  }
	}
	for (ee = 1; ee <= open->no2_e1; ee++) {
	  e = open->obc_e1[ee];
	  c1 = window->m2d[window->e2c[e][0]];
	  c2 = window->m2d[window->e2c[e][1]];
	  windat->depth_e1[e] = (max(windat->eta[c1], windat->eta[c2]) - 
				 window->botzu1[e]);
	  windat->u1av[e] = windat->nu1av[e] = sum[e] / 
	    max(windat->depth_e1[e] * wincon->mdx[window->m2de[e]], wincon->hmin);
	}
      }

      /* Tangential depth averaged velocity                          */
      if (open->bcond_tan2d & VERTIN) {
	double *sum = wincon->d1;
	int e1, c1, c2, zm1;
	memset(sum, 0, window->szeS * sizeof(double));	  
	for (ee = open->no3_e1 + 1; ee <= open->to2_e1; ee++) {
	  e = open->obc_e1[ee];
	  e1 = window->m2de[e];
	  zm1 = window->zm1e[e];
	  while(e != zm1) {
	    sum[e1] += windat->u1[e] * windat->dzu1[e] * wincon->mdx[e1];
	    e = zm1;
	    zm1 = window->zm1e[e];
	  }
	}
	for (ee = open->no3_e1 + 1; ee <= open->to2_e1; ee++) {
	  e = open->obc_e1[ee];
	  c1 = window->m2d[window->e2c[e][0]];
	  c2 = window->m2d[window->e2c[e][1]];
	  windat->depth_e1[e] = (max(windat->eta[c1], windat->eta[c2]) - 
				 window->botzu1[e]);
	  windat->u1av[e] = windat->nu1av[e] = sum[e] / 
	    max(windat->depth_e1[e] * wincon->mdx[window->m2de[e]], wincon->hmin);
	}
      }
    }
  }
}

/* END bdry_tiled_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the u1 velocity from the u1 velocity transport */
/*-------------------------------------------------------------------*/
void extract_u1_3d(geometry_t *window, /* Window geometry            */
                   window_t *windat,   /* Window data                */
                   win_priv_t *wincon  /* Window constants           */
  )
{
  int e, ee, es, n;             /* Local sparse coordinate / counter */

  /* Extract the u1 velocity from velocity transport.                */
  /* Note : the water depths at the forward time steps have been     */
  /* calculated in the last 2D iteration in velocity_update2D()      */
  /* and can be used here.                                           */
  if (wincon->sigma) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      windat->nu1[e] /= wincon->Hn1[es];
    }
  }
}

/* END extract_u1_3d()                                               */
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
  int e, ee, n;                 /* Local sparse coordinate / counter */
  double aconst = 0.1;
  double *u1 = wincon->w1;

  memcpy(u1, windat->u1, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Apply the Asselin filter to remove the computational mode.      */
  /* u1 velocity.                                                    */
  /* For tiled coupling this is done in velocity_adjust(), so that   */
  /* the adjusted nu1 is used in the filtering.                      */
  if (!(master->obcf & DF_TILE)) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      windat->u1[e] = windat->u1[e] + 0.5 * aconst *
	(windat->u1b[e] - 2.0 * windat->u1[e] + windat->nu1[e]);
    }
  }

  /* If the open boundary condition for 2D velocity is FILEIN, then   */
  /* do not perform time filtering (restore boundary velocities to    */
  /* those before filtering is performed).                            */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->bcond_nor & (CUSTOM|FILEIN)) {
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
	windat->u1[e] = u1[e];
	/* For tiled coupling, save the backward velocity for time    */
	/* filtering at the start of the next time-step.              */
	if (open->options & OP_TILED) open->d3[ee] = windat->u1b[e];
      }
    }
    if (open->bcond_tan & (CUSTOM|FILEIN)) {
      for (ee = open->no3_e1+1; ee <= open->to3_e1; ee++) {
	e = open->obc_e1[ee];
	windat->u1[e] = u1[e];
	if (open->options & OP_TILED) open->d3[ee] = windat->u1b[e];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the velocities for the new time level                       */
  memcpy(wincon->w10, windat->u1b, window->sze * sizeof(double));
  memcpy(windat->u1b, windat->u1, window->sze * sizeof(double));
  memcpy(windat->u2b, windat->u2, window->sze * sizeof(double));
  memcpy(windat->u1, windat->nu1, window->sze * sizeof(double));
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
  int e, es, ee;

  for (ee = sb; ee <= eb; ee++) {
    e = obc[ee];
    es = window->m2de[e];
    vel[e] *= Hn(window, windat, wincon, es);
  }
}

/* END reset_sigma_OBC()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to pre-condition vertical velocity for the advection      */
/* computation. Note: wvel has been initialized prior to calling     */
/* this routine.                                                     */
/*-------------------------------------------------------------------*/
void set_wvel(geometry_t *window,           /* Window geometry       */
	      window_t *windat,             /* Window data           */
	      win_priv_t *wincon,           /* Window constants      */
	      double *wvel,                 /* Conditioned w         */
	      int mode                      /* Redundant             */
	      )
{
  int e, ee, es, eb, zm1, zp1;
  int c1, c2, c1s, c2s;
  double top;

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get edge centered vertical        */
  /* velocity. Note that w is cell centered, hence the centered      */
  /* bottom coordinate must be used. Auxiliary cells are not defined */
  /* for cell centers on western and southern boundaries, so the _t  */
  /* vectors cannot be used here and the cell center must be found   */
  /* by mapping downwards from bot_e1 until self-mapping is located. */
  memset(wvel, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->a2_e1; ee++) {

    /* Get the 2D and 3D sparse coordinates                          */
    es = e = window->w2_e1[ee]; /* 2D coordinate                     */
    eb = window->bot_e1[ee];    /* Bottom coordinate                 */
    c1 = window->e2c[e][0];
    c1s = window->m2d[c1];
    top = windat->eta[c1s];
    while (e <= eb && window->gridz[c1] > top) {
      wvel[e] = windat->wtop[c1s];
      e = window->zm1e[e];
      c1 = window->e2c[e][0];
    }

    /* Set the bottom boundary condition. Find the cell centered     */
    /* bottom coordinate and set w to wbot.                          */
    e = eb;                       /* Eedge centered bottom cell      */
    while (e != window->zm1e[e])
      e = window->zm1e[e];        /* Edge centered sediment location */
    eb = window->zp1e[e];
    wvel[eb] = windat->wbot[c1s];
    wvel[e] = wvel[eb];
  }

  /* Set wvel at the surface for ghost cells. If window partitions   */
  /* lie next to ghost cells then these are not captured in the loop */
  /* above but are required for edge centered vertical velocity.     */
  /* Note that wvel[cb] always = 0 and is not explicitly set.        */
  for (ee = window->a2_e1 + 1; ee <= window->n2_e1; ee++) {
    es = e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    c1s = window->m2d[c1];
    top = windat->eta[c1s];
    while (e != window->zm1e[e] && window->gridz[c1] > top) {
      wvel[e] = windat->wtop[c1s];
      e = window->zm1e[e];
      c1 = window->e2c[e][0];
    }
  }

  /* Get the edge average velocity (mean of cell edges either side  */
  /* of e). This is staggered so that the average at layer e is     */
  /* stored in the layer below. This allows values at every         */
  /* vertical edge to be represented, with the bottom edge stored   */
  /* in the sediment cell and sea surface edge in the surface cell. */
  /* Loop from the bottom upwards and set the staggered average.    */
  for (ee = window->b2_e1+1; ee <= window->a2_e1; ee++) {
    e = window->bot_e1[ee];    /* Bottom coordinate                 */
    es = window->sur_e1[ee];   /* Surface coordinate                */
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    zm1 = window->zm1e[e];
    wvel[zm1] = 0.5 * (windat->w[c1] + windat->w[c2]);  /* Bottom   */
    while (e != es && e != window->zp1e[e]) { /* Mid water          */
      zm1 = e;
      e = window->zp1e[e];
      wvel[zm1] = wvel[e];
    }
    c1s = window->m2d[window->e2c[es][0]];
    wvel[e] = windat->wtop[c1s];
  }
  /* Set the surface staggered average                               */
  /*
  for (ee = 1; ee <= wincon->vcs; ee++) {
    e = wincon->s1[ee];
    wvel[e] = wvel[e];
  }
  */
}

/* END set_wvel()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_tr_3d(geometry_t *window,  /* Window geometry       */
                       window_t *windat,    /* Window data           */
                       double trem, /* Time remaining in sub-step    */
                       double *vel, /* Velocity array                */
                       double *as,  /* Multi-dt velocity at time t   */
                       double *ae   /* Multi-dt velocity at time t+1 */
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
/* Routine to set velocities above the free surface equal to zero on */
/* open boundaries.                                                  */
/*-------------------------------------------------------------------*/
void set_dry_bdry(geometry_t *window,             /* Window geometry */
		  int nb, int *obc, int *ceni, double *vel)
{
  int e, ee, j, c1, c1s;
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;

  if (wincon->sigma)
    return;

  for (ee = 1; ee <= nb; ee++) {
    e = obc[ee];
    j = ceni[ee];
    c1 = window->e2c[e][j];
    c1s = window->m2d[c1];
    vel[e] = (windat->eta[c1s] > window->gridz[c1] * wincon->Ds[c1s]) ? 
      vel[e] : 0.0;
  }
}

/* END set_dry_bdry()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the 3D fluxes through the cell edges. This is      */
/* performed over all wet cells only; auxiliary cell 3D fluxes are   */
/* set through a transfer from the master. This is preferrable to    */
/* calculating the 3D fluxes over all wet + auxiliary cells in the   */
/* window and dispensing with the need for transfers of 3D fluxes,   */
/* because complications can arise in the cells to process vectors   */
/* for bottom auxiliary cells where the bathymetry is stepped at the */
/* corresponding wet cell. In this case the wet cell may be a ghost  */
/* cell and thus not included in the cells to process vector v3_e1,  */
/* hence no corresponding auxiliary cell exists in a3_e1 and 3D      */
/* fluxes will not be set at the bottom, even though there may be a  */
/* non-zero velocity associated with the bottom auxiliary cell. This */
/* upsets the vertical velocity calculation at the wet cell which    */
/* subsequently results in incorrect surface tracer calculations     */
/* (continuity is violated).                                         */
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
  int e, ee, es;                /* Sparse coodinate / counter        */
  double *u;
  int n;

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  if (mode & TRANSPORT) {
    u = windat->u1;
  } else {
    u = windat->u1;
  }
  
  memset(windat->u1flux3d, 0, window->sze * sizeof(double));

  /* If using flux-form advection scheme in transport mode, transfer */
  /* u1vm to u1flux3d and u2vm to u2flux3d.                          */

  if (wincon->tmode & SP_U1VM) {

    /* Set the velocities at auxillary cells to zero */

    for (ee = window->b3_e1 + 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      windat->u1[e] = 0.0;
      windat->u1vm[e] = 0.0;
    }

    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      windat->u1flux3d[e] = windat->u1vm[e];
    }
    return;

  }

  /*-----------------------------------------------------------------*/
  /* Fluxes through the edges                                        */
    for (ee = 1; ee <= window->b3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      windat->u1flux3d[e] = u[e] * windat->dzu1[e] * window->h1au1[es] *
	wincon->mdx[es];
    }

    /* Set the fluxes for cells one cell deep                        */
    for (ee = 1; ee <= wincon->ncbot_e1; ee++) {
      e = wincon->cbot_e1[ee];
      es = window->m2de[e];
      windat->u1flux3d[e] = windat->u1flux[es]/windat->dt;
    }
    /* Zero fluxes above the surface                                 */
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      es = window->sur_e1[ee];
      while (e < es) {
        windat->u1flux3d[e] = 0.0;
        e = window->zm1e[e];
      }
    }

    /* Get the mean u1 velocity if required                          */
    if (wincon->means & VEL3D && windat->u1m) {
      double t = windat->dtf;

      for (ee = 1; ee <= window->b3_e1; ee++) {
	e = window->w3_e1[ee];
	es = window->m2de[e];
	c = window->e2c[e][0];
	if (window->wgst[c]) c = window->e2c[e][1];
	cs = window->m2d[c];
	/*if (windat->u1flux3d[c])
	  Averaging when layers are wet can lead to time normanlization  */
	/* issues when free surface drops out of layers; use the     */
	/* no-gradient on u1 to get a surface mean in the top layer. */
	windat->ume[e] = (windat->ume[e] * windat->meanc[cs] + 
			  windat->u1[e] * t) / (windat->meanc[cs] + t);
      }
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	windat->u1m[c] = (windat->u1m[c] * windat->meanc[cs] + 
			  windat->u[c] * t) / (windat->meanc[cs] + t);
	windat->u2m[c] = (windat->u2m[c] * windat->meanc[cs] + 
			  windat->v[c] * t) / (windat->meanc[cs] + t);
      }
    }

    /* Get the mean u1 volume flux if required                       */
    if (wincon->means & VOLFLUX && windat->u1vm) {
      double t = windat->dtf;
      if (wincon->means & TRANSPORT) windat->meanc[0] += 1.0;
      for (ee = 1; ee <= window->b3_e1; ee++) {
	e = window->w3_e1[ee];
	es = window->m2de[e];
	cs = window->e2c[es][0];
	c=window->e2c[e][0];
	if (window->wgst[cs]) cs = window->e2c[es][1];
	windat->u1vm[e] = (windat->u1vm[e] * windat->meanc[cs] + 
			   windat->u1flux3d[e] * t) / (windat->meanc[cs] + t);
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
  int e, ee, es, eb;            /* Sparse coodinate / counter        */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *midx;                 /* SIGMA : depth at the edge         */
  double *midy;                 /* SIGMA : depth at the edge         */
  double *depth;                /* Water depth at the edge           */
  double *adjust;               /* Velocity adjustment               */
  int bn;                       /* OBC counter                       */
  open_bdrys_t **open = window->open; /* Window OBC structure        */

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the edges.                                    */
  /* Set pointers and initialise                                     */
  sum = wincon->d1;
  midx = wincon->d2;
  depth = wincon->d5;           /* Set in precalc_u1()               */
  adjust = wincon->d3;

  memset(sum, 0, window->szeS * sizeof(double));
  /* The  number of surface  u1 cells to process (vcs1) may decrease */
  /* if a cell dries  during the 2D mode - hence adjust[] may not be */
  /* set  but u1 is still updated (with an old value residing in d3. */
  /* Therefore, zero adjust[] before use here.                       */
  memset(adjust, 0, window->szeS * sizeof(double));

  /* Integrate the velocities through the water column. Note: normal */
  /* velocity open boundary sparse coordinates are not included in   */
  /* the cells to process vectors (whereas tangential boundary       */
  /* velocities are) and must be processed separately.               */
  /* Wet cells to process.                                           */
  for (ee = 1; ee <= wincon->vc1; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    midx[es] = wincon->mdx[es];
    sum[es] += windat->u1[e] * windat->dzu1[e] * midx[es];
  }

  /* Open boundary normal velocities. Velocities above the surface   */
  /* are set to zero, so it is possible to integrate over the whole  */
  /* water column rather than just to the free surface.              */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (ee = 1; ee <= open[bn]->no3_e1; ee++) {
	e = open[bn]->obc_e1[ee];
	es = window->m2de[e];
	midx[es] = wincon->mdx[es];
	sum[es] += windat->u1[e] * windat->dzu1[e] * midx[es];
      }
    }
  }

  /* Calculate the velocity adjustment                               */
  for (ee = 1; ee <= wincon->vcs1; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    adjust[es] =
      (windat->u1flux[es] -
       sum[es] * windat->dt * window->h1au1[es]) / (depth[es] *
                                                    windat->dt * midx[es] *
                                                    window->h1au1[es]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (ee = 1; ee <= open[bn]->no2_e1; ee++) {
	es = open[bn]->obc_e1[ee];
	depth[es] -= window->botzu1[es];
	adjust[es] =
	  (windat->u1flux[es] -
	   sum[es] * windat->dt * window->h1au1[es]) / (depth[es] *
							windat->dt *
							midx[es] *
							window->h1au1[es]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (ee = 1; ee <= wincon->vc1; ee++) {
    e = wincon->s1[ee];
    es = window->m2de[e];
    windat->u1[e] += adjust[es];
  }

  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (ee = 1; ee <= open[bn]->no3_e1; ee++) {
	e = open[bn]->obc_e1[ee];
	es = window->m2de[e];
	windat->u1[e] += adjust[es];
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (ee = 1; ee <= wincon->ncdry_e1; ee++) {
    e = es = wincon->cdry_e1[ee];
    windat->u1av[es] = 0.0;
    windat->u1flux[es] = 0.0;
    while (e != window->zm1e[e]) {
      windat->u1[e] = 0.0;
      e = window->zm1e[e];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u1bot, 0, window->szeS * sizeof(double));
  for (ee = 1; ee <= wincon->vcs1; ee++) {
    e = wincon->s1[ee];         /* 3D surface coordinate             */
    es = window->m2de[e];       /* 2D coordinate                     */
    eb = wincon->i5[ee];        /* 3D bottom coordinate              */
    windat->u1bot[es] = windat->u1b[eb] - windat->u1avb[es];
  }

  /* Set the velocity for cells one cell deep                        */
  for (ee = 1; ee <= wincon->ncbot_e1; ee++) {
    e = wincon->cbot_e1[ee];
    es = window->m2de[e];
    windat->u1[e] = windat->u1av[es];
    windat->u1bot[es] = 0.0;
  }

  /* Apply the Asselin filter to remove the computational mode for   */
  /* tiled coupling.                                                 */
  if (master->obcf & DF_TILE) {
    for (ee = 1; ee <= window->v3_e1; ee++) {
      e = window->w3_e1[ee];
      windat->u1b[e] = windat->u1b[e] + 0.5 * 0.1 *
	(wincon->w10[e] - 2.0 * windat->u1b[e] + windat->u1[e]);
    }
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
  int e, ee;                /* Sparse coodinate / counter            */
  int en, eo;               /* New and old surface sparse coodinates */
  int *oldsur;              /* Sparse coordinate of old surface      */

  /* Set pointers and initialise                                     */
  oldsur = wincon->i7;

  memcpy(oldsur, window->sur_e1, window->szeS * sizeof(int));

  /* Get the location of the new surface. Note; by including thin    */
  /* layers in sur_e1 it appears as if elevation has risen if the    */
  /* layer was thin at the start of the time-step. The thin layer    */
  /* now will receive a copy of the (2D adjusted) velocity from the  */
  /* layer below.                                                    */
  set_sur_u1(window, windat, wincon);

  /* Set the velocities at newly wetted u1 cells                     */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];

    /* Set a no-gradient  velocity condition  above the new surface. */
    /* When sea  level rises  into new layers  it is possible that a */
    /* cell centre at c (upstream) is  wet but at xm1 (downstream)   */
    /* is dry. This makes the difference between  cell c and xm1 in  */
    /* tensor momentum advection (which is not invoked as we use the */
    /* vector invariant approach) large, which  can lead to          */
    /* velocities that  are  unrealistically  large. A  no-gradient  */
    /* condition above  the  surface  minimises  this   effect  and  */
    /* maintains stability.                                          */
    en = window->sur_e1[ee];
    while (e < en) {
      windat->u1[e] = windat->u1[oldsur[ee]];
      e = window->zm1e[e];
    }
    /* Set newly wetted cells equal to velocity at the old surface   */
    e = eo = oldsur[ee];
    if (e > en)
      e = window->zp1e[e];
    while (e > en) {
      windat->u1[e] = windat->u1[eo];
      e = window->zp1e[e];
    }
    windat->u1[e] = windat->u1[eo];
  }

  /* Set the velocity at solid edges above explicit maps equal zero  */
  if (window->sm_e1) {
    int eb;
    for (ee = 1; ee <= window->b2_e1; ee++) {
      e = window->w2_e1[ee];
      eb = window->sm_e1[e];
      if (eb) {
	while (e < eb) {
	  windat->u1[e] = 0.0;
	  e = window->zm1e[e];
	}
      }
    }
  }
}

/* END set_new_cells_u1()                                            */
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
  int e, ee, es, j;   /* Edge coodinate / counter                    */
  int *bottom;        /* Bottom sparse coordinate                    */
  int *sur;           /* Minimum of surface coordinate               */
  int *nsur;          /* Surface sparse coordinate after the 2D mode */
  double fctop;       /* Flux at the top edge of the cell            */
  double fcbot;       /* Flux at the bottom edge of the cell         */
  int zp1, zm1;       /* Sparse coordinate below cell c              */
  double hf, vf, cnt; /* Diagnostics for continuity                  */
  int dc = 0;         /* Diagnostic continuity surface coordinate    */
  int dw = 1;         /* Window for diagnostic continuity            */
  double detadt;      /* Rate of change of elevation                 */
  double *deta;       /* Elevation rate of change for SIGMA          */
  double d1;          /* Dummy                                       */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  memset(windat->w, 0, window->szc * sizeof(double));
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
        zp1 = window->zp1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	fctop = fcbot + windat->waterss[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fctop -= window->eSc[j][cs] * windat->u1flux3d[e];
	}

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
	  d1 = 0.0;
	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][cs];
	    d1 += window->eSc[j][cs] * windat->u1flux[e];
	  }
          d1 = wincon->oldeta[cs] - (d1 - windat->waterss[cs] * windat->dt) /
	    window->cellarea[cs];

          /* Elevation difference shound be around 1e-17             */
          emstag(LDEBUG,"hd:vel3d:vel_w_update", "ele %f %e\n", windat->t / 86400,
                  windat->eta[cs] - d1);
          printf("ele %f %d %e\n", windat->t / 86400, cs, windat->eta[cs] - d1);

          /* Check that the vertical integral of the 3D fluxes is    */
          /* equal to the 2D fluxes (dc must be a 2D coordinate).    */
          check_flux(window, windat, wincon, dw, window->m2d[dc]);

          /* Calculate and print the continuity balance if required, */
	  /* i.e. the sum of horizontal divergence in the surface    */
	  /* layer plus the vertical flux into the surface layer     */
	  /* should equal the change in elevation over the timestep. */
          /* Output value should be around 1e-8.                     */
	  hf = 0.0;
	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][zp1];
	    hf += (double)window->eSc[j][cs] * windat->u1flux3d[e];
	    /*printf("u1flux3d %d %d %f\n",j,window->wse[e],windat->u1flux3d[e]);*/
	  }
          hf *= windat->dt;
          vf =
            (0.0 - windat->w[zp1]) * windat->dt * window->cellarea[cs] -
	    windat->waterss[zp1] * windat->dt;
          cnt =
            window->cellarea[cs] * (wincon->oldeta[cs] - windat->eta[cs]) -
            (hf + vf);
	  /*printf("fluxes %d %f %f %f %f %f\n",zp1,wincon->oldeta[cs],windat->eta[cs], hf, vf, windat->w[zp1]);*/
          emstag(LDEBUG,"hd:vel3d:vel_w_update", "cons %f : %e\n", windat->t / 86400, cnt);
          printf("cons %f : %e\n", windat->t / 86400, cnt);
        }
        /* Sub-surface diagnostic checks.                            */
        if (window->wn == dw && c == dc) {
          /* Calculate and print the continuity balance if required. */
          /* Output value should be around 1e-8.                     */
          double wtop = (c == cs) ? windat->wtop[cs] : windat->w[zp1];
	  hf = 0.0;
	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][c];
	    hf += window->eSc[j][cs] * windat->u1flux3d[e];
	  }
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
        zm1 = window->zm1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	fcbot = fctop + deta[c] + windat->waterss[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fcbot += window->eSc[j][cs] * windat->u1flux3d[e];
	}
        /* Velocity at cell top                                      */
        windat->w[c] = fcbot / window->cellarea[cs];

        /* Diagnostic checks.                                        */
        if (window->wn == dw && cs == dc) {
          /* Calculate and print the continuity balance if required. */
          /* Output value should be around 1e-8.                     */
          double wtop =
            (c == cs) ? windat->wtop[cs] : windat->w[window->zp1[c]];
	  hf = deta[c];
	  for (j = 1; j <= window->npe[cs]; j++) {
	    e = window->c2e[j][c];
	    hf += window->eSc[j][cs] * windat->u1flux3d[e];
	  }
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

  /* Set the vertical velocity at dry edges                          */
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
      if (!(wincon->means & TRANSPORT) && wincon->means & (VOLFLUX|PSSFLUX)) {
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
  double fctop;       /* Flux at the top edge of the cell            */

  double fcbot;       /* Flux at the bottom edge of the cell         */
  int j, e;           /* Edge direction and index                    */
  int zp1, zm1;       /* Sparse coordinate below cell c              */
  int dw = 0;         /* Window for diagnostic continuity            */
  double detadt;      /* Rate of change of elevation                 */
  double *deta;       /* Elevation rate of change for SIGMA          */
  double *oldw;       /* w at start of step                          */
  double *Fx;

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  nsur = wincon->i2;            /* Set in set_dz()                   */
  sur = wincon->i3;             /* Set in set_dz()                   */
  deta = wincon->w1;
  Fx = windat->u1flux3d;
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
        zp1 = window->zp1[c];
        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	fctop = fcbot + (windat->waterss[c] * window->cellarea[cs]);
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fctop -= window->eSc[j][cs] * Fx[e];
	}
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
	fctop = fcbot + (windat->waterss[c] * window->cellarea[cs]);
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fctop -= window->eSc[j][cs] * Fx[e];
	}
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
        zp1 = window->zp1[c];
        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	cons = (windat->w[zp1] - windat->w[c] - windat->waterss[c]) * 
	  window->cellarea[cs];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  cons -= window->eSc[j][cs] * Fx[e];
	}
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
        zm1 = window->zm1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	fcbot = fctop + deta[c] + windat->waterss[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fcbot -= window->eSc[j][cs] * windat->u1flux3d[e];
	}

        /* Velocity at cell top                                      */
        windat->w[c] = fcbot / window->cellarea[cs];

        /* Transfer top flux to bottom for next cell                 */
        fctop = fcbot;
        c = zm1;
      }
    }
  }

  /* Set the vertical velocity at dry edges                          */
  /* set_w_dry(window, windat, wincon); */

  /* Set the open boundary condition if required                     */
  /* bdry_w(window, windat, wincon); */

}

/* END ff_sl_w_update()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the vertical velocity through continuity.              */
/*-------------------------------------------------------------------*/
void vel_w_trans(geometry_t *window, /* Window geometry              */
		 window_t *windat,   /* Window data                  */
		 win_priv_t *wincon  /* Window constants             */
  )
{
  int c, cc, cs;      /* Sparse coodinate / counter                  */
  int e, ee, es, j;   /* Edge coodinate / counter                    */
  int *bottom;        /* Bottom sparse coordinate                    */
  int *sur;           /* Minimum of surface coordinate               */
  int *nsur;          /* Surface sparse coordinate after the 2D mode */
  double fctop;       /* Flux at the top edge of the cell            */
  double fcbot;       /* Flux at the bottom edge of the cell         */
  int zp1, zm1;       /* Sparse coordinate below cell c              */
  double hf, vf, cnt; /* Diagnostics for continuity                  */
  int dc = 0;         /* Diagnostic continuity surface coordinate    */
  int dw = 1;         /* Window for diagnostic continuity            */
  double detadt;      /* Rate of change of elevation                 */
  double *deta;       /* Elevation rate of change for SIGMA          */
  double d1;          /* Dummy                                       */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  memset(windat->wm, 0, window->szc * sizeof(double));
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
  vel_w_bounds_tran(window, windat, wincon);

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
      windat->wm[c] = windat->wbot[cs];
      fcbot = windat->wm[c] * window->cellarea[cs];

      while (c > sur[cc]) {
        zp1 = window->zp1[c];

        /* Get the flux at cell top due to inflow from bottom, sides */
        /* and any additional water from sources or sinks.           */
	fctop = fcbot + windat->waterss[c];
	for (j = 1; j <= window->npe[cs]; j++) {
	  e = window->c2e[j][c];
	  fctop -= window->eSc[j][cs] * windat->u1vm[e];
	}

        /* Velocity at cell top                                      */
        windat->wm[zp1] = fctop / window->cellarea[cs];

        /* Transfer top flux to bottom for next cell                 */
        fcbot = fctop;
        c = zp1;
      }

      /* Calculations for those cells through which the surface      */
      /* moves; use value from wtop calculated above.                */
      while (c > nsur[cc]) {
        zp1 = window->zp1[c];
        windat->wm[zp1] = windat->wtop[cs];
        c = zp1;
      }
    }
  }
}

/* END vel_w_trans()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets vertical velocity above the surface to wtop. This aids       */
/* stability in the vertical momentum advection when cells have one  */
/* edge dry.                                                         */
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
  int j, jo;                    /* Edge counter                      */
  int c, cc;                    /* Sparse coodinate / counter        */
  int cs;                       /* 2D sparse coordinate              */
  int cb;                       /* Bottom sparse coordinate          */
  int ee, e, es, eb;            /* Edge coordinates                  */
  int c1, c2;                   /* Centre coordinates                */
  int *bottom;                  /* Bottom sparse coordinate          */
  double eta_l;                 /* Upstream surface elevation        */
  double eta_r;                 /* Downstream surface elevation      */
  double *u1top;                /* Surface e1 velocity               */
  double *u1bot;                /* Bottom e1 velocity                */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  u1top = wincon->d1;
  u1bot = wincon->d3;
  memset(windat->wtop, 0, window->szcS * sizeof(double));
  memset(windat->wbot, 0, window->szcS * sizeof(double));
  memset(u1top, 0, window->szeS * sizeof(double));
  memset(u1bot, 0, window->szeS * sizeof(double));

  /* SIGMA : zero velocity at the sigma boundaries                   */
  if (wincon->sigma)
    return;

  /* Set the surface velocity. Note: this cannot be retreived from   */
  /* the u1 array since neighbours are required and with a cell      */
  /* drying the neighbour may not be the surface.                    */
  for (ee = 1; ee <= window->a2_t; ee++) {
    e = window->w2_t[ee];
    es = window->m2de[e];
    u1top[es] = windat->u1[e];
  }
  for (ee = 1; ee <= window->b2_e1; ee++) {
    eb = window->bot_e1[ee];  /* 3D bottom coordinate                */
    es = window->m2de[eb];
    u1bot[es] = windat->u1[eb];
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
  memcpy(windat->wtop, windat->detadt, window->szcS * sizeof(double));
  memset(windat->wbot, 0.0, window->szcS * sizeof(double));

  if (!wincon->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.                                   */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->i3[cc];
    cs = window->m2d[c]; 

    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][cs];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      /* Calculate the surface vertical velocity                     */
      eta_l = windat->eta[c1];
      eta_r = windat->eta[c2];
      windat->wtop[cs] += (u1top[es] * (eta_r - eta_l) / 
			   (window->h2au1[es] * (double)window->npe[cs]));

      /* Calculate the bottom vertical velocity                      */
      windat->wbot[cs] -= (u1bot[es] *	window->dHde1[es] / 
			   (window->h2au1[es]* (double)window->npe[cs]));
    }
  }
}

/* END vel_w_bounds()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the surface and bottom vertical velocity     */
/*-------------------------------------------------------------------*/
void vel_w_bounds_tran(geometry_t *window, /* Window geometry        */
		       window_t *windat,   /* Window data            */
		       win_priv_t *wincon  /* Window constants       */
  )
{
  int j, jo;                    /* Edge counter                      */
  int c, cc;                    /* Sparse coodinate / counter        */
  int cs;                       /* 2D sparse coordinate              */
  int cb;                       /* Bottom sparse coordinate          */
  int ee, e, es, eb;            /* Edge coordinates                  */
  int c1, c2;                   /* Centre coordinates                */
  int *bottom;                  /* Bottom sparse coordinate          */
  double eta_l;                 /* Upstream surface elevation        */
  double eta_r;                 /* Downstream surface elevation      */
  double *u1top;                /* Surface e1 velocity               */
  double *u1bot;                /* Bottom e1 velocity                */

  /*-----------------------------------------------------------------*/
  /* Set pointers and initialise                                     */
  bottom = wincon->i1;          /* Set in set_dz()                   */
  u1top = wincon->d1;
  u1bot = wincon->d3;
  memset(windat->wtop, 0, window->szcS * sizeof(double));
  memset(windat->wbot, 0, window->szcS * sizeof(double));
  memset(u1top, 0, window->szeS * sizeof(double));
  memset(u1bot, 0, window->szeS * sizeof(double));

  /* SIGMA : zero velocity at the sigma boundaries                   */
  if (wincon->sigma)
    return;

  /* Set the surface velocity. Note: this cannot be retreived from   */
  /* the u1 array since neighbours are required and with a cell      */
  /* drying the neighbour may not be the surface.                    */
  for (ee = 1; ee <= window->a2_t; ee++) {
    e = window->w2_t[ee];
    es = window->m2de[e];
    u1top[es] = windat->ume[e];
  }
  for (ee = 1; ee <= window->b2_e1; ee++) {
    eb = window->bot_e1[ee];  /* 3D bottom coordinate                */
    es = window->m2de[eb];
    u1bot[es] = windat->ume[eb];
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
  memcpy(windat->wtop, windat->detadt, window->szcS * sizeof(double));
  memset(windat->wbot, 0.0, window->szcS * sizeof(double));

  if (!wincon->nonlinear)
    return;

  /* Non-linear case - have to calculate the surface slope terms and */
  /* add or subtract as necessary.                                   */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->i3[cc];
    cs = window->m2d[c]; 

    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][cs];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      /* Calculate the surface vertical velocity                     */
      eta_l = windat->eta[c1];
      eta_r = windat->eta[c2];
      windat->wtop[cs] += (u1top[es] * (eta_r - eta_l) / 
			   (window->h2au1[es] * (double)window->npe[cs]));

      /* Calculate the bottom vertical velocity                      */
      windat->wbot[cs] -= (u1bot[es] *	window->dHde1[es] / 
			   (window->h2au1[es]* (double)window->npe[cs]));
    }
  }
}

/* END vel_w_bounds_tran()                                           */
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
int implicit_vdiff(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon, /* Window constants           */
		   double *var,   /* Array to diffuse                */
		   double *nvar,  /* Output array                    */
		   double *Kzin,  /* Mixing coefficient values       */
		   double *dzin,  /* Cell thicknesses                */
		   double *fb,    /* Flux out of bottom              */
		   double *ft,    /* Flux out of top                 */
		   int *ctp,      /* Surface cells to process        */
		   int *cbm,      /* Bottom cells to process         */
		   int vcs,       /* Last index of surface cells     */
		   double *depth, /* SIGMA : Water depth             */
		   double *scale, /* SIGMA : (depth) scaling         */
		   double *w
  )
{
  int e, c, k;               /* Sparse coordinate                    */
  int es, cs, ks, ks1;       /* Surface sparse coordinate            */
  int eb, cb, kb;            /* Bottom sparse coordinate             */
  int zm1;                   /* Sparse cell below c                  */
  int zp1;                   /* Sparse cell above c                  */
  int ee;                    /* Sparse coordinate counter            */
  int e2;                    /* 2D cell corresponding to 3D location */
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
  for (ee = 1; ee <= vcs; ee++) {
    es = e = ctp[ee];     /* Set cs to the surface sparse coordinate */
    eb = cbm[ee];         /* Bottom sparse coordinate                */
    e2 = window->m2de[e]; /* 2D sparse location corresponding to e   */
    zm1 = window->zm1e[e];
    /*
    c = window->e2c[e][0];
    cs = window->e2c[es][0];
    cb = window->e2c[eb][0];
    */
    c = cs = window->e2ijk[e];
    cb = window->e2ijk[eb];

    dzs = dzin[e];
    ks = ks1 = window->s2k[cs];
    kb = window->s2k[cb];
    thinf = 0;

    /*---------------------------------------------------------------*/
    /* Merge thin  layers with layer  below if  required by  setting */
    /* velocity in the  layers es and zm1 to a mass weighted mean of */
    /* these  two layers, then  assuming zm1  is the layer below the */
    /* thin layer.                                                   */
    if (wincon->thin_merge && dzin[es] < wincon->hmin) {
      thinf = 1;
      es = zm1;
      ks--;
      zm1 = window->zm1e[zm1];
      vel1 = var[e];
      vel2 = var[es];
      /* Thin layer lies in a single layer -> continue               */
      if (es == zm1)
        continue;
      var[e] = var[es] =
        (var[e] * dzs + var[es] * dzin[es]) / (dzin[es] + dzs);
      dzin[es] += dzs;
    }

    /*---------------------------------------------------------------*/
    /* Single layer  case  (i.e. surface lies  in the  bottom layer) */
    if (zm1 == window->zm1e[zm1] || !dzin[zm1]) {
      nvar[es] +=
        windat->dt * (fb[e2] - ft[e2]) / max(dzin[es], wincon->hmin);
      if (e != es) {
        var[e] = var[es];
        dzin[es] -= dzs;
      }
      continue;
    }

    /*---------------------------------------------------------------*/
    /* Calculate constants  required from  the surface  to  bottom.  */
    /* The sparse coordinate, c, lies in  the sediment at the end of */
    /* the loop.                                                     */
    e = es;
    for (k = ks; k >= kb; k--) {
      zm1 = window->zm1e[e];
      dzface[k] = 0.5 * (dzin[zm1] + dzin[e]);
      /* SIGMA : Precondition mixing coefficient and the cell for    */
      /* the sigma calculation. Not yet impemented or tested.        */
      /*
      dz[k] = dzin[e] * depth[e2];
      Kz[k] = Kzin[e] / depth[e2];
      */
      dz[k] = dzin[e];
      Kz[k] = Kzin[e];
      dzf[k] = dzin[zm1] / (dzin[zm1] + dzin[e]);
      e = zm1;
    }

    /*---------------------------------------------------------------*/
    /* Set up tri-diagonal set of equations.                         */
    /* Bottom layer.                                                 */
    zm1 = e;
    e = eb;
    zp1 = window->zp1e[e];
    dzdt = dz[kb] / dt;

    Cm1[kb] = 0.0;
    Cp1[kb] = -Kz[kb + 1] / dzface[kb + 1] + w[e] * dzf[kb+1];
    C[kb] = dzdt - Cp1[kb] + w[e] - w[zm1];
    rhs[kb] = dzdt * var[eb] + fb[e2];

    /* Mid-water layers                                              */
    zm1 = e;
    e = zp1;
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / dt;
      Cm1[k] = -Kz[k] / dzface[k] + w[zm1] * dzf[k];
      Cp1[k] = -Kz[k + 1] / dzface[k + 1] + w[e] * dzf[k+1];
      C[k] = dzdt - Cm1[k] - Cp1[k] + w[e];
      Cm1[k] -= w[zm1];
      rhs[k] = dzdt * var[e];
      zm1 = e;
      e = window->zp1e[e];
    }

    /* Surface layer                                                 */
    dzdt = dz[ks] / dt;
    Cm1[ks] = -Kz[ks] / dzface[ks] + w[zm1] * dzf[ks];
    Cp1[ks] = 0.0;
    C[ks] = dzdt - Cm1[ks] + w[e];
    rhs[ks] = dzdt * var[es] - ft[e2];
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
    dvar[ks] = sol[ks] - var[es];
    e = es;
    for (k = ks - 1; k >= kb; k--) {
      e = window->zm1e[e];
      sol[k] -= ud[k + 1] * sol[k + 1];
      dvar[k] = sol[k] - var[e];
    }

    /*---------------------------------------------------------------*/
    /* Reset the surface dz for thin layers                          */
    e = ctp[ee];
    if (e != es) {
      dzin[es] -= dzs;
      dvar[ks1] = dvar[ks];
    }

    /*---------------------------------------------------------------*/
    /* Update the variable                                           */
    e = ctp[ee];
    /* Restore the backward velocities to their un-merged values if  */
    /* merging occurred. If the bottom lies in the layer below the   */
    /* surface, the u1b values are used in vdiff_u2() to get the     */
    /* bottom stress, and the un-merged values are to be used.       */
    if (thinf) {
      var[e] = vel1;
      var[window->zm1e[e]] = vel2;
    }
    for (k = ks1; k >= kb; k--) {
      /*nvar[e] += (dvar[k] * scale[e2]);*/
      nvar[e] += (dvar[k]);
      e = window->zm1e[e];
    }
  }
  return(0);
}

/* END implicit_vdiff()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes geostrophically balanced velocity                        */
/*-------------------------------------------------------------------*/
void calc_geostrophic(geometry_t *window, /* Window geometry         */
		      window_t *windat,   /* Window data             */
		      win_priv_t *wincon  /* Window constants        */
		      )
{
  int e, es, ee, eoe, n;
  double nut;

  reset_map_t_all(window);

  cells2process_e1(window, windat, wincon);
  precalc_u1(window, windat, wincon);
  set_map_e1(window);
  pressure_u1(window, windat, wincon);
  extract_u1_3d(window, windat, wincon);

  for (ee = 1; ee < window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    nut = 0.0;
    for (n = 1; n <= window->nee[es]; n++) {
      eoe = window->eSe[n][e];
      if (!eoe) continue;
      nut += window->wAe[n][e] * windat->nu1[eoe];
    }
    windat->u1[e] = nut / (windat->dt * wincon->u1c5[es]);
  }

  /* Save the velocities for boundary ramping                        */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->to3_e1; ee++) {
      e = open->obc_e1[ee];
      open->u1d[ee] = windat->u1[e];
    }
  }
  set_map_t(window);
}

/* END calc_geostrophic()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs a momentum balance from tendencies                        */
/*-------------------------------------------------------------------*/
void mom_balance(geometry_t *window, /* Window geometry              */
		 window_t *windat,   /* Window data                  */
		 win_priv_t *wincon  /* Window constants             */
		 )
{
  int e, ee;
  double v, v1, v2, v3, v4, v5, v6, vm, vt;

  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    v1 = fabs(windat->u1_adv[e]);
    v2 = fabs(windat->u1_hdif[e]);
    v3 = fabs(windat->u1_vdif[e]);
    v4 = fabs(windat->u1_cor[e]);
    v5 = fabs(windat->u1_btp[e]);
    v6 = fabs(windat->u1_bcp[e]);

    vt = v1 + v2 + v3 + v4 + v5 + v6;
    v = 1.0;                    /* 1 = advection                     */
    vm = v1;
    v = (v2 > vm) ? 2.0 : v;    /* 2 = horizontal diffusion          */
    vm = (v2 > vm) ? v2 : vm;
    v = (v3 > vm) ? 4.0 : v;    /* 4 = vertical diffusion            */
    vm = (v3 > vm) ? v3 : vm;
    v = (v4 > vm) ? 8.0 : v;    /* 8 = Coriolis                      */
    vm = (v4 > vm) ? v4 : vm;
    v = (v5 > vm) ? 16.0 : v;   /* 16 = barotropic pressure          */
    vm = (v5 > vm) ? v5 : vm;
    v = (v6 > vm) ? 32.0 : v;   /* 32 = baroclinic pressure          */
    vm = (v6 > vm) ? v6 : vm;
    windat->mom_bal[e] = v;
  }
}

/* END mom_balance()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the stokes velocity profile at edges                      */
/*-------------------------------------------------------------------*/
void get_sdc_e1(geometry_t *window, /* Window geometry               */
		window_t *windat,   /* Window data                   */
		win_priv_t *wincon  /* Window constants              */
		)
{
  int e, ee, es;
  int c1, c2;
  double w;                                      /* Angular freq     */
  double k;                                      /* Wave number      */
  double uso, u1, u2, depth, period;
  double ramp = (wincon->rampf & STOKES) ? windat->rampval : 1.0;

  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[es][0];
    c2 = window->e2c[es][1];
    period = 0.5 * (windat->wave_period[c1] + windat->wave_period[c2]);
    w = (period) ? 2.0 * PI / period : 0.0;
    k = w * w / wincon->g;
    uso = ramp * vel_c2e(window, windat->wave_ste1, windat->wave_ste2, e);
    while (e != window->zm1e[e]) {
      depth = min(0.0, window->cellz[c1]);
      wincon->w1[e] = uso * exp(2.0 * k * depth);
      e = window->zm1e[e];
    }
  }
}

/* END get_sdc_e1()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Tidal energt extraction due to tide turbines                      */
/*-------------------------------------------------------------------*/
void turbines(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int n, j, c, cs, e, ed, es;
  int jmin;
  double cext;
  double r0 = 1024.0;
  double vel, fe;
  double *nvel, *u1, *u2;
  double dir, d, dmin;

  if (wincon->mode2d) {
    nvel = wincon->w4;
    u1 = windat->u1av;
    u2 = windat->u2av;
  } else {
    nvel = windat->nu1;
    u1 = windat->u1;
    u2 = windat->u2;
  }

  /* Loop over all turbines in the window                            */
  for (n = 0; n < wincon->nturb; n++) {
    c = wincon->turb[n];
    cs = window->m2d[c];
    cext = wincon->cturb[n];

    /* Get the direction of the cell centered velocity vector        */
    dir =  fmod(atan2(windat->u[c], windat->v[c]) * 180.0 / M_PI + 180.0, 360.0);

    /* Get the edge which aligns most closely with the cell centered */
    /* vector.                                                       */
    dmin = HUGE;
    for (j = 1; j <= window->npe[cs]; j++) {
      e = window->c2e[j][c];
      es = window->m2de[e];
      d = fabs(window->thetau1[es] - dir);
      if (d < dmin) {
	dmin = d;
	jmin = j;
      }
    }

    /* Remove energy from all edges surrounding cell c               */
    for (j = 1; j <= window->npe[cs]; j++) {
      /*if (j != jmin) continue;*/
      e = window->c2e[j][c];
      es = window->m2de[e];
      ed = (wincon->mode2d) ? es : e;
      
      /* Get the current speed at the edge                           */
      vel = sqrt(u1[ed] * u1[ed] * u2[ed] * u2[ed]);

      /* Get the energy loss                                         */
      fe = -0.5 * cext * r0 * vel * u1[ed];

      /* Update the velocity                                         */
      nvel[e] += windat->dt * fe;
      wincon->u1inter[es] += fe * windat->dzu1[e];
    }
  }
}

/* END turbines                                                      */
/*-------------------------------------------------------------------*/


double gridze(geometry_t *geom, int e)
{
  int c1 = geom->e2c[e][0];
  int c2 = geom->e2c[e][1];
  double g1 = geom->gridz[c1];
  double g2 = geom->gridz[c2];
  return(max(g1, g2));
}
