/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/master/hd_step.c
 *  
 *  Description:
 *  hydrodynamic step.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: hd_step.c 6469 2020-02-18 23:49:45Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "hd.h"


void custom_step(hd_data_t* hdata);

/*-------------------------------------------------------------------*/
/* Step the grid forward in time until tstop is reached. This may    */
/* take several grid timesteps.                                      */
/*-------------------------------------------------------------------*/
void hd_step(hd_data_t *hd_data, double tstop)
{
  geometry_t *geom = hd_data->geom;
  master_t *master = hd_data->master;
  geometry_t **window = hd_data->window;
  window_t **windat = hd_data->windat;
  win_priv_t **wincon = hd_data->wincon;
  int n, nc = 0;                /* Counter                           */
  double time_left, time_left_tr;
  int trc = master->tratio;     /* Counter for tracer timestep       */
  trc = master->nstep;
  /* Get the current time                                            */
  master->t = schedule->t;
  master->days = master->t / 86400;

  /* Calculate the time remaining                                    */
  time_left = time_left_tr = tstop - master->t;

  /* Loop over the grid timesteps                                    */
  while (!killed && time_left > DT_EPS) {
    /* Calculate maximum possible dt for this step                   */
    double max_dt;
    double max_dt2d;
    int iratio;

    master->is_filled = 0;
    max_dt = min(master->grid_dt, time_left);
    max_dt2d = master->grid_dt / master->iratio;

    /* Initialise the diagnostics file for this timestep             */
    monitor(master, window, 0);

    PRINT_TIMING_HEADER("hd_step", master);

    /* Set the timesteps                                             */
    master->dt = max_dt;
    master->dtu1 = master->dtu2 = master->dt;
    master->dtb = master->dtf;
    master->dtf = master->dt;
    /* Set the tracer time-step                                      */
    /* dttr = tratio * grid_dt if the schedule step (tstop - t) is   */
    /* long enough (i.e. time_left > dttr). When this is the case    */
    /* dttr is set at tratio intervals or the first step of the      */
    /* schedule, otherwise dttr = 0. If time_left becomes less than  */
    /* the tracer step during schedule, then if this occurs on the   */
    /* first step of the schedule dttr = time_left (the length of    */
    /* the schedule) otherwise dttr is the time remaining in the     */
    /* schedule.                                                     */
    if (master->tratio < 1.0) {
      double dttr = master->tratio * master->grid_dt;
      if(time_left >= dttr) {
        if (!nc || trc % (int)master->tratio == 0)
          master->dttr = min(dttr, time_left);
        else
          master->dttr = 0.0;
      } else {
        if (time_left_tr && master->dt < master->grid_dt)
          master->dttr = master->dt;
        else
          master->dttr = !nc ? time_left : 0.0;
      }
      time_left_tr -= master->dttr;
      trc += 1;
    } else if (master->tratio > 1.0) {
      master->dttc += master->dt;
      if (trc % (int)master->tratio == (int)master->tratio - 1) {
	master->dttr = master->dttc;
	master->dttc = 0.0;
      } else
	master->dttr = 0.0;
      trc += 1;
    } else
      master->dttr = master->dt;
    time_left -= master->dt;

    /*---------------------------------------------------------------*/
    /* Calculate appropriate 2d time step.  iratio is the smallest   */
    /* integer which, when divided into dt, gives a value less than  */
    /* max_dt2d.                                                     */
    /* iratio = (int)(master->dt)/max_dt2d + 1; */
    iratio = (int)ceil(master->dt / max_dt2d);
    /* Only even iratio's allowed with the leapfrog scheme so that   */
    /* fluxes summed on odd timesteps add integrally to the 3d step  */
    /* (note : ic=0; ic<iratio).  */
    if (((iratio + 1) % 2 == 0))
      iratio++;
    master->dt2d = master->dt / iratio;
    /* Input time for OBC data must be at the forward time           */
    master->t3d = master->t + master->dt;
    /* Evaluate ramp value                                           */
    if (master->t >= master->rampend ||
        master->rampstart >= master->rampend)
      master->rampval = 1.0;
    else if (master->t <= master->rampstart)
      master->rampval = 0.0;
    else
      master->rampval = (1.0 - cos(PI * (master->t - master->rampstart) /
                                   (master->rampend -
                                    master->rampstart))) / 2;

    /*---------------------------------------------------------------*/
    /* Sources and sinks of water                                    */
    TIMING_SET;
    sourcesink(master);
    TIMING_DUMP(1, " sourcesink");

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.              */
#if GLOB_BC
    set_lateral_BC_tr(master->tr_wc, master->ntr, geom->nbpt,
                      geom->bpt, geom->bin);
#endif

    /*---------------------------------------------------------------*/
    /* Calculate a heat and salt flux if required                    */
    /*
     * Now done on the windows - see mode3d_step_window_p1
     * calc_heatf(master);
     * calc_saltf(master);
     */

    /*---------------------------------------------------------------*/
    /* Solve the 3D mode in each window                              */
    TIMING_SET;
    mode3d_step(geom, master, window, windat, wincon, master->nwindows);
    TIMING_DUMP(1," mode3d_step");
    if (master->crf == RS_RESTART) return;
#ifdef HAVE_MPI
    /* Compare Vz solution between single and multiwindows */
    if (mpi_check_multi_windows_Vz(geom, master, window, windat))
      hd_quit("Single and multli-windows Vz mismatch\n");
    
    /* Compare velocity solution between single and multiwindows */
    if (mpi_check_multi_windows_velocity(geom, master, window, windat,
					 VEL3D))
      hd_quit("Single and multli-windows u1 mismatch\n");
#endif

    /*---------------------------------------------------------------*/
    /* Solve the 2D mode in each window                              */
    for (n = 1; n <= master->nwindows; n++)
      windat[n]->iratio = iratio;

    TIMING_SET;
    mode2d_step(geom, master, window, windat, wincon, master->nwindows);
    TIMING_DUMP(1," mode2d_step");
    if (master->crf == RS_RESTART) return;
#ifdef HAVE_MPI
    /* Compare velocity solution between single and multiwindows */
    if (mpi_check_multi_windows_velocity(geom, master, window, windat,
					 VEL2D))
      hd_quit("Single and multli-windows u1av mismatch\n");
#endif

    /*---------------------------------------------------------------*/
    /* Do the post 3D mode calculations                              */
    mode3d_post(geom, master, window, windat, wincon, master->nwindows);
#ifdef HAVE_MPI
    /* Compare velocity solution between single and multiwindows */
    if (mpi_check_multi_windows_velocity(geom, master, window, windat,
					 FLUX))
      hd_quit("Single and multli-windows u1flux3d post mismatch\n");
#endif

    /*---------------------------------------------------------------*/
    /* Set the tracer lateral boundary conditions (no-flux)          */
#if GLOB_BC
    set_lateral_BC_vel(master->u1flux3d, geom->nbpt,
                       geom->bpte1, geom->bine1);
    set_lateral_BC_vel(master->u2flux3d, geom->nbpt,
                       geom->bpte2, geom->bine2);
#endif

    /*---------------------------------------------------------------*/
    /* Solve the tracer equation in each window                      */
    TIMING_SET;
    tracer_step(master, window, windat, wincon, master->nwindows);
    TIMING_DUMP(1," tracer_step");
#ifdef HAVE_MPI
    /* Compare single tracer solution between single and multiwindows */
    if (mpi_check_multi_windows_tracer(geom, master, window, windat))
      hd_quit("Single and multli-windows tracer mismatch\n");
#endif
    if (master->crf == RS_RESTART) return;

    /* call any cyustom steps - always at the end!*/
    custom_step(hd_data);

    /*---------------------------------------------------------------*/
    /* Do particle tracking if required                              */
    pt_update(master, window, windat, wincon);

    master->nstep++;
    master->nstep2d += iratio;
    nc++;

    /* Print the run diagnostics */
    monitor(master, window, 1);
    alerts(master);
    /* Re-configure the window partitioning */
    if (master->win_reset && master->nstep % abs(master->win_reset) == 0) {
      reset_windows(hd_data);
      window = hd_data->window;
      windat = hd_data->windat;
      wincon = hd_data->wincon;
      if (master->win_reset < 0)
	master->win_reset = 0;
    }
  }
}

/* END hd_step()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Main time loop for a 2D depth averaged model                      */
/*-------------------------------------------------------------------*/
void hd_step_2d(hd_data_t *hd_data, double tstop)
{
  geometry_t *geom = hd_data->geom;
  master_t *master = hd_data->master;
  geometry_t **window = hd_data->window;
  window_t **windat = hd_data->windat;
  win_priv_t **wincon = hd_data->wincon;
  int n;
  double time_left;

  /* Get the current time */
  master->t = schedule->t;
  master->days = master->t / 86400;

  /* Calculate the time remaining */
  time_left = tstop - master->t;

  /* Loop over the grid timesteps */
  while (!killed && time_left > DT_EPS) {

    /* Calculate maximum possible dt for this step */
    double max_dt = min(master->grid_dt, time_left);
    double max_dt2d = master->grid_dt / master->iratio;
    int iratio;

    master->is_filled = 0;

    /* Initialise the diagnostics file for this timestep */
    monitor(master, window, 0);

    /* Get the 2D timestep */
    /* Set the timesteps */
    master->dt = max_dt;
    master->dtu1 = master->dtu2 = master->dt;
    time_left -= master->dt;
    /* iratio = (int)(master->dt/max_dt2d) + 1; */
    iratio = (int)ceil(master->dt / max_dt2d);
    if (((iratio + 1) % 2 == 0))
      iratio++;
    master->dt2d = master->dt / iratio;

    /* Evaluate ramp value */
    if (master->t >= master->rampend ||
        master->rampstart >= master->rampend)
      master->rampval = 1.0;
    else if (master->t <= master->rampstart)
      master->rampval = 0.0;
    else
      master->rampval = (1.0 - cos(PI * (master->t - master->rampstart) /
                                   (master->rampend -
                                    master->rampstart))) / 2;

    /*---------------------------------------------------------------*/
    /* Do the bio-chem if required */
    /* bmstep_s(master); */

    /*---------------------------------------------------------------*/
    /* Calculate a heatflux if required */
    /*
     * now done on the windows - see window loop in mode3d_prep
     * calc_heatf(master);
     */

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.  */
#if GLOB_BC
    set_lateral_BC_tr(master->tr_wc, master->ntr, geom->nbpt,
                      geom->bpt, geom->bin);
#endif

    /*---------------------------------------------------------------*/
    /* Solve the 2D mode in each window */
    mode3d_prep(geom, master, window, windat, wincon, master->nwindows);
    for (n = 1; n <= master->nwindows; n++)
      windat[n]->iratio = iratio;
    mode2d_step(geom, master, window, windat, wincon, master->nwindows);

    /*---------------------------------------------------------------*/
    /* Set the 3D velocities and do the tracers */
    tracer_step(master, window, windat, wincon, master->nwindows);

    /*UR-ADDED execute any custome exchanges */
    custom_step(hd_data);

    /* Update density field */
    /* if(master->calc_dens)density(master); */
    master->nstep++;
    master->nstep2d += iratio;

    /* Print the run diagnostics */
    monitor(master, window, 1);
    alerts(master);
  }
}

/* END hd_step_2d()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Step the grid forward in time until tstop is reached. This may    */
/* take several grid timesteps.                                      */
/*-------------------------------------------------------------------*/
void hd_step_trans(hd_data_t *hd_data, double tstop)
{
  master_t *master = hd_data->master;
  geometry_t **window = hd_data->window;
  window_t **windat = hd_data->windat;
  win_priv_t **wincon = hd_data->wincon;
  int nc = 0;                          /* Counter                    */
  double time_left, time_left_tr;
  int trc = master->tratio;     /* Counter for tracer timestep       */
  trc = master->nstep;
  /* Get the current time                                            */
  master->t = schedule->t;
  master->days = master->t / 86400;

  /* Calculate the time remaining                                    */
  time_left = time_left_tr = tstop - master->t;

  /* Loop over the grid timesteps                                    */
  while (!killed && time_left > DT_EPS) {
    /* Calculate maximum possible dt for this step                   */
    double max_dt;
    double max_dt2d;
    int iratio;

    max_dt = min(master->grid_dt, time_left);
    max_dt2d = master->grid_dt / master->iratio;

    /* Initialise the diagnostics file for this timestep             */
    monitor(master, window, 0);

    PRINT_TIMING_HEADER("hd_step_trans", master);

    /* Set the timesteps                                             */
    master->dt = max_dt;
    master->dtu1 = master->dtu2 = master->dt;
    master->dtb = master->dtf;
    master->dtf = master->dt;
    /* Set the tracer time-step                                      */
    /* dttr = tratio * grid_dt if the schedule step (tstop - t) is   */
    /* long enough (i.e. time_left > dttr). When this is the case    */
    /* dttr is set at tratio intervals or the first step of the      */
    /* schedule, otherwise dttr = 0. If time_left becomes less than  */
    /* the tracer step during schedule, then if this occurs on the   */
    /* first step of the schedule dttr = time_left (the length of    */
    /* the schedule) otherwise dttr is the time remaining in the     */
    /* schedule.                                                     */
    if (master->tratio > 1.0) {
      double dttr = master->tratio * master->grid_dt;
      if(time_left >= dttr) {
        if (!nc || trc % (int)master->tratio == 0)
          master->dttr = min(dttr, time_left);
        else
          master->dttr = 0.0;
      } else {
        if (time_left_tr && master->dt < master->grid_dt)
          master->dttr = master->dt;
        else
          master->dttr = !nc ? time_left : 0.0;
      }
      time_left_tr -= master->dttr;
      trc += 1;
    } else if (master->tratio < 1.0) {
      master->dttr = master->dt = master->tratio * master->grid_dt;
    } else
      master->dttr = master->dt;
    time_left -= master->dt;

    /*---------------------------------------------------------------*/
    /* Calculate appropriate 2d time step.  iratio is the smallest   */
    /* integer which, when divided into dt, gives a value less than  */
    /* max_dt2d.                                                     */
    /* iratio = (int)(master->dt)/max_dt2d + 1; */
    iratio = (int)ceil(master->dt / max_dt2d);
    /* Only even iratio's allowed with the leapfrog scheme so that   */
    /* fluxes summed on odd timesteps add integrally to the 3d step  */
    /* (note : ic=0; ic<iratio).  */
    if (((iratio + 1) % 2 == 0))
      iratio++;
    master->dt2d = master->dt / iratio;
    master->t3d = master->t;
    /* Evaluate ramp value                                           */
    if (master->t >= master->rampend ||
        master->rampstart >= master->rampend)
      master->rampval = 1.0;
    else if (master->t <= master->rampstart)
      master->rampval = 0.0;
    else
      master->rampval = (1.0 - cos(PI * (master->t - master->rampstart) /
                                   (master->rampend -
                                    master->rampstart))) / 2;

    /*---------------------------------------------------------------*/
    if (master->tmode & SP_CHECK) {
      trans_data_check(master, window, windat, wincon);
      master->nstep++;
      master->nstep2d += iratio;
      nc++;
      monitor(master, window, 1);
      continue;
    }
    if (master->tmode & TR_CHECK)
      trans_data_check(master, window, windat, wincon);

    /*---------------------------------------------------------------*/
    /* Sources and sinks of water                                    */
    TIMING_SET;
    sourcesink(master);
    TIMING_DUMP(1, " sourcesink");

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.              */
#if GLOB_BC
    set_lateral_BC_tr(master->tr_wc, master->ntr, geom->nbpt,
                      geom->bpt, geom->bin);
#endif

    /*---------------------------------------------------------------*/
    /* Calculate a heat and salt flux if required                    */
    /*
     * Now done on the windows - see code in transport_step after
     *                           win_data_fill_3d
     * calc_heatf(master);
     * calc_saltf(master);
     */

    /*---------------------------------------------------------------*/
    /* Set the tracer lateral boundary conditions (no-flux)          */
#if GLOB_BC
    set_lateral_BC_vel(master->u1flux3d, geom->nbpt,
                       geom->bpte1, geom->bine1);
    set_lateral_BC_vel(master->u2flux3d, geom->nbpt,
                       geom->bpte2, geom->bine2);
#endif

    /*---------------------------------------------------------------*/
    /* Set up the model for transport                                */
    TIMING_SET;
    transport_step(master, window, windat, wincon);
    TIMING_DUMP(1," transport_step");

    /*---------------------------------------------------------------*/
    /* Solve the tracer equation in each window.                     */
    /* Note: Using STREAMLINE mode, the first dump of (p,q,r) is     */
    /* usually zero since the first dump is scheduled before the     */
    /* first hd_step (where (p,q,r) are set). Therefore, global      */
    /* filling may use inaccurate total masses, which may affect the */
    /* remainder of the simulation. In this case, skip the first     */
    /* tracer_step, so that non-zero reset values are read from      */
    /* file.                                                         */
    TIMING_SET;

    if (master->tmode & SP_ORIGIN) {
      /*if (master->t > master->tstart + master->dttr)*/
      tracer_step(master, window, windat, wincon, master->nwindows);
    } else
      tracer_step(master, window, windat, wincon, master->nwindows);

    TIMING_DUMP(1," tracer_step");
    /*---------------------------------------------------------------*/
    /* Finalize the model transport                                  */
    transport_post(master, window, windat, wincon);

    /* call any custom steps - always at the end!*/
    custom_step(hd_data);

    /*---------------------------------------------------------------*/
    /* Do particle tracking if required                              */
    pt_update(master, window, windat, wincon);

    master->nstep++;
    master->nstep2d += iratio;
    nc++;

    /* Update simulation time */
    master->t += master->dt;

    /* Print the run diagnostics */
    monitor(master, window, 1);
    alerts(master);

    /* Re-configure the window partitioning */
    if (master->win_reset && master->nstep % master->win_reset == 0) {
      reset_windows(hd_data);
      window = hd_data->window;
      windat = hd_data->windat;
      wincon = hd_data->wincon;
    }
  }
}

/* END hd_step_trans()                                               */
/*-------------------------------------------------------------------*/


/*UR-ADDED-----------------------------------------------------------*/
/* step the cusstom list of exchanges -------------------------------*/
void custom_step(hd_data_t* hdata)
{
  custom_function_t* custom_fnc;
  if(hdata->master->custom_fnstack != NULL )
  {
    custom_fnc = hdata->master->custom_fnstack->functions;
    while(custom_fnc != NULL)
    {
      if(custom_fnc->gather != NULL)
        custom_fnc->gather(custom_fnc, hdata);

      if(custom_fnc->scatter != NULL)
        custom_fnc->scatter(custom_fnc, hdata);
      custom_fnc = custom_fnc->next;
    }
  }

}
