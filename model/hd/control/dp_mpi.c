/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/control/dp_mpi.c
 *  
 *  Description:
 *  Manage distributed processing using MPI.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dp_mpi.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "hd.h"

#if defined(HAVE_MPI)
#include <mpi.h>

/*
 * The basic initialisation and cleanup is done in main
 *
 */
void dp_mpi_init(dp_window_t *dpw) {
  
  int mpi_sz = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_sz);
  /* Store the window number of this process */
  // rank is a global set in main.c and on master in master_build
  if (dpw->geom->nwindows != mpi_sz)
    hd_quit("MPI_SIZE (%d) does not match the number of windows specified (%d)\n",
	    mpi_sz, dpw->geom->nwindows);
}

void dp_mpi_cleanup(dp_window_t *dpw) {
  /* Finalize command in main */
}

/*
 * step   - effectively scatter, send to other windows
 * gather - receive data and only run one window per process/rank
 *
 * A better way would be to install the calling loop functions in the
 * dp details rather than the body to avoid the if's but its okay for
 * now, given enough work for each window
 */
void dp_mpi_vel2d_step_p1(dp_window_t *dpw)
{
  win_data_buffered_send_fill_2d(dpw->windata, dpw->geom->wn);
}

void dp_mpi_vel2d_gather_step_p1(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode2d_step_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_vel2d_step_p2(dp_window_t *dpw)
{
  win_data_buffered_send_refill_2d(dpw->windata, dpw->geom->wn, NVELOCITY);
}

void dp_mpi_vel2d_gather_step_p2(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode2d_step_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_vel2d_step_p3(dp_window_t *dpw)
{
  win_data_buffered_send_refill_2d(dpw->windata, dpw->geom->wn, ETA_A);
}

void dp_mpi_vel2d_gather_step_p3(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode2d_step_window_p3(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_vel3d_step_p1(dp_window_t *dpw)
{
  win_data_buffered_send_fill_3d(dpw->windata, dpw->geom->wn);
}

void dp_mpi_vel3d_gather_step_p1(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode3d_step_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_vel3d_step_p2(dp_window_t *dpw)
{
  win_data_buffered_send_refill_3d(dpw->windata, dpw->geom->wn, MIXING);
}

void dp_mpi_vel3d_gather_step_p2(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode3d_step_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_vel3d_post_p1(dp_window_t *dpw)
{

}

void dp_mpi_vel3d_gather_post_p1(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode3d_post_window_p1(dpw->master, dpw->geom, dpw->windata, dpw->wincon);

  if (dpw->master->mpi_rank == 0)
    win_data_buffered_recv_empty_3d(dpw->windata, dpw->geom->wn, CFL);
}

void dp_mpi_vel3d_post_p2(dp_window_t *dpw)
{
  win_data_buffered_send_refill_3d(dpw->windata, dpw->geom->wn, VELOCITY);
}

void dp_mpi_vel3d_gather_post_p2(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    mode3d_post_window_p2(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_tracer_step(dp_window_t *dpw)
{
  if (dpw->geom->wn == dpw->master->mpi_rank+1)
    tracer_step_window(dpw->master, dpw->geom, dpw->windata, dpw->wincon);
}

void dp_mpi_tracer_gather_step(dp_window_t *dpw)
{
  if (dpw->master->mpi_rank == 0)
    win_data_buffered_recv_empty_3d(dpw->windata, dpw->geom->wn, TRACERS);
}

#endif
