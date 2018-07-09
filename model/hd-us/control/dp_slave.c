/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/control/dp_mpi.c
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
 *  $Id: dp_slave.c 5873 2018-07-06 07:23:48Z riz008 $
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


geometry_t mpi_window;


void kill_signal_handler(int i)
{
  killed = 1;
}

int main(int argc, char *argv[])
{

  killed = 0;
  model_running = 0;

  /* set signal to allow clean up on SIGTERM */
  signal(SIGTERM, kill_signal_handler);
  mpi_run();

  exit(0);
}



void mpi_run()
{
  /*
   * Initialise, wait to receive the window to be run
   */
  mpi_init();
  
  
  
}



void mpi_init()
{
  
  
  
  
  
}


void mpi_cleanup()
{
  
}


void mpi_vel2d_step()
{
  
}


void mpi_vel2d_gather_step(dp_window_t *dpw)
{
  
}


void mpi_vel3d_step_p1(dp_window_t *dpw)
{
  
}


void mpi_vel3d_gather_step_p1(dp_window_t *dpw)
{
  
}


void mpi_vel3d_step_p2(dp_window_t *dpw)
{
  
}


void mpi_vel3d_gather_step_p2(dp_window_t *dpw)
{
  
}


void mpi_vel3d_post_p1(dp_window_t *dpw)
{
  
}


void mpi_vel3d_gather_post_p1(dp_window_t *dpw)
{
  
}


void mpi_vel3d_post_p2(dp_window_t *dpw)
{
  
}


void mpi_vel3d_gather_post_p2(dp_window_t *dpw)
{
  
}


void mpi_tracer_step(dp_window_t *dpw)
{
  
}


void mpi_tracer_gather_step(dp_window_t *dpw)
{
  
}

#endif
