/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/dp.h
 *  
 *  Description:
 *  Manage distributed processing.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dp.h 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#ifndef _DP_H
#define _DP_H

typedef struct dp_window {
  int window_id;
  master_t *master;
  geometry_t *geom;
  window_t *windata;
  win_priv_t *wincon;
  void *private_data;
} dp_window_t;


typedef struct dp_details {
  int nwindows;                 /* Number of windows */
  dp_window_t *dp_windows;      /* Array of windows */

  void (*init) (dp_window_t *dpw);
  void (*cleanup) (dp_window_t *dpw);

  void (*vel2d_step_p1) (dp_window_t *dpw);
  void (*vel2d_gather_step_p1) (dp_window_t *dpw);
  void (*vel2d_step_p2) (dp_window_t *dpw);
  void (*vel2d_gather_step_p2) (dp_window_t *dpw);
  void (*vel3d_step_p1) (dp_window_t *dpw);
  void (*vel3d_gather_step_p1) (dp_window_t *dpw);
  void (*vel3d_step_p2) (dp_window_t *dpw);
  void (*vel3d_gather_step_p2) (dp_window_t *dpw);
  void (*vel3d_post_p1) (dp_window_t *dpw);
  void (*vel3d_gather_post_p1) (dp_window_t *dpw);
  void (*vel3d_post_p2) (dp_window_t *dpw);
  void (*vel3d_gather_post_p2) (dp_window_t *dpw);
  void (*tracer_step) (dp_window_t *dpw);
  void (*tracer_gather_step) (dp_window_t *dpw);

} dp_details_t;

/* Distributed processing prototypes */
void dp_init(master_t *master,
             geometry_t *geom[],
             window_t *windata[], win_priv_t *wincon[], int nwindows);
double dp_clock();
void dp_vel2d_step_p1();
void dp_vel2d_step_p2();
void dp_vel3d_step_p1();
void dp_vel3d_step_p2();
void dp_vel3d_post_p1();
void dp_vel3d_post_p2();
void dp_tracer_step();
void dp_cleanup();

/* Protected methods. To be used by dp.c only! */
void dp_none_init(dp_window_t *dpw);
void dp_none_cleanup(dp_window_t *dpw);
void dp_none_vel2d_step_p1(dp_window_t *dpw);
void dp_none_vel2d_gather_step_p1(dp_window_t *dpw);
void dp_none_vel2d_step_p2(dp_window_t *dpw);
void dp_none_vel2d_gather_step_p2(dp_window_t *dpw);
void dp_none_vel3d_step_p1(dp_window_t *dpw);
void dp_none_vel3d_gather_step_p1(dp_window_t *dpw);
void dp_none_vel3d_step_p2(dp_window_t *dpw);
void dp_none_vel3d_gather_step_p2(dp_window_t *dpw);
void dp_none_vel3d_post_p1(dp_window_t *dpw);
void dp_none_vel3d_gather_post_p1(dp_window_t *dpw);
void dp_none_vel3d_post_p2(dp_window_t *dpw);
void dp_none_vel3d_gather_post_p2(dp_window_t *dpw);
void dp_none_tracer_step(dp_window_t *dpw);
void dp_none_tracer_gather_step(dp_window_t *dpw);

#if defined(HAVE_PTHREADS)
void dp_pthread_init(dp_window_t *dpw);
void dp_pthread_cleanup(dp_window_t *dpw);
void dp_pthread_vel2d_step_p1(dp_window_t *dpw);
void dp_pthread_vel2d_gather_step_p1(dp_window_t *dpw);
void dp_pthread_vel2d_step_p2(dp_window_t *dpw);
void dp_pthread_vel2d_gather_step_p2(dp_window_t *dpw);
void dp_pthread_vel3d_step_p1(dp_window_t *dpw);
void dp_pthread_vel3d_gather_step_p1(dp_window_t *dpw);
void dp_pthread_vel3d_step_p2(dp_window_t *dpw);
void dp_pthread_vel3d_gather_step_p2(dp_window_t *dpw);
void dp_pthread_vel3d_post_p1(dp_window_t *dpw);
void dp_pthread_vel3d_gather_post_p1(dp_window_t *dpw);
void dp_pthread_vel3d_post_p2(dp_window_t *dpw);
void dp_pthread_vel3d_gather_post_p2(dp_window_t *dpw);
void dp_pthread_tracer_step(dp_window_t *dpw);
void dp_pthread_tracer_gather_step(dp_window_t *dpw);
#endif

#if defined(HAVE_MPI)
void dp_mpi_init(dp_window_t *dpw);
void dp_mpi_cleanup(dp_window_t *dpw);
void dp_mpi_vel2d_step_p1(dp_window_t *dpw);
void dp_mpi_vel2d_gather_step_p1(dp_window_t *dpw);
void dp_mpi_vel2d_step_p2(dp_window_t *dpw);
void dp_mpi_vel2d_gather_step_p2(dp_window_t *dpw);
void dp_mpi_vel3d_step_p1(dp_window_t *dpw);
void dp_mpi_vel3d_gather_step_p1(dp_window_t *dpw);
void dp_mpi_vel3d_step_p2(dp_window_t *dpw);
void dp_mpi_vel3d_gather_step_p2(dp_window_t *dpw);
void dp_mpi_vel3d_post_p1(dp_window_t *dpw);
void dp_mpi_vel3d_gather_post_p1(dp_window_t *dpw);
void dp_mpi_vel3d_post_p2(dp_window_t *dpw);
void dp_mpi_vel3d_gather_post_p2(dp_window_t *dpw);
void dp_mpi_tracer_step(dp_window_t *dpw);
void dp_mpi_tracer_gather_step(dp_window_t *dpw);
#endif

#endif                          /* _DP_H */
