/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/include/dp.h
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
 *  $Id: dp.h 5841 2018-06-28 06:51:55Z riz008 $
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


/* Various processing modes */
#define DP_NONE        0x0000
#define DP_PTHREADS    0x0001
#define DP_OPENMP      0x0002
#define DP_MPI         0x0004
/* Use buffered transfers */
#define DP_BUFFEREDv1 0x0100
#define DP_BUFFEREDv2 0x0200
#define DP_BUFFEREDv3 0x0400
#define DP_BUFFERED   (DP_BUFFEREDv1|DP_BUFFEREDv2|DP_BUFFEREDv3)

typedef struct dp_details {
  int nwindows;                 /* Number of windows */
  dp_window_t *dp_windows;      /* Array of windows */

  void (*init) (dp_window_t *dpw);
  void (*cleanup) (dp_window_t *dpw);

  void (*vel2d_step_p1) (dp_window_t *dpw);
  void (*vel2d_gather_step_p1) (dp_window_t *dpw);
  void (*vel2d_step_p2) (dp_window_t *dpw);
  void (*vel2d_gather_step_p2) (dp_window_t *dpw);
  void (*vel2d_step_p3) (dp_window_t *dpw);
  void (*vel2d_gather_step_p3) (dp_window_t *dpw);
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
void dp_vel2d_step_p3();
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
void dp_none_vel2d_step_p3(dp_window_t *dpw);
void dp_none_vel2d_gather_step_p3(dp_window_t *dpw);
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
void dp_pthread_vel2d_step_p3(dp_window_t *dpw);
void dp_pthread_vel2d_gather_step_p3(dp_window_t *dpw);
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
void dp_mpi_vel2d_step_p3(dp_window_t *dpw);
void dp_mpi_vel2d_gather_step_p3(dp_window_t *dpw);
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
