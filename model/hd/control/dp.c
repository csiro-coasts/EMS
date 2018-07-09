/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/control/dp.c
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
 *  $Id: dp.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <sys/times.h>
#include "hd.h"

#if defined(HAVE_OMP)
#include <omp.h>
#endif

/* Root node for the distributed processing xytoij_tree_t */
dp_details_t dp;

/* Initialise the distributed processing */
void dp_init(master_t *master,
             geometry_t *geom[],
             window_t *windata[], win_priv_t *wincon[], int nwindows)
{

  int nn, n;

  /* Setup the threading strategy (e.g. none, pthreads, etc). */
#if defined(HAVE_PTHREADS)
  if (strcasecmp(master->params->dp_mode, "pthreads") == 0) {
    dp.init = dp_pthread_init;
    dp.cleanup = dp_pthread_cleanup;
    dp.vel2d_step_p1 = dp_pthread_vel2d_step_p1;
    dp.vel2d_gather_step_p1 = dp_pthread_vel2d_gather_step_p1;
    dp.vel2d_step_p2 = dp_pthread_vel2d_step_p2;
    dp.vel2d_gather_step_p2 = dp_pthread_vel2d_gather_step_p2;
    dp.vel2d_step_p3 = dp_pthread_vel2d_step_p3;
    dp.vel2d_gather_step_p3 = dp_pthread_vel2d_gather_step_p3;
    dp.vel3d_step_p1 = dp_pthread_vel3d_step_p1;
    dp.vel3d_gather_step_p1 = dp_pthread_vel3d_gather_step_p1;
    dp.vel3d_step_p2 = dp_pthread_vel3d_step_p2;
    dp.vel3d_gather_step_p2 = dp_pthread_vel3d_gather_step_p2;
    dp.vel3d_post_p1 = dp_pthread_vel3d_post_p1;
    dp.vel3d_gather_post_p1 = dp_pthread_vel3d_gather_post_p1;
    dp.vel3d_post_p2 = dp_pthread_vel3d_post_p2;
    dp.vel3d_gather_post_p2 = dp_pthread_vel3d_gather_post_p2;
    dp.tracer_step = dp_pthread_tracer_step;
    dp.tracer_gather_step = dp_pthread_tracer_gather_step;

    master->dp_mode = DP_PTHREADS;

  } else
#else
  if (strcasecmp(master->params->dp_mode, "pthreads") == 0) {
    hd_quit("DP_MODE of PTHREADS has been specified but this executable has not been built with PTHREADS support. Please reconfigure, rebuild and try again.");
  } else
#endif // PTHREADS
#if defined(HAVE_MPI)
  if (strcasecmp(master->params->dp_mode, "mpi") == 0) {
    dp.init    = dp_mpi_init;
    dp.cleanup = dp_mpi_cleanup;
    /* All gathers point to none */
    dp.vel2d_step_p1        = dp_mpi_vel2d_step_p1;
    dp.vel2d_gather_step_p1 = dp_mpi_vel2d_gather_step_p1;
    dp.vel2d_step_p2        = dp_mpi_vel2d_step_p2;
    dp.vel2d_gather_step_p2 = dp_mpi_vel2d_gather_step_p2;
    dp.vel2d_step_p3        = dp_mpi_vel2d_step_p3;
    dp.vel2d_gather_step_p3 = dp_mpi_vel2d_gather_step_p3;
    dp.vel3d_step_p1        = dp_mpi_vel3d_step_p1;
    dp.vel3d_gather_step_p1 = dp_mpi_vel3d_gather_step_p1;
    dp.vel3d_step_p2        = dp_mpi_vel3d_step_p2;
    dp.vel3d_gather_step_p2 = dp_mpi_vel3d_gather_step_p2;
    dp.vel3d_post_p1        = dp_mpi_vel3d_post_p1;
    dp.vel3d_gather_post_p1 = dp_mpi_vel3d_gather_post_p1;
    dp.vel3d_post_p2        = dp_mpi_vel3d_post_p2;
    dp.vel3d_gather_post_p2 = dp_mpi_vel3d_gather_post_p2;
    dp.tracer_step          = dp_mpi_tracer_step;
    dp.tracer_gather_step   = dp_mpi_tracer_gather_step;

    /* Buffered mode 2 is mandatory for MPI */
    master->dp_mode = (DP_MPI|DP_BUFFEREDv2);

  } else
#else
  if (strcasecmp(master->params->dp_mode, "mpi") == 0) {
    hd_quit("DP_MODE of MPI has been specified but this executable has not been built with MPI support. Please reconfigure, rebuild and try again.");
  } else    
#endif // MPI
#if defined(HAVE_OMP)
  if (strncasecmp(master->params->dp_mode, "openmp", 6) == 0) {
    dp.init = dp_none_init;
    dp.cleanup = dp_none_cleanup;
    dp.vel2d_step_p1 = dp_none_vel2d_step_p1;
    dp.vel2d_gather_step_p1 = dp_none_vel2d_gather_step_p1;
    dp.vel2d_step_p2 = dp_none_vel2d_step_p2;
    dp.vel2d_gather_step_p2 = dp_none_vel2d_gather_step_p2;
    dp.vel2d_step_p3 = dp_none_vel2d_step_p3;
    dp.vel2d_gather_step_p3 = dp_none_vel2d_gather_step_p3;
    dp.vel3d_step_p1 = dp_none_vel3d_step_p1;
    dp.vel3d_gather_step_p1 = dp_none_vel3d_gather_step_p1;
    dp.vel3d_step_p2 = dp_none_vel3d_step_p2;
    dp.vel3d_gather_step_p2 = dp_none_vel3d_gather_step_p2;
    dp.vel3d_post_p1 = dp_none_vel3d_post_p1;
    dp.vel3d_gather_post_p1 = dp_none_vel3d_gather_post_p1;
    dp.vel3d_post_p2 = dp_none_vel3d_post_p2;
    dp.vel3d_gather_post_p2 = dp_none_vel3d_gather_post_p2;
    dp.tracer_step = dp_none_tracer_step;
    dp.tracer_gather_step = dp_none_tracer_gather_step;
    /*
     * Note: All of the above functions are indentical to the NONE
     *       case. There are additional directives in the gateway
     *       functions below to account for OpenMP
     */
    omp_set_num_threads(nwindows);

    master->dp_mode = DP_OPENMP;

    /* Check buffered mode */
    if (strcasecmp(master->params->dp_mode, "openmp_bufv1") == 0)
      master->dp_mode |= DP_BUFFEREDv1;
    else
      if (strcasecmp(master->params->dp_mode, "openmp_bufv2") == 0)
	master->dp_mode |= DP_BUFFEREDv2;

  } else
#else
  if (strcasecmp(master->params->dp_mode, "openmp") == 0) {
    hd_quit("DP_MODE of OPENMP has been specified but this executable has not been built with OpenMP support. Please reconfigure, rebuild and try again.");
  } else
#endif // OMP
  {
#if defined(HAVE_OMP)
    /* 
     * If we're here, then the executable has been built with OpenMP
     * support but called with the NONE DP_MODE flag in the prm file,
     * in which case, we setup for sequential mode, unless in TRANSPORT mode
     */
    if (master->runmode & TRANS) {
      omp_set_num_threads(master->params->trans_num_omp);
      omp_set_nested(1);
    } else {
      /* Fully coupled */
      omp_set_num_threads(1);
    }
#endif // OMP

    /* 
     * NONE mode
     *
     * No distributed processing and not using buffered transfer mode 
     */
    master->dp_mode = DP_NONE;

    dp.init = dp_none_init;
    dp.cleanup = dp_none_cleanup;
    dp.vel2d_step_p1 = dp_none_vel2d_step_p1;
    dp.vel2d_gather_step_p1 = dp_none_vel2d_gather_step_p1;
    dp.vel2d_step_p2 = dp_none_vel2d_step_p2;
    dp.vel2d_gather_step_p2 = dp_none_vel2d_gather_step_p2;
    dp.vel2d_step_p3 = dp_none_vel2d_step_p3;
    dp.vel2d_gather_step_p3 = dp_none_vel2d_gather_step_p3;
    dp.vel3d_step_p1 = dp_none_vel3d_step_p1;
    dp.vel3d_gather_step_p1 = dp_none_vel3d_gather_step_p1;
    dp.vel3d_step_p2 = dp_none_vel3d_step_p2;
    dp.vel3d_gather_step_p2 = dp_none_vel3d_gather_step_p2;
    dp.vel3d_post_p1 = dp_none_vel3d_post_p1;
    dp.vel3d_gather_post_p1 = dp_none_vel3d_gather_post_p1;
    dp.vel3d_post_p2 = dp_none_vel3d_post_p2;
    dp.vel3d_gather_post_p2 = dp_none_vel3d_gather_post_p2;
    dp.tracer_step = dp_none_tracer_step;
    dp.tracer_gather_step = dp_none_tracer_gather_step;
  }

  /* Setup each window and initialise. */
  dp.nwindows = nwindows;
  dp.dp_windows =
    (dp_window_t *)malloc(sizeof(dp_window_t) * (nwindows + 1));
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = wincon[1]->twin[nn];
    dp.dp_windows[n].window_id = n;
    dp.dp_windows[n].master = master;
    dp.dp_windows[n].geom = geom[n];
    dp.dp_windows[n].windata = windata[n];
    dp.dp_windows[n].wincon = wincon[n];
    dp.init(&dp.dp_windows[n]);
  }
  /* SCHED_MODE */
  master->thIO = master->params->thIO;
}

double dp_clock(void)
{
  struct tms tms;
  (void)times(&tms);
  return ((double)(tms.tms_utime + tms.tms_cutime)) / sysconf(_SC_CLK_TCK);
}

/*
 * These are all the gateway functions that SHOC calls. The
 * appropriate functions are installed in the dp struct which are
 * dispatched below
 */

void dp_vel2d_step_p1(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_step_p1(&dp.dp_windows[n]);
  }

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_gather_step_p1(&dp.dp_windows[n]);
  }
}

void dp_vel2d_step_p2(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_step_p2(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_gather_step_p2(&dp.dp_windows[n]);
  }
}

void dp_vel2d_step_p3(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_step_p3(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel2d_gather_step_p3(&dp.dp_windows[n]);
  }
}

void dp_vel3d_step_p1(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_step_p1(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_gather_step_p1(&dp.dp_windows[n]);
  }
}

void dp_vel3d_step_p2(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_step_p2(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_gather_step_p2(&dp.dp_windows[n]);
  }
}

void dp_vel3d_post_p1(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_post_p1(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_gather_post_p1(&dp.dp_windows[n]);
  }
}

void dp_vel3d_post_p2(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_post_p2(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.vel3d_gather_post_p2(&dp.dp_windows[n]);
  }
}

void dp_tracer_step(void)
{
  int nn, n;

#if defined(HAVE_OMP)
#pragma omp parallel for private(nn, n)
#endif
  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.tracer_step(&dp.dp_windows[n]);
  }

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.tracer_gather_step(&dp.dp_windows[n]);
  }
}


/* Cleanup distributed processing */
void dp_cleanup(void)
{
  int nn, n;

  for (nn = 1; nn <= dp.nwindows; nn++) {
    n = dp.dp_windows[1].wincon->twin[nn];
    dp.cleanup(&dp.dp_windows[n]);
  }

  free(dp.dp_windows);
}
