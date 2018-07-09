/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/control/dp_pthreads.c
 *  
 *  Description:
 *  Manage distributed processing using pthreads.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dp_pthreads.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "hd.h"

#if defined(HAVE_PTHREADS)
#include <pthread.h>
#include <semaphore.h>


/*UR this is necessary to run -pg at compile time for multiple threads
 * otherwise the model run will fail.
 */
#ifdef HAVE_GPROF
#include <sys/time.h>
#endif



typedef enum {DP_CMD_VEL2D_STEP_P1=1,
  DP_CMD_VEL2D_STEP_P3,
  DP_CMD_VEL2D_STEP_P2,
  DP_CMD_VEL3D_STEP_P1,
  DP_CMD_VEL3D_STEP_P2,
  DP_CMD_VEL3D_POST_P1,
  DP_CMD_VEL3D_POST_P2,
  DP_CMD_TRACER_STEP,
  DP_CMD_EXIT = -1
} dp_pthread_cmd_t;

typedef struct {
  pthread_t thread;
  pthread_attr_t attr;
  pthread_mutex_t mutex;
  sem_t sem_request;
  sem_t sem_response;
  dp_pthread_cmd_t cmd;
#ifdef HAVE_GPROF
  struct itimerval itimer;
#endif
} dp_pthread_t;


static void *dp_window_thread(void *p)
{
  dp_window_t *dpw = (dp_window_t *)p;
  dp_pthread_t *td = (dp_pthread_t *)dpw->private_data;
  int finished = 0;

  while (!finished) {

    dp_pthread_cmd_t cmd;
    /*printf("Thread %ld - Waiting for request from %i\n", pthread_self(),td);*/
    /* Wait for request */
    sem_wait(&td->sem_request);

    /*printf("Thread %ld - Received request from %i\n", pthread_self(),td); */

    /* Process request */
    pthread_mutex_lock(&td->mutex);
    cmd = td->cmd;
    pthread_mutex_unlock(&td->mutex);

    /*printf("Thread %ld - Processing request %d from %i\n", pthread_self(),cmd,td); */


    switch (cmd) {

    case DP_CMD_VEL2D_STEP_P1:
      dp_none_vel2d_step_p1(dpw);
      break;

    case DP_CMD_VEL2D_STEP_P2:
      dp_none_vel2d_step_p2(dpw);
      break;

    case DP_CMD_VEL2D_STEP_P3:
      dp_none_vel2d_step_p3(dpw);
      break;

    case DP_CMD_VEL3D_STEP_P1:
      dp_none_vel3d_step_p1(dpw);
      break;

    case DP_CMD_VEL3D_STEP_P2:
      dp_none_vel3d_step_p2(dpw);
      break;

    case DP_CMD_VEL3D_POST_P1:
      dp_none_vel3d_post_p1(dpw);
      break;

    case DP_CMD_VEL3D_POST_P2:
      dp_none_vel3d_post_p2(dpw);
      break;

    case DP_CMD_TRACER_STEP:
      dp_none_tracer_step(dpw);
      break;

    case DP_CMD_EXIT:
      finished = 1;
      break;
    }

    /* Acknowledge completion */
    sem_post(&td->sem_response);
    /*printf("Thread %ld - Finished request %d from %i\n", pthread_self(),cmd,td); */
  }

  pthread_exit((void *)NULL);

  return (void*)NULL;
}

static void dp_pthread_send_cmd(dp_window_t *dpw, dp_pthread_cmd_t cmd)
{
  dp_pthread_t *td = (dp_pthread_t *)dpw->private_data;

#ifdef HAVE_GPROF
  setitimer(ITIMER_PROF, &(td->itimer), NULL);
#endif

  pthread_mutex_lock(&td->mutex);
  td->cmd = cmd;
  pthread_mutex_unlock(&td->mutex);

  /* release the qeue to process the command */
  sem_post(&td->sem_request);
}

static void dp_pthread_gather(dp_window_t *dpw)
{
  dp_pthread_t *td = (dp_pthread_t *)dpw->private_data;
/*printf("Thread %ld - wait for response %i \n", pthread_self(),td); */
  /* wait for the qeue to finish the command */
  sem_wait(&td->sem_response);
/* printf("Thread %ld - past response %i \n", pthread_self(),td);*/
}


/* Protected PTHREAD methods */
void dp_pthread_init(dp_window_t *dpw)
{
  dp_pthread_t *td = (dp_pthread_t *)malloc(sizeof(dp_pthread_t));
  pthread_attr_init(&td->attr);
  pthread_attr_setdetachstate(&td->attr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&td->attr, SCHED_OTHER);
  pthread_mutex_init(&td->mutex, NULL);
  sem_init(&td->sem_request, 0,0);
  sem_init(&td->sem_response, 0,0);
#ifdef HAVE_GPROF
  getitimer(ITIMER_PROF, &(td->itimer));
#endif
  dpw->private_data = td;
  pthread_create(&td->thread, &td->attr, dp_window_thread, dpw);
}

void dp_pthread_cleanup(dp_window_t *dpw)
{
  dp_pthread_t *td = (dp_pthread_t *)dpw->private_data;

  dp_pthread_send_cmd(dpw, DP_CMD_EXIT);
  dp_pthread_gather(dpw);

  sem_destroy(&td->sem_response);
  sem_destroy(&td->sem_request);
  pthread_mutex_destroy(&td->mutex);
  free(td);
}

void dp_pthread_vel2d_step_p1(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL2D_STEP_P1);
}

void dp_pthread_vel2d_step_p2(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL2D_STEP_P2);
}

void dp_pthread_vel2d_step_p3(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL2D_STEP_P3);
}

void dp_pthread_vel2d_gather_step_p1(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel2d_gather_step_p2(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel2d_gather_step_p3(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel3d_step_p1(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL3D_STEP_P1);
}

void dp_pthread_vel3d_gather_step_p1(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel3d_step_p2(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL3D_STEP_P2);
}

void dp_pthread_vel3d_gather_step_p2(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel3d_post_p1(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL3D_POST_P1);
}

void dp_pthread_vel3d_gather_post_p1(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_vel3d_post_p2(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_VEL3D_POST_P2);
}

void dp_pthread_vel3d_gather_post_p2(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

void dp_pthread_tracer_step(dp_window_t *dpw)
{
  dp_pthread_send_cmd(dpw, DP_CMD_TRACER_STEP);
}

void dp_pthread_tracer_gather_step(dp_window_t *dpw)
{
  dp_pthread_gather(dpw);
}

#endif
