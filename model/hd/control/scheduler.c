/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/control/scheduler.c
 *  
 *  Description:
 *  Schedules when model events should be
 *  dispatched.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: scheduler.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include "hd.h"


#ifdef HAVE_PTHREADS

#include <pthread.h>
#include <time.h>
#include <unistd.h>

typedef struct {
  pthread_t thread;
  pthread_attr_t attr;
#ifdef HAVE_GPROF
  struct itimerval itimer;
#endif
} sched_pthread_t;

double schedDPdispatch(sched_event_t* dispatch, double t)
{
  sched_pthread_t *td = (sched_pthread_t *)malloc(sizeof(sched_pthread_t));
  pthread_attr_init(&td->attr);
  pthread_attr_setschedpolicy(&td->attr, SCHED_OTHER);

#ifdef HAVE_GPROF
  getitimer(ITIMER_PROF, &(td->itimer));
#endif

  emstag(LDEBUG,"hd:scheduler:schedDPdispatch","dispatching event '%s', at %1.1d until %1.1d (sec)",dispatch->name,t,schedule->stop_time);

  /* I'm not sure if this if statement is ever hit or is needed */
  if(t >= schedule->stop_time) {
    emstag(LWARN,"hd:scheduler:schedDPdispatch","Last call, joining thread...");
    pthread_attr_setdetachstate(&td->attr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&td->thread, &td->attr, dispatch->dispatch, dispatch);
    /* we are finished, join the thread to wait */
    pthread_join(td->thread,NULL);
  } else { 
    /* normal case */
    pthread_attr_setdetachstate(&td->attr, PTHREAD_CREATE_DETACHED);
    pthread_create(&td->thread, &td->attr, dispatch->dispatch, dispatch);
  }
  return 0;
}

#endif


double schedDispatch(sched_event_t* dispatch, double t)
{
  dispatch->dispatch(dispatch);
  return 0;
}


/* Create a time schedular structure.
 * The prmfd variable is the FILE id of the parameter file
 * from which the TIME UNITS, MODEL TIME, etc will be read.
 */
scheduler_t *sched_init(FILE * fp,  time_t exec_start_time)
{
  scheduler_t *sched = (scheduler_t *)malloc(sizeof(scheduler_t));
  char buf[MAXSTRLEN];
  sched->head = NULL;

/* Read the time scheduler parameters from the Model parameter file.
 */
  prm_set_errfn(hd_silent_warn);
  if (prm_read_char(fp, "TIMEUNIT", sched->units) == 0) {
    prm_set_errfn(hd_quit);
    prm_read_char(fp, "INTERNAL_TIMEUNIT", sched->units);
  }
  prm_set_errfn(hd_quit);
  if (forced_restart) {
    if (!(prm_read_char(fp, "restart_name", buf)))
      strcpy(buf, "restart.nc");
     sched->start_time = get_restart_time(buf, sched->units);
     prm_get_time_in_secs(fp, "STOP_TIME", &sched->stop_time);
  } else if (nrt_restart) {
    prm_read_char(fp, "INPUT_FILE", buf);
    sched->start_time = get_restart_time(buf, sched->units);
    prm_get_time_in_secs(fp, "STOP_TIME", &sched->stop_time);
    sched->stop_time += sched->start_time;
  } else {
     prm_get_time_in_secs(fp, "START_TIME", &sched->start_time);
     prm_get_time_in_secs(fp, "STOP_TIME", &sched->stop_time);
     if (one_step) {
       double dt;
       prm_get_time_in_secs(fp, "DT", &dt);
       sched->stop_time = sched->start_time + dt;
     }
  }
  /*UR moved to later to get a real time
  sched->exec_start_time = exec_start_time; */
  sched->t = sched->start_time;

 #ifdef HAVE_PTHREADS
  if (prm_read_char(fp, "SCHED_MODE", buf)) {
    if (strcasecmp(buf, "pthreads") == 0)
      sched->dispatch = schedDPdispatch;
    else
      sched->dispatch = schedDispatch;
  } else
#endif
    {
      sched->dispatch = schedDispatch;
    }

  return sched;
}


/* Cleanup a time scheduler structure, including calling the
 * event cleanup routines, and releasing memory.
 */
void sched_end(scheduler_t *sched)
{
#pragma omp master
  if (sched != NULL) {
    sched_event_t *cur = sched->head;
    sched_event_t *next;
    
    /* 
     * Cascade through the list of time scheduler events, cleaning up
     * and deallocating the memory as we go.
     */
    while (cur != NULL) {
      next = cur->next;
      if (cur->cleanup != NULL)
        cur->cleanup(cur, sched->t);
#ifdef HAVE_PTHREADS
      if (cur->sm_fill != NULL)	{
	sem_destroy(cur->sm_fill);
	free(cur->sm_fill);
      }
      if (cur->sm_dump != NULL)	{
	sem_destroy(cur->sm_dump);
	free(cur->sm_dump);
      }
#endif
      free(cur);
      cur = next;
    }
    
    free(sched);
    sched = NULL;
  }
}


/** Functions to register/deregister an interest in an event occuring.
 * The init and event functions return the time when next this
 * 'Event' should be serviced.
 *
 * @param sched - the main data object
 * @param name - the name under which the event can be accessed
 * @param init - the initialisation function, called once here
 * @param event - the function called every time this event is due
 * @param cleanup - the function called at the end of the run
 * @param data - event specific data, accessible via schedGetPublicData
 * @param dispatch - may be NULL, function called after the event has finished
 *                   the function is placed in a separate thread
 */
void sched_register(scheduler_t *sched,
                    char *name,
                    SchedInitFunc init, SchedEventFunc event,
                    SchedCleanupFunc cleanup, void *data,
                    SchedDispatchFunc dispatch, SchedInitFunc in_progress)
{
  sched_event_t *e = sched_get_even_by_name(sched, name);

  if (e != NULL) {
    hd_quit
      ("sched_register: Attempt to override the scheduled event '%s'\n",
       name);
    return;
  }

  e = (sched_event_t *)malloc(sizeof(sched_event_t));
  strcpy(e->name, name);
  e->next = sched->head;
  e->event = event;
  e->dispatch = dispatch;
  e->in_progress = in_progress;
  e->cleanup = cleanup;
  e->public_data = data;
  e->private_data = NULL;
  e->next_event = sched->start_time;
  e->sm_fill = NULL;
  e->sm_dump = NULL;


/* Initialise the routines, and get the next event.
 */
  if (init(e)) {
    e->next_event = event(e, sched->t);
    sched->head = e;
    if(e->dispatch != NULL)
    	e->dispatch(e);
    /* printf("%s initialised\n",e->name); */
  } else {
    /* printf("%s freed\n",e->name); */
    free(e);
  }
}


void sched_deregister(scheduler_t *sched, char *name)
{

  sched_event_t *event = sched_get_even_by_name(sched, name);

  if ((sched != NULL) && (event != NULL)) {
    sched_event_t *prev = NULL;
    sched_event_t *next = NULL;
    sched_event_t *cure = sched->head;

/* Search for the dump, and remove it from the list.
 */
    while (cure != NULL) {
      next = cure->next;
      if (cure == event) {
        if (prev == NULL)
          sched->head = next;
        else
          prev->next = next;
        break;
      }

      prev = cure;
      cure = next;
    }

/* Cleanup the event - force a dump, restore memory, etc.
 */
    if (event->cleanup != NULL)
      event->cleanup(event, sched->t);
    free(event);
  }
}


/* Calculates the time until the next event should occur.
 */
double sched_get_next_event(scheduler_t *sched)
{
  sched_event_t *next = sched->head;
  double earliest = sched->stop_time + 1.0; /* Infinity */

/* All of the pending events have been updated, now we
 * must check to see when the next event should go off.
 */
  while (next != NULL) {
    if ((next->event != NULL) && (next->next_event < earliest)) {
      earliest = next->next_event;
    }
    next = next->next;
  }
  return earliest;
}


void sched_set_time(scheduler_t *sched, double t)
{
  sched_event_t *next = sched->head;
  double time = sched->t;

  /* Crash recovery : reset the model to restart data */
  if (master->crf == RS_RESTART) {
    crash_event(sched_get_even_by_name(sched, "crash"), t);
    t = sched->t;
  }
  /* Check the events to see if the next_event time has been passed.
   * if so, then update the event.
   */
  while (next != NULL) {
    if (next->event != NULL) {
      if (t >= (next->next_event - DT_EPS)) {
	/*printf("Event %s %f\n",next->name,next->next_event);*/
	TIMING_SET;
	next->next_event = next->event(next, t);
	TIMING_PRINT(" ");
	TIMING_DUMP(1, next->name);
	/*printf("Schedule %s\n",next->name);*/
	if(next->dispatch != NULL) {
	  TIMING_SET;
	  sched->dispatch(next,t);
	  TIMING_DUMP(1, "dispatch");
	}
      }
    }
    next = next->next;
  }

  /*
   * If the model time is reset to a past time, then re-initialize
   * all the next_events to this past time, except data assimilation.
   */
  if (sched->t < time) {
    t = sched->t;
    next = sched->head;
    while (next != NULL) {
      if (strcmp(next->name, "dassim") != 0) {
	next->next_event = t;
      }
      next = next->next;
    }
  }

  sched->t = t;
}



/* Search the scheduler list to find an event with the name as specified.
 * Return the event if located or alternatively return NULL.
 */
sched_event_t *sched_get_even_by_name(scheduler_t *sched, char *name)
{
  sched_event_t *event = sched->head;

  while (event != NULL) {
    if (strcmp(event->name, name) == 0)
      return event;
    event = event->next;
  }

  return NULL;
}
