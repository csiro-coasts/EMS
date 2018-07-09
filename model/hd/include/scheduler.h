/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/include/scheduler.h
 *  
 *  Description:
 *  Schedules when model events should be dispatched.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: scheduler.h 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#ifndef _SCHEDULER_H
#define _SCHEDULER_H

#define SEPS 1e-5               /* Allowable error in time scheduling */

/**** Structures, functions, etc.
 ****/
typedef struct sched_event sched_event_t;


typedef int (*SchedInitFunc) (sched_event_t *);
typedef double (*SchedEventFunc) (sched_event_t *, double t);
typedef void (*SchedCleanupFunc) (sched_event_t *, double t);
typedef void* (*SchedDispatchFunc) (void* data);


/* typedef void (*SchedDispatchFunc) (sched_event_t* dispatch, SchedInitFunc in_progress); */

struct sched_event {
  char name[MAXSTRLEN];         /* Name of the event */
  sched_event_t *next;          /* Next link in sched_event_t linked list */
  SchedEventFunc event;         /* Function called when event time reached
                                 */
  SchedCleanupFunc cleanup;     /* Cleaup function */
  SchedDispatchFunc dispatch;   /* trigger asynchroniuous workload, data maybe NULL*/
  SchedInitFunc    in_progress; /* verify that a dispatched function has returned */
  void *public_data;            /* Public data used by event functions */
  void *private_data;           /* Private data used in event functions */
  double next_event;            /* Time of next event */

#ifdef HAVE_PTHREADS
  /* 
   * For synchronisation of threaded dispatches - currently only dumps
   */
  sem_t *sm_fill;  // buffer is filled
  sem_t *sm_dump;  // buffer has dumped
#endif
};


typedef struct {
  sched_event_t *head;          /* Head of the sched_event_t linked list */
  double t;                     /* Current time in units of ts_units */
  double start_time;            /* Start time in units of ts_units */
  double stop_time;             /* Finish time in units of ts_units */
  char units[MAXSTRLEN];        /* ISO units for time */
  time_t exec_start_time;       /* Execution start time */
  SchedEventFunc dispatch;      /* function to dispatch an events dispatch function */
} scheduler_t;



/**** Prototypes for time scheduler functions.
 ****/

/* Create and destroy a time schedular structure.
 * The prmfd variable is the FILE id of the parameter file
 * from which the TIME UNITS, MODEL TIME, etc will be read.
 */
extern scheduler_t *sched_init(FILE * prmfd, time_t exec_start_time);
extern void sched_end(scheduler_t *sched);


/* Functions to register/deregister an interest in an event occuring.
 * The init and event functions return the time when next this
 * 'Event' should be serviced.
 */
extern void sched_register(scheduler_t *sched, char *name,
                           SchedInitFunc init, SchedEventFunc event,
                           SchedCleanupFunc cleanup, void *data,
                           SchedDispatchFunc dispatch, SchedInitFunc in_progress);
extern void sched_deregister(scheduler_t *sched, char *name);

/* Process the events list to determine at which time the
 * next stop should occur. If that time has already been
 * reached, then call the event procedure for each pending
 * event.
 */
extern double sched_get_next_event(scheduler_t *sched);
extern void sched_set_time(scheduler_t *sched, double t);
extern sched_event_t *sched_get_even_by_name(scheduler_t *sched,
                                             char *name);

/* Macro functions. */
#define schedGetName(event) (event->name)
#define schedGetNext(event) (event->next_event)
#define schedSetPublicData(event, data) (event->public_data = data)
#define schedGetPublicData(event) (event->public_data)
#define schedSetPrivateData(event, data) (event->private_data = data)
#define schedGetPrivateData(event) (event->private_data)


#endif                          /* _SCHEDULER_H */
