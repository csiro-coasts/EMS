/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/control/globals.c
 *  
 *  Description:
 *  Definitions of global variables for
 *  3-d non-linear hydrodynamic model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: globals.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#define IN_DECLS
#include <stdio.h>
#include "hd.h"


/* Time variables */
scheduler_t *schedule;          /* Time scheduler */


/* PRIVATE: timeseries.c */
sched_event_t *ts_events;       /* Time scheduler_t timeseries events */
ts_point_t *tslist;             /* time series point list */
int nts;                        /* number of time series points */

/* descriptions */
char executable[MAXSTRLEN];     /*UR executable filename */
char codeheader[MAXSTRLEN];     /* grid type description */
char parameterheader[MAXSTRLEN];  /* parameter description */
char prmname[MAXSTRLEN];        /* parameter file name */
char oname[MAXSTRLEN];          /* output dumpfile name */

int killed = 0;                 /* Flag to indicate if the model should be 
                                   stopped */
int model_running = 0;          /* Is the model running */
double air_dens;                /* density of air */
double ambpress;                /* Ambient atmospheric pressure */
double g;                       /* gravity */
double spec_heat;               /* Specific heat of water */

hd_data_t *hd_data;
geometry_t *geom;
geometry_t **window;
open_bdrys_t **Opens;
window_t **windat;
win_priv_t **wincon;
master_t *master;
parameters_t *params;
dump_data_t *dumpdata;
ts_point_t tsflush;
ts_point_t tsphist;
ts_point_t tsalert;
ts_point_t tstotal;
ts_point_t tstrans;

int debug = 0;                  /* enable/disable debug */
int warnings = 1;               /* enable/disable warning message */
int diag_log = 1;               /* enable/disable diagnostic logging */
char diag_logfile[MAXSTRLEN];
int windows_log = 0;            /* enable/disable warning message */
char window_geom_logfile[MAXSTRLEN];
char window_map_logfile[MAXSTRLEN];
int setup_log = 1;              /* enable/disable warning message */
char setup_logfile[MAXSTRLEN];
char projection[MAXSTRLEN];     /* Projection information */
int autof = 0;                  /* automated initialization */
int forced_restart = 0;         /* force the model restart */
int nrt_restart = 0;            /* Near real-time model restart */
int one_step = 0;               /* Run for exactly one timestep */
int crash_restart = 0;          /* Crash recovery restarts */
char version[MAXSTRLEN];        /* version string of this executable */
/*UR-ADDED */
int stat_log = 0;               /* enable/disable raw runtime statistics */
char stat_logfile[MAXSTRLEN];

/*
 * The following are related to the timing profile
 */
#ifdef DO_TIMING
char timing_logfile[MAXSTRLEN];
FILE *tfp = NULL;
long tfp_pos = 0;
int timing_counter = 1;
int timing_level = TIMING_LEVEL_DEF;
#endif // DO_TIMING

// EOF

