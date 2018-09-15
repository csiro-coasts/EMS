/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/globals.h
 *  
 *  Description:
 *  Global variables.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: globals.h 5915 2018-09-05 03:30:40Z riz008 $
 *
 */

#ifndef _HD_GLOBALS_H
#define _HD_GLOBALS_H

#ifndef IN_DECLS


/* Time variables */
extern scheduler_t *schedule;   /* scheduler_t */

/* version string */
extern char version[];           /* version identifier */

/* descriptions */
extern char executable[MAXSTRLEN];     /*UR executable filename */
extern char codeheader[MAXSTRLEN];  /* grid type description */
extern char parameterheader[MAXSTRLEN]; /* parameter description */
extern char prmname[MAXSTRLEN]; /* parameter file name */
extern char oname[MAXSTRLEN];   /* output dumpfile name */

extern ts_point_t *tslist;      /* time series point list */
extern int nts;                 /* number of time series points */
extern sched_event_t *ts_events;  /* Time scheduler_t timeseries events */
extern ts_point_t tsflush;
extern ts_point_t tsphist;
extern ts_point_t tsalert;
extern ts_point_t tstotal;
extern ts_point_t tstrans;


/* Fundamental values */
extern double g;                /* gravity */
extern double air_dens;         /* density of air */
extern double ambpress;         /* Ambient atmospheric pressure */
extern double spec_heat;        /* specific heat of water */


extern int killed;
extern int model_running;
extern int debug;               /* enable/disable debug messages */
extern int warnings;            /* enable/disable warning messages */
extern int diag_log;            /* enable/disable warning message */
extern char diag_logfile[MAXSTRLEN];
extern int windows_log;         /* enable/disable warning message */
extern char window_geom_logfile[MAXSTRLEN];
extern char window_map_logfile[MAXSTRLEN];
extern int setup_log;           /* enable/disable warning message */
extern char setup_logfile[MAXSTRLEN];
/*UR-ADDED */
extern int stat_log;            /* enable/disable raw runtime statistics */
extern char stat_logfile[MAXSTRLEN];

extern hd_data_t *hd_data;
extern geometry_t *geom;
extern geometry_t **window;
extern open_bdrys_t **Opens;
extern window_t **windat;
extern win_priv_t **wincon;
extern master_t *master;
extern parameters_t *params;
extern dump_data_t *dumpdata;


extern char projection[MAXSTRLEN];  /* Projection information */

extern int forced_restart;         /* force the model restart */
extern int nrt_restart;            /* near real-time model restart */
extern int one_step;               /* one model time step */
extern int crash_restart;          /* restart when model crashes */

extern int threaded;
extern int autof;               /* automated initialization */
extern int wnnum;
extern FILE *op;

extern int mpi_rank;
extern int mpi_size;

/*
 * The following are related to the timing profile
 */
extern char timing_logfile[];
extern FILE *tfp;
extern long tfp_pos;
extern int timing_counter;
extern int timing_level;

#endif

#endif                          /* _HD_GLOBALS_H */
