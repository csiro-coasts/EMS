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
 *  $Id: globals.h 7231 2022-10-25 00:23:13Z her127 $
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
extern ts_point_t tserror;


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

/*-------------------------------------------------------------------*/
/* Valid SST import URLs                                             */
static char *sst_url[7][7] = {
  {"GHRSST_L4a", "GHRSST L4 SST", "1 day", "https://data.nodc.noaa.gov/thredds/dodsC/ghrsst/L4/GLOB/UKMO/OSTIA/", "-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2", "120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc", NULL},
  {"GHRSST_L4b", "GHRSST L4 SST", "1 day", "https://www.ncei.noaa.gov/thredds-ocean/dodsC/ghrsst/L4/GLOB/UKMO/OSTIA/", "-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2", "120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc", NULL},
  {"BOM_L3", "BoM L3S SST", "1 day", "http://rs-data1-mel.csiro.au/thredds/dodsC/imos-srs/sst/ghrsst/L3S-6d/ngt/", "032000-ABOM-L3S_GHRSST-SSTskin-AVHRR_D-6d_night-v02.0-fv02.0.nc", NULL, "NODAY VARIABLES (ghrsst=sea_surface_temperature)(ghrsst_error=quality_level)"},
  {"HIM_L3a", "HIMAWARI L3C SST", "1 hour", "https://www.star.nesdis.noaa.gov/thredds/dodsC/gridH08AHINRTL3CWW00/", "00-STAR-L3C_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc", NULL, "VARIABLES (ghrsst=sea_surface_temperature)(ghrsst_error=quality_level)"},
  {"HIM_L3b", "HIMAWARI L3C SST", "1 hour", "https://www.ncei.noaa.gov/thredds-ocean/dodsC/ghrsst/L3C/PACIFIC/H08/STAR/", "00-STAR-L3C_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.70-v02.0-fv01.0.nc", NULL, "VARIABLES (ghrsst=sea_surface_temperature)(ghrsst_error=quality_level)"},
  {"HIM_L3c", "HIMAWARI BOM L3C SST", "1 hour", "https://dapds00.nci.org.au/thredds/dodsC/qm43/ghrsst/v02.0fv02/Continental/L3C-01hour/ABOM-L3C_GHRSST-SSTskin-AHI_H08/", "00-ABOM-L3C_GHRSST-SSTskin-AHI_H08-1h-v02.0-fv02.0.nc", NULL, "DOMON VARIABLES (ghrsst=sea_surface_temperature)(ghrsst_error=quality_level)"},
  {NULL, NULL, NULL, NULL, NULL, NULL, NULL}
};

#endif

#endif                          /* _HD_GLOBALS_H */
