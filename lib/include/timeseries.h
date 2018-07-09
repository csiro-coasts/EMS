/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/timeseries.h
 *
 *  \brief Include file for routines which deal with time series data
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: timeseries.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _TIMESERIES_H
#define _TIMESERIES_H

#include <stdio.h>
#include "datafile.h"

/* Time series types */
#define        TS_ASCII     0   /* Ascii file */
#define        TS_NC0       1   /* netCDF with 0 spatial dimensions */
#define        TS_NC1       2   /* netCDF with 1 spatial dimension */
#define        TS_NC2       3   /* netCDF with 2 spatial dimensions */
#define        TS_NC3       4   /* netCDF with 3 spatial dimensions */

/* Maximum number of ts files allowed on a single line in the prm file */
#define MAX_TS_FILES (20)

/* Time series modulus */
#define MOD_NONE	0
#define MOD_YEAR	1
#define MOD_MONTH	2
#define MOD_WEEK	3
#define MOD_DAY		4
#define MOD_HOUR	5
#define MOD_MINUTE	6
#define MOD_SECOND	7

/* In converting timeseries to datafile, only a subset of the more
 * commonly accessed structure variables have been preserved.
 * this facilitates some backwards compatibility */
typedef struct {
  int type;                     /* Timeseries type */
  char *name;                   /* Time series name */
  long nt;                      /* Number of time values */
  int ti;                       /* Index of time variable */
  char *t_units;                /* Time units string */
  double *t;                    /* time values */
  int nv;                       /* Number of variables */
  char **varname;               /* df_variable_t names */
  char **varunit;               /* df_variable_t units */
  int t_mod_type;               /* Modulo type */
  double t_mod_scale;           /* Modulo scale */
  double t_unit_scalef;
  double t_base;
  datafile_t *df;               /* Pointer to the data file */
} timeseries_t;


/* Prototypes */
void ts_read(char *name, timeseries_t *ts);
void ts_set_proj_type(timeseries_t *ts, char *type);
void ts_set_default_proj_type(char *type);
void ts_set_default_hashtable_size(int hsize);
double ts_eval(timeseries_t *ts, int varid, double t);
double ts_eval_xy(timeseries_t *ts, int varid, double t, double x,
                  double y);
double ts_eval_xyz(timeseries_t *ts, int varid, double t, double x,
                   double y, double z);
int ts_eval_xy_flag(timeseries_t *ts, int id, double t, double x, double y,
		    double *out);
int ts_multifile_check(int ntsfiles, timeseries_t *tsfiles, char *var,
                       double tstart, double tstop);
double ts_multifile_eval_xyz(int ntsfiles, timeseries_t *tsfiles,
                             char *var, double t, double x, double y,
                             double z);
int ts_get_index(timeseries_t *ts, char *varname);
void ts_convert_time_units(timeseries_t *ts, char *newunits);
int ts_has_time(timeseries_t *ts, double t);
void ts_flush_variable(timeseries_t *ts, int varid);
void ts_free(timeseries_t *ts);
void ts_is_time_monotonic(timeseries_t *ts);
void ts_print_info(timeseries_t *ts, FILE * fp);
double get_file_time(timeseries_t *ts, double t);
void ts_write_ascii(FILE * fp, timeseries_t *ts);
double ts_get_att_value(timeseries_t *ts, char *name, char *attname);

int ts_multifile_read(char *fnames, timeseries_t *ts);
int ts_multifile_get_index(int nfiles, timeseries_t *ts, char *varname);
double ts_multifile_eval(int nfiles,timeseries_t *ts,int multi_varid, double t);
int ts_multifile_get_nv(int nfiles, timeseries_t *ts);
int ts_eval_runcode(timeseries_t *ts);
const char* ts_get_varname(timeseries_t *ts, int vid);
int ts_is_modulo(timeseries_t *ts);
int ts_is_recom(timeseries_t *ts);
int ts_var_z_is_depth(timeseries_t *ts, int id);

#endif                          /* _TIMESERIES_H */
