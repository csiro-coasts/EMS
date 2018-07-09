/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/ems.h
 *
 *  \brief Top-level EMS header
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: ems.h 5871 2018-07-06 07:09:44Z riz008 $
 */

#ifndef	_EMS_H
#define _EMS_H

/* For Window 32 if required */
#ifdef _WIN32
int strcasecmp(const char *s1, const char *s2);
int strncasecmp(const char *s1, const char *s2, int n);
#endif

/* Maximum length of an input line */
#define MAXLINELEN (40000)  /* NOTE: check stack size if increasing this value */

/* Maximum string length */
#define MAXSTRLEN (2048)

/* Maximum length of a file name */
#define MAXFNAMELEN (200)

/* define an empty String */
#define EMPTY ""

#include <stdio.h>
#include <unistd.h>
#include "ems_conf.h"
#include "emsmath.h"
#include "integrator.h"
#include "emsalloc.h"
#include "string_utils.h"
/*UR work in progress #include "num_utils.h" */
#include "prmfile.h"
#include "errfn.h"
#include "poly.h"
#include "colourtable.h"
#include "mapproj.h"
#include "datafile.h"
#include "timeseries.h"
#include "grid.h"
#include "hash.h"
#include "hqueue.h"
#include "ptrack.h"
#include "ncw.h"
#include "nan.h"
#include "topo.h"
#include "coast.h"
#include "tracer.h"
#include "time_utils.h"
#include "emslogger.h"
#include "solar.h"
#include "interp.h"
#include "stringtable.h"
#include "cstmesh.h"
/*
  #include "dyn_loading.h"
*/
#include "underwater.h"

#include "svn_rev.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#if !defined(INLINE)
#if defined(__GNUC__) || defined(__cplusplus__)
#define INLINE inline
#else
#define INLINE
#endif
#endif

/* Release verions and getters */
#define EMSLIB_MAJOR_VERSION 1
#define EMSLIB_MINOR_VERSION 0
#define EMSLIB_PATCH_VERSION 0

int get_emslib_major_vers(void);
int get_emslib_minor_vers(void);
int get_emslib_patch_vers(void);

/* Prototypes */

void    spline(double *x, double *y, long n, int derivspec,
            double start_deriv, double end_deriv, double *ydd);
void    spline_interp(double *xa, double *ya, double *ydd, long n, double x,
                   double *y);
void    contour(double **val, double **x, double **y, long nx, long ny,
             char *cval, double badval);
void    cfft(double data[], int ndata, int dirn);
double  tm_to_juldays(int y, int mo, int d, int h, int mi, int s);
void    tm_to_julsecs(double j, int *y, int *mo, int *d, int *h, int *mi,
                   int *s);
int     tm_scale_to_secs(char *str, double *sec);
char*   tm_time_to_datestr(double t, char *u);
double  tm_datestr_to_julsecs(char *d, char *u);
double  tm_time_to_julsecs(char *d);
void    tm_change_time_units(char *oepoch, char *nepoch, double *times,
                            int n);
double  tm_tz_offset(char *epoch);

double clock_step();
/* ADDED UR 3/11/2004
 */
double  tm_unit_to_sec(const char* units);
char*   tm_units_format( char* units);

/*ADDED UR 08/2010 */
char* tm_extract_unit(char *epoch);
double tm_extract_sec_scale(char *epoch);

int tm_time_to_ymd(double t, char *tunit, int *yr, int *mon, int *day);

int     overwrite_output();
/*UR work in progress
int     large_output();
*/
/* level of logging statements to stdout
 * (will change to  some general logging in the future)
 */
int     verbosity();

void    ems_init(FILE* f, int skip);
void    ems_clear();
void    ems_usage();
/* UR */
#endif                          /* _EMS_H */
