/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/inputs/choosedump.c
 *  
 *  Description:
 *  Allow the user to choose which dump
 *  in the input file is to be used to initialise
 *  the model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: choosedump.c 6608 2020-09-07 03:28:20Z her127 $
 *
 */

#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include "hd.h"

#define START_EPS (0.001)

int dump_choose(dump_data_t *dumpdata, int fid)
{
  unsigned int i;
  size_t n;
  size_t start[2];
  size_t count[2];
  double *tvals;
  int dnum;
  char timeunits[MAXSTRLEN];

  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &n);

  if (DEBUG("dump"))
    dlog("dump", "%d dumps in input file.\n", n);

  if (n < 1)
    hd_quit_and_dump("No dumps in input file!\n");
  if (n == 1)
    return (0L);

  if ((tvals = (double *)malloc(sizeof(double) * n)) == NULL)
    hd_quit_and_dump
      ("dump_choose: Can't allocate memory for time values\n");

  /* read the time values and let the user choose */
  start[0] = 0;
  count[0] = n;
  memset(timeunits, 0, MAXSTRLEN);
  nc_get_vara_double(fid, ncw_var_id(fid, "t"), start, count, tvals);
  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
  emstag(LTRACE,"hd:choosedump:dump_choose","Input file contains the following dumps:");
  for (i = 0; i < n; i++)
    emstag(LTRACE,"hd:choosedump:dump_choose", "%d    %10g %ss\n", i, tvals[i], timeunits);
  emstag(LTRACE,"hd:choosedump:dump_choose", "Enter dump number to load >");
  if ((scanf("%d", &dnum) != 1) || (dnum < 0) || (dnum >= (int)n))
    hd_quit_and_dump("dump_choose: Invalid dump number entered\n");

  free((char *)tvals);
  return (dnum);
}


/* This routine searches through the input dump file for a
 * record with the same time as the model start time specified
 * in the parameter file. The parameters are as follows:
 *
 * dumpdata:    The model grid
 * fid:      The netCDF id of the input dump file
 * t:        The model start time as specified in the parameter file
 *           (but converted to seconds).
 *          
 * The first dump which has a time within START_EPS of t is used
 * as the initialisation dump for the model.
 */
int dump_choose_by_time(parameters_t *params, int fid, double t)
{
  unsigned int i;
  size_t n;
  size_t start;
  size_t count;
  double tvals;
  char timeunits[MAXSTRLEN];

/* This is a little naughty. The times should have been mapped
 * to a common base, and then compared.
 */
  if (strcmp(schedule->units, params->timeunit) != 0)
    hd_quit("dump_choose_by_time: Time units are incompatible.\n");

  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &n);

  if (DEBUG("dump"))
    dlog("dump", "%d dumps in input file.\n", n);

  if (n < 1)
    hd_quit("dump_choose_by_time: No dumps in input file '%s'!\n",
            params->idumpname);
  memset(timeunits, 0, MAXSTRLEN);
  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
  for (i = 0; i < n; ++i) {
    start = i;
    count = 1;
    nc_get_vara_double(fid, ncw_var_id(fid, "t"), &start, &count, &tvals);
    tm_change_time_units(timeunits, params->timeunit, &tvals, 1);
    if (fabs(tvals - t) < START_EPS)
      return i;
  }

  hd_quit
    ("dump_choose_by_time: The dump file '%s' does not contain the time '%.2f %s'.\n",
     params->idumpname, t, params->timeunit);

  return -1;
}

int dump_choose_by_time_m(master_t *master, int fid, double t)
{
  unsigned int i;
  size_t n = 0;
  size_t start;
  size_t count;
  double tvals;
  char timeunits[MAXSTRLEN];
  if (nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &n) < 0)
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, "time"), &n) < 0)
      hd_quit("Can't find 'record' dimension in netCDF file.\n");

  if (n < 1)
    hd_quit("dump_choose_by_time: No dumps in input file!\n");
  memset(timeunits, 0, MAXSTRLEN);
  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
  for (i = 0; i < n; ++i) {
    start = i;
    count = 1;
    nc_get_vara_double(fid, ncw_var_id(fid, "t"), &start, &count, &tvals);
    tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
    if (i == 0 && t < tvals) return i;
    if (i == n-1 && t > tvals) return n-1;
    if (fabs(tvals - t) < START_EPS)
      return i;
  }

  hd_quit
    ("dump_choose_by_time: Input file does not contain the time '%.2f %s'.\n",
     t, master->timeunit);

  return -1;
}

int dump_choose_by_time_p(parameters_t *params, int fid, double t)
{
  unsigned int i;
  size_t n;
  size_t start;
  size_t count;
  double tvals;
  char timeunits[MAXSTRLEN];

  nc_inq_dimlen(fid, ncw_dim_id(fid, "record"), &n);

  if (n < 1)
    hd_quit("dump_choose_by_time: No dumps in input file '%s'!\n",
            params->idumpname);
  memset(timeunits, 0, MAXSTRLEN);
  nc_get_att_text(fid, ncw_var_id(fid, "t"), "units", timeunits);
  for (i = 0; i < n; ++i) {
    start = i;
    count = 1;
    nc_get_vara_double(fid, ncw_var_id(fid, "t"), &start, &count, &tvals);
    tm_change_time_units(timeunits, params->timeunit, &tvals, 1);
    if (fabs(tvals - t) < START_EPS)
      return i;
  }

  hd_quit
    ("dump_choose_by_time: The dump file '%s' does not contain the time '%.2f %s'.\n",
     params->idumpname, t, params->timeunit);

  return -1;
}

int dump_choose_by_time_mom(master_t *master, int fid, double t)
{
  unsigned int i;
  size_t n;
  size_t start;
  size_t count;
  double tvals, pt;
  char timeunits[MAXSTRLEN];

  nc_inq_dimlen(fid, ncw_dim_id(fid, "time"), &n);

  if (n < 1)
    hd_quit("dump_choose_by_time: No dumps in input file!\n");
  memset(timeunits, 0, MAXSTRLEN);
  nc_get_att_text(fid, ncw_var_id(fid, "Time_bounds"), "units", timeunits);

  for (i = 0; i < n; ++i) {
    start = i;
    count = 1;
    nc_get_vara_double(fid, ncw_var_id(fid, "time"), &start, &count, &tvals);
    tm_change_time_units(timeunits, master->timeunit, &tvals, 1);
    if (i == 0 && t < tvals) return i;
    if (i == n-1 && t > tvals) return n-1;
    if (fabs(tvals - t) < START_EPS)
      return i;
    if (i > 0 && t > pt && t < tvals)
      return i;
    pt = tvals;
  }

  hd_quit
    ("dump_choose_by_time: Input file does not contain the time '%.2f %s'.\n",
     t, master->timeunit);

  return -1;
}


int dump_choose_by_time_ts(master_t *master, char *fname, double t)
{
  unsigned int i;
  double tvals;
  double r = t;
  timeseries_t *ts = (timeseries_t *)malloc(sizeof(timeseries_t));
  datafile_t *df;

  if (ts == NULL)
    hd_quit("interp_data_s: No memory available.\n");
  memset(ts, 0, sizeof(timeseries_t));
    
  /* Read the time series                                            */
  ts_read(fname, ts);
  df = ts->df;
  tm_change_time_units(master->timeunit, ts->t_units, &r, 1);

  r = get_file_time(ts, r);
  if (df->rec_modulus) {
    while (r < 0)
      r += df->rec_mod_scale;
    r = fmod(r, df->rec_mod_scale);
  }

  for (i = 0; i < df->nrecords; ++i) {
    tvals = df->records[i];
    if (i == 0 && r < tvals) return i;
    if (i == df->nrecords-1 && r > tvals) return df->nrecords-1;
    if (i < df->nrecords-1 && r >= tvals && r < df->records[i+1]) {
      ts_free(ts);
      return i;
    }
  }

  hd_quit
    ("dump_choose_by_time: Input file does not contain the time '%.2f %s'.\n",
     t, master->timeunit);

  ts_free(ts);
  return(-1);
}



