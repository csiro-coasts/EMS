/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/timeseries.c
 *
 *  \brief Interface to timeseries data files
 *
 *  Routines which deal with time series data. There are several ways
 *  of defining a time series.
 *
 *       1. Ascii file of times and values
 *          If the file is an ascii file, it must have the
 *          following format:
 *
 *          # Comments
 *          ## COLUMNS n
 *          ##
 *          ## COLUMN1.name XXXX
 *          ## COLUMN1.long_name XXXX
 *          ## COLUMN1.units XXXX
 *          ## COLUMN1.missing_value XXXX
 *          ##
 *          ## COLUMN2.name XXXX
 *          ## COLUMN2.long_name XXXX
 *          ## COLUMN2.units XXXX
 *          ## COLUMN2.missing_value XXXX
 *          ##
 *               .
 *               .
 *           .
 *          ##
 *          v   v   v   v   ...
 *          v   v   v   v   ...
 *          v   v   v   v   ...
 *           .
 *           .
 *           .
 *
 *       2. netCDF file with time and 0, 1 or 2 spatial
 *          dimensions.
 *
 *       In either case, the timeseries will contain all
 *       the variables in the specified file.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: timeseries.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "ems.h"

void ascii_write(FILE * fp, datafile_t *df);
static void check_modulus(timeseries_t *ts, int chk);

#define TEPS (0.001)

/** Read and set up a time series.
  * 
  * @param name File containing the time series data.
  * @param ts pointer to time series structure (assumed
  * 		 not previously initialised)
  * 
  * This routine calls quit() if anything goes wrong
  */
void ts_read(char *name, timeseries_t *ts)
{
  int i;
  int index;
  datafile_t *df;

  if ((df = (datafile_t *)malloc(sizeof(datafile_t))) == NULL)
    quit("ts_read: Failed to allocate memory for Datafile.\n");
  memset(df, 0, sizeof(datafile_t));
  emslog(LTRACE,"Reading timeseries from file: %s \n",name);
  
  df_read(name, df);
  if (df->records == NULL) {
    index = df_get_index(df, "t");
    if (index >= 0)
      df_set_record(df, index);
    else if ((index = df_get_index(df, "time")) >= 0)
      df_set_record(df, index);
    else if ((index = df_get_index(df, "Time")) >= 0)
      df_set_record(df, index);
    else if ((index = df_get_index(df, "TIME")) >= 0)
      df_set_record(df, index);
  }

  /* Copy the convenience variables over */
  if (df->type != DFT_ASCII) {
    switch (df->nd) {
    case 0:
      ts->type = TS_NC0;
      break;

    case 1:
      ts->type = TS_NC1;
      break;

    case 2:
      ts->type = TS_NC2;
      break;

    case 3:
      ts->type = TS_NC3;
      break;

    default:
      break;
    }
  } else
    ts->type = TS_ASCII;

  ts->name = df->name;
  ts->nt = df->nrecords;
  ts->ti = df->ri;
  ts->t_units = df->rec_units;
  ts->t = df->records;
  ts->nv = df->nv;

  /* Copy variables */
  if ((ts->varname = (char **)malloc(sizeof(char *) * df->nv)) == NULL)
    quit("ts_read: Failed to allocate memory for variable names.\n");
  memset(ts->varname, 0, sizeof(char *) * df->nv);
  for (i = 0; i < df->nv; ++i) {
    ts->varname[i] = (char *)malloc(sizeof(char)
                                    * (strlen(df->variables[i].name) + 1));
    if (ts->varname[i] == NULL)
      quit("ts_read: Failed to allocate memory.\n");
    strcpy(ts->varname[i], df->variables[i].name);
  }

  if ((ts->varunit = (char **)malloc(sizeof(char *) * df->nv)) == NULL)
    quit("ts_read: Failed to allocate memory for variable units.\n");
  memset(ts->varunit, 0, sizeof(char *) * df->nv);
  for (i = 0; i < df->nv; ++i) {
    if (df->variables[i].units != NULL) {
      ts->varunit[i] = (char *)malloc(sizeof(char)
                                      * (strlen(df->variables[i].units) +
                                         1));
      if (ts->varunit[i] == NULL)
        quit("ts_read: Failed to allocate memory.\n");
      strcpy(ts->varunit[i], df->variables[i].units);
    }
  }

  ts->df = df;
  check_modulus(ts,1);
}

/** Set the form in which all geographic data should be
  * returned.
  * @param ts pointer to time series structure
  * @param type character string containing the geotype.
  * NULL - no prefered output form, no transformation occurs,
  * data same as file format.
  * "geographic" - All 'geographic' related data should be
  * returned as latitudes and longitudes.
  * "proj=xxx .." - Standard map_proj_t format.
  */
void ts_set_proj_type(timeseries_t *ts, char *type)
{
  df_set_proj_type(ts->df, type);
}

/** Set the default geographic data type.
  * @param type character string containing the geotype.
  * NULL - no prefered output form, no transformation occurs,
  * data same as file format.
  * "geographic" - All 'geographic' related data should be
  * returned as latitudes and longitudes.
  * "proj=xxx .." - Standard map_proj_t format.
  */
void ts_set_default_proj_type(char *type)
{
  df_set_default_proj_type(type);
}

/** Set the default hash_table_t size for evaluating numeric grids.
  * @param hsize Hashtable size.
  */
void ts_set_default_hashtable_size(int hsize)
{
  df_set_default_hashtable_size(hsize);
}


/** See if a specific time point exists
  * 
  * @param ts pointer to time series structure
  * @param t specified sample time
  * @return 1 for yes, 0 otherwise
  */
int ts_has_time(timeseries_t *ts, double t)
{
  int b,a;
  double f;
  return(df_find_record(ts->df, t, &b, &a, &f));
}

/** Function that returns the runcode for a timeseries file
 *
 */
int ts_eval_runcode(timeseries_t *ts)
{
  datafile_t   *df = ts->df;
  if (df->runcode == DF_FAIL)
    return 1;
  return 0;
}


/** Evaluate a time series at a specified time.
  * 
  * @param ts pointer to time series structure
  * @param id index of variable to evaluate
  * @param t time value
  * @return value at time t.
  *
  * Calls quit() if something goes wrong.
  */
double ts_eval(timeseries_t *ts, int id, double t)
{
  df_variable_t *v;

  if (ts == NULL)
    quit("ts_eval: NULL TimeSeries pointer\n");

  v = df_get_variable(ts->df, id);
  if (v == NULL)
    quit("ts_eval: No variable corresponding to id = %d\n", id);

  return df_eval(ts->df, v, get_file_time(ts, t));
}


/** Wrapper function that returns a flag for missing values
 *
 */
int ts_eval_xy_flag(timeseries_t *ts, int id, double t, double x, double y,
		    double *out)
{
  datafile_t   *df = ts->df;
  df_variable_t *v = df_get_variable(df, id);

  double val = ts_eval_xy(ts, id, t, x, y);
  /* Ignore all invalid values */
  if (
      fabs(val - v->lflag)     < 1e-5 || 
      fabs(val - v->cflag)     < 1e-5 ||
      fabs(val - v->missing)   < 1e-5 || 
      fabs(val - v->fillvalue) < 1e-5
      ) {
    return(0);
  } else {
    *out = val;
    return(1);
  }
}

/** Evaluate a time series at a specified time and place.
  * 
  * @param ts pointer to time series structure
  * @param id index of variable to evaluate
  * @param t time value
  * @param x x value
  * @param y y value
  * @return value at (t,x,y).
  *
  * Calls quit() if something goes wrong.
  */
double ts_eval_xy(timeseries_t *ts, int id, double t, double x, double y)
{
  datafile_t *df = ts->df;
  df_variable_t *v = df_get_variable(df, id);
  double coords[MAXNUMCOORDS];
  int i;
  int nc = 0;
  int *coordtypes = NULL;

  if (v == NULL)
    quit("ts_eval_xy: Invalid variable id specified (%d).\n", id);

  /* If 0 dimension, then the best we can do is give the interpolated time 
     value */
  if (v->nd == 0)
    return df_eval(df, v, get_file_time(ts, t));

  /* Define the coordinate system if one is not specified. */
  if (v->csystem == NULL) {
    /* Guess by looking at the attributes */
    if (!df_infer_coord_system(df, v)) {
      /* Failed to inter a coordinate system. To be backwardly compatible, 
         we should look for the variables 'x' and 'y'. */
      int x_id = df_get_index(df, "x");
      int y_id = df_get_index(df, "y");

      if ((x_id != -1) && (y_id != -1)) {
        df_coord_mapping_t map;
        int coordids[MAXNUMCOORDS];
        int coordtypes[MAXNUMCOORDS];
        df_variable_t *vx = df_get_variable(df, x_id);
        df_variable_t *vy = df_get_variable(df, y_id);

        df->variables[x_id].type = VT_X;
        df->variables[y_id].type = VT_Y;

        /* If the x and y coordinate have a single dimension only, then we 
           assume the data represents point data. */
        if ((vx->nd < 2) && (vy->nd < 2)) {
        }


        /* If both variables have 2 dimensions we can shall assume a grid
           coordinate system. We shall also try to provide analytic
           attributes if these already exist in the file as global
           attributes. */
        else if ((vx->nd == 2) && (vy->nd == 2)) {
          int dim1size = df->dimensions[vx->dimids[0]].size;
          int dim2size = df->dimensions[vx->dimids[1]].size;

          df_attribute_t *gt = df_get_global_attribute(df, "gridtype");
          if (gt != NULL) {
            char analytic[MAXLINELEN];

            analytic[0] = 0;
            /* We have a analytic expression. Determine the type and add
               this as as an attribute to the x and y variables in the
               datafile structure. Always assume cell edge. */
            if (strcasecmp(ATT_TEXT(gt), "rectangular") == 0) {
              df_attribute_t *x0 = df_get_global_attribute(df, "xorigin");
              df_attribute_t *y0 = df_get_global_attribute(df, "yorigin");
              df_attribute_t *dx = df_get_global_attribute(df, "dx");
              df_attribute_t *dy = df_get_global_attribute(df, "dy");
              df_attribute_t *rot =
                df_get_global_attribute(df, "rotation");
              sprintf(analytic, "rectangular 0.0 0.0 %d %d %g %g %g %g %g",
                      dim1size, dim2size, ATT_DOUBLE(x0, 0), ATT_DOUBLE(y0,
                                                                        0),
                      ATT_DOUBLE(dx, 0), ATT_DOUBLE(dy, 0), ATT_DOUBLE(rot,
                                                                       0));
            }

            else if (strcasecmp(ATT_TEXT(gt), "polar") == 0) {
              df_attribute_t *x0 = df_get_global_attribute(df, "xorigin");
              df_attribute_t *y0 = df_get_global_attribute(df, "yorigin");
              df_attribute_t *arc = df_get_global_attribute(df, "arc");
              df_attribute_t *rmin = df_get_global_attribute(df, "rmin");
              df_attribute_t *rot =
                df_get_global_attribute(df, "rotation");
              sprintf(analytic, "polar 0.0 0.0 %d %d %g %g %g %g %g",
                      dim1size, dim2size, ATT_DOUBLE(x0, 0), ATT_DOUBLE(y0,
                                                                        0),
                      ATT_DOUBLE(arc, 0), ATT_DOUBLE(rmin, 0),
                      ATT_DOUBLE(rot, 0));
            }

            if (analytic[0] != 0) {
              df_add_text_attribute(df, vx, "analytic", analytic);
              df_add_text_attribute(df, vy, "analytic", analytic);
            }

          }
        }

        else
          quit
            ("ts_eval_xy: Unable to locate any appropriate coordinate variables for variable '%s'.\n",
             v->name);



        /* Build the coordinate map and try to set the coordinate. */
        map.nc = 2;
        map.coordids = coordids;
        map.coordtypes = coordtypes;
        map.coordids[0] = x_id;
        map.coordtypes[0] = VT_X;
        map.coordids[1] = y_id;
        map.coordtypes[1] = VT_Y;
        if (df_set_coord_system(df, v, 1, &map) == 0)
          quit
            ("ts_eval_xy: Unable to locate any appropriate coordinate variables for variable '%s'.\n",
             v->name);
      } else
        quit
          ("ts_eval_xy: Unable to locate any appropriate coordinate variables for variable '%s'.\n",
           v->name);
    }
  }

  nc = df_get_num_coords(df, v);
  coordtypes = df_get_coord_types(df, v);
  if (nc != 2)
    quit("ts_eval_xy: Expected only 2 coordinates (X,Y) but had %d.\n",
         nc);

  for (i = 0; i < nc; ++i) {
    if (coordtypes[i] & (VT_X | VT_LONGITUDE))
      coords[i] = x;
    else if (coordtypes[i] & (VT_Y | VT_LATITUDE))
      coords[i] = y;
    else
      quit
        ("ts_eval_xy: The coordinate maps do not contain XY coordinate types.\n");
  }

  return df_eval_coords(df, v, get_file_time(ts, t), coords);
}


/** Evaluate a time series at a specified time and 3d place.
  * 
  * @param ts pointer to time series structure
  * @param id index of variable to evaluate
  * @param t time value
  * @param x x value
  * @param y y value
  * @param z z value
  * @return value at (t,x,y,z).
  *
  * Calls quit() if something goes wrong.
  */
double ts_eval_xyz(timeseries_t *ts, int id, double t, double x, double y,
                   double z)
{
  datafile_t *df = ts->df;
  df_variable_t *v = df_get_variable(df, id);
  double coords[MAXNUMCOORDS];
  int i;
  int nc = 0;
  int *coordtypes = NULL;

  if (v == NULL)
    quit("ts_eval_xyz: Invalid variable id specified (%d).\n", id);

  /* If 0 dimension, then the best we can do is give the interpolated time 
     value */
  if (v->nd == 0)
    return df_eval(df, v, get_file_time(ts, t));

  else if (v->nc == 2)
    return ts_eval_xy(ts, id, t, x, y);


  /* Define the coordinate system if one is not specified. */
  if (v->csystem == NULL) {

    /* Guess by looking at the attributes */
    if (!df_infer_coord_system(df, v))
      quit("ts_eval_xyz: Unable to infer the coordinate system.");
  }

  nc = df_get_num_coords(df, v);
  coordtypes = df_get_coord_types(df, v);
  if (nc == 2)
    return ts_eval_xy(ts, id, t, x, y);
  else if (nc != 3)
    quit("ts_eval_xyz: Expected only 3 coordinates (X,Y,Z) but had %d.\n",
         nc);

  for (i = 0; i < nc; ++i) {
    if (coordtypes[i] & (VT_X | VT_LONGITUDE))
      coords[i] = x;
    else if (coordtypes[i] & (VT_Y | VT_LATITUDE))
      coords[i] = y;
    else if (coordtypes[i] & (VT_Z))
      coords[i] = z;
    else
      quit
        ("ts_eval_xyz: The coordinate maps do not contain XYZ coordinate types.\n");
  }

  return df_eval_coords(df, v, get_file_time(ts, t), coords);
}


/** Check whether the specified variable located in a series
  * of timeseries file. Also check whether the variable spans
  * the range.
  *
  * The time check is still fairly simple minded. It does not
  * look for holes in the data.
  *
  * @param ntsfiles Number of timeseries files in array.
  * @param tsfiles Array of timeseries files.
  * @param var Name of variable to check.
  * @param tstart Lower end of range to check.
  * @param tstop Upper end of range to check.
  * @return non-zero value if successful.
  */
int ts_multifile_check(int ntsfiles, timeseries_t *tsfiles, char *var,
                       double tstart, double tstop)
{
  double f_tstart = 1e300;
  double f_tstop = -1e300;
  int found_var = 0;

  if (ntsfiles > 0) {
    int i;
    for (i = 0; i < ntsfiles; ++i) {
      timeseries_t *ts = &tsfiles[i];
      int varid = ts_get_index(ts, var);

      if (varid >= 0) {
        found_var = 1;
        if (ts->t_mod_type == MOD_NONE) {
          f_tstart = (ts->t[0] < f_tstart) ? ts->t[0] : f_tstart;
          f_tstop =
            (ts->t[ts->nt - 1] > f_tstop) ? ts->t[ts->nt - 1] : f_tstop;
        } else {
          f_tstart = (tstart < f_tstart) ? tstart : f_tstart;
          f_tstop = (tstop > f_tstop) ? tstop : f_tstop;
        }
      }
    }

    if (!found_var) {
      quit
        ("ts_multifile_check: Unable to located the variable '%s' in the timeseries files.\n",
         var);
      return 0;
    } else {
      if (!((tstart >= f_tstart) && (tstop <= f_tstop)))
        warn
          ("ts_multifile_check: The variable '%s' does not span the\ntimeseries files for the time range %.10g to %.10g (%.10g to %.10g).\n",
           var, tstart, tstop, f_tstart, f_tstop);
    }
  } else {
    quit("ts_multifile_check: No timeseries files specified.\n");
    return 0;
  }

  return 1;

}


/** Evaluate in the first appropriate file the specified variable
  * at the given point and time.
  *
  * @param ntsfiles Number of timeseries files in array.
  * @param tsfiles Array of timeseries files.
  * @param var Name of variable to check.
  * @param t time value
  * @param x x value
  * @param y y value
  * @param z z value
  * @return non-zero value if successful.
  */
double ts_multifile_eval_xyz(int ntsfiles, timeseries_t *tsfiles,
                             char *var, double t, double x, double y,
                             double z)
{
  if (ntsfiles > 0) {
    int i;
    int lastvalidindex = -1;
    for (i = 0; i < ntsfiles; ++i) {
      timeseries_t *ts = &tsfiles[i];
      int varid = ts_get_index(ts, var);

      if (varid >= 0)
        lastvalidindex = i;

      /* If the variable exists in the file and the current time is within 
         the bounds of the timeseries files domain, then evaluate the
         point. */
      if ((varid >= 0)
          && ((ts->t_mod_type == MOD_NONE) ||
              ((t >= ts->t[0]) && (t <= ts->t[ts->nt - 1])))) {
        return ts_eval_xyz(ts, varid, t, x, y, z);
      }

      /* If this is the last timeseries file, then evaluate the point. */
      else if (i == (ntsfiles - 1)) {
        if (lastvalidindex >= 0) {
          ts = &tsfiles[lastvalidindex];
          varid = ts_get_index(ts, var);
          return ts_eval_xyz(ts, varid, t, x, y, z);
        } else
          quit
            ("ts_multifile_eval_xyz: Unable to evaluate the variable '%s'.\n",
             var);
      }
    }
  } else
    quit("ts_multifile_eval_xyz: No timeseries files specified.\n");

  return 0.0;
}


/** Find the index of a variable in a time series
  * 
  * @param ts pointer to time series structure
  * @param name name of variable to find
  * @return index value (>= 0) if sucessful, -1 otherwise.
  */
int ts_get_index(timeseries_t *ts, char *name)
{
  return df_get_index(ts->df, name);
}

/** Gives the name of the variable with the given index
  * 
  * @param ts pointer to time series structure
  * @param vid id of variable
  */
const char* ts_get_varname(timeseries_t *ts, int vid)
{
  return ts->varname[vid];
}

/** Find the value of an attribute associated with a variable in a time series
  * 
  * @param ts pointer to time series structure
  * @param name name of variable to find
  * @param attname name of attribute to find
  * @return index value (>= 0) if sucessful, -1 otherwise.
  */
double ts_get_att_value(timeseries_t *ts, char *name, char *attname)
{
  int j;
  double val;
  datafile_t *df = ts->df;
  df_variable_t *v = df_get_variable_by_name(df, name);
  for (j = 0; j < v->na; ++j) {
    df_attribute_t *a = &v->attributes[j];
    if (strcasecmp(a->name, attname) == 0) {
      if (CHK_TYPE(a, AT_BYTE))
	val = ATT_BYTE(a, 0);
      else
	val = ATT_DOUBLE(a, 0);
      return(val * v->scale_factor + v->add_offset);
    }
  }
  return(-1.0);
}


/** Change time units in a time series
  *
  * @param ts pointer to time series structure
  * @param newunits new units string
  *
  * This routine calls quit() if anything goes wrong
  */
void ts_convert_time_units(timeseries_t *ts, char *newunits)
{
  datafile_t *df = ts->df;

  tm_change_time_units(ts->t_units, newunits, df->records, df->nrecords);

  /* Store new units in time series */
  if ((ts->t_units =
       (char *)realloc((void *)ts->t_units, strlen(newunits) + 1)) == NULL)
    quit("ts_convert_time_units: Can't allocate memory for new units\n");
  strcpy(ts->t_units, newunits);

  /* Propogate the change into the DataFile */
  df->rec_units = ts->t_units;
  df->variables[df->ri].units = ts->t_units;
  check_modulus(ts,0);
}

/**
  * Free the memory associated with the specified variable.
  *
  * @param ts pointer to timeseries structure.
  * @param varid The variable id.
  */
void ts_flush_variable(timeseries_t *ts, int varid) {
  df_flush_variable(ts->df, varid);
}

/** Free memory associated with a time series.
  *
  * @param ts pointer to time series structure
  *
  * This routine calls quit() if anything goes wrong.
  */
void ts_free(timeseries_t *ts)
{

  if (ts->varname != NULL) {
    int i;
    for (i = 0; i < ts->nv; ++i) {

      if (ts->varname[i] != NULL)
        free(ts->varname[i]);
      if (ts->varunit[i] != NULL)
        free(ts->varunit[i]);

    }
    free(ts->varname);
    free(ts->varunit);
  }

  df_free(ts->df);
}

/** Check that time values in a time series are
  * monotonic increasing.
  * 
  * @param ts pointer to time series structure
  * 
  * This routine calls quit() if anything goes wrong.
  */
void ts_is_time_monotonic(timeseries_t *ts)
{
  df_check_records(ts->df);
}

/*
 * Check RECOM attribute
 */
int ts_is_recom(timeseries_t *ts)
{
  return(df_is_recom(ts->df));
}

/** Check sign of variable
 */
int ts_var_z_is_depth(timeseries_t *ts, int id)
{
  return(ts->df->variables[id].z_is_depth);
}

/** Print summary information about a time series.
  *
  * @param ts pointer to time series structure
  * @param fp output FILE pointer
  *
  * This routine calls quit() if anything goes wrong.
  */
void ts_print_info(timeseries_t *ts, FILE * fp)
{
  df_print_info(ts->df, fp, 0);
}

/** Wrapper around ts_read but handles multiple files
 *
 * @param fnames whitespace delimited string containing all the filenames
 * @param ts  output pointer to multi-timeseries struct passed in
 * @return number of files found
 */
int ts_multifile_read(char *fnames, timeseries_t *ts)
{
  char *files[MAXSTRLEN * MAX_TS_FILES];
  int ntsfiles, i;
  
  /* Decode the string */
  ntsfiles = parseline(fnames, files, MAX_TS_FILES);

  /* Loop over an initialize eas struct */
  for (i=0; i<ntsfiles; i++) {
    ts_read(files[i], &ts[i]);
  }
  
  return(ntsfiles);
}


/** Wrapper around ts_get_index but handles multiple files
 *
 * @param  nfiles number of ts files
 * @param  ts pointer to the multi-timeseries struct
 * @param  varname name of variable to get the index of
 * @return the linear index of -1 if not found
 */
int ts_multifile_get_index(int nfiles, timeseries_t *ts, char *varname)
{
  int i, index = 0;
  
  for (i=0; i<nfiles; i++) {
    int n = ts_get_index(&ts[i], varname);
    
    /* Keep track of the index */
    if (n >= 0) {
      index += n;
      return(index);
    } else {
      /* Keep counting */
      index += ts[i].nv;
    }
  }

  /* If we're here, then we never found this variable */
  return(-1);
}

/** Wrapper around ts_eval but handles multiple files
 *
 * @param nfiles number of files in the ts array
 * @param ts  output pointer to multi-timeseries struct passed in
 * @param multi_varid variable id
 * @param t time point
 * @return number of files found
 */
double ts_multifile_eval(int nfiles, timeseries_t *ts, int multi_varid, double t)
{

  int v,i, index = 0;

  /* Loop over to recover the local index */
  for (i=0; i<nfiles; i++) {
    for (v=0; v<ts[i].nv; v++) {
      if (multi_varid == index)
	return(ts_eval(&ts[i], v, t));
      /* Increment global index */
      index++;
    }
  }
  
  /* This probably should be an error */
  return(NaN);
}


/** Gets the total number of variables across all of the ts files
 *
 * @param nfiles number of files in the ts array
 * @param ts  output pointer to multi-timeseries struct passed in
 */
int ts_multifile_get_nv(int nfiles, timeseries_t *ts)
{
  int i;
  int nv = 0;

  for (i=0; i<nfiles; i++) {
    nv += ts[i].nv;
  }

  return(nv);
}

/** If a modulus attribute was associated with the time, then
  * convert a 'user' time to a 'file' time.
  *
  * @param ts pointer to time series structure.
  * @param t 'user' time.
  * @return 'file' time.
  */
double get_file_time(timeseries_t *ts, double t)
{

  if (ts->t_mod_type == MOD_NONE)
    return t;

  else {
    double ms = ts->df->rec_mod_scale * ts->t_unit_scalef;
    int y = 0, mo = 0, d = 0, h = 0, mi = 0, sec = 0;
    double s = 0.0, e = 0.0, nt = 0.0;
    double mt = 0.0;

    /* Count the number of extra days for leap years */
    if(ts->df->rec_mod_scale == 365 * 86400) {
      int sy = 0, ey = 0;
      tm_to_julsecs(ts->t_base, &sy, &mo, &d, &h, &mi, &sec);
      tm_to_julsecs(t / 86400 + ts->t_base, &ey, &mo, &d, &h, &mi, &sec);
      mt = ts->t_base;
      for (y = sy; y < ey; y++)
	mt += ((y % 4 == 0 && (y % 100 != 0 || y % 4 == 0)) ? 366.0 : 365.0);
    } else
      mt = ((int)(t * ts->t_unit_scalef / (ms))) * ms / 86400.0 + ts->t_base;
    tm_to_julsecs(mt, &y, &mo, &d, &h, &mi, &sec);

    switch (ts->t_mod_type) {
    case MOD_YEAR:
      s = tm_to_juldays(y, 1, 1, 0, 0, 0);
      e = tm_to_juldays(y + 1, 1, 1, 0, 0, 0);
      break;

    case MOD_MONTH:
      s = tm_to_juldays(y, mo, 1, 0, 0, 0);
      if (mo < 12)
        e = tm_to_juldays(y, mo + 1, 1, 0, 0, 0);
      else
        e = tm_to_juldays(y + 1, 1, 1, 0, 0, 0);
      break;

    case MOD_WEEK:
      s = tm_to_juldays(y, mo, d, 0, 0, 0);
      e = s + 7;
      break;

    case MOD_DAY:
      s = tm_to_juldays(y, mo, d, 0, 0, 0);
      e = s + 1;
      break;

    case MOD_HOUR:
      s = tm_to_juldays(y, mo, d, h, 0, 0);
      e = s + 1 / 24.0;
      break;

    case MOD_MINUTE:
      s = tm_to_juldays(y, mo, d, h, mi, 0);
      e = s + 1 / 1440.0;
      break;

    case MOD_SECOND:
      tm_to_julsecs(mt, &y, &mo, &d, &h, &mi, &sec);
      s = tm_to_juldays(y, mo, d, h, mi, sec);
      e = s + 1 / 86400.0;
      break;
    }


    nt =
      fmod((t * ts->t_unit_scalef) / 86400.0 + ts->t_base - mt, ms) + mt;
    t = ms * ((nt - s) / (ts->t_mod_scale * (e - s))) / ts->t_unit_scalef;
  }

  return t;
}


/** Write the timeseries data as an ASCII file to the specified
  * file stream.
  *
  * @param fp output FILE pointer
  * @param ts pointer to file series structure
  */
void ts_write_ascii(FILE * fp, timeseries_t *ts)
{
  ascii_write(fp, ts->df);
}

/** Checks modulo flag for timeseries file
 *
 */
int ts_is_modulo(timeseries_t *ts)
{
  return(ts->t_mod_type > MOD_NONE);
}

static void check_modulus(timeseries_t *ts, int chk)
{
  int i, j;
  datafile_t *df = ts->df;

  /* Check for modulo */
  ts->t_mod_type = MOD_NONE;

  if (df->ri >= 0) {
    df_variable_t *tv = &df->variables[df->ri];
    for (i = 0; i < tv->na; ++i) {
      df_attribute_t *a = &tv->attributes[i];
      if ((strcasecmp(a->name, "modulo") == 0) && CHK_TYPE(a, AT_TEXT)) {
        char *mt = ATT_TEXT(a);
        char buf[128];

        if (isalpha((int)mt[0]))
          ts->t_mod_scale = 1.0;
        else {
          char buf[128];
          if (sscanf(mt, "%lf %s", &ts->t_mod_scale, buf) != 2)
            quit("ts_read: Invalid modulus attribute (%s).\n",
                 ATT_TEXT(a));
          strcpy(mt, buf);
        }

        df->rec_mod_scale = 86400.0 * ts->t_mod_scale;
        if (strcasecmp(mt, "year") == 0) {
          ts->t_mod_type = MOD_YEAR;
          df->rec_mod_scale *= 365.0;
        } else if (strcasecmp(mt, "month") == 0) {
          ts->t_mod_type = MOD_MONTH;
          df->rec_mod_scale *= 30.0;
        } else if (strcasecmp(mt, "week") == 0) {
          ts->t_mod_type = MOD_WEEK;
          df->rec_mod_scale *= 7.0;
        } else if (strcasecmp(mt, "day") == 0) {
          ts->t_mod_type = MOD_DAY;
          df->rec_mod_scale *= 1.0;
        } else if (strcasecmp(mt, "hour") == 0) {
          ts->t_mod_type = MOD_HOUR;
          df->rec_mod_scale *= 1.0 / 24.0;
        } else if (strcasecmp(mt, "minute") == 0) {
          ts->t_mod_type = MOD_MINUTE;
          df->rec_mod_scale *= 1.0 / 1440.0;
        } else if (strcasecmp(mt, "second") == 0) {
          ts->t_mod_type = MOD_SECOND;
          df->rec_mod_scale *= 1.0 / 86400.0;
        } else
          quit("ts_read: Unknown modulus attribute (%s).\n", ATT_TEXT(a));

        /* Determine the scale factor to convert to julian seconds. */
        sprintf(buf, "1 %s", ts->t_units);
        tm_scale_to_secs(buf, &ts->t_unit_scalef);
        df->rec_mod_scale /= ts->t_unit_scalef;
        ts->t_base = tm_time_to_julsecs(ts->t_units);

        /* 
	 * Check record values
	 *
	 * First time this function is called from ts_read and we do
	 * the time check
	 * When its called from ts_convert_time_units, we've already
	 * done the check and so now we need to time shift to start
	 * from zero
	 */
	if (chk)
	  for (j = 0; j < ts->nt; ++j) {
	    if ((ts->t[j] < 0) || (ts->t[j] > df->rec_mod_scale))
	      quit
		("df_check_modulus: Invalid record(%d) time(%lf) for modulo\ntime-series file. Should be between 0 and %lf.\n",
		 j, ts->t[j], df->rec_mod_scale);
	  }
	else {
	  // Rebase
	  double t0 = ts->t[0];
	  for (j = 0; j < ts->nt; ++j)
	    ts->t[j] -= t0;
	}
        break;
      }
    }
  }
}
