/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/ts_serial.cpp
 *
 *  \brief A serial buffer for 2D timeseries data in netCDF files
 *
 *  This can easily be generalised to nD
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id:$
 */

#include <stdlib.h>
#include <string.h>
#include <vector>
#include "ts_serial.hpp"

/* Default constructor */
tsSerial::tsSerial(char *fname)
{
  /*
   * Initialise timeseries object
   */ 
  ts = new timeseries_t;
  ts_read(fname, ts);

}

/* Default destructor */
tsSerial::~tsSerial(void) {
  ts_free(ts);
  delete ts;
}

/*
 * Runs through the data in the order of the dimensions and serialises
 * xxx assumes 2d lat/lons but it can be made more generic
 */
void tsSerial::serialise(char *vname)
{
  int i, nc, c0, c1, nc0, nc1;
  df_variable_t *var;
  df_variable_t *var_x = NULL, *var_y = NULL;
  datafile_t *df = ts->df;
  int *ctypes, *coordids;
  
  var = df_get_variable_by_name(ts->df, vname);
  if (var->nd != 2)
    quit("tsSerial: Expecting 2D timeseries variable '%s' in file '%s'\n",
	 vname);

  /* Assume all variables are the same coords, dimensions etc .. */
  if (!s2x.empty())
    return;

  /* Let EMS sort out the coordinate system */
  if (!df_infer_coord_system(df, var))
    quit("tsSerial: Unable to infer coordinate system for '%s'\n", var->name);

  nc = df_get_num_coords(df, var);
  if (nc != 2)
    quit("tsSerial: Expected only 2 coordinates (X,Y) but had %d.\n", nc);

  /* Grab x and y id's */
  ctypes   = df_get_coord_types(df, var);
  coordids = df_get_coord_ids(df, var);
  for (i = 0; i < nc; ++i) {
    if (ctypes[i] & (VT_X | VT_LONGITUDE))
      var_x = df_get_variable(df, coordids[i]);
    else if (ctypes[i] & (VT_Y | VT_LATITUDE))
      var_y = df_get_variable(df, coordids[i]);
  }
  
  if (var_x == NULL || var_y == NULL)
    quit("tsSerial: X and Y coordinates not found for '%s'\n", var->name);

  /* Get dimension lengths */
  nc0 = df->dimensions[var->dimids[0]].size;
  nc1 = df->dimensions[var->dimids[1]].size;

  /* Read the first and only records */
  df_read_record(df, var_x, 0);
  df_read_record(df, var_y, 0);

  /* Loop and build serial arrays */
  for (c0=0; c0<nc0; c0++) {
    for (c1=0; c1<nc1; c1++) {
      int coords[] = {c0, c1};
      /* Index space */
      s2c0.push_back(c0);
      s2c1.push_back(c1);
      /* Geographical */
      s2x.push_back(df_get_data_value(df, var_x, 0, coords));
      s2y.push_back(df_get_data_value(df, var_y, 0, coords));
    }
  }
}
