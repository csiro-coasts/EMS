/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/dfeval.c
 *
 *  \brief Evaluate datafile interpolation
 *
 *  Routines which deal interpolating datafile data
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: dfeval.c 7423 2023-10-05 02:15:55Z her127 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "ems.h"

/* local prototypes */
static double* get_dist_from_hash(hash_table_t** ht, int *coordids,
				  double coords[], int ni,
				  df_variable_t* vars);

static double find_close(datafile_t *df, df_variable_t *v, int r, int *is);
static double find_close_bathy(datafile_t *df, df_variable_t *v, int r, int *is);

/** Evaluate a zero dimension datafile variable for
  * a particular record value.
  *
  * @param df pointer to data file structure
  * @param v pointer to variable structure.
  * @param r record value
  * @return Interpolated value.
  *
  * @see quit() If anything goes wrong.
  */
double df_eval(datafile_t *df, df_variable_t *v, double r)
{
  double frac;
  int before;
  int after;
  double val;
  double v_before, v_after;
  
  /* Sanity checks */
  if (df == NULL)
    quit("df_eval: NULL Datafile pointer\n");
  if (v == NULL)
    quit("df_eval: NULL Variable pointer\n");
  if (v->nd > 0)
    quit("df_eval: Interpolation is ambiguous because data \
                 is multi-dimensional, use dfEvalCoords.\n");

  if ((df->records == NULL) || (!v->dim_as_record))
    return v->data[0];

  /* Find nearest records in table */
  if (df_find_record(df, r, &before, &after, &frac) == 0)
    /* Requested record not in range of those in file */
    /* NULL STATEMENT - Should take some action here */ ;

  if ((df->rec_modulus) && (before > after))
     after += df->nrecords;

/*  if (before > after)
    df_read_records(df, v, before, after + df->nrecords - before + 1);
  else */
    df_read_records(df, v, before, after - before + 1);

  val = v->missing;

  v_before = v->data[before - v->start_record];
  v_after  = v->data[after  - v->start_record];

  /* Value exactly at before record */
  if (frac == 0.0 && (fabs(v_before - v->missing) > 1e-10))
    val = v_before;
  /* Value exactly at after record */
  else if (frac == 1.0 && (fabs(v_after - v->missing) > 1e-10))
    val = v_after;
  /* Value inbetween - need to interpolate */
  else if ((fabs(v_after - v->missing) > 1e-10) &&
           (fabs(v_before - v->missing) > 1e-10))
    val = v_before * (1.0 - frac) + v_after * frac;

  return (val);
}

/**
 * Fixes up the sign of the depth variable
 *
 */
void df_fixup_variable_sign(datafile_t *df, int vid)
{
  df_variable_t *v = &df->variables[vid];
  if (v->z_is_depth) {
    int i;
    for(i=0; i<v->nrecords; i++)
      v->data[i] = (v->data[i] > 0) ? -v->data[i] : v->data[i];
  }
}


/** Evaluate a zero dimension datafile variable for
  * a period of record values. It DOES NOT interpolte and returns only
  * the bracketed values about the given time
  *
  * @param df pointer to data file structure
  * @param vid index of variable
  * @param r record value
  * @param dr half of this on either side
  * @param ptr pointer for the begin location
  * @return Number of records found
  *
  * @see quit() If anything goes wrong.
  */
int df_eval_period(datafile_t *df, int vid, double r, double dr, double **ptr)
{
  int before0, before1;
  int after0, after1;
  double frac;
  double *vals;
  double r0 = r - dr/2;
  double r1 = r + dr/2;
  int numrecs = 0;

  df_variable_t *v = &df->variables[vid];

  /* Sanity checks */
  if (df == NULL)
    quit("df_eval: NULL Datafile pointer\n");
  if (v == NULL)
    quit("df_eval: NULL Variable pointer\n");

  if (v->nd > 0)
    quit("df_eval_period: Data is multi-dimensional\n");

  if ((df->records == NULL) || (!v->dim_as_record))
    return 0;

  /* 
   * Find the beginning and end
   * If neither end is found then data is incomplete so bail out
   */
  df_find_record(df, r0, &before0, &after0, &frac);
  df_find_record(df, r1, &before1, &after1, &frac);

  /* How many */
  numrecs = after1 - before0 + 1;

  if (numrecs > 0) {
    /* Go ahead and read */
    df_read_records(df, v, before0, after1 - before0 + 1);
    
    /* 
     * Set output 
     *
     * ASCII files usually suck in the entire data while NETCDF will
     * only have the chunk we've read, hence the offset logic below
     */
    *ptr = &v->data[before0 - v->start_record];
  }

  return (numrecs);
}

/** Interpolate the variable at known coordinates and record.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to variable structure.
  * @param r record value.
  * @param coords coordinates values.
  * @return Interpolated value.
  *
  * @see dfGetNumCoordinates
  * @see dfGetCoordIds
  * @see quit() If anything goes wrong.
  */
double df_eval_coords(datafile_t *df, df_variable_t *v, double r,
                      double coords[])
{
  double val = 0.0;

  /* Sanity checks */
  if (df == NULL)
    quit("df_eval_coords: NULL Datafile pointer\n");
  if (v == NULL)
    quit("df_eval_coords: NULL Variable pointer\n");

  if (v->nd == 0) {
    val = df_eval(df, v, r);
  } else {
    int i;
    int r0, r1;
    double rfrac;
    double recvals[2];
    df_coord_system_t *cs = v->csystem;
    int nc = 0;
    int nd = 0;
    int in[2];

    recvals[0] = 0.0;
    recvals[1] = 0.0;
    in[0] = 0;
    in[1] = 0;

    if (cs == NULL)
      quit
        ("df_eval_coords: No coordinate system specified for variable '%s'\n",
         v->name);

    nc = df_get_num_coords(df, v);
    nd = df_get_num_dims(df, v);
    /* Find nearest records in table */
    if ((v->dim_as_record) && (df->records != NULL))
      df_find_record(df, r, &r0, &r1, &rfrac);
    else {
      r0 = r1 = 0;
      rfrac = 0.0;
    }

    if ((df->rec_modulus) && (r0 > r1))
      r1 += df->nrecords;

    /* Read in the data for each record and interpolate the coordinates
       for each record. Store the interpolated value for each record in
       the recvals array. */
    df_read_records(df, v, r0, r1 - r0 + 1);

    for (i = r0; i <= r1; ++i) {
      int j = i - r0;
      recvals[j] = 0.0;
      v->rid = j;

      /*
       * The interpolation function should've already been defined in
       * df_infer_coord_system using the same exact rules as those
       * that used to be here
       */
      if (nc >= nd) {
	if (v->interp != NULL)
	  recvals[j] = v->interp(df, v, i, coords);
	else
          quit("df_eval_coords: Unable to interpolate data with %d coordinate and %d dimensions. v->interp not defined.\n", nc, nd);
      } else
        quit("df_eval_coords: Less coordinates than dimensions.\n");
    }

    /* Interpolate the record in time */
    val = recvals[0] * (1.0 - rfrac) + recvals[1] * rfrac;
  }

  return val;
}

/*-------------------------------------------------------------------*/
/* PRIVATE and PROTECTED functions                                   */
/*-------------------------------------------------------------------*/
/* Interpolation routines assume that the order of variable          */
/* dimensions is fixed, with the vertical dimension being the        */
/* slowest varying (i.e. defined first), such that typically:        */
/* nz = df->dimensions[v->dimids[0]].size;                           */
/* nj = df->dimensions[v->dimids[1]].size;                           */
/* ni = df->dimensions[v->dimids[2]].size;                           */
/* The exception to this is ROMS station data interpolated with      */
/* interp3d_nearest(), where k varies the fastest. In this instance  */
/* data is accessed with d[n][k] rather than the ususal d[k][n].     */
/* Coordinate ids are retrieved using the df_get_coord_ids();        */
/* note that the order of these coordinate ids is different to that  */
/* accessed using v->coordids[], so do not use the latter.           */
/* Generally the coordids are arranged so that:                      */
/* x = longitude = coordids[0]                                       */
/* y = latitude = coordids[1]                                        */
/* z = depth = coordids[2] = coordids[nc-1] for nc = # coordinates   */
/* Note that this may be input file dependenent, and should be       */
/* checked if data from a new source is used.                        */
/*                                                                   */
/* Interpolation functions supported are:                            */
/* interp_1d_inv_weight : Inverse weighting for 2D point array data  */
/* interp_2d_inv_weight : Inverse weighting for 3D point array data  */
/* interp_1d_inv_weight_hashed : As for interp_1d_inv_weight but     */
/*     hashed weights.                                               */
/* interp_2d_inv_weight_hashed : As for interp_2d_inv_weight but     */
/*     hashed weights.                                               */
/* interp_us_2d() : Unstructured interpolation for 2D UGRID data     */
/* interp_us_3d() : Unstructured interpolation for 3D UGRID data     */
/* interp_us_2d_c() : Unstructured interpolation for 2D gridded      */
/*     Cartesian data.                                               */
/* interp_us_3d_c() : Unstructured interpolation for 3D gridded      */
/*     Cartesian data.                                               */
/* interp_us_2d_i() : Unstructured interpolation for 2D gridded      */
/*     Cartesian data of type VT_INFERRED.                           */
/* interp_us_3d_i() : Unstructured interpolation for 3D gridded      */
/*     Cartesian data of type VT_INFERRED.                           */
/* interp_nearest_within_eps() : Nearest neighbour with a tolerance  */
/*     for 2D point array data.                                      */
/* interp2d_nearest_within_eps() : Nearest neighbour with tolerance  */
/*     for 2D gridded data.                                          */
/* interp3d_nearest_within_eps() : Nearest neighbour with tolerance  */
/*     for 3D gridded data.                                          */
/* interp2d_nearest() : Nearest neighbour for 2D point array data.   */
/* interp3d_nearest() : Nearest neighbour for 3D (ROMS) station data */
/* interp_linear() : Bi/tri-linear interpolation for 2D/3D gridded   */
/*     data.                                                         */
/* interp_linear_degrees() : Same as interp_linear for data within   */
/*     the range 0-360 (degrees), e.g. tidal phases.                 */
/* interp_linear_filled() : Same as interp_linear but fills missing  */
/*     values with the nearest valid value.                          */
/* interp_linear_flagged() : Same as interp_linear but sets a        */
/*     no-gradient across land and doesn't interpolate over cloud    */
/*     values.                                                       */
/* interp_linear_bathy() : Same as interp_linear for bathymetry,     */
/*     fills missing or land (99) values with the nearest valid      */
/* value.                                                            */
/* See dfcoords.c for the rules as to how these functions are        */
/* applied.                                                          */
/*-------------------------------------------------------------------*/

/*
Routine to interpolate 1d spatial data.

NOTE: This is quite simple minded, using an inverse squared distance
weighting scheme. Because this routine is likely to be called twice
for each time series evaluation, the weights do not really need
to be re-calculated the second time. As well, you could use a hash
table mechanism to remember the weights for a large number of points.
This would be a good idea for use with numerical models, which tend
to evaluate time series repeatedly at a fairly small set of points.
For now, I haven't time to do this.
*/
/*
 * From 12/09 this version is no longer used as they have been
 * superseded by the _hashed functions that give considerable
 * performance gains. However I've left them here in case there is a
 * bug we can quickly swap them back. -FR
 */
double interp_1d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[])
{
  int i, j;
  int dsize = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double div = 0.0;
  double sum = 0.0;
  double *d = VAR_1D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;

  /* Calculate weighted sum of all data points */
  for (i = 0; i < dsize; ++i) {
    double dist = 0.0;

    for (j = 0; j < nc; ++j) {
      double dx = VAR_1D(&vars[coordids[j]])[0][i] - coords[j];
      dist += dx * dx;
    }

    dist = sqrt(dist);
    if (dist < 1e-10)
      return d[i];
    div += 1 / dist;
    sum += d[i] / dist;
  }

  return (sum / div);
}

/*
 * See comments for interp_1d_inv_weight
 */
double interp_2d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[])
{
  int i, k, n;
  int nk = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double **d = VAR_2D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  double *zgrid = VAR_1D(&vars[coordids[nc - 1]])[0];
  double z = coords[nc - 1];
  double tdata[2];
  double kindex = 0.0;
  int ks = 0;
  int kt = 0;

  /*UR  */
  tdata[0] = NaN;
  tdata[1] = NaN;
  
  /* Locate the k layer */
  if (z < zgrid[0]) {
    kindex = 0.0;
  } else if (z >= zgrid[nk - 1]) {
    kindex = nk - 1;
  } else {
    for (k = 0; k < nk - 1; ++k) {
      if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
        kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
        break;
      }
    }
  }

  ks = (int)kindex;
  kt = min(ks + 1, nk - 1);
  kindex = kindex - ks;

  for (k = ks; k <= kt; ++k) {
    double sum = 0;
    double div = 0.0;

    /* Calculate weighted sum of all data points */
    for (i = 0; i < ni; ++i) {
      double dist = 0.0;

      if (isnan(d[k][i]))
        continue;

      for (n = 0; n < nc - 1; ++n) {
        double dx = VAR_1D(&vars[coordids[n]])[0][i] - coords[n];
        dist += dx * dx;
      }

      dist = sqrt(dist);
      if (dist < 1e-10) {
        sum = d[k][i];
        div = 1.0;
        break;
      } else {
        div += 1 / dist;
        sum += d[k][i] / dist;
      }
    }
    emstag(LTRACE,"lib:dfeval:interp_2d_inv_weight","tdata index %d",(k - ks));
    tdata[k - ks] = sum / div;
  }

  if (isnan(tdata[0]))
    return tdata[1];
  else if (isnan(tdata[1]))
    return tdata[0];

  return (1 - kindex) * tdata[0] + (kindex) * tdata[1];
}


/*
 * This is the hashed version of the above
 * The optimisation saves quite a bit of time as:
 *   a) The grid boundary doesn't change within a model run so we save
 *      the distance from each of the points to each in the point
 *      array
 *   b) As several timesteps could fall in between 2 in the boundary
 *      file, there is no need to calculate each time
 *
 * The setup is as follows:
 *   There are 2 independent hashes:
 *      ht_dist           : hash of distances for each (x,y)
 *      ht_k (2d only)    : hash of kindex as a function of z
 *      ht_master_[1d|2d] : These are (x,y,z) or (x,y) hashes of
 *                          hashes{varid} of hash-queues{record #}
 *                          which gives the final interpolated value
 *
 *  See the code below for the algorithm
 */
double interp_1d_inv_weight_hashed(datafile_t *df, df_variable_t *v, int record,
				   double coords[])
{
  int ni = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);

  hash_table_t *ht_master = NULL;
  hash_table_t *ht_varid  = NULL;
  hqueue       *hq_val    = NULL;

  double *dists   = NULL;
  double *val_ptr = NULL;
  double  val     = 0.0;


  /*
   * Create the master (x,y) hash, if needed
   */
  if (df->ht_master_1d == NULL) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 2D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 2)
      quit("df_eval_coords : 1d_interp : nc is not 2 but %d!\n", nc);
    // This really should be the number of coords but we don't have
    // that information so this is the best we can do
    df->ht_master_1d = ht_create_d2(ni);
  }

  ht_master = df->ht_master_1d;

  /*
   * See if the value hash for this coord (x,y) exists
   */
  if ( (ht_varid=(hash_table_t *)ht_find(ht_master, &coords[0])) == NULL ) {
    /*
     * Create and add the hash table for this (x,y,z) coordinate
     * There are 2 entries for each variables
     */
    ht_varid = ht_create_i1(df->nv);
    ht_add(ht_master, &coords[0], ht_varid);
  }

  /*
   * Now check if the table for this varid exists
   */
  if ( (hq_val = (hqueue *)ht_find(ht_varid, &v->varid)) == NULL ) {
    /*
     * Only need two rows at the most
     */
    hq_val = hq_create(2);
    ht_add(ht_varid, &v->varid, hq_val);
  }
  
  /*
   * Finally, check to see if we have already calculated a value for
   * this time slot
   */
  if ( (val_ptr = (double *)hq_get_value_by_key(hq_val, record)) == NULL ) {
    int *coordids = df_get_coord_ids(df, v);
    double *data  = VAR_1D(v)[record - v->start_record];
    df_variable_t *vars = df->variables;

    void   *old_ptr  = NULL;
    int i;

    double sum = 0.0;
    double div = 0.0;

    /*
     * Get the distances vector
     */
    dists = get_dist_from_hash(&df->ht_dist, coordids, coords, ni, vars);
    
    /*
     * do the interpolation
     */
    for (i = 0; i < ni; ++i) {
      double dist = dists[i];
	
      if (isnan(data[i]))
	continue;
      
      if (dist < 1e-10) {
	sum = data[i];
	div = 1.0;
	break;
      } else {
	div += 1 / dist;
	sum += data[i] / dist;
      }
    }
    // Final value
    val = sum / div;
    
    // Assign and add to hash-queue
    val_ptr  = (double *)malloc(sizeof(double));
    *val_ptr = val;
    
    // Push returns the bumped value so free if not null
    old_ptr = hq_push(hq_val, record, val_ptr);
    if (old_ptr != NULL)
      free(old_ptr);
    
  } else {
    // We're done
    val = *val_ptr;
  }

  return val;
}


/*
 * This is the same as the 1D except we hash and interpolate in z
 */
double interp_2d_inv_weight_hashed(datafile_t *df, df_variable_t *v, int record,
				   double coords[])
{
  int nk = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);

  hash_table_t *ht_master = NULL;
  hash_table_t *ht_varid  = NULL;
  hqueue       *hq_val    = NULL;

  double  *val_ptr  = NULL;
  double  val       = 0.0;

  /*
   * Create the master (x,y,z) hash, if needed
   */
  if ( df->ht_master_2d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 3)
      quit("df_eval_coords : 2d_interp : nc is not 3 but %d!\n", nc);
    // This really should be the number of coords but we don't have
    // that information so this is the best we can do
    df->ht_master_2d = ht_create_d3(ni*nk);
  }

  ht_master = df->ht_master_2d;

  /*
   * See if the value hash for this coord (x,y,z) exists
   */
  if ( (ht_varid=(hash_table_t *)ht_find(ht_master, &coords[0])) == NULL ) {
    /*
     * Create and add the hash table for this (x,y,z) coordinate
     * There are 2 entries for each variables
     */
    ht_varid = ht_create_i1(df->nv);
    ht_add(ht_master, &coords[0], ht_varid);
  }

  /*
   * Now check if the table for this varid exists
   */
  if ( (hq_val = (hqueue *)ht_find(ht_varid, &v->varid)) == NULL ) {
    /*
     * Only need two rows at the most
     */
    hq_val = hq_create(2);
    ht_add(ht_varid, &v->varid, hq_val);
  }

  /*
   * Finally, check to see if we have already calculated a value for
   * this time slot
   */
  if ( (val_ptr = (double *)hq_get_value_by_key(hq_val, record)) == NULL ) {
    int *coordids = df_get_coord_ids(df, v);
    double **data = VAR_2D(v)[record - v->start_record];
    df_variable_t *vars = df->variables;

    double *kindex_ptr = NULL;
    double *dists      = NULL;
    double  kindex   = 0.0;
    double  tdata[2] = {NaN, NaN};
    void   *old_ptr  = NULL;
    int i,k,ks,kt;

    /*
     * Make sure the depth hash exists
     */
    if (df->ht_k == NULL) {
      df->ht_k = ht_create_d1(nk);
    }

    /*
     * Find the k layer
     */
    if ( (kindex_ptr = (double *)ht_find(df->ht_k,&coords[nc-1])) == NULL) {
      double *zgrid = VAR_1D(&vars[coordids[nc - 1]])[0];
      double  z = coords[nc - 1];
      
      if (z < zgrid[0]) {
	kindex = 0.0;
      } else if (z >= zgrid[nk - 1]) {
	kindex = nk - 1;
      } else {
	for (k = 0; k < nk - 1; ++k) {
	  if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
	    kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
	    break;
	  }
	}
      }
      // Add to hash
      kindex_ptr  = (double *)malloc(sizeof(double));
      *kindex_ptr = kindex;
      ht_add(df->ht_k, &coords[nc-1], kindex_ptr);
      
    } else {
      // This z coordinate is already in the hash
      kindex = *kindex_ptr;
    }
    
    ks = (int)kindex;
    kt = min(ks + 1, nk - 1);
    kindex = kindex - ks;

    /*
     * Get the distances vector
     */
    dists = get_dist_from_hash(&df->ht_dist, coordids, coords, ni, vars);
    
    /*
     * do the interpolation
     */
    for (k = ks; k <= kt; ++k) {
      double sum = 0.0;
      double div = 0.0;
      
      /* Calculate weighted sum of all data points */
      for (i = 0; i < ni; ++i) {
	double dist = dists[i];
	
	if (isnan(data[k][i]))
	  continue;
	
	if (dist < 1e-10) {
	  sum = data[k][i];
	  div = 1.0;
	  break;
	} else {
	  div += 1 / dist;
	  sum += data[k][i] / dist;
	}
      }
      emstag(LTRACE,"lib:dfeval:interp_2d_inv_weight_hashed","tdata index %d",
	     (k - ks));
      tdata[k - ks] = sum / div;
    }
    
    // interpolate between the 2 depth layers
    if (isnan(tdata[0]))
      val = tdata[1];
    else if (isnan(tdata[1]))
      val = tdata[0];
    else 
      val = (1 - kindex) * tdata[0] + (kindex) * tdata[1];
    
    // Assign and add to hash-queue
    val_ptr  = (double *)malloc(sizeof(double));
    *val_ptr = val;

    // Push returns the bumped value so free if not null
    old_ptr = hq_push(hq_val, record, val_ptr);
    if (old_ptr != NULL)
      free(old_ptr);
	   
  } else {

    // We're done
    val = *val_ptr;
  }

  return val;
}



/* Unstructured interpolation using Grid Spec libraries

 */
double interp_us_2d(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int dsize = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double *vd = VAR_1D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS *gs0 = NULL;
  GRID_SPECS *gs1 = NULL;
  int i, j;
  point *p;
  delaunay *d;
  double val;
  int rid = v->rid;
  GRID_SPECS *gs;
  /*char *i_rule = "linear";*/
  char i_rule[MAXSTRLEN];

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Create the grid spec function if required */
  if ( v->gs0_master_2d == NULL && v->gs1_master_2d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 2)
      quit("df_eval_coords : interp_us_2d : nc is not 2 but %d!\n", nc);

    /* Allocate and set the points structure */
    p = (point *)alloc_1d(dsize, sizeof(point));
    for (i = 0; i < dsize; ++i) {
      p[i].x = VAR_1D(&vars[coordids[0]])[0][i];
      p[i].y = VAR_1D(&vars[coordids[1]])[0][i];
    }

    /* Allocate and initialize the Delaunay structure */
    d = delaunay_build(dsize, p, 0, NULL, 0, NULL);
    for (i = 0; i < dsize; ++i) {
      d->points[i].v = d_alloc_1d(1);
      d->points[i].z = d->points[i].v[0] = vd[i];
    }

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = grid_spec_create();
    gs1 = grid_spec_create();
    grid_interp_init_t(gs0, d, i_rule, 0);
    gs0->d = d;
    v->gs0_master_2d = gs0;
    grid_interp_init_t(gs1, d, i_rule, 0);
    gs1->d = d;
    v->gs1_master_2d = gs1;
    v->rv[0] = v->rv[1] = -9999;
  }

 /* Reinitialize the weights if the data has changed */
  gs = (rid == 0) ? v->gs0_master_2d : v->gs1_master_2d;
  if (record != v->rv[rid]) {
    delaunay *d = gs->d;
    v->rv[rid] = record;
    for (i = 0; i < dsize; ++i)
      d->points[i].z = vd[i];
    if (strcmp(i_rule, "nn_sibson") == 0)
      for (i = 0; i < dsize; i++)
	gs->rebuild(gs->interpolator, &d->points[i]);
    else
      gs->rebuild(gs->interpolator, d->points);
  }
  val = grid_interp_on_point(gs, coords[0], coords[1]);

  return(val);
}

/* 2D Cartesian data. Missing values are stripped before */
/* computing the grid_spec.                              */
/* Note: this makes assumptions about the order of       */
/* dimensions and coordinates.                           */
double interp_us_2d_c(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int dsize;
  int nj = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double **vd = VAR_2D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS *gs0 = NULL;
  GRID_SPECS *gs1 = NULL;
  int n, i, j;
  point *p;
  delaunay *d;
  double val;
  int rid = v->rid;
  GRID_SPECS *gs;
  int hasmis = (v->missing) ? 1 : 0;
  int hasfil = (v->fillvalue) ? 1 : 0;
  char i_rule[MAXSTRLEN];
  double **gx = VAR_2D(&vars[coordids[0]])[0];
  double **gy = VAR_2D(&vars[coordids[1]])[0];

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Create the grid spec function if required */
  if ( v->gs0_master_2d == NULL && v->gs1_master_2d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 2)
      quit("df_eval_coords : interp_us_2d : nc is not 2 but %d!\n", nc);

    /* Allocate and set the points structure */
    dsize = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	dsize++;
      }
    }
    p = (point *)alloc_1d(dsize, sizeof(point));
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	p[n].x = gx[j][i];
	p[n++].y = gy[j][i];
      }
    }

    /* Allocate and initialize the Delaunay structure */
    d = delaunay_build(dsize, p, 0, NULL, 0, NULL);
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	d->points[n].v = d_alloc_1d(1);
	d->points[n].z = d->points[n].v[0] = vd[j][i];
	n++;
      }
    }

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = grid_spec_create();
    gs1 = grid_spec_create();
    grid_interp_init_t(gs0, d, i_rule, 0);
    gs0->d = d;
    v->gs0_master_2d = gs0;
    grid_interp_init_t(gs1, d, i_rule, 0);
    gs1->d = d;
    v->gs1_master_2d = gs1;
    v->rv[0] = v->rv[1] = -9999;
  }

 /* Reinitialize the weights if the data has changed */
  gs = (rid == 0) ? v->gs0_master_2d : v->gs1_master_2d;
  if (record != v->rv[rid]) {
    delaunay *d = gs->d;
    v->rv[rid] = record;
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	d->points[n++].z = vd[j][i];
      }
    }
    if (strcmp(i_rule, "nn_sibson") == 0)
      for (n = 0; n < dsize; n++)
	gs->rebuild(gs->interpolator, &d->points[n]);
    else
      gs->rebuild(gs->interpolator, d->points);
  }
  val = grid_interp_on_point(gs, coords[0], coords[1]);
  return(val);
}


/* Unstructured interpolation for VT_INFERRED files */
double interp_us_2d_i(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int dsize;
  int nj = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  int *coordtype = df_get_coord_types(df, v);
  double **vd = VAR_2D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS *gs0 = NULL;
  GRID_SPECS *gs1 = NULL;
  int n, i, j;
  point *p;
  delaunay *d;
  double val;
  double x, y;
  int rid = v->rid;
  GRID_SPECS *gs;
  int hasmis = (v->missing) ? 1 : 0;
  int hasfil = (v->fillvalue) ? 1 : 0;
  char i_rule[MAXSTRLEN];
  double *gx, *gy;

  /* Find which coordinates are lat, long */
  for (n = 0; n < nc; n++) {
    if (coordtype[n] & VT_LONGITUDE) {
      gx = VAR_1D(&vars[coordids[n]])[0];
      x = coords[n];
    }
    if (coordtype[n] & VT_LATITUDE) {
      gy = VAR_1D(&vars[coordids[n]])[0];
      y = coords[n];
    }
  }

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Create the grid spec function if required */
  if ( v->gs0_master_2d == NULL && v->gs1_master_2d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 2)
      quit("df_eval_coords : interp_us_2d : nc is not 2 but %d!\n", nc);

    /* Allocate and set the points structure */
    dsize = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	dsize++;
      }
    }
    p = (point *)alloc_1d(dsize, sizeof(point));
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	p[n].x = gx[i];
	p[n++].y = gy[j];
      }
    }

    /* Allocate and initialize the Delaunay structure */
    d = delaunay_build(dsize, p, 0, NULL, 0, NULL);
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	d->points[n].v = d_alloc_1d(1);
	d->points[n].z = d->points[n].v[0] = vd[j][i];
	n++;
      }
    }

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = grid_spec_create();
    gs1 = grid_spec_create();
    grid_interp_init_t(gs0, d, i_rule, 0);
    gs0->d = d;
    v->gs0_master_2d = gs0;
    grid_interp_init_t(gs1, d, i_rule, 0);
    gs1->d = d;
    v->gs1_master_2d = gs1;
    v->rv[0] = v->rv[1] = -9999;
  }

 /* Reinitialize the weights if the data has changed */
  gs = (rid == 0) ? v->gs0_master_2d : v->gs1_master_2d;
  if (record != v->rv[rid]) {
    delaunay *d = gs->d;
    v->rv[rid] = record;
    n = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[j][i] == v->missing)) continue;
	if ((hasfil && vd[j][i] == v->fillvalue)) continue;
	if (isnan(vd[j][i])) continue;
	d->points[n++].z = vd[j][i];
      }
    }
    if (strcmp(i_rule, "nn_sibson") == 0)
      for (n = 0; n < dsize; n++)
	gs->rebuild(gs->interpolator, &d->points[n]);
    else
      gs->rebuild(gs->interpolator, d->points);
  }
  val = grid_interp_on_point(gs, x, y);

  return(val);
}


/*------------------------------------------------------------------*/
/* Unstructured interpolation for UGRID files                       */
/*------------------------------------------------------------------*/
double interp_us_3d(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int nz = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double **vd = VAR_2D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS **gs0 = NULL;
  GRID_SPECS **gs1 = NULL;
  int i, j, k;
  point **p;
  delaunay **d;
  double val;
  double *zgrid = VAR_1D(&vars[coordids[nc - 1]])[0];
  double z = coords[nc - 1];
  double kindex   = 0.0;
  double  tdata[2] = {NaN, NaN};
  int ks, kt;
  /*char *i_rule = "linear";*/
  int rid = v->rid;
  GRID_SPECS **gs;
  int hasmis = (v->missing) ? 1 : 0;
  char i_rule[MAXSTRLEN];

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Create the grid spec function if required */
  if (v->gs0_master_3d == NULL && v->gs1_master_3d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 3)
      quit("df_eval_coords : interp_us_2d : nc is not 2 but %d!\n", nc);

    /* Allocate the private data */
    v->nz = nz;
    v->kn = i_alloc_1d(nz);
    v->kmap = i_alloc_2d(ni, nz);
    v->kflag = i_alloc_2d(2, nz);
    memset(v->kn, 0, nz * sizeof(int));

    /* Allocate and set the points */
    p = (point **)calloc(nz, sizeof(point *));
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ni; ++i) {
	if ((hasmis && vd[k][i] == v->missing)) continue;
	if (isnan(vd[k][i])) continue;
	  if (p[k] == NULL)
	    p[k] = (point *)alloc_1d(ni, sizeof(point));
	  v->kmap[k][v->kn[k]] = i;
	  p[k][v->kn[k]].x = VAR_1D(&vars[coordids[0]])[0][i];
	  p[k][v->kn[k]++].y = VAR_1D(&vars[coordids[1]])[0][i];

      }
    }

    /* Allocate and set the Delaunay triangulation */
    d = (delaunay **)calloc(nz+1, sizeof(delaunay *));
    for (k = 0; k < nz; k++) {
      d[k] = NULL;
      if (v->kn[k] == 0) continue;
      d[k] = delaunay_build(v->kn[k], p[k], 0, NULL, 0, NULL);
      d[k]->vid = 0;
      d[k]->ptf = 1;

      for (i = 0; i < v->kn[k]; i++) {
	j = v->kmap[k][i];
	d[k]->points[i].v = d_alloc_1d(1);
	d[k]->points[i].z = d[k]->points[i].v[0] = vd[k][j];
      }
    }

    /* Finished with temp array */
    free(p);

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs0_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    gs1 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs1_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    for (k = 0; k < nz; k++) {

      v->gs0_master_3d[k] = v->gs1_master_3d[k] = NULL;
      if (d[k] == NULL) continue;

      gs0[k] = grid_spec_create();
      gs1[k] = grid_spec_create();

      grid_interp_init_t(gs0[k], d[k], i_rule, 0);
      gs0[k]->d = d[k];
      v->gs0_master_3d[k] = gs0[k];

      grid_interp_init_t(gs1[k], d[k], i_rule, 0);
      gs1[k]->d = d[k];
      v->gs1_master_3d[k] = gs1[k];
    }
    v->rv[0] = v->rv[1] = -9999;
  }

  /* Find the k layer */      
  if (z < zgrid[0]) {
    kindex = 0.0;
  } else if (z >= zgrid[nz - 1]) {
    kindex = nz - 1;
  } else {
    for (k = 0; k < nz - 1; ++k) {
      if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
	kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
	break;
      }
    }
  }
  ks = (int)kindex;
  kt = min(ks + 1, nz - 1);
  /* No-gradient below bottom layer */
  if (!v->kn[ks]) ks = kt;
  kindex = kindex - ks;

  /* Reinitialize the weights if the data has changed */
  if (record != v->rv[rid]) {
    for (k = 0; k < nz; k++) v->kflag[rid][k] = 0;
    v->rv[rid] = record;
  }

  /* Rebuild the weights for the required layers */
  gs = (rid == 0) ? v->gs0_master_3d : v->gs1_master_3d;
  for (k = ks; k <= kt; ++k) {
    if (!v->kflag[rid][k]) {
      delaunay *d = gs[k]->d;
      for (i = 0; i < v->kn[k]; ++i) {
	j = v->kmap[k][i];
	d->points[i].z = vd[k][j];
      }
      v->kflag[rid][k] = 1;

      if (strcmp(i_rule, "nn_sibson") == 0)
	for (i = 0; i < ni; i++)
	  gs[k]->rebuild(gs[k]->interpolator, &d->points[i]);
      else
	gs[k]->rebuild(gs[k]->interpolator, d->points);
    }
  }

  /* Do the horizontal interpolation on layers bracketing the depth */
  for (k = ks; k <= kt; ++k) {
    tdata[k - ks] = grid_interp_on_point(gs[k], coords[0], coords[1]);
  }

  /* Do the vertical interpolation */
  if (isnan(tdata[0]))
    val = tdata[1];
  else if (isnan(tdata[1]))
    val = tdata[0];
  else 
    val = (1 - kindex) * tdata[0] + (kindex) * tdata[1];

  return(val);
}

/*------------------------------------------------------------------*/
/* Unstructured interpolation for gridded files (standard type)     */
/*------------------------------------------------------------------*/
double interp_us_3d_c(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int nz = df->dimensions[v->dimids[0]].size;
  int nj = df->dimensions[v->dimids[1]].size;
  int ni = df->dimensions[v->dimids[2]].size;
  int nc = df_get_num_coords(df, v);
  int dsize;
  int *coordids = df_get_coord_ids(df, v);
  int *coordtype = df_get_coord_types(df, v);
  double ***vd = VAR_3D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS **gs0 = NULL;
  GRID_SPECS **gs1 = NULL;
  int n, i, j, k;
  point **p;
  delaunay **d;
  double val;
  double x, y, z = coords[0];
  double kindex   = 0.0;
  double  tdata[2] = {NaN, NaN};
  int ks, kt;
  int rid = v->rid;
  GRID_SPECS **gs;
  int hasmis = (v->missing) ? 1 : 0;
  int hasfil = (v->fillvalue) ? 1 : 0;
  char i_rule[MAXSTRLEN];
  /*double *zgrid = VAR_1D(&vars[coordids[0]])[0];*/
  double zgrid[nz];
  double **gx, **gy;
  int zs, zb;
  int sigma = 0;

  /* Find which coordinates are lat, long, depth */
  for (n = 0; n < nc; n++) {
    if (coordtype[n] & VT_LONGITUDE) {
      gx = VAR_2D(&vars[coordids[n]])[0];
      /*if (nc == v->nd) ni = df->dimensions[coordids[n]].size;*/
      x = coords[n];
    }
    if (coordtype[n] & VT_LATITUDE) {
      gy = VAR_2D(&vars[coordids[n]])[0];
      /*if (nc == v->nd) nj = df->dimensions[coordids[n]].size;*/
      y = coords[n];
    }
    if (coordtype[n] & VT_Z) {
      df_variable_t *vc = &vars[coordids[n]];
      memcpy(zgrid, VAR_1D(&vars[coordids[n]])[0], nz * sizeof(double));
      /*if (nc == v->nd) nz = df->dimensions[coordids[n]].size;*/
      if (vc->type & VT_SIGMA) sigma = 1; 
      z = coords[n];
    }
  }

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Orient the vertical coordinate */
  /*memcpy(zgrid, VAR_1D(&vars[coordids[0]])[0], nz * sizeof(double));*/
  if (!sigma && vars[coordids[0]].z_is_depth)
      for (k = 0; k < nz; k++) zgrid[k] *= -1.0;
  zs = 0; zb = nz - 1;
  if (zgrid[0] < zgrid[nz - 1]) {
    zs = nz - 1;
    zb = 0;
  }

  /* Create the grid spec function if required */
  if (v->gs0_master_3d == NULL && v->gs1_master_3d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 3)
      quit("df_eval_coords : interp_us_3d : nc is not 3 but %d!\n", nc);

    /* Allocate the private data */
    v->nz = nz;
    v->kn = i_alloc_1d(nz);
    v->kflag = i_alloc_2d(2, nz);
    memset(v->kn, 0, nz * sizeof(int));

    /* Allocate and set the points */
    dsize = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[zs][j][i] == v->missing)) continue;
	if ((hasfil && vd[zs][j][i] == v->fillvalue)) continue;
	if (isnan(vd[zs][j][i])) continue;
	dsize++;
      }
    }

    v->kmapi = i_alloc_2d(dsize, nz);
    v->kmapj = i_alloc_2d(dsize, nz);
    p = (point **)calloc(nz, sizeof(point *));
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ni; ++i) {
	for (j = 0; j < nj; ++j) {
	  if ((hasmis && vd[k][j][i] == v->missing)) continue;
	  if ((hasfil && vd[k][j][i] == v->fillvalue)) continue;
	  if (isnan(vd[k][j][i])) continue;
	  if (p[k] == NULL)
	    p[k] = (point *)alloc_1d(dsize, sizeof(point));
	  v->kmapi[k][v->kn[k]] = i;
	  v->kmapj[k][v->kn[k]] = j;
	  p[k][v->kn[k]].x = gx[j][i];
	  p[k][v->kn[k]++].y = gy[j][i];
	}
      }
    }

    /* Allocate and set the Delaunay triangulation */
    d = (delaunay **)calloc(nz+1, sizeof(delaunay *));
    for (k = 0; k < nz; k++) {
      d[k] = NULL;
      if (v->kn[k] == 0) continue;
      d[k] = delaunay_build(v->kn[k], p[k], 0, NULL, 0, NULL);
      d[k]->vid = 0;
      d[k]->ptf = 1;

      for (n = 0; n < v->kn[k]; n++) {
	i = v->kmapi[k][n];
	j = v->kmapj[k][n];
	d[k]->points[n].v = d_alloc_1d(1);
	d[k]->points[n].z = d[k]->points[n].v[0] = vd[k][j][i];
      }
    }

    /* Finished with temp array */
    free(p);

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs0_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    gs1 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs1_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    for (k = 0; k < nz; k++) {

      v->gs0_master_3d[k] = v->gs1_master_3d[k] = NULL;
      if (d[k] == NULL) continue;

      gs0[k] = grid_spec_create();
      gs1[k] = grid_spec_create();

      grid_interp_init_t(gs0[k], d[k], i_rule, 0);
      gs0[k]->d = d[k];
      v->gs0_master_3d[k] = gs0[k];

      grid_interp_init_t(gs1[k], d[k], i_rule, 0);
      gs1[k]->d = d[k];
      v->gs1_master_3d[k] = gs1[k];
    }
    v->rv[0] = v->rv[1] = -9999;
  }

  /* Find the k layer */
  if (z >= zgrid[zs]) {
    kindex = (double)zs;
  } else if (z <= zgrid[zb]) {
    kindex = (double)zb;
  } else {
    for (k = 0; k < nz - 1; ++k) {
      if (zb == 0) {
	if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
	  kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
	  break;
	}
      } else {
	if ((z <= zgrid[k]) && (z > zgrid[k + 1])) {
	  kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
	  break;
	}
      }
    }
  }
  ks = (int)kindex;
  if (zb == 0)
    kt = min(ks + 1, zs);
  else
    kt = min(ks + 1, zb);

  /* No-gradient below bottom layer */
  if (!v->kn[ks]) ks = kt;
  kindex = kindex - ks;

  /* Reinitialize the weights if the data has changed */
  if (record != v->rv[rid]) {
    for (k = 0; k < nz; k++) v->kflag[rid][k] = 0;
    v->rv[rid] = record;
  }

  /* Rebuild the weights for the required layers */
  gs = (rid == 0) ? v->gs0_master_3d : v->gs1_master_3d;
  for (k = ks; k <= kt; ++k) {
    if (!v->kflag[rid][k]) {
      delaunay *d = gs[k]->d;
      for (n = 0; n < v->kn[k]; ++n) {
	i = v->kmapi[k][n];
	j = v->kmapj[k][n];
	d->points[n].z = vd[k][j][i];
      }
      v->kflag[rid][k] = 1;

      if (strcmp(i_rule, "nn_sibson") == 0)
	for (n = 0; n < v->kn[k]; n++)
	  gs[k]->rebuild(gs[k]->interpolator, &d->points[n]);
      else
	gs[k]->rebuild(gs[k]->interpolator, d->points);
    }
  }

  /* Do the horizontal interpolation on layers bracketing the depth */
  for (k = ks; k <= kt; ++k) {
    tdata[k - ks] = grid_interp_on_point(gs[k], x, y);
  }

  /* Do the vertical interpolation */
  if (isnan(tdata[0]))
    val = tdata[1];
  else if (isnan(tdata[1]))
    val = tdata[0];
  else 
    val = (1 - kindex) * tdata[0] + (kindex) * tdata[1];

  return(val);
}


#define Z_PS 0x001
#define Z_PB 0x002
#define Z_NS 0x004
#define Z_NB 0x008

/*------------------------------------------------------------------*/
/* Unstructured interpolation for VT_INFERRED files                 */
/*------------------------------------------------------------------*/
double interp_us_3d_i(datafile_t *df, df_variable_t *v, int record,
		    double coords[])
{
  int nz = df->dimensions[v->dimids[0]].size;
  int nj = df->dimensions[v->dimids[1]].size;
  int ni = df->dimensions[v->dimids[2]].size;
  int nc = df_get_num_coords(df, v);
  int dsize;
  int *coordids = df_get_coord_ids(df, v);
  int *coordtype = df_get_coord_types(df, v);
  double ***vd = VAR_3D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  GRID_SPECS **gs0 = NULL;
  GRID_SPECS **gs1 = NULL;
  int n, i, j, k;
  point **p;
  delaunay **d;
  double val;
  double x, y, z;
  double kindex   = 0.0;
  double  tdata[2] = {NaN, NaN};
  int ks, kt;
  int rid = v->rid;
  GRID_SPECS **gs;
  int hasmis = (v->missing) ? 1 : 0;
  int hasfil = (v->fillvalue) ? 1 : 0;
  char i_rule[MAXSTRLEN];
  double *gx, *gy, *zgrid;
  int zc, zs, zb;
  double sgn = 1.0;

  /* Find which coordinates are lat, long, depth */
  for (n = 0; n < nc; n++) {
    if (coordtype[n] & VT_LONGITUDE) {
      gx = VAR_1D(&vars[coordids[n]])[0];
      /*if (nc == v->nd) ni = df->dimensions[coordids[n]].size;*/
      x = coords[n];
    }
    if (coordtype[n] & VT_LATITUDE) {
      gy = VAR_1D(&vars[coordids[n]])[0];
      /*if (nc == v->nd) nj = df->dimensions[coordids[n]].size;*/
      y = coords[n];
    }
    if (coordtype[n] & VT_Z) {
      memcpy(zgrid, VAR_1D(&vars[coordids[n]])[0], nz * sizeof(double));
      /*if (nc == v->nd) nz = df->dimensions[coordids[n]].size;*/
      z = coords[n];
    }
  }

  /* Figure out the direction of the vertical coordinate */
  if (zgrid[0] >= 0.0 && zgrid[nz-1] >= 0.0 && zgrid[0] < zgrid[nz-1]) 
    zc = Z_PS;
  if (zgrid[0] >= 0.0 && zgrid[nz-1] >= 0.0 && zgrid[0] > zgrid[nz-1]) 
    zc = Z_PB;
  if (zgrid[0] <= 0.0 && zgrid[nz-1] <= 0.0 && zgrid[0] < zgrid[nz-1]) 
    zc = Z_NB;
  if (zgrid[0] <= 0.0 && zgrid[nz-1] <= 0.0 && zgrid[0] > zgrid[nz-1]) 
    zc = Z_NS;
  zs = 0; zb = nz - 1;
  if (zc & (Z_NB|Z_PB)) {
    zs = nz - 1;
    zb = 0;
  }

  if (strlen(df->i_rule))
    strcpy(i_rule, df->i_rule);
  else
    strcpy(i_rule, "linear");

  /* Create the grid spec function if required */
  if (v->gs0_master_3d == NULL && v->gs1_master_3d == NULL ) {
    /*
     * All this hash mechanism relies on the fact that we're dealing
     * with 3D data, so I'm putting in an explicit check just to make
     * sure
     */
    if (nc != 3)
      quit("df_eval_coords : interp_us_3d : nc is not 3 but %d!\n", nc);

    /* Allocate the private data */
    v->nz = nz;
    v->kn = i_alloc_1d(nz);
    v->kflag = i_alloc_2d(2, nz);
    memset(v->kn, 0, nz * sizeof(int));

    /* Allocate and set the points */
    dsize = 0;
    for (i = 0; i < ni; ++i) {
      for (j = 0; j < nj; ++j) {
	if ((hasmis && vd[zs][j][i] == v->missing)) continue;
	if ((hasfil && vd[zs][j][i] == v->fillvalue)) continue;
	if (isnan(vd[zs][j][i])) continue;
	dsize++;
      }
    }

    v->kmapi = i_alloc_2d(dsize, nz);
    v->kmapj = i_alloc_2d(dsize, nz);
    p = (point **)calloc(nz, sizeof(point *));
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ni; ++i) {
	for (j = 0; j < nj; ++j) {
	  if ((hasmis && vd[k][j][i] == v->missing)) continue;
	  if ((hasfil && vd[k][j][i] == v->fillvalue)) continue;
	  if (isnan(vd[k][j][i])) continue;
	  if (p[k] == NULL)
	    p[k] = (point *)alloc_1d(dsize, sizeof(point));
	  v->kmapi[k][v->kn[k]] = i;
	  v->kmapj[k][v->kn[k]] = j;
	  p[k][v->kn[k]].x = gx[i];
	  p[k][v->kn[k]++].y = gy[j];
	}
      }
    }

    /* Allocate and set the Delaunay triangulation */
    d = (delaunay **)calloc(nz+1, sizeof(delaunay *));
    for (k = 0; k < nz; k++) {
      d[k] = NULL;
      if (v->kn[k] == 0) continue;
      d[k] = delaunay_build(v->kn[k], p[k], 0, NULL, 0, NULL);
      d[k]->vid = 0;
      d[k]->ptf = 1;
      for (n = 0; n < v->kn[k]; n++) {
	i = v->kmapi[k][n];
	j = v->kmapj[k][n];
	d[k]->points[n].v = d_alloc_1d(1);
	d[k]->points[n].z = d[k]->points[n].v[0] = vd[k][j][i];
      }
    }

    /* Finished with temp array */
    free(p);

    /* Allocate and set the grid spec structure. Note there are */
    /* two grid_spec structures for the two time levels that */
    /* bracket the target time. These both point to the same */
    /* Delaunay structure, but contain different weights due to */
    /* the different data at the time levels. */
    gs0 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs0_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    gs1 = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    v->gs1_master_3d = (GRID_SPECS **)calloc(nz, sizeof(GRID_SPECS *));
    for (k = 0; k < nz; k++) {

      v->gs0_master_3d[k] = v->gs1_master_3d[k] = NULL;
      if (d[k] == NULL) continue;

      gs0[k] = grid_spec_create();
      gs1[k] = grid_spec_create();

      grid_interp_init_t(gs0[k], d[k], i_rule, 0);
      gs0[k]->d = d[k];
      v->gs0_master_3d[k] = gs0[k];

      grid_interp_init_t(gs1[k], d[k], i_rule, 0);
      gs1[k]->d = d[k];
      v->gs1_master_3d[k] = gs1[k];
    }
    v->rv[0] = v->rv[1] = -9999;
  }

  /* Find the k layer */      
  if (zc & (Z_NS|Z_NB)) sgn = -1.0;

  if (fabs(z) < sgn * zgrid[zs]) {
    kindex = (double)zs;
  } else if (fabs(z) >= sgn * zgrid[zb]) {
    kindex = (double)zb;
  } else {
    for (k = 0; k < nz - 1; ++k) {
      if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
	kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
	break;
      }
    }
  }
  ks = (int)kindex;
  if (zb == 0)
    kt = min(ks + 1, zs);
  else
    kt = min(ks + 1, zb);

  /* No-gradient below bottom layer */
  if (!v->kn[ks]) ks = kt;
  kindex = kindex - ks;

  /* Reinitialize the weights if the data has changed */
  if (record != v->rv[rid]) {
    for (k = 0; k < nz; k++) v->kflag[rid][k] = 0;
    v->rv[rid] = record;
  }

  /* Rebuild the weights for the required layers */
  gs = (rid == 0) ? v->gs0_master_3d : v->gs1_master_3d;
  for (k = ks; k <= kt; ++k) {
    if (!v->kflag[rid][k]) {
      delaunay *d = gs[k]->d;
      for (n = 0; n < v->kn[k]; ++n) {
	i = v->kmapi[k][n];
	j = v->kmapj[k][n];
	d->points[n].z = vd[k][j][i];
      }
      v->kflag[rid][k] = 1;

      if (strcmp(i_rule, "nn_sibson") == 0)
	for (n = 0; n < v->kn[k]; n++)
	  gs[k]->rebuild(gs[k]->interpolator, &d->points[n]);
      else
	gs[k]->rebuild(gs[k]->interpolator, d->points);
    }
  }

  /* Do the horizontal interpolation on layers bracketing the depth */
  for (k = ks; k <= kt; ++k) {
    tdata[k - ks] = grid_interp_on_point(gs[k], x, y);
  }

  /* Do the vertical interpolation */
  if (isnan(tdata[0]))
    val = tdata[1];
  else if (isnan(tdata[1]))
    val = tdata[0];
  else 
    val = (1 - kindex) * tdata[0] + (kindex) * tdata[1];

  return(val);
}


/*
 * Nearest neighbour interpolation within limits
 *
 * Notes:
 *  o) 2D only so will return the same number for every depth
 *  o) Can optimise performance by introducing a hashing mechanism
 */
double interp_nearest_within_eps(datafile_t *df, df_variable_t *v, int record,
				 double coords[])
{
  int dsize = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  double *d = VAR_1D(v)[record - v->start_record];
  df_variable_t *vars = df->variables;
  int i, j;
  float idist;
  double d_eps = 0.0005; /* in degrees ~ 50m */
  
  /*
   * Loop through all points and latch on to the the first one that
   * falls within the limit
   */

  for (i = 0; i < dsize; ++i) {
    double dist = 0.0;

    for (j = 0; j < nc; ++j) {
      double dx = VAR_1D(&vars[coordids[j]])[0][i] - coords[j];
      dist += dx * dx;
    }

    dist = sqrt(dist);

    if (dist <= d_eps)
      return(d[i]);
  }

  /* Nothing found - missing value?? */
  /*return(NaN);*/
  return(0.0);
}


double interp2d_nearest_within_eps(datafile_t *df, df_variable_t *v, int record,
				 double coords[])
{
  int nj = df->dimensions[v->dimids[0]].size;
  int ni = df->dimensions[v->dimids[1]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  df_variable_t *vars = df->variables;
  double **d = VAR_2D(v)[record - v->start_record];
  double dx, dy;
  int i, j, n;
  double d_eps = 0.00001;
  double **gx = VAR_2D(&vars[coordids[0]])[0];
  double **gy = VAR_2D(&vars[coordids[1]])[0];

  /*
   * Loop through all points and latch on to the the first one that
   * falls within the limit
   */

  for (j = 0; j < nj; ++j) {
    for (i = 0; i < ni; ++i) {
      double dist = 0.0;
      double x = gx[j][i];
      double y = gy[j][i];
      if (isnan(x) || isnan(y)) continue;

      dx = x - coords[0];
      dy = y - coords[1];
      dist = sqrt(dx * dx + dy * dy);
      if (dist <= d_eps) {
	return(d[j][i]);
      }

    }
  }

  /* Nothing found - missing value?? */
  /*return(NaN);*/
  return(0.0);
}

double interp3d_nearest_within_eps(datafile_t *df, df_variable_t *v, int record,
				 double coords[])
{
  int nz = df->dimensions[v->dimids[0]].size;
  int nj = df->dimensions[v->dimids[1]].size;
  int ni = df->dimensions[v->dimids[2]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  df_variable_t *vars = df->variables;
  double ***d = VAR_3D(v)[record - v->start_record];
  double *zgrid = VAR_1D(&vars[coordids[nc - 1]])[0];
  double dx, dy, z = coords[nc - 1];
  int i, j, k, n;
  double d_eps = 0.00001;
  double **gx = VAR_2D(&vars[coordids[0]])[0];
  double **gy = VAR_2D(&vars[coordids[1]])[0];

  /*
   * Loop through all points and latch on to the the first one that
   * falls within the limit
   */
  /*
  /* Find the k layer */      
  if (z < zgrid[0]) {
    k = 0;
  } else if (z >= zgrid[nz - 1]) {
    k = nz - 1;
  } else {
    for (k = 0; k < nz - 1; ++k) {
      if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
	break;
      }
    }
  }

  for (j = 0; j < nj; ++j) {
    for (i = 0; i < ni; ++i) {
      double dist = 0.0;
      double x = gx[j][i];
      double y = gy[j][i];
      if (isnan(x) || isnan(y)) continue;

      dx = x - coords[1];
      dy = y - coords[2];
      dist = sqrt(dx * dx + dy * dy);
      if (dist <= d_eps)
	return(d[k][j][i]);
    }
  }

  /* Nothing found - missing value?? */
  /*return(NaN);*/
  return(0.0);
}

/* Nearest neighbour for 2D point array files */
double interp2d_nearest(datafile_t *df, df_variable_t *v, int record,
			    double coords[])
{
  int ni = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);
  int *coordids = df_get_coord_ids(df, v);
  df_variable_t *vars = df->variables;
  double *d = VAR_1D(v)[record - v->start_record];
  double dx, dy;
  int i, mi, n;
  double dmin = HUGE;
  double *gx = VAR_1D(&vars[coordids[0]])[0];
  double *gy = VAR_1D(&vars[coordids[1]])[0];

  /*
   * Loop through all points and latch on to the the first one that
   * falls within the limit
   */

  for (i = 0; i < ni; ++i) {
    double dist = 0.0;
    double x = gx[i];
    double y = gy[i];
    if (isnan(x) || isnan(y)) continue;

    dx = x - coords[0];
    dy = y - coords[1];
    dist = sqrt(dx * dx + dy * dy);
    if (dist <= dmin) {
      dmin = dist;
      mi = i;
    }
  }
  return(d[mi]);
}


/* Nearest neighbour for 3D point array files */
/* Formulated for ROMS station data: note data d[ni][nz] */
/* rather than the usual d[nz][ni]. */
double interp3d_nearest(datafile_t *df, df_variable_t *v, int record,
			  double coords[])
{
  int ni = df->dimensions[v->dimids[0]].size;
  int nc = df_get_num_coords(df, v);
  int nz = df->dimensions[v->dimids[nc-1]].size;
  int *coordids = df_get_coord_ids(df, v);
  df_variable_t *vars = df->variables;
  double **d = VAR_2D(v)[record - v->start_record];
  double *zgrid = VAR_1D(&vars[coordids[nc - 1]])[0];
  double val, dx, dy, z = coords[nc - 1];
  int i, mi, k, ks, kt, n;
  double kindex, dmin = HUGE;
  double x = coords[0];
  double y = coords[1];
  double *gx = VAR_1D(&vars[coordids[0]])[0];
  double *gy = VAR_1D(&vars[coordids[1]])[0];

  /*
   * Loop through all points and latch on to the the first one that
   * falls within the limit
   */
  /*
  /* Locate the k layer */
  if (z < zgrid[0]) {
    kindex = 0.0;
  } else if (z >= zgrid[nz - 1]) {
    kindex = nz - 1;
  } else {
    for (k = 0; k < nz - 1; ++k) {
      if ((z >= zgrid[k]) && (z < zgrid[k + 1])) {
        kindex = k + (z - zgrid[k]) / (zgrid[k + 1] - zgrid[k]);
        break;
      }
    }
  }
  ks = (int)kindex;
  kt = min(ks + 1, nz - 1);
  kindex = kindex - ks;

  for (i = 0; i < ni; ++i) {
    double dist = 0.0;
    dx = gx[i];
    dy = gy[i];
    if (isnan(x) || isnan(y)) continue;

    dx -= x;
    dy -= y;
    dist = sqrt(dx * dx + dy * dy);
    if (dist <= dmin) {
      dmin = dist;
      mi = i;
    }
  }
  val = (1 - kindex) * d[mi][ks] + (kindex) * d[mi][kt];
  return(val);
}
  
/*
Routine to interpolate grid spatial data using a simple bilinear-like
interpolation scheme.
*/
double interp_linear(datafile_t *df, df_variable_t *v, int record,
                     double coords[])
{
  int i, j;
  double val = 0.0;
  double indices[MAXNUMDIMS];
  int nd = df_get_num_dims(df, v);
  int *dimids = df_get_dim_ids(df, v);

  if (df_ctoi(df, v, coords, indices) == nd) {
    double findices[MAXNUMDIMS];
    int iindices[MAXNUMDIMS];
    int corner[MAXNUMDIMS];
    int coffsets[MAXNUMDIMS * MAXNUMDIMS - 1][MAXNUMDIMS];
    int dsize[MAXNUMDIMS];
    int ncorners = 1 << nd;

#if 0
    fprintf(stderr, "Coords =");
    for (i = 0; i < v->csystem->nc; ++i)
      fprintf(stderr, " %g", coords[i]);
    fprintf(stderr, "\n");
    
    fprintf(stderr, "Indices =");
    for (i = 0; i < nd; ++i)
      fprintf(stderr, " %g", indices[i]);
    fprintf(stderr, "\n");
#endif


    /* Compute the fraction indices and integer indicies. Trim at
       boundaries if necessary. */
    for (i = 0; i < nd; ++i) {
      dsize[i] = df->dimensions[dimids[i]].size - 1;
      findices[i] = indices[i];
      if (findices[i] < 0)
        findices[i] = 0.0;
      if (findices[i] > dsize[i])
        findices[i] = (double)dsize[i];
      iindices[i] = (int)floor(findices[i]);
      findices[i] -= iindices[i];
    }
    
    /* Create the corner indice offsets (0 or 1) */
    for (i = 0; i < ncorners; ++i)
      for (j = 0; j < nd; ++j)
        coffsets[i][j] = (i >> j) & 0x01;
    
    /* Now step though all of the corners and perform the linear
       interpolation */
    for (j = 0; j < ncorners; ++j) {
      double term = 1.0;
      for (i = 0; i < nd; ++i) {
        corner[i] = iindices[i] + coffsets[j][i];
        if (corner[i] > dsize[i])
          corner[i] = dsize[i];
        if (coffsets[j][i])
          term *= findices[i];
        else
          term *= (1 - findices[i]);
      }
      
#if 0
      fprintf(stderr, "Corner(%d) =", j);
      for (i = 0; i < nd; ++i)
        fprintf(stderr, " %d", corner[i]);
      fprintf(stderr, "\n");
#endif
      
      /* Round term if very close to a corner. */
      if (fabs(term) < 1e-5)
        term = 0.0;
      else if (fabs(1.0 - term) < 1e-5)
        term = 1.0;
      
      if (term > 0.0)
        val += term * df_get_data_value(df, v, record, corner);
    }
  }
  return (val);
}


/*
MH: Routine to interpolate grid spatial data using a simple bilinear-like
interpolation scheme for variables with a range 0-360.
Values separated by more than 180 (e.g. across the transition 0/360), are 
wrapped (i.e. 360 is added / subtracted to one value).
*/
double interp_linear_degrees(datafile_t *df, df_variable_t *v, int record,
			     double coords[])
{
  int i, j;
  double val = 0.0;
  double indices[MAXNUMDIMS];
  int nd = df_get_num_dims(df, v);
  int *dimids = df_get_dim_ids(df, v);
  double pi2 = 360.0;
 
  if (df_ctoi(df, v, coords, indices) == nd) {
    double findices[MAXNUMDIMS];
    int iindices[MAXNUMDIMS];
    int corner[MAXNUMDIMS];
    int coffsets[MAXNUMDIMS * MAXNUMDIMS - 1][MAXNUMDIMS];
    int dsize[MAXNUMDIMS];
    int ncorners = 1 << nd;
    double dval[ncorners];
    double term[ncorners];

#if 0
    fprintf(stderr, "Coords =");
    for (i = 0; i < v->csystem->nc; ++i)
      fprintf(stderr, " %g", coords[i]);
    fprintf(stderr, "\n");
    
    fprintf(stderr, "Indices =");
    for (i = 0; i < nd; ++i)
      fprintf(stderr, " %g", indices[i]);
    fprintf(stderr, "\n");
#endif


    /* Compute the fraction indices and integer indicies. Trim at
       boundaries if necessary. */
    for (i = 0; i < nd; ++i) {
      dsize[i] = df->dimensions[dimids[i]].size - 1;
      findices[i] = indices[i];
      if (findices[i] < 0)
        findices[i] = 0.0;
      if (findices[i] > dsize[i])
        findices[i] = (double)dsize[i];
      iindices[i] = (int)floor(findices[i]);
      findices[i] -= iindices[i];
    }
    
    /* Create the corner indice offsets (0 or 1) */
    for (i = 0; i < ncorners; ++i)
      for (j = 0; j < nd; ++j)
        coffsets[i][j] = (i >> j) & 0x01;

    /* Now step though all of the corners and perform the linear
       interpolation */
    for (j = 0; j < ncorners; ++j) {
      term[j] = 1.0;
      for (i = 0; i < nd; ++i) {
        corner[i] = iindices[i] + coffsets[j][i];
        if (corner[i] > dsize[i])
          corner[i] = dsize[i];
        if (coffsets[j][i])
          term[j] *= findices[i];
        else
          term[j] *= (1 - findices[i]);
      }
      
#if 0
      fprintf(stderr, "Corner(%d) =", j);
      for (i = 0; i < nd; ++i)
        fprintf(stderr, " %d", corner[i]);
      fprintf(stderr, "\n");
#endif
      
      /* Round term if very close to a corner. */
      if (fabs(term[j]) < 1e-5)
        term[j] = 0.0;
      else if (fabs(1.0 - term[j]) < 1e-5)
        term[j] = 1.0;
    }
    for (j = 0; j < ncorners; ++j) {
      dval[j] = df_get_data_value(df, v, record, corner);
      if (j > 0 && dval[0] - dval[j] > 180.0) dval[j] += pi2;
      if (j > 0 && dval[0] - dval[j] < -180.0) dval[j] -= pi2;
    } 
    for (j = 0; j < ncorners; ++j) {
      if (term[j] > 0.0)
        val += term[j] * dval[j];
    }
    val = (val < 0.0) ? val + pi2 : val;
    val = fmod(val, pi2);
  }
  return (val);
}


#define CF  1  /* Cloud flag    */
#define LF  2  /* Land flag     */
#define FV  4  /* Fill value    */
#define MV  8  /* Missing value */


/*
MH: Routine to interpolate grid spatial data using a simple bilinear-like
interpolation scheme. Sets a no-gradient across land and doesn't interpolate
over cloud values.
*/
double interp_linear_flagged(datafile_t *df, df_variable_t *v, int record,
			     double coords[])
{
  int i, j, ii=0;
  double val = 0.0, dval;
  double indices[MAXNUMDIMS];
  int nd = df_get_num_dims(df, v);
  int *dimids = df_get_dim_ids(df, v);
  int ci, cj;

  if (df_ctoi(df, v, coords, indices) == nd) {
    double findices[MAXNUMDIMS];
    int iindices[MAXNUMDIMS];
    int corner[MAXNUMDIMS];
    int coffsets[MAXNUMDIMS * MAXNUMDIMS - 1][MAXNUMDIMS];
    int dsize[MAXNUMDIMS];
    int ncorners = 1 << nd;
    double term[ncorners];
    double dval[ncorners];
    int mask[ncorners];

#if 0
    fprintf(stderr, "Coords =");
    for (i = 0; i < v->csystem->nc; ++i)
      fprintf(stderr, " %g", coords[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "Indices =");
    for (i = 0; i < nd; ++i)
      fprintf(stderr, " %g", indices[i]);
    fprintf(stderr, "\n");
#endif


    /* Compute the fraction indices and integer indicies. Trim at
       boundaries if necessary. */
    for (i = 0; i < nd; ++i) {
      dsize[i] = df->dimensions[dimids[i]].size - 1;
      findices[i] = indices[i];
      if (findices[i] < 0)
        findices[i] = 0.0;
      if (findices[i] > dsize[i])
        findices[i] = (double)dsize[i];
      iindices[i] = (int)floor(findices[i]);
      findices[i] -= iindices[i];
    }

    /* Create the corner indice offsets (0 or 1) */
    for (i = 0; i < ncorners; ++i)
      for (j = 0; j < nd; ++j)
        coffsets[i][j] = (i >> j) & 0x01;

    /* Now step though all of the corners and perform the linear
       interpolation */
    for (j = 0; j < ncorners; ++j) {
      mask[j] = 0;
      term[j] = 1.0;
      for (i = 0; i < nd; ++i) {
        corner[i] = iindices[i] + coffsets[j][i];
        if (corner[i] > dsize[i])
          corner[i] = dsize[i];
        if (coffsets[j][i])
          term[j] *= findices[i];
        else
          term[j] *= (1 - findices[i]);
      }
      dval[j] = df_get_data_value(df, v, record, corner);
      if (dval[j] == v->lflag) mask[j] |= LF; 
      if (dval[j] == v->cflag) mask[j] |= CF; 
      if (dval[j] == v->missing) mask[j] |= MV;
      if (dval[j] == v->fillvalue) mask[j] |= FV;

#if 0
      fprintf(stderr, "Corner(%d) =", j);
      for (i = 0; i < nd; ++i)
        fprintf(stderr, " %d (%d)", corner[i], mask[i]);
      fprintf(stderr, "\n");
#endif
    }

    /* Find the first valid value */
    for (j = 0; j < ncorners; ++j)
      if (!mask[j]) break;
    /* Fill the non-valid corners up to this corner */
    if (j < 4) {
      for (i = 0; i < j; ++i)
	dval[i] = dval[j];
      /* Fill subsequent non-valid corners */
      for (; j < ncorners; ++j) {
	for (i = 0; i < ncorners, i != j; ++i) 
	  if (!mask[i]) dval[j] = dval[i];
      }
    } else {
      /* Set all non-valid values the same */
      for (j = 0; j < ncorners; ++j) {
	/* If cloud is found, set all the corners to cloud */
	if (mask[j] & CF) {
	  for (i = 0; i < ncorners, i != j; ++i)
	    dval[j] = v->cflag;
	  break;
	}
      }
      /* Otherwise, set all corners to land */
      if (dval[0] != v->cflag)
	for (j = 0; j < ncorners; ++j)
	  dval[j] = v->lflag;
    }

    for (j = 0; j < ncorners; ++j) {
      /* Round term if very close to a corner. */
      if (fabs(term[j]) < 1e-5)
        term[j] = 0.0;
      else if (fabs(1.0 - term[j]) < 1e-5)
        term[j] = 1.0;

      if (term[j] > 0.0) {
        val += term[j] * dval[j];
      }
    }
  }

  return (val);
}


/*
MH: Routine to interpolate grid spatial data using a simple bilinear-like
interpolation scheme. Fills missing values with the nearest valid value.
*/
double interp_linear_filled(datafile_t *df, df_variable_t *v, int record,
			    double coords[])
{
  int i, j, ii=0;
  double val = 0.0, dval;
  double indices[MAXNUMDIMS];
  int nd = df_get_num_dims(df, v);
  int *dimids = df_get_dim_ids(df, v);

  if (df_ctoi(df, v, coords, indices) == nd) {
    double findices[MAXNUMDIMS];
    int iindices[MAXNUMDIMS];
    int corner[MAXNUMDIMS];
    int coffsets[MAXNUMDIMS * MAXNUMDIMS - 1][MAXNUMDIMS];
    int dsize[MAXNUMDIMS];
    int ncorners = 1 << nd;
    double term[ncorners];
    double dval[ncorners];
    int mask[ncorners];

#if 0
    fprintf(stderr, "Coords =");
    for (i = 0; i < v->csystem->nc; ++i)
      fprintf(stderr, " %g", coords[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "Indices =");
    for (i = 0; i < nd; ++i)
      fprintf(stderr, " %g", indices[i]);
    fprintf(stderr, "\n");
#endif


    /* Compute the fraction indices and integer indicies. Trim at
       boundaries if necessary. */
    for (i = 0; i < nd; ++i) {
      dsize[i] = df->dimensions[dimids[i]].size - 1;
      findices[i] = indices[i];
      if (findices[i] < 0)
        findices[i] = 0.0;
      if (findices[i] > dsize[i])
        findices[i] = (double)dsize[i];
      iindices[i] = (int)floor(findices[i]);
      findices[i] -= iindices[i];
    }

    /* Create the corner indice offsets (0 or 1) */
    for (i = 0; i < ncorners; ++i)
      for (j = 0; j < nd; ++j)
        coffsets[i][j] = (i >> j) & 0x01;

    /* Now step though all of the corners and perform the linear
       interpolation */
    for (j = 0; j < ncorners; ++j) {
      mask[j] = 0;
      term[j] = 1.0;
      for (i = 0; i < nd; ++i) {
        corner[i] = iindices[i] + coffsets[j][i];
        if (corner[i] > dsize[i])
          corner[i] = dsize[i];
        if (coffsets[j][i])
          term[j] *= findices[i];
        else
          term[j] *= (1 - findices[i]);
      }
      dval[j] = df_get_data_value(df, v, record, corner);
      if (v->type & VT_MV && dval[j] == (double)v->missing) {
	mask[j] |= MV;
	dval[j] = find_close(df, v, record, corner);
      }
      if (v->type & VT_FV && dval[j] == (double)v->fillvalue) {
	mask[j] |= FV;
	dval[j] = find_close(df, v, record, corner);
      }

#if 0
      fprintf(stderr, "Corner(%d) =", j);
      for (i = 0; i < nd; ++i)
        fprintf(stderr, " %d (%d)", corner[i], mask[i]);
      fprintf(stderr, "\n");
#endif
    }

    for (j = 0; j < ncorners; ++j) {
      /* Round term if very close to a corner. */
      if (fabs(term[j]) < 1e-5)
        term[j] = 0.0;
      else if (fabs(1.0 - term[j]) < 1e-5)
        term[j] = 1.0;

      if (term[j] > 0.0) {
        val += term[j] * dval[j];
      }
    }
  }
  return (val);
}

#define LANDCELL 99
/*
MH: Routine to interpolate grid spatial data using a simple bilinear-like
interpolation scheme. Fills missing or land (99) values with the nearest 
valid value.
*/
double interp_linear_bathy(datafile_t *df, df_variable_t *v, int record,
			    double coords[])
{
  int i, j, ii=0;
  double val = (v->type & VT_BATHY) ? NaN : 0.0;
  double dval;
  double indices[MAXNUMDIMS];
  int nd = df_get_num_dims(df, v);
  int *dimids = df_get_dim_ids(df, v);

  if (df_ctoi(df, v, coords, indices) == nd) {
    double findices[MAXNUMDIMS];
    int iindices[MAXNUMDIMS];
    int corner[MAXNUMDIMS];
    int coffsets[MAXNUMDIMS * MAXNUMDIMS - 1][MAXNUMDIMS];
    int dsize[MAXNUMDIMS];
    int ncorners = 1 << nd;
    double term[ncorners];
    double dval[ncorners];
    int mask[ncorners];
    val = 0.0;
#if 0
    fprintf(stderr, "Coords =");
    for (i = 0; i < v->csystem->nc; ++i)
      fprintf(stderr, " %g", coords[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "Indices =");
    for (i = 0; i < nd; ++i)
      fprintf(stderr, " %g", indices[i]);
    fprintf(stderr, "\n");
#endif


    /* Compute the fraction indices and integer indicies. Trim at
       boundaries if necessary. */
    for (i = 0; i < nd; ++i) {
      dsize[i] = df->dimensions[dimids[i]].size - 1;
      findices[i] = indices[i];
      if (findices[i] < 0)
        findices[i] = 0.0;
      if (findices[i] > dsize[i])
        findices[i] = (double)dsize[i];
      iindices[i] = (int)floor(findices[i]);
      findices[i] -= iindices[i];
    }

    /* Create the corner indice offsets (0 or 1) */
    for (i = 0; i < ncorners; ++i)
      for (j = 0; j < nd; ++j)
        coffsets[i][j] = (i >> j) & 0x01;

    /* Now step though all of the corners and perform the linear
       interpolation */
    for (j = 0; j < ncorners; ++j) {
      mask[j] = 0;
      term[j] = 1.0;
      for (i = 0; i < nd; ++i) {
        corner[i] = iindices[i] + coffsets[j][i];
        if (corner[i] > dsize[i])
          corner[i] = dsize[i];
        if (coffsets[j][i])
          term[j] *= findices[i];
        else
          term[j] *= (1 - findices[i]);
      }
      dval[j] = df_get_data_value(df, v, record, corner);
      if (isnan(dval[j]) || dval[j] == LANDCELL) {
	mask[j] |= MV;
	dval[j] = find_close_bathy(df, v, record, corner);
      }

#if 0
      fprintf(stderr, "Corner(%d) =", j);
      for (i = 0; i < nd; ++i)
        fprintf(stderr, " %d (%d)", corner[i], mask[i]);
      fprintf(stderr, "\n");
#endif
    }

    for (j = 0; j < ncorners; ++j) {
      /* Round term if very close to a corner. */
      if (fabs(term[j]) < 1e-5)
        term[j] = 0.0;
      else if (fabs(1.0 - term[j]) < 1e-5)
        term[j] = 1.0;

      if (term[j] > 0.0) {
        val += term[j] * dval[j];
      }
    }
  }
  
  return (val);
}


/* Locate the closest cell to is[1], is[0] that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
static double find_close(datafile_t *df, df_variable_t *v, int r, int *is)
                 
{
  int i, j, ci, cj;
  double mindist;
  int level;
  int xlim, ylim;
  int ic, jc, kc;
  double val, nval;

  /* Get the indexes */
  if (v->nd == 2) {
    kc = 0;
    jc = is[0];
    ic = is[1];
  } else if (v->nd == 3) {
    kc = is[0];
    jc = is[1];
    ic = is[2];
  } else {
    quit("Landfill only possible for 2D or 3D datafiles.\n");
  }

  /* Get the horizontal grid size */
  for (i = 0; i < v->nd; i++) {
    df_variable_t *dv = &df->variables[v->dimids[i]];
    if (dv->type & VT_LONGITUDE) xlim = df->dimensions[v->dimids[i]].size;
    if (dv->type & VT_LATITUDE) ylim = df->dimensions[v->dimids[i]].size;
  }

  mindist = 1e38;
  level = 1;
  ci = -1;
  cj = -1;
  while (ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > xlim && level > ylim)
      return(v->fillvalue);

    while (!finishedLevel) {

      int jfrom = jc - level;
      int jto = jc + level;
      int ifrom = ic - level;
      int ito = ic + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= xlim) || (j >= ylim))
            continue;
	  if (v->nd == 2) {
	    is[0] = j;
	    is[1] = i;
	  } else if (v->nd == 3) {
	    is[0] = kc;
	    is[1] = j;
	    is[2] = i;
	  }
          if ((val = df_get_data_value(df, v, r, is)) != v->missing) {
	    if (val != v->fillvalue) {
	      double  dist = (double)((ic - i) * (ic - i) + (jc - j) * (jc - j));
	      if (dist < mindist) {
		ci = i;
		cj = j;
		mindist = dist;
		nval = val;
	      }
            }
          }
        }
      }
    }
    ++level;
  }
  if (v->nd == 2) {
    is[0] = cj;
    is[1] = ci;
  } else if (v->nd == 3) {
    is[0] = kc;
    is[1] = cj;
    is[2] = ci;
  }
  return(nval);
}

/* Locate the closest cell to is[1], is[0] that contains valid data.
 * This is not very efficient but so what, we don't run this
 * program very often.
 */
static double find_close_bathy(datafile_t *df, df_variable_t *v, int r, int *is)
                 
{
  int i, j, ci, cj;
  double mindist;
  int level;
  int xlim, ylim;
  int ic, jc, kc;
  double val, nval;

  /* Get the indexes */
  if (v->nd == 2) {
    kc = 0;
    jc = is[0];
    ic = is[1];
  } else if (v->nd == 3) {
    kc = is[0];
    jc = is[1];
    ic = is[2];
  } else {
    quit("Landfill only possible for 2D or 3D datafiles.\n");
  }

  /* Get the horizontal grid size */
  for (i = 0; i < v->nd; i++) {
    df_variable_t *dv = &df->variables[v->dimids[i]];
    if (dv->type & VT_LONGITUDE) xlim = df->dimensions[v->dimids[i]].size;
    if (dv->type & VT_LATITUDE) ylim = df->dimensions[v->dimids[i]].size;
  }

  mindist = 1e38;
  level = 1;
  ci = -1;
  cj = -1;
  while (ci < 0) {
    int finishedLevel = 0;
    int edge = 0;

    /* Scanned the whole grid, must be solid every where. */
    if (level > xlim && level > ylim)
      return(v->fillvalue);

    while (!finishedLevel) {

      int jfrom = jc - level;
      int jto = jc + level;
      int ifrom = ic - level;
      int ito = ic + level;

      switch (edge) {
      case 0:                  /* Left edge */
        ito = ifrom;
        break;

      case 1:                  /* Bottom edge */
        ++ifrom;
        jto = jfrom;
        break;

      case 2:                  /* Right edge */
        ifrom = ito;
        ++jfrom;
        break;

      case 3:                  /* Top edge */
        ++ifrom;
        jfrom = jto;
        finishedLevel = 1;
        break;
      }
      ++edge;

      for (j = jfrom; j <= jto; ++j) {
        for (i = ifrom; i <= ito; ++i) {
          if ((i < 0) || (j < 0) || (i >= xlim) || (j >= ylim))
            continue;
	  if (v->nd == 2) {
	    is[0] = j;
	    is[1] = i;
	  } else if (v->nd == 3) {
	    is[0] = kc;
	    is[1] = j;
	    is[2] = i;
	  }
	  val = df_get_data_value(df, v, r, is);
          if (!isnan(val) && val != LANDCELL) {
	    double  dist = (double)((ic - i) * (ic - i) + (jc - j) * (jc - j));
            if (dist < mindist) {
              ci = i;
              cj = j;
              mindist = dist;
	      nval = val;
            }
          }
        }
      }
    }
    ++level;
  }
  if (v->nd == 2) {
    is[0] = cj;
    is[1] = ci;
  } else if (v->nd == 3) {
    is[0] = kc;
    is[1] = cj;
    is[2] = ci;
  }
  return(nval);
}

// local helper functions

/*
 * Finds and returns the distance pointer from the given hash
 */
static double* get_dist_from_hash(hash_table_t** ht, 
				  int *coordids, double coords[], 
				  int ni, df_variable_t *vars)
{
  double *dists = NULL;
  int i,n;

  // First see if the dist hash exists
  if (*ht == NULL)
    *ht = ht_create_d2(ni);

  if ( (dists = (double *)ht_find(*ht, &coords[0])) == NULL ) {
    dists = (double *)malloc(ni*sizeof(double));
    for (i = 0; i < ni; ++i) {
      double dist = 0.0;
      // Note: This is always a 2D distance
      for (n = 0; n < 2; ++n) {
	double dx = VAR_1D(&vars[coordids[n]])[0][i] - coords[n];
	dist += dx * dx;
      }
      dists[i] = sqrt(dist);
    }
    ht_add(*ht, &coords[0], dists);
  }
  
  return(dists);
}


// EOF
