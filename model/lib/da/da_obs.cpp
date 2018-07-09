/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/da_obs.cpp
 *  
 *  Description:
 *  Data assimilation Observation classes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_obs.cpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#include "da_utils.hpp"
#include "da_obs.hpp"
#include <map>
#include <cstdlib>
#include <algorithm>
#include <set>

/**************/
/* BASE CLASS */
/**************/

/**
 * Constructor
 */
da_obs::da_obs(void)
{
  f_df = df_alloc();
}


/**
 * Destructor
 */
da_obs::~da_obs(void)
{
  /* This frees the pointer itself as well */
  df_free(f_df);
}


/** Init method
 *
 */
void da_obs::init(char *name, 
		  char *file,
		  char *tunits,
		  double err,
		  double dt)
{
  /* Initialise file I/0 */
  init_file(file, tunits);

  /* Set various properties */
  set_name(name);
  set_error(err);
  set_dt(dt);
}


/** Initialises the file
 */
void da_obs::init_file(char *file, char *tunits)
{
  int index = -1;

  /* Initialise the I/O struct */
  df_read(file, f_df);

  /* Set time as the record coordinate */
  index = df_get_index(f_df, "t");
  if (index >= 0)
    df_set_record(f_df, index);
  else if ((index = df_get_index(f_df, "time")) >= 0)
    df_set_record(f_df, index);
  else if ((index = df_get_index(f_df, "Time")) >= 0)
    df_set_record(f_df, index);
  else {
    /* error out */
    quit("No time variable found in '%s'\n", file);
  }

  /* Ensure time units consistency */
  tm_change_time_units(f_df->rec_units, tunits, f_df->records, f_df->nrecords);
}


/** Add a variable to read from file
 */
void da_obs::add_var(char *var, char *ref_var, int offset, int is3D)
{
  int id = df_get_index(get_df(), var);
  if (id < 0) {
    /* silently ignore */
    return;
  }
  /* All good */
  f_vars[id] = offset;

  /* Also cache the variable to state mapping */
  f_var_state[string(var)] = string(ref_var);

  /*
   * Cache whether this is a 3D variable or not - needed for when the
   * variable is 3D but the observations are 2D, e.g. remotely sensed
   * SST
   */
  f_is3D = (is3D != 0) ? true : false;
}


/** Variable string
 */
const char *da_obs::get_var_str(void)
{
  if (f_var_str.empty()) {
    map<string,string>::iterator pos;
    for (pos = f_var_state.begin(); pos != f_var_state.end(); pos++) {
      /* Construct string of the form (salt=salt2) etc ... */
      string var  = pos->first;
      string rvar = pos->second;

      f_var_str += "(" + rvar;
      if (var != rvar)
	f_var_str += "=" + var;
      f_var_str += ")";
    }
  }
  
  return(f_var_str.c_str());
}

/***************/
/* FIXED CLASS */
/***************/

/**
 * set location of this fixed
 * Error out if location not found in domain
 */
void da_obs_fixed::set_location(da_maps *maps, double lon, double lat, double dep)
{
  char buf[1024];
  int i,j,k;
  f_s = maps->getS(lon, lat, dep);
  if (f_s < 0)
    quit("Unable to find location (lon,lat)=(%f,%f) for '%s'\n", lon, lat, 
	                                                          get_name_str());
  /* Set location string */
  maps->getIJK(f_s, &i, &j, &k);
  sprintf(buf, "(%d, %d, %d)", i, j, k);
  f_loc_str = "Fixed 3D coord = " + string(buf);
}


/** Location string
 */
const char *da_obs_fixed::get_location_str(void)
{
  return(f_loc_str.c_str());
}

/** Read method
 *
 */
list<da_obs_exec> &da_obs_fixed::read_obs(double t, da_maps *maps)
{
  map<int,int>::iterator pos;

  /* Start with a clean slate */
  f_obs_list.clear();

  /* Loop over each variable */
  for (pos = f_vars.begin(); pos != f_vars.end(); pos++) {
    int     vid   = pos->first; 
    int     offst = pos->second;
    double *data  = NULL;
    double  m_val = 0;
    double  s_val = 0;
    int     nrecs = 0;
    int i,j,k;

    /* Read the records */
    nrecs = df_eval_period(f_df, vid, t, f_dt, &data);
    if (nrecs > 1) {
      da_obs_exec obs;

      /* 
       * Calculate the mean
       * Note for noisy data, it might be better to use median instead
       */
      m_val = gsl_stats_mean(data, 1, nrecs);

      /* Calculate the standard deviation */
      s_val = gsl_stats_sd_m(data, 1, nrecs, m_val);
      
      /* Set values */
      obs.f_val = m_val;
      obs.f_err = s_val + f_err;
      obs.f_gs  = f_s + offst;
      /* Set lat/lon */
      maps->getIJK(f_s, &i, &j, &k);
      maps->getXY(i, j, &obs.f_x, &obs.f_y);

      /* Add to list */
      f_obs_list.push_back(obs);
    }
  }

  /* Return a reference to the list */
  return(f_obs_list);
}


/******************/
/* 2D-FIELD CLASS */
/*****************/
/**
 * This init method uses ts_read
 */
void da_obs_2d_field::init_file(char *file, char *tunits)
{
  tS = new tsSerial(file);
  
  /* Make timeunits consistent */
  ts_convert_time_units(tS->get_ts(), tunits);

  /* Seed the random number generator */
  srand(time(NULL));
}

/** Add a variable to read from file for 2d-field
 */
void da_obs_2d_field::add_var(char *var, char *ref_var, int offset, int is3D)
{
  /* Call the base class method */
  da_obs::add_var(var, ref_var, offset, is3D);

  /* Perform serialisation of the timeseries data */
  tS->serialise(var);
}

/**
 * o) create a list of all valid points within the domain
 * o) randomly shuffle and pick the first MAX_POINTS worth
 */
list<da_obs_exec> &da_obs_2d_field::read_obs(double t, da_maps *maps)
{
  map<int,int>::iterator pos;
  int i, r0, r1;
  double frac;
  datafile_t *df = tS->get_ts()->df;
  int MAX_POINTS = 2000; /* Or make it some fraction? */
  double mtime, otime;

  /* Start with a clean slate */
  f_obs_list.clear();

  /* We apply daily values, regardless of the actual time in the file */
  df_find_record(df, t, &r0, &r1, &frac);

  /* Check that we're on the same day: not sure if this is the best way? */
  mtime = t/86400.0;
  otime = df->records[r1]/86400.0; // xxx r1 is wrong general - maybe snap to nearest?

  if ( ((int)(mtime)) == ((int)(otime)) ) {
    /* Loop over each variable */
    for (pos = f_vars.begin(); pos != f_vars.end(); pos++) {
      int vid   = pos->first; 
      int offst = pos->second;
      /* sparse array */
      vector<int> wcells;
      /* Keep track of duplicates */
      set<int,greater<int> > locs;

      /* Read data for this variable */
      df_read_records(df, df_get_variable(df, vid), r1, 1);

      warn("DA: Processing %s at obs time = %.4f days for model time = %.4f days\n",
	   get_name_str(), otime, mtime);
      /* 
       * Initialise the sparse array for this file
       * These are all the valid points in the data file that fall
       * within the model domain
       */
      for (i=0; i<tS->getSize(); i++) {
	double X = tS->getX(i);
	double Y = tS->getY(i);
	int coords[] = {tS->getC0(i), tS->getC1(i)};

	/* Get the actual data value for this coord */
	double val = df_get_data_value(df, df_get_variable(df, vid), r1, coords);
	
	/* What about missing values and other flags? */
	if (!isnan(val) && maps->inDomain(X, Y))
	  wcells.push_back(i);
      }

      /* 
       * Now we take a random subset of these points
       */
      random_shuffle(wcells.begin(), wcells.end());

      for (i=0; i<wcells.size(); i++) {
	da_obs_exec obs;
	int    s = wcells[i];
	double X = tS->getX(s);
	double Y = tS->getY(s);
	int coords[] = {tS->getC0(s), tS->getC1(s)};
	double val = df_get_data_value(df, df_get_variable(df, vid), r1, coords);
	int c = maps->getS(X, Y);
	
	/* Insert into the set of empty model grid cells */
	pair<set<int,greater<int> >::iterator,bool> status = locs.insert(c);

	if (status.second) {
	  /* Set values */
	  obs.f_val = val;
	  obs.f_err = f_err;
	  obs.f_x   = X;
	  obs.f_y   = Y;

	  // Really need to sort out the 2D vs 3D sparse coord
	  if (!f_is3D)
	    quit("DA:2d_field read error - not 3D\n");
	  obs.f_gs  = c + offst;

	  /* Add to list */
	  f_obs_list.push_back(obs);
	  
	  if (f_obs_list.size() == MAX_POINTS)
	    break;
	}
      }
      warn("DA: Using %d out of %d SST obs within the domain\n",
	   f_obs_list.size(), wcells.size());
    }
  } else
    warn("DA: No obs in %s found for model time = %.4f days\n",get_name_str(), mtime);
  
  /* Return a reference to the list */
  return(f_obs_list);
}

/******************/
/* 3D-FIELD CLASS */
/******************/

/* not yet implemented */

/****************/
/* MOBILE CLASS */
/****************/
/** Read method for mobile observations
 *
 * Here we simply average out the 3D location
 */
list<da_obs_exec> &da_obs_mobile::read_obs(double t, da_maps *maps)
{
  int   nrecs, nrecs0, nrecs1, nrecs2;
  double *lat = NULL;
  double *lon = NULL;
  double *dep = NULL;
  double x,y,z;
  int loc_i, loc_j;
  int s;
  map<int,int>::iterator pos;

  f_obs_list.clear();

  /* Read the variables from file */
  nrecs0 = df_eval_period(f_df, f_lat_id, t, f_dt, &lat);
  nrecs1 = df_eval_period(f_df, f_lon_id, t, f_dt, &lon);
  nrecs2 = df_eval_period(f_df, f_dep_id, t, f_dt, &dep);

  /* Fix up the sign of the depth variable */
  df_fixup_variable_sign(f_df, f_dep_id);

  /* Sanity check */
  if ( (nrecs0 != nrecs1) && (nrecs1 != nrecs2) )
    quit("Mobile position data are not all of equal length for observation '%s'\n",
	                                                   get_name_str());
  nrecs = nrecs0;

  if (nrecs > 0) {
    /* Get the 2D position */
    x = gsl_stats_mean(lon, 1, nrecs);
    y = gsl_stats_mean(lat, 1, nrecs);
    z = gsl_stats_mean(dep, 1, nrecs);
    
    /* Get the sparse coordinate from the map */
    s = maps->getS(x, y, z);
    
    /* See if we're in the domain */
    if (s >= 0) {
      /* Loop over each variable */
      for (pos = f_vars.begin(); pos != f_vars.end(); pos++) {
	int     vid   = pos->first; 
	int     offst = pos->second;
	double *vdata;
	int nvals;
	da_obs_exec obs;
	
	if ((nvals=df_eval_period(f_df, vid, t, f_dt, &vdata)) != nrecs)
	  quit("Variable '%s' is of length %d but should be %d\n",
	       df_get_variable_name_by_id(f_df, vid), nvals, nrecs);
	
	/* Compute statistics and observation */
	obs.f_val = gsl_stats_mean(vdata, 1, nrecs);
	obs.f_err = gsl_stats_sd_m(vdata, 1, nrecs, obs.f_val) + f_err;
	obs.f_gs  = s + offst;
	obs.f_x   = x;
	obs.f_y   = y;
	
	/* Add this obs to the list */
	f_obs_list.push_back(obs);
      }
    }
  }
  
  /* Return a reference to the list */
  return(f_obs_list);
}
  

/** Sets up the location
 */
void da_obs_mobile::set_location(char *lon, char *lat, char *dep)
{
  /* Find the longitude varid */
  if ( (f_lon_id = df_get_index(f_df, lon)) < 0 )
    quit("Unable to locate longitude variable '%s' in file '%s' for observation '%s'\n", 
	 lon, f_df->name, get_name_str());
  
  /* Find the latitude varid */
  if ( (f_lat_id = df_get_index(f_df, lat)) < 0 )
    quit("Unable to locate latitude variable '%s' in file '%s' for observation '%s'\n", 
	 lat, f_df->name, get_name_str());
  
  /* Find the depth varid */
  if ( (f_dep_id = df_get_index(f_df, dep)) < 0 )
    quit("Unable to locate depth variable '%s' in file '%s' for observation '%s'\n", 
	 dep, f_df->name, get_name_str());
}


/** Location string
 */
const char *da_obs_mobile::get_location_str(void)
{
  if (f_loc_str.empty()) {
    string xstr = string(df_get_variable_name_by_id(f_df, f_lon_id));
    string ystr = string(df_get_variable_name_by_id(f_df, f_lat_id));
    string zstr = string(df_get_variable_name_by_id(f_df, f_dep_id));

    string xs = string("X");
    string ys = string("Y");
    string zs = string("Z");

    f_loc_str  = string("Reading location info from file ");
    f_loc_str += "(" + xs;
    if (xs != xstr)
      f_loc_str += "=" + xstr;
    f_loc_str += ")";
    f_loc_str += "(" + ys;
    if (ys != ystr)
      f_loc_str += "=" + ystr;
    f_loc_str += ")";
    f_loc_str += "(" + zs;
    if (zs != zstr)
      f_loc_str += "=" + zstr;
    f_loc_str += ")";
  }
  
  return(f_loc_str.c_str());
}

/****************/
/* GLIDER CLASS */
/****************/
/** Read method for glider observations
 * 
 * Here we average out the 2D location but bin in the vertical
 *
 */
list<da_obs_exec> &da_obs_glider::read_obs(double t, da_maps *maps)
{
  int   nrecs, nrecs0, nrecs1, nrecs2;
  double *lat   = NULL;
  double *lon   = NULL;
  double *depth = NULL;
  double x,y;
  int loc_i, loc_j;

  f_obs_list.clear();

  /* Read the variables from file */
  nrecs0 = df_eval_period(f_df, f_lat_id, t, f_dt, &lat);
  nrecs1 = df_eval_period(f_df, f_lon_id, t, f_dt, &lon);
  nrecs2 = df_eval_period(f_df, f_dep_id, t, f_dt, &depth);

  /* Fix up the sign of the depth variable */
  df_fixup_variable_sign(f_df, f_dep_id);

  /* Sanity check */
  if ( (nrecs0 != nrecs1) && (nrecs1 != nrecs2) )
    quit("Glider position data are not all of equal length for observation '%s'\n",
	                                                   get_name_str());
  nrecs = nrecs0;

  /* Get the 2D position */
  x = gsl_stats_mean(lon, 1, nrecs);
  y = gsl_stats_mean(lat, 1, nrecs);

  /* See if we're in the domain */
  if (maps->getIJ(x, y, &loc_i, &loc_j)) {
    int n;
    map<int,int>::iterator p_var;
    multimap<int,int> deps;

    /* Bin the indicies with the same sparse coord */
    for (n=0; n<nrecs; n++) {
      int S = maps->getS(loc_i, loc_j, depth[n]);
      if (S >= 0) /* Only work with valid sparse coords */
	deps.insert(make_pair(S, n));
    }

    /* For each variable */
    for (p_var = f_vars.begin(); p_var != f_vars.end(); p_var++) {
      int vid   = p_var->first; 
      int offst = p_var->second;
      int nvals;
      double *vdata;
      multimap<int,int>::iterator p_bin;
      /* 
       * Read in all values and make sure the length is the same as
       * above
       */
      if ((nvals=df_eval_period(f_df, vid, t, f_dt, &vdata)) != nrecs)
	quit("Variable '%s' is of length %d but should be %d\n",
	     df_get_variable_name_by_id(f_df, vid), nvals, nrecs);
      
      /* For each bin */
      for (p_bin = deps.begin(); p_bin != deps.end();) {
	int key      = p_bin->first;   /* The sparse coord associated with this bin */
	int num_vals = deps.count(key);/* The number of values in this bin */

	if (num_vals > 1) {
	  da_obs_exec obs;
	  vector<double> vals;
	  int i;
	  /* For each value in this bin */
	  for (i=0; i<num_vals; i++,p_bin++) {
	    double val = vdata[p_bin->second];
	    /* Skip over any nan's */
	    if (isfinite(val))
	      vals.push_back(val);
	  }
	  /* Compute statistics and observation */
	  obs.f_val = gsl_stats_mean(&vals[0], 1, vals.size());
	  obs.f_err = gsl_stats_sd_m(&vals[0], 1, vals.size(), obs.f_val) + f_err;
	  obs.f_gs  = key + offst;
	  obs.f_x   = x;
	  obs.f_y   = y;

	  /* Add this obs to the list */
	  f_obs_list.push_back(obs);
	} else {
	  /* Move on */
	  p_bin++;
	}
      }
    }
  }

  /* Return a reference to the list */
  return(f_obs_list);
}

/* end da_obs.cpp */
