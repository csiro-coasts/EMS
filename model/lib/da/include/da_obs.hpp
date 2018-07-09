/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/include/da_obs.hpp
 *  
 *  Description:
 *  Class definition for DA observation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_obs.hpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#if !defined(_DA_OBS_HPP)
#define _DA_OBS_HPP

#include <gsl/gsl_statistics.h>
#include <string>
#include <map>
#include <stdio.h>
#include "ts_serial.hpp"

class da_obs_exec;

/**
 * \brief Abstract base class for Observation files
 */
class da_obs
{
public:
  da_obs(void);
  ~da_obs(void);
  
  /* Make this a pure-virtual class */
  virtual list<da_obs_exec> & read_obs(double t, da_maps *maps) = 0;

  /* Add variable */
  virtual void add_var(char *var, char *ref_var, int offset, int is3D);

  /* Initialising function */
  void init(char *name, 
	    char *file,
	    char *tunits,
	    double err,
	    double dt);

  /** Getters - mostly for informational purposes
   */
  inline const char *get_name_str(void) {
    return(f_name.c_str());
  }
  inline double get_error(void) {
    return(f_err);
  }
  inline double get_dt(void) {
    return(f_dt);
  }
  virtual inline const char* get_type(void) {
    return("unknown type");
  };

  inline const char *get_fname(void) {
    return(get_df()->name);
  }
  
  virtual const char *get_location_str(void) {
    return("unknown location");
  }

  const char *get_var_str(void);

protected:
  /********/
  /* DATA */
  /********/
  
  /** Name of this observation
   */
  string f_name;
  
  /** Variable id to offset map
   */
  map<int, int> f_vars;
 
  /** Var to state map - useful for debugging
   */
  map<string, string> f_var_state;

  /** Any static error associated with this
   */
  double f_err;
  
  /** Observation dt, in seconds
   */
  double f_dt;
  
  /** Interface to the I/O library
   */
  datafile_t *f_df;
  
  /** Temporary storage
   */
  list <da_obs_exec> f_obs_list;

  /** Location as a string
   */
  string f_loc_str;

  /** Variable names as a string
   */
  string f_var_str;

  /** Whether this variable is 3D
   */
  bool f_is3D;

  /***********/
  /* METHODS */
  /***********/
  virtual void init_file(char *file, char *tunits);

  inline void set_name(char *name) {
    f_name = string(name);
  }

  inline void set_error(double err) {
    f_err = err;
  }

  inline void set_dt(double dt) {
    f_dt = dt;
  }
  virtual datafile_t *get_df(void) {
    return(f_df);
  }
};

/**
 * \brief Subclass for fixed data at a single location. i.e 0D field
 * 
 * This is designed for moored instruments that do not move in
 * space over time
 */
class da_obs_fixed : public da_obs
{
public:
  /** Default constructor
   */
  da_obs_fixed(void) {
    f_s = -1;
  }

  virtual inline const char* get_type(void) {
    return("fixed");
  }
  virtual const char *get_location_str(void);

  /* This function is mandatory */
  virtual list<da_obs_exec> & read_obs(double t, da_maps *maps);

  /* fixed accepts only an actual position */
  void set_location(da_maps *maps, double lon, double lat, double dep);

protected:
  /* The sparse locations */
  int f_s;
};


/**
 * \brief Subclass for multiple fixed data - low cost nodes,
 * thermister string etc.
 *
 */
class da_obs_multi_fixed : public da_obs_fixed
{
public:
  virtual inline const char* get_type(void) {
    return("multi-fixed");
  }
};

/**
 * \brief Subclass for fixed data of 2D horizontal field
 *
 * This is suitable for 2D variables whose fields do not vary
 * in space over time. Presently this is limited to the surface
 * (i.e. depth=0) therefore appropriate for satellite data. It can
 * easily be modified to handle arbitrary (but constant) depth.
 *
 * No time averaging of data
 */
class da_obs_2d_field : public da_obs_fixed
{
public:
  /**
   * Destructor
   */
  ~da_obs_2d_field(void) {
    delete tS;
  }
  
  virtual inline const char* get_type(void) {
    return("2d-field");
  }

  virtual const char *get_location_str(void) {
    return("Locations interpolated onto model grid");
  }

  virtual list<da_obs_exec> & read_obs(double t, da_maps *maps);
  void init_file(char *file, char *tunits);
  datafile_t *get_df(void) {
    return(tS->get_ts()->df);
  }

  void add_var(char *var, char *ref_var, int offset, int is3D);

private:
  /* Serialised timeseries object */
  tsSerial *tS;
};


/**
 * \brief Subclass for fixed data of full 3D field
 *
 * Not yet implemented.
 */
class da_obs_3d_field : public da_obs_2d_field
{
public:
  virtual inline const char* get_type(void) {
    return("3d-field");
  }
};


/**
 * \brief Subclass for mobile data platform
 *
 * This observation may vary in space over time. The location is read
 * from file at each DA time step
 */
class da_obs_mobile : public da_obs
{
public:
  virtual list<da_obs_exec> & read_obs(double t, da_maps *maps);
  virtual const char *get_location_str(void);
  
  virtual inline const char* get_type(void) {
    return("mobile");
  }

  virtual void set_location(char *lon, char *lat, char *dep);

protected:
  /**
   * Datafile id's of lon, lat and depth, (x,y,z), respectively
   */
  int f_lon_id;
  int f_lat_id;
  int f_dep_id;
};


/**
 * \brief Subclass for glider data
 *
 * The 2D glider location is first averaged over its horizontal extent
 * over the time period and then binned into each of the layers in the model
 */
class da_obs_glider : public da_obs_mobile
{
public:
  virtual list<da_obs_exec> & read_obs(double t, da_maps *maps);
  virtual inline const char* get_type(void) {
    return("glider");
  }
};


/*********************/
/* EXECUTION CONTEXT */
/*********************/
/**
 * \brief Simple class to hold the execution context of each
 * observation
 */
class da_obs_exec
{
  /**
   * Rely on the default constructors from the compiler
   */
public:
  /** The mean value */
  double f_val;

  /** The standard deviation */
  double f_err;

  /** The value of the sparse coordinate + the variable offset */
  int f_gs;

  /** The longitude */
  double f_x;

  /** The latitude */
  double f_y;
};


#endif

/* end of da_obs.hpp */
