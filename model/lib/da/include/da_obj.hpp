/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/include/da_obj.hpp
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
 *  $Id: da_obj.hpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#if !defined(_DA_OBJ_HPP)
#define _DA_OBJ_HPP

extern "C" {
  #include "da_interface.h"
}
#include "da_obs.hpp"
#include <map>

/**
 * \brief This is the main Data Assimilation object
 *
 * This serves as the main gateway for the C interface as well
 */
class da_obj
{
 public:
  da_obj(double ls);
  ~da_obj(void);
  int allocA(int nrows, int ncols);
  void fillA(int row, int col, double val);
  void fillA(int col, double *vals);
  int checkA(void);
  void writeA(const char *fname);
  void dumpMatrix(const char *fname, const gsl_matrix *M);
  void add_obs(da_obs *obs);

  /* Set/Get maps */
  inline void set_maps(da_maps *maps) {
    f_maps = maps;
  }
  inline da_maps *get_maps(void) {
    return(f_maps);
  }

  /* Fill in the maps */
  inline void fill_map(int ls, int i, int j, int k) {
    if (k<0)
      f_maps->fill_map(ls, i, j);
    else
      f_maps->fill_map(ls, i, j, k);
  }    

  inline int get_num_obs(void) {
    return(f_obs_arr.size());
  }

  /* Gateway's to specific obs
   */
  inline void *get_obs_ptr(int i) {
    return((void *)f_obs_arr[i]);
  }

  /* Read all available data */
  int read_all_obs(double t);

  /* Adds the state name and offset */
  void add_state(const char *str, int off);

  /* Get the offset given variable name */
  int get_state_offset(char *ref_var);

  /* Sets the background state vector */
  void set_wb(int index, double data);

  /* Fills in the location vectors */
  void set_xy(int index, double x, double y);

  /* Get a single value of the analysis field */
  double get_wa(int index);

  /* Get a particular observation */
  double get_wo_val(int index);
  
  /* Get a particular observation err */
  double get_wo_err(int index);

  /* Main work horse function */
  void do_analysis(void);

  /* Calculate and apply the localisation factor */
  void applyLocalisaton(gsl_matrix *AHA);

  /* For when there are no observations */
  void copy_wb_wa(void);

private:
  /**
   * The localisation spread, in degrees
   */
  double f_ls;

  /**
   * Holds the anomaly field
   */
  gsl_matrix *f_A;

  /**
   * Error covariance matrix
   */
  gsl_matrix *f_R;

  /**
   * Matrix made up of the rows of A that have observations
   */
  gsl_matrix *f_HA;

  /**
   * Kalman gain matrix
   */
  gsl_matrix *f_K;

  /**
   * The sparse observations vector
   */
  gsl_vector *f_wo;
  
  /**
   * The full NaN padded observations vector
   */
  gsl_vector *f_full_wo_val;
  
  /**
   * The full NaN padded observations vector of standard deviations
   */
  gsl_vector *f_full_wo_err;
  
  /**
   * The full state vector
   */
  gsl_vector *f_wb;

  /**
   * The background for where we have observations
   */
  gsl_vector *f_Hwb;

  /**
   * The final analysis vector
   */
  gsl_vector *f_wa;

  /**
   * holds longitudes
   */
  gsl_vector *f_x;

  /**
   * holds latitudes
   */
  gsl_vector *f_y;

  /**
   * Pointer to maps
   */
  da_maps *f_maps;

  /**
   * Map of offsets with state name as the key
   */
  map<string, int> f_states;

  /** array of observations
   */
  vector <da_obs *> f_obs_arr;

  /** list of actual observations
   */
  list <da_obs_exec> f_obs_exec;

  /*
   * METHODS
   */
  void constructK(void);
  void constructHA_Hwb(void);
  void constructR(void);
  void constructWo(void);
  void construct_full_wo_val(void);
  void construct_full_wo_err(void);
};

#endif

/* end of da_obj.hpp */

