/** \file da_interface.cpp
 *  \brief This is the C interface to the DA library
 *
 *  The DA library is written in C++ which can be used directly. This file
 *  provides a handy C interface
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_interface.cpp 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#include <list>
#include <vector>
#include <string>
using namespace std;

#include <gsl/gsl_matrix.h>
extern "C" {
  #include "errfn.h"
  double **d_alloc_2d(int n1, int n2);
  void d_free_2d(double **pp);
}
#include "da_utils.hpp"
#include "da_obj.hpp"

/**
 * Use this function when invoking the library on its own
 */
extern "C"
void da_warn_on_quit(void)
{
  errfn_set_quit(errfn_warn_default);
}

/** Constructor
 * This object serves as the main gateway for the DA library
 *
 * @return void pointer to the newly allocated da_obj object
 */
extern "C"
void *da_create_object(double da_ls)
{
  da_obj *obj = new da_obj(da_ls);
  
  if (obj == NULL)
    quit("DA: Memory allocation error\n");

  return( (void *) obj);
}


/** Destructor
 * Deletes memory associated with the da_obj object
 */
extern "C"
void da_destroy_object(void *obj)
{
  if (obj != NULL)
    delete( ((da_obj*)obj) );
}


/** Allocate memory for A
 * @param daobj da_obj pointer
 * @param nrows number of states x number of cells
 * @param ncols number of members
 * @return non-zero exit status on error
 */
extern "C"
int da_allocA(void *daobj, int nrows, int ncols)
{
  da_obj *obj = (da_obj *)daobj;
  return(obj->allocA(nrows, ncols));
}


/** Fills one element of A
 * @param daobj da_obj pointer
 * @param row row number
 * @param col column number
 * @param val data value
 *
 */
extern "C"
void da_fill_A(void *daobj, int row, int col, double val)
{
  da_obj *obj = (da_obj *)daobj;
  obj->fillA(row, col, val);
}

 
/** Fills one column of A
 * @param daobj da_obj pointer
 * @param col column number
 * @param vals pointer to first element of column
 *
 */
extern "C"
void da_fill_A_col(void *daobj, int col, double *vals)
{
  da_obj *obj = (da_obj *)daobj;
  obj->fillA(col, vals);
}

 
/** Checks A to make sure there is valid data throughout - useful for
 * debugging
 * @param daobj da_obj pointer
 * @return 1 on success, 0 error
 *
 */
extern "C"
int da_check_A(void *daobj)
{
  return( ((da_obj*)daobj)->checkA() );
}


/** Writes A to file - useful for debugging
 * @param daobj da_obj pointer
 * @param fname name of file to write to
 *
 */
extern "C"
void da_writeA(void *daobj, const char *fname)
{
  ((da_obj*)daobj)->writeA(fname);
}


/** Creates the maps needed to go from cartesian to sparse
 * @param daobj da_obj pointer
 * @param gx (nx+1) x (ny+1) number of x coordinates
 * @param gy (nx+1) x (ny+1) number of y coordinates
 * @param gz (nz+1) vertical layer coordinates
 * @param nx number of x cell centers
 * @param ny number of y cell centers
 * @param nz number of z layers
 */
extern "C"
void da_create_maps(void *daobj, double **gx, double **gy, double *gz,
		                                             int nx, int ny, int nz)
{
  da_maps *maps = new da_maps(gx, gy, gz, nx, ny, nz);
  ((da_obj*)daobj)->set_maps(maps);
}


/** Creates the maps needed to go from cartesian to sparse
 * Special case for Matlab. I don't yet know how to pass 2D pointers
 * from Matlab so here we allocate arrays specifically.
 *
 * @param daobj da_obj pointer
 * @param igx 1D (nx+1) x (ny+1) number of x coordinates - assume col maj
 * @param igy 1D (nx+1) x (ny+1) number of y coordinates - assume col maj
 * @param gz (nz+1) vertical layer coordinates
 * @param nx number of x cell centers
 * @param ny number of y cell centers
 * @param nz number of z layers
 */
extern "C"
void da_create_maps_ml(void *daobj, double *igx, double *igy, double *gz,
		       int nx, int ny, int nz)
{
  da_maps *maps;
  static double **gx = NULL;
  static double **gy = NULL;
  int j, i;

  if (gx != NULL)
    d_free_2d(gx);
  if (gy != NULL)
    d_free_2d(gy);
  
  gx = d_alloc_2d(nx+1, ny+1);
  gy = d_alloc_2d(nx+1, ny+1);

  for (i=0; i<=nx; i++)
    for (j=0; j<=ny; j++) {
      gx[j][i] = igx[i*(ny+1)+j];
      gy[j][i] = igy[i*(ny+1)+j];
    }

  maps = new da_maps(gx, gy, gz, nx, ny, nz);
  ((da_obj*)daobj)->set_maps(maps);
}


/** Fill in the cart to sparse maps. This map is filled in so that
 * observation values can find the appropriate location in the state vectors
 * @param daobj da_obj pointer
 * @param ls local sparse index (i.e. counter value)
 * @param i i index value
 * @param j j index value
 * @param k k index value (specify -1 for 2D)
 */
extern "C"
void da_fill_map(void *daobj, int ls, int i, int j, int k)
{
  ((da_obj*)daobj)->fill_map(ls, i, j, k);
}

/** Fill in the cart to sparse maps. This map is filled in so that
 * observation values can find the appropriate location in the state
 * vectors
 * Vectorised form
 * @param daobj da_obj pointer
 * @param ns number of values
 * @param i i index values
 * @param j j index values
 * @param k k index values (specify -1 for 2D)
 */
extern "C"
void da_fill_map_ml(void *daobj, int ns, int *i, int *j, int *k)
{
  int ls;
  for (ls=0; ls<ns; ls++)
    ((da_obj*)daobj)->fill_map(ls, i[ls], j[ls], k[ls]);
}

/** Add anomaly state name
 * @param daobj da_obj pointer
 * @param str name of state
 * @param off row offset
 */
extern "C"
void da_add_state(void *daobj, const char *str, int off)
{
  ((da_obj*)daobj)->add_state(str, off);
}


/** Creates fixed observation object
 * @param name observation name
 * @param file name of file to read observations from
 * @param tunits units to time to convert to. i.e. the time units with
 * which read will be called
 * @param err instrument error
 * @param dt time interval over which to calculate statistics
 * @return the fixed observation object (da_obs_fixed) as void *
 *
 */
extern "C"
void *da_create_fixed_obs(char *name, 
			    char *file, 
			    char *tunits,
			    double err,
			    double dt)
{
  da_obs_fixed *obs = new da_obs_fixed();

  obs->init(name, file, tunits, err, dt);

  return( (void *)obs );
}

/** Creates glider observation object
 * @param name observation name
 * @param file name of file to read observations from
 * @param tunits units to time to convert to. i.e. the time units with
 * which read will be called
 * @param err instrument error
 * @param dt time interval over which to calculate statistics
 * @return the glider observation object (da_obs_glider) as void *
 *
 */
extern "C"
void *da_create_glider_obs(char *name, 
			   char *file, 
			   char *tunits,
			   double err,
			   double dt)
{
  da_obs_glider *obs = new da_obs_glider();
  
  obs->init(name, file, tunits, err, dt);
  
  return( (void *)obs );
}


/** Creates 2D field observation object
 * @param name observation name
 * @param file name of file to read observations from - usually a multi-netcdf file
 * @param tunits units to time to convert to. i.e. the time units with
 * which read will be called
 * @param err instrument error
 * @return the 2d-field observation object (da_obs_2d_field) as void *
 *
 */
extern "C"
void *da_create_2d_field_obs(char *name, 
			     char *file, 
			     char *tunits,
			     double err)
{
  da_obs_2d_field *obs = new da_obs_2d_field();
  
  obs->init(name, file, tunits, err, -1);
  
  return( (void *)obs );
}


/** Set fixed location
 * @param obs observation object
 * @param obj da object
 * @param lon longitude value
 * @param lat latitude value
 * @param dep depth value
 */
extern "C"
void da_fixed_obs_set_location(void *obs, void *obj, double lon, double lat, double dep)
{
  ((da_obs_fixed *)obs)->set_location(((da_obj *)obj)->get_maps(), lon, lat, dep);
}


/** Sets glider location info
 * @param obs observation object
 * @param lon name of variable associated with longitude
 * @param lat name of variable associated with latitude
 * @param dep name of variable associated with depth
 */
extern "C"
void da_glider_obs_set_location(void *obs, char *lon, char *lat, char *dep)
{
  ((da_obs_glider *)obs)->set_location(lon, lat, dep);
}


/** Add a variable to observation
 * @param daobj da_obj pointer
 * @param daobs observation object
 * @param var of variable
 * @param ref_var of variable this refers to in the state vector
 */
extern "C"
void da_obs_add_var(void *daobj, void *daobs, char *var, char *ref_var,
		    int is3D)
{
  int offset = ((da_obj *)daobj)->get_state_offset(ref_var);
  ((da_obs *)daobs)->add_var(var, ref_var, offset, is3D);
}


/** Adds the observation object to the main da object
 * @param daobj da_obj pointer
 * @param obs observation object
 *
 */
extern "C"
void da_add_obs(void *daobj, void *obs)
{
  ((da_obj*)daobj)->add_obs((da_obs *)obs);
}


/** Geteway to the main work horse function
 * @param daobj da_obj pointer
 *
 */
extern "C"
void da_do_analysis(void *daobj)
{
  ((da_obj*)daobj)->do_analysis();
}


/** Gateway to set the current background state vector
 * @param daobj da_obj pointer
 * @param index index value (including offset) in the big state vector
 * @param data the background value
 */
extern "C"
void da_set_wb(void *daobj, int index, double data)
{
  ((da_obj*)daobj)->set_wb(index, data);
}

/** Gateway to set the locations
 * @param daobj da_obj pointer
 * @param index index value (including offset) in the big state vector
 * @param x longitude value
 * @param y latitude value
 */
extern "C"
void da_set_xy(void *daobj, int index, double x, double y)
{
  ((da_obj*)daobj)->set_xy(index, x, y);
}

/** Gateway to set the current background state vector as a column
 * @param daobj da_obj pointer
 * @param ndata number of data values
 * @param data the background value
 */
extern "C"
void da_set_wb_col(void *daobj, int ndata, double *data)
{
  int index;
  for (index=0; index<ndata; index++)
    ((da_obj*)daobj)->set_wb(index, data[index]);
}

/** Gateway to retrieve the analysis vector
 * @param daobj da_obj pointer
 * @param index index value
 * @return the analysis value
 */
extern "C"
double da_get_wa(void *daobj, int index)
{
  double val = ((da_obj*)daobj)->get_wa(index);
  return(val);
}

/** Gateway to fill in the analysis vector
 * @param daobj da_obj pointer
 * @param ndata number or elements in data
 * @param data allocated data array
 */
extern "C"
void da_fill_wa(void *daobj, int ndata, double *data)
{
  int index;
  for (index=0; index<ndata; index++)
    data[index] = ((da_obj*)daobj)->get_wa(index);
}


/** Gateway to fill in the observation values
 * @param daobj da_obj pointer
 * @param ndata number or elements in data
 * @param data allocated data array
 */
extern "C"
void da_fill_wo_val(void *daobj, int ndata, double *data)
{
  int index;
  for (index=0; index<ndata; index++)
    data[index] = ((da_obj*)daobj)->get_wo_val(index);
}

/** Gateway to fill in the observation errors
 * @param daobj da_obj pointer
 * @param ndata number or elements in data
 * @param data allocated data array
 */
extern "C"
void da_fill_wo_err(void *daobj, int ndata, double *data)
{
  int index;
  for (index=0; index<ndata; index++)
    data[index] = ((da_obj*)daobj)->get_wo_err(index);
}


/** Gateway to retrieve the observations vector
 * @param daobj da_obj pointer 
 * @param index index value
 * @return the observation value
*/
extern "C"
double da_get_wo_val(void *daobj, int index)
{
  double val = ((da_obj*)daobj)->get_wo_val(index);
  return(val);
}


/** Gateway to retrieve the observations errors vector
 * @param daobj da_obj pointer
 * @param index index value
 * @return the observation error value
 */
extern "C"
double da_get_wo_err(void *daobj, int index)
{
  double err = ((da_obj*)daobj)->get_wo_err(index);
  return(err);
}


/**
 * Directive to read in the observations from file for the given time
 * @param daobj da_obj pointer
 * @param t time value to read. Must be in the units specified at
 * observation object creation time
 * @return the number of observations found
 */
extern "C"
int da_read_all_obs(void *daobj, double t)
{
  return(((da_obj*)daobj)->read_all_obs(t));
}


/** Gets the total number of observations specified
 * @param daobj da_obj pointer
 * @return the number of observations
 */
extern "C"
int da_get_num_obs(void *daobj)
{
  return(((da_obj*)daobj)->get_num_obs());
}


/** Gets a particular observation as a void pointer
 * @param daobj da_obj pointer
 * @param i index of observation
 * @return observation object as void pointer
 */
extern "C"
void *da_get_obs_ptr(void *daobj, int i)
{
  return(((da_obj*)daobj)->get_obs_ptr(i));
}


/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
const char *da_obs_get_name(void *daobs) {
  return(((da_obs *)daobs)->get_name_str());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
double da_obs_get_error(void *daobs) {
  return(((da_obs *)daobs)->get_error());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
double da_obs_get_dt(void *daobs) {
  return(((da_obs *)daobs)->get_dt());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
const char* da_obs_get_type(void *daobs) {
  return(((da_obs *)daobs)->get_type());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
const char *da_obs_get_fname(void* daobs) {
  return(((da_obs *)daobs)->get_fname());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
const char *da_obs_get_vars(void *daobs) {
  return(((da_obs *)daobs)->get_var_str());
}

/** Get name of observation
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
const char *da_obs_get_location(void *daobs) {
  return(((da_obs *)daobs)->get_location_str());
}

/** Copies the background field into the analysis
 * @param daobs daobs pointer (from a da_get_obs_ptr call)
 */
extern "C"
void da_copy_wb_wa(void *daobj) {
  ((da_obj *)daobj)->copy_wb_wa();
}
/* end of da_interface.cpp */
