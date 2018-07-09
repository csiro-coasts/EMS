/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/da/include/da_interface.h
 *  
 *  Description:
 *  Data assimilation library C interface header file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: da_interface.h 5850 2018-06-29 05:18:59Z riz008 $
 *
 */

#if !defined(_DA_INTERFACE_H)
#define _DA_INTERFACE_H

void da_warn_on_quit(void);
void *da_create_object(double da_ls);
void da_destroy_object(void *obj);
int da_allocA(void *daobj, int nrows, int ncols);
void da_fill_A(void *daobj, int row, int col, double val);
void da_fill_A_col(void *daobj, int col, double *vals);
int da_check_A(void *daobj);
void da_writeA(void *daobj, const char *fname);
void da_create_maps(void *daobj, double **gx, double **gy, double *gz,
 		                                         int nx, int ny, int nz);
void da_create_maps_ml(void *daobj, double *igx, double *igy, double *gz,
		       int nx, int ny, int nz);
void da_fill_map(void *daobj, int ls, int i, int j, int k);
void da_fill_map_ml(void *daobj, int ns, int *i, int *j, int *k);
void da_add_state(void *daobj, const char *str, int off);
void *da_create_fixed_obs(char *name, 
			    char *file, 
			    char *tunits,
			    double err,
			    double dt);
void *da_create_glider_obs(char *name, 
			   char *file, 
			   char *tunits,
			   double err,
			   double dt);
void *da_create_2d_field_obs(char *name, 
			     char *file, 
			     char *tunits,
			     double err);
void da_fixed_obs_set_location(void *obs, void *obj,
				 double lon, double lat, double dep);
void da_glider_obs_set_location(void *obs, char *lon, char *lat, char *dep);
void da_do_analysis(void *daobj);
void da_set_wb(void *daobj, int index, double data);
void da_set_xy(void *daobj, int index, double x, double y);
void da_set_wb_col(void *daobj, int ndata, double *data);
double da_get_wa(void *daobj, int index);
double da_get_wo_val(void *daobj, int index);
double da_get_wo_err(void *daobj, int index);
void da_fill_wa(void *daobj, int ndata, double *data);
void da_fill_wo_val(void *daobj, int ndata, double *data);
void da_fill_wo_err(void *daobj, int ndata, double *data);
int da_read_all_obs(void *daobj, double t);
void *da_get_obs_ptr(void *daobj, int i);
const char *da_obs_get_name(void *daobs);
double da_obs_get_error(void *daobs);
double da_obs_get_dt(void *daobs);
const char* da_obs_get_type(void *daobs);
const char *da_obs_get_fname(void* daobs);
const char *da_obs_get_location(void *daobs);
const char *da_obs_get_vars(void *daobs);
void da_obs_add_var(void *daobj, void *daobs, char *var, char *ref_var,
		                                             int is3D);
void da_add_obs(void *daobj, void *obs);
void da_copy_wb_wa(void *daobj);

#endif

/* end of da_interface.h */
