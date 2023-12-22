/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/dfcoords.c
 *
 *  \brief Datafile (df) coordinate routines
 *
 *  Builds the association between variables and coordinates.
 *  Only a limited set of coordinate mappings are currently
 *  supported
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: dfcoords.c 7422 2023-10-05 02:15:37Z her127 $
 */

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "ems.h"

#define DEG2RAD (M_PI/180.0)


/* Prototypes for local routines */
void decode_coords(datafile_t *df, df_variable_t *v, char *coords);
void decode_coord_type(char *ctype, VariableType * type,
                       char **coord_domain);
void get_1d_coords(datafile_t *df, int nc, int *coordids, int *nc_1d,
                   int **cids_1d, int *nc_m, int **cids_m);
int search_for_coord(datafile_t *df, VariableType type, int nc,
                     int *coordids);
int df_map_dependent_coords(datafile_t *df, df_coord_mapping_t *cm, int nc,
                            int *coordids);
int are_vars_same_dims(datafile_t *df, int vid1, int vid2);
char *get_text_attribute(datafile_t *df, int varid, const char *attrib);
void free_coord_system(datafile_t *df, df_coord_system_t *csystem);
void set_1d_coords(datafile_t *df, df_variable_t *v,
                   df_coord_mapping_t *cm);
void set_1d_multi_coords(datafile_t *df, df_variable_t *v,
                         df_coord_mapping_t *cm);
void set_2d_coords(datafile_t *df, df_variable_t *v,
                   df_coord_mapping_t *cm);
void set_geo_transform(datafile_t *df, df_coord_mapping_t *cm);
df_coord_transform_t *create_geo_transform(int int_geotype,
                                           char *int_proj, int ext_geotype,
                                           char *ext_proj);

void proj_forward(void *d, int n, double from[], double to[]);
void proj_inverse(void *d, int n, double from[], double to[]);
void proj_free(void *d);


void cm_1d_free(datafile_t *df, df_coord_mapping_t *cm);
int cm_1d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
               const int depindicies[], const double coords[],
               double indices[]);
int cm_1d_itoc(datafile_t *df, df_coord_mapping_t *cm,
               const int depindicies[], const double indices[],
               double coords[]);
void cm_rect_2d_free(datafile_t *df, df_coord_mapping_t *cm);
int cm_rect_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                    const int depindicies[], const double coords[],
                    double indices[]);
int cm_rect_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                    const int depindicies[], const double indices[],
                    double coords[]);
void cm_polar_2d_free(datafile_t *df, df_coord_mapping_t *cm);
int cm_polar_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double coords[],
                     double indices[]);
int cm_polar_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double indices[],
                     double coords[]);
void cm_bilinear_2d_free(datafile_t *df, df_coord_mapping_t *cm);
int cm_bilinear_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                        const int depindicies[], const double coords[],
                        double indices[]);
int cm_bilinear_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                        const int depindicies[], const double indices[],
                        double coords[]);
void cm_1d_multi_free(datafile_t *df, df_coord_mapping_t *cm);
int cm_1d_multi_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double coords[],
                     double indices[]);
int cm_1d_multi_itoc(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double indices[],
                     double coords[]);

double interp_1d_inv_weight_hashed(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_1d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_2d_inv_weight_hashed(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_2d_inv_weight(datafile_t *df, df_variable_t *v, int record,
                            double coords[]);
double interp_linear(datafile_t *df, df_variable_t *v, int record,
                     double coords[]);
double interp_linear_flagged(datafile_t *df, df_variable_t *v, int record,
			     double coords[]);
double interp_linear_filled(datafile_t *df, df_variable_t *v, int record,
			    double coords[]);
double interp_linear_bathy(datafile_t *df, df_variable_t *v, int record,
			   double coords[]);
double interp_linear_degrees(datafile_t *df, df_variable_t *v, int record,
			     double coords[]);


int df_default_hashtable_size = 0;
int df_default_geotype = GT_NONE;
char *df_default_projection = NULL;

/** Associates a coordinate system with a particular variable.
  *
  * @param df pointer to datafile structure (assumed
  * previously initialised)
  * @param v datafile variable pointer
  * @param ncm Number of coordinate mapping requests.
  * @param cmaps An array of coordinate mapping requests.
  * @return non-zero if successful.
  *
  * @see warn() User warning on non-fatal error.
  * @see quit() User warning on fatal error.
  */
int df_set_coord_system(datafile_t *df, df_variable_t *v,
                        int ncm, df_coord_mapping_t *cmaps)
{
  int i, j, k, n;
  df_variable_t *cv = NULL;
  int *coordids = NULL;
  int *coordtypes = NULL;
  int nc = 0;

  if (v->csystem != NULL) {
    free_coord_system(df, v->csystem);
    v->csystem = NULL;
  }

  if (ncm < 0)
    quit("df_set_coord_system: Invalid number of coordinates requests.\n");

  else if (ncm == 0)            /* Clear the coordinate and exit */
    return 1;

  /* Count the number of coordinates and verify they are appropriately
     dimensioned */
  for (i = 0; i < ncm; ++i) {
    if (cmaps[i].nc <= 0)
      quit
        ("df_set_coord_system: Invalid number of coordinates specified.\n");
    nc += cmaps[i].nc;
  }

  /* Recast the coordinate identifiers and types into single arrays. This
     makes it easier to check for duplicates, check for existance of
     variables, and defines the order */
  coordids = (int *)malloc(sizeof(int) * nc);
  coordtypes = (int *)malloc(sizeof(int) * nc);
  for (i = 0, k = 0; i < ncm; ++i) {
    for (j = 0; j < cmaps[i].nc; ++j, ++k) {
      coordids[k] = cmaps[i].coordids[j];
      coordtypes[k] = cmaps[i].coordtypes[j];
    }
  }

  /* Check for duplicate coordinates */
  for (i = 0; i < nc; ++i) {
    for (j = i + 1; j < nc; ++j) {
      if (coordids[i] == coordids[j]) {
        warn
          ("df_set_coord_system: Duplicate coordinate variables specified.\n");
        return 0;
      }

      if (coordtypes[i] == coordtypes[j]) {
        warn
          ("df_set_coord_system: Duplicate coordinate types specified.\n");
        return 0;
      }
    }
  }


  /* Verify that the coord variable exists, and that it's dimensions are a
     subset of the data variable coordinates. All remaining dimensions are
     non-coordinate dimensions. */
  for (i = 0; i < nc; ++i) {

    /* Range check */
    if ((coordids[i] < 0) || (coordids[i] >= df->nv))
      quit("df_set_coord_system: Invalid coordinate id.\n");
    cv = &df->variables[coordids[i]];

    /* Check the dimensions */
    for (j = 0; j < cv->nd; ++j) {
      int not_found = 1;
      for (k = 0; (k < v->nd) && (not_found); ++k) {
        if (cv->dimids[j] == v->dimids[k])
          not_found = 0;
      }

      if (not_found) {
        warn
          ("df_set_coord_system: The dimension '%s' was not found in variable '%s'\n",
           df->dimensions[cv->dimids[j]].name, v->name);
        return 0;
      }
    }

  }

  /* Find the non-coordinate dimensions. Non-coordinate dimensions are
     currently NOT supported. */
  v->nncd = 0;
  v->ncdimids = (int *)malloc(v->nd * sizeof(int));

  for (i = 0; i < v->nd; ++i) {
    int found = 0;
    for (j = 0; (j < nc) && (!found); ++j) {
      cv = &df->variables[coordids[j]];
      for (k = 0; (k < cv->nd) && (!found); ++k) {
        if (v->dimids[i] == cv->dimids[k])
          found = 1;
      }
    }
    if (!found) {
      v->ncdimids[v->nncd++] = v->dimids[i];
      warn
        ("df_set_coord_system: Non-coordinate dimension '%s' for '%s'.\n",
         df->dimensions[v->dimids[i]].name, v->name);
    }
  }

  if (v->nncd != 0) {
    warn
      ("df_set_coord_system: Non-coordinate dimensions are not permitted.\n");
    return 0;
  }

  /* Confirm that for each coordinate map requests, the coordinates have
     the same number of dimensions and that all dimensions are same and in
     the same order. */
  for (i = 0; i < ncm; ++i) {
    df_coord_mapping_t *mp = &cmaps[i];
    df_variable_t *fv = &df->variables[mp->coordids[0]];

    for (j = 1; j < mp->nc; ++j) {
      df_variable_t *sv = &df->variables[mp->coordids[j]];
      if (fv->nd == sv->nd) {
        for (k = 0; k < fv->nd; ++k) {
          if (fv->dimids[k] != sv->dimids[k])
            quit
              ("df_set_coord_system: Coordinate variables have inconsistant dimensions.\n");
        }
      } else
        quit
          ("df_set_coord_system: Inconsistant number of coordinates in mapping request.\n");
    }

  }


  /* For each coordinate map determine the number of independent
     dimensions (i.e. can be evaluted from the coordinates provided) and
     the number of dependent dimensions (must be provided from another
     coordinate map). */
  for (i = 0; i < ncm; ++i) {
    df_coord_mapping_t *mp = &cmaps[i];
    df_variable_t *fv = &df->variables[mp->coordids[0]];
    mp->nd = 0;
    mp->dimids = (int *)malloc(sizeof(int) * fv->nd);
    mp->ndepd = 0;
    mp->depdimids = (int *)malloc(sizeof(int) * fv->nd);

    /* For each dimension find the coordinate map that can resolve this
       dimension. */
    for (j = 0; j < fv->nd; ++j) {
      int found = 0;
      for (k = 0; (k < ncm) && (!found); ++k) {
        if (k != i) {
          df_variable_t *vk = &df->variables[cmaps[k].coordids[0]];
          for (n = 0; (n < vk->nd) && (!found); ++n) {
            if (fv->dimids[j] == vk->dimids[n]) {
              found = 1;
            }
          }
        }
      }

      if (found)
        mp->depdimids[mp->ndepd++] = fv->dimids[j];
      else
        mp->dimids[mp->nd++] = fv->dimids[j];

    }

    /* The number of dimensions that can be evaluated by the coordinates
       must be the same as the number of coordinates */
    if (mp->nd > mp->nc)
      quit
        ("df_set_coord_system: More dimensions than coordinates specified.\n");
  }


  /* Create the coordinate maps from the requests, and map to the
     appropriate functions if supported. */

  /* Allocate the memory for the coordinate system */
  v->csystem = (df_coord_system_t *)malloc(sizeof(df_coord_system_t));
  memset(v->csystem, 0, sizeof(df_coord_system_t));
  v->csystem->nc = nc;
  v->csystem->coordids = coordids;
  v->csystem->coordtypes = coordtypes;
  v->csystem->ncm = ncm;
  /* Duplicate the requests */
  v->csystem->cmaps =
    (df_coord_mapping_t *)malloc(sizeof(df_coord_mapping_t) * ncm);
  memset(v->csystem->cmaps, 0, sizeof(df_coord_mapping_t) * ncm);
  for (i = 0; i < ncm; ++i) {
    df_coord_mapping_t *cm = &v->csystem->cmaps[i];
    cm->nc = cmaps[i].nc;
    cm->coordids = (int *)malloc(sizeof(int) * cm->nc);
    cm->coordtypes = (int *)malloc(sizeof(VariableType) * cm->nc);
    for (j = 0; j < cm->nc; ++j) {
      cm->coordids[j] = cmaps[i].coordids[j];
      cm->coordtypes[j] = cmaps[i].coordtypes[j];
    }


    cm->nd = cmaps[i].nd;
    cm->dimids = cmaps[i].dimids;

    cm->ndepd = cmaps[i].ndepd;
    cm->depdimids = cmaps[i].depdimids;
    cm->free = NULL;
    cm->coords_to_indices = NULL;
    cm->indices_to_coords = NULL;
    cm->special_data = NULL;

    /* If there are the same number of coordinates as dimensions, * then
       we assume that we have a grid of some sort. */
    if (cm->nc == cm->nd) {
      switch (cm->nc) {
      case 1:
        set_1d_coords(df, v, cm);
        break;

      case 2:
	set_2d_coords(df, v, cm);
        break;

      default:
        quit
          ("df_set_coord_system: Multi dimension coordinate system not currently supported.\n");
        break;
      }
    }

    else if ((cm->nc > cm->nd)) {
      if (cm->nd == 1)
        set_1d_multi_coords(df, v, cm);
      else
        quit
          ("df_set_coord_system: Multi coordinate, single dimension coordinate system not currently supported\n");
    }

    else {
      warn("df_set_coord_system: Less coordinates than dimensions!\n");
      return 0;
    }
  }

  /* Define the order in which the coordinate mappings should be
     evaluated. */
  v->csystem->cmorder = (int *)malloc(sizeof(int) * ncm);
  j = 0;
  for (i = 0; i < ncm; ++i) {
    if (!v->csystem->cmaps[i].ndepd)
      v->csystem->cmorder[j++] = i;
  }
  for (i = 0; i < ncm; ++i) {
    if (v->csystem->cmaps[i].ndepd)
      v->csystem->cmorder[j++] = i;
  }

  return 1;
}


/** From attribute information encoded in the datafile, attempts to
  * infer the coordinate mappings and set the coordinate system.
  *
  * @param df pointer to the datafile structure.
  * @param v pointer to datafile variable.
  * @return non-zero value if coordinate system was set.
  */
int df_infer_coord_system(datafile_t *df, df_variable_t *v)
{
  int i, j, k;
  int failed = 0;
  df_variable_t *cv = NULL;
  int nncd = 0;
  int *ncdimids = NULL;
  int nc_1d;
  int *cids_1d;
  int nc_m;
  int *cids_m;
  int ncm;
  int nc, nd;
  df_coord_mapping_t cmaps[MAXNUMCMAPS];
  df_coord_mapping_t *cm;
  int coordids[MAXNUMCOORDS];
  int inferred = 0;

  /* No coordinate system required if 0 dimensional */
  if (v->nd == 0)
    return 1;

  if (v->nc <= 0) {
    warn
      ("df_infer_coord_system: No coordinate information was provided\n");
    return 0;
  }

  /* Check for duplicate coordinates */
  for (i = 0; i < v->nc; ++i) {
    for (j = i + 1; j < v->nc; ++j) {
      if (v->coordids[i] == v->coordids[j]) {
        warn
          ("df_infer_coord_system: Duplicate coordinate variables specified\n");
        return 0;
      }
    }
  }

  /* Make a copy of the coordinate ids. If one of the coordinate id's is
     the record id, then strip it out. The record id cannot be a
     coordinate because it already is. It is implied in the definition of
     the record. */
  nc = 0;
  for (i = 0; i < v->nc; ++i) {
    if (!((v->dim_as_record) && (v->coordids[i] == df->ri)))
      coordids[nc++] = v->coordids[i];

  }

  /* Verify that the coord variable exists, that the variable is a known
     coordinate type, that it's dimensions are a subset of the data
     variable coordinates, and determine which dimensions are free. */
  for (i = 0; i < nc; ++i) {

    /* Range check */
    if ((coordids[i] < 0) || (coordids[i] >= df->nv)) {
      warn("df_infer_coord_system: Invalid coordinate id specified.\n");
      return 0;
    }
    cv = &df->variables[coordids[i]];

    /* Coordinate type check */
    switch (cv->type) {
    case VT_UNKNOWN_COORD:
    case VT_DATA:
      warn
        ("df_infer_coord_system: '%s' is a data variable or an unsupport coordinate variable: variable=%s file=%s\n",
         cv->name, v->name, df->name);
      return 0;

    default:
      break;
    }
    if (cv->type & VT_INFERRED) inferred += 1;

    /* Check the dimensions */
    for (j = 0; j < cv->nd; ++j) {
      int not_found = 1;
      for (k = 0; (k < v->nd) && (not_found); ++k) {
        if (cv->dimids[j] == v->dimids[k])
          not_found = 0;
      }

      if (not_found) {
        warn
          ("df_infer_coord_system: The dimension '%s' was not found in variable '%s'\n",
           df->dimensions[cv->dimids[j]].name, v->name);
        return 0;
      }
    }

  }

  /* Determine which dimensions do not have a coordinate associated with
     them. These are called non-coordinate dimensions and are currently
     not supported. */
  nncd = 0;
  ncdimids = (int *)malloc(v->nd * sizeof(int));
  for (i = 0; i < v->nd; ++i) {
    int found = 0;
    for (j = 0; (j < nc) && (!found); ++j) {
      cv = &df->variables[coordids[j]];
      for (k = 0; (k < cv->nd) && (!found); ++k) {
        if (v->dimids[i] == cv->dimids[k]) {
          found = 1;
	}
      }
    }
    if (!found)
      ncdimids[nncd++] = v->dimids[i];
  }
  free(ncdimids);
  if (nncd != 0) {
    warn
      ("df_infer_coord_system: Some of '%s' dimensions have NO coordinates associated with them.\nThis is not currently supported.\n",
       v->name);
    return 0;
  }


  /* This is the fun bit. Check whether the coordinates form a known
     coordinate system type. */

  /* Create a list of coordinates that depend on one dimension only. */
  get_1d_coords(df, nc, coordids, &nc_1d, &cids_1d, &nc_m, &cids_m);

  /* Build the one dimensional coordinate mappings.  */
  ncm = 0;
  if (nc_1d > 0) {
    df_coord_mapping_t *depcm = &cmaps[ncm];
    if (df_map_dependent_coords(df, depcm, nc_1d, cids_1d))
      ++ncm;
    else
      depcm = NULL;

    for (i = 0; i < nc_1d; ++i) {
      /* Check if this the id was one used in the depcm */
      int used = 0;

      if (depcm != NULL) {
        for (j = 0; j < depcm->nc; ++j) {
          if (cids_1d[i] == depcm->coordids[j]) {
            used = 1;
            break;
          }
        }
      }

      if (!used) {
        cm = &cmaps[ncm];
        cm->nc = 1;
        cm->coordids = (int *)malloc(sizeof(int));
        cm->coordtypes = (int *)malloc(sizeof(int));
        cm->coordids[0] = cids_1d[i];
        cm->coordtypes[0] = df->variables[cids_1d[i]].type;
        ++ncm;
      }
    }
  }

  if (nc_m > 0) {
    if (df_map_dependent_coords(df, &cmaps[ncm], nc_m, cids_m)) {
      ++ncm;
    } else
      failed = 1;
  }

  free(cids_1d);
  free(cids_m);

  if (!failed) {
    /* Now set the coordinate system */
    df_set_coord_system(df, v, ncm, cmaps);

    /* Free memory associated with coordinate mapping */
    for (i = 0; i < ncm; ++i) {
      free(cmaps[i].coordtypes);
      free(cmaps[i].coordids);
      cmaps[i].coordids = NULL;
      cmaps[i].coordtypes = NULL;
    }
  }

  /* Unsupported */
  else {
    warn("dfInferCoordSystem: The coordinate system is unsupported.\n");
    return 0;
  }

  /* Set the interpolation method */
  nd = df_get_num_dims(df, v);
  nc = df_get_num_coords(df, v);
  /* Are there the same number of coordinates as dimensions. We shall
     assume a regular grid. */
  if (nc == nd) {
    if (df_is_irule(df)) {
      if (df_get_num_coords(df, v) == 3 ) {
	if (df_is_sigma(df)) df_set_sigma(df, v);
	if (df_is_zstar(df)) df_set_zstar(df, v);
      }
      if (strcmp(df->i_rule, "nearest_eps") == 0) {
	/* Nearest neighbour for 'simple' files */
	if (df_get_num_coords(df, v) == 2 )
	  v->interp = interp2d_nearest_within_eps;
	else if (df_get_num_coords(df, v) == 3 )
	  v->interp = interp3d_nearest_within_eps;
      } else {
	if (df_get_num_coords(df, v) == 2 ) {	
	  v->gs0_master_2d = v->gs1_master_2d = NULL;
	  if (v->type & VT_INFERRED)
	    v->interp = interp_us_2d_i;
	  else
	    v->interp = interp_us_2d_c;
	} else if (df_get_num_coords(df, v) == 3 ) {
	  v->gs0_master_3d = v->gs1_master_3d = NULL;
	  if (v->type & VT_INFERRED)
	    v->interp = interp_us_3d_i;
	  else
	    v->interp = interp_us_3d_c;
	}
      }
    } else if ((v->lflag && v->lflag != v->add_offset) || 
	       (v->cflag && v->cflag != v->add_offset))
      /* land or cloud flags exist */  
      v->interp = interp_linear_flagged;
    else if (v->type & VT_BATHY)
      v->interp = interp_linear_bathy;
    /*else if (v->type & VT_INFERRED && v->missing && v->missing != v->add_offset)*/
    else if (v->type & VT_INFERRED && v->type & (VT_MV|VT_FV))
      /*else if (inferred == nc && v->missing && v->missing != v->add_offset)*/
      /* missing or fill_value exists */
      v->interp = interp_linear_filled;
    else if (v->type & VT_INFERRED && v->fillvalue && v->fillvalue != v->add_offset)
      /*else if (inferred == nc && v->missing && v->missing != v->add_offset)*/
      /* missing_value exists */
    v->interp = interp_linear_filled;
    else if (v->type & VT_COORD && v->fillvalue && v->fillvalue != v->add_offset)
      /* FillValue exists */
      v->interp = interp_linear_flagged;
    /*v->interp = interp_linear;*/
    else if (v->units && strcmp(v->units, "degrees") == 0)
      v->interp = interp_linear_degrees;
    else {
      if (df_is_sigma(df)) df_set_sigma(df, v);
      if (df_is_zstar(df)) df_set_zstar(df, v);
      v->interp = interp_linear;
    }
  }
  /* 
   * Are there more coordinate than dimensions, and is there only one
   * dimension. If so, then we can assume than we shall interpolate
   * using an inverse weighting scheme. 
   *
   * Note 01/10 (FR) : Changed these to use the _hashed versions of
   *                   the inverse weigthing schemes. The original
   *                   versions are still left in the code for comparision
   */
  else if (nc > nd) {
    if (df->type == DFT_MEMPACK) {
      v->interp = interp_linear;
    }
    else if (df_is_ugrid(df)) {
      /* Only one method supported */
      int *coordtypes = df_get_coord_types(df, v);
      /*
      if (df_get_num_coords(df, v) == 2 &&
	  (coordtypes[0]|coordtypes[1]) == (VT_LONGITUDE|VT_LATITUDE) ) {
	v->interp = interp_nearest_within_eps;
      */

      if (df_get_num_coords(df, v) == 2 ) {
	    v->gs0_master_2d = v->gs1_master_2d = NULL;
	    v->interp = interp_us_2d;
      } else if (df_get_num_coords(df, v) == 3 ) {
	    v->gs0_master_3d = v->gs1_master_3d = NULL;
	    v->interp = interp_us_3d;
      } else
	quit("df_infer_coord_system: Can't infer coordinates for %s in UGRID file %s\n", v->name, df->name);
    }
    else {
      if (nd == 1)
	if (df_is_sigma(df)) {
	  df_set_sigma(df, v);
          v->interp = interp2d_nearest;
	} else
          v->interp = interp_1d_inv_weight_hashed;
      else if (nd == 2) {
	/*if (df_is_zstar(df)) df_set_zstar(df, v);*/
	if (df_is_sigma(df)) {
	  df_set_sigma(df, v);
          v->interp = interp3d_nearest;
	} else
          v->interp = interp_2d_inv_weight_hashed;
      } else {
          warn("df_infer_coord_system: File %s, variable %s : \n",
	       df->name, v->name);
          quit("Unable to interpolate data with %d coordinate and %d dimensions.\n",
	       nc, nd);
      }
    }
  } else {
    warn("df_infer_coord_system: File %s, variable %s : \n",
	 df->name, v->name);
    quit("Less coordinates than dimensions.\n");
  }

  return 1;
}


/** Set the default geographic data type.
  *
  * @param type character string containing the geotype.
  * NULL - no prefered output form, no transformation occurs,
  * data same as file format.
  * "geographic" - All 'geographic' related data should be
  * returned as latitudes and longitudes.
  * "proj=xxx .." - Standard map_proj_t format.
  */
void df_set_default_proj_type(char *type)
{

  if (type == NULL)
    df_default_geotype = GT_NONE;
  else if (strncasecmp(type, GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) == 0)
    df_default_geotype = GT_GEOGRAPHIC;
  else
    df_default_geotype = GT_PROJECTION;

  df_default_projection = type;
}

/** Set the default hash_table_t size for evaluating numeric grids.
  * @param hsize Hashtable size.
  */
void df_set_default_hashtable_size(int hsize)
{

  /* set the value to next power of two */
 /* if(hsize > 0)
  {
    int i=0;
    int s=hsize;
    while(s > 0)
    {
      s = s >> 1;
      i++;
    }
    df_default_hashtable_size = ldexp(1,i+1);

    emstag(LTRACE,"lib:io:dfcoords:df_set_default_hashtable_size","Setting default hashtable size to %d as %d ",&hsize,&df_default_hashtable_size );
  }
  else
   df_default_hashtable_size = 0;
  */
   df_default_hashtable_size = hsize;
   emstag(LTRACE,"lib:io:dfcoords:df_set_default_hashtable_size","Setting default hashtable size to %d ",&hsize );

}


/** Set the geographic data type for a specific datafile.
  *
  * @param df pointer to datafile structure.
  * @param type character string containing the geotype.
  * NULL - no prefered output form, no transformation occurs,
  * data same as file format.
  * "geographic" - All 'geographic' related data should be
  * returned as latitudes and longitudes.
  * "proj=xxx .." - Standard map_proj_t format.
  */
void df_set_proj_type(datafile_t *df, char *type)
{

  if (type == NULL)
    df->geotype = GT_NONE;
  else if (strcasecmp(type, "geographic") == 0)
    df->geotype = GT_GEOGRAPHIC;
  else
    df->geotype = GT_PROJECTION;

  df->projection = type;
}

/** Get the number of coordinates for the specified variable,
  * as defined in the coordinate system structure.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to datafile variable.
  * @return Number of coordiantes or 0 if no coordinate system..
  * @see dfGetCoordIds
  * @see dfGetCoordTypes
  */
int df_get_num_coords(datafile_t *df, df_variable_t *v)
{
  return (v->csystem != NULL) ? v->csystem->nc : 0;
}


/** Get the coordinate ids for the specified variable.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to datafile variable.
  * @return pointer to array of coordinates ids. NULL if none.
  * @see dfGetNumCoordinates
  * @see dfGetCoordTypes
  */
int *df_get_coord_ids(datafile_t *df, df_variable_t *v)
{
  return (v->csystem != NULL) ? v->csystem->coordids : NULL;
}


/** Get the coordinate types for the specified variable.
  *
  * @param df pointer to datafile structure.
  * @param v pointer to datafile variable.
  * @return pointer to array of coordinates types. NULL if none.
  * @see dfGetNumCoordinates
  * @see dfGetCoordIds
  */
int *df_get_coord_types(datafile_t *df, df_variable_t *v)
{
  return (v->csystem != NULL) ? v->csystem->coordtypes : NULL;
}



/** Convert the specified coordinates to indice values.
  *
  * @param df pointer to the datafile structure.
  * @param v pointer to the datafile variable.
  * @param coords Coordinate values in order of original coord system.
  * @param indices Index values in order of dimensions.
  * @returns Zero if failed, else number of indices.
  */
int df_ctoi(datafile_t *df, df_variable_t *v,
            const double coords[], double indices[])
{
  int i, j, k;
  int depindices[MAXNUMDIMS];
  double tcoords[MAXNUMDIMS];
  double ttcoords[MAXNUMDIMS];
  double tindices[MAXNUMDIMS];
  df_coord_system_t *cs = v->csystem;

  if (cs == NULL)
    quit("df_ctoi: No coordinate system specified for variable '%s'\n",
         v->name);

  /* Step through each of the coordinate mappings. Evaluate according to
     the order specified by the coordinate system. */
  for (i = 0; i < cs->ncm; ++i) {
    df_coord_mapping_t *cm = &cs->cmaps[cs->cmorder[i]];
    /* If this coordinate mapping is dependent on dimensions previously
       computed, then populate the depindices array with their values */
    k = 0;
    for (j = 0; (j < v->nd) && (k < cm->ndepd); ++j) {
      if (v->dimids[j] == cm->depdimids[k]) {
        depindices[k++] = (int)indices[j];
      }
    }

    /* Copy the coordinate system coords into the correct positions in the
       coordinate mapping coord array */
    k = 0;
    for (j = 0; (j < cs->nc) && (k < cm->nc); ++j) {
      if (cs->coordids[j] == cm->coordids[k])
        tcoords[k++] = coords[j];
    }

    /* Perform a foward coordinate transform (to file) if required. */
    if (cm->transform != NULL) {
      cm->transform->forward(cm->transform->private_data,
                             cm->nc, tcoords, ttcoords);
      for (j = 0; j < cm->nc; ++j) {
        tcoords[j] = ttcoords[j];
      }
    }
    if (cm->coords_to_indices(df, cm, depindices, tcoords, tindices)) {
      /* Copy the coordinate map indices into the correct positions in the
         coordinate system indice array */
      k = 0;
      for (j = 0; (j < v->nd) && (k < cm->nd); ++j) {
        if (v->dimids[j] == cm->dimids[k])
          indices[j] = tindices[k++];
      }
    } else
      return 0;
  }

  return v->nd;
}


/** Convert the indice values to coordinate values.
  *
  * @param df pointer to the datafile structure.
  * @param v pointer to the datafile variable.
  * @param indices Index values in order of dimensions.
  * @param coords Coordinate values in order of original coord system.
  * @returns Zero if failed, else number of coordinates.
  */
int df_itoc(datafile_t *df, df_variable_t *v,
            const double indices[], double coords[])
{
  int i, j, k;
  int depindices[MAXNUMDIMS];
  double tcoords[MAXNUMDIMS];
  double ttcoords[MAXNUMDIMS];
  double tindices[MAXNUMDIMS];
  df_coord_system_t *cs = v->csystem;

  if (cs == NULL)
    quit("df_itoc: No coordinate system specified for variable '%s'\n",
         v->name);

  /* Step through each of the coordinate mappings. Evaluate according to
     the order specified by the coordinate system. */
  for (i = 0; i < cs->ncm; ++i) {
    df_coord_mapping_t *cm = &cs->cmaps[cs->cmorder[i]];

    /* If this coordinate mapping is dependent on dimensions previously
       computed, then populate the depindices array with their values */
    k = 0;
    for (j = 0; (j < v->nd) && (k < cm->ndepd); ++j) {
      if (v->dimids[j] == cm->depdimids[k])
        depindices[k++] = (int)indices[j];
    }


    /* Copy the coordinate map indices into the correct positions in the
       coordinate system indice array */
    k = 0;
    for (j = 0; (j < v->nd) && (k < cm->nd); ++j) {
      if (v->dimids[j] == cm->dimids[k])
        tindices[k++] = indices[j];
    }

    if (cm->indices_to_coords(df, cm, depindices, tindices, tcoords)) {
      /* Perform an inverse coordinate transform (to user) if required. */
      if (cm->transform != NULL) {
        cm->transform->inverse(cm->transform->private_data,
                               cm->nc, tcoords, ttcoords);
        for (j = 0; j < cm->nc; ++j)
          tcoords[j] = ttcoords[j];
      }

      /* Copy the coordinate system coords into the correct positions in
         the coordinate mapping coord array */
      k = 0;
      for (j = 0; (j < cs->nc) && (k < cm->nc); ++j) {
        if (cs->coordids[j] == cm->coordids[k])
          coords[j] = tcoords[k++];
      }
    } else
      return 0;

  }

  return v->nd;
}


/***** PRIVATE or PROTECTED functions */

/*********************************************************************
Decode a string of coordinates, and populates the coordids array
within the variable structure. Also checks to see if there are any
one dimensional mappings.

Arguments:
    df		-   datafile_t structure.
    v           -   df_variable_t which coords are associated with.
    coords	-   A string of coordinate separated by commas,
		    space or :.

This routine calls quit() if anything goes wrong
*********************************************************************/
void decode_coords(datafile_t *df, df_variable_t *v, char *coords)
{
  int i, j, k;
  char line[MAXLINELEN];
  int maxnc = v->nd + 1;        /* Maximum number of coordinates */
  int *coordids;
  int nc = 0;
  int length = 0;


  /* Copy the coordinate line, nullify all punctuation, and count the
     number of expected names. */
  if (coords != NULL) {
    length = strlen(coords);
    strcpy(line, coords);
    for (i = 0; i < length; ++i) {
      /* If space,tab or comma then delimits a variable */
      if (isspace((int)line[i]) || (line[i] == ',') || (line[i] == ':')) {
        /* Count only when a previous NULL was not evident */
        if ((i > 0) && (line[i - 1] != 0))
          ++maxnc;
        line[i] = 0;
      }
    }
  }

  /* Allocate the array, and copy the charaters */
  coordids = (int *)malloc(maxnc * sizeof(int));
  memset(coordids, 0, maxnc * sizeof(int));

  /* Now check and match the coordinates specified in the * string against
     variable names. Warn if the variable is not found.  Remove
     duplicates. */
  if (coords != NULL) {
    for (i = 0; i < length;) {

      /* Skip nulls to first valid character */
      while ((i < length) && !line[i])
        ++i;

      if (i < length) {
        int not_found = 1;

        /* Check the name against ALL variable names */
        for (j = 0; (j < df->nv) && (not_found); ++j) {
          if (strcasecmp(df->variables[j].name, &line[i]) == 0) {
            int duplicated = 0;
            /* Now check this isn't already duplicated */
            for (k = 0; (k < nc) && (!duplicated); ++k) {
              if (coordids[k] == j)
                duplicated = 1;
            }

            if (!duplicated)
              coordids[nc++] = j;
            not_found = 0;
          }
        }

        if (not_found) {
          warn
            ("decode_coords: Unable to find coordinate '%s' for variable '%s'.\n",
             &line[i], v->name);
        }

        /* Skip to next null */
        while ((i < length) && line[i])
          ++i;
      }
    }
  }

  /* copy to the arguments */
  v->nc = nc;
  if (nc == 0) {
    free(coordids);
    v->coordids = NULL;
  } else
    v->coordids = coordids;

}



/*********************************************************************
Routine to read an attribute from the netCDF file.

Arguments:
    ctype	-    String of coordinate type.
    type	-    df_variable_t type.
    coord_domain-    What's left after the string is parsed.

This routine calls quit() if anything goes wrong
*********************************************************************/
void decode_coord_type(char *ctype, VariableType * type,
                       char **coord_domain)
{
  int i;
  int sep = -1;
  int length;

  /* Search for the first separator. If found then split the string into
     two */
  for (i = 0; ctype[i] && (sep < 0); ++i)
    if (ctype[i] == ':')
      sep = i;

  *type = VT_UNKNOWN_COORD;
  length = (sep > 0) ? sep : (int)strlen(ctype);

  if (strncasecmp(ctype, "time", length) == 0)
    *type = VT_TIME;

  else if (strncasecmp(ctype, "t", length) == 0)
    *type = VT_TIME;

  else if (strncasecmp(ctype, "x", length) == 0)
    *type = VT_X;

  else if (strncasecmp(ctype, "y", length) == 0)
    *type = VT_Y;

  else if (strncasecmp(ctype, "z", length) == 0)
    *type = VT_Z;

  else if (strncasecmp(ctype, "latitude", length) == 0)
    *type = VT_LATITUDE;

  else if (strncasecmp(ctype, "longitude", length) == 0)
    *type = VT_LONGITUDE;

  if (*type == VT_UNKNOWN_COORD)
    warn("decode_coord_type: Unknown coordinate type '%s'\n", ctype);

  if (sep >= 0) {
    *coord_domain = (char *)malloc(sizeof(char *) * MAXLINELEN);
    strcpy(*coord_domain, &ctype[sep + 1]);
  } else
    *coord_domain = NULL;

}


/* Free the data associated with the coordinate system.
  *
  * @param df pointer to a datafile.
  * @param csystem pointer to a coordinate system.
  */
void free_coord_system(datafile_t *df, df_coord_system_t *csystem)
{
  if (csystem != NULL) {
    int i;
    for (i = 0; i < csystem->ncm; ++i) {
      if (csystem->cmaps[i].free != NULL)
        csystem->cmaps[i].free(df, &csystem->cmaps[i]);

      if (csystem->cmaps[i].depdimids != NULL)
        free(csystem->cmaps[i].depdimids);

      if (csystem->cmaps[i].dimids != NULL)
        free(csystem->cmaps[i].dimids);

      if (csystem->cmaps[i].coordids != NULL)
        free(csystem->cmaps[i].coordids);
    }

    if (csystem->ncm > 0) {
      free(csystem->cmaps);
      free(csystem->cmorder);
    }
    free(csystem);
  }
}



/*
Locate the coordinate system that contains the specified indices in
the coordinate system table.

Arguments:
    df          -    Data file.
    nc          -    Number of coordinates.
    coordids    -    list of coordinate identifiers.
    nc_1d       -    number of 1d coord ids.
    cids_1d     -    list of 1d coordinates.
    nc_m        -    number of coord ids with multiple coords.
    cids_m      -    list of multiple coordinates.

This routine calls quit() if anything goes wrong
*/
void get_1d_coords(datafile_t *df, int nc, int *coordids,
                   int *nc_1d, int **cids_1d, int *nc_m, int **cids_m)
{
  int i;
  int *oned, *multi;
  int nc1d = 0, ncm = 0;

  oned = (int *)malloc(nc * sizeof(int));
  multi = (int *)malloc(nc * sizeof(int));
  *nc_1d = *nc_m = 0;

  for (i = 0; i < nc; ++i) {
    if (df->variables[coordids[i]].nd == 1)
      oned[nc1d++] = coordids[i];
    else
      multi[ncm++] = coordids[i];
  }

  *nc_1d = nc1d;
  *cids_1d = oned;

  *nc_m = ncm;
  *cids_m = multi;
}


int search_for_coord(datafile_t *df, VariableType type, int nc,
                     int *coordids)
{
  int i;

  for (i = 0; i < nc; ++i) {
    if (df->variables[coordids[i]].type & type)
      return coordids[i];
  }

  return -1;
}


/* Check whether the two variables have the same number of
 * dimensions and that they are equal.
 */
int are_vars_same_dims(datafile_t *df, int vid1, int vid2)
{
  int i;
  df_variable_t *v1 = &df->variables[vid1];
  df_variable_t *v2 = &df->variables[vid2];

  if (v1->nd != v2->nd)
    return 0;

  for (i = 0; i < v1->nd; ++i)
    if (v1->dimids[i] != v2->dimids[i])
      return 0;

  return 1;
}


int df_map_dependent_coords(datafile_t *df, df_coord_mapping_t *cm,
                            int nc, int *coordids)
{
  int i, x_id, y_id, z_id, t_id;

  /* For the 'dependent' coordinates we must find the coordinates 'pairs'.
     We do this here in a VERY simply minded way, only the following
     coordinates are supported: XY, XYZ, XYT, XYZT. */
  x_id = search_for_coord(df, VT_X, nc, coordids);
  if (x_id < 0)
    x_id = search_for_coord(df, VT_LONGITUDE, nc, coordids);
  y_id = search_for_coord(df, VT_Y, nc, coordids);
  if (y_id < 0)
    y_id = search_for_coord(df, VT_LATITUDE, nc, coordids);
  z_id = search_for_coord(df, VT_Z, nc, coordids);
  t_id = search_for_coord(df, VT_TIME, nc, coordids);

  /* Check if XY */
  if ((x_id >= 0) && (y_id >= 0) && are_vars_same_dims(df, x_id, y_id)) {

    /* Check if XYZ */
    if ((z_id >= 0) && are_vars_same_dims(df, x_id, z_id)) {
      /* Check if XYZT */
      if ((t_id >= 0) && are_vars_same_dims(df, x_id, t_id)) {
        cm->nc = 4;
        cm->coordids = (int *)malloc(sizeof(int) * cm->nc);
        cm->coordtypes = (int *)malloc(sizeof(int) * cm->nc);
        cm->coordids[0] = x_id;
        cm->coordids[1] = y_id;
        cm->coordids[2] = z_id;
        cm->coordids[3] = t_id;
        for (i = 0; i < cm->nc; ++i)
          cm->coordtypes[i] = df->variables[cm->coordids[i]].type;
        return 1;
      } else {
        cm->nc = 3;
        cm->coordids = (int *)malloc(sizeof(int) * cm->nc);
        cm->coordtypes = (int *)malloc(sizeof(int) * cm->nc);
        cm->coordids[0] = x_id;
        cm->coordids[1] = y_id;
        cm->coordids[2] = z_id;
        for (i = 0; i < cm->nc; ++i)
          cm->coordtypes[i] = df->variables[cm->coordids[i]].type;
        return 1;
      }
    }

    /* Check of XYT */
    else if ((t_id >= 0) && are_vars_same_dims(df, x_id, t_id)) {
      cm->nc = 3;
      cm->coordids = (int *)malloc(sizeof(int) * cm->nc);
      cm->coordtypes = (int *)malloc(sizeof(int) * cm->nc);
      cm->coordids[0] = x_id;
      cm->coordids[1] = y_id;
      cm->coordids[2] = t_id;
      for (i = 0; i < cm->nc; ++i)
        cm->coordtypes[i] = df->variables[cm->coordids[i]].type;
      return 1;
    }

    /* Just XY */
    else {
      cm->nc = 2;
      cm->coordids = (int *)malloc(sizeof(int) * cm->nc);
      cm->coordtypes = (int *)malloc(sizeof(int) * cm->nc);
      cm->coordids[0] = x_id;
      cm->coordids[1] = y_id;
      for (i = 0; i < cm->nc; ++i)
        cm->coordtypes[i] = df->variables[cm->coordids[i]].type;
      return 1;
    }
  }

  return 0;
}



/* Locates the specified text attribute in the datafile_t.
 * Returns NULL if not found.
 */
char *get_text_attribute(datafile_t *df, int varid, const char *attrib)
{
  int i;

  df_variable_t *v = &df->variables[varid];

  for (i = 0; i < v->na; ++i) {
    df_attribute_t *a = &v->attributes[i];
    switch (a->type) {
    case AT_TEXT:
      if (strcasecmp(a->name, attrib) == 0)
        return (char *)a->value;
      break;

    default:
      break;
    }
  }

  return NULL;
}


/* 1D COORDINATE MAPPING */

/* Defines an XtoI and ItoX conversion. It is assumed that the
 * coordinate values are monotonic.
 *
 * NOTE: Take care to handle the dependent dimensions separately.
 */
void set_1d_coords(datafile_t *df, df_variable_t *v,
                   df_coord_mapping_t *cm)
{
  if (cm->ndepd == 0) {
    cm->free = cm_1d_free;
    cm->coords_to_indices = cm_1d_ctoi;
    cm->indices_to_coords = cm_1d_itoc;
  } else
    quit
      ("set_1d_coords: Dependent dimensions is not currently supported\n");
}



void cm_1d_free(datafile_t *df, df_coord_mapping_t *cm)
{
}


int cm_1d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
               const int depindicies[], const double coords[],
               double indices[])
{
  df_variable_t *v = &df->variables[cm->coordids[0]];
  double *data = NULL;
  int asc = 1;
  int imid;
  int ilow = 0;
  int ihigh = -1, itemp;
  double value = (v->z_is_depth) ? -coords[0] : coords[0];

  if (v->dim_as_record) {
    ihigh = df->nrecords - 1;
    data = v->data;
  } else {
    ihigh = df->dimensions[v->dimids[0]].size - 1;
    data = VAR_1D(v)[0];
  }

  /* Reverse order if descending */
  if ( (data[ihigh] - data[ilow]) < 0 ) {
    itemp = ihigh;
    ihigh = ilow;
    ilow  = itemp;
  }
  
  if (value < data[ilow]) {
    indices[0] = ilow;
    return 1;
  }
  else if (value > data[ihigh]) {
    indices[0] = ihigh;
    return 1;
  }

  /* perform binary chop to determine values either side of r */
  while (abs(ihigh - ilow) > 1) {
    imid = (ilow + ihigh) / 2;
    if (value >= data[imid])
      ilow = imid;
    else
      ihigh = imid;
  }

  if (ihigh == ilow)
    indices[0] = ilow;
  else
    if (ihigh > ilow)
      indices[0] = ilow + (value - data[ilow]) / (data[ihigh] - data[ilow]);
    else
      indices[0] = ihigh + (value - data[ihigh]) / (data[ilow] - data[ihigh]);
  
  return 1;
}


int cm_1d_itoc(datafile_t *df, df_coord_mapping_t *cm,
               const int depindicies[], const double indices[],
               double coords[])
{
  df_variable_t *v = &df->variables[cm->coordids[0]];
  double *data = NULL;
  int i = 0;
  int ihigh = -1;

  if (v->dim_as_record) {
    ihigh = df->nrecords - 1;
    data = v->data;
  } else {
    ihigh = df->dimensions[v->dimids[0]].size - 1;
    data = VAR_1D(v)[0];
  }

  i = (int)indices[0];
  if ((i < 0) || (i > ihigh))
    return 0;

  coords[0] = (indices[0] - i) * (data[i + 1] - data[i]) + data[i];
  if (v->z_is_depth)
    coords[0] = -coords[0];

  return 1;
}


/* 1D MULTI COORDINATE MAPPING */
/* Defines an X{Y{Z}}toI and ItoX{Y{Z}} conversion.
 * There are many possible ways that a 1d to multiple coordinate
 * map could be handled. The scheme choosen here will truncate all
 * fractional indicies to integers, and will find the index
 * closest to the specified coordinate when performing the inverse
 * mapping.
 */
void set_1d_multi_coords(datafile_t *df, df_variable_t *v,
                         df_coord_mapping_t *cm)
{
  if (cm->ndepd == 0) {
    cm->free = cm_1d_multi_free;
    cm->coords_to_indices = cm_1d_multi_ctoi;
    cm->indices_to_coords = cm_1d_multi_itoc;
  } else
    quit
      ("set_1d_multi_coords: Dependent dimensions is not currently supported\n");
}



void cm_1d_multi_free(datafile_t *df, df_coord_mapping_t *cm)
{
}


int cm_1d_multi_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double coords[],
                     double indices[])
{
  int dsize = df->dimensions[cm->dimids[cm->nd - 1]].size;
  double smalldist = 0.0;
  df_variable_t *vars = df->variables;
  int i, j;

  indices[0] = 0.0;
  /* Calculate weighted sum of all data points */
  for (i = 0; i < dsize; ++i) {
    double dist = 0.0;

    for (j = 0; j < cm->nc; ++j) {
      double dx = VAR_1D(&vars[cm->coordids[j]])[0][i] - coords[j];
      dist += dx * dx;
    }

    if ((i == 0) || (dist < smalldist)) {
      smalldist = dist;
      indices[0] = i;
    }
  }

  return 1;
}


int cm_1d_multi_itoc(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double indices[],
                     double coords[])
{
  int dsize = df->dimensions[cm->dimids[cm->nd - 1]].size;
  df_variable_t *vars = df->variables;
  int i, j;

  i = (int)(indices[0] + 0.5);
  i = (i < 0) ? 0 : i;
  i = (i >= 0) ? dsize - 1 : i;

  for (j = 0; j < cm->nc; ++j)
    coords[j] = VAR_1D(&vars[cm->coordids[j]])[0][i];

  return 1;
}


/* 2D COORDINATE MAPPING */

/* Defines an XYtoIJ and IJtoXY conversion.
 * There are three possible methods supported for managing XY
 * coordinates.
 *    - Rectangular.
 *    - Polar
 *    - Irregular/Numeric.
 * All of these assume a topologically rectangular piece-wise
 * linear grid (I think that's the terminology).
 */
void set_2d_coords(datafile_t *df, df_variable_t *v,
                   df_coord_mapping_t *cm)
{

  if (cm->ndepd == 0) {
    df_variable_t *coordvar0 = &df->variables[cm->coordids[0]];
    df_variable_t *coordvar1 = &df->variables[cm->coordids[1]];

    /* Do each of the coordinate contain an 'analytic' attribute. If so,
       are they the same ?. */
    char *xa = get_text_attribute(df, cm->coordids[0], "analytic");
    char *ya = get_text_attribute(df, cm->coordids[1], "analytic");
    df_attribute_t *ta = NULL;

    if (xa && ya && (strcasecmp(xa, ya) == 0)) {
      char atype[MAXLINELEN];
      if (sscanf(xa, "%s", atype) == 1) {

        /* An analytic solution exists. Is this rectangular or polar. */
        if (strcasecmp(atype, "rectangular") == 0) {
          gi_rect_t *rgi = (gi_rect_t *)malloc(sizeof(gi_rect_t));
          if (sscanf(xa, "%s %lf %lf %d %d %lf %lf %lf %lf %lf",
                     atype, &rgi->ioffset, &rgi->joffset,
                     &rgi->ni, &rgi->nj,
                     &rgi->x0, &rgi->y0, &rgi->dx,
                     &rgi->dy, &rgi->th) == 10) {
            rgi->sinth = sin(rgi->th * M_PI / 180.0);
            rgi->costh = cos(rgi->th * M_PI / 180.0);

            /* Shift origin to account for i/j offset */
            rgi->x0 = rgi->ioffset * rgi->dx * rgi->costh
              - rgi->joffset * rgi->dy * rgi->sinth + rgi->x0;
            rgi->y0 = rgi->ioffset * rgi->dx * rgi->sinth
              + rgi->joffset * rgi->dy * rgi->costh + rgi->y0;
            cm->free = cm_rect_2d_free;
            cm->coords_to_indices = cm_rect_2d_ctoi;
            cm->indices_to_coords = cm_rect_2d_itoc;
            cm->special_data = rgi;
            set_geo_transform(df, cm);
            return;
          }
        }

        else if (strcasecmp(atype, "polar") == 0) {
          gi_polar_t *pgi = (gi_polar_t *)malloc(sizeof(gi_polar_t));
          if (sscanf(xa, "%s %lf %lf %d %d %lf %lf %lf %lf %lf",
                     atype, &pgi->ioffset, &pgi->joffset,
                     &pgi->ni, &pgi->nj,
                     &pgi->x0, &pgi->y0, &pgi->arc,
                     &pgi->rmin, &pgi->rotation) == 10) {

            pgi->rotation = pgi->rotation * DEG2RAD;
            if (pgi->rotation >= 2.0 * M_PI)
              pgi->rotation -= 2.0 * M_PI;
            if (pgi->rotation < 0.0)
              pgi->rotation += 2.0 * M_PI;
            pgi->dth = (pgi->arc / (pgi->ni - 1)) * DEG2RAD;
            pgi->logrmin = log(pgi->rmin);

            cm->free = cm_polar_2d_free;
            cm->coords_to_indices = cm_polar_2d_ctoi;
            cm->indices_to_coords = cm_polar_2d_itoc;
            cm->special_data = pgi;
            set_geo_transform(df, cm);
            return;
          }
        }

        else
          warn
            ("set_2d_coords: '%s' has unknown (%s) analytic type, using bilinear 2d coordinate mapping.\n",
             v->name, atype);
      } else
        warn
          ("set_2d_coords: '%s' has malformed analytic attribute, using bilinear 2d coordinate mapping.\n",
           v->name);
    }

    /* 
     * Use generic bilinear mapping if the dimensions of each coordinate
     * variable are the same.
     */
    if ((coordvar0->nd == 2) && (coordvar1->nd == 2) &&
        (df->dimensions[coordvar0->dimids[0]].size ==
         df->dimensions[coordvar1->dimids[0]].size) &&
        (df->dimensions[coordvar0->dimids[1]].size ==
         df->dimensions[coordvar1->dimids[1]].size)) {
      double **gx, **gy;
      int nx, ny;
      /* Set up transformation function pointer */
      cm->coords_to_indices = cm_bilinear_2d_ctoi;
      cm->indices_to_coords = cm_bilinear_2d_itoc;

      /*
       * 29/11/11 FR : Special handling for SHOC gridded files
       *
       * The original code (in the else section)
       * initialises grid_xytoij with the coordinate variables of the
       * data variable - this is all you can do in general. If this is
       * a simple (eg. rectangular) grid and its domain wholly
       * encompasses the model grid then, there is no problem.
       * The problem occurs when you have very thin branches - the
       * xytoij mechanism gets "stuck" when there is a river only one
       * cell wide, for example. The boundary algorithm gets to the end
       * of the river and has nowhere to go and you get the
       * "create_leaf, boundary vertex not found" error which is
       * a misleading message anyway. What I've tried to do below is
       * to always initialise xytoij with the grid corners in the case
       * of cell centred variables. Care must be taken *not* to 
       * have NaN in land cells though - use zero_fill or cascade_search in
       * the file to interpolate from
       */
      // 11/06/2019 : This test isn't valid anymore not that title has swtiched
      /*
	if ( (ta = df_get_global_attribute(df, "title")) != NULL &&
	strncmp(ATT_TEXT(ta), "SHOC", 4) == 0  && This is a SHOC file
	strcmp(coordvar0->name, "x_centre") == 0 && and we're dealing
	strcmp(coordvar1->name, "y_centre") == 0 )  with cell centred
      */
      /* Infer SHOC standard file */
      if ( strcmp(coordvar0->name, "x_centre") == 0 &&
	   strcmp(coordvar1->name, "y_centre") == 0 &&
	   df_get_index(df, "x_grid") > -1 &&
	   df_get_index(df, "y_grid") > -1 )
	{
	  /* 
	   * Always default to the grid corners
	   */
	  df_variable_t *vid_x = df_get_variable_by_name(df, "x_grid");
	  df_variable_t *vid_y = df_get_variable_by_name(df, "y_grid");
	  
	  df_read_record(df, vid_x, 0);
	  df_read_record(df, vid_y, 0);
	  
	  gx = VAR_2D(vid_x)[0];
	  gy = VAR_2D(vid_y)[0];
	  nx = df->dimensions[vid_x->dimids[1]].size - 1;
	  ny = df->dimensions[vid_x->dimids[0]].size - 1;

	} else {
	/* Assume regular grid. This may fail for complex curvilinear grids */
	df_read_record(df, coordvar0, 0);
	df_read_record(df, coordvar1, 0);
	
	gx = VAR_2D(coordvar0)[0];
	gy = VAR_2D(coordvar1)[0];
	nx = df->dimensions[coordvar0->dimids[1]].size - 1;
	ny = df->dimensions[coordvar0->dimids[0]].size - 1;
      }
      
      /* Initialise the xytoij mapping */
      if (!df_is_irule(df))
	cm->special_data = grid_xytoij_init_hash(gx,gy,nx,ny,df_default_hashtable_size);
      cm->free = cm_bilinear_2d_free;
      set_geo_transform(df, cm);
    } else
      quit
        ("set_2d_coords: Coordinate variables must have the same number and length of dimensions.\n");
  } else
    quit
      ("set_2d_coords: 2d mapping with dependent dimensions is not currently supported.\n");
}


/* GEOGRAPHIC TRANSFORMATIONS */

/* Infer from the 'user' geotype and any coordinate or projection attributes
 * within the file whether a geo transform is required.
 */
void set_geo_transform(datafile_t *df, df_coord_mapping_t *cm)
{
  df_variable_t *cvar0 = &df->variables[cm->coordids[0]];
  df_variable_t *cvar1 = &df->variables[cm->coordids[1]];

  if (cvar0->type == cvar1->type)
    return;

  /* Check if the coordinate variables are X and Y. If so, then search for
     the 'projection' attribute */
  if ((cvar0->type & (VT_X | VT_Y)) && (cvar1->type & (VT_X | VT_Y))) {
    char *xproj = get_text_attribute(df, cm->coordids[0], "projection");
    char *yproj = get_text_attribute(df, cm->coordids[1], "projection");
    if ((xproj != NULL) && (yproj != NULL) &&
        (strcasecmp(xproj, yproj) == 0)) {
      if (strcasecmp(xproj, GEOGRAPHIC_TAG) == 0)
        cm->transform = create_geo_transform(GT_GEOGRAPHIC, NULL,
                                             df->geotype, df->projection);
      else
        cm->transform = create_geo_transform(GT_PROJECTION, xproj,
                                             df->geotype, df->projection);
    }
  }

  /* Check if the coordinate variables are lat and long. */
  else if ((cvar0->type & (VT_LATITUDE | VT_LONGITUDE))
           && (cvar1->type & (VT_LATITUDE | VT_LONGITUDE))) {
    cm->transform = create_geo_transform(GT_GEOGRAPHIC, NULL,
                                         df->geotype, df->projection);
  }
}


/* Build a transformation structure for geographic projections.
 */
df_coord_transform_t *create_geo_transform(int int_geotype,
                                           char *int_proj, int ext_geotype,
                                           char *ext_proj)
{
  df_coord_transform_t *trans = NULL;
  df_geo_transform_t *data = NULL;
  char *args[256];

  if ((int_geotype == GT_NONE) || (ext_geotype == GT_NONE))
    return NULL;

  if ((int_geotype == GT_GEOGRAPHIC) && (ext_geotype == GT_GEOGRAPHIC))
    return NULL;

  if ((int_geotype == GT_PROJECTION) && (ext_geotype == GT_PROJECTION)
      && (strcasecmp(int_proj, ext_proj) == 0))
    return NULL;


  /* Build the geo transformation data */
  data = (df_geo_transform_t *)malloc(sizeof(df_geo_transform_t));
  memset(data, 0, sizeof(df_geo_transform_t));

  data->int_geotype = int_geotype;
  if (int_geotype == GT_PROJECTION) {
    int nargs = parseline(strdup(int_proj), args, 256);
    data->int_mp = mp_init(nargs, args);
  }

  data->ext_geotype = ext_geotype;
  if (ext_geotype == GT_PROJECTION) {
    int nargs = parseline(strdup(ext_proj), args, 256);
    data->ext_mp = mp_init(nargs, args);
  }

  /* Build the transform */
  trans = (df_coord_transform_t *)malloc(sizeof(df_coord_transform_t));
  memset(trans, 0, sizeof(df_coord_transform_t));

  trans->private_data = data;
  trans->forward = proj_forward;  /* User to File */
  trans->inverse = proj_inverse;  /* File to User */
  trans->free = proj_free;

  return trans;
}

/* Geographic transformations.
 */
void proj_forward(void *d, int n, double from[], double to[])
{
  df_geo_transform_t *data = (df_geo_transform_t *)d;
  double lat, lon;

  /* User to lat/long */
  switch (data->ext_geotype) {
  case GT_GEOGRAPHIC:
    lon = from[0];
    lat = from[1];
    break;

  case GT_PROJECTION:
    mp_inverse(data->ext_mp, from[0], from[1], &lat, &lon);
    break;
  }

  /* lat/long to file. */
  switch (data->int_geotype) {
  case GT_GEOGRAPHIC:
    to[0] = lon;
    to[1] = lat;
    break;

  case GT_PROJECTION:
    mp_forward(data->int_mp, lat, lon, &to[0], &to[1]);
    break;
  }
}

void proj_inverse(void *d, int n, double from[], double to[])
{
  df_geo_transform_t *data = (df_geo_transform_t *)d;
  double lat, lon;

  /* File to lat/long */
  switch (data->int_geotype) {
  case GT_GEOGRAPHIC:
    lon = from[0];
    lat = from[1];
    break;

  case GT_PROJECTION:
    mp_inverse(data->int_mp, from[0], from[1], &lat, &lon);
    break;
  }

  /* lat/long to User. */
  switch (data->ext_geotype) {
  case GT_GEOGRAPHIC:
    to[0] = lon;
    to[1] = lat;
    break;

  case GT_PROJECTION:
    mp_forward(data->ext_mp, lat, lon, &to[0], &to[1]);
    break;
  }
}

void proj_free(void *d)
{
  if (d)
    free(d);
}


/* Rectangular mapping */
void cm_rect_2d_free(datafile_t *df, df_coord_mapping_t *cm)
{
  free(cm->special_data);
}


int cm_rect_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                    const int depindicies[], const double coords[],
                    double indices[])
{
  gi_rect_t *rgi = (gi_rect_t *)cm->special_data;
  double i, j, x, y;

  x = coords[0];
  y = coords[1];

  i = ((x - rgi->x0) * rgi->costh + (y - rgi->y0) * rgi->sinth) / rgi->dx;
  j = ((y - rgi->y0) * rgi->costh - (x - rgi->x0) * rgi->sinth) / rgi->dy;

  indices[0] = j;
  indices[1] = i;

  return 1;
}

int cm_rect_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                    const int depindicies[], const double indices[],
                    double coords[])
{
  gi_rect_t *rgi = (gi_rect_t *)cm->special_data;
  double i, j, x, y;

  i = indices[1];
  j = indices[0];

  x = i * rgi->dx * rgi->costh - j * rgi->dy * rgi->sinth + rgi->x0;
  y = i * rgi->dx * rgi->sinth + j * rgi->dy * rgi->costh + rgi->y0;

  coords[0] = x;
  coords[1] = y;

  return 1;
}


/* Square cell polar grid */
void cm_polar_2d_free(datafile_t *df, df_coord_mapping_t *cm)
{
  free(cm->special_data);
}


int cm_polar_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double coords[],
                     double indices[])
{
  gi_polar_t *pgi = (gi_polar_t *)cm->special_data;
  double i, j, x, y;
  double r2, th, xdiff, ydiff;

  x = coords[0];
  y = coords[1];

  xdiff = (x - pgi->x0);
  ydiff = (y - pgi->y0);
  r2 = xdiff * xdiff + ydiff * ydiff;

  th = atan2(ydiff, xdiff);
  while (th > pgi->rotation)
    th -= 2.0 * M_PI;
  i = (pgi->rotation - th) / pgi->dth;
  j = (0.5 * log(r2) - pgi->logrmin) / pgi->dth;

  i -= pgi->ioffset;
  j -= pgi->joffset;

  indices[0] = j;
  indices[1] = i;

  return 1;
}

int cm_polar_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                     const int depindicies[], const double indices[],
                     double coords[])
{
  gi_polar_t *pgi = (gi_polar_t *)cm->special_data;
  double i, j, x, y, r, th;

  i = indices[1] + pgi->ioffset;
  j = indices[0] + pgi->joffset;

  r = pgi->rmin * exp(pgi->dth * j);
  th = pgi->rotation - i * pgi->dth;
  x = pgi->x0 + r * cos(th);
  y = pgi->y0 + r * sin(th);

  coords[0] = x;
  coords[1] = y;

  return 1;
}


/* Generic time-invariant bilinear mapping.  */
void cm_bilinear_2d_free(datafile_t *df, df_coord_mapping_t *cm)
{
  /*UR-CHANGED free the whole data */
  if (cm->special_data != NULL)
    tree_destroy((xytoij_tree_t *)cm->special_data);
  cm->special_data = NULL;
}

int cm_bilinear_2d_ctoi(datafile_t *df, df_coord_mapping_t *cm,
                        const int depindicies[], const double coords[],
                        double indices[])
{
  xytoij_tree_t *partition = (xytoij_tree_t *)cm->special_data;

  return grid_xytofij(partition, coords[0], coords[1], &indices[1],
                      &indices[0]);
}

int cm_bilinear_2d_itoc(datafile_t *df, df_coord_mapping_t *cm,
                        const int depindicies[], const double indices[],
                        double coords[])
{
  xytoij_tree_t *partition = (xytoij_tree_t *)cm->special_data;

  return grid_fgrid_ijtoxy(partition, indices[1], indices[0], &coords[0],
                           &coords[1]);
}
