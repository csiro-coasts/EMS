/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/topo_factory.c
 *
 *  \brief Function to associate bathy/topo file readers with topo API
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: topo_factory.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ems.h"

/* NetCDF CMR bathyfile prototypes. */
int tf_nc_cmr_is_valid(char *fname);
topo_details_t *tf_nc_cmr_open(char *fname);
topo_extent_t *tf_nc_cmr_extent(topo_details_t *td);
int tf_nc_cmr_getz(topo_details_t *td,
      double x[], double y[], double z[],
      int nx, int ny);
int tf_nc_cmr_close(topo_t *tf);

/* NetCDF CMR tiled topofile prototypes. */
int tf_nc_cmr_tiled_is_valid(char *fname);
topo_details_t *tf_nc_cmr_tiled_open(char *fname);
topo_extent_t *tf_nc_cmr_tiled_extent(topo_details_t *td);
int tf_nc_cmr_tiled_getz(topo_details_t *td,
      double x[], double y[], double z[],
      int nx, int ny);
int tf_nc_cmr_tiled_close(topo_t *tf);

/**
 * NOTE: only basic NetCDF topo files supported at present
 * returns a pointer to the topo_file structure
 * @param name the name of the file to use for topo data
 * @return a pointer to a topo_t structure for the file or NULL
 */
topo_t *topo_get_file_from_factory(char *name) {
   topo_t *tf = NULL;

   if(tf_nc_cmr_is_valid(name)) {
      /* OK, we have a legitamite CMR nc_topo file */
      tf = (topo_t*)malloc(sizeof(topo_t));
      memset(tf,0,sizeof(topo_t));
      tf->open = tf_nc_cmr_open;
      tf->close = tf_nc_cmr_close;
      tf->getz = tf_nc_cmr_getz;
      tf->get_extent = tf_nc_cmr_extent;

   } else if(tf_nc_cmr_tiled_is_valid(name)) {
      /* OK, John A has a tiled version of these so... */
      tf = (topo_t*)malloc(sizeof(topo_t));
      memset(tf,0,sizeof(topo_t));
      tf->open = tf_nc_cmr_tiled_open;
      tf->close = tf_nc_cmr_tiled_close;
      tf->getz = tf_nc_cmr_tiled_getz;
      tf->get_extent = tf_nc_cmr_tiled_extent;
   } else {
      quit("topo_get_file_from_factory: don't know how to handle file format for \"%s\"\n", name);
   }

   return(tf);
}
