/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/tf_nc_cmr.c
 *
 *  \brief Functions to access CMR style netCDF bathy/topo files
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: tf_nc_cmr.c 5833 2018-06-27 00:21:35Z riz008 $
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "ems.h"

/* external prototypes */
/* - topo */
int get_lat_index(topo_details_t *td, double lat);

int get_lon_index(topo_details_t *td, double lon);

/* local functions */
/**
 * check to see if file is a valid nc file 
 */
int tf_nc_cmr_is_valid(char *fname) {
   int is_valid = 0;
   int ncid   = -1;
   int latVid = -1;
   int lonVid = -1;
   int htVid  = -1;

   /* check to see if the file exists */
   if(nc_open(fname, NC_NOWRITE, &ncid) == NC_NOERR) {
      /* check to see if it contain latitude, longitude, and height */
      if((nc_inq_varid(ncid, "latitude", &latVid) == NC_NOERR)
         && (nc_inq_varid(ncid, "longitude", &lonVid) == NC_NOERR)
         && ((nc_inq_varid(ncid, "height", &htVid) == NC_NOERR)
         || (nc_inq_varid(ncid, "depth", &htVid) == NC_NOERR))) {
         /* file is a valid CMR nc_topo file */
         is_valid = 1;
         nc_close(ncid);
      }
   }

   return(is_valid);
} /* end tf_nc_cmr_isvalid */

topo_details_t *tf_nc_cmr_open(char *fname) {
   topo_details_t *td;
   int i;
   int ncid = -1;
   int latid = -1;
   size_t nlats = 0;
   int latvid = -1;
   int lonid = -1;
   size_t nlons = 0;
   int lonvid = -1;
   size_t start[2];
   size_t end[2];
   double tmp_del;

  td = (topo_details_t*)malloc(sizeof(topo_details_t));
  memset(td, 0, sizeof(topo_details_t));
  td->extent = (topo_extent_t*)malloc(sizeof(topo_extent_t));
  memset(td->extent, 0, sizeof(topo_extent_t));

   strcpy(td->name, fname);
   /* already checked if it was openable */
   nc_open(fname, NC_NOWRITE, &ncid);
   td->ncid = ncid;
   td->topo_type = GRID_POINT;

   nc_inq_dimid(ncid, "nlines", &latid);
   nc_inq_dimlen(ncid, latid, &nlats);
   td->nlats = nlats;
   td->lats = d_alloc_1d(nlats);
   nc_inq_varid(ncid, "latitude", &latvid);
   start[0] = 0;
   end[0] = nlats;
   nc_get_vara_double(ncid, latvid, start, end, td->lats);

   /* assume that the file is rectilinear */
   td->dely = fabs(td->lats[1] - td->lats[0]);
   /* check to see if it is uniform */
   for(i=2;i<nlats;++i) {
      tmp_del = fabs(td->lats[i] - td->lats[i-1]);
      /* cast to float for proper compare */
      if((float)tmp_del != (float)td->dely) {
         td->dely = -1;
         break;
      }
   }

   nc_inq_dimid(ncid, "ncells", &lonid);
   nc_inq_dimlen(ncid, lonid, &nlons);
   td->nlons = nlons;
   td->lons = d_alloc_1d(nlons);
   nc_inq_varid(ncid, "longitude", &lonvid);
   start[0] = 0;
   end[0] = nlons;
   nc_get_vara_double(ncid, lonvid, start, end, td->lons);

   /* assume that the file is rectilinear */
   td->delx = fabs(td->lons[1] - td->lons[0]);
   /* check to see if it is uniform */
   for(i=2;i<nlons;++i) {
      tmp_del = fabs(td->lons[i] - td->lons[i-1]);
      /* cast to float for proper compare */
      if((float)tmp_del != (float)td->delx) {
         td->delx = -1;
         break;
      }
   }

   if(td->delx < 0 || td->dely < 0) {
      quit("tf_nc_cmr_open: inconsistent delxs and delys are not supported at this time");
   }

   /* set the unrelated values to "non values" */
   strcpy(td->basedir, "");
   td->stats = NULL;
   td->num_file_tiles = 0;
   td->num_file_rows = 0;
   td->num_file_cols = 0;

   return(td);
} /* end tf_nc_cmr_open */

int tf_nc_cmr_close(topo_t *tf) {
   int close_ok = 0;

   if(nc_close(tf->td->ncid) == NC_NOERR) {
      close_ok = 1;
   }

   return(close_ok);
} /* tf_nc_cmr_close */

/*
 * get the z values from the nc topo file
 */
int tf_nc_cmr_getz(topo_details_t *td, double *x, double *y, double *z, 
   int nx, int ny) {
   int got_ok = 0;
   int i, j, k;
   int htvid = -1;
   double make_topo = -1.0;
   size_t start[2];
   size_t count[2];
   double **zz;
   int err;
   int latindex, lonindex;

  if(nc_inq_varid(td->ncid, "height", &htvid) == NC_NOERR)
    make_topo = 1.0;
  else {
    /* make bathy into topo */
    nc_inq_varid(td->ncid, "depth", &htvid);
    make_topo = -1.0;
  }

   /* else there is an nx X ny matrix to read
    * - find the start latindex
    * note: that the X and Ys have already been computed "above" so that nx
    * and ny would be known. */
   latindex = get_lat_index(td, y[0]);
   lonindex = get_lon_index(td, x[0]);

  if(latindex < 0) {
    warn("tf_nc_cmr_getz: lat index not found %.2f\n", y[0]);
    return(0);
  }
  if(lonindex < 0) {
    warn("tf_nc_cmr_getz: lon index not found %.2f\n", x[0]);
    return(0);
  }
  start[0] = (size_t)latindex;
  start[1] = (size_t)lonindex;

   /* set the lon count */
   count[0] = nx;
   /* set the lat count */
   count[1] = ny;
   /* allocate the zz array */
   zz = d_alloc_2d(nx,ny);
   err = nc_get_vara_double(td->ncid, htvid, start, count, zz[0]);
   if(err != NC_NOERR) {
      got_ok = 0;
   } else {
      k = 0;
      for(i=0;i<ny;++i) {
         for(j=0;j<nx;++j) {
           /* NOTE: zz[y][x] */
            z[k] = zz[i][j] * make_topo;
            ++k;
         }
      }
      got_ok = 1;
   } /* end if err != NC_NOERR */
   d_free_2d(zz);

   return(got_ok);
} /* end tf_nc_cmr_getz */

topo_extent_t *tf_nc_cmr_extent(topo_details_t *td) {
   topo_extent_t *te;

  te = (topo_extent_t*)malloc(sizeof(topo_extent_t));
  memset(te, 0, sizeof(topo_extent_t));

  if(td->topo_type == GRID_POINT) {
    te->minLat = td->lats[0] - (td->dely / 2.0);
    te->maxLat = td->lats[td->nlats - 1] + (td->dely / 2.0);
    te->minLon = td->lons[0] - (td->delx / 2.0);
    te->maxLon = td->lons[td->nlons - 1] + (td->delx / 2.0);
  } else if(td->topo_type == GRID_POLY) {
    te->minLat = td->lats[0];
    te->maxLat = td->lats[td->nlats - 1];
    te->minLon = td->lons[0];
    te->maxLon = td->lons[td->nlons - 1];
  } else
    quit("tf_nc_cmr_extent: unsupported topo_type\n");

   return(te);
} /* end tf_nc_cmr_extent */

topo_tile_t *tf_nc_cmr_get_tile(topo_t *tf, double xcoord,  double ycoord) {
   topo_tile_t *tile;
   int lat_index, lon_index;
   int i;
   double darray1[1] = {0.0};
   double darray2[1] = {0.0};

  tile = (topo_tile_t*)malloc(sizeof(topo_tile_t));
  memset(tile, 0, sizeof(topo_tile_t));
  tile->td = (topo_details_t*)malloc(sizeof(topo_details_t));
  memset(tile->td, 0, sizeof(topo_details_t)); 
  tile->td->extent = (topo_extent_t*)malloc(sizeof(topo_extent_t));
  memset(tile->td->extent, 0, sizeof(topo_extent_t));
  tile->neighbours = (topo_neighbours_t*)malloc(sizeof(topo_neighbours_t));
  memset(tile->neighbours, 0, sizeof(topo_neighbours_t)); 

  /* tile inherits topo_type from parent file */
  tile->td->topo_type = tf->td->topo_type;

   /* set the flag to indicate it is a nontiled nc file */
   tile->td->ncid = -1;

   /* set the lats */
   lat_index = get_lat_index(tf->td, ycoord);
   /* make sure we do not go beyond file boundary */
   tile->td->nlats = tf->td->nlats - lat_index;
   /* make sure we do not go beyond the max tile height */
   if(tile->td->nlats > MAXTILEHEIGHT) {
      tile->td->nlats = MAXTILEHEIGHT;
   }
   tile->td->lats = d_alloc_1d(tile->td->nlats);
   for(i=0;i<tile->td->nlats;++i) {
      tile->td->lats[i] = tf->td->lats[lat_index + i];
   }
   tile->td->dely = fabs(tile->td->lats[1] - tile->td->lats[0]);
  if(tile->td->topo_type == GRID_POLY) {
    tile->td->extent->minLat = tile->td->lats[0];
    tile->td->extent->maxLat = tile->td->lats[tile->td->nlats - 1];
  } else if(tile->td->topo_type == GRID_POINT) {
    tile->td->extent->minLat = tile->td->lats[0] - (tile->td->dely / 2.0);
    tile->td->extent->maxLat = tile->td->lats[tile->td->nlats - 1]
      + (tile->td->dely / 2.0);
  } else 
    quit("tf_nc_cmr_get_tile: unsupported topo_type\n");

   /* set the lons */
   lon_index = get_lon_index(tf->td, xcoord);
   /* make sure we do not go beyond file boundary */
   tile->td->nlons = tf->td->nlons - lon_index;
   /* make sure we do not go beyond the max tile height */
   if(tile->td->nlons > MAXTILEWIDTH) {
      tile->td->nlons = MAXTILEWIDTH;
   }
   tile->td->lons = d_alloc_1d(tile->td->nlons);
   for(i=0;i<tile->td->nlons;++i) {
      tile->td->lons[i] = tf->td->lons[lon_index + i];
   }
   tile->td->delx = fabs(tile->td->lons[1] - tile->td->lons[0]);
  if(tile->td->topo_type == GRID_POLY) {
    tile->td->extent->minLon = tile->td->lons[0];
    tile->td->extent->maxLon = tile->td->lons[tile->td->nlons - 1];
  } else if(tile->td->topo_type == GRID_POINT) {
    tile->td->extent->minLon = tile->td->lons[0] - (tile->td->delx / 2.0);
    tile->td->extent->maxLon = tile->td->lons[tile->td->nlons - 1]
      + (tile->td->delx / 2.0);
  } else 
    quit("tf_nc_cmr_get_tile: unsupported topo_type\n");

   /* at this stage the tile is populatable
    * allocate space for the Z data */
   tile->td->z = d_alloc_1d(tile->td->nlats * tile->td->nlons);
   darray1[0] = xcoord;
   darray2[0] = ycoord;

   tf->getz(tf->td, darray1, darray2, tile->td->z, tile->td->nlons, 
      tile->td->nlats);

   /* set the neighbours to unknown */
   tile->neighbours->topLeft  = -1;
   tile->neighbours->top      = -1;
   tile->neighbours->topRight = -1;
   tile->neighbours->bot      = -1;
   tile->neighbours->botLeft  = -1;
   tile->neighbours->botRight = -1;
   tile->neighbours->left     = -1;
   tile->neighbours->right    = -1;

   return(tile);
} /* end tf_nc_cmr_get_tile */
