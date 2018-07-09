/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/coast.c
 *
 *  \brief Frontend for reading coastline files
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: coast.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ems.h"

/* External references */
extern int tiled_coast_isvalid(const char *filename);
extern coast_file_t* tiled_coast_open(const char *filename);

extern int poly_coast_isvalid(const char *filename);
extern coast_file_t* poly_coast_open(const char *filename);


/** Open the coastline file for reading. The file shall be queried first
  * to determine whether it is one of the supported coastline file formats.
  *
  * @param filename The coastline file to open.
  * @return The file details. NULL if the file format is not supported.
  */
coast_file_t* coast_open(const char *filename) {

   if (tiled_coast_isvalid(filename))
      return tiled_coast_open(filename);
   else if (poly_coast_isvalid(filename))
      return poly_coast_open(filename);

   return NULL;
}

/** Close the coastline file.
  *
  * @param cf Pointer to coastline file details.
  */
void coast_close(coast_file_t *cf) {
   assert(cf != NULL);
   cf->close(cf);
}

/** Get the set of 'coastline tiles' that fall within the specified
  * geographic region specified.
  *
  * The tiles can be read one after the other using an iterator pattern.
  *
  * coast_tile_t *tile;
  * coast_tile_t *ti = coast_get_tile_iterator(cf, ...);
  *
  * while ((tile = ti->next_tile(ti)) != NULL) {
  *    ...
  *    ti->free_tile(tile);
  * }
  * coast_free_tile_iterator(ti);
  *
  * @param cf Pointer to coastline file details.
  * @param minlon Miminum longitude of region of interest.
  * @param maxlon Miminum longitude of region of interest.
  * @param minlat Miminum latitude of region of interest.
  * @param maxlat Miminum latitude of region of interest.
  * @return Pointer to the tile iterator.
  */
coast_tile_iterator_t *coast_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat) {
   
   assert(cf != NULL);
   return cf->get_tile_iterator(cf, minlon, maxlon, minlat, maxlat);
}

/** Release memory associated with the tile iterator.
  *
  * @param ti Pointer to the tile iterator.
  */
void coast_free_tile_iterator(coast_tile_iterator_t *ti) {
   assert(ti != NULL);
   return ti->file->free_tile_iterator(ti);
    
}

/** For a given set of geographic locations, determine whether the point is
  * on land or water.
  *
  * @param cf Pointer to coastline file details.
  * @param lon Array of longitudes.
  * @param lat Array of latitudes.
  * @param n Number of points.
  * @return Array of boolean land values (non-zero is land, zero is not).
  */
int *coast_is_land(coast_file_t *cf, double *lon, double *lat, int n) {

   double minlon = 360, maxlon = -360;
   double minlat = 90, maxlat = -90;
   int i=0, j=0;
   int *is_land = NULL;
   coast_tile_iterator_t *ti;
   coast_tile_t *tile;

   if (n <= 0)
       return NULL;

   /* Find the min/max latitude and longitude for the set of points */
   for (i=0; i<n; ++i) {
      if (lon[i] < minlon) minlon = lon[i];
      if (lon[i] > maxlon) maxlon = lon[i];
      if (lat[i] < minlat) minlat = lat[i];
      if (lat[i] > maxlat) maxlat = lat[i];
   }

   is_land = (int *)malloc(sizeof(int)*n);
   memset(is_land, 0, sizeof(int)*n);

   /* Load a tile at a time, and compare the lat/lons */
   ti = coast_get_tile_iterator(cf, minlon, maxlon, minlat, maxlat);
   while((tile = ti->next_tile(ti)) != NULL) {
      for (i=0; i<tile->npolys; ++i) {
         poly_t *p = tile->polys[i];
         if (tile->level[i] != 1) continue;
         for (j=0; j<n; ++j) {
            if (is_land[j]) continue;
            if (poly_contains_point(p, lon[j], lat[j]))
               is_land[j] = 1;
         }
      }
      
      ti->free_tile(tile);
   }
   coast_free_tile_iterator(ti);

   return is_land;
}
