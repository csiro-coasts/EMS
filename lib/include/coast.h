/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/coast.h
 *
 *  \brief Header file for reading coastline files
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: coast.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _COAST_H
#define _COAST_H

/* structures and typedefs */
typedef struct coast_file coast_file_t;
typedef struct coast_tile coast_tile_t;
typedef struct coast_tile_iterator coast_tile_iterator_t;

struct coast_file {
   char filename[MAXFNAMELEN];

   void (*close)(coast_file_t *cf);
   coast_tile_iterator_t *(*get_tile_iterator)(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat);
   void (*free_tile_iterator)(coast_tile_iterator_t *ti);
   void *private_data;
};


struct coast_tile {
   int npolys;
   poly_t **polys;
   int *level;
};

struct coast_tile_iterator {
   coast_file_t *file;
   int start_lon_index;
   int start_lat_index;
   int end_lon_index;
   int end_lat_index;
   int allow_wrap;

   int cur_lon_index;
   int cur_lat_index;

   coast_tile_t *(*next_tile)(coast_tile_iterator_t *it);
   void (*free_tile)(coast_tile_t *tile);

   void *private_data;
};

/* Prototypes */
coast_file_t* coast_open(const char *filename);
void coast_close(coast_file_t *cf);
coast_tile_iterator_t *coast_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat);
void coast_free_tile_iterator(coast_tile_iterator_t *ti);
int *coast_is_land(coast_file_t *cf, double *lon, double *lat, int n);

#endif
