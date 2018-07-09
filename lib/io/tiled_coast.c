/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/tiled_coast.c
 *
 *  \brief A coastline reader for the CMR tiled coastline file format
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: tiled_coast.c 5838 2018-06-27 03:33:01Z riz008 $
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ems.h"

/* definitions */
#define FLIP 1
#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) (((i4) >> 24) + (((i4) >> 8) & 65280) + (((i4) & 65280) << 8) + (((i4) & 255) << 24))

/* structures */
typedef struct tiled_file_data {
   char version[5];
   double minlon;
   double maxlon;
   double minlat;
   double maxlat;
   double dx;
   double dy;
   int nlons;
   int nlats;
   long** indicies;
} tiled_file_data_t;

/* public prototypes */
int tiled_coast_isvalid(const char *filename);
coast_file_t* tiled_coast_open(const char *filename);

/* private prototypes */
static void tc_close(coast_file_t *cf);
static coast_tile_iterator_t *tc_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat);
static void tc_free_tile_iterator(coast_tile_iterator_t *ti);
static coast_tile_t *tc_next_tile(coast_tile_iterator_t *ti);
static void tc_free_tile(coast_tile_t *tile);
static int read_int(FILE* fp);

int tiled_coast_isvalid(const char *filename) {
   char tag[5];
   int valid = 0;

   FILE *fp = fopen(filename, "rb");
   if (fp != NULL) {
      memset(tag, 0, sizeof(tag));
      if (fread(tag, 1, 4, fp) == 4)
         valid = (strcmp(tag, "TILE") == 0);
      fclose(fp);
   }

   return valid;
}

coast_file_t* tiled_coast_open(const char *filename) {

   coast_file_t* cf = NULL;
   char tag[5];

   FILE *fp = fopen(filename, "rb");
   if (fp != NULL) {
      memset(tag, 0, sizeof(tag));
      if (fread(tag, 1, 4, fp) == 4) {
         if (strcmp(tag, "TILE") == 0) {
            tiled_file_data_t *cfp;
            int i=0, j=0;
            
            cf = (coast_file_t *)malloc(sizeof(coast_file_t));
            memset(cf, 0, sizeof(coast_file_t));
            strcpy(cf->filename, filename);
            cf->close = tc_close;
            cf->get_tile_iterator = tc_get_tile_iterator;
            cf->free_tile_iterator = tc_free_tile_iterator;

            cfp = (tiled_file_data_t *)malloc(sizeof(tiled_file_data_t));
            memset(cfp, 0, sizeof(tiled_file_data_t));
            cf->private_data = cfp;

            if(!fread(cfp->version, 1, 4, fp))
	      quit("tiled_coast_open: incorrect format, version not found in %s\n",
		   filename);
            cfp->minlon = read_int(fp) * 1.0e-6;
            cfp->minlat = read_int(fp) * 1.0e-6;
            cfp->maxlon = read_int(fp) * 1.0e-6;
            cfp->maxlat = read_int(fp) * 1.0e-6;
            cfp->nlons = read_int(fp);
            cfp->nlats = read_int(fp);
            cfp->dx = (cfp->maxlon - cfp->minlon) / cfp->nlons;
            cfp->dy = (cfp->maxlat - cfp->minlat) / cfp->nlats;
            cfp->indicies = l_alloc_2d(cfp->nlons, cfp->nlats);
            for (i=0; i<cfp->nlons; ++i)
               for (j=0; j<cfp->nlats; ++j)
                   cfp->indicies[j][i] = (long)read_int(fp);
         }
      }
      fclose(fp);
   }

   return cf;
}

static void tc_close(coast_file_t *cf) {
   if (cf != NULL) {
      free(cf);
   }
}


static coast_tile_iterator_t *tc_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat) {

   tiled_file_data_t *cfp = (tiled_file_data_t*)cf->private_data;
   coast_tile_iterator_t *ti = NULL;

   /* Map the rectangular region into the grid space. */
   int lon_start = (int)floor((minlon - cfp->minlon) / cfp->dx);
   int lon_stop = (int)floor((maxlon - cfp->minlon) / cfp->dx);
   int lat_start = (int)floor((minlat - cfp->minlat) / cfp->dy);
   int lat_stop = (int)floor((maxlat - cfp->minlat) / cfp->dy);
   int allow_wrap = (cfp->minlon <= -180) || (cfp->maxlon >= 180.0);

   /* Trim to the extents. */
   if (!allow_wrap) {
      lon_start = min(max(0, lon_start), cfp->nlons - 1);
      lon_stop = min(max(0, lon_stop), cfp->nlons - 1);
   }

   lat_start = min(max(0, lat_start), cfp->nlats - 1);
   lat_stop = min(max(0, lat_stop), cfp->nlats - 1);

   /* Build tile iterator */
   ti = (coast_tile_iterator_t *)malloc(sizeof(coast_tile_iterator_t));
   ti->file = cf;
   ti->allow_wrap = allow_wrap;
   ti->start_lon_index = lon_start;
   ti->end_lon_index = lon_stop;
   ti->start_lat_index = lat_start;
   ti->end_lat_index = lat_stop;
   ti->cur_lon_index = lon_start;
   ti->cur_lat_index = lat_start - 1;
   ti->next_tile = tc_next_tile;
   ti->free_tile = tc_free_tile;

   return ti;
}

static void tc_free_tile_iterator(coast_tile_iterator_t *ti) {
   if (ti != NULL)
      free(ti);
}

static coast_tile_t *tc_next_tile(coast_tile_iterator_t *ti) {

   coast_tile_t *tile = NULL;
   coast_file_t *cf = ti->file;
   tiled_file_data_t *cfp = (tiled_file_data_t*)cf->private_data;
   double xoffset = 0.0;
   double yoffset = 0.0;
   int i = 0, j = 0;

   /* Increment to the next tile, return NULL if none is found */
   ++(ti->cur_lat_index);
   if (ti->cur_lat_index > ti->end_lat_index) {
      ti->cur_lat_index = ti->start_lat_index;
      ++(ti->cur_lon_index);

      if (ti->cur_lon_index > ti->end_lon_index)
         return NULL;
   }

   i = ti->cur_lon_index;
   j = ti->cur_lat_index;


   /* Wrap the current index if required. */
   if (ti->allow_wrap) {
      int ni = cfp->nlons;
      int nj = cfp->nlats;

      if (i < 0)
         xoffset = -360.0 * ((abs(i) / ni) + 1);

      if (i >= ni)
         xoffset = 360.0 * (i / ni);

      if (j < 0)
         yoffset = -180.0 * ((abs(j) / nj) + 1);

      if (j >= nj)
         yoffset = 180.0 * (j / nj);

      i = i % ni;

      while (i < 0)
         i += ni;

      j = j % nj;

      while (j < 0)
         j += nj;
   }

   if (cfp->indicies[j][i] > 0) {
      FILE *fp = fopen(cf->filename, "rb");
      if (fp != NULL) {
         int ii;

         /* Create the tile and populate with the filled polygons */
         tile = (coast_tile_t *)malloc(sizeof(coast_tile_t));

         fseek(fp, cfp->indicies[j][i], SEEK_SET);
         tile->npolys = read_int(fp);
         tile->polys = (poly_t **)malloc(sizeof(poly_t *)*tile->npolys);
         tile->level = (int *)malloc(sizeof(int)*tile->npolys);

         for (ii=0; ii< tile->npolys; ++ii) {
            int n, jj;

            /* Read the polygon level (coast, lake, island in lake,
             * pond in island in lake, ... */
            tile->level[ii] = read_int(fp);

            /* Build the polygon */
            tile->polys[ii] = poly_create();
            n = read_int(fp);
            for (jj=0; jj<n; ++jj) {
               double x = read_int(fp)*1.0e-6 + xoffset;
               double y = read_int(fp)*1.0e-6 + yoffset;
               poly_add_point(tile->polys[ii], x, y);
            }

            /* Skip shoreline edges. */
            n = read_int(fp);
            if (n > 0)
               fseek(fp, 2*4*n, SEEK_CUR);   /* 4 = sizeof(int) when written */
         }
      }
      fclose(fp);
   }

   return tile;
}


static void tc_free_tile(coast_tile_t *tile) {

   if (tile != NULL) {

      if (tile->level != NULL)
         free(tile->level);

      if (tile->polys != NULL) {
         int i;
         for (i=0; i<tile->npolys; ++i)
            poly_destroy(tile->polys[i]);
         free(tile->polys);
      }

      free(tile);
   }
}

static int read_int(FILE* fp) {
   int value = 0;

   if (!fread(&value, 1, 4, fp))
     quit("error in read_int");
   
#if defined(FLIP)
   value = swabi4 ((unsigned int)value);
#endif

   return value;
}
