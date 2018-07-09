/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/poly_coast.c
 *
 *  \brief A coastline reader for the CMR ascii polygon coastline file format
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: poly_coast.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ems.h"

/* public prototypes */
int poly_coast_isvalid(const char *filename);
coast_file_t* poly_coast_open(const char *filename);

/* private prototypes */
static void pc_close(coast_file_t *cf);
static coast_tile_iterator_t *pc_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat);
static void pc_free_tile_iterator(coast_tile_iterator_t *ti);
static coast_tile_t *pc_next_tile(coast_tile_iterator_t *ti);
static void pc_free_tile(coast_tile_t *tile);

int poly_coast_isvalid(const char *filename) {
   int valid = 0;

   FILE *fp = fopen(filename, "r");
   if (fp != NULL) {
      valid = fgetc(fp) == '#';
      fclose(fp);
   }

   return valid;
}

coast_file_t* poly_coast_open(const char *filename) {

   coast_file_t* cf = NULL;

   cf = (coast_file_t *)malloc(sizeof(coast_file_t));
   memset(cf, 0, sizeof(coast_file_t));
   strcpy(cf->filename, filename);
   cf->close = pc_close;
   cf->get_tile_iterator = pc_get_tile_iterator;
   cf->free_tile_iterator = pc_free_tile_iterator;

   return cf;
}

static void pc_close(coast_file_t *cf) {
   if (cf != NULL) {
      free(cf);
   }
}


static coast_tile_iterator_t *pc_get_tile_iterator(coast_file_t *cf,
                    double minlon, double maxlon,
                    double minlat, double maxlat) {

   coast_tile_iterator_t *ti = NULL;

   /* Build tile iterator */
   ti = (coast_tile_iterator_t *)malloc(sizeof(coast_tile_iterator_t));
   memset(ti, 0, sizeof(coast_tile_iterator_t));
   ti->file = cf;
   ti->next_tile = pc_next_tile;
   ti->free_tile = pc_free_tile;
   ti->private_data = (void *)1;

   return ti;
}

static void pc_free_tile_iterator(coast_tile_iterator_t *ti) {
   if (ti != NULL)
      free(ti);
}

typedef struct poly_list poly_list_t;
struct poly_list {
   poly_t *polygon;
   poly_list_t *next;
};

static coast_tile_t *pc_next_tile(coast_tile_iterator_t *ti) {
   coast_tile_t *tile = NULL;

   if (ti->private_data != NULL) {
      FILE *fp = fopen(ti->file->filename, "r");
      if (fp != NULL) {
         int n = 0;
         int npolys = 0;
         poly_t *p = poly_create();
         poly_list_t *head = NULL;
         poly_list_t *last = NULL;

         /* Build a link list of polygons */
         while((n = poly_read(p, fp)) > 0) {
            poly_list_t *list = (poly_list_t *)malloc(sizeof(poly_list_t));
            list->polygon = p;
            list->next = NULL;
            if (head == NULL) {
               head = list;
               last = list;
            } else {
               last->next = list;
               last = list;
            }
            p = poly_create();
            ++npolys;
         }
         poly_destroy(p);
         fclose(fp);

         /* Build an array of polygons */
         if ((npolys > 0) && (head != NULL)) {
            poly_t **polys = (poly_t **)malloc(sizeof(poly_t *)*npolys);
            int *level = (int *)malloc(sizeof(int)*npolys);

            n = 0;
            do {
               poly_list_t *next = head->next;
               polys[n] = head->polygon;
               level[n] = 1;
               free(head);

               head = next;
               ++n;
            } while (head != NULL);


            /* Create and polulate the tile */
            tile = (coast_tile_t *)malloc(sizeof(coast_tile_t));
            memset(tile, 0, sizeof(coast_tile_t));
            tile->npolys = npolys;
            tile->polys = polys;
            tile->level = level;
            ti->private_data = NULL;
         }
      }
   }

   return tile;
}


static void pc_free_tile(coast_tile_t *tile) {

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
