/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/tf_xyz_cmr.c
 *
 *  \brief Functions to access ASCII XYZ bathy/topo files
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: tf_xyz_cmr.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include "ems.h"

/* external prototypes */
/* - topo */
int get_index(topo_details_t *td, double x, double y);

/* local functions */
/**
 * check to see if file is a valid XYZ ASCII file
 */
int tf_xyz_cmr_isvalid(char *fname) {
   int isvalid = 1;
   FILE *fp;
   char keyword[MAXSTRLEN];
   char vers[MAXSTRLEN];

   fp = fopen(fname, "r");

   if(fp != NULL) {
      /* OK the file is there */
      strcpy(keyword, "VERSION");
      if(prm_read_char(fp, keyword, vers)) {
         /* VERSION keyword is there */
         if(strcmp("xyz_ascii", vers) == 0) {
            /* it is an xyz tuple ASCII file */
         }
      }

      fclose(fp);
   }

   return isvalid;
}

/*
void my_swap_2d(double *vals1, double *vals2, int first, int last) {
   double tmp;

   tmp = vals1[last];
   vals1[last] = vals1[first];
   vals1[first] = tmp;
   tmp = vals2[last];
   vals2[last] = vals2[first];
   vals2[first] = tmp;
}

int my_split_2d(double *vals1, double *vals2, int first, int last) {
   int splitpoint;
   int i;
   double checkval;

   splitpoint = first;
   checkval = vals1[first];
   for(i=first+1;i<=last;++i)
      if(vals1[i] < checkval) {
         ++splitpoint;
         my_swap_2d(vals1, vals2, i, splitpoint);
      }
   my_swap_2d(vals1, vals2, first, splitpoint);

   return splitpoint;
}

void my_q_sort_2d(double *vals1, double *vals2, int first, int last) {
   int splitpoint;
   if(first <last) {
      splitpoint = my_split_2d(vals1, vals2, first, last);
      my_q_sort_2d(vals1, vals2, first, splitpoint-1);
      my_q_sort_2d(vals1, vals2, splitpoint+1, last);
   }
}
*/

topo_details_t *tf_xyz_cmr_open(char *fname) {
   topo_details_t *td;
   FILE *fp;
   char keyword[MAXSTRLEN];
   int npoints, i;
   float x, y, z;

  td = (topo_details_t*)malloc(sizeof(topo_details_t));
  memset(td, 0, sizeof(topo_details_t));
  td->extent = (topo_extent_t*)malloc(sizeof(topo_extent_t));  
  memset(td->extent, 0, sizeof(topo_extent_t));

   fp = fopen(fname, "r");
   strcpy(td->name, fname);
   td->fp = fp;

   strcpy(keyword, "NPOINTS");
   (void)prm_read_int(fp, keyword, &npoints);
   td->nlats = npoints;
   td->nlons = npoints;

   td->lats = d_alloc_1d(npoints);
   td->lons = d_alloc_1d(npoints);
   td->z = d_alloc_1d(npoints);

   /* read in the data */
   for(i=0;i<npoints;++i) {
      if(fgets(keyword, 80, fp) != NULL) {
         sscanf(keyword,"%g%g%g", &x, &y, &z);
         td->lons[i] = (double) x;
         td->lats[i] = (double) y;
         td->z[i] = (double) z;
      } else {
         fclose(fp);
         quit("tf_xyz_cmr_open: NULL occurred whilst reading", fname);
      }
   }
   /* now sort it */
/* no reason to sort!
   my_q_sort_2d(td.lons, td.lats, 0, npoints-1);
*/

   /* set the unrelated values to "non values" */
   td->topo_type = POINT;
   td->delx = -1;
   td->dely = -1;
   td->ncid = -1;
   td->stats = NULL;
   strcpy(td->basedir, "");
   td->stats = NULL;
   td->num_file_tiles = 0;
   td->num_file_rows = 0;
   td->num_file_cols = 0;

   return td;
}

int tf_xyz_cmr_close(topo_t *tf) {
   if(!fclose(tf->td->fp)) {
      /* fclose rets 0 if close OK */
      return 1;
   }

   return 0;
}

int tf_xyz_cmr_getz(topo_details_t *td, double x[], double y[], double z[],
   int nx, int ny) {
   int i;

   /* since we already read in the data all we have to do is find 
    * the nearest points */

   for(i=0;i<nx;++i)
      z[i] = td->z[get_index(td, x[i], y[i])];
      
   return 1;
}

topo_extent_t *tf_xyz_cmr_extent(topo_details_t *td) {
   topo_extent_t *te;
   int i;

  te = (topo_extent_t*)malloc(sizeof(topo_extent_t));  
  memset(te, 0, sizeof(topo_extent_t));

   te->minLat = td->lats[0];
   te->maxLat = td->lats[0];
   te->minLon = td->lons[0];
   te->maxLon = td->lons[0];

   for(i=1;i<td->nlons;++i) {
      te->minLat = min(te->minLat, td->lats[i]);
      te->maxLat = max(te->maxLat, td->lats[i]);
      te->minLon = min(te->minLon, td->lons[i]);
      te->maxLon = max(te->maxLon, td->lons[i]);
   }

   return te;
}

topo_tile_t *tf_xyz_cmr_get_tile(topo_t *tf, double xcoord,
   double ycoord) {
   topo_tile_t *tile;
   int i;
   int nvals, *tmpindxs;

  tile = (topo_tile_t*)malloc(sizeof(topo_tile_t));
  memset(tile, 0, sizeof(topo_tile_t));
  tile->td = (topo_details_t*)malloc(sizeof(topo_details_t));
  memset(tile->td, 0, sizeof(topo_details_t));
  tile->td->extent = (topo_extent_t*)malloc(sizeof(topo_extent_t));  
  memset(tile->td->extent, 0, sizeof(topo_extent_t));
  tile->neighbours = (topo_neighbours_t*)malloc(sizeof(topo_neighbours_t));
  memset(tile->neighbours, 0, sizeof(topo_neighbours_t)); 

   /* set the flag to indicate it is a nonnc file */
   tile->td->ncid = -1;
   tile->td->topo_type = tf->td->topo_type;

   tile->td->nodata = tf->td->nodata;

   /* calc the extent
    * - assume that xcoord, ycoord have been properly calced above */
   tile->td->extent->minLon = xcoord;
   tile->td->extent->minLat = ycoord;
   tile->td->extent->maxLon = xcoord + (tf->td->extent->maxLon 
      - tf->td->extent->minLon) / 5;
   tile->td->extent->maxLat = ycoord + (tf->td->extent->maxLat 
      - tf->td->extent->minLat) / 5;

   /* now find and set the local vals */
   tmpindxs = i_alloc_1d(tf->td->nlons);
   nvals = 0;
   for(i=0;i<tf->td->nlons;++i) {
      /* lets tier this for "faster" finds */
      if(tf->td->lons[i] >= tile->td->extent->minLon)
         if(tf->td->lons[i] < tile->td->extent->maxLon)
            if(tf->td->lats[i] >= tile->td->extent->minLat)
               if(tf->td->lats[i] < tile->td->extent->maxLat) {
                  tmpindxs[nvals] = i;
                  ++nvals;
               }
   }
   tile->td->nlons = nvals+1;
   tile->td->lons = d_alloc_1d(tile->td->nlons);
   tile->td->nlats = nvals+1;
   tile->td->lats = d_alloc_1d(tile->td->nlats);
   tile->td->z = d_alloc_1d(tile->td->nlons);
   for(i=0;i<nvals;++i) {
      tile->td->lons[i] = tf->td->lons[tmpindxs[i]];
      tile->td->lats[i] = tf->td->lats[tmpindxs[i]];
      tile->td->z[i] = tf->td->z[tmpindxs[i]];
   }
   i_free_1d(tmpindxs);

   /* set some non info to non values */
   tile->td->delx = -1;
   tile->td->dely = -1;
   strcpy(tile->td->name, "");
   strcpy(tile->td->basedir,"");
   tile->td->fp = NULL;
   tile->td->fpos = -1;
   tile->td->num_file_tiles = -1;
   tile->td->num_file_rows = -1;
   tile->td->num_file_cols = -1;

   /* set the neighbors to unknown */
   tile->neighbours->topLeft  = -1;
   tile->neighbours->top      = -1;
   tile->neighbours->topRight = -1;
   tile->neighbours->bot      = -1;
   tile->neighbours->botLeft  = -1;
   tile->neighbours->botRight = -1;
   tile->neighbours->left     = -1;
   tile->neighbours->right    = -1;

   return tile;
}



