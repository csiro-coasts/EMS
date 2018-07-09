/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/topo_tile_misc.c
 *
 *  \brief Functions to manipulate tiles of bathy/topo data
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: topo_tile_misc.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "ems.h"

/* external prototypes */
/* - topo */
int get_lat_index(topo_details_t *td, double lat);
int get_lon_index(topo_details_t *td, double lon);
int get_index(topo_details_t *td, double x, double y);
int in_extent(topo_extent_t *te, double xcoord, double ycoord);

/* - tf_nc_cmr */
topo_tile_t *tf_nc_cmr_get_tile(topo_t *tf,
   double xcoord,  double ycoord);

/* - tf_nc_cmr_tiled */
topo_tile_t *tf_nc_cmr_tiled_get_tile(topo_t *tf, int col, int row);
int tf_nc_cmr_tiled_getz(topo_details_t *td,
   double x[], double y[], double z[], int nx, int ny);

/* - tf_xyz_cmr */
topo_tile_t *tf_xyz_cmr_get_tile(topo_t *tf, double xcoord,
   double ycoord);

/* local functions */
double get_tile_z(topo_details_t *td, double x, double y) {
   if(td->topo_type == GRID_POINT || td->topo_type == GRID_POLY)
      if(td->ncid < 0) {
         /* non tiled file */
         return td->z[get_lat_index(td, y) * td->nlons 
            + get_lon_index(td, x)];
      } else {
         /* tiled file */
         double z = NaN;
         if(!tf_nc_cmr_tiled_getz(td, &x, &y, &z, 1, 1))
            quit("get_tile_z: error while reading z from tiled NC file\n");
         return z;
      }
   else if(td->topo_type == POINT)
      return td->z[get_index(td, x, y)];
   else
      quit("get_tile_z: unknown topo_type\n");

   return NaN;
}

/*
 * routine to get a tile from a file
 */
topo_tile_t *get_tile(topo_t *tf, double xcoord, double ycoord) {

   if(tf->td->topo_type == GRID_POINT || tf->td->topo_type == GRID_POLY)
      if(tf->td->num_file_tiles > 0) {
         /* if it is a tiled file then get the tile file info */
         return tf_nc_cmr_tiled_get_tile(tf, (int)xcoord, (int)ycoord);
      } else {
         /* a "normal" nc file */
         return tf_nc_cmr_get_tile(tf, xcoord, ycoord);
      }
   else if(tf->td->topo_type == POINT)
      return tf_xyz_cmr_get_tile(tf, xcoord, ycoord);
   else
      quit("get_tile: unknown topo_type"); 

   return NULL;
}

/*
 *  find the tile based on lat/lon
 */
int get_tile_index(topo_tiles_t *tts, double lon, double lat) {
   int tile = -1;
   int i;

   /* make sure there are actually tiles */
   if(tts->ntiles <= 0) {
      return(-1);
   }

   /* double check lastTile */
   if(tts->last_tile < 0) {
      tts->last_tile = 0;
   } else if(tts->last_tile >= tts->ntiles) {
      tts->last_tile = 0;
   }

   /* 1st check the "current" ie lastTile tile */
   if(in_extent(tts->tiles[tts->last_tile]->td->extent, lon, lat)) {
      return(tts->last_tile);
   }

   /* next check the neighbours */
   /* - right tile */
   if(tts->tiles[tts->last_tile]->neighbours->right >= 0) { 
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile
         ]->neighbours->right]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->right);
      }
   }
   /* - botRight tile */
   if(tts->tiles[tts->last_tile]->neighbours->botRight >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile
         ]->neighbours->botRight]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->botRight);
      }
   }
   /* - bot tile */
   if(tts->tiles[tts->last_tile]->neighbours->bot >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile]->neighbours->bot
         ]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->bot);
      }
   }
   /* - botLeft tile */
   if(tts->tiles[tts->last_tile]->neighbours->botLeft >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile
         ]->neighbours->botLeft]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->botLeft);
      }
   }
   /* - left tile */
   if(tts->tiles[tts->last_tile]->neighbours->left >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile]->neighbours->left
         ]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->left);
      }
   }
   /* - topLeft tile */
   if(tts->tiles[tts->last_tile]->neighbours->topLeft >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile
         ]->neighbours->topLeft]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->topLeft);
      }
   }
   /* - top tile */
   if(tts->tiles[tts->last_tile]->neighbours->top >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile]->neighbours->top
         ]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->top);
      }
   }
   /* - topRight tile */
   if(tts->tiles[tts->last_tile]->neighbours->topRight >= 0) {
      if(in_extent(tts->tiles[tts->tiles[tts->last_tile
         ]->neighbours->topRight]->td->extent, lon, lat)) {
         return(tts->tiles[tts->last_tile]->neighbours->topRight);
      }
   }

   /* else check all the others, but 1st increas lastTile by two as
    * we already checked the current tile and it's neighbours
    * - sub one from ntiles to adjust for the 0 to n-1 indexing */
   tile = min((tts->ntiles-1), (tts->last_tile + 2));
   for(i=0;i<tts->ntiles;++i) {
      if(in_extent(tts->tiles[tile]->td->extent, lon, lat)) {
         return(tile);
      }
      ++tile;
      if(tile >= tts->ntiles) {
         tile = 0;
      }
   }

   /* else return no tile found */
   return(-1);
}

/*
 * find the tile neighbours
 */
int get_tile_neighbours(topo_tiles_t *tts) {
   int got_ok = 1;
   int j;
   double xcoord, ycoord;


   /* 1st check to see if we actually have tiles */
   if(tts->ntiles <= 0) {
      return(0);
   }

   /* now set the neighbours
    *
    *      -1      +1
    *     +---+---+---+     x = tile being processed
    *  +1 |   |   |   |
    *     +---+---+---+
    *     |   | X |   |
    *     +---+---+---+
    *  -1 |   |   |   |
    *     +---+---+---+
    */
   for(j=0;j<tts->ntiles;++j) {
      tts->last_tile = j;
      /* find the left one (x - 1)
       * - make it one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->minLon - tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely in the tile */
      ycoord = tts->tiles[j]->td->extent->minLat + tts->tiles[j]->td->dely;
      /* OK, this is a little circular as how can you get the index if
       * you do not know the neighbours->  As long as all the neighbours 
       * were initially set to -1 then it'll fall through to the for loop */
      tts->tiles[j]->neighbours->left = get_tile_index(tts, xcoord, ycoord);

      /* find the topLeft one (x - 1, y + 1)
       * - make sure we are at least one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->minLon - tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->maxLat + tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->topLeft = get_tile_index(tts, xcoord, ycoord);

      /* find the top one (y + 1)
       * - make sure we are at least one delx in the tile */
      xcoord = tts->tiles[j]->td->extent->minLon + tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->maxLat + tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->top = get_tile_index(tts, xcoord, ycoord);

      /* find the topRight one (x + 1, y + 1)
       * - make sure we are at least one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->maxLon + tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->maxLat + tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->topRight = get_tile_index(tts, xcoord,
         ycoord);

      /* find the right one (x + 1)
       * - make sure we are at least one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->maxLon + tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely in the tile */
      ycoord = tts->tiles[j]->td->extent->minLat + tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->right = get_tile_index(tts, xcoord, ycoord);

      /* find the botRight one (x + 1, y - 1)
       * - make sure we are at least one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->maxLon + tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->minLat - tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->botRight = get_tile_index(tts, xcoord,
         ycoord);

      /* find the bot one (y - 1)
       * - make sure we are at least one delx in the tile */
      xcoord = tts->tiles[j]->td->extent->minLon + tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->minLat - tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->bot = get_tile_index(tts, xcoord, ycoord);

      /* find the botLeft one (x - 1, y - 1)
       * - make sure we are at least one delx out of the tile */
      xcoord = tts->tiles[j]->td->extent->minLon - tts->tiles[j]->td->delx;
      /* - make sure we are at least one dely out of the tile */
      ycoord = tts->tiles[j]->td->extent->minLat - tts->tiles[j]->td->dely;
      tts->tiles[j]->neighbours->botLeft = get_tile_index(tts, xcoord, ycoord);

      /* if no erros yet then check neighbours for errors */
      if(got_ok) {
         if(tts->tiles[j]->neighbours->left == -1
            || tts->tiles[j]->neighbours->topLeft == -1
            || tts->tiles[j]->neighbours->top == -1
            || tts->tiles[j]->neighbours->topRight == -1
            || tts->tiles[j]->neighbours->right == -1
            || tts->tiles[j]->neighbours->botRight == -1
            || tts->tiles[j]->neighbours->bot == -1
            || tts->tiles[j]->neighbours->botLeft == -1) {
            /* one or more of the tile neighbours could not be set. */
            got_ok = 0;
         } /* neighbours == -1 */
      }
   }

   return(got_ok);
}

