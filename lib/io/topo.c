/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/io/topo.c
 *
 *  \brief API and functions to access bathy/topo data
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: topo.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include "ems.h"

/* external prototypes */
/* - topo_tile_misc */
topo_tile_t *get_tile(topo_t *tf, double xcoord, double ycoord);
int get_tile_neighbours(topo_tiles_t *tts);
double get_tile_z(topo_details_t *td, double lon, double lat);
int get_tile_index(topo_tiles_t *tts, double lon, double lat);

/* local functions */

/*
 * find the nearest index in the latitude array
 * return the index
 */
int get_lat_index(topo_details_t *td, double lat) {
   int i;
   double delv, minv, maxv;

  if(td->dely <= 0.0)
    quit("get_lat_index: dely is zero or negative\n");

  if(td->topo_type == GRID_POLY) {
    for(i=0;i<(td->nlats-1);++i)
      if(lat >= td->lats[i] && lat <= td->lats[i+1])
        return i;
  } else if(td->topo_type == GRID_POINT) {
    delv = fabs(td->lats[0] - td->lats[1]);
    minv = td->lats[0] - delv;
    maxv = (td->lats[0] + td->lats[1]) / 2.0;
    if(lat >= minv && lat < maxv)
      return 0;

    for(i=1;i<(td->nlats-1);++i) {
      minv = (td->lats[i-1] + td->lats[i]) / 2.0;
      maxv = (td->lats[i] + td->lats[i+1]) / 2.0;
      if(lat >= minv && lat < maxv)
        return i;
    }

    delv = fabs(td->lats[td->nlats-1] - td->lats[td->nlats-2]);
    minv = (td->lats[td->nlats-1] + td->lats[td->nlats-2]) / 2.0;
    maxv = td->lats[td->nlats-1] + delv;
    if(lat >= minv && lat < maxv)
      return td->nlats-1;
  } else
    quit("get_lat_index: non rectlinear data not supported at present\n");

  return -1;
}

/*
 * find the nearest index in the longitude array
 * return the index
 */
int get_lon_index(topo_details_t *td, double lon) {
   int i;
   double delv, minv, maxv;

  if(td->delx <= 0.0)
    quit("get_lat_index: delx is zero or negative\n");

  if(td->topo_type == GRID_POLY) {
    for(i=0;i<(td->nlons-1);++i)
      if(lon >= td->lons[i] && lon <= td->lons[i+1])
        return i;
  } else if(td->topo_type == GRID_POINT) {
    delv = fabs(td->lons[0] - td->lons[1]);
    minv = td->lons[0] - delv;
    maxv = (td->lons[0] + td->lons[1]) / 2.0;
    if(lon >= minv && lon < maxv)
      return 0;

    for(i=1;i<(td->nlons-1);++i) {
      minv = (td->lons[i-1] + td->lons[i]) / 2.0;
      maxv = (td->lons[i] + td->lons[i+1]) / 2.0;
      if(lon >= minv && lon < maxv)
        return i;
    }

    delv = fabs(td->lons[td->nlons-1] - td->lons[td->nlons-2]);
    minv = (td->lons[td->nlons-1] + td->lons[td->nlons-2]) / 2.0;
    maxv = td->lons[td->nlons-1] + delv;
    if(lon >= minv && lon < maxv)
      return td->nlons-1;
  } else
    quit("get_lon_index: non rectlinear data not supported at present\n");

   return -1;
}

/*
 * find the nearest index to a lat/lon pair
 * this assumes that lats and lons are a pair in topo_details_t
 * return the index
 */
int get_index(topo_details_t *td, double x, double y) {
   int i, i_tmp;
   double dist, min_dist;

   min_dist = sqrt(pow(x - td->lons[0], 2.0) + pow(y - td->lats[0], 2.0));
   i_tmp = 0;
   for(i=1;i<td->nlons;++i) {
      dist = sqrt(pow(x - td->lons[i], 2.0) + pow(y - td->lats[i], 2.0));
      if(dist < min_dist) {
         min_dist = dist;
         i_tmp = i;
      }
   }

   return i;
}

/*
 * see if a point is in the defined extent
 * returns 1 if point is in, 0 if out
 */
int in_extent(topo_extent_t *te, double xcoord, double ycoord) {

  if(ycoord < te->minLat)
    return(0);

  if(ycoord > te->maxLat)
    return(0);

  if(xcoord < te->minLon)
    return(0);

  if(xcoord > te->maxLon)
    return(0);

  return(1);
}

/*
 *  find the file based on lat/lon
 */
int get_file_index(topo_files_t *tfs, double x, double y) {
   int file = 0;
   int i;

   /* make sure there are actually files */
   if(tfs->nfiles <= 0)
      return -1;

   /* double check tfs->last_file */
   if(tfs->last_file < 0)
      tfs->last_file = 0;
   else if(tfs->last_file >= tfs->nfiles)
      tfs->last_file = 0;

   /* 1st check the "current" ie tfs->last_file file */
   if(in_extent(tfs->files[tfs->last_file]->td->extent, x, y))
      return tfs->last_file;

   /* now check others */
   for(i=0;i<tfs->nfiles;++i) {
      if(in_extent(tfs->files[file]->td->extent, x, y))
         return file ;
      ++file;
      if(file >= tfs->nfiles)
         file = 0;
   }

   /* else return no file found */
   return -1;
}

/** Initialize the topo_files based on contents of the parameter structure.
 *
 * @param tp a topo_params file list
 * @return a topo_files structure
 */
topo_files_t *topo_init(topo_params_t *tp) {
   topo_files_t *tfs;
   int ntiles = 0;
   int i, j, r, c;
   double xcoord = -91;
   double ycoord = -91;

  /* HACK HACK */
  if(tp->nfiles > 1)
    /* for the time being we only want to allow one file */
    warn("topo_init: only the 1st file will be used\n");

  tfs = (topo_files_t*)malloc(sizeof(topo_files_t));
  memset(tfs, 0, sizeof(topo_files_t));
   tfs->nfiles = tp->nfiles;
   tfs->last_file = 0;
   tfs->files = (topo_t**)malloc(sizeof(topo_t*) * tp->nfiles);
   memset(tfs->files, 0, sizeof(topo_t*) * tp->nfiles);

   for(i=0;i<tp->nfiles;++i) {
      tfs->files[i] = topo_get_file_from_factory(tp->files[i].filename);
      tfs->files[i]->tts.ntiles = tp->files[i].ntiles;
      tfs->files[i]->tts.last_tile = 0;

      tfs->files[i]->td = tfs->files[i]->open(tp->files[i].filename);
      tfs->files[i]->td->extent = tfs->files[i]->get_extent(tfs->files[i]->td);

      /* in case it is a tiled topo file */
      if(tfs->files[i]->td->num_file_tiles > 0) {
         tfs->files[i]->tts.ntiles = tfs->files[i]->td->num_file_tiles;
      }

      ntiles = tfs->files[i]->tts.ntiles;
      tfs->files[i]->tts.tiles = (topo_tile_t **)malloc(
         sizeof(topo_tile_t*) * ntiles);
      memset(tfs->files[i]->tts.tiles, 0, sizeof(topo_tile_t*) * ntiles);

      /* corner not specified so default to lower left corner */
      xcoord = tfs->files[i]->td->extent->minLon;
      ycoord = tfs->files[i]->td->extent->minLat;

      if(tfs->files[i]->td->num_file_tiles > 0) {
         /* it is a tiled file */
         j = 0;
         for(c=0;c<tfs->files[i]->td->num_file_cols;++c) {
            for(r=0;r<tfs->files[i]->td->num_file_rows;++r) {
               tfs->files[i]->tts.tiles[j] = get_tile(tfs->files[i], c, r);
               ++j;
            }
         }
         /* NOTE: neighbours are set during read in */
      } else {
         /* one file now lets get tiles */
         for(j=0;j<ntiles;++j) {
            tfs->files[i]->tts.tiles[j] = get_tile(tfs->files[i], xcoord, 
               ycoord);
            xcoord = tfs->files[i]->td->lons[get_lon_index(
               tfs->files[i]->td, xcoord)
               + tfs->files[i]->tts.tiles[j]->td->nlons];
            if(xcoord >= tfs->files[i]->td->extent->maxLon) {
               xcoord = tfs->files[i]->td->extent->minLon;
               ycoord = tfs->files[i]->td->lats[
                  get_lat_index(tfs->files[i]->td, ycoord)
                  + tfs->files[i]->tts.tiles[j]->td->nlats];
               if(ycoord >= tfs->files[i]->td->extent->maxLat) {
                  ycoord = tfs->files[i]->td->extent->minLat;
               }
            }
         }

         /* now set neighbours */
         if(ntiles > 0)
            if(!get_tile_neighbours(&tfs->files[i]->tts))
               warn(
                 "topo_init: one or more tile neighbours could not be set\n");
      }
   }

   return(tfs);
}

/** Initialize the topo_files structure for the specified topography file.
 *          
 * @param fname The topography file.
 * @return a topo_files structure
 */ 
topo_files_t *topo_simple_init(char *fname) {
  topo_params_t *tp;

  tp = (topo_params_t*)malloc(sizeof(topo_params_t));
  memset(tp, 0, sizeof(topo_params_t));

  tp->nfiles = 1;
  tp->files = (topo_file_params_t *)malloc(
    sizeof(topo_file_params_t));
  memset(tp->files, 0, sizeof(topo_file_params_t));

  strcpy(tp->files[0].filename, fname);
  tp->files[0].ntiles = 1; /* default to one for now */

  return topo_init(tp);
}

/** Read in the parameter file and populate the data datastructure
 *
 * @param prmfp the file pointer to the parameter file
 * @return a populated topo_params data structure or NULL
 */
topo_params_t *topo_params_read(FILE *prmfp) {
   topo_params_t *tp;
   char keyword[MAXSTRLEN];
   int i_tmp;
   int nfiles = 0;
   int i;

   prm_set_errfn(warn);

  tp = (topo_params_t*)malloc(sizeof(topo_params_t));
  memset(tp, 0, sizeof(topo_params_t));

   if(prmfp != NULL) {
      strcpy(keyword, "NBATHYSRC");
      if(prm_read_int(prmfp, keyword, &nfiles) == 0) {
         quit("topo_params_read: Keyword \"%s\" not found\n", keyword);
      }

      if(nfiles <= 0)
         /* in case someone accidently puts in zero or a negative value */
         quit("topo_params_read: No topofiles to load\n");
      

      tp->nfiles = nfiles;
      tp->files = (topo_file_params_t *)malloc(
         sizeof(topo_file_params_t)*nfiles);
      memset(tp->files, 0, sizeof(topo_file_params_t)*nfiles);
      for(i=0;i<nfiles;++i) {
         sprintf(keyword, "BATHYSRC%i.NAME", i);
         if(!prm_read_char(prmfp, keyword, tp->files[i].filename))
           quit("topo_params_read: Keyword \"%s\" not found\n", keyword);

         /* now we do the tiles */
         tp->files[i].ntiles = 1;
         sprintf(keyword, "BATHYSRC%i.MAXTILES", i);
         if(prm_read_int(prmfp, keyword, &i_tmp))
            tp->files[i].ntiles = i_tmp; /* number of tiles defined */
      }
   }

   return tp;
}

/*
 * Given an nx, ny, delx, dely this rountine sets values in the
 * x and y arrays.
 * returns 1 if no errors or 0 if errors occured
 */
int set_x_y(topo_details_t *td, int nx, int ny, double delx, double dely,
   double *x, double *y) {
   int i, j, k;

   if(delx > 0.0 && dely > 0.0) {
      /* a delx and a dely are given */
      k = 0;
      for(i=0;i<ny;++i)
         for(j=0;j<nx;++j) {
            x[k] = x[0] + delx * (double)j;
            y[k] = y[0] + dely * (double)i;
            ++k;
         }
      return 1;
   };

   if(td->delx > 0.0 || td->dely > 0.0)
      quit("set_x_y: DELX and DELY cannot be defined\n");

   /* data in available file */
   k = 0;
   for(i=0;i<ny;++i)
      for(j=0;j<nx;++j) {
         x[k] = x[0] + td->delx * (double)j;
         y[k] = y[0] + td->dely * (double)i;
         ++k;
      }
   return 1;
}

/*
 * routine to read the X, Y, and Zs
 * returns the number of values read
 */
int topo_param_read_data(FILE *prmfp, topo_files_t *tfs, double **x, double **y,
   double **z, sample_stats_t **ss, int do_stats, int do_debug) {
   char keyword[MAXSTRLEN];
   double xcoord = -91;
   double ycoord = -91;
   int do_interp = 0;
   int nx = -1;
   int ny = -1;
   int n = -1;
   double delx = -1;
   double dely = -1;
   int i;
   double xx[4], yy[4];
   int ff, tt;

   /* now let us get some Z data */
   strcpy(keyword, "XCOORD");
   if(!prm_read_double(prmfp, keyword, &xcoord)) {
      /* ok so could be a XCOORDS file */
      strcpy(keyword, "XCOORDS");
      if(!prm_read_darray(prmfp, keyword, &(*x), &nx))
         quit("keyword \"XCOORD\" nor \"%s\" was found\n", keyword);

      strcpy(keyword, "DELX");
      if(!prm_read_double(prmfp, keyword, &delx)) {
         /* in this case it is used for interps */
      }
   } else {
      strcpy(keyword, "NX");
      if(!prm_read_int(prmfp, keyword, &nx))
         quit("Keyword \"%s\" not found\n", keyword);

      strcpy(keyword, "DELX");
      if(!prm_read_double(prmfp, keyword, &delx)) {
         /* at present no DELX means read NX indexs of lat, default of -1
          * already accounts for that. */
      }
   }

   strcpy(keyword, "YCOORD");
   if(!prm_read_double(prmfp, keyword, &ycoord)) {
      /* ok so could be a YCOORDS file */
      strcpy(keyword, "YCOORDS");
      if(!prm_read_darray(prmfp, keyword, &(*y), &ny))
         quit("Keyword \"YCOORD\" nor \"%s\" was found\n", keyword);

      /* at this stage we know the NX and the NY read in */
      if(nx != ny)
         quit("XCOORDS and YCOORDS: NX must equal NY\n");

      strcpy(keyword, "DELY");
      if(!prm_read_double(prmfp, keyword, &dely)) {
         /* in this case it is used for interps */
      }
      (*z) = d_alloc_1d(nx);
   } else {
      strcpy(keyword, "NY");
      if(!prm_read_int(prmfp, keyword, &ny))
         quit("Keyword \"%s\" not found\n", keyword);

      strcpy(keyword, "DELY");
      if(prm_read_double(prmfp, keyword, &dely)) {
         if(delx < 0)
            quit("\"%s\" present without \"DELX\"\n", keyword);
      } else {
         if(delx > 0)
            quit("\"DELX\" present without \"%s\"\n", keyword);
      }
      /* at this stage we know the NX and NY for reading */
      (*x) = d_alloc_1d(nx*ny);
      (*x)[0] = xcoord;
      (*y) = d_alloc_1d(nx*ny);
      (*y)[0] = ycoord;
      (*z) = d_alloc_1d(nx*ny);
   }

   /* check to see if do an interp or not, default is not to */
   strcpy(keyword, "INTERP");
   if(prm_read_int(prmfp, keyword, &do_interp))
      if(do_interp)
         if(delx < 0 || dely < 0) {
            quit("DELX and DELY must be defined when doing interpolations\n");
         }

   if(xcoord > -91) {
      /* xcoord/ycoord set */
      n = nx*ny;
      tfs->last_file = 0;
      ff= get_file_index(tfs, xcoord, ycoord);
      if(ff < 0)
         quit("\"corner\" request not within available extent(s)\n");
      tt = get_tile_index(&tfs->files[ff]->tts, xcoord, ycoord);
      if(tt < 0)
         set_x_y(tfs->files[ff]->td, nx, ny, delx, dely, (*x), (*y));
      else
         set_x_y(tfs->files[ff]->tts.tiles[tt]->td, nx, ny, delx, dely, (*x),
            (*y));
      if(do_interp) {
         if(do_stats) {
            (*ss) = (sample_stats_t *)malloc(sizeof(sample_stats_t) * n);
            memset((*ss), 0, sizeof(sample_stats_t) * n);
         }
         for(i=0;i<n;i++) {
            xx[0] = (*x)[i] - delx/2.0;
            yy[0] = (*y)[i] - dely/2.0;
            xx[2] = (*x)[i] + delx/2.0;
            yy[2] = (*y)[i] + dely/2.0;
            xx[1] = xx[2];
            yy[1] = yy[0];
            xx[3] = xx[0];
            yy[3] = yy[2];
            (*z)[i] = topo_get_z_under_area(tfs, xx, yy, 4, AVERAGE);
         }
      } else {
         (*z) = topo_get_z_at_points(tfs, (*x), (*y), n, NEAREST);
      }
      if(do_debug)
         /* write out the grid */
         draw_map_ascii(nx, ny, (*z));
   } else {
      /* xcoords/ycoords set */
      n = nx;
      tfs->last_file = 0;
      for(i=0;i<n;++i)
         if(get_file_index(tfs, (*x)[i], (*y)[i]) < 0)
            quit("Coord (%.2f, %.2f) not within available extent(s)\n",
               (*x)[i], (*y)[i]);

      if(do_interp) {
         if(do_stats) {
            (*ss) = (sample_stats_t*)malloc(sizeof(sample_stats_t) * n);
            memset((*ss), 0, sizeof(sample_stats_t) * n);
         }
         for(i=0;i<n;i++) {
            xx[0] = (*x)[i] - delx/2.0;
            yy[0] = (*y)[i] - dely/2.0;
            xx[2] = (*x)[i] + delx/2.0;
            yy[2] = (*y)[i] + dely/2.0;
            xx[1] = xx[2];
            yy[1] = yy[0];
            xx[3] = xx[0];
            yy[3] = yy[2];
            (*z)[i] = topo_get_z_under_area(tfs, xx, yy, 4, AVERAGE);
         }
      } else {
         (*z) = topo_get_z_at_points(tfs, (*x), (*y), n, NEAREST);
      }

      /* write out the value triplets */
      if(do_debug || is_log_enabled(LDEBUG))
      {
         for(i=0;i<n;++i)
            emslog(LDEBUG,"x[%i],y[%i],z[%i]: %.2f, %.2f, %.2f\n",
               i,i,i,(*x)[i],(*y)[i],(*z)[i]);
      }
   }

   return n;
}

/** Close the open topo_files
 *
 * @param tfs the active topo_files structure
 */
void topo_destroy(topo_files_t *tfs) {
  int i,j;

  for(i=0;i<tfs->nfiles;++i) {
    if(!tfs->files[i]->close(tfs->files[i])) {
      warn("topo_destroy: Error while closing file \"%s\"\n",
        tfs->files[i]->td->name);
    }

    if(tfs->files[i]->td->lats != NULL) {
      /* assume if one then the other exists as well */
      d_free_1d(tfs->files[i]->td->lats);
      d_free_1d(tfs->files[i]->td->lons);
    }

    if(tfs->files[i]->td->z != NULL)
      d_free_1d(tfs->files[i]->td->z);

    if(tfs->files[i]->td->extent != NULL)
      free(tfs->files[i]->td->extent);

    if(tfs->files[i]->td->stats != NULL)
      free(tfs->files[i]->td->stats);

    free(tfs->files[i]->td);

    for(j=0;j<tfs->files[i]->tts.ntiles;++j) {
      if(tfs->files[i]->tts.tiles[j]->td->lats != NULL) {
        /* assume if one then the other exists as well */
        d_free_1d(tfs->files[i]->tts.tiles[j]->td->lats);
        d_free_1d(tfs->files[i]->tts.tiles[j]->td->lons);
      }

      if(tfs->files[i]->tts.tiles[j]->td->z != NULL)
        d_free_1d(tfs->files[i]->tts.tiles[j]->td->z);

      if(tfs->files[i]->tts.tiles[j]->td->stats != NULL)
        free(tfs->files[i]->tts.tiles[j]->td->stats);

      free(tfs->files[i]->tts.tiles[j]->td);

      free(tfs->files[i]->tts.tiles[j]->neighbours);

      free(tfs->files[i]->tts.tiles[j]);
    }

    free(tfs->files[i]);
  }
}

/*
 * bilinear interp w/ 4 points
 */
double bilinear(double xcoord, double ycoord, double *x, double *y, double *z) {
   double val = NaN;
   double xSpan, ySpan, x2ll, y2ll;

   /*  2 --- 3
    *  |     |
    *  | *   |   * = point of interest
    *  0 --- 1
    * assume the 0 index is the lower left and the 3 index the upper right
    * calc the spans(assming ~rectilinear relationship) */
   xSpan = fabs(x[0] - x[3]);
   ySpan = fabs(y[0] - y[3]);

   /* calc the distance from the point to the lower left corner */
   x2ll = fabs(xcoord - x[0]);
   y2ll = fabs(ycoord - y[0]);

   val = z[0] * (1.0 - (x2ll/xSpan)) * (1.0 - (y2ll/ySpan))
      + z[1] * (x2ll/xSpan) * (1.0 - (y2ll/ySpan)) 
      + z[2] * (1.0 - (x2ll/xSpan)) * (y2ll/ySpan)
      + z[3] * (x2ll/xSpan) * (y2ll/ySpan);

   return(val);
}

/*
 * a simple weighted distance interp for n points
 * xcoord and ycoord represents the point of interest
 */
double inverse_weighted_distance_interp(double xcoord, double ycoord, 
   double *x, double *y, double *z, int n) {
   double val = 0.0;
   double sigmaRh = 0.0;
   double h[n];
   double R = 0.0;
   double w;
   int i;

   /* This is not the normal "Shepards's method" (Shepard 1968) but one by
    * (Franke & Nielson, 1980) as documented in Environmental Modeling Systems,
    * Inc's SMS package.  It is claimed to provide "superior" results.  At least
    * it is not based on an arbitrary power law scaling.  Unfortunately it is
    * of order 3N vs the original of orrder N.
    * The weight for a given point is given by:
    *
    *            /  R - hi  \ 2           R  = the maximum distance a point is
    *            | -------- |                  from the point of interest
    *            \    Rhi   /
    *  wi = ----------------------        hi = the distance of point i from the
    *          n                               point of interest
    *        . -.  /  R - hj  \ 2
    *         >    | -------- |
    *        ` -'  \    Rhj   /
    *         j=1
    *
    * calc the distances */
   for(i=0;i<n;++i) {
      h[i] = sqrt(pow((xcoord - x[i]), 2.0) + pow((ycoord - y[i]), 2.0));
      R = max(R, h[i]);
   }

   /* calc the normalization factor */
   for(i=0;i<n;++i)
      sigmaRh = sigmaRh + pow((R - h[i]) / (R * h[i]), 2.0);

   /* calc the weighted value */
   for(i=0;i<n;++i) {
      w = pow(((R - h[i]) / (R * h[i])), 2.0) / sigmaRh;
      val = val + z[i] * w;
   }

   return val;
}

/*
 * find the nearest n points
 */
int find_nearest_n_neighbours(topo_files_t *tfs, double x, double y,
   int nmax, double *xx, double *yy, double *zz) {
   int n = -1;

   if(n < 0)
      quit("find_nearest_n_neighbours: no neighbours found/not implemented yet\n");

   return n;
}

/*
 * calc the sample stats
 */
sample_stats_t calc_sample_stats(double *vals, int n, double mean) {
   sample_stats_t ss;
   int i;

   ss.n = n;
   ss.mean = mean;
   ss.min = vals[0];
   ss.max = vals[0];
   /* sample variance:
    *         sigma[1->n](Xi - Xmean)^2
    *  S^2 = ---------------------------
    *              nsamples - 1                  */
   for(i=0;i<n;++i) {
      ss.std = ss.std + pow((vals[i] - mean), 2.0);
      ss.min = min(ss.min, vals[i]);
      ss.max = max(ss.max, vals[i]);
   }
   ss.std = ss.std / (double)(ss.n - 1);
   /* now make it the standard deviation */
   ss.std = sqrt(ss.std);

   return(ss);
}

/** Get z value at a particular point
 *
 * @param tfs datastructure that holds topo file data
 * @param x the x/long value
 * @param y the y/lat value
 * @param hint the hint to be used for calculating the z value
 * @return the z value
 */
double topo_get_z_at_point(topo_files_t *tfs, double x, double y,
   topo_hint_t hint) {
   int ff = -1, tt = -1;
   int lon_index = -1, lat_index = -1, n = -1;
   int max_neighbours = 10;
   double xx[max_neighbours], yy[max_neighbours], zz[max_neighbours];

   ff = get_file_index(tfs, x, y);
    /* if the point is not in a file return NaN*/
   if(ff < 0)
     return NaN;

   if(hint == AUTO) {
      /* try to guess */

      /* add more later */
      if(tfs->files[ff]->td->topo_type == GRID_POINT
         || tfs->files[ff]->td->topo_type == GRID_POLY)
         return topo_get_z_at_point(tfs, x, y, NEAREST);
   } else if(hint == NEAREST) {
      /* use nearest neighbour */

      /* check if in a tile */
      tt = get_tile_index(&tfs->files[ff]->tts, x, y);
      if(tt >= 0) {
         /* point is in a tile */
         return get_tile_z(tfs->files[ff]->tts.tiles[tt]->td, x, y);
      } else {
         /* read point from the file */
         double z = NaN;
         if(!tfs->files[ff]->getz(tfs->files[ff]->td, &x, &y, &z,
            1, 1))
            warn("topo_get_z_at_point: file read error\n");
         return z;
      }
   } else if(hint == LINEAR) {
      /* use bilinear interpolation */
      /* at present use the 4 points "around" the point of interest */
      if(tfs->files[ff]->td->topo_type != GRID_POINT
         && tfs->files[ff]->td->topo_type != GRID_POLY)
         quit("topo_get_z_at_point: LINEAR called for non GRID topo file\n");

      /* check if in a tile */
      tt = get_tile_index(&tfs->files[ff]->tts, x, y);
      if(tt >= 0) {
         /* point is in a tile */
         lat_index = get_lat_index(tfs->files[ff]->tts.tiles[tt]->td, y);
         lon_index = get_lon_index(tfs->files[ff]->tts.tiles[tt]->td, x);
         xx[0] = tfs->files[ff]->tts.tiles[tt]->td->lons[lon_index];
         yy[0] = tfs->files[ff]->tts.tiles[tt]->td->lats[lat_index];
         zz[0] = get_tile_z(tfs->files[ff]->tts.tiles[tt]->td, xx[0],
            yy[0]);
         xx[1] = tfs->files[ff]->tts.tiles[tt]->td->lons[lon_index + 1];
         yy[1] = yy[0];
         zz[1] = get_tile_z(tfs->files[ff]->tts.tiles[tt]->td, xx[1],
            yy[1]);
         xx[2] = xx[0];
         yy[2] = tfs->files[ff]->tts.tiles[tt]->td->lats[lat_index + 1];
         zz[2] = get_tile_z(tfs->files[ff]->tts.tiles[tt]->td, xx[2],
            yy[2]);
         xx[3] = xx[1];
         yy[3] = yy[2];
         zz[3] = get_tile_z(tfs->files[ff]->tts.tiles[tt]->td, xx[3],
            yy[3]);
      } else {
         /* read point from the file */
         xx[0] = x;
         yy[0] = y;
         if(!tfs->files[ff]->getz(tfs->files[ff]->td, &xx[0], &yy[0], &zz[0],
            2, 2)) {
            warn("topo_get_z_at_point: file read error\n");
            return NaN;
         }
      }

      return bilinear(x, y, &xx[0], &yy[0], &zz[0]);
   } else if(hint == NATURAL) {
      /* use natural neighbour */
      quit("topo_get_z_at_point: NATURAL NEIGHBOR not yet added\n");

      /* check if in a tile */
      tt = get_tile_index(&tfs->files[ff]->tts, x, y);
      if(tt >= 0) {
         /* point is in a tile */
      } else {
         /* read point from the file */
      }
   } else if(hint == INVERSE) {
      /* use inverse weighted distance */
      quit("topo_get_z_at_point: INVERSE WEIGHTED DISTANCE not yet added\n");

      n = find_nearest_n_neighbours(tfs, x, y, max_neighbours, 
           &xx[0], &yy[0], &zz[0]);

      return inverse_weighted_distance_interp(x, y, &xx[0], &yy[0], &zz[0], n);
   } else {
      /* unsupported hint */
      quit("topo_get_z_at_point: unsupported hint\n");
   }

   return NaN;
}

/** Get z values for a set of points
 *
 * @param tfs datastructure that holds topo file data
 * @param x pointer to the array of the x/long values
 * @param y pointer to the array of the y/lat values
 * @param hint the hint to be used for calculating the z value
 * @param n number of z levels
 * @return pointer to the array of z values
 */
double* topo_get_z_at_points(topo_files_t *tfs, double *x, double *y, int n,
   topo_hint_t hint) {
   int i;
   double *z;

   z = (double*)malloc(sizeof(double)*n);

   for(i=0;i<n;++i)
      z[i] = topo_get_z_at_point(tfs, x[i], y[i], hint);

   return z;
}

/** Get z value under an area (poly)
 * NOTE: at present this assumes an unrotated rectangle
 *
 * @param tfs datastructure that holds topo file data
 * @param x pointer to the array of the x/lon values
 * @param y pointer to the array of the y/lat values
 * @param hint the hint to be used for calculating the z value
 * @param n number of z levels
 * @return the z value
 */
double topo_get_z_under_area(topo_files_t *tfs, double *x, double *y, int n,
   topo_hint_t hint) {
   double z = 0, z_tmp = NaN;
   double *z4stats = NULL;
   double *x_tmp, *y_tmp;
   double ulx, uly, lrx, lry;
   int ff, tt, ttt, i, j, r, c;
   int num_vals, num_down, num_over;
   int nx, ny, ul_lat_index, ul_lon_index, lr_lat_index, lr_lon_index;
   /*
   sample_stats_t ss;
   */

   /*
    * assume for now that the shape is a rectangle and topo file is a grid
    */
   if(n != 4)
      quit("topo_get_z_under_area: rectangular area not specified\n");

   /* assume all points in file */
   ff = get_file_index(tfs, x[0], y[0]);
   /* if the point is not in a file return NaN*/
   if(ff < 0)
      quit("topo_get_z_under_area: area starting at point(x,y): %.2f, %.2f not in file\n",
         x[0], y[0]);

   if(tfs->files[ff]->td->topo_type == POINT
      || tfs->files[ff]->td->topo_type == COVERAGE)
      quit("topo_get_z_under_area: not a GRID type of topo file\n");

   ulx = min(min(min(x[0], x[1]), x[2]), x[3]);
   uly = max(max(max(y[0], y[1]), y[2]), y[3]);
   lrx = max(max(max(x[0], x[1]), x[2]), x[3]);
   lry = min(min(min(y[0], y[1]), y[2]), y[3]);

   /* OK, it is in the file, now check if in a tile */
   tt = get_tile_index(&tfs->files[ff]->tts, ulx, uly);
   if(tt >= 0) {
      /* ul point is in a tile */
      ul_lat_index = get_lat_index(tfs->files[ff]->tts.tiles[tt]->td, uly);
      ul_lon_index = get_lon_index(tfs->files[ff]->tts.tiles[tt]->td, ulx);
      /* now check for lr in same one */
      if(tt == get_tile_index(&tfs->files[ff]->tts, lrx, lry)) {
         /* whew!, in the same tile */
         lr_lat_index = get_lat_index(tfs->files[ff]->tts.tiles[tt]->td,
            lry);
         lr_lon_index = get_lon_index(tfs->files[ff]->tts.tiles[tt]->td,
            lrx);
         nx = lr_lon_index - ul_lon_index;
         ny = ul_lat_index - lr_lat_index;

         z4stats = d_alloc_1d(nx*ny);
         x_tmp = d_alloc_1d(nx*ny);
         y_tmp = d_alloc_1d(nx*ny);
         num_vals = 0;
         for(i=0;i<nx;++i)
            for(j=0;j<ny;++j) {
               z_tmp = get_tile_z(tfs->files[ff]->tts.tiles[tt]->td,
                  tfs->files[ff]->tts.tiles[tt]->td->lons[lr_lon_index + i],
                  tfs->files[ff]->tts.tiles[tt]->td->lats[ul_lat_index + j]);
               if(z_tmp != NaN) {
                  z = z + z_tmp;
                  z4stats[num_vals] = z_tmp;
                  x_tmp[num_vals] = tfs->files[ff]->tts.tiles[tt]->td->lons[
                     lr_lon_index + i];
                  y_tmp[num_vals] =  tfs->files[ff]->tts.tiles[tt]->td->lats[
                     ul_lat_index + j];
                  ++num_vals;
               }
            }

         if(hint == INVERSE) {
            z = inverse_weighted_distance_interp((lrx-ulx)/2.0,
               (uly-lry)/2.0, x_tmp, y_tmp, z4stats, num_vals);
         } else {
            /* default AVERAGE for now */
            z = z / (double)num_vals;
            /* check on do_stats later
            ss = calc_sample_stats(z4stats, num_vals, z/(double)num_vals);
            */
         }

         d_free_1d(z4stats);
         d_free_1d(x_tmp);
         d_free_1d(y_tmp);

         return z;
      } else if(get_tile_index(&tfs->files[ff]->tts, lrx, lry) >= 0){
         /* assume that if the lower right corner is in a tile then the
          * whole area is available via tiles */

         num_down = 0;
         num_over = 0;
         /* check how far "down" from current tile */
         ttt = tfs->files[ff]->tts.tiles[tt]->neighbours->bot;
         for(i=0;i<=tfs->files[ff]->tts.ntiles;++i) {
            if(in_extent(tfs->files[ff]->tts.tiles[ttt]->td->extent,
               ulx, lry))
               break;
            else {
               ++num_down;
               ttt = tfs->files[ff]->tts.tiles[ttt]->neighbours->bot;
               if(ttt < 0)
                  quit("topo_get_z_under_area: all in tile assumption invalid\n");
            }
         }
         /* check how far "over" from current tile */
         ttt = tfs->files[ff]->tts.tiles[tt]->neighbours->right;
         for(i=0;i<=tfs->files[ff]->tts.ntiles;++i) {
            if(in_extent(tfs->files[ff]->tts.tiles[ttt]->td->extent,
               lrx, uly))
               break;
            else {
               ++num_over;
               ttt = tfs->files[ff]->tts.tiles[ttt]->neighbours->right;
               if(ttt < 0)
                  quit("topo_get_z_under_area: all in tile assumption invalid\n");
            }
         }
         /* now we know the tile range */
         num_vals = 0;
         ttt = tt;
         /* assume that all the tiles are the same size */
         nx = tfs->files[ff]->tts.tiles[ttt]->td->nlons * num_over * num_down;
         ny = tfs->files[ff]->tts.tiles[ttt]->td->nlats * num_over * num_down;
         z4stats = d_alloc_1d(nx*ny);
         x_tmp = d_alloc_1d(nx*ny);
         y_tmp = d_alloc_1d(nx*ny);
         nx = 0;
         ny = 0;
         for(c=0;c<=num_over;++c) {
            ul_lat_index = get_lat_index(tfs->files[ff]->tts.tiles[ttt]->td,
               uly);
            /* if uly is not in tile then end at "top" */
            if(ul_lat_index < 0)
               ul_lat_index = tfs->files[ff]->tts.tiles[ttt]->td->nlats;
            lr_lat_index = get_lat_index(tfs->files[ff]->tts.tiles[ttt]->td,
               lry);
            /* if lry not in tile then start at "bot" */
            if(lr_lat_index < 0)
               lr_lat_index = 0;
            ny = ul_lat_index - lr_lat_index;
            for(r=0;r<=num_down;++r) {
               ul_lon_index = get_lon_index(tfs->files[ff]->tts.tiles[
                  ttt]->td, ulx);
               /* if ulx is not in tile then start at "left" */
               if(ul_lon_index < 0)
                  ul_lon_index = 0;
               lr_lon_index = get_lon_index(tfs->files[ff]->tts.tiles[
                  ttt]->td, lrx);
               /* if lrx is not in tile then end at "right" */
               if(lr_lon_index < 0)
                  lr_lon_index = tfs->files[ff]->tts.tiles[ttt]->td->nlons;
               nx = lr_lon_index - ul_lon_index;
               for(i=0;i<nx;++i)
                  for(j=0;j<ny;++j) {
                     z_tmp = get_tile_z(tfs->files[ff]->tts.tiles[ttt]->td,
                        tfs->files[ff]->tts.tiles[ttt]->td->lons[lr_lon_index 
                           + i],
                        tfs->files[ff]->tts.tiles[tt]->td->lats[ul_lat_index
                           + j]);
                     if(z_tmp != NaN) {
                        z = z + z_tmp;
                        z4stats[num_vals] = z_tmp;
                        x_tmp[num_vals] = tfs->files[ff]->tts.tiles[
                           ttt]->td->lons[lr_lon_index + i];
                        y_tmp[num_vals] =  tfs->files[ff]->tts.tiles[
                           ttt]->td->lats[ul_lat_index + j];
                        ++num_vals;
                     }
                  }
               /* now go "down" one */
               ttt = tfs->files[ff]->tts.tiles[ttt]->neighbours->bot;
            }
            /* now go "over" one */
            ttt = tfs->files[ff]->tts.tiles[tt]->neighbours->right;
         }

         if(hint == INVERSE) {
            z = inverse_weighted_distance_interp((lrx-ulx)/2.0,
               (uly-lry)/2.0, x_tmp, y_tmp, z4stats, num_vals);
         } else {
            /* default AVERAGE for now */
            /* check on do_stats later
            ss = calc_sample_stats(z4stats, num_vals, z/(double)num_vals);
            */
            z = z / (double)num_vals;
         }

         d_free_1d(z4stats);
         d_free_1d(x_tmp);
         d_free_1d(y_tmp);

         return z;
      }
   }

   /* not in a tile so read from file */
   ul_lat_index = get_lat_index(tfs->files[ff]->td, uly);
   ul_lon_index = get_lon_index(tfs->files[ff]->td, ulx);
   lr_lat_index = get_lat_index(tfs->files[ff]->td, lry);
   lr_lon_index = get_lon_index(tfs->files[ff]->td, lrx);

   nx = lr_lon_index - ul_lon_index;
   ny = ul_lat_index - lr_lat_index;

   z4stats = d_alloc_1d(nx*ny);
   x_tmp = d_alloc_1d(nx*ny);
   y_tmp = d_alloc_1d(nx*ny);

   num_vals = 0;
   for(i=0;i<nx;++i)
      for(j=0;j<ny;++j) {
         z_tmp = topo_get_z_at_point(tfs,
            tfs->files[ff]->td->lons[lr_lon_index + i],
            tfs->files[ff]->td->lats[ul_lat_index + j],
            AUTO);
         if(z_tmp != NaN) {
            z = z + z_tmp;
            z4stats[num_vals] = z_tmp;
            x_tmp[num_vals] = tfs->files[ff]->td->lons[
               lr_lon_index + i];
            y_tmp[num_vals] =  tfs->files[ff]->td->lats[
               ul_lat_index + j];
            ++num_vals;
         }
      }

   if(hint == INVERSE) {
      z = inverse_weighted_distance_interp((lrx-ulx)/2.0,
         (uly-lry)/2.0, x_tmp, y_tmp, z4stats, num_vals);
   } else {
      /* default AVERAGE for now */
      /* check on do_stats later
      ss = calc_sample_stats(z4stats, num_vals, z/(double)num_vals);
      */
      z = z / (double)num_vals;
   }

   d_free_1d(z4stats);
   d_free_1d(x_tmp);
   d_free_1d(y_tmp);

   return z;
}

/** Get z values for a grid
 *
 * @param tfs datastructure that holds topo file data
 * @param gx pointer to the 2D array of the nodal x/long values
 *           gx.size = [nce2 + 1] x [nce1 +1]
 * @param gy pointer to the 2D array of the nodal y/lat values
 *           gy.size = [nce2 + 1] x [nce1 +1]
 * @param nce1 number of cellular "rows"
 * @param nce2 number of cellular "columns"
 * @param hint the hint to be used for calculating the z values
 * @return pointer to 2D array to house the z values aka topo
 *         topo.size = [nce2] x [nce1]
 */
double** topo_get_z_for_grid(topo_files_t *tfs, double **gx, double **gy,
   int nce1, int nce2, topo_hint_t hint) {
   int got_topo_ok = 1;
   int i, j;
   double x_ce, y_ce;
   double **topo, *x, *y, *z;

   x = d_alloc_1d(4);
   y = d_alloc_1d(4);
   z = d_alloc_1d(4);
   topo = d_alloc_2d(nce1, nce2);
   for(i=0;i<nce1;++i)
      for(j=0;j<nce2;++j) {
         /* get the centroid of the cell, a simple approach */
         x_ce = (((gx[j][i] + gx[j][i+1]) / 2.0) 
           + ((gx[j+1][i] + gx[j+1][i+1])/ 2.0)) / 2.0;
         y_ce = (((gy[j][i] + gy[j][i+1]) / 2.0)
           + ((gy[j+1][i] + gy[j+1][i+1])/ 2.0)) / 2.0;

         /* get the data based on the hint */
         if(hint == NEAREST) {
            /* now get the value nearest the centroid */
            topo[j][i] = topo_get_z_at_point(tfs, x_ce, y_ce, hint);
         } else if(hint == AVERAGE || hint == AUTO || hint == INVERSE) {
            /* build the 4 point poly */
            x[0] = gx[j][i];
            y[0] = gy[j][i];
            x[1] = gx[j+1][i];
            y[1] = gy[j+1][i];
            x[2] = gx[j][i+1];
            y[2] = gy[j][i+1];
            x[3] = gx[j+1][i+1];
            y[3] = gy[j+1][i+1];
            topo[j][i] = topo_get_z_under_area(tfs, x, y, 4, hint);
         } else if(hint == LINEAR) {
            /* assume that interpolation is to be based on the corners */
            /* NOTE: only bilinear for now */
            /* NOTE: using NEAREST for now, may want others later */
            x[0] = gx[j][i];
            y[0] = gy[j][i];
            z[0] = topo_get_z_at_point(tfs, x[0], y[0], NEAREST);
            x[1] = gx[j+1][i];
            y[1] = gy[j+1][i];
            z[1] = topo_get_z_at_point(tfs, x[1], y[1], NEAREST);
            x[2] = gx[j][i+1];
            y[2] = gy[j][i+1];
            z[2] = topo_get_z_at_point(tfs, x[2], y[2], NEAREST);
            x[3] = gx[j+1][i+1];
            y[3] = gy[j+1][i+1];
            z[3] = topo_get_z_at_point(tfs, x[3], y[3], NEAREST);
            topo[j][i] = bilinear(x_ce, y_ce, x, y, z);
         } else {
            /* else just pass the cell center on to "at_point" for now */
            topo[j][i] = topo_get_z_at_point(tfs, x_ce, y_ce, hint);
         }
         /* if any reads resulted in a NaN the read was not OK */
         if(topo[j][i] == NaN)
            got_topo_ok = 0;
      }

  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(z);

  return topo;
}

/** Debug routine to draw an ascii land vs water map
 *
 * @param nx number of data in the x direction
 * @param ny number of data in the y direction
 * @param z pointer to the 1D array of the z values
 */
void draw_map_ascii(int nx, int ny, double *z) {
   int i, j;
/* <DS/>: draw update */
   int k, nLevels;
   double zMin, zMax, zSpan;

   //zMin = -999999.0;
   zMin = 0;
   zMax = -zMin;
   nLevels = 10;
   for(i=0;i<nx*ny;++i) {
      //zMin = min(zMin, z[i]);
      zMax = max(zMax, z[i]);
   }
   zSpan = (zMax - zMin) / nLevels;

   printf("%0.0f %0.0f %0.0f\n", zMin, zMax, zSpan);
/* </DS> */

   /* the top of the border */
   printf(".");
   for(i=0;i<nx;++i) {
      if((nx*2)>(80-2)) {
         printf("-");
      } else {
         printf("--");
      }
   }
   printf(".\n");

   /* the map and side borders */
   for(i=(ny-1);i>=0;--i) {
      printf("|");
      for(j=0;j<nx;++j) {
         if(isnan(z[i*nx + j])) {
            /* z = NaN */
            if((nx*2)>(80-2)) {
               printf("X");
            } else {
               printf("><");
            }
         } else if(z[i*nx + j] <= -2) {
            /* z <= 0 */
            if((nx*2)>(80-2)) {
               printf(" ");
            } else {
               printf("  ");
            }
         } else if(z[i*nx + j] <= 2) {
            /* z <= 0 */
            if((nx*2)>(80-2)) {
               printf(".");
            } else {
               printf("..");
            }
         } else {
/* <DS/>: draw update */
            //if((nx*2)>(80-2)) {
            //   printf("*");
            //} else {
            //   printf("[]");
            //}
            for(k=0;k<nLevels;++k) {
               if(z[i*nx + j] < (zMin + (double)(k+1) * zSpan)) {
                  if((nx*2)>(80-2)) {
                     printf("%i", k);
                  } else {
                     printf("%i%i", k, k);
                  }
                  break;
               }
            }
/* </DS> */
         }
      }
      printf("|\n");
   }

   /* the bottom of the border */
   printf("`");
   for(i=0;i<nx;++i) {
      if((nx*2)>(80-2)) {
         printf("-");
      } else {
         printf("--");
      }
   }
   printf("'\n");

   /* legend */
   for(i=0;i<(nx-40)/2;++i) {
      printf(" ");
   }
   if((nx*2)>(80-2)) {
     printf("\" \": <= -2, \".\": <= 2, \"*\": > 2, \"X\": NaN\n");
   } else {
     printf("\"  \": <= -2, \"..\": <= 2, \"[]\": > 2, \"><\": NaN\n");
   }
}
