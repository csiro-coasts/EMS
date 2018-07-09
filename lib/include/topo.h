/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/topo.h
 *
 *  \brief Include file for bathymetric reading and interpolation
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: topo.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#include <stdio.h>

#define DODEBUG 1
#define NODEBUG 0
#if !defined(DEBUGMODE)
#define DEBUGMODE DODEBUG
#endif

#define DOSTATS 1
#define NOSTATS 0
#if !defined(STATSMODE)
#define STATSMODE NOSTATS
#endif

#define MAXTILEHEIGHT 128
#define MAXTILEWIDTH 128

/* enumeration types */
typedef enum {
   AUTO, NEAREST, LINEAR, NATURAL, INVERSE, AVERAGE
} topo_hint_t;

typedef enum {
   GRID_POINT, GRID_POLY, POINT, COVERAGE
} topo_type_t;

/* typedefs for topo */
typedef struct sample_stats sample_stats_t;
typedef struct topo_neighbours topo_neighbours_t;
typedef struct topo_extent topo_extent_t;
typedef struct topo_details topo_details_t;
typedef struct topo_tile topo_tile_t;
typedef struct topo_tiles topo_tiles_t;
typedef struct topo topo_t;
typedef struct topo_files topo_files_t;
typedef struct topo_params topo_params_t;
typedef struct topo_file_params topo_file_params_t;

/* structures for topo */
struct sample_stats {
   int n;
   double mean;
   double std;
   double min;
   double max;
};

struct topo_neighbours {
   int topLeft;
   int top;
   int topRight;
   int botLeft;
   int bot;
   int botRight;
   int left;
   int right;
};

struct topo_extent {
   double minLon;
   double minLat;
   double maxLon;
   double maxLat;
};

struct topo_details {
   /* data related members */
   int nlats;
   int nlons;
   double delx;
   double dely;
   double *lats;
   double *lons;
   double *z;
   topo_extent_t *extent;
   sample_stats_t *stats;
   double nodata;
   topo_type_t topo_type;

   /* file related members */
   char name[MAXSTRLEN];
   char basedir[MAXSTRLEN];
   /* - NetCDF uses a file index */
   int ncid;
   /* - others use a file pointer */
   FILE *fp;
   /* - fpos is useful at times */
   long fpos;
   /* - some formats use tiled files */
   int num_file_tiles;
   int num_file_rows;
   int num_file_cols;
};

struct topo_tile {
   topo_details_t *td;
   topo_neighbours_t *neighbours;
};

struct topo_tiles {
   int ntiles;
   topo_tile_t **tiles;
   int last_tile;
};

struct topo {
   topo_details_t *td;
   topo_details_t *(*open)(char *fname);
   topo_extent_t *(*get_extent)(topo_details_t *td);
   int (*getz)(topo_details_t *td, double x[], double y[], double z[], int nx,
      int ny);
   int (*close)(topo_t *tf);
   topo_tiles_t tts;
};

struct topo_files {
   int nfiles;
   topo_t **files;
   int last_file;
};

struct topo_file_params {
   char filename[MAXSTRLEN];

   /* optional, depends on file type */
   int ntiles;
};

struct topo_params {
   int nfiles;
   topo_file_params_t *files;
};

/* prototypes */
/* topo_factory */
topo_t *topo_get_file_from_factory(char *name);

/* topo */
/* for file access use these */
topo_params_t *topo_params_read(FILE *prmfp);

topo_files_t* topo_init(topo_params_t *tp);
topo_files_t* topo_simple_init(char *fnames);
void topo_destroy(topo_files_t *tfs);

/* for data access uses these */
double topo_get_z_at_point(topo_files_t *tfs, double x, double y,
   topo_hint_t hint);
double* topo_get_z_at_points(topo_files_t *tfs, double *x, double *y, int n,
   topo_hint_t hint);
double** topo_get_z_for_grid(topo_files_t *tfs, double **gx, double **gy,
   int nce1, int nce2, topo_hint_t hint);
double topo_get_z_under_area(topo_files_t *tfs, double *x, double *y, int n,
   topo_hint_t hint);

/* for the odd debug */
void draw_map_ascii(int nx, int ny, double *z);


