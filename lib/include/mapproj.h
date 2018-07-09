/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/mapproj.h
 *
 *  \brief Include file for map projection routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: mapproj.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef MAP_PROJ_H
#define MAP_PROJ_H

typedef struct map_proj map_proj_t;

/* This structure provides the data and functions required to
 * transform from eastings/northings to latitude/longitude and
 * vice versa, for any number of arbitary projections.
 */
struct map_proj {
  double ellip_major;           /* Ellipsoidal major axis. */
  double ellip_flat;            /* Ellipsoidal flattening. */
  double falsex;                /* False eastings. */
  double falsey;                /* False northings. */
  void *private_data;           /* Data private to the map projection */
  void *(*init) (map_proj_t *mp, int nargs, char *args[]);
  /* Projection initialisation */
  void (*free) (void *data);    /* Free the private data */
  void (*forward) (map_proj_t *mp, double lat, double lon,
                   double *east, double *north);
  /* Convert lat/lon to eastings/northings */
  void (*inverse) (map_proj_t *mp, double east, double north,
                   double *lat, double *lon);
  /* Convert easting/northings to lat/lon */
};


/* Map projection prototypes */
map_proj_t *mp_init(int nargs, char *args[]);
void mp_cleanup(map_proj_t *mp);
void mp_forward(map_proj_t *mp, double lat, double lon,
                double *east, double *north);
void mp_inverse(map_proj_t *mp, double east, double north,
                double *lat, double *lon);

/* Geodesy prototypes */
double geod_inv_robbins(double x1, double y1, double z1,
                        double x2, double y2, double z2,
                        double a, double ecs);
double geod_inv_geod_fwd_sodanos(double x1, double y1, double x2,
                                 double y2, double a, double ecs);
double geod_inv_sodanos_angles(double x1, double y1, double x2, double y2,
                               double a, double e2, double *raz);
void geod_fwd_sodanos(double x1, double y1, double az1, double s, double a,
                      double e2, double *x2, double *y2);
double geod_gc_distance(double x1, double y1, double x2, double y2);
void geod_fwd_spherical_rot(double elon, double elat, double plon,
                            double plat, double *alon, double *alat);
void geod_inv_spherical_rot(double alon, double alat, double plon,
                            double plat, double *elon, double *elat);

#endif                          /* MAP_PROJ_H */
