/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/pointsourcesink.h
 *
 *  \brief Include file for point source/sink code
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: pointsourcesink.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#ifndef _POINTSOURCESINK_H
#define _POINTSOURCESINK_H

#include "ems.h"
#include "timeseries.h"

/*********************************************************************
The point source/sink structure
*********************************************************************/
typedef struct {
  char name[MAXLINELEN];        /* Name of source/sink */
  double x;                     /* Location x coordinate */
  double y;                     /* Location y coordinate */
  double z;                     /* Location z coordinate */
  double zlow;                  /* Low range of z coordinate */
  double zhigh;                 /* High range of z coordinate */
  timeseries_t *loc;            /* Pointer to location time series */
  int v_offset;                 /* Vertical offset flag */
  int flag;                     /* General purpose flag */
  int x_id;                     /* Index of x variable in location time
                                   series */
  int y_id;                     /* Index of y variable in location time
                                   series */
  int zl_id;                    /* Index of z variable in location time
                                   series */
  int zh_id;                    /* Index of z variable in location time
                                   series */
  int *e1;                      /* Integer model grid coordinate */
  int *e2;                      /* Integer model grid coordinate */
  int *e3;                      /* Integer model grid coordinate */
  int vc;                       /* Number model grid coordinates */
  int *iloc;                    /* i horizontal location */
  int *jloc;                    /* j horizontal location */

  /* Pointer to index routine */
  int (*xyzijk) (void *, double, double, double, int *, int *, int *);
  /* Pointer to data needed by index routines */
  void *model_data;

  /* Number of tsfiles associated with this pss */
  int ntsfiles;

  /* Timeseries data */
  timeseries_t ts[MAX_TS_FILES];

  /* Total number of variables across all the ts files */
  int nv;

  /* Tracer map to the index in master */
  int *tmap;

  int watertsid;       /* Time series multi index of water variable */
  int u1tsid;          /* Time series multi index of u1 velocity */
  int u2tsid;          /* Time series multi index of u2 velocity */
  int s2m;             /* Slave to master map */
} pss_t;

/* Prototypes */
void pss_read(char *name, char *t_units, pss_t **pss,
              int *np, void *data,
              int (*xyzijk) (void *, double, double, double, int *, int *,
                             int *), int (*trI) (void *, char *));

#endif                          /* _POINTSOURCESINK_H */
