/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/cstmesh.h
 *
 *  \brief Header file for cstmesh.c
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: cstmesh.h 6233 2019-05-29 03:31:45Z her127 $
 */

#ifndef _CSTMESH_H
#define _CSTMESH_H

#define IFILE 2
#define OFILE 4
#define REVSE 8
#define LINK_M 1
#define LINK_P 2
#define LINK_N 4
#define LINK_D 8
#define CM_FILE 2
#define CM_LIB  4
#define S_COAST  1
#define S_OBC  2
#define S_LINK  4

typedef struct {
  int *slink;     /* Segment number                                    */
  int *ilink;     /* Index number in segment                           */
  double *llon;   /* Link start/end/mid longitudes                     */
  double *llat;   /* Link start/end/mid latitudes                      */
  double *alon;   /* Acutal segment coordinate                         */
  double *alat;   /* Acutal segment coordinate                         */
  double *segx;   /* Longitude of path for link (optional)             */
  double *segy;   /* Latitude of path for link (optional)              */
  int nseg;       /* Size of link path                                 */
  int dir;        /* Direction to cycle the segment                    */
  int flag;       /* Type of link                                      */
  int id;         /* Link number                                       */
  char text[MAXSTRLEN]; /* Link description (optional)                 */
} link_t;

typedef struct {
  int nobc;       /* Number of points in this obc                      */
  double *olon;   /* Longitude of obc                                  */
  double *olat;   /* latitude of obc                                   */
} cmobc_t;

typedef struct {
  int type;
  int ndims;
  int npoint;
  double **coords;
  int *flag;
  int edge2;
  int **edges;
  int tria3;
  int **triangles;
} msh_t;

typedef struct {
  char geom_file[MAXSTRLEN];
  char mesh_file[MAXSTRLEN];
  int geom_seed;
  int geom_feat;
  int geom_eta1;
  int geom_eta2;
  char hfun_kern[MAXSTRLEN];
  char hfun_scal[MAXSTRLEN];
  char hfun_file[MAXSTRLEN];
  double hfun_hmax;
  double hfun_hmin;
  double hfun_grad;
  int mesh_dims;
  char mesh_kern[MAXSTRLEN];
  int mesh_iter;
  int mesh_top1;
  int mesh_top2;
  double mesh_rad2;
  double mesh_rad3;
  double mesh_eps1;
  double mesh_eps2;
  double mesh_vol3;
  int verbosity;
} jig_t;

typedef struct {
  char cofile[MAXSTRLEN]; /* Coastline file                            */
  char ofile[MAXSTRLEN];  /* JIGSAW .msh geometry input                */
  char pname[MAXSTRLEN];  /* Coastmesh diagnostics                     */
  char hfunf[MAXSTRLEN];  /* JIGSAW mesh size file                     */
  int ismooth;        /* Pre-smoothing passes                          */
  int ismoothz;       /* Pre-smoothing window                          */
  int osmooth;        /* Post-smoothing passes                         */
  int osmoothz;       /* Post-smoothing window                         */
  int nss, **ss;      /* Smoothing individual segments                 */
  int nrs, *sr;       /* Segment removal                               */
  int resample;       /* Resampling interval                           */
  int cutoff;         /* Percentile of ordered sedment sizes below which to cut-off */
  int np;             /* Number of points in the largest segment       */
  double *x;          /* x coordinate of the largest segment           */
  double *y;          /* y coordinate of the largest segment           */
  double slat, slon;  /* Start lat/lon for the largest coastline segment */
  double elat, elon;  /* End lat/lon for the largest coastline segment */
  double mlat, mlon;  /* Mid lat/lon for the largest coastline segment */
  double radius;      /* Minimum radius for segments to include        */
  double length;      /* Minimum length for segments to include        */
  double bslon, bslat, belon, belat; /* Bounding box                   */
  double hfun_min;    /* Minimum mesh size value                       */
  double hfun_max;    /* Maximum mesh size value                       */
  int nobc;           /* Number of open boundaries                     */
  cmobc_t *obc;       /* Open boundary structure                       */
  double auto_l;      /* Auto link length                              */
  int auto_t;         /* Auto link threshold                           */
  double auto_f;      /* Auto link fraction                            */
  int nlink;
  link_t *link;
  msh_t *msh;
  jig_t *jig;
} coamsh_t;

coamsh_t *cm_alloc(void);
void cm_free(coamsh_t *cm);
int coastmesh(coamsh_t *cm, int mode);
void cm_read(coamsh_t *cm, FILE *ip);
void write_jigsaw(jig_t *jig);

#endif
