/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/lagrange.h
 *  
 *  Description: Header for semi-Lagrange interpolation
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: lagrange.h 6595 2020-09-03 03:36:52Z her127 $
 *
 */

#if !defined(_LAGRANGE_H)
#define _LAGRANGE_H

#define I_GLOB  0x001
#define I_FIRST 0x002
#define I_HEX   0x004
#define I_QUAD  0x008
#define I_STICK 0x010
#define I_NRST  0x020
#define I_2D    0x040

#define D_RESET       0x001
#define D_ZERO        0x002
#define D_FILL        0x004

typedef struct {
  int *npt;                /* Number of triangles i-th point belongs to */
  int **pt;                /* Index of j-th triangle i-th point belongs to */
} npts_t;


typedef struct {
  int osl;         /* Interpolation code */
  char trasr[MAXSTRLEN];   /* Interpolation name       */
  int mosl;         /* Interpolation code */
  char momsr[MAXSTRLEN];   /* Interpolation name       */
  int nz;           /* Number of vertical layers */
  int nvec;         /* Number of tracer points */
  delaunay **d;     /* Delaunay structure for velocity */
  delaunay **td;    /* Delaunay structure for tracers */
  int *bot;         /* Bottom index */
  int **u2d;        /* Map from ugrid data to Delaunay index */
  int *d2u;         /* Map from Delaunay surface index to ugrid data */
  int **t2d;        /* Map from ugrid data to Delaunay index */
  int **b2d;        /* Map from Delaunay perimeter index to ugrid index */
  int *dol;         /* Valid layers */
  double *layers;   /* Vertical layer structure */
  double **dzz;     /* Layer thickness */
  double **cx;      /* x streamline location */
  double **cy;      /* y streamline location */
  double **cz;      /* z streamline location */
  int **cl;
  int **ck;
  GRID_SPECS **gsx;
  GRID_SPECS **gsy;
  GRID_SPECS **gsz;
  GRID_SPECS **gst;
  npts_t **nptt;    /* Tracer point triangle structure */
  npts_t **nptm;    /* Momentum point triangle structure */
  int options;
} lagrange_t;

typedef struct {
  meshs_t *meshs;
  lagrange_t *l;
  char vers[MAXSTRLEN];
  char chead[MAXSTRLEN];
  char phead[MAXSTRLEN];
  char prmname[MAXSTRLEN];
  char rcode[MAXSTRLEN];
  char timeunit[MAXSTRLEN];
  char lenunit[MAXSTRLEN];
  char projection[MAXSTRLEN];
  size_t kcentresize;
  size_t kgridsize;
  size_t nMesh2_node;
  size_t nMesh2_edge;
  size_t nMesh2_face;
  size_t nMaxMesh2_face_nodes;
  int oset;
  int si;
  int ns2;
  int npem;
  int type;
  int options;
  double *cellz;
  double *gridz;
  int *npe2;
  int **c2v;
  int *index;
  double *cellx;
  double *celly;
  double *gridx;
  double *gridy;
  double *u1x;
  double *u1y;
  double **x;
  double **y;
  double *bathy;
  double bfill;
} ugrid_t;

ugrid_t *ugrid_open(char *name, int options);
void ugrid_free(ugrid_t *u);
meshs_t *meshs_init(ugrid_t *ugrid);
void create_d_xytoc(ugrid_t *u);
void lagrange_init(ugrid_t *u, int nvec, int *vec);
void lagrange_free(lagrange_t *l);
void delaunay_fill(lagrange_t *l, int vid, double **v, int mode);
void semi_lagrange_source(ugrid_t *ug, double **u, double **v,
			  double **w, double *eta, double dt, int nk,
			  int *kvec, int *nvec, int **vec, int mode);
int semi_lagrange_atxy(ugrid_t *ug, double **u, double **v,
		       double **w, double *eta, double dtin,
		       double *x, double *y, double *z, int mode);
double semi_lagrange_tr(ugrid_t *ug, double **tr, int ki, int ci, int mode);
void test_lagrange(char *infile, char *trname);
#endif
