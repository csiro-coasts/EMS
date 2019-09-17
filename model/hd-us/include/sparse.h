/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/sparse.h
 *  
 *  Description:
 *  Primary sparse structures.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sparse.h 6131 2019-03-04 00:58:32Z her127 $
 *
 */

#if !defined(_SPARSE_H)
#define _SPARSE_H

/*------------------------------------------------------------------*/
/* It is necessary to define reference to some structures in        */
/* advance of it's definition, because there are function pointers  */
/* within the structure that use these structures as an argument.   */

typedef struct win_priv win_priv_t;
typedef struct geometry geometry_t;
typedef struct window window_t;
typedef struct master master_t;
typedef struct transport transport_t;
typedef struct ts_point ts_point_t;

#include <stdio.h>

#include "netcdf.h"
#include "xytoij.h"
#include "timeseries.h"
#include "pointsourcesink.h"
#include "tracer.h"
#include "debug.h"
#include "dumpfile.h"
#include "boundary.h"

#define GLOB_BC         0
#define TR_CK           0
#define NUMWIN          1

#define LO_ORDER        1
#define HI_ORDER        2

#define MAP_A     1
#define GRID_A    2
#define ANGLES_A  4
#define WINDOW_A  8
#define MASTER_A  16
#define CENTRE_A  32
#define EDGE_A    64
#define VERTEX_A  128


/*------------------------------------------------------------------*/
/*UR-ADDED late cleanup                                             */
struct free_stack;
typedef struct free_stack free_stack_t;

/*UR-ADDED create a generic stack of objects which selfsufficiently organize
 * themselfs to exchange data, perform cleanup and shutdown
 * the stack can be added to by calling the generic function
 * custom_stack_add (custom_function_t*)
 * and custom_stack_remove(char* name)
 */
struct custom_stack;
typedef struct custom_stack custom_stack_t;

struct custom_function;
typedef struct custom_function custom_function_t;

struct hd_data;
typedef struct hd_data hd_data_t;

struct custom_data;
typedef struct custom_data custom_data_t;

typedef void (*custom_init_funtion) (hd_data_t* hdata, custom_function_t* func);
typedef void (*custom_scatter_funtion) (custom_function_t* func, hd_data_t* hdata);
typedef void (*custom_gather_function) (custom_function_t* func, hd_data_t* hdata);
typedef void (*custom_destroy_function) (void* d);
typedef void (*custom_function_do_gather) (int n, void* target, void* origin);
typedef void (*custom_function_do_scatter) (int n, void* origin, void* target);
typedef void* (*custom_create_data) (void* d);


/*------------------------------------------------------------------*/
/* Used in pt.c                                                     */
/* Structure to allow vertical interpolation of velocity profiles   */
typedef struct {
  double *z;
  double *val;
  int n;
  int kb;
} profile_t;

/*------------------------------------------------------------------*/
/* Global to local sparse map                                       */
/*------------------------------------------------------------------*/
typedef struct {
  int wn;            /* Window corresponding to global cell centre  */
  int we;            /* Window corresponding to global cell edge    */
  int sc;            /* Local cell centre coordinate in window      */
  int ec;            /* Local edge coordinate in window             */
  int vc;            /* Local vertex coordinate in window           */
  int ac;            /* Auxiliary cell centre in window             */
} global_map_t;



typedef struct {
  char name[MAXSTRLEN];         /* Region name */
  int id;                       /* Region identifier */
  int nregions;                 /* Number of regions */
  int nfaces;                   /* Total number of faces connecting regions */
  int nboundaries;              /* Actual number of exchange boundaries */
  int mode;                     /* Operational flags */
  int *bmap;                    /* Exchange boundary map */
  int *rmap;                    /* Reverse exchange boundary map */
  int *gmap;                    /* Window to global boundary map */
  int obcz;                     /* Open boundary zone */
  double t;                     /* Time */
  double pt;                    /* Previous time */
  double dt;                    /* TS output time interval */
  double next;                  /* Time of next output */
  double pdt;                   /* Previous time interval for averaging */
  double step;                  /* Step counter */
  FILE *fp;                     /* Pointer to the timeseries file */
  int nvec;                     /* Number of cells in this region */
  int *vec;                     /* Cells in this region */
  int *nbe1;                    /* Number of boundary cells in e1 direction*/
  int **be1;                    /* Horizontal boundary cells */
  int *nbz;                     /* Number of horizontal boundary cells */
  int **bz;                     /* Boundary cells in z direction */
  int nvar;                     /* Number of 3D tracer output variables */
  int nvarS;                    /* Number of 2D tracer output variables */
  int *var;                     /* 3D tracer output variable indices */
  int *varS;                    /* 2D tracer output variable indices */
  double *mmass;                /* Mean total mass */
  double *stdev;                /* Mass standard deviation */
  double *smass;                /* Mass at start of interval */
  double *emass;                /* Mass at end of interval */
  double *pssf;                 /* Point source flux */
  double *merr;                 /* Mass error (transport mode) */
  double *gfer;                 /* Global fill error (transport mode) */
  double *sedtr;                /* Sediment transport contribution */
  double *eco;                  /* Biogeochemistry contribution */
  double **trflux;              /* Tracer boundary fluxes */
  double **w1;                  /* Work array */
  double **sgne1;               /* Sign for e1 fluxes */
  double **sgne2;               /* Sign for e2 fluxes */
  double *sgnz;                 /* Sign for z fluxes */
  double *dz;                   /* Cell thickness */
  double vfluxp;                /* Sum of positive volume fluxes */
  double vfluxn;                /* Sum of negative volume fluxes */
  double residence;             /* Residence time using max flow */
  double residencen;            /* Residence time using net flow */
  char **fluxname;              /* Name of flux exchange */
  int transferf;                /* Transfer flag */
} region_t;

/* Function to calculate horizontal mixing of momentum.  */
/* Table of mixing scheme names, and pointers to routines */
/* which implement them. */
typedef struct {
  void (*pre) (geometry_t *window, window_t *windat, win_priv_t *wincon);
  void (*setup) (geometry_t *window, window_t *windat, win_priv_t *wincon);
  void (*u1) (geometry_t *window, window_t *windat, win_priv_t *wincon);
  void (*u2) (geometry_t *window, window_t *windat, win_priv_t *wincon);
  void (*setup_av) (geometry_t *window, window_t *windat, 
		    win_priv_t *wincon, double *tzp);
  void (*u1av) (geometry_t *window, window_t *windat,
		win_priv_t *wincon, double *tzp);
  void (*u2av) (geometry_t *window, window_t *windat,
		win_priv_t *wincon, double *tzp);
} mix_method_t;


/* Particle timeseries input for velocity / mortality */
typedef struct {
  double *val;                  /* Velocity or mortality */
  timeseries_t *ts;             /* Timeseries file */
  int id;                       /* TS id for the variable */
  double data;                  /* Private data */
} pt_ts_t;

typedef struct {
  int rlx;                 /* Relaxation flag */
  char rlxn[MAXSTRLEN];    /* Relaxation filename */
  char rlxtc[MAXSTRLEN];   /* Time constant string */
  int tctype;              /* Time constant type */
  double rlxdt;            /* Relaxation input time */
  double rate;             /* Relaxation rate used */
  double dv0, dv1;         /* Adaptive relaxation endpoints */
  double tc0, tc1;         /* Adaptive relaxation endpoints */
  double slope;            /* Adaptive relaxation slope */
  double *val1, *val2;     /* Relaxation arrays */
  int size;                /* Size of relaxation arrays */
} relax_info_t;

typedef struct {
  int ns;                  /* # coordinates */
  int ns2;                 /* # cells */
  int mnpe;                /* Maximum vertices / cell */
  int *npe;                /* Vertices / cell array */
  double *xloc;            /* x coordinates */
  double *yloc;            /* y coordinates */
  int ***eloc;             /* Index to coordinate map */
  int nobc;                /* # Boundaries */
  int *npts;               /* # OBC cells */
  int **loc;               /* OBC index # */
  int ***obc;              /* OBC locations */  
  int nce1, nce2;          /* Cartesian grid size (optional) */
  int *iloc;               /* i location (optional) */
  int *jloc;               /* j location (optional) */
  int **neic;              /* Cell neighbour map */
  int **neij;              /* Cell neighbor index map */
} mesh_t;

typedef struct {
  int *npt;                /* Number of triangles i-th point belongs to */
  int **pt;                /* Index of j-th triangle i-th point belongs to */
} npt_t;


/*------------------------------------------------------------------*/
/* Window geometry data structure. Contains all variables which are */
/* not dependent on or can be derived from information contained in */
/* another window.                                                  */
/*------------------------------------------------------------------*/
struct win_priv {
  int ntr;                      /* Number of 3D tracers */
  int ntrS;                     /* Number of 2D tracers */
  int nsed;                     /* Number of sediment tracers */
  int ntre;                     /* Number of 3D edge tracers */
  int ntreS;                    /* Number of 2D edge tracers */
  int ntbdy;                    /* Number of tracers to set OBC's */
  int ic;                       /* 2D loop counter */
  int *advect;                  /* Flag to advect tracers */
  int *diffuse;                 /* Flag to diffuse tracers */
  int ntdif_h;                  /* No. of tracers to horiz. diffuse */
  int *tdif_h;                  /* Array of tracers to horiz. diffuse */
  int ntdif_v;                  /* No. of tracers to vert. diffuse */
  int *tdif_v;                  /* Array of tracers to vert. diffuse */
  int ntdec;                    /* Number of tracers to decay */
  int *tdec;                    /* Array of tracers to decay */
  int ndia;                     /* Number of diagnostic tracers */
  int *diagn;                   /* Diagnostic tracers */
  int *sflux;                   /* Surface flux indicies */
  double *mintr;                /* Minimum tracer values */
  double *maxtr;                /* Maximum tracer values */
  double *dectr;                /* Tracer decay constants */
  char **trname;                /* Name of the tracer */
  int *tbdy;                    /* Array of tracers to set OBC's */
  int *twin;                    /* Order to process the windows */
  int *relax;                   /* Tracers to undergo relaxing */
  int nrlx;                     /* Number of tracers to relax */
  int *reset;                   /* Tracers to undergo resetting */
  int nres;                     /* Number of tracers to reset */
  int *reset2d;                 /* 2D Tracers to undergo resetting */
  int nres2d;                   /* Number of 2D tracers to reset */
  double *neweta;               /* Height of free surface at t+1 (m) */
  double *oldeta;               /* Height of free surface at t-1 (m) */
  double *dz;                   /* Layer thickness at cell center (m) */
  double *Ds;                   /* SIGMA scaling factor (m) */
  double *Hs;                   /* SIGMA water depth (m) */
  double *Hn1;                  /* SIGMA scaling factor (m) */
  double *Hn2;                  /* SIGMA scaling factor (m) */
  double *rdens;                /* Horizontal reference density */
  double *Cd;                   /* Drag coefficient */
  double *u1c1;                 /* e1 advection term constant */
  double *u1c3;                 /* e1 advection metric constant */
  double *u1c4;                 /* e1 advection metric constant */
  double *u1c5;                 /* e1 Coriolis constant */
  double *u1c6;                 /* e1 Coriolis constant */
  double *u2c1;                 /* e2 advection term constant */
  double *u2c3;                 /* e2 advection metric constant */
  double *u2c4;                 /* e2 advection metric constant */
  double *u2c5;                 /* e2 Coriolis constant */
  double *u2c6;                 /* e2 Coriolis constant */
  double *topdensu1;            /* Surface density on e1 faces */
  double *topdensu2;            /* Surface density on e2 faces */
  double *densavu1;             /* Mean density on e1 faces */
  double *densavu2;             /* Mean density on e2 faces */
  double *u1adv;                /* e1 dispersion terms */
  double *u2adv;                /* e2 dispersion terms */
  double *u1inter;              /* Vertically averaged u1 terms */
  double *u2inter;              /* Vertically averaged u1 terms */
  double *u1kh;                 /* Horizontal e1 diffusivity (m2s-1) */
  double *u2kh;                 /* Horizontal e2 diffusivity (m2s-1) */
  double *u1vh;                 /* Horizontal e1 viscosity (m2s-1) */
  double *u2vh;                 /* Horizontal e2 viscosity (m2s-1) */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  double *mdx;                  /* Scaling e1 favtor for sigma */
  double *mdy;                  /* Scaling e2 favtor for sigma */
  double *q;                    /* Turbulent kinetic energy */
  double zs;                    /* Surface roughness length scale (m) */
  double *z0;                   /* Bottom roughness length scale (m) */
  double vz0;                   /* Background vertical viscosity (m2s-2) */
  double kz0;                   /* Background vertical diffusivity  */
  double *eta_rlx3d;            /* Eta relaxation flux for 3d mode */

  double tstart;                /* Start time of the simulation */
  double rampstart;             /* Start time of ramp period */
  double rampend;               /* End time of ramp period */
  double g;                     /* Acceleration due to gravity (ms-2) */
  double tz;                    /* Time zone (seconds) */
  double ambpress;              /* Ambient atmospheric pressure */
  double hmin;                  /* Minimum cell thickness (m) */
  double uf;                    /* Background friction velocity */
  double quad_bfc;              /* Quadratic bottom friction coeff. */
  double min_tke;               /* k-epsilon minimum TKE */
  double min_diss;              /* k-epsilon minimum dissipation */
  double Lmin;                  /* Minimum stratified mixing length */
  double eparam;                /* Tuning parameter for stratification */
  double kz_alpha;              /* Csanady coefficient */
  double vz_alpha;              /* Csanady coefficient */
  double wave_alpha;            /* Alpha parameter for waves */
  double wave_hf;               /* Scaling factor for significant wave height */
  double wave_b1;               /* b1 parameter for waves */
  int smooth_VzKz;              /* Shuman smoothing of Vz and Kz */
  int fcf;                      /* Flag for k-w Wilcox (1988)/(1998) models */
  int trasc;                    /* Advection scheme type flag (tracers) */
  int trasf;                    /* Flages for tracer advection */
  char trasr[MAXSTRLEN];        /* Lagrange interpolation rule for tracers */
  char momsr[MAXSTRLEN];        /* Lagrange interpolation rule for momentum */
  int momsc;                    /* Advection scheme type flag (velocity) */
  int momsc2d;                  /* Advection scheme type flag (velocity) */
  int ultimate;                 /* ULTIMATE filter flag */
  int osl;                      /* Order of semi-Lagrange */
  int nosl;                     /* Number of elements in weight array */
  int mosl;                     /* Middle elements in weight array */
  int visc_method;              /* Horizontal diffusion method */
  double smagorinsky;           /* Smagorinsky horizontal diffusion flag */
  double sue1, sue2;            /* Momentum Smagorinsky coefficients */
  double kue1, kue2;            /* Tracer Smagorinsky coefficients */
  double bsue1, bsue2;          /* Base momentum mixing */
  double bkue1, bkue2;          /* Base tracer mixing */
  int smagcode;                 /* Variables using Smagorinsky */
  int diff_scale;               /* Horizontal mixing scaling */
  int stab;                     /* Stability compensation method */
  int thin_merge;               /* Thin layer merging flag */
  int sigma;                    /* Sigma coordinate flag */
  int nonlinear;                /* Non-linearity flag */
  int calc_dens;                /* Calculate density flag */
  int mode2d;                   /* Run in 2D mode */
  int heatflux;                 /* Type of heatflux specification */
  int saltflux;                 /* Type of saltflux specification */
  int cfl;                      /* CFL time-step diagnostic */
  int rampf;                    /* Variables to ramp flag */
  int totals;                   /* Plots integrated totals */
  int robust;                   /* Robust scale for ROAM */
  int fatal;                    /* Instability check flag */
  int smag_smooth;              /* Smagorinsky smoothing */
  int porusplate;               /* Include porus plate sub-gridscale parameterisation */
  double cfl_dt;                /* Time for active cfl specification */
  double lnm;                   /* Level of no motion for steric height */
  int mixlayer;                 /* Mixed layer depth diagnostic */
  int trflux;                   /* Tracer number for flux calculation */
  int trfd1;                    /* Face #1 for flux calculation */
  int trfd2;                    /* Face #2 for flux calculation */
  int trperc;                   /* Tracer no. for percentile calculation */
  int trflsh;                   /* Flushing tracer flag */
  char trage[MAXSTRLEN];        /* Age tracer flag */
  int tendf;                    /* Momentum tendency flag */
  int trtend;                   /* Tracer tendency flag */
  int trsplit;                  /* Split advection using tratio */
  int means;                    /* Mean velocity diagnostic */
  int da;                       /* Data assimilation */
  int nprof;                    /* Normalized profile flag */
  double means_dt;              /* Mean velocity averaging interval */
  double means_next;            /* Next time for zeroing means */
  double means_os;              /* Offset for restarts */
  int means_tra;                /* Mean tracer number */
  int vorticity;                /* Vorticity diagnostic */
  int numbers, numbers1;        /* Numbers diagnostic */
  int save_force;               /* Save input forcing */
  int csr_tide;                 /* Tide constants computed for grid */
  int water_type;               /* Water type */
  int vinit;                    /* Velocity initialisation flag */
  int bulkf;                    /* Bulk scheme for heatflux */
  double zref;                  /* Reference height */
  double hf_ramp;               /* Heat flux ramp time */
  double hftc;                  /* Time constant for heatflux temp */
  double slip;                  /* Slip condition */
  double etamax;                /* Maximum elevation (m) */
  double velmax;                /* Maximum horizontal velocity (ms-1) */
  double velmax2d;              /* Maximum horizontal 2D velocity (ms-1) */
  double etadiff;               /* Maximum difference in mean sea level (m) */
  double wmax;                  /* Maximum vertical velocity (ms-1) */
  double imass;                 /* Initial mass in flushing region */
  short *mask;                  /* Flushing mask */
  int u1_f;                     /* Flag to ommit u1 momentum terms */
  int u1av_f;                   /* Flag to ommit u1av momentum terms */
  int orbital;                  /* Compute orbital wave speed */
  int waves;                    /* Include wave enhanced bottom friction */
  double hmean1;                /* Mean grid spacing between centres (m) */
  double hmean2;                /* Mean grid spacing of edges (m) */
  double amean;                 /* Mean area (m^2) */
  double u1vh0;                 /* Horizontal e1 viscosity (m2s-1) */
  double u2vh0;                 /* Horizontal e2 viscosity (m2s-1) */
  double u1kh0;                 /* Horizontal e1 diffusivity (m2s-1) */
  double u2kh0;                 /* Horizontal e2 diffusivity (m2s-1) */
  int etarlx;                   /* Eta relaxation flag */
  int velrlx;                   /* Velocity relaxation flag */
  int exmapf;                   /* Explicit map flag */
  int alertf;                   /* Perform alert tracking */
  int fillf;                    /* Fill method                        */
  int pssinput;                 /* Point ss input method              */
  int conserve;                 /* Volume conservation */
  int do_closure;               /* Do vertical mixing (transport mode) */
  int do_pt;                    /* Do particle tracking */
  int compatible;               /* Backwards compatible flag */
  int sh_f;                     /* Data input type for specific humidity */
  int filter;                   /* Filtering options */
  int trout;                    /* Transport file output flag */
  int swr_type;                 /* Type of attenuation */
  int dhwf;                     /* Degree heating diagnostic */
  double albedo;                /* Albedo for swr */

  /* Alert thresholds */
  double amax;
  double hmax;
  double vmax;
  double btmax;
  double bcmax;
  double cmax;
  double detamax;
  double dwmax;
  double dtmax;
  double dsmax;
  double smax;
  double emean;
  /* Alert flags */
  int eta_f;
  int vel2d_f;
  int vel3d_f;
  int wvel_f;
  int tend_f;
  int div2d_f;
  int div3d_f;
  int cfl_f;
  int ts_f;
  int tide_r;
  int shear_f;
  int hdiff_f;

  /* Debugging */
  int dbc;   /* Sparse coordinate to debug */
  int dbw;   /* Window containing dbc */
  int dbj;   /* Edge direction to debug */
  FILE *dbf; /* Debugging file handle */
  int dbgf;                     /* Debug flag */
  double dbgtime;               /* Time to commence debugging */

  char timeunit[MAXSTRLEN];     /* Time units */
  tracer_info_t *trinfo_3d;     /* 3D Tracer info data structure */
  tracer_info_t *trinfo_2d;     /* 2D Tracer info data structure */
  tracer_info_t *trinfo_sed;    /* Sediment Tracer info data structure */

  int gint_errfcn;              /* Generic interface error handling flag */
  int *gint_errorf;             /* Generic interface library error flag */
  char **gint_error;            /* Generic interface library error text */

#if defined(HAVE_TRACERSTATS_MODULE)
  void *trs;
#endif

  /* Ecology */
#if defined(HAVE_ECOLOGY_MODULE)
  void *e;                      /* Pointer to ecology structure */
  int do_eco;                   /* Ecology model flag */
  double ecodt;                 /* Ecology timestep */
  int eco_timestep_ratio;       /* Ecological: physical timestep */
#endif

  /* Sediments */
#if defined(HAVE_SEDIMENT_MODULE)
  int do_sed;                   /* Sediment flag */
  sediment_t *sediment;         /* Pointer to sediment structure */
#endif

  /* Waves */
#if defined(HAVE_WAVE_MODULE)
  int do_wave;                  /* Wave flag */
  double wavedt;                /* Wave timestep */
  wave_t *wave;                 /* Pointer to wave structure */
#endif
  double **fetch;               /* Fetch (m) */
  double *coriolis;             /* Coriolis parameter */

  /* Number of OMP threads to use in transport mode */
#ifdef HAVE_OMP
  int trans_num_omp;
#endif

  /* 
   * Work arrays 
   * 
   * The n versions are used in some of the OpenMP parallel routines. eg. advect_diffuse
   */
  double *w1, **w1n;            /* 3D work array #1 */
  double *w2, **w2n;            /* 3D work array #2 */
  double *w3, **w3n;            /* 3D work array #3 */
  double *w4, **w4n;            /* 3D work array #4 */
  double *w5;                   /* 3D work array #5 */
  double *w6;                   /* 3D work array #6 */
  double *w7;                   /* 3D work array #7 */
  double *w8;                   /* 3D work array #8 */
  double *w9;                   /* 3D work array #9 */
  double *w10;                  /* 3D work array #10 */
  double *d1;                   /* 2D work array #1 */
  double *d2;                   /* 2D work array #2 */
  double *d3;                   /* 2D work array #3 */
  double *d4;                   /* 2D work array #4 */
  double *d5;                   /* 2D work array #5 */
  double *d6;                   /* 2D work array #6 */
  double *d7;                   /* 2D work array #7 */
  double *d8;                   /* 2D work array #8 */
  double *v1;                   /* 1D work array #1 */
  double *v2;                   /* 1D work array #2 */
  double *v3;                   /* 1D work array #3 */
  double *v4;                   /* 1D work array #4 */
  double *v5;                   /* 1D work array #5 */
  double *v6;                   /* 1D work array #6 */
  double *v7;                   /* 1D work array #7 */
  double *v8;                   /* 1D work array #8 */
  double *v9;                   /* 1D work array #9 */
  double *v10;                  /* 1D work array #10 */
  double *v11;                  /* 1D work array #11 */
  double *v12;                  /* 1D work array #12 */
  double *sd1;                  /* 1D sediment work array #1 */
  double *one;                  /* 2D work array set to 1.0 */
  double *tendency;             /* Buffer array for tendency diagnostics */
  double **wgt;                 /* Semi-Lagrange weights */
  double *dzz;                  /* cellz thickness for Lagrange */
  double *cellz;                /* cellz levels for Lagrange */
  int **lmap;                   /* Semi-Lagrange map */
  int *nlmap;                   /* Size of lmap */
  double *tmass;                /* Tracer total wc mass */
  double *tsmass;               /* Tracer total sediment mass */
  double *p1, *p2;              /* Dummy pointers */
  double b1, b2, b3, b4;        /* Floating dummies */
  double a1, a2, a3, a4;        /* Floating dummies */
  double *ba;                   /* Floating dummy array */
  int *s1;                      /* 3D integer work array */
  int *s2;                      /* 3D integer work array */
  int *s3;                      /* 3D integer work array */
  int *s4;                      /* 3D integer work array */
  int *s5;                      /* 3D integer work array */
  char *c1;                     /* 3D byte work array */
  char *c2;                     /* 2D byte work array */
  int *m2d;                     /* 3D to 2D map */
  int *gmap;                    /* Ghost cell map */
  int *i1;                      /* 2D integer work array */
  int *i2;                      /* 2D integer work array */
  int *i3;                      /* 2D integer work array */
  int *i4;                      /* 2D integer work array */
  int *i5;                      /* 2D integer work array */
  int *i6;                      /* 2D integer work array */
  int *i7;                      /* 2D integer work array */
  int *kth_e1;                  /* Thin layer vector for e1 velocity */
  int nkth_e1;                  /* Number of cells in kthin, e1 */
  int *kth_e2;                  /* Thin layer vector for e2 velocity */
  int nkth_e2;                  /* Number of cells in kthin, e2 */
  int *cdry_e1;                 /* Dry cell vector for e1 velocity */
  int ncdry_e1;                 /* Number of cells in cdry, e1 */
  int *cdry_e2;                 /* Dry cell vector for e2 velocity */
  int ncdry_e2;                 /* Number of cells in cdry, e2 */
  int *cbot_e1;                 /* Bottom cell vector for e1 velocity */
  int ncbot_e1;                 /* Number of cells in cbot, e1 */
  int *cbot_e2;                 /* Bottom cell vector for e2 velocity */
  int ncbot_e2;                 /* Number of cells in cbot, e2 */
  int vc;                       /* Scalar dummy */
  int vcs;                      /* Scalar dummy */
  int vca;                      /* Scalar dummy */
  int vc1, vcs1, vca1;          /* Scalar dummy */
  int vc2, vcs2, vca2, vci2;    /* Scalar dummy */
  int ncl;                      /* Thin layer num of cells to process */
  int acl;                      /* Non-advective num cells to process */
  int aclS;                     /* Non-advective surface cells to process */
  int dolin_u1;                 /* Linear boundary cells required for u1 */
  int dolin_u2;                 /* Linear boundary cells required for u2 */
  int dobdry_u1;                /* Bdry cells to process required for u1 */
  int dobdry_u2;                /* Bdry cells to process required for u2 */
  int *linmask_u1;              /* Mask for linear u1 calculation */
  int *linmask_u2;              /* Mask for linear u1 calculation */
  int *obcmap;                  /* Maps wet cells to their OBCs */
  int *obctr;                   /* Amalgamated tracers OBCs */
  int tn;                       /* Tracer number */
  pss_t *pss;
  int npss;
  short *agemsk;                /* Age tracer mask */
  short *percmsk;                /* ATracer percentile tracer mask */

  /* FCT work arrays */
  double *Fxh, *Fzh;
  double *Ax, *Az;
  double **ac0, **ac1;
  /* High order advection work arrays */
  double ***B;
  int **Bcell, *nBcell;
  double *trp, *trm;
  /* Linear least squares */
  double ***V;
  int **Vcell, *nVcell;
  qweights *lw;

  /* FFSL work arrays */
  double *crfxf, *crfyf, *crfzf;   /* Fractional factors of total trajectories for faces */
  double *crfxc, *crfyc, *crfzc;   /* Fractional factors of total trajectories for cell centres */
  int *clxf, *clyf, *clzf;         /* Cell counter at head of trajectory, for faces */
  int *clxc, *clyc, *clzc;         /* Cell counter at head of trajectory, for cell centres */
  double *nu, *nv, *nw;            /* Cell-centered velocities */  
  // n versions are for openmp parallel functions
  double *tr_mod, **tr_modn;       /* Modified tracer concentrations */
  double *tr_mod_x, **tr_modn_x;   /* Modified tracer concentrations */
  double *tr_mod_y, **tr_modn_y;   /* Modified tracer concentrations */
  double *tr_mod_z, **tr_modn_z;   /* Modified tracer concentrations */

  /* Function to calculate Vz and Kz.  */
  void (*calc_closure) (geometry_t *, window_t *, win_priv_t *);
  void (*s_func) (double aN, double aM, double *cmu, double *cmu_d);

  /* Function to calculate horizontal mixing of momentum.  */
  mix_method_t *hor_mix;

  /* Function to calculate vorticity for momentum advection */
  double (*pv_calc) (window_t *, int e, int eoe, int v1, int v2);

  /* Tidal harmonic structure */
  tidal_consts_t tc;

  /* Multi-grid transport variables */
  int tmode;                /* Transport mode */
  int tgrid;                /* = 1 if taget & source grids differ       */
  int togn;                 /* Origin for tri-linear interpolation */
  double *xinit;            /* x offset of target origin in source cell */
  double *yinit;            /* y offset of target origin in source cell */

};


/*------------------------------------------------------------------*/
/* Sparse array geometry and navigation structure                   */
/*------------------------------------------------------------------*/
struct geometry {
  /* Structure pointers */
  window_t *windat;             /* Window data structure */
  win_priv_t *wincon;           /* Private window data */

  /* Sparse array dimensions */
  int sgnum;                    /* Number of 3D sparse grid cells */
  int sgnumS;                   /* Number of 2D sparse grid cells */
  int sgsiz;                    /* Size of the 3D sparse array = sgnum+1 */
  int sgsizS;                   /* Size of the 2D sparse array = sgnum2D+1
                                 */
  int nobc;                     /* Number of open boundaries */
  int nce1;                     /* Grid dimension in the e1 direction */
  int nfe1;                     /* nce1 + 1 */
  int nce2;                     /* Grid dimension in the e2 direction */
  int nfe2;                     /* nce2 + 1 */
  int nz;                       /* Number of layers */
  double *layers;               /* Depths of the layer faces */
  int nwindows;                 /* Number of windows */
  double *win_size;             /* Window partition */
  int show_win;                 /* Create plot of windowing */
  int wn;                       /* Window number */
  int *nwn;                     /* Number of cells in each window */
  int **wnx, **wny;             /* List of cells for each window */
  int cmx;                      /* Surface sparse coordinate of max depth */

  /* Re-ordered organisation of the sparse grid */
  int ewet;                     /* End location of wet 3D cells */
  int snon;                     /* Start location of non-wet 3D cells */
  int enon;                     /* End location of non-wet 3D cells */
  int ewetS;                    /* End location of wet 2D cells */
  int snonS;                    /* Start location of non-wet 2D cells */
  int enonS;                    /* End location of non-wet 2D cells */
  int ngsed;                    /* Start location of sedimet ghosts */
  int bdry;                     /* Number of dry cells beneath the sea bed */

  /* Sparse to Cartesian maps. These arrays are currently defined */
  /* for the global sparse array only.  */
  int *s2i;                     /* Sparse to i Cartesian map */
  int *s2j;                     /* Sparse to j Cartesian map */
  int *s2k;                     /* Sparse to k Cartesian map */
  unsigned long ***map;         /* Map from Cartesian to sparse space */
  int *mgc;                     /* Linked list of multiple ghost cells */

  /* Global - window map */
  /* For the global sparse array this contains all wet cells ordered */
  /* consecutively (used to define the window partitions). For windows this */
  /* contains the local to global map. */
  int *wsa;                     /* Array of global 3D cells to process. */
  int *wse;                     /* Array of global 3D edges to process. */
  int *wsv;                     /* Array of global 3D vertices to process. */
  int **g2we;                   /* Global to window map for edges */
  int **g2wv;                   /* Global to window map for vertices */
  global_map_t *fm;             /* Global to local window map. For windows
                                   fm.sc */
                                /* contains the global to local map its window. */

  /* Unstructured */
  int us_type;                  /* Unstructured type */
  int *npe;                     /* Number of nodes per element */
  int npem;                     /* Maximum number of nodes per element */
  int *nve;                     /* Number of nodes surrounding a vertex */
  int nvem;                     /* Maximum number of edges surrounding a vertex */
  int *nvc;                     /* Number of centres surrounding a vertex (2d) */
  int nvcm;                     /* Maximum number of centres surrounding a vertex */
  int *nee;                     /* Number of unique edges surrounding an edge */
  int neem;                     /* Maximum number of unique edges surrounding an edge */
  int **c2c;                    /* Centre to centre spatial map */
  int **c2e;                    /* Centre to edge map */
  int **e2c;                    /* Edge to centre map */
  int **e2e;                    /* Edge to edge map */
  int *ep;                      /* Edge to edge map in e2e[1] direction */
  int *em;                      /* Edge to edge map in e2e[0] direction */
  int **v2e;                    /* Vertex to edge map */
  int **e2v;                    /* Edge to vertex map */
  int **c2v;                    /* Centre to vertex map */
  int **v2c;                    /* Vertex to centre map */
  int *zm1e, *zp1e;             /* Vertical edge maps */
  int *zm1v, *zp1v;             /* Vertical vertex maps */
  int *m2de;                    /* Surface map for edges */
  int *m2dv;                    /* Surface map for vertices */
  int **eSc;                    /* Sign of edge relative to centre */
  int **eSv;                    /* Sign of edge relative to vertex */
  int **eSe;                    /* Unique edges surrounding an edge */
  double **wAe;                 /* Weights of unique edges surrounding an edge */
  double **wSe;                 /* Sign of cross product */
  int **vIc;                    /* Index of 'kite' */
  double *dualarea;             /* Area of the dual cell */
  double **dualareap;           /* Partial area of the dual cell */
  int *e2ijk;                   /* Edge to (i,j,k) map for structured grids */
  int *v2ijk;                   /* Vertex to (i,j,k) map for structured grids */
  int *v2i;                     /* Vertex to i map for structured grids */
  int *v2j;                     /* Vertex to j map for structured grids */
  int *e2k;                     /* Edge to vertical layer map */
  int *v2k;                     /* Vertex to vertical layer map */

  int szc;                      /* Cell centred array size */
  int sze;                      /* Edge centred array size */
  int szv;                      /* Vertex centred array size */
  int szcS;                     /* 2D Cell centred array size */
  int szeS;                     /* 2D Edge centred array size */
  int szvS;                     /* 2D Vertex centred array size */
  int szm;                      /* Maximum 3D size */
  int szmS;                     /* Maximum 2D size */
  int *tri2c;                   /* Triangulation to unstructured map */
  int *c2tri;                   /* Unstructured map to triangulation */
  delaunay *d;                  /* Delaunay data structure for xytoi */
  char i_rule[MAXSTRLEN];       /* Unstructured interpolation method */

  /* 3D cells to process arrays for tracers, u1 and u2 velocities.  */
  /* The sizes v3_* are the limits of cells in w3_* which are to be */
  /* updated in the sparse array. For windows, the sizes n3_* are */
  /* sizes of the vectors and include the first lateral auxiliary */
  /* cell neighbouring wet cells which are required to specify */
  /* fluxes. For the global sparse array a3_t is the total number */
  /* of wet cells (including e1 and e2 boundary faces).  */
  int *w3_t;                    /* 3D local array of tracer cells to
                                   process */
  int *w3_e1;                   /* 3D local array of u1 vel. cells to
                                   process */
  int *w3_e2;                   /* 3D local array of u2 vel. cells to
                                   process */
  int *gsed_t;                  /* Ghost cells for the sediment layer */
  int *ised_t;                  /* Corresponding interior cells to gsed_t */
  int n3_t;                     /* Number of cells in w3_t */
  int n3_e1;                    /* Number of cells in w3_e1 */
  int n3_e2;                    /* Number of cells in w3_e1 */
  int v3_t;                     /* Number of wet cells in w3_t */
  int v3_e1;                    /* Number of wet cells in w3_e1 */
  int v3_e2;                    /* Number of wet cells in w3_e1 */
  int b3_t;                     /* Number of wet + open boundary cells in
                                   w3_t */
  int b3_e1;                    /* Number of wet + open boundary cells in
                                   w3_e1 */
  int b3_e2;                    /* Number of wet + open boundary cells in
                                   w3_e1 */
  int a3_t;                     /* Number of wet + auxiliary cells in w3_t
                                 */
  int a3_e1;                    /* Number of wet + auxiliary cells in
                                   w3_e1 */
  int a3_e2;                    /* Number of wet + auxiliary cells in
                                   w3_e1 */
  int x3_e1;                    /* Number of wet + all auxiliary cells */
  int x3_e2;                    /* Number of wet + all auxiliary cells */

  /* 2D cells to process arrays for elevation, u1 and u2 velocity.  */
  /* The sizes v2_* are the limits of cells in w2_* which are to be */
  /* updated in the window. For windows, the sizes n2_* are sizes */
  /* of the vectors and include the first lateral auxiliary cell */
  /* neighbouring wet cells which are required to specify fluxes.  */
  /* For the global sparse array a2_t is the total number of wet */
  /* cells (including e1 and e2 boundary faces).  */
  int *w2_t;                    /* 2D local array of tracer cells to
                                   process */
  int *w2_e1;                   /* 2D local array of u1 vel. cells to
                                   process */
  int *w2_e2;                   /* 2D local array of u2 vel. cells to
                                   process */
  int *c2cc;                    /* Reverse sparse coordinate to index map */
  int *cc2s;                    /* Input mesh index to sparse coordinate */
  int n2_t;                     /* Number of cells in w2_t */
  int n2_e1;                    /* Number of cells in w2_e1 */
  int n2_e2;                    /* Number of cells in w2_e2 */
  int v2_t;                     /* Number of wet cells in w2_t */
  int v2_e1;                    /* Number of wet cells in w2_e1 */
  int v2_e2;                    /* Number of wet cells in w2_e1 */
  int b2_t;                     /* Number of wet + open boundary cells in
                                   w2_t */
  int b2_e1;                    /* Number of wet + open boundary cells in
                                   w2_e1 */
  int b2_e2;                    /* Number of wet + open boundary cells in
                                   w2_e1 */
  int a2_t;                     /* Number of wet + auxiliary cells in w2_t
                                 */
  int a2_e1;                    /* Number of wet + auxiliary cells in
                                   w2_e1 */
  int a2_e2;                    /* Number of wet + auxiliary cells in
                                   w2_e1 */
  int x2_e1;                    /* Number of wet + all auxiliary cells */
  int x2_e2;                    /* Number of wet + all auxiliary cells */

  /* Ghost cell data. These are currently only defined for the */
  /* global sparse array.  */
  /* Size of the ghost cell arrays */
  int ns2;                      /* 2D sparse size */
  int ns3;                      /* 3D sparse size */
  int nbpt;                     /* Number of 3D lateral land boundary
                                   cells */
  int nbptS;                    /* Number of 2D lateral land boundary
                                   cells */
  int nbpte1;                   /* Number of 3D lateral e1 boundary cells */
  int nbpte1S;                  /* Number of 2D lateral e1 boundary cells */
  int nbpte2;                   /* Number of 3D lateral e2 boundary cells */
  int nbpte2S;                  /* Number of 2D lateral e2 boundary cells */
  int nbe1;                     /* Number of normal velocity e1 3D cells */
  int nbe1S;                    /* Number of normal velocity e1 2D cells */
  int nbe2;                     /* Number of normal velocity e2 3D cells */
  int nbe2S;                    /* Number of normal velocity e2 2D cells */
  /* Solid lateral boundary masks */
  int *bpt;                     /* Solid boundary sparse location */
  int *bin;                     /* Interior cells to solid boundaries */
  int *bin2;                    /* Two interior cells to solid boundaries */
  int *dbin;                    /* Direction of the interior cell */
  int *dbpt;                    /* Direction of the ghost cell */
  int *bpte1;                   /* u1 face centered lateral boundary mask */
  int *bine1;                   /* One point into the interior from bpte1 */
  int *bpte2;                   /* u2 face centered lateral boundary mask */
  int *bine2;                   /* One point into the interior from bpte2 */
  int *bpte1S;                  /* As for bpte1 but for 2D arrays */
  int *bine1S;                  /* As for bine1 but for 2D arrays */
  int *bpte2S;                  /* As for bpte2 but for 2D arrays */
  int *bine2S;                  /* As for bine2 but for 2D arrays */
  int *brsm;                    /* Boundary radiation stress mask */ 
  int *wgst;                    /* Ghost cells for windows */
  int *cask;                    /* Cell centre codes */
  int *eask;                    /* Edge codes */

  /* Cells which need to be transferred from master to slave.  */
  /* These arrays are only defined for windows.  */
  int *m2s;                     /* Master to slave transfer vector */
  int nm2s;                     /* Size of m2s */
  int nm2se1;                   /* Size of m2se1 */
  int nm2sS;                    /* Size of transfer vector for 2D arrays */
  int nm2se1S;                  /* Size of transfer vector for 2D arrays */
  int *m2se1;                   /* Master to slave transfer vector, e1
                                   velocity */
  int *m2se2;                   /* Master to slave transfer vector, e2
                                   velocity */
  int *s2me1;                   /* Slave to master transfer vector, e1
                                   velocity */
  int *s2me2;                   /* Slave to master transfer vector, e2
                                   velocity */
  int *s2m;                     /* Slave to master transfer vector */
  int *s2m_ts;                  /* Time series mapping vector */
  int ns2m;                     /* Size of s2m */
  int ns2me1;                   /* Size of s2me1 */
  int ns2mS;                    /* Size of transfer vector for 2D arrays */
  int ns2me1S;                  /* Size of transfer vector for 2D arrays */
  int ns2m_ts;                  /* Size of ns2m_ts */
  short **owc;                  /* Mapping of slave - master OBC's */

  /* Surface and bottom vectors. These have the sizes of v2_* */
  int *sur_t;                   /* Tracer local vector of surface cells at
                                   t */
  int *nsur_t;                  /* Tracer local vector of surface cells at
                                   t+1 */
  int *bot_t;                   /* Local vector of bottom cells for
                                   tracers */
#if defined(HAVE_SEDIMENT_MODULE)
  int *sed_t;                   /* Sparse locations of the sediment ghosts
                                 */
#endif
  int *sur_e1;                  /* Local vector of surface cells for
                                   tracers */
  int *bot_e1;                  /* Local vector of bottom cells for
                                   tracers */
  int *sur_e2;                  /* Local vector of surface cells for
                                   tracers */
  int *bot_e2;                  /* Local vector of bottom cells for
                                   tracers */
  int *m2d;                     /* 3D to 2D map for windows */
  int cdeep;                    /* Deepest location in the domain */

  /* Local spatial maps. Note: vertical maps for auxiliary cells */
  /* are self-mapping.  */
  int *xp1;                     /* Local array of i+1 locations */
  int *xm1;                     /* Local array of i-1 locations */
  int *yp1;                     /* Local array of j+1 locations */
  int *ym1;                     /* Local array of j-1 locations */
  int *zp1;                     /* Local array of k+1 locations */
  int *zm1;                     /* Local array of k-1 locations */
  int *xpym1;                   /* Local array of i+1,j-1 locations */
  int *xmyp1;                   /* Local array of i-1,j+1 locations */

  /* Open boundary structures */
  open_bdrys_t **open;
  int *obce1;                   /* OUTSIDE e1 OBC cells */
  int nobce1;                   /* Size of obce1 */
  int *obce2;                   /* OUTSIDE e2 OBC cells */
  int nobce2;                   /* Size of obce2 */

  /* Horizontal grid geometry */
  double *h1acell;              /* Cell center e1 grid spacing (m) */
  double *h2acell;              /* Cell center e2 grid spacing (m) */
  double **hacell;              /* Cell center e2 grid spacing (m) */
  double *h1au1;                /* e1 centered e1 grid spacing (m) */
  double *h2au2;                /* e2 centered e2 grid spacing (m) */
  double *h2au1;                /* e1 centered e2 grid spacing (m) */
  double *h1au2;                /* e2 centered e1 grid spacing (m) */
  double *cellarea;             /* Cell area (m2) */
  double *dHde1;                /* Bottom slope in the x direction */
  double *dHde2;                /* Bottom slope in the y direction */
  double *cellx;                /* x co-ordinates of eta points (m) */
  double *celly;                /* x co-ordinates of eta points (m) */
  double *gridx;                /* x co-ordinates of corner points (m) */
  double *gridy;                /* x co-ordinates of corner points (m) */
  double *u1x;                  /* x co-ordinates of u1 points (m) */
  double *u1y;                  /* y co-ordinates of u1 points (m) */
  double *u2x;                  /* x co-ordinates of u2 points (m) */
  double *u2y;                  /* y co-ordinates of u2 points (m) */
  double *thetau1;              /* angle between e1 & x axis at u1 point */
  double *thetau2;              /* angle between e2 & y axis at u2 point */
  double *sinthcell;            /* sin(angle) at cell centre */
  double *costhcell;            /* cos(angle) at cell centre */
  double *sinthu1;              /* sin(angle between e1 & x axes) at u1 */
  double *costhu1;              /* cos(angle between e1 & x axes) at u1 */
  double *sinthu2;              /* sin(angle between e1 & x axes) at u2 */
  double *costhu2;              /* cos(angle between e1 & x axes) at u2 */
  double totarea;               /* Total area */

  /* Vertical grid geometry */
  double *botz;                 /* Cell center height of sea floor (m) */
  double *botzu1;               /* e1 face height of sea floor (m) */
  double *botzu2;               /* e2 face height of sea floor (m) */
  double *botzgrid;             /* Grid corner height of sea floor (m) */
  double *gridz;                /* Vertical level position at face (m) */
  double *cellz;                /* Vertical level position at center (m) */
  double topgrid;               /* Uppermost layerface */

  /* Transport */
  transport_t *trans;

  int is_geog;                  /* Grid specified in geographical projection */
  
  /* Multiple time-stepping arrays for tracers */
  int naux_t;
  int *aux_t;
  double *taux_t;

#if defined(HAVE_ECOLOGY_MODULE)
  /* Ecology */
  int *emap;
#endif

  /* Sediments */
  int sednz;                    /* Number of sediment layers */
  double **cellz_sed;           /* z co-ordinates of sediemnt grid cell */
  double **gridz_sed;           /* z co-ordinates of sediemnt grid level */

  /* Sub-surface explicit maps */
  int neim, neimS;
  int *eims;
  int *eimd;
  int *emmask;
  int *sm_e1;
  int *sm_e2;

  /* Regions */
  int nregions;
  region_t **region;

  /* Process exclusion */
  int ncbgc, *cbgc;
  int ncsed, *csed;
  int ncwave, *cwave;
  int nctran, *ctran;
  int nctrst, *ctrst;

  /* Multi-grid transport maps */
  int *mgm;           /* Maps source cell centre to the target origin */
  GRID_SPECS **gs;
  int **c2p;          /* Cell to Delaunay index map */
  int **gcm;          /* Ghost cell mask */
  xytoij_tree_t *xyij_tree;     /* A xytoij_tree_t for XYtoIJ routines */

  int compatible;               /* Backwards compatible flag */
  int *sask;                    /* Stencil mask */
  double b1, b2, b3, b4;
};


/*------------------------------------------------------------------*/
/* Input parameters data structure                                  */
/*------------------------------------------------------------------*/
typedef struct {
  FILE *prmfd;                  /* Parameter file pointer */

  /* Grid data */
  int nce1;                     /* Grid dimension in the e1 direction */
  int nce2;                     /* Grid dimension in the e2 direction */
  int nfe1;                     /* Grid face dimension in the e1 direction
                                 */
  int nfe2;                     /* Grid face dimension in the e2 direction
                                 */
  int nz;                       /* Number of layers */
  int ns2;                      /* 2D sparse size */
  int ns3;                      /* 3D sparse size */
  double *layers;               /* Depths of the layer faces */
  char gridtype[MAXSTRLEN];     /* Grid name */
  char start_time[MAXSTRLEN];   /* Start time as written in input */
  char stop_time[MAXSTRLEN];    /* Stop time as written in input */
  char oname[MAXSTRLEN];        /* Output parameter file name */
  char opath[MAXSTRLEN];        /* Output path for files */
  char trkey[MAXSTRLEN];        /* Transport keyname */
  double runno;                 /* Unique run identification number */
  char rev[MAXSTRLEN];          /* Version number for parameter file */
  int gridcode;                 /* Code to specify grid type */
  grid_rect_t *rg;              /* Rectangular grid information */
  grid_polar_t *pg;             /* Polar grid information */
  double flon, flat;            /* False pole coordinates */
  double etamax;                /* Maximum elevation (m) */
  double velmax;                /* Maximum horizontal velocity (ms-1) */
  double velmax2d;              /* Maximum horizontal 2D velocity (ms-1) */
  double etadiff;               /* Maximum difference in mean sea level (m) */
  double wmax;                  /* Maximum vertical velocity (ms-1) */
  double bmin;                  /* Minimum bathymetry (m) */
  double bmax;                  /* Maximum bathymetry (m) */
  char mct[MAXSTRLEN];          /* Minimum cell thickness */
  double *bathy;                /* Bathymetry array */
  double **topo;                /* Bathymetry array */
  int nvals;                    /* Number of bathymetry values */
  int reset_bathyf;             /* Reset bathy in -p mode */
  double **x;                   /* Grid e1 coordinate array */
  double **y;                   /* Grid e2 coordinate array */
  double **h1;                  /* Grid e1 metric array */
  double **h2;                  /* Grid e2 metric array */
  double **a1;                  /* Grid e1 angle array */
  double **a2;                  /* Grid e2 angle array */
  int *nwn;
  int **wnx, **wny;

  /* Unstructured */
  int npe;                      /* Number of nodes per element */
  int *npe2;                    /* 2D number of nodes per element */
  int nbu;                      /* Number of open boundaries */
  int *nptsu;                   /* Number of cells per OBC */
  int **locu;                   /* OBC location list */
  double ***posx;               /* OBC x location list */
  double ***posy;               /* OBC y location list */
  int *ic, *jc;                 /* Gridded (i,j) list */
  int ncx, ncy;                 /* Gridded dimensions */
  delaunay *d;                  /* Delaunay data structure */
  char i_rule[MAXSTRLEN];       /* Unstructured interpolation method */
  point *pin;                   /* Points for triangulation */
  int np;                       /* Number of triangulation points */
  int *sin;                     /* Triangulation segments */
  int ns;                       /* Number of segments */
  int uscf;                     /* Unstructured conversion flag */
  mesh_t *mesh;                 /* Input mesh */
  int oset;                     /* Input file indexing offset */
  double mlat, mlon;            /* Interior point for expansion */
  int **neic;                   /* Neighbour cell mapping */
  int **neij;                   /* Neighbour face mapping */
  unsigned long **flag;

  /* Flags */
  int runmode;                  /* Operation mode */
  int nwindows;                 /* Number of windows */
  double *win_size;             /* Window partition */
  int win_reset;                /* Number of steps to reset window loads */
  int show_win;                 /* Create plot of windowing */
  int win_type;                 /* Type of partitioning */
  int win_block;                /* Blocking dimension */
  char win_file[MAXSTRLEN];     /* Window map file; read */
  char wind_file[MAXSTRLEN];    /* Window map file; write */
  char dp_mode[MAXSTRLEN];      /* Distributed processing mode */
  int trasc;                    /* Advection scheme type flag (tracers) */
  int momsc;                    /* Advection scheme type flag (velocity) */
  int momsc2d;                  /* Advection scheme type flag (velocity) */
  int otras;                    /* Order of the advection scheme */
  int osl;                      /* Order of semi-Lagrange */
  int ultimate;                 /* ULTIMATE filter flag */
  char smag[MAXSTRLEN];         /* Smagorisnky input */
  double smagorinsky;           /* Smagorinsky horizontal diffusion flag */
  double sue1, sue2;            /* Momentum Smagorinsky coefficients */
  double kue1, kue2;            /* Tracer Smagorinsky coefficients */
  double bsue1, bsue2;          /* Base momentum mixing */
  double bkue1, bkue2;          /* Base tracer mixing */
  int diff_scale;               /* Horizontal diffusion scaling method */
  int visc_method;              /* Horizontal diffusion method */
  int stab;                     /* Stability compensation method */
  int thin_merge;               /* Thin layer merging flag */
  int sigma;                    /* Sigma coordinate flag */
  int nonlinear;                /* Non-linearity flag */
  int calc_dens;                /* Calculate density flag */
  char densname[MAXSTRLEN];     /* Name of density tracer */
  int mode2d;                   /* Run in 2D mode */
  int cfl;                      /* CFL time-step diagnostic */
  char cfl_dt[MAXSTRLEN];       /* Time for active cfl specification */
  int mixlayer;                 /* Mixed layer depth diagnostic */
  int show_layers;              /* Create plot of sigma layers */
  char trflux[MAXSTRLEN];       /* Tracer number for flux calculation */
  int trfd1;                    /* Face #1 for flux calculation */
  int trfd2;                    /* Face #2 for flux calculation */
  char trperc[MAXSTRLEN];       /* Input for percentile calc. */
  char trpercn[MAXSTRLEN];      /* Tracer name for percentile calc. */
  char trpercr[MAXSTRLEN];      /* Region name for percentile calc. */
  int trflsh;                   /* Flushing tracer flag */
  char trage[MAXSTRLEN];        /* Age tracer flag */
  char dhw[MAXSTRLEN];          /* Degree heating day diagnostic */
  char dhdf[MAXSTRLEN];         /* Degree heating day file */
  double dhw_dt;                /* DHW update increment */
  int dhwf;                     /* Degree heating diagnostic flag */
  int tendf;                    /* Momentum tendency flag */
  char trtend[MAXSTRLEN];       /* Tracer tendency flag */
  int means;                    /* Mean velocity diagnostic */
  char means_dt[MAXSTRLEN];     /* Mean velocity averaging interval */
  char means_os[MAXSTRLEN];     /* Offset for restarts */
  char means_tra[MAXSTRLEN];    /* Offset for restarts */
  char regions[MAXSTRLEN];      /* Name of regions file */
  char region_dt[MAXSTRLEN];    /* Regions interval */
  char region_vars[MAXSTRLEN];  /* Region variables */
  char region_mode[MAXSTRLEN];  /* Region mode */
  int region_obc;               /* Region OBC zone */
  int decorr;                   /* Decorrelation length scale diagnostic */
  int decf;                     /* Decorrelation variable flag */
  char decv[MAXSTRLEN];         /* Decorrelation variable */
  char decs[MAXSTRLEN];         /* Decorrelation scaling */
  int sharp_pyc;                /* Pycnocline sharpening for ROAM */
  int vorticity;                /* Vorticity diagnostic */
  int numbers, numbers1;        /* Numbers diagnostic */
  int save_force;               /* Save input forcing */
  int smooth;                   /* Smoothing flag for topography */
  char smooth_v[MAXSTRLEN];     /* Smoothing flag for other variables */
  char scale_v[MAXSTRLEN];      /* Scaling flag for other variables */
  double slipprm;               /* Slip condition */
  int u1_f;                     /* Flag to ommit u1 momentum terms */
  int u1av_f;                   /* Flag to ommit u1av momentum terms */
  int orbital;                  /* Compute orbital wave speed */
  int waves;                    /* Include wave enhanced bottom friction */
  int rampf;                    /* Variables to ramp flag */
  int totals;                   /* Plots integrated totals */
  double totals_dt;             /* Output interval for totals */
  int ntot;                     /* Number of additional total tracers */
  char **totname;               /* Name of additional totals tracers */
  int roammode;                 /* ROAM configuration version */
  int robust;                   /* Robust scale for ROAM */
  char robusttext[MAXSTRLEN];   /* Description of robust scale */
  int speed;                    /* Speed scale for ROAM */
  int fatal;                    /* Instability check flag */
  int smag_smooth;              /* Smagorinsky smoothing */
  int avhrr;                    /* Include AVHRR SST */
  char avhrr_path[MAXSTRLEN];   /* AVHRR SST data path */
  int ghrsst;                   /* Include GHRSST SST */
  char ghrsst_path[MAXSTRLEN];  /* GHRSST SST data path */
  char alert[MAXSTRLEN];        /* Create alert log */
  char alert_dt[MAXSTRLEN];     /* Time step for alert ts file */
  int eta_f;                    /* Alert action on elevation          */
  int vel2d_f;                  /* Alert action on 2d velocity        */
  int vel3d_f;                  /* Alert action on 3d velocity        */
  int wvel_f;                   /* Alert action on vertical velocity  */
  int tend_f;                   /* Alert action on tendencies         */
  int div2d_f;                  /* Alert action on 2d divergence      */
  int div3d_f;                  /* Alert action on 3d divergence      */
  int cfl_f;                    /* Alert action on CFL                */
  int ts_f;                     /* Alert action on T/S                */
  char alertc[MAXSTRLEN];       /* Alert codes */
  int tide_r;                   /* Active alert tide removal method   */
  int shear_f;                  /* Alert action on velocity shear     */
  int hdiff_f;                  /* Perform enhanced horizontal viscosity */
  double amax;                  /* Maximum advection tendency         */
  double hmax;                  /* Maximum horizontal diffusion tend. */
  double vmax;                  /* Maximum vertical diffusion tend.   */
  double btmax;                 /* Maximum barotropic pressure tend.  */
  double bcmax;                 /* Maximum baroclinic pressure tend.  */
  double cmax;                  /* Maximum coriolis tendency          */
  double detamax;               /* Maximum 2D divergence              */
  double dwmax;                 /* Maximum 3D divergence              */
  double dtmax, dsmax;          /* Maximum T/S rates of change        */
  double smax;                  /* Maximum shear (ms-1)               */
  double emean;                 /* Mean eta mean                      */
  double maxgrad;               /* Maximum bathymetry gradient        */
  int bathyfill;                /* Unstructured land filling method   */
  double bvel;                  /* Velocity for bathymetry check      */
  double lnm;                   /* Level of no motion for steric height */
  int trout;                    /* Transport file output flag */
  int tmode;                    /* Transport mode                     */
  int fillf;                    /* Fill method                        */
  int pssinput;                 /* Point ss input method              */
  int lyear;                    /* Remove leap years */
  int conserve;                 /* Volume conservation                */
  int us_type;                  /* Unstructured type                  */
  int do_closure;               /* Do vertical mixing                 */
  char trans_data[MAXSTRLEN];   /* Transport input filename           */
  char sourcefile[MAXSTRLEN];   /* Transport source grid filename     */
  char trvars[MAXSTRLEN];       /* Transport variables                */
  int parray_inter_botz;        /* Interpolate botz in PARRAY      UR */
  char sequence[MAXSTRLEN];     /* Run sequencing */
  int compatible;               /* Backwards compatible flag */
  int filter;                   /* Filtering options */
  int trfilter;                 /* Tracer filtering options */
  int porusplate;               /* Include porus plate sub-gridscale parameterisation */
  char reef_frac[MAXSTRLEN];    /* Cell area blocked by reef */
  int dbi, dbj, dbk;            /* Sparse coordinate to debug */
  int dbgf;                     /* Debug flag */
  double dbgtime;               /* Time to commence debugging */
  char bdry_file[MAXSTRLEN];    /* Point array boundary file name */
  int do_pt;                    /* Particle tracking flag */
  char ptinname[MAXSTRLEN];     /* Particle input file */
  int gint_errfcn;              /* Generic interface error handling flag */
  int riverflow;                /* Include river flow diagnostic tracer */
  char nprof[MAXSTRLEN];        /* Normalized profile flag */
  /* DATA ASSIM */
  int da;                       /* Data assimilation */
  double da_dt;                 /* Data assimilation time step */
  double da_fcst_dt;            /* Data assimilation forcast time step */
  char restart_name[MAXSTRLEN]; /* Restart file name */
  double restart_dt;            /* Restart file write increment */
  int nland;                    /* Number of cells to redefine as land */
  int *lande1;                  /* e1 list of defined land cells */
  int *lande2;                  /* e2 list of defined land cells */
  char bathystats[MAXSTRLEN];   /* Bathy file for bathymetry statistics */
  int data_infill;              /* Use cascade search on input file data */
  char *da_anom_file;           /* File name for the anomaly fields */
  char *da_anom_states;         /* State names to read from the anomaly fields */
  double da_ls;                 /* DA localisation spread */
  int da_nobs;                  /* Number of Data assimilating observations */
  char **da_obs_names;
  char **da_obs_files;
  char **da_obs_types;
  char **da_obs_locs;
  double *da_obs_dt;
  double *da_obs_errs;
  
  /* Constants */
  char prmname[MAXSTRLEN];      /* Parameter file name */
  char idumpname[MAXSTRLEN];    /* Input dump file name */
  char lenunit[MAXSTRLEN];      /* Length units */
  char codeheader[MAXSTRLEN];
  char parameterheader[MAXSTRLEN];
  char projection[MAXSTRLEN];
  char grid_desc[MAXSTRLEN];
  char grid_name[MAXSTRLEN];
  double g;                     /* Acceleration due to gravity (ms-2) */
  double spec_heat;             /* Specific heat at costant pressure */
  double air_dens;              /* Air density */
  double ambpress;              /* Ambient air pressure */
  double hmin;                  /* Minimum cell thickness (m) */
  double uf;                    /* Background friction velocity */
  char quad_bfc[MAXSTRLEN];     /* Quadratic bottom friction coeff.  */
  double rampstart;             /* Start time of ramp period */
  double rampend;               /* End time of ramp period */
  double rampval;               /* Ramp value for forcing */
  double u1vh;                  /* Horizontal viscosity, e1 direction */
  double u2vh;                  /* Horizontal viscosity, e2 direction */
  double u1kh;                  /* Horizontal diffusivity, e1 direction */
  double u2kh;                  /* Horizontal diffusivity, e2 direction */
  double z0;                    /* Bottom roughness */
  double *z0s;                  /* Spatial bottom roughness */
  double *coriolis;             /* Coriolis parameter */
  double *surface;              /* Initial surface height */

  /* Time */
  double t;                     /* Model time */
  double grid_dt;               /* Time step for the 3D mode (s) */
  double dt;                    /* Time step for the 3D mode (s) */
  double dt2d;                  /* Time step for the 2D mode (s) */
  double dtu1;                  /* Time step for the 3D u1 velocity (s) */
  double dtu2;                  /* Time step for the 3D u2 velocity (s) */
  char timeunit[MAXSTRLEN];     /* Time units */
  char output_tunit[MAXSTRLEN]; /* Output time units */
  int iratio;                   /* Ratio of 3D /2D time-steps */
  double tratio;                /* Ratio of 3D / tracer time-steps */
  int trsplit;                  /* Split advection using tratio */
  /* Transport DT */
  char trans_dt[MAXSTRLEN];

  /* Mixing scheme */
  char mixsc[MAXSTRLEN];        /* Scheme type */
  char s_func[MAXSTRLEN];       /* Stability function */
  double kz0;                   /* Csanady coefficient */
  double kz_alpha;              /* Csanady coefficient */
  double zs;                    /* Mellor-Yamada coefficient */
  double Lmin;                  /* Minimum stratified mixing length */
  double eparam;                /* Tuning parameter for stratification */
  double vz0;                   /* Csanady coefficient */
  double vz_alpha;              /* Csanady coefficient */
  double wave_alpha;            /* Alpha parameter for waves */
  double wave_hf;               /* Scaling factor for significant wave height */
  double wave_b1;               /* b1 parameter for waves */
  double min_tke;               /* k-epsilon minimum TKE */
  double min_diss;              /* k-epsilon minimum dissipation */
  int smooth_VzKz;              /* Shuman smoothing of Vz and Kz */
  int fcf;                      /* Flag for k-w Wilcox (1988)/(1998) models */

  /* Tracer constants and variables */
  int ntr;                      /* Number of tracers */
  int atr;                      /* Number of automatically defined tracers */
  int ntrS;                     /* Number of 2D tracers */
  int atrS;                     /* Number of automatically defined tracers */
  int nsed;                     /* Number of sediment variables */
  int ntre;                     /* Number of 3D edge tracers */
  int ntreS;                    /* Number of 2D edge tracers */
  char tracerdata[MAXSTRLEN];   /* Tracer data */
  double *fill_value_wc;        /* Water column tracer fill value */
  int *advect;                  /* Flag to advect tracers */
  int *diffuse;                 /* Flag to diffuse tracers */
  int ntdif_h;                  /* Number of tracers to horiz. diffuse */
  int *tdif_h;                  /* Array of tracers to horiz. diffuse */
  int ntdif_v;                  /* Number of tracers to vert. diffuse */
  int *tdif_v;                  /* Array of tracers to vert. diffuse */
  double *mintr;                /* Minimum tracer values */
  double *maxtr;                /* Maximum tracer values */
  char **trname;                /* Name of the tracer */
  char **trrlxn;                /* Tracer relaxation filename */
  double *trrlxdt;              /* Tracer relaxation input time */
  char **trrlxtc;               /* Tracer relaxation time constant */
  int rtemp;                    /* Relax temperature flag */
  int rsalt;                    /* Relax salinity flag */
  char **trrest;                /* Tracer reset filename */
  double *trrestdt;             /* Tracer reset input time */
  int *trinc;                   /* 3D Tracer increment flag */
  int *trincS;                  /* 2D Tracer increment flag */
  int *tbdy;                    /* Array of tracers to set OBC's */
  tracer_info_t *trinfo_3d;     /* 3D Tracer info data structure */
  tracer_info_t *trinfo_2d;     /* 2D Tracer info data structure */
  tracer_info_t *trinfo_sed;    /* Sediment Tracer info data structure */

  /* Number of OMP threads to use in transport mode */
#ifdef HAVE_OMP
  int trans_num_omp;
#endif

  /* IO */
  int thIO;                     /* Whether to thread I/O */

  /* timeseries_t file cahcing */
  int tsfile_caching;           /* Flag for tsfile caching */

  /* MOM conversion */
  momgrid_t *momgrid;           /* MOM format grid structure */
  char momfile[MAXSTRLEN];      /* MOM conversion file name           */
  int domom;                    /* Create MOM grid structure */

  /* ROMS conversion */
  romsgrid_t *romsgrid;         /* ROMS format grid structure */
  char romsfile[MAXSTRLEN];     /* ROMS conversion file name           */
  char roms_grid[MAXSTRLEN];    /* ROMS grid file for input */
  double roms_z2s;              /* ROMS Z coord to Sigma num layers factor */
  int doroms;                   /* Create ROMS grid structure */
  int roms_grid_opts;           /* ROMS grid options */

  /* Open boundary structures */
  int nobc;                     /* Number of open boundaries */
  open_bdrys_t **open;          /* Open boundary data structure */
  char orthoweights[MAXSTRLEN]; /* csr tide model orthoweights path */
  char nodal_dir[MAXSTRLEN];    /* Tidal nodal corrections path */
  char tide_con_file[MAXSTRLEN];/* Tidal constituent file */
  char bdrypath[MAXSTRLEN];     /* Path for boundary files */

#if defined(HAVE_WAVE_MODULE)
  int do_wave;                  /* Wave flag */
  double wavedt;                /* wave timestep */
#endif
  char fetch_init[MAXSTRLEN];   /* Name of fetch initialisation file */

#if defined(HAVE_ECOLOGY_MODULE)
  /* Ecology */
  int do_eco;                   /* Ecology model flag */
  double ecodt;                 /* Ecology timestep */
  char eco_vars[MAXSTRLEN];     /* Ecology tracers */
  char eco_defs[MAXSTRLEN];     /* Ecology attribute defaults */
  void *pre_eco;                /* Pointer to ecology pre_built structure */
#endif

#if defined(HAVE_SEDIMENT_MODULE)
  /* Sediments */
  int do_sed;                   /* Sediment flag */
  int do_bm;                    /* Box model transport flag */
  int do_bmphys;                /* Box model physics flag */
  char sed_vars[MAXSTRLEN];     /* Sediment size classes */
  char sed_defs[MAXSTRLEN];     /* Sediment attribute defaults */
#endif
  int sednz;                    /* Number of sediment layers */
  double *gridz_sed;            /* Depths of the sediment layer faces */

  /* Forcing files */
  /* Surface elevation */
  char eta_init[MAXSTRLEN];     /* Name of eta initialisation file */
  int etarlx;                   /* Eta relaxation flag */
  char etarlxn[MAXSTRLEN];      /* Elevation relaxation filename */
  double etarlxdt;              /* Elevation relaxation input time */
  double etarlxtc;              /* Elevation relaxation time constant */
  char etarlxtcs[MAXSTRLEN];    /* Elevation time constant string */

  /* Velocity */
  char vel_init[MAXSTRLEN];     /* Name of velocity initialisation file */
  int velrlx;                   /* Velocity relaxation flag */
  char velrlxn[MAXSTRLEN];      /* Velocity relaxation filename */
  char velrlxtcs[MAXSTRLEN];    /* Velocity time constant string */
  double velrlxdt;              /* Velocity relaxation input time */

  /* Wind */
  char wind[MAXSTRLEN];         /* Name of wind input file */
  double wind_dt;               /* Wind input time-step */
  double wind_scale;            /* Wind stress scaling factor */
  double dlv0;                  /* Drag law v0 */
  double dlv1;                  /* Drag law v1 */
  double dlc0;                  /* Drag law c0 */
  double dlc1;                  /* Drag law c1 */
  int wind_type;                /* Flag to read wind speed or stress */
  int stress_fn;                /* Wind stress function */
  int neutral;                  /* Use neutral drag */
  int nstorm;                   /* Number of storm systems */
  double storm_dt;              /* Storm input time-step */
  /* Heatflux */
  int heatflux;                 /* Type of heatflux specification */
  int saltflux;                 /* Type of saltflux specification */
  int water_type;               /* Water type */
  double hf_ramp;               /* Heat flux ramp time */
  char swr_attn[MAXSTRLEN];     /* Short wave attenuation (for red) */
  char swr_attn1[MAXSTRLEN];    /* Short wave attenuation for blue-green */
  int swr_type;                 /* Type of attenuation */
  char swr_tran[MAXSTRLEN];     /* Short wave surface transmission */
  char swr_babs[MAXSTRLEN];     /* Short wave bottom absorption */
  double zref;                  /* Reference height */
  double albedo;                /* Albedo */
  double albedo_l;              /* Albedo for light */
  int hfadf;                    /* Advection flag */
  int bulkf;                    /* Bulk scheme for heatflux */
  char hftemp[MAXSTRLEN];       /* Heat flux temperature filename */
  double hftemp_dt;             /* Heat flux temperature input time-step */
  char hf[MAXSTRLEN];           /* Name of heat flux file */
  double hf_dt;                 /* Heat flux input time-step */
  double hftc;                  /* Time constant for heatflux temp */
  char airtemp[MAXSTRLEN];      /* Name of air temperature file */
  double airtemp_dt;            /* Air temperature input time-step */
  char evap[MAXSTRLEN];         /* Name of evaporation file */
  double evap_dt;               /* Evaporation input time-step */
  char cloud[MAXSTRLEN];        /* Name of cloud cover file */
  double cloud_dt;              /* Cloud cover input time-step */
  char patm[MAXSTRLEN];         /* Name of atmospheric pressure file */
  double patm_dt;               /* Atmospheric pressure input time-step */
  char precip[MAXSTRLEN];       /* Name of precipitation file */
  double precip_dt;             /* Precipitation input time-step */
  char rh[MAXSTRLEN];           /* Name of relative humidity file */
  double rh_dt;                 /* Relative humidity input time-step */
  char swr[MAXSTRLEN];          /* Name of short wave radiation file */
  double swr_dt;                /* Short wave radiation input time-step */
  char webf[MAXSTRLEN];         /* Name of wave enhanced friction file */
  double webf_dt;               /* Wave friction input time-step */
  char webf_interp[MAXSTRLEN];  /* Interpolation method for webf vars */
  char wetb[MAXSTRLEN];         /* Name of wet bulb temperature file */
  double wetb_dt;               /* Wet bulb temperature input time-step */
  char light[MAXSTRLEN];        /* Name of light radiation file */
  double light_dt;              /* Light input time-step */
  char regulate[MAXSTRLEN];     /* Name of user regulation file */
  double regulate_dt;           /* User regulation time-step */

  /* ROAM */
  char odata[MAXSTRLEN];        /* OFAM forcing data */
  char rdata[MAXSTRLEN];        /* RAMS forcing data */
  char idata[MAXSTRLEN];        /* Initialisation data */
  char tdata[MAXSTRLEN];        /* Temp initialisation data */
  char sdata[MAXSTRLEN];        /* Salt initialisation data */
  char edata[MAXSTRLEN];        /* Eta initialisation data */
  char vdata[MAXSTRLEN];        /* Velocity data */

  /* Explicit maps */
  int exmapf;                   /* Explicit map flag */
  int nemx, *ktx, *kbx;
  int *emisx, *emidx;
  int *emjsx, *emjdx;
  int nemy, *kty, *kby;
  int *emisy, *emidy;
  int *emjsy, *emjdy;

  /* Process exclude */
  int prex, *prxi, *prxj, *prxf;

  /* RIVER loads directory for RECOM */
  char rivldir[MAXSTRLEN];

  int ndf;                      /* Number of dump files */
  char **d_name;
  char **d_filetype;
  char **d_tstart;
  char **d_tinc;
  char **d_tstop;
  char **d_vars;
  char **d_sync;
  char **d_tunit;
} parameters_t;


/*------------------------------------------------------------------*/
/* Master data structure                                            */
/*------------------------------------------------------------------*/
struct master {

  /* Parameter file */
  geometry_t *geom;             /* Global sparse geometry structure */
  geometry_t *sgrid;            /* Unstructured geometry structure */
  dump_data_t *dumpdata;        /* Output dump data structure */
  parameters_t *params;         /* Input parameter data structure */
  FILE *prmfd;                  /* Handle to parameter file */
  char prmname[MAXSTRLEN];      /* Parameter file name */
  char idumpname[MAXSTRLEN];    /* Input dump file name */
  char opath[MAXSTRLEN];        /* Output path for files */

  /* Flags */
  int nwindows;                 /* Number of windows */
  int win_reset;                /* Number of steps to reset window loads */
  int show_win;                 /* Create plot of windowing */
  int runmode;                  /* Type of run (3D, 2D, transport) */
  int trasc;                    /* Advection scheme type flag (tracers) */
  int momsc;                    /* Advection scheme type flag (velocity) */
  int momsc2d;                  /* Advection scheme type flag (velocity) */
  int otras;                    /* Order of the advection scheme */
  int osl;                      /* Order of semi-Lagrange */
  int ultimate;                 /* ULTIMATE filter flag */
  int diff_scale;               /* Horizontal diffusion scaling method */
  int visc_method;              /* Horizontal diffusion method */
  int stab;                     /* Stability compensation method */
  int thin_merge;               /* Thin layer merging flag */
  int sigma;                    /* Sigma coordinate flag */
  int nonlinear;                /* Non-linearity flag */
  int calc_dens;                /* Calculate density flag */
  int mode2d;                   /* Run in 2D mode */
  int cfl;                      /* CFL time-step diagnostic */
  double cfl_dt;                /* Time for active cfl specification */
  double lnm;                   /* Level of no motion for steric height */
  int mixlayer;                 /* Mixed layer depth diagnostic */
  int show_layers;              /* Create plot of sigma layers */
  int trflux;                   /* Tracer number for flux calculation */
  int trfd1;                    /* Face #1 for flux calculation */
  int trfd2;                    /* Face #2 for flux calculation */
  int trperc;                   /* Tracer no. for percentile calculation */
  char trpercr[MAXSTRLEN];      /* Region name for percentile calc. */
  int trflsh;                   /* Flushing tracer flag */
  char trage[MAXSTRLEN];        /* Age tracer flag */
  int tendf;                    /* Momentum tendency flag */
  int trtend;                   /* Tracer tendency flag */
  int means;                    /* Mean velocity diagnostic */
  double means_dt;              /* Mean velocity averaging interval */
  double means_next;            /* Next time for zeroing means */
  double means_os;              /* Offset for restarts */
  double *meanc;                /* Iteration counter for means */
  double *meancs;               /* Seasonal mean iteration counter */
  int means_tra;                /* Mean tracer number */
  double *odeta;                /* detadt at previous timestep */
  int vorticity;                /* Vorticity diagnostic */
  int numbers, numbers1;        /* Numbers diagnostic */
  int save_force;               /* Save input forcing */
  double slip;                  /* Slip condition */
  double flt;                   /* Flushing time */
  double flts;                  /* Flushing time start */
  double imass;                 /* Initial mass in flushing region */
  short *mask;                  /* Flushing mask */
  int u1_f;                     /* Flag to ommit u1 momentum terms */
  int u1av_f;                   /* Flag to ommit u1av momentum terms */
  int orbital;                  /* Compute orbital wave speed */
  int waves;                    /* Include wave enhanced bottom friction */
  int rampf;                    /* Variables to ramp flag */
  int roammode;                 /* ROAM configuration version */
  int robust;                   /* Robust scale for ROAM */
  int fatal;                    /* Instability check flag */
  int smag_smooth;              /* Smagorinsky smoothing */
  int trout;                    /* Transport file output flag */
  int tmode;                    /* Transport mode                     */
  int fillf;                    /* Fill method                        */
  int pssinput;                 /* Point ss input method              */
  int lyear;                    /* Remove leap years */
  int conserve;                 /* Volume conservation                */
  int sh_f;                     /* Data input type for specific humidity */
  int compatible;               /* Backwards compatible flag */
  int filter;                   /* Filtering options */
  int trfilter;                 /* Tracer filtering options */
  int porusplate;               /* Include porus plate sub-gridscale parameterisation */
  char reef_frac[MAXSTRLEN];    /* Cell area blocked by reef */
  double *reefe1;               /* Pointer to e1 reef fraction tracer */
  double *reefe2;               /* Pointer to e2 reef fraction tracer */
  int dhwf;                     /* Degree heating diagnostic */
  int dbc;                      /* Sparse coordinate to debug */
  int dbj;                      /* Edge direction to debug */
  int dbgf;                     /* Debug flag */
  double dbgtime;               /* Time to commence debugging */
  int da;                       /* Data assimilation */
  double da_dt;                 /* Data assimilation time step */
  double da_fcst_dt;            /* Data assimilation forcast time step */
  char restart_name[MAXSTRLEN]; /* Restart file name */
  double restart_dt;            /* Restart file write increment */
  int crf;                      /* Crash recovery flag */
  int regf;                     /* Run regulation flag */
  int data_infill;              /* Use cascade search on input file data */
  char runerror[MAXSTRLEN];     /* Run error message */
  int decorr;                   /* Decorrelation length scale diagnostic */
  double decs;                  /* Decorrelation length scaling */
  int decf;                     /* Decorrelation length flag */
  int decn;                     /* Decorrelation length tracer number */
  int gint_errfcn;              /* Generic interface error handling flag */
  int ntrmap_s2m_3d;            /* Number of tracers to transfer from
				   slave to master. This is ntr minus
				   some number */
  int *trmap_s2m_3d;            /* The actual map for the above */

  /* Constants */
  double g;                     /* Acceleration due to gravity (ms-2) */
  double ambpress;              /* Ambient atmospheric pressure */
  double hmin;                  /* Minimum cell thickness (m) */
  double hmean1;                /* Mean grid spacing between centres (m) */
  double hmean2;                /* Mean grid spacing of edges (m) */
  double amean;                 /* Mean area (m^2) */
  double minres;                /* Minimum resolution (m) */
  double maxres;                /* Maximum resolution (m) */
  int minrese;                  /* Location of minimum resolution */
  int maxrese;                  /* Location of maximum resolution */
  double u1vh0;                 /* Horizontal e1 viscosity (m2s-1) */
  double u2vh0;                 /* Horizontal e2 viscosity (m2s-1) */
  double u1kh0;                 /* Horizontal e1 diffusivity (m2s-1) */
  double u2kh0;                 /* Horizontal e2 diffusivity (m2s-1) */
  double uf;                    /* Background friction velocity */
  double quad_bfc;              /* Quadratic bottom friction coeff.  */
  double mclk;                  /* Mean CPU time */
  double iclk;                  /* Total CPU time for this time-step */
  double mrtr;                  /* Mean run time ratio */
  double *wclk;                 /* CPU time for each window */
  double etamax;                /* Maximum elevation (m) */
  double velmax;                /* Maximum horizontal velocity (ms-1) */
  double velmax2d;              /* Maximum horizontal 2D velocity (ms-1) */
  double etadiff;               /* Maximum difference in mean sea level (m) */
  double wmax;                  /* Maximum vertical velocity (ms-1) */
  double rampstart;             /* Start time of ramp period */
  double rampend;               /* End time of ramp period */
  double rampval;               /* Ramp value for forcing */
  double *Cd;                   /* Drag coefficient */
  double *coriolis;             /* Coriolis parameter */
  double *u1c1;                 /* e1 advection term constant */
  double *u1c3;                 /* e1 advection metric constant */
  double *u1c4;                 /* e1 advection metric constant */
  double *u1c5;                 /* e1 Coriolis constant */
  double *u1c6;                 /* e1 Coriolis constant */
  double *u2c1;                 /* e2 advection term constant */
  double *u2c3;                 /* e2 advection metric constant */
  double *u2c4;                 /* e2 advection metric constant */
  double *u2c5;                 /* e2 Coriolis constant */
  double *u2c6;                 /* e2 Coriolis constant */
  double *one;                  /* 2D work array set to 1.0 */

  /* Diagnostics */
  double mcfl2d;                /* Mean minimum 2D cfl time-step */
  double mcfl3d;                /* Mean minimum 3D cfl time-step */
  int cflc;                     /* Location of min 3D cfl for ACTIVE3D */
  int totals;                   /* Plots integrated totals */
  double totals_dt;             /* Output interval for totals */
  double tvol;                  /* Total volume */
  double tmass;                 /* Total mass */
  int ntot;                     /* Number of additional total tracers */
  char **totname;               /* Name of additional totals tracers */
  double *trtot;                /* Total tracers */
  double *cfl2d;                /* 2D cfl time-step */
  double *cfl3d;                /* 3D cfl time-step */
  double *cour;                 /* Courant stability */
  double *lips;                 /* Lipshitz stability */
  double *ahsb;                 /* Horizontal diffusion stability */
  double *mixl;                 /* Pointer to mixed layer depth */
  double *fluxe1;               /* Advective flux in the e1 direction */
  double *fluxe2;               /* Advective flux in the e1 direction */
  double *fluxw;                /* Vertical advective flux */
  double *fluxkz;               /* Vertical diffusive flux */
  double *u1m;                  /* Pointer to u mean velocity */
  double *u2m;                  /* Pointer to v mean velocity */
  double *ume;                  /* Mean 3D velocity */
  double *wm;                   /* Pointer to w mean velocity */
  double *u1vm;                 /* Pointer to u1 mean volume flux */
  double *u2vm;                 /* Pointer to u2 mean volume flux */
  double *Kzm;                  /* Pointer to mean Kz */
  double *etam;                 /* Pointer to mean elevation */
  double *u1am;                 /* Pointer to uav mean velocity */
  double *u2am;                 /* Pointer to vav mean velocity */
  double *uame;                 /* Mean 2D velocity */
  double *w1m;                  /* Pointer to eastward mean wind */
  double *w2m;                  /* Pointer to northward mean wind */
  double *tempm;                /* Pointer to mean temperature */
  double *saltm;                /* Pointer to mean salinity */
  double *tram;                 /* Pointer to mean tracer */
  double *perc;                 /* Pointer to percentile diagnostic */
  double *fltr;                 /* Pointer to flushing tracer */
  double *agetr;                /* Pointer to age tracer */
  double *steric;               /* Steric height diagnostic */
  double *rv;                   /* Relative vorticity diagnostic */
  double *av;                   /* Absolute vorticity diagnostic */
  double *pv;                   /* Potential vorticity diagnostic */
  double *rv_drvdt;             /* Rate of change rv tendency */
  double *rv_nonlin;            /* Nonlinear rv tendency (adv + diff) */
  double *rv_beta;              /* Planetary vorticity rv tendency */
  double *rv_strch;             /* Topographic stretching rv tendency */
  double *rv_jebar;             /* JEBAR rv tendency */
  double *rv_wsc;               /* Wind stress curl rv tendency */
  double *rv_bsc;               /* Bottom stress curl rv tendency */
  double *riverflow;            /* River flow diagnostic */
  double *riverdepth;           /* River depth diagnostic */
  double *riversalt;            /* River ghost salinity diagnostic */
  double *stream;               /* 2D streamfunction */
  double *brunt;                /* Brunt Vaisala (buoyancy) frequency (s-1) */
  double *int_wave;             /* Internal wave speed (ms-1) */
  double *rich_gr;              /* Gradient Richardson number */
  double *rich_fl;              /* Flux Richardson number */
  double *reynolds;             /* Reynolds number */
  double *froude;               /* Froude number */
  double *sep;                  /* Surface Ekman pumping */
  double *bep;                  /* Bottom Ekman pumping */
  double *tfront;               /* Simpson-Hunter tidal front */
  double *sigma_t;              /* Sigma_t */
  double *rossby_in;            /* Internal Rossby radius (m) */
  double *rossby_ex;            /* External Rossby radius (m) */
  double *speed_2d;             /* 2D current speed */
  double *speed_3d;             /* 3D current speed */
  double *speed_sq;             /* Bottom speed squared */
  double *energy;               /* Total energy (Jm-2) */
  double *kenergy;              /* Kinetic energy (Jm-2) */
  double *obc_phase;            /* OBC phase speed (m/s) */
  double *nprof;                /* Normalized profile */
  int nprofn;                   /* Tracer to normalize */
  double *sound;                /* Speed of sound */
  double *schan;                /* Sound channel depth */
  double *sonic;                /* Sonic layer depth */
  double *wetcell;              /* Wet cell diagnostic */
  double *slope_x;              /* Surface slope (e1 direction) diagnostic */
  double *slope_y;              /* Surface slope (e2 direction) diagnostic */
  double *surfz;                /* Surface layer diagnostic */
  double *wind_Cd;              /* Wind drag */
  double *shear_v;              /* Vertical shear (s-1) */
  double *b_prod;               /* Buoyancy production (m2s-2) */
  double *s_prod;               /* Shear production (m2s-2) */
  double *otemp;                /* OFAM temperature */
  double *osalt;                /* OFAM salinity */
  double *rtemp;                /* Relaxation temperature */
  double *rsalt;                /* Relaxation salinity */
  double *temp_tc;              /* Relaxation temperature time constant */
  double *salt_tc;              /* Relaxation salinity time constant */
  double *eta_tc;               /* Relaxation eta time constant */
  double *eta_inc;              /* Relaxation eta increment */
  double *avhrr;                /* AVHRR SST */
  double *ghrsst;               /* GHRSST SST */
  double *shwin;                /* Window partitioning */
  double *alert_a;              /* Actual alert diagnostic */
  double *alert_c;              /* Cumulative alert diagnostic */
  double *u1vhin;               /* Initial e1 horizontal viscosity */
  double *u2vhin;               /* Initial e2 horizontal viscosity */
  double *layth;                /* Layer thickness for sigma */
  double *u1_rad;               /* u1 velocity radiation stress tend. */
  double *u2_rad;               /* u2 velocity radiation stress tend. */
  double *tau_be1;              /* Bottom stress in e1 direction */
  double *tau_be2;              /* Bottom stress in e2 direction */
  double *tau_bm;               /* Bottom stress magnitude */
  double *dum1, *dum2, *dum3;   /* Dummy arrays */
  double *vcorr, *acorr;        /* Local transport fill corrections */
  double *Vi;                   /* Volume error */
  double *unit;                 /* Unit tracer */
  double *glider;               /* Glider density */
  double *u1vhc;                /* Cell centered horizontal viscosity */
  double *regionid;             /* Region ids */
  double *regres;               /* Region residence time */
  double *sederr;               /* Error percentage map for sediments */
  double *ecoerr;               /* Error percentage map for ecology */
  double *decv;                 /* Decorrelation length scale variable */
  double *decv1;                /* Decorrelation length scale */
  double *dhd;                  /* Degree heating day */
  double *dhwc;                 /* Degree heating week climatology */
  double *dhw;                  /* Degree heating week */
  double *cellres;              /* Mean cell resolution */
  char bathystats[MAXSTRLEN];   /* Bathy file for bathymetry statistics */
  double *bathy_range_max;
  double *bathy_range_min;
  double *bathy_grad_max;
  double *bathy_grad_min;

  /* Vertical grid geometry */
  double *dz;                   /* Layer thickness at cell center (m) */
  double *dzu1;                 /* Layer thickness at e1 face (m) */
  double *dzu2;                 /* Layer thickness at e2 face (m) */
  double *depth_e1;             /* Total depth at e1 face (m) */
  double *depth_e2;             /* Total depth at e2 face (m) */
  double *topz;                 /* Height of the surface (m) */
  int *topk;                    /* Vertical index of top layer */
  int *botk;                    /* Vertical index of bottom layer */

  /* SIGMA variables */
  double *Ds;                   /* SIGMA scaling factor (m) */
  double *Hs;                   /* SIGMA water depth (m) */
  double *Hn1;                  /* SIGMA scaling factor (m) */
  double *Hn2;                  /* SIGMA scaling factor (m) */

  /* e1 velocity variables */
  int vinit;                    /* Velocity initialisation flag */
  double *u1;                   /* Normal velocity to edge (ms-1) */
  double *u;                    /* Eastward velocity (ms-1) */
  double *nu1;                  /* Updated 3D e1 velocity (ms-1) */
  double *u1bot;                /* Bottom u1 velocity (ms-1) */
  double *u1flux3d;             /* Flux in the e1 direction (m3s-1) */
  double *u1adv;                /* e1 dispersion terms */
  double *u1inter;              /* Vertically averaged u1 terms */
  double *u1av;                 /* Normal 2D velocity to edge (ms-1) */
  double *uav;                  /* Eastward 2D velocity (ms-1) */
  double *nu1av;                /* Updated e1 2D velocity (ms-1) */
  double *u1flux;               /* Volume flux through e1 2D faces (m3) */

  /* e2 velocity variables */
  double *u2;                   /* Tangential velocity to edge (ms-1) */
  double *v;                    /* Northward velocity (ms-1) */
  double *nu2;                  /* Updated 3D e1 velocity (ms-1) */
  double *u2bot;                /* Bottom u2 velocity (ms-1) */
  double *u2flux3d;             /* Flux in the e2 direction (m3s-1) */
  double *u2adv;                /* e2 dispersion terms */
  double *u2inter;              /* Vertically averaged u1 terms */
  double *u2av;                 /* Tangential 2D velocity to edge (ms-1) */
  double *vav;                  /* Northward 2D velocity (ms-1) */
  double *nu2av;                /* Updated e2 2D velocity (ms-1) */
  double *u2flux;               /* Volume flux through e2 2D faces (m3) */

  /* Vertical velocity variables */
  double *w;                    /* Vertical velocity (ms-1) */
  double *wtop;                 /* Surface vertical velocity (ms-1) */
  double *wbot;                 /* Bottom vertical velocity (ms-1) */

  double *etab;                 /* Height of the free surface at t-1 (m) */
  double *u1avb;                /* 2D u1 velocity at previous timestep */
  double *u2avb;                /* 2D u2 velocity at previous timestep */
  double *u1b;                  /* 3D u1 velocity at previous timestep */
  double *u2b;                  /* 3D u2 velocity at previous timestep */

  /* Time stepping constants and variables */
  double t;                     /* Simulation time */
  double tstart;                /* Start time of the simulation */
  double grid_dt;               /* Time step for the 3D mode (s) */
  double ogrid_dt;              /* Original step for the 3D mode (s) */
  double dt;                    /* Time step for the 3D mode (s) */
  double dt2d;                  /* Time step for the 2D mode (s) */
  double dtu1;                  /* Time step for the 3D u1 velocity (s) */
  double dtu2;                  /* Time step for the 3D u2 velocity (s) */
  double dttr;                  /* Time step for the tracers (s) */
  double t3d;                   /* Time at the end of 3D step */
  double dtf;                   /* Forward part of leapfrog dt */
  double dtb;                   /* Backward part of leapfrog dt */
  char timeunit[MAXSTRLEN];     /* Time units */
  char output_tunit[MAXSTRLEN]; /* Output time units */
  int iratio;                   /* Ratio of 3D /2D time-steps */
  double tratio;                /* Ratio of 3D / tracer time-steps */
  int trsplit;                  /* Split advection using tratio */
  double days;                  /* Time in days (t/86400) */

  /* Tracer constants and variables */
  int ntr;                      /* Number of 3D tracers */
  int atr;                      /* Number of automatically defined tracers */
  int ntrS;                     /* Number of 2D tracers */
  int atrS;                     /* Number of automatically defined tracers */
  int nsed;                     /* Number of sediment variables */
  int ntre;                     /* Number of 3D edge tracers */
  int ntreS;                    /* Number of 2D edge tracers */
  char tracerdata[MAXSTRLEN];   /* Tracer data */
  double **tr_wc;               /* 3D Tracers */
  double **tr_wcS;              /* 2D Tracers */
  double ***tr_sed;             /* Sediment tracer array */
  int ntbdy;                    /* Number of tracers to set OBC's */
  double *fill_value_wc;        /* Water column tracer fill value */
  int *advect;                  /* Flag to advect tracers */
  int *diffuse;                 /* Flag to diffuse tracers */
  int ntdif_h;                  /* Number of tracers to horiz. diffuse */
  int *tdif_h;                  /* Array of tracers to horiz. diffuse */
  int ntdif_v;                  /* Number of tracers to vert. diffuse */
  int *tdif_v;                  /* Array of tracers to vert. diffuse */
  int ntm_2d;                   /* Number of 2D mean tracers */
  int *tm_2d;                   /* Array of 2D mean tracers */
  int ntm_3d;                   /* Number of 2D mean tracers */
  int *tm_3d;                   /* Array of 2D mean tracers */
  int *relax;                   /* Tracers to undergo relaxing */
  int nrlx;                     /* Number of tracers to relax */
  int *reset;                   /* Tracers to undergo resetting */
  int nres;                     /* Number of tracers to reset */
  int *reset2d;                 /* 2D Tracers to undergo resetting */
  int nres2d;                   /* Number of 2D tracers to reset */
  double *mintr;                /* Minimum tracer values */
  double *maxtr;                /* Maximum tracer values */
  char **trname;                /* Name of the tracer */
  int *tbdy;                    /* Array of tracers to set OBC's */
  tracer_info_t *trinfo_3d;     /* 3D Tracer info data structure */
  tracer_info_t *trinfo_2d;     /* 2D Tracer info data structure */
  tracer_info_t *trinfo_sed;    /* Sediment Tracer info data structure */
  double *sal;                  /* Pointer to salinity tracer */
  double *temp;                 /* Pointer to temperature tracer */
  int tno;                      /* Tracer number for temperature */
  int sno;                      /* Tracer number for salinity */

  /* Mixing constants and variables */
  double *Kz;                   /* Vertical eddy diffusivity (m2s-1) */
  double *Vz;                   /* Vertical eddy viscosity (m2s-1) */
  double vz0;                   /* Background vertical viscosity (m2s-2) */
  double kz0;                   /* Background vertical diffusivity */
  double min_tke;               /* k-epsilon minimum TKE */
  double min_diss;              /* k-epsilon minimum dissipation */
  double Lmin;                  /* Minimum stratified mixing length */
  double eparam;                /* Tuning parameter for stratification */
  double kz_alpha;              /* Csanady coefficient */
  double vz_alpha;              /* Csanady coefficient */
  double wave_alpha;            /* Alpha parameter for waves */
  double wave_hf;               /* Scaling factor for significant wave height */
  double wave_b1;               /* b1 parameter for waves */
  double *q;                    /* Turbulent kinetic energy */
  double *tke;                  /* Turbulent kinetic energy */
  double *diss;                 /* Dissipation */
  double *omega;                /* Turbulent frequency */
  double *Kq;                   /* Turbulence intensity diffusivity */
  double *L;                    /* Turbulence length scale */
  double *Q2;                   /* Turbulence intensity */
  double *Q2L;                  /* Turbulence intensity length scale */
  double zs;                    /* Surface roughness length scale (m) */
  double *z0;                   /* Bottom roughness length scale (m) */
  int smooth_VzKz;              /* Shuman smoothing of Vz and Kz */
  int fcf;                      /* Flag for k-w Wilcox (1988)/(1998) models */
  int do_closure;               /* Do vertical mixing (transport mode) */

  /* Surface elevation variables */
  double *eta;                  /* Height of the free surface at t+1 (m) */
  double *waterss2d;            /* Volume source / sink term */
  double *waterss;              /* Volume source / sink term */
  double *detadt;               /* Rate of change of elevation (ms-1) */
  double *wdiff2d;              /* Surface slope term */
  double *oldeta;               /* Height of the free surface at t (m) */
  /*double *eta_rlx;*/              /* Relaxation  source/sink */
  char etarlxn[MAXSTRLEN];      /* Elevation relaxation filename */
  double etarlxdt;              /* Elevation relaxation input time */
  double etarlxtc;              /* Elevation relaxation time constant */
  char etarlxtcs[MAXSTRLEN];    /* Elevation time constant string */

  /* Density variables */
  double *dens;                 /* Density (kgm-3) */
  double *dens_0;               /* Density at zero pressure (kgm-3) */
  double *topdensu1;            /* Surface density on e1 faces */
  double *topdensu2;            /* Surface density on e2 faces */
  double *densavu1;             /* Mean density on e1 faces */
  double *densavu2;             /* Mean density on e2 faces */

  /* Atmospheric variables */
  double *wind1;                /* Wind stress in the x direction (Nm-2) */
  double *wind2;                /* Wind stress in the y direction (Nm-2) */
  double wind_dt;               /* Wind input time-step */
  int nstorm;                   /* Number of storm systems */
  double storm_dt;              /* Storm input time-step */
  double *swind1;               /* Wind stress in the x direction (Nm-2) */
  double *swind2;               /* Wind stress in the y direction (Nm-2) */
  double *windspeed;            /* Wind speed (ms-1) */
  double *winddir;              /* Wind direction (deg T) */
  double *patm;                 /* Atmospheric pressure (HPa) */
  double *precip;               /* Precipitation (mm day-1) */
  double *evap;                 /* Evaporation (mm day-1) */
  double *airtemp;              /* Air temperature (deg C) */
  double *rh;                   /* Relative humidity (%) */
  double *cloud;                /* Cloud amount (oktas) */
  double *swr;                  /* Short wave radiation (Wm-2) */
  double *light;                /* Daily mean short wave radiation */
  double *wetb;                 /* Wet bulb temperature (deg C) */
  double *hftemp;               /* Sea surface temperature (deg C) */

  /* Radiation albedo and attenuation coefficient */
  double albedo;                /* Albedo for swr */
  double albedo_l;              /* Albedo for light */
  int swr_type;                 /* Type of attenuation */
  double *swr_attn;             /* Short wave attenuation (for red) */
  double *swr_attn1;            /* Short wave attenuation for blue-green */

  /* Heat flux variables */
  int heatflux;                 /* Type of heatflux specification */
  int water_type;               /* Water type */
  double hf_ramp;               /* Heat flux ramp time */
  double *heatf;                /* Surface heat flux (Wm-2) */
  double *swr_tran;             /* Short wave surface transmission */
  double *swr_babs;             /* Short wave bottom absorption */
  double zref;                  /* Reference height */
  int hfadf;                    /* Advection flag */
  int bulkf;                    /* Bulk scheme for heatflux */
  double hftc;                  /* Time constant for heatflux temp */
  double **fetch;               /* Fetch (m) */
  double *nhfd;                 /* Net heat flux diagnostic */
  double *swrd;                 /* Short wave radiation diagnostic */
  double *lwrd;                 /* Long wave radiation diagnostic */
  double *shfd;                 /* Sensible heat flux diagnostic */
  double *lhfd;                 /* Latent heat flux diagnostic */
  double *lwro;                 /* Long wave output radiation */
  int lwrn;                     /* Tracer number for lwr */
  int lhfn;                     /* Tracer number for lhf */
  int shfn;                     /* Tracer number for shf */

  /* Salt flux variables */
  int saltflux;                 /* Type of saltflux specification */
  double *nsfd;                 /* Net salt flux diagnostic */
  int evapn;                    /* Tracer number for evaporation */
  int precipn;                  /* Tracer number for precipitation */

  /* Horizontal diffusion variables */
  double *u1vh;                 /* Horizontal e1 viscosity (m2s-1) */
  double *u2vh;                 /* Horizontal e2 viscosity (m2s-1) */
  double *u1kh;                 /* Horizontal e1 diffusivity (m2s-1) */
  double *u2kh;                 /* Horizontal e2 diffusivity (m2s-1) */
  double *sdc;                  /* Smagorinsky horizontal diffusion */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  int smagcode;                 /* Variables using Smagorinsky */
  double smagorinsky;           /* Smagorinsky horizontal diffusion flag */
  double sue1, sue2;            /* Momentum Smagorinsky coefficients */
  double kue1, kue2;            /* Tracer Smagorinsky coefficients */
  double bsue1, bsue2;          /* Base momentum mixing */
  double bkue1, bkue2;          /* Base tracer mixing */

  /* Relaxation */
  int etarlx;                   /* Eta relaxation flag */
  relax_info_t *eta_rlx;        /* Elevation relaxation */
  int velrlx;                   /* Velocity relaxation flag */
  relax_info_t *vel_rlx;        /* Velocity relaxation */

  /* Unstructured */
  double *circ;                 /* Circulation */
  double *kec;                  /* Kinetic energy at cell centres */
  double *div;                  /* Divergence at cell centres */
  double *rvor;                 /* Relative vorticity */
  double *fv;                   /* Coriolis parameter at vertices */
  double *nrvor;                /* Normalized relative vorticity */
  double *npvor;                /* Normalized potential vorticity */
  double *nrvore;               /* Normalized relative vorticity on edges */
  double *npvore;               /* Normalized potential vorticity on edges */
  double *nrvorc;               /* Normalized relative vorticity on centres */

  /* particle variables */
  int do_pt;                    /* Particle tracking flag */
  long ptn;                     /* Total number of particles */
  particle_t *ptp;              /* Array of particles */
  char ptinname[MAXLINELEN];    /* particle input file name */
  int ptinrec;                  /* particle input record */
  char ptoutname[MAXLINELEN];   /* particle output name */
  double ptstart;               /* particle output start time */
  double ptend;                 /* particle output stop time */
  double ptstep;                /* particle step interval */
  double ptoutinc;              /* particle output step interval */
  double ptreset;               /* particle reset interval */
  int ptfid;
  int ptrec;
  double ptout_t;
  double ptreset_t;
  double ptnext_t;
  double pt_kh;                 /* Horizontal diffusion coefficient */
  double pt_kz_mult;            /* vert diff coeff multiplier */
  double *ptconc;               /* particle concentrations */
  int pt_stickybdry;            /* Is the boundary sticky ? */
  int pt_dumpf;                 /* Flag of particle attributes to dump */
  int pt_fvel;                  /* Flag to invoke particle settling */
  double pt_agelim;             /* Limit for stretching age to short int */
  double pt_sizelim;            /* Limit for stretching size to short int */
  double mage;                  /* Average age of particles */
  double magec;                 /* Counter for mean age */
  int *phist;                   /* Histogram array */
  int shist;                    /* Size of phist */
  double pt_mass;               /* particle mass */
  int pt_nsource;               /* Number of particle sources */
  double *pt_rate;              /* Rate of release from continous source */
  short *pt_colour;             /* Colour of the particle */
  double *pt_size;              /* Size of the particle */
  double *pt_decay;             /* Growth/decay rate of the particle */
  double *pt_svel;              /* Settling velocity of the particle */
  double *pt_sper;              /* Settling period of the particle */
  double *pt_dens;              /* Density of the particle */
  int *svel_type;               /* Type of particle settling */
  pt_ts_t **wvel_s;             /* Vertical settling for sources */
  pt_ts_t *wvel_i;              /* Vertical settling for initial release */
  pt_ts_t *uvel_i;              /* Horizontal e1 swimming for initial release */
  pt_ts_t *vvel_i;              /* Horizontal e2 swimming for initial release */
  pt_ts_t *mort_i;              /* Mortality initial release */
  double *pt_x1;                /* Source line location */
  double *pt_x2;
  double *pt_y1;
  double *pt_y2;
  double *pt_z1;
  double *pt_z2;
  double *pt_accum;
  short *pt_sm;                 /* Particle source map */
  sched_event_t *ts_parts_events; /* Time scheduler for particles */
  float *maxdiffw;              /* Maximum vert diffusion */
  profile_t **u1prof;
  profile_t **u2prof;
  short *ptmsk;                 /* Particle age mask */

  /* Sources and sinks */
  pss_t *pss;
  int npss;
  double *wflux;
  double **tflux;

  /* timeseries_t file cahcing */
  int tsfile_caching;           /* flag for tsfile caching */
  int ntscached;                /* Number of cached ts files.  */
  timeseries_t **tscache;       /* Cached tsfiles */

  /* Explicit maps */
  int exmapf;                   /* Explicit map flag */

  /* Waves */
  double *ustrcw;               /* Wave-current friction velocity at
                                   bottom */
  double *wave_ub;              /* Near bottom wave orbital velocity
                                   [ms-1] */
  double *wave_period;          /* Surface wave period [s] */
  double *wave_dir;             /* Surface wave direction [radian] */
  double *wave_amp;             /* Surface wave amplitude [m] */
  double *wave_Sxy;             /* e1 tangential radiation stress */
  double *wave_Syx;             /* e2 tangential radiation stress */
  double *wave_Fx;              /* e1 wave-induced force */
  double *wave_Fy;              /* e2 wave-induced force */
  double *wave_Cd;              /* Wave enhanced bottom drag */
  double *wave_ste1;            /* Stokes surface velocity, e1 */
  double *wave_ste2;            /* Stokes surface velocity, e2 */
  double *tau_w1;               /* Wave supported wind e1 stress */
  double *tau_w2;               /* Wave supported wind e2 stress */
  double *tau_diss1;            /* Wave to ocean e1 stress */
  double *tau_diss2;            /* Wave to ocean e2 stress */
  double *wave_stke1;           /* Stokes sub-surface velocity, e1 */
  double *wave_stke2;           /* Stokes sub-surface velocity, e2 */
  double *freq;
  int nsfr;

#if defined(HAVE_WAVE_MODULE)
  int do_wave;                  /* Wave flag */
  double wavedt;                /* wave timestep */
#endif

#if defined(HAVE_SEDIMENT_MODULE)
  /* Sediments */
  int do_sed;                   /* Sediment flag */
  int do_bm;                    /* Box model transport flag */
  int do_bmphys;                /* Box model physics flag */

  /* 2D sediment variable pointers */
  char sed_vars[MAXSTRLEN];     /* Sediment size classes */
  char sed_defs[MAXSTRLEN];     /* Sediment attribute defaults */
  double *resuspension;         /* Resuspension at current step */
  double *resuspension_ac;      /* Accumulated resuspension */
  double *deposition_ac;        /* Accumulated deposition */
  double *dzactive;             /* Thickness of sed. active layer [m] */
  double *ustrcw_skin;          /* Skin friction velocity [m/s] */
  double *thickness_sed;        /* Total sediment thickness */
  double *hripples;             /* Ripples height [m] */
  double *lripples;             /* Ripples length [m] */
  double *u1ref, *u2ref;        /* Reference velocity [m/s] */
  /* 3D sediment variable pointers */
  double **por_sed;             /* Sediment bed porosity [nondim] */
  double **cohsedcontent;       /* Mud content in sed bed [Percent] */
  double *svel_floc;            /* Settling velocity of sed. flocs [m/s] */
  double *tss;                  /* Total suspended solids */
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  /* Ecology */
  double ecodt;                 /* Ecology timestep */
  int eco_timestep_ratio;       /* Ecological: physical timestep */
  int do_eco;                   /* Ecology model flag */
  void *e;                      /* Pointer to ecology structure */
#endif

  /* IO */
  int thIO;                     /* Whether to thread I/O */

  /* Number of OMP threads to use in transport mode */
#ifdef HAVE_OMP
  int trans_num_omp;
#endif

  /* Regions */
  int region_nvar;
  int region_nvarS;
  int *region_var;
  int *region_varS;
  double region_dt;
  double region_next;
  int region_mode;
  int region_obcz;

  /* Means for alert tracking */
  int alertf;
  char alertname[MAXSTRLEN];
  double alert_dt;
  int *nalert;
  /* Maximums */
  double meta, seta;
  double mw, sw;
  double mu1, su1;
  double mu1a, su1a;
  double mu2, su2;
  double mu2a, su2a;
  double mdeta, sdeta;
  double mdw, sdw;
  double mdt, sdt;
  double mds, sds;
  double msmin, msmax, ssmin, ssmax;
  double mtmin, mtmax, stmin, stmax;
  double mshear, sshear, memean, semean;
  double mcfl, scfl;
  double ma1, ma2, sa1, sa2;
  double mh1, mh2, sh1, sh2;
  double mv1, mv2, sv1, sv2;
  double mb1, mb2, sb1, sb2;
  double md1, md2, sd1, sd2;
  double mc1, mc2, sc1, sc2;
  /* Thresholds */
  double amax;
  double hmax;
  double vmax;
  double btmax;
  double bcmax;
  double cmax;
  double detamax;
  double dwmax;
  double dtmax;
  double dsmax;
  double smax;
  /* Flags */
  int eta_f;
  int vel2d_f;
  int vel3d_f;
  int wvel_f;
  int tend_f;
  int div2d_f;
  int div3d_f;
  int cfl_f;
  int ts_f;
  int tide_r;
  int shear_f;
  int hdiff_f;
  /* Sparse locations */
  int ceta;
  int cw;
  int ccfl, ncfl;
  int cu1;
  int cu1a;
  int cu2;
  int cu2a;
  int cdeta;
  int cdw;
  int cdt, cds;
  int csmin, csmax;
  int ctmin, ctmax;
  int ca1, ca2;
  int ch1, ch2;
  int cv1, cv2;
  int cb1, cb2;
  int cd1, cd2;
  int cc1, cc2;
  double em;                    /* Excess mass (mean sea level, m) */
  double me;                    /* Mechanical energy (J/m2) */
  double ke;                    /* Kinetic energy (J/m2) */
  double *ef;                   /* OBC energy flux (W/m2) */
  double *vf;                   /* OBC volume flux (m3/s) */

  /* Transport model */
  double *origin;               /* Streamline origin */
  double *pc, *qc, *rc;         /* Streamline origin Courant numbers */
  cstring *trvars;              /* Transport variables */
  int *trvm;                    /* Transport variables map */
  int ntrvars;                  /* Number of transport variables */
  cstring *trvarsS;             /* 2D Transport variables */
  int *trvmS;                   /* 2D Transport variables map */
  int ntrvarsS;                 /* Number of 2D transport variables */
  double *vol_cons;             /* Volume conservation diagnostic */
  int togn;                     /* Origin for tri-linear interpolation */
  
  /* Miscillaneous */
  xytoij_tree_t *xyij_tree;     /* A xytoij_tree_t for XYtoIJ routines */
  void *custdata;               /* Point to a custom data handle */
  int nstep;                    /* Number of 3D model time-steps */
  int nstep2d;                  /* Number of 2D model time-steps */
  int *kth_e1;                  /* Thin layer vector for e1 velocity */
  int nkth_e1;                  /* Number of cells in kthin, e1 */
  int *kth_e2;                  /* Thin layer vector for e2 velocity */
  int nkth_e2;                  /* Number of cells in kthin, e2 */
  int ndf;                      /* Number of dump files */
  dump_file_t *dumplist;        /* Table of dump files */
  int is_filled;                /* Is the master filled? */
  char orthoweights[MAXSTRLEN]; /* csr tide model orthoweights path */
  char nodal_dir[MAXSTRLEN];    /* Tidal nodal corrections path */
  char tide_con_file[MAXSTRLEN];/* Tidal constituent file */
  char bdrypath[MAXSTRLEN];     /* Path for boundary files */
  int *trinc;                   /* 3D Tracer increment flag */
  int *trincS;                  /* 2D Tracer increment flag */
  int df_diagn_set;             /* Flag if the `diagn' flag has been set */
                                /* for one of the output files. */
  int ***i1, ***i2;
  double *d2, *d3;

  /* Function to calculate Vz and Kz.  */
  void (*calc_closure) (geometry_t *, window_t *, win_priv_t *);
  void (*s_func) (double aN, double aM, double *cmu, double *cmu_d);

  /*UR-ADDED late clean up pointer, remains NULL until run is finished */
  free_stack_t* free_stack;

  /*UR-ADDED a list of functions which selfsufficiently initialise and exchange data per timestep*/ 
 custom_stack_t* custom_fnstack;

  /* The MPI Id of this process, 0 for non-mpi */
  int mpi_rank;

  /* function pointers for transfers */
  void (*win_data_fill_3d)  (master_t *master, geometry_t *window, window_t *windat, int nwindows);
  void (*win_data_fill_2d)  (master_t *master, geometry_t *window, window_t *windat, int nwindows);
  void (*win_data_refill_3d) 
         (master_t *master, geometry_t *window, window_t *windat, int nwindow, int mode);
  void (*win_data_refill_2d) 
         (master_t *master, geometry_t *window, window_t *windat, int nwindow, int mode);
  void (*win_data_empty_3d) (master_t *master, geometry_t *window, window_t *windat, int mode);
  void (*win_data_empty_2d) (master_t *master, geometry_t *window, window_t *windat, int mode);
  void (*update_master) (master_t *master, window_t **windat, int mode);
  void (*master_fill) (master_t *master, geometry_t **window, 
		       window_t **windat, win_priv_t **wincon);
  void (*master_fill_ts) (master_t *master, geometry_t **window, 
			  window_t **windat, win_priv_t **wincon);
  void (*master_fill_glider) (master_t *master, geometry_t **window, 
			      window_t **windat, win_priv_t **wincon, ts_point_t *ts, double t);
};


/*------------------------------------------------------------------*/
/* Window data structure                                            */
/*------------------------------------------------------------------*/
struct window {

  int ntr;                      /* Number of 3D tracers */
  int ntrS;                     /* Number of 2D tracers */
  int nsed;                     /* Number of sediment variables */
  int ntre;                     /* Number of 3D edge tracers */
  int ntreS;                    /* Number of 2D edge tracers */
  int nstep;                    /* Iteration step */
  double **tr_wc;               /* 3D Tracers */
  double **tr_wcS;              /* 2D Tracers */
  double ***tr_sed;             /* Sediment Tracers */
  double **tr_wc_as;            /* Multi-dt cells tracer time t values */
  double **tr_wc_ae;            /* Multi-dt cells tracer time t+1 values */
  double *temp;                 /* Temperature */
  double *sal;                  /* Salinity */
  int tno;                      /* Tracer number for temperature */
  int sno;                      /* Tracer number for salinity */
  double *sdc;                  /* Smagorinsky horizontal diffusion */
  double t;                     /* Simulation time (s) */
  double dt;                    /* Time step for the 3D mode (s) */
  double dtu1;                  /* Time step for the 3D u1 velocity (s) */
  double dtu2;                  /* Time step for the 3D u2 velocity (s) */
  double dt2d;                  /* Time step for the 2D mode (s) */
  double dts;                   /* Window 3D sub-time step (s) */
  double dt2s;                  /* Window 2D sub-time step (s) */
  double dtf;                   /* Forward part of leapfrog dt */
  double dtb;                   /* Backward part of leapfrog dt */
  double dtf2;                  /* Forward part of leapfrog dt2d */
  double dtb2;                  /* Backward part of leapfrog dt2d */
  double dttr;                  /* Time step for the tracers (s) */
  double days;                  /* Time in days (t/86400) */
  double wclk;                  /* Mean CPU time for this window */
  double rampval;               /* Ramp value for forcing */
  short iratio;                 /* Ratio of 3D /2D time-steps */
  double mcfl2d;                /* Mean minimum 2D cfl time-step */
  double mcfl3d;                /* Mean minimum 3D cfl time-step */
  int cflc;                     /* Location of min 3D cfl for ACTIVE3D */
  double *eta_as;               /* Multi-dt cells eta time t values */
  double *eta_ae;               /* Multi-dt cells eta time t+1 values */
  double *u1av_as;              /* Multi-dt cells u1av time t values */
  double *u1av_ae;              /* Multi-dt cells u1av time t+1 values */
  double *u2av_as;              /* Multi-dt cells u2av time t values */
  double *u2av_ae;              /* Multi-dt cells u2av time t+1 values */
  double *wflux;                /* Water flux for sourcesinks */
  double **tflux;               /* Tracer flux for sourcesinks */
  double *u1;                   /* Normal velocity to edge (ms-1) */
  double *u;                    /* Eastward velocity (ms-1) */
  double *u2;                   /* Tangential velocity to edge (ms-1) */
  double *v;                    /* Northward velocity (ms-1) */
  double *nu1;                  /* Updated 3D e1 velocity (ms-1) */
  double *nu2;                  /* Updated 3D e2 velocity (ms-1) */
  double *w;                    /* Vertical velocity (ms-1) */
  double *wtop;                 /* Surface vertical velocity (ms-1) */
  double *wbot;                 /* Bottom vertical velocity (ms-1) */
  double *u1av;                 /* e1 2D velocity (ms-1) */
  double *u2av;                 /* e2 2D velocity (ms-1) */
  double *nu1av;                /* Updated 2D e1 velocity (ms-1) */
  double *nu2av;                /* Updated 2D e2 velocity (ms-1) */
  double *uav;                  /* Eastward 2D velocity (ms-1) */
  double *vav;                  /* Northward 2D velocity (ms-1) */
  double *eta;                  /* Height of the free surface at t+1 (m) */
  double *dzu1;                 /* Layer thickness at e1 face (m) */
  double *dzu2;                 /* Layer thickness at e2 face (m) */
  double *depth_e1;             /* Total depth at e1 face (m) */
  double *depth_e2;             /* Total depth at e2 face (m) */
  double *Kz;                   /* Vertical eddy diffusivity (m2s-1) */
  double *Vz;                   /* Vertical eddy viscosity (m2s-1) */
  double *Kq;                   /* Turbulence intensity diffusivity */
  double *tke;                  /* Turbulent kinetic energy */
  double *diss;                 /* Dissipation */
  double *omega;                /* Turbulent frequency */
  double *L;                    /* Turbulence length scale */
  double *Q2;                   /* Turbulence intensity */
  double *Q2L;                  /* Turbulence intensity length scale */
  double *mom_bal;              /* Momentum balance maximums */
  double *u1_adv;               /* u1 velocity advective tendency */
  double *u1_hdif;              /* u1 velocity horizontal diff. tend.  */
  double *u1_vdif;              /* u1 velocity vertical diff. tend.  */
  double *u1_cor;               /* u1 velocity Coriolis tendency */
  double *u1_btp;               /* u1 velocity barotropic pressure tend. */
  double *u1_bcp;               /* u1 velocity baroclinic pressure tend. */
  double *u1_sto;               /* u1 velocity Stokes drift tendency */
  double *u2_adv;               /* u1 velocity advective tendency */
  double *u2_hdif;              /* u1 velocity horizontal diff. tend.  */
  double *u2_vdif;              /* u1 velocity vertical diff. tend.  */
  double *u2_cor;               /* u1 velocity Coriolis tendency */
  double *u2_btp;               /* u1 velocity barotropic pressure tend. */
  double *u2_bcp;               /* u1 velocity baroclinic pressure tend. */
  double *u2_sto;               /* u2 velocity Stokes drift tendency */
  double *tr_adv;               /* Tracers advective tendency */
  double *tr_hdif;              /* Tracers horizontal diffusive tendency */
  double *tr_vdif;              /* Tracers vertical diffusive tendency */
  double *tr_ncon;              /* Tracers non-conservative tendency */
  double *u1flux3d;             /* Flux in the e1 direction (m3s-1) */
  double *u2flux3d;             /* Flux in the e2 direction (m3s-1) */
  double *u1flux;               /* Volume flux through e1 2D faces (m3) */
  double *u2flux;               /* Volume flux through e2 2D faces (m3) */
  double *topz;                 /* Height of the surface (m) */
  double *u1bot;                /* Bottom u1 velocity (ms-1) */
  double *u2bot;                /* Bottom u2 velocity (ms-1) */
  double *wdiff2d;              /* Surface slope term */
  double *detadt;               /* Rate of change of elevation (ms-1) */
  double *dens;                 /* Density (kgm-3) */
  double *dens_0;               /* Density at zero pressure (kgm-3) */
  double *wind1;                /* Wind spped in the x direction (ms-1) */
  double *wind2;                /* Wind spped in the y direction (ms-1) */
  double *windspeed;            /* Wind speed (ms-1) */
  double *winddir;              /* Wind direction (deg T) */
  double *cloud;                /* Cloud amount (oktas) */
  double *airtemp;              /* Air temperature (deg C) */
  double *patm;                 /* Atmospheric pressure (HPa) */
  double *precip;               /* Precipitation (mm day-1) */
  double *evap;                 /* Evaporation (mm day-1) */
  double *heatf;                /* Surface heat flux (Wm-2) */
  double *swr;                  /* Short wave radiation (Wm-2) */
  double *light;                /* Daily mean short wave radiation */
  double *hftemp;               /* Sea surface temperature (deg C) */
  double *wetb;                 /* Wet bulb temperature (deg C) */
  double *rh;                   /* Relative humidity (%) */
  double *waterss;              /* Volume sourse / sink term */
  double *waterss2d;            /* Volume sourse / sink term */
  double *vd;                   /* Density profile work array */
  /*double *eta_rlx;*/              /* Relaxation  source/sink */
  /*double etarlxtc;*/              /* Elevation relaxation time constant */
  int *trinc;                   /* 3D Tracer increment flag */
  int *trincS;                  /* 2D Tracer increment flag */
  int *sur_e1;                  /* e1 surface locations - windows */
  int *sur_e2;                  /* e2 surface locations - windows */

  /* 2D diagnostic variables */
  double *cfl2d;                /* 2D cfl time-step */
  double *cfl3d;                /* 3D cfl time-step */
  double *cour;                 /* Courant stability */
  double *lips;                 /* Lipshitz stability */
  double *ahsb;                 /* Horizontal diffusion stability */
  double *mixl;                 /* Pointer to mixed layer depth */
  double *nsfd;                 /* Net salt flux diagnostic */
  double *nhfd;                 /* Net heat flux diagnostic */
  double *swrd;                 /* Short wave radiation diagnostic */
  double *lwrd;                 /* Long wave radiation diagnostic */
  double *shfd;                 /* Sensible heat flux diagnostic */
  double *lhfd;                 /* Latent heat flux diagnostic */
  double *lwro;                 /* Long wave output radiation */
  double *avhrr;                /* AVHRR SST */
  double *ghrsst;               /* GHRSST SST */
  double *shwin;                /* Window partitioning */
  double *alert_a;              /* Actual alert diagnostic */
  double *alert_c;              /* Cumulative alert diagnostic */
  double *u1vhin;               /* Initial e1 horizontal viscosity */
  double *u2vhin;               /* Initial e2 horizontal viscosity */
  double *swr_attn;             /* Short wave attenuation (for red) */
  double *swr_attn1;            /* Short wave attenuation for blue-green */
  double *swr_tran;             /* Short wave surface transmission */
  double *swr_babs;             /* Short wave bottom absorption */

  /* To keep track of heatflux diagnostics */
  int lwrn;                     /* Tracer number for lwr */
  int lhfn;                     /* Tracer number for lhf */
  int shfn;                     /* Tracer number for shf */

  /* Salt flux diagnostics */
  int evapn;                    /* Tracer number for evaporation */
  int precipn;                  /* Tracer number for precipitation */

  /* Unstructured */
  double *circ;                 /* Circulation */
  double *kec;                  /* Kinetic energy at cell centres */
  double *div;                  /* Divergence at cell centres */
  double *rvor;                 /* Relative vorticity */
  double *fv;                   /* Coriolis parameter at vertices */
  double *nrvor;                /* Normalized relative vorticity */
  double *npvor;                /* Normalized potential vorticity */
  double *nrvore;               /* Normalized relative vorticity on edges */
  double *npvore;               /* Normalized potential vorticity on edges */
  double *nrvorc;               /* Normalized relative vorticity on centres */

  /* 3D diagnostic variables */
  double *fluxe1;               /* Advective flux in the e1 direction */
  double *fluxe2;               /* Advective flux in the e1 direction */
  double *fluxw;                /* Vertical advective flux */
  double *fluxkz;               /* Vertical diffusive flux */
  double *u1m;                  /* Pointer to u mean velocity */
  double *u2m;                  /* Pointer to v mean velocity */
  double *ume;                  /* Mean 3D velocity */
  double *wm;                   /* Pointer to w mean velocity */
  double *u1vm;                 /* Pointer to u1 mean volume flux */
  double *u2vm;                 /* Pointer to u2 mean volume flux */
  double *Kzm;                  /* Pointer to mean Kz */
  double *etam;                 /* Pointer to mean elevation */
  double *u1am;                 /* Pointer to uav mean velocity */
  double *u2am;                 /* Pointer to vav mean velocity */
  double *uame;                 /* Mean 2D velocity */
  double *w1m;                  /* Pointer to eastward mean wind */
  double *w2m;                  /* Pointer to northward mean wind */
  double *tempm;                /* Pointer to mean temperature */
  double *saltm;                /* Pointer to mean salinity */
  double *tram;                 /* Pointer to mean tracer */
  double *meanc;                /* Iteration counter for means */
  double *perc;                 /* Pointer to percentile diagnostic */
  double *odeta;                /* detadt at previous timestep */
  double *rv;                   /* Relative vorticity diagnostic */
  double *av;                   /* Absolute vorticity diagnostic */
  double *pv;                   /* Potential vorticity diagnostic */
  double *rv_drvdt;             /* Rate of change rv tendency */
  double *rv_nonlin;            /* Nonlinear rv tendency (adv + diff) */
  double *rv_beta;              /* Planetary vorticity rv tendency */
  double *rv_strch;             /* Topographic stretching rv tendency */
  double *rv_jebar;             /* JEBAR rv tendency */
  double *rv_wsc;               /* Wind stress curl rv tendency */
  double *rv_bsc;               /* Bottom stress curl rv tendency */
  double *riverflow;            /* River flow diagnostic */
  double *riverdepth;           /* River depth diagnostic */
  double *riversalt;            /* River ghost salinity diagnostic */
  double *steric;               /* Steric height diagnostic */
  double *fltr;                 /* Pointer to flushing tracer */
  double *agetr;                /* Pointer to age tracer */
  double *stream;               /* 2D streamfunction */
  double *brunt;                /* Brunt Vaisala (buoyancy) frequency (s-1) */
  double *int_wave;             /* Internal wave speed (ms-1) */
  double *rich_gr;              /* Gradient Richardson number */
  double *rich_fl;              /* Flux Richardson number */
  double *reynolds;             /* Reynolds number */
  double *froude;               /* Froude number */
  double *sep;                  /* Surface Ekman pumping */
  double *bep;                  /* Bottom Ekman pumping */
  double *tfront;               /* Simpson-Hunter tidal front */
  double *sigma_t;              /* Sigma_t */
  double *rossby_in;            /* Internal Rossby radius (m) */
  double *rossby_ex;            /* External Rossby radius (m) */
  double *speed_2d;             /* 2D current speed */
  double *speed_3d;             /* 3D current speed */
  double *speed_sq;             /* Bottom speed squared */
  double *energy;               /* Total energy (Jm-2) */
  double *kenergy;              /* Kinetic energy (Jm-2) */
  double *obc_phase;            /* OBC phase speed (m/s) */
  double *nprof;                /* Normalized profile */
  double *sound;                /* Speed of sound */
  double *schan;                /* Sound channel depth */
  double *sonic;                /* Sonic layer depth */
  double *wetcell;              /* Wet cell diagnostic */
  double *slope_x;              /* Surface slope (e1 direction) diagnostic */
  double *slope_y;              /* Surface slope (e2 direction) diagnostic */
  double *surfz;                /* Surface layer diagnostic */
  double *wind_Cd;              /* Wind drag */
  double *shear_v;              /* Vertical shear (s-1) */
  double *b_prod;               /* Buoyancy production (m2s-2) */
  double *s_prod;               /* Shear production (m2s-2) */
  double *otemp;                /* OFAM temperature */
  double *osalt;                /* OFAM salinity */
  double *rtemp;                /* Relaxation temperature */
  double *rsalt;                /* Relaxation salinity */
  double *temp_tc;              /* Relaxation temperature time constant */
  double *salt_tc;              /* Relaxation salinity time constant */
  double *eta_tc;               /* Relaxation eta time constant */
  double *eta_inc;              /* Relaxation eta increment */
  double *layth;                /* Layer thickness for sigma */
  double *u1_rad;               /* u1 velocity radiation stress tend. */
  double *u2_rad;               /* u2 velocity radiation stress tend. */
  double *tau_be1;              /* Bottom stress in e1 direction */
  double *tau_be2;              /* Bottom stress in e2 direction */
  double *tau_bm;               /* Bottom stress magnitude */
  double *regionid;             /* Region ids */
  double *regres;               /* Region residence time */
  double tvol;                  /* Total volume */
  double tarea;                 /* Total area */
  double tmass;                 /* Total mass */
  double *trtot;                /* Total tracers */
  double *dum1, *dum2, *dum3;   /* Dummy arrays */
  double *vcorr, *acorr;        /* Local transport fill corrections */
  double *Vi;                   /* Volume error */
  double *unit;                 /* Unit tracer */
  double *glider;               /* Glider density */
  double *u1vhc;                /* Cell centered horizontal viscosity */
  double *reefe1;               /* Pointer to e1 reef fraction tracer */
  double *reefe2;               /* Pointer to e2 reef fraction tracer */
  int *totid;                   /* Tracer numbers for 3D totals */
  int *totidS;                  /* Tracer numbers for 2D totals */
  int *totid_sed;               /* Sediment tracer numbers for totals */
  double *sederr;               /* Error percentage map for sediments */
  double *ecoerr;               /* Error percentage map for ecology */
  double *decv1;                /* Decorrelation length scale */
  double *dhd;                  /* Degree heating day */
  double *dhwc;                 /* Offset degree heating day */
  double *dhw;                  /* Degree heating day */
  double sederrstep;            /* Sediment error step */
  double ecoerrstep;            /* Ecology error step */
  int ntot;                     /* Number of additional total tracers */

  /* Waves */
  double *ustrcw;               /* Wave-current friction velocity at
                                   bottom */
  double *wave_ub;              /* Near bottom wave orbital velocity
                                   [ms-1] */
  double *wave_period;          /* Surface wave period [s] */
  double *wave_dir;             /* Surface wave direction [radian] */
  double *wave_amp;             /* Surface wave amplitude [m] */
  double *wave_Sxy;             /* e1 tangential radiation stress */
  double *wave_Syx;             /* e2 tangential radiation stress */
  double *wave_Fx;              /* e1 wave-induced force */
  double *wave_Fy;              /* e2 wave-induced force */
  double *wave_Cd;              /* Wave enhanced bottom drag */
  double *wave_ste1;            /* Stokes surface velocity, e1 */
  double *wave_ste2;            /* Stokes surface velocity, e2 */
  double *tau_w1;               /* Wave supported wind e1 stress */
  double *tau_w2;               /* Wave supported wind e2 stress */
  double *tau_diss1;            /* Wave to ocean e1 stress */
  double *tau_diss2;            /* Wave to ocean e2 stress */
  double *wave_stke1;           /* Stokes sub-surface velocity, e1 */
  double *wave_stke2;           /* Stokes sub-surface velocity, e2 */
  double *freq;
  int nsfr;

#if defined(HAVE_SEDIMENT_MODULE)
  /* Sediments */
  /* 2D sediment variable pointers */
  double *resuspension;         /* Resuspension at current step */
  double *resuspension_ac;      /* Accumulated resuspension */
  double *deposition_ac;        /* Accumulated deposition */
  double *dzactive;             /* Thickness of sed. active layer [m] */
  double *ustrcw_skin;          /* Skin friction velocity [m/s] */
  double *thickness_sed;        /* Total sediment thickness */
  double *hripples;             /* Ripples height [m] */
  double *lripples;             /* Ripples length [m] */
  double *u1ref, *u2ref;        /* Reference velocity [m/s] */
  /* 3D sediment variable pointers */
  double **por_sed;             /* Sediment bed porosity [nondim] */
  double **cohsedcontent;       /* Mud content in sed bed [Percent] */
  double *svel_floc;            /* Settling velocity of sed. flocs [m/s] */
  double *tss;                  /* Total suspended solids */
#endif

  /* Leapfrog variables */
  double *etab;                 /* Height of the free surface at t-1 (m) */
  double *u1avb;                /* 2D u1 velocity at previous timestep */
  double *u2avb;                /* 2D u2 velocity at previous timestep */
  double *u1b;                  /* 3D u1 velocity at previous timestep */
  double *u2b;                  /* 3D u2 velocity at previous timestep */

  /* Relaxation */
  relax_info_t *vel_rlx;       /* Velocity relaxation */
  relax_info_t *eta_rlx;       /* Elevation relaxation */

  /* Alert tracking */
  int *nalert;
  /* Maximums */
  double meta;
  double mw;
  double mu1;
  double mu1a;
  double mu2;
  double mu2a;
  double mdeta;
  double mdw;
  double mdt;
  double mds;
  double msmin, msmax, ssmin, ssmax;
  double mtmin, mtmax, stmin, stmax;
  double mshear, sshear, memean, semean;
  double ma1, ma2;
  double mh1, mh2;
  double mv1, mv2;
  double mb1, mb2;
  double md1, md2;
  double mc1, mc2;
  double mcfl;
  /* Sparse locations */
  int *ceta;
  int cw;
  int ccfl;
  int *cu1;
  int *cu1a;
  int *cu2;
  int *cu2a;
  int cdeta;
  int cdw;
  int cdt, cds;
  int csmin, csmax;
  int ctmin, ctmax;
  int ca1, ca2;
  int ch1, ch2;
  int cv1, cv2;
  int cb1, cb2;
  int cd1, cd2;
  int cc1, cc2;
  double *tempb;
  double *salb;
  double em;                    /* Excess mass (mean sea level, m) */
  double me;                    /* Mechanical energy (J/m2) */
  double ke;                    /* Kinetic energy (J/m2) */
  double *ef;                   /* OBC energy flux (W/m2) */
  double *vf;                   /* OBC volume flux (m3/s) */

  /* Transport model */
  double *origin;               /* Streamline origin */
  double *pc, *qc, *rc;         /* Streamline origin Courant numbers */
  double *vol_cons;             /* Volume conservation diagnostic */
  delaunay **d;                 /* Delaunay data structure for transport */
  GRID_SPECS **gsx, **gsy, **gsz;  /* u, v interpolation structures */
  GRID_SPECS **gst;              /* Tracer interpolation structures */
  npt_t **nptt;                  /* Tracer point triangle structure */
  npt_t **nptm;                  /* Momentum point triangle structure */

  /* Miscillaneous */
  int df_diagn_set;             /* Flag if the `diagn' flag has been set */
  int regf;                     /* Run regulation flag */
  /* for one of the output files */
};


/* timeseries */
struct ts_point{
  char pname[MAXSTRLEN];
  double x;                     /* x coord of point */
  double y;                     /* y coord of point */
  double z;                     /* z coord of point */
  double tsdt;                  /* TS output time interval */
  double tsout;                 /* Time of next output */
  char tsunits[MAXSTRLEN];
  FILE *fp;                     /* Pointer to the timeseries file */
  master_t *master;             /* Grid in which the point is located */
  int i;                        /* i value of cell containing point */
  int j;                        /* j value of cell containing point */
  int k;                        /* k value of cell containing point */
  int c;                        /* Unstructured value of cell */
  long stepnum;                 /* step number */
  double sinth;                 /* sin value for velocity calculations */
  double costh;                 /* cos value for velocity calculations */
  int v_offset;                 /* Flag to set the vertical offset */
  int type;                     /* Type of output */
  int nvars;                    /* Number of tracer output variables */
  int vars[MAXNUMVARS];         /* Tracer output variable indices */
  int var_type[MAXNUMVARS];     /* Tracer 2D or 3D flag */
  int da_cycle;                 /* DA */
  /*
   * These are for when the output is asked for at varying times
   * and locations, such as glider data
   */
  timeseries_t ts; /* Input ts file */
  int varids[3];   /* lon,lat,depth */
  /* Inline data comparison */
  int ndata;                  /* Number of forcing files              */
  timeseries_t **tsdata;      /* Forcing files                        */
  int **dvarids;              /* Forcing files variable id's          */
  int dnvars;                 /* Number of forcing variables          */
  char **dvars;               /* Name of forcing variables            */
  char **dunits;              /* Units of forcing variables           */
  double **data;              /* Pointer to master data               */
  int *ddim;                  /* 0=2D, 1=3D                           */
  int metric;                 /* Type of metric                       */
  int kernal;                 /* Kernal size                          */
  int last;                   /* Last location in grid                */
  double *val;                /* Data-model comparison value          */
  double *obs;                /* Observation value                    */
  double *minv;               /* Minimum value                        */
  double *maxv;               /* Maximum value                        */
  double *nvals;              /* Number of values read                */
  double **thresh;            /* Threshold value                      */
};

/*------------------------------------------------------------------*/
/* Transport mode structure.
/*------------------------------------------------------------------*/
struct transport {
  geometry_t *tpg;
  win_priv_t *tpc;
  window_t *tpd;
  double *w1;               /* Dummy array for input reads              */
  int ns2;                  /* 2D sparse size                           */
  int ns3;                  /* 3D sparse size                           */
  xytoij_tree_t *xyij_tree; /* A xytoij_tree_t for XYtoIJ routines      */
};

struct hd_data{
  geometry_t *geom;
  master_t *master;
  dump_data_t *dumpdata;
  parameters_t *params;
  geometry_t **window;
  window_t **windat;
  win_priv_t **wincon;
};


/*UR-ADDED  data structure for late destroy of memory
 *  this meant to be executed after everything else but just before the
 * master. If the free function is not NULL, free is called.
 */


struct free_stack {
  free_stack_t* next;
  free_stack_t* last;
  void* object;
  void (*free) (void *object);
};


struct custom_function {
  custom_data_t* data;
  char* name;
  custom_init_funtion init;
  custom_scatter_funtion scatter;
  custom_gather_function gather;
  custom_destroy_function destroy;
  custom_function_do_scatter do_scatter;
  custom_function_do_gather do_gather;
  custom_create_data  create;
  custom_function_t* next;
};


struct custom_stack {
    custom_function_t* functions;
    int nfunctions;
};

struct custom_data {
  void* data;
  int index;
  custom_data_t* next;
};

#endif                          /* _SPARSE_H */
