/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/sediments/include/sediments.h
 *  
 *  Description:
 *  Header file for sediment model
 *  
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sediments.h 7484 2024-02-13 23:53:02Z mar644 $
 *
 */

#ifndef _SEDIMENTS_H
#define _SEDIMENTS_H

/* include stdio for all files */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include "ems.h"

/*UR adjust log redirects */

typedef int (*sedlogfn) (int level, char *format, ...);
typedef int (*sedlogtag) (int level, char *tag, char *format, ...);

extern sedlogfn  sedlog;
extern sedlogtag sedtag;

/* Release verions and getters */
#define SEDIMENTS_MAJOR_VERSION 1
#define SEDIMENTS_MINOR_VERSION 1
#define SEDIMENTS_PATCH_VERSION 3

int get_sediments_major_vers(void);
int get_sediments_minor_vers(void);
int get_sediments_patch_vers(void);

/* set to 0 to run standalone, 1 otherwise */
#define AS_EMSLIB_MODULE 1

#if AS_EMSLIB_MODULE
/* emslib related stuff
 *
 * include loglevels etc.
 */
#include <emslogger.h>

/* this include will also take care of HAVE_SEDIMENT_MODULE */
#include <ems_conf.h>


#else

#define HAVE_SEDIMENT_MODULE 1

/* Log levels  */
#define LTRACE 40
#define LMAIN 0
#define LDEBUG 30
#define LWARN 20
#define LINFO 10
#define LMETRIC 50
#define LERROR -10
#define LFATAL -20
#define LPANIC -30

#endif

/*UR
 *
 * NM - for consistency this should move in above redirects and be fixed up
 *
 */
#include "emsalloc.h"

/*UR - end log redirects
 */

#if !defined(PI)
#define PI (3.14159265358979323846)
#endif

#if !defined(g_grav)
#define g_grav (9.81)
#endif

#if !defined(max)
#define max(x,y) ((x)>(y) ? (x) : (y) )
#endif

#if !defined(min)
#define min(x,y) ((x)<(y) ? (x) : (y) )
#endif


#if !defined(MINVAL)
#define MINVAL (1.e-11)
#endif

/* pre-defined structures and typedefs */
struct sediment;
struct sed_params;
struct sed_column;
struct sed_tracer;
struct sed_spatial;

typedef struct sediment sediment_t;
typedef struct sed_params sed_params_t;
typedef struct sed_column sed_column_t;
typedef struct sed_tracer sed_tracer_t;
typedef struct sed_spatial sed_spatial_t;

/* prototypes */
void bbl(sediment_t *sediment, sed_column_t *sm);
void vert_transport(sediment_t *sediment, sed_column_t *sm);
void hd2sed(sediment_t *sediment, sed_column_t *sm, int c);
void sed2hd(sediment_t *sediment, sed_column_t *sm, int c);
void sed_limits(sediment_t *sediment, sed_column_t *sm, int c);
sed_column_t *alloc_sed_column(sediment_t *sediment);

/* Interface functions */
extern void i_set_error(void* hmodel, int col, int errorf, char *text);
extern double i_get_psize(void* model, char *name);

/* define structures */
/*NMYfix*/
struct sediment {
  void *hmodel;
  sed_params_t *msparam;        /* grid independent parameters */
  sed_tracer_t *mstracers;      /* tracer array */
  sed_column_t *sm;             /* grid dependent parameters */
  sed_column_t **smn;           /* parallel version */
  sed_spatial_t *spatial;
};

/* Sediment transport tracer attribute private data structure */
typedef struct {
  int type;                     /* Private data dype */
  char name[MAXSTRLEN];         /* Name of the default list */
  char trname[MAXSTRLEN];       /* Name of the tracer */
  int cohesive;                 /* Cohesive flag */
  int calcvol;                  
  int floc;                     /* Flocculation flag */
  int resuspend;                /* Resuspension flag */
  int deposit;                  /* Deposition flag */
  int adsorb;                   /* Adsorbtion flag */
  int obc;                      /* Open boundary condition */
  double psize;                 /* Particle size */
  double b_dens;                /* Particle bulk density */
  double i_conc;                /* Initial deposit concentration */
  double f_conc;                /* Compacted deposit concentration; 
				   default = i_conc */
  //2019                                   default = i_conc */
  double css_erosion;
  double css_deposition;

  double svel;                  /* Constant settling velocity */
  char svel_name[MAXSTRLEN];    /* Settling velocity name */
  double adsorbkd;              /* Adsorbtion Kd */
  double adsorbrate;            /* Adsorbtion rate */
  char carriername[MAXSTRLEN];  /* particulate tracer carrying sediment reactive tracer */
  char dissolvedname[MAXSTRLEN];/* dissolved coutrepart of the sediment reactive tracer */
} trinfo_priv_sed_t;

/* (i,j) independent data */
struct sed_params {
  double t;                     /* model time in seconds */
  double dt;                    /* current 3d step - time step */
  int nstep;                    /* currently not used */
  int nz, sednz;                /* number of cells in z direction in wc
                                   and sedibed */
  int ncol;                     /* number of columns (ncol>1 for 3-d
                                   applications) */
  double max_thick_sed;    //2010         
  double *css_er_val, *css_er_depth;


  int ntr;                      /* number of tracers */

  int n_hripple, n_lripple, n_ustrcw_skin, n_depth_sedi, n_dzactive;
  int n_erdepflux_total, n_erdepflux_total_ac, n_tss, n_svel_floc;
  int n_erdepflux_oxygen, n_erdepflux_oxygen_ac;
  int n_por_sed, n_coh_sed;

  /* sediments */
  int geomorph;                 /* calculation mode conserving (1) or not
                                   coserving (0) liquid volume */
  int calcvol_wc;               /* number of volumetric tracers in water
                                   col */
  int calcvol_sed;              /* number of volumetric tracers in sedi
                                   bed */
  double maxseddz;              /* max thickness of sed layer[m]
                                   (currently not used) */
  double minseddz;              /* min thickness of sed layer [m] */
  double css;                   /* Common (default) crititical stress
                                   value N/m*m */
  int cssmode;                  /* Erosion calculation mode (0-3) */
  int flocmode;                 /* Flocculation calculation mode (0, 1) */
  double flocprm1, flocprm2;    /* Flocculation parameters */

  // NMY 2024 floc Amanda
  //char flocfile[MAXSTRLEN];
  int    floc_Nlines;
  int    floc_imax,   floc_jmax,   floc_kmax;
  double floc_salmax, floc_turmax, floc_gggmax;
  double floc_salmin, floc_turmin, floc_gggmin;
  double *floc_d50_1d, ***floc_d50_3d;


  double css_dep, css_scale;    /* Deposition critical stress N/m*m */
  double consolrate;            /* Consolidation time [day] */
  int consolidate;              /* Activators */
  double physriph, physripl;    /* Initial physical ripples height and
                                   length [m] */
  double bioriph, bioripl;      /* Initial bioripples height and length
                                   [m] */
  /* Bio-mixing */
  double biodens;               /* Biological activity (standard animals
                                   per square metre) */
  double maxbiodepth;           /* Max depth for biological activity [m] */
  /* Bio-irrigation rates */
  double bi_dissol_kz;          /* Dissolved diffusion cfft [m*m/s] */
  /* Bio-turbation rates */
  double bt_partic_kz;          /* Particulate diffusion cfft [m*m/s] */

  double bi_dissol_kz_i;        /* surface-water bioirrig exchange cfft [m*m/s] */
  double bt_partic_kz_i;

  /* Spatially varying versions of some of the above */
  double *css_scale_spv;
  double *physriph_spv;
  double *physripl_spv;
  double *maxbiodepth_spv;
  double *biodens_spv;
  double *bi_dissol_kz_spv;
  double *bt_partic_kz_spv;
  double *bi_dissol_kz_i_spv;
  double *bt_partic_kz_i_spv;

  char biosedprofile;           /* Profile of activity in sediments */
  /* Variables used by waves routines */
  int wave_enhanced_bf;         /* not used */
  int bbl_nonlinear;            /* nonlinear (1), linear (0) */
  int hindered_svel_patch;      /* patch to hindered settling */
  int hindered_svel;            /* key to invoke a hindered settling */
  double reef_scale_depth;
  int calc_ripples;             /* Activator */
  /* friction */
  double quad_bfc;              /* minimum quadratic friction coeff */
  double uf;                    /* background friction velocity [m/s],
                                   currently not used */
  /* miscellaneous values */
  double minpor_sed, minpor_wc, finpor_sed; /* Minimal and final
                                               porosities */
  double mindepth_sedi;          /* Minimal sediment thickness [m] */
  char gridkey_buffer[MAXSTRLEN]; /* Buffer for gridkey func. currently
                                     not used */
    int verbose_sed;

  /*UR runtime optimizations avoiding conditional statements */
  /* indices for tracers which are partic not adsorb and floculate */
  int* vtransp_floc;
  int n_vtransp_floc;

  /* indices for tracers which are not diagnostic and decay */
  int* vtransp_decay;
  int n_vtransp_decay;

  /* indices for tracers which are not diagnostic and are in sediment and adsorb */
  int* vtransp_adsorb;
  int n_vtransp_adsorb;

  /* indices for tracers which are partic, not adsorb, calcvol_wc and are in sediment */
  int* vtransp_por_wc;
  int n_vtransp_por_wc;

  /* indices for tracers which are partic, not adsorb, calcvol_sed and are in sediment */
  int* vtransp_por_sed;
  int n_vtransp_por_sed;

  /* indices for tracers which are not diagnostic and are in sediment and are dissolved */
  int* vtransp_dissolv;
  int n_vtransp_dissolv;

  /* indices for tracers which are partic, not adsorb and are in sediment */
  int* vtransp_insed_partic;
  int n_vtransp_insed_partic;

  /* indices for tracers which are partic, not adsorb, not diagnostic and are in sediment */
  int* vtransp_nd_insed_partic;
  int n_vtransp_nd_insed_partic;


  /* indices for tracers which are partic, not adsorb */
  int* vtransp_partic;
  int n_vtransp_partic;

 //2012
  char **prmnameS; // spatially varying prm names
  double ***prmpointS; // pointers to spatial prm
  int nprmS, *prmindexS, *trindexS; //number of spatial prm and prm,tracer indexes
  int ntrB; //number of benthic (ie 2d) tracers
  char  **trnameB; //names of benthis tracers

  int hardsub_numb; //hard substrate tracer number

  int ship_N; //average number of ships per day
  int ship_C; //number of grid-cells on a ship track
  double ship_T; //average hours a ship spends on track
  int *ship_i, *ship_j; // arrays of ship-track i,j cells
  double ship_Kz; // ship-induced vertical diffusion coefficient

//2016
 double erflux_scale; // erision rate scaling parameter, default 1;

 int *fluxsedimap_inst, *fluxsedimap_ac;                /* 2018 map from 2D to 3D variables */

};

/* (i,j) independent tracer data */
struct sed_tracer {
/* header */
  //2012 int prmspatial;
  int n;                        /* index */
  int n_hd_wc, n_hd_sed, n_file;
  char name[MAXSTRLEN];
  char units[MAXSTRLEN];
/* ranges and default values */
  double fill_value_wc;         /* watercolumn fill value */
  double fill_value_sed;        /* sediment fill value */
  int diagn;                    /* flag non-zero if diagnostic; default 0.
                                   if non-zero, 1 if "flux", other if
                                   "value". */
  int inwc;                     /* flag 1 if exists in water; default 1 */
  int insed;                    /* flag 1 if exists in sediments; default
                                   1 */
  int dissol;                   /* flag 1 if dissolved; default !partic
                                   (if present) or 1 */
  int partic;                   /* flag 1 if particulate; default !dissol
                                   (if present) or 0 */
  int resuspend;
  int deposit;

  int calcvol_wc;               /* flag 1 if volumetric; default 0 */
  int calcvol_sed;              /* flag 1 if volumetric; default 0 */
  int adsorb;                   /* flag 1 if adsorbed; default !adsorbed */
  int dissolvednum;             /* number of dissolved fraction for
                                   adsorbed tracer */
  int carriernum;               /* number of sediment, carrying adsorbed
                                   tracer */
  double adsorbkd;              /* sorption-desorption Kd [cub.m/kg] */
  double adsorbrate;            /* adsorption rate [day] */
  int diffuse;                  /* flag 1 if tracer to be diffused;
                                   default 1 */
  int cohesive;                 /* flag 1 if tracer is cohesive; default
                                   cohesive */
  int floc;                     /* flag 1 if tracer flocculates, default 0
                                 */
  double decay_days;                 /* decay e-folding time in days */
/* particulate tracers only (dissol = 0) */
  double psize;                 /* particle size - optional */
  double b_dens;                /* particle bulk density - mandatory
                                   [kg/cub.b] */
  double i_conc;                /* initial deposit concentration -
                                   mandatory [kg/cub.m] */
//2019
  double css_erosion;           /* tracer specific css erosion */
  double css_deposition;

  /* Column based, in case it varies spatially */
  double *svel;                 /* settling velocity - mandatory [m/s] */

  double mass_start_wc, mass_start_sed,mass_end_wc, mass_end_sed; /*mass balance check before and after the sediment cycle*/

  /*UR-ADDED optimization -
   * the scale applied to produce mg m-3 as unit for this tracer  */
  double u_scale;
};


/**** (i,j) dependent data ****/
struct sed_column {
  int i;                        /* Column index */
  int j;                        /* Column index */
  int col_number;               /* column number (1<=col_number<=ncol) */
  int hd_col_index;
  int sed_start;                /*=1 at the beginning of the sed cycle */

  /* these parameters must be specified as an input data at each time step
   */
  double area, *tss_wc,*tss_sed;             /* cell area [m*m] */
  double *gridz_sed, *gridz_wc; /* Layer interface z-coords [m] */
  double **tr_sed, **tr_wc;     /* concentrations [kg/cub.m] */
  double ***ptr_sed, ***ptr_wc; /* pointer tp concentrations [kg/cub.m] */
  double *u1_wc;                /* horizontal velocity [m/s] */
  double *u2_wc;                /* horizontal velocity [m/s] */
  double *Kz_wc;                /* vertical diffusion cfft in watcol
                                   [m*m/s] */
  double *diss_wc;              /* tke dissipation rate */
  double wind1, wind2;          /* wind stress components [N/m*m] */
  double z0_skin;               /* grain roughness [m] */
  double hripples, lripples;    /* ripple height and length [m] */
  double theta;                 /* angle between e1 and x axes [rad] */
  double dHde1, dHde2;          /* Slope of bottom at cell centres in e1 /
                                   e2 dirn [rad] */
  double wave_ub, wave_period, wave_dir;  /* wave data [m/s, s, grad] */
  double *erdepflux_ac;         /* erosion/deposition flux kg m-2 */
  double *cbnm1, *cbnm2, *cbnm3, *cbnm4, *cbfilt; /* HF fileted
                                                     concentrations in
                                                     active layer
                                                     [kg/cub.m] */
  /* can be recalculated at the beginning of each time step */
  int topk_sed, topk_wc;        /* Index of top layer */
  int botk_sed, botk_wc;        /* Index of bottom layer */
  double topz_sed, topz_wc;     /* Top layer interface [m] */
  double botz_sed, botz_wc;     /* Bottom layer interface [m] */
  double depth_sedi, depth_wc;   /* Sediment and water depths [m] */
  double *dz_sed, *dz_wc;       /* Layer thicknesses [m] */
  double *dzface_sed, *dzface_wc;/* dz values at interfaces (recalculated every timestep) */
  double *cellz_sed, *cellz_wc; /* Layer centre z-coords [m] */
  double *sigma_sed;            /* Sigma interfaces in sediments */
  double *dissol_kz;            /* Diffusion cfft for dissolved things in
                                   sedi [m*m/s] */
  double *partic_kz;            /* Diffusion cfft for particulate things
                                   in sedi [m*m/s] */

  double dissol_kz_i;           /* Sediment-water bioirrig exchange cfft [m*m/s] */
  double partic_kz_i;

  double *por_sed, *por_wc;     /* Porosity */
  double windspeed;             /* windspeed [m/s] */
  double z0_phys;               /* physical roughness [m] */


  /* selfrenewable at the beginning of each time step or used for 2-3 d
     visualisation */
  double *dzold_sed, *dzold_wc; /* Old layer thicknesses [m] */
  double *porold_sed, *porold_wc; /* Old porosity */
  double *svel_consolid, *svel_floc;  /* Velocity of consolidating
                                         particles or flocs [m/s] */
  double **svel_wc;
  double *gridvel_sed, *gridvel_wc; /* Grid interface velocity [m/s] */
  double *watvel_sed, *watvel_wc; /* Water vertical velocity [m/s] */
  double *css;                  /* Critical shear stress for erosion
                                   [N/m*m] */
  double css_dep;               /* Critical shear stress for deposition
                                   [N/m*m] */
  double ustrcw, ustrcw_skin;   /* Total and skin friction velocity [m/s] */
  double Cd;                    /* Drag cffnt */
  double dzactive;              /* Thickness of an active layer [m] */
  double *coh_wc, *coh_sed;     /* Cohesive sediment content in % */
  double *erflux;               /* Erosion flux kg m-2 s-1 (ntr) */
  double *depflux;              /* Deposition flux kg m-2 s-1 (ntr) */
  double *erdepflux;            /* Erosion+Deposition kg m-2 s-1 (ntr) */
  double erdepflux_total;       /* flux of sediment particles */
  double erdepflux_oxygen;
  double *bedloadflux;          /* Bedload flux kg m-2 s-1 (ntr) */
  double *ref_c_eq;             /* Equilibrium reference concebtration kg
                                   m-3 (ntr) */
  double *ref_c;                /* Actual concentration in reference point
                                   kg m-3 (ntr) */
  double temp1, temp2, temp3;   /* tmp values */
  double *tmp_sed;              /* tmp array (sednz+1) NMY*/
  double *tmp1000a, *tmp1000b;              /* tmp array (nz+1) NMY*/
  double *tmp1000c;

  double hardsub_scale; /* hard substrate scaling coeff*/

  double *tr_srf_flux;  /* Array of surface tracer fluxes */
};

struct sed_spatial {

  double *theta;
  double *hripples;
  double *lripples;

  double **cbnm1;               /* HF filtered sed. concentration at (n-1)
                                   time step [kg m-3] */
  double **cbnm2;               /* HF filtered sed. concentration at (n-2)
                                   time step [kg m-3] */
  double **cbnm3;               /* HF filtered sed. concentration at (n-3)
                                   time step [kg m-3] */
  double **cbnm4;               /* HF filtered sed. concentration at (n-4)
                                   time step [kg m-3] */
  double **erdeprate;
};

#endif




