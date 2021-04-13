/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/tracerstats/include/tracerstats.h
 *  
 *  Description:
 *  Tracerstats header file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracerstats.h 6762 2021-04-13 02:12:03Z riz008 $
 *
 */

/* include stdio for all files */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* emslib related stuff messages*/
#include "ems.h"

/* Release verions and getters */
#define TRACERSTATS_MAJOR_VERSION 1
#define TRACERSTATS_MINOR_VERSION 0
#define TRACERSTATS_PATCH_VERSION 1

int get_tracerstats_major_vers(void);
int get_tracerstats_minor_vers(void);
int get_tracerstats_patch_vers(void);

#if !defined(MINVAL)
#define MINVAL (1.e-13)
#endif

#define NOT   -0x000001
#define ETA   -0x000002
#define KZ    -0x000004
#define VZ    -0x000008
#define U1VH  -0x000010
#define U2VH  -0x000020
#define U1KH  -0x000040
#define U2KH  -0x000080

#define MEAN   0x000001
#define VINT   0x000002
#define DIFF   0x000004
#define STDEV  0x000008
#define VAR    0x000010
#define COV    0x000020
#define CORR   0x000040
#define FLUXE1 0x000080
#define FLUXE2 0x000100
#define FLUXW  0x000200
#define MEANFLUXE1 0x000400
#define MEANFLUXE2 0x000800
#define MEANFLUXW  0x002000
#define SECTFLUX   0x001000
#define COPY       0x004000
#define SUM        0x008000
#define TRVDIFF    0x010000
#define VMEAN      0x020000
#define SED        0x800000
#define MAX        0x040000
#define MIN        0x080000
#define VMAX       0x100000
#define VMIN       0x200000
#define DHW        0x400000
#define SUM_DT     0x800000
#define RMSE       0x1000000
#define EXPOSURE   0x2000000
#define COPYLAYER  0x4000000
#define REEFTEMP   0x8000000
#define RUN_MEAN   0x10000000
#define VPROF      0x20000000

#define SSECTIONV "w"
#define SSECTION1 "u1"
#define SSECTION2 "u2"

#define SECTIONW  0
#define SECTIONE1 1
#define SECTIONE2 2

#define VDIFFMODNORMAL  0
#define VDIFFMODVOL     1
#define VDIFFMODDIV     2
#define VDIFFMODVOLDIV  3
#define VDIFFSEDFLUX1   4
#define VDIFFSEDFLUX2   5
#define VDIFFSEDFLUX3   6

#define EX_GT  2  /* Greater than exposure threshold flag */
#define EX_LT  4  /* Less than exposure threshold flag */
#define EX_TR  8  /* Exposure threshold is a tracer */
#define EX_VA  16 /* Exposure threshold is a value */

#define MAXSTEPS 10
#define LASTSTEP 0

struct trs;
typedef struct trs trs_t;
/* Define the struct prototypes for 3d and 2d stat functions
 */
struct fluxfn;
typedef struct fluxfn fluxfn_t;

struct statfn2d;
typedef struct statfn2d statfn2d_t;

struct lstatfn2d;
typedef struct lstatfn2d lstatfn2d_t;

struct statfn3d;
typedef struct statfn3d statfn3d_t;

struct lstatfn3d;
typedef struct lstatfn3d lstatfn3d_t;

struct domainfn3d;
typedef struct domainfn3d domainfn3d_t;

/* function specific structs
 * */
struct int_dummy;
typedef struct int_dummy int_dummy_t;

struct sectioncoord;
typedef struct sectioncoord sectioncoord_t;

struct vdiff_range;
typedef struct vdiff_range vdiff_range_t;

struct dependency;
typedef struct dependency dependency_t;

/*
 * Note: The following working arrays are now also used for the tracer
 *       max & min functions. See the respective functions in
 *       statistics.c for details. -FR
 *
 * These have been separted here in order to distinguish between 
 * water column and sediments 
 */
typedef struct {
  double *w1; /* 3d mean working arrays, contains the dt over whichto average */
  double *w2; /* 3d mean working arrays, runnning total of dt sofar */
  double *w3; /* 3d mean working arrays, */
  int    *w4; /* 3d mean working arrays, flag if a mean has been reset */
  int    *w5; /* 3d working arrays */
  double *w6; /* 3d working array */
} stat3d_work_arrays;


/* Define the function prototypes for 3d and 2d stt functions
 */
typedef void (*tracer_statfn2d) (trs_t* t, int* n, void* data);
typedef void (*tracer_statfn3d) (trs_t* t, int* n, void* data, int c);
typedef void (*tracer_statfn3d_new) (trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);

typedef void (*tracer_fluxfn) (trs_t* t, double *flux, void* data, int* n, int c);
typedef void (*tracer_domainfn3d) (void*model, trs_t* t, int* n, void* data, int c);

/* function prototypes for custom step functions */
typedef void* (*tracer_domain_create)(void* d);
typedef void (*tracer_domain_scatter)(int n, void* o, void* t);
typedef void (*tracer_domain_gather)(int n, void* t, void* o);
typedef void* (*tracer_domain_destroy)(void* d);

struct lstatfn2d{
  statfn2d_t* statfn;
  int  nfn;
};

struct statfn2d{
  tracer_statfn2d statfn;
  int  n;
  void* data;
  int step;
  int type;
};


struct lstatfn3d{
  statfn3d_t* statfn;
  int  nfn;
};

struct statfn3d{
  tracer_statfn3d statfn;
  tracer_statfn3d_new statfn_new;
  int  n;
  void* data;
  int step;
  int type;
};

struct fluxfn{
  tracer_fluxfn fluxfn;
  double* flux;
  int n;
  int step;
  void* data;
  int type;
};

struct domainfn3d{
  tracer_domainfn3d domainfn;
  int  n;
  void* data;
  int step;
  int type;
};

struct int_dummy{
  int  d;
};

#define SECTION_MODE_PAD 5

typedef enum {
  SECTION_INTEGRATE = 0,
  SECTION_AVERAGE
} SECTION_MODE;

struct sectioncoord {
  int** data;
  int  nc;
  int  ntrs;
  int* trsn;
  int  dir;
  char* fileout;
  char** tnames;
  char* name;
  char** format;
  double dt;
  double tscale;
  double* outscales;
  char* tunit;
  char* outunit;
  double lastdump;
  double startt;
  double time;
  double* results; /* 
		    * Original concentration flux integrated (mode == 0) 
		    * for mode == 1, this is averaged over space & time
		    */
  double* results_c;   // average concentration
  double* results_sq;  // results squared
  double* results_cusgn_a;  // u-c covariance, part a
  double* results_cusgn_b;  // u-c covariance, part b
  double  u;       /* mean of the normal velocity */
  double  uabs;    /* mean of the absolute value of the normal velocity */
  double  usq;     /* mean of the velocity squared */

  double  Tot_area;    /* total area of the corresponding face */
  int nsteps;
  SECTION_MODE mode; // whether to average or integrate
};

struct vdiff_range{
  int* toprange;
  int* botrange;
  int  nt;
  int  nb;
  int  ntr3d;
  int  strict;
  int  mode;
  char* pairname;
};



/* define the depencies between functions
 */
struct dependency{
  int func;
  int* depends;
};

struct trs {
  int do_tracerstats;
  int cols;
  int nz;
  int sednz;
  int ntr;
  int ntrS;
  int atr;
  int atrS;
  int atrsed;
  int nsed;
  int *e_nstep;
  int *o_nstep;
  int *news;
  int new_step;
  int use_eta;
  int use_etamean;
  int use_Kz;
  int use_Vz;
  int use_fluxe1;
  int use_fluxe2;
  int use_w;
  int use_u1;
  int use_u2;
  int use_area_e1_e2;
  int *nprestep;
  int *nstep;
  int *stat_type;
  int *stat_type_sed;
  int *stat_typeS;
  /*
   * UR define an array of functions to be executed
   * where statfn_t is an array in which the first element is  executed first
   * than the fluxfn_t array,
   * than the domainfn_t functions
   * and finally the elements 1 to nstatfn-1 in ascending order
   * */
  lstatfn3d_t* statfn3d;
  lstatfn3d_t* statfn_sed;
  lstatfn2d_t* statfn2d;
  fluxfn_t* fluxfn;
  domainfn3d_t* domainfn3d;
  int nstatfn2d;
  int nstatfn3d;
  int nstatfn_sed;
  int nfluxfn;
  int ndomainfn3d;
  int **map_2d;
  int **map_3d;
  int **map_sed;
  int topk_wc;
  int botk_wc;
  int topk_sed;
  int botk_sed;
  double dt;
  double time;
  double depth_wc;
  double eta;
  double *dz_wc;
  double *dz_sed;
  double *Kz;
  double *Vz;
  double *u1;
  double *u2;
  double *w;
  double *u1flux3d;
  double *u2flux3d;
  double *area_e1;
  double *area_e2;
  double area_w;
  double h1au2;
  double h2au1;
  double ***tr_wc; /* these are pointers to the host model memory */
  double **tr_e1wc; /* these are derived values and not meant for assignment */
  double **tr_e2wc; /* these are derived values and not meant for assignment */
  double **tr_in;
  double ***tr_sed;
  char **trname_3d;
  char **trname_2d;
  char **trname_sed;
  int *tmap_3d;
  int *tmap_2d;
  int *tmap_sed;
  double d1;
  double d2;
  double *w1S;	/* 2d mean working arrays, */
  double *w2S;	/* 2d mean working arrays, */
  double *w3S;	/* 2d mean working arrays, */

  stat3d_work_arrays w_wc;
  stat3d_work_arrays w_sed;

  void* model;
};

/* Internal routines */
trs_t* trs_create();
trs_t* trs_build(void* model, FILE *fp);
void trs_init(void* model, trs_t *inter, FILE *fp);
void trs_step(void* model, trs_t *inter, int c,int step);

/* 3d stat functions, functions requiring the tracer index n,
 * custom data data and the column identifer c
 */
// The new definitions
void tracer_mean_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_run_mean_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_rmse_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_max_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_min_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_sum_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_stdev_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);
void tracer_var_3d(trs_t *trs, stat3d_work_arrays *w_arr, double ***tr_vals, int n, int botk, int topk, int **map, void* data);

// Still the old ones
void tracer_diff_3d(trs_t *trs, int* n, void* data, int c);
void tracer_cov_3d(trs_t *trs, int* n, void* data, int c);
void tracer_corr_3d_pre(trs_t *trs, int* n, void* data, int c);
void tracer_corr_3d_post(trs_t *trs, int* n, void* data, int c);
void tracer_flux3d_e1(trs_t *trs, int* n, void* data, int c);
void tracer_flux3d_e2(trs_t *trs, int* n, void* data, int c);
void tracer_flux3d_w(trs_t *trs, int* n, void* data, int c);
void tracer_dhw(trs_t *trs, int* n, void* data, int c);
void tracer_sum_dt_3d(trs_t *trs, int* n, void* data, int c);
void tracer_exposure(trs_t *trs, int* n, void* data, int c);
void tracer_reeftemp(trs_t *trs, int* n, void* data, int c);
void tracer_vprof(trs_t *trs, int* n, void* data, int c);

/* 3d stat functions, functions requiring the tracer index 'n',
 * the flux 'flux' and the column identifer 'c'
 */
void tracer_meanflux_3d(trs_t *trs, double *flux, void* data, int* n, int c);

/* domain functions, functions which require knowledge about more than just
 * the column they operate on.
 */
void tracer_section(void* model, trs_t *trs, int* n, void* data, int c);


/* 2d stat functions, functions requiring the tracer index 'n',
 */
void tracer_mean_2d(trs_t *trs, int* n, void* data);
void tracer_run_mean_2d(trs_t *trs, int* n, void* data);
void tracer_rmse_2d(trs_t *trs, int* n, void* data);
void tracer_max_2d(trs_t *trs, int* n, void* data);
void tracer_min_2d(trs_t *trs, int* n, void* data);
void tracer_corr_2d_pre(trs_t *trs, int* n, void* data);
void tracer_corr_2d_post(trs_t *trs, int* n, void* data);
void tracer_cov_2d(trs_t *trs, int* n, void* data);
void tracer_diff_2d(trs_t *trs, int* n, void* data);
void tracer_sum_2d(trs_t *trs, int* n, void* data);
void tracer_stdev_2d(trs_t *trs, int* n, void* data);
void tracer_var_2d(trs_t *trs, int* n, void* data);

void vertical_max(trs_t *trs, int* n, void* data);
void vertical_integrals(trs_t *trs, int* n, void* data);
void vertical_mean(trs_t *trs, int* n, void* data);
void copy_layer(trs_t *trs, int* n, void* data);
void tracer_vertical_strat(trs_t *trs, int* n, void* data);
void tracer_vertical_strat_vol(trs_t *trs, int* n, void* data);
void tracer_vertical_strat_div(trs_t *trs, int* n, void* data);
void tracer_vertical_strat_voldiv(trs_t *trs, int* n, void* data);
void tracer_vertical_sedflux1(trs_t *trs, int* n, void* data);
void tracer_vertical_sedflux2(trs_t *trs, int* n, void* data);


/* Interface routines */
extern double i_get_model_time(void* hmodel);
extern char *i_get_model_timeunit(void* hmodel);
extern double i_get_model_timestep(void* hmodel);
extern int i_is_valid_tracer(int itr);
extern int i_is_valid_sed_tracer(int itr);
extern int  i_get_num_wclayers(void* hmodel);
extern int  i_get_num_columns(void* hmodel);
extern int  i_get_nstep(void* hmodel);
extern int  i_get_num_sedlayers(void* hmodel);
extern int i_get_num_tracers_3d(void* hmodel, int *atr);
extern int i_get_num_tracers_2d(void* hmodel, int *atrS);
extern int i_get_num_sedtracers(void* hmodel);
extern void i_get_names_tracers_3d(void* hmodel, char **trname);
extern void i_get_names_tracers_2d(void* hmodel, char **trname);
extern void i_get_names_tracers_sed(void* hmodel, char **trname);
extern int i_get_host_counter(void* hmodel, int c);
extern int i_get_interface_counter(void* hmodel, int c);
extern double i_get_cellarea_w(void* hmodel, int c);
extern void i_get_cellarea_e1(void* hmodel, int c, double *area);
extern void i_get_cellarea_e2(void* hmodel, int c, double *area);
extern int i_get_topk_wc(void* hmodel, int c);
extern int i_get_botk_wc(void* hmodel, int c);
extern int i_get_topk_sed(void* hmodel, int c);
extern int i_get_botk_sed(void* hmodel, int c);
extern double i_get_botz_wc(void* hmodel, int c);
extern double i_get_topz_wc(void* hmodel, int c);
extern double i_get_depth(void* hmodel, int c);
extern double i_get_eta(void* hmodel, int c);
extern void i_get_fluxe1_wc(void* hmodel, int c, double *u1flux3d);
extern void i_get_fluxe2_wc(void* hmodel, int c, double *u2flux3d);
extern void i_get_u1_wc(void* hmodel, int c, double *u1);
extern void i_get_u2_wc(void* hmodel, int c, double *u2);
extern void i_get_w_wc(void* hmodel, int c, double *w);
extern void i_get_dz_wc(void* hmodel, int c, double *dz_wc);
extern void i_get_dz_sed(void* hmodel, int c, double *dz_sed);
extern void i_get_Kz_wc(void* hmodel, int c, double *Kz_wc);
extern void i_get_Vz_wc(void* hmodel, int c, double *Vz_wc);
extern void i_get_gridz_sed(void* hmodel, int c, double *gridz_sed);
extern void i_get_tracer_wc(void* hmodel, int c, int ntr, int *tmap, double ***tr_wc);
extern void i_get_tracer_e1(void* hmodel, int c, int ntr, int *tmap, double **tr_wc);
extern void i_get_tracer_e2(void* hmodel, int c, int ntr, int *tmap, double **tr_wc);
extern void i_get_tracer_sed(void* hmodel, int c,  int ntr, int *tmap, double ***tr_sed);
extern void i_get_tracer_2d(void* hmodel, int c, int ntr, int *tmap, double **tr_in);
extern int *i_get_tmap_3d(void* hmodel, int ntr, char *trname[]);
extern int *i_get_tmap_sed(void* hmodel, int ntr, char *trname[]);
extern int *i_get_tmap_2d(void* hmodel, int ntr, char *trname[]);
extern int i_get_param_map_2d(void* hmodel, int n);
extern int i_get_param_map_3d(void* hmodel, int n);
extern int i_get_param_map_sed(void* hmodel, int n);
extern int i_tracername_exists(void* model, char*name);
extern int i_tracername_exists_sed(void* model, char*name);

extern int i_is_ijk(void* hmodel, int col, int* ij, int do_k);
extern int i_in_window(void* hmodel, int col, int* i, int* j,int nn);
extern int* i_get_ijk(void* hmodel, int col, int* ij);
extern tracer_info_t* i_get_tracer(void* hmodel, char* trname);
extern int i_get_step(void* model, char* st);
extern void trstat_remove_domain_function(char* name);
extern void trstat_add_domain_function(char* name, int n,int scattertype, tracer_domain_create create, tracer_domain_gather gather, tracer_domain_scatter scatter);/*, tracer_domain_destroy destroy);*/
extern void get_section_area(sectioncoord_t* section, int topk_wc, int botk_wc);
extern void i_set_error(void* hmodel, int col, int errorf, char *text);



