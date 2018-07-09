/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/waves/waves.h
 *  
 *  Description:
 *  Include file for waves library.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: waves.h 5849 2018-06-29 05:03:57Z riz008 $
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

#define VON_KAR 0.408
#define TOBA    0x020   /* Toba (1978) wind waves                    */
#define USAC    0x040   /* US Army Core Engineers wind waves         */
#define WNONE   0x000   /* No action taken                           */
#define WFILE   0x020   /* Wave variables read from file             */
#define WCOMP   0x040   /* Wave variables computed                   */
#define WWIND   0x080   /* Wave variables = wind waves               */

struct wave;
typedef struct wave wave_t;

/* Names of derived quantities in the host model */

struct wave{
  int do_waves;        /* Use the wave module */
  int do_orbital;      /* Wave parameters read from file OK */
  int do_amp;          /* Wave amplitude read from file OK */
  int do_per;          /* Wave period read from file OK */
  int do_dir;          /* Wave direction read from file OK */
  int do_ub;           /* Wave orbital velocity read from file OK */
  int do_wif;          /* Do wave-induced forcing read OK from file */
  int do_wcf;          /* Do wave current friction */
  int do_rs;           /* Do radiation stresses */
  int do_Cd;           /* Do bottom drag modification */
  int do_fetch;        /* Compute fetch */
  int do_stokes;       /* Stokes velocity read from file OK */
  int windwave;        /* Wind wave method */
  int c, cc, ij[3];
  int cols;
  int dt_ratio;
  int nz;
  int sednz;
  int ntr;
  int ntrS;
  int atr;
  int atrS;
  int atrsed;
  int nsed;
  int topk_wc;
  int botk_wc;
  int topk_sed;
  int botk_sed;
  double dt;
  double time;
  double top;
  double bot;
  double bot_l;
  double depth_wc;
  double *dz_wc;
  double *dz_sed;
  double *Kz;
  double *Vz;
  double *u1;
  double *u2;
  double *w;
  double u1bot;
  double u2bot;
  double eta;
  double z0;
  double quad_bfc;
  double h1au2;
  double h2au1;
  double *thetau1;
  double *thetau2;
  double *sinthcell;
  double *costhcell;
  double area;
  double **tr_in;
  char **trname_3d;
  char **trname_2d;
  char **trname_sed;
  int *tmap_3d;
  int *tmap_2d;
  int *tmap_sed;
  int *brsm;
  double **fetch;
  double *amp;
  double *dir;
  double *period;
  double *ub;
  double *Sxy;
  double *Syx;
  double *Fx;
  double *Fy;
  double *ustrcw;
  double *Cd;
  double *ste1;
  double *ste2;
  double wx;
  double wy;
  char amp_name[MAXSTRLEN];
  char period_name[MAXSTRLEN];
  char dir_name[MAXSTRLEN];
  char ub_name[MAXSTRLEN];
  char Sxy_name[MAXSTRLEN];
  char Syx_name[MAXSTRLEN];
  char Fx_name[MAXSTRLEN];
  char Fy_name[MAXSTRLEN];
  char ustrcw_name[MAXSTRLEN];
  char Cd_name[MAXSTRLEN];
  char ste1_name[MAXSTRLEN];
  char ste2_name[MAXSTRLEN];
  int ampid;
  int perid;
  int dirid;
  int ubid;
  int Sxyid;
  int Syxid;
  int Fxid;
  int Fyid;
  int ustid;
  int Cdid;
  int ste1id;
  int ste2id;
  double d1;
  double d2;

  void* model;
} ;

/* Release verions and getters */
#define WAVES_MAJOR_VERSION 1
#define WAVES_MINOR_VERSION 0
#define WAVES_PATCH_VERSION 0

int get_waves_major_vers(void);
int get_waves_minor_vers(void);
int get_waves_patch_vers(void);

/* Internal routines */
wave_t* wave_create();
wave_t* wave_build(void* model, FILE *fp);
void wave_init(void* model, wave_t *wave, FILE *fp);
void wave_step(void* model, wave_t *wave, int c);
void wave_run_setup(FILE *fp, wave_t *wave);

double wavedir_estim(double wy, double wx);
double wavenumb_estim(double period, double depth);
double waveub_estim(double wamp, double period, double wavenumber,
                    double depth);
double wamp_estim(double wx, double wy, double wavenumber,
                  double period);
void wavecomp_estim(double a, double dir, double *wx, double *wy);
void radiation_stress(wave_t *wave, int c);
void ustrcw_estim(wave_t *wave, int c);

/* Generic interface routines */
extern int* i_get_ijk(void* hmodel, int col, int* ij);
extern double i_get_model_time(void* hmodel);
extern double i_get_model_timestep(void* hmodel);
extern int i_is_valid_tracer(int itr);
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
extern int i_get_interface_counter(void* hmodel, int c);
extern double i_get_cellarea(void* hmodel, int c);
extern double i_get_cellarea_w(void* hmodel, int c);
extern int i_get_topk_wc(void* hmodel, int c);
extern int i_get_botk_wc(void* hmodel, int c);
extern int i_get_topk_sed(void* hmodel, int c);
extern int i_get_botk_sed(void* hmodel, int c);
extern double i_get_botz_wc(void* hmodel, int c);
extern double i_get_topz_wc(void* hmodel, int c);
extern double i_get_depth(void* hmodel, int c);
extern double i_get_eta(void* hmodel, int c);
extern void i_get_u1_wc(void* hmodel, int c, double *u1);
extern void i_get_u2_wc(void* hmodel, int c, double *u2);
extern void i_get_w_wc(void* hmodel, int c, double *w);
extern void i_get_dz_wc(void* hmodel, int c, double *dz_wc);
extern void i_get_dz_sed(void* hmodel, int c, double *dz_sed);
extern void i_get_Kz_wc(void* hmodel, int c, double *Kz_wc);
extern void i_get_Vz_wc(void* hmodel, int c, double *Vz_wc);
extern void i_get_gridz_sed(void* hmodel, int c, double *gridz_sed);
extern void i_get_tracer_wc(void* hmodel, int c, int ntr, int *tmap, double ***tr_wc);
extern void i_get_tracer_2d(void* hmodel, int c, int ntr, int *tmap, double **tr_in);
extern int *i_get_tmap_3d(void* hmodel, int ntr, char *trname[]);
extern int *i_get_tmap_sed(void* hmodel, int ntr, char *trname[]);
extern int *i_get_tmap_2d(void* hmodel, int ntr, char *trname[]);

/* Wave interface routines */
extern int w_check_orbital_file(void *hmodel);
extern void w_check_wave_data(int *wa, int *wp, int *wd, int *wu, int *rs, int *stv);
extern double w_get_dt(void *hmodel);
extern int w_check_wave_period(void *hmodel);
extern double w_get_wave_period(void *hmodel, int c);
extern int w_check_wave_amp(void *hmodel);
extern double w_get_wave_amp(void *hmodel, int c);
extern int w_check_wave_dir(void *hmodel);
extern double w_get_wave_dir(void *hmodel, int c);
extern int w_check_wave_ub(void *hmodel);
extern double w_get_wave_ub(void *hmodel, int c);
extern int w_check_wave_Fx(void *hmodel);
extern double w_get_wave_Fx(void *hmodel, int c);
extern int w_check_wave_Fy(void *hmodel);
extern double w_get_wave_Fy(void *hmodel, int c);
extern int w_check_wind(void *hmodel);
extern double w_get_wave_wind1(void *hmodel, int c);
extern double w_get_wave_wind2(void *hmodel, int c);
extern double w_get_thetau1(void *hmodel, int c);
extern double w_get_thetau2(void *hmodel, int c);
extern int w_get_winsize(void *hmodel);
extern void w_get_brsm(void *hmodel, int *brsm);
extern void w_get_bot_vel(void *hmodel, double sinthcell, double costhcell,
			  double *u1bot, double *u2bot, double *botz, int c);
extern double w_get_quad_bfc(void *hmodel);
extern double w_get_z0(void *hmodel, int c); 
extern double w_get_sinthcell(void *hmodel, int cc);
extern double w_get_costhcell(void *hmodel, int cc);
extern double w_get_fetch(void *hmodel, int cc, int n);
extern double w_check_fetch(void *hmodel);
extern void i_set_error(void* hmodel, int col, int errorf, char *text);
extern int w_check_wave_ste1(void *hmodel);
extern double w_get_wave_ste1(void *hmodel, int c);
extern int w_check_wave_ste2(void *hmodel);
extern double w_get_wave_ste2(void *hmodel, int c);
