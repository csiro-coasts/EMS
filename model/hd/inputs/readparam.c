/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/inputs/readparam.c
 *  
 *  Description:
 *  Routine to read model parameters
 *  from a file.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: readparam.c 5901 2018-08-28 02:10:22Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0

#define MAXBDRY    200
#define MAXDEPTH   9000
#define LAYERMIN   0.5

#define RECT       1
#define GEO_RECT   2
#define FPOLE      4
#define CART       8

#define N_BDRY    1
#define S_BDRY    2
#define E_BDRY    4
#define W_BDRY    8

int nce1;
int nce2;

void geog_false_pole_coord(double **x, double **y, double **h1,
                           double **h2, double **a1, double **a2,
                           long int nce1, long int nce2, double x00,
                           double y00, double flon, double flat,
                           double xinc, double yinc);
void geog_dreckon_coord(double **x, double **y, double **h1, double **h2,
                        double **a1, double **a2, long int nce1,
                        long int nce2, double x00, double y00, double rotn,
                        double xinc, double yinc);
void f_pole(double ilat, double ilon, double ang, double *olat,
            double *olon);
static int read_bathy_from_file(FILE *fp, parameters_t *params);
static double *read_bathy_from_db(FILE *fp, parameters_t *params);
static void mask_bathy_with_coast(const char *cfname, parameters_t *params);
void reset_bathy(parameters_t *params, int io, int jo, double val);
int is_wet(parameters_t *params, double bathy);
int is_land(parameters_t *params, double bathy);
int is_outside(parameters_t *params, double bathy);
int read2darray(FILE * fp, char *label, double **array, long n1, long n2);
void shift_grid(double **x, double **y, double **h1, double **h2,
                double **a1, double **a2, int nce1, int nce2, int sft);
void eta_init(geometry_t *geom, parameters_t *params, master_t *master);
void neighbour_none(dump_data_t *dumpdata);
void read_zoom(parameters_t *params, FILE *fp);
void read_eta_relax(parameters_t *params, FILE *fp);
void read_vel_relax(parameters_t *params, FILE *fp);
void read_win_type(parameters_t *params, FILE *fp);
void check_TS_relax(parameters_t *params, FILE *fp);
void baro_vortex(master_t *master, geometry_t *geom, int io, int jo);
int check_river_bdry(parameters_t *params, FILE *fp);
int check_bdry_options(parameters_t *params, FILE *fp);
char *bdryname(int m);
char *mixname(int m);
char *cflname(int m);
char *substepname(int m);
char *heatfluxname(int m);
char *meansname(int m, char *buf);
char *btname(int m);
char *adname(int m);
char *tf(int m);
void write_grid_specs(FILE *op, parameters_t *params);


/*-------------------------------------------------------------------*/
/* Routine to set default parameters                                 */
/*-------------------------------------------------------------------*/
void set_default_param(parameters_t *params)
{
  sprintf(params->codeheader, "%c", '\0');
  sprintf(params->parameterheader, "%c", '\0');
  sprintf(params->opath, "%c", '\0');
  sprintf(params->trkey, "%c", '\0');
  sprintf(params->sequence, "%c", '\0');
  sprintf(params->bdrypath, "%c", '\0');
  sprintf(params->tracerdata, "%c", '\0');
  sprintf(params->rivldir, "%c", '\0');
  sprintf(params->trans_dt, "%c", '\0');
  params->runno = 0;
  sprintf(params->rev, "%c", '\0');
  params->momsc = ORDER2;
  params->trasc = VANLEER;
  params->ultimate = 0;
  params->osl = 0;
  params->smagorinsky = 0.0;
  params->sue1 = params->sue2 = 0.0;
  params->kue1 = params->kue2 = 0.0;
  params->bsue1 = params->bsue2 = 0.0;
  params->bkue1 = params->bkue2 = 0.0;
  sprintf(params->smag, "%c", '\0');
  params->smag_smooth = 0;
  params->diff_scale = LINEAR;
  params->visc_method = LAPLACIAN;
  params->stab = NONE;
  params->atr = 0;
  params->ntr = 0;
  params->ntrS = 0;
  params->nsed = 0;
  params->sednz = 0;
  params->tratio = 1.0;
  params->trsplit = 0;
  /* params->sliding = 0; / *UR*/
  params->cfl = NONE;
  memset(params->cfl_dt, 0, sizeof(params->cfl_dt));
  params->mixlayer = NONE;
  params->show_layers = 0;
  params->means = NONE;
  sprintf(params->means_os, "%c", '\0');
  sprintf(params->bathystats, "%c", '\0');
  params->reset_bathyf = 0;
  params->data_infill = 0;
  params->vorticity = NONE;
  params->numbers = NONE;
  params->u1_f = params->u1av_f = params->u2_f = params->u2av_f = NONE;
  params->save_force = NONE;
  params->lnm = 0.0;
  strcpy(params->trflux, "NONE");
  strcpy(params->trperc, "NONE");
  sprintf(params->trpercr, "%c", '\0');
  params->trflsh = 0;
  sprintf(params->trage, "%c", '\0');
  params->dhwf = 0;
  sprintf(params->dhw, "%c", '\0');
  params->tendf = 0;
  sprintf(params->trtend, "%c", '\0');
  params->thin_merge = 1;
  params->sigma = 0;
  params->maxgrad = -1.0;
  params->bvel = 0.0;
  params->nonlinear = 1;
  params->calc_dens = 1;
  params->mode2d = 0;
  params->smooth = 0;
  params->slipprm = 1.0;
  params->nwindows = 1;
  params->win_reset = 0;
  params->win_type = STRIPE_E1;
  params->win_block = 0;
  params->win_size = NULL;
  params->nwn = NULL;
  params->wnx = NULL;
  params->wny = NULL;
  params->show_win = 0;
  params->restart_dt = 0.0;
  params->da = 0;
  params->da_dt = 0.0;
  params->da_fcst_dt = 0.0;
  params->prex = 0;
  params->ndf = 0;
  params->riverflow = 0;
  sprintf(params->restart_name, "%c", '\0');
  sprintf(params->win_file, "%c", '\0');
  sprintf(params->wind_file, "%c", '\0');
  sprintf(params->bdry_file, "%c", '\0');
  sprintf(params->ptinname, "%c", '\0');
  sprintf(params->swr_babs, "%c", '\0');
  sprintf(params->swr_attn, "%c", '\0');
  sprintf(params->swr_attn1, "%c", '\0');
  sprintf(params->swr_tran, "%c", '\0');
  sprintf(params->densname, "%c", '\0');
  sprintf(params->regions, "%c", '\0');
  sprintf(params->region_dt, "%c", '\0');
  sprintf(params->region_vars, "%c", '\0');
  sprintf(params->region_mode, "%c", '\0');
  params->region_obc = 0;
  params->nbl1 = params->nbl2 = 0;
  params->do_pt = 0;
  strcpy(params->dp_mode, "none");
  params->edgef = 0;
  params->waves = 0;
  params->decorr = 0;
  params->decf = NONE;
  params->orbital = 0;
  params->rampf = 0;
  params->wind_scale = 1.0;
  params->totals = 0;
  params->roammode = A_ROAM_R1;
  params->robust = 1;
  params->speed = 1;
  params->fatal = ETA_A | NANF;
  params->avhrr = 0;
  params->ghrsst = 0;
  params->dozoom = NOZOOM;
  params->zoomf = 1;
  params->zmfe1 = 1;
  params->zmfe2 = 1;
  params->exmapf = 0;
  params->heatflux = NONE;
  params->saltflux = NONE;
  params->water_type = NONE;
  params->tsfile_caching = 1;
  params->etamax = 10.0;
  params->velmax = 3.87;
  params->velmax2d = 2.55;
  params->etadiff = 0.0;
  params->wmax = 100.0;
  params->trsplit = 0;
  params->etarlx = NONE;
  params->velrlx = NONE;
  params->albedo = -9999.0;
  params->domom = 0;
  params->doroms = 0;
  params->zref = 10.0;
  params->neutral = 0;
  params->wave_alpha = 100.0;
  params->wave_hf = 1.0;
  params->wave_b1 = 0.0014;
  params->fillf = NONE;
  params->filter = NONE;
  params->trfilter = NONE;
  params->conserve = NONE;
  params->lyear = 0;
  params->do_closure = 0;
  params->rsalt = 0;
  params->rtemp = 0;
  params->smooth_VzKz = 0;
  params->albedo_l = -1;
  params->trout = 0;
  params->noutside = 0;
  params->porusplate = 0;
  params->sharp_pyc = 0;
  sprintf(params->reef_frac, "%c", '\0');
  params->dbi = params->dbj = params->dbk = -1;
  params->dbgf = NONE;
  params->dbgtime = 0.0;
  memset(params->momfile, 0, sizeof(params->momfile));
  memset(params->avhrr_path, 0, sizeof(params->avhrr_path));
  memset(params->ghrsst_path, 0, sizeof(params->ghrsst_path));
  params->rampf = WIND|FILEIN|CUSTOM|TIDALH|TIDALC|TIDEBC|FLATHR|ETA_RELAX;
  memset(params->trvars, 0, sizeof(params->trvars));
  sprintf(params->smooth_v, "%c", '\0');
  sprintf(params->scale_v, "%c", '\0');
  sprintf(params->orthoweights, "%c", '\0');
  sprintf(params->nodal_dir, "%c", '\0');
  sprintf(params->tide_con_file, "%c", '\0');
  sprintf(params->webf, "%c", '\0');
  sprintf(params->regulate, "%c", '\0');
  params->regulate_dt = 0.0;
  memset(params->alert, 0, sizeof(params->alert));
  memset(params->alertc, 0, sizeof(params->alertc));
  memset(params->alert_dt, 0, sizeof(params->alert_dt));
  /* Set which alerts are activated                                  */
  params->eta_f = 1;      /* Implement alert action on elevation     */
  params->vel2d_f = 1;    /* Implement alert action on 2d velocity   */
  params->vel3d_f = 1;    /* Implement alert action on 3d velocity   */
  params->wvel_f = 0;     /* Implement alert action on w velocity    */
  params->tend_f = 0;     /* Implement alert action on tendencies    */
  params->div2d_f = 0;    /* Implement alert action on 2d divergence */
  params->div3d_f = 0;    /* Implement alert action on 3d divergence */
  params->cfl_f = 0;      /* Implement alert action on CFL           */
  params->ts_f = 0;       /* Implement alert action on T/S rates     */
  params->shear_f = 0;    /* Implement alert action on shear         */
  params->hdiff_f = 1;    /* Perform enhanced horizontal viscosity   */
  params->tide_r = NONE;  /* Tide removal method                    */
  params->parray_inter_botz = 0; /* Interpolate botz in PARRAY files */
  params->compatible = 0; /* Backwards compatibility                 */
  params->gint_errfcn = LWARN;

  /* Set the alert thresholds. Note: velocities are set via          */
  /* params->velmax and params->wmax.                                */
  /* Statistically derived maximums :                                */
  /*
  params->amax = 0.30;
  params->hmax = 0.15;
  params->vmax = 0.35;
  params->btmax = 0.30;
  params->bcmax = 0.20;
  params->cmax = 0.04;
  params->detamax = 3.5e-3;
  params->dwmax = 1.8e-3;
  */
  params->amax = 0.30;         /* Maximum advection tendency         */
  params->hmax = 5.15;         /* Maximum horizontal diffusion tend. */
  params->vmax = 5.35;         /* Maximum vertical diffusion tend.   */
  params->btmax = 5.30;        /* Maximum barotropic pressure tend.  */
  params->bcmax = 5.20;        /* Maximum baroclinic pressure tend.  */
  params->cmax = 5.04;         /* Maximum coriolis tendency          */
  params->detamax = 3.5e-3;    /* Maximum 2D divergence              */
  params->dwmax = 1.8e-3;      /* Maximum 3D divergence              */
  params->dtmax = 1.0;         /* Maximum T rate of change           */
  params->dsmax = 1.0;         /* Maximum S rate of change           */
  params->smax = 2.03e-4;      /* Maximum shear, dv/dx (s-1)         */
  /*params->smax = 1.3e-4;*/       /* Maximum shear, dv/dx (s-1)         */


#if defined(HAVE_SEDIMENT_MODULE)
  params->do_sed = 0;
  sprintf(params->sed_vars, "%c", '\0');
  sprintf(params->sed_defs, "%c", '\0');
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  params->do_eco = 0;
  sprintf(params->eco_vars, "%c", '\0');
  sprintf(params->eco_defs, "%c", '\0');
#endif

#if defined(HAVE_WAVE_MODULE)
  params->do_wave = 0;
#endif

  /* ROMS */
  params->roms_grid[0] = '\0';
}

/* END set_default_param()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the input parameters                              */
/*-------------------------------------------------------------------*/
parameters_t *params_read(FILE *fp)
{
  parameters_t *params;
  char buf[MAXSTRLEN];
  char buf1[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  int m, n;
  int cc;                       /* Sparse counters */
  double d1;

  /* Allocate memory for the parameter data structure */
  params = params_alloc();

  /* Read in the name of the input file */
  params->prmfd = fp;
  prm_set_errfn(hd_quit);

  if (forced_restart) {
    if (prm_read_char(fp, "restart_name", buf))
      strcpy(params->idumpname, buf);
    else
      strcpy(params->idumpname, "restart.nc");
     params->t = get_restart_time(params->idumpname, schedule->units);
     sprintf(params->start_time, "%g days", params->t / 86400);
  } else {
     sprintf(keyword, "INPUT_FILE");
     prm_read_char(fp, keyword, params->idumpname);

     /* Read the start and stop time (used only for diagnostics) */
     sprintf(keyword, "START_TIME");
     prm_read_char(fp, keyword, params->start_time);
     prm_get_time_in_secs(fp, "START_TIME", &params->t);
  }

  sprintf(keyword, "STOP_TIME");
  prm_read_char(fp, keyword, params->stop_time);
  sprintf(keyword, "LENUNIT");
  prm_read_char(fp, keyword, params->lenunit);
  sprintf(keyword, "TIMEUNIT");
  prm_read_char(fp, keyword, params->timeunit);

  /* Near real-time definitions */
  if (nrt_restart) {
    /* Check the model TIMEUNIT is seconds */
    if (contains_token(schedule->units, "seconds") == NULL)
      hd_quit("get_nrt_time: TIMEUNIT must be 'seconds since ....'.\n");
    params->t = get_restart_time(params->idumpname, schedule->units);
    sprintf(params->start_time, "%g days", params->t / 86400.0);
    params->rampstart = params->rampend = params->t;
    tm_scale_to_secs(params->stop_time, &d1);
    sprintf(params->stop_time, "%g days", (d1 + params->t) / 86400.0);
  }

  /* Grid dimensions */
  sprintf(keyword, "NCE1");
  prm_read_int(fp, keyword, &params->nce1);
  params->nfe1 = params->nce1 + 1;
  sprintf(keyword, "NCE2");
  prm_read_int(fp, keyword, &params->nce2);
  params->nfe2 = params->nce2 + 1;
  sprintf(keyword, "LAYERFACES");
  prm_read_darray(fp, keyword, &params->layers, &params->nz);
  if (--params->nz < 1)
    hd_quit("Number of layers in vertical must be 1 or more\n");
  params->bathy = NULL;

  read_bathy_from_file(fp,params);
  
  sprintf(keyword, "GRIDTYPE");
  prm_read_char(fp, keyword, params->gridtype);

  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "DESCRIPTION");
  prm_read_char(fp, keyword, params->grid_desc);
  sprintf(keyword, "NAME");
  prm_read_char(fp, keyword, params->grid_name);
  /* gridhoriz(params,geom,fp); */

  /* Flags */
  /* Set defaults */
  params->runmode = MANUAL;
  set_default_param(params);
  sprintf(keyword, "ETAMAX");
  prm_read_double(fp, keyword, &params->etamax);
  sprintf(keyword, "MIN_CELL_THICKNESS");
  prm_read_char(fp, keyword, params->mct);

  /* Reset bathymetry */
  sprintf(keyword, "RESET_BATHY");
  if (prm_read_char(fp, keyword, buf))
    params->reset_bathyf = is_true(buf);

  /* Bathymetry statistics */
  sprintf(keyword, "BATHY_STATS");
  if (prm_read_char(fp, keyword, params->bathystats))
    params->ntrS += 4;

  /* Use cascade search on input data */
  sprintf(keyword, "DATA_INFILL");
  if (prm_read_char(fp, keyword, buf))
    params->data_infill = is_true(buf);

  /* User regulation */
  sprintf(keyword, "REGULATE_FILE");
  prm_read_char(fp, keyword, params->regulate);
  if (prm_read_char(fp, "REGULATE_DT", buf))
    tm_scale_to_secs(buf, &params->regulate_dt);

  /* ID number and revision */
  sprintf(keyword, "ID_NUMBER");
  prm_read_double(fp, keyword, &params->runno);
  sprintf(keyword, "REVISION");
  prm_read_char(fp, keyword, params->rev);
  /* Fatal instability checking */
  sprintf(keyword, "FATAL");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->fatal = NONE;
    } else {
      params->fatal = 0;
      if (contains_token(buf, "ETA") != NULL)
	params->fatal |= ETA_A;
      if (contains_token(buf, "VEL2D") != NULL)
	params->fatal |= VEL2D;
      if (contains_token(buf, "VEL3D") != NULL)
	params->fatal |= VEL3D;
      if (contains_token(buf, "WIND") != NULL)
	params->fatal |= WIND;
      if (contains_token(buf, "T/S") != NULL)
	params->fatal |= TS;
      if (contains_token(buf, "NAN") != NULL)
	params->fatal |= NANF;
    }
  }

  /* Optional parameters */
  sprintf(keyword, "MOM_SCHEME");
  if (prm_read_char(fp, keyword, buf)) {
    params->momsc = 0;
    if (contains_token(buf, "ORDER1") != NULL)
      params->momsc |= ORDER1;
    if (contains_token(buf, "ORDER2") != NULL)
      params->momsc |= ORDER2;
    if (contains_token(buf, "VANLEER") != NULL)
      params->momsc |= VANLEER;
    if (contains_token(buf, "ANGULAR") != NULL)
      params->momsc |= ANGULAR;
    if (contains_token(buf, "ANGULAR3D") != NULL)
      params->momsc |= ANGULAR3D;
    if (contains_token(buf, "WIMPLICIT") != NULL)
      params->momsc |= WIMPLICIT;
    if (contains_token(buf, "EXPLICIT") != NULL)
      params->momsc |= EXPLICIT;
    if (contains_token(buf, "ADVECT_FORM") != NULL)
      params->momsc |= ADVECT_FORM;
    if (contains_token(buf, "WTOP_O4") != NULL)
      params->momsc |= WTOP_O4;
    if (contains_token(buf, "WTOP_O2") != NULL)
      params->momsc |= WTOP_O2;
    if (contains_token(buf, "ZERO_DRYK") != NULL)
      params->momsc |= ZERO_DRYK;
    if (contains_token(buf, "SHAPIRO") != NULL)
      params->momsc |= SHAPIRO;
    if (contains_token(buf, "LAGRANGE") != NULL)
      params->momsc |= LAGRANGE;
  }
  sprintf(keyword, "ORDER_SL");
  prm_read_int(fp, keyword, &params->osl);

  /* Use 4th order approximations for wtop by default */
  if (!(params->momsc & WTOP_O2))
    params->momsc |= WTOP_O4;
  if (params->momsc & SHAPIRO)
    params->filter = ADVECT;

  sprintf(keyword, "TRA_SCHEME");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "ORIGINAL") == 0)
      params->trasc = ORIGINAL;
    if (strcmp(buf, "ORDER1") == 0)
      params->trasc = ORDER1;
    if (strcmp(buf, "ORDER2") == 0)
      params->trasc = ORDER2;
    if (strcmp(buf, "QUICKEST_AD") == 0)
      params->trasc = QUICKEST_AD;
    if (strcmp(buf, "ORDER4") == 0)
      params->trasc = ORDER4;
    if (strcmp(buf, "QUICKEST") == 0)
      params->trasc = QUICKEST;
    if (strcmp(buf, "QUICKEST_CO") == 0)
      params->trasc = QUICKEST_CO;
    if (strcmp(buf, "VANLEER") == 0)
      params->trasc = VANLEER;
    if (strcmp(buf, "FFSL") == 0)
      params->trasc = FFSL;
    if (strcmp(buf, "ORDER2_UPWIND") == 0)
      params->trasc = ORDER2_UW;
    if (strcmp(buf, "LAGRANGE") == 0)
      params->trasc = LAGRANGE;
    if (strcmp(buf, "LAGRANGE|VANLEER") == 0) {
      params->trasc = LAGRANGE|VANLEER;
      params->trsplit = 1;
    }
  }
  sprintf(keyword, "ULTIMATE");
  if (prm_read_char(fp, keyword, buf))
    params->ultimate = is_true(buf);

  /* Setting robust > 1 will cap and smooth smagorinsky */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);

  sprintf(keyword, "STABILITY");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->stab = NONE;
    if (strcmp(buf, "ORIGINAL") == 0)
      params->stab = ORIGINAL;
    if (strcmp(buf, "SUB-STEP") == 0)
      params->stab = SUB_STEP;
    if (strcmp(buf, "SUB-STEP-NOSURF") == 0)
      params->stab = SUB_STEP_NOSURF;
    if (strcmp(buf, "SUB-STEP-TRACER") == 0)
      params->stab = SUB_STEP_TRACER;
  }
  sprintf(keyword, "MERGE_THIN");
  if (prm_read_char(fp, keyword, buf))
    params->thin_merge = is_true(buf);
  sprintf(keyword, "SIGMA");
  if (prm_read_char(fp, keyword, buf))
    params->sigma = is_true(buf);
  if (prm_read_char(fp, "SHOW_LAYERS", buf)) {
    params->show_layers = is_true(buf);
    if (params->show_layers)
      params->ntr++;
  }
  sprintf(keyword, "FILTERING");
  if (prm_read_char(fp, keyword, buf)) {
    params->filter = 0;
    if (strcmp(buf, "NONE") == 0)
      params->filter = NONE;
    if (strcmp(buf, "ADVECT") == 0)
      params->filter |= ADVECT;
  }
  read_trfilter(params, fp);

  sprintf(keyword, "PORUS_PLATE");
  if (prm_read_char(fp, keyword, params->reef_frac)) {
    params->porusplate = 1;
    params->ntr += 2;
  }

  sprintf(keyword, "TRANS_OUTPUT");
  if (prm_read_char(fp, keyword, buf))
    params->trout = is_true(buf);
  params->smooth = 0;
  sprintf(keyword, "SMOOTHING");
  prm_read_int(fp, keyword, &params->smooth);
  sprintf(keyword, "SMOOTH_VARS");
  prm_read_char(fp, keyword, params->smooth_v);
  sprintf(keyword, "SCALE_VARS");
  prm_read_char(fp, keyword, params->scale_v);
  sprintf(keyword, "NONLINEAR");
  if (prm_read_char(fp, keyword, buf))
    params->nonlinear = is_true(buf);
  if (prm_read_char(fp, "CALCDENS", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    params->calc_dens = is_true(fields[0]);
    if (n > 1) strcpy(params->densname, fields[1]);      
  }
  sprintf(keyword, "2D-MODE");
  if (prm_read_char(fp, keyword, buf))
    params->mode2d = is_true(buf);
  sprintf(keyword, "SLIP");
  prm_read_double(fp, keyword, &params->slipprm);

  read_means(params, fp, 0);

  /* Window sizes */
  read_window_info(params, fp);

  sprintf(keyword, "PARRAY_INTPL_BOTZ");
  if(prm_read_char(fp, keyword, buf))
  	params->parray_inter_botz = is_true(buf);

  sprintf(keyword, "GRID_EDGES");
  prm_read_int(fp, keyword, &params->edgef);

  /* Compatibility */
  read_compatible(params, fp);

  sprintf(keyword, "CFL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->cfl = NONE;
    if (strcmp(buf, "PASSIVE") == 0)
      params->cfl = PASSIVE;
    if (strcmp(buf, "PASSIVE|WVEL") == 0)
      params->cfl = PASSIVE | WVEL;
    if (strcmp(buf, "ACTIVE3D") == 0)
      params->cfl = ACTIVE3D;
    if (strcmp(buf, "ACTIVE3D|WVEL") == 0)
      params->cfl = ACTIVE3D | WVEL;
    if (strcmp(buf, "ACTIVE") == 0)
      params->cfl = ACTIVE;
    if (strcmp(buf, "ACTIVE|WVEL") == 0)
      params->cfl = ACTIVE | WVEL;
    if (params->cfl & (ACTIVE | ACTIVE3D))
      prm_set_errfn(hd_quit);
    sprintf(keyword, "CFL_DT");
    prm_read_char(fp, keyword, params->cfl_dt);
    prm_set_errfn(hd_silent_warn);
  }
  if (!(params->cfl & NONE))
    params->ntrS += 5;

  params->vorticity = 0;
  sprintf(keyword, "VORTICITY");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ABSOLUTE") != NULL) {
      params->vorticity |= ABSOLUTE;
      params->ntrS++;
    }
    if (contains_token(buf, "RELATIVE") != NULL) {
      params->vorticity |= RELATIVE;
      params->ntrS++;
    }
    if (contains_token(buf, "POTENTIAL") != NULL) {
      params->vorticity |= POTENTIAL;
      params->ntrS++;
    }
    if (contains_token(buf, "TENDENCY") != NULL) {
      params->vorticity |= TENDENCY;
      params->ntrS += 7;
      if (!(params->vorticity & RELATIVE)) {
        params->vorticity |= RELATIVE;
        params->ntrS++;
      }
    }
  } else
    params->vorticity = NONE;

  sprintf(keyword, "MIX_LAYER");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->mixlayer = NONE;
    if (strcmp(buf, "DENS_MIX") == 0)
      params->mixlayer = DENS_MIX;
    if (strcmp(buf, "TKE_MIX") == 0)
      params->mixlayer = TKE_MIX;
    if (strcmp(buf, "TEMP_MIX") == 0)
      params->mixlayer = TEMP_MIX;
    if (!(params->mixlayer & NONE))
      params->ntrS++;
  }
  sprintf(keyword, "CALC_FLUXES");
  if (prm_read_char(fp, keyword, params->trflux)) {
    if (strcmp(params->trflux, "NONE") != 0) {
      params->ntr += 4;
    }
  } else 
    sprintf(params->trflux, "NONE");
  sprintf(keyword, "CALC_PERCS");
  if (prm_read_char(fp, keyword, params->trperc)) {
    if (strcmp(params->trperc, "NONE") != 0) {
      params->ntr += 1;
      sprintf(keyword, "PERC_REGION");
      prm_read_char(fp, keyword, params->trpercr);
    }
  } else 
    sprintf(params->trperc, "NONE");
  sprintf(keyword, "MOM_TEND");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->tendf = is_true(buf))) {
      params->ntr += 13;
    }
  }
  sprintf(keyword, "TRA_TEND");
  if (prm_read_char(fp, keyword, params->trtend)) {
    if (strcmp(params->trtend, "NONE") != 0)
      params->ntr += 3;
    else
      sprintf(params->trtend, "%c", '\0');
  }
  sprintf(keyword, "FLUSHING_TR");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->trflsh = is_true(buf)))
      params->ntr += 1;
  }
  sprintf(keyword, "AGE_TR");
  if (prm_read_char(fp, keyword, params->trage)) {
    params->ntr += 1;
  }
  sprintf(keyword, "STERIC_HEIGHT");
  prm_read_double(fp, keyword, &params->lnm);
  if (params->lnm != 0.0)
    params->ntrS++;
  read_debug(params, fp);

  read_decorr(params, fp, 0);

  /* Totals diagnostics */
  read_totals(params, fp);

  /* Regions */
  read_region_info(params, fp);

  /* AVHRR SST */
  if (prm_read_char(fp, "AVHRR", params->avhrr_path)) {
    create_avhrr_list(params);
    params->avhrr = 1;
    params->ntrS++;
  }

  /* GHRSST SST */
  if (prm_read_char(fp, "GHRSST", params->ghrsst_path)) {
    create_ghrsst_list(params);
    params->ntrS++;
  }

  /* Particle tracking */
  if(prm_read_char(fp, "PT_InputFile", params->ptinname)) {
    params->do_pt = 1;
    params->ntr++;
  }

  /* Transport mode files */
  read_tmode_params(fp, params);

  /* Data assimilation */
  params->da = 0;
  sprintf(keyword, "DATA_ASSIMILATION");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->da = is_true(buf))) {
      /* Global dt and err for obs */
      double obs_dt  = -1;
      double obs_err = -1;
#ifndef HAVE_DA
      hd_quit("DATA ASSIMILATION has been invoked in the prm file but this version of SHOC is not built with DA support. Please either disable DA or reconfigure with the --enable-da option and rebuild");
#endif
      if (prm_read_char(fp, "DA_DT", buf))
	tm_scale_to_secs(buf, &params->da_dt);
      if (prm_read_char(fp, "DA_FCST_DT", buf))
	tm_scale_to_secs(buf, &params->da_fcst_dt);
      if (prm_read_char(fp, "DA_ANOMALY_FILE", buf)) {
	params->da_anom_file = (char *)malloc((strlen(buf)+1)*sizeof(char));
	strcpy(params->da_anom_file, buf);
      }
      if (prm_read_char(fp, "DA_ANOMALY_STATES", buf)) {
	params->da_anom_states = (char *)malloc((strlen(buf)+1)*sizeof(char));
	strcpy(params->da_anom_states, buf);
      }

      /* Always convert DT to seconds */
      prm_get_time_in_secs(fp, "DA_OBS_DT", &obs_dt);
      prm_read_double(fp, "DA_OBS_ERR", &obs_err);
      
      /* See if we're using localisation */
      params->da_ls = 0.0;
      prm_read_double(fp, "DA_LS", &params->da_ls);

      if (prm_read_int(fp, "DA_NOBS", &params->da_nobs)) {
	params->da_obs_names = (char **)malloc(params->da_nobs*sizeof(char *));
	params->da_obs_files = (char **)malloc(params->da_nobs*sizeof(char *));
	params->da_obs_types = (char **)malloc(params->da_nobs*sizeof(char *));
	params->da_obs_locs  = (char **)malloc(params->da_nobs*sizeof(char *));
	params->da_obs_dt    = (double *)malloc(params->da_nobs*sizeof(double));
	params->da_obs_errs  = (double *)malloc(params->da_nobs*sizeof(double));
	for (n = 0; n < params->da_nobs; n++) {
	  /* Name of this observation */
	  sprintf(keyword, "DA_OBS%1.1d.name", n);
	  params->da_obs_names[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_names[n] = strdup(buf);
	  /* File name for the data */
	  sprintf(keyword, "DA_OBS%1.1d.data", n);
	  params->da_obs_files[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_files[n] = strdup(buf);
	  /* What the platform type is */
	  sprintf(keyword, "DA_OBS%1.1d.type", n);
	  params->da_obs_types[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_types[n] = strdup(buf);
	  /* The location of this observation */
	  sprintf(keyword, "DA_OBS%1.1d.location", n);
	  params->da_obs_locs[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_locs[n] = strdup(buf);
	  /* Any DT  associated with this */
	  sprintf(keyword, "DA_OBS%1.1d.dt", n);
	  params->da_obs_dt[n] = obs_dt;
	  prm_read_double(fp, keyword, &params->da_obs_dt[n]);
	  /* Any errors associated with this */
	  sprintf(keyword, "DA_OBS%1.1d.error", n);
	  params->da_obs_errs[n] = obs_err;
	  prm_read_double(fp, keyword, &params->da_obs_errs[n]);
	}
      }
    }
  }

  /* Diagnistic numbers */
  params->ntr += numbers_init(params);

  /* Degree heating days */
  params->ntr += read_dhw(params, fp);

  /* River flow diagnostic */
  if (check_river_bdry(params, fp)) {
    params->riverflow = 1;
    params->ntrS++;
    if (check_bdry_options(params, fp) & OP_DYNAHC) {
      params->riverflow = 2;
      params->ntrS++;
      params->ntr++;
    }
  }

  params->u1_f = params->u1av_f = params->u2_f = params->u2av_f = 0;
  sprintf(keyword, "U1_OMIT");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ADVECT") != NULL)
      params->u1_f |= ADVECT;
    if (contains_token(buf, "HDIFF") != NULL)
      params->u1_f |= HDIFF;
    if (contains_token(buf, "VDIFF") != NULL)
      params->u1_f |= VDIFF;
    if (contains_token(buf, "PRESS_BT") != NULL)
      params->u1_f |= PRESS_BT;
    if (contains_token(buf, "PRESS_BC") != NULL)
      params->u1_f |= PRESS_BC;
    if (contains_token(buf, "CORIOLIS") != NULL)
      params->u1_f |= CORIOLIS;
  } else
    params->u1_f = NONE;
  sprintf(keyword, "U2_OMIT");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ADVECT") != NULL)
      params->u2_f |= ADVECT;
    if (contains_token(buf, "HDIFF") != NULL)
      params->u2_f |= HDIFF;
    if (contains_token(buf, "VDIFF") != NULL)
      params->u2_f |= VDIFF;
    if (contains_token(buf, "PRESS_BT") != NULL)
      params->u2_f |= PRESS_BT;
    if (contains_token(buf, "PRESS_BC") != NULL)
      params->u2_f |= PRESS_BC;
    if (contains_token(buf, "CORIOLIS") != NULL)
      params->u2_f |= CORIOLIS;
  } else
    params->u2_f = NONE;
  sprintf(keyword, "U1AV_OMIT");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ADVECT") != NULL)
      params->u1av_f |= ADVECT;
    if (contains_token(buf, "HDIFF") != NULL)
      params->u1av_f |= HDIFF;
    if (contains_token(buf, "VDIFF") != NULL)
      params->u1av_f |= VDIFF;
    if (contains_token(buf, "PRESS_BT") != NULL)
      params->u1av_f |= PRESS_BT;
    if (contains_token(buf, "PRESS_BC") != NULL)
      params->u1av_f |= PRESS_BC;
    if (contains_token(buf, "CORIOLIS") != NULL)
      params->u1av_f |= CORIOLIS;
  } else
    params->u1av_f = NONE;
  sprintf(keyword, "U2AV_OMIT");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "ADVECT") != NULL)
      params->u2av_f |= ADVECT;
    if (contains_token(buf, "HDIFF") != NULL)
      params->u2av_f |= HDIFF;
    if (contains_token(buf, "VDIFF") != NULL)
      params->u2av_f |= VDIFF;
    if (contains_token(buf, "PRESS_BT") != NULL)
      params->u2av_f |= PRESS_BT;
    if (contains_token(buf, "PRESS_BC") != NULL)
      params->u2av_f |= PRESS_BC;
    if (contains_token(buf, "CORIOLIS") != NULL)
      params->u2av_f |= CORIOLIS;
  } else
    params->u2av_f = NONE;

  sprintf(keyword, "ALERT");
  if (prm_read_char(fp, keyword, params->alert)) {
    if (contains_token(params->alert, "NONE") != NULL)
      memset(params->alert, 0, sizeof(params->alert));
    if (contains_token(params->alert, "ACTIVE") != NULL)
      params->ntrS += 4;
    sprintf(keyword, "ALERT_CODE");
    if (prm_read_char(fp, keyword, params->alertc)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      n = parseline(params->alertc, fields, MAXNUMARGS);
      if (n == 11) {
	params->eta_f = atoi(fields[0]);
	params->vel2d_f = atoi(fields[1]);
	params->vel3d_f = atoi(fields[2]);
	params->wvel_f = atoi(fields[3]);
	params->tend_f = atoi(fields[4]);
	params->div2d_f = atoi(fields[5]);
	params->div3d_f = atoi(fields[6]);
	params->cfl_f = atoi(fields[7]);
	params->ts_f = atoi(fields[8]);
	params->shear_f = atoi(fields[9]);
	params->hdiff_f = atoi(fields[10]);
      }
    }
  }
  sprintf(keyword, "ALERT_DT");
  prm_read_char(fp, keyword, params->alert_dt);

  params->rampf = 0;
  sprintf(keyword, "RAMPVARS");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "WIND") != NULL)
      params->rampf |= WIND;
    if (contains_token(buf, "ETA_RELAX") != NULL)
      params->rampf |= ETA_RELAX;
    if (contains_token(buf, "FILEIN") != NULL)
      params->rampf |= FILEIN;
    if (contains_token(buf, "CUSTOM") != NULL)
      params->rampf |= CUSTOM;
    if (contains_token(buf, "TIDALH") != NULL)
      params->rampf |= TIDALH;
    if (contains_token(buf, "TIDALC") != NULL)
      params->rampf |= TIDALC;
    if (contains_token(buf, "TIDEBC") != NULL)
      params->rampf |= TIDEBC;
    if (contains_token(buf, "FLATHR") != NULL)
      params->rampf |= FLATHR;
    if (contains_token(buf, "INV_BARO") != NULL)
      params->rampf |= INV_BARO;
    if (contains_token(buf, "FLUX_ADJUST") != NULL)
      params->rampf |= FLUX_ADJUST;
    if (contains_token(buf, "STOKES") != NULL)
      params->rampf |= STOKES;
  } else
    params->rampf = WIND|FILEIN|CUSTOM|TIDALH|TIDALC|TIDEBC|FLATHR|INV_BARO|ETA_RELAX|STOKES;

  sprintf(keyword, "VELMAX");
  prm_read_double(fp, keyword, &params->velmax);
  sprintf(keyword, "VELMAX_2D");
  if(!params->velmax2d && !prm_read_double(fp, keyword, &params->velmax2d))
    params->velmax2d = params->velmax;
  sprintf(keyword, "ETA_DIFF");
  prm_read_double(fp, keyword, &params->etadiff);
  sprintf(keyword, "WMAX");
  prm_read_double(fp, keyword, &params->wmax);

  /* Mandatory parameters */
  prm_set_errfn(quit);
  sprintf(keyword, "HMIN");
  prm_read_double(fp, keyword, &params->hmin);
  sprintf(keyword, "UF");
  prm_read_double(fp, keyword, &params->uf);
  sprintf(keyword, "QBFC");
  prm_read_double(fp, keyword, &params->quad_bfc);
  sprintf(keyword, "RAMPSTART");
  prm_get_time_in_secs(fp, keyword, &params->rampstart);
  sprintf(keyword, "RAMPEND");
  prm_get_time_in_secs(fp, keyword, &params->rampend);
  if (params->rampend < params->rampstart)
    hd_quit("readparam: ramp end time before ramp start time\n");
  if(params->rampstart > params->t)
    hd_warn("readparam: ramp start time after model start time\n");

  /* Time */
  sprintf(keyword, "DT");
  prm_get_time_in_secs(fp, keyword, &params->grid_dt);
  sprintf(keyword, "IRATIO");
  prm_read_int(fp, keyword, &params->iratio);
  if ((params->iratio + 1) % 2 == 0) {
    params->iratio++;
    hd_warn("Only even iratio values allowed : setting iratio=%d\n",
            params->iratio);
  }
  sprintf(keyword, "TRATIO");
  prm_read_double(fp, keyword, &params->tratio);
  sprintf(keyword, "OUTPUT_TIMEUNIT");
  if (prm_read_char(fp, keyword, params->output_tunit)) {
    char *p = strstr(params->output_tunit, "since");
    if (p == NULL)
      sprintf(params->output_tunit, "%s %s",
              strdup(params->output_tunit),
              strstr(params->timeunit, "since"));
  } else
    strcpy(params->output_tunit, params->timeunit);

  /* Constants */
  prm_set_errfn(hd_quit);
  strcpy(params->prmname, prmname);
  sprintf(keyword, "CODEHEADER");
  prm_read_char(fp, keyword, params->codeheader);
  strcpy(codeheader, params->codeheader);
  sprintf(keyword, "PARAMETERHEADER");
  prm_read_char(fp, keyword, params->parameterheader);
  strcpy(parameterheader, params->parameterheader);
  sprintf(keyword, "G");
  prm_read_double(fp, keyword, &params->g);
  sprintf(keyword, "SPECHEAT");
  prm_read_double(fp, keyword, &params->spec_heat);
  sprintf(keyword, "AIRDENS");
  prm_read_double(fp, keyword, &params->air_dens);
  air_dens = params->air_dens;
  sprintf(keyword, "AMBIENT_AIR_PRESSURE");
  prm_read_double(fp, keyword, &params->ambpress);
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "PROJECTION");
  if (prm_read_char(fp, keyword, params->projection)) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }
  get_output_path(params, fp);
  prm_read_char(fp, "restart_name", params->restart_name);
  if (prm_read_char(fp, "restart_dt", buf)) {
     tm_scale_to_secs(buf, &params->restart_dt);
  }
  if (!prm_read_char(fp, "OutputTransport", params->trkey)) {
    if (params->trout) {
      strcpy(buf, params->prmname);
      stripend(buf);
      strcpy(params->trkey, buf);
    }
  }

  /* Horizontal diffusion */
  read_hdiff(params, fp, 0);
  /*
  if (params->smagorinsky > 0.0)
    params->ntr++;
  */

  z0_init(params);
  /*
  prm_read_double(fp, keyword, &params->z0);
  if (params->z0 > 1.0)
    hd_quit_and_dump
      ("params_read() : Bottom roughness parameter too large, %f\n",
       params->z0);
  */

  /* Mixing scheme */
  sprintf(params->mixsc, "%c", '\0');
  sprintf(params->s_func, "%c", '\0');
  params->min_tke = 0.0;
  params->min_diss = 0.0;
  params->vz0 = 0.0;
  params->kz0 = 0.0;
  prm_set_errfn(hd_quit);
  sprintf(keyword, "MIXING_SCHEME");
  prm_read_char(fp, keyword, params->mixsc);
  if (strcmp(params->mixsc, "constant") == 0 &&
      params->mixlayer == TKE_MIX) {
    params->mixlayer = NONE;
    hd_warn("MIX_LAYER = TKE_MIX not functional with MIXING_SCHEME = constant.\n");
  }
  sprintf(keyword, "KZ0");
  prm_read_double(fp, keyword, &params->kz0);
  sprintf(keyword, "VZ0");
  prm_read_double(fp, keyword, &params->vz0);
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "KZ_ALPHA");
  prm_read_double(fp, keyword, &params->kz_alpha);
  sprintf(keyword, "VZ_ALPHA");
  prm_read_double(fp, keyword, &params->vz_alpha);
  sprintf(keyword, "ZS");
  if (!(prm_read_double(fp, keyword, &params->zs))) {
    if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
        strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
        strcmp(params->mixsc, "harcourt") == 0 ||
        strcmp(params->mixsc, "k-e") == 0 ||
        strcmp(params->mixsc, "k-w") == 0 ||
        strcmp(params->mixsc, "W88") == 0 ||
        strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
      hd_quit("%s requires surface length scale, ZS.\n", params->mixsc);
  }
  sprintf(keyword, "STABILITY_FUNC");
  prm_read_char(fp, keyword, params->s_func);
  if (prm_read_char(fp, "SMOOTH_VzKz", buf))
    params->smooth_VzKz = is_true(buf);
  sprintf(keyword, "LMIN");
  prm_read_double(fp, keyword, &params->Lmin);
  sprintf(keyword, "E");
  prm_read_double(fp, keyword, &params->eparam);
  sprintf(keyword, "WAVE_ALPHA");
  prm_read_double(fp, keyword, &params->wave_alpha);
  sprintf(keyword, "WAVE_B1");
  prm_read_double(fp, keyword, &params->wave_b1);
  sprintf(keyword, "WAVE_HEIGHT_FACT");
  prm_read_double(fp, keyword, &params->wave_hf);
  sprintf(keyword, "MIN_TKE");
  prm_read_double(fp, keyword, &params->min_tke);
  sprintf(keyword, "MIN_DISS");
  prm_read_double(fp, keyword, &params->min_diss);
  if (strcmp(params->mixsc, "k-e") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "k-w") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "W88") == 0)
    params->ntr += 2;
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0)
    params->ntr += 4;
  if (strcmp(params->mixsc, "harcourt") == 0)
    params->ntr += 4;

  /* Coriolis */
  prm_set_errfn(hd_silent_warn);
  params->coriolis = NULL;
  sprintf(keyword, "CORIOLIS");
  prm_read_darray(fp, keyword, &params->coriolis, &params->nvals);

  /* Surface height */
  sprintf(params->eta_init, "%c", '\0');
  sprintf(keyword, "SURFACE");
  prm_read_char(fp, keyword, params->eta_init);
  sprintf(params->vel_init, "%c", '\0');
  sprintf(keyword, "VELOCITY");
  prm_read_char(fp, keyword, params->vel_init);

  /* Grid refinement */
  if (prm_read_char(fp, "GRID_REFINEMENT", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    if((params->dozoom = is_true(fields[0]))) {
      read_zoom(params, fp);
      if (n > 1 && strcmp(fields[1], "PRECONDITIONED") == 0)
	params->dozoom |= PRECOND;
    }
  }

  /* Blending */
  if (prm_read_int(fp, "NE1_BLEND", &params->nbl1)) {
    params->e1_blend = (char **)malloc(params->nbl1 * sizeof(char *));
    for (n = 0; n < params->nbl1; n++) {
      params->e1_blend[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      sprintf(keyword, "E1_BLEND%1d", n);
      prm_read_char(fp, keyword, params->e1_blend[n]);
    }
  }
  if (prm_read_int(fp, "NE2_BLEND", &params->nbl2)) {
    params->e2_blend = (char **)malloc(params->nbl2 * sizeof(char *));
    for (n = 0; n < params->nbl2; n++) {
      params->e2_blend[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      sprintf(keyword, "E2_BLEND%1d", n);
      prm_read_char(fp, keyword, params->e2_blend[n]);
    }
  }

  /* Cells explicitly defined as OUTSIDE */
  read_blocks(fp, "NOUTSIDE", &params->noutside, &params->oute1, &params->oute2);
  read_blocks(fp, "NLAND", &params->nland, &params->lande1, &params->lande2);

  /* Surface relaxation */
  read_eta_relax(params, fp);
  read_vel_relax(params, fp);

  /* Wind */
  sprintf(params->wind, "%c", '\0');
  params->wind_dt = params->storm_dt = 0.0;
  prm_set_errfn(hd_silent_warn);
  params->wind_type = SPEED;
  if (prm_read_char(fp, "WIND_TS", params->wind)) {
    prm_set_errfn(hd_quit);
    prm_get_time_in_secs(fp, "WIND_INPUT_DT", &params->wind_dt);
    prm_read_double(fp, "WIND_SPEED_SCALE", &params->wind_scale);
    prm_set_errfn(hd_silent_warn);
    if (prm_read_char(fp, "WIND_STRESS_FCTN", buf)) {
      prm_set_errfn(hd_quit);
      if(strcmp(buf, "L&P") == 0) {
	params->stress_fn = LARGEPOND;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else if(strcmp(buf, "B") == 0) {
	params->stress_fn = BUNKER;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else if(strcmp(buf, "K/W") == 0) {
	params->stress_fn = KITIAG;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
      } else if(strcmp(buf, "Ko") == 0) {
	params->stress_fn = KONDO;
	prm_read_double(fp, "WIND_STRESS_REFH", &params->dlv0);
	if (prm_read_char(fp, "WIND_STRESS_NEUTRAL", buf) > 0)
	  params->neutral = is_true(buf);
      } else {
	prm_set_errfn(hd_quit);
	params->stress_fn = ORIGINAL;
	prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
	prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
	prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
	prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
      }
    } else {
      prm_set_errfn(hd_quit);
      params->stress_fn = ORIGINAL;
      prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
      prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
      prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
      prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
    }
    if (prm_read_char(fp, "WIND_TYPE", buf)) {
      if(strcmp(buf, "STRESS") == 0)
	params->wind_type = STRESS;
      else
	params->wind_type = SPEED;
    }
  }
  if (prm_read_int(fp, "NSTORM", &params->nstorm)) {
    prm_set_errfn(hd_quit);
    prm_get_time_in_secs(fp, "STORM_INPUT_DT", &params->storm_dt);
    if (!params->wind_scale || !params->dlv0) {
      prm_read_double(fp, "WIND_SPEED_SCALE", &params->wind_scale);
      prm_read_double(fp, "DRAG_LAW_V0", &params->dlv0);
      prm_read_double(fp, "DRAG_LAW_V1", &params->dlv1);
      prm_read_double(fp, "DRAG_LAW_CD0", &params->dlc0);
      prm_read_double(fp, "DRAG_LAW_CD1", &params->dlc1);
    }
  }

  /* Atmospherics */
  prm_set_errfn(hd_silent_warn);
  sprintf(params->patm, "%c", '\0');
  sprintf(keyword, "PRESSURE");
  prm_read_char(fp, keyword, params->patm);
  sprintf(keyword, "PRESSURE_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->patm_dt);
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  prm_read_char(fp, keyword, params->precip);
  sprintf(keyword, "PRECIPITATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  prm_read_char(fp, keyword, params->evap);
  sprintf(keyword, "EVAPORATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  sprintf(params->airtemp, "%c", '\0');
  sprintf(keyword, "AIRTEMP");
  prm_read_char(fp, keyword, params->airtemp);
  sprintf(keyword, "AIRTEMP_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->airtemp_dt);
  sprintf(params->rh, "%c", '\0');
  sprintf(keyword, "HUMIDITY");
  prm_read_char(fp, keyword, params->rh);
  sprintf(keyword, "HUMIDITY_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->rh_dt);
  sprintf(params->cloud, "%c", '\0');
  sprintf(keyword, "CLOUD");
  prm_read_char(fp, keyword, params->cloud);
  sprintf(keyword, "CLOUD_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->cloud_dt);
  sprintf(params->swr, "%c", '\0');
  sprintf(keyword, "RADIATION");
  prm_read_char(fp, keyword, params->swr);
  sprintf(keyword, "RADIATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->swr_dt);
  prm_read_double(fp, "ALBEDO", &params->albedo);
  sprintf(params->light, "%c", '\0');
  sprintf(keyword, "LIGHT");
  prm_read_char(fp, keyword, params->light);
  sprintf(keyword, "LIGHT_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->light_dt);
  prm_read_double(fp, "ALBEDO_LIGHT", &params->albedo_l);
  sprintf(params->wetb, "%c", '\0');
  sprintf(keyword, "WET_BULB");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "WET_BULB_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);
  sprintf(keyword, "DEW_POINT");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "DEW_POINT_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);

  /* Heatflux */
  prm_set_errfn(hd_silent_warn);
  params->heatflux = NONE;
  sprintf(params->hftemp, "%c", '\0');
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->heatflux = NONE;
    if (strcmp(buf, "ORIGINAL") == 0) {
      hd_warn
        ("ORIGINAL HEATFLUX formulation no longer supported : using BULK formulation\n");
      params->heatflux = ADVANCED;
    }
    if ((strcmp(buf, "NET_HEAT") == 0) || 
	(strcmp(buf, "COMP_HEAT_MOM") == 0) || 
	(strcmp(buf, "COMP_HEAT") == 0)) {
    if (strcmp(buf, "COMP_HEAT") == 0)
      params->heatflux = COMP_HEAT;
    if (strcmp(buf, "COMP_HEAT_MOM") == 0)
      params->heatflux = COMP_HEAT_MOM;
    if (strcmp(buf, "NET_HEAT") == 0)
      params->heatflux = NET_HEAT;

      /* Read the swr parameters                                     */
    read_swr(params, fp, 0);
    }

    if (strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0 || 
	params->heatflux & ADVANCED) {
      params->heatflux = ADVANCED;

      /* Read the swr parameters                                     */
      read_swr(params, fp, 0);
      read_hf_bulk(params, fp);
      params->ntrS += 5;
    }
    if (strcmp(buf, "SURF_RELAX") == 0) {
      params->heatflux = SURF_RELAX;
      sprintf(keyword, "HEATFLUX_TEMP");
      if (!(prm_read_char(fp, keyword, params->hftemp)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TEMP file.\n");
      sprintf(keyword, "HEATFLUX_TEMP_DT");
      if (!(prm_get_time_in_secs(fp, keyword, &params->hftemp_dt)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TEMP_DT constant.\n");
      sprintf(keyword, "HEATFLUX_TC");
      if (!(prm_get_time_in_secs(fp, keyword, &params->hftc)))
	hd_warn
	  ("params_read() : SURF_RELAX heatflux requires HEATFLUX_TC constant.\n");
    }
    if (strcmp(buf, "AVHRR") == 0) {
      params->heatflux = AVHRR;
      if (!strlen(params->avhrr_path)) {
	hd_warn
	  ("params_read() : AVHRR heatflux requires an AVHRR diagnostic path name.\n");
	params->heatflux = NONE;
      } else {
	sprintf(keyword, "HEATFLUX_TC");
	if (!(prm_get_time_in_secs(fp, keyword, &params->hftc))) {
	  hd_warn
	    ("params_read() : AVHRR heatflux requires HEATFLUX_TC constant.\n");
	  params->hftc = 86400.0;
	}
      }
    }

    if (strcmp(buf, "INVERSE") == 0) {
      params->heatflux = INVERSE;
      sprintf(keyword, "HEATFLUX_TEMP");
      prm_read_char(fp, keyword, params->hftemp);
      sprintf(keyword, "HEATFLUX_TEMP_DT");
      prm_get_time_in_secs(fp, keyword, &params->hftemp_dt);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      if(params->mixlayer & NONE) {
	hd_quit_and_dump
	  ("params_read() : INVERSE heatflux requires MIX_LAYER to be set.\n");
      }
      params->ntrS += 1;
    }
    if (strcmp(buf, "NET_HEAT") == 0) {
      params->heatflux = NET_HEAT;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 2;
    }
    if (strcmp(buf, "COMP_HEAT") == 0) {
      params->heatflux = COMP_HEAT;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 5;
    }
    if (strcmp(buf, "COMP_HEAT_MOM") == 0) {
      params->heatflux = COMP_HEAT_MOM;
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 5;
    }
    sprintf(keyword, "HEATFLUX_RAMP");
    if (prm_read_char(fp, keyword, buf))
      tm_scale_to_secs(buf, &params->hf_ramp);
    else
      params->hf_ramp = params->t;
  }

  /* Saltflux */
  prm_set_errfn(hd_silent_warn);
  params->saltflux = NONE;
  sprintf(keyword, "SALTFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->saltflux = NONE;
    if (strcmp(buf, "ORIGINAL") == 0) {
      params->saltflux = ORIGINAL;
      params->ntrS += 3;
    }
    if (strcmp(buf, "ADVANCED") == 0) {
      params->saltflux = ADVANCED;
      params->ntrS += 3;
    }
    if (strcmp(buf, "BULK") == 0) {
      params->saltflux = BULK;
      params->ntrS += 1;
    }
  }

  prm_set_errfn(hd_silent_warn);

  /* Sediment geometry */
  params->sednz = 0;
  /* Note : sediment layers are read in assuming the first layer is */
  /* closest to the water column. This is reversed when copying to */
  /* the master so that the sediment origin is the deepest layer.  */
  if (prm_read_darray(fp, "LAYERFACES_SED", &params->gridz_sed,
                      &params->sednz)) {
    if (--params->sednz < 1)
      hd_quit("Number of sediment layers must be 2 or more\n");
  } else {
    double *dz_sed = NULL;
    if (prm_read_int(fp, "NSEDLAYERS", &params->sednz)) {
      dz_sed = d_alloc_1d(params->sednz);
      if (prm_read_darray(fp, "NSEDLAYERS", &dz_sed, &params->sednz)) {
        if (params->sednz) {
          params->gridz_sed = d_alloc_1d(params->sednz + 1);

          params->gridz_sed[params->sednz] = 0.0;
          for (m = params->sednz - 1; m >= 0; m--) {
            params->gridz_sed[m] = params->gridz_sed[m + 1] -
              dz_sed[params->sednz - m - 1];
         }
        }
      }
      if(dz_sed)
        d_free_1d(dz_sed);
    }
  }

#if defined(HAVE_SEDIMENT_MODULE)
  read_sediments(params, fp, &params->ntr);
#endif

#if defined(HAVE_WAVE_MODULE)
  params->do_wave = 0;
  sprintf(keyword, "DO_WAVES");
  if (prm_read_char(fp, keyword, buf) && is_true(buf)) {
    params->do_wave = is_true(buf);
    prm_get_time_in_secs(fp, "WAVES_DT", &params->wavedt);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  read_ecology(params, fp, &params->ntr);
#endif
  /* Waves */
  read_waves(params, fp, 0);
  if (params->waves & STOKES_DRIFT && params->tendf) params->ntr+=2;

  /* Library error handling */
  sprintf(keyword, "LIB_ERROR_FCN");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "QUIT") == 0)
      params->gint_errfcn = LFATAL;
    if (strcmp(buf, "WARN") == 0)
      params->gint_errfcn = LWARN;
  }

  /* 3D Tracer constants and variables.  */
  /* Tracers relating to certain diagnostics etc. are automatically */
  /* generated if the relevant flag is set. If these tracers are */
  /* manually defined in the input parameter file, then the manual */
  /* definition over-rides the auto generation, and the defined */
  /* tracer_info is used.  */
  /* Always assume salt and temp are auto tracers */
  params->ntr += 2;
  params->atr = params->ntr;
  tracer_setup(params, fp);
  create_tracer_3d(params);
  for (n = 0; n < params->atr; n++) {
    params->trinfo_3d[n].m = -1;
    params->trinfo_3d[n].n = n;
  }

  /* Get the tracers which undergo diffusion (horizontal and */
  /* vertical)  and store in tdif_h[] and tdif_v[]. */
  params->ntdif_h = params->ntdif_v = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->ntdif_h++;

      if (params->compatible & V5342) {
#if defined(HAVE_SEDIMENT_MODULE)
        if (params->do_sed) {
          if (strcmp(tracer->name, "salt") == 0 ||
              strcmp(tracer->name, "temp") == 0)
            params->ntdif_v++;
        } else
          params->ntdif_v++;
#else
        params->ntdif_v++;
#endif
      } else {
#if defined(HAVE_SEDIMENT_MODULE)
        if (params->do_sed) {
	  if ( (tracer->type & WATER && !(tracer->type & CLOSURE) && !(tracer->type & SEDIMENT)) ||
	       strcmp(tracer->name, "salt") == 0 || strcmp(tracer->name, "temp") == 0 )
	    params->ntdif_v++;
	}  else
	  if(!(tracer->type & CLOSURE))
	    params->ntdif_v++;
#else
	if(!(tracer->type & CLOSURE)) {
	  params->ntdif_v++;
	}
#endif
      }
    }
  }
  if (params->ntdif_h)
    params->tdif_h = i_alloc_1d(params->ntdif_h);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->tdif_h[m] = n;
      m++;
    }
  }
  if (params->ntdif_v)
    params->tdif_v = i_alloc_1d(params->ntdif_v);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      if (params->compatible & V5342) {
#if defined(HAVE_SEDIMENT_MODULE)
	if (params->do_sed) {
	  if (strcmp(tracer->name, "salt") == 0 ||
	      strcmp(tracer->name, "temp") == 0) {
	    params->tdif_v[m] = n;
	    m++;
	  }
	} else {
	  params->tdif_v[m] = n;
	  m++;
	}
#else
	params->tdif_v[m] = n;
	m++;
#endif
      } else {
#if defined(HAVE_SEDIMENT_MODULE)
        if (params->do_sed) {
	  if ( (tracer->type & WATER && !(tracer->type & CLOSURE) && !(tracer->type & SEDIMENT)) ||
	       strcmp(tracer->name, "salt") == 0 || strcmp(tracer->name, "temp") == 0 ) {
	    params->tdif_v[m] = n;
            m++;
	  }
        } else {
	  if(!(tracer->type & CLOSURE)) {
	    params->tdif_v[m] = n;
            m++;
	  }
	}
#else
        if(!(tracer->type & CLOSURE)) {
          params->tdif_v[m] = n;
          m++;
	}
#endif
      }
    }
  }

  /* Sediment tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  tracer_read(fp, NULL, SEDIM, hd_quit, hd_warn, hd_silent_warn,
        &params->nsed, &params->trinfo_sed);
  if (params->nsed && !params->sednz)
    hd_quit("Sediment tracers requires NSEDLAYERS > 0\n");

  /* 2D Tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  params->atrS = params->ntrS;
  tracer_read(fp, NULL, INTER, hd_quit, hd_warn, hd_silent_warn,
              &params->ntrS, &params->trinfo_2d);

  for (n = 0; n < params->ntrS; n++)
    if (params->trinfo_2d[n].increment & ETA_M) {
      if (params->etarlx & NONE)
	params->etarlx = INCREMENT;
      else
	params->etarlx |= INCREMENT;
    }

  /* Timeseries file caching */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "CACHE_TSFILES");
  if (prm_read_char(fp, keyword, buf) > 0)
    params->tsfile_caching = is_true(buf);
  else
    params->tsfile_caching = 1;
  emstag(LDEBUG,"hd:readparam:params_read","Setting ts_file_cachinf to: %s",(params->tsfile_caching?"true":"false"));

  /* Explicit mappings */
  read_explicit_maps(params, fp);

  /* Process exclusion */
  read_exclude_points(params, fp);

  /* Open boudaries */
  get_bdry_params(params, fp);

  /* Allocate memory for the open boundary data structure */
  params->open =
    (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);

  for (n = 0; n < params->nobc; n++) {

    /* Allocate memory for the boundary structures */
    params->open[n] = OBC_alloc();

    /* Read the open boudary conditions */
    params->open[n]->ntr = params->ntr;
    params->open[n]->atr = params->atr;
    get_OBC_conds(params, params->open[n], fp, n, params->trinfo_3d);
    get_OBC_relax(params, params->open[n], fp, n);
    get_obc_list(params->open[n], fp, n, "BOUNDARY");
    prm_set_errfn(hd_quit);
  }

  /* Read the csr tide model paths if required */
  for (n = 0; n < params->nobc; n++) {
    if(params->open[n]->bcond_ele & TIDALH) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TIDE_CSR_ORTHOWEIGHTS");
      prm_read_char(fp, keyword, params->orthoweights);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      prm_set_errfn(hd_silent_warn);
    }
    if(params->open[n]->bcond_ele & TIDALC ||
       params->open[n]->bcond_nor & TIDALC ||
       params->open[n]->bcond_nor2d & TIDALC ||
       params->open[n]->bcond_tan & TIDALC ||
       params->open[n]->bcond_tan2d & TIDALC) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      sprintf(keyword, "TIDE_CONSTITUENTS");
      prm_read_char(fp, keyword, params->tide_con_file);
    }
  }

  sprintf(keyword, "TIDAL_REMOVAL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "CSR") == 0) params->tide_r = CSR_R;
    if (strcmp(buf, "MEAN") == 0) params->tide_r = MEAN_R;
  }

  /* Whether we read do threaded I/O on some forcing */
  params->thIO = 0;
  if (prm_read_char(fp, "SCHED_MODE", buf)) {
    if (strcasecmp(buf, "pthreads") == 0) {
#ifdef HDF5_THREADSAFE
      params->thIO = 1;
#else
      hd_quit("SCHED_MODE error: This executable is not built with the correct HDF5 threadsafe library version via NETCDF4");
#endif
    }
  }

  /* Output files to write to setup.txt */
  sprintf(keyword, "OutputFiles");
  if (prm_read_char(fp, keyword, buf)) {
    if (fopen(buf, "r") == NULL) {
      params->ndf = atoi(buf);
      params->d_name = (char **)malloc(params->ndf * sizeof(char *));
      params->d_tinc = (char **)malloc(params->ndf * sizeof(char *));
      params->d_sync = (char **)malloc(params->ndf * sizeof(char *));
      params->d_filetype = (char **)malloc(params->ndf * sizeof(char *));
      for (n = 0; n < params->ndf; n++) {
	params->d_name[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.name",n);
	prm_read_char(fp, keyword, params->d_name[n]);
	if(!endswith(params->d_name[n],".nc"))
	  strcat(params->d_name[n],".nc");
	params->d_filetype[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.filetype",n);
	prm_read_char(fp, keyword, params->d_filetype[n]);
	prm_read_char(fp, keyword, buf);
	if (strcmp(params->d_filetype[n], "mom") == 0) params->domom |= RDGRID;
	if (strcmp(params->d_filetype[n], "roms") == 0) params->doroms |= RDGRID;
	
	params->d_tinc[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.tinc",n);
	prm_read_char(fp, keyword, params->d_tinc[n]);
	
	params->d_sync[n] = (char *)malloc(sizeof(char)*MAXLINELEN);
	strcpy(params->d_sync[n], params->d_tinc[n]);
	sprintf(keyword, "file%1.1d.sync_dt",n);
	prm_read_char(fp, keyword, params->d_sync[n]);
      }
    }
  }

  sprintf(keyword, "REREAD_INPUT");/*UR*/
  if(prm_read_char(fp, keyword, buf1))
    {/*
       if(strcmp(buf1,"YES") == 0 || strcmp(buf1,"1") == 0)
       params->sliding = 1;*/
    }

  if (DEBUG("init_m"))
    dlog("init_m", "Input parameters read OK\n");
  return (params);
}

/* END params_read()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set input parameters to standard values                */
/*-------------------------------------------------------------------*/
parameters_t *auto_params(FILE * fp, int autof)
{
  double **bathy;               /* Bathymetry array */
  int i, j, m, n, c, nz;        /* Counters */
  int ii, jj, nr, rs, tn;       /* Counters */
  int bdry = 0;                 /* Tracks open boundaries */
  open_bdrys_t *open;           /* Pointer to boundary structure */
  char keyword[MAXSTRLEN];      /* Input dummy */
  char buf[MAXSTRLEN];          /* Input dummy */
  size_t sz;                    /* Size of tracer structure */
  double d1;                    /* Dummy */
  int dataf;                    /* Flag for boundary data */
  int is, ie, js, je;           /* Limits of grid */
  int *ib, *jb, *rb, *rf;
  int **mask, *nib;
  int ic = 0, jc = 0;
  double *rhc;
  FILE* fptr = fp;

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the parameter data structure */
  params = params_alloc();

  /* Read in the name of the input file */
  params->prmfd = fp;
  switch (autof) {
  case 1: /* -a, auto mode */
    params->runmode = AUTO;
    break;

  case 2: /* -g, input file generation */
    params->runmode = DUMP | EXIT;
    break;

  case 3: /* -ag, auto mode with initialisation termination */
    params->runmode = AUTO | EXIT;
    break;

  case 4: /* -r, ROAM mode */
    params->runmode = AUTO | ROAM;
    break;

  case 5: /* -rg, ROAM mode with initialisation termination */
    params->runmode = AUTO | ROAM | EXIT;
    break;

  case 6: /* -rs, ROAM mode with initialisation from previous run */
    params->runmode = AUTO | ROAM | RE_ROAM;
    break;

  case 7: /* -rso, ROAM mode with initialisation from previous run, T/S from OFAM */
    params->runmode = AUTO | ROAM | RE_ROAM | RE_NOTS;
    break;

  case 10:/* -rt, ROAM PRE_MARVL mode - same as -rg plus generates transport file */
    params->runmode = AUTO | ROAM | PRE_MARVL | EXIT;
    break;

  }

  /*-----------------------------------------------------------------*/
  /* Read mandatory parameters */
  set_default_param(params);
  prm_set_errfn(hd_quit);
  sprintf(keyword, "START_TIME");
  prm_read_char(fp, keyword, params->start_time);
  sprintf(keyword, "STOP_TIME");
  prm_read_char(fp, keyword, params->stop_time);
  prm_get_time_in_secs(fp, "START_TIME", &params->t);
  sprintf(keyword, "NCE1");
  prm_read_int(fp, keyword, &nce1);
  sprintf(keyword, "NCE2");
  prm_read_int(fp, keyword, &nce2);
  params->nce1 = nce1;
  params->nce2 = nce2;
  prm_set_errfn(hd_silent_warn);

  /* edgef>0 sets a land boundary around the entire grid (to */
  /* accomodate ghost cells).  */
  sprintf(keyword, "GRID_EDGES");
  if (prm_read_int(fp, keyword, &params->edgef)) {
    if (params->edgef) {
      params->nce1 = nce1 + 2 * params->edgef;
      params->nce2 = nce2 + 2 * params->edgef;
    }
  }

  sprintf(keyword, "PORUS_PLATE");
  if (prm_read_char(fp, keyword, params->reef_frac)) {
    params->porusplate = 1;
    params->atr += 2;
  }

  /* Window sizes */
  read_window_info(params, fp);

  /* Compatibility */
  read_compatible(params, fp);

  /* Debugging */
  read_debug(params, fp);

  /* Used for writing transport files */
  sprintf(keyword, "INPUT_FILE");
  prm_read_char(fp, keyword, params->idumpname);
  prm_read_char(fp, "TRANS_DATA", params->trans_data);

  /* Time unit (optional). Note; this is read in as a mandatory  */
  /* parameter in sched_init(). */
  sprintf(keyword, "TIMEUNIT");
  if (!prm_read_char(fp, keyword, params->timeunit))
    strcpy(params->timeunit, "seconds since 1990-01-01 00:00:00 +08");

  /* ID number and revision */
  sprintf(keyword, "ID_NUMBER");
  prm_read_double(fp, keyword, &params->runno);
  sprintf(keyword, "REVISION");
  prm_read_char(fp, keyword, params->rev);

  /* Output files (optional) */
  get_output_path(params, fp);
  prm_read_char(fp, "OutputTransport", params->trkey);
  sprintf(keyword, "OutputFiles");
  if (prm_read_char(fp, keyword, buf)) {
    if (fopen(buf, "r") == NULL) {
      params->ndf = atoi(buf);
      params->d_name = (char **)malloc(params->ndf * sizeof(char *));
      params->d_filetype = (char **)malloc(params->ndf * sizeof(char *));
      params->d_tstart = (char **)malloc(params->ndf * sizeof(char *));
      params->d_tinc = (char **)malloc(params->ndf * sizeof(char *));
      params->d_tstop = (char **)malloc(params->ndf * sizeof(char *));
      params->d_vars = (char **)malloc(params->ndf * sizeof(char *));
      params->d_sync = (char **)malloc(params->ndf * sizeof(char *));
      params->d_tunit = (char **)malloc(params->ndf * sizeof(char *));
      for (n = 0; n < params->ndf; n++) {
	params->d_name[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.name",n);
	prm_read_char(fp, keyword, params->d_name[n]);
	
	params->d_filetype[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.filetype",n);
	prm_read_char(fp, keyword, params->d_filetype[n]);
	if (strcmp(params->d_filetype[n], "mom") == 0) params->domom |= RDGRID;
	if (strcmp(params->d_filetype[n], "roms") == 0) params->doroms |= RDGRID;
	
	params->d_tstart[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.tstart",n);
	prm_read_char(fp, keyword, params->d_tstart[n]);
	
	params->d_tinc[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.tinc",n);
	prm_read_char(fp, keyword, params->d_tinc[n]);
	
	params->d_tstop[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
	sprintf(keyword, "file%1.1d.tstop",n);
	prm_read_char(fp, keyword, params->d_tstop[n]);
	
	params->d_vars[n] = (char *)malloc(sizeof(char)*MAXLINELEN);
	sprintf(keyword, "file%1.1d.vars",n);
	prm_read_char(fp, keyword, params->d_vars[n]);
	
	params->d_sync[n] = (char *)malloc(sizeof(char)*MAXLINELEN);
	sprintf(keyword, "file%1.1d.sync",n);
	prm_read_char(fp, keyword, params->d_sync[n]);

	params->d_tunit[n] = (char *)malloc(sizeof(char)*MAXLINELEN);
	sprintf(keyword, "file%1.1d.timeunit",n);
	prm_read_char(fp, keyword, params->d_tunit[n]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the flags to default values */
  params->grid_dt = 0.0;
  params->iratio = 0;
  params->stab = SUB_STEP_NOSURF;
  params->uf = 1e-4;
  params->quad_bfc = 0.003;
  params->rampstart = params->t;
  params->rampend = params->rampstart + 86400;
  params->rampf = WIND|FILEIN|CUSTOM|TIDALH|TIDALC|TIDEBC|FLATHR|INV_BARO;
  params->etamax = 10.0;
  strcpy(params->mct, "25%");
  params->g = 9.81;
  params->spec_heat = 3990;
  air_dens = params->air_dens = 1.225;
  params->ambpress = 100800;
  params->u1vh = 0.0;
  params->u2vh = 0.0;
  params->u1kh = 0.0;
  params->u2kh = 0.0;
  params->hmin = 0.0;
  params->z0 = 0.0;
  params->nfe1 = params->nce1 + 1;
  params->nfe2 = params->nce2 + 1;
  strcpy(params->lenunit, "metre");
  ts_set_default_hashtable_size((params->nce1 + 1) * (params->nce2 +
                                                      1) * 4);
  /* Headers for the input dump */
  strcpy(codeheader, "SHOC default version");
  strcpy(params->prmname, prmname);
  strcpy(params->codeheader, codeheader);
  strcpy(params->grid_name, "SHOC grid");
  sprintf(params->grid_desc, "Automated grid from %s", prmname);
  /* Parameter header (optional) */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader)))
    strcpy(params->parameterheader, "Auto grid");
  strcpy(parameterheader, params->parameterheader);

  /* Mixing scheme (optional) */
  sprintf(keyword, "MIXING_SCHEME");
  if (prm_read_char(fp, keyword, params->mixsc)) {
    params->min_tke = 0.0;
    params->min_diss = 0.0;
    params->vz0 = 0.0;
    params->kz0 = 0.0;
    prm_set_errfn(hd_silent_warn);
    sprintf(keyword, "KZ0");
    prm_read_double(fp, keyword, &params->kz0);
    sprintf(keyword, "VZ0");
    prm_read_double(fp, keyword, &params->vz0);
    prm_set_errfn(hd_silent_warn);
    sprintf(keyword, "KZ_ALPHA");
    prm_read_double(fp, keyword, &params->kz_alpha);
    sprintf(keyword, "VZ_ALPHA");
    prm_read_double(fp, keyword, &params->vz_alpha);
    sprintf(keyword, "ZS");
    if (!(prm_read_double(fp, keyword, &params->zs))) {
      if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
          strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
          strcmp(params->mixsc, "harcourt") == 0 ||
          strcmp(params->mixsc, "k-e") == 0 ||
          strcmp(params->mixsc, "k-w") == 0 ||
          strcmp(params->mixsc, "W88") == 0 ||
          strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
        hd_quit("Vertical mixing requires surface length scale, ZS.\n");
    }
    sprintf(keyword, "LMIN");
    prm_read_double(fp, keyword, &params->Lmin);
    sprintf(keyword, "E");
    prm_read_double(fp, keyword, &params->eparam);
    sprintf(keyword, "WAVE_ALPHA");
    prm_read_double(fp, keyword, &params->wave_alpha);
    sprintf(keyword, "WAVE_B1");
    prm_read_double(fp, keyword, &params->wave_b1);
    sprintf(keyword, "WAVE_HEIGHT_FACT");
    prm_read_double(fp, keyword, &params->wave_hf);
    sprintf(keyword, "MIN_TKE");
    prm_read_double(fp, keyword, &params->min_tke);
    sprintf(keyword, "MIN_DISS");
    prm_read_double(fp, keyword, &params->min_diss);
    if (strcmp(params->mixsc, "k-e") == 0)
      params->atr += 2;
    if (strcmp(params->mixsc, "k-w") == 0)
      params->atr += 2;
    if (strcmp(params->mixsc, "W88") == 0)
      params->atr += 2;
    if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0)
      params->atr += 4;
    if (strcmp(params->mixsc, "harcourt") == 0)
      params->atr += 4;
  } else {
    strcpy(params->mixsc, "k-e");
    params->min_tke = 7.6e-6;
    params->min_diss = 5e-10;
    params->vz0 = 1e-5;
    params->kz0 = 1e-5;
    params->zs = 0.01;
    params->atr += 2;
  }

  /* Density (optional) */
  if (prm_read_char(fp, "CALCDENS", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    params->calc_dens = is_true(fields[0]);
    if (n > 1) strcpy(params->densname, fields[1]);      
  }

  /* Surface height */
  sprintf(params->eta_init, "%c", '\0');
  sprintf(keyword, "SURFACE");
  prm_read_char(fp, keyword, params->eta_init);
  sprintf(params->vel_init, "%c", '\0');
  sprintf(keyword, "VELOCITY");
  prm_read_char(fp, keyword, params->vel_init);

  /* Grid refinement */
  if (prm_read_char(fp, "GRID_REFINEMENT", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    if((params->dozoom = is_true(fields[0]))) {
      read_zoom(params, fp);
      if (n > 1 && strcmp(fields[1], "PRECONDITIONED") == 0)
	params->dozoom |= PRECOND;
    }
  }

  /* Cells explicitly defined as OUTSIDE */
  read_blocks(fp, "NOUTSIDE", &params->noutside, &params->oute1, &params->oute2);
  read_blocks(fp, "NLAND", &params->nland, &params->lande1, &params->lande2);

  /* MOM grid conversion */
  if (prm_read_char(fp, "MOM_CONVERT", params->momfile))
    params->domom |= CDF;

  /* ROMS grid conversion */
  if (prm_read_char(fp, "ROMS_CONVERT", params->romsfile))
    params->doroms |= CDF;

  /* Get the Z to Sigma levels scaling factor for ROMS */
  if (params->romsfile) {
    params->roms_z2s = 1.0;
    prm_read_double(fp, "ROMS_Z2S_FACTOR", &params->roms_z2s);
  }

  /* Forcing files */
  prm_set_errfn(hd_silent_warn);
  sprintf(params->wind, "%c", '\0');
  prm_read_char(fp, "WIND_TS", params->wind);
  prm_get_time_in_secs(fp, "WIND_INPUT_DT", &params->wind_dt);
  if (prm_read_int(fp, "NSTORM", &params->nstorm)) {
    prm_set_errfn(hd_quit);
    prm_get_time_in_secs(fp, "STORM_INPUT_DT", &params->storm_dt);
  }
  params->wind_scale = 1.0;
  params->dlv0 = 10.0;
  params->dlv1 = 26.0;
  params->dlc0 = 0.00114;
  params->dlc1 = 0.00218;
  params->wind_type = SPEED;
  if (prm_read_char(fp, "WIND_TYPE", buf)) {
    if(strcmp(buf, "STRESS") == 0)
      params->wind_type = STRESS;
    else
      params->wind_type = SPEED;
  }
  params->stress_fn = ORIGINAL;

  sprintf(params->patm, "%c", '\0');
  sprintf(keyword, "PRESSURE");
  prm_read_char(fp, keyword, params->patm);
  sprintf(keyword, "PRESSURE_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->patm_dt);

  /* Heatflux (optional) */
  sprintf(keyword, "HEATFLUX");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NET_HEAT") == 0) {
      params->heatflux = NET_HEAT;
      sprintf(params->swr, "%c", '\0');
      params->albedo = -1;
      sprintf(keyword, "RADIATION");
      prm_read_char(fp, keyword, params->swr);
      sprintf(keyword, "RADIATION_INPUT_DT");
      prm_get_time_in_secs(fp, keyword, &params->swr_dt);
      prm_read_double(fp, "ALBEDO", &params->albedo);
      sprintf(params->hf, "%c", '\0');
      prm_set_errfn(hd_quit);
      sprintf(keyword, "HEATFLUX_FILE");
      prm_read_char(fp, keyword, params->hf);
      sprintf(keyword, "HEATFLUX_DT");
      prm_get_time_in_secs(fp, keyword, &params->hf_dt);
      prm_set_errfn(hd_silent_warn);
      sprintf(keyword, "HEATFLUX_TC");
      prm_get_time_in_secs(fp, keyword, &params->hftc);
      params->ntrS += 2;
    }
    sprintf(keyword, "HEATFLUX_RAMP");
    if (prm_read_char(fp, keyword, buf)) {
      tm_scale_to_secs(buf, &params->hf_ramp);
    }
    else
      params->hf_ramp = params->t;
  }

  /* SWR parameters (optional) */
  read_swr(params, fp, 1);

  params->coriolis = NULL;
  sprintf(keyword, "CORIOLIS");
  prm_read_darray(fp, keyword, &params->coriolis, &params->nvals);

  /* Sigma (optional) */
  sprintf(keyword, "SIGMA");
  if (prm_read_char(fp, keyword, buf)) {
    params->sigma = is_true(buf);
    if (params->sigma) {
      params->smagorinsky = 0.1;
      params->atr++;
    }
  }

  /* Robustness (optional) */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);
  if (params->robust == 0) {
    params->roammode = A_ROAM_R2;
    /*
    params->nwindows = 4;
    read_window_info(params, fp);
    strcpy(params->dp_mode, "openmp");
    */
  }

  /* Pycnocline sharpening */
  sprintf(keyword, "SHARP_PYC");
  if (prm_read_char(fp, keyword, buf))
    params->sharp_pyc = is_true(buf);

  /* Smagorinsky (optional)
  sprintf(keyword, "SMAGORINSKY");
  prm_read_double(fp, keyword, &params->smagorinsky);
  if (params->smagorinsky > 0.0)
    params->atr++;
  */
  if (params->runmode & DUMP)
    params->sigma=0;
  else {
    if (prm_read_char(fp, "SHOW_LAYERS", buf)) {
      params->show_layers = is_true(buf);
      if (params->show_layers)
	params->atr++;
    }
  }
  sprintf(keyword, "SMAG_SMOOTH");
  prm_read_int(fp, keyword, &params->smag_smooth);

  /* Fatal instability checking (optional) */
  sprintf(keyword, "FATAL");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->fatal = NONE;
    } else {
      params->fatal = 0;
      if (contains_token(buf, "ETA") != NULL)
	params->fatal |= ETA_A;
      if (contains_token(buf, "VEL2D") != NULL)
	params->fatal |= VEL2D;
      if (contains_token(buf, "VEL3D") != NULL)
	params->fatal |= VEL3D;
      if (contains_token(buf, "WIND") != NULL)
	params->fatal |= WIND;
      if (contains_token(buf, "T/S") != NULL)
	params->fatal |= TS;
      if (contains_token(buf, "NAN") != NULL)
	params->fatal |= NANF;
    }
  }

  /* Diagnostics (optional)                                          */
  prm_set_errfn(hd_silent_warn);
  /* Alert tracking                                                  */
  sprintf(keyword, "ALERT");
  if (prm_read_char(fp, keyword, params->alert)) {
    if (contains_token(params->alert, "NONE") != NULL)
      memset(params->alert, 0, sizeof(params->alert));
    if (contains_token(params->alert, "ACTIVE") != NULL)
      params->ntrS += 4;
  }

  sprintf(keyword, "VELMAX");
  prm_read_double(fp, keyword, &params->velmax);
  sprintf(keyword, "VELMAX_2D");
  if(!prm_read_double(fp, keyword, &params->velmax2d))
    params->velmax2d = params->velmax;
  sprintf(keyword, "ETA_DIFF");
  prm_read_double(fp, keyword, &params->etadiff);
  sprintf(keyword, "WMAX");
  prm_read_double(fp, keyword, &params->wmax);
  sprintf(keyword, "ALERT_DT");
  prm_read_char(fp, keyword, params->alert_dt);
  /* CFL diagnostic */
  sprintf(keyword, "CFL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->cfl = NONE;
    if (strcmp(buf, "PASSIVE") == 0)
      params->cfl = PASSIVE;
    if (strcmp(buf, "PASSIVE|WVEL") == 0)
      params->cfl = PASSIVE | WVEL;
    if (strcmp(buf, "ACTIVE3D") == 0)
      params->cfl = ACTIVE3D;
    if (strcmp(buf, "ACTIVE3D|WVEL") == 0)
      params->cfl = ACTIVE3D | WVEL;
    if (strcmp(buf, "ACTIVE") == 0)
      params->cfl = ACTIVE;
    if (strcmp(buf, "ACTIVE|WVEL") == 0)
      params->cfl = ACTIVE | WVEL;
    if (params->cfl & (ACTIVE | ACTIVE3D)) {
      sprintf(keyword, "CFL_DT");
      if (!(prm_read_char(fp, keyword, params->cfl_dt)))
	strcpy(params->cfl_dt, params->stop_time);
    }
  }
  /* CFL */
  if (!(params->cfl & NONE))
    params->ntrS += 5;
  /* Momentum tendency (optional) */
  sprintf(keyword, "MOM_TEND");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->tendf = is_true(buf))) {
      params->atr += 13;
    }
  }
  /* Tracer tendency (optional) */
  sprintf(keyword, "TRA_TEND");
  if (prm_read_char(fp, keyword, params->trtend)) {
    if (strcmp(params->trtend, "NONE") != 0)
      params->atr += 3;
    else
      sprintf(params->trtend, "%c", '\0');
  }
  /* AVHRR SST */
  if (prm_read_char(fp, "AVHRR", params->avhrr_path)) {
    create_avhrr_list(params);
    params->avhrr = 1;
    params->ntrS++;
  }
  /* GHRSST SST */
  if (prm_read_char(fp, "GHRSST", params->ghrsst_path)) {
    create_ghrsst_list(params);
    params->ntrS++;
  }

  /* Transport mode files (optional) */
  read_tmode_params(fp, params);

  /* Means (optional) */
  read_means(params, fp, 1);

  /* Fluxes (optional) */
  sprintf(keyword, "CALC_FLUXES");
  if (prm_read_char(fp, keyword, params->trflux)) {
    if (strcmp(params->trflux, "NONE") != 0) {
      params->atr += 4;
    }
  } else 
    sprintf(params->trflux, "NONE");

  /* Particle tracking */
  if (prm_read_char(fp, "PT_InputFile", params->ptinname)) {
    params->do_pt = 1;
    params->atr++;
  }

  /* Diagnistic numbers (optional) */
  params->atr += numbers_init(params);

  params->atr += read_dhw(params, fp);

  /* River flow diagnostic */
  if (check_river_bdry(params, fp)) {
    params->riverflow = 1;
    params->ntrS++;
    if (check_bdry_options(params, fp) & OP_DYNAHC) {
      params->riverflow = 2;
      params->ntrS++;
      params->atr++;
    }
  }

  /* Bathymetry checks (optional) */
  sprintf(keyword, "MAXGRAD");
  prm_read_double(fp, keyword, &params->maxgrad);
  sprintf(keyword, "VEL_BATHY");
  prm_read_double(fp, keyword, &params->bvel);
  sprintf(keyword, "SMOOTHING");
  prm_read_int(fp, keyword, &params->smooth);
  sprintf(keyword, "SMOOTH_VARS");
  prm_read_char(fp, keyword, params->smooth_v);
  sprintf(keyword, "SCALE_VARS");
  prm_read_char(fp, keyword, params->scale_v);
  /* Bathymetry statistics */
  sprintf(keyword, "BATHY_STATS");
  if (prm_read_char(fp, keyword, params->bathystats))
    params->ntrS += 4;

  /* Time steps (optional) */
  sprintf(keyword, "DT");
  prm_get_time_in_secs(fp, keyword, &params->grid_dt);
  sprintf(keyword, "IRATIO");
  prm_read_int(fp, keyword, &params->iratio);
  if ((params->iratio + 1) % 2 == 0) {
    params->iratio++;
    hd_warn("Only even iratio values allowed : setting iratio=%d\n",
            params->iratio);
  }

  /* Horizontal mixing (optional) */
  read_hdiff(params, fp, 1);

  /* Minimum surface layer thickness (optional) */
  sprintf(keyword, "HMIN");
  prm_read_double(fp, keyword, &params->hmin);

  /* Bottom roughness (optional) */
  z0_init(params);
  /*
  sprintf(keyword, "Z0");
  prm_read_double(fp, keyword, &params->z0);
  */

  /* Surface relaxation (optional) */
  read_eta_relax(params, fp);
  read_vel_relax(params, fp);

  /* Sediment layers (optional) */
  /* Note : sediment layers are read in assuming the first layer is */
  /* closest to the water column. This is reversed when copying to */
  /* the master so that the sediment origin is the deepest layer.  */
  params->sednz = 0;
  if (prm_read_darray(fp, "LAYERFACES_SED", &params->gridz_sed,
                      &params->sednz)) {
    if (--params->sednz < 1)
      hd_quit("Number of sediment layers must be 2 or more\n");
  } else {
    double *dz_sed = NULL;
    if (prm_read_int(fp, "NSEDLAYERS", &params->sednz)) {
      dz_sed = d_alloc_1d(params->sednz);
      if (prm_read_darray(fp, "NSEDLAYERS", &dz_sed, &params->sednz)) {
        if (params->sednz) {
          params->gridz_sed = d_alloc_1d(params->sednz + 1);

          params->gridz_sed[params->sednz] = 0.0;
          for (m = params->sednz - 1; m >= 0; m--) {
            params->gridz_sed[m] = params->gridz_sed[m + 1] -
              dz_sed[params->sednz - m - 1];
	  }
        }
      }
      if(dz_sed)
  d_free_1d(dz_sed);
    }
  }

#if defined(HAVE_SEDIMENT_MODULE)
  read_sediments(params, fp, &params->atr);
#endif

#if defined(HAVE_WAVE_MODULE)
  params->do_wave = 0;
  sprintf(keyword, "DO_WAVES");
  if (prm_read_char(fp, keyword, buf) && is_true(buf)) {
    params->do_wave = is_true(buf);
    prm_get_time_in_secs(fp, "WAVES_DT", &params->wavedt);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  read_ecology(params, fp, &params->atr);
#endif

  /* Waves */
  read_waves(params, fp, 1);
  if (params->waves & STOKES_DRIFT && params->tendf) params->atr+=2;

  /* netCDF dumpfile name */
  if (params->runmode & AUTO) {
    sprintf(keyword, "INPUT_FILE");
    if (!prm_read_char(fp, keyword, params->oname))
      strcpy(params->oname, "auto");
  } else {
    strcpy(params->oname, oname);
  }

  /* Output time unit (optional) */
  sprintf(keyword, "OUTPUT_TIMEUNIT");
  if (prm_read_char(fp, keyword, params->output_tunit)) {
    char *p = strstr(params->output_tunit, "since");
    if (p == NULL)
      sprintf(params->output_tunit, "%s %s",
              strdup(params->output_tunit),
              strstr(params->timeunit, "since"));
  } else
    strcpy(params->output_tunit, params->timeunit);

  /*-----------------------------------------------------------------*/
  /* Tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "NTRACERS");
  /*UR added*/
  if(!prm_read_int(fp, keyword, &n)) {
    if(prm_read_char(fp, "TRACERFILE", buf)) {
      fptr  = fopen(buf,"r");
      if(fptr == NULL) {
	hd_warn("External tracer list cannot be read: %s ",buf);
	fptr = fp;
      }
    }else
      fptr = fp;
  }

  /* Reset the pre-tracer ROAM parameterisation if required */
  if (params->runmode & ROAM) {
    if (prm_read_char(fp, "ROAMv", buf)) {
      if (strcmp(buf, "1") == 0 || strcmp(buf, "CPD") == 0)
	params->roammode = A_ROAM_CPD1;
      else if (strcmp(buf, "2") == 0 || strcmp(buf, "RMD") == 0)
	params->roammode = A_ROAM_CPD2;
      else if (strcmp(buf, "3") == 0 || strcmp(buf, "FLA") == 0)
	params->roammode = A_ROAM_FLA;
      else if (strcmp(buf, "4") == 0 || strcmp(buf, "ROAMv1") == 0)
	params->roammode = A_ROAM_R1;
      else if (strcmp(buf, "5") == 0 || strcmp(buf, "ROAMv2") == 0)
	params->roammode = A_ROAM_R2;
      else if (strcmp(buf, "ROAMv3") == 0)
	params->roammode = A_ROAM_R3;
      else if (strcmp(buf, "6") == 0 || strcmp(buf, "RECOMv1") == 0)
	params->roammode = A_RECOM_R1;
      else if (strcmp(buf, "7") == 0 || strcmp(buf, "RECOMv2") == 0)
	params->roammode = A_RECOM_R2;
    }

    if (params->roammode == A_RECOM_R2)
      auto_params_recom_pre2(fp, params);
    else if (params->roammode == A_RECOM_R1)
      auto_params_recom_pre1(fp, params);
    else if (params->roammode == A_ROAM_R2)
      auto_params_roam_pre2(fp, params);
    else if (params->roammode == A_ROAM_R3)
      auto_params_roam_pre3(fp, params);
    else
      auto_params_roam_pre1(fp, params);
  }
  read_trfilter(params, fp);

  /* Set up tracers */
  sprintf(keyword, "NTRACERS");
  if (prm_read_int(fptr, keyword, &n)) {
    params->atr += 2;   /* Temp and salt always auto */
    params->ntr = params->atr;
    tracer_setup(params, fptr);
    create_tracer_3d(params);
  } else {
    params->ntr = 2 + params->atr;
    params->atr = 0;
    sz = sizeof(tracer_info_t) * params->ntr;
    params->trinfo_3d = (tracer_info_t *)malloc(sz);
    memset(params->trinfo_3d, 0, sz);
    create_tracer_3d(params);
  }
  for (n = 0; n < params->atr; n++) {
    params->trinfo_3d[n].m = -1;
    params->trinfo_3d[n].n = n;
  }

  /* Tracer relaxation */
  params->trrlxn = malloc(params->ntr * sizeof(char *));
  params->trrlxdt = malloc(params->ntr * sizeof(double));
  params->trrlxtc = malloc(params->ntr * sizeof(char *));
  for (n = 0; n < params->ntr; n++) {
    char rdt[MAXSTRLEN];
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    params->trrlxn[n] = malloc(MAXSTRLEN * sizeof(char));
    params->trrlxtc[n] = malloc(MAXSTRLEN * sizeof(char));
    sprintf(params->trrlxn[n], "%c", '\0');
    sprintf(params->trrlxtc[n], "%c", '\0');
    sprintf(keyword, "TRACER%1.1d.relaxation_file", tm);
    if (prm_read_char(fptr, keyword, params->trrlxn[n]) > 0) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TRACER%1.1d.relaxation_input_dt", tm);
      if (prm_read_char(fptr, keyword, rdt) > 0) {
	tm_scale_to_secs(rdt, &params->trrlxdt[n]);
	sprintf(keyword, "TRACER%1.1d.relaxation_time_constant", tm);
	prm_read_char(fptr, keyword, params->trrlxtc[n]);
      }
      prm_set_errfn(hd_silent_warn);
    }
    sprintf(keyword, "relax_%s", tracer->name);
    if (prm_read_char(fptr, keyword, buf)) {
      char fname[MAXSTRLEN], dt[MAXSTRLEN], unit[MAXSTRLEN];
      double in_dt;
      if (sscanf(buf,"%s %s %lf %s", fname, dt, &in_dt, unit) == 4) {
	strcpy(params->trrlxn[n], fname);
	sprintf(rdt,"%f %s\n",in_dt, unit);
	tm_scale_to_secs(rdt, &params->trrlxdt[n]);
	strcpy(params->trrlxtc[n], dt);
      } else
	hd_quit("relax_%s format = filename.nc dt.ts n days\n");      
    }
    strcpy(params->trinfo_3d[n].relax_file, params->trrlxn[n]);
    strcpy(params->trinfo_3d[n].r_rate, params->trrlxtc[n]);
    strcpy(params->trinfo_3d[n].relax_dt, rdt);
  }
  /* Tracer reset */
  params->trrest = malloc(params->ntr * sizeof(char *));
  params->trrestdt = malloc(params->ntr * sizeof(double));
  for (n = 0; n < params->ntr; n++) {
    int tm = params->trinfo_3d[n].m;
    params->trrest[n] = malloc(MAXSTRLEN * sizeof(char));
    memset(params->trrest[n], 0, sizeof(params->trrest[n]));
    sprintf(keyword, "TRACER%1.1d.reset_file", tm);
    if (prm_read_char(fptr, keyword, params->trrest[n]) > 0) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TRACER%1.1d.reset_dt", tm);
      if (prm_read_char(fptr, keyword, buf) > 0)
	tm_scale_to_secs(buf, &params->trrestdt[n]);
      prm_set_errfn(hd_silent_warn);
    }
  }
  prm_set_errfn(hd_quit);

  params->ntdif_h = params->ntdif_v = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->ntdif_h++;
      if (params->compatible & V5342) {
#if defined(HAVE_SEDIMENT_MODULE)
	if (params->do_sed) {
	  if (strcmp(tracer->name, "salt") == 0 ||
	      strcmp(tracer->name, "temp") == 0)
	    params->ntdif_v++;
	} else
	  params->ntdif_v++;
#else
	params->ntdif_v++;
#endif
      } else {
#if defined(HAVE_SEDIMENT_MODULE)
        if (params->do_sed) {
	  if ((tracer->type & WATER && !(tracer->type & CLOSURE) && !(tracer->type & SEDIMENT)) ||
	      strcmp(tracer->name, "salt") == 0 || strcmp(tracer->name, "temp") == 0 )
	    params->ntdif_v++;
        } else
          if(!(tracer->type & CLOSURE))
	    params->ntdif_v++;
#else
        if(!(tracer->type & CLOSURE)) {
	  params->ntdif_v++;
	}
#endif
      }
    }
  }
  if (params->ntdif_h)
    params->tdif_h = i_alloc_1d(params->ntdif_h);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      params->tdif_h[m] = n;
      m++;
    }
  }
  if (params->ntdif_v)
    params->tdif_v = i_alloc_1d(params->ntdif_v);
  m = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];
    if (!tracer->diagn && tracer->type&WATER && tracer->diffuse) {
      if (params->compatible & V5342) {
#if defined(HAVE_SEDIMENT_MODULE)
	if (params->do_sed) {
	  if (strcmp(tracer->name, "salt") == 0 ||
	      strcmp(tracer->name, "temp") == 0) {
	    params->tdif_v[m] = n;
	    m++;
	  }
	} else {
	  params->tdif_v[m] = n;
	  m++;
	}
#else
	params->tdif_v[m] = n;
	m++;
#endif
      } else {
#if defined(HAVE_SEDIMENT_MODULE)
        if (params->do_sed) {
	  if ((tracer->type & WATER && !(tracer->type & CLOSURE) && !(tracer->type & SEDIMENT) ||
	       strcmp(tracer->name, "salt") == 0 || strcmp(tracer->name, "temp") == 0 )) {
            params->tdif_v[m] = n;
            m++;
	  }
        } else {
	  if(!(tracer->type & CLOSURE)) {
            params->tdif_v[m] = n;
            m++;
	  }
	}
#else
        if(!(tracer->type & CLOSURE)) {
          params->tdif_v[m] = n;
          m++;
	}
#endif
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sediment tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  tracer_read(fp, NULL, SEDIM, hd_quit, hd_warn, hd_silent_warn,
        &params->nsed, &params->trinfo_sed);
  if (params->nsed && !params->sednz)
    hd_quit("Sediment tracers requires NSEDLAYERS > 0\n");

  /*-----------------------------------------------------------------*/
  /* 2D Tracer constants and variables */
  prm_set_errfn(hd_silent_warn);
  params->atrS = params->ntrS;
  memset(keyword, 0, sizeof(keyword));
  tracer_read(fp, NULL, INTER, hd_quit, hd_warn, hd_silent_warn,
              &params->ntrS, &params->trinfo_2d);

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  read_grid(params);

  /*-----------------------------------------------------------------*/
  /* Read the bathymetry */
  /* Note : */
  /* gridx[j][i] = params->x[j*2][i*2]; */
  /* u1x[j][i] = params->x[j*2+1][i*2]; */
  /* u2x[j][i] = params->x[j*2][i*2+1]; */
  /* cellx[j][i] = params->x[j*2+1][i*2+1]; */
  /* h1acell[j][i] = 2.0*params->h1[j*2+1][i*2+1]; */
  /* h1au1[j][i] = 2.0*params->h1[j*2+1][i*2]; */
  /* h2au1[j][i] = 2.0*params->h2[j*2+1][i*2]; */
  /* h1au2[j][i] = 2.0*params->h1[j*2][i*2+1]; */
  /* h2au2[j][i] = 2.0*params->h2[j*2][i*2+1]; */
  /* thetau1[j][i] = params->a1[j*2+1][i*2]; */
  /* thetau2[j][i] = params->a2[j*2][i*2+1] - PI/2.0; */
  params->bathy = NULL;

  /* First check to see if the BATHY parameter exists. If so, then use
   * this in preference to the BATHYFILE parameters. */
  prm_set_errfn(hd_silent_warn);
  /*
  if (prm_read_darray(fp, "BATHY", &params->bathy, &params->nvals) ) {
     if (params->nvals != nce1 * nce2)
       hd_quit("auto_params: incorrect BATHY data : %d\n", params->nvals);
  }
  */

  if(!read_bathy_from_file(fp, params)) {
    if (params->roms_grid_opts & ROMS_GRID_BATHY) {
      params->bathy = read_bathy_from_roms(params->roms_grid, nce1, nce2, NULL);
      params->nvals = nce1 * nce2;
    } else if (prm_read_char(fp, "BATHYFILE", buf)) {
      params->bathy = read_bathy_from_db(fp, params);
      if (params->bathy == NULL)
	hd_quit("auto_params: bathymetry data could not be read!");
      else
      	params->nvals = nce1 * nce2;
    } else if (prm_read_char(fp, "MOM_GRID", buf)) {
      params->bathy = read_bathy_from_mom(buf, nce1, nce2);
      params->nvals = nce1 * nce2;
    } else
      hd_quit("auto_params: No bathymetry data was supplied.");
  }

  /* Mask the batymetry with a coastline if required.*/
  prm_set_errfn(hd_silent_warn);

  /* Read from ROMS first */
  if (params->roms_grid_opts & ROMS_GRID_MASK)
    mask_bathy_from_roms(params->roms_grid, params, NULL);
  else 
    if (prm_read_char(fp, "COASTFILE", buf)) 
      mask_bathy_with_coast(buf, params);
  
  /* Interpolate the bathymetry onto the grid */
  bathy = d_alloc_2d(params->nce1, params->nce2);
  n = 0;
  is = js = 0;
  ie = params->nce1;
  je = params->nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = params->nce1 - params->edgef;
    je = params->nce2 - params->edgef;
  }

  for (j = 0; j < params->nce2; j++)
    for (i = 0; i < params->nce1; i++)
      bathy[j][i] = LANDCELL;
  for (j = js; j < je; j++)
    for (i = is; i < ie; i++)
      bathy[j][i] = -params->bathy[n++];

  prm_set_errfn(hd_silent_warn);
  if (!prm_read_double(fp, "BATHYMAX", &params->bmax)) {
    params->bmax = 0;
    for (c = 0; c < params->nvals; c++) {
      if (params->bathy[c] > params->bmax && params->bathy[c] != -NOTVALID
          && params->bathy[c] <= MAXDEPTH)
        params->bmax = params->bathy[c];
    }
    if (params->bmax < 10.0) {
      n = (int)params->bmax / 2;
      params->bmax = (double)n *2.0 + 1.0;
    } else {
      n = (int)params->bmax / 10;
      params->bmax = (double)n *10.0 + 10.0;
    }
    if (params->bmax > MAXDEPTH)
      params->bmax = MAXDEPTH;
  } else {
    for (c = 0; c < params->nvals; c++) {
      if (params->bathy[c] <= MAXDEPTH && params->bathy[c] > params->bmax)
        params->bathy[c] = params->bmax;
    }
  }
  /*
    for (c = 0; c < params->nvals; c++)
      printf("%d %f\n",c,params->bathy[c]);
  */
  if (!prm_read_double(fp, "BATHYMIN", &params->bmin)) {
    params->bmin = 1e10;
    for (c = 0; c < params->nvals; c++) {
      if (params->bathy[c] < params->bmin && params->bathy[c] > 0)
        params->bmin = params->bathy[c];
    }
  }
  prm_read_char(fp, "MIN_CELL_THICKNESS", params->mct);

  /* Read the vertical grid from the input file if required */
  if (prm_read_darray(fp, "LAYERFACES", &params->layers, &params->nz)) {
    if (--params->nz < 1)
      hd_quit("Number of layers in vertical must be 1 or more\n");
  } else if (prm_read_char(fp, "MOM_GRID", buf)) {
    params->layers = read_mom_layers(buf, &params->nz);
  } else {
    /* Set the vertical grid spacing to the default */
    set_sigma_distrib(params, params->bmin, params->bmax);

    for (n = 0; n < params->nz; n++) {
      i = (int)(params->layers[n] * params->bmax);
      j = (int)(params->layers[n + 1] * params->bmax);
      if (j - i) {
	       if (params->layers[n] * params->bmax <= -100.0) {
	         m = (int)(0.1 * params->layers[n] * params->bmax);
	         d1 = 0.1;
	       }
	       else if (params->layers[n] * params->bmax <= -10.0) {
	         m = (int)(params->layers[n] * params->bmax);
	         d1 = 1.0;
	       }
	       else {
          m = (int)(10 * params->layers[n] * params->bmax);
          d1 = 10.0;
         }
      } else {
        i = (int)(10 * params->layers[n] * params->bmax);
        j = (int)(10 * params->layers[n + 1] * params->bmax);
        if (j - i) {
          m = (int)(10 * params->layers[n] * params->bmax);
          d1 = 10.0;
        } else {
          m = (int)(100 * params->layers[n] * params->bmax);
          d1 = 100.0;
        }
      }
      params->layers[n] = (double)m / d1;
    }
  }

  if (params->hmin == 0.0)
    params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
                                    params->layers[params->nz - 1]));

  /*-----------------------------------------------------------------*/
  /* Explicit mappings                                               */
  read_explicit_maps(params, fp);

  /*-----------------------------------------------------------------*/
  /* Open boundaries */
  /* If boundary info is included in the input file, read and return */
  get_bdry_params(params, fp);
  
  if (prm_read_int(fp, "NBOUNDARIES", &params->nobc)) {
    if (params->nobc == 0) {
      d_free_2d(bathy);
      return (params);
    }
    /* Allocate memory for the open boundary data structure */
    params->open =
      (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);

    for (n = 0; n < params->nobc; n++) {

      /* Allocate memory for the boundary structures */
      params->open[n] = OBC_alloc();
      /* Read the open boudary conditions */
      params->open[n]->ntr = params->ntr;
      params->open[n]->atr = params->atr;
      get_OBC_conds(params, params->open[n], fp, n, params->trinfo_3d);

      /* Set the flags for this open boundary */
      get_obc_list(params->open[n], fp, n, "BOUNDARY");
      prm_set_errfn(hd_quit);
    }
    d_free_2d(bathy);
    /* Read the csr tide model paths if required */
    for (n = 0; n < params->nobc; n++) {
      if(params->open[n]->bcond_ele & TIDALH) {
	sprintf(keyword, "TIDE_CSR_ORTHOWEIGHTS");
	prm_read_char(fp, keyword, params->orthoweights);
	sprintf(keyword, "TIDE_CSR_CON_DIR");
	prm_read_char(fp, keyword, params->nodal_dir);
      }
      if(params->open[n]->bcond_ele & TIDALC ||
	 params->open[n]->bcond_nor & TIDALC ||
	 params->open[n]->bcond_nor2d & TIDALC ||
	 params->open[n]->bcond_tan & TIDALC ||
	 params->open[n]->bcond_tan2d & TIDALC) {
	prm_set_errfn(hd_quit);
	sprintf(keyword, "TIDE_CSR_CON_DIR");
	prm_read_char(fp, keyword, params->nodal_dir);
	sprintf(keyword, "TIDE_CONSTITUENTS");
	prm_read_char(fp, keyword, params->tide_con_file);
      }
    }
    return (params);
  }

  /* Open boundary numbers are assigned in the following order: */
  /* Edges : west (i=0), east (i=nce1), south (j=0) north (j=nce2) */
  /* Interior : checks west, east, south north for each cell over */
  /* all cells in the domain.  */
  /* The numbering in the auto input file should reflect this */
  params->nobc = 0;
  nz = params->nz;

  i = 0;
  for (j = 0; j < params->nce2; j++) {
    if (is_wet(params, bathy[j][i])) {
      params->nobc++;
      bdry |= W_BDRY;
      break;
    }
  }
  i = params->nce1 - 1;
  for (j = 0; j < params->nce2; j++) {
    if (is_wet(params, bathy[j][i])) {
      params->nobc++;
      bdry |= E_BDRY;
      break;
    }
  }
  j = 0;
  for (i = 0; i < params->nce1; i++) {
    if (is_wet(params, bathy[j][i])) {
      params->nobc++;
      bdry |= S_BDRY;
      break;
    }
  }
  j = params->nce2 - 1;
  for (i = 0; i < params->nce1; i++) {
    if (is_wet(params, bathy[j][i])) {
      params->nobc++;
      bdry |= N_BDRY;
      break;
    }
  }
  if (params->nobc > MAXBDRY) 
    hd_quit("Number of boundaries found > allowed maximum (%d).\n", MAXBDRY);

  /* Get open boundaries in the interior */
  nib = i_alloc_1d(MAXBDRY);
  memset(nib, 0, MAXBDRY * sizeof(int));
  mask = i_alloc_2d(params->nce1, params->nce2);
  for (i = 0; i < params->nce1; i++)
    for (j = 0; j < params->nce2; j++)
      mask[j][i] = 0;
  for (i = 1; i < params->nce1 - 1; i++) {
    for (j = 1; j < params->nce2 - 1; j++) {
      ii = i; jj = j;
      while (is_wet(params, bathy[jj][i]) && 
	     is_outside(params, bathy[jj][i-1]) && 
	     !(mask[jj][i-1] & DRY) && jj < params->nce2-1) {
	mask[jj][i-1] |= (DRY|L_EDGE);
	nib[params->nobc]++;
	jj++;
      }
      if (nib[params->nobc]) {
        params->nobc++;
	continue;
      }
      while (is_wet(params, bathy[jj][i]) 
	     && is_outside(params, bathy[jj][i+1]) && 
	     !(mask[jj][i+1] & DRY) && jj < params->nce2-1) {
	mask[jj][i+1] |= (DRY|R_EDGE);
	nib[params->nobc]++;
	jj++;
      }
      if (nib[params->nobc]) {
        params->nobc++;
	continue;
      }
      while (is_wet(params, bathy[j][ii]) &&
	     is_outside(params, bathy[j-1][ii]) && 
	     !(mask[j-1][ii] & DRY) && ii < params->nce1-1) {
	mask[j-1][ii] |= (DRY|B_EDGE);
	nib[params->nobc]++;
	ii++;
      }
      if (nib[params->nobc]) {
        params->nobc++;
	continue;
      }
      while (is_wet(params, bathy[j][ii]) &&
	     is_outside(params, bathy[j+1][ii]) && 
	     !(mask[j+1][ii] & DRY) && ii < params->nce1-1) {
	mask[j+1][ii] |= (DRY|F_EDGE);
	nib[params->nobc]++;
	ii++;
      }
      if (nib[params->nobc]) {
        params->nobc++;
	continue;
      }
    }
  }
  if (params->nobc > MAXBDRY) 
    hd_quit("Number of boundaries found > allowed maximum (%d).\n", MAXBDRY);
  /* Get the river boundaries if explicitly specified using RIVER */
  ib = i_alloc_1d(MAXBDRY);
  jb = i_alloc_1d(MAXBDRY);
  rb = i_alloc_1d(MAXBDRY);
  memset(rb, 0, MAXBDRY * sizeof(int));
  rf = i_alloc_1d(MAXBDRY);
  memset(rf, 0, MAXBDRY * sizeof(int));
  rhc = d_alloc_1d(MAXBDRY);
  nr = m = 0;
  rs = params->nobc;
  sprintf(keyword, "RIVER%1d", nr);
  while(prm_read_char(fp, keyword, buf)) {
    int found = 0;
    double lat, lon;
    char name[MAXSTRLEN], buf1[MAXSTRLEN];
    m++;
    sprintf(keyword, "RIVER%1d", m);
    if (sscanf(buf, "%s %lf %lf %s", name, &lat, &lon, buf1) == 4) {
      int rmode = (params->roammode & (A_RECOM_R1|A_RECOM_R2)) ? 1 : 0;
      if (lat > 90.0 || lat < -90.0) { /* Swap */
	d1 = lon;
	lon = lat;
	lat = d1;
      }
      /*
       * If RECOM, force mode=1 in find_closest_b below as placement
       * of RIVER cell is strictly enforced on the actual grid
       */
      if (find_closest_b(params, bathy, lat, lon, &ii, &jj, rmode)) {
	hd_warn("Can't locate river %s in the grid.\n", name);
	continue;
      } else {
	if (rmode)
	  hd_warn("River %s wet cell located at cell (%d %d).\n", name, ii, jj);
	else
	  hd_warn("River %s dry cell located at cell (%d %d).\n", name, ii, jj);
      }
      /* Set the cell as OUTSIDE if right two cell centres are wet */
      if (is_wet(params, bathy[jj][ii+1]) && 
	  is_wet(params, bathy[jj][ii+2]) && !(mask[jj][ii] & DRY)) {
	rf[nr] |= (U1BDRY|L_EDGE);
	bathy[jj][ii] = NOTVALID;
	reset_bathy(params, ii, jj, NOTVALID);
	rhc[nr] = bathy[jj][ii+1];
	ib[nr] = ii+1; jb[nr] = jj;
	found = 1;
      }
      /* Set the cell as OUTSIDE if left two cell centres are wet */
      else if (is_wet(params, bathy[jj][ii-1]) && 
	       is_wet(params, bathy[jj][ii-2]) && !(mask[jj][ii] & DRY)) {
	rf[nr] |= (U1BDRY|R_EDGE);
	bathy[jj][ii] = NOTVALID;
	reset_bathy(params, ii, jj, NOTVALID);
	rhc[nr] = bathy[jj][ii-1];
	ib[nr] = ii; jb[nr] = jj;
	found = 1;
      }
      /* Set the cell as OUTSIDE if top two cell centres are wet */
      else if (is_wet(params, bathy[jj+1][ii]) &&
	       is_wet(params, bathy[jj+2][ii]) && !(mask[jj][ii] & DRY)) {
	rf[nr] |= (U2BDRY|B_EDGE);
	bathy[jj][ii] = NOTVALID;
	reset_bathy(params, ii, jj, NOTVALID);
	rhc[nr] = bathy[jj+1][ii];
	ib[nr] = ii; jb[nr] = jj+1;
	found = 1;
      }
      /* Set the cell as OUTSIDE if bottom two cell centres are wet */
      else if (is_wet(params, bathy[jj-1][ii]) &&
	       is_wet(params, bathy[jj-2][ii]) && !(mask[jj][ii] & DRY)) {
	rf[nr] |= (U2BDRY|F_EDGE);
	bathy[jj][ii] = NOTVALID;
	reset_bathy(params, ii, jj, NOTVALID);
	rhc[nr] = bathy[jj-1][ii];
	ib[nr] = ii; jb[nr] = jj;
	found = 1;
      }
      if (!found) {
	if (find_closest_b(params, bathy, lat, lon, &ii, &jj, 1)) {
	  hd_warn("Can't locate river %s in the grid.\n", name);
	  continue;
	} else
	  hd_warn("River %s wet cell located at cell (%d %d).\n", name, ii, jj);
	if (is_wet(params, bathy[jj][ii])) {
	  if (ii > 0 && ii < params->nce1-1 && 
	      is_wet(params, bathy[jj][ii+1]) && is_land(params, bathy[jj][ii-1])) {
	    rf[nr] |= (U1BDRY|L_EDGE);
	    bathy[jj][ii-1] = NOTVALID;
	    reset_bathy(params, ii-1, jj, NOTVALID);
	    rhc[nr] = bathy[jj][ii];
	    ib[nr] = ii; jb[nr] = jj;
	    found = 1;
	  } else if (ii > 0 && ii < params->nce1-1 && 
		     is_wet(params, bathy[jj][ii-1]) && is_land(params, bathy[jj][ii+1])) {
	    rf[nr] |= (U1BDRY|R_EDGE);
	    bathy[jj][ii+1] = NOTVALID;
	    reset_bathy(params, ii+1, jj, NOTVALID);
	    rhc[nr] = bathy[jj][ii];
	    ib[nr] = ii; jb[nr] = jj;
	    found = 1;
	  } else if (jj > 0 && jj < params->nce2-1 &&
		     is_wet(params, bathy[jj+1][ii]) && is_land(params, bathy[jj-1][ii])) {
	    rf[nr] |= (U2BDRY|B_EDGE);
	    bathy[jj-1][ii] = NOTVALID;
	    reset_bathy(params, ii, jj-1, NOTVALID);
	    rhc[nr] = bathy[jj][ii];
	    ib[nr] = ii; jb[nr] = jj;
	    found = 1;
	  } else if (jj > 0 && jj < params->nce2-1 && 
		     is_wet(params, bathy[jj-1][ii]) && is_land(params, bathy[jj+1][ii])) {
	    rf[nr] |= (U2BDRY|F_EDGE);
	    bathy[jj+1][ii] = NOTVALID;
	    reset_bathy(params, ii, jj+1, NOTVALID);
	    rhc[nr] = bathy[jj][ii];
	    ib[nr] = ii; jb[nr] = jj;
	    found = 1;
	  }
	}
      }
      if (!found) {
	hd_warn("Can't position river %s in the grid.\n", name);
	continue;
      }
      rb[nr] = m-1;
      nib[params->nobc] = 1;
      params->nobc++;
      nr++;
    }
  }

  if (params->nobc > MAXBDRY) 
    hd_quit("Number of boundaries found > allowed maximum (%d).\n", MAXBDRY);

  /* Allocate memory for the open boundary data structure */
  params->open =
    (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);
  for (n = 0; n < params->nobc; n++) {
    /* Allocate memory for the boundary structures */
    params->open[n] = OBC_alloc();
    open = params->open[n];
    open->id = n;
    open->ntr = params->ntr;
    open->atr = params->atr;

    /* Read the open boudary conditions */
    dataf = 0;
    init_OBC_conds(params, open);

    open->relax_zone_tan = 0;
    open->relax_zone_ele = 0;
    open->bathycon = 0;
    open->relax_time = 0;
    open->relax_timei = 0;
    open->sponge_zone = 0;
    open->sponge_zone_h = 0;
    open->sponge_f = 0;
    open->inverse_barometer = 0;
    open->tidemem_depth = 0.0;
    sprintf(open->tsfn, "%c", '\0');
    open->nbstd = 0;
    open->bstd = (char **)malloc(NBSTD * sizeof(char *));
    for (i = 0; i < NBSTD; i++) {
      open->bstd[i] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      sprintf(open->bstd[i], "%c", '\0');
    }

    sprintf(keyword, "BOUNDARY%1d.DATA", n);
    if (prm_read_char(fp, keyword, open->tsfn))
      dataf = 1;

    open->ntr = params->ntr;
    open->bcond_nor = NOGRAD;
    open->bcond_nor2d = NOGRAD;
    open->relax_zone_nor = 0;
    open->bcond_tan = CLAMPD;
    open->bcond_tan2d = CLAMPD;
    open->bcond_w = NOTHIN;
    if (dataf) {
      open->bcond_ele = FILEIN;
      /*open->relax_time=900.0;*/
    } else {
      open->bcond_ele = NOTHIN;
      /*open->relax_time=900.0;*/
    }
    open->bcond_Vz = NOGRAD;
    open->bcond_Kz = NOGRAD;
    open->bcond_tra = i_alloc_1d(open->ntr);
    open->clampv = d_alloc_1d(open->ntr);
    open->trpc = d_alloc_1d(open->ntr);
    open->relax_zone_tra = i_alloc_1d(open->ntr);
    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if (dataf && ((strcmp(tracer->name, "salt") == 0) ||
		    (strcmp(tracer->name, "temp") == 0)))
        open->bcond_tra[tn] = UPSTRM | FILEIN;
      else
        open->bcond_tra[tn] = NOGRAD;
      open->relax_zone_tra[tn] = 0;
      open->clampv[tn] = 0;
      open->trpc[tn] = 0.0;
    }
  }
  /* Read the csr tide model paths if required */
  for (n = 0; n < params->nobc; n++) {
    if(params->open[n]->bcond_ele & TIDALH) {
      sprintf(keyword, "TIDE_CSR_ORTHOWEIGHTS");
      prm_read_char(fp, keyword, params->orthoweights);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      break;
    }
    if(params->open[n]->bcond_ele & TIDALC) {
      prm_set_errfn(hd_quit);
      sprintf(keyword, "TIDE_CSR_CON_DIR");
      prm_read_char(fp, keyword, params->nodal_dir);
      sprintf(keyword, "TIDE_CONSTITUENTS");
      prm_read_char(fp, keyword, params->tide_con_file);
      break;
    }
  }

  /* Allocate OBC location vectors */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    if (nib[n]) {
      open->iloc = i_alloc_1d(abs(nib[n]));
      open->jloc = i_alloc_1d(abs(nib[n]));
    }
  }

  /* Get open boundaries on the grid edges */
  n = 0;
  if (bdry & W_BDRY) {
    i = 0;
    open = params->open[n];
    open->type = U1BDRY;
    strcpy(open->name, "West");
    open->npts = 0;
    for (j = 0; j < params->nce2; j++) {
      if (is_wet(params, bathy[j][i])) {
        open->npts++;
      }
    }
    open->iloc = i_alloc_1d(open->npts);
    open->jloc = i_alloc_1d(open->npts);
    open->npts = 0;
    for (j = 0; j < params->nce2; j++) {
      if (is_wet(params, bathy[j][i])) {
        open->iloc[open->npts] = i;
        open->jloc[open->npts] = j;
        open->npts++;
      }
    }
    n++;
  }
  if (bdry & E_BDRY) {
    i = params->nce1 - 1;
    open = params->open[n];
    open->type = U1BDRY;
    strcpy(open->name, "East");
    open->npts = 0;
    for (j = 0; j < params->nce2; j++) {
      if (is_wet(params, bathy[j][i])) {
        open->npts++;
      }
    }
    open->iloc = i_alloc_1d(open->npts);
    open->jloc = i_alloc_1d(open->npts);
    open->npts = 0;
    for (j = 0; j < params->nce2; j++) {
      if (is_wet(params, bathy[j][i])) {
        open->iloc[open->npts] = i + 1;
        open->jloc[open->npts] = j;
        open->npts++;
      }
    }
    n++;
  }
  if (bdry & S_BDRY) {
    j = 0;
    open = params->open[n];
    open->type = U2BDRY;
    strcpy(open->name, "South");
    open->npts = 0;
    for (i = 0; i < params->nce1; i++) {
      if (is_wet(params, bathy[j][i])) {
        open->npts++;
      }
    }
    open->iloc = i_alloc_1d(open->npts);
    open->jloc = i_alloc_1d(open->npts);
    open->npts = 0;
    for (i = 0; i < params->nce1; i++) {
      if (is_wet(params, bathy[j][i])) {
        open->iloc[open->npts] = i;
        open->jloc[open->npts] = j;
        open->npts++;
      }
    }
    n++;
  }
  if (bdry & N_BDRY) {
    j = params->nce2 - 1;
    open = params->open[n];
    open->type = U2BDRY;
    strcpy(open->name, "North");
    open->npts = 0;
    for (i = 0; i < params->nce1; i++) {
      if (is_wet(params, bathy[j][i])) {
        open->npts++;
      }
    }
    open->iloc = i_alloc_1d(open->npts);
    open->jloc = i_alloc_1d(open->npts);
    open->npts = 0;
    for (i = 0; i < params->nce1; i++) {
      if (is_wet(params, bathy[j][i])) {
        open->iloc[open->npts] = i;
        open->jloc[open->npts] = j + 1;
        open->npts++;
      }
    }
    n++;
  }

  /* Get open boundaries in the interior */
  m = 0;
  for (i = 0; i < params->nce1; i++)
    for (j = 0; j < params->nce2; j++)
      if (mask[j][i] & DRY) mask[j][i] &= ~DRY;
  for (i = 1; i < params->nce1 - 1; i++) {
    if (n == params->nobc) break;
    for (j = 1; j < params->nce2 - 1; j++) {
      if (n == params->nobc) break;
      open = params->open[n];
      open->npts = 0;
      ii = i; jj = j;
      /* L_EDGE */
      while (is_wet(params, bathy[jj][i]) && 
	     is_outside(params, bathy[jj][i-1]) && 
	     mask[jj][i-1] & L_EDGE && !(mask[jj][i-1] & DRY) &&
	     jj < params->nce2-1) {
	mask[jj][i-1] |= DRY;
        open->type = U1BDRY;
	open->iloc[open->npts] = i;
	open->jloc[open->npts] = jj;
	open->npts++;
        sprintf(open->name, "bdry%1.1d", n);
	jj++;
      }
      if (open->npts) {
	n++;
	continue;
      }
      /* R_EDGE */
      while (is_wet(params, bathy[jj][i]) &&
	     is_outside(params, bathy[jj][i+1]) &&  
	     mask[jj][i+1] & R_EDGE && !(mask[jj][i+1] & DRY) &&
	     jj < params->nce2-1) {
	mask[jj][i+1] |= DRY;
        open->type = U1BDRY;
	open->iloc[open->npts] = i+1;
	open->jloc[open->npts] = jj;
	open->npts++;
        sprintf(open->name, "bdry%1.1d", n);
	jj++;
      }
      if (open->npts) {
	n++;
	continue;
      }
      /* B_EDGE */
      while (is_wet(params, bathy[j][ii]) && 
	     is_outside(params, bathy[j-1][ii]) && 
	     mask[j-1][ii] & B_EDGE && !(mask[j-1][ii] & DRY) &&
	     ii < params->nce1-1) {
	mask[j-1][ii] |= DRY;
        open->type = U2BDRY;
	open->iloc[open->npts] = ii;
	open->jloc[open->npts] = j;
	open->npts++;
        sprintf(open->name, "bdry%1.1d", n);
	ii++;
      }
      if (open->npts) {
	n++;
	continue;
      }
      /* F_EDGE */
      while (is_wet(params, bathy[j][ii]) &&
	     is_outside(params, bathy[j+1][ii]) &&  
	     mask[j+1][ii] & F_EDGE && !(mask[j+1][ii] & DRY) &&
	     ii < params->nce1-1) {
	mask[j+1][ii] |= DRY;
        open->type = U2BDRY;
	open->iloc[open->npts] = ii;
	open->jloc[open->npts] = j+1;
	open->npts++;
        sprintf(open->name, "bdry%1.1d", n);
	ii++;
      }
      if (open->npts) {
	n++;
	continue;
      }
    }
  }

  /* Get the RIVER open boundaries */
  for (m = 0; m < nr; m++) {
    open = params->open[m + rs];
    open->npts = 0;
    if (rf[m] & U1BDRY)
      open->type = U1BDRY;
    else
      open->type = U2BDRY;
    open->iloc[open->npts] = ib[m];
    open->jloc[open->npts] = jb[m];
    open->npts++;
    sprintf(open->name, "river%1.1d", m);
  }

  /* Read the flow data file if explicitly specified */
  for (n = 0; n < params->nobc; n++) {
    int bdf = 0;
    char name[MAXSTRLEN];
    double lat, lon;
    open = params->open[n];
    sprintf(open->bflow, "%c", '\0');
    open->bhc = 0.0;
    sprintf(open->cusname_u1, "%c", '\0');
    sprintf(open->cusname_u2, "%c", '\0');
    for (i = 0; i < open->ntr; i++) {
      open->cusname_t[i] = (char *)malloc(MAXSTRLEN);
      sprintf(open->cusname_t[i], "%c", '\0');
    }	    
    open->custype = 0;

    switch (open->type) {
    case U1BDRY:
      open->custype = U1BDRY;
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u1");
      prm_read_char(fp, keyword, open->cusname_u1);
      break;

    case U2BDRY:
      open->custype = U2BDRY;
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u2");
      prm_read_char(fp, keyword, open->cusname_u2);
      break;

    case (U1BDRY | U2BDRY):
      open->custype = U1BDRY | U2BDRY;
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u1");
      if (prm_read_char(fp, keyword, open->cusname_u1))
        bdf = 1;
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n, "u2");
      if (prm_read_char(fp, keyword, open->cusname_u2)) {
        if (bdf)
          hd_quit
            ("Multiple custom routines (u1 and u2) specified for BOUNDARY%1d.\n",
             n);
        else {
          bdf = 1;
        }
      }
      break;
    }
    if (strcmp(open->cusname_u1, "u1flowbdry") == 0) {
      sprintf(keyword, "BOUNDARY%1d.U1_FLOW", n);
      prm_read_char(fp, keyword, open->bflow);
      sprintf(keyword, "BOUNDARY%1d.U1_HC", n);
      prm_read_double(fp, keyword, &open->bhc);
      open->custype = U1BDRY;
    }
    if (strcmp(open->cusname_u2, "u2flowbdry") == 0) {
      sprintf(keyword, "BOUNDARY%1d.U2_FLOW", n);
      prm_read_char(fp, keyword, open->bflow);
      sprintf(keyword, "BOUNDARY%1d.U2_HC", n);
      prm_read_double(fp, keyword, &open->bhc);
      open->custype = U2BDRY;
    }
    if (strlen(open->cusname_u1) || strlen(open->cusname_u2)) {
      open->bcond_ele = NOTHIN;
      open->bcond_nor = CUSTOM;
      open->bcond_nor2d = VERTIN;
      open->bcond_tra[0] = UPSTRM | FILEIN;
      open->bcond_tra[1] = UPSTRM | FILEIN;
    }

    /* Set the river OBCs if specified via RIVER */
    nr = -1;
    for (m = 0; m < params->nobc - rs; m++) {
      if (ib[m] == open->iloc[0] && jb[m] == open->jloc[0]) {
	nr = m;
	break;
      }
    }
    if (nr >= 0) {
      m = nr;
      nr = rb[m];
      sprintf(keyword, "RIVER%1d", nr);
      if (prm_read_char(fp, keyword, buf)) {
	if (sscanf(buf, "%s %lf %lf %s", name, &lat, &lon, open->bflow) == 4) {
	  strcpy(open->name, name);
	  if (rf[m] & U1BDRY) {
	    strcpy(open->cusname_u1, "u1flowbdry");
	    open->custype = U1BDRY;
	  }
	  if (rf[m] & U2BDRY) {
	    strcpy(open->cusname_u2, "u2flowbdry");
	    open->custype = U2BDRY;
	  }
	  open->bhc = rhc[m];
	  open->bcond_ele = NOTHIN;
	  open->bcond_nor = CUSTOM;
	  open->bcond_nor2d = VERTIN;
	  open->bcond_tan = CLAMPD;
	  open->bcond_tan2d = CLAMPD;
	  for (tn = 0; tn < params->ntr; tn++) {
	    tracer_info_t *tracer = &params->trinfo_3d[tn];
	    if (strcmp(tracer->name, "salt") == 0) {
	      open->bcond_tra[tn] = TRFLUX|CUSTOM;
	      if (!strlen(open->cusname_t[tn]))
		sprintf(open->cusname_t[tn],"0.0");
	    }
	    if (strcmp(tracer->name, "temp") == 0)
	      open->bcond_tra[tn] = TRCONC|FILEIN;
	  }
	  open->bgz = laux;
	}
      }
    }

    /* Tracers */
    for (i = 0; i < open->ntr; i++) {
      /* 
       * cusname_t should already be set to empty upon malloc (see above)
       * Therefore if not found it won't be overwritten - otherwise
       * salt 0.0 was being nulled out
       */
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.%s", n,
              params->trinfo_3d[i].name);
      prm_read_char(fp, keyword, open->cusname_t[i]);

      open->scale_s[i] = (char *)malloc(MAXSTRLEN);
      sprintf(open->scale_s[i], "%c", '\0');
      sprintf(keyword, "BOUNDARY%1d.SCALE_S.%s", n, params->trinfo_3d[i].name);
      prm_read_char(fp, keyword, open->scale_s[i]);

      open->scale_p[i] = (char *)malloc(MAXSTRLEN);
      sprintf(open->scale_p[i], "%c", '\0');
      sprintf(keyword, "BOUNDARY%1d.SCALE_P.%s", n, params->trinfo_3d[i].name);
      prm_read_char(fp, keyword, open->scale_p[i]);
    }
  }

  /* Reset the post-tracer ROAM parameterisation if required */
  if (params->runmode & ROAM) {
    if (params->roammode == A_ROAM_CPD1)   /* Standard ROAM */
      auto_params_roam_post1(fp, params);
    if (params->roammode == A_ROAM_CPD2)   /* ETA_RELAX ramp, RAYMND bcond_nor */
      auto_params_roam_post2(fp, params);
    if (params->roammode == A_ROAM_FLA)   /* ETA_RELAX ramp, FLATHR */
      auto_params_roam_post3(fp, params);
    if (params->roammode == A_ROAM_R1)   /* ETA_RELAX ramp, velocity forced */
      auto_params_roam_post4(fp, params);
    if (params->roammode == A_ROAM_R2)   /* A_ROAM_R2 with alternative robust parameterisations */
      auto_params_roam_post5(fp, params);
    if (params->roammode == A_ROAM_R3)   /* A_ROAM_R3 with alternative robust parameterisations */
      auto_params_roam_post5(fp, params);
    if (params->roammode == A_RECOM_R1)   /* RECOM */
      auto_params_recom_post1(fp, params);
    if (params->roammode == A_RECOM_R2)   /* RECOM + ROBUST */
      auto_params_recom_post2(fp, params);
  }

  d_free_2d(bathy);
  i_free_2d(mask);
  i_free_1d(ib);
  i_free_1d(jb); 
  i_free_1d(rb);
  i_free_1d(rf);
  i_free_1d(nib);
  return (params);
}

/* END auto_params()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void autoset(parameters_t *params, master_t *master)
{
  geometry_t *geom = master->geom;
  int n, cc, c, lc;             /* Sparse counters */
  int i, j, k;                  /* Cartesian counters */
  double spd = 2.0;             /* Maximum advective velocity expected */
  double cif = 2.0;             /* Scaling factor for internal wave speed */
  double sf = 0.8;              /* Safety factor for cfl calculations */
  double hf = 0.05;             /* Factor for horizontal diffusion */
  double sfr = 1.0;             /* Depth dependent safety factor */
  double eps = 1e5;             /* Large value for cfl calculations */
  double cfl2d, cfl3d;          /* CFL timesteps */
  double hmax;                  /* Maximum horizontal diffusion */
  double lat, lon;              /* Latitude and longitude for Coriolis */
  double d1, d2, d3;            /* Dummies */
  dump_data_t *dumpdata = master->dumpdata;

  /*-----------------------------------------------------------------*/
  /* Time */
  master->t = dumpdata->t = params->t;

  /*-----------------------------------------------------------------*/
  /* Vertical grid geometry. The bathymetry is set in build_sparse_map() */
  /* and overwritten in dump_read() if no automation is performed.  */
  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;                 /* Continue on ghosts */
    if (i < params->nfe1 && j < params->nfe2) {
      geom->gridz[c] = params->layers[k];
      geom->cellz[c] = 0.5 * (params->layers[k] + params->layers[k + 1]);
    }
  }

  /* Sediments */
  /* Note : sediment layers are read in assuming the first layer is */
  /* closest to the water column. This is reversed here so that the */
  /* sediment origin is the deepest layer.  */
  if (geom->sednz) {
    for (c = 1; c <= geom->sgnumS; c++) {
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
        continue;               /* Continue on ghosts */
      if (i < params->nce1 && j < params->nce2) {
        for (lc = 0; lc <= geom->sednz; lc++) {
          geom->gridz_sed[lc][c] = params->gridz_sed[lc];
        }
        for (lc = geom->sednz - 1; lc >= 0; lc--) {
          geom->cellz_sed[lc][c] = 0.5 * (geom->gridz_sed[lc][c] +
                                          geom->gridz_sed[lc + 1][c]);
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Cell mappings for output */
  neighbour_none(dumpdata);

  /*-----------------------------------------------------------------*/
  /* Horizontal grid geometry */
  memset(geom->h1au1, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2au2, 0, geom->sgsizS * sizeof(double));
  memset(geom->h1acell, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2acell, 0, geom->sgsizS * sizeof(double));
  memset(geom->h2au1, 0, geom->sgsizS * sizeof(double));
  memset(geom->h1au2, 0, geom->sgsizS * sizeof(double));
  memset(geom->thetau1, 0, geom->sgsizS * sizeof(double));
  memset(geom->thetau2, 0, geom->sgsizS * sizeof(double));
  for (j = 0; j < geom->nce2 + 1; j++) {
    /* Corner points */
    for (i = 0; i < geom->nce1 + 1; i++) {
      dumpdata->gridx[j][i] = params->x[j * 2][i * 2];
      dumpdata->gridy[j][i] = params->y[j * 2][i * 2];
      dumpdata->h1agrid[j][i] = 2.0 * params->h1[j * 2][i * 2];
      dumpdata->h2agrid[j][i] = 2.0 * params->h2[j * 2][i * 2];
      c = geom->map[params->nz - 1][j][i];
      if (c) {
        c = geom->m2d[c];
        geom->gridx[c] = params->x[j * 2][i * 2];
        geom->gridy[c] = params->y[j * 2][i * 2];
      }
    }
    /* u2 points */
    for (i = 0; i < params->nce1; i++) {
      dumpdata->u2x[j][i] = params->x[j * 2][i * 2 + 1];
      dumpdata->u2y[j][i] = params->y[j * 2][i * 2 + 1];
      dumpdata->h1au2[j][i] = 2.0 * params->h1[j * 2][i * 2 + 1];
      dumpdata->h2au2[j][i] = 2.0 * params->h2[j * 2][i * 2 + 1];
      dumpdata->thetau2[j][i] = params->a2[j * 2][i * 2 + 1] - PI / 2.0;
      c = geom->map[params->nz - 1][j][i];
      if (c) {
        c = geom->m2d[c];
        geom->u2x[c] = params->x[j * 2][i * 2 + 1];
        geom->u2y[c] = params->y[j * 2][i * 2 + 1];
        geom->h1au2[c] = 2.0 * params->h1[j * 2][i * 2 + 1];
        geom->h2au2[c] = 2.0 * params->h2[j * 2][i * 2 + 1];
        geom->thetau2[c] = params->a2[j * 2][i * 2 + 1] - PI / 2.0;
      }
    }
  }
  for (j = 0; j < params->nce2; j++) {
    /* u1 points */
    for (i = 0; i < params->nce1 + 1; i++) {
      dumpdata->u1x[j][i] = params->x[j * 2 + 1][i * 2];
      dumpdata->u1y[j][i] = params->y[j * 2 + 1][i * 2];
      dumpdata->h1au1[j][i] = 2.0 * params->h1[j * 2 + 1][i * 2];
      dumpdata->h2au1[j][i] = 2.0 * params->h2[j * 2 + 1][i * 2];
      dumpdata->thetau1[j][i] = params->a1[j * 2 + 1][i * 2];
      c = geom->map[params->nz - 1][j][i];
      if (c) {
        c = geom->m2d[c];
        geom->u1x[c] = params->x[j * 2 + 1][i * 2];
        geom->u1y[c] = params->y[j * 2 + 1][i * 2];
        geom->h1au1[c] = 2.0 * params->h1[j * 2 + 1][i * 2];
        geom->h2au1[c] = 2.0 * params->h2[j * 2 + 1][i * 2];
        geom->thetau1[c] = params->a1[j * 2 + 1][i * 2];
      }
    }
    /* Cell points */
    for (i = 0; i < params->nce1; i++) {
      dumpdata->cellx[j][i] = params->x[j * 2 + 1][i * 2 + 1];
      dumpdata->celly[j][i] = params->y[j * 2 + 1][i * 2 + 1];
      dumpdata->h1acell[j][i] = 2.0 * params->h1[j * 2 + 1][i * 2 + 1];
      dumpdata->h2acell[j][i] = 2.0 * params->h2[j * 2 + 1][i * 2 + 1];
      c = geom->map[params->nz - 1][j][i];
      if (c) {
        c = geom->m2d[c];
        geom->cellx[c] = params->x[j * 2 + 1][i * 2 + 1];
        geom->celly[c] = params->y[j * 2 + 1][i * 2 + 1];
        geom->h1acell[c] = 2.0 * params->h1[j * 2 + 1][i * 2 + 1];
        geom->h2acell[c] = 2.0 * params->h2[j * 2 + 1][i * 2 + 1];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Load the tracers. Note hd_ts_read() is used to read the tracers */
  /* from file and this requires time variables to be set on the */
  /* master. These are only set in compute_constants(), called after */
  /* this routine, so preset them here.  */
  master->t = master->tstart = params->t;
  strcpy(master->timeunit, params->timeunit);
  strcpy(master->output_tunit, params->output_tunit);
  load_tracer_step_3d(params, master, params->prmfd);
  load_tracer_step_2d(params, master, params->prmfd);

  /*-----------------------------------------------------------------*/
  /* Get the CFL conditions */
  /* Get the initial density distribution */
  density_m(master);

  if (params->grid_dt == 0.0) {
    double btws, btwsm = 0.0;
    double bcws, bcwsm = 0.0;

    cfl2d = cfl3d = 1e10;
    if (params->robust >= 6)
      sf *= (1.0 - 0.5 * ((double)(params->robust) - 6.0) / 4.0);

    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = lc = geom->w2_t[cc];

      d1 = sqrt(1.0 / (geom->h1acell[c] * geom->h1acell[c]) +
                1.0 / (geom->h2acell[c] * geom->h2acell[c]));
      bcws = int_wave_spd_cont_m(master, c, cc);
      bcwsm = (bcws > bcwsm) ? bcws : bcwsm;
      btws = sqrt(-params->g * geom->botz[c]);
      btwsm = (btws > btwsm) ? btws : btwsm;
      d2 = 2.0 * cif * bcws + spd;
      d3 = 2.0 * btws + spd;
      d3 = (d3 > 0.0) ? 1.0 / (d1 * d3) : eps;
      if (d3 < cfl2d)
        cfl2d = d3;
      d2 = (d2 > 0.0) ? 1.0 / (d1 * d2) : eps;
      if (d2 < cfl3d)
        cfl3d = d2;
    }

    if (DEBUG("init_m")) {
      dlog("init_m : auto dt","Maximum advective velocity expected = %3.1f",spd);
      dlog("init_m : auto dt","Internal wave speed scaling factor= %3.1f",cif);
      dlog("init_m : auto dt","Safety factor for cfl calculations = %3.1f",sf);
      dlog("init_m : auto dt","Maximum barotropic wave speed = %3.1f",btwsm);
      dlog("init_m : auto dt","Maximum baroclinic wave speed = %3.1f",bcwsm);
      dlog("init_m : auto dt","cfl3d = %5.1f, cfl2d = %3.1f",cfl3d,cfl2d);
    }

    if (params->bmax > 3000)
      sfr = 0.35;
    else if (params->bmax > 1000 && params->bmax <= 3000)
      sfr = 0.5;
    sf *= sfr;
    if (sf * cfl3d > 300)
      d1 = 30;
    else
      d1 = 5;
    c = (int)cfl3d *sf / d1;
    params->grid_dt = d1 * (double)c;
    c = (int)cfl2d *sf;
    while (c > 1 && (c % 2 != 0 && c % 3 != 0 && c % 5 != 0))
      c--;
    if (c)
      cfl2d = (double)c;
    else {
      if (params->roammode == A_RECOM_R2) {
	c = (int)(1000 * cfl2d);
	cfl2d = (double)c / 1000.0;
      } else
	cfl2d = floor(cfl2d);
    }
    if (cfl2d <= 1)
      hd_warn("2D time-step is less than 1 (%f)\n", cfl2d);
    params->iratio = (int)(params->grid_dt / cfl2d);
    if (params->iratio % 2 != 0)
      params->iratio++;
    params->grid_dt = (double)params->iratio * cfl2d;
    if (DEBUG("init_m")) {
      dlog("init_m : auto dt","Scaled cfl3d = %5.1f",params->grid_dt);
      dlog("init_m : auto dt","Scaled cfl2d = %5.1f",cfl2d);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Horizontal mixing coefficients */
  n = 0;
  d1 = d2 = 1e-10;
  hmax = 1e10;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    d1 += geom->h1acell[c];
    d2 += geom->h2acell[c];
    d3 = 1.0 / ((1.0 / (geom->h1acell[c] * geom->h1acell[c]) +
                 1.0 / (geom->h2acell[c] * geom->h2acell[c])) * 4.0 *
                params->grid_dt);
    if (d3 < hmax)
      hmax = d3;
    n += 1;
  }
  if (hmax > 5.0) {
    c = (int)hmax / 5;
    hmax = 5.0 * (double)c;
  }

  if (n) {
    int u1khf = 0, u1vhf = 0, u2khf = 0, u2vhf = 0;
    if (params->u1kh == 0.0)
      u1khf = 1;
    if (params->u1vh == 0.0)
      u1vhf = 1;
    if (params->u2kh == 0.0)
      u2khf = 1;
    if (params->u2vh == 0.0)
      u2vhf = 1;
    d1 /= (double)n;
    d2 /= (double)n;
    d1 = 0.01 * d1 * d1 / params->grid_dt;
    d2 = 0.01 * d2 * d2 / params->grid_dt;
    if (d1 > 5.0) {
      c = (int)d1 / 5;
      if (u1khf)
	params->u1kh = 5.0 * (double)c;
      if (u1vhf)
	params->u1vh = 5.0 * (double)c;
    } else {
      if (u1khf)
	params->u1kh = d1;
      if (u1vhf)
	params->u1vh = d1;
    }
    if (d2 > 5.0) {
      c = (int)d2 / 5;
      if (u2khf)
	params->u2kh = 5.0 * (double)c;
      if (u2vhf)
	params->u2vh = 5.0 * (double)c;
    } else {
      if (u2khf)
	params->u2kh = d2;
      if (u2vhf)
	params->u2vh = d2;
    }
    /* Set limits */
    if (u1khf) {
      if (params->u1kh > hmax)
        params->u1kh = hmax;
      if (hf * hmax > 5.0 && params->u1kh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u1kh = 5.0 * (double)c;
      }
    }
    if (u1vhf) {
      /* Linear from u1vh *= 0.5 for robust = 6, u1vh *= 1 for       */
      /* robust = 10                                                 */
      /*
      if (params->robust >= 6)
	      params->u1vh *= (((double)params->robust - 10.0) / 8.0 + 1.0);
      */
      if (params->u1vh > hmax)
        params->u1vh = hmax;
      if (hf * hmax > 5.0 && params->u1vh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u1vh = 5.0 * (double)c;
      }
    }
    if (u2khf) {
      if (params->u2kh > hmax)
        params->u2kh = hmax;
      if (hf * hmax > 5.0 && params->u2kh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u2kh = 5.0 * (double)c;
      }
    }
    if (u2vhf) {
      /* Linear from u2vh *= 0.5 for robust = 6, u2vh *= 1 for       */
      /* robust = 10                                                 */
      /*
      if (params->robust >= 6)
	      params->u2vh *= (((double)params->robust - 10.0) / 8.0 + 1.0);
      */
      if (params->u2vh > hmax)
        params->u2vh = hmax;
      if (hf * hmax > 5.0 && params->u2vh < hf * hmax) {
        c = (int)(hf * hmax) / 5;
        params->u2vh = 5.0 * (double)c;
      }
    }
  } else
    hd_warn("Horizontal mixing is set to zero\n");
  if(params->smagorinsky > 0.0) {
    params->u1vh *= -1.0;
    params->u2vh *= -1.0;
  }

  /*-----------------------------------------------------------------*/
  /* Bottom roughness */
  if (params->z0 == 0.0) {
    n = 0;
    d1 = 0.0;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->zp1[geom->bot_t[cc]];
      lc = geom->m2d[c];
      if (params->sigma)
        d1 += 0.5 * ((-geom->gridz[c] * geom->botz[lc]) - geom->botz[lc]);
      else
        d1 += 0.5 * (geom->gridz[c] - geom->botz[lc]);
      n++;
    }
    if (n)
      d1 /= (double)n;
    n = 1e4 * d1 / (exp(0.4 * sqrt(1.0 / params->quad_bfc)) - 1);
    params->z0 = (double)n / 1e4;
  }

  /*-----------------------------------------------------------------*/
  /* Coriolis */
  memset(master->coriolis, 0, geom->sgsizS * sizeof(double));
  if (params->coriolis) {
    n = 0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++) {
        dumpdata->coriolis[j][i] = params->coriolis[n++];
       c = geom->map[params->nz - 1][j][i];
        if (c) {
          c = geom->m2d[c];
          master->coriolis[c] = dumpdata->coriolis[j][i];
        }
      }
  } else if (strlen(projection) > 0) {
    map_proj_t *mp = NULL;
    prm_set_errfn(hd_silent_warn);
    if (strncasecmp(projection, GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) !=
        0) {
      char *args[256];
      int nargs = parseline(strdup(projection), args, 256);
      mp = mp_init(nargs, args);
    }
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++) {
        lon = dumpdata->cellx[j][i];
        lat = dumpdata->celly[j][i];
        if (mp != NULL)
          mp_inverse(mp, lon, lat, &lat, &lon);
        dumpdata->coriolis[j][i] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
        c = geom->map[params->nz - 1][j][i];
        if (c) {
          c = geom->m2d[c];
          master->coriolis[c] = dumpdata->coriolis[j][i];
        }
      }
  } else {
    lat = -30.0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++) {
        dumpdata->coriolis[j][i] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
        c = geom->map[params->nz - 1][j][i];
        if (c) {
          c = geom->m2d[c];
          master->coriolis[c] = dumpdata->coriolis[j][i];
        }
      }
  }

  /*-----------------------------------------------------------------*/
  /* Initialize */
  memset(master->u1av, 0, geom->sgsizS * sizeof(double));
  memset(master->u2av, 0, geom->sgsizS * sizeof(double));
  memset(master->u1bot, 0, geom->sgsizS * sizeof(double));
  memset(master->u2bot, 0, geom->sgsizS * sizeof(double));
  memset(master->wind1, 0, geom->sgsizS * sizeof(double));
  memset(master->wind2, 0, geom->sgsizS * sizeof(double));
  memset(master->wtop, 0, geom->sgsizS * sizeof(double));
  memset(master->topz, 0, geom->sgsizS * sizeof(double));
  memset(master->eta, 0, geom->sgsizS * sizeof(double));
  memset(master->patm, 0, geom->sgsizS * sizeof(double));
  memset(master->u1, 0, geom->sgsiz * sizeof(double));
  memset(master->u2, 0, geom->sgsiz * sizeof(double));
  memset(master->w, 0, geom->sgsiz * sizeof(double));
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    master->Vz[c] = params->vz0;
    master->Kz[c] = params->kz0;
  }
  /* Initialise the surface elevation if required */
  eta_init(geom, params, master);
  /* Initialise the velocity if required */
  vel_init(geom, params, master);

  /* Get the roms grid structure if required                         */
  dumpdata->romsgrid = convert_roms_grid(params, dumpdata);

}

/* END autoset()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract bathymetry from the database bathy file        */
/*-------------------------------------------------------------------*/
static double *read_bathy_from_db(FILE *fp, parameters_t *params) {

  char bfname[MAXFNAMELEN];
  double **src_topo = NULL;
  double *dst_bathy = NULL;
  int has_proj = (strlen(params->projection) > 0);
  int is_geog = has_proj && (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
  int i, j, n;
  int is, js, ie, je, jm, im;
  int nce1 = params->nce1;
  int nce2 = params->nce2;
  int hint = 1;        /* Interpolation type : 0 = AUTO    */
                       /*                      1 = NEAREST */
                       /*                      2 = LINEAR  */
                       /*                      3 = NATURAL */
                       /*                      4 = INVERSE */
                       /*                      5 = AVERAGE */

  /* Adjust for grid edges */
  is = js = 0;
  ie = nce1;
  je = nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = nce1 - params->edgef;
    je = nce2 - params->edgef;
    nce1 -= 2 * params->edgef;
    nce2 -= 2 * params->edgef;
  }

  prm_set_errfn(hd_quit);

  /* Read the bathymetry data. */
  if (prm_read_char(fp, "BATHYFILE", bfname)) {
    topo_files_t *tf;
    double **x = d_alloc_2d(nce1+1, nce2+1);
    double **y = d_alloc_2d(nce1+1, nce2+1);

    /* Reduce double grid to single grid corners */

    if (is_geog) {
      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i - params->edgef;
          jm = j - params->edgef;
          x[jm][im] = params->x[2*j][2*i];
          y[jm][im] = params->y[2*j][2*i];
        }
      }
    } else if (has_proj) {
      char *args[256];
      int nargs = parseline(strdup(params->projection), args, 256);
      map_proj_t *mp = mp_init(nargs, args);

      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i - params->edgef;
          jm = j - params->edgef;
          mp_inverse(mp, params->x[2*j][2*i], params->y[2*j][2*i],
              &y[jm][im], &x[jm][im]);
        }
      }
    }

    // FR: 28-05-2018 Disabling to allow for models to run across 0 deg longitude
    // Reverted 14-06-2018
    /* Check ranges; must be 0 - 360 */
    for (j = 0; j < nce2 + 1; j++) {
      for (i = 0; i < nce1 + 1; i++) {
	while (x[j][i] < 0.0 || x[j][i] > 360.0) {
	  if (x[j][i] < 0.0) x[j][i] += 360.0;
	  if (x[j][i] > 360.0) x[j][i] -= 360.0;
	}
      }
    }

    /* Populate the grid */
    tf = topo_simple_init(bfname);

    src_topo = topo_get_z_for_grid(tf, x, y, nce1, nce2, (topo_hint_t)hint);

    d_free_2d(y);
    d_free_2d(x);
    topo_destroy(tf);

    if (src_topo == NULL)
      hd_quit("read_bathy_from_db: unable to populate bathymetric grid.");
  }

  /* Convert 2D bathy grid to a 1d bathymetry array. */
  dst_bathy = d_alloc_1d(nce1 * nce2);
  n = 0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++, n++) {
      dst_bathy[n] = -src_topo[j][i];
      if (dst_bathy[n] <= 0)
         dst_bathy[n] = -LANDCELL;
    }
  }

  if (src_topo != NULL)
    d_free_2d(src_topo);

  return dst_bathy;
}

/* END read_bathy_from_db()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract bathymetry from the database bathy file        */
/*-------------------------------------------------------------------*/
static int read_bathy_from_file(FILE *fp, parameters_t *params) {
  FILE* fb = fp;
  char buf[MAXLINELEN];
  int is, js, ie, je, jm, im;
  int nce1 = params->nce1;
  int nce2 = params->nce2;

  /* Adjust for grid edges */
  is = js = 0;
  ie = nce1;
  je = nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = nce1 - params->edgef;
    je = nce2 - params->edgef;
    nce1 -= 2 * params->edgef;
    nce2 -= 2 * params->edgef;
  }

  if(!prm_read_char(fb, "BATHY", buf)) {
    if(prm_read_char(fb, "BATHYEXTFILE", buf))
      fb=fopen(buf,"r");
    else {
      hd_warn("No bathymetry provided!?");
      return 0;
    }
  }	
  prm_read_darray(fb, "BATHY", &params->bathy, &params->nvals);
  params->bmin = 0;
  prm_read_double(fb, "BATHYMIN", &params->bmin);
  params->bmax = 0;
  prm_read_double(fb, "BATHYMAX", &params->bmax);
	
  if(fb != fp)
    fclose(fb);

  if (params->nvals != nce1 * nce2)
    hd_quit("auto_params: incorrect BATHY data : %d\n", params->nvals);

  return 1;
}
	

/*-------------------------------------------------------------------*/
/* Routine to extract bathymetry from the database bathy file        */
/*-------------------------------------------------------------------*/
static void mask_bathy_with_coast(const char *cfname, parameters_t *params) {

  coast_file_t* cf = coast_open(cfname);
  if (cf != NULL) {
    double minlon = 360, maxlon = -360;
    double minlat = 90, maxlat = -90;
    double **x, **y;
    int has_proj = (strlen(params->projection) > 0);
    int is_geog = has_proj && (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
    int i, j, n;
    int is, js, ie, je, jm, im;
    int nce1 = params->nce1;
    int nce2 = params->nce2;
    coast_tile_iterator_t *ti;
    coast_tile_t *tile;

    /* Adjust for grid edges */
    is = js = 0;
    ie = nce1;
    je = nce2;
    if (params->edgef) {
      is = js = params->edgef;
      ie = nce1 - params->edgef;
      je = nce2 - params->edgef;
      nce1 -= 2 * params->edgef;
      nce2 -= 2 * params->edgef;
    }

    /* Reduce double grid to single grid corners */
    x = d_alloc_2d(nce1+1, nce2+1);
    y = d_alloc_2d(nce1+1, nce2+1);

    if (is_geog) {
      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i - params->edgef;
          jm = j - params->edgef;
          x[jm][im] = params->x[2*j][2*i];
          y[jm][im] = params->y[2*j][2*i];
        }
      }
    } else if (has_proj) {
      char *args[256];
      int nargs = parseline(strdup(params->projection), args, 256);
      map_proj_t *mp = mp_init(nargs, args);

      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i - params->edgef;
          jm = j - params->edgef;
          mp_inverse(mp, params->x[2*j][2*i], params->y[2*j][2*i],
                 &y[jm][im], &x[jm][im]);
        }
      }
    }


    /* Read the coastline data if available */

    for (j=0; j<nce2; j++) {
      for (i=0; i<nce1; i++, ++n) {
        if (x[j][i] < minlon) minlon = x[j][i];
        if (x[j][i] > maxlon) maxlon = x[j][i];
        if (y[j][i] < minlat) minlat = y[j][i];
        if (y[j][i] > maxlat) maxlat = y[j][i];
      }
    }

    /* Load a tile at a time, and compare the lat/lons */
    ti = coast_get_tile_iterator(cf, minlon, maxlon, minlat, maxlat);
    while((tile = ti->next_tile(ti)) != NULL) {
      for (i=0; i<tile->npolys; ++i) {
        poly_t *p = tile->polys[i];
        if (tile->level[i] != 1) continue;
        n = 0;
        for (j=0; j<nce2; j++) {
          for (i=0; i<nce1; i++, ++n) {
            int ii, jj, land = 0;
            if (params->bathy[n] == -LANDCELL) continue;

            /* Build a 9 x 9 grid. If more than more than 4 cells land,
             * then assume the cell is LAND. */
            for (jj=0; jj<=2; ++jj) {
              for (ii=0; ii<=2; ++ii) {
                double nx =
                    (x[j][i]*(1.0-ii/2.0) + x[j][i+1]*(ii/2.0))
                            * (1.0-jj/2.0)
                  + (x[j+1][i]*(1.0-ii/2.0) + x[j+1][i+1]*(ii/2.0))
                            * (jj/2.0);
                double ny =
                       (y[j][i]*(1.0-jj/2.0) + y[j+1][i]*(jj/2.0))
                            * (1.0-ii/2.0)
                  + (y[j][i+1]*(1.0-jj/2.0) + y[j+1][i+1]*(jj/2.0))
                            * (ii/2.0);
                if (poly_contains_point(p, nx, ny))
                   ++land;
              }
            }

            if (land > 4)
               params->bathy[n] = -LANDCELL;
          }
        }
      }

      ti->free_tile(tile);
    }

    coast_free_tile_iterator(ti);

    d_free_2d(y);
    d_free_2d(x);

    coast_close(cf);
  }
}

/* END mask_bathy_with_coast()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise surface elevation.                          */
/*-------------------------------------------------------------------*/
void eta_init(geometry_t *geom,      /* Global geometry              */
	      parameters_t *params,  /* Input parameter data         */
	      master_t *master       /* Master data                  */
  )
{
  int c, cc;
  int i, j, n;
  int number;
  char buf[MAXSTRLEN];
  char key[MAXSTRLEN];
  double val;

  if (strlen(params->eta_init)) {

    strcpy(buf, params->eta_init);
    if (sscanf(buf,"%s %d %d", key, &i, &j) == 3) {
      if (strcmp(key, "BARO_VORTEX") == 0)
	baro_vortex(master, geom, i, j);
      return;
    }

    value_init_2d(master, master->eta, params->prmfd, params->eta_init,
		  "eta", "SURFACE", 0.0, 0); 

    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      val = DRY_FRAC * params->hmin + geom->botz[c];
      master->eta[c] = max(val, master->eta[c]);
    }
    if (params->means & ETA_M)
      memcpy(master->etam, master->eta, geom->sgsizS * sizeof(double));
  }

  sprintf(buf, "SURFACE_POINTS");
  if (prm_read_int(params->prmfd, buf, &n)) {
    prm_flush_line(params->prmfd);
    for (cc = 0; cc < n; cc++) {
      if (fscanf(params->prmfd, "%d %d : %lf", &i, &j, &val) != 3)
	hd_quit("eta_init: Can't read i j : val in SURFACE_POINTS list.\n", n);
      /* Flush the remainder of the line */
      prm_flush_line(params->prmfd);
      c = geom->map[geom->nz - 1][j][i];
      if (c > 0 && c < geom->sgnum)
	master->eta[c] = val;
    }
  }

  /* Check for NaNs                                                  */
  s2c_2d(geom, master->eta, dumpdata->eta, nce1, nce2);
  for (cc = 1; cc <= geom->b2_t; cc++) {
    int ci, cj;
    c = geom->w2_t[cc];
    if (isnan(master->eta[c])) {
      find_closest_nonnan(dumpdata, dumpdata->eta, geom->s2i[c], geom->s2j[c],
			  geom->nz-1, &ci, &cj);
      master->eta[c] = dumpdata->eta[cj][ci];
      hd_warn("eta_init: Found NaN at (%d %d); replaced with (%d %d).\n", 
	      geom->s2i[c], geom->s2j[c], ci, cj);
    }
  }
}

/* END eta_init()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialises the model with a baroclinic vortex around (io,jo)     */
/* From Penven et al (2006), Ocean Modelling, 12, 157-187            */
/*-------------------------------------------------------------------*/
void baro_vortex(master_t *master, geometry_t *geom, int io, int jo)
{
  int c, cs, cc, i, j;
  int co;
  double rho0 = 1024.0;
  double H1 = 2500.0; /* Depth of no motion                          */
  double ef = 60.0;   /* e-folding length scale                      */
  double umax = 1.0;  /* Maximum surface geostrophic velocity        */
  double N = 0.003, N2;;
  double g = 9.81;
  double Po, P;
  double z, dx, dy, f;
  double t1, t2, t3, d1, d2, d3;
  double rho;
  double fo = 2.0 * M_PI / 86400.0;
 
  co = geom->map[geom->nz-1][jo][io];
  ef *= 1e3;
  N2 = N * N;
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    z = geom->cellz[c];
    i = geom->s2i[c];
    j = geom->s2j[c];
    dx = (double)abs(i - io) * geom->h1acell[cs];
    dy = (double)abs(j - jo) * geom->h2acell[cs];
    f = 2 * fo * sin(geom->celly[cs] * M_PI / 180.0);
    Po = fabs(rho0 * f * umax * ef * sqrt(exp(1.0)));

    t1 = rho0 * N2 * z / g;
    t3 = exp(-(dx * dx + dy * dy) / (2.0 * ef * ef));
    if (z > -H1) {
      t2 = (1.0 - exp(-(z+H1))) / (H1 - 1.0 + exp(-H1));
      rho = rho0 - t1 - Po * t2 * t3 / g;
      P = Po * t3 * (H1 - 1.0 + z + exp(-(z+H1))) / (H1 - 1.0 + exp(-H1));
    } else {
      rho = rho0 - t1;
      P = 0.0;
    }
    master->temp[c] = (1030 - rho) / 0.28;
    master->dens[c] = rho;
    if(geom->s2k[c] == geom->nz-1)
      master->eta[cs] = Po * t3 / (master->dens[cs] * g);
  }
  /*calc_density(master);*/
  density_m(master);
}

/* END baro_vortex()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise velocity.                                   */
/*-------------------------------------------------------------------*/
void vel_init(geometry_t *geom,      /* Global geometry              */
	      parameters_t *params,  /* Input parameter data         */
	      master_t *master       /* Master data                  */
  )
{
  int c, cc, cs;
  int tu1, tu2;
  char buf[MAXSTRLEN];
  double u, v;
  double thetau1, thetau2;

  if (strlen(params->vel_init)) {
    timeseries_t *ts;
    int idu, idv;

    if (strcmp(params->vel_init, "GEOSTROPHIC") == 0) return;
    /* Read from file                                            */
    /* hd_ts_read() is used to read eta from file and this       */
    /* requires time variables to be set on the master. These    */
    /* are only set in compute_constants(), called after this    */
    /* routine, so preset them here.                             */
    master->t = params->t;
    strcpy(master->timeunit, params->timeunit);
    strcpy(master->output_tunit, params->output_tunit);

    ts = hd_ts_read(master, params->vel_init, 0);
    idu = ts_get_index(ts, fv_get_varname(params->vel_init, "u", buf));
    tu1 = 1;
    if (idu < 0) {
      idu = ts_get_index(ts, fv_get_varname(params->vel_init, "u1", buf));
      tu1 = 2;
    }
    if (idu < 0)
      hd_quit("vel_init: The file '%s' does not contain the velocity 'u'.\n",
	      params->vel_init);
    idv = ts_get_index(ts, fv_get_varname(params->vel_init, "v", buf));
    tu2 = 1;
    if (idv < 0) {
      idv = ts_get_index(ts, fv_get_varname(params->vel_init, "u2", buf));
      tu2 = 2;
    }
    if (idv < 0)
      hd_quit("vel_init: The file '%s' does not contain the velocity 'v'.\n",
	      params->vel_init);

    if (tu1 == 1 && tu2 == 1) {
      for (cc = 1; cc <= geom->b3_e1; cc++) {
	c = geom->w3_e1[cc];
	cs = geom->m2d[c];
	thetau1 = geom->thetau1[cs];
	u = ts_eval_xyz(ts, idu, master->t, geom->u1x[cs], geom->u1y[cs],
			geom->cellz[c]);
	v = ts_eval_xyz(ts, idv, master->t, geom->u1x[cs], geom->u1y[cs],
			geom->cellz[c]);
	master->u1[c] = cos(thetau1) * u + sin(thetau1) * v;
      }
      for (cc = 1; cc <= geom->b3_e2; cc++) {
	c = geom->w3_e2[cc];
	cs = geom->m2d[c];
	thetau2 = geom->thetau2[cs];
	u = ts_eval_xyz(ts, idu, master->t, geom->u2x[cs], geom->u2y[cs],
			geom->cellz[c]);
	v = ts_eval_xyz(ts, idv, master->t, geom->u2x[cs], geom->u2y[cs],
			geom->cellz[c]);
	master->u2[c] = cos(thetau2) * v - sin(thetau2) * u;
      }
    } else if (tu1 == 2 && tu2 == 2) {
      for (cc = 1; cc <= geom->b3_e1; cc++) {
	c = geom->w3_e1[cc];
	cs = geom->m2d[c];
	master->u1[c] = ts_eval_xyz(ts, idu, master->t, geom->u1x[cs], geom->u1y[cs],
				    geom->cellz[c]);
      }
      for (cc = 1; cc <= geom->b3_e2; cc++) {
	c = geom->w3_e2[cc];
	cs = geom->m2d[c];
	master->u2[c] = ts_eval_xyz(ts, idv, master->t, geom->u2x[cs], geom->u2y[cs],
				    geom->cellz[c]);
      }
    } else
      hd_quit("vel_init: Incompatible velocity components.\n");
  }
}

/* END vel_init()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise diagnostic numbers                          */
/*-------------------------------------------------------------------*/
int numbers_init(parameters_t *params      /* Input parameter data   */
  )
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  int ntr = 0;

  params->numbers = 0;
  sprintf(keyword, "NUMBERS");
  if (prm_read_char(params->prmfd, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->numbers = NONE;
    } else {
      if (contains_token(buf, "BRUNT") != NULL) {
	params->numbers |= BRUNT;
	ntr++;
      }
      if (contains_token(buf, "INT_WAVE") != NULL) {
	params->numbers |= INT_WAVE;
	ntr++;
      }
      if (contains_token(buf, "RICHARDSON_GR") != NULL) {
	params->numbers |= RICHARD_GR;
	ntr++;
      }
      if (contains_token(buf, "RICHARDSON_FL") != NULL) {
	params->numbers |= RICHARD_FL;
	ntr++;
      }
      if (contains_token(buf, "REYNOLDS") != NULL) {
	params->numbers |= REYNOLDS;
	ntr++;
      }
      if (contains_token(buf, "FROUDE") != NULL) {
	params->numbers |= FROUDE;
	ntr++;
      }
      if (contains_token(buf, "SIGMA_T") != NULL) {
	params->numbers |= SIGMA_T;
	ntr++;
      }
      if (contains_token(buf, "SOUND") != NULL) {
	params->numbers |= SOUND;
	ntr+=2;
	params->ntrS++;
      }
      if (contains_token(buf, "SHEAR_V") != NULL) {
	params->numbers |= SHEAR_V;
	ntr++;
      }
      if (contains_token(buf, "BUOY_PROD") != NULL) {
	params->numbers |= BUOY_PROD;
	ntr++;
      }
      if (contains_token(buf, "SHEAR_PROD") != NULL) {
	params->numbers |= SHEAR_PROD;
	ntr++;
      }
      if (contains_token(buf, "ROSSBY_IN") != NULL) {
	params->numbers |= ROSSBY_IN;
	ntr++;
      }
      if (contains_token(buf, "ROSSBY_EX") != NULL) {
	params->numbers |= ROSSBY_EX;
	params->ntrS++;
      }
      if (contains_token(buf, "SPEED_3D") != NULL) {
	params->numbers |= SPEED_3D;
	ntr += 2;
      }
      if (contains_token(buf, "SPEED_2D") != NULL) {
	params->numbers |= SPEED_2D;
	params->ntrS += 2;
      }
      if (contains_token(buf, "SPEED_SQ") != NULL) {
	params->numbers |= SPEED_SQ;
	params->ntrS++;
      }
      if (contains_token(buf, "ENERGY") != NULL) {
	params->numbers |= ENERGY;
	ntr++;
      }
      if (contains_token(buf, "KINETIC") != NULL) {
	params->numbers |= KINETIC;
	ntr++;
      }
      if (contains_token(buf, "WET_CELLS") != NULL) {
	params->numbers |= WET_CELLS;
	params->ntrS++;
      }
      if (contains_token(buf, "SLOPE") != NULL) {
	params->numbers |= SLOPE;
	params->ntrS+=2;
      }
      if (contains_token(buf, "BOTSTRESS") != NULL) {
	params->numbers |= BOTSTRESS;
	params->ntrS+=3;
      }
      if (contains_token(buf, "SURF_LAYER") != NULL) {
	params->numbers |= SURF_LAYER;
	params->ntrS++;
      }
      if (contains_token(buf, "WIND_CD") != NULL) {
	params->numbers |= WIND_CD;
	params->ntrS++;
      }
      if (contains_token(buf, "OBC_PHASE") != NULL) {
	params->numbers |= OBC_PHASE;
	params->ntrS++;
      }
      if (contains_token(buf, "DUMMIES") != NULL) {
	params->numbers |= DUMMIES;
	ntr+=3;
      }
      if (contains_token(buf, "UNIT") != NULL) {
	params->numbers |= UNIT;
	ntr++;
      }
      if (contains_token(buf, "EKMAN_PUMP") != NULL) {
	params->numbers |= EKPUMP;
	params->ntrS+=2;
      }
      if (contains_token(buf, "ALL_NUMBERS") != NULL) {
	params->numbers |= (BRUNT|INT_WAVE|RICHARD_GR|RICHARD_FL|REYNOLDS|
			    FROUDE|SOUND|ROSSBY_IN|ROSSBY_EX|SHEAR_V|
			    BUOY_PROD|SHEAR_PROD|SPEED_3D|SPEED_2D|
			    OBC_PHASE|SPEED_SQ|SIGMA_T|ENERGY|WET_CELLS|
			    SLOPE|SURF_LAYER|BOTSTRESS|UNIT|EKPUMP);
	ntr+=17;
	params->ntrS += 15;
      }
    }
  } else {
    params->numbers = NONE;
  }
  return(ntr);
}

/* END numbers_init()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initializes the degree heating diagnostic                         */
/*-------------------------------------------------------------------*/
int read_dhw(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n, ret;

  sprintf(keyword, "DHW");

  if (prm_read_char(fp, keyword, buf)) {
    n = parseline(buf, fields, MAXNUMARGS);
    if (n == 1) {
      strcpy(params->dhw, fields[0]);
      ret = 3;
      params->dhwf = DHW_RT;
    }
    if (n == 2) {
      strcpy(params->dhw, fields[0]);
      strcpy(params->dhdf, fields[1]);
      sprintf(keyword, "DHW_DT");
      if (prm_read_char(fp, keyword, buf))
	tm_scale_to_secs(buf, &params->dhw_dt);
      else
	params->dhw_dt = 86400.0;
      params->dhwf = DHW_NOAA;
      ret = 3;
    }
    return(ret);
  }
  return(0);
}

/* END init_dhw()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise bottom roughness.                           */
/*-------------------------------------------------------------------*/
void z0_init(parameters_t *params    /* Input parameter data         */
  )
{
  int n1;

  prm_set_errfn(hd_quit);
  params->z0s = NULL;
  prm_read_double(params->prmfd, "Z0", &params->z0);

  /* Check if it is a list of values in the parameter file          */
  if ((int)params->z0 == params->nce1 * params->nce2) {
    if (prm_read_darray(params->prmfd, "Z0", &params->z0s, &n1) > 0) {
      if (n1 != params->nce1 * params->nce2)
	hd_quit("z0_init : incorrect number specified for Z0 parameter.\n");
    }
  }
}

/* END z0_init()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise a value from a number, array or file input  */
/*-------------------------------------------------------------------*/
double *value_init(master_t *master,      /* Master data             */
		   parameters_t *params,  /* Input parameter data    */
		   char *keyword,         /* Variable keyword        */
		   char *vname            /* Variable name           */
  )
{
  geometry_t *geom = master->geom;
  int c, cc;
  int i, j, n;
  int number, nvals;
  char buf[MAXSTRLEN];
  double val, *d1, **d2;
  double *ret = NULL;

  if (strlen(keyword)) {

    ret = d_alloc_1d(geom->sgsizS);

    /* Check if it is a list of values in the parameter file         */
    if (sscanf(keyword, "%lf", &val) == 1) {
      val = atof(keyword);
      /* Check if it is a list of values in the parameter file          */
      if ((int)val == params->nce1 * params->nce2) {
	if (prm_read_darray(params->prmfd, vname, &d1, &nvals) > 0) {
	  if (nvals == params->nce1 * params->nce2) {
	    double **d2 = d_alloc_2d(params->nce1, params->nce2);
	    n = 0;
	    for (j = 0; j < params->nce2; j++)
	      for (i = 0; i < params->nce1; i++)
		d2[j][i] = d1[n++];
	    c2s_2d(geom, ret, d2, geom->nce1, geom->nce2);
	    d_free_2d(d2);
	    d_free_1d(d1);
	  } else {
	    hd_quit("value_init : incorrect number specified for %s parameter.\n", keyword);
	  }
	}
      } else {
	for (cc = 1; cc < geom->sgnumS; cc++) {
	  ret[cc] = val;
	}
      }
    } else {
      /* Read from file                                            */
      timeseries_t *ts;
      int id;

      /* hd_ts_read() is used to read eta from file and this       */
      /* requires time variables to be set on the master. These    */
      /* are only set in compute_constants(), called after this    */
      /* routine, so preset them here.                             */
      master->t = params->t;
      strcpy(master->timeunit, params->timeunit);
      strcpy(master->output_tunit, params->output_tunit);

      ts = hd_ts_read(master, keyword, 0);
      id = ts_get_index(ts, fv_get_varname(keyword, vname, buf));
      if (id < 0)
	hd_quit("value_init: The file '%s' does not contain the tracer %s.\n",
		keyword, vname);

      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	ret[c] = ts_eval_xy(ts, id, master->t, geom->cellx[c], geom->celly[c]);
      }
    }
  }
  return(ret);
}

/* END value_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the geographic metrics on the sphere using a false pole   */
/*-------------------------------------------------------------------*/
void geog_false_pole_coord(double **x,  /* where to store grid x values */
                           double **y,  /* where to store grid y values */
                           double **h1, /* where to store h1 metric values
                                         */
                           double **h2, /* where to store h2 metric values
                                         */
                           double **a1, /* where to store a1 angle values */
                           double **a2, /* where to store a2 angle values */
                           long int nce1, /* number of cells in e1
                                             direction */
                           long int nce2, /* number of cells in e2
                                             direction */
                           double x00,  /* x origin offset */
                           double y00,  /* y origin offset */
                           double flon, /* False longitude */
                           double flat, /* False latitude */
                           double xinc, /* cell size in x direction */
                           double yinc  /* cell size in y direction */
  )
{
  long i, j;
  double fx00, fy00, xval, yval;
  double rflon = DEG2RAD(flon);
  double rflat = DEG2RAD(flat);

  geod_fwd_spherical_rot(DEG2RAD(x00), DEG2RAD(y00), rflon, rflat, &fx00,
                         &fy00);
  xinc = DEG2RAD(xinc);
  yinc = DEG2RAD(yinc);

  for (j = 0; j < nce2 + 1; j++) {
    yval = j * yinc;
    for (i = 0; i < nce1 + 1; i++) {
      xval = i * xinc;
      geod_inv_spherical_rot(fx00 + xval, fy00 + yval, rflon, rflat,
                             &x[j][i], &y[j][i]);
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }

  /* Calculate h1 and h2 numerically */
  grid_get_geog_metrics(x, y, nce1, nce2, h1, h2);

  /* calculate a1, a2 numerically */
  grid_get_geog_angle(x, y, nce1, nce2, a1, a2);

  /* Check ranges; must be 0 to 360 */
  // FR: 28-05-2018 Disabling to allow for models to run across 0 deg longitude
  // Note: mask_bathy_from_coast may need fixing as well
  // Reverted 14-06-2018
  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      while (x[j][i] < 0.0 || x[j][i] > 360.0) {
	if (x[j][i] < 0.0) x[j][i] += 360.0;
	if (x[j][i] > 360.0) x[j][i] -= 360.0;
      }
    }
  }
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Compute the geographic metrics on the sphere by computing the     */
/* latitude and longitude by dead reckoning. We use the Sodano's     */
/* direct formulation to compute the latitude/longitude given a      */
/* range and bearing.                                                */
/*-------------------------------------------------------------------*/
void geog_dreckon_coord(double **x, /* where to store grid x values */
                        double **y, /* where to store grid y values */
                        double **h1,  /* where to store h1 metric values */
                        double **h2,  /* where to store h2 metric values */
                        double **a1,  /* where to store a1 angle values */
                        double **a2,  /* where to store a2 angle values */
                        long int nce1,  /* number of cells in e1 direction
                                         */
                        long int nce2,  /* number of cells in e2 direction
                                         */
                        double x00, /* x origin offset */
                        double y00, /* y origin offset */
                        double rotn,  /* Rotation */
                        double xinc,  /* cell size in x direction */
                        double yinc /* cell size in y direction */
  )
{
  long i, j;
  double rx00 = DEG2RAD(x00);
  double ry00 = DEG2RAD(y00);
  double xaz = M_PI / 2 - DEG2RAD(rotn);
  double yaz = xaz - M_PI / 2;

#if 0
  for (j = 0; j < nce2 + 1; j++) {
    double slat;                /* Start latitude of line */
    double slon;                /* Start longitude of line */
    geod_fwd_sodanos(rx00, ry00, yaz, j * yinc, RADIUS, ECC, &slon, &slat);
    for (i = 0; i < nce1 + 1; i++) {
      geod_fwd_sodanos(slon, slat, xaz, i * xinc, RADIUS, ECC, &x[j][i],
                       &y[j][i]);
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }
#else
  x[0][0] = rx00;
  y[0][0] = ry00;
  for (j = 0; j < nce2 + 1; j++) {
    if (j > 0)
      geod_fwd_sodanos(x[j - 1][0], y[j - 1][0], yaz, yinc, RADIUS, ECC,
                       &x[j][0], &y[j][0]);
    for (i = 1; i < nce1 + 1; i++)
      geod_fwd_sodanos(x[j][i - 1], y[j][i - 1], xaz, xinc, RADIUS, ECC,
                       &x[j][i], &y[j][i]);
  }

  for (j = 0; j < nce2 + 1; j++) {
    for (i = 0; i < nce1 + 1; i++) {
      x[j][i] = RAD2DEG(x[j][i]);
      y[j][i] = RAD2DEG(y[j][i]);
    }
  }
#endif

  /* Calculate h1 and h2 numerically */
  grid_get_geog_metrics(x, y, nce1, nce2, h1, h2);

  /* calculate a1, a2 numerically */
  grid_get_geog_angle(x, y, nce1, nce2, a1, a2);
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
int read2darray(FILE * fp, char *label, double **array, long n1, long n2)
{
  double *vals;
  int numvals;
  long n;
  long i;
  long j;

  vals = NULL;
  if (!prm_read_darray(fp, label, &vals, &numvals))
    return (0);
  if (numvals != n1 * n2)
    hd_quit("read2darray: incorrect number of values\n");
  n = 0;
  for (j = 0; j < n2; j++)
    for (i = 0; i < n1; i++)
      array[j][i] = vals[n++];
  free(vals);
  return 1;
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void params_write(parameters_t *params, dump_data_t *dumpdata)
{

  FILE *op, *fopen();
  int n, tn, i, j, ntr;
  double d1;
  char tag[MAXSTRLEN];
  char key[MAXSTRLEN];
  char bname[MAXSTRLEN];

  sprintf(tag, "%s.prm", params->oname);
  op = fopen(tag, "w");
  fprintf(op, "# SHOC parameter file\n");
  fprintf(op, "CODEHEADER           %s\n", params->codeheader);
  fprintf(op, "PARAMETERHEADER      %s\n", params->parameterheader);
  fprintf(op, "DESCRIPTION          %s\n", params->grid_desc);
  fprintf(op, "NAME                 %s\n", params->grid_name);
  if (strlen(params->projection))
    fprintf(op, "PROJECTION           %s\n", params->projection);
  fprintf(op, "TIMEUNIT             %s\n", params->timeunit);
  fprintf(op, "OUTPUT_TIMEUNIT      %s\n", params->output_tunit);
  fprintf(op, "LENUNIT              %s\n", params->lenunit);
  fprintf(op, "START_TIME           %s\n", params->start_time);
  fprintf(op, "STOP_TIME            %s\n\n", params->stop_time);
  fprintf(op, "INPUT_FILE           %s.nc\n\n", params->oname);

  if (params->trout) {
    if (strlen(params->trans_dt))
      fprintf(op,"TRANS_DT             %s\n", params->trans_dt);
    fprintf(op,"TRANS_OUTPUT         YES\n");
    if (params->tmode == SP_FFSL)
      fprintf(op,"TRANS_MODE           SP_FFSL\n");
    if (strlen(params->trkey))
      fprintf(op,"OutputTransport        %s\n",params->trkey);
    fprintf(op,"\n");
  }
    

  if(params->ndf) {
    fprintf(op,"# Output files\n");
    if (strlen(params->opath))
      fprintf(op,"OutputPath             %s\n\n",params->opath);
    fprintf(op,"OutputFiles            %d\n\n",params->ndf);
    for(n = 0; n < params->ndf; n++) {
      fprintf(op,"file%1.1d.name           %s\n",n,params->d_name[n]);
      fprintf(op,"file%1.1d.filetype       %s\n",n,params->d_filetype[n]);
      if (strlen(params->d_tstart[n]))
	fprintf(op,"file%1.1d.tstart         %s\n",n,params->d_tstart[n]);
      fprintf(op,"file%1.1d.tinc           %s\n",n,params->d_tinc[n]);
      if (strlen(params->d_tstop[n]))
	fprintf(op,"file%1.1d.tstop          %s\n",n,params->d_tstop[n]);
      fprintf(op,"file%1.1d.bytespervalue  4\n",n);
      fprintf(op,"file%1.1d.vars           %s\n",n,params->d_vars[n]);
      if (strlen(params->d_tunit[n]))
	fprintf(op,"file%1.1d.timeunit       %s\n",n,params->d_tunit[n]);
      fprintf(op,"\n");
    }
  }
  else {
    fprintf(op, "# Output files\n");
    fprintf(op, "OutputFiles 1\n\n");
    fprintf(op, "file0.name           out1\n");
    fprintf(op, "file0.filetype       standard\n");
    fprintf(op, "file0.tstart         %s\n", params->start_time);
    fprintf(op, "file0.tinc           5 days\n");
    fprintf(op, "file0.tstop          %s\n", params->stop_time);
    fprintf(op, "file0.bytespervalue  4\n");
    fprintf(op, "file0.vars           ALL\n\n");
  }

  fprintf(op, "# Flags\n");
  fprintf(op, "WINDOWS              %d\n", params->nwindows);
  if (params->nwindows > 1) {
    /* Default form RECOM */
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
      fprintf(op, "DP_MODE              OPENMP\n");
    else
      fprintf(op, "DP_MODE              %s\n", params->dp_mode);
  }
  fprintf(op, "NONLINEAR            %s\n", tf(params->nonlinear));
  fprintf(op, "CALCDENS             %s\n", tf(params->calc_dens));
  fprintf(op, "HEATFLUX             %s\n", heatfluxname(params->heatflux));
  fprintf(op, "SALTFLUX             %s\n", heatfluxname(params->saltflux));
  fprintf(op, "2D-MODE              %s\n", tf(params->mode2d));
  fprintf(op, "STABILITY            %s\n", substepname(params->stab));
  fprintf(op, "RAMPSTART            %d days\n", (int)params->t / 86400);
  fprintf(op, "RAMPEND              %d days\n",
          (int)params->t / 86400 + 1);
  if (params->rampf != 0) {
    fprintf(op, "RAMPVARS            ");
    if (params->rampf & WIND) fprintf(op, " WIND");
    if (params->rampf & TIDALH) fprintf(op, " TIDALH");
    if (params->rampf & TIDALC) fprintf(op, " TIDALC");
    if (params->rampf & FILEIN) fprintf(op, " FILEIN");
    if (params->rampf & CUSTOM) fprintf(op, " CUSTOM");
    if (params->rampf & FLATHR) fprintf(op, " FLATHR");
    if (params->rampf & ETA_RELAX) fprintf(op, " ETA_RELAX");
    if (params->rampf & INV_BARO) fprintf(op, " INV_BARO");
    if (params->rampf & FLUX_ADJUST) fprintf(op, " FLUX_ADJUST");
    if (params->rampf & STOKES) fprintf(op, " STOKES");
    fprintf(op, "\n");
  }

  fprintf(op, "MERGE_THIN           %s\n", tf(params->thin_merge));
  fprintf(op, "HMIN                 %-6.4f\n", params->hmin);
  fprintf(op, "SLIP                 %-6.1f\n", params->slipprm);
  fprintf(op, "SIGMA                %s\n", tf(params->sigma));
  fprintf(op, "\n");

  fprintf(op, "# Time steps\n");
  fprintf(op, "DT                   %-8.2f seconds\n", params->grid_dt);
  fprintf(op, "IRATIO               %d\n", params->iratio);
  fprintf(op, "TRATIO               %5.1f\n\n", params->tratio);

  fprintf(op, "# Advection\n");
  fprintf(op, "MOM_SCHEME           %s\n", adname(params->momsc));
  fprintf(op, "TRA_SCHEME           %s\n", adname(params->trasc));
  fprintf(op, "ULTIMATE             %s\n\n", tf(params->ultimate));

  fprintf(op, "# Horizontal mixing\n");
  if (params->u1vh > 1.0)
    fprintf(op, "U1VH                 %-6.1f\n", params->u1vh);
  else
    fprintf(op, "U1VH                 %-8.3f\n", params->u1vh);
  if (params->u2vh > 1.0)
    fprintf(op, "U2VH                 %-6.1f\n", params->u2vh);
  else
    fprintf(op, "U2VH                 %-8.3f\n", params->u2vh);
  if (params->u1kh > 1.0)
    fprintf(op, "U1KH                 %-6.1f\n", params->u1kh);
  else
    fprintf(op, "U1KH                 %-8.3f\n", params->u1kh);
  if (params->u2kh > 1.0)
    fprintf(op, "U2KH                 %-6.1f\n", params->u2kh);
  else
    fprintf(op, "U2KH                 %-8.3f\n", params->u2kh);
  /*fprintf(op, "SMAGORINSKY          %-6.4f\n", params->smagorinsky);*/
  fprintf(op, "SMAGORINSKY          %s\n", params->smag);
  if (params->smag_smooth)
    fprintf(op, "SMAG_SMOOTH          %d\n", params->smag_smooth);
  if (params->visc_method == NONE)
    fprintf(op, "VISC_METHOD          NONE\n");
  else if (params->visc_method == LAPLACIAN)
    fprintf(op, "VISC_METHOD          LAPLACIAN\n");
  else if (params->visc_method == SIMPLE)
    fprintf(op, "VISC_METHOD          SIMPLE\n");
  else if (params->visc_method == (LAPLACIAN|PRE794))
    fprintf(op, "VISC_METHOD          PRE_V794\n");
  if (params->robust > 1)
    fprintf(op, "ROBUST               %d\n", params->robust);
  /*
  if (params->speed > 1)
    fprintf(op, "SPEED                %d\n", params->speed);
  */
  fprintf(op, "\n");

  fprintf(op, "# Vertical mixing\n");
  fprintf(op, "MIXING_SCHEME        %s\n", params->mixsc);
  fprintf(op, "VZ0                  %-6.4e\n", params->vz0);
  fprintf(op, "KZ0                  %-6.4e\n", params->kz0);
  if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
      strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
    fprintf(op, "ZS                   %-6.4f\n\n", params->zs);
  else {
    fprintf(op, "MIN_TKE              %-6.4e\n\n", params->min_tke);
    fprintf(op, "ZS                   %-6.4f\n\n", params->zs);
  }

  fprintf(op, "# Bottom friction\n");
  fprintf(op, "QBFC                 %-6.4f\n", params->quad_bfc);
  fprintf(op, "UF                   %-6.5f\n", params->uf);
  fprintf(op, "Z0                   %-6.4f\n\n", params->z0);
  if (params->z0s)
    fprintf(op, "%-6.4e\n", params->z0s[0]);

  fprintf(op, "# Constants\n");
  fprintf(op, "G                    %-6.4f\n", params->g);
  fprintf(op, "SPECHEAT             %-6.1f\n", params->spec_heat);
  fprintf(op, "AIRDENS              %-6.4f\n", params->air_dens);
  fprintf(op, "AMBIENT_AIR_PRESSURE %-6.4f\n", params->ambpress);
  if (params->coriolis) {
    fprintf(op, "CORIOLIS             %d\n", params->nce1 * params->nce2);
    fprintf(op, "%-6.4e\n", params->coriolis[0]);
  }
  fprintf(op, "\n");

  fprintf(op, "# Diagnostics\n");
  fprintf(op, "CFL                  %s\n", cflname(params->cfl));
  if (params->cfl & (ACTIVE | ACTIVE3D))
    fprintf(op, "CFL_DT               %s\n", params->cfl_dt);
  fprintf(op, "MIX_LAYER            %s\n", mixname(params->mixlayer));
  fprintf(op, "MEAN                 %s\n", meansname(params->means, key));
  if (strlen(params->means_dt))
    fprintf(op, "MEAN_DT              %s\n", params->means_dt);
  if (strlen(params->means_os))
    fprintf(op, "MEAN_OFFSET          %s\n", params->means_os);
  fprintf(op, "VORTICITY            NONE\n");
  fprintf(op, "NUMBERS              NONE\n");
  if (strlen(params->alert)) {
    fprintf(op, "ALERT                %s\n", params->alert);
    if (contains_token(params->alert, "ACTIVE") != NULL) {
      fprintf(op, "ALERT_CODE           %d %d %d %d %d %d %d %d %d %d %d\n",
	      params->eta_f, params->vel2d_f, params->vel3d_f,
	      params->wvel_f, params->tend_f, params->div2d_f,
	      params->div3d_f, params->cfl_f, params->ts_f,
	      params->shear_f, params->hdiff_f);
    }
    if (params->eta_f && params->etarlx & ALERT) {
      fprintf(op, "eta_relaxation_file  %s\n", params->etarlxn);
      if (params->compatible & V1670)
	fprintf(op, "COMPATIBLE           V1670\n");
    }
    if (strlen(params->alert_dt))
      fprintf(op, "ALERT_DT             %s\n", params->alert_dt);
  }
  else {
    fprintf(op, "ALERT                NONE\n");
  }
  fprintf(op, "MOM_TEND             %s\n", tf(params->tendf));
  if (strlen(params->trtend))
    fprintf(op, "TRA_TEND             %s\n", params->trtend);
  if (strlen(params->trflux))
    fprintf(op, "CALC_FLUXES          %s\n", params->trflux);
  if (strlen(params->trperc))
    fprintf(op, "CALC_PERCS           %s\n", params->trperc);
  fprintf(op, "FLUSHING_TR          %s\n", tf(params->trflsh));
  if (strlen(params->trage))
    fprintf(op, "AGE_TR               %s\n", params->trage);
  fprintf(op, "TOTALS               %s\n", tf(params->totals));
  if (params->avhrr)
    fprintf(op, "AVHRR                %s\n", params->avhrr_path);
  if (params->ghrsst)
    fprintf(op, "GHRSST               %s\n", params->ghrsst_path);
  fprintf(op, "STERIC_HEIGHT        %-6.2f\n", params->lnm);

  fprintf(op, "\n");
  
  write_grid_specs(op, params);
  
  fprintf(op, "# Vertical grid spacing\n");
  fprintf(op, "LAYERFACES           %d\n", params->nz + 1);
  d1 = 1.0;
  if (params->sigma)
    d1 = params->bmax;
  for (n = 0; n < params->nz; n++)
    fprintf(op, "%-8.2f\n", params->layers[n] * d1);
  fprintf(op, "%-8.2f\n\n", params->layers[params->nz]);

  fprintf(op, "# Bathymetry limits\n");
  fprintf(op, "BATHYMIN             %-6.1f\n", params->bmin);
  fprintf(op, "BATHYMAX             %-6.1f\n", params->bmax);
  fprintf(op, "ETAMAX               %-6.1f\n", params->etamax);
  if (params->etadiff)
    fprintf(op, "ETA_DIFF             %-6.1f\n", params->etadiff);
  fprintf(op, "MIN_CELL_THICKNESS   %s\n", params->mct);
  if (params->smooth)
    fprintf(op, "SMOOTHING            %d\n", params->smooth);
  if (params->maxgrad > 0.0)
    fprintf(op, "MAXGRAD              %-6.2f\n", params->maxgrad);

  fprintf(op, "\n# Tracers\n");
  ntr = params->ntr;
  if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
    ntr = 2; if (strlen(params->dhw)) ntr += 3;
  }
  fprintf(op, "NTRACERS             %d\n\n", ntr);
  tn = 0;
  for (n = 0; n < params->ntr; n++) {
    tracer_info_t *tracer = &params->trinfo_3d[n];

    /* Special handling for Degree Heating weeks in RECOM */
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2) && n > 1)
      if (strncmp(tracer->name, "dhw", 3) != 0 &&
	  strcmp(tracer->name, "salt") != 0    &&
	  strcmp(tracer->name, "temp") != 0)
	continue;

    fprintf(op, "TRACER%1.1d.name         %s\n", tn, tracer->name);
    fprintf(op, "TRACER%1.1d.long_name    %s\n", tn, tracer->long_name);
    fprintf(op, "TRACER%1.1d.units        %s\n", tn, tracer->units);
    if (tracer->fill_value_wc < 1)
      fprintf(op, "TRACER%1.1d.fill_value   %-4.2e\n", tn,
              tracer->fill_value_wc);
    else
      fprintf(op, "TRACER%1.1d.fill_value   %-6.1f\n", tn,
              tracer->fill_value_wc);
    if (n > 1)
      fprintf(op, "TRACER%1.1d.valid_range  %-4.2e %-4.2e\n", tn,
	      tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
    else
      fprintf(op, "TRACER%1.1d.valid_range  %-6.1f %-6.1f\n", tn,
	      tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
    fprintf(op, "TRACER%1.1d.advect       %d\n", tn, tracer->advect);
    fprintf(op, "TRACER%1.1d.diffuse      %d\n", tn, tracer->diffuse);
    fprintf(op, "TRACER%1.1d.diagn        %d\n", tn, tracer->diagn);
    if(!(params->roammode & (A_RECOM_R1|A_RECOM_R2)) && (tracer->type & (WATER|SEDIM) || tracer->type & INTER))
      fprintf(op, "TRACER%1.1d.type         %s\n", tn, trtypename(tracer->type, key));
    if (strlen(tracer->tracerstat)) {
      fprintf(op, "TRACER%1.1d.tracerstat   %s\n", tn, tracer->tracerstat);
      /* Special handling for scale factor for DHW in RECOM */
      if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
	if (strcmp(tracer->name, "dhw") == 0) {
	  fprintf(op, "TRACER%1.1d.dt           7 days\n", tn);
	  fprintf(op, "TRACER%1.1d.scale_factor 7 days\n", tn);
	}
    }
    if (strlen(params->trinfo_3d[n].data))
      fprintf(op, "TRACER%1.1d.data         %s\n", tn, params->trinfo_3d[n].data);
    if (strlen(params->trrlxn[n])) {
      fprintf(op, "TRACER%1.1d.relaxation_file     %s\n",
              tn, params->trrlxn[n]);
      fprintf(op, "TRACER%1.1d.relaxation_input_dt %s\n",
              tn, otime(params->trrlxdt[n], tag));
      fprintf(op, "TRACER%1.1d.relaxation_time_constant  %s\n",
              tn, params->trrlxtc[n]);
    }
    if (strlen(params->trrest[n])) {
      fprintf(op, "TRACER%1.1d.reset_file   %s\n",
              tn, params->trrest[n]);
      if (params->save_force && (strcmp(tracer->name, "otemp") == 0 ||
				 strcmp(tracer->name, "osalt") == 0 ||
				 strcmp(tracer->name, "ovelu") == 0 ||
				 strcmp(tracer->name, "ovelv") == 0))
	fprintf(op, "TRACER%1.1d.reset_dt     1 day\n", tn);
	/*
	fprintf(op, "TRACER%1.1d.reset_dt     %-8.2f seconds\n",
		tn, params->grid_dt);
	*/
      else
	fprintf(op, "TRACER%1.1d.reset_dt     %s\n",
		tn, otime(params->trrestdt[n], tag));
    }
    tn++;
    fprintf(op, "\n");
  }

  fprintf(op, "# Forcing\n");
  if (strlen(params->wind)) {
    fprintf(op, "WIND_TS              %s\n", params->wind);
    fprintf(op, "WIND_INPUT_DT        %s\n", otime(params->wind_dt, tag));
    fprintf(op, "WIND_SPEED_SCALE     %-6.1f\n", params->wind_scale);
    if (params->wind_type == SPEED) {
      fprintf(op, "DRAG_LAW_V0          %-6.1f\n", params->dlv0);
      fprintf(op, "DRAG_LAW_V1          %-6.1f\n", params->dlv1);
      fprintf(op, "DRAG_LAW_CD0         %-8.5f\n", params->dlc0);
      fprintf(op, "DRAG_LAW_CD1         %-8.5f\n", params->dlc1);
    } else
      fprintf(op, "WIND_TYPE            STRESS\n");
    fprintf(op, "\n");
  }
  if (strlen(params->patm)) {
    fprintf(op, "PRESSURE             %s\n", params->patm);
    fprintf(op, "PRESSURE_INPUT_DT    %s\n", otime(params->patm_dt, tag));
    fprintf(op, "\n");
  }
  if (strlen(params->eta_init))
    fprintf(op, "SURFACE              %s\n\n", params->eta_init);
  if (params->etarlx & RELAX) {
    fprintf(op, "eta_relaxation_file           %s\n", params->etarlxn);
    fprintf(op, "eta_relaxation_input_dt       %s\n",
	    otime(params->etarlxdt, tag));
    fprintf(op, "eta_relaxation_time_constant  %s\n\n",
	    otime(params->etarlxtc, tag));
  }
  if (strlen(params->vel_init))
    fprintf(op, "VELOCITY             %s\n\n", params->vel_init);
  if (params->heatflux & COMP_HEAT) {
    fprintf(op, "HEATFLUX             COMP_HEAT\n");
    fprintf(op, "HEATFLUX_FILE	      %s\n", params->hf);
    fprintf(op, "HEATFLUX_DT          %s\n", otime(params->hf_dt, tag));
    fprintf(op, "RADIATION            %s\n", params->swr);
    fprintf(op, "RADIATION_INPUT_DT   %s\n", otime(params->swr_dt, tag));
    fprintf(op, "SWR_ATTENUATION      %s\n", params->swr_attn);
    fprintf(op, "SWR_TRANSMISSION     %s\n", params->swr_tran);
    fprintf(op, "SWR_BOT_ABSORB       %s\n", params->swr_babs);
    fprintf(op, "ALBEDO               %-6.2f\n", params->albedo);
    fprintf(op, "\n");
  }
  if (params->heatflux & ADVANCED) {
    fprintf(op, "HEATFLUX             BULK\n");
    if (master->airtemp) {
      fprintf(op, "AIRTEMP              %s\n", params->airtemp);
      fprintf(op, "AIRTEMP_INPUT_DT     %s\n", otime(params->airtemp_dt, tag));
    }
    if (master->wetb && master->sh_f == WETBULB) {
      fprintf(op, "WET_BULB             %s\n", params->wetb);
      fprintf(op, "WET_BULB_INPUT_DT    %s\n", otime(params->wetb_dt, tag));
    }
    if (master->wetb && master->sh_f == DEWPOINT) {
      fprintf(op, "DEW_POINT            %s\n", params->wetb);
      fprintf(op, "DEW_POINT_INPUT_DT   %s\n", otime(params->wetb_dt, tag));
    }
    if (master->rh) {
      fprintf(op, "HUMIDITY             %s\n", params->rh);
      fprintf(op, "HUMIDITY_INPUT_DT    %s\n", otime(params->rh_dt, tag));
    }
    if (master->cloud) {
      fprintf(op, "CLOUD                %s\n", params->cloud);
      fprintf(op, "CLOUD_INPUT_DT       %s\n", otime(params->cloud_dt, tag));
    }
    fprintf(op, "SWR_ATTENUATION      %s\n", params->swr_attn);
    fprintf(op, "SWR_TRANSMISSION     %s\n", params->swr_tran);
    fprintf(op, "SWR_BOT_ABSORB       %s\n", params->swr_babs);
    fprintf(op, "\n");
  }
  if (params->saltflux & BULK) {
    fprintf(op, "SALTFLUX             BULK\n");
    if (master->cloud) {
      fprintf(op, "PRECIPITATION        %s\n", params->precip);
      fprintf(op, "PRECIPITATION_INPUT_DT %s\n", otime(params->precip_dt, tag));
    }
    fprintf(op, "\n");
  }
  if (params->heatflux & COMP_HEAT_MOM) {
    fprintf(op, "HEATFLUX             COMP_HEAT_MOM\n");
    fprintf(op, "HEATFLUX_FILE	      %s\n", params->hf);
    fprintf(op, "HEATFLUX_DT          %s\n", otime(params->hf_dt, tag));
    fprintf(op, "RADIATION            %s\n", params->swr);
    fprintf(op, "RADIATION_INPUT_DT   %s\n", otime(params->swr_dt, tag));
    fprintf(op, "SWR_ATTENUATION      %s\n", params->swr_attn);
    fprintf(op, "SWR_TRANSMISSION     %s\n", params->swr_tran);
    fprintf(op, "SWR_BOT_ABSORB       %s\n", params->swr_babs);
    fprintf(op, "ALBEDO               %-6.2f\n", params->albedo_l);
    fprintf(op, "\n");
  }
  fprintf(op, "# Time series\n");
  prm_set_errfn(hd_silent_warn);
  if (prm_read_int(params->prmfd, "TSPOINTS", &tn)) {
    fprintf(op, "TSPOINTS             %d\n\n", tn);
    for (n = 0; n < tn; n++) {
      sprintf(key, "TS%1d.name", n);
      prm_read_char(params->prmfd, key, tag);
      fprintf(op, "TS%1d.name              %s\n", n, tag);
      sprintf(key, "TS%1d.location", n);
      prm_read_char(params->prmfd, key, tag);
      fprintf(op, "TS%1d.location          %s\n", n, tag);
      sprintf(key, "TS%1d.dt", n);
      prm_read_char(params->prmfd, key, tag);
      fprintf(op, "TS%1d.dt                %s\n", n, tag);
      sprintf(key, "TS%1d.reference", n);
      if (prm_read_char(params->prmfd, key, tag))
        fprintf(op, "TS%1d.reference         %s\n", n, tag);
      else
        fprintf(op, "TS%1d.reference         msl\n", n);
      fprintf(op, "\n");
    }
  } else {
    fprintf(op, "TSPOINTS             1\n\n");
    fprintf(op, "TS0.name             loc1.ts\n");
    fprintf(op, "TS0.location         %-10.4f %-10.4f 0\n",
            dumpdata->cellx[params->nce2 / 2][params->nce1 / 2],
            dumpdata->celly[params->nce2 / 2][params->nce1 / 2]);
    fprintf(op, "TS0.dt               1 hour\n");
    fprintf(op, "TS0.reference        msl\n\n");
  }

  fprintf(op, "# Open boundaries\n");
  fprintf(op, "NBOUNDARIES           %d\n\n", params->nobc);
  if (strlen(params->orthoweights))
    fprintf(op, "TIDE_CSR_CON_DIR %s\n", params->nodal_dir);
  if (strlen(params->nodal_dir))
    fprintf(op, "TIDE_CSR_ORTHOWEIGHTS %s\n\n", params->orthoweights);
  if (strlen(params->tide_con_file))
    fprintf(op, "TIDE_CONSTITUENTS %s\n\n", params->tide_con_file);
  if (params->tide_r == CSR_R)
    fprintf(op, "TIDAL_REMOVAL       CSR\n");
  if (params->tide_r == MEAN_R)
    fprintf(op, "TIDAL_REMOVAL       MEAN\n");
  for (n = 0; n < params->nobc; n++) {
    open_bdrys_t *open = params->open[n];
    fprintf(op, "BOUNDARY%1.1d.NAME          %s\n", n, open->name);
    fprintf(op, "BOUNDARY%1.1d.TYPE          %s\n", n, btname(open->type));
    bcname(open->bcond_nor, bname);
    fprintf(op, "BOUNDARY%1.1d.BCOND_NOR     %s\n", n, bname);
    if (open->bcond_nor2d != open->bcond_nor) {
      bcname(open->bcond_nor2d, bname);
      fprintf(op, "BOUNDARY%1.1d.BCOND_NOR2D   %s\n", n, bname);
    }
    bcname(open->bcond_ele, bname);
    fprintf(op, "BOUNDARY%1.1d.BCOND_ELE     %s\n", n, bname);
    bcname(open->bcond_tan, bname);
    fprintf(op, "BOUNDARY%1.1d.BCOND_TAN     %s\n", n, bname);
    if (open->bcond_tan2d != open->bcond_tan) {
      bcname(open->bcond_tan2d, bname);
      fprintf(op, "BOUNDARY%1.1d.BCOND_TAN2D   %s\n", n, bname);
    }
    fprintf(op, "BOUNDARY%1.1d.BCOND_TRA_ALL NOGRAD\n", n);
    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      bcname(open->bcond_tra[tn], bname);
      if (open->bcond_tra[tn] != NOGRAD)
	fprintf(op, "BOUNDARY%1.1d.BCOND_TRA%1.1d    %s\n", n, tn, bname);
      if (open->bcond_tra[tn] & CUSTOM)
	fprintf(op, "BOUNDARY%1.1d.CUSTOM.%s   %s\n", n, tracer->name, open->cusname_t[tn]);
    }
    bcname(open->bcond_Vz, bname);
    fprintf(op, "BOUNDARY%1.1d.BCOND_VZ      %s\n", n, bname);
    bcname(open->bcond_Kz, bname);
    fprintf(op, "BOUNDARY%1.1d.BCOND_KZ      %s\n", n, bname);
    fprintf(op, "BOUNDARY%1.1d.INVERSE_BAROMETER %s\n", n,
            tf(open->inverse_barometer));
    if (params->roammode & (A_RECOM_R1|A_RECOM_R2)) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      if (strlen(open->tsfn)) {
	strcpy(tag, open->tsfn);
	i = parseline(tag, fields, MAXNUMARGS);
	fprintf(op, "BOUNDARY%1.1d.DATA          %s\n", n, fields[0]);
      }
    } else {
      if (strlen(open->tsfn))
	fprintf(op, "BOUNDARY%1.1d.DATA          %s\n", n, open->tsfn);
    }
    if (strlen(open->cusname_u1))
      fprintf(op, "BOUNDARY%1.1d.CUSTOM.u1     %s\n", n, open->cusname_u1);
    if (strlen(open->cusname_u2))
      fprintf(op, "BOUNDARY%1.1d.CUSTOM.u2     %s\n", n, open->cusname_u2);
    if (strlen(open->cusname_u1av))
      fprintf(op, "BOUNDARY%1.1d.CUSTOM.u1av   %s\n", n, open->cusname_u1av);
    if (strlen(open->cusname_u2av))
      fprintf(op, "BOUNDARY%1.1d.CUSTOM.u2av   %s\n", n, open->cusname_u2av);
    if (strlen(open->bflow)) {
      if (open->type & U1BDRY) {
        fprintf(op, "BOUNDARY%1.1d.U1_FLOW       %s\n", n, open->bflow);
        fprintf(op, "BOUNDARY%1.1d.U1_HC         %-6.1f\n", n,
                open->bhc);
      }
      if (open->type & U2BDRY) {
        fprintf(op, "BOUNDARY%1.1d.U2_FLOW       %s\n", n, open->bflow);
        fprintf(op, "BOUNDARY%1.1d.U2_HC         %-6.1f\n", n,
                open->bhc);
      }
    }
    if (open->sponge_zone_h)
      fprintf(op, "BOUNDARY%1.1d.NSPONGE_HORZ  %d\n", n,
              open->sponge_zone_h);
    if (open->sponge_f)
      fprintf(op, "BOUNDARY%1.1d.SPONGE_FACT  %-6.1f\n", n,
              open->sponge_f);
    if (open->sponge_zone)
      fprintf(op, "BOUNDARY%1.1d.NSPONGE_VERT  %d\n", n, open->sponge_zone);
    if (open->relax_timei) {
      fprintf(op, "BOUNDARY%1.1d.RELAX_OUT     %s\n", n,
	      otime(open->relax_time, tag));
      fprintf(op, "BOUNDARY%1.1d.RELAX_IN      %s\n", n,
	      otime(open->relax_timei, tag));
    }
    else if (open->relax_time)
      fprintf(op, "BOUNDARY%1.1d.RELAX_TIME    %s\n", n,
	      otime(open->relax_time, tag));
    if (open->rele_b && open->rele_i)
      fprintf(op, "BOUNDARY%1.1d.RELAX_ELE     %d %d %d\n", n,
	      open->relax_ele, (int)open->rele_b, (int)open->rele_i);
    if (open->adjust_flux)
      fprintf(op, "BOUNDARY%1.1d.ADJUST_FLUX   %s\n", n,
	      otime(open->adjust_flux, tag));      
    if (open->relax_zone_nor)
      fprintf(op, "BOUNDARY%1.1d.RELAX_ZONE_NOR %d\n", n, open->relax_zone_nor);
    if (open->relax_zone_tan)
      fprintf(op, "BOUNDARY%1.1d.RELAX_ZONE_TAN %d\n", n, open->relax_zone_tan);
    if (open->stagger & INFACE)
      fprintf(op, "BOUNDARY%1.1d.STAGGER       INFACE\n", n);
    if (open->bathycon)
      fprintf(op, "BOUNDARY%1.1d.BATHY_CON     %d\n", n, open->bathycon);
    fprintf(op, "BOUNDARY%1.1d.POINTS        %d\n", n, open->npts);
    for (tn = 0; tn < open->npts; tn++)
      fprintf(op, "%d %d\n", open->iloc[tn], open->jloc[tn]);
    fprintf(op, "\n");
  }

  fprintf(op, "# Bathymetry\n");
  fprintf(op, "BATHY    %d\n", params->nce1 * params->nce2);
  for (j = 0; j < params->nce2; j++)
    for (i = 0; i < params->nce1; i++)
      if (isnan(dumpdata->botz[j][i]))
        fprintf(op, "%8.3f\n", -(double)NOTVALID);
      else
        fprintf(op, "%8.3f\n", -dumpdata->botz[j][i]);

  fprintf(op, "\n");
  if (strcmp(params->gridtype, "NUMERICAL") == 0) {
    fprintf(op, "XCOORDS              %d\n", (2 * params->nce1 + 1) *
            (2 * params->nce2 + 1));
    for (j = 0; j < 2 * params->nce2 + 1; j++)
      for (i = 0; i < 2 * params->nce1 + 1; i++)
        fprintf(op, "%f\n", params->x[j][i]);
    fprintf(op, "\n");
    fprintf(op, "YCOORDS              %d\n", (2 * params->nce1 + 1) *
            (2 * params->nce2 + 1));
    for (j = 0; j < 2 * params->nce2 + 1; j++)
      for (i = 0; i < 2 * params->nce1 + 1; i++)
        fprintf(op, "%f\n", params->y[j][i]);
  }

  fclose(op);
  if (params->coriolis)
    d_free_1d(params->coriolis);
  if (params->surface)
    d_free_1d(params->surface);
  d_free_2d(params->x);
  d_free_2d(params->y);
  d_free_2d(params->h1);
  d_free_2d(params->h2);
  d_free_2d(params->a1);
  d_free_2d(params->a2);
}

/* END params_write()                                                */
/*-------------------------------------------------------------------*/

/* Common to both param/trans_write */
void write_grid_specs(FILE *op, parameters_t *params)
{
  fprintf(op, "# Grid\n");
  if (strlen(params->roms_grid))
    fprintf(op, "ROMS_GRID            %s\n", params->roms_grid);
  else
    fprintf(op, "GRIDTYPE             %s\n", params->gridtype);
  fprintf(op, "NCE1                 %d\n", params->nce1);
  fprintf(op, "NCE2                 %d\n", params->nce2);
  /*fprintf(op, "PERIODIC             NONE\n");*/
  if (params->gridcode & (GRECT_DXY_ROT|GRECT_DLATLON_ROT|GRECT_DLATLON_FP)) {
    fprintf(op, "X00                  %-8.5f\n", params->rg->x0);
    fprintf(op, "Y00                  %-8.5f\n", params->rg->y0);
    if (params->gridcode & (GRECT_DLATLON_ROT|GRECT_DLATLON_FP)) {
      fprintf(op, "DLAMBDA              %-8.5f\n", params->rg->dx);
      fprintf(op, "DPHI                 %-8.5f\n", params->rg->dy);
    } else {
      fprintf(op, "DX                   %-8.5f\n", params->rg->dx);
      fprintf(op, "DY                   %-8.5f\n", params->rg->dy);
    }
    if (params->gridcode & GRECT_DLATLON_FP) {
      fprintf(op, "POLE_LATITUDE        %-8.5f\n", params->flat);
      fprintf(op, "POLE_LONGITUDE       %-8.5f\n", params->flon);
    } else {
      fprintf(op, "ROTATION             %-8.5f\n", params->rg->th);
    }
  } else if (params->gridcode & RECT_DXY_ROT) {
    fprintf(op, "X00                  %-8.5f\n", params->rg->x0);
    fprintf(op, "Y00                  %-8.5f\n", params->rg->y0);
    fprintf(op, "DX                   %-8.5f\n", params->rg->dx);
    fprintf(op, "DY                   %-8.5f\n", params->rg->dy);
    fprintf(op, "ROTATION             %-8.5f\n", params->rg->th);
  } else if (params->gridcode & POLAR_DARC_ROT) {
    fprintf(op, "X00                  %-8.5f\n", params->pg->x0);
    fprintf(op, "Y00                  %-8.5f\n", params->pg->y0);
    fprintf(op, "ARC                  %-8.5f\n", params->pg->arc);
    fprintf(op, "R0                   %-8.5f\n", params->pg->rmin);
    fprintf(op, "ROTATION             %-8.5f\n", params->pg->rotation);
  }
  fprintf(op, "\n");
}


/*-------------------------------------------------------------------*/
/* calculate cell and face neighbours for non-periodic grid          */
/*-------------------------------------------------------------------*/
void neighbour_none(dump_data_t *dumpdata)
{
  long i, j;

  /* cells and faces left and right from cell i */
  for (i = 0; i < dumpdata->nce1; i++) {
    dumpdata->crci[i] =
      (short)((i < dumpdata->nce1 - 1) ? i + 1 : NOTVALID);
    dumpdata->clci[i] = (short)((i > 0) ? i - 1 : NOTVALID);
    dumpdata->frci[i] = (short)(i + 1);
    dumpdata->flci[i] = (short)(i);
  }
  /* cells and faces left and right from face i */
  for (i = 0; i < dumpdata->nfe1; i++) {
    dumpdata->crfi[i] = (short)((i < dumpdata->nfe1 - 1) ? i : NOTVALID);
    dumpdata->clfi[i] = (short)((i > 0) ? i - 1 : NOTVALID);
    dumpdata->frfi[i] =
      (short)((i < dumpdata->nfe1 - 1) ? i + 1 : NOTVALID);
    dumpdata->flfi[i] = (short)((i > 0) ? i - 1 : NOTVALID);
  }
  /* cells and faces forward and back from cell j */
  for (j = 0; j < dumpdata->nce2; j++) {
    dumpdata->cfcj[j] =
      (short)((j < dumpdata->nce2 - 1) ? j + 1 : NOTVALID);
    dumpdata->cbcj[j] = (short)((j > 0) ? j - 1 : NOTVALID);
    dumpdata->ffcj[j] = (short)(j + 1);
    dumpdata->fbcj[j] = (short)(j);
  }
  /* cells and faces forward and back from face j */
  for (j = 0; j < dumpdata->nfe2; j++) {
    dumpdata->cffj[j] = (short)((j < dumpdata->nfe2 - 1) ? j : NOTVALID);
    dumpdata->cbfj[j] = (short)((j > 0) ? j - 1 : NOTVALID);
    dumpdata->fffj[j] =
      (short)((j < dumpdata->nfe2 - 1) ? j + 1 : NOTVALID);
    dumpdata->fbfj[j] = (short)((j > 0) ? j - 1 : NOTVALID);
  }
}

/* END neighbour_none()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the grid specification parameters                 */
/*-------------------------------------------------------------------*/
void read_grid(parameters_t *params)
{
  FILE *fp;                     /* Input parameter file */
  char buf[MAXSTRLEN];          /* Buffer */
  double slat, slon;            /* North-western domain coordinates */
  double elat, elon;            /* South-eastern domain coordinates */
  double x00;                   /* Grid e1 offset */
  double y00;                   /* Grid e2 offset */
  double rotn;                  /* Grid rotation */
  double xinc;                  /* e1 grid spacing */
  double yinc;                  /* e2 grid spacing */
  double lx, ly;                /* Size of the grid */
  double arc, rmin;             /* Polar grid attributes */
  double ella;                  /* Elliptic grid attributes */
  double taumax, taumin;        /* Elliptic grid attributes */
  int nsimm;                    /* Elliptic grid attributes */
  int nce1;                     /* e1 grid dimension */
  int nce2;                     /* e2 grid dimension */
  /* double nm=1852.0; *//* Meters in a nautical mile */
  /* double dm=60.0; *//* Minutes in 1 degree */
  int i, j;

  fp =params->prmfd;
  nce1 = params->nce1;
  nce2 = params->nce2;

  if (!nce1 || !nce2)
    hd_quit("Error reading grid size : nce1 and nce2 must be supplied\n");

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  /* Allocate memory for the coordinate and metric arrays which are  */
  /* at twice the final resolution.                                  */
  params->x = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  params->y = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  params->h1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  params->h2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  params->a1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  params->a2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  xinc = yinc = rotn = params->flat = params->flon = 0.0;
  prm_set_errfn(hd_silent_warn);

  memset(projection, 0, sizeof(projection));
  prm_read_char(fp, "PROJECTION", params->projection);

  /* Generate or read the coordinates, metrics and angles            */
  if (prm_read_double(fp, "SLAT", &slat) &&
      prm_read_double(fp, "ELAT", &elat) &&
      prm_read_double(fp, "SLON", &slon) &&
      prm_read_double(fp, "ELON", &elon)) {
    double mlat, mlon;          /* Lat and long of the mid-point of the
                                   grid */

    strcpy(params->gridtype, "GEOGRAPHIC_RECTANGULAR");
    strcpy(params->projection, "geographic");   /* Override projection */
    strcpy(projection, params->projection);
    mlat = 0.5 * (slat + elat);
    mlon = 0.5 * (slon + elon);

    ly = fabs(elat - slat);
    lx = fabs(elon - slon);

    xinc = cos(PI*mlat/180.0) * lx / (double)nce1;
    yinc = ly / (double)nce2;
    x00 = min(slon, elon);
    y00 = min(elat, slat);
    prm_read_double(fp, "ROTATION", &rotn);

    f_pole(mlat, mlon, rotn, &params->flon, &params->flat);
    geog_false_pole_coord(params->x, params->y, params->h1, params->h2,
                          params->a1, params->a2,
                          2 * nce1, 2 * nce2,
                          x00, y00, params->flon, params->flat, xinc / 2.0,
                          yinc / 2.0);

    params->rg = (grid_rect_t *)malloc(sizeof(grid_rect_t));
    params->rg->x0 = x00;
    params->rg->y0 = y00;
    params->rg->dx = xinc;
    params->rg->dy = yinc;
    params->rg->th = rotn;
    params->rg->sinth = sin(params->rg->th * M_PI / 180.0);
    params->rg->costh = cos(params->rg->th * M_PI / 180.0);
    params->gridcode = GRECT_DXY_ROT;

  } else if (prm_read_char(fp, "ROMS_GRID", buf)) {
    /* Read ROMS grid ahead of other types of grid specification */
    sprintf(params->roms_grid, "%s", buf);
    read_roms_grid(params->roms_grid, nce1, nce2, params->x, params->y, NULL);
    strcpy(params->gridtype, "NUMERICAL");
    grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
		     params->h1, params->h2);
    grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
		   params->a1, params->a2);
    params->gridcode = NUMERICAL;
    /* Read the various options */
    params->roms_grid_opts = ROMS_GRID_2D;
    if (prm_read_char(fp, "ROMS_GRID_OPTS", buf)) {
      if (contains_token(buf, "MASK") != NULL)
	params->roms_grid_opts |= ROMS_GRID_MASK;
      if (contains_token(buf, "BATHY") != NULL)
	params->roms_grid_opts |= ROMS_GRID_BATHY;
    }	  
  } else if (prm_read_double(fp, "DX", &xinc) &&
             prm_read_double(fp, "DY", &yinc)) {
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);
    prm_read_double(fp, "ROTATION", &rotn);
    params->rg = (grid_rect_t *)malloc(sizeof(grid_rect_t));
    params->rg->x0 = x00;
    params->rg->y0 = y00;
    params->rg->dx = xinc;
    params->rg->dy = yinc;
    params->rg->th = rotn;
    if (strlen(params->projection) > 0 && strcmp(params->projection,
"geographic") == 0) {
      strcpy(params->gridtype, "GEOGRAPHIC_RECTANGULAR");
      geog_dreckon_coord(params->x, params->y, params->h1, params->h2,
                         params->a1, params->a2, 2 * nce1, 2 * nce2,
                         x00, y00, rotn, xinc / 2.0, yinc / 2.0);
      params->gridcode = GRECT_DXY_ROT;
      if (prm_read_char(fp, "FALSE_POLE", buf)) {
	if (is_true(buf)) {
	  double mlat, mlon;
	  mlon = 0.5 * (params->x[0][0] + params->x[0][2*nce1-1]);
	  mlat = 0.5 * (params->y[0][0] + params->y[2*nce2-1][0]);
	  f_pole(mlat, mlon, rotn, &params->flon, &params->flat);
	  xinc = yinc = params->y[nce2/2+2][nce1/2] - params->y[nce2/2][nce1/2];
	  geog_false_pole_coord(params->x, params->y, params->h1, params->h2,
				params->a1, params->a2,
				2 * nce1, 2 * nce2,
				x00, y00, params->flon, params->flat,
				xinc / 2.0, yinc / 2.0);
	  params->gridcode = GRECT_DLATLON_FP;
	}
      }
    } else {
      strcpy(params->gridtype, "RECTANGULAR");
      grid_gen_rect_coord(params->x, params->y, params->h1, params->h2,
                          params->a1, params->a2, 2 * nce1, 2 * nce2,
                          x00, y00, rotn, xinc / 2.0, yinc / 2.0);
      params->rg->sinth = sin(params->rg->th * M_PI / 180.0);
      params->rg->costh = cos(params->rg->th * M_PI / 180.0);
      params->gridcode = RECT_DXY_ROT;
    }
  } else if (prm_read_double(fp, "DLAMBDA", &xinc) &&
             prm_read_double(fp, "DPHI", &yinc)) {

    strcpy(params->gridtype, "GEOGRAPHIC_RECTANGULAR");
    strcpy(params->projection, "geographic");  /* Override projection */
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);

    if (prm_read_double(fp, "POLE_LATITUDE", &params->flat) &&
        prm_read_double(fp, "POLE_LONGITUDE", &params->flon)) {
      /* Set the longitude increment relative to the grid centre rather than grid origin. */
      /* Note: ROAM compensates for this. */
      /*
      if(!(params->runmode & ROAM))
	xinc *= cos(PI*0.5*(2 * y00 + nce2 * yinc)/180.0);
      */

      geog_false_pole_coord(params->x, params->y, params->h1, params->h2,
                            params->a1, params->a2,
                            2 * nce1, 2 * nce2,
                            x00, y00, params->flon, params->flat,
                            xinc / 2.0, yinc / 2.0);

      params->gridcode = GRECT_DLATLON_FP;
    } else if (prm_read_double(fp, "ROTATION", &rotn)) {
      grid_gen_rect_coord(params->x, params->y, params->h1, params->h2,
                          params->a1, params->a2, 2 * nce1, 2 * nce2,
                          x00, y00, rotn, xinc / 2.0, yinc / 2.0);

      params->gridcode = GRECT_DLATLON_ROT;
    }
    else
      hd_quit
        ("POLE_LATITUDE/PLOT_LONGITUDE or ROTATION must be specified.\n");
    params->rg = (grid_rect_t *)malloc(sizeof(grid_rect_t));
    params->rg->x0 = x00;
    params->rg->y0 = y00;
    params->rg->dx = xinc;
    params->rg->dy = yinc;
    params->rg->th = rotn;

  } else if (prm_read_double(fp, "ARC", &arc) &&
             prm_read_double(fp, "R0", &rmin)) {
    strcpy(params->gridtype, "POLAR");
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);
    prm_read_double(fp, "ROTATION", &rotn);
    grid_gen_polar_coord(params->x, params->y, params->h1, params->h2,
                         params->a1, params->a2, 2 * nce1,
                         2 * nce2, x00, y00, rotn, arc, rmin);
    params->pg = (grid_polar_t *)malloc(sizeof(grid_polar_t));
    params->pg->x0 = x00;
    params->pg->y0 = y00;
    params->pg->arc = arc;
    params->pg->rmin = rmin;
    params->pg->rotation = rotn;
    params->gridcode = POLAR_DARC_ROT;
  } else if (prm_read_double(fp, "ELLA", &ella) &&
             prm_read_double(fp, "TAUMAX", &taumax) &&
             prm_read_double(fp, "TAUMIN", &taumin) &&
             prm_read_int(fp, "NSIMM", &nsimm)) {
    strcpy(params->gridtype, "ELLIPTIC");
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);
    prm_read_double(fp, "ROTATION", &rotn);
    grid_gen_elliptic_coord(params->x, params->y, params->h1, params->h2,
                            params->a1, params->a2, 2 * nce1, 2 * nce2,
                            x00, y00, rotn, ella, taumax, taumin,
                            2L * nsimm);
    params->gridcode = ELLIPTIC_DTAU_ROT;
  } else if (prm_read_char(fp, "XCOORDS", buf) && prm_read_char(fp, "YCOORDS", buf) ) {
    read2darray(fp, "XCOORDS", params->x, 2 * nce1 + 1, 2 * nce2 + 1);
    read2darray(fp, "YCOORDS", params->y, 2 * nce1 + 1, 2 * nce2 + 1);
    strcpy(params->gridtype, "NUMERICAL");
    grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
                     params->h1, params->h2);
    grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
                   params->a1, params->a2);
    params->gridcode = NUMERICAL;
  } else if (prm_read_char(fp, "COORDFILE", buf)) {
  	FILE* fcoord = fopen(buf,"r");
  	if(fcoord == NULL)
  		hd_quit("readparam:read_grid; Cannot read coordinates file: %s ",buf);
  	
  	read2darray(fcoord, "XCOORDS", params->x, 2 * nce1 + 1, 2 * nce2 + 1);
    read2darray(fcoord, "YCOORDS", params->y, 2 * nce1 + 1, 2 * nce2 + 1);
    strcpy(params->gridtype, "NUMERICAL");
    grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
                     params->h1, params->h2);
    grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
                   params->a1, params->a2);
    params->gridcode = NUMERICAL;
  	fclose(fcoord);
  
  } else if (prm_read_char(fp, "MOM_GRID", buf)) {
      read_mom_grid(buf, nce1, nce2, params->x, params->y);
      strcpy(params->gridtype, "NUMERICAL");
      grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
		       params->h1, params->h2);
      grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
		     params->a1, params->a2);
      params->gridcode = NUMERICAL;
  } else
    hd_quit("Cannot read grid information\n");

/* Set the time-series file projection */
  if (strlen(params->projection) > 0) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }

  /* Input parameter files are sometimes created with non-boundary */
  /* water cells located at the edge of the grid. The sparse array */
  /* requires a ghost cell adjacent to these, thus the input array */
  /* should ideally have a border of land cells. When this is not */
  /* the case (i.e. edgef>0) then the grid is extended by two cells */
  /* which are set to land and have metrics copied from the first */
  /* valid cells (e.g. i=1,nce1-2 etc).  */
  if (params->edgef)
    shift_grid(params->x, params->y, params->h1,
               params->h2, params->a1, params->a2,
               2 * nce1 + 1, 2 * nce2 + 1, params->edgef);

  /* If the the projection is specified and is geographic, then */
  /* compute the metrics and angles for the sphere.  */
  if (strncasecmp(params->projection,
                  GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) == 0) {
    grid_get_geog_metrics(params->x, params->y, 2 * nce1,
                          2 * nce2, params->h1, params->h2);
    grid_get_geog_angle(params->x, params->y, 2 * nce1,
                        2 * nce2, params->a1, params->a2);
  }
}

/* END read_grid()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add strips of land around an existing grid by linearly */
/* interpolating the grid metrics. Note that metric arrays are of    */
/* size (2*nce1)(2*nce2) hence the interpolation must assign values  */
/* over a range of 2*sft-1 cells around the grid, where sft is the   */
/* extra number of cells required.                                   */
/*-------------------------------------------------------------------*/
void
shift_grid(double **x, double **y, double **h1, double **h2, double **a1,
           double **a2, int nce1, int nce2, int sft)
{
  int i, j;
  int n, np1, np2, nm4, nm3, nm2;
  int ns = 2 * sft - 1;

  /* Shift the array values by (i+sft,j+sft) */
  for (j = nce2 - 1 - 2 * sft; j > ns; j--) {
    for (i = nce1 - 1 - 2 * sft; i > ns; i--) {
      x[j][i] = x[j - 1 - ns][i - 1 - ns];
      y[j][i] = y[j - 1 - ns][i - 1 - ns];
      h1[j][i] = h1[j - 1 - ns][i - 1 - ns];
      h2[j][i] = h2[j - 1 - ns][i - 1 - ns];
      a1[j][i] = a1[j - 1 - ns][i - 1 - ns];
      a2[j][i] = a2[j - 1 - ns][i - 1 - ns];
    }
  }

  /* Set values on the edges */
  for (j = 0; j < nce2; j++) {
    for (n = ns; n >= 0; n--) {
      np1 = n + 1;
      np2 = n + 2;
      nm4 = nce1 - n - 3;
      nm3 = nce1 - n - 2;
      nm2 = nce1 - n - 1;
      x[j][n] = intp(x[j][np1], x[j][np2], np1, np2, n);
      x[j][nm2] = intp(x[j][nm4], x[j][nm3], nm4, nm3, nm2);

      y[j][n] = intp(y[j][np1], y[j][np2], np1, np2, n);
      y[j][nm2] = intp(y[j][nm4], y[j][nm3], nm4, nm3, nm2);

      h1[j][n] = intp(h1[j][np1], h1[j][np2], np1, np2, n);
      h1[j][nm2] = intp(h1[j][nm4], h1[j][nm3], nm4, nm3, nm2);

      h2[j][n] = intp(h2[j][np1], h2[j][np2], np1, np2, n);
      h2[j][nm2] = intp(h2[j][nm4], h2[j][nm3], nm4, nm3, nm2);

      a1[j][n] = intp(a1[j][np1], a1[j][np2], np1, np2, n);
      a1[j][nm2] = intp(a1[j][nm4], a1[j][nm3], nm4, nm3, nm2);

      a1[j][n] = intp(a1[j][np1], a1[j][np2], np1, np2, n);
      a2[j][nm2] = intp(a2[j][nm4], a2[j][nm3], nm4, nm3, nm2);
    }
  }

  for (i = 0; i < nce1; i++) {
    for (n = ns; n >= 0; n--) {
      np1 = n + 1;
      np2 = n + 2;
      nm4 = nce2 - n - 3;
      nm3 = nce2 - n - 2;
      nm2 = nce2 - n - 1;

      x[n][i] = intp(x[np1][i], x[np2][i], np1, np2, n);
      x[nm2][i] = intp(x[nm4][i], x[nm3][i], nm4, nm3, nm2);

      y[n][i] = intp(y[np1][i], y[np2][i], np1, np2, n);
      y[nm2][i] = intp(y[nm4][i], y[nm3][i], nm4, nm3, nm2);

      h1[n][i] = intp(h1[np1][i], h1[np2][i], np1, np2, n);
      h1[nm2][i] = intp(h1[nm4][i], h1[nm3][i], nm4, nm3, nm2);

      h2[n][i] = intp(h2[np1][i], h2[np2][i], np1, np2, n);
      h2[nm2][i] = intp(h2[nm4][i], h2[nm3][i], nm4, nm3, nm2);

      a1[n][i] = intp(a1[np1][i], a1[np2][i], np1, np2, n);
      a1[nm2][i] = intp(a1[nm4][i], a1[nm3][i], nm4, nm3, nm2);

      a2[n][i] = intp(a2[np1][i], a2[np2][i], np1, np2, n);
      a2[nm2][i] = intp(a2[nm4][i], a2[nm3][i], nm4, nm3, nm2);
    }
  }
}

/* END shift_grid()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the tracer filtering options                      */
/*-------------------------------------------------------------------*/
void read_trfilter(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];

  sprintf(keyword, "TRACER_FILTER");
  if (prm_read_char(fp, keyword, buf)) {
    params->trfilter = 0;
    if (contains_token(buf, "NONE") != NULL)
      params->trfilter |= NONE;
    if (contains_token(buf, "FILL") != NULL)
      params->trfilter |= (TRF_FILL|TRF_FILL3D|TRF_FILL2D|TRF_FILLSED);
    if (contains_token(buf, "SMOOTH") != NULL)
      params->trfilter |= TRF_SMOO;
    if (contains_token(buf, "SHAPIRO") != NULL)
      params->trfilter |= TRF_SHAP;
    if (contains_token(buf, "MEDIAN") != NULL)
      params->trfilter |= TRF_MEDI;
    if (contains_token(buf, "SHUMAN") != NULL)
      params->trfilter |= TRF_SHUM;
    if (contains_token(buf, "FILL3D") != NULL)
      params->trfilter |= (TRF_FILL3D|TRF_FILL);
    if (contains_token(buf, "FILL2D") != NULL)
      params->trfilter |= (TRF_FILL2D|TRF_FILL);
    if (contains_token(buf, "FILLSED") != NULL)
      params->trfilter |= (TRF_FILLSED|TRF_FILL);
  }
}


/* END read_trfilter()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the windowing parameters                          */
/*-------------------------------------------------------------------*/
void read_window_info(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];
  int n, m, cc;

  sprintf(keyword, "WINDOWS");
  prm_read_int(fp, keyword, &params->nwindows);
  sprintf(keyword, "DUMP_WIN_MAP");
  prm_read_char(fp, keyword, params->wind_file);
  sprintf(keyword, "READ_WIN_MAP");
  prm_read_char(fp, keyword, params->win_file);
  if (params->nwindows > 1 && params->trasc & LAGRANGE)
    hd_quit("Multiple windows not supported with advection scheme 'LAGRANGE'\n");
  if (params->tmode & (SP_CHECK|TR_CHECK)) params->nwindows = 1;

  if (params->nwindows > 1) {
    read_win_type(params, fp);
    sprintf(keyword, "SHOW_WINDOWS");
    if (prm_read_char(fp, keyword, buf)) {
      if (is_true(buf)) {
	if (params->show_win == 0)
	  params->ntrS += 1;
	params->show_win = 1;
      }
    }
    sprintf(keyword, "WINDOW_RESET");
    prm_read_int(fp, keyword, &params->win_reset);
    params->win_size = NULL;
    sprintf(keyword, "WINDOW_SIZE");
    if (prm_read_char(fp, keyword, buf)) {
      if (params->win_size)
	d_free_1d(params->win_size);
      params->win_size = d_alloc_1d(params->nwindows);
      if (strcmp(buf, "default") == 0) {
	for (n = 0; n < params->nwindows; n++)
	  params->win_size[n] = 1.0 / (double)params->nwindows;
      } else {
	prm_skip_to_end_of_key(fp, keyword);
	for (n = 0; n < params->nwindows; n++)
	  fscanf(fp, "%lf", &params->win_size[n]);
      }
    }
    /* Read the explicit window partitions */
    if (params->nwn != NULL)
      i_free_1d(params->nwn);
    params->nwn = i_alloc_1d(params->nwindows);
    memset(params->nwn, 0, params->nwindows * sizeof(int));
    if (params->wnx) i_free_2d(params->wnx);
    if (params->wny) i_free_2d(params->wny);
    params->wnx = params->wny = NULL; m = 0;    
    /*UR to be replaced */
    cc = 0;
    for (n = 0; n < params->nwindows; n++) {
      sprintf(keyword, "WINDOW%d_POINTS", n + 1);
      if (prm_read_int(fp, keyword, &params->nwn[n])) {
	if (params->nwn[n] > m)
	  m = params->nwn[n];
	cc++;
      }
    }
    if (params->nwindows > 1 && cc > 0 && cc != params->nwindows)
      hd_warn("WINDOWS = %d but only %d WINDOW_POINTS specified. This may lead to error.\n", params->nwindows, cc);
    if (m) {
      params->wnx = i_alloc_2d(m, params->nwindows);
      params->wny = i_alloc_2d(m, params->nwindows);
      for (n = 0; n < params->nwindows; n++) {
	sprintf(keyword, "WINDOW%d_POINTS", n + 1);
	if (prm_skip_to_end_of_key(fp, keyword)) {
	  prm_flush_line(fp);
	  for (cc = 0; cc < params->nwn[n]; cc++) {
	    if (fscanf(fp, "%d %d", &params->wnx[n][cc], &params->wny[n][cc]) != 2)
	      hd_quit("params_read: Can't read i j in WINDOW%d_POINTS list.\n", n);
	    /* Flush the remainder of the line */
	    prm_flush_line(fp);
	  }
	}
      }
    }
    /*UR end to be replaced */
  }

  strcpy(params->dp_mode, "none");
  sprintf(keyword, "DP_MODE");
  if(!prm_read_char(fp, keyword, params->dp_mode)) {
    sprintf(keyword, "DPMODE");
    prm_read_char(fp, keyword, params->dp_mode);
  }
}

/* END read_window_info()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads the window type for partitioning                            */
/*-------------------------------------------------------------------*/
void read_win_type(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN];
  
  if (prm_read_char(fp, "WINDOW_TYPE", buf)) {
    int n;    
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    if (strcmp(fields[0], "STRIPE_E1") == 0)
      params->win_type = STRIPE_E1;
    if (strcmp(fields[0], "STRIPE_E2") == 0)
      params->win_type = STRIPE_E2;
    if (strcmp(fields[0], "BLOCK_E1") == 0) {
      params->win_type = BLOCK_E1;
      if (n > 1) params->win_block = atoi(fields[1]);
    }
    if (strcmp(fields[0], "BLOCK_E2") == 0) {
      params->win_type = BLOCK_E2;
      if (n > 1) params->win_block = atoi(fields[1]);
    }
    if (strcmp(fields[0], "EXPLICIT") == 0)
      params->win_type = WIN_EXP;
  }
}

/* END read_win_type()                                               */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Routine to read horizontal diffusion parameters                   */
/*-------------------------------------------------------------------*/
void read_hdiff(parameters_t *params, FILE *fp, int mode)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  prm_set_errfn(hd_silent_warn);

  params->u1vh = params->u2vh = 0.0;
  params->u1kh = params->u2kh = 0.0;

  sprintf(keyword, "U1VH");
  prm_read_double(fp, keyword, &params->u1vh);
  sprintf(keyword, "U2VH");
  prm_read_double(fp, keyword, &params->u2vh);
  sprintf(keyword, "U1KH");
  prm_read_double(fp, keyword, &params->u1kh);
  sprintf(keyword, "U2KH");
  prm_read_double(fp, keyword, &params->u2kh);

  sprintf(keyword, "DIFF_SCALE");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->diff_scale = NONE;
    if (strcmp(buf, "LINEAR") == 0)
      params->diff_scale = LINEAR;
    if (strcmp(buf, "NONLIN") == 0)
      params->diff_scale = NONLIN;
    if (strcmp(buf, "AUTO") == 0) {
      params->diff_scale = AUTO;
      /*
      if (params->u1vh > 0.0) params->u1vh *= -1.0;
      if (params->u2vh > 0.0) params->u2vh *= -1.0;
      if (params->u1kh > 0.0) params->u1kh *= -1.0;
      if (params->u2kh > 0.0) params->u2kh *= -1.0;
      */
    }
    if (strcmp(buf, "SMAG") == 0) {
      params->diff_scale = SMAG;
      params->smagorinsky = 0.1;
      if (params->u1vh > 0.0) 
	params->u1vh *= -1.0;
      else if (params->u1vh == 0.0) 
	params->u1vh = -1.0;
      if (params->u2vh > 0.0) 
	params->u2vh *= -1.0;
      else if (params->u2vh == 0.0) 
	params->u2vh = -1.0;
      if (params->u1kh > 0.0) 
	params->u1kh *= -1.0;
      else if (params->u1kh == 0.0) 
	params->u1kh = -1.0;
      if (params->u2kh > 0.0) 
	params->u2kh *= -1.0;
      else if (params->u2kh == 0.0) 
	params->u2kh = -1.0;
    }
  }

  sprintf(keyword, "SMAGORINSKY");
  /*prm_read_double(fp, keyword, &params->smagorinsky);*/
  if (prm_read_char(fp, keyword, params->smag) ) {
    int n;    
    char *fields[MAXSTRLEN * MAXNUMARGS];
    strcpy(buf, params->smag);
    n = parseline(buf, fields, MAXNUMARGS);
    if (n == 1) {
      params->smagorinsky = atof(fields[0]);
      params->sue1 = params->sue2 = 1.0;
      params->kue1 = params->kue2 = 1.0;
      params->bsue1 = params->bsue2 = 0.0;
      params->bkue1 = params->bkue2 = 0.0;
    }
    if (n == 4) {
      params->bsue1 = floor(atof(fields[0]));
      params->bsue2 = floor(atof(fields[1]));
      params->bkue1 = floor(atof(fields[2]));
      params->bkue2 = floor(atof(fields[3]));
      params->sue1 = atof(fields[0]) - params->bsue1;
      params->sue2 = atof(fields[1]) - params->bsue2;
      params->kue1 = atof(fields[2]) - params->bkue1;
      params->kue2 = atof(fields[3]) - params->bkue2;
      params->smagorinsky = 1.0;
    }
  }
  if(params->smagorinsky > 0.0) {
    if (mode)
      params->atr += 1;
    else
      params->ntr += 1;
  }

  sprintf(keyword, "SMAG_SMOOTH");
  prm_read_int(fp, keyword, &params->smag_smooth);

  sprintf(keyword, "VISC_METHOD");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->visc_method = NONE;
    if (strcmp(buf, "LAPLACIAN") == 0)
      params->visc_method = LAPLACIAN;
    if (strcmp(buf, "SIMPLE") == 0)
      params->visc_method = SIMPLE;
    if (strcmp(buf, "PRE_V794") == 0)
      params->visc_method = (LAPLACIAN|PRE794);
  }
  if (params->compatible & V794)
    params->visc_method |= PRE794;
}

/* END read_hdiff()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read wave parameters                                   */
/*-------------------------------------------------------------------*/
void read_waves(parameters_t *params, FILE *fp, int mode)
{
  char buf[MAXSTRLEN];

  /* Read the file name containing the wave information              */
  if (prm_read_char(fp, "WAVE_VARS", params->webf)) {
    if (prm_get_time_in_secs(fp, "WAVE_VARS_INPUT_DT", &params->webf_dt)) {
      params->orbital = 1;
      params->waves |= ORBITAL;
      params->ntrS += 5;
    } else {
      if (is_true(params->webf)) {
	params->orbital = 1;
	params->webf_dt = 0.0;
	params->waves |= ORBITAL;
	params->ntrS += 5;
      }
    }
    if (prm_read_char(fp, "WAVE_VARS_INTERP_TYPE", buf))
      strcpy(params->webf_interp, buf);
    else
      sprintf(params->webf_interp, "%c", '\0');
  }

  /* Read the wave feedback options                                  */
  if (prm_read_char(fp, "WAVES", buf)) {
    if (params->orbital) {
      if (contains_token(buf, "WAVE_FORCING") != NULL) {
	params->waves |= WAVE_FOR;
	params->ntrS += 2;
	if (params->tendf) params->ntrS += 2;
      }
      if (contains_token(buf, "TAN_RADIATION") != NULL) {
	params->waves |= TAN_RAD;
	params->ntrS += 2;
	if (params->tendf) params->ntrS += 2;
      }
      if (contains_token(buf, "STOKES") != NULL) {
	params->waves |= (STOKES|STOKES_DRIFT|STOKES_MIX|STOKES_WIND);
	params->ntrS += 6;
      }
      if (contains_token(buf, "STOKES_DRIFT") != NULL) {
	params->waves |= (STOKES|STOKES_DRIFT);
	params->ntrS += 6;
      }
      if (contains_token(buf, "STOKES_MIX") != NULL) {
	params->waves |= (STOKES|STOKES_MIX);
	params->ntrS += 6;
      }
      if (contains_token(buf, "STOKES_WIND") != NULL) {
	params->waves |= (STOKES|STOKES_WIND);
	params->ntrS += 6;
      }
      if (contains_token(buf, "SPECTRAL") != NULL) {
	params->waves |= SPECTRAL;
	if (mode)
	  params->atr += 2;
	else
	  params->ntr += 2;
      }
      if (contains_token(buf, "BOT_STRESS") != NULL) {
	params->waves |= BOT_STR;
	params->ntrS += 1;
      }
      if (contains_token(buf, "VERT_MIX") != NULL) {
	params->waves |= VERTMIX;
	if (contains_token(buf, "MIX_JONES") != NULL)
	  params->waves |= JONES;
	else if (contains_token(buf, "MIX_WOM") != NULL)
	  params->waves |= WOM;
	else if (contains_token(buf, "MIX_BVM") != NULL)
	  params->waves |= BVM;
      }
    } else
      hd_warn("params_read: WAVES requires WAVE_VARS file. Ignoring WAVES...\n");
  }

  /* Read the fetch file                                             */
  prm_read_char(fp, "FETCH", params->fetch_init);
}

/* END read_waves()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read short wave radiation parameters                   */
/*-------------------------------------------------------------------*/
void read_swr(parameters_t *params, FILE *fp, int mode)
{
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];
  /* Extinction and transmission for water types, from POM User's    */
  /* Guide, Mellor (1992).                                           */
  double swr_type[2][6] = {
    {0.0, 0.037, 0.042, 0.056, 0.073, 0.127},
    {0.0, 0.320, 0.310, 0.290, 0.260, 0.240}
  };

  /* Check for swr bottom absorption; default = 1                    */
  params->swr_type = 0;
  sprintf(keyword, "SWR_BOT_ABSORB");
  if (prm_read_char(fp, keyword, params->swr_babs))
    params->ntrS++;
  /* Check for swr partitioning; default = none                      */
  sprintf(keyword, "SWR_ATTENUATION");
  if (prm_read_char(fp, keyword, params->swr_attn)) {
    /* Check for 2 layer swr attenuation                             */
    sprintf(keyword, "SWR_ATTENUATION_DEEP");
    if (prm_read_char(fp, keyword, params->swr_attn1)) {
      /* Get the partitioning fraction                               */
      sprintf(keyword, "SWR_FRACTION");
      if (!(prm_read_char(fp, keyword, params->swr_tran))) {
	hd_quit_and_dump
	  ("params_read() : Split swr attenuation requires SWR_FRACTION set.\n");
      }
      params->swr_type = (SWR_2D|SWR_SPLIT);
      params->ntrS+=3;
    } else {
      /* Check the transmission coefficient : default = 1            */
      sprintf(keyword, "SWR_TRANSMISSION");
      if (!(prm_read_char(fp, keyword, params->swr_tran))) {
	hd_warn("params_read() : SWR_TRANSMISSION set to 1.0.\n");
	strcpy(params->swr_tran, "1.0");
      }
      params->swr_type = SWR_2D;
      params->ntrS+=2;
    }
  }
  sprintf(keyword, "SWR_ATTENUATION3D");
  if (prm_read_char(fp, keyword, params->swr_attn)) {
    /* Check the transmission coefficient : default = 1            */
    sprintf(keyword, "SWR_TRANSMISSION");
    if (!(prm_read_char(fp, keyword, params->swr_tran))) {
      hd_warn("params_read() : SWR_TRANSMISSION set to 1.0.\n");
      strcpy(params->swr_tran, "1.0");
    }
    params->ntrS+=1;
    if (mode)
      params->atr += 1;
    else
      params->ntr += 1;
    params->swr_type = SWR_3D;
  }

  if (strlen(params->swr_attn) && !strlen(params->swr_tran))
    strcpy(params->swr_tran, "1.0");
  sprintf(keyword, "WATER_TYPE");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "TYPE_I") == 0)
      params->water_type = TYPE_I;
    if (strcmp(buf, "TYPE_IA") == 0)
      params->water_type = TYPE_IA;
    if (strcmp(buf, "TYPE_IB") == 0)
      params->water_type = TYPE_IB;
    if (strcmp(buf, "TYPE_II") == 0)
      params->water_type = TYPE_II;
    if (strcmp(buf, "TYPE_III") == 0)
      params->water_type = TYPE_III;
  }
  if (params->water_type != NONE) {
    sprintf(params->swr_attn, "%f", swr_type[0][params->water_type]);
    sprintf(params->swr_tran, "%f", swr_type[1][params->water_type]);
    params->swr_type = SWR_2D;
    params->ntrS+=2;
  }
  if (strlen(params->swr_attn) && !strlen(params->swr_babs)) {
    strcpy(params->swr_babs, "1.0");
    params->ntrS++;
  }
}

/* END read_swr()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read bulk heatflux formulation parameters              */
/*-------------------------------------------------------------------*/
void read_hf_bulk(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN], buf[MAXSTRLEN], buf1[MAXSTRLEN];

  sprintf(keyword, "HEATFLUX_REFH");
  if (!(prm_read_double(fp, keyword, &params->zref)))
    params->zref = 10.0;
  /* Note : the input codes 0 - 4 are included for backwards */
  /* compatibility.                                          */
  if (prm_read_char(fp, "BULK_SCHEME", buf)) {
    if((strcmp(buf, "L&P") == 0) || (strcmp(buf, "1") == 0))
      params->bulkf = LARGEPOND;
    else if((strcmp(buf, "B") == 0) || (strcmp(buf, "4") == 0))
      params->bulkf = BUNKER;
    else if((strcmp(buf, "K/W") == 0) || (strcmp(buf, "3") == 0))
      params->bulkf = KITIAG;
    else if((strcmp(buf, "Ko") == 0) || (strcmp(buf, "0") == 0))
      params->bulkf = KONDO;
    else if((strcmp(buf, "M") == 0) || (strcmp(buf, "2") == 0))
      params->bulkf = MASAG;
  } else
    params->bulkf = KONDO;
  params->hfadf = 0;
  sprintf(keyword, "HEATFLUX_ADVECT");
  if (prm_read_char(fp, keyword, buf1))
    params->hfadf = is_true(buf1);
  sprintf(keyword, "HEATFLUX_TEMP");
  prm_read_char(fp, keyword, params->hftemp);
  sprintf(keyword, "HEATFLUX_TEMP_DT");
  prm_get_time_in_secs(fp, keyword, &params->hftemp_dt);
  sprintf(keyword, "HEATFLUX_TC");
  prm_get_time_in_secs(fp, keyword, &params->hftc);
}

/* END read_hf_bulk()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read backwards compatibility parameter                 */
/*-------------------------------------------------------------------*/
void read_compatible(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  sprintf(keyword, "COMPATIBLE");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "V794") != NULL)
      params->compatible |= V794;
    if (contains_token(buf, "V1246") != NULL)
      params->compatible |= V1246;
    if (contains_token(buf, "V1562") != NULL)
      params->compatible |= V1562;
    if (contains_token(buf, "V1598") != NULL) {
      params->compatible |= V1598;
      if (params->momsc & WTOP_O4) {
	params->momsc &= ~WTOP_O4;
	params->momsc |= WTOP_O2;
      }
      if (params->momsc & ZERO_DRYK)
	params->momsc &= ~ZERO_DRYK;
      if (params->momsc & SHAPIRO)
	params->momsc &= ~SHAPIRO;
    }
    if (contains_token(buf, "V1652") != NULL)
      params->compatible |= V1652;
    if (contains_token(buf, "V1670") != NULL)
      params->compatible |= V1670;
    if (contains_token(buf, "V1957") != NULL)
      params->compatible |= V1957;
    if (contains_token(buf, "V4201") != NULL)
      params->compatible |= V4201;
    if (contains_token(buf, "V5342") != NULL)
      params->compatible |= V5342;
    if (contains_token(buf, "V5895") != NULL)
      params->compatible |= V5895;
  }
}

/* END read_compatible()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read cell centers for grid refinement                  */
/*-------------------------------------------------------------------*/
void read_zoom(parameters_t *params, FILE *fp)
{
  int n;                     /* Number of times to map to get a face */
  char buf[MAXSTRLEN];
  /*
  if (params->nwindows == 1)
    hd_quit("GRID_REFINEMENT requires more than 1 window %d\n", params->nwindows);
  */

  /* Read the locations of the cell centers in (i,j) space for the   */
  /* zoom zone.  */
  prm_set_errfn(quit);
  if (prm_read_int(fp, "ZOOM_FACTOR", &params->zoomf))
    params->zmfe1 = params->zmfe2 = params->zoomf;
  else {
    prm_read_int(fp, "ZOOM_FACTOR_E1", &params->zmfe1);
    prm_read_int(fp, "ZOOM_FACTOR_E2", &params->zmfe2);
  }
  if (prm_read_int(fp, "ZOOM_POINTS", &params->nzoom)) {
    params->zci = i_alloc_1d(params->nzoom + 1);
    params->zcj = i_alloc_1d(params->nzoom + 1);
    for (n = 1; n <= params->nzoom; n++) {
      if (fscanf(fp, "%d %d", &params->zci[n], &params->zcj[n]) != 2)
	hd_quit("params_read: Can't read i j in 'zoom_cells' list.\n");
    }
  }
  if (prm_read_char(fp, "ZOOM_HR_ZONE", buf)) {
    int i, j, is, ie, js, je;
    if (sscanf(buf, "(%d,%d)-(%d,%d)",&is, &js, &ie, &je) == 4) {
      params->nzoom = 0;
      for (j = params->zoomf / 2; j < params->nce2; j += params->zoomf)
	for (i = params->zoomf / 2; i < params->nce1; i += params->zoomf)
	  if (i < is || i > ie || j < js || j > je)
	    params->nzoom++;
      params->zci = i_alloc_1d(params->nzoom + 1);
      params->zcj = i_alloc_1d(params->nzoom + 1);
      n = 1;
      for (j = params->zoomf / 2; j < params->nce2; j += params->zoomf)
	for (i = params->zoomf / 2; i < params->nce1; i += params->zoomf)
	  if (i < is || i > ie || j < js || j > je) {
	    params->zci[n] = i;
	    params->zcj[n] = j;
	    n++;
	  }
    } else
      params->zoomf = 1;
  }
  if (prm_read_char(fp, "ZOOM_CR_ZONE", buf)) {
    int i, j, is, ie, js, je;
    if (sscanf(buf, "(%d,%d)-(%d,%d)",&is, &js, &ie, &je) == 4) {
      params->nzoom = 0;
      for (j = params->zoomf / 2; j < params->nce2; j += params->zoomf)
	for (i = params->zoomf / 2; i < params->nce1; i += params->zoomf)
	  if (i >= is && i <= ie && j >= js && j <= je)
	    params->nzoom++;
      params->zci = i_alloc_1d(params->nzoom + 1);
      params->zcj = i_alloc_1d(params->nzoom + 1);
      n = 1;
      for (j = params->zoomf / 2; j < params->nce2; j += params->zoomf)
	for (i = params->zoomf / 2; i < params->nce1; i += params->zoomf)
	  if (i >= is && i <= ie && j >= js && j <= je) {
	    params->zci[n] = i;
	    params->zcj[n] = j;
	    n++;
	  }
    } else
      params->zoomf = 1;
  }
  prm_set_errfn(hd_silent_warn);
}

/* END read_zoom()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in the TOTALS diagnostic                          */
/*-------------------------------------------------------------------*/
void read_totals(parameters_t *params, FILE *fp)
{
  int n, m, ntr;
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];

  if (prm_read_char(fp, "TOTALS", buf)) {
    n = parseline(buf, fields, MAXNUMARGS);
    params->totals = is_true(fields[0]);
    ntr = params->ntot = n + 1; /* Always include T & S */
    if (ntr) {
      params->totname = (char **)malloc(ntr * sizeof(char *));      
      for (m = 0; m < ntr; m++) {
	params->totname[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);	
	if (m == 0)
	  strcpy(params->totname[m], "temp");
	else if (m == 1)
	  strcpy(params->totname[m], "salt");
	else
	  strcpy(params->totname[m], fields[m - 1]);
      }
    }
    if (prm_read_char(fp, "TOTALS_DT", buf))
      tm_scale_to_secs(buf, &params->totals_dt);
    else
      params->totals_dt = 3600.0;
  }
}

/* END read_totals()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in eta relaxation                                 */
/*-------------------------------------------------------------------*/
void read_eta_relax(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], tu0[MAXSTRLEN], tu1[MAXSTRLEN];
  char fname[] = "eta_relaxation_file";
  char dt_name[] = "eta_relaxation_input_dt";
  char tc_name[] = "eta_relaxation_time_constant";

  sprintf(params->etarlxn, "%c", '\0');
  if (prm_read_char(fp, fname, params->etarlxn) > 0) {
    if (prm_read_char(fp, dt_name, buf) > 0 &&
	prm_read_char(fp, tc_name, params->etarlxtcs) > 0) {
      double d0, d1, r0, r1;
      tm_scale_to_secs(buf, &params->etarlxdt);
      if (sscanf(params->etarlxtcs, "%lf %s", &params->etarlxtc, buf) == 2)
	tm_scale_to_secs(params->etarlxtcs, &params->etarlxtc);
      params->etarlx = RELAX;
      params->ntrS++;
      if (sscanf(params->etarlxtcs, "%s %lf %lf %s %lf %lf %s", 
		 buf, &d0, &r0, tu0, &d1, &r1, tu1) == 7) {
	params->etarlx |= ETA_ADPT;
	params->ntrS += 2;
      }
    } else if (prm_read_char(fp, dt_name, buf) > 0) {
      tm_scale_to_secs(buf, &params->etarlxdt);
      params->ntrS++;
    } else {
      params->etarlxdt = 3600.0;
      params->etarlxtc = 0.0;
      params->etarlx = ALERT;
      params->ntrS++;
    }
  }
}
/* END read_eta_relax()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read in eta relaxation                                 */
/*-------------------------------------------------------------------*/
void read_vel_relax(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN];
  char fname[] = "vel_relaxation_file";
  char dt_name[] = "vel_relaxation_input_dt";
  char tc_name[] = "vel_relaxation_time_constant";

  params->velrlx = NONE;
  sprintf(params->velrlxn, "%c", '\0');
  if (prm_read_char(fp, fname, params->velrlxn) > 0) {
    if (prm_read_char(fp, dt_name, buf) > 0 &&
	prm_read_char(fp, tc_name, params->velrlxtcs) > 0) {
      tm_scale_to_secs(buf, &params->velrlxdt);
      params->velrlx = RELAX;
    }
  }
}

/* END read_vel_relax()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read debugging location                                */
/*-------------------------------------------------------------------*/
void read_debug(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  int n, m = 3;

  if (prm_read_char(fp, "DEBUG_LOC", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    params->dbi = atoi(fields[0]);
    params->dbj = atoi(fields[1]);
    params->dbk = atoi(fields[2]);
    if (m != n) params->dbgf = 0;
    while (m != n) {
      if (strcmp(fields[m], "append") == 0) {
	params->dbgf |= D_APPEND;
	m++;
      } else if(strcmp(fields[m], "after") == 0) {
	params->dbgf |= D_AFTER;
	sprintf(buf, "%s %s", fields[m+1], fields[m+2]);
	tm_scale_to_secs(buf, &params->dbgtime);
	m += 3;
      } else if (strcmp(fields[m], "step") == 0) {
	params->dbgf |= D_STEP;
	m++;
      } else {
	hd_warn("read_debug: can't understand keyword %s\n", fields[m]);
	m++;
      }
    }
  }
}

/* END read_debug()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read decorrelation diagnostic                          */
/*-------------------------------------------------------------------*/
void read_decorr(parameters_t *params, FILE *fp, int mode)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];

  sprintf(keyword, "DECORR_LENGTH");
  if (prm_read_char(fp, keyword, buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int n = parseline(buf, fields, MAXNUMARGS);
    if (n >= 2) {
      strcpy(params->decv, fields[0]);
      if (strcmp(params->decv, "eta") == 0)
	params->decf = DEC_ETA;
      else if (strcmp(params->decv, "u1") == 0)
	params->decf = DEC_U1;
      else if (strcmp(params->decv, "u2") == 0)
	params->decf = DEC_U2;
      else
	params->decf = DEC_TRA;
      if (params->decf & (DEC_ETA))
	params->ntrS += 2;
      else {
	if (mode)
	  params->atr += 2;
	else
	  params->ntr += 2;
      }
      params->decorr = atoi(fields[1]);
      if (n == 3)
	strcpy(params->decs, fields[2]);
      else
	strcpy(params->decs, "m");
    }
  }
}

/* END read_decorr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create the AVHRR file list                             */
/*-------------------------------------------------------------------*/
void create_avhrr_list(parameters_t *params)
{
  FILE *ap;
  char infile[MAXSTRLEN], path[MAXSTRLEN], date[MAXSTRLEN];
  double start;            /* Model start time                       */
  double stop;             /* Model stop time                        */
  int nfiles, nf = 0;      /* Number of SST files                    */
  int fid;                 /* netcdf file handle                     */
  int ncerr;               /* netcdf error code                      */
  int i, n;                /* Counter                                */
  int product;             /* Number of days in composite            */
  int ys, mos, ds, hs, mis, ss;  /* Start year, month, day           */
  int ye, moe, de, he, mie, se;  /* End year, month, day             */
  int y, m, d, lp;               /* Year, month, day                 */
  double ep;                     /* Epoch                            */
  int dayno[2][13] = {           /* Days in the month                */
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
  };
  int lyr[2] = {365, 366};
  char files[MAXNUMTSFILES][MAXSTRLEN];

  /* Open the file list                                              */
  strcpy(path, params->avhrr_path);
  if(strcmp(path,"avhrr_list.txt") == 0) {
    if ((ap = fopen(path, "r")) == NULL)
      hd_warn("Can't open AVHRR file '%s'\n", path);
    else {
      fclose(ap);
      return;
    }
  }
  ap = fopen("avhrr_list.txt", "w");

  /* Get start and end year, month, day                              */
  ep = tm_time_to_julsecs(params->timeunit);
  tm_scale_to_secs(params->start_time, &start);
  tm_scale_to_secs(params->stop_time, &stop);
  nfiles = (ceil(stop) - floor(start)) / 86400;
  start = start / 86400 + ep;;
  tm_to_julsecs(start, &ys, &mos, &ds, &hs, &mis, &ss);
  stop = stop / 86400.0 + ep;
  tm_to_julsecs(stop, &ye, &moe, &de, &he, &mie, &se);

  /* If the path points to a netcdf file, include the file in the    */
  /* list.                                                           */
  n = parseline(path, (char **)files, MAXNUMTSFILES);
  i = strlen(path);
  if (nf > 1 ||
      (path[i - 3] == '.' && path[i - 2] == 'n' && path[i - 1] == 'c')) {
    nfiles = n;
    nf = 1;
  }

  /* Print the file header                                           */
  fprintf(ap, "multi-netcdf-version 1.0\n\n");
  fprintf(ap, "nfiles %d\n", nfiles);

  /* Print the file list                                             */
  if (nf) {
    for (i = 0; i < nfiles; i++) {
      strcpy(infile, ((char **)files)[i]);
      /* Check the file exists, if so include in the list            */
      if ((ncerr = nc_open(infile, NC_NOWRITE, &fid)) != NC_NOERR) {
	hd_warn("Can't find AVHRR SST input file :\n");
	hd_warn("  File = %s\n", infile);
      } else {
	fprintf(ap, "file0.filename %s\n", infile);
	nc_close(fid);
      }
    }
  } else {
    /* Make a list of files for each day.                            */
    /* Get the composite length (days).                              */
    for(i = 0 ; i < strlen(path) - 4; i++)
      if(path[i] == 'c' && path[i + 1] == 'o' &&
	 path[i + 2] == 'm' && path[i + 3] == 'p')
	break;
    if (path[i + 5] == 'd')
      sprintf(infile, "%c", path[i + 4]);
    else
      sprintf(infile, "%c%c", path[i + 4], path[i + 5]);
    product = atoi(infile);

    /* Print the list of files                                       */
    y = ys; m = mos; d = ds; n = 0;
    for (i = 0; i < nfiles; i++) {
      sprintf(date, "%d", y);
      if (m < 10)
	sprintf(date, "%s0%d", date, m);
      else
	sprintf(date, "%s%d", date, m);
      if (d < 10)
	sprintf(date, "%s0%d", date, d);
      else
	sprintf(date, "%s%d", date, d);
      lp = y % 4 == 0 && (y % 100 != 0 || y % 400 == 0);
      d += 1;
      if (d > dayno[lp][m]) {
	m++;
	d = 1;
      }
      if (d > lyr[lp] || m > 12) {
	y++;
	m = 1;
	d = 1;
      }
      if (path[strlen(path) - 1] == '/')
	sprintf(infile, "%sSST%dd_%s.nc", path, product, date);
      else
	sprintf(infile, "%s/SST%dd_%s.nc", path, product, date);

      /* Check the file exists, if so include in the list            */
      if ((ncerr = nc_open(infile, NC_NOWRITE, &fid)) != NC_NOERR) {
	hd_warn("Can't find AVHRR SST input file for date %s\n", date);
	hd_warn("  File = %s\n", infile);
      } else {
	fprintf(ap, "file%d.filename %s\n", n, infile);
	nc_close(fid);
	n++;
      }
    }
    if(!n) {
      hd_quit("No AVHRR SST files found from %d/%d/%d to %d/%d/%d\n",
	      ds, mos, ys, de, moe, ye);
    } else {
      rewind(ap);
      fprintf(ap, "multi-netcdf-version 1.0\n\n");
      fprintf(ap, "nfiles %d\n", n);
    }
  }
  fclose(ap);
}

/* END create_avhrr_list()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create the GHRSST file list                            */
/*-------------------------------------------------------------------*/
void create_ghrsst_list(parameters_t *params)
{
  FILE *ap;
  char infile[MAXSTRLEN], path[MAXSTRLEN], date[MAXSTRLEN], fname[MAXSTRLEN];
  char fname2[MAXSTRLEN], daystr[MAXSTRLEN];
  double start;            /* Model start time                       */
  double stop;             /* Model stop time                        */
  int nfiles, nf = 0;      /* Number of SST files                    */
  int fid;                 /* netcdf file handle                     */
  int ncerr;               /* netcdf error code                      */
  int i, n;                /* Counter                                */
  int product;             /* Number of days in composite            */
  int ys, mos, ds, hs, mis, ss;  /* Start year, month, day           */
  int ye, moe, de, he, mie, se;  /* End year, month, day             */
  int y, m, d, day, lp;          /* Year, month, day                 */
  int is_mnc;
  double ep;                     /* Epoch                            */
  int tday;                      /* Transition day                   */
  int mode = 0;                  /* Transition mode                  */
  int dayno[2][13] = {           /* Days in the month                */
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
  };
  int lyr[2] = {365, 366};
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];

  /* Get the input arguments                                         */
  nfiles = parseline(params->ghrsst_path, fields, MAXNUMARGS);

  /* Check if all names are .mnc files                               */
  is_mnc = 0;
  for (n = 0; n < nfiles; n++) {
    strcpy(path, fields[n]);
    m = strlen(path);
    if (path[m-4] == '.' && path[m-3] == 'm' &&
	path[m-2] == 'n' && path[m-1] == 'c')
      is_mnc += 1;
  }
  is_mnc = (is_mnc == nfiles) ? 1 : 0;

  if (is_mnc) {
    for (n = 0; n < nfiles; n++) {
      strcpy(path, fields[n]);
      if (n == nfiles-1)
	sprintf(fname, "%s(ghrsst=analysed_sst)", fields[n]);
      else
	sprintf(fname, "%s(ghrsst=analysed_sst) ", fields[n]);
      if ((ap = fopen(path, "r")) == NULL)
	hd_warn("Can't open GHRSST file '%s'\n", path);
      else {
	fclose(ap);
      }
    }
    params->ghrsst = 2;
    strcpy(params->ghrsst_path, fname);
    return;
  } else {
    if (nfiles == 2) {
      strcpy(path, fields[0]);
      strcpy(fname, fields[1]);
      mode = 1;
    } else if (nfiles == 4) {
      strcpy(path, fields[0]);
      strcpy(fname, fields[1]);
      strcpy(fname2, fields[2]);
      tday = atoi(fields[3]);
      mode = 2;
    } else {
      hd_warn("create_ghrsst_list: GHRSST input format = <path> <filename>\n");
      return;
    }
  }

  /* Open the file list                                              */
  ap = fopen("ghrsst_list.mnc", "w");
  params->ghrsst = 1;

  /* Get start and end year, month, day                              */
  ep = tm_time_to_julsecs(params->timeunit);
  tm_scale_to_secs(params->start_time, &start);
  tm_scale_to_secs(params->stop_time, &stop);
  nfiles = (ceil(stop) - floor(start)) / 86400;
  start = start / 86400 + ep;;
  tm_to_julsecs(start, &ys, &mos, &ds, &hs, &mis, &ss);
  stop = stop / 86400.0 + ep;
  tm_to_julsecs(stop, &ye, &moe, &de, &he, &mie, &se);

  /* Print the file header                                           */
  fprintf(ap, "multi-netcdf-version 1.0\n\n");
  fprintf(ap, "nfiles %d\n", nfiles);

  /* Make a list of files for each day.                              */
  /* Print the list of files                                         */
  y = ys; m = mos; d = ds; n = 0;
  for (i = 0; i < nfiles; i++) {
    day = yrday(y, m, d);
    if (day < 10)
      sprintf(daystr, "00%d", day);
    else if (day < 100)
      sprintf(daystr, "0%d", day);
    else
      sprintf(daystr, "%d", day);
    if (path[strlen(path) - 1] == '/')
      sprintf(infile, "%s%d/%s", path, y, daystr);
    else
      sprintf(infile, "%s/%d/%s", path, y, daystr);
    
    sprintf(date, "%d", y);
    if (m < 10)
      sprintf(date, "%s0%d", date, m);
    else
      sprintf(date, "%s%d", date, m);
    if (d < 10)
      sprintf(date, "%s0%d", date, d);
    else
      sprintf(date, "%s%d", date, d);

    lp = y % 4 == 0 && (y % 100 != 0 || y % 400 == 0);
    d += 1;
    if (d > dayno[lp][m]) {
      m++;
      d = 1;
    }
    if (d > lyr[lp] || m > 12) {
      y++;
      m = 1;
      d = 1;
    }
    if (mode == 1)
      sprintf(infile, "%s/%s%s", infile, date, fname);
    else {
      if (day <= tday)
	sprintf(infile, "%s/%s%s", infile, date, fname);
      else
	sprintf(infile, "%s/%s%s", infile, date, fname2);
    }

    /* Check the file exists, if so include in the list            */
    /*
    if ((ncerr = nc_open(infile, NC_NOWRITE, &fid)) != NC_NOERR) {
      hd_warn("Can't find GHRSST SST input file for date %s\n", date);
      hd_warn("  File = %s\n", infile);
    } else {
    */
      fprintf(ap, "file%d.filename %s\n", n, infile);
      nc_close(fid);
      n++;

  }
  if(!n) {
    hd_quit("No GHRSST SST files found from %d/%d/%d to %d/%d/%d\n",
	    ds, mos, ys, de, moe, ye);
  } else {
    rewind(ap);
    fprintf(ap, "multi-netcdf-version 1.0\n");
    fprintf(ap, "nfiles %d\n", n);
  }
  fclose(ap);
}

/* END create_ghrsst_list()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read explicit maps from the parameter file             */
/*-------------------------------------------------------------------*/
void read_explicit_maps(parameters_t *params, FILE *fp)
{
  int n, cc;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "MAP_POINTS_E1");
  if (prm_read_int(fp, keyword, &params->nemx)) {
    params->emisx = i_alloc_1d(params->nemx);
    params->emjsx = i_alloc_1d(params->nemx);
    params->emidx = i_alloc_1d(params->nemx);
    params->emjdx = i_alloc_1d(params->nemx);
    params->ktx = i_alloc_1d(params->nemx);
    params->kbx = i_alloc_1d(params->nemx);
    for (cc = 0; cc < params->nemx; cc++) {
      prm_next_line(buf, MAXSTRLEN, fp);
      n = sscanf(buf, "%d %d : %d %d : %d %d", &params->emisx[cc],
		 &params->emjsx[cc],&params->emidx[cc], &params->emjdx[cc],
		 &params->ktx[cc], &params->kbx[cc]);
      if (n == 4) {
	params->ktx[cc] = params->nz - 1;
	params->kbx[cc] = 0;
      }
      params->ktx[cc] = min(params->nz - 1, params->ktx[cc]);
      params->kbx[cc] = max(0, params->kbx[cc]);
      if (n !=4 && n != 6) {
	hd_warn("params_read: line %d, format = is js : ie je : kt kb.\n", cc);
	hd_quit("params_read: Can't read i j in MAP_POINTS_E1 list.\n");
      }
    }
  }
  sprintf(keyword, "MAP_POINTS_E2");
  if (prm_read_int(fp, keyword, &params->nemy)) {
    params->emisy = i_alloc_1d(params->nemy);
    params->emjsy = i_alloc_1d(params->nemy);
    params->emidy = i_alloc_1d(params->nemy);
    params->emjdy = i_alloc_1d(params->nemy);
    params->kty = i_alloc_1d(params->nemy);
    params->kby = i_alloc_1d(params->nemy);
    for (cc = 0; cc < params->nemy; cc++) {
      prm_next_line(buf, MAXSTRLEN, fp);
      n = sscanf(buf, "%d %d : %d %d : %d %d", &params->emisy[cc],
		 &params->emjsy[cc],&params->emidy[cc], &params->emjdy[cc],
		 &params->kty[cc], &params->kby[cc]);
      if (n == 4) {
	params->kty[cc] = params->nz - 1;
	params->kby[cc] = 0;
      }
      params->kty[cc] = min(params->nz - 1, params->kty[cc]);
      params->kby[cc] = max(0, params->kby[cc]);
      if (n !=4 && n != 6) {
	hd_warn("params_read: line %d, format = is js : ie je : kt kb.\n", cc);
	hd_quit("params_read: Can't read i j in MAP_POINTS_E2 list.\n");
      }
    }
  }
}

/* END read_explicit_maps()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read exclusion points from library processes           */
/*-------------------------------------------------------------------*/
void read_exclude_points(parameters_t *params, FILE *fp)
{
  int n, cc;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf;

  read_blocks(fp, "EXCLUDE_BGCSED", &params->prex, &params->prxi, &params->prxj);
  if (params->prex) {
    params->prxf = i_alloc_1d(params->prex + 1);
    for (n = 1; n <= params->prex; n++) {
      params->prxf[n] |= (EX_BGC|EX_SED);
    }
  }

  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "EXCLUDE_PROCESS_POINTS");
  if (prm_read_int(fp, keyword, &params->prex)) {
    params->prxi = i_alloc_1d(params->prex);
    params->prxj = i_alloc_1d(params->prex);
    params->prxf = i_alloc_1d(params->prex);
    for (cc = 1; cc <= params->prex; cc++) {
      params->prxf[cc] = 0;
      prm_next_line(buf, MAXSTRLEN, fp);
      n = sscanf(buf, "%d %d", &params->prxi[cc], &params->prxj[cc]);
      if (n < 2)
	hd_quit("read_exclude_points: Can't read i j in EXCLUDE_PROCESS_POINTS list.\n");
      nf = parseline(buf, fields, MAXNUMARGS);
      for (n = 2; n < nf; n++) {
	if (strcmp(fields[n], "EX_BGC") == 0) params->prxf[cc] |= EX_BGC;
	if (strcmp(fields[n], "EX_SED") == 0) params->prxf[cc] |= EX_SED;
	if (strcmp(fields[n], "EX_TRAN") == 0) params->prxf[cc] |= EX_TRAN;
	if (strcmp(fields[n], "EX_WAVE") == 0) params->prxf[cc] |= EX_WAVE;
	if (strcmp(fields[n], "EX_TRST") == 0) params->prxf[cc] |= EX_TRST;
      }
    }
  }
}

/* END read_exclude_points()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the mean diagnostics                              */
/*-------------------------------------------------------------------*/
void read_means(parameters_t *params, FILE *fp, int mode)
{
  int n, cc, ntr = 0;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char trname[MAXSTRLEN];

  prm_set_errfn(hd_silent_warn);
  prm_read_int(fp, "NTRACERS", &ntr);

  sprintf(keyword, "MEAN");
  params->means = 0;
  sprintf(params->means_dt, "%c", '\0');
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL)
      params->means |= NONE;
    if (contains_token(buf, "VEL3D") != NULL)
      params->means |= VEL3D;
    if (contains_token(buf, "VEL2D") != NULL)
      params->means |= VEL2D;
    if (contains_token(buf, "FLUX") != NULL)
      params->means |= FLUX;
    if (contains_token(buf, "TENDENCY") != NULL)
      params->means |= TENDENCY;
    if (contains_token(buf, "WIND") != NULL)
      params->means |= WIND;
    if (contains_token(buf, "ETA") != NULL)
      params->means |= ETA_M;
    if (contains_token(buf, "KZ_M") != NULL)
      params->means |= KZ_M;
    if (contains_token(buf, "TS") != NULL)
      params->means |= TS;
    if (contains_token(buf, "TS|MMM") != NULL)
      params->means |= (TS|MMM);
    if (contains_token(buf, "TIDAL") != NULL)
      params->means |= TIDAL;
    if (contains_token(buf, "TRANSPORT") != NULL)
      params->means |= TRANSPORT;
    if (contains_token(buf, "VOLFLUX") != NULL)
      params->means |= VOLFLUX;
    for (n = 0; n < ntr; n++) {
      sprintf(keyword, "TRACER%1.1d.name", n);
      if (prm_read_char(fp, keyword, trname) > 0) {
	if (contains_token(buf, trname) != NULL) {
	  strcpy(params->means_tra, trname);
	  sprintf(keyword, "TRACER%1.1d.type", n);
	  if (prm_read_char(fp, keyword, trname))
	    if (strcmp(trname, "WC2D") == 0)
	      params->means |= MTRA2D;
	    else
	      params->means |= MTRA3D;
	  else
	    params->means |= MTRA3D;
	}
      }
    }

    if (params->means & ETA_M)
      params->ntrS += 1;
    if (params->means & WIND)
      params->ntrS += 2;
    if (params->means & VEL2D)
      params->ntrS += 2;
    if (params->means & MTRA2D)
      params->ntrS += 1;
    if (params->means & VEL3D) {
      if(mode)
	params->atr += 3;
      else
	params->ntr += 3;
    }
    if (params->means & TS) {
      if(mode)
	params->atr += 2;
      else
	params->ntr += 2;
    }
    if (params->means & MTRA3D) {
      if(mode)
	params->atr += 1;
      else
	params->ntr += 1;
    }
    if (params->means & KZ_M) {
      if(mode)
	params->atr++;
      else
	params->ntr++;
    }
    if (params->means & VOLFLUX) {
      if(mode)
	params->atr += 2;
      else
	params->ntr += 2;
    }
    strcpy(params->means_dt, params->stop_time);
    sprintf(keyword, "MEAN_DT");
    if (!prm_read_char(fp, keyword, params->means_dt))
      sprintf(params->means_dt, "%f days", get_run_length(fp) + 1);
    sprintf(keyword, "MEAN_OFFSET");
    prm_read_char(fp, keyword, params->means_os);
  } else
    params->means = NONE;

  if (params->trout) {
    if (!prm_read_char(fp, "MEAN_DT", params->means_dt) &&
	!prm_read_char(fp, "TRANS_DT", params->means_dt))
      strcpy(params->means_dt, "1 hour");
    if (params->means & NONE) params->means &= ~NONE;
    if (!(params->means & VEL3D)) {
      params->means |= VEL3D;
      if(mode)
	params->atr += 3;
      else
	params->ntr += 3;
    }
    if (!(params->means & KZ_M)) {
      params->means |= KZ_M;
      if(mode)
	params->atr++;
      else
	params->ntr++;
    }
  }
  if (params->trout) {
    if (params->means & NONE) params->means &= ~NONE;
    if (!(params->means & VOLFLUX)) {
      params->means |= VOLFLUX;
      if(mode)
	params->atr += 2;
      else
	params->ntr += 2;
    }
  }
  if (params->tmode & SP_FFSL) {
    if(mode)
      params->atr += 2;
    else
      params->ntr += 2;
  }
}

/* END read_means()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads sediment transport information                              */
/*-------------------------------------------------------------------*/
#if defined(HAVE_SEDIMENT_MODULE)
void read_sediments(parameters_t *params, FILE *fp, int *ntr)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];

  params->do_sed = 0;
  sprintf(keyword, "DO_SEDIMENTS");
  if (prm_read_char(fp, keyword, buf)) {
    if (!params->sednz)
      hd_quit("Sediments requires NSEDLAYERS > 0\n");
    if (strcmp(buf, "DO") == 0)
      params->do_sed |= LIB_DO;
    else if (strcmp(buf, "WRITE") == 0)
      params->do_sed |= LIB_WRITE;
    else if (is_true(buf))
      params->do_sed |= LIB_DO;
    if (params->do_sed) {
      params->ntrS += 1;
      if (prm_read_char(fp, "SED_VARS", params->sed_vars)) {
	int n = count_sed_classes(params->sed_vars); /* Sed classes + tss */
	*ntr += n;  
	params->nsed += n + 2; /* sed + salt and temp */
	params->ntrS += count_sed_2d();    /* ustrcw_skin and depth_sed */
      }
      prm_read_char(fp, "SED_VARS_ATTS", params->sed_defs);
    }
  }
}
#endif
/* END read_sediments()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads sediment transport information                              */
/*-------------------------------------------------------------------*/
#if defined(HAVE_ECOLOGY_MODULE)
void read_ecology(parameters_t *params, FILE *fp, int *ntr)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];

  prm_set_errfn(hd_silent_warn);
  params->do_eco = 0;
  if (prm_read_char(fp, "DO_ECOLOGY", buf)) {
    if (!params->sednz)
      hd_quit("Ecology requires NSEDLAYERS > 0\n");
    if (strcmp(buf, "DO") == 0)
      params->do_eco |= LIB_DO;
    else if (strcmp(buf, "WRITE") == 0)
      params->do_eco |= LIB_WRITE;
    else if (is_true(buf))
      params->do_eco |= LIB_DO;
    prm_get_time_in_secs(fp, "ECOLOGY_DT", &params->ecodt);
    if (params->do_eco) {
      prm_read_char(fp, "ECO_VARS_ATTS", params->eco_defs);
      params->ntrS += 1;
      // Build temp ecology, if needed
      if (strlen(params->eco_defs))
	params->pre_eco = ecology_pre_build(params->eco_vars, 
					    params->eco_defs, params->prmfd);
      if (strlen(params->eco_vars)) {
	int n = count_eco_3d(params->eco_vars);
	*ntr += n;  
	params->nsed += n;
	params->ntrS += count_eco_2d(params->eco_vars);
      }
    }
  }
}
#endif

/* END read_ecology()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check if T/S is subject to relaxation, and increment   */
/* the auto tracers if so.                                           */
/*-------------------------------------------------------------------*/
void check_TS_relax(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char name[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n, ntr;

  prm_set_errfn(hd_silent_warn);
  if (prm_read_int(fp, "NTRACERS", &ntr)) {
    for (n = 0; n < ntr; n++) {
      sprintf(keyword, "TRACER%1.1d.name", n);
      if (prm_read_char(fp, keyword, name) > 0) {
	sprintf(keyword, "TRACER%1.1d.relaxation_file", n);
	if (prm_read_char(fp, keyword, buf) > 0) {
	  sprintf(keyword, "TRACER%1.1d.relaxation_input_dt", n);
	  if (prm_read_char(fp, keyword, buf) > 0) {
	    sprintf(keyword, "TRACER%1.1d.relaxation_time_constant", n);
	    if (prm_read_char(fp, keyword, buf) > 0) {
	      int nf = parseline(buf, fields, MAXNUMARGS);
	      if (strcmp(name, "salt") == 0) {
		if (nf <= 2) {
		  params->rsalt = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		} else { /*Adaptive relaxation */
		  params->rsalt = RLX_ADPT;
		  params->atr += 2;
		  params->ntr += 2;
		}
	      }
	      if (strcmp(name, "temp") == 0) {
		if (nf <= 2) {
		  params->rtemp = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		} else { /*Adaptive relaxation */
		  params->rtemp = RLX_ADPT;
		  params->atr += 2;
		  params->ntr += 2;
		}
	      }
	    }
	  }
	} else {
	  sprintf(keyword, "relax_%s", name);
	  if (prm_read_char(fp, keyword, buf)) {
	    char fname[MAXSTRLEN], dt[MAXSTRLEN], unit[MAXSTRLEN];
	    double in_dt;
	    if (sscanf(buf,"%s %s %lf %s", fname, dt, &in_dt, unit) == 4) {
	      if (strcmp(name, "salt") == 0) {
		params->rsalt = 1;
		params->atr += 1;
		params->ntr += 1;
	      }
	      if (strcmp(name, "temp") == 0) {
		params->rtemp = 1;
		params->atr += 1;
		params->ntr += 1;
	      }
	    }
	  }
	}
      }
    }
  }
}

/* END check_TS_relax()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Tracers relating to certain diagnostics etc. are automatically    */
/* generated if the relevant flag is set. If these tracers are       */
/* manually defined in the input parameter file, then the manual     */
/* definition over-rides the auto generation, and the defined        */
/* tracer_info is used.                                              */
/*-------------------------------------------------------------------*/
void tracer_setup(parameters_t *params, FILE *fp)
{
  int ntr, m, n;
  char buf[MAXSTRLEN];

  prm_read_char(fp, "TRACER_DATA", params->tracerdata);

  check_TS_relax(params, fp);

  tracer_read(fp, NULL, WATER, hd_quit, hd_warn, hd_silent_warn,
              &params->ntr, &params->trinfo_3d);

  /* Check for duplicate tracers (i.e. automatically generated */
  /* tracers that are also defined in the input parameter file).  */
  ntr = 0;
  if ((tracer_find_index("temp", params->ntr, params->trinfo_3d)) >= 0)
    ntr++;
  if ((tracer_find_index("salt", params->ntr, params->trinfo_3d)) >= 0)
    ntr++;

  if (params->means & VEL3D) {
    if ((tracer_find_index("u1mean", params->ntr, params->trinfo_3d)) >= 0)
      ntr++;
    if (tracer_find_index("u2mean", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("wmean", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("Kzmean", params->ntr, params->trinfo_3d) >= 0 &&
	(params->means & KZ_M || params->trout))
      ntr++;
  }
  if (params->means & TS) {
    if (tracer_find_index("temp_mean", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("salt_mean", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->means & VOLFLUX) {
    if ((tracer_find_index("u1vmean", params->ntr, params->trinfo_3d)) >= 0)
      ntr++;
    if (tracer_find_index("u2vmean", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->trflux, "NONE") != 0) {
    if (tracer_find_index("flux_e1", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("flux_e2", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("flux_w", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("flux_kz", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->trperc, "NONE") != 0) {
    sprintf(buf, "percentile_%s", params->trperc);
    if (tracer_find_index(buf, params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->numbers & SOUND) {
    if (tracer_find_index("sound", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("sound_channel", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->trflsh) {
    if (tracer_find_index("flush", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strlen(params->trage)) {
    if (tracer_find_index("age", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->smagorinsky > 0.0) {
    if (tracer_find_index("smagorinsky", params->ntr, params->trinfo_3d) >=
        0)
      ntr++;
  }
  if (strlen(params->trtend)) {
    if (tracer_find_index("tra_adv", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("tra_vdif", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("tra_ncon", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->tendf) {
    if (tracer_find_index("mom_balance", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_adv", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_hdif", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_vdif", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_btp", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_bcp", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u1_cor", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_adv", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_hdif", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_vdif", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_btp", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_bcp", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("u2_cor", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->mixsc, "k-e") == 0) {
    if (tracer_find_index("tke", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("diss", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & OTEMP) {
    if (tracer_find_index("otemp", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & OSALT) {
    if (tracer_find_index("osalt", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & OVELU) {
    if (tracer_find_index("ovelu", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & OVELV) {
    if (tracer_find_index("ovelv", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->mixsc, "k-w") == 0 || strcmp(params->mixsc, "W88") == 0) {
    if (tracer_find_index("tke", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("omega", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
      strcmp(params->mixsc, "harcourt") == 0) {
    if (tracer_find_index("tki", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("tki_l", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("lscale", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("Kq", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->rtemp)
    if (tracer_find_index("rtemp", params->ntr, params->trinfo_3d) >= 0) {
      params->rtemp = 0;
      ntr++;
    }
  if (params->rsalt)
    if (tracer_find_index("rsalt", params->ntr, params->trinfo_3d) >= 0) {
      params->rsalt = 0;
      ntr++;
    }
  if (params->dhwf & DHW_NOAA) {
    if (tracer_find_index("dhw", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("dhd", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("dhwc", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }

  /* Re-order tracers defined in the input parameter file to account */
  /* for duplicate tracers.  */
  m = params->atr - ntr;
  for (n = params->atr; n < params->ntr; n++) {
    tracer_copy(&params->trinfo_3d[m], &params->trinfo_3d[n]);
    params->trinfo_3d[m].n = m;
    if (params->trinfo_3d[m].flag & (SC_ST|SC_PT))
      params->trinfo_3d[m].scale -= (double)ntr;
    m++;
  }
  params->atr -= ntr;
  params->ntr -= ntr;

  /*
  for (n = 0; n < params->ntr; n++) {
    printf("%d %s\n",n,params->trinfo_3d[n].name);
  }
  if (params->trflux != -1) {
    for (n = params->atr; n < params->ntr; n++) {
      if(params->trflux == params->trinfo_3d[n].m) {
	params->trflux = n;
	break;
      }
    }
  }
  */
}

/* END tracer_setup()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads in a list of (i,j) as a series of blocks                    */
/*-------------------------------------------------------------------*/
void read_blocks(FILE *fp, char *key, int *nb, int **listi, int **listj)
{
  char buf[MAXSTRLEN];
  int nlist, n, nn;
  int i, j, is, js, ie, je;
  int *xi, *xj;

  if (prm_read_int(fp, key, &nlist)) {
    *nb = 0;
    for (n = 0; n < nlist; n++) {
      prm_next_line(buf, MAXSTRLEN, fp);
      if (sscanf(buf, "(%d,%d)-(%d,%d)",&is, &js, &ie, &je) == 4) {
	for (j = js; j <= je; j++)
	  for (i = is; i <= ie; i++)
	    *nb += 1;
      } else if (sscanf(buf, "%d %d",&is, &js) == 2) {
	*nb += 1;
      }
    }
    xi = i_alloc_1d(*nb+1);
    xj = i_alloc_1d(*nb+1);
    nn = 1;
    prm_read_int(fp, key, &nlist);
    for (n = 0; n < nlist; n++) {
      prm_next_line(buf, MAXSTRLEN, fp);
      if (sscanf(buf, "(%d,%d)-(%d,%d)",&is, &js, &ie, &je) == 4) {
	for (j = js; j <= je; j++)
	  for (i = is; i <= ie; i++) {
	    xi[nn] = i;
	    xj[nn] = j;
	    nn++;
	  }
      } else if (sscanf(buf, "%d %d",&is, &js) == 2) {
	xi[nn] = is;
	xj[nn] = js;
	nn++;
      }
    }
    *listi = xi;
    *listj = xj;
  } else
    *nb = 0;
}

/* END read_blocks()                                                 */
/*-------------------------------------------------------------------*/


void f_pole(double ilat, double ilon, double ang, double *olat,
            double *olon)
{
  double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, az = 0.0;
  double qdist = (2.0 * M_PI * RADIUS) / 4.0;

  x1 = DEG2RAD(ilon);
  y1 = DEG2RAD(ilat);
  az = DEG2RAD(ang);

  geod_fwd_sodanos(x1, y1, az, qdist, RADIUS, ECC, &x2, &y2);
  *olat = RAD2DEG(x2);
  *olon = RAD2DEG(y2);
}


double intp(double a, double b, int xs, int xe, int x)
{
  return ((double)(x - xs) * (b - a) / ((double)(xe - xs)) + a);
}

double intpf(double a, double b, double xs, double xe, double x)
{
  return ((x - xs) * (b - a) / ((xe - xs)) + a);
}

/*-------------------------------------------------------------------*/
char *otime(double dt, char *tag)
{
  if (dt >= 86400)
    sprintf(tag, "%-6.1f days", dt / 86400);
  else if (dt >= 3600 && dt < 86400)
    sprintf(tag, "%-6.1f hours", dt / 3600);
  else if (dt >= 60 && dt < 3600)
    sprintf(tag, "%-6.1f minutes", dt / 60);
  else if (dt < 60)
    sprintf(tag, "%-6.1f seconds", dt);
  return (tag);
}


char *tf(int m)
{
  switch (m) {
  case 0:
    return ("NO");
  case 1:
    return ("YES");
  }
  return ("NO");
}

char *btname(int m)
{
  switch (m) {
  case U1BDRY:
    return ("u1");
  case U2BDRY:
    return ("u2");
  case U1BDRY | U2BDRY:
    return ("velocity");
  }
  return ("NULL");
}

char *bdryname(int m)
{
  switch (m) {
  case 1:
    return ("FILEIN");
  case 2:
    return ("CUSTOM");
  case 3:
    return ("LINEAR");
  case 4:
    return ("NOGRAD");
  case 8:
    return ("LINEXT");
  case 16:
    return ("POLEXT");
  case 32:
    return ("GRAVTY");
  case 64:
    return ("ORLANS");
  case 128:
    return ("CAMOBR");
  case 256:
    return ("MILLER");
  case 512:
    return ("VERTIN");
  case 1024:
    return ("CYCLIC");
  case 2048:
    return ("NOTHIN");
  case 4096:
    return ("RAYMND");
  case 8192:
    return ("TIDEBC");
  case 16384:
    return ("CLAMPD");
  case 32768:
    return ("STATIS");
  case 65536:
    return ("UPSTRM");
  case 65536 | 1:
    return ("UPSTRM");
  case 131072:
    return ("TIDALM");
  case 262144:
    return ("CYCLED");
  }
  return ("NULL");
}

char *substepname(int m)
{
  switch (m) {
  case ORIGINAL:
    return ("ORIGINAL");
  case SUB_STEP:
    return ("SUB-STEP");
  case SUB_STEP_NOSURF:
    return ("SUB-STEP-NOSURF");
  case SUB_STEP_TRACER:
    return ("SUB-STEP-TRACER");
  case MULTI_DT:
    return ("MULTI-DT");
  case NONE:
    return ("NONE");
  }
  return ("NULL");
}

char *mixname(int m)
{
  switch (m) {
  case DENS_MIX:
    return ("DENS_MIX");
  case TKE_MIX:
    return ("TKE_MIX");
  case TEMP_MIX:
    return ("TEMP_MIX");
  case NONE:
    return ("NONE");
  }
  return ("NULL");
}

char *heatfluxname(int m)
{
  switch (m) {
  case ORIGINAL:
    return ("ORIGINAL");
  case ADVANCED:
    return ("BULK");
  case INVERSE:
    return ("INVERSE");
  case NET_HEAT:
    return ("NET_HEAT");
  case COMP_HEAT:
    return ("COMP_HEAT");
  case COMP_HEAT_NONE:
    return ("COMP_HEAT_NONE");
  case SURF_RELAX:
    return ("SURF_RELAX");
  case NONE:
    return ("NONE");
  }
  return ("NULL");
}

char *meansname(int m, char *buf)
{
  sprintf(buf, "%c", '\0');
  if (m & FLUX)
    sprintf(buf, "%s FLUX ", buf);
  if (m & VEL3D)
    sprintf(buf, "%s VEL3D ", buf);
  if (m & VEL2D)
    sprintf(buf, "%s VEL2D ", buf);
  if (m & TS)
    sprintf(buf, "%s TS ", buf);
  if (m & WIND)
    sprintf(buf, "%s WIND ", buf);
  if (m & TENDENCY)
    sprintf(buf, "%s TENDENCY ", buf);
  if (m & ETA_M)
    sprintf(buf, "%s ETA ", buf);
  if (m & NONE)
    sprintf(buf, "NONE");
  return (buf);
}

char *trtypename(int m, char *buf)
{
  sprintf(buf, "%c", '\0');
  if (m == WATER)
    return("WC3D");
  if (m == INTER)
    return("WC2D");
  if (m == SEDIM)
    return("SED");
  if (m & WATER)
    sprintf(buf, "%s WC3D ", buf);
  if (m & INTER)
    sprintf(buf, "%s WC2D ", buf);
  if (m & SEDIM)
    sprintf(buf, "%s SED ", buf);
  if (m & HYDRO)
    sprintf(buf, "%s HYDRO ", buf);
  if (m & SEDIMENT)
    sprintf(buf, "%s SEDIMENT ", buf);
  if (m & ECOLOGY)
    sprintf(buf, "%s ECOLOGY ", buf);
  if (m & WAVE)
    sprintf(buf, "%s WAVE ", buf);
  if (m & TRACERSTAT)
    sprintf(buf, "%s TRACERSTAT ", buf);
  if (m & PROGNOSTIC)
    sprintf(buf, "%s PROGNOSTIC ", buf);
  if (m & DIAGNOSTIC)
    sprintf(buf, "%s DIAGNOSTIC ", buf);
  if (m & PARAMETER)
    sprintf(buf, "%s PARAMETER ", buf);
  if (m & FORCING)
    sprintf(buf, "%s FORCING ", buf);
  return (buf);
}

char *cflname(int m)
{
  switch (m) {
  case PASSIVE:
    return ("PASSIVE");
  case PASSIVE|WVEL:
    return ("PASSIVE|WVEL");
  case ACTIVE:
    return ("ACTIVE");
  case ACTIVE|WVEL:
    return ("ACTIVE|WVEL");
  case ACTIVE3D:
    return ("ACTIVE3D");
  case ACTIVE3D|WVEL:
    return ("ACTIVE3D|WVEL");
  case NONE:
    return ("NONE");
  }
  return ("NULL");
}

char *adname(int m)
{
  switch (m) {
  case ORDER1:
    return ("ORDER1");
  case ORDER2:
    return ("ORDER2");
  case QUICKEST_AD:
    return ("QUICKEST_AD");
  case ORDER4:
    return ("ORDER4");
  case QUICKEST:
    return ("QUICKEST");
  case QUICKEST_CO:
    return ("QUICKEST_CO");
  case VANLEER:
    return ("VANLEER");
  case ORDER2_UW:
    return ("ORDER2_UPWIND");
  case LAGRANGE:
    return ("LAGRANGE");
  case FFSL:
      return ("FFSL");
  case ANGULAR:
    return ("ANGULAR");
  case WIMPLICIT:
    return ("WIMPLICIT");
  case LAGRANGE|VANLEER:
    return ("LAGRANGE|VANLEER");
  }
  return ("NULL");
}

double get_restart_time(char *filename, char *itimeunits) {
  double restart_time = NaN;
  int cdfid = 0;
  int ndims = 0;
  int vid;

  ncw_open(filename, NC_NOWRITE, &cdfid);
  vid = ncw_var_id(cdfid, "t");
  if ((nc_inq_varndims(cdfid, vid, &ndims) == NC_NOERR) && (ndims == 1)) {
     int dimid = 0;
     size_t start = 0;
     size_t nt = 0;

     nc_inq_vardimid(cdfid, vid, &dimid);
     nc_inq_dimlen(cdfid, dimid, &nt);

     if (nt > 0) {
        double *t = d_alloc_1d(nt);
        nc_get_vara_double(cdfid,vid,&start,&nt,t);
        if (itimeunits) {
           char ftimeunits[MAXSTRLEN];
           memset(ftimeunits,0,MAXSTRLEN);
           nc_get_att_text(cdfid,vid,"units",ftimeunits);
           tm_change_time_units(ftimeunits, itimeunits, t, nt);
        }

        restart_time = t[nt-1];
        d_free_1d(t);
     }
  } else
     hd_quit("get_restart_time: Expected restart file time variable to have single dimension.\n");

  return restart_time;
}


double get_nrt_time(char *filename, char *itimeunits) {
  double nrt_time = NaN;
  int cdfid = 0;
  int ndims = 0;
  int vid;


  ncw_open(filename, NC_NOWRITE, &cdfid);
  vid = ncw_var_id(cdfid, "t");
  if ((nc_inq_varndims(cdfid, vid, &ndims) == NC_NOERR) && (ndims == 1)) {
     int dimid = 0;
     size_t start = 0;
     size_t nt = 0;

     nc_inq_vardimid(cdfid, vid, &dimid);
     nc_inq_dimlen(cdfid, dimid, &nt);

     if (nt > 0) {
        double *t = d_alloc_1d(nt);
        nc_get_vara_double(cdfid,vid,&start,&nt,t);
        if (itimeunits) {
           char ftimeunits[MAXSTRLEN];
           memset(ftimeunits,0,MAXSTRLEN);
           nc_get_att_text(cdfid,vid,"units",ftimeunits);
           tm_change_time_units(ftimeunits, itimeunits, t, nt);
        }

        nrt_time = t[nt-1];
        d_free_1d(t);
     }
  } else
     hd_quit("get_nrt_time: Expected nrt file time variable to have single dimension.\n");

  return nrt_time;
}


void get_bdry_params(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  sprintf(keyword, "WRITE_BDRY");
  prm_read_char(fp, keyword, params->bdry_file);
  prm_set_errfn(hd_quit);
  sprintf(keyword, "NBOUNDARIES");
  prm_read_int(fp, keyword, &params->nobc);

  sprintf(keyword, "BDRY_PATH");
  if (prm_read_char(fp, keyword, buf)) {
    if (endswith(buf,"/"))
      sprintf(keyword, "%s", buf);
    else
      sprintf(keyword, "%s/", buf);
    strcpy(params->bdrypath, keyword);
  }

  /* Read river loads directory location for RECOM */
  if (params->roammode & (A_RECOM_R1|A_RECOM_R2))
    prm_read_char(fp, "RIVER_LOADS_DIR", params->rivldir);

}


int get_smoothing(char *list, char *tag) {
  int n;
  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  int nf;

  if (!strlen(list)) return(0);

  strcpy(buf, list);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (nf <= 0) {
    hd_warn("get_smoothing : SMOOTH_VARS incorrectly specified.\n");
    return(0);
  }

  for (n = 0; n < nf; n++) {
    char *tok;
    tok = strtok(fields[n], ":");
    if (strcmp(tok, tag) == 0) {
      tok = strtok(NULL, ":");
      return(atoi(tok));
    }
  }
  return(0);
}

double get_scaling(char *list, char *tag) {
  int n;
  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  int nf;

  if (!strlen(list)) return(1.0);

  strcpy(buf, list);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (nf <= 0) {
    hd_warn("get_scaling : SMOOTH_VARS incorrectly specified.\n");
    return(1.0);
  }

  for (n = 0; n < nf; n++) {
    char *tok;
    tok = strtok(fields[n], ":");
    if (strcmp(tok, tag) == 0) {
      tok = strtok(NULL, ":");
      return(atof(tok));
    }
  }
  return(1.0);
}

int decode_tag(char *list, char *tag, char *value) {
  int n;
  char *fields[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  int nf;

  sprintf(value, "%c", '\0');
  if (!strlen(list)) return(0);

  strcpy(buf, list);
  nf = parseline(buf, fields, MAXNUMARGS);
  if (nf <= 0) {
    hd_warn("decode_tag : tag incorrectly specified.\n");
    return(0);
  }

  for (n = 0; n < nf; n++) {
    char *tok;
    tok = strtok(fields[n], ":");
    if (strcmp(tok, tag) == 0) {
      tok = strtok(NULL, ":");
      strcpy(value, tok);
      return(1);
    }
  }
  return(0);
}

double get_run_length(FILE *fp)
{
  double start, stop, len;

  prm_get_time_in_secs(fp, "START_TIME", &start);
  prm_get_time_in_secs(fp, "STOP_TIME", &stop);
  return((stop - start) / 86400.0);
}

void get_output_path(parameters_t *params, FILE *fp) {

  char buf[MAXSTRLEN];
  char path[MAXSTRLEN];
  int sequence = -1;
  int n;

  if (prm_read_char(fp, "OutputPath", path)) {
    if (prm_read_char(fp, "SEQUENCE", params->sequence)) {
      if (sscanf(params->sequence, "%d", &n) == 1)
	sequence = n;
      else {
	FILE *sp;
	sp = fopen(params->sequence, "r");
	prm_read_int(sp, "Run #", &n);
	sequence = n + 1;
	sprintf(params->sequence, "%d", n+1);
	fclose(sp);
      }
      if (!params->runno) params->runno = (double)sequence;
      if (endswith(path,"/"))
	sprintf(buf, "%srun%d/", path, sequence);
      else
	sprintf(buf, "%s/run%d/", path, sequence);
      strcpy(path, buf);
    } else {
      if (endswith(path,"/"))
	sprintf(buf, "%s", path);
      else
	sprintf(buf, "%s/", path);
    }
    strcpy(params->opath, buf);
  }
}

int is_wet(parameters_t *params, double bathy) {
  /* Note; LANDCELL has not been yet set, hence use etamax/layers[nz] */
  if (bathy < params->layers[params->nz] && bathy < params->etamax &&
      bathy != NOTVALID)
    return(1);
  else
    return(0);
}

int is_land(parameters_t *params, double bathy) {
  if (bathy != NOTVALID && 
      (bathy >= params->layers[params->nz] || bathy >= params->etamax))      
    return(1);
  else
    return(0);
}

int is_outside(parameters_t *params, double bathy) {
  if (bathy == NOTVALID)
    return(1);
  else
    return(0);
}

void reset_bathy(parameters_t *params, int io, int jo, double val)
{
  int i, j, n = 0;
  int is = 0;
  int js = 0;
  int ie = params->nce1;
  int je = params->nce2;
  if (params->edgef) {
    is = js = params->edgef;
    ie = params->nce1 - params->edgef;
    je = params->nce2 - params->edgef;
  }
  for (j = js; j < je; j++) {
    for (i = is; i < ie; i++) {
      if (i == io && j == jo)
	params->bathy[n] = -val;
      n++;
    }
  }
}

int check_river_bdry(parameters_t *params, FILE *fp)
{
  int n, nobc;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  sprintf(keyword, "NBOUNDARIES");
  if (prm_read_int(fp, keyword, &nobc)) {
    for (n = 0; n < nobc; n++) {
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.u1", n);
      if (prm_read_char(fp, keyword, buf)) {
	if (strcmp(buf, "u1flowbdry") == 0) return 1;
      }
      sprintf(keyword, "BOUNDARY%1d.CUSTOM.u2", n);
      if (prm_read_char(fp, keyword, buf)) {
	if (strcmp(buf, "u2flowbdry") == 0) return 1;
      }
      sprintf(keyword, "BOUNDARY%1d.BCOND0", n);
      if (prm_read_char(fp, keyword, buf)) {
	char *fields[MAXSTRLEN * MAXNUMARGS];
	(void)parseline(buf, fields, MAXNUMARGS);
	if (strcmp(fields[0], "RIVER") == 0) return 1;
      }
    }
  }
  sprintf(keyword, "RIVER0");
  if (prm_read_char(fp, keyword, buf)) {
    return 1;
  }
  return 0;
}

int check_bdry_options(parameters_t *params, FILE *fp)
{
  int n, nobc;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  int ret = 0;

  sprintf(keyword, "NBOUNDARIES");
  if (prm_read_int(fp, keyword, &nobc)) {
    for (n = 0; n < nobc; n++) {
      sprintf(keyword, "BOUNDARY%1d.OPTIONS", n);
      if (prm_read_char(fp, keyword, buf)) {
	if (contains_token(buf, "NONE") != NULL) {
	  ret = NONE;
	} else {
	  ret = 0;
	  if (contains_token(buf, "UPSTRM") != NULL)
	    ret |= OP_UPSTRM;
	  if (contains_token(buf, "GEOSTR") != NULL)
	    ret |= OP_GEOSTR;
	  if (contains_token(buf, "RLEN") != NULL)
	    ret |= OP_RLEN;
	  if (contains_token(buf, "DYNAMIC_HC") != NULL)
	    ret |= OP_DYNAHC;
	  if (contains_token(buf, "NO_HDIFF") != NULL)
	    ret |= OP_NOHDIF;
	}
      }
    }
  }
  return ret;
}

