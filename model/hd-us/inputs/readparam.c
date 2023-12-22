/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/readparam.c
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
 *  $Id: readparam.c 7463 2023-12-13 03:51:41Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "hd.h"
/* JIGSAW grid generation */
#ifdef HAVE_JIGSAWLIB
#include "lib_jigsaw.h"
#endif

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
#define O_BDRY    16
#define I_BDRY    32

/*-------------------------------------------------------------------*/
/* Valid bathymetry netCDF dimension names                           */
static char *bathy_dims[6][6] = {
  {"botz", "i_centre", "j_centre", "x_centre", "y_centre", "standard"},
  {"height", "lon", "lat", "lon", "lat", "nc_bathy"},
  {"height", "longitude", "latitude", "longitude", "latitude", "nc_bathy"},
  {"height", "nlon", "nlat", "lon", "lat", "nc_bathy"},
  {"elevation", "lon", "lat", "lon", "lat", "nc_bathy"},
  {NULL, NULL, NULL, NULL, NULL, NULL}
};

int nce1;
int nce2;

static int read_bathy_from_file(FILE *fp, parameters_t *params);
static double *read_bathy_from_db(FILE *fp, parameters_t *params);
static void mask_bathy_with_coast_s(const char *cfname, parameters_t *params);
static void mask_bathy_with_coast_us(const char *cfname, parameters_t *params);
static int read_bathy_from_nc(parameters_t *params, char *fname);
static void read_bathy_from_sparse_nc(parameters_t *params, char *fname);
static void read_bathy_from_sparse_bty(parameters_t *params, char *fname);
void reset_bathy(parameters_t *params, int io, int jo, double val);
int read2darray(FILE * fp, char *label, double **array, long n1, long n2);
void shift_grid(double **x, double **y, double **h1, double **h2,
                double **a1, double **a2, int nce1, int nce2, int sft);
void eta_init(geometry_t *geom, parameters_t *params, master_t *master);
void vel_init(geometry_t *geom, parameters_t *params, master_t *master);
void read_explicit_maps(parameters_t *params, FILE *fp);
void read_hdiff(parameters_t *params, FILE *fp, int mode);
void read_vdiff(parameters_t *params, FILE *fp, int mode);
void read_eta_relax(parameters_t *params, FILE *fp);
void read_vel_relax(parameters_t *params, FILE *fp);
void read_compatible(parameters_t *params, FILE *fp);
void read_means(parameters_t *params, FILE *fp, int mode);
void read_debug(parameters_t *params, FILE *fp);
void read_win_type(parameters_t *params, FILE *fp);
void read_profile(parameters_t *params, FILE *fp);
void read_advect(parameters_t *params, FILE *fp);
void read_turb(parameters_t *params, FILE *fp);
void read_history(parameters_t *params, FILE *fp);
void check_TS_relax(parameters_t *params, FILE *fp);
void get_output_path(parameters_t *params, FILE *fp);
void baro_vortex(master_t *master, geometry_t *geom, int io, int jo);
int check_river_bdry(parameters_t *params, FILE *fp);
int check_bdry_options(parameters_t *params, FILE *fp);
int count_auto_OBC_us(parameters_t *params, int **nib, int ***jc, int **mask);
int set_auto_OBC_us(parameters_t *params, int *nib, int **jc, int *mask);
int count_auto_OBC(parameters_t *params, double **bathy, int **nib, 
		    int ***mask, int *onr, int *ors, 
		    int **ib, int **jb, int **rb, int **rf, double **rhc);
void set_auto_OBC(parameters_t *params, double **bathy, int bdry, int *nib, 
		  int **mask, int *onr, int *ors, 
		  int *ib, int *jb, int *rb, int *rf, double *rhc);
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
static double *read_bathy_from_db_us(FILE *fp, parameters_t *params);
double* topo_get_z_for_grid_us(topo_files_t *tfs, double *gx, double *gy,
			       int ns2, topo_hint_t hint);
void remove_OBC_corner(parameters_t *params);
void bathy_interp_us(parameters_t *params, char *fname, char *i_rule, int mode);
void bathy_interp_s(parameters_t *params, char *fname, int mode);
void write_tile_dump(master_t *master, FILE *op, char *fname, char *vars, int tile, int obc,
		     double *regionid, int nreg, int *reg, int **flag, int nobc, int fm, int mode);
void decode_id(parameters_t *params, char *ifile, char *name, double *grd_id,
	       double *hyd_id, double *sed_id, double *bgc_id);
void decode_idc(parameters_t *params, char *ifile, char *name, char *grd_id,
		char *hyd_id, char *sed_id, char *bgc_id);

#ifdef HAVE_JIGSAWLIB
void create_jigsaw_mesh(coamsh_t *cm, jigsaw_msh_t *J_mesh, jigsaw_msh_t *J_hfun,
			int powf, int stproj, int filef);
void hfun_from_bathy(parameters_t *params, char *fname, coamsh_t *cm, jigsaw_msh_t *J_hfun, int mode);
void hfun_from_coast(parameters_t *params, coamsh_t *cm, jigsaw_msh_t *J_hfun, int mode);
void create_hex_radius(double crad, double x00, double y00,
		       double rmin, double rmax, double gscale,
		       jigsaw_msh_t *J_mesh, int powf);
void create_ellipse(parameters_t *params, long int nce1, long int nce2, 
		    double x00, double y00, double flon, double flat,
		    double xinc, double yinc, 
		    double elf, double ores,
		    jigsaw_msh_t *J_mesh, jigsaw_msh_t *J_hfun);
void neighbour_finder_j(jigsaw_msh_t *msh, int ***neic);
#endif

/*-------------------------------------------------------------------*/
/* Routine to set default parameters                                 */
/*-------------------------------------------------------------------*/
void set_default_param(parameters_t *params)
{
  int n;
  sprintf(params->codeheader, "%c", '\0');
  sprintf(params->parameterheader, "%c", '\0');
  sprintf(params->reference, "%c", '\0');
  sprintf(params->notes, "%c", '\0');
  sprintf(params->trl, "%c", '\0');
  sprintf(params->opath, "%c", '\0');
  sprintf(params->trkey, "%c", '\0');
  sprintf(params->sequence, "%c", '\0');
  sprintf(params->bdrypath, "%c", '\0');
  sprintf(params->tracerdata, "%c", '\0');
  sprintf(params->rivldir, "%c", '\0');
  sprintf(params->trans_dt, "%c", '\0');
  sprintf(params->runnoc, "%c", '\0');
  sprintf(params->mesh_reorder, "%c", '\0');
  params->mrf = 0;
  params->runno = 0;
  params->history = NONE;
  sprintf(params->runcode, "%c", '\0');
  sprintf(params->rev, "%c", '\0');
  params->momsc = RINGLER|WTOP_O2|WIMPLICIT|PV_ENEUT;
  /*params->trasc = VANLEER|HIORDER;*/
  params->trasc = VANLEER;
  params->ultimate = 0;
  params->kinetic = ORDER2;
  params->kfact = 1.0;
  params->osl = L_LSLIN;
  params->rkstage = 1;
  params->smagorinsky = 0.0;
  params->sue1 = 0.0;
  params->kue1 = 0.0;
  params->bsue1 = 0.0;
  params->bkue1 = 0.0;
  sprintf(params->smag, "%c", '\0');
  params->smag_smooth = 0;
  params->diff_scale = NONLIN;
  params->visc_method = US_LAPLACIAN;
  params->visc_fact = 0.0;
  params->stab = NONE;
  params->atr = 0;
  params->ntr = 0;
  params->ntrS = 0;
  params->nsed = 0;
  params->sednz = 0;
  params->tratio = 1.0;
  params->trsplit = 0;
  sprintf(params->autotrpath, "%c", '\0');
  /* params->sliding = 0; / *UR*/
  params->cfl = NONE;
  memset(params->cfl_dt, 0, sizeof(params->cfl_dt));
  params->mixlayer = NONE;
  params->show_layers = 0;
  params->means = NONE;
  sprintf(params->means_os, "%c", '\0');
  sprintf(params->bathystats, "%c", '\0');
  sprintf(params->particles, "%c", '\0');
  sprintf(params->addquad, "%c", '\0');
  params->vorticity = NONE;
  params->numbers = params->numbers1 = NONE;
  params->u1_f = params->u1av_f = NONE;
  params->save_force = NONE;
  params->lnm = 0.0;
  strcpy(params->trflux, "NONE");
  params->trfd1 = 1;
  params->trfd2 = 4;
  strcpy(params->trperc, "NONE");
  sprintf(params->trpercr, "%c", '\0');
  params->trflsh = 0;
  sprintf(params->trage, "%c", '\0');
  params->ndhw = 0;
  params->tendf = 0;
  sprintf(params->trtend, "%c", '\0');
  params->thin_merge = 1;
  params->sigma = 0;
  params->maxgrad = -1.0;
  sprintf(params->maxdiff, "%c", '\0');
  params->bathyfill = B_REMOVE;
  params->bvel = 0.0;
  params->nonlinear = 1;
  params->calc_dens = 1;
  params->mode2d = 0;
  params->smooth = 0;
  params->slipprm = 1.0;
  params->map_type = 0;
  params->nwindows = 1;
  params->win_reset = 0;
  params->win_type = GROUPED;
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
  params->meshinfo = 1;
  params->swr_type = NONE;
  for (n = 0; n < 6; n++) params->swr_ens[n] = -1.0;
  sprintf(params->u1vhc, "%c", '\0');
  sprintf(params->u1khc, "%c", '\0');
  sprintf(params->restart_name, "%c", '\0');
  sprintf(params->win_file, "%c", '\0');
  sprintf(params->wind_file, "%c", '\0');
  sprintf(params->geom_file, "%c", '\0');
  sprintf(params->bdry_file, "%c", '\0');
  sprintf(params->ptinname, "%c", '\0');
  sprintf(params->wind, "%c", '\0');
  sprintf(params->wind_interp, "%c", '\0');
  sprintf(params->patm, "%c", '\0');
  sprintf(params->patm_interp, "%c", '\0');
  sprintf(params->precip, "%c", '\0');
  sprintf(params->precip_interp, "%c", '\0');
  sprintf(params->evap, "%c", '\0');
  sprintf(params->evap_interp, "%c", '\0');
  sprintf(params->airtemp, "%c", '\0');
  sprintf(params->airtemp_interp, "%c", '\0');
  sprintf(params->rh, "%c", '\0');
  sprintf(params->cloud, "%c", '\0');
  sprintf(params->cloud_interp, "%c", '\0');
  sprintf(params->swr, "%c", '\0');
  sprintf(params->light, "%c", '\0');
  sprintf(params->wetb, "%c", '\0');
  sprintf(params->wetb_interp, "%c", '\0');
  sprintf(params->hftemp, "%c", '\0');
  sprintf(params->hf, "%c", '\0');
  sprintf(params->swr_babs, "%c", '\0');
  sprintf(params->swr_attn, "%c", '\0');
  sprintf(params->swr_attn1, "%c", '\0');
  sprintf(params->swr_tran, "%c", '\0');
  sprintf(params->swr_regions, "%c", '\0');
  params->swreg_dt = 86400.0;
  strcpy(params->swr_data, "GHRSST");
  sprintf(params->densname, "%c", '\0');
  sprintf(params->regions, "%c", '\0');
  sprintf(params->region_dt, "%c", '\0');
  sprintf(params->region_vars, "%c", '\0');
  sprintf(params->region_mode, "%c", '\0');
  sprintf(params->i_rule, "%c", '\0');
  sprintf(params->cookiecut, "%c", '\0');
  sprintf(params->cmacrop, "%c", '\0');
  sprintf(params->cmicrop, "%c", '\0');
  params->region_obc = 0;
  sprintf(params->imp2df, "%c", '\0');
  sprintf(params->imp3df, "%c", '\0');
  params->do_pt = 0;
  params->do_lag = 0;
  strcpy(params->dp_mode, "none");
  params->waves = 0;
  params->decorr = 0;
  params->decf = NONE;
  sprintf(params->monotr, "%c", '\0');
  params->monomn = 0.0;
  params->monomx = 0.0;
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
  params->eta_ib = 0;
  params->etarlx = NONE;
  params->velrlx = NONE;
  params->albedo = -9999.0;
  params->domom = 0;
  params->doroms = 0;
  params->doswan = 0;
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
  params->closf = 0;
  sprintf(params->mixsc, "%c", '\0');
  sprintf(params->s_func, "%c", '\0');
  params->min_tke = 0.0;
  params->min_diss = 0.0;
  params->vz0 = 0.0;
  params->kz0 = 0.0;
  sprintf(params->vz0c, "%c", '\0');
  sprintf(params->kz0c, "%c", '\0');
  params->smooth_VzKz = 0;
  params->albedo_l = -1;
  params->trout = 0;
  params->porusplate = 0;
  params->sharp_pyc = 0;
  params->uscf = 0;
  params->tidef = 0;
  params->tidep = 0;
  params->eqt_alpha = 0.948;
  params->eqt_beta = 0.7;
  params->mlat = params->mlon = NOTVALID;
  sprintf(params->nprof, "%c", '\0');
  sprintf(params->nprof2d, "%c", '\0');
  sprintf(params->reef_frac, "%c", '\0');
  sprintf(params->crf, "%c", '\0');
  params->dbi = params->dbj = params->dbk = -1;
  params->dbgf = NONE;
  params->dbgtime = 0.0;
  memset(params->momfile, 0, sizeof(params->momfile));
  memset(params->avhrr_path, 0, sizeof(params->avhrr_path));
  params->ghrsst_type = G_OSTIA;
  sprintf(params->ghrsst_path, "%c", '\0');
  sprintf(params->ghrsst_opt, "%c", '\0');
  sprintf(params->ghrsst_irule, "%c", '\0');
  strcpy(params->ghrsst_name, "GHRSST L4 SST");
  strcpy(params->ghrsst_dt, "1 day");
  sprintf(params->errornorm, "%c", '\0');
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
  sprintf(params->feta_input_dt, "%c", '\0');
  sprintf(params->feta_interp, "%c", '\0');
  sprintf(params->fsalt_input_dt, "%c", '\0');
  sprintf(params->fsalt_interp, "%c", '\0');
  sprintf(params->ftemp_input_dt, "%c", '\0');
  sprintf(params->ftemp_interp, "%c", '\0');
  sprintf(params->fvelu_input_dt, "%c", '\0');
  sprintf(params->fvelu_interp, "%c", '\0');
  sprintf(params->fvelv_input_dt, "%c", '\0');
  sprintf(params->fvelv_interp, "%c", '\0');
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
  params->nptsu = NULL;
  params->nland = 0;
  params->nturb = 0;

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

  sprintf(params->rendername, "%c", '\0');
  sprintf(params->renderpath, "%c", '\0');
  sprintf(params->renderdesc, "%c", '\0');
  sprintf(params->renderrem, "%c", '\0');
  sprintf(params->rendertype, "%c", '\0');
  sprintf(params->renderopts, "%c", '\0');
  sprintf(params->ecosedconfig, "%c", '\0');
#if defined(HAVE_SEDIMENT_MODULE)
  params->do_sed = 0.0;
  sprintf(params->sed_vars, "%c", '\0');
  sprintf(params->sed_defs, "%c", '\0');
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  params->do_eco = 0;
  sprintf(params->eco_vars, "%c", '\0');
  sprintf(params->eco_defs, "%c", '\0');
#endif

#if defined(HAVE_WAVE_MODULE)
  params->do_wave = NONE;
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
  int cc;
  double d1;
  struct timeval tm1;
  int rconf = 0;

  /* Allocate memory for the parameter data structure                */
  params = params_alloc();

  /* Read in the name of the input file                              */
  params->prmfd = fp;
  prm_set_errfn(hd_quit);
  gettimeofday(&tm1, NULL);
  params->initt = tm1.tv_sec + tm1.tv_usec * 1e-6;

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

  /* Near real-time definitions                                      */
  if (nrt_restart) {
    /* Check the model TIMEUNIT is seconds                           */
    if (contains_token(schedule->units, "seconds") == NULL)
      hd_quit("get_nrt_time: TIMEUNIT must be 'seconds since ....'.\n");
    params->t = get_restart_time(params->idumpname, schedule->units);
    sprintf(params->start_time, "%g days", params->t / 86400.0);
    params->rampstart = params->rampend = params->t;
    tm_scale_to_secs(params->stop_time, &d1);
    sprintf(params->stop_time, "%g days", (d1 + params->t) / 86400.0);
  }

  /* Input file type                                                 */
  params->us_type = (US_RUS|US_WUS|US_P);
  sprintf(keyword, "INPUT_FILE_TYPE");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->us_type = (US_RUS|US_WUS|US_P);
    } else {
      if (contains_token(buf, "STRUCTURED") != NULL)
	params->us_type |= (US_RS|US_WS|US_P);
      if (contains_token(buf, "UNSTRUCTURED") != NULL)
	params->us_type |= (US_RUS|US_WUS|US_P);
    }
  }
  if (prm_read_char(fp, "POWER_MESH", buf)) {
    if (is_true(buf)) params->us_type |= US_POW;
  }
  if (prm_read_char(fp, "STEREOGRAPHIC_MESH", buf)) {
    if (is_true(buf)) params->us_type |= US_STER;
  }
  sprintf(keyword, "INTERP_RULE");
  if (prm_read_char(fp, keyword, params->i_rule));

  /* Grid dimensions                                                 */
  prm_read_int(fp, "nMaxMesh2_face_nodes", &params->npe);
  prm_read_int(fp, "nMesh2_face", &params->ns2);
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

  if (params->us_type & US_RS)
    read_bathy_from_file(fp,params);

  sprintf(keyword, "GRIDTYPE");
  prm_read_char(fp, keyword, params->gridtype);

  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "DESCRIPTION");
  prm_read_char(fp, keyword, params->grid_desc);
  sprintf(keyword, "NAME");
  prm_read_char(fp, keyword, params->grid_name);
  /* gridhoriz(params,geom,fp); */

  /* Flags                                                           */
  /* Set defaults                                                    */
  params->runmode = MANUAL;
  set_default_param(params);
  if (prm_read_char(fp, "HYDRO_CONFIG", buf)) {
    if (get_rendered_hydroparams(buf, params))
      hd_warn("Can't find hydrodynamic configuration '%s'\n", buf);
    else {
      hd_warn("Hydrodynamic parameters configured to '%s.h'\n", buf);
      rconf = 1;
    }
  }

  sprintf(keyword, "MESH_INFO");
  if (prm_read_char(fp, keyword, buf))
    params->meshinfo = is_true(buf);
  if (prm_read_char(fp, "REORDER_WRITE", params->mesh_reorder))
    params->mrf |= MR_WRITE;
  if (prm_read_char(fp, "REORDER_WRITEX", params->mesh_reorder))
    params->mrf |= MR_WRITEX;
  if (prm_read_char(fp, "REORDER_READ", params->mesh_reorder))
    params->mrf |= MR_READ;
  sprintf(keyword, "ETAMAX");
  prm_read_double(fp, keyword, &params->etamax);
  sprintf(keyword, "MIN_CELL_THICKNESS");
  prm_read_char(fp, keyword, params->mct);
  read_history(params, fp);

  /* Bathymetry statistics                                           */
  sprintf(keyword, "BATHY_STATS");
  if (prm_read_char(fp, keyword, params->bathystats))
    params->ntrS += 4;

  /* User regulation                                                 */
  sprintf(keyword, "REGULATE_FILE");
  prm_read_char(fp, keyword, params->regulate);
  if (prm_read_char(fp, "REGULATE_DT", buf))
    tm_scale_to_secs(buf, &params->regulate_dt);
  sprintf(keyword, "CRASH_RECOVERY");
  prm_read_char(fp, keyword, params->crf);

  /* ID number and revision                                          */
  sprintf(keyword, "ID_NUMBER");
  prm_read_char(fp, keyword, params->runnoc);
  params->runno = atof(params->runnoc);
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

  /*-----------------------------------------------------------------*/
  /* Advection schemes                                               */
  read_advect(params, fp);

  /* Setting robust > 1 will cap and smooth smagorinsky              */
  sprintf(keyword, "ROBUST");
  prm_read_int(fp, keyword, &params->robust);

  /* Stability sub-stepping in advection schemes                     */
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

  /* Thin layer merging                                              */
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

  /* Variable filtering                                              */
  sprintf(keyword, "FILTERING");
  if (prm_read_char(fp, keyword, buf)) {
    params->filter = 0;
    if (strcmp(buf, "NONE") == 0)
      params->filter = NONE;
    if (strcmp(buf, "ADVECT") == 0)
      params->filter |= ADVECT;
  }
  read_trfilter(params, fp);

  /* Porus plate sub-gridscale parameterisation                      */
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
  sprintf(keyword, "TIDE_POTENTIAL");
  if (prm_read_char(fp, keyword, buf)) {
    if (params->tidep = is_true(buf)) {
      prm_read_double(fp, "EQT_ALPHA", &params->eqt_alpha);
      prm_read_double(fp, "EQT_BETA", &params->eqt_beta);
      params->ntrS += 1;
      params->tidef |= TD_EQT;
    }
  }
  read_means(params, fp, 0);

  /* Window sizes                                                    */
  read_window_info(params, fp);

  sprintf(keyword, "PARRAY_INTPL_BOTZ");
  if(prm_read_char(fp, keyword, buf))
    params->parray_inter_botz = is_true(buf);

  /* Compatibility                                                   */
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
    params->ntrS += 6;

  params->vorticity = 0;
  sprintf(keyword, "VORTICITY");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->vorticity = NONE;
    }
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

  read_trflux(params, fp);
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
  read_profile(params, fp);

  read_decorr(params, fp, 0);
  read_monotone(params, fp, 0);

  /* Totals diagnostics                                              */
  read_totals(params, fp);

  /* Regions                                                         */
  read_region_info(params, fp);

  /* AVHRR SST */
  if (prm_read_char(fp, "AVHRR", params->avhrr_path)) {
    not_included("AVHRR import");
    hd_warn("Use GHRSST!");
    /*
    create_avhrr_list(params);
    params->avhrr = 1;
    params->ntrS++;
    */
  }

  /* GHRSST SST                                                      */
  if (prm_read_char(fp, "GHRSST", params->ghrsst_path)) {
    prm_read_char(fp, "GHRSST_OPTIONS", params->ghrsst_opt);
    prm_read_char(fp, "GHRSST_INTERP", params->ghrsst_irule);
    create_ghrsst_list(params);
    params->ntrS+=2;
  }

  /* Error norm diagnostic */
  if (prm_read_char(fp, "ERROR_NORM", params->errornorm)) {
    if (prm_read_char(fp, "ERROR_NORM_DT", buf))
      tm_scale_to_secs(buf, &params->enorm_dt);
    else
      params->enorm_dt = 3600.0;
  }

  sprintf(keyword, "MACRO_PLASTIC");
  prm_read_char(fp, keyword, params->cmacrop);
  sprintf(keyword, "MICRO_PLASTIC");
  prm_read_char(fp, keyword, params->cmicrop);

  /* Particle tracking                                               */
  if(prm_read_char(fp, "PT_InputFile", params->ptinname)) {
    params->do_pt = PT_DO;
    params->ntr++;
    /* Auto particle source                                          */
    if (prm_read_char(fp, "particles", params->particles))
      params->do_pt = PT_AUTO;
    params->do_lag = params->do_pt;
  }

  /* Transport mode files                                            */
  read_tmode_params(fp, params);

  /* Data assimilation                                               */
  params->da = 0;
  sprintf(keyword, "DATA_ASSIMILATION");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->da = is_true(buf))) {
      /* Global dt and err for obs                                   */
      double obs_dt  = -1;
      double obs_err = -1;
#ifndef HAVE_DA
      hd_quit("DATA ASSIMILATION has been invoked in the prm file but this version of COMPAS is not built with DA support. Please either disable DA or reconfigure with the --enable-da option and rebuild");
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

      /* Always convert DT to seconds                                */
      prm_get_time_in_secs(fp, "DA_OBS_DT", &obs_dt);
      prm_read_double(fp, "DA_OBS_ERR", &obs_err);
      
      /* See if we're using localisation                             */
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
	  /* Name of this observation                                */
	  sprintf(keyword, "DA_OBS%1.1d.name", n);
	  params->da_obs_names[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_names[n] = strdup(buf);
	  /* File name for the data                                  */
	  sprintf(keyword, "DA_OBS%1.1d.data", n);
	  params->da_obs_files[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_files[n] = strdup(buf);
	  /* What the platform type is                               */
	  sprintf(keyword, "DA_OBS%1.1d.type", n);
	  params->da_obs_types[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_types[n] = strdup(buf);
	  /* The location of this observation                        */
	  sprintf(keyword, "DA_OBS%1.1d.location", n);
	  params->da_obs_locs[n] = NULL;
	  if (prm_read_char(fp, keyword, buf))
	    params->da_obs_locs[n] = strdup(buf);
	  /* Any DT  associated with this                            */
	  sprintf(keyword, "DA_OBS%1.1d.dt", n);
	  params->da_obs_dt[n] = obs_dt;
	  prm_read_double(fp, keyword, &params->da_obs_dt[n]);
	  /* Any errors associated with this                         */
	  sprintf(keyword, "DA_OBS%1.1d.error", n);
	  params->da_obs_errs[n] = obs_err;
	  prm_read_double(fp, keyword, &params->da_obs_errs[n]);
	}
      }
    }
  }

  /* Diagnistic numbers                                              */
  params->ntr += numbers_init(params);
  params->ntr += import_init(params, fp);

  /* Degree heating days                                             */
  params->ntr += read_dhw(params, fp);

  /* Auto point source                                               */
  if (prm_read_char(fp, "pss", buf)) {
    params->numbers |= PASS;
    params->ntr += 1;
  }

  /* River flow diagnostic                                           */
  if (check_river_bdry(params, fp)) {
    params->riverflow = 1;
    params->ntrS++;
    if (check_bdry_options(params, fp) & OP_DYNAHC) {
      params->riverflow = 2;
      params->ntrS++;
      params->ntr++;
    }
  }

  /* Processes to omit                                               */
  params->u1_f = params->u1av_f = 0;
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

  /* Runtime alerts and constraint                                   */
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

  /* Ramping                                                         */
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

  /* Variable maximum values                                         */
  sprintf(keyword, "VELMAX");
  prm_read_double(fp, keyword, &params->velmax);
  sprintf(keyword, "VELMAX_2D");
  if(!params->velmax2d && !prm_read_double(fp, keyword, &params->velmax2d))
    params->velmax2d = params->velmax;
  sprintf(keyword, "ETA_DIFF");
  prm_read_double(fp, keyword, &params->etadiff);
  sprintf(keyword, "WMAX");
  prm_read_double(fp, keyword, &params->wmax);

  /* General mandatory parameters                                    */
  prm_set_errfn(quit);
  sprintf(keyword, "HMIN");
  prm_read_double(fp, keyword, &params->hmin);
  sprintf(keyword, "UF");
  prm_read_double(fp, keyword, &params->uf);
  sprintf(keyword, "QBFC");
  prm_read_char(fp, keyword, params->quad_bfc);
  sprintf(keyword, "RAMPSTART");
  prm_get_time_in_secs(fp, keyword, &params->rampstart);
  sprintf(keyword, "RAMPEND");
  prm_get_time_in_secs(fp, keyword, &params->rampend);
  if (params->rampend < params->rampstart)
    hd_quit("readparam: ramp end time before ramp start time\n");
  if(params->rampstart > params->t)
    hd_warn("readparam: ramp start time after model start time\n");

  /* Time                                                            */
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

  /* Constants                                                       */
  prm_set_errfn(hd_quit);
  strcpy(params->prmname, prmname);
  sprintf(keyword, "CODEHEADER");
  prm_read_char(fp, keyword, params->codeheader);
  strcpy(codeheader, params->codeheader);
  sprintf(keyword, "PARAMETERHEADER");
  prm_read_char(fp, keyword, params->parameterheader);
  strcpy(parameterheader, params->parameterheader);
  sprintf(keyword, "REFERENCE");
  prm_read_char(fp, keyword, params->reference);
  sprintf(keyword, "NOTES");
  prm_read_char(fp, keyword, params->notes);
  sprintf(keyword, "TECHNOLOGY_READINESS_LEVEL");
  prm_read_char(fp, keyword, params->trl);
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

  /* Horizontal diffusion                                            */
  read_hdiff(params, fp, 0);

  /* Bottom roughness                                                */
  z0_init(params);

  /* Vertical mixing schemes                                         */
  prm_set_errfn(hd_quit);
  sprintf(keyword, "MIXING_SCHEME");
  prm_read_char(fp, keyword, params->mixsc);
  if (strcmp(params->mixsc, "constant") == 0 &&
      params->mixlayer == TKE_MIX) {
    params->mixlayer = NONE;
    hd_warn("MIX_LAYER = TKE_MIX not functional with MIXING_SCHEME = constant.\n");
  }
  read_vdiff(params, fp, 0);
  /*
  sprintf(keyword, "KZ0");
  prm_read_double(fp, keyword, &params->kz0);
  sprintf(keyword, "VZ0");
  prm_read_double(fp, keyword, &params->vz0);
  sprintf(keyword, "KZ0");
  prm_read_char(fp, keyword, params->kz0c);
  sprintf(keyword, "VZ0");
  prm_read_char(fp, keyword, params->vz0c);
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "KZ_ALPHA");
  prm_read_double(fp, keyword, &params->kz_alpha);
  sprintf(keyword, "VZ_ALPHA");
  prm_read_double(fp, keyword, &params->vz_alpha);
  */
  sprintf(keyword, "ZS");
  prm_read_double(fp, keyword, &params->zs);
  if (params->zs == 0.0) {
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
  if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
      strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
    params->ntr += 1;
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0)
    params->ntr += 4;
  if (strcmp(params->mixsc, "harcourt") == 0)
    params->ntr += 4;

  /* Coriolis                                                        */
  prm_set_errfn(hd_silent_warn);
  params->coriolis = NULL;
  sprintf(keyword, "CORIOLIS");
  prm_read_darray(fp, keyword, &params->coriolis, &params->nvals);
  if (prm_read_darray(fp, keyword, &params->coriolis, &params->nvals))
    if (params->nce1 && params->nce2 && params->nvals != params->nce1 * params->nce2)
      hd_quit(" Number of values for CORIOLIS incorrect: %d != %d\n",
	      params->nvals, params->nce1 * params->nce2);

  /* Surface height                                                  */
  sprintf(params->eta_init, "%c", '\0');
  sprintf(keyword, "SURFACE");
  prm_read_char(fp, keyword, params->eta_init);
  /* Inverse barometer initial condition only for -g option
  sprintf(keyword, "SURFACE_INV_BAR");
  if (prm_read_char(fp, keyword, buf))
    params->eta_ib = is_true(buf);
  */
  sprintf(params->vel_init, "%c", '\0');
  sprintf(keyword, "VELOCITY");
  prm_read_char(fp, keyword, params->vel_init);

  /* Cells explicitly defined as OUTSIDE                             */
  params->polyland = (char **)malloc(MAXNUMVARS * sizeof(char *));
  for (n = 0; n < MAXNUMVARS; n++) {
    params->polyland[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    sprintf(params->polyland[n], "%c", '\0');
  }
  read_blocks(fp, "NLAND", &params->nland, &params->lande1, &params->lande2, params->polyland);

  /* Forcing data */
  params->save_force = 0;
  if(prm_read_char(fp, "ETA_DATA", params->edata)) {
    if (!(prm_read_char(fp, "ETA_INPUT_DT", params->feta_input_dt)))
      strcpy(params->feta_input_dt, "1 day");
    prm_read_char(fp, "ETA_IRULE", params->feta_interp);
    params->save_force |= FETA;
    params->ntrS += 1;
  }
  if(prm_read_char(fp, "SALT_DATA", params->sdata)) {
    if (!(prm_read_char(fp, "TEMP_INPUT_DT", params->fsalt_input_dt)))
      strcpy(params->fsalt_input_dt, "1 day");
    prm_read_char(fp, "TEMP_IRULE", params->fsalt_interp);
    params->save_force |= FSALT;
    params->ntr += 1;
  }
  if(prm_read_char(fp, "TEMP_DATA", params->tdata)) {
    if (!(prm_read_char(fp, "TEMP_INPUT_DT", params->ftemp_input_dt)))
      strcpy(params->ftemp_input_dt, "1 day");
    prm_read_char(fp, "TEMP_IRULE", params->ftemp_interp);
    params->save_force |= FTEMP;
    params->ntr += 1;
  }
  if(prm_read_char(fp, "VELOCITY_DATA", params->vdata)) {
    if (prm_read_char(fp, "VELOCITY_INPUT_DT", params->fvelu_input_dt))
      strcpy(params->fvelv_input_dt, params->fvelu_input_dt);
    else {
      strcpy(params->fvelu_input_dt, "1 day");
      strcpy(params->fvelv_input_dt, "1 day");
    }
    if (prm_read_char(fp, "VELOCITY_IRULE", params->fvelu_interp))
      strcpy(params->fvelv_interp, params->fvelu_interp);
    params->save_force |= (FVELU|FVELV);
    params->ntr += 2;
  }
  if (params->save_force == 0) params->save_force = NONE;

  /* Surface relaxation                                              */
  read_eta_relax(params, fp);
  read_vel_relax(params, fp);

  /* Wind                                                            */
  sprintf(params->wind, "%c", '\0');  sprintf(params->wind, "%c", '\0');  sprintf(params->wind, "%c", '\0');
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
    prm_read_char(fp, "WIND_IRULE", params->wind_interp);
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

  /* Atmospherics                                                    */
  prm_set_errfn(hd_silent_warn);
  sprintf(params->patm, "%c", '\0');
  sprintf(keyword, "PRESSURE");
  prm_read_char(fp, keyword, params->patm);
  sprintf(keyword, "PRESSURE_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->patm_dt);
  sprintf(keyword, "PRESSURE_IRULE");
  prm_read_char(fp, keyword, params->patm_interp);
  sprintf(params->precip, "%c", '\0');
  sprintf(keyword, "PRECIPITATION");
  prm_read_char(fp, keyword, params->precip);
  sprintf(keyword, "PRECIPITATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->precip_dt);
  sprintf(keyword, "PRECIPITATION_IRULE");
  prm_read_char(fp, keyword, params->precip_interp);
  sprintf(params->evap, "%c", '\0');
  sprintf(keyword, "EVAPORATION");
  prm_read_char(fp, keyword, params->evap);
  sprintf(keyword, "EVAPORATION_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->evap_dt);
  sprintf(keyword, "EVAPORATION_IRULE");
  prm_read_char(fp, keyword, params->evap_interp);
  sprintf(params->airtemp, "%c", '\0');
  sprintf(keyword, "AIRTEMP");
  prm_read_char(fp, keyword, params->airtemp);
  sprintf(keyword, "AIRTEMP_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->airtemp_dt);
  sprintf(keyword, "AIRTEMP_IRULE");
  prm_read_char(fp, keyword, params->airtemp_interp);
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
  sprintf(keyword, "CLOUD_IRULE");
  prm_read_char(fp, keyword, params->cloud_interp);
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
  sprintf(params->lwri, "%c", '\0');
  sprintf(keyword, "LONGWAVE_IN");
  prm_read_char(fp, keyword, params->lwri);
  sprintf(keyword, "LONGWAVE_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->lwri_dt);
  sprintf(params->wetb, "%c", '\0');
  sprintf(keyword, "WET_BULB");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "WET_BULB_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);
  sprintf(keyword, "WET_BULB_IRULE");
  prm_read_char(fp, keyword, params->wetb_interp);
  sprintf(keyword, "DEW_POINT");
  prm_read_char(fp, keyword, params->wetb);
  sprintf(keyword, "DEW_POINT_INPUT_DT");
  prm_get_time_in_secs(fp, keyword, &params->wetb_dt);
  sprintf(keyword, "DEW_POINT_IRULE");
  prm_read_char(fp, keyword, params->wetb_interp);

  /* Heatflux                                                        */
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

    if (contains_token(buf, "ADVANCED") != NULL) params->heatflux = ADVANCED;
    if (contains_token(buf, "BULK") != NULL) params->heatflux = ADVANCED;
    if (params->heatflux & ADVANCED) {
      /*
    if (strcmp(buf, "ADVANCED") == 0 || strcmp(buf, "BULK") == 0 || 
	params->heatflux & ADVANCED) {
      params->heatflux = ADVANCED;
    */

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

    if (contains_token(buf, "GHRSST") != NULL) {
      params->heatflux |= GHRSST;
      if (!strlen(params->ghrsst_path)) {
	hd_warn
	  ("params_read() : GHRSST heatflux requires a GHRSST diagnostic path name.\n");
	params->heatflux = NONE;
      } else {
	sprintf(keyword, "HEATFLUX_TC");
	if (!(prm_get_time_in_secs(fp, keyword, &params->hftc))) {
	  hd_warn
	    ("params_read() : GHRSST heatflux requires HEATFLUX_TC constant.\n");
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

  /* Saltflux                                                        */
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
      if (!(params->heatflux & ADVANCED))
	hd_warn("params_read(): BULK saltflux must also use BULK heatflux. No saltflux included.\n");
      else {
	params->saltflux = BULK;
	params->ntrS += 1;
      }
    }
  }

  /* Tidal energy extraction                                         */
  read_turb(params, fp);

  prm_set_errfn(hd_silent_warn);

  /* Sediment geometry                                               */
  params->sednz = 0;

  /* Note : sediment layers are read in assuming the first layer is  */
  /* closest to the water column. This is reversed when copying to   */
  /* the master so that the sediment origin is the deepest layer.    */
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

  sprintf(keyword, "RENDER_NAME");
  if (prm_read_char(fp, keyword, params->rendername)) {
    sprintf(keyword, "RENDER_DESC");
    prm_read_char(fp, keyword, params->renderdesc);
  }
  sprintf(keyword, "RENDER_PATH");
  prm_read_char(fp, keyword, params->renderpath);
  sprintf(keyword, "RENDER_TYPE");
  prm_read_char(fp, keyword, params->rendertype);
  sprintf(keyword, "RENDER_OPTIONS");
  prm_read_char(fp, keyword, params->renderopts);
  sprintf(keyword, "RENDER_REMOVE");
  prm_read_char(fp, keyword, params->renderrem);
  sprintf(keyword, "ECOSED_CONFIG");
  prm_read_char(fp, keyword, params->ecosedconfig);
#if defined(HAVE_SEDIMENT_MODULE)
  /* Parameters for coupling to sediment module                      */
  read_sediments(params, fp, &params->ntr);
#endif

#if defined(HAVE_WAVE_MODULE)
  /* Parameters for coupling to wave module                          */
  params->do_wave = NONE;
  sprintf(keyword, "DO_WAVES");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->do_wave = NONE;
    if (strcmp(buf, "FILE") == 0)
      params->do_wave = W_FILE;
    if (strcmp(buf, "COMP") == 0)
      params->do_wave = W_COMP;
    if (strcmp(buf, "SWAN") == 0)
      params->do_wave = (W_SWAN|W_SWANM);
    if (strcmp(buf, "SWAN_W") == 0)
      params->do_wave = (W_SWAN|W_SWANW);
    prm_get_time_in_secs(fp, "WAVES_DT", &params->wavedt);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  /* Parameters for coupling to biogeochemical                       */
  read_ecology(params, fp, &params->ntr);
#endif

  /* Wave module parameters                                          */
  read_waves(params, fp, 0);
  if (params->waves & STOKES_DRIFT && params->tendf) {
    params->ntr+=2;
  }

  /* Library error handling                                          */
  sprintf(keyword, "LIB_ERROR_FCN");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "QUIT") == 0)
      params->gint_errfcn = LFATAL;
    if (strcmp(buf, "WARN") == 0)
      params->gint_errfcn = LWARN;
  }

  /*-----------------------------------------------------------------*/
  /* Tracers                                                         */
  /* 3D Tracer constants and variables.                              */
  /* Tracers relating to certain diagnostics etc. are automatically  */
  /* generated if the relevant flag is set. If these tracers are     */
  /* manually defined in the input parameter file, then the manual   */
  /* definition over-rides the auto generation, and the defined      */
  /* tracer_info is used.                                            */
  /* Always assume salt and temp are auto tracers.                   */
  params->ntr += 2;
  params->atr = params->ntr;
  tracer_setup(params, fp);
  create_tracer_3d(params);
  for (n = 0; n < params->atr; n++) {
    params->trinfo_3d[n].m = -1;
    params->trinfo_3d[n].n = n;
  }
  sprintf(keyword, "AUTOTRACERPATH");
  prm_read_char(fp, keyword, params->autotrpath);

  /* Get the tracers which undergo diffusion (horizontal and         */
  /* vertical)  and store in tdif_h[] and tdif_v[].                  */
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

  /* Reset the arrays for split tracer transport                     */
  if (params->trsplit) {
    int nvec, *vec;
    /* Horizontal diffusion                                          */
    nvec = params->ntdif_h;
    vec = i_alloc_1d(nvec);
    memcpy(vec, params->tdif_h, nvec * sizeof(int));
    if (params->ntdif_h) {
      params->ntdif_hs = params->ntdif_h - 2;
      if (params->ntdif_hs) {
	params->tdif_hs = i_alloc_1d(params->ntdif_hs);
	params->ntdif_h = params->ntdif_hs = 0;
	for (n = 0; n < nvec; n++) {
	  tracer_info_t *tracer;
	  int tn;
	  tn = vec[n];
	  tracer = &params->trinfo_3d[tn];
	  if (strcmp(tracer->name, "salt") == 0)
	    params->tdif_h[params->ntdif_h++] = tn;
	  else if (strcmp(tracer->name, "temp") == 0)
	    params->tdif_h[params->ntdif_h++] = tn;
	  else
	    params->tdif_hs[params->ntdif_hs++] = tn;
	}
      }
    }
    i_free_1d(vec);
    /* Vertical diffusion                                            */
    nvec = params->ntdif_v;
    vec = i_alloc_1d(nvec);
    memcpy(vec, params->tdif_v, nvec * sizeof(int));
    if (params->ntdif_v) {
      params->ntdif_vs = params->ntdif_v - 2;
      if (params->ntdif_vs) {
	params->tdif_vs = i_alloc_1d(params->ntdif_vs);
	params->ntdif_v = params->ntdif_vs = 0;
	for (n = 0; n < nvec; n++) {
	  tracer_info_t *tracer;
	  int tn;
	  tn = vec[n];
	  tracer = &params->trinfo_3d[tn];
	  if (strcmp(tracer->name, "salt") == 0)
	    params->tdif_v[params->ntdif_v++] = tn;
	  else if (strcmp(tracer->name, "temp") == 0)
	    params->tdif_v[params->ntdif_v++] = tn;
	  else
	    params->tdif_vs[params->ntdif_vs++] = tn;
	}
      }
    }
    i_free_1d(vec);
  }

  /* Sediment tracer constants and variables                         */
  prm_set_errfn(hd_silent_warn);
  tracer_read(fp, NULL, SEDIM, hd_quit, hd_warn, hd_silent_warn,
        &params->nsed, &params->trinfo_sed);
  if (params->nsed && !params->sednz)
    hd_quit("Sediment tracers requires NSEDLAYERS > 0\n");

  /* 2D Tracer constants and variables                               */
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

  /* Timeseries file caching                                         */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "CACHE_TSFILES");
  if (prm_read_char(fp, keyword, buf) > 0)
    params->tsfile_caching = is_true(buf);
  else
    params->tsfile_caching = 1;
  emstag(LDEBUG,"hd:readparam:params_read","Setting ts_file_cachinf to: %s",(params->tsfile_caching?"true":"false"));

  /* Explicit mappings                                               */
  read_explicit_maps(params, fp);

  /* Process exclusion                                               */
  read_exclude_points(params, fp);

  /*-----------------------------------------------------------------*/
  /* Open boudaries                                                  */
  get_bdry_params(params, fp);

  /* Allocate memory for the open boundary data structure            */
  if (params->open == NULL) {
    params->open =
      (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);
  }

  for (n = 0; n < params->nobc; n++) {

    /* Allocate memory for the boundary structures                   */
    if (!rconf) params->open[n] = OBC_alloc();

    /* Read the open boudary conditions                              */
    params->open[n]->ntr = params->ntr;
    params->open[n]->atr = params->atr;
    get_OBC_conds(params, params->open[n], fp, n, params->trinfo_3d);
    get_OBC_relax(params, params->open[n], fp, n);
    get_obc_list(params->open[n], fp, n, "BOUNDARY");
    prm_set_errfn(hd_quit);
  }

  /* Read the csr tide model paths if required                       */
  sprintf(keyword, "TIDE_CSR_ORTHOWEIGHTS");
  prm_read_char(fp, keyword, params->orthoweights);
  sprintf(keyword, "TIDE_CSR_CON_DIR");
  prm_read_char(fp, keyword, params->nodal_dir);
  sprintf(keyword, "TIDE_CONSTITUENTS");
  prm_read_char(fp, keyword, params->tide_con_file);
  sprintf(keyword, "TIDE_CONSTITUENTS_TRAN");
  if (prm_read_char(fp, keyword, params->tide_con_file)) 
    params->tidef |= TD_TRAN;
  sprintf(keyword, "TIDE_CONSTITUENTS_VEL");
  if (prm_read_char(fp, keyword, params->tide_con_file)) 
    params->tidef |= TD_VEL;
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

  /* Whether we read do threaded I/O on some forcing                 */
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

  /* Output files to write to setup.txt                              */
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

  sprintf(keyword, "REREAD_INPUT");
  if(prm_read_char(fp, keyword, buf1))
    not_included("REREAD_INPUT");

  if (DEBUG("init_m"))
    dlog("init_m", "Input parameters read OK\n");
  return (params);
}

/* END params_read()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set input parameters to standard values                */
/*-------------------------------------------------------------------*/
parameters_t *auto_params(FILE * fp, int autof)
{
  double **bathy = NULL;        /* Bathymetry array                  */
  int i, j, m, n, c, nz;        /* Counters                          */
  int ii, jj, tn;               /* Counters                          */
  int nr = 0, rs = 0;           /* River counters                    */
  int bdry = 0;                 /* Tracks open boundaries            */
  open_bdrys_t *open;           /* Pointer to boundary structure     */
  char keyword[MAXSTRLEN];      /* Input dummy                       */
  char buf[MAXSTRLEN];          /* Input dummy                       */
  size_t sz;                    /* Size of tracer structure          */
  double d1;                    /* Dummy                             */
  int dataf;                    /* Flag for boundary data            */
  int is, ie, js, je;           /* Limits of structured grid         */
  double bsgn;                  /* Bathymetry sign                   */
  int *ib = NULL, *jb = NULL, *rb = NULL, *rf = NULL;
  int **jci = NULL, *masku = NULL;
  int *nib = NULL, **mask = NULL; 
  int ic = 0, jc = 0;
  double *rhc;
  struct timeval tm1;
  FILE* fptr = fp;

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the parameter data structure                */
  params = params_alloc();
  gettimeofday(&tm1, NULL);
  params->initt = tm1.tv_sec + tm1.tv_usec * 1e-6;

  /* Read in the name of the input file                              */
  params->prmfd = fp;
  switch (autof) {
  case 1: /* -a, auto mode                                           */
    params->runmode = AUTO;
    break;

  case 2: /* -g, input file generation                               */
    params->runmode = DUMP | EXIT;
    break;

  case 3: /* -ag, auto mode with initialisation termination          */
    params->runmode = AUTO | EXIT;
    break;

  case 4: /* -r, ROAM mode                                           */
    params->runmode = AUTO | ROAM;
    break;

  case 5: /* -rg, ROAM mode with initialisation termination          */
    params->runmode = AUTO | ROAM | EXIT;
    break;

  case 6: /* -rs, ROAM mode with initialisation from previous run    */
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
  /* Read mandatory parameters                                       */
  set_default_param(params);
  prm_set_errfn(hd_quit);
  sprintf(keyword, "START_TIME");
  prm_read_char(fp, keyword, params->start_time);
  sprintf(keyword, "STOP_TIME");
  prm_read_char(fp, keyword, params->stop_time);
  prm_get_time_in_secs(fp, "START_TIME", &params->t);
  prm_read_int(fp, "nMaxMesh2_face_nodes", &params->npe);
  prm_read_int(fp, "nMesh2_face", &params->ns2);
  sprintf(keyword, "NCE1");
  prm_read_int(fp, keyword, &nce1);
  sprintf(keyword, "NCE2");
  prm_read_int(fp, keyword, &nce2);
  params->nce1 = nce1;
  params->nce2 = nce2;
  if (prm_read_char(fp, "REORDER_WRITE", params->mesh_reorder))
    params->mrf |= MR_WRITE;
  if (prm_read_char(fp, "REORDER_WRITEX", params->mesh_reorder))
    params->mrf |= MR_WRITEX;
  if (prm_read_char(fp, "REORDER_READ", params->mesh_reorder))
    params->mrf |= MR_READ;
  prm_set_errfn(hd_silent_warn);

  /* Add a quad grid to the structured mesh                          */
  sprintf(keyword, "ADD_QUAD");
  prm_read_char(fp, keyword, params->addquad);

  /* Porus plate sub-gridscale parameterisation                      */
  sprintf(keyword, "PORUS_PLATE");
  if (prm_read_char(fp, keyword, params->reef_frac)) {
    params->porusplate = 1;
    params->atr += 2;
  }

  /* Window sizes                                                    */
  read_window_info(params, fp);

  /* Compatibility                                                   */
  read_compatible(params, fp);

  /* Debugging                                                       */
  read_debug(params, fp);

  /* Used for writing transport files                                */
  sprintf(keyword, "INPUT_FILE");
  prm_read_char(fp, keyword, params->idumpname);
  prm_read_char(fp, "TRANS_DATA", params->trans_data);

  /* Time unit (optional). Note; this is read in as a mandatory      */
  /* parameter in sched_init().                                      */
  sprintf(keyword, "TIMEUNIT");
  if (!prm_read_char(fp, keyword, params->timeunit))
    strcpy(params->timeunit, "seconds since 1990-01-01 00:00:00 +08");

  /* ID number and revision                                          */
  sprintf(keyword, "ID_NUMBER");
  prm_read_char(fp, keyword, params->runnoc);
  params->runno = atof(params->runnoc);
  sprintf(keyword, "REVISION");
  prm_read_char(fp, keyword, params->rev);

  /* Input mesh types                                                */
  sprintf(keyword, "INPUT_FILE_TYPE");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->us_type = (US_RUS|US_WUS);
    } else {
      if (contains_token(buf, "STRUCTURED") != NULL)
	params->us_type |= (US_RS|US_WS|US_G);
      if (contains_token(buf, "UNSTRUCTURED") != NULL)
	params->us_type |= (US_RUS|US_WUS);
    }
  }
  if (prm_read_char(fp, "POWER_MESH", buf)) {
    if (is_true(buf)) params->us_type |= US_POW;
  }
  if (prm_read_char(fp, "STEREOGRAPHIC_MESH", buf)) {
    if (is_true(buf)) params->us_type |= US_STER;
  }
  sprintf(keyword, "CONVERSION");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->uscf = 0;
    } else {
      if (contains_token(buf, "TRI") != NULL)
	params->uscf = US_TRI;
      if (contains_token(buf, "QUAD") != NULL)
	params->uscf = US_QUAD;
      if (contains_token(buf, "HEX") != NULL)
	params->uscf = US_HEX;
    }
  }
  sprintf(keyword, "INTERP_RULE");
  if (prm_read_char(fp, keyword, params->i_rule));
  prm_read_char(fp, "COOKIE_CUT", params->cookiecut);

  /* Advection (optional)                                            */
  read_advect(params, fp);

  /* Output files (optional)                                         */
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
  /* Set the flags to default values                                 */
  params->grid_dt = 0.0;
  params->iratio = 0;
  params->stab = SUB_STEP_NOSURF;
  params->uf = 1e-4;
  strcpy(params->quad_bfc, "0.003");
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
  ts_set_default_hashtable_size((params->nce1 + 1) * (params->nce2 + 1) * 4);
  /* Headers for the input dump                                      */
  strcpy(codeheader, "COMPAS default version");
  strcpy(params->prmname, prmname);
  strcpy(params->codeheader, codeheader);
  sprintf(keyword, "NAME");
  if (!prm_read_char(fp, keyword, params->grid_name))
    strcpy(params->grid_name, "COMPAS_grid");
  sprintf(params->grid_desc, "Automated grid from %s", prmname);
  /* Parameter header (optional)                                     */
  sprintf(keyword, "PARAMETERHEADER");
  if (!(prm_read_char(fp, keyword, params->parameterheader)))
    strcpy(params->parameterheader, "Auto grid");
  strcpy(parameterheader, params->parameterheader);
  sprintf(keyword, "REFERENCE");
  prm_read_char(fp, keyword, params->reference);
  sprintf(keyword, "NOTES");
  prm_read_char(fp, keyword, params->notes);
  sprintf(keyword, "TECHNOLOGY_READINESS_LEVEL");
  prm_read_char(fp, keyword, params->trl);
  prm_read_char(fp, "restart_name", params->restart_name);

  /* Mixing scheme (optional)                                        */
  sprintf(keyword, "MIXING_SCHEME");
  if (prm_read_char(fp, keyword, params->mixsc)) {
    params->min_tke = 0.0;
    params->min_diss = 0.0;
    params->vz0 = 0.0;
    params->kz0 = 0.0;
    prm_set_errfn(hd_silent_warn);
    read_vdiff(params, fp, 1);
    /*
    sprintf(keyword, "KZ0");
    prm_read_double(fp, keyword, &params->kz0);
    sprintf(keyword, "VZ0");
    prm_read_double(fp, keyword, &params->vz0);
    sprintf(keyword, "KZ0");
    prm_read_char(fp, keyword, params->kz0c);
    sprintf(keyword, "VZ0");
    prm_read_char(fp, keyword, params->vz0c);
    prm_set_errfn(hd_silent_warn);
    sprintf(keyword, "KZ_ALPHA");
    prm_read_double(fp, keyword, &params->kz_alpha);
    sprintf(keyword, "VZ_ALPHA");
    prm_read_double(fp, keyword, &params->vz_alpha);
    */
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
    if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
	strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0)
      params->atr += 1;
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

  /* Density (optional)                                              */
  if (prm_read_char(fp, "CALCDENS", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    params->calc_dens = is_true(fields[0]);
    if (n > 1) strcpy(params->densname, fields[1]);      
  }

  /* Surface height                                                  */
  sprintf(params->eta_init, "%c", '\0');
  sprintf(keyword, "SURFACE");
  prm_read_char(fp, keyword, params->eta_init);
  sprintf(keyword, "SURFACE_INV_BAR");
  if (prm_read_char(fp, keyword, buf))
    params->eta_ib = is_true(buf);
  sprintf(params->vel_init, "%c", '\0');
  sprintf(keyword, "VELOCITY");
  prm_read_char(fp, keyword, params->vel_init);

  /* Cells explicitly defined as OUTSIDE                             */
  params->polyland = (char **)malloc(MAXNUMVARS * sizeof(char *));
  for (n = 0; n < MAXNUMVARS; n++) {
    params->polyland[n] = (char *)malloc(sizeof(char)*MAXSTRLEN);
    sprintf(params->polyland[n], "%c", '\0');
  }
  read_blocks(fp, "NLAND", &params->nland, &params->lande1, &params->lande2, params->polyland);

  /* MOM grid conversion                                             */
  if (prm_read_char(fp, "MOM_CONVERT", params->momfile)) {
    not_included("MOM grid conversion");
    /*
    params->domom |= CDF;
    */
  }

  /* ROMS grid conversion                                            */
  if (prm_read_char(fp, "ROMS_CONVERT", params->romsfile)) {
    not_included("ROMS grid conversion");
    /*
    params->doroms |= CDF;
    */
  }

  /* SWAN conversion                                                 */
  if (prm_read_char(fp, "SWAN_CONVERT", buf)) {
    if (is_true(buf)) params->doswan |= 1;
  }

  /* Get the Z to Sigma levels scaling factor for ROMS               */
  /*
  if (params->romsfile) {
    params->roms_z2s = 1.0;
    prm_read_double(fp, "ROMS_Z2S_FACTOR", &params->roms_z2s);
  }
  */

  /* Forcing files                                                   */
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

  /* Heatflux (optional)                                             */
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

  /* SWR parameters (optional)                                       */
  read_swr(params, fp, 1);

  /* Coriolis (optional)                                             */
  params->coriolis = NULL;
  sprintf(keyword, "CORIOLIS");
  if (prm_read_darray(fp, keyword, &params->coriolis, &params->nvals))
    if (params->nce1 && params->nce2 && params->nvals != params->nce1 * params->nce2)
      hd_quit(" Number of values for CORIOLIS incorrect: %d != %d\n",
	      params->nvals, params->nce1 * params->nce2);

  /* Sigma (optional)                                                */
  sprintf(keyword, "SIGMA");
  if (prm_read_char(fp, keyword, buf)) {
    params->sigma = is_true(buf);
    if (params->sigma) {
      params->smagorinsky = 0.1;
      params->atr++;
    }
  }

  /* Robustness (optional)                                           */
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

  /* Pycnocline sharpening                                           */
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

  /* Fatal instability checking (optional)                           */
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

  /* CFL                                                             */
  if (!(params->cfl & NONE))
    params->ntrS += 6;

  /* Momentum tendency (optional)                                    */
  sprintf(keyword, "MOM_TEND");
  if (prm_read_char(fp, keyword, buf)) {
    if ((params->tendf = is_true(buf))) {
      params->atr += 13;
    }
  }

  /* Tracer tendency (optional)                                      */
  sprintf(keyword, "TRA_TEND");
  if (prm_read_char(fp, keyword, params->trtend)) {
    if (strcmp(params->trtend, "NONE") != 0)
      params->atr += 3;
    else
      sprintf(params->trtend, "%c", '\0');
  }

  /* AVHRR SST                                                       */
  if (prm_read_char(fp, "AVHRR", params->avhrr_path)) {
    not_included("AVHRR import");
    hd_warn("Use GHRSST!");
    /*
    create_avhrr_list(params);
    params->avhrr = 1;
    params->ntrS++;
    */
  }

  /* GHRSST SST                                                      */
  if (prm_read_char(fp, "GHRSST", params->ghrsst_path)) {
    prm_read_char(fp, "GHRSST_OPTIONS", params->ghrsst_opt);
    prm_read_char(fp, "GHRSST_INTERP", params->ghrsst_irule);
    create_ghrsst_list(params);
    params->ntrS+=2;
  }

  /* Transport mode files (optional)                                 */
  read_tmode_params(fp, params);
  /* Transport output                                                */
  sprintf(keyword, "TRANS_OUTPUT");
  if (prm_read_char(fp, keyword, buf))
    params->trout = is_true(buf);
  if (!prm_read_char(fp, "OutputTransport", params->trkey)) {
    if (params->trout) {
      strcpy(buf, params->prmname);
      stripend(buf);
      strcpy(params->trkey, buf);
    }
  }

  /* Means (optional)                                                */
  read_means(params, fp, 1);

  /* Fluxes (optional)                                               */
  read_trflux(params, fp);

  /* Particle tracking      */
  if (prm_read_char(fp, "PT_InputFile", params->ptinname)) {
    params->do_pt = PT_DO;
    params->atr++;
    /* Auto particle source                                          */
    if (prm_read_char(fp, "particles", params->particles))
      params->do_pt = PT_AUTO;
    params->do_lag = params->do_pt;
  }

  /* Diagnistic numbers (optional)                                   */
  params->atr += numbers_init(params);
  params->atr += import_init(params, fp);

  params->atr += read_dhw(params, fp);

  /* Auto point source                                               */
  if (prm_read_char(fp, "pss", buf)) {
    params->numbers |= PASS;
    params->ntr += 1;
  }

  /* River flow diagnostic                                           */
  if (check_river_bdry(params, fp)) {
    params->riverflow = 1;
    params->ntrS++;
    if (check_bdry_options(params, fp) & OP_DYNAHC) {
      params->riverflow = 2;
      params->ntrS++;
      params->atr++;
    }
  }

  /* Bathymetry checks (optional)                                    */
  sprintf(keyword, "BATHYFILL");
  if (prm_read_char(fp, "BATHYFILL", buf)) {
    if (strcmp(buf, "REMOVE") == 0)
      params->bathyfill = B_REMOVE;
    if (strcmp(buf, "MINIMUM") == 0)
      params->bathyfill = B_MINIMUM;
    if (strcmp(buf, "AVERAGE") == 0)
      params->bathyfill = B_AVERAGE;
  }
  sprintf(keyword, "MAXGRAD");
  prm_read_double(fp, keyword, &params->maxgrad);
  sprintf(keyword, "MAXDIFF");
  prm_read_char(fp, keyword, params->maxdiff);
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

  /* Time steps (optional)                                           */
  sprintf(keyword, "DT");
  prm_get_time_in_secs(fp, keyword, &params->grid_dt);
  sprintf(keyword, "IRATIO");
  prm_read_int(fp, keyword, &params->iratio);
  if ((params->iratio + 1) % 2 == 0) {
    params->iratio++;
    hd_warn("Only even iratio values allowed : setting iratio=%d\n",
            params->iratio);
  }

  /* Horizontal mixing (optional)                                    */
  read_hdiff(params, fp, 1);

  /* Minimum surface layer thickness (optional)                      */
  sprintf(keyword, "HMIN");
  prm_read_double(fp, keyword, &params->hmin);

  /* Bottom roughness (optional)                                     */
  sprintf(keyword, "QBFC");
  prm_read_char(fp, keyword, params->quad_bfc);
  z0_init(params);

  /* Surface relaxation (optional)                                   */
  read_eta_relax(params, fp);
  read_vel_relax(params, fp);

  /* Sediment layers (optional)                                      */
  /* Note : sediment layers are read in assuming the first layer is  */
  /* closest to the water column. This is reversed when copying to   */
  /* the master so that the sediment origin is the deepest layer.    */
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
  params->do_wave = NONE;
  sprintf(keyword, "DO_WAVES");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "NONE") == 0)
      params->do_wave = NONE;
    if (strcmp(buf, "FILE") == 0)
      params->do_wave = W_FILE;
    if (strcmp(buf, "COMP") == 0)
      params->do_wave = W_COMP;
    if (strcmp(buf, "SWAN") == 0)
      params->do_wave = (W_SWAN|W_SWANM);
    if (strcmp(buf, "SWAN_W") == 0)
      params->do_wave = (W_SWAN|W_SWANW);
    prm_get_time_in_secs(fp, "WAVES_DT", &params->wavedt);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  read_ecology(params, fp, &params->atr);
#endif

  /* Waves                                                           */
  read_waves(params, fp, 1);
  if (params->waves & STOKES_DRIFT && params->tendf) {
    params->atr+=2;
  }

  /* netCDF dumpfile name                                            */
  if (params->runmode & AUTO) {
    sprintf(keyword, "INPUT_FILE");
    if (!prm_read_char(fp, keyword, params->oname))
      strcpy(params->oname, "auto");
  } else {
    strcpy(params->oname, oname);
  }

  /* Output time unit (optional)                                     */
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
  /* Tracer constants and variables                                  */
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "NTRACERS");

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
  sprintf(keyword, "AUTOTRACERPATH");
  prm_read_char(fp, keyword, params->autotrpath);

  /* Reset the pre-tracer ROAM parameterisation if required          */
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
      else if (strcmp(buf, "ROAMv4") == 0)
	params->roammode = A_ROAM_R4;
      else if (strcmp(buf, "ROAMv5") == 0)
	params->roammode = A_ROAM_R5;
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
    else if (params->roammode == A_ROAM_R4)
      auto_params_roam_pre3(fp, params);
    else if (params->roammode == A_ROAM_R5)
      auto_params_roam_pre4(fp, params);
    else
      auto_params_roam_pre1(fp, params);
  }
  read_trfilter(params, fp);

  /* Set up tracers                                                  */
  sprintf(keyword, "NTRACERS");
  if (prm_read_int(fptr, keyword, &n)) {
    params->atr += 2;   /* Temp and salt always auto */
    params->ntr = params->atr;
    tracer_setup(params, fptr);
    create_tracer_3d(params);
    for (n = 0; n < params->atr; n++) {
      params->trinfo_3d[n].m = -1;
      params->trinfo_3d[n].n = n;
    }
  } else {
    params->ntr = 2 + params->atr;
    params->atr = 0;
    sz = sizeof(tracer_info_t) * params->ntr;
    params->trinfo_3d = (tracer_info_t *)malloc(sz);
    memset(params->trinfo_3d, 0, sz);
    create_tracer_3d(params);
    for (n = 0; n < params->ntr; n++)
      params->trinfo_3d[n].n = params->trinfo_3d[n].m;
  }

  /* Tracer relaxation                                               */
  params->trrlxn = malloc(params->ntr * sizeof(char *));
  params->trrlxdt = malloc(params->ntr * sizeof(double));
  params->trrlxtc = malloc(params->ntr * sizeof(char *));
  for (n = 0; n < params->ntr; n++) {
    char rdt[MAXSTRLEN];
    int tm = params->trinfo_3d[n].m;
    tracer_info_t *tracer = &params->trinfo_3d[n];
    sprintf(rdt, "%c", '\0');
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
  /* Tracer reset                                                    */
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
  /* Sediment tracer constants and variables                         */
  prm_set_errfn(hd_silent_warn);
  tracer_read(fp, NULL, SEDIM, hd_quit, hd_warn, hd_silent_warn,
        &params->nsed, &params->trinfo_sed);
  if (params->nsed && !params->sednz)
    hd_quit("Sediment tracers requires NSEDLAYERS > 0\n");

  /*-----------------------------------------------------------------*/
  /* 2D Tracer constants and variables                               */
  prm_set_errfn(hd_silent_warn);
  params->atrS = params->ntrS;
  memset(keyword, 0, sizeof(keyword));
  tracer_read(fp, NULL, INTER, hd_quit, hd_warn, hd_silent_warn,
              &params->ntrS, &params->trinfo_2d);

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  if (prm_read_char(params->prmfd, "INTERIOR_COORD", buf)) {
    sscanf(buf, "%lf %lf", &params->mlon, &params->mlat);
  }
  read_grid(params);

  /*-----------------------------------------------------------------*/
  /* Read the bathymetry                                             */
  /* Note for structured grids:                                      */
  /* gridx[j][i] = params->x[j*2][i*2];                              */
  /* u1x[j][i] = params->x[j*2+1][i*2];                              */
  /* u2x[j][i] = params->x[j*2][i*2+1];                              */
  /* cellx[j][i] = params->x[j*2+1][i*2+1];                          */
  /* h1acell[j][i] = 2.0*params->h1[j*2+1][i*2+1];                   */
  /* h1au1[j][i] = 2.0*params->h1[j*2+1][i*2];                       */
  /* h2au1[j][i] = 2.0*params->h2[j*2+1][i*2];                       */
  /* h1au2[j][i] = 2.0*params->h1[j*2][i*2+1];                       */
  /* h2au2[j][i] = 2.0*params->h2[j*2][i*2+1];                       */
  /* thetau1[j][i] = params->a1[j*2+1][i*2];                         */
  /* thetau2[j][i] = params->a2[j*2][i*2+1] - PI/2.0;                */
  params->bathy = NULL;
  prm_read_double(fp, "BATHYMAX", &params->bmax);
  prm_read_double(fp, "BATHYMIN", &params->bmin);

  /* First check to see if the BATHY parameter exists. If so, then   */
  /* use this in preference to the BATHYFILE parameters.             */
  prm_set_errfn(hd_silent_warn);
  if(!read_bathy_from_file(fp, params)) {
    if (params->us_type & US_IJ && params->roms_grid_opts & ROMS_GRID_BATHY) {
      params->bathy = read_bathy_from_roms(params->roms_grid, nce1, nce2, NULL);
      params->nvals = nce1 * nce2;
    } else if (prm_read_char(fp, "BATHYFILE", buf)) {
      if (endswith(buf, ".nc") || endswith(buf, ".mnc"))
	read_bathy_from_sparse_nc(params, buf);
      else if (endswith(buf, ".bty"))
	read_bathy_from_sparse_bty(params, buf);
      else {
	if (params->us_type & US_IJ) {
	  params->bathy = read_bathy_from_db(fp, params);
	  params->nvals = nce1 * nce2;
	} else {
	  params->bathy = read_bathy_from_db_us(fp, params);
	  params->nvals = params->ns2 + 1;
	}
	if (params->bathy == NULL)
	  hd_quit("auto_params: bathymetry data could not be read!");
      }
    } else if (params->us_type & US_IJ && prm_read_char(fp, "MOM_GRID", buf)) {
      params->bathy = read_bathy_from_mom(buf, nce1, nce2);
      params->nvals = nce1 * nce2;
    } else
      hd_quit("auto_params: No bathymetry data was supplied.");
  }

  /* Mask the batymetry with a coastline if required.                */
  prm_set_errfn(hd_silent_warn);

  /* Read from ROMS first                                            */
  if (params->roms_grid_opts & ROMS_GRID_MASK)
    mask_bathy_from_roms(params->roms_grid, params, NULL);
  else {
    if (prm_read_char(fp, "COASTFILE", buf)) {
      if (params->us_type & US_IJ)
	mask_bathy_with_coast_s(buf, params);
      else
	mask_bathy_with_coast_us(buf, params);
    }
  }

  /* Figure out if bathymetry is entered as positive or negative     */
  /* values.                                                         */
  i = j = 0;
  for (c = 0; c < params->nvals; c++) {
    if (fabs(params->bathy[c]) != fabs(LANDCELL) && fabs(params->bathy[c]) != fabs(NOTVALID)) {
      if (params->bathy[c] >= 0.0)
	i++;
      else
	j++;
    }
  }
  bsgn = (i > j) ? 1.0 : -1.0;

  /* Insert land cells for unstructured meshes if required           */  
  /*
  if (!(params->us_type & US_IJ)) {
    for (n = 1; n <= params->nland; n++) {
      c = params->lande1[n];
      params->bathy[c] = LANDCELL;
    }
  }
  */

  /* Interpolate the bathymetry onto the grid                        */
  if (params->us_type & US_IJ) {
    bathy = d_alloc_2d(params->nce1, params->nce2);
    n = 0;
    is = js = 0;
    ie = params->nce1;
    je = params->nce2;

    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	bathy[j][i] = LANDCELL;
    for (j = js; j < je; j++)
      for (i = is; i < ie; i++)
	bathy[j][i] = -bsgn * params->bathy[n++];
  }

  prm_set_errfn(hd_silent_warn);
  if (!prm_read_double(fp, "BATHYMAX", &params->bmax)) {
    params->bmax = 0;
    for (c = 0; c < params->nvals; c++) {
      if (bsgn * params->bathy[c] > params->bmax && 
	  bsgn * params->bathy[c] != NOTVALID && 
	  bsgn * params->bathy[c] <= MAXDEPTH)
        params->bmax = bsgn * params->bathy[c];
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
      if (bsgn * params->bathy[c] <= MAXDEPTH && bsgn * params->bathy[c] > params->bmax)
        params->bathy[c] = bsgn * params->bmax;
    }
  }

  if (!prm_read_double(fp, "BATHYMIN", &params->bmin)) {
    params->bmin = 1e10;
    for (c = 0; c < params->nvals; c++) {
      if (bsgn * params->bathy[c] < params->bmin && bsgn * params->bathy[c] > 0)
        params->bmin = bsgn * params->bathy[c];
    }
  }
  prm_read_char(fp, "MIN_CELL_THICKNESS", params->mct);

  /* Read the vertical grid from the input file if required          */
  if (prm_read_darray(fp, "LAYERFACES", &params->layers, &params->nz)) {
    if (--params->nz < 1)
      hd_quit("Number of layers in vertical must be 1 or more\n");
  } else if (prm_read_char(fp, "MOM_GRID", buf)) {
    params->layers = read_mom_layers(buf, &params->nz);
  } else {
    /* Set the vertical grid spacing to the default                  */
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
  if (!prm_skip_to_end_of_key(fp, "HMIN")) {
    if (params->hmin == 0.0)
      params->hmin = min(0.1, 0.07 * (params->layers[params->nz] -
				      params->layers[params->nz - 1]));
  }

  if (DEBUG("init_m"))
    dlog("init_m", "\nBathymetry processed OK\n");

  /*-----------------------------------------------------------------*/
  /* Explicit mappings                                               */
  read_explicit_maps(params, fp);

  /*-----------------------------------------------------------------*/
  /* Tidal energy extraction                                         */
  read_turb(params, fp);

  /*-----------------------------------------------------------------*/
  /* Grid conversion                                                 */
  if (params->uscf & (US_TRI|US_HEX)) {
    convert_quad_mesh(params);
    write_mesh_desc(params, NULL, NULL, M_MESHSTRUCT);
  }

  /*-----------------------------------------------------------------*/
  /* Open boundaries                                                 */
  /* If boundary info is included in the input file, read and return */
  get_bdry_params(params, fp);
  
  if (prm_read_int(fp, "NBOUNDARIES", &params->nobc)) {
    if (params->nobc == 0) {
      if (bathy != NULL) d_free_2d(bathy);
      return (params);
    }
    /* Allocate memory for the open boundary data structure          */
    params->open =
      (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);

    for (n = 0; n < params->nobc; n++) {

      /* Allocate memory for the boundary structures                 */
      params->open[n] = OBC_alloc();
      /* Read the open boudary conditions                            */
      params->open[n]->ntr = params->ntr;
      params->open[n]->atr = params->atr;
      get_OBC_conds(params, params->open[n], fp, n, params->trinfo_3d);
      /* Set the flags for this open boundary                        */
      get_obc_list(params->open[n], fp, n, "BOUNDARY");
      prm_set_errfn(hd_quit);
    }
    if (bathy != NULL) d_free_2d(bathy);
    /* Read the csr tide model paths if required                     */
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
	sprintf(keyword, "TIDE_CONSTITUENTS_TRAN");
	if (prm_read_char(fp, keyword, params->tide_con_file)) 
	  params->tidef |= TD_TRAN;
	sprintf(keyword, "TIDE_CONSTITUENTS_VEL");
	if (prm_read_char(fp, keyword, params->tide_con_file)) 
	  params->tidef |= TD_VEL;
      }
    }
    return (params);
  }

  /* Open boundary numbers are assigned in the following order:      */
  /* Edges : west (i=0), east (i=nce1), south (j=0) north (j=nce2)   */
  /* Interior : checks west, east, south north for each cell over    */
  /* all cells in the domain.                                        */
  /* The numbering in the auto input file should reflect this        */
  /* Count the number of OBCs                                        */
  if (params->us_type & US_IJ)
    bdry = count_auto_OBC(params, bathy, &nib, &mask, &nr, &rs, 
		   &ib, &jb, &rb, &rf, &rhc);
  else
    bdry = count_auto_OBC_us(params, &nib, &jci, &masku);

  if (params->nobc > MAXBDRY) 
    hd_quit("Number of boundaries found > allowed maximum (%d).\n", MAXBDRY);

  /* Allocate memory for the open boundary data structure            */
  params->open =
    (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * params->nobc);
  for (n = 0; n < params->nobc; n++) {
    /* Allocate memory for the boundary structures                   */
    params->open[n] = OBC_alloc();
    open = params->open[n];
    open->id = n;
    open->ntr = params->ntr;
    open->atr = params->atr;

    /* Read the open boudary conditions                              */
    dataf = 0;
    init_OBC_conds(params, open);

    open->relax_zone_tan = 0;
    open->relax_zone_ele = 0;
    open->bathycon = 0;
    open->relax_time = 0;
    open->relax_timei = 0;
    open->sponge_zone = 0;
    open->sponge_zone_h = 0.0;
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
  /* Read the csr tide model paths if required                       */
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
      sprintf(keyword, "TIDE_CONSTITUENTS_TRAN");
      if (prm_read_char(fp, keyword, params->tide_con_file)) 
	params->tidef |= TD_TRAN;
      sprintf(keyword, "TIDE_CONSTITUENTS_VEL");
      if (prm_read_char(fp, keyword, params->tide_con_file)) 
	params->tidef |= TD_VEL;
      break;
    }
  }

  /* Set the OBC locations                                           */
  if (params->us_type & US_IJ)
    set_auto_OBC(params, bathy, bdry, nib, mask, &nr, &rs, 
		 ib, jb, rb, rf, rhc);
  else
    set_auto_OBC_us(params, nib, jci, masku);

  /* Read the flow data file if explicitly specified                 */
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

    /* Set the river OBCs if specified via RIVER (for structured     */
    /* grids only).                                                  */
    if (params->us_type & US_IJ) {
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
    }

    if (DEBUG("init_m"))
      dlog("init_m", "\nOBCs processed OK\n");

    /*---------------------------------------------------------------*/
    /* Tracers                                                       */
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

  /* Reset the post-tracer ROAM parameterisation if required         */
  if (params->runmode & ROAM) {
    if (params->roammode == A_ROAM_CPD1)   /* Standard ROAM          */
      auto_params_roam_post1(fp, params);
    if (params->roammode == A_ROAM_CPD2)   /* ETA_RELAX ramp, RAYMND bcond_nor */
      auto_params_roam_post2(fp, params);
    if (params->roammode == A_ROAM_FLA)    /* ETA_RELAX ramp, FLATHR */
      auto_params_roam_post3(fp, params);
    if (params->roammode == A_ROAM_R1)     /* ETA_RELAX ramp, velocity forced */
      auto_params_roam_post4(fp, params);
    if (params->roammode == A_ROAM_R2)     /* A_ROAM_R1 with alternative robust parameterisations */
      auto_params_roam_post5(fp, params);
    if (params->roammode == A_ROAM_R3)     /* A_ROAM_R2 with alternative robust parameterisations */
      auto_params_roam_post5(fp, params);
    if (params->roammode == A_ROAM_R4)     /* A_ROAM_R3 with TPXO tide                            */
      auto_params_roam_post6(fp, params);
    if (params->roammode == A_ROAM_R5)   /* A_ROAM_R5 with TPXO tide and dual relaxation */
      auto_params_roam_post7(fp, params);
    if (params->roammode == A_RECOM_R1)    /* RECOM                  */
      auto_params_recom_post1(fp, params);
    if (params->roammode == A_RECOM_R2)    /* RECOM + ROBUST         */
      auto_params_recom_post2(fp, params);
  }

  if (bathy) d_free_2d(bathy);
  if (mask) i_free_2d(mask);
  if (masku) i_free_1d(masku);
  if (nib) i_free_1d(nib);
  if (jci) i_free_2d(jci);
  if (ib) i_free_1d(ib);
  if (jb) i_free_1d(jb); 
  if (rb) i_free_1d(rb);
  if (rf) i_free_1d(rf);
  return (params);
}

/* END auto_params()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Automatically sets mandatory model parameters                     */
/*-------------------------------------------------------------------*/
void autoset(parameters_t *params, master_t *master, geometry_t *geom)
{
  int n, cc, c, lc, np, m;      /* Sparse counters                   */
  int vv, v, ee, e;             /* Sparse counters                   */
  int i, j, k;                  /* Cartesian counters                */
  double spd = 2.0;             /* Maximum advective velocity expected */
  double cif = 2.0;             /* Scaling factor for internal wave speed */
  double sf = 0.8;              /* Safety factor for cfl calculations */
  double hf = 0.05;             /* Factor for horizontal diffusion   */
  double sfr = 1.0;             /* Depth dependent safety factor     */
  double eps = 1e5;             /* Large value for cfl calculations  */
  double cfl2d, cfl3d;          /* CFL timesteps                     */
  double hmax;                  /* Maximum horizontal diffusion      */
  double lat, lon;              /* Latitude and longitude for Coriolis */
  double d1, d2, d3;            /* Dummies                           */
  dump_data_t *dumpdata = master->dumpdata;

  /*-----------------------------------------------------------------*/
  /* Time                                                            */
  master->t = dumpdata->t = params->t;

  /*-----------------------------------------------------------------*/
  /* Vertical grid geometry. The bathymetry is set in build_sparse_map() */
  /* and overwritten in dump_read() if no automation is performed.  */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    k = geom->s2k[c];
    geom->gridz[c] = params->layers[k];
    geom->cellz[c] = 0.5 * (params->layers[k] + params->layers[k + 1]);
  }
  for (cc = geom->b3_t + 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    geom->gridz[c] = geom->gridz[geom->wgst[c]];
    geom->cellz[c] = geom->cellz[geom->wgst[c]];
  }    

  /* Sediments                                                       */
  /* Note : sediment layers are read in assuming the first layer is  */
  /* closest to the water column. This is reversed here so that the  */
  /* sediment origin is the deepest layer.                           */
  if (geom->sednz) {
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[c];
      if (geom->wgst[c])
        continue;               /* Continue on ghosts                */
      for (k = 0; k < geom->sednz; k++) {
	geom->gridz_sed[k][c] = params->gridz_sed[k];
      }
      for (k = geom->sednz - 1; k >= 0; k--) {
	geom->cellz_sed[k][c] = 0.5 * (geom->gridz_sed[k][c] +
					geom->gridz_sed[k + 1][c]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Cell mappings for output                                        */
  neighbour_none(dumpdata);

  /*-----------------------------------------------------------------*/
  /* Horizontal grid geometry                                        */
  dumpdata_init_geom(params, geom, dumpdata);

  /*-----------------------------------------------------------------*/
  /* Load the tracers. Note hd_ts_read() is used to read the tracers */
  /* from file and this requires time variables to be set on the     */
  /* master. These are only set in compute_constants(), called after */
  /* this routine, so preset them here.                              */
  master->t = master->tstart = params->t;
  strcpy(master->timeunit, params->timeunit);
  strcpy(master->output_tunit, params->output_tunit);
  load_tracer_step_3d(params, master, params->prmfd);
  load_tracer_step_2d(params, master, params->prmfd);

  /*-----------------------------------------------------------------*/
  /* Get the initial density distribution                            */
  density_m(master);

  /*-----------------------------------------------------------------*/
  /* Get the CFL conditions                                          */
  if (params->grid_dt == 0.0) {
    double btws, btwsm = 0.0;
    double bcws, bcwsm = 0.0;

    cfl2d = cfl3d = 1e10;
    if (params->robust >= 6)
      sf *= (1.0 - 0.5 * ((double)(params->robust) - 6.0) / 4.0);

    for (cc = 1; cc <= geom->v2_t; cc++) {
      c = lc = geom->w2_t[cc];

      if (params->us_type & US_IJ) {
	d1 = 0.0;
	for (j = 1; j <= geom->npe[c]; j++) {
	  e = geom->c2e[j][c];
	  if (c == geom->e2c[e][0]) {
	    if (geom->h1acell[e])
	      d1 += 1.0 / (geom->h1acell[e] * geom->h1acell[e]);
	  }
	}
	d1 = sqrt(d1);
      } else
	d1 = sqrt((double)geom->npe[c] / (2.0 * geom->cellarea[c]));

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
    else if (sf * cfl3d < 20)
      d1 = 1;
    else
      d1 = 5;
    c = (int)cfl3d *sf / d1;
    params->grid_dt = d1 * (double)c;
    c = (int)cfl2d *sf;
    if (c > 1.0) {
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
    } else {
      c = (int)(cfl2d *10) * sf;
      cfl2d = (double)c / 10.0;
      hd_warn("floor(2D time-step) <= 1 (%f)\n", cfl2d);
    }
    params->iratio = max(1, (int)(params->grid_dt / cfl2d));
    if (params->iratio % 2 != 0)
      params->iratio++;
    params->grid_dt = (double)params->iratio * cfl2d;
    if (DEBUG("init_m")) {
      dlog("init_m : auto dt","Scaled cfl3d = %5.1f",params->grid_dt);
      dlog("init_m : auto dt","Scaled cfl2d = %5.1f",cfl2d);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Horizontal mixing coefficients                                  */
  /* Note: Stability criterion for a quad is:                        */
  /* 1/[(1/h1*h1 + 1/h2*h2) * 4 * dt]. For arbitary polgons assume   */
  /* h1 = h2 = sqrt(area), then the stability criterion              */
  /* is 1/[2/area] * 4 * dt.                                         */
  n = 0;
  d1 = d2 = 1e-10;
  hmax = 1e10;

  for (cc = 1; cc <= geom->v2_t; cc++) {
    c = geom->w2_t[cc];

    if (params->us_type & US_IJ) {
      d2 = 0.0;
      for (j = 1; j <= geom->npe[c]; j++) {
	e = geom->c2e[j][c];
	if (c == geom->e2c[e][0]) {
	  if (geom->h1acell[e])
	    d2 += (geom->h1acell[e] * geom->h1acell[e]);
	}
      }
    } else {
      d2 = geom->cellarea[c];
    }
    d1 += d2;
    d3 = (d2) ? d2 / (8.0 * params->grid_dt) : 1e10;
    if (d3 < hmax)
      hmax = d3;
    n += 1;
  }
  d1 /= (double)n;
  hmax = (d1) ? d1 / (8.0 * params->grid_dt) : 1e10;
  if (hmax > 5.0) {
    c = (int)hmax / 5;
    hmax = 5.0 * (double)c;
  }

  if (n) {
    int u1khf = 0, u1vhf = 0;
    /*if (params->u1kh == 0.0)*/
    if (!strlen(params->u1khc) && !strlen(params->smag))
      u1khf = 1;
    /*if (params->u1vh == 0.0)*/
    if (!strlen(params->u1khc) && !strlen(params->smag))
      u1vhf = 1;
    d1 = 0.01 * d1 / params->grid_dt;

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
    /* Set limits                                                    */
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
  } else
    hd_warn("Horizontal mixing is set to zero\n");

  if(params->smagorinsky > 0.0) {
    params->u1vh *= -1.0;
    params->u2vh *= -1.0;
  }

  /*-----------------------------------------------------------------*/
  /* Bottom roughness                                                */
  if (params->z0 == 0.0) {
    double qbfc = atof(params->quad_bfc);
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
    n = 1e4 * d1 / (exp(0.4 * sqrt(1.0 / fabs(qbfc))) - 1);
    params->z0 = (double)n / 1e4;
  }

  /*-----------------------------------------------------------------*/
  /* Coriolis                                                        */
  memset(master->coriolis, 0, geom->sgsizS * sizeof(double));
  memset(master->fv, 0, geom->szvS * sizeof(double));
  if (params->coriolis) {
    if (params->us_type & US_IJ) {
      n = 0;
      for (j = 0; j < geom->nce2; j++)
	for (i = 0; i < geom->nce1; i++) {
	  dumpdata->coriolis[j][i] = params->coriolis[n++];
	  c = geom->map[params->nz - 1][j][i];
	  if (c) {
	    c = geom->m2d[c];
	    master->coriolis[c] = dumpdata->coriolis[j][i];
	  }
	}
    } else {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	n = geom->c2cc[c] - 1;
	master->coriolis[c] = params->coriolis[n];
      }
    }
    for (vv = 1; vv <= geom->n2_e2; vv++) {
      v = geom->w2_e2[vv];
      d1 = d2 = 0.0;
      for (n = 1; n <= geom->nvc[v]; n++) {
	c = geom->v2c[v][n];
	if (c && !geom->wgst[c]) {
	  d2 += master->coriolis[c];
	  d1 += 1.0;
	}
      }
      if (d1) master->fv[v] = d2 / d1;
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
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      lon = geom->cellx[c];
      lat = geom->celly[c];
      if (mp != NULL)
	mp_inverse(mp, lon, lat, &lat, &lon);
      master->coriolis[c] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
    }
    for (vv = 1; vv <= geom->n2_e2; vv++) {
      v = geom->w2_e2[vv];
      lon = geom->gridx[v];
      lat = geom->gridy[v];
      if (mp != NULL)
	mp_inverse(mp, lon, lat, &lat, &lon);
      master->fv[v] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
    }
  } else {
    lat = -30.0;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->coriolis[c] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
      if (params->us_type & US_IJ) {
	i = geom->s2i[c];
	j = geom->s2j[c];
	dumpdata->coriolis[j][i] = master->coriolis[c];
      }
    }
    for (vv = 1; vv <= geom->n2_e2; vv++) {
      v = geom->w2_e2[vv];
      master->fv[v] = 2 * 7.29e-05 * sin(lat * M_PI / 180.0);
    }
  }
  dumpdata_init(dumpdata, geom, master);

  /*-----------------------------------------------------------------*/
  /* Initialize                                                      */
  memset(master->u1av, 0, geom->szeS * sizeof(double));
  memset(master->u1bot, 0, geom->szeS * sizeof(double));
  memset(master->wind1, 0, geom->szeS * sizeof(double));
  memset(master->wtop, 0, geom->sgsizS * sizeof(double));
  memset(master->topz, 0, geom->sgsizS * sizeof(double));
  memset(master->eta, 0, geom->sgsizS * sizeof(double));
  memset(master->patm, 0, geom->sgsizS * sizeof(double));
  memset(master->u1, 0, geom->sze * sizeof(double));
  memset(master->w, 0, geom->sgsiz * sizeof(double));

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    master->Vz[c] = params->vz0;
    master->Kz[c] = params->kz0;
  }

  /* Initialise the surface elevation if required                    */
  eta_init(geom, params, master);

  /* Initialise the velocity if required                             */
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
  int hint = 1;        /* Interpolation type : 0 = AUTO              */
                       /*                      1 = NEAREST           */
                       /*                      2 = LINEAR            */
                       /*                      3 = NATURAL NEIGHBOURS*/
                       /*                      4 = INVERSE           */
                       /*                      5 = AVERAGE           */

  /* Adjust for grid edges */
  is = js = 0;
  ie = nce1;
  je = nce2;

  prm_set_errfn(hd_quit);

  /* Read the bathymetry data.                                       */
  if (prm_read_char(fp, "BATHYFILE", bfname)) {
    topo_files_t *tf;
    double **x = d_alloc_2d(nce1+1, nce2+1);
    double **y = d_alloc_2d(nce1+1, nce2+1);

    /* Reduce double grid to single grid corners                     */
    if (is_geog) {
      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i;
          jm = j;
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
          im = i;
          jm = j;
          mp_inverse(mp, params->x[2*j][2*i], params->y[2*j][2*i],
              &y[jm][im], &x[jm][im]);
        }
      }
    }

    /* Check ranges; must be 0 - 360                                 */
    for (j = 0; j < nce2 + 1; j++) {
      for (i = 0; i < nce1 + 1; i++) {
	while (x[j][i] < 0.0 || x[j][i] > 360.0) {
	  if (x[j][i] < 0.0) x[j][i] += 360.0;
	  if (x[j][i] > 360.0) x[j][i] -= 360.0;
	}
      }
    }

    /* Populate the grid                                             */
    tf = topo_simple_init(bfname);

    src_topo = topo_get_z_for_grid(tf, x, y, nce1, nce2, (topo_hint_t)hint);

    d_free_2d(y);
    d_free_2d(x);
    topo_destroy(tf);

    if (src_topo == NULL)
      hd_quit("read_bathy_from_db: unable to populate bathymetric grid.");
  }

  /* Convert 2D bathy grid to a 1d bathymetry array.                 */
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


static double *read_bathy_from_db_us(FILE *fp, parameters_t *params) {

  char bfname[MAXFNAMELEN];
  double *topo = NULL;
  double *dst_bathy = NULL;
  int has_proj = (strlen(params->projection) > 0);
  int is_geog = has_proj && (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
  int i, j, n;
  int is, js, ie, je, jm, im;
  int ns2 = params->ns2;
  int hint = 1;        /* Interpolation type : 0 = AUTO              */
                       /*                      1 = NEAREST           */
                       /*                      2 = LINEAR            */
                       /*                      3 = NATURAL NEIGHBOURS*/
                       /*                      4 = INVERSE           */
                       /*                      5 = AVERAGE           */

  prm_set_errfn(hd_quit);

  /* Read the bathymetry data.                                       */
  if (prm_read_char(fp, "BATHYFILE", bfname)) {
    topo_files_t *tf;
    double *x = d_alloc_1d(ns2+1);
    double *y = d_alloc_1d(ns2+1);

    /* Reduce double grid to single grid corners                     */
    if (is_geog) {
      for (i=1; i<=ns2; i++) {
	x[i] = params->x[i][0];
	y[i] = params->y[i][0];
      }
    }
    /* Check ranges; must be 0 - 360                                 */
    for (i=1; i<=ns2; i++) {
      while (x[i] < 0.0 || x[i] > 360.0) {
	if (x[i] < 0.0) x[i] += 360.0;
	if (x[i] > 360.0) x[i] -= 360.0;
      }
    }
    
    /* Populate the grid                                             */
    tf = topo_simple_init(bfname);

    topo = topo_get_z_for_grid_us(tf, x, y, ns2, (topo_hint_t)hint);

    d_free_1d(y);
    d_free_1d(x);
    topo_destroy(tf);

    if (topo == NULL)
      hd_quit("read_bathy_from_db_us: unable to populate bathymetric grid.");
  }

  /* Convert 2D bathy grid to a 1d bathymetry array.                 */
  dst_bathy = d_alloc_1d(ns2+1);
  for (i=1; i<=ns2; i++) {
    dst_bathy[i] = topo[i];
    if (topo[i] > 0)
      dst_bathy[i] = LANDCELL;
  }

  if (topo != NULL)
    d_free_1d(topo);

  return dst_bathy;
}

/* END read_bathy_from_db()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract bathymetry from the database bathy file        */
/*-------------------------------------------------------------------*/
static int read_bathy_from_file(FILE *fp, parameters_t *params) {
  FILE* fb = fp;
  char buf[MAXLINELEN];
  int cc, cn, is, js, ie, je, jm, im;
  int nce1 = params->nce1;
  int nce2 = params->nce2;

  /* Adjust for grid edges                                           */
  is = js = 0;
  ie = nce1;
  je = nce2;

  if(!prm_read_char(fb, "BATHY", buf)) {
    if(prm_read_char(fb, "BATHYEXTFILE", buf))
      fb=fopen(buf,"r");
    else {
      hd_warn("No bathymetry provided in parameter file.");
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

  if (params->us_type & US_RS) {
    if (params->nvals != nce1 * nce2)
      hd_quit("read_bathy: incorrect BATHY data : %d  (requires %d)\n", 
	      params->nvals, nce1 * nce2);
  } else {
    if (params->nvals != params->ns2 + params->nland)
      hd_quit("read_bathy: incorrect BATHY data : %d  (requires %d)\n", 
	      params->nvals, params->ns2);
  }

  if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
    double *bathy;
    int *mask = i_alloc_1d(params->nvals);
    /* Set the land override mask                                    */
    memset(mask, 0, params->nvals * sizeof(int));
    for (cc = 1; cc <= params->nland; cc++) {
      mask[params->lande1[cc]] = cc;
    }
    /* Re-order the bathymetry array                                 */
    bathy = d_alloc_1d(params->nvals);
    memcpy(bathy, params->bathy, params->nvals * sizeof(double));
    d_free_1d(params->bathy);
    params->bathy = d_alloc_1d(params->ns2 + 1);
    cn = 1;
    for (cc = 0; cc < params->ns2; cc++) {
      if (!mask[cc]) params->bathy[cn++] = bathy[cc];
    }
    d_free_1d(bathy);
    i_free_1d(mask);
    if (params->mrf & MR_READ) reorder_bathy(params);
  }
  return 1;
}

/* END read_bathy_from_file()                                        */	
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates bathymetry from netcdf file for unstructured grids   */
/*-------------------------------------------------------------------*/
static int read_bathy_from_nc(parameters_t *params, char *fname)
{
  char buf[MAXSTRLEN];
  int cc, i, j;
  timeseries_t *ts = NULL;
  int idb;

  /* Initialize the timeseries file                                  */
  ts = (timeseries_t *)malloc(sizeof(timeseries_t));
  if (ts == NULL)
    hd_quit("hd_ts_read: No memory available.\n");
  memset(ts, 0, sizeof(timeseries_t));

  /* Read the time series                                            */
  ts_read(fname, ts);
  if ((idb = ts_get_index(ts, fv_get_varname(fname, "botz", buf))) == -1)
    hd_quit("read_bathy_from_nc: Can't find variable botz in file %s\n", fname);

  /* Read the values                                                 */
  if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
    params->bathy = d_alloc_1d(params->ns2 + 1);
    for (cc = 1; cc <= params->ns2; cc++) {
      params->bathy[cc] = ts_eval_xy(ts, idb, 0.0, params->x[cc][0], params->y[cc][0]);
    }
    params->nvals = cc-1;
  } else {
    params->bathy = d_alloc_1d(params->nce1 * params->nce2);
    cc = 0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++) {
	params->bathy[cc] = -ts_eval_xy(ts, idb, 0.0, params->x[2*j][2*i], params->y[2*j][2*i]);
	cc++;
      }
    params->nvals = params->nce1 * params->nce2;
  }

  /*ts_free((timeseries_t*)ts);
  free(ts);*/
}

/* END read_bathy_from_nc()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured bathymetry from a netCDF file, packs into a      */
/* vector and interpolates onto an unstructured mesh.                */
/*-------------------------------------------------------------------*/
void read_bathy_from_sparse_nc(parameters_t *params, char *fname)
{
  char buf[MAXSTRLEN], fnamei[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  char *rules[MAXSTRLEN * MAXNUMARGS];
  GRID_SPECS *gs = NULL;
  int nbath;
  double *x, *y, *b, **botz, **cellx, **celly, *celx, *cely;
  double bmean;
  int fid;
  int ncerr;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  int n, m, i, j, cc, c;
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int bverbose = 0;
  int intype = -1;
  int dof = 0;
  int spf = 0;           /* Sparse data interpolation                */
  int nf = 0;            /* Number of bathymetry files               */
  int nr = 0;            /* Number of i_rules                        */
  int idb;
  poly_t *pl;

  sprintf(i_rule, "%c", '\0');
  if (prm_read_char(params->prmfd, "BATHY_INTERP_RULE", i_rule)) {
    spf = 1;
    nr = parseline(i_rule, rules, MAXNUMARGS);
  }
  if (endswith(fname, ".mnc")) {
    FILE *fp;
    char key[MAXSTRLEN];
    double vers;
    if ((fp = fopen(fname, "r")) == NULL)
      hd_quit("Can't find bathy file %s\n", fname);
    sprintf(fnamei, "%c", '\0');
    sprintf(i_rule, "%c", '\0');
    prm_read_double(fp, "multi-netcdf-version", &vers);
    prm_read_int(fp, "nfiles", &n);
    for (i = 0; i < n; i++) {
      sprintf(key, "file%d.filename", i);
      prm_read_char(fp, key, buf);
      m = parseline(buf, fields, MAXNUMARGS);
      sprintf(fnamei, "%s %s", fnamei, fields[0]);
      if (m == 2) sprintf(i_rule, "%s %s", i_rule, fields[1]);
    }
    spf = 1;
    nr = parseline(i_rule, rules, MAXNUMARGS);
    fclose(fp);
  } else
    strcpy(fnamei, fname);
  nf = parseline(fnamei, fields, MAXNUMARGS);

  if (nr > 1 && nr != nf)
    hd_quit("read_bathy_from_sparse_nc: Number of INTERP_RULE must equal 1 or number of BATHYFILES.\n");
  m = 0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
    params->bathy = d_alloc_1d(params->ns2 + 1);
    for (cc = 1; cc <= params->ns2; cc++) params->bathy[cc] = NaN;
  } else {
    params->bathy = d_alloc_1d(params->nce1 * params->nce2);
    cc = 0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	params->bathy[cc++] = NaN;
  }

  /*-----------------------------------------------------------------*/
  /* Interpolate the bathymetry                                      */
  /*-----------------------------------------------------------------*/

  /*-----------------------------------------------------------------*/
  /* Interpolation of the first nf-1 files.                          */
  /* If nf = 1 and no i_rule is supplied then use structured         */
  /* interpolation.                                                  */
  /* If nf > 1 and one i_rule is set then use structured             */
  /* interpolation.                                                  */
  /* If nf > 1 and  multiple i_rules are supplied then use the       */
  /* corresponding rule. If i_rule is 'none' then also use           */
  /* structured interpolation.                                       */
  /* i_rule options:                                                 */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (nf > 1 || !spf) {

    if (DEBUG("init_m"))
      dlog("init_m", "\nStructured bathymetry reading from file\n");

    /* If nf > 1 and i_rule is set, then interpolate file 0:nf-2     */
    /* with ts_read() and file nf-1 using i_rule.                    */
    n = (nf > 1 && spf) ? nf - 1 : nf;

    for (m = 0; m < n; m++) {
      if (spf && nr > 1 && strcmp(rules[m], "none") != 0)
	bathy_interp_us(params, fields[m], rules[m], 1);
      else
	bathy_interp_s(params, fields[m], nf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Interpolate the last file with unstructured methods if an       */
  /* i_rule is supplied.                                             */
  if (spf) {
    bathy_interp_us(params, fields[m], rules[m], 0);
  }

  /* Set any unfilled NaNs to land                                   */
  for (n = 0; n < params->nvals; n++) {
    if (isnan(params->bathy[n])) {
      if (params->us_type & US_RS)
	params->bathy[n] = -LANDCELL;
      else
	params->bathy[n] = LANDCELL;
    }
  }

  /* Print to file if required                                       */
  if (prm_read_char(params->prmfd, "BATHY_PRINT", buf)) {
    FILE *fp1, *fp2;
    sprintf(fnamei, "%s.bath", buf);
    fp1 = fopen(fnamei, "w");
    sprintf(fnamei, "%s.bty", buf);
    fp2 = fopen(fnamei, "w");
    fprintf(fp1, "BATHY %d\n", params->nvals);
    if (params->us_type & US_RS) {
      c = 0;
      for (j = 0; j < params->nce2; j++) {
	for (i = 0; i < params->nce1; i++) {
	  fprintf(fp1, "%f\n", params->bathy[c]);
	  if (params->bathy[c] != -LANDCELL && !isnan(params->x[2*j][2*i]) && !isnan(params->y[2*j][2*i]))
	    fprintf(fp2, "%f %f %f\n", params->x[2*j][2*i], params->y[2*j][2*i], params->bathy[c]);
	  c++;
	}
      }
    } else {
      for (c = 1; c <= params->ns2; c++) {
	fprintf(fp1, "%f\n", params->bathy[c]);
	fprintf(fp2, "%f %f %f\n", params->x[c][0], params->y[c][0], params->bathy[c]);
      }
    }
    fclose(fp1);
    fclose(fp2);
  }

  if (DEBUG("init_m"))
    dlog("init_m", "\nBathymetry read from file OK\n");
}


void read_bathy_from_sparse_nc_o(parameters_t *params, char *fname)
{
  char buf[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  GRID_SPECS *gs = NULL;
  int nbath;
  double *x, *y, *b, **botz, **cellx, **celly, *celx, *cely;
  double bmean;
  int fid;
  int ncerr;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  int n, m, i, j, cc, c;
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int bverbose = 0;
  int intype = -1;
  int dof = 0;
  int spf = 0;           /* Sparse data interpolation                */
  int nf = 0;            /* Number of bathymetry files               */
  int idb;
  poly_t *pl;

  sprintf(i_rule, "%c", '\0');
  if (prm_read_char(params->prmfd, "BATHY_INTERP_RULE", i_rule)) spf = 1;
  nf = parseline(fname, fields, MAXNUMARGS);
  m = 0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
    params->bathy = d_alloc_1d(params->ns2 + 1);
    for (cc = 1; cc <= params->ns2; cc++) params->bathy[cc] = NaN;
  } else {
    params->bathy = d_alloc_1d(params->nce1 * params->nce2);
    cc = 0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	params->bathy[cc++] = NaN;
  }

  /*-----------------------------------------------------------------*/
  /* Interpolate the bathymetry                                      */
  /*-----------------------------------------------------------------*/

  /*-----------------------------------------------------------------*/
  /* Interpolation from structured input using ts_read().            */
  if (nf > 1 || !spf) {
    timeseries_t *ts = NULL;

    if (DEBUG("init_m"))
      dlog("init_m", "\nStructured bathymetry reading from file\n");

    /* If nf > 1 and i_rule is set, then interpolate file 0:nf-2     */
    /* with ts_read() and file nf-1 using i_rule.                    */
    n = (nf > 1 && spf) ? nf - 1 : nf;

    for (m = 0; m < n; m++) {

      /* Open the file and get dimensions                            */
      if ((ncerr = nc_open(fields[m], NC_NOWRITE, &fid)) != NC_NOERR) {
	printf("Can't find input bathymetry file %s\n", fields[m]);
	hd_quit((char *)nc_strerror(ncerr));
      }
      i = 0;
      while (bathy_dims[i][0] != NULL) {
	if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &nce2) == 0 &&
	    nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &nce1) == 0) {
	  intype = i;
	  break;
	}
	i++;
      }
      if (intype < 0)      
	hd_quit("read_bathy_from_sparse_nc: Can't find bathymetry attributes in file %s.\n", fields[m]);

      /* Initialize the timeseries file                              */
      ts = (timeseries_t *)malloc(sizeof(timeseries_t));
      if (ts == NULL)
	hd_quit("read_bathy_from_sparse_nc: No memory available.\n");
      memset(ts, 0, sizeof(timeseries_t));
    
      /* Read the time series                                        */
      ts_read(fields[m], ts);
      if ((idb = ts_get_index(ts, fv_get_varname(fields[m], bathy_dims[intype][0], buf))) == -1)
	hd_quit("read_bathy_from_sparse_nc: Can't find variable %s in file %s\n", 
		bathy_dims[intype][0], fields[m]);

      /* Get a polgon of the perimeter of the bathy file             */
      pl = nc2poly(fid, nce1, nce2, bathy_dims[intype][3], bathy_dims[intype][4], fields[m], ts);

      /* Read the values                                             */
      if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
	/* Unstructured.                                             */
	for (cc = 1; cc <= params->ns2; cc++) {
	  if (isnan(params->bathy[cc])) {
	    /* Note: ts_eval_xy() for bathymetry files (idb = 'botz' */
	    /* or 'height') will return NaN if the xytoij() mapping  */
	    /* does not return valid indices (i.e. (x,y) lies        */
	    /* outside the (i,j) bounds of the file).                */
	    /* This seems not entirely reliable though, so we use a  */
	    /* bounding polygon for the perimeter.                   */
	    /* Check if the location lies within the file's bounding */
	    /* perimeter.                                            */
	    if (!poly_contains_point(pl, params->x[cc][0], params->y[cc][0])) continue;

	    params->bathy[cc] = ts_eval_xy(ts, idb, 0.0, params->x[cc][0], params->y[cc][0]);
	    
	    /*printf("%d %f %d %f : %f %f\n",m,(double)cc/(double)params->ns2,cc,params->bathy[cc],params->x[cc][0], params->y[cc][0]);*/
	    if (limf && !isnan(params->bathy[cc])) {
	      if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
		params->bathy[cc] = -params->bmax;
	      if (params->bmin && fabs(params->bathy[cc]) < params->bmin) 
		params->bathy[c] = -params->bmin;
	      if (params->bmin && params->bathy[cc] > params->bmin) {
		if (nf > 1)
		  params->bathy[cc] = NaN;
		else
		  params->bathy[cc] = LANDCELL;
	      }
	    }
	  }
	}
	params->nvals = cc-1;
      } else {
	/* Structured                                                */
	cc = 0;
	for (j = 0; j < params->nce2; j++)
	  for (i = 0; i < params->nce1; i++) {
	    params->bathy[cc] = -ts_eval_xy(ts, idb, 0.0, params->x[2*j][2*i], params->y[2*j][2*i]);
	    cc++;
	  }
	params->nvals = params->nce1 * params->nce2;
      }
      ts_free((timeseries_t*)ts);
      free(ts);
      poly_destroy(pl);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sparse interpolation method. Options:                           */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (spf) {
    int lond, latd;
    int varid;

    if (DEBUG("init_m"))
      dlog("init_m", "\nUnstructured bathymetry reading from file\n");

    /*---------------------------------------------------------------*/
    /* Open the file and get dimensions                              */
    if ((ncerr = nc_open(fields[m], NC_NOWRITE, &fid)) != NC_NOERR) {
      printf("Can't find input bathymetry file %s\n", fields[m]);
      hd_quit((char *)nc_strerror(ncerr));
    }

    /* Get dimensions                                                  */
    i = 0;
    while (bathy_dims[i][0] != NULL) {
      if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &nce2) == 0 &&
	  nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &nce1) == 0) {
	intype = i;
	break;
      }
      i++;
    }

    varid = ncw_var_id(fid, bathy_dims[intype][3]);
    nc_inq_varndims(fid, varid, &lond);
    varid = ncw_var_id(fid, bathy_dims[intype][4]);
    nc_inq_varndims(fid, varid, &latd);

    if (intype < 0)      
      hd_quit("read_bathy_from_sparse_nc: Can't find bathymetry attributes in file %s.\n", fields[m]);

    /*---------------------------------------------------------------*/
    /* Interpolation from unstructured input using a triangulation.  */
    /* The input bathymetry is structured, but is triangulated then  */
    /* interpolated using i_rule.                                    */
    /* Allocate and read.                                            */
    if (lond == 1)
      celx = d_alloc_1d(nce1);
    else if (lond == 2)
      cellx = d_alloc_2d(nce1, nce2);
    if (latd == 1)
      cely = d_alloc_1d(nce2);
    else if (latd == 2)
      celly = d_alloc_2d(nce1, nce2);
    botz = d_alloc_2d(nce1, nce2);
    start[0] = 0L;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = nce2;
    count[1] = nce1;
    count[2] = 0;
    count[3] = 0;
    if (lond == 1) {
      count[0] = nce1;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, celx);
    } else if (lond == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, cellx[0]);
    }
    if (latd == 1) {
      count[0] = nce2;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, cely);
    } else if (latd == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, celly[0]);
    }
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][0]), start, count, botz[0]);
    nc_close(fid);

    /*---------------------------------------------------------------*/
    /* Set the wet bathymetry vector (to interpolate from)           */
    nbath = n = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) nbath++;
      }
    }

    x = d_alloc_1d(nbath);
    y = d_alloc_1d(nbath);
    b = d_alloc_1d(nbath);
    bmean = 0.0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) {
	  if (lond == 1)
	    x[n] = celx[i];
	  else if (lond == 2)
	    x[n] = cellx[j][i];
	  if (latd == 1)
	    y[n] = cely[j]; 
	  else if (latd == 2)
	    y[n] = celly[j][i]; 
	  b[n] = botz[j][i];
	  bmean += b[n];
	  //if(b[n] > 0.0)printf("%d(%d %d) %f %f %f %f\n",n,i,j,x[n],y[n],b[n],OUTSIDE);
	  n++;
	}
      }
    }
    if (n) bmean /= (double)n;
    if (lond == 1)
      d_free_1d(celx);
    else if (lond == 2)
      d_free_2d(cellx);
    if (latd == 1)
      d_free_1d(cely);
    else if (latd == 2)
      d_free_2d(celly);
    d_free_2d(botz);

    /* Interpolate from a triangulation                              */
    gs = grid_interp_init(x, y, b, nbath, i_rule);

    for (c = 1; c <= params->ns2; c++) {
      if(isnan(params->bathy[c])) {
	params->bathy[c] = grid_interp_on_point(gs, params->x[c][0], params->y[c][0]);
	if (limf) {
	  if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	    params->bathy[c] = -params->bmax;
	  if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
	    params->bathy[c] = -params->bmin;
	  if (params->bmin && params->bathy[c] > params->bmin)
	    params->bathy[c] = LANDCELL;
	}
      }
      if (isnan(params->bathy[c])) params->bathy[c] = bmean;
      if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[c][0], params->y[c][0]);
      /*printf("%d %f %d %f\n",m,(double)c/(double)params->ns2,c,params->bathy[c]);*/
    }
    grid_specs_destroy(gs);
    params->nvals = params->ns2-1;
  }
  if (DEBUG("init_m"))
    dlog("init_m", "\nBathymetry read  from file OK\n");
}

/* END read_bathy_from_sparse_nc()                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured bathymetry from a netCDF file, packs into a      */
/* vector and interpolates onto an unstructured mesh.                */
/*-------------------------------------------------------------------*/
void read_bathy_from_sparse_bty(parameters_t *params, char *fname)
{
  FILE *bf;
  char buf[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nbath;
  double *x, *y, *b;
  double bmean;
  int n, i, j, cc, c;
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int bverbose = 0;

  /* Read the bathymetry values                                    */
  if ((bf = fopen(fname, "r")) == NULL)
    hd_quit("read_bathy_from_sparse_bty: Can't open bathymetry file %s.\n", fname);
  nbath = 0;
  while (fgets(buf, MAXSTRLEN, bf) != NULL) {
    nbath++;
  }
  rewind(bf);
  x = d_alloc_1d(nbath);
  y = d_alloc_1d(nbath);
  b = d_alloc_1d(nbath);
  n = 0;
  bmean = 0.0;
  while (fgets(buf, MAXSTRLEN, bf) != NULL) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    i = parseline(buf, fields, MAXNUMARGS);
    x[n] = atof(fields[0]);
    y[n] = atof(fields[1]);
    b[n] = atof(fields[2]);
    bmean += b[n];
    if (b[n] > 0) b[n] *= -1.0;
    if (bverbose) printf("%d %f %f : %f\n",n, x[n], y[n], b[n]);
    n++;
  }
  bmean /= (double)n;
  fclose(bf);

  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  sprintf(i_rule, "%c", '\0');
  if (!(prm_read_char(params->prmfd, "BATHY_INTERP_RULE", i_rule)))
    strcpy(i_rule, "linear");

  /* Interpolate from a triangulation                                */
  gs = grid_interp_init(x, y, b, nbath, i_rule);
  params->bathy = d_alloc_1d(params->ns2+1);
  for (c = 1; c <= params->ns2; c++) {
    params->bathy[c] = grid_interp_on_point(gs, params->x[c][0], params->y[c][0]);
    if (limf) {
      if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	params->bathy[c] = -params->bmax;
      if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
	params->bathy[c] = -params->bmin;
      if (params->bmin && params->bathy[c] > params->bmin)
	params->bathy[c] = LANDCELL;
    }
    if (isnan(params->bathy[c])) params->bathy[c] = bmean;
    if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[c][0], params->y[c][0]);
  }
  grid_specs_destroy(gs);
  params->nvals = params->ns2-1;

  if (DEBUG("init_m"))
    dlog("init_m", "\nBathymetry read  from file OK\n");
}

/* END read_bathy_from_sparse_bty()                                   */
/*-------------------------------------------------------------------*/

#define TF_USE  0x001
#define TF_TEST 0x002

/*-------------------------------------------------------------------*/
/* Interpolates bathymetry using unstructured method i_rule from the */
/* GRID_SPEC libraries.                                              */
/* Large input files may exhaust memory during the triangulation     */
/* process. This function sub-sections a file if its size id greater */
/* than ms points, and triangulates the sub-sections. If tf = 1 the  */
/* points triangulated are reduced further by only using those that  */
/* are coincident with the convex hull of the mesh. Sub-sections are */
/* overlapped ol points, and the convex hull is expanded by fact to  */
/* ensure that sufficient points encompass the mesh cells for        */
/* accurate interpolation.                                           */
/*-------------------------------------------------------------------*/
void bathy_interp_us(parameters_t *params, char *fname, char *i_rule, int mode)
{
  char buf[MAXSTRLEN];
  int n, m, i, j, cc, c;
  int is, ie, ni, ns;
  int fid;
  int ncerr;
  int lond, latd;
  int varid;
  int nbath;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  double *x, *y, *b, **botz, **cellx, **celly, *celx, *cely, *bz;
  double latx, latn, lonx, lonn;
  int ftype = 0;    /* 0 for netCDF, 1 for .bty                      */
  int intype = -1;  /* Bathymetry attributes type                    */
  double bmean;     /* Mean bathymetry                               */
  double sgn = 1.0; /* Sign multiplier for bathymetry                */
  int ms = 60e6;    /* Filesize limit for sub-sectioning (points)    */
  int bf = 1;       /* Omit bathy above msl (> 0)                    */
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int dm = 1;       /* Input decimation (not used)                   */
  int tf = 0;       /* Bounding box for input mesh                   */
  int bs = 1;       /* Sub-sectioning of bathymetry file             */
  int bverbose = 0; /* Print bathymetry output for debugging         */
  int sverbose = 0; /* Print control flow for debugging              */
  int ol = 5;       /* Overlap (points) for file sub-sections        */  
  double fact = 1.15; /* Expansion factor for mesh convex hull       */
  poly_t *pl;       /* Polygon of bathymetry convex hull             */
  poly_t *plb;      /* Polygon of mesh domain convex hull            */
  GRID_SPECS *gs = NULL;

  if (prm_read_char(params->prmfd, "BATHY_REPORT", buf))
    if (is_true(buf)) sverbose = 1;
  if (DEBUG("init_m"))
    dlog("init_m", "\nUnstructured bathymetry reading from file\n");
  hd_warn("Unstructured interpolation of bathy file %s using rule %s\n", fname, i_rule);
  if (sverbose) printf("\nUnstructured interpolation of bathy file %s using rule %s\n", fname, i_rule);

  /*-----------------------------------------------------------------*/
  /* Open the file and get dimensions                              */
  if (endswith(fname, ".nc")) {
    if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
      printf("Can't find input bathymetry file %s\n", fname);
      hd_quit((char *)nc_strerror(ncerr));
    }

    /* Get dimensions                                                  */
    i = 0;
    while (bathy_dims[i][0] != NULL) {
      if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &nce2) == 0 &&
	  nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &nce1) == 0 &&
	  ncw_var_id(fid, bathy_dims[i][0]) >= 0) {
	intype = i;
	break;
      }
      i++;
    }
    if (intype < 0)      
      hd_quit("bathy_interp_us: Can't find bathymetry attributes in file %s.\n", fname);
    
    varid = ncw_var_id(fid, bathy_dims[intype][3]);
    nc_inq_varndims(fid, varid, &lond);
    varid = ncw_var_id(fid, bathy_dims[intype][4]);
    nc_inq_varndims(fid, varid, &latd);
    if (sverbose) printf("Found attribues of %s, type = %d\n", fname, intype);

    /*-----------------------------------------------------------------*/
    /* Interpolation from unstructured input using a triangulation.    */
    /* The input bathymetry is structured, but is triangulated then    */
    /* interpolated using i_rule.                                      */
    /* Allocate and read.                                              */
    if (lond == 1)
      celx = d_alloc_1d(nce1);
    else if (lond == 2)
      cellx = d_alloc_2d(nce1, nce2);
    if (latd == 1)
      cely = d_alloc_1d(nce2);
    else if (latd == 2)
      celly = d_alloc_2d(nce1, nce2);
    botz = d_alloc_2d(nce1, nce2);
    start[0] = 0L;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = nce2;
    count[1] = nce1;
    count[2] = 0;
    count[3] = 0;
    if (lond == 1) {
      count[0] = nce1;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, celx);
    } else if (lond == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, cellx[0]);
    }
    if (sverbose) printf("Read longitude, dimension = %d\n", lond);
    if (latd == 1) {
      count[0] = nce2;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, cely);
    } else if (latd == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, celly[0]);
    }
    if (sverbose) printf("Read latitude, dimension = %d\n", lond);
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][0]), start, count, botz[0]);
    nc_close(fid);

    /*-----------------------------------------------------------------*/
    /* Set the wet bathymetry vector (to interpolate from)             */
    nbath = n = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) nbath++;
      }
    }
    if (sverbose) printf("Read bathy, %d points\n", nbath);
  }

  if (endswith(fname, ".bty") || endswith(fname, ".xyz")) {
    FILE *bf;
    char buf[MAXSTRLEN];
    int nn = 0, np = 0;
    ftype = 1;
    lond = latd = 2;
    nbath = n = 0;
    latn = lonn = HUGE;
    latx = lonx = -HUGE;
    if ((bf = fopen(fname, "r")) == NULL)
      hd_quit("bathy_interp_us: Can't open bathymetry file %s.\n", fname);
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      nbath++;
    }
    rewind(bf);
    cellx = d_alloc_2d(nbath, 1);
    celly = d_alloc_2d(nbath, 1);
    botz = d_alloc_2d(nbath, 1);
    n = 0;
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      i = parseline(buf, fields, MAXNUMARGS);
      cellx[0][n] = atof(fields[0]);
      celly[0][n] = atof(fields[1]);
      botz[0][n] = atof(fields[2]);
      if (botz[0][n] >= 0.0) np++;
      if (botz[0][n] < 0.0) nn++;
      latx = max(celly[0][n], latx);
      lonx = max(cellx[0][n], lonx);
      latn = min(celly[0][n], latn);
      lonn = min(cellx[0][n], lonn);
      n++;
    }
    fclose(bf);
    nce1 = nbath;
    nce2 = 1;
    n = 0;
    /* Bathy should be negative. If there are more positive bathys     */
    /* than negative, then assume the positive values are to be used.  */
    sgn = (np > nn) ? -1.0 : 1.0;
  }

  /*-------------------------------------------------------------------*/
  /* Check if the input file is large and needs subsectioning          */
  i = nce1 * nce2;
  if (i > ms) {
    bs = i /ms;
    hd_warn("Large bathy filesize: subdividing into %d subsections\n", bs);
    if (sverbose) printf("Large bathy filesize: subdividing into %d subsections\n", bs);
    /* Only interpolate within the bounds of the subsection (mode=1)   */
    /* within the bounds of the mesh (tf != 0).                        */
    mode = 1;
    tf |= TF_TEST;
  }
  /* Explicitly set the number of sections if required                 */
  if (prm_read_int(params->prmfd, "BATHY_SECTION", &bs)) {
    hd_warn("Large bathy filesize: subdividing into %d subsections\n", bs);
    if (sverbose) printf("Large bathy filesize: subdividing into %d subsections\n", bs);
    /* Only interpolate within the bounds of the subsection (mode=1)   */
    /* within the bounds of the mesh (tf != 0).                        */
    mode = 1;
    tf |= TF_TEST;
  }

  /*-------------------------------------------------------------------*/
  /* Get a convex hull of the model domain. Only use bathy data in     */
  /* this hull so as to reduce the size of the triangulation of bathy  */
  /* data.                                                             */
  if (prm_read_char(params->prmfd, "BATHY_TRUNCATE", buf))
    if (is_true(buf)) tf |= (TF_USE|TF_TEST);
  if (tf) {
    point *pin;
    point *pout;
    double xc = 0.0; 
    double yc = 0.0;

    /* Make a polygon of the convex hull of the input domain           */
    if (params->us_type & US_RS) {
      /* Input grid is a structured grid                               */
      pin = malloc(params->nce1 * params->nce2 * sizeof(point));
      cc = 0;
      for (j = 0; j < params->nce2; j++)
	for (i = 0; i < params->nce1; i++) {
	  if (!isnan(params->x[2*j][2*i]) && !isnan(params->y[2*j][2*i])) {
	    pin[cc].x = params->x[2*j][2*i];
	    pin[cc].y = params->y[2*j][2*i];
	    cc++;
	  }
	}
      n = cc;
    } else {
      /* Input grid is an unsructured mesh                             */
      pin = malloc(params->ns2 * sizeof(point));
      for (cc = 1; cc <= params->ns2; cc++) {
	pin[cc-1].x = params->x[cc][0];
	pin[cc-1].y = params->y[cc][0];
      }
      n = params->ns2;
    }
    pout = convex_hull(pin, &n);
    /* Get the centre of mass                                      */
    for (cc = 0; cc < n; cc++) {
      xc += pout[cc].x;
      yc += pout[cc].y;
    }
    xc /= (double)n;
    yc /= (double)n;
    /* Expand the convex hull by fact so as to encompas the        */
    /* perimeter cells (without this there may be issues on the    */
    /* boundaries with some interpolation schemes).                */
    plb = poly_create();
    for (cc = 0; cc < n; cc++) {
      double x = pout[cc].x - xc;
      double y = pout[cc].y - yc;
      double dist = fact * sqrt(x * x + y * y);
      double dir =  atan2(y, x);
      double sinth = sin(dir);
      double costh = cos(dir);
      x = xc + dist * costh;
      y = yc + dist * sinth;
      poly_add_point(plb, x, y);
    }
    poly_add_point(plb, pout[0].x, pout[0].y);
    free((point *)pout);
    if (sverbose) printf("Generated domain convex hull, %d points\n", plb->n);
  }

  /* Decimation - not used
  if (prm_read_int(params->prmfd, "BATHY_DECIMATE", &dm))
    hd_warn("Bathy file %s decimated by %d\n", fname, dm);
  */

  /*---------------------------------------------------------------*/
  /* Extract the bathymetry points used in the triangulation       */
  x = d_alloc_1d(nbath);
  y = d_alloc_1d(nbath);
  b = d_alloc_1d(nbath);
  is = ie = 0;
  ni = ceil(nce1 / bs);
  /* Loop over sections                                            */
  for (ns = 1; ns <= bs; ns++) {
    int fillf = 0;
    is = ie;
    ie = min(ns * ni, nce1);
    if (sverbose) {
      printf("--------------------------------------------\n");
      printf("Pass %d (i = %d to %d) of %d (interval = %d)\n", ns, is, ie, (int)nce1, ni);
    }
    n = 0;
    bmean = 0.0;
    /* Loop over points in this section                            */
    for (j = 0; j < nce2; j++) {
      int iss = max(is-ol, 0);
      int iee = min(ie+ol, nce1);
      for (i = iss; i < iee; i++) {
	double xin = (lond == 1) ? celx[i] : cellx[j][i];
	double yin = (latd == 1) ? cely[j] : celly[j][i];
	if (isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (tf) {
	  m = poly_contains_point(plb, xin, yin);
	  if (tf & TF_TEST && m) fillf = 1;
	  if (tf & TF_USE && !m) continue;
	} else
	  fillf = 1;
	/*if (dm > 1 && (j * nce1 + i)%dm != 0) continue;*/
	if (bf && sgn * botz[j][i] > 0.0) continue;
	if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) {
	  if (lond == 1)
	    x[n] = celx[i];
	  else if (lond == 2)
	    x[n] = cellx[j][i];
	  if (latd == 1)
	    y[n] = cely[j]; 
	  else if (latd == 2)
	    y[n] = celly[j][i]; 
	  b[n] = sgn * botz[j][i];
	  bmean += b[n];
	  n++;
	}
      }
    }
    if (!n || !fillf) {
      if (sverbose) printf("No bathy points found in mesh domain for this pass\n");
      continue;
    } else
      if (sverbose) printf("Converted bathy to 1D vector, %d points\n", n);

    if (n) bmean /= (double)n;
    if (bs == 1 || (bs > 1 && ns == bs)) {
      if (lond == 1)
	d_free_1d(celx);
      else if (lond == 2)
	d_free_2d(cellx);
      if (latd == 1)
	d_free_1d(cely);
      else if (latd == 2)
	d_free_2d(celly);
      d_free_2d(botz);
    }
    nbath = n;

    /*---------------------------------------------------------------*/
    /* Interpolate from a triangulation                              */
    gs = grid_interp_init(x, y, b, nbath, i_rule);
    if (sverbose) printf("Created GRID_SPEC using %s\n", i_rule);

    /*---------------------------------------------------------------*/
    /* If mode = 1 only interpolate within the bound of the supplied */
    /* file. Get a polygon of the perimeter of the bathy file.       */
    if (mode) {
      /* Get the bounding convex hull of the bathy points            */
      point *pin= malloc(nbath * sizeof(point));
      point *pout;
      for (n = 0; n < nbath; n++) {
	pin[n].x = x[n];
	pin[n].y = y[n];
      }
      n = nbath;
      pout = convex_hull(pin, &n);
      /* Make a polygon of the convex hull                           */
      pl = poly_create();
      for (i = 0; i < n; i++) {
	poly_add_point(pl, pout[i].x, pout[i].y);
      }
      poly_add_point(pl, pout[0].x, pout[0].y);
      free((point *)pin);
      free((point *)pout);
      if (sverbose) printf("Created convex hull of bathy file, %d points\n", pl->n);
    }

    /*---------------------------------------------------------------*/
    /* Do the interpolation                                          */
    if (params->us_type & US_RS) {
      n = 0;
      c = 0;
      for (j = 0; j < params->nce2; j++) {
	for (i = 0; i < params->nce1; i++) {
	  if(isnan(params->bathy[c])) {
	    if (isnan(params->x[2*j][2*i]) || isnan(params->y[2*j][2*i])) {
	      params->bathy[c] = -NOTVALID;
	    } else {
	      if (mode && !poly_contains_point(pl, params->x[2*j][2*i], params->y[2*j][2*i])) {
		c++;
		continue;
	      }
	      params->bathy[c] = -grid_interp_on_point(gs, params->x[2*j][2*i], params->y[2*j][2*i]);
	      if (limf) {
		if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
		  params->bathy[c] = params->bmax;
		if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
		  params->bathy[c] = params->bmin;
		/*
		if (params->bmin && params->bathy[c] > params->bmin)
		  params->bathy[c] = LANDCELL;
		*/
	      }
	      n++;
	      if (isnan(params->bathy[c])) params->bathy[c] = -bmean;
	      if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[2*j][2*i], params->y[2*j][2*i]);
	    }
	  }
	  c++;
	}
      }
      params->nvals = params->nce1 * params->nce2;
    } else {
      n = 0;
      for (c = 1; c <= params->ns2; c++) {
	if(isnan(params->bathy[c])) {
	  if (mode && !poly_contains_point(pl, params->x[c][0], params->y[c][0])) continue;
	  params->bathy[c] = grid_interp_on_point(gs, params->x[c][0], params->y[c][0]);
	  if (limf) {
	    if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	      params->bathy[c] = -params->bmax;
	    if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
	      params->bathy[c] = -params->bmin;
	    if (params->bmin && params->bathy[c] > params->bmin)
	      params->bathy[c] = LANDCELL;
	  }
	  n++;
	  if (isnan(params->bathy[c])) params->bathy[c] = bmean;
	  if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[c][0], params->y[c][0]);
	}
      }
      params->nvals = params->ns2-1;
    }
    if (sverbose) printf("Interpolted %d bathy points onto mesh\n", n);
    grid_specs_destroy(gs);
    if (mode) poly_destroy(pl);
  }

  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(b);
  if (tf) poly_destroy(plb);
}

/* END bathy_interp_us()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates bathymetry using unstructured method i_rule from the */
/* GRID_SPEC libraries.                                              */
/*-------------------------------------------------------------------*/
void bathy_interp_us_o(parameters_t *params, char *fname, char *i_rule, int mode)
{
  char buf[MAXSTRLEN];
  int n, m, i, j, cc, c;
  int intype = -1;
  int ftype = 0;
  int fid;
  int ncerr;
  int lond, latd;
  int varid;
  int nbath;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  double *x, *y, *b, **botz, **cellx, **celly, *celx, *cely, *bz;
  double latx, latn, lonx, lonn;
  double bmean;
  double sgn = 1.0;
  int bf = 1;       /* Omit bathy above msl (> 0)                    */
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  int dm = 1;       /* Input decimation */
  int tf = 0;
  GRID_SPECS *gs = NULL;
  int bverbose = 0;
  int sverbose = 1;
  poly_t *pl, *plb;

  if (DEBUG("init_m"))
    dlog("init_m", "\nUnstructured bathymetry reading from file\n");
  hd_warn("Unstructured interpolation of bathy file %s using rule %s\n", fname, i_rule);
  if (sverbose) printf("Unstructured interpolation of bathy file %s using rule %s\n", fname, i_rule);

  /*-----------------------------------------------------------------*/
  /* Open the file and get dimensions                              */
  if (endswith(fname, ".nc")) {
    if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
      printf("Can't find input bathymetry file %s\n", fname);
      hd_quit((char *)nc_strerror(ncerr));
    }

    /* Get dimensions                                                  */
    i = 0;
    while (bathy_dims[i][0] != NULL) {
      if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &nce2) == 0 &&
	  nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &nce1) == 0) {
	intype = i;
	break;
      }
      i++;
    }
    if (intype < 0)      
      hd_quit("bathy_interp_us: Can't find bathymetry attributes in file %s.\n", fname);
    
    varid = ncw_var_id(fid, bathy_dims[intype][3]);
    nc_inq_varndims(fid, varid, &lond);
    varid = ncw_var_id(fid, bathy_dims[intype][4]);
    nc_inq_varndims(fid, varid, &latd);
    if (sverbose) printf("Found attribues of %s, type = %d\n", fname, intype);

    /*-----------------------------------------------------------------*/
    /* Interpolation from unstructured input using a triangulation.    */
    /* The input bathymetry is structured, but is triangulated then    */
    /* interpolated using i_rule.                                      */
    /* Allocate and read.                                              */
    if (lond == 1)
      celx = d_alloc_1d(nce1);
    else if (lond == 2)
      cellx = d_alloc_2d(nce1, nce2);
    if (latd == 1)
      cely = d_alloc_1d(nce2);
    else if (latd == 2)
      celly = d_alloc_2d(nce1, nce2);
    botz = d_alloc_2d(nce1, nce2);
    start[0] = 0L;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = nce2;
    count[1] = nce1;
    count[2] = 0;
    count[3] = 0;
    if (lond == 1) {
      count[0] = nce1;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, celx);
    } else if (lond == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][3]), start, count, cellx[0]);
    }
    if (sverbose) printf("Read longitude, dimension = %d\n", lond);
    if (latd == 1) {
      count[0] = nce2;
      count[1] = 0;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, cely);
    } else if (latd == 2) {
      count[0] = nce2;
      count[1] = nce1;
      nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][4]), start, count, celly[0]);
    }
    if (sverbose) printf("Read latitude, dimension = %d\n", lond);
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, bathy_dims[intype][0]), start, count, botz[0]);
    nc_close(fid);

    /*-----------------------------------------------------------------*/
    /* Set the wet bathymetry vector (to interpolate from)             */
    nbath = n = 0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) nbath++;
      }
    }
    if (sverbose) printf("Read bathy, %d points\n", nbath);
  }

  if (endswith(fname, ".bty") || endswith(fname, ".xyz")) {
    FILE *bf;
    char buf[MAXSTRLEN];
    int nn = 0, np = 0;
    ftype = 1;
    lond = latd = 2;
    nbath = n = 0;
    latn = lonn = HUGE;
    latx = lonx = -HUGE;
    if ((bf = fopen(fname, "r")) == NULL)
      hd_quit("bathy_interp_us: Can't open bathymetry file %s.\n", fname);
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      nbath++;
    }
    rewind(bf);
    cellx = d_alloc_2d(nbath, 1);
    celly = d_alloc_2d(nbath, 1);
    botz = d_alloc_2d(nbath, 1);
    n = 0;
    while (fgets(buf, MAXSTRLEN, bf) != NULL) {
      char *fields[MAXSTRLEN * MAXNUMARGS];
      i = parseline(buf, fields, MAXNUMARGS);
      cellx[0][n] = atof(fields[0]);
      celly[0][n] = atof(fields[1]);
      botz[0][n] = atof(fields[2]);
      if (botz[0][n] >= 0.0) np++;
      if (botz[0][n] < 0.0) nn++;
      latx = max(celly[0][n], latx);
      lonx = max(cellx[0][n], lonx);
      latn = min(celly[0][n], latn);
      lonn = min(cellx[0][n], lonn);
      n++;
    }
    fclose(bf);
    nce1 = nbath;
    nce2 = 1;
    n = 0;
    /* Bathy should be negative. If there are more positive bathys     */
    /* than negative, then assume the positive values are to be used.  */
    sgn = (np > nn) ? -1.0 : 1.0;
  }

  /* Get a convex hull of the model domain. Only use bathy data in     */
  /* this hull so as to reduce the size of the triangulation of bathy  */
  /* data.                                                             */
  if (prm_read_char(params->prmfd, "BATHY_TRUNCATE", buf))
    tf = is_true(buf);
  if (tf) {
    point *pin= malloc(params->ns2 * sizeof(point));
    point *pout;
    double xc = 0.0; 
    double yc = 0.0;
    double fact = 1.02;
    for (cc = 1; cc <= params->ns2; cc++) {
      pin[cc-1].x = params->x[cc][0];
      pin[cc-1].y = params->y[cc][0];
    }
    n = params->ns2;
    pout = convex_hull(pin, &n);
    /* Make a polygon of the convex hull                           */
    plb = poly_create();
    for (i = 0; i < n; i++) {
      poly_add_point(plb, pout[i].x, pout[i].y);
    }
    poly_add_point(plb, pout[0].x, pout[0].y);
    free((point *)pin);
    free((point *)pout);

    /* Get the centre of mass                                      */
    for (cc = 0; cc < plb->n; cc++) {
      xc += plb->x[cc];
      yc += plb->y[cc];
    }
    xc /= (double)plb->n;
    yc /= (double)plb->n;
    /* Expand the convex hull by fact so as to encompas the        */
    /* perimeter cells (without this there may be issues on the    */
    /* boundaries with some interpolation schemes).                */
    for (cc = 0; cc < plb->n; cc++) {
      double x = plb->x[cc] - xc;
      double y = plb->y[cc] - yc;
      double dist = fact * sqrt(x * x + y * y);
      double dir =  atan2(y, x);
      double sinth = sin(dir);
      double costh = cos(dir);
      plb->x[cc] = xc + dist * costh;
      plb->y[cc] = yc + dist * sinth;
    }
    if (sverbose) printf("Generated domain convex hull, %d points\n", plb->n);
  }

  /* Decimation                                                        */
  if (prm_read_int(params->prmfd, "BATHY_DECIMATE", &dm))
    hd_warn("Bathy file %s decimated by %d\n", fname, dm);

  x = d_alloc_1d(nbath);
  y = d_alloc_1d(nbath);
  b = d_alloc_1d(nbath);
  bmean = 0.0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      double xin = (lond == 1) ? celx[i] : cellx[j][i];
      double yin = (latd == 1) ? cely[j] : celly[j][i];
      if (isnan(botz[j][i])) continue;
      if (lond == 1 && isnan(celx[i])) continue;
      if (lond == 2 && isnan(cellx[j][i])) continue;
      if (latd == 1 && isnan(cely[j])) continue;
      if (latd == 2 && isnan(celly[j][i])) continue;
      if (tf && !poly_contains_point(plb, xin, yin)) continue;
      if (dm > 1 && (j * nce1 + i)%dm != 0) continue;
      if (bf && sgn * botz[j][i] > 0.0) continue;
      if (botz[j][i] != LANDCELL && fabs(botz[j][i]) != fabs(NOTVALID)) {
	if (lond == 1)
	  x[n] = celx[i];
	else if (lond == 2)
	  x[n] = cellx[j][i];
	if (latd == 1)
	  y[n] = cely[j]; 
	else if (latd == 2)
	    y[n] = celly[j][i]; 
	b[n] = sgn * botz[j][i];
	bmean += b[n];
	n++;
      }
    }
  }
  if (sverbose) printf("Converted bathy to 1D vector, %d points\n", n);

  if (n) bmean /= (double)n;
  if (lond == 1)
    d_free_1d(celx);
  else if (lond == 2)
    d_free_2d(cellx);
  if (latd == 1)
    d_free_1d(cely);
  else if (latd == 2)
    d_free_2d(celly);
  d_free_2d(botz);
  nbath = n;

  /* Interpolate from a triangulation                                */
  gs = grid_interp_init(x, y, b, nbath, i_rule);
  if (sverbose) printf("Created GRID_SPEC using %s\n", i_rule);

  /* If mode = 1 only interpolate within the bound of the supplied   */
  /* file. Get a polgon of the perimeter of the bathy file.          */
  if (mode) {
    if (ftype == 0) {
      timeseries_t *ts = NULL;
      ts = (timeseries_t *)malloc(sizeof(timeseries_t));
      if (ts == NULL)
	hd_quit("read_bathy_from_sparse_nc: No memory available.\n");
      memset(ts, 0, sizeof(timeseries_t));
      ts_read(fname, ts);
      pl = nc2poly(fid, nce1, nce2, bathy_dims[intype][3], bathy_dims[intype][4], fname, ts);
      ts_free((timeseries_t*)ts);
      free(ts);
      if (sverbose) printf("Created poly of bathy file bounds, %d points\n", pl->n);
    }
    if (ftype == 1) {
      /* Get the bounding convex hull of the bathy points            */
      point *pin= malloc(nbath * sizeof(point));
      point *pout;
      for (n = 0; n < nbath; n++) {
	pin[n].x = x[n];
	pin[n].y = y[n];
      }
      n = nbath;
      pout = convex_hull(pin, &n);
      /* Make a polygon of the convex hull                           */
      pl = poly_create();
      for (i = 0; i < n; i++) {
	poly_add_point(pl, pout[i].x, pout[i].y);
      }
      poly_add_point(pl, pout[0].x, pout[0].y);
      free((point *)pin);
      free((point *)pout);
      if (sverbose) printf("Created convex hull of bathy file, %d points\n", pl->n);
      /*
      poly_add_point(pl, lonx, latx);
      poly_add_point(pl, lonx, latn);
      poly_add_point(pl, lonn, latn);
      poly_add_point(pl, lonn, latx);
      */
    }
  }

  n = 0;
  for (c = 1; c <= params->ns2; c++) {
    if(isnan(params->bathy[c])) {
      if (mode && !poly_contains_point(pl, params->x[c][0], params->y[c][0])) continue;
      params->bathy[c] = grid_interp_on_point(gs, params->x[c][0], params->y[c][0]);
      if (limf) {
	if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	  params->bathy[c] = -params->bmax;
	if (params->bmin && fabs(params->bathy[c]) < params->bmin) 
	  params->bathy[c] = -params->bmin;
	if (params->bmin && params->bathy[c] > params->bmin)
	  params->bathy[c] = LANDCELL;
      }
    }
    if (isnan(params->bathy[c])) params->bathy[c] = bmean;
    if (bverbose) printf("%d %f : %f %f\n",c, params->bathy[c], params->x[c][0], params->y[c][0]);
    n++;
    /*printf("%d %f %d %f\n",m,(double)c/(double)params->ns2,c,params->bathy[c]);*/
  }
  if (sverbose) printf("Interpolted %d bathy points onto mesh\n", n);

  grid_specs_destroy(gs);
  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(b);
  if (mode) poly_destroy(pl);
  if (tf) poly_destroy(plb);
  params->nvals = params->ns2-1;
}

/* END bathy_interp_us()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates bathymetry using bilinear methods via timeseries     */
/* libraries.                                                        */
/*-------------------------------------------------------------------*/
void bathy_interp_s(parameters_t *params, char *fname, int mode)
{
  char buf[MAXSTRLEN];
  int n, m, i, j, cc, c;
  int intype = -1;
  int fid;
  int ncerr;
  int lond, latd;
  int varid;
  int nbath;
  int idb;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  double *x, *y, *b, **botz, **cellx, **celly, *celx, *cely;
  double bmean;
  int limf = 1;     /* Limit bathymetry to bmin and bmax             */
  timeseries_t *ts = NULL;
  poly_t *pl;
  int bverbose = 0;

  hd_warn("Structured interpolation of bathy file %s\n", fname);

  /* Open the file and get dimensions                                */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input bathymetry file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }
  i = 0;
  while (bathy_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, bathy_dims[i][1]), &nce1) == 0) {
      intype = i;
      break;
    }
    i++;
  }
  if (intype < 0)      
    hd_quit("read_bathy_from_sparse_nc: Can't find bathymetry attributes in file %s.\n", fname);

  /* Initialize the timeseries file                                  */
  ts = (timeseries_t *)malloc(sizeof(timeseries_t));
  if (ts == NULL)
    hd_quit("read_bathy_from_sparse_nc: No memory available.\n");
  memset(ts, 0, sizeof(timeseries_t));
    
  /* Read the time series                                            */
  ts_read(fname, ts);
  if ((idb = ts_get_index(ts, fv_get_varname(fname, bathy_dims[intype][0], buf))) == -1)
    hd_quit("read_bathy_from_sparse_nc: Can't find variable %s in file %s\n", 
	    bathy_dims[intype][0], fname);

  /* Get a polgon of the perimeter of the bathy file                 */
  pl = nc2poly(fid, nce1, nce2, bathy_dims[intype][3], bathy_dims[intype][4], fname, ts);

  /* Read the values                                                 */
  if (params->us_type & (US_RUS|US_G) && !(params->us_type & US_RS)) {
    /* Unstructured.                                                 */
    for (cc = 1; cc <= params->ns2; cc++) {
      if (isnan(params->bathy[cc])) {
	/* Note: ts_eval_xy() for bathymetry files (idb = 'botz'     */
	/* or 'height') will return NaN if the xytoij() mapping does */
	/* not return valid indices (i.e. (x,y) lies outside the     */
	/* (i,j) bounds of the file).                                */
	/* This seems not entirely reliable though, so we use a      */
	/* bounding polygon for the perimeter.                       */
	/* Check if the location lies within the file's bounding     */
	/* perimeter.                                                */
	if (!poly_contains_point(pl, params->x[cc][0], params->y[cc][0])) continue;

	params->bathy[cc] = ts_eval_xy(ts, idb, 0.0, params->x[cc][0], params->y[cc][0]);

	if (limf && !isnan(params->bathy[cc])) {
	  if (params->bmax && fabs(params->bathy[c]) > params->bmax) 
	    params->bathy[cc] = -params->bmax;
	  if (params->bmin && fabs(params->bathy[cc]) < params->bmin) 
	    params->bathy[c] = -params->bmin;
	  if (params->bmin && params->bathy[cc] > params->bmin) {
	    if (mode > 1)
	      params->bathy[cc] = NaN;
	    else
	      params->bathy[cc] = LANDCELL;
	  }
	}
      }
    }
    params->nvals = cc-1;
  } else {
    /* Structured                                                    */
    cc = 0;
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++) {
	if (!isnan(params->x[2*j][2*i]) && !isnan(params->y[2*j][2*i])) {

	  if (!poly_contains_point(pl, params->x[2*j][2*i], params->y[2*j][2*i])) continue;
	  params->bathy[cc] = -ts_eval_xy(ts, idb, 0.0, params->x[2*j][2*i], params->y[2*j][2*i]);
	} else
	  params->bathy[cc] = NOTVALID;
	cc++;
      }
    params->nvals = params->nce1 * params->nce2;
  }
  ts_free((timeseries_t*)ts);
  free(ts);
  poly_destroy(pl);
  nc_close(fid);
}

/* END bathy_interp_s()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract bathymetry from the database bathy file        */
/*-------------------------------------------------------------------*/
static void mask_bathy_with_coast_us(const char *cfname, parameters_t *params) {

  coast_file_t* cf = coast_open(cfname);
  if (cf != NULL) {
    double minlon = 360, maxlon = -360;
    double minlat = 90, maxlat = -90;
    int cc, c, j, m;
    coast_tile_iterator_t *ti;
    coast_tile_t *tile;
    int ns2 = params->ns2;
    int filef = 1;
    FILE *fp;

    /* Find limits of the domain                                     */
    for (cc = 1; cc <= ns2; cc++) {
      for (j = 0; j <= params->npe2[cc]; j++) {
	if (params->x[cc][j] < minlon) minlon = params->x[cc][j];
	if (params->x[cc][j] > maxlon) maxlon = params->x[cc][j];
	if (params->y[cc][j] < minlat) minlat = params->y[cc][j];
	if (params->y[cc][j] > maxlat) maxlat = params->y[cc][j];
      }
    }
    if(filef) {
      /*printf("coast minmax %f %f %f %f\n",minlon,maxlon,minlat,maxlat);*/
      fp = fopen("coast.cst", "w");
    }
    
    /* Load a tile at a time, and compare the lat/lons               */
    ti = coast_get_tile_iterator(cf, minlon, maxlon, minlat, maxlat);
    while((tile = ti->next_tile(ti)) != NULL) {
      for (m=0; m<tile->npolys; ++m) {
        poly_t *p = tile->polys[m];
        if (tile->level[m] != 1) continue;
	if(filef) {
	  for (cc = 0; cc < p->n; cc++) {
	    fprintf(fp,"%f %f\n", p->x[cc], p->y[cc]);
	  }
	  fprintf(fp, "NaN NaN\n");
	}
	for (cc = 1; cc <= ns2; cc++) {
	  int land = 0;
	  int npe = params->npe2[cc];
	  int lmin = npe + 1;
	  if (params->bathy[cc] == LANDCELL) continue;

	  /* Check the cell centre and vertices of the polygon.      */
	  /* If more than more than 4 cells land, then assume the    */
	  /* cell is LAND.                                           */
	  for (j = 0; j <= npe; j++) {
	    double nx = params->x[cc][j];
	    double ny = params->y[cc][j];
	    if (poly_contains_point(p, nx, ny))
	      ++land;
	  }
	  if (land > lmin / 2)
	    params->bathy[cc] = LANDCELL;
	}
      }
      ti->free_tile(tile);
    }
    if(filef) fclose(fp);
    coast_free_tile_iterator(ti);
    coast_close(cf);
  }
}

static void mask_bathy_with_coast_s(const char *cfname, parameters_t *params) {

  coast_file_t* cf = coast_open(cfname);
  if (cf != NULL) {
    double minlon = 360, maxlon = -360;
    double minlat = 90, maxlat = -90;
    double **x, **y;
    int has_proj = (strlen(params->projection) > 0);
    int is_geog = has_proj && (strcasecmp(params->projection, GEOGRAPHIC_TAG) == 0);
    int i, j, n, m;
    int is, js, ie, je, jm, im;
    int nce1 = params->nce1;
    int nce2 = params->nce2;
    coast_tile_iterator_t *ti;
    coast_tile_t *tile;

    /* Adjust for grid edges                                         */
    is = js = 0;
    ie = nce1;
    je = nce2;

    /* Reduce double grid to single grid corners                     */
    x = d_alloc_2d(nce1+1, nce2+1);
    y = d_alloc_2d(nce1+1, nce2+1);

    if (is_geog) {
      for (j=js; j<je+1; j++) {
        for (i=is; i<ie+1; i++) {
          im = i;
          jm = j;
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
          im = i;
          jm = j;
          mp_inverse(mp, params->x[2*j][2*i], params->y[2*j][2*i],
                 &y[jm][im], &x[jm][im]);
        }
      }
    }


    /* Read the coastline data if available                          */

    for (j=0; j<nce2; j++) {
      for (i=0; i<nce1; i++, ++n) {
        if (x[j][i] < minlon) minlon = x[j][i];
        if (x[j][i] > maxlon) maxlon = x[j][i];
        if (y[j][i] < minlat) minlat = y[j][i];
        if (y[j][i] > maxlat) maxlat = y[j][i];
      }
    }

    /* Load a tile at a time, and compare the lat/lons               */
    ti = coast_get_tile_iterator(cf, minlon, maxlon, minlat, maxlat);
    while((tile = ti->next_tile(ti)) != NULL) {
      for (m=0; m<tile->npolys; ++m) {
        poly_t *p = tile->polys[m];
        if (tile->level[m] != 1) continue;
        n = 0;
        for (j=0; j<nce2; j++) {
          for (i=0; i<nce1; i++, ++n) {
            int ii, jj, land = 0;
            if (params->bathy[n] == -LANDCELL) continue;

            /* Build a 9 x 9 grid. If more than more than 4 cells land,
             * then assume the cell is LAND.                         */
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
/* Note inverse barometer is added to the initial condition in       */
/* pre_run_setup() because patm and dens aren't initialised yet.     */
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
  char i_rule[MAXSTRLEN];
  double val;
  char *fields[MAXSTRLEN * MAXNUMARGS];

  if (strlen(params->eta_init)) {

    strcpy(buf, params->eta_init);
    if (sscanf(buf,"%s %d %d", key, &i, &j) == 3) {
      if (strcmp(key, "BARO_VORTEX") == 0)
	baro_vortex(master, geom, i, j);
      return;
    }

    sprintf(i_rule, "%c", '\0');
    if ((i = parseline(buf, fields, MAXNUMARGS)) == 2) {
      strcpy(i_rule, fields[1]);
    }

    value_init_2d(master, master->eta, params->prmfd, params->eta_init,
		  "eta", "SURFACE", 0.0, i_rule); 

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
      if (fscanf(params->prmfd, "%d : %lf", &i, &val) != 2)
	hd_quit("eta_init: Can't read cc : val in SURFACE_POINTS list.\n", n);
      /* Flush the remainder of the line */
      prm_flush_line(params->prmfd);
      c = geom->cc2s[i];
      if (c > 0 && c < geom->szcS)
	master->eta[c] = val;
    }
  }

  /* Scale                                                           */
  if ((val = get_scaling(params->scale_v, "ETA")) != 1.0) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      master->eta[c] += val;
    }
  }

  /* Check for NaNs                                                  */
  if (params->us_type & US_IJ) {
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
  int e, ee, es, c;
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
     for (ee = 1; ee <= geom->b3_e1; ee++) {
	e = geom->w3_e1[ee];
	es = geom->m2de[e];
	c = geom->e2c[e][0];
	if (geom->wgst[c]) c = geom->e2c[e][1];
	thetau1 = geom->thetau1[es];
	u = ts_eval_xyz(ts, idu, master->t, geom->u1x[es], geom->u1y[es],
			geom->cellz[c]);
	v = ts_eval_xyz(ts, idv, master->t, geom->u1x[es], geom->u1y[es],
			geom->cellz[c]);
	master->u1[e] = cos(thetau1) * u + sin(thetau1) * v;
      }
    } else if (tu1 == 2 && tu2 == 2) {
      for (ee = 1; ee <= geom->b3_e1; ee++) {
	e = geom->w3_e1[ee];
	es = geom->m2d[e];
	c = geom->e2c[e][0];
	if (geom->wgst[c]) c = geom->e2c[e][1];
	master->u1[e] = ts_eval_xyz(ts, idu, master->t, geom->u1x[es], geom->u1y[es],
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

  params->numbers = params->numbers1 = 0;
  sprintf(keyword, "NUMBERS");
  if (prm_read_char(params->prmfd, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->numbers = params->numbers1 = NONE;
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
	ntr++;
      }
      if (contains_token(buf, "SPEED_2D") != NULL) {
	params->numbers |= SPEED_2D;
	params->ntrS++;
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
	params->ntrS+=1;
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
      if (contains_token(buf, "PASSIVE") != NULL) {
	params->numbers |= PASS;
	ntr++;
      }
      if (contains_token(buf, "GLIDER") != NULL) {
	params->numbers |= GLIDER;
	ntr++;
      }
      if (contains_token(buf, "RESOLUTION") != NULL) {
	params->numbers |= CELLRES;
	params->ntrS+=3;
      }
      if (contains_token(buf, "AREA") != NULL) {
	params->numbers1 |= CELLAREA;
	params->ntrS+=2;
      }
      if (contains_token(buf, "EKMAN_PUMP") != NULL) {
	params->numbers |= EKPUMP;
	params->ntrS+=1;
      }
      if (contains_token(buf, "TIDE_FRONT") != NULL) {
	params->numbers |= TIDEFR;
	params->ntrS+=1;
      }
      if (contains_token(buf, "U1VH") != NULL) {
	params->numbers1 |= U1VHC;
	ntr++;
      }
      if (contains_token(buf, "CONTINUITY") != NULL) {
	params->numbers1 |= VOLCONT;
	ntr++;
      }
      if (contains_token(buf, "CELL_INDEX") != NULL) {
	params->numbers1 |= CENTI;
	ntr++;
      }
      if (contains_token(buf, "UNIFORMITY") != NULL) {
	params->numbers1 |= MESHUN;
	params->ntrS+=1;
      }
      if (contains_token(buf, "TPXO") != NULL) {
	params->numbers1 |= TPXO;
	params->tidef |= TD_TPXOE;
	params->ntrS++;
      }
      if (contains_token(buf, "TPXO_VEL") != NULL) {
	params->numbers1 |= TPXOV;
	params->tidef |= TD_TPXOV;
	params->ntrS+=2;
      }
      if (contains_token(buf, "TPXO_TRAN") != NULL) {
	params->numbers1 |= TPXOT;
	params->ntrS+=2;
      }
      if (contains_token(buf, "TRAN2D") != NULL) {
	params->numbers1 |= TRAN2D;
	params->ntrS+=2;
      }
      if (contains_token(buf, "WIND") != NULL) {
	params->numbers1 |= WINDSPDI;
	params->ntrS+=2;
      }
      if (contains_token(buf, "ALL_NUMBERS") != NULL) {
	params->numbers |= (BRUNT|INT_WAVE|RICHARD_GR|RICHARD_FL|REYNOLDS|
			    FROUDE|SOUND|ROSSBY_IN|ROSSBY_EX|SHEAR_V|
			    BUOY_PROD|SHEAR_PROD|SPEED_3D|SPEED_2D|
			    OBC_PHASE|SPEED_SQ|SIGMA_T|ENERGY|WET_CELLS|
			    SLOPE|SURF_LAYER|BOTSTRESS|UNIT|EKPUMP|TIDEFR|
			    CELLRES|WINDSPDI);
	params->numbers1 |= (U1VHC);
	ntr+=18;
	params->ntrS += 17;
      }
      if (params->numbers == 0) params->numbers = NONE;
      if (params->numbers1 == 0) params->numbers1 = NONE;
    }
  } else {
    params->numbers = params->numbers1 = NONE;
  }
  return(ntr);
}

/* END numbers_init()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads automated 2D or 3D file import                              */
/*-------------------------------------------------------------------*/
int import_init(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n;

  sprintf(keyword, "IMPORT2D");
  if (prm_read_char(fp, keyword, buf)) {
    n = parseline(buf, fields, MAXNUMARGS);
    strcpy(params->imp2dn, fields[0]);
    strcpy(params->imp2du, fields[1]);
    strcpy(params->imp2df, fields[2]);
    if (n > 3) {
      sprintf(params->imp2dt, "%s%s", fields[3], fields[4]);
    }
    params->ntrS += 1;
  }
  sprintf(keyword, "IMPORT3D");
  if (prm_read_char(fp, keyword, buf)) {
    n = parseline(buf, fields, MAXNUMARGS);
    strcpy(params->imp3dn, fields[0]);
    strcpy(params->imp3du, fields[1]);
    strcpy(params->imp3df, fields[2]);
    if (n > 3) {
      sprintf(params->imp3dt, "%s%s", fields[3], fields[4]);
    }
    return(1);
  }
  return(0);
}

/* END import_init()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initializes the degree heating diagnostic                         */
/*-------------------------------------------------------------------*/
int read_dhw(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN], keyword[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n, m, ret = 0;

  if (prm_read_int(fp, "NDHW", &params->ndhw)) {
    
    params->dhw = (char **)malloc(params->ndhw * sizeof(char *));
    params->dhdf = (char **)malloc(params->ndhw * sizeof(char *));
    params->dhwt = (char **)malloc(params->ndhw * sizeof(char *));
    params->dhwf = i_alloc_1d(params->ndhw);
    params->dhw_dt = d_alloc_1d(params->ndhw);
    params->dhwh = d_alloc_1d(params->ndhw);

    for (m = 0; m < params->ndhw; m++) {
      params->dhw[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      params->dhdf[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      params->dhwt[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      params->dhwf[m] = 0;
      sprintf(params->dhw[m], "%c", '\0');
      sprintf(params->dhwt[m], "%c", '\0');
      
      sprintf(keyword, "DHW%d", m);
      if (prm_read_char(fp, keyword, buf)) {
	sprintf(keyword, "DHW%d.text", m);
	prm_read_char(fp, keyword, params->dhwt[m]);
	n = parseline(buf, fields, MAXNUMARGS);
	if (n == 1) {
	  strcpy(params->dhw[m], fields[0]);
	  ret += 3;
	  params->dhwf[m] = DHW_RT;
	}
	if (n >= 2) {
	  strcpy(params->dhw[m], fields[0]);
	  strcpy(params->dhdf[m], fields[1]);
	  sprintf(keyword, "DHW_DT%d", m);
	  if (prm_read_char(fp, keyword, buf))
	    tm_scale_to_secs(buf, &params->dhw_dt[m]);
	  else
	    params->dhw_dt[m] = 86400.0;
	  params->dhwf[m] = DHW_NOAA;
	  if (n == 3) {
	    if (strcmp(fields[2], "mean") == 0 || strcmp(fields[2], "MEAN") == 0) 
	      params->dhwf[m] |= DHW_MEAN;
	    else {
	      params->dhwf[m] |= DHW_SNAP;
	      params->dhwh[m] = atof(fields[2]);
	    }
	  } else
	    params->dhwf[m] |= DHW_INT;
	  ret += 3;
	}
      }
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
  int n1, ns2 = 0;
  int doread = 0;

  prm_set_errfn(hd_quit);
  params->z0s = NULL;
  prm_read_double(params->prmfd, "Z0", &params->z0);
  if (params->us_type & US_IJ) {
    if ((int)params->z0 == params->nce1 * params->nce2)
      doread = 1;
  } else {
    if (prm_read_int(params->prmfd, "nMesh2_face", &ns2))
      if ((int)params->z0 == ns2) doread = 1;
  }

  /* Check if it is a list of values in the parameter file          */
  if (doread) {
    if (prm_read_darray(params->prmfd, "Z0", &params->z0s, &n1) > 0) {
      if (n1 != params->ns2)
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
      /* Check if it is a list of values in the parameter file       */
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
/* Reads a 2D structured array                                       */
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
/* Creates a polygon of the perimeter of data within a netCDF file   */
/*-------------------------------------------------------------------*/
poly_t *nc2poly(int fid, int nce1, int nce2, char *xname, char *yname, char *bname, timeseries_t *ts)
{ 
  int i, j;                                     /* Cartesian couters */
  size_t start[4];
  size_t count[4];
  double **lon, **lat;
  poly_t *pl;
  int n;
  df_variable_t *vlon, *vlat;
  pl = poly_create();
  datafile_t *df = ts->df;

  /* Get variables                                                   */
  if ((vlon = df_get_variable_by_name(df, xname)) == NULL)
    hd_quit("nc2poly: Can't find variable %s in file %s\n", xname, bname);
  if ((vlat = df_get_variable_by_name(df, yname)) == NULL)
    hd_quit("nc2poly: Can't find variable %s in file %s\n", yname, bname);

  /* Allocate                                                        */
  lon = d_alloc_2d(nce1, nce2);
  lat = d_alloc_2d(nce1, nce2);

  /* Read in the coordinate information                              */
  if (vlon->nd == 1) {
    double *x = d_alloc_1d(nce1);
    start[0] = 0;
    count[0] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, xname), start, count, x);
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
	lon[j][i] = x[i];
    d_free_1d(x);
  } else if (vlon->nd == 2) {
    start[0] = 0;
    start[1] = 0;
    count[0] = nce2;
    count[0] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, xname), start, count, lon[0]);
  }
  if (vlat->nd == 1) {
    double *y = d_alloc_1d(nce2);
    start[0] = 0;
    count[0] = nce2;
    nc_get_vara_double(fid, ncw_var_id(fid, yname), start, count, y);
    for (j = 0; j < nce2; j++)
      for (i = 0; i < nce1; i++)
	lat[j][i] = y[j];
    d_free_1d(y);
  } else if (vlon->nd == 2) {
    start[0] = 0;
    start[1] = 0;
    count[0] = nce2;
    count[0] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, xname), start, count, lat[0]);
  }

  /* Create the polygon from the perimeter                           */
  i = 0;
  for (j = 0; j < nce2; j++) {
    poly_add_point(pl, lon[j][0], lat[j][0]);
  }
  for (i = 0; i < nce1; i++) {
    poly_add_point(pl, lon[nce2-1][i], lat[nce2-1][i]);
  }
  for (j = nce2-1; j >= 0; j--) {
    poly_add_point(pl, lon[j][nce1-1], lat[j][nce1-1]);
  }
  for (i = nce1-1; i >= 0; i--) {
    poly_add_point(pl, lon[0][i], lat[0][i]);
  }
  poly_add_point(pl, lon[0][0], lat[0][0]);

  d_free_2d(lon);
  d_free_2d(lat);
  return(pl);
}

/* END nc2poly()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes parameter structure specification to file                  */
/*-------------------------------------------------------------------*/
void params_write(parameters_t *params, dump_data_t *dumpdata, char *name)
{
  FILE *op, *fopen();
  int n, tn, i, j, ntr, atr;
  double d1;
  char tag[MAXSTRLEN];
  char key[MAXSTRLEN];
  char bname[MAXSTRLEN];

  sprintf(tag, "%s.prm", name);
  op = fopen(tag, "w");

  fprintf(op, "# COMPAS parameter file\n");
  fprintf(op, "CODEHEADER           %s\n", params->codeheader);
  fprintf(op, "PARAMETERHEADER      %s\n", params->parameterheader);
  if (strlen(params->reference))
    fprintf(op, "REFERENCE            %s\n", params->reference);
  if (strlen(params->notes))
    fprintf(op, "NOTES                %s\n", params->notes);
  fprintf(op, "DESCRIPTION          %s\n", params->grid_desc);
  fprintf(op, "NAME                 %s\n", params->grid_name);
  if (params->runno > 0)
    fprintf(op, "ID_NUMBER            %s\n", params->runnoc);
  else
    fprintf(op, "ID_NUMBER            1.0\n");
  if (strlen(params->runcode))
    fprintf(op, "ID_CODE              %s\n", params->runcode);
  else
    fprintf(op, "ID_CODE              %s|G1.0|H0.0|S0.0|B0.0\n", params->grid_name);
  if (strlen(params->rev))
  fprintf(op, "REVISION             %s\n", params->rev);
  if (strlen(params->projection))
    fprintf(op, "PROJECTION           %s\n", params->projection);
  fprintf(op, "TIMEUNIT             %s\n", params->timeunit);
  fprintf(op, "OUTPUT_TIMEUNIT      %s\n", params->output_tunit);
  fprintf(op, "LENUNIT              %s\n", params->lenunit);
  fprintf(op, "START_TIME           %s\n", params->start_time);
  fprintf(op, "STOP_TIME            %s\n", params->stop_time);
  if (!(params->history & NONE)) {
    fprintf(op, "HISTORY             ");
    if (params->history & HST_LOG) fprintf(op, " LOG");
    if (params->history & HST_DIF) fprintf(op, " DIF");
    if (params->history & HST_RESET) fprintf(op, " RESET");
    if (params->history & HST_NOTES) fprintf(op, " NOTES");
    fprintf(op, "\n");
  }
  if (strlen(params->trl))
    fprintf(op, "TECHNOLOGY_READINESS_LEVEL %s\n", params->trl);
  else
    fprintf(op, "TECHNOLOGY_READINESS_LEVEL TR6\n");
  fprintf(op, "\n");
  
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
    fprintf(op,"OutputFiles          %d\n\n",params->ndf);
    for(n = 0; n < params->ndf; n++) {
      fprintf(op,"file%1.1d.name           %s\n",n,params->d_name[n]);
      fprintf(op,"file%1.1d.filetype       %s\n",n,params->d_filetype[n]);
      if (params->d_tstart != NULL && strlen(params->d_tstart[n]))
	fprintf(op,"file%1.1d.tstart         %s\n",n,params->d_tstart[n]);
      fprintf(op,"file%1.1d.tinc           %s\n",n,params->d_tinc[n]);
      if (params->d_tstop != NULL && strlen(params->d_tstop[n]))
	fprintf(op,"file%1.1d.tstop          %s\n",n,params->d_tstop[n]);
      fprintf(op,"file%1.1d.bytespervalue  4\n",n);
      if (params->d_vars != NULL)
	fprintf(op,"file%1.1d.vars           %s\n",n,params->d_vars[n]);
      if (params->d_tunit != NULL && strlen(params->d_tunit[n]))
	fprintf(op,"file%1.1d.timeunit       %s\n",n,params->d_tunit[n]);
      fprintf(op,"\n");
    }
  } else {
    if (params->us_type & US_IJ)
      strcpy(tag, "standard");
    else
      strcpy(tag, "ugrid");
    fprintf(op, "# Output files\n");
    fprintf(op, "OutputFiles          1\n\n");
    fprintf(op, "file0.name           out1\n");
    fprintf(op, "file0.filetype       %s\n", tag);
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
    if (params->win_type & GROUPED)
      fprintf(op, "DP_MODE              %s\n", params->dp_mode);
    if (params->win_type & GROUPED)
      fprintf(op, "WINDOW_TYPE          GROUPED\n");
    else if (params->win_type & STRIPE_E1)
      fprintf(op, "WINDOW_TYPE          STRIPE_E1\n");
    else if (params->win_type & STRIPE_E2)
      fprintf(op, "WINDOW_TYPE          STRIPE_E2\n");
    else if (params->win_type & BLOCK_E1)
      fprintf(op, "WINDOW_TYPE          BLOCK_E1\n");
    else if (params->win_type & BLOCK_E2)
      fprintf(op, "WINDOW_TYPE          BLOCK_E2\n");
    else if (params->win_type & WIN_METIS)
      fprintf(op, "WINDOW_TYPE          METIS\n");
    if (params->show_win)
      fprintf(op, "SHOW_WINDOWS         YES\n");
    if (params->map_type & (WIN_DUMP|GEOM_DUMP)) {
      fprintf(op, "DUMP_MAP            ");
      if (params->map_type & WIN_DUMP && strlen(params->wind_file))
	fprintf(op, " WIN:%s", params->wind_file);
      if (params->map_type & GEOM_DUMP && strlen(params->geom_file))
	fprintf(op, " GEOM:%s", params->geom_file);
      fprintf(op,"\n");
    }
    if (params->map_type & (WIN_READ|GEOM_READ|WIN_CHECK|GEOM_CHECK)) {
      fprintf(op, "READ_MAP            ");
      if (params->map_type & WIN_READ && strlen(params->win_file))
	fprintf(op, " WIN:%s", params->win_file);
      if (params->map_type & GEOM_READ && strlen(params->geom_file))
	fprintf(op, " GEOM:%s", params->geom_file);
      if (params->map_type & GEOM_CHECK)
	fprintf(op, " CHECK_G:%s", params->geom_file);
      if (params->map_type & WIN_CHECK)
	fprintf(op, " CHECK_W:%s", params->win_file);
      fprintf(op,"\n");
    }
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
  fprintf(op, "TRATIO               %-8.1f\n\n", params->tratio);

  fprintf(op, "# Advection\n");
  fprintf(op, "MOM_SCHEME           %s\n", adname(params->momsc));
  fprintf(op, "TRA_SCHEME           %s\n", adname(params->trasc));
  fprintf(op, "ULTIMATE             %s\n\n", tf(params->ultimate));

  fprintf(op, "# Horizontal mixing\n");
  if (params->u1vh > 1.0)
    fprintf(op, "U1VH                 %-6.1f\n", params->u1vh);
  else
    fprintf(op, "U1VH                 %-8.3f\n", params->u1vh);
  if (params->u1kh > 1.0)
    fprintf(op, "U1KH                 %-6.1f\n", params->u1kh);
  else
    fprintf(op, "U1KH                 %-8.3f\n", params->u1kh);
  /*fprintf(op, "SMAGORINSKY          %-6.4f\n", params->smagorinsky);*/
  if (strlen(params->smag))
    fprintf(op, "SMAGORINSKY          %s\n", params->smag);
  if (params->smag_smooth)
    fprintf(op, "SMAG_SMOOTH          %d\n", params->smag_smooth);
  if (params->diff_scale == NONE)
    fprintf(op, "DIFF_SCALE           NONE\n");
  if (params->diff_scale & LINEAR)
    fprintf(op, "DIFF_SCALE           LINEAR ");
  if (params->diff_scale & NONLIN)
    fprintf(op, "DIFF_SCALE           NONLIN ");
  if (params->diff_scale & CUBIC)
    fprintf(op, "DIFF_SCALE           CUBIC ");
  if (params->diff_scale & AREAL)
    fprintf(op, "DIFF_SCALE           AREA ");
  if (params->diff_scale & AUTO)
    fprintf(op, "DIFF_SCALE           AUTO ");
  if (params->diff_scale & SCALEBI)
    fprintf(op, "SCALEBI ");
  if (params->diff_scale & SCALE2D)
    fprintf(op, "SCALE2D ");
  if (params->diff_scale & H_LEN)
    fprintf(op, " LENGTH ");
  if (params->diff_scale & E_LEN)
    fprintf(op, " E_AREA ");
  if (params->diff_scale & C_LEN)
    fprintf(op, " C_AREA ");
  fprintf(op, "\n");

  if (params->visc_method == NONE)
    fprintf(op, "VISC_METHOD          NONE\n");
  else if (params->visc_method == LAPLACIAN)
    fprintf(op, "VISC_METHOD          LAPLACIAN\n");
  else if (params->visc_method == US_LAPLACIAN)
    fprintf(op, "VISC_METHOD          US_LAPLACIAN\n");
  else if (params->visc_method == US_BIHARMONIC) {
    if (params->visc_fact == 0.0)
      fprintf(op, "VISC_METHOD          US_BIHARMONIC\n");
    else
      fprintf(op, "VISC_METHOD          US_LAPLACIAN US_BIHARMONIC %2.1f\n", params->visc_fact);
  }
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
    if (params->min_tke)
      fprintf(op, "MIN_TKE              %-6.4e\n\n", params->min_tke);
    fprintf(op, "ZS                   %-6.4f\n\n", params->zs);
  }

  fprintf(op, "# Bottom friction\n");
  fprintf(op, "QBFC                 %s\n", params->quad_bfc);
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
    if (params->eta_f && params->etarlx & ALERT &&
	!(params->etarlx & RELAX)) {
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
  if (params->compatible & V4201)
    fprintf(op, "COMPATIBLE           V4201\n");
  fprintf(op, "\n");
  
  write_grid_specs(op, params);
  fprintf(op, "\n# Mean horizontal edge length = %8.2f m\n", master->hmean1);
  fprintf(op, "# Mean horizontal distance between centres = %8.2f m\n", master->hmean2);
  fprintf(op, "# Minimum horizontal distance between centres = %8.2f m [%5.3f %5.3f]\n", 
	  master->minres, geom->u1x[master->minrese], geom->u1y[master->minrese]);
  fprintf(op, "# Maximum horizontal distance between centres = %8.2f m [%5.3f %5.3f]\n", 
	  master->maxres, geom->u1x[master->maxrese], geom->u1y[master->maxrese]);
  fprintf(op, "# Mean cell area = %5.4e m^2\n\n", master->amean);

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
  if (strlen(params->maxdiff))
    fprintf(op, "MAXDIFF              %s\n", params->maxdiff);

  fprintf(op, "\n# Tracers\n");
  ntr = (params->roammode & (A_RECOM_R1|A_RECOM_R2)) ? 2 : params->ntr;
  atr = (params->roammode & (A_RECOM_R1|A_RECOM_R2)) ? 0 : params->atr;
  fprintf(op, "NTRACERS             %d\n\n", ntr - atr);
  n = 0;
  for (tn = atr; tn < ntr; tn++) {
    tracer_info_t *tracer = &params->trinfo_3d[tn];
    fprintf(op, "TRACER%1.1d.name         %s\n", n, tracer->name);
    fprintf(op, "TRACER%1.1d.long_name    %s\n", n, tracer->long_name);
    fprintf(op, "TRACER%1.1d.units        %s\n", n, tracer->units);
    if (tracer->fill_value_wc < 1)
      fprintf(op, "TRACER%1.1d.fill_value   %-4.2e\n", n,
              tracer->fill_value_wc);
    else
      fprintf(op, "TRACER%1.1d.fill_value   %-6.1f\n", n,
              tracer->fill_value_wc);
    if (tn > 1)
      fprintf(op, "TRACER%1.1d.valid_range  %-4.2e %-4.2e\n", n,
	      tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
    else
      fprintf(op, "TRACER%1.1d.valid_range  %-6.1f %-6.1f\n", n,
	      tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
    fprintf(op, "TRACER%1.1d.advect       %d\n", n, tracer->advect);
    fprintf(op, "TRACER%1.1d.diffuse      %d\n", n, tracer->diffuse);
    fprintf(op, "TRACER%1.1d.diagn        %d\n", n, tracer->diagn);
    if(!(params->roammode & (A_RECOM_R1|A_RECOM_R2)) && (tracer->type & (WATER|SEDIM) || tracer->type & INTER))
      fprintf(op, "TRACER%1.1d.type         %s\n", n, trtypename(tracer->type, key));
    /*
    fprintf(op, "TRACER%1.1d.decay        %s\n", n, tracer->decay);
    fprintf(op, "TRACER%1.1d.svel         %-6.1f\n", n, tracer->svel);
    fprintf(op, "TRACER%1.1d.resusp_rate  %-6.1f\n", n,
            tracer->resusp_rate);
    fprintf(op, "TRACER%1.1d.crit_stress  %-6.1f\n", n,
            tracer->crit_stress);
    */
    if (strlen(tracer->tracerstat))
      fprintf(op, "TRACER%1.1d.tracerstat   %s\n", n, tracer->tracerstat);
    if (strlen(params->trinfo_3d[tn].data))
      fprintf(op, "TRACER%1.1d.data         %s\n", n, params->trinfo_3d[tn].data);
    if (params->trrlxn != NULL && strlen(params->trrlxn[tn])) {
      fprintf(op, "TRACER%1.1d.relaxation_file     %s\n",
              n, params->trrlxn[tn]);
      fprintf(op, "TRACER%1.1d.relaxation_input_dt %s\n",
              n, otime(params->trrlxdt[tn], tag));
      fprintf(op, "TRACER%1.1d.relaxation_time_constant  %s\n",
              n, params->trrlxtc[tn]);
    }
    if (params->trrest != NULL && strlen(params->trrest[tn])) {
      fprintf(op, "TRACER%1.1d.reset_file   %s\n",
              n, params->trrest[tn]);
      if (params->save_force && (strcmp(tracer->name, "otemp") == 0 ||
				 strcmp(tracer->name, "osalt") == 0 ||
				 strcmp(tracer->name, "ovelu") == 0 ||
				 strcmp(tracer->name, "ovelv") == 0))
	fprintf(op, "TRACER%1.1d.reset_dt     1 day\n", n);
	/*
	fprintf(op, "TRACER%1.1d.reset_dt     %-8.2f seconds\n",
		n, params->grid_dt);
	*/
      else
	fprintf(op, "TRACER%1.1d.reset_dt     %s\n",
		n, otime(params->trrestdt[tn], tag));
    }
    n++;
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
    fprintf(op, "eta_relaxation_time_constant  %s\n\n", params->etarlxtcs);
    /*	    otime(params->etarlxtc, tag));*/
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
    fprintf(op, "SWR_TRANSMISION     %s\n", params->swr_tran);
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
    if (params->us_type & US_IJ) {
      fprintf(op, "TS0.location         %-10.4f %-10.4f 0\n",
	      dumpdata->cellx[params->nce2 / 2][params->nce1 / 2],
	      dumpdata->celly[params->nce2 / 2][params->nce1 / 2]);
    } else {
      master_t *master = dumpdata->master;
      geometry_t *geom = master->geom;
      int cc = geom->b2_t / 2;
      int c = geom->w2_t[cc];
      fprintf(op, "TS0.location         %-10.4f %-10.4f 0\n",
	      geom->cellx[c], geom->celly[c]);
    }
    fprintf(op, "TS0.dt               1 hour\n");
    fprintf(op, "TS0.reference        msl\n\n");
  }

  if (strlen(params->cookiecut)) {
    fprintf(op, "NBOUNDARIES           0\n\n");
    fclose(op);
    return;
  }

  fprintf(op, "# Open boundaries\n");
  fprintf(op, "NBOUNDARIES           %d\n\n", params->nobc);
  if (strlen(params->nodal_dir))
    fprintf(op, "TIDE_CSR_CON_DIR %s\n", params->nodal_dir);
  if (strlen(params->orthoweights))
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
      fprintf(op, "BOUNDARY%1.1d.NSPONGE_HORZ  %5.1f\n", n,
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
    if (open->adjust_flux_s)
      fprintf(op, "BOUNDARY%1.1d.ADJUST_TIDE   %s\n", n,
	      otime(open->adjust_flux_s, tag));      
    if (open->relax_zone_nor)
      fprintf(op, "BOUNDARY%1.1d.RELAX_ZONE_NOR %d\n", n, open->relax_zone_nor);
    if (open->relax_zone_tan)
      fprintf(op, "BOUNDARY%1.1d.RELAX_ZONE_TAN %d\n", n, open->relax_zone_tan);
    if (strlen(open->scale_e))
      fprintf(op, "BOUNDARY%1.1d.SCALE_ETA      %s\n", n, open->scale_e);
    if (strlen(open->tide_con))
      fprintf(op, "BOUNDARY%1.1d.T_CONSTITUENTS  %s\n", n, open->tide_con);
    if (open->stagger & INFACE)
      fprintf(op, "BOUNDARY%1.1d.STAGGER       INFACE\n", n);
    if (open->bathycon)
      fprintf(op, "BOUNDARY%1.1d.BATHY_CON     %d\n", n, open->bathycon);
    if(params->us_type & US_IJ) {
      fprintf(op, "BOUNDARY%1.1d.POINTS        %d\n", n, open->npts);
      for (tn = 0; tn < open->npts; tn++)
	fprintf(op, "%d %d\n", open->iloc[tn], open->jloc[tn]);
    } else {
      mesh_t *m = params->mesh;
      fprintf(op, "BOUNDARY%1.1d.UPOINTS        %d\n", n, m->npts[n]);
      for (i = 1; i <= m->npts[n]; i++) {
	/*	
	fprintf(op, "%d (%f,%f)-(%f,%f)\n", m->loc[n][i],
		m->xloc[m->obc[n][i][0]], m->yloc[m->obc[n][i][0]],
		m->xloc[m->obc[n][i][1]], m->yloc[m->obc[n][i][1]]);
	*/
	fprintf(op, "%d (%d %d)\n", m->loc[n][i], m->obc[n][i][0], m->obc[n][i][1]);
      }
    }
    fprintf(op, "\n");
  }

  if (params->us_type & US_IJ) {
    fprintf(op, "# Bathymetry\n");
    fprintf(op, "BATHY    %d\n", params->nce1 * params->nce2);
    for (j = 0; j < params->nce2; j++)
      for (i = 0; i < params->nce1; i++)
	if (isnan(dumpdata->botz[j][i])) {
	  fprintf(op, "%8.3f\n", -(double)LANDCELL);
	}
	else {
	  if (dumpdata->botz[j][i] > -params->bmin) dumpdata->botz[j][i] = (double)LANDCELL;
	  fprintf(op, "%8.3f\n", -dumpdata->botz[j][i]);
	}
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
  } else {
    write_mesh_us(params, params->bathy, op, 2);
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


/*-------------------------------------------------------------------*/
/* Writes a template for an auto file                                */
/*-------------------------------------------------------------------*/
void auto_write(char *name)
{
  FILE *fp;
  char buf[MAXSTRLEN];

  sprintf(buf, "%s.prm", name);
  fp = fopen(buf, "w");

  fprintf(fp,"###################################\n");
  fprintf(fp,"# auto file for %s\n\n", name);
  fprintf(fp,"TIMEUNIT             seconds since 1990-01-01 00:00:00\n");
  fprintf(fp,"OUTPUT_TIMEUNIT      days since 1990-01-01 00:00:00\n");
  fprintf(fp,"LENUNIT              metre\n");
  fprintf(fp,"START_TIME           0 days\n");
  fprintf(fp,"STOP_TIME            1 days\n");
  fprintf(fp,"NAME                 %s\n", name);
  fprintf(fp,"ID_NUMBER            1.0\n");
  fprintf(fp,"PARAMETERHEADER      Auto %s\n\n", name);
  fprintf(fp,"# Input file\n");
  fprintf(fp,"INPUT_FILE           in\n\n");
  fprintf(fp,"# Grid\n");
  fprintf(fp,"projection           geographic\n");
  fprintf(fp,"GRIDTYPE             UNSTRUCTURED\n");
  fprintf(fp,"INTERIOR_COORD       0.0 0.0\n");
  fprintf(fp,"NUMBERS              RESOLUTION\n");
  fprintf(fp,"\n###################################\n");
  fprintf(fp,"# JIGSAW options\n");
  fprintf(fp,"JIG_GEOM_FILE        %s_geom.msh\n", name);
  fprintf(fp,"HFUN_FILE            %s_hfun.msh\n", name);
  fprintf(fp,"POWER_MESH           YES\n");
  fprintf(fp,"STEREOGRAPHIC_MESH   NO\n");
  fprintf(fp,"\n###################################\n");
  fprintf(fp,"# Weighting function\n");
  fprintf(fp,"# HMIN/HMAX = resolution\n");
  fprintf(fp,"# BMIN/BMAX = bathymetry / distance from coast\n");
  fprintf(fp,"HFUN_BATHY_FILE      COAST\n");
  fprintf(fp,"HFUN_GRID            yes\n");
  fprintf(fp,"HFUN_SMOOTH          1\n\n");
  fprintf(fp,"NHFUN                2\n");
  fprintf(fp,"HMIN0                1000.0\n");
  fprintf(fp,"HMAX0                1500.0\n");
  fprintf(fp,"BMIN0                200.0\n");
  fprintf(fp,"BMAX0                25000.0\n");
  fprintf(fp,"TYPE0                0.0\n\n");
  fprintf(fp,"HMIN1                1500.0\n");
  fprintf(fp,"HMAX1                3000.0\n");
  fprintf(fp,"BMIN1                25000.0\n");
  fprintf(fp,"BMAX1                150000.0\n");
  fprintf(fp,"TYPE1                0.0\n\n");
  fprintf(fp,"#HFUN_EXCLUDE         0\n");
  fprintf(fp,"#HFUN_OVERRIDE        0\n");
  fprintf(fp,"\n###################################\n");
  fprintf(fp,"# Coastmesh\n");
  fprintf(fp,"COASTFILE            coastfile.cst\n");
  fprintf(fp,"PLOTFILE             %s_p\n", name);
  fprintf(fp,"CUTOFF               50\n");
  fprintf(fp,"#RADIUS               1000\n");
  fprintf(fp,"#LENGTH               1000\n");
  fprintf(fp,"MINLON               0.0\n");
  fprintf(fp,"MINLAT               0.0\n");
  fprintf(fp,"MAXLON               0.0\n");
  fprintf(fp,"MAXLAT               0.0\n");
  fprintf(fp,"START_COORD          0.0 0.0\n");
  fprintf(fp,"MID_COORD            0.0 0.0\n");
  fprintf(fp,"END_COORD            0.0 0.0\n");
  fprintf(fp,"SMOOTH               1\n");
  fprintf(fp,"SMOOTHZ              2\n");
  fprintf(fp,"#RESAMPLE             2\n");
  fprintf(fp,"NOBC                 1\n");
  fprintf(fp,"OBC0                 obc.xy\n\n");
  fprintf(fp,"#LINKS		    0\n");
  fprintf(fp,"\n###################################\n");
  fprintf(fp,"# Bathymetry\n");
  fprintf(fp,"BATHYFILE            bathy.nc\n");
  fprintf(fp,"#BATHY_INTERP_RULE    nn_non_sibson\n");
  fprintf(fp,"BATHY_TRUNCATE       NO\n");
  fprintf(fp,"BATHY_REPORT         NO\n");
  fprintf(fp,"BATHYFILL            AVERAGE\n");
  fprintf(fp,"BATHYMIN             1.0\n");
  fprintf(fp,"BATHYMAX             6000.0\n");
  fprintf(fp,"MIN_CELL_THICKNESS   20%\n");
  fprintf(fp,"SMOOTHING            1\n");
  fprintf(fp,"MAXDIFF              median\n");
  fprintf(fp,"MAXGRAD              0.04\n");
  fprintf(fp,"\n###################################\n");
  fprintf(fp,"# Boundaries and forcing\n");
  fprintf(fp,"NBOUNDARIES          0\n\n");
  fprintf(fp,"outputfiles          0\n\n");
  fclose(fp);
}

/* END auto_write()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Common to both param/trans_write                                  */
/*-------------------------------------------------------------------*/
void write_grid_specs(FILE *op, parameters_t *params)
{
  mesh_t *mesh = params->mesh;
  fprintf(op, "# Grid\n");
  if (strlen(params->roms_grid))
    fprintf(op, "ROMS_GRID            %s\n", params->roms_grid);
  else
    fprintf(op, "GRIDTYPE             %s\n", params->gridtype);

  if (params->us_type & (US_RS|US_WS|US_G))
    fprintf(op, "INPUT_FILE_TYPE      STRUCTURED\n");
  if (params->us_type & (US_RUS|US_WUS))
    fprintf(op, "INPUT_FILE_TYPE      UNSTRUCTURED\n");

  if (params->us_type & US_POW)
    fprintf(op, "POWER_MESH           YES\n");
  if (params->us_type & US_STER)
    fprintf(op, "STEREOGRAPHIC_MESH   YES\n");

  if (params->us_type & US_IJ) {
    fprintf(op, "NCE1                 %d\n", params->nce1);
    fprintf(op, "NCE2                 %d\n", params->nce2);
  } else {
    fprintf(op, "Mesh2 unstructured   v1.0\n");
    fprintf(op, "nMaxMesh2_face_nodes %d\n", mesh->mnpe);
    fprintf(op, "nMesh2_face_indices  %d\n", mesh->ns);
    fprintf(op, "nMesh2_face          %d\n",mesh->ns2);
  }

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

/* END write_grid_specs()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculate cell and face neighbours for non-periodic structured    */
/* grids.                                                            */
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
  for (j = 0; j < dumpdata->nce2; j++)
    for (i = 0; i < dumpdata->nce1; i++) {
      dumpdata->cellx[j][i] = NaN;
      dumpdata->celly[j][i] = NaN;
    }
  for (j = 0; j < dumpdata->nce2; j++)
    for (i = 0; i < dumpdata->nfe1; i++) {
      dumpdata->u1x[j][i] = NaN;
      dumpdata->u1y[j][i] = NaN;
    }
  for (j = 0; j < dumpdata->nfe2; j++)
    for (i = 0; i < dumpdata->nce1; i++) {
      dumpdata->u2x[j][i] = NaN;
      dumpdata->u2y[j][i] = NaN;
    }
  for (j = 0; j < dumpdata->nfe2; j++)
    for (i = 0; i < dumpdata->nfe1; i++) {
      dumpdata->gridx[j][i] = NaN;
      dumpdata->gridy[j][i] = NaN;
    }
}

/* END neighbour_none()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the grid specification parameters                 */
/*-------------------------------------------------------------------*/
void read_grid(parameters_t *params)
{
  FILE *fp;                     /* Input parameter file              */
  char buf[MAXSTRLEN];          /* Buffer                            */
  double slat, slon;            /* North-western domain coordinates  */
  double elat, elon;            /* South-eastern domain coordinates  */
  double x00;                   /* Grid e1 offset                    */
  double y00;                   /* Grid e2 offset                    */
  double rotn;                  /* Grid rotation                     */
  double xinc;                  /* e1 grid spacing                   */
  double yinc;                  /* e2 grid spacing                   */
  double lx, ly;                /* Size of the grid                  */
  double arc, rmin;             /* Polar grid attributes             */
  double ella;                  /* Elliptic grid attributes          */
  double taumax, taumin;        /* Elliptic grid attributes          */
  double rmax, crad;            /* Circular hex attributes           */
  int nsimm;                    /* Elliptic grid attributes          */
  int nce1;                     /* e1 grid dimension                 */
  int nce2;                     /* e2 grid dimension                 */
  /* double nm=1852.0; */       /* Meters in a nautical mile         */
  /* double dm=60.0; */         /* Minutes in 1 degree               */
  int i, j;

  fp =params->prmfd;
  nce1 = params->nce1;
  nce2 = params->nce2;

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  /* Allocate memory for the coordinate and metric arrays for        */
  /* structured grids which are at twice the final resolution.       */
  if (nce1 && nce2) {
    params->x = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
    params->y = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
    params->h1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
    params->h2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
    params->a1 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
    params->a2 = d_alloc_2d(2 * nce1 + 1, 2 * nce2 + 1);
  }
  xinc = yinc = rotn = params->flat = params->flon = 0.0;
  prm_set_errfn(hd_silent_warn);

  memset(projection, 0, sizeof(projection));
  prm_read_char(fp, "PROJECTION", params->projection);

  prm_read_char(fp, "TESTCASE", params->testcase);

  /*-----------------------------------------------------------------*/
  /* Specified unstructured meshes                                   */
  if (prm_read_int(fp, "nMaxMesh2_face_nodes", &params->npe) &&
      prm_read_int(fp, "nMesh2_face_indices", &params->ns) &&
      prm_read_int(fp, "nMesh2_face", &params->ns2)) {
    i = read_mesh_us(params);
    strcpy(params->gridtype, "UNSTRUCTURED");
    params->gridcode = UNSTRUCTURED;
    params->us_type |= US_G;
    if (i > 0) return;

  } else if (prm_read_double(fp, "HEX_RADIUS", &crad)) {
#ifdef HAVE_JIGSAWLIB

    /*---------------------------------------------------------------*/
    /* Circular unstructured meshes                                  */
    jigsaw_msh_t J_mesh;
    double gscale = 0.25;
    int powf = (params->us_type & US_POW) ? 1 : 0;
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);
    prm_read_double(fp, "DX_MIN", &rmin);
    prm_read_double(fp, "DX_MAX", &rmax);
    prm_read_double(fp, "HEX_GSCALE", &gscale);
    strcpy(params->gridtype, "UNSTRUCTURED");
    params->gridcode = UNSTRUCTURED;
    params->us_type |= US_G;

    /* Initialise JISGSAW output mesh                                */
    jigsaw_init_msh_t(&J_mesh);
    
    /* Creates circular coastline, hfun and runs JIGSAW              */
    create_hex_radius(crad, x00, y00, rmin, rmax, gscale, &J_mesh, powf);

    /* Copy out JIGSAW mesh into COMPAS data structures              */
    convert_jigsaw_msh(params, NULL, &J_mesh);

    /* No longer need this */
    jigsaw_free_msh_t(&J_mesh);

    write_mesh_desc(params, NULL, NULL, (M_CONVERT_JIG|M_MESHSTRUCT));
#else
    hd_quit("read_grid: HEX_RADIUS requires the JIGSAW library\n");
#endif
  } else if (prm_read_char(fp, "JIGSAW_FILE", buf)) {
    /*---------------------------------------------------------------*/
    /* JIGSAW file input unstructured meshes                         */
    convert_jigsaw_msh(params, buf, NULL); 
    write_mesh_desc(params, NULL, NULL, (M_CREATE_JIG|M_CONVERT_JIG|M_MESHSTRUCT));
    strcpy(params->gridtype, "UNSTRUCTURED");
    params->gridcode = UNSTRUCTURED;
    params->us_type |= US_G;
  } else if (prm_read_char(fp, "HFUN_BATHY_FILE", buf) ||
	     prm_read_char(fp, "WEIGHTING", buf)) {
#ifdef HAVE_JIGSAWLIB
    char *fields[MAXSTRLEN * MAXNUMARGS];
    jigsaw_msh_t J_mesh;
    jigsaw_msh_t J_hfun;
    coamsh_t *cm;
    int powf = (params->us_type & US_POW) ? 1 : 0;
    int stproj = (params->us_type & US_STER) ? 1 : 0;
    int mode;

    mode = 0;
    if(endswith(buf,".nc"))
      mode |= (H_BATHY|H_NC);
    if(endswith(buf,".bty"))
      mode |= (H_BATHY|H_BTY);
    if(endswith(buf,".msh"))
      mode |= (H_BATHY|H_MSH);
    if (contains_token(buf, "GWS") != NULL) {
      mode |= H_GWS;
      i = parseline(buf, fields, MAXNUMTSFILES);
      strcpy(buf, fields[i-1]);
    }
    if (contains_token(buf, "GRAD") != NULL) {
      mode |= H_GRAD;
      i = parseline(buf, fields, MAXNUMTSFILES);
      strcpy(buf, fields[i-1]);
    }
    if (contains_token(buf, "COAST") != NULL)
      mode |= (H_COAST|H_CST);
    if (contains_token(buf, "POINT") != NULL)
      mode |= (H_COAST|H_POINT);
    if (contains_token(buf, "POLY") != NULL)
      mode |= (H_COAST|H_POLY);
    if (!(mode & (H_BATHY|H_COAST)))
      mode |= (H_COAST|H_CONST);

    /* Set up the bounding perimeters                                */ 
    cm = cm_alloc();
    cm_read(cm, params->prmfd);
    if (prm_skip_to_end_of_key(fp, "JIG_GEOM_FILE"))
      coastmesh(cm, (CM_LIB|CM_FILE));
    else
      coastmesh(cm, CM_LIB);

    /* Build the hfun file                                           */
    jigsaw_init_msh_t (&J_hfun);
    if (mode & H_COAST)
      hfun_from_coast(params, cm, &J_hfun, mode);
    if (mode & H_BATHY)
      hfun_from_bathy(params, buf, cm, &J_hfun, mode);
    /*
    if (strcmp(buf, "COAST") == 0)
      hfun_from_coast(params, cm, &J_hfun, 0);
    else if (strcmp(buf, "POINT") == 0)
      hfun_from_coast(params, cm, &J_hfun, 1);
    else if (strcmp(buf, "POLY") == 0)
      hfun_from_coast(params, cm, &J_hfun, 2);
    else
      hfun_from_bathy(params, buf, cm, &J_hfun);
    */

    /* Initialise JISGSAW output mesh                                */
    jigsaw_init_msh_t(&J_mesh);

    /* Create the mesh                                               */
    create_jigsaw_mesh(cm, &J_mesh, &J_hfun, powf, stproj, params->meshinfo);

    /* Copy out JIGSAW mesh into COMPAS data structures              */
    convert_jigsaw_msh(params, NULL, &J_mesh);

    write_mesh_desc(params, cm, NULL, (M_BOUND_HFUN|M_CREATE_JIG|M_CONVERT_JIG|M_MESHSTRUCT));

    /* Clean                                                         */
    cm_free(cm);
    jigsaw_free_msh_t(&J_mesh);
    strcpy(params->gridtype, "UNSTRUCTURED");
    params->gridcode = UNSTRUCTURED;
    params->us_type |= US_G;
#else
    hd_quit("read_grid: HFUN_BATHY_FILE requires the JIGSAW library\n");
#endif
  } else if (prm_read_double(fp, "SLAT", &slat) &&
      prm_read_double(fp, "ELAT", &elat) &&
      prm_read_double(fp, "SLON", &slon) &&
      prm_read_double(fp, "ELON", &elon)) {
    /*---------------------------------------------------------------*/
    /* Structured defined by a bounding box                          */
    double mlat, mlon;          /* Lat and long of the mid-point of the
                                   grid */
    strcpy(params->gridtype, "GEOGRAPHIC_RECTANGULAR");
    strcpy(params->projection, "geographic"); /* Override projection */
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
    /*---------------------------------------------------------------*/
    /* Structured ROMS grid conversion                               */
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
    /*---------------------------------------------------------------*/
    /* Structured defined by a corner location and (m) resolution    */
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
    /*---------------------------------------------------------------*/
    /* Structured defined by a corner location and (deg) resolution  */
    strcpy(params->gridtype, "GEOGRAPHIC_RECTANGULAR");
    strcpy(params->projection, "geographic"); /* Override projection */
    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);

    if (prm_read_double(fp, "POLE_LATITUDE", &params->flat) &&
        prm_read_double(fp, "POLE_LONGITUDE", &params->flon)) {
      double elf = 0.0;

      /* Set the longitude increment relative to the grid centre rather than grid origin. */
      /* Note: ROAM compensates for this. */
      /*
      if(!(params->runmode & ROAM))
	xinc *= cos(PI*0.5*(2 * y00 + nce2 * yinc)/180.0);
      */

      /* Hex ellipse option */
      if (prm_read_double(fp, "ELLIPSE", &elf)) {
	jigsaw_msh_t J_mesh;
	jigsaw_msh_t J_hfun;
	double deg2m = 60.0 * 1852.0;
	double ores;       /* Outer ellipse resolution                 */
	strcpy(params->gridtype, "UNSTRUCTURED");
	params->gridcode = UNSTRUCTURED;
	params->us_type |= US_G;

	/* Initialise JISGSAW output mesh                              */
	jigsaw_init_msh_t(&J_mesh);
	jigsaw_init_msh_t (&J_hfun);

	if (!(prm_read_double(fp, "OUTER_RES", &ores))) 
	  ores = 3000.0 / deg2m;

	create_ellipse(params, nce1, nce2,
		       x00, y00, params->flon, params->flat,
		       xinc, yinc, elf, ores, &J_mesh, &J_hfun);

	/* Copy out JIGSAW mesh into COMPAS data structures          */
	convert_jigsaw_msh(params, NULL, &J_mesh);

	write_mesh_desc(params, NULL, NULL, (M_CONVERT_JIG|M_MESHSTRUCT));

	/* No longer need this */
	jigsaw_free_msh_t(&J_mesh);
	jigsaw_free_msh_t(&J_hfun);
      } else {

	geog_false_pole_coord(params->x, params->y, params->h1, params->h2,
			      params->a1, params->a2,
			      2 * nce1, 2 * nce2,
			      x00, y00, params->flon, params->flat,
			      xinc / 2.0, yinc / 2.0);
	
	params->gridcode = GRECT_DLATLON_FP;
      }

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
    /*---------------------------------------------------------------*/
    /* Structured polar grid                                         */
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
    /*---------------------------------------------------------------*/
    /* Structured elliptic grid                                      */
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
    /*---------------------------------------------------------------*/
    /* Specified structured grid                                     */
    read2darray(fp, "XCOORDS", params->x, 2 * nce1 + 1, 2 * nce2 + 1);
    read2darray(fp, "YCOORDS", params->y, 2 * nce1 + 1, 2 * nce2 + 1);
    strcpy(params->gridtype, "NUMERICAL");
    grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
                     params->h1, params->h2);
    grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
                   params->a1, params->a2);
    params->gridcode = NUMERICAL;

  } else if (prm_read_char(fp, "COORDFILE", buf)) {
    /*---------------------------------------------------------------*/
    /* Specified structured grid from file                           */
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
    /*---------------------------------------------------------------*/
    /* Structured MOM grid conversion                                */
    read_mom_grid(buf, nce1, nce2, params->x, params->y);
    strcpy(params->gridtype, "NUMERICAL");
    grid_get_metrics(params->x, params->y, 2 * nce1, 2 * nce2,
		     params->h1, params->h2);
    grid_get_angle(params->x, params->y, 2 * nce1, 2 * nce2,
		   params->a1, params->a2);
    params->gridcode = NUMERICAL;

  } else if (prm_read_double(fp, "TRI_DX", &xinc) ||
	     prm_read_double(fp, "HEX_DX", &xinc)) {
    /*---------------------------------------------------------------*/
    /* Untructured defined by a corner location and resolution       */
    int np;
    double *x, *y;
    point *pin;
    char code[64];

    if (prm_skip_to_end_of_key(fp, "TRI_DX"))
      params->us_type |= US_TRI;
    if (prm_skip_to_end_of_key(fp, "HEX_DX"))
      params->us_type |= US_HEX;

    prm_read_double(fp, "X00", &x00);
    prm_read_double(fp, "Y00", &y00);
    prm_read_double(fp, "ROTATION", &rotn);
    strcpy(params->gridtype, "UNSTRUCTURED");
    if (strlen(params->projection) > 0 && strcmp(params->projection,
"geographic") == 0) {
      strcpy(code, "eznCD");
      x = d_alloc_1d((nce1+1)*(nce2+1) - nce2/2);
      y = d_alloc_1d((nce1+1)*(nce2+1) - nce2/2);
      np = delaunay_dreckon_coord(x, y, nce1, nce2, x00, y00, rotn, xinc);
      params->gridcode = UNSTRUCTURED|GTRI_DXY_ROT;
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
	  params->gridcode = UNSTRUCTURED|GTRI_DLATLON_FP;
	}
      }
    } else {
      strcpy(code, "eznCDq");
      x = d_alloc_1d((nce1+1)*(nce2+1) - nce2/2);
      y = d_alloc_1d((nce1+1)*(nce2+1) - nce2/2);
      np = delaunay_rect_coord(x, y, nce1, nce2, x00, y00, rotn, xinc);
      params->gridcode = UNSTRUCTURED|TRI_DXY_ROT;
    }

    pin = malloc(np * sizeof(point));
    for (i = 0; i < np; i++) {
      pin[i].x = x[i] - x00;
      pin[i].y = y[i] - y00;
    }
    params->d = create_tri_mesh(np, pin, 0, NULL, 0, NULL, code);

    /* Note: sometimes using regadj=1 in convert_hex_mesh() with     */
    /* regular meshes the mesh is truncated. Try setting regadj=0 in */
    /* this case.                                                    */
    if (params->us_type & US_HEX)
      convert_hex_mesh(params, params->d, 0);
    if (params->us_type & US_TRI)
      convert_tri_mesh(params, params->d);

    for (i = 1; i <= params->ns2; i++) {
      params->x[i][0] += x00;
      params->y[i][0] += y00;
      for (j = 1; j <= params->npe2[i]; j++) {
	params->x[i][j] += x00;
	params->y[i][j] += y00;
      }
    }

    /*free((delaunay *)params->d);*/

  } else if (read_mesh_us(params)) {
    /*---------------------------------------------------------------*/
    /* Untructured. For safety: should be captured at the start of   */
    /* the routine.                                                  */
    strcpy(params->gridtype, "UNSTRUCTURED");
    params->gridcode = UNSTRUCTURED;
    params->us_type |= US_G;
  } else
    hd_quit("Cannot read grid information\n");

  /*---------------------------------------------------------------*/
  /* If a structured grid definition is to be converted to         */
  /* an unstructured grid then set the type to a structured read   */
  /* flag, US_RS.                                                  */
  if (!(params->gridcode & UNSTRUCTURED)) params->us_type |= (US_RS|US_IJ);

  /* Set the time-series file projection */
  if (strlen(params->projection) > 0) {
    strcpy(projection, params->projection);
    ts_set_default_proj_type(projection);
  }

  /* If the the projection is specified and is geographic, then */
  /* compute the metrics and angles for the sphere.  */
  if (params->us_type & US_IJ && 
      strncasecmp(params->projection,
		  GEOGRAPHIC_TAG, strlen(GEOGRAPHIC_TAG)) == 0) {
    grid_get_geog_metrics(params->x, params->y, 2 * nce1,
                          2 * nce2, params->h1, params->h2);
    grid_get_geog_angle(params->x, params->y, 2 * nce1,
                        2 * nce2, params->a1, params->a2);

  }
  if (params->osl & L_BILIN & !(params->us_type & US_IJ))
    hd_quit("read_grid: Can only use BILINEAR with quad meshes.\n");
}
/* END read_grid()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Counts the number of OBCs present in a structured grid            */
/*-------------------------------------------------------------------*/
int count_auto_OBC(parameters_t *params, 
		    double **bathy,
		    int **nibi,
		    int ***maski, 
		    int *onr,
		    int *ors,
		    int **ibi, 
		    int **jbi, 
		    int **rbi, 
		    int **rfi,
		    double **rhci
		    )
{
  FILE *fp = params->prmfd;
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];
  int i, j, ii, jj;
  int m;
  int is, ie, js, je;           /* Limits of grid */
  int bdry = 0;
  int nr, rs;
  double d1;
  int *nib, *ib, *jb, *rb, *rf;
  int **mask;
  double *rhc;

  params->nobc = 0;

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
      if (lat > 90.0 || lat < -90.0) { /* Swap */
	d1 = lon;
	lon = lat;
	lat = d1;
      }
      if (find_closest_b(params, bathy, lat, lon, &ii, &jj, 0)) {
	hd_warn("Can't locate river %s in the grid.\n", name);
	continue;
      } else
	hd_warn("River %s dry cell located at cell (%d %d).\n", name, ii, jj);
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
  *onr = nr;
  *ors = rs;
  *nibi = nib;
  *maski = mask;
  *ibi = ib;
  *jbi = jb;
  *rbi = rb;
  *rfi = rf;
  *rhci = rhc;
  return(bdry);
}

/* END count_auto_OBC()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Counts the number of OBCs present in an unstructured grid         */
/*-------------------------------------------------------------------*/
int count_auto_OBC_us(parameters_t *params,
		      int **nibi,
		      int ***jci,
		      int **maski
		      )
{
  int **neic, **neij;
  int cc, cs, cn, i, j, jj;
  int ns2, npe, npem;
  int bdry;
  double *bathy = params->bathy;
  int *nib;
  int **jc;
  int *mask;
  int verbose = 0;
  int bathycon = 1;  /* Set a no-gradient for bathymetry on the boundary */

  remove_OBC_corner(params);

  params->nobc = 0;
  ns2 = params->ns2;
  npem = params->npe;

  /* Get the neighbour mappings                                      */
  neic = i_alloc_2d(ns2+1, npem+1);
  neij = i_alloc_2d(ns2+1, npem+1);
  jc = i_alloc_2d(ns2+1, npem+1);
  mask = i_alloc_1d(ns2+1);
  nib = i_alloc_1d(MAXBDRY);
  memset(nib, 0, MAXBDRY * sizeof(int));

  for (cc = 1; cc <= ns2; cc++) {
    mask[cc] = -1;
    npe = params->npe2[cc];
    for (j = 1; j <= npe; j++) {
      jj = j;
      cn = find_neighbour(cc, params->x, params->y, params->npe2, ns2, &jj);
      neic[j][cc] = cn;
      neij[j][cc] = jj;
    }
  }

  /*-----------------------------------------------------------------*/
  /* OBCs on the limits of the grid (no neighbour). This follows the */
  /* edge of the domain and increments the number of OBCs if there   */
  /* is a discontinuity (e.g. land cell) interrupting the OBC.       */
  for (cc = 1; cc <= ns2; cc++) {
    if (mask[cc] >= 0) continue;
    npe = params->npe2[cc];
    for (j = 1; j <= npe; j++) {
      jc[j][cc] = 0;
      cn = neic[j][cc];
      if (is_wet(params, bathy[cc]) && !cn) {
	/* found=1 if a neighbour of cc is wet, and itself does not  */
	/* have a neighbour in one direction (i.e. is on the grid    */
	/* perimeter).                                               */
	int cno, cso, found = 0;
	bdry |= O_BDRY;
	/* Look for a wet neighbour                                  */
	for (jj = 1; jj <= npe; jj++) {
	  cn = neic[jj][cc];
	  if (cn && is_wet(params, bathy[cn])) {
	    /* Check if cn is on the grid edge                       */
	    for (i = 1; i <= params->npe2[cn]; i++) {
	      cs = neic[i][cn];
	      if (!cs) {
		found = 1;
		cno = cn;
		cso = cs;
	      }
	    }
	    if (found) {
	      if (mask[cn] >= 0) {
		mask[cc] = mask[cn];
	      }
	    }
	  }
	}
	if (mask[cc] == -1) {
	  mask[cc] = params->nobc;
	  params->nobc++;
	  if (verbose) printf("nobc=%d cc=%d j=%d found=%d(%d %d) %f\n",
			      params->nobc, cc, j, found, cno, cso, bathy[cc]);
	}
	nib[mask[cc]]++;
	jc[j][cc] = 1;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a no-gradient for bathymetry on the boundary if required    */
  if (bathycon) {
    for (cc = 1; cc <= ns2; cc++) {
      if (mask[cc] >= 0) {
	npe = params->npe2[cc];
	for (j = 1; j <= npe; j++) {
	  if (jc[j][cc]) {
	    jj = jo(j, npe);
	    cn = neic[jj][cc];
	    if (mask[cn] >= 0)
	      cn = neic[jj][cn];
	    if (cn && is_wet(params, bathy[cn])) {
	      params->bathy[cc] = params->bathy[cn];
	      break;
	    }
	  }
	}
      }
    }
  }

  if (params->nobc > MAXBDRY) 
    hd_quit("Number of boundaries found > allowed maximum (%d).\n", MAXBDRY);

  /*-----------------------------------------------------------------*/
  /* Get open boundaries in the interior (dry neighbour)             */
  for (cc = 1; cc <= ns2; cc++) {
    if (mask[cc] >= 0) continue;
    npe = params->npe2[cc];
    for (j = 1; j <= npe; j++) {
      cn = neic[j][cc];
      if (is_wet(params, bathy[cc]) && is_outside(params, bathy[cn])) {
	/* found=1 if a neighbour of cc is wet, and itself has an    */
	/* OUTSIDE neighbour in any direction.                       */
	int cno, cso, found = 0;
	bdry |= O_BDRY;
	/* Look for a wet neighbour                                  */
	for (jj = 1; jj <= npe; jj++) {
	  cn = neic[jj][cc];
	  if (cn && is_wet(params, bathy[cn])) {
	    /* Check if cn is on the grid edge                       */
	    for (i = 1; i <= params->npe2[cn]; i++) {
	      cs = neic[i][cn];
	      if (is_outside(params, bathy[cs])) {
		found = 1;
		cno = cn;
		cso = cs;
	      }
	    }
	    if (found) {
	      if (mask[cn] >= 0) {
		mask[cc] = mask[cn];
	      }
	    }
	  }
	}
	if (mask[cc] == -1) {
	  mask[cc] = params->nobc;
	  params->nobc++;
	  if (verbose) printf("nobc=%d cc=%d j=%d found=%d(%d %d) %f\n",
			      params->nobc, cc, j, found, cno, cso, bathy[cc]);
	}
	nib[mask[cc]]++;
	jc[j][cc] = 1;
      }
    }
  }

  *nibi = nib;
  *maski = mask;
  *jci = jc;
  i_free_2d(neic);
  i_free_2d(neij);
  return(bdry);
}

/* END count_auto_OBC_us()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the OBCs present in a structured grid                        */
/*-------------------------------------------------------------------*/
void set_auto_OBC(parameters_t *params, 
		  double **bathy,
		  int bdry,
		  int *nib,
		  int **mask, 
		  int *onr,
		  int *ors,
		  int *ib, 
		  int *jb, 
		  int *rb, 
		  int *rf,
		  double *rhc
		  )
{
  open_bdrys_t *open;           /* Pointer to boundary structure     */
  FILE *fp = params->prmfd;
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];
  int i, j, ii, jj;
  int n, m;
  int is, ie, js, je;           /* Limits of grid                    */
  int nr, rs;
  double d1;

  /* Allocate OBC location vectors                                   */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    if (nib[n]) {
      open->iloc = i_alloc_1d(abs(nib[n]));
      open->jloc = i_alloc_1d(abs(nib[n]));
    }
  }

  /* Get open boundaries on the grid edges                           */
  n = 0;
  if (bdry & W_BDRY) {
    i = 0;
    open = params->open[n];
    open->type = U1BDRY;
    open->intype = O_POI;
    strcpy(open->name, "West");
    open->intype = O_POI;
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
    open->intype = O_POI;
    strcpy(open->name, "East");
    open->intype = O_POI;
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
    open->intype = O_POI;
    strcpy(open->name, "South");
    open->intype = O_POI;
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
    open->intype = O_POI;
    strcpy(open->name, "North");
    open->intype = O_POI;
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

  /* Get open boundaries in the interior                             */
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
      /* L_EDGE                                                      */
      while (is_wet(params, bathy[jj][i]) && 
	     is_outside(params, bathy[jj][i-1]) && 
	     mask[jj][i-1] & L_EDGE && !(mask[jj][i-1] & DRY) &&
	     jj < params->nce2-1) {
	mask[jj][i-1] |= DRY;
        open->type = U1BDRY;
	open->intype = O_POI;
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
      /* R_EDGE                                                      */
      while (is_wet(params, bathy[jj][i]) &&
	     is_outside(params, bathy[jj][i+1]) &&  
	     mask[jj][i+1] & R_EDGE && !(mask[jj][i+1] & DRY) &&
	     jj < params->nce2-1) {
	mask[jj][i+1] |= DRY;
        open->type = U1BDRY;
	open->intype = O_POI;
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
      /* B_EDGE                                                      */
      while (is_wet(params, bathy[j][ii]) && 
	     is_outside(params, bathy[j-1][ii]) && 
	     mask[j-1][ii] & B_EDGE && !(mask[j-1][ii] & DRY) &&
	     ii < params->nce1-1) {
	mask[j-1][ii] |= DRY;
        open->type = U2BDRY;
	open->intype = O_POI;
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
      /* F_EDGE                                                      */
      while (is_wet(params, bathy[j][ii]) &&
	     is_outside(params, bathy[j+1][ii]) &&  
	     mask[j+1][ii] & F_EDGE && !(mask[j+1][ii] & DRY) &&
	     ii < params->nce1-1) {
	mask[j+1][ii] |= DRY;
        open->type = U2BDRY;
	open->intype = O_POI;
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

  /* Get the RIVER open boundaries                                   */
  for (m = 0; m < nr; m++) {
    open = params->open[n];
    open->npts = 0;
    if (rf[m] & U1BDRY)
      open->type = U1BDRY;
    else
      open->type = U2BDRY;
    open->intype = O_POI;
    open->iloc[open->npts] = ib[m];
    open->jloc[open->npts] = jb[m];
    open->npts++;
    sprintf(open->name, "river%1.1d", m);
  }
}

/* END set_auto_OBC()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the OBCs present in an unstructured grid                     */
/*-------------------------------------------------------------------*/
int set_auto_OBC_us(parameters_t *params,
		    int *nib,
		    int **jc,   /* = 1 for boundary edges            */ 
		    int *mask   /* = n for boundary n                */
		    )
{
  open_bdrys_t *open;               /* Pointer to boundary structure */
  int cc, cs, cn, i, j, jj, n;
  int verbose = 0;

  /* Allocate OBC location vectors                                   */
  for (n = 0; n < params->nobc; n++) {
    open = params->open[n];
    if (nib[n]) {
      open->npts = nib[n];
      open->locu = i_alloc_1d(open->npts);
      open->posx = d_alloc_2d(2 ,open->npts);
      open->posy = d_alloc_2d(2 ,open->npts);
    }
  }

  /* Get open boundaries on the grid edges                           */
  memset(nib, 0, MAXBDRY * sizeof(int));
  for (cc = 1; cc <= params->ns2; cc++) {
    n = mask[cc];
    if (n < 0) continue;
    open = params->open[n];
    sprintf(open->name, "OBC%d", n);
    open->type = U1BDRY;
    open->intype = O_UPC;
    for (j = 1; j <= params->npe2[cc]; j++) {
      if (jc[j][cc]) {
	open->locu[nib[n]] = cc;
	jj = (j == params->npe2[cc]) ? 1 : j + 1;
	open->posx[nib[n]][0] = params->x[cc][j];
	open->posx[nib[n]][1] = params->x[cc][jj];
	open->posy[nib[n]][0] = params->y[cc][j];
	open->posy[nib[n]][1] = params->y[cc][jj];
	nib[n]++;
      }
    }
  }
}

/* END set_auto_OBC_us()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Remove cells from the grid that are boundary cells (on the edge   */
/* of the grid) and that don't have a neighbour that is completely   */
/* wet (not on the edge of the grid).                                */
/*-------------------------------------------------------------------*/
void remove_OBC_corner(parameters_t *params)
{
  int **neic;
  int c, cc, cs, cn, i, j, jj;
  int ns2, npe, npem;
  double *bathy = params->bathy;
  int found, wetf;
  int verbose = 0;

  ns2 = params->ns2;
  npem = params->npe;

  /* Get the neighbour mappings                                      */
  neic = i_alloc_2d(ns2+1, npem+1);
  for (cc = 1; cc <= ns2; cc++) {
    npe = params->npe2[cc];
    for (j = 1; j <= npe; j++) {
      jj = j;
      cn = find_neighbour(cc, params->x, params->y, params->npe2, ns2, &jj);
      neic[j][cc] = cn;
    }
  }

  /* Remove corner cells                                             */
  c = 1;
  for (cc = 1; cc <= ns2; cc++) {
    npe = params->npe2[cc];
    for (j = 1; j <= npe; j++) {
      wetf = 1;
      cn = neic[j][cc];
      if (is_wet(params, bathy[cc]) && !cn) {
	wetf = 0;
	/* Look for a wet neighbour                                  */
	for (jj = 1; jj <= npe; jj++) {
	  found = 0;
	  cn = neic[jj][cc];
	  if (cn && is_wet(params, bathy[cn])) {
	    /* Check if cn is on the grid edge                       */
	    for (i = 1; i <= params->npe2[cn]; i++) {
	      cs = neic[i][cn];
	      if (!cs) {
		found = 1;
	      }
	    }
	    if (!found) wetf = 1;
	  }
	}
	if (!wetf) {
	  if (verbose) printf("Removing corner cell at cc=%d (bathy=NOTVALID)\n", cc);
	  break;
	} 
      }
    }
    if (wetf) {
      for (j = 0; j <= npe; j++) {
	params->x[c][j] = params->x[cc][j];
	params->y[c][j] = params->y[cc][j];
      }
      params->bathy[c] = params->bathy[cc];
      params->npe2[c] = params->npe2[cc];
      c++;
    }
  }
  params->ns2 = c - 1;
  i_free_2d(neic);
}

/* END remove_OBC_corner()                                           */
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
  char buf[MAXSTRLEN], keyword[MAXSTRLEN], val[MAXSTRLEN];
  int n, m, cc;

  sprintf(keyword, "WINDOWS");
  prm_read_int(fp, keyword, &params->nwindows);

  /* Window maps: original specification                             */
  sprintf(keyword, "DUMP_WIN_MAP");
  if (prm_read_char(fp, keyword, params->wind_file))
    params->map_type |= WIN_DUMP;
  sprintf(keyword, "READ_WIN_MAP");
  if (prm_read_char(fp, keyword, params->win_file)) {
    params->win_type = WIN_FILE;
    params->map_type |= WIN_READ;
  }
  /* New specification                                               */
  sprintf(keyword, "DUMP_MAP");
  if (prm_read_char(fp, keyword, buf)) {
    if (decode_tag(buf, "WIN", val)) {
      strcpy(params->wind_file, val);
      params->map_type |= WIN_DUMP;
    }
    if (decode_tag(buf, "GEOM", val)) {
      strcpy(params->geom_file, val);
      params->map_type |= GEOM_DUMP;
    }
  }
  sprintf(keyword, "READ_MAP");
  if (prm_read_char(fp, keyword, buf)) {
    if (decode_tag(buf, "WIN", val)) {
      strcpy(params->win_file, val);
      params->map_type |= WIN_READ;
    }
    if (decode_tag(buf, "GEOM", val)) {
      strcpy(params->geom_file, val);
      params->map_type |= GEOM_READ;
    }
    if (decode_tag(buf, "CHECK_G", val)) {
      strcpy(params->geom_file, val);
      params->map_type |= GEOM_CHECK;
    }
    if (decode_tag(buf, "CHECK_W", val)) {
      strcpy(params->win_file, val);
      params->win_type |= WIN_CHECK;
      params->map_type |= WIN_CHECK;
    }
  }

  if (params->nwindows > 1 && params->trasc & LAGRANGE)
    hd_quit("Multiple windows not supported with advection scheme 'LAGRANGE'\n");
  if (params->tmode & (SP_CHECK|TR_CHECK)) params->nwindows = 1;

  if (params->nwindows > 1) {
    read_win_type(params, fp);
    if (params->win_type == WIN_METIS) {
      /* Read options */
      sprintf(keyword, "METIS_OPTIONS");
      params->metis_opts = 0;
      if (prm_read_char(fp, keyword, buf)) {
	if (contains_token(buf, "VOLUME_WEIGHTED") != NULL)
	  params->metis_opts |= METIS_VOLUME_WEIGHTED;
      }
    }
    sprintf(keyword, "CHECK_WIN_MAP");
    if (prm_read_char(fp, keyword, params->win_file)) params->win_type |= WIN_CHECK;

    sprintf(keyword, "SHOW_WINDOWS");
    if (prm_read_char(fp, keyword, buf)) {
      if (is_true(buf)) {
	/*if (params->show_win == 0)*/
	  params->ntrS += 2;
	params->show_win = 1;
      }
    }
    if (params->show_win) params->ntrS += 2;
      
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
    /* Read the explicit window partitions                           */
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
  }

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
    if (strcmp(fields[0], "GROUPED") == 0)
      params->win_type = GROUPED;
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
    if (strcmp(fields[0], "METIS") == 0)
      params->win_type = WIN_METIS;
    if (strcmp(fields[0], "REGION") == 0) {
      params->win_type = WIN_REG;
      strcpy(params->win_file, fields[1]);
    }
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
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n1, n2;
  int  is_vhreg = 0;
  int  is_khreg = 0;

  prm_set_errfn(hd_silent_warn);

  params->u1vh = 0.0;
  params->u1kh = 0.0;

  sprintf(keyword, "U1VH");
  if (prm_read_char(fp, keyword, buf))
    strcpy(params->u1vhc, buf);
  else 
    strcpy(buf, params->u1vhc);
  n1 = parseline(buf, fields, MAXNUMARGS);
  if (n1 == 1) params->u1vh = atof(fields[0]);
  if (n1 > 1 && strcmp(fields[0], "region") == 0) is_vhreg = 1;

  sprintf(keyword, "U1KH");
  if (prm_read_char(fp, keyword, buf))
    strcpy(params->u1khc, buf);
  else 
    strcpy(buf, params->u1khc);
  n2 = parseline(buf, fields, MAXNUMARGS);
  if (n2 == 1) params->u1kh = atof(fields[0]);
  if (n2 > 1 && strcmp(fields[0], "region") == 0) is_khreg = 1;

  sprintf(keyword, "DIFF_SCALE");
  if (prm_read_char(fp, keyword, buf)) {
    params->diff_scale = 0;
    if (contains_token(buf, "NONE") != NULL)
      params->diff_scale |= NONE;
    if (contains_token(buf, "LINEAR") != NULL)
      params->diff_scale |= LINEAR;
    if (contains_token(buf, "NONLIN") != NULL)
      params->diff_scale |= NONLIN;
    if (contains_token(buf, "CUBIC") != NULL)
      params->diff_scale |= CUBIC;
    if (contains_token(buf, "AREA") != NULL)
      params->diff_scale |= AREAL;
    if (contains_token(buf, "AUTO") != NULL)
      params->diff_scale |= (AUTO|NONE);
    if (contains_token(buf, "LENGTH") != NULL)
      params->diff_scale |= H_LEN;
    if (contains_token(buf, "E_AREA") != NULL)
      params->diff_scale |= E_LEN;
    if (contains_token(buf, "C_AREA") != NULL)
      params->diff_scale |= C_LEN;
    if (contains_token(buf, "SMAG") != NULL) {
      params->diff_scale |= SMAG;
      params->smagorinsky = 0.1;
      if (params->u1vh > 0.0) 
	params->u1vh *= -1.0;
      else if (params->u1vh == 0.0) 
	params->u1vh = -1.0;
      if (params->u1kh > 0.0) 
	params->u1kh *= -1.0;
      else if (params->u1kh == 0.0) 
	params->u1kh = -1.0;
    }
    if (contains_token(buf, "SCALE2D") != NULL)
      params->diff_scale |= SCALE2D;
    if (contains_token(buf, "SCALEBI") != NULL)
      params->diff_scale |= SCALEBI;
  }
  if (is_vhreg) {
    params->diff_scale |= VH_REG;
    params->smagorinsky = 1.0;
    params->ntrS += 1;
  }
  if (is_khreg) {
    params->diff_scale |= KH_REG;
    params->smagorinsky = 1.0;
  }

  sprintf(keyword, "SMAGORINSKY");
  if (prm_read_char(fp, keyword, params->smag) ) {
    int n;    
    char *fields[MAXSTRLEN * MAXNUMARGS];
    strcpy(buf, params->smag);
    n = parseline(buf, fields, MAXNUMARGS);
    if (n == 1) {
      params->smagorinsky = atof(fields[0]);
      params->sue1 = 1.0;
      params->kue1 = 1.0;
      params->bsue1 = 0.0;
      params->bkue1 = 0.0;
    }
    if (n == 2) {
      params->bsue1 = floor(atof(fields[0]));
      params->bkue1 = floor(atof(fields[1]));
      params->sue1 = atof(fields[0]) - params->bsue1;
      params->kue1 = atof(fields[1]) - params->bkue1;
      params->smagorinsky = 1.0;
    }
    if (n == 4) {      /* COMPAS compatibility */
      params->bsue1 = floor(atof(fields[0]));
      params->bkue1 = floor(atof(fields[3]));
      params->sue1 = atof(fields[0]) - params->bsue1;
      params->kue1 = atof(fields[3]) - params->bkue1;
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
    int n;    
    char *fields[MAXSTRLEN * MAXNUMARGS];
    if (strcmp(buf, "NONE") == 0)
      params->visc_method = NONE;
    if (strcmp(buf, "LAPLACIAN") == 0)
      params->visc_method = LAPLACIAN;
    if (strcmp(buf, "US_LAPLACIAN") == 0)
      params->visc_method = US_LAPLACIAN;
    if (strcmp(buf, "US_BIHARMONIC") == 0)
      params->visc_method = US_BIHARMONIC;
    if (strcmp(buf, "SIMPLE") == 0)
      params->visc_method = SIMPLE;
    if (strcmp(buf, "PRE_V794") == 0)
      params->visc_method = (LAPLACIAN|PRE794);
    n = parseline(buf, fields, MAXNUMARGS);
    if (n == 3 && strcmp(fields[0], "US_LAPLACIAN") == 0 && strcmp(fields[1], "US_BIHARMONIC") == 0) {
      params->visc_method = US_BIHARMONIC;
      params->visc_fact = atof(fields[2]);
    }
  }
  if (params->compatible & V794)
    params->visc_method |= PRE794;

}

/* END read_hdiff()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read vertical diffusion parameters                     */
/*-------------------------------------------------------------------*/
void read_vdiff(parameters_t *params, FILE *fp, int mode)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int n;
  
  /* Background viscosity                                            */
  sprintf(keyword, "VZ0");
  prm_read_char(fp, keyword, params->vz0c);
  strcpy(buf, params->vz0c);
  n = parseline(buf, fields, MAXNUMARGS);
  if (n == 1) {
    params->vz0 = atof(params->vz0c);
    params->closf |= VZ_C;
  }
  if (n > 1 && strcmp(fields[0], "region") == 0) {
    if (mode)
      params->atr += 1;
    else
      params->ntr += 1;
    params->closf |= VZ_R;
  }

  /* Background diffusivity                                          */
  sprintf(keyword, "KZ0");
  prm_read_char(fp, keyword, params->kz0c);
  strcpy(buf, params->kz0c);
  n = parseline(buf, fields, MAXNUMARGS);
  if (n == 1) {
    params->kz0 = atof(params->kz0c);
    params->closf |= KZ_C;
  }
  if (n > 1 && strcmp(fields[0], "region") == 0) {
    if (mode)
      params->atr += 1;
    else
      params->ntr += 1;
    params->closf |= KZ_R;
  }
    
  prm_set_errfn(hd_silent_warn);
  sprintf(keyword, "KZ_ALPHA");
  prm_read_double(fp, keyword, &params->kz_alpha);
  sprintf(keyword, "VZ_ALPHA");
  prm_read_double(fp, keyword, &params->vz_alpha);
}


/* END read_vdiff()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read wave parameters                                   */
/*-------------------------------------------------------------------*/
void read_waves(parameters_t *params, FILE *fp, int mode)
{
  char buf[MAXSTRLEN];
  int dok = 0;

  /* Read the file name containing the wave information              */
  if (params->do_wave & W_FILE) {
    if (prm_read_char(fp, "WAVE_VARS", params->webf)) {
      if (prm_get_time_in_secs(fp, "WAVE_VARS_INPUT_DT", &params->webf_dt)) {
	params->waves |= ORBITAL;
	params->ntrS += 5;
      }
      if (prm_read_char(fp, "WAVE_VARS_INTERP_TYPE", buf))
	strcpy(params->webf_interp, buf);
      else
	sprintf(params->webf_interp, "%c", '\0');
    }
  }
  if (params->do_wave & W_COMP) {
    params->webf_dt = 0.0;
    params->waves |= ORBITAL;
    params->ntrS += 5;
  }
  if (params->do_wave & W_SWAN) {
    params->webf_dt = 0.0;
    params->waves |= ORBITAL;
    params->ntrS += 5;
  }
  if (params->do_wave & NONE) return;

  /* Read the wave feedback options                                  */
  if (prm_read_char(fp, "WAVES", buf)) {
    if (contains_token(buf, "WAVE_FORCING") != NULL) {
      params->waves |= WAVE_FOR;
      params->ntrS += 2;
      if (params->tendf) params->ntrS += 2;
    }
    if (contains_token(buf, "TAN_RADIATION") != NULL) {
      not_included("tangential radiation stress");
      hd_warn("Use 'WAVE_FORCING' instead!");
      /*
      params->waves |= TAN_RAD;
      params->ntrS += 2;
      if (params->tendf) params->ntrS += 2;
      */
    }
    if (contains_token(buf, "STOKES") != NULL) {
      params->waves |= (STOKES|STOKES_DRIFT|STOKES_MIX);
      params->ntrS += 6;
      /* Note: nearshore and Stokes drift need wavenumber */
      if (!dok) {
	params->ntrS += 1;
	dok = 1;
      }
    }
    if (contains_token(buf, "STOKES_DRIFT") != NULL) {
      params->waves |= (STOKES|STOKES_DRIFT);
      params->ntrS += 6;
      if (!dok) {
	params->ntrS += 1;
	dok = 1;
      }
    }
    if (contains_token(buf, "NEARSHORE") != NULL) {
      params->waves |= NEARSHORE;
      params->ntrS += 16;
      if (!dok) {
	params->ntrS += 1;
	dok = 1;
      }
    }
    if (contains_token(buf, "STOKES_MIX") != NULL) {
      params->waves |= (STOKES|STOKES_MIX);
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
    hd_warn("params_read: WAVES requires DO_WAVES specification. Ignoring WAVES...\n");

  /* Read the fetch file                                             */
  prm_read_char(fp, "FETCH", params->fetch_init);
}

/* END read_waves()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read tracer flux diagnostic                            */
/*-------------------------------------------------------------------*/
void read_trflux(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int i;

  sprintf(keyword, "CALC_FLUXES");
  if (prm_read_char(fp, keyword, buf)) {
    i = parseline(buf, fields, MAXNUMARGS);
    if (strcmp(fields[0], "NONE") != 0) {
      strcpy(params->trflux, fields[0]);
      params->ntr += 4;
      if (i == 3) {
	params->trfd1 = atoi(fields[1]);
	params->trfd2 = atoi(fields[2]);
      }
    }
  } else 
    sprintf(params->trflux, "NONE");
}

/* END read_trflux()                                                 */
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
    params->ntrS+=1;
    params->swr_type |= (SWR_2D|SWR_ATTN);
    /* Check for 2 layer swr attenuation                             */
    sprintf(keyword, "SWR_ATTENUATION_DEEP");
    if (prm_read_char(fp, keyword, params->swr_attn1)) {
      /* Get the partitioning fraction                               */
      sprintf(keyword, "SWR_FRACTION");
      if (!(prm_read_char(fp, keyword, params->swr_tran))) {
	hd_quit_and_dump
	  ("params_read() : Split swr attenuation requires SWR_FRACTION set.\n");
      }
      params->swr_type |= SWR_SPLIT;
      params->ntrS+=2;
    } else {
      /* Check the transmission coefficient                          */
      sprintf(keyword, "SWR_TRANSMISSION");
      if (prm_read_char(fp, keyword, params->swr_tran))
	params->swr_type |= SWR_TRAN;
      else
	strcpy(params->swr_tran, "1.0");
      params->ntrS+=1;
    }
  }

  sprintf(keyword, "SWR_ATTENUATION3D");
  if (prm_read_char(fp, keyword, params->swr_attn)) {
    /* Check the transmission coefficient : default = 1            */
    sprintf(keyword, "SWR_TRANSMISSION");
    if (prm_read_char(fp, keyword, params->swr_tran))
      params->swr_type |= SWR_TRAN;
    else
      strcpy(params->swr_tran, "1.0");
    params->ntrS+=1;
    if (mode)
      params->atr += 1;
    else
      params->ntr += 1;
    params->swr_type |= (SWR_3D|SWR_ATTN);
  }

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
    params->swr_type = (SWR_2D|SWR_ATTN|SWR_TRAN);
    params->ntrS+=2;
  }

  /* SWR parameter estimation                                        */
  sprintf(keyword, "SWR_REGIONS");
  if (prm_read_char(fp, keyword, params->swr_regions)) {
    sprintf(keyword, "SWR_DT");
    prm_get_time_in_secs(fp, keyword, &params->swreg_dt);
    prm_read_char(fp, "SWR_DATA", params->swr_data);
    params->ntrS+=4;
    params->do_lag = 1;
    if (!(params->swr_type & SWR_ATTN)) {
      strcpy(params->swr_attn, "0.0");
      params->swr_type |= SWR_2D;
      params->ntrS+=1;
    }
    if (!(params->swr_type & SWR_TRAN)) {
      if (prm_read_char(fp, "SWR_TRANSMISSION", params->swr_tran))
	params->swr_type |= SWR_TRAN;
      else
	strcpy(params->swr_tran, "0.0");
      if (!(params->swr_type & SWR_SPLIT)) params->ntrS+=1;
    }

    sprintf(keyword, "SWR_ENSEMBLE");
    if (prm_read_char(fp, keyword, params->swr_ensemble)) {
      int i;
      char *fields[MAXSTRLEN * MAXNUMARGS];
      strcpy(buf, params->swr_ensemble);
      i = parseline(buf, fields, MAXNUMARGS);
      if (i == 4) {
	params->swr_ens[0] = atof(fields[0]);
	params->swr_ens[1] = atof(fields[1]);
	params->swr_ens[2] = 1.0;
	params->swr_ens[3] = atof(fields[2]);
	params->swr_ens[4] = atof(fields[3]);
	params->swr_ens[5] = 1.0;
      }
      if (i == 6) {
	params->swr_ens[0] = atof(fields[0]);
	params->swr_ens[1] = atof(fields[1]);
	params->swr_ens[2] = atof(fields[2]);
	params->swr_ens[3] = atof(fields[3]);
	params->swr_ens[4] = atof(fields[4]);
	params->swr_ens[5] = atof(fields[5]);
      }
    }
  } else {
    if (params->swr_type & SWR_ATTN && !(params->swr_type & SWR_TRAN) &&
	!(params->swr_type & SWR_SPLIT)) {
      strcpy(params->swr_tran, "1.0");
      params->ntrS++;
    }
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
  /* Note : the input codes 0 - 4 are included for backwards         */
  /* compatibility.                                                  */
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
/* Routine to read advection scheme parameters                       */
/*-------------------------------------------------------------------*/
void read_advect(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];

  /* Momentum                                                        */
  sprintf(keyword, "MOM_SCHEME");
  if (prm_read_char(fp, keyword, buf)) {
    params->momsc = 0;
    if (contains_token(buf, "ORDER1") != NULL)
      params->momsc |= ORDER1;
    if (contains_token(buf, "ORDER2") != NULL)
      params->momsc |= ORDER2;
    if (contains_token(buf, "VANLEER") != NULL)
      params->momsc |= VANLEER;
    if (contains_token(buf, "WIMPLICIT") != NULL)
      params->momsc |= WIMPLICIT;
    if (contains_token(buf, "EXPLICIT") != NULL)
      params->momsc |= EXPLICIT;
    if (contains_token(buf, "ADVECT_FORM") != NULL)
      params->momsc |= ADVECT_FORM;
    if (contains_token(buf, "WTOP_O4") != NULL) {
      not_included("High order surface vertical velocity");
      hd_quit("Exiting");
      params->momsc |= WTOP_O4;
    }
    if (contains_token(buf, "WTOP_O2") != NULL)
      params->momsc |= WTOP_O2;
    if (contains_token(buf, "ZERO_DRYK") != NULL)
      params->momsc |= ZERO_DRYK;
    if (contains_token(buf, "SHAPIRO") != NULL)
      params->momsc |= SHAPIRO;
    if (contains_token(buf, "LAGRANGE") != NULL)
      params->momsc |= LAGRANGE;
    if (contains_token(buf, "RINGLER") != NULL)
      params->momsc |= RINGLER;
    if (contains_token(buf, "VECINVAR") != NULL)
      params->momsc |= RINGLER;
    if (contains_token(buf, "NEUTRAL") != NULL)
      params->momsc |= PV_ENEUT;
    if (contains_token(buf, "CONSERVE") != NULL)
      params->momsc |= PV_ENSCO;
    if (contains_token(buf, "DISSIPATE") != NULL)
      params->momsc |= PV_ENSDS;
    if (contains_token(buf, "APVM") != NULL)
      params->momsc |= PV_APVM;
    if (contains_token(buf, "LUST") != NULL)
      params->momsc |= PV_LUST;
    if (contains_token(buf, "CLUST") != NULL)
      params->momsc |= PV_CLUST;
    if (!(params->momsc & (PV_ENEUT|PV_ENSCO|PV_ENSDS)))
      params->momsc |= PV_ENEUT;
  }

  /* Use 4th order approximations for wtop by default                */
  if (!(params->momsc & WTOP_O2))
    params->momsc |= WTOP_O4;
  if (params->momsc & SHAPIRO)
    params->filter = ADVECT;

  /* Tracers                                                         */
  sprintf(keyword, "TRA_SCHEME");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "ORIGINAL") == 0)
      params->trasc = ORIGINAL;
    if (strcmp(buf, "ORDER1") == 0)
      params->trasc = ORDER1;
    if (strcmp(buf, "ORDER2") == 0)
      params->trasc = ORDER2;
    if (strcmp(buf, "ORDER4") == 0)
      params->trasc = ORDER4;
    if (strcmp(buf, "QUICKEST") == 0)
      params->trasc = QUICKEST;
    if (strcmp(buf, "QUICKEST|US") == 0)
      params->trasc = QUICKEST|HIORDER;
    if (strcmp(buf, "VANLEER") == 0)
      params->trasc = VANLEER;
    if (strcmp(buf, "VANLEER|US") == 0)
      params->trasc = VANLEER|HIORDER;
    if (strcmp(buf, "FFSL") == 0)
      params->trasc = FFSL;
    if (strcmp(buf, "FCT") == 0 || strcmp(buf, "FCT|ORDER2") == 0)
      params->trasc = FCT|ORDER2;
    if (strcmp(buf, "FCT|ORDER3US") == 0)
      params->trasc = FCT|ORDER3US;
    if (strcmp(buf, "FCT|ORDER4US") == 0)
      params->trasc = FCT|ORDER4US;
    if (strcmp(buf, "ORDER3US") == 0)
      params->trasc = ORDER3US;
    if (strcmp(buf, "ORDER4US") == 0)
      params->trasc = ORDER4US;
    if (strcmp(buf, "ORDER2_UPWIND") == 0)
      params->trasc = ORDER2_UW;
    if (strcmp(buf, "LAGRANGE") == 0)
      params->trasc = LAGRANGE;
    if (strcmp(buf, "LAGRANGE|VANLEER") == 0) {
      params->trasc = LAGRANGE|VANLEER;
      params->trsplit = 1;
    }
    if (strcmp(buf, "FFSL|VANLEER") == 0) {
      params->trasc = FFSL|VANLEER|HIORDER;
      params->trsplit = 1;
    }
  }
  /* ULTIMATE limiter                                                */
  sprintf(keyword, "ULTIMATE");
  if (prm_read_char(fp, keyword, buf))
    params->ultimate = is_true(buf);

  /* Kinetic enegy formulation                                       */
  sprintf(keyword, "KINETIC");
  if (prm_read_char(fp, keyword, buf)) {
    char val[MAXSTRLEN];
    params->kinetic = 0;
    if (decode_tag(buf, "GASSMANN", val)) {
      params->kinetic |= K_GASS;
      params->kfact = atof(val);
    }
    if (decode_tag(buf, "YU", val)) {
      params->kinetic |= K_YU;
      params->kfact = atof(val);
    }
    if (decode_tag(buf, "ORDER2", val))
      params->kinetic |= ORDER2;
    else if (decode_tag(buf, "ORDER1", val))
      params->kinetic |= ORDER1;
    else if (decode_tag(buf, "ORDER4", val))
      params->kinetic |= ORDER4;
    else if (decode_tag(buf, "ORDER2_UW", val))
      params->kinetic |= ORDER2_UW;
    else
      params->kinetic |= ORDER2;
  }

  /* Runge-Kutta stages                                              */
  sprintf(keyword, "RUNGE-KUTTA");
  prm_read_int(fp, keyword, &params->rkstage);
  /* Semi-Lagrange interpolations                                    */
  sprintf(keyword, "ORDER_SL");
  if (prm_read_char(fp, keyword, buf)) {
    if (strcmp(buf, "LINEAR") == 0)
      params->osl = L_LINEAR;
    else if (strcmp(buf, "NN_SIBSON") == 0)
      params->osl = L_SIB;
    else if (strcmp(buf, "NN_NON_SIBSON") == 0)
      params->osl = L_NONSIB;
    else if (strcmp(buf, "CUBIC") == 0)
      params->osl = L_CUBIC;
    else if (strcmp(buf, "QUADRATIC") == 0)
      params->osl = L_LSQUAD;
    else if (strcmp(buf, "LINEARLSQ") == 0)
      params->osl = L_LSLIN;
    else if (strcmp(buf, "BILINEAR") == 0)
      params->osl = L_BILIN;
    else if (strcmp(buf, "BAYCENTRIC") == 0)
      params->osl = L_BAYLIN;
  }
  if (params->trasc & FFSL && params->osl & (L_BILIN|L_BAYLIN))
    hd_quit("Cannot use FFSL advection with bilinear or baycentric interpolation, due to streamline tracking from edges.\n");

  sprintf(keyword, "FILL_METHOD");
  if (prm_read_char(fp, keyword, buf)) {
    if (contains_token(buf, "NONE") != NULL) {
      params->fillf = NONE;
    } else {
      params->fillf = 0;
      if (contains_token(buf, "MONOTONIC") != NULL)
	params->fillf |= MONOTONIC;
      if (contains_token(buf, "CLIP") != NULL)
	params->fillf |= CLIP;
    }
  }
  if (params->fillf & MONOTONIC)
    params->ntrS += 1;
}

/* END read_advect()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read n tidal energy extraction                         */
/*-------------------------------------------------------------------*/
void read_turb(parameters_t *params, FILE *fp)
{
  int n, m;
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];

  if(prm_read_int(fp, "NTURBINES", &params->nturb)) {
    params->turbv = d_alloc_2d(params->nturb, 4);

    params->turbs = (char **)malloc(params->nturb * sizeof(char *));
    for (m = 0; m < params->nturb; m++) {
      params->turbs[m] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      sprintf(keyword, "TURB%d", m);
      prm_read_char(fp, keyword, params->turbs[m]);
      strcpy(buf, params->turbs[m]);
      n = parseline(buf, fields, MAXNUMARGS);
      if (n != 4) hd_quit("read_turb: Format is TURB# x y z v\n");
      params->turbv[0][m] = atof(fields[0]);
      params->turbv[1][m] = atof(fields[1]);
      params->turbv[2][m] = atof(fields[2]);
      params->turbv[3][m] = atof(fields[3]);
    }
  }
}

/* END read_turb()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read histort log specification                         */
/*-------------------------------------------------------------------*/
void read_history(parameters_t *params, FILE *fp)
{
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];
  char hist[MAXSTRLEN];

  if (prm_read_char(fp, "HISTORY", hist)) {
    strcpy(buf, hist);
    if (contains_token(buf, "NONE") != NULL) {
      params->history = NONE;
    } else {
      params->history = 0;
      if (contains_token(buf, "LOG") != NULL)
	params->history |= HST_LOG;
      if (contains_token(buf, "DIFF") != NULL)
	params->history |= HST_DIF;
      if (contains_token(buf, "RESET") != NULL)
	params->history |= HST_RESET;
      if (contains_token(buf, "NOTES") != NULL)
	params->history |= HST_NOTES;
      if (contains_token(buf, "MASTER") != NULL) {
	int n, i;
	char *fields[MAXSTRLEN * MAXNUMARGS];
	strcpy(keyword, hist);
	params->history |= HST_MASTER;
	n = parseline(keyword, fields, MAXNUMARGS);
	for (i = 0; i < n; i++) {
	  if (strcmp(fields[i], "MASTER") == 0) {
	    strcpy(params->histnamem, fields[i+1]);
	    break;
	  }
	}
      }
    }
  }
}

/* END read_history()                                                */
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
    if (contains_token(buf, "V6257") != NULL)
      params->compatible |= V6257;
    if (contains_token(buf, "V6898") != NULL)
      params->compatible |= V6898;
    if (contains_token(buf, "V7367") != NULL)
      params->compatible |= V7367;
    if (contains_token(buf, "V7367") != NULL)
      params->compatible |= V7367;
  }
}

/* END read_compatible()                                             */
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
    ntr = params->ntot = n + 1; /* Always include T & S              */
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
    } else if (strcmp(params->etarlxn, "TPXO") == 0) {
      prm_read_char(fp, tc_name, params->etarlxtcs);
      if (sscanf(params->etarlxtcs, "%lf %s", &params->etarlxtc, buf) == 2)
	tm_scale_to_secs(params->etarlxtcs, &params->etarlxtc);
      if (prm_read_char(fp, dt_name, buf) > 0)
	tm_scale_to_secs(buf, &params->etarlxdt);
      else
	params->etarlxdt = 86400.0;
      params->ntrS++;
      params->etarlx = ETA_TPXO;
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
  int n, m;

  if (prm_read_char(fp, "DEBUG_LOC", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    if (strcmp(fields[0], "us") == 0) {
      params->dbi = atoi(fields[1]);
      params->dbj = atoi(fields[2]);
      m = 3;
    } else {
      if (!(params->us_type & US_IJ)) 
	hd_quit("Use unstructured format for DEBUG_LOC: DEBUG_LOC  us <c> <j>.\n");
      params->dbi = atoi(fields[0]);
      params->dbj = atoi(fields[1]);
      params->dbk = atoi(fields[2]);
      m = 3;
    }
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
/* Routine to read profile generation                                */
/*-------------------------------------------------------------------*/
void read_profile(parameters_t *params, FILE *fp)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  int n, m;

  sprintf(keyword, "PROFILE");
  if (prm_read_char(fp, keyword, buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(buf, fields, MAXNUMARGS);
    if (n == 1) {
      strcpy(params->nprof, fields[0]);
      params->ntr += 1;
    }
    if (n == 2) {
      strcpy(params->nprof, fields[0]);
      strcpy(params->nprof2d, fields[1]);
      params->ntr += 1;
    }
  }
}

/* END read_profile()                                                */
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
/* Routine to read monotonicity diagnostic                           */
/*-------------------------------------------------------------------*/
void read_monotone(parameters_t *params, FILE *fp, int mode)
{
  char buf[MAXSTRLEN];
  char keyword[MAXSTRLEN];
  sprintf(keyword, "MONOTONE");
  if (prm_read_char(fp, keyword, buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int n = parseline(buf, fields, MAXNUMARGS);
    if (n >= 3) {
      strcpy(params->monotr, fields[0]);
      params->monomn = atof(fields[1]);
      params->monomx = atof(fields[2]);
      if (mode)
	params->atr += 1;
      else
	params->ntr += 1;
    } else {
      hd_warn("format: MONOTONE <tracer> <min> <max>\n");
    }
  }
}

/* END read_monotone()                                               */
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
  char fname2[MAXSTRLEN], daystr[MAXSTRLEN], infile2[MAXSTRLEN], buf[MAXSTRLEN];
  double start;            /* Model start time                       */
  double stop;             /* Model stop time                        */
  int nfiles, nf = 0;      /* Number of SST files                    */
  int fid;                 /* netcdf file handle                     */
  int ncerr;               /* netcdf error code                      */
  int i, n, nc;            /* Counter                                */
  int product;             /* Number of days in composite            */
  int ys, mos, ds, hs, mis, ss;  /* Start year, month, day           */
  int ye, moe, de, he, mie, se;  /* End year, month, day             */
  int y, m, d, day, lp, py;      /* Year, month, day                 */
  int is_mnc;
  double ep;                     /* Epoch                            */
  int tday;                      /* Transition day                   */
  int mode = 0;                  /* Transition mode                  */
  int dayno[2][13] = {           /* Days in the month                */
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
  };
  char *hra[24] = {"0000", "0100", "0200", "0300", "0400", "0500", "0600",
		   "0700", "0800", "0900", "1000", "1100", "1200", "1300",
		   "1400", "1500", "1600", "1700", "1800", "1900", "2000",
		   "2100", "2200", "2300"};
  char *monc[13] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
		    "10", "11", "12"};
  int hrc = 0;
  int lyr[2] = {365, 366};
  char files[MAXNUMTSFILES][MAXSTRLEN];
  char fields[MAXNUMTSFILES][MAXSTRLEN];
  char *opts[MAXSTRLEN * MAXNUMARGS];
  char vars[MAXSTRLEN], i_rule[MAXSTRLEN];
  int domon = 0;
  int noday = 0;
  int yrmnc = 0;
  int is_him = 0;
  int urlt = -1;
  int intype;

  /* Get the input arguments                                         */
  nfiles = parseline(params->ghrsst_path, opts, MAXNUMARGS);

  /* Check if all names are .mnc files                               */
  is_mnc = 0;
  for (n = 0; n < nfiles; n++) {
    strcpy(fields[n], opts[n]);
    strcpy(path, fields[n]);
    m = strlen(path);
    if (path[m-4] == '.' && path[m-3] == 'm' &&
	path[m-2] == 'n' && path[m-1] == 'c') {
      is_mnc += 1;
      intype = 1;
    } else { /* File may contain variable substitution */
      for (i = 0; i < m-6; i++) {
	if (path[i] == '.' && path[i+1] == 'm' &&
	    path[i+2] == 'n' && path[i+3] == 'c' &&
	    path[i+4] == '(' && path[m-1] == ')') {
	  is_mnc += 1;
	  intype = 2;
	}
      }
    }
  }
  is_mnc = (is_mnc == nfiles) ? 1 : 0;

  /* Check for shortcuts (sst_url in globals.h)                      */
  if (!is_mnc && nfiles == 1) {
    i = 0;
    while (sst_url[i][0] != NULL) {
      if (strcmp(fields[0], sst_url[i][0]) == 0) {
	strcpy(params->ghrsst_name, sst_url[i][1]);
	strcpy(params->ghrsst_dt, sst_url[i][2]);
	strcpy(files[0], sst_url[i][3]);
	strcpy(fields[0], files[0]);
	strcpy(files[1], sst_url[i][4]);
	strcpy(fields[1], files[1]);
	nfiles = 2;
	if (sst_url[i][5] != NULL)  {
	  strcpy(files[2], sst_url[i][5]);
	  strcpy(fields[2], files[2]);
	  nfiles += 1;
	}
	if (sst_url[i][6] != NULL)  {
	  strcpy(params->ghrsst_opt, sst_url[i][6]);
	}
	if (i == 4 || i == 5) is_him = 1;
	urlt = i;
	if (i == 0 || i == 1) params->ghrsst_type = G_OSTIA;
	if (i == 2) params->ghrsst_type = G_BOM;
	if (i == 3 || i == 4 || i == 5) params->ghrsst_type = G_HIM;
	hd_warn("Satellite import from %s\n", sst_url[i][2]);
	break;
      }
      i++;
    }
  }

  /* Get any options                                                 */
  sprintf(vars, "%c", '\0');
  if (strlen(params->ghrsst_opt)) {
    nc = parseline(params->ghrsst_opt, opts, MAXNUMARGS);
    for (n = 0; n < nc; n++) {
      strcpy(path, opts[n]);
      if (strcmp(path, "DOMON") == 0) domon = 1;
      if (strcmp(path, "NODAY") == 0) noday = 1;
      if (strcmp(path, "YRMNC") == 0) yrmnc = 1;
      if (strcmp(path, "VARIABLES") == 0) {
	strcpy(vars, opts[n+1]);
	n++;
      }
    }
  }

  if (is_mnc) {
    for (n = 0; n < nfiles; n++) {
      strcpy(path, fields[n]);
      if (is_mnc == 1) {
	if (intype == 1) {    /* Assume it is a GHRSST file          */
	  if (n == nfiles-1)
	    sprintf(fname, "%s(ghrsst=analysed_sst)(ghrsst_error=analysis_error)", fields[n]);
	  else
	    sprintf(fname, "%s(ghrsst=analysed_sst)(ghrsst_error=analysis_error) ", fields[n]);
	  if ((ap = fopen(path, "r")) == NULL)
	    hd_warn("Can't open GHRSST file '%s'\n", path);
	  else {
	    fclose(ap);
	  }
	} else if (intype == 2) {
	  /* Filename includes variable substitution                 */
	  if (n == 0) 
	    sprintf(fname, "%s", fields[n]);
	  else
	    sprintf(fname, "%s%s", fname, fields[n]);
	}
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
    } else if (nfiles == 3) {
      strcpy(path, fields[0]);
      strcpy(fname, fields[1]);
      strcpy(fname2, fields[2]);
      mode = 3;
    } else {
      hd_warn("create_ghrsst_list: GHRSST input format = <path> <filename>\n");
      return;
    }
  }

  /* Open the file list                                              */
  params->ghrsst = 1;

  /* Get start and end year, month, day                              */
  ep = tm_time_to_julsecs(params->timeunit);
  tm_scale_to_secs(params->start_time, &start);
  tm_scale_to_secs(params->stop_time, &stop);
  nfiles = (ceil(stop) - floor(start)) / 86400;
  if (is_him) nfiles *=24;
  start = start / 86400 + ep;;
  tm_to_julsecs(start, &ys, &mos, &ds, &hs, &mis, &ss);
  stop = stop / 86400.0 + ep;
  tm_to_julsecs(stop, &ye, &moe, &de, &he, &mie, &se);
  y = py = ys; m = mos; d = ds; n = 0;

  /* Print the file header                                           */
  if (yrmnc) {
    sprintf(buf, "ghrsst_list-%d.mnc", y);
    strcpy(files[yrmnc], buf);
  } else
    sprintf(buf, "ghrsst_list.mnc");
  ap = fopen(buf, "w");
  fprintf(ap, "multi-netcdf-version 1.0\n");
  fprintf(ap, "nfiles %d\n\n", nfiles);

  /* Make a list of files for each day.                              */
  /* Print the list of files                                         */
  for (i = 0; i < nfiles; i++) {
    day = yrday(y, m, d);
    /* Separate .mnc files for each year if required                 */
    if (yrmnc && y != py) {
      py = y;
      rewind(ap);
      fprintf(ap, "multi-netcdf-version 1.0\n");
      fprintf(ap, "nfiles %d\n", n);
      fclose(ap);
      sprintf(buf, "ghrsst_list-%d.mnc", y);
      ap = fopen(buf, "w");
      fprintf(ap, "multi-netcdf-version 1.0\n");
      fprintf(ap, "nfiles %d\n\n", nfiles);
      yrmnc++;
      strcpy(files[yrmnc], buf);
      n = 0;
    }
    if (day < 10)
      sprintf(daystr, "00%d", day);
    else if (day < 100)
      sprintf(daystr, "0%d", day);
    else
      sprintf(daystr, "%d", day);
    if (!noday) {
      if (domon) {
	if (path[strlen(path) - 1] == '/')
	  sprintf(infile, "%s%d/%s", path, y, monc[m]);
	else
	  sprintf(infile, "%s/%d/%s/", path, y, monc[m]);
      } else {
	if (path[strlen(path) - 1] == '/')
	  sprintf(infile, "%s%d/%s", path, y, daystr);
	else
	  sprintf(infile, "%s/%d/%s", path, y, daystr);
      }
    } else {
      if (path[strlen(path) - 1] == '/')
	sprintf(infile, "%s%d", path, y);
      else
	sprintf(infile, "%s/%d/", path, y);
    }

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
    if (!is_him) d += 1;
    if (d > dayno[lp][m]) {
      m++;
      d = 1;
    }
    if (d > lyr[lp] || m > 12) {
      y++;
      m = 1;
      d = 1;
    }
    if (mode == 1 && !is_him)
      sprintf(infile, "%s/%s%s", infile, date, fname);
    else if (mode == 3 || (mode == 1 && is_him)) {
      /*
      timeseries_t *ts ;
      ts = (timeseries_t *)malloc(sizeof(timeseries_t));
      */
      strcpy(infile2, infile);
      if (is_him) {
	sprintf(infile, "%s/%s%s%s", infile, date, hra[hrc], fname);
	if (hrc == 23) {
	  d += 1;
	  hrc = 0;
	} else
	  hrc += 1;
      } else
	sprintf(infile, "%s/%s%s", infile, date, fname);
      if ((ncerr = nc_open(infile, NC_NOWRITE, &fid)) != NC_NOERR) {
	printf("file%d %s not found\n",n, infile);
	sprintf(infile2, "%s/%s%s", infile2, date, fname2);
	if ((ncerr = nc_open(infile2, NC_NOWRITE, &fid)) != NC_NOERR) {
	  continue;
	  /*hd_quit("Can't open file%d %s\n", n, infile2);*/
	}
	strcpy(infile, infile2);
      }
      nc_close(fid);
    } else {
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

  /* Write the .mnc files to params->ghrsst_path                     */
  if (yrmnc) {
    if (strlen(vars))
      sprintf(params->ghrsst_path, "%s%s", files[1], vars);
    else
      sprintf(params->ghrsst_path, "%s(ghrsst=analysed_sst)(ghrsst_error=analysis_error)", files[1]);
    for (n = 2; n <= yrmnc; n++) {
      if (strlen(vars))
	sprintf(params->ghrsst_path, "%s %s%s", params->ghrsst_path, files[n], vars);
      else
	sprintf(params->ghrsst_path, "%s %s(ghrsst=analysed_sst)(ghrsst_error=analysis_error)", 
	       params->ghrsst_path, files[n]);
    }
  } else {
    if (strlen(vars))
      sprintf(params->ghrsst_path, "ghrsst_list.mnc%s", vars);
    else
      strcpy(params->ghrsst_path, "ghrsst_list.mnc(ghrsst=analysed_sst)(ghrsst_error=analysis_error)");
  }
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

    not_implemented("explicit mapping");
    return;

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

    not_implemented("explicit mapping");
    return;

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

  read_blocks(fp, "EXCLUDE_BGCSED", &params->prex, &params->prxi, &params->prxj, NULL);
  if (params->prex) {
    params->prxf = i_alloc_1d(params->prex);
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
  char meanbuf[MAXSTRLEN];

  sprintf(keyword, "MEAN");
  params->means = 0;
  sprintf(params->means_dt, "%c", '\0');
  if (prm_read_char(fp, keyword, buf)) {
    strcpy(meanbuf, buf);
    if (contains_token(buf, "NONE") != NULL)
      params->means |= NONE;
    if (contains_token(buf, "VEL3D") != NULL) {
      params->means |= VEL3D;
      strip(meanbuf, "VEL3D");
    }
    if (contains_token(buf, "VEL2D") != NULL) {
      params->means |= VEL2D;
      strip(meanbuf, "VEL2D");
    }
    if (contains_token(buf, "FLUX") != NULL) {
      params->means |= FLUX;
      strip(meanbuf, "FLUX");
    }
    if (contains_token(buf, "TENDENCY") != NULL) {
      params->means |= TENDENCY;
      strip(meanbuf, "TENDENCY");
    }
    if (contains_token(buf, "WIND") != NULL) {
      params->means |= WIND;
      strip(meanbuf, "WIND");
    }
    if (contains_token(buf, "ETA") != NULL) {
      params->means |= ETA_M;
      strip(meanbuf, "ETA");
    }
    if (contains_token(buf, "KZ_M") != NULL) {
      params->means |= KZ_M;
      strip(meanbuf, "KZ_M");
    }
    if (contains_token(buf, "TS") != NULL) {
      params->means |= TS;
      strip(meanbuf, "TS");
    }
    if (contains_token(buf, "TS|MMM") != NULL) {
      params->means |= (TS|MMM);
      strip(meanbuf, "TS|MMM");
    }
    if (contains_token(buf, "TIDAL") != NULL) {
      params->means |= TIDAL;
      strip(meanbuf, "TIDAL");
    }
    if (contains_token(buf, "TRANSPORT") != NULL) {
      params->means |= (TRANSPORT|VEL3D|VOLFLUX|KZ_M);
      strip(meanbuf, "TRANSPORT");
    }
    if (contains_token(buf, "VOLFLUX") != NULL) {
      params->means |= VOLFLUX;
      strip(meanbuf, "VOLFLUX");
    }
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
    /* Generic (auto) tracers                                        */
    if (contains_token(buf, "TRA3D") != NULL) {
      params->means |= MTRA3D;
      strip(meanbuf, "TRA3D");
      strcpy(params->means_tra, meanbuf);
    }
    if (contains_token(buf, "TRA2D") != NULL) {
      params->means |= MTRA2D;
      strip(meanbuf, "TRA2D");
      strcpy(params->means_tra, meanbuf);
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
    /* u1vmean is now a state variable rather than a tracer
    if (params->means & VOLFLUX) {
      if(mode)
	params->atr += 2;
      else
	params->ntr += 2;
    }
    */
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
      /*
      if(mode)
	params->atr += 2;
      else
	params->ntr += 2;
      */
    }
  }
  /*
  if (params->tmode & SP_FFSL) {
    if(mode)
      params->atr += 2;
    else
      params->ntr += 2;
  }
  */
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
      if (strlen(params->ecosedconfig))
	strcpy(params->sed_defs, params->ecosedconfig);
      else
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
      if (strlen(params->ecosedconfig))
	strcpy(params->eco_defs, params->ecosedconfig);
      else
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
		/*
		if (nf <= 2) {
		  params->rsalt = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		*/
		if (strcmp(fields[0], "obc") == 0) {
		  params->rsalt = RLX_OBC;
		  params->atr += 2;
		  params->ntr += 2;
		} else if (strcmp(fields[0], "region") == 0) {
		  params->rsalt = RLX_REG;
		  params->atr += 2;
		  params->ntr += 2;
		} else if (endswith(fields[0], ".nc")) {
		  params->rsalt = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		} else { /* Adaptive relaxation */
		  params->rsalt = RLX_ADPT;
		  params->atr += 2;
		  params->ntr += 2;
		}
	      }
	      if (strcmp(name, "temp") == 0) {
		/*
		if (nf <= 2) {
		  params->rtemp = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		*/
		if (strcmp(fields[0], "obc") == 0) {
		  params->rtemp = RLX_OBC;
		  params->atr += 2;
		  params->ntr += 2;
		} else if (strcmp(fields[0], "region") == 0) {
		  params->rtemp = RLX_REG;
		  params->atr += 2;
		  params->ntr += 2;
		} else if (endswith(fields[0], ".nc")) {
		  params->rtemp = RLX_FILE;
		  params->atr += 1;
		  params->ntr += 1;
		} else { /* Adaptive relaxation */
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

  /* Check for duplicate tracers (i.e. automatically generated       */
  /* tracers that are also defined in the input parameter file).     */
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
    if (params->waves & STOKES_DRIFT) {
      if (tracer_find_index("u1_sto", params->ntr, params->trinfo_3d) >= 0)
	ntr++;
      if (tracer_find_index("u2_sto", params->ntr, params->trinfo_3d) >= 0)
	ntr++;
    }
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
  if (params->save_force & FTEMP) {
    if (tracer_find_index("temp_force", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & FSALT) {
    if (tracer_find_index("salt_force", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & FVELU) {
    if (tracer_find_index("velu_force", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (params->save_force & FVELV) {
    if (tracer_find_index("velv_force", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->mixsc, "k-w") == 0 || strcmp(params->mixsc, "W88") == 0) {
    if (tracer_find_index("tke", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
    if (tracer_find_index("omega", params->ntr, params->trinfo_3d) >= 0)
      ntr++;
  }
  if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
      strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0) {
    if (tracer_find_index("lscale", params->ntr, params->trinfo_3d) >= 0)
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
  /* Re-order tracers defined in the input parameter file to account */
  /* for duplicate tracers.                                          */

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
/* Reads in a list of mesh indicies as a series of blocks            */
/*-------------------------------------------------------------------*/
int read_blocks(FILE *fp, char *key, int *nb, int **listi, int **listj, char **pname)
{
  char buf[MAXSTRLEN], pbuf[MAXSTRLEN];;
  int nlist, n, nn;
  int i, j, is, js, ie, je;
  int *xi, *xj;
  int ret = NONE;

  if (prm_read_int(fp, key, &nlist)) {
    *nb = 0;
    for (n = 0; n < nlist; n++) {
      prm_next_line(buf, MAXSTRLEN, fp);
      if (sscanf(buf, "(%d,%d)-(%d,%d)",&is, &js, &ie, &je) == 4) {
	ret = B_BLOCKIJ;
	for (j = js; j <= je; j++)
	  for (i = is; i <= ie; i++)
	    *nb += 1;
      } else if (sscanf(buf, "(%d)-(%d)",&is, &ie) == 2) {
	ret = B_BLOCKC;
	for (i = is; i <= ie; i++)
	  *nb += 1;
      } else if (sscanf(buf, "%d %d",&is, &js) == 2) {
	ret = B_LISTIJ;
	*nb += 1;
      } else if (sscanf(buf, "%d",&is) == 1) {
	ret = B_LISTC;
	*nb += 1;
      } else {
	if (pname != NULL) {
	  FILE *fp;
	  sscanf(buf,"%s",pbuf);
	  if ((fp = fopen(pbuf, "r")) == NULL)
	    hd_quit("Can't open file '%s' for %s\n", pbuf, key);
	  strcpy(pname[n+1], pbuf);
	  ret = B_POLY;
	  *nb += 1;
	}
      }
    }

    if (ret == B_POLY) return(0);
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
      } else if (sscanf(buf, "(%d)-(%d)",&is, &ie) == 2) {
	for (i = is; i <= ie; i++) {
	  xi[nn] = i;
	  nn++;
	}
      } else if (sscanf(buf, "%d %d",&is, &js) == 2) {
	xi[nn] = is;
	xj[nn] = js;
	nn++;
      } else if (sscanf(buf, "%d",&is) == 1) {
	xi[nn] = is;
	nn++;
      }
    }
    *listi = xi;
    *listj = xj;
  } else
    *nb = 0;
  return(ret);
}

/* END read_blocks()                                                 */
/*-------------------------------------------------------------------*/

#define R_OUT  1
#define R_INT  2
#define R_OBC  4
#define R_NOR  8
#define R_TAN  16

/*-------------------------------------------------------------------*/
/* Cookie cuts a mesh given a .bncc region file.                     */
/* Format:                                                           */
/* COOKIE_CUT file.bncc t1 t2 ... tn tname                           */
/* makes .prm file tnamet1.prm, tnamet2.prm etc.                     */
/* file.bncc is the region file specifying the tiles.                */
/* t1 t2 ... tn are the regions in file.bncc to tile.                */
/* For 2 way nesting:                                                */
/* COOKIE_CUT file.bncc t1 t2 ... tn tname sf bf 2way                */
/* sf = separation interface (cells)                                 */
/* bf = 0 (NO) no barotropic coupling, bf = 1 (YES) for barotropic   */
/*      coupling.                                                    */
/*-------------------------------------------------------------------*/
void cookie_cut(master_t *master, parameters_t *params)
{
  FILE *op, *bp;
  geometry_t *geom = master->geom;
  mesh_t *mesh = params->mesh;
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int c, cc, cn, ee, e, v;
  int j, jp, n, nn, nb, m, mr, nr, fn, nv, rgn, ns2, ns;
  int nreg, *reg, *mask, *maskb, *mapc, *mapv, **flag, *rask, *eask;
  int *nobc, *nobn, *nobt, *obc, *obe, tobc;
  int found[geom->nobc];
  double *regionid;
  int barof = 0;          /* Barotropic coupling for 2 way           */
  int sf = 0;             /* Interface separation for 2 way          */
  int verbose = 0;

  if (!strlen(params->cookiecut)) return;
  n = parseline(params->cookiecut, fields, MAXNUMARGS);
  if (strcmp(fields[n-1], "2way") == 0) {
    sf = 2;
    mr = n - 5;           /* Number of tiles to produce              */
    fn = n - 4;           /* Index of file name argument             */
    /*barof = atoi(fields[n-2]);*/
    barof = is_true(fields[n-2]);
    sf = atoi(fields[n-3]);
  } else {
    mr = n - 2;
    fn = n - 1;
  }

  if (n < 1) hd_quit("COOKIE_CUT format = region.bncc region_no output.prm\n");
  if (!(endswith(fields[0], ".bncc")))
    hd_quit("COOKIE_CUT requires a .bncc (COMPAS) region file\n");
  regionid = d_alloc_1d(geom->szc);
  nr = read_regioni(master, fields[0], regionid);
  if (!nr) {
    hd_warn("cookie_cut(): Region file %s contains no valid boxnos\n", fields[0]);
    return;
  }
  rask = i_alloc_1d(nr);
  memset(rask, 0, nr * sizeof(int));
  for (m = 0; m < mr; m++) {
    rgn = atoi(fields[m+1]);
    rask[rgn] = 1;
  }

  /* Loop over supplied regions and create .prm files. Parameters    */
  /* are written according to current .prm settings.                 */
  for (m = 0; m < mr; m++) {
    /* Open the output file                                          */
    sprintf(params->oname, "%s%s", fields[fn], fields[m+1]);
    sprintf(buf, "%s.prm", params->oname);

    /* Get the region to cookie cut and save cell centres            */
    reg = i_alloc_1d(geom->szcS);
    mask = i_alloc_1d(geom->szcS);
    rgn = atoi(fields[m+1]);
    memset(mask, 0, geom->szcS * sizeof(int));
    nreg = 0;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      if (rgn == (int)regionid[c]) {
	reg[nreg++] = c;
	mask[c] = 1;
	if (sf) {
	  for (j = 1; j <= geom->npe[c]; j++) {
	    cn = geom->c2c[j][c];
	    n = (int)regionid[cn];
	    if (!mask[cn] && n != NOTVALID && n != rgn) {
	      reg[nreg++] = cn;
	      mask[cn] = 1;
	    }
	  }
	}
      }
    }

    /* Extend the interface separation                               */
    for (ns = 1; ns < sf; ns++) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	if (rgn != (int)regionid[c] && mask[c] == ns) {
	  for (j = 1; j <= geom->npe[c]; j++) {
	    cn = geom->c2c[j][c];
	    n = (int)regionid[cn];
	    if (!mask[cn] && n != NOTVALID && n != rgn && !geom->wgst[cn]) {
	      reg[nreg++] = cn;
	      mask[cn] = ns+1;
	    }
	  }
	}
      }
    }
    i_free_1d(mask);

    /* Get the edges of the cell centres                             */
    nv = 0;
    mask = i_alloc_1d(geom->szvS);
    memset(mask, 0, geom->szvS * sizeof(int));
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      for (j = 1; j <= geom->npe[c]; j++) {
	v = geom->c2v[j][c];
	if (!mask[v]) {
	  mask[v] = v;
	  nv++;
	}
      }
    }
  
    /* Reset mesh size and write the parameter file                  */
    ns2 = mesh->ns2;
    ns = mesh->ns;
    mesh->ns2 = nreg;
    mesh->ns = nreg + nv;
    params_write(params, master->dumpdata, params->oname);
    if ((op = fopen(buf, "a")) == 0)
      hd_quit("cookie_cut(): Can't open file %s\n", buf);
    fseek(op, 0, SEEK_END);
    if (params->us_type & US_IJ) {
      fprintf(op, "\nGRIDTYPE             UNSTRUCTURED\n");
      fprintf(op, "Mesh2 unstructured   v1.0\n");
      fprintf(op, "nMaxMesh2_face_nodes %d\n", mesh->mnpe);
      fprintf(op, "nMesh2_face_indices  %d\n", mesh->ns);
      fprintf(op, "nMesh2_face          %d\n\n", mesh->ns2);
    }
    mesh->ns2 = ns2;
    mesh->ns = ns;

    /* Make maps from vertices and centres to mesh indices           */
    n = 1;
    mapc = i_alloc_1d(geom->szcS);
    mapv = i_alloc_1d(geom->szvS);
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      mapc[c] = n;
      n++;
    }
    for (cc = 1; cc <= geom->n2_e2; cc++) {
      v = geom->w2_e2[cc];
      if (v == mask[v]) {
	mapv[v] = n;
	n++;
      }
    }
    if(n - 1 != nv + nreg) hd_quit("cookie cut error.\n");

    /* Write the dumpfiles for data transfer across tiles.           */
    /* First make a flag array describing the status of cells and    */
    /* count the R_OBC cells for each boundary.                      */
    fprintf(op, "OutputFiles           df_tile%d.txt\n\n",rgn);
    flag = i_alloc_2d(nreg, geom->npem+1);
    maskb = i_alloc_1d(geom->szcS);
    memset(maskb, 0, geom->szcS * sizeof(int));
    eask = i_alloc_1d(geom->szeS);
    memset(eask, 0, geom->szeS * sizeof(int));
    nobc = i_alloc_1d(nr);
    nobn = i_alloc_1d(nr);
    nobt = i_alloc_1d(nr);
    memset(nobc, 0, nr * sizeof(int));
    memset(nobn, 0, nr * sizeof(int));
    memset(nobt, 0, nr * sizeof(int));
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      n = (int)regionid[c];
      maskb[c] = n+1;
      for (j = 1; j <= geom->npe[c]; j++) flag[j][cc] = 0;
    }

    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      n = (int)regionid[c];
      /* Centres                                                     */
      if (n == rgn)
	flag[0][cc] = R_INT;
      else {
	flag[0][cc] = R_OBC;
	nobc[n]++;
      }
      /* Edges                                                       */
      for (j = 1; j <= geom->npe[c]; j++) {
	cn = geom->c2c[j][c];
	if (geom->wgst[cn]) continue;
	nn = (int)regionid[cn];
	if (nn == NOTVALID) 
	  flag[j][cc] = R_OUT;
	else {
	  if (nn == rgn) {
	    flag[j][cc] = R_INT;
	  } else {
	    if ((sf && maskb[cn] && maskb[c]-1 != rgn && maskb[cn]-1 != rgn)) {
	      int e = geom->c2e[j][c];
	      if (!eask[e]) {
		flag[j][cc] = (R_OBC|R_TAN);
		eask[e] = 1;
		nobt[nn]++;
	      }
	    } 
	    if (!maskb[cn]) {
	      flag[j][cc] = (R_OBC|R_NOR);
	      nobn[nn]++;
	    }
	  }
	}
      }
    }

    tobc = 0;
    for (n = 0; n < nr; n++) {
      if (!rask[n]) nobc[n] = nobn[n] = nobt[n] = 0;
      if (nobc[n]) {
	if (verbose) printf("OBC%d-%d = %d centres\n",rgn,n,nobc[n]);
	tobc++;
      }
    }

    if (sf) {
      /* Print the output from this tile used by other tiles         */
      sprintf(key, "df_tile%d.txt",rgn);
      bp = fopen(key, "w");
      fprintf(bp, "multi-dumpfile-version 1.0\n");
      fprintf(bp, "nfiles                %d\n",tobc);
      c = 0;
      for (n = 0; n < nr; n++) {
	if (nobc[n]) fprintf(bp, "file%d df_tile%d-%d.txt\n",c++, rgn, n);
      }
      fclose(bp);

      for (n = 0; n < geom->nobc; n++) {
	open_bdrys_t *open = geom->open[n];
	found[n] = 0;
	for (ee = 1; ee <= open->no2_e1; ee++) {
	  c = open->obc_e2[ee];
	  ns = (int)regionid[c];	
	  if (ns == rgn) found[n]++;
	}
      }
      for (n = 0; n < geom->nobc; n++) {
	if (found[n]) tobc++;
      }
    } else {
      tobc = 0;
      memset(nobc, 0, nr * sizeof(int));
    }

    /* Print the open boundaries for this tile                       */
    fprintf(op, "#TIDE_CSR_CON_DIR        /home/cem/tides/jhtide2/\n");
    fprintf(op, "#TIDE_CSR_ORTHOWEIGHTS   /home/cem/tides/csr4.0/ortho_csr_4.0\n");
    fprintf(op, "#TIDE_CONSTITUENTS       /home/cem/tides/otps/consts_otps_zUV_AUS.nc\n\n");

    fprintf(op, "NOTE: Cross-shelf boundary may require specification.\n");
    fprintf(op, "NBOUNDARIES           %d\n",tobc);
    nb = 0;
    /* Open boundaries that fall within this tile                    */
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (found[n]) {
	fprintf(op, "BOUNDARY%d.NAME        %s\n", nb, open->name);
	fprintf(op, "BOUNDARY%d.TYPE        u1\n", nb);
	c = open->obc_t[1];
	fprintf(op, "#BOUNDARY%d.START_LOC   %f %f\n",nb, geom->cellx[c], geom->celly[c]);
	c = open->obc_t[open->no2_t];
	fprintf(op, "#BOUNDARY%d.END_LOC     %f %f\n",nb, geom->cellx[c], geom->celly[c]);
	c = open->obc_t[open->no2_t/2];
	fprintf(op, "#BOUNDARY%d.MID_LOC     %f %f\n",nb, geom->cellx[c], geom->celly[c]);
	fprintf(op, "BOUNDARY%d.UPOINTS     %d\n",nb, found[n]);
	for (ee = 1; ee <= open->no2_e1; ee++) {
	  c = open->obc_e2[ee];
	  e = open->obc_e1[ee];
	  ns = (int)regionid[c];	
	  if (ns == rgn) {
	    int v1 = geom->e2v[e][0];
	    int v2 = geom->e2v[e][1];
	    fprintf(op, "%d (%5.3f,%5.3f)-(%5.3f,%5.3f)\n",mapc[c],
		    geom->gridx[v1], geom->gridy[v1],
		    geom->gridx[v2], geom->gridy[v2]);
	  }
	}
	nb++;
      }
      fprintf(op, "\n");
    }
    /* Open boundaries from overlapping tiles                        */
    for (n = 0; n < nr; n++) {
      if (nobc[n]) {
	int obc = 0;
	double sx, sy, ex, ey, mx, my;
	sprintf(key, "boundary%d-%d.xy", n, rgn);
	if (verbose) bp = fopen(key, "w");
	fprintf(op, "BOUNDARY%d.NAME        tile%d\n", nb, n);
	fprintf(op, "BOUNDARY%d.TYPE        u1\n", nb);
	sprintf(key, "bdry%d-%d_", n, rgn);
	fprintf(op, "BOUNDARY%d.BCOND0      SOLID\n", nb);
	if (barof) {
	  fprintf(op, "#BOUNDARY%d.BCOND0      NEST2WAY %seta.mpk %sts.mpk %suv_nor.mpk %suv_tan.mpk\n", 
		  nb, key, key, key, key);
	  fprintf(op, "BOUNDARY%d.BCOND_nor     FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_nor2d   FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_tan     FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_tan2d   FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_ele     FILEIN\n", nb);
	  fprintf(op, "#BOUNDARY%d.ADJUST_RATIO  1.2\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_salt    FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_temp    FILEIN\n", nb);
	  fprintf(op, "#BOUNDARY%d.FILEIN_DT     %5.2f seconds\n", nb, master->grid_dt);
	  fprintf(op, "BOUNDARY%d.OPTIONS       TILED\n", nb);
	  fprintf(op, "BOUNDARY%d.INVERSE_BAROMETER NO\n", nb);
	  fprintf(op, "BOUNDARY%d.data          %seta.mpk %sts.mpk %suv_nor.mpk %suv_tan.mpk %suvav_nor.mpk %suvav_tan.mpk\n", nb, key, key, key, key, key, key);
	}
	else {
	  fprintf(op, "#BOUNDARY%d.BCOND0      NEST2WAY %sets.mpk %suv_nor.mpk %suv_tan.mpk\n", 
		  nb, key, key, key);
	  fprintf(op, "BOUNDARY%d.BCOND_nor     FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_tan     FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_ele     NOTHIN|FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.ADJUST_RATIO  1.2\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_salt    FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.BCOND_temp    FILEIN\n", nb);
	  fprintf(op, "BOUNDARY%d.FILEIN_DT     %5.2f seconds\n", nb, master->grid_dt);
	  fprintf(op, "BOUNDARY%d.OPTIONS       TILED\n", nb);
	  fprintf(op, "BOUNDARY%d.INVERSE_BAROMETER NO\n", nb);
	  fprintf(op, "BOUNDARY%d.data          %seta.mpk %sts.mpk %suv_nor.mpk %suv_tan.mpk\n", 
		  nb, key, key, key, key);
	}
	fprintf(op, "BOUNDARY%d.UPOINTS     %d\n",nb, nobn[n]);
	for (cc = 0; cc < nreg; cc++) {
	  c = reg[cc];
	  ns = (int)regionid[c];
	  if (ns == n) {
	    for (j = 1; j <= geom->npe[c]; j++) {
	      int e = geom->c2e[j][c];
	      int v1 = geom->e2v[e][0];
	      int v2 = geom->e2v[e][1];
	      cn = geom->c2c[j][c];
	      if (flag[j][cc] & R_NOR) {
		fprintf(op, "%d (%5.3f,%5.3f)-(%5.3f,%5.3f)\n",mapc[c],
			geom->gridx[v1], geom->gridy[v1],
			geom->gridx[v2], geom->gridy[v2]);
		if (obc == 0) {
		  sx = geom->cellx[c];
		  sy = geom->celly[c];
		}
		if (obc == 1) {
		  mx = geom->cellx[c];
		  my = geom->celly[c];
		}
		ex = geom->cellx[c];
		ey = geom->celly[c];
		obc++;
		if (verbose) {
		  fprintf(bp, "%f %f\n",geom->gridx[v1], geom->gridy[v1]);
		  fprintf(bp, "%f %f\n",geom->gridx[v2], geom->gridy[v2]);
		}
	      }
	    }
	  }
	}
	fprintf(op, "#BOUNDARY%d.START_LOC   %f %f\n",nb, sx, sy);
	fprintf(op, "#BOUNDARY%d.END_LOC     %f %f\n",nb, ex, ey);
	fprintf(op, "#BOUNDARY%d.MID_LOC     %f %f\n",nb, mx, my);
	fprintf(op, "\n");
	if (verbose) {
	  fprintf(bp, "NaN NaN\n");
	  fclose(bp);
	}

	/* Print the output file specification from tile n used by   */
	/* this tile.                                                */
	sprintf(key, "df_tile%d-%d.txt",n,rgn);
	bp = fopen(key, "w");
	ns = 0;
	if (barof) {
	  fprintf(bp, "OutputFiles           6\n\n");
	  write_tile_dump(master, bp, "eta", "eta", n, rgn, regionid, nreg, reg, flag, nobc[n], ns++, R_OBC);
	  sprintf(key, "temp salt");
	  write_tile_dump(master, bp, "ts", key, n, rgn, regionid, nreg, reg, flag, nobc[n], ns++, R_OBC);
	} else {
	  fprintf(bp, "OutputFiles           3\n\n");
	  sprintf(key, "eta temp salt");
	  write_tile_dump(master, bp, "ets", key, n, rgn, regionid, nreg, reg, flag, nobc[n], ns++, R_OBC);
	}
	sprintf(key, "u1");
	write_tile_dump(master, bp, "uv_nor", key, n, rgn, regionid, nreg, reg, flag, nobn[n], ns++, R_NOR);
	write_tile_dump(master, bp, "uv_tan", key, n, rgn, regionid, nreg, reg, flag, nobt[n], ns++, R_TAN);
	if (barof) {
	  sprintf(key, "u1av");
	  write_tile_dump(master, bp, "uvav_nor", key, n, rgn, regionid, nreg, reg, flag, nobn[n], ns++, R_NOR);
	  write_tile_dump(master, bp, "uvav_tan", key, n, rgn, regionid, nreg, reg, flag, nobt[n], ns++, R_TAN);
	}
	fclose(bp);
      }
    }

    /* Write the coordinates and save mapping from centres and edges */
    /* to coordinate indices.                                        */
    n = 1;
    fprintf(op, "\nCoordinates\n");
    for (cc = 1; cc <= geom->n2_e2; cc++) {
      v = geom->w2_e2[cc];
      if (v == mask[v]) {
	fprintf(op, "%d %f %f\n", n, geom->gridx[v], geom->gridy[v]);
	mapv[v] = n;
	n++;
      }
    }
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      fprintf(op, "%d %f %f\n", n, geom->cellx[c], geom->celly[c]);
      mapc[c] = n;
      n++;
    }
    if(n - 1 != nv + nreg) hd_quit("cookie cut error.\n");

    /* Write the mesh indices                                        */
    n = 1;
    fprintf(op, "\nIndices\n");
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      fprintf(op, "%d %d %d\n",n, geom->npe[c], mapc[c]);
      for (j = 1; j <= geom->npe[c]; j++) {
	v = geom->c2v[j][c];
	fprintf(op, "%d %d ", j, mapv[v]);
	jp = (j == geom->npe[c]) ? 1 : j + 1;
	v = geom->c2v[jp][c];
	fprintf(op, "%d\n", mapv[v]);
      }
      n++;
    }

    /* Write the bathymetry                                          */
    fprintf(op, "\nBATHY    %d\n", nreg);
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      fprintf(op, "%f\n", geom->botz[c]);
    }

    fclose(op);
    i_free_1d(mapc);
    i_free_1d(mapv);
    i_free_1d(mask);
    i_free_1d(maskb);
    i_free_1d(eask);
    i_free_1d(reg);
    i_free_1d(nobc);
    i_free_1d(nobn);
    i_free_1d(nobt);
  }
  i_free_1d(rask);
  d_free_1d(regionid);
}

/* END cookie_cut()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes boundary (2way) dumpfile specifications into a df.txt file */
/*-------------------------------------------------------------------*/
void write_tile_dump(master_t *master,   /* Master data              */ 
		     FILE *op,           /* Dump file pointer        */
		     char *fname,        /* File name                */
		     char *vars,         /* Variable names           */
		     int tile,           /* Tile to dump from        */
		     int obc,            /* Tile files are used by   */
		     double *regionid,   /* Region array             */
		     int nreg,           /* # points in tile obc     */
		     int *reg,           /* Indices of tile obc      */
		     int **flag,         /* Flag array for tile obc  */
		     int nobc,           /* Number of points         */
		     int fn,             /* File number              */
		     int mode            /* Centers or edges         */
		     )
{
  geometry_t *geom = master->geom;
  char buf[MAXSTRLEN];
  int c, cc, e, j, n;

  sprintf(buf, "bdry%d-%d_%s.mpk", tile, obc, fname);
  fprintf(op, "file%d.name            %s\n", fn, buf);
  fprintf(op, "file%d.filetype        memory\n", fn);
  fprintf(op, "file%d.tinc            %5.2f seconds\n", fn, master->grid_dt);
  fprintf(op, "#file%d.sync_dt         0 seconds\n", fn);
  fprintf(op, "#file%d.filter          copy\n", fn);
  fprintf(op, "file%d.fill_rule       cascade_search\n", fn);
  fprintf(op, "file%d.bytespervalue   4\n", fn);
  fprintf(op, "file%d.vars            %s\n", fn, vars);
  if (mode & R_NOR) fprintf(op, "file%d.obc            NOR\n", fn);
  if (mode & R_TAN) fprintf(op, "file%d.obc            TAN\n", fn);
  fprintf(op, "file%d.points          %d\n", fn, nobc);

  if (mode == R_OBC) {
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      n = (int)regionid[c];
      if (n == tile) {
	if (flag[0][cc] & R_OBC)
	  fprintf(op, "%f %f\n", geom->cellx[c], geom->celly[c]);
      }
    }
  }
  if (mode & (R_NOR|R_TAN)) {
    for (cc = 0; cc < nreg; cc++) {
      c = reg[cc];
      n = (int)regionid[c];
      if (n == tile) {
	for (j = 1; j <= geom->npe[c]; j++) {
	  if (flag[j][cc] & mode) {
	    e = geom->c2e[j][c];
	    fprintf(op, "%f %f\n", geom->u1x[e], geom->u1y[e]);
	  }
	}
      }
    }
  }
  fprintf(op, "\n");
}

/* END write_tile_dump()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Constructs a run code from previous codes in output files, and    */
/* run identifiers in the parameter file.                            */
/* The code is of the form:                                          */
/* NAME|Ggrd_id|Hhyd_id|Ssed_id|Bbgc_id                              */
/* NAME is a text string, grd_id, hyd_id, sed_id and bgc_id are      */
/* floating point, e.g.                                              */
/* GBR|G2.3|H5.2|S1.0|B0.0                                           */
/* The individual ids can be in any order, but must be separated by  */
/* '|'.                                                              */
/* An id of 0.0 indicates that model component has not been invoked. */
/*-------------------------------------------------------------------*/
void create_code(parameters_t *params)
{
  char name[MAXSTRLEN], buf[MAXSTRLEN];
  char grd_c[MAXSTRLEN], hyd_c[MAXSTRLEN], sed_c[MAXSTRLEN], bgc_c[MAXSTRLEN];
  double grd_id, hyd_id, sed_id, bgc_id;
  FILE *fp;
  long t;
  int sedf = 0;

  /* Define if sediments is invoked stand alone or in ecology        */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) sedf = SEDIM;
#endif 
#if defined(HAVE_ECOLOGY_MODULE)
  sedf = 0;
  if (params->do_eco) sedf = ECOLOGY; 
  if (params->do_sed) {
    if (params->do_eco) 
      sedf = ECOLOGY;
    else
      sedf = SEDIM;
  }
#endif 

  /* Open the parameter file                                         */
  fp = fopen(params->prmname, "r+");
  strcpy(name, params->grid_name);

  /* Grid generation from AUTO with -a or -r. grd_id always set to   */
  /* 1.0 unless specified otherwise.                                 */
  if (params->runmode & AUTO) {
    strcpy(name, params->grid_name);
    if (params->runno == 0.0) {
      grd_id = 1.0;
      sprintf(grd_c, "1.0");
    } else {
      grd_id = params->runno;
      strcpy(grd_c, params->runnoc);
    }
  }
  /* netCDF INPUT file generated by -g. grd_id inherits the current  */
  /* ID_NUMBER.  If no number is specified, set to 1.0.               */
  if (params->runmode & DUMP) {
    strcpy(name, params->grid_name);
    if (params->runno == 0.0) {
      grd_id = 1.0;
      sprintf(grd_c, "1.0");
    } else {
      grd_id = params->runno;
      strcpy(grd_c, params->runnoc);
    }
  }
  /* Run using -p. Read the GRD ID from input file and set the HYD   */
  /* ID to the ID_NUMBER.                                            */
  if (params->runmode & (MANUAL | RE_ROAM)) {
    /* Read the NAME and GRD ID from INPUT_FILE                      */
    decode_id(params, params->idumpname, name, 
	      &grd_id, &hyd_id, &sed_id, &bgc_id);
    decode_idc(params, params->idumpname, name, 
	       grd_c, hyd_c, sed_c, bgc_c);

    /* Set the HYD ID                                                */
    if (params->runno == 0.0) {
      hyd_id = 1.0;
      sprintf(hyd_c, "1.0");
    } else {
      hyd_id = params->runno;
      strcpy(hyd_c, params->runnoc);
    }
  }
  /* Run using -t. Read GRD and hyd_id from the transport file, and  */
  /* set the sed_id or bgc_id to the ID_NUMBER.                      */
  if (params->runmode & TRANS) {
    /* Read the NAME, GRD ID and HYD ID from INPUT_FILE              */
    char fname[MAXSTRLEN];
    strcpy(fname, params->trans_data);
    decode_id(params, fname, name, 
	      &grd_id, &hyd_id, &sed_id, &bgc_id);
    decode_idc(params, params->idumpname, name, 
	       grd_c, hyd_c, sed_c, bgc_c);
    /* If the model component is SED, then record the runno in the   */
    /* parameter file. This is used downstream to include SED IDs    */
    /* in BGC runs.                                                  */
    if (sedf == SEDIM) {
      sed_id = params->runno;
      strcpy(sed_c, params->runnoc);
    }
    /* If the model component is BGC, then read the sediment ID from */
    /* the parameter file.                                           */
    if (sedf == ECOLOGY) {
      prm_read_char(fp, "ID_CODE", buf);
      sed_id = get_idcode(buf, "S");
      get_idcodec(buf, "S", sed_c);
      if (sed_id == 0.0)
	hd_warn("Can't find sediment transport ID in file %s\n", params->prmname);
      bgc_id = params->runno;
      strcpy(bgc_c, params->runnoc);
    } 
  }
  /* Make the id code                                                */
  sprintf(params->runcode, "%s|G%3.2f|H%3.2f|S%3.2f|B%3.2f", name,
	  grd_id, hyd_id, sed_id, bgc_id);
  sprintf(params->runcode, "%s|G%s|H%s|S%s|B%s", name,
	  grd_c, hyd_c, sed_c, bgc_c);

  if (sedf == SEDIM) {
    /*
    if (prm_skip_to_start_of_key(fp, "ID_CODE")) {
      fprintf(fp, "ID_CODE %s\n", params->runcode);
    } else {
      fseek(fp, 0, SEEK_END);
      fprintf(fp, "\nID_CODE %s\n", params->runcode);
    }
    */
  }
  /*
  time(&t);
  if (prm_skip_to_start_of_key(fp, "Last_run")) {
    fprintf(fp, "Last_run %s\n", ctime(&t));
  } else {
    fseek(fp, 0, SEEK_END);
    fprintf(fp, "Last_run %s\n", ctime(&t));
  }
  */
  fclose(fp);
}

/* END create_code()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Decodes a runcode from netCDF output into its constituent IDs.    */
/*-------------------------------------------------------------------*/
void decode_id(parameters_t *params,
	       char *ifile,
	       char *name,
	       double *grd_id,
	       double *hyd_id,
	       double *sed_id,
	       double *bgc_id
	       )
{
  int i, n, fid, ncerr;
  char rcode[MAXSTRLEN], buf[MAXSTRLEN], fname[MAXSTRLEN];
  char *tok;
  char *fields[MAXSTRLEN * MAXNUMARGS];

  n = parseline(ifile, fields, MAXNUMARGS);
  if (n == 1)
    strcpy(fname, fields[0]);
  if (endswith(fname, ".mnc")) {
    FILE *fp;
    fp = fopen(fname, "r");
    prm_read_char(fp, "file0.filename", buf);
    strcpy(fname, buf);
    fclose(fp);
  }
  /* Strip out any variable substitution                             */
  if (!(endswith(fname, ".nc"))) {
    strcpy(buf, fname);
    tok = strtok(buf, ".");
    sprintf(fname, "%s.nc", buf);
  }

  for (i = 0; i < MAXSTRLEN; i++) {
    rcode[i] = 0;
  }

  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("Can't find input file %s\n", fname);
  }
  nc_get_att_text(fid, NC_GLOBAL, "Run_code", rcode);
  strcpy(buf, rcode);
  tok = strtok(buf, "|");
  if (tok != NULL) strcpy(name, tok);
  *grd_id = get_idcode(rcode, "G");
  *hyd_id = get_idcode(rcode, "H");
  *sed_id = get_idcode(rcode, "S");
  *bgc_id = get_idcode(rcode, "B");
  nc_close(fid);
}

void decode_idc(parameters_t *params,
		char *ifile,
		char *name,
		char *grd_id,
	        char *hyd_id,
		char *sed_id,
		char *bgc_id
	       )
{
  int i, n, fid, ncerr;
  char rcode[MAXSTRLEN], buf[MAXSTRLEN], fname[MAXSTRLEN];
  char *tok;
  char *fields[MAXSTRLEN * MAXNUMARGS];

  n = parseline(ifile, fields, MAXNUMARGS);
  if (n == 1)
    strcpy(fname, fields[0]);
  if (endswith(fname, ".mnc")) {
    FILE *fp;
    fp = fopen(fname, "r");
    prm_read_char(fp, "file0.filename", buf);
    strcpy(fname, buf);
    fclose(fp);
  }
  /* Strip out any variable substitution                             */
  if (!(endswith(fname, ".nc"))) {
    strcpy(buf, fname);
    tok = strtok(buf, ".");
    sprintf(fname, "%s.nc", buf);
  }

  for (i = 0; i < MAXSTRLEN; i++) {
    rcode[i] = 0;
  }

  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("Can't find input file %s\n", fname);
  }
  nc_get_att_text(fid, NC_GLOBAL, "Run_code", rcode);
  strcpy(buf, rcode);
  tok = strtok(buf, "|");
  if (tok != NULL) strcpy(name, tok);
  get_idcodec(rcode, "G", grd_id);
  get_idcodec(rcode, "H", hyd_id);
  get_idcodec(rcode, "S", sed_id);
  get_idcodec(rcode, "B", bgc_id);
  nc_close(fid);
}

/* END decode_id()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets an id from a run code                                        */
/*-------------------------------------------------------------------*/
double get_idcode(char *code, char *id)
{
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char *tok;
  int i, nn, n;
  int j = 0;

  sprintf(key, "%c", '\0');
  strcpy(buf, code);
  tok = strtok(buf, "|");
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "G") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "H") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "S") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "B") == 0) strcpy(key, tok);
  if (strlen(key)) {
    for (n = 0; n < strlen(key) - 1; n++)
      buf[n] = key[n+1];
    buf[n] = '\0';
    return(atof(buf));
  } else {
    return(0.0);
  }
}

void get_idcodec(char *code, char *id, char *ret)
{
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char *tok;
  int i, nn, n;
  int j = 0;

  sprintf(key, "%c", '\0');
  strcpy(buf, code);
  tok = strtok(buf, "|");
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "G") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "H") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "S") == 0) strcpy(key, tok);
  tok = strtok(NULL, "|");
  if (tok != NULL && strcmp(id, "B") == 0) strcpy(key, tok);
  if (strlen(key)) {
    for (n = 0; n < strlen(key) - 1; n++)
      buf[n] = key[n+1];
    buf[n] = '\0';
  } else {
    strcpy(buf, "0.0");
  }
  strcpy(ret, buf);
}

/* END get_idcode()                                                  */
/*-------------------------------------------------------------------*/


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
  if (m & U1BDRY)
    return ("u1");
  else if (m & U2BDRY)
    return ("u2");
  else if (m & (U1BDRY|U2BDRY))
    return ("velocity");
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
    /*
  case 131072:
    return ("TIDALM");
    */
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
  case ORDER4:
    return ("ORDER4");
  case QUICKEST:
    return ("QUICKEST");
  case QUICKEST|HIORDER:
    return ("QUICKEST|US");
  case VANLEER:
    return ("VANLEER");
  case VANLEER|HIORDER:
    return ("VANLEER|US");
  case ORDER2_UW:
    return ("ORDER2_UPWIND");
  case LAGRANGE:
    return ("LAGRANGE");
  case FFSL:
      return ("FFSL");
  case WIMPLICIT:
    return ("WIMPLICIT");
  case RINGLER:
    return ("VECINVAR");
  case RINGLER|WTOP_O2|WIMPLICIT:
    return ("VECINVAR WTOP_O2 WIMPLICIT");
  case RINGLER|WTOP_O2|WIMPLICIT|PV_ENEUT:
    return ("VECINVAR WTOP_O2 WIMPLICIT NEUTRAL");
  case RINGLER|WTOP_O2|WIMPLICIT|PV_ENEUT|PV_APVM:
    return ("VECINVAR WTOP_O2 WIMPLICIT NEUTRAL APVM");
  case RINGLER|WTOP_O2|WIMPLICIT|PV_ENEUT|PV_LUST:
    return ("VECINVAR WTOP_O2 WIMPLICIT NEUTRAL LUST");
  case RINGLER|WTOP_O2|WIMPLICIT|PV_ENEUT|PV_CLUST:
    return ("VECINVAR WTOP_O2 WIMPLICIT NEUTRAL CLUST");
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
  params->nbu = params->nobc;

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

  /* Check for tokens without a value */
  for (n = 0; n < nf; n++) {
    if (strcmp(fields[n], tag) == 0)
      return(2);
  }
  /* Check for tokens with a value separated by ':' */
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

int has_neighbour(int c, int npe, int **neic)
{
  int j;
  for (j = 1; j <= npe; j++) {
    if (neic[j][c] == 0)
      return(j);
  }
  return(0);
}

int is_edge(parameters_t *params, int cc, double *bathy, int **neic) {
  int cn, j;
  int found = 0;

  for (j = 1; j <= params->npe2[cc]; j++) {
    cn = neic[j][cc];
    if (!cn)
      return(1);
  }
  return(0);
}

int has_outside_neighbour(parameters_t *params, int cc, double *bathy, int **neic) {
  int cn, j;
  int found = 0;

  for (j = 1; j <= params->npe2[cc]; j++) {
    cn = neic[j][cc];
    if (is_outside(params, bathy[cn]))
      return(cn);
  }
  return(0);
}

int has_wet_neighbour(parameters_t *params, int cc, double *bathy, int **neic) {
  int cn, j;
  int found = 0;

  for (j = 1; j <= params->npe2[cc]; j++) {
    cn = neic[j][cc];
    if (is_wet(params, bathy[cn]))
      return(cn);
  }
  return(0);
}


void reset_bathy(parameters_t *params, int io, int jo, double val)
{
  int i, j, n = 0;
  int is = 0;
  int js = 0;
  int ie = params->nce1;
  int je = params->nce2;
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
  int n, nobc = params->nobc;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  if (params->riverflow > 0) return 1;

  sprintf(keyword, "NBOUNDARIES");
  if (nobc || prm_read_int(fp, keyword, &nobc)) {
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
      /*
      if (params->open[n] != NULL) {
	if (strlen(params->open[n]->bflow) &&
	    (strcmp(params->open[n]->cusname_u1, "u1flowbdry") == 0 ||
	     strcmp(params->open[n]->cusname_u2, "u2flowbdry") == 0))
	  return(1);
      }
      */
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



/** Get z values for a grid
 *
 * @param tfs datastructure that holds topo file data
 * @param gx pointer to the 2D array of the centre x/long values
 *           gx.size = [ns2 + 1]
 * @param gy pointer to the 2D array of the centre y/lat values
 *           gy.size = [ns2 + 1]
 * @param ns2 number of cellular centres
 * @param hint the hint to be used for calculating the z values
 * @return pointer to 2D array to house the z values aka topo
 *         topo.size = [nce2] x [nce1]
 */
double* topo_get_z_for_grid_us(topo_files_t *tfs, double *gx, double *gy,
   int ns2, topo_hint_t hint) {
   int got_topo_ok = 1;
   int i, j;
   double x_ce, y_ce;
   double *topo, *x, *y, *z;

   x = d_alloc_1d(4);
   y = d_alloc_1d(4);
   z = d_alloc_1d(4);
   topo = d_alloc_1d(ns2+1);
   for(i=1;i<=ns2;++i) {
     /* get the centroid of the cell, a simple approach */
     x_ce = gx[i];
     y_ce = gy[i];

     /* get the data based on the hint */
     if(hint == NEAREST) {
       /* now get the value nearest the centroid */
       topo[i] = topo_get_z_at_point(tfs, x_ce, y_ce, hint);
     } else if(hint == AVERAGE || hint == AUTO || hint == INVERSE) {
       /* build the 4 point poly */
       /*
	 x[0] = gx[j][i];
	 y[0] = gy[j][i];
	 x[1] = gx[j+1][i];
	 y[1] = gy[j+1][i];
	 x[2] = gx[j][i+1];
	 y[2] = gy[j][i+1];
	 x[3] = gx[j+1][i+1];
	 y[3] = gy[j+1][i+1];
	 topo[j][i] = topo_get_z_under_area(tfs, x, y, 4, hint);
       */
     } else if(hint == LINEAR) {
       /* assume that interpolation is to be based on the corners */
       /* NOTE: only bilinear for now */
       /* NOTE: using NEAREST for now, may want others later */
       /*
	 x[0] = gx[j][i];
	 y[0] = gy[j][i];
	 z[0] = topo_get_z_at_point(tfs, x[0], y[0], NEAREST);
	 x[1] = gx[j+1][i];
	 y[1] = gy[j+1][i];
	 z[1] = topo_get_z_at_point(tfs, x[1], y[1], NEAREST);
	 x[2] = gx[j][i+1];
	 y[2] = gy[j][i+1];
	 z[2] = topo_get_z_at_point(tfs, x[2], y[2], NEAREST);
	 x[3] = gx[j+1][i+1];
	 y[3] = gy[j+1][i+1];
	 z[3] = topo_get_z_at_point(tfs, x[3], y[3], NEAREST);
	 topo[j][i] = bilinear(x_ce, y_ce, x, y, z);
       */
     } else {
       /* else just pass the cell center on to "at_point" for now */
       topo[i] = topo_get_z_at_point(tfs, x_ce, y_ce, hint);
     }
     /* if any reads resulted in a NaN the read was not OK */
     if(topo[i] == NaN)
       got_topo_ok = 0;
   }

  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(z);

  return topo;
}

void not_implemented(char *text)
{
  hd_warn("The function %s is not yet implemented in COMPAS.\n", text);
}

void not_included(char *text)
{
  hd_warn("Functionality %s is not included in COMPAS.\n", text);
}

