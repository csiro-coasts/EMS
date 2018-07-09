/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/lib/waves/waves.c
 *  
 *  Description:
 *  Wave initialisation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: waves.c 5849 2018-06-29 05:03:57Z riz008 $
 *
 */

#include "waves.h"
#include <ctype.h>

#define NOTVALID -9999
#define U1BDRY 0x040   /* u1 point open boundary                         */
#define L_EDGE 0x080   /* u1 bdry point on left side of interior cell    */
#define R_EDGE 0x100   /* u1 bdry point on right side of interior cell   */
#define U2BDRY 0x200   /* u2 point open boundary                         */
#define B_EDGE 0x400   /* u2 bdry point on back side of interior cell    */
#define F_EDGE 0x800   /* u2 bdry point on forward side of interior cell */

double g = 9.81;

static double get_fetch(double *fetch, double dir);
static void wind_wave_toba(wave_t *wave, double *f);
static void wind_wave_usac(wave_t *wave, double *f);
static double fwc94(double cmu, double cukw);
static void stress2wind(double *wx, double *wy);
static void madsen94o(int i, int j, double ubr, double wr, double ucr,
                      double zr, double phiwc, double zo,
                      double *pustrc, double *pustrwm,
                      double *pustrr, double *pfwc, double *pzoa);

/* Version information */
int get_waves_major_vers(void)
{
  return(WAVES_MAJOR_VERSION);
}

int get_waves_minor_vers(void)
{
  return(WAVES_MINOR_VERSION);
}

int get_waves_patch_vers(void)
{
  return(WAVES_PATCH_VERSION);
}

/*-------------------------------------------------------------------*/
/* Builds and initialises the interface                              */
/*-------------------------------------------------------------------*/
wave_t* wave_build(void* model, FILE *fp)
{
  wave_t* wave = wave_create();

  wave_init(model, wave, fp);
  return wave;
}
/* END wave_build()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocates memory for the interface structure                      */
/*-------------------------------------------------------------------*/
wave_t* wave_create() {
  wave_t *wave = (wave_t *)malloc(sizeof(wave_t));

  wave->do_waves = 0;
  wave->do_amp = WNONE;
  wave->do_per = WNONE;
  wave->do_dir = WNONE;
  wave->do_ub = WNONE;
  wave->do_wif = WNONE;
  wave->do_stokes = WNONE;
  wave->do_rs = 0;
  wave->do_Cd = 0;
  wave->do_fetch = 0;
  wave->windwave = TOBA;
  wave->cols = -1;
  wave->nz = -1;
  wave->sednz = -1;
  wave->ntr = 0;
  wave->ntrS = 0;
  wave->atr = 0;
  wave->atrS = 0;
  wave->nsed = 0;
  wave->dt = 0.0;
  wave->dt_ratio = 0;
  wave->topk_wc = -1;
  wave->botk_wc = -1;
  wave->topk_sed = -1;
  wave->botk_sed = -1;
  wave->top = 0;
  wave->bot = 0;
  wave->bot_l = 0;
  wave->dz_wc = NULL;
  wave->Kz = NULL;
  wave->Vz = NULL;
  wave->u1 = NULL;
  wave->u2 = NULL;
  wave->w = NULL;
  wave->u1bot = 0;
  wave->u2bot = 0;
  wave->z0 = 0;
  wave->quad_bfc;
  strcpy(wave->amp_name, "wave_amp");
  wave->amp = NULL;
  wave->ampid = -1;
  strcpy(wave->period_name, "wave_period");
  wave->period = NULL;
  wave->perid = -1;
  strcpy(wave->dir_name, "wave_dir");
  wave->dir = NULL;
  wave->dirid = -1;
  strcpy(wave->ub_name, "wave_ub");
  wave->ub = NULL;
  wave->ubid = -1;
  strcpy(wave->Sxy_name, "wave_Sxy");
  wave->Sxy = NULL;
  wave->Sxyid = -1;
  strcpy(wave->Syx_name, "wave_Syx");
  wave->Syx = NULL;
  wave->Syxid = -1;
  strcpy(wave->Fx_name, "wave_Fx");
  wave->Fx = NULL;
  wave->Fxid = -1;
  strcpy(wave->Fy_name, "wave_Fy");
  wave->Fy = NULL;
  wave->Fyid = -1;
  strcpy(wave->ste1_name, "stx");
  wave->ste1 = NULL;
  wave->ste1id = -1;
  strcpy(wave->ste2_name, "sty");
  wave->ste2 = NULL;
  wave->ste2id = -1;


  strcpy(wave->ustrcw_name, "ustrcw");
  wave->ustrcw = NULL;
  wave->ustid = -1;
  strcpy(wave->Cd_name, "wave_Cd");
  wave->Cd = NULL;
  wave->Cdid = -1;


  wave->eta = 0.0;
  wave->area = 0.0;
  wave->h1au2 = 0.0;
  wave->h2au1 = 0.0;
  wave->fetch = NULL;
  wave->thetau1 = NULL;
  wave->thetau2 = NULL;
  wave->sinthcell = NULL;
  wave->costhcell = NULL;
  wave->trname_3d = NULL;
  wave->trname_2d = NULL;
  wave->trname_sed = NULL;
  wave->tmap_3d = NULL;
  wave->tmap_2d = NULL;
  wave->tmap_sed = NULL;
  wave->brsm = NULL;
  wave->depth_wc = -1;
  wave->eta = -1;
  wave->model = NULL;

  return wave;
}
/* END wave_create()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialises the tracer statistics structure using information     */
/* from the host model.                                              */
/*-------------------------------------------------------------------*/
void wave_init(void* model, wave_t *wave, FILE *fp) {
  int n, nn, m, i, j, k, ntr_tmp;
  int ntr, atr;           /* Number of 3D tracers in host model      */
  int ntrS, atrS;         /* Number of 2D tracers in host model      */
  char **trname_3d;       /* 3D tracer names                         */
  char **trname_2d;       /* 2D tracer names                         */
  char **trname_sed;      /* Sediment tracer names                   */
  int size;               /* Size of the grid                        */
  double dt;
  char buf[MAXSTRLEN];

  /*-----------------------------------------------------------------*/
  /* Set up the grid                                                 */
  /* Get the wave timestep                                           */
  wave->dt = w_get_dt(model);
  dt = i_get_model_timestep(model);
  wave->dt_ratio = ceil(wave->dt / dt);
  if (wave->dt_ratio == 0) {
    wave->dt_ratio = 1;
    emstag(LWARN,"waves:init","Wave time ratio = 0 : setting to 1\n");
  }
  emstag(LTRACE,"waves:init","Wave time ratio = %d \n", wave->dt_ratio);

  /* Get the number of water column layers                           */
  wave->nz = i_get_num_wclayers(model);

  /* Allocate memory for water column layer thickness                */
  wave->dz_wc = d_alloc_1d(wave->nz);

  /* Get the number of sediment layers                               */
  wave->sednz = i_get_num_sedlayers(model);

  /* Get the number of water columns                                 */
  wave->cols = i_get_num_columns(model);

  /* Get the timestep                                                */
  wave->dt = i_get_model_timestep(model);

  /* Get the model start time                                        */
  wave->time = i_get_model_time(model);

  /* Assign a reference to the parent model                          */
  wave->model = model;

  /* Get the quadratic friction coefficient                          */
  wave->quad_bfc = w_get_quad_bfc(model);

  /* Set the method of computing the wave variables                  */
  wave->do_amp = wave->do_per = wave->do_dir = wave->do_ub = WCOMP;
  w_check_wave_data(&n, &m, &i, &j, &k, &nn);
  if (n) wave->do_amp = WFILE;
  if (m) wave->do_per = WFILE;
  if (i) wave->do_dir = WFILE;
  if (j) wave->do_ub = WFILE;
  if (k) wave->do_wif = WFILE;
  if (nn) wave->do_stokes = WFILE;

  /* Check that wave variables exist */
  if(!w_check_orbital_file(model)) {
    emstag(LFATAL,"waves:init"," Wave input file is undefined. Use ORBITAL_VEL in the parameter file.\n");
    exit(0);
  }
  if (!w_check_wave_period(model)) {
    emstag(LFATAL,"waves:init"," Wave period variable, '%s', undefined.\n", wave->period_name );
    exit(0);
  }
  if (!w_check_wave_amp(model)) {
    emstag(LFATAL,"waves:init"," Wave amplitude variable, '%s', undefined.\n", wave->amp_name);
    exit(0);
  }
  if (!w_check_wave_dir(model)) {
    emstag(LFATAL,"waves:init"," Wave direction variable, '%s', undefined.\n", wave->dir_name);
    exit(0);
  }
  if (!w_check_wave_ub(model)) {
    emstag(LFATAL,"waves:init"," Wave orbital velocity variable, '%s', undefined.\n", wave->ub_name);
    exit(0);
  }
  if (wave->do_wif & WFILE) {
    if (!w_check_wave_Fx(model)) {
      emstag(LFATAL,"waves:init"," Wave-induced force e1 variable, '%s', undefined.\n", wave->Fx_name);
      exit(0);
    }
    if (!w_check_wave_Fy(model)) {
      emstag(LFATAL,"waves:init"," Wave-induced force e2 variable, '%s', undefined.\n", wave->Fy_name);
      exit(0);
    }
  }
  if (!w_check_wind(model)) {
    emstag(LFATAL,"waves:init"," Wind variables undefined.\n");
    exit(0);
  }

  if (wave->do_stokes & WFILE) {
    if (!w_check_wave_ste1(model)) {
      emstag(LFATAL,"waves:init"," Stokes velocity x variable, '%s', undefined.\n", wave->ste1_name);
      exit(0);
    }
    if (!w_check_wave_ste2(model)) {
      emstag(LFATAL,"waves:init"," Stokes velocity y variable, '%s', undefined.\n", wave->ste2_name);
      exit(0);
    }
  }

  /* Get the 3d tracers names                                        */
  ntr = i_get_num_tracers_3d(model, &atr);
  trname_3d = malloc(ntr * sizeof(char *));
  for (n = 0; n < ntr; n++)
    trname_3d[n] = malloc(MAXSTRLEN * sizeof(char));
  i_get_names_tracers_3d(model, trname_3d);

  /* Get the 2d tracers names                                        */
  ntrS = i_get_num_tracers_2d(model, &atrS);
  trname_2d = malloc(ntrS * sizeof(char *));
  for (n = 0; n < ntrS; n++)
    trname_2d[n] = malloc(MAXSTRLEN * sizeof(char));
  i_get_names_tracers_2d(model, trname_2d);

  /* Get the number of 2D tracers required by the wave model         */
  wave->ntrS = 0;
  for (n = 0; n < ntrS; n++){
    if (strcmp(wave->period_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->amp_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->dir_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->ub_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->Sxy_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->Syx_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->Fx_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->Fy_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->ste1_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->ste2_name, trname_2d[n]) == 0)
      wave->ntrS++;

    if (strcmp(wave->ustrcw_name, trname_2d[n]) == 0)
      wave->ntrS++;
    if (strcmp(wave->Cd_name, trname_2d[n]) == 0)
      wave->ntrS++;

  }

  /* Set up the derived wave variables                               */
  /* Get the names of the wave variables                             */
  if(wave->ntrS) {
    wave->trname_2d = malloc(wave->ntrS * sizeof(char *));
    for (n = 0; n < wave->ntrS; n++) {
      wave->trname_2d[n] = malloc(MAXSTRLEN * sizeof(char));
    }
  }
  /* Allocate memory for wave variables & get the wave map function  */
  wave->tr_in = (double**)p_alloc_1d(wave->ntrS);

  /* Assign the wave pointers                                       */
  m = 0;
  for (n = 0; n < ntrS; n++){
    if (strcmp(wave->period_name, trname_2d[n]) == 0) {
      wave->perid = m;
      strcpy(wave->trname_2d[m], wave->period_name);
      m++;
    }
    if (strcmp(wave->amp_name, trname_2d[n]) == 0) {
      wave->ampid = m;
      strcpy(wave->trname_2d[m], wave->amp_name);
      m++;
    }
    if (strcmp(wave->dir_name, trname_2d[n]) == 0) {
      wave->dirid = m;
      strcpy(wave->trname_2d[m], wave->dir_name);
      m++;
    }
    if (strcmp(wave->ub_name, trname_2d[n]) == 0) {
      wave->ubid = m;
      strcpy(wave->trname_2d[m], wave->ub_name);
      m++;
    }
    if (strcmp(wave->Sxy_name, trname_2d[n]) == 0) {
      wave->Sxyid = m;
      strcpy(wave->trname_2d[m], wave->Sxy_name);
      m++;
    }
    if (strcmp(wave->Syx_name, trname_2d[n]) == 0) {
      wave->Syxid = m;
      strcpy(wave->trname_2d[m], wave->Syx_name);
      m++;
    }
    if (strcmp(wave->Fx_name, trname_2d[n]) == 0) {
      wave->Fxid = m;
      strcpy(wave->trname_2d[m], wave->Fx_name);
      wave->do_wif = 1;
      m++;
    }
    if (strcmp(wave->Fy_name, trname_2d[n]) == 0) {
      wave->Fyid = m;
      strcpy(wave->trname_2d[m], wave->Fy_name);
      m++;
    }
    if (strcmp(wave->ste1_name, trname_2d[n]) == 0) {
      wave->ste1id = m;
      strcpy(wave->trname_2d[m], wave->ste1_name);
      m++;
    }
    if (strcmp(wave->ste2_name, trname_2d[n]) == 0) {
      wave->ste2id= m;
      strcpy(wave->trname_2d[m], wave->ste2_name);
      m++;
    }

    if (strcmp(wave->ustrcw_name, trname_2d[n]) == 0) {
      wave->ustid = m;
      strcpy(wave->trname_2d[m], wave->ustrcw_name);
      m++;
    }
    if (strcmp(wave->Cd_name, trname_2d[n]) == 0) {
      wave->Cdid = m;
      strcpy(wave->trname_2d[m], wave->Cd_name);
      wave->do_Cd = 1;
      m++;
    }


  }

  if (wave->Sxyid >= 0 && wave->Syxid >= 0)
     wave->do_rs = 1;
  if (wave->ste1id >= 0 && wave->ste2id >= 0)
    wave->do_stokes = 1;

  /* Get the tracer map and check if all tracers have been found    */
  wave->tmap_2d = i_get_tmap_2d(model, wave->ntrS, wave->trname_2d);
  for(n = 0; n < wave->ntrS; n++) {
    if(wave->tmap_2d[n] < 0 ) {	
      emstag(LPANIC,"waves:init","Attempting to reference a 2d tracer which doesn't exist in host model [%s] exiting...",wave->trname_2d[n]);
      exit(1);
    } else
      emstag(LTRACE,"waves:init","Wave tracer %d = %s.\n", n, wave->trname_2d[n]);
  }

  /* Get the boundary radiation stress mask                         */
  size = w_get_winsize(model);
  wave->brsm = i_alloc_1d(size);
  memset(wave->brsm, 0, n * sizeof(int));
  w_get_brsm(model, wave->brsm);

  /* Get the grid angles                                            */
  wave->thetau1 = d_alloc_1d(size);
  wave->thetau2 = d_alloc_1d(size);
  wave->sinthcell = d_alloc_1d(size);
  wave->costhcell = d_alloc_1d(size);
  for (i = 1; i <= wave->cols; i++) {
    wave->thetau1[i] = w_get_thetau1(model, i);
    wave->thetau2[i] = w_get_thetau2(model, i);
    wave->sinthcell[i] = w_get_sinthcell(model, i);
    wave->sinthcell[i] = w_get_costhcell(model, i);
  }

  /* Get the fetch                                                  */
  if (w_check_fetch(model)) {
    wave->do_fetch = 1;
    wave->fetch = d_alloc_2d(8, size);
    for (i = 1; i <= wave->cols; i++) {
      for (n = 0; n < 8; n++)
	wave->fetch[i][n] = w_get_fetch(model, i, n);
    }
  }

  /* Get the wind wave method                                       */
  if (!(wave->do_per & WFILE)) {
    if (prm_read_char(fp, "WIND_WAVE", buf)) {
      if (strcmp(buf, "TOBA") == 0) {
	wave->windwave = TOBA;
	emstag(LTRACE,"waves:init","Using Toba (1978) for wind waves.\n");
      }
      if (strcmp(buf, "USAC") == 0) {
	wave->windwave = USAC;
	emstag(LTRACE,"waves:init","Using US Army Core of Engineers for wind waves.\n");
      }
    } else {
      wave->windwave = TOBA;
      emstag(LWARN,"waves:init","Can't file wave period file input. Using Toba (1978) for wind waves.\n");
    }
    wave->do_per = WWIND;
    if (wave->do_amp & WFILE)
      emstag(LWARN,"waves:init","Overwriting wave amplitude file input with wind wave amplitude.\n");
    wave->do_amp = WWIND;
    if (wave->do_dir & WFILE)
      emstag(LWARN,"waves:init","Overwriting wave direction file input with wind wave direction.\n");
    wave->do_dir = WWIND;
  }

  /* Free memory                                                    */
  free((char **)trname_2d);
  free((char **)trname_3d);

  wave->do_waves = 1;
  emstag(LTRACE,"waves:init","Waves initialised OK.\n");
}


/*-------------------------------------------------------------------*/
/* Does the interface work.                                          */
/*-------------------------------------------------------------------*/
void wave_step(void* model, wave_t *wave, int c) {
  int cc, n, i;
  double wamp;
  double dir;
  double wavenumb;
  double ub;
  double ustrcw;
  double wx, wy;
  double thetau1;     /* angle between e1 & x axis at u1 point       */ 
  double thetau2;     /* angle between e2 & y axis at u2 point       */ 

  /* Return if nothing to do for this step                           */
  if(!wave->do_waves)
    return;

  /*-----------------------------------------------------------------*/
  /* Get the model attributes                                        */
  /* Get the interface column index                                  */
  i = i_get_interface_counter(model, c);
  wave->c = c;
  wave->cc = i + 1;
  i_get_ijk(model, c, wave->ij);

  /* Get the top index of the water column                           */
  wave->topk_wc = i_get_topk_wc(model, c);

  /* Get the bottom index of the water column                        */
  wave->botk_wc = i_get_botk_wc(model, c);

  /* Get the height of the top layer                                 */
  wave->top = i_get_topz_wc(model, c);

  /* Get the depth of the sea bed                                    */
  wave->bot = i_get_botz_wc(model, c);

  /* Get the timestep for this step                                  */
  wave->dt = i_get_model_timestep(model);

  /* Get the model time at this step                                 */
  wave->time = i_get_model_time(model);

  /* Get the water column layer thickness                            */
  i_get_dz_wc(model, c, wave->dz_wc);

  /* Get the water depth                                             */
  wave->depth_wc = i_get_depth(model, c);

  /* Get the cell area                                               */
  wave->area = i_get_cellarea_w(model, c);

  /* Get the borrom roughness parameter                              */
  wave->z0 = w_get_z0(model, c);

  /* Get the surface elevation                                       */
  wave->eta = i_get_eta(model, c);

  /* Get the bottom velocities                                       */
  w_get_bot_vel(model, wave->sinthcell[wave->cc], wave->costhcell[wave->cc],
		&wave->u1bot, &wave->u2bot, & wave->bot_l, c);

  /* Get the wind components                                         */
  wave->wx = w_get_wave_wind1(model, c);
  wave->wy = w_get_wave_wind2(model, c);

  /* Read in the 2D tracers                                          */
  i_get_tracer_2d(model, c, wave->ntrS, wave->tmap_2d, wave->tr_in);

  /* Assign the pointers                                             */
  wave->amp = wave->tr_in[wave->ampid];
  wave->period = wave->tr_in[wave->perid];
  wave->dir = wave->tr_in[wave->dirid];
  wave->ub = wave->tr_in[wave->ubid];
  wave->ustrcw = wave->tr_in[wave->ustid];
  if (wave->do_rs) {
    wave->Sxy = wave->tr_in[wave->Sxyid];
    wave->Syx = wave->tr_in[wave->Syxid];
  }
  if (wave->do_wif) {
    wave->Fx = wave->tr_in[wave->Fxid];
    wave->Fy = wave->tr_in[wave->Fyid];
  }
  if (wave->do_stokes) {
    wave->ste1 = wave->tr_in[wave->ste1id];
    wave->ste2 = wave->tr_in[wave->ste2id];
  }
  if (wave->do_Cd)
    wave->Cd = wave->tr_in[wave->Cdid];

  /* Estimate wave period from wind if not provided                  */
  if (wave->do_per & WWIND && wave->do_fetch) {
    if (wave->windwave & TOBA)
      wind_wave_toba(wave, wave->fetch[wave->cc]);
    if (wave->windwave & USAC)
      wind_wave_usac(wave, wave->fetch[wave->cc]);
  }

  /* Estimate wave direction from wind if not provided               */
  if (wave->do_dir & WCOMP)
    *wave->dir = wavedir_estim(wave->wy, wave->wx);

  /* Estimate wave amplitude if not provided                         */
  if (wave->do_amp & WCOMP) {
    wavenumb = wavenumb_estim(*wave->period, wave->depth_wc);
    *wave->amp = wamp_estim(wave->wx, wave->wy, wavenumb, *wave->period);
  }

  /* Estimate wave orbital velocity if not provided                  */
  if (wave->do_ub & (WCOMP|WWIND)) {
    wavenumb = wavenumb_estim(*wave->period, wave->depth_wc);
    *wave->ub = waveub_estim(*wave->amp, *wave->period, wavenumb, 
			     wave->depth_wc);
  }

  /* Get the wave current friction                                   */
  if (wave->do_Cd) {
     ustrcw_estim(wave, c);
  }

  /* Estimate the radiation stresses if not provided                 */
  if (wave->do_rs) {
    radiation_stress(wave, c);
  }
}

/*-------------------------------------------------------------------*/
/* Writes wave info into the setup.txt file                          */
/*-------------------------------------------------------------------*/
void wave_run_setup(FILE *fp, wave_t *wave)
{
  fprintf(fp, "\n");
  fprintf(fp, "Waves library invoked\n");
  
  if (wave->do_per & WWIND && wave->do_fetch) {
    fprintf(fp, "  wave amplitude, period & direction computed from winds\n");
  }
  if (wave->do_dir & WCOMP)
    fprintf(fp, "  wave direction computed\n");

  if (wave->do_amp & WCOMP)
    fprintf(fp, "  wave amplitude computed\n");

  if (wave->do_ub & (WCOMP|WWIND))
    fprintf(fp, "  wave orbital velocity computed\n");

  if (wave->do_Cd)
    fprintf(fp, "  wave current friction computed\n");

  if (wave->do_rs)
    fprintf(fp, "  wave radiation stresses computed\n");

  if (wave->do_wif & WFILE)
    fprintf(fp, "  wave induced forces read from file\n");
  
  fprintf(fp, "\n");
}

/*-------------------------------------------------------------------*/
/* Computes radiations stresses in deep water and at the coast.      */
/* Uses the methodology of Bye (1977) The flow series of Thallasso - */
/* models. The FLOWM model, supplement to the FLOWC model. Selected  */
/* topics in atmospheric and marine sciences, No 6, Flinders Uni. of */
/* South Australia.                                                  */
/*-------------------------------------------------------------------*/
void radiation_stress(wave_t *wave,      /* Wave data                */
		      int c              /* Sparse location          */
		      )
{
  int cc = wave->cc;        /* Array index                           */
  double dir = *wave->dir;  /* Wave direction (deg T)                */
  double wamp = *wave->amp; /* Wave amplitude (m)                    */
  double wx, wy;            /* Wave components in Cartesian coords   */
  double wcx, wcy;          /* Wave components in local coords       */
  double dir1;              /* Wave direction relative to e1 axis    */
  double dir2;              /* Wave direction relative to e2 axis    */
  double ldir;              /* Wave direction in local coordinates   */
  double Sxy, Syx;          /* Radiation stresses at the coast       */
  double S;                 /* Dummy                                 */
  double costhu1;
  double sinthu1;
  double costhu2;
  double sinthu2;
  double d1,d2;

  /* Get the grid angles                                             */
  costhu1 = cos(wave->thetau1[cc]);
  sinthu1 = sin(wave->thetau1[cc]);
  costhu2 = cos(wave->thetau2[cc]);
  sinthu2 = sin(wave->thetau2[cc]);

  /* Get the radiation stresses in deep water                        */
  /* Note; dir is the direction (deg T) the wave train is heading    */
  d1 = dir * PI / 180.;
  *wave->Syx = 0.25 * g * wamp * wamp * sin(d1) * cos(d1);
  *wave->Sxy = *wave->Syx;

  /* Get the wave components in the cartesian system                 */
  wavecomp_estim(wamp, dir, &wx, &wy);

  /* Rotate to the curvilinear system (local coordinates)            */
  wcx = wx * costhu1 + wy * sinthu1;
  wcy = -wx * sinthu2 + wy * costhu2;

  /* Get the direction relative to local coordinates                 */
  ldir = d1 = (wamp > 0.) ? acos(wcy / wamp) : 0.0;
  ldir = d2 = (wcx < 0.) ? -ldir * 180. / PI : ldir * 180. / PI;
  ldir = (ldir < 0.) ? ldir + 360. : ldir;

  /* Get the direction relative to e1 and e2 local axes              */
  /* |dir1| < 90 and |dir2| < 90                                     */
  dir1 = (ldir > 180.) ? ldir - 270. : 90. - ldir;
  dir2 = (ldir < 90.) ? ldir : (ldir > 270) ? ldir - 360. : 180. - ldir;
  /*if (ldir <= 180.) dir1 *= -1.0;*/
  /*if (ldir >= 270. || ldir <= 90.) dir2 *= -1.0;*/

  /* Get the radiation stresses at the coast                         */
  Syx = g * wamp * wamp  * 0.1 * dir1 / 90.0;
  Sxy = g * wamp * wamp  * 0.1 * dir2 / 90.0;

  /* Set the radiation stresses across lateral boundaries            */
  /* Northern and southern boundaries                                */
  if (wave->brsm[c] & U2BDRY) {
    if (wave->brsm[c] & F_EDGE) {
      Sxy = (ldir <= 90. || ldir >= 270.) ? Sxy : 0.;
    } else if (wave->brsm[c] & B_EDGE) {
      Sxy = (ldir >= 90. && ldir <= 270.) ? Sxy : 0.;
    }
    *wave->Sxy = Sxy;
  }
  /* Western and eastern boundaries                                  */
  if (wave->brsm[c] & U1BDRY) {
    if (wave->brsm[c] & R_EDGE) {
      Syx = (ldir >= 0. && ldir <= 180.) ? Syx : 0.;
    } else if (wave->brsm[c] & L_EDGE) {
      Syx = (ldir >= 180. || ldir <= 0.) ? Syx : 0.;
    }
    *wave->Syx = Syx;
  }
}


/* The minimum reference level set to 1 m so that in a typical application 
   it is above the wave bbl */
/*-------------------------------------------------------------------*/
/* Compute the friction due to waves at the sea bed                  */
/*-------------------------------------------------------------------*/
void ustrcw_estim(wave_t *wave, int c)
{
  double uval = wave->u1bot; /* Cell centre e1 velocity magnitude    */
  double vval = wave->u2bot; /* Cell centre e2 velocity magnitude    */
  double udir;               /* Cell centre velocity direction       */
  double u;                  /* Cell centre speed                    */
  double period = *wave->period;
  double dir = *wave->dir;
  double ub = *wave->ub;
  double z0 = wave->z0; //???????
  double ustrc, ustrwm, ustrr, fwc, zoa, zr, zc, v;

  /* Compute the reference height                                    */
  if (wave->topk_wc > wave->botk_wc)
    zc = (wave->bot_l - wave->bot) / 2.;
  else
    zc = (wave->top - wave->bot) / 2.;


 /* NMY If the reference height is less than 1m use log-profile to 
    interpolate velocities to 1m hight, so that in a typical application 
    zr exceeds the thickness of the wave bbl */
  zr=zc;
  if (zc<1.) zr=1.;
  if (zc < z0) {
    emstag(LPANIC,"waves:ustrcw_estim",
	   "Cell centre below the bottom roughness height\n");
    return;
  }
  // Reference velocity 
  uval = uval * log(zr / z0) / log(zc / z0);
  vval = vval * log(zr / z0) / log(zc / z0);
  /* end NMY */


  if ((uval == 0.0) && (vval == 0.0)) {
    u = 0.0;
    udir = dir;
  } else {
    u = sqrt(uval * uval + vval * vval);
    udir = fmod(450.0 - atan2(vval, uval) * 180.0 / PI, 180.0);
  }

  /* Compute the apparent Z0.                                        */
  madsen94o(wave->ij[0], wave->ij[1], ub, 2 * M_PI / period, u, zr, 
	    dir - udir, wave->z0, &ustrc, &ustrwm, &ustrr, &fwc, &zoa);

  /* Store bottom friction velocity FIX - Is this the right one???   */
  if (!isnan(ustrr))
    *wave->ustrcw = ustrr;
  else
    *wave->ustrcw = 0.0;

  v = log(zc / zoa) / VON_KAR;

  if (wave->do_Cd) {
    if(zc>zoa) {
      v = log(zc / zoa) / VON_KAR;
      *wave->Cd = max(wave->quad_bfc, 1.0 / (v * v));
    } else
      *wave->Cd = wave->quad_bfc;
  }
}



/*-------------------------------------------------------------------*/
/* Compute the friction due to waves at the sea bed                  */
/*-------------------------------------------------------------------*/
void ustrcw_estim_original(wave_t *wave, int c)
{
  double uval = wave->u1bot; /* Cell centre e1 velocity magnitude    */
  double vval = wave->u2bot; /* Cell centre e2 velocity magnitude    */
  double udir;               /* Cell centre velocity direction       */
  double u;                  /* Cell centre speed                    */
  double period = *wave->period;
  double dir = *wave->dir;
  double ub = *wave->ub;
  double ustrc, ustrwm, ustrr, fwc, zoa, zr, v;

  if ((uval == 0.0) && (vval == 0.0)) {
    u = 0.0;
    udir = dir;
  } else {
    u = sqrt(uval * uval + vval * vval);
    udir = fmod(450.0 - atan2(vval, uval) * 180.0 / PI, 180.0);
  }

  /* Compute the reference height                                    */
  if (wave->topk_wc > wave->botk_wc)
    zr = (wave->bot_l - wave->bot) / 2.;
  else
    zr = (wave->top - wave->bot) / 2.;

  /* Compute the apparent Z0.                                        */
  madsen94o(wave->ij[0], wave->ij[1], ub, 2 * M_PI / period, u, zr, 
	    dir - udir, wave->z0, &ustrc, &ustrwm, &ustrr, &fwc, &zoa);

  /* Store bottom friction velocity FIX - Is this the right one???   */
  if (!isnan(ustrr))
    *wave->ustrcw = ustrr;
  else
    *wave->ustrcw = 0.0;

  v = log(zr / zoa) / VON_KAR;

  if (wave->do_Cd)
    *wave->Cd = max(wave->quad_bfc, 1.0 / (v * v));
}


/*-------------------------------------------------------------------*/
/* Evaluate wave components given an amplitude & direction in deg T  */
/*-------------------------------------------------------------------*/
void wavecomp_estim(double a, double dir, double *wx, double *wy)
{
  dir = (dir < 0.) ? dir + 360. : dir;
  dir = (dir > 360.) ? dir - 360. : dir;
  *wy = a * cos(dir * PI / 180.);
  *wx = sqrt(a * a - (*wy) * (*wy));
  *wx = (dir > 180.) ? -*wx : *wx;
}


/*-------------------------------------------------------------------*/
/* Evaluate surface wave direction.                                  */
/*   Input: wind stress components                                   */
/*   Output: wave direction                                          */
/* Author NMY                                                        */
/*-------------------------------------------------------------------*/
double wavedir_estim(double wy_str, double wx_str)
{
  double wavedirection;
  wavedirection = fmod(450.0 - atan2(wy_str, wx_str) * 180.0 / PI, 180.0);
  return (wavedirection);
}


/*-------------------------------------------------------------------*/
/* Find surface wave number using shallow water approximation.       */
/*   Input: wave period and water depth                              */
/*   Output: wave number                                             */
/* Author NMY                                                        */
/*-------------------------------------------------------------------*/
double wavenumb_estim(double period, double depth)
{
  double wavenumber, wn1, wn2, wn3, er1,er2,er3;
  /* wavenumber = (2. * PI / period) / sqrt(g * depth); shallow water approx*/
  if (period == 0.0) return(0.0);
  wn1 = (2. * PI / period) / sqrt(g * depth);
  wn2 = (2. * PI / period)*(2. * PI / period) / g ;
  wn3 = 0.5*(wn1+wn2);
  er1 = fabs(atanh(wn1*depth)-(2. * PI / period)*(2. * PI / period) / (g*wn1));
  er2 = fabs(atanh(wn2*depth)-(2. * PI / period)*(2. * PI / period) / (g*wn2));
  er3 = fabs(atanh(wn3*depth)-(2. * PI / period)*(2. * PI / period) / (g*wn3));

  if (er2<er1) {wavenumber= wn2; er1=er2;}
  else {wavenumber= wn1; er2=er1;}
  if (er3<er1) wavenumber = wn3;
 
  return (wavenumber);
}


/*-------------------------------------------------------------------*/
/* Find near bottom orbital velocity using linear wave theory.       */
/*   Input: wave amplitude, period, number, and water depth          */
/*   Output: Near-bottom wave orbital velocity                       */
/* Author NMY                                                        */
/*-------------------------------------------------------------------*/
double waveub_estim(double wamp, double period, double wavenumber,
                    double depth)
{
  double waveub;
  if (wavenumber == 0.0) return(0.0);
  waveub = (1.0 * wamp * (PI / period)) / sinh(wavenumber * depth);
  return (waveub);
}


/*-------------------------------------------------------------------*/
/* Evaluate surface wave amplitude.                                  */
/* (Y.M.Tang,R.Grimshaw 1996, JGR Vol.101,No.C10, pp. 22,705-22,714) */
/*   Input: wind stress components, wave number, period, and depth   */
/*   Output: Wave amplitude                                          */
/* Author NMY                                                        */
/*-------------------------------------------------------------------*/
double wamp_estim(double wx_str, double wy_str, double wavenumber,
                  double period)
{
  double wamp;
  double wa_c, wa_diss;
  double w_str, w_ustr;
  wa_c = (2. * 3.14 / period) / wavenumber; /* wave phase speed */
  wa_diss = 100.;               /* Dissipation constant */
  w_str = sqrt(wx_str * wx_str + wy_str * wy_str);  /* wind stress module */
  w_ustr = sqrt(w_str);         /* wind shear velocity */
  if (28. * w_ustr > wa_c) {
    wamp = 0.2 * 0.001 * fabs(28. * (w_ustr / wa_c) - 1.) / wa_diss;
    wamp = sqrt(sqrt(wamp));
    wamp = wamp / wavenumber;
  } else
    wamp = 0.;
  return (wamp);                /* wave amplitude */
}


/*-------------------------------------------------------------------*/
/* Wind waves according to Toba (1978) Stochastic form of the growth */
/* of wind waves in a single-parameter representation with physical  */
/* implications. J. Phys. Oceanogr., 8, 494-507.                     */
/*-------------------------------------------------------------------*/
void wind_wave_toba(wave_t *wave, double *f)
{
  double Hs;           /* Significant wave height (m)                */
  double period;       /* Wave period (sec)                          */
  double ustar;        /* Wind friction velocity                     */
  double dir;          /* Wind direction (deg T)                     */
  double fetch;        /* Fetch (m)                                  */
  double H, T;         /* Dummies                                    */

  dir = wavedir_estim(wave->wy, wave->wx);
  fetch = get_fetch(f, dir);
  ustar =  sqrt(sqrt(wave->wx * wave->wx + wave->wy * wave->wy) / 1.225);
  H = (ustar) ? 0.05 * sqrt(g * fetch / (ustar * ustar)) : 0.0;
  Hs = H * ustar * ustar / g;
  *wave->amp = 0.5 * Hs;
  T = pow(6.4 * H, 0.666667);
  *wave->period = T * ustar / g;
  *wave->dir = dir;
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Wind waves according to US Army Core of Engineers.                */
/* This has been documented in Hipsey et al (2006) Computational     */
/* Aquatic Ecosystem Dynamics Model, v2.3 Science Manual. CWR, UWA.  */
/* http://www.cwr.uwa.edu.au/services/models/legacy/model/caedym/    */
/*-------------------------------------------------------------------*/
void wind_wave_usac(wave_t *wave, double *f)
{
  double Hs;           /* Significant wave height (m)                */
  double period;       /* Wave period (sec)                          */
  double ustar;        /* Wind friction velocity                     */
  double dir;          /* Wind direction (deg T)                     */
  double fetch;        /* Fetch (m)                                  */
  double H, T;         /* Dummies                                    */
  double w1, w2;       /* Wind speed (ms-1)                          */

  dir = wavedir_estim(wave->wy, wave->wx);
  fetch = get_fetch(f, dir);
  w1 = wave->wx;
  w2 = wave->wy;
  stress2wind(&w1, &w2);
  w2 = w1 * w1 + w2 * w2;

  w1 = 0.530 * pow(g * wave->depth_wc / w2, 0.75);
  Hs = w2 * 0.283 * tanh(w1) * 
    tanh(0.00565 * sqrt(g * fetch / w2) / tanh(w1)) / g;

  w1 = 0.833 * pow(g * wave->depth_wc / w2, 0.375);
  period = 2.0 * PI * sqrt(w2) * 1.20 * tanh(w1) *
    tanh(0.0379 * pow(g * fetch / w2, 0.333) / tanh(w1)) / g;

  *wave->amp = 0.5 * Hs;
  *wave->period = period;
  *wave->dir = dir;
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Wave-current friction factor.                                     */
/* (Equations 32 and 33 in Madsen, 1994).                            */
/* @param cmu Relative strength of currents from Eqn 27              */
/* @param cukw Relative roughness = cmu*ubr/(kn*wr)                  */
/* @return Wave-current friction factor                              */
/* @author Chris Sherwood, CSIRO                                     */
/* @version 11 July 1997                                             */
/*-------------------------------------------------------------------*/
static double fwc94(double cmu, double cukw)
{
  if (cukw <= 100.)
    return (cmu * exp(7.02 * pow(cukw, -0.078) - 8.82));
  else                          /* ( cukw > 100. ) */
    return (cmu * exp(5.61 * pow(cukw, -0.109) - 7.30));
}

/*-------------------------------------------------------------------*/
/* Grant-Madsen wave-current calculations from Madsen (1994).        */
/* This version supercedes madsen94 (it is better optimized) but the */
/* argument list is differant.                                       */
/*  (Note that this formulation does not address angular difference  */
/*  between stress and velocity).                                    */
/* @param c Cell number.                                             */
/* @param ub Wave-orbital velocity outside wbl [m/s]                 */
/* @param T Wave period [s]                                          */
/* @param ur Current velocity at height zr [m/s]                     */
/* @param zr Reference height for current velocity [m]               */
/* @param phiwc Angle between currents and waves at zr (degrees)     */
/* @param zo Bottom roughness height [m]                             */
/* @param pustrc Current shear velocity u*c [m/s]                    */
/* @param pustrr Wave shear velocity u*w [m/s]                       */
/* @param pphicw Angle between u*cw and waves at bed (degrees)       */
/* @param pustrcw Wave-current combined shear velocity u*cw [m/s]    */
/* @param pfwc Wave friction factor [ ]                              */
/* @param pzoa Apparent bottom roughness [m]                         */
/* @author Chris Sherwood, CSIRO                                     */
/* @version 5 June 1997                                              */
/*-------------------------------------------------------------------*/
static void madsen94o(int ii, int jj, double ubr, double wr, double ucr,
		      double zr, double phiwc, double zo,
		      double *pustrc, double *pustrwm,
		      double *pustrr, double *pfwc, double *pzoa)
{
#define MAXIT 20
#define LOWU (0.01)
  double rmu[MAXIT], Cmu[MAXIT], fwc[MAXIT], dwc[MAXIT];
  double ustrwm2[MAXIT], ustrr2[MAXIT], ustrc[MAXIT];
  double kN, cosphiwc, ustrr, lnzr, lndw, lnln, bigsqr, diff, zoa;
  int i, nit, errflg;

  kN = 30. * zo;
  zoa = zo;

  /* junk return values */
  *pustrc = 0.0;
  *pustrwm = 0.0;
  *pustrr = 0.0;
  *pfwc = .4;
  *pzoa = zoa;

  /* some data checks */
  if (wr <= 0.) {
    emstag(LPANIC,"waves:madsen94o","Bad value for ang. freq. in Madsen94o at (%d,%d): wr = %g\n", ii, jj, wr);
    return;
  }

  if (ubr < 0.) {
    emstag(LPANIC,"waves:madsen94o","Bad value for orbital vel. in Madsen94o at (%d,%d): ub = %g\n",ii, jj, ubr);
    return;
  }
  if (kN <= 0.) {
    emstag(LPANIC,"waves:madsen94o","Negative roughness in Madsen94o at (%d,%d): kN = %g\n", ii, jj, kN);
    return;
  }
  if (zr < 5. * kN) {
    emstag(LPANIC,"waves:madsen94o","Low value for ref. level in Madsen94o at (%d,%d): zr = %g\n", ii, jj, zr);
    return;
  }
  if (ubr <= LOWU) {
    if (ucr <= LOWU) {
      *pustrc = 0.;
      *pustrwm = 0.;
      *pustrr = 0.;
      *pzoa = zo;

      return;
    }
    ustrc[0] = ucr * VON_KAR / log(zr / zo);
    *pustrc = ustrc[0];
    *pustrwm = 0.;
    *pustrr = ustrc[0];
    *pzoa = zo;
    return;
  }

  cosphiwc = fabs(cos(phiwc * (M_PI / 180.)));
  rmu[0] = 0.;
  Cmu[0] = 1.;
  fwc[0] = fwc94(Cmu[0], (Cmu[0] * ubr / (kN * wr))); /* Eqn. 32 or 33 */
  ustrwm2[0] = 0.5 * fwc[0] * ubr * ubr;  /* Eqn. 29 */
  ustrr2[0] = Cmu[0] * ustrwm2[0];  /* Eqn. 26 */
  ustrr = sqrt(ustrr2[0]);
  dwc[0] = (Cmu[0] * ubr / (kN * wr)) >= 8. ? 2. * VON_KAR * ustrr / wr : kN;
  lnzr = log(zr / dwc[0]);
  lndw = log(dwc[0] / zo);
  lnln = lnzr / lndw;
  bigsqr =
    (-1. + sqrt(1 + ((4. * VON_KAR * lndw) / (lnzr * lnzr)) * ucr / ustrr));
  ustrc[0] = 0.5 * ustrr * lnln * bigsqr;
  errflg = 0;
  for (i = 1, nit = 1, diff = 1.; diff > 0.0005 && i < MAXIT; i++, nit++) {
    rmu[i] = ustrc[i - 1] * ustrc[i - 1] / ustrwm2[i - 1];  /* Eqn. 28 */
    Cmu[i] = sqrt(1. + 2. * rmu[i] * cosphiwc + rmu[i] * rmu[i]); /* Eqn
                                                                     27 */
    fwc[i] = fwc94(Cmu[i], (Cmu[i] * ubr / (kN * wr))); /* Eqn. 32 or 33 */
    ustrwm2[i] = 0.5 * fwc[i] * ubr * ubr;  /* Eqn. 29 */
    ustrr2[i] = Cmu[i] * ustrwm2[i];  /* Eqn. 26 */
    ustrr = sqrt(ustrr2[i]);
    dwc[i] = (Cmu[i] * ubr / (kN * wr)) >= 8. ? 2. * VON_KAR * ustrr / wr : kN;  /* Eqn. 
                                                                               36 
                                                                             */
    if (dwc[i] > 0.8 * zr) {
      errflg = 1;
      dwc[i] = 0.8 * zr;        /* clumsy fix for ill posed cases */
    }
    lnzr = log(zr / dwc[i]);
    lndw = log(dwc[i] / zo);
    lnln = lnzr / lndw;
    bigsqr =
      (-1. + sqrt(1 + ((4. * VON_KAR * lndw) / (lnzr * lnzr)) * ucr / ustrr));
    ustrc[i] = 0.5 * ustrr * lnln * bigsqr; /* Eqn. 38 */
    diff = fabs((fwc[i] - fwc[i - 1]) / fwc[i]);
  }
  if (errflg)
    emstag(LPANIC,"waves:madsen94o","dwc > 0.8*zr\n");

  *pustrwm = sqrt(ustrwm2[nit - 1]);
  *pustrc = ustrc[nit - 1];
  *pustrr = ustrr;
  *pfwc = fwc[nit - 1];
  zoa =
    exp(log(dwc[nit - 1]) -
        (ustrc[nit - 1] / ustrr) * log(dwc[nit - 1] / zo));
  *pzoa = zoa;
  return;
}


/*-------------------------------------------------------------------*/
/* Routine to return the fetch given a wind direction                */
/*-------------------------------------------------------------------*/
double get_fetch(double *fetch, double dir)
{

  if (dir < 0.0) dir += 360.0;

  if (dir < 30 || dir >= 330)
    return (fetch[4]);
  else if (dir >= 30 && dir < 60)
    return (fetch[5]);
  else if (dir >= 60 && dir < 120)
    return (fetch[6]);
  else if (dir >= 120 && dir < 150)
    return (fetch[7]);
  else if (dir >= 150 && dir < 210)
    return (fetch[0]);
  else if (dir >= 210 && dir < 240)
    return (fetch[1]);
  else if (dir >= 240 && dir < 300)
    return (fetch[2]);
  else if (dir >= 300 && dir < 330)
    return (fetch[3]);
  else
    return (0.0);
}

/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate wind speed given wind stress                 */
/*-------------------------------------------------------------------*/
void stress2wind(double *wx, double *wy)
		 
{
  double cd = 1.4e-3, w0;
  double v0 = 10.0;
  double v1 = 26.0;
  double c0 = 0.00114;
  double c1 = 0.00218;
  double air_dens = 1.225;
  double sx = *wx < 0. ? -1.0 : 1.0;
  double sy = *wy < 0. ? -1.0 : 1.0;

  w0 = sqrt(fabs(*wx) / (cd * air_dens));
  if (w0 <= v0)
    cd = c0;
  else if (w0 >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (w0 - v0) / (v1 - v0);
  *wx = sx * sqrt(fabs(*wx) / (cd * air_dens));

  w0 = sqrt(fabs(*wy) / (cd * air_dens));
  if (w0 <= v0)
    cd = c0;
  else if (w0 >= v1)
    cd = c1;
  else
    cd = c0 + (c1 - c0) * (w0 - v0) / (v1 - v0);
  *wy = sy * sqrt(fabs(*wy) / (cd * air_dens));
}

/*-------------------------------------------------------------------*/
