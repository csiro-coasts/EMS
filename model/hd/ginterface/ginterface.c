/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hds/ginterface/ginterface.c
 *  
 *  Description:
 *  Generic interface
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ginterface.c 7475 2023-12-16 11:17:08Z bai155 $
 *
 */

#include <stdio.h>
#include "hd.h"

#define NOT   -0x000001
#define ETA   -0x000002
#define KZ    -0x000004
#define VZ    -0x000008
#define U1VH  -0x000010
#define U2VH  -0x000020
#define U1KH  -0x000040
#define U2KH  -0x000080
/*!!!  add above definitions to below */
int LOCAL_TRACERS[] = {NOT, ETA, KZ, VZ,U1VH, U2VH, U1KH, U2KH};
int NLOCAL_TRACERS = sizeof(LOCAL_TRACERS) / sizeof(int);

char* ginterface_get2Dtracername(void* model, int i);
int  ginterface_get_max_numbercolumns(void* model);

static int e_ntr;
static int tr_map[MAXNUMVARS];
static int sed_map[MAXNUMVARS];
static int e_nepi;
static int epi_map[MAXNUMVARS];

int i_get_c(void *model, int b);

#if defined(HAVE_TRACERSTATS_MODULE)

trs_t* trs_create();

#define WIN_SUFFIX "wn"
/*-------------------------------------------------------------------*/
/* Tracer statistics interface step                                  */
/*-------------------------------------------------------------------*/
void tracerstats_step(geometry_t *window)
{
  win_priv_t *wincon = window->wincon;

  /*-----------------------------------------------------------------*/
  /* Do the interface routines */
  int cc, c;

  process_cell_mask(window, window->ctrst, window->nctrst);
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    /* Exclude process points if required */
    if (wincon->c2[c]) continue;
    trs_step(window, wincon->trs, c, LASTSTEP);
  }
}


/*-------------------------------------------------------------------*/
/* Tracer statistics interface step                                  */
/*-------------------------------------------------------------------*/
void tracerstats_prestep(geometry_t *window,int step)
{
  win_priv_t *wincon = window->wincon;

  /*-----------------------------------------------------------------*/
  /* Do the interface routines */
  int cc, c;
  if(step >= MAXSTEPS )
  	  hd_warn("ginterface : tracerstats-prestep index out of bounds %u - IGNORING!",step);

  emstag(LTRACE,"hd:ginterface:tracerstats_prestep","Executing tracerstats prestep %u",step);
  process_cell_mask(window, window->ctrst, window->nctrst);
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    /* Exclude process points if required */
    if (wincon->c2[c]) continue;
    trs_step(window, wincon->trs, c, step);
  }
}

#endif

/* Ecolgy routines
int einterface_get_eco_flag(void* model, char* name)
double sinterface_gettracersvel(void* model, char* name)
*/

/* Sediment routines
int sinterface_getcss_er_val(FILE* prmfd, double *css_er_val)
int sinterface_getcss_er_depth(FILE* prmfd, double *css_er_depth)
void sinterface_putgridz_sed(void* hmodel, int c, double *gridz_sed)
void sinterface_putcellz_sed(void* hmodel, int c, double *cellz_sed)
void sinterface_putgriddz_wc(void* hmodel, int c, double *gridz_wc)
void sinterface_putcellz_wc(void* hmodel, int c, double *cellz_wc)
void sinterface_putdz_wc(void* hmodel, int c, double *dz_wc)
void sinterface_puttopz_wc(void* hmodel, int c, double topz_wc)
void sinterface_putbotz_wc(void* hmodel, int c, double botz_wc)
int sinterface_put_ustrcw(void* hmodel, int c, double ustrcw)
int sinterface_getverbose_sed(FILE* prmfd)
int sinterface_getgeomorph(FILE* prmfd)
int sinterface_getconsolidate(FILE* prmfd)
double sinterface_getfinpor_sed(FILE* prmfd)
double sinterface_getconsolrate(FILE* prmfd)
double sinterface_getmaxthicksed(FILE* prmfd)
int sinterface_getcssmode(FILE* prmfd)
double sinterface_getcss(FILE* prmfd)
double sinterface_getcssdep(FILE* prmfd)
int sinterface_getflocmode(FILE* prmfd)
double sinterface_getflocprm1(FILE* prmfd)
double sinterface_getflocprm2(FILE* prmfd)
int sinterface_getbblnonlinear(FILE* prmfd)
int sinterface_gethindered_svel_patch(FILE* prmfd)
int sinterface_gethindered_svel(FILE* prmfd)
double sinterface_reef_scale_depth(FILE* prmfd)
int sinterface_getcalcripples(FILE* prmfd)
double sinterface_getcssscale(FILE* prmfd)
double sinterface_getphysriph(FILE* prmfd)
double sinterface_getphysripl(FILE* prmfd)
double sinterface_getbioriph(FILE* prmfd)
double sinterface_getbioripl(FILE* prmfd)
double sinterface_getbiodens(FILE* prmfd)
double sinterface_getmaxbiodepth(FILE* prmfd)
double sinterface_getbi_dissol_kz(FILE* prmfd)
double sinterface_getbt_partic_kz(FILE* prmfd)
double sinterface_getbi_dissol_kz_i(FILE* prmfd)
double sinterface_getbt_partic_kz_i(FILE* prmfd)
char sinterface_getbiosedprofile(FILE* prmfd)
double sinterface_getz0(FILE* prmfd)
double sinterface_getquad_bfc(FILE* prmfd)
int sinterface_getshipfile(FILE* prmfd, char *shipfile)
double sinterface_erflux_scale(FILE* prmfd)
FILE* si_getparamfile_tracer(FILE* fp)
FILE* si_getparamfile_sed(FILE *fp)
double sinterface_get_psize(void* model, char *name)
double sinterface_get_b_dens(void* model, char *name)
double sinterface_get_i_conc(void* model, char *name)
double sinterface_get_svel(void* model, char *name)
void sinterface_get_svel_name(void* model, char *name, char *sname)
int sinterface_get_cohesive(void* model, char *name)
int sinterface_get_resuspend(void* model, char *name)
int sinterface_get_deposit(void* model, char *name)
int sinterface_get_floc(void* model, char *name)
int sinterface_get_calcvol(void* model, char *name)
int sinterface_get_adsorb(void* model, char *name)
void sinterface_get_carriername(void* model, char *name,char *carriername )
void sinterface_get_dissolvedname(void* model, char *name, char *dissolvedname)
double sinterface_get_adsorb_kd(void* model, char *name)
double sinterface_get_adsorb_rate(void* model, char *name)
*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Generic interface functions. General functions are prefixed with  */
/* i_. Routines prefixed by ginterface_ are typically used by        */
/* ecology or sediments. There exist some wrapper functions in this  */
/* class for convenience that call the more generic i_ functions.    */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

/* Routines specific to SHOC
double i_get_sinthcell(void *hmodel, int cc)
double i_get_costhcell(void *hmodel, int cc)
double i_get_wave_wind1(void *hmodel, int c)
double i_get_wave_wind2(void *hmodel, int c)
void i_get_bot_vel(void *hmodel, double sinthcell, double costhcell,
		   double *u1bot, double *u2bot, double *botz, int c)
void i_get_tracer_e1(void* hmodel, int c, int ntr, int *tmap, double **tr_e1)
void i_get_tracer_e2(void* hmodel, int c, int ntr, int *tmap, double **tr_e1)
*/

/*-------------------------------------------------------------------*/
/* 3D (water column) tracers                                         */
/*-------------------------------------------------------------------*/
/* Returns 1 if itr is a valid tracer or state variable defined in   */
/* the array LOCAL_TRACERS.                                          */
int i_is_valid_tracer(int itr)
{
  int i;
  if(itr >= 0)
    return (1);
  
  for(i = 0;i< NLOCAL_TRACERS;i++ )
    {
      if(LOCAL_TRACERS[i] == itr )
	return (1);
    }
  
  return (0);
}

/* Returns 1 if 3D tracer 'name' exists                              */
int i_tracername_exists(void* model, char*name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn;

  for(tn=0; tn<wincon->ntr; tn++) {
    /* Water column */
    if((tracer_find_index(name, wincon->ntr,
              wincon->trinfo_3d) >= 0))
      return 1;
  }
  return 0;
}
/* Wrapper to 3D tracer name exists                                  */
int ginterface_tracername_exists(void* model, char*name)
{
  return i_tracername_exists(model, name);
}

/* Returns 1 if itr is a valid sediment tracer                       */
int i_is_valid_sed_tracer(int itr)
{
  int i;
  if(itr >= 0)
    return (1);
  
  return (0);
}

/* Gets the number of explicit and auto 3D tracers                   */
int i_get_num_tracers_3d(void* hmodel, int *atr)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    *atr = master->atr;
    return wincon->ntr;
}

/* Returns 1 if the model grid is geographic                         */
int i_is_geographic(void* hmodel)
{
  geometry_t* window = (geometry_t*) hmodel;
  return(window->is_geog);
}

/* Returns the number of 3D tracers                                  */
int ginterface_getntracers(void* model)
{
    return 0;
}

/* Gets an array of 3D tracer names                                  */
void i_get_names_tracers_3d(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->ntr; n++) {
      strcpy(trname[n], wincon->trinfo_3d[n].name);
    }
}

/* Gets an array of 3D tracer names                                  */
void i_get_tracerstats_3d(void* hmodel, char *trname, int tn) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    strcpy(trname, wincon->trinfo_3d[tn].tracerstat);
}

/* Returns the name of a 3D tracer i                                 */
char* ginterface_gettracername(void* model, int i)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  if(i >= 0 && i < wincon->ntr)
    return wincon->trinfo_3d[i].name;

  return NULL;
}

/* Gets the index in the tracer info list of 3D tracer n             */
int i_get_param_map_3d(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    return wincon->trinfo_3d[n].m;
}

/* Returns the number of 3D tracer with name 'name'                  */
/* Can also use tracer_find_index().                                 */
int gi_getmap2hdwctracer(void *hmodel, int ns, char *name)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon=window->wincon;
    int n;
    /*int index = tracer_find_index(name, wincon->ntr, wincon->trinfo_3d);*/
    for(n = 0; n < windat->ntr; n++)
	if(strcmp(name, wincon->trinfo_3d[n].name) == 0)
	    return n;
    return (-1);
}

/* Returns the tracer info for a 3D tracer with name 'trname'.       */
/* Returns NULL if the tracer can't be found.                        */
tracer_info_t* i_get_tracer(void* hmodel, char* trname)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int n;
  emstag(LTRACE,"hd:ginterface:i_get_tracer","looking for tracer with name: %s",trname);


  for (n = 0; n < wincon->ntr; n++) {
    if(strcmp(trname,wincon->trinfo_3d[n].name) == 0)
      return &(wincon->trinfo_3d[n]);
  }
  return NULL;
}

/* Returns the tracer info for a epi with name 'trname'.       */
/* Returns NULL if the tracer can't be found.                        */
tracer_info_t* i_get_tracer2d(void* hmodel, char* trname)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int n;
  emstag(LTRACE,"hd:ginterface:i_get_tracer2d","looking for epi with name: %s",trname);


  for (n = 0; n < wincon->ntrS; n++) {
    if(strcmp(trname,wincon->trinfo_2d[n].name) == 0)
      return &(wincon->trinfo_2d[n]);
  }
  return NULL;
}

/* Returns a list of 3D tracer indexes corresponding to tracer names */
/* trname[].                                                         */
int *i_get_tmap_3d(void* model, int ntr, char *trname[])
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn;
  int *tmap_3d = i_alloc_1d(ntr);

  for(tn=0; tn<ntr; tn++) {
    if(strcmp(trname[tn], "Kz") == 0)
      tmap_3d[tn] = KZ;
    else if(strcmp(trname[tn], "Vz") == 0)
      tmap_3d[tn] = VZ;
    else if(strcmp(trname[tn], "u1vh") == 0)
      tmap_3d[tn] = U1VH;
    else if(strcmp(trname[tn], "u2vh") == 0)
      tmap_3d[tn] = U2VH;
    else if(strcmp(trname[tn], "u1kh") == 0)
      tmap_3d[tn] = U1VH;
    else if(strcmp(trname[tn], "u2kh") == 0)
      tmap_3d[tn] = U2VH;
    else {
      if((tmap_3d[tn] = tracer_find_index(trname[tn], wincon->ntr,
            wincon->trinfo_3d)) < 0)
	hd_warn("ginterface : Can't find type WC3D tracer '%s' in parameter file at %d.\n", trname[tn],tn);
    }
  }
  return(tmap_3d);
}

/* Returns the fill value for 3D tracer 'name'                       */
double ginterface_get_fillvalue_wc(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    
    double v = 0;
    v = tr->fill_value_wc;
    return v;
}

/* Returns the diagnostic flag of 3D tracer 'name'                   */
int ginterface_gettracerdiagnflag(void* model, char* name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int index = tracer_find_index(name, wincon->ntr, wincon->trinfo_3d);
  return wincon->trinfo_3d[index].diagn;
}
/* Wrapper for diagnostic flag                                       */
int ginterface_get_diagn(void* model, char *name)
{
    return ginterface_gettracerdiagnflag(model, name);
}

/* Returns the dissolved flag of a 3D tracer                         */
int ginterface_get_dissol(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    
    int v = 1;
    v = tr->dissol;
    return v;
}

/* Returns the particle flag of 3D tracers                           */
int ginterface_gettracerparticflag(void* model, char* name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int index = tracer_find_index(name, wincon->ntr, wincon->trinfo_3d);
  return wincon->trinfo_3d[index].partic;
}
/* Wrapper for particulate flag                                      */
int ginterface_get_partic(void* model, char *name)
{
    return ginterface_gettracerparticflag(model, name);
}

/* Returns the diffuse flag of 3D tracers                            */
int ginterface_get_diffuse(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    
    int v = 1;
    v = tr->diffuse;
    return v;
}

/* Returns the decay flag of 3D tracers                              */
double ginterface_get_decay(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    
    double v = 0;
    v = (tr->flag & (DE_TR2|DE_TR3)) ? 0.0 : atof(tr->decay);
    return v;
}

/* Returns the units flag of 3D tracers                              */
void ginterface_get_tracerunits(void* model, char *name, char *units)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    strcpy(units, tr->units);
}

/* Gets 2D pointers to an array of 3D tracers at coordinate c,       */
/* including 3D state variable tracers.                              */
void i_get_tracer_wc(void* hmodel, int c, int ntr, int *tmap, double ***tr_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon = window->wincon;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k, n, m;
    for(k = top_k; k >= bot_k; k--) {
      for(n = 0; n < ntr; n++) {
        m = tmap[n];
        if(m == NOT)
          continue;
        else if(m == KZ)
          tr_wc[n][k] = &windat->Kz[c3];
        else if(m == VZ)
          tr_wc[n][k] = &windat->Vz[c3];
        else if(m == U1VH)
          tr_wc[n][k] = &wincon->u1vh[c3];
        else if(m == U2VH)
          tr_wc[n][k] = &wincon->u2vh[c3];
        else if(m == U1KH)
          tr_wc[n][k] = &wincon->u1kh[c3];
        else if(m == U2KH)
          tr_wc[n][k] = &wincon->u2kh[c3];
        else {
          tr_wc[n][k] = &windat->tr_wc[m][c3];
	}
      }
      c3 = window->zm1[c3];
    }
}
/* Returns 1D pointers to an array of tracers at index b.            */
double**  ginterface_getwctracers(void* model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nz = window->nz;
  int ntr = e_ntr;
  double **wctr = (double **)calloc(ntr * nz, sizeof(double*));
  int cs2 = window->wincon->s2[b+1];
  int c2 = window->m2d[cs2];
  int cc = window->c2cc[c2];
  int cs = window->nsur_t[cc];
  int cb = window->bot_t[cc];
  int c, k = window->s2k[cb], n;

  /*
   * Loop from the bottom cell up
   */
  for (c = cb; c != cs; c = window->zp1[c]) {
    for (n = 0; n < ntr; ++n) {
      wctr[k * ntr + n] = &windat->tr_wc[tr_map[n]][c];
    }
    k++;
  }
  // The surface cell
  for (n = 0; n < ntr; ++n) {
    wctr[k * ntr + n] = &windat->tr_wc[tr_map[n]][c];
  }

  return wctr;
}

/* Gets 2D pointers to an array of 3D tracers at coordinate c on the */
/* e1 face, including 3D cell centered state variable tracers.       */
/* Only valid for the structured model.                              */
void i_get_tracer_e1(void* hmodel, int c, int ntr, int *tmap, double **tr_e1)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int xm1 = window->xm1[c3];
  int k, n, m;
  for(k = top_k; k >= bot_k; k--) {
    for(n = 0; n < ntr; n++) {
      m = tmap[n];
      if(m == NOT)
	continue;
      else if(m == KZ)
	tr_e1[n][k] = 0.5 * (windat->Kz[c3] + windat->Kz[xm1]);
      else if(m == VZ)
	tr_e1[n][k] = 0.5 * (windat->Vz[c3] + windat->Vz[xm1]);
      else if(m == U1VH)
	tr_e1[n][k] = wincon->u1vh[c3];
      else if(m == U2VH)
	tr_e1[n][k] = wincon->u2vh[c3];
      else if(m == U1KH)
	tr_e1[n][k] = wincon->u1kh[c3];
      else if(m == U2KH)
	tr_e1[n][k] = wincon->u2kh[c3];
      else
	tr_e1[n][k] = 0.5 * (windat->tr_wc[m][c3] + windat->tr_wc[m][xm1]);
    }
    c3 = window->zm1[c3];
    xm1 = window->xm1[c3];
  }
}

/* Gets 2D pointers to an array of 3D tracers at coordinate c on the */
/* e2 face, including 3D cell centered state variable tracers.       */
/* Only valid for the structured model.                              */
void i_get_tracer_e2(void* hmodel, int c, int ntr, int *tmap, double **tr_e2)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int ym1 = window->ym1[c3];
  int k, n, m;
  for(k = top_k; k >= bot_k; k--) {
    for(n = 0; n < ntr; n++) {
      m = tmap[n];
      if(m == NOT)
	continue;
      else if(m == KZ)
	tr_e2[n][k] = 0.5 * (windat->Kz[c3] + windat->Kz[ym1]);
      else if(m == VZ)
	tr_e2[n][k] = 0.5 * (windat->Vz[c3] + windat->Vz[ym1]);
      else if(m == U1VH)
	tr_e2[n][k] = wincon->u1vh[c3];
      else if(m == U2VH)
	tr_e2[n][k] = wincon->u2vh[c3];
      else if(m == U1KH)
	tr_e2[n][k] = wincon->u1kh[c3];
      else if(m == U2KH)
	tr_e2[n][k] = wincon->u2kh[c3];
      else
	tr_e2[n][k] = 0.5 * (windat->tr_wc[m][c3] + windat->tr_wc[m][ym1]);
    }
    c3 = window->zm1[c3];
    ym1 = window->ym1[c3];
  }
}

/* Gets 2D pointers to an array of 3D tracers at coordinate c on the */
/* gridz face, including 3D cell centered state variable tracers.    */
void i_get_tracer_w(void* hmodel, int c, int ntr, int *tmap, double **tr_w)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int zm1 = window->zm1[c3];
  int k, n, m;
  for(k = top_k; k >= bot_k; k--) {
    for(n = 0; n < ntr; n++) {
      m = tmap[n];
      if(m == NOT)
	continue;
      else if(m == KZ)
	tr_w[n][k] = 0.5 * (windat->Kz[c3] + windat->Kz[zm1]);
      else if(m == VZ)
	tr_w[n][k] = 0.5 * (windat->Vz[c3] + windat->Vz[zm1]);
      else if(m == U1VH)
	tr_w[n][k] = 0.5 * (wincon->u1vh[c3] + wincon->u1vh[zm1]);
      else if(m == U2VH)
	tr_w[n][k] = 0.5 * (wincon->u2vh[c3] + wincon->u2vh[zm1]);
      else if(m == U1KH)
	tr_w[n][k] = 0.5 * (wincon->u1kh[c3] + wincon->u1kh[zm1]);
      else if(m == U2KH)
	tr_w[n][k] = 0.5 * (wincon->u2kh[c3] + wincon->u2kh[zm1]);
      else
	tr_w[n][k] = 0.5 * (windat->tr_wc[m][c3] + windat->tr_wc[m][zm1]);
    }
    c3 = window->zm1[c3];
    zm1 = window->zm1[c3];
  }
}

/* Returns the surface flux of 3D tracer 'name' at coordinate c      */
double ginterface_get_srf_flux(void* hmodel, char *name, int c)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  win_priv_t *wincon=window->wincon;
  tracer_info_t *tr = i_get_tracer(hmodel, name);
  int n = tr->n;
  int c2;
  double v;
  if(wincon->sflux[n]>=0) {
  	c2=window->m2d[c];
  	v = windat->tr_wcS[wincon->sflux[n]][c2];
	//	fprintf(stderr,"get_srf_flux=%lf \n", v);
	v = (windat->eta[c2] - window->botz[c2] <= wincon->hmin) ? 0.0 : v;
  	return v;
  }
  else
    return 0;
}

/* Returns the size of tracers to transfer to SWAN                   */
int i_get_swan_size(void* hmodel)
{
  return (0);
}
int i_get_swan_size_m(void* hmodel)
{
  return (0);
}

int i_get_swan_nobc_m(void* hmodel)
{
  return (0);
}
int *i_get_swan_obc_m(void* hmodel)
{
  return (NULL);
}

/* Checks the status of the hotstart flag */
int i_check_swan_hs_m(void* hmodel)
{
  return (0);
}

int i_check_swan_hs(void* hmodel)
{
  return (0);
}

/* Returns a poiner to the sea level tracer to transfer to SWAN      */
double *i_get_swan_eta(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_eta_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the UAV tracer to transfer to SWAN            */
double *i_get_swan_uav(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_uav_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the VAV tracer to transfer to SWAN            */
double *i_get_swan_vav(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_vav_m(void* hmodel)
{
  return (NULL);
}

/* Returns a pointer to the east component of cell centered wind     */
/* stress. Unstructured meshes only; assumes the centered u wind     */
/* stress is in wincon->d2 via vel_cen().                            */
double *i_get_swan_wx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wx_m(void* hmodel)
{
  return (NULL);
}

/* Returns a pointer to the north component of cell centered wind    */
/* stress. Unstructured meshes only; assumes the centered v wind     */
/* stress is in wincon->d3 via vel_cen().                            */
double *i_get_swan_wy(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wy_m(void* hmodel)
{
  return (NULL);
}


/* Returns a poiner to the bathymetry to transfer to SWAN            */
double *i_get_swan_dep(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_dep_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave amplitude tracer to get from SWAN    */
double *i_get_swan_amp(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_amp_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave period tracer to get from SWAN       */
double *i_get_swan_per(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_per_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave direction tracer to get from SWAN    */
double *i_get_swan_dir(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_dir_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave orbital velocity tracer to get from  */
/* SWAN.                                                             */
double *i_get_swan_ub(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_ub_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the x wave induced force tracer to get from   */
/* SWAN.                                                             */
double *i_get_swan_Fx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_Fx_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the y wave induced force tracer to get from   */
/* SWAN.                                                             */
double *i_get_swan_Fy(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_Fy_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the x Stokes velocity tracer to get from      */
/* SWAN.                                                             */
double *i_get_swan_ste1(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_ste1_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the y Stokes velocity tracer to get from      */
/* SWAN.                                                             */
double *i_get_swan_ste2(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_ste2_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wavenumber tracer to get from SWAN         */
double *i_get_swan_k(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_k_m(void* hmodel)
{
  return (NULL);
}
 /* Returns a poiner to the wave Bernoulli tracer to get from SWAN    */
double *i_get_swan_Kb(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_Kb_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave whitecapping tracer to get from SWAN    */
double *i_get_swan_fwcapx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fwcapx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fwcapy(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fwcapy_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave depth-induced breaking tracer to get from SWAN    */
double *i_get_swan_fbrex(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fbrex_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fbrey(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fbrey_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave bottom friction dissapation tracer to get from SWAN    */
double *i_get_swan_fbotx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fbotx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fboty(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fboty_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave surface streaming tracer to get from SWAN    */
double *i_get_swan_fsurx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fsurx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fsury(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_fsury_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave form drag tracer to get from SWAN    */
double *i_get_swan_wfdx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wfdx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wfdy(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wfdy_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave ocean viscous stress tracer to get from SWAN    */
double *i_get_swan_wovsx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wovsx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wovsy(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_wovsy_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the wave roller tracer to get from SWAN    */
double *i_get_swan_frolx(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_frolx_m(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_froly(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_froly_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the dum1 tracer 1 to get from SWAN            */
double *i_get_swan_dum1(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_dum1_m(void* hmodel)
{
  return (NULL);
}

/* Returns a poiner to the dum2 tracer 2 to get from SWAN            */
double *i_get_swan_dum2(void* hmodel)
{
  return (NULL);
}
double *i_get_swan_dum2_m(void* hmodel)
{
  return (NULL);
}

/*-------------------------------------------------------------------*/
/* 2D (benthic) tracers                                              */
/*-------------------------------------------------------------------*/

/* Returns 1 if 2D tracer 'name' exists                              */
int i_tracername_exists_2d(void* model, char*name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn;

  for(tn=0; tn<wincon->ntrS; tn++) {
    /* Water column */
    if((tracer_find_index(name, wincon->ntrS,
              wincon->trinfo_2d) >= 0))
      return 1;
  }
  return 0;
}
/* Wrapper to 2D tracer name exists                                  */
int ginterface_tracername_exists_epi(void* model, char*name)
{
  return i_tracername_exists_2d(model, name);
}

/* Gets the number of explicit and auto 2D tracers                   */
int i_get_num_tracers_2d(void* hmodel, int *atr)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    *atr = master->atrS;
    return wincon->ntrS;
}
/* Wrapper for number of 2D tracers                                  */
int ginterface_getnumberofBtracer(void* hmodel)
{
    int atr;
    return i_get_num_tracers_2d(hmodel, &atr);
}

/* Returns the number of 2D tracers                                  */
int ginterface_getnepis(void* model)
{
    return 0;
}

/* Gets an array of 2D tracer names                                  */
void i_get_names_tracers_2d(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->ntrS; n++) {
      strcpy(trname[n], wincon->trinfo_2d[n].name);
    }
}

/* Gets an 2D tracerstat name                                        */
void i_get_tracerstats_2d(void* hmodel, char *trname, int tn) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    strcpy(trname, wincon->trinfo_2d[tn].tracerstat);
}

/* Returns the name of a 2D tracer i                                 */
char* ginterface_get2Dtracername(void* model, int i)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  if(i >= 0 && i < wincon->ntrS)
    return wincon->trinfo_2d[i].name;

  return NULL;
}
/* Wrapper for tracer name                                           */
char* ginterface_getepiname(void* model, int i)
{
  return ginterface_get2Dtracername(model,i);
  /*  return NULL; */
}
void ginterface_getnameofBtracer(void* hmodel, int n, char *tracername)
{
    strcpy(tracername, ginterface_get2Dtracername(hmodel, n));
}

/* Gets the index in the tracer info list of 2D tracer n             */
int i_get_param_map_2d(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    if(wincon->trinfo_2d != NULL && wincon->ntrS > 0)
      return wincon->trinfo_2d[n].m;
    else
      return 0;
}

/* Returns the value of tracer n at coordinate c                     */
double ginterface_getvalueofBtracer(void* hmodel, int n, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2=window->m2d[c];
    // fprintf(stderr,"sediments:tarcer number %d, val %lf \n", n, windat->tr_wcS[n][c2] );
  return windat->tr_wcS[n][c2];
    //    return windat->tr_wcS[n][c];
}

/* Returns a list of 2D tracer indexes corresponding to tracer names */
/* trname[].                                                         */
int *i_get_tmap_2d(void* model, int ntr, char *trname[])
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn;
  int *tmap_2d = i_alloc_1d(ntr);
  int n;
  for(tn=0; tn<ntr; tn++) {
    if(strcmp(trname[tn], "eta") == 0) {
      tmap_2d[tn] = ETA;
    }
    else {
      if((tmap_2d[tn] = tracer_find_index(trname[tn], wincon->ntrS,
            wincon->trinfo_2d)) < 0)
  hd_warn("ginterface : Can't find type WC2D tracer '%s' in parameter file.\n", trname[tn]);
    }
  }
  return(tmap_2d);
}

/* Returns the number of 2D column tracers whos name is prefixed by  */
/* 'R_'.                                                             */
int ginterface_get_num_rsr_tracers(void* model)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn, n;

  n = 0;
  for(tn = 0; tn < wincon->ntrS; tn++) {
    char *trname = wincon->trinfo_2d[tn].name;
    if (strncmp(trname, "R_", 2) == 0)
      n++;
  }
  return(n);
}

/* Returns an array of 1D pointers to tracers whos name is prefixed  */
/* by 'R_'.                                                          */
void ginterface_get_rsr_tracers(void* model, int *rtns)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn, n;

  n = 0;
  for(tn = 0; tn < wincon->ntrS; tn++) {
    char *trname = wincon->trinfo_2d[tn].name;
    if (strncmp(trname, "R_", 2) == 0)
      rtns[n++] = tn;
  }
}

/* Returns the number of 3D column tracers whos name is prefixed by  */
 /* 'Ed_'.                                                             */
 int ginterface_get_num_ed_tracers(void* model)
 {
   geometry_t *window = (geometry_t *)model;
   win_priv_t *wincon = window->wincon;
   int tn, n;
 
   n = 0;
   for(tn = 0; tn < wincon->ntr; tn++) {
     char *trname = wincon->trinfo_3d[tn].name;
     if (strncmp(trname, "Ed_", 3) == 0)
       n++;
   }
   return(n);
 }
 
 /* Returns an array of 1D pointers to tracers whos name is prefixed  */
 /* by 'Ed_'.                                                          */
 void ginterface_get_ed_tracers(void* model, int *rtns)
 {
   geometry_t *window = (geometry_t *)model;
   win_priv_t *wincon = window->wincon;
   int tn, n;
 
   n = 0;
   for(tn = 0; tn < wincon->ntr; tn++) {
     char *trname = wincon->trinfo_3d[tn].name;
     if (strncmp(trname, "Ed_", 3) == 0)
       rtns[n++] = tn;
   }
 }


/* Returns the diagnostic flag of 2D tracer 'name'                   */
int ginterface_getepidiagnflag(void* model, char* name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int index = tracer_find_index(name, wincon->ntrS, wincon->trinfo_2d);
  return wincon->trinfo_2d[index].diagn;
}

/* Returns the tracer info of 2D tracers                             */
tracer_info_t* ginterface_getepiinfo(void* model, int* n)
{
    return NULL;
}

/* Gets 1D pointers to an array of 2D tracers at coordinate c.        */
void i_get_tracer_2d(void* hmodel, int c, int ntr, int *tmap, double **tr)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int n, m;
    for(n = 0; n < ntr; n++) {
      m = tmap[n];
      if(m == NOT)
        continue;
      else if(m == ETA)
        tr[n] = &windat->eta[c2];
      else
        tr[n] = &windat->tr_wcS[m][c2];
    }
}

/* Returns 1D pointers to an array of 2D tracers at index b.         */
double **ginterface_getepivars(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nepi = e_nepi;
  double **epivar = (double **)calloc(nepi, sizeof(double*));
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];
  int n;

  for (n = 0; n < nepi; ++n) {
    epivar[n] = &windat->tr_wcS[epi_map[n]][c2];
  }
  return epivar;
}

/* Returns a pointer to tracer 'name' at coordinate c                */
double *ginterface_getpointerBtracer(void* hmodel, int n, int c)
{
    double *a;
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    a = &windat->tr_wcS[n][c2];
    return a;
}

/*-------------------------------------------------------------------*/
/* Sediment tracers                                                  */
/*-------------------------------------------------------------------*/

/* Returns 1 if sediment tracer 'name' exists                        */
int i_tracername_exists_sed(void* model, char*name)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn;

  for(tn=0; tn<wincon->nsed; tn++) {
    /* SEDIMENT */
    if((tracer_find_index(name, wincon->nsed,
              wincon->trinfo_sed) >= 0))
      return 1;
  }
  return 0;
}

/* Returns the number of 3d tracers that have a sediment component   */
int ginterface_getnumberoftracers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon=window->wincon;
    int ns,n,m;
    char name_wc[MAXSTRLEN];
    char name_sed[MAXSTRLEN];
    ns = 0;
    for(n = 0; n < windat->ntr; n++) {
      strcpy(name_wc, wincon->trinfo_3d[n].name);
      for(m = 0; m < windat->nsed; m++) {
	strcpy(name_sed, wincon->trinfo_sed[m].name);
	if(strcmp(name_sed,name_wc) == 0) {
	  ns = ns + 1;
	}
      }
    }
    return ns;
}

/* Gets the number of 3d tracers that have a sediment component      */
int gi_gettracernames(void* hmodel, char **tracername)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon=window->wincon;
    int ns,n,m;
    char name_wc[MAXSTRLEN];
    char name_sed[MAXSTRLEN];
    ns=0;
    for(n = 0; n < windat->ntr; n++) {
	strcpy(name_wc, wincon->trinfo_3d[n].name);
	for(m = 0; m < windat->nsed; m++) {
	    strcpy(name_sed, wincon->trinfo_sed[m].name);
	    if(strcmp(name_sed,name_wc) == 0) {
		strcpy(tracername[ns], name_wc);
		ns = ns + 1;
	    }
	}
    }
    return ns;
}

/* Gets the number of sediment tracers                               */
int i_get_num_sedtracers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    return wincon->nsed;
}

/* Gets an array of sediment tracer names                            */
void i_get_names_tracers_sed(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->nsed; n++) {
      strcpy(trname[n], wincon->trinfo_sed[n].name);
    }
}

/* Gets the index in the tracer info list of sediment tracer n       */
int i_get_param_map_sed(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    if(wincon->trinfo_sed != NULL && wincon->nsed > 0)
      return wincon->trinfo_sed[n].m;
    else
      return 0;
}

/* Returns the number of sediment tracer with name 'name'            */
/* Can also use tracer_find_index().                                 */
int gi_getmap2hdsedtracer(void *hmodel, int ns, char *name)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    win_priv_t *wincon=window->wincon;
    int n;
    /*int index = tracer_find_index(name, windat->nsed, wincon->trinfo_sed);*/
    for(n = 0; n < windat->nsed; n++)
	if(strcmp(name, wincon->trinfo_sed[n].name )==0)
		 return n;
    return (-1);
}

/* Returns a list of sediment tracer indexes corresponding to tracer */
/* names trname[].                                                   */
int *i_get_tmap_sed(void* model, int ntr, char *trname[])
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  int tn,n;
  int *tmap_sed = i_alloc_1d(ntr);
  
  for(tn=0; tn<ntr; tn++) {
    n = tracer_find_index(trname[tn], wincon->nsed, wincon->trinfo_sed);
    if(n < 0) {
      hd_warn("ginterface : Can't find type SED tracer '%s' in parameter file at %d.\n", trname[tn],tn);
    }
    tmap_sed[tn] = n;
  }

  return(tmap_sed);
}

/* Returns the fill value for sediment tracer 'name'                 */
double ginterface_get_fillvalue_sed(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    
    double v=0;
    v = tr->fill_value_sed;
    return v;
}

/* Gets 2D pointers to an array of sediment tracers at coordinate c  */
void i_get_tracer_sed(void* hmodel, int c,  int ntr, int *tmap, double ***tr_sed)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  int bot_k = 0;
  int top_k = window->sednz-1;
  int k, n, m;
  
  for(k = bot_k; k <= top_k; k++) {
    for(n = 0; n < ntr; n++) {
      if((m = tmap[n]) != -1)
	tr_sed[n][k] = &windat->tr_sed[m][k][c2];
    }
  }
}
/* Returns 1D pointers to an array of sediment tracers at index b    */
double** ginterface_getsedtracers(void* model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int nz = window->sednz;
  int ntr = e_ntr;
  double **sedtr = (double **)calloc(ntr * nz, sizeof(double*));
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];
  int k, n;

  for (k = 0; k < nz; ++k)
    for (n = 0; n < ntr; ++n) {
      sedtr[k * ntr + n] = &windat->tr_sed[sed_map[n]][k][c2];
    }
  return sedtr;
}

double ginterface_get_svel(void* model, char *name)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    double v=0;
    if (data)
      v = data->svel;
    return v;
}

void ginterface_get_svel_name(void* model, char *name, char *sname)
{
    geometry_t* window = (geometry_t*) model;
    win_priv_t *wincon=window->wincon;
    tracer_info_t *tr = i_get_tracer(model, name);
    trinfo_priv_sed_t *data = tr->private_data[TR_PRV_DATA_SED];
    
    if (data)
      strcpy(sname, data->svel_name);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

/* Returns the model time                                            */
double i_get_model_time(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->t;
}
/* Wrapper for model time                                            */
double ginterface_getmodeltime(void* model)
{
  return i_get_model_time(model);
}

/* Returns the master timeunit                                       */
char *i_get_model_timeunit(void* hmodel)
{
  return(master->timeunit);
}

/* Returns the window timeunit                                       */
char *ginterface_gettimeunits(void *model)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  return wincon->timeunit;
}

/* Returns the output timeunit                                       */
char *ginterface_getoutputtimeunits(void *model)
{
  return master->output_tunit;
}

/* Returns the model 3D time-step                                    */
double i_get_model_timestep(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->dt;
}
/* Wrapper for 3D time-step                                          */
double ginterface_getmodeltimestep(void* hmodel)
{
    return i_get_model_timestep(hmodel);
}

/* Returns 1 if the driver is the transport model                    */
int ginterface_transport_mode(void)
{
  return(master->runmode & TRANS);
}

char * ginterface_get_output_path(void)
{
  char *opath = master->opath;
  if (strlen(opath))
    return(opath);
  else
    return(NULL);
}

/* Returns the number of processors for transport                    */
#if defined(HAVE_OMP)
int ginterface_get_trans_num_omp(void *model)
{
  geometry_t* window = (geometry_t*) model;
  win_priv_t *wincon=window->wincon;
  return(wincon->trans_num_omp);
}
#endif

/* Returns the window number                                         */
int ginterface_get_win_num(void *model)
{
  return(((geometry_t *)model)->wn);
}

/* Sets the error function to warning                                */
void gi_set_errfn_warn(void *hmodel)
{
 prm_set_errfn(warn);
}

/* Returns the number of vertical layers                             */
int  i_get_num_wclayers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->nz;
}
/* Wrapper for model layers                                          */
int  ginterface_getnumberwclayers(void* model)
{
  return i_get_num_wclayers(model);
}

/* Returns the size of the 2D vector                                 */
int i_get_winsize(void *hmodel)
{
  geometry_t* window = (geometry_t*) hmodel;
  return (window->sgsizS);
}

/* Returns the number of 2D columns                                  */
int  i_get_num_columns(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->b2_t;
}
/* Wrapper for 2D columns                                            */
int  ginterface_get_max_numbercolumns(void* model)
{
  return i_get_num_columns(model);
}
int  ginterface_getnumbercolumns(void* hmodel) 
{
  return i_get_num_columns(hmodel);
}

/* Returns the number of columns in the cells-to-process vector (wet */
/* cells only; i.e. excluding OBC cells). Note: this has to be       */
/* called after setup is complete. Use the max_columns function if   */
/* you need to pre-allocate stuff.                                   */
int ginterface_getnumbercolumns_e(void* model)
{
  geometry_t *window = (geometry_t *)model;
  process_cell_mask(window, window->cbgc, window->ncbgc);
  return window->wincon->vcs2;
}

/* Returns 1 if cell index b is a boundary cell                      */
int ginterface_isboundarycolumn(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, b);
  
  /* Exclude process points if required */
  if (window->wincon->c2[window->m2d[c]]) return 1;

  /* Do ecology on boundary cells */
  // FR (11/09) : By definition now that we've changed over to vca2,
  //              we wont ever have boundary or dry columns
  return 0;
  /* Don't do ecology on boundary cells */
  /*return b + 1 > window->v2_t ? 1 : 0;*/
}

/* Returns the model simulation step                                 */
int i_get_nstep(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->nstep;
}

/* Returns the number of sediment layers                             */
int  i_get_num_sedlayers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->sednz;
}
/* Wrapper to sediment layers                                        */
int  ginterface_getnumbersedlayers(void* model)
{
  return i_get_num_sedlayers(model);
}

/* Returns the coordinate corresponding to the index in the          */
/* cells-to-process list.                                            */
int i_get_c(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  return window->wincon->s2[b+1];
}

/* Returns the counter index (referenced to zero) corresponding to   */
/* the coordinate c.                                                 */
int i_get_interface_counter(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int sc = cc - 1;
    return sc;
}

/* Returns the counter index referenced to one.                      */
int i_get_host_counter(void* hmodel, int c)
{
    return c+1;
}

/* Returns the k layer of the surface for coordinate c.              */
int i_get_topk_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int v = window->s2k[window->nsur_t[cc]];
    return v;
}
/* Wrapper for surface layer                                         */
int ginterface_gettopk_wc(void* hmodel, int c)
{
  return i_get_topk_wc(hmodel, c);
}
/* Returns the k layer of the surface for cells-to-process index b.  */
int ginterface_getwctopk(void *model, int b)
{
  int c = i_get_c(model, b);
  return i_get_topk_wc(model, c);
}

/* Returns the k layer of the bottom for coordinate c.               */
int i_get_botk_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int v = window->s2k[window->bot_t[cc]];
    return v;
}
/* Wrapper for bottom layer                                          */
/* If in hd model cell numbering increases downward, return the      */
/* number of the top hd cell.                                        */
int ginterface_getbotk_wc(void* hmodel, int c)
{
  return i_get_botk_wc(hmodel, c);
}
/* Returns the k layer of the bottom for cells-to-process index b.   */
int ginterface_getwcbotk(void *model, int b)
{
  int c = i_get_c(model, b);
  return i_get_botk_wc(model, c);
}

/* Returns the k layer of the sediment surface. c not used.          */
int i_get_topk_sed(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->sednz;
}
/* Returns the k layer of the sediment surface. b not used.          */
int ginterface_getsedtopk(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  if (window != NULL)
    return i_get_topk_sed(model, b) - 1;
  else {
    // This happens for ecology_pre_build
    return 0;
  }
}

/* Returns the k layer of the sediment bottom. c not used.           */
int i_get_botk_sed(void* hmodel, int c)
{
  return 0;
}
/* Returns the k layer of the sediment bottom. b not used.           */
int ginterface_getsedbotk(void *model, int b)
{
  return i_get_botk_sed(model, b);
}

/* Returns cell thickness in column c from surface to bottom. The dz */
/* array is supplied to the function. The index c is that for the    */
/* full 3D array.                                                    */
/* Note:                                                             */
/* k index changes from botk_wc to topk_wc, where botk_wc and        */
/* topk_wc are numbers of the bottom and top cells respectively,     */
/* botk_wc < topk_wc <= (nwclayers-1).                               */
/* dz_wc[n][botk_wc] is thickness of the bottom cell                 */
/* dz_wc[n][topk_wc] is thickness of the top cell                    */
/* If in the hd model cell numbering increases downward              */
/* (topk_wc_hd < botk_wc_hd), make sure that in dz_wc the cell       */
/* numbering increases upward, starting from topk_wc_hd in the       */
/* bottom cell and to  botk_wc_hd in the top cell.                   */
void i_get_dz_wc(void* hmodel, int c, double *dz_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;

    for(k = top_k; k >= bot_k; k--) {
      dz_wc[k] = wincon->dz[c3] * wincon->Hn1[c2];
      c3 = window->zm1[c3];
    }
}
/* Wrapper for cell thickness                                        */
void ginterface_getdz_wc(void* hmodel, int c, double *dz_wc)
{
  i_get_dz_wc(hmodel, c, dz_wc);
}
/* Returns cell thickness in column b from surface to bottom. The dz */
/* array is returned from the function. The index b is that for the  */
/* reduced cells-to-process vector.                                  */
double *ginterface_getwccellthicknesses(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, b);
  double *dz = calloc(window->nz, sizeof(double));
  i_get_dz_wc(model, c, dz);
  return dz;
}

/* Gets the sediment layer depths                                    */
void i_get_gridz_sed(void* hmodel, int c, double *gridz_sed)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int bot_k = 0;
    int top_k = window->sednz-1;
    int k;
    for(k = bot_k; k <= top_k+1; k++)
      gridz_sed[k] = window->gridz_sed[k][c2];
}
/* Wrapper for sediment layer depths                                 */
void ginterface_getgridz_sed(void* hmodel, int c, double *gridz_sed)
{
  return i_get_gridz_sed(hmodel, c, gridz_sed);
}

/* Returns sediment cell thickness in sediment column c. The dz      */
/* array is supplied to the function. The index c is that for the    */
/* full 3D array.                                                    */
void i_get_dz_sed(void* model, int c, double *dz_sed)
{
  geometry_t *window = (geometry_t *)model;
  int c2 = window->m2d[c];
  int nz = window->sednz;
  int k;
  if(nz <= 0)
  {
    emstag(LPANIC,"hd:ginterface:i_get_dz_sed","Requesting sdiment cell thickness without sediment layers,exiting ...");
    exit(0);
  }
  for (k = 0; k < nz; ++k) {
    dz_sed[k] = window->gridz_sed[k + 1][c2] - window->gridz_sed[k][c2];
  }
}
/* Returns sediment cell thickness in sediment column b. The dz      */
/* array is supplied to the function. The index b is that for the    */
/* reduced cells-to-process vector.                                  */
double *ginterface_getsedcellthicknesses(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, b);
  int nz = window->sednz;
  double *dz = calloc(nz, sizeof(double));

  i_get_dz_sed(model, c, dz);
  return dz;
}

/* Returns the cell area of the host grid at coordinate c            */
double i_get_cellarea_w(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    double v = window->cellarea[c2];
    return v;
}
/* Wrapper for cell area                                             */
double ginterface_getcellarea(void* hmodel, int c)
{
  return i_get_cellarea_w(hmodel, c);
}
/* Returns the cell area of the host grid at index b                 */
double ginterface_cellarea(void* hmodel, int b)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c = i_get_c(hmodel, b);
    return i_get_cellarea_w(hmodel, c);
}

/* Gets the e1 layer thicknesses for the column at c                 */
void i_get_cellarea_e1(void* hmodel, int c, double *area)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int k;
  for(k = top_k; k >= bot_k; k--) {
    /*
     * Note that the dz used here is the dz at the end of the hydro
     * timestep which could be different to the one at the start which
     * was used to calculate u1flux3d. see vel3d for more details
     */
    area[k] = windat->dzu1[c3] * window->h2au1[c2];
    c3 = window->zm1[c3];
  }
}

/* Gets the e2 layer thicknesses for the column at c                 */
void i_get_cellarea_e2(void* hmodel, int c, double *area)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int k;
  for(k = top_k; k >= bot_k; k--) {
    /*
     * Note that the dz used here is the dz at the end of the hydro
     * timestep which could be different to the one at the start which
     * was used to calculate u2flux3d. see vel3d for more details
     */
    area[k] = windat->dzu2[c3] * window->h1au2[c2];
    c3 = window->zm1[c3];
  }
}

/* Returns the depth referenced to msl at coordinate c               */
double i_get_botz_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon=window->wincon;
    int c2 = window->m2d[c];
    double v =  window->botz[c2] * wincon->Hn1[c2];
    return v;
}
/* Wrapper for depth.                                                */
/* If in hd model OZ is directed downward (ie. topz<botz) return (   */
/* -1)*botz.                                                         */
double ginterface_getbotz_wc(void* hmodel, int c)
{
  return  i_get_botz_wc(hmodel, c);
}

/* Returns the surface at coordinate c                               */
double i_get_topz_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    double v = windat->topz[c2];
    /* no need to multiply by depth as it
       is zero by definition if sigma grid is applied */
    return v;
}

/* Returns cell centred depth for index b, layer k in the            */
/* cells-to-process vector.                                          */
double ginterface_getcellz(void *model, int b, int k)
{
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int cs = window->nsur_t[cc];
  int cb = window->bot_t[cc];
  int kk = window->s2k[cb];
  
  /* Search for the correct level */
  for (c = cb; c != cs; c = window->zp1[c])
    if (k == kk++)
      break;
  
  return(window->cellz[c]);
}

/* Returns the depth referenced to the surface at coordinate c       */
double i_get_depth(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    double v = windat->eta[c2] - window->botz[c2];
    return v;
}

/* Returns the sea surface height                                    */
double i_get_eta(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    double v = windat->eta[c2];
    return v;
}
/* Wrapper for surface height.                                       */
/* If in hd model OZ is directed downward (ei. topz<botz) return     */
/* (-1)*topz.                                                        */
double ginterface_gettopz_wc(void* hmodel, int c)
{
  return i_get_eta(hmodel, c);
}

/* Gets the (i,j) index corresponding to index col. Only suitable    */
/* for structured grids.                                             */
void ginterface_get_ij(void* model, int col, int *ij)
{
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, col);
  int cs = window->m2d[c];

  ij[0] = window->s2i[cs];
  ij[1] = window->s2j[cs];

}

/* Returns the bottom stress magnitude for index b. Added for use by */
/* benthic plants.                                                   */
double ginterface_botstress(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];

  if (windat->tau_bm != NULL) 
    return windat->tau_bm[c2];
  else
    return 0.0;
}

/* Returns the wind speed magnitude for index b. Added to calculate  */
/* gas exchange.                                                     */
double ginterface_get_windspeed(void *model, int b)
{
   geometry_t *window = (geometry_t *)model;
   int c = i_get_c(model, b);
   int c2 = window->m2d[c];
   return window->windat->windspeed[c2];
}
/* Wrapper for windspeed                                             */
double ginterface_getwindspeed(void* hmodel, int c)
{
  return ginterface_get_windspeed(hmodel, c);
}

/* Returns the surface light for index b.                            */
double ginterface_getlighttop(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];

  if (windat->light != NULL) 
    return windat->light[c2];
  else
    return 0.0;
}

/* Returns the porosity in sediment layers for index b.              */
double *ginterface_getporosity(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  double por_def = 0.4;
  double *por = wincon->sd1;
  int nz = window->sednz;
  int ntr = windat->nsed;
  int c = i_get_c(model, b);
  int c2=window->m2d[c];
  int n, k;
  for(k = 0; k < nz; k++)
    por[k] = por_def;

#if defined(HAVE_SEDIMENT_MODULE)
  if (wincon->do_sed) {
    for(n = 0; n < ntr; n++) {
      // FR : This is to make sure we are backwards compatible
      if( (strcmp(wincon->trinfo_sed[n].name,"por_sed") == 0) ||
	  (strcmp(wincon->trinfo_sed[n].name,"porosity") == 0)) {
        for(k = 0; k < nz; k++)
          por[k] = windat->tr_sed[n][k][c2];
	break; // end for loop
      }
    }
  }
#endif

  return por;
}

/* Returns the erosion rate in sediment layers for index b.         */
double ginterface_geterosionrate(void *model, int b) {
}

/* Returns the bottom wave stress for index b.                       */
double ginterface_getustrcw(void *model, int b)
{
  geometry_t *window = (geometry_t *)model;
  window_t *windat = window->windat;
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];
  return windat->ustrcw[c2];
}

/* Calculates the Zenith at index b using the library function       */
double ginterface_calc_zenith(void *model, double t, int b)
{
  double lat;
  double elev;
  geometry_t *window = (geometry_t *)model;
  int c = i_get_c(model, b);
  int c2 = window->m2d[c];
  double ang = 7.29e-5;    /* Earth's angular velocity (s-1)         */
  char *tunit = window->wincon->timeunit;
  char *ounit = master->params->output_tunit;

  /* Latitude from Coriolis                                          */
  lat = asin(window->wincon->coriolis[c2] / (2.0 * ang));

  /* Call the library function to calculate the solar elevation      */
  if (window->is_geog)
    elev = calc_solar_elevation(NULL, tunit, t, lat, NULL, &window->cellx[c2]);
  else
    elev = calc_solar_elevation(ounit, tunit, t, lat, NULL, NULL);

  /* Zenith                                                          */
  return ( (PI/2.0) - elev);
}

/* Returns the angle of the e1 normal vector                         */
double i_get_thetau1(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  double v = window->thetau1[c];
  return v;
}

/* Returns the angle of the e2 normal vector                         */
double i_get_thetau2(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  double v = window->thetau2[c];
  return v;
}

/* Gets an array of anges between e1 direction of the grid and W-E   */
/* direction. Note: c index changes from 0 to ncol-1, the angle      */
/* units are rad and e1 direction of the grid is cell centered.      */
void ginterface_gettheta(void* hmodel, double *theta, int ncol)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon=window->wincon;
    int c, cc;
    if (window->b2_t > ncol)
	hd_quit("Wrong number of columns in ginterface_gettheta");
    for(cc=1; cc<=window->b2_t; cc++) {
	c=window->w2_t[cc];
	theta[cc-1] = wincon->d1[c];
    }
}

/* Returns cell centered sin() of the angle of the e1 normal vector. */
/* Structured meshes only.                                           */
double i_get_sinthcell(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  int xp1 = window->xp1[c];
  int yp1 = window->xp1[c];
  double v = (sin(window->thetau1[c]) +
	      sin(window->thetau1[xp1]) +
	      sin(window->thetau2[c]) +
	      sin(window->thetau2[yp1])) / 4.0;
  return v;
}


/* Returns cell centered cos() of the angle of the e2 normal vector. */
/* Structured meshes only.                                           */
double i_get_costhcell(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  int xp1 = window->xp1[c];
  int yp1 = window->xp1[c];
  double v = (cos(window->thetau1[c]) +
	      cos(window->thetau1[xp1]) +
	      cos(window->thetau2[c]) +
	      cos(window->thetau2[yp1])) / 4.0;
  return v;
}

/* Reads initial thicknesses of the sediment grid layers from the    */
/* parameter file.                                                   */
/* prmfd is the pointer to the parameter file                        */
/* sednzc is number of the sed layers                                */
/* returns sednz - number of sediment layers                         */
/* updates dz_sed[k] - a 1d array of sediment thicknesses            */
/* Note:                                                             */
/* k index changes from 0 to (sednz-1)                               */
/* dz_sed[0] is thickness of the bottom layer                        */
/* dz_sed[sednz-1] is the thickness of the top layer.                */
int ginterface_getdzinit_sed(FILE* prmfd, double *dz_sed, int sednzc)
{
    int k;
    int sednz=0;
    double *tmp;
    prm_read_darray(prmfd, "NSEDLAYERS", &dz_sed, &sednz);
    /* change order of numbering in dz_sed[k] */
    tmp = d_alloc_1d(sednz);
    for (k = 0; k < sednz; k++)
      tmp[k] = dz_sed[k];
    for (k = 0; k < sednz; k++)
      dz_sed[k] = tmp[sednz-1-k];
    d_free_1d(tmp);
    return sednz;
}


/* Returns the verbosity flag                                        */
int ginterface_getverbosity(void *model)
{
  extern int debug;

  return debug;
}

/* Gets a column of e1 volume fluxes at coordinate c                 */
void i_get_fluxe1_wc(void* hmodel, int c, double *u1flux3d)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int bot_k = window->s2k[window->bot_t[cc]];
  int top_k = window->s2k[window->nsur_t[cc]];
  int c3 = window->nsur_t[cc];
  int k;
  for(k = top_k; k >= bot_k; k--) {
    u1flux3d[k] = windat->u1flux3d[c3];
    c3 = window->zm1[c3];
  }
}

/* Gets a column of e2 volume fluxes at coordinate c                 */
void i_get_fluxe2_wc(void* hmodel, int c, double *u2flux3d)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      u2flux3d[k] = windat->u2flux3d[c3];
      c3 = window->zm1[c3];
    }
}

/* Gets a column of e1 velocities at coordinate c                    */
void i_get_u1_wc(void* hmodel, int c, double *u1)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      u1[k] = windat->u1[c3];
      c3 = window->zm1[c3];
    }
}
/* Get the water column horizontal velocities at coordinate c        */
/* updates u1_wc[k] - u1 component of water velocity.                */
/* Note:                                                             */
/* k index changes from botk_wc to topk_wc, where botk_wc and        */
/* topk_wc are numbers of the bottom and top cells, respectively,    */
/* and botk_wc < topk_wc <= (nwclayers-1).                           */
/* u1_wc[n][botk_wc] is velocity in the bottom cell                  */
/* u1_wc[n][topk_wc] is velocity in the top cell                     */
/* If in the hd model cell numbering increases downward              */
/* (topk_wc_hd < botk_wc_hd), make sure that in u1_wc the cell       */
/* numbering increases upward, starting from topk_wc_hd in the       */
/* bottom cell and to  botk_wc_hd in the top cell.                   */
/* u1 velocity must be copied to dummy array wincon->w4.             */
void ginterface_getu1_wc(void* hmodel, int c, double *u1_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon=window->wincon;
    int c2=window->m2d[c];
    int cc=window->c2cc[c2];
    int bot_k= window->s2k[window->bot_t[cc]];
    int top_k= window->s2k[window->nsur_t[cc]];
    int c3=window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
    	u1_wc[k] = wincon->w4[c3];
	//tmp 2013
	//	if(fabs(u1_wc[k]) > 1e-9) fprintf(stderr, "ginterface_get_u1: u1=%lf \n", u1_wc[k]); 
    	/*u1_wc[k] = windat->u1[c3];*/
    	c3 = window->zm1[c3];
    }
}

/* Gets a column of e2 velocities at coordinate c                    */
void i_get_u2_wc(void* hmodel, int c, double *u2)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      u2[k] = windat->u2[c3];
      c3 = window->zm1[c3];
    }
}

/* u2 velocity must be copied to dummy array wincon->w5.             */
void ginterface_getu2_wc(void* hmodel, int c, double *u2_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon=window->wincon;
    int c2=window->m2d[c];
    int cc=window->c2cc[c2];
    int bot_k= window->s2k[window->bot_t[cc]];
    int top_k= window->s2k[window->nsur_t[cc]];
    int c3=window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
    	u2_wc[k] = wincon->w5[c3];
    	/*u2_wc[k] = windat->u2[c3];*/
    	c3 = window->zm1[c3];
    }
}

/* Gets a column of vertical velocities at coordinate c              */
void i_get_w_wc(void* hmodel, int c, double *w)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      w[k] = windat->w[c3];
      c3 = window->zm1[c3];
    }
}

/* Read water column vertical diffusion coefficients at coordinate   */
/* c. Updates Kz_wc[k] - vertical diffusion coefficient in water.    */
/* Note:                                                             */
/* k index changes from botk_wc to topk_wc+1, where botk_wc and      */
/* topk_wc are numbers of the bottom and top cells, respectively,    */
/* and botk_wc < topk_wc <= (nwclayers-1).                           */
/* Kz_wc[n][botk_wc] is Kz in the bottom cell                        */
/* Kz_wc[n][topk_wc] is Kz in the top cell                           */
/* If in the hd model cell numbering increases downward              */
/* (topk_wc_hd < botk_wc_hd), make sure that in Kz_wc the cell       */
/* numbering increases upward, starting from topk_wc_hd in the       */
/* bottom cell and to  (botk_wc_hd+1) in the top cell.               */
void i_get_Kz_wc(void* hmodel, int c, double *Kz_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      Kz_wc[k] = windat->Kz[c3];
      c3 = window->zm1[c3];
    }
    Kz_wc[top_k+1] = 0.;
}
/* Wrapper for vertical diffusion                                    */
void ginterface_getkz_wc(void* hmodel, int c, double *Kz_wc)
{
  i_get_Kz_wc(hmodel, c, Kz_wc);
}

/* Gets a column of vertical viscosities at coordinate c             */
void i_get_Vz_wc(void* hmodel, int c, double *Vz_wc)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int bot_k = window->s2k[window->bot_t[cc]];
    int top_k = window->s2k[window->nsur_t[cc]];
    int c3 = window->nsur_t[cc];
    int k;
    for(k = top_k; k >= bot_k; k--) {
      Vz_wc[k] = windat->Vz[c3];
      c3 = window->zm1[c3];
    }
}

/** Verify if the current sparse index represents a  valid cartesian coordinate
 * given by the two int arrays i and j. Note that the cartesian coordinates must
 * match by index as well (must be at the position in the respective array).
 *
 * @param hmodel - the reference to the host model
 * @param col - the sparse index
 * @param i - i coordiante of the cartesian map
 * @param j - j coordinate of the cartesian map
 * @param nn - length of i and j
 * @return the position in i and j if found, -1 otherwise
 */
int i_in_window(void* hmodel, int col, int* i, int* j,int nn)
{
  int ni;
  int ci,cj;
  geometry_t* window = (geometry_t*) hmodel;
  ci = window->s2i[col];
  cj = window->s2j[col];
  for(ni = 0 ; ni < nn ; ni++)
  {
    if(i[ni] == ci && j[ni] == cj)
      return ni;
  }
  return -1;
}

/* Returns the tracerstat step given a text input                    */
int i_get_step(void* model, char* st)
{
    if(strcmp(st,"PRE_DECAY") == 0)
        return 1;
    else if(strcmp(st,"POST_DIAG") == 0)
        return 1;
    else if(strcmp(st,"PRE_WAVE") == 0)
        return 2;
    else if(strcmp(st,"POST_DECAY") == 0)
        return 2;
    else if(strcmp(st,"PRE_SED") == 0)
        return 6;
    else if(strcmp(st,"POST_WAVE") == 0)
        return 6;
    else if(strcmp(st,"POST_SED") == 0)
        return 3;
    else if(strcmp(st,"PRE_ECO") == 0)
        return 3;
    else if(strcmp(st,"POST_ECO") == 0)
        return 4;
    else if(strcmp(st,"PRE_RTSTAT") == 0)
        return 4;
    else if(strcmp(st, "POST_RTSTAT") == 0)
        return 5;

    return 0;
}

/** Retrieve the cartesian coordiantes ij and depth k for sparse index col
 *
 * @param hmodel - the reference to the host model
 * @param col- the sparse index
 * @param ij - a valid int array of length 3
 * @return ij
 *
 */
int* i_get_ijk(void* hmodel, int col, int* ij)
{
  geometry_t* window = (geometry_t*) hmodel;
  ij[0] = window->s2i[col];
  ij[1] = window->s2j[col];
  ij[2] = window->s2k[col];
  emstag(LTRACE,"hd:ginterface","got ijk %u %u %u ",window->s2i[col],window->s2j[col],window->s2k[col]);

  return ij;
}

/* Returns the i Cartesian coordinate corresponding to coordinate c. */
/* For structured grids only.                                        */
int ginterface_getcelli(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->s2i[c];
}

/* Returns the j Cartesian coordinate corresponding to coordinate c. */
/* For structured grids only.                                        */
int ginterface_getcellj(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->s2j[c];
}

/** Test the cartesian coordiantes ij and depth k to match sparse index col
 *
 * @param hmodel - the reference to the host model
 * @param col- the sparse index
 * @param ij - a valid int array of length 3 or 2
 * @param do_k if k should be tested (that requires ij to be of length 3)
 * @return 1 ifthe cartesian and sparse inde match, 0 otherwise
 *
 */
int i_is_ijk(void* hmodel, int col, int* ij, int do_k)
{
  geometry_t* window = (geometry_t*) hmodel;
  emstag(LTRACE,"hd:ginterface","testing c %u [ij %u %u] against %u %u",col,window->s2i[col],window->s2j[col],ij[0],ij[1]);

  if(ij[0] == window->s2i[col] && ij[1] == window->s2j[col])
  {
  	if(do_k)
  	{
  		if(ij[2] == window->s2k[col])
  			return 1;
  	  else
  	    return 0;
  	}else
  		return 1;
  }
  return 0;
}

/* Insert text into the error buffer                                 */
void i_set_error(void* hmodel, int col, int errorf, char *text)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  wincon->gint_errorf[col] = errorf;
  sprintf(wincon->gint_error[col], "%c", '\0');
  if (text != NULL) // In case of a reset
    strcpy(wincon->gint_error[col], text);
}

/* Returns the error flag for coordinate col                         */
int i_get_error(void* hmodel, int col)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  return(wincon->gint_errorf[col]);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Tracerstat functions                                              */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

/*UR in lack of a better place to put it ...                         */
/* custom_function add and remove                                    */
custom_function_t* custom_stack_remove(char* fnc)
{
  /* this function uses the global structs master */
  custom_function_t* custom_fnc1;
  custom_function_t* custom_fnc2 = NULL;
  if(fnc == NULL)
  {
    emstag(LWARN,"hd:ginterface:custom_stack_add","Cannot remove NULL to function stack!");
    return NULL;
  }

  if(master->custom_fnstack == NULL || master->custom_fnstack->nfunctions == 0)
    return NULL;

  custom_fnc1 = hd_data->master->custom_fnstack->functions;
  while(custom_fnc1 != NULL)
  {
    if(strcmp(custom_fnc1->name,fnc) == 0)
    {
      if(custom_fnc2 == NULL)
        master->custom_fnstack->functions = custom_fnc1->next;
      else
        custom_fnc2->next = custom_fnc1->next;
      master->custom_fnstack->nfunctions--;
      emstag(LTRACE,"hd:ginterface:custom_stack_add","removed custom function '%s' from function stack!",fnc);

      return custom_fnc1;
    }
    custom_fnc2 = custom_fnc1;
    custom_fnc1 = custom_fnc1->next;
  }
  return NULL;
}


char* custom_stack_add(custom_function_t* fnc)
{
  /* this function uses the global struct master */
  custom_function_t* custom_fnc;
  if(fnc == NULL)
  {
    emstag(LWARN,"hd:ginterface:custom_stack_add","Cannot add NULL to function stack!");
    return NULL;
  }

  if(master->custom_fnstack == NULL)
  {
    emstag(LTRACE,"hd:ginterface:custom_stack_add","creating function stack!");
    master->custom_fnstack = (custom_stack_t*)malloc(sizeof(custom_stack_t));
    if(master->custom_fnstack == NULL)
      hd_quit("custom_stack_add - Failed to assign memory for custom_funtion stack!");
    master->custom_fnstack->nfunctions = 0;
    master->custom_fnstack->functions= NULL;
  }

  custom_fnc = master->custom_fnstack->functions;
  if(custom_fnc == NULL)
    master->custom_fnstack->functions = fnc;
  else
  {
    while(custom_fnc->next != NULL)
    {
      if(strcmp(custom_fnc->name,fnc->name) == 0)
      {
        emstag(LTRACE,"hd:ginterface:trstat_domain_function_init","Attempted to add the custom function '%s' twice! - Ignored at %d ",fnc->name,master->custom_fnstack->nfunctions);
        return NULL;
      }
      custom_fnc = custom_fnc->next;
    }
    if(strcmp(custom_fnc->name,fnc->name) == 0)
    {
      emstag(LWARN,"hd:ginterface:custom_stack_add","Attempted to add the custom function '%s' twice! - Ignored at %d ",fnc->name,master->custom_fnstack->nfunctions);
      return NULL;
    }
    custom_fnc->next = fnc;
  }
  emstag(LDEBUG,"hd:ginterface:custom_stack_add","adding custom function '%s' to function stack!",fnc->name);
  master->custom_fnstack->nfunctions++;
  /* fnc->init(hd_data,fnc); */
  return fnc->name;
}



#if defined(HAVE_TRACERSTATS_MODULE)
/* generic domain-functions to add to allow runtime collection of data */
#if defined(HAVE_MPI)
void trstat_domain_function_init(hd_data_t* hdata, custom_function_t* fnc)
{
}

void trstat_domain_function_gather(custom_function_t* func, hd_data_t* hdata)
{
}

void trstat_domain_function_scatter1(custom_function_t* func, hd_data_t* hdata)
{
}

void trstat_domain_function_scatter2(custom_function_t* func, hd_data_t* hdata)
{
}

#else

void trstat_domain_function_init(hd_data_t* hdata, custom_function_t* fnc)
{
  emstag(LTRACE,"hd:ginterface:trstat_domain_function_init","attempting Init Domain Function, as: %s ",fnc->name);
  
  fnc->data->data = fnc->create(((trs_t*)hdata->wincon[1]->trs)->domainfn3d[fnc->data->index].data);
  
}


void trstat_domain_function_gather(custom_function_t* func, hd_data_t* hdata)
{
  master_t* master = hdata->master;
  geometry_t** windows = hdata->window;
  int n = func->data->index;
  int i;
  for(i = 1; i <= master->nwindows;i++ )
    func->do_gather(n,func->data->data,windows[i]->wincon->trs);
}

/* a domain functionwhich calls do_scatter per window */
void trstat_domain_function_scatter1(custom_function_t* func, hd_data_t* hdata)
{
  master_t* master = hdata->master;
  geometry_t** windows = hdata->window;
  int n = func->data->index;

  int i;
  for(i = 1; i <= master->nwindows;i++ )
    func->do_scatter(n,func->data->data,windows[i]->wincon->trs);
}


/*a domain funtion which calls do_scatter only once with the master i.e.for i/o out put */
void trstat_domain_function_scatter2(custom_function_t* func, hd_data_t* hdata)
{
  int n = func->data->index;
  func->do_scatter(n,func->data->data,NULL);
}



void trstat_add_domain_function(char* name, int n,int scattertype, custom_create_data create , custom_function_do_gather gather, custom_function_do_scatter scatter)/* , custom_destroy_function destroy )*/
{
  custom_function_t* fnc = (custom_function_t*) malloc(sizeof(custom_function_t));
  emstag(LTRACE,"hd:ginterface:trstat_add_domain_function","Starting to create Domain Function, adding as: %s",name);
  if(fnc == NULL)
  {
    emstag(LFATAL,"hd:ginterface:trstat_add_domain_function","Failed to assign memory for custom function!- Exiting");
    exit(1);
  }
  fnc->data = (custom_data_t*) malloc(sizeof(custom_data_t));
  if(fnc->data == NULL)
  {
    emstag(LFATAL,"hd:ginterface:trstat_add_domain_function","Failed to assign memory for custom function data!- Exiting");
    exit(1);
  }
  fnc->do_gather = gather;
  fnc->do_scatter= scatter;
  fnc->data->index = n;
  fnc->gather = trstat_domain_function_gather;
  if(scattertype == 2)
    fnc->scatter = trstat_domain_function_scatter2;
  else
    fnc->scatter = trstat_domain_function_scatter1;
  fnc->init =  trstat_domain_function_init;
  fnc->create = create;
  fnc->name = (char*) malloc(sizeof(char) * (strlen(name)+1));
  strcpy(fnc->name,name);
  fnc->next = NULL;
  fnc->destroy = NULL;/*destroy; */
  emstag(LTRACE,"hd:ginterface:trstat_add_domain_function","Created Domain Function, adding as: %s",name);
  custom_stack_add(fnc);
}

void trstat_remove_domain_function(char* name)
{
  custom_function_t* fnc =  custom_stack_remove(name);
  if(fnc != NULL)
  {
    if(fnc->destroy != NULL)
      fnc->destroy(fnc->data);
    free(fnc);
  }
}
#endif
#endif

#if defined(HAVE_WAVE_MODULE)

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Wave functions                                                    */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Wave step                                                         */
/*-------------------------------------------------------------------*/
void wave_interface_step(geometry_t *window)
{
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  int cc, c;

  /*-----------------------------------------------------------------*/
  /* Do the ecology if required */
#if defined(HAVE_WAVE_MODULE)
  if ((wincon->do_wave) && !(windat->nstep % wincon->wave->dt_ratio)) {
    /*
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
    */
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = wincon->s2[cc];
      /* Exclude process points if required */
      if (window->ncwave)
	if (ANY(c, window->cwave, window->ncwave)) continue;
      wave_step(window, wincon->wave, c);
    }
  }
#endif

}

/* END wave_step()                                                   */
/*-------------------------------------------------------------------*/

/* Returns the wave forcing method                                   */
int w_do_waves(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  if (wincon->do_wave & W_FILE)
    return 1;
  if (wincon->do_wave & W_COMP)
    return 2;
  if (wincon->do_wave & W_SWAN)
    return 4;
  return 0;
}

/* Returns the wave time-step                                        */
double w_get_dt(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  double v = wincon->wavedt;
  return v;
}

/* Returns the bottom drag limit                                     */
double i_get_quad_bfc(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  double v = wincon->quad_bfc;
  return v;
}

/* Returns the bottom roughness                                      */
double i_get_z0(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int c2 = window->m2d[c];
  double v = wincon->z0[c2];
  return v;
}

/* Checks if orbital velocities are required to be computed          */
int i_check_orbital_file(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  if (wincon->orbital)
    return(1);
  else
    return(0);
}

/* Checks if the wave period vector is available                     */
int i_check_wave_period(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_period == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_period_m(void *hmodel) 
{
  return(0);
}

/* Returns the wave period at coordinate c                           */
double i_get_wave_period(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_period[c2];
  return v;
}
/* Wrapper for wave period                                           */
double ginterface_getwave_period(void* hmodel, int c)
{
  if (i_check_wave_period(hmodel))
    return i_get_wave_period(hmodel, c) ;
  else
    return 1;
}

/* Checks if the wave direction vector is available                  */
int i_check_wave_dir(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_dir == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_dir_m(void *hmodel) 
{
  return(0);
}

/* Returns the wave direction at coordinate c                        */
double i_get_wave_dir(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_dir[c2];
  return v;
}
/* Wrapper for wave direction                                        */
double ginterface_getwave_dir(void* hmodel, int c)
{
  if (i_check_wave_dir(hmodel))
    return i_get_wave_dir(hmodel, c) ;
  else
    return 0;
}

/* Checks if the wave amplitude vector is available                  */
int i_check_wave_amp(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_amp == NULL)
    return(0);
  else
    return(1);
}

/* Returns the wave amplitude at coordinate c                        */
double i_get_wave_amp(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  return NOTVALID;
}

int i_check_wave_amp_m(void *hmodel) 
{
  return(0);
}

/* Checks if the bottom orbital velociity vector is available        */
int i_check_wave_ub(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ub == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_ub_m(void *hmodel) 
{
  return(0);
}

/* Returns the bottom orbital velocity                               */
double i_get_wave_ub(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ub[c2];
  return v;
}
/* Wrapper for orbital velocity                                      */
double ginterface_getwave_ub(void* hmodel, int c)
{
  if (i_check_wave_ub(hmodel))
    return i_get_wave_ub(hmodel, c);
  else
    return 0;
}

/* Checks if the e1 radiation stress vector is available             */
int i_check_wave_Fx(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_Fx == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_Fx_m(void *hmodel) 
{
  return(0);
}

/* Returns the e1 component of radiation stress at coordinate c      */
double i_get_wave_Fx(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_Fx[c2];
  return v;
}

/* Checks if the e2 radiation stress vector is available             */
int i_check_wave_Fy(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_Fy == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_Fy_m(void *hmodel) 
{
  return(0);
}

/* Returns the e2 component of radiation stress at coordinate c      */
double i_get_wave_Fy(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_Fy[c2];
  return v;
}

/* Checks if surface Stokes e1 velocities are available              */
int i_check_wave_ste1(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ste1 == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_ste1_m(void *hmodel) 
{
  return(0);
}

/* Checks if surface Stokes e2 velocities are available              */
int i_check_wave_ste2(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ste2 == NULL)
    return(0);
  else
    return(1);
}

int i_check_wave_ste2_m(void *hmodel) 
{
  return(0);
}

/* Returns the e1 component of surface Stokes drift at coordinate c  */
double i_get_wave_ste1(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ste1[c2];
  return v;
}

/* Returns the e2 component of surface Stokes drift at coordinate c  */
double i_get_wave_ste2(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ste2[c2];
  return v;
}

/* Checks the wind stress vectors are available                      */
int i_check_wind(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  if (windat->wind1 == NULL || windat->wind2 == NULL)
    return(0);
  else
    return(1);
}

/* Returns the e1 component of wind stress at coordinate c           */
double i_get_wave_wind1(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wind1[c2];
  return v;
}
/* Wrapper for e1 wind stress                                        */
double ginterface_getwind1(void* hmodel, int c)
{
  return i_get_wave_wind1(hmodel, c);
}

/* Returns the e2 component of wind stress at coordinate c           */
double i_get_wave_wind2(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wind2[c2];
  return v;
}
/* Wrapper for e2 wind stress                                        */
double ginterface_getwind2(void* hmodel, int c)
{
  return i_get_wave_wind2(hmodel, c);
}

/* Checks if the wavenumber is available                                  */
int i_check_wave_k(void *hmodel) 
{
  return(0);
}
int i_check_wave_k_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave Bernoulli head vector is available                  */
int i_check_wave_Kb(void *hmodel) 
{
  return(0);
}
int i_check_wave_Kb_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave whitecapping vector is available                  */
int i_check_wave_fwcapx(void *hmodel) 
{
  return(0);
}
int i_check_wave_fwcapx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_fwcapy(void *hmodel) 
{
  return(0);
}
int i_check_wave_fwcapy_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave depth-induced breaking vector is available                  */
int i_check_wave_fbrex(void *hmodel) 
{
  return(0);
}
int i_check_wave_fbrex_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_fbrey(void *hmodel) 
{
  return(0);
}
int i_check_wave_fbrey_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave bottom friction dissapation vector is available                  */
int i_check_wave_fbotx(void *hmodel) 
{
  return(0);
}
int i_check_wave_fbotx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_fboty(void *hmodel) 
{
  return(0);
}
int i_check_wave_fboty_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave surface streaming vector is available                  */
int i_check_wave_fsurx(void *hmodel) 
{
  return(0);
}
int i_check_wave_fsurx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_fsury(void *hmodel) 
{
  return(0);
}
int i_check_wave_fsury_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave form drag vector is available                  */
int i_check_wave_wfdx(void *hmodel) 
{
  return(0);
}
int i_check_wave_wfdx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_wfdy(void *hmodel) 
{
  return(0);
}
int i_check_wave_wfdy_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave ocean viscous stress vector is available                  */
int i_check_wave_wovsx(void *hmodel) 
{
  return(0);
}
int i_check_wave_wovsx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_wovsy(void *hmodel) 
{
  return(0);
}
int i_check_wave_wovsy_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave roller vector is available                  */
int i_check_wave_frolx(void *hmodel) 
{
  return(0);
}
int i_check_wave_frolx_m(void *hmodel) 
{
  return(0);
}
int i_check_wave_froly(void *hmodel) 
{
  return(0);
}
int i_check_wave_froly_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave dummy 1 tracer is available                    */
int i_check_wave_dum1(void *hmodel) 
{
  return(0);
}
int i_check_wave_dum1_m(void *hmodel) 
{
  return(0);
}

/* Checks if the wave dummy 2 tracer is available                    */
int i_check_wave_dum2(void *hmodel) 
{
  return(0);
}
int i_check_wave_dum2_m(void *hmodel) 
{
  return(0);
}

/* Checks if fetch has been computed                                 */
double i_check_fetch(void *hmodel)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  if (wincon->fetch)
    return (1);
  else
    return (0);
}

/* Returns the fetch for index cc                                    */
double i_get_fetch(void *hmodel, int cc, int n) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int c = window->w2_t[cc];
  double v = wincon->fetch[c][n];
  return v;
}

/* Gets the edge normal and tangential volocities given the cell     */
/* centered velocities.                                              */
void i_get_bot_vel(void *hmodel, double sinthcell, double costhcell,
		   double *u1bot, double *u2bot, double *botz, int c)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  int cc = window->c2cc[c2];
  int cb = window->bot_t[cc];

  double ui = (windat->u1[cb] + windat->u1[window->xp1[cb]]) / 2.;
  double uj = (windat->u2[cb] + windat->u2[window->yp1[cb]]) / 2.;
  *u1bot = ui * costhcell - uj * sinthcell;
  *u2bot = ui * sinthcell + uj * costhcell;
  *botz = window->gridz[window->zp1[cb]];
}

void w_get_brsm(void *hmodel, int *brsm) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int cc, c, lc;
  return;
  /*
  for (cc = window->nbe1S + 1; cc <= window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    lc = window->bine1S[cc];
    if (lc == window->ym1[c]) brsm[lc] |= (U2BDRY|F_EDGE);
    if (lc == window->yp1[c]) brsm[lc] |= (U2BDRY|B_EDGE);
  }
  for (cc = window->nbe2S + 1; cc <= window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    lc = window->bine2S[cc];
    if (lc == window->xm1[c]) brsm[lc] |= (U1BDRY|R_EDGE);
    if (lc == window->xp1[c]) brsm[lc] |= (U1BDRY|L_EDGE);
  }
  */
}

/*
 * Moon angle interface
 */
void ginterface_moonvars(void *hmodel, int b,
			 double *mlon, double *mlat,
			 double *dist_earth_sun, double *dist_moon_earth,
			 double *lunar_angle, double *sun_angle,
			 double *moon_phase,
			 double *lunar_dec)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c = i_get_c(hmodel, b);
  int cs = window->m2d[c];
  double Al, dec;
  int yr, mon, nday, jd; /* Year, month, day, Julian date            */
  double day;            /* Day of year                              */
  double hrs;            /* Hours at longitude lon                   */
  double d2r = PI/180.0;  /* Degrees to radians                      */
  double Rl = 3844e5;    /* Mean distance of moon from earth         */
  double Rs = 1496e8;    /* Mean distance of sun from earth          */
  double es = 0.0167086; /* Solar eccentricity                       */
  double el = 0.0549;    /* Lunar eccentricity                       */
  double as = 149.60e9;  /* Sun semi-maor axis (m)                   */
  double lat;            /* Latitude (radians)                       */ 
  double lon = window->cellx[cs];  /* Longitude (degrees)            */
  double gmst;           /* Greenwhich mean sidereal time (hours)    */
  double jt;             /* Julian date                              */
  double hrang[2];       /* Lunar and solar hour angle (rad)         */
  double radius[2];      /* Lunar and solar distance (m)             */

  /* Get the Julian date                                             */
  tm_time_to_ymd(windat->t, master->timeunit, &yr, &mon, &nday);
  jd = date_to_jul(mon, nday, yr);

  /* Julian date starts at midday - adjust half a day                */
  jt = (double)jd - 0.5 ;

  /* Get the day of the year                                         */
  if (window->is_geog) {
    dtime(NULL, master->timeunit, windat->t, &yr, &day, &lon);
  } else {
    hd_quit("Can't compute moon angles for non-geographic mesh.\n");
    return;
  }

  /* Get the Greenwhich mean sidereal time in hours                  */
  /* see https://aa.usno.navy.mil/faq/docs/GAST.php                  */
  jt += windat->days - floor(windat->days);

  // appliedOrbitalHW5.m has 18.697374458 !

  jt -= tm_tz_offset(ginterface_gettimeunits(hmodel)) / 86400.0;  // moonvars is expecting UTC.
  
  gmst = fmod(18.697374558 + 24.06570982441908 * (jt - 2451545.0), 24.0) * PI / 12.0;
  
  /* Call EMS library function */

  moonvars(jt, &Al, &dec, &mlon[0], &mlat[0], &radius[0]);

  radius[0] *= 1e3;   /* Convert to metres */

  /* Get the lunar hour angle (Pugh Eq. 3.20a)                       */
  hrang[0] = lon * d2r + gmst - Al;
  
  /* Get the solar hour angle                                        */
  nday = (int)day;
  hrs = 24 * (day - (double)nday);
  hrang[1] = (hrs - 12.0) * PI / 12.0;
  radius[1] = d2r * day * 360.0 / 365.25;
  radius[1] = as * (1.0 - es * es) / (1.0 + es *cos(radius[1]));

  double moonphase = 0.0; // calculated elsewhere.

  /*
   * Assign outputs
   */
  *dist_moon_earth = radius[0];
  *dist_earth_sun  = radius[1];
  *lunar_angle = hrang[0]; // radians
  *sun_angle   = hrang[1]; // radians
  *moon_phase = moonphase*M_PI; // radians
  *lunar_dec = dec;  //radians

  return;
}

/*
 * Fractional cloud cover
 */
double ginterface_get_cloud(void *hmodel, int b)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c = i_get_c(hmodel, b);
  int cs = window->m2d[c];
  double tcld;

  tcld = min(max((windat->cloud) ? windat->cloud[c] : 0.0, 0.0), 1.0);

  return(tcld);
}

void ginterface_get_windspeed_and_cloud_and_mslp_and_rh(void *hmodel, int b, double *windspeed, double *tcld, double *mslp, double *rh)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c = i_get_c(hmodel, b);
  int cs = window->m2d[c];

  *windspeed = windat->windspeed[cs];
  *tcld = min(max((windat->cloud) ? windat->cloud[cs] : 0.0, 0.0), 1.0);
  *mslp = (windat->patm) ? windat->patm[cs] : 0.0;
  *rh =   (windat->rh) ? windat->rh[cs] : 0.0;
  
  return;
}

void ginterface_get_cloud_and_wind(void *hmodel, int b, double *tcld, double *wind)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c = i_get_c(hmodel, b);
  int cs = window->m2d[c];

  *tcld = min(max((windat->cloud) ? windat->cloud[cs] : 0.0, 0.0), 1.0);

  double wind1 = windat->wind1[cs];
  double wind2 = windat->wind2[cs];

  *wind = 1000.0 * sqrt(wind1*wind1 + wind2*wind2);
    
  return;
}

void ginterface_get_lat_lon(void *hmodel, int b, double *lat, double *lon)
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c = i_get_c(hmodel, b);
  int cs = window->m2d[c];
  
  *lon = window->cellx[cs];  /* Longitude (degrees) */
  *lat = window->celly[cs];  /* Latitude  (degrees) */
  return;
}
#endif

/* Returns the sea surface height                                    */
double ginterface_get_eta(void* hmodel, int b)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c = i_get_c(hmodel, b);
    int c2 = window->m2d[c];
    double v = windat->eta[c2];
    return v;
}
