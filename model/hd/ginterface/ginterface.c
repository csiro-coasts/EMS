/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/ginterface/ginterface.c
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
 *  $Id: ginterface.c 5841 2018-06-28 06:51:55Z riz008 $
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

int tr_map[MAXNUMVARS];
int epi_map[MAXNUMVARS];


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

int i_is_valid_sed_tracer(int itr)
{
  int i;
  if(itr >= 0)
    return (1);
  
  return (0);
}


/* Read parameters from hydromodel */

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

double i_get_model_time(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->t;
}

char *i_get_model_timeunit(void* hmodel)
{
  return(master->timeunit);
}

double i_get_model_timestep(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->dt;
}

int  i_get_num_wclayers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->nz;
}

int  i_get_num_columns(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->b2_t;
}

int i_get_nstep(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    return windat->nstep;
}

int  i_get_num_sedlayers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->sednz;
}

int i_get_num_tracers_3d(void* hmodel, int *atr)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    *atr = master->atr;
    return wincon->ntr;
}

void i_get_names_tracers_3d(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->ntr; n++) {
      strcpy(trname[n], wincon->trinfo_3d[n].name);
    }
}

int i_get_num_tracers_2d(void* hmodel, int *atr)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    *atr = master->atrS;
    return wincon->ntrS;
}

void i_get_names_tracers_2d(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->ntrS; n++) {
      strcpy(trname[n], wincon->trinfo_2d[n].name);
    }
}

int i_get_param_map_3d(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    return wincon->trinfo_3d[n].m;
}

int i_get_param_map_2d(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    if(wincon->trinfo_2d != NULL && wincon->ntrS > 0)
      return wincon->trinfo_2d[n].m;
    else
      return 0;
}

int i_get_param_map_sed(void* hmodel, int n) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    if(wincon->trinfo_sed != NULL && wincon->nsed > 0)
      return wincon->trinfo_sed[n].m;
    else
      return 0;
}

int i_get_num_sedtracers(void* hmodel)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    return wincon->nsed;
}

void i_get_names_tracers_sed(void* hmodel, char **trname) {
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon = window->wincon;
    int n;
    for (n = 0; n < wincon->nsed; n++) {
      strcpy(trname[n], wincon->trinfo_sed[n].name);
    }
}

int i_get_interface_counter(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int sc = cc - 1;
    return sc;
}

int i_get_host_counter(void* hmodel, int c)
{
    return c+1;
}

double i_get_cellarea_w(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    double v = window->cellarea[c2];
    return v;
}

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

int i_get_topk_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int v = window->s2k[window->nsur_t[cc]];
    return v;
}
int i_get_botk_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    int c2 = window->m2d[c];
    int cc = window->c2cc[c2];
    int v = window->s2k[window->bot_t[cc]];
    return v;
}

int i_get_topk_sed(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    return window->sednz;
}

int i_get_botk_sed(void* hmodel, int c)
{
  return 0;
}

double i_get_botz_wc(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    win_priv_t *wincon=window->wincon;
    int c2 = window->m2d[c];
    double v =  window->botz[c2] * wincon->Hn1[c2];
    return v;
}

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

double i_get_depth(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    double v = windat->eta[c2] - window->botz[c2];
    return v;
}

double i_get_eta(void* hmodel, int c)
{
    geometry_t* window = (geometry_t*) hmodel;
    window_t *windat = window->windat;
    int c2 = window->m2d[c];
    double v = windat->eta[c2];
    return v;
}

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
}

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

/*UR-ADDED */
/** Retrive tracer information given a name
 *
 * @param hmodel - the reference to the host model
 * @param trname - the string representing the tracer name
 * @return a tracer_info_t data structure or NULL if the tracer can't be found
 */
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


void i_set_error(void* hmodel, int col, int errorf, char *text)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  wincon->gint_errorf[col] = errorf;
  sprintf(wincon->gint_error[col], "%c", '\0');
  if (text != NULL) // In case of a reset
    strcpy(wincon->gint_error[col], text);
}


int i_get_error(void* hmodel, int col)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  return(wincon->gint_errorf[col]);
}


/*UR in lack of a better place to put it ... */

/* custom_function add and remove */
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
// xxx should these be quits?
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

void trstat_add_domain_function(char* name, int n,int scattertype, custom_create_data create , custom_function_do_gather gather, custom_function_do_scatter scatter)
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
    for (cc = 1; cc <= wincon->vca2; cc++) {
      c = wincon->s2[cc];
      /* Exclude process points if required */
      if (window->ncwave)
	if (ANY(c, window->cwave, window->ncwave)) continue;
      wave_step(window, wincon->wave, c);
    }
    if (wincon->waves & TAN_RAD)
      set_lateral_BC_rad(window, windat, wincon);
    else if (wincon->waves & WAVE_FOR)
       set_lateral_BC_waf(window, windat, wincon);
  }
#endif

}

/* END wave_step()                                                   */
/*-------------------------------------------------------------------*/

double w_get_dt(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  double v = wincon->wavedt;
  return v;
}
double w_get_quad_bfc(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  double v = wincon->quad_bfc;
  return v;
}
double w_get_z0(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int c2 = window->m2d[c];
  double v = wincon->z0[c2];
  return v;
}
int w_check_orbital_file(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  if (wincon->orbital)
    return(1);
  else
    return(0);
}
int w_check_wave_period(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_period == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_period(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_period[c2];
  return v;
}

int w_check_wave_dir(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_dir == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_dir(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_dir[c2];
  return v;
}

int w_check_wave_amp(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_amp == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_amp(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  return NOTVALID;
}
int w_check_wave_ub(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ub == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_ub(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ub[c2];
  return v;
}
int w_check_wave_Fx(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_Fx == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_Fx(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_Fx[c2];
  return v;
}
int w_check_wave_Fy(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_Fy == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_Fy(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_Fy[c2];
  return v;
}
int w_check_wave_ste1(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ste1 == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_ste1(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ste1[c2];
  return v;
}
int w_check_wave_ste2(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  window_t *windat = window->windat;
  if (wincon->waves & NONE || windat->wave_ste2 == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_ste2(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wave_ste2[c2];
  return v;
}

int w_check_wind(void *hmodel) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  if (windat->wind1 == NULL || windat->wind2 == NULL)
    return(0);
  else
    return(1);
}
double w_get_wave_wind1(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wind1[c2];
  return v;
}
double w_get_wave_wind2(void *hmodel, int c) 
{
  geometry_t* window = (geometry_t*) hmodel;
  window_t *windat = window->windat;
  int c2 = window->m2d[c];
  double v = windat->wind2[c2];
  return v;
}
double w_check_fetch(void *hmodel)
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  if (wincon->fetch)
    return (1);
  else
    return (0);
}
double w_get_fetch(void *hmodel, int cc, int n) 
{
  geometry_t* window = (geometry_t*) hmodel;
  win_priv_t *wincon = window->wincon;
  int c = window->w2_t[cc];
  double v = wincon->fetch[c][n];
  return v;
}
double w_get_thetau1(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  double v = window->thetau1[c];
  return v;
}
double w_get_thetau2(void *hmodel, int cc) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int c = window->w2_t[cc];
  double v = window->thetau2[c];
  return v;
}
double w_get_sinthcell(void *hmodel, int cc) 
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
double w_get_costhcell(void *hmodel, int cc) 
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
void w_get_bot_vel(void *hmodel, double sinthcell, double costhcell,
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

int w_get_winsize(void *hmodel)
{
  geometry_t* window = (geometry_t*) hmodel;
  return (window->sgsizS);
}
void w_get_brsm(void *hmodel, int *brsm) 
{
  geometry_t* window = (geometry_t*) hmodel;
  int cc, c, lc;

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
}

#endif
