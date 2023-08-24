/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/tracers/load_tracer.c
 *  
 *  Description:
 *  Initialise tracers and density values.
 *  Used exclusivelt by gendump.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: load_tracer.c 7380 2023-07-26 04:36:43Z her127 $
 *
 */

#include <stdio.h>
#include <math.h>
#include "hd.h"
#include "tracer.h"

void load_wc_tracer_step_3d(parameters_t *params, master_t *master,
                            FILE * fp);
void load_wc_tracer_step_2d(parameters_t *params, master_t *master,
                            FILE * fp);
void load_wc_tracer_defaults_3d(master_t *master);
void load_wc_tracer_defaults_2d(master_t *master);
void calc_density(master_t *master);
void load_sed_tracer_step_3d(master_t *master, FILE * fp);
void load_sed_tracer_defaults(master_t *master);
void read_tsfile(timeseries_t *ts, char *name, char *t_units);
void addatt(int cdfid, int varid, const char*name, const char *text);
void calc_scaling(parameters_t *params, master_t *master, char *mapname);
int face_dist(double *layers, int nz, double depth);
void tr_dataset(char *data, tracer_info_t *tr, double def);
void trn_dataset(char *data, tracer_info_t *trinfo, int tn, int ntr, int atr, double **tr, double def);
void trf_dataset(char *data, tracer_info_t *trinfo, int tn, int ntr, int atr, double **tr, double def);
void dens_profile(master_t *master, geometry_t *geom, double *ret, double v1, double v2);
void dens_grad(master_t *master, geometry_t *geom, double *ret, double v1, double v2);
void dens_scale(master_t *master, geometry_t *geom, double *ret, double v1, double v2, char *sgn);
void dens_scale_file(master_t *master, geometry_t *geom, double *ret, 
		     char *vname, char *infile, double v2, char *sgn);
void dens_grad_file(master_t *master, geometry_t *geom, double *ret, 
		    char *vname, char *infile);
void dens_scale_mid(master_t *master, geometry_t *geom, double *ret, 
		    char *vname, char *infile, double v1, double v2, char *sgn) ;
void tracer_fill(master_t *master, timeseries_t *ts, char *fname, double *ret, int sz, int *vec, int vc, 
		 tracer_info_t *tr, double fill);
int value_init_sparse2d_o(master_t *master, double *ret, char *fname,
			char *vname, char *i_rule);
int value_init_sparse3d_o(master_t *master, double *ret, char *fname,
			char *vname, char *i_rule);
int value_init_sparse2d(master_t *master, double *ret, char *fname, char *vname,
			char *in_rule, double t, int *rask, timeseries_t *ts);
int value_init_sparse3d(master_t *master, double *ret, char *fname, char *vname,
			char *in_rule, double t, int *rask, timeseries_t *ts);
void decode_indata(char *in, char *fname, char *iname, double *t);
void interp_data_s(master_t *master, char *fname, char *vname, double *ret, int mode,
		   int *mask, double t);
void filter_tracer(tracer_info_t *trinfo, geometry_t *window, double **tr, int n);
void glider_scaling(parameters_t *params, master_t *master, char *mapname, int tm);
int find_autotracer_by_name(char *name);
void copy_autotracer_by_name(char *name, tracer_info_t tr[], int ntr, int *n, 
			     double **tra, double **trp);
int set_autotracer_by_name(char *name, tracer_info_t tr[], int ntr, int *n, 
			   double **tra, double **trp, char *buf);
void find_autotracer_by_groupkey(char *name, int *tra, int *n);
int copy_autotracer_by_groupkey(char *name, tracer_info_t tr[], int ntr, int *n);
int set_autotracer_by_groupkey(char *name, tracer_info_t tr[], int ntr, 
			       int *n, char *key);
void update_autotracer(tracer_info_t tr[], int ntr, int *atr, int *mtr);
void typename(int code, char *name);
int set_tracer_3d(parameters_t *params, master_t *master, int ntr, 
		  tracer_info_t *trinfo, double **tr);
int set_tracer_2d(parameters_t *params, master_t *master, int ntr, 
		  tracer_info_t *trinfo, double **tr);
void duplicate_error(char *name, int tn);

extern int NAUTOTR;
extern tracer_info_t autotracerlist[];

/*------------------------------------------------------------------*/
/* Loads 3D tracer values                                           */
/*------------------------------------------------------------------*/
void load_tracer_step_3d(parameters_t *params, master_t *master, FILE * fp)
{

  load_wc_tracer_step_3d(params, master, fp);

  if (geom->sednz > 0)
    load_sed_tracer_step_3d(master, fp);
}

/* END load_tracer_step_3d()                                        */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Loads 2D tracer values                                           */
/*------------------------------------------------------------------*/
void load_tracer_step_2d(parameters_t *params, master_t *master, FILE * fp)
{
  load_wc_tracer_step_2d(params, master, fp);
}

/* END load_tracer_step_2d()                                        */
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/* Read in the tracer data for the water column                     */
/*------------------------------------------------------------------*/
void load_wc_tracer_step_3d(parameters_t *params, master_t *master,
                            FILE * fp)
{
  int c, cs, cc;
  int i, j, k;
  int t;
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  char tag[MAXSTRLEN];
  geometry_t *geom = master->sgrid;

  /*----------------------------------------------------------------*/
  /* Create a scaling file if required                              */
  for (t = master->atr; t < master->ntr; ++t) {
    int tm = master->trinfo_3d[t].m;
    prm_set_errfn(hd_silent_warn);
    sprintf(tag, "TRACER%1.1d.create_scale", tm);
    if (prm_read_char(fp, tag, buf)) {
      calc_scaling(params, master, buf);
    }
  }

  /*----------------------------------------------------------------*/
  /* Read the initial condition for AUTO | DUMP mode                */
  /* Load up the default tracer values */
  load_wc_tracer_defaults_3d(master);
  /*  for (t = master->atr; t < master->ntr; ++t) {*/
  for (t = 0; t < master->ntr; ++t) {
    int tm = master->trinfo_3d[t].m;
    /* Only initialize explicit tracers. Exception is swr_attn.     */
    if (t < master->atr && !strlen(master->trinfo_3d[t].data)) continue;
    prm_set_errfn(hd_silent_warn);
    /*
    sprintf(tag, "TRACER%1.1d.interp_type", tm);
    if (!(prm_read_char(fp, tag, buf2)))
      buf2[0] = '\0';
    */
    sprintf(tag, "TRACER%1.1d.data", tm);
    value_init_3d(master, master->tr_wc[t], fp, master->trinfo_3d[t].data, 
		  master->trinfo_3d[t].name, tag, 
		  master->trinfo_3d[t].fill_value_wc,
		  master->trinfo_3d[t].i_rule);
  }

  /* Load temp and salt for ROAM */
  /* Problem for RECOM using standard file */
  if ( !(params->roammode & (A_RECOM_R1|A_RECOM_R2)))
    temp_salt_init(params, master);

  /* Glider scaling */
  for (t = master->atr; t < master->ntr; ++t) {
    int tm = master->trinfo_3d[t].m;
    sprintf(tag, "TRACER%1.1d.glider_scale", tm);
    if (prm_read_char(fp, tag, buf))
      glider_scaling(params, master, buf, t);
  }

  /* Scaling */
  for (t = master->atr; t < master->ntr; ++t)
    scale_tracer(master->trinfo_3d, geom, master->tr_wc, t);

  /* Filtering */
  for (t = master->atr; t < master->ntr; ++t)
    filter_tracer(master->trinfo_3d, geom, master->tr_wc, t);

  prm_set_errfn(hd_quit);

  calc_density(master);
}

/* END load_wc_tracer_step_3d()                                     */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Scales a tracer                                                  */
/*------------------------------------------------------------------*/
void scale_tracer(tracer_info_t *trinfo, geometry_t *window, double **tr, int n)
{
  int c, cc, id;
  int *vec, nvec;
  int flag = trinfo[n].flag;
  double scale = trinfo[n].scale;

  if (trinfo[n].type & WATER) {
    vec = window->w3_t;
    nvec = window->b3_t;
  }
  if (trinfo[n].type & INTER) {
    vec = window->w2_t;
    nvec = window->b2_t;
  }

  if (flag & SC_SC) {
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] += scale;
    }
  } else if (flag & SC_ST) {
    int id = (int)scale;
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] += tr[id][c];
    }
  } else if (flag & SC_PC) {
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] *= scale;
    }
  } else if (flag & SC_PT) {
    int id = (int)scale;
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] *= tr[id][c];
    }
  }
  /* Clip                                                           */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    tr[n][c] = min(max(trinfo[n].valid_range_wc[0], tr[n][c]), 
		   trinfo[n].valid_range_wc[1]);
  }
}

/* END scale_tracer()                                               */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Filters a tracer                                                 */
/*------------------------------------------------------------------*/
void filter_tracer(tracer_info_t *trinfo, geometry_t *geom, double **tr, int n)
{
  int c, cc, cp, cm, id, i, m;
  int *vec, nvec;
  int flag = trinfo[n].flag;
  double scale = 0.25;
  double cval = trinfo[n].scale;

  if (trinfo[n].type & WATER) {
    vec = geom->w3_t;
    nvec = geom->b3_t;
  }
  if (trinfo[n].type & INTER) {
    vec = geom->w2_t;
    nvec = geom->b2_t;
  }

  if (flag & V_HP) {
    double f[9];
    double *b = d_alloc_1d(geom->sgsiz);

    /* Set the boundary conditions (no flux) */
    for (cc = 1; cc <= geom->nbpt; cc++) {
      c = geom->bpt[cc];
      cp = geom->bin[cc];
      tr[n][c] = tr[n][cp];
    }
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->bot_t[cc];
      tr[n][geom->zm1[c]] = tr[n][c];
    }
    for (m = 0; m < geom->nobc; m++) {
      open_bdrys_t *open = geom->open[m];
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	cm = (open->type & U1BDRY) ? open->obc_e1[cc] : open->obc_e2[cc];
	tr[n][cm] = tr[n][c];
      }
    }

    /* Set the filter kernal */
    f[0] = f[2] = -1.0;
    f[3] = f[4] = f[5] = 0.0;
    f[6] = f[8] = 1.0;
    f[1] = -cval;
    f[7] = cval;

    for (cc = 1; cc <= nvec; cc++) {
      double a[9];
      c = vec[cc];
      cp = geom->zp1[c];
      cm = geom->zm1[c];
      /*
      a[4] = tr[n][c];
      a[5] = tr[n][geom->xp1[c]];
      a[3] = tr[n][geom->xm1[c]];
      a[1] = tr[n][cp];
      a[2] = tr[n][geom->xp1[cp]];
      a[0] = tr[n][geom->xm1[cp]];
      a[7] = tr[n][cm];
      a[8] = tr[n][geom->xp1[cm]];
      a[6] = tr[n][geom->xm1[cm]];
      */
      b[c] = tr[n][c];
      for (i = 0; i < 9; i++)
	b[c] += scale * (a[i] * f[i]);
    }
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] = b[c];
    }

    for (cc = 1; cc <= nvec; cc++) {
      double a[9];
      c = vec[cc];
      cp = geom->zp1[c];
      cm = geom->zm1[c];
      a[4] = tr[n][c];
      /*
      a[5] = tr[n][geom->yp1[c]];
      a[3] = tr[n][geom->ym1[c]];
      a[1] = tr[n][cp];
      a[2] = tr[n][geom->yp1[cp]];
      a[0] = tr[n][geom->ym1[cp]];
      a[7] = tr[n][cm];
      a[8] = tr[n][geom->yp1[cm]];
      a[6] = tr[n][geom->ym1[cm]];
      */
      b[c] = tr[n][c];
      for (i = 0; i < 9; i++)
	b[c] += scale * (a[i] * f[i]);
    }
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      tr[n][c] = b[c];
    }
     d_free_1d(b);
  }
}

/* END filter_tracer()                                              */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialises T & S from input file                                */
/*------------------------------------------------------------------*/
void temp_salt_init(parameters_t *params, master_t *master)
{
  int c, cs, cc;
  int t;
  int tns = 1, tnt = 1;
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  geometry_t *geom = master->sgrid;

  if (!strlen(params->tdata)) {
    hd_warn("temp_salt_init: No temp forcing data specified for TEMP_DATA\n");
    tnt = 0;
  }
  if (!strlen(params->sdata)) {
    hd_warn("temp_salt_init: No salt forcing data specified for SALT_DATA\n");
    tns = 0;
  }

  /* Load temp and salt for ROAM */
  for (t = 0; t < master->ntr; ++t) {
    timeseries_t *ts;
    int id;
    
    if (tns && strcmp(master->trinfo_3d[t].name, "salt") == 0)
      strcpy(buf, params->sdata);
    else if(tnt && strcmp(master->trinfo_3d[t].name, "temp") == 0)
      strcpy(buf, params->tdata);
    else
      continue;
    
    ts = hd_ts_read(master, buf, 0);
    if (value_init_sparse3d(master, master->tr_wc[t], buf, 
			    master->trinfo_3d[t].name, 
			    master->trinfo_3d[t].i_rule,
			    schedule->start_time, NULL, ts) == 0)
      continue;

    id = ts_get_index(ts, fv_get_varname(buf,
					 master->trinfo_3d[t].name, buf2));
    if (id >= 0) {
      strcpy(master->trinfo_3d[t].data, buf);
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	master->tr_wc[t][c] = ts_eval_xyz(ts, id, master->t,
					  geom->cellx[cs],
					  geom->celly[cs],
					  geom->cellz[c]);
      }
    } else {
      hd_warn
	("temp_salt_init: Using default values for '%s' tracer.\n",
	 master->trinfo_3d[t].name);
    }
  }
  calc_density(master);
}
/* END temp_salt_init()                                             */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Read in the tracer data for the 2D variables                     */
/*------------------------------------------------------------------*/
void load_wc_tracer_step_2d(parameters_t *params, master_t *master,
                            FILE * fp)
{
  int c, cs, cc;
  int i, j, k;
  int t;
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  char tag[MAXSTRLEN];
  geometry_t *geom = master->sgrid;

  /* Load up the default tracer values */
  load_wc_tracer_defaults_2d(master);

  prm_set_errfn(hd_silent_warn);
  for (t = 0; t < master->ntrS; ++t) {
    int tm = master->trinfo_2d[t].m;
    /*
    sprintf(tag, "TRACER%1.1d.interp_type", tm);
    if (!(prm_read_char(fp, tag, buf2)))
      buf2[0] = '\0';
    */
    strcpy(buf2, master->trinfo_2d[t].i_rule);
    sprintf(tag, "TRACER%1.1d.data", tm);

    value_init_2d(master, master->tr_wcS[t], fp, master->trinfo_2d[t].data, 
		  master->trinfo_2d[t].name, tag, 
		  master->trinfo_2d[t].fill_value_wc, buf2);

    /* Scaling */
    scale_tracer(master->trinfo_2d, geom, master->tr_wcS, t);
  }

  prm_set_errfn(hd_quit);
}

/* END load_wc_tracer_step_2d()                                     */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialise the tracers to the values specified in the tracer     */
/* fill_value.                                                      */
/*------------------------------------------------------------------*/
void load_wc_tracer_defaults_3d(master_t *master)
{
  int n, c, cs, cc;
  geometry_t *geom = master->sgrid;

  /* Fill watercolumn values */
  for (n = 0; n < master->ntr; n++) {
    memset(master->tr_wc[n],0,geom->sgsiz*sizeof(double));
    if(master->trinfo_3d[n].fill_value_wc != 0.0) {
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
        master->tr_wc[n][c] = master->trinfo_3d[n].fill_value_wc;
      }
    }
  }
}

/* END load_wc_tracer_defaults_3d()                                    */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialise the 2D tracers to the values specified in the tracer  */
/* fill_value.                                                      */
/*------------------------------------------------------------------*/
void load_wc_tracer_defaults_2d(master_t *master)
{
  int n, c, cc;
  geometry_t *geom = master->sgrid;

  /* Fill 2D tracer values */
  for (n = 0; n < master->ntrS; n++) {
    memset(master->tr_wcS[n],0,geom->sgsizS*sizeof(double));
    if(master->trinfo_2d[n].fill_value_wc != 0.0) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
        master->tr_wcS[n][c] = master->trinfo_2d[n].fill_value_wc;
      }
    }
  }
}

/* END load_wc_tracer_defaults_2d()                                 */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Read in the tracer data for a specific tracer 'trname'           */
/*------------------------------------------------------------------*/
void load_wc_tracer_name(master_t *master, FILE * fp, char *trname, int dim)
			  
{
  int t, tm, ntr, trn = -1;
  double fill;
  char fname[MAXSTRLEN];
  char vname[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  char tag[MAXSTRLEN];
  tracer_info_t *tr;

  if (dim == WATER) {
    tr = master->trinfo_3d;
    strcpy(vname, trname);
    ntr = master->ntr;
  } else if (dim == INTER) {
    tr = master->trinfo_2d;
    strcpy(vname, trname);
    ntr = master->ntrS;
  } else if (dim == SEDIM) {
    tr = master->trinfo_sed;
    strcpy(vname, trname);
    strcat(vname, "_sed");
    ntr = master->nsed;
  }

  /*----------------------------------------------------------------*/
  /* Save the name of the data file in the parameter structure */
  for (t = 0; t < ntr; t++) {
    tm = tr[t].m;
    if (strcmp(tr[t].name, trname) == 0) {
      trn = t;
      break;
    }
  }
  if (trn == -1)
    return;

  if (dim == WATER) {
    fill = master->trinfo_3d[trn].fill_value_wc;
    strcpy(i_rule, master->trinfo_3d[trn].i_rule);
  } else if (dim == INTER) {
    fill = master->trinfo_2d[trn].fill_value_wc;
    strcpy(i_rule, master->trinfo_2d[trn].i_rule);
  } else if (dim == SEDIM) {
    fill = master->trinfo_sed[trn].fill_value_sed;
    strcpy(i_rule, master->trinfo_sed[trn].i_rule);
  }

  /* Search for a data file for this variable in the prm file     */
  /*
  sprintf(tag, "TRACER%1.1d.interp_type", tm);
  if (!(prm_read_char(fp, tag, i_rule)))
    i_rule[0] = '\0';
  */
  if (dim == SEDIM)
    sprintf(tag, "TRACER%1.1d.data_sed", tm);
  else
    sprintf(tag, "TRACER%1.1d.data", tm);
  /*
  if (!(prm_read_char(fp, tag, fname)))
    fname[0] = '\0';
  */
  strcpy(fname, tr[trn].data);

  if (dim == WATER) {
    value_init_3d(master, master->tr_wc[trn], fp, fname, vname, tag, fill, i_rule);
    /* Glider scaling */
    for (t = master->atr; t < master->ntr; ++t) {
      int tm = master->trinfo_3d[t].m;
      char buf[MAXSTRLEN];
      sprintf(tag, "TRACER%1.1d.glider_scale", tm);
      if (prm_read_char(fp, tag, buf))
	glider_scaling(params, master, buf, t);
    }
    /* Scaling */
    for (t = master->atr; t < master->ntr; ++t)
      scale_tracer(master->trinfo_3d, geom, master->tr_wc, t);
  } else if (dim == INTER)
    value_init_2d(master, master->tr_wcS[trn], fp, fname, vname, tag, fill, i_rule);
  else if (dim == SEDIM)
    value_init_sed(master, master->tr_sed[trn], fp, fname, vname, tag, fill, i_rule);

  prm_set_errfn(hd_quit);
}

/* END load_wc_tracer_name()                                        */
/*------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/* To add a new autotracer:                                          */
/* 1) Define the tracer attributes in autotracer.c                   */
/* 2) Copy the attributes to trinfo_3d when the keyword is invoked   */
/*    in set_tracer_2d(), assign pointers, and perform any custom    */
/*    initialisation in this routine.                                */
/* 3) Assign pointers for the windows in win_data_init().            */
/*-------------------------------------------------------------------*/
void init_tracer_2d(parameters_t *params, /* Input parameters data   */
                    master_t *master      /* Master data             */
  )
{
  int tn;
  geometry_t *geom = master->geom;
  tn = 0;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  master->alert_a = master->alert_c = NULL;
  master->avhrr = master->ghrsst = master->ghrsste = master->shwin = master->shinx = master->layth = NULL;
  master->cfl2d = master->cfl3d = master->cour = master->courn = master->lips = master->ahsb = NULL;
  master->mixl = master->steric = NULL;
  master->av = master->rv = master->pv = NULL;
  master->rv_nonlin = master->rv_beta = master->rv_strch = NULL;
  master->rv_drvdt = master->rv_jebar = master->rv_wsc = master->rv_bsc =
    NULL;
  master->w1m = master->w2m = master->u1am = master->u2am = master->etam = NULL;
  master->nhfd = master->swrd = master->lwrd = master->lhfd =
    master->shfd = NULL;
  master->nsfd = NULL;
  master->rossby_ex = master->speed_2d = master->wind_Cd = NULL;
  master->obc_phase = master->speed_sq = NULL;
  master->u1_rad = master->u2_rad = master->sonic = master->wetcell = NULL;
  master->surfz = master->slope_x = NULL;
  master->tau_be1 = master->tau_be2 = master->tau_bm =  NULL;
  master->sederr = master->ecoerr = master->riverflow = master->riverdepth = NULL;
  master->bathy_range_max = master->bathy_range_min = NULL;
  master->bathy_grad_max = master->bathy_grad_min = master->eta_tc = master->eta_inc = NULL;
  master->cellres = master->equitide =   master->tpxotide = master->vhreg = NULL;
  master->tpxovelu = master->tpxovelv = master->tpxotranu = master->tpxotranv = NULL;
  master->uat = master->vat = master->meshun = NULL;
  master->carea = master->earea = master->sarea = master->searea = NULL;
  master->feta = NULL;

  /* Waves                                                           */
  master->ustrcw = master->wave_amp = master->wave_Cd = NULL;
  master->wave_ub = master->wave_period = master->wave_dir = NULL;
  master->wave_Sxy = master->wave_Syx = master->vol_cons = NULL;
  master->wave_Fx = master->wave_Fy = NULL;
  master->wave_ste1 = master->wave_ste2 = NULL;
  master->wave_P = master->wave_Kb = master->wave_k = NULL;
  master->wave_fwcapx = master->wave_fwcapy = master->wave_fbrex = master->wave_fbrey = NULL;
  master->wave_fbotx = master->wave_fboty = master->wave_fsurx = master->wave_fsury = NULL;
  master->wave_wfdx = master->wave_wfdy = master->wave_wovsx = master->wave_wovsy = NULL;
  master->wave_frolx = master->wave_froly = NULL;
  master->tau_w1 = master->tau_w2 = master->tau_diss1 = master->tau_diss2 = NULL;
  master->sep = master->bep = master->tfront = NULL;

  /* SWR                                                             */
  master->swr_attn1 = master->swr_tran = NULL;
  master->swr_babs = master->swreg = master->swrms = master->attn_mean = master->tran_mean = NULL;
  if (params->swr_type & SWR_2D) master->swr_attn = NULL;

  /*-----------------------------------------------------------------*/
  /* Assign the automtically generated 2D pointers and set up the    */
  /* tracer_info structure. If any of these tracers have been        */
  /* explicitly defined in the input parameter file, then continue.  */
  tn = set_tracer_2d(params, master, master->ntrS, master->trinfo_2d, master->tr_wcS);
  /* Sediment error maps                                             */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    master->sederr = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "sed_error");
    strcpy(master->trinfo_2d[tn].long_name, "Sediment error");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|SEDIMENT|DIAGNOSTIC;
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    memset(master->sederr, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
    /* Set up sediment tracers if required                           */
    tn = sediment_autotracer_2d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, master->trinfo_2d, master->ntrS, tn);
  }
#endif
  /* Ecology error maps                                              */
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    master->ecoerr = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "eco_error");
    strcpy(master->trinfo_2d[tn].long_name, "Ecology error");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|ECOLOGY|DIAGNOSTIC;
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    memset(master->ecoerr, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
    /* Set up sediment tracers if required                           */
    tn = ecology_autotracer_2d(params->prmfd, params->do_eco, params->eco_vars, 
			       params->eco_defs, params->pre_eco, 
			       master->trinfo_2d, master->ntrS, tn);
  }
#endif

}

/* END init_tracer_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the sediment tracers in the master          */
/*-------------------------------------------------------------------*/
void init_tracer_sed(parameters_t *params, /* Input parameters data  */
		     master_t *master      /* Master data            */
  )
{
  int tn;
  geometry_t *geom = master->sgrid;
  char buf[MAXSTRLEN];
  tn = 0;

  /* Initialise */

  /* Assign the 2D pointers.  */
  sprintf(buf, "%c", '\0');

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    /* Set up sediment tracers if required */
    tn = sediment_autotracer_sed(params->prmfd, params->do_sed, params->sed_vars, 
				 params->sed_defs, master->trinfo_sed, master->nsed, tn);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    /* Set up ecology tracers if required */
    tn = ecology_autotracer_sed(params->prmfd, params->do_eco, 
				params->eco_vars, params->eco_defs, 
				params->pre_eco, 
				master->trinfo_sed, master->nsed, tn);
  }
#endif

}

/* END init_tracer_sed()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 3D tracers in the master.               */
/* Note: 3D tracers are first defined in params->trinfo_3d in the    */
/* routine create_tracer_3d(), copied to master->trinfo_3d in        */
/* hd_init(), then pointers are assigned here.                       */
/* To add a new autotracer:                                          */
/* 1) Define the tracer attributes in autotracer.c                   */
/* 2) Copy the attributes to trinfo_3d when the keyword is invoked   */
/*    in set_tracer_3d(). Note; capability exists for this routine   */
/*    to assign pointers, but we don't use this here.                */
/* 3) Assign ponters to the tracer and define any custom             */
/*    initialisation in this routine.                                */
/* 4) Assign pointers for the windows in win_data_build().           */
/*-------------------------------------------------------------------*/
void init_tracer_3d(parameters_t *params, /* Input parameters data   */
                    master_t *master      /* Master data             */
  )
{
  int n, tn;
  geometry_t *geom = master->sgrid;
  char buf[MAXSTRLEN], name[MAXSTRLEN];
  if (geom == NULL) geom = master->geom;

  /*-----------------------------------------------------------------*/
  /* Initialise pointers                                             */
  master->sal = master->temp = NULL;
  master->tke = master->diss = master->L = master->omega = NULL;
  master->Q2 = master->Q2L = master->Kq = master->sdc = NULL;
  master->u1m = master->u2m = master->wm = master->Kzm = NULL;
  master->fluxe1 = master->fluxe2 = master->tempm = master->saltm = NULL;
  master->tram = master->fluxw = master->fluxkz = NULL;
  master->brunt = master->int_wave = master->rich_gr = master->rich_fl = NULL;
  master->froude = master->reynolds = master->rossby_in = master->sound = NULL;
  master->shear_v = master->b_prod = master->s_prod = master->speed_3d = NULL;
  master->otemp = master->osalt = master->perc = NULL;
  master->fltr = master->agetr = NULL;
  master->rtemp = master->rsalt = master->schan = master->sigma_t = NULL;
  master->ptconc = master->energy = master->vcorr = master->acorr = NULL;
  master->dum1 = master->dum2 = master->dum3 = master->kenergy = NULL;
  master->riversalt = master->regionid = master->regres = NULL;
  master->Vi = master->reefe1 = master->reefe2 = NULL;
  master->temp_tc = master->salt_tc = master->unit = master->mono = NULL;
  master->wave_stke1 = master->wave_stke2 = NULL;
  master->fsalt = master->ftemp = master->fvelu = master->fvelv = NULL;
  if (params->swr_type & SWR_3D) master->swr_attn = NULL;
  master->glider = master->nprof = master->u1vhc = NULL;
  master->volcont = master->centi = NULL;
  if (params->ndhw) {
    master->dhw = (double **)p_alloc_1d(params->ndhw);
    master->dhd = (double **)p_alloc_1d(params->ndhw);
    master->dhwc = (double **)p_alloc_1d(params->ndhw);
    for (n = 0; n < params->ndhw; n++) {
      master->dhw[n] = NULL;
      master->dhd[n] = NULL;
      master->dhwc[n] = NULL;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Assign the water column tracer pointers to any tracers defined  */
  /* in the input parameter file.                                    */
  for (tn = 0; tn < master->ntr; tn++) {
    strcpy(name, master->trinfo_3d[tn].name);

    if (strcmp("salt", name) == 0) {
      master->sal = master->tr_wc[tn];
      master->sno = tn;
    }
    else if (strcmp("temp", name) == 0) {
      master->temp = master->tr_wc[tn];
      master->tno = tn;
    }
    else if (strcmp("tke", name) == 0) {
      master->tke = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("diss", name) == 0) {
      master->diss = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("omega", name) == 0) {
      master->omega = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("tki", name) == 0) {
      master->Q2 = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("tki_l", name) == 0) {
      master->Q2L = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("Kq", name) == 0) {
      master->Kq = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    } else if (strcmp("lscale", name) == 0) {
      master->L = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    } else if (strcmp("u1mean", name) == 0) {
      master->u1m = master->tr_wc[tn];
      memset(master->u1m, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("u2mean", name) == 0) {
      master->u2m = master->tr_wc[tn];
      memset(master->u2m, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("wmean", name) == 0) {
      master->wm = master->tr_wc[tn];
      memset(master->wm, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("temp_mean", name) == 0) {
      master->tempm = master->tr_wc[tn];
      memset(master->tempm, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("salt_mean", name) == 0) {
      master->saltm = master->tr_wc[tn];
      memset(master->saltm, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("tracer_mean", name) == 0) {
      master->tram = master->tr_wc[tn];
      memset(master->tram, 0, geom->sgsiz * sizeof(double));
      sprintf(buf, "Mean %s", params->means_tra);
      strcpy(master->trinfo_3d[tn].long_name, buf);
    } else if (strcmp("Kzmean", name) == 0) {
      master->Kzm = master->tr_wc[tn];
      memset(master->Kzm, 0, geom->sgsiz * sizeof(double));
    } else if (strcmp("flux_e1", name) == 0) {
      master->fluxe1 = master->tr_wc[tn];
      sprintf(master->trinfo_3d[tn].long_name, "Advective flux through face %d", 
	      params->trfd1);
    } else if (strcmp("flux_e2", name) == 0) {
      master->fluxe2 = master->tr_wc[tn];
      sprintf(master->trinfo_3d[tn].long_name, "Advective flux through face %d", 
	      params->trfd2);
    } else if (strcmp("flux_w", name) == 0)
      master->fluxw = master->tr_wc[tn];
    else if (strcmp("flux_kz", name) == 0)
      master->fluxkz = master->tr_wc[tn];
    else if (strcmp("flush", name) == 0)
      master->fltr = master->tr_wc[tn];
    else if (strcmp("age", name) == 0)
      master->agetr = master->tr_wc[tn];
    else if (strcmp("smagorinsky", name) == 0) {
      master->sdc = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    }
    else if (strcmp("brunt_vaisala", name) == 0)
      master->brunt = master->tr_wc[tn];
    else if (strcmp("int_wave_speed", name) == 0)
      master->int_wave = master->tr_wc[tn];
    else if (strcmp("richardson_gr", name) == 0)
      master->rich_gr = master->tr_wc[tn];
    else if (strcmp("richardson_fl", name) == 0)
      master->rich_fl = master->tr_wc[tn];
    else if (strcmp("reynolds", name) == 0)
      master->reynolds = master->tr_wc[tn];
    else if (strcmp("froude", name) == 0)
      master->froude = master->tr_wc[tn];
    else if (strcmp("sigma_t", name) == 0)
      master->sigma_t = master->tr_wc[tn];
    else if (strcmp("energy", name) == 0)
      master->energy = master->tr_wc[tn];
    else if (strcmp("kenergy", name) == 0)
      master->kenergy = master->tr_wc[tn];
    else if (strcmp("ptconc", name) == 0)
      master->ptconc = master->tr_wc[tn];
    else if (strcmp("sound", name) == 0)
      master->sound = master->tr_wc[tn];
    else if (strcmp("sound_channel", name) == 0)
      master->schan = master->tr_wc[tn];
    else if (strcmp("shear_vert", name) == 0)
      master->shear_v = master->tr_wc[tn];
    else if (strcmp("buoy_prod", name) == 0)
      master->b_prod = master->tr_wc[tn];
    else if (strcmp("shear_prod", name) == 0)
      master->s_prod = master->tr_wc[tn];
    else if (strcmp("rossby_in", name) == 0)
      master->rossby_in = master->tr_wc[tn];
    else if (strcmp("current_speed_3d", name) == 0)
      master->speed_3d = master->tr_wc[tn];
    else if (strcmp("otemp", name) == 0) {
      master->otemp = master->tr_wc[tn];
      if (strlen(params->odata)) {
	if (params->save_force & ROAM)
	  sprintf(buf, "%s(otemp=temp)", params->tdata);
	else
	  sprintf(buf, "%s", params->odata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
      }
    } else if (strcmp("osalt", name) == 0) {
      master->osalt = master->tr_wc[tn];
      if (strlen(params->odata)) {
	if (params->save_force & ROAM)
	  sprintf(buf, "%s(osalt=salt)", params->sdata);
	else
	  sprintf(buf, "%s", params->odata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
      }
    } else if (strcmp("rtemp", name) == 0) {
      master->rtemp = master->tr_wc[tn];
      if (strlen(params->tdata)) {
	char buf1[MAXSTRLEN];
	if (params->save_force & ROAM) {
	  sprintf(buf, "%s(rtemp=temp)", params->tdata);
	  strcpy(buf1, "1 day");
	}
	else {
	  sprintf(buf, "%s", params->tdata);
	  strcpy(buf1, "1 second");
	}
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, buf1);
      }
    } else if (strcmp("rsalt", name) == 0) {
      master->rsalt = master->tr_wc[tn];
      if (strlen(params->sdata)) {
	char buf1[MAXSTRLEN];
	if (params->save_force & ROAM) {
	  sprintf(buf, "%s(rsalt=salt)", params->sdata);
	  strcpy(buf1, "1 day");
	}
	else {
	  sprintf(buf, "%s", params->sdata);
	  strcpy(buf1, "1 second");
	}
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, buf1);
      }
    } else if (strcmp("ovelu", name) == 0) {
      if (strlen(params->vdata)) {
	if (params->save_force & ROAM)
	  sprintf(buf, "%s(ovelu=u)", params->vdata);
	else
	  sprintf(buf, "%s", params->vdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
      }
    } else if (strcmp("ovelv", name) == 0) {
      if (strlen(params->vdata)) {
	if (params->save_force & ROAM)
	  sprintf(buf, "%s(ovelv=v)", params->vdata);
	else
	  sprintf(buf, "%s", params->vdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
      }
    } else if (strcmp("salt_force", name) == 0) {
      master->fsalt = master->tr_wc[tn];
      if (params->save_force & FSALT) {
	sprintf(buf, "%s(salt_force=salt)", params->sdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, params->fsalt_input_dt);
	strcpy(master->trinfo_3d[tn].reset_interp, params->fsalt_interp);
      }
    } else if (strcmp("temp_force", name) == 0) {
      master->ftemp = master->tr_wc[tn];
      if (params->save_force & FTEMP) {
	sprintf(buf, "%s(temp_force=temp)", params->tdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, params->ftemp_input_dt);
	strcpy(master->trinfo_3d[tn].reset_interp, params->ftemp_interp);
      }
    } else if (strcmp("velu_force", name) == 0) {
      master->fvelu = master->tr_wc[tn];
      if (params->save_force & FVELU) {
	sprintf(buf, "%s(velu_force=u)", params->vdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, params->fvelu_input_dt);
	strcpy(master->trinfo_3d[tn].reset_interp, params->fvelu_interp);
      }
    } else if (strcmp("velv_force", name) == 0) {
      master->fvelv = master->tr_wc[tn];
      if (params->save_force & FVELV) {
	sprintf(buf, "%s(velv_force=v)", params->vdata);
	strcpy(master->trinfo_3d[tn].reset_file, buf);
	strcpy(master->trinfo_3d[tn].reset_dt, params->fvelv_input_dt);
	strcpy(master->trinfo_3d[tn].reset_interp, params->fvelv_interp);
      }
    } else if (strcmp("imp3df", name) == 0) {
      strcpy(master->trinfo_3d[tn].name, params->imp3dn);
      strcpy(master->trinfo_3d[tn].long_name, params->imp3dn);
      strcpy(master->trinfo_3d[tn].units, params->imp3du);
      sprintf(buf, "%c", '\0');
      tr_dataset(buf, &master->trinfo_3d[tn], 0.0);
      strcpy(master->trinfo_3d[tn].i_rule, "nn_sibson");
      if (strlen(params->imp3dt))
	sprintf(master->trinfo_3d[tn].data, "[data=%s(t=%s)]", params->imp3df, 
		params->imp3dt);
      else
	sprintf(master->trinfo_3d[tn].data, "[data=%s]", params->imp3df);
    } else if (strcmp("temp_tc", name) == 0)
      master->temp_tc = master->tr_wc[tn];
    else if (strcmp("salt_tc", name) == 0)
      master->salt_tc = master->tr_wc[tn];
    else if (strcmp("tracer1", name) == 0)
      master->dum1 = master->tr_wc[tn];
    else if (strcmp("tracer2", name) == 0)
      master->dum2 = master->tr_wc[tn];
    else if (strcmp("tracer3", name) == 0)
      master->dum3 = master->tr_wc[tn];
    else if (strcmp("reef_fraction_e1", name) == 0)
      master->reefe1 = master->tr_wc[tn];
    else if (strcmp("reef_fraction_e2", name) == 0)
      master->reefe2 = master->tr_wc[tn];
    else if (strcmp("Vcorr", name) == 0)
      master->vcorr = master->tr_wc[tn];
    else if (strcmp("Acorr", name) == 0)
      master->acorr = master->tr_wc[tn];
    else if (strcmp("Vi", name) == 0)
      master->Vi = master->tr_wc[tn];
    else if (strcmp("unit", name) == 0)
      master->unit = master->tr_wc[tn];
    else if (strcmp("regionid", name) == 0)
      master->regionid = master->tr_wc[tn];
    else if (strcmp("residence", name) == 0)
      master->regres = master->tr_wc[tn];
    else if (strcmp("decorr_e1", name) == 0)
      master->decv1 = master->tr_wc[tn];
    else if (strcmp("wave_stke1", name) == 0)
      master->wave_stke1 = master->tr_wc[tn];
    else if (strcmp("wave_stke2", name) == 0)
      master->wave_stke2 = master->tr_wc[tn];
    else if (strcmp("glider", name) == 0)
      master->glider = master->tr_wc[tn];
    else if (strcmp("u1vhc", name) == 0)
      master->u1vhc = master->tr_wc[tn];
    else if (strcmp("vol_cont", name) == 0)
      master->volcont = master->tr_wc[tn];
    else if (strcmp("cell_index", name) == 0)
      master->centi = master->tr_wc[tn];
    else if (strcmp("nprof", name) == 0)
      master->nprof = master->tr_wc[tn];
    else if (strcmp("mono", name) == 0) {
      master->mono = master->tr_wc[tn];
      sprintf(master->trinfo_3d[tn].long_name, "Monotinicity of %s", 
	      params->monotr);
    } else if (params->swr_type & SWR_3D && strcmp("swr_attenuation", name) == 0) {
      master->swr_attn = master->tr_wc[tn];
      trn_dataset(params->swr_attn, master->trinfo_3d, tn, master->ntr, 
		  master->atr, master->tr_wc, 0.073);
    } else if (params->ndhw) {
      for (n = 0; n < params->ndhw; n++) {
	sprintf(buf, "dhw%d", n);
	if (strcmp(buf, name) == 0)
	  master->dhw[n] = master->tr_wc[tn];
	sprintf(buf, "dhd%d", n);
	if (strcmp(buf, name) == 0)
	  master->dhd[n] = master->tr_wc[tn];
	sprintf(buf, "dhwc%d", n);
	if (strcmp(buf, name) == 0)
	  master->dhwc[n] = master->tr_wc[tn];
      }
    }
  }
  sprintf(buf, "percentile_%s", params->trperc);
  for (tn = 0; tn < master->ntr; tn++) {
    if (strcmp(buf, master->trinfo_3d[tn].name) == 0) {
      master->perc = master->tr_wc[tn];
      strcpy(master->trinfo_3d[n].name, buf);
      sprintf(buf, "Percentile for %s", params->trperc);
      strcpy(master->trinfo_3d[n].long_name, buf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Assign the automtically generated 3D pointers and set up the    */
  /* tracer_info structure. If any of these tracers have been        */
  /* explicitly defined in the input parameter file, then continue.  */
  /*
  set_tracer_3d(params, master, master->ntr, master->trinfo_3d, master->tr_wc);
  */
  /*
  for(tn=0; tn<master->ntr; tn++)
    printf("%d %s\n",tn,master->trinfo_3d[n].name);
  */
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    /* Set up sediment tracers if required                           */
    tn = sediment_autotracer_3d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, master->trinfo_3d, master->ntr, tn);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    /* Set up ecology tracers if required                            */
    tn = ecology_autotracer_3d(params->prmfd, params->do_eco, params->eco_vars,
			       params->eco_defs, params->pre_eco,
			       master->trinfo_3d, master->ntr, tn);
  }
#endif

}

/* END init_tracer_3d()                                              */
/*-------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Read in the tracer data for the specified tracer                 */
/*------------------------------------------------------------------*/
void load_sed_tracer_step_3d(master_t *master, FILE * fp)
{
  int k, c, cc;
  int t;
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  char tag[MAXSTRLEN];
  geometry_t *geom = master->sgrid;

  /* Load up the default tracer values */
  load_sed_tracer_defaults(master);

  prm_set_errfn(hd_silent_warn);

  for (t = 0; t < master->nsed; ++t) {
    int tm = master->trinfo_sed[t].m;

    /*
    sprintf(tag, "TRACER%1.1d.interp_type", tm);
    if (!(prm_read_char(fp, tag, buf2)))
      buf2[0] = '\0';
    */
    strcpy(buf2, master->trinfo_sed[t].i_rule);
    sprintf(tag, "TRACER%1.1d.data_sed", tm);
    strcpy(buf, master->trinfo_sed[t].name);
    strcat(buf, "_sed");
    value_init_sed(master, master->tr_sed[t], fp, master->trinfo_sed[t].data, 
		   buf, tag, master->trinfo_sed[t].fill_value_sed, buf2);
  }
  prm_set_errfn(hd_quit);
}

/* END load_sed_tracer_step_3d()                                    */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Initialise the tracers to the values specified in the tracer     */
/* fill_value. Also compute the density.                            */
/*------------------------------------------------------------------*/
void load_sed_tracer_defaults(master_t *master)
{
  int k, c, cc, n;
  geometry_t *geom = master->sgrid;

  /* Fill sediment values */
  if (geom->sednz == 0)
    return;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    for (k = 0; k < geom->sednz; k++) {
      for (n = 0; n < master->nsed; ++n) {
	if (master->trinfo_sed[n].fill_value_sed) {
	  master->tr_sed[n][k][c] = master->trinfo_sed[n].fill_value_sed;
	} else {
	  master->tr_sed[n][k][c] = 0.0;
	}
      }
    }
  }
}

/* END load_sed_tracer_defaults()                                   */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Compute the density from the salt and temp tracers               */
/*------------------------------------------------------------------*/
void calc_density(master_t *master)
{
  int cc, c;
  geometry_t *geom = master->sgrid;
  int ndens = -master->calc_dens;

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    if (ndens > 0)
      master->dens[c] = master->tr_wc[ndens][c];
    else if (master->sal != NULL && master->temp != NULL)
      master->dens[c] = eos(master->sal[c],
                            master->temp[c],
                            -geom->cellz[c] * 9.81 * 1025.0);
    else
      master->dens[c] = 1025.0;
  }
}

/* END calc_density()                                               */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Sets up the autotracers from the global set in autotracer.c,     */
/* and sets the pointers for the master to appropriate tracers.     */
/* This is called to populate the master->trinfo_3d and             */
/* params->trinfo_3d datastructures. If a tracer array (**tr) is    */
/* not NULL, then pointers from the master can also be set to the   */
/* relevant tracer variables (this is not done for params trinfo).  */
/* Note: when using a groupkey, then set pointers to variables in   */
/* the same order they are found in the autotracerlist.             */
/*------------------------------------------------------------------*/
int set_tracer_3d(parameters_t *params, 
		  master_t *master,
		  int ntr, 
		  tracer_info_t *trinfo, 
		  double **tr)
{
  geometry_t *geom;
  char buf[MAXSTRLEN];
  int n, m, ntrp, tn = 0;

  if (master != NULL) geom= master->geom;

  copy_autotracer_by_name("salt", trinfo, ntr, &tn, tr, &master->sal);
  copy_autotracer_by_name("temp", trinfo, ntr, &tn, tr, &master->temp);

  if (params->means & VEL3D) {
    n = copy_autotracer_by_groupkey("vel3d_mean", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->u1m = tr[n++];
      master->u2m = tr[n++];
      master->wm = tr[n++];
      memset(master->u1m, 0, geom->sgsiz * sizeof(double));
      memset(master->u2m, 0, geom->sgsiz * sizeof(double));
      memset(master->wm, 0, geom->sgsiz * sizeof(double));
    }
  }
  if (params->means & KZ_M) {
    copy_autotracer_by_name("Kzmean", trinfo, ntr, &tn, tr, &master->Kzm);
    if (tr != NULL)
      memset(master->Kzm, 0, geom->sgsiz * sizeof(double));
  }
  if (params->means & TS) {
    n = copy_autotracer_by_groupkey("ts_mean", trinfo, ntr, &tn);	
    if (tr != NULL) {
      master->tempm = tr[n++];
      master->saltm = tr[n++];
      memset(master->tempm, 0, geom->sgsiz * sizeof(double));
      memset(master->saltm, 0, geom->sgsiz * sizeof(double));
    }
  }
  if (params->means & MTRA3D) {
    n = tn;
    copy_autotracer_by_name("tracer_mean", trinfo, ntr, &tn, 
			    tr, &master->tram);
    sprintf(buf, "Mean %s", params->means_tra);
    strcpy(master->trinfo_3d[n].long_name, buf);
    memset(master->tram, 0, geom->sgsiz * sizeof(double));
  }
  if (strcmp(params->trflux, "NONE") != 0) {
    n = m = copy_autotracer_by_groupkey("tr_flux", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->fluxe1 = tr[n++]; 
      master->fluxe2 = tr[n++]; 
      master->fluxw = tr[n++];
      master->fluxkz = tr[n++];
    }
    sprintf(master->trinfo_3d[m++].long_name, "Advective flux through face %d", 
	    params->trfd1);
    sprintf(master->trinfo_3d[m++].long_name, "Advective flux through face %d", 
	    params->trfd2);
  }
  if (strlen(params->regions)) {
    n = copy_autotracer_by_groupkey("regions", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->regionid = tr[n++];
      master->regres = tr[n++];
    }
  }
  if (params->trflsh) {
    copy_autotracer_by_name("flush", trinfo, ntr, &tn, tr, &master->fltr);
  }
  if (strlen(params->trage)) {
    copy_autotracer_by_name("age", trinfo, ntr, &tn, tr, &master->agetr);
  }
  if (strcmp(params->trperc, "NONE") != 0) {
    n = tn;
    sprintf(buf, "percentile_%s", params->trperc);
    copy_autotracer_by_name(buf, trinfo, ntr, &tn, tr, &master->perc);
    strcpy(master->trinfo_3d[n].name, buf);
    sprintf(buf, "Percentile for %s", params->trperc);
    strcpy(trinfo[n].long_name, buf);
  }
  if (strcmp(params->mixsc, "k-e") == 0) {
    n = copy_autotracer_by_groupkey("k-e", trinfo, ntr, &tn); 
    if (tr != NULL) {
      master->tke = tr[n++];
      master->diss = tr[n++];
    }
  }
  if (strcmp(params->mixsc, "k-w") == 0) {
    n = copy_autotracer_by_groupkey("k-w", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->tke = tr[n++];
      master->omega = tr[n++];
    }
  }
  if (strcmp(params->mixsc, "W88") == 0) {
    n = copy_autotracer_by_groupkey("k-w", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->tke = tr[n++];
      master->omega = tr[n++];
    }
  }
  if (strcmp(params->mixsc, "mellor_yamada_2_0") == 0 ||
      strcmp(params->mixsc, "mellor_yamada_2_0_estuarine") == 0) {
    n = copy_autotracer_by_groupkey("my2.0", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->L = tr[n++];
    }
  }
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0 || 
      strcmp(params->mixsc, "harcourt") == 0) {
    n = copy_autotracer_by_groupkey("my2.5", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->Q2 = tr[n++];
      master->Q2L = tr[n++];
      master->L = tr[n++];
      master->Kq = tr[n++];
    }
  }
  if (params->smagorinsky > 0.0) {
    copy_autotracer_by_name("smagorinsky", trinfo, ntr, &tn, 
			    tr, &master->sdc);
  }
  if (params->show_layers) {
    copy_autotracer_by_name("layer_thick", trinfo, ntr, &tn, tr, &master->layth);
  }
  if (params->save_force & OTEMP) {
    n = tn;
    copy_autotracer_by_name("otemp", trinfo, ntr, &tn, tr, &master->otemp);
    if (strlen(params->odata)) {
      if (params->save_force & ROAM)
	sprintf(buf, "%s(otemp=temp)", params->tdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(trinfo[n].reset_file, buf);
      strcpy(trinfo[n].reset_dt, "1 second");
    }
  }
  if (params->save_force & OSALT) {
    n = tn;
    copy_autotracer_by_name("osalt", trinfo, ntr, &tn, tr, &master->osalt);
   if (strlen(params->odata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(osalt=salt)", params->sdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(trinfo[n].reset_file, buf);
      strcpy(trinfo[n].reset_dt, "1 second");
    }
  }
  if (params->rtemp) {
    n = tn;
    copy_autotracer_by_name("rtemp", trinfo, ntr, &tn, tr, &master->rtemp);
    if (strlen(params->tdata)) {
      char buf1[MAXSTRLEN];
      if (params->save_force & ROAM) {
	sprintf(buf, "%s(rtemp=temp)", params->tdata);
	strcpy(buf1, "1 day");
      }
      else {
	sprintf(buf, "%s", params->tdata);
	strcpy(buf1, "1 second");
      }
      strcpy(trinfo[n].reset_file, buf);
      strcpy(trinfo[n].reset_dt, buf1);
    }
  }
  if (params->rsalt) {
    n = tn;
    copy_autotracer_by_name("rsalt", trinfo, ntr, &tn, tr, &master->rsalt);
    if (strlen(params->sdata)) {
      char buf[MAXSTRLEN];
      char buf1[MAXSTRLEN];
      if (params->save_force & ROAM) {
	sprintf(buf, "%s(rsalt=salt)", params->sdata);
	strcpy(buf1, "1 day");
      }
      else {
	sprintf(buf, "%s", params->sdata);
	strcpy(buf1, "1 second");
      }
      strcpy(trinfo[n].reset_file, buf);
      strcpy(trinfo[n].reset_dt, buf1);
    }
  }
  if (params->save_force & FTEMP) {
    n = tn;
    copy_autotracer_by_name("temp_force", trinfo, ntr, &tn, tr, &master->ftemp);
    sprintf(buf, "%s(temp_force=temp)", params->tdata);
    strcpy(trinfo[n].reset_file, buf);
    strcpy(trinfo[n].reset_dt, params->ftemp_input_dt);
    strcpy(trinfo[n].reset_interp, params->ftemp_interp);
  }
  if (params->save_force & FSALT) {
    n = tn;
    copy_autotracer_by_name("salt_force", trinfo, ntr, &tn, tr, &master->fsalt);
    sprintf(buf, "%s(salt_force=salt)", params->sdata);
    strcpy(trinfo[n].reset_file, buf);
    strcpy(trinfo[n].reset_dt, params->fsalt_input_dt);
    strcpy(trinfo[n].reset_interp, params->fsalt_interp);
  }
  if (params->save_force & FVELU) {
    n = tn;
    copy_autotracer_by_name("velu_force", trinfo, ntr, &tn, tr, &master->fvelu);
    sprintf(buf, "%s(velu_force=u)", params->vdata);
    strcpy(trinfo[n].reset_file, buf);
    strcpy(trinfo[n].reset_dt, params->fvelu_input_dt);
    strcpy(trinfo[n].reset_interp, params->fvelu_interp);
  }
  if (params->save_force & FVELV) {
    n = tn;
    copy_autotracer_by_name("velv_force", trinfo, ntr, &tn, tr, &master->fvelv);
    sprintf(buf, "%s(velv_force=v)", params->vdata);
    strcpy(trinfo[n].reset_file, buf);
    strcpy(trinfo[n].reset_dt, params->fvelv_input_dt);
    strcpy(trinfo[n].reset_interp, params->fvelv_interp);
  }
  if (params->rtemp & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    copy_autotracer_by_name("temp_tc", trinfo, ntr, &tn, tr, &master->temp_tc);
  }  
  if (params->rsalt & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    copy_autotracer_by_name("salt_tc", trinfo, ntr, &tn, tr, &master->salt_tc);
  }
  if (params->save_force & OVELU) {
    copy_autotracer_by_name("ovelu", trinfo, ntr, &tn, tr, NULL);
    if (strlen(params->vdata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(ovelu=u)", params->vdata);
      else
	sprintf(buf, "%s", params->vdata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
  }
  if (params->save_force & OVELV) {
    copy_autotracer_by_name("ovelv", trinfo, ntr, &tn, tr, NULL);
    if (strlen(params->vdata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(ovelv=v)", params->vdata);
      else
	sprintf(buf, "%s", params->vdata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
  }



  if (params->tendf) {
    copy_autotracer_by_groupkey("tend", trinfo, ntr, &tn);
    if (params->waves & STOKES_DRIFT) {
      copy_autotracer_by_groupkey("tend_wave", trinfo, ntr, &tn);
    }
  }
  if (strlen(params->trtend)) {
    copy_autotracer_by_groupkey("tr_tend", trinfo, ntr, &tn);
  }
  if (params->waves & SPECTRAL) {
    n = copy_autotracer_by_groupkey("spec_wave", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->wave_stke1 = tr[n++]; 
      master->wave_stke2 = tr[n++]; 
    }
  }
  if (params->numbers & BRUNT) {
    copy_autotracer_by_name("brunt", trinfo, ntr, &tn, 
			    tr, &master->brunt);
  }
  if (params->numbers & INT_WAVE) {
    copy_autotracer_by_name("int_wave_speed", trinfo, ntr, &tn, 
			    tr, &master->int_wave);
  }
  if (params->numbers & RICHARD_GR) {
    copy_autotracer_by_name("richardson_gr", trinfo, ntr, &tn, 
			    tr, &master->rich_gr);
  }
  if (params->numbers & RICHARD_FL) {
    copy_autotracer_by_name("richardson_fl", trinfo, ntr, &tn, 
			    tr, &master->rich_fl);
  }
  if (params->numbers & REYNOLDS) {
    copy_autotracer_by_name("reynolds", trinfo, ntr, &tn, 
			    tr, &master->reynolds);
  }
  if (params->numbers & FROUDE) {
    copy_autotracer_by_name("froude", trinfo, ntr, &tn, 
			    tr, &master->froude);
  }
  if (params->numbers & SIGMA_T) {
    copy_autotracer_by_name("sigma_t", trinfo, ntr, &tn, 
			    tr, &master->sigma_t);
  }
  if (params->numbers & ENERGY) {
    copy_autotracer_by_name("energy", trinfo, ntr, &tn, 
			    tr, &master->energy);
  }
  if (params->numbers & KINETIC) {
    copy_autotracer_by_name("kenergy", trinfo, ntr, &tn, 
			    tr, &master->kenergy);
  }
  if (params->do_pt) {
    copy_autotracer_by_name("ptconc", trinfo, ntr, &tn, 
			    tr, &master->ptconc);
  }
  if (params->numbers & SOUND) {
    n = copy_autotracer_by_groupkey("sound", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->sound = tr[n++];
      master->schan = tr[n++];
    }
  }
  if (params->numbers & ROSSBY_IN) {
    copy_autotracer_by_name("rossby_internal", trinfo, ntr, &tn, 
			    tr, &master->rossby_in);
  }
  if (params->numbers & SPEED_3D) {
    copy_autotracer_by_name("current_speed_3d", trinfo, ntr, &tn, 
			    tr, &master->speed_3d);
  }
  if (params->numbers & SHEAR_V) {
    copy_autotracer_by_name("shear_vert", trinfo, ntr, &tn, 
			    tr, &master->shear_v);
  }
  if (params->numbers & BUOY_PROD) {
    copy_autotracer_by_name("buoy_prod", trinfo, ntr, &tn, 
			    tr, &master->b_prod);
  }
  if (params->numbers & SHEAR_PROD) {
    copy_autotracer_by_name("shear_prod", trinfo, ntr, &tn, 
			    tr, &master->s_prod);
  }
  if (!(params->decf & (NONE|DEC_ETA))) {
    copy_autotracer_by_name("decorr_e1", trinfo, ntr, &tn, 
			    tr, &master->decv1);
  }
  if (strlen(params->imp3df)) {
    n = tn;
    copy_autotracer_by_name("imp3df", trinfo, ntr, &tn, tr, NULL);
    strcpy(trinfo[n].name, params->imp3dn);
    strcpy(trinfo[n].long_name, params->imp3dn);
    strcpy(trinfo[n].units, params->imp3du);
    tr_dataset(buf, &trinfo[n], 0.0);
    strcpy(trinfo[n].i_rule, "nn_sibson");
    if (strlen(params->imp3dt))
      sprintf(trinfo[n].data, "[data=%s(t=%s)]", params->imp3df, 
	      params->imp3dt);
    else
      sprintf(trinfo[n].data, "[data=%s]", params->imp3df);
  }
  if (params->numbers & DUMMIES) {
    n = copy_autotracer_by_groupkey("dummies", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->dum1 = tr[n++];
      master->dum2 = tr[n++];
      master->dum3 = tr[n++];
    }
  }
  if (params->numbers & UNIT) {
    copy_autotracer_by_name("unit", trinfo, ntr, &tn, tr, &master->unit);
  }
  if (params->numbers & PASS) {
    copy_autotracer_by_name("passive", trinfo, ntr, &tn, tr, NULL);
  }
  if (params->numbers & GLIDER) {
    copy_autotracer_by_name("glider", trinfo, ntr, &tn, tr, &master->glider);
  }
  if (params->numbers1 & U1VHC) {
    copy_autotracer_by_name("u1vhc", trinfo, ntr, &tn, tr, &master->u1vhc);
  }
  if (params->numbers1 & VOLCONT) {
    copy_autotracer_by_name("vol_cont", trinfo, ntr, &tn, tr, &master->volcont);
  }
  if (params->numbers1 & CENTI) {
    copy_autotracer_by_name("cell_index", trinfo, ntr, &tn, tr, &master->centi);
  }
  if (strlen(params->nprof)) {
    copy_autotracer_by_name("nprof", trinfo, ntr, &tn, tr, &master->nprof);
  }
  if (strlen(params->monotr)) {
    n = tn;
    copy_autotracer_by_name("mono", trinfo, ntr, &tn, tr, &master->mono);
    sprintf(trinfo[n].long_name, "Monotinicity of %s", params->monotr);
  }
  if (params->porusplate) {
    n = copy_autotracer_by_groupkey("porusplate", trinfo, ntr, &tn);
    if (tr != NULL) {
      master->reefe1 =  tr[n++];
      master->reefe2 =  tr[n++];
    }
  }
  if (params->riverflow == 2) {
    copy_autotracer_by_name("flow_salt", trinfo, ntr, &tn, 
			    tr, &master->riversalt);
  }
  if (params->swr_type & SWR_3D && strlen(params->swr_attn)) {
    n = tn;
    if (tracer_find_index("swr_attenuation", ntr, trinfo) == -1) {
      copy_autotracer_by_name("swr_attenuation", trinfo, ntr, &tn, 
			      tr, &master->swr_attn);
      if (tr != NULL)
	trn_dataset(params->swr_attn, trinfo, n, master->ntr, 
		    master->atr, tr, 0.073);
    } else {
      int tm = tracer_find_index("swr_attenuation", ntr, trinfo);
      if (tr != NULL)
	trn_dataset(params->swr_attn, trinfo, tm, master->ntr, 
		    master->atr, tr, 0.073);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Degree heating weeks - special treatment as we don't know how   */
  /* many of these there will be.                                    */
  for (n = 0; n < params->ndhw; n++) {
    sprintf(buf, "dhw%d", n);    
    if (tracer_find_index(buf, ntr, trinfo) == -1) {
      if (tr != NULL) master->dhw[n] = master->tr_wc[tn];
      strcpy(trinfo[tn].name, buf);
      if (strlen(params->dhwt[n]))
	sprintf(buf, "Degree heating week #%d %s", n, params->dhwt[n]);
      else
	sprintf(buf, "Degree heating week #%d", n);
      strcpy(trinfo[tn].long_name, buf);
      strcpy(trinfo[tn].units, "DegC-week");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      /* Set DA_ tag for dhw so that slave-master transfers don't occur. */
      strcpy(trinfo[tn].tag, "DA_");
      trinfo[tn].valid_range_wc[0] = -1e10;
      trinfo[tn].valid_range_wc[1] = 1e10;
      if (params->dhwf[n] & DHW_RT)
	sprintf(trinfo[tn].tracerstat, "exposure(temp:dhwc%d:dhwt%d)",n,n);
      trinfo[tn].m = -1;
      trinfo[tn].n = tn;
      tn++;
    }
    if (params->dhwf[n] & DHW_RT)
      sprintf(buf, "dhwt%d", n);
    if (params->dhwf[n] & DHW_NOAA)
      sprintf(buf, "dhd%d", n);
    if (tracer_find_index(buf, ntr, trinfo) == -1) {
      if (tr != NULL) master->dhd[n] = master->tr_wc[tn];
      strcpy(trinfo[tn].name, buf);

      if (params->dhwf[n] & DHW_RT)
	sprintf(buf, "Degree heating exposure time #%d", n);
      if (params->dhwf[n] & DHW_NOAA)
	sprintf(buf, "Degree heating day #%d",n);
      strcpy(trinfo[tn].long_name, buf);

      if (params->dhwf[n] & DHW_RT)
	strcpy(buf, "week");
      if (params->dhwf[n] & DHW_NOAA)
	strcpy(buf, "DegCday");
      strcpy(trinfo[tn].units, buf);

      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 7;
      trinfo[tn].m = -1;
      trinfo[tn].n = tn;
      tn++;
    }
    sprintf(buf, "dhwc%d", n);    
    if (tracer_find_index(buf, ntr, trinfo) == -1) {
      if (tr != NULL) master->dhwc[n] = master->tr_wc[tn];
      strcpy(trinfo[tn].name, buf);
      sprintf(buf, "Degree heating threshold #%d", n);
      strcpy(trinfo[tn].long_name, buf);
      strcpy(trinfo[tn].units, "Degrees C");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = -1e10;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = -1;
      trinfo[tn].n = tn;
      strcpy(trinfo[tn].reset_file, params->dhw[n]);
      strcpy(trinfo[tn].reset_dt, "1 day");
      tn++;
    }
  }
  params->ntr = ntr;
  return(tn);
}

/* END set_tracer_3d()                                               */
/*-------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Sets up the autotracers from the global set in autotracer.c,     */
/* and sets the pointers for the master to appropriate tracers.     */
/* This is called to populate the master->trinfo_2d.                */
/*------------------------------------------------------------------*/
int set_tracer_2d(parameters_t *params, 
		   master_t *master,
		   int ntr, 
		   tracer_info_t *trinfo, 
		   double **tr)
{
  geometry_t *geom;
  char buf[MAXSTRLEN];
  int n, m, ntrp, tn = 0;

  if (master != NULL) geom= master->geom;
  sprintf(buf, "%c", '\0');

  if (!(params->cfl & NONE)) {
    n = set_autotracer_by_groupkey("cfl", trinfo, ntr, &tn, buf);
    master->cfl2d = tr[n++];
    master->cfl3d = tr[n++];
    master->cour = tr[n++];
    master->lips = tr[n++];
    master->ahsb = tr[n++];
    master->courn = tr[n++];
  }
  if (!(params->mixlayer & NONE))
    set_autotracer_by_name("mixed_layer", trinfo, ntr, &tn, tr, &master->mixl, buf);
  if (params->lnm != 0.0)
    set_autotracer_by_name("steric", trinfo, ntr, &tn, tr, &master->steric, buf);
  if (params->vorticity & ABSOLUTE)
    set_autotracer_by_name("abs_vor", trinfo, ntr, &tn, tr, &master->av, buf);
  if (params->vorticity & RELATIVE)
    set_autotracer_by_name("rel_vor", trinfo, ntr, &tn, tr, &master->rv, buf);
  if (params->vorticity & POTENTIAL)
    set_autotracer_by_name("pot_vor", trinfo, ntr, &tn, tr, &master->pv, buf);
  if (params->vorticity & TENDENCY) {
    n = set_autotracer_by_groupkey("vorticity", trinfo, ntr, &tn, buf);
    master->rv_drvdt = tr[n++];
    master->rv_nonlin = tr[n++];
    master->rv_beta = tr[n++];
    master->rv_strch = tr[n++];
    master->rv_jebar = tr[n++];
    master->rv_wsc = tr[n++];
    master->rv_bsc = tr[n++];
  }
  if (params->diff_scale & VH_REG)
    set_autotracer_by_name("u1vh_region", trinfo, ntr, &tn, tr, &master->vhreg, buf);
  if (params->numbers & ROSSBY_EX)
    set_autotracer_by_name("rossby_external", trinfo, ntr, &tn, tr, &master->rossby_ex, buf);
  if (params->numbers & SPEED_2D)
    set_autotracer_by_name("current_speed_2d", trinfo, ntr, &tn, tr, &master->speed_2d, buf);
  if (params->numbers & SPEED_SQ)
    set_autotracer_by_name("speed_sq", trinfo, ntr, &tn, tr, &master->speed_sq, buf);
  if (params->numbers & OBC_PHASE)
    set_autotracer_by_name("obc_phase", trinfo, ntr, &tn, tr, &master->obc_phase, buf);
  if (params->numbers & WIND_CD)
    set_autotracer_by_name("wind_Cd", trinfo, ntr, &tn, tr, &master->wind_Cd, buf);
  if (params->numbers & CELLRES) {
    n = set_autotracer_by_groupkey("resolution", trinfo, ntr, &tn, buf);
    master->cellres = tr[n++];
    master->sarea = tr[n++];
    master->searea = tr[n++];
  }
  if (params->numbers1 & CELLAREA) {
    n = set_autotracer_by_groupkey("area", trinfo, ntr, &tn, buf);
    master->carea = tr[n++];
    master->earea = tr[n++];
  }
  if (params->numbers1 & MESHUN)
    set_autotracer_by_name("mesh_uniformity", trinfo, ntr, &tn, tr, &master->meshun, buf);
  if (params->waves & (TAN_RAD|WAVE_FOR) && params->tendf) {
    n = set_autotracer_by_groupkey("rad_stress", trinfo, ntr, &tn, buf);
    master->u1_rad = tr[n++];
    master->u2_rad = tr[n++];
  }
  if (params->means & ETA_M)
    set_autotracer_by_name("eta_mean", trinfo, ntr, &tn, tr, &master->etam, buf);
  if (params->save_force & FETA) {
    n = set_autotracer_by_name("eta_force", trinfo, ntr, &tn, tr, &master->feta, buf);
    sprintf(buf, "%s(eta_force=eta)", params->edata);
    strcpy(master->trinfo_2d[n].reset_file, buf);
    strcpy(master->trinfo_2d[n].reset_dt, params->feta_input_dt);
    strcpy(master->trinfo_2d[n].reset_interp, params->feta_interp);
  }
  if (params->means & WIND) {
    n = set_autotracer_by_groupkey("wind_mean", trinfo, ntr, &tn, buf);
    master->w1m = tr[n++];
    master->w2m = tr[n++];
  }
  if (params->means & VEL2D) {
    n = set_autotracer_by_groupkey("vel2d_mean", trinfo, ntr, &tn, buf);
    master->u1am = tr[n++];
    master->u2am = tr[n++];
  }
  if (params->means & MTRA2D) {
    char key[MAXSTRLEN];
    n = set_autotracer_by_name("tracer_mean_2d", trinfo, ntr, &tn, tr, &master->tram, buf);
    strcpy(master->trinfo_2d[n].name, "tracer_mean");
    sprintf(key, "Mean %s", params->means_tra);
    strcpy(master->trinfo_2d[n].long_name, key);
    memset(master->tram, 0, geom->sgsizS * sizeof(double));
  }
  if (params->heatflux & ADVANCED) {
    n = set_autotracer_by_groupkey("heatflux", trinfo, ntr, &tn, buf);
    master->nhfd = tr[n++];
    master->swrd = tr[n++];
    master->lwrd = tr[n++];
    master->lhfd = tr[n++];
    master->shfd = tr[n++];
  }
  if (params->heatflux & (INVERSE|COMP_HEAT|COMP_HEAT_MOM))
    set_autotracer_by_name("nhf", trinfo, ntr, &tn, tr, &master->nhfd, buf);
  if (params->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE)) {
    n = set_autotracer_by_groupkey("heatcomp", trinfo, ntr, &tn, buf);
    if (params->heatflux & (COMP_HEAT_MOM | COMP_HEAT_NONE)) {
      /* See logic in heatflux.c:comp_heat_mom */
      master->swr = tr[n++];
    } else 
      master->swrd = tr[n++];
    master->lwrn = n++;
    master->lhfn = n++;
    master->shfn = n++;
    if (params->heatflux & COMP_HEAT_NONE) {
      if (strlen(params->precip)) {
	n = set_autotracer_by_name("precip", trinfo, ntr, &tn, tr, NULL, buf);
	master->precipn = n;
      }
      if (strlen(params->evap)) {
	n = set_autotracer_by_name("evap", trinfo, ntr, &tn, tr, NULL, buf);
	master->evapn = n;
      }
    }
  }
  if (params->heatflux & NET_HEAT) {
    n = set_autotracer_by_groupkey("netheat", trinfo, ntr, &tn, buf);
    master->nhfd = tr[n++];
    master->swrd = tr[n++];
  }
  if (params->saltflux & (ADVANCED | BULK | ORIGINAL))
    set_autotracer_by_name("nsf", trinfo, ntr, &tn, tr, &master->nsfd, buf);
  if (params->saltflux & (ADVANCED | ORIGINAL)) {
    n = set_autotracer_by_groupkey("saltflux", trinfo, ntr, &tn, buf);
    master->precipn = n++;
    master->evapn = n++;
  }
  if (params->waves & BOT_STR)
    set_autotracer_by_name("wave_Cd", trinfo, ntr, &tn, tr, &master->wave_Cd, buf);
  if (params->waves & TAN_RAD) {
    n = set_autotracer_by_groupkey("rad_force", trinfo, ntr, &tn, buf);
    master->wave_Sxy = tr[n++];
    master->wave_Syx = tr[n++];
  }
  if (params->waves & WAVE_FOR) {
    n = set_autotracer_by_groupkey("wave_force", trinfo, ntr, &tn, buf);
    master->wave_Fx = tr[n++];
    master->wave_Fy = tr[n++];
  }
  if (params->waves & (STOKES|SPECTRAL)) {
    n = set_autotracer_by_groupkey("wave_stokes", trinfo, ntr, &tn, buf);
    master->wave_ste1 = tr[n++];
    master->wave_ste2 = tr[n++];
    if (params->waves & STOKES_DRIFT) {
      n = set_autotracer_by_groupkey("wave_stress", trinfo, ntr, &tn, buf);
      master->tau_w1 = tr[n++];
      master->tau_w2 = tr[n++];
      master->tau_diss1 = tr[n++];
      master->tau_diss2 = tr[n++];
    }
  }
  if (params->waves & NEARSHORE) {
    n = set_autotracer_by_groupkey("wave_nearshore", trinfo, ntr, &tn, buf);
    master->wave_Kb = tr[n++];
    master->wave_k = tr[n++];
    master->wave_P = tr[n++];
    master->wave_fwcapx = tr[n++];
    master->wave_fbrex = tr[n++];
    master->wave_fbotx = tr[n++];
    master->wave_fsurx = tr[n++];
    master->wave_wfdx = tr[n++];
    master->wave_wovsx = tr[n++];
    master->wave_frolx = tr[n++];
    master->wave_fwcapy = tr[n++];
    master->wave_fbrey = tr[n++];
    master->wave_fboty = tr[n++];
    master->wave_fsury = tr[n++];
    master->wave_wfdy = tr[n++];
    master->wave_wovsy = tr[n++];
    master->wave_froly = tr[n++];
  }
  if (!(params->do_wave & NONE)) {
    n = set_autotracer_by_groupkey("waves", trinfo, ntr, &tn, buf);
    master->ustrcw = tr[n++];
    master->wave_ub = tr[n++];
    master->wave_period = tr[n++];
    master->wave_dir = tr[n++];
    master->wave_amp = tr[n++];
  }
  if (params->etarlx & (RELAX|ALERT|BOUNDARY))
    set_autotracer_by_name("oeta", trinfo, ntr, &tn, tr, NULL, buf);
  if (params->etarlx & ETA_TPXO || params->etarlx & ETA_ADPT)
    set_autotracer_by_name("eta_tc", trinfo, ntr, &tn, tr, &master->eta_tc, buf);
  if (params->etarlx & ETA_ADPT)
    set_autotracer_by_name("eta_inc", trinfo, ntr, &tn, tr, &master->eta_inc, buf);
  if (params->avhrr)
    set_autotracer_by_name("AVHRR", trinfo, ntr, &tn, tr, &master->avhrr, buf);
  if (params->ghrsst) {
    n = set_autotracer_by_groupkey("ghrsst", trinfo, ntr, &tn, buf);
    strcpy(master->trinfo_2d[n].reset_file, params->ghrsst_path);
    strcpy(master->trinfo_2d[n].reset_dt, params->ghrsst_dt);
    if (strlen(params->ghrsst_irule))
      strcpy(master->trinfo_2d[n].reset_interp, params->ghrsst_irule);
    master->ghrsst = tr[n++];
    master->ghrsste = tr[n++];
  }
  if (params->show_win) {
    /* Window index tracer currently not used, as it's the same as  */
    /* cell_index. In the surface layer cc and w2_t[cc] are the     */
    /* same for windows.                                            */
    n = set_autotracer_by_groupkey("windiag", trinfo, ntr, &tn, buf);
    master->shwin = tr[n++];
    master->shinx= tr[n++];
    /*set_autotracer_by_name("windows", trinfo, ntr, &tn, tr, &master->shwin, buf);*/
  }
  if (strlen(params->bathystats)) {
    n = set_autotracer_by_groupkey("bathystat", trinfo, ntr, &tn, buf);
    master->bathy_range_min = tr[n++];
    master->bathy_range_max = tr[n++];
    master->bathy_grad_min = tr[n++];
    master->bathy_grad_max = tr[n++];
  }
  if (contains_token(params->alert, "ACTIVE") != NULL) {
    n = set_autotracer_by_groupkey("alerts", trinfo, ntr, &tn, buf);
    master->alert_a = tr[n++];
    master->alert_c = tr[n++];
    master->u1vhin = tr[n++];
    master->u2vhin = tr[n++];
  }
  if (params->fillf & (WEIGHTED|MONOTONIC) || (params->tmode & SP_FFSL))
    set_autotracer_by_name("vol_cons", trinfo, ntr, &tn, tr, &master->vol_cons, buf);
  if (params->numbers & SOUND)
    set_autotracer_by_name("sonic_depth", trinfo, ntr, &tn, tr, &master->sonic, buf);
  if (params->numbers & EKPUMP) {
    n = set_autotracer_by_groupkey("ekman", trinfo, ntr, &tn, buf);
    master->sep = tr[n++];
    master->bep = tr[n++];
  }
  if (params->numbers & TIDEFR)
    set_autotracer_by_name("SH_tide_front", trinfo, ntr, &tn, tr, &master->tfront, buf);
  if (params->numbers & WET_CELLS)
    set_autotracer_by_name("wet_cells", trinfo, ntr, &tn, tr, &master->wetcell, buf);
  if (params->numbers & SURF_LAYER)
    set_autotracer_by_name("surf_layer", trinfo, ntr, &tn, tr, &master->surfz, buf);
  if (params->numbers & SLOPE)
    set_autotracer_by_name("surf_slope", trinfo, ntr, &tn, tr, &master->slope_x, buf);
  if (params->numbers & BOTSTRESS) {
    n = set_autotracer_by_groupkey("botstress", trinfo, ntr, &tn, buf);
    master->tau_be1 = tr[n++];
    master->tau_be2 = tr[n++];
    master->tau_bm = tr[n++];
  }
  if (strlen(params->swr_babs)) {
    n = tn;
    copy_autotracer_by_name("swr_bot_absorb", trinfo, ntr, &tn, tr, &master->swr_babs);
    trn_dataset(params->swr_babs, trinfo, n, params->ntrS, params->atrS, tr, 1.0);
  }
  if (params->swr_type & SWR_2D && strlen(params->swr_attn)) {
    n = tn;
    copy_autotracer_by_name("swr_attenuation", trinfo, ntr, &tn, tr, &master->swr_attn);
    if (strlen(params->swr_regions))
      trf_dataset(params->swr_attn, trinfo, n, params->ntrS, params->atrS, tr, 0.073);
    else
      trn_dataset(params->swr_attn, trinfo, n, params->ntrS, params->atrS, tr, 0.073);
  }
  if (strlen(params->swr_attn1)) {
    n = tn;
    copy_autotracer_by_name("swr_deep_attenuation", trinfo, ntr, &tn, tr, &master->swr_attn1);
    tr_dataset(params->swr_attn1, &trinfo[n], 0.073);
  }
  if (strlen(params->swr_tran)) {
    n = tn;
    copy_autotracer_by_name("swr_transmission", trinfo, ntr, &tn, tr, &master->swr_tran);
    if (strlen(params->swr_regions))
      trf_dataset(params->swr_tran, trinfo, n, params->ntrS, params->atrS, tr, 0.26);
    else
      trn_dataset(params->swr_tran, trinfo, n, params->ntrS, params->atrS, tr, 0.26);
  }
  if (strlen(params->swr_regions)) {
    n = set_autotracer_by_groupkey("swr_regions", trinfo, ntr, &tn, buf);
    master->swreg = tr[n++];
    master->swrms = tr[n++];
    master->attn_mean = tr[n++];
    master->tran_mean = tr[n++];
  }
  if (params->riverflow) {
    set_autotracer_by_name("flow", trinfo, ntr, &tn, tr, &master->riverflow, buf);
    if (params->riverflow == 2)
      set_autotracer_by_name("flow_depth", trinfo, ntr, &tn, tr, &master->riverdepth, buf);
  }
  if (params->tidep)
    set_autotracer_by_name("equitide", trinfo, ntr, &tn, tr, &master->equitide, buf);
  if (params->numbers1 & TPXO)
    set_autotracer_by_name("tpxotide", trinfo, ntr, &tn, tr, &master->tpxotide, buf);
  if (params->numbers1 & TPXOV) {
    n = set_autotracer_by_groupkey("tpxo_vel", trinfo, ntr, &tn, buf);
    master->tpxovelu = tr[n++];
    master->tpxovelv = tr[n++];
  }
  if (params->numbers1 & TPXOT) {
    n = set_autotracer_by_groupkey("tpxo_tran", trinfo, ntr, &tn, buf);
    master->tpxotranu = tr[n++];
    master->tpxotranv = tr[n++];
  }
  if (params->numbers1 & TRAN2D) {
    n = set_autotracer_by_groupkey("transport_2d", trinfo, ntr, &tn, buf);
    master->uat = tr[n++];
    master->vat = tr[n++];
  }
  if (params->decf & DEC_ETA)
    set_autotracer_by_name("decorr_e1", trinfo, ntr, &tn, tr, &master->decv1, buf);
  if (strlen(params->imp2df)) {
    n = tn;
    set_autotracer_by_name("imp2df", trinfo, ntr, &tn, tr, NULL, buf);
    strcpy(trinfo[n].name, params->imp2dn);
    strcpy(trinfo[n].long_name, params->imp2dn);
    strcpy(trinfo[n].units, params->imp2du);
    strcpy(master->trinfo_2d[tn].i_rule, "nn_sibson");
    if (strlen(params->imp2dt))
      sprintf(trinfo[n].data, "[data=%s(t=%s)]", params->imp2df, params->imp2dt);
    else
      sprintf(trinfo[n].data, "[data=%s]", params->imp2df);
  }
  return(tn);
}

/* END set_tracer_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create the 3D tracers in the params tracer_info.       */
/* Defaults for auto runmodes are set in this routine.               */
/* This information is copied to master->trinfo_3d in hd_init().     */
/*-------------------------------------------------------------------*/
void create_tracer_3d(parameters_t *params)   /* Input parameters    */

{
  tracer_info_t *trinfo = params->trinfo_3d;
  char buf[MAXSTRLEN];
  int n, tn = 0;

  tn = set_tracer_3d(params, NULL, params->ntr, params->trinfo_3d, NULL);

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    /* Set up sediment tracers if required                           */
    tn = sediment_autotracer_3d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, trinfo, params->ntr, tn);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    /* Set up ecology tracers if required                            */
    tn = ecology_autotracer_3d(params->prmfd, params->do_eco, params->eco_vars, 
			       params->eco_defs, params->pre_eco,
			       trinfo, params->ntr, tn);
  }
#endif

}

/* END create_tracer_3d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create a scaing netCDF file (spatially and temporally  */
/* dependent) given a list of time series files at a single (x,y,z)  */
/* location, or a profile at a (x,y) location and a forcing data     */
/* file. The scaling is calculated via:                              */
/* scale(x,y,z,t) = time series(x,y,z,t) - forcing(x,y,z,t)          */
/*-------------------------------------------------------------------*/
void calc_scaling(parameters_t *params, /* Input parameters data     */
		  master_t *master,     /* Master data               */
		  char *mapname         /* Name of input data file   */
		  )
{
  FILE *ip = NULL;           /* Input data file pointer              */
  FILE *op = NULL;           /* Raw oputput file                     */
  int i, m, j, k, nf;        /* Counters                             */
  int kk, k1, k2;            /* Interpolation indicies               */
  int cdfid;                 /* netCDF output file handle            */
  char key[MAXSTRLEN];       /* Dummy buffer                         */
  char buf[MAXSTRLEN];       /* Dummy buffer                         */
  char t_units[MAXSTRLEN];   /* Time units                           */
  char v_units[MAXSTRLEN];   /* Units for scaling variable           */
  char v_name[MAXSTRLEN];    /* Scaling variable name                */
  char f_name[MAXSTRLEN];    /* Forcing file name                    */
  char o_name[MAXSTRLEN];    /* Output variable name                 */
  char r_name[MAXSTRLEN];    /* Reference file name (optional)       */
  char b_name[MAXSTRLEN];    /* Bottom file name (optional)          */
  int ntsfiles;              /* Number of time series files          */
  timeseries_t *tsfiles;     /* Time series files                    */
  timeseries_t rforce;       /* Reference file                       */
  timeseries_t bforce;       /* Bottom file                          */
  int nffiles;               /* Number of forcing files              */
  /*timeseries_t tsforce;*/  /* Forcing file                         */
  timeseries_t **tsforce;    /* Forcing files                        */
  int varids[MAXNUMTSFILES]; /* Forcing files variable id's          */
  cstring *tsvars;           /* Timeseries files scaling name        */
  char outfile[MAXSTRLEN];   /* Name of netCDF output file           */
  double start;              /* Start time of output file            */
  double stop;               /* Stop time of output file             */
  double step;               /* Data interval of output file         */
  double t, pt, ta;          /* Time counter                         */
  double zp = -9999.0;       /* No-data value                        */
  double zd;                 /* Cell centered depth                  */
  int ti;                    /* Time counter                         */
  int recdimid;              /* netCDF Record id                     */
  int pointsid;              /* netCDF #timeseries id                */
  int dims[10];              /* netCDF dimension id                  */
  size_t sindex[3];          /* netCDF array                         */
  size_t count[3];           /* netCDF array                         */
  int npoints;               /* Number of timeseries (x,y) locations */
  int  *mask;                /* Mask for (x,y) timeseries locations  */
  int tid;                   /* netCDF time id                       */
  int kid;                   /* netCDF (x,y) location id             */
  int xid;                   /* netCDF longitude id                  */
  int yid;                   /* netCDF latitude id                   */
  int zid;                   /* netCDF depth id                      */
  int etaid;                 /* netcdf variable id                   */
  int *fof;                  /* Flag for bottom cell forcing data    */
  int **map;                 /* ntsfiles to npoints mapping          */
  int *kbot;                 /* Bottom layer of (x,y) locations      */
  int *cs, *cb;              /* Surface / bottom sparse locations    */
  double *lat, *lon;         /* Latitude / longitude of (x,y) locs   */
  double *x, *y, *z;         /* Lat, lon, depth of timeseries files  */
  double **zz;               /* Layer structure to interpolate onto  */
  double *fv, *sc, *scda;    /* Forcing value, scaled value          */
  double refv;               /* Reference scaling value              */
  double zref;               /* Reference scaling depth              */
  double botv;               /* Bottom scaling value                 */
  double minval, maxval;     /* Tracer bounds                        */
  int kref = zp;             /* Reference k level                    */
  int kbof = params->nz;     /* Bottom k level                       */
  int plotif = 1;            /* Print individual files               */
  int daf = 0;               /* Depth average scaling                */
  geometry_t *geom = master->sgrid;

  /*-----------------------------------------------------------------*/
  /* Open the input scaling data file                                */
  if( (ip=fopen(mapname,"r")) == NULL ) {
    hd_warn("calc_scaling: Can't open input scaling data file '%s'\n",mapname);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Read the scaling variable and units                             */
  prm_read_char(ip, "VAR_NAME", v_name);
  prm_read_char(ip, "VAR_UNITS", v_units);
  prm_read_char(ip, "OUT_FILE", outfile);
  prm_read_char(ip, "OUT_NAME", o_name);
  memset(f_name, 0, sizeof(char) * MAXSTRLEN);
  prm_read_char(ip, "FORCING", f_name);
  if (prm_read_char(ip, "DEPTH_AVERAGE", buf))
    daf = is_true(buf);

  /*-----------------------------------------------------------------*/
  /* Set start / stop and step times                                 */
  sprintf(t_units, "%s", params->timeunit);
  tm_scale_to_secs(params->start_time, &start);
  tm_scale_to_secs(params->stop_time, &stop);  
  step = 86400.0;   /* Default of 1 day                              */
  if (prm_read_char(ip, "OUT_DT", buf))
    tm_scale_to_secs(buf, &step);

  /*-----------------------------------------------------------------*/
  /* Read the reference value                                        */
  memset(r_name, 0, sizeof(char) * MAXSTRLEN);
  if (prm_read_char(ip, "REF_VALUE", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int nf = parseline(buf, fields, MAXNUMARGS);
    if (nf <= 0) {
      hd_warn("calc_scaling : REF_VALUE incorrectly specified.\n");
      return;
    }
    if (!(prm_read_double(ip, "REF_DEPTH", &zref))) {
      hd_warn("calc_scaling : REF_DEPTH incorrectly specified.\n");
      return;
    }
    if (zref > 0.0)zref = -zref;
    if (nf == 1) {
      if (strcasecmp(fields[0], "default") == 0)  /* Default value   */
	refv = 0.0;
      else if (sscanf(fields[0], "%lf", &refv) != 1) { /* A number   */
	refv = zp;
	sprintf(r_name,"%s",fields[0]);           /* A filename      */
	read_tsfile(&rforce, r_name, t_units);
      }
    }
    botv = refv;
    kref = face_dist(params->layers, params->nz, zref);
  }

  /*-----------------------------------------------------------------*/
  /* Read the bottom value                                           */
  memset(b_name, 0, sizeof(char) * MAXSTRLEN);
  if (prm_read_char(ip, "BOT_VALUE", buf)) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    int nf = parseline(buf, fields, MAXNUMARGS);
    if (nf <= 0) {
      hd_warn("calc_scaling : BOT_VALUE incorrectly specified.\n");
      return;
    }
    if (nf == 1) {
      if (strcasecmp(fields[0], "default") == 0)  /* Default value   */
	botv = 0.0;
      else if (sscanf(fields[0], "%lf", &botv) != 1) { /* A number   */
	botv = zp;
	sprintf(b_name,"%s",fields[0]);           /* A filename      */
	read_tsfile(&bforce, b_name, t_units);
      }
    }
    kbof = params->nz;
    for (i = 1; i <= geom->v2_t;  i++) {
      j = geom->w2_t[i];
      if (geom->s2k[geom->bot_t[j]] < kbof)
	kbof = geom->s2k[geom->bot_t[j]];
    }
  }

  /*-----------------------------------------------------------------*/	
  /* Get the number of files inputs and allocate memory              */
  prm_read_int(ip,"nfiles",&ntsfiles);
  if (ntsfiles <= 0)
    hd_quit("calc_scaling: No files specified in parameter file.\n");
  tsfiles = (timeseries_t *)malloc(sizeof(timeseries_t) * ntsfiles);
  tsvars = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  memset(tsfiles, 0, sizeof(timeseries_t) * ntsfiles);
  memset(tsvars, 0, sizeof(cstring) * ntsfiles);
  x = d_alloc_1d(ntsfiles);
  y = d_alloc_1d(ntsfiles);
  z = d_alloc_1d(ntsfiles);

  /*-----------------------------------------------------------------*/	
  /* Read the time-series and forcing files                          */
  for (i=0; i<ntsfiles; ++i) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    sprintf(key, "file%d", i);
    prm_read_char(ip, key, buf);
    nf = parseline(buf, fields, MAXNUMARGS);
    if (nf < 3) hd_quit("File %d must include : filename xloc yloc depth var_name.\n", i);
    read_tsfile(&tsfiles[i], strdup(fields[0]), t_units);
    x[i] = atof(fields[1]);
    y[i] = atof(fields[2]);
    if (strcmp(fields[3],"profile")==0) {
      z[i] = zp;
    } else {
      z[i] = atof(fields[3]);
    }
    if (nf > 1)
      strcpy(tsvars[i], fields[4]);
    else
      strcpy(tsvars[i], v_name);
  }

  /*-----------------------------------------------------------------*/	
  /* Open the forcing file                                           */
  if(strlen(f_name)) {
    /*
    read_tsfile(&tsforce, f_name, t_units);
    */
    char files[MAXNUMTSFILES][MAXSTRLEN];
    cstring *filenames;

    nffiles = parseline(f_name, (char **)files, MAXNUMTSFILES);
    filenames = (cstring *)malloc(sizeof(cstring) * nffiles);
    for (i = 0; i < nffiles; ++i)
     strcpy(filenames[i], ((char **)files)[i]);

    tsforce = hd_ts_multifile_read(master, nffiles, filenames);
    if (tsforce == NULL) {
      free(filenames);
      free(tsforce);
      hd_quit("calc_scaling: Can't open forcing file %s.\n", f_name);
    }
    if (hd_ts_multifile_get_index(nffiles, tsforce, filenames, v_name, varids) > 0) {
      for(i = 0; i < nffiles; i++) {
	if (varids[i] < 0)
	  hd_warn("calc_scaling: Can't find variable '%s' in timeseries file %s.\n", 
		  v_name, filenames[i]);
      }
    } else
      hd_quit("calc_scaling: Can't find variable '%s' in any timeseries file %s.\n", 
	      v_name, f_name);
  }

  /*-----------------------------------------------------------------*/	
  /* Count the number of points (i.e. (x,y) locations)               */
  mask = i_alloc_1d(ntsfiles);
  for (i=0; i<ntsfiles; ++i)
    mask[i] = 1;
  for (i=0; i<ntsfiles; ++i) {
    if (!mask[i]) continue;
    for (j=i+1; j<ntsfiles; ++j) {
      if (x[i] == x[j] && y[i] == y[j])
	mask[j] = 0;
    }
  }
  npoints = 0;
  for (i=0; i<ntsfiles; ++i)
    npoints += mask[i];

  lon = d_alloc_1d(npoints);
  lat = d_alloc_1d(npoints);
  zz = d_alloc_2d(params->nz + 1, npoints);
  map = i_alloc_2d(params->nz + 1, npoints);
  cs = i_alloc_1d(npoints);
  cb = i_alloc_1d(npoints);
  kbot = i_alloc_1d(npoints);
  fof = i_alloc_1d(npoints);
  fv = d_alloc_1d(params->nz + 1);
  sc = d_alloc_1d(params->nz + 1);
  if (daf) scda = d_alloc_1d(npoints);
  memset(scda, 0, sizeof(double) * npoints);
  memset(fof, 0, sizeof(int) * npoints);

  /*-----------------------------------------------------------------*/	
  /* Assign the point lat, lon and depth                             */
  j = 0;
  for (i=0; i<ntsfiles; ++i) {
    if (mask[i]) {
      for (k=0; k<params->nz; k++) {
	if (z[i] == zp)
	  zz[j][k] = 0.5 * (params->layers[k] + params->layers[k+1]);
	else
	  zz[j][k] = zp;
      }
      lon[j] = x[i];
      lat[j] = y[i];
      /* Get the bottom k value */
      hd_xyztoindex_m(geom, lon[j], lat[j], 0.0, &kk, &cs[j], &cb[j]);
      kbot[j] = geom->s2k[cb[j]];
      for (m = 0; m < ntsfiles; m++) {
	if (x[m] != lon[j] && y[m] != lat[j]) continue;
 	/* If the mooring depth is below the sea floor, truncate its */
	/* depth to the sea floor.                                   */
	if (z[m] < params->layers[kbot[j]])
	  z[m] = params->layers[kbot[j]];
	/* Set the cell centered depth corresponding to z[m]         */
	kk = face_dist(params->layers, params->nz, z[m]);
	zz[j][kk] = 0.5 * (params->layers[kk] + params->layers[kk+1]);
	map[j][kk] = m;
      }
      /* Set the bottom layer depth */
      kk = kbot[j];
      if (strlen(f_name)) {
	if (zz[j][kk] == zp) fof[j] = 1;
	zz[j][kk] = 0.5 * (params->layers[kk] + params->layers[kk + 1]);
      }
      j++;
    }
  }
  fclose(ip);

  /*-----------------------------------------------------------------*/
  /* Now that all the data has been parameters have been read,       */
  /* Create the netCDF file.                                         */
  if (nc_create(outfile,NC_NOCLOBBER, &cdfid) != NC_NOERR) {
    hd_warn("Couldn't create output dump file %s\n", outfile);
    return;
  }

  /* define dimensions */
  nc_def_dim(cdfid,"record",NC_UNLIMITED,&recdimid);
  nc_def_dim(cdfid,"points", npoints,&pointsid);
  nc_def_dim(cdfid,"k",params->nz,&kid);

  dims[0] = kid;
  nc_def_var(cdfid,"z",NC_DOUBLE,1,dims,&zid);
  addatt(cdfid,zid,"units", "metre");
  addatt(cdfid,zid,"long_name", "Z coordinate");
  addatt(cdfid,zid,"coordinate_type", "Z");

  dims[0] = pointsid;
  nc_def_var(cdfid,"x",NC_DOUBLE,1,dims,&xid);
  addatt(cdfid,xid,"units", "degrees_east");
  addatt(cdfid,xid,"long_name", "Longitude");
  addatt(cdfid,xid,"coordinate_type", "longitude");

  nc_def_var(cdfid,"y",NC_DOUBLE,1,dims,&yid);
  addatt(cdfid,yid,"units", "degrees_north");
  addatt(cdfid,yid,"long_name", "Latitude");
  addatt(cdfid,yid,"coordinate_type", "latitude");

  /* time dependent variables */
  dims[0] = recdimid;
  nc_def_var(cdfid,"t",NC_DOUBLE,1,dims,&tid);
  addatt(cdfid,tid,"units", t_units);
  addatt(cdfid,tid,"long_name", "Time");
  addatt(cdfid,tid,"coordinate_type", "time");

  dims[1] = kid;
  dims[2] = pointsid;
  nc_def_var(cdfid,o_name,NC_DOUBLE,3,dims,&etaid);
  addatt(cdfid,etaid,"units", v_units);
  addatt(cdfid,etaid,"long_name", o_name);
  addatt(cdfid,etaid,"coordinates", "t, z, y, x");

  nc_enddef(cdfid);

  /*-----------------------------------------------------------------*/	
  /* Write the Position information                                  */
  for (i=0; i<npoints; ++i) {
    sindex[0] = i;
    count[0] = 1;
    nc_put_vara_double(cdfid,xid,sindex,count,&lon[i]);
    nc_put_vara_double(cdfid,yid,sindex,count,&lat[i]);
  }
  for (k=0; k<params->nz; ++k) {
    sindex[0] = k;
    count[0] = 1;
    zd = 0.5 * (params->layers[k] + params->layers[k+1]);
    nc_put_vara_double(cdfid,zid,sindex,count,&zd);
  }

  /*-----------------------------------------------------------------*/	
  /* Step through all time, evaluate and write the output value      */
  if (daf) {
    if (endswith(o_name, ".nc")) {
      m = strlen(o_name);
      for (i = 0; i < m-3; i++)
	buf[i] = o_name[i];
      buf[i] = '\0';
    } else
      strcpy(buf, o_name);
    strcat(buf, ".txt");
    op = fopen(buf, "w");
    for (i=0; i<ntsfiles; ++i) {
      fprintf(op, "## file%d %s at (%4.2f %4.2f %2.2f)\n", i, tsfiles[i].name, x[i], y[i], z[i]);
    }
    fprintf(op, "## COLUMNS 3\n");
    fprintf(op, "##\n");
    fprintf(op, "## COLUMN1.name  lon\n");
    fprintf(op, "## COLUMN1.long_name  Longitude\n");
    fprintf(op, "## COLUMN1.units  degrees_east\n");
    fprintf(op, "##\n");
    fprintf(op, "## COLUMN2.name  lat\n");
    fprintf(op, "## COLUMN2.long_name  Latitude\n");
    fprintf(op, "## COLUMN2.units  degrees_north\n");
    fprintf(op, "##\n");
    fprintf(op, "## COLUMN3.name  %s\n", o_name);
    fprintf(op, "## COLUMN3.long_name  %s\n", o_name);
    fprintf(op, "## COLUMN3.units  %s\n", v_units);
    fprintf(op, "##\n");
  }
  t = start;
  ti = 0;
  ta = 0.0;
  while (t <= stop) {
    double val = 0.0, sref = 0.0, bref = 0.0;
    sindex[0] = ti;
    count[0] = 1L;
    nc_put_vara_double(cdfid,tid,sindex,count,&t);
    ta += 1.0;

    for (i=0; i<npoints; ++i) {

      for (k=0; k<params->nz; k++) {

        /*-----------------------------------------------------------*/	
	/* Compute the new value                                     */
	val = 0.0;
	int varid = ts_get_index(&tsfiles[i], tsvars[i]);
	if (strlen(f_name) && k == kbot[i] && fof[i]) {
	  /*
	  varid = ts_get_index(&tsforce, v_name);
	  val = ts_eval_xyz(&tsforce, varid, t, lon[i], lat[i], zz[i][k]);
	  */
	  val = hd_ts_multifile_eval_xyz(nffiles, tsforce, varids,
					t, lon[i], lat[i], zz[i][k]);
	} else if (zz[i][k] != zp) {
	  j = map[i][k];
	  varid = ts_get_index(&tsfiles[j], tsvars[j]);
	  val = ts_eval_xyz(&tsfiles[j], varid, t, lon[i], lat[i], zz[i][k]);
	} else {
	  val = zp;
	}
	fv[k] = val;

	/* Overwrite with reference value if required                */
	if (k == kref) {
	  if (k < kbot[i])kbot[i] = k;
	  if (strlen(r_name)) {
	    varid = ts_get_index(&rforce, v_name);
	    fv[k] = ts_eval_xyz(&rforce, varid, t, lon[i], lat[i], zref);
	  } else
	    fv[k] = refv;
	}
	if (k == kbof) {
	  if (strlen(b_name)) {
	    zd = 0.5 * (params->layers[kbof] + params->layers[kbof+1]);
	    varid = ts_get_index(&bforce, v_name);
	    fv[k] = ts_eval_xyz(&bforce, varid, t, lon[i], lat[i], zd);
	  } else
	    fv[k] = botv;
	}
      }

      /*-------------------------------------------------------------*/
      /* Interpolate                                                 */
      for (k=0; k<kbot[i]; k++)
	sc[k] = 0.0;

      k1 = k2 = kbot[i];
      for (k=kbot[i]; k<params->nz; k++) {
	/*int varid = ts_get_index(&tsforce, v_name);*/
	zd = 0.5 * (params->layers[k] + params->layers[k+1]);
	/*val = ts_eval_xyz(&tsforce, varid, t, lon[i], lat[i], zd);*/
	val = hd_ts_multifile_eval_xyz(nffiles, tsforce, varids,
				       t, lon[i], lat[i], zd);
	if (fv[k] == zp) {
	  for (kk = k; kk < params->nz; kk++)
	    if (fv[kk] != zp) {
	      k2 = kk;
	      break;
	    }
	  if (k1 == k2)
	    fv[k] = fv[k1];
	  else
	    fv[k] = (fv[k2] - fv[k1]) * (double)(k - k1) / 
	      (double)(k2 - k1) + fv[k1];
	}
	if (k == k2)
	  k1 = k2;
	sc[k] = fv[k] -  val;
	if (k == kref) sref = sc[k];
	if (k == kbot[i]) bref = sc[k];
      }
      if (kref != zp) {
	for (k = kbof; k < kref; k++)
	  sc[k] = (sref - fv[kbof]) * (double)(k - kbof) / 
	    (double)(kref - kbof) + fv[kbof];
      } else {
	for (k = kbof; k < kbot[i]; k++)
	  sc[k] = (bref - fv[kbof]) * (double)(k - kbof) / 
	    (double)(kbot[i] - kbof) + fv[kbof];
      }

      /*-------------------------------------------------------------*/
      /* Depth average                                               */
      if (daf) {
	double depth = 0.0;
	val = 0.0;
	for (k=0; k<params->nz; k++) {
	  zd = params->layers[k+1] - params->layers[k];
	  val += sc[k] * zd;
	  depth += zd;
	}
	scda[i] += val / depth;
	for (k=0; k<params->nz; k++)
	  sc[k] = val / depth;
      }

      /*-------------------------------------------------------------*/
      /* Write to file                                               */
      for (k=0; k<params->nz; k++) {
	sindex[0] = ti;
	sindex[1] = k;
	sindex[2] = i;
	count[0] = 1;
	count[1] = 1;
	count[2] = 1;
	nc_put_vara_double(cdfid,etaid,sindex,count,&sc[k]);
      }
    }

    nc_sync(cdfid);

    t+=step;
    ++ti;
  }
  if (daf) {
    for (i=0; i<npoints; ++i)
      fprintf(op, "%f %f %f\n", lon[i], lat[i], scda[i] / ta);
    fclose(op);
    d_free_1d(scda);
  }

  /*-----------------------------------------------------------------*/
  /* Find the valid ranges for the tracer                            */
  minval = -1e30;
  maxval = 1e30;
  for (i = 0; i < master->ntr; i++) {
    tracer_info_t *tr = &master->trinfo_3d[i];
    if (strcmp(v_name, tr->name) == 0) {
      minval = tr->valid_range_wc[0];
      maxval = tr->valid_range_wc[1];
      break;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Write the raw data to file at the times and locations specified */
  /* in the observation files.                                       */
  if (endswith(o_name, ".nc")) {
    m = strlen(o_name);
    for (i = 0; i < m-3; i++)
      buf[i] = o_name[i];
    buf[i] = '\0';
  } else
    strcpy(buf, o_name);
  strcat(buf, ".ts");
  op = fopen(buf, "w");
  for (i=0; i<ntsfiles; ++i) {
    fprintf(op, "## file%d %s at (%4.2f %4.2f %2.2f)\n", i, tsfiles[i].name, x[i], y[i], z[i]);
  }
  fprintf(op, "## COLUMNS 4\n");
  fprintf(op, "##\n");
  fprintf(op, "## COLUMN1.name  Time\n");
  fprintf(op, "## COLUMN1.long_name  Time\n");
  fprintf(op, "## COLUMN1.units  %s\n", master->output_tunit);
  fprintf(op, "##\n");
  m = 1;
  fprintf(op, "## COLUMN%1.1d.name  file_no\n", ++m);
  fprintf(op, "## COLUMN%1.1d.long_name File number\n",m);
  fprintf(op, "## COLUMN%1.1d.units  \n",m);
  fprintf(op, "##\n");
  fprintf(op, "## COLUMN%1.1d.name  observation\n", ++m);
  fprintf(op, "## COLUMN%1.1d.long_name Measured data value\n",m);
  fprintf(op, "## COLUMN%1.1d.units %s\n",m, v_units);
  fprintf(op, "##\n");
  fprintf(op, "## COLUMN%1.1d.name  forcing\n", ++m);
  fprintf(op, "## COLUMN%1.1d.long_name Forcing data value\n",m);
  fprintf(op, "## COLUMN%1.1d.units %s\n",m, v_units);
  fprintf(op, "##\n");

  for (i=0; i<ntsfiles; ++i) {
    int varid;
    double v1, v2;
    timeseries_t *ts;
    FILE *fp;
    ts = &tsfiles[i];

    if (plotif) {
      strcpy(key, ts->name);
      for (m = strlen(key)-1; m >= 0; m--)
	if (key[m] == '/') {
	  m++;
	  break;
	}
      for (k = m; k < strlen(key); k++)
	buf[k-m] = key[k];
      buf[k-m] = '\0';
      sprintf(key, "SC_%s", buf);
      fp = fopen(key, "w");
      fprintf(fp, "## COLUMNS 3\n");
      fprintf(fp, "##\n");
      fprintf(fp, "## COLUMN1.name  Time\n");
      fprintf(fp, "## COLUMN1.long_name  Time\n");
      fprintf(fp, "## COLUMN1.units  %s\n", master->output_tunit);
      fprintf(fp, "##\n");
      fprintf(fp, "## COLUMN2.name  observation\n");
      fprintf(fp, "## COLUMN2.long_name Measured data value\n");
      fprintf(fp, "## COLUMN2.units %s\n", v_units);
      fprintf(fp, "##\n");
      fprintf(fp, "## COLUMN3.name  forcing\n");
      fprintf(fp, "## COLUMN3.long_name Forcing data value\n");
      fprintf(fp, "## COLUMN3.units %s\n", v_units);
      fprintf(fp, "##\n");
    }
    pt = ts->t[0];
    for (m = 0; m < ts->nt; m++) {
      t = ts->t[m];
      if (t < start || t > stop) continue;
      if (m > 0 && t - pt < params->grid_dt) continue;
      varid = ts_get_index(&tsfiles[i], tsvars[i]);
      v1 = ts_eval_xyz(&tsfiles[i], varid, t, x[i], y[i], z[i]);
      /*
      varid = ts_get_index(&tsforce, v_name);
      v2 = ts_eval_xyz(&tsforce, varid, t, x[i], y[i], z[i]);
      */
      v2 = hd_ts_multifile_eval_xyz(nffiles, tsforce, varids,
				    t, x[i], y[i], z[i]);
      pt = t;
      tm_change_time_units(master->timeunit, master->output_tunit, &t, 1);
      v2 = (v2 < minval) ? NaN : v2;
      v2 = (v2 > maxval) ? NaN : v2;
      fprintf(op, "%f %d %f %f\n",t, i, v1, v2);
      if (plotif) fprintf(fp, "%f %f %f\n",t, v1, v2);
    }
    if (plotif) fclose(fp);
  }
  fclose(op);

  nc_close(cdfid);
  /*ts_free(&tsforce);*/
  for (i=0; i<nffiles; ++i)
    hd_ts_free(master, tsforce[i]);
  if (strlen(r_name))ts_free(&rforce);
  if (strlen(b_name))ts_free(&bforce);
  for (i = 0; i < ntsfiles; ++i)
    ts_free(&tsfiles[i]);

  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(z);
  d_free_1d(lon);
  d_free_1d(lat);
  d_free_1d(fv);
  d_free_1d(sc);
  d_free_2d(zz);
  i_free_2d(map);
  i_free_1d(cs);
  i_free_1d(cb);
  i_free_1d(kbot);
  i_free_1d(fof);

}

/* END calc_scaling()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create a scaing tracer from glider data.               */
/* scale(x,y,z,t) = time series(x,y,z,t) - forcing(x,y,z,t)          */
/* The scaling can be done for individual glider missions, applied   */
/* to different forcing data. The spatial distribution can be        */
/* uniform from a mean scaling profile, or interpolated using        */
/* grid_specs interpolation.                                         */
/*-------------------------------------------------------------------*/
void glider_scaling(parameters_t *params, /* Input parameters data   */
		    master_t *master,     /* Master data             */
		    char *mapname,        /* Name of input data file */
		    int tm                /* Tracer number           */
		    )
{
  geometry_t *geom = master->geom;
  FILE *ip = NULL;           /* Input data file pointer              */
  int i, m, j, jj, k, nf;    /* Counters                             */
  int c, cc, cs, cb, dj, c2; /* Counters                             */
  int **ck, *nck;            /* Layer maps                           */
  char i_rule[MAXSTRLEN];    /* Interpolation rule                   */
  char key[MAXSTRLEN];       /* Dummy buffer                         */
  char buf[MAXSTRLEN];       /* Dummy buffer                         */
  char t_units[MAXSTRLEN];   /* Time units                           */
  char v_units[MAXSTRLEN];   /* Units for scaling variable           */
  char v_name[MAXSTRLEN];    /* Scaling variable name                */
  int ntsfiles;              /* Number of time series files          */
  timeseries_t *tsfiles;     /* Time series files                    */
  int nffiles;               /* Number of forcing files              */
  timeseries_t **tsforce;    /* Forcing files                        */
  int varids[MAXNUMTSFILES]; /* Forcing files variable id's          */
  int dot[MAXNUMTSFILES];    /* Forcing files time flag              */
  double fstart[MAXNUMTSFILES]; /* Forcing files start time          */
  double fstop[MAXNUMTSFILES]; /* Forcing files stop time            */
  cstring *tsvars;           /* Timeseries files scaling name        */
  cstring *xname;            /* x coordinate name                    */
  cstring *yname;            /* y coordinate name                    */
  cstring *zname;            /* z coordinate name                    */
  cstring *fname;            /* z coordinate name                    */
  double start;              /* Start time of output file            */
  double stop;               /* Stop time of output file             */
  double step;               /* Data interval of output file         */
  double t;                  /* Time counter                         */
  int npoints;               /* Number of timeseries (x,y) locations */
  int  *points;              /* Points in each layer                 */
  double **lat, **lon;       /* Latitude / longitude of (x,y) locs   */
  double x, y, z;            /* Geographic location                  */
  double fv, val, **sc;      /* Forcing value, scaled value          */
  int intf = 0;              /* Interpolation method flag            */
  GRID_SPECS **gs;           /* Interpolation data structure         */
  int nz = geom->nz;         /* Vertical dimension                   */
  double vprof[nz];          /* Mean vertical profile                */
  tracer_info_t *tr = &master->trinfo_3d[tm]; /* Tracer info         */
  int found = 0;             /* Found some valid data                */
  int verbose = 0;

  /*-----------------------------------------------------------------*/
  /* Open the input scaling data file                                */
  if( (ip=fopen(mapname,"r")) == NULL ) {
    hd_warn("glider_scaling: Can't open input scaling data file '%s'\n",mapname);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Read the scaling variable and units                             */
  sprintf(t_units, "%s", params->timeunit);
  prm_read_char(ip, "VAR_NAME", v_name);
  prm_read_char(ip, "VAR_UNITS", v_units);
  prm_get_time_in_secs(ip, "START", &start);
  prm_get_time_in_secs(ip, "STOP", &stop);
  prm_get_time_in_secs(ip, "STEP", &step);
  if (prm_read_char(ip, "METHOD", i_rule)) {
    if (strcmp(i_rule, "profile") == 0) 
      intf = 0;
    else if (strcmp(i_rule, "cascade") == 0) {
      intf = 1;
      if (params->us_type & US_IJ)
	hd_quit("Glider scaling with 'cascade' filling not optional with structured grids.\n");
    } else
      intf = 2;
  }

  /*-----------------------------------------------------------------*/	
  /* Get the number of files inputs and allocate memory              */
  prm_read_int(ip,"nfiles",&ntsfiles);
  if (ntsfiles <= 0)
    hd_quit("calc_scaling: No files specified in parameter file.\n");
  tsfiles = (timeseries_t *)malloc(sizeof(timeseries_t) * ntsfiles);
  tsforce = (timeseries_t **)malloc(sizeof(timeseries_t) * ntsfiles);
  tsvars = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  xname = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  yname = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  zname = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  fname = (cstring *)malloc(sizeof(cstring) * ntsfiles);
  memset(tsfiles, 0, sizeof(timeseries_t) * ntsfiles);
  memset(tsforce, 0, sizeof(timeseries_t) * ntsfiles);
  memset(tsvars, 0, sizeof(cstring) * ntsfiles);
  memset(xname, 0, sizeof(cstring) * ntsfiles);
  memset(yname, 0, sizeof(cstring) * ntsfiles);
  memset(zname, 0, sizeof(cstring) * ntsfiles);
  memset(fname, 0, sizeof(cstring) * ntsfiles);
  memset(vprof, 0, nz * sizeof(double));
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    master->tr_wc[tm][c] = NOTVALID;
  }

  /*-----------------------------------------------------------------*/	
  /* Read the time-series and forcing files                          */
  for (i=0; i<ntsfiles; ++i) {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    sprintf(key, "file%d", i);
    prm_read_char(ip, key, buf);
    nf = parseline(buf, fields, MAXNUMARGS);
    if (nf < 7) hd_quit("File %d must include : filename xloc yloc depth var_name forcing_file forcing_name.\n", i);
    read_tsfile(&tsfiles[i], strdup(fields[0]), t_units);
    strcpy(xname[i], fields[1]);
    strcpy(yname[i], fields[2]);
    strcpy(zname[i], fields[3]);
    if (nf > 1)
      strcpy(tsvars[i], fields[4]);
    else
      strcpy(tsvars[i], v_name);
    strcpy(fname[i], fields[5]);
    tsforce[i] = hd_ts_read(master, fname[i], 0);
    strcpy(v_name, fields[6]);
    if ((varids[i] = ts_get_index(tsforce[i], v_name)) < 0) {
      hd_warn("glider_scaling: Can't find variable '%s' in timeseries file %s.\n", 
	      v_name, fname[i]);
    }
    dot[i] = 0;
    fstart[i] = start;
    fstop[i] = stop;
    if (nf == 10) {
      sprintf(buf, "%s %s", fields[7], fields[9]);
      tm_scale_to_secs(buf, &fstart[i]);
      sprintf(buf, "%s %s", fields[8], fields[9]);
      tm_scale_to_secs(buf, &fstop[i]);
      dot[i] = 1;
    }
    if (verbose) {
      printf("file%d data: %s\n", i, fields[0]);
      printf("file%d forcing: %s\n", i, fname[i]);
    }
  }

  /*-----------------------------------------------------------------*/	
  /* Count the number of points (i.e. (x,y) locations)               */
  npoints  = ntsfiles * ((stop - start) / step + 1);
  points = i_alloc_1d(nz + 1);
  memset(points, 0, (nz+1) * sizeof(int));
  lon = d_alloc_2d(npoints, nz);
  lat = d_alloc_2d(npoints, nz);
  sc = d_alloc_2d(npoints, nz);

  /*-----------------------------------------------------------------*/	
  /* Step through all time, evaluate and write the output value      */
  t = start;
  while (t <= stop) {
    int varid, id;

    for (i = 0; i < ntsfiles; i++) {
      timeseries_t *ts = &tsfiles[i];
      timeseries_t *tsf = tsforce[i];
      double t_g = t;
      double t_f = t;

      /* Continue if time is out of range                            */
      tm_change_time_units(master->timeunit, ts->t_units, &t_g, 1);
      m = ts_has_time(ts, t_g);
      if (!m) continue;
      tm_change_time_units(master->timeunit, tsf->t_units, &t_f, 1);
      m = ts_has_time(tsf, t_f);
      if (!m) continue;
      if (dot[i] && t < fstart[i]) continue;
      if (dot[i] && t > fstop[i]) continue;

      varid = ts_get_index(ts, zname[i]);
      z = ts_eval(ts, varid, t);
      for (k = nz - 1; k >= 0; k--) {
	if (params->layers[k] < z) break;
      }
      if (k < 0) continue;
      varid = ts_get_index(ts, xname[i]);
      x = ts_eval(ts, varid, t);
      varid = ts_get_index(ts, yname[i]);
      y = ts_eval(ts, varid, t);
      varid = ts_get_index(ts, tsvars[i]);
      fv = ts_eval(ts, varid, t);
      val = ts_eval_xyz(tsf, varids[i], t, x, y, z);
      if (isnan(fv) || isnan(val)) continue;

      /* Clip                                                        */
      fv = min(max(tr->valid_range_wc[0], fv), tr->valid_range_wc[1]);
      val = min(max(tr->valid_range_wc[0], val), tr->valid_range_wc[1]);

      /* Use all glider data to make a 1D mean profile               */
      if (intf == 0) {
	if (vprof[k] == 0.0)
	  vprof[k] = fv -  val;
	else
	  vprof[k] = 0.5 * (sc[k][dj] + (fv -  val));
	continue;
      }

      /* Check for duplicates                                        */
      dj = -1;
      for (jj = 0; jj < points[k]; jj++) {
	if (x == lon[k][jj] && y == lat[k][jj]) {
	  dj = jj;
	  break;
	}
      }

      /* Save position and scaling for interpolation                 */
      if (dj >= 0) {
	j = dj;
	sc[k][dj] = 0.5 * (sc[k][dj] + (fv -  val));
	if (verbose) printf("%f : %f %f %f (%d) : %f %f %f dup\n",t, x, y, z, k,
			    fv,val,sc[k][dj]);

      } else {
	j = points[k];
	lon[k][j] = x;
	lat[k][j] = y;
	sc[k][j] = fv -  val;
	points[k] += 1;
	if (verbose) printf("%f : %f %f %f (%d) : %f %f %f\n",t, x, y, z, k,
			    fv,val,sc[k][j]);
      }

      /* Nearest neighbour spatial interpolation                     */
      if (intf == 1) {
	if ((c = hd_grid_xyztoc_m(master, lon[k][j], lat[k][j], z)) != -1)
	  master->tr_wc[tm][c] = sc[k][j];
      }
      found = 1;
    }
    t += step;
  }
  if (!found) return;

  /* Interpolate onto the scaling tracer                             */
  if (intf == 0) {
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      k = geom->s2k[c];
      master->tr_wc[tm][c] = vprof[k];
    }
  }
  if (intf == 1) {
    int *mask = i_alloc_1d(geom->szc);
    nck = i_alloc_1d(geom->nz);
    ck = i_alloc_2d(geom->szcS, geom->nz);
    memset(nck, 0, geom->nz * sizeof(int));
    memset(mask, 0, geom->sgsiz * sizeof(int));
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      k = geom->s2k[c];
      ck[k][nck[k]++] = c;
      if ( master->tr_wc[tm][c] != NOTVALID) mask[c] = 1;
    }

    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      c2 = geom->m2d[c];
      k = geom->s2k[c];
      if (master->tr_wc[tm][c] == NOTVALID) {
	x = geom->cellx[c2];
	y = geom->celly[c2];
	if (points[k]) {
	  double dist, x1, y1;
	  double dm = HUGE;
	  for (i = 0; i < nck[k]; i++) {
	    cs = ck[k][i];
	    if (cs != c && mask[cs]) {
	      x1 = geom->cellx[geom->m2d[cs]] - x;
	      y1 = geom->celly[geom->m2d[cs]] - y;
	      dist = sqrt(x1 * x1 + y1 * y1);
	      if (dist < dm) {
		master->tr_wc[tm][c] = master->tr_wc[tm][cs];
		dm = dist;
	      }
	    }
	  }
	}
      }
    }
    i_free_1d(nck);
    i_free_2d(ck);
    i_free_1d(mask);
  }
 
  if (intf == 2) {
    gs = (GRID_SPECS **)calloc(nz+1, sizeof(GRID_SPECS *));
    for (k = nz - 1; k >= 0; k--) {
      gs[k] = NULL;
      if (points[k] > 2)
	gs[k] = grid_interp_init(lon[k], lat[k], sc[k], points[k], i_rule);
    }
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      master->tr_wc[tm][c] = NOTVALID;
      
      x = geom->cellx[cs];
      y = geom->celly[cs];
      k = geom->s2k[c];

      if (gs[k])
	master->tr_wc[tm][c] = grid_interp_on_point(gs[k], x, y);
    }

    for (k = nz - 1; k >= 0; k--) {
      if (gs[k])
	grid_specs_destroy(gs[k]);
    }
  }

  /* Fill any not valid layers                                       */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    double v1, v2;
    int i1 = 0, i2 = 0;
    val = 0.0;
    c = cs = geom->w3_t[cc];
    if (master->tr_wc[tm][c] == NOTVALID) {
      while ((v1 = master->tr_wc[tm][c]) == NOTVALID)
	c = geom->zp1[c];
      if (c == geom->zp1[c] && master->tr_wc[tm][c] == NOTVALID) i1 = 1;
      c = cs;
      while ((v2 = master->tr_wc[tm][c]) == NOTVALID)
	c = geom->zm1[c];
      if (c == geom->zm1[c]) i2 = 1;
      if (!i1) val += v1;
      if (!i2) val += v2;
      if (!i1 && !i2) val *= 0.5;
      master->tr_wc[tm][cs] = val;
    }
  }

  for (i = 0; i < ntsfiles; ++i) {
    ts_free(&tsfiles[i]);
    hd_ts_free(master, tsforce[i]);
  }

  d_free_2d(lon);
  d_free_2d(lat);
  d_free_2d(sc);

}

/* END glider_scaling()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Populates an array with values by region.                         */
/*-------------------------------------------------------------------*/
int value_init_regions(master_t *master, char *dname, double *tr, int mode)
{
  geometry_t *geom = master->sgrid;
  char *files[MAXSTRLEN * MAXNUMARGS];
  int n, nf, nr, rgn, cc, c;
  double *regionid, rgv;
  nf = parseline(dname, files, MAXNUMARGS);
  if(strcmp(files[0], "region") == 0) {
    regionid = d_alloc_1d(geom->sgsiz);
    nr = read_regioni(master, files[1], regionid);
    if (nr) {
      for (n = 2; n < nf; n++) {
	sscanf(files[n], "%d:%lf", &rgn, &rgv);
	if (mode == 3) {
	  for (cc = 1; cc <= geom->b3_t; cc++) {
	    c = geom->w3_t[cc];
	    if (rgn == (int)regionid[c])
	      tr[c] = rgv;
	  }
	}
	if (mode == 2) {
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c])
	      tr[c] = rgv;
	  }     
	}
      }
      d_free_1d(regionid);
      return(1);
    }
    d_free_1d(regionid);
    return(0);
  }
  return(0);
}

/* END value_init_regions()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Populates a 2D array with a regionalisation.                      */
/* Note that the input to read_regioni() is a 3D array, so we copy   */
/* the surface layer into the returning array.                       */
/*-------------------------------------------------------------------*/
int set_regions(master_t *master, char *dname, double *regionid)
{
  geometry_t *geom = master->geom;
  char *files[MAXSTRLEN * MAXNUMARGS];
  int cc, c, nf, nr;
  double *reg;

  reg = d_alloc_1d(geom->szc);
  nf = parseline(dname, files, MAXNUMARGS);
  if(strcmp(files[0], "region") == 0) {
    nr = read_regioni(master, files[1], reg);
    if (nr) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	regionid[c] = reg[c];
      }
      d_free_1d(reg);
      return(1);
    } else
      return(0);
  }
  return(0);
}

/* END get_regions()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise a 3d array from a number, array or file     */
/* input.                                                            */
/*-------------------------------------------------------------------*/
void value_init_3d(master_t *master,     /* Master data              */
		   double *ret,          /* 3D array                 */
		   FILE *fp,             /* File pointer             */
		   char *fname,          /* File name                */
		   char *vname,          /* Variable name            */
		   char *tag,            /* Input file tag           */
		   double fill,          /* Fill value               */
		   char *i_rule          /* Interpolation type       */
		   )
{
  geometry_t *geom = master->sgrid;
  int nce1 = geom->nce1;
  int nce2 = geom->nce2;
  int nz = geom->nz;
  int c, cc, cs;
  int i, j, k, n;
  int number, nvals;
  char buf[MAXSTRLEN], key[MAXSTRLEN], sgn[MAXSTRLEN], line[MAXSTRLEN];
  double val, v1, v2, *d1 = NULL, ***d2;
  int *vec, nvec;
  int size;

  size = geom->sgsiz;
  nvec = geom->b3_t;
  vec = geom->w3_t;

  /*-----------------------------------------------------------------*/
  /* Initialise to the fill value                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    ret[c] = fill;
  }

  /*-----------------------------------------------------------------*/
  /* First check if this variable is contained in tracerdata         */
  if (strlen(master->tracerdata)) {
    timeseries_t *ts = NULL;
    int id;
    ts = hd_ts_read(master, master->tracerdata, 0);
    id = ts_get_index(ts, fv_get_varname(master->tracerdata, vname, buf));
    if (id >= 0) {
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	cs = geom->m2d[c];
	ret[c] = ts_eval_xyz(ts, id, master->t,
			     geom->cellx[cs],
			     geom->celly[cs],
			     geom->cellz[c]);
      }
      if (master->trfilter & TRF_FILL3D) {
	/* 
	 * Note: For hydro autotracers tracer_find_index does not work
	 */
	if ((n = tracer_find_index(vname, master->ntr, master->trinfo_3d)) >= 0)
	  tracer_fill(master, ts, fname, ret, geom->sgsiz, vec, nvec, &master->trinfo_3d[n], fill);
      }
      hd_ts_free(master, ts);
      return;
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }

  if (!strlen(fname)) {
    hd_warn("value_init_3d: No 3D information provided to populate array %s; using defaults\n", vname);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check for regionalized formats                                  */
  if (i_rule && strlen(i_rule))
    sprintf(buf, "[var=%s] [i_rule=%s] %s", vname, i_rule, fname);
  else
    sprintf(buf, "[var=%s] %s\n", vname, fname);
  if ((n = set_variable(master, buf, ret, NULL))) {
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check if it is a region file                                    */
  strcpy(buf, fname);
  if (value_init_regions(master, buf, ret, 3))
    return;

  /*-----------------------------------------------------------------*/
  /* Check if it is a profile specification                          */
  sprintf(sgn, "%c", '\0');
  strcpy(line, fname);
  if (sscanf(line, "%s %s %lf %lf %s", buf, key, &v1, &v2, sgn) == 5) {
    /* dens_scale fits a value to the surface and bottom and         */
    /* a profile by incrementing this value by a scaled vertical     */
    /* density difference.                                           */ 
    if (strcmp(buf, "dens_scale") == 0) {
      dens_scale_mid(master, geom, ret, vname, key, v1, v2, sgn);
    }
    n = tracer_find_index(vname, master->ntr, master->trinfo_3d);
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      ret[c] = max(ret[c], master->trinfo_3d[n].valid_range_wc[0]);
      ret[c] = min(ret[c], master->trinfo_3d[n].valid_range_wc[1]);
    }
    return;
  }

  sprintf(sgn, "%c", '\0');
  strcpy(line, fname);
  if (sscanf(line, "%s %s %lf", buf, key, &v2) == 3 ||
      sscanf(line, "%s %s %lf %s", buf, key, &v2, sgn) == 4) {
    int mf = 0;
    if (v1 = atof(key)) mf = 1;
    /* dens_profile normalizes the density profile to top and       */
    /* densities and creates a new profile by fitting top and       */
    /* bottom tracer endpoints.                                     */
    if (strcmp(buf, "dens_profile") == 0) {
      dens_profile(master, geom, ret, v1, v2);
    }
    /* dens_grad normalizes the density profile to the sum of       */
    /* density gradients and creates a new profile by normalizing   */
    /* tracer to the sum of gradients to a particular depth.        */
    if (strcmp(buf, "dens_grad") == 0) {
      dens_grad(master, geom, ret, v1, v2);
    }
    /* dens_scale fits a value to the surface and bottom and        */
    /* a profile by incrementing this value by a scaled vertical    */
    /* density difference.                                          */ 
    if (strcmp(buf, "dens_scale") == 0) {
      if (mf)
	dens_scale(master, geom, ret, v1, v2, sgn);
      else
	dens_scale_file(master, geom, ret, vname, key, v2, sgn);
    }
    n = tracer_find_index(vname, master->ntr, master->trinfo_3d);
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      ret[c] = max(ret[c], master->trinfo_3d[n].valid_range_wc[0]);
      ret[c] = min(ret[c], master->trinfo_3d[n].valid_range_wc[1]);
    }
    return;
  }

  strcpy(line, fname);
  if (sscanf(line, "%s %s", buf, key) == 2) {
    int mf = 0;
    if (v1 = atof(key)) mf = 1;
    /* dens_grad normalizes the density profile to the sum of       */
    /* density gradients and creates a new profile by normalizing   */
    /* tracer to the sum of gradients to a particular depth.        */
    if (strcmp(buf, "dens_grad") == 0) {
      dens_grad_file(master, geom, ret, vname, key);
    }
    n = tracer_find_index(vname, master->ntr, master->trinfo_3d);
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      ret[c] = max(ret[c], master->trinfo_3d[n].valid_range_wc[0]);
      ret[c] = min(ret[c], master->trinfo_3d[n].valid_range_wc[1]);
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check if it is a list of values in the parameter file           */
  if (sscanf(fname, "%lf", &val) == 1) {
    val = atof(fname);
    /* Check if it is a list of values in the parameter file         */
    if ((int)val == nce1 * nce2 * nz) {
      if (prm_read_darray(fp, tag, &d1, &nvals) > 0) {
	if (nvals == nce1 * nce2 * nz) {
	  double ***d2 = d_alloc_3d(nce1, nce2, nz);
	  n = 0;
	  for (k = 0; k < nz; k++)
	    for (j = 0; j < nce2; j++)
	      for (i = 0; i < nce1; i++)
		d2[k][j][i] = d1[n++];
	  c2s_3d(geom, ret, d2, nce1, nce2, nz);
	  d_free_3d(d2);
	  d_free_1d(d1);
	} else {
	  hd_warn("value_init_3d: incorrect number specified for %s variable.\n", vname);
	}
      }
    } else {
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	ret[c] = val;
      }
    }
  } else {
    /* Not a number nor so assume it is a file                       */
    timeseries_t *ts;
    int id; 

    if (i_rule && strlen(i_rule))
      ts = hd_ts_read_us(master, fname, 0, i_rule);
    else
      ts = hd_ts_read(master, fname, 0);

    /* Infer if this is a plain ascii file (as opposed to a          */
    /* timeseries file) by the absence of the allocation of the time */
    /* field. This interpolation is depth invariant.                 */
    if (ts->t == NULL) {    /* Plain file                            */	
      GRID_SPECS *gs = NULL;
      double *x = NULL, *y = NULL, *z = NULL;
      int ii;
	  
      /* Find and assign the needed variables                        */
      for (ii=0; ii<ts->df->nv; ii++) {
	if (strcmp("lon", ts->df->variables[ii].name) == 0)
	  x = ts->df->variables[ii].data;
	else if (strcmp("lat", ts->df->variables[ii].name) == 0)
	  y = ts->df->variables[ii].data;
	else if (strcmp(vname, ts->df->variables[ii].name) == 0)
	  z = ts->df->variables[ii].data;
      }

      /* Error checking                                              */
      if (x == NULL)
	hd_quit("val_init: 'lon' column not found in '%s'\n", fname);
      if (y == NULL)
	hd_quit("val_init: 'lat' column not found in '%s'\n", fname);
      if (z == NULL)
	hd_quit("val_init: '%s' column not found in '%s'\n", vname, fname);
      
      /* Initiliase the grid specs struct                            */
      gs = grid_interp_init(x, y, z, ts->df->dimensions[0].size, i_rule);
	  
      /* Do the interpolation                                        */
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	cs = geom->m2d[c];
	ret[c] = grid_interp_on_point(gs, geom->cellx[cs], geom->celly[cs]);
	    
	/* Check for nan's                                           */
	if (isnan(ret[c])) ret[c] = fill;
      }
	  
      /* Cleanup                                                     */
      grid_specs_destroy(gs);
    } else {
      /* Sanity check to see if the variable is in this datafile     */
      id = ts_get_index(ts, fv_get_varname(fname, vname, buf));

      if (id < 0)
	hd_warn("val_init_3d: The file '%s' does not contain the tracer '%s'\n", fname, vname);
      else {
	/* Use a sparse interpolation for 'standard' files           */
	/*if (value_init_sparse3d(master, ret, fname, vname, i_rule,
	  schedule->start_time, NULL, ts) == 1) {*/
	  /* Otherwise interpolate as a gridded file                 */
	  for (cc = 1; cc <= nvec; cc++) {
	    c = vec[cc];
	    cs = geom->m2d[c];
	    ret[c] = ts_eval_xyz(ts, id, master->t,
				 geom->cellx[cs],
				 geom->celly[cs],
				 geom->cellz[c]);
	  }
	
      }
      if (master->trfilter & TRF_FILL3D) {
	if ((n = tracer_find_index(vname, master->ntr, master->trinfo_3d)) >= 0)
	  tracer_fill(master, ts, fname, ret, geom->sgsiz, vec, nvec, &master->trinfo_3d[n], fill);
	else
	  tracer_fill(master, ts, fname, ret, geom->sgsiz, vec, nvec, NULL, fill);
      }
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }

  /* Limits */
  /*
  n = tracer_find_index(vname, master->ntr, master->trinfo_3d);
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    if (isnan(ret[c])) ret[c] = fill;
    ret[c] = max(ret[c], master->trinfo_3d[n].valid_range_wc[0]);
    ret[c] = min(ret[c], master->trinfo_3d[n].valid_range_wc[1]);
  }
  */
}

/* END value_init_3d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Valid bathymetry netCDF dimension names                           */
static char *in_dims[5][8] = {
  {"botz", "i_centre", "j_centre", "k_centre", "x_centre", "y_centre", "z_centre", "t"},
  {"botz", "i", "j", "k", "longitude", "latitude", "zc", "time"},
  {"height", "longitude", "latitude", "zc", "lon", "lat", "zc", "time"},
  {"height", "ni", "nj", "nk", "x", "y", "z", "t"},
  {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}
};


/*-------------------------------------------------------------------*/
/* Sets vales of a tracer. Usage:                                    */
/* (var=temp) (i_rule=linear) (region=region.bnc) (0:data0) (1:val1) */
/*-------------------------------------------------------------------*/
int set_variable(master_t *master, char *tag, double *ret, double *tin)
{
  geometry_t *geom = master->geom;
  FILE *fp = master->prmfd;
  char sname[MAXSTRLEN];
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char vname[MAXSTRLEN];
  char rfile[MAXSTRLEN];
  char dfile[MAXSTRLEN];
  char i_rule[MAXSTRLEN];
  char buf1[MAXSTRLEN], key1[MAXSTRLEN], iname[MAXSTRLEN];
  char *files[MAXSTRLEN * MAXNUMARGS];
  int n, m, nf, nr, rgn, cc, c, tn, mode;
  int rf = 0;
  int nv, *vec;
  double *regionid, rgv;
  double *tr;
  double t, fill;

  /* Read the variable specification                                 */
  nf = parseline(tag, files, MAXNUMARGS);

  /* Find the variable                                        */
  sprintf(vname, "%c", '\0');
  for (n = 0; n < nf; n++) {
    strcpy(key1, "var=");
    if (find_token(files[n], key1, buf, ']')) {
      strcpy(vname, buf);
      break;
    }
  }
  if (strlen(vname) == 0) return(0);
  if (ret != NULL)
    tr = ret;
  else
    tr = NULL;
  if ((tn = tracer_find_index(vname, master->ntr, master->trinfo_3d)) >= 0) {
    mode = 3;
    nv = geom->b3_t;
    vec = geom->w3_t;
    if (tr == NULL) tr = master->tr_wc[tn];
  }
  if ((tn = tracer_find_index(vname, master->ntrS, master->trinfo_2d)) >= 0) {
    mode = 2;
    nv = geom->b2_t;
    vec = geom->w2_t;
    if (tr == NULL) tr = master->tr_wcS[tn];
  }
  if (tr == NULL) return(0);

  /* Get the interpolation rule                               */
  sprintf(i_rule, "%c", '\0');
  for (n = 0; n < nf; n++) {
    if (find_token(files[n], "i_rule=", buf, ']')) {
      strcpy(i_rule, buf);
      break;
    }
  }

  /* Check for data                                           */
  for (n = 0; n < nf; n++) {
    if (find_token(files[n], "data=", buf, ']')) {
      strcpy(dfile, buf);
      rf = 0;
      sscanf(dfile, "%lf", &fill);
      strcpy(iname, vname);
      if (tin == NULL)  {
	t = NOTVALID;
	t = schedule->start_time;
      } else {
	t =  *tin;
      }
      decode_indata(buf, buf1, iname, &t);
      /* Check for variable substitution                      
      sprintf(buf1, "(%s%s)", key, buf);
      sprintf(key1, "(%s=", vname);
      if (find_token(buf1, key1, iname, ')')) {
	find_token(buf1, key, buf, '(');
      }*/
      if (strlen(i_rule)) {
	/*interp_data_us(master, buf1, iname, tr, mode, i_rule, NULL, t);*/
	if (mode == 2)
	  value_init_sparse2d(master, tr, buf1, iname, i_rule, t, NULL, NULL);
	if (mode == 3)
	  value_init_sparse3d(master, tr, buf1, iname, i_rule, t, NULL, NULL);
      } else
	interp_data_s(master, buf1, iname, tr, mode, NULL, t);
      return(1);
    }
  }

  /* Check for regionalization                                */
  for (n = 0; n < nf; n++) {
    if (find_token(files[n], "region=", buf, ']')) {
      int nr, nn;
      double *regionid = d_alloc_1d(geom->sgsiz);
      int *mask = i_alloc_1d(geom->sgsiz);
      double t;
      strcpy(rfile, buf);
      rf = 1;
      nr = read_regioni(master, rfile, regionid);
      /* Loop through regions                                 */
      for (m = 0; m < nr; m++) {
	sprintf(key,"%d:", m);
	/* Find the specification for this region             */
	for (nn = 0; nn < nf; nn++) {
	  if (find_token(files[nn], key, buf, ']')) {
	    /* Region value is a number                       */
	    if (sscanf(buf, "%lf", &fill) == 1) {
	      fill = atof(buf);
	      for (cc = 1; cc <= nv; cc++) {
		c = vec[cc];
		if (m == (int)regionid[c])
		  tr[c] = fill;
	      }
	      continue;
	    } else {
	      memset(mask, 0, geom->sgsiz * sizeof(int));
	      for (cc = 1; cc <= nv; cc++) {
		c = vec[cc];
		if (m == (int)regionid[c]) mask[c] = 1;
	      }
	      strcpy(iname, vname);
	      if (tin == NULL)
		t = NOTVALID;
	      else
		t =  *tin;
	      decode_indata(buf, buf1, iname, &t);
	      /* Check for variable substitution              */
	      /*
	      sprintf(buf1, "(%s%s)", key, buf);
	      sprintf(key1, "(%s=", vname);
	      if (find_token(buf1, key1, iname, ')')) {
		find_token(buf1, key, buf, '(');
	      }
	      */
	      if (strlen(i_rule)) {
		/*interp_data_us(master, buf1, iname, tr, mode, i_rule, mask, t);*/
		if (mode == 2)
		  value_init_sparse2d(master, tr, buf1, iname, i_rule, t, mask, NULL);
		if (mode == 3)
		  value_init_sparse3d(master, tr, buf1, iname, i_rule, t, mask, NULL);
	      } else
		interp_data_s(master, buf1, iname, tr, mode, mask, t);
	      continue;
	    }
	  }
	}
      }
      return(1);
    }
  }
  return(0);
}

/* END set_variable()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Interpolates bathymetry using bilinear methods via timeseries     */
/* libraries.                                                        */
/*-------------------------------------------------------------------*/
void interp_data_s(master_t *master,  /* Master data                 */
		   char *fname,       /* File name containing data   */
		   char *vname,       /* Variable name               */
		   double *ret,       /* Variable to interpolate     */
		   int mode,          /* 2D or 3D                    */
		   int *mask,         /* Mask for vec (optional)     */
		   double t           /* Time for temporal data      */
		   )
{
  char buf[MAXSTRLEN];
  int n, m, i, j, cc, c, cs;
  int nvec;          /* Size of vec                 */ 
  int *vec;          /* Indices to interpolate onto */
  int intype = -1;
  int fid;
  int ncerr;
  int idb;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  timeseries_t *ts = NULL;
  int bverbose = 0;
  double tin = (t == NOTVALID) ? 0.0 : t;
  poly_t *pl;

  if (mode == 3) {
    nvec = geom->b3_t;
    vec = geom->w3_t;
  }
  if (mode == 2) {
    nvec = geom->b2_t;
    vec = geom->w2_t;
  }
  /* Open the file and get dimensions                                */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    printf("Can't find input file %s\n", fname);
    hd_quit((char *)nc_strerror(ncerr));
  }
  if (ncw_var_id(fid, vname) < 0)
    hd_quit("Can't find variable %s in file %s\n", vname, fname);
  i = 0;
  while (in_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][1]), &nce1) == 0) {
      intype = i;
      break;
    }
    i++;
  }
  if (intype < 0)      
    hd_quit("interp_data_s: Can't find attributes in file %s for variable %s.\n", fname, vname);

  /* Initialize the timeseries file                                  */
  ts = (timeseries_t *)malloc(sizeof(timeseries_t));
  if (ts == NULL)
    hd_quit("interp_data_s: No memory available.\n");
  memset(ts, 0, sizeof(timeseries_t));
    
  /* Read the time series                                            */
  ts_read(fname, ts);
  if ((idb = ts_get_index(ts, fv_get_varname(fname, vname, buf))) == -1)
    hd_quit("interp_data_s: Can't find variable %s in file %s\n", vname, fname);

  /* Get a polgon of the perimeter of the bathy file                 */
  pl = nc2poly(fid, nce1, nce2, in_dims[intype][4], in_dims[intype][5], fname, ts);

  /* Read the values                                                 */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cs = geom->m2d[c];
    if (mask != NULL && !mask[c]) continue;
    /* Note: ts_eval_xy() for bathymetry files (idb = 'botz'     */
    /* or 'height') will return NaN if the xytoij() mapping does */
    /* not return valid indices (i.e. (x,y) lies outside the     */
    /* (i,j) bounds of the file).                                */
    /* This seems not entirely reliable though, so we use a      */
    /* bounding polygon for the perimeter.                       */
    /* Check if the location lies within the file's bounding     */
    /* perimeter.                                                */
    /*if (!poly_contains_point(pl, geom->cellx[cs], geom->celly[cs])) continue;*/
    if (mode == 2)
      ret[c] = ts_eval_xy(ts, idb, tin, geom->cellx[cs], geom->celly[cs]);
    if (mode == 3)
      ret[c] = ts_eval_xyz(ts, idb, tin, geom->cellx[cs], geom->celly[cs], geom->cellz[c]);
  }
  ts_free((timeseries_t*)ts);
  free(ts);
  poly_destroy(pl);
}

/* END interp_data_s()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks if a set_variable() specification is set.                  */
/*-------------------------------------------------------------------*/
int is_set_variable(char *fnames)
{
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  int nf = parseline(fnames, fields, MAXNUMARGS);
  int n;

  for (n = 0; n < nf; n++) {
    if (find_token(fields[n], "data=", buf, ']'))
      return(1);
    if (find_token(fields[n], "region=", buf, ']'))
      return(1);
  }
  return(0);
}

/* END is_set_variable()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Decodes a string for tracer data interpolation                    */
/*-------------------------------------------------------------------*/
void decode_indata(char *in, char *fname, char *iname, double *t)
{
  char buf[MAXSTRLEN], key[MAXSTRLEN];
  char inbuf[MAXSTRLEN], infile[MAXSTRLEN];
  int found = 0;

  sprintf(infile, "[%s]", in);

  /* Get the data name                                               */
  sprintf(inbuf, "[%s]", in);
  sprintf(key, "[");
  find_token(inbuf, key, fname, ']');

  /* Check for variable substitution              */
  strcpy(inbuf, in);
  sprintf(key, "(%s=", iname);
  if (find_token(inbuf, key, buf, ')')) {
    strcpy(iname, buf);
    sprintf(key, "[");
    find_token(infile, key, fname, '(');
    found = 1;
  }

  /* Check for a time input                                          */
  strcpy(inbuf, in);
  sprintf(key, "(t=");
  if (find_token(inbuf, key, buf, ')')) {
    if (!found) {
      sprintf(key, "[");
      find_token(infile, key, fname, '(');
    }
    tm_scale_to_secs(buf, t);
  }
}

/* END decode_indata()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured 2D tracer values from a netCDF file, packs into  */
/* a vector and interpolates onto an unstructured mesh.              */
/*-------------------------------------------------------------------*/
int value_init_sparse2d(master_t *master, /* Master data             */
			double *ret,      /* 2D array                */
			char *fname,      /* File name               */
			char *vname,      /* Variable name           */
			char *in_rule,    /* Interpolation type      */
			double t,         /* Time for temporal data  */
			int *rask,        /* Mask for vec (optional) */
			timeseries_t *ts  /* Timeseries structure    */
			)
{
  geometry_t *geom = master->geom;
  char i_rule[MAXSTRLEN];
  char buf[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nvar;
  double *x, *y, *v, **botz, **cellx, **celly, *celx, *cely, **var;
  double vmean, mv;
  int fid;
  int ncerr;
  int lond, latd;
  int varid;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  int ti;
  int n, i, j, c, cc;
  int bverbose = 0;
  int intype = -1;    /* 0 = EMS, > 0 = OFAM                         */
  int intime = 0;
  int baf = 1;

  /*-----------------------------------------------------------------*/
  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (strlen(in_rule) == 0) return(1);
  if (in_rule == NULL) {
    strcpy(i_rule, "nn_non_sibson");
  } else
    strcpy(i_rule, in_rule);

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("value_init_sparse_2d: Can't find input file %s\n", fname);
    return(1);
  }

  /* Get dimensions                                                  */
  i = 0;
  while (in_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][1]), &nce1) == 0) {
      intype = i;
      if (strcmp(in_dims[i][7], "time") == 0) intime = 1;
      break;
    }
    i++;
  }
  if (intype < 0) {
    hd_warn("value_init_sparse2d: Can't find %s attributes in file %s. Using gridded interpolation.\n", vname, fname);
    return(1);
  }
  if (ncw_var_id(fid, in_dims[intype][0]) < 0)
    baf = 0;

  /* Get the time index                                              */
  if (t != NOTVALID) {
    double rfrac;
    if (ts != NULL) 
      df_find_record(ts->df, t, &ti, &i, &rfrac);
    else
      ti = dump_choose_by_time_m(master, fid, t);
    if (ti == -1) return(1);
  }

  varid = ncw_var_id(fid, in_dims[intype][4]);
  nc_inq_varndims(fid, varid, &lond);
  varid = ncw_var_id(fid, in_dims[intype][5]);
  nc_inq_varndims(fid, varid, &latd);

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
  var = d_alloc_2d(nce1, nce2);

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
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, celx);
  } else if (lond == 2) {
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, cellx[0]);
  }
  if (latd == 1) {
    count[0] = nce2;
    count[1] = 0;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, cely);
  } else if (latd == 2) {
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, celly[0]);
  }
  count[0] = nce2;
  count[1] = nce1;
  if (baf)
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][0]), start, count, botz[0]);
  if (t != NOTVALID) {
    start[0] = ti;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = 1L;
    count[1] = nce2;
    count[2] = nce1;
    count[3] = 0;
 
 }
  nc_get_vara_double(fid, ncw_var_id(fid, vname), start, count, var[0]);

  /*-----------------------------------------------------------------*/
  /* Set the wet tracer vector (to interpolate from)                 */
  nvar = n = 0;
  c = nc_get_att_double(fid, ncw_var_id(fid, vname), "missing_value", &mv);
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (c >= 0 && var[j][i] == mv) continue;
      if (isnan(var[j][i])) continue;
      if (baf && isnan(botz[j][i])) continue;
      if (lond == 1 && isnan(celx[i])) continue;
      if (lond == 2 && isnan(cellx[j][i])) continue;
      if (latd == 1 && isnan(cely[j])) continue;
      if (latd == 2 && isnan(celly[j][i])) continue;
      if (baf && (botz[j][i] == LANDCELL || fabs(botz[j][i]) == fabs(NOTVALID))) continue;
      nvar++;
    }
  }

  if (nvar) {
    x = d_alloc_1d(nvar);
    y = d_alloc_1d(nvar);
    v = d_alloc_1d(nvar);
  } else
    hd_quit("value_init_sparse2d: Can't find valid %s values in file %s.\n", vname, fname);

  vmean = 0.0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (c >= 0 && var[j][i] == mv) continue;
      if (isnan(var[j][i])) continue;
      if (baf && isnan(botz[j][i])) continue;
      if (lond == 1 && isnan(celx[i])) continue;
      if (lond == 2 && isnan(cellx[j][i])) continue;
      if (latd == 1 && isnan(cely[j])) continue;
      if (latd == 2 && isnan(celly[j][i])) continue;
      if (baf && (botz[j][i] == LANDCELL || fabs(botz[j][i]) == fabs(NOTVALID))) continue;
      if (lond == 1)
	x[n] = celx[i];
      else if (lond == 2)
	x[n] = cellx[j][i];
      if (latd == 1)
	y[n] = cely[j]; 
      else if (latd == 2)
	y[n] = celly[j][i]; 
      v[n] = var[j][i];
      vmean += v[n];
      n++;
    }
  }

  if (n) vmean /= (double)n;
  if (lond == 1)
    d_free_1d(celx);
  else if (lond == 2)
    d_free_2d(cellx);
  if (latd == 1)
    d_free_1d(cely);
  else if (latd == 2)
    d_free_2d(celly);
  d_free_2d(botz);
  d_free_2d(var);
  nc_close(fid);

  /*-----------------------------------------------------------------*/
  /* Interpolate the tracer value                                    */
  gs = grid_interp_init(x, y, v, nvar, i_rule);
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (rask != NULL && !rask[c]) continue;
    ret[c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
    /*printf("%d %f : %f %f\n",c, ret[c], geom->cellx[c], geom->celly[c]);*/
    if (isnan(ret[c])) ret[c] = vmean;
    if (bverbose) printf("%d %f : %f %f\n",c, ret[c], geom->cellx[c], geom->celly[c]);
  }

  grid_specs_destroy(gs);
  hd_warn("2D Variable %s interpolated from file %s using %s.\n", vname, fname, i_rule);
  d_free_1d(x);
  d_free_1d(y);
  d_free_1d(v);
  return(0);
}

/* END value_init_sparse2d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured 3D tracer values from a netCDF file, packs into  */
/* a vector and interpolates onto an unstructured mesh.              */
/*-------------------------------------------------------------------*/
int value_init_sparse3d(master_t *master, /* Master data             */
			double *ret,      /* 3D array                */
			char *fname,      /* File name               */
			char *vname,      /* Variable name           */
			char *in_rule,    /* Interpolation type      */
			double t,         /* Time for temporal data  */
			int *rask,        /* Mask for vec (optional) */
			timeseries_t *ts  /* Timeseries structure    */
			)
{
  geometry_t *geom = master->geom;
  char i_rule[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nvar;
  double *x, *y, *v, **botz, **cellx, **celly, *celx, *cely;
  double ***var, **var2, *layers;
  int **mask;
  double *vmean, mv;
  int fid;
  int ncerr;
  int lond, latd;
  int varid;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  size_t nz;
  int ti;
  int n, i, j, k, c, cs, cc;
  int bverbose = 0;
  int intype = -1;    /* 0 = EMS, > 0 = OFAM                         */
  int intime = 0;
  int baf = 1;

  /*-----------------------------------------------------------------*/
  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (strlen(in_rule) == 0) return(1);
  if (in_rule == NULL) {
    strcpy(i_rule, "nn_non_sibson");
  } else
    strcpy(i_rule, in_rule);

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("value_init_sparse3d: Can't find input file %s\n", fname);
    return(1);
  }

  /* Get dimensions                                                  */
  i = 0;
  while (in_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][1]), &nce1) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][3]), &nz) == 0
	) {
      intype = i;
      if (strcmp(in_dims[i][7], "time") == 0) intime = 1;
      break;
    }
    i++;
  }
  if (intype < 0) {
    hd_warn("value_init_sparse3d: Can't find %s attributes in file %s. Using gridded interpolation.\n", vname, fname);
    return(1);
  }
  if (ncw_var_id(fid, in_dims[intype][0]) < 0)
    baf = 0;

  /* Get the time index                                              */
  if (t != NOTVALID) {
    double rfrac;
    if (intime == 0) {
      if (ts != NULL) 
	df_find_record(ts->df, t, &ti, &i, &rfrac);
      else
	ti = dump_choose_by_time_m(master, fid, t);
    } else {
      if (ts != NULL) 
	df_find_record(ts->df, t, &ti, &i, &rfrac);
      else
	ti = dump_choose_by_time_mom(master, fid, t);
    }
    if (ti == -1) return(1);
  }

  varid = ncw_var_id(fid, in_dims[intype][4]);
  nc_inq_varndims(fid, varid, &lond);
  varid = ncw_var_id(fid, in_dims[intype][5]);
  nc_inq_varndims(fid, varid, &latd);

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
  layers = d_alloc_1d(nz);
  botz = d_alloc_2d(nce1, nce2);
  mask = i_alloc_2d(nce1, nce2);
  var = d_alloc_3d(nce1, nce2, nz);
  start[0] = 0L;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][6]), start, count, layers);
  count[0] = nce2;
  count[1] = nce1;
  if (lond == 1) {
    count[0] = nce1;
    count[1] = 0;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, celx);
  } else if (lond == 2) {
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, cellx[0]);
  }
  if (latd == 1) {
    count[0] = nce2;
    count[1] = 0;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, cely);
  } else if (latd == 2) {
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, celly[0]);
  }
  count[0] = nce2;
  count[1] = nce1;
  if (baf)
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][0]), start, count, botz[0]);
  if (t != NOTVALID) {
    start[0] = ti;
    start[1] = 0L;
    start[2] = 0L;
    start[3] = 0L;
    count[0] = 1L;
    count[1] = nz;
    count[2] = nce2;
    count[3] = nce1;
  }
  nc_get_vara_double(fid, ncw_var_id(fid, vname), start, count, var[0][0]);

  /*-----------------------------------------------------------------*/
  /* Set the wet tracer vector (to interpolate from)                 */
  nvar = n = 0;
  c = nc_get_att_double(fid, ncw_var_id(fid, vname), "missing_value", &mv);
  for (k = 0; k < nz; k++) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (c >= 0 && var[nz-1][j][i] == mv) continue;
	if (isnan(var[nz-1][j][i])) continue;
	if (baf && isnan(botz[j][i])) continue;
	if (lond == 1 && isnan(celx[i])) continue;
	if (lond == 2 && isnan(cellx[j][i])) continue;
	if (latd == 1 && isnan(cely[j])) continue;
	if (latd == 2 && isnan(celly[j][i])) continue;
	if (baf && (botz[j][i] == LANDCELL || fabs(botz[j][i]) == fabs(NOTVALID))) continue;
	nvar++;
      }
    }
  }

  if (nvar) {
    x = d_alloc_1d(nvar);
    y = d_alloc_1d(nvar);
    v = d_alloc_1d(nvar);
  } else
    hd_quit("value_init_sparse3d: Can't find valid %s values in file %s.\n", vname, fname);

  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      mask[j][i] = -1;
      if (var[nz-1][j][i] == mv) continue;
      if (isnan(var[nz-1][j][i])) continue;
      if (baf && isnan(botz[j][i])) continue;
      if (lond == 1 && isnan(celx[i])) continue;
      if (lond == 2 && isnan(cellx[j][i])) continue;
      if (latd == 1 && isnan(cely[j])) continue;
      if (latd == 2 && isnan(celly[j][i])) continue;
      if (baf && (botz[j][i] == LANDCELL || fabs(botz[j][i]) == fabs(NOTVALID))) continue;
      if (lond == 1)
	x[n] = celx[i];
      else if (lond == 2)
	x[n] = cellx[j][i];
      if (latd == 1)
	y[n] = cely[j]; 
      else if (latd == 2)
	y[n] = celly[j][i]; 
      n++;
      mask[j][i] = 0;
      if (baf) {
	for (k = nz-1; k >= 0; k--) {
	  if (layers[k] > botz[j][i])
	    mask[j][i] = k;
	}
      }
    }
  }
  if (lond == 1)
    d_free_1d(celx);
  else if (lond == 2)
    d_free_2d(cellx);
  if (latd == 1)
    d_free_1d(cely);
  else if (latd == 2)
    d_free_2d(celly);
  d_free_2d(botz);
  nc_close(fid);

  /*-----------------------------------------------------------------*/
  /* Interpolate the tracer value                                    */
  var2 = d_alloc_2d(geom->sgsizS, nz);
  vmean = d_alloc_1d(nz);
  for (k = 0; k < nz; k++) {
    n = 0;
    vmean[k] = 0.0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (mask[j][i] >= 0) {
	  v[n] = var[k][j][i];
	  /* Set a no-gradient below the bottom                     */
	  if (k < mask[j][i]) v[n] = var[mask[j][i]][j][i];
	  /* Get the mean                                           */
	  vmean[k] += v[n];
	  n++;
	}
      }
    }
    /* Interpolate horizontally                                    */
    if (n) vmean[k] /= (double)n;
    gs = grid_interp_init(x, y, v, nvar, i_rule);
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      var2[k][c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
      if (isnan(var2[k][c])) var2[k][c] = vmean[k];
      if (bverbose) printf("%d %d %f : %f %f : %f\n",c, k, var2[k][c], geom->cellx[c], geom->celly[c], vmean[k]);
    }
    grid_specs_destroy(gs);
  }
  /* Interpolate vertically                                          */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    if (rask != NULL && !rask[cs]) continue;
    if (geom->cellz[c] >= layers[nz-1]) 
      ret[c] = var2[nz-1][cs];
    else if (geom->cellz[c] <= layers[0]) 
      ret[c] = var2[0][cs];
    else {
      k = nz-1;
      while(layers[k] > geom->cellz[c])
	k--;
      ret[c] = (var2[k][cs] - var2[k+1][cs]) * (geom->cellz[c] - layers[k]) /
	(layers[k] - layers[k+1]) + var2[k][cs];
    }
  }

  d_free_3d(var);
  d_free_2d(var2);
  d_free_1d(v);
  d_free_1d(vmean);
  i_free_2d(mask);
  nc_close(fid);
  hd_warn("3D Variable %s interpolated from file %s using %s.\n", vname, fname, i_rule);
  return(0);
}

/* END value_init_sparse3d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured 2D tracer values from a netCDF file, packs into  */
/* a vector and interpolates onto an unstructured mesh.              */
/*-------------------------------------------------------------------*/
int value_init_sparse2d_o(master_t *master,   /* Master data           */
			double *ret,        /* 2D array              */
			char *fname,        /* File name             */
			char *vname,        /* Variable name         */
			char *in_rule       /* Interpolation type    */
			)
{
  geometry_t *geom = master->geom;
  char i_rule[MAXSTRLEN];
  char buf[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nvar;
  double *x, *y, *v, **botz, **cellx, **celly, **var;
  double vmean, mv;
  int fid;
  int ncerr;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  int ti;
  int n, i, j, c, cc;
  int bverbose = 0;
  int intype = -1;    /* 0 = EMS, > 0 = OFAM                         */
  int intime = 0;

  /*-----------------------------------------------------------------*/
  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (strlen(in_rule) == 0) return(1);
  if (in_rule == NULL) {
    strcpy(i_rule, "nn_non_sibson");
  } else
    strcpy(i_rule, in_rule);

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("value_init_sparse_2d: Can't find input file %s\n", fname);
    return(1);
  }

  /* Get dimensions                                                  */
  i = 0;
  while (in_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][1]), &nce1) == 0) {
      intype = i;
      if (strcmp(in_dims[i][7], "time") == 0) intime = 1;
      break;
    }
    i++;
  }
  if (intype < 0) {
    hd_warn("value_init_sparse2d: Can't find %s attributes in file %s. Using gridded interpolation.\n", vname, fname);
    return(1);
  }

  /* Get the time index                                              */
  ti = dump_choose_by_time_m(master, fid, schedule->start_time);
  if (ti == -1) return(1);

  /* Allocate and read                                               */
  cellx = d_alloc_2d(nce1, nce2);
  celly = d_alloc_2d(nce1, nce2);
  botz = d_alloc_2d(nce1, nce2);
  var = d_alloc_2d(nce1, nce2);
  start[0] = 0L;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = nce2;
  count[1] = nce1;
  count[2] = 0;
  count[3] = 0;
  nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, cellx[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, celly[0]);
  nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][0]), start, count, botz[0]);
  start[0] = ti;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nce2;
  count[2] = nce1;
  count[3] = 0;
  nc_get_vara_double(fid, ncw_var_id(fid, vname), start, count, var[0]);

  /*-----------------------------------------------------------------*/
  /* Set the wet tracer vector (to interpolate from)                 */
  nvar = n = 0;
  c = nc_get_att_double(fid, ncw_var_id(fid, vname), "missing_value", &mv);
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (c >= 0 && var[j][i] == mv) continue;
      if (isnan(var[j][i])) continue;
      if (isnan(botz[j][i])) continue;
      if (isnan(cellx[j][i])) continue;
      if (isnan(celly[j][i])) continue;
      if (botz[j][i] != LANDCELL && botz[j][i] != NOTVALID) nvar++;
    }
  }

  if (nvar) {
    x = d_alloc_1d(nvar);
    y = d_alloc_1d(nvar);
    v = d_alloc_1d(nvar);
  } else
    hd_quit("value_init_sparse2d: Can't find valid %s values in file %s.\n", vname, fname);

  vmean = 0.0;
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      if (var[j][i] == mv) continue;
      if (isnan(botz[j][i])) continue;
      if (isnan(cellx[j][i])) continue;
      if (isnan(celly[j][i])) continue;
      if (botz[j][i] != LANDCELL && botz[j][i] != NOTVALID) {
	x[n] = cellx[j][i];
	y[n] = celly[j][i]; 
	v[n] = var[j][i];
	vmean += v[n];
	n++;
      }
    }
  }

  if (n) vmean /= (double)n;
  d_free_2d(cellx);
  d_free_2d(celly);
  d_free_2d(botz);
  d_free_2d(var);
  nc_close(fid);

  /*-----------------------------------------------------------------*/
  /* Interpolate the tracer value                                    */
  gs = grid_interp_init(x, y, v, nvar, i_rule);
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    ret[c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
    if (isnan(ret[c])) ret[c] = vmean;
    if (bverbose) printf("%d %f : %f %f\n",c, ret[c], geom->cellx[c], geom->celly[c]);
  }
  grid_specs_destroy(gs);
  hd_warn("2D Variable %s interpolated from file %s using %s.\n", vname, fname, i_rule);
  return(0);
}

/* END value_init_sparse2d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Reads structured 3D tracer values from a netCDF file, packs into  */
/* a vector and interpolates onto an unstructured mesh.              */
/*-------------------------------------------------------------------*/
int value_init_sparse3d_o(master_t *master,   /* Master data           */
			double *ret,        /* 3D array              */
			char *fname,        /* File name             */
			char *vname,        /* Variable name         */
			char *in_rule       /* Interpolation type    */
			)
{
  geometry_t *geom = master->geom;
  char i_rule[MAXSTRLEN];
  GRID_SPECS *gs = NULL;
  int nvar;
  double *x, *y, *v, **botz, **cellx, **celly;
  double ***var, **var2, *layers;
  int **mask;
  double *vmean, mv;
  int fid;
  int ncerr;
  size_t start[4];
  size_t count[4];
  size_t nce1;
  size_t nce2;
  size_t nz;
  int ti;
  int n, i, j, k, c, cs, cc;
  int bverbose = 0;
  int intype = -1;    /* 0 = EMS, > 0 = OFAM                         */
  int intime = 0;

  /*-----------------------------------------------------------------*/
  /* Interpolation method. Options:                                  */
  /* linear, cubic, nn_sibson, nn_non_sibson, average                */
  if (strlen(in_rule) == 0) return(1);
  if (in_rule == NULL) {
    strcpy(i_rule, "nn_non_sibson");
  } else
    strcpy(i_rule, in_rule);

  /* Open the dump file for reading                                  */
  if ((ncerr = nc_open(fname, NC_NOWRITE, &fid)) != NC_NOERR) {
    hd_warn("value_init_sparse3d: Can't find input file %s\n", fname);
    return(1);
  }

  /* Get dimensions                                                  */
  i = 0;
  while (in_dims[i][0] != NULL) {
    if (nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][2]), &nce2) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][1]), &nce1) == 0 &&
	nc_inq_dimlen(fid, ncw_dim_id(fid, in_dims[i][3]), &nz) == 0
	) {
      intype = i;
      if (strcmp(in_dims[i][7], "time") == 0) intime = 1;
      break;
    }
    i++;
  }

  if (intype < 0) {
    hd_warn("value_init_sparse3d: Can't find %s attributes in file %s. Using gridded interpolation.\n", vname, fname);
    return(1);
  }

  /* Get the time index                                              */
  if (intime == 0)
    ti = dump_choose_by_time_m(master, fid, schedule->start_time);
  else
    ti = dump_choose_by_time_mom(master, fid, schedule->start_time);
  if (ti == -1) return(1);

  /* Allocate and read                                               */
  cellx = d_alloc_2d(nce1, nce2);
  celly = d_alloc_2d(nce1, nce2);
  layers = d_alloc_1d(nz);
  botz = d_alloc_2d(nce1, nce2);
  mask = i_alloc_2d(nce1, nce2);
  var = d_alloc_3d(nce1, nce2, nz);
  start[0] = 0L;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = nz;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;
  if (intype == 0) {
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][6]), start, count, layers);
    count[0] = nce2;
    count[1] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, cellx[0]);
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, celly[0]);
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][0]), start, count, botz[0]);
  } else {
    double *lon = d_alloc_1d(nce1);
    double *lat = d_alloc_1d(nce2);
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][6]), start, count, layers);
    count[0] = nce1;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][4]), start, count, lon);
    count[1] = nce2;
    nc_get_vara_double(fid, ncw_var_id(fid, in_dims[intype][5]), start, count, lat);
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	cellx[j][i] = lon[i];
	celly[j][i] = lon[j];
	botz[j][i] = -5000.0;
      }
    }
    d_free_1d(lon);
    d_free_1d(lat);
  }
  start[0] = ti;
  start[1] = 0L;
  start[2] = 0L;
  start[3] = 0L;
  count[0] = 1L;
  count[1] = nz;
  count[2] = nce2;
  count[3] = nce1;
  nc_get_vara_double(fid, ncw_var_id(fid, vname), start, count, var[0][0]);

  /*-----------------------------------------------------------------*/
  /* Set the wet tracer vector (to interpolate from)                 */
  nvar = n = 0;
  c = nc_get_att_double(fid, ncw_var_id(fid, vname), "missing_value", &mv);
  for (k = 0; k < nz; k++) {
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (c >= 0 && var[nz-1][j][i] == mv) continue;
	if (isnan(var[nz-1][j][i])) continue;
	if (isnan(botz[j][i])) continue;
	if (isnan(cellx[j][i])) continue;
	if (isnan(celly[j][i])) continue;
	if (k == 0 && botz[j][i] != LANDCELL && botz[j][i] != NOTVALID) nvar++;
      }
    }
  }

  x = d_alloc_1d(nvar);
  y = d_alloc_1d(nvar);
  v = d_alloc_1d(nvar);
  for (j = 0; j < nce2; j++) {
    for (i = 0; i < nce1; i++) {
      mask[j][i] = -1;
      if (var[nz-1][j][i] == mv) continue;
      if (isnan(var[nz-1][j][i])) continue;
      if (isnan(botz[j][i])) continue;
      if (isnan(cellx[j][i])) continue;
      if (isnan(celly[j][i])) continue;
      if (botz[j][i] != LANDCELL && botz[j][i] != NOTVALID) {
	x[n] = cellx[j][i];
	y[n] = celly[j][i]; 
	n++;
	for (k = nz-1; k >= 0; k--) {
	  if (layers[k] > botz[j][i])
	    mask[j][i] = k;
	}
      }
    }
  }
  d_free_2d(cellx);
  d_free_2d(celly);
  d_free_2d(botz);
  nc_close(fid);

  /*-----------------------------------------------------------------*/
  /* Interpolate the tracer value                                    */
  var2 = d_alloc_2d(geom->szcS, nz);
  vmean = d_alloc_1d(nz);
  for (k = 0; k < nz; k++) {
    n = 0;
    vmean[k] = 0.0;
    for (j = 0; j < nce2; j++) {
      for (i = 0; i < nce1; i++) {
	if (mask[j][i] >= 0) {
	  v[n] = var[k][j][i];
	  /* Set a no-gradient below the bottom                     */
	  if (k < mask[j][i]) v[n] = var[mask[j][i]][j][i];
	  /* Get the mean                                           */
	  vmean[k] += v[n];
	  n++;
	}
      }
    }
    /* Interpolate horizontally                                    */
    if (n) vmean[k] /= (double)n;
    gs = grid_interp_init(x, y, v, nvar, i_rule);
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      var2[k][c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
      if (isnan(var2[k][c])) var2[k][c] = vmean[k];
      if (bverbose) printf("%d %d %f : %f %f : %f\n",c, k, var2[k][c], geom->cellx[c], geom->celly[c], vmean[k]);
    }
    grid_specs_destroy(gs);
  }
  /* Interpolate vertically                                          */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    if (geom->cellz[c] >= layers[nz-1]) 
      ret[c] = var2[nz-1][cs];
    else if (geom->cellz[c] <= layers[0]) 
      ret[c] = var2[0][cs];
    else {
      k = nz-1;
      while(layers[k] > geom->cellz[c])
	k--;
      ret[c] = (var2[k][cs] - var2[k+1][cs]) * (geom->cellz[c] - layers[k]) /
	(layers[k] - layers[k+1]) + var2[k][cs];
    }
  }

  d_free_3d(var);
  d_free_2d(var2);
  d_free_1d(v);
  d_free_1d(vmean);
  i_free_2d(mask);
  nc_close(fid);
  hd_warn("2D Variable %s interpolated from file %s using %s.\n", vname, fname, i_rule);
  return(0);
}

/* END value_init_sparse3d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Fills cells in the grid that are non-wet cells in the data that   */
/* being interpolated from. Also performs filtering of the input     */
/* data.                                                             */
/*-------------------------------------------------------------------*/
void tracer_fill(master_t *master,   /* Master structure             */
		 timeseries_t *ts,   /* File to interpolate from     */
		 char *fname,        /* Name of interpolation file   */
		 double *ret,        /* Variable to fill             */
		 int sz,             /* Size of vec                  */
		 int *vec,           /* Cells to process             */
		 int vc,             /* Size of vec                  */
		 tracer_info_t *tr,  /* Tracer info                  */
		 double fill         /* Tracer fill value            */
		 )
{
  geometry_t *geom = master->sgrid;
  dump_data_t *dumpdata = master->dumpdata;
  char buf[MAXSTRLEN];
  int n, m, k, c, cm, cc, cs, id, *mask, *map, *done;
  double val, *tmp, **a2, ***a3;

  strcpy(buf, fname);

  /*
   * special handling for RECOM
   * By default, the RECOM initialisation files are non-cascade
   * searched simple_cf files which use 1e35 instead of NaN's
   */
  if (ts_is_recom(ts)) {
    for (cc = 1; cc <= vc; cc++) {
      c = vec[cc];
      if (ret[c] > 1e20) ret[c] = NaN;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Fill data outside the grid being interpolated from.             */
  if (master->trfilter & TRF_FILL) {
    /* Get the land mask from the bathymetry variable                */
      /* Good values (i.e. ones that don't contain a NaN or land in  */
      /* the interpolation) are set to 1.                            */
    mask = i_alloc_1d(geom->sgsizS);
    done = i_alloc_1d(sz);
    id = ts_get_index(ts, fv_get_varname(fname, "botz", buf));
    if (id < 0) {
      hd_warn("tracerfill: The file '%s' does not contain the tracer 'bot_z'\n");
      return;
    } else {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	val = ts_eval_xy(ts, id, master->t,
			 geom->cellx[c],
			 geom->celly[c]);
	/* Check depth sign */
	if (ts_var_z_is_depth(ts, id)) val = -val;
	mask[c] = (!isnan(val) && (val < 0.0 && val != NOTVALID)) ? 1 : 0;
      }
    }

    /* Set the wet cell mask in the tracer                           */
    memset(done, 0, sz * sizeof(int));
    for (cc = 1; cc <= vc; cc++) {
      c = vec[cc];
      cs = geom->m2d[c];
      if (mask[cs]) done[c] = 1;
    }

    /* Fill the tracer variable                                      */
    for (cc = 1; cc <= vc; cc++) {
      int cs;
      c = vec[cc];
      cs = geom->m2d[c];
      for (n = 1; n <= geom->npe[cs]; n++) {
	cm = c;
	cs = geom->m2d[cm];
	map = geom->c2c[n];
	if (mask[cs] && !mask[map[cs]] || (!isnan(ret[cm]) && isnan(ret[map[cm]]))) {
	  while (cm != map[cm] || (!isnan(ret[cm]) && isnan(ret[map[cm]]))) {
	    if (cm == map[cm]) break;
	    if (!mask[map[cs]] || isnan(ret[map[cm]])) {
	      if (!done[map[cm]])
		ret[map[cm]] = ret[cm];
	      else
		ret[map[cm]] = 0.5 * (ret[cm] + ret[map[cm]]);
	      done[map[cm]] = 1;
	    }
	    cm = map[cm];
	    cs = geom->m2d[cm];
	  }
	}
      }
    }
    i_free_1d(mask);
    i_free_1d(done);
  }

  /*-----------------------------------------------------------------*/
  /* Perform filtering if required                                   */
  /* Set a no-gradient beyond the open boundary                      */
  for (n = 0; n < geom->nobc; n++) {
    int e, ee;
    open_bdrys_t *open = geom->open[n];
    m = (sz == geom->sgsiz) ? open->no3_e1 : open->no2_e1;
    for (ee = 1; ee <= m; ee++) {
      c = cm = open->obc_e2[ee];
      e = open->obc_e1[ee];
      while (cm != open->omape[e][cm]) {
	ret[cm] = ret[c];
	cm = open->omape[e][cm];
      }
    }
  }
  n = (sz == geom->sgsiz) ? geom->nbpt : geom->nbptS;
  tmp = d_alloc_1d(sz);

  /* Smoothing filter                                                */
  if (master->trfilter & TRF_SMOO) {
    memcpy(tmp, ret, sz * sizeof(double));
    for (cc = 1; cc <= n; cc++) {
      c = geom->bpt[cc];
      cm = geom->bin[cc];
      ret[c] = ret[cm];
    }
    for (cc = 1; cc <= vc; cc++) {
      c = vec[cc];
      tmp[c] = cvol1(master, ret, c, 0);
    }
    memcpy(ret, tmp, sz * sizeof(double));
  }

  /* Shuman filter                                                   */
  if (master->trfilter & TRF_SHUM) {
    memcpy(tmp, ret, sz * sizeof(double));
    for (cc = 1; cc <= n; cc++) {
      c = geom->bpt[cc];
      cm = geom->bin[cc];
      ret[c] = ret[cm];
    }
    for (cc = 1; cc <= vc; cc++) {
      c = vec[cc];
      tmp[c] = shuman(geom, ret, c);
    }
    memcpy(ret, tmp, sz * sizeof(double));
  }

  /* Shapiro filter                                                  */
  if (master->trfilter & TRF_SHAP) {
    memcpy(tmp, ret, sz * sizeof(double));
    for (cc = 1; cc <= n; cc++) {
      c = geom->bpt[cc];
      cm = geom->bin[cc];
      tmp[c] = tmp[cm];
    }
    for (cc = 1; cc <= vc; cc++) {
      int j, cn, cs;
      double val = ret[c];
      c = vec[cc];
      cs = geom->m2d[c];
      tmp[c] = 0.0;
      for (j = 1; j <= geom->npe[cs]; j++) {
	shapiro_smooth(geom, ret, c, j);
	tmp[c] += ret[c];
	ret[c] = val;
      }
      tmp[c] /= (double)geom->npe[cs];
    }
    memcpy(ret, tmp, sz * sizeof(double));
  }

  /* Median filter                                                   */
  if (master->trfilter & TRF_MEDI) {
    memcpy(tmp, ret, sz * sizeof(double));
    for (cc = 1; cc <= n; cc++) {
      c = geom->bpt[cc];
      cm = geom->bin[cc];
      ret[c] = ret[cm];
    }
    for (cc = 1; cc <= vc; cc++) {
      c = vec[cc];
      tmp[c] = con_median(geom, ret, tr->valid_range_wc[1], c, ST_SQ3);
    }
    memcpy(ret, tmp, sz * sizeof(double));
  }
  d_free_1d(tmp);

  /*-----------------------------------------------------------------*/
  /* Check for NaNs and replace with nearest non-NaN                 */
  if (tr != NULL && strncasecmp(tr->tag, "DA_OBS", 6) == 0) return;
  a2 = NULL;
  a3 = NULL;
  if (sz == geom->sgsiz) {
    a3 = d_alloc_3d(geom->nce1, geom->nce2, geom->nz);
    s2c_3d(geom, ret, a3, geom->nce1, geom->nce2, geom->nz);
  } else {
    a2 = d_alloc_2d(geom->nce1, geom->nce2);
    s2c_2d(geom, ret, a2, geom->nce1, geom->nce2);
  }

  for (cc = 1; cc <= vc; cc++) {
    c = vec[cc];
    if (isnan(ret[c])) {
      int ci, cj;	
      double **a;
      if (sz == geom->sgsiz) {
	k = geom->s2k[c];
	a = a3[k];
      } else {
	k = geom->nz - 1;
	a = a2;
      }
      if (find_closest_nonnan(dumpdata, a, geom->s2i[c], geom->s2j[c], k, &ci, &cj))
	ret[c] = a[cj][ci];
      else {
	/* Search upwards for a valid value                          */
	ci = c;
	while (isnan(ret[ci]) && ci != geom->zp1[ci])
	  ci = geom->zp1[ci];
	ret[c] = (isnan(ret[ci])) ? fill : ret[ci];
      }
    }
  }

  if (a2) d_free_2d(a2);
  if (a3) d_free_3d(a3);
}

/* END tracer_fill()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_profile normalizes the density profile to top and densities  */
/* and creates a new profile by fitting top and bottom tracer        */
/* endpoints.                                                        */
/*-------------------------------------------------------------------*/
void dens_profile(master_t *master, geometry_t *geom, double *ret, double v1, double v2)
{
  double *d1 = NULL, dt, db;
  int cc, c, cs, cb, k;

  density_m(master);
  d1 = d_alloc_1d(geom->nz + 1);
  /* Get the maximum bottom density                              */
  db = 0.0;
  for (cc = 1; cc <= geom->b2_t; cc++) {
    db = max(db, master->dens_0[geom->bot_t[cc]]);
  }
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = cs = geom->w2_t[cc];
    /* Get surface densities                                     */
    dt = master->dens_0[cs];
    /* Get the normalized density profile                        */
    while (c != geom->zm1[c]) {
      k = geom->s2k[c];
      d1[k] = (db == dt) ? 1.0 : (master->dens_0[c] - dt) / (db - dt);
      c = geom->zm1[c];
    }
    for(k--; k >= 0; k--)
      d1[k] = d1[k+1];
    /* Reconstruct the profile.                                  */
    c = cs;
    while (c != geom->zm1[c]) {
      k = geom->s2k[c];
      ret[c] = v1 + d1[geom->s2k[c]] * (v2 - v1);
      c = geom->zm1[c];
    }
  }
  d_free_1d(d1);
}

/* END dens_profile()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_grad normalizes the density profile to the sum of density    */
/* gradients and creates a new profile by normalizing tracer to the  */
/* sum of gradients to a particular depth.                           */
/*-------------------------------------------------------------------*/
void dens_grad(master_t *master, geometry_t *geom, double *ret, double v1, double v2) 
{
  double *d1 = NULL, dt, db;
  int cc, c, cs, cb;
  int zm1, zp1;
  double ddendz, dtdz, vb;
  density_m(master);
  /* Set a no-gradient at the bottom                            */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->bot_t[cc];
    master->dens[geom->zm1[c]] = master->dens[c];
  }
  /* Find the maximum depth                                     */
  dt = 0.0;
  d1 = d_alloc_1d(geom->nz + 1);
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (geom->botz[c] < dt) {
      dt = geom->botz[c];
      cb = geom->bot_t[cc];
      cs = c;
    }
  }
  /* Get the profile at the deepest point                        */
  ddendz = 0.0;
  c = cs;
  zm1 = geom->zm1[c];
  while (c != zm1 ) {
    ddendz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
    c = zm1;
    zm1 = geom->zm1[c];
  }
  c = cb = geom->zp1[c];
  zp1 = geom->zp1[c];
  dtdz = 0.0;
  while (c != geom->zp1[c]) {
    d1[geom->s2k[c]] = v2 + (v1 - v2) * dtdz / ddendz;
    zm1 = c;
    c = zp1;
    zp1 = geom->zp1[zp1];
    dtdz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
  }
  /* Get the profile everywhere else                               */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = cs = geom->w2_t[cc];
    /* Get the sum of density gradients                          */
    ddendz = 0.0;
    zm1 = geom->zm1[c];
    while (c != zm1 ) {
      ddendz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
      c = zm1;
      zm1 = geom->zm1[c];
    }
    c = cb = geom->zp1[c];
    zp1 = geom->zp1[c];
    dtdz = 0.0;
    vb = d1[geom->s2k[cb]];
    while (c != geom->zp1[c]) {
      ret[c] = vb + (v1 - vb) * dtdz / ddendz;
      zm1 = c;
      c = zp1;
      zp1 = geom->zp1[zp1];
      dtdz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
    }
  }
  d_free_1d(d1);
}

/* END dens_grad()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_scale fits a value to the surface and bottom and a profile   */
/* by incrementing this value by a scaled vertical density           */
/* difference.                                                       */
/*-------------------------------------------------------------------*/
void dens_scale(master_t *master, geometry_t *geom, double *ret, double v1, double v2, char *sgn) 
{
  double *d1 = NULL, dt, db, vs;
  int cc, c, cs, cb;
  int zm1, zp1;
  double s = 1.0;

  if (strlen(sgn))
    if ((strcmp(sgn, "n") == 0) || (strcmp(sgn, "inverse") == 0))
      s = -1.0;
  density_m(master);
  if (v2 >= 0.0) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = cs = geom->w2_t[cc];
      /* Get the normalized density profile                      */
      vs = v1;
      ret[c] = vs;
      zp1 = c;
      c = geom->zm1[c];
      while (c != geom->zm1[c]) {
	dt = master->dens_0[zp1];
	db = master->dens_0[c];
	ret[c] = (db == 0.0) ? vs : vs + s * vs * v2 * 2.0 * (db - dt) / (dt + db);
	vs = ret[c];
	zp1 = c;
	c = geom->zm1[c];
      }
    }
  } else {
    /* Find the maximum depth                                    */
    dt = 0.0;
    d1 = d_alloc_1d(geom->nz + 1);
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      if (geom->botz[c] < dt) {
	dt = geom->botz[c];
	cb = geom->bot_t[cc];
      }
    }
    /* Get the profile at the deepest point                      */
    c = cb;
    vs = v1;
    d1[geom->s2k[c]] = vs;
    zm1 = c;
    c = geom->zp1[c];
    while (c != geom->zp1[c]) {
      dt = master->dens_0[c];
      db = master->dens_0[zm1];
      d1[geom->s2k[c]] = (db == 0.0) ? vs : vs + s * vs * v2 * 2.0 * (db - dt) / (dt + db);
      vs = d1[geom->s2k[c]];
      zm1 = c;
      c = geom->zp1[c];
    }
    /* Get the normalized density profile                        */
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->bot_t[cc];
      ret[c] = vs = d1[geom->s2k[c]];
      zm1 = c;
      c = geom->zp1[c];
      while (c != geom->zp1[c]) {
	dt = master->dens_0[c];
	db = master->dens_0[zm1];
	ret[c] = vs + s * vs * v2 * 2.0 * (db - dt) / (dt + db);
	vs = ret[c];
	zm1 = c;
	c = geom->zp1[c];
      }
    }
    d_free_1d(d1);
  }
}

/* END dens_scale()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_scale fits a value to the surface and bottom and a profile   */
/* by incrementing this value by a scaled vertical density           */
/* difference.                                                       */
/*-------------------------------------------------------------------*/
void dens_scale_file(master_t *master, 
		     geometry_t *geom, 
		     double *ret, 
		     char *vname, 
		     char *infile, 
		     double v2, 
		     char *sgn) 
{
  double *d1 = NULL, dt, db, vs;
  int cc, c, cs, cb;
  int zm1, zp1;
  double s = 1.0;
  timeseries_t *ts;
  int id; 
  char buf[MAXSTRLEN];

  if (strlen(sgn))
    if ((strcmp(sgn, "n") == 0) || (strcmp(sgn, "inverse") == 0))
      s = -1.0;
  ts = hd_ts_read(master, infile, 0);
  if (ts->t != NULL) {
    double *indat = d_alloc_1d(geom->sgsizS);
    id = ts_get_index(ts, fv_get_varname(infile, vname, buf));
    if (id < 0)
      hd_warn("val_init_3d: The file '%s' does not contain the tracer '%s'\n", buf, vname);
    else {
      /* Interpolate                                             */
      for (cc = 1; cc <= geom->b2_t; cc++) {
	if (v2 >= 0.0)
	  c = geom->w2_t[cc];
	else
	  c = geom->bot_t[cc];
	cs = geom->m2d[c];
	indat[cc] = ts_eval_xyz(ts, id, master->t,
				geom->cellx[cs],
				geom->celly[cs],
				geom->cellz[c]);
      }
      hd_ts_free(master, ts);
    }
    density_m(master);
    if (v2 >= 0.0) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = cs = geom->w2_t[cc];
	/* Get the normalized density profile                      */
	ret[c] = vs = indat[cc];
	zp1 = c;
	c = geom->zm1[c];
	while (c != geom->zm1[c]) {
	  dt = master->dens_0[zp1];
	  db = master->dens_0[c];
	  ret[c] = (db == 0.0) ? vs : vs + s * vs * v2 * 2.0 * (db - dt) / (dt + db);
	  vs = ret[c];
	  zp1 = c;
	  c = geom->zm1[c];
	}
      }
    } else {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->bot_t[cc];
	/* Get the normalized density profile                        */
	ret[c] = vs = indat[cc];
	zm1 = c;
	c = geom->zp1[c];
	while (c != geom->zp1[c]) {
	  dt = master->dens_0[c];
	  db = master->dens_0[zm1];
	  ret[c] = vs + s * vs * v2 * 2.0 * (db - dt) / (dt + db);
	  vs = ret[c];
	  zm1 = c;
	  c = geom->zp1[c];
	}
      }
    }
    d_free_1d(indat);
  }
}

/* ENS dens_scale_file()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_grad normalizes the density profile to the sum of density    */
/* gradients and creates a new profile by normalizing tracer to the  */
/* sum of gradients to a particular depth.                           */
/*-------------------------------------------------------------------*/
void dens_grad_file(master_t *master, 
		    geometry_t *geom, 
		    double *ret, 
		    char *vname, 
		    char *infile) 
{
  double *d1 = NULL, dt, db, vs;
  int cc, c, cs, cb;
  int zm1, zp1;
  timeseries_t *ts;
  int id; 
  double ddendz, dtdz;
  char buf[MAXSTRLEN];

  /* Read in the data                                            */
  ts = hd_ts_read(master, infile, 0);
  if (ts->t != NULL) {
    double *indat = d_alloc_1d(geom->sgsiz);
    id = ts_get_index(ts, fv_get_varname(infile, vname, buf));
    if (id < 0)
      hd_warn("val_init_3d: The file '%s' does not contain the tracer '%s'\n", buf, vname);
    else {
      /* Interpolate                                             */
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	indat[c] = ts_eval_xyz(ts, id, master->t,
			       geom->cellx[cs],
			       geom->celly[cs],
			       geom->cellz[c]);
      }
      hd_ts_free(master, ts);
    }
    density_m(master);
    /* Set a no-gradient at the bottom                            */
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->bot_t[cc];
      master->dens[geom->zm1[c]] = master->dens[c];
    }
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = cs = geom->w2_t[cc];
      /* Get the sum of density gradients                          */
      ddendz = 0.0;
      zm1 = geom->zm1[c];
      while (c != zm1 ) {
	ddendz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
	c = zm1;
	zm1 = geom->zm1[c];
      }
      c = cb = geom->zp1[c];
      zp1 = geom->zp1[c];
      dtdz = 0.0;
      while (c != geom->zp1[c]) {
	ret[c] = indat[cb] + (indat[cs] - indat[cb]) * dtdz / ddendz;
	zm1 = c;
	c = zp1;
	zp1 = geom->zp1[zp1];
	dtdz += (master->dens[c] - master->dens[zm1]) / (geom->cellz[c] - geom->cellz[zm1]);
      }
    }
    d_free_1d(indat);
  }
}

/* END dens_grad_file()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* dens_scale fits a value to the surface and bottom and a profile   */
/* by incrementing this value by a scaled vertical density           */
/* difference. Inverse profile can be specified using sgn='n' and    */
/* a copy below depth v1 can be specified using sgn='b'.             */ 
/*-------------------------------------------------------------------*/
void dens_scale_mid(master_t *master, 
		    geometry_t *geom, 
		    double *ret,
		    char *vname, 
		    char *infile, 
		    double v1, 
		    double v2, 
		    char *sgn) 
{
  double *d1 = NULL, dt, db, vs;
  int cc, c, cs, cb, cm;
  int zm1, zp1;
  timeseries_t *ts;
  int id; 
  char buf[MAXSTRLEN];
  int ff = 0;

  if (strlen(sgn)) {
    if ((strcmp(sgn, "c") == 0))
      ff = 1;;
    if ((strcmp(sgn, "t") == 0))
      ff = 2;
  }
  ts = hd_ts_read(master, infile, 0);
  if (ts->t != NULL) {
    double *indat = d_alloc_1d(geom->sgsiz);
    id = ts_get_index(ts, fv_get_varname(infile, vname, buf));
    if (id < 0)
      hd_warn("val_init_3d: The file '%s' does not contain the tracer '%s'\n", buf, vname);
    else {
      /* Interpolate                                             */
      for (cc = 1; cc <= geom->b3_t; cc++) {
	c = geom->w3_t[cc];
	cs = geom->m2d[c];
	indat[c] = ts_eval_xyz(ts, id, master->t,
				geom->cellx[cs],
				geom->celly[cs],
				geom->cellz[c]);
      }
      hd_ts_free(master, ts);
    }
    density_m(master);
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = cs = geom->w2_t[cc];
      /* Find the layer corresponding to v1                      */
      cb = geom->bot_t[cc];
      while (c != geom->zm1[c] && fabs(geom->cellz[c]) < fabs(v1)) {
	c = geom->zm1[c];
      }
      c = cm = geom->zp1[c];
      /* Get the normalized density profile upwards              */
      ret[c] = vs = indat[c];
      zm1 = c;
      c = geom->zp1[c];
      while (c != geom->zp1[c]) {
	dt = master->dens_0[c];
	db = master->dens_0[zm1];
	ret[c] = vs - vs * v2 * 2.0 * (db - dt) / (dt + db);
	vs = ret[c];
	zm1 = c;
	c = geom->zp1[c];
      }
      dt = master->dens_0[c];
      db = master->dens_0[zm1];
      ret[c] = vs - vs * v2 * 2.0 * (db - dt) / (dt + db);
      /* Get the normalized density profile downwards            */
      zp1 = cm;
      c = geom->zm1[cm];
      vs = indat[cm];
      while (c != geom->zm1[c]) {
	if (ff == 1)
	  ret[c] = indat[c];
	else {
	  dt = master->dens_0[zp1];
	  db = master->dens_0[c];
	  ret[c] = (db == 0.0) ? vs : vs + vs * v2 * 2.0 * (db - dt) / (dt + db);
	  if (ff == 2 && ret[c] > indat[c]) ret[c] = indat[c];
	  vs = ret[c];
	}
	zp1 = c;
	c = geom->zm1[c];
      }
    }
    d_free_1d(indat);
  }
}

/* ENS dens_scale_mid()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise a 2d array from a number, array or file     */
/* input.                                                            */
/*-------------------------------------------------------------------*/
void value_init_2d(master_t *master,     /* Master data              */
		   double *ret,          /* 2D tracer to fill        */
		   FILE *fp,             /* File pointer             */
		   char *fname,          /* File name                */
		   char *vname,          /* Variable name            */
		   char *tag,            /* Input file tag           */
		   double fill,          /* Fill value               */
		   char *i_rule          /* Interpolation type       */
		 )
{
  geometry_t *geom = master->sgrid;
  int nce1 = geom->nce1;
  int nce2 = geom->nce2;
  int c, cc;
  int i, j, n;
  int number, nvals;
  char buf[MAXSTRLEN];
  double val, **d2;
  double *d1 = NULL;
  int *vec, nvec;
  int size;
  parameters_t *params=master->params;

  size = geom->sgsizS;
  nvec = geom->b2_t;
  vec = geom->w2_t;

  /*-----------------------------------------------------------------*/
  /* Initialise to the fill value                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    ret[c] = fill;
  }

  /*-----------------------------------------------------------------*/
  /* First check if this variable is contained in tracerdata         */
  if (strlen(fname) == 0 && strlen(master->tracerdata)) {
    timeseries_t *ts = NULL;
    int id;
    ts = hd_ts_read(master, master->tracerdata, 0);
    id = ts_get_index(ts, fv_get_varname(master->tracerdata, vname, buf));
    if (id >= 0) {
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	ret[c] = ts_eval_xy(ts, id, master->t,
			    geom->cellx[c],
			    geom->celly[c]);
      }
      if (master->trfilter & TRF_FILL2D)
	if ((n = tracer_find_index(vname, master->ntrS, master->trinfo_2d)) >= 0)
	  tracer_fill(master, ts, fname, ret, geom->sgsizS, vec, nvec, &master->trinfo_2d[n], fill);
      hd_ts_free(master, ts);
      return;
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }

  if (!strlen(fname)) {
    hd_warn("value_init_2d: No information provided to populate array %s; using defaults\n", vname);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check for regionalized formats                                  */
  if (i_rule && strlen(i_rule))
    sprintf(buf, "[var=%s] [i_rule=%s] %s", vname, i_rule, fname);
  else
    sprintf(buf, "[var=%s] %s\n", vname, fname);
  if ((n = set_variable(master, buf, ret, NULL))) {
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check if it is a region file                                    */
  if (value_init_regions(master, fname, ret, 2))
    return;

  /*-----------------------------------------------------------------*/
  /* Check if it is a list of values in the parameter file           */
  if (sscanf(fname, "%lf", &val) == 1) {
    val = atof(fname);
    /* Check if it is a list of values in the parameter file         */
    if (params->us_type & US_IJ && (int)val == nce1 * nce2) {
      if (prm_read_darray(fp, tag, &d1, &nvals) > 0) {
	if (nvals == nce1 * nce2) {
	  d2 = d_alloc_2d(nce1, nce2);
	  n = 0;
	  for (j = 0; j < nce2; j++)
	    for (i = 0; i < nce1; i++)
	      d2[j][i] = d1[n++];
	  c2s_2d(geom, ret, d2, nce1, nce2);
	  d_free_2d(d2);
	  d_free_1d(d1);
	} else {
	  hd_warn("value_init_2d: incorrect number specified for %s variable.\n", vname);
	}
      }
    } else if ((int)val == params->ns2) {
      if (prm_read_darray(fp, tag, &d1, &nvals) > 0) {
	if (nvals == params->ns2) {
	  for (i = 1; i <= params->ns2; i++)
	    ret[i] = d1[i-1];
	  d_free_1d(d1);
	} else {
	  hd_warn("value_init_2d: incorrect number specified for %s variable.\n", vname);
	}
      }
    }
  } else {
    /* Not a number nor so assume it is a file                       */
    timeseries_t *ts;
    int id; 

    if (strlen(i_rule))
      ts = hd_ts_read_us(master, fname, 0, i_rule);
    else
      ts = hd_ts_read(master, fname, 0);

    /* Infer if this is a plain ascii file (as opposed to a          */
    /* timeseries file) by the absence of the allocation of the time */
    /* field.                                                        */
    if (ts->t == NULL) {    /* Plain file                            */	
      GRID_SPECS *gs = NULL;
      double *x = NULL, *y = NULL, *z = NULL;
      int ii;

      /* Find and assign the needed variables                        */
      for (ii=0; ii<ts->df->nv; ii++) {
	if (strcmp("lon", ts->df->variables[ii].name) == 0)
	  x = ts->df->variables[ii].data;
	else if (strcmp("lat", ts->df->variables[ii].name) == 0)
	  y = ts->df->variables[ii].data;
	else if (strcmp(vname, ts->df->variables[ii].name) == 0)
	  z = ts->df->variables[ii].data;
      }

      /* Error checking                                              */
      if (x == NULL)
	hd_quit("val_init: 'lon' column not found in '%s'\n", fname);
      if (y == NULL)
	hd_quit("val_init: 'lat' column not found in '%s'\n", fname);
      if (z == NULL)
	hd_quit("val_init: '%s' column not found in '%s'\n", vname, fname);
      
      /* Initiliase the grid specs struct                            */
      gs = grid_interp_init(x, y, z, ts->df->dimensions[0].size, i_rule);
	  
      /* Do the interpolation                                        */
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	ret[c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
	    
	/* Check for nan's                                           */
	if (isnan(ret[c])) ret[c] = fill;
      }
	  
      /* Cleanup                                                     */
      grid_specs_destroy(gs);
    } else {
      /* Sanity check to see if the variable is in this datafile     */
      id = ts_get_index(ts, fv_get_varname(fname, vname, buf));

      if (id < 0)
	hd_warn("val_init_3d: The file '%s' does not contain the tracer '%s'\n", fname, vname);
      else {

	/* Use a sparse interpolation for 'standard' files           */
	if (value_init_sparse2d(master, ret, fname, vname, i_rule,
				schedule->start_time, NULL, ts) == 1) {
	  /* Otherwise interpolate as a gridded file                 */
	  for (cc = 1; cc <= nvec; cc++) {
	    c = vec[cc];
	    ret[c] = ts_eval_xy(ts, id, master->t,
				geom->cellx[c],
				geom->celly[c]);
	  }
	}
      }
      if (master->trfilter & TRF_FILL2D) {
	if ((n = tracer_find_index(vname, master->ntrS, master->trinfo_2d)) >= 0)
	  tracer_fill(master, ts, fname, ret, geom->sgsizS, vec, nvec, &master->trinfo_2d[n], fill);
	else
	  tracer_fill(master, ts, fname, ret, geom->sgsizS, vec, nvec, NULL, fill);
      }
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }
}

/* END value_init_2d()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise a sediment array from a number, array or    */
/* file input.                                                       */
/*-------------------------------------------------------------------*/
void value_init_sed(master_t *master,      /* Master data            */
		    double **ret,          /* 2D tracer to fill      */
		    FILE *fp,              /* File pointer           */
		    char *fname,           /* File name              */
		    char *vname,           /* Variable name          */
		    char *tag,             /* Input file tag         */
		    double fill,           /* Fill value             */
		    char *i_rule           /* Interpolation type     */
		    )
{
  geometry_t *geom = master->sgrid;
  int nce1 = geom->nce1;
  int nce2 = geom->nce2;
  int nz = geom->sednz;
  int c, cc;
  int i, j, k, n;
  int number, nvals;
  char buf[MAXSTRLEN];
  double val, ***d2;
  double *d1 = NULL;
  int *vec, nvec;
  int size;

  size = geom->sgsizS;
  nvec = geom->b2_t;
  vec = geom->w2_t;

  /*-----------------------------------------------------------------*/
  /* Initialise to the fill value                                    */
  for (k = 0; k < nz; k++)
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];
      ret[k][c] = fill;
    }

  /*-----------------------------------------------------------------*/
  /* First check if this variable is contained in tracerdata         */
  if (strlen(master->tracerdata)) {
    timeseries_t *ts = NULL;
    int id;
    ts = hd_ts_read(master, master->tracerdata, 0);
    id = ts_get_index(ts, fv_get_varname(master->tracerdata, vname, buf));
    if (id >= 0) {
      for (k = 0; k < nz; k++) {
	for (cc = 1; cc <= nvec; cc++) {
	  c = vec[cc];
	  ret[k][c] = ts_eval_xyz(ts, id, master->t,
				  geom->cellx[c],
				  geom->celly[c],
				  geom->gridz_sed[k][c]);
	}
	if (master->trfilter & TRF_FILLSED) {
	  /* Strip off the _sed suffix */
	  memset(buf, 0, MAXSTRLEN);
	  strncpy(buf, vname, strlen(vname)-4);
	  if ((n = tracer_find_index(buf, master->nsed, master->trinfo_sed)) >= 0)
	    tracer_fill(master, ts, fname, ret[k], geom->sgsizS, vec, nvec, &master->trinfo_sed[n], fill);
	}
      }
      hd_ts_free(master, ts);
      return;
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }

  if (!strlen(fname)) {
    hd_warn("value_init_sed: No sediment information provided to populate array %s; using defaults\n", vname);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Check if it is a region file                                    */
  n = 1;
  for (k = 0; k < nz; k++) {
    if (value_init_regions(master, fname, ret[k], 2) == 0) {
      n = 0;
      break;
    }
  }
  if (n) return;

  /*-----------------------------------------------------------------*/
  /* Check if it is a list of values in the parameter file           */
  if (sscanf(fname, "%lf", &val) == 1) {
    val = atof(fname);
    /* Check if it is a list of values in the parameter file         */
    if ((int)val == nce1 * nce2 * nz) {
      if (prm_read_darray(fp, tag, &d1, &nvals) > 0) {
	if (nvals == nce1 * nce2 * nz) {
	  double ***d2 = d_alloc_3d(nce1, nce2, nz);
	  n = 0;
	  for (k = 0; k < nz; k++)
	    for (j = 0; j < nce2; j++)
	      for (i = 0; i < nce1; i++)
		d2[k][j][i] = d1[n++];
	  for (k = 0; k < nz; k++)
	    c2s_2d(geom, ret[k], d2[k], nce1, nce2);
	  d_free_3d(d2);
	  d_free_1d(d1);
	} else {
	  hd_quit("value_init_sed: incorrect number specified for %s variable.\n", vname);
	}
      }
    }
  } else {
    /* Not a number nor so assume it is a file                       */
    timeseries_t *ts;
    int id; 

    ts = hd_ts_read(master, fname, 0);
	
    /* Infer if this is a plain ascii file (as opposed to a          */
    /* timeseries file) by the absence of the allocation of the time */
    /* field.                                                        */
    if (ts->t == NULL) {    /* Plain file                            */	
      GRID_SPECS *gs = NULL;
      double *x = NULL, *y = NULL, *z = NULL;
      int ii;

      /* Find and assign the needed variables                        */
      for (ii = 0; ii < ts->df->nv; ii++) {
	if (strcmp("lon", ts->df->variables[ii].name) == 0)
	  x = ts->df->variables[ii].data;
	else if (strcmp("lat", ts->df->variables[ii].name) == 0)
	  y = ts->df->variables[ii].data;
	else if (strcmp(vname, ts->df->variables[ii].name) == 0)
	  z = ts->df->variables[ii].data;
      }

      /* Error checking                                              */
      if (x == NULL)
	hd_quit("value_init_sed: 'lon' column not found in '%s'\n", fname);
      if (y == NULL)
	hd_quit("value_init_sed: 'lat' column not found in '%s'\n", fname);
      if (z == NULL)
	hd_quit("value_init_sed: '%s' column not found in '%s'\n", vname, fname);

      /* Initiliase the grid specs struct                            */
      gs = grid_interp_init(x, y, z, ts->df->dimensions[0].size, i_rule);

      /* Do the interpolation                                        */
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	for (k = 0; k < nz; k++) {
	  ret[k][c] = grid_interp_on_point(gs, geom->cellx[c], geom->celly[c]);
	  
	  /* Check for nan's                                         */
	  if (isnan(ret[k][c])) ret[k][c] = fill;
	}
      }

      /* Cleanup                                                     */
      grid_specs_destroy(gs);
    } else {
      /* Sanity check to see if the variable is in this datafile     */
      id = ts_get_index(ts, fv_get_varname(fname, vname, buf));
	
      if (id < 0)
	hd_quit("val_init_sed: The file '%s' does not contain the tracer '%s'\n", fname, vname);
	  
      /* Interpolate                                                 */
      for (k = 0; k < nz; k++) {
	for (cc = 1; cc <= nvec; cc++) {
	  c = vec[cc];
	  ret[k][c] = ts_eval_xyz(ts, id, master->t,
				  geom->cellx[c],
				  geom->celly[c],
				  geom->gridz_sed[k][c]);
	}
	if (master->trfilter & TRF_FILLSED) {
	  /* Strip off the _sed suffix */
	  memset(buf, 0, MAXSTRLEN);
	  strncpy(buf, vname, strlen(vname)-4);
	  if ((n = tracer_find_index(vname, master->nsed, master->trinfo_sed)) >= 0)
	    tracer_fill(master, ts, fname, ret[k], geom->sgsizS, vec, nvec, &master->trinfo_sed[n], fill);
	  else
	    tracer_fill(master, ts, fname, ret[k], geom->sgsizS, vec, nvec, NULL, fill);
	}
      }
    }
    if (ts != NULL) hd_ts_free(master, ts);
  }
}

/* END value_init_sed()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds an autotracer in the global list by its name                */
/*-------------------------------------------------------------------*/
int find_autotracer_by_name(char *name)
{
  int i;
  for (i = 0; i < NAUTOTR; i++)
    if (strcmp(autotracerlist[i].name, name) == 0) break;
  return(i);
}

/* END find_autotracer_by_name()                                     */
/*-------------------------------------------------------------------*


/*-------------------------------------------------------------------*/
/* Copies an autotracer attributes to a tracer info structure given  */
/* a tracer name.                                                    */
/*-------------------------------------------------------------------*/
void copy_autotracer_by_name(char *name, 
			     tracer_info_t tr[], 
			     int ntr, 
			     int *n, 
			     double **tra, 
			     double **trp)
{
  int tn;
  if ((tn = tracer_find_index(name, ntr, tr)) < 0) {
    int i;
    for (i = 0; i < NAUTOTR; i++)
      if (strcmp(autotracerlist[i].name, name) == 0) break;
    tracer_copy(&tr[*n], &autotracerlist[i]);
    tr[*n].n = tr[*n].m = *n;
    if (tra != NULL && trp != NULL) *trp = tra[*n];
    *n += 1;
  } else
    duplicate_error(name, tn);
}

/* END copy_autotracer_by_name()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Same as copy_autotracer_by_name() except includes tracer          */
/* initialisation.                                                   */
/*-------------------------------------------------------------------*/
int set_autotracer_by_name(char *name, tracer_info_t tr[], int ntr, 
			    int *n, double **tra, double **trp, char *buf)
{
  int tn, sn = *n;
  if ((tn = tracer_find_index(name, ntr, tr) < 0)) {
    int i;
    for (i = 0; i < NAUTOTR; i++)
      if (strcmp(autotracerlist[i].name, name) == 0) break;
    tracer_copy(&tr[*n], &autotracerlist[i]);
    tr[*n].n = tr[*n].m = *n;
    tr_dataset(buf, &tr[*n], 0.0);
    if (tra != NULL && trp != NULL) *trp = tra[*n];
    *n += 1;
  } else
    duplicate_error(name, tn);
  return(sn);
}

/* END set_autotracer_by_name()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds an autotracer in the global list by its group key           */
/*-------------------------------------------------------------------*/
void find_autotracer_by_groupkey(char *name, int *tra, int *n)
{
  int i;

  *n = 0;
  for (i = 0; i < NAUTOTR; i++) {
    if (strcmp(autotracerlist[i].groupkey, name) == 0) {
      tra[*n] = i;
      *n += 1;
    }
  }
}

/* END find_autotracer_by_groupkey()                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Copies an autotracer attributes to a tracer info structure given  */
/* a tracer groupkey.                                                */
/*-------------------------------------------------------------------*/
int copy_autotracer_by_groupkey(char *name, tracer_info_t tr[], int ntr, int *n)
{
  int i, j;
  int nf, nt = 0;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  char buf[MAXSTRLEN];
  int tn, sn = *n;

  for (i = 0; i < NAUTOTR; i++) {
    if (strlen(autotracerlist[i].groupkey) == 0) continue;
    if (strcmp(autotracerlist[i].groupkey, "NONE") == 0) continue;
    strcpy(buf, autotracerlist[i].groupkey); 
    nf = parseline(buf, fields, MAXNUMARGS);
    for (j = 0; j < nf; j++) {
      if (strcmp(fields[j], name) == 0) {
	if ((tn = tracer_find_index(autotracerlist[i].name, ntr, tr)) < 0) {
	  tracer_copy(&tr[*n], &autotracerlist[i]);
	  tr[*n].n = tr[*n].m = *n;
	  *n += 1;
	} else
	  duplicate_error(name, tn);
      }
    }
  }
  return(sn);
}

int set_autotracer_by_groupkey(char *name, tracer_info_t tr[], int ntr, 
			       int *n, char *key)
{
  int i, j;
  int nf, nt = 0;
  char *fields[MAXSTRLEN * MAXNUMARGS];
  char buf[MAXSTRLEN];
  int tn, sn = *n;

  for (i = 0; i < NAUTOTR; i++) {
    if (strlen(autotracerlist[i].groupkey) == 0) continue;
    strcpy(buf, autotracerlist[i].groupkey); 
    nf = parseline(buf, fields, MAXNUMARGS);
    for (j = 0; j < nf; j++) {
      if (strcmp(fields[j], name) == 0) {
	if ((tn = tracer_find_index(autotracerlist[i].name, ntr, tr)) < 0) {
	  tracer_copy(&tr[*n], &autotracerlist[i]);
	  tr_dataset(key, &tr[*n], tr[*n].fill_value_wc);
	  tr[*n].n = tr[*n].m = *n;
	  *n += 1;
	} else
	  duplicate_error(name, tn);
      }
    }
  }
  return(sn);
}

/* END find_autotracer_by_groupkey()                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints an error if an autotracer is also explicitly defined.      */
/* Note: this can be over-ridden if the keyword and tracer are       */
/* included in tracer_setup().                                       */
/*-------------------------------------------------------------------*/
void duplicate_error(char *name, int tn)
{
  if (strcmp(name, "temp") == 0 || strcmp(name, "salt") == 0) return;
  hd_quit("ERROR autotracers: duplicate TRACER %s: remove from tracer list or remove autotracer keyword.\n",
	  name);
}

/* END duplicate_error()                                             */
/*-------------------------------------------------------------------*/

void read_tsfile(timeseries_t *ts, char *name, char *t_units)
{
  ts_read(name,ts);
  /* Check data time units */
  if( strcmp(ts->t_units,t_units) != 0 )
    ts_convert_time_units(ts,t_units);
}

void addatt(int cdfid, int varid, const char*name, const char *text)
{
  nc_put_att_text(cdfid, varid, name, strlen(text), text);
}

int face_dist(double *layers, int nz, double depth)
{
  int k;
  for (k=0; k<=nz; k++) {
    if (layers[k] > depth)
      break;
  }
  return(k-1);
}


int centre_dist_min(double *layers, int nz, double depth)
{
  double zd;
  double zmin;               /* Min / max variables                  */
  int k, kk;

  kk = 0;
  zmin = 1e10;
  for (k=0; k<nz; k++) {
    zd = 0.5 * (layers[k] + params->layers[k+1]);
    if (fabs(zd - depth) < zmin) {
      zmin = fabs(zd - depth);
      kk = k;
    }
  }
  return(k);
}

int face_dist_min(double *layers, int nz, double depth)
{
  double zd;
  double zmin;               /* Min / max variables                  */ 
  int k, kk;

  kk = 0;
  zmin = 1e10;
  for (k=0; k<nz; k++) {
    zd = layers[k];
    if (fabs(zd - depth) < zmin) {
      zmin = fabs(zd - depth);
      kk = k;
    }
  }
  return(k);
}


void trn_dataset(char *data, tracer_info_t *trinfo, int tn, int ntr, int atr, double **tr, double def)
{
  double val;
  int n;
  if (strlen(data)) {
    for (n = atr; n < ntr; n++) {              /* data is a tracer */
      if (strcmp(data, trinfo[n].name) == 0) {
	tr[tn] = tr[n];
	return;
      }
    }
    if (sscanf(data, "%lf", &val) == 1) {      /* data is a value */
      trinfo[tn].fill_value_wc = atof(data);
      sprintf(trinfo[tn].data, "%c", '\0');
    } else
      strcpy(trinfo[tn].data, data);           /* data is a file */
  } else {
    trinfo[tn].fill_value_wc = def;            /* Default value */
    sprintf(trinfo[tn].data, "%c", '\0');
  }
}

void trf_dataset(char *data, tracer_info_t *trinfo, int tn, int ntr, int atr, double **tr, double def)
{
  double val;
  int n;

  if (strlen(data)) {
    trinfo[tn].flag = -1;
    for (n = atr; n < ntr; n++) {              /* data is a tracer */
      if (strcmp(data, trinfo[n].name) == 0) {
	trinfo[tn].flag = n;
	return;
      }
    }
    if (sscanf(data, "%lf", &val) == 1) {      /* data is a value */
      trinfo[tn].fill_value_wc = atof(data);
      sprintf(trinfo[tn].data, "%c", '\0');
    } else
      strcpy(trinfo[tn].data, data);           /* data is a file */
  } else {
    trinfo[tn].fill_value_wc = def;            /* Default value */
    sprintf(trinfo[tn].data, "%c", '\0');
  }
  return;
}

void tr_dataset(char *data, tracer_info_t *tr, double def)
{
  double val;
  if (strlen(data)) {
    if (sscanf(data, "%lf", &val) == 1)
      tr->fill_value_wc = atof(data);
    else
      strcpy(tr->data, data);
  } else {
    tr->fill_value_wc = def;
    sprintf(tr->data, "%c", '\0');
  }
}

#define AUTO_REM 0x002
#define AUTO_ADD 0x004
#define AUTO_UPD 0x008

/*-------------------------------------------------------------------*/
/* Updates the autotracer list                                       */
/*-------------------------------------------------------------------*/
void update_autotracer(tracer_info_t tr[], int ntr, int *atr, int *mtr)
{
  int n, tn;
  char buf[MAXSTRLEN];
  int found, flag;

  for (n = 0; n < ntr; n++) {
    found = 0;
    strcpy(buf, tr[n].tag);
    if (strlen(buf)) {
      if (contains_token(buf, "auto_update") != NULL) flag = AUTO_UPD;
      if (contains_token(buf, "auto_remove") != NULL) flag = AUTO_REM;
      if (contains_token(buf, "auto_add") != NULL) flag = AUTO_ADD;
      for (tn = 0; tn < NAUTOTR; tn++) {
	if (strcmp(tr[n].name, autotracerlist[tn].name) == 0) {
	  if (flag & AUTO_UPD) {
	    tracer_copy(&autotracerlist[tn], &tr[n]);
	    atr[tn] |= AUTO_UPD;
	    hd_warn("update_autotracer: Updating auto tracer %s\n", tr[n].name);
	    found = AUTO_UPD;
	  }
	  if (flag & AUTO_REM) {
	    atr[tn] |= AUTO_REM;
	    found = AUTO_REM;
	    hd_warn("update_autotracer: Removing auto tracer %s\n", tr[n].name);
	  }
	}
      }
      if (flag & AUTO_UPD && !(found & AUTO_UPD))
	hd_warn("update_autotracer: Can't find auto tracer %s to update - 'auto_add' first.\n", tr[n].name);
      if (flag & AUTO_REM && !(found & AUTO_REM))
	hd_warn("update_autotracer: Can't find auto tracer %s to remove\n", tr[n].name);
      if (flag & AUTO_ADD) {
	mtr[n] |= AUTO_ADD;
	hd_warn("update_autotracer: Adding tracer %s to autotracers\n", tr[n].name);
      }
    }
  }
}

/* END update_autotracer()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes a C module containing the autotracer list                  */
/*-------------------------------------------------------------------*/
void write_autotracer(master_t *master)
{
  FILE *fp;
  int n, m, trt, tc;
  int ntr = NAUTOTR;
  long t;
  int *atr, *m3d, *m2d;
  char *autoname = "autotracer.c";
  char buf[MAXSTRLEN];

  /* Allocate                                                        */
  atr = i_alloc_1d(ntr);
  m3d = i_alloc_1d(master->ntr);
  m2d = i_alloc_1d(master->ntrS);
  memset(atr, 0, ntr * sizeof(int));
  memset(m3d, 0, master->ntr * sizeof(int));
  memset(m2d, 0, master->ntrS * sizeof(int));

  /* Update the autotracers and count the total                      */
  update_autotracer(master->trinfo_3d, master->ntr, atr, m3d);
  update_autotracer(master->trinfo_2d, master->ntrS, atr, m2d);
  trt = ntr;
  for (n = master->ntrS-1; n >= 0; n--) {
    if (m2d[n] & AUTO_ADD) trt += 1;
  }
  for (n = master->ntr-1; n >= 0; n--) {
    if (m3d[n] & AUTO_ADD) trt += 1;
  }
  for (n = ntr-1; n >= 0; n--) {
    if (atr[n] & AUTO_REM) trt -= 1;
  }

  /* Open the autotracer module                                      */
  if (strlen(master->autotrpath)) {
    if (endswith(master->autotrpath, "/"))
      sprintf(buf, "%s%s", master->autotrpath, autoname);
    else
      sprintf(buf, "%s/%s", master->autotrpath, autoname);
  } else
    sprintf(buf, "%s", autoname);
  fp = fopen(buf, "w");
  hd_warn("write_autotracer: writing autotracer module to %s\n", buf);

  /* Write the autotracer list                                       */
  time(&t);
  fprintf(fp,"/*\n");
  fprintf(fp," *\n");
  fprintf(fp," *  ENVIRONMENTAL MODELLING SUITE (EMS)\n");
  fprintf(fp," *\n");  
  fprintf(fp," *  File: model/hd-us/tracers/autotracer.c\n");
  fprintf(fp," *\n");
  fprintf(fp," *  Description:\n");
  fprintf(fp," *  Autotracer list.\n");
  fprintf(fp," *  \n");
  fprintf(fp," *  Copyright:\n");
  fprintf(fp," *  Copyright (c) 2018. Commonwealth Scientific and Industrial\n");
  fprintf(fp," *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights\n");
  fprintf(fp," *  reserved. See the license file for disclaimer and full\n");
  fprintf(fp," *  use/redistribution conditions.\n");
  fprintf(fp," *  \n");
  fprintf(fp," *  $Id: load_tracer.c 7380 2023-07-26 04:36:43Z her127 $\n", version, ctime(&t));
  fprintf(fp," *\n");
  fprintf(fp," */\n\n");

  fprintf(fp,"#include <stdio.h>\n");
  fprintf(fp,"#include <math.h>\n");
  fprintf(fp,"#include \"hd.h\"\n");
  fprintf(fp,"#include \"tracer.h\"\n\n");

  fprintf(fp, "int NAUTOTR = %d;\n", trt);
  fprintf(fp, "tracer_info_t autotracerlist[] = {\n");
  tc = 1;
  for (n = 0; n < ntr; n++) {
    if (atr[n] & AUTO_REM) continue;
    m = (tc == trt) ? 0 : 1;
    write_auto_atts(fp, &autotracerlist[n], m);
    tc++;
  }
  for (n = 0; n < master->ntr; n++) {
    if (m3d[n] & AUTO_ADD) {
      m = (tc == trt) ? 0 : 1;
      write_auto_atts(fp, &master->trinfo_3d[n], m);
      tc++;
    }
  }
  for (n = 0; n < master->ntrS; n++) {
    if (m2d[n] & AUTO_ADD) {
      m = (tc == trt) ? 0 : 1;
      write_auto_atts(fp, &master->trinfo_2d[n], m);
      tc++;
    }
  }
  fprintf(fp, "};\n");
  i_free_1d(atr);
  i_free_1d(m3d);
  i_free_1d(m2d);
  fclose(fp);
}

/* END write_autotracer()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes tracer attributes to the autotracerlist file               */
/*-------------------------------------------------------------------*/
void write_auto_atts(FILE *fp, tracer_info_t *tr, int mode)
{
  char tname[MAXSTRLEN];

  fprintf(fp, "  {\n");
  fprintf(fp, "    .name = \"%s\",\n", tr->name);
  fprintf(fp, "    .long_name = \"%s\",\n", tr->long_name);
  if (strlen(tr->std_name))
    fprintf(fp, "    .std_name = \"%s\",\n", tr->std_name);
  fprintf(fp, "    .units = \"%s\",\n", tr->units);
  fprintf(fp, "    .valid_range_wc[0] = %4.2e,\n", tr->valid_range_wc[0]);
  fprintf(fp, "    .valid_range_wc[1] = %4.2e,\n", tr->valid_range_wc[1]);
  fprintf(fp, "    .fill_value_wc = %4.2f,\n", tr->fill_value_wc);
  if (tr->type & SEDIMENT) {
    fprintf(fp, "    .fill_value_sed = %4.2f,\n", tr->fill_value_sed);
    fprintf(fp, "    .valid_range_sed[0] = %4.2e,\n", tr->valid_range_sed[0]);
    fprintf(fp, "    .valid_range_sed[1] = %4.2e,\n", tr->valid_range_sed[1]);
  }
  typename(tr->type, tname);
  fprintf(fp, "    .type = %s,\n", tname);
  /*fprintf(fp, "    .type = 0x%x,\n", tr->type);*/
  fprintf(fp, "    .inwc = %d,\n", tr->inwc);
  if (tr->insed)
    fprintf(fp, "    .insed = %d,\n", tr->insed);
  if (tr->dissol)
    fprintf(fp, "    .dissol = %d,\n", tr->dissol);
  if (tr->partic)
    fprintf(fp, "    .partic = %d,\n", tr->partic);
  fprintf(fp, "    .advect = %d,\n", tr->advect);
  fprintf(fp, "    .diffuse = %d,\n", tr->diffuse);
  fprintf(fp, "    .diagn = %d,\n", tr->diagn);
  if (strlen(tr->decay))
    fprintf(fp, "    .decay = \"%s\",\n", tr->decay);
  /*fprintf(fp, "    .m = %d,\n", tr->m);*/
  fprintf(fp, "    .m = -1,\n");
  if (strlen(tr->tag) && strcmp(tr->tag, "auto_add") != 0 && strcmp(tr->tag, "auto_update") != 0)
    fprintf(fp, "    .tag = \"%s\",\n", tr->tag);
  if (tr->flag)
    fprintf(fp, "    .flag = %d,\n", tr->flag);
  if (strlen(tr->data))
    fprintf(fp, "    .data = \"%s\",\n", tr->data);
  if (strlen(tr->i_rule))
    fprintf(fp, "    .i_rule = \"%s\",\n", tr->i_rule);
  if (tr->scale)
    fprintf(fp, "    .scale = %5.2f,\n", tr->scale);
  if (strlen(tr->relax_file)) {
    fprintf(fp, "    .relax_file = \"%s\",\n", tr->relax_file);
    fprintf(fp, "    .relax_dt = \"%s\",\n", tr->relax_dt);
    fprintf(fp, "    .r_rate = \"%s\",\n", tr->r_rate);
    fprintf(fp, "    .tctype = \"%d\",\n", tr->tctype);
  }
  if (strlen(tr->reset_file)) {
    fprintf(fp, "    .reset_file = \"%s\",\n", tr->reset_file);
    fprintf(fp, "    .reset_dt = \"%s\",\n", tr->reset_dt);
  }
  if (strlen(tr->reset_interp))
    fprintf(fp, "    .reset_interp = \"%s\",\n", tr->reset_interp);
  if (strlen(tr->vector_name))
    fprintf(fp, "    .vector_name = \"%s\",\n", tr->vector_name);
  if (strlen(tr->vector_components))
    fprintf(fp, "    .vector_components = \"%s\",\n", tr->vector_components);
  if (strlen(tr->tracerstat))
    fprintf(fp, "    .tracerstat = \"%s\",\n", tr->tracerstat);
  if (strlen(tr->trstat_tag))
    fprintf(fp, "    .dt = \"%s\",\n", tr->trstat_tag);
  if (strlen(tr->groupkey))
    fprintf(fp, "    .groupkey = \"%s\"\n", tr->groupkey);
  else
    fprintf(fp, "    .groupkey = \"NONE\"\n");
  if (mode)
    fprintf(fp, "  },\n");
  else
    fprintf(fp, "  }\n");
}

/* END write_auto_atts()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Converts a tracer type to a string                                */
/*-------------------------------------------------------------------*/
void typename(int code, char *name)
{
  int found = 0;

  if (code & WATER) {
    strcpy(name, "WATER");
    found = 1;
  }
  if (code & SEDIM) {
    if (found)
      sprintf(name, "%s|SEDIM", name);
    else {
      strcpy(name, "SEDIM");
      found = 1;
    }
  }
  if (code & INTER) {
    if (found)
      sprintf(name, "%s|INTER", name);
    else {
      strcpy(name, "INTER");
      found = 1;
    }
  }
  if (code & HYDRO)
    sprintf(name, "%s|HYDRO", name);
  if (code & SEDIMENT)
    sprintf(name, "%s|SEDIMENT", name);
  if (code & ECOLOGY)
    sprintf(name, "%s|ECOLOGY", name);
  if (code & WAVE)
    sprintf(name, "%s|WAVE", name);
  if (code & TRACERSTAT)
    sprintf(name, "%s|TRACERSTAT", name);
  if (code & PROGNOSTIC)
    sprintf(name, "%s|PROGNOSTIC", name);
  if (code & DIAGNOSTIC)
    sprintf(name, "%s|DIAGNOSTIC", name);
  if (code & PARAMETER)
    sprintf(name, "%s|PARAMETER", name);
  if (code & FORCING)
    sprintf(name, "%s|FORCING", name);
  if (code & E1VAR)
    sprintf(name, "%s|E1VAR", name);
  if (code & E2VAR)
    sprintf(name, "%s|E2VAR", name);
  if (code & CLOSURE)
    sprintf(name, "%s|CLOSURE", name);
  if (code & OPTICAL)
    sprintf(name, "%s|OPTICAL", name);

}

/* END typename()                                                    */
/*-------------------------------------------------------------------*/
