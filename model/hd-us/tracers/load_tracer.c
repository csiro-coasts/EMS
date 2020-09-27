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
 *  $Id: load_tracer.c 6476 2020-02-18 23:53:31Z her127 $
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
int value_init_sparse2d(master_t *master, double *ret, char *fname,
			char *vname, char *i_rule);
int value_init_sparse3d(master_t *master, double *ret, char *fname,
			char *vname, char *i_rule);
void filter_tracer(tracer_info_t *trinfo, geometry_t *window, double **tr, int n);

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


  int t;
  char buf[MAXSTRLEN];

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
  int c, cc;
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
}

/* END scale_tracer()                                               */
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Filters a tracer                                                 */
/*------------------------------------------------------------------*/
void filter_tracer(tracer_info_t *trinfo, geometry_t *geom, double **tr, int n)
{
  int c, cc, cp, cm, i, m;
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
    
    if (value_init_sparse3d(master, master->tr_wc[t], buf, 
			    master->trinfo_3d[t].name, 
			    master->trinfo_3d[t].i_rule) == 0)
      continue;

    ts = hd_ts_read(master, buf, 0);
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
 

  int t;

  char buf2[MAXSTRLEN];
  char tag[MAXSTRLEN];
  geometry_t *geom = master->sgrid;

  /* Load up the default tracer values */
  load_wc_tracer_defaults_2d(master);

  prm_set_errfn(hd_silent_warn);
  for (t = 0; t < master->ntrS; ++t) {
    int tm = master->trinfo_2d[t].m;

    sprintf(tag, "TRACER%1.1d.interp_type", tm);
    if (!(prm_read_char(fp, tag, buf2)))
      buf2[0] = '\0';
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

  if (dim == WATER)
    value_init_3d(master, master->tr_wc[trn], fp, fname, vname, tag, fill, i_rule);
  else if (dim == INTER)
    value_init_2d(master, master->tr_wcS[trn], fp, fname, vname, tag, fill, i_rule);
  else if (dim == SEDIM)
    value_init_sed(master, master->tr_sed[trn], fp, fname, vname, tag, fill, i_rule);

  prm_set_errfn(hd_quit);
}

/* END load_wc_tracer_name()                                        */
/*------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 2D tracers in the master                */
/*-------------------------------------------------------------------*/
void init_tracer_2d(parameters_t *params, /* Input parameters data   */
                    master_t *master      /* Master data             */
  )
{
  int tn;
  geometry_t *geom = master->geom;
  char buf[MAXSTRLEN];
  tn = 0;

  /* Initialise */
  master->alert_a = master->alert_c = NULL;
  master->avhrr = master->ghrsst = master->ghrsste = master->shwin = master->layth = NULL;
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
  master->cellres = master->equitide =   master->tpxotide = NULL;
  master->tpxovelu = master->tpxovelv = master->tpxotranu = master->tpxotranv = NULL;
  master->uat = master->vat = master->meshun = NULL;

  /* Waves */
  master->ustrcw = master->wave_amp = master->wave_Cd = NULL;
  master->wave_ub = master->wave_period = master->wave_dir = NULL;
  master->wave_Sxy = master->wave_Syx = master->vol_cons = NULL;
  master->wave_Fx = master->wave_Fy = NULL;
  master->wave_ste1 = master->wave_ste2 = NULL;
  master->tau_w1 = master->tau_w2 = master->tau_diss1 = master->tau_diss2 = NULL;
  master->sep = master->bep = master->tfront = NULL;

  /* SWR */
  master->swr_attn1 = master->swr_tran = NULL;
  master->swr_babs = master->swreg = master->swrms = master->attn_mean = master->tran_mean = NULL;
  if (params->swr_type & SWR_2D) master->swr_attn = NULL;

  /* Assign the 2D pointers.  */
  sprintf(buf, "%c", '\0');
  /* CFL time-step diagnostic.  */
  if (!(params->cfl & NONE)) {
    master->cfl2d = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "cfl2d");
    strcpy(master->trinfo_2d[tn].long_name, "2D CFL timestep");
    strcpy(master->trinfo_2d[tn].units, "s");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = 0;
    master->cfl3d = master->tr_wcS[tn + 1];
    strcpy(master->trinfo_2d[tn + 1].name, "cfl3d");
    strcpy(master->trinfo_2d[tn + 1].long_name, "3D CFL timestep");
    strcpy(master->trinfo_2d[tn + 1].units, "s");
    master->trinfo_2d[tn + 1].type = INTER;
    tr_dataset(buf, &master->trinfo_2d[tn + 1], 0.0);
    master->trinfo_2d[tn + 1].valid_range_wc[0] = 0;
    master->trinfo_2d[tn + 1].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn + 1].n = 1;
    master->cour = master->tr_wcS[tn + 2];
    strcpy(master->trinfo_2d[tn + 2].name, "courant");
    strcpy(master->trinfo_2d[tn + 2].long_name, "Courant timestep");
    strcpy(master->trinfo_2d[tn + 2].units, "s");
    master->trinfo_2d[tn + 2].type = INTER;
    tr_dataset(buf, &master->trinfo_2d[tn + 2], 0.0);
    master->trinfo_2d[tn + 2].valid_range_wc[0] = 0;
    master->trinfo_2d[tn + 2].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn + 2].n = 1;
    master->lips = master->tr_wcS[tn + 3];
    strcpy(master->trinfo_2d[tn + 3].name, "lipschitz");
    strcpy(master->trinfo_2d[tn + 3].long_name, "Lipschitz timestep");
    strcpy(master->trinfo_2d[tn + 3].units, "s");
    master->trinfo_2d[tn + 3].type = INTER;
    tr_dataset(buf, &master->trinfo_2d[tn + 3], 0.0);
    master->trinfo_2d[tn + 3].valid_range_wc[0] = 0;
    master->trinfo_2d[tn + 3].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn + 3].n = 1;
    master->ahsb = master->tr_wcS[tn + 4];
    strcpy(master->trinfo_2d[tn + 4].name, "diffstab");
    strcpy(master->trinfo_2d[tn + 4].long_name, "Horizontal diffusion timestep");
    strcpy(master->trinfo_2d[tn + 4].units, "s");
    master->trinfo_2d[tn + 4].type = INTER;
    tr_dataset(buf, &master->trinfo_2d[tn + 4], 0.0);
    master->trinfo_2d[tn + 4].valid_range_wc[0] = 0;
    master->trinfo_2d[tn + 4].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn + 4].n = 1;
    master->courn = master->tr_wcS[tn + 5];
    strcpy(master->trinfo_2d[tn + 5].name, "courant_no");
    strcpy(master->trinfo_2d[tn + 5].long_name, "Courant Number");
    strcpy(master->trinfo_2d[tn + 5].units, "");
    master->trinfo_2d[tn + 5].type = INTER;
    tr_dataset(buf, &master->trinfo_2d[tn + 5], 0.0);
    master->trinfo_2d[tn + 5].valid_range_wc[0] = 0;
    master->trinfo_2d[tn + 5].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn + 5].n = 1;
    tn = 6;
  }
  /* Mixed layer depth diagnostic */
  if (!(params->mixlayer & NONE)) {
    master->mixl = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "mixed_layer");
    strcpy(master->trinfo_2d[tn].long_name, "Mixed layer depth");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;

    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* Steric height diagnostic */
  if (params->lnm != 0.0) {
    master->steric = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "steric");
    strcpy(master->trinfo_2d[tn].long_name, "Steric height");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;

    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* Vorticity diagnostic */
  if (params->vorticity & ABSOLUTE) {
    master->av = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "abs_vor");
    strcpy(master->trinfo_2d[tn].long_name, "Absolute vorticity");
    strcpy(master->trinfo_2d[tn].units, "s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;

    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->vorticity & RELATIVE) {
    master->rv = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rel_vor");
    strcpy(master->trinfo_2d[tn].long_name, "Relative vorticity");
    strcpy(master->trinfo_2d[tn].units, "s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;

    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->vorticity & POTENTIAL) {
    master->pv = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "pot_vor");
    strcpy(master->trinfo_2d[tn].long_name, "Potential vorticity");
    strcpy(master->trinfo_2d[tn].units, "m-1s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->vorticity & TENDENCY) {
    master->rv_drvdt = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_drvdt");
    strcpy(master->trinfo_2d[tn].long_name, "RV rate of change tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_nonlin = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_nonlin");
    strcpy(master->trinfo_2d[tn].long_name, "RV nonlinear tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_beta = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_beta");
    strcpy(master->trinfo_2d[tn].long_name,
           "RV planetary vorticity tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_strch = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_strch");
    strcpy(master->trinfo_2d[tn].long_name,
           "RV topographic vorticity tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_jebar = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_jebar");
    strcpy(master->trinfo_2d[tn].long_name, "RV JEBAR tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_wsc = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_wsc");
    strcpy(master->trinfo_2d[tn].long_name, "RV wind stress curl tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->rv_bsc = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rv_bsc");
    strcpy(master->trinfo_2d[tn].long_name,
           "RV bottom stress curl tendency");
    strcpy(master->trinfo_2d[tn].units, "s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* Numbers diagnostic */
  if (params->numbers & ROSSBY_EX) {
    master->rossby_ex = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "rossby_external");
    strcpy(master->trinfo_2d[tn].long_name,
           "External Rossby radius");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SPEED_2D) {
    master->speed_2d = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "current_speed_2d");
    strcpy(master->trinfo_2d[tn].long_name,
           "Current Speed 2D");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SPEED_SQ) {
    master->speed_sq = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "speed_sq");
    strcpy(master->trinfo_2d[tn].long_name,
           "Bottom Current^2");
    strcpy(master->trinfo_2d[tn].units, "m2s-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & OBC_PHASE) {
    master->obc_phase = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "obc_phase");
    strcpy(master->trinfo_2d[tn].long_name,
           "OBC phase speed");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & WIND_CD) {
    master->wind_Cd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wind_Cd");
    strcpy(master->trinfo_2d[tn].long_name, 
           "Momentum drag coefficient");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    tn++;
  }
  if (params->numbers & CELLRES) {
    master->cellres = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "cell_resolution");
    strcpy(master->trinfo_2d[tn].long_name, "Mean Cell Resolution");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers1 & MESHUN) {
    master->meshun = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "mesh_uniformity");
    strcpy(master->trinfo_2d[tn].long_name, "Mesh uniformity index");
    strcpy(master->trinfo_2d[tn].units, "%");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->waves & TAN_RAD && params->tendf) {
    master->u1_rad = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "u1_rad");
    strcpy(master->trinfo_2d[tn].long_name,
           "u1 radiation stress tendency");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|WAVE|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->u2_rad = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "u2_rad");
    strcpy(master->trinfo_2d[tn].long_name,
           "u2 radiation stress tendency");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|WAVE|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->means & ETA_M) {
    master->etam = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "eta_mean");
    strcpy(master->trinfo_2d[tn].long_name, "Mean elevation");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    memset(master->etam, 0, geom->szcS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->means & WIND) {
    master->w1m = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "w1mean");
    strcpy(master->trinfo_2d[tn].long_name, "Mean eastward wind stress");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    memset(master->w1m, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->w2m = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "w2mean");
    strcpy(master->trinfo_2d[tn].long_name, "Mean northward wind stress");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    memset(master->w2m, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->means & VEL2D) {
    master->u1am = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "u1av_mean");
    strcpy(master->trinfo_2d[tn].long_name, "Mean 2D east velocity ");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    strcpy(master->trinfo_2d[tn].vector_name, "Mean 2D current");
    strcpy(master->trinfo_2d[tn].vector_components, "u1av_mean u2av_mean");
    memset(master->u1am, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->u2am = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "u2av_mean");
    strcpy(master->trinfo_2d[tn].long_name, "Mean 2D north velocity");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    strcpy(master->trinfo_2d[tn].vector_name, "Mean 2D current");
    strcpy(master->trinfo_2d[tn].vector_components, "u1av_mean u2av_mean");
    memset(master->u2am, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->means & MTRA2D) {
    master->tram = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tracer_mean");
    sprintf(buf, "Mean %s", params->means_tra);
    strcpy(master->trinfo_2d[tn].long_name, buf);
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    memset(master->tram, 0, geom->sgsizS * sizeof(double));
    master->trinfo_2d[tn].n = tn;
    tn++;
  }

  if (params->heatflux & ADVANCED) {
    master->nhfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "nhf");
    strcpy(master->trinfo_2d[tn].long_name, "Net heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->swrd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr");
    strcpy(master->trinfo_2d[tn].long_name, "Short wave radiation");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;

    master->lwrd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "lwr");
    strcpy(master->trinfo_2d[tn].long_name, "Long wave radiation");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->lhfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "lhf");
    strcpy(master->trinfo_2d[tn].long_name, "Latent heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->shfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "shf");
    strcpy(master->trinfo_2d[tn].long_name, "Sensible heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }

  if (params->heatflux & INVERSE) {
    master->nhfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "nhf");
    strcpy(master->trinfo_2d[tn].long_name, "Net heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->heatflux & (COMP_HEAT | COMP_HEAT_MOM)) {
    master->nhfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "nhf");
    strcpy(master->trinfo_2d[tn].long_name, "Net heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE)) {
    if (params->heatflux & (COMP_HEAT_MOM | COMP_HEAT_NONE)) {
      /* See logic in heatflux.c:comp_heat_mom */
      master->swr = master->tr_wcS[tn];
    } else 
      master->swrd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr");
    strcpy(master->trinfo_2d[tn].long_name, "Short wave radiation");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->lwrn = tn;
    strcpy(master->trinfo_2d[tn].name, "lwr");
    strcpy(master->trinfo_2d[tn].long_name, "Long wave radiation");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->lhfn = tn;
    strcpy(master->trinfo_2d[tn].name, "lhf");
    strcpy(master->trinfo_2d[tn].long_name, "Latent heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->shfn = tn;
    strcpy(master->trinfo_2d[tn].name, "shf");
    strcpy(master->trinfo_2d[tn].long_name, "Sensible heat flux");
    strcpy(master->trinfo_2d[tn].units, "W m-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 1000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    if (params->heatflux & COMP_HEAT_NONE) {
      if (strlen(params->precip)) {
	master->precipn = tn;
	strcpy(master->trinfo_2d[tn].name, "precip");
	strcpy(master->trinfo_2d[tn].long_name, "Precipitation rate");
	strcpy(master->trinfo_2d[tn].units, "ms-1");
	master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
	tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
	master->trinfo_2d[tn].valid_range_wc[0] = 0;
	master->trinfo_2d[tn].valid_range_wc[1] = 2000;
	master->trinfo_2d[tn].n = tn;
	tn++;
      }
      if (strlen(params->evap)) {
	master->evapn = tn;
	strcpy(master->trinfo_2d[tn].name, "evap");
	strcpy(master->trinfo_2d[tn].long_name, "Evaporation rate");
	strcpy(master->trinfo_2d[tn].units, "ms-1");
	master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
	tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
	master->trinfo_2d[tn].valid_range_wc[0] = 0;
	master->trinfo_2d[tn].valid_range_wc[1] = 2000;
	master->trinfo_2d[tn].n = tn;
	tn++;
      }
    }
  }
  if (params->heatflux & NET_HEAT) {
    master->nhfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "nhf");
    strcpy(master->trinfo_2d[tn].long_name, "Net heat flux");
    strcpy(master->trinfo_2d[tn].units, "Wm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1000;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->swrd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr");
    strcpy(master->trinfo_2d[tn].long_name, "Short wave radiation");
    strcpy(master->trinfo_2d[tn].units, "Wm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->saltflux & (ADVANCED | BULK | ORIGINAL)) {
    master->nsfd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "nsf");
    strcpy(master->trinfo_2d[tn].long_name, "Net salt flux");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->saltflux & (ADVANCED | ORIGINAL)) {
    master->evapn = tn;
    strcpy(master->trinfo_2d[tn].name, "evap");
    strcpy(master->trinfo_2d[tn].long_name, "Evaporation rate");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->precipn = tn;
    strcpy(master->trinfo_2d[tn].name, "precip");
    strcpy(master->trinfo_2d[tn].long_name, "Precipitation rate");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 2000;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->waves & BOT_STR) {
    master->wave_Cd = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_Cd");
    strcpy(master->trinfo_2d[tn].long_name, "Wave enhanced bottom drag");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->waves & TAN_RAD) {
    master->wave_Sxy = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_Sxy");
    strcpy(master->trinfo_2d[tn].long_name, "Radiation stress, x tangential");
    strcpy(master->trinfo_2d[tn].units, "m2s-2");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_Syx = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_Syx");
    strcpy(master->trinfo_2d[tn].long_name, "Radiation stress, y tangential");
    strcpy(master->trinfo_2d[tn].units, "m2s-2");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->waves & WAVE_FOR) {
    master->wave_Fx = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_Fx");
    strcpy(master->trinfo_2d[tn].long_name, "Wave-induced force along e1");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_Fy = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_Fy");
    strcpy(master->trinfo_2d[tn].long_name, "Wave-induced force along e2");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->waves & (STOKES|SPECTRAL)) {
    master->wave_ste1 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_ste1");
    strcpy(master->trinfo_2d[tn].long_name, "Stokes velocity along e1");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -100;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_ste2 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_ste2");
    strcpy(master->trinfo_2d[tn].long_name, "Stokes velocity along e2");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -100;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    if (params->waves & STOKES_DRIFT) {
      master->tau_w1 = master->tr_wcS[tn];
      strcpy(master->trinfo_2d[tn].name, "tau_w1");
      strcpy(master->trinfo_2d[tn].long_name, "Wave supported wind stress e1 direction");
      strcpy(master->trinfo_2d[tn].units, "Nm-2");
      master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
      tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
      master->trinfo_2d[tn].valid_range_wc[0] = -100;
      master->trinfo_2d[tn].valid_range_wc[1] = 100;
      master->trinfo_2d[tn].n = tn;
      tn++;
      master->tau_w2 = master->tr_wcS[tn];
      strcpy(master->trinfo_2d[tn].name, "tau_w2");
      strcpy(master->trinfo_2d[tn].long_name, "Wave supported wind stress e2 direction");
      strcpy(master->trinfo_2d[tn].units, "Nm-2");
      master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
      tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
      master->trinfo_2d[tn].valid_range_wc[0] = -100;
      master->trinfo_2d[tn].valid_range_wc[1] = 100;
      master->trinfo_2d[tn].n = tn;
      tn++;
      master->tau_diss1 = master->tr_wcS[tn];
      strcpy(master->trinfo_2d[tn].name, "tau_diss1");
      strcpy(master->trinfo_2d[tn].long_name, "Wave to ocean stress e1 direction");
      strcpy(master->trinfo_2d[tn].units, "Nm-2");
      master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
      tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
      master->trinfo_2d[tn].valid_range_wc[0] = -100;
      master->trinfo_2d[tn].valid_range_wc[1] = 100;
      master->trinfo_2d[tn].n = tn;
      tn++;
      master->tau_diss2 = master->tr_wcS[tn];
      strcpy(master->trinfo_2d[tn].name, "tau_diss2");
      strcpy(master->trinfo_2d[tn].long_name, "Wave to ocean stress e2 direction");
      strcpy(master->trinfo_2d[tn].units, "Nm-2");
      master->trinfo_2d[tn].type = INTER|WAVE|FORCING;
      tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
      master->trinfo_2d[tn].valid_range_wc[0] = -100;
      master->trinfo_2d[tn].valid_range_wc[1] = 100;
      master->trinfo_2d[tn].n = tn;
      tn++;
    }
  }
  if (params->orbital) {
    master->ustrcw = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "ustrcw");
    strcpy(master->trinfo_2d[tn].long_name, "Wave current friction");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|WAVE|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_ub = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_ub");
    strcpy(master->trinfo_2d[tn].long_name, "Wave orbital velocity");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|WAVE|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_period = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_period");
    strcpy(master->trinfo_2d[tn].long_name, "Wave period");
    strcpy(master->trinfo_2d[tn].units, "s");
    master->trinfo_2d[tn].type = INTER|WAVE|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_dir = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_dir");
    strcpy(master->trinfo_2d[tn].long_name, "Wave direction");
    strcpy(master->trinfo_2d[tn].units, "deg");
    master->trinfo_2d[tn].type = INTER|WAVE|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->wave_amp = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wave_amp");
    strcpy(master->trinfo_2d[tn].long_name, "Wave amplitude");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].type = INTER|WAVE|PROGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->etarlx & (RELAX|ALERT|BOUNDARY)) {
    /*master->eta_rlx->val1 = master->tr_wcS[tn];*/
    strcpy(master->trinfo_2d[tn].name, "oeta");
    strcpy(master->trinfo_2d[tn].long_name, "Relaxation elevation");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].fill_value_wc = 0.0;
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].diagn = 1;
    master->trinfo_2d[tn].valid_range_wc[0] = -1e35;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e35;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->etarlx & ETA_TPXO) {
    /*
    strcpy(master->trinfo_2d[tn].name, "oeta");
    strcpy(master->trinfo_2d[tn].long_name, "Relaxation elevation");
    strcpy(master->trinfo_2d[tn].units, "metre");
    master->trinfo_2d[tn].fill_value_wc = 0.0;
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].diagn = 1;
    master->trinfo_2d[tn].valid_range_wc[0] = -1e35;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e35;
    master->trinfo_2d[tn].n = tn;
    tn++;
    */
    master->eta_tc = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "eta_tc");
    strcpy(master->trinfo_2d[tn].long_name, "Relaxation eta time constant");
    strcpy(master->trinfo_2d[tn].units, "days");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* AVHRR SST */
  if (params->avhrr) {
    master->avhrr = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "AVHRR");
    strcpy(master->trinfo_2d[tn].long_name, "AVHRR SST");
    strcpy(master->trinfo_2d[tn].units, "degrees C");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    strcpy(master->trinfo_2d[tn].reset_file, "avhrr_list.txt(AVHRR=sst)");
    strcpy(master->trinfo_2d[tn].reset_dt, "1 day");
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* GHRSST SST */
  if (params->ghrsst) {
    master->ghrsst = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "ghrsst");
    strcpy(master->trinfo_2d[tn].long_name, "GHRSST L4 SST");
    strcpy(master->trinfo_2d[tn].units, "degrees C");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    strcpy(master->trinfo_2d[tn].reset_file, params->ghrsst_path);
    strcpy(master->trinfo_2d[tn].reset_dt, "1 day");
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->ghrsste = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "ghrsst_error");
    strcpy(master->trinfo_2d[tn].long_name, "GHRSST L4 SST error");
    strcpy(master->trinfo_2d[tn].units, "degrees C");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    strcpy(master->trinfo_2d[tn].reset_file, params->ghrsst_path);
    strcpy(master->trinfo_2d[tn].reset_dt, "1 day");
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->show_win) {
    master->shwin = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "windows");
    strcpy(master->trinfo_2d[tn].long_name, "window partitions");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (strlen(params->bathystats)) {
    master->bathy_range_min = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "bathy_range_min");
    strcpy(master->trinfo_2d[tn].long_name, "Bathymetry deviation above");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->bathy_range_max = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "bathy_range_max");
    strcpy(master->trinfo_2d[tn].long_name, "Bathymetry deviation below");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->bathy_grad_min = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "bathy_grad");
    strcpy(master->trinfo_2d[tn].long_name, "Bathymetry gradient");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->bathy_grad_max = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "bathy_grad_max");
    strcpy(master->trinfo_2d[tn].long_name, "Bathymetry maximum gradient deviation");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (contains_token(params->alert, "ACTIVE") != NULL) {
    master->alert_a = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "alerts_actual");
    strcpy(master->trinfo_2d[tn].long_name, "alerts; actual");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->alert_c = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "alerts_cumulative");
    strcpy(master->trinfo_2d[tn].long_name, "alerts; cumulative");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e36;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->u1vhin = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "U1VH0");
    strcpy(master->trinfo_2d[tn].long_name, "Initial e1 horizontal viscosity");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e36;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->u2vhin = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "U2VH0");
    strcpy(master->trinfo_2d[tn].long_name, "Initial e2 horizontal viscosity");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e36;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->fillf & (WEIGHTED|MONOTONIC) || (params->tmode & SP_FFSL)) {
    master->vol_cons = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "vol_cons");
    strcpy(master->trinfo_2d[tn].long_name, "Volume conservation");
    strcpy(master->trinfo_2d[tn].units, "%");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SOUND) {
    master->sonic = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "sonic_depth");
    strcpy(master->trinfo_2d[tn].long_name, "Sonic depth");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 100;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & EKPUMP) {
    master->sep = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "sep");
    strcpy(master->trinfo_2d[tn].long_name, "Surface Ekman pumping");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->bep = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "bep");
    strcpy(master->trinfo_2d[tn].long_name, "Bottom Ekman pumping");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & TIDEFR) {
    master->tfront = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tide_front");
    strcpy(master->trinfo_2d[tn].long_name, "Tide front");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & WET_CELLS) {
    master->wetcell = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "wet_cells");
    strcpy(master->trinfo_2d[tn].long_name, "Wet cell diagnostic");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SURF_LAYER) {
    master->surfz = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "surf_layer");
    strcpy(master->trinfo_2d[tn].long_name, "k index of surface layer");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SLOPE) {
    master->slope_x = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "surf_slope");
    strcpy(master->trinfo_2d[tn].long_name, "Mean surface slope");
    strcpy(master->trinfo_2d[tn].units, "m/m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers & BOTSTRESS) {
    master->tau_be1 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tau_be1");
    strcpy(master->trinfo_2d[tn].long_name, "Bottom stress in east direction");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    strcpy(master->trinfo_2d[tn].vector_name, "Bottom stress");
    strcpy(master->trinfo_2d[tn].vector_components, "tau_be1 tau_be2");
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->tau_be2 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tau_be2");
    strcpy(master->trinfo_2d[tn].long_name, "Bottom stress in north direction");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    strcpy(master->trinfo_2d[tn].vector_name, "Bottom stress");
    strcpy(master->trinfo_2d[tn].vector_components, "tau_be1 tau_be2");
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->tau_bm = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tau_bm");
    strcpy(master->trinfo_2d[tn].long_name, "Bottom stress magnitude");
    strcpy(master->trinfo_2d[tn].units, "Nm-2");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e4;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (strlen(params->swr_babs)) {
    master->swr_babs = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr_bot_absorb");
    strcpy(master->trinfo_2d[tn].long_name, "SWR bottom absorption");
    strcpy(master->trinfo_2d[tn].units, "");
    trn_dataset(params->swr_babs, master->trinfo_2d, tn, params->ntrS, params->atrS, master->tr_wcS, 1.0);
    /*tr_dataset(params->swr_babs, &master->trinfo_2d[tn], 1.0);*/
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->swr_type & SWR_2D && strlen(params->swr_attn)) {
    master->swr_attn = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr_attenuation");
    strcpy(master->trinfo_2d[tn].long_name, "SWR attenuation");
    strcpy(master->trinfo_2d[tn].units, "m-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    /*tr_dataset(params->swr_attn, &master->trinfo_2d[tn], 0.073);*/
    trn_dataset(params->swr_attn, master->trinfo_2d, tn, params->ntrS, params->atrS, master->tr_wcS, 0.073);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (strlen(params->swr_attn1)) {
    master->swr_attn1 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr_deep_attenuation");
    strcpy(master->trinfo_2d[tn].long_name, "SWR deep attenuation");
    strcpy(master->trinfo_2d[tn].units, "m-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(params->swr_attn1, &master->trinfo_2d[tn], 0.073);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (strlen(params->swr_tran)) {
    master->swr_tran = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swr_transmission");
    strcpy(master->trinfo_2d[tn].long_name, "SWR transmission");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    trn_dataset(params->swr_tran, master->trinfo_2d, tn, params->ntrS, params->atrS, master->tr_wcS, 0.26);
    /*tr_dataset(params->swr_tran, &master->trinfo_2d[tn], 0.26);*/
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
 if (strlen(params->swr_regions)) {
    master->swreg = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swreg");
    strcpy(master->trinfo_2d[tn].long_name, "SWR param estimation regions");
    strcpy(master->trinfo_2d[tn].units, "");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->swrms = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "swrms");
    strcpy(master->trinfo_2d[tn].long_name, "SWR param estimation RMSE");
    strcpy(master->trinfo_2d[tn].units, "degrees C");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->swrms = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "attn_mean");
    strcpy(master->trinfo_2d[tn].long_name, "SWR mean attenuation");
    strcpy(master->trinfo_2d[tn].units, "m-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->swrms = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tran_mean");
    strcpy(master->trinfo_2d[tn].long_name, "SWR mean transmission");
    strcpy(master->trinfo_2d[tn].units, " ");
    master->trinfo_2d[tn].type = INTER|HYDRO|PARAMETER;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->riverflow) {
    master->riverflow = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "flow");
    strcpy(master->trinfo_2d[tn].long_name, "River flow");
    strcpy(master->trinfo_2d[tn].units, "m3s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    if (params->riverflow == 2) {
      master->riverdepth = master->tr_wcS[tn];
      strcpy(master->trinfo_2d[tn].name, "flow_depth");
      strcpy(master->trinfo_2d[tn].long_name, "River flow depth");
      strcpy(master->trinfo_2d[tn].units, "m");
      master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
      tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
      master->trinfo_2d[tn].valid_range_wc[0] = 0;
      master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_2d[tn].n = tn;
      tn++;
    }
  }
  if (params->tidep) {
    master->equitide = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "equitide");
    strcpy(master->trinfo_2d[tn].long_name, "Equilibrium tide");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers1 & TPXO) {
    master->tpxotide = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tpxotide");
    strcpy(master->trinfo_2d[tn].long_name, "TPXO tide");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers1 & TPXOV) {
    master->tpxovelu = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tpxou");
    strcpy(master->trinfo_2d[tn].long_name, "TPXO x velocity");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->tpxovelv = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tpxov");
    strcpy(master->trinfo_2d[tn].long_name, "TPXO v velocity");
    strcpy(master->trinfo_2d[tn].units, "ms-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers1 & TPXOT) {
    master->tpxotranu = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tpxoU");
    strcpy(master->trinfo_2d[tn].long_name, "TPXO x transport");
    strcpy(master->trinfo_2d[tn].units, "m2s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->tpxotranv = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "tpxoV");
    strcpy(master->trinfo_2d[tn].long_name, "TPXO v transport");
    strcpy(master->trinfo_2d[tn].units, "m2s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->numbers1 & TRAN2D) {
    master->uat = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "UAT");
    strcpy(master->trinfo_2d[tn].long_name, "Eastward 2D transport");
    strcpy(master->trinfo_2d[tn].units, "m2s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->vat = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "VAT");
    strcpy(master->trinfo_2d[tn].long_name, "Northward 2D transport");
    strcpy(master->trinfo_2d[tn].units, "m2s-1");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->decf & DEC_ETA) {
    master->decv1 = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "decorr_e1");
    strcpy(master->trinfo_2d[tn].long_name, "Decorrelation length scale");
    strcpy(master->trinfo_2d[tn].units, params->decs);
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  if (params->etarlx & ETA_ADPT) {
    master->eta_tc = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "eta_tc");
    strcpy(master->trinfo_2d[tn].long_name, "Relaxation eta time constant");
    strcpy(master->trinfo_2d[tn].units, "days");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = 0;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
    master->eta_inc = master->tr_wcS[tn];
    strcpy(master->trinfo_2d[tn].name, "eta_inc");
    strcpy(master->trinfo_2d[tn].long_name, "Relaxation eta increment");
    strcpy(master->trinfo_2d[tn].units, "m");
    master->trinfo_2d[tn].type = INTER|HYDRO|DIAGNOSTIC;
    tr_dataset(buf, &master->trinfo_2d[tn], 0.0);
    master->trinfo_2d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_2d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_2d[tn].n = tn;
    tn++;
  }
  /* Sediment error maps */
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
    /* Set up sediment tracers if required */
    tn = sediment_autotracer_2d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, master->trinfo_2d, master->ntrS, tn);
  }
#endif
  /* Ecology error maps */
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
    /* Set up sediment tracers if required */
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

/* END init_tracer_sed()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the 3D tracers in the master                */
/*-------------------------------------------------------------------*/
void init_tracer_3d(parameters_t *params, /* Input parameters data   */
                    master_t *master      /* Master data             */
  )
{
  int tn;
  geometry_t *geom = master->sgrid;
  char buf[MAXSTRLEN];
  if (geom == NULL) geom = master->geom;

  /* Initialise pointers */
  master->sal = master->temp = NULL;
  master->tke = master->diss = master->L = master->omega = NULL;
  master->Q2 = master->Q2L = master->Kq = master->sdc = NULL;
  master->u1m = master->u2m = master->wm = master->Kzm = NULL;
  master->fluxe1 = master->fluxe2 = master->tempm = master->saltm = master->tram = NULL;
  master->fluxw = master->fluxkz = NULL;
  master->brunt = master->int_wave = master->rich_gr = master->rich_fl = NULL;
  master->froude = master->reynolds = master->rossby_in = master->sound = NULL;
  master->shear_v = master->b_prod = master->s_prod = master->speed_3d = NULL;
  master->otemp = master->osalt = master->perc = master->fltr = master->agetr = NULL;
  master->rtemp = master->rsalt = master->schan = master->sigma_t = NULL;
  master->ptconc = master->energy = master->vcorr = master->acorr = NULL;
  master->dum1 = master->dum2 = master->dum3 = master->kenergy = master->riversalt = NULL;
  master->regionid = master->regres = master->Vi = master->reefe1 = master->reefe2 = NULL;
  master->temp_tc = master->salt_tc = master->unit = master->mono = NULL;
  master->wave_stke1 = master->wave_stke2 = NULL;
  if (params->swr_type & SWR_3D) master->swr_attn = NULL;
  master->dhw = master->dhd = master->dhwc = master->glider = master->nprof = master->u1vhc = NULL;

  /*-----------------------------------------------------------------*/
  /* Assign the water column tracer pointers to any tracers defined */
  /* in the input parameter file.  */
  for (tn = 0; tn < master->ntr; tn++) {
    if (strcmp("salt", master->trinfo_3d[tn].name) == 0) {
      master->sal = master->tr_wc[tn];
      master->sno = tn;
    }
    else if (strcmp("temp", master->trinfo_3d[tn].name) == 0) {
      master->temp = master->tr_wc[tn];
      master->tno = tn;
    }
    else if (strcmp("tke", master->trinfo_3d[tn].name) == 0) {
      master->tke = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("diss", master->trinfo_3d[tn].name) == 0) {
      master->diss = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("omega", master->trinfo_3d[tn].name) == 0) {
      master->omega = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("tki", master->trinfo_3d[tn].name) == 0) {
      master->Q2 = master->tr_wc[tn];
      master->trinfo_3d[tn].diffuse = 0;
    } else if (strcmp("tki_l", master->trinfo_3d[tn].name) == 0) {
      master->trinfo_3d[tn].diffuse = 0;
      master->Q2L = master->tr_wc[tn];
    } else if (strcmp("Kq", master->trinfo_3d[tn].name) == 0) {
      master->Kq = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    } else if (strcmp("lscale", master->trinfo_3d[tn].name) == 0) {
      master->L = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    } else if (strcmp("u1mean", master->trinfo_3d[tn].name) == 0)
      master->u1m = master->tr_wc[tn];
    else if (strcmp("u2mean", master->trinfo_3d[tn].name) == 0)
      master->u2m = master->tr_wc[tn];
    else if (strcmp("wmean", master->trinfo_3d[tn].name) == 0)
      master->wm = master->tr_wc[tn];
    else if (strcmp("temp_mean", master->trinfo_3d[tn].name) == 0)
      master->tempm = master->tr_wc[tn];
    else if (strcmp("salt_mean", master->trinfo_3d[tn].name) == 0)
      master->saltm = master->tr_wc[tn];
    else if (strcmp("tracer_mean", master->trinfo_3d[tn].name) == 0)
      master->tram = master->tr_wc[tn];
    else if (strcmp("Kzmean", master->trinfo_3d[tn].name) == 0)
      master->Kzm = master->tr_wc[tn];
    else if (strcmp("flux_e1", master->trinfo_3d[tn].name) == 0)
      master->fluxe1 = master->tr_wc[tn];
    else if (strcmp("flux_e2", master->trinfo_3d[tn].name) == 0)
      master->fluxe2 = master->tr_wc[tn];
    else if (strcmp("flux_w", master->trinfo_3d[tn].name) == 0)
      master->fluxw = master->tr_wc[tn];
    else if (strcmp("flux_kz", master->trinfo_3d[tn].name) == 0)
      master->fluxkz = master->tr_wc[tn];
    else if (strcmp("flush", master->trinfo_3d[tn].name) == 0)
      master->fltr = master->tr_wc[tn];
    else if (strcmp("age", master->trinfo_3d[tn].name) == 0)
      master->agetr = master->tr_wc[tn];
    else if (strcmp("smagorinsky", master->trinfo_3d[tn].name) == 0) {
      master->sdc = master->tr_wc[tn];
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].diagn = 0;
    }
    else if (strcmp("brunt_vaisala", master->trinfo_3d[tn].name) == 0)
      master->brunt = master->tr_wc[tn];
    else if (strcmp("int_wave_speed", master->trinfo_3d[tn].name) == 0)
      master->int_wave = master->tr_wc[tn];
    else if (strcmp("richardson_gr", master->trinfo_3d[tn].name) == 0)
      master->rich_gr = master->tr_wc[tn];
    else if (strcmp("richardson_fl", master->trinfo_3d[tn].name) == 0)
      master->rich_fl = master->tr_wc[tn];
    else if (strcmp("reynolds", master->trinfo_3d[tn].name) == 0)
      master->reynolds = master->tr_wc[tn];
    else if (strcmp("froude", master->trinfo_3d[tn].name) == 0)
      master->froude = master->tr_wc[tn];
    else if (strcmp("sigma_t", master->trinfo_3d[tn].name) == 0)
      master->sigma_t = master->tr_wc[tn];
    else if (strcmp("energy", master->trinfo_3d[tn].name) == 0)
      master->energy = master->tr_wc[tn];
    else if (strcmp("kenergy", master->trinfo_3d[tn].name) == 0)
      master->kenergy = master->tr_wc[tn];
    else if (strcmp("ptconc", master->trinfo_3d[tn].name) == 0)
      master->ptconc = master->tr_wc[tn];
    else if (strcmp("sound", master->trinfo_3d[tn].name) == 0)
      master->sound = master->tr_wc[tn];
    else if (strcmp("sound_channel", master->trinfo_3d[tn].name) == 0)
      master->schan = master->tr_wc[tn];
    else if (strcmp("shear_vert", master->trinfo_3d[tn].name) == 0)
      master->shear_v = master->tr_wc[tn];
    else if (strcmp("buoy_prod", master->trinfo_3d[tn].name) == 0)
      master->b_prod = master->tr_wc[tn];
    else if (strcmp("shear_prod", master->trinfo_3d[tn].name) == 0)
      master->s_prod = master->tr_wc[tn];
    else if (strcmp("rossby_in", master->trinfo_3d[tn].name) == 0)
      master->rossby_in = master->tr_wc[tn];
    else if (strcmp("current_speed_3d", master->trinfo_3d[tn].name) == 0)
      master->speed_3d = master->tr_wc[tn];
    else if (strcmp("otemp", master->trinfo_3d[tn].name) == 0)
      master->otemp = master->tr_wc[tn];
    else if (strcmp("osalt", master->trinfo_3d[tn].name) == 0)
      master->osalt = master->tr_wc[tn];
    else if (strcmp("rtemp", master->trinfo_3d[tn].name) == 0)
      master->rtemp = master->tr_wc[tn];
    else if (strcmp("rsalt", master->trinfo_3d[tn].name) == 0)
      master->rsalt = master->tr_wc[tn];
    else if (strcmp("temp_tc", master->trinfo_3d[tn].name) == 0)
      master->temp_tc = master->tr_wc[tn];
    else if (strcmp("salt_tc", master->trinfo_3d[tn].name) == 0)
      master->salt_tc = master->tr_wc[tn];
    else if (strcmp("tracer1", master->trinfo_3d[tn].name) == 0)
      master->dum1 = master->tr_wc[tn];
    else if (strcmp("tracer2", master->trinfo_3d[tn].name) == 0)
      master->dum2 = master->tr_wc[tn];
    else if (strcmp("tracer3", master->trinfo_3d[tn].name) == 0)
      master->dum3 = master->tr_wc[tn];
    else if (strcmp("reef_fraction_e1", master->trinfo_3d[tn].name) == 0)
      master->reefe1 = master->tr_wc[tn];
    else if (strcmp("reef_fraction_e2", master->trinfo_3d[tn].name) == 0)
      master->reefe2 = master->tr_wc[tn];
    else if (strcmp("Vcorr", master->trinfo_3d[tn].name) == 0)
      master->vcorr = master->tr_wc[tn];
    else if (strcmp("Acorr", master->trinfo_3d[tn].name) == 0)
      master->acorr = master->tr_wc[tn];
    else if (strcmp("Vi", master->trinfo_3d[tn].name) == 0)
      master->Vi = master->tr_wc[tn];
    else if (strcmp("unit", master->trinfo_3d[tn].name) == 0)
      master->unit = master->tr_wc[tn];
    else if (strcmp("regionid", master->trinfo_3d[tn].name) == 0)
      master->regionid = master->tr_wc[tn];
    else if (strcmp("residence", master->trinfo_3d[tn].name) == 0)
      master->regres = master->tr_wc[tn];
    else if (strcmp("decorr_e1", master->trinfo_3d[tn].name) == 0)
      master->decv1 = master->tr_wc[tn];
    else if (strcmp("wave_stke1", master->trinfo_3d[tn].name) == 0)
      master->wave_stke1 = master->tr_wc[tn];
    else if (strcmp("wave_stke2", master->trinfo_3d[tn].name) == 0)
      master->wave_stke2 = master->tr_wc[tn];
    else if (strcmp("dhw", master->trinfo_3d[tn].name) == 0)
      master->dhw = master->tr_wc[tn];
    else if (strcmp("dhd", master->trinfo_3d[tn].name) == 0)
      master->dhd = master->tr_wc[tn];
    else if (strcmp("dhwc", master->trinfo_3d[tn].name) == 0)
      master->dhwc = master->tr_wc[tn];
    else if (strcmp("glider", master->trinfo_3d[tn].name) == 0)
      master->glider = master->tr_wc[tn];
    else if (strcmp("u1vhc", master->trinfo_3d[tn].name) == 0)
      master->u1vhc = master->tr_wc[tn];
    else if (strcmp("nprof", master->trinfo_3d[tn].name) == 0)
      master->nprof = master->tr_wc[tn];
    else if (strcmp("mono", master->trinfo_3d[tn].name) == 0)
      master->mono = master->tr_wc[tn];
    else if (params->swr_type & SWR_3D && strcmp("swr_attenuation", master->trinfo_3d[tn].name) == 0)
      master->swr_attn = master->tr_wc[tn];
  }

  sprintf(buf, "percentile_%s", params->trperc);
  for (tn = 0; tn < master->ntr; tn++) {
    if (strcmp(buf, master->trinfo_3d[tn].name) == 0)
      master->perc = master->tr_wc[tn];
  }

  /* Assign the automtically generated 3D pointers and set up the */
  /* tracer_info structure. If any of these tracers have been */
  /* explicitly defined in the input parameter file, then continue.  */
  tn = 0;
  /* Temperature and salinity */
  if ((tracer_find_index("salt", master->ntr, master->trinfo_3d)) < 0) {
    master->sal = master->tr_wc[tn];
    master->sno = tn;
    strcpy(master->trinfo_3d[tn].name, "salt");
    strcpy(master->trinfo_3d[tn].long_name, "Salinity");
    strcpy(master->trinfo_3d[tn].units, "PSU");
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 40;
    master->trinfo_3d[tn].fill_value_wc = 35.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].dissol = 1;
    master->trinfo_3d[tn].advect = 1;
    master->trinfo_3d[tn].diffuse = 1;
    strcpy(master->trinfo_3d[tn].decay, "0.0");
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].m = tn;
    tn++;
  }
  if ((tracer_find_index("temp", master->ntr, master->trinfo_3d)) < 0) {
    master->temp = master->tr_wc[tn];
    master->tno = tn;
    strcpy(master->trinfo_3d[tn].name, "temp");
    strcpy(master->trinfo_3d[tn].long_name, "Temperature");
    strcpy(master->trinfo_3d[tn].units, "degrees C");
    master->trinfo_3d[tn].valid_range_wc[0] = -4;
    master->trinfo_3d[tn].valid_range_wc[1] = 40;
    master->trinfo_3d[tn].fill_value_wc = 20.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].dissol = 1;
    master->trinfo_3d[tn].advect = 1;
    master->trinfo_3d[tn].diffuse = 1;
    strcpy(master->trinfo_3d[tn].decay, "0.0");
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].m = tn;
    tn++;
  }
  /* 3D mean velocity */
  if (params->means & VEL3D) {
    if (tracer_find_index("u1mean", master->ntr, master->trinfo_3d) == -1) {
      master->u1m = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "u1mean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean 3D east velocity");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC|E1VAR;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -100;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      strcpy(master->trinfo_3d[tn].vector_name, "Mean 3D current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1mean u2mean");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->u1m, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
    if (tracer_find_index("u2mean", master->ntr, master->trinfo_3d) == -1) {
      master->u2m = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "u2mean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean 3D north velocity");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC|E2VAR;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -100;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      strcpy(master->trinfo_3d[tn].vector_name, "Mean 3D current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1mean u2mean");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->u2m, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
    if (tracer_find_index("wmean", master->ntr, master->trinfo_3d) == -1) {
      master->wm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "wmean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean w velocity");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -100;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->wm, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
  }
    /* Temperature and salinity */
  if (params->means & TS) {
    if (tracer_find_index("temp_mean", master->ntr, master->trinfo_3d) == -1) {
      master->tempm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "temp_mean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean temperature");
      strcpy(master->trinfo_3d[tn].units, "deg C");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->tempm, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
    if (tracer_find_index("salt_mean", master->ntr, master->trinfo_3d) == -1) {
      master->saltm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "salt_mean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean salinity");
      strcpy(master->trinfo_3d[tn].units, "psu");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->saltm, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
  }
  /* Tracer mean */
  if (params->means & MTRA3D) {
    if (tracer_find_index("tracer_mean", master->ntr, master->trinfo_3d) == -1) {
      master->tram = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tracer_mean");
      sprintf(buf, "Mean %s", params->means_tra);
      strcpy(master->trinfo_3d[tn].long_name, buf);
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->tram, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
  }
  /* 3D mean volume flux */
  /*
  if ((params->means & VOLFLUX) || (params->tmode & SP_FFSL)) {
    if (tracer_find_index("u1vmean", master->ntr, master->trinfo_3d) == -1) {
      master->u1vm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "u1vmean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean u1 volume flux");
      strcpy(master->trinfo_3d[tn].units, "m3s-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC|E1VAR;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e30;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e30;
      strcpy(master->trinfo_3d[tn].vector_name, "Mean volume flux");
      strcpy(master->trinfo_3d[tn].vector_components, "u1vmean u2vmean");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->u1vm, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
    if (tracer_find_index("u2vmean", master->ntr, master->trinfo_3d) == -1) {
      master->u2vm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "u2vmean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean u2 volume flux");
      strcpy(master->trinfo_3d[tn].units, "m3s-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC|E2VAR;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e30;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e30;
      strcpy(master->trinfo_3d[tn].vector_name, "Mean volume flux");
      strcpy(master->trinfo_3d[tn].vector_components, "u1vmean u2vmean");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->u2vm, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
  }
  */
  /* Mean Kz */
  if (params->means & KZ_M) {
    if (tracer_find_index("Kzmean", master->ntr, master->trinfo_3d) == -1) {
      master->Kzm = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "Kzmean");
      strcpy(master->trinfo_3d[tn].long_name, "Mean Vertical Diffusivity");
      strcpy(master->trinfo_3d[tn].units, "m2s-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Tracer fluxes */
  if (strcmp(params->trflux, "NONE") != 0) {
    if (tracer_find_index("flux_e1", master->ntr, master->trinfo_3d) == -1) {
      master->fluxe1 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flux_e1");
      sprintf(master->trinfo_3d[tn].long_name, "Advective flux through face %d", params->trfd1);
      strcpy(master->trinfo_3d[tn].units, "kgs-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("flux_e2", master->ntr, master->trinfo_3d) == -1) {
      master->fluxe2 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flux_e2");
      sprintf(master->trinfo_3d[tn].long_name, "Advective flux through face %d", params->trfd2);
      strcpy(master->trinfo_3d[tn].units, "kgs-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("flux_w", master->ntr, master->trinfo_3d) == -1) {
      master->fluxw = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flux_w");
      strcpy(master->trinfo_3d[tn].long_name, "Vertical advective flux");
      strcpy(master->trinfo_3d[tn].units, "kgs-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("flux_kz", master->ntr, master->trinfo_3d) == -1) {
      master->fluxkz = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flux_kz");
      strcpy(master->trinfo_3d[tn].long_name, "Vertical diffusive flux");
      strcpy(master->trinfo_3d[tn].units, "kgs-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Regions */
  if (strlen(params->regions)) {
    if (tracer_find_index("regionid", master->ntr, master->trinfo_3d) == -1) {
      master->regionid = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "regionid");
      strcpy(master->trinfo_3d[tn].long_name, "Region identifier");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("residence", master->ntr, master->trinfo_3d) == -1) {
      master->regres = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "residence");
      strcpy(master->trinfo_3d[tn].long_name, "Residence time");
      strcpy(master->trinfo_3d[tn].units, "days");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
 /* Flushing tracer */
  if (params->trflsh) {
    if (tracer_find_index("flush", master->ntr, master->trinfo_3d) == -1) {
      master->fltr = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flush");
      strcpy(master->trinfo_3d[tn].long_name, "Flushing tracer");
      strcpy(master->trinfo_3d[tn].units, "mgL-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
 /* Age tracer */
  if (strlen(params->trage)) {
    if (tracer_find_index("age", master->ntr, master->trinfo_3d) == -1) {
      master->agetr = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "age");
      strcpy(master->trinfo_3d[tn].long_name, "Age tracer");
      strcpy(master->trinfo_3d[tn].units, "days");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Tracer percentiles */
  if (strcmp(params->trperc, "NONE") != 0) {
    sprintf(buf, "percentile_%s", params->trperc);
    if (tracer_find_index(buf, master->ntr, master->trinfo_3d) == -1) {
      master->perc = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, buf);
      sprintf(buf, "Percentile for %s", params->trperc);
      strcpy(master->trinfo_3d[tn].long_name, buf);
      strcpy(master->trinfo_3d[tn].units, "%");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* k-e closure */
  if (strcmp(params->mixsc, "k-e") == 0) {
    if (tracer_find_index("tke", master->ntr, master->trinfo_3d) == -1) {
      master->tke = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tke");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulent Kinetic Energy");
      strcpy(master->trinfo_3d[tn].units, "m2s-2");
      master->trinfo_3d[tn].fill_value_wc = 7.6e-6;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("diss", master->ntr, master->trinfo_3d) == -1) {
      master->diss = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "diss");
      strcpy(master->trinfo_3d[tn].long_name, "Dissipation");
      strcpy(master->trinfo_3d[tn].units, "m2s-3");
      master->trinfo_3d[tn].fill_value_wc = 5e-10;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* k-w closure */
  if (strcmp(params->mixsc, "k-w") == 0 || strcmp(params->mixsc, "W88") == 0) {
    if (tracer_find_index("tke", master->ntr, master->trinfo_3d) == -1) {
      master->tke = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tke");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulent Kinetic Energy");
      strcpy(master->trinfo_3d[tn].units, "m2s-2");
      master->trinfo_3d[tn].fill_value_wc = 7.6e-6;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("omega", master->ntr, master->trinfo_3d) == -1) {
      master->omega = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "omega");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulence frequency");
      strcpy(master->trinfo_3d[tn].units, "s-1");
      master->trinfo_3d[tn].fill_value_wc = 1.0e-12;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Mellor-Yamada 2.5 closure */
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
      strcmp(params->mixsc, "harcourt") == 0) {
    if (tracer_find_index("tki", master->ntr, master->trinfo_3d) == -1) {
      master->Q2 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tki");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulent Kinetic Intensity");
      strcpy(master->trinfo_3d[tn].units, "m2s-2");
      master->trinfo_3d[tn].fill_value_wc = 2.0e-8;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("tki_l", master->ntr, master->trinfo_3d) == -1) {
      master->Q2L = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tki_l");
      strcpy(master->trinfo_3d[tn].long_name,
             "Turbulent Kinetic Intensity Length Scale");
      strcpy(master->trinfo_3d[tn].units, "m2s-1");
      master->trinfo_3d[tn].fill_value_wc = 3.4e-9;
      master->trinfo_3d[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("lscale", master->ntr, master->trinfo_3d) == -1) {
      master->L = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "lscale");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulence length scale");
      strcpy(master->trinfo_3d[tn].units, "metre");
      master->trinfo_3d[tn].fill_value_wc = 0.17;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e4;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("Kq", master->ntr, master->trinfo_3d) == -1) {
      master->Kq = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "Kq");
      strcpy(master->trinfo_3d[tn].long_name, "Turbulence mixing");
      strcpy(master->trinfo_3d[tn].units, "m2s-1");
      master->trinfo_3d[tn].fill_value_wc = 1e-5;
      master->trinfo_3d[tn].type = WATER|HYDRO|PARAMETER;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Smagorinsky diffusion */
  if (params->smagorinsky > 0.0) {
    if (tracer_find_index("smagorinsky", master->ntr, master->trinfo_3d) ==
        -1) {
      master->sdc = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "smagorinsky");
      strcpy(master->trinfo_3d[tn].long_name, "Smagorinsky diffusion");
      strcpy(master->trinfo_3d[tn].units, "m2s-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|PARAMETER;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      memset(master->sdc, 0, geom->sgsiz * sizeof(double));
      tn++;
    }
  }
  /* Layer thickness */
  if (params->show_layers) {
    master->layth = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "layer_thick");
    strcpy(master->trinfo_3d[tn].long_name, "Layer thickness");
    strcpy(master->trinfo_3d[tn].units, "metre");
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e4;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  /* Forcing data saved to file */
  if (params->save_force & OTEMP &&
      tracer_find_index("otemp", master->ntr, master->trinfo_3d) == -1) {
    master->otemp = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "otemp");
    strcpy(master->trinfo_3d[tn].long_name, "OFAM temperature");
    strcpy(master->trinfo_3d[tn].units, "degrees C");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -4;
    master->trinfo_3d[tn].valid_range_wc[1] = 40;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    if (strlen(params->odata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(otemp=temp)", params->tdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(master->trinfo_3d[tn].reset_file, buf);
      strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->save_force & OSALT &&
      tracer_find_index("osalt", master->ntr, master->trinfo_3d) == -1) {
    master->osalt = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "osalt");
    strcpy(master->trinfo_3d[tn].long_name, "OFAM salinity");
    strcpy(master->trinfo_3d[tn].units, "PSU");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 40;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    if (strlen(params->odata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(osalt=salt)", params->sdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(master->trinfo_3d[tn].reset_file, buf);
      strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
    }
    tn++;
  }
  /* Relaxation for temp and salt */
  if (params->rtemp &&
      tracer_find_index("rtemp", master->ntr, master->trinfo_3d) == -1) {
    master->rtemp = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "rtemp");
    strcpy(master->trinfo_3d[tn].long_name, "Relaxation temperature");
    strcpy(master->trinfo_3d[tn].units, "degrees C");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 50;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    if (strlen(params->tdata)) {
      char buf[MAXSTRLEN];
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
    tn++;
  }
  if (params->rsalt &&
      tracer_find_index("rsalt", master->ntr, master->trinfo_3d) == -1) {
    master->rsalt = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "rsalt");
    strcpy(master->trinfo_3d[tn].long_name, "Relaxation salinity");
    strcpy(master->trinfo_3d[tn].units, "PSU");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 40;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
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
      strcpy(master->trinfo_3d[tn].reset_file, buf);
      strcpy(master->trinfo_3d[tn].reset_dt, buf1);
    }
    tn++;
  }
  if (params->rtemp & (RLX_ADPT|RLX_REG|RLX_OBC) &&
      tracer_find_index("temp_tc", master->ntr, master->trinfo_3d) == -1) {
    master->temp_tc = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "temp_tc");
    strcpy(master->trinfo_3d[tn].long_name, "Relaxation temperature time constant");
    strcpy(master->trinfo_3d[tn].units, "days");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->rsalt & (RLX_ADPT|RLX_REG|RLX_OBC) &&
      tracer_find_index("salt_tc", master->ntr, master->trinfo_3d) == -1) {
    master->salt_tc = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "salt_tc");
    strcpy(master->trinfo_3d[tn].long_name, "Relaxation salinity time constant");
    strcpy(master->trinfo_3d[tn].units, "days");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->save_force & OVELU &&
      tracer_find_index("ovelu", master->ntr, master->trinfo_3d) == -1) {
    strcpy(master->trinfo_3d[tn].name, "ovelu");
    strcpy(master->trinfo_3d[tn].long_name, "OFAM east velocity");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -100;
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    strcpy(master->trinfo_3d[tn].vector_name, "Global 3D current");
    strcpy(master->trinfo_3d[tn].vector_components, "ovelu ovelv");
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    if (strlen(params->vdata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(ovelu=u)", params->vdata);
      else
	sprintf(buf, "%s", params->vdata);
      strcpy(master->trinfo_3d[tn].reset_file, buf);
      strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->save_force & OVELV &&
      tracer_find_index("ovelv", master->ntr, master->trinfo_3d) == -1) {
    strcpy(master->trinfo_3d[tn].name, "ovelv");
    strcpy(master->trinfo_3d[tn].long_name, "OFAM north velocity");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -100;
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    strcpy(master->trinfo_3d[tn].vector_name, "Global 3D current");
    strcpy(master->trinfo_3d[tn].vector_components, "ovelu ovelv");
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    if (strlen(params->vdata)) {
      char buf[MAXSTRLEN];
      if (params->save_force & ROAM)
	sprintf(buf, "%s(ovelv=v)", params->vdata);
      else
	sprintf(buf, "%s", params->vdata);
      strcpy(master->trinfo_3d[tn].reset_file, buf);
      strcpy(master->trinfo_3d[tn].reset_dt, "1 second");
    }
    tn++;
  }
  /* Momentum tendencies */
  if (params->tendf) {
    if (tracer_find_index("mom_balance", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "mom_balance");
      strcpy(master->trinfo_3d[tn].long_name, "Momentum balance maximum");
      strcpy(master->trinfo_3d[tn].units, "%");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_adv", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_adv");
      strcpy(master->trinfo_3d[tn].long_name, "East advective tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Advective current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_adv u2_adv");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_hdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_hdif");
      strcpy(master->trinfo_3d[tn].long_name,
             "East horizontal diffusion tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Hdiffusion current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_hdif u2_hdif");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_vdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_vdif");
      strcpy(master->trinfo_3d[tn].long_name,
             "East vertical diffusion tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Vdiffusion current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_vdif u2_vdif");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_btp", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_btp");
      strcpy(master->trinfo_3d[tn].long_name,
             "East barotropic pressure gradient tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Barotropic current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_btp u2_btp");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_bcp", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_bcp");
      strcpy(master->trinfo_3d[tn].long_name,
             "East baroclinic pressure tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Baroclinic current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_bcp u2_bcp");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u1_cor", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u1_cor");
      strcpy(master->trinfo_3d[tn].long_name, "East Coriolis tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Coriolis surrent");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_cor u2_cor");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (params->waves & STOKES_DRIFT) {
      if (tracer_find_index("u1_sto", master->ntr, master->trinfo_3d) == -1) {
	strcpy(master->trinfo_3d[tn].name, "u1_sto");
	strcpy(master->trinfo_3d[tn].long_name, "East Stokes tendency");
	strcpy(master->trinfo_3d[tn].units, "ms-1");
	master->trinfo_3d[tn].fill_value_wc = 0.0;
	master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
	master->trinfo_3d[tn].diagn = 0;
	master->trinfo_3d[tn].advect = 0;
	master->trinfo_3d[tn].diffuse = 0;
	master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
	master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
	strcpy(master->trinfo_3d[tn].vector_name, "Stokes current");
	strcpy(master->trinfo_3d[tn].vector_components, "u1_sto u2_sto");
	master->trinfo_3d[tn].m = -1;
	master->trinfo_3d[tn].n = tn;
	tn++;
      }
    }
    if (tracer_find_index("u2_adv", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_adv");
      strcpy(master->trinfo_3d[tn].long_name, "North advective tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Advective current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_adv u2_adv");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u2_hdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_hdif");
      strcpy(master->trinfo_3d[tn].long_name,
             "North horizontal diffusion tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Hdiffusion current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_hdif u2_hdif");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u2_vdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_vdif");
      strcpy(master->trinfo_3d[tn].long_name,
             "North vertical diffusion tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Visffusion current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_vdif u2_vdif");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u2_btp", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_btp");
      strcpy(master->trinfo_3d[tn].long_name,
             "North barotropic pressure gradient tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Barotropic currents");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_btp u2_btp");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u2_bcp", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_bcp");
      strcpy(master->trinfo_3d[tn].long_name,
             "North baroclinic pressure tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Baroclinic current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_bcp u2_bcp");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("u2_cor", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "u2_cor");
      strcpy(master->trinfo_3d[tn].long_name, "North Coriolis tendency");
      strcpy(master->trinfo_3d[tn].units, "ms-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      strcpy(master->trinfo_3d[tn].vector_name, "Coriolis current");
      strcpy(master->trinfo_3d[tn].vector_components, "u1_cor u2_cor");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (params->waves & STOKES_DRIFT) {
      if (tracer_find_index("u2_sto", master->ntr, master->trinfo_3d) == -1) {
	strcpy(master->trinfo_3d[tn].name, "u2_sto");
	strcpy(master->trinfo_3d[tn].long_name, "North Stokes tendency");
	strcpy(master->trinfo_3d[tn].units, "ms-1");
	master->trinfo_3d[tn].fill_value_wc = 0.0;
	master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
	master->trinfo_3d[tn].diagn = 0;
	master->trinfo_3d[tn].advect = 0;
	master->trinfo_3d[tn].diffuse = 0;
	master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
	master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
	strcpy(master->trinfo_3d[tn].vector_name, "Stokes current");
	strcpy(master->trinfo_3d[tn].vector_components, "u1_sto u2_sto");
	master->trinfo_3d[tn].m = -1;
	master->trinfo_3d[tn].n = tn;
	tn++;
      }
    }
  }
  /* Tracer tendencies */
  if (strlen(params->trtend)) {
    if (tracer_find_index("tra_adv", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "tra_adv");
      strcpy(master->trinfo_3d[tn].long_name, "Tracer advective tendency");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    /*
    if (tracer_find_index("tra_hdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "tra_hdif");
      strcpy(master->trinfo_3d[tn].long_name, "Tracer horizontal diffusive tendency");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    */
    if (tracer_find_index("tra_vdif", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "tra_vdif");
      strcpy(master->trinfo_3d[tn].long_name, "Tracer vertical diffusive tendency");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("tra_ncon", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "tra_ncon");
      strcpy(master->trinfo_3d[tn].long_name, "Tracer non-conservative tendency");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e35;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e35;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  /* Waves */
  if (params->waves & SPECTRAL) {
    master->wave_stke1 = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "wave_stke1");
    strcpy(master->trinfo_3d[tn].long_name, "Stokes sub-surface velocity along e1");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|WAVE|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -100;
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
    master->wave_stke2 = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "wave_stke2");
    strcpy(master->trinfo_3d[tn].long_name, "Stokes sub-surface velocity along e2");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|WAVE|FORCING;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -100;
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  /* Diagnostic numbers */
  if (params->numbers & BRUNT &&
      tracer_find_index("brunt_vaisala", master->ntr, master->trinfo_3d) == -1) {
    master->brunt = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "brunt_vaisala");
    strcpy(master->trinfo_3d[tn].long_name, "Brunt Vaisala Frequency");
    strcpy(master->trinfo_3d[tn].units, "s-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & INT_WAVE &&
      tracer_find_index("int_wave_speed", master->ntr, master->trinfo_3d) == -1) {
    master->int_wave = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "int_wave_speed");
    strcpy(master->trinfo_3d[tn].long_name, "Internal wave speed");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & RICHARD_GR &&
      tracer_find_index("richardson_gr", master->ntr, master->trinfo_3d) == -1) {
    master->rich_gr = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "richardson_gr");
    strcpy(master->trinfo_3d[tn].long_name, "Gradient Richardson number");
    strcpy(master->trinfo_3d[tn].units, "");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e308;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e308;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & RICHARD_FL &&
      tracer_find_index("richardson_fl", master->ntr, master->trinfo_3d) == -1) {
    master->rich_fl = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "richardson_fl");
    strcpy(master->trinfo_3d[tn].long_name, "Flux Richardson number");
    strcpy(master->trinfo_3d[tn].units, "");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e308;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e308;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & REYNOLDS &&
      tracer_find_index("reynolds", master->ntr, master->trinfo_3d) == -1) {
    master->reynolds = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "reynolds");
    strcpy(master->trinfo_3d[tn].long_name, "Reynolds number");
    strcpy(master->trinfo_3d[tn].units, "");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & FROUDE &&
      tracer_find_index("froude", master->ntr, master->trinfo_3d) == -1) {
    master->froude = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "froude");
    strcpy(master->trinfo_3d[tn].long_name, "Froude number");
    strcpy(master->trinfo_3d[tn].units, "");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SIGMA_T &&
      tracer_find_index("sigma_t", master->ntr, master->trinfo_3d) == -1) {
    master->sigma_t = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "sigma_t");
    strcpy(master->trinfo_3d[tn].long_name, "Sigma_t");
    strcpy(master->trinfo_3d[tn].units, "kgm-3");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -100;
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & ENERGY &&
      tracer_find_index("energy", master->ntr, master->trinfo_3d) == -1) {
    master->energy = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "energy");
    strcpy(master->trinfo_3d[tn].long_name, "Mechanical energy");
    strcpy(master->trinfo_3d[tn].units, "Jm-3");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & KINETIC &&
      tracer_find_index("kenergy", master->ntr, master->trinfo_3d) == -1) {
    master->kenergy = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "kenergy");
    strcpy(master->trinfo_3d[tn].long_name, "Kinetic energy");
    strcpy(master->trinfo_3d[tn].units, "Jm-3");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->do_pt && 
      tracer_find_index("ptconc", master->ntr, master->trinfo_3d) == -1) {
    master->ptconc = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "ptconc");
    strcpy(master->trinfo_3d[tn].long_name, "Particle concentration");
    strcpy(master->trinfo_3d[tn].units, "kgm-3");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 0;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SOUND &&
      tracer_find_index("sound", master->ntr, master->trinfo_3d) == -1) {
    master->sound = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "sound");
    strcpy(master->trinfo_3d[tn].long_name, "Speed of sound");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 1000;
    master->trinfo_3d[tn].valid_range_wc[1] = 2000;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SOUND &&
      tracer_find_index("sound_channel", master->ntr, master->trinfo_3d) == -1) {
    master->sound = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "sound_channel");
    strcpy(master->trinfo_3d[tn].long_name, "Sound channel depth");
    strcpy(master->trinfo_3d[tn].units, "m");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
/*UR    master->trinfo_3d[tn].positive = down; */
    master->trinfo_3d[tn].valid_range_wc[0] = -1e4;
/*UR made 100 as not-realistic but theoretically possible limit */
    master->trinfo_3d[tn].valid_range_wc[1] = 100;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & ROSSBY_IN &&
      tracer_find_index("rossby_internal", master->ntr, master->trinfo_3d) == -1) {
    master->rossby_in = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "rossby_internal");
    strcpy(master->trinfo_3d[tn].long_name, "Internal Rossby radius");
    strcpy(master->trinfo_3d[tn].units, "metre");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = 0;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SPEED_3D &&
      tracer_find_index("current_speed_3d", master->ntr, master->trinfo_3d) == -1) {
    master->speed_3d = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "current_speed_3d");
    strcpy(master->trinfo_3d[tn].long_name, "Current Speed 3D");
    strcpy(master->trinfo_3d[tn].units, "ms-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SHEAR_V &&
      tracer_find_index("shear_vert", master->ntr, master->trinfo_3d) == -1) {
    master->shear_v = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "shear_vert");
    strcpy(master->trinfo_3d[tn].long_name, "Vertical velocity shear");
    strcpy(master->trinfo_3d[tn].units, "s-1");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & BUOY_PROD &&
      tracer_find_index("buoy_prod", master->ntr, master->trinfo_3d) == -1) {
    master->b_prod = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "buoy_prod");
    strcpy(master->trinfo_3d[tn].long_name, "Buoyancy production");
    strcpy(master->trinfo_3d[tn].units, "m2s-2");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (params->numbers & SHEAR_PROD &&
      tracer_find_index("shear_prod", master->ntr, master->trinfo_3d) == -1) {
    master->s_prod = master->tr_wc[tn];
    strcpy(master->trinfo_3d[tn].name, "shear_prod");
    strcpy(master->trinfo_3d[tn].long_name, "Shear production");
    strcpy(master->trinfo_3d[tn].units, "m2s-2");
    master->trinfo_3d[tn].fill_value_wc = 0.0;
    master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
    master->trinfo_3d[tn].diagn = 1;
    master->trinfo_3d[tn].advect = 0;
    master->trinfo_3d[tn].diffuse = 0;
    master->trinfo_3d[tn].inwc = 1;
    master->trinfo_3d[tn].insed = 0;
    master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
    master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
    master->trinfo_3d[tn].m = -1;
    master->trinfo_3d[tn].n = tn;
    tn++;
  }
  if (!(params->decf & (NONE|DEC_ETA))) {
    if (tracer_find_index("decorr_e1", master->ntr, master->trinfo_3d) == -1) {
      master->decv1 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "decorr_e1");
      strcpy(master->trinfo_3d[tn].long_name, "Decorrelation length scale");
      strcpy(master->trinfo_3d[tn].units, params->decs);
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0.0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->numbers & DUMMIES) {
    if (tracer_find_index("tracer1", master->ntr, master->trinfo_3d) == -1) {
      master->dum1 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tracer1");
      strcpy(master->trinfo_3d[tn].long_name, "Dummy tracer 1");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("tracer2", master->ntr, master->trinfo_3d) == -1) {
      master->dum2 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tracer2");
      strcpy(master->trinfo_3d[tn].long_name, "Dummy tracer 2");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("tracer3", master->ntr, master->trinfo_3d) == -1) {
      master->dum3 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "tracer3");
      strcpy(master->trinfo_3d[tn].long_name, "Dummy tracer 3");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->numbers & UNIT) {
    if (tracer_find_index("unit", master->ntr, master->trinfo_3d) == -1) {
      master->unit = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "unit");
      strcpy(master->trinfo_3d[tn].long_name, "Unit tracer");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 1.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->numbers & PASS) {
    if (tracer_find_index("passive", master->ntr, master->trinfo_3d) == -1) {
      strcpy(master->trinfo_3d[tn].name, "passive");
      strcpy(master->trinfo_3d[tn].long_name, "Passive tracer");
      strcpy(master->trinfo_3d[tn].units, "kgm-3");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 1;
      master->trinfo_3d[tn].diffuse = 1;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->numbers & GLIDER) {
    if (tracer_find_index("glider", master->ntr, master->trinfo_3d) == -1) {
      master->glider = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "glider");
      strcpy(master->trinfo_3d[tn].long_name, "Glider density");
      strcpy(master->trinfo_3d[tn].units, "kgm-3");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->numbers1 & U1VHC) {
    if (tracer_find_index("u1vhc", master->ntr, master->trinfo_3d) == -1) {
      master->u1vhc = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "u1vhc");
      strcpy(master->trinfo_3d[tn].long_name, "Horizontal viscosity cell centered");
      strcpy(master->trinfo_3d[tn].units, "m2s-1");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (strlen(params->nprof)) {
    if (tracer_find_index("nprof", master->ntr, master->trinfo_3d) == -1) {
      master->nprof = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "nprof");
      sprintf(master->trinfo_3d[tn].long_name, "Normalized  vertical profile of %s", params->nprof);
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (strlen(params->monotr)) {
    if (tracer_find_index("mono", master->ntr, master->trinfo_3d) == -1) {
      master->mono = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "mono");
      sprintf(master->trinfo_3d[tn].long_name, "Monotinicity of %s", params->monotr);
      strcpy(master->trinfo_3d[tn].units, " ");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 100;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->dhwf) {
    if (tracer_find_index("dhw", master->ntr, master->trinfo_3d) == -1) {
      master->dhw = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "dhw");
      strcpy(master->trinfo_3d[tn].long_name, "Degree heating week");
      strcpy(master->trinfo_3d[tn].units, "DegC-week");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      if (params->dhwf & DHW_RT)
	strcpy(master->trinfo_3d[tn].tracerstat, "exposure(temp:dhwc:dhwt)");
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (params->dhwf & DHW_RT)
      strcpy(buf, "dhwt");
    if (params->dhwf & DHW_NOAA)
      strcpy(buf, "dhd");
    if (tracer_find_index(buf, master->ntr, master->trinfo_3d) == -1) {

      master->dhd = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, buf);

      if (params->dhwf & DHW_RT)
	strcpy(buf, "Degree heating exposure time");
      if (params->dhwf & DHW_NOAA)
	strcpy(buf, "Degree heating day");
      strcpy(master->trinfo_3d[tn].long_name, buf);

      if (params->dhwf & DHW_RT)
	strcpy(buf, "week");
      if (params->dhwf & DHW_NOAA)
	strcpy(buf, "DegCday");
      strcpy(master->trinfo_3d[tn].units, buf);

      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 7;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("dhwc", master->ntr, master->trinfo_3d) == -1) {
      master->dhwc = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "dhwc");
      strcpy(master->trinfo_3d[tn].long_name, "Degree heating threshold");
      strcpy(master->trinfo_3d[tn].units, "Degrees C");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      strcpy(master->trinfo_3d[tn].reset_file, params->dhw);
      strcpy(master->trinfo_3d[tn].reset_dt, "1 day");
      tn++;
    }
  }
  if (params->porusplate) {
    if (tracer_find_index("reef_fraction_e1", master->ntr, master->trinfo_3d) == -1) {
      master->reefe1 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "reef_fraction_e1");
      strcpy(master->trinfo_3d[tn].long_name, "Fraction of reef in cell e1 direction");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      /*tr_dataset(params->reef_frac, &master->trinfo_3d[tn], 0.0);*/
      tn++;
    }
    if (tracer_find_index("reef_fraction_e2", master->ntr, master->trinfo_3d) == -1) {
      master->reefe2 = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "reef_fraction_e2");
      strcpy(master->trinfo_3d[tn].long_name, "Fraction of reef in cell e2 direction");
      strcpy(master->trinfo_3d[tn].units, "");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      /*tr_dataset(params->reef_frac, &master->trinfo_3d[tn], 0.0);*/
      tn++;
    }
  }
  if (params->riverflow == 2) {
    if (tracer_find_index("flow_salt", master->ntr, master->trinfo_3d) == -1) {
      master->riversalt = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "flow_salt");
      strcpy(master->trinfo_3d[tn].long_name, "River flow salinity");
      strcpy(master->trinfo_3d[tn].units, "psu");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 1;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }
  if (params->swr_type & SWR_3D && strlen(params->swr_attn)) {
    if (tracer_find_index("swr_attenuation", master->ntr, master->trinfo_3d) == -1) {
      master->swr_attn = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "swr_attenuation");
      strcpy(master->trinfo_3d[tn].long_name, "SWR attenuation");
      strcpy(master->trinfo_3d[tn].units, "m-1");
      master->trinfo_3d[tn].type = WATER|HYDRO|PARAMETER;
      trn_dataset(params->swr_attn, master->trinfo_3d, tn, master->ntr, master->atr, master->tr_wc, 0.073);
      master->trinfo_3d[tn].valid_range_wc[0] = 0;
      master->trinfo_3d[tn].valid_range_wc[1] = 10;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].n = tn;
      tn++;
    } else {
      int tm = tracer_find_index("swr_attenuation", master->ntr, master->trinfo_3d);
      trn_dataset(params->swr_attn, master->trinfo_3d, tm, master->ntr, master->atr, master->tr_wc, 0.073);
    }
  }
  if (params->trasc & LAGRANGE) {
    /*
    if (tracer_find_index("Vi", master->ntr, master->trinfo_3d) == -1) {
      master->Vi = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "Vi");
      strcpy(master->trinfo_3d[tn].long_name, "Volume error");
      strcpy(master->trinfo_3d[tn].units, "m^3");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    */
  }
  if (params->fillf & LOCAL && params->fillf & DIAGNOSE) {
    if (tracer_find_index("Vcorr", master->ntr, master->trinfo_3d) == -1) {
      master->vcorr = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "Vcorr");
      strcpy(master->trinfo_3d[tn].long_name, "Volume correction");
      strcpy(master->trinfo_3d[tn].units, "%");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
    if (tracer_find_index("Acorr", master->ntr, master->trinfo_3d) == -1) {
      master->acorr = master->tr_wc[tn];
      strcpy(master->trinfo_3d[tn].name, "Acorr");
      strcpy(master->trinfo_3d[tn].long_name, "Aij correction");
      strcpy(master->trinfo_3d[tn].units, "%");
      master->trinfo_3d[tn].fill_value_wc = 0.0;
      master->trinfo_3d[tn].type = WATER|HYDRO|DIAGNOSTIC;
      master->trinfo_3d[tn].diagn = 0;
      master->trinfo_3d[tn].advect = 0;
      master->trinfo_3d[tn].diffuse = 0;
      master->trinfo_3d[tn].inwc = 1;
      master->trinfo_3d[tn].insed = 0;
      master->trinfo_3d[tn].valid_range_wc[0] = -1e10;
      master->trinfo_3d[tn].valid_range_wc[1] = 1e10;
      master->trinfo_3d[tn].m = -1;
      master->trinfo_3d[tn].n = tn;
      tn++;
    }
  }

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    /* Set up sediment tracers if required */
    tn = sediment_autotracer_3d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, master->trinfo_3d, master->ntr, tn);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    /* Set up ecology tracers if required */
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

  int t;
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];
  char tag[MAXSTRLEN];


  /* Load up the default tracer values */
  load_sed_tracer_defaults(master);

  prm_set_errfn(hd_silent_warn);

  for (t = 0; t < master->nsed; ++t) {
    int tm = master->trinfo_sed[t].m;

    sprintf(tag, "TRACER%1.1d.interp_type", tm);
    if (!(prm_read_char(fp, tag, buf2)))
      buf2[0] = '\0';
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


/*-------------------------------------------------------------------*/
/* Routine to create the 3D tracers in the params tracer_info.       */
/* Defaults for auto runmodes are set in this routine.               */
/*-------------------------------------------------------------------*/
void create_tracer_3d(parameters_t *params)   /* Input parameters    */

{
  tracer_info_t *trinfo = params->trinfo_3d;
  char buf[MAXSTRLEN];
  int tn = 0;

  /* Set up the tracer_info structure with mandatory and optional   */
  /* tracers.                                                       */
  if ((tracer_find_index("salt", params->ntr, params->trinfo_3d)) < 0) {
    strcpy(trinfo[tn].name, "salt");
    strcpy(trinfo[tn].long_name, "Salinity");
    strcpy(trinfo[tn].units, "PSU");
    /*UR added */
    strcpy(trinfo[tn].std_name, "sea_water_salinity");
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].fill_value_wc = 35.0;
    trinfo[tn].type = WATER;
    trinfo[tn].inwc = 1;
    trinfo[tn].dissol = 1;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 1;
    strcpy(trinfo[tn].decay, "0.0");
    trinfo[tn].diagn = 0;
    trinfo[tn].m = tn;
    tn++;
  }
  if ((tracer_find_index("temp", params->ntr, params->trinfo_3d)) < 0) {
    strcpy(trinfo[tn].name, "temp");
    strcpy(trinfo[tn].long_name, "Temperature");
    strcpy(trinfo[tn].units, "degrees C");
    /*UR added */
    strcpy(trinfo[tn].std_name, "sea_water_temperature");
    trinfo[tn].valid_range_wc[0] = -4;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].fill_value_wc = 20.0;
    trinfo[tn].type = WATER;
    trinfo[tn].inwc = 1;
    trinfo[tn].dissol = 1;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 1;
    strcpy(trinfo[tn].decay, "0.0");
    trinfo[tn].diagn = 0;
    trinfo[tn].m = tn;
    tn++;
  }
  /* 3D mean velocity */
  if (params->means & VEL3D) {
    if ((tracer_find_index("u1mean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1mean");
      strcpy(trinfo[tn].long_name, "Mean 3D east velocity");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC|E1VAR;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -100;
      trinfo[tn].valid_range_wc[1] = 100;
      strcpy(trinfo[tn].vector_name, "Mean 3D current");
      strcpy(trinfo[tn].vector_components, "u1mean u2mean");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2mean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2mean");
      strcpy(trinfo[tn].long_name, "Mean 3D north velocity");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC|E2VAR;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -100;
      trinfo[tn].valid_range_wc[1] = 100;
      strcpy(trinfo[tn].vector_name, "Mean 3D current");
      strcpy(trinfo[tn].vector_components, "u1mean u2mean");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("wmean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "wmean");
      strcpy(trinfo[tn].long_name, "Mean w velocity");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -100;
      trinfo[tn].valid_range_wc[1] = 100;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  if (params->means & TS) {
    if ((tracer_find_index("temp_mean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "temp_mean");
      strcpy(trinfo[tn].long_name, "Mean temperature");
      strcpy(trinfo[tn].units, "deg C");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e10;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("salt_mean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "salt_mean");
      strcpy(trinfo[tn].long_name, "Mean salinity");
      strcpy(trinfo[tn].units, "psu");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e10;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  if (params->means & MTRA3D) {
    if ((tracer_find_index("tracer_mean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "tracer_mean");
      sprintf(buf, "Mean %s", params->means_tra);
      strcpy(trinfo[tn].long_name, buf);
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e10;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* 3D mean volume flux */
  /*
  if ((params->means & VOLFLUX) || (params->tmode & SP_FFSL)) {
    if ((tracer_find_index("u1vmean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1vmean");
      strcpy(trinfo[tn].long_name, "Mean u1 volume flux");
      strcpy(trinfo[tn].units, "m3s-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC|E1VAR;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e30;
      trinfo[tn].valid_range_wc[1] = 1e30;
      strcpy(trinfo[tn].vector_name, "Mean volume flux");
      strcpy(trinfo[tn].vector_components, "u1vmean u2vmean");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2vmean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2vmean");
      strcpy(trinfo[tn].long_name, "Mean u2 volume flux");
      strcpy(trinfo[tn].units, "m3s-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC|E2VAR;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e30;
      trinfo[tn].valid_range_wc[1] = 1e30;
      strcpy(trinfo[tn].vector_name, "Mean volume flux");
      strcpy(trinfo[tn].vector_components, "u1vmean u2vmean");
      trinfo[tn].m = tn;
      tn++;
    }
  }
  */
  /* Mean Kz */
  if (params->means & KZ_M) {
    if ((tracer_find_index("Kzmean", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "Kzmean");
      strcpy(trinfo[tn].long_name, "Mean Vertical Diffusivity");
      strcpy(trinfo[tn].units, "m2s-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Tracer fluxes */
  if (strcmp(params->trflux, "NONE") != 0) {
    if ((tracer_find_index("flux_e1", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "flux_e1");
      sprintf(trinfo[tn].long_name, "Advective flux through face %d", params->trfd1);
      strcpy(trinfo[tn].units, "kgs-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("flux_e2", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "flux_e2");
      sprintf(trinfo[tn].long_name, "Advective flux through face %d", params->trfd2);
      strcpy(trinfo[tn].units, "kgs-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("flux_w", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "flux_w");
      strcpy(trinfo[tn].long_name, "Vertical advective flux");
      strcpy(trinfo[tn].units, "kgs-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("flux_kz", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "flux_kz");
      strcpy(trinfo[tn].long_name, "Vertical diffusive flux");
      strcpy(trinfo[tn].units, "kgs-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Regions */
  if (strlen(params->regions)) {
    if ((tracer_find_index("regionid", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "regionid");
      strcpy(trinfo[tn].long_name, "Region identifier");
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("residence", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "residence");
      strcpy(trinfo[tn].long_name, "Residence time");
      strcpy(trinfo[tn].units, "days");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Flushing tracer */
  if (params->trflsh) {
    if ((tracer_find_index("flush", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "flush");
      strcpy(trinfo[tn].long_name, "Flushing tracer");
      strcpy(trinfo[tn].units, "mgL-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Age tracer */
  if (strlen(params->trage)) {
    if ((tracer_find_index("age", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "age");
      strcpy(trinfo[tn].long_name, "Age tracer");
      strcpy(trinfo[tn].units, "days");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1e10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Tracer percentiles */
  if (strcmp(params->trperc, "NONE") != 0) {
    char buf[MAXSTRLEN];
    sprintf(buf, "percentile_%s", params->trperc);
    if (tracer_find_index(buf, params->ntr, params->trinfo_3d) < 0) {
      strcpy(trinfo[tn].name, buf);
      sprintf(buf, "Percentile for %s", params->trperc);
      strcpy(trinfo[tn].long_name, buf);
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 100;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* k-e closure */
  if (strcmp(params->mixsc, "k-e") == 0) {
    if ((tracer_find_index("tke", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "tke");
      strcpy(trinfo[tn].long_name, "Turbulent Kinetic Energy");
      strcpy(trinfo[tn].units, "m2s-2");
      trinfo[tn].fill_value_wc = 7.6e-6;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 10;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("diss", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "diss");
      strcpy(trinfo[tn].long_name, "Dissipation");
      strcpy(trinfo[tn].units, "m2s-3");
      trinfo[tn].fill_value_wc = 5e-10;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 10;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* k-w closure */
  if (strcmp(params->mixsc, "k-w") == 0 || strcmp(params->mixsc, "W88") == 0) {
    if ((tracer_find_index("tke", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "tke");
      strcpy(trinfo[tn].long_name, "Turbulent Kinetic Energy");
      strcpy(trinfo[tn].units, "m2s-2");
      trinfo[tn].fill_value_wc = 7.6e-6;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 100;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("omega", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "omega");
      strcpy(trinfo[tn].long_name, "Turbulence frequency");
      strcpy(trinfo[tn].units, "s-1");
      trinfo[tn].fill_value_wc = 1.0e-12;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 100;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Mellor-Yamada 2.5 closure */
  if (strcmp(params->mixsc, "mellor_yamada_2_5") == 0 ||
      strcmp(params->mixsc, "harcourt") == 0) {
    if ((tracer_find_index("tki", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "tki");
      strcpy(trinfo[tn].long_name, "Turbulent Kinetic Intensity");
      strcpy(trinfo[tn].units, "m2s-2");
      trinfo[tn].fill_value_wc = 2.0e-8;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 10;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("tki_l", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "tki_l");
      strcpy(trinfo[tn].long_name,
	     "Turbulent Kinetic Intensity Length Scale");
      strcpy(trinfo[tn].units, "m2s-1");
      trinfo[tn].fill_value_wc = 3.4e-9;
      trinfo[tn].type = WATER|HYDRO|PROGNOSTIC|CLOSURE;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 1;
      trinfo[tn].diffuse = 1;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 100;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("lscale", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "lscale");
      strcpy(trinfo[tn].long_name, "Turbulence length scale");
      strcpy(trinfo[tn].units, "metre");
      trinfo[tn].fill_value_wc = 0.17;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1e4;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("Kq", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "Kq");
      strcpy(trinfo[tn].long_name, "Turbulence mixing");
      strcpy(trinfo[tn].units, "m2s-1");
      trinfo[tn].fill_value_wc = 1e-5;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].inwc = 1;
      trinfo[tn].insed = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1; 
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Smagorinsky diffusion */
  if (params->smagorinsky > 0.0) {
    if ((tracer_find_index("smagorinsky", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "smagorinsky");
      strcpy(trinfo[tn].long_name, "Smagorinsky diffusion");
      strcpy(trinfo[tn].units, "m2s-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Layer thickness */
  if (params->show_layers) {
    if ((tracer_find_index("layer_thick", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "layer_thick");
      strcpy(trinfo[tn].long_name, "Layer thickness");
      strcpy(trinfo[tn].units, "metre");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 1;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = 0;
      trinfo[tn].valid_range_wc[1] = 5000;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Momentum tendencies */
  if (params->tendf) {
    if ((tracer_find_index("mom_balance", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "mom_balance");
      strcpy(trinfo[tn].long_name, "Momentum balance maximum");
      strcpy(trinfo[tn].units, "%");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_adv", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_adv");
      strcpy(trinfo[tn].long_name, "u1 advective tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Advective current");
      strcpy(trinfo[tn].vector_components, "u1_adv u2_adv");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_hdif", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_hdif");
      strcpy(trinfo[tn].long_name, "u1 horizontal diffusion tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Hdiffusion current");
      strcpy(trinfo[tn].vector_components, "u1_hdif u2_hdif");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_vdif", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_vdif");
      strcpy(trinfo[tn].long_name, "u1 vertical diffusion tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Vdiffusion current");
      strcpy(trinfo[tn].vector_components, "u1_vdif u2_vdif");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_btp", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_btp");
      strcpy(trinfo[tn].long_name, "u1 barotropic pressure gradient tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Barotropic current");
      strcpy(trinfo[tn].vector_components, "u1_btp u2_btp");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_bcp", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_bcp");
      strcpy(trinfo[tn].long_name, "u1 baroclinic pressure tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Baroclinic current");
      strcpy(trinfo[tn].vector_components, "u1_bcp u2_bcp");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u1_cor", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u1_cor");
      strcpy(trinfo[tn].long_name, "u1 Coriolis tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Coriolis current");
      strcpy(trinfo[tn].vector_components, "u1_cor u2_cor");
      trinfo[tn].m = tn;
      tn++;
    }
    if (params->waves & STOKES_DRIFT) {
      if ((tracer_find_index("u1_sto", params->ntr, params->trinfo_3d)) < 0) {
	strcpy(trinfo[tn].name, "u1_sto");
	strcpy(trinfo[tn].long_name, "u1 Stokes tendency");
	strcpy(trinfo[tn].units, "ms-1");
	trinfo[tn].fill_value_wc = 0.0;
	trinfo[tn].type = WATER;
	trinfo[tn].diagn = 0;
	trinfo[tn].advect = 0;
	trinfo[tn].diffuse = 0;
	trinfo[tn].valid_range_wc[0] = -1e35;
	trinfo[tn].valid_range_wc[1] = 1e35;
	strcpy(trinfo[tn].vector_name, "Stokes current");
	strcpy(trinfo[tn].vector_components, "u1_sto u2_sto");
	trinfo[tn].m = tn;
	tn++;
      }
    }
    if ((tracer_find_index("u2_adv", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_adv");
      strcpy(trinfo[tn].long_name, "u2 advective tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Advective current");
      strcpy(trinfo[tn].vector_components, "u1_adv u2_adv");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2_hdif", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_hdif");
      strcpy(trinfo[tn].long_name,
	     "u2 horizontal diffusion tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Hdiffusion current");
      strcpy(trinfo[tn].vector_components, "u1_hdif u2_hdif");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2_vdif", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_vdif");
      strcpy(trinfo[tn].long_name,
	     "u2 vertical diffusion tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Vdiffusion current");
      strcpy(trinfo[tn].vector_components, "u1_vdif u2_vdif");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2_btp", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_btp");
      strcpy(trinfo[tn].long_name,
	     "u2 barotropic pressure gradient tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Barotropic current");
      strcpy(trinfo[tn].vector_components, "u1_btp u2_btp");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2_bcp", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_bcp");
      strcpy(trinfo[tn].long_name,
	     "u2 baroclinic pressure tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Baroclinic current");
      strcpy(trinfo[tn].vector_components, "u1_bcp u2_bcp");
      trinfo[tn].m = tn;
      tn++;
    }
    if ((tracer_find_index("u2_cor", params->ntr, params->trinfo_3d)) < 0) {
      strcpy(trinfo[tn].name, "u2_cor");
      strcpy(trinfo[tn].long_name, "u2 Coriolis tendency");
      strcpy(trinfo[tn].units, "ms-1");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      strcpy(trinfo[tn].vector_name, "Coriolis current");
      strcpy(trinfo[tn].vector_components, "u1_cor u2_cor");
      trinfo[tn].m = tn;
      tn++;
    }
    if (params->waves & STOKES_DRIFT) {
      if ((tracer_find_index("u2_sto", params->ntr, params->trinfo_3d)) < 0) {
	strcpy(trinfo[tn].name, "u2_sto");
	strcpy(trinfo[tn].long_name, "u2 Stokes tendency");
	strcpy(trinfo[tn].units, "ms-1");
	trinfo[tn].fill_value_wc = 0.0;
	trinfo[tn].type = WATER;
	trinfo[tn].diagn = 0;
	trinfo[tn].advect = 0;
	trinfo[tn].diffuse = 0;
	trinfo[tn].valid_range_wc[0] = -1e35;
	trinfo[tn].valid_range_wc[1] = 1e35;
	strcpy(trinfo[tn].vector_name, "Stokes current");
	strcpy(trinfo[tn].vector_components, "u1_sto u2_sto");
	trinfo[tn].m = tn;
	tn++;
      }
    }
  }
  /* Tracer tendencies */
  if (strlen(params->trtend)) {
    if (tracer_find_index("tra_adv", params->ntr, params->trinfo_3d) == -1) {
      strcpy(trinfo[tn].name, "tra_adv");
      strcpy(trinfo[tn].long_name, "Tracer advective tendency");
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    /*
    if (tracer_find_index("tra_hdif", params->ntr, params->trinfo_3d) == -1) {
      strcpy(trinfo[tn].name, "tra_hdif");
      strcpy(trinfo[tn].long_name, "Tracer horizontal diffusive tendency");
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    */
    if (tracer_find_index("tra_vdif", params->ntr, params->trinfo_3d) == -1) {
      strcpy(trinfo[tn].name, "tra_vdif");
      strcpy(trinfo[tn].long_name, "Tracer vertical diffusive tendency");
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
    if (tracer_find_index("tra_ncon", params->ntr, params->trinfo_3d) == -1) {
      strcpy(trinfo[tn].name, "tra_ncon");
      strcpy(trinfo[tn].long_name, "Tracer non-conservative tendency");
      strcpy(trinfo[tn].units, "");
      trinfo[tn].fill_value_wc = 0.0;
      trinfo[tn].type = WATER;
      trinfo[tn].diagn = 0;
      trinfo[tn].advect = 0;
      trinfo[tn].diffuse = 0;
      trinfo[tn].valid_range_wc[0] = -1e35;
      trinfo[tn].valid_range_wc[1] = 1e35;
      trinfo[tn].m = tn;
      tn++;
    }
  }
  /* Waves */
  if (params->waves & SPECTRAL) {
    strcpy(trinfo[tn].name, "wave_stke1");
    strcpy(trinfo[tn].long_name, "Stokes sub-surface velocity along e1");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER|WAVE|FORCING;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -100.0;
    trinfo[tn].valid_range_wc[1] = 100.0;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "wave_stke2");
    strcpy(trinfo[tn].long_name, "Stokes sub-surface velocity along e2");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER|WAVE|FORCING;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -100.0;
    trinfo[tn].valid_range_wc[1] = 100.0;
    trinfo[tn].m = tn;
    tn++;
  }
  /* Diagnostic numbers */
  if (params->numbers & BRUNT) {
    strcpy(trinfo[tn].name, "brunt_vaisala");
    strcpy(trinfo[tn].long_name, "Brunt Vaisala Frequency");
    strcpy(trinfo[tn].units, "s-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & INT_WAVE) {
    strcpy(trinfo[tn].name, "int_wave_speed");
    strcpy(trinfo[tn].long_name, "Internal wave speed");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & RICHARD_GR) {
    strcpy(trinfo[tn].name, "richardson_gr");
    strcpy(trinfo[tn].long_name, "Gradient Richardson number");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e308;
    trinfo[tn].valid_range_wc[1] = 1e308;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & RICHARD_FL) {
    strcpy(trinfo[tn].name, "richardson_fl");
    strcpy(trinfo[tn].long_name, "Flux Richardson number");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e308;
    trinfo[tn].valid_range_wc[1] = 1e308;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & REYNOLDS) {
    strcpy(trinfo[tn].name, "reynolds");
    strcpy(trinfo[tn].long_name, "Reynolds number");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & FROUDE) {
    strcpy(trinfo[tn].name, "froude");
    strcpy(trinfo[tn].long_name, "Froude number");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & SIGMA_T) {
    strcpy(trinfo[tn].name, "sigma_t");
    strcpy(trinfo[tn].long_name, "Sigma_t");
    strcpy(trinfo[tn].units, "kgm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & ENERGY) {
    strcpy(trinfo[tn].name, "energy");
    strcpy(trinfo[tn].long_name, "Mechanical energy");
    strcpy(trinfo[tn].units, "Jm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & KINETIC) {
    strcpy(trinfo[tn].name, "kenergy");
    strcpy(trinfo[tn].long_name, "Kinetic energy");
    strcpy(trinfo[tn].units, "Jm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->do_pt) {
    strcpy(trinfo[tn].name, "ptconc");
    strcpy(trinfo[tn].long_name, "Particle concentration");
    strcpy(trinfo[tn].units, "kgm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & SOUND) {
    strcpy(trinfo[tn].name, "sound");
    strcpy(trinfo[tn].long_name, "Speed of sound");
    strcpy(trinfo[tn].units, "ms-1");
    /*UR added */
    strcpy(trinfo[tn].std_name, "speed_of_sound_in_sea_water");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 1000.;
    trinfo[tn].valid_range_wc[1] = 2000.;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "sound_channel");
    strcpy(trinfo[tn].long_name, "Sound channel depth");
    strcpy(trinfo[tn].units, "m");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e4;
/*UR made 100 as not-realistic but theoretically possible limit */
    trinfo[tn].valid_range_wc[1] = 100;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & ROSSBY_IN) {
    strcpy(trinfo[tn].name, "rossby_internal");
    strcpy(trinfo[tn].long_name, "Internal Rossby radius");
    strcpy(trinfo[tn].units, "metre");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & SPEED_3D) {
    strcpy(trinfo[tn].name, "current_speed_3d");
    strcpy(trinfo[tn].long_name, "Current Speed 3D");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -100;
    trinfo[tn].valid_range_wc[1] = 100;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & SHEAR_V) {
    strcpy(trinfo[tn].name, "shear_vert");
    strcpy(trinfo[tn].long_name, "Vertical velocity shear");
    strcpy(trinfo[tn].units, "s-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & BUOY_PROD) {
    strcpy(trinfo[tn].name, "buoy_prod");
    strcpy(trinfo[tn].long_name, "Buoyancy production");
    strcpy(trinfo[tn].units, "m2s-2");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & SHEAR_PROD) {
    strcpy(trinfo[tn].name, "shear_prod");
    strcpy(trinfo[tn].long_name, "Shear production");
    strcpy(trinfo[tn].units, "m2s-2");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (!(params->decf & (NONE|DEC_ETA))) {
    strcpy(trinfo[tn].name, "decorr_e1");
    strcpy(trinfo[tn].long_name, "Decorrelation length scale e1");
    strcpy(trinfo[tn].units, params->decs);
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "decorr_e2");
    strcpy(trinfo[tn].long_name, "Decorrelation length scale e2");
    strcpy(trinfo[tn].units, params->decs);
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->save_force & OTEMP) {
    strcpy(trinfo[tn].name, "otemp");
    strcpy(trinfo[tn].long_name, "OFAM temperature");
    strcpy(trinfo[tn].units, "degrees C");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -4;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].m = tn;
    if (strlen(params->odata)) {
      if (params->save_force & ROAM)
	sprintf(buf, "%s(otemp=temp)", params->tdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->save_force & OSALT) {
    strcpy(trinfo[tn].name, "osalt");
    strcpy(trinfo[tn].long_name, "OFAM salinity");
    strcpy(trinfo[tn].units, "PSU");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].m = tn;
    if (strlen(params->odata)) {
      if (params->save_force & ROAM)
	sprintf(buf, "%s(osalt=salt)", params->sdata);
      else
	sprintf(buf, "%s", params->odata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->rtemp) {
    strcpy(trinfo[tn].name, "rtemp");
    strcpy(trinfo[tn].long_name, "Relaxation temperature");
    strcpy(trinfo[tn].units, "degrees C");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -4;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].m = tn;
    if (strlen(params->tdata)) {
      char buf[MAXSTRLEN];
      char buf1[MAXSTRLEN];
      char buf2[MAXSTRLEN];
      if (params->save_force & ROAM) {
	/* See if temp is being variable substituted */
	if (find_token(params->tdata, "temp=", buf2, ')'))
	  sprintf(buf, "%s(rtemp=%s)", params->tdata, buf2);
	else 
	  sprintf(buf, "%s(rtemp=temp)", params->tdata);
	strcpy(buf1, "1 day");
      }
      else {
	sprintf(buf, "%s", params->tdata);
	strcpy(buf1, "1 second");
      }
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, buf1);
    }
    tn++;
  }
  if (params->rsalt) {
    strcpy(trinfo[tn].name, "rsalt");
    strcpy(trinfo[tn].long_name, "Relaxation salinity");
    strcpy(trinfo[tn].units, "PSU");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 40;
    trinfo[tn].m = tn;
    if (strlen(params->sdata)) {
      char buf[MAXSTRLEN];
      char buf1[MAXSTRLEN];
      char buf2[MAXSTRLEN];
      if (params->save_force & ROAM) {
	/* See if salt is being variable substituted */
	if (find_token(params->sdata, "salt=", buf2, ')'))
	  sprintf(buf, "%s(rsalt=%s)", params->sdata, buf2);
	else
	  sprintf(buf, "%s(rsalt=salt)", params->sdata);
	strcpy(buf1, "1 day");
      }
      else {
	sprintf(buf, "%s", params->sdata);
	strcpy(buf1, "1 second");
      }
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, buf1);
    }
    tn++;
  }
  if (params->rtemp & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    strcpy(trinfo[tn].name, "temp_tc");
    strcpy(trinfo[tn].long_name, "Relaxation temperature time constant");
    strcpy(trinfo[tn].units, "days");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->rsalt & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    strcpy(trinfo[tn].name, "salt_tc");
    strcpy(trinfo[tn].long_name, "Relaxation salinity time constant");
    strcpy(trinfo[tn].units, "days");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->save_force & OVELU) {
    strcpy(trinfo[tn].name, "ovelu");
    strcpy(trinfo[tn].long_name, "OFAM east velocity");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -100;
    trinfo[tn].valid_range_wc[1] = 100;
    strcpy(trinfo[tn].vector_name, "Global 3D current");
    strcpy(trinfo[tn].vector_components, "ovelu ovelv");
    trinfo[tn].m = tn;
    if (strlen(params->vdata)) {
      if (params->save_force & ROAM) {
	char buf2[MAXSTRLEN];
	/* See if u is being variable substituted */
	if (find_token(params->vdata, "u=", buf2, ')'))
	  sprintf(buf, "%s(ovelu=%s)", params->vdata, buf2);
	else
	  sprintf(buf, "%s(ovelu=u)", params->vdata);
      } else
	sprintf(buf, "%s", params->vdata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->save_force & OVELV) {
    strcpy(trinfo[tn].name, "ovelv");
    strcpy(trinfo[tn].long_name, "OFAM north velocity");
    strcpy(trinfo[tn].units, "ms-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -100;
    trinfo[tn].valid_range_wc[1] = 100;
    strcpy(trinfo[tn].vector_name, "Global 3D current");
    strcpy(trinfo[tn].vector_components, "ovelu ovelv");
    trinfo[tn].m = tn;
    if (strlen(params->vdata)) {
      if (params->save_force & ROAM) {
	char buf2[MAXSTRLEN];
	/* See if v is being variable substituted */
	if (find_token(params->vdata, "v=", buf2, ')'))
	  sprintf(buf, "%s(ovelv=%s)", params->vdata, buf2);
	else
	  sprintf(buf, "%s(ovelv=v)", params->vdata);
      } else
	sprintf(buf, "%s", params->vdata);
      strcpy(trinfo[tn].reset_file, buf);
      strcpy(trinfo[tn].reset_dt, "1 second");
    }
    tn++;
  }
  if (params->porusplate) {
    strcpy(trinfo[tn].name, "reef_fraction_e1");
    strcpy(trinfo[tn].long_name, "Fraction of reef in cell e1 direction");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    /*tr_dataset(params->reef_frac, &params->trinfo_3d[tn], 0.0);*/
    tn++;
    strcpy(trinfo[tn].name, "reef_fraction_e2");
    strcpy(trinfo[tn].long_name, "Fraction of reef in cell e2 direction");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    /*tr_dataset(params->reef_frac, &params->trinfo_3d[tn], 0.0);*/
    tn++;
  }
  if (params->numbers & DUMMIES) {
    strcpy(trinfo[tn].name, "tracer1");
    strcpy(trinfo[tn].long_name, "Dummy tracer 1");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "tracer2");
    strcpy(trinfo[tn].long_name, "Dummy tracer 2");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "tracer3");
    strcpy(trinfo[tn].long_name, "Dummy tracer 3");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & UNIT) {
    strcpy(trinfo[tn].name, "unit");
    strcpy(trinfo[tn].long_name, "Unit tracer");
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 1.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & PASS) {
    strcpy(trinfo[tn].name, "passive");
    strcpy(trinfo[tn].long_name, "Passive tracer");
    strcpy(trinfo[tn].units, "kgm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 1;
    trinfo[tn].diffuse = 1;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers & GLIDER) {
    strcpy(trinfo[tn].name, "glider");
    strcpy(trinfo[tn].long_name, "Glider density");
    strcpy(trinfo[tn].units, "kgm-3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->numbers1 & U1VHC) {
    strcpy(trinfo[tn].name, "u1vhc");
    strcpy(trinfo[tn].long_name, "Horizontal viscosity cell centered");
    strcpy(trinfo[tn].units, "m2s-1");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  if (strlen(params->nprof)) {
    strcpy(trinfo[tn].name, "nprof");
    sprintf(trinfo[tn].long_name, "Normalized  vertical profile of %s", params->nprof);
    strcpy(trinfo[tn].units, "");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 100;
    trinfo[tn].m = tn;
    tn++;
  }
  if (strlen(params->monotr)) {
    strcpy(trinfo[tn].name, "mono");
    sprintf(trinfo[tn].long_name, "Monotinicity of %s", params->monotr);
    strcpy(trinfo[tn].units, " ");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 100;
    trinfo[tn].m = tn;
    tn++;
  }
  if (params->dhwf) {
    strcpy(trinfo[tn].name, "dhw");
    strcpy(trinfo[tn].long_name, "Degree heating week");
    strcpy(trinfo[tn].units, "DegC-week");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    if (params->dhwf & DHW_RT)
      strcpy(trinfo[tn].tracerstat, "exposure(temp:dhwc:dhwt)");
    trinfo[tn].m = -1;
    trinfo[tn].n = tn;
    tn++;

    if (params->dhwf & DHW_RT)
      strcpy(buf, "dhwt");
    if (params->dhwf & DHW_NOAA)
      strcpy(buf, "dhd");
    strcpy(trinfo[tn].name, buf);
    if (params->dhwf & DHW_RT)
      strcpy(buf, "Degree heating exposure time");
    if (params->dhwf & DHW_NOAA)
      strcpy(buf, "Degree heating day");
    strcpy(trinfo[tn].long_name, buf);
    if (params->dhwf & DHW_RT)
      strcpy(buf, "week");
    if (params->dhwf & DHW_NOAA)
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

    strcpy(trinfo[tn].name, "dhwc");
    strcpy(trinfo[tn].long_name, "Degree heating threshold");
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
    strcpy(trinfo[tn].reset_file, params->dhw);
    strcpy(trinfo[tn].reset_dt, "1 day");
    tn++;
  }
  if (params->riverflow == 2) {
    strcpy(trinfo[tn].name, "flow_salt");
    strcpy(trinfo[tn].long_name, "River flow salinity");
    strcpy(trinfo[tn].units, "psu");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER|HYDRO|DIAGNOSTIC;
    trinfo[tn].diagn = 1;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].n = tn;
    tn++;
  }
  if (params->swr_type & SWR_3D && strlen(params->swr_attn)) {
    strcpy(trinfo[tn].name, "swr_attenuation");
    strcpy(trinfo[tn].long_name, "SWR attenuation");
    strcpy(trinfo[tn].units, "m-1");
    trinfo[tn].type = WATER|HYDRO|PARAMETER;
    tr_dataset(params->swr_attn, &params->trinfo_3d[tn], 0.073);
    trinfo[tn].valid_range_wc[0] = 0;
    trinfo[tn].valid_range_wc[1] = 10;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].n = tn;
    tn++;
  }
  /*
  if (params->trasc & LAGRANGE) {
    strcpy(trinfo[tn].name, "Vi");
    strcpy(trinfo[tn].long_name, "Volume error");
    strcpy(trinfo[tn].units, "m^3");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }
  */
  if (params->fillf & LOCAL && params->fillf & DIAGNOSE) {
    strcpy(trinfo[tn].name, "Vcorr");
    strcpy(trinfo[tn].long_name, "Volume correction");
    strcpy(trinfo[tn].units, "%");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
    strcpy(trinfo[tn].name, "Acorr");
    strcpy(trinfo[tn].long_name, "Aij correction");
    strcpy(trinfo[tn].units, "%");
    trinfo[tn].fill_value_wc = 0.0;
    trinfo[tn].type = WATER;
    trinfo[tn].diagn = 0;
    trinfo[tn].advect = 0;
    trinfo[tn].diffuse = 0;
    trinfo[tn].inwc = 1;
    trinfo[tn].insed = 0;
    trinfo[tn].valid_range_wc[0] = -1e10;
    trinfo[tn].valid_range_wc[1] = 1e10;
    trinfo[tn].m = tn;
    tn++;
  }

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed) {
    /* Set up sediment tracers if required */
    tn = sediment_autotracer_3d(params->prmfd, params->do_sed, params->sed_vars, 
				params->sed_defs, trinfo, params->ntr, tn);
  }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco) {
    /* Set up ecology tracers if required */
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
  int nvals;
  char buf[MAXSTRLEN], key[MAXSTRLEN], sgn[MAXSTRLEN], line[MAXSTRLEN];
  double val, v1, v2, *d1 = NULL;
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
	if (value_init_sparse3d(master, ret, fname, vname, i_rule) == 1) {
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
static char *in_dims[4][8] = {
  {"botz", "i_centre", "j_centre", "k_centre", "x_centre", "y_centre", "z_centre", "t"},
  {"height", "longitude", "latitude", "zc", "lon", "lat", "zc", "time"},
  {"height", "ni", "nj", "nk", "x", "y", "z", "t"},
  {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}
};


/*-------------------------------------------------------------------*/
/* Reads structured 2D tracer values from a netCDF file, packs into  */
/* a vector and interpolates onto an unstructured mesh.              */
/*-------------------------------------------------------------------*/
int value_init_sparse2d(master_t *master,   /* Master data           */
			double *ret,        /* 2D array              */
			char *fname,        /* File name             */
			char *vname,        /* Variable name         */
			char *in_rule       /* Interpolation type    */
			)
{
  geometry_t *geom = master->geom;
  char i_rule[MAXSTRLEN];

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
int value_init_sparse3d(master_t *master,   /* Master data           */
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
      int j, cs;
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
  int cc, c, cs, k;

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
  double *d1 = NULL, dt;
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
  double dt, db, vs;
  int cc, c, cs;
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
  double dt, db, vs;
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
  int nvals;
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
	if (value_init_sparse2d(master, ret, fname, vname, i_rule) == 1) {
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
  int nvals;
  char buf[MAXSTRLEN];
  double val;
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

