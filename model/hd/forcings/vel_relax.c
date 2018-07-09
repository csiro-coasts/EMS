/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/forcings/vel_relax.c
 *  
 *  Description:
 *  Relaxation of velocity.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: vel_relax.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

double get_relax_rate(relax_info_t *rlx, double u1, double u2, double v1, double v2);

typedef struct {
  timeseries_t *ts;             /* Time-Series file */
  int vid1;                     /* df_variable_t id into TS file. */
  int vid2;                     /* df_variable_t id into TS file. */
  double dt;                    /* Relaxation time step */
  double rate;                  /* Relaxation rate */
  double tnext;                 /* Next time step */
  timeseries_t *tc;             /* Relaxation rate timeseries */
  int tcid;                     /* Relaxation rate id */
} vel_relax_data_t;


int vel_relax_init(sched_event_t *event)
{
  char buf[MAXSTRLEN];
  char fname[MAXSTRLEN];
  double dt, tconst;
  vel_relax_data_t *relax = NULL;
  master_t *master = (master_t *)schedGetPublicData(event);
  char tc_name[] = "vel_relaxation_time_constant";
  
  prm_set_errfn(hd_silent_warn);
  
  if (master->velrlx & RELAX) {
    relax_info_t *rlx = master->vel_rlx;
    
    /* Set the input time step, and time constant */
    strcpy(fname, rlx->rlxn);
    dt = rlx->rlxdt;
    
    /* Allocate memory for relaxation structure */
    relax = (vel_relax_data_t *)malloc(sizeof(vel_relax_data_t));
    memset(relax, 0, sizeof(vel_relax_data_t));
    
    /* Assign the time constant                                */
    relax->tc = NULL;
    rlx->tctype = NONE;
    if (strlen(rlx->rlxtc)) {
      char tu0[MAXSTRLEN], tu1[MAXSTRLEN];
      double r0, r1, d0, d1;
      if (sscanf(rlx->rlxtc, "%lf %s", &tconst, buf) == 2) {
	tm_scale_to_secs(rlx->rlxtc, &tconst);
	relax->rate = rlx->rate = tconst;
	rlx->tctype = RLX_CONS;
      } else if (sscanf(rlx->rlxtc, "%s %lf %lf %s %lf %lf %s", 
			buf, &d0, &r0, tu0, &d1, &r1, tu1) == 7) {
	if(strcmp(buf, "linear") == 0) {	
	  hd_warn("Linear adaptive relaxation performed on velocity.\n");
	  tconst = 1.0;
	  relax->rate = 0.0;
	  rlx->dv0 = d0;
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  rlx->tc0 = r0;
	  rlx->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	  rlx->tctype = (RLX_ADPT|RLX_LINR);
	} else
	  hd_quit("Linear adaptive relaxation format, e.g; 'linear 10 1 day 0.1 10 days\n");
      } else if (sscanf(rlx->rlxtc, "%s %lf %lf %s", buf, &d0, &r0, tu0) == 4) {
	if (strcmp(buf, "exponential") == 0) {
	  hd_warn("Exponential adaptive relaxation performed on elevation.\n");
	  tconst = 1.0;
	  relax->rate = 0.0;
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  rlx->dv0 = d0 * log(r0);
	  rlx->slope = d0 * log(r0 / 86400.0);
	  rlx->tc0 = r0;
	  rlx->tctype = (RLX_ADPT|RLX_EXP);
	} else
	  hd_quit("Exponential adaptive relaxation format, e.g; 'exponential 10 1 day'\n");
      } else {
	/* Relaxation constant is in a file */
	relax->tc = hd_ts_read(master, rlx->rlxtc, 0);
	relax->tcid = ts_get_index(relax->tc, 
				   fv_get_varname(rlx->rlxtc, tc_name, buf));
	if (relax->tcid < 0)
	  hd_quit("%s does not contain '%s'.\n", rlx->rlxtc, tc_name);
	tconst = ts_eval(relax->tc, relax->tcid, master->t);
	sprintf(buf, "%f %s\n",tconst, relax->tc->varunit[relax->tcid]);
	tm_scale_to_secs(buf, &tconst);
	rlx->rate = tconst;
	rlx->tctype = RLX_FILE;
      }
    }
    if (master->velrlx & RELAX && tconst <= 0)
      hd_quit
        ("vel_relax_init: Invalid velocity relaxation time constant specified.\n");

    /* Populate and register the scheduler events. */
    relax->dt = dt;
    if (tconst)relax->rate = rlx->rate = tconst;
    relax->tnext = schedule->start_time;
    relax->ts = hd_ts_read(master, fname, 1);
    relax->vid1 =
      ts_get_index(relax->ts, fv_get_varname(fname, "u", buf));
    if (relax->vid1 < 0) {
      hd_quit("vel_relax_init: The file '%s' does not contain 'u'\n",
              fname);
      return 0;
    }
    relax->vid2 =
      ts_get_index(relax->ts, fv_get_varname(fname, "v", buf));
    if (relax->vid2 < 0) {
      hd_quit("vel_relax_init: The file '%s' does not contain 'v'\n",
              fname);
      return 0;
    }

    schedSetPrivateData(event, relax);
    return 1;
  }

  return 0;
}


/* Read the relaxation data in. */
double vel_relax_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  vel_relax_data_t *relax = (vel_relax_data_t *)schedGetPrivateData(event);

  /* Check whether the current time exceeds the next update time */
  if (t >= (relax->tnext - SEPS)) {
    relax_info_t *rlx = master->vel_rlx;
    int cc, cs, c = 0;
    double rr, u, v;
    
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      u = v = 0.0;
      u = ts_eval_xyz(relax->ts, relax->vid1, t,
		      geom->u1x[cs], geom->u1y[cs], 
		      geom->cellz[c] * master->Ds[cs]);
      v = ts_eval_xyz(relax->ts, relax->vid2, t,
		      geom->u1x[cs], geom->u1y[cs],
		      geom->cellz[c] * master->Ds[cs]);
      rlx->val1[c] = cos(geom->thetau1[cs]) * u + sin(geom->thetau1[cs]) * v;
      u = v = 0.0;
      u = ts_eval_xyz(relax->ts, relax->vid1, t,
		      geom->u2x[cs], geom->u2y[cs], 
		      geom->cellz[c] * master->Ds[cs]);
      v = ts_eval_xyz(relax->ts, relax->vid2, t,
		      geom->u2x[cs], geom->u2y[cs],
		      geom->cellz[c] * master->Ds[cs]);
      rlx->val2[c] = cos(geom->thetau2[cs]) * v - sin(geom->thetau2[cs]) * u;
    }
    /* Get the relaxation time constant */
    if (rlx->tctype & RLX_CONS)
      rlx->rate = relax->rate;
    else if (rlx->tctype & RLX_FILE) {
    /* Relaxation rate read from file. Temporally variable only. */
      char buf[MAXSTRLEN];
      rr = ts_eval(relax->tc, relax->tcid, master->t);
      sprintf(buf, "%f %s\n",rr, relax->tc->varunit[relax->tcid]);
      tm_scale_to_secs(buf, &rr);
      rlx->rate = rr;
    }
    
    /* Set the next update time */
    relax->tnext += relax->dt;
  }

  return relax->tnext;
}


/* Cleanup routines */
void vel_relax_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  vel_relax_data_t *relax = (vel_relax_data_t *)schedGetPrivateData(event);
  relax_info_t *rlx = master->vel_rlx;

  hd_ts_free(master, relax->ts);
  if (relax->tc != NULL)
    hd_ts_free(master, relax->tc);
  if (rlx->size && rlx->val1) d_free_1d(rlx->val1);
  if (rlx->size && rlx->val2) d_free_1d(rlx->val2);
  free(relax);
}


/*-------------------------------------------------------------------*/
/* Relax velocity towards specified values in the window.            */
/* The change in elevation is calculated from the vel_rlx passed     */
/* from the master every scheduled timestep, and applied to u1/u2    */
/* on the 3D step.                                                   */
/*-------------------------------------------------------------------*/
void do_vel_relax(geometry_t *window,  /* Window geometry            */
		  window_t *windat,    /* Window data                */		  
		  win_priv_t *wincon,  /* Window constants           */
		  int mode             /* 2D / 3D mode flag          */
		  )
{
  relax_info_t *rlx = windat->vel_rlx;
  double *depth;
  double rr;
  int bn, c, cs, cc;

  if (wincon->velrlx & NONE) return;

  if (mode == VEL3D) {
    double *u1 = windat->u1;
    double *u2 = windat->u2;
    if (wincon->tide_r == MEAN_R) {
      if (windat->u1m) u1 = windat->u1m;
      if (windat->u2m) u2 = windat->u2m;
    }
    for (cc = 1; cc <= wincon->vc1; cc++) {
      c = wincon->s1[cc];
      rr = get_relax_rate(rlx, u1[c], u2[c], rlx->val1[c], rlx->val2[c]);
      windat->u1[c] -= windat->dt2d * (u1[c] - rlx->val1[c]) / rr;
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	rr = get_relax_rate(rlx, u1[c], u2[c], rlx->val1[c], rlx->val2[c]);
	windat->u1[c] -= windat->dt2d * (u1[c] - rlx->val1[c]) / rr;
      }
    }
    for (cc = 1; cc <= wincon->vc2; cc++) {
      c = wincon->s2[cc];
      rr = get_relax_rate(rlx, u1[c], u2[c], rlx->val1[c], rlx->val2[c]);
      windat->u2[c] -= windat->dt2d * (u2[c] - rlx->val2[c]) / rr;
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	rr = get_relax_rate(rlx, u1[c], u2[c], rlx->val1[c], rlx->val2[c]);
	windat->u2[c] -= windat->dt2d * (u2[c] - rlx->val2[c]) / rr;
      }
    }
  } else {
    double *v1 = wincon->d1;
    double *v2 = wincon->d2;
    double *u1av = windat->u1av;
    double *u2av = windat->u2av;
    if (wincon->tide_r == MEAN_R ) {
      if (windat->u1am) u1av = windat->u1am;
      if (windat->u2am) u2av = windat->u2am;
    }
    depth = wincon->d5;           /* Set in precalc_u1()               */
    memset(v1, 0, window->sgsizS * sizeof(double));
    memset(v2, 0, window->sgsizS * sizeof(double));

    /* Wet cells to process.                                           */
    for (cc = 1; cc <= wincon->vc1; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      v1[cs] += rlx->val1[c] * windat->dzu1[c] * wincon->mdx[cs];
    }
    /* Open boundary normal velocities. */
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	cs = window->m2d[c];
	v1[cs] += rlx->val1[c] * windat->dzu1[c] * wincon->mdx[cs];
      }
    }
    /* Calculate the depth averaged relaxation.                        */
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      v1[cs] /= depth[cs];
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no2_e1; cc++) {
	cs = open->obc_e1[cc];
	v1[cs] /= (depth[cs] - window->botzu1[cs]);
      }
    }

    depth = wincon->d6;
    memset(v2, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= wincon->vc2; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      v2[cs] += rlx->val2[c] * windat->dzu2[c] * wincon->mdy[cs];
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	cs = window->m2d[c];
	v2[cs] += rlx->val2[c] * windat->dzu2[c] * wincon->mdy[cs];
      }
    }
    /* Calculate the depth averaged relaxation.                        */
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      v2[cs] /= depth[cs];
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no2_e2; cc++) {
	cs = open->obc_e2[cc];
	v2[cs] /= (depth[cs] - window->botzu2[cs]);
      }
    }

    /* Adjust the 2D velocities.                                       */
    for (cc = 1; cc <= wincon->vc1; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      rr = get_relax_rate(rlx, u1av[cs], u2av[cs], v1[cs], v2[cs]);
      windat->u1av[cs] -= windat->dt2d * (u1av[cs] - v1[cs]) / rr;
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	cs = window->m2d[c];
	rr = get_relax_rate(rlx, u1av[cs], u2av[cs], v1[cs], v2[cs]);
	windat->u1av[cs] -= windat->dt2d * (u1av[cs] - v1[cs]) / rr;
      }
    }
    for (cc = 1; cc <= wincon->vc2; cc++) {
      c = wincon->s2[cc];
      cs = window->m2d[c];
      rr = get_relax_rate(rlx, u1av[cs], u2av[cs], v1[cs], v2[cs]);
      windat->u2av[cs] -= windat->dt2d * (u2av[cs] - v2[cs]) / rr;
    }
    for (bn = 0; bn < window->nobc; bn++) {
      open_bdrys_t *open = window->open[bn];
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	cs = window->m2d[c];
	rr = get_relax_rate(rlx, u1av[cs], u2av[cs], v1[cs], v2[cs]);
	windat->u2av[cs] -= windat->dt2d * (u2av[cs] - v2[cs]) / rr;
      }
    }
  }
}

/* END do_vel_relax_()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the relaxation time constant                       */
/*-------------------------------------------------------------------*/
double get_relax_rate(relax_info_t *rlx, double u1, double u2, double v1, double v2)
{
  double dv;

  if (rlx->tctype & (RLX_CONS|RLX_FILE))
    return(rlx->rate);
  else if (rlx->tctype & RLX_LINR) {
    dv = sqrt(u1 * u1 + u2 * u2) - sqrt(v1 * v1 + v2 * v2);
    return((fabs(dv) - rlx->dv0) * rlx->slope + rlx->tc0);
  } else if (rlx->tctype & RLX_EXP) {
    dv = sqrt(u1 * u1 + u2 * u2) - sqrt(v1 * v1 + v2 * v2);
    return(exp(rlx->dv0 / fabs(dv)));
  }
  return 0.0;
}

/* END get_relax_rate()                                              */
/*-------------------------------------------------------------------*/
