/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/eta_relax.c
 *  
 *  Description:
 *  Relaxation of surface elevation. At present,
 *  this code is simple minded, and probably should
 *  not be used in an application which includes
 *  tidal forcing. Tracer concentrations are not
 *  considered here - whether this code affects
 *  them or not depends on when it is called during
 *  the model time step.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: eta_relax.c 6440 2019-12-04 03:00:23Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

typedef struct {
  timeseries_t *ts;             /* Time-Series file */
  int vid;                      /* df_variable_t id into TS file. */
  double dt;                    /* Relaxation time step */
  double rate;                  /* Relaxation rate */
  timeseries_t *tc;             /* Relaxation rate timeseries */
  int tcid;                     /* Relaxation rate id */
} eta_relax_data_t;

int eta_relax_init(sched_event_t *event)
{
  char buf[MAXSTRLEN];
  char fname[MAXSTRLEN];
  double dt, tconst;
  eta_relax_data_t *relax = NULL;
  master_t *master = (master_t *)schedGetPublicData(event);
  char tc_name[] = "eta_relaxation_time_constant";
  relax_info_t *rlx = master->eta_rlx;

  prm_set_errfn(hd_silent_warn);

  if (master->etarlx & (RELAX|ALERT|BOUNDARY|ETA_TPXO)) {

    /* Set the input time step, and time constant */
    strcpy(fname, rlx->rlxn);
    dt = rlx->rlxdt;

    /* Allocate memory for relaxation structure */
    relax = (eta_relax_data_t *)malloc(sizeof(eta_relax_data_t));
    memset(relax, 0, sizeof(eta_relax_data_t));

    /* Assign the time constant                                */
    relax->tc = NULL;
    if (strlen(rlx->rlxtc)) {
      char tu0[MAXSTRLEN], tu1[MAXSTRLEN];
      double r0, r1, d0, d1, s1, s2;
      if (sscanf(rlx->rlxtc, "%lf %s", &tconst, buf) == 2) {
	tm_scale_to_secs(rlx->rlxtc, &tconst);
	rlx->tctype = RLX_CONS;
      } else if (sscanf(rlx->rlxtc, "%s %lf %lf %s %lf %lf %s", 
			buf, &d0, &r0, tu0, &d1, &r1, tu1) == 7) {
	if(strcmp(buf, "linear") == 0) {
	  hd_warn("Linear adaptive relaxation performed on elevation.\n");
	  if(d1 > d0) {
	    s1 = d1; s2 = r1;
	    d1 = d0; r1 = r0;
	    d0 = s1; r0 = s2;
	  }
	  tconst = 1.0;
	  rlx->dv0 = d0;
	  rlx->dv1 = d1;
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  rlx->tc0 = r0;
	  rlx->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	  rlx->tctype = (RLX_ADPT|RLX_LINR);
	} else if(strcmp(buf, "temporal") == 0) {
	  hd_warn("Temporal linear relaxation performed on elevation.\n");
	  tconst = 1.0;
	  rlx->dv0 = d0;
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  rlx->tc0 = r0;
	  rlx->tc1 = r1;
	  rlx->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	  rlx->tctype = (RLX_TIME);
	} else if(strcmp(buf, "depth") == 0) {
	  /* Depth scaled linear relaxation */
	  hd_warn("Depth scaled linear relaxation performed on elevation\n");
	  tconst = 1.0;
	  rlx->dv0 = d0;
	  rlx->dv1 = d1;
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  rlx->tc0 = r0;
	  rlx->tc1 = r1;
	  rlx->slope = (d1 - d0) ? (r1 - r0) / (d1 - d0) : 0.0;
	  rlx->tctype = (RLX_ADPT|RLX_DEP);
	} else if(strcmp(buf, "exp_depth") == 0) {
	  /* Depth scaled exponential relaxation                    */
	  /* rate = (r0 - r1)exp(-depth/d0) + (r1 - r0.exp(-d1/d0)) */
	  /* d0 = scaling factor (e.g. 500)                         */
	  /* r0 = rate at depth equal 0                             */
	  /* d1 = maximum depth                                     */
	  /* r1 = rate at depth d1                                  */
	  /* e.g. exp_depth 500 100 days 4000 7 days                */
	  /* d0 ~ 40 for d1 = 200                                   */
	  /* d0 ~ 500 for d1 = 4000                                 */
	  hd_warn("Depth scaled exponential relaxation performed on elevation\n");
	  tconst = 1.0;
	  rlx->dv0 = fabs(d0);
	  rlx->dv1 = fabs(d1);
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  rlx->tc0 = r0;
	  rlx->tc1 = r1;
	  rlx->slope = r0 * exp(-d1 / d0);
	  rlx->tctype = (RLX_ADPT|RLX_EDEP);
	} else if(strcmp(buf, "cos_depth") == 0) {
	  double d, r;
	  /* Depth scaled cosine relaxation                                */
	  /* rate = 0.5*((r0-r1)*cos(d*pi/(d1-d0)-d0*PI/(d1-d0))+(t0+t1)); */
	  /* This is designed for d0 < d1                                  */
	  hd_warn("Depth scaled cosine relaxation performed on elevation\n");
	  tconst = 1.0;
	  d0 = fabs(d0);
	  d1 = fabs(d1);
	  sprintf(buf, "%f %s", r0, tu0);
	  tm_scale_to_secs(buf, &r0);
	  sprintf(buf, "%f %s", r1, tu1);
	  tm_scale_to_secs(buf, &r1);
	  if (d0 > d1) {
	    d = d0;
	    d0 = d1;
	    d1 = d;
	    r = r0;
	    r0 = r1;
	    r1 = r;
	  }
	  rlx->dv0 = fabs(d0);
	  rlx->dv1 = fabs(d1);
	  rlx->tc0 = r0;
	  rlx->tc1 = r1;
	  rlx->slope = d0 * PI / (d1 - d0);
	  rlx->tctype = (RLX_ADPT|RLX_CDEP);
	} else
	  hd_quit("Linear adaptive relaxation format, e.g; 'linear 10 1 day 0.1 10 days'\n");
      } else if (sscanf(rlx->rlxtc, "%s %lf %lf %s", buf, &d0, &r0, tu0) == 4) {
	if (strcmp(buf, "exponential") == 0) {
	  hd_warn("Exponential adaptive relaxation performed on elevation.\n");
	  tconst = 1.0;
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
	rlx->tctype = RLX_FILE;
      }
    }

    if (master->etarlx & RELAX && tconst <= 0)
      hd_quit
        ("eta_relax_init: Invalid surface elevation relaxation time constant specified.\n");

    /* Populate and register the scheduler events. */
    relax->dt = dt;
    if (tconst)relax->rate = rlx->rate = tconst;
    if (master->etarlx & ETA_TPXO) {
      int tn;
      if ((tn = tracer_find_index("tpxotide", master->ntrS, master->trinfo_2d)) >= 0)
	rlx->val1 = master->tr_wcS[tn];
      else
	hd_quit("eta_relax_init: Can't find tpxotide: use NUMBERS TPXO.\n");
    } else {
      relax->ts = hd_ts_read(master, fname, 1);
      relax->vid =
	ts_get_index(relax->ts, fv_get_varname(fname, "eta", buf));
      if (relax->vid < 0) {
	hd_quit("eta_relax_init: The file '%s' does not contain 'eta'\n",
		fname);
      }
      return 0;
    }

    schedSetPrivateData(event, relax);
    return 1;
  }

  return 0;
}


/* Read the relaxation data in. */
double eta_relax_event(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  geometry_t *geom = master->geom;
  eta_relax_data_t *relax = (eta_relax_data_t *)schedGetPrivateData(event);
  relax_info_t *rlx = master->eta_rlx;

  /*
   * this causes a compilter warning but not needed anyway as the init 
   * has a quit for this case
   */
  /* if (relax->vid < 0) return; */

  /* Check whether the current time exceeds the next update time */
  if (t >= (event->next_event - SEPS)) {
    int cc, c = 0;
    if (!(master->etarlx & ETA_TPXO)) {
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	rlx->val1[c] = ts_eval_xy(relax->ts, relax->vid, t,
				  geom->cellx[c], geom->celly[c]);
      }
    }

    /* Get the relaxation time constant */
    if (rlx->tctype & RLX_CONS)
      rlx->rate = relax->rate;
    else if (rlx->tctype & RLX_FILE) {
    /* Relaxation rate read from file. Temporally variable only. */
      char buf[MAXSTRLEN];
      double rr;
      rr = ts_eval(relax->tc, relax->tcid, master->t);
      sprintf(buf, "%f %s\n",rr, relax->tc->varunit[relax->tcid]);
      tm_scale_to_secs(buf, &rr);
      rlx->rate = rr;
    }

    /* Set the next update time */
    event->next_event += relax->dt;
  }

  return event->next_event;
}


/* Cleanup routines */
void eta_relax_cleanup(sched_event_t *event, double t)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  eta_relax_data_t *relax = (eta_relax_data_t *)schedGetPrivateData(event);
  relax_info_t *rlx = master->eta_rlx;

  hd_ts_free(master, relax->ts);
  if (rlx->size && rlx->val1) d_free_1d(rlx->val1);
  free(relax);
}


/*-------------------------------------------------------------------*/
/* Relax surface elevation towards specified values in the window.   */
/* The change in elevation is calculated from the eta_rlx passed     */
/* from the master every 3D timestep, and this is converted to an    */
/* elevation (m) and applied to neweta on the 2D step. Note that     */
/* fluxes for the 3D mode are integrated over the 2D step and        */
/* applied to the tracer fluxes at the surface.                      */
/*-------------------------------------------------------------------*/
void do_eta_relax(geometry_t *window,  /* Window geometry            */
		  int c,               /* Coordinate to relax        */
		  double relax_rate,   /* Relaxation rate            */
		  int tidef            /* Flag to remove tide        */
		  )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  relax_info_t *rlx = windat->eta_rlx;
  double *eta = wincon->neweta;
  double tide = 0.0;
  double colflux;

  if (isnan(rlx->val1[c])) return;

  /* Get an estimate of tidal height via csr tide model              */
  if (tidef == CSR_R && wincon->tc.nt) {
    double gmt = windat->t - wincon->tz;
    double ramp = (wincon->rampf & TIDALH) ? windat->rampval : 1.0;
    tide = ramp * csr_tide_eval(&wincon->tc, window->c2cc[c], gmt);
  }
  /* Use the mean elevation to remove relaxation elevation           */
  if (tidef == MEAN_R && windat->etam) {
    if (!windat->meanc[c]) 
      return;
    else
      eta = windat->etam;
  }

  /* Get the difference between eta and the relaxation value         */
  colflux = (eta[c] - tide - rlx->val1[c]);
  /* Adaptive relaxation                                             */
  if (rlx->tctype & RLX_LINR) {
    double d = max(min(fabs(colflux), rlx->dv0), rlx->dv1);
    relax_rate = (d - rlx->dv0) * rlx->slope + rlx->tc0;
  }
  if (rlx->tctype & RLX_EXP) {
    double d = 1000.0 * 86400.0;
    relax_rate = min(exp(rlx->dv0 / fabs(colflux)), d);
  }
  if (rlx->tctype & RLX_TIME)
    relax_rate = min((windat->days - rlx->dv0) * rlx->slope + rlx->tc0, rlx->tc1);
  if (rlx->tctype & RLX_DEP) {
    double depth = windat->eta[c] - window->botz[c];
    double dmin = min(rlx->dv0, rlx->dv1);
    double dmax = max(rlx->dv0, rlx->dv1);
    depth = min(max(depth, dmin), dmax);
    relax_rate = (depth - rlx->dv0) * rlx->slope + rlx->tc0;
  } else if (rlx->tctype & RLX_EDEP) {
    double depth = windat->eta[c] - window->botz[c];
    double rmin = min(rlx->tc0, rlx->tc1);
    double rmax = max(rlx->tc0, rlx->tc1);
    relax_rate = max(rmin,
		     (rmax - rmin) * exp(-fabs(depth) / rlx->dv0) + (rmin - rlx->slope));
  } else if (rlx->tctype & RLX_CDEP) {
    double depth = windat->eta[c] - window->botz[c];
    double r0 = rlx->tc0;
    double r1 = rlx->tc1;
    double d = fabs(depth);
    double dd = rlx->dv1 - rlx->dv0;
    if (d < rlx->dv0)
      relax_rate = r0;
    else if (d > rlx->dv1)
      relax_rate = r1;
    else
      relax_rate = 0.5 * ((r0 - r1) * cos(d * PI / dd - rlx->slope) + (r0 + r1));
  }
  if (windat->eta_tc) windat->eta_tc[c] = relax_rate / 86400.0;

  /* Get the elevation to add to the predicted elevation to converge */
  /* to the relaxation elevation.                                    */
  colflux = windat->dt2d * colflux / relax_rate;
  if (windat->eta_inc) windat->eta_inc[c] = colflux;

  /* Calculate the new elevation                                     */
  wincon->neweta[c] = max(wincon->neweta[c] - colflux, window->botz[c]);

  /* Rate of change of eta                                           */
  windat->detadt[c] = (wincon->neweta[c] - windat->etab[c]) / windat->dt2d;

  /* Integrate the volume (m3) for the 3D mode (on odd substeps of   */
  /* the 2D mode as is done with 2D momentum fluxes).                */
  if (((wincon->ic + 1) % 2) == 0)
    wincon->eta_rlx3d[c] += colflux * window->cellarea[c];

}

/* END do_eta_relax_()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the tracer increment to elevation. The elevation   */
/* increment is spread over the 2D step.                             */
/*-------------------------------------------------------------------*/
void do_eta_increment(geometry_t *window   /* Window geometry        */
		      )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double einc;

  int n, c, cc;

  for (n = 0; n < windat->ntrS; n++) {
    if (wincon->trinfo_2d[n].increment & ETA_M && windat->trincS[n]) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	einc = windat->tr_wcS[n][c] / (double)windat->iratio;
	wincon->neweta[c] += einc;
	/* Rate of change of eta                                     */
	windat->detadt[c] = (wincon->neweta[c] - windat->etab[c]) / 
	  windat->dt2d;
	/* Integrate the volume (m3) for the 3D mode (on odd         */
	/* substeps of the 2D mode as is done with 2D momentum       */
	/* fluxes).                                                  */
	if (((wincon->ic + 1) % 2) == 0)
	  wincon->eta_rlx3d[c] -= (einc * window->cellarea[c]);
      }
      if (wincon->ic + 1 == windat->iratio) windat->trincS[n] = 0;      
    }
  }
}

/* END do_eta_increment()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the tracer increment to velocity. The velocity     */
/* increment is added at the end of the 3D step.                     */
/*-------------------------------------------------------------------*/
void do_vel_increment_3d(geometry_t *window   /* Window geometry     */
		      )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */

  int n, c, cc;

  for (n = 0; n < windat->ntr; n++) {
    if (wincon->trinfo_3d[n].increment & U1VEL && windat->trinc[n]) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	windat->u1[c] += windat->tr_wcS[n][c];
      }
      windat->trinc[n] = 0;
    }
    if (wincon->trinfo_3d[n].increment & U2VEL && windat->trinc[n]) {
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	windat->u2[c] += windat->tr_wcS[n][c];
      }
      windat->trinc[n] = 0;
    }
  }
}

/* END do_vel_increment_3d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the tracer increment to 2D velocity. The velocity  */
/* increment is spread over the 2D step.                             */
/*-------------------------------------------------------------------*/
void do_vel_increment_2d(geometry_t *window   /* Window geometry     */
		      )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  open_bdrys_t **open = window->open;   /* Window OBC structure      */
  double *vinc = wincon->d1;
  double *depth;
  int n, c, cs, cc, bn;

  for (n = 0; n < windat->ntr; n++) {
    if (wincon->trinfo_3d[n].increment & U1VEL && windat->trinc[n]) {
      depth = wincon->d5;      /* Set in precalc_u1()                */
      memset(vinc, 0, window->sgsizS * sizeof(double));
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	vinc[cs] += windat->tr_wcS[n][c] * windat->dzu1[c] * 
	  wincon->mdx[cs];
      }
      for (cc = 1; cc <= wincon->vcs1; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	windat->u1av[cs] += (vinc[cs] / 
			     ((double)windat->iratio * depth[cs]));
      }
      if(!wincon->dobdry_u1) {
	for (bn = 0; bn < window->nobc; bn++) {
	  for (cc = 1; cc <= open[bn]->no2_e1; cc++) {
	    cs = open[bn]->obc_e1[cc];
	    windat->u1av[cs] += (vinc[cs] / ((double)windat->iratio * 
				      (depth[cs] - window->botzu1[cs])));
	  }
	}
      }
    }
    if (wincon->trinfo_3d[n].increment & U2VEL && windat->trinc[n]) {
      depth = wincon->d6;      /* Set in precalc_u2()                */
      memset(vinc, 0, window->sgsizS * sizeof(double));
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	vinc[cs] += windat->tr_wcS[n][c] * windat->dzu2[c] * 
	  wincon->mdy[cs];
      }
      for (cc = 1; cc <= wincon->vcs1; cc++) {
	c = wincon->s1[cc];
	cs = window->m2d[c];
	windat->u2av[cs] += (vinc[cs] / 
			     ((double)windat->iratio * depth[cs]));
      }
      if(!wincon->dobdry_u1) {
	for (bn = 0; bn < window->nobc; bn++) {
	  for (cc = 1; cc <= open[bn]->no2_e2; cc++) {
	    cs = open[bn]->obc_e2[cc];
	    windat->u2av[cs] += (vinc[cs] / ((double)windat->iratio * 
				      (depth[cs] - window->botzu2[cs])));
	  }
	}
      }
    }
  }
}

/* END do_vel_increment_2d()                                         */
/*-------------------------------------------------------------------*/
