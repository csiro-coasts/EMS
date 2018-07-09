/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/inputs/sourcesink.c
 *  
 *  Description:
 *  Routines dealing with sources and sinks
 *  of water and tracers in the meco model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sourcesink.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hd.h>
#include <tracer.h>

static int parse_ss_location(pss_t *p, char *line);
static void pss_locate(pss_t *p, double t);

void ref_depth(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	       pss_t *pss, double *zlow, double *zhigh, int cs);
int hd_pss_vec_m(void *data, pss_t *pss);
int hd_pss_vec_w(void *data, pss_t *pss);
int pss_vec_m(master_t *master, pss_t *pss);
int pss_vec_w(geometry_t *window, pss_t *pss);
int hd_pss_region_m(void *data, char *dname, pss_t *pss);
int pss_init_region_m(master_t *master, char *dname, pss_t *pss);
int hd_pss_region_w(void *data, char *dname, pss_t *pss);
int pss_init_region_w(geometry_t *window, char *dname, pss_t *pss);
void hd_pss_read(char *name, char *t_units, pss_t **pss,
		 int *np, void *data,
		 int (*xyzijk) (void *, double, double, double, int *, int *,
				int *), int (*trI) (void *, char *), int (*pss_vec) (void *, pss_t *),
		 int (*pss_region) (void *, char *, pss_t *));

/*-------------------------------------------------------------------*/
/* Initialisation routine for sources and sinks - called once at     */
/* model start-up time to read source/sink lists, and surface heat   */
/* flux parameters.                                                  */
/*-------------------------------------------------------------------*/
void sourcesink_init(parameters_t *params,  /* Input parameters      */
                     master_t *master,      /* Master data           */
                     geometry_t **window,   /* Window geometry       */
                     window_t **windat,     /* Window data           */
                     win_priv_t **wincon    /* Window constants      */
  )
{
  geometry_t *geom = master->geom;          /* Global geometry       */
  int n, s, m = 0;
  int t;

  prm_set_errfn(hd_silent_warn);

  /* Initialise point source/sinks on the master                     */
  hd_pss_read(params->prmname, master->timeunit, &master->pss,
           &master->npss, master, hd_xyztoindex_m, hd_get_tracer_index_m, 
	      hd_pss_vec_m, hd_pss_region_m);
  
  for (s = 0; s < master->npss; s++) {
    pss_t *pss = &master->pss[s];

    for (t = 0; t < pss->ntsfiles; t++) {
      hd_ts_check(master, &pss->ts[t]);
    }
    if (master->pss[s].loc)
      hd_ts_check(master, pss->loc);
    if (pss->e1[0] == NOTVALID && pss->e1[0] ==NOTVALID && pss->e1[0] ==NOTVALID)
      hd_quit
        ("sourcesink_init: %s is in solid or outside cell; (i,j,k) = (%d,%d,%d).\n",
         pss->name, geom->s2i[pss->e1[0]], geom->s2j[pss->e1[0]],
         geom->s2k[pss->e1[0]]);
  }

  if (master->npss) {
    master->wflux = d_alloc_1d(master->npss);
    master->tflux = d_alloc_2d(master->npss, master->ntr);
    memset(master->wflux, 0, master->npss * sizeof(double));
    /* If flow is used in the point sources, then this must be saved */
    /* if transport files are written, and read if FFSL is used in   */
    /* the transport mode. The variable w is used to store/read the  */
    /* waterss array, hence waterss must be copied to wmean for      */
    /* transport dumps and in transport mode copied back to waterss. */
    for (s = 0; s < master->npss; s++) {
      pss_t *p = &master->pss[s];
      if (p->watertsid >= 0) {
	/* If transport files are written and point sources contain  */
	/* flow, then set a flag to copy waterss to wmean.           */
	if (master->means & VOLFLUX) {
	  master->means |= PSSFLUX;
	  for (n = 1; n <= master->nwindows; n++) 
	    wincon[n]->means |= PSSFLUX;
	}
	/* If transport mode is SP_FFSL and point sources contain    */
	/* flow, then set a flag to copy w to waterss and set the    */
	/* conservation flag to re-compute w.                        */
	if (master->tmode & SP_FFSL) {
	  if (master->conserve == NONE)
	    master->conserve = (CONS_PSS|CONS_W);
	  else {
	    master->conserve |= (CONS_PSS|CONS_W);
	    if (master->conserve & CONS_WS) master->conserve &= ~CONS_WS;
	  }
	  for (n = 1; n <= master->nwindows; n++) 
	    wincon[n]->conserve = master->conserve;
	}
      }
    }
  }
  /* Initialise point source/sinks in each window                    */
  t = 1;
  for (n = 1; n <= master->nwindows; n++) {
    hd_pss_read(params->prmname, master->timeunit, &wincon[n]->pss,
             &wincon[n]->npss, window[n],
             hd_xyztoindex_w, hd_get_tracer_index_w, hd_pss_vec_w,
	     hd_pss_region_w);

    if (wincon[n]->npss) {
      m += wincon[n]->npss;
      for (s = 0; s < wincon[n]->npss; s++) {
	pss_t *p = &wincon[n]->pss[s];
	hd_warn("Sourcesink%d (%s) located in window%d : (%d %d)\n",
		s, wincon[n]->pss[s].name, n,
		window[n]->s2i[wincon[n]->pss[s].e1[0]],
		window[n]->s2j[wincon[n]->pss[s].e1[0]]);
	if (p->vc > 1) t = 0;
      }
      windat[n]->wflux = d_alloc_1d(wincon[n]->npss);
      windat[n]->tflux = d_alloc_2d(wincon[n]->npss, windat[n]->ntr);
    }
  }
  if (t && m != master->npss)
    hd_quit("Mismatch between point sourcesinks on master(%d) and slave(%d). Check the runlog.\n", master->npss, m);
  /* Find the master - salve pss map                                 */
  for (n = 1; n <= master->nwindows; n++) {
    for (s = 0; s < wincon[n]->npss; s++) {
      pss_t *p = &wincon[n]->pss[s];
      int m;
      for (m = 0; m < master->npss; m++) {
	pss_t *pss = &master->pss[m];
	if (strcmp(pss->name, p->name) == 0)
	  p->s2m = m;
      }
    }
  }
}

/* END sourcesink_init()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to evaluate additional water fluxes due to sources or     */
/* sinks in the window in the master.                                */
/*-------------------------------------------------------------------*/
void sourcesink(master_t *master  /* Master data                     */
  )
{
  int s, n, nn;

  /* Loop over the pointsourcesinks                                  */
  for (s = 0; s < master->npss; s++) {
    pss_t *p = &master->pss[s];
    double f = 1.0 / (double)p->vc;
    if (p->watertsid >= 0)
      master->wflux[s] = 
	ts_multifile_eval(p->ntsfiles, p->ts, p->watertsid, master->t) * f;
    /* Get tracer index for this sourcesink                          */
    for (n = 0; n < master->ntr; n++) {
      int ts_id = -1;
      master->tflux[n][s] = NOTVALID;
      for (nn = 0; nn < p->nv && ts_id < 0; nn++) {
        if (p->tmap[nn] == n) {
          ts_id = nn;
	  break;
	}
      }

      if (ts_id < 0)
        continue;

      /* Calculate flux of tracer to add (kg/s)                      */
      if (p->watertsid >= 0) {
        /* In this case, we have already calculated wflux above, so  */
        /* only need to calculate tracer concentration in source.    */
        double conc = ts_multifile_eval(p->ntsfiles, p->ts, ts_id, master->t);
	                                                    /* kg/m3 */
        master->tflux[n][s] = master->wflux[s] * conc;
      } else {
        /* No water associated with this source, so tracer mass flux */
        /* flux is directly specified by input time series.          */
        master->tflux[n][s] = 
	  ts_multifile_eval(p->ntsfiles, p->ts, ts_id, master->t) * f;
      }
    }
  }

  /* Clear the 2-d and 3-d water source/sink arrays                  */
  memset(master->waterss2d, 0, geom->szcS);
  memset(master->waterss, 0, geom->szc);
}

/* END sourcesink()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to evaluate additional water fluxes due to sources or     */
/* sinks in the window.                                              */
/*-------------------------------------------------------------------*/
void ss_water(geometry_t *window, /* Window geometry                 */
              window_t *windat,   /* Window data                     */
              win_priv_t *wincon  /* Window constants                */
  )
{
  int s, cc, c, nc;
  double wflux;

  /* Clear the 2-d and 3-d water source/sink arrays                  */
  memset(windat->waterss, 0, window->szc * sizeof(double));
  memset(windat->waterss2d, 0, window->szcS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* For the transport mode using FFSL, if point source volume       */
  /* fluxes are included in the transport file (saved to wmean) then */
  /* copy these into waterss and reverse engineer the original flow  */
  /* rates and tracer fluxes.                                        */
  if (master->trasc == FFSL && master->conserve & CONS_PSS && 
      !(master->conserve & CONS_NOF)) {
    int n, nn, cs;
    /* Only comput fluxes when w is read from file (in reset.c). For */
    /* sub-stepping using TRATIO < 1.0, only compute on the first    */
    /* step.                                                         */
    if(windat->w[0] == 0.0) return;

    /* Initialise                                                    */
    memcpy(windat->waterss, windat->w, window->sgsiz * sizeof(double));
    memset(windat->waterss2d, 0, window->szcS * sizeof(double));
    /* Get the vertical integral of the inflow                       */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      windat->waterss2d[cs] += windat->waterss[c] * window->cellarea[cs];
    }
    /* Reset the inflow and tracer flux                              */
    for (s = 0; s < wincon->npss; s++) {
      pss_t *p = &wincon->pss[s];
      if (p->watertsid >= 0) {
	/* Get the tracer sourcesink concentration as read from file */
	for (n = 0; n < windat->ntr; n++) {
	  int ts_id = -1;
	  for (nn = 0; nn < p->nv && ts_id < 0; nn++) {
	    if (p->tmap[nn] == n) {
	      ts_id = nn;
	      break;
	    }
	  }
	  if (ts_id < 0)
	    continue;
	  if (windat->wflux[s]) 
	    windat->tflux[n][s] /= windat->wflux[s];
	  else
	    windat->tflux[n][s] = 0.0;
	}

	/* Reset the volume fluxes                                     */
	for (nc = 0; nc < p->vc; nc++) {
	  /* Find position of source/sink horizontally and vertically, */
	  /* and truncate to available water column.                   */
	  pss_locate(p, windat->t);
	  c = p->e1[nc];
	  cs = window->m2d[c];
	  windat->wflux[s] = windat->waterss2d[cs];
	}
	/* Reset the tracer fluxes using saved wflux and above tflux.  */
	/* Not used : we now use waterss for volume fluxes directly.   
	for (n = 0; n < windat->ntr; n++) {
	  int ts_id = -1;
	  for (nn = 0; nn < p->nv && ts_id < 0; nn++)
	    if (p->tmap[nn] == n) {
	      ts_id = nn;
	      break;
	    }
	  if (ts_id < 0)
	    continue;
	  windat->tflux[n][s] *= windat->wflux[s];
	}
	*/
      }
    }
    windat->w[0] = 0.0;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Loop over the pointsourcesinks                                  */
  for (s = 0; s < wincon->npss; s++) {
    pss_t *p = &wincon->pss[s];
    /* Check that this pointsourcesink has water                     */
    if (p->watertsid >= 0) {
      int c, cb, cs, c2, zp1;
      double pdz;
      double zhigh_i;
      double zlow_i;

      /* Flow of water                                               */
      wflux = windat->wflux[s];

      /* Loop over horizontal locations                              */
      for (nc = 0; nc < p->vc; nc++) {
	/* Find position of source/sink horizontally and vertically, */
	/* and truncate to available water column.                   */
	pss_locate(p, windat->t);
	zhigh_i = p->zhigh;
	zlow_i  = p->zlow;
	c = p->e1[nc];
	cs = p->e2[nc];
	cb = p->e3[nc];
	ref_depth(window, windat, wincon, p, &zlow_i, &zhigh_i, cs);
	c2 = window->m2d[c];
	zp1 = window->zp1[c];

	/* Not in the interior of the model, so ignore it!           */
	if (c == 0)
	  continue;
	
	zlow_i = max(zlow_i, window->botz[c2]);
	zhigh_i = min(zhigh_i, windat->topz[c2]);
	pdz = zhigh_i - zlow_i;

	/* Add the water, either distributed vertically or into one  */
	/* layer.                                                    */
	if (pdz > 0.0) {
	  /* Distributed - loop over the affected layers             */
	  c = c2;
	  while (c != window->zm1[c]) {
	    double ctop = (c == cs) ? windat->topz[c2] : window->gridz[zp1];
	    double cbot = (c == cb) ? window->botz[c2] : window->gridz[c];
	    double zlow;
	    double zhigh;
	    double frac;

	    /* Range of source/sink in this cell                     */
	    if (zlow_i <= cbot)
	      zlow = cbot;
	    else if (zlow_i >= ctop)
	      zlow = ctop;
	    else
	      zlow = zlow_i;
	    if (zhigh_i <= cbot)
	      zhigh = cbot;
	    else if (zhigh_i >= ctop)
	      zhigh = ctop;
	    else
	      zhigh = zhigh_i;
	    
	    /* Fraction of flow in this cell                         */
	    frac = (zhigh - zlow) / pdz;
	    windat->waterss[c] += wflux * frac;
	    zp1 = c;
	    c = window->zm1[c];
	  }
	} else {
	  /* Source/sink is not distributed vertically               */
	  c = p->e1[nc];
	  windat->waterss[c] += wflux;
	}
	/* Add flow to 2-d array                                     */
	c2 = window->m2d[c];
	windat->waterss2d[c2] += wflux;
      }
    }
  }

  /* Add precipitation and evaporation. This is added to the 2D part */
  /* here so that elevation can be set in the 2D mode. The 3D part   */
  /* is set in vel_w_update() when the sparse coordinates of the     */
  /* surface are known (i.e. after set_dz() is called).              */
  if (wincon->saltflux & ORIGINAL) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      /* Convert from m s-1 to m3 s-1                                */
      windat->waterss2d[c] += (windat->nsfd[c] * window->cellarea[c]);
    }
  }
}

/* END ss_water()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to evaluate additional water fluxes due to sources or     */
/* sinks in the window.                                              */
/*-------------------------------------------------------------------*/
void ss_momentum(geometry_t *window,  /* Window geometry             */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon,  /* Window constants            */
		 int mode             /* 2D or 3D mode               */
		 )
{
  int s, cc, c, cs, cb, c2, zp1, nc;
  double u1, u2;
  double u1av, u2av;
  double *depth_e1 = wincon->d5;           /* Set in precalc_u1()    */
  double *depth_e2 = wincon->d6;           /* Set in precalc_u2()    */

  if (wincon->npss) {
    int s;
    for (s = 0; s < wincon->npss; s++) {
      pss_t *p = &wincon->pss[s];
      double dens0 = 0.0;
      double zhigh_i;
      double zlow_i;
      int zp1;

      if (p->u1tsid < 0 && p->u2tsid < 0) continue;

      /*-------------------------------------------------------------*/
      /* Get the source velocity                                     */
      u1 = u2 = u1av = u2av = 0.0;
      if (p->u1tsid >= 0)
	u1 = ts_multifile_eval(p->ntsfiles, p->ts, p->u1tsid, windat->t);
      if (p->u2tsid >= 0)
	u2 = ts_multifile_eval(p->ntsfiles, p->ts, p->u2tsid, windat->t);

      for (nc = 0; nc < p->vc; nc++) {
	/*-----------------------------------------------------------*/
	/* Get the source location                                   */
	pss_locate(p, windat->t);
	zhigh_i = p->zhigh;
	zlow_i  = p->zlow;
	c = p->e1[nc];
	cs = p->e2[nc];
	cb = p->e3[nc];
	ref_depth(window, windat, wincon, p, &zlow_i, &zhigh_i, cs);
	c2 = window->m2d[c];
	zp1 = window->zp1[c];

	/* Not in the interior of the model, so ignore it! */
	if (c == 0)
	  continue;

	/*-----------------------------------------------------------*/
	/* Momentum line source.                                     */
	/* Force = (source density)(source flow)(source velocity) /  */
	/*         (water density)(cell volume)                      */
	/*       = rate of momentum input per unit mass              */
	if (p->watertsid >= 0) {
	  double dfrac, frac, pdz, d1;
	  double wflux = windat->wflux[s];
	  double a1 = window->h1au1[c2] * window->h2au1[c2];
	  double a2 = window->h1au2[c2] * window->h2au2[c2];

	  /* Source density                                          */
	  if (windat->tflux[windat->tno][s] != NOTVALID) {
	    double temp = windat->tflux[windat->tno][s] / wflux;
	    if (windat->tflux[windat->sno][s] != NOTVALID) {
	      double salt = windat->tflux[windat->sno][s] / wflux;
	      eos2(salt, temp, 0.0, &d1, &dens0);
	    }
	  }

	  /* Add momentum to the 2D mode.                            */
	  if (mode & VEL2D) {
	    double dens = wincon->densavu1[c] /window->h1au1[c2];
	    dfrac = (dens0 > 0.0) ? dens0 / dens : 1.0;
	    d1 = windat->rampval * windat->dt2d * dfrac * wflux;
	    if (u1) windat->nu1av[c2] += (d1 * u1 / (depth_e1[c2] * wincon->mdx[c2] * a1));
	    if (u2) windat->nu2av[c2] += (d1 * u2 / (depth_e2[c2] * wincon->mdy[c2] * a2));
	    return;
	  }

	  zlow_i = max(zlow_i, window->botz[c2]);
	  zhigh_i = min(zhigh_i, windat->topz[c2]);
	  pdz = zhigh_i - zlow_i;

	  /* Add momentum to the 3D mode.                            */
	  /* Add the water, either distributed vertically or into    */
	  /* one layer.                                              */
	  if (pdz > 0.0) {
	    /* Distributed - loop over the affected layers           */
	    c = c2;
	    while (c != window->zm1[c]) {
	      double ctop = (c == cs) ? windat->topz[c2] : window->gridz[zp1];
	      double cbot = (c == cb) ? window->botz[c2] : window->gridz[c];
	      double zlow;
	      double zhigh;
	      double frac;
	      
	      /* Range of source/sink in this cell                     */
	      if (zlow_i <= cbot)
		zlow = cbot;
	      else if (zlow_i >= ctop)
		zlow = ctop;
	      else
		zlow = zlow_i;
	      if (zhigh_i <= cbot)
		zhigh = cbot;
	      else if (zhigh_i >= ctop)
		zhigh = ctop;
	      else
		zhigh = zhigh_i;
	      /* Fraction of flow in this cell                         */
	      frac = (zhigh - zlow) / pdz;
	      
	      if (frac != 0.0) {
		dfrac = (dens0 > 0.0) ? dens0 / windat->dens[c] : 1.0;
		d1 = windat->rampval * windat->dt * dfrac * frac * wflux / pdz;
		/* Update the 3D mode                                  */
		if (u1) windat->nu1[c] += (d1 * u1 / a1);
		if (u2) windat->nu2[c] += (d1 * u2 / a2);
	      }
	      zp1 = c;
	      c = window->zm1[c];
	    }
	  } else {
	    /* Input into a single layer.                              */
	    dfrac = (dens0 > 0.0) ? dens0 / windat->dens[c] : 1.0;
	    d1 = windat->rampval * windat->dt * dfrac * wflux;
	    if (u1) windat->nu1[c] += (d1 * u1 / (windat->dzu1[c] * a1));
	    if (u2) windat->nu2[c] += (d1 * u2 / (windat->dzu2[c] * a2));
	  }
	} else {
	  /*---------------------------------------------------------*/
	  /* Add source velocity directly.                           */
	  zp1 = window->zp1[c];
	  /* Add momentum to the 3D mode.                            */
	  while (c != window->zm1[c]) {
	    if(window->cellz[c] <= zhigh_i && window->cellz[c] >= zlow_i) {
	      if (u1) {
		if (mode & VEL3D) windat->nu1[c] += u1;
		u1av += u1 * windat->dzu1[c] * wincon->mdx[c2];
	      }
	      if (u2) {
		if (mode & VEL3D) windat->nu2[c] += u2;
		u1av += u2 * windat->dzu2[c] * wincon->mdy[c2];
	      }
	    }
	    zp1 = c;
	    c = window->zm1[c];
	  }
	  if (mode & VEL2D) {
	    /* Add momentum to the 2D mode.                            */
	    if (u1) windat->nu1av[c2] += (u1av / (depth_e1[c2] * wincon->mdx[c2]));
	    if (u2) windat->nu2av[c2] += (u2av / (depth_e2[c2] * wincon->mdy[c2]));
	  }
	}
      }
    }
  }
}

/* END ss_momentum()                                                 */
/*-------------------------------------------------------------------*/

#define FLUX_TR   0x01   /* Flux input: tracer has flux specified    */
#define FLUX_NOTR 0x02   /* Flux input: no flux associated with tracer  */
#define CONC_TR   0x04   /* Flow input: tracer has conc. specified   */
#define CONC_NOTR 0x08   /* Flow input: no conc. specified for T/S   */
#define NOVAL     0x10   /* Flow input: no conc. associated with tracer */
#define CONC_UNIT 0x20   /* Flow input: no conc. specified unit tracer  */

/*-------------------------------------------------------------------*/
/* Routine to evaluate any additional sources/sinks of tracer.       */
/* NOTE: the values in dtracer are the fluxes (kg/s) of tracer OUT   */
/* of each cell;                                                     */
/*-------------------------------------------------------------------*/
void ss_tracer(geometry_t *window,  /* Window geometry               */
               window_t *windat,    /* Window data                   */
               win_priv_t *wincon,  /* Window constants              */
               int n,               /* Tracer number                 */
	       double *dtracer,     /* Flux divergence               */
	       double dt            /* Time step                     */
	       )
{
  int s;

  /*-----------------------------------------------------------------*/
  /* Loop over the pointsourcesinks                                  */
  for (s = 0; s < wincon->npss; s++) {
    pss_t *p = &wincon->pss[s];
    int c, c2, cb, cs, nc;
    int zp1;
    double pdz;
    double zhigh_i;
    double zlow_i;
    double sf;
    double val;
    int trf = NOVAL;

    /* Set the type of action to take on this tracer                 */
    if (p->watertsid >= 0) {
      if (windat->tflux[n][s] == NOTVALID) {
	if (n == windat->tno || n == windat->sno)
	  trf = CONC_NOTR;
        else if (strcmp(wincon->trname[n], "unit") == 0)
	  trf = CONC_UNIT;
	else
	  trf = NOVAL;
      } else
	trf = CONC_TR;
    } else {
      if (windat->tflux[n][s] == NOTVALID) 
	trf = FLUX_NOTR;
      else
	trf = FLUX_TR;
    }

    /*---------------------------------------------------------------*/
    /* Loop over horizontal locations                                */
    for (nc = 0; nc < p->vc; nc++) {
      /* Find position of source/sink horizontally and vertically,   */
      /* and truncate to available water column. It would be nice if */
      /* this didn't have to be done repeatedly for each tracer      */
      /* (each time this routine is called).                         */
      pss_locate(p, windat->t);
      zhigh_i = p->zhigh;
      zlow_i  = p->zlow;
      c = p->e1[nc];
      cs = p->e2[nc];
      cb = p->e3[nc];
      ref_depth(window, windat, wincon, p, &zlow_i, &zhigh_i, cs);
      c2 = window->m2d[c];
      zp1 = window->zp1[c];

      /* Not in the  interior of the model, so ignore it             */
      if (c == 0)
	continue;

      zlow_i = max(zlow_i, window->botz[c2]);
      zhigh_i = min(zhigh_i, windat->topz[c2]);
      pdz = zhigh_i - zlow_i;

      /* If water flows out, we need to remove tracer from model     */
      /* regardless of whether or not the source/sink explicitly     */
      /* involves this tracer.                                       */
      if (p->watertsid >= 0 && windat->wflux[s] < 0.0) {
	/* Is sink distributed vertically?                           */
	if (pdz > 0.0) {
	  /* Loop over the affected layers                           */
	  c = c2;
	  while (c != window->zm1[c]) {
	    /* Top and bottom of this cell                           */
	    double ctop = (c == cs) ? windat->topz[c2] : window->gridz[zp1];
	    double cbot = (c == cb) ? window->botz[c2] : window->gridz[c];
	    double zlow;
	    double zhigh;
	    double frac;
	    
	    /* Range of source/sink in this cell                     */
	    if (zlow_i <= cbot)
	      zlow = cbot;
	    else if (zlow_i >= ctop)
	      zlow = ctop;
	    else
	      zlow = zlow_i;
	    if (zhigh_i <= cbot)
	      zhigh = cbot;
	    else if (zhigh_i >= ctop)
	      zhigh = ctop;
	    else
	      zhigh = zhigh_i;

	    /* Fraction of flow in this cell                         */
	    frac = (zhigh - zlow) / pdz;
	    dtracer[c] -= frac * windat->wflux[s] * windat->tr_wc[n][c] * dt ;
	    zp1 = c;
	    c = window->zm1[c];
	  }
	} else {
	  /* Single layer                                            */
	  dtracer[p->e1[nc]] -= windat->wflux[s] * windat->tr_wc[n][c] * dt;
	}
      } else {
	/* Here we know that we don't have a sink for water, so we   */
	/* need to check whether this source/sink involves this      */
	/* tracer.                                                   */
	if (trf & (FLUX_NOTR|NOVAL)) continue;
	/* Is this distributed vertically?                           */
	if (pdz > 0.0) {
	  /* Loop over the affected layers                           */
	  c = c2;
	  while (c != window->zm1[c]) {
	    /* Top and bottom of this cell                           */
	    double ctop = (c == cs) ? windat->topz[c2] : window->gridz[zp1];
	    double cbot = (c == cb) ? window->botz[c2] : window->gridz[c];
	    double zlow;
	    double zhigh;
	    double frac;

	    /* Scaling factor                                        */
	    sf = 1.0;
	    if (p->flag & PSS_AW)
	      sf = window->cellarea[cs];
	    if (p->flag & PSS_VW)
	      sf = (window->cellarea[cs] * wincon->dz[c]);

	    /* Range of source/sink in this cell                     */
	    if (zlow_i <= cbot)
	      zlow = cbot;
	    else if (zlow_i >= ctop) 
	      zlow = ctop;
	    else
	      zlow = zlow_i;
	    if (zhigh_i <= cbot)
	      zhigh = cbot;
	    else if (zhigh_i >= ctop)
	      zhigh = ctop;
	    else
	      zhigh = zhigh_i;

	    /* Fraction of flow in this cell                         */
	    frac = (zhigh - zlow) / pdz;

	    /* If this source/sink is not associated with temp or    */
	    /* sal and a volume influx is specified (trf=CONC_NOTR)  */
	    /* then assume tracer is input with the ambient water    */
	    /* column concentration. If this is not done, then these */
	    /* tracer concentrations will change due to violations   */
	    /* of volume conservation, which subsequently impacts    */
	    /* density and flow.                                     */

	    /* For the FFSL transport mode with volume inputs, use   */
	    /* the volume fluxes in each cell saved to wmean and     */
	    /* copied to waterss. Note that frac uses ctop based on  */
	    /* topz, which is the oldeta value and is out of sync in */
	    /* transport mode. Using the saved waterss directly      */
	    /* ensures exact volume flux input in transport mode.    */
	    if (p->watertsid >= 0) {
	      if (master->trasc == FFSL && master->conserve & CONS_PSS)
		val = windat->waterss[c] * window->cellarea[c2];
	      else
		val = (trf & CONC_UNIT) ? frac * windat->wflux[s] : frac;
	    } else
	      val = frac;

	    /* Update the tracer divergence                          */
	    if (trf & CONC_NOTR)
	      dtracer[c] -= frac * sf * windat->wflux[s] * windat->tr_wc[n][c] * dt;
	    else if (trf & CONC_UNIT)
	      dtracer[c] -= val * sf * dt;
	    else {
	      dtracer[c] -= windat->tflux[n][s] * val * dt;
	    }
	    zp1 = c;
	    c = window->zm1[c];
	  }
	} else {
	  if (p->watertsid >= 0) {
	    if (master->trasc == FFSL && master->conserve & CONS_PSS)
	      val = windat->waterss[c] * window->cellarea[c2];
	    else
	      val = (trf & CONC_UNIT) ? windat->wflux[s] : 1.0;
	  } else
	    val = 1.0;

	  /* Scaling factor                                          */
	  sf = 1.0;
	  if (p->flag & PSS_AW)
	    sf = window->cellarea[p->e2[nc]];
	  if (p->flag & PSS_VW)
	    sf = (window->cellarea[p->e2[nc]] * wincon->dz[p->e1[nc]]);

	  /* Single layer                                            */
	  if (trf & CONC_NOTR)
	    dtracer[p->e1[nc]] -= sf * windat->wflux[s] * windat->tr_wc[n][c] * dt;
	  else if (trf & CONC_UNIT)
	    dtracer[p->e1[nc]] -= sf * val * dt;
	  else
	    dtracer[p->e1[nc]] -= windat->tflux[n][s] * val * dt;
	}
      }
    }
  }
  

  /* Add tracers fluxes due to eta relaxation.                       */
  /* Note : dt is the subbstepping timestep, hence for substepping   */
  /* only add the fraction (dt/windat->dt) of total mass.            */
  /* Note : eta_rlx3d (m3) is the flow out of the cell (negative for */
  /* flow into the cell).                                            */
  if (!(wincon->etarlx & NONE)) {
    int c, cc, cs;
    double sgn = (wincon->sigma) ? -1.0 : 1.0;
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      dtracer[c] += sgn * (wincon->eta_rlx3d[cs] * windat->tr_wc[n][c] * dt /
		     windat->dt);
    }
  }
}

/* END ss_tracer()                                                   */
/*-------------------------------------------------------------------*/




/** Reads a list of point sources and sinks from an ascii file.
  * The sources sinks are specified as shown below.
  * pssnn means pss followed by the integer point source/sink number nn
  * (with the minimum number of digits needed).
  * Here, S is a string (which must not contain whitespace),
  * N is an integer and X, Y, Z and A are floating point
  * numbers.
  *
  * This routine depends on the existence of a global variable:
  *
\begin{verbatim}
# Number of point source/sinks
npointss     N

# Parameters for each point source/sink
# Point source/sink name
pssnn.name    S

# Location ( x y z)
pssnn.location X Y Z

# Data - the next line is a time series definition
# as used by my timeseries routines in sjwlib. The example
# below assumes that the data is in an ascii or netCDF file.
pssnn.data  filename
\end{verbatim}
  *
  * @param name Name of ascii file containing point source/sink list
  * @param t_units Time units to be used for time series
  * @param pss Returned pointer to point source/sink list
  * @param np Returned number of point sources/sinks
  * @param xyzijk Pointer to function which converts (x,y,z) to (i,j,k)
  * @param trI Pointer to function which finds model index of tracer name
  */
void hd_pss_read(char *name,        /* File name                      */
		 char *t_units,     /* Time units for time series     */
		 pss_t **pss,       /* Pointer to source/sink list    */
		 int *np,           /* Number of point source/sinks   */
		 void *model_data,  /* Data for xyzijk and trI        */
		 int (*xyzijk) (void *, double, double, double, int *, int *, int *),
		 int (*trI) (void *, char *), 
		 int (*pss_vec) (void *, pss_t *),
		 int (*pss_region) (void *, char *, pss_t *)
		 )
{
  FILE* fp;
  FILE* fpp;
  int npss = 0;
  pss_t *p = NULL;
  int i, j;
  int n;
  int t, id_counter;
  char buf[MAXLINELEN];
  char key[MAXLINELEN];

  /*-----------------------------------------------------------------*/
  /* Open the file                                                   */
  if ((fp = fopen(name, "r")) == NULL)
    hd_quit("pss_read: Can't open %s\n", name);
  if(!prm_read_char(fp, "npointss", buf)) {
    if(prm_read_char(fp, "pointssfile", buf)) {
      fclose(fp);
      if ((fpp = fopen(buf, "r")) == NULL)
	hd_quit("pss_read: Can't open %s\n", buf);
    } else {
      npss =0;
      *pss = p;
      *np = npss;
      return;
    }
  } else
    fpp = fp;

  /*-----------------------------------------------------------------*/
  /* Get the number of point inputs                                  */
  prm_read_int(fpp, "npointss", &npss);

  /*-----------------------------------------------------------------*/
  /* Allocate memory for list of point inputs                        */
  if (npss > 0 && ((p = (pss_t *)malloc(npss * sizeof(pss_t))) == NULL))
    hd_quit
      ("pss_read: Can't allocate memory for point source/sink list\n");

  /*-----------------------------------------------------------------*/
  /* Read each point input                                           */
  for (i = 0, j = 0; i < npss; i++) {

    p[j].vc = 0;

    /* Store index routine pointer and data needed by it             */
    p[j].xyzijk = xyzijk;
    p[j].model_data = model_data;

    /* Name                                                          */
    sprintf(key, "pss%d.name", i);
    prm_read_char(fpp, key, p[j].name);

    /* Location                                                      */
    sprintf(key, "pss%d.location", i);
    prm_read_char(fpp, key, buf);
    if (!parse_ss_location(&p[j], buf))
      continue;

    /* Regions or blocks                                             */
    if (p[j].x == NOTVALID && p[j].y == NOTVALID) {     
      sprintf(key, "pss%d.region", i);
      if (prm_read_char(fpp, key, buf)) {
        if (!((*pss_region) (model_data, buf, &p[j])))
          continue;
      }
      sprintf(key, "pss%d.ncells", i);
      if (prm_skip_to_end_of_key(fpp, key)) {
	read_blocks(fpp, key, &p[j].vc, &p[j].iloc, &p[j].jloc);
	if (!((*pss_vec) (model_data, &p[j])))
	  continue;
      }
    }

    /* Data                                                          */
    sprintf(key, "pss%d.data", i);
    prm_read_char(fpp, key, buf);
    p[j].ntsfiles = ts_multifile_read(buf, p[j].ts);

    /* Check data time units                                         */
    for (n=0; n<p[j].ntsfiles; n++) {
      if (strcmp(p[j].ts[n].t_units, t_units) != 0) {
	ts_convert_time_units(&p[j].ts[n], t_units);
      }
    }

    /* Flags                                                         */
    p[j].flag = 0;
    sprintf(key, "pss%1d.flag", i);
    if (prm_read_char(fpp, key, buf)) {
      if (contains_token(buf, "NONE") != NULL) {
	p[j].flag = NONE;
      } else {
	if (contains_token(buf, "AREA_WEIGHTED") != NULL)
	  p[j].flag |= PSS_AW;
	if (contains_token(buf, "VOL_WEIGHTED") != NULL)
	  p[j].flag |= PSS_VW;
      }
    }

    /* Include variable vertical references for the time series.     */
    /* 0 = referenced to mean sea level                              */
    /* 1 = referenced to the bottom                                  */
    /* 2 = referenced to the free surface                            */
    p[j].v_offset = 0;
    sprintf(key, "pss%1d.reference", i);
    if (prm_read_char(fpp, key, buf)) {
      if (strcmp(buf, "msl") == 0 || strcmp(buf, "MSL") == 0)
        p[j].v_offset = 0;
      if (strcmp(buf, "bottom") == 0 || strcmp(buf, "BOTTOM") == 0)
        p[j].v_offset = 1;
      if (strcmp(buf, "surface") == 0 || strcmp(buf, "SURFACE") == 0)
        p[j].v_offset = 2;
    }

    /* Find variable indices for each of the time series variables   */
    p[j].nv  = ts_multifile_get_nv(p[j].ntsfiles, p[j].ts);
    p[j].tmap = i_alloc_1d(p[j].nv);
    p[j].watertsid = -1;
    p[j].u1tsid = -1;
    p[j].u2tsid = -1;

    id_counter = 0;
    for (t = 0; t < p[j].ntsfiles; t++) {
      for (n = 0; n < p[j].ts[t].nv; n++) {
	if (strcasecmp(p[j].ts[t].varname[n], "water") == 0 ||
	    strcasecmp(p[j].ts[t].varname[n], "flow") == 0) {
	  p[j].watertsid = id_counter;
	}
	if (strcasecmp(p[j].ts[t].varname[n], "u1") == 0 ||
	    strcasecmp(p[j].ts[t].varname[n], "U1") == 0) {
	  p[j].u1tsid = id_counter;
	}
	if (strcasecmp(p[j].ts[t].varname[n], "u2") == 0 ||
	    strcasecmp(p[j].ts[t].varname[n], "U2") == 0) {
	  p[j].u2tsid = id_counter;
	}

	/* Cache the tracer index to the master                       */
	p[j].tmap[id_counter] = (*trI) (model_data, p[j].ts[t].varname[n]);
	if (
	    p[j].tmap[id_counter] < 0 && /* invalid flux id           */
	    n != p[j].ts[t].ti        && /* and not time variable     */
	    id_counter != p[j].watertsid /* and not the water/flow id */
	    )
	  emstag(LWARN,"pss_read","%s in %s in file %s isn't a model tracer!\n",
		 p[j].ts[t].varname[n], p[j].name, p[j].ts[t].name);
	/*
	  warn("pss_read: %s in %s in file %s isn't a model tracer!\n",
	       p[j].ts[t].varname[n], p[j].name, p[j].ts[t].name);
	*/
	/* Increment the global (across tsfiles) id counter */
	id_counter++;
      }
    }
    j++;
  }

  npss = j;

  /* Close the file                                                  */
  fclose(fpp);

  /* Store pointer to list of point source/sinks                     */
  *pss = p;
  *np = npss;
}

static void pss_locate(pss_t *p, double t)
{
  /* Check if the source/sink has a time varying location            */
  if (p->loc == NULL)
    return;

  /* Get x, y and z values                                           */
  p->x = ts_eval(p->loc, p->x_id, t);
  p->y = ts_eval(p->loc, p->y_id, t);
  p->zlow = ts_eval(p->loc, p->zl_id, t);
  p->zhigh = ts_eval(p->loc, p->zh_id, t);
  p->z = (p->zlow + p->zhigh) / 2;

  /* Convert to model indices                                        */
  if ((*(p->xyzijk))
      (p->model_data, p->x, p->y, p->z, &p->e1[0], &p->e2[0], &p->e3[0]) < 0)
    hd_quit
      ("pss_locate: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",
       p->name, p->x, p->y, p->z);
}

static int parse_ss_location(pss_t *p, char *line)
{
  int e1, e2, e3;
  FILE *fp;

  /* Clear location time series pointer                              */
  p->loc = NULL;
  p->vc = 0;

  /* Time independent, range of z values                             */
  if (sscanf(line, "%lf %lf %lf %lf", &p->x, &p->y, &p->zlow, &p->zhigh) ==
      4) {
    if (p->zlow >= p->zhigh)
      hd_quit("pss_read: %s has bad z range (must be low then high)\n",
              p->name);
    p->z = (p->zlow + p->zhigh) / 2;
    if ((*(p->xyzijk)) (p->model_data, p->x, p->y, p->z, &e1, &e2, &e3) <=
        0) {
      emstag(LWARN,"pss_read","%s location (%.10g,%.10g,%.10g) can't be converted to indices\n",
         p->name, p->x, p->y, p->z);
      /*
      warn
        ("pss_read: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",
         p->name, p->x, p->y, p->z);
      */
      return 0;
    }
    p->vc = 1;
    p->e1 = i_alloc_1d(p->vc);
    p->e2 = i_alloc_1d(p->vc);
    p->e3 = i_alloc_1d(p->vc);
    p->e1[0] = e1;
    p->e2[0] = e2;
    p->e3[0] = e3;
  }
  /* Time independent, single z value                                */
  else if (sscanf(line, "%lf %lf %lf", &p->x, &p->y, &p->z) == 3) {
    if ((*(p->xyzijk)) (p->model_data, p->x, p->y, p->z, &e1, &e2, &e3) <=
        0) {
      emstag(LTRACE,"pss_read","%s location (%.10g,%.10g,%.10g) can't be converted to indices\n",
         p->name, p->x, p->y, p->z);
      /*
      warn
        ("pss_read: %s location (%.10g,%.10g,%.10g) can't be converted to indices\n",
         p->name, p->x, p->y, p->z);
      */
      return 0;
    }
    p->vc = 1;
    p->e1 = i_alloc_1d(p->vc);
    p->e2 = i_alloc_1d(p->vc);
    p->e3 = i_alloc_1d(p->vc);
    p->e1[0] = e1;
    p->e2[0] = e2;
    p->e3[0] = e3;
    p->zlow = p->z;
    p->zhigh = p->z;
  }
  else if (sscanf(line, "%lf %lf", &p->zlow, &p->zhigh) == 2) {
    if (p->zlow >= p->zhigh)
      hd_quit("pss_read: %s has bad z range (must be low then high)\n",
              p->name);
    p->z = (p->zlow + p->zhigh) / 2;
    p->x = p->y = NOTVALID;
  }
  /* Time dependent, all coords from a time series                   */
  else if ((fp = fopen(line, "r+")) != NULL) {
    fclose(fp);

    if ((p->loc = (timeseries_t *)malloc(sizeof(timeseries_t))) == NULL)
      hd_quit("pss_read: No memory for location time series\n");
    ts_read(line, p->loc);
    /* Get coordinate indices */
    if ((p->x_id = ts_get_index(p->loc, "x")) < 0 &&
        (p->x_id = ts_get_index(p->loc, "X")) < 0)
      hd_quit
        ("pss_read: Location timeseries must have a variable named x (or X)\n");
    if ((p->y_id = ts_get_index(p->loc, "y")) < 0 &&
        (p->y_id = ts_get_index(p->loc, "Y")) < 0)
      hd_quit
        ("pss_read: Location timeseries must have a variable named y (or Y)\n");
    if ((p->zl_id = ts_get_index(p->loc, "z_low")) < 0 ||
        (p->zh_id = ts_get_index(p->loc, "z_high")) < 0)
      hd_quit
        ("pss_read: Location timeseries must have variables named z_low and z_high\n");

    p->vc = 1;
    p->e1 = i_alloc_1d(p->vc);
    p->e2 = i_alloc_1d(p->vc);
    p->e3 = i_alloc_1d(p->vc);
    /* Set location at some arbitrary initial time. Probably don't really
       need to do this here, as any routines using the structure should
       call the location routine themselves. */
    pss_locate(p, 0.0);
  } else {
    hd_quit("pss_read: Can't understand location definition for %s\n",
            p->name);
  }
  return 1;
}


/*-------------------------------------------------------------------*/
/* Free memory associted with source/sinks                           */
/*-------------------------------------------------------------------*/
void sourcesink_free(master_t *master, 
		     geometry_t **window, 
		     window_t **windat, 
		     win_priv_t **wincon)
{
  int n, s, t;

  if (master->wflux) d_free_1d(master->wflux);
  if (master->tflux) d_free_2d(master->tflux);

  for (s = 0; s < master->npss; s++) {
    pss_t *pss = &master->pss[s];
    for (t = 0; t < pss->ntsfiles; t++) {
      hd_ts_free(master, &pss->ts[t]);
    }
    if (master->pss[s].loc)
      hd_ts_free(master, pss->loc);
    free((pss_t *)pss);
  }
  /*free((pss_t **)master->pss);*/
  for (n = 1; n <= master->nwindows; n++) {
    if (windat[n]->wflux) d_free_1d(windat[n]->wflux);
    if (windat[n]->tflux) d_free_2d(windat[n]->tflux);
    for (s = 0; s < wincon[n]->npss; s++) {
      pss_t *p = &wincon[n]->pss[s];
      if (p->vc) {
	i_free_1d(p->e1);
	i_free_1d(p->e2);
	i_free_1d(p->e3);
      }
      free((pss_t *)p);
    }
    /*free((pss_t **)wincon[n]->pss);*/
  }
}

/* END sourcesink_free()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-initialise source/sinks                                        */
/*-------------------------------------------------------------------*/
void sourcesink_reinit(master_t *master,     /* Master data          */
		       geometry_t **window,  /* Window geometry      */
		       window_t **windat,    /* Window data          */
		       win_priv_t **wincon   /* Window constants     */
		       )
{
  parameters_t *params = master->params;

  sourcesink_free(master, window, windat, wincon);
  sourcesink_init(params, master, window, windat, wincon);
  master->regf &= ~RS_PSSSET;
}

/* END sourcesink_reinit()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the vertical cell location given a depth.                 */
/*-------------------------------------------------------------------*/
void ref_depth(geometry_t *window,
	       window_t *windat,
	       win_priv_t *wincon,
	       pss_t *pss,         /* Pointsource data structure     */
	       double *zlow,       /* Lower depth                    */
	       double *zhigh,      /* Upper depth                    */
	       int cs              /* Surface coordinate             */
	       )

{
  double d1;

  if (pss->v_offset == 0) {
    *zlow = pss->zlow;
    *zhigh = pss->zhigh;
  }  else if (pss->v_offset == 1) {
    *zlow = window->botz[cs] * wincon->Ds[cs] + fabs(pss->zlow);
    *zhigh = window->botz[cs] * wincon->Ds[cs] + fabs(pss->zhigh);
  } else {
    *zlow += windat->eta[cs];
    *zhigh += windat->eta[cs];
  }
  if (*zlow > *zhigh) {
    d1 = *zlow;
    *zlow = *zhigh;
    *zhigh = d1;
  }
}

void ref_depth_m(master_t *master,
		 pss_t *pss,         /* Pointsource data structure   */
		 double *zlow,       /* Lower depth                  */
		 double *zhigh,      /* Upper depth                  */
		 int cs              /* Surface coordinate           */
		 )

{
  geometry_t *geom = master->geom;
  double d1;

  if (pss->v_offset == 0) {
    *zlow = pss->zlow;
    *zhigh = pss->zhigh;
  }  else if (pss->v_offset == 1) {
    *zlow = geom->botz[cs] * master->Ds[cs] + fabs(pss->zlow);
    *zhigh = geom->botz[cs] * master->Ds[cs] + fabs(pss->zhigh);
  } else {
    *zlow += master->eta[cs];
    *zhigh += master->eta[cs];
  }
  if (*zlow > *zhigh) {
    d1 = *zlow;
    *zlow = *zhigh;
    *zhigh = d1;
  }
}

/* END ref_depth()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets horizontal locations of soursinks on the master              */
/*-------------------------------------------------------------------*/
int hd_pss_vec_m(void *data, pss_t *pss)
{
  return pss_vec_m((master_t *)data, pss);
}

int pss_vec_m(master_t *master, pss_t *pss)
{
  geometry_t *geom = master->geom;
  int c, cs, cc, n, i, ci;

  if (pss->vc) {
    n = pss->vc;
    pss->vc = 0;
    /* Count valid locations                                         */
    for(cc = 1; cc <= n; cc++) {
      ci = geom->cc2s[pss->iloc[cc]];
      c = min(ci, geom->szc - 1);
      if (c > 0 && c < geom->szc) pss->vc++;
    }
    /* Allocate and populate                                         */
    pss->e1 = i_alloc_1d(pss->vc);
    pss->e2 = i_alloc_1d(pss->vc);
    pss->e3 = i_alloc_1d(pss->vc);
    pss->vc = 0;
    for(cc = 1; cc <= n; cc++) {
      ci = geom->cc2s[pss->iloc[cc]];
      c = min(ci, geom->szc - 1);
      if (c > 0 && c < geom->szc){
	cs = geom->m2d[c];
	pss->e1[pss->vc] = c;
	pss->e2[pss->vc] = cs;
	pss->e3[pss->vc] = geom->bot_t[geom->c2cc[cs]];
	pss->vc++;
      }
    }

    return 1;
  }
  return 0;
}

/* END pss_vec_m()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets horizontal locations of soursinks in each window             */
/*-------------------------------------------------------------------*/
int hd_pss_vec_w(void *data, pss_t *pss)
{
  return pss_vec_w((geometry_t *)data, pss);
}

int pss_vec_w(geometry_t *window, pss_t *pss)
{
  int c, cs, cc, i, ci;
  int nn, n = pss->vc;
  int *mask;

  if (n != 0) {
    mask = i_alloc_1d(n);
    pss->vc = 0;
    for(nn = 0; nn < n; nn++) {
      mask[nn] = 0;
      ci = geom->cc2s[pss->iloc[nn+1]];
      for(cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	if (window->wsa[c] == ci) {
	  mask[pss->vc] = c; 
	  pss->vc++;	  
	  break;
	}
      }
    }
    if (pss->vc) {
      pss->e1 = i_alloc_1d(pss->vc);
      pss->e2 = i_alloc_1d(pss->vc);
      pss->e3 = i_alloc_1d(pss->vc);
      for(cc = 0; cc < pss->vc; cc++) {
	c = mask[cc];
	cs = window->m2d[c];
	pss->e1[cc] = c;
	pss->e2[cc] = cs;
	pss->e3[cc] = window->bot_t[window->c2cc[cs]];
      }
      i_free_1d(mask);
      return 1;
    }
    i_free_1d(mask);
    return 0;
  }
  return 0;
}

/* END pss_vec_w()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets sourcesink locations for a region in the master              */
/*-------------------------------------------------------------------*/
int hd_pss_region_m(void *data, char *dname, pss_t *pss)
{
  return pss_init_region_m((master_t *)data, dname, pss);
}


int pss_init_region_m(master_t *master, char *dname, pss_t *pss)
{
  geometry_t *geom = master->geom;
  char *files[MAXSTRLEN * MAXNUMARGS];
  int n, nf, nr, rgn, cc, c;
  double *regionid;

  nf = parseline(dname, files, MAXNUMARGS);
  regionid = d_alloc_1d(geom->sgsiz);
  nr = read_regioni(master, files[0], regionid);
  if (nr) {
    pss->vc = 0;
    for (n = 1; n < nf; n++) {
      sscanf(files[n], "%d", &rgn);
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	if (rgn == (int)regionid[c])
	  pss->vc++;
      }
    }
    if (pss->vc) {
      pss->e1 = i_alloc_1d(pss->vc);
      pss->e2 = i_alloc_1d(pss->vc);
      pss->e3 = i_alloc_1d(pss->vc);
      pss->vc = 0;
      for (n = 1; n < nf; n++) {
	sscanf(files[n], "%d", &rgn);
	for (cc = 1; cc <= geom->b2_t; cc++) {
	  c = geom->w2_t[cc];
	  if (rgn == (int)regionid[c]) {
	    pss->e1[pss->vc] = c;
	    pss->e2[pss->vc] = c;
	    pss->e3[pss->vc] = geom->bot_t[geom->c2cc[c]];
	    pss->vc++;
	  }
	}
      }
    } else {
      hd_warn("pss_init_regions: Can't find valid regions in file %s\n", files[0]);
      d_free_1d(regionid);
      return(0);
    }
    d_free_1d(regionid);
    return(1);
  }
  d_free_1d(regionid);
  return(0);
}

/* END pss_init_regions_m()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets sourcesink locations for a region in the master              */
/*-------------------------------------------------------------------*/
int hd_pss_region_w(void *data, char *dname, pss_t *pss)
{
  return pss_init_region_w((geometry_t *)data, dname, pss);
}


int pss_init_region_w(geometry_t *window, char *dname, pss_t *pss)
{
  char *files[MAXSTRLEN * MAXNUMARGS];
  int n, nf, nr, rgn, cc, c;
  double *regionid;
  nf = parseline(dname, files, MAXNUMARGS);
  regionid = d_alloc_1d(master->geom->sgsiz);
  nr = read_regioni(master, files[0], regionid);
  if (nr) {
    pss->vc = 0;
    for (n = 1; n < nf; n++) {
      sscanf(files[n], "%d", &rgn);
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	if (rgn == (int)regionid[window->wsa[c]])
	  pss->vc++;
      }
    }
    if (pss->vc) {
      pss->e1 = i_alloc_1d(pss->vc);
      pss->e2 = i_alloc_1d(pss->vc);
      pss->e3 = i_alloc_1d(pss->vc);
      pss->vc = 0;
      for (n = 1; n < nf; n++) {
	sscanf(files[n], "%d", &rgn);
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  if (rgn == (int)regionid[window->wsa[c]]) {
	    pss->e1[pss->vc] = c;
	    pss->e2[pss->vc] = c;
	    pss->e3[pss->vc] = window->bot_t[window->c2cc[c]];
	    pss->vc++;
	  }
	}
      }
    } else {
      d_free_1d(regionid);
      return(0);
    }
    d_free_1d(regionid);
    return(1);
  }
  d_free_1d(regionid);
  return(0);
}

/* END pss_init_regions_w()                                          */
/*-------------------------------------------------------------------*/
