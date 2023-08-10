/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/momentum/mommisc.c
 *  
 *  Description:
 *  Calculate boundary u2 values
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: mommisc.c 6803 2021-06-23 01:36:08Z her127 $
 *
 */

#include "hd.h"

/*-------------------------------------------------------------------*/
/* Get the cells and weights within the sponge zones                 */
/*-------------------------------------------------------------------*/
void set_sponge_cells(geometry_t *window)
{
  int nn, n, cc, c;
  int ee, e, eb, ep;
  int i1, i2, cn, cb, j;
  double d1, d2, dist, dmin;
  int nscm, *scm;
  int *mask, *eask;
  double deg2m = 60.0 * 1852.0;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int verbose = 0;

  mask = i_alloc_1d(window->szm);
  eask = i_alloc_1d(window->szm);
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->sponge_zone_h) {
      d1 = 0.0;
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	d1 += sqrt(window->cellarea[c]);
      }
      d1 /= (double)open->no2_t;
      if (d1 > (double)open->sponge_zone_h) 
	hd_quit("Sponge zone for %s (%5.2f m) is less than mean grid size (%5.2f m). Increase sponge zone.\n",
		open->name, open->sponge_zone_h, d1);

      /*-------------------------------------------------------------*/
      /* Count the cell centres in the zone                          */
      open->nspc = 0;
      memset(mask, 0, window->szm * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	for (i1 = 1; i1 <= open->no2_t; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (!mask[c] && dist <= (double)open->sponge_zone_h) {
	    open->nspc++;
	    mask[c] = 1;
	  }
	}
      }
      
      /*-------------------------------------------------------------*/
      /* Get the cells in the zone                                   */
      open->spc = i_alloc_1d(open->nspc + 1);
      open->swc = d_alloc_1d(open->nspc + 1);
      open->snc = i_alloc_1d(open->nspc + 1);
      open->smc = i_alloc_1d(open->nspc + 1);
      open->nspc = 1;
      memset(mask, 0, window->szm * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	for (i1 = 1; i1 <= open->no2_t; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (!mask[c] && dist <= (double)open->sponge_zone_h) {
	    open->spc[open->nspc++] = c;
	    mask[c] = 1;
	  }
	}
      }
      open->nspc--;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the mask                                                 */
  memset(mask, 0, window->szm * sizeof(int));
  for (cc = window->v2_t+1; cc <= window->n2_t; cc++) {
    c = window->w2_t[cc];
    mask[c] = 1;
  }
  for (cc = window->b2_t+1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    mask[c] = 1;
  }
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      c = open->obc_e2[ee];
      cn = open->ogc_t[ee];
      mask[c] = mask[cn] = 1;
    }
    if (open->sponge_zone_h) {
      /*-------------------------------------------------------------*/
      /* Get the closest boundary cell to sponge cells               */
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	if (window->cask) window->cask[c] |= W_SOBC;
	dmin = HUGE;
	for (i1 = 1; i1 <= open->no2_t; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (dist <= dmin) {
	    dmin = dist;
	    open->swc[cc] = dist;
	    open->snc[cc] = i2;
	  }
	}
	mask[c] = 1;
      }
    }
  }

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->sponge_zone_h) {
      /*-------------------------------------------------------------*/
      /* Find the cells on the outer limit of the zone               */
      nscm = 0;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  cn = window->c2c[j][c];
	  if (!mask[cn]) {
	    nscm++;
	    break;
	  }
	}
      }
      scm = i_alloc_1d(nscm + 1);
      memset(scm, 0, (nscm + 1) * sizeof(int));
      nscm = 1;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  cn = window->c2c[j][c];
	  if (!mask[cn]) {
	    scm[nscm++] = c;
	    break;
	  }
	}
      }
      nscm--;

      /*-------------------------------------------------------------*/
      /* Get the closest cell on the outer perimeter, its distance   */
      /* to the closest cell on the boundary, and the ratio of       */
      /* distances (weight). The value in the sponge zone is then:   */
      /* v = swc * (v(perimeter) - v(boundary))  + v(boundary)       */
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	dmin = HUGE;
	for (i1 = 1; i1 <= nscm; i1++) {
	  i2 = scm[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (dist <= dmin) {
	    dmin = dist;
	    open->smc[cc] = i2;

	    cn = open->snc[cc];
	    d1 = window->cellx[cn] - window->cellx[c];
	    d2 = window->celly[cn] - window->celly[c];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swc[cc] = dist;

	    d1 = window->cellx[cn] - window->cellx[i2];
	    d2 = window->celly[cn] - window->celly[i2];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swc[cc] = (dist) ? open->swc[cc] / dist : 0.0;
	  }
	}
      }
      if (verbose) {
	for (cc = 1; cc <= open->nspc; cc++) {
	  c = open->spc[cc];
	  printf("wn=%d c=(%d %d) mn=(%d %d) mx=(%d %d) %f\n",window->wn,window->s2i[c],window->s2j[c],
		 window->s2i[open->snc[cc]],window->s2j[open->snc[cc]],
		 window->s2i[open->smc[cc]],window->s2j[open->smc[cc]],open->swc[cc]);
	}
      }

      /*-------------------------------------------------------------*/
      /* Count the edges in the zone                                 */
      memset(eask, 0, window->szm * sizeof(int));
      open->nspe1 = 0;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  e = window->c2e[j][c];
	  if (!eask[e]) {
	    open->nspe1++;
	    eask[e] = 1;
	  }
	}
      }

      open->spe1 = i_alloc_1d(open->nspe1 + 1);
      open->swe1 = d_alloc_1d(open->nspe1 + 1);
      open->sne1 = i_alloc_1d(open->nspe1 + 1);
      open->sme1 = i_alloc_1d(open->nspe1 + 1);
      memset(eask, 0, window->szm * sizeof(int));
      open->nspe1 = 1;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  e = window->c2e[j][c];
	  if (!eask[e]) {
	    cn = open->snc[cc];

	    /* Get the edge corresponding to the closest cell on the */
	    /* boundary (i.e. the closest edge on the boundary).     */
	    for (ee = 1; ee <= open->no2_e1; ee++) {
	      eb = open->obc_e1[ee];
	      cb = open->obc_e2[ee];
	      if (cb == cn) {
		open->sne1[open->nspe1] = eb;
		break;
	      }
	    }

	    /* Distance to the closest boundary edge                 */
	    d1 = window->u1x[e] - window->u1x[eb];
	    d2 = window->u1y[e] - window->u1y[eb];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swe1[open->nspe1] = dist;
	    /* Get the edge corresponding to the closest cell on the */
	    /* outer perimeter (i.e. the closest edge on the outer   */
	    /* perimeter).                                           */
	    dmin = HUGE;
	    for (i1 = 1; i1 <= nscm; i1++) {
	      i2 = scm[i1];
	      for (ee = 1; ee <= window->npe[i2]; ee++) {
		ep = window->c2e[ee][i2];
		d1 = window->u1x[e] - window->u1x[ep];
		d2 = window->u1y[e] - window->u1y[ep];
		dist = sqrt(d1 * d1 + d2 * d2);
		if (is_geog) dist *= deg2m;
		if (dist < dmin) {
		  dmin = dist;
		  open->sme1[open->nspe1] = ep;
		}
	      }
	    }
	    if(n==1&&e==2696) {
	    }
	    ep = open->sme1[open->nspe1];
	    d1 = window->u1x[ep] - window->u1x[eb];
	    d2 = window->u1y[ep] - window->u1y[eb];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swe1[open->nspe1] = (dist) ? open->swe1[open->nspe1] / dist : 0.0;
	    if (window->eask) window->eask[e] |= W_SOBC;
	    open->spe1[open->nspe1++] = e;
	    eask[e] = 1;
	  }
	}
      }
      open->nspe1--;
      i_free_1d(scm);
    }
  }
  i_free_1d(mask);
  i_free_1d(eask);
}

/* END set_sponge_cells()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone at cell centres where friction increases       */
/* linearly towards the boundary.                                    */
/*-------------------------------------------------------------------*/
void set_sponge_c(geometry_t *window, /* Window geometry             */
		  double *AH,         /* Viscosity                   */
		  double dt           /* 3D time step                */
		  )
{
  open_bdrys_t *open;           /* Open boudary data structure       */
  int c, cc, cs, cb, cp;        /* Centre counters                   */
  double vb, vp, vm;            /* Diffusivities                     */
  double sfact = 0.9;           /* Safety factor                     */
  int  n;                       /* Boundary counter                  */

  for (n = 0; n < window->nobc; ++n) {
    int cb, cp;
    open = window->open[n];
    if (open->sponge_zone_h) {
      for (cc = 1; cc <= open->nspc; cc++) {
	c = cs = open->spc[cc];
	cb = open->snc[cc];
	cp = open->smc[cc];
	while (c != window->zm1[c]) {
	  vp = AH[cp];
	  vm = sfact * window->cellarea[cs] / (4.0 * dt);
	  vm = 1.0 / ((2.0 / window->cellarea[cs]) * 4.0 * dt);
	  if (open->sponge_f)
	    vb = open->sponge_f * vp;
	  else {
	    vb = vm;
	  }
	  AH[c] = min(vm, open->swc[cc] * (vp - vb) + vb);
	  c = window->zm1[c];
	  cb = window->zm1[cb];
	  if (cb == window->zm1[cb]) cb = window->zp1[cb];
	  cp = window->zm1[cp];
	  if (cp == window->zm1[cp]) cp = window->zp1[cp];
	}
      }
    }
  }
}

/* END set_sponge_c()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone at cell edges where friction increases         */
/* linearly towards the boundary.                                    */
/*-------------------------------------------------------------------*/
void set_sponge_e(geometry_t *window, /* Window geometry             */
		  double *AH,         /* Viscosity                   */
		  double dt           /* 3D time step                */
 	  )
{
  open_bdrys_t *open;           /* Open boudary data structure       */
  int e, ee, es, eb, ep;        /* Edge counters                     */
  double vb, vp, vm;            /* Diffusivities                     */
  double sfact = 0.9;           /* Safety factor                     */
  int  n;                       /* Boundary counter                  */
  double bif, a;

  for (n = 0; n < window->nobc; ++n) {
    int cb, cp;
    open = window->open[n];
    if (open->sponge_zone_h) {
      for (ee = 1; ee <= open->nspe1; ee++) {
	e = es = open->spe1[ee];
	a = 1.0 / (0.125 * window->edgearea[es]);
	eb = open->sne1[ee];
	ep = open->sme1[ee];
	bif = (window->wincon->diff_scale & SCALEBI) ? 0.125 * window->edgearea[es] : 1.0;
	while (e != window->zm1e[e]) {
	  vp = AH[ep];

	  vm = 1.0 / (window->h1au1[es] * window->h1au1[es]) +
	    1.0 / (window->h2au1[es] * window->h2au1[es]);
	  vm = sfact * bif / (4.0 * vm * dt);

	  if (open->sponge_f)
	    vb = open->sponge_f * vp;
	  else
	    vb = vm;

	  /* v = swc * (v(perimeter) - v(boundary))  + v(boundary)   */
	  AH[e] = min(vm, open->swe1[ee] * (vp - vb) + vb);
	  e = window->zm1e[e];
	  eb = window->zm1e[eb];
	  if (eb == window->zm1e[eb]) eb = window->zp1e[eb];
	  ep = window->zm1e[ep];
	  if (ep == window->zm1e[ep]) ep = window->zp1e[ep];
	}
      }
    }
  }
}

/* END set_sponge_e()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets the sponge zone to original mixing                         */
/*-------------------------------------------------------------------*/
void reset_sponge_zone(geometry_t *window)   /* Window geometry      */
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  open_bdrys_t *open;           /* Open boudary data structure       */ 
  int n, bn, spgn;
  int ee, e, c, es;                     /* Counters, cell locations  */
  int *map;

  if (windat->u1vhin == NULL) return;

  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];
    if (open->sponge_zone_h) {
      for (ee = 1; ee <= open->nspe1; ee++) {
	e = es = open->spe1[ee];
	while (e != window->zm1e[e]) {
	  if (wincon->sue1 > 0.0 && wincon->u1vh0 < 0.0) 
	    wincon->u1vh[e] = wincon->sue1 * windat->sdc[c];
	  else
	    wincon->u1vh[e] = windat->u1vhin[es];
	  e = window->zm1e[e];
	}
      }
    }
  }
}

/* END reset_sponge_zone()                                           */
/*-------------------------------------------------------------------*/
