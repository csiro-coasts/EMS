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
 *  $Id: mommisc.c 6738 2021-03-30 00:40:45Z her127 $
 *
 */

#include "hd.h"

/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the boundary  */
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


/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the boundary  */
/*-------------------------------------------------------------------*/
void set_sponge(geometry_t *window,   /* Window geometry             */
		double *u1vh,         /* e1 viscosity                */
		double dt,            /* 3D time step                */
		double *mask          /* 3D dummy array              */
		)
{
  open_bdrys_t *open;           /* Open boudary data structure       */
  int c, cc, cs;                /* Centre counters                   */
  int e, es, eb, ee;            /* Edge counters                     */
  double sfact = 0.9;           /* Safety factor                     */
  double tfact = 0.3;           /* tanh() multiplier : small = long  */
  int bn, sc, n;                /* Boundary counter                  */
  int j, jo, ji;                /* Edge counter                      */
  double rm, vh;                /* Dummies                           */
  int spgn;                     /* Sponge zone                       */
  int *map;                     /* Interior map                      */
  int mode = 4;                 /* 1: normal to boundary only        */
                                /* 2: all directions from boundary   */
                                /* 4: as for 2 with lateral fill     */

  memset(mask, 0, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Increase in horizontal friction in the sponge: u1vh             */
  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];
    if (open->sponge_zone_h) {
      if (mode == 1) {
	/*-----------------------------------------------------------*/
	/* Blend in a direction normal to the boundary edge only     */
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = eb = open->obc_e1[ee];
	  es = window->m2de[e];
	  if (mask[e]) continue;
	  spgn = open->sponge_zone_h;
	  map = open->ceni[ee] ? window->em : window->ep;
	  if (open->sponge_f)
	    rm = open->sponge_f * u1vh[e];
	  else {
	    rm = 1.0 / (window->h1au1[es] * window->h1au1[es]) +
	      1.0 / (window->h2au1[es] * window->h2au1[es]);
	    rm = sfact / (4.0 * rm * dt);
	  }
	  for (bn = 1; bn <= spgn; bn++) {
	    vh = u1vh[e];
	    /* Linear ramp to rm on the boundary                     */
	    u1vh[e] = (fabs(u1vh[e]) - rm) * 
	      (double)(bn - 1) / (double)spgn + rm;
	    /* Hyperbolic tangent ramp to rm on the boundary         */
	    /*
	    u1vh[e] =  vh + (rm - vh) * 
	      (1 - tanh(tfact * (double)(bn - 1)));
	    */
	    mask[e] = 1.0;
	    e = map[e];
	  }
	}
      } else {
	/*-----------------------------------------------------------*/
	/* Blend in all directions away from the boundary            */
	int cb;     
	spgn = open->sponge_zone_h;
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = eb = open->obc_e1[ee];
	  mask[e] = 1.0;
	}
	/* First pass; all directions from the boundary cell         */
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = eb = open->obc_e1[ee];
	  es = window->m2de[e];
	  jo = open->outi[ee];
	  if (open->sponge_f)
	    rm = open->sponge_f * u1vh[e];
	  else {
	    rm = 1.0 / (window->h1au1[es] * window->h1au1[es]) +
	      1.0 / (window->h2au1[es] * window->h2au1[es]);
	    rm = sfact / (4.0 * rm * dt);
	  }
	  cb = open->obc_e2[ee];
	  cs = window->m2d[cb];
	  for (j = 1; j <= window->npe[cs]; j++) {
	    c = cb;
	    for (bn = 1; bn <= spgn; bn++) {
	      e = window->c2e[jo][c];
	      if (!e) continue;
	      if (e == eb || (!mask[e] && !window->wgst[c])) {
		/* Linear ramp to rm on the boundary                 */
		u1vh[e] = (fabs(u1vh[e]) - rm) * 
		  (double)(bn - 1) / (double)spgn + rm;
		/* Hyperbolic tangent ramp to rm on the boundary     */
		/*	
		  vh = u1vh[e];
		  u1vh[e] =  vh + (rm - vh) * 
		  (1 - tanh(tfact * (double)(bn - 1)));
		*/
		mask[e] = (double)jo;
	      }
	      c = window->c2c[j][c];
	    }
	  }
	}

	/* Second pass; all directions except j from current cell    */
	if (mode == 4) {
	  int cs1, bn1, j1;
	  for (ee = 1; ee <= open->no3_e1; ee++) {
	    e = eb = open->obc_e1[ee];
	    es = window->m2de[e];
	    jo = open->outi[ee];
	    ji = open->ini[ee];
	    cb = open->obc_e2[ee];
	    for (bn = 1; bn < spgn; bn++) {
	      cb = window->c2c[ji][cb];
	      cs = window->m2d[cb];
	      for (j = 1; j <= window->npe[cs]; j++) {
		if (j == ji || j == jo) continue;
		c = cb;
		rm = u1vh[window->c2e[jo][c]];
		for (bn1 = 1; bn1 <= spgn; bn1++) {
		  e = window->c2e[jo][c];
		  if (!e) continue;
		  if (!mask[e] && !window->wgst[c]) {
		    u1vh[e] = (fabs(u1vh[e]) - rm) * 
		      (double)(bn1 - 1) / (double)spgn + rm;
		  }
		  c = window->c2c[j][c];
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/* END set_sponge()                                                  */
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
      spgn = open->sponge_zone_h;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
	es = window->m2de[e];
	c = open->obc_e2[ee];
	map = open->ceni[ee] ? window->em : window->ep;
	for (bn = 1; bn <= spgn; bn++) {
	  if (wincon->sue1 > 0.0 && wincon->u1vh0 < 0.0) 
	    wincon->u1vh[e] = wincon->sue1 * windat->sdc[c];
	  else
	    wincon->u1vh[e] = windat->u1vhin[es];
	  e = map[e];
	  es = window->m2de[e];
	}
      }
    }
  }
}

/* END reset_sponge_zone()                                           */
/*-------------------------------------------------------------------*/
