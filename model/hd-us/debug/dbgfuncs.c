/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/debug/dbgfuncs.c
 *  
 *  Description:
 *  Convienience routines for outputing debug
 *  information.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: dbgfuncs.c 6484 2020-03-26 01:14:36Z her127 $
 *
 */

#include <math.h>
#include <stdio.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Returns 1 if the sparse coordinate corresponds to the mesh index  */
/*-------------------------------------------------------------------*/
int is_index(geometry_t *geom,     /* Window or global geometry      */
	     int cc,               /* Mesh index                     */
	     int c                 /* Sparse cell centre coordinate  */
	     )
{
  geometry_t *ggeom = master->geom;
  int ret = 0;

  if (geom == ggeom) {
    if (c == ggeom->cc2s[cc]) ret = 1;
  } else {
    if (geom->wsa[c] == ggeom->cc2s[cc]) ret = 1;
  }
  return(ret);
}

/* END is_index()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to print the i,j,k location given a global sparse         */
/* coordinate.                                                       */
/*-------------------------------------------------------------------*/
void psg(int c)
{
  int i = geom->s2i[c];
  int j = geom->s2j[c];
  int k = geom->s2k[c];
  printf("%d (%d %d %d)\n", c, i, j, k);
}

/* END psg()                                                         */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to print the i,j,k location given a global sparse         */
/* coordinate.                                                       */
/*-------------------------------------------------------------------*/
void psl(geometry_t *window, int cl)
{
  int c = window->wsa[cl];
  int i = geom->s2i[c];
  int j = geom->s2j[c];
  int k = geom->s2k[c];
  hd_warn("%d (%d %d %d)\n", cl, i, j, k);
}

/* END psl()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks to see if a NaN exists in an array                         */
/*-------------------------------------------------------------------*/
void check_nan(geometry_t *window, double *A, int as, int ae,
               int *vec, char *tag)
{
  int c, cc, cg;                /* Local sparse coordinate / counter */
  int i, j, k;                  /* Cartesian coordinates */

  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (isnan(A[c])) {
      cg = window->wsa[c];
      i = geom->s2i[cg];
      j = geom->s2j[cg];
      k = geom->s2k[cg];
      hd_warn("%s NaN at %d (%d %d %d) : %d\n", tag, c, i, j, k, cg);
    }
  }
}

/* END check_nan()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Prints the maximum value in an array                              */
/*-------------------------------------------------------------------*/
void check_max(geometry_t *window, double *A, int as, int ae,
               int *vec, char *tag)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  int i, j, k;                  /* Cartesian coordinates */
  double amax = -1e10;

  cm = cg = i = j = k = 0;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (A[c] > amax) {
      amax = A[c];
      cm = c;
      cg = window->wsa[c];
      i = geom->s2i[cg];
      j = geom->s2j[cg];
      k = geom->s2k[cg];
    }
  }
  hd_warn("%s max = %f at %d (%d %d %d) : %d\n", tag, amax, cm, i, j, k,
          cg);
}


/*-------------------------------------------------------------------*/
/* Prints the absolute maximum value in an array to file fp          */
/*-------------------------------------------------------------------*/
double print_max(geometry_t *window, double *A, int as, int ae,
		 int *vec, char *tag, FILE *fp, int *cl)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  int i, j, k;                  /* Cartesian coordinates */
  double amax = -1e10;

  cm = cg = i = j = k = 0;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (fabs(A[c]) > amax) {
      amax = fabs(A[c]);
      cm = c;
      cg = window->wsa[c];
      i = geom->s2i[cg];
      j = geom->s2j[cg];
      k = geom->s2k[cg];
    }
  }
  if (A[cm] >= 0.0)
    fprintf(fp,"%s %f at %d (%d %d %d) : %d\n", tag, amax, cm, i, j, k,
	    cg);
  else
    fprintf(fp,"%s -%f at %d (%d %d %d) : %d\n", tag, amax, cm, i, j, k,
	    cg);
  *cl = cm;
  return(amax);
}


/*-------------------------------------------------------------------*/
/* Returns the absolute maximum value in an array                     */
/*-------------------------------------------------------------------*/
double get_max(geometry_t *window, double *A, int as, int ae,
	       int *vec, int *cl)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  double amax = -1e10;

  cm = cg = 0;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (fabs(A[c]) > amax) {
      amax = fabs(A[c]);
      *cl = c;
    }
  }
  return(amax);
}


/*-------------------------------------------------------------------*/
/* Returns the minimum value in an array                             */
/*-------------------------------------------------------------------*/
double get_min(geometry_t *window, double *A, int as, int ae,
	       int *vec, int *cl)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  double amin = 1e10;

  cm = cg = 0;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (A[c] < amin) {
      amin = A[c];
      *cl = c;
    }
  }
  return(amin);
}


/*-------------------------------------------------------------------*/
/* Returns the absolute maximum value in an array                     */
/*-------------------------------------------------------------------*/
double get_maxs(geometry_t *window, double *A, int as, int ae,
	       int *vec, int *cl, double vmax)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  double amax = -1e10;

  cm = cg = 1;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (fabs(A[c]) > vmax) {
      cl[cg++] = c;
    }
    if (fabs(A[c]) > amax) {
      amax = fabs(A[c]);
      cl[0] = c;
    }
  }
  return(amax);
}


/*-------------------------------------------------------------------*/
/* Prints the minimum value in an array                              */
/*-------------------------------------------------------------------*/
void check_min(geometry_t *window, double *A, int as, int ae,
               int *vec, char *tag)
{
  int c, cc, cg, cm;            /* Local sparse coordinate / counter */
  int i, j, k;                  /* Cartesian coordinates */
  double amax = 1e10;

  cm = cg = i = j = k = 0;
  for (cc = as; cc <= ae; cc++) {
    c = vec[cc];
    if (A[c] < amax) {
      amax = A[c];
      cm = c;
      cg = window->wsa[c];
      i = geom->s2i[cg];
      j = geom->s2j[cg];
      k = geom->s2k[cg];
    }
  }
  hd_warn("%s min = %f at %d (%d %d %d) : %d\n", tag, amax, cm, i, j, k,
          cg);
}

/* END check_min()                                                   */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to test the continuity equation. The output value should  */
/* be around 1e-8.                                                   */
/*-------------------------------------------------------------------*/
double check_continuity(geometry_t *window, int c)
{
  window_t *windat = window->windat;
  double hf, vf;
  int j, e, cs = window->m2d[c];
  int zp1 = window->zp1[c];

  /* Calculate and print the continuity balance if required. */
  double wtop = (c == cs) ? windat->wtop[cs] : windat->w[zp1];
  hf = 0.0;
  for (j = 1; j <= window->npe[cs]; j++) {
    hf -= window->eSc[j][cs] * windat->u1flux3d[e];
  }
  vf = ((wtop - windat->w[c])) * window->cellarea[cs];
  fprintf(stderr, "continuity %f %d : %e\n", windat->t / 86400,
          geom->s2k[c], hf + vf);
  return (hf + vf);
}

/* END check_cont()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks the continuity (volume) balance for correct volume fluxes  */
/* and sea levels. Should be ~1e-10.                                 */
/*-------------------------------------------------------------------*/
void check_volcons(geometry_t *window,       /* Window geometry      */
		   window_t *windat,         /* Window data          */
		   win_priv_t *wincon,       /* Window constants     */
		   double *oeta,             /* eta at t-1           */
		   double *eta,              /* eta at t             */
		   double dt                 /* Timestep             */
		)
{
  double colflux;               /* Flux divergence                   */
  int c, ci, cc, c2, ee, e, es; /* Sparse coodinate / counter        */
  int co, cs, cb, zp1;
  double *u1, *w, *u1flux;
  double volo, vol, d1, flux, cin, cout;

  if (wincon->means & TRANSPORT) {
    u1 = windat->ume;
    w = windat->wm;
    u1flux = windat->u1vm;
  } else {
    u1 = windat->u1;
    w = windat->w;
    u1flux = windat->u1flux3d;
  }

  /* Loop over surface layers                                        */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    co = wincon->i3[cc];
    cs = wincon->i2[cc];
    cb = wincon->i1[cc];
    c2 = window->m2d[c];
    while (c != window->zm1[c]) {
      cin = cout = volo = 0.0;
      if (c == cs) {
	double watertop;
	double waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
	double dz = eta[c2] - waterbot;
	ci = c;
	vol = dz * window->cellarea[c2];
	/* First increment the old volume from the current layer to  */
	/* the layer below that containing the old surface. If       */
	/* elevation has risen this loop is skipped.                 */
	while (ci != co) {
	  zp1 = window->zp1[ci];

	  /* Horizontal fluxes                                       */
	  for (ee = 1; ee <= window->npe[c2]; ee++) {
	    e = window->c2e[ee][ci];
	    es = window->m2de[e];
	    flux = window->eSc[ee][c2] * u1flux[e] * dt;
	    if (flux > 0.0) {
	      cout += flux;
	    } else if (flux < 0.0) {
	      cin -= flux;
	    }
	  }

	  /* Old volume                                              */
	  watertop = window->gridz[zp1];
	  d1 = window->cellarea[c2] * (watertop - waterbot);
	  volo += d1;

	  ci = zp1;
	  waterbot = window->gridz[ci];
	}
	/* Now increment the old volume for the layer containing the */
	/* old surface.                                              */
	waterbot = (ci == cb) ? window->botz[c2] : window->gridz[ci];
	d1 = (oeta[c2] - waterbot) * window->cellarea[c2];
	volo += d1;

	/* Horizontal fluxes                                         */
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][ci];
	  es = window->m2de[e];
	  flux = window->eSc[ee][c2] * u1flux[e] * dt;
	  if (flux > 0.0)
	    cout += flux;
	  else if (flux < 0.0)
	    cin -= flux;
	}

	/* Increment the fluxes from above the old surface elevation */
	/* to the top of the vertical grid. If the surface has risen */
	/* then these cells now contain water and horizontal         */
	/* divergence may alter the total amount of tracer. If       */
	/* elevation has dropped then there is no water (and hence   */
	/* no divergence) in these cells.                            */
	if (ci != c2) {
	  do {
	    ci = window->zp1[ci];
	    for (ee = 1; ee <= window->npe[c2]; ee++) {
	      e = window->c2e[ee][ci];
	      es = window->m2de[e];
	      flux = window->eSc[ee][c2] * u1flux[e] * dt;
	      if (flux > 0.0)
		cout += flux;
	      else if (flux < 0.0)
		cin -= flux;
	    }
	  } while (ci != c2);
	}
      } else {
	vol = volo = window->cellarea[c2] * wincon->dz[c];
	/* Horizontal fluxes                                         */
	for (ee = 1; ee <= window->npe[c2]; ee++) {
	  e = window->c2e[ee][c];
	  es = window->m2de[e];
	  flux = window->eSc[ee][c2] * u1flux[e] * dt;
	  if (flux > 0.0)
	    cout += flux;
	  else if (flux < 0.0)
	    cin -= flux;
	}
	/* Vertical fluxes in the upper layer                        */
	zp1 = window->zp1[c];
	flux = w[zp1] * window->cellarea[c2] * dt;
	if (flux < 0.0)
	  cin -= flux;
	else if (flux > 0.0)
	  cout += flux;
      }
      /* Vertical fluxes in the lower layer                          */
      flux = w[c] * window->cellarea[c2] * dt;
      if (flux > 0.0)
	cin += flux;
      else if (flux < 0.0)
	cout -= flux;
      if (windat->volcont) windat->volcont[c] = (volo - vol) - ((cout-cin));
      c = window->zm1[c];
    }
  }
}

/* END check_volcons()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void check_bathy(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int cc, c, zp1, c2;
  double d1, d2;
  double minval = 1e-10;        /* Minimum value for velocity */

  memset(windat->u1, 0, window->sgsiz * sizeof(double));
  memset(windat->u2, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->v3_e1; cc++) {
    c = window->w3_e1[cc];
    windat->u1[c] = 0.5;
  }
  set_dz_at_u1(window, windat, wincon);
  set_flux_3d(window, windat, wincon, VEL3D);
  vel_w_update(window, windat, wincon);

  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    zp1 = window->zp1[c];
    c2 = window->m2d[c];
    d1 =
      (cc <= window->v2_t) ? 0.0 : 0.5 * (windat->w[c] + windat->w[zp1]);
    if (fabs(d1) < minval)
      d1 = minval;
    d2 = fabs(d1 * windat->dt / (wincon->dz[c] * wincon->Ds[c2]));
    if (d2 > 2) {
      hd_warn("cfl at (%d %d %d) %f %f\n", geom->s2i[c], geom->s2j[c],
              geom->s2k[c], d2, d1);
    }
  }
}

/*-------------------------------------------------------------------*/



void printflags(FILE * out, unsigned long f)
{
  if (f & U1SOLID)
    fprintf(out, "U1SOLID ");
  if (f & U2SOLID)
    fprintf(out, "U2SOLID ");
  if (f & U1OUTSIDE)
    fprintf(out, "U1OUTSIDE ");
  if (f & U2OUTSIDE)
    fprintf(out, "U2OUTSIDE ");
  if (f & U1BDRY)
    fprintf(out, "U1BDRY ");
  if (f & U2BDRY)
    fprintf(out, "U2BDRY ");
  if (f & L_EDGE)
    fprintf(out, "L_EDGE ");
  if (f & R_EDGE)
    fprintf(out, "R_EDGE ");
  if (f & B_EDGE)
    fprintf(out, "B_EDGE ");
  if (f & F_EDGE)
    fprintf(out, "F_EDGE ");
  if (f & ETASPEC)
    fprintf(out, "ETASPEC ");
  if (f & DRY)
    fprintf(out, "DRY (No longer used) ");
  if (f & SOLID)
    fprintf(out, "SOLID ");
  if (f & OUTSIDE)
    fprintf(out, "OUTSIDE ");
  if (f & ALLWATER)
    fprintf(out, "ALLWATER ");
  if (f & TRACERSPEC)
    fprintf(out, "TRACERSPEC ");
  if (f & U1TFLUX)
    fprintf(out, "U1TFLUX ");
  if (f & U2TFLUX)
    fprintf(out, "U2TFLUX ");
  if (f & U1AVZERO)
    fprintf(out, "U1AVZERO ");
  if (f & U2AVZERO)
    fprintf(out, "U2AVZERO ");
  fprintf(out, "\n");
}
