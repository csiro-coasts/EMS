/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/debug/dbgfuncs.c
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
 *  $Id: dbgfuncs.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <math.h>
#include <stdio.h>
#include "hd.h"

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
/* Routine to get the 2D streamfunction                              */
/*-------------------------------------------------------------------*/
void streamfunction(geometry_t *window, /* Processing window */
                    window_t *windat, /* Window data structure */
                    win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc;

  if (wincon->vorticity & RELATIVE) {
    memset(windat->rv, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->nbe1S; cc++) {
      c = window->bpte1S[cc];
      if (c == window->xp1[c])
        while (c != window->xm1[c]) {
          c = window->xm1[c];
          windat->rv[c] += windat->u2av[c];
        }
    }
  }
  if (wincon->vorticity & POTENTIAL) {
    memset(windat->pv, 0, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->nbe2S; cc++) {
      c = window->bpte2S[cc];
      if (c == window->yp1[c])
        while (c != window->ym1[c]) {
          c = window->ym1[c];
          windat->pv[c] -= windat->u1av[c];
        }
    }
  }
}

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
  int cs = window->m2d[c];
  int zp1 = window->zp1[c];

  /* Calculate and print the continuity balance if required. */
  double wtop = (c == cs) ? windat->wtop[cs] : windat->w[zp1];
  hf = windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
    windat->u2flux3d[window->yp1[c]] - windat->u2flux3d[c];
  vf = ((wtop - windat->w[c])) * window->cellarea[cs];
  fprintf(stderr, "continuity %f %d : %e\n", windat->t / 86400,
          geom->s2k[c], hf + vf);
  return (hf + vf);
}

/* END check_cont()                                                  */
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
  set_dz_at_u2(window, windat, wincon);
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
