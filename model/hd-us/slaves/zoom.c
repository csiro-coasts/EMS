/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/slaves/zoom.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: zoom.c 6746 2021-03-30 00:44:44Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void smooth(master_t *master, double *A, int *ctp, int nctp);
double cvolt(geometry_t *window, double *a, int c, int edge);

/*-------------------------------------------------------------------*/
/* Routine to sum e1 fluxes over a grid so as to be applicable to a  */
/* zoomed cell.                                                      */
/*-------------------------------------------------------------------*/
void zflux_e1(geometry_t *geom, /* Sparse global geometery */
              geometry_t *window, /* Processing window */
              double *Ag,       /* Global flux array */
              double *Al,       /* Local flux array */
              int *vec,         /* Cells on which Al is defined */
              int *evec,        /* Cells on which Ag is defined */
              int nvec          /* Size of vec */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */
  int c1, c2;

  /* Transfer the flux in the master to the slave.                  */
  /* zoomed window then sum the fluxes.                             */
  for (cc = 1; cc <= nvec; cc++) {
    lc = vec[cc];
    c = evec[lc];
    Al[lc] = Ag[c];
  }
}

/* END zflux_e1()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths an array using a convolution kernal                       */
/*-------------------------------------------------------------------*/
void smooth(master_t *master, /* Model data structure                */
	    double *A,        /* Array to smooth                     */
	    int *ctp,         /* Cells to process array              */
	    int nctp          /* Size of ctp                         */
  )
{
  geometry_t *geom = master->geom;
  int c, cc;
  double* aa;

  aa = d_alloc_1d(geom->sgsizS);
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    aa[c] = cvol1(master, A, c, 0);
  }
  memcpy(A, aa, geom->sgsizS * sizeof(double));

  d_free_1d(aa);
}

/* END smooth()                                                      */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Smooths an array using a convolution kernal                       */
/*-------------------------------------------------------------------*/
void smooth3(master_t *master, /* Model data structure               */
	     double *A,        /* Array to smooth                    */
	     int *ctp,         /* Cells to process array             */
	     int nctp,         /* Size of ctp                        */
	     int sz,           /* Size of A                          */
	     int edge          /* Edge to smooth                     */
  )
{
  geometry_t *geom = master->geom;
  int c, ci, cc, j;
  double* aa;

  aa = d_alloc_1d(sz);
  memcpy(aa, A, sz * sizeof(double));
  for (cc = 1; cc <= nctp; cc++) {
    ci = c = ctp[cc];
    if (edge) {
      for (j = 1; j <= geom->npe[geom->m2d[ci]]; j++) {
	c = geom->c2e[j][ci];
	aa[c] = cvol1(master, A, ci, j);
      }
    } else
      aa[c] = cvol1(master, A, ci, 0);
  }
  memcpy(A, aa, sz * sizeof(double));

  d_free_1d(aa);
}

/* END smooth3()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths an edge array using area weighted eSe[]                   */
/*-------------------------------------------------------------------*/
void smooth3e(master_t *master, /* Model data structure              */
	      double *A,        /* Array to smooth                   */
	      int *ctp,         /* Cells to process array            */
	      int nctp,         /* Size of ctp                       */
	      int sz            /* Size of A                         */
  )
{
  geometry_t *geom = master->geom;
  int e, es, ee, j, eoe;
  double *aa, d1;

  aa = d_alloc_1d(sz);
  memcpy(aa, A, sz * sizeof(double));
  for (ee = 1; ee <= nctp; ee++) {
    e = ctp[ee];
    es = geom->m2de[e];
    d1 = 0.0;
    for (j = 1; j <= geom->nee[e]; j++) {
      eoe = geom->eSe[j][e];
      es = geom->m2de[eoe];
      if (eoe) {
	aa[e] += A[eoe] * geom->edgearea[es];
	d1 += geom->edgearea[es];
      }
    }
    aa[e] = (d1) ? aa[e]/ d1 : A[e];
  }
  memcpy(A, aa, sz * sizeof(double));
  d_free_1d(aa);
}

/* END smooth3e()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths an array using a convolution kernal                       */
/*-------------------------------------------------------------------*/
void smooth_w(geometry_t *window,  /* Window data structure          */
	      double *A,           /* Array to smooth                */
	      double *B,           /* Dummy array                    */
	      int *ctp,            /* Cells to process array         */
	      int nctp,            /* Size of ctp                    */
	      int n,               /* Smoothing passes               */
	      double lim,          /* Upper limit                    */
	      int edgef            /* 0=centres, 1=edges             */
	      )
{
  int c, ci, cc, nn;
  int en = 0;

  for (nn = 0; nn < n; nn++) {

    memcpy(B, A, window->szm * sizeof(double));
    for (cc = 1; cc <= nctp; cc++) {
      c = ci = ctp[cc];
      if (edgef) {
	en = window->e2e[c][0];
	ci = window->e2c[c][0];
      }
      B[c] = min(lim, cvolt(window, A, ci, en));
    }
    memcpy(A, B, window->szm * sizeof(double));
  }
}

/* END smooth_w()                                                    */
/*-------------------------------------------------------------------*/

double cvolt(geometry_t *window,  /* Model data structure */
             double *a,         /* Array to smooth */
             int c,             /* Cell centre location */
	     int edge           /* Edge to smooth */
  )
{
  double *kk;
  double fi, nf;
  int cc;
  int *st, sz = 3;
  int type = ST_SQ3;

  if( edge > 0) type |= ST_EDGE;

  st = stencil(window, c, &sz, type, edge);

  fi = nf = 0.0;
  for (cc = 0; cc < sz; cc++) {
    if (a[st[cc]] != 0.0) {
      fi += a[st[cc]];
      nf += 1.0;
    }
  }
  if (nf > 0.0) fi /= nf;
  i_free_1d(st);
  return (fi);
}

/* END cvol1()                                                       */
/*-------------------------------------------------------------------*/
