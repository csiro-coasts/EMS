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
 *  $Id: zoom.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void smooth(master_t *master, double *A, int *ctp, int nctp);

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
  int c, ci, cc;
  double* aa;

  aa = d_alloc_1d(sz);
  memcpy(aa, A, sz * sizeof(double));
  for (cc = 1; cc <= nctp; cc++) {
    ci = ctp[cc];
    c = (!edge) ? ci : geom->c2e[edge][ci];
    aa[c] = cvol1(master, A, ci, edge);
  }
  memcpy(A, aa, sz * sizeof(double));

  d_free_1d(aa);
}

/* END smooth3()                                                     */
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
	      double lim           /* Upper limit                    */
	      )
{
  int c, cc, nn;

  for (nn = 0; nn < n; nn++) {

    memcpy(B, A, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= nctp; cc++) {
      c = ctp[cc];
      B[c] = min(lim, cvol1(master, A, c, 0));
    }
    memcpy(A, B, window->sgsiz * sizeof(double));
  }
}

/* END smooth_w()                                                    */
/*-------------------------------------------------------------------*/
