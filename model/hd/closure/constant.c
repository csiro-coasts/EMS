/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/closure/constant.c
 *  
 *  Description:
 *  Implements very simple (constant) mixing scheme
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: constant.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

/*-------------------------------------------------------------------*/
/* Routine to initialize the constant mixing scheme                  */
/*-------------------------------------------------------------------*/
void closure_constant_init(geometry_t *geom,  /* Sparse global geometry
                                                 structure */
                           master_t *master /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coordinates */

  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    master->Vz[c] = master->vz0;
    master->Kz[c] = master->kz0;
  }
}

/* END closure_constant_init()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* constant mixing coefficients.                                     */
/*-------------------------------------------------------------------*/
void closure_constant(geometry_t *window, /* Processing window */
                      window_t *windat, /* Window data structure */
                      win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc, c3;                /* Sparse coodinate / counter */
  int cb;                       /* Bottom sparse coordinate */
  int zm1;                      /* Cell below cell c */

  for (cc = 1; cc <= window->b2_t; cc++) {
    c3 = c = window->nsur_t[cc];
    cb = window->bot_t[cc];

    /*---------------------------------------------------------------*/
    /* Compute cell centred viscosity, Vz, and diffusion, Kz. Use */
    /* the minimum of surface,bottom and q-moment length scales.  */
    c = c3;
    zm1 = window->zm1[c];

    while (c < cb) {
      windat->Vz[c] = wincon->vz0;
      windat->Kz[c] = wincon->kz0;
      c = zm1;
      zm1 = window->zm1[c];
    }
  }
}

/* END closure_constant()                                            */
/*-------------------------------------------------------------------*/
