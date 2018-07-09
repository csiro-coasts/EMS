/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/closure/csanady.c
 *  
 *  Description:
 *  Implements Csanady style mixing scheme
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: csanady.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hd.h"

#define RI_EPS (1e-6)

#if defined(OS_IRIX) || defined(OS_SOLARIS)
#define HAS_CUBE_ROOT
extern double cbrt(double);     /* Not always proto-typed */
#endif


/*-------------------------------------------------------------------*/
/* Routine to initialize the Csanady mixing scheme                   */
/*-------------------------------------------------------------------*/
void closure_csanady_init(geometry_t *geom, /* Sparse global geometry
                                               structure */
                          master_t *master  /* Master data structure */
  )
{
  int c, cc;                    /* Sparse coordinates */

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    master->Vz[c] = master->vz0;
    master->Kz[c] = master->kz0;
  }
}

/* END closure_csanady_init()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the vertical viscosity and diffusivity using */
/* the Csanady closure scheme.                                       */
/*-------------------------------------------------------------------*/
void closure_csanady(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int cc, zm1;                  /* Sparse counter */
  double Ri;                    /* Richardson number */
  double *vel;
  double tmp, zu, zl, dz;

  /* Get the cell centered current speed */
  vel = wincon->w1;
  for (cc = 1; cc <= window->b3_t; cc++) {
    int c = window->w3_t[cc];
    double u1val = windat->u[c];
    double u2val = windat->v[c];
    vel[c] = sqrt(u1val * u1val + u2val * u2val);
  }

  /* Calculate the mixing coefficients */
  for (cc = 1; cc <= window->b2_t; cc++) {
    int c = window->w2_t[cc];
    int cs = c;
    int cb = window->bot_t[cc];
    double top = windat->topz[cs];
    double bot = window->botz[cs];
    double H = top - bot;
    double ustar;
    double wind_ustar = sqrt(windat->wind1[cs] * windat->wind1[cs] +
                             windat->wind2[cs] * windat->wind2[cs]);
    wind_ustar = sqrt(wind_ustar / 1025.0);

    /* Limit depth - FIX - THIS IS AN AWFUL FUDGE */
    if (H < 10)
      H = 10.0;
    if (H > 100)
      H = 100.0;

    while (c != cb) {

      zm1 = window->zm1[c];

      /* Use max of wind ustar, bottom ustar */
      ustar = max(sqrt(wincon->Cd[cs]) * vel[c], wind_ustar);

      /* Z-coord of cell centre, this layer */
      zu = (c == cs) ? (top + window->gridz[cs]) / 2 : window->cellz[c];
      /* Z-coord of cell centre, next layer down */
      zl = (c == window->zp1[cb]) ? (bot + window->gridz[c]) / 2 :
        window->cellz[zm1];
      dz = zu - zl;

      /* Calculate Vz/Kz at cell interfaces */
      tmp = vel[c] - vel[zm1];
      tmp = tmp * tmp * (windat->dens[c] + windat->dens[zm1]) / 2 + RI_EPS;
      Ri = g * (windat->dens_0[zm1] - windat->dens_0[c]) * dz / tmp;
      if (Ri < -0.0001) {
        /* We have a density inversion */
        windat->Kz[c] = 0.1;
        windat->Vz[c] = 0.1;
      } else {
        windat->Kz[c] = wincon->kz0 + wincon->kz_alpha * ustar
          * H / ((1 + 3.33 * Ri) * sqrt(1 + 3.33 * Ri));
        windat->Vz[c] = wincon->vz0 + wincon->vz_alpha * ustar
          * H / sqrt(1.0 + 10.0 * Ri);
      }
      c = zm1;
    }

    /* Store Vz/Kz at bottom (may be used in particle tracking) */
    windat->Kz[cb] = windat->Kz[window->zp1[cb]];
    windat->Vz[cb] = windat->Vz[window->zp1[cb]];
  }
}

/* END closure_csanady()                                             */
/*-------------------------------------------------------------------*/
