/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/debug/totals.c
 *  
 *  Description:
 *  Calculates total water, salt and heat in a grid.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: totals.c 5845 2018-06-29 04:08:51Z riz008 $
 *
 */

#include "hd.h"

void totals(geometry_t *window, /* Processing window */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon, /* Window geometry / constants */
            char *s)
{
  int c, cc, c2, cb;            /* Sparse coodinate / counter */
  double vol;
  double totsalt = 0.0;         /* total salt */
  double tottemp = 0.0;         /* total temp */
  double totwater = 0.0;        /* total water */


  totsalt = 0.0;
  tottemp = 0.0;
  totwater = 0.0;
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cb = wincon->i1[cc];
    vol = (windat->topz[c2] - window->botz[c2]) * window->cellarea[c2];
    totwater += vol;
    while (c <= cb) {
      vol = window->cellarea[c2] * wincon->dz[c];
      if (windat->sal)
        totsalt += vol * windat->sal[c];
      if (windat->temp)
        tottemp += vol * windat->temp[c];
      c = window->zm1[c];
    }
  }

  if (DEBUG("totals")) {
    if (s)
      dlog("totals", "Totals at %s\n", s);
    dlog("totals", "Total water %f\n", totwater);
    dlog("totals", "Total salt %f\n", totsalt);
    dlog("totals", "Total heat %f\n", tottemp);
  }
}

void print_tracer_mass(geometry_t *window,  /* Processing window */
                       window_t *windat,  /* Window data structure */
                       win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, cc, c2, cb;            /* Sparse coodinate / counter */
  int nn;
  double area;

  dlog("totals", "Window %d: tracer mass = ", window->wn);
  for (nn = 0; nn < wincon->ntbdy; ++nn) {
    int n = wincon->tbdy[nn];
    double mwc = 0.0;
    /* double msed = 0.0; */
    int first_NaN = 1;
    dlog("totals", "\n  %s:", wincon->trname[n]);
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      cb = wincon->i1[cc];
      area = window->cellarea[c2];
      while (c <= cb) {
        double m, vol;
        vol = wincon->dz[c] * area;
        m = windat->tr_wc[n][c];
        if (isnan(m) || m < 0.0) {
          if (first_NaN) {
            first_NaN = 0;
            dlog("totals", " NaN or < 0 in (%d,%d,%d)", geom->s2i[c],
                 geom->s2j[c], geom->s2k[c]);
          } else
            dlog("totals", ",(%d,%d,%d)", geom->s2i[c], geom->s2j[c],
                 geom->s2k[c]);
        } else
          mwc += m * vol;
      }
      c = window->zm1[c];
    }
    /* 
       if( cgrid->sednz &&
       !(cgrid->flag[cgrid->nz-1][j][i]&(SOLID|OUTSIDE|JOINT)) ) {
       sediment_t* sm = cgrid->sm[j][i]; for (k = 0; k < cgrid->sednz;
       ++k) { double m = cgrid->tr_sed[n][k][j][i]; if( isnan(m) || m <
       0.0 ) { if ( first_NaN ) { first_NaN = 0; dlog("totals", " NaN or < 
       0 in (%d,%d,%d)", i, j, k); } else dlog("totals", ",(%d,%d,%d)", i, 
       j, k); } else { if (cgrid->tracers[n].dissol) msed += m * area *
       sm->dz[k] * sm->porosity[k]; else if (cgrid->tracers[n].partic)
       msed += m * area * sm->dz[k]; } } } */

    if (!first_NaN)
      dlog("totals", ", value =");
    /* dlog("totals", " (%g, %g, %g)", mwc, msed, mwc+msed); */
    dlog("totals", " (%g)", mwc);
  }
  dlog("totals", "\n");
}

void check_tracers_positiveness(geometry_t *window, /* Processing window */
                                window_t *windat, /* Window data structure 
                                                   */
                                win_priv_t *wincon  /* Window geometry /
                                                       constants */
  )
{
  int c, cc, c2;                /* Sparse coodinate / counter */
  int n, nn;
  for (nn = 0; nn < wincon->ntbdy; ++nn) {
    n = wincon->tbdy[nn];
    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      if (windat->tr_wc[n][c] < 0.0)
        hd_warn
          ("Negative tracer value for %s: tr_wc[%d][%d][%d] = %f, time = %d.\n",
           wincon->trname[n], geom->s2k[c], geom->s2j[c], geom->s2i[c],
           windat->tr_wc[n][c], (int)windat->t);
    }
  }
}
