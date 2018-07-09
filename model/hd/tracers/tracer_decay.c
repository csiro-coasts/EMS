/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/tracers/tracer_decay.c
 *  
 *  Description:
 *  Calculate tracer decay
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracer_decay.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include "hd.h"
#include "tracer.h"

/*------------------------------------------------------------------*/
/* Decays tracer values                                             */
/*------------------------------------------------------------------*/
void tracer_decay(geometry_t *window,  /* Window geometry           */
                  window_t *windat,    /* Window data               */
                  win_priv_t *wincon,  /* Window constants          */
                  double decay,        /* Decay constant            */
                  double *tracer,      /* Tracer to decay           */
		  tracer_info_t *tr    /* Tracer info               */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  double d;
  double sc;
  int tn = -1;
  char units[MAXSTRLEN];

  if (decay == 0.0)
    return;

  sc = 1.0;
  tn = (int)decay;
  if (tr->flag & DE_TR2)
    strcpy(units, wincon->trinfo_2d[tn].units);
  if (tr->flag & DE_TR3)
    strcpy(units, wincon->trinfo_3d[tn].units);
  if (strcmp(units, "days") == 0 || strcmp(units, "day") == 0)       
    sc = 86400.0;
  if (strcmp(units, "hours") == 0 || strcmp(units, "hour") == 0) 
    sc = 3600.0;
  if (strcmp(units, "minutes") == 0 || strcmp(units, "minute") == 0) 
    sc = 60.0;

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    d = decay;
    if (tr->flag & DE_TR2) d = windat->tr_wcS[tn][window->m2d[c]];
    if (tr->flag & DE_TR3) d = windat->tr_wc[tn][c];
    tracer[c] -= windat->dt * tracer[c] / (d * sc);
  }
}

/* END tracer_decay()                                               */
/*------------------------------------------------------------------*/
