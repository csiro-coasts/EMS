/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/mommisc.c
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
 *  $Id: mommisc.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include "hd.h"


void blend_var(open_bdrys_t *open, double *var, double *var2, double *mask, int spgn, int c);
void get_blend_val(double *val, int c, int spgn, int *map, double *mask);

/*-------------------------------------------------------------------*/
/* Routine to adjust the 3D velocity so that the vertical integral   */
/* equals the 2D velocity. This routine assumes water depth and      */
/* grid spacings correspond to the start of the time-step, t.        */
/* The adjustment is weighted so that larger adjustments are made to */
/* those velocities having larger absolute values. This is so that   */
/* the adjustment does not overwhelm the solutions to the momentum   */
/* equations if adjustments are un-weighted. Note that if the weight */
/* w[c] = 1 is equivalent to the un-weighted case.                   */
/*-------------------------------------------------------------------*/
void velocity_adjust_weighted(geometry_t *window, /* Window geometry */
			      window_t *windat,   /* Window data     */
			      win_priv_t *wincon /* Window constants */
  )
{
  int c, cc, cs, cb;            /* Sparse coodinate / counter        */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *midx;                 /* SIGMA : depth at the e1 face      */
  double *midy;                 /* SIGMA : depth at the e1 face      */
  double *depth;                /* Water depth at the cell face      */
  double *adjust;               /* Velocity adjustment               */
  int bn;                       /* OBC counter                       */
  open_bdrys_t **open = window->open; /* Window OBC structure        */
  double *vmax, *vmin;
  double *w, *wtot;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e1 faces.                                 */
  /* Set pointers and initialise                                     */
  sum = wincon->d1;
  midx = wincon->d2;
  depth = wincon->d5;           /* Set in precalc_u1()               */
  adjust = wincon->d3;
  vmax = wincon->d4;
  vmin = wincon->d7;
  w = wincon->w4;
  wtot = wincon->w5;

  memset(sum, 0, window->sgsizS * sizeof(double));
  memset(wtot, 0, window->sgsizS * sizeof(double));
  /* The  number of surface  u1 cells to process (vcs1) may decrease */
  /* if a cell dries  during the 2D mode - hence adjust[] may not be */
  /* set  but u1 is still updated (with an old value residing in d3. */
  /* Therefore, zero adjust[] before use here.                       */
  memset(adjust, 0, window->sgsizS * sizeof(double));

  /* Get the maximum and minimum absolute velocities                 */
  memset(vmax, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    vmin[c] = wincon->velmax;
  }
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    vmax[cs] = (fabs(windat->u1[c]) > vmax[cs]) ? 
		fabs(windat->u1[c]) : vmax[cs];
    vmin[cs] = (fabs(windat->u1[c]) < vmin[cs]) ? 
		fabs(windat->u1[c]) : vmin[cs];
  }

  /* Integrate the velocities through the water column. Note: normal */
  /* velocity open boundary sparse coordinates are not included in   */
  /* the cells to process vectors (whereas tangential boundary       */
  /* velocities are) and must be processed separately.               */
  /* Wet cells to process.                                           */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    midx[cs] = wincon->mdx[cs];
    sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
    /* Get the weights                                               */
    /* Weighted to maximum velocities                                */
    w[c] = (vmax[cs]) ? ((fabs(windat->u1[c]) - vmin[cs])) / vmax[cs] : 1.0;
    /* Weighted to minimum velocities                                */
    /*w[c] = (vmin[cs] - vmax[cs]) ? (windat->u1[c] - vmin[cs]) / (vmin[cs] - vmax[cs]) + 1.0 : 1.0;*/
    /* No weighting                                                  */
    /*w[c] = 1.0;*/
    wtot[cs] += w[c] * windat->dzu1[c] * midx[cs];
  }

  /* Open boundary normal velocities. Velocities above the surface   */
  /* are set to zero, so it is possible to integrate over the whole  */
  /* water column rather than just to the free surface.              */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	midx[cs] = wincon->mdx[cs];
	sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
	/* Get the weights                                           */
	w[c] = (vmax[cs]) ? 
	  ((fabs(windat->u1[c]) - vmin[cs])) / vmax[cs] : 1.0;
	/*w[c] = (vmin[cs] - vmax[cs]) ? (windat->u1[c] - vmin[cs]) / (vmin[cs] - vmax[cs]) + 1.0 : 1.0;*/
	wtot[cs] += w[c] * windat->dzu1[c] * midx[cs];
      }
    }
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u1flux[cs] - sum[cs] * windat->dt * window->h2au1[cs]) /
      (wtot[cs] * windat->dt * window->h2au1[cs]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e1; cc++) {
	cs = open[bn]->obc_e1[cc];
	depth[cs] -= window->botzu1[cs];
	adjust[cs] = 
	  (windat->u1flux[cs] - sum[cs] * windat->dt * window->h2au1[cs]) /
	  (wtot[cs] * windat->dt * window->h2au1[cs]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    windat->u1[c] += (adjust[cs] * w[c]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	windat->u1[c] += (adjust[cs] * w[c]);
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e1; cc++) {
    c = cs = wincon->cdry_e1[cc];
    windat->u1av[cs] = 0.0;
    windat->u1flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u1[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u1bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i5[cc];        /* 3D bottom coordinate              */
    windat->u1bot[cs] = windat->u1b[c] - windat->u1avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e1; cc++) {
    c = wincon->cbot_e1[cc];
    cs = window->m2d[c];
    windat->u1[c] = windat->u1av[cs];
    windat->u1bot[cs] = 0.0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e2 faces.                                 */
  /* Set pointers and initialise                                     */
  midy = wincon->d2;
  depth = wincon->d6;           /* Set in precalc_u2()               */
  memset(sum, 0, window->sgsizS * sizeof(double));
  memset(adjust, 0, window->sgsizS * sizeof(double));
  memset(wtot, 0, window->sgsizS * sizeof(double));

  /* Get the maximum and minimum absolute velocities                 */
  memset(vmax, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    vmin[c] = wincon->velmax;
  }
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    vmax[cs] = (fabs(windat->u2[c]) > vmax[cs]) ? 
		fabs(windat->u2[c]) : vmax[cs];
    vmin[cs] = (fabs(windat->u2[c]) < vmin[cs]) ? 
		fabs(windat->u2[c]) : vmin[cs];
  }

  /* Integrate the velocities through the water column.              */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    midy[cs] = wincon->mdy[cs];
    sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
    /* Get the weights                                               */
    w[c] = (vmax[cs]) ? ((fabs(windat->u2[c]) - vmin[cs])) / vmax[cs] : 1.0;
    /*w[c] = (vmax[cs]) ? fabs(windat->u2[c]) / vmax[cs] : 1.0;*/
    /*w[c] = (vmin[cs] - vmax[cs]) ? (windat->u2[c] - vmin[cs]) / (vmin[cs] - vmax[cs]) + 1.0 : 1.0;*/
    /*w[c] = 1.0;*/
    wtot[cs] += w[c] * windat->dzu2[c] * midy[cs];
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	midy[cs] = wincon->mdy[cs];
	sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
	/* Get the weights                                               */
	w[c] = (vmax[cs]) ? 
	  ((fabs(windat->u2[c]) - vmin[cs])) / vmax[cs] : 1.0;
	/*w[c] = (vmax[cs]) ? fabs(windat->u2[c]) / vmax[cs] : 1.0;*/
	/*w[c] = (vmin[cs] - vmax[cs]) ? (windat->u2[c] - vmin[cs]) / (vmin[cs] - vmax[cs]) + 1.0 : 1.0;*/
	wtot[cs] += w[c] * windat->dzu2[c] * midy[cs];
      }
    }
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u2flux[cs] - sum[cs] * windat->dt * window->h1au2[cs]) / 
      (wtot[cs] * windat->dt * window->h1au2[cs]);
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e2; cc++) {
	cs = open[bn]->obc_e2[cc];
	depth[cs] -= window->botzu2[cs];
	adjust[cs] =
	  (windat->u2flux[cs] - sum[cs] * windat->dt * window->h1au2[cs]) /
	  (wtot[cs] * windat->dt * window->h1au2[cs]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->u2[c] += (adjust[cs] * w[c]);
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	windat->u2[c] += (adjust[cs] * w[c]);
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e2; cc++) {
    c = cs = wincon->cdry_e2[cc];
    windat->u2av[cs] = 0.0;
    windat->u2flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u2[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u2bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i6[cc];        /* 3D bottom coordinate              */
    windat->u2bot[cs] = windat->u2b[cb] - windat->u2avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e2; cc++) {
    c = wincon->cbot_e2[cc];
    cs = window->m2d[c];
    windat->u2[c] = windat->u2av[cs];
    windat->u2bot[cs] = 0.0;
  }
}

/* END velocity_adjust_weighted()                                    */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to adjust the 3D velocity so that the vertical integral   */
/* equals the 2D velocity. This routine assumes water depth and      */
/* grid spacings correspond to the start of the time-step, t.        */
/* The adjustment is weighted so that larger adjustments are made to */
/* those velocities having larger differences between the unweighted */
/* adjusted velocity and the un-adjusted velocity. This prevents     */
/* large un-adjusted velocities from becomming larger.               */       
/*-------------------------------------------------------------------*/
void velocity_adjust_weighted_r(geometry_t *window, /* Window geometry */
			      window_t *windat,   /* Window data     */
			      win_priv_t *wincon /* Window constants */
  )
{
  int c, cc, cs, cb;            /* Sparse coodinate / counter        */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *midx;                 /* SIGMA : depth at the e1 face      */
  double *midy;                 /* SIGMA : depth at the e1 face      */
  double *depth;                /* Water depth at the cell face      */
  double *adjust;               /* Velocity adjustment               */
  int bn;                       /* OBC counter                       */
  open_bdrys_t **open = window->open; /* Window OBC structure        */
  double *vmax, *vmin;
  double *w, *wtot, d1, d2;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e1 faces.                                 */
  /* Set pointers and initialise                                     */
  sum = wincon->d1;
  midx = wincon->d2;
  depth = wincon->d5;           /* Set in precalc_u1()               */
  adjust = wincon->d3;
  vmax = wincon->d4;
  vmin = wincon->d7;
  w = wincon->w4;
  wtot = wincon->w5;

  memset(sum, 0, window->sgsizS * sizeof(double));
  memset(wtot, 0, window->sgsizS * sizeof(double));
  /* The  number of surface  u1 cells to process (vcs1) may decrease */
  /* if a cell dries  during the 2D mode - hence adjust[] may not be */
  /* set  but u1 is still updated (with an old value residing in d3. */
  /* Therefore, zero adjust[] before use here.                       */
  memset(adjust, 0, window->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column. Note: normal */
  /* velocity open boundary sparse coordinates are not included in   */
  /* the cells to process vectors (whereas tangential boundary       */
  /* velocities are) and must be processed separately.               */
  /* Wet cells to process.                                           */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    midx[cs] = wincon->mdx[cs];
    sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
  }

  /* Open boundary normal velocities. Velocities above the surface   */
  /* are set to zero, so it is possible to integrate over the whole  */
  /* water column rather than just to the free surface.              */
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	midx[cs] = wincon->mdx[cs];
	sum[cs] += windat->u1[c] * windat->dzu1[c] * midx[cs];
      }
    }
  }

  /* Estimate the velocity adjustment based on no weighting          */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    cs = window->w3_e1[cc];
    adjust[cs] =
      (windat->u1flux[cs] - sum[cs] * windat->dt * window->h2au1[cs]) /
      (depth[cs] * midx[cs] * windat->dt * window->h2au1[cs]);
  }
  /* Get the maximum and minimum absolute velocities differences      */
  memset(vmax, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    vmin[c] = wincon->velmax;
  }
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    d1 = fabs(windat->u1[c] + adjust[cs]);
    vmax[cs] = (d1 > vmax[cs]) ? d1 : vmax[cs];
    vmin[cs] = (d1 < vmin[cs]) ? d1 : vmin[cs];
  }
  /* Get the weights                                                 */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    d1 = fabs(windat->u1[c] + adjust[cs]);
    d2 = vmin[cs] - vmax[cs];
    w[c] = (d2) ? (d1 - vmax[cs]) / d2 : 1.0;
    wtot[cs] += w[c] * windat->dzu1[c] * midx[cs];
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u1flux[cs] - sum[cs] * windat->dt * window->h2au1[cs]) /
      (wtot[cs] * windat->dt * window->h2au1[cs]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e1; cc++) {
	cs = open[bn]->obc_e1[cc];
	depth[cs] -= window->botzu1[cs];
	adjust[cs] = 
	  (windat->u1flux[cs] - sum[cs] * windat->dt * window->h2au1[cs]) /
	  (wtot[cs] * windat->dt * window->h2au1[cs]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    windat->u1[c] += (adjust[cs] * w[c]);
  }
  if(!wincon->dobdry_u1) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e1; cc++) {
	c = open[bn]->obc_e1[cc];
	cs = window->m2d[c];
	windat->u1[c] += (adjust[cs] * w[c]);
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e1; cc++) {
    c = cs = wincon->cdry_e1[cc];
    windat->u1av[cs] = 0.0;
    windat->u1flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u1[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u1bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = wincon->s1[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i5[cc];        /* 3D bottom coordinate              */
    windat->u1bot[cs] = windat->u1b[c] - windat->u1avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e1; cc++) {
    c = wincon->cbot_e1[cc];
    cs = window->m2d[c];
    windat->u1[c] = windat->u1av[cs];
    windat->u1bot[cs] = 0.0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e2 faces.                                 */
  /* Set pointers and initialise                                     */
  midy = wincon->d2;
  depth = wincon->d6;           /* Set in precalc_u2()               */
  memset(sum, 0, window->sgsizS * sizeof(double));
  memset(adjust, 0, window->sgsizS * sizeof(double));
  memset(wtot, 0, window->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column.              */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    midy[cs] = wincon->mdy[cs];
    sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	midy[cs] = wincon->mdy[cs];
	sum[cs] += windat->u2[c] * windat->dzu2[c] * midy[cs];
      }
    }
  }

  /* Estimate the velocity adjustment based on no weighting          */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    cs = window->w3_e2[cc];
    adjust[cs] =
      (windat->u2flux[cs] - sum[cs] * windat->dt * window->h1au2[cs]) /
      (depth[cs] * midy[cs] * windat->dt * window->h1au2[cs]);
  }
  /* Get the maximum and minimum absolute velocities differences      */
  memset(vmax, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    vmin[c] = wincon->velmax;
  }
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    d1 = fabs(windat->u2[c] + adjust[cs]);
    vmax[cs] = (d1 > vmax[cs]) ? d1 : vmax[cs];
    vmin[cs] = (d1 < vmin[cs]) ? d1 : vmin[cs];
  }
  /* Get the weights                                                 */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    d1 = fabs(windat->u2[c] + adjust[cs]);
    d2 = vmin[cs] - vmax[cs];
    w[c] = (d2) ? (d1 - vmax[cs]) / d2 : 1.0;
    wtot[cs] += w[c] * windat->dzu2[c] * midy[cs];
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    adjust[cs] =
      (windat->u2flux[cs] - sum[cs] * windat->dt * window->h1au2[cs]) / 
      (wtot[cs] * windat->dt * window->h1au2[cs]);
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no2_e2; cc++) {
	cs = open[bn]->obc_e2[cc];
	depth[cs] -= window->botzu2[cs];
	adjust[cs] =
	  (windat->u2flux[cs] - sum[cs] * windat->dt * window->h1au2[cs]) /
	  (wtot[cs] * windat->dt * window->h1au2[cs]);
      }
    }
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs = window->m2d[c];
    windat->u2[c] += (adjust[cs] * w[c]);
  }
  if(!wincon->dobdry_u2) {
    for (bn = 0; bn < window->nobc; bn++) {
      for (cc = 1; cc <= open[bn]->no3_e2; cc++) {
	c = open[bn]->obc_e2[cc];
	cs = window->m2d[c];
	windat->u2[c] += (adjust[cs] * w[c]);
      }
    }
  }

  /* Set the velocity at dry water columns equal to zero             */
  for (cc=1; cc <= wincon->ncdry_e2; cc++) {
    c = cs = wincon->cdry_e2[cc];
    windat->u2av[cs] = 0.0;
    windat->u2flux[cs] = 0.0;
    while (c != window->zm1[c]) {
      windat->u2[c] = 0.0;
      c = window->zm1[c];
    }
  }

  /* Save the bottom velocity for the 2D mode                        */
  memset(windat->u2bot, 0, window->sgsizS * sizeof(double));
  for (cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];         /* 3D surface coordinate             */
    cs = window->m2d[c];        /* 2D coordinate                     */
    cb = wincon->i6[cc];        /* 3D bottom coordinate              */
    windat->u2bot[cs] = windat->u2b[cb] - windat->u2avb[cs];
  }

  /* Set the velocity for cells one cell deep                        */
  for (cc = 1; cc <= wincon->ncbot_e2; cc++) {
    c = wincon->cbot_e2[cc];
    cs = window->m2d[c];
    windat->u2[c] = windat->u2av[cs];
    windat->u2bot[cs] = 0.0;
  }
}

/* END velocity_adjust_weighted()                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the boundary  */
/*-------------------------------------------------------------------*/
void set_sponge_m(geometry_t *geom, /* Global geometry             */
		  master_t *master  /* Master data structure */
  )
{
  open_bdrys_t *open;           /* Open boudary data structure */
  int c, cb, cc;                /* Sparse counters */
  double scale = 4.0;           /* Max bottom friction multiplier */
  double sfact = 0.9;           /* Safety factor */
  int bn, n;                    /* Boundary counter */
  double rm, spgn;              /* Dummies */

  for (n = 0; n < geom->nobc; ++n) {
    open = geom->open[n];

    /*---------------------------------------------------------------*/
    /* Increase in bottom friction in the sponge */
    if (open->sponge_zone) {
      for (cc = 1; cc <= open->no2_t; cc++) {
        c = open->obc_t[cc];

        spgn = open->sponge_zone;
        rm = scale * master->Cd[c];
        for (bn = 1; bn <= spgn; bn++) {
          master->Cd[c] = (rm - master->Cd[c]) * (double)(bn - spgn) /
            (double)(1 - spgn) + master->Cd[c];
          c = open->nmap[c];
        }
      }
    }

    /*---------------------------------------------------------------*/
    /* Increase in horizontal friction in the sponge. Note if the    */
    /* normal component has Smagorinsky diffusion then no sponge     */
    /* is set since viscosity is set every timestep.                 */
    if (open->sponge_zone_h) {
      if (master->u1vh != master->sdc) {
	for (cc = 1; cc <= open->no2_e1; cc++) {
	  c = cb = open->obc_e1[cc];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    if (open->sponge_f)
	      rm = open->sponge_f * master->u1vh[c];
	    else {
	      rm = 1.0 / (geom->h1au1[c] * geom->h1au1[c]) +
		1.0 / (geom->h2au1[c] * geom->h2au1[c]);
	      rm = sfact / (4.0 * rm * master->grid_dt);
	    }
	    master->u1vh[c] =
	      (fabs(master->u1vh[c]) - rm) * (double)(bn - 1) / 
	      (double)spgn + rm;
	    c = open->nmap[c];
	  }
	  if (open->ocodex & R_EDGE)
	    master->u1vh[geom->xp1[cb]] = master->u1vh[cb];
	  if (open->ocodex & L_EDGE)
	    master->u1vh[geom->xm1[cb]] = master->u1vh[cb];
	}
      }
      if (master->u2vh != master->sdc) {
	for (cc = 1; cc <= open->no2_e2; cc++) {
	  c = cb = open->obc_e2[cc];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    if (open->sponge_f)
	      rm = open->sponge_f * master->u2vh[c];
	    else {
	      rm = 1.0 / (geom->h1au2[c] * geom->h1au2[c]) +
		1.0 / (geom->h2au2[c] * geom->h2au2[c]);
	      rm = sfact / (4.0 * rm * master->grid_dt);
	    }
	    master->u2vh[c] = 
	      (fabs(master->u2vh[c]) - rm) * (double)(bn - 1) / 
	      (double)spgn + rm;
	    c = open->nmap[c];
	  }
	  if (open->ocodey & F_EDGE)
	    master->u2vh[geom->yp1[cb]] = master->u2vh[cb];
	  if (open->ocodey & B_EDGE)
	    master->u2vh[geom->ym1[cb]] = master->u2vh[cb];
	}
      }
    }
  }
}

/* END set_sponge_m()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the boundary  */
/*-------------------------------------------------------------------*/
void set_sponge_w(geometry_t *window,   /* Window geometry           */
		  window_t *windat,     /* Window data               */
		  win_priv_t *wincon    /* Window constants          */
  )
{
  open_bdrys_t *open;           /* Open boudary data structure       */
  int c, cs, cc;                /* Sparse counters                   */
  double sfact = 0.9;           /* Safety factor                     */
  double tfact = 0.3;           /* tanh() multiplier : small = long  */
  int bn, ln, lc, sc, n;        /* Boundary counter                  */
  double rm, spgn, vh;          /* Dummies                           */
  double *mask = wincon->w1;
  int blendf = 1;

  memset(mask, 0, window->sgsiz * sizeof(double));
  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];

    /*---------------------------------------------------------------*/
    /* Increase in horizontal friction in the sponge.                */
    if (open->sponge_zone_h) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	cs = window->m2d[c];
	spgn = open->sponge_zone_h;
	for (bn = 1; bn <= spgn; bn++) {
	  if (open->sponge_f)
	    rm = open->sponge_f * wincon->u1vh[c];
	  else {
	    rm = 1.0 / (window->h1au1[cs] * window->h1au1[cs]) +
	      1.0 / (window->h2au1[cs] * window->h2au1[cs]);
	    rm = sfact / (4.0 * rm * windat->dt);
	  }
	  vh = wincon->u1vh[c];
	  /* Linear ramp to rm on the boundary                      */
	  wincon->u1vh[c] = (fabs(wincon->u1vh[c]) - rm) * 
	    (double)(bn - 1) / (double)spgn + rm;
	  /* Hyperbolic tangent ramp to rm on the boundary          */

	  wincon->u1vh[c] =  vh + (rm - vh) * 
	    (1 - tanh(tfact * (double)(bn - 1)));

	  mask[c] = 1;
	  c = open->nmap[c];
	}
      }

      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	cs = window->m2d[c];
	spgn = open->sponge_zone_h;
	for (bn = 1; bn <= spgn; bn++) {
	  if (open->sponge_f)
	    rm = open->sponge_f * wincon->u2vh[c];
	  else {
	    rm = 1.0 / (window->h1au2[cs] * window->h1au2[cs]) +
	      1.0 / (window->h2au2[cs] * window->h2au2[cs]);
	    rm = sfact / (4.0 * rm * windat->dt);
	  }
	  vh = wincon->u2vh[c];
	  /* Linear ramp to rm on the boundary                      */
	  wincon->u2vh[c] = (fabs(wincon->u2vh[c]) - rm) * 
	    (double)(bn - 1) / (double)spgn + rm;
	  /* Hyperbolic tangent ramp to rm on the boundary          */

	  wincon->u2vh[c] =  vh + (rm - vh) * 
	    (1 - tanh(tfact * (double)(bn - 1)));

	  mask[c] = 1;
	  c = open->nmap[c];
	}
      }

      /* Blend the sponge laterally to remove large discontinuties  */
      /* in viscosity.                                              */
      if (blendf) {
	for (cc = 1; cc <= open->no3_e1; cc++) {
	  c = open->obc_e1[cc];
	  cs = window->m2d[c];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    vh = wincon->u1vh[c];
	    /* Blend u1vh to the front (yp1) of current position    */
	    lc = sc = window->yp1[c];
	    if (!mask[lc]) {
	      /* Get viscosity spgn cells from the current cell     */
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->yp1[lc];
	      rm = wincon->u1vh[lc];
	      /* Blend to this viscosity                            */
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		wincon->u1vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		wincon->u1vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		mask[sc] = 1;
		sc = window->yp1[sc];
	      }
	    }
	    /* Blend u2vh to the back (ym1) of current position     */
	    lc = sc = window->ym1[c];
	    if (!mask[lc]) {
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->ym1[lc];
	      rm = wincon->u1vh[lc];
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		wincon->u1vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		wincon->u1vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		mask[sc] = 1;
		sc = window->ym1[sc];
	      }
	    }
	    c = open->nmap[c];
	  }
	}
	for (cc = 1; cc <= open->no3_e2; cc++) {
	  c = open->obc_e2[cc];
	  cs = window->m2d[c];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    vh = wincon->u2vh[c];
	    /* Blend u2vh to the right (xp1) of current position    */
	    lc = sc = window->xp1[c];
	    if (!mask[lc]) {
	      /* Get viscosity spgn cells from the current cell     */
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->xp1[lc];
	      rm = wincon->u2vh[lc];
	      /* Blend to this viscosity                            */
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		wincon->u2vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		wincon->u2vh[sc] = rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		mask[sc] = 1;
		sc = window->xp1[sc];
	      }
	    }
	    /* Blend u2vh to the left (xm1) of current position     */
	    lc = sc = window->xm1[c];
	    if (!mask[lc]) {
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->xm1[lc];
	      rm = wincon->u2vh[lc];
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		wincon->u2vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		wincon->u2vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		mask[sc] = 1;
		sc = window->xm1[sc];
	      }
	    }
	    c = open->nmap[c];
	  }
	}
      }
    }

  }
}

/* END set_sponge_w()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the boundary  */
/*-------------------------------------------------------------------*/
void set_sponge(geometry_t *window,   /* Window geometry             */
		double *u1vh,         /* e1 viscosity                */
		double *u2vh,         /* e2 viscosity                */
		double dt,            /* 3D time step                */
		double *mask          /* 3D dummy array              */
  )
{
  open_bdrys_t *open;           /* Open boudary data structure       */
  int c, cs, cb, cc;            /* Sparse counters                   */
  double sfact = 0.9;           /* Safety factor                     */
  double tfact = 0.3;           /* tanh() multiplier : small = long  */
  int bn, ln, lc, sc, n;        /* Boundary counter                  */
  double rm, spgn, vh;          /* Dummies                           */
  int blendf = 2;

  /*-----------------------------------------------------------------*/
  /* Increase in horizontal friction in the sponge: u1vh             */
  memset(mask, 0, window->sgsiz * sizeof(double));
  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];

    if (open->sponge_zone_h) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = cb = open->obc_e1[cc];
	cs = window->m2d[c];
	spgn = open->sponge_zone_h;
	for (bn = 1; bn <= spgn; bn++) {
	  if (open->sponge_f)
	    rm = open->sponge_f * u1vh[c];
	  else {
	    rm = 1.0 / (window->h1au1[cs] * window->h1au1[cs]) +
	      1.0 / (window->h2au1[cs] * window->h2au1[cs]);
	    rm = sfact / (4.0 * rm * dt);
	  }
	  vh = u1vh[c];
	  /* Linear ramp to rm on the boundary                      */
	  u1vh[c] = (fabs(u1vh[c]) - rm) * 
	    (double)(bn - 1) / (double)spgn + rm;
	  /* Hyperbolic tangent ramp to rm on the boundary          */

	  u1vh[c] =  vh + (rm - vh) * 
	    (1 - tanh(tfact * (double)(bn - 1)));

	  if (open->options & OP_ISPNG) u2vh[c] = u1vh[c];
	  mask[c] = 1;
	  c = open->nmap[c];
	}
	if (open->ocodex & R_EDGE)
	  u1vh[window->xp1[cb]] = u1vh[cb];
	if (open->ocodex & L_EDGE)
	  u1vh[window->xm1[cb]] = u1vh[cb];
      }

      /*------------------------------------------------------------*/
      /* Blend the sponge laterally to remove large discontinuties  */
      /* in viscosity.                                              */
      if (blendf == 1) {
	for (cc = 1; cc <= open->no3_e1; cc++) {
	  c = open->obc_e1[cc];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    vh = u1vh[c];
	    /* Blend u1vh to the front (yp1) of current position    */
	    lc = sc = window->yp1[c];
	    if (!mask[lc]) {
	      /* Get viscosity spgn cells from the current cell     */
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->yp1[lc];
	      rm = u1vh[lc];
	      /* Blend to this viscosity                            */
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		u1vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		u1vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		if (open->options & OP_ISPNG) u2vh[sc] = u1vh[sc];
		mask[sc] = 1;
		sc = window->yp1[sc];
	      }
	    }
	    /* Blend u1vh to the back (ym1) of current position     */
	    lc = sc = window->ym1[c];
	    if (!mask[lc]) {
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->ym1[lc];
	      rm = u1vh[lc];
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		u1vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		u1vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		if (open->options & OP_ISPNG) u2vh[sc] = u1vh[sc];
		mask[sc] = 1;
		sc = window->ym1[sc];
	      }
	    }
	    c = open->nmap[c];
	  }
	} 
      } else if (blendf == 2) { /* Unidirectional blending           */
	for (cc = 1; cc <= open->no3_e1; cc++) {
	  c = open->obc_e1[cc];
	  spgn = open->sponge_zone_h;
	  blend_var(open, u1vh, u2vh, mask, spgn, c);
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Increase in horizontal friction in the sponge: u2vh             */
  memset(mask, 0, window->sgsiz * sizeof(double));
  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];

    if (open->sponge_zone_h) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
      c = cb = open->obc_e2[cc];
	cs = window->m2d[c];
	spgn = open->sponge_zone_h;
	for (bn = 1; bn <= spgn; bn++) {
	  if (open->sponge_f)
	    rm = open->sponge_f * u2vh[c];
	  else {
	    rm = 1.0 / (window->h1au2[cs] * window->h1au2[cs]) +
	      1.0 / (window->h2au2[cs] * window->h2au2[cs]);
	    rm = sfact / (4.0 * rm * dt);
	  }
	  vh = u2vh[c];
	  /* Linear ramp to rm on the boundary                      */
	  u2vh[c] = (fabs(u2vh[c]) - rm) * 
	    (double)(bn - 1) / (double)spgn + rm;
	  /* Hyperbolic tangent ramp to rm on the boundary          */

	  u2vh[c] =  vh + (rm - vh) * 
	    (1 - tanh(tfact * (double)(bn - 1)));

	  if (open->options & OP_ISPNG) u1vh[c] = u2vh[c];
	  mask[c] = 1;
	  c = open->nmap[c];
	}
	if (open->ocodey & F_EDGE)
	  u2vh[window->yp1[cb]] = u2vh[cb];
	if (open->ocodey & B_EDGE)
	  u2vh[window->ym1[cb]] = u2vh[cb];
      }

      /*------------------------------------------------------------*/
      /* Blend the sponge laterally to remove large discontinuties  */
      /* in viscosity.                                              */
      if (blendf == 1) {
	for (cc = 1; cc <= open->no3_e2; cc++) {
	  c = open->obc_e2[cc];
	  spgn = open->sponge_zone_h;
	  for (bn = 1; bn <= spgn; bn++) {
	    vh = u2vh[c];
	    /* Blend u2vh to the right (xp1) of current position    */
	    lc = sc = window->xp1[c];
	    if (!mask[lc]) {
	      /* Get viscosity spgn cells from the current cell     */
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->xp1[lc];
	      rm = u2vh[lc];
	      /* Blend to this viscosity                            */
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		u2vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		u2vh[sc] = rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		if (open->options & OP_ISPNG) u1vh[sc] = u2vh[sc];
		mask[sc] = 1;
		sc = window->xp1[sc];
	      }
	    }
	    /* Blend u2vh to the left (xm1) of current position     */
	    lc = sc = window->xm1[c];
	    if (!mask[lc]) {
	      for (ln = 2; ln <= spgn && !mask[lc]; ln++)
		lc = window->xm1[lc];
	      rm = u2vh[lc];
	      for (ln = 2; ln <= spgn; ln++) {
		/* Linear                                           */
		u2vh[sc] = (vh - rm) * (double)(ln - 1) /
		  (double)(1 - spgn) + vh;
		/* Hyperbolic tangent                               */
		/*
		u2vh[sc] =  rm + (vh - rm) * 
		  (1 - tanh(tfact * (double)(ln - 1)));
		*/
		if (open->options & OP_ISPNG) u1vh[sc] = u2vh[sc];
		mask[sc] = 1;
		sc = window->xm1[sc];
	      }
	    }
	    c = open->nmap[c];
	  }
	}
      }  else if (blendf == 2) { /* Unidirectional blending           */
	for (cc = 1; cc <= open->no3_e2; cc++) {
	  c = open->obc_e2[cc];
	  spgn = open->sponge_zone_h;
	  blend_var(open, u2vh, u1vh, mask, spgn, c);
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
  int cc, c, cs;                        /* Counters, cell locations  */

  if (windat->u1vhin == NULL || windat->u2vhin == NULL) return;

  for (n = 0; n < window->nobc; ++n) {
    open = window->open[n];

    if (open->sponge_zone_h) {
      spgn = open->sponge_zone_h;
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	cs = window->m2d[c];
	for (bn = 1; bn <= spgn; bn++) {
	  if (wincon->sue1 > 0.0 && wincon->u1vh0 < 0.0) {
	    wincon->u1vh[c] = wincon->sue1 * windat->sdc[c];
	    if (wincon->bsue1) { 
	      if (wincon->diff_scale == LINEAR)
		wincon->u1vh[c] += fabs(wincon->bsue1 * window->h1au1[cs] / wincon->hmean1);
	      if (wincon->diff_scale == NONLIN)
		wincon->u1vh[c] += fabs(wincon->bsue1 * window->h1au1[cs] * window->h1au1[cs] / 
					(wincon->hmean1 * wincon->hmean1));
	    }
	  } else {
	    wincon->u1vh[c] = windat->u1vhin[cs];
	  }
	  c = open->nmap[c];
	  cs = window->m2d[c];
	}
      }
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	cs = window->m2d[c];
	for (bn = 1; bn <= spgn; bn++) {
	  if (wincon->sue2 > 0.0 && wincon->u2vh0 < 0.0) {
	    wincon->u2vh[c] = wincon->sue2 * windat->sdc[c];
	    if (wincon->bsue2) { 
	      if (wincon->diff_scale == LINEAR)
		wincon->u2vh[c] += fabs(wincon->bsue2 * window->h2au2[cs] / wincon->hmean2);
	      if (wincon->diff_scale == NONLIN)
		wincon->u2vh[c] += fabs(wincon->bsue2 * window->h2au2[cs] * window->h2au2[cs] / 
					(wincon->hmean2 * wincon->hmean2));
	    }
	  } else
	    wincon->u2vh[c] = windat->u2vhin[cs];
	  c = open->nmap[c];
	  cs = window->m2d[c];
	}
      }
    }
  }
}

/* END reset_sponge_zone()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to a variable blend unidirectionally                      */
/*-------------------------------------------------------------------*/
void blend_var(open_bdrys_t *open,
	       double *var,
	       double *var2,
	       double *mask,
	       int spgn,
	       int c
	       )
{
  int lc, bn;


  lc = c;
  get_blend_val(var, lc, spgn, open->tmpp, mask);
  if (open->options & OP_ISPNG)
    get_blend_val(var2, lc, spgn, open->tmpp, mask);
  if (!mask[open->tmpp[c]]) {
    for (bn = 1; bn <= spgn; bn++) {
      lc = open->tmpp[lc];
      get_blend_val(var, lc, spgn, open->nmap, mask);
      get_blend_val(var, lc, spgn, open->omap, mask);
      if (open->options & OP_ISPNG) {
	get_blend_val(var2, lc, spgn, open->nmap, mask);
	get_blend_val(var2, lc, spgn, open->omap, mask);
      }
      mask[lc] = 1;
    }
  }

  get_blend_val(var, lc, spgn, open->tmpm, mask);
  if (open->options & OP_ISPNG)
    get_blend_val(var2, lc, spgn, open->tmpm, mask);
  if (!mask[open->tmpm[c]]) {
    for (bn = 2; bn <= spgn; bn++) {
      lc = open->tmpm[lc];
      get_blend_val(var, lc, spgn, open->nmap, mask);
      get_blend_val(var, lc, spgn, open->omap, mask);
      if (open->options & OP_ISPNG) {
	get_blend_val(var2, lc, spgn, open->nmap, mask);
	get_blend_val(var2, lc, spgn, open->omap, mask);
      }
      mask[lc] = 1;
    }
  }
}

/* END blend_var()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the blended values                                       */
/*-------------------------------------------------------------------*/
void get_blend_val(double *val, int c, int spgn, int *map, double *mask)
{
  int ln, lc;
  int bflag = 1;
  double vh, rm;
  double tfact = 0.3;           /* tanh() multiplier : small = long  */

  vh = val[c];
  lc = map[c];

  if (!mask[lc]) {
    for (ln = 1; ln <= spgn && !mask[lc]; ln++)
      lc = map[lc];
    rm = val[lc];
    lc = map[c];
    for (ln = 2; ln <= spgn; ln++) {
      if (bflag) {
	/* Linear                                                    */
	val[lc] = (vh - rm) * (double)(ln - 1) /
	  (double)(1 - spgn) + vh;
      } else {
	/* Hyperbolic tangent                                        */
	val[lc] = rm + (vh - rm) * 
	  (1 - tanh(tfact * (double)(ln - 1)));
      }
      lc = map[lc];
    }
  }  
}

/* END get_blend_val()                                               */
/*-------------------------------------------------------------------*/

