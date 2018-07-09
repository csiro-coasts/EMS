/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/mommix.c
 *  
 *  Description: Horizontal mixing of momentum
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: mommix.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void scale_hdiff(geometry_t *window, double *AH, double smag, int mode);
void scale_hdiff_m(master_t *master, double *AH, double smag, double AH0, int mode);

/*-------------------------------------------------------------------*/
/* This routine reads the parameter file to determine the            */
/* horizontal mixing method to be used, and sets the appropriate     */
/* pointers.                                                         */
/*-------------------------------------------------------------------*/
void hvisc_init(master_t *master,    /* Master data                  */
		win_priv_t **wincon  /* Window private data          */
  )
{
  int i, n;

  /* Set the horizontal mixing method                              */
  for (n = 1; n <= master->nwindows; n++) {
    wincon[n]->hor_mix = (mix_method_t *)malloc(sizeof(mix_method_t));
    wincon[n]->hor_mix->pre = hvisc_setup_pre;
    wincon[n]->hor_mix->setup = hvisc_setup;
    wincon[n]->hor_mix->setup_av = hvisc_setup_2d;
    if (wincon[n]->visc_method & PRE794) {
      wincon[n]->hor_mix->pre = hvisc_null;
      wincon[n]->hor_mix->setup = hvisc_setup_old;
      wincon[n]->hor_mix->setup_av = hvisc_setup_2d_old;
    }
    if (wincon[n]->visc_method & LAPLACIAN) {
      wincon[n]->hor_mix->u1 = hvisc_u1_3d;
      wincon[n]->hor_mix->u2 = hvisc_u2_3d;
      wincon[n]->hor_mix->u1av = hvisc_u1_2d;
      wincon[n]->hor_mix->u2av = hvisc_u2_2d;
    } else if (wincon[n]->visc_method & SIMPLE) {
      if (wincon[n]->u1vh0 >= 0.0 && wincon[n]->u2vh0 >= 0.0) {
	wincon[n]->hor_mix->setup = hvisc_null;
	wincon[n]->hor_mix->setup_av = hvisc_2d_null;
      }
      wincon[n]->hor_mix->u1 = hvisc_u1_3d_simple;
      wincon[n]->hor_mix->u2 = hvisc_u2_3d_simple;
      wincon[n]->hor_mix->u1av = hvisc_u1_2d_simple;
      wincon[n]->hor_mix->u2av = hvisc_u2_2d_simple;
    }
    if (wincon[n]->u1_f & HDIFF)
      wincon[n]->hor_mix->u1 = hvisc_null;
    if (wincon[n]->u2_f & HDIFF)
      wincon[n]->hor_mix->u2 = hvisc_null;
    if (wincon[n]->u1av_f & HDIFF)
      wincon[n]->hor_mix->u1av = hvisc_2d_null;
    if (wincon[n]->u2av_f & HDIFF)
      wincon[n]->hor_mix->u2av = hvisc_2d_null;
  }
}

/* END hvisc_init()                                                  */
/*-------------------------------------------------------------------*/


void hvisc_null(geometry_t *window,   /* Window geometry             */
		window_t *windat,     /* Window data                 */
		win_priv_t *wincon    /* Window constants            */
		)
{
}

void hvisc_2d_null(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   double *tzp           /* Pointer to surface level */
		   )
{
}


/*-------------------------------------------------------------------*/
/* Computs the stress tensors and Smagorinsky diffusion. This is     */
/* done before the momentum calculations so that the Smagorinsky     */
/* values can be transferred, rather than computed over the whole    */
/* wet+auxilliary+ghost array. The latter method can cause problems  */
/* using the t12[window->xp1[yp1]] value.                            */
/* Blumberg and Herring (1987) multiply the first term in t11 by h2. */
/* This seems to be in error, as for uniform h the units of the      */
/* mixing term don't equate to ms-2. Also Apel (1987) derives the    */
/* mixing tensors in cartesian coordinates and multiplies the first  */
/* term in the stress tensor by 2, rather than h2. This latter       */
/* formulation is considered correct, and is formulated in this      */
/* routine. The formulation of Blumberg & Herring, without the       */
/* multiplication by h2 so that units are consistent, is retained in */
/* the routine hvisc_setup_old().                                    */
/*-------------------------------------------------------------------*/
void hvisc_setup_pre(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
  )
{
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface sparse coordinates          */
  int zm1;                    /* 3D sparse coordinate at k-1         */
  int xp1, xm1;               /* 3D sparse coordinate at i+1, i-1    */
  int yp1, ym1;               /* 3D sparse coordinate at j+1, j-1    */
  int xp1s, xm1s;             /* 2D sparse coordinate at i+1, i-1    */
  int yp1s, ym1s;             /* 2D sparse coordinate at j+1, j-1    */
  int xmym1s;                 /* Sparse coordinate at (i-1,j-1)      */
  int xmyp1s;                 /* Sparse coordinate at (i-1,j+1)      */
  int xmyp1;                  /* Sparse coordinate at (i-1,j+1)      */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  double c1;                  /* Constant for x diffusion            */
  double c2;                  /* Constant for y diffusion            */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(t11, 0, window->sgsiz * sizeof(double));
  memset(t12, 0, window->sgsiz * sizeof(double));
  memset(t22, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate the stress tensors.                                   */
  /* This is performed over all sparse cells (wet, auxiliary and     */
  /* ghost). The tensors t11 and t22 are cell centered, t12 is cell  */
  /* cornered thus the centered bottom coordinate must be used (this */
  /* is different to the face centered e1 bottom coordinate on       */
  /* western edges due to the stagger of the velocity cells).        */
  /* Auxiliary cells are not defined for cell centers on western and */
  /* southern boundaries, so the _t vectors cannot be used here.     */
  /* Simply loop over the entire sparse vector.                      */
  for (cc = 1; cc <= window->enonS; cc++) {
    c = cc;
    zm1 = window->zm1[c];
    
    while (c != zm1) {

      cs = window->m2d[c];
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      xm1 = window->xm1[c];
      ym1 = window->ym1[c];

      xp1s = window->xp1[cs];
      yp1s = window->yp1[cs];
      xm1s = window->xm1[cs];
      ym1s = window->ym1[cs];
      xmym1s = window->ym1[xm1s];

      /* Save the (x,x) direction stress in t11.                     */
      /* Cell centered. Note: factor of 2 cancels 0.5 from u1b mean  */
      /* in the 2nd term. Also mulptiply all terms by Ds here for    */
      /* sigma.                                                      */
      t11[c] = wincon->Ds[cs] *
        (2.0*(windat->u1b[xp1] - windat->u1b[c]) / window->h1acell[cs] +
         (windat->u2b[yp1] + windat->u2b[c]) *
         (window->h1au2[yp1s] - window->h1au2[cs]) / window->cellarea[cs]);

      /* Save the (y,y) direction stress in t22.                     */
      /* Cell centered.                                              */
      t22[c] = wincon->Ds[cs] *
        (2.0*(windat->u2b[yp1] - windat->u2b[c]) / window->h2acell[cs] +
	 (windat->u1b[xp1] + windat->u1b[c]) *
         (window->h2au1[xp1s] - window->h2au1[cs]) / window->cellarea[cs]);

      /* Save the (x,y) direction stress in t12                      */
      /* Cell cornered.                                              */
      c1 = 0.25 * (wincon->Ds[cs] + wincon->Ds[xm1s] + wincon->Ds[ym1s] +
                   wincon->Ds[xmym1s]);

      t12[c] = c1 * (windat->u1b[c] / window->h1au1[cs] -
                     windat->u1b[ym1] / window->h1au1[ym1s]) *
        (window->h1au1[cs] + window->h1au1[ym1s]) /
        (window->h2au1[cs] + window->h2au1[ym1s]);

      t12[c] += c1 * (windat->u2b[c] / window->h2au2[cs] -
                      windat->u2b[xm1] / window->h2au2[xm1s]) *
        (window->h2au2[cs] + window->h2au2[xm1s]) /
        (window->h1au2[cs] + window->h1au2[xm1s]);
      c = zm1;
      zm1 = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the Smagorinsky diffusivity if required                     */
  if (wincon->smagorinsky != 0.0) {
    if (!windat->sdc)
      hd_quit("Smagorinsky diffusion requires tracer called smag\n");
    memset(windat->sdc, 0, window->sgsiz*sizeof(double));
    /* Calculate the Smagorinsky diffusion coefficient               */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];

      zm1 = window->zm1[c];
      
      while (c != zm1) {

	xp1 = window->xp1[c];
	yp1 = window->yp1[c];
	cs = window->m2d[c];
	c1 = 2.0 * wincon->Ds[cs] * wincon->Ds[cs];
	
	/* Cell centered diffusivity                                   */
	c2 = 0.25 * (t12[c] + t12[xp1] + t12[yp1] + t12[window->xp1[yp1]]);
	windat->sdc[c] = sqrt(t11[c] * t11[c] / c1 +
			      c2 * c2 / c1 +
			      t22[c] * t22[c] / c1) *
	  window->h1acell[cs] * window->h2acell[cs] * wincon->smagorinsky;
  	c = zm1;
	zm1 = window->zm1[c];
      }
    }
  }
}

/* END hvisc_setup_pre()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Implements smoothing, local increases and sponge zones on         */
/* horizontal mixing. Note; smoothing and local increases will not   */
/* result in identical field for 1 and multiple windows. To do this  */
/* sdc should be gathered, the smoothing/increases performed on the  */
/* master and then scattered. Smoothing with a 9 point kernal will   */
/* require an extra transfer after this routine is invoked.          */
/*-------------------------------------------------------------------*/
void hvisc_setup(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface sparse coordinates          */
  int sf = 0;                 /* Sponge flag                         */

  /*-----------------------------------------------------------------*/
  /* Set the Smagorinsky diffusivity if required                     */
  if (wincon->smagorinsky != 0.0) {

    /* Adjust the Smagorinsky coefficient for active alert actions   */
    /* and horizontal sponge zones.                                  */

    /* If smagcode & U1_A|U2_A|U1_AK|U2_AK then the relevant mixing  */
    /* array is allocated and filled with (grid weighted) constant   */
    /* coefficients in hd_init.c, and this code is not invoked.      */

    /* If smagcode & U1_SP|U2_SP|U1_SPK|U2_SPK then the relevant     */
    /* mixing array points to sdc, set in windows.c.                 */

    /* Set the flag for re-defining sponge zones                     */
    if (wincon->smagcode & (U1_SA|U2_SA)) sf = 1;
    if (wincon->smagcode & (U1_A|U2_A)) {
      if (wincon->alertf & (LEVEL1_U|LEVEL1_V|LEVEL1_UA|LEVEL1_VA))
	sf = 1;
    }

    Set_lateral_BC_density_w(windat->sdc, window->nbpt, window->bpt,
                             window->bin);

    /* Reset viscosity for U1_A/U2_A prior to shear alerts           */
    if (wincon->u1vh[0] == 1) {
      set_hdiff(window, wincon->u1vh, wincon->u1vh0, 0);
      sf = 1;
    }
    if (wincon->u2vh[0] == 1) {
      set_hdiff(window, wincon->u2vh, wincon->u2vh0, 1);
      sf = 1;
    }
    wincon->u1vh[0] = wincon->u2vh[0] = 0;

    /* Copy Smagorinsky diffusion to horizontal diffusion. This is   */
    /* not subject to LEVEL1 alerts or sponge zones.                 */
    if (wincon->smagcode & U1_SAK) {
      memcpy(wincon->u1kh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u1kh, wincon->kue1, 3);
    }
    if (wincon->smagcode & U2_SAK) {
      memcpy(wincon->u2kh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u2kh, wincon->kue2, 4);
    }

    /* Copy Smagorinsky diffusion to horizontal viscosity. This is   */
    /* subject to LEVEL1 alerts and separate sponge zones in the e1  */
    /* and e2 directions.                                            */
    /* Uncomment this code if separate sponges apply in the e1 and   */
    /* e2 directions.                                                */
    if (wincon->smagcode & U1_SA) {
      memcpy(wincon->u1vh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u1vh, wincon->sue1, 1);
    }
    if (wincon->smagcode & U2_SA) {
      memcpy(wincon->u2vh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u2vh, wincon->sue2, 2);
    }
    /* Reset diffusion throughout the water column for active alerts */
    /* Alerts at 3d locations.                                       */
    if (wincon->alertf & LEVEL1_U) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_U;
      /*dump_snapshot(sched_get_even_by_name(schedule, "dumps"), windat->t);*/
    }
    if (wincon->alertf & LEVEL1_V) {
      if (wincon->smagcode & U2_A) wincon->u2vh[0] = 1;
      for(cc = 1; (c = windat->cu2[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_V;
      /*dump_snapshot(sched_get_even_by_name(schedule, "dumps"), windat->t);*/
    }
    if (wincon->alertf & LEVEL1_W) {
      c = windat->cw;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL1_W;
    }
    if (wincon->alertf & LEVEL2w) {
      c = windat->cdw;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL2w;
    }
    /* Alerts at 2d locations.                                       */
    if (wincon->alertf & LEVEL1_UA) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1a[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_UA;
    }
    if (wincon->alertf & LEVEL1_VA) {
      if (wincon->smagcode & U2_A) wincon->u2vh[0] = 1;
      for(cc = 1; (c = windat->cu2a[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_VA;
    }
    if (wincon->alertf & LEVEL2) {
      c = windat->cdeta;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL2;
    }

    /* Increase diffusion if required */
    /*
    if (wincon->alertf & LEVEL1_U1VH) {
      for (cc = 1; cc <= window->n3_e1; cc++) {
	c = window->w3_e1[cc];
	wincon->u1vh[c] = min(-wincon->u1vh0, wincon->u1vh[c] / wincon->smagorinsky);
      }
      wincon->alertf &= ~LEVEL1_U1VH;
    }
    if (wincon->alertf & LEVEL1_U2VH) {
      for (cc = 1; cc <= window->n3_e2; cc++) {
	c = window->w3_e2[cc];
	wincon->u2vh[c] = min(-wincon->u2vh0, wincon->u2vh[c] / wincon->smagorinsky);
      }
      wincon->alertf &= ~LEVEL1_U2VH;
    }
    */

    /* Set the sponge zone                                           */
    /*set_sponge_w(window, windat, wincon);*/
    if (sf) {
      reset_sponge_zone(window);
      set_sponge(window, wincon->u1vh, wincon->u2vh, windat->dt,  wincon->w1);
    }

    /* Ensure the upper limit of the smagorinsky diffusion is the    */
    /* constant viscosity, and smooth.                               */
    if (wincon->smag_smooth && wincon->smagorinsky > 0.0) {
      if (wincon->u1vh0 < 0.0)
	smooth_w(window, wincon->u1vh, wincon->w6, window->w3_e1,
		 window->n3_e1, wincon->smag_smooth, -wincon->u1vh0);
      if (wincon->u2vh0 < 0.0)
	smooth_w(window, wincon->u2vh, wincon->w6, window->w3_e2,
		 window->n3_e2, wincon->smag_smooth, -wincon->u2vh0);
      if (wincon->u1kh0 < 0.0)
	smooth_w(window, wincon->u1kh, wincon->w6, window->w3_t,
		 window->n3_t, wincon->smag_smooth, -wincon->u1kh0);
      if (wincon->u2kh0 < 0.0)
	smooth_w(window, wincon->u2kh, wincon->w6, window->w3_t,
		 window->n3_t, wincon->smag_smooth, -wincon->u2kh0);
    }
  }

  /* Create a sponge zone on coarse-fine boundaries for zoomed grids */
  if (wincon->dozoom == DOZOOM) 
    set_sponge_zoom(window, windat, wincon, wincon->w2, wincon->w3,
		    window->nzbe1, window->nzbe2);

}

/* END hvisc_setup()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Original pre v794 formulation : factor of 2 missing on stresses.  */
/* Sets the tensors, Smagorinsky coefficient (if required) and       */
/* sponges for 3D horizontal mixing of momentum.                     */
/*-------------------------------------------------------------------*/
void hvisc_setup_old(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
  )
{
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface sparse coordinates          */
  int zm1;                    /* 3D sparse coordinate at k-1         */
  int xp1, xm1;               /* 3D sparse coordinate at i+1, i-1    */
  int yp1, ym1;               /* 3D sparse coordinate at j+1, j-1    */
  int xp1s, xm1s;             /* 2D sparse coordinate at i+1, i-1    */
  int yp1s, ym1s;             /* 2D sparse coordinate at j+1, j-1    */
  int xmym1s;                 /* Sparse coordinate at (i-1,j-1)      */
  int xmyp1s;                 /* Sparse coordinate at (i-1,j+1)      */
  int xmyp1;                  /* Sparse coordinate at (i-1,j+1)      */
  int sf = 0;                 /* Sponge flag                         */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  double c1;                  /* Constant for x diffusion            */
  double c2;                  /* Constant for y diffusion            */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(t11, 0, window->sgsiz * sizeof(double));
  memset(t12, 0, window->sgsiz * sizeof(double));
  memset(t22, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate the stress tensors.                                   */
  /* This is performed over all sparse cells (wet, auxiliary and     */
  /* ghost). The tensors t11 and t22 are cell centered, t12 is cell  */
  /* cornered thus the centered bottom coordinate must be used (this */
  /* is different to the face centered e1 bottom coordinate on       */
  /* western edges due to the stagger of the velocity cells).        */
  /* Auxiliary cells are not defined for cell centers on western and */
  /* southern boundaries, so the _t vectors cannot be used here and  */
  /* the cell center must be found by mapping downwards from the     */
  /* surface until self-mapping is located.                          */
  for (cc = 1; cc <= window->enonS; cc++) {
    c = cc;
    zm1 = window->zm1[c];
    while (c != zm1) {

      cs = window->m2d[c];
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      xm1 = window->xm1[c];
      ym1 = window->ym1[c];

      xp1s = window->xp1[cs];
      yp1s = window->yp1[cs];
      xm1s = window->xm1[cs];
      ym1s = window->ym1[cs];
      xmym1s = window->ym1[xm1s];

      /* Save the (x,x) direction stress in t11.                     */
      /* Cell centered.                                              */
      t11[c] = wincon->Ds[cs] *
        ((windat->u1b[xp1] - windat->u1b[c]) / window->h1acell[cs] +
         (windat->u2b[yp1] + windat->u2b[c]) *
         (window->h1au2[yp1s] - window->h1au2[cs]) / window->cellarea[cs]);

      /* Save the (y,y) direction stress in t22.                     */
      /* Cell centered.                                              */
      t22[c] = wincon->Ds[cs] *
        ((windat->u2b[yp1] - windat->u2b[c]) / window->h2acell[cs] +
         (windat->u1b[xp1] + windat->u1b[c]) *
         (window->h2au1[xp1s] - window->h2au1[cs]) / window->cellarea[cs]);

      /* Save the (x,y) direction stress in t12                      */
      /* Cell cornered.                                              */
      c1 = 0.25 * (wincon->Ds[cs] + wincon->Ds[xm1s] + wincon->Ds[ym1s] +
                   wincon->Ds[xmym1s]);

      t12[c] = c1 * (windat->u1b[c] / window->h1au1[cs] -
                     windat->u1b[ym1] / window->h1au1[ym1s]) *
        (window->h1au1[cs] + window->h1au1[ym1s]) /
        (window->h2au1[cs] + window->h2au1[ym1s]);

      t12[c] += c1 * (windat->u2b[c] / window->h2au2[cs] -
                      windat->u2b[xm1] / window->h2au2[xm1s]) *
        (window->h2au2[cs] + window->h2au2[xm1s]) /
        (window->h1au2[cs] + window->h1au2[xm1s]);

      c = zm1;
      zm1 = window->zm1[c];

    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the Smagorinsky diffusivity if required                     */
  if (wincon->smagorinsky != 0.0) {
    if (!windat->sdc)
      hd_quit("Smagorinsky diffusion requires tracer called smag\n");
    memset(windat->sdc, 0, window->sgsiz*sizeof(double));
    /* Calculate the Smagorinsky diffusion coefficient               */
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      zm1 = window->zm1[c];
      
      while (c != zm1) {

	xp1 = window->xp1[c];
	yp1 = window->yp1[c];
	cs = window->m2d[c];
	c1 = wincon->Ds[cs] * wincon->Ds[cs];

	/* Cell centered diffusivity                                   */
	c2 = 0.25 * (t12[c] + t12[xp1] + t12[yp1] + t12[window->yp1[xp1]]);
	windat->sdc[c] = sqrt(t11[c] * t11[c] / c1 +
			      0.5 * c2 * c2 / c1 +
			      t22[c] * t22[c] / c1) *
	  window->h1acell[cs] * window->h2acell[cs] * wincon->smagorinsky;
  	c = zm1;
	zm1 = window->zm1[c];
      }
    }

    /* Adjust the Smagorinsky coefficient for active alert actions   */
    /* and horizontal sponge zones.                                  */

    /* If smagcode & U1_A|U2_A|U1_AK|U2_AK then the relevant mixing  */
    /* array is allocated and filled with (grid weighted) constant   */
    /* coefficients in hd_init.c, and this code is not invoked.      */

    /* If smagcode & U1_SP|U2_SP|U1_SPK|U2_SPK then the relevant     */
    /* mixing array points to sdc, set in windows.c.                 */

    /* Set the flag for re-defining sponge zones                     */
    if (wincon->smagcode & (U1_SA|U2_SA)) sf = 1;
    if (wincon->smagcode & (U1_A|U2_A)) {
      if (wincon->alertf & (LEVEL1_U|LEVEL1_V|LEVEL1_UA|LEVEL1_VA))
	sf = 1;
    }

    /* Reset viscosity for U1_A/U2_A prior to shear alerts           */
    if (wincon->u1vh[0] == 1) {
      set_hdiff(window, wincon->u1vh, wincon->u1vh0, 0);
      sf = 1;
    }
    if (wincon->u2vh[0] == 1) {
      set_hdiff(window, wincon->u2vh, wincon->u2vh0, 1);
      sf = 1;
    }
    wincon->u1vh[0] = wincon->u2vh[0] = 0;

    /* Copy Smagorinsky diffusion to horizontal diffusion. This is   */
    /* not subject to LEVEL1 alerts or sponge zones.                 */
    if (wincon->smagcode & U1_SAK) {
      memcpy(wincon->u1kh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u1kh, wincon->kue1, 3);
    }
    if (wincon->smagcode & U2_SAK) {
      memcpy(wincon->u2kh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u2kh, wincon->kue2, 4);
    }

    /* Copy Smagorinsky diffusion to horizontal viscosity. This is   */
    /* subject to LEVEL1 alerts and separate sponge zones in the e1  */
    /* and e2 directions.                                            */
    if (wincon->smagcode & U1_SA) {
      memcpy(wincon->u1vh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u1vh, wincon->sue1, 1);
    }
    if (wincon->smagcode & U2_SA) {
      memcpy(wincon->u2vh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u2vh, wincon->sue2, 2);
    }

    /* Ensure the upper limit of the smagorinsky diffusion is the    */
    /* constant viscosity, and smooth.                               */
    if (wincon->smag_smooth && wincon->smagorinsky > 0.0) {
      if (wincon->u1vh0 < 0.0)
	smooth_w(window, wincon->u1vh, wincon->w6, window->w3_e1,
		 window->n3_e1, wincon->smag_smooth, -wincon->u1vh0);
      if (wincon->u2vh0 < 0.0)
	smooth_w(window, wincon->u2vh, wincon->w6, window->w3_e2,
		 window->n3_e2, wincon->smag_smooth, -wincon->u2vh0);
      if (wincon->u1kh0 < 0.0)
	smooth_w(window, wincon->u1kh, wincon->w6, window->w3_t,
		 window->n3_t, wincon->smag_smooth, -wincon->u1kh0);
      if (wincon->u2kh0 < 0.0)
	smooth_w(window, wincon->u2kh, wincon->w6, window->w3_t,
		 window->n3_t, wincon->smag_smooth, -wincon->u2kh0);
    }

    /* Reset diffusion throughout the water column for active alerts */
    /* Alerts at 3d locations.                                       */
    if (wincon->alertf & LEVEL1_U) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_U;
      /*dump_snapshot(sched_get_even_by_name(schedule, "dumps"), windat->t);*/
    }
    if (wincon->alertf & LEVEL1_V) {
      if (wincon->smagcode & U2_A) wincon->u2vh[0] = 1;
      for(cc = 1; (c = windat->cu2[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_V;
      /*dump_snapshot(sched_get_even_by_name(schedule, "dumps"), windat->t);*/
    }
    if (wincon->alertf & LEVEL1_W) {
      c = windat->cw;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL1_W;
    }
    if (wincon->alertf & LEVEL2w) {
      c = windat->cdw;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL2w;
    }
    /* Alerts at 2d locations.                                       */
    if (wincon->alertf & LEVEL1_UA) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1a[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_UA;
    }
    if (wincon->alertf & LEVEL1_VA) {
      if (wincon->smagcode & U2_A) wincon->u2vh[0] = 1;
      for(cc = 1; (c = windat->cu2a[cc]); cc++) {
	reset_hdiff(window, c, 0);
	reset_hdiff(window, c, 1);
      }
      wincon->alertf &= ~LEVEL1_VA;
    }
    if (wincon->alertf & LEVEL2) {
      c = windat->cdeta;
      reset_hdiff(window, c, 0);
      reset_hdiff(window, c, 1);
      wincon->alertf &= ~LEVEL2;
    }

    /* Increase diffusion if required */
    /*
    if (wincon->alertf & LEVEL1_U1VH) {
      for (cc = 1; cc <= window->n3_e1; cc++) {
	c = window->w3_e1[cc];
	wincon->u1vh[c] = min(-wincon->u1vh0, wincon->u1vh[c] / wincon->smagorinsky);
      }
      wincon->alertf &= ~LEVEL1_U1VH;
    }
    if (wincon->alertf & LEVEL1_U2VH) {
      for (cc = 1; cc <= window->n3_e2; cc++) {
	c = window->w3_e2[cc];
	wincon->u2vh[c] = min(-wincon->u2vh0, wincon->u2vh[c] / wincon->smagorinsky);
      }
      wincon->alertf &= ~LEVEL1_U2VH;
    }
    */

    /* Set the sponge zone                                           */
    /*set_sponge_w(window, windat, wincon);*/
    if (sf) set_sponge(window, wincon->u1vh, wincon->u2vh, windat->dt, 
		       wincon->w1);
  }

  /* Create a sponge zone on coarse-fine boundaries for zoomed grids */
  if (wincon->dozoom == DOZOOM) 
    set_sponge_zoom(window, windat, wincon, wincon->w2, wincon->w3,
		    window->nzbe1, window->nzbe2);

}

/* END hvisc_setup_old()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/*-------------------------------------------------------------------*/
void hvisc_u1_3d(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface sparse coordinates          */
  int zm1;                    /* 3D sparse coordinate at k-1         */
  int xm1;                    /* 3D sparse coordinate at i-1         */
  int yp1;                    /* 3D sparse coordinate at j+1         */
  int xm1s;                   /* 2D sparse coordinate at i-1         */
  int yp1s;                   /* 2D sparse coordinate at j+1         */
  int xmyp1s;                 /* Sparse coordinate at (i-1,j+1)      */
  int xmyp1;                  /* Sparse coordinate at (i-1,j+1)      */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double c1;                  /* Constant for x diffusion            */
  double c2;                  /* Constant for y diffusion            */
  double d1;                  /* Dummy                               */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  int *cells;                 /* Cells to process vector             */
  int vc;                     /* Size of cells[]                     */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;
  AH = wincon->u1vh;

  /*-----------------------------------------------------------------*/
  /* Use the thin layer cells to process if required                 */
  if (wincon->dolin_u1) {
    cells = wincon->s4;
    vc = wincon->acl;
  } else if (wincon->thin_merge) {
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion                           */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    
    xm1s = window->xm1[cs];
    yp1s = window->yp1[cs];
    xmyp1s = window->xmyp1[cs];
    
    d1 = window->h1au1[cs] * window->h2au1[cs];
    
    c1 = window->h2acell[cs] * AH[c] * t11[c] -
      window->h2acell[xm1s] * AH[xm1] * t11[xm1] +
      0.5 * (window->h1au2[yp1s] +
	     window->h1au2[xmyp1s]) * AH[yp1] * t12[yp1] -
      0.5 * (window->h1au2[cs] + window->h1au2[xm1s]) * AH[c] * t12[c];

    c2 = (0.5 * (window->h1au2[yp1s] + window->h1au2[xmyp1s]) -
	  0.5 * (window->h1au2[cs] + window->h1au2[xm1s])) *
      0.5 * (AH[c] * t12[c] + AH[yp1] * t12[yp1]) -
      (window->h2acell[cs] - window->h2acell[xm1s]) *
      0.5 * (AH[c] * t22[c] + AH[xm1] * t22[xm1]);
    
    /* Step forward in time                                          */
    windat->nu1[c] += windat->dt * (c1 + c2) / d1;
  }
}

/* END hvisc_u1_3d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hvisc_u1_3d_simple(geometry_t *window, /* Window geometry       */
                        window_t *windat,   /* Window data           */
                        win_priv_t *wincon  /* Window constants      */
  )
{
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface sparse coordinates          */
  int xp1, xm1;               /* 3D sparse coordinate at i+1, i-1    */
  int yp1, ym1;               /* 3D sparse coordinate at j+1, j-1    */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double c1;                  /* Constant for x diffusion            */
  double c2;                  /* Constant for y diffusion            */
  int *cells;                 /* Cells to process vector             */
  int vc;                     /* Size of cells[]                     */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  AH = wincon->u1vh;

  /*-----------------------------------------------------------------*/
  /* Use the thin layer cells to process if required                 */
  if (wincon->dolin_u1) {
    cells = wincon->s4;
    vc = wincon->acl;
  } else if (wincon->thin_merge) {
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }


  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion                           */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    c1 =
      AH[c] * (windat->u1b[xp1] + windat->u1b[xm1] -
               2.0 * windat->u1b[c]) / (window->h1au1[cs] *
                                        window->h1au1[cs]);

    c2 =
      AH[c] * (windat->u1b[yp1] + windat->u1b[ym1] -
               2.0 * windat->u1b[c]) / (window->h2au1[cs] *
                                      window->h2au1[cs]);

    windat->nu1[c] += windat->dt * (c1 + c2);
  }
}

/* END hvisc_u1_3d_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/*-------------------------------------------------------------------*/
void hvisc_u2_3d(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int cs;                       /* Surface / bottom sparse coordinates */
  int xp1;                      /* 3D sparse coordinate at i+1 */
  int ym1;                      /* 3D sparse coordinate at j-1 */
  int xp1s;                     /* 2D sparse coordinate at i+1 */
  int ym1s;                     /* 2D sparse coordinate at j-1 */
  int xpym1, xpym1s;            /* Sparse coordinate at (i-1,j+1) */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  double d1;                    /* Dummy */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  int *cells;                   /* Cells to process vector */
  int vc;                       /* Size of cells[] */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;
  AH = wincon->u2vh;

  if (wincon->dolin_u2) {
    cells = wincon->s4;
    vc = wincon->acl;
  } else if (wincon->thin_merge) {
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s2;
    vc = wincon->vc;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e2 horizontal diffusion */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    /* xpym1=window->ym1[xp1]; */
    xpym1 = window->xpym1[c];

    ym1s = window->xm1[cs];
    xp1s = window->yp1[cs];
    /* xpym1s=window->ym1[xp1s]; */
    xpym1s = window->xpym1[cs];

    d1 = window->h2au2[cs] * window->h1au2[cs];

    c1 = window->h1acell[cs] * AH[c] * t22[c] -
      window->h1acell[ym1s] * AH[ym1] * t22[ym1] +
      0.5 * (window->h2au1[xp1s] +
             window->h2au1[xpym1s]) * AH[xp1] * t12[xp1] -
      0.5 * (window->h2au1[cs] + window->h2au1[ym1s]) * AH[c] * t12[c];

    c2 = (0.5 * (window->h2au1[xp1s] + window->h2au1[xpym1s]) -
          0.5 * (window->h2au1[cs] + window->h2au1[ym1s])) *
      0.5 * (AH[c] * t12[c] + AH[xp1] * t12[xp1]) -
      (window->h1acell[cs] - window->h1acell[ym1s]) *
      0.5 * (AH[c] * t11[c] + AH[ym1] * t11[ym1]);

    /* Step forward in time */
    windat->nu2[c] += windat->dt * (c1 + c2) / d1;
  }
}

/* END hvisc_u2_3d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hvisc_u2_3d_simple(geometry_t *window, /* Window geometry */
                        window_t *windat, /* Window data           */
                        win_priv_t *wincon  /* Window constants 
                                             */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int cs;                       /* Surface sparse coordinates */
  int xp1, xm1;                 /* 3D sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* 3D sparse coordinate at j+1, j-1 */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  int *cells;                   /* Cells to process vector */
  int vc;                       /* Size of cells[] */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  AH = wincon->u2vh;

  if (wincon->dolin_u2) {
    cells = wincon->s4;
    vc = wincon->acl;
  } else if (wincon->thin_merge) {
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s2;
    vc = wincon->vc;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e2 horizontal diffusion */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    c1 =
      AH[c] * (windat->u2b[xp1] + windat->u2b[xm1] -
               2.0 * windat->u2b[c]) / (window->h1au2[cs] *
                                        window->h1au2[cs]);

    c2 =
      AH[c] * (windat->u2b[yp1] + windat->u2b[ym1] -
               2.0 * windat->u2b[c]) / (window->h2au2[cs] *
                                        window->h2au2[cs]);

    windat->nu2[c] += windat->dt * (c1 + c2);
  }
}

/* END hvisc_u2_3d_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Sets the mixing tensors for the 2D mode.                          */
/*-------------------------------------------------------------------*/
void hvisc_setup_2d(geometry_t *window,  /* Window geometry          */
		    window_t *windat,    /* Window data              */
		    win_priv_t *wincon,  /* Window constants         */
		    double *tzp          /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1 */
  int xmym1;                    /* Sparse coordinate at (i-1,j-1) */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;

  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(t11, 0, window->sgsizS * sizeof(double));
  memset(t12, 0, window->sgsizS * sizeof(double));
  memset(t22, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate the stress tensors.  */
  /* This is performed over all surface sparse cells                 */
  for (c = 1; c <= window->enonS; c++) {
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    xmym1 = window->ym1[xm1];

    /* Save the (x,x) direction stress in t11.  */
    /* Cell centered.  */

    t11[c] = wincon->Ds[c] * (tzp[c] - window->botz[c]) *
      (2.0*(windat->u1avb[xp1] - windat->u1avb[c]) / window->h1acell[c] +
       (windat->u2avb[yp1] + windat->u2avb[c]) *
       (window->h1au2[yp1] - window->h1au2[c]) / window->cellarea[c]);

    /* Save the (y,y) direction stress in t22.  */
    /* Cell centered.  */

    t22[c] = wincon->Ds[c] * (tzp[c] - window->botz[c]) *
      (2.0*(windat->u2avb[yp1] - windat->u2avb[c]) / window->h2acell[c] +
       (windat->u1avb[xp1] + windat->u1avb[c]) *
       (window->h2au1[xp1] - window->h2au1[c]) / window->cellarea[c]);

    /* Save the (x,y) direction stress in t12 */
    /* Cell cornered.  */
    c1 = 0.25 * (tzp[c] - window->botzgrid[c]) *
      (wincon->Ds[c] + wincon->Ds[xm1] + wincon->Ds[ym1] +
       wincon->Ds[xmym1]);

    t12[c] = c1 * (windat->u1avb[c] / window->h1au1[c] -
                   windat->u1avb[ym1] / window->h1au1[ym1]) *
      (window->h1au1[c] + window->h1au1[ym1]) /
      (window->h2au1[c] + window->h2au1[ym1]);

    t12[c] += c1 * (windat->u2avb[c] / window->h2au2[c] -
                    windat->u2avb[xm1] / window->h2au2[xm1]) *
      (window->h2au2[c] + window->h2au2[xm1]) /
      (window->h1au2[c] + window->h1au2[xm1]);
  }
}

/* END hvisc_setup_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Original pre v794 formulation : factor of 2 missing on stresses.  */
/* Sets the mixing tensors for the 2D mode.                          */
/*-------------------------------------------------------------------*/
void hvisc_setup_2d_old(geometry_t *window,  /* Window geometry      */
			window_t *windat,    /* Window data          */
			win_priv_t *wincon,  /* Window constants     */
			double *tzp      /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1 */
  int xmym1;                    /* Sparse coordinate at (i-1,j-1) */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;

  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(t11, 0, window->sgsizS * sizeof(double));
  memset(t12, 0, window->sgsizS * sizeof(double));
  memset(t22, 0, window->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate the stress tensors.  */
  /* This is performed over all surface sparse cells                 */
  for (c = 1; c <= window->enonS; c++) {
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    xmym1 = window->ym1[xm1];

    /* Save the (x,x) direction stress in t11.  */
    /* Cell centered.  */

    t11[c] = wincon->Ds[c] * (tzp[c] - window->botz[c]) *
      ((windat->u1avb[xp1] - windat->u1avb[c]) / window->h1acell[c] +
       (windat->u2avb[yp1] + windat->u2avb[c]) *
       (window->h1au2[yp1] - window->h1au2[c]) / window->cellarea[c]);

    /* Save the (y,y) direction stress in t22.  */
    /* Cell centered.  */
    t22[c] = wincon->Ds[c] * (tzp[c] - window->botz[c]) *
      ((windat->u2avb[yp1] - windat->u2avb[c]) / window->h2acell[c] +
       (windat->u1avb[xp1] + windat->u1avb[c]) *
       (window->h2au1[xp1] - window->h2au1[c]) / window->cellarea[c]);

    /* Save the (x,y) direction stress in t12 */
    /* Cell cornered.  */
    c1 = 0.25 * (tzp[c] - window->botzgrid[c]) *
      (wincon->Ds[c] + wincon->Ds[xm1] + wincon->Ds[ym1] +
       wincon->Ds[xmym1]);

    t12[c] = c1 * (windat->u1avb[c] / window->h1au1[c] -
                   windat->u1avb[ym1] / window->h1au1[ym1]) *
      (window->h1au1[c] + window->h1au1[ym1]) /
      (window->h2au1[c] + window->h2au1[ym1]);

    t12[c] += c1 * (windat->u2avb[c] / window->h2au2[c] -
                    windat->u2avb[xm1] / window->h2au2[xm1]) *
      (window->h2au2[c] + window->h2au2[xm1]) /
      (window->h1au2[c] + window->h1au2[xm1]);
  }
}

/* END hvisc_setup_2d_old()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the viscosity from the Smagorinsky mixing calculated in the  */
/* 3d mode. The 2d viscosity is the vertical integral of the 3d      */
/* viscosity and remains invariant over the 2d step. e1 viscosity is */
/* stored in wincon->w2 and e2 in viscosity wincon->w3.              */
/*-------------------------------------------------------------------*/
void set_viscosity_2d(geometry_t *window) /* Window geometry         */
{
  window_t *windat = window->windat;      /* Window data             */
  win_priv_t *wincon = window->wincon;    /* Window geometry         */
  int c, cc;                  /* Sparse coordinate / counter         */
  int cs;                     /* Surface / bottom sparse coordinates */
  int zm1;                    /* Sparse coordinate at k-1            */
  double depth;               /* Water depth                         */
  double dz;                  /* Cell thickness                      */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double *tzp = (wincon->nonlinear && !wincon->sigma) ? windat->eta : 
                                                        windat->topz;

  /* Set the e1 viscosity                                            */
  if (wincon->smagcode & (U1_SP|U1_SA)) {
    AH = wincon->w2;
    memset(AH, 0, window->sgsizS * sizeof(double));
    /* Vertically integrate the Smagorinsky diffusivity              */
    for (cc = 1; cc <= window->n2_e1; cc++) {
      c = window->w2_e1[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      while (c != zm1) {
	dz = 0.5 * (windat->dzu1[c] + windat->dzu1[window->xp1[c]]);
	AH[cs] += wincon->u1vh[c] * dz;
        c = zm1;
        zm1 = window->zm1[c];
      }
      depth = max(tzp[cs] - window->botz[cs], wincon->hmin);
      AH[cs] /= depth;
    }
  } else {
    memcpy(wincon->w2, wincon->u1vh, window->sgsizS * sizeof(double));
  }

  if (wincon->smagcode & (U2_SP|U2_SA)) {
    AH = wincon->w3;
    memset(AH, 0, window->sgsizS * sizeof(double));
    /* Vertically integrate the Smagorinsky diffusivity              */
    for (cc = 1; cc <= window->n2_e2; cc++) {
      c = window->w2_e2[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      while (c != zm1) {
	dz = 0.5 * (windat->dzu2[c] + windat->dzu2[window->yp1[c]]);
	AH[cs] += wincon->u2vh[c] * dz;
        c = zm1;
        zm1 = window->zm1[c];
      }
      depth = max(tzp[cs] - window->botz[cs], wincon->hmin);
      AH[cs] /= depth;
    }
  } else {
    memcpy(wincon->w3, wincon->u2vh, window->sgsizS * sizeof(double));
  }

  /* Create a sponge zone on coarse-fine boundaries for zoomed grids */
  if (wincon->dozoom == DOZOOM && wincon->mode2d) {
    set_sponge_zoom(window, windat, wincon, wincon->w2, wincon->w3,
		    window->nzbe1S, window->nzbe2S);
  }
}

/* END set_viscosity_2d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/*-------------------------------------------------------------------*/
void hvisc_u1_2d(geometry_t *window,  /* Processing window */
                 window_t *windat,  /* Window data structure */
                 win_priv_t *wincon,  /* Window geometry / constants */
                 double *tzp    /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, xm1;                 /* Sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* Sparse coordinate at j+1, j-1 */
  int xmym1;                    /* Sparse coordinate at (i-1,j-1) */
  int xmyp1;                    /* Sparse coordinate at (i-1,j+1) */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double depth;                 /* Water depth at u1 face */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  double d1;                    /* Dummy */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;
  AH = wincon->w2;

  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion */
  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];

    d1 = window->h1au1[c] * window->h2au1[c];
    depth = windat->depth_e1[c];

    c1 = window->h2acell[c] * AH[c] * t11[c] -
      window->h2acell[xm1] * AH[xm1] * t11[xm1] +
      0.5 * (window->h1au2[yp1] +
             window->h1au2[xmyp1]) * AH[yp1] * t12[yp1] -
      0.5 * (window->h1au2[c] + window->h1au2[xm1]) * AH[c] * t12[c];
    c2 = (0.5 * (window->h1au2[yp1] + window->h1au2[xmyp1]) -
          0.5 * (window->h1au2[c] + window->h1au2[xm1])) *
      0.5 * (AH[c] * t12[c] + AH[yp1] * t12[yp1]) -
      (window->h2acell[c] - window->h2acell[xm1]) *
      0.5 * (AH[c] * t22[c] + AH[xm1] * t22[xm1]);

    /* Step forward in time */
    windat->nu1av[c] += windat->dt2d * (c1 + c2) / (d1 * depth);
  }
  debug_c(window, D_UA, D_HDIFF);
}

/* END hvisc_u1_2d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hvisc_u1_2d_simple(geometry_t *window, /* Processing window */
                        window_t *windat, /* Window data structure */
                        win_priv_t *wincon, /* Window geometry / constants 
                                             */
                        double *tzp /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, xm1;                 /* 3D sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* 3D sparse coordinate at j+1, j-1 */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  AH = wincon->w2;

  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion */
  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    c1 = AH[c] / (window->h1au1[c] * window->h1au1[c]) *
      (windat->u1avb[xp1] + windat->u1avb[xm1] - 2 * windat->u1avb[c]);
    c2 = AH[c] / (window->h2au1[c] * window->h2au1[c]) *
      (windat->u1avb[yp1] + windat->u1avb[ym1] - 2 * windat->u1avb[c]);
    windat->nu1av[c] += windat->dt2d * (c1 + c2);
  }
  debug_c(window, D_UA, D_HDIFF);
}

/* END hvisc_u1_2d_simple()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/*-------------------------------------------------------------------*/
void hvisc_u2_2d(geometry_t *window,  /* Processing window */
                 window_t *windat,  /* Window data structure */
                 win_priv_t *wincon,  /* Window geometry / constants */
                 double *tzp    /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1;                      /* Sparse coordinate at i+1 */
  int ym1;                      /* Sparse coordinate at j-1 */
  int xpym1;                    /* Sparse coordinate at (i+1,j-1) */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */
  double *t11;                  /* Horizontal stress tensor, (x,x) */
  double *t12;                  /* Horizontal stress tensor, (x,y) */
  double *t22;                  /* Horizontal stress tensor, (y,y) */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double depth;                 /* Water depth at u1 face */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  double d1;                    /* Dummy */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;
  AH = wincon->w3;

  if(wincon->dolin_u2) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion */
  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];

    d1 = window->h2au2[c] * window->h1au2[c];
    depth = windat->depth_e2[c];

    c1 = window->h1acell[c] * AH[c] * t22[c] -
      window->h1acell[ym1] * AH[ym1] * t22[ym1] +
      0.5 * (window->h2au1[xp1] +
             window->h2au1[xpym1]) * AH[xp1] * t12[xp1] -
      0.5 * (window->h2au1[c] + window->h2au1[ym1]) * AH[c] * t12[c];

    c2 = (0.5 * (window->h2au1[xp1] + window->h2au1[xpym1]) -
          0.5 * (window->h2au1[c] + window->h2au1[ym1])) *
      0.5 * (AH[c] * t12[c] + AH[xp1] * t12[xp1]) -
      (window->h1acell[c] - window->h1acell[ym1]) *
      0.5 * (AH[c] * t11[c] + AH[ym1] * t11[ym1]);

    /* Step forward in time */
    windat->nu2av[c] += windat->dt2d * (c1 + c2) / (d1 * depth);
  }
  debug_c(window, D_VA, D_HDIFF);
}

/* END hvisc_u2_2d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hvisc_u2_2d_simple(geometry_t *window, /* Processing window */
                        window_t *windat, /* Window data structure */
                        win_priv_t *wincon, /* Window geometry / constants 
                                             */
                        double *tzp /* Pointer to surface level */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, xm1;                 /* 3D sparse coordinate at i+1, i-1 */
  int yp1, ym1;                 /* 3D sparse coordinate at j+1, j-1 */
  double *AH;                   /* e1 horizontal diffusion coefficient */
  double c1;                    /* Constant for x diffusion */
  double c2;                    /* Constant for y diffusion */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  AH = wincon->w3;

  if(wincon->dolin_u2) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion */
  for (cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    c1 = AH[c] / (window->h1au2[c] * window->h1au2[c]) *
      (windat->u2avb[xp1] + windat->u2avb[xm1] - 2 * windat->u2avb[c]);
    c2 = AH[c] / (window->h2au2[c] * window->h2au2[c]) *
      (windat->u2avb[yp1] + windat->u2avb[ym1] - 2 * windat->u2avb[c]);
    windat->nu2av[c] += windat->dt2d * (c1 + c2);
  }
  debug_c(window, D_VA, D_HDIFF);
}

/* END hvisc_u2_2d_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the scaled horizontal mixing                       */
/*-------------------------------------------------------------------*/
void set_hdiff(geometry_t *window,      /* Window geometry           */
	       double *AH,              /* Mixing coefficient        */
	       double AH0,              /* Constant mixing           */
	       int mode                 /* 0 = e1, 1 = e2            */
	       )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double *h1, hm;                       /* Cell size                 */
  int cc, c, c2;                        /* Counters, cell locations  */
  int nv;                               /* Stencil size              */
  int *vec;                             /* Cells to process          */

  /* Assign pointers                                                 */
  if (mode) {
    h1 = window->h2au1;
    hm = wincon->hmean2;
    nv = window->n3_e2;
    vec = window->w3_e2;
  } else {
    h1 = window->h1au1;
    hm = wincon->hmean1;
    nv = window->n3_e1;
    vec = window->w3_e1;
  }

  for (cc = 1; cc <= nv; cc++) {
    c = vec[cc];
    c2 = window->m2d[c];
    AH[c] = fabs(AH0 * h1[c2] / hm);
    /*AH[c] = fabs(AH0 * h1[c2] * h1[c2] / (hm * hm));*/
  }
}

/* END set_hdiff()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the scaled horizontal Smagorinsky mixing           */
/*-------------------------------------------------------------------*/
void scale_hdiff(geometry_t *window,      /* Window geometry         */
		 double *AH,              /* Mixing coefficient      */
		 double smag,             /* Smagorinsky coefficient */
		 int mode                 /* 0 = e1, 1 = e2          */
		 )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double *h1, hm;                       /* Cell size                 */
  int cc, c, c2;                        /* Counters, cell locations  */
  int nv;                               /* Stencil size              */
  int *vec;                             /* Cells to process          */
  double AH0;                           /* Constant mixing           */
  int scalef;                           /* Constant mixing scaling   */

  if (smag == 1.0 && wincon->smagorinsky != 0.0) return;

  /* Assign pointers                                                 */
  scalef = wincon->diff_scale;
  if (mode == 1) {
    h1 = window->h1au1;
    hm = wincon->hmean1;
    nv = window->n3_t;
    vec = window->w3_t;
    AH0 = wincon->bsue1;
  } else if (mode == 2) {
    h1 = window->h2au2;
    hm = wincon->hmean2;
    nv = window->n3_t;
    vec = window->w3_t;
    AH0 = wincon->bsue2;
  } else if (mode == 3) {
    h1 = window->h1acell;
    hm = wincon->hmean1;
    nv = window->n3_t;
    vec = window->w3_t;
    AH0 = wincon->bkue1;
  } else if (mode == 4) {
    h1 = window->h2acell;
    hm = wincon->hmean2;
    nv = window->n3_t;
    vec = window->w3_t;
    AH0 = wincon->bkue2;
  }

  for (cc = 1; cc <= nv; cc++) {
    c = vec[cc];
    c2 = window->m2d[c];
    AH[c] = smag * AH[c];
    if (scalef == LINEAR)
      AH[c] += fabs(AH0 * h1[c2] / hm);
    if (scalef == NONLIN)
      AH[c] += fabs(AH0 * h1[c2] * h1[c2] / (hm * hm));
  }
}


void scale_hdiff_m(master_t *master,   /* Master data                */
		   double *AH,         /* Mixing coefficient         */
		   double smag,        /* Smagorinsky coefficient    */
		   double AH0,         /* Base constant  mixing      */
		   int mode            /* 0 = e1, 1 = e2             */
		   )
{
  geometry_t *geom = master->geom;
  double *h1, hm;                       /* Cell size                 */
  int cc, c, c2;                        /* Counters, cell locations  */
  int nv;                               /* Stencil size              */
  int *vec;                             /* Cells to process          */
  int scalef;                           /* Constant mixing scaling   */

  /* Smagorinsky only                                                */
  if (smag == 1.0 && master->smagorinsky != 0.0) return;

  /* Assign pointers. Note; master->smagorinsky = 1.0                */
  scalef = master->diff_scale;

  if (mode == 1) {
    h1 = geom->h1au1;
    hm = master->hmean1;
    nv = geom->n3_t;
    vec = geom->w3_t;
  } else if (mode == 2) {
    h1 = geom->h2au2;
    hm = master->hmean2;
    nv = geom->n3_t;
    vec = geom->w3_t;
  } else if (mode == 3) {
    h1 = geom->h1acell;
    hm = master->hmean1;
    nv = geom->n3_t;
    vec = geom->w3_t;
  } else if (mode == 4) {
    h1 = geom->h2acell;
    hm = master->hmean2;
    nv = geom->n3_t;
    vec = geom->w3_t;
  }

  for (cc = 1; cc <= nv; cc++) {
    c = vec[cc];
    c2 = geom->m2d[c];

    AH[c] = smag * AH[c];
    if (scalef == LINEAR)
      AH[c] += fabs(AH0 * h1[c2] / hm);
    if (scalef == NONLIN)
      AH[c] += fabs(AH0 * h1[c2] * h1[c2] / (hm * hm));
  }
}

/* END scale_hdiff()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to increase horizontal diffusion over a predefined area.  */
/* The stencil area is assumed to be star shaped.                    */
/*-------------------------------------------------------------------*/
void reset_hdiff(geometry_t *window,    /* Window geometry           */
		 int cl,                /* Central coordinate        */
		 int mode               /* 0 = e1, 1 = e2            */
		 )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double *VH, *VH0, vh;                 /* Modified diffusion        */
  double *h1, *h2;                      /* Cell size                 */
  double rm, rmax;                      /* Maximum diffusion         */
  double sfact = 0.9;                   /* Safety factor             */
  double mfact = 4.0;                   /* Multiplication factor     */
  int cc, c, cs, n;                     /* Counters, cell locations  */
  int ns;                               /* Stencil size              */
  int *dist, sz;                        /* Distance array            */
  int *st, size = 15;                   /* Stencil diameter          */
  int smag = 0;

  if (!wincon->hdiff_f) return;

  /* Assign pointers                                                 */
  cs = window->m2d[cl];
  if (mode) {
    h1 = window->h2au1;
    h2 = window->h2au2;
    VH = wincon->u2vh;
    VH0 = windat->u2vhin;
    rm = fabs(wincon->u2vh0 * window->h2au2[cs] / wincon->hmean2);
    if (wincon->u2vh0 < 0.0) smag = 1;
  } else {
    h1 = window->h1au1;
    h2 = window->h1au2;
    VH = wincon->u1vh;
    rm = fabs(wincon->u1vh0 * window->h1au1[cs] / wincon->hmean1);
    VH0 = windat->u1vhin;
    if (wincon->u1vh0 < 0.0) smag = 1;
  }
  rmax = 1.0 / (h1[cs] * h1[cs]) + 1.0 / (h2[cs] * h2[cs]);
  rmax = sfact / (4.0 * rmax * window->windat->dt);
  rm = min(mfact * rm, rmax);

  /* Get the stencil                                                 */
  ns = size;
  if (wincon->visc_method & ROAM_WIN) {
    size = ns = 5;
    st = stencil(window, cl, &ns, 0);
    dist = i_alloc_1d(ns + 1);
    for (cc = 0; cc < ns; cc++) dist[cc] = 2;
    dist[12] = 0;
    dist[6] = dist[7] = dist[8] = 1;
    dist[11] = dist[13] = 1;
    dist[16] = dist[17] = dist[18] = 1;
  } else {
    st = stencil(window, cl, &ns, 1);
    /* Get the distance of stencil locations from the center           */
    sz = floor(size / 2);
    dist = i_alloc_1d(ns + 1);
    for (cc = 0; cc < 4; cc++)
      dist[cc] = sz;
    for (n = sz - 1; n > 0; n--) {
      for (c = 1; c <= 4 * (n + 1); c++) {
	dist[cc++] = n;
      }
    }
  }
  for (cc = 0; cc < ns; cc++) {
    c = window->m2d[st[cc]];
    while (c != window->zm1[c]) {
      cs = window->m2d[c];
      vh = (smag) ? VH[c] : VH0[cs];
      /* Get the maximum stable diffusion                              */
      VH[c] = rm + (double)dist[cc] * (vh - rm) / (double)(sz + 1);
      /*VH[c] = h1[cs] * h1[cs] / (4.0 * windat->dt);*/
      c = window->zm1[c];
    }
  }
  i_free_1d(dist);
  i_free_1d(st);
  windat->nalert[7]++;
  /*
  if (mode == 0)
    hd_warn("U1VH reset at %f days : (%d %d)\n", windat->t/86400,
	    window->s2i[cl], window->s2j[cl]);
  if (mode == 0)
    hd_warn("U2VH reset at %f days : (%d %d)\n", windat->t/86400,
	    window->s2i[cl], window->s2j[cl]);
  */
}

/* END reset_hdiff()                                                 */
/*-------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/* Optimised horizontal mixing                                      */
/*------------------------------------------------------------------*/
void reset_hor_diff(master_t *master, double u1vh, double u2vh, int flag)
{
  geometry_t *geom = master->geom;
  double vh1, vh2;
  int c, cc, c2;

  master->u1vh0 = u1vh;
  master->u2vh0 = u2vh;
  if(flag == AUTO) {
    double hf = 0.05;             /* Factor for horizontal diffusion */
    double step = 1;              /* Integral step of diffusion > 1  */
    double hmax = 1e10;
    double d1, d2;
    int u1khf = 0, u1vhf = 0, u2khf = 0, u2vhf = 0, i1, cs;
    if (u1vh <= 0.0) u1vhf = 1;      
    if (u2vh <= 0.0) u2vhf = 1;
    for (cc = 1; cc <= geom->n3_t; cc++) {
      c = geom->w3_t[cc];
      cs = geom->m2d[c];
      hmax = 1.0 / ((1.0 / (geom->h1acell[cs] * geom->h1acell[cs]) +
		     1.0 / (geom->h2acell[cs] * geom->h2acell[cs])) * 4.0 *
		    master->grid_dt);
      i1 = (int)hmax / (int)step;
      hmax = step * (double)i1;
      d1 = 0.01 * geom->h1acell[cs] * geom->h1acell[cs] / master->grid_dt;
      d2 = 0.01 * geom->h2acell[cs] * geom->h2acell[cs] / master->grid_dt;
      i1 = (int)d1 / (int)step;
      if (u1vhf)
	master->u1vh[c] = step * (double)i1;
      i1 = (int)d2 / (int)step;
      if (u2vhf)
	master->u2vh[c] = step * (double)i1;
      /* Set limits */
      if (u1vhf) {
	if (master->u1vh[c] > hmax)
	  master->u1vh[c] = hmax;
	if (master->u1vh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
        master->u1vh[c] = step * (double)i1;
	}
      }
      if (u2vhf) {
	if (master->u2vh[c] > hmax)
	  master->u2vh[c] = hmax;
	if (master->u2vh[c] < hf * hmax) {
	  i1 = (int)(hf * hmax) / (int)step;
	  master->u2vh[c] = step * (double)i1;
	}
      }
    }
    smooth3(master, master->u1vh, geom->w3_e1, geom->n3_e1);
    smooth3(master, master->u2vh, geom->w3_e2, geom->n3_e2);
  } else {
    memset(master->u1vh, 0.0, geom->sgsiz * sizeof(double));
    if (u1vh > 0.0) {
      vh1 = (master->bsue1) ? master->bsue1 : u1vh;
      if (master->smagcode & U1_SA)
	memcpy(master->u1vh, master->sdc, geom->sgsiz * sizeof(double));
    } else
      vh1 = 0.01 * master->hmean1 * master->hmean1 / master->grid_dt;
    scale_hdiff_m(master, master->u1vh, master->sue1, vh1, 1);

    memset(master->u2vh, 0.0, geom->sgsiz * sizeof(double));
    if (u2vh > 0.0) {
      vh2 = (master->bsue2) ? master->bsue2 : u2vh;
      if (master->smagcode & U2_SA)
	memcpy(master->u2vh, master->sdc, geom->sgsiz * sizeof(double));
    } else
      vh2 = 0.01 * master->hmean2 * master->hmean2 / master->grid_dt;
    scale_hdiff_m(master, master->u2vh, master->sue2, vh2, 2);
  }
}

/* END reset_hor_diff()                                             */
/*------------------------------------------------------------------*/
