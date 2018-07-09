/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/momentum/mommix.c
 *  
 *  Description: Horizontal mixing of momentum
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: mommix.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void scale_hdiff(geometry_t *window, double *AH, double smag, int mode);

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

  /* Set the horizontal mixing method                                */
  for (n = 1; n <= master->nwindows; n++) {
    wincon[n]->hor_mix = (mix_method_t *)malloc(sizeof(mix_method_t));
    wincon[n]->hor_mix->pre = hvisc_setup_pre;
    wincon[n]->hor_mix->setup = hvisc_setup;
    wincon[n]->hor_mix->setup_av = hvisc_setup_2d;
    if (wincon[n]->visc_method & LAPLACIAN) {
      wincon[n]->hor_mix->u1 = hvisc_u1_3d;
      wincon[n]->hor_mix->u1av = hvisc_u1_2d;
    } else if (wincon[n]->visc_method & US_LAPLACIAN) {
      wincon[n]->hor_mix->u1 = hvisc_u1_3dus;
      wincon[n]->hor_mix->u1av = hvisc_u1_2dus;
    } else if (wincon[n]->visc_method & SIMPLE) {
      if (wincon[n]->u1vh0 >= 0.0 && wincon[n]->u2vh0 >= 0.0) {
	wincon[n]->hor_mix->setup = hvisc_null;
	wincon[n]->hor_mix->setup_av = hvisc_2d_null;
      }
      wincon[n]->hor_mix->u1 = hvisc_u1_3d_simple;
      wincon[n]->hor_mix->u1av = hvisc_u1_2d_simple;
    }
    if (wincon[n]->u1_f & HDIFF) {
      wincon[n]->hor_mix->u1 = hvisc_null;
    }
    if (wincon[n]->u1av_f & HDIFF) {
      wincon[n]->hor_mix->u1av = hvisc_2d_null;
    }
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
/* wet+auxilliary+ghost array.                                       */
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
  int c, cv, cs, cvs, cc;     /* Cell coordinate / counter           */
  int e, ee, es;              /* Edge coordinates                    */
  int vv, v, vs;              /* Vertex coordinate : c2v[1]          */
  int n;                      /* Counter                             */
  int zm1;                    /* 3D cell coordinate at k-1           */
  int e1, e2, e3, e4;         /* Edges surrounding cell c            */
  int e1s, e2s, e3s, e4s;     /* Surface edges surrounding c         */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  double c1;                  /* Constant for x diffusion            */
  double c2;                  /* Constant for y diffusion            */
  double h1, h2, h3, h4;      /* Vertex mean of h2au1 and h1au1      */

  if (wincon->smagorinsky == 0.0) return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  t11 = wincon->t11;
  t12 = wincon->t12;
  t22 = wincon->t22;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  memset(t11, 0, window->sze * sizeof(double));
  memset(t12, 0, window->sze * sizeof(double));
  memset(t22, 0, window->sze * sizeof(double));

  vel_grad(window, windat, wincon, windat->u1b, t11, 
	   wincon->w5, GRAD_3D|GRAD_NOR);
  vel_grad(window, windat, wincon, windat->u1b, wincon->w6, 
	   t22, GRAD_3D|GRAD_TAN);

  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    e1 = window->m2d[window->e2c[e][0]];
    e2 = window->m2d[window->e2c[e][1]];
    c1 = 0.5 * (wincon->Ds[e1] + wincon->Ds[e2]);

    t11[e] *= c1;
    t22[e] *= c1;
    t12[e] = c1 * (wincon->w5[e] + wincon->w6[e]);

  }

  /*-----------------------------------------------------------------*/
  /* Set the Smagorinsky diffusivity if required                     */
  if (wincon->smagorinsky != 0.0) {
    if (!windat->sdc)
      hd_quit("Smagorinsky diffusion requires tracer called smag\n");
    memset(windat->sdc, 0, window->sgsiz*sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      for (n = 1; n <= window->npe[cs]; n++) {
	e = window->c2e[n][c];
	es = window->m2de[e];
	windat->sdc[c] += (window->h1au1[es] * window->h1au1[es] *
			   sqrt(t11[e] * t11[e] +
				0.5 * t12[e] * t12[e] +
				t22[e] * t22[e]));
      }
      windat->sdc[c] *= (wincon->smagorinsky / (double)window->npe[cs]);
    }
  }

  /* Average onto the cell centres and vertices                      */
  memcpy(wincon->w5, t11, window->sze * sizeof(double));
  memcpy(wincon->w6, t22, window->sze * sizeof(double));
  memset(t11, 0, window->szc * sizeof(double));
  memset(t22, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    c1 = 1.0 / (double)window->npe[cs];
    for (ee = 1; ee <= window->npe[cs]; ee++) {
      e = window->c2e[ee][c];
      t11[c] += wincon->w5[e] * c1;
      t22[c] += wincon->w6[e] * c1;
    }
  }
  memcpy(wincon->w5, t12, window->sze * sizeof(double));
  memset(t12, 0, window->szv * sizeof(double));
  for (vv = 1; vv <= window->b3_e2; vv++) {
    v = window->w3_e2[vv];
    vs = window->m2dv[v];
    c1 = 1.0 / window->nve[vs];
    for (ee = 1; ee <= window->nve[vs]; ee++) {
      e = window->v2e[v][ee];
      if (e)
	t12[v] += wincon->w5[e] *c1;
    }
  }

  if (wincon->smagcode & (U1_SP|U1_SA)) {
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      wincon->u1vh[e] = 0.5 * (windat->sdc[window->e2c[e][0]] + 
			       windat->sdc[window->e2c[e][1]]);
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
  int c, cc, ee, e;           /* Sparse coordinate / counter         */
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
      set_hdiff(window, wincon->u1vh, wincon->u1vh0);
      sf = 1;
    }
    wincon->u1vh[0] = 0;

    /* Copy Smagorinsky diffusion to horizontal diffusion. This is   */
    /* not subject to LEVEL1 alerts or sponge zones.                 */
    if (wincon->smagcode & U1_SAK) {
      memcpy(wincon->u1kh, windat->sdc, window->sgsiz * sizeof(double));
      scale_hdiff(window, wincon->u1kh, wincon->kue1, 3);
    }

    /* Copy Smagorinsky diffusion to horizontal viscosity. This is   */
    /* subject to LEVEL1 alerts and separate sponge zones in the e1  */
    /* and e2 directions.                                            */
    /* Uncomment this code if separate sponges apply in the e1 and   */
    /* e2 directions.                                                */
    if (wincon->smagcode & U1_SA) {
      for (ee = 1; ee <= window->n3_e1; ee++) {
	e = window->w3_e1[ee];
	wincon->u1vh[e] = 0.5 * (windat->sdc[window->e2c[e][0]] +
				 windat->sdc[window->e2c[e][1]]);
      }
      scale_hdiff(window, wincon->u1vh, wincon->sue1, 1);
    }
    /* Reset diffusion throughout the water column for active alerts */
    /* Alerts at 3d locations.                                       */
    if (wincon->alertf & LEVEL1_U) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1[cc]); cc++) {
	reset_hdiff(window, c);
      }
      wincon->alertf &= ~LEVEL1_U;
    }
    if (wincon->alertf & LEVEL1_W) {
      c = windat->cw;
      reset_hdiff(window, c);
      wincon->alertf &= ~LEVEL1_W;
    }
    if (wincon->alertf & LEVEL2w) {
      c = windat->cdw;
      reset_hdiff(window, c);
      wincon->alertf &= ~LEVEL2w;
    }
    /* Alerts at 2d locations.                                       */
    if (wincon->alertf & LEVEL1_UA) {
      if (wincon->smagcode & U1_A) wincon->u1vh[0] = 1;
      for(cc = 1; (c = windat->cu1a[cc]); cc++) {
	reset_hdiff(window, c);
      }
      wincon->alertf &= ~LEVEL1_UA;
    }
    if (wincon->alertf & LEVEL2) {
      c = windat->cdeta;
      reset_hdiff(window, c);
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
    if (sf) {
      reset_sponge_zone(window);
      set_sponge(window, wincon->u1vh, windat->dt, wincon->w1);
    }

    /* Ensure the upper limit of the smagorinsky diffusion is the    */
    /* constant viscosity, and smooth.                               */
    if (wincon->smag_smooth && wincon->smagorinsky > 0.0) {
      if (wincon->u1vh0 < 0.0)
	smooth_w(window, wincon->u1vh, wincon->w6, window->w3_e1,
		 window->n3_e1, wincon->smag_smooth, -wincon->u1vh0);
      if (wincon->u1kh0 < 0.0)
	smooth_w(window, wincon->u1kh, wincon->w6, window->w3_t,
		 window->n3_t, wincon->smag_smooth, -wincon->u1kh0);
    }
  }
}

/* END hvisc_setup()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion for unstructured grids and variable eddy     */
/* viscosity using Laplacian mixing. Eq. A.8 and A.39 in             */
/* Ringler et al, 2013.                                              */
/* Note: divergence and relative vorticity are computed in           */
/* nonlin_coriolis_3d().                                             */
/* Note: (rvor[v2] - rvor[v1]) / h2au1[e] is rvor pointing from      */
/* vertex 2 to vertex 1, or equivalently + k times rvor pointing     */
/* from c1 to c2.                                                    */
/* Note also the tangent velocity points to v1 here, whereas it      */
/* poinrs to v2 in Ringler et al, 2010, Fig3.                        */
/*-------------------------------------------------------------------*/
void hvisc_u1_3dus(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
		   )
{
  int e, ee;                  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c1, c2;                 /* Cells adjacent to e                 */
  int v1, v2;                 /* Vertices adjacent to e              */
  double d1;                  /* Dummy                               */
  double dtu;
  int cc, c, cs;
  int vv, v, vs;
  double d2, iarea;

  dtu = windat->dtf + windat->dtb;

  /*-----------------------------------------------------------------*/
  /* Divergence is given by Eq. 21, Ringler et al, 2010.             */
  /*if (!(wincon->momsc & RINGLER)) {*/
    memset(windat->div, 0, window->sgsiz * sizeof(double));
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      d1 = 1.0 / window->cellarea[cs];
      for (ee = 1; ee <= window->npe[cs]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	d2 = window->h1au1[es] * windat->u1b[e] * d1;
	windat->div[c] += window->eSc[ee][cs] * d2;
      }
    }
    memset(windat->rvor, 0, window->szv * sizeof(double));
    for (vv = 1; vv <= window->a3_e2; vv++) {
      v = window->w3_e2[vv];
      vs = window->m2dv[v];
      iarea = 1.0 / window->dualarea[vs];
      for (ee = 1; ee <= window->nve[vs]; ee++) {
	e = window->v2e[v][ee];
	if (e) {
	  es = window->m2de[e];
	  d2 = window->h2au1[es] * windat->u1b[e];
	  windat->rvor[v] += window->eSv[ee][v] * iarea * d2;
	}
      }
    }
    /*-----------------------------------------------------------------*/
    /* Set relative vorticity along solid boudaries equal to zero.     */
    for (vv = window->a3_e2 + 1; vv <= window->n3_e2; vv++) {
      v = window->w3_e2[vv];
      windat->rvor[v] = 0.0;
    }


  /*-----------------------------------------------------------------*/
  /* Compute the horizontal mixing tendency                          */
  for (ee = 1; ee <= window->v3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    v1 = window->e2v[e][0];
    v2 = window->e2v[e][1];
    d1 = (windat->div[c1] - windat->div[c2]) / window->h2au1[es] +
      (windat->rvor[v2] - windat->rvor[v1]) / window->h1au1[es] ;
    windat->nu1[e] += windat->dt * wincon->u1vh[e] * d1;
  }
}

/* END hvisc_u1_3dus()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/* Old code - no longer used for unstructured grids since we use     */
/* hvisc_u1_3dus() instead. We retain this code for a possible       */
/* future tensor approach.                                           */
/*-------------------------------------------------------------------*/
void hvisc_u1_3d(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon   /* Window constants            */
  )
{
  int e, ee;                  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c, c1, c2;              /* Cells adjacent to e                 */
  int c1s, c2s, cs;           /* Surface cells                       */
  int v;                      /* Vertex coordinate                   */
  int vc;                     /* Size of cells[]                     */
  int n, m;                   /* Counters                            */
  int j1, j2;                 /* Directions from edge                */
  int zm1;                    /* 3D sparse coordinate at k-1         */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double d1;                  /* Constant for x diffusion            */
  double d2;                  /* Constant for y diffusion            */
  double d3;                  /* Dummy                               */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  int *cells;                 /* Cells to process vector             */
  double nvc;                 /* Cells surrounding vertices          */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  t12 = wincon->t12;
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
  for (ee = 1; ee <= vc; ee++) {
    e = cells[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    c1s = window->m2d[c1];
    c2s = window->m2d[c2];

    j1 = window->e2e[e][0];
    j2 = jp(j1, window->npe[c1s]);
    t11 = (j1 == 1) ? wincon->t11 : wincon->t22;
    t22 = (j1 == 1) ? wincon->t22 : wincon->t11;

    d3 = window->h1au1[es] * window->h2au1[es];

    d2 = 0.0;
    d1 = AH[e] * (2.0 * window->hacell[j2][c1s] * t11[c1] -
		  2.0 * window->hacell[j2][c2s] * t11[c2]);

    for (m = 0; m <= 1; m++) {
      int vs;
      double sgn = (m) ? -1.0 : 1.0;
      v = window->e2v[e][m];
      vs = window->m2dv[v];
      nvc = window->nvc[vs];
      for (n = 1; n <= nvc; n++) {
	c = window->v2c[v][n];
	cs = window->m2d[c];
	d1 += (sgn * t12[v] * AH[e] * 2.0 * window->hacell[j1][cs]) / (double)nvc;
	d2 += (sgn * AH[e] * 2.0 * window->hacell[j1][cs]) / (double)nvc;
      }
      d2 *= 0.5 * t12[v];
    }

    d2 -=  (2.0 * (window->hacell[j2][c1s] - window->hacell[j2][c2s])) * 
      0.5 * AH[e] * (t22[c1] + t22[c2]);

    /* Step forward in time                                          */
    windat->nu1[e] += windat->dt * (d1 + d2) / d3;
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
  int e, ep, em, eu, ed, ee;  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c1, c2;                 /* Cells adjacent to e                 */
  int j1, j2;
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double d1;                  /* Constant for x diffusion            */
  double d2;                  /* Constant for y diffusion            */
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
  for (ee = 1; ee <= vc; ee++) {
    e = cells[ee];
    ep = window->ep[e];
    em = window->em[e];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];

    j1 = window->e2e[e][0];

    eu = (j1 == 1) ? e2e(window, e, 2) : e2e(window, e, 3);
    ed = (j1 == 1) ? e2e(window, e, 4) : e2e(window, e, 1);

    c1 = AH[e] * (windat->u1b[ep] + windat->u1b[em] -      
		  2.0 * windat->u1b[e]) / (window->h2au1[es] *
					   window->h2au1[es]);

    c2 = AH[e] * (windat->u1b[eu] + windat->u1b[ed] -
		  2.0 * windat->u1b[e]) / (window->h1au1[es] *
					   window->h1au1[es]);

    windat->nu1[e] += windat->dt * (c1 + c2);
  }
}

/* END hvisc_u1_3d_simple()                                          */
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
  int c, cs, cv, cc;            /* Cell coordinate / counter         */
  int e, es;                    /* Edge coordinates                  */
  int vv, v;                    /* Vertex coordinate : c2v[1]        */
  int n;                        /* Counter                           */
  int e1, e2, e3, e4;           /* Edges surrounding cell c          */
  int *ctp;                     /* Cells to process vector           */
  int vcs;                      /* Number of cells to process        */
  double *t11;                  /* Horizontal stress tensor, (x,x)   */
  double *t12;                  /* Horizontal stress tensor, (x,y)   */
  double *t22;                  /* Horizontal stress tensor, (y,y)   */
  double c1;                    /* Constant for x diffusion          */
  double c2;                    /* Constant for y diffusion          */
  double h1, h2, h3, h4;        /* Vertex mean of h2au1 and h1au1    */

  /*-----------------------------------------------------------------*/
  /* These tensors are no longer required if momentum is solved      */
  /* for unstructured grids using hvisc_u1_2dus().                   */
  /* The routine is retained for a possible future tensor approach.  */
}

/* END hvisc_setup_2d()                                              */
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
  int e, ee;                  /* Sparse coordinate / counter         */
  int es;                     /* Surface / bottom sparse coordinates */
  int zm1;                    /* Sparse coordinate at k-1            */
  double depth;               /* Water depth                         */
  double dz;                  /* Cell thickness                      */
  double *AH;                 /* e1 horizontal diffusion coefficient */

  /* Set the e1 viscosity                                            */
  if (wincon->smagcode & (U1_SP|U1_SA)) {
    AH = wincon->w2;
    memset(AH, 0, window->szeS * sizeof(double));
    depth = 0.0;
    /* Vertically integrate the Smagorinsky diffusivity              */
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      es = window->m2de[e];
      zm1 = window->zm1e[e];
      while (e != zm1) {
	AH[es] += wincon->u1vh[e] * windat->dzu1[e];
	depth += windat->dzu1[e];
        e = zm1;
        zm1 = window->zm1e[e];
      }
      AH[es] /= depth;
    }
  } else {
    memcpy(wincon->w2, wincon->u1vh, window->szeS * sizeof(double));
  }
}

/* END set_viscosity_2d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion for unstructured grids and variable eddy     */
/* viscosity using Laplacian mixing. Eq. A.39 Ringler et al, 2013.   */
/* Note: divergence and relative vorticity are computed in           */
/* nonlin_coriolis_2d().                                             */
/*-------------------------------------------------------------------*/
void hvisc_u1_2dus(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon,  /* Window constants          */
		   double *tzp          /* Pointer to surface level  */
		   )
{
  int e, ee;                  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c1, c2;                 /* Cells adjacent to e                 */
  int v1, v2;                 /* Vertices adjacent to e              */
  double d1;                  /* Dummy                               */
  double dtu;
  int cc, c, cs;
  int vv, v, vs;
  double d2, iarea;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  dtu = windat->dtu1 / d1;          /* 2D forward time-step          */

  /*-----------------------------------------------------------------*/
  /* Divergence is given by Eq. 21, Ringler et al, 2010.             */
  /*if (!(wincon->momsc & RINGLER)) {*/
    memset(windat->div, 0, window->szcS * sizeof(double));
    for(cc = 1; cc <= window->a2_t; cc++) {
      c = window->w2_t[cc];
      cs = window->m2d[c];
      d1 = 1.0 / window->cellarea[cs];
      for (ee = 1; ee <= window->npe[cs]; ee++) {
	e = window->c2e[ee][c];
	es = window->m2de[e];
	d2 = window->h1au1[es] * windat->u1avb[e] * d1;
	windat->div[c] += window->eSc[ee][cs] * d2;
      }
    }
    memset(windat->rvor, 0, window->szvS * sizeof(double));
    for (vv = 1; vv <= window->a2_e2; vv++) {
      v = window->w2_e2[vv];
      vs = window->m2dv[v];
      iarea = 1.0 / window->dualarea[vs];
      for (ee = 1; ee <= window->nve[vs]; ee++) {
	e = window->v2e[v][ee];
	if (e) {
	  d2 = window->h2au1[e] * windat->u1avb[e];
	  windat->rvor[v] += window->eSv[ee][v] * iarea * d2;
	}
      }
    }
    /*-----------------------------------------------------------------*/
    /* Set relative vorticity along solid boudaries equal to zero.     */
    for (vv = window->a2_e2 + 1; vv <= window->n2_e2; vv++) {
      v = window->w2_e2[vv];
      windat->rvor[v] = 0.0;
    }


  /*-----------------------------------------------------------------*/
  /* Compute the horizontal mixing tendency                          */
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    v1 = window->e2v[e][0];
    v2 = window->e2v[e][1];
    d1 = (windat->div[c1] - windat->div[c2]) / window->h2au1[es] +
      (windat->rvor[v2] - windat->rvor[v1]) / window->h1au1[es] ;
    windat->nu1av[e] += windat->dt2d * wincon->u1vh[e] * d1;
  }
}

/* END hvisc_u1_2dus()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion for non-uniform grids and variable eddy      */
/* viscosity. Uses the full terms in the horizontal diffusion        */
/* equation.                                                         */
/* Old code - no longer used for unstructured grids since we use     */
/* hvisc_u1_2dus() instead. We retain this code for a possible       */
/* future tensor approach.                                           */
/*-------------------------------------------------------------------*/
void hvisc_u1_2d(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data                 */
                 win_priv_t *wincon,  /* Window constants            */
                 double *tzp          /* Pointer to surface level    */
  )
{
  int e, ee;                  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c, c1, c2;              /* Cells adjacent to e                 */
  int v;                      /* Vertex coordinate                   */
  int vc;                     /* Size of cells[]                     */
  int n, m;                   /* Counters                            */
  int j1, j2;                 /* Directions from edge                */
  int *ctp;                   /* Cells to process vector             */
  int vcs;                    /* Number of cells to process          */
  double *t11;                /* Horizontal stress tensor, (x,x)     */
  double *t12;                /* Horizontal stress tensor, (x,y)     */
  double *t22;                /* Horizontal stress tensor, (y,y)     */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double depth;               /* Water depth at u1 face              */
  double d1;                  /* Constant for x diffusion            */
  double d2;                  /* Constant for y diffusion            */
  double d3;                  /* Dummy                               */
  double nvc;                 /* Cells surrounding vertices          */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
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

  vcs = window->v2_t;
  ctp = window->w2_t;

  /*-----------------------------------------------------------------*/
  /* Calculate the e1 horizontal diffusion                           */
  for (ee = 1; ee <= vcs; ee++) {
    e = ctp[ee];
    es = window->m2de[e];
    c1 = window->m2d[window->e2c[e][0]];
    c2 = window->m2d[window->e2c[e][1]];

    j1 = window->e2e[e][0];
    j2 = jp(j1, window->npe[c1]);
    t11 = (j1 == 1) ? wincon->t11 : wincon->t22;
    t22 = (j1 == 1) ? wincon->t22 : wincon->t11;

    d3 = window->h1au1[es] * window->h2au1[es];
    depth = windat->depth_e1[es];
    d2 = 0.0;
    d1 = AH[es] * (2.0 * window->hacell[j2][c1] * t11[c1] -
		   2.0 * window->hacell[j2][c2] * t11[c2]);

    for (m = 0; m <= 1; m++) {
      double sgn = (m) ? -1.0 : 1.0;
      v = window->m2dv[window->e2v[e][m]];
      nvc = window->nvc[v];
      for (n = 1; n <= nvc; n++) {
	c = window->v2c[v][n];
	d1 += (sgn * t12[v] * AH[es] * 2.0 * window->hacell[j1][c]) / (double)nvc;
	d2 += (sgn * AH[es] * 2.0 * window->hacell[j1][c]) / (double)nvc;
      }
      d2 *= 0.5 * t12[v];
    }

    d2 -=  (2.0 * (window->hacell[j2][c1] - window->hacell[j2][c2])) * 
      0.5 * AH[es] * (t22[c1] + t22[c2]);

    /* Step forward in time                                          */
    windat->nu1av[es] += windat->dt2d * (d1 + d2) / (d3 * depth);
  }
  debug_c(window, D_UA, D_HDIFF);
}

/* END hvisc_u1_2d()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hvisc_u1_2d_simple(geometry_t *window, /* Window geometry       */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon, /* Window constants      */
			double *tzp   /* Pointer to surface level    */
  )
{
  int e, ep, em, eu, ed, ee;  /* Edge coordinate / counter           */
  int es;                     /* Surface edge coordinates            */
  int c1, c2;                 /* Cells adjacent to e                 */
  int j1, j2;
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double d1;                  /* Constant for x diffusion            */
  double d2;                  /* Constant for y diffusion            */
  int *ctp;                   /* Cells to process vector             */
  int vcs;                    /* Number of cells to process          */

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
  /* Calculate the e1 horizontal diffusion                           */
  for (ee = 1; ee <= vcs; ee++) {
    e = ctp[ee];
    ep = window->ep[e];
    em = window->em[e];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];

    j1 = window->e2e[e][0];

    eu = (j1 == 1) ? e2e(window, e, 2) : e2e(window, e, 3);
    ed = (j1 == 1) ? e2e(window, e, 4) : e2e(window, e, 1);

    d1 = AH[e] * (windat->u1avb[ep] + windat->u1avb[em] -      
		  2.0 * windat->u1avb[e]) / (window->h2au1[es] *
					     window->h2au1[es]);

    d2 = AH[e] * (windat->u1avb[eu] + windat->u1avb[ed] -
		  2.0 * windat->u1avb[e]) / (window->h1au1[es] *
					     window->h1au1[es]);

    windat->nu1av[e] += windat->dt2d * (d1 + d2); 
  }

  debug_c(window, D_UA, D_HDIFF);
}

/* END hvisc_u1_2d_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the scaled horizontal mixing                       */
/*-------------------------------------------------------------------*/
void set_hdiff(geometry_t *window,      /* Window geometry           */
	       double *AH,              /* Mixing coefficient        */
	       double AH0               /* Constant mixing           */
	       )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double *h1, hm;                       /* Cell size                 */
  int cc, c, c2;                        /* Counters, cell locations  */
  int nv;                               /* Stencil size              */
  int *vec;                             /* Cells to process          */

  /* Assign pointers                                                 */
  h1 = window->h1au1;
  hm = wincon->hmean1;
  nv = window->n3_e1;
  vec = window->w3_e1;

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
		 int mode                 /* 1 = edges, 3 = centres  */
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
  int *m2d;

  if (smag == 1.0 && wincon->smagorinsky != 0.0) return;

  /* Assign pointers                                                 */
  scalef = wincon->diff_scale;
  if (mode == 1) {
    h1 = window->h1au1;
    hm = wincon->hmean1;
    nv = window->n3_e1;
    vec = window->w3_e1;
    AH0 = wincon->bsue1;
    m2d = window->m2de;
  } else if (mode == 3) {
    h1 = window->cellarea;
    hm = wincon->amean;
    nv = window->n3_t;
    vec = window->w3_t;
    AH0 = wincon->bkue1;
    m2d = window->m2d;
  }

  for (cc = 1; cc <= nv; cc++) {
    c = vec[cc];
    c2 = m2d[c];

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
		 int el                 /* Central coordinate        */
		 )
{
  window_t *windat = window->windat;    /* Window data               */
  win_priv_t *wincon = window->wincon;  /* Window constants          */
  double vh;                            /* Modified diffusion        */
  double rm, rmax, dmax;                /* Maximum diffusion         */
  double x0, y0, x1, y1;                /* Edge coordinates          */
  double sfact = 0.9;                   /* Safety factor             */
  double mfact = 4.0;                   /* Multiplication factor     */
  double *dist;                         /* Distance array            */
  int e, ee, es, j, cl;                 /* Counters, cell locations  */
  int ns;                               /* Stencil size              */
  int *st, sz, size = 3;                /* Stencil diameter          */
  int smag = 0;

  int code = ST_SQ3|ST_EDGE;

  if (!wincon->hdiff_f) return;

  /* Assign pointers                                                 */
  es = window->m2de[el];
  if (wincon->u1vh0 < 0.0) smag = 1;
  rm = fabs(wincon->u1vh0 * window->h1au1[es] / wincon->hmean1);
  rmax = 1.0 / (window->h1au1[es] * window->h1au1[es]);
  rmax = sfact / (4.0 * rmax * window->windat->dt);
  rm = min(mfact * rm, rmax);

  /* Get the stencil                                                 */
  cl = window->c2e[e][0];
  j = window->e2e[e][0];
  st = stencil(window, cl, &ns, code, j);

  /* Get the distance of stencil locations from the center           */
  sz = floor(size / 2);
  dist = d_alloc_1d(ns + 1);
  x0 = window->u1x[window->m2de[el]];
  y0 = window->u1y[window->m2de[el]];
  dmax = 0.0;
  for (ee = 0; ee < ns; ee++) {
    e = window->m2de[st[ee]];
    x1 = window->u1x[e];
    y1 = window->u1y[e];
    dist[ee] = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
    dmax = max(dmax, dist[ee]);
  }
  for (ee = 0; ee < ns; ee++) {
    e = es = window->m2d[st[ee]];
    while (e != window->zm1e[e]) {
      vh = (smag) ? wincon->u1vh[e] : windat->u1vhin[es];
      /* Get the maximum stable diffusion                            */
      wincon->u1vh[e] = rm + dist[ee] * (vh - rm) / (dmax + 1.0);
      e = window->zm1e[e];
    }
  }
  d_free_1d(dist);
  i_free_1d(st);
  windat->nalert[7]++;
}

/* END reset_hdiff()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Optimised horizontal mixing                                       */
/*-------------------------------------------------------------------*/
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
    smooth3(master, master->u1vh, geom->w3_e1, geom->n3_e1, geom->sze, -1);
  } else {
    if (u1vh > 0.0)
      vh1= u1vh;
    else
      vh1 = 0.01 * master->hmean1 * master->hmean1 / master->grid_dt;
    for (cc = 1; cc <= geom->n3_e1; cc++) {
      c = geom->w3_e1[cc];
      c2 = geom->m2d[c];
      if (flag == LINEAR)
	master->u1vh[c] = fabs(vh1 * geom->h1au1[c2] / master->hmean1);			       
      if (flag == NONLIN)
	master->u1vh[c] = fabs(vh1 * geom->h1au1[c2] * geom->h1au1[c2] / 			       
			       (master->hmean1 * master->hmean1));
    }
  }
}

/* END reset_hor_diff()                                             */
/*------------------------------------------------------------------*/
