/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/momentum/advect.c
 *  
 *  Description: Momentum advection schemes for 2D and 3D modes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: advect.c 6075 2019-02-08 04:11:42Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

double porus_plate_e1(geometry_t *window, window_t *windat, win_priv_t *wincon, int c);
double porus_plate_e2(geometry_t *window, window_t *windat, win_priv_t *wincon, int c);
int cg;
#define SMALL 1e-10             /* Minimum value for velocity bounds */

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U1 ADVECTION ROUTINES                                             */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Computes relative vorticity on vertices, normalized relative and  */
/* planetary vorticity on vertices and edges.                        */
/*-------------------------------------------------------------------*/
int nonlin_coriolis_3d(geometry_t *window,  /* Window geometry       */
		       window_t *windat,    /* Window data           */
		       win_priv_t *wincon   /* Window constants      */
		       )  
{
  int cc, c, cs, c1, c2;    /* Cell counters                         */
  int ee, e, es, eoe;       /* Edge counters                         */
  int vv, v, vs, v1, v2;    /* Vertex counters                       */
  int n, j;                 /* Genereal counters                     */
  double d1, d2, d3;        /* Dummies                               */
  double iarea;             /* 1 / area                              */
  double *dzv = wincon->w1;
  double *vel;              /* Velocity to use in spatial gradients  */
  double dtu, dtm;
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  int itermax = 20;         /* Maximum number of substeps            */
  char tp = 'n';            /* Sub step code                         */
  double mvel;              /* Maximum velocity                      */
  int ii, jj, kk;
  int vc, *cells;

  /*-----------------------------------------------------------------*/
  /* Initialize                                                      */
  windat->dtu1 = windat->dtf;
  vel = windat->u1;

  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on  */
  /* the u1 face to calculate the maximum timesteps.                 */
  dtm = dtu = windat->dtf;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF)) {
    for (ee = 1; ee <= vc; ee++) {
      int c1s, c2s;
      double vlc;
      double sf = 0.8;          /* Fraction of sub-timestep used     */
      double minval = 1e-10;    /* Minimum value for velocity        */
      double *u2au1 = wincon->w7;
      e = cells[ee];            /* Wet cell to process               */
      es = window->m2de[e];     /* 2D cell corresponding to 3D cell  */
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];

      /* Minimum time-step due to x velocity                         */
      vlc = vel[e];
      if (fabs(vlc) < minval)
        vlc = minval;
      d1 = sf * fabs(window->h2au1[es] / vlc);

     if (d1 < dtm) {
        dtm = d1;
        tp = 'u';
        mvel = vlc;
        ii = c1s;
        jj = es;
        kk = window->s2k[c1];
      }
    }
  }
  if (dtm != dtu) {
    e = ceil((windat->dtb + dtu) / (2 * dtm));
    dtm = (windat->dtb + dtu) / (2.0 * (double)e);
    dtu = 2.0 * dtm;
    hd_warn
      ("Sub-time stepping in u1 (%c=%6.3f, h1=%f) at %8.3f days (cs=%3d es=%3d %3d): dt=%5.2f [%f %f]\n",
       tp, mvel, window->h2au1[es], windat->t / 86400.0, ii, jj, kk, dtm,
       window->u1x[es], window->u1y[es]);
    if (dtm == 0.0 || (int)(dtu / dtm) > itermax) {
      hd_quit_and_dump
	("advect_u1_3d: maximum number of sub-steps (%d) exceeded.\n",
	 itermax);
      return(1);
    }
  } else
    dtu = dtm + windat->dtb;
  windat->dtu1 = dtm;

  trem = windat->dt;
  tremf = windat->dtf;

  /* Loop over the time interval until no time remains               */
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Compute the relative vorticity at vertices. This is         */
      /* achieved by computing the circulation around the perimeter  */
      /* of the dual grid normalized by the area of the dual grid    */
      /* (e.g. Apel Eq. 6.53). This corresponds to Eq. 23 in Ringler */
      /* et al (2010).                                               */
      memset(windat->circ, 0, window->szv * sizeof(double));
      memset(windat->rvor, 0, window->szv * sizeof(double));
      for (vv = 1; vv <= window->a3_e2; vv++) {
	v = window->w3_e2[vv];
	vs = window->m2dv[v];
	iarea = 1.0 / window->dualarea[vs];
	for (ee = 1; ee <= window->nve[vs]; ee++) {
	  e = window->v2e[v][ee];
	  if (e) {
	    es = window->m2de[e];
	    d2 = window->h2au1[es] * vel[e];
	    windat->circ[v] += window->eSv[ee][v] * d2;
	    windat->rvor[v] += window->eSv[ee][v] * iarea * d2;
	  }
	}
      }

      /*-------------------------------------------------------------*/
      /* Set relative vorticity along solid boudaries equal to zero. */
      /* This equates to a free slip boundary condition (Ringler et  */
      /* al, (2013) pA.39, Sect. A.4.                                */
      for (vv = window->a3_e2 + 1; vv <= window->n3_e2; vv++) {
	v = window->w3_e2[vv];
	windat->rvor[v] = 0.0;
      }
      
      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* vertices (e.g. Eq. A.44 Ringler et al, 2013). The sum of    */
      /* these is the potential vorticity.                           */
      memset(dzv, 0, window->szv * sizeof(double));
      for (vv = 1; vv <= window->n3_e2; vv++) {
	v = window->w3_e2[vv];
	vs = window->m2dv[v];
	d1 = 1.0 / window->dualarea[vs];
	for (cc = 1; cc <= window->nvc[vs]; cc++) {
	  dzv[v] += wincon->dz[window->v2c[v][cc]] * window->dualareap[vs][cc];
	}
	dzv[v] *= d1;
	windat->nrvor[v] = windat->rvor[v] / dzv[v];
	windat->npvor[v] = windat->fv[vs] / dzv[v];
      }

      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* edges. This is the quantity used in Eq. A.39, Ringler et    */
      /* al, (2013).                                                 */
      memset(windat->nrvore, 0, window->sze * sizeof(double));
      memset(windat->npvore, 0, window->sze * sizeof(double));
      for (ee = 1; ee <= window->n3_e1; ee++) {
	e = window->w3_e1[ee];
	v1 = window->e2v[e][0];
	v2 = window->e2v[e][1];
	windat->nrvore[e] = 0.5 * (windat->nrvor[v1] + windat->nrvor[v2]);
	windat->npvore[e] = 0.5 * (windat->npvor[v1] + windat->npvor[v2]);
      }
      
      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* cell centres.                                               */
      memset(windat->nrvorc, 0, window->sgsiz * sizeof(double));
      for(cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	iarea = 1.0 / window->cellarea[cs];
	for (n = 1; n <= window->npe[cs]; n++) {
	  j = window->vIc[n][cs];
	  v = window->c2v[n][c];
	  vs = window->m2dv[v];
	  windat->nrvorc[c] += window->dualareap[vs][j] * windat->nrvor[v] * iarea;
	}
      }

      /*-------------------------------------------------------------*/
      /* Compute the kinetic energy at cell centres. Used in the     */
      /* decomposition of momentum advection, e.g. Eq. 4, Ringler et */
      /* al, 2010.                                                   */
      /* Divergence is given by Eq. 21, Ringler et al, 2010.         */
      memset(windat->kec, 0, window->szc * sizeof(double));
      memset(windat->div, 0, window->szc * sizeof(double));
      for(cc = 1; cc <= window->a3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	d1 = 1.0 / window->cellarea[cs];
	for (ee = 1; ee <= window->npe[cs]; ee++) {
	  e = window->c2e[ee][c];
	  es = window->m2de[e];
	  d2 = window->h1au1[es] * vel[e] * d1;
	  windat->div[c] += window->eSc[ee][cs] * d2;
	  windat->kec[c] += 0.25 * d2 * window->h2au1[es] * vel[e];
	}
      }

      /*-------------------------------------------------------------*/
      /* Compute the nonlinear Coriolis tendency. The tangential     */
      /* velocity is computed using the primal to dual mesh mapping  */
      /* (Eq. 24, Ringler et al, 2010) and a mean of absolute        */
      /* vorticity at the edge in question. The gradient of kinetic  */
      /* energey is also added (see Eq. A.39, Ringler et al, 2013).  */
      if (wincon->u1_f & ADVECT) {
	memset(windat->nrvore, 0, window->sze * sizeof(double));
	memset(windat->kec, 0, window->sgsiz * sizeof(double));
      }
      if (wincon->u1_f & CORIOLIS)
	memset(windat->npvore, 0, window->sze * sizeof(double));

      for (ee = 1; ee <= vc; ee++) {
	e = cells[ee];
	es = window->m2de[e];
	c1 = window->e2c[e][0];
	c2 = window->e2c[e][1];
	d1 = 1.0 / window->h2au1[es];
	d2 = 0.0; 
	for (n = 1; n <= window->nee[es]; n++) {
	  eoe = window->eSe[n][e];
	  if (!eoe) continue;
	  /* Ringler et, al. (2010) Eq. 49                           */
	  /*d3 = 0.5 * (windat->nrvore[e] + windat->npvore[e] + windat->nrvore[eoe] + windat->npvore[eoe]);*/
	  d3 = wincon->pv_calc(windat, e, eoe, window->e2v[e][0], window->e2v[e][1]);
	  d2 += window->wAe[n][e] * vel[eoe] * d3 * windat->dzu1[eoe];
	}
	d3 = d2 - (windat->kec[c1] - windat->kec[c2]) * d1;
	wincon->u1inter[es] += (d3 * dtm * windat->dzu1[e]);
	windat->nu1[e] += d3 * dtu;
      }
      if (wincon->u1_f & ADVECT)
	memset(wincon->u1inter, 0, window->szeS * sizeof(double));

      /* Get the timestep for the next iteration.                    */
      /* Decrement the time remaining. This must be computed for     */
      /* both the leapfrog step (trem) used to update velocity, and  */
      /* the forward step (tremf) used to integrate u1inter[] over   */
      /* forward time-step windat->dtf.                              */
      trem -= dtu;
      tremf -= dtm;
      if (dtm)
	dtu = dtm * windat->dt / windat->dtf;
      if (trem < dtu)
        dtu = trem;
      if (tremf < dtm)
        dtm = tremf;
      if(!wincon->sigma)
	vel = windat->nu1;
    }
  }
  return(0);
}

/* END nonlin_coriolis_3d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Energetically neutral formulation of nonlinear Coriolis force:    */
/* Ringler et al, (2010) Eq.49.                                      */
/*-------------------------------------------------------------------*/
double pv_energy_neutral(window_t *windat,    /* Window data         */
			 int e,               /* Edge coordinate     */
			 int eoe,             /* Edge coordinate     */
			 int v1,              /* Vertex coordinate   */
			 int v2               /* Vertex coordinate   */
			 )
{
  double v;
  v = 0.5 * (windat->nrvore[e] + windat->npvore[e] + windat->nrvore[eoe] + windat->npvore[eoe]);
  return(v);
}

/* END pv_energy_neutral()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Potential enstrophy conserving formulation of nonlinear Coriolis  */
/* force: Ringler et al, (2010) Eq.71.                               */
/*-------------------------------------------------------------------*/
double pv_enstrophy_conserve(window_t *windat, /* Window data        */
			     int e,            /* Edge coordinate    */
			     int eoe,          /* Edge coordinate    */
			     int v1,           /* Vertex coordinate  */
			     int v2            /* Vertex coordinate  */
			     )
{

  return(windat->nrvore[e] + windat->npvore[e]);
}

/* END pv_enstrophy_conserve()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Potential enstrophy dissapating formulation of nonlinear Coriolis */
/* force: Ringler et al, (2010) Eq.73.                               */
/*-------------------------------------------------------------------*/
double pv_enstrophy_dissipate(window_t *windat, /* Window data       */
			      int e,            /* Edge coordinate   */
			      int eoe,          /* Edge coordinate   */
			      int v1,           /* Vertex coordinate */
			      int v2            /* Vertex coordinate */
			      )
{
  double ret;

  ret = (windat->u2[e] >=0) ? (windat->nrvor[v2] + windat->npvor[v2]) :
    (windat->nrvor[v1] + windat->npvor[v1]);
  return(ret);
}

/* END pv_enstrophy_dissipate()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/* Old code - no longer used for unstructured grids since we use     */
/* nonlin_coriolis_3d() instead. We retain this code for a possible  */
/* future tensor approach.                                           */
/*-------------------------------------------------------------------*/
int advect_u1_3d(geometry_t *window,  /* Processing window           */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon   /* Window constants            */
		 )
{
  int e, ee, es, ep, em;    /* Edge coordinate / counter             */
  int e1, e2, e1s, e2s;     /* Edge coordinate / counter             */
  int c, cc;                /* Cell coordinate / counter             */
  int c1, c2, c1s, c2s;     /* Cell coordinate / counter             */
  int cs, cb, ems;          /* Surface / bottom sparse coordinate    */
  int j, npe;
  int zp1, zm1, zp2;        /* Sparse coordinate at k+1, k-1         */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double u1val;             /* Cell centered u1 value                */
  double u2val;             /* Cell cornered u2 value                */
  double wval;              /* Face centered vertical velocity       */
  double *u2au1;            /* u2 value at the e1 face               */
  double h1;                /* Face centered grid spacing            */
  double h1av, h2av;        /* Mean cell spacings x and y directions */
  double dzav;              /* Mean cell spacings in the z direction */
  double m, m1, m2;         /* Dummy variables                       */
  double d1, d2, d3;        /* Dummy variables                       */
  double *vel;              /* Velocity to use in spatial gradients  */
  double *area;             /* Cell area                             */
  double *ush1h2;           /* u1 * u1 * h1 * h2                     */
  double *u1u2hs;           /* u1 * u2 * h1 * h1                     */
  double *wu1;              /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *vface;            /* wvel centered on the grid face        */
  double *cn;               /* Vertical courant number               */
  double *velbuf;           /* u1 buffer for vertical flux           */
  double *dzu1;             /* Pointer for thickness                 */
  double *div;              /* Momentum divergence                   */
  double midx;              /* Depth at the e1 face                  */
  double top;               /* Surface level at the cell face        */
  double sf = 0.8;          /* Fraction of sub-timestep used         */
  double wsf = 0.8;         /* Fraction of vertical subtimestep used */
  double minval = 1e-10;    /* Minimum value for velocity            */
  int *cells;               /* Cells to process                      */
  int itermax = 20;         /* Maximum number of substeps            */
  char tp = 'n';            /* Sub step code                         */
  int ii = 0, jj = 0, kk = 0;   /* Sub step coordinates              */
  double vs = 0, vlc = 0;       /* Sub step velocity                 */
  double *wq;

  if (!wincon->nonlinear)
    return(0);

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  ush1h2 = wincon->w1;
  u1u2hs = wincon->w2;
  wu1 = wincon->w3;
  wvel = wincon->w4;
  velbuf = wincon->w5;
  div = wincon->w6;
  area = window->cellarea;
  u2au1 = wincon->w7;
  vel = windat->u1;
  wq = wincon->d1;
  if (wincon->momsc & ADVECT_FORM) {
    dzu1 = wincon->w8;
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      dzu1[e] = 1.0;
    }
  } else
    dzu1 = windat->dzu1;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required.                             */
  /* The advective terms can create large velocities due to momentum */
  /* divergence if the sea level lies just above a layer (i.e. a     */
  /* layer) in a cell and just below the layer in the adjacent cell. */
  /* Here dzu1 is thin at cellc, ush1h2 may be large at c and small  */
  /* at the neighboring, hence (ush1h2[c]-ush1h2[cn])/dzu1[c] can be */
  /* large.                                                          */
  vcs = wincon->vcs;
  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }
  if(wincon->dolin_u1) {
    linear_bdry_cell(window, windat, wincon, cells, vcs,
		     wincon->linmask_u1, wincon->i5);
    cells = wincon->s4;
    vc = wincon->acl;
    vcs = wincon->aclS;
  }
  set_surf_cells(window, windat, wincon, 1);

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on  */
  /* the u1 face to calculate the maximum timesteps. If wincon->stab */
  /* = SUB_STEP_NOSURF then the maximum sub-step is not calculated   */
  /* in the surface layer for vertical velocity. This avoids sub-    */
  /* stepping when the surface elevation is only slightly greater    */
  /* than a layer level (thin surface layers) and moderately large   */
  /* vertical velocities exist. This condition may occur often and   */
  /* increase the run time ratio (slower execution) due to frequent  */
  /* sub-stepping. However, the model may go unstable in the surface */
  /* layer in this case, and wincon->hmin may need to be increased.  */
  dtm = dtu = windat->dtf;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF)) {
    for (ee = 1; ee <= vc; ee++) {
      e = cells[ee];            /* Wet cell to process               */
      es = window->m2de[e];      /* 2D cell corresponding to 3D cell  */
      zp1 = window->zp1e[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      c1s = window->m2d[c1];
      c2s = window->m2d[c2];
      zp1 = window->zp1[c1];
      zp2 = window->zp1[c2];

      /* Minimum time-step due to x velocity                         */
      vlc = windat->u1[e];
      if (ee <= vcs) wq[es] = dtu * vlc / window->h2au1[es];
      if (fabs(vlc) < minval)
        vlc = minval;
      d1 = sf * fabs(window->h2au1[es] / vlc);
      if (d1 < dtm) {
        dtm = d1;
        tp = 'u';
        vs = vlc;
        ii = window->s2i[c1];
        jj = window->s2j[c1];
        kk = window->s2k[c1];
      }
      /* Minimum time-step due to y velocity                         */
      vlc = u2au1[e];
      if (ee <= vcs) wq[es] = max(wq[es], dtu * vlc / window->h1au1[es]);
      if (fabs(vlc) < minval)
        vlc = minval;
      d2 = sf * fabs(window->h1au1[es] / vlc);
      if (d2 < dtm) {
        dtm = d2;
        tp = 'v';
        vs = vlc;
        ii = window->s2i[c1];
        jj = window->s2j[c1];
        kk = window->s2k[c1];
      }
      /* Minimum time-step due to z velocity                         */
      vlc = (ee <= vcs) ? windat->wtop[c1s] + windat->wtop[c2s] : windat->w[zp1] + windat->w[zp2];
      vlc = 0.25 * (vlc + windat->w[c1] + windat->w[c2]);
      if (ee <= vcs) wq[es] = max(wq[es], dtu * vlc / (windat->dzu1[e] * wincon->mdx[es]));
      if (wincon->stab & SUB_STEP_NOSURF && ee <= vcs)
        vlc = 0.0;
      if (fabs(vlc) < minval)
        vlc = minval;
      /* SIGMA : Multiply by depth */
      d3 =
        wsf * fabs(max(windat->dzu1[e] * wincon->mdx[es], wincon->hmin) / vlc);
      if (d3 < dtm) {
        dtm = d3;
        tp = 'w';
        vs = vlc;
        ii = window->s2i[c1];
        jj = window->s2j[c1];
        kk = window->s2k[c1];
      }
    }
  }
  if (dtm != dtu) {
    e = ceil((windat->dtb + dtu) / (2 * dtm));
    dtm = (windat->dtb + dtu) / (2.0 * (double)e);
    dtu = 2.0 * dtm;

    hd_warn
      ("Sub-time stepping in u1 (%c=%6.3f) at %8.3f days (%3d %3d %3d): dt=%5.2f\n",
       tp, vs, windat->t / 86400.0, ii, jj, kk, dtm);
    if ((int)(dtu / dtm) > itermax) {
      hd_quit_and_dump
	("advect_u1_3d: maximum number of sub-steps (%d) exceeded.\n",
	 itermax);
      return(1);
    }
  }
  else
    dtu = dtm + windat->dtb;
  windat->dtu1 = dtm;

  if (wincon->u1_f & ADVECT)
    return(0);

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity. Note that w is cell centered, hence the centered      */
  /* bottom coordinate must be used. Auxiliary cells are not defined */
  /* for cell centers on western and southern boundaries, so the _t  */
  /* vectors cannot be used here and the cell center must be found   */
  /* by mapping downwards from bot_e1 until self-mapping is located. */
  memcpy(wincon->w4, windat->w, window->szc * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_t[cc]; /* 2D coordinate                     */
    cb = window->bot_t[cc];    /* Bottom e1 coordinate              */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs]; /* w : top set to wtop             */
      c = window->zm1[c];
    }
    /* Set the bottom boundary condition. Find the cell centered     */
    /* bottom coordinate and set w to wbot.                          */
    c = cb;                       /* e1 face centered bottom cell    */
    while (c != window->zm1[c])
      c = window->zm1[c];         /* Cell centered sediment location */
    cb = window->zp1[c];          /* cell centered bottom location   */
    wvel[cb] = windat->wbot[cs];  /* w : bottom set to wbot          */
  }
  /* Set wvel at the surface for ghost cells. If window partitions   */
  /* lie next to ghost cells then these are not captured in the loop */
  /* above but are required for face centered vertical velocity.     */
  /* Note that wvel[cb] always = 0 and is not explicitly set.        */
  for (cc = window->a2_t + 1; cc <= window->n2_t; cc++) {
    cs = c = window->w2_t[cc];
    top = windat->eta[cs];
    while (c != window->zm1[c] && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs];
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the first timestep for leapfrog                             */
  trem = windat->dt;
  tremf = windat->dtf;
  /* Loop over the time interval until no time remains               */
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Initialise                                                  */
      memset(ush1h2, 0, window->sze * sizeof(double));
      memset(u1u2hs, 0, window->sze * sizeof(double));
      memset(wu1, 0, window->sze * sizeof(double));
      if (wincon->momsc & ZERO_DRYK)
	set_surf_cond(window, vel, 0.0, window->sur_e1, window->b2_e1);

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres (u1*u1*h1*h2*dz).    */
      /* This must be calculated at ghost cells also since these     */
      /* locations will give a non-zero cell centered flux.          */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  cs = window->e2c[es][0];
	  j = window->e2e[e][1];
	  ep = e2e(window, e, j);

          u1val = 0.5 * (vel[e] + vel[ep]);
          /* SIGMA : Multiply by the cell centre depth               */
          if (u1val > 0)
            m = dzu1[e] * vel[e] * wincon->mdx[es];
          else
            m =
              dzu1[ep] * vel[ep] * wincon->mdx[window->m2de[ep]];
          ush1h2[e] = m * u1val * area[cs];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  cs = window->e2c[es][0];
	  j = window->e2e[e][1];
	  ep = e2e(window, e, j);
          u1val = 0.5 * (vel[e] + vel[ep]);

          /* SIGMA : Multiply by the cell centre depth               */
          m = 0.5 * (dzu1[ep] * vel[ep] *
                   wincon->mdx[window->m2de[ep]] +
                   dzu1[e] * vel[e] * wincon->mdx[es]);
          ush1h2[e] = m * u1val * area[cs];
	}
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  cs = window->e2c[es][0];
	  j = window->e2e[e][1];
	  ep = e2e(window, e, j);
          vface[c] = 0.5 * (vel[e] + vel[ep]);
          cn[e] = vface[e] * dtu / window->h1acell[cs];
          /* The interpolated velocity must be shifted by one cell   */
          /* when using van_leer_do() since this routine is set up   */
          /* to interpolate onto the cell face and we need to        */
          /* interpolate to the cell center. Use wu1 as a dummy.     */
          wu1[e] =
            dzu1[ep] * vel[ep] * wincon->mdx[window->m2de[ep]];
        }
	/* Ghost cells 2 cells in */
        for (ee = window->b3_e1 + 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  j = window->e2e[e][0];
	  em = e2e(window, e, j);
          wu1[em] = dzu1[e] * vel[e] *
            wincon->mdx[window->m2de[e]];
	}

        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
	for (ee = 1; ee <= window->n3_e1; ee++) {
	  int cp, cm;
	  e = window->w3_e1[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
	  j = window->e2e[e][0];
	  cm = window->c2c[j][c2];
	  j = window->e2e[e][1];
	  cp = window->c2c[j][c1];
	  van_leer_do(ush1h2, wu1, vface, cn, e, c1, cp, c2, cm);
	}

        /* Multiply by centered velocity and shift one cell down     */
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  cs = window->e2c[es][0];
          ush1h2[e] *= vface[e] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners (u1*u2*h1*h1*dz).    */
      /* This must be calculated at ghost cells also since these     */
      /* locations will give a non-zero flux at ouside corners.      */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
	  c = window->m2d[window->e2c[e][0]];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ete(window, e, &e1, &e2);

          h1 = 0.5 * (window->h1au1[window->m2de[e1]] + window->h1au1[window->m2de[e2]]);
          u2val = 0.5 * (windat->u1[e1] * wincon->mdx[window->m2de[e1]] +
                         windat->u1[e2] * wincon->mdx[window->m2de[e2]]);
          if (u2val > 0)
            m = dzu1[em] * vel[em];
          else
            m = dzu1[e] * vel[e];
          u1u2hs[e] = m * h1 * h1 * u2val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  c = window->e2c[es][0];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ems = window->m2de[em];
	  ete(window, e, &e1, &e2);
	  e1s = window->m2de[e1];
	  e2s = window->m2de[e2];

          h1 = 0.5 * (window->h1au1[e1s] + window->h1au1[e2s]);
          h1av = window->h1au1[e2s] / (window->h1au1[e1s] + window->h1au1[e2s]);
          h2av = window->h1au1[ems] / (window->h1au1[ems] + window->h1au1[es]);
          /* SIGMA : Multiply by the cell corner depth               */
          u2val = h1av * (windat->u1[e1] * wincon->mdx[e1s] -
                          windat->u1[e2] * wincon->mdx[e2s]) +
            windat->u1[e2] * wincon->mdx[e2s];
          m1 = dzu1[em] * vel[em];
          m2 = dzu1[e] * vel[e];
          m = h2av * (m2 - m1) + m1;
          u1u2hs[e] = m * h1 * h1 * u2val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sze * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  c = window->e2c[es][0];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ems = window->m2de[em];
	  ete(window, e, &e1, &e2);
	  e1s = window->m2de[e1];
	  e2s = window->m2de[e2];
          vface[e] = 0.5 * (windat->u1[em] + windat->u1[e]);
          cn[e] = 2.0 * vface[e] * dtu / (window->h1au1[e1s] + window->h1au1[e2s]);
          wu1[e] = dzu1[e] * vel[e] * wincon->mdx[es];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
	for (ee = 1; ee <= window->n3_e1; ee++) {
	  int cp, cm;
	  e = window->w3_e1[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
	  j = window->e2e[e][0];
	  cm = window->c2c[j][c2];
	  j = window->e2e[e][1];
	  cp = window->c2c[j][c1];
	  van_leer_do(u1u2hs, wu1, vface, cn, e, c1, cp, c2, cm);
	}

        /* Multiply by velocity and shift one cell down              */
        for (ee = 1; ee <= window->n3_e1; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  cs = window->e2c[es][0];
          u1u2hs[e] *= vface[e] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the vertical advective fluxes (u1*w).                   */
      /* First set the bottom boundary condition for u1 and w. Due   */
      /* to the stagger on velocity cells, there may not exist a     */
      /* separate sediment ghost cell for velocity at locations      */
      /* where there is a step in the bathymetry on eastern faces.   */
      /* In these situations the cell center is wet and the cell     */
      /* face is dry, and the the 'sediment' cell for u1 velocity    */
      /* is associated with a wet cell whose face is a lateral ghost */
      /* cell for u1 velocity. Setting a vertical no gradient at     */
      /* this cell will actually also set a non-zero normal velocity */
      /* at this lateral ghost cell. A buffer is used here for u1    */
      /* velocity to avoid this.                                     */
      if (wincon->momsc & ZERO_DRYK)
	set_surf_cond(window, vel, -1.0, window->sur_e1, window->b2_e1);
      memcpy(velbuf, vel, window->sze * sizeof(double));
      for (ee = 1; ee <= window->a2_e1; ee++) {
        e = window->bot_e1[ee]; /* 3D bottom coordinate              */
        zm1 = window->zm1e[e];  /* Sediment coordinate               */
        velbuf[zm1] = vel[e];   /* u1 : no-gradient across bottom    */
      }

      /* Surface vertical advective flux                             */
      for (ee = 1; ee <= wincon->vcs; ee++) {
        e = cells[ee];
        es = window->m2de[e];
	c1 = window->e2c[es][0];
	c2 = window->e2c[es][1];
        wu1[e] = 0.5 * vel[e] * (windat->wtop[c1] + windat->wtop[c2]);
      }

      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        /* Interior layers. Note that wincon->w4 has been set up to  */
        /* contain vertical velocity at all wet cells and the        */
        /* surface vertical velocity, wtop[], at all cells above the */
        /* free surface. The vertical flux is stored staggered in    */
        /* wu1 so that the the any given sparse coordinate contains  */
        /* the vertical advective flux at the UPPER FACE of that     */
        /* cell (as opposed to vertical velocity which resides at    */
        /* lower face.                                               */
        for (ee = 1; ee <= vc; ee++) {
          e = cells[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
          zm1 = window->zm1e[e];
          wval = 0.5 * (wvel[c1] + wvel[c2]);
          if (wval > 0)
            wu1[zm1] = velbuf[zm1] * wval;
          else
            wu1[zm1] = velbuf[e] * wval;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (ee = 1; ee <= vc; ee++) {
          e = cells[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
          zm1 = window->zm1e[e];
          wval = 0.5 * (wvel[c1] + wvel[c2]);
	  dzav = windat->dzu1[zm1] / (windat->dzu1[zm1] + windat->dzu1[e]);
	  wu1[zm1] = wval * (dzav * (velbuf[e] - velbuf[zm1]) + velbuf[zm1]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sze * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (ee = 1; ee <= window->n3_t; ee++) {
          e = window->w3_e1[ee];
          es = window->m2de[e];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
          dzav = 0.5 * (windat->dzu1[e] + windat->dzu1[window->zm1e[e]]);
          vface[e] = 0.5 * (wvel[c1] + wvel[c2]);
          cn[e] = vface[e] * dtu / (dzav * wincon->mdx[es]);
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wvel.                                         */
	for (cc = 1; cc <= wincon->vc; cc++) {
	  int zm2;
	  c = cells[cc];
	  zp1 = window->zp1[c];
	  zm1 = window->zm1[c];
	  zm2 = window->zm1[zm1];
	  van_leer_do(wvel, velbuf, vface, cn, c, c, zp1, zm1, zm2);
	}

        /* Multiply by vertical velocity and shift one cell down     */
        for (ee = 1; ee <= vc; ee++) {
          e = cells[ee];
          zm1 = window->zm1e[e];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
          wval = 0.5 * (wvel[c1] + wvel[c2]); 
	  wu1[zm1] = vface[e] * wval;
        }
      }

      /* Bottom vertical advective flux. Note: This formulation is   */
      /* consistent with the original meco, however, taking this     */
      /* loop out makes the formulation consistent with the surface  */
      /* flux calculation, i.e. if many layers separate the bottom   */
      /* cells either side of a face (the bottom is much deeper on   */
      /* one side of a face than the other), then to set the bottom  */
      /* velocity at the face equal to the mean wbot's either side   */
      /* of the face is not a good approximation - better to take a  */
      /* mean in the horizontal.                                     */
      for (ee = 1; ee <= window->v2_e1; ee++) {
        e = window->bot_e1[ee]; /* Bottom coordinate                 */
        es = window->w2_e1[ee];
	c1 = window->e2c[es][0];
	c2 = window->e2c[es][1];
        zm1 = window->zm1[e];
        wu1[zm1] =
          0.5 * velbuf[e] * (windat->wbot[c1] + windat->wbot[c2]);
      }

      /*-------------------------------------------------------------*/
      /* Get the momentum divergence and metric terms.               */
      for (ee = 1; ee <= vc; ee++) {
        e = cells[ee];
        es = window->m2de[e];
	j = window->e2e[e][0];
	em = e2e(window, e, j);
	j = window->e2e[e][0];
	j = (j <= npe/2) ? jp(j, npe) : jm(j, npe);
	ep = e2e(window, e, j);
        zm1 = window->zm1e[e];
        midx = wincon->mdx[es];

        /* SIGMA : Multiply by cell depth                            */
        d1 =  wincon->u1c1[es] * (ush1h2[e] - ush1h2[em] + 
				  u1u2hs[ep] - u1u2hs[e]) /
	  max(dzu1[e], wincon->hmin) +
          wincon->u1c3[es] * vel[e] * vel[e] * midx +
          wincon->u1c4[es] * u2au1[e] * u2au1[e] * midx;

	if (wincon->porusplate)
	  d1 += porus_plate_e1(window, windat, wincon, e);

        /* Integrate the e1 dispersion term vertically and with time */
        /* (i.e. over the sub-timestepping; windat->dtf).            */
        wincon->u1inter[es] += (d1 * dtm * windat->dzu1[e]);

        /* Add the vertical advective flux. Note : the difference of */
        /* this flux is wu1[c]-wu1[zm1] due to the stagger imposed,  */
        /* not wu1[zp1]-wu1[c].                                      */
	if (!(wincon->momsc & WIMPLICIT))
	  d1 -= (wu1[e] - wu1[zm1]) / max(windat->dzu1[e], wincon->hmin);

        /* Add to new u1 value */
	div[e] = d1 * dtu;
        /*windat->nu1[c] += d1 * dtu;*/
      }

      /* Filter the divergence if required                           */
      if (wincon->filter & ADVECT) {
	/*not_included("momentum tensor Shapiro filtering");*/
	vel2D_lbc(div, window->nbpte1, window->nbe1,
		  window->bpte1, window->bine1, wincon->slip);
	shapiro(window, div, velbuf, cells, vc, 1, 1, XDIR);
      }

      /* Get the u1 updated solution                                 */
      for (ee = 1; ee <= vc; ee++) {
        e = cells[ee];
        windat->nu1[e] += div[e];
      }

      /* Get the timestep for the next iteration.                    */
      /* Decrement the time remaining. This must be computed for     */
      /* both the leapfrog step (trem) used to update velocity, and  */
      /* the forward step (tremf) used to integrate u1inter[] over   */
      /* forward time-step windat->dtf.                              */
      trem -= dtu;
      tremf -= dtm;
      if (dtm)
	dtu = dtm * windat->dt / windat->dtf;
      if (trem < dtu)
        dtu = trem;
      if (tremf < dtm)
        dtm = tremf;
      if(!wincon->sigma)
	vel = windat->nu1;
    }
  }
  return(0);
}

/* END advect_u1_3d()                                          (slf) */
/*-------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U1AV ADVECTION ROUTINES                                           */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Computes relative vorticity on vertices, normalized relative and  */
/* planetary vorticity on vertices and edges.                        */
/*-------------------------------------------------------------------*/
int nonlin_coriolis_2d(geometry_t *window,  /* Window geometry       */
		       window_t *windat,    /* Window data           */
		       win_priv_t *wincon   /* Window constants      */
		       )
{
  int cc, c, cs, c1, c2;    /* Cell counters                         */
  int ee, e, es, eoe;       /* Edge counters                         */
  int vv, v, vs, v1, v2;    /* Vertex counters                       */
  int n, j;                 /* Genereal counters                     */
  double d1, d2, d3;        /* Dummies                               */
  double iarea;             /* 1 / area                              */
  double *dzv = wincon->w1;
  double dtu, dtm, lr;
  double *vel;              /* Velocity to use in spatial gradients  */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double *prv = windat->rv_drvdt;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  dtu = windat->dtu1 / d1;          /* 2D forward time-step          */
  vel = windat->u1av;

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  if(windat->dtu1 == windat->dtf)
    dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping  */
  else
    dtu *= 2.0;           /* 2D leapfrog time-step with sub-stepping */
  lr = windat->dtf2 / windat->dt2d;
  dtm = dtu * lr;
  trem = windat->dt2d;
  tremf = windat->dtf2;
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Compute the relative vorticity at vertices. This is         */
      /* achieved by computing the circulation around the perimeter  */
      /* of the dual grid normalized by the area of the dual grid    */
      /* (e.g. Apel Eq. 6.53). This corresponds to Eq. 23 in Ringler */
      /* et al (2010).                                               */
      memset(windat->circ, 0, window->szvS * sizeof(double));
      memset(windat->rvor, 0, window->szvS * sizeof(double));
      for (vv = 1; vv <= window->a2_e2; vv++) {
	v = window->w2_e2[vv];
	vs = window->m2dv[v];
	iarea = 1.0 / window->dualarea[vs];
	for (ee = 1; ee <= window->nve[vs]; ee++) {
	  e = window->v2e[v][ee];
	  if (e) {
	    d2 = window->h2au1[e] * windat->u1av[e];
	    windat->circ[v] += window->eSv[ee][vs] * d2;
	    windat->rvor[v] += window->eSv[ee][vs] * iarea * d2;
	  }
	}
      }

      /*-------------------------------------------------------------*/
      /* Set relative vorticity along solid boudaries equal to zero. */
      /* This equates to a free slip boundary condition (Ringler et  */
      /* al, (2013) pA.39, Sect. A.4.                                */
      for (vv = window->a2_e2 + 1; vv <= window->n2_e2; vv++) {
	v = window->w2_e2[vv];
	windat->rvor[v] = 0.0;
      }

      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* vertices (e.g. Eq. A.44 Ringler et al, 2013). The sum of    */
      /* these is the potential vorticity.                           */
      memset(dzv, 0, window->szvS * sizeof(double));
      for (vv = 1; vv <= window->n2_e2; vv++) {
	v = window->w2_e2[vv];
	vs = window->m2dv[v];
	d1 = 1.0 / window->dualarea[vs];
	for (cc = 1; cc <= window->nvc[vs]; cc++) {
	  cs = window->v2c[v][cc];
	  d2 = windat->eta[cs] - window->botz[cs];
	  dzv[v] += d2 * window->dualareap[vs][cc];
	}
	dzv[v] /= window->dualarea[vs];
	windat->nrvor[v] = windat->rvor[v] / dzv[v];
	windat->npvor[v] = windat->fv[vs] / dzv[v];
      }

      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* edges. This is the quantity used in Eq. A.39, Ringler et al */
      /* (2013).                                                     */
      memset(windat->nrvore, 0, window->szeS * sizeof(double));
      memset(windat->npvore, 0, window->szeS * sizeof(double));
      for (ee = 1; ee <= window->n2_e1; ee++) {
	e = window->w2_e1[ee];
	v1 = window->e2v[e][0];
	v2 = window->e2v[e][1];
	windat->nrvore[e] = 0.5 * (windat->nrvor[v1] + windat->nrvor[v2]);
	windat->npvore[e] = 0.5 * (windat->npvor[v1] + windat->npvor[v2]);
      }

      /*-------------------------------------------------------------*/
      /* Compute the normalized relative and planetary vorticity at  */
      /* cell centres.                                               */
      memset(windat->nrvorc, 0, window->szcS * sizeof(double));
      if (wincon->vorticity & POTENTIAL)
	memset(windat->pv, 0, window->szcS * sizeof(double));
      if (wincon->vorticity & ABSOLUTE)
	memset(windat->av, 0, window->szcS * sizeof(double));
      for(cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	cs = window->m2d[c];
	iarea = 1.0 / window->cellarea[cs];
	for (n = 1; n <= window->npe[cs]; n++) {
	  j = window->vIc[n][cs];
	  v = window->c2v[n][c];
	  vs = window->m2dv[v];
	  windat->nrvorc[c] += window->dualareap[vs][j] * windat->nrvor[v] * iarea;
	  if (wincon->vorticity & POTENTIAL)
	    windat->pv[c] += window->dualareap[vs][j] * iarea *
	      (windat->nrvor[v] + windat->npvor[v]);
	  if (wincon->vorticity & ABSOLUTE) {
	    windat->av[c] += window->dualareap[vs][j] * iarea *
	      (windat->rvor[v] + windat->fv[v]);
	  }
	}
      }

      /* Save the relative vorticity                                   */
      if (wincon->vorticity & (RELATIVE | TENDENCY)) {
	/* Save the relative vorticity at the previous timestep        */
	if (wincon->vorticity & TENDENCY)
	  memcpy(prv, windat->rv, window->szcS * sizeof(double));
	memcpy(windat->rv, windat->nrvorc, window->szcS * sizeof(double));
      }

      /*-------------------------------------------------------------*/
      /* Compute the kinetic energy at cell centres. Used in the     */
      /* decomposition of momentum advection, e.g. Eq. 4, Ringler et */
      /* al, 2010.                                                   */
      memset(windat->kec, 0, window->szcS * sizeof(double));
      memset(windat->div, 0, window->szcS * sizeof(double));
      for(cc = 1; cc <= window->a2_t; cc++) {
	c = window->w2_t[cc];
	cs = window->m2d[c];
	d1 = 1.0 / window->cellarea[cs];
	for (ee = 1; ee <= window->npe[cs]; ee++) {
	  e = window->c2e[ee][c];
	  es = window->m2de[e];
	  d2 = window->h1au1[es] * windat->u1av[e] * d1;
	  windat->div[c] += window->eSc[ee][cs] * d2;
	  windat->kec[c] += 0.25 * d2 * window->h2au1[es] * windat->u1av[e];
	}
      }
      /*
      for (cc = 1; cc <= window->nbptS; cc++) {
	c = window->bpt[cc];
	c1 = window->bin[cc];
	windat->kec[c] = windat->kec[c1];
      }
      */

      /*-------------------------------------------------------------*/
      /* Compute the nonlinear Coriolis tendency. The tangential     */
      /* velocity is computed using the primal to dual mesh mapping  */
      /* (Eq. 24, Ringler et al, 2010) and a mean of absolute        */
      /* vorticity at the edge in question. The gradient of kinetic  */
      /* energey is also added (see Eq. A.39, Ringler et al, 2013).  */
      if (wincon->u1av_f & ADVECT) {
	memset(windat->nrvore, 0, window->szeS * sizeof(double));
	memset(windat->kec, 0, window->szcS * sizeof(double));
      }
      if (wincon->u1av_f & CORIOLIS)
	memset(windat->npvore, 0, window->szeS * sizeof(double));
      for (ee = 1; ee <= window->v2_e1; ee++) {
	e = window->w2_e1[ee];
	es = window->m2de[e];
	c1 = window->e2c[e][0];
	c2 = window->e2c[e][1];
	
	d1 = 1.0 / window->h2au1[es];
	d2 = 0.0; 
	for (n = 1; n <= window->nee[es]; n++) {
	  eoe = window->eSe[n][e];
	  if (!eoe) continue;
	  /*d3 = 0.5 * (windat->nrvore[e] + windat->npvore[e] + windat->nrvore[eoe] + windat->npvore[eoe]);*/
	  d3 = wincon->pv_calc(windat, e, eoe, window->e2v[e][0], window->e2v[e][1]);
	  d2 += window->wAe[n][e] * windat->u1av[eoe] * d3 * windat->depth_e1[eoe];
	}
	d3 = d2 - (windat->kec[c1] - windat->kec[c2]) * d1;
	wincon->u1adv[e] -= (dtm * d3);
	windat->nu1av[e] += (dtu * d3);
      }
      if (wincon->u1_f & ADVECT)
	memset(wincon->u1adv, 0, window->szeS * sizeof(double));      
    }  /* End if(dtu > 0.0)                                          */

    /* Get the timestep for the next iteration.                      */
    /* Decrement the time remaining. This must be computed for       */
    /* both the leapfrog step (trem) used to update velocity, amd    */
    /* the forward step (tremf) used to integrate u1inter[] over     */
    /* forward time-step windat->dtf.                                */
    trem -= dtu;
    tremf -= dtm;
    dtu = dtm / lr;
    if (trem < dtu)
      dtu = trem;
    if (tremf < dtm)
      dtm = tremf;
    if(!wincon->sigma)
      vel = windat->nu1av;
  }  /* End while (trem > 0) */

  return(0);
}

/* END nonlin_coriolis_2d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/* Old code - no longer used for unstructured grids since we use     */
/* nonlin_coriolis_2d() instead. We retain this code for a possible  */
/* future tensor approach.                                           */
/*-------------------------------------------------------------------*/
void advect_u1_2d(geometry_t *window, /* Processing window           */
                  window_t *windat,   /* Window data structure       */
                  win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int e, ee, es, ep, em;    /* Edge coordinate / counter             */
  int e1, e2, e1s, e2s;     /* Edge coordinate / counter             */
  int c, cc;                /* Cell coordinate / counter             */
  int c1, c2, c1s, c2s;     /* Cell coordinate / counter             */
  int cs, cb, ems;          /* Surface / bottom sparse coordinate    */
  int j, npe;
  int *ctp;                 /* Cells to process vector               */
  int vcs;                  /* Number of cells to process            */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double u1val;             /* Cell centered u1 value                */
  double u2val;             /* Cell cornered u2 value                */
  double h1;                /* Face centered grid spacing            */
  double h1av, h2av;        /* Mean cell spacings x and y directions */
  double m, d1;             /* Dummy variables                       */
  double *Du1;              /* Dummy for e1 advective fluxes         */
  double *vel;              /* Velocity to use in spatial gradients  */
  double *area;             /* Cell area                             */
  double *u1sh1h2_av;       /* u1 * u1 * h1 * h2                     */
  double *u1u2h1s_av;       /* u1 * u2 * h1 * h1                     */
  double *depth;            /* Depth of the water column             */
  double *dvel;             /* Velocity to interpolate               */
  double *vface;            /* Velocity centered on the grid face    */
  double *cn;               /* Courant number                        */
  double midx;              /* Depth at the e1 face                  */
  double *u2au1;            /* u2 value at the e1 face               */
  double lr;                /* Ratio of dtf : dt2d                   */

  if (!wincon->nonlinear)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  Du1 = wincon->d1;
  u1sh1h2_av = wincon->d2;
  u1u2h1s_av = wincon->d3;
  area = window->cellarea;
  depth = windat->depth_e1;
  vel = windat->u1av;
  u2au1 = wincon->w7;
  if(wincon->dolin_u1) {
    ctp = wincon->s4;
    vcs = wincon->aclS;
  } else {
    ctp = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  if(windat->dtu1 == windat->dtf)
    dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping  */
  else
    dtu *= 2.0;           /* 2D leapfrog time-step with sub-stepping */
  lr = windat->dtf2 / windat->dt2d;
  dtm = dtu * lr;
  trem = windat->dt2d;
  tremf = windat->dtf2;
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Initialise                                                  */
      /* memset(depth,0,window->sgsizS*sizeof(double)); */
      memset(Du1, 0, window->szeS * sizeof(double));
      memset(u1sh1h2_av, 0, window->szeS * sizeof(double));
      memset(u1u2h1s_av, 0, window->szeS * sizeof(double));

      /*-------------------------------------------------------------*/
      /* Get the tracer values in multi-dt auxiliary cells. This is  */
      /* a linear interpolation between the value at the start of    */
      /* the timestep and the end of longer timesteps.               */
      set_multidt_tr_2d(window, windat, trem, vel, windat->u1av_as,
                        windat->u1av_ae);

      /*-------------------------------------------------------------*/
      /* Loop over the u1 points to calculate depths. Note: Du1 at   */
      /* the coordinate ym1 is required for the cell corner fluxes,  */
      /* and this cell is not included in the e1 cells to process    */
      /* vector, hence Du1 must be calculated over auxiliary cells   */
      /* on all sides in order to accomodate this.                   */
      for (ee = 1; ee <= window->n2_e1; ee++) {
        e = window->w2_e1[ee];
        Du1[e] = depth[e] * vel[e];
      }
      if (wincon->momsc2d & VANLEER) {
	for (ee = window->b2_e1 + 1; ee <= window->n2_e1; ee++) {
	  ee = window->w2_e1[ee];
	  em = window->em[e];
	  Du1[em] = depth[em] * vel[em];
	}
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres. This must be        */
      /* calculated at ghost cells also since these locations will   */
      /* give a non-zero cell centered flux.                         */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  cs = window->e2c[e][0];
	  ep = window->ep[e];

          u1val = 0.5 * (vel[e] + vel[ep]);
          /* SIGMA : Multiply by the cell centre depth               */
          if (u1val > 0)
            m = Du1[e] * wincon->mdx[e];
          else
            m = Du1[ep] * wincon->mdx[ep];
          u1sh1h2_av[e] = m * u1val * area[cs];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  cs = window->e2c[e][0];
	  ep = window->ep[e];

          u1val = 0.5 * (vel[e] + vel[ep]);
          /* SIGMA : Multiply by the cell centre depth               */
          u1sh1h2_av[e] =
            u1val * area[cs] * 0.5 * (Du1[e] * wincon->mdx[e] +
				      Du1[ep] * wincon->mdx[window->m2de[ep]]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  cs = window->e2c[e][0];
	  ep = window->ep[e];
          vface[e] = 0.5 * (vel[e] + vel[ep]);
          cn[c] = vface[c] * dtu / window->h1acell[c];
          /* The interpolated velocity must be shifted by one cell   */
          /* when using van_leer_do() since this routine is set up   */
          /* to interpolate onto the cell face and we need to        */
          /* interpolate to the cell center.                         */
          dvel[e] = Du1[ep] * wincon->mdx[ep];
        }
	/* Ghost cells 2 cells in */
        for (ee = window->b2_e1 + 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  em = window->em[e];
          dvel[em] = Du1[e] * wincon->mdx[e];
	}
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array u2h1h2.                                       */
	for (ee = 1; ee <= window->n2_e1; ee++) {
	  int cp, cm;
	  e = window->w2_e1[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
	  j = window->e2e[e][0];
	  cm = window->c2c[j][c2];
	  j = window->e2e[e][1];
	  cp = window->c2c[j][c1];
	  van_leer_do(u1sh1h2_av, dvel, vface, cn, e, c1, cp, c2, cm);
	}

        /* Multiply by vertical velocity and shift one cell down     */
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  cs = window->e2c[e][0];
          u1sh1h2_av[e] *= vface[e] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners. This must be        */
      /* calculated at ghost cells also since these locations will   */
      /* give a non-zero flux at ouside corners.                     */
      /* First order upstream scheme                                 */
      if (wincon->momsc2d & ORDER1) {
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w3_e1[ee];
	  c = window->m2d[window->e2c[e][0]];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ete(window, e, &e1, &e2);
          h1 = 0.5 * (window->h1au1[e1] + window->h1au1[e2]);
          u2val = 0.5 * (windat->u1av[e1] * wincon->mdy[e1] +
                         windat->u1av[e2] * wincon->mdy[e2]);
          if (u2val > 0)
            m = Du1[em];
          else
            m = Du1[e];
          u1u2h1s_av[e] = m * h1 * h1 * u2val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc2d & ORDER2) {
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w3_e1[ee];
	  c = window->m2d[window->e2c[e][0]];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ete(window, e, &e1, &e2);

          h1 = 0.5 * (window->h1au1[e1] + window->h1au1[e2]);
          u2val = 0.5 * (windat->u1av[e1] * wincon->mdy[e1] +
                         windat->u1av[e2] * wincon->mdy[e2]);
          h1av = window->h1au1[e2] / (window->h1au1[e1] + window->h1au1[e2]);
          h2av = window->h2au1[em] / (window->h2au1[em] + window->h2au1[e]);
          /* SIGMA : Multiply by the cell corner depth               */

          u2val = h1av * (windat->u1av[e1] * wincon->mdx[e1] -
                          windat->u1av[e2] * wincon->mdx[e2]) +
            windat->u1av[e2] * wincon->mdx[e2];
          m = h2av * (Du1[e] - Du1[em]) + Du1[em];

          u1u2h1s_av[e] = m * h1 * h1 * u2val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc2d & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        dvel = wincon->w4;
        memset(dvel, 0, window->sgsizS * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w3_e1[ee];
	  c = window->m2d[window->e2c[e][0]];
	  npe = window->npe[c];
	  j = window->e2e[e][0];
	  j = (j <= npe/2) ? jm(j, npe) : jp(j, npe);
	  em = e2e(window, e, j);
	  ete(window, e, &e1, &e2);
          vface[e] = 0.5 * (windat->u1av[e1] + windat->u1av[e2]);
          cn[e] =
            2.0 * vface[e] * dtu / (window->h1au1[e1] + window->h1au1[e2]);
          dvel[e] = Du1[e] * wincon->mdx[e];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array dvel.                                         */
	for (ee = 1; ee <= window->n2_e1; ee++) {
	  int cp, cm;
	  e = window->w2_e1[ee];
	  c1 = window->e2c[e][0];
	  c2 = window->e2c[e][1];
	  j = window->e2e[e][0];
	  cm = window->c2c[j][c2];
	  j = window->e2e[e][1];
	  cp = window->c2c[j][c1];
	  van_leer_do(u1u2h1s_av, dvel, vface, cn, e, c1, cp, c2, cm);
	}

        /* Multiply by vertical velocity and shift one cell down     */
        for (ee = 1; ee <= window->n2_e1; ee++) {
          e = window->w2_e1[ee];
	  cs = window->e2c[e][0];
          u1u2h1s_av[e] *= vface[e] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the u1 updated solution                                 */
      for (ee = 1; ee <= vcs; ee++) {
        e = ctp[ee];
	j = window->e2e[e][0];
	em = e2e(window, e, j);
	j = window->e2e[e][0];
	j = (j <= npe/2) ? jp(j, npe) : jm(j, npe);
	ep = e2e(window, e, j);

        /* SIGMA : Multiply by cell depth                            */
        d1 = wincon->u1c1[e] * (u1sh1h2_av[e] - u1sh1h2_av[em] + 
				u1u2h1s_av[ep] - u1u2h1s_av[e]) / 
	  max(depth[e], wincon->hmin) +
          wincon->u1c3[e] * vel[e] * vel[e] * wincon->mdx[e] +
          wincon->u1c4[e] * u2au1[e] * u2au1[e] * wincon->mdx[e];

        windat->nu1av[e] += d1 * dtu;

        /* Integrate the dispersion terms over the sub-timestep.     */
        /* u1adv is divided by dt to get the mean over the 2D mode   */
        /* in vel_u1_update().                                       */
        wincon->u1adv[e] -= (d1 * dtm);
      }
    }

    /* Get the timestep for the next iteration.                      */
    /* Decrement the time remaining. This must be computed for       */
    /* both the leapfrog step (trem) used to update velocity, amd    */
    /* the forward step (tremf) used to integrate u1inter[] over     */
    /* forward time-step windat->dtf.                                */
    trem -= dtu;
    tremf -= dtm;
    dtu = dtm / lr;
    if (trem < dtu)
      dtu = trem;
    if (tremf < dtm)
      dtm = tremf;
    if(!wincon->sigma)
      vel = windat->nu1av;
  }
  debug_c(window, D_UA, D_ADVECT);
}

/* END advect_u1_2d()                                          (slf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implement the porus plate algorithm in the e1 direction           */
/*-------------------------------------------------------------------*/
double porus_plate_e1(geometry_t *window,    /* Processing window    */
		      window_t *windat,      /* Window data          */
		      win_priv_t *wincon,    /* Window constants     */
		      int e
		      )
{
  int c1,c2,e2;
  double Cdrag = 1.5;        /* Drag coefficient for coral reefs     */
  double eff = 0.35;         /* Effective wet cross sectional area   */
  double closs;              /* Energy loss coefficient              */
  double d1, reef;

  e2 = window->m2de[e];
  c1 = window->e2c[e][0];
  c2 = window->e2c[e][1];
  c1 = (window->wgst[c1]) ? window->wgst[c1] : c1;
  c2 = (window->wgst[c2]) ? window->wgst[c2] : c2;
  reef = (windat->reefe1[c1] > 0.0) ? windat->reefe1[c1] : 0.0;
  if (windat->reefe1[c2] > 0.0) reef = 0.5 * (reef + windat->reefe1[c2]);

  if (reef) {
    /* Get the energy loss coefficient                             */
    d1 = 1.0 / eff;
    closs = Cdrag * reef * d1 * d1 / 2.0;
    /* Get the reef friction term                                  */
    d1 = windat->u1[e] * windat->u1[e] + windat->u2[e] * windat->u2[e];
    return(-closs * windat->u1[e] * sqrt(d1) / window->h2au1[e2]);
  }
  return(0.0);
}

/* END porus_plate_e1()                                              */
/*-------------------------------------------------------------------*/




