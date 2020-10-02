/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/momentum/advect.c
 *  
 *  Description: Momentum advection schemes for 2D and 3D modes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: advect.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void set_dzzv(geometry_t *window, double *dzz, int *cells, int vc, int vcs, int mode);
void vel_face(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	      double *nu, double *nv, double *nw, 
	      int *cells, int vcs, int vc, int mode);
int get_posv(geometry_t *window, int c, double nu, double nv, double nw,
	     double *cx, double *cy, double *cz, double dt, int *nc, 
	     int dir, double *pp, double *qq);
void set_eov(geometry_t *window, win_priv_t *wincon, double h1, double h2, double *dzz,
	     int *c, double *cx, double *cy, double *cz);
void weights_oov(geometry_t *window, int c, int co, double pin, double qin, double rin,
		 double *h1, double *h2, double *dzz);
void weights_eov(geometry_t *window, int c, int co, double pin, double qin, double rin,
		 double *h1, double *h2, double *dzz);
int get_posvS(geometry_t *window, int c, double nu, double nv, 
	      double *cx, double *cy, double dt, int dir, double *p, double *q);
void set_eovS(geometry_t *window, win_priv_t *wincon, double h1, double h2,
	      int *c, double *cx, double *cy);
void weights_oovS(geometry_t *window, int c, int co, double pin, double qin, 
		  double *h1, double *h2);
void weights_eovS(geometry_t *window, int c, int co, double pin, double qin,
		  double *h1, double *h2);
double geth(double *h, double f, int c, int *m1);
double porus_plate_e1(geometry_t *window, window_t *windat, win_priv_t *wincon, int c);
double porus_plate_e2(geometry_t *window, window_t *windat, win_priv_t *wincon, int c);
int cg;
#define SMALL 1e-10             /* Minimum value for velocity bounds */

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U1 ADVECTION ROUTINES                                             */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/*-------------------------------------------------------------------*/
int advect_u1_3d(geometry_t *window, /* Processing window           */
		 window_t *windat,   /* Window data                 */
		 win_priv_t *wincon  /* Window constants            */
  )
{
  int c, cc;                /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int xmyp1;                /* Sparse coordinate at (i-1,j+1)        */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int xmzp1;                /* Sparse coordinate at (i-1,k-1)        */
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
    for (cc = 1; cc <= window->n3_e1; cc++) {
      c = window->w3_e1[cc];
      dzu1[c] = 1.0;
    }
  } else
    dzu1 = windat->dzu1;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required.                             */
  /* The advective terms can create large velocities due to momentum */
  /* divergence if the sea level lies just above a layer (i.e. a     */
  /* layer) in a cell and just below the layer in the adjacent cell. */
  /* Here dzu1 is thin at cellc, ush1h2 may be large at c and small  */
  /* at xm1, hence (ush1h2[c]-ush1h2[xm1])/dzu1[c] can be large.     */
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
    for (cc = 1; cc <= vc; cc++) {
      c = cells[cc];            /* Wet cell to process               */
      cs = window->m2d[c];      /* 2D cell corresponding to 3D cell  */
      xm1 = window->xm1[c];
      yp1 = window->yp1[c];
      /* xmyp1=window->yp1[xm1]; */
      xmyp1 = window->xmyp1[c];
      zp1 = window->zp1[c];
      xmzp1 = window->zp1[xm1];

      /* Minimum time-step due to x velocity                         */
      vlc = windat->u1[c];
      if (cc <= vcs) wq[cs] = dtu * vlc / window->h1au1[cs];
      if (fabs(vlc) < minval)
        vlc = minval;
      d1 = sf * fabs(window->h1au1[cs] / vlc);
      if (d1 < dtm) {
        dtm = d1;
        tp = 'u';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
      /* Minimum time-step due to y velocity                         */
      vlc = u2au1[c];
      if (cc <= vcs) wq[cs] = max(wq[cs], dtu * vlc / window->h2au1[cs]);
      if (fabs(vlc) < minval)
        vlc = minval;
      d2 = sf * fabs(window->h2au1[cs] / vlc);
      if (d2 < dtm) {
        dtm = d2;
        tp = 'v';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
      /* Minimum time-step due to z velocity                         */
      vlc =
        (cc <=
         vcs) ? windat->wtop[cs] +
        windat->wtop[window->xm1[cs]] : windat->w[zp1] + windat->w[xmzp1];
      vlc = 0.25 * (vlc + windat->w[c] + windat->w[xm1]);
      if (cc <= vcs) wq[cs] = max(wq[cs], dtu * vlc / (windat->dzu1[c] * wincon->mdx[cs]));
      if (wincon->stab & SUB_STEP_NOSURF && cc <= vcs)
        vlc = 0.0;
      if (fabs(vlc) < minval)
        vlc = minval;
      /* SIGMA : Multiply by depth */
      d3 =
        wsf * fabs(max(windat->dzu1[c] * wincon->mdx[cs], wincon->hmin) /
                   vlc);
      if (d3 < dtm) {
        dtm = d3;
        tp = 'w';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
    }
  }
  if (dtm != dtu) {
    c = ceil((windat->dtb + dtu) / (2 * dtm));
    dtm = (windat->dtb + dtu) / (2.0 * (double)c);
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

  if (wincon->momsc & LAGRANGE) {
    semi_lagrange_u1(window, windat, wincon, cells, vcs, vc);
    if (wincon->nbl1) {
      vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
		window->bpte1, window->bine1, wincon->slip);
      blend_vel(window, windat, wincon, U13D, vel);
    }
    return(0);
  }

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
  memcpy(wincon->w4, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e1; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_e1[cc]; /* 2D coordinate                     */
    cb = window->bot_e1[cc];    /* Bottom e1 coordinate              */

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
  for (cc = window->a2_e1 + 1; cc <= window->n2_e1; cc++) {
    cs = c = window->w2_e1[cc];
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
      memset(ush1h2, 0, window->sgsiz * sizeof(double));
      memset(u1u2hs, 0, window->sgsiz * sizeof(double));
      memset(wu1, 0, window->sgsiz * sizeof(double));
      if (wincon->momsc & ZERO_DRYK)
	set_surf_cond(window, vel, 0.0, window->sur_e1, window->b2_e1);

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres (u1*u1*h1*h2*dz).    */
      /* This must be calculated at ghost cells also since these     */
      /* locations will give a non-zero cell centered flux.          */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xp1 = window->xp1[c];

          u1val = 0.5 * (vel[xp1] + vel[c]);
          /* SIGMA : Multiply by the cell centre depth               */
          if (u1val > 0)
            m = dzu1[c] * vel[c] * wincon->mdx[cs];
          else
            m =
              dzu1[xp1] * vel[xp1] * wincon->mdx[window->m2d[xp1]];
          ush1h2[c] = m * u1val * area[cs];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xp1 = window->xp1[c];
          u1val = 0.5 * (vel[xp1] + vel[c]);
          /* SIGMA : Multiply by the cell centre depth               */
          m =
            0.5 * (dzu1[xp1] * vel[xp1] *
                   wincon->mdx[window->m2d[xp1]] +
                   dzu1[c] * vel[c] * wincon->mdx[cs]);
          ush1h2[c] = m * u1val * area[cs];
	}
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xp1 = window->xp1[c];
          vface[c] = 0.5 * (vel[xp1] + vel[c]);
          cn[c] = vface[c] * dtu / window->h1acell[cs];
          /* The interpolated velocity must be shifted by xp1 when   */
          /* using van_leer_do() since this routine is set up to     */
          /* interpolate onto the cell face and we need to           */
          /* interpolate to the cell center. Use wu1 as a dummy.     */
          wu1[c] =
            dzu1[xp1] * vel[xp1] * wincon->mdx[window->m2d[xp1]];
        }
	/* Ghost cells 2 cells in */
        for (cc = window->b3_e1 + 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xm1 = window->xm1[c];
          wu1[xm1] = dzu1[c] * vel[c] *
            wincon->mdx[window->m2d[c]];
	}

        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
        van_leer_do(ush1h2, wu1, vface, cn, 1, window->n3_e1,
                    window->w3_e1, window->xp1, window->xm1);

        /* Multiply by centered velocity and shift one cell down     */
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xp1 = window->xp1[c];
          ush1h2[c] *= vface[c] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners (u1*u2*h1*h1*dz).    */
      /* This must be calculated at ghost cells also since these     */
      /* locations will give a non-zero flux at ouside corners.      */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          cs = window->m2d[c];
          h1 = 0.5 * (window->h1au2[window->xm1[cs]] + window->h1au2[cs]);
          u2val = 0.5 * (windat->u2[xm1] * wincon->mdy[window->m2d[xm1]] +
                         windat->u2[c] * wincon->mdy[cs]);
          if (u2val > 0)
            m = dzu1[ym1] * vel[ym1];
          else
            m = dzu1[c] * vel[c];
          u1u2hs[c] = m * h1 * h1 * u2val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          cs = window->m2d[c];

          xm1s = window->xm1[cs];
          ym1s = window->ym1[cs];
          h1 = 0.5 * (window->h1au2[xm1s] + window->h1au2[cs]);
          h1av =
            window->h1au2[xm1s] / (window->h1au2[xm1s] +
                                   window->h1au2[cs]);
          h2av =
            window->h2au1[ym1s] / (window->h2au1[ym1s] +
                                   window->h2au1[cs]);
          /* SIGMA : Multiply by the cell corner depth               */
          u2val = h1av * (windat->u2[c] * wincon->Ds[cs] -
                          windat->u2[xm1] * wincon->Ds[xm1s]) +
            windat->u2[xm1] * wincon->Ds[xm1s];
          m1 = dzu1[ym1] * vel[ym1];
          m2 = dzu1[c] * vel[c];
          m = h2av * (m2 - m1) + m1;
          u1u2hs[c] = m * h1 * h1 * u2val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sgsiz * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          xm1 = window->xm1[c];
          vface[c] = 0.5 * (windat->u2[xm1] + windat->u2[c]);
          cn[c] = 2.0 * vface[c] * dtu /
            (window->h2au2[cs] + window->h2au2[window->xm1[cs]]);
          wu1[c] = dzu1[c] * vel[c] * wincon->mdx[cs];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
        van_leer_do(u1u2hs, wu1, vface, cn, 1, window->n3_e1,
                    window->w3_e1, window->yp1, window->ym1);

        /* Multiply by velocity and shift one cell down              */
        for (cc = 1; cc <= window->n3_e1; cc++) {
          c = window->w3_e1[cc];
          cs = window->m2d[c];
          u1u2hs[c] *= vface[c] * area[cs];
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
      memcpy(velbuf, vel, window->sgsiz * sizeof(double));
      for (cc = 1; cc <= window->a2_e1; cc++) {
        c = window->bot_e1[cc]; /* 3D bottom coordinate              */
        zm1 = window->zm1[c];   /* Sediment coordinate               */
        velbuf[zm1] = vel[c];   /* u1 : no-gradient across bottom    */
      }

      /* Surface vertical advective flux                             */
      for (cc = 1; cc <= wincon->vcs; cc++) {
        c = cells[cc];
        cs = window->m2d[c];
        xm1 = window->xm1[cs];
        wu1[c] = 0.5 * vel[c] * (windat->wtop[xm1] + windat->wtop[cs]);
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
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          xm1 = window->xm1[c];
          zm1 = window->zm1[c];
          wval = 0.5 * (wvel[xm1] + wvel[c]);
          if (wval > 0)
            wu1[zm1] = velbuf[zm1] * wval;
          else
            wu1[zm1] = velbuf[c] * wval;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          xm1 = window->xm1[c];
          zm1 = window->zm1[c];
          wval = 0.5 * (wvel[xm1] + wvel[c]);
	  dzav = windat->dzu1[zm1] / (windat->dzu1[zm1] + windat->dzu1[c]);
	  wu1[zm1] = wval * (dzav * (velbuf[c] - velbuf[zm1]) + velbuf[zm1]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sgsiz * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_t; cc++) {
          c = window->w3_t[cc];
          cs = window->m2d[c];
          xm1 = window->xm1[c];
          dzav = 0.5 * (windat->dzu1[c] + windat->dzu1[window->zm1[c]]);
          vface[c] = 0.5 * (wvel[xm1] + wvel[c]);
          cn[c] = vface[c] * dtu / (dzav * wincon->Ds[cs]);
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wvel.                                         */
        van_leer_do(wvel, velbuf, vface, cn, 1, vc, cells,
                    window->zp1, window->zm1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          zm1 = window->zm1[c];
          wu1[zm1] = vface[c] * wvel[c];
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
      /* mean in the horizontal, 0.5(wbot[c]+w[xm1]) or              */
      /* 0.5(w[c]+wbot[xm1]).                                        */
      for (cc = 1; cc <= window->v2_e1; cc++) {
        c = window->bot_e1[cc]; /* Bottom coordinate                 */
        cs = window->w2_e1[cc];
        xm1 = window->xm1[cs];
        zm1 = window->zm1[c];
        wu1[zm1] =
          0.5 * velbuf[c] * (windat->wbot[cs] + windat->wbot[xm1]);
      }

      /*-------------------------------------------------------------*/
      /* Get the momentum divergence and metric terms.               */
      for (cc = 1; cc <= vc; cc++) {
        c = cells[cc];
        cs = window->m2d[c];
        xm1 = window->xm1[c];
        yp1 = window->yp1[c];
        zm1 = window->zm1[c];
        midx = wincon->mdx[cs];

        /* SIGMA : Multiply by cell depth                            */
        d1 =  wincon->u1c1[cs] * (ush1h2[c] - ush1h2[xm1] + 
				  u1u2hs[yp1] - u1u2hs[c]) /
	  max(dzu1[c], wincon->hmin) +
          wincon->u1c3[cs] * vel[c] * vel[c] * midx +
          wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx;

	if (wincon->porusplate)
	  d1 += porus_plate_e1(window, windat, wincon, c);

        /* Integrate the e1 dispersion term vertically and with time */
        /* (i.e. over the sub-timestepping; windat->dtf).            */
        wincon->u1inter[cs] += (d1 * dtm * windat->dzu1[c]);

        /* Add the vertical advective flux. Note : the difference of */
        /* this flux is wu1[c]-wu1[zm1] due to the stagger imposed,  */
        /* not wu1[zp1]-wu1[c].                                      */
	if (!(wincon->momsc & WIMPLICIT))
	  d1 -= (wu1[c] - wu1[zm1]) / max(windat->dzu1[c], wincon->hmin);

        /* Add to new u1 value */
	div[c] = d1 * dtu;
        /*windat->nu1[c] += d1 * dtu;*/
      }

      /* Filter the divergence if required                           */
      if (wincon->filter & ADVECT) {
	vel2D_lbc(div, window->nbpte1, window->nbe1,
		  window->bpte1, window->bine1, wincon->slip);
	shapiro(window, div, velbuf, cells, vc, 1, XDIR);
	shapiro(window, div, velbuf, cells, vc, 1, YDIR);
      }

      /* No advection at the blend-zone boundaries for stability     */
      for (ii = 0; ii < wincon->nbl2; ii++) {
	blend_t *blend = wincon->ble2[ii];
	for (cc = 1; cc <= blend->nlim; cc++) {
	  c = blend->lim[cc];
	  div[c] = 0.0;
	}
	for (cc = 1; cc <= blend->nlimS; cc++) {
	  c = blend->lim[cc];
	  wincon->u1inter[c] = 0.0;
	}
      }

      /* Get the u1 updated solution                                 */
      for (cc = 1; cc <= vc; cc++) {
        c = cells[cc];
        windat->nu1[c] += div[c];
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

      if (wincon->nbl1) {
	vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
		  window->bpte1, window->bine1, wincon->slip);
	blend_vel(window, windat, wincon, U13D, vel);
      }

    }
  }
  return(0);
}

/* END advect_u1_3d()                                          (slf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* backward direction. Formulated using the flux form.               */
/*-------------------------------------------------------------------*/
void advect_u1_3d_ang_adv_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int xmyp, xmzm;
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h1au1;            /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midx;             /* Depth at the e1 face                  */
  double dz;                /* Cell thickness at cell face           */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int *ctp;                 /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w5;
  u2au1 = wincon->w7;
  h1au1 = window->h1au1;
  vel = windat->u1;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu1;
  cells = wincon->s4;
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required                              */
  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    ctp = wincon->s3;
    vc = wincon->ncl;
  } else {
    ctp = wincon->s1;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  reorder_cells(window, cells, ctp, vc, vcs, wincon->i5, 1);
  dtu = windat->dtf + windat->dtb;
  windat->dtu1 = windat->dt;

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
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_t[cc];  /* 2D coordinate                     */
    cb = window->bot_t[cc];     /* Bottom e1 coordinate              */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];         /* w : top set to wtop             */
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xm1s = window->m2d[xm1];
    umid[c] = 0.25 * (2.0 * vel[c] + vel[xp1] + vel[xm1]) *
              window->h2au1[cs];
    vmid[c] = windat->u2[c] * window->h1au2[cs];
  }
  /* Stagger wmid so that wtop occupies the surface cell             */
  for(cc = 1; cc <= vc; cc++) {
    c = ctp[cc];
    zm1 = window->zm1[c];
    wmid[zm1] = wvel[c];
  }
  /* Surface and bottom conditions for wmid                          */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cb = window->zm1[wincon->i5[cc]];
    cs = window->m2d[c];
    wmid[c] = wtop[cs];
    wmid[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the cells above the surface                                 */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    zp1 = window->zp1[c];
    bmask[c] = 4.0;
    while (c != zp1) {
      bmask[zp1] = 2.0;
      c = zp1;
      zp1 = window->zp1[c];
    }
  }
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1 + 1; cc <=  window->nbpte1; cc++) {
    c = window->bpte1[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
	      open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
	      open[n]->cyc_e1, windat->nu1, windat->u1, 
	      windat->u1b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
	      windat->u1, windat->u1b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u1;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration : nu1[zm1] is known.                         */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = vc; cc >= 1; cc--) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    xmyp = window->yp1[xm1];
    xmzm = window->zm1[xm1];

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xp1] == 1) {
      ush1h2 = 0.5 * wincon->u1c1[cs] * umid[c] * 
	(vel[c] * h1au1[cs] - vel[xm1] * h1au1[xm1s]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h1au1[cs];
      d2 = h1au1[xp1s];
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[cs] * umid[c] * (d1 * vel[c] - d2 * vel[xp1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[cs] * umid[c] * 
	(vel[c] * h1au1[cs] + nvel[xp1] * h1au1[xp1s] -
	 vel[xm1] * h1au1[xm1s]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c] * h1au1[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[yp1] = nu1[c] due to the no-slip condition.    */
    if(bmask[yp1] == 1) {
      u1u2hs = 0.25 * wincon->u1c1[cs] * 
	(vmid[c] + vmid[xm1]) * (vel[c] * h1au1[cs] - vel[ym1] * h1au1[ym1s]);
				 
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u1c1[cs] * 
	((vmid[yp1] + vmid[xmyp]) * nvel[yp1] * h1au1[yp1s] +
	 (vmid[c] + vmid[xm1]) * (vel[c] * h1au1[cs] - vel[ym1] * h1au1[ym1s]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u1c1[cs] * (vmid[yp1] + vmid[xmyp]) * h1au1[cs];
      exe2 = ime2 * vel[c];
    }

    dz = max(windat->dzu1[c], wincon->hmin);
    /* Get the contribution from w.du/dz (note the wmid stagger)     */
    if(bmask[c] == 4) {
      d1 = wmid[c] + wmid[xm1] + wmid[zm1] + wmid[xmzm];
      wu1 = -0.125 * d1 * (vel[c] - (vel[c] + nvel[zm1])) / dz;

      /* Get the implicit multiplier for v.du/dy                       */
      imz = 0.125 * d1 /dz;
      exz = imz * vel[c];
    }
    else {
      wu1 = -0.25 * ((wmid[c] + wmid[xm1]) * (vel[zp1] - vel[c]) +
		     (wmid[zm1] + wmid[xmzm]) * nvel[zm1]) / dz;
    
      /* Get the implicit multiplier for v.du/dy                       */
      imz = 0.25 * (wmid[zm1] + wmid[xmzm]) / dz;
      exz = imz * vel[c];
    }
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2 + imz);
    adv = (ush1h2 + u1u2hs + wu1);

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv -= (exe1 + exe2 + exz);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[cs] * vel[c] * vel[c] * midx[cs] +
           wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx[cs];
    /* Note : nu1 is a copy of u1b, but may be adjusted for thin     */
    /* layers, hence use nu1 for u1b here.                           */
    velb = windat->nu1[c];
    windat->nu1[c] = (velb + dtu * adv ) / immp;
    
    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = (u1b-dtu.w.du/dz) - nu1       */
    if(wincon->momsc & EXPLICIT) {
      adv -= (wu1 - exz);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
    }
    else {
      immp = 1.0 + dtu * imz;
      adv = (velb + dtu * wu1) / immp;
      wincon->u1inter[cs] += ((adv - windat->nu1[c]) * windat->dzu1[c] *
			      windat->dtf / dtu);

      /* Fill cells above the free surface to the top of the grid    */
      if(bmask[zp1] == 2) {
	cs = c;
	while (c != zp1) {
	  windat->nu1[zp1] = windat->nu1[cs];
	  c = zp1;
	  zp1 = window->zp1[c];
	}
      }
    }
  }
}

/* END advect_u1_3d_ang_adv_b()                                (sbf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* forward direction. Formulated using the advective form.           */
/*-------------------------------------------------------------------*/
void advect_u1_3d_ang_adv_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int xmyp, xmzm;
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h1au1;            /* Cell width at e2 face                 */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midx;             /* Depth at the e1 face                  */
  double dz;                /* Cell thickness at cell face           */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u2au1 = wincon->w7;
  vel = wincon->w8;
  h1au1 = window->h1au1;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu1;
  memcpy(vel, windat->u1, window->sgsiz * sizeof(double));
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required                              */
  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  dtu = windat->dtf + windat->dtb;
  windat->dtu1 = windat->dt;

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
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e1; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_e1[cc]; /* 2D coordinate                     */
    cb = window->bot_e1[cc];    /* Bottom e1 coordinate              */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];         /* w : top set to wtop             */
      c = window->zm1[c];
    }
    /* Set the bottom boundary condition. Find the cell centered     */
    /* bottom coordinate and set w to wbot.                          */
    c = cb;                      /* e1 face centered bottom cell     */
    while (c != window->zm1[c])
      c = window->zm1[c];        /* Cell centered sediment location  */
    cb = window->zp1[c];         /* cell centered bottom location    */
    wvel[cb] = windat->wbot[cs]; /* w : bottom set to wbot           */
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    umid[c] = 0.25 * (2.0 * vel[c] + vel[xp1] + vel[xm1]) *
              window->h2au1[cs];
    vmid[c] = windat->u2[c] * window->h1au2[cs];
  }
  /* Stagger wmid so that wtop occupies the surface cell             */
  for(cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    zm1 = window->zm1[c];
    wmid[zm1] = wvel[c];
  }
  /* Surface and bottom conditions for wmid                          */
  for(cc = 1; cc <= vcs; cc++) {
    c = cells[cc];
    cb = window->zm1[wincon->i5[cc]];
    cs = window->m2d[c];
    wmid[c] = wtop[cs];
    wmid[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1 + 1; cc <=  window->nbpte1; cc++) {
    c = window->bpte1[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
	      open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
	      open[n]->cyc_e1, windat->nu1, windat->u1, 
	      windat->u1b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
	      windat->u1, windat->u1b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u1;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    xmyp = window->yp1[xm1];
    xmzm = window->zm1[xm1];

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xm1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xm1]) {
      ush1h2 = 0.5 * wincon->u1c1[cs] * umid[c] * 
	(vel[xp1] * h1au1[xp1s] - vel[c] * h1au1[cs]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h1au1[c];
      d2 = h1au1[xm1s];
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[cs] * umid[c] * (d1 * vel[c] - d2 * vel[xm1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[cs] * umid[c] * 
	(vel[xp1] * h1au1[xp1s] -
	 vel[c] * h1au1[cs] + nvel[xm1] * h1au1[xm1s]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c] * h1au1[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[ym1] = nu1[c] due to the no-slip condition.    */
    if(bmask[ym1]) {
      u1u2hs = 0.25 * wincon->u1c1[cs] * 
	(vmid[yp1] + vmid[xmyp]) *
	(vel[yp1] * h1au1[yp1s] - vel[c] * h1au1[cs]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u1c1[cs] * 
	((vmid[yp1] + vmid[xmyp]) *
	 (vel[yp1] * h1au1[yp1s] - vel[c] * h1au1[cs]) +
	 (vmid[c] + vmid[xm1]) * nvel[ym1] * h1au1[ym1s]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u1c1[cs] * (vmid[c] + vmid[xm1]) * h1au1[cs];
      exe2 = ime2 * vel[c];
    }

    /* Get the contribution from w.du/dz                             */
    dz = max(windat->dzu1[c], wincon->hmin);
    wu1 = -0.25 * ((wmid[c] + wmid[xm1]) * nvel[zp1] +
		   (wmid[zm1] + wmid[xmzm]) * (vel[c] - vel[zm1])) / dz;

    /* Get the implicit multiplier for w.du/dz                       */
    imz = 0.25 * (wmid[c] + wmid[xm1]) / dz;
    exz = imz * vel[c];
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz);
    adv = ush1h2 + u1u2hs + wu1;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += (exe1 + exe2 + exz);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[cs] * vel[c] * vel[c] * midx[cs] +
           wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx[cs];
    velb = windat->nu1[c];
    windat->nu1[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = (u1b-dtu.w.du/dz) - nu1       */
    if(wincon->momsc & EXPLICIT) {
      adv -= (exz + wu1);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
    }
    else {
      immp = 1.0 - dtu * (imz);
      adv = (velb + dtu * wu1) / immp;
      wincon->u1inter[cs] += ((adv - windat->nu1[c]) * windat->dzu1[c] *
			      windat->dtf / dtu);
    }
  }
}

/* END advect_u1_3d_ang_adv_f()                                (saf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* backward direction. Formulated using the flux form.               */
/*-------------------------------------------------------------------*/
void advect_u1_3d_ang_flux_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double fluxc;             /* Flux through cell face                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midx;             /* Depth at the e1 face                  */
  double *dz;               /* Cell thickness at cell face           */
  double dzd;               /* Cell thickness divisor                */
  double h2mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int *ctp;                 /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w5;
  u2au1 = wincon->w7;
  area = window->cellarea;
  vel = windat->u1;
  midx = wincon->mdx;
  dz = windat->dzu1;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu1;
  cells = wincon->s4;
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required                              */
  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    ctp = wincon->s3;
    vc = wincon->ncl;
    /*
    wtop = wincon->d3;
    reset_thin_wtop(window, windat, wincon, wincon->kth_e1, 
		    wincon->nkth_e1, wtop);
    */
  } else {
    ctp = wincon->s1;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  reorder_cells(window, cells, ctp, vc, vcs, wincon->i5, 1);
  dtu = windat->dtf + windat->dtb;
  windat->dtu1 = windat->dt;

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
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_t[cc];  /* 2D coordinate                     */
    cb = window->bot_t[cc];     /* Bottom e1 coordinate              */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];         /* w : top set to wtop             */
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xm1s = window->m2d[xm1];
    h2mid = 0.5 * (window->h1au2[cs] + window->h1au2[xm1s]);
    umid[c] = 0.5 * area[cs] * (vel[c] + vel[xp1]);
    vmid[c] = 0.5 * h2mid * h2mid * 
             (windat->u2[c] * Ds[cs] + windat->u2[xm1] * Ds[xm1s]);
  }
  /* Stagger wmid so that wtop occupies the surface cell             */
  for(cc = 1; cc <= vc; cc++) {
    c = ctp[cc];
    zm1 = window->zm1[c];
    xm1 = window->xm1[c];
    wmid[zm1] = 0.5 * (wvel[c] + wvel[xm1]);
  }
  /* Surface and bottom conditions for wmid                          */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cb = window->zm1[wincon->i5[cc]];
    cs = window->m2d[c];
    xm1 = window->xm1[cs];
    wmid[c] = 0.5 * (wtop[cs] + wtop[xm1]);
    wmid[cb] = 0.5 * (windat->wbot[cs] + windat->wbot[xm1]);
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the cells above the surface                                 */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    zp1 = window->zp1[c];
    bmask[c] = 4.0;
    while (c != zp1) {
      bmask[zp1] = 2.0;
      c = zp1;
      zp1 = window->zp1[c];
    }
  }
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1 + 1; cc <=  window->nbpte1; cc++) {
    c = window->bpte1[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
	      open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
	      open[n]->cyc_e1, windat->nu1, windat->u1, 
	      windat->u1b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
	      windat->u1, windat->u1b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u1;
  }
  if(wincon->momsc & ADVECT_FORM) {
    dz = wincon->w5;
    for(cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      dz[c] = 1.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration : nu1[zm1] is known.                         */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = vc; cc >= 1; cc--) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu1[c], wincon->hmin);

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midx[cs] * vel[c];
    if(bmask[xp1] == 1) {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
              (umid[c] * fluxc -
               umid[xm1] * dz[xm1] * midx[xm1s] * vel[xm1]);

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = umid[xm1] * dz[c] * midx[cs];
      d2 = umid[c] * dz[xp1] * midx[xp1s];
      ime1 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[cs] * (d1 * vel[c] - d2 * vel[xp1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
  	      (umid[c] * (dz[xp1] * midx[xp1s] * nvel[xp1] + fluxc) -
               umid[xm1] * dz[xm1] * midx[xm1s] * vel[xm1]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[cs] * umid[xm1] * dz[c] * midx[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[yp1] = nu1[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[yp1] == 1) {
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * fluxc - vmid[c] * dz[ym1] * vel[ym1]);

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = vmid[c] * dz[c];
      d2 = vmid[yp1] * dz[yp1];
      ime2 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
      ime2 = 0.5 * wincon->u1c1[cs] * (d1 * vel[c] - d2 * vel[yp1]);
    }
    else { 
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * (fluxc + dz[yp1] * nvel[yp1]) -
	       vmid[c] * dz[ym1] * vel[ym1]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u1c1[cs] * vmid[c] * dz[c];
      exe2 = ime2 * vel[c];
    }
    if(bmask[c] == 4.0) {
      /* Get the contribution from w.du/dz (note the wmid stagger)   */
      wu1 = 0.5 * wmid[zm1] * (vel[c] + nvel[zm1]) - wmid[c] * vel[c];
      imz = 0.0;
    }
    else {
      /* Get the contribution from w.du/dz (note the wmid stagger)   */
      wu1 = -0.5 * (wmid[c] * vel[zp1] -
		    wmid[zm1] * (vel[c] + nvel[zm1]));

      /* Get the implicit multiplier for w.du/dz                     */
      imz = 0.5 * wmid[c];
      exz = imz * vel[c];
    }
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2 + imz) / dzd;
    adv = (ush1h2 + u1u2hs + wu1) / dzd;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv -= ((exe1 + exe2 + exz) / dzd);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[cs] * vel[c] * vel[c] * midx[cs] +
           wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx[cs];
    /* Note : nu1 is a copy of u1b, but may be adjusted for thin     */
    /* layers, hence use nu1 for u1b here.                           */
    velb = windat->nu1[c];
    windat->nu1[c] = (velb + dtu * adv ) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = (u1b-dtu.w.du/dz) - nu1       */
    if(wincon->momsc & EXPLICIT) {
      adv -= ((wu1 - exz) / dzd);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
    }
    else {
      immp = 1.0 + dtu * imz / dzd;
      adv = (velb + dtu * wu1 / dzd) / immp;
      wincon->u1inter[cs] += ((adv - windat->nu1[c]) * windat->dzu1[c] *
			      windat->dtf / dtu);
      /*
      adv -= ((exe1 + exe2 + wu1 - exz) / dzd);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
      */

      /* Fill cells above the free surface to the top of the grid    */
      if(bmask[zp1] == 2) {
	cs = c;
	while (c != zp1) {
	  windat->nu1[zp1] = windat->nu1[cs];
	  c = zp1;
	  zp1 = window->zp1[c];
	}
      }
    }
  }
}

/* END advect_u1_3d_ang_flux_b()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* forward direction. Formulated using the flux form.                */
/*-------------------------------------------------------------------*/
void advect_u1_3d_ang_flux_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double fluxc;             /* Flux through cell face                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midx;             /* Depth at the e1 face                  */
  double *dz;               /* Cell thickness at cell face           */
  double dzd;               /* Cell thickness divisor                */
  double h2mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u2au1 = wincon->w7;
  vel = wincon->w8;
  area = window->cellarea;
  midx = wincon->mdx;
  dz = windat->dzu1;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu1;
  memcpy(vel, windat->u1, window->sgsiz * sizeof(double));
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required                              */
  if (wincon->thin_merge) {
    set_thin_u1(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
    /*
    wtop = wincon->d3;
    reset_thin_wtop(window, windat, wincon, wincon->kth_e1, 
		    wincon->nkth_e1, wtop);
    */
  } else {
    cells = wincon->s1;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  dtu = windat->dtf + windat->dtb;
  windat->dtu1 = windat->dt;

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
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e1; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_e1[cc]; /* 2D coordinate                     */
    cb = window->bot_e1[cc];    /* Bottom e1 coordinate              */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];         /* w : top set to wtop             */
      c = window->zm1[c];
    }
    /* Set the bottom boundary condition. Find the cell centered     */
    /* bottom coordinate and set w to wbot.                          */
    c = cb;                      /* e1 face centered bottom cell     */
    while (c != window->zm1[c])
      c = window->zm1[c];        /* Cell centered sediment location  */
    cb = window->zp1[c];         /* cell centered bottom location    */
    wvel[cb] = windat->wbot[cs]; /* w : bottom set to wbot           */
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e1; cc++) {
    c = window->w3_e1[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xm1s = window->m2d[xm1];
    h2mid = 0.5 * (window->h1au2[cs] + window->h1au2[xm1s]);
    umid[c] = 0.5 * area[cs] * (vel[c] + vel[xp1]);
    vmid[c] = 0.5 * h2mid * h2mid * 
             (windat->u2[c] * Ds[cs] + windat->u2[xm1] * Ds[xm1s]);
    wmid[c] = 0.5 * (wvel[c] + wvel[xm1]);
  }
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xm1 = window->xm1[c];
    wsur[c] = 0.5 * (wtop[c] + wtop[xm1]);
  }
  /* Set a no-gradient on u1 at the bottom and wmid at the bottom    */
  /* equal to the mean of wbot[].                                    */
  for(cc = 1; cc <= vcs; cc++) {
    c = wincon->i5[cc];
    cs = window->m2d[c];
    xm1 = window->xm1[cs];
    zm1 = window->zm1[c];
    wmid[c] = 0.5 * (windat->wbot[cs] + windat->wbot[xm1]);
    vel[zm1] = vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1 + 1; cc <=  window->nbpte1; cc++) {
    c = window->bpte1[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e1,
	      open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
	      open[n]->cyc_e1, windat->nu1, windat->u1, 
	      windat->u1b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to3_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1,
	      windat->u1, windat->u1b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u1;
  }
  if(wincon->momsc & ADVECT_FORM) {
    dz = wincon->w5;
    for(cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      dz[c] = 1.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  /* Surface layers : nu1[zp1] = nu1[c] from no-slip condition.      */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = 1; cc <= vcs; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu1[c], wincon->hmin);

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xm1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midx[cs] * vel[c];
    if(bmask[xm1]) {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
              (umid[c] * dz[xp1] * midx[xp1s] * vel[xp1] -
               umid[xm1] * fluxc);

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = umid[c] * dz[c] * midx[cs];
      d2 = umid[xm1] * dz[xm1] * midx[xm1s];
      ime1 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
	      (umid[c] * dz[xp1] * midx[xp1s] * vel[xp1] -
               umid[xm1] * (fluxc + dz[xm1] * midx[xm1s] * nvel[xm1]));

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c] * dz[c] * midx[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[ym1] = nu1[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[ym1]) {
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * dz[yp1] * vel[yp1] - vmid[c] * fluxc);

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = vmid[yp1] * dz[c];
      d2 = vmid[c] * dz[ym1];
      ime2 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * dz[yp1] * vel[yp1] -
	       vmid[c] * (fluxc + dz[ym1] * nvel[ym1]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u1c1[cs] * vmid[yp1] * dz[c];
      exe2 = ime2 * vel[c];
    }

    /* Get the contribution from w.du/dz                             */
    /*wu1 = -0.5 * (wsur[cs] * vel[c] - wmid[c] * vel[zm1]);*/
    /*wu1 = 0.5 * wmid[c] * vel[zm1] - wsur[cs] * vel[c];*/
    wu1 = 0.5 * wmid[c] * (vel[zm1] + vel[c]) - wsur[cs] * vel[c];

    /* Get the implicit multiplier for w.du/dz                       */
    /*imz = 0.5 * (wmid[c] - wsur[cs]);*/
    /*imz = 0.5 * wmid[c];*/
    imz = 0.0;
    exz = 0.5 * wmid[c] * vel[c];
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz) / dzd;
    adv = (ush1h2 + u1u2hs + wu1) / dzd;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += (exe1 + exe2 + exz) / dzd;
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[cs] * vel[c] * vel[c] * midx[cs] +
           wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx[cs];
    /* Note : nu1 is a copy of u1b, but may be adjusted for thin     */
    /* layers, hence use nu1 for u1b here.                           */
    velb = windat->nu1[c];
    windat->nu1[c] = (velb + dtu * adv ) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = (u1b-dtu.w.du/dz) - nu1       */
    if(wincon->momsc & EXPLICIT) {
      adv -= ((exz + wu1) / dzd);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
    }
    else {
      immp = 1.0 - dtu * imz / dzd;
      adv = (velb + dtu * wu1 / dzd) / immp;
      wincon->u1inter[cs] += ((adv - windat->nu1[c]) * windat->dzu1[c] *
			      windat->dtf / dtu);
      
      /* Fill the cells above the free surface to the top of grid    */
      cs = c;
      while (c != zp1) {
	windat->nu1[zp1] = windat->nu1[cs];
	c = zp1;
	zp1 = window->zp1[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sub - surface cells : nu1[zp1] is known                         */
  for (cc = vcs + 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu1[c], wincon->hmin);

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xm1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midx[cs] * vel[c];
    if(bmask[xm1]) {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
              (umid[c] * dz[xp1] * midx[xp1s] * vel[xp1] -
               umid[xm1] * fluxc) / dzd;

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = umid[c];
      d2 = umid[xm1] * dz[xm1] * midx[xm1s] / dzd;
      ime1 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[cs] * 
              (umid[c] * dz[xp1] * midx[xp1s] * vel[xp1] -
               umid[xm1] * (fluxc + dz[xm1] * midx[xm1s] * nvel[xm1])) / dzd;

      /* Get the implicit multiplier for u.du/dx                       */
      ime1 = 0.5 * wincon->u1c1[cs] * umid[c];
      exe1 = ime1 * fluxc / dzd;
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[ym1] = nu1[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[ym1]) {
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * dz[yp1] * vel[yp1] - vmid[c] * fluxc) / dzd;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = vmid[yp1];
      d2 = vmid[c] * dz[ym1] / dzd;
      ime2 = 0.5 * wincon->u1c1[cs] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u1c1[cs] *
	      (vmid[yp1] * dz[yp1] * vel[yp1] -
	       vmid[c] * (fluxc + dz[ym1] * nvel[ym1])) / dzd;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u1c1[cs] * vmid[yp1];
      exe2 = ime2 * fluxc / dzd;
    }

    /* Get the contribution from w.du/dz                             */
    wu1 = -0.5 * (wmid[zp1] * (nvel[zp1] + vel[c]) -
		  wmid[c] * vel[zm1]) / dzd;

    /* Get the implicit multiplier for w.du/dz                       */
    imz = 0.5 * wmid[c] / dzd;
    exz = imz * vel[c];

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz);
    adv = ush1h2 + u1u2hs + wu1;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += (exe1 + exe2 + exz);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[cs] * vel[c] * vel[c] * midx[cs] +
           wincon->u1c4[cs] * u2au1[c] * u2au1[c] * midx[cs];
    velb = windat->nu1[c];
    windat->nu1[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = (u1b-dtu.w.du/dz) - nu1       */
    if(wincon->momsc & EXPLICIT) {
      adv -= (exz + wu1);
      wincon->u1inter[cs] += (windat->dtf * adv * windat->dzu1[c]);
    }
    else {
      immp = 1.0 - dtu * (imz);
      adv = (velb + dtu * wu1) / immp;
      wincon->u1inter[cs] += ((adv - windat->nu1[c]) * windat->dzu1[c] *
			      windat->dtf / dtu);
    }
  }
}

/* END advect_u1_3d_ang_flux_f()                                     */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U2 ADVECTION ROUTINES                                             */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 3D mode      */
/* using sub-time stepping to render the scheme unconditionally      */
/* stable.                                                           */
/*-------------------------------------------------------------------*/
int advect_u2_3d(geometry_t *window, /* Window geometry             */
		 window_t *windat,   /* Window data                 */
		 win_priv_t *wincon  /* Window constants            */
  )
{
  int c, cc;                /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int xpym1;                /* Sparse coordinate at (i+1,j-1)        */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int ymzp1;                /* Sparse coordinate at (j-1,k-1)        */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double trem;              /* Time remaining in leapfrog step       */
  double tremf;             /* Time remaining in forward step        */
  double u1val;             /* Cell centered u1 value                */
  double u2val;             /* Cell cornered u2 value                */
  double wval;              /* Face centered vertical velocity       */
  double *u1au2;            /* u1 value at the e2 face               */
  double h2;                /* Face centered grid spacing            */
  double h1av, h2av;        /* Mean cell spacings x and y directions */
  double dzav;              /* Mean cell spacings in the z direction */
  double m, m1, m2;         /* Dummy variables                       */
  double d1, d2, d3;        /* Dummy variables                       */
  double *vel;              /* Velocity to use in spatial gradients  */
  double *area;             /* Cell area                             */
  double *ush1h2;           /* u2 * u2 * h1 * h2                     */
  double *u1u2hs;           /* u1 * u2 * h2 * h2                     */
  double *wu1;              /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *vface;            /* wvel centered on the grid face        */
  double *cn;               /* Vertical courant number               */
  double *velbuf;           /* u2 buffer for vertical flux           */
  double *dzu2;             /* Pointer for thickness                 */
  double *div;              /* Momentum divergence                   */
  double midy;              /* Depth at the e2 face                  */
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
  /* Assign pointers */
  ush1h2 = wincon->w1;
  u1u2hs = wincon->w2;
  wu1 = wincon->w3;
  wvel = wincon->w4;
  velbuf = wincon->w5;
  div = wincon->w6;
  area = window->cellarea;
  u1au2 = wincon->w7;
  vel = windat->u2;
  wq = wincon->d1;
  if (wincon->momsc & ADVECT_FORM) {
    dzu2 = wincon->w8;
    for (cc = 1; cc <= window->n3_e2; cc++) {
      c = window->w3_e2[cc];
      dzu2[c] = 1.0;
    }
  } else
    dzu2 = windat->dzu2;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required                              */
  vcs = wincon->vcs;
  if (wincon->thin_merge) {
    set_thin_u2(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s2;
    vc = wincon->vc;
  }
  if(wincon->dolin_u2) {
    linear_bdry_cell(window, windat, wincon, cells, vcs,
		     wincon->linmask_u2, wincon->i6);
    cells = wincon->s4;
    vc = wincon->acl;
    vcs = wincon->aclS;
  }
  set_surf_cells(window, windat, wincon, 0);

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
    for (cc = 1; cc <= vc; cc++) {
      c = cells[cc];            /* Wet cell to process               */
      cs = window->m2d[c];      /* 2D cell corresponding to 3D cell  */
      xp1 = window->xp1[c];
      ym1 = window->ym1[c];
      /* xpym1=window->ym1[xp1]; */
      xpym1 = window->xpym1[c];
      zp1 = window->zp1[c];
      ymzp1 = window->zp1[ym1];

      /* Minimum time-step due to x velocity                         */
      vlc = 0.25 * (windat->u1[ym1] + windat->u1[xpym1] +
                    windat->u1[c] + windat->u1[xp1]);
      if (fabs(vlc) < minval)
        vlc = minval;
      d1 = sf * fabs(window->h1au2[cs] / vlc);
      if (d1 < dtm) {
        dtm = d1;
        tp = 'u';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }

      /* Minimum time-step due to y velocity                         */
      vlc = windat->u2[c];
      if (fabs(vlc) < minval)
        vlc = minval;
      d2 = sf * fabs(window->h2au2[cs] / vlc);
      if (d2 < dtm) {
        dtm = d2;
        tp = 'v';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
      /* Minimum time-step due to z velocity                         */
      vlc =
        (cc <=
         vcs) ? windat->wtop[cs] +
        windat->wtop[window->ym1[cs]] : windat->w[zp1] + windat->w[ymzp1];
      vlc = 0.25 * (vlc + windat->w[c] + windat->w[ym1]);
      if (cc <= vcs) wq[cs] = dtu * vlc / (windat->dzu2[c] * wincon->mdx[cs]);
      if (wincon->stab & SUB_STEP_NOSURF && cc <= vcs)
        vlc = 0.0;
      if (fabs(vlc) < minval)
        vlc = minval;
      d3 =
        wsf * fabs(max(windat->dzu2[c] * wincon->mdy[cs], wincon->hmin) /
                   vlc);
      if (d3 < dtm) {
        dtm = d3;
        tp = 'w';
        vs = vlc;
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
    }
  }
  if (dtm != dtu) {
    c = ceil((windat->dtb + dtu) / (2 * dtm));
    dtm = (windat->dtb + dtu) / (2.0 * (double)c);
    dtu = 2.0 * dtm;
    hd_warn
      ("Sub-time stepping in u2 (%c=%6.3f) at %8.3f days (%3d %3d %3d): dt=%5.2f\n",
       tp, vs, windat->t / 86400.0, ii, jj, kk, dtm);
    if ((int)(dtu / dtm) > itermax) {
      hd_quit_and_dump
	("advect_u2_3d: maximum number of sub-steps (%d) exceeded.\n",
	 itermax);
      return(1);
    }
  }
  else
    dtu = dtm + windat->dtb;
  windat->dtu2 = dtm;
  if (wincon->u2_f & ADVECT)
    return(0);

  if (wincon->momsc & LAGRANGE) {
    semi_lagrange_u2(window, windat, wincon, cells, vcs, vc);
    if (wincon->nbl2) {
      vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
		window->bpte2, window->bine2, wincon->slip);
      blend_vel(window, windat, wincon, U23D, vel);
    }
    return(0);
  }

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity.                                                       */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {

    /* Get the 2D and 3D sparse coordinates                          */
    cs = c = window->w2_e2[cc]; /* 2D coordinate                     */
    cb = window->bot_e2[cc];    /* Bottom coordinate                 */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs];
      c = window->zm1[c];
    }

    c = cb;
    while (c != window->zm1[c])
      c = window->zm1[c];
    cb = window->zp1[c];
    wvel[cb] = windat->wbot[cs];
  }
  /* Set wvel at the surface for ghost cells. If window partitions   */
  /* lie next to ghost cells then these are not captured in the loop */
  /* above but are required for face centered vertical velocity.     */
  /* Note that wvel[cb] always = 0 and is not explicitly set.        */
  for (cc = window->a2_e2 + 1; cc <= window->n2_e2; cc++) {
    cs = c = window->w2_e2[cc];
    top = windat->eta[cs];
    while (c != window->zm1[c] && window->gridz[c] > top) {
      wvel[c] = windat->wtop[cs];
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  /* Set the first timestep for leapfrog                             */
  trem = windat->dt;
  tremf = windat->dtf;
  /* Loop over the time interval until no time remains               */
  while (trem > 0) {
    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Initialise                                                  */
      memset(ush1h2, 0, window->sgsiz * sizeof(double));
      memset(u1u2hs, 0, window->sgsiz * sizeof(double));
      memset(wu1, 0, window->sgsiz * sizeof(double));
      if (wincon->momsc & ZERO_DRYK)
	set_surf_cond(window, vel, 0.0, window->sur_e2, window->b2_e2);

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell centres (u2*u2*h1*h2*dz)     */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          yp1 = window->yp1[c];
          u2val = 0.5 * (vel[yp1] + vel[c]);
          /* SIGMA : Multiply by the cell centre depth               */
          if (u2val > 0)
            m = dzu2[c] * vel[c] * wincon->mdy[cs];
          else
            m =
              dzu2[yp1] * vel[yp1] * wincon->mdy[window->m2d[yp1]];
          ush1h2[c] = m * u2val * area[cs];
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];

          cs = window->m2d[c];
          yp1 = window->yp1[c];
          u2val = 0.5 * (vel[yp1] + vel[c]);

          /* SIGMA : Multiply by the cell centre depth               */
          m = 0.5 * (dzu2[yp1] * vel[yp1] * 
                     wincon->mdy[window->m2d[yp1]] +
                     dzu2[c] * vel[c] * wincon->mdy[cs]);
          ush1h2[c] = m * u2val * area[cs];	
	}
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          yp1 = window->yp1[c];
          vface[c] = 0.5 * (vel[yp1] + vel[c]);
          cn[c] = vface[c] * dtu / window->h2acell[cs];
          wu1[c] = dzu2[yp1] * vel[yp1] *
            wincon->mdy[window->m2d[yp1]];
        }
	/* Ghost cells 2 cells in */
        for (cc = window->b3_e2 + 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          ym1 = window->ym1[c];
          wu1[ym1] = dzu2[c] * vel[c] *
            wincon->mdy[window->m2d[c]];
	}

        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
        van_leer_do(ush1h2, wu1, vface, cn, 1, window->n3_e2,
                    window->w3_e2, window->yp1, window->ym1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          ush1h2[c] *= vface[c] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the advective flux at cell corners (u1*u2*h2*h2*dz)     */
      /* First order upstream scheme                                 */
      if (wincon->momsc & ORDER1) {
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          cs = window->m2d[c];
          h2 = 0.5 * (window->h2au1[window->ym1[cs]] + window->h2au1[cs]);
          u1val = 0.5 * (windat->u1[ym1] * wincon->mdx[ym1] +
                         windat->u1[c] * wincon->mdx[cs]);

          if (u1val > 0)
            m = dzu2[xm1] * vel[xm1];
          else
            m = dzu2[c] * vel[c];
          u1u2hs[c] = m * h2 * h2 * u1val;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];

          xm1 = window->xm1[c];
          ym1 = window->ym1[c];
          cs = window->m2d[c];
          xm1s = window->xm1[cs];
          ym1s = window->ym1[cs];
          h2 = 0.5 * (window->h2au1[ym1s] + window->h2au1[cs]);
          h2av =
            window->h2au1[ym1s] / (window->h2au1[ym1s] +
                                   window->h2au1[cs]);
          h1av =
            window->h1au2[xm1s] / (window->h1au2[xm1s] +
                                   window->h1au2[cs]);
          /* SIGMA : Multiply by the cell corner depth               */
          u1val = h2av * (windat->u1[c] * wincon->Ds[cs] -
                          windat->u1[ym1] * wincon->Ds[ym1s]) +
            windat->u1[ym1] * wincon->Ds[ym1s];
          m1 = dzu2[xm1] * vel[xm1];
          m2 = dzu2[c] * vel[c];
          m = h1av * (m2 - m1) + m1;
          u1u2hs[c] = m * h2 * h2 * u1val;
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sgsiz * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          ym1 = window->ym1[c];
          vface[c] = 0.5 * (windat->u1[ym1] + windat->u1[c]);
          cn[c] = 2.0 * vface[c] * dtu /
            (window->h1au1[cs] + window->h1au1[window->ym1[cs]]);
          wu1[c] = dzu2[c] * vel[c] * wincon->mdy[cs];
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wu1.                                          */
        van_leer_do(u1u2hs, wu1, vface, cn, 1, window->n3_e2,
                    window->w3_e2, window->xp1, window->xm1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= window->n3_e2; cc++) {
          c = window->w3_e2[cc];
          cs = window->m2d[c];
          u1u2hs[c] *= vface[c] * area[cs];
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the vertical advective fluxes (u2*w)                    */
      /* solution.                                                   */
      /* First set the bottom boundary condition for u2 and w        */
      if (wincon->momsc & ZERO_DRYK)
	set_surf_cond(window, vel, -1.0, window->sur_e2, window->b2_e2);

      memcpy(velbuf, vel, window->sgsiz * sizeof(double));
      for (cc = 1; cc <= window->a2_e2; cc++) {
        c = window->bot_e2[cc]; /* 3D bottom coordinate              */
        zm1 = window->zm1[c];   /* Sediment coordinate               */
        velbuf[zm1] = vel[c];   /* u2 : no-gradient across bottom    */
      }

      /* Surface vertical advective flux                             */
      for (cc = 1; cc <= vcs; cc++) {
        c = cells[cc];
        cs = window->m2d[c];
        ym1 = window->ym1[cs];
        wu1[c] = 0.5 * vel[c] * (windat->wtop[ym1] + windat->wtop[cs]);
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
        /* lower face).                                              */
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          ym1 = window->ym1[c];
          zm1 = window->zm1[c];
          wval = 0.5 * (wvel[ym1] + wvel[c]);
          if (wval > 0)
            wu1[zm1] = velbuf[zm1] * wval;
          else
            wu1[zm1] = velbuf[c] * wval;
        }
      }
      /* Second order centered scheme                                */
      else if (wincon->momsc & ORDER2) {
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          ym1 = window->ym1[c];
          zm1 = window->zm1[c];
          wval = 0.5 * (wvel[ym1] + wvel[c]);
          dzav = windat->dzu2[zm1] / (windat->dzu2[zm1] + windat->dzu2[c]);
          wu1[zm1] =
            wval * (dzav * (velbuf[c] - velbuf[zm1]) + velbuf[zm1]);
        }
      }
      /* Van Leer's scheme                                           */
      else if (wincon->momsc & VANLEER) {
        cn = wincon->w8;
        vface = wincon->w6;
        memset(wu1, 0, window->sgsiz * sizeof(double));
        /* Get the vertical Courant number and vertical velocity at  */
        /* the e1 face.                                              */
        for (cc = 1; cc <= window->n3_t; cc++) {
          c = window->w3_t[cc];
          cs = window->m2d[c];
          ym1 = window->ym1[c];
          dzav = 0.5 * (windat->dzu2[c] + windat->dzu2[window->zm1[c]]);
          vface[c] = 0.5 * (wvel[ym1] + wvel[c]);
          cn[c] = vface[c] * dtu / (dzav * wincon->Ds[cs]);
        }
        /* Get the u1 velocity at the grid face and store in the     */
        /* dummy array wvel.                                         */
        van_leer_do(wvel, velbuf, vface, cn, 1, vc, cells, window->zp1,
                    window->zm1);

        /* Multiply by vertical velocity and shift one cell down     */
        for (cc = 1; cc <= vc; cc++) {
          c = cells[cc];
          zm1 = window->zm1[c];

          wu1[zm1] = vface[c] * wvel[c];
        }
      }

      /* Bottom vertical advective flux                              */
      for (cc = 1; cc <= window->v2_e2; cc++) {
        c = window->bot_e2[cc]; /* Bottom coordinate                 */
        cs = window->w2_e2[cc];
        ym1 = window->ym1[cs];
        zm1 = window->zm1[c];
        wu1[zm1] =
          0.5 * velbuf[c] * (windat->wbot[cs] + windat->wbot[ym1]);
      }

      /*-------------------------------------------------------------*/
      /* Get the momentum divergence and metric terms.               */
      for (cc = 1; cc <= vc; cc++) {
        c = cells[cc];
        cs = window->m2d[c];
        ym1 = window->ym1[c];
        xp1 = window->xp1[c];
        zm1 = window->zm1[c];
        midy = wincon->mdy[cs];

        /* SIGMA : Multiply by cell depth                            */
        d1 = 
          wincon->u2c1[cs] * (ush1h2[c] - ush1h2[ym1] + u1u2hs[xp1] -
                              u1u2hs[c]) / max(dzu2[c], wincon->hmin) +
          wincon->u2c3[cs] * vel[c] * vel[c] * midy +
          wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy;

	if (wincon->porusplate)
	  d1 += porus_plate_e2(window, windat, wincon, c);

        /* Integrate the e1 dispersion term vertically and with time */
        /* (i.e. over the sub-timestepping; windat->dtf).            */
	wincon->u2inter[cs] += (d1 * dtm * windat->dzu2[c]);

        /* Add the vertical advective flux. Note : the difference of */
        /* this flux is wu1[c]-wu1[zm1] due to the stagger imposed,  */
        /* not wu1[zp1]-wu1[c].                                      */
	if (!(wincon->momsc & WIMPLICIT))
	  d1 -= (wu1[c] - wu1[zm1]) / max(windat->dzu2[c], wincon->hmin);

        /* Add to new u1 value                                       */
	div[c] = d1 * dtu;
        /*windat->nu2[c] += d1 * dtu;*/
      }

      /* Filter the divergence if required                           */
      if (wincon->filter & ADVECT) {
	vel2D_lbc(div, window->nbpte2, window->nbe2,
		  window->bpte2, window->bine2, wincon->slip);
	shapiro(window, div, velbuf, cells, vc, 1, XDIR);
	shapiro(window, div, velbuf, cells, vc, 1, YDIR);
      }

      /* No advection at the blend-zone boundaries for stability     */
      for (ii = 0; ii < wincon->nbl1; ii++) {
	blend_t *blend = wincon->ble1[ii];
	for (cc = 1; cc <= blend->nlim; cc++) {
	  c = blend->lim[cc];
	  div[c] = 0.0;
	}
	for (cc = 1; cc <= blend->nlimS; cc++) {
	  c = blend->lim[cc];
	  wincon->u2inter[c] = 0.0;
	}
      }

      /* Get the u1 updated solution                                 */
      for (cc = 1; cc <= vc; cc++) {
        c = cells[cc];
        windat->nu2[c] += div[c];
      }

      /* Decrement the time remaining                                */
      trem -= dtu;
      tremf -= dtm;
      if (dtm)
	dtu = dtm * windat->dt / windat->dtf;
      if (trem < dtu)
        dtu = trem;
      if (tremf < dtm)
        dtm = tremf;
       if(!wincon->sigma)
	vel = windat->nu2;

       if (wincon->nbl2) {
	 vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
		   window->bpte2, window->bine2, wincon->slip);
	 blend_vel(window, windat, wincon, U23D, vel);
       }

     }
  }
  return(0);
}

/* END advect_u2_3d()                                          (slf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* backward direction. Formulated using the advective form.          */
/*-------------------------------------------------------------------*/
void advect_u2_3d_ang_adv_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
		  	    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int ymxp, ymzm;
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h2au2;            /* Cell width at e2 face                 */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midy;             /* Depth at the e1 face                  */
  double dz;                /* Cell thickness at cell face           */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int *ctp;                 /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u1au2 = wincon->w7;
  h2au2 = window->h2au2;
  vel = windat->u2;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu2;
  cells = wincon->s4;
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required */
  if (wincon->thin_merge) {
    set_thin_u2(window, windat, wincon);
    ctp = wincon->s3;
    vc = wincon->ncl;
  } else {
    vc = wincon->vc;
    ctp = wincon->s2;
  }
  vcs = wincon->vcs;
  reorder_cells(window, cells, ctp, vc, vcs, wincon->i6, 1);
  dtu = windat->dtf + windat->dtb;
  windat->dtu2 = windat->dt;

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity.                                                       */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {

    /* Get the 2D and 3D sparse coordinates */
    cs = c = window->w2_e2[cc]; /* 2D coordinate */
    cb = window->bot_e2[cc];    /* Bottom coordinate */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    vmid[c] = 0.25 * (2.0 * vel[c] + vel[yp1] + vel[ym1]) *
              window->h1au2[cs];
    umid[c] = windat->u1[c] * window->h2au1[cs];
  }
  /* Stagger wmid so that wtop occupies the surface cell             */
  for(cc = 1; cc <= vc; cc++) {
    c = ctp[cc];
    zm1 = window->zm1[c];
    wmid[zm1] = wvel[c];
  }
  /* Surface and bottom conditions for wmid                          */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cb = window->zm1[wincon->i6[cc]];
    cs = window->m2d[c];
    wmid[c] = wtop[cs];
    wmid[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the cells above the surface                                 */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    zp1 = window->zp1[c];
    bmask[c] = 4.0;
    while (c != zp1) {
      bmask[zp1] = 2.0;
      c = zp1;
      zp1 = window->zp1[c];
    }
  }
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2 + 1; cc <=  window->nbpte2; cc++) {
    c = window->bpte2[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e2,
	      open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
	      open[n]->cyc_e2, windat->nu2, windat->u2, 
	      windat->u2b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to3_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2,
	      windat->u2, windat->u2b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u2;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration : nu2[zp1] is known.                         */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = vc; cc >= 1; cc--) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    ymxp = window->xp1[ym1];
    ymzm = window->zm1[ym1];

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[yp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[yp1] == 1) {
      ush1h2 = 0.5 * wincon->u2c1[cs] * vmid[c] * 
	      (vel[c] * h2au2[cs] - vel[ym1] * h2au2[ym1s]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h2au2[cs];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = h2au2[yp1s];
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u2c1[cs] * vmid[c] * (d1 * vel[c] - d2 * vel[yp1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[cs] * vmid[c] * 
	(vel[c] * h2au2[cs] + nvel[yp1] * h2au2[yp1s] -
	 vel[ym1] * h2au2[ym1s]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c] * h2au2[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land or open      */
    /* boundaries nu2[xp1] = nu2[c] due to the no-slip condition.    */
    if(bmask[xp1] == 1) {
      u1u2hs = 0.25 * wincon->u2c1[cs] * 
	(umid[c] + umid[ym1]) * (vel[c] * h2au2[cs] - vel[xm1] * h2au2[xm1s]);
					   
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u2c1[cs] * 
	((umid[xp1] + umid[ymxp]) * nvel[xp1] * h2au2[xp1s] +
	 (umid[c] + umid[ym1]) *
	 (vel[c] * h2au2[cs] - vel[xm1] * h2au2[xm1s]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u2c1[cs] * (umid[xp1] + umid[ymxp]) * h2au2[cs];
      exe2 = ime2 * vel[c];
    }

    dz = max(windat->dzu2[c], wincon->hmin);
    /* Get the contribution from w.du/dz (note the wmid stagger)     */
    if(bmask[c] == 4) {
      d1 = wmid[c] + wmid[ym1] + wmid[zm1] + wmid[ymzm];
      wu1 = -0.125 * d1 * (vel[c] - (vel[c] + nvel[zm1])) / dz;

      /* Get the implicit multiplier for v.du/dy                       */
      imz = 0.125 * d1 /dz;
      exz = imz * vel[c];
    }
    else {
      wu1 = -0.25 * ((wmid[c] + wmid[ym1]) * (vel[zp1] - vel[c]) +
		     (wmid[zm1] + wmid[ymzm]) * nvel[zm1]) / dz;
    
      /* Get the implicit multiplier for v.du/dy                       */
      imz = 0.25 * (wmid[zm1] + wmid[ymzm]) / dz;
      exz = imz * vel[c];
    }
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2 + imz);
    adv = ush1h2 + u1u2hs + wu1;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv -= (exe1 + exe2 + exz);
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[cs] * vel[c] * vel[c] * midy[cs] +
           wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy[cs];
    /* Note : nu2 is a copy of u2b, but may be adjusted for thin     */
    /* layers, hence use nu2 for u2b here.                           */
    velb = windat->nu2[c];
    windat->nu2[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.dv/dx + v.dv/dy) = (u2b-dtu.w.dv/dz) - nu2       */
    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      adv -= (wu1 - exz);
      wincon->u2inter[cs] += (windat->dtf * adv * windat->dzu2[c]);
    }
    else {
      immp = 1.0 + dtu * imz;
      adv = (velb + dtu * wu1) / immp;
      wincon->u2inter[cs] += ((adv - windat->nu2[c]) * windat->dzu2[c] * 
		  	      windat->dtf / dtu);

      /* Fill cells above the free surface to the top of the grid    */
      if(bmask[zp1] == 2) {
	cs = c;
	while (c != zp1) {
	  windat->nu2[zp1] = windat->nu2[cs];
	  c = zp1;
	  zp1 = window->zp1[c];
	}
      }
    }
  }
}

/* END advect_u2_3d_ang_adv_b()                                (sbf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* forward direction. Formulated using the advective form.           */
/*-------------------------------------------------------------------*/
void advect_u2_3d_ang_adv_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int ymxp, ymzm;
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h2au2;            /* Cell width at e2 face                 */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midy;             /* Depth at the e1 face                  */
  double dz;                /* Cell thickness at cell face           */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u1au2 = wincon->w7;
  vel = wincon->w8;
  h2au2 = window->h2au2;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu2;
  memcpy(vel, windat->u2, window->sgsiz * sizeof(double));
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required */
  vcs = wincon->vcs;
  if (wincon->thin_merge) {
    set_thin_u2(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
  } else {
    cells = wincon->s2;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  dtu = windat->dtf + windat->dtb;
  windat->dtu2 = windat->dt;

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity.                                                       */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {

    /* Get the 2D and 3D sparse coordinates */
    cs = c = window->w2_e2[cc]; /* 2D coordinate */
    cb = window->bot_e2[cc];    /* Bottom coordinate */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];
      c = window->zm1[c];
    }

    c = cb;
    while (c != window->zm1[c])
      c = window->zm1[c];
    cb = window->zp1[c];
    wvel[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    vmid[c] = 0.25 * (2.0 * vel[c] + vel[yp1] + vel[ym1]) *
              window->h1au2[cs];
    umid[c] = windat->u1[c] * window->h2au1[cs];
  }
  /* Stagger wmid so that wtop occupies the surface cell             */
  for(cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    zm1 = window->zm1[c];
    wmid[zm1] = wvel[c];
  }
  /* Surface and bottom conditions for wmid                          */
  for(cc = 1; cc <= vcs; cc++) {
    c = cells[cc];
    cb = window->zm1[wincon->i6[cc]];
    cs = window->m2d[c];
    wmid[c] = wtop[cs];
    wmid[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2 + 1; cc <=  window->nbpte2; cc++) {
    c = window->bpte2[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e2,
	      open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
	      open[n]->cyc_e2, windat->nu2, windat->u2, 
	      windat->u2b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to3_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2,
	      windat->u2, windat->u2b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u2;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = vcs + 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    ymxp = window->xp1[ym1];
    ymzm = window->zm1[ym1];

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[ym1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[ym1]) {
      ush1h2 = 0.5 * wincon->u2c1[cs] * vmid[c] * 
	      (vel[yp1] * h2au2[yp1s] - vel[c] * h2au2[cs]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h2au2[cs];
      d2 = h2au2[yp1s];
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u2c1[cs] * vmid[c] * (d1 * vel[c] - d2 * vel[ym1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[cs] * vmid[c] * 
	(vel[yp1] * h2au2[yp1s] - 
	 (vel[c] * h2au2[cs] + nvel[ym1] * h2au2[ym1s]));

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c] * h2au2[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land or open      */
    /* boundaries nu2[xm1] = nu2[c] due to the no-slip condition.    */
    if(bmask[xm1]) {
      u1u2hs = 0.25 * wincon->u2c1[cs] * 
	(umid[xp1] + umid[ymxp]) * (vel[xp1] * h2au2[xp1s] - 
				    vel[c] * h2au2[cs]);
					   
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u2c1[cs] * 
	((umid[xp1] + umid[ymxp]) * 
	 (vel[xp1] * h2au2[xp1s] - vel[c] * h2au2[cs]) -
	 (umid[c] + umid[ym1]) * nvel[xm1] * h2au2[xm1s]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u2c1[cs] * (umid[c] + umid[ym1]) * h2au2[cs];
      exe2 = ime2 * vel[c];
    }

    /* Get the contribution from w.dv/dz                             */
    dz = max(windat->dzu2[c], wincon->hmin);
    wu1 = -0.25 * ((wmid[c] + wmid[ym1]) * nvel[zp1] +
		   (wmid[zm1] + wmid[ymzm]) * (vel[c] - vel[zm1])) / dz;
    
    /* Get the implicit multiplier for v.dv/dz                       */
    imz = 0.25 * (wmid[c] + wmid[ym1]) / dz;
    exz = imz * vel[c];
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz);
    adv = ush1h2 + u1u2hs + wu1;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2 + exz;
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[cs] * vel[c] * vel[c] * midy[cs] +
           wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy[cs];
    velb = windat->nu2[c];
    windat->nu2[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.dv/dx + v.dv/dy) = (u2b-dtu.w.dv/dz) - nu2       */
    if(wincon->momsc & EXPLICIT) {
      adv -= (exz + wu1);
      wincon->u2inter[cs] += (windat->dtf * adv * windat->dzu2[c]);
    }
    else {
      immp = 1.0 - dtu * (imz);
      adv = (velb + dtu * wu1) / immp;
      wincon->u2inter[cs] += ((adv - windat->nu2[c]) * windat->dzu2[c] * 
			      windat->dtf / dtu);
    }
  }
}

/* END advect_u2_3d_ang_adv_f()                                (saf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* backward direction. Formulated using the flux form.               */
/*-------------------------------------------------------------------*/
void advect_u2_3d_ang_flux_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
		  	    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double fluxc;             /* Flux through cell face                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midy;             /* Depth at the e1 face                  */
  double *dz;               /* Cell thickness at cell face           */
  double dzd;               /* Cell thickness divisor                */
  double h1mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int *ctp;                 /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u1au2 = wincon->w7;
  area = window->cellarea;
  vel = windat->u2;
  midy = wincon->mdy;
  dz = windat->dzu2;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu2;
  cells = wincon->s4;
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required */
  if (wincon->thin_merge) {
    set_thin_u2(window, windat, wincon);
    ctp = wincon->s3;
    vc = wincon->ncl;
    /*
    wtop = wincon->d3;
    reset_thin_wtop(window, windat, wincon, wincon->kth_e2, 
		    wincon->nkth_e2, wtop);
    */
  } else {
    vc = wincon->vc;
    ctp = wincon->s2;
  }
  vcs = wincon->vcs;
  reorder_cells(window, cells, ctp, vc, vcs, wincon->i6, 1);
  dtu = windat->dtf + windat->dtb;
  windat->dtu2 = windat->dt;

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity.                                                       */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {

    /* Get the 2D and 3D sparse coordinates */
    cs = c = window->w2_e2[cc]; /* 2D coordinate */
    cb = window->bot_e2[cc];    /* Bottom coordinate */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];
      c = window->zm1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    ym1s = window->m2d[ym1];
    zm1 = window->zm1[c];
    h1mid = 0.5 * (window->h2au1[cs] + window->h2au1[ym1s]);
    umid[c] = 0.5 * h1mid * h1mid * 
             (windat->u1[c] * Ds[cs] + windat->u1[ym1] * Ds[ym1s]);
    vmid[c] = 0.5 * area[cs] * (vel[c] + vel[yp1]);
    /* Stagger wmid so that wtop occupies the surface cell           */
    wmid[zm1] = 0.5 * (wvel[c] + wvel[ym1]);
  }
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    cb = window->zm1[wincon->i6[cc]];
    cs = window->m2d[c];
    ym1 = window->ym1[cs];
    wmid[c] = 0.5 * (wtop[cs] + wtop[ym1]);
    wmid[cb] = 0.5 * (windat->wbot[cs] + windat->wbot[ym1]);
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the cells above the surface                                 */
  for(cc = 1; cc <= vcs; cc++) {
    c = ctp[cc];
    zp1 = window->zp1[c];
    bmask[c] = 4.0;
    while (c != zp1) {
      bmask[zp1] = 2.0;
      c = zp1;
      zp1 = window->zp1[c];
    }
  }
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2 + 1; cc <=  window->nbpte2; cc++) {
    c = window->bpte2[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e2,
	      open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
	      open[n]->cyc_e2, windat->nu2, windat->u2,
	      windat->u2b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to3_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2,
	      windat->u2, windat->u2b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u2;
  }
  if(wincon->momsc & ADVECT_FORM) {
    dz = wincon->w5;
    for(cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      dz[c] = 1.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration : nu2[zp1] is known.                         */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = vc; cc >= 1; cc--) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu2[c], wincon->hmin);

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[yp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midy[cs] * vel[c];
    if(bmask[yp1] == 1) {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
              (vmid[c] * fluxc -
               vmid[ym1] * dz[ym1] * midy[ym1s] * vel[ym1]);

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = vmid[ym1] * dz[c] * midy[cs];
      d2 = vmid[c] * dz[yp1] * midy[yp1s];
      ime1 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
 	      (vmid[c] * (dz[yp1] * midy[yp1s] * nvel[yp1] + fluxc) -
               vmid[ym1] * dz[ym1] * midy[ym1s] * vel[ym1]);

      /* Get the implicit multiplier for v.dv/dy                     */
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[ym1] * dz[c] * midy[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land or open      */
    /* boundaries nu2[xp1] = nu2[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[xp1] == 1) {
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * fluxc - umid[c] * dz[xm1] * vel[xm1]);

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = umid[c] * dz[c];
      d2 = umid[xp1] * dz[xp1];
      ime2 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * (fluxc + dz[xp1] * nvel[xp1]) -
	       umid[c] * dz[xm1] * vel[xm1]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u2c1[cs] * umid[c] * dz[c];
      exe2 = ime2 * vel[c];
    }

    /* Get the contribution from w.dv/dz                             */
    if(bmask[c] == 4.0) {
      /* Get the contribution from w.du/dz (note the wmid stagger)   */
      wu1 = 0.5 * wmid[zm1] * (vel[c] + nvel[zm1]) - wmid[c] * vel[c];
      imz = 0.0;
    }
    else {
      wu1 = -0.5 * (wmid[c] * vel[zp1] -
		    wmid[zm1] * (vel[c] + nvel[zm1]));

      /* Get the implicit multiplier for w.dv/dz                     */
      imz = 0.5 * wmid[c];
      exz = imz * vel[c];
    }
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2 + imz) / dzd;
    adv = (ush1h2 + u1u2hs + wu1) / dzd;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv -= ((exe1 + exe2 + exz) / dzd);
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[cs] * vel[c] * vel[c] * midy[cs] +
           wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy[cs];
    /* Note : nu2 is a copy of u2b, but may be adjusted for thin     */
    /* layers, hence use nu2 for u2b here.                           */
    velb = windat->nu2[c];
    windat->nu2[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.dv/dx + v.dv/dy) = (u2b-dtu.w.dv/dz) - nu2       */
    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      adv -= ((wu1 - exz) / dzd);
      wincon->u2inter[cs] += (windat->dtf * adv * windat->dzu2[c]);
    }
    else {
      immp = 1.0 + dtu * imz / dzd;
      adv = (velb + dtu * wu1 / dzd) / immp;
      wincon->u2inter[cs] += ((adv - windat->nu2[c]) * windat->dzu2[c] * 
		  	      windat->dtf / dtu);

      /* Fill cells above the free surface to the top of the grid    */
      if(bmask[zp1] == 2) {
	cs = c;
	while (c != zp1) {
	  windat->nu2[zp1] = windat->nu2[cs];
	  c = zp1;
	  zp1 = window->zp1[c];
	}
      }
    }
  }
}

/* END advect_u2_3d_ang_flux_b()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 3D mode      */
/* using an unconditionally stable angular centered derivative (see  */
/* Kowalik and Murty (1993) p51, 186) stepping through the grid in a */
/* forward direction. Formulated using the flux form.                */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void advect_u2_3d_ang_flux_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int cs, cb;               /* Surface / bottom sparse coordinate    */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xp1s, yp1s;           /* 2D sparse coordinate at i+1, j+1      */
  int xm1s, ym1s;           /* 2D sparse coordinate at i-1, j-1      */
  int zp1, zm1;             /* Sparse coordinate at k+1, k-1         */
  int vc, vcs;              /* Cells to process counters             */
  double dtu;               /* Leapfrog sub-time step to use         */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double imz;               /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double exz;               /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double fluxc;             /* Flux through cell face                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double wu1;               /* w * u1                                */
  double *wvel;             /* Vertical velocity                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *wmid;             /* Corner centered w velocity            */
  double *wsur;             /* Corner centered surface w velocity    */
  double *wtop;             /* Surface vertical velocity             */
  double *midy;             /* Depth at the e1 face                  */
  double *dz;               /* Cell thickness at cell face           */
  double dzd;               /* Cell thickness divisor                */
  double h1mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double top;               /* Surface level at the cell face        */
  double d1, d2;            /* Dummies                               */
  int *cells;               /* Cells to process                      */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  wmid = wincon->w3;
  wvel = wincon->w4;
  u1au2 = wincon->w7;
  vel = wincon->w8;
  area = window->cellarea;
  midy = wincon->mdy;
  dz = windat->dzu2;
  Ds = wincon->Ds;
  wsur = wincon->d2;
  wtop = windat->wtop;
  nvel = windat->nu2;
  memcpy(vel, windat->u2, window->sgsiz * sizeof(double));
  exe1 = exe2 = exz = 0.0;

  /*-----------------------------------------------------------------*/
  /* Adjust the thin layers if required */
  vcs = wincon->vcs;
  if (wincon->thin_merge) {
    set_thin_u2(window, windat, wincon);
    cells = wincon->s3;
    vc = wincon->ncl;
    /*
    wtop = wincon->d3;
    reset_thin_wtop(window, windat, wincon, wincon->kth_e2, 
		    wincon->nkth_e2, wtop);
    */
  } else {
    cells = wincon->s2;
    vc = wincon->vc;
  }
  vcs = wincon->vcs;
  dtu = windat->dtf + windat->dtb;
  windat->dtu2 = windat->dt;

  /*-----------------------------------------------------------------*/
  /* Precondition the vertical velocity array by setting all         */
  /* velocities above the surface coordinate equal to the surface    */
  /* vertical velocity, wtop. This must be performed over auxilliary */
  /* cells since these are used to get face centered vertical        */
  /* velocity.                                                       */
  memcpy(wvel, windat->w, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_e2; cc++) {

    /* Get the 2D and 3D sparse coordinates */
    cs = c = window->w2_e2[cc]; /* 2D coordinate */
    cb = window->bot_e2[cc];    /* Bottom coordinate */

    top = windat->eta[cs];
    while (c <= cb && window->gridz[c] > top) {
      wvel[c] = wtop[cs];
      c = window->zm1[c];
    }

    c = cb;
    while (c != window->zm1[c])
      c = window->zm1[c];
    cb = window->zp1[c];
    wvel[cb] = windat->wbot[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  memset(wmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->n3_e2; cc++) {
    c = window->w3_e2[cc];
    cs = window->m2d[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    ym1s = window->m2d[ym1];
    h1mid = 0.5 * (window->h2au1[cs] + window->h2au1[ym1s]);
    umid[c] = 0.5 * h1mid * h1mid * 
             (windat->u1[c] * Ds[cs] + windat->u1[ym1] * Ds[ym1s]);
    vmid[c] = 0.5 * area[cs] * (vel[c] + vel[yp1]);
    wmid[c] = 0.5 * (wvel[c] + wvel[ym1]);
  }
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    ym1 = window->ym1[c];
    wsur[c] = 0.5 * (wtop[c] + wtop[ym1]);
  }
  for(cc = 1; cc <= vcs; cc++) {
    c = wincon->i6[cc];
    cs = window->m2d[c];
    ym1 = window->ym1[cs];
    zm1 = window->zm1[c];
    wmid[c] = 0.5 * (windat->wbot[cs] + windat->wbot[ym1]);
    vel[zm1] = vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bmask = wincon->w4;
  bcond = (FILEIN | CUSTOM | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsiz * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2 + 1; cc <=  window->nbpte2; cc++) {
    c = window->bpte2[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor & bcond)) {
      for(cc = 1; cc <=  open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no3_e2,
	      open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
	      open[n]->cyc_e2, windat->nu2, windat->u2, 
	      windat->u2b, open[n]->bcond_nor,
	      windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, 0);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to3_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2,
	      windat->u2, windat->u2b, 
	      open[n]->bcond_tan, windat->dtf, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, 0);
  }
  if(wincon->momsc & EXPLICIT) {
    memset(bmask, 0, window->sgsiz * sizeof(double));
    nvel = windat->u2;
  }
  if(wincon->momsc & ADVECT_FORM) {
    dz = wincon->w5;
    for(cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      dz[c] = 1.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  /* Surface layers : nu1[zp1] = nu1[c] from no-slip condition.      */
  /* Note : over thin layers dz[c] != dzd hence multiplication by    */
  /* dz must be explicitly included in the flux terms. Below the     */
  /* surface this need not be done since dz[c] = dzd always.         */
  for (cc = 1; cc <= vcs; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu2[c], wincon->hmin);

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[ym1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midy[cs] * vel[c];
    if(bmask[ym1]) {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
              (vmid[c] * dz[yp1] * midy[yp1s] * vel[yp1] -
               vmid[ym1] * fluxc);

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = vmid[c] * dz[c] * midy[cs];
      d2 = vmid[ym1] * dz[ym1] * midy[ym1s];
      ime1 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
              (vmid[c] * dz[yp1] * midy[yp1s] * vel[yp1] -
               vmid[ym1] * (fluxc + dz[ym1] * midy[ym1s] * nvel[ym1]));

      /* Get the implicit multiplier for v.dv/dy                     */
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c] * dz[c] * midy[cs];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land or open      */
    /* boundaries nu2[xm1] = nu2[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[xm1]) {
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * dz[xp1] * vel[xp1] - umid[c] * fluxc);

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = umid[xp1] * dz[c];
      d2 = umid[c] * dz[xm1];
      ime2 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * dz[xp1] * vel[xp1] -
	       umid[c] * (fluxc + dz[xm1] * nvel[xm1]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u2c1[cs] * umid[xp1] * dz[c];
      exe2 = ime2 * vel[c];
    }

    /* Get the contribution from w.dv/dz                             */
    /*wu1 = -0.5 * (wsur[cs] * vel[c] - wmid[c] * vel[zm1]);*/
    /*wu1 = 0.5 * wmid[c] * vel[zm1] - wsur[cs] * vel[c];*/
    wu1 = 0.5 * wmid[c] * (vel[zm1] + vel[c]) - wsur[cs] * vel[c];

    /* Get the implicit multiplier for w.dv/dz                       */
    /*imz = 0.5 * (wmid[c] - wsur[cs]);*/
    /*imz = 0.5 * wmid[c];*/
    imz = 0.0;
    exz = 0.5 * wmid[c] * vel[c];
    if (wincon->momsc & WIMPLICIT) wu1 = imz = exz = 0.0;

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz) / dzd;
    adv = (ush1h2 + u1u2hs + wu1) / dzd;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += (exe1 + exe2 + exz) / dzd;
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[cs] * vel[c] * vel[c] * midy[cs] +
           wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy[cs];
    /* Note : nu2 is a copy of u2b, but may be adjusted for thin     */
    /* layers, hence use nu2 for u2b here.                           */
    velb = windat->nu2[c];
    windat->nu2[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.dv/dx + v.dv/dy) = (u2b-dtu.w.dv/dz) - nu2       */
    if(wincon->momsc & EXPLICIT) {
      adv -= ((exz + wu1) / dzd);
      wincon->u2inter[cs] += (windat->dtf * adv * windat->dzu2[c]);
    }
    else {
      immp = 1.0 - dtu * imz / dzd;
      adv = (velb + dtu * wu1 / dzd) / immp;
      wincon->u2inter[cs] += ((adv - windat->nu2[c]) * windat->dzu2[c] * 
			      windat->dtf / dtu);

      /* Fill the cells aove the free surface to the top of the grid */
      cs = c;
      while (c != zp1) {
	windat->nu2[zp1] = windat->nu2[cs];
	c = zp1;
	zp1 = window->zp1[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Sub - surface cells : nu2[zp1] is known                         */
  for (cc = vcs + 1; cc <= vc; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    xp1s = window->xp1[cs];
    xm1s = window->xm1[cs];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    yp1s = window->yp1[cs];
    ym1s = window->ym1[cs];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    dzd = max(windat->dzu2[c], wincon->hmin);

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[ym1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    fluxc = dz[c] * midy[cs] * vel[c];
    if(bmask[ym1]) {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
              (vmid[c] * dz[yp1] * midy[yp1s] * vel[yp1] -
               vmid[ym1] * fluxc) / dzd;

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = vmid[c];
      d2 = vmid[ym1] * dz[ym1] * midy[ym1s] / dzd;
      ime1 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[cs] * 
              (vmid[c] * dz[yp1] * midy[yp1s] * vel[yp1] -
               vmid[ym1] * (fluxc + dz[ym1] * midy[ym1s] * nvel[ym1])) / dzd;

      /* Get the implicit multiplier for v.dv/dy                     */
      ime1 = 0.5 * wincon->u2c1[cs] * vmid[c];
      exe1 = ime1 * fluxc / dzd;
    }

    /* Get the contribution from u.dv/dx. Note; at land or open      */
    /* boundaries nu2[xm1] = nu2[c] due to the no-slip condition.    */
    fluxc = dz[c] * vel[c];
    if(bmask[xm1]) {
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * dz[xp1] * vel[xp1] - umid[c] * fluxc) / dzd;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = umid[xp1];
      d2 = umid[c] * dz[xm1] / dzd;
      ime2 = 0.5 * wincon->u2c1[cs] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u2c1[cs] *
	      (umid[xp1] * dz[xp1] * vel[xp1] -
	       umid[c] * (fluxc + dz[xm1] * nvel[xm1])) / dzd;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u2c1[cs] * umid[xp1];
      exe2 = ime2 * fluxc / dzd;
    }

    /* Get the contribution from w.dv/dz                             */
    wu1 = -0.5 * (wmid[zp1] * (nvel[zp1] + vel[c]) -
		  wmid[c] * vel[zm1]) / dzd;

    /* Get the implicit multiplier for w.dv/dz                       */
    imz = 0.5 * wmid[c] / dzd;
    exz = imz * vel[c];

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2 + imz);
    adv = ush1h2 + u1u2hs + wu1;

    /* Explicit                                                      */
    if(wincon->momsc & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2 + exz;
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[cs] * vel[c] * vel[c] * midy[cs] +
           wincon->u2c4[cs] * u1au2[c] * u1au2[c] * midy[cs];
    velb = windat->nu2[c];
    windat->nu2[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.dv/dx + v.dv/dy) = (u2b-dtu.w.dv/dz) - nu2       */
    if(wincon->momsc & EXPLICIT) {
      adv -= (exz + wu1);
      wincon->u2inter[cs] += (windat->dtf * adv * windat->dzu2[c]);
    }
    else {
      immp = 1.0 - dtu * (imz);
      adv = (velb + dtu * wu1) / immp;
      wincon->u2inter[cs] += ((adv - windat->nu2[c]) * windat->dzu2[c] * 
			      windat->dtf / dtu);
    }
  }
}

/* END advect_u2_3d_ang_flux_f()                                     */
/*-------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U1AV ADVECTION ROUTINES                                           */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a backward direction. Formulated     */
/* using the advective form.                                         */
/*-------------------------------------------------------------------*/
void advect_u1_2d_ang_adv_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xmyp;
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h1au1;            /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midx;             /* Depth at the e1 face                  */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u2au1 = wincon->w7;
  h1au1 = window->h1au1;
  vel = windat->u1av;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  bmask = wincon->d3;
  nvel = windat->nu1av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    umid[c] = 0.25 * (2.0 * vel[c] + vel[xp1] + vel[xm1]) *
              window->h2au1[c];
    vmid[c] = windat->u2av[c] * window->h1au2[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1S + 1; cc <=  window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    bmask[c] = 1.0;
  }

  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e1,
              open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
              open[n]->cyc_e1, windat->nu1av, windat->u1av, 
	      windat->u1avb, open[n]->bcond_nor2d,
	      windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, U1BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	    open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	      windat->u1av, windat->u1avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, U1BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u1av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = wincon->vcs; cc >=1 ; cc--) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    xmyp = window->yp1[xm1];

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xp1]) {
      ush1h2 = 0.5 * wincon->u1c1[c] * umid[c] * 
	(vel[c] * h1au1[c] - vel[xm1] * h1au1[xm1]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h1au1[c];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = h1au1[xp1] * window->h2au1[c] / 
	  (midx[xp1] * window->h2au1[xp1]);
      ime1 = 0.5 * wincon->u1c1[c] * umid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[c] * umid[c] * (d1 * vel[c] - d2 * vel[xp1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[c] * umid[c] * 
	(vel[c] * h1au1[c] + nvel[xp1] * h1au1[xp1] -
	 vel[xm1] * h1au1[xm1]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[c] * umid[c] * h1au1[c];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[yp1] = nu1[c] due to the no-slip condition.    */
    if(bmask[yp1]) {
      u1u2hs = 0.25 * wincon->u1c1[c] * 
	(vmid[c] + vmid[xm1]) * (vel[c] * h1au1[c] - vel[ym1] * h1au1[ym1]);
				 
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u1c1[c] * 
	((vmid[yp1] + vmid[xmyp]) * nvel[yp1] * h1au1[yp1] +
	 (vmid[c] + vmid[xm1]) * (vel[c] * h1au1[c] - vel[ym1] * h1au1[ym1]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u1c1[c] * (vmid[yp1] + vmid[xmyp]) * h1au1[c];
      exe2 = ime2 * vel[c];
    }

    /* Sum the advective terms and implicit multipliers. Note the    */
    /* multipliers are added to 1.0 for backward stepping.           */
    immp = 1.0 + dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1.0;
      adv -= (exe1 + exe2);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[c] * vel[c] * vel[c] * midx[c] +
           wincon->u1c4[c] * u2au1[c] * u2au1[c] * midx[c];
    /* Note : nu1av is a copy of u1avb                               */
    velb = windat->nu1av[c];
    windat->nu1av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else {
      /*
      adv -= (exe1 + exe2);
      wincon->u1adv[c] += (adv * dtm);
      */
      wincon->u1adv[c] += ((windat->nu1av[c] - velb) * dtm / dtu);
    }
  }
}

/* END advect_u1_2d_ang_adv_b()                                (sbf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a forward direction. Formulated      */
/* using the advective form.                                         */
/*-------------------------------------------------------------------*/
void advect_u1_2d_ang_adv_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int xmyp;
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midx;             /* Depth at the e1 face                  */
  double *h1au1;            /* Cell area                             */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u2au1 = wincon->w7;
  h1au1 = window->h1au1;
  vel = windat->u1av;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  bmask = wincon->d3;
  nvel = windat->nu1av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    umid[c] = 0.25 * (2.0 * vel[c] + vel[xp1] + vel[xm1]) *
              window->h2au1[c];
    vmid[c] = windat->u2av[c] * window->h1au2[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1S + 1; cc <=  window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e1,
              open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
              open[n]->cyc_e1, windat->nu1av, windat->u1av, 
	      windat->u1avb, open[n]->bcond_nor2d,
               windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, U1BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	      windat->u1av, windat->u1avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, U1BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u1av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    xmyp = window->yp1[xm1];

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xm1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xm1]) {
      ush1h2 = 0.5 * wincon->u1c1[c] * umid[c] * 
	(vel[xp1] * h1au1[xp1] - vel[c] * h1au1[c]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h1au1[c];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = h1au1[xm1] * window->h2au1[c] / 
	  (midx[xm1] * window->h2au1[xm1]);
      ime1 = 0.5 * wincon->u1c1[c] * umid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[c] * umid[c] * (d1 * vel[c] - d2 * vel[xm1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[c] * umid[c] * 
	(vel[xp1] * h1au1[xp1] -
	 vel[c] * h1au1[c] + nvel[xm1] * h1au1[xm1]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[c] * umid[c] * h1au1[c];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[ym1] = nu1[c] due to the no-slip condition.    */
    if(bmask[ym1]) {
      u1u2hs = 0.25 * wincon->u1c1[c] * 
	(vmid[yp1] + vmid[xmyp]) *
	(vel[yp1] * h1au1[yp1] - vel[c] * h1au1[c]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u1c1[c] * 
	((vmid[yp1] + vmid[xmyp]) *
	 (vel[yp1] * h1au1[yp1] - vel[c] * h1au1[c]) +
	 (vmid[c] + vmid[xm1]) * nvel[ym1] * h1au1[ym1]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u1c1[c] * (vmid[c] + vmid[xm1]) * h1au1[c];
      exe2 = ime2 * vel[c];
    }

    /* Sum the advective terms and implicit multipliers. Note the    */
    /* multipliers are subtracted to 1.0 for backward stepping.      */
    immp = 1.0 - dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2;
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[c] * vel[c] * vel[c] * midx[c] +
           wincon->u1c4[c] * u2au1[c] * u2au1[c] * midx[c];
    /* Note : nu1av is a copy of u1avb                               */
    velb = windat->nu1av[c];
    windat->nu1av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else
      wincon->u1adv[c] += ((windat->nu1av[c] - velb) * dtm / dtu);
  }
}

/* END advect_u1_2d_ang_adv_f()                                (saf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a backward direction. Formulated     */
/* using the flux form.                                              */
/*-------------------------------------------------------------------*/
void advect_u1_2d_ang_flux_b(geometry_t *window,/* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midx;             /* Depth at the e1 face                  */
  double h2mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *Du1;              /* Dummy for e1 advective fluxes         */
  double *depth;            /* Depth of the water column             */
  double dth;               /* Depth divisor                         */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u2au1 = wincon->w7;
  area = window->cellarea;
  vel = windat->u1av;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  Du1 = wincon->d2;
  bmask = wincon->d3;
  if(!wincon->momsc2d & ADVECT_FORM)
    depth = windat->depth_e1;
  else
    depth = wincon->one;
  nvel = windat->nu1av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(Du1, 0, window->sgsizS * sizeof(double));
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    h2mid = 0.5 * (window->h1au2[c] + window->h1au2[xm1]);
    umid[c] = 0.5 * area[c] * (vel[c] + vel[xp1]);
    vmid[c] = 0.5 * h2mid * h2mid * 
              (windat->u2av[c] * Ds[c] + windat->u2av[xm1] * Ds[xm1]);
  }
  for (cc = 1; cc <= window->x2_e1; cc++) {
    c = window->w2_e1[cc];
    Du1[c] = depth[c] * vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1S + 1; cc <=  window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    bmask[c] = 1.0;
  }

  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e1,
              open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
              open[n]->cyc_e1, windat->nu1av, windat->u1av, 
	      windat->u1avb, open[n]->bcond_nor2d,
	      windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, U1BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	      windat->u1av, windat->u1avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, U1BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u1av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = wincon->vcs; cc >=1 ; cc--) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    dth = max(windat->depth_e1[c], wincon->hmin);

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xp1]) {
      ush1h2 = 0.5 * wincon->u1c1[c] * 
              (umid[c] * Du1[c] * midx[c] -
               umid[xm1] * Du1[xm1] * midx[xm1]) / dth;
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = umid[xm1];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = umid[c] * midx[c] * window->h2au1[c] / 
	  (midx[xp1] * window->h2au1[xp1]);
      ime1 = 0.5 * wincon->u1c1[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u1c1[c] * (umid[xm1] * Du1[c] * midx[c] -
				      umid[c] * Du1[xp1] * midx[xp1]) / dth;
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[c] * 
	      (umid[c] * (Du1[c] * midx[c] + 
			  depth[xp1] * midx[xp1] * nvel[xp1]) -
	       umid[xm1] * (Du1[c] * midx[c] + Du1[xm1] * midx[xm1])) / dth;

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[c] * umid[xm1];
      exe1 = ime1 * Du1[c] * midx[c] / dth;
      ime1 = exe1 = 0.0;
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[yp1] = nu1[c] due to the no-slip condition.    */
    if(bmask[yp1]) {
      u1u2hs = 0.5 * wincon->u1c1[c] *
	      (vmid[yp1] * Du1[c] - vmid[c] * Du1[ym1]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = vmid[c];
      d2 = vmid[yp1] * depth[yp1] / dth;
      ime2 = 0.5 * wincon->u1c1[c] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u1c1[c] *
	      (vmid[yp1] * (Du1[c] + depth[yp1] * nvel[yp1]) -
	       vmid[c] * Du1[ym1]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u1c1[c] * vmid[c];
      exe2 = ime2 * Du1[c] / dth;
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

      immp = 1.0 + dtu * (ime1);
      adv -= (exe2);

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1.0;
      adv -= (exe1 + exe2);
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[c] * vel[c] * vel[c] * midx[c] +
           wincon->u1c4[c] * u2au1[c] * u2au1[c] * midx[c];
    /* Note : nu1av is a copy of u1avb                               */
    velb = windat->nu1av[c];
    windat->nu1av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else {
      adv -= (exe1);
      wincon->u1adv[c] += (adv * dtm);
    /*wincon->u1adv[c] += ((windat->nu1av[c] - velb) * dtm / dtu);*/
    }
  }
}

/* END advect_u1_2d_ang_flux_b()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the e1 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a forward direction. Formulated      */
/* using the flux form.                                              */
/*-------------------------------------------------------------------*/
void advect_u1_2d_ang_flux_f(geometry_t *window,/* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u2au1;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midx;             /* Depth at the e1 face                  */
  double h2mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *Du1;              /* Dummy for e1 advective fluxes         */
  double *depth;            /* Depth of the water column             */
  double dth;               /* Depth divisor                         */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u1av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u2au1 = wincon->w7;
  area = window->cellarea;
  vel = windat->u1av;
  midx = wincon->mdx;
  Ds = wincon->Ds;
  Du1 = wincon->d2;
  bmask = wincon->d3;
  if(!wincon->momsc2d & ADVECT_FORM)
    depth = windat->depth_e1;
  else
    depth = wincon->one;
  nvel = windat->nu1av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(Du1, 0, window->sgsizS * sizeof(double));
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    h2mid = 0.5 * (window->h1au2[c] + window->h1au2[xm1]);
    umid[c] = 0.5 * area[c] * (vel[c] + vel[xp1]);
    vmid[c] = 0.5 * h2mid * h2mid * 
              (windat->u2av[c] * Ds[c] + windat->u2av[xm1] * Ds[xm1]);
  }
  for (cc = 1; cc <= window->x2_e1; cc++) {
    c = window->w2_e1[cc];
    Du1[c] = depth[c] * vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe1S + 1; cc <=  window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e1,
              open[n]->obc_e1, open[n]->oi1_e1, open[n]->oi2_e1,
              open[n]->cyc_e1, windat->nu1av, windat->u1av, 
	      windat->u1avb, open[n]->bcond_nor2d,
	      windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_nor, U1BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
	c = open[n]->obc_e1[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e1 + 1,
	      open[n]->to2_e1, open[n]->obc_e1, open[n]->oi1_e1,
	      open[n]->oi2_e1, open[n]->cyc_e1, windat->nu1av,
	      windat->u1av, windat->u1avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau1, 
	      open[n]->transfer_u1, open[n]->relax_zone_tan, U1BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u1av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    dth = max(windat->depth_e1[c], wincon->hmin);

    /* Get the contribution from u.du/dx. Note; at land boundaries   */
    /* nu1[xm1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[xm1]) {
      ush1h2 = 0.5 * wincon->u1c1[c] * 
              (umid[c] * Du1[xp1] * midx[xp1] -
               umid[xm1] * Du1[c] * midx[c]) / dth;

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = umid[c];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = umid[xm1] * midx[c] * window->h2au1[c] / 
	  (midx[xm1] * window->h2au1[xm1]);
      ime1 = 0.5 * wincon->u1c1[c] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u1c1[c] * 
              (umid[c] * Du1[xp1] * midx[xp1] -
               umid[xm1] * (Du1[c] * midx[c] +
			    depth[xm1] * midx[xm1] * nvel[xm1])) / dth;

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u1c1[c] * umid[c];
      exe1 = ime1 * Du1[c] * midx[c] / dth;
    }

    /* Get the contribution from v.du/dy. Note; at land or open      */
    /* boundaries nu1[ym1] = nu1[c] due to the no-slip condition.    */
    if(bmask[ym1]) {

      u1u2hs = 0.5 * wincon->u1c1[c] *
	      (vmid[yp1] * Du1[yp1] - vmid[c] * Du1[c]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = vmid[yp1];
      d2 = vmid[c] * depth[ym1] / dth;
      ime2 = 0.5 * wincon->u1c1[c] * (d1 - d2);
    }
    else { 

      u1u2hs = 0.5 * wincon->u1c1[c] *
	      (vmid[yp1] * Du1[yp1] -
	       vmid[c] * (depth[ym1] * nvel[ym1] + Du1[c])) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u1c1[c] * vmid[yp1];
      exe2 = ime2 * Du1[c] / dth;
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2;
    }

    /* Update the velocity                                           */
    adv += wincon->u1c3[c] * vel[c] * vel[c] * midx[c] +
           wincon->u1c4[c] * u2au1[c] * u2au1[c] * midx[c];
    /* Note : nu1av is a copy of u1avb                               */
    velb = windat->nu1av[c];
    windat->nu1av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else
      wincon->u1adv[c] += ((windat->nu1av[c] - velb) * dtm / dtu);
  }
}

/* END advect_u1_2d_ang_flux_f()                                     */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* U2AV ADVECTION ROUTINES                                           */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a backward direction. Formulated     */
/* using the advective form.                                         */
/*-------------------------------------------------------------------*/
void advect_u2_2d_ang_adv_b(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int ymxp;
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h2au2;            /* Cell width at e2 location             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midy;             /* Depth at the e1 face                  */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u1au2 = wincon->w7;
  h2au2 = window->h2au2;
  vel = windat->u2av;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  bmask = wincon->d3;
  nvel = windat->nu2av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    vmid[c] = 0.25 * (2.0 * vel[c] + vel[yp1] + vel[ym1]) *
              window->h1au2[c];
    umid[c] = windat->u1av[c] * window->h2au1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2S + 1; cc <=  window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e2,
              open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
              open[n]->cyc_e2, windat->nu2av, windat->u2av, 
	      windat->u2avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, U2BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to2_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2av,
	      windat->u2av, windat->u2avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, U2BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u2av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = wincon->vcs; cc >= 1; cc--) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    ymxp = window->xp1[ym1];

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[yp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[yp1]) {
      ush1h2 = 0.5 * wincon->u2c1[c] * vmid[c] * 
	      (vel[c] * h2au2[c] - vel[ym1] * h2au2[ym1]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h2au2[c];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = h2au2[yp1] * window->h1au2[c] / 
	  (midy[yp1] * window->h1au2[yp1]);
      ime1 = 0.5 * wincon->u2c1[c] * vmid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u2c1[c] * vmid[c] * (d1 * vel[c] - d2 * vel[yp1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[c] * vmid[c] * 
	(vel[c] * h2au2[c] + nvel[yp1] * h2au2[yp1] -
	 vel[ym1] * h2au2[ym1]);

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u2c1[c] * vmid[c] * h2au2[c];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land boundaries   */
    /* nu2[xp1] = nu2[c] due to the no-slip condition.               */
    if(bmask[xp1]) {
      u1u2hs = 0.25 * wincon->u2c1[c] * 
	(umid[c] + umid[ym1]) * (vel[c] * h2au2[c] - vel[xm1] * h2au2[xm1]);
					   
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u2c1[c] * 
	((umid[xp1] + umid[ymxp]) * nvel[xp1] * h2au2[xp1] +
	 (umid[c] + umid[ym1]) *
	 (vel[c] * h2au2[c] - vel[xm1] * h2au2[xm1]));

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u2c1[c] * (umid[xp1] + umid[ymxp]) * h2au2[c];
      exe2 = ime2 * vel[c];
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv -= (exe1 + exe2);
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[c] * vel[c] * vel[c] * midy[c] +
           wincon->u2c4[c] * u1au2[c] * u1au2[c] * midy[c];
    /* Note : nu2av is a copy of u2avb                               */
    velb = windat->nu2av[c];
    windat->nu2av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u2adv[c] += (adv * dtm);
    else {
      /*
      adv -= (exe1 + exe2);
      wincon->u2adv[c] += (adv * dtm);
      */
      wincon->u2adv[c] += ((windat->nu2av[c] - velb) * dtm / dtu);
    }
  }
}

/* END advect_u2_2d_ang_adv_b()                                (sbf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a forwward direction. Formulated     */
/* using the advective form.                                         */
/*-------------------------------------------------------------------*/
void advect_u2_2d_ang_adv_f(geometry_t *window, /* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  int ymxp;
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *h2au2;            /* Cell width at e2 location             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midy;             /* Depth at the e1 face                  */
  double *Ds;               /* Total depth for sigma                 */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u1au2 = wincon->w7;
  h2au2 = window->h2au2;
  vel = windat->u2av;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  bmask = wincon->d3;
  nvel = windat->nu2av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    vmid[c] = 0.25 * (2.0 * vel[c] + vel[yp1] + vel[ym1]) *
              window->h1au2[c];
    umid[c] = windat->u1av[c] * window->h2au1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2S + 1; cc <=  window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e2,
              open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
              open[n]->cyc_e2, windat->nu2av, windat->u2av, 
	      windat->u2avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, U2BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to2_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2av,
	      windat->u2av, windat->u2avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, U2BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u2av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    ymxp = window->xp1[ym1];

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[ym1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[ym1]) {
      ush1h2 = 0.5 * wincon->u2c1[c] * vmid[c] * 
	      (vel[yp1] * h2au2[yp1] - vel[c] * h2au2[c]);
	       
      /* Get the implicit multiplier for u.du/dx                     */
      d1 = h2au2[c];
      /* See bdry_u1_2d() for 2D NOGRAD formulation                  */
      d2 = h2au2[yp1] * window->h1au2[c] / 
	  (midy[yp1] * window->h1au2[yp1]);
      ime1 = 0.5 * wincon->u2c1[c] * vmid[c] * (d1 - d2);
      exe1 = 0.5 * wincon->u2c1[c] * vmid[c] * (d1 * vel[c] - d2 * vel[ym1]);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[c] * vmid[c] * 
	(vel[yp1] * h2au2[yp1] - 
	 (vel[c] * h2au2[c] + nvel[ym1] * h2au2[ym1]));

      /* Get the implicit multiplier for u.du/dx                     */
      ime1 = 0.5 * wincon->u2c1[c] * vmid[c] * h2au2[c];
      exe1 = ime1 * vel[c];
    }

    /* Get the contribution from u.dv/dx. Note; at land boundaries   */
    /* nu2[xm1] = nu2[c] due to the no-slip condition.               */
    if(bmask[xm1]) {
      u1u2hs = 0.25 * wincon->u2c1[c] * 
	(umid[xp1] + umid[ymxp]) * (vel[xp1] * h2au2[xp1] - vel[c] * h2au2[c]);
					   
      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.0;
    }
    else { 
      u1u2hs = 0.25 * wincon->u2c1[c] * 
	((umid[xp1] + umid[ymxp]) * 
	 (vel[xp1] * h2au2[xp1] - vel[c] * h2au2[c]) -
	 (umid[c] + umid[ym1]) * nvel[xm1] * h2au2[xm1]);

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.25 * wincon->u2c1[c] * (umid[c] + umid[ym1]) * h2au2[c];
      exe2 = ime2 * vel[c];
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2;
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[c] * vel[c] * vel[c] * midy[c] +
           wincon->u2c4[c] * u1au2[c] * u1au2[c] * midy[c];
    /* Note : nu2av is a copy of u2avb                               */
    velb = windat->nu2av[c];
    windat->nu2av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else
      wincon->u2adv[c] += ((windat->nu2av[c] - velb) * dtm / dtu);
  }
}

/* END advect_u2_2d_ang_adv_f()                                (saf) */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative       */
/* stepping through the grid in a forwward direction. Formulated     */
/* using the flux form.                                              */
/*-------------------------------------------------------------------*/
void advect_u2_2d_ang_flux_f(geometry_t *window,/* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midy;             /* Depth at the e1 face                  */
  double h1mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *Du2;              /* Dummy for e1 advective fluxes         */
  double *depth;            /* Depth of the water column             */
  double dth;               /* Depth divisor                         */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u1au2 = wincon->w7;
  area = window->cellarea;
  vel = windat->u2av;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  Du2 = wincon->d2;
  bmask = wincon->d3;
  if(!wincon->momsc2d & ADVECT_FORM)
    depth = windat->depth_e2;
  else
    depth = wincon->one;
  nvel = windat->nu2av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(Du2, 0, window->sgsizS * sizeof(double));
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    h1mid = 0.5 * (window->h2au1[c] + window->h2au1[ym1]);
    vmid[c] = 0.5 * area[c] * (vel[c] + vel[yp1]);
    umid[c] = 0.5 * h1mid * h1mid * 
             (windat->u1av[c] * Ds[c] + windat->u1av[ym1] * Ds[ym1]);
  }
  for (cc = 1; cc <= window->x2_e2; cc++) {
    c = window->w2_e2[cc];
    Du2[c] = depth[c] * vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2S + 1; cc <=  window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & R_EDGE || open[n]->ocodey & F_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e2,
              open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
              open[n]->cyc_e2, windat->nu2av, windat->u2av, 
	      windat->u2avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, U2BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to2_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2av,
	      windat->u2av, windat->u2avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, U2BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u2av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    dth = max(windat->depth_e2[c], wincon->hmin);

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[ym1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[ym1]) {
      ush1h2 = 0.5 * wincon->u2c1[c] * 
              (vmid[c] * Du2[yp1] * midy[yp1] -
               vmid[ym1] * Du2[c] * midy[c]) / dth;

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = vmid[c];
      /* See bdry_u2_2d() for 2D NOGRAD formulation                  */
      d2 = vmid[ym1] * midy[c] * window->h1au2[c]/ 
 	  (midy[ym1] * window->h1au2[ym1]);
      ime1 = 0.5 * wincon->u2c1[c] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[c] * 
              (vmid[c] * Du2[yp1] * midy[yp1] -
               vmid[ym1] * (Du2[c] * midy[c] +
			    depth[ym1] * midy[ym1] * nvel[ym1])) / dth;

      /* Get the implicit multiplier for u.du/dx                       */
      ime1 = 0.5 * wincon->u2c1[c] * vmid[c];
      exe1 = ime1 * Du2[c] * midy[c] / dth;
    }

    /* Get the contribution from u.dv/dx. Note; at land boundaries   */
    /* nu2[xm1] = nu2[c] due to the no-slip condition.               */
    if(bmask[xm1]) {

      u1u2hs = 0.5 * wincon->u2c1[c] *
	      (umid[xp1] * Du2[xp1] - umid[c] * Du2[c]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = umid[xp1];
      d2 = umid[c] * depth[xm1] / dth;
      ime2 = 0.5 * wincon->u2c1[c] * (d1 - d2);
    }
    else { 

      u1u2hs = 0.5 * wincon->u2c1[c] *
	      (umid[xp1] * Du2[xp1] -
	       umid[c] * (depth[xm1] * nvel[xm1] + Du2[c])) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u2c1[c] * umid[xp1];
      exe2 = ime2 * Du2[c] / dth;
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 - dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv += exe1 + exe2;
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[c] * vel[c] * vel[c] * midy[c] +
           wincon->u2c4[c] * u1au2[c] * u1au2[c] * midy[c];
    /* Note : nu2av is a copy of u2avb                               */
    velb = windat->nu2av[c];
    windat->nu2av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u1adv[c] += (adv * dtm);
    else
      wincon->u2adv[c] += ((windat->nu2av[c] - velb) * dtm / dtu);
  }
}

/* END advect_u2_2d_ang_flux_f()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the e2 advective fluxes for the 2D mode      */
/* using an unconditionally stable angular centered derivative.      */
/*-------------------------------------------------------------------*/
void advect_u2_2d_ang_flux_b(geometry_t *window,/* Processing window */
			    window_t *windat,   /* Window data       */
			    win_priv_t *wincon  /* Window constants  */
  )
{
  int c, cc, n;             /* Sparse coordinate / counter           */
  int xp1, xm1;             /* 3D sparse coordinate at i+1, i-1      */
  int yp1, ym1;             /* 3D sparse coordinate at j+1, j-1      */
  double dtu;               /* Leapfrog sub-time step to use         */
  double dtm;               /* Minimum forward time step             */
  double ime1;              /* Implicit multiplier                   */
  double ime2;              /* Implicit multiplier                   */
  double immp;              /* Implicit multiplier                   */
  double exe1;              /* Explicit multiplier                   */
  double exe2;              /* Explicit multiplier                   */
  double adv;               /* Sum of advective terms                */
  double *u1au2;            /* u2 value at the e1 face               */
  double *vel;              /* Velocity to use in spatial gradients  */
  double velb;              /* Velocity at previous timestep         */
  double *nvel;             /* Updated velocity                      */
  double *area;             /* Cell area                             */
  double ush1h2;            /* u1 * u1 * h1 * h2                     */
  double u1u2hs;            /* u1 * u2 * h1 * h1                     */
  double *umid;             /* Cell centered u velocity              */
  double *vmid;             /* Corner centered v velocity            */
  double *midy;             /* Depth at the e1 face                  */
  double h1mid;             /* Corner centered h1au2                 */
  double *Ds;               /* Total depth for sigma                 */
  double *Du2;              /* Dummy for e1 advective fluxes         */
  double *depth;            /* Depth of the water column             */
  double dth;               /* Depth divisor                         */
  double *bmask;            /* Open boundary mask                    */
  double d1, d2;            /* Dummies                               */
  int bcond;                /* Normal boundary condition             */
  open_bdrys_t **open = window->open;

  if (!wincon->nonlinear || wincon->u2av_f & ADVECT)
    return;

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  umid = wincon->w1;
  vmid = wincon->w2;
  u1au2 = wincon->w7;
  area = window->cellarea;
  vel = windat->u2av;
  midy = wincon->mdy;
  Ds = wincon->Ds;
  Du2 = wincon->d2;
  bmask = wincon->d3;
  if(!wincon->momsc2d & ADVECT_FORM)
    depth = windat->depth_e2;
  else
    depth = wincon->one;
  nvel = windat->nu2av;
  exe1 = exe2 = 0.0;

  /*-----------------------------------------------------------------*/
  /* Get the sub-timestep                                            */
  d1 = windat->dt / windat->dtf2;   /* iratio for this step          */
  if (wincon->stab & (MULTI_DT))
    dtu = windat->dt2s;
  else
    dtu = windat->dtu1 / d1;        /* 2D forward time-step          */
  dtm = dtu;
  dtu += windat->dtb2;  /* 2D leapfrog time-step; no sub-stepping    */

  /*-----------------------------------------------------------------*/
  /* Get the velocities at the cell faces / corners                  */
  memset(Du2, 0, window->sgsizS * sizeof(double));
  memset(umid, 0, window->sgsiz * sizeof(double));
  memset(vmid, 0, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    h1mid = 0.5 * (window->h2au1[c] + window->h2au1[ym1]);
    vmid[c] = 0.5 * area[c] * (vel[c] + vel[yp1]);
    umid[c] = 0.5 * h1mid * h1mid * 
             (windat->u1av[c] * Ds[c] + windat->u1av[ym1] * Ds[ym1]);
  }
  for (cc = 1; cc <= window->x2_e2; cc++) {
    c = window->w2_e2[cc];
    Du2[c] = depth[c] * vel[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary mask                                           */
  bcond = (FILEIN | CUSTOM | VERTIN | CLAMPD | CYCLIC);
  memset(bmask, 0, window->sgsizS * sizeof(double));
  /* Set the tangential land boundaries                              */
  for(cc = window->nbe2S + 1; cc <=  window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    bmask[c] = 1.0;
  }
  /* Set the open boundaries                                         */
  for (n = 0; n < window->nobc; n++) {

    if (open[n]->ocodex & L_EDGE || open[n]->ocodey & B_EDGE)
      continue;

    /* Set the normal open boundaries if !bcond_nor&bcond. These     */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_nor2d & bcond)) {
      for(cc = 1; cc <=  open[n]->no2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      /* Update the normal velocity on the boundary                  */
      set_OBC(window, windat, wincon, open[n], 1, open[n]->no2_e2,
              open[n]->obc_e2, open[n]->oi1_e2, open[n]->oi2_e2,
              open[n]->cyc_e2, windat->nu2av, windat->u2av, 
	      windat->u2avb, open[n]->bcond_nor2d,
              windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_nor, U2BDRY);

    /* Set the tangential open boundaries if !bcond_tan&bcond. These */
    /* boundaries are then treated with a no-gradient condition.     */
    if(!(open[n]->bcond_tan2d & bcond)) {
      for(cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
	c = open[n]->obc_e2[cc];
	bmask[c] = 1.0;
      }
    }
    else
      set_OBC(window, windat, wincon, open[n], open[n]->no3_e2 + 1,
	      open[n]->to2_e2, open[n]->obc_e2, open[n]->oi1_e2,
	      open[n]->oi2_e2, open[n]->cyc_e2, windat->nu2av,
	      windat->u2av, windat->u2avb, 
	      open[n]->bcond_tan2d, windat->dtf2, &open[n]->datau2, 
	      open[n]->transfer_u2, open[n]->relax_zone_tan, U2BDRY);
  }
  if(wincon->momsc2d & EXPLICIT) {
    memset(bmask, 0, window->sgsizS * sizeof(double));
    nvel = windat->u2av;
  }

  /*-----------------------------------------------------------------*/
  /* Do the integration                                              */
  for (cc = wincon->vcs; cc >= 1; cc--) {
    c = wincon->s3[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    dth = max(windat->depth_e2[c], wincon->hmin);

    /* Get the contribution from v.dv/dy. Note; at land boundaries   */
    /* nu2[yp1] = 0 due to the zero flux condition.                  */
    /* Set a no-gradient condition over open boundaries              */
    if(bmask[yp1]==-1) {
      ush1h2 = 0.5 * wincon->u2c1[c] * 
              (vmid[c] * Du2[c] * midy[c] -
               vmid[ym1] * Du2[ym1] * midy[ym1]) / dth;

      /* Get the implicit multiplier for u.du/dx                     */
      d1 = vmid[ym1];
      /* See bdry_u2_2d() for 2D NOGRAD formulation                  */
      d2 = vmid[c] * midy[c] * window->h1au2[c]/ 
 	   (midy[yp1] * window->h1au2[yp1]);
      ime1 = 0.5 * wincon->u2c1[c] * (d1 - d2);
    }
    else {
      ush1h2 = 0.5 * wincon->u2c1[c] * 
              (vmid[c] * (Du2[c] * midy[c] +
			  depth[yp1] * midy[yp1] * nvel[yp1]) -
               vmid[ym1] * Du2[ym1] * midy[ym1]) / dth;

      /* Get the implicit multiplier for u.du/dx                       */
      ime1 = 0.5 * wincon->u2c1[c] * vmid[ym1];
      exe1 = ime1 * Du2[c] * midy[c] / dth;
    }

    /* Get the contribution from u.dv/dx. Note; at land boundaries   */
    /* nu2[xp1] = nu2[c] due to the no-slip condition.               */
    if(bmask[xp1]) {
      u1u2hs = 0.5 * wincon->u2c1[c] *
	      (umid[xp1] * Du2[c] - umid[c] * Du2[xm1]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      d1 = umid[c];
      d2 = umid[xp1] * depth[xp1] / dth;
      ime2 = 0.5 * wincon->u2c1[c] * (d1 - d2);
    }
    else { 
      u1u2hs = 0.5 * wincon->u2c1[c] *
	      (umid[xp1] * (Du2[c] + depth[xp1] * nvel[xp1]) -
	       umid[c] * Du2[xm1]) / dth;

      /* Get the implicit multiplier for v.du/dy                     */
      ime2 = 0.5 * wincon->u2c1[c] * umid[c];
      exe2 = ime2 * Du2[c] / dth;
    }

    /* Sum the advective terms and implicit multipliers              */
    immp = 1.0 + dtu * (ime1 + ime2);
    adv = ush1h2 + u1u2hs;

      immp = 1.0 + dtu * (ime2);
      adv -= (exe1);

    /* Explicit                                                      */
    if(wincon->momsc2d & EXPLICIT) {
      immp = 1;
      adv -= (exe1 + exe2);
    }

    /* Update the velocity                                           */
    adv += wincon->u2c3[c] * vel[c] * vel[c] * midy[c] +
           wincon->u2c4[c] * u1au2[c] * u1au2[c] * midy[c];
    /* Note : nu2av is a copy of u2avb                               */
    velb = windat->nu2av[c];
    windat->nu2av[c] = (velb + dtu * adv) / immp;

    /* Integrate the e1 dispersion term vertically and with time     */
    /* (i.e. over the forward time-step; windat->dtf).               */
    /* Note : dtu(u.du/dx + v.du/dy) = u1avb - nu1av                 */
    if(wincon->momsc2d & EXPLICIT)
      wincon->u2adv[c] += (adv * dtm);
    else {
      /*wincon->u2adv[c] += ((windat->nu2av[c] - velb) * dtm / dtu);*/
      adv -= (exe2);
      wincon->u2adv[c] += (adv * dtm);
    }
  }
}

/* END advect_u2_2d_ang_flux_b()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-Lagrangian advection for u1 component of velocity.           */
/*-------------------------------------------------------------------*/
void semi_lagrange_u1(geometry_t *window,  /* Processing window      */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int *cells,          /* Cells to process       */
		      int vcs,             /* No. 2D cells           */
		      int vc               /* No. 3D cells           */
  )
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->w1;
  double *nv = wincon->w2;
  double *nw = wincon->w3;
  double *cx = wincon->w4;
  double *cy = wincon->w5;
  double *cz = wincon->w6;
  double *dzz = wincon->w8;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  double u, v, w;
  double p, q, r, pp, qq, he1, he2;
  int cc, c, ci, c2, zm1;
  int *cl = wincon->s5;        /* Streamline origin */
  int *c2cc = wincon->s4;      /* c to cc map */
  int nc, n;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(cx, 0, window->sgsiz * sizeof(double));
  memset(cy, 0, window->sgsiz * sizeof(double));
  memset(cz, 0, window->sgsiz * sizeof(double));
  memcpy(nu, windat->u1, window->sgsiz * sizeof(double));
  memcpy(nv, wincon->w7, window->sgsiz * sizeof(double));
  memset(nw, 0, window->sgsiz * sizeof(double));
  memset(cl, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing */
  set_dzzv(window, dzz, cells, vc, vcs, U13D);
  for (c = 1; c <= window->enonS; c++) {
    h1[c] = window->h1acell[window->xm1[c]];
    h2[c] = 0.5 * (window->h2au2[c] + window->h2au2[window->xm1[c]]);
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face */
  vel_face(window, windat, wincon, nu, nv, nw, cells, vcs, vc, U13D);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    nu[c] = 0.0;
    nv[c] = nv[ci];
    windat->u1b[c] = 0.0;
  }
  for (cc = 1; cc <= window->nbe1; cc++) {
    c = window->bpte1[cc];
    ci = window->bine1[cc];
    nu[c] = 0.0;
    windat->u1b[c] = 0.0;
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }
  /* Set the tangential velocity boundary conditions (slip)          */
  for (; cc <= window->nbpte1; cc++) {
    c = window->bpte1[cc];
    ci = window->bine1[cc];
    nu[c] = wincon->slip * nu[ci];
    nv[c] = -nv[ci];
    windat->u1b[c] = wincon->slip * windat->u1b[ci];
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    zm1 = window->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
    windat->u1b[zm1] = -windat->u1b[c];
  }
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    ci = window->zp1[c];
    nu[c] = nu[ci];
    nv[c] = nv[ci];
    nw[c] = nw[ci];
    windat->u1b[c] = windat->u1b[ci];
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  for (cc = 1; cc <= vc; cc++) {
    c = cg = cl[cc] = cells[cc];
    time_left = windat->dt;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    pp = qq = 0.0;
    n=0;
    while (time_left > 0) {
      int ins = (cl[cc] == wincon->m2d[cl[cc]]) ? 1 : 0;
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      c2 = window->m2d[cl[cc]];
      p = (nu[cl[cc]] > 0.0) ? h1[c2] : h1[window->xp1[c2]];
      q = (nv[cl[cc]] > 0.0) ? h2[c2] : h2[window->yp1[c2]];
      r = (nw[cl[cc]] < 0.0) ? dzz[cl[cc]] : dzz[window->zm1[cl[cc]]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = (r) ? ((ins && nw[cl[cc]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[cc]]))) : dt;
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position.  */
      cl[cc] =
        get_posv(window, cl[cc], u, v, w, &cx[c], &cy[c], &cz[c], dt, &nc, U1GEN, &pp, &qq);
      if (nc > 2)hd_warn("Streamline distance > 2 cells @ (%d %d) : %4.2f days\n",
			 window->s2i[c], window->s2j[c], windat->days);

      c2 = window->m2d[cl[cc]];
      p = cx[c] / geth(h1, qq, c2, window->ym1);
      q = cy[c] / geth(h2, pp, c2, window->xm1);
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;

      /* Since the velocity field may be non-conservative, it is     */
      /* possible that the streamline may arrive above the surface   */
      /* or below the sea bed, hence if the origin is in the top or  */
      /* bottom layer, the distance should be limited to eta-cellz   */
      /* or cellz-botz (i.e. dzz). Limiting 0<r<1 is equivalent to   */
      /* this.                                                       */
      /* This is currently controlled with a limit on *cz in         */
      /* get_pos().                                                  */
      /*r = min(max(0.0, r), 1.0);*/

      wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - r) * (1.0 - p) * q;
      wgt[2] = (1.0 - r) * p * (1.0 - q);
      wgt[3] = (1.0 - r) * p * q;
      wgt[4] = r * (1.0 - p) * (1.0 - q);
      wgt[5] = r * (1.0 - p) * q;
      wgt[6] = r * p * (1.0 - q);
      wgt[7] = r * p * q;

      /* Interpolate x,y,z velocities */
      u = int_val(window, cl[cc], wgt, nu);
      v = int_val(window, cl[cc], wgt, nv);
      w = int_val(window, cl[cc], wgt, nw);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      n++;

      /*
      if(p > 1.0 || p < 0.0)
	hd_warn("x cell distance outside of bounds at %d : %f\n", c, p);
      if(q > 1.0 || q < 0.0)
	hd_warn("y cell distance outside of bounds at %d : %f\n", c, q);
      if(r > 1.0 || r < 0.0)
	hd_warn("r cell distance outside of bounds at %d : %f\n", c, r);
      */
    }

    sl_check(window, &cl[cc], &cx[c], &cy[c], &cz[c]);
    c2 = window->m2d[cl[cc]];
    he1 = geth(h1, qq, c2, window->ym1);;
    he2 = geth(h2, pp, c2, window->xm1);

    if (wincon->osl == 0) {

      p = cx[c] / he1;
      q = cy[c] / he2;
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;
      /*r = min(max(0.0, r), 1.0);*/

      /* Get the interpolation weights */
      wincon->wgt[c][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][1] = (1.0 - r) * (1.0 - p) * q;
      wincon->wgt[c][2] = (1.0 - r) * p * (1.0 - q);
      wincon->wgt[c][3] = (1.0 - r) * p * q;
      wincon->wgt[c][4] = r * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][5] = r * (1.0 - p) * q;
      wincon->wgt[c][6] = r * p * (1.0 - q);
      wincon->wgt[c][7] = r * p * q;
    } else {
      if (wincon->osl == 2 || wincon->osl == 4) {
	set_eov(window, wincon, he1, he2, dzz, &cl[cc], &cx[c], &cy[c], &cz[c]);
	weights_eov(window, c, cl[cc], cx[c], cy[c], cz[c], h1, h2, dzz);
      } else
	weights_oov(window, c, cl[cc], cx[c], cy[c], cz[c], h1, h2, dzz);
    }
  }

  /* Do the final interpolation */
  memcpy(windat->nu1, windat->u1b, window->sgsiz * sizeof(double));
  if (wincon->osl == 0) {
    if(wincon->togn & TOPRIGHT) {
      for (cc = 1; cc <= vc; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	/*check_tridiagonal(window,windat->u1b,"before",c,cl[cc]);*/
	windat->nu1[c] = int_val(window, cl[cc], wincon->wgt[c], windat->u1b);
	/*
	if(isnan(windat->nu1[c])) {
	  print_tridiagonal(window,windat->u1b,c,cl[cc]);
	  s2ijk(window,c);
	  hd_quit("NaN at %f days\n",windat->days);
	}
	*/
      }
    } else {
      for (cc = 1; cc <= vc; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	windat->nu1[c] = int_val_bl(window, cl[cc], wincon->wgt[c], windat->u1b);
      }
    }
  } else {
    for (cc = 1; cc <= vc; cc++) {
      c = cg = cells[cc];
      c2cc[c] = cc;
      windat->nu1[c] = int_valo(window, wincon->wgt[c], wincon->lmap[cl[cc]], 
				windat->u1b, wincon->osl);
    }
  }
}

/* END semi_lagrange_u1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-Lagrangian advection for u2 component of velocity.           */
/*-------------------------------------------------------------------*/
void semi_lagrange_u2(geometry_t *window,  /* Processing window      */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int *cells,          /* Cells to process       */
		      int vcs,             /* No. 2D cells           */
		      int vc               /* No. 3D cells           */
  )
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->w1;
  double *nv = wincon->w2;
  double *nw = wincon->w3;
  double *cx = wincon->w4;
  double *cy = wincon->w5;
  double *cz = wincon->w6;
  double *dzz = wincon->w8;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  double u, v, w;
  double p, q, r, pp, qq, he1, he2;
  int cc, c, ci, c2, zm1;
  int *cl = wincon->s5;        /* Streamline origin */
  int *c2cc = wincon->s4;      /* c to cc map */
  int nc, n;

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(cx, 0, window->sgsiz * sizeof(double));
  memset(cy, 0, window->sgsiz * sizeof(double));
  memset(cz, 0, window->sgsiz * sizeof(double));
  memcpy(nu, wincon->w7, window->sgsiz * sizeof(double));
  memcpy(nv, windat->u2, window->sgsiz * sizeof(double));
  memset(nw, 0, window->sgsiz * sizeof(double));
  memset(cl, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the vertical cell spacing */
  set_dzzv(window, dzz, cells, vc, vcs, U23D);
  for (c = 1; c <= window->enonS; c++) {
    h1[c] = 0.5 * (window->h1au1[c] + window->h1au1[window->ym1[c]]);
    h2[c] = window->h2acell[window->ym1[c]];
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face */
  vel_face(window, windat, wincon, nu, nv, nw, cells, vcs, vc, U23D);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    ci = window->bine2[cc];
    nv[c] = 0.0;
    nu[c] = nu[ci];
    windat->u2b[c] = 0.0;
  }
  for (cc = 1; cc <= window->nbe2; cc++) {
    c = window->bpte2[cc];
    ci = window->bine2[cc];
    nv[c] = 0.0;
    windat->u2b[c] = 0.0;
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }
  /* Set the tangential velocity boundary conditions (slip)          */
  for (; cc <= window->nbpte2; cc++) {
    c = window->bpte2[cc];
    ci = window->bine2[cc];
    nu[c] = -nu[ci];
    nv[c] = wincon->slip * nv[ci];
    windat->u2b[c] = wincon->slip * windat->u2b[ci];
    nw[c] = nw[ci];
    dzz[c] = dzz[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /* Set the sediment value */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    zm1 = window->zm1[c];
    nu[zm1] = -nu[c];
    nv[zm1] = -nv[c];
    nw[zm1] = -nw[c];
    windat->u2b[zm1] = -windat->u2b[c];
  }
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    ci = window->zp1[c];
    nu[c] = nu[ci];
    nv[c] = nv[ci];
    nw[c] = nw[ci];
    windat->u2b[c] = windat->u2b[ci];
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  for (cc = 1; cc <= vc; cc++) {
    c = cl[cc] = cells[cc];
    time_left = windat->dt;
    u = nu[c];
    v = nv[c];
    w = nw[c];
    pp = qq = 0.0;
    n=0;

    while (time_left > 0) {

      int ins = (cl[cc] == wincon->m2d[cl[cc]]) ? 1 : 0;
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      c2 = window->m2d[cl[cc]];
      p = (nu[cl[cc]] > 0.0) ? h1[c2] : h1[window->xp1[c2]];
      q = (nv[cl[cc]] > 0.0) ? h2[c2] : h2[window->yp1[c2]];
      r = (nw[cl[cc]] < 0.0) ? dzz[cl[cc]] : dzz[window->zm1[cl[cc]]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = (r) ? ((ins && nw[cl[cc]] < 0.0) ? dt : min(dt, SCALE * fabs(r / nw[cl[cc]]))) : dt;
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position.  */
      cl[cc] =
        get_posv(window, cl[cc], u, v, w, &cx[c], &cy[c], &cz[c], dt, &nc, U2GEN, &pp, &qq);
      if (nc > 2)hd_warn("Streamline distance > 2 cells @ (%d %d) : %4.2f days\n",
			 window->s2i[c], window->s2j[c], windat->days);

      c2 = window->m2d[cl[cc]];
      p = cx[c] / geth(h1, qq, c2, window->ym1);
      q = cy[c] / geth(h2, pp, c2, window->xm1);
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;

      /* Since the velocity field may be non-conservative, it is     */
      /* possible that the streamline may arrive above the surface   */
      /* or below the sea bed, hence if the origin is in the top or  */
      /* bottom layer, the distance should be limited to eta-cellz   */
      /* or cellz-botz (i.e. dzz). Limiting 0<r<1 is equivalent to   */
      /* this.                                                       */
      /* This is currently controlled with a limit on *cz in         */
      /* get_pos().                                                  */
      /*r = min(max(0.0, r), 1.0);*/

      wgt[0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - r) * (1.0 - p) * q;
      wgt[2] = (1.0 - r) * p * (1.0 - q);
      wgt[3] = (1.0 - r) * p * q;
      wgt[4] = r * (1.0 - p) * (1.0 - q);
      wgt[5] = r * (1.0 - p) * q;
      wgt[6] = r * p * (1.0 - q);
      wgt[7] = r * p * q;

      /* Interpolate x,y,z velocities */
      u = int_val(window, cl[cc], wgt, nu);
      v = int_val(window, cl[cc], wgt, nv);
      w = int_val(window, cl[cc], wgt, nw);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      n++;

      /*
      if(p > 1.0 || p < 0.0)
	hd_warn("x cell distance outside of bounds at %d : %f\n", c, p);
      if(q > 1.0 || q < 0.0)
	hd_warn("y cell distance outside of bounds at %d : %f\n", c, q);
      if(r > 1.0 || r < 0.0)
	hd_warn("r cell distance outside of bounds at %d : %f\n", c, r);
      */
    }

    sl_check(window, &cl[cc], &cx[c], &cy[c], &cz[c]);
    c2 = window->m2d[cl[cc]];
    he1 = geth(h1, qq, c2, window->ym1);
    he2 = geth(h2, pp, c2, window->xm1);
    if (wincon->osl == 0) {

      p = cx[c] / he1;
      q = cy[c] / he2;
      r = (dzz[cl[cc]]) ? cz[c] / dzz[cl[cc]] : 0.0;
      /*r = min(max(0.0, r), 1.0);*/

      /* Get the interpolation weights */
      wincon->wgt[c][0] = (1.0 - r) * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][1] = (1.0 - r) * (1.0 - p) * q;
      wincon->wgt[c][2] = (1.0 - r) * p * (1.0 - q);
      wincon->wgt[c][3] = (1.0 - r) * p * q;
      wincon->wgt[c][4] = r * (1.0 - p) * (1.0 - q);
      wincon->wgt[c][5] = r * (1.0 - p) * q;
      wincon->wgt[c][6] = r * p * (1.0 - q);
      wincon->wgt[c][7] = r * p * q;
    } else {
      if (wincon->osl == 2 || wincon->osl == 4) {
	set_eov(window, wincon, he1, he2, dzz, &cl[cc], &cx[c], &cy[c], &cz[c]);
	weights_eov(window, c, cl[cc], cx[c], cy[c], cz[c], h1, h2, dzz);
      } else
	weights_oov(window, c, cl[cc], cx[c], cy[c], cz[c], h1, h2, dzz);
    }
  }

  /* Do the final interpolation */
  memcpy(windat->nu2, windat->u2b, window->sgsiz * sizeof(double));
  if (wincon->osl == 0) {
    if(wincon->togn & TOPRIGHT) {
      for (cc = 1; cc <= vc; cc++) {
	c = cells[cc];
	c2cc[c] = cc;
	windat->nu2[c] = int_val(window, cl[cc], wincon->wgt[c], windat->u2b);
	/*
	if(isnan(windat->nu1[c])) {
	  print_tridiagonal(window,windat->u1b,c,cl[cc]);
	  s2ijk(window,c);
	  hd_quit("NaN at %f days\n",windat->days);
	}
	*/
      }
    } else {
      for (cc = 1; cc <= vc; cc++) {
	c = cells[cc];
	c2cc[c] = cc;
	windat->nu2[c] = int_val_bl(window, cl[cc], wincon->wgt[c], windat->u2b);
      }
    }
  } else {
    for (cc = 1; cc <= vc; cc++) {
      c = cells[cc];
      c2cc[c] = cc;
      windat->nu2[c] = int_valo(window, wincon->wgt[c], wincon->lmap[cl[cc]], 
				windat->u2b, wincon->osl);
    }
  }
}

/* END semi_lagrange_u2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the vertical cell spacing centred on the grid layers          */
/*-------------------------------------------------------------------*/
void set_dzzv(geometry_t *window, double *dzz, int *cells, int vc, int vcs, int mode)
{
  window_t *windat = window->windat;
  win_priv_t *wincon = window->wincon;
  int c, cc, c2, zm1, nn;
  int *bot, *dry, ndry;
  int sc, ec, *aux;
  double *dz;

  /* Set pointers */
  if (mode & U13D) {
    dz = windat->dzu1;
    bot = wincon->i5;
    dry = wincon->cdry_e1;
    ndry = wincon->ncdry_e1;
    sc = window->b2_e1;
    ec = window->n2_e1;
    aux = window->bot_e1;
  } else if (mode & U23D) {
    dz = windat->dzu2;
    bot = wincon->i6;
    dry = wincon->cdry_e2;
    ndry = wincon->ncdry_e2;
    sc = window->b2_e2;
    ec = window->n2_e2;
    aux = window->w2_e2;
  }

  /* Initialize in case a streamline lands in a dry cell, a valid    */
  /* (i.e. non inf) increment is computed.                           */
  memcpy(dzz, dz, window->sgsiz * sizeof(double));
  /* Overwrite with valid wet cell thicknesses                       */
  for (cc = 1; cc <= vc; cc++) {
    c = cells[cc];
    zm1 = window->zm1[c];
    dzz[zm1] = 0.5 * (dz[c] + dz[zm1]);
  }
  /* Set the surface layer                                           */
  for (cc = 1; cc <= vcs; cc++) {
    c = c2 = cells[cc];
    dzz[c] = 0.5 * dz[c];
    /* Set the surface and bottom boundary conditions                */
    while (c2 != window->zp1[c2]) {
      c2 = window->zp1[c2];
      dzz[c2] = dzz[c];
    }
    c2 = bot[cc];
    dzz[window->zm1[c2]] = 0.5 * dz[c2];
  }
  /* FRONT and RIGHT open boundaries                                 */
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn];
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	c2 = open->oi1_e1[cc];
	dzz[c] = dzz[c2];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	c2 = open->oi1_e2[cc];
	dzz[c] = dzz[c2];
      }
    }
  }
  /* Surface and bottom conditions over dry cells */
  for (cc = 1; cc <= ndry; cc++) {
    c = c2 = dry[cc];
    dzz[window->zm1[c]] = dzz[c];
    while (c2 != window->zp1[c2]) {
      c2 = window->zp1[c2];
      dzz[c2] = dzz[c];
    }
  }
  /* Set ghost cells                                                 */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    c2 = window->bin[cc];
    dzz[c] = dzz[c2];
  }
  for (cc = 1; cc <= window->ngsed; cc++) {
    c = window->gsed_t[cc];
    c2 = window->zp1[c];
    if(!dzz[c]) dzz[c] = dzz[c2];
  }
}

/* END set_dzzv()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to centre a velocity on the u1 or u2 point                */
/*-------------------------------------------------------------------*/
void vel_face(geometry_t *window, /* Processing window */
	      window_t *windat, /* Window data structure */
	      win_priv_t *wincon, /* Window geometry / constants */
	      double *nu,     /* Centered u1 velocity */
	      double *nv,     /* Centered u2 velocity */
	      double *nw,     /* Centered z velocity */
	      int *cells,
	      int vcs,
	      int vc,
	      int mode
  )
{
  int cc, c, n, ncells;
  double w, wtop;
  int *m1, *obc, nobc;

  if (mode & (U13D|U12D)) {
    m1 = window->xm1;
  } else if (mode & (U23D|U22D)) {
    m1 = window->ym1;
  }
  if (mode & (U13D|U23D))
    ncells = vc;
  else
    ncells = vcs;

  /* Face centre the velocities */
  for (cc = 1; cc <= ncells; cc++) {
    c = cells[cc];
    if (fabs(nu[c]) < SMALL) nu[c] = SMALL;
    if (fabs(nv[c]) < SMALL) nv[c] = SMALL;
  }
  if (mode & (U13D|U23D)) {
    for (cc = 1; cc <= vc; cc++) {
      if (cc <= vcs)
	wtop = 0.5 * (windat->wtop[window->m2d[c]] + windat->wtop[window->m2d[m1[c]]]);
      else
	wtop = 0.5 * (windat->w[window->zp1[c]] + windat->w[m1[window->zp1[c]]]);
      w = 0.5 * (windat->w[c] + windat->w[m1[c]]);
      nw[c] = 0.5 * (w + wtop);
      if (fabs(nw[c]) < SMALL)
	nw[c] = SMALL;
    }
  }

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (mode & U13D) {
      nobc = open->no3_e1;
      obc = open->obc_e1;
    } else if (mode & U23D) {
      nobc = open->no3_e2;
      obc = open->obc_e2;
    } else if (mode & U12D) {
      nobc = open->no2_e1;
      obc = open->obc_e1;
    } else if (mode & U22D) {
      nobc = open->no2_e2;
      obc = open->obc_e2;
    }
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= nobc; cc++) {
	c = obc[cc];
	nu[c] = (nu[c] < 0.0) ? 0.0 : nu[c];
      }
    }
    if (open->ocodex & L_EDGE) {
      for (cc = 1; cc <= nobc; cc++) {
	c = obc[cc];
	nu[c] = (nu[c] > 0.0) ? 0.0 : nu[c];
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= nobc; cc++) {
	c = obc[cc];
	nv[c] = (nv[c] < 0.0) ? 0.0 : nv[c];
      }
    }
    if (open->ocodey & B_EDGE) {
      for (cc = 1; cc <= nobc; cc++) {
	c = obc[cc];
	nv[c] = (nv[c] > 0.0) ? 0.0 : nv[c];
      }
    }
  }

  if (mode & (U12D|U22D)) return;

  /* Check for undefined velocities */
    for (cc = 1; cc <= ncells; cc++) {
      c = cells[cc];
      if (fabs(nu[c]) > wincon->velmax || isnan(nu[c])) {
	emstag(LWARN,"Lagrange:","Undefined u1 velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nu[c]);
	nu[c] = SMALL;
      }
      if (fabs(nv[c]) > wincon->velmax || isnan(nv[c])) {
	emstag(LWARN,"Lagrange:","Undefined u2 velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nv[c]);
	nv[c] = SMALL;
      }
      if (fabs(nw[c]) > wincon->velmax || isnan(nw[c])) {
	emstag(LWARN,"Lagrange:","Undefined w velocity at %f (%d %d %d) : %f\n", windat->t/86400, window->s2i[c], window->s2j[c], window->s2k[c], nw[c]);
	nw[c] = SMALL;
      }
    }
}

/* END vel_face()                                                    */
/*-------------------------------------------------------------------*/

#define TINY_INSIDE (1e-8)

/*-------------------------------------------------------------------*/
/* Locates the position of the streamline in sparse coordinates.     */
/* The origin for interpolations is in the top right horizontal and  */
/* bottom right vertical corner of the cell. Due to the orthogonal   */
/* grid (non-uniform grid spacing) the streamline distance must be   */
/* decremented for each individual cell.                             */
/*-------------------------------------------------------------------*/
int get_posv(geometry_t *window, int c, double nu, double nv, double nw,
	     double *cx, double *cy, double *cz, double dt, int *nc, 
	     int dir, double *p, double *q)
{
  win_priv_t *wincon = window->wincon;
  double *dzz = wincon->w8;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  int co, cm;
  double d, h;
  int ncx, ncy;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;
  if (nw == SMALL)
    nw = 0.0;
  ncx = ncy = 0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  d = nu * dt + (*cx);
  if (d < 0.0) {
    co = window->xp1[co];
    h = geth(h1, *q, window->m2d[co], window->ym1);
    ncx = 1;
    while (fabs(d) >= h) {
      d += h;
      co = window->xp1[co];
      h = geth(h1, *q, window->m2d[co], window->ym1);
      ncx++;
    }
    d += h;
  }
  if (d > 0.0) {
    h = geth(h1, *q, window->m2d[co], window->ym1);
    while (fabs(d) >= h) {
      /* Prevent for intrusions into land                            */
      if (dir & U1GEN && wincon->gmap[co] & L_EDGE) {
	co = window->xp1[co];
	h = geth(h1, *q, window->m2d[co], window->ym1);
	d = h - TINY_INSIDE;
	break;
      }
      d -= h;
      co = window->xm1[co];
      ncx--;
    }
  }
  *cx = fabs(d);
  *p = *cx / geth(h1, *q, window->m2d[co], window->ym1);

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = nv * dt + (*cy);
  if (d < 0.0) {
    co = window->yp1[co];
    h = geth(h2, *p, window->m2d[co], window->xm1);
    ncy = 1;
    while (fabs(d) >= h) {
      d += h;
      co = window->yp1[co];
      h = geth(h2, *p, window->m2d[co], window->xm1);
      ncy++;
    }
    d += h;
  }
  if (d > 0.0) {
    h = geth(h2, *p, window->m2d[co], window->xm1);
    while (fabs(d) >= h) {
      /* Prevent for intrusions into land                            */
      if (dir & U2GEN && wincon->gmap[co] & B_EDGE) {
	co = window->yp1[co];
	h = geth(h2, *p, window->m2d[co], window->xm1);
	d = h - TINY_INSIDE;
	break;
      }
      d -= h;
      co = window->ym1[co];
      ncy--;
    }
  }
  *cy = fabs(d);
  *q = *cy / geth(h2, *p, window->m2d[co], window->xm1);
  *nc  = (fabs(ncx) > fabs(ncy)) ? ncx : ncy;

  /*-----------------------------------------------------------------*/
  /* z direction.                                                    */
  d = -nw * dt + (*cz);
  if (d < 0.0) {
    co = window->zm1[co];
    while (fabs(d) >= dzz[co]) {
      if (co == window->zm1[co]) {
	*cz = fabs(fmod(d, dzz[co]));
	d += dzz[co];
	return(co);
      }
      d += dzz[co];
      co = window->zm1[co];
    }
    d += dzz[co];
  }
  if (d > 0.0) {
    while (d >= dzz[co] && co != window->zp1[co]) {
      if (co == wincon->m2d[co]) {
	*cz = fmod(d, dzz[co]);
	return(co);
      }
      d -= dzz[co];
      co = window->zp1[co];
    }
  }
  *cz = fabs(d);
  /* Sanity check: *cz must be < dzz */
  if (*cz > dzz[co]) *cz = fmod(*cz, dzz[co]);

  return (co);
}

int get_posvS(geometry_t *window, int c, double nu, double nv, 
	      double *cx, double *cy, double dt, int dir, double *p, double *q)
{
  win_priv_t *wincon = window->wincon;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  int co, cm;
  double d, h;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  d = nu * dt + (*cx);
  if (d < 0.0) {
    co = window->xp1[co];
    h = geth(h1, *q, window->m2d[co], window->ym1);
    while (fabs(d) >= h) {
      d += h;
      co = window->xp1[co];
      h = geth(h1, *q, window->m2d[co], window->ym1);
    }
    d += h;
  }
  if (d > 0.0) {
    h = geth(h1, *q, window->m2d[co], window->ym1);
    while (fabs(d) >= h) {
      /* Prevent for intrusions into land                            */
      if (dir & U1GEN && wincon->gmap[co] & L_EDGE) {
	co = window->xp1[co];
	h = geth(h1, *q, window->m2d[co], window->ym1);
	d = h - TINY_INSIDE;
	break;
      }
      d -= h;
      co = window->xm1[co];
    }
  }
  *cx = fabs(d);
  *p = *cx / geth(h1, *q, window->m2d[co], window->ym1);

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = nv * dt + (*cy);
  if (d < 0.0) {
    co = window->yp1[co];
    h = geth(h2, *p, window->m2d[co], window->xm1);
    while (fabs(d) >= h) {
      d += h;
      co = window->yp1[co];
      h = geth(h2, *p, window->m2d[co], window->xm1);
    }
    d += h;
  }
  if (d > 0.0) {
    h = geth(h2, *p, window->m2d[co], window->xm1);
    while (fabs(d) >= h) {
      /* Prevent for intrusions into land                            */
      if (dir & U2GEN && wincon->gmap[co] & B_EDGE) {
	co = window->yp1[co];
	h = geth(h2, *p, window->m2d[co], window->xm1);
	d = h - TINY_INSIDE;
	break;
      }
      d -= h;
      co = window->ym1[co];
    }
  }
  *cy = fabs(d);
  *q = *cy / geth(h2, *p, window->m2d[co], window->xm1);

  return (co);
}

double geth(double *h, double f, int c, int *m1)
{
  int cm = m1[c];
  double h1 = h[c];
  double h2 = h[m1[c]];
  return(f * (h2 - h1) + h1);
}

int get_posvo(geometry_t *window, int c, double nu, double nv, double nw,
	     double *cx, double *cy, double *cz, double dt, int *nc, int dir)
{
  win_priv_t *wincon = window->wincon;
  double *dzz = wincon->w8;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  int co, cm;
  double d;
  int ncx, ncy;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;
  if (nw == SMALL)
    nw = 0.0;
  ncx = ncy = 0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  d = nu * dt + (*cx);
  if (d < 0.0) {
    co = window->xp1[co];
    ncx = 1;
    while (fabs(d) >= h1[window->m2d[co]]) {
      d += h1[window->m2d[co]];
      co = window->xp1[co];
      ncx++;
    }
    d += h1[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= h1[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      if (dir & U1GEN && wincon->gmap[co] & L_EDGE) {
	co = window->xp1[co];
	d = h1[window->m2d[co]] - TINY_INSIDE;
	break;
      }
      d -= h1[window->m2d[co]];
      co = window->xm1[co];
      ncx--;
    }
  }
  *cx = fabs(d);

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = nv * dt + (*cy);
  if (d < 0.0) {
    co = window->yp1[co];
    ncy = 1;
    while (fabs(d) >= h2[window->m2d[co]]) {
      d += h2[window->m2d[co]];
      co = window->yp1[co];
      ncy++;
    }
    d += h2[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= h2[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      if (dir & U2GEN && wincon->gmap[co] & B_EDGE) {
	co = window->yp1[co];
	d = h2[window->m2d[co]] - TINY_INSIDE;
	break;
      }
      d -= h2[window->m2d[co]];
      co = window->ym1[co];
      ncy--;
    }
  }
  *cy = fabs(d);

  /*-----------------------------------------------------------------*/
  /* Adjust for converging curvilinear grids                         */
  while (*cx >= h1[window->m2d[co]] ||
	 *cy >= h2[window->m2d[co]]) {
    if (*cx >= h1[window->m2d[co]]) {
      cm = window->xm1[co];
      ncx--;
      if (cm != window->xm1[cm]) {
	*cx -= h1[window->m2d[co]];
	co = cm;
      } else /* Streamline heading into land - keep it wet           */
	*cx = fabs(*cx - h1[window->m2d[co]]);
    }
    if (*cy >= h2[window->m2d[co]]) {
      cm = window->ym1[co];
      ncy--;
      if (cm != window->ym1[cm]) {
	*cy -= h2[window->m2d[co]];
	co = window->ym1[co];
      } else
	*cy = fabs(*cy - h2[window->m2d[co]]);
    }
  }
  *nc  = (fabs(ncx) > fabs(ncy)) ? ncx : ncy;

/*-----------------------------------------------------------------*/
  /* z direction.                                                    */
  d = -nw * dt + (*cz);
  if (d < 0.0) {
    co = window->zm1[co];
    while (fabs(d) >= dzz[co]) {
      if (co == window->zm1[co]) {
	*cz = fabs(fmod(d, dzz[co]));
	d += dzz[co];
	return(co);
      }
      d += dzz[co];
      co = window->zm1[co];
    }
    d += dzz[co];
  }
  if (d > 0.0) {
    while (d >= dzz[co] && co != window->zp1[co]) {
      if (co == wincon->m2d[co]) {
	*cz = fmod(d, dzz[co]);
	return(co);
      }
      d -= dzz[co];
      co = window->zp1[co];
    }
  }
  *cz = fabs(d);
  /* Sanity check: *cz must be < dzz */
  if (*cz > dzz[co]) *cz = fmod(*cz, dzz[co]);

  return (co);
}

int get_posvSo(geometry_t *window, int c, double nu, double nv, 
	      double *cx, double *cy, double dt, int dir)
{
  win_priv_t *wincon = window->wincon;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  int co, cm;
  double d;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  co = c;
  if (nu == SMALL)
    nu = 0.0;
  if (nv == SMALL)
    nv = 0.0;

  /*-----------------------------------------------------------------*/
  /* x direction.                                                    */
  d = nu * dt + (*cx);
  if (d < 0.0) {
    co = window->xp1[co];
    while (fabs(d) >= h1[window->m2d[co]]) {
      d += h1[window->m2d[co]];
      co = window->xp1[co];
    }
    d += h1[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= h1[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      if (dir & U1GEN && wincon->gmap[co] & L_EDGE) {
	co = window->xp1[co];
	d = h1[window->m2d[co]] - TINY_INSIDE;
	break;
      }
      d -= h1[window->m2d[co]];
      co = window->xm1[co];
    }
  }
  *cx = fabs(d);

  /*-----------------------------------------------------------------*/
  /* y direction.                                                    */
  d = nv * dt + (*cy);
  if (d < 0.0) {
    co = window->yp1[co];
    while (fabs(d) >= h2[window->m2d[co]]) {
      d += h2[window->m2d[co]];
      co = window->yp1[co];
    }
    d += h2[window->m2d[co]];
  }
  if (d > 0.0) {
    while (fabs(d) >= h2[window->m2d[co]]) {
      /* Prevent for intrusions into land                            */
      if (dir & U2GEN && wincon->gmap[co] & B_EDGE) {
	co = window->yp1[co];
	d = h2[window->m2d[co]] - TINY_INSIDE;
	break;
      }
      d -= h2[window->m2d[co]];
      co = window->ym1[co];
    }
  }
  *cy = fabs(d);

  /*-----------------------------------------------------------------*/
  /* Adjust for converging curvilinear grids                         */
  while (*cx >= h1[window->m2d[co]] ||
	 *cy >= h2[window->m2d[co]]) {
    if (*cx >= h1[window->m2d[co]]) {
      cm = window->xm1[co];
      if (cm != window->xm1[cm]) {
	*cx -= h1[window->m2d[co]];
	co = window->xm1[co];
      } else
	*cx = fabs(*cx - h1[window->m2d[co]]);
    }
    if (*cy >= h2[window->m2d[co]]) {
      cm = window->ym1[co];
      if (cm != window->ym1[cm]) {
	*cy -= h2[window->m2d[co]];
	co = window->ym1[co];
      } else
	*cy = fabs(*cy - h2[window->m2d[co]]);
    }
  }
  return (co);
}

/* END get_posv()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the weights for interpolations of order 1 & 3. Uses Eqn. 5   */ 
/* McDonald (1984) Mon. Wea. Rev. 112, 1267-1275. The distances are  */
/* re-mapped to local distances, where h[0] = 0. Note that (Eqn. 7)  */
/* hx[i-1] < p < h[i], hy[j-1] < q < hy[j] and hz[i] < r < hz[i+1].  */
/*-------------------------------------------------------------------*/
void weights_oov(geometry_t *window, 
		 int c,            /* Destination cell               */
		 int co,           /* Source cell                    */
		 double pin,       /* x distance from dest origin    */
		 double qin,       /* y distance from dest origin    */
		 double rin,       /* z distance from dest origin    */
		 double *h1,       /* e1 distance between faces      */
		 double *h2,       /* e2 distance between faces      */
		 double *dzz)      /* Distance between layer centres */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj, ck;
  double p, q, r;
  double hx[5], hy[5], hz[5];

  /* Set up the normalized distances */
  ck = co;
  ci = cj = window->m2d[co];
  p = h1[ci] - pin;
  q = h2[cj] - qin;
  r = rin;
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    ck = window->zm1[ck];
    p += h1[ci];
    q += h2[cj];
    r += dzz[ck];
  }
  hx[0] = hy[0] = hz[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + h1[ci];
    hy[m] = hy[m-1] + h2[cj];
    hz[m] = hz[m-1] + dzz[ck];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
    ck = window->zp1[ck];
  }
  /* Loop to get the weights */
  in = 0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	wincon->wgt[c][in] = 1.0;
	for (m = 0; m <= n; m++) {
	  if (m != k)
	    wincon->wgt[c][in] *= (r - hz[m]) / (hz[k] - hz[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != j)
	    wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != i)
	    wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	}
	in++;
      }
  }
}

void weights_oovS(geometry_t *window, 
		  int c,            /* Destination cell               */
		  int co,           /* Source cell                    */
		  double pin,       /* x distance from dest origin    */
		  double qin,       /* y distance from dest origin    */
		  double *h1,       /* e1 distance between faces      */
		  double *h2)       /* e2 distance between faces      */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj;
  double p, q;
  double hx[5], hy[5];

  for (m = 0; m < wincon->nosl; m++) wincon->wgt[c][m] = 0.0;

  /* Set up the normalized distances */
  ci = cj = co;
  p = h1[ci] - pin;
  q = h2[cj] - qin;
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    p += h1[ci];
    q += h2[cj];
  }
  hx[0] = hy[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + h1[ci];
    hy[m] = hy[m-1] + h2[cj];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
  }

  /* Loop to get the weights */
  in = 0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	if (k == wincon->mosl) {
	  wincon->wgt[c][in] = 1.0;
	  for (m = 0; m <= n; m++) {
	    if (m != j)
	      wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	  }
	  for (m = 0; m <= n; m++) {
	    if (m != i)
	      wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	  }
	}
	in++;
      }
  }

}

/* END weights_oov()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the weights for interpolations of order 2 & 4. Uses Eqn. 5   */ 
/* McDonald (1984) Mon. Wea. Rev. 112, 1267-1275. The distances are  */
/* re-mapped to local distances, where h[0] = 0. Note that (Eqn. 8)  */
/* hx[i-1/2] < p < h[i+1/2], hy[j-1/2] < q < hy[j+1/2] and           */
/* hz[i-1/2] < r < hz[i+1/2].                                        */
/* Note: weights can be negative.                                    */
/*-------------------------------------------------------------------*/
void weights_eov(geometry_t *window, 
		 int c,            /* Destination cell               */
		 int co,           /* Source cell                    */
		 double pin,       /* x distance from dest origin    */
		 double qin,       /* y distance from dest origin    */
		 double rin,       /* z distance from dest origin    */
		 double *h1,       /* e1 distance between faces      */
		 double *h2,       /* e2 distance between faces      */
		 double *dzz)      /* Distance between layer faces   */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj, ck;
  double p, q, r, d;
  double hx[5], hy[5], hz[5];

  /* Set up the normalized distances */
  ck = co;
  ci = cj = window->m2d[co];
  p = h1[ci] - pin;
  q = h2[ci] - qin;
  r = dzz[ck] + rin;
  ci = window->xm1[ci];
  cj = window->ym1[cj];
  ck = window->zm1[ck];
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    ck = window->zm1[ck];
    p += h1[ci];
    q += h2[cj];
    r += dzz[ck];
  }
  hx[0] = hy[0] = hz[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + h1[ci];
    hy[m] = hy[m-1] + h2[cj];
    hz[m] = hz[m-1] + dzz[ck];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
    ck = window->zp1[ck];
  }

  /* Loop to get the weights */
  in = 0;
  d = 0.0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	wincon->wgt[c][in] = 1.0;
	for (m = 0; m <= n; m++) {
	  if (m != k)
	    wincon->wgt[c][in] *= (r - hz[m]) / (hz[k] - hz[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != j)
	    wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	}
	for (m = 0; m <= n; m++) {
	  if (m != i)
	    wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	}
	d += wincon->wgt[c][in];
	in++;
      }
  }
  /*
  if (fabs(d - 1.0) > SMALL) hd_warn("Weights do not sum to 1 at %d : %e\n",c, d);
  */
}


void weights_eovS(geometry_t *window, 
		  int c,            /* Destination cell               */
		  int co,           /* Source cell                    */
		  double pin,       /* x distance from dest origin    */
		  double qin,       /* y distance from dest origin    */
		  double *h1,       /* e1 distance between faces      */
		  double *h2)       /* e2 distance between faces      */
{
  win_priv_t *wincon = window->wincon;
  int n = wincon->osl;
  int i, j, k, m, in;
  int ci, cj;
  double p, q;
  double hx[5], hy[5];

  for (m = 0; m < wincon->nosl; m++) wincon->wgt[c][m] = 0.0;

  /* Set up the normalized distances */
  ci = cj = co;
  p = h1[ci] - pin;
  q = h2[ci] - qin;
  ci = window->xm1[ci];
  cj = window->ym1[cj];
  if (n > 2) {
    ci = window->xm1[ci];
    cj = window->ym1[cj];
    p += h1[ci];
    q += h2[cj];
  }
  hx[0] = hy[0] = 0.0;
  for (m = 1; m <= n; m++) {
    hx[m] = hx[m-1] + h1[ci];
    hy[m] = hy[m-1] + h2[cj];
    ci = window->xp1[ci];
    cj = window->yp1[cj];
  }

  /* Loop to get the weights */
  in = 0;
  for (i = 0; i <= n; i++) {
    for (j = 0; j <= n; j++)
      for (k = 0; k <= n; k++) {
	if (k == wincon->mosl) {
	  wincon->wgt[c][in] = 1.0;
	  for (m = 0; m <= n; m++) {
	    if (m != j)
	      wincon->wgt[c][in] *= (q - hy[m]) / (hy[j] - hy[m]);
	  }
	  for (m = 0; m <= n; m++) {
	    if (m != i)
	      wincon->wgt[c][in] *= (p - hx[m]) / (hx[i] - hx[m]);
	  }
	}
	in++;
      }
  }
}

/* END weights_eov()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the origin  for tri-quadratic and tri-quartic    */
/* interpolations. See McDonald (1984) Mon. Wea. Rev. 112, 1267-1275 */
/* Eqn. 8.                                                           */ 
/*-------------------------------------------------------------------*/
void set_eov(geometry_t *window, win_priv_t *wincon, double h1, double h2, double *dzz,
	     int *c, double *cx, double *cy, double *cz)
{
  int co = *c;

  if (*cx > 0.5 * h1) {
    *cx -= h1;
    co = window->xm1[co];
  }
  if (*cy > 0.5 * h2) {
    *cy -= h2;
    co = window->ym1[co];
  }
  if (*cz > 0.5 * wincon->dz[co]) {
    *cz -= dzz[co];
    co = window->zp1[co];
  }
  *c = co;
}

void set_eovS(geometry_t *window, win_priv_t *wincon, double h1, double h2, 
	      int *c, double *cx, double *cy)
{
  int co = *c;
  double h;

  if (*cx > 0.5 * h1) {
    *cx -= h1;
    co = window->xm1[co];
  }
  if (*cy > 0.5 * h2) {
    *cy -= h2;
    co = window->ym1[co];
  }
  *c = co;
}

/* END set_eo()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-Lagrangian advection for u1av component of velocity.         */
/*-------------------------------------------------------------------*/
void semi_lagrange_u1av(geometry_t *window,  /* Processing window    */
			window_t *windat,    /* Window data          */
			win_priv_t *wincon   /* Window constants     */
			)
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->w1;
  double *nv = wincon->w2;
  double *cx = wincon->w4;
  double *cy = wincon->w5;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  double u, v, he1, he2;
  double p, q, pp, qq;
  int *cells, vcs;
  int cc, c, ci, c2;
  int *cl = wincon->s5;        /* Streamline origin */
  int *c2cc = wincon->i7;      /* c to cc map */
  int n;

  if(wincon->dolin_u1) {
    cells = wincon->s4;
    vcs = wincon->aclS;
  } else {
    cells = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(cx, 0, window->sgsizS * sizeof(double));
  memset(cy, 0, window->sgsizS * sizeof(double));
  memcpy(nu, windat->u1av, window->sgsizS * sizeof(double));
  memcpy(nv, wincon->w7, window->sgsizS * sizeof(double));
  memset(cl, 0, window->sgsizS * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the cell spacing */
  for (c = 1; c <= window->enonS; c++) {
    h1[c] = window->h1acell[window->xm1[c]];
    h2[c] = 0.5 * (window->h2au2[c] + window->h2au2[window->xm1[c]]);
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face */
  vel_face(window, windat, wincon, nu, nv, NULL, cells, vcs, 0, U12D);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    ci = window->bin[cc];
    nu[c] = 0.0;
    nv[c] = nv[ci];
    windat->u1avb[c] = 0.0;
  }
  for (cc = 1; cc <= window->nbe1S; cc++) {
    c = window->bpte1S[cc];
    ci = window->bine1S[cc];
    nu[c] = 0.0;
    /*nv[c] = wincon->slip * nv[ci];*/
    windat->u1avb[c] = 0.0;
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }
  /* Set the tangential velocity boundary conditions (slip)          */
  for (; cc <= window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    ci = window->bine1S[cc];
    nu[c] = wincon->slip * nu[ci];
    nv[c] = -nv[ci];
    windat->u1avb[c] = wincon->slip * windat->u1avb[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  for (cc = 1; cc <= vcs; cc++) {
    c = cg = cl[cc] = cells[cc];
    time_left = windat->dt2d;
    u = nu[c];
    v = nv[c];
    pp = qq = 0.0;
    n=0;

    while (time_left > 0) {
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      p = (nu[cl[cc]] > 0.0) ? h1[c] : h1[window->xp1[c]];
      q = (nv[cl[cc]] > 0.0) ? h2[c] : h2[window->yp1[c]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position.  */
      cl[cc] = c2 = get_posvS(window, cl[cc], u, v, &cx[c], &cy[c], dt, U1GEN, &pp, &qq);

      p = cx[c] / geth(h1, qq, c2, window->ym1);
      q = cy[c] / geth(h2, pp, c2, window->xm1);

      wgt[0] = (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - p) * q;
      wgt[2] = p * (1.0 - q);
      wgt[3] = p * q;
      wgt[4] = 0.0;
      wgt[5] = 0.0;
      wgt[6] = 0.0;
      wgt[7] = 0.0;

      /* Interpolate x,y,z velocities */
      u = int_val(window, cl[cc], wgt, nu);
      v = int_val(window, cl[cc], wgt, nv);
      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      n++;
    }

    he1 = geth(h1, qq, c2, window->ym1);
    he2 = geth(h2, pp, c2, window->xm1);
    if (wincon->osl == 0) {
      p = cx[c] / he1;
      q = cy[c] / he2;

      /* Get the interpolation weights */
      wincon->wgt[c][0] = (1.0 - p) * (1.0 - q);
      wincon->wgt[c][1] = (1.0 - p) * q;
      wincon->wgt[c][2] = p * (1.0 - q);
      wincon->wgt[c][3] = p * q;
      wincon->wgt[c][4] = 0.0;
      wincon->wgt[c][5] = 0.0;
      wincon->wgt[c][6] = 0.0;
      wincon->wgt[c][7] = 0.0;
    } else {
      if (wincon->osl == 2 || wincon->osl == 4) {
	set_eovS(window, wincon, he1, he2, &cl[cc], &cx[c], &cy[c]);
	weights_eovS(window, c, cl[cc], cx[c], cy[c], h1, h2);
      } else
	weights_oovS(window, c, cl[cc], cx[c], cy[c], h1, h2);
    }
  }

  /* Do the final interpolation */
  memcpy(windat->nu1av, windat->u1avb, window->sgsizS * sizeof(double));
  if (wincon->osl == 0) {
    if(wincon->togn & TOPRIGHT) {
      for (cc = 1; cc <= vcs; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	/*check_tridiagonal(window,windat->u1avb,"before",c,cl[cc]);*/
	windat->nu1av[c] = int_val(window, cl[cc], wincon->wgt[c], windat->u1avb);
	/*
	if(isnan(windat->nu1av[c])) {
	  print_tridiagonal(window,windat->u1avb,c,cl[cc]);
	  s2ijk(window,c);
	  hd_quit("NaN at %f days\n",windat->days);
	}
	*/
      }
    } else {
      for (cc = 1; cc <= vcs; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	windat->nu1av[c] = int_val_bl(window, cl[cc], wincon->wgt[c], windat->u1avb);
      }
    }
  } else {
    for (cc = 1; cc <= vcs; cc++) {
      c = cg = cells[cc];
      c2cc[c] = cc;
      windat->nu1av[c] = int_valo(window, wincon->wgt[c], wincon->lmap[cl[cc]], 
				  windat->u1avb, wincon->osl);
    }
  }
}

/* END semi_lagrange_u1av()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Semi-Lagrangian advection for u1av component of velocity.         */
/*-------------------------------------------------------------------*/
void semi_lagrange_u2av(geometry_t *window,  /* Processing window    */
			window_t *windat,    /* Window data          */
			win_priv_t *wincon   /* Window constants     */
			)
{
  double dt;                    /* Sub-time step for Euler velocity
                                   interp.  */
  double time_left;             /* Number of sub-timesteps per dt */
  double SCALE = 0.95;          /* Use 95% of allowable sub-timestep */
  double wgt[8];                /* Weights for trilinear interpolation */
  double *nu = wincon->w1;
  double *nv = wincon->w2;
  double *cx = wincon->w4;
  double *cy = wincon->w5;
  double *h1 = wincon->d3;
  double *h2 = wincon->d4;
  double u, v, he1, he2;
  double p, q, pp, qq;
  int *cells, vcs;
  int cc, c, ci, c2;
  int *cl = wincon->s5;        /* Streamline origin */
  int *c2cc = wincon->i7;      /* c to cc map */
  int n;

  if(wincon->dolin_u2) {
    cells = wincon->s4;
    vcs = wincon->aclS;
  } else {
    cells = wincon->s3;
    vcs = wincon->vcs;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise */
  memset(cx, 0, window->sgsizS * sizeof(double));
  memset(cy, 0, window->sgsizS * sizeof(double));
  memcpy(nu, wincon->w7, window->sgsizS * sizeof(double));
  memcpy(nv, windat->u2av, window->sgsizS * sizeof(double));
  memset(cl, 0, window->sgsizS * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Set the cell spacing */
  for (c = 1; c <= window->enonS; c++) {
    h1[c] = 0.5 * (window->h1au1[c] + window->h1au1[window->ym1[c]]);
    h2[c] = window->h2acell[window->ym1[c]];
  }

  /*-----------------------------------------------------------------*/
  /* Center the velocities on the appropriate face */
  vel_face(window, windat, wincon, nu, nv, NULL, cells, vcs, 0, U22D);

  /* Set the ghost cells */
  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    nv[c] = 0.0;
    nu[c] = nu[ci];
    windat->u2avb[c] = 0.0;
  }
  for (cc = 1; cc <= window->nbe2S; cc++) {
    c = window->bpte2S[cc];
    ci = window->bine2S[cc];
    nv[c] = 0.0;
    windat->u2avb[c] = 0.0;
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }
  /* Set the tangential velocity boundary conditions (slip)          */
  for (; cc <= window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    ci = window->bine2S[cc];
    nu[c] = -nu[ci];
    nv[c] = wincon->slip * nv[ci];
    windat->u2avb[c] = wincon->slip * windat->u2avb[ci];
    wincon->m2d[c] = (ci == wincon->m2d[ci]) ? c : 0;
  }

  /*-----------------------------------------------------------------*/
  /* Trace the streamline back to the origin */
  for (cc = 1; cc <= vcs; cc++) {
    c = cg = cl[cc] = cells[cc];
    time_left = windat->dt2d;
    u = nu[c];
    v = nv[c];
    pp = qq = 0.0;
    n=0;

    while (time_left > 0) {
      /* Get the sub-timestep for the next integration. Streamlines */
      /* cannot cross more than one cell in one sub-timestep.  */
      p = (nu[cl[cc]] > 0.0) ? h1[c] : h1[window->xp1[c]];
      q = (nv[cl[cc]] > 0.0) ? h2[c] : h2[window->yp1[c]];
      dt = min(SCALE * fabs(p / nu[cl[cc]]),
               SCALE * fabs(q / nv[cl[cc]]));
      dt = min(dt, time_left);

      /* Get the bottom left grid coordinate defining the streamline */
      /* position.  */
      cl[cc] = c2 = get_posvS(window, cl[cc], u, v, &cx[c], &cy[c], dt, U2GEN, &pp, &qq);

      p = cx[c] / geth(h1, qq, c2, window->ym1);
      q = cy[c] / geth(h2, pp, c2, window->xm1);

      wgt[0] = (1.0 - p) * (1.0 - q);
      wgt[1] = (1.0 - p) * q;
      wgt[2] = p * (1.0 - q);
      wgt[3] = p * q;
      wgt[4] = 0.0;
      wgt[5] = 0.0;
      wgt[6] = 0.0;
      wgt[7] = 0.0;

      /* Interpolate x,y,z velocities */
      u = int_val(window, cl[cc], wgt, nu);
      v = int_val(window, cl[cc], wgt, nv);

      /* Get the number of sub-timesteps used */
      time_left = time_left - dt;
      n++;
    }
    he1 = geth(h1, qq, c2, window->ym1);
    he2 = geth(h2, pp, c2, window->xm1);

    if (wincon->osl == 0) {
      p = cx[c] / he1;
      q = cy[c] / he2;

      /* Get the interpolation weights */
      wincon->wgt[c][0] = (1.0 - p) * (1.0 - q);
      wincon->wgt[c][1] = (1.0 - p) * q;
      wincon->wgt[c][2] = p * (1.0 - q);
      wincon->wgt[c][3] = p * q;
      wincon->wgt[c][4] = 0.0;
      wincon->wgt[c][5] = 0.0;
      wincon->wgt[c][6] = 0.0;
      wincon->wgt[c][7] = 0.0;
    } else {
      if (wincon->osl == 2 || wincon->osl == 4) {
	set_eovS(window, wincon, he1, he2, &cl[cc], &cx[c], &cy[c]);
	weights_eovS(window, c, cl[cc], cx[c], cy[c], h1, h2);
      } else
	weights_oovS(window, c, cl[cc], cx[c], cy[c], h1, h2);
    }
  }

  /* Do the final interpolation */
  memcpy(windat->nu2av, windat->u2avb, window->sgsizS * sizeof(double));
  if (wincon->osl == 0) {
    if(wincon->togn & TOPRIGHT) {
      for (cc = 1; cc <= vcs; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	/*check_tridiagonal(window,windat->u2avb,"before",c,cl[cc]);*/
	windat->nu2av[c] = int_val(window, cl[cc], wincon->wgt[c], windat->u2avb);
	/*
	if(isnan(windat->nu2av[c])) {
	  print_tridiagonal(window,windat->u2avb,c,cl[cc]);
	  s2ijk(window,c);
	  hd_quit("NaN at %f days\n",windat->days);
	}
	*/
      }
    } else {
      for (cc = 1; cc <= vcs; cc++) {
	c = cg = cells[cc];
	c2cc[c] = cc;
	windat->nu2av[c] = int_val_bl(window, cl[cc], wincon->wgt[c], windat->u2avb);
      }
    }
  } else {
    for (cc = 1; cc <= vcs; cc++) {
      c = cg = cells[cc];
      c2cc[c] = cc;
      windat->nu2av[c] = int_valo(window, wincon->wgt[c], wincon->lmap[cl[cc]], 
				  windat->u2avb, wincon->osl);
    }
  }
}

/* END semi_lagrange_u2av()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implement the porus plate algorithm in the e1 direction           */
/*-------------------------------------------------------------------*/
double porus_plate_e1(geometry_t *window,    /* Processing window    */
		      window_t *windat,      /* Window data          */
		      win_priv_t *wincon,    /* Window constants     */
		      int c
		      )
{
  int c2;
  double Cdrag = 1.5;        /* Drag coefficient for coral reefs     */
  double eff = 0.35;         /* Effective wet cross sectional area   */
  double closs;              /* Energy loss coefficient              */
  double d1;

  c2 = window->m2d[c];
  if (windat->reefe1[c] > 0.0) {
    /* Get the energy loss coefficient                             */
    d1 = 1.0 / eff;
    closs = Cdrag * windat->reefe1[c] * d1 * d1 / 2.0;
    /* Get the reef friction term                                  */
    d1 = windat->u1[c] * windat->u1[c] + windat->u2[c] * windat->u2[c];
    return(-closs * windat->u1[c] * sqrt(d1) / window->h1au1[c2]);
  }
  return(0.0);
}

/* END porus_plate_e1()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implement the porus plate algorithm in the e2 direction           */
/*-------------------------------------------------------------------*/
double porus_plate_e2(geometry_t *window,    /* Processing window    */
		      window_t *windat,      /* Window data          */
		      win_priv_t *wincon,    /* Window constants     */
		      int c
		      )
{
  int c2;
  double Cdrag = 1.5;        /* Drag coefficient for coral reefs     */
  double eff = 0.35;         /* Effective wet cross sectional area   */
  double closs;              /* Energy loss coefficient              */
  double d1;

  c2 = window->m2d[c];
  if (windat->reefe2[c] > 0.0) {
    /* Get the energy loss coefficient                             */
    d1 = 1.0 / eff;
    closs = Cdrag * windat->reefe2[c] * d1 * d1 / 2.0;
    /* Get the reef friction term                                  */
    d1 = windat->u1[c] * windat->u1[c] + windat->u2[c] * windat->u2[c];
    return(-closs * windat->u2[c] * sqrt(d1) / window->h2au2[c2]);
  }
  return(0.0);
}

/* END porus_plate_e2()                                              */
/*-------------------------------------------------------------------*/


