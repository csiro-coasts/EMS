/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/slaves/zoom.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: zoom.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void geom_sinterp_e1(geometry_t *geom, geometry_t **window);
void geom_sinterp_e2(geometry_t *geom, geometry_t **window);
void geom_bdry_e1(geometry_t *geom, geometry_t **window);
void geom_bdry_e2(geometry_t *geom, geometry_t **window);
void geom_bdry_he1(geometry_t *geom, geometry_t **window);
void geom_bdry_he2(geometry_t *geom, geometry_t **window);
void geom_m2s(geometry_t *geom, geometry_t **window);
void smooth(master_t *master, double *A, int *ctp, int nctp);
void smooth3(master_t *master, double *A, int *ctp, int nctp);
void print_zoom_bdry(char *tag, int c, int cnr, int ce, int cn, int cne,
		     double x, double y);

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* zoom_shift       : shifts zoom cell faces to cell centers         */
/* global_interp    : interpolates global zoomed grid centers        */
/* global_interp_e1 : interpolates global zoomed grid e1 faces       */
/* global_interp_e2 : interpolates global zoomed grid e2 faces       */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Routine to shift zoomed cell faces to cell centers                */
/*-------------------------------------------------------------------*/
void zoom_shift(geometry_t *geom, /* Sparse global geometery */
                geometry_t *window  /* Processing window */
  )
{
  int c, cc;
  int zc;
  int zs1 = (int)(window->zmfe1 / 2);
  int zs2 = (int)(window->zmfe2 / 2);
  int c1, cs = 0, cb = 0;

  if (window->zoomf == 1) {
    window->nzin = window->nzinS = 0;
    window->zoomc = i_alloc_1d(window->enonS + 1);
    for (cc = 1; cc <= window->enonS; cc++)
      window->zoomc[cc] = (ZC | ZE1 | ZE2);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Reset the lateral maps for ghost cells and set the zoom codes */
  /* for this window.  */
  reset_zoom_ghosts(geom, window, window->w3_t, window->b3_t);
  window->zoomc = i_alloc_1d(window->enon + 1);
  for (cc = 1; cc <= window->enon; cc++)
    window->zoomc[cc] = ZN;
  set_window_codes(geom, window, window->w3_e1, window->b3_e1, ZE1);
  set_window_codes(geom, window, window->w3_e2, window->b3_e2, ZE2);
  set_window_codes(geom, window, window->w3_t, window->b3_t, ZC);

  /*-----------------------------------------------------------------*/
  /* Get the correct maps from ghost cells to auxiliary cell centers */
  /* (all ghost cells map to zoom centers).  */
  for (cc = 1; cc <= window->enon; cc++) {  /* cc = local sparse coord.  */
    c = window->wsa[cc];        /* c = global sparse coord.  */
    if (!geom->fm[c].wn) {
      /* Left edge : auxiliary cells still map zoomf cells into the */
      /* interior. Reduce this by zoomf/2.  */
      c1 = window->xp1[cc];     /* c1 = cell to east of c */
      /* If cc is west self pointing (ghost cell on west edge) and */
      /* the cell to the east is not a zoom center, then map the */
      /* ghost cell to the cell center (zs global xm1's to the west */
      /* of c1).  */
      if (window->xm1[cc] == cc && !(window->zoomc[c1] & ZC)) {
        c1 = window->wsa[c1];
        for (zc = 1; zc <= zs1; zc++)
          c1 = geom->xm1[c1];
        window->xp1[cc] = window->fm[c1].sc;
      }
      /* Right edge : auxiliary cells still map one cell into the */
      /* interior. Increase this by zoomf/2.  */
      c1 = window->xm1[cc];
      if (window->xp1[cc] == cc && !(window->zoomc[c1] & ZC)) {
        c1 = window->wsa[c1];
        for (zc = 1; zc <= zs1; zc++)
          c1 = geom->xp1[c1];
        window->xm1[cc] = window->fm[c].sc;
      }
      c1 = window->yp1[cc];
      if (window->ym1[cc] == cc && !(window->zoomc[c1] & ZC)) {
        c1 = window->wsa[c1];
        for (zc = 1; zc <= zs2; zc++)
          c1 = geom->ym1[c1];
        window->yp1[cc] = window->fm[c1].sc;
      }
      c1 = window->ym1[cc];
      if (window->yp1[cc] == cc && !(window->zoomc[c1] & ZC)) {
        c1 = window->wsa[c1];
        for (zc = 1; zc <= zs2; zc++)
          c1 = geom->yp1[c1];
        window->ym1[cc] = window->fm[c].sc;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Shift the e1 arrays */
  for (cc = 1; cc <= window->n2_e1; cc++) {
    c = window->w2_e1[cc];
    c1 = geom->xp1[window->wsa[c]];
    if (cc <= window->x2_e1) {
      cs = window->sur_e1[cc];
      cs = geom->xp1[window->wsa[cs]];
      cb = window->bot_e1[cc];
      cb = geom->xp1[window->wsa[cb]];
    }
    /* Note : boundary cells are not shifted since the OBC arrays    */
    /* were shifted when created from the global arrays, and the     */
    /* cells to process boundary cells were mapped directly from the */
    /* (shifted) OBC arrays.                                         */
    if (cc <= window->v2_e1 || cc > window->b2_e1) {
      for (zc = 1; zc <= zs1; zc++) {
	window->w2_e1[cc] = window->fm[c1].sc;
	if (cc <= window->x2_e1) {
	  window->sur_e1[cc] = window->fm[cs].sc;
	  window->bot_e1[cc] = window->fm[cb].sc;
	  cs = geom->xp1[cs];
	  cb = geom->xp1[cb];
	}
	c1 = geom->xp1[c1];
      }
    }    
  }
  for (cc = 1; cc <= window->n3_e1; cc++) {
    c = window->w3_e1[cc];
    c1 = geom->xp1[window->wsa[c]];
    if (cc <= window->v3_e1 || cc > window->b3_e1) {
      for (zc = 1; zc <= zs1; zc++) {
	window->w3_e1[cc] = window->fm[c1].sc;
	c1 = geom->xp1[c1];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Shift the e2 arrays */
  for (cc = 1; cc <= window->n2_e2; cc++) {
    c = window->w2_e2[cc];
    c1 = geom->yp1[window->wsa[c]];
    if (cc <= window->x2_e2) {
      cs = window->sur_e2[cc];
      cs = geom->yp1[window->wsa[cs]];
      cb = window->bot_e2[cc];
      cb = geom->yp1[window->wsa[cb]];
    }
    if (cc <= window->v2_e2 || cc > window->b2_e2) {
      for (zc = 1; zc <= zs2; zc++) {
	window->w2_e2[cc] = window->fm[c1].sc;
	if (cc <= window->x2_e2) {
	  window->sur_e2[cc] = window->fm[cs].sc;
	  window->bot_e2[cc] = window->fm[cb].sc;
	  cs = geom->yp1[cs];
	  cb = geom->yp1[cb];
	}
	c1 = geom->yp1[c1];
      }
    }
  }
  for (cc = 1; cc <= window->n3_e2; cc++) {
    c = window->w3_e2[cc];
    c1 = geom->yp1[window->wsa[c]];
    if (cc <= window->v3_e2 || cc > window->b3_e2) {
      for (zc = 1; zc <= zs2; zc++) {
	window->w3_e2[cc] = window->fm[c1].sc;
	c1 = geom->yp1[c1];
      }
    }
  }
}

/* END zoom_shift()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the lateral maps for ghost cells. Ghost cells    */
/* may not map to zoom cell centers at this stage, since the maps    */
/* are created for all cells in the zoom window and ghost cells map  */
/* zoomf cells into the interior for east/south boundaries and one   */
/* cell into the interior for west/north boundaries rather than the  */
/* 1+zoomf/2 cells needed to map to a cell center.                   */
/*-------------------------------------------------------------------*/
void reset_zoom_ghosts(geometry_t *geom,  /* Sparse global geometery */
                       geometry_t *window,  /* Processing window */
                       int *vec,  /* Cells to process vector */
                       int nvec /* Size of vec */
  )
{
  int c, cc, c1;

  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    c1 = window->xm1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->xp1[c1] = c;
    c1 = window->xp1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->xm1[c1] = c;
    c1 = window->ym1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->yp1[c1] = c;
    c1 = window->yp1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->ym1[c1] = c;
    c1 = window->xmyp1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->xpym1[c1] = c;
    c1 = window->xpym1[c];
    if (!geom->fm[window->wsa[c1]].wn)
      window->xmyp1[c1] = c;
  }
}

/* END reset_zoom_ghosts()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computs the scaling factor to apply to coarse fluxes to derive    */
/* the fluxes for fine cells. The scaling is based on the fine       */
/* cell's width comparsed to the coarse cell. Note that if some of   */
/* the coarse cell is associated with solid boundaries then the      */
/* scaling must be larger to take this into account.                 */
/*-------------------------------------------------------------------*/
void set_flux_scale(geometry_t *geom,
		    geometry_t *window)
{
  win_priv_t *wincon = window->wincon;
  int c, cs, cg, cp, cm, cc, cps, cms;
  int z1, z2;
  int zs1 = (int)(window->zmfe1 / 2);
  int zs2 = (int)(window->zmfe2 / 2);

  if (window->zoomf == 1)
    return;

  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    cg = cp = cm = window->wsa[c];
    wincon->w4[c] = 0.0;

    for (z1 = 1; z1 <= zs1; z1++)
      cg = geom->xm1[cg];
    cs = geom->m2d[cg];
    cg = geom->xm1[cg]; /* Fine cell adjacent to coarse */
    if (geom->xm1[cg] != cg)  /* Check for solid edge */
      wincon->w4[c] += geom->h2au1[cs];

    for (z2 = 1; z2 <= zs2; z2++) {
      cp = geom->yp1[cp];
      cm = geom->ym1[cm];
      for (z1 = 1; z1 <= zs1; z1++) {
	cp = geom->xm1[cp];
	cm = geom->xm1[cm];
      }
      cps = geom->m2d[cp];
      cms = geom->m2d[cm];
      cp = geom->xm1[cp];   /* Fine grid cell adjacent to coarse cell */
      cm = geom->xm1[cm];
      if (geom->xm1[cp] != cp)
	wincon->w4[c] += geom->h2au1[cps];
      if (geom->xm1[cm] != cm)
	wincon->w4[c] += geom->h2au1[cms];
    }
  }

  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    cg = cp = cm = window->wsa[c];
    wincon->w5[c] = 0.0;

    for (z2 = 1; z2 <= zs2; z2++)
      cg = geom->ym1[cg];
    cs = geom->m2d[cg];
    cg = geom->ym1[cg]; /* Fine cell adjacent to coarse */
    if (geom->ym1[cg] != cg)  /* Check for solid edge */
      wincon->w5[c] += geom->h1au2[cs];

    for (z1 = 1; z1 <= zs1; z1++) {
      cp = geom->xp1[cp];
      cm = geom->xm1[cm];
      for (z2 = 1; z2 <= zs1; z2++) {
	cp = geom->ym1[cp];
	cm = geom->ym1[cm];
      }
      cps = geom->m2d[cp];
      cms = geom->m2d[cm];
      cp = geom->ym1[cp];   /* Fine grid cell adjacent to coarse cell */
      cm = geom->ym1[cm];
      if (geom->ym1[cp] != cp)
	wincon->w5[c] += geom->h1au2[cps];
      if (geom->ym1[cm] != cm)
	wincon->w5[c] += geom->h1au2[cms];
    }
  }
}

/* END set_flux_scale()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the zoom codes within a window. These are defined  */
/* over all wet, auxiliary and ghost cells for the full 3D array.    */
/*-------------------------------------------------------------------*/
void set_window_codes(geometry_t *geom, /* Sparse global geometery */
                      geometry_t *window, /* Processing window */
                      int *vec, /* Cells to process vector */
                      int nvec, /* Size of vec */
                      int zc    /* Zoom code */
  )
{
  int c, cc, c1;
  int *wm;                      /* Window / ghost cell mask */
  int *mask;

  wm = i_alloc_1d(geom->sgsiz);
  mask = i_alloc_1d(geom->sgsiz);
  for (c = 1; c <= geom->sgnum; c++)
    wm[c] = geom->fm[c].wn;
  if (zc == ZE1) {
    /* Western edge ghost cells for e1 cells to process. These cells */
    /* are actually wet and defined to be in a window but due to the */
    /* stagger for e1 velocities at western boundaries they are */
    /* ghost cells for u1. These are identified if the west[c] is */
    /* non-zero (i.e. a wet cell in a window) and west[west[c]] is */
    /* zero (i.e. a ghost cell).  */
    for (cc = 1; cc <= geom->v3_e1; cc++) {
      c = geom->w3_e1[cc];      /* Global cell to process */
      c1 = geom->xm1[c];        /* Global cell to the west of c */
      if (geom->fm[c1].wn && !(geom->fm[geom->xm1[c1]].wn))
        wm[c1] = 0;
    }
  }
  if (zc == ZE2) {
    /* Southern edge ghost cells for e2 cells to process.  */
    for (cc = 1; cc <= geom->v3_e2; cc++) {
      c = geom->w3_e2[cc];      /* Global cell to process */
      c1 = geom->ym1[c];        /* Global cell to the south of c */
      if (geom->fm[c1].wn && !(geom->fm[geom->ym1[c1]].wn))
        wm[c1] = 0;
    }
  }

  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];

    /* Set the zoom code at the wet cell to process */
    if (window->zoomc[c] == ZN)
      window->zoomc[c] = zc;
    else
      window->zoomc[c] |= zc;

    /* Set the zoom code at auxiliary cells */
    if (wm[window->wsa[window->xp1[c]]] != window->wn) {
      c1 = c;
      while (c1 != window->xp1[c1]) {
        c1 = window->xp1[c1];
        window->zoomc[c1] |= zc;
      }
    }
    if (wm[window->wsa[window->xm1[c]]] != window->wn) {
      c1 = c;
      while (c1 != window->xm1[c1]) {
        c1 = window->xm1[c1];
        window->zoomc[c1] |= zc;
      }
    }
    if (wm[window->wsa[window->yp1[c]]] != window->wn) {
      c1 = c;
      while (c1 != window->yp1[c1]) {
        c1 = window->yp1[c1];
        window->zoomc[c1] |= zc;
      }
    }
    if (wm[window->wsa[window->ym1[c]]] != window->wn) {
      c1 = c;
      while (c1 != window->ym1[c1]) {
        c1 = window->ym1[c1];
        window->zoomc[c1] |= zc;
      }
    }
    if (wm[window->wsa[window->xm1[window->ym1[c]]]] != window->wn) {
      c1 = c;
      while (c1 != window->xm1[window->ym1[c1]]) {
        c1 = window->xm1[window->ym1[c1]];
        window->zoomc[c1] |= zc;
      }
    }
    if (wm[window->wsa[window->xm1[window->yp1[c]]]] != window->wn) {
      c1 = c;
      memset(mask, 0, geom->sgsiz * sizeof(int));
      while (c1 != window->xm1[window->yp1[c1]] && !(mask[c1])) {
        c1 = window->xm1[window->yp1[c1]];
        window->zoomc[c1] |= zc;
	mask[c1] = 1;
      }
    }
    if (wm[window->wsa[window->xp1[window->ym1[c]]]] != window->wn) {
      c1 = c;
      memset(mask, 0, geom->sgsiz * sizeof(int));
      while (c1 != window->xp1[window->ym1[c1]] && !(mask[c1])) {
        c1 = window->xp1[window->ym1[c1]];
        window->zoomc[c1] |= zc;
	mask[c1] = 1;
      }
    }
    if (wm[window->wsa[window->xp1[window->yp1[c]]]] != window->wn) {
      c1 = c;
      memset(mask, 0, geom->sgsiz * sizeof(int));
      while (c1 != window->xp1[window->yp1[c1]] && !(mask[c1])) {
        c1 = window->xp1[window->yp1[c1]];
        window->zoomc[c1] |= zc;
	mask[c1] = 1;
      }
    }
  }
  i_free_1d(wm);
  i_free_1d(mask);
}

/* END set_window_codes()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate a variable over the zoomed grid using      */
/* bilinear interpolation.                                           */
/*-------------------------------------------------------------------*/
void global_interp(geometry_t *geom,  /* Sparse geometry */
                   double *A,   /* Array to interpolate */
                   int nz       /* Number of interpolation cells */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int ix;                       /* Weight counter */
  double wt[4], vals[4];        /* Weights, corner values */

  if (geom->zoom & PRECOND) {
    for (cc = 1; cc <= nz; cc++) {
      c = geom->zin[cc].c;
      A[c] = A[geom->zin[cc].cnr];
      /*fill_zoom_cell(geom, A, c, geom->zin[cc].ce, geom->zin[cc].cn);*/
    }
    return;
  }

  for (cc = 1; cc <= nz; cc++) {
    /* Get the weights */
    c = geom->zin[cc].c;
    wt[0] = (1.0 - geom->zin[cc].x) * (1.0 - geom->zin[cc].y);
    wt[1] = (1.0 - geom->zin[cc].x) * geom->zin[cc].y;
    wt[2] = geom->zin[cc].x * (1.0 - geom->zin[cc].y);
    wt[3] = geom->zin[cc].x * geom->zin[cc].y;
    /* Get the four integral neighbours of the non-integral point */
    vals[0] = A[geom->zin[cc].cnr];
    vals[1] = A[geom->zin[cc].cn];
    vals[2] = A[geom->zin[cc].ce];
    vals[3] = A[geom->zin[cc].cne];
    /* Do the bilinear interpolation */
    A[c] = 0.0;
    /*if(A==master->u1av&&c==448)printf("%f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);*/
    for (ix = 0; ix < 4; ix++)
      A[c] += wt[ix] * vals[ix];
  }
}

/* END global_interp()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate a variable on e1 faces over the zoomed     */
/* grid using bilinear interpolation.                                */
/*-------------------------------------------------------------------*/
void global_interp_e1(geometry_t *geom, /* Sparse geometry */
                      double *A,  /* Array to interpolate */
                      int nz    /* Number of interpolation cells */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int ix;                       /* Weight counter */
  double wt[4], vals[4];        /* Weights, corner values */

  if (geom->zoom & PRECOND) {
    for (cc = 1; cc <= nz; cc++) {
      c = geom->zin[cc].c;
      A[c] = A[geom->zin[cc].cnr];
      /*fill_zoom_cell(geom, A, c, geom->zin[cc].ce, geom->zin[cc].cn);*/
    }
    return;
  }

  for (cc = 1; cc <= nz; cc++) {
    /* Get the weights */
    c = geom->zine1[cc].c;

    wt[0] = (1.0 - geom->zine1[cc].x) * (1.0 - geom->zine1[cc].y);
    wt[1] = (1.0 - geom->zine1[cc].x) * geom->zine1[cc].y;
    wt[2] = geom->zine1[cc].x * (1.0 - geom->zine1[cc].y);
    wt[3] = geom->zine1[cc].x * geom->zine1[cc].y;
    /* Get the four integral neighbours of the non-integral point */
    vals[0] = A[geom->zine1[cc].cnr];
    vals[1] = A[geom->zine1[cc].cn];
    vals[2] = A[geom->zine1[cc].ce];
    vals[3] = A[geom->zine1[cc].cne];

    /* Do the bilinear interpolation */
    A[c] = 0.0;
    for (ix = 0; ix < 4; ix++)
      A[c] += wt[ix] * vals[ix];
  }
}

/* END global_interp_e1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate a variable on e2 faces over the zoomed     */
/* grid using bilinear interpolation.                                */
/*-------------------------------------------------------------------*/
void global_interp_e2(geometry_t *geom, /* Sparse geometry */
                      double *A,  /* Array to interpolate */
                      int nz    /* Number of interpolation cells */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int ix;                       /* Weight counter */
  double wt[4], vals[4];        /* Weights, corner values */

  if (geom->zoom & PRECOND) {
    for (cc = 1; cc <= nz; cc++) {
      c = geom->zin[cc].c;
      A[c] = A[geom->zin[cc].cnr];
      /*fill_zoom_cell(geom, A, c, geom->zin[cc].ce, geom->zin[cc].cn);*/
    }
    return;
  }

  for (cc = 1; cc <= nz; cc++) {
    /* Get the weights */
    c = geom->zine2[cc].c;

    wt[0] = (1.0 - geom->zine2[cc].x) * (1.0 - geom->zine2[cc].y);
    wt[1] = (1.0 - geom->zine2[cc].x) * geom->zine2[cc].y;
    wt[2] = geom->zine2[cc].x * (1.0 - geom->zine2[cc].y);
    wt[3] = geom->zine2[cc].x * geom->zine2[cc].y;
    /* Get the four integral neighbours of the non-integral point */

    vals[0] = A[geom->zine2[cc].cnr];
    vals[1] = A[geom->zine2[cc].cn];
    vals[2] = A[geom->zine2[cc].ce];
    vals[3] = A[geom->zine2[cc].cne];
    /* Do the bilinear interpolation */
     A[c] = 0.0;
    for (ix = 0; ix < 4; ix++)
      A[c] += wt[ix] * vals[ix];
  }
}

/* END global_interp_e2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate a variable on e1 faces over the zoomed     */
/* grid using linear interpolation in e1 followed by e2 directions.  */
/*-------------------------------------------------------------------*/
void global_sinterp_e1(geometry_t *geom, /* Sparse geometry          */
                      double *A,    /* Array to interpolate          */
                      int nzp,      /* Number of interpolation cells */
                      int nzc       /* Number of interpolation cells */
  )
{
  int c, cc;                    /* Sparse coordinate / counter       */
  int ix;                       /* Weight counter                    */
  double wt[2], vals[2];        /* Weights, corner values            */

  /* Interpolate in the parallel (e2) direction                      */
  for (cc = 1; cc <= nzp; cc++) {
    /* Get the weights */
    c = geom->zipe1[cc].c;
    wt[0] = 1.0 - geom->zipe1[cc].x;
    wt[1] = geom->zipe1[cc].x;
    vals[0] = A[geom->zipe1[cc].cnr];
    vals[1] = A[geom->zipe1[cc].ce];
    /* Do the linear interpolation                                   */
    A[c] = 0.0;
    for (ix = 0; ix < 2; ix++)
      A[c] += wt[ix] * vals[ix];
  }
  /* Interpolate in the cross (e1) direction                         */
  for (cc = 1; cc <= nzc; cc++) {
    /* Get the weights */
    c = geom->zice1[cc].c;
    wt[0] = 1.0 - geom->zice1[cc].y;
    wt[1] = geom->zice1[cc].y;
    vals[0] = A[geom->zice1[cc].cnr];
    vals[1] = A[geom->zice1[cc].cn];
    /* Do the linear interpolation                                   */
    A[c] = 0.0;
    for (ix = 0; ix < 2; ix++)
      A[c] += wt[ix] * vals[ix];
  }
}

/* END global_sinterp_e1()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate a variable on e2 faces over the zoomed     */
/* grid using linear interpolation in e2 followed by e1 directions.  */
/*-------------------------------------------------------------------*/
void global_sinterp_e2(geometry_t *geom, /* Sparse geometry          */
                      double *A,    /* Array to interpolate          */
                      int nzp,      /* Number of interpolation cells */
                      int nzc       /* Number of interpolation cells */
  )
{
  int c, cc;                    /* Sparse coordinate / counter       */
  int ix;                       /* Weight counter                    */
  double wt[2], vals[2];        /* Weights, corner values            */

  /* Interpolate in the parallel (e2) direction                      */
  for (cc = 1; cc <= nzp; cc++) {
    /* Get the weights */
    c = geom->zipe2[cc].c;
    wt[0] = 1.0 - geom->zipe2[cc].y;
    wt[1] = geom->zipe2[cc].y ;
    vals[0] = A[geom->zipe2[cc].cnr];
    vals[1] = A[geom->zipe2[cc].cn];
    /* Do the linear interpolation                                   */
    A[c] = 0.0;
    for (ix = 0; ix < 2; ix++)
      A[c] += wt[ix] * vals[ix];
  }
  /* Interpolate in the cross (e1) direction                         */
  for (cc = 1; cc <= nzc; cc++) {
    /* Get the weights */
    c = geom->zice2[cc].c;
    wt[0] = 1.0 - geom->zice2[cc].x;
    wt[1] = geom->zice2[cc].x ;
    vals[0] = A[geom->zice2[cc].cnr];
    vals[1] = A[geom->zice2[cc].ce];
    /* Do the linear interpolation                                   */
    A[c] = 0.0;
    for (ix = 0; ix < 2; ix++)
      A[c] += wt[ix] * vals[ix];
  }
}

/* END global_sinterp_e2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate 2D velocity variables                      */
/*-------------------------------------------------------------------*/
void global_interp_vel2d(master_t *master, /* Master data            */
			 geometry_t *geom  /* Global geometry        */
			 )
{
  int smoothf = 0;

  /* Interpolate the updated velocities onto the coarse grid         */
  global_interp_e1(geom, master->nu1av, geom->nzine1S);
  global_interp_e2(geom, master->nu2av, geom->nzine2S);

  /*
  global_sinterp_e1(geom, master->nu1av, geom->nipe1S, geom->nice1S);
  global_sinterp_e2(geom, master->nu2av, geom->nipe2S, geom->nice2S);
  */
  memcpy(master->d3, master->u1av, geom->sgsizS * sizeof(double));
  memcpy(master->d4, master->u2av, geom->sgsizS * sizeof(double));

  /* Perform smoothing on interpolated velocities                    */
  /* Note; the whole window must be transferred to the master if     */
  /* this is to be performed.                                        */
  if (smoothf) {
    int c, cc, cv;
    memcpy(master->d6, master->nu1av, geom->sgsizS * sizeof(double));
    for (cc = 1; cc <= geom->nipe1S; cc++) {
      c = geom->zipe1[cc].c;
      cv = geom->zipe1[cc].ce;
      if (geom->zipe1[cc].x >= 0.0)
	master->nu1av[cv] = shuman(geom, master->d6, c);
    }
    memcpy(master->d6, master->nu2av, geom->sgsizS * sizeof(double));
    for (cc = 1; cc <= geom->nipe1S; cc++) {
      c = geom->zipe1[cc].c;
      cv = geom->zipe1[cc].cn;
      if (geom->zipe1[cc].y >= 0.0) {
	master->nu2av[cv] = shuman(geom, master->d6, c);
      }
    }
  }
}

/* END global_interp_vel2d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate eta at on the master coarse grid           */
/*-------------------------------------------------------------------*/
void global_interp_pre2d(master_t *master, /* Master data            */
			 geometry_t *geom, /* Global geometry        */
			 int ic,
			 double dt2d
			 )
{
  global_interp_flux(master, geom, ic, dt2d);
}

/* END global_interp_pre2d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate total depth                                  */
/*-------------------------------------------------------------------*/
void global_interp_post2d(master_t *master, /* Master data           */
			  geometry_t *geom  /* Global geometry       */
			  )
{
  int smoothf = 1;

  global_interp_depth(master, geom);

  /* Perform smoothing on interpolated velocities                    */
  /* Note; the whole window must be transferred to the master if     */
  /* this is to be performed.                                        */
  if (smoothf) {
    int c, cc, cv;
    memcpy(master->d6, master->u1av, geom->sgsizS * sizeof(double));
    for (cc = 1; cc <= geom->nipe1S; cc++) {
      c = geom->zipe1[cc].c;
      cv = geom->zipe1[cc].ce;
      if (geom->zipe1[cc].x >= 0.0)
	master->u1av[c] = shuman(geom, master->d6, c);
    }
    memcpy(master->d6, master->u2av, geom->sgsizS * sizeof(double));
    for (cc = 1; cc <= geom->nipe1S; cc++) {
      c = geom->zipe1[cc].c;
      cv = geom->zipe1[cc].cn;
      if (geom->zipe1[cc].y >= 0.0) {
	master->u2av[c] = shuman(geom, master->d6, c);
      }
    }
  }
}

/* END global_interp_post2d()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to step forward in time.                                  */
/*-------------------------------------------------------------------*/
void global_interp_pre(master_t *master, /* Master data              */
		       geometry_t *geom  /* Global geometry          */
		       )
{
  int c, cc;

  /* Step forward in time                                            */
  for (cc = 1; cc <= geom->nzine1; cc++) {
    c = geom->zine1[cc].c;
    master->u1b[c] = master->u1[c];
  }
  for (cc = 1; cc <= geom->nzine2; cc++) {
    c = geom->zine2[cc].c;
    master->u2b[c] = master->u2[c];
  }
}

/* END global_interp_pre()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to interpolate vertical velocity variables and step       */
/* in time.                                                          */
/*-------------------------------------------------------------------*/
void global_interp_post(master_t *master, /* Master data             */
			geometry_t *geom  /* Global geometry         */
			)
{

  /* Interpolate vertical velocity variables on the master           */
  global_interp(geom, master->w, geom->nzin);
  global_interp(geom, master->wtop, geom->nzinS);
  global_interp(geom, master->wbot, geom->nzinS);
}

/* END global_interp_post()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the 3D interpolated e1 velocities on the zoomed */
/* grid.                                                             */
/*-------------------------------------------------------------------*/
void interp_adjust_e1(master_t *master, /* Master data               */
		      geometry_t *geom  /* Global geometry           */
  )
{
  int c, cc, cs, cb, zm1;       /* Sparse coordinate / counter       */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *adjust;               /* Velocity adjustment               */
  double *maxeta;                             
  double depth;
  int nzS = geom->nzine1S;
  int nz = geom->nzine1;

  /*-----------------------------------------------------------------*/
  /* Adjust the e1 velocities.                                       */
  /* Set pointers and initialise                                     */
  sum = d_alloc_1d(geom->sgsizS);
  adjust = d_alloc_1d(geom->sgsizS);
  maxeta = master->d1;
  memset(sum, 0, geom->sgsizS * sizeof(double));
  memset(adjust, 0, geom->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column               */
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine1[cc].c;
    cs = geom->m2d[c];
    sum[cs] += master->u1[c] * master->dzu1[c];
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine1[cc].c;
    cs = geom->m2d[c];
    depth = maxeta[cs] - geom->botzu1[cs];
    adjust[cs] = (master->u1flux[cs] - 
		  sum[cs] * master->dt * geom->h2au1[cs]) / 
      (depth * master->dt * geom->h2au1[cs]);
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine1[cc].c;
    cs = geom->m2d[c];
    master->u1[c] += adjust[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Fluxes through the e1 faces                                     */
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine1[cc].c;
    cs = geom->m2d[c];
    master->u1flux3d[c] = master->u1[c] * master->dzu1[c] * geom->h2au1[cs];
  }
  /* Zero fluxes above the surface                                   */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine1[cc].c;
    cb = geom->zine1[cc].cb;
    cs = geom->m2d[c];
    zm1 = geom->zm1[c];
    /* Loop down the water column                                    */
    while (c != zm1 && geom->gridz[c] > maxeta[cs]) {
      master->u1flux3d[c] = 0.0;
      c = zm1;
      zm1 = geom->zm1[c];
    }
    master->u1bot[cs] = master->u1b[cb] - master->u1avb[cs];
  }
  d_free_1d(sum);
  d_free_1d(adjust);
}

/* END interp_adjust_e1()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the 3D interpolated e2 velocities on the zoomed */
/* grid.                                                             */
/*-------------------------------------------------------------------*/
void interp_adjust_e2(master_t *master, /* Master data               */
		      geometry_t *geom  /* Global geometry           */
  )
{
  int c, cc, cs, cb, zm1;       /* Sparse coordinate / counter       */
  double *sum;                  /* Vertically integrated 3D velocity */
  double *adjust;               /* Velocity adjustment               */
  double *maxeta;                             
  double depth;
  int nzS = geom->nzine2S;
  int nz = geom->nzine2;

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e2 faces.                                 */
  /* Set pointers and initialise                                     */
  sum = d_alloc_1d(geom->sgsizS);
  adjust = d_alloc_1d(geom->sgsizS);
  maxeta = master->d2;
  memset(sum, 0, geom->sgsizS * sizeof(double));
  memset(adjust, 0, geom->sgsizS * sizeof(double));

  /* Integrate the velocities through the water column               */
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    if (c != geom->zm1[c])
      sum[cs] += master->u2[c] * master->dzu2[c];
  }

  /* Calculate the velocity adjustment                               */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    depth = maxeta[cs] - geom->botzu2[cs];
    adjust[cs] = (master->u2flux[cs] - 
		  sum[cs] * master->dt * geom->h1au2[cs]) / 
      (depth * master->dt * geom->h1au2[cs]);
  }

  /* Adjust the 3D velocities                                        */
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    master->u2[c] += adjust[cs];
  }

  /*-----------------------------------------------------------------*/
  /* Fluxes through the e2 faces                                     */
  memset(sum, 0, geom->sgsizS * sizeof(double));
  for (cc = 1; cc <= nz; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    master->u2flux3d[c] = master->u2[c] * master->dzu2[c] * geom->h1au2[cs];
    sum[cs]+=master->dzu2[c];
  }

  /* Zero fluxes above the surface                                   */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    cb = geom->zine2[cc].cb;
    zm1 = geom->zm1[c];
    /* Loop down the water column                                    */
    while (c != zm1 && geom->gridz[c] > maxeta[cs]) {
      master->u2flux3d[c] = 0.0;
      c = zm1;
      zm1 = geom->zm1[c];
    }
    master->u2bot[cs] = master->u2b[cb] - master->u2avb[cs];
  }
  d_free_1d(sum);
  d_free_1d(adjust);
}

/* END interp_adjust_e2()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the cell thickness at e1 faces.              */
/*-------------------------------------------------------------------*/
void interp_dz_at_u1(master_t *master,  /* Master data               */
		     geometry_t *geom   /* Global geometry           */
  )
{
  int c, cc, cs, c3, cb, zm1;         /* Sparse coordinate / counter */
  double *maxeta;              
  int *sur_e1;
  double top, bot;
  int nzS = geom->nzine1S;

  /*-----------------------------------------------------------------*/
  /* Calculate the cell spacings.                                    */
  /* Set pointers and initialise.                                    */
  maxeta = master->d1;
  sur_e1 = i_alloc_1d(nzS + 1);

  /* Get the elevation at the e1 face                                */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine1[cc].c;
    cs = geom->m2d[c];
    maxeta[cs] = max(master->eta[c], master->eta[geom->xm1[c]]);
    maxeta[cs] = max(maxeta[cs], geom->botzu1[c]);
    zm1 = geom->zm1[c];
    /* Loop down the water column                                    */
    while (c != zm1 && geom->gridz[c] > maxeta[cs]) {
      c = zm1;
      zm1 = geom->zm1[c];
    }
    sur_e1[cc] = c;
  }

  for (cc = 1; cc <= nzS; cc++) {
    /* Get the 2D and 3D sparse coordinates of the surface         */
    c = c3 = sur_e1[cc];
    cs = geom->m2d[c];          /* 2D cell corresponding to c      */
    cb = geom->zine1[cc].cb;    /* 3D bottom coordinate            */

    /* Set the cell thickness from the surface to the layer above  */
    /* the bottom.                                                 */
    top = maxeta[cs];
    while (c != cb) {
      bot = geom->gridz[c];
      master->dzu1[c] = top - bot;
      top = bot;
      c = geom->zm1[c];
    }
    /*
    zm1 = geom->zm1[c];
    while (zm1 != geom->zm1[zm1]) {
      bot = geom->gridz[c];
      master->dzu1[c] = top - bot;
      top = bot;
      c = zm1;
      zm1 = geom->zm1[c];
    }
    */
    /* Set the cell thickness at the bottom                          */
    master->dzu1[c] = top - geom->botzu1[cs];

    /* Set the cell thickness above the old surface equal to zero. */
    c = c3;
    top = master->dzu1[c];
    while (c != cs) {
      c = geom->zp1[c];
      master->dzu1[c] = top;
    }
  }
  i_free_1d(sur_e1);
}

/* END interp_dz_at_u1()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the cell thickness at e1 faces.              */
/*-------------------------------------------------------------------*/
void interp_dz_at_u2(master_t *master,  /* Master data               */
		     geometry_t *geom   /* Global geometry           */
  )
{
  int c, cc, cs, c3, cb, zm1;         /* Sparse coordinate / counter */
  double *maxeta;              
  int *sur_e2;
  double top, bot;
  int nzS = geom->nzine2S;

  /*-----------------------------------------------------------------*/
  /* Calculate the cell spacings.                                    */
  /* Set pointers and initialise.                                    */
  maxeta = master->d2;
  sur_e2 = i_alloc_1d(nzS + 1);

  /* Get the elevation at the e2 face                                */
  for (cc = 1; cc <= nzS; cc++) {
    c = geom->zine2[cc].c;
    cs = geom->m2d[c];
    maxeta[cs] = max(master->eta[c], master->eta[geom->ym1[c]]);
    maxeta[cs] = max(maxeta[cs], geom->botzu2[c]);
    zm1 = geom->zm1[c];
    /* Loop down the water column                                    */
    while (c != zm1 && geom->gridz[c] > maxeta[cs]) {
      c = zm1;
      zm1 = geom->zm1[c];
    }
    sur_e2[cc] = c;
  }

  for (cc = 1; cc <= nzS; cc++) {
    /* Get the 2D and 3D sparse coordinates of the surface         */
    c = c3 = sur_e2[cc];
    cs = geom->m2d[c];          /* 2D cell corresponding to c      */
    cb = geom->zine2[cc].cb;    /* 3D bottom coordinate            */

    /* Set the cell thickness from the surface to the layer above  */
    /* the bottom.                                                 */
    top = maxeta[cs];
    while (c != cb) {
      bot = geom->gridz[c];
      master->dzu2[c] = top - bot;
      top = bot;
      c = geom->zm1[c];
    }

    /* Set the cell thickness at the bottom                          */
    master->dzu2[c] = top - geom->botzu2[cs];

    /* Set the cell thickness above the old surface equal to zero. */
    c = c3;
    top = master->dzu2[c];
    while (c != cs) {
      c = geom->zp1[c];
      master->dzu2[c] = top;
    }
  }
  i_free_1d(sur_e2);
}

/* END interp_dz_at_u2()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the depths at the e1 and e2 faces.           */
/*-------------------------------------------------------------------*/
void global_interp_depth(master_t *master,    /* Master data         */
			 geometry_t *geom     /* Global geometry     */
			 )
{
  int c, cc;            /* Sparse coordinate / counter               */
  double *tzp;          /* Pointer to sutface elevation              */

  /*-----------------------------------------------------------------*/
  /* Set pointers                                                    */
  tzp = (master->nonlinear && !master->sigma) ? master->eta : master->topz;

  /*-----------------------------------------------------------------*/
  /* Get the depths at the e1 and e2 faces.                          */
  for (cc = 1; cc <= geom->nzine1S; cc++) {
    c = geom->zine1[cc].c;
    master->depth_e1[c] = max(tzp[geom->xm1[c]], tzp[c]) - geom->botzu1[c];
  }
  for (cc = 1; cc <= geom->nzine2S; cc++) {
    c = geom->zine2[cc].c;
    master->depth_e2[c] = max(tzp[geom->ym1[c]], tzp[c]) - geom->botzu2[c];
  }
}

/* END global_interp_depth()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the interpolated 2D fluxes on the zoomed     */
/* grid.                                                             */
/*-------------------------------------------------------------------*/
void global_interp_flux(master_t *master,    /* Master data          */
			geometry_t *geom,    /* Global geometry      */
			int ic,              /* 2D mode counter      */
			double dt2d          /* 2D leapfrog timestep */
  )
{
  int c, cc;            /* Sparse coordinate / counter               */
  double aconst = 0.1;  /* Constant for Asselin time filer           */
  double *tzp;          /* Pointer to sutface elevation              */
  double *u1av, *u2av;  /* Pointer to unfiltered velocity            */
  double *etab;         /* Pointer to eta(t-1) prior to s2m          */

  /*-----------------------------------------------------------------*/
  /* Set pointers.                                                   */
  /* Note that when s2m transfers are made the master is filled with */
  /* zero values since the cells requiring interpolation on the      */
  /* windows are not processed. These (auxiliary) cells in the       */
  /* windows only receive non-zero values from the master if m2s     */
  /* transfers are made after master interpolation. If s2m transfer  */
  /* fills zeros into the master and the variable is required (e.g.  */
  /* master->etab for calculation of master->detadt) then a dummy    */
  /* must be used (e.g. master->d5).                                 */
  tzp = (master->nonlinear && !master->sigma) ? master->eta : master->topz;
  u1av = master->d3;
  u2av = master->d4;
  etab = master->d5;

  /*-----------------------------------------------------------------*/
  /* Apply the Asselin time fileter to interpolated velocities       */
  for (cc = 1; cc <= geom->nzine1S; cc++) {
    c = geom->zine1[cc].c;
    master->u1av[c] = u1av[c] + 0.5 * aconst *
      (master->u1avb[c] - 2.0 * u1av[c] + master->nu1av[c]);
  }
  for (cc = 1; cc <= geom->nzine2S; cc++) {
    c = geom->zine2[cc].c;
    master->u2av[c] = u2av[c] + 0.5 * aconst *
      (master->u2avb[c] - 2.0 * u2av[c] + master->nu2av[c]);
  }

  /*-----------------------------------------------------------------*/
  /* Interpolate the elevation to the zoomed grid                    */
  global_interp(geom, tzp, geom->nzinS);       /* Bi-linear          */
  /*global_interp_nn(geom, tzp, geom->nzinS);*//* Natural neighbours */
  /* Set eta on the master using continuity                          */
  /*global_eta(geom, dt2d);*/

  /* Smooth the elevation on the zoomed grid.                        */
  /* Note; the whole window must be transferred to the master if     */
  /* this is to be performed.                                        */
  /*
  memcpy(master->d6, master->nu1av, geom->sgsizS * sizeof(double));
  for (cc = 1; cc <= geom->nipe1S; cc++) {
    c = geom->zipe1[cc].c;
    master->eta[c] = shuman(geom, master->d6, c);
  }
  for (cc = 1; cc <= geom->nzinS; cc++) {
    c = geom->zin[cc].c;
    master->eta[c] = shuman(geom, master->d6, c);
  }
  for (cc = 1; cc <= geom->nice1S; cc++) {
    c = geom->zice1[cc].c;
    master->eta[c] = shuman(geom, master->d6, c);
  }
  for (cc = 1; cc <= geom->nice2S; cc++) {
    c = geom->zice2[cc].c;
    master->eta[c] = shuman(geom, master->d6, c);
  }
  */

  /* Calculate the rate of change of elevation                       */
  for (cc = 1; cc <= geom->nzinS; cc++) {
    c = geom->zin[cc].c;
    master->detadt[c] = (master->eta[c] - etab[c]) / dt2d;
  }      
  memcpy(master->etab, master->eta, geom->sgsizS * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Get the fluxes at the e1 and e2 faces and sum the fluxes on odd */
  /* timesteps. Note the time filtered velocity must be used here.   */
  if (((ic + 1) % 2) == 0) {
    for (cc = 1; cc <= geom->nzine1S; cc++) {
      c = geom->zine1[cc].c;
      master->u1flux[c] += master->u1av[c] * master->depth_e1[c] * 
	geom->h2au1[c] * dt2d;
    }
    for (cc = 1; cc <= geom->nzine2S; cc++) {
      c = geom->zine2[cc].c;
      master->u2flux[c] += master->u2av[c] * master->depth_e2[c] * 
	geom->h1au2[c] * dt2d;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Step forward in time at interpolated velocity cells. These      */
  /* velocities at other time levels are required to apply the time  */
  /* filter at the current timestep prior to flux calculation.       */
  /* Also, master->etab is required for calculation of detadt.       */
  for (cc = 1; cc <= geom->nzine1S; cc++) {
    c = geom->zine1[cc].c;
    master->u1avb[c] = master->u1av[c];
    master->u1av[c] = master->nu1av[c];
  }
  for (cc = 1; cc <= geom->nzine2S; cc++) {
    c = geom->zine2[cc].c;
    master->u2avb[c] = master->u2av[c];
    master->u2av[c] = master->nu2av[c];
  }
  for (cc = 1; cc <= geom->nzinS; cc++) {
    c = geom->zin[cc].c;
    master->d5[c] = master->oldeta[c];
    master->oldeta[c] = master->eta[c];
  }
}

/* END global_interp_flux()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the 2D fluxes on the fine grid that are      */
/* required on the coarse grid prior to transfer (and summation) on  */
/* the coarse grid. This must be performed before the NVELOCITY      */
/* transfer is performed in the 1st part (velocity) of the 2D mode;  */
/* 2D fluxes are not transferred and are 'in situ' calculated.       */
/* Note that these fluxes cannot be summed in the coarse grid in     */
/* eta_step because the coarse maps map zoomf cells and maps to the  */
/* inter-zoom fine cell fluxes required for summation are not        */
/* available. The dummy u1flux_z, u2flux_z is used so that fine and  */
/* summmed coarse fluxes are kept seperate.                          */
/*-------------------------------------------------------------------*/
void set_zoom_flux(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon    /* Window constants         */
		   )
{
  double aconst = 0.1;          /* Time filter constant */
  double vfil;
  int c, cc;

  if(wincon->dozoom != NOZOOM && window->zoomf == 1) {

    for (cc = 1; cc <= window->ns2mS; cc++) {
      c = window->s2m[cc];
      /* Calculate the flux at e1 wet and boundary cells             */
      vfil = windat->u1av[c] + 0.5 * aconst *
        (windat->u1avb[c] - 2.0 * windat->u1av[c] + windat->nu1av[c]);
      windat->u1flux_z[c] = vfil * windat->depth_e1[c] * 
	window->h2au1[c] * wincon->mdx[c] * windat->dt2d;      
      /* Calculate the flux at e2 wet and boundary cells             */
      vfil = windat->u2av[c] + 0.5 * aconst *
        (windat->u2avb[c] - 2.0 * windat->u2av[c] + windat->nu2av[c]);
      windat->u2flux_z[c] = vfil * windat->depth_e2[c] * 
	window->h1au2[c] * wincon->mdy[c] * windat->dt2d;      
    }
  }
}

/* END set_zoom_flux()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate eta over the zoomed grid using the           */
/* interpolated 2D fluxes.                                           */
/*-------------------------------------------------------------------*/
void global_eta(geometry_t *geom,  /* Global geometry */
		double dt2d        /* 2D time-step */
		)
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xp1, yp1;
  double u1p, u1, u2p, u2, colflux;

  for (cc = 1; cc <= geom->nzinS; cc++) {
    c = geom->zin[cc].c;
    xp1 = geom->xp1[c];
    yp1 = geom->yp1[c];

    u1 = master->u1av[c] * master->depth_e1[c] * geom->h2au1[c] * dt2d;
    u1p = master->u1av[xp1] * master->depth_e1[xp1] * geom->h2au1[xp1] *
      dt2d;
    u2 = master->u2av[c] * master->depth_e2[c] * geom->h1au2[c] * dt2d;
    u2p = master->u2av[yp1] * master->depth_e2[yp1] * geom->h1au2[yp1] * 
      dt2d;
    colflux = u1p - u1 + u2p - u2;
    master->eta[c] = max(master->etab[c] - colflux / geom->cellarea[c],
			 geom->botz[c]);
  }
}

/* END global_eta()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the horizontal grid spacings for zoomed zones   */
/*-------------------------------------------------------------------*/
void hgrid_zoom(geometry_t *geom,       /* Global geometry           */
		geometry_t *window      /* Window geometry           */
  )
{
  int c, cg, cc;
  int wn, zc;
  int zs1 = (int)(window->zmfe1 / 2);
  int zs2 = (int)(window->zmfe2 / 2);
  int x1, x2, y1, y2;

  if (window->zoomf == 1)
    return;

  for (cc = 1; cc <= window->enonS; cc++) {
    c = window->w2_t[cc];
    c = cc;
    cg = x1 = x2 = y1 = y2 = window->wsa[c];
    wn = geom->fm[c].wn;

    for (zc = 1; zc <= zs1; zc++) {
      x1 = geom->xp1[x1];
      x2 = geom->xm1[x2];
      window->h1acell[c] += (geom->h1acell[x1] + geom->h1acell[x2]);
      window->h1au1[c] += (geom->h1au1[x1] + geom->h1au1[x2]);
      window->h1au2[c] += (geom->h1au2[x1] + geom->h1au2[x2]);

    }
    for (zc = 1; zc <= zs2; zc++) {
      y1 = geom->yp1[y1];
      y2 = geom->ym1[y2];
      window->h2acell[c] += (geom->h2acell[y1] + geom->h2acell[y2]);
      window->h2au2[c] += (geom->h2au2[y1] + geom->h2au2[y2]);
      window->h2au1[c] += (geom->h2au1[y1] + geom->h2au1[y2]);
    }
    window->cellarea[c] = window->h1acell[c] * window->h2acell[c];
  }

  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    x1 = window->wsa[window->xm1[c]];
    window->botzu1[c] = max(geom->botz[x1], window->botz[c]);
  }
  for (cc = 1; cc <= window->nbpte1S; cc++) {
    c = window->bpte1S[cc];
    window->botzu1[c] = HUGE_VAL;
  }
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    y1 = window->wsa[window->ym1[c]];
    window->botzu2[c] = max(geom->botz[y1], window->botz[c]);
  }
  for (cc = 1; cc <= window->nbpte2S; cc++) {
    c = window->bpte2S[cc];
    window->botzu2[c] = HUGE_VAL;
  }
}

/* END hgrid_zoom()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the horizontal grid spacings for zoomed zones   */
/*-------------------------------------------------------------------*/
void hgrid_zoom_m(geometry_t *geom,       /* Global geometry         */
		  geometry_t **window,    /* Window geometry         */
		  dump_data_t *dumpdata   /* Dump data               */
		  )
{
  int c, cc;
  int n, zc;
  int zs1, zs2; 
  int cp, cm;
  double *d1, *d2, *d3;
  double h1, h2;
  int *mask;

  if (geom->zoomf == 1)
    return;

 /* Make the mask for zoomed cells */
  mask = i_alloc_1d(geom->sgsizS);
  d1 = d_alloc_1d(geom->sgsizS);
  d2 = d_alloc_1d(geom->sgsizS);
  d3 = d_alloc_1d(geom->sgsizS);
  memset(mask, 0, geom->sgsizS * sizeof(int));

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    n = geom->fm[c].wn;
    zs1 = window[n]->zmfe1;
    zs2 = window[n]->zmfe2;
    if (zs1 > 1 || zs2 > 1) mask[c] = 1;
  }

  /* Reset the grid metrics in the e1 direction for zoomed cells */
  memcpy(d1, geom->h1acell, geom->sgsizS * sizeof(double));
  memcpy(d2, geom->h1au1, geom->sgsizS * sizeof(double));
  memcpy(d3, geom->h1au2, geom->sgsizS * sizeof(double));
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    n = geom->fm[c].wn;
    zs1 = (int)(window[n]->zmfe1 / 2);

    h1 = h2 = 0.0;
    cp = cm = c;
    for (zc = 1; zc <= zs1; zc++) {
      cp = geom->xp1[cp];
      cm = geom->xm1[cm];
      if (geom->zoomc[c] & (ZC|ZCB)) {
	geom->h1acell[c] += (d1[cp] + d1[cm]);
	h1 += d2[cp];
      }
      if (geom->zoomc[c] & (ZE1|ZE1B)) {
	geom->h1au1[c] += (d2[cp] + d2[cm]);
	h2 += d2[cm];
      }
      if (geom->zoomc[c] & (ZE2|ZE2B))
	geom->h1au2[c] += (d3[cp] + d3[cm]);
    }
    if (geom->zoomc[c] & (ZC|ZCB)) {
      zc = geom->xp1[cp];
      if (mask[cp] && !mask[zc])
	geom->h1au1[zc] += h1;
    }
    if (geom->zoomc[c] & (ZE1|ZE1B)) {
      zc = geom->xm1[c];
      if (mask[c] && !mask[zc])
	geom->h1au1[c] -= h2;
    }
  }

  /* Reset the grid metrics in the e2 direction for zoomed cells */
  memcpy(d1, geom->h2acell, geom->sgsizS * sizeof(double));
  memcpy(d2, geom->h2au2, geom->sgsizS * sizeof(double));
  memcpy(d3, geom->h2au1, geom->sgsizS * sizeof(double));
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    n = geom->fm[c].wn;
    zs2 = (int)(window[n]->zmfe2 / 2);

    h1 = h2 = 0.0;
    cp = cm = c;
    for (zc = 1; zc <= zs2; zc++) {
      cp = geom->yp1[cp];
      cm = geom->ym1[cm];
      if (geom->zoomc[c] & (ZC|ZCB)) {
	geom->h2acell[c] += (d1[cp] + d1[cm]);
	h1 += d2[cp];
      }
      if (geom->zoomc[c] & (ZE2|ZE2B)) {
	geom->h2au2[c] += (d2[cp] + d2[cm]);
	h2 += d2[cm];
      }
      if (geom->zoomc[c] & (ZE1|ZE1B))
	geom->h2au1[c] += (d3[cp] + d3[cm]);
    }
    if (geom->zoomc[c] & (ZC|ZCB)) {
      zc = geom->yp1[cp];
      if (mask[cp] && !mask[zc])
	geom->h2au2[zc] += h1;
    }
    if (geom->zoomc[c] & (ZE2|ZE2B)) {
      zc = geom->ym1[c];
      if (mask[c] && !mask[zc])
	geom->h2au2[c] -= h2;
    }
    geom->cellarea[c] = geom->h1acell[c] * d1[c];
  }

  /* Populate the new metrics throughout the whole zoomed cell */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = cm = cp = geom->w2_t[cc];
    n = geom->fm[c].wn;
    zs1 = window[n]->zmfe1;
    zs2 = window[n]->zmfe2;

    if (geom->zoomc[c] & (ZC|ZCB)) {
      fill_zoom_cell(geom, geom->h1acell, c, zs1, zs2);
      fill_zoom_cell(geom, geom->h2acell, c, zs1, zs2);
      fill_zoom_cell(geom, geom->cellarea, c, zs1, zs2);
    }
    if (geom->zoomc[c] & (ZE1|ZE1B)) {
      fill_zoom_cell(geom, geom->h1au1, c, zs1, zs2);
      fill_zoom_cell(geom, geom->h2au1, c, zs1, zs2);
    }
    if (geom->zoomc[c] & (ZE2|ZE2B)) {
      fill_zoom_cell(geom, geom->h1au2, c, zs1, zs2);
      fill_zoom_cell(geom, geom->h2au2, c, zs1, zs2);
    }
  }

  /* Reset the dumpdata structure */
  s2c_2d(geom, geom->h1au1, dumpdata->h1au1, geom->nfe1, geom->nce2);
  s2c_2d(geom, geom->h1au2, dumpdata->h1au2, geom->nce1, geom->nfe2);
  s2c_2d(geom, geom->h1acell, dumpdata->h1acell, geom->nce1, geom->nce2);
  s2c_2d(geom, geom->h2au1, dumpdata->h2au1, geom->nfe1, geom->nce2);
  s2c_2d(geom, geom->h2au2, dumpdata->h2au2, geom->nce1, geom->nfe2);
  s2c_2d(geom, geom->h2acell, dumpdata->h2acell, geom->nce1, geom->nce2);

  i_free_1d(mask);
  d_free_1d(d1);
  d_free_1d(d2);
  d_free_1d(d3);
}



void hgrid_zoom_p(geometry_t *geom,       /* Global geometry         */
		  geometry_t **window     /* Window geometry         */
		  )
{
  int c, lc, cc;
  int n, zc, nn, ll;
  int zs1, zs2; 
  int cp, cm;
  double h1, h2;

  if (geom->zoomf == 1)
    return;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    n = geom->fm[c].wn;
    zs1 = (int)(window[n]->zmfe1 / 2);
    zs2 = (int)(window[n]->zmfe2 / 2);
    if (geom->zoomc[c] & ZC) {
      lc = geom->fm[c].sc;
      if (zs1) window[n]->h1acell[lc] = geom->h1acell[c];
      if (zs2) window[n]->h2acell[lc] = geom->h2acell[c];
      window[n]->cellarea[lc] = geom->cellarea[c];
    }
    if (zs1 && geom->zoomc[c] & ZE1) {
      lc = c;
      for (zc = 1; zc <= zs1; zc++) lc = geom->xp1[lc];
      lc = geom->fm[lc].sc;
      window[n]->h1au1[lc] = geom->h1au1[c];
      window[n]->h2au1[lc] = geom->h2au1[c];
      window[n]->botzu1[lc] = geom->botzu1[c];
    }
    if (zs2 && geom->zoomc[c] & ZE2) {
      lc = c;
      for (zc = 1; zc <= zs2; zc++) lc = geom->yp1[lc];
      lc = geom->fm[lc].sc;
      window[n]->h1au2[lc] = geom->h1au2[c];
      window[n]->h2au2[lc] = geom->h2au2[c];
      window[n]->botzu2[lc] = geom->botzu2[c];
    } 

  }
  /* Reset auxiliary metrics at fine/coarse boundaries */
  for (nn = 1; nn <= geom->nwindows; nn++) {
    for (cc = 1; cc <= window[nn]->b2_t; cc++) {
      lc = window[nn]->w2_t[cc];
      /* Front boundaries */
      ll = window[nn]->yp1[lc];
      c = window[nn]->wsa[ll];
      n = geom->fm[c].wn;
      h1 = (double)window[n]->zmfe1;
      if (n != nn && h1 > 1.0) {
	window[nn]->h1acell[ll] /= h1;
	window[nn]->h1au1[ll] /= h1;
	window[nn]->h1au2[ll] /= h1;
	window[nn]->cellarea[ll] /= h1;
      }
      /* Back boundaries */
      ll = window[nn]->ym1[lc];
      c = window[nn]->wsa[ll];
      n = geom->fm[c].wn;
      h1 = (double)window[n]->zmfe1;
      if (n != nn && h1 > 1.0) {
	window[nn]->h1acell[ll] /= h1;
	window[nn]->h1au1[ll] /= h1;
	window[nn]->h1au2[ll] /= h1;
	window[nn]->cellarea[ll] /= h1;
      }
      /* Right boundaries */
      ll = window[nn]->xp1[lc];
      c = window[nn]->wsa[ll];
      n = geom->fm[c].wn;
      h2 = (double)window[n]->zmfe2;
      if (n != nn && h2 > 1.0) {
	window[nn]->h2acell[ll] /= h2;
	window[nn]->h2au1[ll] /= h2;
	window[nn]->h2au2[ll] /= h2;
	window[nn]->cellarea[ll] /= h2;
      }
      /* Left boundaries */
      ll = window[nn]->xm1[lc];
      c = window[nn]->wsa[ll];
      n = geom->fm[c].wn;
      h2 = (double)window[n]->zmfe2;
      if (n != nn && h2 > 1.0) {
	window[nn]->h2acell[ll] /= h2;
	window[nn]->h2au1[ll] /= h2;
	window[nn]->h2au2[ll] /= h2;
	window[nn]->cellarea[ll] /= h2;
      }
    }
  }
}

/* END hgrid_zoom_p()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to adjust the horizontal grid spacings for zoomed zones   */
/*-------------------------------------------------------------------*/
void consts_zoom(master_t *master,       /* Master data              */
		 geometry_t *geom,       /* Global geometry           */
		 geometry_t *window,     /* Window geometry          */
		 win_priv_t *wincon      /* Window constants         */
		 )
{
  int c, cc;
  int xp1, xm1, yp1, ym1;
  int zs1 = (int)(window->zmfe1 / 2);
  int zs2 = (int)(window->zmfe2 / 2);
  int x1, x2, y1, y2;

  if (window->zoomf == 1)
    return;

  /*-----------------------------------------------------------------*/
  /* Reset the precalculated constants.                              */
  /* u1 momentum equations.                                          */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->w2_e1[cc];
    xp1 = window->xp1[c];
    xm1 = window->xm1[c];
    wincon->u1c1[c] =
      -1.0 / (window->h1au1[c] * window->h1au1[c] * window->h2au1[c]);
    wincon->u1c3[c] =
      1.0 * (window->h1acell[c] -
             window->h1acell[xm1]) / (window->h1au1[c] * window->h1au1[c]);
    wincon->u1c4[c] =
      1.0 * (window->h2acell[c] -
             window->h2acell[xm1]) / (window->h1au1[c] * window->h2au1[c]);
    wincon->u1c6[c] = -1.0 / window->h1au1[c];
    x1 = window->wsa[c];
    if (zs1)
      x1 = geom->zse1[x1];
    x2 = geom->xm1[x1];
    wincon->u1c5[c] =
      1.0 * (master->coriolis[x1] + master->coriolis[x2]) / 2;
  }

  /* u2 momentum equations.                                          */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->w2_e2[cc];
    yp1 = window->yp1[c];
    ym1 = window->ym1[c];
    wincon->u2c1[c] =
      -1.0 / (window->h1au2[c] * window->h2au2[c] * window->h2au2[c]);
    wincon->u2c3[c] =
      1.0 * (window->h1acell[c] -
             window->h1acell[ym1]) / (window->h1au2[c] * window->h2au2[c]);
    wincon->u2c4[c] =
      1.0 * (window->h2acell[c] -
             window->h2acell[ym1]) / (window->h2au2[c] * window->h2au2[c]);
    wincon->u2c6[c] = -1.0 / window->h2au2[c];
    y1 = window->wsa[c];
    if (zs2)
      y1 = geom->zse2[y1];
    y2 = geom->ym1[y1];
    wincon->u2c5[c] =
      -1.0 * (master->coriolis[y1] + master->coriolis[y2]) / 2;
  }
}

/* END consts_zoom()                                                 */
/*-------------------------------------------------------------------*/


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

  /* Transfer the flux in the master to the slave. If the slave is a */
  /* zoomed window then sum the fluxes.  */
  if (window->zoomf == 1) {
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Al[lc] = Ag[c];
    }
  } else {
    int zc, zs = window->zmee2;
    int zm1;
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = c1 = c2 = evec[cc];
      Al[lc] = Ag[c];
      for (zc = 1; zc <= zs; zc++) {
        c1 = geom->ym1[c1];
        c2 = geom->yp1[c2];
        Al[lc] += (Ag[c1] + Ag[c2]);
	/* Code for non-uniform zoomed bathymetry
	zm1 = window->zm1[lc];
	if (zm1 == window->zm1[zm1]) {
	  c1 = geom->zm1[c1];
	  while(c1 != geom->zm1[c1]) {
	    Al[lc] += Ag[c1];
	    c1 = geom->zm1[c1];
	  }
	  c2 = geom->zm1[c2];
	  while(c2 != geom->zm1[c2]) {
	    Al[lc] += Ag[c2];
	    c2 = geom->zm1[c2];
	  }
	}
	*/
      }
    }
  }
}

/* END zflux_e1()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to sum e1 fluxes over a grid so as to be applicable to a  */
/* zoomed cell.                                                      */
/*-------------------------------------------------------------------*/
void zflux_e2(geometry_t *geom, /* Sparse global geometery */
              geometry_t *window, /* Processing window */
              double *Ag,       /* Global flux array */
              double *Al,       /* Local flux array */
              int *vec,         /* Cells on which A is defined */
              int *evec,        /* Cells on which Ag is defined */
              int nvec          /* Size of vec */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */
  int c1, c2;

  /* Transfer the flux in the master to the slave. If the slave is a */
  /* zoomed window then sum the fluxes.  */
  if (window->zoomf == 1) {
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Al[lc] = Ag[c];
    }
  } else {
    int zc, zs = window->zmee1;
    int zm1;
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = c1 = c2 = evec[cc];
      Al[lc] = Ag[c];
      for (zc = 1; zc <= zs; zc++) {
        c1 = geom->xm1[c1];
        c2 = geom->xp1[c2];
        Al[lc] += (Ag[c1] + Ag[c2]);
	/*if(Ag==master->u2flux3d)printf("%d %d %d %f\n",c,c1,c2,Al[lc]);*/
	/* Code for non-uniform zoomed bathymetry
	zm1 = window->zm1[lc];
	if (zm1 == window->zm1[zm1]) {
	  c1 = geom->zm1[c1];
	  while(c1 != geom->zm1[c1]) {
	    Al[lc] += Ag[c1];
	    c1 = geom->zm1[c1];
	  }
	  c2 = geom->zm1[c2];
	  while(c2 != geom->zm1[c2]) {
	    Al[lc] += Ag[c2];
	    c2 = geom->zm1[c2];
	  }
	}
	*/
      }
    }
  }
}

/* END zflux_e2()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to sum e1 velocity over a grid so as to be applicable to  */
/* a zoomed cell.                                                    */
/*-------------------------------------------------------------------*/
void zvel_filter(geometry_t *geom, /* Sparse global geometery */
		 geometry_t *window, /* Processing window */
		 double *Ag,       /* Global velocity array */
		 double *Al,       /* Local velocity array */
		 double zs,        /* zmee2 for e1, zmee1 for e2 */
		 int *vec,         /* Cells on which Al is defined */
		 int *evec,        /* Cells on which Ag is defined */
		 int nvec          /* Size of vec */		 

  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */
  int c1, c2;

  /* Transfer the flux in the master to the slave. If the slave is a */
  /* zoomed window then sum the fluxes.  */
  if (window->zoomf == 1) {
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Al[lc] = Ag[c];
    }
  } else {
    int zc;
    int zm1;
    double area, wgt;
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = c1 = c2 = evec[cc];
      /* Area weighted */
      Al[lc] = Ag[c];
      area = geom->cellarea[geom->m2d[c]];
      /* Shapiro zs == 1
      Al[lc] = 2.0 * Ag[c]; 
      area = 4.0;
      */
      /* Full weighting zs == 2
      Al[lc] = 3.0 * Ag[c]; 
      area = 9.0;
      */
      for (zc = 1; zc <= zs; zc++) {
        c1 = geom->ym1[c1];
	/* Area weighted */
        c2 = geom->yp1[c2];
        Al[lc] += (geom->cellarea[geom->m2d[c1]]*Ag[c1] + 
		   geom->cellarea[geom->m2d[c2]]*Ag[c2]);
	area += (geom->cellarea[geom->m2d[c1]] + geom->cellarea[geom->m2d[c2]]);
	/* Shapiro
        Al[lc] += (Ag[c1] + Ag[c2]);
	*/
	/* Full weighting
	wgt = (zc == 1) ? 2.0 : 1.0;
        Al[lc] += (wgt * Ag[c1] + wgt * Ag[c2]);
	*/
      }
      Al[lc] /= area;
    }
  }
}

/* END zvel_e2()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to sum e1 velocity over a grid so as to be applicable to  */
/* a zoomed cell.                                                    */
/*-------------------------------------------------------------------*/
void zvel_e2(geometry_t *geom, /* Sparse global geometery */
	     geometry_t *window, /* Processing window */
	     double *Ag,       /* Global velocity array */
	     double *Al,       /* Local velocity array */
	     int *vec,         /* Cells on which Al is defined */
	     int *evec,        /* Cells on which Ag is defined */
	     int nvec          /* Size of vec */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */
  int c1, c2;

  /* Transfer the flux in the master to the slave. If the slave is a */
  /* zoomed window then sum the fluxes.  */
  if (window->zoomf == 1) {
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Al[lc] = Ag[c];
    }
  } else {
    int zc, zs = window->zmee2;
    int zm1;
    double area, wgt;
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = c1 = c2 = evec[cc];
      /* Area weighted */
      Al[lc] = Ag[c];
      area = geom->cellarea[geom->m2d[c]];
      /* Shapiro zs == 1
      Al[lc] = 2.0 * Ag[c]; 
      area = 4.0;
      */
      /* Full weighting zs == 2
      Al[lc] = 3.0 * Ag[c]; 
      area = 9.0;
      */
      for (zc = 1; zc <= zs; zc++) {
        c1 = geom->ym1[c1];
	/* Area weighted */
        c2 = geom->yp1[c2];
        Al[lc] += (geom->cellarea[geom->m2d[c1]]*Ag[c1] + 
		   geom->cellarea[geom->m2d[c2]]*Ag[c2]);
	area += (geom->cellarea[geom->m2d[c1]] + geom->cellarea[geom->m2d[c2]]);
	/* Shapiro
        Al[lc] += (Ag[c1] + Ag[c2]);
	*/
	/* Full weighting
	wgt = (zc == 1) ? 2.0 : 1.0;
        Al[lc] += (wgt * Ag[c1] + wgt * Ag[c2]);
	*/
      }
      Al[lc] /= area;
    }
  }
}

/* END zvel_e2()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for interpolation of zoomed grid        */
/* variables on the master.                                          */
/*-------------------------------------------------------------------*/
void geom_interp(geometry_t *geom,  /* Sparse global geometery */
                 geometry_t **window  /* Processing window */
  )
{
  int c, c1, c2 = 0, cc;        /* Sparse coordinate / counter */
  int n, m, n1, n2;             /* Window counter */
  int zc;                       /* Interpolation structure counter */
  int *mask;                    /* Duumy mask */
  int i, j;
  int cp, c3;
  int lc = 0;                   /* Local wet and auxiliary cells */
  int zs, zmf, zf, zmn, zmt;    /* Zoom counter, zoom factor */
  int ze1 = 0, ze2 = 0;
  int cnr = 0, ce = 0, cn = 0, cne = 0;
  double x = 0.0, x0 = 0.0, xinc = 0.0;
  double y = 0.0, y0 = 0.0, yinc = 0.0;
  int wn;
  int *mx = NULL;
  int *my = NULL;

  if (!geom->zoomf)
    return;

  /*-----------------------------------------------------------------*/
  /* For PRECONDITIONED grids the zoom cell centre is saved and the  */
  /* interpolation over the zoom cell sets all values to the centre  */
  /* value.                                                          */
  if (geom->zoom & PRECOND) {
    geom->nzin = geom->nzinS = 0;
    /* Count the zoom cell centres */
    for (n = 1; n <= geom->nwindows; n++) {
      if (window[n]->zoomf == 1)
	continue;
      for (cc = 1; cc <= window[n]->b3_t; cc++) { 
	c = window[n]->w3_t[cc];
	c2 = window[n]->m2d[c];
	c3 = window[n]->wsa[c];
	if (window[n]->zoomc[c2] & ZC) {
	  geom->nzin += (window[n]->zmfe1 * window[n]->zmfe2 - 1);
	  if (c3 <= geom->enonS)
	    geom->nzinS += (window[n]->zmfe1 * window[n]->zmfe2 - 1);
	}
      }
    }
    geom->zin =
      (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzin + 1));
    geom->nzin = geom->nzinS;
    geom->nzinS = 0;
    /* Fill the interpolation array */
    for (n = 1; n <= geom->nwindows; n++) {
      if (window[n]->zoomf == 1)
	continue;
      for (cc = 1; cc <= window[n]->b3_t; cc++) {
	c = c1 = window[n]->w3_t[cc];
	c2 = window[n]->m2d[c2];
	c3 = cnr = window[n]->wsa[c];
	if (window[n]->zoomc[c2] & ZC) {
	  for (n2 = 1; n2 <= window[n]->zmee2; n2++)
	    c3 = geom->ym1[c3];
	  for (n1 = 1; n1 <= window[n]->zmee1; n1++)
	    c3 = geom->xm1[c3];
	  /* Fill the array A                                                */
	  for (n2 = 1; n2 <= window[n]->zmfe2; n2++) {
	    c1 = c3;
	    for (n1 = 1; n1 <= window[n]->zmfe1; n1++) {
	      if (c1 != cnr) {
		if (c1 <= geom->enonS) {
		  geom->nzinS++;
		  i = geom->nzinS;
		} else {
		  geom->nzin++;
		  i = geom->nzin;
		}
		geom->zin[i].c = c1;
		geom->zin[i].cnr = cnr;
		geom->zin[i].ce = window[n]->zmfe1;
		geom->zin[i].cn = window[n]->zmfe2;
	      }
	      c1 = geom->xp1[c1];
	    }
	    c3 = geom->yp1[c3];
	  }
	}
      }
    }
    return;
  }

  mask = i_alloc_1d(geom->enon + 1);
  /*-----------------------------------------------------------------*/
  /* Get the local interpolation information */
  /* Count the number of non-center zoom wet cells */
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nzin = window[n]->nzinS = 0;
    if (window[n]->zoomf == 1)
      continue;
    for (cc = 1; cc <= window[n]->enon; cc++) { /* local sparse coord.  */
      c = window[n]->wsa[cc];   /* global sparse coord.  */
      if (window[n]->zoomc[cc] & ZC) {
        c2 = c;                 /* global zoom center */
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;              /* global interzoom cell */
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;  /* local interzoom cell */
            if (!(window[n]->zoomc[lc] & ZC) &&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {
              window[n]->nzin++;
              if (cc <= window[n]->enonS)
                window[n]->nzinS++;
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }
    window[n]->zin = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                             (window[n]->nzin + 1));

    /* Get the interpolation corner locations */
    zc = 1;
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      if (window[n]->zoomc[cc] & ZC) {
        c2 = c;
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;

            if (!(window[n]->zoomc[lc] & ZC) &&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {

              window[n]->zin[zc].c = lc;

              /* If cc is an ghost cell find the interior wet cell */
              cnr = cp = cc;
              if (!geom->fm[c].wn) {
                c3 = ANY(cc, window[n]->bpt, window[n]->nbpt);
                cp = window[n]->bin[c3];
              }
              window[n]->zin[zc].cnr = cp;

              /* Corner cells are self pointing over ghost cells */
              cp = window[n]->xp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpt, window[n]->nbpt);
                cp = window[n]->bin[c3];
              }
              window[n]->zin[zc].ce = cp;

              cp = window[n]->yp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpt, window[n]->nbpt);
                cp = window[n]->bin[c3];
              }
              window[n]->zin[zc].cn = cp;

              cp = window[n]->yp1[window[n]->xp1[cnr]];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpt, window[n]->nbpt);
                cp = window[n]->bin[c3];
              }
              window[n]->zin[zc].cne = cp;

              window[n]->zin[zc].x =
                (double)(i - 1) / (double)window[n]->zmfe1;
              window[n]->zin[zc].y =
                (double)(j - 1) / (double)window[n]->zmfe2;
              zc++;
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }
    if (zc - 1 != window[n]->nzin)
      printf
        ("WARNING : Incorrect number of zoom interpolation cells located, window %d : %d %d\n",
         window[n]->wn, zc - 1, window[n]->nzin);
  }

  /*-----------------------------------------------------------------*/
  /* Get the global interpolation information */
  /* Count the number of interpolation cells */
  geom->nzin = geom->nzinS = 0;
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  /* Set the mask to zero on global OBC's and count global OBC       */
  /* interpolation cells in the zoom window.                         */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      wn = geom->fm[c].wn;
      zmf = (open->type & U1BDRY) ? window[wn]->zmfe1 : window[wn]->zmfe2;
      if (zmf > 1) {
	for (zs = 1; zs <= (int)zmf / 2; zs++) {
	  geom->nzin++;
	  if (c <= geom->enonS)
	    geom->nzinS++;
	  mask[c] = 1;
	  c = open->nmap[c];
	}
	if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) {
	  if (!(window[wn]->zoomc[window[wn]->fm[c].sc] & ZC)) {
	    geom->nzin++;
	    if (c <= geom->enonS)
	      geom->nzinS++;
	  }
	  mask[c] = 1;
	}
      }
    }
  }

  /* Count remaining interpolation cells                             */
  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzin; zc++) {
      c = window[n]->wsa[window[n]->zin[zc].c];
      if (!mask[c]) {
        geom->nzin++;
        if (c <= geom->enonS)
          geom->nzinS++;
        mask[c] = 1;

      }
    }
  }

  geom->zin =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzin + 1));
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  geom->nzin = geom->nzinS;
  geom->nzinS = 0;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      wn = geom->fm[c].wn;
      zmn = (open->type & U1BDRY) ? window[wn]->zmfe1 : window[wn]->zmfe2;
      zmt = (open->type & U1BDRY) ? window[wn]->zmfe2 : window[wn]->zmfe1;
      zf = (int)zmn / 2;
      /* Only process cells in zoomed windows that have not been         */
      /* processed (i.e. mask[c] == 0).                                  */
      if (!mask[c] && zmn > 1) {
	/* Find the zoomed grid cell center interior to the boundary     */
	for (zs = 1; zs <= zf; zs++)
	  c = open->nmap[c];
	m = 0;
	for (zc = 1; zc <= zmt; zc++) {
	  c = open->tmpp[c];
	  lc = window[wn]->fm[c].sc;
	  m++;
	  if (window[wn]->zoomc[lc] & ZC) break;
	}
	/* Find the corner points and map pointers                       */
	if (open->ocodex & L_EDGE) {	    
	  cn = cne = c;
	  if (geom->fm[cn].wn > 0 && geom->fm[cne].wn != wn)
	    cn = window[wn]->wsa[window[wn]->xm1[lc]];
	  ce = cnr = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  if (geom->fm[cn].wn > 0 && geom->fm[ce].wn != wn)
	    cnr = window[wn]->wsa[window[wn]->xm1[lc]];
	  c2 = open->obc_t[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 1.0 / (double)(zf + 1.0); 
	  xinc = 1.0 / (double)(zf + 1.0);
	  y0 = (double)(zmt - m) / (double)zmt;
	  yinc = 1.0 / (double)zmt;
	  ze1 = zf;
	  ze2 = zmt;
	  ze2 = m;
	}
	if (open->ocodex & R_EDGE) {
	  cn = cne = c;
	  if (geom->fm[cn].wn > 0 && geom->fm[cn].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      cne = geom->xp1[cne];
	  }
	  cnr = ce = window[wn]->wsa[(lc = window[wn]->ym1[lc])];

	  if (geom->fm[cn].wn > 0 && geom->fm[cnr].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      ce = geom->xp1[ce];
	  }
	  c2 = open->obc_t[cc];
	  for (zc = 1; zc <= zf; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0;
	  xinc = 1.0 / (double)(zf + 1.0);
	  y0 = (double)(zmt - m) / (double)zmt;
	  yinc = 1.0 / (double)zmt;
	  ze1 = zf + 1;
	  ze2 = m;
	}
	if (open->ocodey & B_EDGE) {
	  cne = ce = c;
	  if (geom->fm[cn].wn > 0 && geom->fm[cne].wn != wn)
	    ce = window[wn]->wsa[window[wn]->ym1[lc]];
	  cn = cnr = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  if (geom->fm[cn].wn > 0 && geom->fm[cne].wn != wn)
	    cnr = window[wn]->wsa[window[wn]->ym1[lc]];
	  c2 = open->obc_t[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = (double)(zmt - m) / (double)zmt;
	  xinc = 1.0 / (double)zmt;
	  y0 = 1.0 / (double)(zf + 1.0); 
	  yinc = 1.0 / (double)(zf + 1.0);
	  ze2 = zf;
	  ze1 = m;
	}
	if (open->ocodey & F_EDGE) {
	  cne = ce = c;
	  if (geom->fm[cn].wn > 0 && geom->fm[ce].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      cne = geom->yp1[cne];
	  }
	  cnr = cn = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  if (geom->fm[cn].wn > 0 && geom->fm[cnr].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      cn = geom->yp1[cn];
	  }
	  c2 = open->obc_t[cc];
	  for (zc = 1; zc <= zf; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->ym1;
	  x0 = (double)(zmt - m) / (double)zmt;
	  xinc = 1.0 / (double)zmt;
	  y0 = 0.0;
	  yinc = 1.0 / (double)(zf + 1.0);
	  ze2 = zf + 1;
	  ze1 = m;
	}

	/* Fill all global cells bounded by the corners and re-set the   */
	/* mask for these cells.                                         */
	x = x0;
	for (zs = 1; zs <= ze1; zs++) {
	  c = c2;
	  y = y0;
	  for (zc = 1; zc <= ze2; zc++) {
	    if (!mask[c] && geom->fm[c].wn == wn &&
		!(window[wn]->zoomc[window[wn]->fm[c].sc] & ZC)) {
	      if (c <= geom->enonS) {
		geom->nzinS++;
		i = geom->nzinS;
		/*
		print_zoom_bdry("cb" c, cnr, ce, cn, cne, x, y);
		*/
	      } else {
		geom->nzin++;
		i = geom->nzin;
	      }
	      geom->zin[i].cnr = cnr;
	      geom->zin[i].ce = ce;
	      geom->zin[i].cn = cn;
	      geom->zin[i].cne = cne;
	      geom->zin[i].c = c;
	      geom->zin[i].x = x;
	      geom->zin[i].y = y;
	      mask[c] = 1;
	    }
	    y += yinc;
	    c = my[c];
	  }
	  x += xinc;
	  c2 = mx[c2];
	}
      }
    }
  }

  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzin; zc++) {
      c = window[n]->wsa[window[n]->zin[zc].c];
      if (!mask[c]) {

        cnr = window[n]->wsa[window[n]->zin[zc].cnr];
        ce = window[n]->wsa[window[n]->zin[zc].ce];
        cn = window[n]->wsa[window[n]->zin[zc].cn];
        cne = window[n]->wsa[window[n]->zin[zc].cne];
        x = window[n]->zin[zc].x;
        y = window[n]->zin[zc].y;

	/* Save the corner global cells                              */
	if (c <= geom->enonS) {
	  geom->nzinS++;
	  i = geom->nzinS;
	} else {
	  geom->nzin++;
	  i = geom->nzin;
	}
	geom->zin[i].cnr = cnr;
	geom->zin[i].ce = ce;
	geom->zin[i].cn = cn;
	geom->zin[i].cne = cne;
	geom->zin[i].c = c;
	geom->zin[i].x = x;
	geom->zin[i].y = y;
        mask[c] = 1;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Deallocate unrequired structure memory */
  i_free_1d(mask);
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->zoomf > 1)
      free((zoom_interp_t *)window[n]->zin);
    window[n]->nzin = window[n]->nzinS = 0;
  }
}

/* END geom_interp()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for interpolation of zoomed variables   */
/* variables on e1 faces in the master. Over ghost cells the corner  */
/* interpolation calls are only self pointing for tangential faces,  */
/* normal e1 faces use the (zero valued) ghost cell.                 */
/*-------------------------------------------------------------------*/
void geom_interp_e1(geometry_t *geom, /* Sparse global geometery */
                    geometry_t **window /* Processing window */
  )
{
  int c, c1, c2 = 0, cc;        /* Sparse coordinate / counter */
  int n, m;                     /* Window counter */
  int zc;                       /* Interpolation structure counter */
  int *mask;                    /* Duumy mask */
  int *c2cc;
  int lc = 0;                   /* Local wet and auxiliary cells */
  int zs, zmf, zf, zmn, zmt;    /* Zoom counter, zoom factor */
  int ze1 = 0, ze2 = 0;
  int cnr = 0, ce = 0, cn = 0, cne = 0;
  double x = 0.0, x0 = 0.0, xinc = 0.0;
  double y = 0.0, y0 = 0.0, yinc = 0.0;
  int wn;
  int *mx = NULL;
  int *my = NULL;
  int i, j;
  int zm1, cp, c3;

  if (!geom->zoomf)
    return;

  /* Set cells to process for increasing horizontal diffusion */
  geom_bdry_he1(geom, window);

  /* Set cells to process for smoothing at coarse-fine boundaries */
  /*geom_bdry_e1(geom, window);*/

  /* Set cells to process for smoothing at fine boundaries */
  geom_m2s(geom, window);

  /* Set cell to process for dual linear interpolations */
  /*geom_sinterp_e1(geom, window);*/

  /* Allocate the c to cc array */
  c2cc = i_alloc_1d(geom->sgsiz);
  for (cc = 1; cc <= geom->b2_e1; cc++) {
    c = geom->w2_e1[cc];
    c2cc[c] = cc;
  }

  if (geom->zoom & PRECOND) {
    geom->nzine1 = geom->nzin;
    geom->nzine1S = geom->nzinS;
    geom->zine1 =
      (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzine1 + 1));
    for (cc = 1; cc <= geom->nzine1; cc++) {
      geom->zine1[cc].c = c = geom->zin[cc].c;
      geom->zine1[cc].cnr = geom->zin[cc].cnr;
      geom->zine1[cc].ce = geom->zin[cc].ce;
      geom->zine1[cc].cn = geom->zin[cc].cn;
      c2 = c2cc[geom->m2d[c]];
      if (c2cc[c2]) {
	geom->zine1[cc].cb = geom->bot_e1[c2];
      } 
      if (!c2cc[c2] || !geom->zine1[cc].cb) {
	while(c != geom->zm1[c])
	  c = geom->zm1[c];
	geom->zine1[cc].cb = geom->zp1[c];
      }
    }
    i_free_1d(c2cc);
    return;
  }

  mask = i_alloc_1d(geom->sgsiz);

  /* Set the mask adjacent to southern solid boundaries */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->v3_e1; cc++) {
    c = geom->w3_e1[cc];      /* Global cell to process */
    c1 = geom->xm1[c];        /* Global cell to the west of c */
    if (geom->fm[c1].wn && !(geom->fm[geom->xm1[c1]].wn))
      mask[c1] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the local interpolation information */
  /* Count the number of non-center zoom wet cells */
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nzine1 = window[n]->nzine1S = 0;
    if (window[n]->zoomf == 1)
      continue;
    for (cc = 1; cc <= window[n]->enon; cc++) { /* local sparse coord.  */
      c = window[n]->wsa[cc];   /* global sparse coord.  */
      if (window[n]->zoomc[cc] & ZE1) {
        c2 = c;                 /* global zoom center */
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;              /* global interzoom cell */
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;  /* local interzoom cell */
            if (!(window[n]->zoomc[lc] & ZE1) &&
		(ANY(c1, geom->w3_e1, geom->b3_e1) >= 0) &&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {
              window[n]->nzine1++;
              if (cc <= window[n]->enonS) {
                window[n]->nzine1S++;
	      }
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }
    window[n]->zine1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nzine1 + 1));

    /* Get the interpolation corner locations */
    zc = 1;
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];

      if (window[n]->zoomc[cc] & ZE1) {
        c2 = c;
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;

            if (!(window[n]->zoomc[lc] & ZE1) &&
		(ANY(c1, geom->w3_e1, geom->b3_e1) >= 0) &&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {

              window[n]->zine1[zc].c = lc;

              /* If cc is an ghost cell find the interior wet cell.  */
              /* Note : an interior cell to a ghost cell, bine1, is */
              /* only used for tangential cells (c3>window->nbe1) so */
              /* as to emulate a no-gradient condition across the */
              /* ghost cell.  */
              /* The interior maps, bine1, always map to zoom cell */
              /* centers in the window; the correct location for */
              /* tangential e1 velocities in the global (shifted) */
              /* grid. The interior map, bine1, cannot be used for */
              /* normal velocities in the global grid, since these */
              /* are shifted zoomf/2 from the cell center: however */
              /* the normal velocities must be set to zero at the */
              /* ghost cell (zero flux condition) hence the zero */
              /* valued ghost cell is used as the corner cell and it */
              /* is not required to find the interior cell via the */
              /* interior map bine1.  */
              cnr = cp = cc;
              if (!geom->fm[c].wn) {
                c3 = ANY(cc, window[n]->bpte1, window[n]->nbpte1);
                if (c3 > window[n]->nbe1)
                  cp = window[n]->bine1[c3];
              }
              window[n]->zine1[zc].cnr = cp;

              cp = window[n]->xp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte1, window[n]->nbpte1);
                if (c3 > window[n]->nbe1)
                  cp = window[n]->bine1[c3];
              }
              window[n]->zine1[zc].ce = cp;

              cp = window[n]->yp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte1, window[n]->nbpte1);
                if (c3 > window[n]->nbe1)
                  cp = window[n]->bine1[c3];
              }
              window[n]->zine1[zc].cn = cp;

              cp = window[n]->yp1[window[n]->xp1[cnr]];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte1, window[n]->nbpte1);
                if (c3 > window[n]->nbe1)
                  cp = window[n]->bine1[c3];
              }
              window[n]->zine1[zc].cne = cp;
              window[n]->zine1[zc].x =
                (double)(i - 1) / (double)window[n]->zmfe1;
              window[n]->zine1[zc].y =
                (double)(j - 1) / (double)window[n]->zmfe2;
              zc++;
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }
    if (zc - 1 != window[n]->nzine1)
      printf
        ("WARNING : Incorrect number of zoom interpolation e1 faces located, window %d : %d %d\n",
         window[n]->wn, zc - 1, window[n]->nzine1);

  }

  /*-----------------------------------------------------------------*/
  /* Get the global interpolation information                        */
  /*-----------------------------------------------------------------*/
  /* Count the number of interpolation cells                         */
  geom->nzine1 = geom->nzine1S = 0;
  /* Set the mask at all valid cells in the master corresponding to  */
  /* the window's zoomed cell to process.                            */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->b3_e1; cc++) {
    c2 = cp = geom->w3_e1[cc];
    wn = geom->fm[c2].wn;
    if (geom->zoomc[geom->m2d[c2]] & (ZE1 | ZE1B)) {
      ze1 = window[wn]->zmfe1;
      ze2 = window[wn]->zmfe2;
      for (zc = 1; zc <= (int)ze2 / 2; zc++)
	c2 = geom->ym1[c2];
      for (zc = 1; zc <= ze2; zc++) {
	c = c2;
	for (zs = 1; zs <= ze1; zs++) {
	  mask[c] = 1;
	  /* Non-uniform zoomed bathymetry 
	     Keep going to the bottom for the last cell
	  zm1 = geom->zm1[cp];
	  if (geom->zm1[zm1] == zm1) {
	    c3 = geom->zm1[c];
	    while(c3 != geom->zm1[c3]) {
	      mask[c3] = 1;
	      c3 = geom->zm1[c3];
	    }
	  }
	  */
	  c = geom->xp1[c];
	}
	c2 = geom->yp1[c2];
      }
    }
  }

  /* Count global OBC interpolation cells in the zoom window         */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->to3_e1; cc++) {
      c = open->obc_e1[cc];
      wn = geom->fm[open->nmap[c]].wn;
      if (wn == 0) 
	wn = geom->fm[c].wn;	
      zmf = (open->type & U1BDRY) ? window[wn]->zmfe1 : window[wn]->zmfe2;
      if (window[wn]->zoomf > 1) {
	if (cc <= open->no3_e1)
	  zf = (int)zmf;      /* Normal cells          */
	else
	  zf = (int)zmf / 2;  /* Tangential cells      */
	/* Normal cells on right edges; include an extra row         */
	if (open->ocodex & R_EDGE) {
	  if (mask[c] && !(geom->zoomc[geom->m2d[c]] & ZE1NB)) {
	    geom->nzine1++;
	    if (c <= geom->enonS) {
	      geom->nzine1S++;
	    }
	    mask[c] = 0;
	  }
	  c = open->nmap[c];
	}
	/* Normal and tangential cells                               */
	for (zs = 1; zs <= zf; zs++) {
	  if (mask[c] && !(geom->zoomc[geom->m2d[c]] & (ZE1|ZE1NB|ZE1TB))) {
	    geom->nzine1++;
	    if (c <= geom->ewetS) {
	      geom->nzine1S++;
	    }
	    mask[c] = 0;
	  }
	  c = open->nmap[c];
	}
	/* Tangential cells on front edges; include an extra row     */
	if (open->ocodey & F_EDGE) {
	  if (mask[c] && !(window[wn]->zoomc[window[wn]->fm[c].sc] & ZE1)) {
	    geom->nzine1++;
	    if (c <= geom->enonS)
	      geom->nzine1S++;
	    mask[c] = 0;
	  }
	}
      }
    }
  }

  /* Count remaining interpolation cells                             */
  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzine1; zc++) {
      c = window[n]->wsa[window[n]->zine1[zc].c];
      if (mask[c]) {
        geom->nzine1++;
        if (c <= geom->enonS)
          geom->nzine1S++;
        mask[c] = 0;
      }
    }
  }

  /* Fill the global interpolation structure */
  geom->zine1 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzine1 + 1));
  /* Set the mask at all valid cells in the master corresponding to  */
  /* the window's zoomed cell to process.                            */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->b3_e1; cc++) {
    c2 = geom->w3_e1[cc];
    wn = geom->fm[c2].wn;
    if (geom->zoomc[geom->m2d[c2]] & (ZE1 | ZE1B)) {
      ze1 = window[wn]->zmfe1;
      ze2 = window[wn]->zmfe2;
      for (zc = 1; zc <= (int)ze2 / 2; zc++)
	c2 = geom->ym1[c2];
      for (zc = 1; zc <= ze2; zc++) {
	c = c2;
	for (zs = 1; zs <= ze1; zs++) {
	  mask[c] = 1;
	  c = geom->xp1[c];
	  /* Non-uniform zoomed bathymetry 
	     Keep going to the bottom for the last cell
	  zm1 = geom->zm1[cp];
	  if (geom->zm1[zm1] == zm1) {
	    c3 = geom->zm1[c];
	    while(c3 != geom->zm1[c3]) {
	      mask[c3] = 1;
	      c3 = geom->zm1[c3];
	    }
	  }
	  */
	}
	c2 = geom->yp1[c2];
      }
    }
  }
  geom->nzine1 = geom->nzine1S;
  geom->nzine1S = 0;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->to3_e1; cc++) {
      c = open->obc_e1[cc];
      wn = geom->fm[open->nmap[c]].wn;
      if (wn == 0) 
	wn = geom->fm[c].wn;	
      zmn = window[wn]->zmfe1;
      zmt = window[wn]->zmfe2;
      if (cc <= open->no3_e1)
	zf = (int)zmn;      /* Normal cells          */
      else
	zf = (int)zmn / 2;  /* Tangential cells      */
      /* Only process cells in zoomed windows that have not been         */
      /* processed (i.e. mask[c] == 1).                                  */
      if (mask[c] && window[wn]->zoomf > 1) {

	/* Normal cells on R_EDGE                                        */
	m = 0;
	if (open->ocodex & R_EDGE) {
	  c2 = c;
	  for (zc = 1; zc <= zmt; zc++) {
	    c2 = open->tmpp[c2];
	    lc = window[wn]->fm[c2].sc;
	    m++;
	    if (geom->fm[c2].wn == wn && geom->zoomc[geom->m2d[c2]] & ZE1NB)
	      break;
	  }
	  cn = cne = c2;
	  cnr = ce = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  c2 = open->obc_e1[cc];
	  x  = 0.0; 
	  y = (double)(zmt - m) / (double)zmt;
	  yinc = 1.0 / (double)zmt;
	  for (zc = 1; zc <= m; zc++) {
	    if (!(geom->zoomc[geom->m2d[c2]] & ZE1NB) && 
		mask[c2] && geom->fm[c2].wn == wn) {
	      if (c2 <= geom->enonS) {
		geom->nzine1S++;
		i = geom->nzine1S;
	      } else {
		geom->nzine1++;
		i = geom->nzine1;
	      }
	      geom->zine1[i].cb = geom->bot_e1[c2cc[geom->m2d[c2]]];
	      geom->zine1[i].cnr = cnr;
	      geom->zine1[i].ce = ce;
	      geom->zine1[i].cn = cn;
	      geom->zine1[i].cne = cne;
	      geom->zine1[i].c = c2;
	      geom->zine1[i].x = x;
	      geom->zine1[i].y = y;
	      mask[c2] = 0;
	    }
	    y += yinc;
	    c2 = geom->yp1[c2];
	  }
	}

	/* Find the zoomed grid cell center interior to the boundary     */
	m = 0;
	for (zs = 1; zs <= zf; zs++)
	  c = open->nmap[c];
	for (zc = 1; zc <= zmt; zc++) {
	  c = open->tmpp[c];
	  lc = window[wn]->fm[c].sc;
	  m++;
	  if (geom->fm[c].wn == wn && window[wn]->zoomc[lc] & ZE1) break;
	}
	/* Find the corner points and map pointers                       */
	if (open->ocodex & L_EDGE) {	    
	  cne = c;
	  cn = window[wn]->wsa[window[wn]->xm1[lc]];
	  ce = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  cnr = window[wn]->wsa[window[wn]->xm1[lc]];
	  c2 = open->obc_e1[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0;
	  xinc = 1.0 / (double)zmn;
	  y0 = (double)(zmt - m) / (double)zmt;
	  yinc = 1.0 / (double)zmt;
	  ze1 = zmn;
	  ze2 = m;
	}
	if (open->ocodex & R_EDGE) {
	  cn = c;
	  cne = window[wn]->wsa[window[wn]->xp1[lc]];
	  if (geom->fm[cne].wn != wn && cn == cne) {
	    for (zc = 1; zc <= zmn; zc++)
	      cne = geom->xp1[cne];
	  }
	  cnr = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  ce = window[wn]->wsa[window[wn]->xp1[lc]];
	  c2 = open->obc_e1[cc];
	  for (zc = 1; zc <= zmn; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0;
	  xinc = 1.0 / (double)zmn;
	  y0 = (double)(zmt - m) / (double)zmt;
	  yinc = 1.0 / (double)zmt;
	  ze1 = zmn;
	  ze2 = m;
	}
	if (open->ocodey & B_EDGE) {
	  cne = ce = c;
	  if (geom->fm[cne].wn != wn)
	    ce = window[wn]->wsa[window[wn]->ym1[lc]];
	  cn = cnr = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  if (geom->fm[cn].wn != wn)
	    cnr = window[wn]->wsa[window[wn]->ym1[lc]];
	  c2 = open->obc_e1[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0; 
	  xinc = 1.0 / (double)zmn;
	  y0 = 1.0 / (double)(zf + 1.0); 
	  yinc = 1.0 / (double)(zf + 1.0);
	  ze2 = zf;
	  ze1 = m;
	}
	if (open->ocodey & F_EDGE) {
	  ce = cne = c;
	  if (geom->fm[ce].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      cne = geom->yp1[cne];
	  }
	  cnr = cn = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  if (geom->fm[cnr].wn != wn) {
	    for (zc = 1; zc <= (int)zf; zc++)
	      cn = geom->yp1[cn];
	  }
	  c2 = open->obc_e1[cc];
	  for (zc = 1; zc <= zf; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0; 
	  xinc = 1.0 / (double)zmn;
	  y0 = 0.0;
	  yinc = 1.0 / (double)(zf + 1.0);
	  ze2 = zf + 1;
	  ze1 = m;
	}

	/* Fill all global cells bounded by the corners and re-set the   */
	/* mask for these cells.                                         */
	x = x0;
	for (zs = 1; zs <= ze1; zs++) {
	  c = c2;
	  y = y0;
	  for (zc = 1; zc <= ze2; zc++) {
	    if (!(geom->zoomc[geom->m2d[c]] & (ZE1 | ZE1NB | ZE1TB)) && 
		mask[c] && geom->fm[c].wn == wn) {
	      if (c <= geom->enonS) {
		geom->nzine1S++;
		i = geom->nzine1S;
		/*
		print_zoom_bdry("u1b", c, cnr, ce, cn, cne, x, y);
		*/
	      } else {
		geom->nzine1++;
		i = geom->nzine1;
	      }
	      geom->zine1[i].cb = geom->bot_e1[c2cc[geom->m2d[c]]];
	      geom->zine1[i].cnr = cnr;
	      geom->zine1[i].ce = ce;
	      geom->zine1[i].cn = cn;
	      geom->zine1[i].cne = cne;
	      geom->zine1[i].c = c;
	      geom->zine1[i].x = x;
	      geom->zine1[i].y = y;
	      mask[c] = 0;
	    }
	    y += yinc;
	    c = my[c];
	  }
	  x += xinc;
	  c2 = mx[c2];
	}
      }
    }
  }

  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzine1; zc++) {
      c = window[n]->wsa[window[n]->zine1[zc].c];

      if (mask[c]) {
        cnr = window[n]->wsa[window[n]->zine1[zc].cnr];
        ce = window[n]->wsa[window[n]->zine1[zc].ce];
        cn = window[n]->wsa[window[n]->zine1[zc].cn];
        cne = window[n]->wsa[window[n]->zine1[zc].cne];
        x = window[n]->zine1[zc].x;
        y = window[n]->zine1[zc].y;

	/* Save the corner global cells                              */
	if (c <= geom->enonS) {
	  geom->nzine1S++;
	  i = geom->nzine1S;
	} else {
	  geom->nzine1++;
	  i = geom->nzine1;
	}
	geom->zine1[i].cb = geom->bot_e1[c2cc[geom->m2d[c]]];
	geom->zine1[i].cnr = cnr;
	geom->zine1[i].ce = ce;
	geom->zine1[i].cn = cn;
	geom->zine1[i].cne = cne;
	geom->zine1[i].c = c;
	geom->zine1[i].x = x;
	geom->zine1[i].y = y;
        mask[c] = 0;
	if (!geom->zine1[i].cb) {
	  c1 = c;
	  while(c1 != geom->zm1[c1])
	    c1 = geom->zm1[c1];
	  geom->zine1[i].cb = geom->zp1[c1];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Deallocate unrequired structure memory */
  i_free_1d(mask);
  i_free_1d(c2cc);
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->zoomf > 1)
      free((zoom_interp_t *)window[n]->zine1);
    window[n]->nzine1 = window[n]->nzine1S = 0;
  }
}

/* END geom_interp_e1()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for interpolation of zoomed variables   */
/* variables on e2 faces in the master. Over ghost cells the corner  */
/* interpolation calls are only self pointing for tangential faces,  */
/* normal e1 faces use the (zero valued) ghost cell.                 */
/*-------------------------------------------------------------------*/
void geom_interp_e2(geometry_t *geom, /* Sparse global geometery */
                    geometry_t **window /* Processing window */
  )
{
  int c, c1, c2 = 0, cc;        /* Sparse coordinate / counter */
  int n, m;                     /* Window counter */
  int zc;                       /* Interpolation structure counter */
  int *mask;                    /* Duumy mask */
  int *c2cc;
  int lc = 0;                   /* Local wet and auxiliary cells */
  int zs, zmf, zf, zmn, zmt;    /* Zoom counter, zoom factor */
  int ze1 = 0, ze2 = 0;
  int cnr = 0, ce = 0, cn = 0, cne = 0;
  double x = 0.0, x0 = 0.0, xinc = 0.0;
  double y = 0.0, y0 = 0.0, yinc = 0.0;
  int wn;
  int *mx = NULL;
  int *my = NULL;
  int i, j;
  int cp, c3;

  if (!geom->zoomf)
    return;

  /* Set cells to process for increasing horizontal diffusion */

  geom_bdry_he2(geom, window);

  /* Set cells to process for smoothing at coarse-fine boundaries */
  /*geom_bdry_e2(geom, window);*/

  /* Set cell to process for dual linear interpolations */
  /*geom_sinterp_e2(geom, window);*/

  /* Allocate the c to cc array */
  c2cc = i_alloc_1d(geom->sgsiz);
  for (cc = 1; cc <= geom->b2_e2; cc++) {
    c = geom->w2_e2[cc];
    c2cc[c] = cc;
  }

  if (geom->zoom & PRECOND) {
    geom->nzine2 = geom->nzin;
    geom->nzine2S = geom->nzinS;
    geom->zine2 =
      (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzine2 + 1));
    for (cc = 1; cc <= geom->nzine2; cc++) {
      c = geom->m2d[geom->zin[cc].c];
      geom->zine2[cc].c = geom->zin[cc].c;
      geom->zine2[cc].cnr = geom->zin[cc].cnr;
      geom->zine2[cc].ce = geom->zin[cc].ce;
      geom->zine2[cc].cn = geom->zin[cc].cn;
      c2 = c2cc[geom->m2d[c]];
      if (c2cc[c2]) {
	geom->zine2[cc].cb = geom->bot_e2[c2];
      }
      if (!c2cc[c2] || !geom->zine2[cc].cb) {
	while(c != geom->zm1[c])
	  c = geom->zm1[c];
	geom->zine2[cc].cb = geom->zp1[c];
      }
    }
    i_free_1d(c2cc);
    return;
  }

  mask = i_alloc_1d(geom->sgsiz);

  /* Set the mask adjacent to southern solid boundaries */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->v3_e2; cc++) {
    c = geom->w3_e2[cc];      /* Global cell to process */
    c1 = geom->ym1[c];        /* Global cell to the south of c */
    if (geom->fm[c1].wn && !(geom->fm[geom->ym1[c1]].wn))
      mask[c1] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the local interpolation information */
  /* Count the number of non-center zoom wet cells */
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nzine2 = window[n]->nzine2S = 0;
    if (window[n]->zoomf == 1)
      continue;
    for (cc = 1; cc <= window[n]->enon; cc++) { /* local sparse coord.  */
      c = window[n]->wsa[cc];   /* global sparse coord.  */

      if (window[n]->zoomc[cc] & ZE2) {
        c2 = c;                 /* global zoom center */
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;              /* global interzoom cell */
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;  /* local interzoom cell */
            if (!(window[n]->zoomc[lc] & ZE2) && 
		(ANY(c1, geom->w3_e2, geom->b3_e2)  >= 0)&&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {
              window[n]->nzine2++;
              if (cc <= window[n]->enonS) {
                window[n]->nzine2S++;
	      }
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }

    window[n]->zine2 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nzine2 + 1));

    /* Get the interpolation corner locations */
    zc = 1;
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];

      if (window[n]->zoomc[cc] & ZE2) {
        c2 = c;
        for (j = 1; j <= window[n]->zmfe2; j++) {
          c1 = c2;
          for (i = 1; i <= window[n]->zmfe1; i++) {
            lc = window[n]->fm[c1].sc;

            if (!(window[n]->zoomc[lc] & ZE2) &&
		(ANY(c1, geom->w3_e2, geom->b3_e2) >= 0) &&
                (!geom->fm[c1].wn || geom->fm[c1].wn == window[n]->wn)) {

              window[n]->zine2[zc].c = lc;

              /* If cc is an ghost cell find the interior wet cell */
              cnr = cp = cc;
              if (!geom->fm[c].wn) {
                c3 = ANY(cc, window[n]->bpte2, window[n]->nbpte2);
                if (c3 > window[n]->nbe2)
                  cp = window[n]->bine2[c3];
              }
              window[n]->zine2[zc].cnr = cp;
              /* Corner cells are self pointing over tangential */
              /* ghost cell locations only.  */
              cp = window[n]->xp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte2, window[n]->nbpte2);
                if (c3 > window[n]->nbe2)
                  cp = window[n]->bine2[c3];
              }
              window[n]->zine2[zc].ce = cp;
              cp = window[n]->yp1[cnr];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte2, window[n]->nbpte2);
                if (c3 > window[n]->nbe2)
                  cp = window[n]->bine2[c3];
              }
              window[n]->zine2[zc].cn = cp;
              cp = window[n]->yp1[window[n]->xp1[cnr]];
              if (!geom->fm[window[n]->wsa[cp]].wn) {
                c3 = ANY(cp, window[n]->bpte2, window[n]->nbpte2);
                if (c3 > window[n]->nbe2)
                  cp = window[n]->bine2[c3];
              }
              window[n]->zine2[zc].cne = cp;
              window[n]->zine2[zc].x =
                (double)(i - 1) / (double)window[n]->zmfe1;
              window[n]->zine2[zc].y =
                (double)(j - 1) / (double)window[n]->zmfe2;
              zc++;
            }
            c1 = geom->xp1[c1];
          }
          c2 = geom->yp1[c2];
        }
      }
    }
    if (zc - 1 != window[n]->nzine2)
      printf
        ("WARNING : Incorrect number of zoom interpolation e2 faces located, window %d : %d %d\n",
         window[n]->wn, zc - 1, window[n]->nzine2);
  }

  /*-----------------------------------------------------------------*/
  /* Get the global interpolation information.                       */
  /*-----------------------------------------------------------------*/
  /* Count the number of interpolation cells.                        */
  geom->nzine2 = geom->nzine2S = 0;
  /* Set the mask at all valid cells in the master corresponding to  */
  /* the window's zoomed cell to process.                            */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->b3_e2; cc++) {
    c2 = geom->w3_e2[cc];
    wn = geom->fm[c2].wn;
    if (geom->zoomc[geom->m2d[c2]] & (ZE2 | ZE2B)) {
      ze1 = window[wn]->zmfe1;
      ze2 = window[wn]->zmfe2;
      for (zc = 1; zc <= (int)ze1 / 2; zc++)
	c2 = geom->xm1[c2];
      for (zc = 1; zc <= ze2; zc++) {
	c = c2;
	for (zs = 1; zs <= ze1; zs++) {
	  mask[c] = 1;
	  c = geom->xp1[c];
	}
	c2 = geom->yp1[c2];
      }
    }
  }
  /* Count global OBC interpolation cells in the zoom window         */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->to3_e2; cc++) {
      c = open->obc_e2[cc];
      wn = geom->fm[open->nmap[c]].wn;
      if (wn == 0) 
	wn = geom->fm[c].wn;
      zmf = (open->type & U1BDRY) ? window[wn]->zmfe1 : window[wn]->zmfe2;
      if (window[wn]->zoomf > 1) {
	lc = geom->fm[c].sc;
	if (cc <= open->no3_e2)
	  zf = (int)zmf;      /* Normal cells          */
	else
	  zf = (int)zmf / 2;  /* Tangential cells      */
	/* Normal cells on front edges; include an extra row        */
	if (open->ocodey & F_EDGE) {
	  if (mask[c] && !(geom->zoomc[geom->m2d[c]] & ZE2NB)) {
	    geom->nzine2++;
	    if (c <= geom->enonS) {
	      geom->nzine2S++;
	    }
	    mask[c] = 0;
	  }
	  c = open->nmap[c];
	}
	/* Normal and tangential cells                              */
	for (zs = 1; zs <= zf; zs++) {
	  if (mask[c] && !(geom->zoomc[geom->m2d[c]] & (ZE2|ZE2NB|ZE2TB))) {
	    geom->nzine2++;
	    if (c <= geom->enonS) {
	      geom->nzine2S++;
	    }
	    mask[c] = 0;
	  }
	  c = open->nmap[c];
	}
	/* Tangential cells on right edges; include an extra row    */
	if (open->ocodex & R_EDGE) {
	  if (mask[c] && !(window[wn]->zoomc[window[wn]->fm[c].sc] & ZE2)) {
	    geom->nzine2++;
	    if (c <= geom->enonS) {
	      geom->nzine2S++;
	    }
	    mask[c] = 0;
	  }
	}
      }
    }
  }

  /* Count remaining interpolation cells                             */
  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzine2; zc++) {
      c = window[n]->wsa[window[n]->zine2[zc].c];
      if (mask[c]) {
        geom->nzine2++;
        if (c <= geom->enonS)
          geom->nzine2S++;
        mask[c] = 0;

      }
    }
  }

  /* Fill the global interpolation structure */
  geom->zine2 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nzine2 + 1));
  /* Set the mask at all valid cells in the master corresponding to  */
  /* the window's zoomed cell to process.                            */
  memset(mask, 0, (geom->sgsiz) * sizeof(int));
  for (cc = 1; cc <= geom->b3_e2; cc++) {
    c2 = geom->w3_e2[cc];
    wn = geom->fm[c2].wn;
    if (geom->zoomc[geom->m2d[c2]] & (ZE2 | ZE2B)) {
      ze1 = window[wn]->zmfe1;
      ze2 = window[wn]->zmfe2;
      for (zc = 1; zc <= (int)ze1 / 2; zc++)
	c2 = geom->xm1[c2];
      for (zc = 1; zc <= ze2; zc++) {
	c = c2;
	for (zs = 1; zs <= ze1; zs++) {
	  mask[c] = 1;
	  c = geom->xp1[c];
	}
	c2 = geom->yp1[c2];
      }
    }
  }
  geom->nzine2 = geom->nzine2S;
  geom->nzine2S = 0;

  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->to3_e2; cc++) {
      c = open->obc_e2[cc];
      wn = geom->fm[open->nmap[c]].wn;
      if (wn == 0) 
	wn = geom->fm[c].wn;
      zmn = window[wn]->zmfe2;
      zmt = window[wn]->zmfe1;
      if (cc <= open->no3_e2)
	zf = (int)zmn;      /* Normal cells          */
      else
	zf = (int)zmn / 2;  /* Tangential cells      */
      /* Only process cells in zoomed windows that have not been         */
      /* processed (i.e. mask[c] == 1).                                  */
      if (mask[c] && window[wn]->zoomf > 1) {

	/* Normal cells on F_EDGE                                        */
	m = 0;
	if (open->ocodey & F_EDGE) {
	  c2 = c;
	  for (zc = 1; zc <= zmt; zc++) {
	    c2 = open->tmpp[c2];
	    lc = window[wn]->fm[c2].sc;
	    m++;
	    if (geom->fm[c2].wn == wn && geom->zoomc[geom->m2d[c2]] & ZE2NB)
	      break;
	  }
	  ce = cne = c2;
	  cnr = cn = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  c2 = open->obc_e2[cc];
	  y = 0.0; 
	  x = (double)(zmt - m) / (double)zmt;
	  xinc = 1.0 / (double)zmt;
	  for (zc = 1; zc <= m; zc++) {
	    if (!(geom->zoomc[geom->m2d[c2]] & ZE2NB) && 
		mask[c2] && geom->fm[c2].wn == wn) {
	      if (c2 <= geom->enonS) {
		geom->nzine2S++;
		i = geom->nzine2S;
	      } else {
		geom->nzine2++;
		i = geom->nzine2;
	      }
	      geom->zine2[i].cb = geom->bot_e2[c2cc[geom->m2d[c2]]];
	      geom->zine2[i].cnr = cnr;
	      geom->zine2[i].ce = ce;
	      geom->zine2[i].cn = cn;
	      geom->zine2[i].cne = cne;
	      geom->zine2[i].c = c2;
	      geom->zine2[i].x = x;
	      geom->zine2[i].y = y;
	      mask[c2] = 0;
	    }
	    x += xinc;
	    c2 = geom->xp1[c2];
	  }
	}

	/* Find the zoomed grid cell center interior to the boundary     */
	m = 0;
	for (zs = 1; zs <= zf; zs++)
	  c = open->nmap[c];
	for (zc = 1; zc <= zmt; zc++) {
	  c = open->tmpp[c];
	  lc = window[wn]->fm[c].sc;
	  m++;
	  if (geom->fm[c].wn == wn && window[wn]->zoomc[lc] & ZE2) break;
	}
	/* Find the corner points and map pointers                       */
	if (open->ocodex & L_EDGE) {	    
	  cne = cn = c;
	  if (geom->fm[cne].wn != wn)
	    cn = window[wn]->wsa[window[wn]->xm1[lc]];
	  ce = cnr = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  if (geom->fm[ce].wn != wn)
	    cnr = window[wn]->wsa[window[wn]->xm1[lc]];
	  c2 = open->obc_e2[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 1.0 / (double)(zf + 1.0); 
	  xinc = 1.0 / (double)(zf + 1.0);
	  y0 = 0.0; 
	  yinc = 1.0 / (double)zmn;
	  ze1 = zf;
	  ze2 = m;
	}
	if (open->ocodex & R_EDGE) {
	  cn = cne = c;
	  if (geom->fm[cn].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      cne = geom->xp1[cne];
	  }
	  cnr = ce = window[wn]->wsa[(lc = window[wn]->ym1[lc])];
	  if (geom->fm[cnr].wn != wn) {
	    for (zc = 1; zc <= zf; zc++)
	      ce = geom->xp1[ce];
	  }
	  c2 = open->obc_e2[cc];
	  for (zc = 1; zc <= zf; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = 0.0;
	  xinc = 1.0 / (double)(zf + 1.0);
	  y0 = 0.0; 
	  yinc = 1.0 / (double)zmn;
	  ze1 = zf + 1;
	  ze2 = m;
	}
	if (open->ocodey & B_EDGE) {	    
	  cne = c;
	  ce = window[wn]->wsa[window[wn]->ym1[lc]];
	  cn = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  cnr = window[wn]->wsa[window[wn]->ym1[lc]];
	  c2 = open->obc_e2[cc];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = (double)(zmt - m) / (double)zmt;
	  xinc = 1.0 / (double)zmt;
	  y0 = 0.0;
	  yinc = 1.0 / (double)zmn;
	  ze1 = m;
	  ze2 = window[wn]->zoomf;
	}
	if (open->ocodey & F_EDGE) {
	  ce = c;
	  cne = window[wn]->wsa[window[wn]->yp1[lc]];
	  if (geom->fm[cne].wn != wn && ce == cne) {
	    for (zc = 1; zc <= zmn; zc++)
	      cne = geom->yp1[cne];
	  }
	  cnr = window[wn]->wsa[(lc = window[wn]->xm1[lc])];
	  cn = window[wn]->wsa[window[wn]->yp1[lc]];
	  c2 = open->obc_e2[cc];
	  for (zc = 1; zc <= zmn; zc++)
	    c2 = open->nmap[c2];
	  mx = geom->xp1;
	  my = geom->yp1;
	  x0 = (double)(zmt - m) / (double)zmt;
	  xinc = 1.0 / (double)zmt;
	  y = 0.0; 
	  yinc = 1.0 / (double)zmn;
	  ze1 = m;
	  ze2 = zmn;
	}
	/* Fill all global cells bounded by the corners and re-set the   */
	/* mask for these cells.                                         */
	x = x0;
	for (zs = 1; zs <= ze1; zs++) {
	  c = c2;
	  y = y0;
	  for (zc = 1; zc <= ze2; zc++) {
	    if (!(geom->zoomc[geom->m2d[c]] & (ZE2NB | ZE2TB)) && 
		mask[c] && geom->fm[c].wn == wn) {
	      if (c <= geom->enonS) {
		geom->nzine2S++;
		i = geom->nzine2S;
		/*
		print_zoom_bdry("u2b" c, cnr, ce, cn, cne, x, y);
		*/
	      } else {
		geom->nzine2++;
		i = geom->nzine2;
	      }
	      geom->zine2[i].cb = geom->bot_e2[c2cc[geom->m2d[c]]];
	      geom->zine2[i].cnr = cnr;
	      geom->zine2[i].ce = ce;
	      geom->zine2[i].cn = cn;
	      geom->zine2[i].cne = cne;
	      geom->zine2[i].c = c;
	      geom->zine2[i].x = x;
	      geom->zine2[i].y = y;
	      mask[c] = 0;
	    }
	    y += yinc;
	    c = my[c];
	  }
	  x += xinc;
	  c2 = mx[c2];
	}
      }
    }
  }

  for (n = 1; n <= geom->nwindows; n++) {
    for (zc = 1; zc <= window[n]->nzine2; zc++) {
      c = window[n]->wsa[window[n]->zine2[zc].c];

      if (mask[c]) {
        cnr = window[n]->wsa[window[n]->zine2[zc].cnr];
        ce = window[n]->wsa[window[n]->zine2[zc].ce];
        cn = window[n]->wsa[window[n]->zine2[zc].cn];
        cne = window[n]->wsa[window[n]->zine2[zc].cne];
        x = window[n]->zine2[zc].x;
        y = window[n]->zine2[zc].y;

	if (c <= geom->enonS) {
	  geom->nzine2S++;
	  i = geom->nzine2S;
	} else {
	  geom->nzine2++;
	  i = geom->nzine2;
	}
	geom->zine2[i].cb = geom->bot_e2[c2cc[geom->m2d[c]]];
	geom->zine2[i].cnr = cnr;
	geom->zine2[i].ce = ce;
	geom->zine2[i].cn = cn;
	geom->zine2[i].cne = cne;
	geom->zine2[i].c = c;
	geom->zine2[i].x = x;
	geom->zine2[i].y = y;
        mask[c] = 0;
	if (!geom->zine2[i].cb) {
	  c1 = c;
	  while(c1 != geom->zm1[c1])
	    c1 = geom->zm1[c1];
	  geom->zine2[i].cb = geom->zp1[c1];
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Deallocate unrequired structure memory */
  i_free_1d(mask);
  i_free_1d(c2cc);
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->zoomf > 1)
      free((zoom_interp_t *)window[n]->zine2);
    window[n]->nzine2 = window[n]->nzine2S = 0;
  }
}

/* END geom_interp_e2()                                              */
/*-------------------------------------------------------------------*/


void print_zoom_bdry(char *tag, int c, int cnr, int ce, int cn, int cne,
		     double x, double y)
{
  printf("%s %d : %d %d : cnr(%d %d) ce(%d %d) cn(%d %d) cne(%d %d) %f %f\n",
	 tag, c, geom->s2i[c], geom->s2j[c], geom->s2i[cnr], geom->s2j[cnr],
	 geom->s2i[ce], geom->s2j[ce], geom->s2i[cn], geom->s2j[cn], 
	 geom->s2i[cne], geom->s2j[cne], x, y);
}

/*-------------------------------------------------------------------*/
/* Routine to set the arrays for increasing horizontal diffusion     */
/* over the e1 region at fine - coarse grid boundaries.              */
/*-------------------------------------------------------------------*/
void geom_bdry_he1(geometry_t *geom,     /* Sparse global geometery  */
		   geometry_t **window   /* Processing window        */
  )
{
  int c, c1, c2, cc;        /* Sparse coordinate / counter */
  int n;                    /* Window counter */
  int zc;                   /* Interpolation structure counter */
  int c3, lc;
  int zse1, zse2, zn, gc, wn;

  if (geom->zbe1) i_free_1d(geom->zbe1);
  if (geom->zme1) i_free_1d(geom->zme1);

  /*-----------------------------------------------------------------*/
  /* Find the cells that form the boundary between zoomed windows.   */
  /* Initialize                                                      */
  geom->nzbe1 = 0;
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->nzbe1 = window[n]->nzbe1S = 0;
  /* Count the cells that border zoomed grids. Cells are counted in  */
  /* both the fine and coarse grid.                                  */
  for (n = 1; n <= geom->nwindows; n++) {

    if (window[n]->zoomf <= 1)
      continue;

    zse1 = (int)(window[n]->zmfe1 / 2);
    for (cc = 1; cc <= window[n]->b3_e1; cc++) { 
      c = window[n]->w3_e1[cc];  /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->xm1[c];    /* Local cell to west               */
      c2 = window[n]->wsa[c1];   /* Global cell to west              */
      wn = geom->fm[c2].wn;      /* Window number to the west        */
      if(window[n]->zmfe1 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nzbe1 += 1 + laux;
	if (c <= window[n]->enonS) window[n]->nzbe1S += 1 + laux;
	/* Increment the boundary cells in window wn                 */
	lc = geom->fm[c2].sc;
	for (zc = 1; zc <= zse1 + 1 + laux; zc++) {
	  window[wn]->nzbe1 += window[n]->zmfe1;
	  if (lc <= window[wn]->enonS) window[wn]->nzbe1S += window[n]->zmfe1;
	  lc = window[wn]->xp1[lc];
	}
	geom->nzbe1 += 1 + window[n]->zmfe1;
      }
      c1 = window[n]->xp1[c];    /* Local cell to east               */
      c2 = window[n]->wsa[c1];   /* Global cell to east              */
      wn = geom->fm[c2].wn;      /* Window number to the east        */
      if(window[n]->zmfe1 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nzbe1 += 1 + laux;
	if (c <= window[n]->enonS ) window[n]->nzbe1S += 1 + laux;
	/* Increment the boundary cells in window wn                 */
	lc = geom->fm[c2].sc;
	for (zc = 1; zc <= zse1 + 1 + laux; zc++) {
	  window[wn]->nzbe1 += window[n]->zmfe1;
	  if (lc <= window[wn]->enonS) window[wn]->nzbe1S += window[n]->zmfe1;
	  lc = window[wn]->xm1[lc];
	}
	geom->nzbe1 += 1 + window[n]->zmfe1;
      }
    }
  }

  /* Allocate memory                                                 */
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->nzbe1) {
      window[n]->zbe1 = i_alloc_1d(window[n]->nzbe1 + 1);
      window[n]->zme1 = i_alloc_1d(window[n]->nzbe1 + 1);
      window[n]->nzbe1 = window[n]->nzbe1S + 1;
      window[n]->nzbe1S = 1;
    }
  }
  if (geom->nzbe1)
    geom->zbe1 = i_alloc_1d(geom->nzbe1 + 1);
  geom->nzbe1 = geom->nzbe1S = 0;

  /* Fill the arrays with bordering cells                            */
  for (n = 1; n <= geom->nwindows; n++) {

    if (window[n]->zoomf <= 1)
      continue;

    zse1 = (int)(window[n]->zmfe1 / 2);
    zse2 = (int)(window[n]->zmfe2 / 2);
    for (cc = 1; cc <= window[n]->b3_e1; cc++) { 
      c = window[n]->w3_e1[cc];  /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->xm1[c];    /* Local cell to east               */
      c2 = window[n]->wsa[c1];   /* Global cell to east              */
      wn = geom->fm[c2].wn;      /* Window number to the east        */
      if(window[n]->zmfe1 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n. Note include   */
	/* auxiliary cells in the bordering cells.                   */
	lc = c;
	for (zc = 1; zc <= 1 + laux; zc++) {
	  if (lc <= window[n]->enonS) {
	    window[n]->zbe1[window[n]->nzbe1S] = lc;
	    window[n]->zme1[window[n]->nzbe1S] = R_EDGE;
	    window[n]->nzbe1S += 1;
	  } else {
	    window[n]->zbe1[window[n]->nzbe1] = lc;
	    window[n]->zme1[window[n]->nzbe1] = R_EDGE;
	    window[n]->nzbe1 += 1;
	  }
	  lc = window[n]->xm1[lc];
	}
	/* Increment the bordering cells in window wn                */
	c3 = geom->fm[c2].sc;
	for (zc = 1; zc <= zse2; zc++)
	  c3 = window[wn]->ym1[c3];
	for (zn = 1; zn <= window[n]->zmfe2; zn++) {
	  lc = c3;
	  for (zc = 1; zc <= zse1 + 1 + laux; zc++) {
	    if (lc <= window[wn]->enonS) {
	      window[wn]->zbe1[window[wn]->nzbe1S] = lc;
	      window[wn]->zme1[window[wn]->nzbe1S] = L_EDGE;
	      window[wn]->nzbe1S += 1;
	    } else {
	      window[wn]->zbe1[window[wn]->nzbe1] = lc;
	      window[wn]->zme1[window[wn]->nzbe1] = L_EDGE;
	      window[wn]->nzbe1 += 1;
	    }
	    lc = window[wn]->xp1[lc];
	  }
	  c3 = window[wn]->yp1[c3];
	}
	/* Increment the cells in the master                         */
	lc = c2;
	for (zc = 1; zc <= 1 + window[n]->zmfe1; zc++) {
	  geom->nzbe1 += 1;
	  geom->zbe1[geom->nzbe1] = lc;
	  if (lc <= geom->enonS) geom->nzbe1S += 1;
	  lc = geom->xp1[lc];
	}
      }
      c1 = window[n]->xp1[c];    /* Local cell to west               */
      c2 = window[n]->wsa[c1];   /* Global cell to west              */
      wn = geom->fm[c2].wn;      /* Window number to the west        */
      if(window[n]->zmfe1 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n                 */
	lc = c;
	for (zc = 1; zc <= 1 + laux; zc++) {
	  if (lc <= window[n]->enonS) {
	    window[n]->zbe1[window[n]->nzbe1S] = lc;
	    window[n]->zme1[window[n]->nzbe1S] = L_EDGE;
	    window[n]->nzbe1S += 1;
	  } else {
	    window[n]->zbe1[window[n]->nzbe1] = lc;
	    window[n]->zme1[window[n]->nzbe1] = L_EDGE;
	    window[n]->nzbe1 += 1;
	  }
	  lc = window[n]->xp1[lc];
	}
	/* Increment the bordering cells in window wn                */
	c3 = geom->fm[c2].sc;
	for (zc = 1; zc <= zse2; zc++)
	  c3 = window[wn]->ym1[c3];
	for (zn = 1; zn <= window[n]->zmfe2; zn++) {
	  lc = c3;
	  for (zc = 1; zc <= zse1 + 1 + laux; zc++) {
	    if (lc <= window[wn]->enonS) {
	      window[wn]->zbe1[window[wn]->nzbe1S] = lc;
	      window[wn]->zme1[window[wn]->nzbe1S] = R_EDGE;
	      window[wn]->nzbe1S += 1;
	    } else {
	      window[wn]->zbe1[window[wn]->nzbe1] = lc;
	      window[wn]->zme1[window[wn]->nzbe1] = R_EDGE;
	      window[wn]->nzbe1 += 1;
	    }
	    lc = window[wn]->xm1[lc];
	  }
	  c3 = window[wn]->yp1[c3];
	}
	/* Increment the cells in the master                         */
	lc = c2;
	for (zc = 1; zc <= 1 + window[n]->zmfe1; zc++) {
	  geom->nzbe1 += 1;
	  geom->zbe1[geom->nzbe1] = lc;
	  if (lc <= geom->enonS) geom->nzbe1S += 1;
	  lc = geom->xm1[lc];
	}
      }
    }
  }
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nzbe1--;
    window[n]->nzbe1S--;
  }
}

/* END geom_bdry_he1()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for increasing horizontal diffusion     */
/* over the e2 region at fine - coarse grid boundaries.              */
/*-------------------------------------------------------------------*/
void geom_bdry_he2(geometry_t *geom,     /* Sparse global geometery  */
		   geometry_t **window   /* Processing window        */
  )
{
  int c, c1, c2, cc;        /* Sparse coordinate / counter */
  int n;                    /* Window counter */
  int zc;                   /* Interpolation structure counter */
  int c3, lc;
  int zse1, zse2, zn, gc, wn;

  if (geom->zbe2) i_free_1d(geom->zbe2);
  if (geom->zme2) i_free_1d(geom->zme2);

  /*-----------------------------------------------------------------*/
  /* Find the cells that form the boundary between zoomed windows.   */
  /* Initialize                                                      */
  geom->nzbe2 = 0;
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->nzbe2 = window[n]->nzbe2S = 0;
  /* Count the cells that border zoomed grids. Cells are counted in  */
  /* both the fine and coarse grid.                                  */
  for (n = 1; n <= geom->nwindows; n++) {

    if (window[n]->zoomf <= 1)
      continue;

    zse1 = (int)(window[n]->zmfe1 / 2);
    zse2 = (int)(window[n]->zmfe2 / 2);
    /* Count the wet cells                                           */
    for (cc = 1; cc <= window[n]->b3_e2; cc++) { 
      c = window[n]->w3_e2[cc];  /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->ym1[c];    /* Local cell to south              */
      c2 = window[n]->wsa[c1];   /* Global cell to south             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      if(window[n]->zmfe2 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nzbe2 += 1 + laux;
	if (c <= window[n]->enonS) window[n]->nzbe2S += 1 + laux;
	/* Increment the boundary cells in window wn                 */
	lc = geom->fm[c2].sc;
	for (zc = 1; zc <= zse2 + 1 + laux; zc++) {
	  window[wn]->nzbe2 += window[n]->zmfe2;
	  if (lc <= window[wn]->enonS) window[wn]->nzbe2S += window[n]->zmfe2;
	  lc = window[wn]->yp1[lc];
	}
	geom->nzbe2 += 1 + window[n]->zmfe2;
      }
      c1 = window[n]->yp1[c];    /* Local cell to north              */
      c2 = window[n]->wsa[c1];   /* Global cell to north             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      if(window[n]->zmfe2 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nzbe2 += 1 + laux;
	if (c <= window[n]->enonS ) window[n]->nzbe2S += 1 + laux;
	/* Increment the boundary cells in window wn                 */
	lc = geom->fm[c2].sc;
	for (zc = 1; zc <= zse2 + 1 + laux; zc++) {
	  window[wn]->nzbe2 += window[n]->zmfe2;
	  if (lc <= window[wn]->enonS) window[wn]->nzbe2S += window[n]->zmfe2;
	  lc = window[wn]->ym1[lc];
	}
	geom->nzbe2 += 1 + window[n]->zmfe2;
      }
    }
  }
  /* Allocate memory                                                 */
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->nzbe2) {
      window[n]->zbe2 = i_alloc_1d(window[n]->nzbe2 + 1);
      window[n]->zme2 = i_alloc_1d(window[n]->nzbe2 + 1);
      window[n]->nzbe2 = window[n]->nzbe2S + 1;
      window[n]->nzbe2S = 1;
    }
  }
  if (geom->nzbe2)
    geom->zbe2 = i_alloc_1d(geom->nzbe2 + 1);
  geom->nzbe2 = 0;

  /* Fill the arrays with bordering cells                            */
  for (n = 1; n <= geom->nwindows; n++) {

    if (window[n]->zoomf <= 1)
      continue;

    zse1 = (int)(window[n]->zmfe1 / 2);
    zse2 = (int)(window[n]->zmfe2 / 2);
    for (cc = 1; cc <= window[n]->b3_e2; cc++) { 
      c = window[n]->w3_e2[cc];  /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->ym1[c];    /* Local cell to south              */
      c2 = window[n]->wsa[c1];   /* Global cell to south             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      if(window[n]->zmfe2 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n. Note include   */
	/* auxiliary cells in the bordering cells.                   */
	lc = c;
	for (zc = 1; zc <= 1 + laux; zc++) {
	  if (lc <= window[n]->enonS) {
	    window[n]->zbe2[window[n]->nzbe2S] = lc;
	    window[n]->zme2[window[n]->nzbe2S] = F_EDGE;
	    window[n]->nzbe2S += 1;
	  } else {
	    window[n]->zbe2[window[n]->nzbe2] = lc;
	    window[n]->zme2[window[n]->nzbe2] = F_EDGE;
	    window[n]->nzbe2 += 1;
	  }
	  lc = window[n]->ym1[lc];
	}
	/* Increment the bordering cells in window wn                */
	c3 = geom->fm[c2].sc;
	for (zc = 1; zc <= zse1; zc++)
	  c3 = window[wn]->xm1[c3];
	for (zn = 1; zn <= window[n]->zmfe1; zn++) {
	  lc = c3;
	  for (zc = 1; zc <= zse2 + 1 + laux; zc++) {
	    if (lc <= window[wn]->enonS) {
	      window[wn]->zbe2[window[wn]->nzbe2S] = lc;
	      window[wn]->zme2[window[wn]->nzbe2S] = B_EDGE;
	      window[wn]->nzbe2S += 1;
	    } else {
	      window[wn]->zbe2[window[wn]->nzbe2] = lc;
	      window[wn]->zme2[window[wn]->nzbe2] = B_EDGE;
	      window[wn]->nzbe2 += 1;
	    }
	    lc = window[wn]->yp1[lc];
	  }
	  c3 = window[wn]->xp1[c3];
	}
	/* Increment the cells in the master                         */
	lc = c2;
	for (zc = 1; zc <= 1 + window[n]->zmfe2; zc++) {
	  geom->nzbe2 += 1;
	  geom->zbe2[geom->nzbe2] = lc;
	  if (lc <= geom->enonS) geom->nzbe2S += 1;
	  lc = geom->yp1[lc];
	}
      }
      c1 = window[n]->yp1[c];    /* Local cell to north              */
      c2 = window[n]->wsa[c1];   /* Global cell to north             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      if(window[n]->zmfe2 > 1 && geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n                 */
	lc = c;
	for (zc = 1; zc <= 1 + laux; zc++) {
	  if (lc <= window[n]->enonS) {
	    window[n]->zbe2[window[n]->nzbe2S] = lc;
	    window[n]->zme2[window[n]->nzbe2S] = B_EDGE;
	    window[n]->nzbe2S += 1;
	  } else {
	    window[n]->zbe2[window[n]->nzbe2] = lc;
	    window[n]->zme2[window[n]->nzbe2] = B_EDGE;
	    window[n]->nzbe2 += 1;
	  }
	  lc = window[n]->yp1[lc];
	}
	/* Increment the bordering cells in window wn                */
	c3 = geom->fm[c2].sc;
	for (zc = 1; zc <= zse1; zc++)
	  c3 = window[wn]->xm1[c3];
	for (zn = 1; zn <= window[n]->zmfe1; zn++) {
	  lc = c3;
	  for (zc = 1; zc <= zse2 + 1 + laux; zc++) {
	    if (lc <= window[wn]->enonS) {
	      window[wn]->zbe2[window[wn]->nzbe2S] = lc;
	      window[wn]->zme2[window[wn]->nzbe2S] = F_EDGE;
	      window[wn]->nzbe2S += 1;
	    } else {
	      window[wn]->zbe2[window[wn]->nzbe2] = lc;
	      window[wn]->zme2[window[wn]->nzbe2] = F_EDGE;
	      window[wn]->nzbe2 += 1;
	    }
	    lc = window[wn]->ym1[lc];
	  }
	  c3 = window[wn]->xp1[c3];
	}
	/* Increment the cells in the master                         */
	lc = c2;
	for (zc = 1; zc <= 1 + window[n]->zmfe2; zc++) {
	  geom->nzbe2 += 1;
	  geom->zbe2[geom->nzbe2] = lc;
	  if (lc <= geom->enonS) geom->nzbe2S += 1;
	  lc = geom->ym1[lc];
	}
      }
    }
  }
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nzbe2--;
    window[n]->nzbe2S--;
  }
}

/* END geom_bdry_he2()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the fine - coarse grid boundary arrays.            */
/* The cell at the interface is stored in zice1[].c                  */
/* The auxiliary cell over the interface is stored in zice1[].cn     */
/* The auxiliary cell in the global array (geom) at e1 and e2 faces  */
/* are stored in zice1[].ce and zice1[].cne respectively.            */
/*-------------------------------------------------------------------*/
void geom_bdry_e1(geometry_t *geom,      /* Sparse global geometery  */
		  geometry_t **window    /* Processing window        */
  )
{
  int c, c1, c2, cc;        /* Sparse coordinate / counter */
  int n;                    /* Window counter */
  int gc, wn, ce1, ce2;
  int zc, zse1, zse2;

  /*-----------------------------------------------------------------*/
  /* Find the cells that form the boundary between zoomed windows.   */
  /* Initialize                                                      */
  geom->nice1 = 0;
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->nice1 = window[n]->nice1S = 0;

  /* Count the cells that border zoomed grids. Cells are counted in  */
  /* both the fine and coarse grid.                                  */
  for (n = 1; n <= geom->nwindows; n++) {
    for (cc = 1; cc <= window[n]->b3_t; cc++) { 
      c = window[n]->w3_t[cc];   /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->xm1[c];    /* Local cell to west               */
      c2 = window[n]->wsa[c1];   /* Global cell to west              */
      wn = geom->fm[c2].wn;      /* Window number to the west        */
      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nice1 += 1;
	geom->nice1 += 1;
	if (c1 <= window[n]->enonS) window[n]->nice1S += 1;
	if (c2 <= geom->enonS) geom->nice1S += 1;
      }
      c1 = window[n]->xp1[c];    /* Local cell to east               */
      c2 = window[n]->wsa[c1];   /* Global cell to east              */
      wn = geom->fm[c2].wn;      /* Window number to the east        */
      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nice1 += 1;
	geom->nice1 += 1;
	if (c1 <= window[n]->enonS ) window[n]->nice1S += 1;
	if (c2 <= geom->enonS) geom->nice1S += 1;
      }
    }
  }

  /* Allocate memory                                                 */
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->nice1)
      window[n]->zice1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
						 (window[n]->nice1 + 1));
    window[n]->nice1 = 1;
  }
  if (geom->nice1) {
    geom->zice1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
					  (geom->nice1 + 1));
  }
  geom->nice1 = geom->nice1S + 1;
  geom->nice1S = 1;

  /* Fill the arrays with bordering cells                            */
  for (n = 1; n <= geom->nwindows; n++) {
    zse1 = (int)(window[n]->zmfe1 / 2);
    zse2 = (int)(window[n]->zmfe2 / 2);
    for (cc = 1; cc <= window[n]->b3_t; cc++) { 
      c = window[n]->w3_t[cc];   /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->xm1[c];    /* Local cell to east               */
      c2 = window[n]->wsa[c1];   /* Global cell to east              */
      wn = geom->fm[c2].wn;      /* Window number to the east        */
      ce1 = ce2 = c2;
      for (zc = 1; zc <= zse1; zc++)
	ce1 = geom->xm1[ce1];
      for (zc = 1; zc <= zse2; zc++)
	ce2 = geom->ym1[ce2];

      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n.                */
	window[n]->zice1[window[n]->nice1].c = c1;
	window[n]->zice1[window[n]->nice1].cn = c;
	window[n]->nice1++;
	if (c2 <= geom->enonS) {
	  geom->zice1[geom->nice1S].c = c2;
	  geom->zice1[geom->nice1S].cn = gc;
	  geom->zice1[geom->nice1S].ce = ce1;
	  geom->zice1[geom->nice1S].cne = ce2;
	  geom->nice1S++;
	} else {
	  geom->zice1[geom->nice1].c = c2;
	  geom->zice1[geom->nice1].cn = gc;
	  geom->zice1[geom->nice1].ce = ce1;
	  geom->zice1[geom->nice1].cne = ce2;
	  geom->nice1++;
	}
      }
      c1 = window[n]->xp1[c];    /* Local cell to west               */
      c2 = window[n]->wsa[c1];   /* Global cell to west              */
      wn = geom->fm[c2].wn;      /* Window number to the west        */
      ce1 = ce2 = c2;
      for (zc = 1; zc <= zse1; zc++)
	ce1 = geom->xm1[ce1];
      for (zc = 1; zc <= zse2; zc++)
	ce2 = geom->ym1[ce2];

      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n                 */
	window[n]->zice1[window[n]->nice1].c = c1;
	window[n]->zice1[window[n]->nice1].cn = c;
	window[n]->nice1++;
	if (c2 <= geom->enonS) {
	  geom->zice1[geom->nice1S].c = c2;
	  geom->zice1[geom->nice1S].cn = gc;
	  geom->zice1[geom->nice1S].ce = ce1;
	  geom->zice1[geom->nice1S].cne = ce2;
	  geom->nice1S++;
	} else {
	  geom->zice1[geom->nice1].c = c2;
	  geom->zice1[geom->nice1].cn = gc;
	  geom->zice1[geom->nice1].ce = ce1;
	  geom->zice1[geom->nice1].cne = ce2;
	  geom->nice1++;
	}
      }
    }
    window[n]->nice1--;
  }
  geom->nice1--;
  geom->nice1S--;
}

/* END geom_bdry_e1()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the fine auxiliary cells used in the coarse grid.  */
/* The auxiliary cell center is stored in zipe1[].c                  */
/* The auxiliary cell e1 face is stored in zipe1[].ce                */
/* The auxiliary cell e2 face is stored in zipe1[].cn                */
/* The variable zipe1[].x contains a code regarding the nature of    */
/* the e1 cell:                                                      */
/* zipe1[].x = 1 -> cell is adjacent to an e2 (normal) boundary      */
/* zipe1[].x = 0 -> cell is adjacent to an e1 (tangential) boundary  */
/* zipe1[].x = -1 -> cell is adjacent to a solid boundary            */
/* The variable zipe1[].y contains a code regarding the nature of    */
/* the e2 cell:                                                      */
/* zipe1[].y = 1 -> cell is adjacent to an e1 (normal) boundary      */
/* zipe1[].y = 0 -> cell is adjacent to an e2 (tangential) boundary  */
/* zipe1[].y = -1 -> cell is adjacent to a solid boundary            */
/*-------------------------------------------------------------------*/
void geom_m2s(geometry_t *geom,          /* Sparse global geometery  */
	      geometry_t **window        /* Processing window        */
  )
{
  int c, lc, cc;            /* Sparse coordinate / counter */
  int n;                    /* Window counter */
  int wn, ce1, ce2;
  int xm1, xp1, yp1, ym1;

  /*-----------------------------------------------------------------*/
  /* Find the cells that form the boundary between zoomed windows.   */
  /* Initialize                                                      */
  geom->nipe1 = 0;
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->nipe1 = window[n]->nipe1S = 0;

  /* Count the cells that border zoomed grids. Cells are counted in  */
  /* the finest grid only.                                           */
  for (n = 1; n <= geom->nwindows; n++) {
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      c = window[n]->wsa[lc];
      wn = geom->fm[c].wn;
      if (window[wn]->zmfe1 < window[n]->zmfe1 || window[wn]->zmfe2 < window[n]->zmfe2) {
	geom->nipe1++;
	if (c <= geom->enonS) geom->nipe1S += 1;
      }
    }
  }

  /* Allocate memory                                                 */
  if (geom->nipe1) {
    geom->zipe1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
					  (geom->nipe1 + 1));
  }
  geom->nipe1 = geom->nipe1S + 1;
  geom->nipe1S = 1;

  /* Fill the cell locations                                         */
  for (n = 1; n <= geom->nwindows; n++) {
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      xm1 = window[n]->xm1[lc];
      xp1 = window[n]->xp1[lc];
      ym1 = window[n]->ym1[lc];
      yp1 = window[n]->yp1[lc];
      c = window[n]->wsa[lc];
      ce1 = window[n]->m2se1[cc];
      ce2 = window[n]->m2se2[cc];

      wn = geom->fm[c].wn;
      if (window[wn]->zmfe1 < window[n]->zmfe1 || window[wn]->zmfe2 < window[n]->zmfe2) {
	if (c <= geom->enonS) {
	  geom->zipe1[geom->nipe1S].c = c;
	  geom->zipe1[geom->nipe1S].ce = ce1;
	  geom->zipe1[geom->nipe1S].cn = ce2;

 	  if (geom->fm[window[n]->wsa[xm1]].wn == n || 
	      geom->fm[window[n]->wsa[xp1]].wn == n) 
	    geom->zipe1[geom->nipe1S].x = 1;
	  else if (geom->fm[window[n]->wsa[xm1]].wn == 0 || 
	      geom->fm[window[n]->wsa[xp1]].wn == 0) 
	    geom->zipe1[geom->nipe1S].x = -1;
	  else 
	    geom->zipe1[geom->nipe1S].x = 0;

	  if (geom->fm[window[n]->wsa[ym1]].wn == n || 
	      geom->fm[window[n]->wsa[yp1]].wn == n) 
	    geom->zipe1[geom->nipe1S].y = 1;
	  else if (geom->fm[window[n]->wsa[ym1]].wn == 0 || 
	      geom->fm[window[n]->wsa[yp1]].wn == 0) 
	    geom->zipe1[geom->nipe1S].y = -1;
	  else 
	    geom->zipe1[geom->nipe1S].y = 0;

	  geom->nipe1S += 1;
	} else {
	  geom->zipe1[geom->nipe1].c = c;
	  geom->zipe1[geom->nipe1].ce = ce1;
	  geom->zipe1[geom->nipe1].cn = ce2;

	  if (geom->fm[window[n]->wsa[xm1]].wn == n || 
	      geom->fm[window[n]->wsa[xp1]].wn == n) 
	    geom->zipe1[geom->nipe1].x = 1;
	  else if (geom->fm[window[n]->wsa[xm1]].wn == 0 || 
	      geom->fm[window[n]->wsa[xp1]].wn == 0) 
	    geom->zipe1[geom->nipe1].x = -1;
	  else 
	    geom->zipe1[geom->nipe1].x = 0;

	  if (geom->fm[window[n]->wsa[ym1]].wn == n || 
	      geom->fm[window[n]->wsa[yp1]].wn == n) 
	    geom->zipe1[geom->nipe1].y = 1;
	  else if (geom->fm[window[n]->wsa[ym1]].wn == 0 || 
	      geom->fm[window[n]->wsa[yp1]].wn == 0) 
	    geom->zipe1[geom->nipe1].y = -1;
	  else 
	    geom->zipe1[geom->nipe1].y = 0;

	  geom->nipe1++;
	}
      }
    }
  }
  geom->nipe1--;
  geom->nipe1S--;
}

/* END geom_m2s() ()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the fine - coarse grid boundary arrays.            */
/* The cell at the interface is stored in zice2[].c                  */
/* The auxiliary cell over the interface is stored in zice2[].cn     */
/*-------------------------------------------------------------------*/
void geom_bdry_e2(geometry_t *geom,      /* Sparse global geometery  */
		   geometry_t **window   /* Processing window        */
  )
{
  int c, c1, c2, cc;        /* Sparse coordinate / counter */
  int n;                    /* Window counter */
  int gc, wn, ce1, ce2;
  int zc, zse1, zse2;

  /*-----------------------------------------------------------------*/
  /* Find the cells that form the boundary between zoomed windows.   */
  /* Initialize                                                      */
  geom->nice2 = 0;
  for (n = 1; n <= geom->nwindows; n++)
    window[n]->nice2 = window[n]->nice2S = 0;

  /* Count the cells that border zoomed grids. Cells are counted in  */
  /* both the fine and coarse grid.                                  */
  for (n = 1; n <= geom->nwindows; n++) {
    for (cc = 1; cc <= window[n]->b3_t; cc++) { 
      c = window[n]->w3_t[cc];   /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->ym1[c];    /* Local cell to south              */
      c2 = window[n]->wsa[c1];   /* Global cell to south             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nice2 += 1;
	geom->nice2 += 1;
	if (c1 <= window[n]->enonS) window[n]->nice2S += 1;
	if (c2 <= geom->enonS) geom->nice2S += 1;
      }
      c1 = window[n]->yp1[c];    /* Local cell to north              */
      c2 = window[n]->wsa[c1];   /* Global cell to north             */
      wn = geom->fm[c2].wn;      /* Window number to the north       */
      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the boundary cells in window n                  */
	window[n]->nice2 += 1;
	geom->nice2 += 1;
	if (c1 <= window[n]->enonS ) window[n]->nice2S += 1;
	if (c2 <= geom->enonS) geom->nice2S += 1;
      }
    }
  }

  /* Allocate memory                                                 */
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->nice2)
      window[n]->zice2 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
						 (window[n]->nice2 + 1));
    window[n]->nice2 = 1;
  }
  if (geom->nice2) {
    geom->zice2 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
					  (geom->nice2 + 1));
  }
  geom->nice2 = geom->nice2S + 1;
  geom->nice2S = 1;

  /* Fill the arrays with bordering cells                            */
  for (n = 1; n <= geom->nwindows; n++) {
    zse1 = (int)(window[n]->zmfe1 / 2);
    zse2 = (int)(window[n]->zmfe2 / 2);
    for (cc = 1; cc <= window[n]->b3_t; cc++) { 
      c = window[n]->w3_t[cc];   /* Local sparse coordinate          */
      gc = window[n]->wsa[c];    /* Global sparse coordinate         */
      c1 = window[n]->ym1[c];    /* Local cell to south              */
      c2 = window[n]->wsa[c1];   /* Global cell to south             */
      wn = geom->fm[c2].wn;      /* Window number to the south       */
      ce1 = ce2 = c2;
      for (zc = 1; zc <= zse1; zc++)
	ce1 = geom->xm1[ce1];
      for (zc = 1; zc <= zse2; zc++)
	ce2 = geom->ym1[ce2];

      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n.                */
	window[n]->zice2[window[n]->nice2].c = c1;
	window[n]->zice2[window[n]->nice2].cn = c;
	window[n]->nice2++;
	if (c2 <= geom->enonS) {
	  geom->zice2[geom->nice2S].c = c2;
	  geom->zice2[geom->nice2S].cn = gc;
	  geom->zice2[geom->nice2S].ce = ce1;
	  geom->zice2[geom->nice2S].cne = ce2;
	  geom->nice2S++;
	} else {
	  geom->zice2[geom->nice2].c = c2;
	  geom->zice2[geom->nice2].cn = gc;
	  geom->zice2[geom->nice2].ce = ce1;
	  geom->zice2[geom->nice2].cne = ce2;
	  geom->nice2++;
	}
      }
      c1 = window[n]->yp1[c];    /* Local cell to north              */
      c2 = window[n]->wsa[c1];   /* Global cell to north             */
      wn = geom->fm[c2].wn;      /* Window number to the north       */
      ce1 = ce2 = c2;
      for (zc = 1; zc <= zse1; zc++)
	ce1 = geom->xm1[ce1];
      for (zc = 1; zc <= zse2; zc++)
	ce2 = geom->ym1[ce2];

      if(geom->fm[gc].wn == n && wn > 0 && wn != n) {
	/* Increment the bordering cells in window n                 */
	window[n]->zice2[window[n]->nice2].c = c1;
	window[n]->zice2[window[n]->nice2].cn = c;
	window[n]->nice2++;
	if (c2 <= geom->enonS) {
	  geom->zice2[geom->nice2S].c = c2;
	  geom->zice2[geom->nice2S].cn = gc;
	  geom->zice2[geom->nice2S].ce = ce1;
	  geom->zice2[geom->nice2S].cne = ce2;
	  geom->nice2S++;
	} else {
	  geom->zice2[geom->nice2].c = c2;
	  geom->zice2[geom->nice2].cn = gc;
	  geom->zice2[geom->nice2].ce = ce1;
	  geom->zice2[geom->nice2].cne = ce2;
	  geom->nice2++;
	}
      }
    }
    window[n]->nice2--;
  }
  geom->nice2--;
  geom->nice2S--;
}

/* END geom_bdry_e2()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for interpolation of zoomed variables   */
/* variables on e1 faces in the master. Over ghost cells the corner  */
/* interpolation calls are only self pointing for tangential faces,  */
/* normal e1 faces use the (zero valued) ghost cell.                 */
/*-------------------------------------------------------------------*/
void geom_sinterp_e1(geometry_t *geom, /* Sparse global geometery */
		     geometry_t **window /* Processing window */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int n;                        /* Window counter */
  int z1, z2;                   /* Interpolation structure counter */
  int *mask;                    /* Duumy mask */
  int *c2cc;
  int i, j;
  int lc;
  int zse1, zse2;
  int cnr, ce, cn, cs;
  double x, y, x1, x2;
  
  if (!geom->zoomf)
    return;
  mask = i_alloc_1d(geom->enon + 1);
  c2cc = i_alloc_1d(geom->sgsizS);
  
  /* Allocate the c to cc array */
  i = 0;
  for (cc = 1; cc <= geom->b2_e1; cc++) {
    if (geom->w2_e1[cc] > i)
      i = geom->w2_e1[cc];
  }
  c2cc = i_alloc_1d(i + 1);
  for (cc = 1; cc <= geom->b2_e1; cc++) {
    c = geom->w2_e1[cc];
    c2cc[c] = cc;
  }

  /* Set the mask adjacent to southern solid boundaries */
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  for (cc = 1; cc <= geom->v3_e1; cc++) {
    c = geom->w3_e1[cc];      /* Global cell to process */
    cn = geom->xm1[c];        /* Global cell to the west of c */
    if (geom->fm[cn].wn && !(geom->fm[geom->xm1[cn]].wn))
      mask[cn] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the local interpolation information */
  /* Count the number of non-center zoom wet cells */
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nipe1 = window[n]->nipe1S = 0;
    window[n]->nice1 = window[n]->nice1S = 0;
    if (window[n]->zoomf == 1)
      continue;
    /* Count the number of parallel cells and allocate the structure */
    zse1 = window[n]->zmfe1 - 1;
    window[n]->nipe1 = window[n]->v3_t * zse1;
    window[n]->nipe1S = window[n]->v2_t * zse1;
    window[n]->zipe1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nipe1 + 1));
    /* Count the number of cross cells and allocate the structure */
    zse1 = window[n]->zmfe1 + window[n]->zmfe2;
    window[n]->nice1 = window[n]->v3_t * zse1;
    window[n]->nice1S = window[n]->v2_t * zse1;
    window[n]->zice1 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nice1 + 1));

    z1 = z2 = 1;
    zse1 = window[n]->zmfe1 / 2;
    zse2 = window[n]->zmfe2 / 2;
    for (cc = 1; cc <= window[n]->v3_t; cc++) {
      lc = window[n]->w3_t[cc];
      cnr = ce = window[n]->wsa[lc];
      for (i = 1; i <= zse1 ; i++) {
	cnr = geom->xm1[cnr];
	ce = geom->xp1[ce];
      }
      ce = geom->xp1[ce];

      /* Get the parallel interpolation locations */
      c = cnr;
      x = 0.0;
      for (i = 1; i < window[n]->zmfe1; i++) {
	c = geom->xp1[c];
	x += (1.0 / (double)window[n]->zmfe1);
	window[n]->zipe1[z1].c = c;
	window[n]->zipe1[z1].cnr = cnr;
	window[n]->zipe1[z1].ce = ce;
	window[n]->zipe1[z1].x = x;
	z1++;
      }

      /* Get the cross interpolation locations */
      for (j = 1; j <= window[n]->zmfe1; j++) {
	cs = cn = cnr;
	x1 = x2 = 0.0;
	for (i = 1; i <= window[n]->zmfe2; i++) {
	  cn = geom->yp1[cn];
	  x1 += 1.0;
	  if (geom->fm[cn].wn != n) break;
	}
	for (i = 1; i <= window[n]->zmfe2; i++) {
	  cs = geom->ym1[cs];
	  x2 += 1.0;
	  if (geom->fm[cs].wn != n) break;
	}
	c = cnr;
	y = 0.0;
	for (i = 1; i <= zse2; i++) {
	  c = geom->yp1[c];
	  y += (1.0 / x1);
	  window[n]->zice1[z2].c = c;
	  window[n]->zice1[z2].cnr = cnr;
	  window[n]->zice1[z2].cn = cn;
	  window[n]->zice1[z2].y = y;
	  z2++;
	}
	c = cnr;
	y = 1.0;
	for (i = 1; i <= zse2; i++) {
	  c = geom->ym1[c];
	  y -= (1.0 / x2);
	  window[n]->zice1[z2].c = c;
	  window[n]->zice1[z2].cnr = cs;
	  window[n]->zice1[z2].cn = cnr;
	  window[n]->zice1[z2].y = y;
	  z2++;
	}
	cnr = geom->xp1[cnr];
      }
    }
    if (z1 - 1 != window[n]->nipe1)
      printf
        ("WARNING : Incorrect number of parallel zoom interpolation e1 faces located, window %d : %d %d\n",
         window[n]->wn, z1 - 1, window[n]->nipe1);
    if (z2 - 1 != window[n]->nice1)
      printf
        ("WARNING : Incorrect number of cross zoom interpolation e1 faces located, window %d : %d %d\n",
         window[n]->wn, z2 - 1, window[n]->nice1);
  }

  /*-----------------------------------------------------------------*/
  /* Get the global interpolation information */
  /* Count the number of interpolation cells */
  geom->nzine1 = geom->nzine1S = 0;
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nipe1; z1++) {
      c = window[n]->zipe1[z1].c;
      if (!mask[c]) {
        geom->nipe1++;
        if (c <= geom->enonS)
          geom->nipe1S++;
        mask[c] = 1;
      }
    }
    for (z2 = 1; z2 <= window[n]->nice1; z2++) {
      c = window[n]->zice1[z2].c;
      if (!mask[c]) {
        geom->nice1++;
        if (c <= geom->enonS)
          geom->nice1S++;
        mask[c] = 1;
      }
    }
  }

  /* Fill the global interpolation structure */
  geom->zipe1 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nipe1 + 1));
  geom->zice1 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nice1 + 1));
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  geom->nipe1 = 1;
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nipe1; z1++) {
      c = window[n]->zipe1[z1].c;
      if (!mask[c]) {
        geom->zipe1[geom->nipe1].c = c;
	geom->zipe1[geom->nipe1].cb = geom->bot_e1[c2cc[c]];
        cnr = window[n]->zipe1[z1].cnr;
        ce = window[n]->zipe1[z1].ce;
        x = window[n]->zipe1[z1].x;

        geom->zipe1[geom->nipe1].cnr = cnr;
        geom->zipe1[geom->nipe1].ce = ce;
        geom->zipe1[geom->nipe1].x = x;
        geom->nipe1++;
        mask[c] = 1;
      }

    }
  }
  geom->nipe1--;

  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  geom->nice1 = 1;
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nice1; z1++) {
      c = window[n]->zice1[z1].c;
      if (!mask[c]) {
        geom->zice1[geom->nice1].c = c;
	geom->zice1[geom->nice1].cb = geom->bot_e1[c2cc[c]];
        cnr = window[n]->zice1[z1].cnr;
        cn = window[n]->zice1[z1].cn;
        y = window[n]->zice1[z1].y;

        geom->zice1[geom->nice1].cnr = cnr;
        geom->zice1[geom->nice1].cn = cn;
        geom->zice1[geom->nice1].y = y;
        geom->nice1++;
        mask[c] = 1;
      }

    }
  }
  geom->nice1--;

  /*-----------------------------------------------------------------*/
  /* Deallocate unrequired structure memory */
  i_free_1d(mask);
  i_free_1d(c2cc);
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->zoomf > 1) {
      free((zoom_interp_t *)window[n]->zipe1);
      free((zoom_interp_t *)window[n]->zice1);
    }
    window[n]->nipe1 = window[n]->nice1S = 0;
  }
}

/* END geom_sinterp_e1()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the arrays for interpolation of zoomed variables   */
/* variables on e2 faces in the master. Over ghost cells the corner  */
/* interpolation calls are only self pointing for tangential faces,  */
/* normal e2 faces use the (zero valued) ghost cell.                 */
/*-------------------------------------------------------------------*/
void geom_sinterp_e2(geometry_t *geom, /* Sparse global geometery */
		     geometry_t **window /* Processing window */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int n;                        /* Window counter */
  int z1, z2;                   /* Interpolation structure counter */
  int *mask;                    /* Duumy mask */
  int *c2cc;
  int i, j;
  int lc;
  int zse1, zse2;
  int cnr, ce, cw, cn;
  double x, y, x1, x2;

  if (!geom->zoomf)
    return;
  mask = i_alloc_1d(geom->enon + 1);
  c2cc = i_alloc_1d(geom->sgsizS);

  /* Allocate the c to cc array */
  i = 0;
  for (cc = 1; cc <= geom->b2_e2; cc++) {
    if (geom->w2_e2[cc] > i)
      i = geom->w2_e2[cc];
  }
  c2cc = i_alloc_1d(i + 1);
  for (cc = 1; cc <= geom->b2_e2; cc++) {
    c = geom->w2_e2[cc];
    c2cc[c] = cc;
  }

  /* Set the mask adjacent to southern solid boundaries */
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  for (cc = 1; cc <= geom->v3_e2; cc++) {
    c = geom->w3_e2[cc];      /* Global cell to process */
    cn = geom->ym1[c];        /* Global cell to the south of c */
    if (geom->fm[cn].wn && !(geom->fm[geom->ym1[cn]].wn))
      mask[cn] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the local interpolation information */
  /* Count the number of non-center zoom wet cells */
  for (n = 1; n <= geom->nwindows; n++) {
    window[n]->nipe2 = window[n]->nipe2S = 0;
    window[n]->nice2 = window[n]->nice2S = 0;
    if (window[n]->zoomf == 1)
      continue;
    /* Count the number of parallel cells and allocate the structure */
    zse2 = window[n]->zmfe2 - 1;
    window[n]->nipe2 = window[n]->v3_t * zse2;
    window[n]->nipe2S = window[n]->v2_t * zse2;
    window[n]->zipe2 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nipe2 + 1));
    /* Count the number of cross cells and allocate the structure */
    zse2 = window[n]->zmfe2 + window[n]->zmfe1;
    window[n]->nice2 = window[n]->v3_t * zse2;
    window[n]->nice2S = window[n]->v2_t * zse2;
    window[n]->zice2 = (zoom_interp_t *)malloc(sizeof(zoom_interp_t) *
                                               (window[n]->nice2 + 1));

    z1 = z2 = 1;
    zse2 = window[n]->zmfe2 / 2;
    zse1 = window[n]->zmfe1 / 2;
    for (cc = 1; cc <= window[n]->v3_t; cc++) {
      lc = window[n]->w3_t[cc];
      cnr = cn = window[n]->wsa[lc];
      for (i = 1; i <= zse2 ; i++) {
	cnr = geom->ym1[cnr];
	cn = geom->yp1[cn];
      }
      cn = geom->yp1[cn];

      /* Get the parallel interpolation locations */
      c = cnr;
      y = 0.0;
      for (i = 1; i < window[n]->zmfe2; i++) {
	c = geom->yp1[c];
	y += (1.0 / (double)window[n]->zmfe2);
	window[n]->zipe2[z1].c = c;
	window[n]->zipe2[z1].cnr = cnr;
	window[n]->zipe2[z1].cn = cn;
	window[n]->zipe2[z1].y = y;
	z1++;
      }

      /* Get the cross interpolation locations */
      for (j = 1; j <= window[n]->zmfe2; j++) {
	ce = cw = cnr;
	x1 = x2 = 0.0;
	for (i = 1; i <= window[n]->zmfe1; i++) {
	  ce = geom->xp1[ce];
	  x1 += 1.0;
	  if (geom->fm[ce].wn != n) break;
	}
	for (i = 1; i <= window[n]->zmfe1; i++) {
	  cw = geom->xm1[cw];
	  x2 += 1.0;
	  if (geom->fm[cw].wn != n) break;
	}
	c = cnr;
	x = 0.0;
	for (i = 1; i <= zse1; i++) {
	  c = geom->xp1[c];
	  x += (1.0 / x1);
	  window[n]->zice2[z2].c = c;
	  window[n]->zice2[z2].cnr = cnr;
	  window[n]->zice2[z2].ce = ce;
	  window[n]->zice2[z2].x = x;
	  z2++;
	}
	c = cnr;
	x = 1.0;
	for (i = 1; i <= zse1; i++) {
	  c = geom->xm1[c];
	  x -= (1.0 / x2);
	  window[n]->zice2[z2].c = c;
	  window[n]->zice2[z2].cnr = cw;
	  window[n]->zice2[z2].ce = cnr;
	  window[n]->zice2[z2].x = x;
	  z2++;
	}
	cnr = geom->yp1[cnr];
      }
    }
    if (z1 - 1 != window[n]->nipe2)
      printf
        ("WARNING : Incorrect number of parallel zoom interpolation e2 faces located, window %d : %d %d\n",
         window[n]->wn, z1 - 1, window[n]->nipe2);
    if (z2 - 1 != window[n]->nice2)
      printf
        ("WARNING : Incorrect number of cross zoom interpolation e2 faces located, window %d : %d %d\n",
         window[n]->wn, z2 - 1, window[n]->nice2);
  }

  /*-----------------------------------------------------------------*/
  /* Get the global interpolation information */
  /* Count the number of interpolation cells */
  geom->nzine2 = geom->nzine2S = 0;
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nipe2; z1++) {
      c = window[n]->zipe2[z1].c;
      if (!mask[c]) {
        geom->nipe2++;
        if (c <= geom->enonS)
          geom->nipe2S++;
        mask[c] = 1;
      }
    }
    for (z2 = 1; z2 <= window[n]->nice2; z2++) {
      c = window[n]->zice2[z2].c;
      if (!mask[c]) {
        geom->nice2++;
        if (c <= geom->enonS)
          geom->nice2S++;
        mask[c] = 1;
      }
    }
  }

  /* Fill the global interpolation structure */
  geom->zipe2 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nipe2 + 1));
  geom->zice2 =
    (zoom_interp_t *)malloc(sizeof(zoom_interp_t) * (geom->nice2 + 1));
  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  geom->nipe2 = 1;
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nipe2; z1++) {
      c = window[n]->zipe2[z1].c;
      if (!mask[c]) {
        geom->zipe2[geom->nipe2].c = c;
	geom->zipe2[geom->nipe2].cb = geom->bot_e2[c2cc[c]];
        cnr = window[n]->zipe2[z1].cnr;
        cn = window[n]->zipe2[z1].cn;
        y = window[n]->zipe2[z1].y;

        geom->zipe2[geom->nipe2].cnr = cnr;
        geom->zipe2[geom->nipe2].cn = cn;
        geom->zipe2[geom->nipe2].y = y;
        geom->nipe2++;
        mask[c] = 1;
      }

    }
  }
  geom->nipe2--;

  memset(mask, 0, (geom->enon + 1) * sizeof(int));
  geom->nice2 = 1;
  for (n = 1; n <= geom->nwindows; n++) {
    for (z1 = 1; z1 <= window[n]->nice2; z1++) {
      c = window[n]->zice2[z1].c;
      if (!mask[c]) {
        geom->zice2[geom->nice2].c = c;
	geom->zice2[geom->nice2].cb = geom->bot_e2[c2cc[c]];
        cnr = window[n]->zice2[z1].cnr;
        ce = window[n]->zice1[z1].ce;
        x = window[n]->zice1[z1].x;

        geom->zice2[geom->nice2].cnr = cnr;
        geom->zice2[geom->nice2].ce = ce;
        geom->zice2[geom->nice2].x = x;
        geom->nice2++;
        mask[c] = 1;
      }

    }
  }
  geom->nice2--;

  /*-----------------------------------------------------------------*/
  /* Deallocate unrequired structure memory */
  i_free_1d(mask);
  i_free_1d(c2cc);
  for (n = 1; n <= geom->nwindows; n++) {
    if (window[n]->zoomf > 1) {
      free((zoom_interp_t *)window[n]->zipe2);
      free((zoom_interp_t *)window[n]->zice2);
    }
    window[n]->nipe2 = window[n]->nice2S = 0;
  }
}

/* END geom_sinterp_e2()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to fill all the values of A in the zoomed cell to A[c]    */
/*-------------------------------------------------------------------*/
void fill_zoom_cell(geometry_t *geom, double *A, int c, int ze1, int ze2)
{
  int cc, n1, n2;
  int c1, cs, cg = c;
  int z1 = (int)(ze1 / 2);
  int z2 = (int)(ze2 / 2);

  cs = geom->m2d[c];
  /* Get the corner location                                         */
  if (geom->zoomc[cs] & (ZC|ZE1))
    for (n2 = 1; n2 <= z2; n2++)
      cg = geom->ym1[cg];
  if (geom->zoomc[cs] & (ZC|ZE2))
    for (n1 = 1; n1 <= z1; n1++)
      cg = geom->xm1[cg];

  /* Fill the array A                                                */
  for (n2 = 1; n2 <= ze2; n2++) {
    c1 = cg;
    for (n1 = 1; n1 <= ze1; n1++) {
      A[c1] = A[c];
      c1 = geom->xp1[c1];
    }
    cg = geom->yp1[cg];
  }
}

/* END fill_zoom_cell()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the 3D - 2D maps for regions where zoomf>1 and   */
/* a step in the bathymetry occurs at a location between zoomed cell */
/* nodes. In this case there exists no corresponding surface cell in */
/* the coarsened grid and this is reset as the nearest coarse cell.  */
/*                                                                   */
/*                        | x | o |   | x |                          */
/*                        |---|---|---|---|   x = wet cell node      */
/*                        |   |   |   |   |  -- = vertical layers    */
/*                        |___|___|   |   |  __ = bottom boundary    */
/*                        |---|---|---|---|                          */
/*                        |   | g |   | x |   g = ghost cell         */
/*                        |   |   |___|___|   o = surface cell above */
/*                        |---|---|---|---|       g.                 */
/*                        |   |   |   |   |                          */
/*                                                                   */
/*-------------------------------------------------------------------*/
void set_zoom_m2d(geometry_t *window  /* Window data structure */
  )
{
  int c, cc;
  int zoomf = window->zoomf;

  if (zoomf == 1)
    return;
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->m2d[cc];

    if ((c <= 0 || c > window->enonS)) {
      /* Get the direction in which a wet cell is adjacent to c */
      c = window->m2d[window->xp1[cc]];
      if (c > 0 && c <= window->enonS)
        window->m2d[cc] = c;
      c = window->m2d[window->xm1[cc]];
      if (c > 0 && c <= window->enonS)
        window->m2d[cc] = c;
      c = window->m2d[window->yp1[cc]];
      if (c > 0 && c <= window->enonS)
        window->m2d[cc] = c;
      c = window->m2d[window->ym1[cc]];
      if (c > 0 && c <= window->enonS)
        window->m2d[cc] = c;
    }
  }
}

/* END set_zoom_m2d()                                                */
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
    aa[c] = cvol1(master, A, c);
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
	     int nctp          /* Size of ctp                        */
  )
{
  geometry_t *geom = master->geom;
  int c, cc;
  double* aa;

  aa = d_alloc_1d(geom->sgsiz);
  for (cc = 1; cc <= nctp; cc++) {
    c = ctp[cc];
    aa[c] = cvol1(master, A, c);
  }
  memcpy(A, aa, geom->sgsiz * sizeof(double));

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
      B[c] = min(lim, cvol1w(window, A, c));
    }
    memcpy(A, B, window->sgsiz * sizeof(double));
    /*
    for (cc = 1; cc <= nctp; cc++) {
      c = ctp[cc];
      A[c] = min(lim, A[c]);
    }
    memcpy(B, A, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= nctp; cc++) {
      c = ctp[cc];
      B[c] = cvol1w(window, A, c);
    }
    memcpy(A, B, window->sgsiz * sizeof(double));
    */
  }
}

/* END smooth_w()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a sponge zone where friction increases towards the coarse    */
/* fine interface.                                                   */
/*-------------------------------------------------------------------*/
void set_sponge_zoom(geometry_t *window,   /* Window geometry        */
		     window_t *windat,     /* Window data            */
		     win_priv_t *wincon,   /* Window constants       */
		     double *u1vh,         /* e1 viscosity           */
		     double *u2vh,         /* e2 viscosity           */
		     int nzbe1,
		     int nzbe2
		     )
{
  int c, cs, cc;                /* Sparse counters                   */
  double sfact = 0.6;           /* Safety factor                     */
  double tfact = 0.3;           /* tanh() multiplier : small = long  */
  int bn, ln, lc, sc;           /* Boundary counter                  */
  double rm, vh;                /* Dummies                           */
  int *mask = wincon->i7;
  int blendf = 1;
  int spgn = 5;
  int *imap;

  memset(mask, 0, window->sgsizS * sizeof(int));

  if(!(wincon->exmapf & E1_INNER)) {
    for (cc = 1; cc <= nzbe1; cc++) {
      c = lc = window->zbe1[cc];
      cs = window->m2d[c];
      imap = (window->zme1[cc] == R_EDGE) ? window->xp1 : window->xm1;
      /* Get viscosity spgn cells from the current cell              */
      for (ln = 1; ln < spgn; ln++)
	lc = imap[lc];
      for (bn = 1; bn <= spgn; bn++) {
	rm = 1.0 / (window->h1au1[cs] * window->h1au1[cs]) +
	  1.0 / (window->h2au1[cs] * window->h2au1[cs]);
	rm = sfact / (4.0 * rm * windat->dt);
	/*rm = 0.01 * window->h1au1[cs] * window->h1au1[cs] / windat->dt;*/
	vh = u1vh[lc];
	/* Linear ramp to rm on the boundary                         */
	/*
	  u1vh[c] = (vh - rm) * (double)(bn - 1) / ((double)spgn - 1.0) + rm;
	*/
	/* Hyperbolic tangent ramp to rm on the boundary             */
	u1vh[c] =  vh + (rm - vh) * (1 - tanh(tfact * (double)(bn - 1)));
	
	mask[window->m2d[c]] |= 1;
	c = imap[c];
      }
    }
  }

  if(!(wincon->exmapf & E2_INNER)) {
    for (cc = 1; cc <= nzbe2; cc++) {
      c = lc = window->zbe2[cc];
      cs = window->m2d[c];
      imap = (window->zme2[cc] == F_EDGE) ? window->yp1 : window->ym1;
      /* Get viscosity spgn cells from the current cell              */
      for (ln = 1; ln < spgn; ln++)
	lc = imap[lc];
      for (bn = 1; bn <= spgn; bn++) {
	rm = 1.0 / (window->h1au2[cs] * window->h1au2[cs]) +
	  1.0 / (window->h2au2[cs] * window->h2au2[cs]);
	rm = sfact / (4.0 * rm * windat->dt);
	/*rm = 0.01 * window->h2au2[cs] * window->h2au2[cs] / windat->dt;*/
	vh = u2vh[c];
	/* Linear ramp to rm on the boundary                    */
	u2vh[c] = (fabs(u2vh[c]) - rm) * 
	  (double)(bn - 1) / (double)spgn + rm;
	/* Hyperbolic tangent ramp to rm on the boundary        */
	u2vh[c] =  vh + (rm - vh) * (1 - tanh(tfact * (double)(bn - 1)));
	mask[window->m2d[c]] |= 2;
	c = imap[c];
      }
    }
  }

  /* Blend the sponge laterally to remove large discontinuties in    */
  /* viscosity.                                                      */
  if (blendf && !(wincon->exmapf & E1_INNER)) {
    for (cc = 1; cc <= nzbe1; cc++) {
      c = window->zbe1[cc];
      cs = window->m2d[c];    
      imap = (window->zme1[cc] == R_EDGE) ? window->xp1 : window->xm1;
      for (bn = 1; bn <= spgn; bn++) {
	vh = u1vh[c];
	/* Blend u1vh to the front (yp1) of current position         */
	lc = sc = window->yp1[c];
	if (!(mask[window->m2d[lc]] & 1)) {
	  /* Get viscosity spgn cells from the current cell          */
	  for (ln = 2; ln <= spgn && !(mask[window->m2d[lc]] & 1); ln++)
	    lc = window->yp1[lc];
	  rm = u1vh[lc];
	  /* Blend to this viscosity                                 */
	  for (ln = 2; ln <= spgn; ln++) {
	    /* Linear                                                */
	    u1vh[sc] = (vh - rm) * (double)(ln - 1) / 
	      (double)(1 - spgn) + vh;
	    /* Hyperbolic tangent                                      */
	    /*
	    u1vh[sc] =  rm + (vh - rm) * (1 - tanh(tfact * (double)(ln - 1)));
	    */
	    mask[window->m2d[sc]] |= 1;
	    sc = window->yp1[sc];
	  }
	}
	/* Blend u1vh to the back (ym1) of current position          */
	lc = sc = window->ym1[c];
	if (!(mask[window->m2d[lc]] & 1)) {
	  for (ln = 2; ln <= spgn && !(mask[window->m2d[lc]] & 1); ln++)
	    lc = window->ym1[lc];
	  rm = u1vh[lc];
	  for (ln = 2; ln <= spgn; ln++) {
	    /* Linear                                                 */
	    u1vh[sc] = (vh - rm) * (double)(ln - 1) / 
	      (double)(1 - spgn) + vh;
	    /* Hyperbolic tangent                                     */
	    /*
	    u1vh[sc] =  rm + (vh - rm) * (1 - tanh(tfact * (double)(ln - 1)));
	    */
	    mask[window->m2d[sc]] |= 1;
	    sc = window->ym1[sc];
	  }
	}
	c = imap[c];
      }
    }
  }

  if (blendf && !(wincon->exmapf & E2_INNER)) {
    for (cc = 1; cc <= nzbe2; cc++) {
      c = window->zbe2[cc];
      cs = window->m2d[c];
      imap = (window->zme2[cc] == F_EDGE) ? window->yp1 : window->ym1;
      for (bn = 1; bn <= spgn; bn++) {
	vh = u2vh[c];
	/* Blend u2vh to the right (xp1) of current position        */
	lc = sc = window->xp1[c];
	if (!(mask[window->m2d[lc]] & 2)) {
	  /* Get viscosity spgn cells from the current cell         */
	  for (ln = 2; ln <= spgn && !(mask[window->m2d[lc]] & 2); ln++)
	    lc = window->xp1[lc];
	  rm = u2vh[lc];
	  /* Blend to this viscosity                                */
	  for (ln = 2; ln <= spgn; ln++) {
	    /* Linear                                               */
	    u2vh[sc] = (vh - rm) * (double)(ln - 1) / 
	      (double)(1 - spgn) + vh;
	    /* Hyperbolic tangent                                   */
	    /*
	    u2vh[sc] = rm + (vh - rm) * (1 - tanh(tfact * (double)(ln - 1)));
	    */
	    mask[window->m2d[sc]] |= 2;
	    sc = window->xp1[sc];
	  }
	}

	/* Blend u2vh to the left (xm1) of current position         */
	lc = sc = window->xm1[c];
	if (!(mask[window->m2d[lc]] & 2)) {
	  for (ln = 2; ln <= spgn && !(mask[window->m2d[lc]] & 2); ln++)
	    lc = window->xm1[lc];
	  rm = u2vh[lc];
	  for (ln = 2; ln <= spgn; ln++) {
	    /* Linear                                               */
	    u2vh[sc] = (vh - rm) * (double)(ln - 1) / 
	      (double)(1 - spgn) + vh;
	    /* Hyperbolic tangent                                   */
	    /*
	    u2vh[sc] =  rm + (vh - rm) * (1 - tanh(tfact * (double)(ln - 1)));
	    */
	    mask[window->m2d[sc]] |= 1;
	    sc = window->xm1[sc];
	  }
	}
	c = imap[c];
      }
    }
  }
}

/* END set_sponge_zoom()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up blending zones                                  */
/*-------------------------------------------------------------------*/
void set_blend_zones(master_t *master, geometry_t *window, win_priv_t *wincon)
{
  int c, c1, cc, tn, i;
  int winsize;

  /*-----------------------------------------------------------------*/
  /* Set the zones in e1 directions                                  */
  if (master->nbl1) {
    wincon->nbl1 = master->nbl1;
    wincon->ble1 = (blend_t **)malloc(sizeof(blend_t *) * master->nbl1);
    winsize = geom->nce2;

    for (tn = 0; tn < wincon->nbl1; tn++) {
      wincon->ble1[tn] = (blend_t *)malloc(sizeof(blend_t));
      blend_t *blend = wincon->ble1[tn];
      if (master->e1bs[tn] >= 0 && master->e1bs[tn] < geom->nce1 && 
	  master->e1be[tn] >= 0 && master->e1be[tn] < geom->nce1) {

	/* Reverse the indicies if required                          */
	if (master->e1bs[tn] > master->e1be[tn]) {
	  c = master->e1be[tn];
	  master->e1be[tn] = master->e1bs[tn];
	  master->e1bs[tn] = c;
	}
	if (master->e2bs[tn] > master->e2be[tn]) {
	  c = master->e2be[tn];
	  master->e2be[tn] = master->e2bs[tn];
	  master->e2bs[tn] = c;
	}

	/* Count the number of cells in the zone                     */
	blend->ne1b = blend->ne1bS = 0;
	blend->nlim = blend->nlimS = 0;
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      window->s2j[c] >= master->e2bs[tn] &&
	      window->s2j[c] <= master->e2be[tn]) {
	    blend->ne1b++;
	    if (c <= window->enonS) blend->ne1bS++;
	  }
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      (window->s2j[c] == master->e2bs[tn] ||
	       window->s2j[c] == master->e2be[tn])) {
	    blend->nlim++;
	    if (c <= window->enonS) blend->nlimS++;
	  }
	}

	/* Allocate memory                                           */
	blend->c = i_alloc_1d(blend->ne1b + 1);
	blend->p = i_alloc_1d(blend->ne1b + 1);
	blend->m = i_alloc_1d(blend->ne1b + 1);
	blend->cs = i_alloc_1d(blend->ne1b + 1);
	blend->ce = i_alloc_1d(blend->ne1b + 1);
	blend->os = i_alloc_1d(blend->ne1b + 1);
	blend->oe = i_alloc_1d(blend->ne1b + 1);
	blend->lim = i_alloc_1d(blend->nlim + 1);

	/* Get the limit (end) locations of the blend zone           */
	i = 1;
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      (window->s2j[c] == master->e2bs[tn] ||
	       window->s2j[c] == master->e2be[tn])) {
	    blend->lim[i] = c;
	    i++;
	  }
	}
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (c == window->zp1[c]) continue;
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      (window->s2j[c] == master->e2bs[tn] ||
	       window->s2j[c] == master->e2be[tn])) {
	    blend->lim[i] = c;
	    i++;
	  }
	}

	/* Get the surface locations in the zone                     */
	i = 1;
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      window->s2j[c] >= master->e2bs[tn] &&
	      window->s2j[c] <= master->e2be[tn]) {
	    blend->c[i] = c;
	    blend->p[i] = window->xp1[c];
	    blend->m[i] = window->xm1[c];
	    c1 = c;
	    while (window->s2i[c1] !=  master->e1bs[tn]) {
	      c1 = window->xm1[c1];
	      if (c1 == window->xm1[c1]) {
		c1 = window->xp1[c1];
		break;
	      }
	    }
	    blend->cs[i] = c1;
	    blend->os[i] = window->xm1[c1];
	    c1 = c;
	    while (window->s2i[c1] !=  master->e1be[tn]) {
	      c1 = window->xp1[c1];
	      if (c1 == window->xp1[c1]) {
		c1 = window->xm1[c1];
		break;
	      }
	    }
	    blend->ce[i] = c1;
	    blend->oe[i] = window->xp1[c1];
	    i++;
	  }
	}

	/* Get the sub-surface locations                             */
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (c == window->zp1[c]) continue;
	  if (window->s2i[c] >= master->e1bs[tn] && 
	      window->s2i[c] <= master->e1be[tn] && 
	      window->s2j[c] >= master->e2bs[tn] &&
	      window->s2j[c] <= master->e2be[tn]) {
	    blend->c[i] = c;
	    blend->p[i] = window->xp1[c];
	    blend->m[i] = window->xm1[c];
	    c1 = c;
	    while (window->s2i[c1] !=  master->e1bs[tn]) {
	      c1 = window->xm1[c1];
	      if (c1 == window->xm1[c1]) {
		c1 = window->xp1[c1];
		break;
	      }
	    }
	    blend->cs[i] = c1;
	    blend->os[i] = window->xm1[c1];
	    c1 = c;
	    while (window->s2i[c1] !=  master->e1be[tn]) {
	      c1 = window->xp1[c1];
	      if (c1 == window->xp1[c1]) {
		c1 = window->xm1[c1];
		break;
	      }
	    }
	    blend->ce[i] = c1;
	    blend->oe[i] = window->xp1[c1];
	    i++;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the zones in e2 directions                                  */
  if (master->nbl2) {
    wincon->nbl2 = master->nbl2;
    wincon->ble2 = (blend_t **)malloc(sizeof(blend_t *) * master->nbl2);
    winsize = geom->nce1;
    for (tn = 0; tn < wincon->nbl2; tn++) {
      wincon->ble2[tn] = (blend_t *)malloc(sizeof(blend_t));
      blend_t *blend = wincon->ble2[tn];
      if (master->e2bs[tn] >= 0 && master->e2bs[tn] < geom->nce2 && 
	  master->e2be[tn] >= 0 && master->e2be[tn] < geom->nce2) {

	/* Reverse the indicies if required                          */
	if (master->e2bs[tn] > master->e2be[tn]) {
	  c = master->e2be[tn];
	  master->e2be[tn] = master->e2bs[tn];
	  master->e2bs[tn] = c;
	}
	if (master->e1bs[tn] > master->e1be[tn]) {
	  c = master->e1be[tn];
	  master->e1be[tn] = master->e1bs[tn];
	  master->e1bs[tn] = c;
	}

	/* Count the number of cells                                 */
	blend->ne2b = blend->ne2bS = 0;
	blend->nlim = blend->nlimS = 0;
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      window->s2i[c] >= master->e1bs[tn] &&
	      window->s2i[c] <= master->e1be[tn]) {
	    blend->ne2b++;
	    if (c <= window->enonS) blend->ne2bS++;
	  }
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      (window->s2i[c] == master->e1bs[tn] ||
	       window->s2i[c] == master->e1be[tn])) {
	    blend->nlim++;
	    if (c <= window->enonS) blend->nlimS++;
	  }
	}

	/* Allocate memory                                           */
	blend->c = i_alloc_1d(blend->ne2b + 1);
	blend->p = i_alloc_1d(blend->ne2b + 1);
	blend->m = i_alloc_1d(blend->ne2b + 1);
	blend->cs = i_alloc_1d(blend->ne2b + 1);
	blend->ce = i_alloc_1d(blend->ne2b + 1);
	blend->os = i_alloc_1d(blend->ne2b + 1);
	blend->oe = i_alloc_1d(blend->ne2b + 1);
	blend->lim = i_alloc_1d(blend->nlim + 1);

	/* Get the limit locations of the blend zone                 */
	i = 1;
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      (window->s2i[c] == master->e1bs[tn] ||
	       window->s2i[c] == master->e1be[tn])) {
	    blend->lim[i] = c;
	    i++;
	  }
	}
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (c == window->zp1[c]) continue;
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      (window->s2i[c] == master->e1bs[tn] ||
	       window->s2i[c] == master->e1be[tn])) {
	    blend->lim[i] = c;
	    i++;
	  }
	}

	/* Get the surface locations                                 */
	i = 1;
	for (cc = 1; cc <= window->b2_t; cc++) {
	  c = window->w2_t[cc];
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      window->s2i[c] >= master->e1bs[tn] &&
	      window->s2i[c] <= master->e1be[tn]) {
	    blend->c[i] = c;
	    blend->p[i] = window->yp1[c];
	    blend->m[i] = window->ym1[c];
	    c1 = c;
	    while (window->s2j[c1] !=  master->e2bs[tn]) {
	      c1 = window->ym1[c1];
	      if (c1 == window->ym1[c1]) {
		c1 = window->yp1[c1];
		break;
	      }
	    }
	    blend->cs[i] = c1;
	    blend->os[i] = window->ym1[c1];
	    c1 = c;
	    while (window->s2j[c1] !=  master->e2be[tn]) {
	      c1 = window->yp1[c1];
	      if (c1 == window->yp1[c1]) {
		c1 = window->ym1[c1];
		break;
	      }
	    }
	    blend->ce[i] = c1;
	    blend->oe[i] = window->yp1[c1];
	    i++;
	  }
	}

	/* Get the sub-surface locations                             */
	for (cc = 1; cc <= window->b3_t; cc++) {
	  c = window->w3_t[cc];
	  if (c == window->zp1[c]) continue;
	  if (window->s2j[c] >= master->e2bs[tn] && 
	      window->s2j[c] <= master->e2be[tn] && 
	      window->s2i[c] >= master->e1bs[tn] &&
	      window->s2i[c] <= master->e1be[tn]) {
	    blend->c[i] = c;
	    blend->p[i] = window->yp1[c];
	    blend->m[i] = window->ym1[c];
	    c1 = c;
	    while (window->s2j[c1] !=  master->e2bs[tn]) {
	      c1 = window->ym1[c1];
	      if (c1 == window->ym1[c1]) {
		c1 = window->yp1[c1];
		break;
	      }
	    }
	    blend->cs[i] = c1;
	    blend->os[i] = window->ym1[c1];
	    c1 = c;
	    while (window->s2j[c1] !=  master->e2be[tn]) {
	      c1 = window->yp1[c1];
	      if (c1 == window->yp1[c1]) {
		c1 = window->ym1[c1];
		break;
	      }
	    }
	    blend->ce[i] = c1;
	    blend->oe[i] = window->yp1[c1];
	    i++;
	  }
	}
      }
    }
  }
}

/* END set_blend_zones()                                             */
/*-------------------------------------------------------------------*/
