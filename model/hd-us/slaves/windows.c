/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/slaves/windows.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: windows.c 5943 2018-09-13 04:39:09Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* window_build()    : generates the window data structures          */
/* get_local_maps() : creates the mapping function in each window    */
/* local_map_build() : sets spatial maps                             */
/* reset_map()       : initialises auxiliary cells                   */
/* get_local_wsa()  : generates the local work sparse arrays         */
/* wsa_cells()      : fills the cells to process vectors             */
/* OBC_build()       : makes the window open boundary structures     */
/* surfbot_build()    : generates the surface and bottom maps        */
/* reorder_gl_map() : re-orders the global to local map              */
/* windat_fill()    : fills the window data with master data         */
/* window_alloc()     : allocates memory for window structures       */
/* OBC_alloc()     : allocates memory for OBC data structures        */
/* win_data_alloc()  : allocates memory for window data structures   */
/* master_alloc()  : allocates memory for the master data structure  */
/* win_data_build()  : makes window data structures for all windows  */
/* win_data_init()  : allocates memory for the window data structure */
/* win_data_clear() : clears memory for window data structures       */
/* get_timesteps()  : sets CFL time-step in each window              */
/* timeaux()        : finds cells common to windows for multi dt     */
/*-------------------------------------------------------------------*/

int *win_vector_build(geometry_t *geom, geometry_t *window, int wn,
		      int *gv, int mode);
int get_zcell(geometry_t *geom, open_bdrys_t *open, int cg, int zfe1, 
	      int zfe2, int mode);
int get_zcell_cy(geometry_t *geom, geometry_t *window, open_bdrys_t *open,
		 int c, int mode);
void get_process_exclude(parameters_t *params, geometry_t *geom, geometry_t *window);
void wave_setup(geometry_t *geom, geometry_t *window, int mode);
void nan_check(geometry_t **window, window_t **windat, 
	       win_priv_t **wincon, int nwindows);
int oedge(int npe, int n);

void get_local_obc_a(geometry_t *geom, int *vec, int nvec, int nvec2D, geometry_t **window, int nwindows); 		     
void set_reef_frac(master_t *master, geometry_t **window, window_t **windat, 
		     win_priv_t **wincon);
void aux_maps(geometry_t *geom, geometry_t *window, int cc);
void alloc_geom_us(geometry_t *geom, unsigned long mask);
void get_inner_exmap(geometry_t *geom, geometry_t *window);
int eiw(geometry_t *geom, int e, int *mode);
int eic(geometry_t *geom, int e, int wn, int *mode);
int jow(geometry_t *geom, int c, int j);
int joc(geometry_t *geom, int c, int cn);
void local_map_build_v2c(geometry_t *geom, int c, int cc, int wn,
			 int **nmap, int *wsa, int *ac, int nsaux);

/* Codes for edges */
#define W_GHOST 2
#define W_SAME  4
#define W_DIFF  8

/* Codes for vertices */
#define V_W  1
#define V_B  2
#define V_G  4

/* Codes for edges */
#define E_W  1
#define E_B  2
#define E_G  4

/*-------------------------------------------------------------------*/
/* Routine to check if an integer is a member of an integer array    */
/*-------------------------------------------------------------------*/
int ANY(int var,                /* Integer variable                  */
        int array[],            /* Integer array                     */
        int ns)
{                               /* Size of array                     */
  int nn;
  
  for (nn = 1; nn <= ns; nn++)
    if (var == array[nn])
      return (nn);
  return (0);
}
  
/*
 * This function is identical to the above ANY except it is zero-based
 * eg, called from write_run_setup
 * NOTE: This returns 1 or 0 depending if there was a match or not -
 *       which is different to above which returns the actual non-zero
 *       index
 */
int ANY0(
	 int var,     /* Integer variable */
	 int array[], /* Integer array    */
	 int ns       /* Size of array    */
	 )
{         
  int nn;
  
  for (nn = 0; nn < ns; nn++)
    if (var == array[nn])
      return (1);
  return (0);
}
  
/*
 * Same routine as above, but executes in parallel
 */
#if defined(HAVE_OMP)
int ANY_OMP(int var,      /* Integer variable */
	    int array[],  /* Integer array    */
	    int ns)       /* Size of array    */
{
  int nn;
  int done = 0; /* flag */
  int ret  = 0; /* This is the return value */

#pragma omp parallel for private(nn), shared(done,ret)
  for (nn = 1; nn <= ns; nn++) {
    if (!done) {
      if (var == array[nn]) {
#pragma omp critical
	{
	  ret  = nn;
	  done = 1;
	}
      }
    }
  }	
  
  return (ret);
}
#endif

/* ANYf stands for ANY function */
#ifdef HAVE_OMP
#define ANYf(x,y,z) ANY_OMP(x,y,z)
#else
#define ANYf(x,y,z) ANY(x,y,z)
#endif

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Window to decompose the global domain into windows. This routine  */
/* will do this by dividing the valid surface cells into the number  */
/* of required windows, and for each set of surface sets the window  */
/* is filled with sparse locations until the bottom is reached, i.e. */
/* each window consists of a number of water columns.                */
/*-------------------------------------------------------------------*/
void window_build(geometry_t *geom,     /* Global geometry           */
		  parameters_t *params  /* Parameter structure       */
		  )
{
  int n, j;                     /* Window counters                   */
  int c, cc, ee, e, c1;         /* Sparse location counters          */
  int vv, v, vs;                /* Vertex counters                   */
  int cs, es;                   /* Surface counters                  */
  int *wsa;                     /* Window sparse array               */
  int **ws2;                    /* Window sparse array               */
  int *wsz, *wsz2D;             /* Size of 3D/2D window arrays       */
  int nwindows;                 /* Number of windows to make         */
  int readwin = 0;              /* =1 : read window map from file    */

  nwindows = geom->nwindows;
  if (DEBUG("init_w"))
    dlog("init_w", "Start making %d windows\n", nwindows);
  if (params->runmode & MANUAL && nwindows > 1 && strlen(params->win_file)) 
    readwin = 1;

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  geom->fm = (global_map_t *)malloc(sizeof(global_map_t) * geom->szm);
  geom->g2we = i_alloc_2d(geom->sze, nwindows + 1);
  geom->g2wv = i_alloc_2d(geom->szv, nwindows + 1);
  window = (geometry_t **)malloc(sizeof(geometry_t *) * (nwindows + 1));

  wsa = i_alloc_1d(geom->szc);
  ws2 = i_alloc_2d(geom->szcS, nwindows + 1);
  wsz = i_alloc_1d(nwindows + 1);
  wsz2D = i_alloc_1d(nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Initialise the global to local map                              */
  for (c = 0; c < geom->szc; c++) {
    geom->fm[c].wn = 0;
    geom->fm[c].sc = 0;
    geom->fm[c].ac = 0;
  }

  for (ee = 0; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    geom->fm[e].ec = 0;
    geom->fm[e].we = 0;
    for (n = 1; n < nwindows; n++)
      geom->g2we[n][e] = 0;
  }
  for (vv = 0; vv <= geom->n3_e2; vv++) {
    v = geom->w3_e2[vv];
    geom->fm[v].vc = 0;
    for (n = 1; n < nwindows; n++)
      geom->g2wv[n][v] = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Single window : in this case window arrays point to the global  */
  /* arrays.                                                         */
  if (nwindows == 1) {
    window[1] = window_alloc();
    /* Get the sparse array sizes                                    */
    window[1]->enon = geom->sgnum;
    window[1]->sgsiz = window[1]->szc = geom->sgnum + 1;
    window[1]->ewet = geom->a3_t;
    window[1]->snon = geom->sgnum - geom->nbpt + 1;
    window[1]->enonS = geom->sgnumS;
    window[1]->sgsizS = window[1]->szcS = geom->sgnumS + 1;
    window[1]->ewetS = geom->a2_t;
    window[1]->snonS = geom->sgnumS - geom->nbptS + 1;
    window[1]->nz = geom->nz;
    window[1]->layers = geom->layers;
    window[1]->sednz = geom->sednz;
    window[1]->nwindows = nwindows;
    window[1]->wn = 1;
    /* Ghost sediment cells are only used in the semi-Lagrange       */
    /* scheme, which currently only runs with 1 window.              */
    window[1]->ngsed = geom->ngsed;
    window[1]->gsed_t = geom->gsed_t;
    window[1]->ised_t = geom->ised_t;
    window[1]->npem = geom->npem;
    window[1]->nvem = geom->nvem;
    window[1]->nvcm = geom->nvcm;
    window[1]->neem = geom->neem;
    window[1]->szc = geom->szc;
    window[1]->szcS = geom->szcS;
    window[1]->sze = geom->sze;
    window[1]->szeS = geom->szeS;
    window[1]->szv = geom->szv;
    window[1]->szvS = geom->szvS;
    window[1]->us_type = geom->us_type;
    window[1]->map = geom->map;

    /* Allocate memory for the window structure                      */
    alloc_geom_us(window[1], WINDOW_A);
    /* Point the window arrays to the global arrays                  */
    point_geom_us(window[1], geom);

    /* Define the diagnol maps from the directional maps             */
    /* Define the local - global maps                                */
    for (c = 1; c <= geom->sgnum; c++) {
      window[1]->wsa[c] = c;
      geom->fm[c].sc = c;
    }

    /* Define the fm.wn map in window #1 for wet cells only          */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      geom->fm[c].wn = 1;
    }
    for (e = 1; e < geom->sze; e++) {
      window[1]->wse[e] = e;
      geom->fm[e].ec = e;
      geom->fm[e].we = 1;
      geom->g2we[1][e] = e;
    }
    for (v = 1; v < geom->szv; v++) {
      window[1]->wsv[v] = v;
      geom->fm[v].vc = v;
      geom->g2wv[1][v] = v;
    }
    if (DEBUG("init_w"))
      dlog("init_w", "  Local maps for window 1 created OK\n");
  }
  /*-----------------------------------------------------------------*/
  /* Multiple windows : here the window sizes, spatial maps and      */
  /* processsing vectors must be defined for each window.            */
  else {

    /* Partition the surface layer into windows                      */
    if (readwin) {
      read_windows_us(geom, window, params->win_file);
    } else {
      if (params->win_type & STRIPE_E2)
	window_cells_linear_e2(geom, nwindows, ws2, wsz2D);
      else if (params->win_type & BLOCK_E1)
	window_cells_block_e1(geom, nwindows, ws2, wsz2D, params->win_block);
      else if (params->win_type & BLOCK_E2)
	window_cells_block_e2(geom, nwindows, ws2, wsz2D, params->win_block);
      else if (params->win_type & STRIPE_E1)
	window_cells_linear_e1(geom, nwindows, ws2, wsz2D);
      else
	window_cells_grouped(geom, nwindows, ws2, wsz2D);

      /* Get the global to local maps                                */
      get_gl_maps(geom, nwindows, ws2, wsz2D);

      if (DEBUG("init_w"))
	dlog("init_w", "  Global local map created OK\n");

      /*-------------------------------------------------------------*/
      /* Get the 3D sparse locations in each window based on the     */
      /* surface partitioning defined above.                         */
      for (n = 1; n <= nwindows; n++) {
	int cellf = 0;
	window[n] = window_alloc();
	window[n]->wn = n;
	window[n]->npem = geom->npem;
	window[n]->nvem = geom->nvem;
	window[n]->nvcm = geom->nvcm;
	window[n]->neem = geom->neem;
	window[n]->nwindows = nwindows;
	window[n]->us_type = geom->us_type;
	window[n]->nz = geom->nz;
	window[n]->layers = d_alloc_1d(window[n]->nz + 1);
	memcpy(window[n]->layers, geom->layers, (window[n]->nz + 1) * sizeof(double));
	window[n]->sednz = geom->sednz;

	/* Get the 3D local cells in window n                        */
	get_window_cells_h(geom, n, wsa, wsz, ws2[n], wsz2D[n]);

	/* Set the local wet array and local maps for window n       */
	get_local_maps(geom, window[n], n, wsa, wsz[n], ws2[n], wsz2D[n], cellf);

	if (DEBUG("init_w"))
	  dlog("init_w", "  Local maps for window %d created OK\n", n);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Create the open boundary structures for each window             */
  OBC_build(geom->open, geom, window, nwindows);

  if (!readwin) {
    /*---------------------------------------------------------------*/
    /* Set up the local processing vectors arrays for window wn      */
    get_local_wsc(geom, geom->w3_t, geom->v3_t, geom->v2_t, window, nwindows);
    get_local_wse(geom, geom->w3_e1, geom->v3_e1, geom->v2_e1, window, nwindows,
		  params->smagorinsky);
    get_local_wsv(geom, geom->w3_e2, geom->v3_e2, geom->v2_e2, window, nwindows);

    if (DEBUG("init_w"))
      dlog("init_w", "  Local work arrays created OK\n\n");
  } else {
    if (!(params->compatible & V1957))
      get_local_obc_a(geom, geom->w3_t, geom->v3_t, geom->v2_t, window, nwindows);
  }

  /*-----------------------------------------------------------------*/
  /* Set the free surface and bottom maps                            */
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,cc,cs,c1,e,ee,es,v,vv,vs,j)
#endif
  for (n = 1; n <= nwindows; n++) {
    int check = 0;

    if (!readwin) {
      window[n]->sur_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->nsur_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->bot_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->sur_e1 = i_alloc_1d(window[n]->szeS);
      window[n]->bot_e1 = i_alloc_1d(window[n]->szeS);
      window[n]->sur_e2 = i_alloc_1d(window[n]->szvS);
      window[n]->bot_e2 = i_alloc_1d(window[n]->szvS);
      window[n]->m2d = i_alloc_1d(window[n]->szc);
      window[n]->m2de = i_alloc_1d(window[n]->sze);
      window[n]->m2dv = i_alloc_1d(window[n]->szv);
      window[n]->s2i = i_alloc_1d(window[n]->szc);
      window[n]->s2j = i_alloc_1d(window[n]->szc);
      window[n]->s2k = i_alloc_1d(window[n]->szc);
      window[n]->e2ijk = i_alloc_1d(window[n]->sze);
      window[n]->v2ijk = i_alloc_1d(window[n]->szv);

      surfbot_build(geom, window[n], n, window[n]->nsur_t, window[n]->bot_t,
		    geom->bot_t, geom->v2_t, window[n]->w2_t,
		    window[n]->v2_t, window[n]->a2_t);
      surfbot_builde(geom, window[n], n, window[n]->sur_e1, window[n]->bot_e1,
		     geom->bot_e1, geom->v2_e1, window[n]->w2_e1,
		     window[n]->v2_e1, window[n]->x2_e1);

      /*-------------------------------------------------------------*/
      /* Get the 3D - 2D map (including the sediment). Note: the     */
      /* bottom ghost cell maps vertically to the interior sediment  */
      /* cell, hence we do not set m2d maps for ghost sediments.     */
      /* Cell centres                                                */
      for (cc = 1; cc <= window[n]->enon; cc++)
	window[n]->m2d[cc] = 0;
      for (cc = 1; cc <= window[n]->n2_t; cc++) {
	c = cs = window[n]->w2_t[cc];
	while (window[n]->zm1[c] && c != window[n]->zm1[c]) {
	  window[n]->m2d[c] = cs;
	  c = window[n]->zm1[c];
	}
	if (cc <= window[n]->ewetS) window[n]->m2d[c] = cs;
      }
      /* Get the 3D - 2D map for sub-surface ghost cells (i.e. ghost */
      /* cells where the cell above is wet).                         */
      for (cc = 1; cc <= window[n]->enon; cc++) {
	c = cc;
	c1 = window[n]->zp1[cc];
	while (window[n]->m2d[cc] == 0 && c != c1) {
	  window[n]->m2d[cc] = window[n]->m2d[c1];
	  c = c1;
	  c1 = window[n]->zp1[c1];
	}
      }

      /* Cell edges                                                  */
      for (ee = 1; ee <= window[n]->n3_e1; ee++) {
	e = window[n]->w3_e1[ee];
	window[n]->m2de[e] = 0;
      }
      for (es = 1; es < window[n]->szeS; es++) {
	e = es;
	while (window[n]->zm1e[e] && e != window[n]->zm1e[e]) {
	  window[n]->m2de[e] = es;
	  e = window[n]->zm1e[e];
	}
	window[n]->m2de[e] = es;
      }
      for (ee = 1; ee <= window[n]->n3_e1; ee++) {
	e = window[n]->w3_e1[ee];
	if (window[n]->m2de[e] == 0) {
	  es = geom->m2de[window[n]->wse[e]];
	  window[n]->m2de[e] = geom->g2we[n][es];
	}
      }

      /* Reset the vertical maps                                     */
      for (ee = 1; ee <= window[n]->n2_e1; ee++) {
	int eb = window[n]->bot_e1[ee];
	e = window[n]->w2_e1[ee];
	while (e != window[n]->zm1e[e])
	  e = window[n]->zm1e[e];
	window[n]->zm1e[eb] = e;
	window[n]->zp1e[e] = eb;
      }

      /* Cell vertexes                                               */
      for (vv = 1; vv <= window[n]->n3_e2; vv++) {
	v = window[n]->w3_e2[vv];
	window[n]->m2dv[vv] = 0;
      }

      for (vv = 1; vv <= window[n]->n2_e2; vv++) {
	v = vs = window[n]->w2_e2[vv];
	while (window[n]->zm1v[v] && v != window[n]->zm1v[v]) {
	  window[n]->m2dv[v] = vs;
	  v = window[n]->zm1v[v];
	}
	window[n]->m2dv[v] = vs;
      }

      if (check) {
	for (vv = 1; vv <= window[n]->n3_e2; vv++) {
	  v = window[n]->w3_e2[vv];
	  if (v) {
	    vs = window[n]->m2dv[v];
	    c = window[n]->wse[v];
	    if (vs == 0) printf("Can't find 2D vertex map at vv=%d, v=%d(%d)[%f %f] wn=%d\n",vv, v, c,
				 geom->gridx[geom->m2dv[c]], geom->gridy[geom->m2dv[c]], n);
	  } else
	    printf("zero at %d\n",vv);
	}
      }

      /* Check for un-mapped cells, and define using lateral maps    */
      for (cc = 1; cc <= window[n]->enon; cc++) {
	cs = window[n]->m2d[cc];
	if(cs == 0) {
	  for (j = 1; j <= window[n]->npem; j++) {
	    if ((c = window[n]->c2c[j][window[n]->m2d[window[n]->c2c[oedge(window[n]->npem,j)][cc]]])) {
	      window[n]->m2d[cc] = c;
	      break;
	    }
	  }
	  if(window[n]->m2d[cc] == 0)
	    emstag(LDEBUG,"hd:windows:window_build",
		   "Can't find surface map for window %d cell %d\n", n, cc);
	}
      }

      /*-------------------------------------------------------------*/
      /* Get the 3D - 1D map                                         */
      memset(window[n]->s2k, 0, window[n]->szc * sizeof(int));
      for (cc = 1; cc <= window[n]->n3_t; cc++) {
	c = window[n]->w3_t[cc];
	c1 = window[n]->wsa[c];
	window[n]->s2i[c] = geom->s2i[c1];
	window[n]->s2j[c] = geom->s2j[c1];
	window[n]->s2k[c] = geom->s2k[c1];
      }

      for (ee = 1; ee <= window[n]->n3_e1; ee++) {
	e = window[n]->w3_e1[ee];
	window[n]->e2ijk[e] = NOTVALID;
      }

      /* Edge maps to cell centre to retrieve (i,j,k)                */
      for (ee = 1; ee <= window[n]->n3_e1; ee++) {
	e = window[n]->w3_e1[ee];
	c = window[n]->e2c[e][0];
	if (geom->wgst[window[n]->wsa[c]])
	  c = window[n]->e2c[e][1];
	window[n]->e2ijk[e] = c;
      }
      for (vv = 1; vv <= window[n]->b3_e2; vv++) {
	v = window[n]->w3_e2[vv];
	window[n]->v2ijk[v] = window[n]->v2c[v][1];
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the coordinate - index map (currently only used with      */
    /* sediments).                                                   */
    c1 = 0;
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (c > c1)
	c1 = c;
    }
    window[n]->c2cc = i_alloc_1d(c1 + 1);
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      window[n]->c2cc[c] = cc;
    }

    if (!readwin) {
      /*-------------------------------------------------------------*/
      /* Get the transfer vectors                                    */
      build_transfer_maps(geom, window, n, nwindows);

      /*-------------------------------------------------------------*/
      /* Get the ghost cell vectors                                  */
      get_local_ghost(geom, window[n], n);
    }

    /*---------------------------------------------------------------*/
    /* Build any local vectors required                              */
    if (window[n]->thetau1 == NULL)
      window[n]->thetau1 = d_alloc_1d(window[n]->szeS);
    if (window[n]->thetau2 == NULL)
      window[n]->thetau2 = d_alloc_1d(window[n]->szeS);
    for (ee = 1; ee <= window[n]->n2_e1; ee++) {
      e = window[n]->w2_e1[ee];
      es = window[n]->wse[e];
      window[n]->thetau1[e] = geom->thetau1[es];
      window[n]->thetau2[e] = geom->thetau2[es];
    }

    if(geom->sm_e1)
      window[n]->sm_e1 = win_vector_build(geom, window[n], n, geom->sm_e1, 1);
    if(geom->sm_e2)
      window[n]->sm_e2 = win_vector_build(geom, window[n], n, geom->sm_e2, 2);

    get_inner_exmap(geom, window[n]);
    get_process_exclude(params, geom, window[n]);
    set_mask(window[n]);
  }

  if (DEBUG("init_w"))
    dlog("init_w", "  Surface and bottom arrays created OK\n");

  /*-----------------------------------------------------------------*/
  /* Free memory                                                     */
  i_free_1d(wsa);
  i_free_2d(ws2);
  i_free_1d(wsz);
  i_free_1d(wsz2D);

  /*-----------------------------------------------------------------*/
  /* Print window information                                        */
  for (n = 1; n <= nwindows; n++) {
    if (DEBUG("init_w")) {
      dlog("init_w", "\n");
      dlog("init_w", "Window #%d\n", n);
      dlog("init_w", "  number of 3D wet cells = %d (%4.1f %%)\n",
           window[n]->ewet, (double)(100 * window[n]->ewet / geom->a3_t));
      dlog("init_w", "  number of 2D wet cells = %d (%4.1f %%)\n",
           window[n]->ewetS,
           (double)(100 * window[n]->ewetS / geom->a2_t));
      c = window[n]->enon - window[n]->snon + 1;
      dlog("init_w", "  number of 3D auxiliary cells = %d (%d %%)\n", c,
           100 * c / window[n]->enon);
      c = window[n]->enonS - window[n]->snonS + 1;
      dlog("init_w", "  number of 2D auxiliary cells = %d (%d %%)\n", c,
           100 * c / window[n]->enonS);
      dlog("init_w", "  number of 3D tracer wet cells to process = %d\n",
           window[n]->v3_t);
      dlog("init_w", "  number of 2D tracer wet cells to process = %d\n",
           window[n]->v2_t);
      dlog("init_w", "  number of 3D u1 wet cells to process = %d\n",
           window[n]->v3_e1);
      dlog("init_w", "  number of 2D u1 wet cells to process = %d\n",
           window[n]->v2_e1);
      dlog("init_w", "  number of 3D u2 wet cells to process = %d\n",
           window[n]->v3_e2);
      dlog("init_w", "  number of 2D u2 wet cells to process = %d\n",
           window[n]->v2_e2);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Free memory not required in future computations                 */
  for (n = 1; n <= nwindows; n++) {
    /* This map is of size geom->szc in each window with only        */
    /* window[n]->enon locations used.                               */
    if (nwindows > 1)
      free((global_map_t *)window[n]->fm);
  }
  i_free_2d(geom->g2we);
  i_free_2d(geom->g2wv);
}

/* END window_build()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create an ascii file of window partitions              */
/*-------------------------------------------------------------------*/
void write_windows(geometry_t *geom,      /* Sparse global geometery */
                   unsigned long ***flag  /* Flags structure         */
  )
{
  int c;                        /* Sparse counters                   */
  int n;                        /* Window counter                    */
  int i, j, k;
  FILE *op, *fopen();
  short **aa;

  if (!windows_log)
    return;
  aa = s_alloc_2d(geom->nfe1, geom->nfe2);
  op = fopen(window_geom_logfile, "w");
  if (op == NULL)
    hd_quit("write_windows: Can't open file %s\n", window_geom_logfile);
  k = geom->nz - 1;
  for (i = 0; i < geom->nfe1; i++)
    for (j = 0; j < geom->nfe2; j++) {
      aa[j][i] = 0;
      if (!(flag[k][j][i] & (SOLID | OUTSIDE))) {
        c = geom->map[k][j][i];
        n = geom->fm[c].wn;
        aa[j][i] = n;
      }
    }
  for (j = geom->nfe2 - 1; j >= 0; j--) {
    for (i = 0; i < geom->nfe1; i++) {
      fprintf(op, "%d ", aa[j][i]);
    }
    fprintf(op, "\n");
  }
  fclose(op);

  s_free_2d(aa);
}

/* END write_windows()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Divides the domain into nwindows regions where each region        */
/* contains adjacent cells.                                          */
/*-------------------------------------------------------------------*/
void window_cells_grouped(geometry_t *geom,   /* Global geometery    */
			  int nwindows,   /* Number of windows       */
			  int **ws2,      /* 2D wet cells in window  */
			  int *wsizeS     /* Number of 2D wet cells  */
		  )
{
  int cc, c, ci, cn, j, n;
  int found, ni, *filla;
  int wi, wn = 1;
  int ns = geom->szcS;
  int *mask, *wm;
  int checkf = 1;
  int verbose = 0;
  int v2_t = geom->v2_t;

  /* Initialise                                                      */
  wm = i_alloc_1d(geom->szcS);
  memset(wm, 0, geom->szcS * sizeof(int));

  /* Set the open boundary mask                                      */
  mask = i_alloc_1d(geom->szcS);
  memset(mask, 0, geom->szcS * sizeof(int));
  for (cc = geom->v2_t + 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    mask[c] = 1;
  }
  /* Sponge zones are allocated to a single window                   */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->sponge_zone_h) {
      v2_t -= (open->nspc - open->no2_t);
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	mask[c] = 1;
      }
    }
  }

  /* Get the number of cells in this window                          */
  c = (int)ceil(v2_t / nwindows);
  for (n = 1; n < nwindows; n++)
    wsizeS[n] = c;
  wsizeS[nwindows] = v2_t - (nwindows - 1) * c;

  /* If window sizes were specified, use these                       */
  if (geom->win_size) {
    cc = 0;
    for (n = 1; n < nwindows; n++) {
      wsizeS[n] = (int)ceil(v2_t * geom->win_size[n]);
      geom->win_size[n] = (double)wsizeS[n] / (double)v2_t;
      cc += wsizeS[n];
    }
    wsizeS[nwindows] = v2_t - cc;
    geom->win_size[nwindows] =
      (double)wsizeS[nwindows] / (double)v2_t;
  }

  /* Set the mask fanning out from the interior coordinate           */
  ci = 1;
  for (ci = 1; ci < geom->szcS; ci++) {
    if (!mask[ci]) break;
  }
  filla = i_alloc_1d(ns);
  memset(filla, 0, ns * sizeof(int));
  found = 1;
  ni = wi = 1;
  wm[ci] = wn;
  ws2[wn][wi++] = ci;
  ni++;
  filla[ci] = mask[ci] = 1;
  while (found) {
    found = 0;
    for (n = 1; n <= wn; n++) {
      for (cc = 1; cc <= ni; cc++) {
	c = ws2[n][cc];
	if (filla[c] == 2) continue;
	for (j = 1; j <= geom->npe[c]; j++) {
	  cn = geom->c2c[j][c];
	  if (geom->wgst[cn] || mask[cn]) continue;
	  if (!filla[cn]) {
	    filla[cn] = 1;
	    filla[c] = 2;
	    if (verbose) printf("wn=%d wi=%d c=%d[%d %d]\n",wn,wi,cn,geom->s2i[cn],geom->s2j[cn]);
	    wm[cn] = wn;
	    ws2[wn][wi++] = cn;
	    /*mask[cn] = 1;*/
	    ni++;
	    if (wi > wsizeS[wn] && wn < nwindows) {
	      wi = 1;
	      wn++;
	    }
	    found = 1;
	  }
	}
      }
    }
    if (found && ni > geom->szcS) 
      hd_warn("window_cells_grouped(): too many cells processed, %d out of %d\n", 
	      ni, geom->szcS);
  }

  /* Assign any 'holes' to the last window                           */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    if (!filla[c] && !mask[c]) {
      ws2[wn][wi++] = c;
      wm[c] = wn;
    }
  }

  /* Assign sponge zones to a single window                          */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (open->sponge_zone_h) {
      found = 0;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= geom->npe[c]; j++) {  
	  cn = geom->c2c[j][c];
	  if ((wn = wm[cn])) {
	    if(wn > 0) {
	      found = 1;
	      break;
	    }
	  }
	}
	if (found) break;
      }
      if (!found) hd_quit("window_cells_grouped: Can't find window for sponge zone OBC %s\n", open->name);
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	wsizeS[wn]++;
	ws2[wn][wsizeS[wn]] = c;
	wm[c] = wn;
	if (geom->win_size)
	  geom->win_size[wn] = (double)wsizeS[wn] / (double)geom->b2_t;
      }
    }
  }

  /* Assign the obc cells to the window corresponding to their       */
  /* interior wet neighbours.                                        */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    if (!open->sponge_zone_h) {
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	ci = open->oi1_t[cc];
	wn = wm[ci];
	wsizeS[wn]++;
	ws2[wn][wsizeS[wn]] = c;
	if (verbose) printf("OBC wn=%d wi=%d c=%d[%d %d]\n",wn,wsizeS[wn],c,geom->s2i[c],geom->s2j[c]);
	wm[c] = wn;
	if (geom->win_size)
	  geom->win_size[wn] = (double)wsizeS[wn] / (double)geom->b2_t;
      }
    }
  }

  /* Check all cells are accounted for                               */
  if (checkf) {
    int wc[nwindows + 1];
    int mw[nwindows + 1];
    int nb;
    for (n = 1; n <= nwindows; n++) wc[n] = 0;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      found = 0;
      for (n = 1; n <= nwindows; n++) {
	for (j = 1; j <= wsizeS[n]; j++) 
	  if (c == ws2[n][j]) {
	    found = 1;
	    wc[n]++;
	    break;
	  }
	if (found) break;
      }
      if (!found) {
      /* Add any cells not found to neighboring windows
	for (n = 1; n <= nwindows; n++) mw[n] = 0;
	for (j = 1; j <= geom->npe[c]; j++) {
	  cn = geom->c2c[j][c];
	  if ((wn = wm[cn])) mw[wn] += 1;
	}
	wn = 1;
	wi = mw[wn];
	for (n = 1; n <= nwindows; n++) {
	  if (mw[n] > wi) {
	    wi = mw[n];
	    wn = n;
	  }
	}
	wsizeS[wn]++;
	ws2[wn][wsizeS[wn]] = c;
	wm[c] = wn;
	if (geom->win_size)
	  geom->win_size[wn] = (double)wsizeS[wn] / (double)geom->b2_t;
	hd_warn("window_cells_grouped: Cell %d[%d %d] not found: allocated to window%d (%f %f)\n", c, geom->s2i[c], geom->s2j[c], wn, geom->cellx[c], geom->celly[c]);
	*/
	hd_warn("window_cells_grouped: Can't find window partition for cell %d[%d %d] %f %f\n", c, geom->s2i[c], geom->s2j[c], geom->cellx[c], geom->celly[c]);
      }
    }
    for (n = 1; n <= nwindows; n++) {
      if (wc[n] != wsizeS[n])
	hd_quit("window_cells_grouped: Incorrect number of cells for window %d: %d != %d\n", n, wc[n], wsizeS[n]);
    }
  }

  i_free_1d(filla);
  i_free_1d(mask);
  i_free_1d(wm);
}

/* END window_cells_grouped()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to partition the surface layer into windows if the domain */
/* is divided into nwindows consecutive rows of cells as presented   */
/* in the work array w2_t. Any open boundary cells are assigned the  */
/* window number corresponding to their interior wet neighbour.      */
/*-------------------------------------------------------------------*/
void window_cells_linear_e1(geometry_t *geom, /* Global geometery    */
			    int nwindows, /* Number of windows       */
			    int **ws2,    /* 2D wet cells in window  */
			    int *wsizeS   /* Number of 2D wet cells  */
  )
{
  int c, cc, n;
  int mode = 0;
  int sm = 0;
  int *wm;
  int verbose = 0;

  /* Get the number of cells in this window                          */
  c = (int)ceil(geom->v2_t / nwindows);
  for (n = 1; n < nwindows; n++)
    wsizeS[n] = c;
  wsizeS[nwindows] = geom->v2_t - (nwindows - 1) * c;

  /* If window sizes were specified, use these                       */
  if (geom->win_size) {
    cc = 0;
    for (n = 1; n < nwindows; n++) {
      wsizeS[n] = (int)ceil(geom->v2_t * geom->win_size[n]);
      geom->win_size[n] = (double)wsizeS[n] / (double)geom->v2_t;
      cc += wsizeS[n];
    }
    wsizeS[nwindows] = geom->v2_t - cc;
    geom->win_size[nwindows] =
      (double)wsizeS[nwindows] / (double)geom->v2_t;
  }

  /* If window sizes are explicitly specified, use these             */
  for (n = 1; n <= nwindows; n++) {
    if (geom->nwn[n - 1]) {
      wsizeS[n] = geom->nwn[n - 1];
      mode = 1;
    }
    sm += wsizeS[n];
  }

  /* Save the surface cell locations in ws2                          */
  if (mode) {
    int cl;
    int *mask = i_alloc_1d(geom->szcS);
    memset(mask, 0, (geom->szcS) * sizeof(int));
    for (n = 1; n < nwindows; n++) {
      for (cc = 1; cc <= wsizeS[n]; cc++) {
	int i = geom->wnx[n - 1][cc - 1];
	int j = geom->wny[n - 1][cc - 1];
	int k = geom->nz - 1;
	c = geom->map[k][j][i];
	ws2[n][cc] = c;
	mask[c] = 1;
      }
    }
    n = nwindows; cc = 1;
    for (c = 1; c <= geom->a2_t; c++) {
      cl = geom->wsa[c];
      if (!mask[cl]) {
	ws2[n][cc] = cl;
	cc++;
	mask[cl] = 1;
      }
    }
    wsizeS[n] = cc - 1;
    i_free_1d(mask);
  } else {
    int nn, ci;
    wm = i_alloc_1d(geom->szcS);
    /* Assign the wet cells in the work array to the windows.        */
    n = 1;
    nn = 1;
    for (cc = 1; cc <= geom->v2_t; cc++) {
      c = geom->w2_t[cc];
      if (nn > wsizeS[n]) {
	n++;
	nn = 1;
      }
      ws2[n][nn++] = c;
      wm[c] = n;
    }
    /* Assign the obc cells to the window corresponding to their     */
    /* interior wet neighbours.                                      */
    for (nn = 0; nn < geom->nobc; nn++) {
      open_bdrys_t *open = geom->open[nn];
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	ci = open->oi1_t[cc];
	n = wm[ci];
	wsizeS[n]++;
	ws2[n][wsizeS[n]] = c;
	wm[c] = n;
	if (geom->win_size)
	  geom->win_size[n] = (double)wsizeS[n] / (double)geom->b2_t;
	if (verbose)
	  printf("%s wn=%d %d(%d %d) %d\n",open->name, n, c, geom->s2i[c], geom->s2j[c], wsizeS[n]);
      }
    }
    i_free_1d(wm);
  }
}

/* END window_cells_linear_e1()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to partition the surface layer into windows if the domain */
/* is divided into nwindows consecutive rows of cells. Currently the */
/* windows are partitioned via a loop through the grid in the k,j,i  */
/* direction and consecutively assigns cells to windows.             */
/*-------------------------------------------------------------------*/
void window_cells_linear_e2(geometry_t *geom, /* Global geometery    */
			    int nwindows, /* Number of windows       */
			    int **ws2,   /*2D wet cells in window wn */
			    int *wsizeS  /* Number of 2D wet cells   */
  )
{
  int c, cc, n, m;
  int c1, i, j;
  int *mask, *wetc;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  wetc = i_alloc_1d(geom->enonS + 1);
  mask = i_alloc_1d(geom->enonS + 1);

  /*-----------------------------------------------------------------*/
  /* Make a mask of wet cells                                        */
  memset(wetc, 0, (geom->enonS + 1) * sizeof(int));
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    wetc[c] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the number of cells in this window                          */
  c = (int)ceil(geom->a2_t / nwindows);
  for (n = 1; n < nwindows; n++)
    wsizeS[n] = c;
  wsizeS[nwindows] = geom->a2_t - (nwindows - 1) * c;

  /* If window sizes were specified, use these                       */
  if (geom->win_size) {
    cc = 0;
    for (n = 1; n < nwindows; n++) {
      wsizeS[n] = (int)ceil(geom->a2_t * geom->win_size[n]);
      geom->win_size[n] = (double)wsizeS[n] / (double)geom->a2_t;
      cc += wsizeS[n];
    }
    wsizeS[nwindows] = geom->a2_t - cc;
    geom->win_size[nwindows] =
      (double)wsizeS[nwindows] / (double)geom->a2_t;
  }

  /*-----------------------------------------------------------------*/
  /* Save the surface cell locations in ws2                          */
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  cc = n = 1;
  for (i = 0; i < geom->nfe1; i++) {
    for (j = 0; j < geom->nfe2; j++) {
      c = geom->map[geom->nz-1][j][i];
      if (mask[c]) continue;
      if (wetc[c]) {
	if (n > nwindows) 
	  hd_quit("Striping error in e2 direction : windows = %d.\n", n);
	ws2[n][cc] = c;
	mask[c] = 1;
	cc++;
        if (cc > wsizeS[n]) {
	  int bf = 1, nn, cb, cbc;
	  /* Don't allow window transition to occur on OBCs          */
	  for (nn = 0; nn < geom->nobc; nn++) {
	    open_bdrys_t *open = geom->open[nn];
	    if (open->type == U2BDRY && ANYf(c, open->obc_t, open->no2_t)) {
	      bf = 0;
	      wsizeS[n]+=1;
	      wsizeS[nwindows]-=1;
	    }
	    if (open->type == U2BDRY && ANYf(c, open->oi1_t, open->no2_t)) {
	      bf = 0;
	      wsizeS[n]+=2;
	      wsizeS[nwindows]-=2;
	    }
	  }
	  if (bf) {
	    n++;
	    cc = 1;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Collect any unallocated cells                                   */
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (n = 1; n <= nwindows; n++) {
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];
      mask[c] = 1;
    }
  }
  m = 0;
  i = wsizeS[nwindows]+1;
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    if (mask[c] == 0) {
      ws2[nwindows][i] = c;
      wsizeS[nwindows]++;
      i++;
      m++;
    }
  }
  if (DEBUG("init_w")) {
    dlog("init_w", "Found %d unallocated cells\n", m);
  }
  i_free_1d(wetc);
  i_free_1d(mask);
}

/* END window_cells_linear_e2()                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to partition the surface layer into windows if the domain */
/* is divided into blosks of size sqrt(nce1 * nce1). The offset, os, */
/* allows domains to be extended in the e1 direction, making the     */
/* blocks rectangular.                                               */
/*-------------------------------------------------------------------*/
void window_cells_block_e1(geometry_t *geom, /* Global geometery     */
			   int nwindows, /* Number of windows        */
			   int **ws2,    /*2D wet cells in window wn */
			   int *wsizeS,  /* Number of 2D wet cells   */
			   int os
			   )
{
  int c, cc, n, m;
  int c1, i, j, ii, jj;
  int *mask, *wetc;
  int wetf;
  double d1;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  wetc = i_alloc_1d(geom->enonS + 1);
  mask = i_alloc_1d(geom->enonS + 1);

  /*-----------------------------------------------------------------*/
  /* Make a mask of wet cells                                        */
  memset(wetc, 0, (geom->enonS + 1) * sizeof(int));
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    wetc[c] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the number of cells in this window                          */
  d1 = (double)(geom->nfe1 * geom->nfe2) / (double)(nwindows);
  d1 = (-(double)os+sqrt((double)(os * os) + 4.0 * d1))/2.0;
  c1 = (d1 > 2.0) ? (int)d1-1 : 1;

  /*-----------------------------------------------------------------*/
  /* Find the size of the block iteratively                          */
  cc = 0;
  while (cc <= geom->a2_t-c1*c1) {
    i = j = cc = 0;
    memset(mask, 0, (geom->enonS + 1) * sizeof(int));
    for (n = 1; n <= nwindows; n++) {
      wetf = 1;
      for (jj = j; jj < j + c1; jj++) {
	if (jj >= geom->nfe2) continue;
	for (ii = i; ii < i + c1 + os; ii++) {
	  if (ii >= geom->nfe1) continue;
	  c = geom->map[geom->nz-1][jj][ii];
	  if (mask[c]) continue;
	  if (wetc[c]) {
	    mask[c] = 1;
	    wetf = 0;
	    cc++;
	  }
	}
      }
      if (ii >= geom->nfe1 && jj >= geom->nfe2) break;
      i = ii;
      if (i >= geom->nfe1) {
	i = 0;
	j = jj;
      }
      if (jj > geom->nfe2) jj = j;
      if (wetf)n--;
    }
    c1++;
  }
  c1-=2;

  /*-----------------------------------------------------------------*/
  /* Save the surface cell locations in ws2                          */
  i = j = 0;
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (n = 1; n <= nwindows; n++) {
    wsizeS[n] = 0;
    cc = 1;
    wetf = 1;
    for (jj = j; jj < j + c1; jj++) {
      if (jj >= geom->nfe2) continue;
      for (ii = i; ii < i + c1 + os; ii++) {
	if (ii >= geom->nfe1) continue;
	c = geom->map[geom->nz-1][jj][ii];
	if (mask[c]) continue;
	if (wetc[c]) {
	  ws2[n][cc] = c;
	  mask[c] = 1;
	  wetf = 0;
	  wsizeS[n]++;
	  cc++;
	}
      }
    }
    if (ii >= geom->nfe1 && jj >= geom->nfe2) break;
    i = ii;
    if (i >= geom->nfe1) {
      i = 0;
      j = jj;
    }
    if (jj > geom->nfe2) jj = j;
    if (wetf) n--;
  }

  /*-----------------------------------------------------------------*/
  /* Collect any unallocated cells                                   */
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (n = 1; n <= nwindows; n++) {
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];
      mask[c] = 1;
    }
  }
  m = 0;
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    if (mask[c] == 0) {
      wsizeS[nwindows]++;
      ws2[nwindows][wsizeS[nwindows]] = c;
      m++;
    }
  }
  if (DEBUG("init_w")) {
    for (n = 1; n <= nwindows; n++)
      dlog("init_w", "Blocking in e1 direction: window%d = %d cells\n", n, wsizeS[n]);
  }
  i_free_1d(wetc);
  i_free_1d(mask);
}

/* END window_cells_block_e1()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to partition the surface layer into windows if the domain */
/* is divided into blosks of size sqrt(nce1 * nce1). The offset, os, */
/* allows domains to be extended in the e2 direction, making the     */
/* blocks rectangular.                                               */
/*-------------------------------------------------------------------*/
void window_cells_block_e2(geometry_t *geom, /* Global geometery     */
			   int nwindows, /* Number of windows        */
			   int **ws2,    /*2D wet cells in window wn */
			   int *wsizeS,  /* Number of 2D wet cells   */
			   int os
			   )
{
  int c, cc, n, m;
  int c1, i, j, ii, jj;
  int *mask, *wetc;
  int wetf;
  double d1;

  /*-----------------------------------------------------------------*/
  /* Initialise                                                      */
  wetc = i_alloc_1d(geom->enonS + 1);
  mask = i_alloc_1d(geom->enonS + 1);

  /*-----------------------------------------------------------------*/
  /* Make a mask of wet cells                                        */
  memset(wetc, 0, (geom->enonS + 1) * sizeof(int));
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    wetc[c] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Get the number of cells in this window                          */
  d1 = (double)(geom->nfe1 * geom->nfe2) / (double)(nwindows);
  d1 = (-(double)os+sqrt((double)(os * os) + 4.0 * d1))/2.0;
  c1 = (d1 > 2.0) ? (int)d1-1 : 1;

  /*-----------------------------------------------------------------*/
  /* Find the size of the block iteratively                          */
  cc = 0;
  while (cc <= geom->a2_t-c1*c1) {
    i = j = cc = 0;
    memset(mask, 0, (geom->enonS + 1) * sizeof(int));
    for (n = 1; n <= nwindows; n++) {
      wetf = 1;
      for (ii = i; ii < i + c1 + os; ii++) {
	if (ii >= geom->nfe1) continue;
	for (jj = j; jj < j + c1; jj++) {
	  if (jj >= geom->nfe2) continue;
	  c = geom->map[geom->nz-1][jj][ii];
	  if (mask[c]) continue;
	  if (wetc[c]) {
	    mask[c] = 1;
	    wetf = 0;
	    cc++;
	  }
	}
      }
      if (ii >= geom->nfe1 && jj >= geom->nfe2) break;
      j = jj;
      if (j >= geom->nfe2) {
	j = 0;
	i = ii;
      }
      if (ii > geom->nfe1) ii = i;
      if (wetf)n--;
    }
    c1++;
  }
  c1-=2;

  /*-----------------------------------------------------------------*/
  /* Save the surface cell locations in ws2                          */
  i = j = 0;
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (n = 1; n <= nwindows; n++) {
    wsizeS[n] = 0;
    cc = 1;
    wetf = 1;
    for (ii = i; ii < i + c1 + os; ii++) {
      if (ii >= geom->nfe1) continue;
      for (jj = j; jj < j + c1; jj++) {
	if (jj >= geom->nfe2) continue;
	c = geom->map[geom->nz-1][jj][ii];
	if (mask[c]) continue;
	if (wetc[c]) {
	  ws2[n][cc] = c;
	  mask[c] = 1;
	  wetf = 0;
	  wsizeS[n]++;
	  cc++;
	}
      }
    }
    if (ii >= geom->nfe1 && jj >= geom->nfe2) break;
    j = jj;
    if (j >= geom->nfe2) {
      j = 0;
      i = ii;
    }
    if (ii > geom->nfe1) ii = i;
    if (wetf) n--;
  }

  /*-----------------------------------------------------------------*/
  /* Collect any unallocated cells                                   */
  memset(mask, 0, (geom->enonS + 1) * sizeof(int));
  for (n = 1; n <= nwindows; n++) {
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];
      mask[c] = 1;
    }
  }
  m = 0;
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    if (mask[c] == 0) {
      wsizeS[nwindows]++;
      ws2[nwindows][wsizeS[nwindows]] = c;
      m++;
    }
  }
  if (DEBUG("init_w")) {
    for (n = 1; n <= nwindows; n++)
      dlog("init_w", "Blocking in e2 direction: window%d = %d cells\n", n, wsizeS[n]);
  }
  i_free_1d(wetc);
  i_free_1d(mask);
}

/* END window_cells_block_e2()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the global to local maps for each window, i.e.     */
/* derive the mapping between a global sparse coordinate and the     */
/* corresponding window, and local sparse coordinate in that window. */
/*-------------------------------------------------------------------*/
void get_gl_maps(geometry_t *geom,  /* Sparse global geometery       */
                 int nwindows,      /* Number of windows             */
                 int **ws2,         /* 2D wet cells in window wn     */
                 int *wsizeS        /* # 2D wet cells in window wn   */
  )
{
  int c, cc, c1, c2, n;         /* Sparse array counters             */
  int *w2;                      /* Dummy 2D surface cell vector      */

  for (n = 1; n <= nwindows; n++) {
    /* Allocate memory */
    w2 = i_alloc_1d(wsizeS[n] + 1);

    /* Get the global to local surface maps for the sub-surface      */
    /* layers filling the array in the (i,j,k) direction.            */
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];           /* Global coordinate in window n     */
      geom->fm[c].sc = cc;      /* Surface global to local map       */
      geom->fm[c].wn = n;       /* Surface global window             */
      w2[cc] = c;               /* Make a copy of ws2                */
    }

    /* Map all wet cells below the surface layer. This creates the   */
    /* map in the (i,j,k) direction, however, the order is not       */
    /* critical since geom->fm gets re-ordered in the routine        */
    /* get_local_maps().                                             */
    c2 = wsizeS[n];
    while (c2) {
      for (c1 = 1; c1 <= wsizeS[n]; c1++) {
        c = w2[c1];        /* Surface global coordinate in window wn */
        if (c) {           /* if(not a sediment cell) : ws2[c1]=0    */
          c = geom->zm1[c];  /* Global location one layer down       */
          w2[c1] = c;      /* Set the 2D local location to this cell */
          if (c != geom->zm1[c]) {  /* if(not a sediment cell)       */
            geom->fm[c].sc = cc;    /* Set global-local map location */
            geom->fm[c].wn = n; /* Set the global-local map window   */
            cc++;               /* Increment the local coordinate    */
          } else {              /* if(sediment cell)                 */
            if (w2[c1])
              c2--;             /* Set this column as done           */
            w2[c1] = 0;         /* Flag the local location sediment  */
	    /*
            geom->fm[c].sc = cc;
            cc++;
	    */
          }
        }
      }
    }
    i_free_1d(w2);
  }
}

/* END get_gl_maps()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the 3D local cells in window n given the surface   */
/* window partitioning. Fills the local sparse array in horizontal   */
/* slabs first.                                                      */
/*-------------------------------------------------------------------*/
void get_window_cells_h(geometry_t *geom, /* Global geometery        */
                        int wn,           /* Window number           */
                        int *wsa,       /* 3D wet cells in window wn */
                        int *wsize,     /* Size of wsa               */
                        int *ws2,       /* 2D wet cells in window wn */
                        int wsizeS      /* Size of ws2               */
			)
{
  int c, cc, c1, zm1;                  /* Cell indices / counters    */
  int *layer;                          /* Vertical layer map         */

  /* Initialise the vertical layer map                               */
  layer = i_alloc_1d(wsizeS + 1);
  memcpy(layer, ws2, (wsizeS + 1) * sizeof(int));
  layer[0] = 0;

  /* Initialise auxiliary maps for this window                        */
  for (c = 1; c <= geom->sgnum; c++)
    geom->fm[c].ac = 0;

  /* Get the global to local surface maps for the sub-surface layers  */
  /* filling the array in the (i,j,k) direction.                      */
  c1 = 1;                       /* 3D cell counter                    */
  wsize[wn] = 0;                /* Initialise number of 3D cells      */
  while (layer[0] < wsizeS) {
    for (cc = 1; cc <= wsizeS; cc++) {
      c = layer[cc];            /* Global coordinate in window wn     */
      zm1 = geom->zm1[c];
      if (c != zm1) {
        wsa[c1] = c;
        c1++;
        wsize[wn]++;
        layer[cc] = zm1;
      } else {
        if (layer[cc])
          layer[0]++;
        layer[cc] = 0;
      }
    }
  }
  i_free_1d(layer);
}

/* END get_window_cells_h()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to generate the local maps in the window structure        */
/*-------------------------------------------------------------------*/
void get_local_maps(geometry_t *geom,   /* Global geometry           */
		    geometry_t *window, /* Window geometry           */
                    int wn,             /* Window number             */
                    int *wsa,         /* 3D wet cells in window wn   */
                    int wsize,        /* # 3D wet cells in window wn */
                    int *ws2,         /* 2D wet cells in window wn   */
                    int wsize2D,      /* # 2D wet cells in window wn */
		    int cellf
  )
{
  int c, cc, i, j, jj;      /* Global coordinate, counter            */
  int e, ee, ae, ge, ges;   /* Edges                                 */
  int v, vv,vs,  av, vv1;   /* Vertices                              */
  int c2D;                  /* Local cell coordinate                 */
  int e2D;                  /* Local edge coordinate                 */
  int ac;                   /* Auxiliary centre                      */
  int cl, cu, cd, cs;       /* Local sparse coordinate               */
  int *wsar;                /* Reordered work sparse array           */
  int **nc2c;               /* Cell to cell spatial buffer           */
  int **ne2c;               /* Edge to cell spatial buffer           */
  int **ne2e;               /* Edge to edge spatial buffer           */
  int **nc2e;               /* Cell to edge spatial buffer           */
  int *nzp1, *nzm1;         /* Spatial buffer arrays for k+1 and k-1 */
  int *e2ee;                /* Global to local edge map              */
  int *v2vv;                /* Global to local vertex map            */
  int *sedm;
  int *mask;
  int verbose = 0;
  int npe;
  int checkf = 1;

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the buffer arrays. Since it is not known at */
  /* this stage how many auxiliary cells are to be included in the   */
  /* window arrays, set the size of these buffers to the sparse grid */
  /* size.                                                           */
  if (verbose) printf("Window %d\n", wn);
  wsar = i_alloc_1d(geom->szc);
  nc2c = i_alloc_2d(geom->szc, geom->npem+1);
  ne2c = i_alloc_2d(2, geom->sze);
  ne2e = i_alloc_2d(2, geom->sze);
  nc2e = i_alloc_2d(geom->szc, geom->npem+1);
  nzp1 = i_alloc_1d(geom->szc);
  nzm1 = i_alloc_1d(geom->szc);

  /*-----------------------------------------------------------------*/
  /* CENTRES                                                         */
  /*-----------------------------------------------------------------*/
  /* The following is defined for centres on each window:            */
  /* window->szc                                                     */
  /* window->szcS                                                    */
  /* window->c2c                                                     */
  /* window->zm1                                                     */
  /* window->zp1                                                     */
  /* Initialise the buffers. The local maps are made self mapping.   */
  /* The local wet work array is filled with global wet sparse       */
  /* locations; auxiliary global locations are appended below.       */
  window->szc = wsize + 1;
  window->szcS = wsize2D + 1;
  window->ewet = wsize;
  window->snon = wsize + 1;
  window->ewetS = wsize2D;
  window->snonS = wsize2D + 1;
  memset(geom->sask, 0, geom->sze * sizeof(int));

  for (c = 0; c < geom->szc; c++) {
    wsar[c] = 0;
    for(j = 1; j <= geom->npem; j++)
      nc2c[j][c] = c;
    nzp1[c] = c;
    nzm1[c] = c;
  }
  for (c = 1; c <= wsize2D; c++)
    wsar[c] = ws2[c];

  /*-----------------------------------------------------------------*/
  /* Get the work sparse array and related horizontal maps for the   */
  /* surface layer. Any sparse locations encountered in the surface  */
  /* layer are tagged as zero.                                       */
  c2D = 1;
  ac = wsize2D;
  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];

    if (c == ws2[c2D]) {
      cs = geom->m2d[c];
      local_map_build_c2c(geom, c, c2D, wn, geom->c2c, nc2c, wsar, &ac, laux, geom->npe[cs]);
      wsa[cc] = 0;
      c2D++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reorder the work array, wsar. This currently includes the wet   */
  /* surface cells followed by the surface auxiliary cells. Append   */
  /* all non-surface wet cells to the end of the surface layer here. */
  window->enonS = ac;
  window->szcS = ac + 1;
  for (cc = 1; cc <= wsize; cc++) {
    if (wsa[cc]) {
      ac++;
      wsar[ac] = wsa[cc];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reset the global-local map to that global sparse locations map  */
  /* correctly to the re-ordered local coordinates. This can only be */
  /* performed once the number of auxiliary cells in the surface     */
  /* layer are known.                                                */
  reorder_gl_map(geom, wn, wsa, wsize, wsar, window->enonS);

  /*-----------------------------------------------------------------*/
  /* Get the vertical maps for the surface layer wet cells. These    */
  /* should all map to other wet cells (downward maps), map to the   */
  /* sediment (downward maps in a water column one cell deep) or be  */
  /* self mapping (upward maps). Only downward maps to the sediment  */
  /* will increment ac; all other local sparse locations are already */
  /* defined.                                                        */
  for (cc = 1; cc <= wsize2D; cc++) {
    c = wsar[cc];
    local_map_build_z(geom, c, cc, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux);
    local_map_build_z(geom, c, cc, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux);
  }

  /*-----------------------------------------------------------------*/
  /* Get the work sparse array and related maps for the non-surface  */
  /* layer.                                                          */

  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];
    if (c) {
      cs = geom->m2d[c];
      cl = geom->fm[c].sc;
      local_map_build_c2c(geom, c, cl, wn, geom->c2c, nc2c, wsar, &ac, laux, geom->npe[cs]);
      local_map_build_z(geom, c, cl, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux);
      local_map_build_z(geom, c, cl, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add sediment locations beneath bottom auxiliary cells           */
  for (cs = 1; cs < geom->szcS; cs++) {
    cc = cs;                    /* cs = global surface location      */
    /* Search down for the sediment : cc = global sediment location  */
    while (cc != geom->zm1[cc]) {
      cc = geom->zm1[cc];
    }
    c = geom->zp1[cc];          /* c = global bottom location       */
    if (geom->fm[c].wn && geom->fm[c].wn != wn && geom->fm[cs].ac) {
      cu = cc;
      /* If a step in the bathymetry occurs at the wet-auxiliary    */
      /* face then the auxiliary sediment and non-sediment cell may */
      /* be isolated (non-contiguous maps). Therefore need to loop  */
      /* up from the sediment cell until the non-sediment auxiliary */
      /* cell located. This is necessary if cell thickness arrays   */
      /* (dz[]) are set at auxiliary locations since the bottom of  */
      /* the water column is located where the downward map is self-*/
      /* mapping in this case (i.e. the auxiliary sediment cell     */
      /* must be the only cell which is self mapping).              */
      while (geom->fm[cu].ac == 0) {
        ac++;
        geom->fm[cu].ac = ac;
        wsar[ac] = cu;
        cu = geom->zp1[cu];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the local wet array and local maps          */
  window->enon = ac;
  window->szc = ac + 1;
  alloc_geom_us(window, CENTRE_A);
  if (verbose) printf("# centres: 2D = %d, 3D = %d\n", window->szcS-1, window->szc-1);

  /*-----------------------------------------------------------------*/
  /* Fill the window structure with sparse locations (global for wsa */
  /* and local for the maps) from the buffer arrays.                 */
  if (ac >= geom->szc)
    hd_quit("get_local_maps: Error at %d %d\n", wn, cc);

  for (cc = 1; cc <= ac; cc++) {
    c = wsar[cc];
    cs = geom->m2d[c];
    npe = geom->npe[cs];

    window->wsa[cc] = c;
    
    if (cc <= window->enonS) window->npe[cc] = npe;
    for(j = 1; j <= npe; j++)
      window->c2c[j][cc] = cc;
    for (j = 1; j <= npe; j++)
      window->c2c[j][cc] = nc2c[j][cc];
    window->zp1[cc] = nzp1[cc];
    window->zm1[cc] = nzm1[cc];
  }

  /*-----------------------------------------------------------------*/
  /* Get the vertical maps for all auxiliary cells. These local      */
  /* locations are derived through the global vertical maps.         */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];

    /* Auxiliary cells and ghost cells                               */
    if ((geom->fm[c].wn != wn && geom->fm[geom->zp1[c]].wn != wn) ||
        (geom->fm[c].wn == 0 && geom->fm[geom->zp1[c]].wn == wn)) {
      cl = geom->fm[c].ac;
      if (cl) {
        /* Upward maps for auxiliary cells and ghost cells           */
        cu = geom->fm[geom->zp1[c]].ac;
        /* Upward maps for sub-surface ghost cells (i.e. ghost       */
        /* cells where the cell above is wet).                       */
        if (!cu)
          cu = geom->fm[geom->zp1[c]].sc;
        cd = geom->fm[geom->zm1[c]].ac;
        if (cu)
          window->zp1[cl] = cu;
        if (cd)
          window->zm1[cl] = cd;
      }
    }
  }

  /* Get the undefined lateral maps for 2D auxiliary cells           */
  for (cc = window->snonS; cc <= window->enonS; cc++)
    aux_maps(geom, window, cc);
  /* Get the undefined lateral maps for 3D auxiliary cells           */
  for (cc = window->snon; cc <= window->enon; cc++)
    aux_maps(geom, window, cc);

  /*-----------------------------------------------------------------*/
  /* Make the global to local map for this window. This differs from */
  /* the global to local map in geom in that it is defined only over */
  /* all global cells in the window, including ghost cells which are */
  /* associated with zero window and local coordinate in geom->fm.   */
  window->fm =
    (global_map_t *)malloc(sizeof(global_map_t) * (geom->enon + 1));
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    window->fm[c].wn = window->wn;
    window->fm[c].sc = cc;
    if (geom->fm[c].wn == 0)
      window->fm[c].ac = 0;
    else if (geom->fm[c].wn == window->wn)
      window->fm[c].ac = 1;
    else
      window->fm[c].ac = 2;
    if (verbose && cc <= window->enonS) {
      if (window->fm[c].ac == 0)
	printf("centres: l=%d(ghost) wn=%d g=%d\n",cc, wn, c);
      if (window->fm[c].ac == 1)
	printf("centres: l=%d(wet) wn=%d g=%d\n",cc, wn, c);
      if (window->fm[c].ac == 2)
	printf("centres: l=%d(aux) wn=%d g=%d\n",cc, wn, c);
    }
  }

  /* Deallocate the buffer arrays                                    */
  i_free_1d(wsar);
  i_free_2d(nc2c);
  i_free_1d(nzp1);
  i_free_1d(nzm1);

  if (checkf) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      cs = geom->m2d[c];
      if (geom->fm[cs].wn == wn) {
	for (j = 1; j <= geom->npe[cs]; j++) {
	  int cl = window->c2c[j][cc];
	  int clg = window->wsa[cl];
	  int cn = geom->c2c[j][cs];
	  if (clg != cn) hd_quit("get_local_maps: window%d local cell %d(%d) neighbour %d(%d) (j=%d) does not map to global cell %d ([%f %f]-[%f %f])\n",
				 wn, cc, cs, cl, cn, j, clg, geom->cellx[cs], geom->celly[cs], geom->cellx[cn], geom->celly[cn]);
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* EDGES                                                           */
  /*-----------------------------------------------------------------*/
  /* The following is defined for edges on each window:              */
  /* window->sze                                                     */
  /* window->szeS                                                    */
  /* window->wse                                                     */
  /* window->c2e                                                     */
  /* window->e2c                                                     */
  /* window->e2e                                                     */
  /* window->eSc                                                     */
  /* window->eSe                                                     */
  /* window->wAe                                                     */
  /* window->ep                                                      */
  /* window->em                                                      */
  /* window->zm1e                                                    */
  /* window->zp1e                                                    */
  /* window->e2k                                                     */
  /* Count the edges in the window.                                  */
  window->szeS = window->sze = 1;
  mask = i_alloc_1d(geom->sze);
  sedm = i_alloc_1d(geom->szeS);
  memset(mask, 0, geom->sze * sizeof(int));
  memset(sedm, 0, geom->szeS * sizeof(int));

  /* First get all edges surrounding wet and auxiliary cells in a    */
  /* window.                                                         */
  /* Note that wet edges in a window are defined as those edges      */
  /* whose global sign vector (eSc) is equal to -1, i.e. those       */
  /* vectors that point in towards the cell centre.                  */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    if (!geom->fm[c].wn) continue;   /* Ghost cell                   */
    for (j = 1; j <= geom->npe[cs]; j++) {
      /*if (geom->eSc[j][cs] == -1) {*/
      e = geom->c2e[j][c];
      if (!mask[e]) {
	if (cc <= window->enonS)
	  window->szeS++;
	window->sze++;
	mask[e] = 1;
      }
    }
  }

  /* Next get and edges that were missed in the above loop that      */
  /* correspond to all those edges connecting each vertex of cells   */
  /* that are adjacent to a wet edge in the window. These are        */
  /* required for voriticity computations in the nonlinear Coriolis  */
  /* momentum advection, and depending on the window arrangement may */
  /* have been missed above. If these edges are not associated with  */
  /* centres in the window, then e2c maps cannot be set.             */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    if (!geom->fm[c].wn) continue;   /* Ghost cell                   */
    if (geom->fm[c].wn != wn) continue; /* Centre not in window      */
    for (j = 1; j <= geom->npe[cs]; j++) {
      if (geom->eSc[j][cs] == -1) {  /* Edge is in the window        */
	ge = geom->c2e[j][c];
	for (i = 0; i <= 1; i++) {
	  v = geom->e2v[ge][i];      /* Vertex of e                  */
	  vs = geom->m2dv[v];
	  for (jj = 1; jj <= geom->nve[vs]; jj++) {
	    e = geom->v2e[v][jj];    /* Edge connected to v          */
	    if (e && !mask[e]) {
	      if (cc <= window->enonS)
		window->szeS++;
	      window->sze++;
	      mask[e] = 1;
	    }
	  }
	}
      }
    }
  }

  ges = window->sze;
  window->sze += window->szeS;  /* Sediments                         */

  /* Allocate                                                        */
  alloc_geom_us(window, EDGE_A);
  window->e2c = i_alloc_2d(2, window->sze);
  window->e2e = i_alloc_2d(2, window->sze);
  window->c2e = i_alloc_2d(window->szc, window->npem+1);
  window->eSc = i_alloc_2d(window->szcS, window->npem + 1);
  window->eSe = i_alloc_2d(window->sze, window->neem + 1);
  window->wAe = d_alloc_2d(window->sze, window->neem + 1);
  window->zp1e = i_alloc_1d(window->sze);
  window->zm1e = i_alloc_1d(window->sze);
  window->ep = i_alloc_1d(window->sze);
  window->em = i_alloc_1d(window->sze);
  window->e2k = i_alloc_1d(window->sze);
  e2ee = i_alloc_1d(geom->sze);
  wsar = i_alloc_1d(geom->sze);
  memset(wsar, 0, geom->sze * sizeof(int));
  if (verbose) printf("# edges: 2D = %d, 3D = %d\n", window->szeS-1, window->sze-1);

  /* Set the maps.                                                   */
  /* First get all edges surrounding wet and auxiliary cells in a    */
  /* window.                                                         */
  ee = 1;
  memset(mask, 0, geom->sze * sizeof(int));
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];

    if (!geom->fm[c].wn) continue;   /* Ghost cell                   */
    npe = geom->npe[cs];
    for (j = 1; j <= npe; j++) {
      int cn = geom->c2c[j][c];     /* Neighbour of c in direction j */
      e = geom->c2e[j][c];          /* Edge in direction j           */

      if (!mask[e]) {
	window->wse[ee] = e;
	e2ee[e] = ee;
	geom->g2we[wn][e] = ee;
	/* Note: fm.ec contains the wet edge only in window wn.      */
	/* These are unique for window wn.                           */
	if (geom->fm[c].wn == wn) {
	  if (geom->fm[cn].wn == 0)       /* Ghost cells included    */
	    geom->fm[e].ec = ee;
	  else if (geom->fm[cn].wn == wn) /* Wet cells included      */
	    geom->fm[e].ec = ee;
	  else {
	    if (geom->eSc[j][cs] == -1)   /* Neighbour is auxiliary  */
	      geom->fm[e].ec = ee;        /* but vector is inward.   */
	  }
	}
	wsar[e] = 1;
	if (geom->eSc[j][cs] == -1) {
	  window->e2c[ee][0] = cc;
	  window->e2c[ee][1] = window->c2c[j][cc];
	  window->e2e[ee][0] = j;
	  window->e2e[ee][1] = geom->e2e[e][1];
	} else {
	  window->e2c[ee][1] = cc;
	  window->e2c[ee][0] = window->c2c[j][cc];
	  window->e2e[ee][0] = geom->e2e[e][0];
	  window->e2e[ee][1] = j;
	}

	window->e2k[ee] = geom->e2k[e];
	if (verbose && cc <= window->enonS)
	  printf("edges (centre): l=%d g=%d c0=%d(%d) c1=%d(%d)\n",ee, window->wse[ee], window->e2c[ee][0], 
		 window->wsa[window->e2c[ee][0]], window->e2c[ee][1], window->wsa[window->e2c[ee][1]]);
	mask[e] = 1;
	ee++;
      }
      window->c2e[j][cc] = e2ee[e];
      if (cc <= window->enonS)
	window->eSc[j][cc] = geom->eSc[j][cs];
    }
  }

  /* Next get and edges that were missed in the above loop that      */
  /* correspond to all those edges connecting each vertex of cells   */
  /* that are adjacent to a wet edge in the window.                  */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    if (!geom->fm[c].wn) continue;   /* Ghost cell                   */
    if (geom->fm[c].wn != wn) continue; /* Centre not in window      */
    npe = geom->npe[cs];
    for (j = 1; j <= npe; j++) {
      if (geom->eSc[j][cs] == -1) {  /* Edge is in the window        */
	int cn = geom->c2c[j][c];    /* Neighbour of c, direction j  */
	ge = geom->c2e[j][c];        /* Global edge                  */
	for (i = 0; i <= 1; i++) {
	  v = geom->e2v[ge][i];      /* Vertex of e                  */
	  vs = geom->m2dv[v];
	  for (jj = 1; jj <= geom->nve[vs]; jj++) {
	    e = geom->v2e[v][jj];   /* Edge connected to v          */
	    
	    if (e && !mask[e]) {
	      window->wse[ee] = e;
	      e2ee[e] = ee;
	      geom->g2we[wn][e] = ee;
	      wsar[e] = 1;

	      window->e2e[ee][0] = j;
	      window->e2e[ee][1] = geom->e2e[e][1];
	      window->e2k[ee] = geom->e2k[e];
	      if (cc <= window->enonS) {
		if (verbose)
		  printf("edges (vertex): l=%d g=%d c0=%d(%d) c1=%d(%d)\n",ee, window->wse[ee], window->e2c[ee][0], 
		       window->wsa[window->e2c[ee][0]], window->e2c[ee][1], window->wsa[window->e2c[ee][1]]);
	      }
	      mask[e] = 1;
	      ee++;
	    }
	  }
	}
      }
    }
  }

  /* Ghost cells, treated the same as the master in pp_us()          */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    if ((cd = geom->wgst[c])) {
      cs = geom->m2d[cd];
      for (j = 1; j <= geom->npe[cs]; j++) {
	if ((e = geom->c2e[j][c]))
	  window->c2e[j][cc] = e2ee[e];
      }
    }
  }

  /* Set self-mapping                                                */
  for (ee = 1; ee < window->sze; ee++) {
    window->zp1e[ee] = ee;
    window->zm1e[ee] = ee;
    window->ep[ee] = ee;
    window->em[ee] = ee;
  }

  /* Get the edges surrounding an edge map                           */
  for (ee = 1; ee < window->szeS; ee++) {
    e = window->wse[ee];
    window->nee[ee] = geom->nee[e];   
  }

  /* Set the edge spatial maps                                       */
  jj = 1;
  memset(mask, 0, geom->sze * sizeof(int));
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    npe = geom->npe[cs];
    for (j = 1; j <= npe; j++) {
      int es;
      e = geom->c2e[j][c];
      es = geom->m2de[e];
      if (!e) continue;
      if (!mask[e]) {
	ee = e2ee[e];

	/* Edge maps and weights                                     */
	for (i = 1; i <= geom->nee[es]; i++) {
	  ae = geom->eSe[i][e];
	  window->eSe[i][ee] = e2ee[ae];
	  window->wAe[i][ee] = geom->wAe[i][e];
	}
	/* Vertical maps                                             */
	ae = geom->zp1e[e];
	window->zp1e[ee] = e2ee[ae];
	ae = geom->zm1e[e];
	window->zm1e[ee] = e2ee[ae];

	/* Sediments under wet cells. When a sediment edge is        */
	/* assigned, set the mask sedm. This prevents additional     */
	/* sediment edges created from access via cells deeper than  */
	/* the shallowest cell sharing an edge.                      */
	if (!sedm[es] && !geom->wgst[c] && geom->fm[c].wn && geom->zm1e[ae] == ae) {
	  window->zm1e[ee] = ges;
	  if (ges>window->sze)hd_quit("get_local_maps: wn%d sed size %d > %d\n",window->wn,ges,window->sze);
	  window->zm1e[ges] = ges;
	  window->zp1e[ges] = ee;
	  sedm[geom->m2de[e]] = ee;
	  ges++;
	  jj++;
	}

	/* Horizontal maps                                           */
	ae = geom->ep[e];
	if (wsar[ae])
	  window->ep[ee] = e2ee[ae];
	ae = geom->em[e];
	if (wsar[ae])
	  window->em[ee] = e2ee[ae];
	/*
	if (verbose && cc <= window->enonS)
	  printf("edges: l=%d g=%d c0=%d(%d) c1=%d(%d)\n",ee, window->wse[ee], window->e2c[ee][0], 
		 window->wsa[window->e2c[ee][0]], window->e2c[ee][1], window->wsa[window->e2c[ee][1]]);
	*/
	mask[e] = 1;
      }
    }
  }

  i_free_1d(wsar);
  i_free_1d(mask);
  i_free_1d(sedm);

  /*-----------------------------------------------------------------*/
  /* VERTICES                                                        */
  /*-----------------------------------------------------------------*/
  /* The following is defined for vertices on each window:           */
  /* window->szv                                                     */
  /* window->szvS                                                    */
  /* window->szm                                                     */
  /* window->szmS                                                    */
  /* window->c2v                                                     */
  /* window->e2v                                                     */
  /* window->v2c                                                     */
  /* window->v2e                                                     */
  /* window->eSv                                                     */
  /* window->vIc                                                     */
  /* window->zm1v                                                    */
  /* window->zp1v                                                    */
  /* window->dualarea                                                */
  /* window->dualareap                                               */
  /* Count the vertices in the window.                               */
  window->szvS = window->szv = 1;
  wsar = i_alloc_1d(geom->szv);
  v2vv = i_alloc_1d(geom->szv);
  memset(wsar, 0, geom->szv * sizeof(int));

  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    /* Get the vertices in the window. Vertices of every cell centre */
    /* in the window are included.                                   */
    for (j = 1; j <= geom->npe[cs]; j++) {
      v = geom->c2v[j][c];
      if (v) {
	if (!wsar[v]) {
	  wsar[v] = 1;
	  if (cc <= window->enonS)
	    window->szvS++;
	  window->szv++;
	}
      }
    }
  }

  /* Allocate                                                        */
  alloc_geom_us(window, VERTEX_A);
  window->c2v = i_alloc_2d(window->szc, window->npem+1);
  window->v2c = i_alloc_2d(window->nvcm+1, window->szv);
  window->v2e = i_alloc_2d(window->nvem+1, window->szv);
  window->e2v = i_alloc_2d(2, window->sze);
  window->zp1v = i_alloc_1d(window->szv);
  window->zm1v = i_alloc_1d(window->szv);
  window->vIc = i_alloc_2d(window->szcS, window->npem + 1);
  window->eSv = i_alloc_2d(window->szv, window->nvem + 1);
  memset(wsar, 0, geom->szv * sizeof(int));
  for (vv = 1; vv < window->szv; vv++)
    for (i = 1; i <= window->nvcm; i++)
      window->v2c[vv][i] = 0;
  if (verbose) printf("# vertices: 2D = %d, 3D = %d\n", window->szvS-1, window->szv-1);

  /* Set the maps                                                    */
  vv = 1;
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];

    for (j = 1; j <= geom->npe[cs]; j++) {
      int ic;
      v = geom->c2v[j][c];
      vs = geom->m2dv[v];

      if (v) {
	/* Vertices can be accessed from multiple centres, hence to  */
	/* ensure the local vertex is unique, set a mask (wsar) of   */
	/* the global vertex (v) to the local vertex (vv), and use   */
	/* this value when the vertex is subsequently accessed.      */
	/* First time v is accessed: vv1 = vv                        */
	/* Subsequent times vv is accessed: vv1 = wsar[v]            */
	vv1 = (wsar[v]) ? wsar[v] : vv;
	window->wsv[vv1] = v;
	v2vv[v] = vv1;
	geom->g2wv[wn][v] = vv1;

	/* Note: fm.vc contains the wet vertices only in window wn.  */
	/* These are unique for window wn.                           */
	if (geom->fm[c].wn == wn) {
	  geom->fm[v].vc = vv1;
	}
	window->c2v[j][cc] = vv1;

	if (cc <= window->enonS) {
	  window->vIc[j][cc] = geom->vIc[j][c];
	  window->dualarea[vv1] = geom->dualarea[v];
	  window->nve[vv1] = geom->nve[v];
	  window->nvc[vv1] = geom->nvc[v];
	  for (i = 1; i <= window->nvc[vv1]; i++)
	    window->dualareap[vv1][i] = geom->dualareap[v][i];
	}

	for (i = 1; i <= geom->nvc[vs]; i++) {
	  /* Note: maps from v to ghost cell centres adjacent to     */
	  /* auxiliary cells may not be set if the ghost cell is not */
	  /* is not directional to a wet cell. These cells should    */
	  /* never be accessed however. They manifest as v2c[cc]=0.  */
	  if (c == geom->v2c[v][i] && !window->v2c[vv1][i]) {
	    window->v2c[vv1][i] = cc;
	    ic = i;
	  }
	}
	for (i = 1; i <= geom->nve[vs]; i++) {
	  e = geom->v2e[v][i];
	  if (e && !wsar[v]) {
	    ee = e2ee[e];
	    window->v2e[vv1][i] = ee;
	    if (v == geom->e2v[e][0]) {
	      window->e2v[ee][0] = vv1;
	      window->eSv[i][vv1] = 1;
	    }
	    if (v == geom->e2v[e][1]) {
	      window->e2v[ee][1] = vv1;
	      window->eSv[i][vv1] = -1;
	    }
	  }
	}
	if (verbose && cc <= window->enonS && window->v2c[vv1][ic])
	  printf("vertices: l=%d g=%d v2c(%d)=%d(%d)\n",vv1, window->wsv[vv1], ic, window->v2c[vv1][ic], c);
	if (!wsar[v]) {
	  wsar[v] = vv;
	  vv++;
	}
      }
    }
  }

  /* Set self-mapping                                                */
  for (vv = 1; vv < window->szv; vv++) {
    window->zp1v[vv] = vv;
    window->zm1v[vv] = vv;
  }
  /* Set the edge spatial maps                                       */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    cs = geom->m2d[c];
    for (j = 1; j <= geom->npe[cs]; j++) {
      v = geom->c2v[j][c];
      vv = v2vv[v];
      
      if (v) {
	/* Vertical maps                                             */
	av = geom->zp1v[v];
	window->zp1v[vv] = v2vv[av];
	av = geom->zm1v[v];
	/* Global vertex maps map downward to the depest edge. For   */
	/* windows the deepest edge may lie in another window (e.g.  */
	/* for vertices on the outmost auxiliary cell) and have no   */
	/* associated coordinate for that window, even though the    */
	/* global coordate does. Do not map downward in these cases. */
	if (v2vv[av])
	  window->zm1v[vv] = v2vv[av];
      }
    }
  }

  i_free_1d(wsar);
  i_free_1d(e2ee);
  i_free_1d(v2vv);

  window->szm = max(max(window->szc, window->sze), window->szv);
  window->szmS = max(max(window->szcS, window->szeS), window->szvS);
  alloc_geom_us(window, MAP_A);
}

/* END get_local_maps()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the local spatial maps. If a wet local cell maps   */
/* to an auxiliary cell, then the auxiliary cell maps to other       */
/* auxiliary cells saux times. This is required to handle the higher */
/* order schemes. The value of saux is dependent on the order of the */
/* scheme used and if the ULTIMATE filter is used (e.g. for a second */
/* order scheme saux=1 and second order+ULTIMATE saux=2). If cells   */
/* are self-mapping (i.e. boundary cells) then the auxiliary mapping */
/* terminates.                                                       */
/*-------------------------------------------------------------------*/
void local_map_build_c2c(geometry_t *geom, /* Global geometry        */
			 int c,     /* Global sparse coordinate      */
			 int cc,    /* Local sparse coordinate       */
			 int wn,    /* Window number                 */
			 int **map, /* Global spatial map            */
			 int **nmap,/* Local spatial map             */
			 int *wsa,  /* Local 3D wet cell array       */
			 int *ac,   /* Local auxiliary cell location */
			 int nsaux, /* # cells to include laterally  */
			 int npe    /* Number of vertices            */
			 )
{
  int n, j, jj, jop;        /* Counters / edge directions            */
  int jv, v, vs;            /* Vertices                              */
  int cgm, cg, cgms;        /* Global locations of mapped coordinate */
  int cax, cl;              /* Local location of mapped coordinate   */
  int e;                    /* Global edge                           */
  int *mask = geom->sask;
  int verbose = 0;

  /*-----------------------------------------------------------------*/
  /* First get the cells in all edge directions                      */
  for (j = 1; j <= npe; j++) {  /* Edge number                       */

    /* Get the global map of the global sparse coordinate            */
    cgm = map[j][c];
    cgms = geom->m2d[cgm];
    e = geom->c2e[j][c];
    if (!mask[e]) {

      /* Get the local map of the local sparse coordinate            */
      nmap[j][cc] = cax = reset_map(geom, cgm, wn, wsa, ac);
      if (cc != cax) {
	jop = jow(geom, c, j);
	nmap[jop][cax] = cc;
      }

      mask[e] = 1;
      if (verbose && c == geom->m2d[c])
	printf("c2c wn=%d c=%d cgm=%d j=%d cl=%d ac=%d\n",wn,c,cgm,j,cax,*ac);

      /* If the local map is an auxiliary cell, then map this cell a */
      /* further nsaux times in the direction of map[].              */
      if (geom->fm[cgm].ac != 0) {
	for (n = 1; n < nsaux; n++) {
	  
	  if (cgm != map[j][cgm]) {   /* Check for self-mapping      */
	    cg = cgm;
	    cgm = map[j][cgm];
	    e = geom->c2e[j][cg];
	    if (!mask[e]) {
	      nmap[j][cax] = cl = reset_map(geom, cgm, wn, wsa, ac);
	      if (cax != nmap[j][cax]) {
		int cr = (geom->wgst[cg]) ? geom->wgst[cg] : cg;
		jop = jow(geom, cr, j);
		nmap[jop][cl] = cax;
	      }
	      cax = cl;
	      mask[e] = 1;
	    }
	    if (verbose && c == geom->m2d[c])
	      printf("aux wn=%d c=%d cgm=%d j=%d cl=%d\n",wn,cg,cgm,j,cax);

	  } else
	  continue;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Next include centres surrounding each vertex of a given cell.   */
  /* These cells are required to normalize the relative vorticity    */
  /* computed at vertices.                                           */
  /* Cells around vertices of c.                                     */
  local_map_build_v2c(geom, c, cc, wn, nmap, wsa, ac, nsaux);

  /* Cells around vertices of mapped centres.                        */
  for (j = 1; j <= npe; j++) {  /* Edge number                       */

    /* Get the global map of the global sparse coordinate            */
    cgm = map[j][c];
    cax = nmap[j][cc];

    /* Cells around vertices of cgm                                  */
    local_map_build_v2c(geom, cgm, cax, wn, nmap, wsa, ac, nsaux);

    /* If the local map is an auxiliary cell, then map this cell a   */
    /* further nsaux times in the direction of map[].                */
    if (geom->fm[cgm].ac != 0) {
      for (n = 1; n < nsaux; n++) {
	  
	if (cgm != map[j][cgm]) {     /* Check for self-mapping      */
	  cgm = map[j][cgm];
	  cax = nmap[j][cax];

	  /* Cells around vertices of cgm                            */
	  local_map_build_v2c(geom, cgm, cax, wn, nmap, wsa, ac, nsaux);
	} else
	  continue;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the cell centres required for high order advection          */
  for (j = 1; j <= npe; j++) {  /* Edge number                       */
    int es, eoe, m;
    int cax1, cax2, j1, j2;
    /* Get the global map of the global sparse coordinate            */
    e = geom->c2e[j][c];        /* Global edge                       */
    es = geom->m2de[e];         /* Surface global edge               */

    for (n = 1; n <= geom->nee[es]; n++) {
      eoe = geom->eSe[n][e];
      if (!eoe || mask[eoe]) continue;      
      cgm = geom->e2c[eoe][0];
      j1 = geom->e2e[eoe][0];
      cax1 = reset_map(geom, cgm, wn, wsa, ac);

      cgm = geom->e2c[eoe][1];
      j2 = geom->e2e[eoe][1];
      cax2 = reset_map(geom, cgm, wn, wsa, ac);

      nmap[j1][cax1] = cax2;
      nmap[j2][cax2] = cax1;
    }
    mask[eoe] = 1;
  }
}

/* END local_map_build_c2c()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* The centres to cell c are found by locating all centres           */
/* surrounding every vertex of c. This is done for an expanding      */
/* circuit of centres, saux times, around c.                         */
/*-------------------------------------------------------------------*/
void local_map_build_v2c(geometry_t *geom, /* Global geometry        */
			 int c,     /* Global sparse coordinate      */
			 int cc,    /* Local sparse coordinate       */
			 int wn,    /* Window number                 */
			 int **nmap,/* Local spatial map             */
			 int *wsa,  /* Local 3D wet cell array       */
			 int *ac,   /* Local auxiliary cell location */
			 int nsaux /* # cells to include laterally  */
			 )
{
  int n, nc, jc, jv;        /* Counters                              */
  int j, jop;               /* Direction of centres                  */
  int v, vs;                /* Vertices                              */
  int cs;                   /* Surface global centre                 */
  int cgm;                  /* Global centres                        */
  int cax;                  /* Local centres                         */
  int e;                    /* Global edge                           */
  int npe;                  /* # cells surrounding c                 */
  int domap = 0;            /* Set c2c mappings between cells        */
  int *mask = geom->sask;
  int verbose = 0;


  cs = geom->m2d[c];
  npe = geom->npe[cs];

  /* Cells of vertices surrounding c                                 */
  for (jc = 1; jc <= npe; jc++) {      /* Vertex counter             */
    v = geom->c2v[jc][c];              /* Global vertex to c         */
    if (v) {
      vs = geom->m2dv[v];              /* Surface global vertex      */
      for (jv = 1; jv <= geom->nvc[vs]; jv++) {
	cgm = geom->v2c[v][jv];        /* Global centre to v         */ 
	if (cgm && cgm != c) {
	  cax = reset_map(geom, cgm, wn, wsa, ac);
	  j = joc(geom, c, cgm);       /* Direction of cgm from c    */
	  if(verbose && c==cs && !j)
	    printf("v2c_um wn=%d c=%d cgm=%d v=%d cl=%d j=%d\n",wn,c,cgm,v,cgm,cax,j);
	  if (j && domap) {
	    e = geom->c2e[j][c];       /* Edge in j direction        */
	    if (!mask[e]) {
	      /* Get the local map of the local sparse coordinate    */
	      if (domap) {
		nmap[j][cc] = cax;
		if (cc != cax) {
		  int cr = (geom->wgst[c]) ? geom->wgst[c] : c;
		  jop = jow(geom, cr, j); /* Direction of c from cgm */
		  nmap[jop][cax] = cc ;
		}
	      }
	      mask[e] = 1;
	      if(verbose && c==cs)
		printf("v2c_m wn=%d c=%d cgm=%d v=%d cl=%d j=%d\n",wn,c,cgm,v,cgm,cax,j);
	    }
	  }
	}
      }
    }
  }
}

/* END local_map_build_v2c()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the opposite edge number to edge j                        */
/*-------------------------------------------------------------------*/
int jow(geometry_t *geom, int c, int j)
{
  int jj;
  int cn = geom->c2c[j][c];

  /* Ghosts are treated consistently as in pp_us()                   */
  if (geom->wgst[cn])
    return(jo(j, geom->npe[geom->m2d[c]]));

  for (jj = 1; jj <= geom->npe[geom->m2d[cn]]; jj++) {
    if (c == geom->c2c[jj][cn]) {
      return(jj);
    }
  }
  return(jo(j, geom->npe[geom->m2d[c]]));
}

/* END jow()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns the direction of cell cn from c                           */
/*-------------------------------------------------------------------*/
int joc(geometry_t *geom, int c, int cn)
{
  int j;
  int cs = geom->m2d[c];

  for(j = 1; j <= geom->npe[cs]; j++)
    if (cn == geom->c2c[j][c])
      return(j);
  return(0);
}

/* END joc()                                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Vertical mappings                                                 */
/*-------------------------------------------------------------------*/
void local_map_build_z(geometry_t *geom,   /* Global geometry        */
		       int c,     /* Global sparse coordinate        */
		       int cc,    /* Local sparse coordinate         */
		       int wn,    /* Window number                   */
		       int *map,  /* Global spatial map              */
		       int *nmap, /* Local spatial map               */
		       int *rmap, /* Reverse local spatial map       */
		       int *wsa,  /* Local 3D wet cell sparse array  */
		       int *ac,   /* Local auxiliary cell location   */
		       int nsaux  /* # cells to include laterally    */
		       )
{
  int n;                    /* Counter                               */
  int cgm;                  /* Global locations of mapped coordinate */
  int cax;                  /* Local location of mapped coordinate   */

  /* Get the global map of the global sparse coordinate              */
  cgm = map[c];

  /* Get the local map of the local sparse coordinate                */
  nmap[cc] = cax = reset_map(geom, cgm, wn, wsa, ac);
  if (cc != cax)
    rmap[cax] = cc;

  /* If the local map is an auxiliary cell, then map this cell a     */
  /* further nsaux times in the direction of map[].                  */
  /*if (cax == *ac || geom->fm[cgm].ac != 0) {*/
  if (geom->fm[cgm].ac != 0) {
    for (n = 1; n < nsaux; n++) {
      if (cgm != map[cgm]) {    /* Check for self-mapping locations  */
        nmap[cax] = reset_map(geom, cgm, wn, wsa, ac);
        if (cax != nmap[cax])
          rmap[nmap[cax]] = cax;
        cax = nmap[cax];
      } else
        return;
    }
  }
}


/* END local_map_build_z()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to check if a given global sparse coordinate is a member  */
/* of a given window wn. If so, then return the local sparse         */
/* coordinate in that window. If not, then add a new coordinate to   */
/* the window, set the reverse map of the new local coordinate to    */
/* the global coordinate and return the new local coordinate.        */
/*-------------------------------------------------------------------*/
int reset_map(geometry_t *geom,  /* Global geometry              */
		  int cm,   /* Mapped global sparse coordinate       */
		  int wn,   /* Window number                         */
		  int *wsa, /* Local 3D wet cell sparse array        */
		  int *al   /* Local location of auxiliary cell      */
	      )
{
  int sc;                 /* Local sparse coordinate of window       */
  int wl;                 /* Window number of cm                     */

  sc = geom->fm[cm].sc;
  wl = geom->fm[cm].wn;

  if (wl != wn) {
    /* Add auxiliary cells to the 3D array                           */
    if (geom->fm[cm].ac) {
      return (geom->fm[cm].ac);
    } else {
      *al += 1;
      wsa[*al] = cm;          /* Include the auxiliary cell in wsa[] */
      geom->fm[cm].ac = *al;  /* Map global cell to auxiliary cell   */
    }
    return (*al);
  } else
    return (sc);
}

/* END reset_map()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to extract the local sparse coordinates for the work      */
/* array of cells to process from the global sparse work array.      */
/* These vectors are arranged so that the first v3_ cells are wet    */
/* cells followed by b3_ minus v3_ normal open boundary cells        */
/* (tangential open boundary cells are included as wet cells)        */
/* followed by a3_ minus  b3_ auxiliary cells (this is followed by a */
/* further x3_ minus a3_ southern edge auxiliary cells for u1        */
/* velocity and western edge auxiliary cells for u2 velocity)        */
/* followed by n3_ minus a3_ ghost cells (for both 3D and 2D         */
/* vectors). Note that special treatment of e1 ghost cells on        */
/* western edges and e2 ghost cells on southern edges must be        */
/* applied since these faces are wet and thus part of the window.    */
/* Note that a cell is a ghost cell if it is associated with window  */
/* = 0.                                                              */
/*-------------------------------------------------------------------*/
void get_local_wsc(geometry_t *geom,    /* Global geometry           */
		   int *vec,            /* Global work sparse array  */
                   int nvec,            /* Size of vec               */
                   int nvec2D,          /* Last 2D cell loc in vec   */
                   geometry_t **window, /* Window data structure     */
                   int nwindows         /* Number of windows         */
		   )
{
  int *wntmp;                   /* Temporary buffer */
  int *wntmp2D;                 /* Temporary buffer */
  int c, cs, cc;                /* Centre counters */
  int ee, e, es;                /* Edge counters */
  int n, nn, wn;                /* Window counters */
  int ac, cgm;                  /* Sparse coordinates */
  int **mask;                   /* Set to 1 when cell is processed */
  int *wm;                      /* Window / ghost cell mask */
  int ms;                       /* Size of mask */
  int sc;                       /* Wet cell in window */
  int *gc, *gc2D;               /* Ghost cell in window */
  int *ec, *ec2D;               /* Extra auxiliary cells for mode=1 & 2 */
  int nf;
  int v, vs, i, j, eoe;
  int mode = 0;
  int checkf = 0;
  int verbose = 0;
  int **g2w;

  /*-----------------------------------------------------------------*/
  /* Initialise the counters */
  wntmp = i_alloc_1d(nwindows + 1);
  wntmp2D = i_alloc_1d(nwindows + 1);
  gc = i_alloc_1d(nwindows + 1);
  gc2D = i_alloc_1d(nwindows + 1);
  ec = i_alloc_1d(nwindows + 1);
  ec2D = i_alloc_1d(nwindows + 1);
  wm = i_alloc_1d(geom->szm);
  g2w = i_alloc_2d(geom->szc, nwindows+1);

  for (n = 1; n <= nwindows; n++) {
    wntmp[n] = 0;
    wntmp2D[n] = 0;
    gc[n] = 0;
    gc2D[n] = 0;
    ec[n] = 0;
    ec2D[n] = 0;
  }
  ms = 0;
  for (n = 1; n <= nwindows; n++) {
    if (window[n]->szm > ms)
      ms = window[n]->szm;
  }
  mask = i_alloc_2d(ms + 1, nwindows + 1);

  for (n = 1; n <= nwindows; n++) {
    for (c = 1; c <= ms; c++)
      mask[n][c] = 0;
    for (c = 1; c < window[n]->szc; c++) {
      ac = window[n]->wsa[c];
      g2w[n][ac] = c;
    }
  }
  for (c = 1; c <= geom->sgnum; c++)
    wm[c] = geom->fm[c].wn;

  /*-----------------------------------------------------------------*/
  /* Count the number of cells to process in each window. This       */
  /* includes the wet cells in each window and the first lateral     */
  /* auxiliary cell neighboring a wet cell. These extra cells are    */
  /* required to calculate the horizontal fluxes through faces (for  */
  /* tracers) or centers (velocity) which are subsequently used to   */
  /* update the wet cell values.                                     */
  /* Count wet cells in window n.                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    sc = geom->fm[c].sc;
    if (n > 0 && mask[n][sc] == 0) {
      wntmp[n]++;
      if (cc <= nvec2D) {
	wntmp2D[n]++;
      }
      mask[n][sc] = 1;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Count the open boundary cells (OBC) in window n.                */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      cgm = window[n]->open[nn]->no3_t;
      for (cc = 1; cc <= cgm; cc++) {
	ac = window[n]->open[nn]->obc_t[cc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          if (ac <= window[n]->enonS)
            wntmp2D[n]++;
          mask[n][ac] = 1;
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Count the first auxiliary cell in window n                      */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];                /* Global cell to process            */
    n = wm[c];                  /* Window associated with cell c     */
    sc = geom->fm[c].sc;        /* Local cell in window n            *

    /* Auxiliary and ghost cells mapped from c2c. These cells are    */
    /* required for the specification of fluxes used to subsequently */
    /* update tracers, and velocity.                                 */
    cs = geom->m2d[c];
    for (nf = 1; nf <= geom->npe[cs]; nf++) {
      cgm = geom->c2c[nf][c];
      if (wm[cgm] != n) {
	ac = window[n]->c2c[nf][sc];
	if (mask[n][ac] == 0) {
	  wntmp[n]++;
	  if (wm[cgm] == 0)  /* Ghost cells                          */
	    gc[n]++;
	  else               /* Auxiliary cells                      */
	    ec[n]++;
	  if (cc <= nvec2D) {
	    wntmp2D[n]++;
	    if (wm[cgm] == 0)
	      gc2D[n]++;
	    else {
	      ec2D[n]++;
	    }
	  }
	  mask[n][ac] = 1;
	}
      }
    }
  }

  /* Include centres surrounding each vertex of a given cell. These  */
  /* cells are required to normalize the relative vorticity computed */
  /* at vertices (dz is required at these centres).                  */
  for (ee = 1; ee <= geom->b3_e1; ee++) {
    e = geom->w3_e1[ee];
    es = geom->m2de[e];
    n = wm[geom->e2c[e][0]];
    if (!n) n = wm[geom->e2c[e][1]];
    for (i = 1; i <= geom->nee[es]; i++) {
      eoe = geom->eSe[i][e];
      for (j = 0; j <= 1; j++) {
	v = geom->e2v[eoe][j];
	vs = geom->m2dv[v];
	for (cc = 1; cc <= geom->nvc[vs]; cc++) {
	  cgm = geom->v2c[v][cc];
	  if (cgm && wm[cgm] != n) {
	    ac = g2w[n][cgm];
	    if (mask[n][ac] == 0) {
	      wntmp[n]++;
	      if (wm[cgm] == 0)   /* Ghost cells                     */
		gc[n]++;
	      else                /* Auxiliary cells                 */
		ec[n]++;
	      if (cgm == geom->m2d[cgm]) {
		wntmp2D[n]++;
		if (wm[cgm] == 0)
		  gc2D[n]++;
		else {
		  ec2D[n]++;
		}
	      }
	      mask[n][ac] = 1;
	    }
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Count the boundary ghost cells in window n                      */
  for (cc = geom->v3_t+1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];         /* Global cell to process            */
    cs = geom->m2d[c];          /* Surface cell                      */
    n = wm[c];                  /* Window associated with cell c     */
    sc = geom->fm[c].sc;        /* Local cell in window n            */

    for (nf = 1; nf <= geom->npe[cs]; nf++) {
      cgm = geom->c2c[nf][c];
      if (wm[cgm] != n) {
	ac = window[n]->c2c[nf][sc];
	if (mask[n][ac] == 0) {
	  if (wm[cgm] == 0) {
	    gc[n]++;
	    wntmp[n]++;
	    if (c == geom->m2d[c]) {
	      wntmp2D[n]++;
	      if (wm[cgm] == 0)
		gc2D[n]++;
	    }
	    mask[n][ac] = 1;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the cells to process vectors                */
  for (n = 1; n <= nwindows; n++) {
    gc[n] = wntmp[n] - (gc[n] + ec[n]);
    gc2D[n] = wntmp2D[n] - (gc2D[n] + ec2D[n]);
    ec[n] += gc[n];
    ec2D[n] += gc2D[n];

    window[n]->n3_t = wntmp[n];   /* Total 3D cells to process       */
    window[n]->n2_t = wntmp2D[n]; /* Total 2D cells to process       */
    window[n]->a3_t = ec[n];      /* 3D wet + auxiliary cells        */
    window[n]->a2_t = ec2D[n];    /* 2D wet + auxiliary cells        */
    window[n]->w3_t = i_alloc_1d(wntmp[n] + 1);
    window[n]->w2_t = i_alloc_1d(wntmp2D[n] + 1);
  }

  /*-----------------------------------------------------------------*/
  /* Fill the local window cells to process arrays                   */
  for (n = 1; n <= nwindows; n++) {
    wntmp[n] = 1;
    wntmp2D[n] = 1;
    c = ec[n];
    cc = ec2D[n];
    ec[n] = gc[n] + 1;
    ec2D[n] = gc2D[n] + 1;
    gc[n] = c + 1;
    gc2D[n] = cc + 1;
    for (c = 1; c <= ms; c++)
      mask[n][c] = 0;
  }

  /* First fill the cells to process vectors with wet cells          */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    if (n > 0) {
      if (verbose == n && cc <= nvec2D) 
	printf("wet wn=%d cc=%d c=%d(%d)\n", n, wntmp2D[n], geom->fm[c].sc, c);
      wsa_cells(window[n], &wntmp[n], &wntmp2D[n], geom->fm[c].sc, mask[n],
                mode);
    }
  }

  /* Save the number of wet cells in each vector                     */
  for (n = 1; n <= nwindows; n++) {
    window[n]->v3_t = wntmp[n] - 1;
    window[n]->v2_t = wntmp2D[n] - 1;
  }

  /* Fill the cells to process vectors with OBC cells                */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      cgm = window[n]->open[nn]->no3_t;
      for (cc = 1; cc <= cgm; cc++) {
	ac = window[n]->open[nn]->obc_t[cc];
	if (verbose == n && ac < window[n]->szcS) 
	  printf("OBC wn=%d cc=%d c=%d(%d)\n", n, wntmp2D[n], ac, window[n]->wsa[ac]);
        wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
      }
    }
  }

  /* Save the number of OBC cells in each vector                     */
  for (n = 1; n <= nwindows; n++) {
    window[n]->b3_t = wntmp[n] - 1;
    window[n]->b2_t = wntmp2D[n] - 1;
  }

  /* Fill the cells to process vectors with the first auxiliary cell */
  /* neighbouring the wet cells.                                     */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cs = geom->m2d[c];
    n = wm[c];
    ac = geom->fm[c].sc;

    /* Auxiliary and ghost cells mapped from c2c.                    */
    for (nf = 1; nf <= geom->npe[cs]; nf++) {
      cgm = geom->c2c[nf][c];
      wn = wm[cgm];
      sc = window[n]->c2c[nf][ac];
      if (wn != n) {
	if (wn == 0) {
	  if (verbose == n && !mask[n][sc] && cc <= nvec2D)
	    printf("ghost wn=%d cc=%d c=%d(%d)\n", n, gc2D[n], sc, cgm);
	  wsa_cells(window[n], &gc[n], &gc2D[n], sc, mask[n], mode);
	} else {
	  if (verbose == n && !mask[n][sc] && cc <= nvec2D)
	    printf("aux (c2c) wn=%d cc=%d c=%d(%d)\n", n, wntmp2D[n], sc, cgm);
	  wsa_cells(window[n], &wntmp[n], &wntmp2D[n], sc, mask[n], mode);
	}
      }
    }
  }

  /* Include centres surrounding each vertex of a given cell. These  */
  /* cells are required to normalize the relative vorticity computed */
  /* at vertices (dz is required at these centres).                  */
  for (ee = 1; ee <= geom->b3_e1; ee++) {
    e = geom->w3_e1[ee];
    es = geom->m2de[e];
    n = wm[geom->e2c[e][0]];
    if (!n) n = wm[geom->e2c[e][1]];
    for (i = 1; i <= geom->nee[es]; i++) {
      eoe = geom->eSe[i][e];
      for (j = 0; j <= 1; j++) {
	v = geom->e2v[eoe][j];
	vs = geom->m2dv[v];
	for (cc = 1; cc <= geom->nvc[vs]; cc++) {
	  cgm = geom->v2c[v][cc];
	  if (cgm && wm[cgm] != n) {
	    sc = g2w[n][cgm];
	    if (wm[cgm] == 0) {
	      if (verbose == n && !mask[n][sc] && e < geom->szeS)
		printf("ghost wn=%d cc=%d c=%d(%d)\n", n, gc2D[n], sc, cgm);
	      wsa_cells(window[n], &gc[n], &gc2D[n], sc, mask[n], mode);
	    } else {
	      if (verbose == n && !mask[n][sc] && e < geom->szeS)
		printf("aux (vertex) wn=%d cc=%d c=%d(%d)\n", n, wntmp2D[n], sc, cgm);
	      wsa_cells(window[n], &wntmp[n], &wntmp2D[n], sc, mask[n], mode);
	    }
	  }
	}
      }
    }
  }

  /*------------------------------------------------------------------*/
  /* Fill with the boundary ghost cells in window n                   */
  for (cc = geom->v3_t+1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];         /* Global cell to process             */
    cs = geom->m2d[c];          /* Surface cell                       */
    n = wm[c];                  /* Window associated with cell c      */
    sc = geom->fm[c].sc;        /* Local cell in window n             */

    for (nf = 1; nf <= geom->npe[cs]; nf++) {
      cgm = geom->c2c[nf][c];
      if (wm[cgm] != n) {
	ac = window[n]->c2c[nf][sc];
	if (mask[n][ac] == 0) {
	  if (wm[cgm] == 0)
	    wsa_cells(window[n], &gc[n], &gc2D[n], ac, mask[n], mode);
	}
      }
    }
  }
  i_free_1d(wntmp);
  i_free_1d(wntmp2D);
  i_free_1d(gc);
  i_free_1d(gc2D);
  i_free_1d(ec);
  i_free_1d(ec2D);
  i_free_1d(wm);
  i_free_2d(g2w);
  i_free_2d(mask);

  /*------------------------------------------------------------------*/
  /* Check that all cells are accounted for                           */
  if (checkf) {
    int ii, nc;
    int sv = 1;
    int ev = geom->v3_t;
    wm = i_alloc_1d(geom->szc);
    memset(wm, 0, geom->szc * sizeof(int));
    for (n = 1; n <= geom->nwindows; n++) {
      for (cc = sv; cc <= ev; cc++) {
	c = geom->w3_t[cc];
	for (ii = 1; ii <= window[n]->v3_t; ii++) {
	  sc = window[n]->w3_t[ii];
	  if (sc >= window[n]->szc) printf("Error in w3_t : lc=%d, size=%d\n", sc, window[n]->szc);
	  if (c == window[n]->wsa[sc]) {
	    wm[c] = n;
	    break;
	  }
	}
      }
      printf("Cell window %d\n", n);
      printf("  wet = %d(2D) %d(3D)\n", window[n]->v2_t, window[n]->v3_t);
      printf("  OBC = %d(2D) %d(3D)\n", window[n]->b2_t, window[n]->b3_t);
      printf("  aux = %d(2D) %d(3D)\n", window[n]->a2_t, window[n]->a3_t);
      printf("  ghost = %d(2D) %d(3D)\n", window[n]->n2_t, window[n]->n3_t);
      printf("  size = %d(2D) %d(3D)\n", window[n]->szcS, window[n]->szc);
    }
    printf("Global\n");
    printf("  wet = %d(2D) %d(3D)\n", geom->v2_t, geom->v3_t);
    printf("  OBC = %d(2D) %d(3D)\n", geom->b2_t, geom->b3_t);
    printf("  ghost = %d(2D) %d(3D)\n", geom->n2_t, geom->n3_t);
    printf("  size = %d(2D) %d(3D)\n", geom->szcS, geom->szc);
    nc = 0;  
    for (cc = sv; cc <= ev; cc++) {
      c = geom->w3_t[cc];
      if(!wm[c]) {
	nc++;
	printf("Centre %d c not found at cc=%d cg=%d (%d %d %d)\n",nc,cc,c,geom->s2i[c],geom->s2j[c],geom->s2k[c]);
      }
    }
    i_free_1d(wm);

    if (geom->nwindows == 1) {
      if (geom->v3_t != window[1]->v3_t)
	printf("Error in v3_t : %d vs %d\n",geom->v3_t, window[1]->v3_t);
      if (geom->n3_t != window[1]->n3_t)
	printf("Error in n3_t : %d vs %d\n",geom->n3_t, window[1]->n3_t);
      for (cc = 1; cc <= geom->n3_t; cc++) {
	c = geom->w3_t[cc];
	sc = window[1]->w3_t[cc];
	if (c != sc)
	  printf("Error in w3_t : cc=%d gc=%d lc=%d\n",cc, c, sc);
      }
    }
  }
}

/* END get_local_wsc()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the local sparse coordinates for the work      */
/* array of edges to process from the global sparse work array.      */
/*-------------------------------------------------------------------*/
void get_local_wse(geometry_t *geom,
		   int *vec,            /* Global work sparse array  */
                   int nvec,            /* Size of vec               */
                   int nvec2D,          /* Last 2D cell loc in vec   */
                   geometry_t **window, /* Window data structure     */
                   int nwindows,        /* Number of windows         */
		   double smag          /* Smagorinsky flag          */
		   )
{
  int *wntmp;                   /* Temporary buffer                  */
  int *wntmp2D;                 /* Temporary buffer                  */
  int c, cs, cc, lc;            /* Centre counters                   */
  int ee, e;                    /* Edge counters                     */
  int n, wn;                    /* Window counters                   */
  int ac, le, ge, wae;          /* Mapped centres, edges             */
  int vv, v, vs, gv, lv;        /* Mapped vertices                   */
  int **mask;                   /* Set to 1 when cell is processed   */
  int *wm;                      /* Window / ghost cell mask          */
  int *msk;                     /* Edge status mask                  */
  int sc;                       /* Local wet cell                    */
  int *gc, *gc2D;               /* Ghost edge                        */
  int *ec, *ec2D;               /* Auxiliary edge                    */
  int *bc, *bc2D;               /* OBC edge                          */
  int i, ii, j, jj;             /* Counters                          */
  int checkf = 0;
  int verbose = 0;
  int dof = 0;

  /*-----------------------------------------------------------------*/
  /* Initialise the counters                                         */
  wntmp = i_alloc_1d(nwindows + 1);
  wntmp2D = i_alloc_1d(nwindows + 1);
  gc = i_alloc_1d(nwindows + 1);
  gc2D = i_alloc_1d(nwindows + 1);
  ec = i_alloc_1d(nwindows + 1);
  ec2D = i_alloc_1d(nwindows + 1);
  bc = i_alloc_1d(nwindows + 1);
  bc2D = i_alloc_1d(nwindows + 1);
  wm = i_alloc_1d(geom->sze);
  msk = i_alloc_1d(geom->sze);
  for (n = 1; n <= nwindows; n++) {
    wntmp[n] = 0;
    wntmp2D[n] = 0;
    gc[n] = 0;
    gc2D[n] = 0;
    ec[n] = 0;
    ec2D[n] = 0;
    bc[n] = 0;
    bc2D[n] = 0;
  }
  mask = i_alloc_2d(geom->sze, nwindows + 1);
  for (n = 1; n <= nwindows; n++) {
    for (e = 1; e < geom->sze; e++)
      mask[n][e] = 0;
    window[n]->v2_e1 = window[n]->b2_e1 = window[n]->a2_e1 = window[n]->n2_e1 = window[n]->x2_e1 = 0;
    window[n]->v3_e1 = window[n]->b3_e1 = window[n]->a3_e1 = window[n]->n3_e1 = window[n]->x3_e1 = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the window map for edges. Here the edges are defined for    */
  /* each window; edges are defined the same a window number as the  */
  /* centre in the direction e2c[e][0], or equivalently those edges  */
  /* where eSc = -1.                                                 */
  for (ee = 1; ee <= geom->v3_e1; ee++) {
    e = geom->w3_e1[ee];
    c = geom->e2c[e][0];
    wm[e] = geom->fm[c].wn;
    msk[e] = E_W;
  }
  for (ee = geom->v3_e1 + 1; ee <= geom->b3_e1; ee++) {
    e = geom->w3_e1[ee];
    c = geom->e2c[e][0];
    wm[e] = geom->fm[c].wn;
    /* If c is a ghost cell, use the other centre adjacent to e      */
    if (!wm[e]) wm[e] = geom->fm[geom->e2c[e][1]].wn;
    msk[e] = E_B;
  }
  for (ee = geom->b3_e1 + 1; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    c = geom->e2c[e][0];
    wm[e] = geom->fm[c].wn;
    if (!wm[e]) wm[e] = geom->fm[geom->e2c[e][1]].wn;
    msk[e] = E_G;
  }

  /*-----------------------------------------------------------------*/
  /* Count the edges in each window.                                 */
  /* Edges are required for all those edges connecting each vertex   */
  /* of cells that are adjacent to a wet edge in the window.         */
  /* First ghost and auxiliary cells surrounding cell centres to     */
  /* compute divergence. These will not be captured in the loop      */
  /* below over stepped bathymetry.                                  */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];            /* Global cell centre             */
    cs = geom->m2d[c];             /* Surface global cell centre     */
    wn = geom->fm[c].wn;           /* Window corresponding to c      */
    sc = geom->fm[c].sc;           /* Local cell centre in window wn */
    for (j = 1; j <= geom->npe[cs]; j++) {
      e = geom->c2e[j][c];         /* Global edge surrounding c      */
      le = window[wn]->c2e[j][sc]; /* Local edge surrounding sc      */
      if (le) {
	wae = wm[e];               /* Window corresponding to e      */
	if (msk[e] & E_G || wae != wn) {
	  if (mask[wn][le] == 0 && le != sc) {
	    if (msk[e] & E_G) {
	      window[wn]->n3_e1++;
	      if (c == cs) {
		window[wn]->n2_e1++;
	      }
	    } else {
	      window[wn]->a3_e1++;
	      if (c == cs) {
		window[wn]->a2_e1++;
	      }
	    }
	    mask[wn][le] = 1;
	  }
	}
      }
    }
  }

  /* Next ghost and auxiliary edges to compute vorticity             */
  for (ee = 1; ee <= nvec; ee++) {
    int gc;
    e = vec[ee];                /* Global edge to process            */
    wn = wm[e];                 /* Window associated with edge e     */
    sc = geom->g2we[wn][e];     /* Local edge in window wn           */
    for (i = 0; i <= 1; i++) {  /* Loop over centres adjacent to e   */
      c = geom->e2c[e][i];      /* Cell centre in e2c[i] direction   */
      cs = geom->m2d[c];        /* Surface cell                      */
      lc = window[wn]->e2c[sc][i];
      for (j = 1; j <= geom->npe[cs]; j++) {
	gv = geom->c2v[j][c];           /* Global vertex of c        */
	vs = geom->m2dv[gv];            /* Surface vertex of gv      */
	lv = window[wn]->c2v[j][lc];    /* Local vertex of c         */
	for (ii = 1; ii <= geom->nve[vs]; ii++) {
	  ge = geom->v2e[gv][ii];      /* Global edge connected to v */
	  le = window[wn]->v2e[lv][ii]; /* Local edge                */
	  if (le) {
	    wae = wm[ge];              /* Window corresponding to ge */
	    if (msk[ge] & E_G || wae != wn) {
	      if (mask[wn][le] == 0 && le != sc) {
		if (msk[ge] & E_G) {
		  window[wn]->n3_e1++;
		  if (ee <= nvec2D) {
		    window[wn]->n2_e1++;
		  }
		} else {
		  window[wn]->a3_e1++;
		  if (ee <= nvec2D) {
		    window[wn]->a2_e1++;
		  }
		}
		mask[wn][le] = 1;
	      }
	    }
	  }
	}
      }
    }
  }

  /* Wet edges                                                       */
  for (ee = 1; ee <= geom->v3_e1; ee++) {
    e = geom->w3_e1[ee];
    wn = wm[e];
    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    if (msk[e] & E_W && !mask[wn][sc]) {
      window[wn]->v3_e1++;
      if (e < geom->szeS) window[wn]->v2_e1++;
      mask[wn][sc] = 1;
    }
  }
  /* Boundary edges                                                  */
  for (ee = geom->v3_e1 + 1; ee <= geom->b3_e1; ee++) {
    e = geom->w3_e1[ee];
    wn = wm[e];
    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    if (msk[e] & E_B && !mask[wn][sc]) {
      window[wn]->b3_e1++;
      if (e < geom->szeS) window[wn]->b2_e1++;
      mask[wn][sc] = 1;
    }
  }

  /* Tangential ghost edges.                                         */
  /* These are currently not used, but are included to make the size */
  /* of the vector the same as the master.                           */
  for (ee = geom->nbpte1 + 1; ee <= geom->nbe1; ee++) {
    e = geom->bpte1[ee];

    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    c = geom->bine1[ee];
    n = geom->fm[c].wn;            /* Window associated with edge e1 */
    if (mask[n][sc] == 0) {
      window[wn]->x3_e1++;
      if (c == geom->m2d[c]) {
	window[wn]->x2_e1++;
      }
      mask[n][sc] = 1;
    }
  }

  /* Get edges used in Smagorinsky for multiple windows missed       */
  /* above and put in auxiliary vectors.                             */
  if (smag != 0.0) {
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];       /* Global cell centre                */
      cs = geom->m2d[c];        /* Surface cell                      */
      wn = geom->fm[c].wn;      /* Window associated with cell c     */ 
      lc = geom->fm[c].sc;      /* Local coordinate                  */
      for (j = 1; j <= geom->npe[cs]; j++) {
	int cn, cns, scn;
	cn = geom->c2c[j][c];   /* Neighbour of c                    */
	cns = geom->m2d[cn];    /* Surface neighbour                 */
	if (!geom->wgst[cn] && wn != geom->fm[cn].wn) {
	  scn = window[wn]->c2c[j][lc]; /* Local neighbour           */

	  for (ee = 1; ee <= geom->npe[cns]; ee++) {
	    e = geom->c2e[ee][cn];         /* Global edge to cn      */
	    if (!e) continue;
	    le = window[wn]->c2e[ee][scn]; /* Local edge to scn      */
	    if (mask[wn][le] == 0) {
	      if (geom->wgst[cn] || msk[e] & E_G) {
		window[wn]->n3_e1++;
		if (c == cs) {
		  window[wn]->n2_e1++;
		}
	      } else {
		window[wn]->a3_e1++;
		if (c == cs) {
		  window[wn]->a2_e1++;
		}
	      }
	      mask[wn][le] = 1;
	    }
	  }

	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the edges to process vectors                */
  for (n = 1; n <= nwindows; n++) {
    window[n]->b3_e1 += window[n]->v3_e1;
    window[n]->a3_e1 += window[n]->b3_e1;
    window[n]->n3_e1 += window[n]->a3_e1;
    window[n]->x3_e1 += window[n]->n3_e1;
    window[n]->b2_e1 += window[n]->v2_e1;
    window[n]->a2_e1 += window[n]->b2_e1;
    window[n]->n2_e1 += window[n]->a2_e1;
    window[n]->x2_e1 += window[n]->n2_e1;
    /*
    if (nwindows == 1) {
      window[n]->sze = window[n]->x3_e1 + window[n]->x2_e1 + 1;
      window[n]->szeS = window[n]->x2_e1 + 1;
    }
    */
    window[n]->w3_e1 = i_alloc_1d(window[n]->sze);
    window[n]->w2_e1 = i_alloc_1d(window[n]->szeS);
    wntmp[n] = 1;
    wntmp2D[n] = 1;
    bc[n] = window[n]->v3_e1 + 1;
    bc2D[n] = window[n]->v2_e1 + 1;
    ec[n] = window[n]->b3_e1 + 1;
    ec2D[n] = window[n]->b2_e1 + 1;
    gc[n] = window[n]->a3_e1 + 1;
    gc2D[n] = window[n]->a2_e1 + 1;
  }
  for (n = 1; n <= nwindows; n++) {
    for (e = 1; e < geom->sze; e++)
      mask[n][e] = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Fill the edges to process vectors.                              */
  /* Ghost and auxiliary edges.                                      */
  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    wn = geom->fm[c].wn;
    sc = geom->fm[c].sc;
    for (j = 1; j <= geom->npe[cs]; j++) {
      e = geom->c2e[j][c];
      le = window[wn]->c2e[j][sc];    /* Local edge of c         */
      if (le) {
	wae = wm[e];              /* Window corresponding to e */

	if (msk[e] & E_G || wae != wn) {
	  if (mask[wn][le] == 0 && le != sc) {
	    if (msk[e] & E_G) {
	      window[wn]->w3_e1[gc[wn]] = le;
	      gc[wn]++;
	      if (c == cs) {
		window[wn]->w2_e1[gc2D[wn]] = le;
		if (verbose) printf("ghost wn=%d ee=%d e=%d(%d)\n",wn, gc2D[wn], le, e);
		gc2D[wn]++;
	      }
	    } else {
	      window[wn]->w3_e1[ec[wn]] = le;
	      ec[wn]++;
	      if (c == cs) {
		window[wn]->w2_e1[ec2D[wn]] = le;
		if (verbose) printf("aux wn=%d ee=%d e=%d(%d)\n",wn, ec2D[wn], le, e);
		ec2D[wn]++;
	      }
	    }
	    mask[wn][le] = 1;
	  }
	}
      }
    }
  }

  for (ee = 1; ee <= nvec; ee++) {
    e = vec[ee];                /* Global edge to process            */
    wn = wm[e];                 /* Window associated with edge e     */
    sc = geom->g2we[wn][e];     /* Local edge in window wn           */

    for (i = 0; i <= 1; i++) {  /* Loop over centres adjacent to e   */
      c = geom->e2c[e][i];      /* Cell centre in e2c[i] direction   */
      cs = geom->m2d[c];        /* Surface cell                      */
      lc = window[wn]->e2c[sc][i];
      for (j = 1; j <= geom->npe[cs]; j++) {
	gv = geom->c2v[j][c];           /* Global vertex of c        */
	vs = geom->m2dv[gv];            /* Surface vertex of gv      */
	lv = window[wn]->c2v[j][lc];    /* Local vertex of c         */
	for (ii = 1; ii <= geom->nve[vs]; ii++) {
	  ge = geom->v2e[gv][ii];      /* Global edge connected to v */
	  le = window[wn]->v2e[lv][ii]; /* Local edge                */
	  if (le) {
	    wae = wm[ge];              /* Window corresponding to ge */
	    if (msk[ge] & E_G || wae != wn) {
	      if (mask[wn][le] == 0 && le != sc) {
		if (msk[ge] & E_G) {
		  window[wn]->w3_e1[gc[wn]] = le;
		  gc[wn]++;
		  if (ee <= nvec2D) {
		    window[wn]->w2_e1[gc2D[wn]] = le;
		    if (verbose) printf("ghost wn=%d ee=%d e=%d(%d)\n",wn, gc2D[wn], le, ge);
		    gc2D[wn]++;
		  }
		} else {
		  window[wn]->w3_e1[ec[wn]] = le;
		  ec[wn]++;

		  if (ee <= nvec2D) {
		    window[wn]->w2_e1[ec2D[wn]] = le;
		    if (verbose) printf("aux wn=%d ee=%d e=%d(%d)\n",wn, ec2D[wn], le, ge);
		    ec2D[wn]++;
		  }
		}
		mask[wn][le] = 1;
	      }
	    }
	  }
	}
      }
    }
  }

  /* Wet edges                                                       */
  for (ee = 1; ee <= geom->v3_e1; ee++) {
    e = geom->w3_e1[ee];
    wn = wm[e];
    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    if (msk[e] & E_W && !mask[wn][sc]) {
      window[wn]->w3_e1[wntmp[wn]] = sc;
      wntmp[wn]++;
      if (e < geom->szeS) {
	window[wn]->w2_e1[wntmp2D[wn]] = sc;
	if (verbose) printf("wet wn=%d ee=%d e=%d(%d)\n",wn, wntmp2D[wn], sc, e);
	wntmp2D[wn]++;
      }
      mask[wn][sc] = 1;
    }
  }

  /* Boundary edges                                                  */
  for (ee = geom->v3_e1 + 1; ee <= geom->b3_e1; ee++) {
    e = geom->w3_e1[ee];
    wn = wm[e];
    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    if (msk[e] & E_B && !mask[wn][sc]) {
      window[wn]->w3_e1[bc[wn]] = sc;
      bc[wn]++;
      if (e < geom->szeS) {
	window[wn]->w2_e1[bc2D[wn]] = sc;
	if (verbose) printf("OBC wn=%d ee=%d e=%d(%d)\n",wn, bc2D[wn], sc, e);
	bc2D[wn]++;
      }
      mask[wn][sc] = 1;
    }
  }

  /* Tangential ghost edges.                                         */
  /* These are currently not used, but are included to make the size */
  /* of the vector the same as the master.                           */

  for (ee = geom->nbpte1 + 1; ee <= geom->nbe1; ee++) {
    e = geom->bpte1[ee];

    /*sc = geom->fm[e].ec;*/
    sc = geom->g2we[wn][e];
    c = geom->bine1[ee];
    n = geom->fm[c].wn;            /* Window associated with edge e1 */
    if (mask[n][sc] == 0) {
      mask[n][sc] = 1;
    }
  }

  /*
  if (nwindows == 1) printf("aa %d\n", window[1]->c2e[1][1245]);
  if (nwindows == 5) printf("aa %d\n", window[4]->c2e[1][350]);
  */
  /* Get edges used in Smagorinsky missed above                      */
  if (smag != 0.0) {
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];       /* Global cell centre                */
      cs = geom->m2d[c];        /* Surface cell                      */
      wn = geom->fm[c].wn;      /* Window associated with cell c     */ 
      lc = geom->fm[c].sc;      /* Local coordinate                  */
      for (j = 1; j <= geom->npe[cs]; j++) {
	int cn, cns, scn;
	cn = geom->c2c[j][c];   /* Neighbour of c                    */
	cns = geom->m2d[cn];    /* Surface neighbour                 */
	if (!geom->wgst[cn] && wn != geom->fm[cn].wn) {
	  scn = window[wn]->c2c[j][lc]; /* Local neighbour           */

	  for (ee = 1; ee <= geom->npe[cns]; ee++) {
	    e = geom->c2e[ee][cn];         /* Global edge to cn      */
	    if (!e) continue;
	    le = window[wn]->c2e[ee][scn]; /* Local edge to scn      */
	    if (mask[wn][le] == 0) {
	      if (geom->wgst[cn] || msk[e] & E_G) {
		window[wn]->w3_e1[gc[wn]] = le;
		gc[wn]++;
		if (c == cs) {
		  window[wn]->w2_e1[gc2D[wn]] = le;
		  if (verbose) printf("Smag ghost wn=%d ee=%d e=%d(%d)\n",wn, gc2D[wn], le, e);
		  gc2D[wn]++;
		  window[wn]->n2_e1++;
		}
	      } else {
		window[wn]->w3_e1[ec[wn]] = le;
		ec[wn]++;
		if (c == cs) {
		  window[wn]->w2_e1[ec2D[wn]] = le;
		  if (verbose) printf("Smag aux wn=%d ee=%d e=%d(%d)\n",wn, ec2D[wn], le, e);
		  ec2D[wn]++;
		}
	      }
	      mask[wn][le] = 1;
	    }
	  }
	}
      }
    }
  }

  for (ee = 1; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    geom->fm[e].we = wm[e];
  }

  i_free_1d(wntmp);
  i_free_1d(wntmp2D);
  i_free_1d(gc);
  i_free_1d(gc2D);
  i_free_1d(ec);
  i_free_1d(ec2D);
  i_free_1d(bc);
  i_free_1d(bc2D);
  i_free_1d(wm);
  i_free_1d(msk);
  i_free_2d(mask);

  /* Check that all edges are accounted for */
  if (checkf) {
    int ii;
    wm = i_alloc_1d(geom->sze);
    memset(wm, 0, geom->sze * sizeof(int));
    for (n = 1; n <= geom->nwindows; n++) {
      for (ee = 1; ee <= geom->n3_e1; ee++) {
	e = geom->w3_e1[ee];
	if (geom->nwindows == 1 && e != window[1]->w3_e1[ee])
	  printf("Edge different at %d %d : %d\n",ee, e, window[1]->w3_e1[ee]);
	for (ii = 1; ii <= window[n]->n3_e1; ii++) {
	  sc = window[n]->w3_e1[ii];
	  if (sc >= window[n]->sze) printf("Error in w3_e1 : le=%d, size=%d\n", sc, window[n]->sze);
	  if (e == window[n]->wse[sc]) {
	    wm[e] = 1;
	    break;
	  }
	}
      }
      printf("Edge window %d\n", n);
      printf("  wet = %d(2D) %d(3D)\n", window[n]->v2_e1, window[n]->v3_e1);
      printf("  OBC = %d(2D) %d(3D)\n", window[n]->b2_e1, window[n]->b3_e1);
      printf("  aux = %d(2D) %d(3D)\n", window[n]->a2_e1, window[n]->a3_e1);
      printf("  ghost = %d(2D) %d(3D)\n", window[n]->n2_e1, window[n]->n3_e1);
      printf("  size = %d(2D) %d(3D)\n", window[n]->szeS, window[n]->sze);
    }
    printf("Global\n");
    printf("  wet = %d(2D) %d(3D)\n", geom->v2_e1, geom->v3_e1);
    printf("  OBC = %d(2D) %d(3D)\n", geom->b2_e1, geom->b3_e1);
    printf("  aux = %d(2D) %d(3D)\n", geom->a2_e1, geom->a3_e1);
    printf("  ghost = %d(2D) %d(3D)\n", geom->n2_e1, geom->n3_e1);
    printf("  size = %d(2D) %d(3D)\n", geom->szeS, geom->sze);
    cc = 0;  
    for (ee = 1; ee <= geom->n3_e1; ee++) {
      e = geom->w3_e1[ee];
      if(!wm[e]) {
	cc++;
	c = geom->e2c[e][0];
	if (geom->wgst[c])
	  c = geom->e2c[e][1];
	printf("Edge %d e not found at ee=%d e=%d (%d %d %d)\n",cc,ee,e,geom->s2i[c],geom->s2j[c],geom->s2k[c]);
      }
    }
    i_free_1d(wm);

    if (geom->nwindows == 1) {
      if (geom->v3_e1 != window[1]->v3_e1)
	printf("Error in v3_e1 : %d vs %d\n",geom->v3_e1, window[1]->v3_e1);
      if (geom->n3_e1 != window[1]->n3_e1)
	printf("Error in n3_e1 : %d vs %d\n",geom->n3_e1, window[1]->n3_e1);
      for (cc = 1; cc <= geom->v3_e1; cc++) {
	c = geom->w3_e1[cc];
	sc = window[1]->w3_e1[cc];
	if (c != sc)
	  printf("Error in w3_e1 : ee=%d ge=%d le=%d\n",cc, c, sc);
      }
    }
  }
}

/* END get_local_wse()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to extract the local sparse coordinates for the work      */
/* array of vertices to process from the global sparse work array.   */
/* Note: the same vertex may be included in multiple windows. To     */
/* find vertices associated with windows, loop over all wet and      */
/* boundary centres, and each vertex of the centre is assigned to    */
/* the window of that centre.                                        */
/*-------------------------------------------------------------------*/
void get_local_wsv(geometry_t *geom,
		   int *vec,            /* Global work sparse array  */
                   int nvec,            /* Size of vec               */
                   int nvec2D,          /* Last 2D cell loc in vec   */
                   geometry_t **window, /* Window data structure     */
                   int nwindows         /* Number of windows         */
		   )
{
  int c, cs, cc, n, wn, j;      /* Counters                          */
  int vv, sv, v, vs;            /* Vertices                          */
  int **mask;                   /* Set to 1 when cell is processed   */
  int *wm;                      /* Window / ghost cell mask          */
  int *msk;                     /* Vertex status mask                */
  int sc;                       /* Local eet cell in window          */
  int checkf = 0;
  /*int nwindows = geom->nwindows;*/

  wm = i_alloc_1d(geom->szc);
  msk = i_alloc_1d(geom->szv);
  mask = i_alloc_2d(geom->szv, nwindows + 1);
  for (n = 1; n <= nwindows; n++) {
    for (v = 1; v < geom->szv; v++)
      mask[n][v] = 0;
    window[n]->v2_e2 = window[n]->b2_e2 = window[n]->a2_e2 = window[n]->n2_e2 = 0;
    window[n]->v3_e2 = window[n]->b3_e2 = window[n]->a3_e2 = window[n]->n3_e2 = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the window map for edges                                    */
  for (cc = 1; cc <= geom->n3_t; cc++) {
    c = geom->w3_t[cc];
    wm[c] = geom->fm[c].wn;
  }
  for (vv = 1; vv <= geom->v3_e2; vv++) {
    v = geom->w3_e2[vv];
    msk[v] = V_W;
  }
  for (vv = geom->v3_e2 + 1; vv <= geom->b3_e2; vv++) {
    v = geom->w3_e2[vv];
    msk[v] = V_B;
  }
  for (vv = geom->b3_e2 + 1; vv <= geom->n3_e2; vv++) {
    v = geom->w3_e2[vv];
    msk[v] = V_G;
  }

  /*-----------------------------------------------------------------*/
  /* Count the vertices in each window.                              */
  for (wn = 1; wn <= nwindows; wn++) {
    for (cc = 1; cc <= window[wn]->enon; cc++) {
      c = window[wn]->wsa[cc];
      cs = geom->m2d[c];
      /* Get the vertices in the window. Vertices of every cell      */
      /* centre in the window are included.                          */

      for (j = 1; j <= geom->npe[cs]; j++) {
	int ic;
	v = geom->c2v[j][c];
	vs = geom->m2dv[v];

	if (v) {
	  sv = geom->g2wv[wn][v];

	  /* Cells within window wn                                  */
	  if (wn == wm[c]) {
	    /* Ghost cells                                           */
	    if(msk[v] & V_G && !mask[wn][sv]) {
	      window[wn]->n3_e2++;
	      if (v < geom->szvS) window[wn]->n2_e2++;
	      mask[wn][sv] = 1;
	    }
	    /* Wet cells                                             */
	    if(msk[v] & V_W && !mask[wn][sv]) {
	      window[wn]->v3_e2++;
	      if (v < geom->szvS) window[wn]->v2_e2++;
	      mask[wn][sv] = 1;
	    }
	    /* OBC cells                                             */
	    if(msk[v] & V_B && !mask[wn][sv]) {
	      window[wn]->b3_e2++;
	      if (v < geom->szvS) window[wn]->b2_e2++;
	      mask[wn][sv] = 1;
	    }
	  } else {
	    /* Auxiliary cells                                       */
	    /* Ghost cells                                           */
	    if(msk[v] & V_G && !mask[wn][sv]) {
	      window[wn]->n3_e2++;
	      if (v < geom->szvS) window[wn]->n2_e2++;
	      mask[wn][sv] = 1;
	    }
	    /* Wet / OBC cells                                       */
	    if(msk[v] & (V_W|V_B|V_G) && !mask[wn][sv]) {
	      window[wn]->a3_e2++;
	      if (v < geom->szvS) window[wn]->a2_e2++;
	      mask[wn][sv] = 1;
	    }
	  }
	}
      }
    }
  }	    

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the vertices to process vectors             */
  for (n = 1; n <= nwindows; n++) {
    if (nwindows == 1) {
      window[n]->szv = window[n]->n3_e2 + window[n]->a3_e2 + window[n]->b3_e2 + window[n]->v3_e2 + 1;
      window[n]->szvS = window[n]->n2_e2 + window[n]->a2_e2 + window[n]->b2_e2 + window[n]->v2_e2 + 1;
    }
    window[n]->n3_e2 = window[n]->a3_e2 + window[n]->b3_e2 + window[n]->v3_e2 + 1;
    window[n]->a3_e2 = window[n]->b3_e2 + window[n]->v3_e2 + 1;
    window[n]->b3_e2 = window[n]->v3_e2 + 1;
    window[n]->v3_e2 = 1;
    window[n]->n2_e2 = window[n]->a2_e2 + window[n]->b2_e2 + window[n]->v2_e2 + 1;
    window[n]->a2_e2 = window[n]->b2_e2 + window[n]->v2_e2 + 1;
    window[n]->b2_e2 = window[n]->v2_e2 + 1;
    window[n]->v2_e2 = 1;
    window[n]->w3_e2 = i_alloc_1d(window[n]->szv);
    window[n]->w2_e2 = i_alloc_1d(window[n]->szvS);
  }
  for (n = 1; n <= nwindows; n++)
    for (v = 1; v < geom->szv; v++)
      mask[n][v] = 0;

  /*-----------------------------------------------------------------*/
  /* Fill the vertices to process vectors.                           */
  for (wn = 1; wn <= nwindows; wn++) {
    for (cc = 1; cc <= window[wn]->enon; cc++) {
      c = window[wn]->wsa[cc];
      cs = geom->m2d[c];
      /* Get the vertices in the window. Vertices of every cell      */
      /* centre in the window are included.                          */

      for (j = 1; j <= geom->npe[cs]; j++) {
	int ic;
	v = geom->c2v[j][c];
	vs = geom->m2dv[v];

	if (v) {
	  sv = geom->g2wv[wn][v];

	  /* Cells within window wn                                  */
	  if (wn == wm[c]) {
	    /* Ghost cells                                           */
	    if(msk[v] & V_G && !mask[wn][sv]) {
	      window[wn]->w3_e2[window[wn]->n3_e2++] = sv;
	      if (v < geom->szvS) window[wn]->w2_e2[window[wn]->n2_e2++] = sv;
	      mask[wn][sv] = 1;
	    }
	    /* Wet cells                                             */
	    if(msk[v] & V_W && !mask[wn][sv]) {
	       window[wn]->w3_e2[window[wn]->v3_e2++] = sv;
	      if (v < geom->szvS)  window[wn]->w2_e2[window[wn]->v2_e2++] = sv;
	      mask[wn][sv] = 1;
	    }
	    /* OBC cells                                             */
	    if(msk[v] & V_B && !mask[wn][sv]) {
	       window[wn]->w3_e2[window[wn]->b3_e2++] = sv;
	      if (v < geom->szvS)  window[wn]->w2_e2[window[wn]->b2_e2++] = sv;
	      mask[wn][sv] = 1;
	    }
	  } else {
	    /* Auxiliary cells                                       */
	    /* Ghost cells                                           */
	    if(msk[v] & V_G && !mask[wn][sv]) {
	      window[wn]->w3_e2[window[wn]->n3_e2++] = sv;
	      if (v < geom->szvS)  window[wn]->w2_e2[window[wn]->n2_e2++] = sv;
	      mask[wn][sv] = 1;
	    }
	    /* Wet / OBC cells                                       */
	    if(msk[v] & (V_W|V_B|V_G) && !mask[wn][sv]) {
	      window[wn]->w3_e2[window[wn]->a3_e2++] = sv;
	      if (v < geom->szvS)  window[wn]->w2_e2[window[wn]->a2_e2++] = sv;
	      mask[wn][sv] = 1;
	    }
	  }
	}
      }
    }
    window[wn]->n3_e2--;
    window[wn]->a3_e2--;
    window[wn]->b3_e2--;
    window[wn]->v3_e2--;
    window[wn]->n2_e2--;
    window[wn]->a2_e2--;
    window[wn]->b2_e2--;
    window[wn]->v2_e2--;
    if (checkf) {
      printf("win%d szv=%d n3_e2=%d\n",window[wn]->szv, window[wn]->n3_e2);
      printf("szvS=%d n2_e2=%d\n",window[wn]->szvS, window[wn]->n2_e2);
    }
  }
}

/* END get_local_wsv()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Fills the cells to process vectors                                */
/*-------------------------------------------------------------------*/
void wsa_cells(geometry_t *window,  /* Window data structure         */
               int *loc3,   /* Location of local 3D cell to process  */
               int *loc2,   /* Location of local 2D cell to process  */
               int c,       /* Local sparse cell to process          */
               int *mask,   /* Mask for cells already processed      */
               int mode     /* mode=0 : cell center cells to process */
	       )
                          /* mode=1 : e1 face cells to process       */
                          /* mode=2 : e2 face cells to process       */
{

  if (mask[c])
    return;
  mask[c] = 1;

  if (mode == 0) {
    window->w3_t[*loc3] = c;
    *loc3 += 1;
    if (c <= window->enonS) {
      window->w2_t[*loc2] = c;
      *loc2 += 1;
    }
  }
  if (mode == 1) {
    window->w3_e1[*loc3] = c;
    *loc3 += 1;
    if (c < window->szeS) {
      window->w2_e1[*loc2] = c;
      *loc2 += 1;
    }
  }
}

/* END wsa_cells()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocate the open boundary ghost cells at R_EDGE and F_EDGE if    */
/* the maps have already been read in from file.                     */
/*-------------------------------------------------------------------*/
void get_local_obc_a(geometry_t *geom,    /* Global geometry         */
		     int *vec,            /* Global work array       */
		     int nvec,            /* Size of vec             */
		     int nvec2D,          /* Last 2D cell loc in vec */
		     geometry_t **window, /* Window data structure   */
		     int nwindows         /* Number of windows       */
		     )
{
  int c, cc, cs, n, nn, wn, j;  /* Counters                          */
  int ac, cgm;                  /* Sparse coordinates                */
  int **mask;                   /* Set to 1 when cell is processed   */
  int *wm;                      /* Window / ghost cell mask          */
  int ms;                       /* Size of mask                      */
  int sc;                       /* Wet cell in window                */

  /*-----------------------------------------------------------------*/
  /* Initialise the counters */
  ms = 0;
  for (n = 1; n <= nwindows; n++) {
    if (window[n]->enon > ms)
      ms = window[n]->enon;
  }
  mask = i_alloc_2d(ms + 1, nwindows + 1);
  for (n = 1; n <= nwindows; n++)
    for (c = 1; c <= ms; c++)
      mask[n][c] = 0;
  wm = i_alloc_1d(geom->sgsiz);
  for (c = 1; c <= geom->sgnum; c++)
    wm[c] = geom->fm[c].wn;


  /*-----------------------------------------------------------------*/
  /* Count wet cells in window n.                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    sc = geom->fm[c].sc;
    if (n > 0 && mask[n][sc] == 0) {
      mask[n][sc] = 1;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Count the open boundary cells (OBC) in window n.                */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      for (cc = 1; cc <= window[n]->open[nn]->no3_t; cc++) {
	ac = window[n]->open[nn]->obc_t[cc];
        if (mask[n][ac] == 0)
          mask[n][ac] = 1;
      }
      window[n]->open[nn]->no2_a = window[n]->open[nn]->no3_a = 0;
      for(cc = 1; cc <= window[n]->open[nn]->no3_t; cc++) {
	c = window[n]->open[nn]->obc_t[cc];
	cs = window[n]->m2d[c];
	for(j = 1; j <= window[n]->npe[cs]; j++) {
	  if (window[n]->wgst[window[n]->c2c[j][c]]) {
	    j--;
	    break;
	  }
	}
	ac = window[n]->c2c[j][c];
	cgm = window[n]->wsa[ac];
	if (wm[cgm] > 0 && wm[cgm] != n && mask[n][ac] == 0) {
	  window[n]->open[nn]->no3_a++;
	  if (ac <= window[n]->enonS) {
	    window[n]->open[nn]->no2_a++;
	  }
	  mask[n][ac] = 1;
	}
      }
    }
  }
  /* Allocate the OBC ghost cells                                    */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      window[n]->open[nn]->obc_a = i_alloc_1d(window[n]->open[nn]->no3_a + 1);
      window[n]->open[nn]->no3_a = 1;
      for(cc = 1; cc <= window[n]->open[nn]->no3_t; cc++) {
	c = window[n]->open[nn]->obc_t[cc];
	cs = window[n]->m2d[c];
	for(j = 1; j <= window[n]->npe[cs]; j++) {
	  if (window[n]->wgst[window[n]->c2c[j][c]]) {
	    j--;
	    break;
	  }
	}
	ac = window[n]->c2c[j][c];
	cgm = window[n]->wsa[ac];
	if (wm[cgm] > 0 && wm[cgm] != n) {
	  window[n]->open[nn]->obc_a[window[n]->open[nn]->no3_a] = ac;
	  window[n]->open[nn]->no3_a++;
	}
      }
    }
  }
  i_free_1d(wm);
  i_free_2d(mask);
}

/* END get_local_obc_a()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to decompose the global open boundarys into OBCs within   */
/* each window.                                                      */
/*-------------------------------------------------------------------*/
void OBC_build(open_bdrys_t **open, /* Global OBC structure          */
	       geometry_t *geom,    /* Global geometry               */
               geometry_t **window, /* Window data structure         */
               int nwindows         /* Number of windows             */
  )
{
  int n, nn, tn, m, j, jj;     /* Counters                           */
  int c, cl, cc, c1, c2;       /* Sparse coordinates / counters      */
  int ee, e;                   /* Edge coordinate, counters          */
  int i1, i2, cy;              /* Interior local sparse coordinates  */
  int g1, g2;                  /* Interior global sparse coordinates */
  int *nobc;                   /* Number of OBC's for each window    */
  int **tti;                   /* Tranfer index counter              */
  short **ff;                  /* Flag for file forcing              */
  int zfe1, zfe2;

  if (!geom->nobc)
    return;
  geom->owc = s_alloc_2d(nwindows + 1, geom->nobc);
  nobc = i_alloc_1d(nwindows + 1);
  ff = s_alloc_2d(geom->nobc, nwindows + 1);
  tti = i_alloc_2d(geom->nobc, nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Count the number of windows containing open boundaries          */
  for (n = 0; n < geom->nobc; n++) {
    for (nn = 1; nn <= nwindows; nn++) {
      geom->owc[n][nn] = -1;
      tti[nn][n] = 1;
    }
  }
  for (n = 0; n < geom->nobc; n++) {
    for (ee = 1; ee <= open[n]->no3_e1; ee++) {
      c = open[n]->obc_e2[ee];  /* Cell corresponding to obc edge    */
      geom->owc[n][geom->fm[c].wn] = n;
    }
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      geom->owc[n][geom->fm[c].wn] = n;

    }
  }
  for (nn = 1; nn <= nwindows; nn++)
    nobc[nn] = 0;
  for (n = 0; n < geom->nobc; n++) {
    for (nn = 1; nn <= nwindows; nn++) {
      if (geom->owc[n][nn] != -1)
        nobc[nn]++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the boundary structures in each window      */
  for (nn = 1; nn <= nwindows; nn++) {
    window[nn]->open = (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * nobc[nn]);
    memset(window[nn]->open, 0, sizeof(open_bdrys_t *) * nobc[nn]);
  }
  for (nn = 1; nn <= nwindows; nn++) {
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n] = OBC_alloc();
      window[nn]->open[n]->bcond_tra = i_alloc_1d(open[0]->ntr);
      ff[nn][n] = 0;
    }
    window[nn]->nobc = nobc[nn];
  }

  /*-----------------------------------------------------------------*/
  /* Fill the window OBC structues with OBC information              */
  for (nn = 1; nn <= nwindows; nn++) {
    c1 = 0;
    for (n = 0; n < geom->nobc; n++) {
      if (geom->owc[n][nn] != -1) {
        strcpy(window[nn]->open[c1]->name, open[n]->name);
        window[nn]->open[c1]->id = open[n]->id;
        window[nn]->open[c1]->type = open[n]->type;
        window[nn]->open[c1]->ocodex = open[n]->ocodex;
        window[nn]->open[c1]->ocodey = open[n]->ocodey;
        window[nn]->open[c1]->ocodec = open[n]->ocodec;
        window[nn]->open[c1]->bcond_nor = open[n]->bcond_nor;
        window[nn]->open[c1]->bcond_nor2d = open[n]->bcond_nor2d;
        window[nn]->open[c1]->bcond_tan = open[n]->bcond_tan;
        window[nn]->open[c1]->bcond_tan2d = open[n]->bcond_tan2d;
        window[nn]->open[c1]->bcond_ele = open[n]->bcond_ele;
        window[nn]->open[c1]->bcond_w = open[n]->bcond_w;
        window[nn]->open[c1]->bcond_Vz = open[n]->bcond_Vz;
        window[nn]->open[c1]->bcond_Kz = open[n]->bcond_Kz;
        window[nn]->open[c1]->clampv = d_alloc_1d(open[n]->ntr);
        window[nn]->open[c1]->relax_time = open[n]->relax_time;
        window[nn]->open[c1]->relax_timei = open[n]->relax_timei;
        window[nn]->open[c1]->relax_zone_nor = open[n]->relax_zone_nor;
        window[nn]->open[c1]->relax_zone_tan = open[n]->relax_zone_tan;
        window[nn]->open[c1]->relax_zone_ele = open[n]->relax_zone_ele;
        window[nn]->open[c1]->linear_zone_nor = open[n]->linear_zone_nor;
        window[nn]->open[c1]->linear_zone_tan = open[n]->linear_zone_tan;
        window[nn]->open[c1]->relax_ele = open[n]->relax_ele;
        window[nn]->open[c1]->rele_b = open[n]->rele_b;
        window[nn]->open[c1]->rele_i = open[n]->rele_i;
        window[nn]->open[c1]->adjust_flux = open[n]->adjust_flux;
        window[nn]->open[c1]->adjust_flux_s = open[n]->adjust_flux_s;
        window[nn]->open[c1]->stagger = open[n]->stagger;
        window[nn]->open[c1]->spf = open[n]->spf;
        window[nn]->open[c1]->meanc = open[n]->meanc;
        window[nn]->open[c1]->bflux_2d = open[n]->bflux_2d;
        window[nn]->open[c1]->bflux_3d = open[n]->bflux_3d;
        window[nn]->open[c1]->inverse_barometer = open[n]->inverse_barometer;
        window[nn]->open[c1]->upmeth = open[n]->upmeth;
        window[nn]->open[c1]->sponge_zone = open[n]->sponge_zone;
        window[nn]->open[c1]->sponge_zone_h = open[n]->sponge_zone_h;
        window[nn]->open[c1]->sponge_f = open[n]->sponge_f;
        window[nn]->open[c1]->ntr = open[n]->ntr;
	window[nn]->open[c1]->bgz = open[n]->bgz;
	window[nn]->open[c1]->rlen = open[n]->rlen;
	window[nn]->open[c1]->options = open[n]->options;
        window[nn]->open[c1]->relax_zone_tra = i_alloc_1d(open[n]->ntr);
        window[nn]->open[c1]->trpc = d_alloc_1d(open[n]->ntr);
        if ((window[nn]->open[c1]->ntt = open[n]->ntt))
          window[nn]->open[c1]->trm = i_alloc_1d(open[n]->ntr);

	window[nn]->open[c1]->ntflx = 0;
        for (tn = 0; tn < open[n]->ntr; tn++) {
          window[nn]->open[c1]->clampv[tn] = open[n]->clampv[tn];
          window[nn]->open[c1]->bcond_tra[tn] = open[n]->bcond_tra[tn];
          window[nn]->open[c1]->relax_zone_tra[tn] =
            open[n]->relax_zone_tra[tn];
          window[nn]->open[c1]->trpc[tn] = open[n]->trpc[tn];
          if (window[nn]->open[c1]->ntt)
            window[nn]->open[c1]->trm[tn] = open[n]->trm[tn];
	  if (open[n]->bcond_tra[tn] & (TRFLUX|TRCONC|TRCONF))
	    window[nn]->open[c1]->ntflx++;
	}

	/* Tracers using a flux OBC                                  */
	if (window[nn]->open[c1]->ntflx)
	  window[nn]->open[c1]->tflx = i_alloc_1d(window[nn]->open[c1]->ntflx);
	window[nn]->open[c1]->ntflx = 0;
        for (tn = 0; tn < open[n]->ntr; tn++) {
	  if (open[n]->bcond_tra[tn] & (TRFLUX|TRCONC|TRCONF)) {
	    window[nn]->open[c1]->tflx[window[nn]->open[c1]->ntflx] = tn;
	    window[nn]->open[c1]->ntflx++;
	  }
	}

        /* Tidal forcing for elevation                               */
	strcpy(window[nn]->open[c1]->tide_con, open[n]->tide_con);
        window[nn]->open[c1]->ntide = open[n]->ntide;
        if (window[nn]->open[c1]->ntide > 0)
          window[nn]->open[c1]->tideforce =
            (tide_details_t *)malloc(sizeof(tide_details_t) *
                                     window[nn]->open[c1]->ntide);
        for (tn = 0; tn < window[nn]->open[c1]->ntide; tn++) {
          window[nn]->open[c1]->tideforce[tn].ntides =
            open[n]->tideforce[tn].ntides;
          strcpy(window[nn]->open[c1]->tideforce[tn].tname,
                 open[n]->tideforce[tn].tname);
          window[nn]->open[c1]->tideforce[tn].ic =
            open[n]->tideforce[tn].ic;
          window[nn]->open[c1]->tideforce[tn].jc =
            open[n]->tideforce[tn].jc;
          window[nn]->open[c1]->tideforce[tn].amp =
            open[n]->tideforce[tn].amp;
          window[nn]->open[c1]->tideforce[tn].per =
            open[n]->tideforce[tn].per;
          window[nn]->open[c1]->tideforce[tn].mta =
            open[n]->tideforce[tn].mta;
          window[nn]->open[c1]->tideforce[tn].dta =
            open[n]->tideforce[tn].dta;
          window[nn]->open[c1]->tideforce[tn].mtp =
            open[n]->tideforce[tn].mtp;
          window[nn]->open[c1]->tideforce[tn].dtp =
            open[n]->tideforce[tn].dtp;
        }

        /* Point the interior cell maps on the boundaries to the     */
        /* correct spatial map.                                      */
	window[nn]->open[c1]->nmap = window[nn]->c2c[geom->open[n]->ocodex];
	window[nn]->open[c1]->omap = window[nn]->c2c[geom->open[n]->ocodey];
	if(geom->npem == 4) {
	  if (geom->open[n]->ocodex == 1 || geom->open[n]->ocodex == 3) {
	    window[nn]->open[c1]->tmpp = window[nn]->c2c[2];
	    window[nn]->open[c1]->tmpm = window[nn]->c2c[4];
	  }
	  if (geom->open[n]->ocodex == 2 || geom->open[n]->ocodex == 4) {
	    window[nn]->open[c1]->tmpp = window[nn]->c2c[3];
	    window[nn]->open[c1]->tmpm = window[nn]->c2c[1];
	  }
	}
        geom->owc[n][nn] = c1;
        c1++;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Fill the open boundary vectors with sparse locations of the     */
  /* open boundaries in that window.                                 */
  /* First count the OBC cells in each window.                       */
  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n]->no3_t = 0;
      window[nn]->open[n]->no3_e1 = 0;
      window[nn]->open[n]->to3_e1 = 0;
      window[nn]->open[n]->no2_t = 0;
      window[nn]->open[n]->no2_e1 = 0;
      window[nn]->open[n]->to2_e1 = 0;
    }

  for (n = 0; n < geom->nobc; n++) {
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      nn = geom->fm[c].wn;
      window[nn]->open[geom->owc[n][nn]]->no3_t++;
      if (geom->zp1[c] == c) {
        window[nn]->open[geom->owc[n][nn]]->no2_t++;
      }
    }

    for (ee = 1; ee <= open[n]->no3_e1; ee++) {
      int wn, ed, mode;
      e = open[n]->obc_e1[ee];
      nn = eiw(geom, e, &mode);
      c = eic(geom, e, nn, &ed);
      window[nn]->open[geom->owc[n][nn]]->no3_e1++;
      if (geom->zp1[c] == c)
        window[nn]->open[geom->owc[n][nn]]->no2_e1++;
    }
    for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to3_e1; ee++) {
      int wn, ed, mode;
      e = open[n]->obc_e1[ee];
      nn = eiw(geom, e, &mode);
      c = eic(geom, e, nn, &ed);
      window[nn]->open[geom->owc[n][nn]]->to3_e1++;
      if (geom->zp1[c] == c)
        window[nn]->open[geom->owc[n][nn]]->to2_e1++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for master - slave transfer vectors             */
  if (nwindows == 1) {
    for (n = 0; n < geom->nobc; n++) {
      window[1]->open[n]->transfer_u1 = open[n]->transfer_u1;
      window[1]->open[n]->transfer_u2 = open[n]->transfer_u2;
      window[1]->open[n]->transfer_u1av = open[n]->transfer_u1av;
      window[1]->open[n]->transfer_u2av = open[n]->transfer_u2av;
      window[1]->open[n]->transfer_eta = open[n]->transfer_eta;
      window[1]->open[n]->t_transfer = open[n]->t_transfer;
      window[1]->open[n]->t_imap = open[n]->t_imap;
      window[1]->open[n]->ttsz = open[n]->ttsz;
    }
  } else {
    for (nn = 1; nn <= nwindows; nn++)
      for (n = 0; n < nobc[nn]; n++) {
	int etasize = window[nn]->open[n]->no2_t + 1;
	int u1size_n = window[nn]->open[n]->no3_e1 + 1;
	int u1size_t = u1size_n + window[nn]->open[n]->to3_e1 + 1;
	int zn = (window[nn]->open[n]->relax_zone_nor) ? window[nn]->open[n]->relax_zone_nor : 1;
	int zt = (window[nn]->open[n]->relax_zone_tan) ? window[nn]->open[n]->relax_zone_tan : 1;

	/* Elevation transfers                                       */
	/*
        if (window[nn]->open[n]->bcond_ele & (FILEIN | CUSTOM) ||
	    window[nn]->open[n]->bcond_nor2d & FLATHR) {
	*/
        if (window[nn]->open[n]->bcond_ele & (FILEIN | CUSTOM)) {
          window[nn]->open[n]->transfer_eta = d_alloc_1d(etasize);
          window[nn]->open[n]->tmap = i_alloc_1d(etasize);
          ff[nn][n] |= 1;
        }

	/* u1 boundary transfers (normal and tangential)             */
        if (window[nn]->open[n]->type & U1BDRY) {
	  g1 = 0;
	  if (window[nn]->open[n]->bcond_nor & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u1 = d_alloc_1d(zn * u1size_n);
	    g1 |= 1;
	  }
	  if (window[nn]->open[n]->bcond_nor2d & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u1av = d_alloc_1d(u1size_n);
	    g1 |= 2;
	  }
	  if (window[nn]->open[n]->bcond_tan & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u2 = d_alloc_1d(zt * u1size_t);
	    g1 |= 4;
	  }
	  if (window[nn]->open[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u2av = d_alloc_1d(u1size_t);
	    g1 |= 4;
	  }
	  if(g1 & (1|2)) {
	    window[nn]->open[n]->tmap_u1 = i_alloc_1d(u1size_n);
	    ff[nn][n] |= 2;
	  }
	  if(g1 & 4) {
	    window[nn]->open[n]->tmap_u2 = i_alloc_1d(u1size_t);
	    ff[nn][n] |= 16;
	  }
        }

	/* Tracer transfers                                          */
        if (window[nn]->open[n]->ntt) {
	  m = (window[nn]->open[n]->bgz) ? window[nn]->open[n]->bgz : 1;
	  window[nn]->open[n]->ttsz = window[nn]->open[n]->no3_t + 
	    m * window[nn]->open[n]->no3_e1 + 1;
          window[nn]->open[n]->t_transfer = 
	    d_alloc_2d(window[nn]->open[n]->ttsz, 
		       window[nn]->open[n]->ntt);
          window[nn]->open[n]->t_tmap =
            i_alloc_1d(window[nn]->open[n]->ttsz);
	  window[nn]->open[n]->t_imap = 
	    i_alloc_2d(m + 1, window[nn]->open[n]->ttsz);
          ff[nn][n] |= 8;
        }
      }
  }

  /* Allocate memory for the open boundary vectors                   */
  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n]->obc_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->oi1_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->oi2_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->cyc_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->nepc = 
        d_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->inc = 
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->bec = 
        i_alloc_2d(window[nn]->open[n]->no3_t + 1, window[nn]->npem + 1);
      window[nn]->open[n]->bcc = 
        i_alloc_2d(window[nn]->open[n]->no3_t + 1, window[nn]->npem + 1);
      window[nn]->open[n]->olap = i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->no3_t = 1;

      c1 = window[nn]->open[n]->no3_e1 + window[nn]->open[n]->to3_e1 + 1;
      window[nn]->open[n]->obc_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->obc_e2 = i_alloc_1d(c1);
      window[nn]->open[n]->oi1_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->oi2_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->cyc_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->e2c_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->ogc_t = i_alloc_1d(window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->ini = i_alloc_1d(window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->outi = i_alloc_1d(window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->ceni = i_alloc_1d(window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->dir = i_alloc_1d(window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->nmape = (int **)p_alloc_1d(window[nn]->open[n]->no3_e1+1);
      window[nn]->open[n]->omape = (int **)p_alloc_1d(window[nn]->open[n]->no3_e1+1);
      window[nn]->open[n]->to3_e1 = window[nn]->open[n]->no3_e1 + 1;
      window[nn]->open[n]->no3_e1 = 1;
    }

  /*-----------------------------------------------------------------*/
  /* Fill the open boundary vectors with local sparse locations      */
  for (n = 0; n < geom->nobc; n++) {

    /* Tracer open boundary vectors                                  */
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];
      window[nn]->open[c1]->mwn = n;
      if (ff[nn][c1] & 1 && cc <= open[n]->no2_t)
        window[nn]->open[c1]->tmap[window[nn]->open[c1]->no3_t] = cc;
      if (ff[nn][c1] & 8 && window[nn]->open[c1]->ntt) {
        window[nn]->open[c1]->t_tmap[tti[nn][c1]] = cc;
        window[nn]->open[c1]->t_imap[window[nn]->open[c1]->no3_t][0] = tti[nn][c1];
	tti[nn][c1]++;
      }
      cl = geom->fm[c].sc;
      cy = open[n]->cyc_t[cc];
      cy = geom->fm[cy].sc;
      /* Cell centre coordinate                                      */
      window[nn]->open[c1]->obc_t[window[nn]->open[c1]->no3_t] = cl;
      /* Number of boundary edges around this centre                 */
      window[nn]->open[c1]->nepc[window[nn]->open[c1]->no3_t] = open[n]->nepc[cc];
      /* Average index mapping into the interior                     */
      window[nn]->open[c1]->inc[window[nn]->open[c1]->no3_t] = open[n]->inc[cc];
      /* Map containing boundary edges around this centre            */
      for (j = 1; j <= geom->npe[geom->m2d[c]]; j++) {
	window[nn]->open[c1]->bec[j][window[nn]->open[c1]->no3_t] = 0;
	window[nn]->open[c1]->bcc[j][window[nn]->open[c1]->no3_t] = 0;
	if (geom->wgst[geom->c2c[j][c]]) {
	  if (open[n]->bec[j][cc]) {
	    window[nn]->open[c1]->bec[j][window[nn]->open[c1]->no3_t] = window[nn]->c2e[j][cl];
	    window[nn]->open[c1]->bcc[j][window[nn]->open[c1]->no3_t] = window[nn]->c2c[j][cl];
	  }
	} else {
	  if (open[n]->bcc[j][cc] < 0)
	    window[nn]->open[c1]->bcc[j][window[nn]->open[c1]->no3_t] = -window[nn]->c2c[j][cl];
	}
      }
      /* Boundary neighbours                                         */
      window[nn]->open[c1]->olap[window[nn]->open[c1]->no3_t] = open[n]->olap[cc];
      /* Cyclic centre                                               */
      window[nn]->open[c1]->cyc_t[window[nn]->open[c1]->no3_t] = cy;
      /* Interior centres                                            */
      g1 = open[n]->oi1_t[cc];
      i1 = get_icell(geom, window[nn], cl, open[n]->obc_t[cc], g1);
      window[nn]->open[c1]->oi1_t[window[nn]->open[c1]->no3_t] = i1;

      g2 = open[n]->oi2_t[cc];
      i2 = get_icell(geom, window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_t[window[nn]->open[c1]->no3_t] = i2;

      window[nn]->open[c1]->no3_t++;
    }

    /* Normal u1 velocity open boundary vectors                      */
    for (ee = 1; ee <= open[n]->no3_e1; ee++) {
      int el, ei, ed, j, mode;
      e = open[n]->obc_e1[ee];
      nn = eiw(geom, e, &mode);
      c = eic(geom, e, nn, &ed);
      c1 = geom->owc[n][nn];

      if (ff[nn][c1] & 2)
        window[nn]->open[c1]->tmap_u1[window[nn]->open[c1]->no3_e1] = ee;

      cl = geom->fm[c].sc;
      el = window[nn]->c2e[ed][cl];
      cy = geom->g2we[nn][open[n]->cyc_e1[ee]];
      /* Normal edge coordinate                                      */
      window[nn]->open[c1]->obc_e1[window[nn]->open[c1]->no3_e1] = el;
      /* Centre corresponding to the normal edge                     */
      window[nn]->open[c1]->obc_e2[window[nn]->open[c1]->no3_e1] = cl;
      /* Cyclic edge                                                 */
      window[nn]->open[c1]->cyc_e1[window[nn]->open[c1]->no3_e1] = cy;
      /* Index pointing to the interior of the domain                */
      window[nn]->open[c1]->ini[window[nn]->open[c1]->no3_e1] = open[n]->ini[ee];
      /* Index pointing out of the domain                            */
      window[nn]->open[c1]->outi[window[nn]->open[c1]->no3_e1] = open[n]->outi[ee];
      /* Edge orientation                                            */
      window[nn]->open[c1]->ceni[window[nn]->open[c1]->no3_e1] = open[n]->ceni[ee];
      window[nn]->open[c1]->dir[window[nn]->open[c1]->no3_e1] = open[n]->dir[ee];
      /* Edge map into the interior                                  */
      j = open[n]->ini[ee];
      window[nn]->open[c1]->nmape[window[nn]->open[c1]->no3_e1] = window[nn]->c2c[j];
      /* Edge map out of the domain                                  */
      j = open[n]->outi[ee];
      window[nn]->open[c1]->omape[window[nn]->open[c1]->no3_e1] = window[nn]->c2c[j];
      /* First boundary ghost cell centre adjacent to this edge      */
      window[nn]->open[c1]->ogc_t[window[nn]->open[c1]->no3_e1] = window[nn]->c2c[j][cl];
      /* Get the centre index (cc) corresponding to the edge index.  */
      /* Note: open->no3_t still needs to be decremented, so use <   */
      cy = window[nn]->e2c[window[nn]->c2e[ed][cl]][open[n]->ini[ee]];
      for (cc = 1; cc < window[nn]->open[c1]->no3_t; cc++) {
	if (cl == window[nn]->open[c1]->obc_t[cc]) {
	  window[nn]->open[c1]->e2c_e1[window[nn]->open[c1]->no3_e1] = cc;
	  break;
	}
      }
      /* Interior edges                                              */
      /*ei = jo(ed, geom->npe[geom->m2d[c]]);*/
      ei = jow(geom, c, ed);
      cl = window[nn]->c2c[ei][cl];
      el = window[nn]->c2e[ei][cl];
      window[nn]->open[c1]->oi1_e1[window[nn]->open[c1]->no3_e1] = el;
      cl = window[nn]->c2c[ei][cl];
      el = window[nn]->c2e[ei][cl];
      window[nn]->open[c1]->oi2_e1[window[nn]->open[c1]->no3_e1] = el;
      /* Tracer transfer map                                         */
      if (ff[nn][c1] & 8 && window[nn]->open[c1]->ntt) {
	for (m = 1; m <= window[nn]->open[c1]->bgz; m++) {
	  window[nn]->open[c1]->t_imap[window[nn]->open[c1]->no3_e1][m] = tti[nn][c1];
	  window[nn]->open[c1]->t_tmap[tti[nn][c1]] = open[n]->t_imap[ee][m];
	  tti[nn][c1]++;
	}
      }

      window[nn]->open[c1]->no3_e1++;
    }

    /* Tangential u1 velocity open boundary vectors                  */
    for (ee = open[n]->no3_e1 + 1; ee <= open[n]->to3_e1; ee++) {
      int ed, mode;
      e = open[n]->obc_e1[ee];
      nn = eiw(geom, e, &mode);
      c = eic(geom, e, nn, &ed);
      c1 = geom->owc[n][nn];

      if (ff[nn][c1] & 16) {
        window[nn]->open[c1]->tmap_u2[window[nn]->open[c1]->to3_e1] = ee;
      }

      cl = geom->fm[c].sc;
      cy = geom->g2we[nn][open[n]->cyc_e1[ee]];

      /* Tangential edge coordinate                                  */
      window[nn]->open[c1]->obc_e1[window[nn]->open[c1]->to3_e1] = window[nn]->c2e[ed][cl];
      /* Centre corresponding to the tangential edge                 */
      window[nn]->open[c1]->obc_e2[window[nn]->open[c1]->to3_e1] = cl;
      /* Cyclic edge                                                 */
      window[nn]->open[c1]->cyc_e1[window[nn]->open[c1]->to3_e1] = cy;
      /* Get the centre index (cc) corresponding to the edge index   */
      for (cc = 1; cc < window[nn]->open[c1]->no3_t; cc++) {
	if (cl == window[nn]->open[c1]->obc_t[cc]) {
	  window[nn]->open[c1]->e2c_e1[window[nn]->open[c1]->no3_e1] = cc;
	  jj = window[nn]->open[c1]->inc[cc];
	  break;
	}
      }
      /* Tangential interiors are found by finding the centre for    */
      /* the edge (cl) and the index the edge has from this centre   */
      /* (j), then propagating the centre into the interior in its   */
      /* the average direction (jj), then recovering the edge j of   */
      /* that cell. For non-recangular grids this doesn't propagate  */
      /* perpendicularly into the interior.                          */
      j = geom->e2e[e][0];
      i1 = window[nn]->c2c[jj][cl];

      window[nn]->open[c1]->oi1_e1[window[nn]->open[c1]->to3_e1] = window[nn]->c2e[j][i1];

      i2 = window[nn]->c2c[jj][i1];
      window[nn]->open[c1]->oi2_e1[window[nn]->open[c1]->to3_e1] = window[nn]->c2e[j][i2];
      window[nn]->open[c1]->to3_e1++;
    }
  }

  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n]->no3_t -= 1;
      window[nn]->open[n]->no3_e1 -= 1;
      window[nn]->open[n]->to2_e1 += window[nn]->open[n]->no3_e1;
      window[nn]->open[n]->to3_e1 -= 1;

      /* Set self-mapping maps over all edges for no OBC ghost zones */
      /*
      if (window[nn]->open[n]->bgz == 0) {
	for (ee = 1; ee <= window[nn]->open[n]->no3_e1; ee++) {
	  int ed, m;
	  e = window[nn]->open[n]->obc_e1[ee];
	  c = eic(geom, e, nn, &ed);
	  if (window[nn]->open[n]->stagger & INFACE) c = window[nn]->c2c[ed][c];
	  window[nn]->c2c[ed][c] = c;
	}
      }
      */
      /* Get the bottom coordinate vector                            */
      window[nn]->open[n]->cdeep = 1;
      window[nn]->open[n]->bot_t = i_alloc_1d(window[nn]->open[n]->no2_t + 1);
      for (cc = 1; cc <= window[nn]->open[n]->no2_t; cc++) {
	c = c2 = window[nn]->open[n]->obc_t[cc];
	while(c != window[nn]->zm1[c])
	  c = window[nn]->zm1[c];
	window[nn]->open[n]->bot_t[cc] = window[nn]->zp1[c];
	cl = window[nn]->open[n]->obc_t[window[nn]->open[n]->cdeep];
	if(window[nn]->botz[c2] < window[nn]->botz[cl])
	  window[nn]->open[n]->cdeep = cc;
      }

      cl = max(window[nn]->open[n]->no3_t + 1, window[nn]->open[n]->to3_e1 + 1);

      window[nn]->open[n]->i1 = i_alloc_1d(cl);
      window[nn]->open[n]->dum = d_alloc_1d(cl);
      window[nn]->open[n]->dum1 = d_alloc_1d(cl);
      window[nn]->open[n]->dumn = d_alloc_2d(cl, window[nn]->open[n]->ntr);

      window[nn]->open[n]->u1d =
	d_alloc_1d(window[nn]->open[n]->to3_e1 + 1);
      window[nn]->open[n]->flow =
	d_alloc_1d(window[nn]->open[n]->no3_t + 1);

      /* Dummy for tracers using a flux OBC                          */
      if (window[nn]->open[n]->ntflx) {
	cl = max(window[nn]->open[n]->no3_t + 1, window[nn]->open[n]->to3_e1 + 1);
	window[nn]->open[n]->dumtr = d_alloc_2d(cl, window[nn]->open[n]->ntflx);
      }
	  
      /* Allocate the phase smoothing array if required */
      if (window[nn]->open[n]->spf) {
	window[nn]->open[n]->sphase =
	  d_alloc_1d(window[nn]->open[n]->no2_t + 1);
	memset(window[nn]->open[n]->sphase, 0, 
	       (window[nn]->open[n]->no2_t + 1) * sizeof(double));
      }
    }

  /* Set the global map geom->fm for OBC ghost cells. This ensures   */
  /* these cells are treated as wet cells when the transfer maps are */
  /* created in build_transfer_maps().                               */
  for (nn = 1; nn <= nwindows; nn++) {
    for (n = 0; n < window[nn]->nobc; n++) {
      open_bdrys_t *open = window[nn]->open[n];
      for (ee = 1; ee <= open->no3_e1; ee++) {
	cl = open->ogc_t[ee];
	c = window[nn]->wsa[cl];
	geom->fm[c].wn = nn;
	geom->fm[c].sc = cl;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Print OBC information                                           */
  if (DEBUG("init_w")) {
    dlog("init_w", "\nWindow open boundaries\n");
    for (nn = 1; nn <= nwindows; nn++) {
      dlog("init_w", "  Window %d = %d open boundaries\n", nn,
           window[nn]->nobc);
      for (n = 0; n < window[nn]->nobc; n++) {
        dlog("init_w", "    boundary #%d : %s\n", n,
             window[nn]->open[n]->name);
        dlog("init_w", "      Tracer = %d cells\n",
             window[nn]->open[n]->no3_t);
        dlog("init_w",
             "      u1 velocity = %d normal cells, %d tangential cells\n",
             window[nn]->open[n]->no3_e1,
             window[nn]->open[n]->to3_e1 - window[nn]->open[n]->no3_e1);
        dlog("init_w", "      Elevation = %d cells\n",
             window[nn]->open[n]->no2_t);
        dlog("init_w",
             "      u1av velocity = normal %d cells, %d tangential cells\n",
             window[nn]->open[n]->no2_e1,
             window[nn]->open[n]->to2_e1 - window[nn]->open[n]->no3_e1);
      }
    }
  }

  i_free_1d(nobc);
  s_free_2d(ff);
  i_free_2d(tti);
}

/* END OBC_build()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Get the cells and weights within the sponge zones                 */
/*-------------------------------------------------------------------*/
void set_sponge_cells(geometry_t *window)
{
  int nn, n, cc, c;
  int ee, e, eb, ep;
  int i1, i2, cn, cb, j;
  double d1, d2, dist, dmin;
  int nscm, *scm;
  int *mask;
  double deg2m = 60.0 * 1852.0;
  int has_proj = (strlen(projection) > 0);
  int is_geog = has_proj && (strcasecmp(projection, GEOGRAPHIC_TAG) == 0);
  int verbose = 0;

  mask = i_alloc_1d(window->szm);
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->sponge_zone_h) {

      /*-------------------------------------------------------------*/
      /* Count the cell centres in the zone                          */
      open->nspc = 0;
      memset(mask, 0, window->szm * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	for (i1 = 1; i1 <= open->no2_t; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (!mask[c] && dist <= (double)open->sponge_zone_h) {
	    open->nspc++;
	    mask[c] = 1;
	  }
	}
      }
      
      /*-------------------------------------------------------------*/
      /* Get the cells in the zone                                   */
      open->spc = i_alloc_1d(open->nspc + 1);
      open->swc = d_alloc_1d(open->nspc + 1);
      open->snc = i_alloc_1d(open->nspc + 1);
      open->smc = i_alloc_1d(open->nspc + 1);
      open->nspc = 1;
      memset(mask, 0, window->szm * sizeof(int));
      for (cc = 1; cc <= window->b2_t; cc++) {
	c = window->w2_t[cc];
	for (i1 = 1; i1 <= open->no2_e1; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (!mask[c] && dist <= (double)open->sponge_zone_h) {
	    open->spc[open->nspc++] = c;
	    mask[c] = 1;
	  }
	}
      }
      open->nspc--;

      /*-------------------------------------------------------------*/
      /* Set up the mask                                             */
      memset(mask, 0, window->szm * sizeof(int));
      for (cc = window->v2_t+1; cc <= window->n2_t; cc++) {
	c = window->w2_t[cc];
	mask[c] = 1;
      }
      for (cc = window->b2_t+1; cc <= window->a2_t; cc++) {
	c = window->w2_t[cc];
	mask[c] = 0;
      }
      for (ee = 1; ee <= open->no2_e1; ee++) {
	c = open->obc_e2[ee];
	cn = open->ogc_t[ee];
	mask[c] = mask[cn] = 1;
      }

      /*-------------------------------------------------------------*/
      /* Get the closest boundary cell to sponge cells               */
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	dmin = HUGE;
	for (i1 = 1; i1 <= open->no2_t; i1++) {
	  i2 = open->obc_t[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (dist <= dmin) {
	    dmin = dist;
	    open->swc[cc] = dist;
	    open->snc[cc] = i2;
	  }
	}
	mask[c] = 1;
      }

      /*-------------------------------------------------------------*/
      /* Find the cells on the outer limit of the zone               */
      nscm = 0;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  cn = window->c2c[j][c];
	  if (!mask[cn]) {
	    nscm++;
	    break;
	  }
	}
      }
      scm = i_alloc_1d(nscm + 1);
      nscm = 1;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  cn = window->c2c[j][c];
	  if (!mask[cn]) {
	    scm[nscm++] = c;
	    break;
	  }
	}
      }
      nscm--;

      /*-------------------------------------------------------------*/
      /* Get the closest cell on the outer perimeter, its distance   */
      /* to the closest cell on the boundary, and the ratio of       */
      /* distances (weight). The value in the sponge zone is then:   */
      /* v = swc * (v(perimeter) - v(boundary))  + v(boundary)       */
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	dmin = HUGE;
	for (i1 = 1; i1 <= nscm; i1++) {
	  i2 = scm[i1];
	  d1 = window->cellx[c] - window->cellx[i2];
	  d2 = window->celly[c] - window->celly[i2];
	  dist = sqrt(d1 * d1 + d2 * d2);
	  if (is_geog) dist *= deg2m;
	  if (dist <= dmin) {
	    dmin = dist;
	    open->smc[cc] = i2;

	    cn = open->snc[cc];
	    d1 = window->cellx[cn] - window->cellx[c];
	    d2 = window->celly[cn] - window->celly[c];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swc[cc] = dist;

	    d1 = window->cellx[cn] - window->cellx[i2];
	    d2 = window->celly[cn] - window->celly[i2];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swc[cc] /= dist;
	  }
	}
      }
      if (verbose) {
	for (cc = 1; cc <= open->nspc; cc++) {
	  c = open->spc[cc];
	  printf("wn=%d c=(%d %d) mn=(%d %d) mx=(%d %d) %f\n",window->wn,window->s2i[c],window->s2j[c],
		 window->s2i[open->snc[cc]],window->s2j[open->snc[cc]],
		 window->s2i[open->smc[cc]],window->s2j[open->smc[cc]],open->swc[cc]);
	}
      }

      /*-------------------------------------------------------------*/
      /* Count the edges in the zone                                 */
      memset(mask, 0, window->szm * sizeof(int));
      open->nspe1 = 0;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  e = window->c2e[j][c];
	  if (!mask[e]) {
	    open->nspe1++;
	    mask[e] = 1;
	  }
	}
      }
      open->spe1 = i_alloc_1d(open->nspe1 + 1);
      open->swe1 = d_alloc_1d(open->nspe1 + 1);
      open->sne1 = i_alloc_1d(open->nspe1 + 1);
      open->sme1 = i_alloc_1d(open->nspe1 + 1);
      memset(mask, 0, window->szm * sizeof(int));
      open->nspe1 = 1;
      for (cc = 1; cc <= open->nspc; cc++) {
	c = open->spc[cc];
	for (j = 1; j <= window->npe[c]; j++) {
	  e = window->c2e[j][c];
	  if (!mask[e]) {
	    cn = open->snc[cc];

	    /* Get the edge corresponding to the closest cell on the */
	    /* boundary (i.e. the closest edge on the boundary).     */
	    for (ee = 1; ee <= open->no2_e1; ee++) {
	      eb = open->obc_e1[ee];
	      cb = open->obc_e2[ee];
	      if (cb == cn) {
		open->sne1[open->nspe1] = eb;
		break;
	      }
	    }

	    /* Distance to the closest boundary edge                 */
	    d1 = window->u1x[e] - window->u1x[eb];
	    d2 = window->u1y[e] - window->u1y[eb];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swe1[open->nspe1] = dist;

	    /* Get the edge corresponding to the closest cell on the */
	    /* outer perimeter (i.e. the closest edge on the outer   */
	    /* perimeter).                                           */
	    dmin = HUGE;
	    for (i1 = 1; i1 <= nscm; i1++) {
	      i2 = scm[i1];
	      for (ee = 1; ee <= window->npe[i2]; ee++) {
		ep = window->c2e[ee][i2];
		d1 = window->u1x[e] - window->u1x[ep];
		d2 = window->u1y[e] - window->u1y[ep];
		dist = sqrt(d1 * d1 + d2 * d2);
		if (is_geog) dist *= deg2m;
		if (dist < dmin) {
		  dmin = dist;
		  open->sme1[open->nspe1] = ep;
		}
	      }
	    }
	    ep = open->sme1[open->nspe1];
	    d1 = window->u1x[ep] - window->u1x[eb];
	    d2 = window->u1y[ep] - window->u1y[eb];
	    dist = sqrt(d1 * d1 + d2 * d2);
	    if (is_geog) dist *= deg2m;
	    open->swe1[open->nspe1] /= dist;
	    /*
	      if(window->e2e[e][0]==1&&window->s2j[c]==10)
	      printf("%d %d %d %d %f\n",window->s2i[c],e,ep,open->sne1[open->nspe1],open->swe1[open->nspe1]);
	    */
	    open->spe1[open->nspe1++] = e;
	    mask[e] = 1;
	  }
	}
      }
      i_free_1d(scm);
    }
  }
  i_free_1d(mask);
}

/* END set_sponge_cells()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return an interior cell in a window to any given cell  */
/* on the boundary. The interior cell is located based on the global */
/* maps.                                                             */
/*-------------------------------------------------------------------*/
int get_icell(geometry_t *geom,   /* Global geometry                 */
	      geometry_t *window, /* Window data structure           */
              int cl,           /* Local cell on the boundary        */
              int cg,           /* Global cell on the boundary       */
              int ci            /* Global cell interior to cg        */
  )
{
  int i1 = 0;                   /* Local interior cell               */
  int n;
  int cs = geom->m2d[cg];

  for (n = 1; n <= geom->npe[cs]; n++) {
    if (ci == geom->c2c[n][cg])
      i1 = window->c2c[n][cl];
  }
  return (i1);
}

/* END get_icell()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the local ghost vectors. The approach used here is */
/* to identify if a ghost cell is included in a window's sparse      */
/* array, and find the corresponding interior cells via the local    */
/* maps and the global lateral boundary directions rather than use a */
/* a direct comparison of the global interior cell with the window's */
/* sparse array. This maps of ghost cells to interior cells which    */
/* outside the window self mapping, and sets nested sparse array     */
/* interior cells correctly.                                         */
/*-------------------------------------------------------------------*/
void get_local_ghost(geometry_t *geom,   /* Global; geometry         */
		     geometry_t *window, /* Local geometry           */
		     int n               /* Window number            */
		     )
{
  int c, cc;                    /* Cell centre / counter             */
  int lc, ic1, ic2;             /* Local sparse coordinates          */
  int e, le, ei, ee;            /* Cell edge / counter               */
  int nb, nc, ncS, j;           /* Counters                          */
  int dr;
  int *map;

  if (DEBUG("init_w")) {
    dlog("init_w", "\n");
    dlog("init_w", "  Start making local ghost vectors : window %d\n", n);
  }

  if (geom->nwindows == 1) {
    window->nbpt = geom->nbpt;
    window->nbptS = geom->nbptS;
    window->bpt = geom->bpt;
    window->bin = geom->bin;
    window->bin2 = geom->bin2;
    window->dbpt = geom->dbpt;
    window->wgst = geom->wgst;

    window->nbpte1 = geom->nbpte1;
    window->nbpte1S = geom->nbpte1S;
    window->bpte1 = geom->bpte1;
    window->bine1 = geom->bine1;
    window->bpte1S = geom->bpte1S;
    window->bine1S = geom->bine1S;
    window->nbe1 = geom->nbe1;
    window->nbe1S = geom->nbe1S;

    window->nbpte2 = geom->nbpte2;
    window->nbpte2S = geom->nbpte2S;
    window->bpte2 = geom->bpte2;
    window->bine2 = geom->bine2;
    window->bpte2S = geom->bpte2S;
    window->bine2S = geom->bine2S;
    window->nbe2 = geom->nbe2;
    window->nbe2S = geom->nbe2S;
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Get the lateral ghost vectors for cell centers.                 */
  /* The lateral ghost cells that are included in a window's sparse  */
  /* array comprise those cells that are east, west, north, south    */
  /* and the corners south-west, north-west and south-east from a    */
  /* wet cell in the window. Note that the corner maps are mapped    */
  /* only once while the other maps are mapped laux times. Also note */
  /* that the map north-east doesn't exist since it is not needed.   */
  if (DEBUG("init_w"))
    dlog("init_w", "    Start cell center local ghost vectors\n");
  /* Count the local ghost cells in window n                         */
  window->nbpt = 0;
  window->nbptS = 0;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];          /* Global boundary cell              */
    /* Note : cannot use ANY(c,window->w3_t) here since the cells to */
    /* process vectors do not contain all ghost cells, only those    */
    /* required to specify the fluxes used to update cell centers,   */
    /* i.e. eastern and northern edge ghost cells.                   */
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell  */
    if (lc) {
      window->nbpt++;
      if (lc <= window->enonS)
        window->nbptS++;
    }
  }

  /* Allocate memory and initialise counters                         */
  window->bpt = i_alloc_1d(window->nbpt + 1);
  memset(window->bpt, 0, window->nbpt * sizeof(int));
  window->bin = i_alloc_1d(window->nbpt + 1);
  memset(window->bin, 0, window->nbpt * sizeof(int));
  window->bin2 = i_alloc_1d(window->nbpt + 1);
  memset(window->bin2, 0, window->nbpt * sizeof(int));
  window->dbpt = i_alloc_1d(window->nbpt + 1);
  memset(window->dbpt, 0, window->nbpt * sizeof(int));
  window->wgst = i_alloc_1d(window->enon + 1);
  window->nbpt = window->nbptS + 1;
  window->nbptS = 1;

  /* Assign the lateral ghost cells to the boundary vectors          */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];            /* Global boundary cell            */
    dr = geom->dbin[cc];          /* Direction of bin from bpt       */
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell  */
    if (lc) {
      ic1 = window->c2c[dr][lc];  /* Interior cell to lc             */
      ic2 = window->c2c[dr][ic1]; /* Interior cell to ic1            */
      if (lc <= window->enonS) {
        window->bpt[window->nbptS] = lc;
        window->bin[window->nbptS] = ic1;
        window->bin2[window->nbptS] = ic2;
        window->dbpt[window->nbptS] = geom->dbpt[cc];
        window->nbptS++;
      } else {
        window->bpt[window->nbpt] = lc;
        window->bin[window->nbpt] = ic1;
        window->bin2[window->nbpt] = ic2;
        window->dbpt[window->nbpt] = geom->dbpt[cc];
        window->nbpt++;
      }
      window->wgst[lc] = ic1;
    }
  }
  window->nbpt--;
  window->nbptS--;

  /* Assign 2D OBC ghost cells to the window ghost array             */
  for (nb = 0; nb < geom->nobc; nb++) {
    open_bdrys_t *open = geom->open[nb];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      c = open->obc_e2[ee];
      lc = ANYf(c, window->wsa, window->enon);
      if (lc) {
	j = open->outi[ee];
	window->wgst[window->c2c[j][lc]] = lc;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the lateral ghost vectors for e1 face centered cells.       */
  /* The lateral ghost cells that are included in a window's sparse  */
  /* array for e1 faces actually include wet and auxilliary cells on */
  /* western edges (southern edges for e2 faces). This means that    */
  /* some auxiliary cells may be identified as ghost cells for e1    */
  /* faces at locations where a solid boundary exists at the limits  */
  /* of the window. Consequently these e1 ghost cells are included   */
  /* in the boundary array, but the corresponding cell center ghost  */
  /* cell is not because the directional maps don't include this     */
  /* ghost in sparse array. This results in different numbers of     */
  /* lateral ghost cells included in the cell and face centered      */
  /* boundary arrays. This only applies for window partitions and    */
  /* does not apply for the global vectors or nwindows=1.            */
  /*                                                                 */
  /*                                                                 */
  /*     e.g. plan view     |   ]   |   |   |                        */
  /*                        |   ]   |   |   |                        */
  /*                        |---]---|---|---|                        */
  /*                        |   ] a | a | a |                        */
  /*     Extra e1 ghost  -> |   ] e |   |   |                        */
  /*     in this layer.     |---]---|---|---|  -- = cells            */
  /*     This ghost is   -> | g ] a | a | a |  __ = window limit     */
  /*     derived from the   |___|_e_|___|___|   ] = solid boundary   */
  /*     xm1[yp1] map.      |---]---|---|---|                        */
  /*                        | g ] w | w | w |   g = centered ghost   */
  /*                        |   ] e |   |   |   w = wet cell         */
  /*                        |---]---|---|---|   a = auxiliary cell   */
  /*                        |   ]   |   |   |   e = e1 ghost         */
  /*                                                                 */
  if (DEBUG("init_w"))
    dlog("init_w", "    Start e1 centered local ghost vectors\n");

  /* Count the local e1 ghost cells in window n                      */
  window->nbpte1 = 0;
  window->nbpte1S = 0;
  for (ee = 1; ee <= geom->nbpte1; ee++) {
    e = geom->bpte1[ee];        /* Global boundary edge              */
    le = ANYf(e, window->wse, window->sze-1); /* Local boundary cell */
    if (le) {
      window->nbpte1++;
      if (le < window->szeS)
        window->nbpte1S++;
    }
  }

  /* Allocate memory and initialise counters                         */
  window->bpte1 = i_alloc_1d(window->nbpte1 + 1);
  memset(window->bpte1, 0, window->nbpte1 * sizeof(int));
  window->bine1 = i_alloc_1d(window->nbpte1 + 1);
  memset(window->bine1, 0, window->nbpte1 * sizeof(int));
  window->bpte1S = i_alloc_1d(window->nbpte1S + 1);
  memset(window->bpte1S, 0, window->nbpte1S * sizeof(int));
  window->bine1S = i_alloc_1d(window->nbpte1S + 1);
  memset(window->bine1S, 0, window->nbpte1S * sizeof(int));

  nc = ncS = 1;
  window->nbe1 = window->nbe1S = 0;

  /* Assign the lateral ghost cells to the boundary vectors          */
  for (ee = 1; ee <= geom->nbpte1; ee++) {
    e = geom->bpte1[ee];        /* Global boundary edge              */
    ei = geom->bin[ee];         /* Global interior edge              */
    map = (e == geom->em[ei]) ? window->ep : window->em;
    le = ANYf(e, window->wse, window->sze-1); /* Local boundary cell */
    if (le) {
      /* It is possible that the lateral ghost cell for e1 faces     */
      /* corresponds to an auxiliary cell in the sparse array (see   */
      /* above) and at locations where the bathymetry is stepped,    */
      /* the interior cell to this ghost cell may be outside the     */
      /* limits of the window. These ghost cell's are never used in  */
      /* the window's computations, so leave these cells self-       */
      /* mapping.                                                    */
      ic1 = map[le];            /* Local interior cell to le         */
      window->bpte1[nc] = le;
      window->bine1[nc] = ic1;
      nc++;
      if (ee <= geom->nbe1)
        window->nbe1++;

      if (le < window->szeS) {
        window->bpte1S[ncS] = le;
        window->bine1S[ncS] = ic1;
        ncS++;
        if (ee <= geom->nbe1S)
          window->nbe1S++;
      }
    }
  }
  if (nc - 1 != window->nbpte1)
    hd_warn
      ("get_local_ghosts: Incorrect number of 3D e1 ghost cells located, window %d : %d %d\n",
       n, nc - 1, window->nbpte1);
  if (ncS - 1 != window->nbpte1S)
    hd_warn
      ("get_local_ghosts: Incorrect number of 2D e1 ghost cells located, window %d : %d %d\n",
       n, ncS - 1, window->nbpte1S);

  /*-----------------------------------------------------------------*/
  /* Sediment ghosts under auxiliary and lateral ghost cells         */
  window->ngsed = 0;
  map = i_alloc_1d(window->szc);
  memset(map, 0, window->szc);
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    map[c] = 1;
  }
  for (c = 1; c <= window->enonS; c++) {
    if (!map[c]) window->ngsed++;
  }
  cc = 1;
  window->gsed_t = i_alloc_1d(window->ngsed + 1);
  for (c = 1; c <= window->enonS; c++) {
    if (!map[c]) {
      lc = c;
      while(lc != window->zm1[lc])
	lc = window->zm1[lc];
      window->gsed_t[cc] = lc;
      cc++;
    }
  }
  i_free_1d(map);

  if (DEBUG("init_w"))
    dlog("init_w", "  Local ghost vectors OK : %d cell centers\n\n",
         window->nbpt);
}

/* END get_local_ghost()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to generate the vectors containing the local sparse       */
/* coordinates of the free surface and bottom layer for each window. */
/* Note that the surface map is dynamic and must be reset at each    */
/* timestep as the free surface rises and falls through layers.      */
/*-------------------------------------------------------------------*/
void surfbot_build(geometry_t *geom,    /* Global geometry           */
		   geometry_t *window,  /* Window data structure     */
                   int wn,      /* Window number                     */
                   int *sur,    /* Surface cell vector               */
                   int *bot,    /* Bottom cell vector                */
                   int *vec,    /* Global bottom cells vector        */
                   int nvec,    /* Size of vec                       */
                   int *vec2D,  /* 2D cells to process vector        */
                   int nvec2D,  /* Number of wet cells in vec2D      */
                   int avec2D   /* Wet + auxiliary cells in vec2D    */
  )
{
  int c;                      /* Local sparse coordinate             */
  int gc, gcb;                /* Global sparse coordinate            */
  int cc, c1, n;              /* Counters                            */
  int back;                   /* Cell in a backward direction from c */
  int *map;                   /* Local backward map                  */
  int *gmap;                  /* Global backward map                 */

  /* Get the auxiliary mapping direction                             */
  map = window->zm1;
  gmap = NULL;

  /* Get the surface map. This is simply a copy of the 2D wet cells  */
  /* to process vector and is only included here for consistency     */
  /* since it is overwritten each time-step.                         */
  for (cc = 1; cc <= avec2D; cc++) {
    sur[cc] = vec2D[cc];
  }

  /* Get the bottom map. This is located via the global bottom maps. */
  for (cc = 1; cc <= nvec2D; cc++)
    bot[cc] = 0;
  c1 = 1;
  for (cc = 1; cc <= nvec; cc++) {
    c = gc = vec[cc];
    if (geom->fm[c].wn == wn) {
      /* back is a land cell => use the original c                   */
      if (c == geom->zp1[c]) c = gc;
 
      bot[c1] = geom->fm[c].sc;

      c1++;
    }
  }

  /* Add auxiliary location bottom coordinates                       */
  for (cc = nvec2D + 1; cc <= avec2D; cc++) {
    c = vec2D[cc];              /* Local cell to process             */
    back = map[c];              /* Backward map of local cell        */
    gc = window->wsa[c];        /* Global cell corresponding to c    */
    gcb = window->wsa[back];    /* Global cell corresponding to back */

    /* Bottom coordinate for auxiliary cells                         */
    while (geom->fm[gcb].wn != 0 && geom->fm[gc].wn != 0) {
      c = window->zm1[c];
      back = map[c];
      gc = window->wsa[c];
      gcb = window->wsa[back];
    }
    bot[cc] = c;
  }
}


void surfbot_builde(geometry_t *geom,   /* Global geometry */
		   geometry_t *window,  /* Window data structure */
                   int wn,      /* Window number */
                   int *sur,    /* Surface cell vector */
                   int *bot,    /* Bottom cell vector */
                   int *vec,    /* Global bottom cells vector */
                   int nvec,    /* Size of vec */
                   int *vec2D,  /* 2D cells to process vector */
                   int nvec2D,  /* Number of wet cells in vec2D */
                   int avec2D   /* Wet + auxiliary cells in vec2D */
  )
{
  int c, e, es, ee, e1;         /* Local sparse coordinate */
  int gc, ge, le;               /* Global sparse coordinate */
  int cc, c1, c2, n, j, k;      /* Counters */
  int *wgst;
  int *e2ee;

  /*-----------------------------------------------------------------*/
  /* Get the surface map. This is simply a copy of the 2D wet cells  */
  /* to process vector and is only included here for consistency     */
  /* since it is overwritten each time-step.                         */
  e2ee = i_alloc_1d(window->szeS);
  for (ee = 1; ee <= window->n2_e1; ee++) {
    e = window->w2_e1[ee];
    window->sur_e1[ee] = e;
    window->bot_e1[ee] = 0;
    e2ee[e] = ee;
  }

  /*-----------------------------------------------------------------*/
  /* Get the bottom map. This is located via the global bottom maps. */
  for (ee = 1; ee <= geom->b2_e1; ee++) {
    ge = geom->bot_e1[ee];

    if (geom->fm[ge].we == wn) {
      es = geom->m2de[ge];
      e1 = e2ee[geom->g2we[wn][es]];
      e = geom->g2we[wn][ge];
      window->bot_e1[e1] = e;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Add auxiliary location bottom coordinates                       */
  for (ee = nvec2D + 1; ee <= avec2D; ee++) {
    e = vec2D[ee];              /* Local edge to process             */
    /* Bottom coordinate for auxiliary cells */
    while (window->zm1e[e] != e)
      e = window->zm1e[e];
    bot[ee] = window->zp1e[e];

    /* Bottom coordinate for boundary cells adjacent to OUTSIDE */
    /* cells in the domain interior.  */
  }

  /*-----------------------------------------------------------------*/
  /* Get the bottom coordinate via the shallowest cell centres       */
  /* adjacent to the edge. This is how bottom edges are defined in   */
  /* build sparse_grid_us().                                         */
  wgst = i_alloc_1d(window->szcS);
  memset(wgst, 0, window->szcS * sizeof(int));
  for (cc = window->a2_t+1; cc <= window->n2_t; cc++) {
    c = window->w2_t[cc];
    wgst[c] = c;
  }
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no2_e1; ee++) {
      c = open->obc_e2[ee];
      c1 = open->ogc_t[ee];
      wgst[c1] = c;
    }
  }
  memset(window->s2k, 0, window->szc * sizeof(int));
  for (cc = 1; cc <= window->n2_t; cc++) {
    c = window->w2_t[cc];
    window->s2k[c] = window->nz-1;
    while (c != window->zm1[c]) {
      c1 = c;
      c = window->zm1[c];
      window->s2k[c] = window->s2k[c1] - 1;
    }
  }
  for (ee = window->v2_e1 + 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    j = window->e2e[e][0];
    c2 = window->e2c[e][1];

    /* Set any ghost cells to wet cells */
    if (wgst[c1]) {
      c1 = c2;
      j = window->e2e[e][1];
    }
    if (wgst[c2]) c2 = c1;

    /* Find the shallowest cell */
    while (c1 != window->zm1[c1]) {
      c1 = window->zm1[c1];
    }
    c1 = window->zp1[c1];

    while (c2 != window->zm1[c2]) {
      c2 = window->zm1[c2];
    }
    c2 = window->zp1[c2];

    if (window->s2k[c2] > window->s2k[c1]) {
      c1 = c2;
      j = window->e2e[e][1];
    }
    window->bot_e1[ee] = window->c2e[j][c1];
  }

  i_free_1d(wgst);
}

/* END surfbot_build()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to build a 2D vector from the global maps                 */
/*-------------------------------------------------------------------*/
int *win_vector_build(geometry_t *geom,   /* Global geometry         */
		      geometry_t *window, /* Local geometry          */
		      int wn,             /* Window number           */
		      int *gv,            /* Global 2D vector        */
		      int mode            /* Cell center,e1,e2 code  */
		      )
{
  int c;              /* Local sparse coordinate                     */
  int cg;             /* Global sparse coordinate                    */
  int cc;             /* Counters                                    */
  int *vec;           /* Local cells to process vector               */
  int nvec;           /* Size of vec                                 */
  int *lv;            /* Local 2D vector to fill                     */

  /* Allocate memory                                                 */
  lv = i_alloc_1d(window->enonS + 1);
  memset(lv, 0, (window->enonS + 1) * sizeof(int));

  /* Get the zoom code and auxiliary mapping direction */
  if (mode == 0) {
    vec = window->w2_t;
    nvec = window->b2_t;
  } else if (mode == 1) {
    vec = window->w2_e1;
    nvec = window->b2_e1;
  } else {
    vec = window->w2_e2;
    nvec = window->b2_e2;
  }

  /* Fill the window vector from the global vector                   */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cg = window->wsa[c];
    if (geom->fm[cg].wn == wn) {
      lv[c] = gv[cg];
    }
  }
  return(lv);
}

/* END win_vector_build()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to reorder the global to local map. Any wet cells in the  */
/* mask (i.e. any cells with a window number the same as the global  */
/* to local window map, fm.wn) are assigned the first local sparse   */
/* locations. Any wet cells not in the mask are assigned local       */
/* locations from position msize onwards. This allows for auxiliary  */
/* cells in the surface layer to follow the wet surface cells so     */
/* the first msize locations in the local array correspond to the    */
/* surface wet cells followed by surface auxiliary cells, allowing   */
/* 2D maps to be extracted from 3D maps.                             */
/* NOTE : this assumes that wsa is filled with sparse locations in   */
/* the (k,i,j) direction and mask is filled in the (i,j) direction.  */
/*-------------------------------------------------------------------*/
void reorder_gl_map(geometry_t *geom, /* Global geometry             */
		    int wn,     /* Window number                     */
                    int *wsa,   /* Local work sparse array           */
                    int wsize,  /* Size of wsa                       */
                    int *mask,  /* Global locations to order first   */
                    int msize   /* Size of mask                      */
  )
{
  int c, cc, c1, c2;            /* Counters                          */

  /* Reset the wet cells in the mask to the first locations in the   */
  /* local sparse array. If mask is the surface layer then this is   */
  /* not necessary since the surface layer wet cells already occupy  */
  /* the first local sparse locations of fm.sc.                      */
  c1 = 1;
  for (cc = 1; cc <= msize; cc++) {
    c = mask[cc];
    if (geom->fm[c].wn == wn) {
      geom->fm[c].sc = c1;
      c1++;
    }
  }

  /* Reset the wet cells not in the mask to the local locations from */
  /* msize onwards.                                                  */
  c2 = 1;
  c1 = msize + 1;
  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];
    if (geom->fm[c].wn == wn && c != mask[c2]) {
      geom->fm[c].sc = c1;
      c1++;
      c2++;
    }
  }
}

/* END reorder_gl_map()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set a mask for cell identification                                */
/*-------------------------------------------------------------------*/
void set_mask(geometry_t *window)
{
  int cc, c, ee, e, ci, n;

  window->cask = i_alloc_1d(window->szc);
  memset(window->cask, 0, window->szc * sizeof(int));
  for (cc = 1; cc <= window->v3_t; cc++) {
    c = window->w3_t[cc];
    window->cask[c] = W_WET;
  }
  for (cc = 1; cc <= window->v2_t; cc++) {
    c = window->w2_t[cc];
    window->cask[c] |= W_SURF;
    c = window->bot_t[cc];
    window->cask[c] |= W_BOT;
    c = window->zm1[c];
    window->cask[c] |= W_SED;
  }
  for (cc = window->v3_t+1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    window->cask[c] |= W_AUX;
  }
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    window->cask[c] = W_GST;
    c = window->bin[cc];
    window->cask[c] |= W_INT;
  }
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      window->cask[c] = W_NOBC;
    }
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      window->cask[c] = W_GOBC;
    }
  }

  window->eask = i_alloc_1d(window->sze);
  memset(window->eask, 0, window->sze * sizeof(int));
  for (ee = 1; ee <= window->v3_e1; ee++) {
    e = window->w3_e1[ee];
    window->eask[e] = W_WET;
  }
  for (ee = 1; ee <= window->v2_e1; ee++) {
    e = window->w2_e1[ee];
    window->eask[e] |= W_SURF;
    e = window->bot_e1[ee];
    window->eask[e] |= W_BOT;
    e = window->zm1e[e];
    window->eask[e] |= W_SED;
  }
  for (ee = window->v3_e1+1; ee <= window->a3_e1; ee++) {
    e = window->w3_e1[ee];
    window->eask[e] |= W_AUX;
  }
  for (ee = 1; ee <= window->nbpte1; ee++) {
    e = window->bpte1[ee];
    window->eask[e] = W_GST;
    e = window->bine1[ee];
    window->eask[e] |= W_INT;
  }
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      window->eask[e] = W_NOBC;
    }
    for (ee = open->no3_e1 + 1; ee <= open->to3_e1; ee++) {
      e = open->obc_e1[ee];
      window->eask[e] = W_TOBC;
    }
  }
}

/* END set_mask()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise the window data structure                   */
/*-------------------------------------------------------------------*/
window_t **win_data_build(master_t *master,   /* Master data         */
                          geometry_t **window /* Window data         */
  )
{
  window_t **windat;            /* Window data                       */
  int nwindows;                 /* Number of windows                 */
  int n, j, tn, k = 0;          /* Counters                          */
  int c, cs, cc;                /* Cell locations                    */
  int e, es, ee;                /* Edge locations                    */

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  nwindows = master->nwindows;
  windat = (window_t **)malloc(sizeof(window_t *) * (nwindows + 1));

  for (n = 1; n <= nwindows; n++) {
    windat[n] = win_data_init(master, window[n]);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise window data                                          */
  for (n = 1; n <= nwindows; n++) {

    /* Constants                                                     */
    windat[n]->dt = master->grid_dt;
    windat[n]->iratio = master->iratio;
    windat[n]->t = master->t;
    /*windat[n]->etarlxtc = master->etarlxtc;*/
    if (!(master->cfl & NONE))
      windat[n]->mcfl2d = windat[n]->mcfl3d = HUGE;

    /* Set the 3D water column tracers pointers                      */
    windat[n]->sal = windat[n]->temp = NULL;
    windat[n]->tke = windat[n]->diss = windat[n]->L = windat[n]->omega =
      NULL;
    windat[n]->Q2 = windat[n]->Q2L = windat[n]->Kq = NULL;
    windat[n]->u1m = windat[n]->u2m = windat[n]->wm = windat[n]->Kzm = NULL;
    windat[n]->tempm = windat[n]->saltm = windat[n]->tram = NULL;
    windat[n]->u1_adv = windat[n]->u1_hdif = windat[n]->u1_vdif = NULL;
    windat[n]->u1_cor = windat[n]->u1_btp = windat[n]->u1_bcp = NULL;
    windat[n]->u2_adv = windat[n]->u2_hdif = windat[n]->u2_vdif = NULL;
    windat[n]->u1_sto = windat[n]->u2_sto = NULL;
    windat[n]->fluxe1 = windat[n]->fluxe2 = windat[n]->mom_bal = NULL;
    windat[n]->fluxw = windat[n]->fluxkz = NULL;
    windat[n]->brunt = windat[n]->int_wave = NULL;
    windat[n]->rich_gr = windat[n]->rich_fl = NULL;
    windat[n]->reynolds = windat[n]->froude = windat[n]->sigma_t = NULL;
    windat[n]->rossby_in = windat[n]->sound = windat[n]->schan = NULL;
    windat[n]->shear_v = windat[n]->b_prod = windat[n]->s_prod = NULL;
    windat[n]->speed_3d = windat[n]->perc = windat[n]->energy = NULL;
    windat[n]->dum1 = windat[n]->dum2 = windat[n]->dum3 = NULL;
    windat[n]->vcorr = windat[n]->acorr = windat[n]->kenergy = NULL;
    windat[n]->regionid = windat[n]->regres = windat[n]->Vi = NULL;
    windat[n]->reefe1 = windat[n]->reefe2 = windat[n]->agetr = NULL;
    windat[n]->tr_adv = windat[n]->tr_hdif = windat[n]->tr_vdif = windat[n]->tr_ncon = NULL;
    windat[n]->wave_stke1 = windat[n]->wave_stke1 = NULL;

    for (tn = 0; tn < windat[n]->ntr; tn++) {
      if (strcmp("salt", master->trname[tn]) == 0) {
        windat[n]->sal = windat[n]->tr_wc[tn];
        windat[n]->sno = tn;
      } else if (strcmp("temp", master->trname[tn]) == 0) {
        windat[n]->temp = windat[n]->tr_wc[tn];
        windat[n]->tno = tn;
      } else if (strcmp("diss", master->trname[tn]) == 0)
        windat[n]->diss = windat[n]->tr_wc[tn];
      else if (strcmp("tke", master->trname[tn]) == 0)
        windat[n]->tke = windat[n]->tr_wc[tn];
      else if (strcmp("omega", master->trname[tn]) == 0)
        windat[n]->omega = windat[n]->tr_wc[tn];
      else if (strcmp("tki", master->trname[tn]) == 0)
        windat[n]->Q2 = windat[n]->tr_wc[tn];
      else if (strcmp("tki_l", master->trname[tn]) == 0)
        windat[n]->Q2L = windat[n]->tr_wc[tn];
      else if (strcmp("Kq", master->trname[tn]) == 0)
        windat[n]->Kq = windat[n]->tr_wc[tn];
      else if (strcmp("lscale", master->trname[tn]) == 0)
        windat[n]->L = windat[n]->tr_wc[tn];
      else if (strcmp("u1mean", master->trname[tn]) == 0)
        windat[n]->u1m = windat[n]->tr_wc[tn];
      else if (strcmp("u2mean", master->trname[tn]) == 0)
        windat[n]->u2m = windat[n]->tr_wc[tn];
      else if (strcmp("wmean", master->trname[tn]) == 0)
        windat[n]->wm = windat[n]->tr_wc[tn];
      else if (strcmp("temp_mean", master->trname[tn]) == 0)
        windat[n]->tempm = windat[n]->tr_wc[tn];
      else if (strcmp("salt_mean", master->trname[tn]) == 0)
        windat[n]->saltm = windat[n]->tr_wc[tn];
      else if (strcmp("tracer_mean", master->trname[tn]) == 0)
        windat[n]->tram = windat[n]->tr_wc[tn];
      else if (strcmp("Kzmean", master->trname[tn]) == 0)
        windat[n]->Kzm = windat[n]->tr_wc[tn];
      else if (strcmp("flux_e1", master->trname[tn]) == 0)
        windat[n]->fluxe1 = windat[n]->tr_wc[tn];
      else if (strcmp("flux_e2", master->trname[tn]) == 0)
        windat[n]->fluxe2 = windat[n]->tr_wc[tn];
      else if (strcmp("flux_w", master->trname[tn]) == 0)
        windat[n]->fluxw = windat[n]->tr_wc[tn];
      else if (strcmp("flux_kz", master->trname[tn]) == 0)
        windat[n]->fluxkz = windat[n]->tr_wc[tn];
      else if (strcmp("smagorinsky", master->trname[tn]) == 0)
        windat[n]->sdc = windat[n]->tr_wc[tn];
      else if (strcmp("mom_balance", master->trname[tn]) == 0) {
        windat[n]->mom_bal = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_adv", master->trname[tn]) == 0) {
        windat[n]->u1_adv = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_cor", master->trname[tn]) == 0) {
        windat[n]->u1_cor = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_sto", master->trname[tn]) == 0) {
        windat[n]->u1_sto = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_hdif", master->trname[tn]) == 0) {
        windat[n]->u1_hdif = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_vdif", master->trname[tn]) == 0) {
        windat[n]->u1_vdif = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_btp", master->trname[tn]) == 0) {
        windat[n]->u1_btp = windat[n]->tr_wc[tn];
      } else if (strcmp("u1_bcp", master->trname[tn]) == 0) {
        windat[n]->u1_bcp = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_adv", master->trname[tn]) == 0) {
        windat[n]->u2_adv = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_cor", master->trname[tn]) == 0) {
        windat[n]->u2_cor = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_sto", master->trname[tn]) == 0) {
        windat[n]->u2_sto = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_hdif", master->trname[tn]) == 0) {
        windat[n]->u2_hdif = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_vdif", master->trname[tn]) == 0) {
        windat[n]->u2_vdif = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_btp", master->trname[tn]) == 0) {
        windat[n]->u2_btp = windat[n]->tr_wc[tn];
      } else if (strcmp("u2_bcp", master->trname[tn]) == 0) {
        windat[n]->u2_bcp = windat[n]->tr_wc[tn];
      } else if (strcmp("wave_stke1", master->trname[tn]) == 0) {
        windat[n]->wave_stke1 = windat[n]->tr_wc[tn];
      } else if (strcmp("wave_stke2", master->trname[tn]) == 0) {
        windat[n]->wave_stke2 = windat[n]->tr_wc[tn];
      } else if (strcmp("brunt_vaisala", master->trname[tn]) == 0) {
        windat[n]->brunt = windat[n]->tr_wc[tn];
      } else if (strcmp("int_wave_speed", master->trname[tn]) == 0) {
        windat[n]->int_wave = windat[n]->tr_wc[tn];
      } else if (strcmp("richardson_gr", master->trname[tn]) == 0) {
        windat[n]->rich_gr = windat[n]->tr_wc[tn];
      } else if (strcmp("richardson_fl", master->trname[tn]) == 0) {
        windat[n]->rich_fl = windat[n]->tr_wc[tn];
      } else if (strcmp("froude", master->trname[tn]) == 0) {
        windat[n]->froude = windat[n]->tr_wc[tn];
      } else if (strcmp("sigma_t", master->trname[tn]) == 0) {
        windat[n]->sigma_t = windat[n]->tr_wc[tn];
      } else if (strcmp("energy", master->trname[tn]) == 0) {
        windat[n]->energy = windat[n]->tr_wc[tn];
      } else if (strcmp("kenergy", master->trname[tn]) == 0) {
        windat[n]->kenergy = windat[n]->tr_wc[tn];
      } else if (strcmp("sound", master->trname[tn]) == 0) {
        windat[n]->sound = windat[n]->tr_wc[tn];
      } else if (strcmp("sound_channel", master->trname[tn]) == 0) {
        windat[n]->schan = windat[n]->tr_wc[tn];
      } else if (strcmp("tracer1", master->trname[tn]) == 0) {
        windat[n]->dum1 = windat[n]->tr_wc[tn];
      } else if (strcmp("tracer2", master->trname[tn]) == 0) {
        windat[n]->dum2 = windat[n]->tr_wc[tn];
      } else if (strcmp("tracer3", master->trname[tn]) == 0) {
        windat[n]->dum3 = windat[n]->tr_wc[tn];
      } else if (strcmp("Vcorr", master->trname[tn]) == 0) {
        windat[n]->vcorr = windat[n]->tr_wc[tn];
      } else if (strcmp("Acorr", master->trname[tn]) == 0) {
        windat[n]->acorr = windat[n]->tr_wc[tn];
      } else if (strcmp("Vi", master->trname[tn]) == 0) {
        windat[n]->Vi = windat[n]->tr_wc[tn];
      } else if (strcmp("flow_salt", master->trname[tn]) == 0) {
        windat[n]->riversalt = windat[n]->tr_wc[tn];
      } else if (strcmp("reef_fraction_e1", master->trname[tn]) == 0) {
        windat[n]->reefe1 = windat[n]->tr_wc[tn];
      } else if (strcmp("reef_fraction_e2", master->trname[tn]) == 0) {
        windat[n]->reefe2 = windat[n]->tr_wc[tn];
      } else if (strcmp("reynolds", master->trname[tn]) == 0) {
        windat[n]->reynolds = windat[n]->tr_wc[tn];
      } else if (strcmp("shear_vert", master->trname[tn]) == 0) {
        windat[n]->shear_v = windat[n]->tr_wc[tn];
      } else if (strcmp("buoy_prod", master->trname[tn]) == 0) {
        windat[n]->b_prod = windat[n]->tr_wc[tn];
      } else if (strcmp("shear_prod", master->trname[tn]) == 0) {
        windat[n]->s_prod = windat[n]->tr_wc[tn];
      } else if (strcmp("rossby_internal", master->trname[tn]) == 0) {
        windat[n]->rossby_in = windat[n]->tr_wc[tn];
      } else if (strcmp("current_speed_3d", master->trname[tn]) == 0) {
        windat[n]->speed_3d = windat[n]->tr_wc[tn];
      } else if (strcmp("regionid", master->trname[tn]) == 0) {
        windat[n]->regionid = windat[n]->tr_wc[tn];
      } else if (strcmp("residence", master->trname[tn]) == 0) {
        windat[n]->regres = windat[n]->tr_wc[tn];
      } else if (strcmp("rtemp", master->trname[tn]) == 0) {
        windat[n]->rtemp = windat[n]->tr_wc[tn];
      } else if (strcmp("rsalt", master->trname[tn]) == 0) {
        windat[n]->rsalt = windat[n]->tr_wc[tn];
      } else if (strcmp("temp_tc", master->trname[tn]) == 0) {
        windat[n]->temp_tc = windat[n]->tr_wc[tn];
      } else if (strcmp("salt_tc", master->trname[tn]) == 0) {
        windat[n]->salt_tc = windat[n]->tr_wc[tn];
      } else if (strcmp("flush", master->trname[tn]) == 0) {
        windat[n]->fltr = windat[n]->tr_wc[tn];
      } else if (strcmp("age", master->trname[tn]) == 0) {
        windat[n]->agetr = windat[n]->tr_wc[tn];
      } else if (strcmp("decorr_e1", master->trname[tn]) == 0) {
        windat[n]->decv1 = windat[n]->tr_wc[tn];
      } else if (master->swr_type & SWR_3D && strcmp("swr_attenuation", master->trname[tn]) == 0) {
        windat[n]->swr_attn = windat[n]->tr_wc[tn];
      } else if (master->trperc >= 0) {
	char buf[MAXSTRLEN];
	sprintf(buf, "percentile_%s", 
		master->trinfo_3d[master->trperc].name);
	if (strcmp(buf, master->trname[tn]) == 0) {
	  windat[n]->perc = windat[n]->tr_wc[tn];
	}
      } else if (master->trtend >= 0) {
	char buf[MAXSTRLEN];
	sprintf(buf, "%s_adv", master->trinfo_3d[master->trtend].name);
	if (strcmp(buf, master->trname[tn]) == 0)
	  windat[n]->tr_adv = windat[n]->tr_wc[tn];
	sprintf(buf, "%s_hdif", master->trinfo_3d[master->trtend].name);
	if (strcmp(buf, master->trname[tn]) == 0)
	  windat[n]->tr_hdif = windat[n]->tr_wc[tn];
	sprintf(buf, "%s_vdif", master->trinfo_3d[master->trtend].name);
	if (strcmp(buf, master->trname[tn]) == 0)
	  windat[n]->tr_vdif = windat[n]->tr_wc[tn];
	sprintf(buf, "%s_ncon", master->trinfo_3d[master->trtend].name);
	if (strcmp(buf, master->trname[tn]) == 0)
	  windat[n]->tr_ncon = windat[n]->tr_wc[tn];
      }
    }

    /* Initialise all 3D and 2D variables required by the window     */
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];

      for (tn = 0; tn < windat[n]->ntr; tn++) {
        windat[n]->tr_wc[tn][cc] = master->tr_wc[tn][c];
      }

      windat[n]->w[cc] = master->w[c];
      windat[n]->Kz[cc] = master->Kz[c];
      windat[n]->Vz[cc] = master->Vz[c];

      if (cc <= window[n]->enonS) {
        windat[n]->eta[cc] = master->eta[c];
        windat[n]->detadt[cc] = master->detadt[c];
        windat[n]->wtop[cc] = master->wtop[c];
        windat[n]->wbot[cc] = master->wbot[c];
        windat[n]->patm[cc] = master->patm[c];

        if (master->ntrS) {
	  for (tn = 0; tn < windat[n]->ntrS; tn++) {
	    windat[n]->tr_wcS[tn][cc] = master->tr_wcS[tn][c];
	  }
	}
        if (master->nsed) {
          for (k = 0; k < window[n]->sednz; k++)
            for (tn = 0; tn < windat[n]->nsed; tn++) {
              windat[n]->tr_sed[tn][k][cc] = master->tr_sed[tn][k][c];
            }
        }
      }
    }

    /* Edge arrays                                                   */
    for (ee = 1; ee <= window[n]->n3_e1; ee++) {
      e = window[n]->wse[ee];
      windat[n]->u1[ee] = master->u1[e];
      if (ee <= window[n]->n2_e1) {
	windat[n]->u1av[ee] = master->u1av[e];
        windat[n]->wind1[ee] = master->wind1[e];
        windat[n]->windspeed[ee] = master->windspeed[e];
        windat[n]->winddir[ee] = master->winddir[e];
      }
    }

    /* Save the bottom velocity for the 2D mode                      */
    for (ee = 1; ee <= window[n]->b2_e1; ee++) {
      es = window[n]->w2_e1[ee];  /* 2D coordiate corresponding to e */
      e = window[n]->bot_e1[ee];  /* 3D bottom coordinate            */
      windat[n]->u1bot[es] = windat[n]->u1[e] - windat[n]->u1av[es];
    }

    window[n]->windat = windat[n];
  }

  if (master->show_win) {
    for (cc = 1; cc <= geom->a2_t; cc++) {
      c = geom->wsa[cc];
      n = geom->fm[c].wn;
      cs = geom->fm[c].sc;
      master->shwin[c] = windat[n]->shwin[cs] = (double)n;
    }
  }

  return (windat);
}

/* END win_data_build()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to initialise the window data structure                   */
/*-------------------------------------------------------------------*/
window_t *win_data_init(master_t *master,   /* Master data structure */
                        geometry_t *window  /* Window data structure */
			)
{
  window_t *windat;             /* Window data structure             */
  int winsize;                  /* Number of sparse cells in window  */
  int n, m;

  windat = win_data_alloc();

  /* Multiple windows. Allocate memory for each window               */
  if (master->nwindows > 1) {

    /*---------------------------------------------------------------*/
    /* 3D arrays                                                     */
    /* Cell centred                                                  */
    winsize = window->szc;
    windat->ntr = master->ntr;
    windat->tr_wc = d_alloc_2d(winsize, windat->ntr);
    windat->w = d_alloc_1d(winsize);
    windat->Kz = d_alloc_1d(winsize);
    windat->Vz = d_alloc_1d(winsize);
    windat->dens = d_alloc_1d(winsize);
    windat->dens_0 = d_alloc_1d(winsize);
    windat->waterss = d_alloc_1d(winsize);
    windat->nrvorc = d_alloc_1d(winsize);
    windat->kec = d_alloc_1d(winsize);
    windat->div = d_alloc_1d(winsize);
    windat->u = d_alloc_1d(winsize);
    windat->v = d_alloc_1d(winsize);
    /* Edges                                                         */
    winsize = window->sze;
    windat->u1 = d_alloc_1d(winsize);
    windat->u2 = d_alloc_1d(winsize);
    windat->dzu1 = d_alloc_1d(winsize);
    windat->u1flux3d = d_alloc_1d(winsize);
    windat->u1b = d_alloc_1d(winsize);
    windat->u2b = d_alloc_1d(winsize);
    windat->nrvore = d_alloc_1d(winsize);
    windat->npvore = d_alloc_1d(winsize);
    /* Vertices                                                      */
    winsize = window->szv;
    windat->fv = d_alloc_1d(winsize);
    windat->circ = d_alloc_1d(winsize);
    windat->rvor = d_alloc_1d(winsize);
    windat->nrvor = d_alloc_1d(winsize);
    windat->npvor = d_alloc_1d(winsize);
    /* Others                                                        */
    if (master->velrlx & RELAX) {
      relax_info_t *relax = master->vel_rlx;
      windat->vel_rlx = relax_info_init(relax->rlxn, relax->rlxtc, 
					relax->rlxdt, winsize, winsize);
      windat->vel_rlx->rlx = relax->rlx;
    }

    /*---------------------------------------------------------------*/
    /* 2D arrays                                                     */
    windat->ntrS = master->ntrS;
    /* Cell centred                                                  */
    winsize = window->szcS;
    if (windat->ntrS)
      windat->tr_wcS = d_alloc_2d(winsize, windat->ntrS);
    windat->eta = d_alloc_1d(winsize);
    windat->topz = d_alloc_1d(winsize);
    windat->wdiff2d = d_alloc_1d(winsize);
    windat->detadt = d_alloc_1d(winsize);
    windat->wtop = d_alloc_1d(winsize);
    windat->wbot = d_alloc_1d(winsize);
    windat->patm = d_alloc_1d(winsize);
    windat->waterss2d = d_alloc_1d(winsize);
    windat->etab = d_alloc_1d(winsize);
    windat->uav = d_alloc_1d(winsize);
    windat->vav = d_alloc_1d(winsize);
    if (!(master->means & NONE)) {
      windat->meanc = d_alloc_1d(winsize);
      if (master->means & TIDAL) {
	windat->odeta = d_alloc_1d(winsize);
      }
      if (master->means & VEL3D) {
	windat->ume = d_alloc_1d(window->sze);
	memset(windat->ume, 0, window->sze * sizeof(double));
      }
      if (master->means & VEL2D) {
	windat->uame = d_alloc_1d(window->szeS);
	memset(windat->uame, 0, window->szeS * sizeof(double));
      }
      if (master->means & VOLFLUX || master->tmode & SP_U1VM) {
	windat->u1vm = d_alloc_1d(window->sze);
	memset(windat->u1vm, 0, window->sze * sizeof(double));
      }
    }
    if (master->swr)
      windat->swr = d_alloc_1d(winsize);
    if (master->heatflux & (NET_HEAT | ADVANCED | INVERSE | COMP_HEAT | COMP_HEAT_MOM)) {
      windat->heatf = d_alloc_1d(winsize);
      if (master->heatflux & ADVANCED)
	windat->cloud = d_alloc_1d(winsize);
      if (master->heatflux & COMP_HEAT)
	windat->lwrd = d_alloc_1d(winsize);
    }
    if (master->heatflux & (COMP_HEAT_MOM | COMP_HEAT_NONE))
      if (master->swrd) windat->swrd = d_alloc_1d(winsize);
    if (master->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE)) {
      // The if's are needed for comp_heat_none
      if (master->lwro) windat->lwro = d_alloc_1d(winsize);
      if (master->lhfd) windat->lhfd = d_alloc_1d(winsize);
      if (master->shfd) windat->shfd = d_alloc_1d(winsize);
    }

    /* Allocate salt flux arrays, if needed                          */
    if (master->airtemp) windat->airtemp = d_alloc_1d(winsize);
    if (master->precip)  windat->precip  = d_alloc_1d(winsize);
    if (master->evap)    windat->evap    = d_alloc_1d(winsize);
    
    if (master->sh_f & (DEWPOINT|WETBULB))
      windat->wetb = d_alloc_1d(winsize);
    else
      if (master->sh_f & RELHUM)
	windat->rh = d_alloc_1d(winsize);
    
    if (master->albedo_l)
      windat->light = d_alloc_1d(winsize);

    if (master->heatflux & SURF_RELAX || master->hftemp)
      windat->hftemp = d_alloc_1d(winsize);

    /* Edges                                                         */
    winsize = window->szeS;
    windat->u1av = d_alloc_1d(winsize);
    windat->u2av = d_alloc_1d(winsize);
    windat->nu1av = d_alloc_1d(winsize);
    windat->u1flux = d_alloc_1d(winsize);
    windat->depth_e1 = d_alloc_1d(winsize);
    windat->u1bot = d_alloc_1d(winsize);
    windat->wind1 = d_alloc_1d(winsize);
    windat->windspeed = d_alloc_1d(winsize);
    windat->winddir   = d_alloc_1d(winsize);
    windat->u1avb = d_alloc_1d(winsize);
    windat->nsed = master->nsed;
    if (windat->nsed)
      windat->tr_sed = d_alloc_3d(winsize, window->sednz, windat->nsed);

    /* Vertices                                                      */

    /* Other                                                         */
    if (master->etarlx & (RELAX|ALERT|BOUNDARY)) {
      relax_info_t *relax = master->eta_rlx;
      windat->eta_rlx = relax_info_init(relax->rlxn, relax->rlxtc, 
					relax->rlxdt, winsize, 0);
      if ((m = tracer_find_index("oeta", master->ntrS, master->trinfo_2d)) >= 0)
	windat->eta_rlx->val1 = windat->tr_wcS[m];
      windat->eta_rlx->rlx = relax->rlx;
      windat->eta_rlx->dv0 = relax->dv0;
      windat->eta_rlx->dv1 = relax->dv1;
      windat->eta_rlx->tc0 = relax->tc0;
      windat->eta_rlx->tc1 = relax->tc1;
      windat->eta_rlx->slope = relax->slope;
    }

    if (master->ntr) windat->trinc = i_alloc_1d(master->ntr);
    if (master->ntrS) windat->trincS = i_alloc_1d(master->ntrS);

    for (m = 0; m < windat->ntrS; m++) {
      if (strcmp("cfl2d", master->trinfo_2d[m].name) == 0)
        windat->cfl2d = windat->tr_wcS[m];
      if (strcmp("cfl3d", master->trinfo_2d[m].name) == 0)
        windat->cfl3d = windat->tr_wcS[m];
      if (strcmp("courant", master->trinfo_2d[m].name) == 0)
        windat->cour = windat->tr_wcS[m];
      if (strcmp("lipschitz", master->trinfo_2d[m].name) == 0)
        windat->lips = windat->tr_wcS[m];
      if (strcmp("diffstab", master->trinfo_2d[m].name) == 0)
        windat->ahsb = windat->tr_wcS[m];
      if (strcmp("mixed_layer", master->trinfo_2d[m].name) == 0)
        windat->mixl = windat->tr_wcS[m];
      if (strcmp("steric", master->trinfo_2d[m].name) == 0)
        windat->steric = windat->tr_wcS[m];
      if (strcmp("abs_vor", master->trinfo_2d[m].name) == 0)
        windat->av = windat->tr_wcS[m];
      if (strcmp("rel_vor", master->trinfo_2d[m].name) == 0)
        windat->rv = windat->tr_wcS[m];
      if (strcmp("pot_vor", master->trinfo_2d[m].name) == 0)
        windat->pv = windat->tr_wcS[m];
      if (strcmp("rv_drvdt", master->trinfo_2d[m].name) == 0)
        windat->rv_drvdt = windat->tr_wcS[m];
      if (strcmp("rv_nonlin", master->trinfo_2d[m].name) == 0)
        windat->rv_nonlin = windat->tr_wcS[m];
      if (strcmp("rv_beta", master->trinfo_2d[m].name) == 0)
        windat->rv_beta = windat->tr_wcS[m];
      if (strcmp("rv_strch", master->trinfo_2d[m].name) == 0)
        windat->rv_strch = windat->tr_wcS[m];
      if (strcmp("rv_jebar", master->trinfo_2d[m].name) == 0)
        windat->rv_jebar = windat->tr_wcS[m];
      if (strcmp("rv_wsc", master->trinfo_2d[m].name) == 0)
        windat->rv_wsc = windat->tr_wcS[m];
      if (strcmp("rv_bsc", master->trinfo_2d[m].name) == 0)
        windat->rv_bsc = windat->tr_wcS[m];
      if (strcmp("w1mean", master->trinfo_2d[m].name) == 0)
        windat->w1m = windat->tr_wcS[m];
      if (strcmp("w2mean", master->trinfo_2d[m].name) == 0)
        windat->w2m = windat->tr_wcS[m];
      if (strcmp("eta_mean", master->trinfo_2d[m].name) == 0)
        windat->etam = windat->tr_wcS[m];
      if (strcmp("u1av_mean", master->trinfo_2d[m].name) == 0)
        windat->u1am = windat->tr_wcS[m];
      if (strcmp("u2av_mean", master->trinfo_2d[m].name) == 0)
        windat->u2am = windat->tr_wcS[m];
      if (strcmp("nhf", master->trinfo_2d[m].name) == 0)
        windat->nhfd = windat->tr_wcS[m];
      if (strcmp("swr", master->trinfo_2d[m].name) == 0) {
	if (master->params->heatflux & (COMP_HEAT_MOM | COMP_HEAT_NONE)) {
	  /* See logic in heatflux.c:comp_heat_mom                   */
	  windat->swr = windat->tr_wcS[m];
	} else 
	  windat->swrd = windat->tr_wcS[m];
      }

      /* These heatflux tracer numbers reflect logic in load_tracer  */
      if (strcmp("lwr", master->trinfo_2d[m].name) == 0) {
	if (master->params->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE))
	  windat->lwrn = m;
	else
	  windat->lwrd = windat->tr_wcS[m];
      }
      if (strcmp("lhf", master->trinfo_2d[m].name) == 0) {
	if (master->params->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE))
	  windat->lhfn = m;
	else
	  windat->lhfd = windat->tr_wcS[m];
      }
      if (strcmp("shf", master->trinfo_2d[m].name) == 0) {
	if (master->params->heatflux & (COMP_HEAT | COMP_HEAT_MOM | COMP_HEAT_NONE))
	  windat->shfn = m;
	else 
	  windat->shfd = windat->tr_wcS[m];
      }
      if (strcmp("precip", master->trinfo_2d[m].name) == 0) {
	if ( (master->params->heatflux & COMP_HEAT_NONE) ||
	     params->saltflux & (ADVANCED | ORIGINAL) )
	  windat->precipn = m;
      }
      if (strcmp("evap", master->trinfo_2d[m].name) == 0) {
	if ( (master->params->heatflux & COMP_HEAT_NONE) ||
	     params->saltflux & (ADVANCED | ORIGINAL) )
	  windat->evapn = m;
      }
      if (strcmp("nsf", master->trinfo_2d[m].name) == 0)
        windat->nsfd = windat->tr_wcS[m];
      if (strcmp("ustrcw", master->trinfo_2d[m].name) == 0)
        windat->ustrcw = windat->tr_wcS[m];
      if (strcmp("wave_ub", master->trinfo_2d[m].name) == 0)
        windat->wave_ub = windat->tr_wcS[m];
      if (strcmp("wave_period", master->trinfo_2d[m].name) == 0)
        windat->wave_period = windat->tr_wcS[m];
      if (strcmp("wave_amp", master->trinfo_2d[m].name) == 0)
        windat->wave_amp = windat->tr_wcS[m];
      if (strcmp("wave_dir", master->trinfo_2d[m].name) == 0)
        windat->wave_dir = windat->tr_wcS[m];
      if (strcmp("wave_Sxy", master->trinfo_2d[m].name) == 0)
        windat->wave_Sxy = windat->tr_wcS[m];
      if (strcmp("wave_Syx", master->trinfo_2d[m].name) == 0)
        windat->wave_Syx = windat->tr_wcS[m];
      if (strcmp("wave_Fx", master->trinfo_2d[m].name) == 0)
        windat->wave_Fx = windat->tr_wcS[m];
      if (strcmp("wave_Fy", master->trinfo_2d[m].name) == 0)
        windat->wave_Fy = windat->tr_wcS[m];
      if (strcmp("wave_ste1", master->trinfo_2d[m].name) == 0)
        windat->wave_ste1 = windat->tr_wcS[m];
      if (strcmp("wave_ste2", master->trinfo_2d[m].name) == 0)
        windat->wave_ste2 = windat->tr_wcS[m];
      if (strcmp("tau_w1", master->trinfo_2d[m].name) == 0)
        windat->tau_w1 = windat->tr_wcS[m];
      if (strcmp("tau_w2", master->trinfo_2d[m].name) == 0)
        windat->tau_w2 = windat->tr_wcS[m];
      if (strcmp("tau_diss1", master->trinfo_2d[m].name) == 0)
        windat->tau_diss1 = windat->tr_wcS[m];
      if (strcmp("tau_diss2", master->trinfo_2d[m].name) == 0)
        windat->tau_diss2 = windat->tr_wcS[m];
      if (strcmp("wave_Cd", master->trinfo_2d[m].name) == 0)
        windat->wave_Cd = windat->tr_wcS[m];
      if (strcmp("rossby_external", master->trinfo_2d[m].name) == 0)
        windat->rossby_ex = windat->tr_wcS[m];
      if (strcmp("current_speed_2d", master->trinfo_2d[m].name) == 0)
        windat->speed_2d = windat->tr_wcS[m];
      if (strcmp("speed_sq", master->trinfo_2d[m].name) == 0)
        windat->speed_sq = windat->tr_wcS[m];
      if (strcmp("obc_phase", master->trinfo_2d[m].name) == 0)
        windat->obc_phase = windat->tr_wcS[m];
      if (strcmp("wind_Cd", master->trinfo_2d[m].name) == 0)
        windat->wind_Cd = windat->tr_wcS[m];
      if (strcmp("AVHRR", master->trinfo_2d[m].name) == 0)
        windat->avhrr = windat->tr_wcS[m];
      if (strcmp("ghrsst", master->trinfo_2d[m].name) == 0)
        windat->ghrsst = windat->tr_wcS[m];
      if (strcmp("u1_rad", master->trinfo_2d[m].name) == 0)
        windat->u1_rad = windat->tr_wcS[m];
      if (strcmp("u2_rad", master->trinfo_2d[m].name) == 0)
        windat->u2_rad = windat->tr_wcS[m];
      if (strcmp("windows", master->trinfo_2d[m].name) == 0)
        windat->shwin = windat->tr_wcS[m];
      if (strcmp("alerts_actual", master->trinfo_2d[m].name) == 0)
        windat->alert_a = windat->tr_wcS[m];
      if (strcmp("alerts_cumulative", master->trinfo_2d[m].name) == 0)
        windat->alert_c = windat->tr_wcS[m];
      if (strcmp("U1VH0", master->trinfo_2d[m].name) == 0)
        windat->u1vhin = windat->tr_wcS[m];
      if (strcmp("U2VH0", master->trinfo_2d[m].name) == 0)
        windat->u2vhin = windat->tr_wcS[m];
      if (strcmp("vol_cons", master->trinfo_2d[m].name) == 0)
        windat->vol_cons = windat->tr_wcS[m];
      if (strcmp("sonic_depth", master->trinfo_2d[m].name) == 0)
        windat->sonic = windat->tr_wcS[m];
      if (strcmp("wet_cells", master->trinfo_2d[m].name) == 0)
        windat->wetcell = windat->tr_wcS[m];
      if (strcmp("surf_layer", master->trinfo_2d[m].name) == 0)
        windat->surfz = windat->tr_wcS[m];
      if (strcmp("surf_slope", master->trinfo_2d[m].name) == 0)
        windat->slope_x = windat->tr_wcS[m];
      if (strcmp("tau_be1", master->trinfo_2d[m].name) == 0)
        windat->tau_be1 = windat->tr_wcS[m];
      if (strcmp("tau_be2", master->trinfo_2d[m].name) == 0)
        windat->tau_be2 = windat->tr_wcS[m];
      if (strcmp("tau_bm", master->trinfo_2d[m].name) == 0)
        windat->tau_bm = windat->tr_wcS[m];
      if (master->swr_type & SWR_2D && strcmp("swr_attenuation", master->trinfo_2d[m].name) == 0)
        windat->swr_attn = windat->tr_wcS[m];
      if (strcmp("swr_deep_attenuation", master->trinfo_2d[m].name) == 0)
        windat->swr_attn1 = windat->tr_wcS[m];
      if (strcmp("swr_transmission", master->trinfo_2d[m].name) == 0)
        windat->swr_tran = windat->tr_wcS[m];
      if (strcmp("swr_bot_absorb", master->trinfo_2d[m].name) == 0)
        windat->swr_babs = windat->tr_wcS[m];
      if (strcmp("flow", master->trinfo_2d[m].name) == 0)
        windat->riverflow = windat->tr_wcS[m];
      if (strcmp("flow_depth", master->trinfo_2d[m].name) == 0)
        windat->riverdepth = windat->tr_wcS[m];
      if (strcmp("eta_tc", master->trinfo_2d[m].name) == 0)
        windat->eta_tc = windat->tr_wcS[m];
      if (strcmp("eta_inc", master->trinfo_2d[m].name) == 0)
        windat->eta_inc = windat->tr_wcS[m];
      if (strcmp("sed_error", master->trinfo_2d[m].name) == 0)
        windat->sederr = windat->tr_wcS[m];
      if (strcmp("eco_error", master->trinfo_2d[m].name) == 0)
        windat->ecoerr = windat->tr_wcS[m];
      if (strcmp("decorr_e1", master->trinfo_2d[m].name) == 0)
        windat->decv1 = windat->tr_wcS[m];
      /*if (strcmp("oeta", master->trinfo_2d[m].name) == 0 && windat->eta_rlx)
        windat->eta_rlx->val1 = windat->tr_wcS[m];*/
    }
  }

  /* Single window. The master and the window are the same - point   */
  /* the window to the master data structure arrays.                 */
  else {
    /* 3D arrays                                                     */
    windat->ntr = master->ntr;
    windat->tr_wc = master->tr_wc;
    windat->nsed = master->nsed;
    if (windat->nsed)
      windat->tr_sed = master->tr_sed;
    windat->u1 = master->u1;
    windat->u2 = master->u2;
    windat->u = master->u;
    windat->v = master->v;
    windat->w = master->w;
    windat->Kz = master->Kz;
    windat->Vz = master->Vz;
    windat->dzu1 = master->dzu1;
    windat->u1flux3d = master->u1flux3d;
    windat->dens = master->dens;
    windat->dens_0 = master->dens_0;
    windat->u1b = master->u1b;
    windat->u2b = master->u2b;
    windat->waterss = master->waterss;
    if (master->velrlx & RELAX) 
      windat->vel_rlx = master->vel_rlx;
    if (master->means & VEL3D)
      windat->ume = master->ume;
    if (master->means & VOLFLUX || master->tmode & SP_U1VM)
      windat->u1vm = master->u1vm;

    /* 2D arrays                                                     */
    windat->ntrS = master->ntrS;
    if (windat->ntrS)
      windat->tr_wcS = master->tr_wcS;
    windat->eta = master->eta;
    windat->u1av = master->u1av;
    windat->u2av = master->u2av;
    windat->nu1av = master->nu1av;
    windat->u1flux = master->u1flux;
    windat->uav = master->uav;
    windat->vav = master->vav;
    windat->depth_e1 = master->depth_e1;
    windat->u1bot = master->u1bot;
    windat->topz = master->topz;
    windat->wdiff2d = master->wdiff2d;
    windat->detadt = master->detadt;
    windat->wtop = master->wtop;
    windat->wbot = master->wbot;
    windat->wind1 = master->wind1;
    windat->windspeed = master->windspeed;
    windat->winddir = master->winddir;
    windat->patm = master->patm;
    windat->airtemp = master->airtemp;
    windat->cloud = master->cloud;
    windat->precip = master->precip;
    windat->evap = master->evap;
    windat->rh = master->rh;
    windat->wetb = master->wetb;
    windat->waterss2d = master->waterss2d;
    /*windat->eta_rlx = master->eta_rlx;*/

    windat->etab = master->etab;
    windat->u1avb = master->u1avb;

    windat->heatf = master->heatf;
    windat->swr = master->swr;
    windat->light = master->light;
    windat->hftemp = master->hftemp;

    windat->nrvore = master->nrvore;
    windat->npvore = master->npvore;
    windat->nrvorc = master->nrvorc;
    windat->fv = master->fv;

    windat->circ = master->circ;
    windat->kec = master->kec;
    windat->div = master->div;
    windat->rvor = master->rvor;
    windat->nrvor = master->nrvor;
    windat->npvor = master->npvor;

    windat->trinc = master->trinc;
    windat->trincS = master->trincS;
    windat->cfl2d = master->cfl2d;
    windat->cfl3d = master->cfl3d;
    windat->cour = master->cour;
    windat->lips = master->lips;
    windat->ahsb = master->ahsb;
    windat->steric = master->steric;
    windat->av = master->av;
    windat->rv = master->rv;
    windat->pv = master->pv;
    windat->rv_drvdt = master->rv_drvdt;
    windat->rv_nonlin = master->rv_nonlin;
    windat->rv_beta = master->rv_beta;
    windat->rv_strch = master->rv_strch;
    windat->rv_jebar = master->rv_jebar;
    windat->rv_wsc = master->rv_wsc;
    windat->rv_bsc = master->rv_bsc;
    windat->rossby_ex = master->rossby_ex;
    windat->speed_2d = master->speed_2d;
    windat->speed_sq = master->speed_sq;
    windat->obc_phase = master->obc_phase;
    windat->wind_Cd = master->wind_Cd;
    windat->etam = master->etam;
    windat->u1am = master->u1am;
    windat->u2am = master->u2am;
    windat->w1m = master->w1m;
    windat->w2m = master->w2m;
    windat->mixl = master->mixl;
    windat->nhfd = master->nhfd;
    windat->swrd = master->swrd;
    windat->lwrd = master->lwrd;
    windat->lwro = master->lwro;
    windat->shfd = master->shfd;
    windat->lhfd = master->lhfd;
    windat->nsfd = master->nsfd;
    windat->meanc = master->meanc;
    windat->odeta = master->odeta;
    windat->avhrr = master->avhrr;
    windat->ghrsst = master->ghrsst;
    windat->shwin = master->shwin;
    windat->alert_a = master->alert_a;
    windat->alert_c = master->alert_c;
    windat->u1vhin = master->u1vhin;
    windat->u2vhin = master->u2vhin;
    windat->ustrcw = master->ustrcw;
    windat->wave_ub = master->wave_ub;
    windat->wave_period = master->wave_period;
    windat->wave_amp = master->wave_amp;
    windat->wave_dir = master->wave_dir;
    windat->wave_Sxy = master->wave_Sxy;
    windat->wave_Syx = master->wave_Syx;
    windat->wave_Fx = master->wave_Fx;
    windat->wave_Fy = master->wave_Fy;
    windat->wave_ste1 = master->wave_ste1;
    windat->wave_ste2 = master->wave_ste2;
    windat->wave_Cd = master->wave_Cd;
    windat->u1_rad = master->u1_rad;
    windat->u2_rad = master->u2_rad;
    windat->vol_cons = master->vol_cons;
    windat->regionid = master->regionid;
    windat->regres = master->regres;
    windat->sonic = master->sonic;
    windat->wetcell = master->wetcell;
    windat->surfz = master->surfz;
    windat->slope_x = master->slope_x;
    windat->tau_be1 = master->tau_be1;
    windat->tau_be2 = master->tau_be2;
    windat->tau_bm = master->tau_bm;
    windat->swr_attn = master->swr_attn;
    windat->swr_attn1 = master->swr_attn1;
    windat->swr_tran = master->swr_tran;
    windat->swr_babs = master->swr_babs;
    windat->riverflow = master->riverflow;
    windat->riverdepth = master->riverdepth;
    windat->eta_tc = master->eta_tc;
    windat->sederr = master->sederr;
    windat->ecoerr = master->ecoerr;
    if (master->etarlx & (RELAX|ALERT|BOUNDARY)) 
      windat->eta_rlx = master->eta_rlx;

    /* Diagnostic indicies                                           */
    windat->precipn = master->precipn;
    windat->evapn   = master->evapn;
    windat->lwrn    = master->lwrn;
    windat->shfn    = master->shfn;
    windat->lhfn    = master->lhfn;
  }
  winsize = window->enon + 1;
  if (windat->mom_bal)
    memset(windat->mom_bal, 0, winsize * sizeof(double));
  if (windat->u1_adv)
    memset(windat->u1_adv, 0, winsize * sizeof(double));
  if (windat->u1_hdif)
    memset(windat->u1_hdif, 0, winsize * sizeof(double));
  if (windat->u1_vdif)
    memset(windat->u1_vdif, 0, winsize * sizeof(double));
  if (windat->u1_btp)
    memset(windat->u1_btp, 0, winsize * sizeof(double));
  if (windat->u1_bcp)
    memset(windat->u1_bcp, 0, winsize * sizeof(double));
  if (windat->u1_cor)
    memset(windat->u1_cor, 0, winsize * sizeof(double));
  if (windat->u2_adv)
    memset(windat->u2_adv, 0, winsize * sizeof(double));
  if (windat->u2_hdif)
    memset(windat->u2_hdif, 0, winsize * sizeof(double));
  if (windat->u2_vdif)
    memset(windat->u2_vdif, 0, winsize * sizeof(double));
  if (windat->u2_btp)
    memset(windat->u2_btp, 0, winsize * sizeof(double));
  if (windat->u2_bcp)
    memset(windat->u2_bcp, 0, winsize * sizeof(double));
  if (windat->u2_cor)
    memset(windat->u2_cor, 0, winsize * sizeof(double));
  if (windat->tr_adv)
    memset(windat->tr_adv, 0, winsize * sizeof(double));
  if (windat->tr_hdif)
    memset(windat->tr_hdif, 0, winsize * sizeof(double));
  if (windat->tr_vdif)
    memset(windat->tr_vdif, 0, winsize * sizeof(double));
  windat->sur_e1 = i_alloc_1d(window->szeS);
  winsize = window->enonS + 1;
  if (windat->sederr)
    memset(windat->sederr, 0, winsize * sizeof(double));
  if (windat->ecoerr)
    memset(windat->ecoerr, 0, winsize * sizeof(double));

  windat->vd = NULL;
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (m = 0; m < windat->ntr; m++) {
      if (open->bcond_tra[m] & DESCAL && windat->vd == NULL)
	windat->vd = d_alloc_1d(window->nz + 1);
    }
  }
  return (windat);
}

/* END win_data_init()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to clear memory for the window data structure             */
/*-------------------------------------------------------------------*/
void win_data_clear(window_t *windat  /* Window data structure       */
  )
{
  return;
  /* 3D arrays                                                       */
  d_free_2d(windat->tr_wc);
  d_free_1d(windat->u1);
  d_free_1d(windat->u2);
  d_free_1d(windat->u);
  d_free_1d(windat->v);
  d_free_1d(windat->w);
  d_free_1d(windat->Kz);
  d_free_1d(windat->Vz);
  d_free_1d(windat->dzu1);
  d_free_1d(windat->u1flux3d);
  d_free_1d(windat->dens);
  d_free_1d(windat->dens_0);
  d_free_1d(windat->u1b);
  d_free_1d(windat->u2b);
  d_free_1d(windat->waterss);

  /* 2D arrays                                                       */
  if (windat->ntrS)
    d_free_2d(windat->tr_wcS);
  d_free_1d(windat->eta);
  d_free_1d(windat->u1av);
  d_free_1d(windat->u2av);
  d_free_1d(windat->uav);
  d_free_1d(windat->vav);
  d_free_1d(windat->etab);
  d_free_1d(windat->u1avb);
  i_free_1d(windat->sur_e1);
  d_free_1d(windat->depth_e1);
  d_free_1d(windat->u1flux);
  d_free_1d(windat->u1bot);
  d_free_1d(windat->topz);
  d_free_1d(windat->wdiff2d);
  d_free_1d(windat->detadt);
  d_free_1d(windat->wtop);
  d_free_1d(windat->wbot);
  d_free_1d(windat->patm);
  if (windat->airtemp)
    d_free_1d(windat->airtemp);
  if (windat->cloud)
    d_free_1d(windat->cloud);
  d_free_1d(windat->wind1);
  if(windat->windspeed != NULL)
    d_free_1d(windat->windspeed);
  if(windat->winddir != NULL)
    d_free_1d(windat->winddir);
  d_free_1d(windat->waterss2d);
  if (windat->nalert)
    i_free_1d(windat->nalert);
  if (windat->ef)
    d_free_1d(windat->ef);
  if (windat->vf)
    d_free_1d(windat->vf);
  if (windat->ceta)
    i_free_1d(windat->ceta);
  if (windat->cu1)
    i_free_1d(windat->cu1);
  if (windat->cu1a)
    i_free_1d(windat->cu1a);
  if (windat->tempb)
    d_free_1d(windat->tempb);
  if (windat->salb)
    d_free_1d(windat->salb);
  if(windat->totid)
    i_free_1d(windat->totid);
  if(windat->trtot)
    d_free_1d(windat->trtot);
  if(windat->trinc)
    i_free_1d(windat->trinc);
  if(windat->trincS)
    i_free_1d(windat->trincS);
  if(windat->vd)
    d_free_1d(windat->vd);
  if(windat->ume)
    d_free_1d(windat->ume);
  if(windat->u1vm)
    d_free_1d(windat->u1vm);
  if(windat->uame)
    d_free_1d(windat->uame);

  /* Space for these variables is deallocated by the master via the  */
  /* scheduler cleanup for 1 window.                                 */
  if (master->nwindows > 1) {
    if (windat->heatf)
      d_free_1d(windat->heatf);
    if (windat->light)
      d_free_1d(windat->light);
    if (windat->hftemp)
      d_free_1d(windat->hftemp);
    if (windat->airtemp)
      d_free_1d(windat->airtemp);
    if (windat->cloud)
      d_free_1d(windat->cloud);
  }
}

/* END win_data_clear()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to initialise the window data structure                   */
/*-------------------------------------------------------------------*/
win_priv_t **win_consts_init(master_t *master,    /* Master data     */
                             geometry_t **window  /* Window data     */
  )
{
  geometry_t *geom = master->geom;  /* Global geometry               */
  win_priv_t **wincon;          /* Window constants structure        */
  window_t *windat;             /* Window data structure             */
  int nwindows;                 /* Number of windows                 */
  int n, cc, tn;                /* Counters                          */
  int c, cs, c1, cb;            /* Sparse locations                  */
  int szm;                      /* Max 3D cell, edge, vertex size    */
  int szmS;                     /* Max 2D cell, edge, vertex size    */
  int winsize;                  /* Number of sparse cells in window  */
  int i, j, k;                  /* Cartesian locations               */
  int ee, e;

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the wincon arrays                           */
  nwindows = master->nwindows;
  wincon = (win_priv_t **)malloc(sizeof(win_priv_t *) * (nwindows + 1));
  for (n = 1; n <= nwindows; n++) {
    windat = window[n]->windat;
    wincon[n] = win_consts_alloc();
    szm = window[n]->szm;
    szmS = window[n]->szmS;
    wincon[n]->tstart = master->tstart;
    wincon[n]->ntr = master->ntr;
    wincon[n]->ntrS = master->ntrS;
    wincon[n]->nsed = master->nsed;
    wincon[n]->ntdif_h = master->ntdif_h;
    wincon[n]->ntdif_v = master->ntdif_v;
    wincon[n]->rampf = master->rampf;
    wincon[n]->calc_closure = master->calc_closure;
    wincon[n]->s_func = master->s_func;
    wincon[n]->tgrid = NONE;
    wincon[n]->togn = master->togn;
    wincon[n]->sflux = i_alloc_1d(master->ntr);

    /* 3D work arrays                                               */
    winsize = window[n]->sgsiz = window[n]->enon + 1;
    wincon[n]->w1 = d_alloc_1d(szm);
    wincon[n]->w2 = d_alloc_1d(szm);
    wincon[n]->w3 = d_alloc_1d(szm);
    wincon[n]->w4 = d_alloc_1d(szm);
    /* Also allocate the parallel versions, if required              */
#ifdef HAVE_OMP
    wincon[n]->trans_num_omp = master->trans_num_omp;
    if (wincon[n]->trans_num_omp > 1) {
      wincon[n]->w1n = d_alloc_2d(szm, wincon[n]->trans_num_omp);
      wincon[n]->w2n = d_alloc_2d(szm, wincon[n]->trans_num_omp);
      wincon[n]->w3n = d_alloc_2d(szm, wincon[n]->trans_num_omp);
      wincon[n]->w4n = d_alloc_2d(szm, wincon[n]->trans_num_omp);
    }
#endif
    wincon[n]->w5 = d_alloc_1d(szm);
    wincon[n]->w6 = d_alloc_1d(szm);
    wincon[n]->w7 = d_alloc_1d(szm);
    wincon[n]->w8 = d_alloc_1d(szm);
    wincon[n]->w9 = d_alloc_1d(szm);
    wincon[n]->w10 = d_alloc_1d(szm);
    wincon[n]->rdens = d_alloc_1d(winsize);
    wincon[n]->s1 = i_alloc_1d(szm);
    wincon[n]->s2 = i_alloc_1d(szm);
    wincon[n]->s3 = i_alloc_1d(szm);
    wincon[n]->s4 = i_alloc_1d(szm);
    wincon[n]->c1 = c_alloc_1d(szm);
    wincon[n]->ba = d_alloc_1d(window[n]->npem+1);
    wincon[n]->gmap = i_alloc_1d(szm);
    memset(wincon[n]->s1, 0, szm * sizeof(int));
    memset(wincon[n]->s2, 0, szm * sizeof(int));
    if (windat->u1_adv || windat->u1_hdif || windat->u1_vdif ||
        windat->u1_cor || windat->u1_btp || windat->u1_bcp ||
        windat->u2_adv || windat->u2_hdif || windat->u2_vdif ||
        windat->u2_cor || windat->u2_btp || windat->u2_bcp ||
	windat->tr_adv || windat->tr_hdif || windat->tr_vdif ||
	windat->tr_ncon)
      wincon[n]->tendency = d_alloc_1d(szm);
    if(master->trasc & (LAGRANGE|FFSL) || 
       master->momsc & LAGRANGE ||
       master->do_pt) {
      wincon[n]->m2d = i_alloc_1d(winsize);
      wincon[n]->dzz = d_alloc_1d(szm);
      wincon[n]->cellz = d_alloc_1d(szm);
      wincon[n]->s5 = i_alloc_1d(szm);
    }
    if(master->trasc & FFSL || master->do_pt) {
      wincon[n]->tr_mod = d_alloc_1d(szm);
      wincon[n]->tr_mod_x = d_alloc_1d(szm);
      wincon[n]->tr_mod_y = d_alloc_1d(szm);
      wincon[n]->tr_mod_z = d_alloc_1d(szm);
    }
    if(master->trasc & FCT) {
      wincon[n]->Fxh = d_alloc_1d(window[n]->sze);
      wincon[n]->Fzh = d_alloc_1d(window[n]->szc);
      wincon[n]->Ax = d_alloc_1d(window[n]->sze);
      wincon[n]->Az = d_alloc_1d(window[n]->szc);
    }
    if(master->trasc & HIORDER) {
      wincon[n]->trp = d_alloc_1d(window[n]->sze);
      wincon[n]->trm = d_alloc_1d(window[n]->sze);
    }
    wincon[n]->nu = d_alloc_1d(winsize);
    wincon[n]->nv = d_alloc_1d(winsize);
    wincon[n]->nw = d_alloc_1d(winsize);
    wincon[n]->clxf = i_alloc_1d(winsize);
    wincon[n]->clyf = i_alloc_1d(winsize);
    wincon[n]->clzf = i_alloc_1d(winsize);
    wincon[n]->clxc = i_alloc_1d(winsize);
    wincon[n]->clyc = i_alloc_1d(winsize);
    wincon[n]->clzc = i_alloc_1d(winsize);
    wincon[n]->crfxf = d_alloc_1d(winsize);
    wincon[n]->crfyf = d_alloc_1d(winsize);
    wincon[n]->crfzf = d_alloc_1d(winsize);
    wincon[n]->crfxc = d_alloc_1d(winsize);
    wincon[n]->crfyc = d_alloc_1d(winsize);
    wincon[n]->crfzc = d_alloc_1d(winsize);
    /* Also allocate the parallel versions, if required              */
#ifdef HAVE_OMP
    if (wincon[n]->trans_num_omp > 1) {
      wincon[n]->tr_modn = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_x = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_y = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_z = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
    }
#endif
    /* 2D work arrays                                                */
    winsize = window[n]->sgsizS = window[n]->enonS + 1;
    wincon[n]->d1 = d_alloc_1d(szmS);
    wincon[n]->d2 = d_alloc_1d(szmS);
    wincon[n]->d3 = d_alloc_1d(szmS);
    wincon[n]->d4 = d_alloc_1d(szmS);
    wincon[n]->d5 = d_alloc_1d(szmS);
    wincon[n]->d6 = d_alloc_1d(szmS);
    wincon[n]->d7 = d_alloc_1d(szmS);
    wincon[n]->c2 = c_alloc_1d(szmS);
    if (master->thin_merge) {
      wincon[n]->kth_e1 = i_alloc_1d(szmS);
      wincon[n]->kth_e2 = i_alloc_1d(szmS);
    }
    wincon[n]->cdry_e1 = i_alloc_1d(szmS);
    wincon[n]->cdry_e2 = i_alloc_1d(szmS);
    wincon[n]->cbot_e1 = i_alloc_1d(szmS);
    wincon[n]->cbot_e2 = i_alloc_1d(szmS);
    wincon[n]->i1 = i_alloc_1d(szmS);
    wincon[n]->i2 = i_alloc_1d(szmS);
    wincon[n]->i3 = i_alloc_1d(szmS);
    wincon[n]->i4 = i_alloc_1d(szmS);
    wincon[n]->i5 = i_alloc_1d(szmS);
    wincon[n]->i6 = i_alloc_1d(szmS);
    wincon[n]->i7 = i_alloc_1d(szmS);
    if (master->eta_rlx || master->etarlx & INCREMENT)
      wincon[n]->eta_rlx3d = d_alloc_1d(winsize);
    if (master->fetch)
      wincon[n]->fetch = d_alloc_2d(8, winsize);

    /*
     * Error handling, strickly speaking these don't need to be
     * winsize, just the number of parallel threads but this way the
     * indexing is more straightforward
     */
    wincon[n]->gint_error  = c_alloc_2d(MAXSTRLEN, winsize);
    wincon[n]->gint_errorf = i_alloc_1d(winsize);

    /* 1D work arrays                                                */
    winsize = window[n]->nz + 1;
    wincon[n]->v1 = d_alloc_1d(winsize);
    wincon[n]->v2 = d_alloc_1d(winsize);
    wincon[n]->v3 = d_alloc_1d(winsize);
    wincon[n]->v4 = d_alloc_1d(winsize);
    wincon[n]->v5 = d_alloc_1d(winsize);
    wincon[n]->v6 = d_alloc_1d(winsize);
    wincon[n]->v7 = d_alloc_1d(winsize);
    wincon[n]->v8 = d_alloc_1d(winsize);
    wincon[n]->v9 = d_alloc_1d(winsize);
    wincon[n]->v10 = d_alloc_1d(winsize);
    wincon[n]->v11 = d_alloc_1d(winsize);
    wincon[n]->v12 = d_alloc_1d(winsize);
    wincon[n]->tmass = d_alloc_1d(master->ntr + 3);
    wincon[n]->tsmass = d_alloc_1d(master->ntr);

    if(window[n]->sednz) {
      winsize = window[n]->sednz + 1;
      wincon[n]->sd1 = d_alloc_1d(winsize);
    }

    /* Sigma arrays                                                  */
    wincon[n]->mdx = d_alloc_1d(window[n]->szeS);
    
    /* Multiple windows. Allocate memory for each window             */
    if (master->nwindows > 1) {
      /* 3D arrays                                                   */
      int szc = window[n]->szc;
      int sze = window[n]->sze;
      int szm = window[n]->szm;
      wincon[n]->dz = d_alloc_1d(szc);
      wincon[n]->t11 = d_alloc_1d(szm);
      wincon[n]->t12 = d_alloc_1d(szm);
      wincon[n]->t22 = d_alloc_1d(szm);
      wincon[n]->q = (master->q != NULL) ? d_alloc_1d(szc) : NULL;
      wincon[n]->u1vh = d_alloc_1d(sze);

      if (master->smagcode & (U1_AK|U1_SAK))
        wincon[n]->u1kh = d_alloc_1d(szc);
      else
	wincon[n]->u1kh = windat->sdc;


      /* 2D arrays                                                   */
      szc = window[n]->szcS;
      sze = window[n]->szeS;
      szm = window[n]->szmS;
      wincon[n]->oldeta = d_alloc_1d(szc);
      wincon[n]->one = d_alloc_1d(window[n]->szm);
      for (c = 1; c < window[n]->szm; c++)
        wincon[n]->one[c] = 1.0;
      if (master->sigma) {
        wincon[n]->Ds = d_alloc_1d(szc);
        wincon[n]->Hs = d_alloc_1d(szc);
        wincon[n]->Hn1 = d_alloc_1d(sze);
        wincon[n]->Hn2 = d_alloc_1d(sze);
      }
      wincon[n]->Cd = d_alloc_1d(szc);
      wincon[n]->z0 = d_alloc_1d(szc);
      wincon[n]->topdensu1 = d_alloc_1d(sze);
      wincon[n]->densavu1 = d_alloc_1d(sze);
      wincon[n]->u1c1 = d_alloc_1d(sze);
      wincon[n]->u1c3 = d_alloc_1d(sze);
      wincon[n]->u1c4 = d_alloc_1d(sze);
      wincon[n]->u1c5 = d_alloc_1d(sze);
      wincon[n]->u1c6 = d_alloc_1d(sze);
      wincon[n]->u1adv = d_alloc_1d(sze);
      wincon[n]->u1inter = d_alloc_1d(sze);
      wincon[n]->coriolis = d_alloc_1d(szc);
    } else {
      /* Single window. The master and the window are the same -     */
      /* point the window to the master data structure arrays.       */
      /* 3D arrays */
      wincon[n]->dz = master->dz;
      wincon[n]->t11 = master->t11;
      wincon[n]->t12 = master->t12;
      wincon[n]->t22 = master->t22;
      wincon[n]->q = master->q;
      wincon[n]->u1vh = master->u1vh;
      wincon[n]->u1kh = master->u1kh;

      /* 2D arrays                                                   */
      wincon[n]->oldeta = master->oldeta;
      wincon[n]->one = master->one;
      if (master->sigma) {
        wincon[n]->Ds = master->Ds;
        wincon[n]->Hs = master->Hs;
        wincon[n]->Hn1 = master->Hn1;
        wincon[n]->Hn2 = master->Hn2;
      }
      wincon[n]->Cd = master->Cd;
      wincon[n]->z0 = master->z0;
      wincon[n]->topdensu1 = master->topdensu1;
      wincon[n]->densavu1 = master->densavu1;
      wincon[n]->u1c1 = master->u1c1;
      wincon[n]->u1c3 = master->u1c3;
      wincon[n]->u1c4 = master->u1c4;
      wincon[n]->u1c5 = master->u1c5;
      wincon[n]->u1c6 = master->u1c6;
      wincon[n]->u1adv = master->u1adv;
      wincon[n]->u1inter = master->u1inter;
      wincon[n]->coriolis = master->coriolis;
    }

    /*---------------------------------------------------------------*/
    /* Initialise the wincon arrays                                  */
    wincon[n]->g = master->g;
    wincon[n]->tz = tm_tz_offset(master->timeunit);
    wincon[n]->ambpress = master->ambpress;
    wincon[n]->trasc = master->trasc;
    wincon[n]->trasf = 0;
    wincon[n]->momsc = master->momsc;
    wincon[n]->momsc2d = master->momsc2d;
    wincon[n]->hmin = master->hmin;
    wincon[n]->uf = master->uf;
    wincon[n]->quad_bfc = master->quad_bfc;
    wincon[n]->slip = master->slip;
    wincon[n]->etamax = master->etamax;
    wincon[n]->velmax = master->velmax;
    wincon[n]->velmax2d = master->velmax2d;
    wincon[n]->etadiff = master->etadiff;
    wincon[n]->wmax = master->wmax;
    wincon[n]->ultimate = master->ultimate;
    wincon[n]->smagorinsky = master->smagorinsky;
    wincon[n]->sue1 = master->sue1;
    wincon[n]->kue1 = master->kue1;
    wincon[n]->bsue1 = master->bsue1;
    wincon[n]->bkue1 = master->bkue1;
    wincon[n]->diff_scale = master->diff_scale;
    wincon[n]->smag_smooth = master->smag_smooth;
    wincon[n]->visc_method = master->visc_method;
    wincon[n]->stab = master->stab;
    wincon[n]->cfl = master->cfl;
    wincon[n]->cfl_dt = master->cfl_dt;
    wincon[n]->lnm = master->lnm;
    wincon[n]->vorticity = master->vorticity;
    wincon[n]->numbers = master->numbers;
    wincon[n]->compatible = master->compatible;
    wincon[n]->filter = master->filter;
    wincon[n]->porusplate = master->porusplate;
    wincon[n]->totals = master->totals;
    wincon[n]->robust = master->robust;
    wincon[n]->fatal = master->fatal;
    wincon[n]->mode2d = master->mode2d;
    wincon[n]->means = master->means;
    wincon[n]->means_dt = master->means_dt;
    wincon[n]->means_next = master->means_next;
    wincon[n]->tendf = master->tendf;
    wincon[n]->trtend = master->trtend;
    wincon[n]->trsplit = master->trsplit;
    wincon[n]->trflux = master->trflux;
    wincon[n]->trfd1 = master->trfd1;
    wincon[n]->trfd2 = master->trfd2;
    wincon[n]->trperc = master->trperc;
    wincon[n]->trflsh = master->trflsh;
    strcpy(wincon[n]->trage, master->trage);
    wincon[n]->mixlayer = master->mixlayer;
    wincon[n]->thin_merge = master->thin_merge;
    wincon[n]->sigma = master->sigma;
    wincon[n]->nonlinear = master->nonlinear;
    wincon[n]->calc_dens = master->calc_dens;
    wincon[n]->heatflux = master->heatflux;
    wincon[n]->sh_f = master->sh_f;
    wincon[n]->albedo = master->albedo;
    wincon[n]->saltflux = master->saltflux;
    wincon[n]->hftc = master->hftc;
    wincon[n]->bulkf = master->bulkf;
    wincon[n]->zref = master->zref;
    wincon[n]->rampstart = master->rampstart;
    wincon[n]->rampend = master->rampend;
    wincon[n]->hf_ramp = master->hf_ramp;
    wincon[n]->u1vh0 = master->u1vh0;
    wincon[n]->u2vh0 = master->u2vh0;
    wincon[n]->u1kh0 = master->u1kh0;
    wincon[n]->u2kh0 = master->u2kh0;
    wincon[n]->hmean1 = master->hmean1;
    wincon[n]->hmean2 = master->hmean2;
    wincon[n]->amean = master->amean;
    wincon[n]->zs = master->zs;
    wincon[n]->vz0 = master->vz0;
    wincon[n]->kz0 = master->kz0;
    wincon[n]->min_tke = master->min_tke;
    wincon[n]->min_diss = master->min_diss;
    wincon[n]->Lmin = master->Lmin;
    wincon[n]->smooth_VzKz = master->smooth_VzKz;
    wincon[n]->eparam = master->eparam;
    wincon[n]->kz_alpha = master->kz_alpha;
    wincon[n]->vz_alpha = master->vz_alpha;
    wincon[n]->fcf = master->fcf;
    wincon[n]->wave_alpha = master->wave_alpha;
    wincon[n]->wave_b1 = master->wave_b1;
    wincon[n]->wave_hf = master->wave_hf;
    wincon[n]->smagcode = master->smagcode;
    wincon[n]->u1_f = master->u1_f;
    wincon[n]->u1av_f = master->u1av_f;
    wincon[n]->save_force = master->save_force;
    wincon[n]->vinit = master->vinit;
    wincon[n]->etarlx = master->etarlx;
    wincon[n]->velrlx = master->velrlx;
    wincon[n]->alertf = master->alertf;
    wincon[n]->u1_f = master->u1_f;
    wincon[n]->u1av_f = master->u1av_f;
    wincon[n]->eta_f = master->eta_f;
    wincon[n]->vel2d_f = master->vel2d_f;
    wincon[n]->vel3d_f = master->vel3d_f;
    wincon[n]->wvel_f = master->wvel_f;
    wincon[n]->tend_f = master->tend_f;
    wincon[n]->div2d_f = master->div2d_f;
    wincon[n]->div3d_f = master->div3d_f;
    wincon[n]->cfl_f = master->cfl_f;
    wincon[n]->ts_f = master->ts_f;
    wincon[n]->shear_f = master->shear_f;
    wincon[n]->hdiff_f = master->hdiff_f;
    wincon[n]->tide_r = master->tide_r;
    wincon[n]->amax = master->amax;
    wincon[n]->hmax = master->hmax;
    wincon[n]->vmax = master->vmax;
    wincon[n]->btmax = master->btmax;
    wincon[n]->bcmax = master->bcmax;
    wincon[n]->cmax = master->cmax;
    wincon[n]->detamax = master->detamax;
    wincon[n]->dwmax = master->dwmax;
    wincon[n]->dtmax = master->dtmax;
    wincon[n]->dsmax = master->dsmax;
    wincon[n]->smax = master->smax;
    wincon[n]->exmapf = master->exmapf;
    wincon[n]->waves = master->waves;
    wincon[n]->orbital = master->orbital;
    wincon[n]->fillf = master->fillf;
    wincon[n]->pssinput = master->pssinput;
    wincon[n]->conserve = master->conserve;
    wincon[n]->do_closure = master->do_closure;
    wincon[n]->do_pt = master->do_pt;
    wincon[n]->tmode = master->tmode;
    wincon[n]->trout = master->trout;
    wincon[n]->gint_errfcn = master->gint_errfcn;
    wincon[n]->da = master->da;
    wincon[n]->swr_type = master->swr_type;

    /* FR: Strictly speaking this is not good
     *     It will obviously work for shared memory and luckily
     *     the MPI version has its own master.
     *     Also, this is probably not the only place we've got
     *     this issue
     */
    wincon[n]->trinfo_3d = master->trinfo_3d;
    wincon[n]->trinfo_2d = master->trinfo_2d;
    wincon[n]->trinfo_sed = master->trinfo_sed;
#if defined(HAVE_ECOLOGY_MODULE)
    wincon[n]->do_eco = master->do_eco;
    wincon[n]->ecodt = master->ecodt;
#endif
#if defined(HAVE_SEDIMENT_MODULE)
    wincon[n]->do_sed = master->do_sed;
#endif
#if defined(HAVE_WAVE_MODULE)
    wincon[n]->do_wave = master->do_wave;
    wincon[n]->wavedt = master->wavedt;
#endif
    strcpy(wincon[n]->timeunit, master->timeunit);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise tracer info                                          */
  for (n = 1; n <= nwindows; n++) {

    if (nwindows > 1) {
      /* 3D Tracers                                                  */
      wincon[n]->trinfo_3d = (tracer_info_t *)malloc(sizeof(tracer_info_t) *
                                                   master->ntr);
      memset(wincon[n]->trinfo_3d, 0, sizeof(tracer_info_t) * master->ntr);
      for (tn = 0; tn < master->ntr; tn++) {
        tracer_copy(&wincon[n]->trinfo_3d[tn], &master->trinfo_3d[tn]);
        strcpy(wincon[n]->trinfo_3d[tn].name, master->trinfo_3d[tn].name);
      }

      /* 2D Tracers                                                  */
      if (wincon[n]->ntrS) {
        wincon[n]->trinfo_2d =
          (tracer_info_t *)malloc(sizeof(tracer_info_t) * wincon[n]->ntrS);
        memset(wincon[n]->trinfo_2d, 0,
               sizeof(tracer_info_t) * wincon[n]->ntrS);
        for (tn = 0; tn < wincon[n]->ntrS; tn++) {
          tracer_copy(&wincon[n]->trinfo_2d[tn], &master->trinfo_2d[tn]);
          strcpy(wincon[n]->trinfo_2d[tn].name, master->trinfo_2d[tn].name);
        }
      }

      /* Sediment tracers                                            */
      if (wincon[n]->nsed) {
        wincon[n]->trinfo_sed =
          (tracer_info_t *)malloc(sizeof(tracer_info_t) * wincon[n]->nsed);
        memset(wincon[n]->trinfo_sed, 0,
               sizeof(tracer_info_t) * wincon[n]->nsed);
        for (tn = 0; tn < wincon[n]->nsed; tn++) {
          tracer_copy(&wincon[n]->trinfo_sed[tn], &master->trinfo_sed[tn]);
          strcpy(wincon[n]->trinfo_sed[tn].name, master->trinfo_sed[tn].name);
        }
      }
    }

    wincon[n]->trname = (char **)malloc(wincon[n]->ntr * sizeof(char *));
    for (tn = 0; tn < master->ntr; tn++) {
      wincon[n]->trname[tn] = (char *)malloc(sizeof(char)*MAXSTRLEN);
      strcpy(wincon[n]->trname[tn], master->trname[tn]);
    }
    wincon[n]->mintr = d_alloc_1d(wincon[n]->ntr);
    wincon[n]->maxtr = d_alloc_1d(wincon[n]->ntr);
    wincon[n]->advect = i_alloc_1d(wincon[n]->ntr);
    wincon[n]->diffuse = i_alloc_1d(wincon[n]->ntr);
    if (wincon[n]->ntdif_h)
      wincon[n]->tdif_h = i_alloc_1d(wincon[n]->ntdif_h);
    if (wincon[n]->ntdif_v)
      wincon[n]->tdif_v = i_alloc_1d(wincon[n]->ntdif_v);
    wincon[n]->ntbdy = 0;
    wincon[n]->nrlx = master->nrlx;
    wincon[n]->nres = master->nres;
    wincon[n]->nres2d = master->nres2d;
    if (wincon[n]->nrlx)
      wincon[n]->relax = i_alloc_1d(wincon[n]->nrlx);
    if (wincon[n]->nres)
      wincon[n]->reset = i_alloc_1d(wincon[n]->nres);
    if (wincon[n]->nres2d)
      wincon[n]->reset2d = i_alloc_1d(wincon[n]->nres2d);

    /* Get the tracers to advect and diffuse for the window          */
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      wincon[n]->advect[tn] = master->advect[tn];
      wincon[n]->diffuse[tn] = master->diffuse[tn];
      wincon[n]->mintr[tn] = master->mintr[tn];
      wincon[n]->maxtr[tn] = master->maxtr[tn];
      if (master->advect[tn])
        wincon[n]->ntbdy++;
      if (master->diffuse[tn] && !master->advect[tn])
        wincon[n]->ntbdy++;
    }
    for (tn = 0; tn < wincon[n]->nrlx; tn++)
      wincon[n]->relax[tn] = master->relax[tn];
    for (tn = 0; tn < wincon[n]->nres; tn++)
      wincon[n]->reset[tn] = master->reset[tn];
    for (tn = 0; tn < wincon[n]->nres2d; tn++)
      wincon[n]->reset2d[tn] = master->reset2d[tn];
    for (tn = 0; tn < wincon[n]->ntdif_h; tn++)
      wincon[n]->tdif_h[tn] = master->tdif_h[tn];
    for (tn = 0; tn < wincon[n]->ntdif_v; tn++)
      wincon[n]->tdif_v[tn] = master->tdif_v[tn];

    /* Get the tracers which need OBC's to be set                    */
    if (wincon[n]->ntbdy)
      wincon[n]->tbdy = i_alloc_1d(wincon[n]->ntbdy);
    wincon[n]->ntbdy = 0;
    /* Tracers that are advected */
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      if (master->advect[tn]) {
        wincon[n]->tbdy[wincon[n]->ntbdy] = tn;
        wincon[n]->ntbdy++;
      }
      if (master->diffuse[tn] && !master->advect[tn]) {
        wincon[n]->tbdy[wincon[n]->ntbdy] = tn;
        wincon[n]->ntbdy++;
      }
    }

    /* Get the tracers that have a surface flux tracer set           */
    for (tn = 0; tn < master->ntr; tn++) {
      char key[MAXSTRLEN], buf[MAXSTRLEN];
      tracer_info_t *tinfo = &master->trinfo_3d[tn];
      wincon[n]->sflux[tn] = -1;
      if (strlen(tinfo->tag) && tinfo->diffuse) {
	strcpy(buf, tinfo->tag);
	if (decode_tag(buf, "surf_flux", key)) {
	  for (j = 0; j < master->ntrS; j++) {
	    if (strcmp(key, master->trinfo_2d[j].name) == 0) {
	      wincon[n]->sflux[tn] = j;
	      break;
	    }
	  }
	  if (wincon[n]->sflux[tn] < 0)
	    hd_warn("win_consts_init: window%d: Can't find 2D surface flux tracer '%s' for tracer '%s'.\n",
		    n, key, tinfo->name);
	}
      }
    }

    /* Get the tracers which are to be decayed & diagnostic tracers  */
    wincon[n]->ntdec = wincon[n]->ndia = 0;
    wincon[n]->dectr = d_alloc_1d(wincon[n]->ntr);
    memset(wincon[n]->dectr, 0, wincon[n]->ntr * sizeof(double));
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      tracer_info_t *tinfo = &master->trinfo_3d[tn];
      if (!tinfo->diagn && tinfo->inwc && strlen(tinfo->decay))
        wincon[n]->ntdec++;
      if (tinfo->diagn)
        wincon[n]->ndia++;
      if ((i = tracer_find_index(tinfo->decay, master->ntr, master->trinfo_3d)) >= 0) {
	wincon[n]->trinfo_3d[tn].flag |= DE_TR3;
	master->trinfo_3d[tn].flag |= DE_TR3;
	params->trinfo_3d[tn].flag |= DE_TR3;
	wincon[n]->dectr[tn] = (double)i;
      } else if ((i = tracer_find_index(tinfo->decay, master->ntrS, master->trinfo_2d)) >= 0) {
	wincon[n]->trinfo_3d[tn].flag |= DE_TR2;
	master->trinfo_3d[tn].flag |= DE_TR2;
	params->trinfo_3d[tn].flag |= DE_TR2;
	wincon[n]->dectr[tn] = (double)i;
      } else
	wincon[n]->dectr[tn] = atof(tinfo->decay);
    }

    if (wincon[n]->ntdec)
      wincon[n]->tdec = i_alloc_1d(wincon[n]->ntdec);
    if (wincon[n]->ndia)
      wincon[n]->diagn = i_alloc_1d(wincon[n]->ndia);
    wincon[n]->ntdec = wincon[n]->ndia = 0;
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      tracer_info_t *tinfo = &master->trinfo_3d[tn];
      if (tinfo->diagn) {
        wincon[n]->diagn[wincon[n]->ndia] = tn;
        wincon[n]->ndia++;
      }
      if (!tinfo->diagn && tinfo->inwc && strlen(tinfo->decay)) {
        wincon[n]->tdec[wincon[n]->ntdec] = tn;
        wincon[n]->ntdec++;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Initialise window data                                          */
  for (n = 1; n <= nwindows; n++) {

    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      k = geom->s2k[c];

      /* If a ghost cell is encountered, save it's wet neighbour in  */
      /* wincon->s2 and continue.                                    */
      wincon[n]->s2[cc] = 0;
      if (geom->fm[c].wn != n && geom->wgst[c]) {
	wincon[n]->s2[cc] = window[n]->wgst[cc];
        continue;
      }

      /* Cell centered 3D arrays                                     */
      if (master->q != NULL)
        wincon[n]->q[cc] = master->q[c];
      wincon[n]->dz[cc] = master->dz[c];

      wincon[n]->u1kh[cc] = master->u1kh[c];

      if (cc <= window[n]->enonS) {
        /* Cell centered 2D arrays                                   */
        wincon[n]->Cd[cc] = master->Cd[c];
        wincon[n]->z0[cc] = master->z0[c];

	/* Coriolis                                                  */
        wincon[n]->coriolis[cc] = master->coriolis[c];

	/* Fetch                                                     */
	if (master->fetch) {
	  for (tn = 0; tn < 8; tn++)
	    wincon[n]->fetch[cc][tn] = master->fetch[c][tn];
	}
      }
    }

    /* Edge arrays                                                   */
    for (ee = 1; ee <= window[n]->n3_e1; ee++) {
      e = window[n]->wse[ee];
      wincon[n]->u1vh[ee] = master->u1vh[e];
      if (ee <= window[n]->n2_e1) {
        wincon[n]->u1c1[ee] = master->u1c1[e];
        wincon[n]->u1c3[ee] = master->u1c3[e];
        wincon[n]->u1c4[ee] = master->u1c4[e];
        wincon[n]->u1c5[ee] = master->u1c5[e];
        wincon[n]->u1c6[ee] = master->u1c6[e];
      }
    }

    /* Initialise the sigma arrays                                   */
    if (wincon[n]->sigma) {
      for (cc = 1; cc <= window[n]->enonS; cc++) {
        c = window[n]->wsa[cc];
        wincon[n]->Ds[cc] = master->Ds[c];
        wincon[n]->Hs[cc] = master->Hs[c];
        wincon[n]->Hn1[cc] = master->Hn1[c];
        wincon[n]->Hn2[cc] = master->Hn2[c];
      }
    } else {
      for (ee = 1; ee <= window[n]->n2_e1; ee++) {
	e = window[n]->w2_e1[ee];
	wincon[n]->mdx[e] = 1.0;
      }
      wincon[n]->Ds = wincon[n]->one;
      wincon[n]->Hs = wincon[n]->one;
      wincon[n]->Hn1 = wincon[n]->one;
      wincon[n]->Hn2 = wincon[n]->one;
    }

    /* Set linear advection flags                                    */
    wincon[n]->dolin_u1 = wincon[n]->dolin_u2 = 0;
    wincon[n]->dobdry_u1 = wincon[n]->dobdry_u2 = 0;
    wincon[n]->linmask_u1 = wincon[n]->linmask_u1 = NULL;
    wincon[n]->obcmap = i_alloc_1d(window[n]->sgsiz);
    memset(wincon[n]->obcmap, 0, window[n]->sgsiz * sizeof(int));
    for (i = 0; i < window[n]->nobc; i++) {
      open_bdrys_t *open = window[n]->open[i];

      /* Map from cells to open boundaries                           */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	wincon[n]->obcmap[c] = i;
      }
    }

    for (i = 0; i < geom->nobc; i++) {
      open_bdrys_t *open = geom->open[i];
      /* Get the linear momentum zones                               */
      /* Note: including CUSTOM for normal velocities in dobdry_     */
      /* allows the thin layer algorithms to be properly implemented */
      /* on the open boundaries.                                     */
      if (open->type & U1BDRY) {
	if (open->bcond_nor & LINEAR || open->linear_zone_nor)
	  wincon[n]->dolin_u1 = 1;
	if (open->bcond_nor & (LINEAR|NOTHIN|CUSTOM))
	  wincon[n]->dobdry_u1 = 1;
	if (open->stagger & INFACE)
	  wincon[n]->dobdry_u1 = 1;
      }
      if (open->type & U2BDRY && (open->bcond_tan & LINEAR ||
				  open->linear_zone_tan))
	wincon[n]->dolin_u1 = 1;
      if (open->type & U2BDRY) {
	if (open->bcond_nor & LINEAR || open->linear_zone_nor)
	  wincon[n]->dolin_u2 = 1;
	if (open->bcond_nor & (LINEAR|NOTHIN|CUSTOM))
	  wincon[n]->dobdry_u2 = 1;
	if (open->stagger & INFACE)
	  wincon[n]->dobdry_u2 = 1;
      }
      if (open->type & U1BDRY && (open->bcond_tan & LINEAR ||
				  open->linear_zone_tan))
	wincon[n]->dolin_u2 = 1;
    }

    /* Set the mask for u1 linear momentum cells                     */
    if (wincon[n]->dolin_u1) {
      wincon[n]->linmask_u1 = i_alloc_1d(window[n]->szeS);
      memset(wincon[n]->linmask_u1, 0, window[n]->sgsizS * sizeof(int));
      for (i = 0; i < window[n]->nobc; i++) {
	open_bdrys_t *open = window[n]->open[i];
	if (open->bcond_nor & LINEAR) {
	  for(cc = 1; cc <= open->no2_e1; cc++) {
	    c = open->obc_e1[cc];
	    wincon[n]->linmask_u1[c] = 1;
	  }
	}
	if (open->linear_zone_nor) {
	  int *imap;
	  for(cc = 1; cc <= open->no2_e1; cc++) {
	    c = open->oi1_e1[cc];
	    imap = (open->ceni[cc]) ? window[n]->em : window[n]->ep;
	    for (j=0; j < open->linear_zone_nor; j++) {
	      wincon[n]->linmask_u1[c] = 1;
	      c = imap[c];
	    }
	  }
	}
      }
    }

    /* Save the amalgamated OBC conditions for each advected tracer  */
    wincon[n]->obctr = i_alloc_1d(wincon[n]->ntr * sizeof(int));
    memset(wincon[n]->obctr, 0, wincon[n]->ntr * sizeof(int));
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      for (i = 0; i < window[n]->nobc; i++) {
	open_bdrys_t *open = window[n]->open[i];
	wincon[n]->obctr[tn] |= open->bcond_tra[tn];
      }
    }

    /* Get the maps for semi-Lagrange advection                      */
    if(master->trasc & (LAGRANGE|FFSL) || 
       master->momsc & LAGRANGE ||
       master->do_pt) {
      wincon[n]->trasf = TR_FIRST;
      wincon[n]->osl = master->osl;

      if (wincon[n]->osl & L_LINEAR)
	strcpy(wincon[n]->trasr, "linear");
      else if (wincon[n]->osl & L_BILIN)
	strcpy(wincon[n]->trasr, "bilinear");
      else if (wincon[n]->osl & L_BAYLIN)
	strcpy(wincon[n]->trasr, "baylinear");
      else if (wincon[n]->osl & L_SIB)
	strcpy(wincon[n]->trasr, "nn_sibson");
      else if (wincon[n]->osl & L_NONSIB)
	strcpy(wincon[n]->trasr, "nn_non_sibson");
      else if (wincon[n]->osl & L_CUBIC)
	strcpy(wincon[n]->trasr, "cubic");
      else if (wincon[n]->osl & L_LSQUAD)
	strcpy(wincon[n]->trasr, "quadratic");
      else
	strcpy(wincon[n]->trasr, "linear");

    }

    /* Set the offset for means                                      */
    if (master->means_os) {
      wincon[n]->means_os = master->means_os;
      for (cc = 1; cc <= window[n]->enonS; cc++)
	windat->meanc[cc] += wincon[n]->means_os;
    }

    if (DEBUG("init_w"))
      dlog("init_w", "Window %d initialised OK\n", n);
  }

  /*-----------------------------------------------------------------*/
  /* Set values in ghost locations                                   */
  for (n = 1; n <= nwindows; n++) {
    for (cc = 1; cc <= window[n]->enonS; cc++) {
      c = c1 = cc;
      cb = wincon[n]->s2[cc];
      if (cb) {
        wincon[n]->Ds[c] = wincon[n]->Ds[cb];
        /* 3D variables                                              */
        while (c1 != window[n]->zm1[c1]) {
          /*wincon[n]->u1vh[c1] = wincon[n]->u1vh[cb];*/
          wincon[n]->u1kh[c1] = wincon[n]->u1kh[cb];
          c1 = window[n]->zm1[c1];
        }
      }
    }

#if defined(HAVE_SEDIMENT_MODULE)
    /* Sediment temporary values                                     */
    if (wincon[n]->do_sed) {
      for (cc = 1; cc <= window[n]->enonS; cc++) {
        c = window[n]->wsa[cc];
        wincon[n]->d1[cc] = 0.5 * (geom->thetau1[c] + geom->thetau2[c]);
      }
    }
#endif

    /* Set the pointers for horizontal mixing                        */
    hvisc_init(master, wincon);

    window[n]->wincon = wincon[n];
  }

  /* Set the ghost cell map                                          */
  for (n = 1; n <= nwindows; n++) {
    int nf;
    winsize = window[n]->sgsiz;
    memset(wincon[n]->gmap, 0, winsize * sizeof(int));
    for (cc = 1; cc <= window[n]->b3_t; cc++) {
      c = window[n]->w3_t[cc];
      cs = window[n]->m2d[c];
      for (nf = 1; nf <= window[n]->npe[cs]; nf++) {
	c1 = window[n]->c2c[nf][c];
	if (window[n]->c2c[nf][c1] == c1) wincon[n]->gmap[c] |= nf;
      }
    }
    /* Copy to the ghost cells                                       */
    for (cc = 1; cc <= window[n]->nbpt; cc++) {
      c = window[n]->bpt[cc];
      c1 = window[n]->bin[cc];
      cs = window[n]->m2d[c1];
      for (nf = 1; nf <= window[n]->npe[cs]; nf++) {
	if (c == window[n]->c2c[nf][c1])
	  wincon[n]->gmap[c] |= (L_EDGE|U1SOLID);
      }
    }
  }

  /* Get the debug coordinate and window                             */
  if (master->dbc) {
    int wn = geom->fm[master->dbc].wn;
    c = geom->fm[master->dbc].sc;
    for (n = 1; n <= nwindows; n++) {
      wincon[n]->dbw = wn;
      wincon[n]->dbc = c;
      wincon[n]->dbj = master->dbj;
      wincon[n]->dbgf = master->dbgf;
      wincon[n]->dbgtime = master->dbgtime;
    }
  }

  /* Set the csr tidal model constituents if required                */
  csr_tide_grid_init(master, window);

  return (wincon);
}

/* END win_consts_init()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to clear the wincon structure                             */
/*-------------------------------------------------------------------*/
void win_consts_clear(geometry_t **window, int nwindows)
{
  int n;

  for (n = 1; n <= nwindows; n++) {
    win_priv_t *wincon = window[n]->wincon;
#if defined(HAVE_ECOLOGY_MODULE)
    if (wincon->do_eco)
      ecology_destroy(wincon->e);
#endif

    if (wincon->maxtr)
      d_free_1d(wincon->maxtr);
    if (wincon->mintr)
      d_free_1d(wincon->mintr);
    if (wincon->dectr)
      d_free_1d(wincon->dectr);
    if (wincon->advect)
      i_free_1d(wincon->advect);
    if (wincon->diffuse)
      i_free_1d(wincon->diffuse);
    if (wincon->tdec)
      i_free_1d(wincon->tdec);
    if (wincon->tdif_h)
      i_free_1d(wincon->tdif_h);
    if (wincon->tdif_v)
      i_free_1d(wincon->tdif_v);
    if (wincon->diagn)
      i_free_1d(wincon->diagn);
    if (wincon->tbdy)
      i_free_1d(wincon->tbdy);
    if (wincon->twin)
      i_free_1d(wincon->twin);
    if (wincon->reset)
      i_free_1d(wincon->reset);
    if (wincon->relax)
      i_free_1d(wincon->relax);
    if (wincon->sflux)
      i_free_1d(wincon->sflux);
    if (wincon->kth_e1)
      i_free_1d(wincon->kth_e1);
    if (wincon->kth_e2)
      i_free_1d(wincon->kth_e2);
    i_free_1d(wincon->cdry_e1);
    i_free_1d(wincon->cdry_e2);
    i_free_1d(wincon->cbot_e1);
    i_free_1d(wincon->cbot_e2);
    d_free_1d(wincon->oldeta);
    d_free_1d(wincon->dz);
    d_free_1d(wincon->one);
    d_free_1d(wincon->mdx);
    if (wincon->sigma) {
      d_free_1d(wincon->Ds);
      d_free_1d(wincon->Hs);
      d_free_1d(wincon->Hn1);
      d_free_1d(wincon->Hn2);
    }
    d_free_1d(wincon->rdens);
    d_free_1d(wincon->Cd);
    d_free_1d(wincon->u1c1);
    d_free_1d(wincon->u1c3);
    d_free_1d(wincon->u1c4);
    d_free_1d(wincon->u1c5);
    d_free_1d(wincon->u1c6);
    d_free_1d(wincon->topdensu1);
    d_free_1d(wincon->densavu1);
    d_free_1d(wincon->u1adv);
    d_free_1d(wincon->u1inter);
    /* Don't free horizontal mixing variables if Smagorinsky         */
    /* diffusion is used: this is handled when tracers (sdc) are     */
    /* freed.                                                        */
    if (wincon->u1vh && !(wincon->smagcode & U1_SP))
      d_free_1d(wincon->u1vh);
    if (wincon->u1kh && !(wincon->smagcode & U1_SPK))
      d_free_1d(wincon->u1kh);
    d_free_1d(wincon->t11);
    d_free_1d(wincon->t12);
    d_free_1d(wincon->t22);
    if (wincon->q)
      d_free_1d(wincon->q);
    d_free_1d(wincon->z0);
    d_free_1d(wincon->w1);
    d_free_1d(wincon->w2);
    d_free_1d(wincon->w3);
    d_free_1d(wincon->w4);
    if (wincon->w1n)
      d_free_2d(wincon->w1n);
    if (wincon->w2n)
      d_free_2d(wincon->w2n);
    if (wincon->w3n)
      d_free_2d(wincon->w3n);
    if (wincon->w4n)
      d_free_2d(wincon->w4n);
    d_free_1d(wincon->w5);
    d_free_1d(wincon->w6);
    d_free_1d(wincon->w7);
    d_free_1d(wincon->w8);
    d_free_1d(wincon->w9);
    d_free_1d(wincon->w10);
    d_free_1d(wincon->d1);
    d_free_1d(wincon->d2);
    d_free_1d(wincon->d3);
    d_free_1d(wincon->d4);
    d_free_1d(wincon->d5);
    d_free_1d(wincon->d6);
    d_free_1d(wincon->d7);
    d_free_1d(wincon->v1);
    d_free_1d(wincon->v2);
    d_free_1d(wincon->v3);
    d_free_1d(wincon->v4);
    d_free_1d(wincon->v5);
    d_free_1d(wincon->v6);
    d_free_1d(wincon->v7);
    d_free_1d(wincon->v8);
    d_free_1d(wincon->v9);
    d_free_1d(wincon->v10);
    d_free_1d(wincon->v11);
    d_free_1d(wincon->v12);
    i_free_1d(wincon->s1);
    i_free_1d(wincon->s2);
    i_free_1d(wincon->s3);
    i_free_1d(wincon->s4);
    if (wincon->s5)
      i_free_1d(wincon->s5);
    c_free_1d(wincon->c1);
    c_free_1d(wincon->c2);
    i_free_1d(wincon->i1);
    i_free_1d(wincon->i2);
    i_free_1d(wincon->i3);
    i_free_1d(wincon->i4);
    i_free_1d(wincon->i5);
    i_free_1d(wincon->i6);
    i_free_1d(wincon->i7);
    i_free_1d(wincon->gmap);
    d_free_1d(wincon->tmass);
    d_free_1d(wincon->tsmass);
    i_free_1d(wincon->obcmap);
    i_free_1d(wincon->obctr);
    if(wincon->linmask_u1)
      i_free_1d(wincon->linmask_u1);
    if(wincon->linmask_u2)
      i_free_1d(wincon->linmask_u2);
    if(window[n]->sednz)
      d_free_1d(wincon->sd1);
    if(wincon->eta_rlx3d)
      d_free_1d(wincon->eta_rlx3d);
    if(wincon->wgt)
      d_free_2d(wincon->wgt);
    if(wincon->lmap)
      i_free_2d(wincon->lmap);
    if(wincon->m2d)
      i_free_1d(wincon->m2d);
    if (wincon->fetch)
      d_free_2d(wincon->fetch);
    if(wincon->agemsk)
      s_free_1d(wincon->agemsk);
    if(wincon->percmsk)
      s_free_1d(wincon->percmsk);
    d_free_1d(wincon->nu);
    d_free_1d(wincon->nv);
    d_free_1d(wincon->nw);
    i_free_1d(wincon->clxc);
    i_free_1d(wincon->clyc);
    i_free_1d(wincon->clzc);
    i_free_1d(wincon->clxf);
    i_free_1d(wincon->clyf);
    i_free_1d(wincon->clzf);
    d_free_1d(wincon->crfxc);
    d_free_1d(wincon->crfyc);
    d_free_1d(wincon->crfzc);
    d_free_1d(wincon->crfxf);
    d_free_1d(wincon->crfyf);
    d_free_1d(wincon->crfzf);
    d_free_1d(wincon->tr_mod);
    d_free_1d(wincon->tr_mod_x);
    d_free_1d(wincon->tr_mod_y);
    d_free_1d(wincon->tr_mod_z);
    if (wincon->tr_modn)
      d_free_2d(wincon->tr_modn);
    if (wincon->tr_modn_x)
      d_free_2d(wincon->tr_modn_x);
    if (wincon->tr_modn_y)
      d_free_2d(wincon->tr_modn_y);
    if (wincon->tr_modn_z)
      d_free_2d(wincon->tr_modn_z);
    c_free_2d(wincon->gint_error);
    if (nwindows > 1) {
      if (wincon->coriolis)
	d_free_1d(wincon->coriolis);
    }
  }
}

/* END win_consts_clear()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise arrays prior to the run                     */
/*-------------------------------------------------------------------*/
void pre_run_setup(master_t *master,    /* Master data structure     */
                   geometry_t **window, /* Window geometry           */
                   window_t **windat,   /* Window data               */
                   win_priv_t **wincon  /* Data private to windows   */
  )
{
  geometry_t *geom = master->geom;
  int nwindows = master->nwindows;
  int cc, c, ee, e, lc, n, m, v, vv;

  for (n = 1; n <= nwindows; n++) {

    win_data_fill_3d(master, window[n], windat[n], master->nwindows);
    win_data_fill_2d(master, window[n], windat[n], master->nwindows);
    set_map_inside(window[n]);
    memcpy(wincon[n]->oldeta, windat[n]->eta, window[n]->sgsizS * sizeof(double));
    set_dz(window[n], windat[n], wincon[n]);
    get_depths(window[n], windat[n], wincon[n]);
    density_w(window[n], windat[n], wincon[n]);
    if (!(master->vinit & NONE)) {
      win_data_fill_3d(master, window[n], windat[n], master->nwindows);
      set_dz_at_u1(window[n], windat[n], wincon[n]);
      if (master->vinit == VINIT_GEO)
	calc_geostrophic(window[n], windat[n], wincon[n]);
      set_flux_3d(window[n], windat[n], wincon[n], VEL3D);
      set_dz(window[n], windat[n], wincon[n]);

      /* Set a no-gradient over inside stagger OBCs                  */
      for (m = 0; m < window[n]->nobc; m++) {
	open_bdrys_t *open = window[n]->open[m];
	int *mi;
	if (open->stagger & INFACE) {
	  for (ee = 1; ee <= open->no3_e1; ee++) {
	    mi = (!open->ceni[ee]) ? window[n]->em : window[n]->ep;
	    e = open->obc_e1[ee];
	    windat[n]->u1flux3d[mi[e]] = windat[n]->u1flux3d[e];
	  }
	}
      }
      win_data_empty_3d(master, window[n], windat[n], VELOCITY);
    }

    /*---------------------------------------------------------------*/
    /* Coriolis at vertices                                          */
    for (vv = 1; vv <= window[n]->n2_e2; vv++) {
      double d1, d2;
      v = window[n]->w2_e2[vv];
      d1 = d2 = 0.0;
      for (m = 1; m <= window[n]->nvc[v]; m++) {
	c = window[n]->v2c[v][m];
	if (c && !window[n]->wgst[c]) {
	  d2 += wincon[n]->coriolis[c];
	  d1 += 1.0;
	}
      }
      if (d1) windat[n]->fv[v] = d2 / d1;
    }    
  }

  /*-----------------------------------------------------------------*/
  /* Set the reef fractions                                          */
  set_reef_frac(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set sponge zones                                                */
  for (n = 1; n <= nwindows; n++) {
    set_sponge_cells(window[n]);
    set_sponge_c(window[n], wincon[n]->u1kh, windat[n]->dt);
    set_sponge_e(window[n], wincon[n]->u1vh, windat[n]->dt);
  }

  /*-----------------------------------------------------------------*/
  /* Vertival velociyu and velocity initialisation                   */
  if (!(master->vinit & NONE)) {
    for (n = 1; n <= nwindows; n++) {
      /* Transfer velocities for w computation                       */
      win_data_refill_3d(master, window[n], windat[n], master->nwindows, VELOCITY);
      /* Vertical velocity computation                               */
      vel_w_update(window[n], windat[n], wincon[n]);
      vint_3d(window[n], windat[n], wincon[n]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Master fill                                                     */
  master->is_filled = 0;
  master_fill(master, window, windat, wincon);
  master->is_filled = 0;

  /*-----------------------------------------------------------------*/
  /* Variable filling and initialisation                             */
  for (n = 1; n <= nwindows; n++) {
    win_data_fill_3d(master, window[n], windat[n], master->nwindows);
    win_data_fill_2d(master, window[n], windat[n], master->nwindows);
    compute_ref_density(master, window[n], windat[n], wincon[n]);
    Set_lateral_BC_density_w(windat[n]->dens, window[n]->nbpt,
                             window[n]->bpt, window[n]->bin);

    windat[n]->dtf2 = master->dt / master->iratio;
    wincon[n]->calc_closure(window[n], windat[n], wincon[n]);

    /* Set a no-gradient over ghost OBCs                             */
    OBC_bgz_nograd(window[n]);

    /* Set the leapfrog arrays for the first iteration               */
    for (ee = 1; ee <= window[n]->b3_e1; ee++) {
      e = window[n]->w3_e1[ee];
      windat[n]->u1b[e] = windat[n]->u1[e];
      windat[n]->u2b[e] = windat[n]->u2[e];
    }
    for (ee = 1; ee <= window[n]->b2_e1; ee++) {
      e = window[n]->w2_e1[ee];
      windat[n]->u1avb[e] = windat[n]->u1av[e];
    }
    for (cc = 1; cc <= window[n]->b2_t; cc ++) {
      c = window[n]->w2_t[cc];
      windat[n]->etab[c] = wincon[n]->oldeta[c] = windat[n]->eta[c];
    }

    /* Set the lateral boundary conditions for velocity.             */
#if !GLOB_BC
    vel2D_lbc(windat[n]->u1, window[n]->nbpte1, window[n]->nbe1,
	      window[n]->bpte1, window[n]->bine1, wincon[n]->slip);
    vel2D_lbc(windat[n]->u1b, window[n]->nbpte1, window[n]->nbe1,
	      window[n]->bpte1, window[n]->bine1, wincon[n]->slip);
    /*
    vel2D_lbc(windat[n]->u1av, window[n]->nbpte1S, window[n]->nbe1S,
	      window[n]->bpte1S, window[n]->bine1S, wincon[n]->slip);
    */
    set_lateral_bc_eta(windat[n]->etab, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
    set_lateral_bc_eta(wincon[n]->oldeta, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
#endif

    if (master->mode2d)
      mode2d_tracer_init(window[n], windat[n], wincon[n]);

    win_data_empty_3d(master, window[n], windat[n], (VELOCITY|WVEL));
    win_data_empty_2d(master, window[n], windat[n], DEPTH);


    /* Transport mode initialisation                                 */
    if (master->runmode & TRANS) {
      calc_tmass(window[n], windat[n], wincon[n], 
		 wincon[n]->tmass, wincon[n]->tsmass);
      /* Use the viscosity tensors to store streamline Courant nos.  */
      windat[n]->origin = windat[n]->Vz;
      windat[n]->pc = wincon[n]->t11;
      windat[n]->qc = wincon[n]->t12;
      windat[n]->rc = wincon[n]->t22;
      wincon[n]->p1 = windat[n]->u1b;
      /* Initialise the streamline origin for dumps                  */
      if (!(master->tmode & SP_ORIGIN)) {
	for (cc = 1; cc < window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  windat[n]->origin[c] = (double)c;
	}
      }
    } else {
      if (wincon[n]->trasc & LAGRANGE)
	wincon[n]->p1 = d_alloc_1d(window[n]->szc);
    }

    /* Semi-Lagrange initialisation                                  */
    if(wincon[n]->trasc & (LAGRANGE|FFSL) || 
       wincon[n]->momsc & LAGRANGE ||
       wincon[n]->do_pt) {
      tran_grid_init(window[n], windat[n], wincon[n]);
    }

    /* Get the weights for the second derivative                     */
    if(wincon[n]->trasc & (ORDER3US|ORDER4US)) {
      build_advect_weights(window[n], windat[n], wincon[n]); 
    }
    if(wincon[n]->trasc & HIORDER) {
      build_quadratic_weights(window[n], windat[n], wincon[n]); 
    }

    /* Initialise the sediments                                      */
#if defined(HAVE_SEDIMENT_MODULE)
    if (wincon[n]->do_sed)
    {
      emstag(LDEBUG,"hd:windows:pre_run_setup","Setup Sediment...");
      wincon[n]->sediment = sed_init(master->prmfd, window[n]);
      emstag(LDEBUG,"hd:windows:pre_run_setup","Finished Sediment.");
    }
#else
      emstag(LDEBUG,"hd:windows:pre_run_setup","Sediment Ignored.");
#endif

    /* Initialise the waves                                          */
#if defined(HAVE_WAVE_MODULE)
    if (wincon[n]->do_wave) {
      wave_setup(geom, window[n], 1);
      emstag(LDEBUG,"hd:windows:pre_run_setup","Setup Waves...");
      wincon[n]->wave = wave_build(window[n], master->prmfd);
      emstag(LDEBUG,"hd:windows:pre_run_setup","Finished Waves.");
      wave_setup(geom, window[n], 0);
    }
#else
      emstag(LDEBUG,"hd:windows:pre_run_setup","Waves Ignored.");
#endif

    /* Initialise the ecology                                        */
#if defined(HAVE_ECOLOGY_MODULE)
    if (wincon[n]->do_eco) {
      emstag(LDEBUG,"hd:windows:pre_run_setup","Setup Ecology...");
      wincon[n]->e = ecology_build(window[n], master->prmname);
      emstag(LDEBUG,"hd:windows:pre_run_setup","Finished Ecology.");
    }
#else
    emstag(LTRACE,"hd:windows:pre_run_setup","Ecology Ignored.");
#endif

#if defined(HAVE_TRACERSTATS_MODULE)
    emstag(LDEBUG,"hd:windows:pre_run_setup","Setup Tracerstats...");
    prm_set_errfn(hd_silent_warn);
    wincon[n]->trs = trs_build(window[n], master->prmfd);
    emstag(LDEBUG,"hd:windows:pre_run_setup","Finished Tracerstats.");
#else
    emstag(LTRACE,"hd:windows:pre_run_setup","Tracerstats Ignored.");
#endif
  }

  /* Transfer the backward arrays to auxiliary cells                 */
  /* Slave to master */
  for (n = 1; n <= nwindows; n++) {
    s2m_vel(master->u1avb, windat[n]->u1avb,
            window[n]->s2me1, window[n]->wse, window[n]->ns2me1S);
  }
  /* Master to slave                                                 */
  for (n = 1; n <= nwindows; n++) {
    int ce1;
    for (ee = 1; ee <= window[n]->nm2se1S; ee++) {
      e = window[n]->m2se1[ee];
      ce1 = window[n]->wse[e];
      windat[n]->u1avb[e] = master->u1avb[ce1];
    }
  }

  get_timesteps(window, windat, wincon, nwindows, master);
  timeaux(window, windat, wincon, nwindows);
  init_flushing(master, window, wincon);
  init_age(master, window, windat, wincon);
  init_trperc(master, window, windat, wincon);
  init_totals(master, window, wincon);
  init_trans(master, window, wincon);
  nan_check(window, windat, wincon, nwindows);
  pt_setup(master, window, windat, wincon);
}

/* END pre_run_setup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the reef fractions for the porus plate algorithm             */
/*-------------------------------------------------------------------*/
void set_reef_frac(master_t *master, 
		   geometry_t **window, 
		   window_t **windat, 
		   win_priv_t **wincon
		   )
{
  geometry_t *geom = master->geom;
  int n, i, j, c, cs, cc, lc, wn;
  FILE *pf;
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  double d1;

  if (master->porusplate) {
    not_included("porus plate");
    master->porusplate = 0;
    return;
  }

  if (master->porusplate && strlen(master->reef_frac)) {
    n = parseline(params->reef_frac, fields, MAXNUMARGS);
    if (n == 2) {
      /* 'Along i' file sets grid ratio from i to i+1, i.e. through  */
      /* e2 faces.                                                   */   
      if ((pf = fopen(fields[0], "r")) == NULL)
	hd_quit("set_reef_frac: Can't open e1 reef fraction file %s\n", fields[0]);
      fgets(buf, 256, pf);
      while (fgets(buf, 256, pf) != NULL) {
	sscanf(buf,"%d %d %lf",&i,&j,&d1);
	if(i >= 0 && i < geom->nce1 && j >=0 && j < geom->nce2) {
	  cs = geom->map[geom->nz-1][j][i];
	  if (cs > 0 && cs < geom->enonS) {
	    c = geom->bot_t[geom->c2cc[cs]];
	    if (c > 0) {
	      wn = geom->fm[c].wn;
	      lc = geom->fm[c].sc;
	      master->reefe1[c] = windat[wn]->reefe1[lc] = d1;
	    }
	  }
	}
      }
      fclose(pf);
      /* 'Along j' file sets grid ratio from j to j+1, i.e. through  */
      /* e1 faces.                                                   */
      if ((pf = fopen(fields[1], "r")) == NULL)
	hd_quit("set_reef_frac: Can't open e2 reef fraction file %s\n", fields[1]);
      fgets(buf, 256, pf);
      while (fgets(buf, 256, pf) != NULL) {
	sscanf(buf,"%d %d %lf",&i,&j,&d1);
	if(i >=0 && i < geom->nce1 && j >=0 && j < geom->nce2) {
	  cs = geom->map[geom->nz-1][j][i];
	  if (cs > 0 && cs < geom->enonS) {
	    c = geom->bot_t[geom->c2cc[cs]];
	    if (c > 0) {
	      wn = geom->fm[c].wn;
	      lc = geom->fm[c].sc;
	      master->reefe2[c] = windat[wn]->reefe2[lc] = d1;
	    }
	  }
	}
      }
      fclose(pf);
    } else if (n == 1) {
      d1 = atof(master->reef_frac);
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->bot_t[cc];
	cs = geom->m2d[c];
	wn = geom->fm[c].wn;
	lc = geom->fm[c].sc;
	master->reefe1[c] = windat[wn]->reefe1[lc] = d1;
	master->reefe2[c] = windat[wn]->reefe2[lc] = d1;
      }
    } else
      hd_warn("porusplate: require 2 reef fraction files.\n");
  }
}

/* END set_reef_frac()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the window structure to initialize the wave modeule       */
/*-------------------------------------------------------------------*/
void wave_setup(geometry_t *geom, geometry_t *window, int mode)
{
  win_priv_t *wincon = window->wincon;
  int c, cc;

  if (window->thetau1 == NULL) {
    window->thetau1 = d_alloc_1d(window->szcS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->thetau1[cc] = geom->thetau1[c];
      }
    }
  }
  if (window->thetau2 == NULL) {
    window->thetau2 = d_alloc_1d(window->szcS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->thetau2[cc] = geom->thetau2[c];
      }
    }
  }
  if (window->sinthcell == NULL) {
    window->sinthcell = d_alloc_1d(window->szcS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->sinthcell[cc] = geom->sinthcell[c];
      }
    }
  }
  if (window->costhcell == NULL) {
    window->costhcell = d_alloc_1d(window->szcS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->costhcell[cc] = geom->costhcell[c];
      }
    }
  }
  if (!mode) {
    if(window->thetau1) d_free_1d(window->thetau1);
    if(window->thetau2) d_free_1d(window->thetau2);
    if(window->sinthcell) d_free_1d(window->sinthcell);
    if(window->costhcell) d_free_1d(window->costhcell);
    if(wincon->fetch) d_free_2d(wincon->fetch);
    window->thetau1 = window->thetau2 = NULL;
    window->sinthcell = window->costhcell = NULL;
    wincon->fetch = NULL;
  }
}

/* END pre_wave_setup()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Free global arrays that are no longer required                    */
/*-------------------------------------------------------------------*/
void geom_free_us(master_t *master,  /* Master data structure        */
		  geometry_t *geom,  /* Global geometry              */
		  int mode           /* Code for arrays to free      */
  )
{
  int n, tn, flg;

  if (mode & ALL) {
    /* Global maps                                                   */
#if !TR_CK
    free((global_map_t *)geom->fm);
#endif

    if (master->nwindows > 1) {
      master_free_nwin(master);
      i_free_2d(geom->c2c);
      i_free_2d(geom->c2e);
      i_free_2d(geom->c2v);
      i_free_2d(geom->e2e);
      i_free_2d(geom->e2c);
      i_free_2d(geom->e2v);
      i_free_2d(geom->v2c);
      i_free_2d(geom->v2e);

      i_free_2d(geom->eSc);
      i_free_2d(geom->eSv);
      i_free_2d(geom->eSe);
      i_free_2d(geom->vIc);
      d_free_2d(geom->wAe);
      d_free_2d(geom->wSe);

      i_free_1d(geom->zp1);
      i_free_1d(geom->zm1);
      i_free_1d(geom->zp1e);
      i_free_1d(geom->zm1e);
    }

    /* Cells to process vectors                                      */
    i_free_1d(geom->w3_t);
    i_free_1d(geom->w3_e1);
    i_free_1d(geom->w3_e2);

    /* Sparse arrays                                                 */
    d_free_1d(geom->gridz);
    d_free_1d(geom->cellz);

    if (geom->sednz) {
      d_free_2d(geom->gridz_sed);
      d_free_2d(geom->cellz_sed);
    }
    if (geom->gsed_t)
      i_free_1d(geom->gsed_t);
    if (geom->ised_t)
      i_free_1d(geom->ised_t);
    d_free_1d(geom->h1acell);
    d_free_1d(geom->sinthcell);
    d_free_1d(geom->costhcell);
    d_free_1d(geom->thetau1);
    d_free_1d(geom->sinthu1);
    d_free_1d(geom->costhu1);
    d_free_1d(geom->thetau2);
    d_free_1d(geom->sinthu2);
    d_free_1d(geom->costhu2);
    d_free_1d(geom->cellarea);
    d_free_1d(geom->botz);
    d_free_1d(geom->botzgrid);
    d_free_1d(geom->cellx);
    d_free_1d(geom->celly);
    d_free_1d(geom->gridx);
    d_free_1d(geom->gridy);
    d_free_1d(geom->u1x);
    d_free_1d(geom->u1y);
    d_free_1d(geom->dHde1);
    d_free_1d(geom->h1au1);
    d_free_1d(geom->h2au1);
    d_free_1d(geom->botzu1);
    if (geom->sm_e1)
      i_free_1d(geom->sm_e1);
    if (geom->sm_e2)
      i_free_1d(geom->sm_e2);
  }
  if (mode & (UNUSED | ALL)) {

    if (geom->mgc) i_free_1d(geom->mgc);

#if defined(HAVE_SEDIMENT_MODULE)
    if (geom->sed_t) i_free_1d(geom->sed_t);
#endif

    return;

    if (geom->nobc)
      s_free_2d(geom->owc);
    for (n = 0; n < geom->nobc; n++) {
      /* Open boundary vectors                                       */
      i_free_1d(geom->open[n]->oi1_t);
      i_free_1d(geom->open[n]->oi1_e1);
      i_free_1d(geom->open[n]->oi2_t);
      i_free_1d(geom->open[n]->oi2_e1);
      /* Some of the OBC vectors may be required for reading data    */
      /* from file on the master; don't deallocate in this case.     */
      flg = 1;
      for (tn = 0; tn < geom->open[n]->ntr; tn++)
        if (geom->open[n]->bcond_tra[tn] & (FILEIN))
          flg = 0;
      if (flg && !(geom->open[n]->bcond_ele & (FILEIN)))
        i_free_1d(geom->open[n]->obc_t);
      if (geom->open[n]->type & U1BDRY &&
          !(geom->open[n]->bcond_nor & (FILEIN | CUSTOM)))
        i_free_1d(geom->open[n]->obc_e1);
      if (geom->open[n]->type & U2BDRY &&
          !(geom->open[n]->bcond_nor & (FILEIN | CUSTOM)))
        i_free_1d(geom->open[n]->obc_e2);
    }
  }
}

/* END geom_free()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to initialise the window data structure                   */
/*-------------------------------------------------------------------*/
void window_init(geometry_t *geom,    /* Global geometry             */
                 geometry_t **window  /* Window geometry             */
  )
{
  int nwindows;                 /* Number of windows                 */
  int n, m, cc, ee, vv;         /* Counters                          */
  int c, e, v, c1, cb, eb;      /* Sparse locations                  */
  int i, j, k;                  /* Cartesian locations               */
  int **s2, **se, **sv;         /* Dummy buffer for ghost locations  */

  nwindows = geom->nwindows;
  c = e = v = 0;
  for (n = 1; n <= nwindows; n++) {
    if (window[n]->enon > c)
      c = window[n]->enon;
    if (window[n]->sze > e)
      e = window[n]->sze;
    if (window[n]->szv > v)
      v = window[n]->szv;
  }
  s2 = i_alloc_2d(c + 1, nwindows + 1);
  se = i_alloc_2d(e, nwindows + 1);
  sv = i_alloc_2d(v, nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Initialise window data                                          */
  for (n = 1; n <= nwindows; n++) {

    window[n]->topgrid = geom->topgrid;
    window[n]->totarea = geom->totarea;

    /* Set geog flag */
    window[n]->is_geog = geom->is_geog;
    
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      k = geom->s2k[c];

      /* Cell centered 3D arrays                                     */
      window[n]->gridz[cc] = geom->gridz[c];
      window[n]->cellz[cc] = geom->cellz[c];

      /* If a ghost cell is encountered, save it's wet neighbour in  */
      /* wincon->s2 and continue.                                    */
      s2[n][cc] = 0;
      if (geom->fm[c].wn != n && geom->wgst[c]) {
	s2[n][cc] = window[n]->wsa[window[n]->wgst[cc]];
        continue;
      }

      if (cc <= window[n]->enonS) {

        /* Cell centered 2D arrays                                   */
        window[n]->cellarea[cc] = geom->cellarea[c];
        window[n]->botz[cc] = geom->botz[c];
        window[n]->cellx[cc] = geom->cellx[c];
        window[n]->celly[cc] = geom->celly[c];

	for (m = 1; m <= window[n]->npe[cc]; m++)
	  window[n]->hacell[m][cc] = geom->hacell[m][c];

        /* Sediments                                                 */
        if (geom->sednz > 0) {
          int kk;
          for (kk = 0; kk < geom->sednz + 1; kk++)
            window[n]->gridz_sed[kk][cc] = geom->gridz_sed[kk][c];
          for (kk = 0; kk < geom->sednz; kk++)
            window[n]->cellz_sed[kk][cc] = geom->cellz_sed[kk][c];
        }
      }
    }

    /* Edge centered arrays 2D arrays                                */
    for (ee = 1; ee < window[n]->szeS; ee++) {
      e = window[n]->wse[ee];

      se[n][ee] = 0;
      c = geom->e2c[e][0];
      cc = window[n]->e2c[ee][0];
      j = geom->e2e[e][0];
      if (geom->fm[c].wn != n && geom->wgst[c]) {
	se[n][ee] = window[n]->wse[window[n]->c2e[j][window[n]->wgst[cc]]];
      }
      window[n]->h1acell[ee] = geom->h1acell[e];
      window[n]->h1au1[ee] = geom->h1au1[e];
      window[n]->h2au1[ee] = geom->h2au1[e];
      window[n]->u1x[ee] = geom->u1x[e];
      window[n]->u1y[ee] = geom->u1y[e];
      window[n]->botzu1[ee] = geom->botzu1[e];
      window[n]->dHde1[ee] = geom->dHde1[e];
      window[n]->costhu1[ee] = geom->costhu1[e];
      window[n]->sinthu1[ee] = geom->sinthu1[e];
      window[n]->costhu2[ee] = geom->costhu2[e];
      window[n]->sinthu2[ee] = geom->sinthu2[e];
    }

    /* Vertex centered arrays 2D arrays                              */
    for (vv = 1; vv < window[n]->szvS; vv++) {
      v = window[n]->wsv[vv];
 
      sv[n][vv] = 0;
      c = geom->v2c[v][1];
      cc = window[n]->v2c[vv][1];
      if (geom->fm[c].wn != n && geom->wgst[c]) {
	sv[n][vv] = window[n]->wsv[window[n]->c2v[1][window[n]->wgst[vv]]];
        continue;
      }
 
      window[n]->gridx[vv] = geom->gridx[v];
      window[n]->gridy[vv] = geom->gridy[v];
      window[n]->botzgrid[vv] = geom->botzgrid[v];
    }
    if (DEBUG("init_w"))
      dlog("init_w", "Window geometry %d initialised OK\n", n);
  }

  /*-----------------------------------------------------------------*/
  /* Set values in ghost locations                                   */
  if (nwindows > 1) {

  for (n = 1; n <= nwindows; n++) {

    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = cc;
      cb = s2[n][cc];
      if (cb) {
        window[n]->gridz[c] = geom->gridz[window[n]->wsa[cc]];
      }
    }
    for (cc = 1; cc <= window[n]->enonS; cc++) {
      c = cc;
      cb = s2[n][cc];

      if (cb) {
        window[n]->cellarea[c] = geom->cellarea[cb];
        window[n]->cellx[c] = geom->cellx[window[n]->wsa[cc]];
        window[n]->celly[c] = geom->celly[window[n]->wsa[cc]];
      }
    }
    for (ee = 1; ee < window[n]->szeS; ee++) {
      cb = se[n][ee];
      if (cb) {
        window[n]->h1au1[ee] = geom->h1au1[cb];
        window[n]->h1acell[ee] = geom->h1acell[cb];
        window[n]->h2au1[ee] = geom->h2au1[cb];
	window[n]->costhu1[ee] = geom->costhu1[cb];
	window[n]->sinthu1[ee] = geom->sinthu1[cb];
	window[n]->costhu2[ee] = geom->costhu2[cb];
	window[n]->sinthu2[ee] = geom->sinthu2[cb];
      }
    }

    for (vv = 1; vv < window[n]->szvS; vv++) {
      cb = sv[n][vv];
      if (cb) {
        window[n]->gridx[vv] = geom->gridx[window[n]->wsv[vv]];
        window[n]->gridy[vv] = geom->gridy[window[n]->wsv[vv]];
        window[n]->botzgrid[vv] = geom->botzgrid[window[n]->wsv[vv]];
      }
    }

    for (i = 0; i < window[n]->nobc; i++) {
      open_bdrys_t *open = window[n]->open[i];
      /*
      for(cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	cb = open->oi1_t[cc];
	window[n]->cellarea[c] = window[n]->cellarea[cb];
      }
      for(ee = 1; ee <= open->no2_e1; ee++) {
	e = open->obc_e1[ee];
	eb = open->oi1_e1[ee];
	window[n]->h1au1[e] = window[n]->h1au1[eb];
	window[n]->h1acell[e] = window[n]->h1acell[eb];
	window[n]->h2au1[e] = window[n]->h2au1[eb];
	window[n]->costhu1[e] = geom->costhu1[window[n]->wse[eb]];
	window[n]->sinthu1[e] = geom->sinthu1[window[n]->wse[eb]];
	window[n]->costhu2[e] = geom->costhu2[window[n]->wse[eb]];
	window[n]->sinthu2[e] = geom->sinthu2[window[n]->wse[eb]];
      }
      */
      /* Metric terms for OBC ghost cells                            */
      for (ee = 1; ee <= open->no2_e1; ee++) {
	c = open->obc_e2[ee];
	cb = window[n]->wsa[c];
	for (j = 0; j < open->bgz; j++) {
	  c = open->omape[ee][c];
	  cb = window[n]->wsa[c];
	  /* Note: cell metrics (e.g. h1au2) could be set here, e.g. */
	  /* window[n]->h1au1[c] = geom->h1au1[cb];                  */
	  /* but these are captured using s2[n] and wgst[] above.    */
	  /*window[n]->botz[c] = geom->botz[cb];*/
	}
	/* Reset the grid size for long rivers                       */
	/*
	if (open->options & OP_RLEN && open->rlen != 0.0) {
	  int p1;
	  c = open->obc_t[cc];
	  p1 = open->nmape[ee][c];
	  c = open->omape[ee][c];

	  if (open->ocodex & R_EDGE)
	  window[n]->h1au1[c] += 0.5 * (window[n]->h1acell[c] + open->rlen - window[n]->h1au1[c]);
	  else
	  window[n]->h1au1[p1] += 0.5 * (window[n]->h1acell[c] + open->rlen - window[n]->h1au1[c]);
	  window[n]->cellarea[c] = open->rlen * window[n]->h2acell[c];
	  window[n]->h1acell[c] += open->rlen;
	}
	*/
      }
    }
  }
  }
  i_free_2d(s2);
}

/* END window_init()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to free the window structure                              */
/*-------------------------------------------------------------------*/
void win_geom_clear(geometry_t **window,  /* Window geometry         */
                    int nwindows          /* Number of windows       */
  )
{
  int n, nn, tn;

  for (n = 1; n <= nwindows; n++) {
    i_free_1d(window[n]->s2i);
    i_free_1d(window[n]->s2j);
    i_free_1d(window[n]->s2k);
    i_free_1d(window[n]->wsa);
    i_free_1d(window[n]->w3_t);
    i_free_1d(window[n]->w3_e1);
    i_free_1d(window[n]->w3_e2);
    i_free_1d(window[n]->w2_t);
    i_free_1d(window[n]->w2_e1);
    i_free_1d(window[n]->w2_e2);
    if (window[n]->c2cc)
      i_free_1d(window[n]->c2cc);
    if (window[n]->mgm)
      i_free_1d(window[n]->mgm);
    if (nwindows > 1) {
      i_free_1d(window[n]->bpt);
      i_free_1d(window[n]->bin);
      i_free_1d(window[n]->bin2);
      i_free_1d(window[n]->bpte1);
      i_free_1d(window[n]->bine1);
      i_free_1d(window[n]->bpte2);
      i_free_1d(window[n]->bine2);
      i_free_1d(window[n]->bpte1S);
      i_free_1d(window[n]->bine1S);
      i_free_1d(window[n]->bpte2S);
      i_free_1d(window[n]->bine2S);
    }
    i_free_1d(window[n]->m2s);
    i_free_1d(window[n]->m2se1);
    i_free_1d(window[n]->s2m);
    i_free_1d(window[n]->s2me1);
    i_free_1d(window[n]->sur_t);
    i_free_1d(window[n]->nsur_t);
    i_free_1d(window[n]->bot_t);
    i_free_1d(window[n]->sur_e1);
    i_free_1d(window[n]->bot_e1);
    i_free_1d(window[n]->sur_e2);
    i_free_1d(window[n]->bot_e2);
    i_free_1d(window[n]->m2d);
    i_free_2d(window[n]->c2c);
    i_free_2d(window[n]->c2e);
    i_free_2d(window[n]->c2v);
    i_free_2d(window[n]->e2e);
    i_free_2d(window[n]->e2c);
    i_free_2d(window[n]->e2v);
    i_free_2d(window[n]->v2c);
    i_free_2d(window[n]->v2e);
    i_free_2d(window[n]->eSc);
    i_free_2d(window[n]->eSv);
    i_free_2d(window[n]->eSe);
    i_free_2d(window[n]->vIc);
    d_free_2d(window[n]->wAe);
    d_free_2d(window[n]->wSe);
    i_free_1d(window[n]->zp1);
    i_free_1d(window[n]->zm1);
    d_free_1d(window[n]->h1acell);
    d_free_1d(window[n]->h2acell);
    d_free_1d(window[n]->h1au1);
    d_free_1d(window[n]->h2au2);
    d_free_1d(window[n]->h2au1);
    d_free_1d(window[n]->h1au2);
    d_free_1d(window[n]->cellarea);
    d_free_1d(window[n]->dHde1);
    d_free_1d(window[n]->dHde2);
    d_free_1d(window[n]->cellx);
    d_free_1d(window[n]->celly);
    d_free_1d(window[n]->gridx);
    d_free_1d(window[n]->gridy);
    d_free_1d(window[n]->u1x);
    d_free_1d(window[n]->u1y);
    d_free_1d(window[n]->u2x);
    d_free_1d(window[n]->u2y);
    d_free_1d(window[n]->botz);
    d_free_1d(window[n]->botzu1);
    d_free_1d(window[n]->botzu2);
    d_free_1d(window[n]->botzgrid);
    d_free_1d(window[n]->gridz);
    d_free_1d(window[n]->cellz);
    if (window[n]->thetau1)
      d_free_1d(window[n]->thetau1);
    if (window[n]->thetau2)
      d_free_1d(window[n]->thetau2);
    if (window[n]->sinthcell)
      d_free_1d(window[n]->sinthcell);
    if (window[n]->costhcell)
      d_free_1d(window[n]->costhcell);
    if (window[n]->sinthu1)
      d_free_1d(window[n]->sinthu1);
    if (window[n]->costhu1)
      d_free_1d(window[n]->costhu1);
    if (window[n]->sinthu2)
      d_free_1d(window[n]->sinthu2);
    if (window[n]->costhu2)
      d_free_1d(window[n]->costhu2);
    if (window[n]->aux_t)
      i_free_1d(window[n]->aux_t);
    if (window[n]->taux_t)
      d_free_1d(window[n]->taux_t);
    if (window[n]->sm_e1)
      i_free_1d(window[n]->sm_e1);
    if (window[n]->sm_e2)
      i_free_1d(window[n]->sm_e2);
    if (window[n]->cbgc)
      i_free_1d(window[n]->cbgc);
    if (window[n]->csed)
      i_free_1d(window[n]->csed);
    if (window[n]->ctran)
      i_free_1d(window[n]->ctran);
    if (window[n]->cwave)
      i_free_1d(window[n]->cwave);
    if (window[n]->ctrst)
      i_free_1d(window[n]->ctrst);
    if (window[n]->cask)
      i_free_1d(window[n]->cask);
    if (window[n]->eask)
      i_free_1d(window[n]->eask);

    /* Open boundary cell vectors */
    for (nn = 0; nn < window[n]->nobc; nn++) {
      i_free_1d(window[n]->open[nn]->obc_t);
      i_free_1d(window[n]->open[nn]->bot_t);
      i_free_1d(window[n]->open[nn]->obc_e1);
      i_free_1d(window[n]->open[nn]->obc_e2);
      i_free_1d(window[n]->open[nn]->oi1_t);
      i_free_1d(window[n]->open[nn]->oi1_e1);
      i_free_1d(window[n]->open[nn]->oi2_t);
      i_free_1d(window[n]->open[nn]->oi2_e1);
      i_free_1d(window[n]->open[nn]->e2c_e1);
      if (window[n]->open[nn]->cyc_t)
        i_free_1d(window[n]->open[nn]->cyc_t);
      if (window[n]->open[nn]->cyc_e1)
        i_free_1d(window[n]->open[nn]->cyc_e1);
      if (window[n]->open[nn]->tmap)
        i_free_1d(window[n]->open[nn]->tmap);
      if (window[n]->open[nn]->tmap_u1)
        i_free_1d(window[n]->open[nn]->tmap_u1);
      if (window[n]->open[nn]->tmap_u2)
        i_free_1d(window[n]->open[nn]->tmap_u2);
      if (window[n]->open[nn]->transfer_eta)
        d_free_1d(window[n]->open[nn]->transfer_eta);
      if (window[n]->open[nn]->transfer_u1)
	d_free_1d(window[n]->open[nn]->transfer_u1);
      if (window[n]->open[nn]->transfer_u2)
        d_free_1d(window[n]->open[nn]->transfer_u2);
      if (window[n]->open[nn]->transfer_u1av)
        d_free_1d(window[n]->open[nn]->transfer_u1av);
      if (window[n]->open[nn]->transfer_u2av)
        d_free_1d(window[n]->open[nn]->transfer_u2av);
      if (window[n]->open[nn]->t_tmap)
        i_free_1d(window[n]->open[nn]->t_tmap);
      if (window[n]->open[nn]->t_imap)
        i_free_2d(window[n]->open[nn]->t_imap);
      if (window[n]->open[nn]->trm)
        i_free_1d(window[n]->open[nn]->trm);
      if (window[n]->open[nn]->t_transfer)
        d_free_2d(window[n]->open[nn]->t_transfer);
      if (window[n]->open[nn]->bcond_tra)
        i_free_1d(window[n]->open[nn]->bcond_tra);
      if (window[n]->open[nn]->clampv)
        d_free_1d(window[n]->open[nn]->clampv);
      if (window[n]->open[nn]->trpc)
        d_free_1d(window[n]->open[nn]->trpc);
      if (window[n]->open[nn]->relax_zone_tra)
        i_free_1d(window[n]->open[nn]->relax_zone_tra);
      if (window[n]->open[nn]->iloc)
        i_free_1d(window[n]->open[nn]->iloc);
      if (window[n]->open[nn]->jloc)
        i_free_1d(window[n]->open[nn]->jloc);
      if (window[n]->open[nn]->ilocc)
        i_free_1d(window[n]->open[nn]->ilocc);
      if (window[n]->open[nn]->jlocc)
        i_free_1d(window[n]->open[nn]->jlocc);
      if (window[n]->open[nn]->etabold)
        free((bdry_old_values_t *)window[n]->open[nn]->etabold);
      if (window[n]->open[nn]->u1_on)
        free((bdry_old_values_t *)window[n]->open[nn]->u1_on);
      if (window[n]->open[nn]->u2_on)
        free((bdry_old_values_t *)window[n]->open[nn]->u2_on);
      if (window[n]->open[nn]->u1_ot)
        free((bdry_old_values_t *)window[n]->open[nn]->u1_ot);
      if (window[n]->open[nn]->u2_ot)
        free((bdry_old_values_t *)window[n]->open[nn]->u2_ot);
      if (window[n]->open[nn]->u1av_on)
        free((bdry_old_values_t *)window[n]->open[nn]->u1av_on);
      if (window[n]->open[nn]->u2av_on)
        free((bdry_old_values_t *)window[n]->open[nn]->u2av_on);
      if (window[n]->open[nn]->u1av_ot)
        free((bdry_old_values_t *)window[n]->open[nn]->u1av_ot);
      if (window[n]->open[nn]->u2av_ot)
        free((bdry_old_values_t *)window[n]->open[nn]->u2av_ot);
      if (window[n]->open[nn]->tideforce)
        free((tide_details_t *)window[n]->open[nn]->tideforce);
      if (window[n]->open[nn]->tide)
        free((tidal_memory_t *)window[n]->open[nn]->tide);
      for (tn = 0; tn < window[n]->open[nn]->ntr; tn++) {
        if (window[n]->open[nn]->cusname_t[tn])
          free((char *)window[n]->open[nn]->cusname_t[tn]);
      }
      if (window[n]->open[nn])
	free((open_bdrys_t **)window[n]->open[nn]);
    }

    if (window[n]->cellz_sed)
      d_free_2d(window[n]->cellz_sed);
    if (window[n]->gridz_sed)
      d_free_2d(window[n]->gridz_sed);
#if defined(HAVE_SEDIMENT_MODULE)
    if(window[n]->sed_t)
      i_free_1d(window[n]->sed_t);
#endif

  }
}

/* END win_geom_clear()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to clear the windows                                      */
/*-------------------------------------------------------------------*/
void windows_clear(hd_data_t *hd_data)
{
  geometry_t *geom = hd_data->geom;
  master_t *master = hd_data->master;
  geometry_t **window = hd_data->window;
  window_t **windat = hd_data->windat;
  win_priv_t **wincon = hd_data->wincon;
  int n;

  /*-----------------------------------------------------------------*/
  /* Copy window data to the master                                  */
  master_fill(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Close existing windows                                          */
  dp_cleanup();

  for (n = 1; n <= geom->nwindows; n++)
    win_data_clear(windat[n]);
  win_consts_clear(window, geom->nwindows);
  win_geom_clear(window, geom->nwindows);
  free((window_t **)windat);
  free((win_priv_t **)wincon);
  free((geometry_t **)window);
}

/* END windows_clear()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the window data structure          */
/*-------------------------------------------------------------------*/
geometry_t *window_alloc(void)
{
  geometry_t *window = (geometry_t *)malloc(sizeof(geometry_t));
  memset(window, 0, sizeof(geometry_t));
  return window;
}

/* END window_alloc()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the window constants structure     */
/*-------------------------------------------------------------------*/
win_priv_t *win_consts_alloc(void)
{
  win_priv_t *wincon = (win_priv_t *)malloc(sizeof(win_priv_t));
  memset(wincon, 0, sizeof(win_priv_t));
  return wincon;
}

/* END win_consts_alloc()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the window data structure          */
/*-------------------------------------------------------------------*/
window_t *win_data_alloc(void)
{
  window_t *windat = (window_t *)malloc(sizeof(window_t));
  memset(windat, 0, sizeof(window_t));
  return windat;
}

/* END win_data_alloc()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to generate auxiliary cell maps in directions that have   */
/* not already been defined. This is accomplished through the global */
/* maps. Any cell that maps outside its window remains self-mapping. */
/*-------------------------------------------------------------------*/
void aux_maps(geometry_t *geom,       /* Global geometry             */
	      geometry_t *window,     /* Local geometry              */
	      int cc                  /* Local sparse coordinate     */
	      )
{
  int j;
  int c, cs;                    /* Global sparse coordinate          */
  int cg;                       /* Mapped global sparse counters     */
  int cl;                       /* Mapped local sparse coordinate    */

  c = window->wsa[cc];
  cs = geom->m2d[c];

  for (j = 1; j <= geom->npe[cs]; j++) {
    if (window->c2c[j][cc] == cc) {
      cg = geom->c2c[j][c];
      cl = geom->fm[cg].ac;
      if (cl)
	window->c2c[j][cc] = cl;
    }
  }

  /* Generate ghost cell maps in directions that have not already    */
  /* been defined.                                                   */
  if (!geom->fm[c].wn) {            /* Window = 0 => ghost cell      */
    for (j = 1; j <= geom->npe[cs]; j++) {
      if (window->c2c[j][cc] == cc) { /* Self pointing ghost cell    */
      cg = geom->c2c[j][c];     /* Global cell in direction j
      cl = geom->fm[cg].ac;     /* Auxiliary cell to the next to cc  */
      if (!cl)
        cl = geom->fm[cg].sc;   /* Wet cell if aux. cell undefined   */
      if (cl)
        window->c2c[j][cc] = cl;
      }
    }
  }
}

/* END aux_maps()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set multi dt auxiliary cells that are associated with  */
/* windows having a longer time-step than the current window to the  */
/* updated values at the end of the long time-step.                  */
/*-------------------------------------------------------------------*/
void fill_multidt_aux(master_t *master,     /* Master data           */
                      geometry_t *window,   /* Window geometry       */
                      window_t *windat      /* Window data           */
  )
{
  int c, cc;                    /* Local sparse coordinate / counter */
  int tn;                       /* Tracer counter                    */

  /* Auxiliary multi-dt cells                                        */
  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    c = window->wsa[c];
    for (tn = 0; tn < windat->ntr; tn++)
      windat->tr_wc_ae[tn][cc] = master->tr_wc[tn][c];
  }
}

/* END fill_multidt_aux()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set multi dt auxiliary cells that are associated with  */
/* windows having a longer time-step than the current window to the  */
/* updated values at the end of the long time-step.                  */
/*-------------------------------------------------------------------*/
void fill_multidt_aux_2d(geometry_t *window,  /* Window geometry     */
                         window_t *windat,    /* Window data         */
                         double *vel,         /* Velocity array      */
                         double *ae /* Multi-dt cells velocity time t+1
                                       values                        */
  )
{
  int c, cc;                    /* Local sparse coordinate / counter */

  /* Auxiliary multi-dt cells                                        */
  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    if (c < window->enonS) {
      c = window->wsa[c];
      ae[cc] = vel[c];
    }
  }
}

/* END fill_multidt_aux_2d()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the time-step in each window based on the CFL      */
/* condition and order the windows from longest to shortest time-    */
/* steps.                                                            */
/*-------------------------------------------------------------------*/
void get_timesteps(geometry_t **window, /* Window geometry           */
                   window_t **windat,   /* Window data               */
                   win_priv_t **wincon, /* Window constants          */
                   int nwindows,        /* Number of windows         */
                   master_t *master     /* Model data structure      */
  )
{
  int wn, n, j, *nn;            /* Window counters                   */
  int c, c2, cc, e;             /* Sparse cponters                   */
  double top;                   /* Surface level                     */
  double d1, *d2, d3, *d4;      /* Dummies                           */
  double maxvel = 0.5;          /* Maximum 3D velocity expected      */
  double int_wave_speed = 2.0;  /* Maximum internal wave speed       */
  double sf = 0.9;              /* Safety factor                     */

  d2 = d_alloc_1d(nwindows + 1);
  d4 = d_alloc_1d(nwindows + 1);
  nn = i_alloc_1d(nwindows + 1);
  if (wincon[1]->stab & MULTI_DT)
    master->dt = 0.0;
  else if (!(wincon[1]->stab & MULTI_DT) && !windat[1]->dt)
    master->dt = 1e10;

  for (n = 1; n <= nwindows; n++) {

    wincon[n]->twin = i_alloc_1d(nwindows + 1);

    /*---------------------------------------------------------------*/
    /* If the time-step is already defined (non-zero) then use this  */
    /* time-step for all windows.                                    */
    if (windat[n]->dt) {
      for (c = 1; c <= nwindows; c++)
        nn[c] = c;
      d2[n] = windat[n]->dt / windat[n]->iratio;
      d4[n] = windat[n]->dt;
    }

    /*---------------------------------------------------------------*/
    /* Get the time-steps in each window from the CFL condition      */
    else {
      d2[n] = d4[n] = 1e10;

      nn[n] = n;
      for (cc = 1; cc <= window[n]->v3_t; cc++) {
        c = window[n]->w3_t[cc];
        c2 = window[n]->m2d[c];

	d1 = 0.0;
	for (j = 1; j <= window[n]->npe[c2]; j++) {
	  e = window[n]->c2e[j][c];
	  if (c == window[n]->e2c[e][0]) {
	    d1 += 1.0 / (window[n]->h1acell[e] * window[n]->h1acell[e]);
	  }
	}
        d3 = (1.0 / sqrt(d1)) / (maxvel + 2.0 * int_wave_speed);
        top = wincon[n]->sigma ? 0.0 : windat[n]->eta[c2];
        d1 =
          (1.0 / sqrt(d1)) / (maxvel +
                              2.0 * sqrt(g *
                                         fabs(top -
                                              window[n]->botz[c2]) *
                                         wincon[n]->Ds[c2]));
        if (d1 < d2[n])
          d2[n] = d1;
        if (d3 < d4[n])
          d4[n] = d3;
      }
      if (n == 4) {
        d4[n] /= 2.0;
        d2[n] /= 2.0;
      }
      if (n == 2) {
        d4[n] /= 4.0;
        d2[n] /= 4.0;
      }

      if (d2[n] > 1.0) {
        c = (int)(sf * d2[n]);
        cc = c / 10;
        if (c >= 10)
          d2[n] = (double)cc *10.0;
        else
          d2[n] = (double)cc;
      }
      if (d4[n] > 1.0) {
        c = (int)(sf * d4[n]);
        cc = c / 10;
        if (c >= 10)
          d4[n] = (double)cc *10.0;
        else
          d4[n] = (double)cc;
      }
      c = d4[n] / d2[n];
      if (c)
        d4[n] = d2[n] * c;

      if (wincon[n]->stab & MULTI_DT) {
        if (d4[n] > master->dt) {
          master->dt = d4[n];
          master->iratio = (int)(d4[n] / d2[n]);
        }
      } else {
        if (d4[n] < master->dt) {
          master->dt = d4[n];
          master->iratio = (int)(d4[n] / d2[n]);
        }
      }
    }
  }
  if (!windat[1]->dt)
    sorts2(d4, d2, nn, nwindows);

  for (cc = 1; cc <= nwindows; cc++) {
    c = nwindows - cc + 1;
    wn = nn[c];
    windat[wn]->dt = master->dt;
    windat[wn]->iratio = master->iratio;
    windat[wn]->dts = d4[c];
    windat[wn]->dt2s = d2[c];
    /* windat[wn]->iratio=d4[c]/d2[c]; */

    for (n = 1; n <= nwindows; n++) {
      c = nwindows - n + 1;
      wincon[wn]->twin[n] = nn[c];
    }

    if (DEBUG("init_w") && (wincon[1]->stab & MULTI_DT)) {
      dlog("init_w", "3D time step, window %d = %4.1f, IRATIO=%d\n",
           wincon[wn]->twin[cc], windat[wn]->dts, windat[wn]->iratio);
      dlog("init_w", "2D time step, window %d = %4.1f\n",
           wincon[wn]->twin[cc], windat[wn]->dt2s);
    }
  }
  d_free_1d(d2);
  d_free_1d(d4);
  i_free_1d(nn);
}

/* END get_timesteps()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to find & save the local sparse coordinates of auxiliary  */
/* cells which are associated with windows having a longer time-step */
/* than a given window. Memory for tracer buffers is also allocated. */
/*-------------------------------------------------------------------*/
void timeaux(geometry_t **window,  /* Window geometry                */
             window_t **windat,    /* Window data                    */
             win_priv_t **wincon,  /* Window constants               */
             int nwindows          /* Number of windows              */
  )
{
  int c, cc;                    /* Sparse coordinate / counter       */
  int n;                        /* Window counter                    */
  int gc;                       /* Global sparse coordinate          */
  int ac;                       /* Auxiliary cell counter            */
  int wn;                       /* Window corresponding to gc        */

  /* Loop over all windows */
  for (n = 1; n <= nwindows; n++)
    window[n]->naux_t = 0;

  if (!(wincon[1]->stab & MULTI_DT))
    return;

  for (n = 1; n <= nwindows; n++) {
    window[n]->aux_t = NULL;
    ac = 0;

    /*---------------------------------------------------------------*/
    /* Count the number of auxiliary cells in window n which are     */
    /* associated with a window having a longer time-step than       */
    /* window n.                                                     */
    for (cc = window[n]->v3_t + 1; cc <= window[n]->n3_t; cc++) {
      c = window[n]->w3_t[cc];  /* Local auxiliary/ghost location    */
      gc = window[n]->wsa[c];   /* Global location                   */
      wn = geom->fm[gc].wn;     /* Window corresponding to gc        */
      /* Note : ghost cells have a corresponding window = 0.         */
      if (wn && wn != n && windat[wn]->dts > windat[n]->dts)
        ac++;
    }
    if (DEBUG("init_w") && ac)
      dlog("init_w", "Window %d = %d multi time-step auxiliary cells\n", n,
           ac);

    /*---------------------------------------------------------------*/
    /* Allocate meory for long time-step auxiliary cells             */
    window[n]->naux_t = ac;
    window[n]->aux_t = i_alloc_1d(ac + 1);
    window[n]->taux_t = d_alloc_1d(ac + 1);
    windat[n]->tr_wc_as = d_alloc_2d(ac + 1, windat[n]->ntr);
    windat[n]->tr_wc_ae = d_alloc_2d(ac + 1, windat[n]->ntr);
    windat[n]->eta_as = d_alloc_1d(ac + 1);
    windat[n]->eta_ae = d_alloc_1d(ac + 1);
    windat[n]->u1av_as = d_alloc_1d(ac + 1);
    windat[n]->u1av_ae = d_alloc_1d(ac + 1);
    windat[n]->u2av_as = d_alloc_1d(ac + 1);
    windat[n]->u2av_ae = d_alloc_1d(ac + 1);

    /*---------------------------------------------------------------*/
    /* Fill the long time-step auxiliary cell vectors with the local */
    /* sparse coordinates in window n and save the longer time-step  */
    /* corresponding to the auxiliary cell window to array.          */
    ac = 1;
    for (cc = window[n]->v3_t + 1; cc <= window[n]->n3_t; cc++) {
      c = window[n]->w3_t[cc];
      gc = window[n]->wsa[c];
      wn = geom->fm[gc].wn;
      if (wn && wn != n && windat[wn]->dts > windat[n]->dts) {
        window[n]->aux_t[ac] = c;
        window[n]->taux_t[ac] = windat[wn]->dts;
        ac++;
      }
    }
  }
}

/* END timeaux()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sorts an array into increasing order                              */
/*-------------------------------------------------------------------*/
void sorts2(double *a, double *b, int *c, int n)
{
  int i, j;
  for (i = 1; i <= n; ++i)
    for (j = n; i < j; --j)
      orders2(&a[j - 1], &a[j], &b[j - 1], &b[j], &c[j - 1], &c[j]);
}

void orders2(double *p, double *q, double *r, double *s, int *i, int *j)
{
  double t1;
  int t2;
  if (*p > *q) {
    t1 = *p;
    *p = *q;
    *q = t1;
    t1 = *r;
    *r = *s;
    *s = t1;
    t2 = *i;
    *i = *j;
    *j = t2;
  }
}

/* END sorts2()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to print an array in ascii format                         */
/*-------------------------------------------------------------------*/
void prints(double *A,          /* Array to print                    */
            char *fname,        /* File name                         */
            int xlim, int ylim, int k, double scale)
{
  int i, j, c;
  FILE *fp;
  double fillv = 0.0;

  fp = fopen(fname, "w");
  if (fp == NULL)
    hd_quit("prints: Can't open file %s\n", fname);

  for (j = ylim - 1; j >= 0; j--) {
    for (i = 0; i < xlim; i++) {
      c = geom->map[k][j][i];
      if (c)
        fprintf(fp, "%8.3e ", A[c] * scale);
      else
        fprintf(fp, "%8.3e ", fillv);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

/* END prints()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-orders a cells to process vector                               */
/* mode = 1 : order from top to bottom                               */
/* mode = 2 : order from bottom to top                               */
/*-------------------------------------------------------------------*/
void reorder_cells(geometry_t *window, /* Processing window          */
		   int *ncells,        /* Re-ordered cells           */
		   int *cells,         /* Cells to process (ctp)     */
		   int vc,             /* Number of ctp              */
		   int vcs,            /* Surface number of ctp      */
		   int *bot,           /* Bottom vector              */
		   int mode            /* Re-ordering method         */
		   )
{
  int c, cc, cs, cb, nc;

  memcpy(ncells, cells, (window->enon + 1) * sizeof(int));
  nc = 1;

  /* Order from the surface to the bottom, for all i then j          */
  if(mode == 1) {
    for (cc = 1; cc <= vcs; cc++) {
      c = cells[cc];
      cb = bot[cc];
      while (c != cb) {
	ncells[nc] = c;
	nc++;
	c = window->zm1[c];
      }
      ncells[nc] = c;
      nc++;
    }
  }
  /* Order from the bottom to the surface, for all i then j          */
  else if(mode == 2) {
    for (cc = 1; cc <= vcs; cc++) {
      c = bot[cc];
      cs = cells[cc];
      while (c != cs) {
	ncells[nc] = c;
	nc++;
	c = window->zp1[c];
      }
      ncells[nc] = c;
      nc++;
    }
  }

  if(nc - 1 != vc)
    hd_warn("reorder_cells : different number of re-ordered cells (%d != %d)\n",
	    vc, nc - 1);
}

/* END reorder_cells()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the cells to exclude for particular processes   */
/*-------------------------------------------------------------------*/
void get_process_exclude(parameters_t *params, /* Params structure   */
			 geometry_t *geom,     /* Global geometry    */
			 geometry_t *window    /* Window geometry    */
			 )
{
  int cc, c, cs, i, j;

  window->ncbgc = 0;
  window->ncsed = 0;
  window->ncwave = 0;
  window->nctran = 0;
  window->nctrst = 0;
  if (params->prex == 0) return;

  for (cc = 1; cc <= params->prex; cc++) {
    i = params->prxi[cc];
    j = params->prxj[cc];
    c = geom->map[geom->nz - 1][j][i];
    if (geom->fm[c].wn == window->wn) {
      if (params->prxf[cc] & EX_BGC) window->ncbgc++;
      if (params->prxf[cc] & EX_SED) window->ncsed++;
      if (params->prxf[cc] & EX_WAVE) window->ncwave++;
      if (params->prxf[cc] & EX_TRAN) window->nctran++;
      if (params->prxf[cc] & EX_TRST) window->nctrst++;
    }
  }
  if (window->ncbgc)
    window->cbgc = i_alloc_1d(window->ncbgc + 1);
  if (window->ncsed)
    window->csed = i_alloc_1d(window->ncsed + 1);
  if (window->ncwave)
    window->cwave = i_alloc_1d(window->ncwave + 1);
  if (window->nctran)
    window->ctran = i_alloc_1d(window->nctran + 1);
  if (window->nctrst)
    window->ctrst = i_alloc_1d(window->nctrst + 1);
  window->ncbgc = window->ncsed = window->ncwave = 1;
  window->nctran = window->nctrst = 1;

  for (cc = 1; cc <= params->prex; cc++) {
    i = params->prxi[cc];
    j = params->prxj[cc];
    c = geom->map[geom->nz - 1][j][i];
    if (geom->fm[c].wn == window->wn) {
      cs = geom->fm[c].sc;
      if (params->prxf[cc] & EX_BGC)
	window->cbgc[window->ncbgc++] = cs;
      if (params->prxf[cc] & EX_SED)
	window->csed[window->ncsed++] = cs;
      if (params->prxf[cc] & EX_WAVE)
	window->cwave[window->ncwave++] = cs;
      if (params->prxf[cc] & EX_TRAN)
	window->ctran[window->nctran++] = cs;
      if (params->prxf[cc] & EX_TRST)
	window->ctrst[window->nctrst++] = cs;
    }
  }
  window->ncbgc--;
  window->ncsed--;
  window->ncwave--;
  window->nctran--;
  window->nctrst--;
}

/* END get_process_exclude()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
 /* Set the explicit masp for the window                             */
/*-------------------------------------------------------------------*/
void get_inner_exmap(geometry_t *geom,  /* Global geometry           */
		     geometry_t *window /* Window geometry           */
		     )
{
  int cc, c, cs, cd;
  int *mask = geom->emmask;

  /* Get the inner explicit maps for this window                     */
  if (geom->emmask) {

    window->neim = window->neimS = 0;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      if (mask[c] && geom->fm[c].wn == window->wn) {
	window->neim++;
	if (c <= geom->enonS)
	  window->neimS++;
      }
    }
    window->eims = i_alloc_1d(window->neim + 1);
    window->eimd = i_alloc_1d(window->neim + 1);
    memset(window->eims, 0, window->neim * sizeof(int));
    memset(window->eims, 0, window->neim * sizeof(int));
    window->neim = 1;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      if (mask[c] && geom->fm[c].wn == window->wn) {
	cs = geom->fm[c].sc;
	cd = geom->fm[mask[c]].sc;
	window->eims[window->neim] = cs;
	window->eimd[window->neim] = cd;
	window->neim++;
      }
    }
    window->neim--;
  }
}

/* END get_inner_exmap()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to make a mask array for cells to exclude                 */
/*-------------------------------------------------------------------*/
void process_cell_mask(geometry_t *window, int *cells, int ncells)
{
  win_priv_t *wincon = window->wincon;
  int c, cc;
  memset(wincon->c2, 0, window->sgsizS * sizeof(char));
  for (cc = 1; cc <= ncells; cc++) {
    c = cells[cc];
    wincon->c2[c] = 1;
  }
}

/* END process_cell_mask()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns a local coordinate of level k given the local surface     */
/* coordinate.                                                       */
/*-------------------------------------------------------------------*/
int get_local_sur(geometry_t *window, int cl, int ks)
{
  int k, c = cl;

  for (k = window->nz - 1; k > ks; k--) c = window->zm1[c];
  return(c);
}

/* END get_local_sur()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check for NaNs prior to the run start                  */
/*-------------------------------------------------------------------*/
void nan_check(geometry_t **window,   /* Window geometry             */
	       window_t **windat,     /* Window data                 */
	       win_priv_t **wincon,   /* Window constants            */
	       int nwindows           /* Number of windows           */
	       )
{
  int n, c, cc;

  for (n = 1; n <= nwindows; n++) {
    /* eta and atmospheric pressure                                  */
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (isnan(windat[n]->eta[c]))
	hd_warn("NaN found in eta, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
      if (isnan(windat[n]->patm[c]))
	hd_warn("NaN found in patm, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
    }
    /* Wind and 2D velocity                                          */
    for (cc = 1; cc <= window[n]->b2_e1; cc++) {
      c = window[n]->w2_e1[cc];
      if (isnan(windat[n]->u1av[c]))
	hd_warn("NaN found in u1av, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
      if (isnan(windat[n]->wind1[c]))
	hd_warn("NaN found in wind1, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
    }
    /* T/S                                                           */
    for (cc = 1; cc <= window[n]->b3_t; cc++) {
      c = window[n]->w3_t[cc];
      if (isnan(windat[n]->temp[c]))
	hd_warn("NaN found in temp, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
      if (isnan(windat[n]->sal[c]))
	hd_warn("NaN found in salt, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    /* 3D velocity                                                   */
    for (cc = 1; cc <= window[n]->b3_e1; cc++) {
      c = window[n]->w3_e1[cc];
      if (isnan(windat[n]->u1[c]))
	hd_warn("NaN found in u1, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    for (cc = window[n]->b3_e1 + 1; cc <= window[n]->a3_e1; cc++) {
      c = window[n]->w3_e1[cc];
      if (isnan(windat[n]->u1[c]))
	hd_warn("NaN found in u1, window %d at auxiliary (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
  }
}

/* END nan_check()                                                   */
/*-------------------------------------------------------------------*/

int eiw(geometry_t *geom, int e, int *mode)
{
  int c1 = geom->e2c[e][0];
  int c2 = geom->e2c[e][1];
  int wn1 = geom->fm[c1].wn;
  int wn2 = geom->fm[c2].wn;

  if (wn1 > 0 && wn2 > 0 && wn1 == wn2) {
    *mode = W_SAME;
    return(wn1);
  } else if (wn1 > 0 && wn2 > 0 && wn1 != wn2) {
    *mode = W_DIFF;
    return(wn1);
  } else if (wn1 == 0 && wn2 > 0) {
    *mode = W_GHOST;
    return(wn2);
  } else if (wn2 == 0 && wn1 > 0) {
    *mode = W_GHOST;
    return(wn1);
  } else
    hd_quit("Can't assign window to edge %d\n", e);
}


int eic(geometry_t *geom, int e, int wn, int *mode)
{
  int c1 = geom->e2c[e][0];
  int c2 = geom->e2c[e][1];
  int wn1 = geom->fm[c1].wn;
  int wn2 = geom->fm[c2].wn;

  if (wn1 > 0 && wn2 > 0 && wn1 == wn2 && wn1 == wn) {
    *mode = geom->e2e[e][0];
    return(c1);
  } else if (wn1 > 0 && wn2 > 0 && wn1 != wn2) {
    if (wn1 == wn) {
      *mode = geom->e2e[e][0];
      return(c1);
    }
    if (wn2 == wn) {
      *mode = geom->e2e[e][1];
      return(c2);
    }
  } else if (wn1 == 0 && wn2 > 0 && wn2 == wn) {
    *mode = geom->e2e[e][1];
    return(c2);
  } else if (wn2 == 0 && wn1 > 0 && wn1 == wn) {
    *mode = geom->e2e[e][0];
    return(c1);
  }
}

