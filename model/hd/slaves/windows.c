/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/slaves/windows.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: windows.c 6110 2019-02-15 05:43:04Z her127 $
 *
 */

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
void set_zoom_OBC(geometry_t *geom, open_bdrys_t **open, geometry_t **window);
int get_zcell(geometry_t *geom, open_bdrys_t *open, int cg, int zfe1, 
	      int zfe2, int mode);
int get_icellz(geometry_t *geom, int zf, int cg, int ci);
int get_zcell_cy(geometry_t *geom, geometry_t *window, open_bdrys_t *open,
		 int c, int mode);
void get_inner_exmap(geometry_t *geom, geometry_t *window);
void get_process_exclude(parameters_t *params, geometry_t *geom, geometry_t *window);
void tidalc_setup(geometry_t *geom, geometry_t *window, open_bdrys_t *open);
void wave_setup(geometry_t *geom, geometry_t *window, int mode);
void nan_check(geometry_t **window, window_t **windat, 
	       win_priv_t **wincon, int nwindows);
void local_map_build_d(int c, int cc, int wn, int *e1map, int *e2map, int *nmap,
		       int *rmap, int *wsa, int *ac, int nsaux, 
		       int *zmfe1, int *zmfe2);
void get_local_obc_a(int *vec, int nvec, int nvec2D, geometry_t **window, int nwindows); 		     
void set_reef_frac(master_t *master, geometry_t **window, window_t **windat, 
		     win_priv_t **wincon);
void window_cells_check(geometry_t *geom, int nwindows, int **ws2, int *wsizeS);

int zmode = 1;                /* zmode=0 : only zoom centers in zone */
                              /* zmode=1 : centers & faces in zone */

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
  int n;                        /* Window counters                   */
  int c, cc, c1;                /* Sparse location counters          */
  int *wsa;                     /* Window sparse array               */
  int **ws2;                    /* Window sparse array               */
  int *wsz, *wsz2D;             /* Size of 3D/2D window arrays       */
  int nwindows;                 /* Number of windows to make         */
  int *zmfe1;                   /* Zoom factor for coarsened grids   */
  int *zmfe2;                   /* Zoom factor for coarsened grids   */
  int readwin = 0;              /* =1 : read window map from file    */

  nwindows = geom->nwindows;
  if (DEBUG("init_w"))
    dlog("init_w", "Start making %d windows\n", nwindows);
  if (params->runmode & MANUAL && nwindows > 1 && strlen(params->win_file)) 
    readwin = 1;

  /*-----------------------------------------------------------------*/
  /* Allocate memory                                                 */
  geom->fm = (global_map_t *)malloc(sizeof(global_map_t) * geom->sgsiz);
  window = (geometry_t **)malloc(sizeof(geometry_t *) * (nwindows + 1));

  wsa = i_alloc_1d(geom->sgsiz);
  ws2 = i_alloc_2d(geom->sgsizS, nwindows + 1);
  wsz = i_alloc_1d(nwindows + 1);
  wsz2D = i_alloc_1d(nwindows + 1);
  zmfe1 = i_alloc_1d(nwindows + 1);
  zmfe2 = i_alloc_1d(nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Initialise the global to local map                              */
  /*UR-FIX even though 0 is not used initialise it from 0
   * some conditional statements depend on geom->fm[c].wn in
   * reorder_gl_map, despite the fact they shouldn't get there.	     */
  for (c = 0; c <= geom->sgnum; c++) {
    geom->fm[c].wn = 0;
    geom->fm[c].sc = 0;
    geom->fm[c].ac = 0;
  }
  geom->zoom = params->dozoom;
  for (n = 1; n <= nwindows; n++) {
    zmfe1[n] = 1;
    zmfe2[n] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Single window : in this case window arrays point to the global  */
  /* arrays.                                                         */
  if (nwindows == 1) {
    window[1] = window_alloc();
    window[1]->zoom = 0;
    window[1]->zmfe1 = 1;
    window[1]->zmfe2 = 1;
    window[1]->zmee1 = window[1]->zmee2 = 0;
    window[1]->zmme1 = window[1]->zmme2 = 0;
    window[1]->zoomf = 1;
    /* Get the sparse array sizes                                    */
    window[1]->enon = geom->sgnum;
    window[1]->ewet = geom->a3_t;
    window[1]->snon = geom->sgnum - geom->nbpt + 1;
    window[1]->enonS = geom->sgnumS;
    window[1]->ewetS = geom->a2_t;
    window[1]->snonS = geom->sgnumS - geom->nbptS + 1;
    window[1]->nz = geom->nz;
    window[1]->sednz = geom->sednz;
    window[1]->nwindows = nwindows;
    window[1]->wn = 1;
    /* Ghost sediment cells are only used in the semi-Lagrange       */
    /* scheme, which currently only runs with 1 window.              */
    window[1]->ngsed = geom->ngsed;
    window[1]->gsed_t = geom->gsed_t;
    window[1]->ised_t = geom->ised_t;

    /* Allocate memory for the window structure                      */
    alloc_geom(window[1], WINDOW_A);
    /* Point the window arrays to the global arrays                  */
    point_geom(window[1], geom);
    /* Define the diagnol maps from the directional maps             */
    for (c = 1; c <= geom->sgnum; c++) {
      window[1]->xmyp1[c] = geom->yp1[geom->xm1[c]];
      window[1]->xpym1[c] = geom->xp1[geom->ym1[c]];
    }
    /* Define the local - global maps                                */
    geom->zoom = geom->zoomf = geom->zmfe1 = geom->zmfe2 = 0;
    geom->zoomc = i_alloc_1d(geom->enonS + 1);
    for (c = 1; c <= geom->sgnum; c++) {
      window[1]->wsa[c] = c;
      geom->fm[c].sc = c;
      if (c <= geom->sgnumS)
        geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
    }
    /* Define the fm.wn map in window #1 for wet cells only          */
    for (cc = 1; cc <= geom->a3_t; cc++) {
      c = geom->w3_t[cc];
      geom->fm[c].wn = 1;
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
      read_windows(geom, window, params->win_file);
    } else {
      if (params->dozoom) {
	window_cells_zoom(geom, params, nwindows, zmfe1, zmfe2, ws2, wsz2D);
      }
      else {
	if (params->win_type & STRIPE_E2)
	  window_cells_linear_e2(geom, nwindows, ws2, wsz2D);
	else if (params->win_type & BLOCK_E1)
	  window_cells_block_e1(geom, nwindows, ws2, wsz2D, params->win_block);
	else if (params->win_type & BLOCK_E2)
	  window_cells_block_e2(geom, nwindows, ws2, wsz2D, params->win_block);
	else
	  window_cells_linear_e1(geom, nwindows, ws2, wsz2D);
      }
      window_cells_check(geom, nwindows, ws2, wsz2D);

      /* Get the global to local maps                                */
      get_gl_maps(geom, nwindows, ws2, wsz2D);

      if (DEBUG("init_w"))
	dlog("init_w", "  Global local map created OK\n");

      /*-------------------------------------------------------------*/
      /* Get the 3D sparse locations in each window based on the     */
      /* surface partitioning defined above.                         */
      for (n = 1; n <= nwindows; n++) {
	int cellf = 0;
	if (params->compatible & V1652 || params->zmfe1 == params->zmfe2)
	  cellf = 1;

	window[n] = window_alloc();
	window[n]->wn = n;
	window[n]->nwindows = nwindows;
	window[n]->zoom = params->dozoom;
	window[n]->zmfe1 = zmfe1[n];
	window[n]->zmfe2 = zmfe2[n];
	window[n]->zoomf = max(zmfe1[n], zmfe2[n]);
	window[n]->zmee1 = (int)(window[n]->zmfe1 / 2);
	window[n]->zmee2 = (int)(window[n]->zmfe2 / 2);
	if (params->dozoom & PRECOND) {
	  window[n]->zmme1 = window[n]->zmee1 + 1;
	  window[n]->zmme2 = window[n]->zmee2 + 1;
	} else {
	  window[n]->zmme1 = window[n]->zmfe1;
	  window[n]->zmme2 = window[n]->zmfe2;
	}
	window[n]->nz = geom->nz;
	window[n]->sednz = geom->sednz;
	
	/* Get the 3D local cells in window n                        */
	get_window_cells_h(geom, n, wsa, wsz, ws2[n], wsz2D[n]);

	/* Set the local wet array and local maps for window n       */
	get_local_maps(window[n], n, wsa, wsz[n], ws2[n], wsz2D[n], cellf);

	if (DEBUG("init_w"))
	  dlog("init_w", "  Local maps for window %d created OK\n", n);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Create the open boundary structures for each window             */
  OBC_build(geom->open, window, nwindows);

  if (!readwin) {
    /*---------------------------------------------------------------*/
    /* Set up the local processing vectors arrays for window wn      */
    get_local_wsa(geom->w3_t, geom->v3_t, geom->v2_t, window, nwindows, 
		  params->compatible, 0);
    get_local_wsa(geom->w3_e1, geom->v3_e1, geom->v2_e1, window, nwindows,
		  params->compatible, 1);
    get_local_wsa(geom->w3_e2, geom->v3_e2, geom->v2_e2, window, nwindows,
		  params->compatible, 2);

    if (DEBUG("init_w"))
      dlog("init_w", "  Local work arrays created OK\n\n");
  } else {
    if (!(params->compatible & V1957))
      get_local_obc_a(geom->w3_t, geom->v3_t, geom->v2_t, window, nwindows);
  }

  /*-----------------------------------------------------------------*/
  /* Get the master e1 and e2 zoomed maps                            */
  build_zoom_maps(geom, window);

  /*-----------------------------------------------------------------*/
  /* Set the free surface and bottom maps                            */
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,cc,c1),shared(window)
#endif
  for (n = 1; n <= nwindows; n++) {
    if (!readwin) {
      window[n]->sur_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->nsur_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->bot_t = i_alloc_1d(window[n]->a2_t + 1);
      window[n]->sur_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
      window[n]->bot_e1 = i_alloc_1d(window[n]->x2_e1 + 1);
      window[n]->sur_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
      window[n]->bot_e2 = i_alloc_1d(window[n]->x2_e2 + 1);
      window[n]->m2d = i_alloc_1d(window[n]->enon + 1);
      window[n]->s2i = i_alloc_1d(window[n]->enon + 1);
      window[n]->s2j = i_alloc_1d(window[n]->enon + 1);
      window[n]->s2k = i_alloc_1d(window[n]->enon + 1);

      surfbot_build(window[n], n, window[n]->nsur_t, window[n]->bot_t,
		    geom->bot_t, geom->v2_t, window[n]->w2_t,
		    window[n]->v2_t, window[n]->a2_t, 0);
      surfbot_build(window[n], n, window[n]->sur_e1, window[n]->bot_e1,
		    geom->bot_e1, geom->v2_e1, window[n]->w2_e1,
		    window[n]->v2_e1, window[n]->x2_e1, 1);
      surfbot_build(window[n], n, window[n]->sur_e2, window[n]->bot_e2,
		    geom->bot_e2, geom->v2_e2, window[n]->w2_e2,
		    window[n]->v2_e2, window[n]->x2_e2, 2);

      /*-------------------------------------------------------------*/
      /* Get the 3D - 2D map (including the sediment)                */
      for (cc = 1; cc <= window[n]->enon; cc++)
	window[n]->m2d[cc] = 0;
      for (cc = 1; cc <= window[n]->enonS; cc++) {
	c = cc;
	while (window[n]->zm1[c] && c != window[n]->zm1[c]) {
	  window[n]->m2d[c] = cc;
	  c = window[n]->zm1[c];
	}
	window[n]->m2d[c] = cc;
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

      /* Check for un-mapped cells, and define using lateral maps    */
      for (cc = 1; cc <= window[n]->enon; cc++) {
	if(window[n]->m2d[cc] == 0) {
	  if ((c = window[n]->xp1[window[n]->m2d[window[n]->xm1[cc]]])) 
	    window[n]->m2d[cc] = c;
	  else if ((c = window[n]->xm1[window[n]->m2d[window[n]->xp1[cc]]]))
	    window[n]->m2d[cc] = c;
	  else if ((c = window[n]->yp1[window[n]->m2d[window[n]->ym1[cc]]]))
	    window[n]->m2d[cc] = c;
	  else if ((c = window[n]->ym1[window[n]->m2d[window[n]->yp1[cc]]]))
	    window[n]->m2d[cc] = c;
	  /*
	  else
	    emstag(LDEBUG,"hd:windows:window_build",
		   "Can't find surface map for window %d cell %d\n", n, cc);
	  */
	}
      }

      /*-------------------------------------------------------------*/
      /* Shift the zoom cell faces to cell centers                   */
      if(zmode)
	zoom_shift(geom, window[n]);
      else
	set_zoom_m2d(window[n]);

      /*-------------------------------------------------------------*/
      /* Get the 3D - 1D map                                         */
      memset(window[n]->s2k, 0, window[n]->enon + 1);
      for (cc = 1; cc <= window[n]->enon; cc++) {
	c = window[n]->wsa[cc];
	window[n]->s2i[cc] = geom->s2i[c];
	window[n]->s2j[cc] = geom->s2j[c];
	window[n]->s2k[cc] = geom->s2k[c];
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the coordinate - index map (currently only used with      */
    /* sediments).                                                   */
    /*#if defined(HAVE_SEDIMENT_MODULE)*/
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
    /*#endif*/

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
    if(geom->sm_e1)
      window[n]->sm_e1 = win_vector_build(geom, window[n], n, geom->sm_e1, 1);
    if(geom->sm_e2)
      window[n]->sm_e2 = win_vector_build(geom, window[n], n, geom->sm_e2, 2);

    get_inner_exmap(geom, window[n]);
    get_process_exclude(params, geom, window[n]);
  }

  if (DEBUG("init_w"))
    dlog("init_w", "  Surface and bottom arrays created OK\n");

  /*-----------------------------------------------------------------*/
  /* Set the zoom interpolation structures in the global array */
  geom_interp(geom, window);
  geom_interp_e1(geom, window);
  geom_interp_e2(geom, window);

  /*-----------------------------------------------------------------*/
  /* Free memory */
  i_free_1d(wsa);
  i_free_2d(ws2);
  i_free_1d(wsz);
  i_free_1d(wsz2D);

  /*-----------------------------------------------------------------*/
  /* Print window information */
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
  /* Free memory not required in future computations */
  for (n = 1; n <= nwindows; n++) {
    /* This map is of size geom->sgsiz in each window with only */
    /* window[n]->enon locations used.  */
    if (nwindows > 1)
      free((global_map_t *)window[n]->fm);
  }
}

/* END window_build()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create an ascii file of window partitions              */
/*-------------------------------------------------------------------*/
void write_windows(geometry_t *geom,  /* Sparse global geometery */
                   unsigned long ***flag  /* Flags structure */
  )
{
  int c;                        /* Sparse counters */
  int n;                        /* Window counter */
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
  /*
     op=fopen("zoom.out","w"); for(i=0; i<geom->nfe1; i++) for(j=0;
     j<geom->nfe2; j++) { aa[j][i]=0; if(!(flag[k][j][i]&(SOLID|OUTSIDE)))
     { c=geom->map[k][j][i]; n=geom->zoomc[c]; aa[j][i]=n; } }
     for(j=geom->nfe2-1; j>=0; j--) { for(i=0; i<geom->nfe1; i++){
     fprintf(op,"%d ",aa[j][i]); } fprintf(op,"\n"); } fclose(op); */
  s_free_2d(aa);
}

/* END write_windows()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to partition the surface layer into windows if the domain */
/* is divided into nwindows consecutive rows of cells. Currently the */
/* windows are partitioned via a loop through the grid in the k,i,j  */
/* direction and consecutively assigns cells to windows.             */
/*-------------------------------------------------------------------*/
void window_cells_linear_e1(geometry_t *geom, /* Global geometery    */
			    int nwindows, /* Number of windows       */
			    int **ws2,   /*2D wet cells in window wn */
			    int *wsizeS  /* Number of 2D wet cells   */
  )
{
  int c, cc, n;
  int mode = 0;
  int sm = 0;

  /* Initialise */
  geom->zmfe1 = geom->zmfe2 = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  geom->zbe1 = i_alloc_1d(geom->enon + 1);
  geom->zme1 = i_alloc_1d(geom->enon + 1);
  geom->zbe2 = i_alloc_1d(geom->enon + 1);
  geom->zme2 = i_alloc_1d(geom->enon + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  for (c = 1; c <= geom->enon; c++) {
    geom->zbe1[c] = 1;
    geom->zme1[c] = 1;
    geom->zbe2[c] = 1;
    geom->zme2[c] = 1;
  }

  /* Get the number of cells in this window */
  c = (int)ceil(geom->a2_t / nwindows);
  for (n = 1; n < nwindows; n++)
    wsizeS[n] = c;
  wsizeS[nwindows] = geom->a2_t - (nwindows - 1) * c;

  /* If window sizes were specified, use these */
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

  /* If window sizes are explicitly specified, use these */
  for (n = 1; n <= nwindows; n++) {
    if (geom->nwn[n - 1]) {
      wsizeS[n] = geom->nwn[n - 1];
      mode = 1;
    }
    sm += wsizeS[n];
  }

  /* Save the surface cell locations in ws2 */
  if (mode) {
    int cl;
    int *mask = i_alloc_1d(geom->enonS + 1);
    memset(mask, 0, (geom->enonS + 1) * sizeof(int));
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
    int nn, cbc, cb;
    c = 1;
    for (n = 1; n <= nwindows; n++) {
      for (cc = 1; cc <= wsizeS[n]; cc++) {
	/* Don't allow window transition to occur on OBCs          */
	if (cc == wsizeS[n]) {
	  for (nn = 0; nn < geom->nobc; nn++) {
	    open_bdrys_t *open = geom->open[nn];
	    if (open->type == U1BDRY && ANYf(geom->wsa[c], open->obc_t, open->no2_t)) {
	      wsizeS[n]+=1;
	      wsizeS[nwindows]-=1;
	    }	    
	    if (open->type == U1BDRY && ANYf(geom->wsa[c], open->oi1_t, open->no2_t)) {
	      wsizeS[n]+=2;
	      wsizeS[nwindows]-=2;
	    }	    
	  }
	}
	ws2[n][cc] = geom->wsa[c];
	c++;
      }
    }
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
  geom->zmfe1 = geom->zmfe2 = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  geom->zbe1 = i_alloc_1d(geom->enon + 1);
  geom->zme1 = i_alloc_1d(geom->enon + 1);
  geom->zbe2 = i_alloc_1d(geom->enon + 1);
  geom->zme2 = i_alloc_1d(geom->enon + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  for (c = 1; c <= geom->enon; c++) {
    geom->zbe1[c] = 1;
    geom->zme1[c] = 1;
    geom->zbe2[c] = 1;
    geom->zme2[c] = 1;
  }

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
  geom->zmfe1 = geom->zmfe2 = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  geom->zbe1 = i_alloc_1d(geom->enon + 1);
  geom->zme1 = i_alloc_1d(geom->enon + 1);
  geom->zbe2 = i_alloc_1d(geom->enon + 1);
  geom->zme2 = i_alloc_1d(geom->enon + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  for (c = 1; c <= geom->enon; c++) {
    geom->zbe1[c] = 1;
    geom->zme1[c] = 1;
    geom->zbe2[c] = 1;
    geom->zme2[c] = 1;
  }

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
  geom->zmfe1 = geom->zmfe2 = 0;
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  geom->zbe1 = i_alloc_1d(geom->enon + 1);
  geom->zme1 = i_alloc_1d(geom->enon + 1);
  geom->zbe2 = i_alloc_1d(geom->enon + 1);
  geom->zme2 = i_alloc_1d(geom->enon + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZN | ZC | ZE1 | ZE2);
  }
  for (c = 1; c <= geom->enon; c++) {
    geom->zbe1[c] = 1;
    geom->zme1[c] = 1;
    geom->zbe2[c] = 1;
    geom->zme2[c] = 1;
  }

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
/* Routine to get the cells in each window if the domain is divided  */
/* into nwindows consecutive rows of cells with a section of the     */
/* domain zoomed (i.e. coarser resolution).                          */
/* The file zoom_cells contains a list of cell centers locations for */
/* the coarse grid. The zoom factor must be an odd number (3, 5 etc) */
/*-------------------------------------------------------------------*/
void window_cells_zoom(geometry_t *geom,     /* Global geometery     */
		       parameters_t *params, /* Parameter structure  */
                       int nwindows,         /* Number of windows    */
                       int *zmfe1,           /* e1 zoom factor       */
                       int *zmfe2,           /* e2 zoom factor       */
                       int **ws2,       /* 2D wet cells in window wn */
                       int *wsizeS      /* Size of ws2               */
  )
{
  int c, cc, zm1;               /* Sparse coordinate / counter */
  int c1, c2, c3, g1, g2;       /* Sparse coordinates */
  int n;                        /* Window counter */
  int wc, zc;                   /* Cell counter */
  int i, j;                     /* Counters */
  int zi, zj;                   /* Cartesian zoom zone coordinates */
  int zoomw = nwindows;         /* Window to use as the zoom zone */
  int zmwfe1;                   /* Zoom factor for the zoom window */
  int zmwfe2;                   /* Zoom factor for the zoom window */
  int nzme1, nzme2;             /* Number of times to map to get a face */
  int *ze1, *ze2;               /* Zoom factors */

  /* Initialize                                                      */
  if (params->zoomf > 1) {
    zmwfe1 = geom->zmfe1 = params->zoomf;
    zmwfe2 = geom->zmfe2 = params->zoomf;
    geom->zoomf = max(geom->zmfe1, geom->zmfe2);
  } else {
    zmwfe1 = geom->zmfe1 = params->zmfe1;
    zmwfe2 = geom->zmfe2 = params->zmfe2;
    geom->zoomf = params->zoomf = max(geom->zmfe1, geom->zmfe2);
  }
  nzme1 = (int)(zmwfe1 / 2);
  nzme2 = (int)(zmwfe2 / 2);

  /* Initialise */
  zmfe1[0] = zmfe2[0] = 0;
  for (n = 1; n <= nwindows; n++) {
    wsizeS[n] = 0;
    zmfe1[n] = 1;
    zmfe2[n] = 1;
    if (n == zoomw) {
      zmfe1[n] = zmwfe1;
      zmfe2[n] = zmwfe2;
    }
  }
  geom->zoomc = i_alloc_1d(geom->enonS + 1);
  geom->zbe1 = i_alloc_1d(geom->enon + 1);
  geom->zme1 = i_alloc_1d(geom->enon + 1);
  geom->zbe2 = i_alloc_1d(geom->enon + 1);
  geom->zme2 = i_alloc_1d(geom->enon + 1);
  ze1 = i_alloc_1d(geom->enon + 1);
  ze2 = i_alloc_1d(geom->enon + 1);
  for (c = 1; c <= geom->enonS; c++) {
    geom->zoomc[c] = (ZC | ZE1 | ZE2);
    ze1[c] = ze2[c] = 1;
    for (n = 1; n <= nwindows; n++)
      ws2[n][c] = 0;
  }
  for (c = 1; c <= geom->enon; c++) {
    geom->zbe1[c] = 1;
    geom->zme1[c] = 1;
    geom->zbe2[c] = 1;
    geom->zme2[c] = 1;
  }

  /* Set the codes for cell centers and faces in the zoomed grid */
  /* cell, count the number of non-zoomed cells in this window and */
  /* save the surface non-zoomed cell locations in ws2.  */
  wc = zc = 1;
  for (cc = 1; cc <= params->nzoom; cc++) {
    zi = params->zci[cc];
    zj = params->zcj[cc];
    c = geom->map[geom->nz - 1][zj][zi];
    if ((c > 0 && c <= geom->ewetS) ||
	(c > geom->enonS && c <= geom->ewet)) {
      g1 = g2 = c3 = c;
      for (n = 1; n <= nzme1; n++) {
        g1 = geom->xm1[g1];     /* Zoom cell e1 location */
        c3 = geom->xm1[c3];     /* Zoom cell corner location */
      }
      for (n = 1; n <= nzme2; n++) {
        g2 = geom->ym1[g2];     /* Zoom cell e2 location */
        c3 = geom->ym1[c3];     /* Zoom cell corner location */
      }
      /* Set the window number for all cells in the zoom zone */
      c2 = c3;
      for (j = 1; j <= zmwfe2; j++) {
        c1 = c2;
        for (i = 1; i <= zmwfe1; i++) {
          if (zmode) {
            ws2[zoomw][wc] = c1;
            wc++;
            geom->fm[c1].wn = zoomw;
          } else {
            geom->fm[c1].wn = ZF;
          }
          geom->zoomc[c1] = ZN;
	  ze1[c1] = zmwfe1;
	  ze2[c1] = zmwfe2;
          c1 = geom->xp1[c1];
          zc++;
        }
        c2 = geom->yp1[c2];
      }
      /* Set the codes for the e1 and e2 faces */
      if (zmode) {
        geom->zoomc[c] = ZC;
        geom->zoomc[g1] = (g1 == c) ? (ZC | ZE1) : ZE1;
        geom->zoomc[g2] = (g2 == c) ? (ZC | ZE2) : ZE2;
      } else {
        ws2[zoomw][wc] = c;
        wc++;
        geom->zoomc[c] = (ZC | ZE1 | ZE2);
        geom->fm[c].wn = zoomw;
      }
    }
    /* Set the zoom codes for R_EDGE and F_EDGE open boundaries     */
    /* R_EDGE boundaries                                            */
    c1 = c;
    for (i = 1; i <= nzme1 + 1; i++)
      c1 = geom->xp1[c1];
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (ANYf(c1, open->obc_e1, open->no2_e1)) {
	geom->zoomc[c1] = ZE1;
	geom->fm[c1].wn = zoomw;
	ws2[zoomw][wc] = c1;
	wc++;
	zc++;
	c2 = c1;
	for (j = 1; j <= nzme2; j++) {
	  c1 = geom->yp1[c1];
	  c2 = geom->ym1[c2];
	  geom->zoomc[c1] = geom->zoomc[c2] = ZN;
	  geom->fm[c1].wn = geom->fm[c2].wn = zoomw;
	  ws2[zoomw][wc] = c1;
	  wc++;
	  ws2[zoomw][wc] = c2;
	  wc++;
	  zc += 2;
	}
      }
    }
    /* F_EDGE boundaries                                            */
    c2 = c;
    for (j = 1; j <= nzme2 + 1; j++)
      c2 = geom->yp1[c2];
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      if (ANYf(c2, open->obc_e2, open->no2_e2)) {
	geom->zoomc[c2] = ZE2;
	geom->fm[c2].wn = zoomw;
	ws2[zoomw][wc] = c2;
	wc++;
	zc++;
	c1 = c2;
	for (i = 1; i <= nzme1; i++) {
	  c1 = geom->xp1[c1];
	  c2 = geom->xm1[c2];
	  geom->zoomc[c1] = geom->zoomc[c2] = ZN;
	  geom->fm[c1].wn = geom->fm[c2].wn = zoomw;
	  ws2[zoomw][wc] = c1;
	  wc++;
	  ws2[zoomw][wc] = c2;
	  wc++;
	  zc += 2;
	}
      }
    }
  }
  wc--;
  zc--;

  if (DEBUG("zoom") && (zc / (zmwfe1 * zmwfe2) != params->nzoom))
    dlog("zoom", "Warning : Non integral cell number in zoom zone\n");

  /* Divide the remaining cells outsize the zoom zone into nwindows */
  /* minus 1 windows. The last window is the zoom zone with wc */
  /* cells.  */
  cc = geom->a2_t - zc;
  c = (int)ceil(cc / (nwindows - 1));
  for (n = 1; n < nwindows - 1; n++)
    wsizeS[n] = c;
  wsizeS[nwindows - 1] = cc - (nwindows - 2) * c;
  wsizeS[zoomw] = wc;

  /* Save the surface cell locations in ws2 */
  n = wc = 1;
  if (n == zoomw)
    n++;
  for (c = 1; c <= geom->a2_t; c++) {
    if (n <= nwindows && !(geom->fm[c].wn)) {
      ws2[n][wc] = c;
      wc++;
      if (wc > wsizeS[n]) {
        n++;
        if (n == zoomw)
          n++;
        wc = 1;
      }
    }
  }

  /* Set the window number of sub-surface unused cells in the zoom */
  /* zone to -1.  */
  if (!zmode) {
    for (cc = 1; cc <= geom->a2_t; cc++) {
      c = cc;
      zm1 = geom->zm1[c];
      if (geom->fm[c].wn == ZF) {
        geom->fm[c].wn = zoomw;
        do {
          geom->fm[c].wn = zoomw;
          c = zm1;
          zm1 = geom->zm1[c];
        } while (c != zm1);
      }
    }
  }

  /* Set the zoom factors for merged cells. The merged cells have    */
  /* zoom factors of zmwfe1/zmwfe2, non-merged cells have factors of */
  /* 1, and the merged-regular boundaries must be adjusted.          */
  for (c = 1; c <= geom->enonS; c++) {
    g1 = g2 = c;
    geom->zbe1[c] = geom->zme1[c] = ze1[c];
    geom->zbe2[c] = geom->zme2[c] = ze2[c];
    if (params->dozoom & PRECOND) {
      for (zc = 1; zc <= ze1[c]; zc++) {
	g1 = geom->xp1[g1];
	g2 = geom->xm1[g2];
      }
      /* Merged - regular boundaries */
      if (ze1[c] > 1 && ze1[g1] == 1) geom->zbe1[c] -= nzme1;
      if (ze1[c] > 1 && ze1[g2] == 1) geom->zme1[c] -= nzme1; 
      /* Regular - merged boundaries */
      g1 = geom->xp1[c];
      g2 = geom->xm1[c];
      if (ze1[c] == 1 && ze1[g1] > 1) geom->zbe1[c] += nzme1;
      if (ze1[c] == 1 && ze1[g2] > 1) geom->zme1[c] += nzme1;

      g1 = g2 = c;
      for (zc = 1; zc <= ze2[c]; zc++) {
	g1 = geom->yp1[g1];
	g2 = geom->ym1[g2];
      }
      if (ze2[c] > 1 && ze2[g1] == 1) geom->zbe2[c] -= nzme2;
      if (ze2[c] > 1 && ze2[g2] == 1) geom->zme2[c] -= nzme2;
      g1 = geom->yp1[c];
      g2 = geom->ym1[c];
      if (ze2[c] == 1 && ze2[g1] > 1) geom->zbe2[c] += nzme2;
      if (ze2[c] == 1 && ze2[g2] > 1) geom->zme2[c] += nzme2;
    }
  }

  /* Set the sub-surface zoom factors */
  for (cc = 1; cc <= geom->enonS; cc++) {
    c = geom->zm1[cc];
    while (c != geom->zm1[c]) {
      geom->zbe1[c] = geom->zbe1[cc];
      geom->zme1[c] = geom->zme1[cc];
      geom->zbe2[c] = geom->zbe2[cc];
      geom->zme2[c] = geom->zme2[cc];
      c = geom->zm1[c];
    }
  }
  for (cc = 1; cc <= params->nzoom; cc++) {
    zi = params->zci[cc];
    zj = params->zcj[cc];
    c = geom->map[geom->nz - 1][zj][zi];
    if ((c > 0 && c <= geom->ewetS) ||
	(c > geom->enonS && c <= geom->ewet)) {
      geom->zoomc[c] |= ZZ;
    }
  }
  i_free_1d(ze1);
  i_free_1d(ze2);
}

/* END window_cells_zoom()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void window_cells_check(geometry_t *geom, /* Global geometery        */
			int nwindows, /* Number of windows           */
			int **ws2,    /* 2D wet cells in window wn   */
			int *wsizeS   /* Number of 2D wet cells      */
			)
{
  int wn, n, m, c, ci, co, cc;
  int *wmap, *cmap, *mask;
  int **ws, *wsz;

  wmap = i_alloc_1d(geom->enonS + 1);
  cmap = i_alloc_1d(geom->enonS + 1);
  ws = i_alloc_2d(geom->sgsizS, nwindows + 1);
  wsz = i_alloc_1d(nwindows + 1);

  for (n = 1; n <= nwindows; n++) {
    wsz[n] = wsizeS[n];
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];
      wmap[c] = n;
      cmap[c] = cc;
      ws[n][cc] = ws2[n][cc];
    }
    wsizeS[n] = 0.0;
  }

  /* Set cells interior and exterior to OBCs to be in the same       */
  /* window.                                                         */
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      ci = open->oi1_t[cc];
      wn = wmap[c];
      if (wmap[c] != wmap[ci]) wmap[ci] = wn;
    }
  }
  for (n = 0; n < geom->nobc; n++) {
    open_bdrys_t *open = geom->open[n];
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      co = open->omap[c];
      wn = wmap[c];
      if (wmap[c] != wmap[co]) wmap[co] = wn;
    }
  }

  /* Count the cells in each window and re-assign the coordinate     */
  for (cc = 1; cc <= geom->a2_t; cc++) {
    c = geom->wsa[cc];
    n = wmap[c];
    wsizeS[n]++;
    ws2[n][wsizeS[n]] = c;
  }

  i_free_1d(wmap);
  i_free_1d(cmap);
}

/* END window_cells_check()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the global to local maps for each window, i.e.     */
/* derive the mapping between a global sparse coordinate and the     */
/* corresponding window, and local sparse coordinate in that window. */
/*-------------------------------------------------------------------*/
void get_gl_maps(geometry_t *geom,  /* Sparse global geometery */
                 int nwindows,  /* Number of windows */
                 int **ws2,     /* Array of 2D wet cells in window wn */
                 int *wsizeS    /* Number of 2D wet cells in window wn */
  )
{
  int c, cc, c1, c2, n;         /* Sparse array counters */
  int *w2;                      /* Dummy 2D surface cell vector */

  for (n = 1; n <= nwindows; n++) {
    /* Allocate memory */
    w2 = i_alloc_1d(wsizeS[n] + 1);

    /* Get the global to local surface maps for the sub-surface */
    /* layers filling the array in the (i,j,k) direction.  */
    for (cc = 1; cc <= wsizeS[n]; cc++) {
      c = ws2[n][cc];           /* Global coordinate in window n */
      geom->fm[c].sc = cc;      /* Surface global to local map */
      geom->fm[c].wn = n;       /* Surface global window */
      w2[cc] = c;               /* Make a copy of ws2 */
    }

    /* Map all wet cells below the surface layer. This creates the */
    /* map in the (i,j,k) direction, however, the order is not */
    /* critical since geom->fm gets re-ordered in the routine */
    /* get_local_maps().  */
    c2 = wsizeS[n];
    while (c2) {
      for (c1 = 1; c1 <= wsizeS[n]; c1++) {
        c = w2[c1];             /* Surface global coordinate in window wn */
        if (c) {                /* if(not a sediment cell) : ws2[c1]=0 */
          c = geom->zm1[c];     /* Global location one layer down */
          w2[c1] = c;           /* Set the 2D local location to this cell */
          if (c != geom->zm1[c]) {  /* if(not a sediment cell) */
            geom->fm[c].sc = cc;  /* Set global-local map location */
            geom->fm[c].wn = n; /* Set the global-local map window */
            cc++;               /* Increment the local coordinate */
          } else {              /* if(sediment cell) */
            if (w2[c1])
              c2--;             /* Set this column as done */
            w2[c1] = 0;         /* Flag the local location sediment */
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
/* window partitioning. Fills the local sparse array in columns      */
/* first.                                                            */
/*-------------------------------------------------------------------*/
void get_window_cells(geometry_t *geom, /* Sparse global geometery */
                      int wn,   /* Window number */
                      int *wsa, /* Array of 3D wet cells in window wn */
                      int *wsize, /* Number of 3D wet cells in window wn */
                      int *ws2, /* Array of 2D wet cells in window wn */
                      int wsizeS  /* Number of 2D wet cells in window wn */
  )
{
  int c, cc, c1;                /* Sparse array counters */

  /* Initialise auxiliary maps for this window */
  for (c = 1; c <= geom->sgnum; c++)
    geom->fm[c].ac = 0;

  /* Get the global to local surface maps for the sub-surface layers */
  /* filling the array in the (i,j,k) direction.  */
  c1 = 1;                       /* 3D cell counter */
  wsize[wn] = 0;                /* Initialise number of 3D cells */
  for (cc = 1; cc <= wsizeS; cc++) {
    c = ws2[cc];                /* Global coordinate in window wn */
    while (c != geom->zm1[c]) {
      wsa[c1] = c;
      c1++;
      wsize[wn]++;
      c = geom->zm1[c];
    }
  }
}

/* END get_window_cells()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the 3D local cells in window n given the surface   */
/* window partitioning. Fills the local sparse array in horizontal   */
/* slabs first.                                                      */
/*-------------------------------------------------------------------*/
void get_window_cells_h(geometry_t *geom, /* Sparse global geometery */
                        int wn, /* Window number */
                        int *wsa, /* Array of 3D wet cells in window wn */
                        int *wsize, /* Number of 3D wet cells in window wn
                                     */
                        int *ws2, /* Array of 2D wet cells in window wn */
                        int wsizeS  /* Number of 2D wet cells in window wn
                                     */
  )
{
  int c, cc, c1, zm1;           /* Sparse array counters */
  int *layer;                   /* Vertical layer map */

  /* Initialise the vertical layer map */
  layer = i_alloc_1d(wsizeS + 1);
  memcpy(layer, ws2, (wsizeS + 1) * sizeof(int));
  layer[0] = 0;

  /* Initialise auxiliary maps for this window */
  for (c = 1; c <= geom->sgnum; c++)
    geom->fm[c].ac = 0;

  /* Get the global to local surface maps for the sub-surface layers */
  /* filling the array in the (i,j,k) direction.  */
  c1 = 1;                       /* 3D cell counter */
  wsize[wn] = 0;                /* Initialise number of 3D cells */
  while (layer[0] < wsizeS) {
    for (cc = 1; cc <= wsizeS; cc++) {
      c = layer[cc];            /* Global coordinate in window wn */
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
void get_local_maps(geometry_t *window, /* Window data structure */
                    int wn,     /* Window number */
                    int *wsa,   /* Array of 3D wet cells in window wn */
                    int wsize,  /* Number of 3D wet cells in window wn */
                    int *ws2,   /* Array of 2D wet cells in window wn */
                    int wsize2D,/* Number of 2D wet cells in window wn */
		    int cellf   /* Use isotropic map generation      */
  )
{
  int c, cc, c2D;               /* Sparse array counters */
  int ac;                       /* Auxiliary cell counter for the window */
  int cl, cu, cd, cs;           /* Local sparse coordinate */
  int *wsar;                    /* Reordered work sparse array */
  int *nxp1, *nxm1;             /* Spatial buffer arrays for i+1 and i-1 */
  int *nyp1, *nym1;             /* Spatial buffer arrays for j+1 and j-1 */
  int *nzp1, *nzm1;             /* Spatial buffer arrays for k+1 and k-1 */
  int *nxpym1, *nxmyp1;         /* Spatial buffer arrays for diagnols */
  int *nxmym1, *nxpyp1;         /* Spatial buffer arrays for diagnols */
  int *xpym1, *xmyp1;           /* Spatial buffer arrays for diagnols */
  int *xmym1, *xpyp1;           /* Spatial buffer arrays for diagnols */
  int *nypxm1, *ypxm1;
  int *one;

  /* int *nyp2,*nym2; *//* Spatial buffer arrays for j-1 and j-2 */
  /* int *nxp2,*nxm2; *//* Spatial buffer arrays for i-1 and i-2 */
  /* int *nzp2,*nzm2; *//* Spatial buffer arrays for k-1 and k-2 */

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the buffer arrays. Since it is not known at */
  /* this stage how many auxiliary cells are to be included in the */
  /* window arrays, set the size of these buffers to the sparse grid */
  /* size.  */
  wsar = i_alloc_1d(geom->sgsiz);
  nxp1 = i_alloc_1d(geom->sgsiz);
  nxm1 = i_alloc_1d(geom->sgsiz);
  nyp1 = i_alloc_1d(geom->sgsiz);
  nym1 = i_alloc_1d(geom->sgsiz);
  nzp1 = i_alloc_1d(geom->sgsiz);
  nzm1 = i_alloc_1d(geom->sgsiz);
  nxpym1 = i_alloc_1d(geom->sgsiz);
  nxmyp1 = i_alloc_1d(geom->sgsiz);
  nxmym1 = i_alloc_1d(geom->sgsiz);
  nxpyp1 = i_alloc_1d(geom->sgsiz);
  nypxm1 = i_alloc_1d(geom->sgsiz);
  one = i_alloc_1d(geom->sgsiz);

  /*
     nxp2 = i_alloc_1d(geom->sgsiz); nxm2 = i_alloc_1d(geom->sgsiz); nyp2
     = i_alloc_1d(geom->sgsiz); nym2 = i_alloc_1d(geom->sgsiz); nzp2 =
     i_alloc_1d(geom->sgsiz); nzm2 = i_alloc_1d(geom->sgsiz); */

  /* Create the global diagonal maps */
  xmyp1 = i_alloc_1d(geom->sgsiz);
  xpym1 = i_alloc_1d(geom->sgsiz);
  xmym1 = i_alloc_1d(geom->sgsiz);
  xpyp1 = i_alloc_1d(geom->sgsiz);
  ypxm1 = i_alloc_1d(geom->sgsiz);
  for (c = 1; c <= geom->sgnum; c++) {
    xmyp1[c] = geom->yp1[geom->xm1[c]];
    xpym1[c] = geom->xp1[geom->ym1[c]];
    xmym1[c] = geom->xm1[geom->ym1[c]];
    xpyp1[c] = geom->xp1[geom->yp1[c]];
    ypxm1[c] = geom->xm1[geom->yp1[c]];
    one[c] = 1;
  }

  /*-----------------------------------------------------------------*/
  /* Initialise the buffers. The local maps are made self mapping.  */
  /* The local wet work array is filled with global wet sparse */
  /* locations; auxiliary global locations are appended below.  */
  window->ewet = wsize;
  window->snon = wsize + 1;
  window->ewetS = wsize2D;
  window->snonS = wsize2D + 1;
  /*UR-FIX initialise from 0 - odd times that is used */
  for (c = 0; c <= geom->sgnum; c++) {
    wsar[c] = 0;
    nxp1[c] = c;
    nxm1[c] = c;
    nyp1[c] = c;
    nym1[c] = c;
    nzp1[c] = c;
    nzm1[c] = c;
    nxpym1[c] = c;
    nxmyp1[c] = c;
    nxmym1[c] = c;
    nxpyp1[c] = c;
    nypxm1[c] = c;
    /*
       nxp2[c] = c; nxm2[c] = c; nyp2[c] = c; nym2[c] = c; nzp2[c] = c;
       nzm2[c] = c; */
  }
  for (c = 1; c <= wsize2D; c++)
    wsar[c] = ws2[c];

  /*-----------------------------------------------------------------*/
  /* Get the work sparse array and related horizontal maps for the */
  /* surface layer. Any sparse locations encountered in the surface */
  /* layer are tagged as zero.  */
  c2D = 1;
  ac = wsize2D;
  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];
    if (c == ws2[c2D]) {
      if (cellf) {
	local_map_build_o(c, c2D, wn, geom->xp1, nxp1, nxm1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, geom->xm1, nxm1, nxp1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, geom->yp1, nyp1, nym1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, geom->ym1, nym1, nyp1, wsar, &ac, laux,
			  window->zoomf);
      } else {
	local_map_build(c, c2D, wn, geom->xp1, nxp1, nxm1, wsar, &ac, laux,
			geom->zbe1);
	local_map_build(c, c2D, wn, geom->xm1, nxm1, nxp1, wsar, &ac, laux,
			geom->zme1);
	local_map_build(c, c2D, wn, geom->yp1, nyp1, nym1, wsar, &ac, laux,
			geom->zbe2);
	local_map_build(c, c2D, wn, geom->ym1, nym1, nyp1, wsar, &ac, laux,
			geom->zme2);
      }
      /* Different ghost cells can be obtained using maps xmyp1 and  */
      /* ypxm1. The lateral maps, bpt, use the innermost map, i.e.   */
      /* xm1 for xmyp1 or yp1 for ypxm1. Ghost cells are not defined */
      /* for auxiliary cells, but are required (e.g for diffusion    */
      /* tensors). These ghost cells may be defined using the        */
      /* diagonals, but maps from both directions (xmyp1 and ypxm1)  */
      /* are required to set a unique lateral boundary condition.    */
      if (cellf) {
	local_map_build_o(c, c2D, wn, ypxm1, nypxm1, nxpym1, wsar, &ac, 1,
			  window->zoomf);
      } else {
	/* Take note of the order of the maps in this function       */
	local_map_build_d(c, c2D, wn, geom->xm1, geom->yp1, nypxm1, nxpym1, 
			  wsar, &ac, 1, geom->zme1, geom->zbe2);
      }
      /* Map the diagonals                                           */
      if (cellf) {
	local_map_build_o(c, c2D, wn, xpym1, nxpym1, nxmyp1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, xmyp1, nxmyp1, nxpym1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, xmym1, nxmym1, nxpyp1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, c2D, wn, xpyp1, nxpyp1, nxmym1, wsar, &ac, 1,
			  window->zoomf);
      } else {
	local_map_build_d(c, c2D, wn, geom->xp1, geom->ym1, nxpym1, nxmyp1, 
			  wsar, &ac, 1, geom->zme2, geom->zbe1);
	local_map_build_d(c, c2D, wn, geom->yp1, geom->xm1, nxmyp1, nxpym1, 
			  wsar, &ac, 1, geom->zbe2, geom->zme1);
	local_map_build_d(c, c2D, wn, geom->xm1, geom->ym1, nxmym1, nxpyp1, 
			  wsar, &ac, 1, geom->zme2, geom->zme1);
	local_map_build_d(c, c2D, wn, geom->xp1, geom->yp1, nxpyp1, nxmym1,
			  wsar, &ac, 1, geom->zbe2, geom->zbe1);
      }
      wsa[cc] = 0;
      c2D++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reorder the work array, wsar. This currently includes the wet */
  /* surface cells followed by the surface auxiliary cells. Append */
  /* all non-surface wet cells to the end of the sutface layer here. */
  window->enonS = ac;
  for (cc = 1; cc <= wsize; cc++) {
    if (wsa[cc]) {
      ac++;
      wsar[ac] = wsa[cc];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Reset the global-local map to that global sparse locations map */
  /* correctly to the re-ordered local coordinates. This can only be */
  /* performed once the number of auxiliary cells in the surface */
  /* layer are known.  */
  reorder_gl_map(wn, wsa, wsize, wsar, window->enonS);

  /*-----------------------------------------------------------------*/
  /* Get the vertical maps for the surface layer wet cells. These */
  /* should all map to other wet cells (downward maps), map to the */
  /* sediment (downward maps in a water column one cell deep) or be */
  /* self mapping (upward maps). Only downward maps to the sediment */
  /* will increment ac; all other local sparse locations are already */
  /* defined.  */
  for (cc = 1; cc <= wsize2D; cc++) {
    c = wsar[cc];
    if (cellf) {
      local_map_build_o(c, cc, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux, 1);
      local_map_build_o(c, cc, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux, 1);
    } else {
      local_map_build(c, cc, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux, one);
      local_map_build(c, cc, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux, one);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the work sparse array and related maps for the non-surface */
  /* layer.  */
  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];
    if (c) {
      cl = geom->fm[c].sc;
      if (cellf) {
	local_map_build_o(c, cl, wn, geom->xp1, nxp1, nxm1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, cl, wn, geom->xm1, nxm1, nxp1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, cl, wn, geom->yp1, nyp1, nym1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, cl, wn, geom->ym1, nym1, nyp1, wsar, &ac, laux,
			  window->zoomf);
	local_map_build_o(c, cl, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux,
			  1);
	local_map_build_o(c, cl, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux,
			  1);
      } else {
	local_map_build(c, cl, wn, geom->xp1, nxp1, nxm1, wsar, &ac, laux,
			geom->zbe1);
	local_map_build(c, cl, wn, geom->xm1, nxm1, nxp1, wsar, &ac, laux,
			geom->zme1);
	local_map_build(c, cl, wn, geom->yp1, nyp1, nym1, wsar, &ac, laux,
			geom->zbe2);
	local_map_build(c, cl, wn, geom->ym1, nym1, nyp1, wsar, &ac, laux,
			geom->zme2);
	local_map_build(c, cl, wn, geom->zp1, nzp1, nzm1, wsar, &ac, laux,
			one);
	local_map_build(c, cl, wn, geom->zm1, nzm1, nzp1, wsar, &ac, laux,
			one);
      }
    }
  }

  /* Get the 3D diagonal maps and sparse locations of diagonals that */
  /* were not captured by the lateral maps (i.e. bottom right and */
  /* top left outside corners). These auxiliary cells only need to */
  /* be included one cell laterally away from the wet cell, since */
  /* no higher order approximations are required with these maps.  */
  /* The xmyp1 and xpym1 maps are used to get u2au1 and u1au2; the */
  /* xmym1 auxiliary cell is used in the horizontal diffusion tensor */
  /* calculation (note: a map of xmym1 is not saved here, the */
  /* auxiliary cell in this direction is only identified and */
  /* included in the window and is accessed via ym1[xm1[c]]).  */
  /* The xpyp1 map is used for Smagorinsky diffusion. */
  for (cc = 1; cc <= wsize; cc++) {
    c = wsa[cc];
    if (c) {
      cl = geom->fm[c].sc;
      if (cellf) {
	local_map_build_o(c, cl, wn, ypxm1, nypxm1, nxpym1, wsar, &ac, 1,
			  window->zoomf);

	local_map_build_o(c, cl, wn, xpym1, nxpym1, nxmyp1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, cl, wn, xmyp1, nxmyp1, nxpym1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, cl, wn, xmym1, nxmym1, nxpyp1, wsar, &ac, 1,
			  window->zoomf);
	local_map_build_o(c, cl, wn, xpyp1, nxpyp1, nxmym1, wsar, &ac, 1,
			  window->zoomf);
      } else {
	local_map_build_d(c, cl, wn, geom->xm1, geom->yp1, nypxm1, nxpym1,
			  wsar, &ac, 1,geom->zme1, geom->zbe2);
	
	local_map_build_d(c, cl, wn, geom->xp1, geom->ym1, nxpym1, nxmyp1,
			  wsar, &ac, 1, geom->zme2, geom->zbe1);
	local_map_build_d(c, cl, wn, geom->yp1, geom->xm1, nxmyp1, nxpym1,
			  wsar, &ac, 1, geom->zbe2, geom->zme1);
	local_map_build_d(c, cl, wn, geom->xm1, geom->ym1, nxmym1, nxpyp1,
			  wsar, &ac, 1, geom->zme2, geom->zme1);
	local_map_build_d(c, cl, wn, geom->xp1, geom->yp1, nxpyp1, nxmym1,
			  wsar, &ac, 1, geom->zbe2, geom->zbe1);
      }
    }
  }

  /* Add sediment locations beneath bottom auxiliary cells */
  for (cs = 1; cs <= geom->sgnumS; cs++) {
    cc = cs;                    /* cs = global surface location */
    /* Search down for the sediment : cc = global sediment location */
    while (cc != geom->zm1[cc]) {
      cc = geom->zm1[cc];
    }
    c = geom->zp1[cc];          /* c = global bottom location */
    /* Used to have && geom->fm[c].ac==0 in the statement below */
    /* but this never created auxiliary sediment cells since the */
    /* cell above the sediment is always assigned an auxilliary.  */
    if (geom->fm[c].wn && geom->fm[c].wn != wn && geom->fm[cs].ac) {
      cu = cc;
      /* If a step in the bathymetry occurs at the wet-auxiliary */
      /* face then the auxiliary sediment and non-sediment cell may */
      /* be isolated (non-contiguous maps). Therefore need to loop */
      /* up from the sediment cell until the non-sediment auxiliary */
      /* cell is located. This is necessary if cell thickness arrays */
      /* (dz[]) are set at auxiliary locations since the bottom of */
      /* the water column is located where the downward map is self- */
      /* mapping in this case (i.e. the auxiliary sediment cell must */
      /* be the only cell which is self mapping).  */
      while (geom->fm[cu].ac == 0) {
        ac++;
        geom->fm[cu].ac = ac;
        wsar[ac] = cu;
        cu = geom->zp1[cu];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the local wet array and local maps */
  window->enon = ac;
  alloc_geom(window, (MAP_A | GRID_A | WINDOW_A));

  /*-----------------------------------------------------------------*/
  /* Fill the window structure with sparse locations (global for wsa */
  /* and local for the maps) from the buffer arrays.  */
  if (ac > geom->sgnum)
    hd_quit("get_local_maps: Error at %d %d\n", wn, cc);

  for (cc = 1; cc <= ac; cc++) {
    window->wsa[cc] = wsar[cc];
    window->xp1[cc] = nxp1[cc];
    window->xm1[cc] = nxm1[cc];
    window->yp1[cc] = nyp1[cc];
    window->ym1[cc] = nym1[cc];
    window->zp1[cc] = nzp1[cc];
    window->zm1[cc] = nzm1[cc];
    window->xmyp1[cc] = nxmyp1[cc];
    window->xpym1[cc] = nxpym1[cc];

    /*
       window->xp2[cc]=nxp2[cc]; window->xm2[cc]=nxm2[cc];
       window->yp2[cc]=nyp2[cc]; window->ym2[cc]=nym2[cc];
       window->zp2[cc]=nzp2[cc]; window->zm2[cc]=nzm2[cc]; */
    if (window->wsa[cc] <= 0 || window->wsa[cc] > geom->sgnum)
      hd_quit
        ("get_local_maps: Error at window %d local location %d = %d\n", wn,
         cc, window->wsa[cc]);
  }

  /*-----------------------------------------------------------------*/
  /* Get the vertical maps for all auxiliary cells. These local */
  /* locations are derived through the global vertical maps.  */
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];

    /* Auxiliary cells and ghost cells */
    if ((geom->fm[c].wn != wn && geom->fm[geom->zp1[c]].wn != wn) ||
        (geom->fm[c].wn == 0 && geom->fm[geom->zp1[c]].wn == wn)) {
      cl = geom->fm[c].ac;
      if (cl) {
        /* Upward maps for auxiliary cells and ghost cells */
        cu = geom->fm[geom->zp1[c]].ac;
        /* Upward maps for sub-surface ghost cells (i.e. ghost */
        /* cells where the cell above is wet).  */
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

  /* Get the undefined lateral maps for 2D auxiliary cells */
  for (cc = window->snonS; cc <= window->enonS; cc++)
    aux_maps(window, cc);
  /* Get the undefined lateral maps for 3D auxiliary cells */
  for (cc = window->snon; cc <= window->enon; cc++)
    aux_maps(window, cc);

  /*-----------------------------------------------------------------*/
  /* Make the global to local map for this window. This differs from */
  /* the global to local map in geom in that it is defined only over */
  /* all global cells in the window, including ghost cells which are */
  /* associated with zero window and local coordinate in geom->fm.  */
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
  }

  /*-----------------------------------------------------------------*/
  /* Reset the window maps for zoomed grids over inner explicit maps */
  if ((geom->zmfe1 || geom->zmfe2) && geom->emmask) {
    int *mask = geom->emmask;
    int zc, cd, mode = 0;
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      if (window->zoomf > 1 && mask[c]) {
	cs = c;
	cd = mask[c];
	if (geom->ym1[cs] == cd) mode = 1;
	for (zc = 1; zc <= window->zmee1; zc++) {
	  if (!mode) {
	    cs = geom->xp1[cs];
	    cd = geom->xp1[cd];
	  }
	}
	for (zc = 1; zc <= window->zmee2; zc++) {
	  if (mode) {
	    cs = geom->yp1[cs];
	    cd = geom->yp1[cd];
	  }
	}
	cs = geom->fm[cs].sc;
	cd = geom->fm[cd].sc;
	if (mode) {
	  window->ym1[cs] = cd;
	  window->xpym1[cs] = window->xp1[cd];
	} else {
	  window->xm1[cs] = cd;
	  window->xmyp1[cs] = window->yp1[cd];
	}
	cs = cd = c;
	for (zc = 1; zc <= window->zmfe1; zc++) {
	  if (!mode) {
	    cd = geom->xp1[cd];
	  }
	}
	for (zc = 1; zc <= window->zmfe2; zc++) {
	  if (mode) {
	    cd = geom->yp1[cd];
	  }
	}
	cs = window->fm[cs].sc;
	cd = window->fm[cd].sc;
	if (mode) {
	  window->yp1[cs] = cd;
	} else {
	  window->xp1[cs] = cd;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Deallocate the buffer arrays */
  i_free_1d(wsar);
  i_free_1d(nxp1);
  i_free_1d(nxm1);
  i_free_1d(nyp1);
  i_free_1d(nym1);
  i_free_1d(nzp1);
  i_free_1d(nzm1);
  i_free_1d(nxmyp1);
  i_free_1d(nxpym1);
  i_free_1d(nxmym1);
  i_free_1d(nxpyp1);
  i_free_1d(nypxm1);
  i_free_1d(xmyp1);
  i_free_1d(xpym1);
  i_free_1d(xmym1);
  i_free_1d(xpyp1);
  i_free_1d(ypxm1);
  i_free_1d(one);
  /*
     i_free_1d(nxp2); i_free_1d(nxm2); i_free_1d(nyp2); i_free_1d(nym2);
     i_free_1d(nzp2); i_free_1d(nzm2); */
}

/* END get_local_maps()                                              */
/*-------------------------------------------------------------------*/


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
/* Original Pre-V1652 mapping                                        */
void local_map_build_o(int c,     /* Global sparse coordinate */
                     int cc,    /* Local sparse coordinate */
                     int wn,    /* Window number */
                     int *map,  /* Global spatial map */
                     int *nmap, /* Local spatial map */
                     int *rmap, /* Reverse local spatial map */
                     int *wsa,  /* Local 3D wet cell sparse array */
                     int *ac,   /* Local sparse location of auxiliary cell
                                 */
                     int nsaux, /* Number of cells to include laterally */
                     int zoomf  /* Zoom factor for coarse grid mappings */
  )
{
  int n;                        /* Counter */
  int cgm;                      /* Global locations of mapped coordinate */
  int cax;                      /* Local location of mapped coordinate */
  int zn;                       /* Counter for zoom mappings */

  /* Get the global map of the global sparse coordinate */
  cgm = map[c];
  for (zn = 1; zn < zoomf; zn++)
    cgm = map[cgm];

  /* Get the local map of the local sparse coordinate */
  nmap[cc] = cax = reset_map(cgm, wn, wsa, ac);
  if (cc != cax)
    rmap[cax] = cc;

  /* If the local map is an auxiliary cell, then map this cell a */
  /* further nsaux times in the direction of map[].  */
  if (cax == *ac || geom->fm[cgm].ac != 0) {
    for (n = 1; n < nsaux; n++) {
      if (cgm != map[cgm]) {    /* Check for self-mapping locations */
        for (zn = 1; zn <= zoomf; zn++)
          cgm = map[cgm];
        nmap[cax] = reset_map(cgm, wn, wsa, ac);
        if (cax != nmap[cax])
          rmap[nmap[cax]] = cax;
        cax = nmap[cax];
      } else
        return;
    }
  }
}

/* Mappings for anisotropic grid refinement                          */
void local_map_build(int c,     /* Global sparse coordinate */
                     int cc,    /* Local sparse coordinate */
                     int wn,    /* Window number */
                     int *map,  /* Global spatial map */
                     int *nmap, /* Local spatial map */
                     int *rmap, /* Reverse local spatial map */
                     int *wsa,  /* Local 3D wet cell sparse array */
                     int *ac,   /* Local sparse location of auxiliary cell
                                 */
                     int nsaux, /* Number of cells to include laterally */
                     int *zoomf /* Zoom factor for coarse grid mappings */
  )
{
  int n;                        /* Counter */
  int cgm;                      /* Global locations of mapped coordinate */
  int cax;                      /* Local location of mapped coordinate */
  int zn;                       /* Counter for zoom mappings */

  /* Get the global map of the global sparse coordinate */
  cgm = map[c];
  for (zn = 1; zn < zoomf[c]; zn++)
    cgm = map[cgm];

  /* Get the local map of the local sparse coordinate */
  nmap[cc] = cax = reset_map(cgm, wn, wsa, ac);
  if (cc != cax)
    rmap[cax] = cc;

  /* If the local map is an auxiliary cell, then map this cell a */
  /* further nsaux times in the direction of map[].  */
  if (cax == *ac || geom->fm[cgm].ac != 0) {
    for (n = 1; n < nsaux; n++) {
      if (cgm != map[cgm]) {    /* Check for self-mapping locations */
        for (zn = 1; zn <= zoomf[c]; zn++)
          cgm = map[cgm];
        nmap[cax] = reset_map(cgm, wn, wsa, ac);
        if (cax != nmap[cax])
          rmap[nmap[cax]] = cax;
        cax = nmap[cax];
      } else
        return;
    }
  }
}

void local_map_build_d(int c,     /* Global sparse coordinate */
		       int cc,    /* Local sparse coordinate */
		       int wn,    /* Window number */
		       int *mape1, /* Global spatial map */
		       int *mape2, /* Global spatial map */
		       int *nmap, /* Local spatial map */
		       int *rmap, /* Reverse local spatial map */
		       int *wsa,  /* Local 3D wet cell sparse array */
		       int *ac,   /* Local sparse location of auxiliary cell
				   */
		       int nsaux, /* Number of cells to include laterally */
		       int *zmfe1, /* Zoom factor for coarse grid mappings */
		       int *zmfe2  /* Zoom factor for coarse grid mappings */
		       )
{
  int n;                        /* Counter */
  int cgm;                      /* Global locations of mapped coordinate */
  int cax;                      /* Local location of mapped coordinate */
  int zn;                       /* Counter for zoom mappings */

  /* Get the global map of the global sparse coordinate */
  cgm = c;
  for (zn = 1; zn <= zmfe2[c]; zn++)
    cgm = mape2[cgm];
  for (zn = 1; zn <= zmfe1[c]; zn++)
    cgm = mape1[cgm];

  /* Get the local map of the local sparse coordinate */
  nmap[cc] = cax = reset_map(cgm, wn, wsa, ac);
  if (cc != cax)
    rmap[cax] = cc;

  /* If the local map is an auxiliary cell, then map this cell a */
  /* further nsaux times in the direction of map[].  */
  if (cax == *ac || geom->fm[cgm].ac != 0) {
    for (n = 1; n < nsaux; n++) {
      if (cgm != mape1[mape2[cgm]]) {    /* Check for self-mapping locations */
        for (zn = 1; zn <= zmfe2[c]; zn++)
          cgm = mape2[cgm];
        for (zn = 1; zn <= zmfe1[c]; zn++)
          cgm = mape1[cgm];
        nmap[cax] = reset_map(cgm, wn, wsa, ac);
        if (cax != nmap[cax])
          rmap[nmap[cax]] = cax;
        cax = nmap[cax];
      } else
        return;
    }
  }
}

/* END local_map_build()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to check if a given global sparse coordinate is a member  */
/* of a given window wn. If so, then return the local sparse         */
/* coordinate in that window. If not, then add a new coordinate to   */
/* the window, set the reverse map of the new local coordinate to    */
/* the global coordinate and return the new local coordinate.        */
/*-------------------------------------------------------------------*/
int reset_map(int cm,           /* Mapped global sparse coordinate */
              int wn,           /* Window number */
              int *wsa,         /* Local 3D wet cell sparse array */
              int *al           /* Local sparse location of auxiliary cell
                                 */
  )
{
  int sc;                       /* Local sparse coordinate of window */
  int wl;                       /* Window number of c */

  sc = geom->fm[cm].sc;
  wl = geom->fm[cm].wn;

  if (wl != wn) {
    /* Add auxiliary cells to the 3D array */
    if (geom->fm[cm].ac) {
      return (geom->fm[cm].ac);
    } else {
      *al += 1;
      wsa[*al] = cm;            /* Include the auxiliary cell in wsa[] */
      geom->fm[cm].ac = *al;    /* Map the global cell to auxiliary cell */
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
void get_local_wsa(int *vec,            /* Global work sparse array  */
                   int nvec,            /* Size of vec               */
                   int nvec2D,          /* Last 2D cell loc in vec   */
                   geometry_t **window, /* Window data structure     */
                   int nwindows,        /* Number of windows         */
		   int compatible,      /* Compatibility flag        */
                   int mode /* mode=0 : cell center cells to process */
  )                         /* mode=1 : e1 face cells to process     */
                            /* mode=2 : e2 face cells to process     */
{
  int *wntmp;                   /* Temporary buffer */
  int *wntmp2D;                 /* Temporary buffer */
  int c, cc, n, nn, wn;         /* Counters */
  int ac, cgm;                  /* Sparse coordinates */
  short **mask;                 /* Set to 1 when cell is processed */
  int *wm;                      /* Window / ghost cell mask */
  int ms;                       /* Size of mask */
  int sc;                       /* Wet cell in window */
  int *gc, *gc2D;               /* Ghost cell in window */
  int *ec, *ec2D;               /* Extra auxiliary cells for mode=1 & 2 */
  int zn;
  int zc;

  /*-----------------------------------------------------------------*/
  /* Initialise the counters */
  wntmp = i_alloc_1d(nwindows + 1);
  wntmp2D = i_alloc_1d(nwindows + 1);
  gc = i_alloc_1d(nwindows + 1);
  gc2D = i_alloc_1d(nwindows + 1);
  ec = i_alloc_1d(nwindows + 1);
  ec2D = i_alloc_1d(nwindows + 1);
  wm = i_alloc_1d(geom->sgsiz);
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
    if (window[n]->enon > ms)
      ms = window[n]->enon;
  }
  mask = s_alloc_2d(ms + 1, nwindows + 1);
  for (n = 1; n <= nwindows; n++)
    for (c = 1; c <= ms; c++)
      mask[n][c] = 0;
  for (c = 1; c <= geom->sgnum; c++)
    wm[c] = geom->fm[c].wn;

  /*-----------------------------------------------------------------*/
  /* Set the ghost cells in the window mask for mode = 1 & 2 */
  zc = ZC;
  if (mode == 1) {
    /* Western edge ghost cells for e1 cells to process. These cells */
    /* are actually wet and defined to be in a window but due to the */
    /* stagger for e1 velocities at western boundaries they are */
    /* ghost cells for u1. These are identified if the west[c] is */
    /* non-zero (i.e. a wet cell in a window) and west[west[c]] is */
    /* zero (i.e. a ghost cell).  */
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];              /* Global cell to process */
      n = geom->fm[c].wn;       /* Window associated with c */
      cgm = geom->xm1[c];       /* Global cell to the west of c */
      /* The zn loop must not be done for zmode=1 */
      if (!zmode)
	for(zn = 1; zn < window[n]->zmme1; zn++)cgm = geom->xm1[cgm];
      n = geom->fm[cgm].wn;     /* Window associated with cgm */
      /* This multiple mapping for zmme1 != 0 allows ghost cells to */
      /* be identified across zoom zones.  */
      sc = 1;
      if (n)
        sc += window[n]->zmee1 + 1;
      ac = geom->xm1[cgm];      /* Global cell to the west of cgm */
      if (!zmode)
	for(zn = 1; zn <= sc; zn++)ac = geom->xm1[ac];
      if (geom->fm[cgm].wn && !(geom->fm[ac].wn))
        wm[cgm] = 0;

      /* For zoomf > 1 subsurface cells, the e1 faces must lie at    */
      /* the vertical level corresponding to the cell center, not    */
      /* the cell face. Hence, for faces which are wet but the       */
      /* center is solid, mask the face.                             */
      /* e.g. side view                                              */
      /*                  cell1   cell2                              */
      /*               |__.__.__|__.__.__|                           */
      /*               |  .  .  y  .  .  |                           */
      /*               |__.b1.__|__.__.__|                           */
      /*               |     ]  x  .  .  |                           */
      /*               |solid]__|__.b2.__|                           */
      /*               |        ]solid   |                           */
      /*                                                             */
      /* Cell1 has a step in the bathymetry at the 3rd subzoom cell. */
      /* The bottom u1 face (marked as x) corresponding to the cell  */
      /* center, b2, is wet in cell 2  but should be located at y,   */
      /* the face corresponding to the cell center b1.               */
      ac = cgm;
      for(zn = 1; zn <= window[n]->zmee1; zn++)
	ac = geom->xm1[ac];
      if (geom->fm[c].wn && !(geom->fm[ac].wn))
        wm[c] = 0;

    }
    /* Mark single wet cells with land to the east and west as ghost */
    /* cells. These cells are not captured with the loop above since */
    /* since there is no global wet cell to process to the east at   */
    /* xp1. These cells are required for the advective terms.        */
    if (!(compatible & V1283)) {
      for (cc = 1; cc < geom->nbpt; cc++) {
	c = geom->bpt[cc];
	cgm = geom->bin[cc];
	if (c == geom->xm1[cgm])
	  if (geom->fm[cgm].wn && !(geom->fm[c].wn))
	    if (wm[cgm]) wm[cgm] = 0;
      }
    }
    zc = ZE1;
  }
  if (mode == 2) {
    /* Southern edge ghost cells for e2 cells to process.  */
    for (cc = 1; cc <= nvec; cc++) {
      c = vec[cc];              /* Global cell to process */
      n = geom->fm[c].wn;       /* Window associated with c */
      cgm = geom->ym1[c];       /* Global cell to the south of c */
      if (!zmode)
	for(zn = 1; zn < window[n]->zmme2; zn++)cgm = geom->ym1[cgm];
      n = geom->fm[cgm].wn;     /* Window associated with cgm */
      sc = 1;
      if (n)
        sc += window[n]->zmee2;
      ac = geom->ym1[cgm];      /* Global cell to the south of cgm */
      if (!zmode)
	for(zn = 1; zn <= sc; zn++)ac = geom->ym1[ac];
      if (geom->fm[cgm].wn && !(geom->fm[ac].wn))
        wm[cgm] = 0;
      /* Mask cell faces whose centers at the same vertical level    */
      /* are solid.                                                  */
      ac = cgm;
      for(zn = 1; zn <= window[n]->zmee2; zn++)
	ac = geom->ym1[ac];
      if (geom->fm[c].wn && !(geom->fm[ac].wn))
        wm[c] = 0;

    }
    if (!(compatible & V1283)) {
      for (cc = 1; cc < geom->nbpt; cc++) {
	c = geom->bpt[cc];
	cgm = geom->bin[cc];
	if (c == geom->ym1[cgm])
	  if (geom->fm[cgm].wn && !(geom->fm[c].wn))
	    if (wm[cgm]) wm[cgm] = 0;
      }
    }
    zc = ZE2;
  }

  /*-----------------------------------------------------------------*/
  /* Count the number of cells to process in each window. This       */
  /* includes the wet cells in each window and the first lateral     */
  /* auxiliary cell neighboring a wet cell. These extra cells are    */
  /* required to calculate the horizontal fluxes through faces (for  */
  /* tracers) or centers (velocity) which are subsequently used to   */
  /* update the wet cell values.                                     */
  /* Note for zoomed cells the zoom codes at locations corresponding */
  /* to normal OBC cells in zoomed grids are turned off in           */
  /* set_zoom_OBC() so they are not included as wet cells here.      */
  /* Count wet cells in window n.                                    */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    sc = geom->fm[c].sc;
    if (!(geom->zoomc[geom->m2d[c]] & zc))
      continue;
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
  /* Note that zoomed grid OBC cells have been set to correspond to  */
  /* cell centers in the zoomed grid in OBC_build() so that they are */
  /* naturally included in the OBC vectors here.                     */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      if (mode == 0)
        cgm = window[n]->open[nn]->no3_t;
      else if (mode == 1)
        cgm = window[n]->open[nn]->no3_e1;
      else
        cgm = window[n]->open[nn]->no3_e2;
      for (cc = 1; cc <= cgm; cc++) {
	if (mode == 0)
          ac = window[n]->open[nn]->obc_t[cc];
        else if (mode == 1)
          ac = window[n]->open[nn]->obc_e1[cc];
        else
          ac = window[n]->open[nn]->obc_e2[cc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          if (ac <= window[n]->enonS)
            wntmp2D[n]++;
          mask[n][ac] = 1;
        }
      }
      /* Add the yp1/xp1 at the last cell for U1BDRY/U2BDRY into the  */
      /* auxiliary cells. This is required in bcond_ele = NOTHIN.     */
      if (mode == 0) {
	window[n]->open[nn]->no2_a = window[n]->open[nn]->no3_a = 0;
	if (compatible & V1957) {
	  cc = window[n]->open[nn]->no2_t;
	  c = window[n]->open[nn]->obc_t[cc];
	  if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	  if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
	  cgm = window[n]->wsa[ac];
	  if (wm[cgm] > 0 && wm[cgm] != n && mask[n][ac] == 0) {
	    wntmp[n]++;
	    if (ac <= window[n]->enonS)
	      wntmp2D[n]++;
	    mask[n][ac] = 1;
	  }
	} else {
	  for(cc = 1; cc <= window[n]->open[nn]->no3_t; cc++) {
	    c = window[n]->open[nn]->obc_t[cc];
	    if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	    if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
	    cgm = window[n]->wsa[ac];
	    if (wm[cgm] > 0 && wm[cgm] != n && mask[n][ac] == 0) {
	      wntmp[n]++;
	      window[n]->open[nn]->no3_a++;
	      if (ac <= window[n]->enonS) {
		wntmp2D[n]++;
		window[n]->open[nn]->no2_a++;
	      }
	      mask[n][ac] = 1;
	    }
	  }
	}
	
	/* Add the e1/e2 open boundary locations associated with */
	/* R_EDGE/F_EDGE boundaries to the auxiliary cells. */
	if (window[n]->open[nn]->ocodex & R_EDGE) {
	  for(cc = 1; cc <= window[n]->open[nn]->no3_e1; cc++) {
	    ac = window[n]->open[nn]->obc_e1[cc];
            if (mask[n][ac] == 0) {
              wntmp[n]++;
              if (ac <= window[n]->enonS)
                wntmp2D[n]++;
              mask[n][ac] = 1;
            }
          }
	}
	if (window[n]->open[nn]->ocodey & F_EDGE) {
	  for(cc = 1; cc <= window[n]->open[nn]->no3_e2; cc++) {
	    ac = window[n]->open[nn]->obc_e2[cc];
            if (mask[n][ac] == 0) {
              wntmp[n]++;
              if (ac <= window[n]->enonS)
                wntmp2D[n]++;
              mask[n][ac] = 1;
            }
          }
	}
      }
    }
  }

  /*------------------------------------------------------------------*/
  /* Count the first auxiliary cell in window n */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];                /* Global cell to process */
    n = wm[c];                  /* Window associated with cell c */
    sc = geom->fm[c].sc;        /* Local cell in window n corresponding to
                                   c */

    /* Note : the test !n is only valid for zoomed windows. The cells */
    /* to process vector, vec[], uses only the wet cells (1 to nvec) */
    /* and these should never contain ghost cells (i.e. cells where */
    /* n=0). For zoomed cells it is possible that the zoomed cell */
    /* center is a member of vec[] and hence not a ghost cell but the */
    /* cell center was set to a ghost cell above in wm[]. Thus n=0 */
    /* which is not allowed in the tests below.  */
    if (!(geom->zoomc[geom->m2d[c]] & zc) || !n)
      continue;

    /* Eastern edge auxiliary and ghost cells. These cells are */
    /* required for the specification of fluxes used to subsequently */
    /* update tracers, e1 and e2 velocities.  */
    cgm = geom->xp1[c];
    for (zn = 1; zn < window[n]->zmme1; zn++)
      cgm = geom->xp1[cgm];      
    if (wm[cgm] != n) {
      ac = window[n]->xp1[sc];
      if (mask[n][ac] == 0) {
        wntmp[n]++;
        if (wm[cgm] == 0)
          gc[n]++;
        if (cc <= nvec2D) {
          wntmp2D[n]++;
          if (wm[cgm] == 0)
            gc2D[n]++;
        }
        mask[n][ac] = 1;
      }
    }

    if (mode == 1) {
      /* Western edge auxiliary and ghost cells. These cells are */
      /* required for the u1 non-linear terms.  */
      cgm = geom->xm1[c];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      if (wm[cgm] != n) {
        ac = window[n]->xm1[sc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          if (wm[cgm] == 0)
            gc[n]++;
          if (cc <= nvec2D) {
            wntmp2D[n]++;
            if (wm[cgm] == 0)
              gc2D[n]++;
          }
          mask[n][ac] = 1;
        }
      }
    }

    /* Northern edge auxiliary and ghost cells */
    cgm = geom->yp1[c];
    for (zn = 1; zn < window[n]->zmme2; zn++)
      cgm = geom->yp1[cgm];
    if (wm[cgm] != n) {
      ac = window[n]->yp1[sc];
      if (mask[n][ac] == 0) {
        wntmp[n]++;
        if (wm[cgm] == 0)
          gc[n]++;
        if (cc <= nvec2D) {
          wntmp2D[n]++;
          if (wm[cgm] == 0)
            gc2D[n]++;
        }
        mask[n][ac] = 1;
      }
    }

    if (mode == 2) {
      /* Southern edge auxiliary and ghost cells. These cells are */
      /* required for the u2 non-linear terms.  */
      cgm = geom->ym1[c];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      if (wm[cgm] != n) {
        ac = window[n]->ym1[sc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          if (wm[cgm] == 0)
            gc[n]++;
          if (cc <= nvec2D) {
            wntmp2D[n]++;
            if (wm[cgm] == 0)
              gc2D[n]++;
          }
          mask[n][ac] = 1;
        }
      }
    }
  }

  /* Get the extra auxiliary cells for u1 and u2 velocity (x3_e1 and */
  /* x3_e2). This must be performed after the main auxiliary cells */
  /* (a3_e1 and a3_e2) are defined above so that where cells may be */
  /* defined as both a3_ and x3_ (e.g. at southern inside corners for */
  /* u1 and western inside corners for u2) the cell is included in */
  /* a3_ vector first.  */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];                /* Global cell to process */
    n = wm[c];                  /* Window associated with cell c */
    sc = geom->fm[c].sc;        /* Local cell in window n corresponding to
                                   c */
    if (!(geom->zoomc[geom->m2d[c]] & zc) || !n)
      continue;

    /* Western edge auxiliary cells for u2 velocity. Cell thickness */
    /* is required at these locations for use in the advective flux */
    /* calculation.  */
    if (mode == 2) {
      cgm = geom->xm1[c];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      if (wm[cgm] != n) {
        ac = window[n]->xm1[sc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          ec[n]++;
          if (cc <= nvec2D) {
            wntmp2D[n]++;
            ec2D[n]++;
          }
          mask[n][ac] = 1;
        }
      }
    }

    /* Southern edge auxiliary cells for u1 velocity.  */
    if (mode == 1) {
      cgm = geom->ym1[c];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      if (wm[cgm] != n) {
        ac = window[n]->ym1[sc];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          ec[n]++;
          if (cc <= nvec2D) {
            wntmp2D[n]++;
            ec2D[n]++;
          }
          mask[n][ac] = 1;
        }
      }

      /* The bottom left diagonal (xmym1) is required in the */
      /* calculation of the stress tensors for horizontal momentum */
      /* diffusion. The _e1 cells to process vectors are used to */
      /* assign the tensors hence the xmym1 cell is included here.  */
      /* This cell may be missed using xm1 or ym1 maps for e1 faces */
      /* in the situation where land lies at xm1 and another window */
      /* is located at ym1.  */
      cgm = geom->xm1[geom->ym1[c]];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      if (wm[cgm] != n) {
        ac = window[n]->xm1[window[n]->ym1[sc]];
        if (mask[n][ac] == 0) {
          wntmp[n]++;
          ec[n]++;
          if (cc <= nvec2D) {
            wntmp2D[n]++;
            ec2D[n]++;
          }
          mask[n][ac] = 1;
        }
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the cells to process vectors */
  for (n = 1; n <= nwindows; n++) {
    gc[n] = wntmp[n] - (gc[n] + ec[n]);
    gc2D[n] = wntmp2D[n] - (gc2D[n] + ec2D[n]);
    ec[n] += gc[n];
    ec2D[n] += gc2D[n];

    if (mode == 0) {
      window[n]->n3_t = wntmp[n]; /* Total 3D cells to process */
      window[n]->n2_t = wntmp2D[n]; /* Total 2D cells to process */
      window[n]->a3_t = gc[n];  /* 3D wet + auxiliary cells */
      window[n]->a2_t = gc2D[n];  /* 2D wet + auxiliary cells */
      window[n]->w3_t = i_alloc_1d(wntmp[n] + 1);
      window[n]->w2_t = i_alloc_1d(wntmp2D[n] + 1);
    }
    if (mode == 1) {
      window[n]->n3_e1 = wntmp[n];
      window[n]->n2_e1 = wntmp2D[n];
      window[n]->a3_e1 = gc[n];
      window[n]->a2_e1 = gc2D[n];
      window[n]->x3_e1 = ec[n];
      window[n]->x2_e1 = ec2D[n];
      window[n]->w3_e1 = i_alloc_1d(wntmp[n] + 1);
      window[n]->w2_e1 = i_alloc_1d(wntmp2D[n] + 1);
      /*
         memset(window[n]->w3_e1,0,sizeof(int)*(wntmp[n]+1));
         memset(window[n]->w2_e1,0,sizeof(int)*(wntmp2D[n]+1)); */
    }
    if (mode == 2) {
      window[n]->n3_e2 = wntmp[n];
      window[n]->a3_e2 = gc[n];
      window[n]->a2_e2 = gc2D[n];
      window[n]->n2_e2 = wntmp2D[n];
      window[n]->x3_e2 = ec[n];
      window[n]->x2_e2 = ec2D[n];
      window[n]->w3_e2 = i_alloc_1d(wntmp[n] + 1);
      window[n]->w2_e2 = i_alloc_1d(wntmp2D[n] + 1);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Fill the local window cells to process arrays */
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

  /* First fill the cells to process vectors with wet cells */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    if (n > 0 && geom->zoomc[geom->m2d[c]] & zc)
      wsa_cells(window[n], &wntmp[n], &wntmp2D[n], geom->fm[c].sc, mask[n],
                mode);
  }

  /* Save the number of wet cells in each vector */
  for (n = 1; n <= nwindows; n++) {
    if (mode == 0) {
      window[n]->v3_t = wntmp[n] - 1;
      window[n]->v2_t = wntmp2D[n] - 1;
    }
    if (mode == 1) {
      window[n]->v3_e1 = wntmp[n] - 1;
      window[n]->v2_e1 = wntmp2D[n] - 1;
    }
    if (mode == 2) {
      window[n]->v3_e2 = wntmp[n] - 1;
      window[n]->v2_e2 = wntmp2D[n] - 1;
    }
  }

  /* Fill the cells to process vectors with OBC cells */
  for (n = 1; n <= nwindows; n++) {
    for (nn = 0; nn < window[n]->nobc; nn++) {
      if (mode == 0)
        cgm = window[n]->open[nn]->no3_t;
      else if (mode == 1)
        cgm = window[n]->open[nn]->no3_e1;
      else
        cgm = window[n]->open[nn]->no3_e2;
      for (cc = 1; cc <= cgm; cc++) {
        if (mode == 0)
          ac = window[n]->open[nn]->obc_t[cc];
        else if (mode == 1)
          ac = window[n]->open[nn]->obc_e1[cc];
        else
          ac = window[n]->open[nn]->obc_e2[cc];
        wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
      }
    }
  }

  /* Save the number of OBC cells in each vector */
  for (n = 1; n <= nwindows; n++) {
    if (mode == 0) {
      window[n]->b3_t = wntmp[n] - 1;
      window[n]->b2_t = wntmp2D[n] - 1;
    }
    if (mode == 1) {
      window[n]->b3_e1 = wntmp[n] - 1;
      window[n]->b2_e1 = wntmp2D[n] - 1;
    }
    if (mode == 2) {
      window[n]->b3_e2 = wntmp[n] - 1;
      window[n]->b2_e2 = wntmp2D[n] - 1;
    }
  }

  /* Add the yp1/xp1 at the last cell for U1BDRY/U2BDRY into the      */
  /* auxiliary cells.                                                 */
  if (mode == 0) {
    for (n = 1; n <= nwindows; n++) {
      for (nn = 0; nn < window[n]->nobc; nn++) {
	if (window[n]->open[nn]->no3_a)
	  window[n]->open[nn]->obc_a = i_alloc_1d(window[n]->open[nn]->no3_a + 1);
	if (compatible & V1957) {
	  cc = window[n]->open[nn]->no2_t;
	  c = window[n]->open[nn]->obc_t[cc];
	  if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	  if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
	  cgm = window[n]->wsa[ac];
	  if (wm[cgm] > 0 && wm[cgm] != n)
	    wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
	} else {
	  window[n]->open[nn]->no3_a = 1;
	  for(cc = 1; cc <= window[n]->open[nn]->no3_t; cc++) {
	    c = window[n]->open[nn]->obc_t[cc];
	    if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	    if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
	    cgm = window[n]->wsa[ac];
	    if (wm[cgm] > 0 && wm[cgm] != n) {
	      wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
	      window[n]->open[nn]->obc_a[window[n]->open[nn]->no3_a] = ac;
	      window[n]->open[nn]->no3_a++;
	    }
	  }
	}

	/* Add the e1/e2 open boundary locations associated with */
	/* R_EDGE/F_EDGE boundaries to the auxiliary cells. */
	if (window[n]->open[nn]->ocodex & R_EDGE) {
	  for(cc = 1; cc <= window[n]->open[nn]->no3_e1; cc++) {
	    ac = window[n]->open[nn]->obc_e1[cc];
	    wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
	  }
	}
	if (window[n]->open[nn]->ocodey & F_EDGE) {
	  for(cc = 1; cc <= window[n]->open[nn]->no3_e2; cc++) {
	    ac = window[n]->open[nn]->obc_e2[cc];
	    wsa_cells(window[n], &wntmp[n], &wntmp2D[n], ac, mask[n], mode);
	  }
	}
      }
    }
  }

  /* Fill the cells to process vectors with the first auxiliary cell */
  /* neighbouring the wet cells.  */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    ac = geom->fm[c].sc;
    if (!(geom->zoomc[geom->m2d[c]] & zc) || !n)
      continue;

    /* Eastern edge auxiliary and ghost cells */
    cgm = geom->xp1[c];
    for (zn = 1; zn < window[n]->zmme1; zn++)
      cgm = geom->xp1[cgm];
    wn = wm[cgm];
    if (wn != n) {
      if (wn == 0)
        wsa_cells(window[n], &gc[n], &gc2D[n], window[n]->xp1[ac], mask[n],
                  mode);
      else
        wsa_cells(window[n], &wntmp[n], &wntmp2D[n], window[n]->xp1[ac],
                  mask[n], mode);
    }

    if (mode == 1) {
      /* Western edge auxiliary/ghost cells for e1 cells to process */
      cgm = geom->xm1[c];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      wn = wm[cgm];
      if (wn != n) {
        if (wn == 0)
          wsa_cells(window[n], &gc[n], &gc2D[n], window[n]->xm1[ac],
                    mask[n], mode);
        else
          wsa_cells(window[n], &wntmp[n], &wntmp2D[n], window[n]->xm1[ac],
                    mask[n], mode);
      }
    }

    /* Northern edge auxiliary and ghost cells */
    cgm = geom->yp1[c];
    for (zn = 1; zn < window[n]->zmme2; zn++)
      cgm = geom->yp1[cgm];
    wn = wm[cgm];
    if (wn != n) {
      if (wn == 0)
        wsa_cells(window[n], &gc[n], &gc2D[n], window[n]->yp1[ac], mask[n],
                  mode);
      else
        wsa_cells(window[n], &wntmp[n], &wntmp2D[n], window[n]->yp1[ac],
                  mask[n], mode);
    }

    if (mode == 2) {
      /* Southern edge auxiliary/ghost cells for e1 cells to process */
      cgm = geom->ym1[c];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      wn = wm[cgm];
      if (wn != n) {
        if (wn == 0)
          wsa_cells(window[n], &gc[n], &gc2D[n], window[n]->ym1[ac],
                    mask[n], mode);
        else
          wsa_cells(window[n], &wntmp[n], &wntmp2D[n], window[n]->ym1[ac],
                    mask[n], mode);
      }
    }
  }

  /* Fill the cells to process vectors with extra auxiliary cells */
  /* neighbouring the wet cells. Again this is performed after the */
  /* above loop so that cells that are may be both a3_ and x3_ are */
  /* included in the a3_ cells to process first.  */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    n = wm[c];
    ac = geom->fm[c].sc;
    if (!(geom->zoomc[geom->m2d[c]] & zc) || !n)
      continue;

    if (mode == 2) {
      /* Western edge ghost cells for e2 cells to process */
      cgm = geom->xm1[c];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      wn = wm[cgm];
      if (wn != n)
        wsa_cells(window[n], &ec[n], &ec2D[n], window[n]->xm1[ac], mask[n],
                  mode);
    }

    if (mode == 1) {
      /* Southern edge ghost cells for e1 cells to process */
      cgm = geom->ym1[c];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      wn = wm[cgm];
      if (wn != n)
        wsa_cells(window[n], &ec[n], &ec2D[n], window[n]->ym1[ac], mask[n],
                  mode);

      /* Get the bottom left diagonal (xmym1) e1 cell to process */
      cgm = geom->xm1[geom->ym1[c]];
      for (zn = 1; zn < window[n]->zmme2; zn++)
        cgm = geom->ym1[cgm];
      for (zn = 1; zn < window[n]->zmme1; zn++)
        cgm = geom->xm1[cgm];
      wn = wm[cgm];
      if (wn != n)
        wsa_cells(window[n], &ec[n], &ec2D[n],
                  window[n]->xm1[window[n]->ym1[ac]], mask[n], mode);
    }
  }

  i_free_1d(wntmp);
  i_free_1d(wntmp2D);
  i_free_1d(gc);
  i_free_1d(gc2D);
  i_free_1d(ec);
  i_free_1d(ec2D);
  i_free_1d(wm);
  s_free_2d(mask);
}

/* END get_local_wsa()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Fills the cells to process vectors                                */
/*-------------------------------------------------------------------*/
void wsa_cells(geometry_t *window,  /* Window data structure */
               int *loc3,       /* Location of local 3D cell to process */
               int *loc2,       /* Location of local 2D cell to process */
               int c,           /* Local sparse cell to process */
               short *mask,     /* Mask for cells already processed */
               int mode         /* mode=0 : cell center cells to process */
  )
                          /* mode=1 : e1 face cells to process */
                          /* mode=2 : e2 face cells to process */
{

  if (mask[c])
    return;
  mask[c] = 1;

  if (mode == 0) {
    window->w3_t[*loc3] = c;
    if (c <= window->enonS) {
      window->w2_t[*loc2] = c;
    }
  }
  if (mode == 1) {
    window->w3_e1[*loc3] = c;
    if (c <= window->enonS)
      window->w2_e1[*loc2] = c;
  }
  if (mode == 2) {
    window->w3_e2[*loc3] = c;
    if (c <= window->enonS)
      window->w2_e2[*loc2] = c;
  }
  *loc3 += 1;
  if (c <= window->enonS)
    *loc2 += 1;
}

/* END wsa_cells()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Allocate the open boundary ghost cells at R_EDGE and F_EDGE if    */
/* the maps have already been read in from file.                     */
/*-------------------------------------------------------------------*/
void get_local_obc_a(int *vec,            /* Global work array       */
		     int nvec,            /* Size of vec             */
		     int nvec2D,          /* Last 2D cell loc in vec */
		     geometry_t **window, /* Window data structure   */
		     int nwindows         /* Number of windows       */
		     )
{
  int c, cc, n, nn, wn;         /* Counters */
  int ac, cgm;                  /* Sparse coordinates */
  short **mask;                 /* Set to 1 when cell is processed */
  int *wm;                      /* Window / ghost cell mask */
  int ms;                       /* Size of mask */
  int sc;                       /* Wet cell in window */
  int zc = ZC;

  /*-----------------------------------------------------------------*/
  /* Initialise the counters */
  ms = 0;
  for (n = 1; n <= nwindows; n++) {
    if (window[n]->enon > ms)
      ms = window[n]->enon;
  }
  mask = s_alloc_2d(ms + 1, nwindows + 1);
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
    if (!(geom->zoomc[geom->m2d[c]] & zc))
      continue;
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
	if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
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
	if (window[n]->open[nn]->type & U1BDRY) ac = window[n]->yp1[c];
	if (window[n]->open[nn]->type & U2BDRY) ac = window[n]->xp1[c];
	cgm = window[n]->wsa[ac];
	if (wm[cgm] > 0 && wm[cgm] != n) {
	  window[n]->open[nn]->obc_a[window[n]->open[nn]->no3_a] = ac;
	  window[n]->open[nn]->no3_a++;
	}
      }
    }
  }
  i_free_1d(wm);
  s_free_2d(mask);
}

/* END get_local_obc_a()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to decompose the global open boundarys into OBCs within   */
/* each window.                                                      */
/*-------------------------------------------------------------------*/
void OBC_build(open_bdrys_t **open, /* Global open boundary structure */
               geometry_t **window, /* Window data structure */
               int nwindows     /* Number of windows */
  )
{
  int n, nn, tn, m;             /* Counters */
  int c, cl, cc, c1, c2;        /* Sparse coordinates / counters */
  int i1, i2, cy;               /* Interior local sparse coordinates */
  int g1, g2;                   /* Interior global sparse coordinates */
  int *nobc;                    /* Number of OBC's for each window */
  short **ff;                   /* Flag for file forcing */
  int zfe1, zfe2;

  if (!geom->nobc)
    return;
  geom->owc = s_alloc_2d(nwindows + 1, geom->nobc);
  nobc = i_alloc_1d(nwindows + 1);
  ff = s_alloc_2d(geom->nobc, nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Count the number of windows containing open boundaries */
  for (n = 0; n < geom->nobc; n++)
    for (nn = 1; nn <= nwindows; nn++)
      geom->owc[n][nn] = -1;
  for (n = 0; n < geom->nobc; n++) {
    for (cc = 1; cc <= open[n]->no3_e1; cc++) {
      c = open[n]->obc_e1[cc];
      geom->owc[n][geom->fm[c].wn] = n;
    }
    for (cc = 1; cc <= open[n]->no3_e2; cc++) {
      c = open[n]->obc_e2[cc];
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
  /* Allocate memory for the boundary structures in each window */
  for (nn = 1; nn <= nwindows; nn++) {
    window[nn]->open =
      (open_bdrys_t **)malloc(sizeof(open_bdrys_t *) * nobc[nn]);
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
  /* Fill the window OBC structues with OBC information */
  for (nn = 1; nn <= nwindows; nn++) {
    c1 = 0;
    for (n = 0; n < geom->nobc; n++) {
      if (geom->owc[n][nn] != -1) {
        strcpy(window[nn]->open[c1]->name, open[n]->name);
        window[nn]->open[c1]->id = open[n]->id;
        window[nn]->open[c1]->type = open[n]->type;
        window[nn]->open[c1]->ocodex = open[n]->ocodex;
        window[nn]->open[c1]->ocodey = open[n]->ocodey;
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
	/* Tracers using a flux OBC */
	if (window[nn]->open[c1]->ntflx)
	  window[nn]->open[c1]->tflx = i_alloc_1d(window[nn]->open[c1]->ntflx);
	window[nn]->open[c1]->ntflx = 0;
        for (tn = 0; tn < open[n]->ntr; tn++) {
	  if (open[n]->bcond_tra[tn] & (TRFLUX|TRCONC|TRCONF)) {
	    window[nn]->open[c1]->tflx[window[nn]->open[c1]->ntflx] = tn;
	    window[nn]->open[c1]->ntflx++;
	  }
	}

        /* Tidal forcing for elevation */
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

        /*
           bdry_data_copy(&window[nn]->open[c1]->bdata,&open[n]->bdata);
           window[nn]->open[c1]->bdata_t=(bdry_details_t *)
           malloc(sizeof(bdry_details_t)*open[n]->ntr); for(tn=0;
           tn<open[n]->ntr; tn++) {
           dlog("init_w","%d\n",open[n]->bdata_t[tn].explct);
           bdry_data_copy(&window[nn]->open[c1]->bdata_t[tn],
           &open[n]->bdata_t[tn]); } */

        /* Point the interior cell maps on the boundaries to the */
        /* correct spatial map.  */
        if (open[n]->ocodex & (L_EDGE|LO_EDGE)) {
          window[nn]->open[c1]->nmap = window[nn]->xp1;
          window[nn]->open[c1]->omap = window[nn]->xm1;
          window[nn]->open[c1]->tmpp = window[nn]->yp1;
          window[nn]->open[c1]->tmpm = window[nn]->ym1;
        }
        if (open[n]->ocodex & R_EDGE) {
          window[nn]->open[c1]->nmap = window[nn]->xm1;
          window[nn]->open[c1]->omap = window[nn]->xp1;
          window[nn]->open[c1]->tmpp = window[nn]->yp1;
          window[nn]->open[c1]->tmpm = window[nn]->ym1;
        }
        if (open[n]->ocodey & (B_EDGE|BO_EDGE)) {
          window[nn]->open[c1]->nmap = window[nn]->yp1;
          window[nn]->open[c1]->omap = window[nn]->ym1;
          window[nn]->open[c1]->tmpp = window[nn]->xp1;
          window[nn]->open[c1]->tmpm = window[nn]->xm1;
        }
        if (open[n]->ocodey & F_EDGE) {
          window[nn]->open[c1]->nmap = window[nn]->ym1;
          window[nn]->open[c1]->omap = window[nn]->yp1;
          window[nn]->open[c1]->tmpp = window[nn]->xp1;
          window[nn]->open[c1]->tmpm = window[nn]->xm1;
        }
        geom->owc[n][nn] = c1;
        c1++;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Re-set the zoom codes at open boundary locations on the master  */
  set_zoom_OBC(geom, open, window);

  /*-----------------------------------------------------------------*/
  /* Fill the open boundary vectors with sparse locations of the */
  /* open boundaries in that window.  */
  /* First count the OBC cells in each window.  */
  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n]->no3_t = 0;
      window[nn]->open[n]->no3_e1 = 0;
      window[nn]->open[n]->no3_e2 = 0;
      window[nn]->open[n]->to3_e1 = 0;
      window[nn]->open[n]->to3_e2 = 0;
      window[nn]->open[n]->no2_t = 0;
      window[nn]->open[n]->to2_e1 = 0;
      window[nn]->open[n]->to2_e2 = 0;
    }

  for (n = 0; n < geom->nobc; n++) {
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      nn = geom->fm[c].wn;
      if (!(geom->zoomc[geom->m2d[c]] & ZCB)) continue;
      window[nn]->open[geom->owc[n][nn]]->no3_t++;
      if (geom->zp1[c] == c) {
        window[nn]->open[geom->owc[n][nn]]->no2_t++;
      }
    }

    for (cc = 1; cc <= open[n]->no3_e1; cc++) {
      c = open[n]->obc_e1[cc];
      nn = geom->fm[c].wn;
      if (!(geom->zoomc[geom->m2d[c]] & ZE1B)) continue;
      window[nn]->open[geom->owc[n][nn]]->no3_e1++;
      if (geom->zp1[c] == c)
        window[nn]->open[geom->owc[n][nn]]->no2_e1++;
    }
    for (cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
      c = open[n]->obc_e1[cc];
      nn = geom->fm[c].wn;
      if (!(geom->zoomc[geom->m2d[c]] & ZE1B)) continue;
      window[nn]->open[geom->owc[n][nn]]->to3_e1++;
      if (geom->zp1[c] == c)
        window[nn]->open[geom->owc[n][nn]]->to2_e1++;
    }
    for (cc = 1; cc <= open[n]->no3_e2; cc++) {
      c = open[n]->obc_e2[cc];
      nn = geom->fm[c].wn;
      if (!(geom->zoomc[geom->m2d[c]] & ZE2B)) continue;
      window[nn]->open[geom->owc[n][nn]]->no3_e2++;
      if (geom->zp1[c] == c) {
        window[nn]->open[geom->owc[n][nn]]->no2_e2++;
      }
    }
    for (cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
      c = open[n]->obc_e2[cc];
      nn = geom->fm[c].wn;
      if (!(geom->zoomc[geom->m2d[c]] & ZE2B)) continue;
      window[nn]->open[geom->owc[n][nn]]->to3_e2++;
      if (geom->zp1[c] == c)
        window[nn]->open[geom->owc[n][nn]]->to2_e2++;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Allocate memory for master - slave transfer vectors */
  if (nwindows == 1) {
    for (n = 0; n < geom->nobc; n++) {
      window[1]->open[n]->transfer_u1 = open[n]->transfer_u1;
      window[1]->open[n]->transfer_u2 = open[n]->transfer_u2;
      window[1]->open[n]->transfer_u1av = open[n]->transfer_u1av;
      window[1]->open[n]->transfer_u2av = open[n]->transfer_u2av;
      window[1]->open[n]->transfer_eta = open[n]->transfer_eta;
      window[1]->open[n]->t_transfer = open[n]->t_transfer;
    }
  } else {
    for (nn = 1; nn <= nwindows; nn++)
      for (n = 0; n < nobc[nn]; n++) {
	int etasize = window[nn]->open[n]->no2_t + 1;
	int u1size_n = window[nn]->open[n]->no3_e1 + 1;
	int u2size_n = window[nn]->open[n]->no3_e2 + 1;
	int u1size_t = window[nn]->open[n]->to3_e1 + 1;
	int u2size_t = window[nn]->open[n]->to3_e2 + 1;
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
	    window[nn]->open[n]->transfer_u2 = d_alloc_1d(zt * u2size_t);
	    g1 |= 4;
	  }
	  if (window[nn]->open[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u2av = d_alloc_1d(u2size_t);
	    g1 |= 4;
	  }
	  if(g1 & (1|2)) {
	    window[nn]->open[n]->tmap_u1 = i_alloc_1d(u1size_n);
	    ff[nn][n] |= 2;
	  }
	  if(g1 & 4) {
	    window[nn]->open[n]->tmap_u2 = i_alloc_1d(u2size_t);
	    ff[nn][n] |= 16;
	  }
        }

	/* u2 boundary transfers (normal and tangential).            */
	/* Note: for VELOCITY boundaries transfer_u2 and tmap_u2 may */
	/* already be allocated above. In this case free and         */
	/* re-allocate if the array size is larger (to3), or take no */
	/* action if the array size is smaller (no3).                */
        if (window[nn]->open[n]->type & U2BDRY ) {
	  g2 = 0;
	  if (window[nn]->open[n]->bcond_nor & (FILEIN | CUSTOM) &&
	      !window[nn]->open[n]->transfer_u2) {
	    window[nn]->open[n]->transfer_u2 = d_alloc_1d(zn * u2size_n);
	    g2 |= 1;
	  }
	  if (window[nn]->open[n]->bcond_nor2d & (FILEIN | CUSTOM)) {
	    window[nn]->open[n]->transfer_u2av = d_alloc_1d(u2size_n);
	    g2 |= 2;
	  }
	  if (window[nn]->open[n]->bcond_tan & (FILEIN | CUSTOM)) {
	    if (window[nn]->open[n]->transfer_u1)
	      d_free_1d(window[nn]->open[n]->transfer_u1);
	    window[nn]->open[n]->transfer_u1 = d_alloc_1d(zt * u1size_t);
	    g2 |= 4;
	  }
	  if (window[nn]->open[n]->bcond_tan2d & (FILEIN | CUSTOM)) {
	    if (window[nn]->open[n]->transfer_u1av)
	      d_free_1d(window[nn]->open[n]->transfer_u1av);
	    window[nn]->open[n]->transfer_u1av = d_alloc_1d(u1size_t);
	    g2 |= 4;
	  }
	  if(g2 & (1|2) && !window[nn]->open[n]->tmap_u2) {
	    window[nn]->open[n]->tmap_u2 = i_alloc_1d(u2size_n);
	    ff[nn][n] |= 4;
	  }
	  if(g2 & 4) {
	    if (window[nn]->open[n]->tmap_u1)
	      i_free_1d(window[nn]->open[n]->tmap_u1);
	    window[nn]->open[n]->tmap_u1 = i_alloc_1d(u1size_t);
	    ff[nn][n] |= 32;
	  }
        }

	/* Tracer transfers                                          */
        if (window[nn]->open[n]->ntt) {
          window[nn]->open[n]->t_tmap =
            i_alloc_1d(window[nn]->open[n]->no3_t + 1);
	  m = (window[nn]->open[n]->bgz) ? window[nn]->open[n]->bgz : 1;
          window[nn]->open[n]->t_transfer =
            d_alloc_2d((m * window[nn]->open[n]->no3_t) + 1, window[nn]->open[n]->ntt);
          ff[nn][n] |= 8;
        }
      }
  }
      
  /* Allocate memory for the open boundary vectors */
  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {

      /* Old values data structures */
      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->no2_t + 1);
      window[nn]->open[n]->etabold = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->etabold, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->no2_e1 + 1);
      window[nn]->open[n]->u1av_on = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u1av_on, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->to2_e1 + 1);
      window[nn]->open[n]->u1av_ot = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u1av_ot, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->no2_e2 + 1);
      window[nn]->open[n]->u2av_on = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u2av_on, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->to2_e2 + 1);
      window[nn]->open[n]->u2av_ot = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u2av_ot, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->no3_e1 + 1);
      window[nn]->open[n]->u1_on = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u1_on, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->to3_e1 + 1);
      window[nn]->open[n]->u1_ot = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u1_ot, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->no3_e2 + 1);
      window[nn]->open[n]->u2_on = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u2_on, 0, c2);

      c2 = sizeof(bdry_old_values_t) * (window[nn]->open[n]->to3_e2 + 1);
      window[nn]->open[n]->u2_ot = (bdry_old_values_t *)malloc(c2);
      memset(window[nn]->open[n]->u2_ot, 0, c2);

      window[nn]->open[n]->obc_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->oi1_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->oi2_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->cyc_t =
        i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      if (window[nn]->open[n]->ocodex & (L_EDGE|LO_EDGE|R_EDGE) || 
	  window[nn]->open[n]->ocodey & (B_EDGE|BO_EDGE|F_EDGE))
	window[nn]->open[n]->ogc_t =
	  i_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->no3_t = 1;

      c1 = window[nn]->open[n]->no3_e1 + window[nn]->open[n]->to3_e1 + 1;
      window[nn]->open[n]->obc_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->oi1_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->oi2_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->cyc_e1 = i_alloc_1d(c1);
      window[nn]->open[n]->to3_e1 = window[nn]->open[n]->no3_e1 + 1;
      window[nn]->open[n]->no3_e1 = 1;

      c1 = window[nn]->open[n]->no3_e2 + window[nn]->open[n]->to3_e2 + 1;
      window[nn]->open[n]->obc_e2 = i_alloc_1d(c1);
      window[nn]->open[n]->oi1_e2 = i_alloc_1d(c1);
      window[nn]->open[n]->oi2_e2 = i_alloc_1d(c1);
      window[nn]->open[n]->cyc_e2 = i_alloc_1d(c1);
      window[nn]->open[n]->to3_e2 = window[nn]->open[n]->no3_e2 + 1;
      window[nn]->open[n]->no3_e2 = 1;
      /*
      if (window[nn]->open[n]->ocodex & R_EDGE)
	window[nn]->open[n]->ogc_t = window[nn]->open[n]->obc_e1;
      if (window[nn]->open[n]->ocodey & F_EDGE)
	window[nn]->open[n]->ogc_t = window[nn]->open[n]->obc_e2;
      */
    }

  /*-----------------------------------------------------------------*/
  /* Fill the open boundary vectors with local sparse locations */
  for (n = 0; n < geom->nobc; n++) {
    /* Tracer open boundary vectors */
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];
      window[nn]->open[c1]->mwn = n;

      /* Re-set the cell center coordinate for zoomed grids */
      zfe1 = window[nn]->zmfe1;
      zfe2 = window[nn]->zmfe2;
      if ((zfe1 > 1 || zfe2 > 1) && !(geom->zoomc[geom->m2d[c]] & ZCB)) continue;
      c = get_zcell(geom, open[n], c, zfe1, zfe2, 0);

      if (ff[nn][c1] & 1 && cc <= open[n]->no2_t)
        window[nn]->open[c1]->tmap[window[nn]->open[c1]->no3_t] = cc;

      if (ff[nn][c1] & 8 && window[nn]->open[c1]->ntt)
        window[nn]->open[c1]->t_tmap[window[nn]->open[c1]->no3_t] = cc;

      cl = geom->fm[c].sc;
      cy = get_zcell_cy(geom, window[nn], open[n], open[n]->cyc_t[cc], 0);
      cy = geom->fm[cy].sc;
      window[nn]->open[c1]->obc_t[window[nn]->open[c1]->no3_t] = cl;
      window[nn]->open[c1]->cyc_t[window[nn]->open[c1]->no3_t] = cy;

      g1 = open[n]->oi1_t[cc];
      i1 = get_icell(window[nn], cl, open[n]->obc_t[cc], g1);
      window[nn]->open[c1]->oi1_t[window[nn]->open[c1]->no3_t] = i1;

      g2 = open[n]->oi2_t[cc];
      i2 = get_icell(window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_t[window[nn]->open[c1]->no3_t] = i2;
      if (window[nn]->open[c1]->ocodex & (L_EDGE|LO_EDGE|R_EDGE) ||
	  window[nn]->open[c1]->ocodey & (B_EDGE|BO_EDGE|F_EDGE))
	window[nn]->open[c1]->ogc_t[window[nn]->open[c1]->no3_t] = 
	  window[nn]->open[c1]->omap[cl];

      window[nn]->open[c1]->no3_t++;
    }

    /* Normal u1 velocity open boundary vectors */
    for (cc = 1; cc <= open[n]->no3_e1; cc++) {
      c = open[n]->obc_e1[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];

      /* Re-set the e1 face coordinate for zoomed grids */
      zfe1 = window[nn]->zmfe1;
      zfe2 = window[nn]->zmfe2;
      if ((zfe1 > 1 || zfe2 > 1) && !(geom->zoomc[geom->m2d[c]] & ZE1B)) continue;
      if (open[n]->ocodex & L_EDGE)
	c = get_zcell(geom, open[n], c, zfe1, zfe2, 0);

      if (ff[nn][c1] & 2)
        window[nn]->open[c1]->tmap_u1[window[nn]->open[c1]->no3_e1] = cc;

      cl = geom->fm[c].sc;
      cy = get_zcell_cy(geom, window[nn], open[n], open[n]->cyc_e1[cc], 1);
      cy = geom->fm[cy].sc;
      window[nn]->open[c1]->obc_e1[window[nn]->open[c1]->no3_e1] = cl;
      window[nn]->open[c1]->cyc_e1[window[nn]->open[c1]->no3_e1] = cy;

      g1 = open[n]->oi1_e1[cc];
      if (zfe1 > 1 && open[n]->ocodex & R_EDGE)
	i1 = get_icellz(geom, zfe1, open[n]->obc_e1[cc], g1);
      else
	i1 = get_icell(window[nn], cl, open[n]->obc_e1[cc], g1);
      window[nn]->open[c1]->oi1_e1[window[nn]->open[c1]->no3_e1] = i1;
      g2 = open[n]->oi2_e1[cc];
      i2 = get_icell(window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_e1[window[nn]->open[c1]->no3_e1] = i2;
      window[nn]->open[c1]->no3_e1++;
    }

    /* Tangential u1 velocity open boundary vectors */
    for (cc = open[n]->no3_e1 + 1; cc <= open[n]->to3_e1; cc++) {
      c = open[n]->obc_e1[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];

      /* Re-set the e1 face coordinate for zoomed grids */
      zfe1 = window[nn]->zmfe1;
      zfe2 = window[nn]->zmfe2;
      if ((zfe1 > 1 || zfe2 > 1) && !(geom->zoomc[geom->m2d[c]] & ZE1B)) continue;
      c = get_zcell(geom, open[n], c, zfe1, zfe2, 0);

      if (ff[nn][c1] & 32)
        window[nn]->open[c1]->tmap_u1[window[nn]->open[c1]->to3_e1] = cc;

      cl = geom->fm[c].sc;
      cy = get_zcell_cy(geom, window[nn], open[n], open[n]->cyc_e1[cc], 0);
      cy = geom->fm[cy].sc;
      window[nn]->open[c1]->obc_e1[window[nn]->open[c1]->to3_e1] = cl;
      window[nn]->open[c1]->cyc_e1[window[nn]->open[c1]->to3_e1] = cy;

      g1 = open[n]->oi1_e1[cc];
      i1 = get_icell(window[nn], cl, open[n]->obc_e1[cc], g1);
      window[nn]->open[c1]->oi1_e1[window[nn]->open[c1]->to3_e1] = i1;

      g2 = open[n]->oi2_e1[cc];
      i2 = get_icell(window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_e1[window[nn]->open[c1]->to3_e1] = i2;
      window[nn]->open[c1]->to3_e1++;
    }

    /* Normal u2 velocity open boundary vectors */
    for (cc = 1; cc <= open[n]->no3_e2; cc++) {
      c = open[n]->obc_e2[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];

      /* Re-set the e1 face coordinate for zoomed grids */
      zfe1 = window[nn]->zmfe1;
      zfe2 = window[nn]->zmfe2;
      if ((zfe1 > 1 || zfe2 > 1) && !(geom->zoomc[geom->m2d[c]] & ZE2B)) continue;
      if (open[n]->ocodey & B_EDGE)
	c = get_zcell(geom, open[n], c, zfe1, zfe2, 0);

      if (ff[nn][c1] & 4)
        window[nn]->open[c1]->tmap_u2[window[nn]->open[c1]->no3_e2] = cc;

      cl = geom->fm[c].sc;
      cy = get_zcell_cy(geom, window[nn], open[n], open[n]->cyc_e2[cc], 1);
      cy = geom->fm[cy].sc;
      window[nn]->open[c1]->obc_e2[window[nn]->open[c1]->no3_e2] = cl;
      window[nn]->open[c1]->cyc_e2[window[nn]->open[c1]->no3_e2] = cy;

      g1 = open[n]->oi1_e2[cc];
      if (zfe2 > 1 && open[n]->ocodey & F_EDGE)
	i1 = get_icellz(geom, zfe2, open[n]->obc_e2[cc], g1);
      else
	i1 = get_icell(window[nn], cl, open[n]->obc_e2[cc], g1);
      window[nn]->open[c1]->oi1_e2[window[nn]->open[c1]->no3_e2] = i1;

      g2 = open[n]->oi2_e2[cc];
      i2 = get_icell(window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_e2[window[nn]->open[c1]->no3_e2] = i2;
      window[nn]->open[c1]->no3_e2++;
    }

    /* Tangential u2 velocity open boundary vectors */
    for (cc = open[n]->no3_e2 + 1; cc <= open[n]->to3_e2; cc++) {
      c = open[n]->obc_e2[cc];
      nn = geom->fm[c].wn;
      c1 = geom->owc[n][nn];

      /* Re-set the e1 face coordinate for zoomed grids */
      zfe1 = window[nn]->zmfe1;
      zfe2 = window[nn]->zmfe2;
      if ((zfe1 > 1 || zfe2 > 1) && !(geom->zoomc[geom->m2d[c]] & ZE2B)) continue;
      c = get_zcell(geom, open[n], c, zfe1, zfe2, 0);

      if (ff[nn][c1] & 16)
        window[nn]->open[c1]->tmap_u2[window[nn]->open[c1]->to3_e2] = cc;

      cl = geom->fm[c].sc;
      cy = get_zcell_cy(geom, window[nn], open[n], open[n]->cyc_e2[cc], 0);
      cy = geom->fm[cy].sc;
      window[nn]->open[c1]->obc_e2[window[nn]->open[c1]->to3_e2] = cl;
      window[nn]->open[c1]->cyc_e2[window[nn]->open[c1]->to3_e2] = cy;

      g1 = open[n]->oi1_e2[cc];
      i1 = get_icell(window[nn], cl, open[n]->obc_e2[cc], g1);
      window[nn]->open[c1]->oi1_e2[window[nn]->open[c1]->to3_e2] = i1;

      g2 = open[n]->oi2_e2[cc];
      i2 = get_icell(window[nn], i1, g1, g2);
      window[nn]->open[c1]->oi2_e2[window[nn]->open[c1]->to3_e2] = i2;
      window[nn]->open[c1]->to3_e2++;
    }
  }

  for (nn = 1; nn <= nwindows; nn++)
    for (n = 0; n < nobc[nn]; n++) {
      window[nn]->open[n]->no3_t -= 1;
      window[nn]->open[n]->no3_e1 -= 1;
      window[nn]->open[n]->no3_e2 -= 1;
      window[nn]->open[n]->to2_e1 += window[nn]->open[n]->no3_e1;
      window[nn]->open[n]->to2_e2 += window[nn]->open[n]->no3_e2;
      window[nn]->open[n]->to3_e1 -= 1;
      window[nn]->open[n]->to3_e2 -= 1;

      /* Set self-mapping maps over all edges for no OBC ghost zones     */
      if (window[nn]->open[n]->bgz == 0) {
	if (window[nn]->open[n]->ocodex & L_EDGE) {
	  for (cc = 1; cc <= window[nn]->open[n]->no3_e1; cc++) {
	    c = window[nn]->open[n]->obc_e1[cc];
	    if (window[nn]->open[n]->stagger & INFACE) c = window[nn]->xm1[c];
	    window[nn]->xm1[c] = c;
	  }
	}
	if (window[nn]->open[n]->ocodey & B_EDGE) {
	  for (cc = 1; cc <= window[nn]->open[n]->no3_e2; cc++) {
	    c = window[nn]->open[n]->obc_e2[cc];
	    if (window[nn]->open[n]->stagger & INFACE) c = window[nn]->xm1[c];
	    window[nn]->ym1[c] = c;
	  }
	}
      }

      /* Get the bottom coordinate vector                                */
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
      window[nn]->open[n]->olap =
	i_alloc_1d(window[nn]->open[n]->no2_t + 1);
      window[nn]->open[n]->dum =
	d_alloc_1d(window[nn]->open[n]->no3_t + 1);
      window[nn]->open[n]->dumn =
	d_alloc_2d(window[nn]->open[n]->no3_t + 1, 
		   window[nn]->open[n]->ntr);
      window[nn]->open[n]->u1d =
	d_alloc_1d(window[nn]->open[n]->to3_e1 + 1);
      window[nn]->open[n]->u2d =
	d_alloc_1d(window[nn]->open[n]->to3_e2 + 1);
      window[nn]->open[n]->flow =
	d_alloc_1d(window[nn]->open[n]->no3_t + 1);
      /* Dummy for tracers using a flux OBC */
      if (window[nn]->open[n]->ntflx)
	window[nn]->open[n]->dumtr =
	  d_alloc_2d(window[nn]->open[n]->no3_t + 1,
		     window[nn]->open[n]->ntflx);
      /* Allocate the phase smoothing array if required */
      if (window[nn]->open[n]->spf) {
	window[nn]->open[n]->sphase =
	  d_alloc_1d(window[nn]->open[n]->no2_t + 1);
	memset(window[nn]->open[n]->sphase, 0, 
	       (window[nn]->open[n]->no2_t + 1) * sizeof(double));
      }
    }
    
  
  /*-----------------------------------------------------------------*/
  /* Print OBC information */
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
        dlog("init_w",
             "      u2 velocity = %d normal cells, %d tangential cells\n",
             window[nn]->open[n]->no3_e2,
             window[nn]->open[n]->to3_e2 - window[nn]->open[n]->no3_e2);
        dlog("init_w", "      Elevation = %d cells\n",
             window[nn]->open[n]->no2_t);
        dlog("init_w",
             "      u1av velocity = normal %d cells, %d tangential cells\n",
             window[nn]->open[n]->no2_e1,
             window[nn]->open[n]->to2_e1 - window[nn]->open[n]->no3_e1);
        dlog("init_w",
             "      u2av velocity = %d normal cells, %d tangential cells\n",
             window[nn]->open[n]->no2_e2,
             window[nn]->open[n]->to2_e2 - window[nn]->open[n]->no3_e2);
      }
    }
  }
  i_free_1d(nobc);
  s_free_2d(ff);
}

/* END OBC_build()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set zoom codes on the master recognising OBC's. The ZB */
/* code is placed on the master at the location corresponding to the */
/* window's cell center.                                             */
/*-------------------------------------------------------------------*/
void set_zoom_OBC(geometry_t *geom,    /* Global geometry */
		  open_bdrys_t **open, /* Global open boundarys */
		  geometry_t **window  /* Local geometry */
		  )
{
  int n, wn;
  int c, cz, cc;
  int zfe1, zfe2;

  for (n = 0; n < geom->nobc; n++) {

    /* Tracer open boundary vectors */
    for (cc = 1; cc <= open[n]->no2_t; cc++) {
      c = open[n]->obc_t[cc];
      wn = geom->fm[c].wn;
      zfe1 = window[wn]->zmfe1;
      zfe2 = window[wn]->zmfe2;
      cz = get_zcell(geom, open[n], c, zfe1, zfe2, 0);
      if (geom->zoomc[cz] & ZC && geom->fm[cz].wn == wn) {
	geom->zoomc[c] |= ZCB;
	if (zfe1 > 1 || zfe2 > 1)
	  geom->zoomc[cz] &= ~ZC;
      }
    }
    /* Normal u1 velocity open boundary vectors */
    for (cc = 1; cc <= open[n]->no2_e1; cc++) {
      c = cz = open[n]->obc_e1[cc];
      wn = geom->fm[c].wn;
      zfe1 = window[wn]->zmfe1;
      if (geom->zoomc[cz] & ZE1 && geom->fm[cz].wn == wn) {
	geom->zoomc[c] |= ZE1B;
	if (zfe1 > 1) {
	  geom->zoomc[cz] &= ~ZE1;
	  geom->zoomc[cz] |= ZE1NB;
	}
      }
    }
    /* Tangential u1 velocity open boundary vectors */
    for (cc = open[n]->no3_e1 + 1; cc <= open[n]->to2_e1; cc++) {
      c = open[n]->obc_e1[cc];
      wn = geom->fm[c].wn;
      zfe1 = window[wn]->zmfe1;
      zfe2 = window[wn]->zmfe2;
      cz = get_zcell(geom, open[n], c, zfe1, zfe2, 1);
      if (geom->zoomc[cz] & ZE1 && geom->fm[cz].wn == wn) {
	/* Don't include cells adjacent to solid walls */
	/*if (geom->xm1[cz] != geom->xm1[geom->ym1[cz]]) {*/
	  geom->zoomc[c] |= ZE1B;
	  if (zfe1 > 1) {
	    /*geom->zoomc[cz] &= ~ZE1;*/
	    geom->zoomc[cz] |= ZE1TB;
	  }

      }
    }
    /* Normal u2 velocity open boundary vectors */
    for (cc = 1; cc <= open[n]->no2_e2; cc++) {
      c = cz = open[n]->obc_e2[cc];
      wn = geom->fm[c].wn;
      zfe2 = window[wn]->zmfe2;
      if (geom->zoomc[cz] & ZE2 && geom->fm[cz].wn == wn) {
	geom->zoomc[c] |= ZE2B;
	if (zfe2 > 1) {
	  geom->zoomc[cz] |= ZE2NB;
	  geom->zoomc[cz] &= ~ZE2;
	}
      }
    }
    /* Tangential u2 velocity open boundary vectors */
    for (cc = open[n]->no3_e2 + 1; cc <= open[n]->to2_e2; cc++) {
      c = open[n]->obc_e2[cc];
      wn = geom->fm[c].wn;
      zfe1 = window[wn]->zmfe1;
      zfe2 = window[wn]->zmfe2;
      cz = get_zcell(geom, open[n], c, zfe1, zfe2, 2);
      if (geom->zoomc[cz] & ZE2 && geom->fm[cz].wn == wn) {
	/* Don't include cells adjacent to solid walls */
	/*if (geom->ym1[cz] != geom->ym1[geom->ym1[cz]]) {*/
	  geom->zoomc[c] |= ZE2B;
	  if (zfe2 > 1) {
	    /*geom->zoomc[cz] &= ~ZE2;*/
	    geom->zoomc[cz] |= ZE2TB;
	  }

      }
    }
  }
}

/* END set_zoom_OBC()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return an interior cell in a window to any given cell  */
/* on the boundary. The interior cell is located based on the global */
/* maps.                                                             */
/*-------------------------------------------------------------------*/
int get_icell(geometry_t *window, /* Window data structure */
              int cl,           /* Local cell on the boundary */
              int cg,           /* Global cell on the boundary */
              int ci            /* Global cell interior to cg */
  )
{
  int i1 = 0;                   /* Local interior cell */
  if (ci == geom->xp1[cg])
    i1 = window->xp1[cl];
  else if (ci == geom->xm1[cg])
    i1 = window->xm1[cl];
  else if (ci == geom->yp1[cg])
    i1 = window->yp1[cl];
  else if (ci == geom->ym1[cg])
    i1 = window->ym1[cl];
  return (i1);
}

/* END get_icell()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return an interior cell in a window to any given cell  */
/* on the boundary for zoomed R_EDGE or F_EDGE boundaries.           */
/*-------------------------------------------------------------------*/
int get_icellz(geometry_t *geom, /* Global geometry                  */
	       int zf,           /* Local cell on the boundary */
	       int cg,           /* Global cell on the boundary */
	       int ci            /* Global cell interior to cg */
  )
{
  int zc;
  if (ci == geom->xm1[cg]) {
    for (zc = 1; zc <= (int)(zf / 2); zc++)
      ci = geom->xm1[ci];
  }
  else if (ci == geom->ym1[cg]) {
    for (zc = 1; zc <= (int)(zf / 2); zc++)
      ci = geom->ym1[ci];
  }
  return (geom->fm[ci].sc);
}

/* END get_icell()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to shift a boundary cell for zoomed grids                 */
/*-------------------------------------------------------------------*/
int get_zcell(geometry_t *geom,    /* Global geometry */
	      open_bdrys_t *open,  /* Global open boundary structure */
              int cg,              /* Global cell on the boundary */
	      int zfe1,            /* e1 zoom flag for the window */
	      int zfe2,            /* e2 zoom flag for the window */
	      int mode
	      )
{
  int n;
  int c = cg;
  int zse1 = (int)(zfe1 / 2);
  int zse2 = (int)(zfe2 / 2);

  if (zse1) {
    for (n = 1; n <= zse1; n++) {
      if (open->ocodex & L_EDGE)
	c = geom->xp1[c];
      if (open->ocodex & R_EDGE)
	c = geom->xm1[c];
    }
  }
  if (zse2) {
    for (n = 1; n <= zse2; n++) {
      if (open->ocodey & B_EDGE)
	c = geom->yp1[c];
      if (open->ocodey & F_EDGE)
	c = geom->ym1[c];
    }
  }
  if (mode == 1)
    for (n = 1; n <= zse1; n++)
      c = geom->xm1[c];
  if (mode == 2)
    for (n = 1; n <= zse2; n++)
      c = geom->ym1[c];

  return(c);
}

/* END get_zcell()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to shift a boundary cyclic cell for zoomed grids          */
/*-------------------------------------------------------------------*/
int get_zcell_cy(geometry_t *geom,    /* Global geometry */
		 geometry_t *window,  /* Local geometry */
		 open_bdrys_t *open,  /* Global open boundary structure */
		 int c,               /* Global cell on the boundary */
		 int mode             /* Shift for normal OBCs */
		 )
{
  int n, lc = 0;
  int zse1 = window->zmee1;
  int zse2 = window->zmee2;

  if (zse1) {
    if (open->ocodex & L_EDGE) {
      if (mode) zse1++;
      for (n = 1; n < zse1; n++)
	c = geom->xm1[c];
      if (mode)
	lc = geom->fm[c].sc;
      else
	lc = window->xm1[geom->fm[c].sc];
    }
    if (open->ocodex & R_EDGE) {
      for (n = 1; n < zse1; n++)
	c = geom->xp1[c];
      lc = window->xp1[geom->fm[c].sc];
    }
  }
  if (zse2) {
    if (open->ocodey & B_EDGE) {
      if (mode) zse2++;
      for (n = 1; n < zse2; n++)
	c = geom->ym1[c];
      if (mode)
	lc = geom->fm[c].sc;
      else
	lc = window->ym1[geom->fm[c].sc];
    }
    if (open->ocodey & F_EDGE) {
      for (n = 1; n < zse2; n++)
	c = geom->yp1[c];
      lc = window->yp1[geom->fm[c].sc];
    }
    c = window->wsa[lc];
  }
  return(c);
}

/* END get_zcell_cy()                                                */
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
void get_local_ghost(geometry_t *geom, geometry_t *window, int n  /* Window
                                                                     number
                                                                   */
  )
{
  int c, cc;                    /* Sparse location/counter */
  int lc, ic1, ic2;             /* Local sparse coordinates */
  int nb, nc, ncS;              /* Counters */
  int zc;
  int zse1 = window->zmee1;
  int zse2 = window->zmee2;
  int dr;
  int **maps;
  int *xpyp1;
  int *xpym1;
  int *xmyp1;
  int *xmym1;
  int *dbine1, *dbine2;
  /*
     int *mask; mask=i_alloc_1d(master->enon+1);
     memset(mask,0,(master->enon+1)*sizof(int)); for(cc=1;
     cc<=master->nbpt; cc++) { mask[master->bpt[cc]]=cc; } */

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

  /* Set up the compass of directional maps */
  xpyp1 = i_alloc_1d(window->enon + 1);
  xpym1 = i_alloc_1d(window->enon + 1);
  xmyp1 = i_alloc_1d(window->enon + 1);
  xmym1 = i_alloc_1d(window->enon + 1);
  for (c = 1; c <= window->enon; c++) {
    xmyp1[c] = window->yp1[window->xm1[c]];
    xpym1[c] = window->xp1[window->ym1[c]];
    xmym1[c] = window->xm1[window->ym1[c]];
    xpyp1[c] = window->xp1[window->yp1[c]];
  }
  c = 8 * sizeof(int *);
  maps = (int **)malloc(c);
  maps[0] = window->yp1;
  maps[1] = xpyp1;
  maps[2] = window->xp1;
  maps[3] = xpym1;
  maps[4] = window->ym1;
  maps[5] = xmym1;
  maps[6] = window->xm1;
  maps[7] = xmyp1;

  /* Get the directions of the interior cells for e1 and e2 vectors */
  dbine1 = i_alloc_1d(geom->nbpt + 1);
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpte1[cc];
    lc = geom->bine1[cc];

    if (geom->xm1[lc] == c)
      dbine1[cc] = 2;
    else if (geom->xp1[lc] == c)
      dbine1[cc] = 6;
    else if (geom->ym1[lc] == c)
      dbine1[cc] = 0;
    else if (geom->yp1[lc] == c)
      dbine1[cc] = 4;
    else if (geom->xm1[geom->ym1[lc]] == c)
      dbine1[cc] = 1;
    else if (geom->xp1[geom->yp1[lc]] == c)
      dbine1[cc] = 5;
    else if (geom->xm1[geom->yp1[lc]] == c)
      dbine1[cc] = 3;
    else if (geom->xp1[geom->ym1[lc]] == c)
      dbine1[cc] = 7;
  }
  dbine2 = i_alloc_1d(geom->nbpt + 1);
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpte2[cc];
    lc = geom->bine2[cc];

    if (geom->xm1[lc] == c)
      dbine2[cc] = 2;
    else if (geom->xp1[lc] == c)
      dbine2[cc] = 6;
    else if (geom->ym1[lc] == c)
      dbine2[cc] = 0;
    else if (geom->yp1[lc] == c)
      dbine2[cc] = 4;
    else if (geom->xm1[geom->ym1[lc]] == c)
      dbine2[cc] = 1;
    else if (geom->xp1[geom->yp1[lc]] == c)
      dbine2[cc] = 5;
    else if (geom->xm1[geom->yp1[lc]] == c)
      dbine2[cc] = 3;
    else if (geom->xp1[geom->ym1[lc]] == c)
      dbine2[cc] = 7;
  }

  /*-----------------------------------------------------------------*/
  /* Get the lateral ghost vectors for cell centers.  */
  /* The lateral ghost cells that are included in a window's sparse */
  /* array comprise those cells that are east, west, north, south */
  /* and the corners south-west, north-west and south-east from a */
  /* wet cell in the window. Note that the corner maps are mapped */
  /* only once while the other maps are mapped laux times. Also note */
  /* that the map north-east doesn't exist since it is not needed.  */
  if (DEBUG("init_w"))
    dlog("init_w", "    Start cell center local ghost vectors\n");
  /* Count the local ghost cells in window n */
  window->nbpt = 0;
  window->nbptS = 0;
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];          /* Global boundary cell */
    /* Note : cannot use ANY(c,window->w3_t) here since the cells to */
    /* process vectors do not contain all ghost cells, only those */
    /* required to specify the fluxes used to update cell centers, */
    /* i.e. eastern and northern edge ghost cells.  */
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      window->nbpt++;
      if (lc <= window->enonS)
        window->nbptS++;
    }
  }

  /* Allocate memory and initialise counters */
  window->bpt = i_alloc_1d(window->nbpt + 1);
  memset(window->bpt, 0, window->nbpt * sizeof(int));
  window->bin = i_alloc_1d(window->nbpt + 1);
  memset(window->bin, 0, window->nbpt * sizeof(int));
  window->bin2 = i_alloc_1d(window->nbpt + 1);
  memset(window->bin2, 0, window->nbpt * sizeof(int));
  window->wgst = i_alloc_1d(window->enon + 1);
  window->nbpt = window->nbptS + 1;
  window->nbptS = 1;

  /* Assign the lateral ghost cells to the boundary vectors */
  for (cc = 1; cc <= geom->nbpt; cc++) {
    c = geom->bpt[cc];          /* Global boundary cell */
    dr = geom->dbin[cc];        /* Direction of bin from bpt */
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      ic1 = maps[dr][lc];       /* Interior cell to lc */
      ic2 = maps[dr][ic1];      /* Interior cell to ic1 */
      if (lc <= window->enonS) {
        window->bpt[window->nbptS] = lc;
        window->bin[window->nbptS] = ic1;
        window->bin2[window->nbptS] = ic2;
        window->nbptS++;
      } else {
        window->bpt[window->nbpt] = lc;
        window->bin[window->nbpt] = ic1;
        window->bin2[window->nbpt] = ic2;
        window->nbpt++;
      }
      window->wgst[lc] = ic1;
    }
  }
  window->nbpt--;
  window->nbptS--;

  /* Assign 2D OBC ghost cells to the window ghost array */
  if (!(geom->compatible & V1957)) {
    for (nb = 0; nb < geom->nobc; nb++) {
      open_bdrys_t *open = geom->open[nb];
      if (open->bgz) {
	int wb, m, mm;
	mm = (open->ocodex & R_EDGE || open->ocodey & F_EDGE) ? 2 : 1;
	for (cc = 1; cc <= open->no2_t; cc++) {
	  c = open->obc_t[cc];
	  for (m = 0; m < mm; m++)
	    c = open->omap[c];
	  ic1 = lc = ANYf(c, window->wsa, window->enon);
	  wb = geom->owc[nb][n];
	  if(ic1 && wb >= 0) {
	    for (m = 0; m < mm; m++)
	      ic1 = window->open[wb]->nmap[ic1];
	    window->wgst[lc] = ic1;
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the lateral ghost vectors for e1 face centered cells.  */
  /* The lateral ghost cells that are included in a window's sparse */
  /* array for e1 faces actually include wet and auxilliary cells on */
  /* western edges (southern edges for e2 faces). This means that */
  /* some auxiliary cells may be identified as ghost cells for e1 */
  /* faces at locations where a solid boundary exists at the limits */
  /* of the window. Consequently these e1 ghost cells are included */
  /* in the boundary array, but the corresponding cell center ghost */
  /* cell is not because the directional maps don't include this */
  /* ghost in sparse array. This results in different numbers of */
  /* lateral ghost cells included in the cell and face centered */
  /* boundary arrays. This only applies for window partitions and */
  /* does not apply for the global vectors or nwindows=1.            */
  /*                                                                 */
  /*                                                                 */
  /*     e.g.               |   ]   |   |   |                        */
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
  /* Count the local e1 ghost cells in window n */

  window->nbpte1 = 0;
  window->nbpte1S = 0;
  for (cc = 1; cc <= geom->nbpte1; cc++) {
    c = geom->bpte1[cc];        /* Global boundary cell */
    dr = dbine1[cc];            /* Direction of bin from bpt */
    /* The ghost cell for e1 faces on western edges is a wet cell.  */
    /* For zoomed windows this cell is identified as a ghost cell, */
    /* but needs to be shifted xp1 by zs so as to correspond to */
    /* the correct cell center.  */
    if (dr == 2)
      for (zc = 1; zc <= zse1; zc++)
        c = geom->xp1[c];
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      window->nbpte1++;
      if (lc <= window->enonS)
        window->nbpte1S++;
    }
  }

  /* Allocate memory and initialise counters */
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

  /* Assign the lateral ghost cells to the boundary vectors */
  for (cc = 1; cc <= geom->nbpte1; cc++) {
    c = geom->bpte1[cc];        /* Global boundary cell */
    dr = dbine1[cc];            /* Direction of bin from bpt */
    if (dr == 2)
      for (zc = 1; zc <= zse1; zc++)
        c = geom->xp1[c];
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      /* It is possible that the lateral ghost cell for e1 faces */
      /* corresponds to an auxiliary cell in the sparse array (see */
      /* above) and at locations where the bathymetry is stepped, */
      /* the interior cell to this ghost cell may be outside the */
      /* limits of the window. These ghost cell's are never used in */
      /* the window's computations, so leave these cells self- */
      /* mapping.  */
      ic1 = maps[dr][lc];       /* Interior cell to lc */
      window->bpte1[nc] = lc;
      window->bine1[nc] = ic1;
      nc++;
      if (cc <= geom->nbe1)
        window->nbe1++;
      if (lc <= window->enonS) {
        window->bpte1S[ncS] = lc;
        window->bine1S[ncS] = ic1;
        ncS++;
        if (cc <= geom->nbe1S)
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
  /* Get the lateral ghost vectors for e2 face centered cells */
  if (DEBUG("init_w"))
    dlog("init_w", "    Start e2 centered local ghost vectors\n");
  /* Count the local e2 ghost cells in window n */
  window->nbpte2 = 0;
  window->nbpte2S = 0;
  for (cc = 1; cc <= geom->nbpte2; cc++) {
    c = geom->bpte2[cc];        /* Global boundary cell */
    dr = dbine2[cc];            /* Direction of bin from bpt */
    if (dr == 0)
      for (zc = 1; zc <= zse2; zc++)
        c = geom->yp1[c];
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      window->nbpte2++;
      if (lc <= window->enonS)
        window->nbpte2S++;
    }
  }

  /* Allocate memory and initialise counters */
  window->bpte2 = i_alloc_1d(window->nbpte2 + 1);
  memset(window->bpte2, 0, window->nbpte2 * sizeof(int));
  window->bine2 = i_alloc_1d(window->nbpte2 + 1);
  memset(window->bine2, 0, window->nbpte2 * sizeof(int));
  window->bpte2S = i_alloc_1d(window->nbpte2S + 1);
  memset(window->bpte2S, 0, window->nbpte2S * sizeof(int));
  window->bine2S = i_alloc_1d(window->nbpte2S + 1);
  memset(window->bine2S, 0, window->nbpte2S * sizeof(int));
  nc = ncS = 1;
  window->nbe2 = window->nbe2S = 0;

  /* Assign the lateral ghost cells to the boundary vectors */
  for (cc = 1; cc <= geom->nbpte2; cc++) {
    int gc;
    c = gc = geom->bpte2[cc];        /* Global boundary cell */
    dr = dbine2[cc];            /* Direction of bin from bpt */
    if (dr == 0)
      for (zc = 1; zc <= zse2; zc++)
        c = geom->yp1[c];
    lc = ANYf(c, window->wsa, window->enon); /* Local boundary cell */
    if (lc) {
      ic1 = maps[dr][lc];       /* Interior cell to lc */
      window->bpte2[nc] = lc;
      window->bine2[nc] = ic1;
      nc++;
      if (cc <= geom->nbe2)
        window->nbe2++;
      if (lc <= window->enonS) {
        window->bpte2S[ncS] = lc;
        window->bine2S[ncS] = ic1;
        ncS++;
        if (cc <= geom->nbe2S)
          window->nbe2S++;
      }
    }
  }
  if (nc - 1 != window->nbpte2)
    hd_warn
      ("get_local_ghosts: Incorrect number of 3D e2 ghost cells located, window %d : %d %d\n",
       n, nc, window->nbpte2);
  if (ncS - 1 != window->nbpte2S)
    hd_warn
      ("get_local_ghosts: Incorrect number of 2D e2 ghost cells located, window %d : %d %d\n",
       n, ncS, window->nbpte2S);

  /* Sediment ghosts under auxiliary and lateral ghost cells */
  window->ngsed = 0;
  memset(xpyp1, 0, window->enon + 1);
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    xpyp1[c] = 1;
  }
  for (c = 1; c <= window->enonS; c++) {
    if (!xpyp1[c]) window->ngsed++;
  }
  cc = 1;
  window->gsed_t = i_alloc_1d(window->ngsed + 1);
  for (c = 1; c <= window->enonS; c++) {
    if (!xpyp1[c]) {
      lc = c;
      while(lc != window->zm1[lc])
	lc = window->zm1[lc];
      window->gsed_t[cc] = lc;
      cc++;
    }
  }

  i_free_1d(dbine1);
  i_free_1d(dbine2);
  i_free_1d(xpyp1);
  i_free_1d(xpym1);
  i_free_1d(xmyp1);
  i_free_1d(xmym1);
  free((void *)maps);

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
void surfbot_build(geometry_t *window,  /* Window data structure */
                   int wn,      /* Window number */
                   int *sur,    /* Surface cell vector */
                   int *bot,    /* Bottom cell vector */
                   int *vec,    /* Global bottom cells vector */
                   int nvec,    /* Size of vec */
                   int *vec2D,  /* 2D cells to process vector */
                   int nvec2D,  /* Number of wet cells in vec2D */
                   int avec2D,  /* Wet + auxiliary cells in vec2D */
                   int mode     /* Cell center, e1 or e2 face code */
  )
{
  int c;                        /* Local sparse coordinate */
  int gc, gcb;                  /* Global sparse coordinate */
  int cc, c1, n;                /* Counters */
  int back;                     /* Cell in a backward direction from c */
  int *map;                     /* Local backward map */
  int *gmap;                    /* Global backward map */
  int zc;                       /* Zoom code for cell center or cell face */
  int zmf, zs;

  /* Get the zoom code and auxiliary mapping direction */
  if (mode == 0) {
    zc = ZC;
    map = window->zm1;
    gmap = NULL;
    zmf = window->zoomf;
  } else if (mode == 1) {
    zc = ZE1;
    map = window->xm1;
    gmap = geom->xm1;
    zmf = window->zmfe1;
  } else {
    zc = ZE2;
    map = window->ym1;
    gmap = geom->ym1;
    zmf = window->zmfe2;
  }
  zs = (int)(zmf / 2);

  /* Get the surface map. This is simply a copy of the 2D wet cells */
  /* to process vector and is only included here for consistency */
  /* since it is overwritten each time-step.  */
  for (cc = 1; cc <= avec2D; cc++) {
    sur[cc] = vec2D[cc];
  }

  /* Get the bottom map. This is located via the global bottom maps. */
  for (cc = 1; cc <= nvec2D; cc++)
    bot[cc] = 0;
  c1 = 1;
  for (cc = 1; cc <= nvec; cc++) {
    c = gc = vec[cc];
    if (geom->zoomc[geom->m2d[c]] & zc && geom->fm[c].wn == wn) {
      /* For zoomf > 1 bottom cells, the faces must lie at the       */
      /* the vertical level corresponding to the cell center, not    */
      /* the cell face. Hence move these cells upwards until they    */
      /* are associated with a wet cell center.                      */
      if (zmf > 1 && mode)
	do {
	  for (n = 1, back = c; n <= zs + 1; n++)
	    back = gmap[back];
	  if (geom->fm[back].wn == 0)
	    c = geom->zp1[c];
	} while(geom->fm[back].wn == 0 && c != geom->zp1[c]);
      /* back is a land cell => use the original c */
      if (c == geom->zp1[c]) c = gc;
 
      bot[c1] = geom->fm[c].sc;
      c1++;
    }
  }

  /* Add auxiliary location bottom coordinates */
  for (cc = nvec2D + 1; cc <= avec2D; cc++) {
    c = vec2D[cc];              /* Local cell to process */
    back = map[c];              /* Backward map of local cell */
    gc = window->wsa[c];        /* Global cell corresponding to c */
    gcb = window->wsa[back];    /* Global cell corresponding to back */

    /* Bottom coordinate for auxiliary cells */
    while (geom->fm[gcb].wn != 0 && geom->fm[gc].wn != 0) {
      c = window->zm1[c];
      back = map[c];
      gc = window->wsa[c];
      gcb = window->wsa[back];
    }

    /* Bottom coordinate for boundary cells adjacent to OUTSIDE */
    /* cells in the domain interior.  */
    while ((mode == 1 || mode == 2) && geom->fm[gc].wn == wn) {
      c = window->zm1[c];
      gc = window->wsa[c];
    }

    if (!mode)
      bot[cc] = c;
    else
      bot[cc] = window->zp1[c];
  }

  /*
     for(cc=nvec2D+1; cc<=avec2D; cc++) { if(!bot[cc]) { c=vec2D[cc];
     zm1=window->zm1[c]; while(c!=zm1) { c=zm1; zm1=window->zm1[zm1]; }
     if(!mode) bot[cc]=c; else bot[cc]=window->zp1[c]; } } */
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
  int zc;             /* Zoom code for cell center or cell face      */
  int *vec;           /* Local cells to process vector               */
  int nvec;           /* Size of vec                                 */
  int *lv;            /* Local 2D vector to fill                     */

  /* Allocate memory                                                 */
  lv = i_alloc_1d(window->enonS + 1);
  memset(lv, 0, (window->enonS + 1) * sizeof(int));

  /* Get the zoom code and auxiliary mapping direction */
  if (mode == 0) {
    zc = ZC;
    vec = window->w2_t;
    nvec = window->b2_t;
  } else if (mode == 1) {
    zc = ZE1;
    vec = window->w2_e1;
    nvec = window->b2_e1;
  } else {
    zc = ZE2;
    vec = window->w2_e2;
    nvec = window->b2_e2;
  }

  /* Fill the window vector from the global vector                   */
  for (cc = 1; cc <= nvec; cc++) {
    c = vec[cc];
    cg = window->wsa[c];
    if (geom->zoomc[geom->m2d[cg]] & zc && geom->fm[cg].wn == wn) {
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
void reorder_gl_map(int wn,     /* Window number */
                    int *wsa,   /* Local work sparse array */
                    int wsize,  /* Size of wsa */
                    int *mask,  /* Array of global locations to order
                                   first */
                    int msize   /* Size of mask */
  )
{
  int c, cc, c1, c2;            /* Counters */

  /* Reset the wet cells in the mask to the first locations in the */
  /* local sparse array. If mask is the surface layer then this is */
  /* not necessary since the surface layer wet cells already occupy */
  /* the first local sparse locations of fm.sc.  */
  c1 = 1;
  for (cc = 1; cc <= msize; cc++) {
    c = mask[cc];
    if (geom->fm[c].wn == wn) {
      geom->fm[c].sc = c1;
      c1++;
    }
  }

  /* Reset the wet cells not in the mask to the local locations from */
  /* msize onwards.  */
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
/* Routine to initialise the window data structure                   */
/*-------------------------------------------------------------------*/
window_t **win_data_build(master_t *master, /* Model data structure */
                          geometry_t **window /* Window data structure */
  )
{
  window_t **windat;            /* Window data structure */
  int nwindows;                 /* Number of windows */
  int n, cc, tn, k = 0;         /* Counters */
  int c, cs;                    /* Sparse locations */

  /*-----------------------------------------------------------------*/
  /* Allocate memory */
  nwindows = master->nwindows;
  windat = (window_t **)malloc(sizeof(window_t *) * (nwindows + 1));

  for (n = 1; n <= nwindows; n++) {
    windat[n] = win_data_init(master, window[n]);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise window data */
  for (n = 1; n <= nwindows; n++) {

    /* Constants */
    windat[n]->dt = master->grid_dt;
    windat[n]->iratio = master->iratio;
    windat[n]->t = master->t;
    /*windat[n]->etarlxtc = master->etarlxtc;*/
    if (!(master->cfl & NONE))
      windat[n]->mcfl2d = windat[n]->mcfl3d = HUGE;

    /* Set the 3D water column tracers pointers */
    windat[n]->sal = windat[n]->temp = NULL;
    windat[n]->tke = windat[n]->diss = windat[n]->L = windat[n]->omega =
      NULL;
    windat[n]->Q2 = windat[n]->Q2L = windat[n]->Kq = NULL;
    windat[n]->u1m = windat[n]->u2m = windat[n]->wm = windat[n]->Kzm = NULL;
    windat[n]->u1vm = windat[n]->u2vm = windat[n]->tempm = windat[n]->saltm = NULL;
    windat[n]->tram = NULL;
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
    windat[n]->speed_3d = windat[n]->speedd_3d = windat[n]->perc = windat[n]->energy = NULL;
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
      else if (strcmp("u1vmean", master->trname[tn]) == 0)
        windat[n]->u1vm = windat[n]->tr_wc[tn];
      else if (strcmp("u2vmean", master->trname[tn]) == 0)
        windat[n]->u2vm = windat[n]->tr_wc[tn];
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
      } else if (strcmp("current_dir_3d", master->trname[tn]) == 0) {
        windat[n]->speedd_3d = windat[n]->tr_wc[tn];
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
      } else if (strcmp("glider", master->trname[tn]) == 0) {
        windat[n]->glider = windat[n]->tr_wc[tn];
      } else if (strcmp("nprof", master->trname[tn]) == 0) {
        windat[n]->nprof = windat[n]->tr_wc[tn];
      } else if (strcmp("unit", master->trname[tn]) == 0) {
        windat[n]->unit = windat[n]->tr_wc[tn];
      } else if (strcmp("decorr_e1", master->trname[tn]) == 0) {
        windat[n]->decv1 = windat[n]->tr_wc[tn];
      } else if (strcmp("decorr_e2", master->trname[tn]) == 0) {
        windat[n]->decv2 = windat[n]->tr_wc[tn];
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
      } else if (master->dhwf & DHW_NOAA) {
	if (strcmp("dhd", master->trname[tn]) == 0)
	  windat[n]->dhd = windat[n]->tr_wc[tn];
	else if (strcmp("dhwc", master->trname[tn]) == 0)
	  windat[n]->dhwc = windat[n]->tr_wc[tn];
	else if (strcmp("dhw", master->trname[tn]) == 0)
	  windat[n]->dhw = windat[n]->tr_wc[tn];
      }
    }

    /* Initialise all 3D and 2D variables required by the window */
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
        windat[n]->wind1[cc] = master->wind1[c];
        windat[n]->wind2[cc] = master->wind2[c];
        windat[n]->windspeed[cc] = master->windspeed[c];
        windat[n]->winddir[cc] = master->winddir[c];
	if (master->meanc) windat[n]->meanc[cc] = master->meanc[c];

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

      /* e1 face centered arrays */
      windat[n]->u1[cc] = master->u1[c];
      if (cc <= window[n]->enonS) {
        windat[n]->u1av[cc] = master->u1av[c];
      }

      /* e2 face centered arrays */
      windat[n]->u2[cc] = master->u2[c];
      if (cc <= window[n]->enonS) {
        windat[n]->u2av[cc] = master->u2av[c];
      }
    }

    /* Save the bottom velocity for the 2D mode */
    for (cc = 1; cc <= window[n]->b2_e1; cc++) {
      c = window[n]->bot_e1[cc];  /* 3D bottom coordinate */
      cs = window[n]->m2d[c];   /* 2D coordiate corresponding to c */
      windat[n]->u1bot[cs] = windat[n]->u1[c] - windat[n]->u1av[cs];
    }
    for (cc = 1; cc <= window[n]->b2_e2; cc++) {
      c = window[n]->bot_e2[cc];  /* 3D bottom coordinate */
      cs = window[n]->m2d[c];   /* 2D coordiate corresponding to c */
      windat[n]->u2bot[cs] = windat[n]->u2[c] - windat[n]->u2av[cs];
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
window_t *win_data_init(master_t *master, /* Master data structure */
                        geometry_t *window  /* Window data structure */
  )
{
  window_t *windat;             /* Window data structure */
  int winsize;                  /* Number of sparse cells in window */
  int n, m;

  windat = win_data_alloc();

  /* Multiple windows. Allocate memory for each window */
  if (master->nwindows > 1) {
    /* 3D arrays */
    winsize = window->enon + 1;
    windat->ntr = master->ntr;
    windat->tr_wc = d_alloc_2d(winsize, windat->ntr);
    windat->u1 = d_alloc_1d(winsize);
    windat->u2 = d_alloc_1d(winsize);
    windat->w = d_alloc_1d(winsize);
    windat->Kz = d_alloc_1d(winsize);
    windat->Vz = d_alloc_1d(winsize);
    windat->dzu1 = d_alloc_1d(winsize);
    windat->dzu2 = d_alloc_1d(winsize);
    windat->u1flux3d = d_alloc_1d(winsize);
    windat->u2flux3d = d_alloc_1d(winsize);
    windat->dens = d_alloc_1d(winsize);
    windat->dens_0 = d_alloc_1d(winsize);
    windat->u1b = d_alloc_1d(winsize);
    windat->u2b = d_alloc_1d(winsize);
    windat->waterss = d_alloc_1d(winsize);
    if (master->velrlx & RELAX) {
      relax_info_t *relax = master->vel_rlx;
      windat->vel_rlx = relax_info_init(relax->rlxn, relax->rlxtc, 
					relax->rlxdt, winsize, winsize);
      windat->vel_rlx->rlx = relax->rlx;
    }

    /* 2D arrays */
    windat->ntrS = master->ntrS;
    winsize = window->enonS + 1;
    if (windat->ntrS)
      windat->tr_wcS = d_alloc_2d(winsize, windat->ntrS);
    windat->eta = d_alloc_1d(winsize);
    windat->u1av = d_alloc_1d(winsize);
    windat->u2av = d_alloc_1d(winsize);
    windat->nu1av = d_alloc_1d(winsize);
    windat->nu2av = d_alloc_1d(winsize);
    windat->u1flux = d_alloc_1d(winsize);
    windat->u2flux = d_alloc_1d(winsize);
    windat->depth_e1 = d_alloc_1d(winsize);
    windat->depth_e2 = d_alloc_1d(winsize);
    windat->u1bot = d_alloc_1d(winsize);
    windat->u2bot = d_alloc_1d(winsize);
    windat->topz = d_alloc_1d(winsize);
    windat->wdiff2d = d_alloc_1d(winsize);
    windat->detadt = d_alloc_1d(winsize);
    windat->wtop = d_alloc_1d(winsize);
    windat->wbot = d_alloc_1d(winsize);
    windat->patm = d_alloc_1d(winsize);
    windat->wind1 = d_alloc_1d(winsize);
    windat->wind2 = d_alloc_1d(winsize);
    windat->windspeed = d_alloc_1d(winsize);
    windat->winddir   = d_alloc_1d(winsize);
    windat->waterss2d = d_alloc_1d(winsize);

    windat->etab = d_alloc_1d(winsize);
    windat->u1avb = d_alloc_1d(winsize);
    windat->u2avb = d_alloc_1d(winsize);

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

    if (master->dozoom) {
      windat->u1flux_z = d_alloc_1d(winsize);
      windat->u2flux_z = d_alloc_1d(winsize);
    }
    if (!(master->means & NONE)) {
      windat->meanc = d_alloc_1d(geom->sgsizS);
      if (master->means & TIDAL) {
	windat->odeta = d_alloc_1d(geom->sgsizS);
      }
    }
    /*UR-TODO why is this overwritten from above ?? 
    windat->meanc = d_alloc_1d(winsize);
    windat->odeta = d_alloc_1d(winsize);
    */
    if (master->ntr) windat->trinc = i_alloc_1d(master->ntr);
    if (master->ntrS) windat->trincS = i_alloc_1d(master->ntrS);

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

    /* Allocate salt flux arrays, if needed */
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

    windat->nsed = master->nsed;
    if (windat->nsed)
      windat->tr_sed = d_alloc_3d(winsize, window->sednz, windat->nsed);

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
	  /* See logic in heatflux.c:comp_heat_mom */
	  windat->swr = windat->tr_wcS[m];
	} else 
	  windat->swrd = windat->tr_wcS[m];
      }
      /*
       * These heatflux tracer numbers reflect logic in load_tracer
       */
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
      if (strcmp("current_dir_2d", master->trinfo_2d[m].name) == 0)
        windat->speedd_2d = windat->tr_wcS[m];
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
      if (strcmp("surf_slope_x", master->trinfo_2d[m].name) == 0)
        windat->slope_x = windat->tr_wcS[m];
      if (strcmp("surf_slope_y", master->trinfo_2d[m].name) == 0)
        windat->slope_y = windat->tr_wcS[m];
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
      if (strcmp("decorr_e2", master->trinfo_2d[m].name) == 0)
        windat->decv2 = windat->tr_wcS[m];
      if (strcmp("sep", master->trinfo_2d[m].name) == 0)
        windat->sep = windat->tr_wcS[m];
      if (strcmp("bep", master->trinfo_2d[m].name) == 0)
        windat->bep = windat->tr_wcS[m];
      /*if (strcmp("oeta", master->trinfo_2d[m].name) == 0 && windat->eta_rlx)
        windat->eta_rlx->val1 = windat->tr_wcS[m];*/
    }
  }

  /* Single window. The master and the window are the same - point */
  /* the window to the master data structure arrays.  */
  else {
    /* 3D arrays */
    windat->ntr = master->ntr;
    windat->tr_wc = master->tr_wc;
    windat->nsed = master->nsed;
    if (windat->nsed)
      windat->tr_sed = master->tr_sed;
    windat->u1 = master->u1;
    windat->u2 = master->u2;
    windat->w = master->w;
    windat->Kz = master->Kz;
    windat->Vz = master->Vz;
    windat->dzu1 = master->dzu1;
    windat->dzu2 = master->dzu2;
    windat->u1flux3d = master->u1flux3d;
    windat->u2flux3d = master->u2flux3d;
    windat->dens = master->dens;
    windat->dens_0 = master->dens_0;
    windat->u1b = master->u1b;
    windat->u2b = master->u2b;
    windat->waterss = master->waterss;
    if (master->velrlx & RELAX) 
      windat->vel_rlx = master->vel_rlx;

    /* 2D arrays */
    windat->ntrS = master->ntrS;
    if (windat->ntrS)
      windat->tr_wcS = master->tr_wcS;
    windat->eta = master->eta;
    windat->u1av = master->u1av;
    windat->u2av = master->u2av;
    windat->nu1av = master->nu1av;
    windat->nu2av = master->nu2av;
    windat->u1flux = master->u1flux;
    windat->u2flux = master->u2flux;
    windat->depth_e1 = master->depth_e1;
    windat->depth_e2 = master->depth_e2;
    windat->u1bot = master->u1bot;
    windat->u2bot = master->u2bot;
    windat->topz = master->topz;
    windat->wdiff2d = master->wdiff2d;
    windat->detadt = master->detadt;
    windat->wtop = master->wtop;
    windat->wbot = master->wbot;
    windat->wind1 = master->wind1;
    windat->wind2 = master->wind2;
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
    windat->u1flux_z = master->u1flux_z;
    windat->u2flux_z = master->u2flux_z;

    windat->etab = master->etab;
    windat->u1avb = master->u1avb;
    windat->u2avb = master->u2avb;

    windat->heatf = master->heatf;
    windat->swr = master->swr;
    windat->light = master->light;
    windat->hftemp = master->hftemp;

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
    windat->speedd_2d = master->speedd_2d;
    windat->speed_3d = master->speed_3d;
    windat->speedd_3d = master->speedd_3d;
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
    windat->slope_y = master->slope_y;
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

    /* Diagnostic indicies */
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
  winsize = window->enonS + 1;
  windat->sur_e1 = i_alloc_1d(winsize);
  windat->sur_e2 = i_alloc_1d(winsize);
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
void win_data_clear(window_t *windat  /* Window data structure */
  )
{

  /* 3D arrays */
  d_free_2d(windat->tr_wc);
  d_free_1d(windat->u1);
  d_free_1d(windat->u2);
  d_free_1d(windat->w);
  d_free_1d(windat->Kz);
  d_free_1d(windat->Vz);
  d_free_1d(windat->dzu1);
  d_free_1d(windat->dzu2);
  d_free_1d(windat->u1flux3d);
  d_free_1d(windat->u2flux3d);
  d_free_1d(windat->dens);
  d_free_1d(windat->dens_0);
  d_free_1d(windat->u1b);
  d_free_1d(windat->u2b);
  d_free_1d(windat->waterss);

  /* 2D arrays */
  if (windat->ntrS)
    d_free_2d(windat->tr_wcS);
  d_free_1d(windat->eta);
  d_free_1d(windat->u1av);
  d_free_1d(windat->u2av);
  d_free_1d(windat->etab);
  d_free_1d(windat->u1avb);
  d_free_1d(windat->u2avb);
  i_free_1d(windat->sur_e1);
  i_free_1d(windat->sur_e2);
  d_free_1d(windat->depth_e1);
  d_free_1d(windat->depth_e2);
  d_free_1d(windat->u1flux);
  d_free_1d(windat->u2flux);
  d_free_1d(windat->u1bot);
  d_free_1d(windat->u2bot);
  d_free_1d(windat->topz);
  d_free_1d(windat->wdiff2d);
  d_free_1d(windat->detadt);
  d_free_1d(windat->wtop);
  d_free_1d(windat->wbot);
  d_free_1d(windat->patm);
  d_free_1d(windat->wind1);
  d_free_1d(windat->wind2);
  /*UR-ADDED; for one window this is also clearerd at wind forcing */
  if(windat->windspeed != NULL)
    d_free_1d(windat->windspeed);
  if(windat->winddir != NULL)
    d_free_1d(windat->winddir);
  d_free_1d(windat->waterss2d);
  if (windat->u1flux_z)
    d_free_1d(windat->u1flux_z);
  if (windat->u2flux_z)
    d_free_1d(windat->u2flux_z);
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
  if (windat->cu2)
    i_free_1d(windat->cu2);
  if (windat->cu1a)
    i_free_1d(windat->cu1a);
  if (windat->cu2a)
    i_free_1d(windat->cu2a);
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
  /* Space for these variables is deallocated by the master */
  /* via the scheduler cleanup for 1 window.                */
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

  win_data_init_transfer_buf_cleanup(windat);

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
  geometry_t *geom = master->geom;  /* Global geometry */
  win_priv_t **wincon;          /* Window constants structure */
  window_t *windat;             /* Window data structure */
  int nwindows;                 /* Number of windows */
  int n, cc, tn;                /* Counters */
  int c, c1, cb;                /* Sparse locations */
  int winsize;                  /* Number of sparse cells in window */
  int i, j, k;                  /* Cartesian locations */

  /*-----------------------------------------------------------------*/
  /* Allocate memory for the wincon arrays */
  nwindows = master->nwindows;
  wincon = (win_priv_t **)malloc(sizeof(win_priv_t *) * (nwindows + 1));
  for (n = 1; n <= nwindows; n++) {
    windat = window[n]->windat;
    wincon[n] = win_consts_alloc();
    /* window[n]->geom = wincon[n]; */
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

    /* 3D work arrays */
    winsize = window[n]->sgsiz = window[n]->enon + 1;
    wincon[n]->w1 = d_alloc_1d(winsize);
    wincon[n]->w2 = d_alloc_1d(winsize);
    wincon[n]->w3 = d_alloc_1d(winsize);
    wincon[n]->w4 = d_alloc_1d(winsize);
    /* Also allocate the parallel versions, if required */
#ifdef HAVE_OMP
    wincon[n]->trans_num_omp = master->trans_num_omp;
    if (wincon[n]->trans_num_omp > 1) {
      wincon[n]->w1n = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->w2n = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->w3n = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->w4n = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
    }
#endif
    wincon[n]->w5 = d_alloc_1d(winsize);
    wincon[n]->w6 = d_alloc_1d(winsize);
    wincon[n]->w7 = d_alloc_1d(winsize);
    wincon[n]->w8 = d_alloc_1d(winsize);
    wincon[n]->w9 = d_alloc_1d(winsize);
    wincon[n]->w10 = d_alloc_1d(winsize);
    wincon[n]->rdens = d_alloc_1d(winsize);
    wincon[n]->s1 = i_alloc_1d(winsize);
    wincon[n]->s2 = i_alloc_1d(winsize);
    wincon[n]->s3 = i_alloc_1d(winsize);
    wincon[n]->s4 = i_alloc_1d(winsize);
    wincon[n]->c1 = c_alloc_1d(winsize);
    wincon[n]->gmap = i_alloc_1d(winsize);
    memset(wincon[n]->s1, 0, winsize * sizeof(int));
    memset(wincon[n]->s2, 0, winsize * sizeof(int));
    if (windat->u1_adv || windat->u1_hdif || windat->u1_vdif ||
        windat->u1_cor || windat->u1_btp || windat->u1_bcp ||
        windat->u2_adv || windat->u2_hdif || windat->u2_vdif ||
        windat->u2_cor || windat->u2_btp || windat->u2_bcp ||
	windat->tr_adv || windat->tr_hdif || windat->tr_vdif ||
	windat->tr_ncon)
      wincon[n]->tendency = d_alloc_1d(winsize);
    if(master->trasc & (LAGRANGE|FFSL) || master->momsc & LAGRANGE) {
      wincon[n]->m2d = i_alloc_1d(winsize);
      wincon[n]->s5 = i_alloc_1d(winsize);
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
    wincon[n]->tr_mod = d_alloc_1d(winsize);
    wincon[n]->tr_mod_x = d_alloc_1d(winsize);
    wincon[n]->tr_mod_y = d_alloc_1d(winsize);
    wincon[n]->tr_mod_z = d_alloc_1d(winsize);
    /* Also allocate the parallel versions, if required */
#ifdef HAVE_OMP
    if (wincon[n]->trans_num_omp > 1) {
      wincon[n]->tr_modn = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_x = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_y = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
      wincon[n]->tr_modn_z = d_alloc_2d(winsize, wincon[n]->trans_num_omp);
    }
#endif
    /* 2D work arrays */
    winsize = window[n]->sgsizS = window[n]->enonS + 1;
    wincon[n]->d1 = d_alloc_1d(winsize);
    wincon[n]->d2 = d_alloc_1d(winsize);
    wincon[n]->d3 = d_alloc_1d(winsize);
    wincon[n]->d4 = d_alloc_1d(winsize);
    wincon[n]->d5 = d_alloc_1d(winsize);
    wincon[n]->d6 = d_alloc_1d(winsize);
    wincon[n]->d7 = d_alloc_1d(winsize);
    wincon[n]->c2 = c_alloc_1d(winsize);
    if (master->thin_merge) {
      wincon[n]->kth_e1 = i_alloc_1d(winsize);
      wincon[n]->kth_e2 = i_alloc_1d(winsize);
    }
    wincon[n]->cdry_e1 = i_alloc_1d(winsize);
    wincon[n]->cdry_e2 = i_alloc_1d(winsize);
    wincon[n]->cbot_e1 = i_alloc_1d(winsize);
    wincon[n]->cbot_e2 = i_alloc_1d(winsize);
    wincon[n]->i1 = i_alloc_1d(winsize);
    wincon[n]->i2 = i_alloc_1d(winsize);
    wincon[n]->i3 = i_alloc_1d(winsize);
    wincon[n]->i4 = i_alloc_1d(winsize);
    wincon[n]->i5 = i_alloc_1d(winsize);
    wincon[n]->i6 = i_alloc_1d(winsize);
    wincon[n]->i7 = i_alloc_1d(winsize);
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

    /* 1D work arrays */
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

    /* Sigma arrays */
    if (master->sigma) {
      winsize = window[n]->sgsizS;
      wincon[n]->mdx = d_alloc_1d(winsize);
      wincon[n]->mdy = d_alloc_1d(winsize);
    }

    /* Multiple windows. Allocate memory for each window */
    if (master->nwindows > 1) {
      /* 3D arrays */
      winsize = window[n]->sgsiz;
      wincon[n]->dz = d_alloc_1d(winsize);
      wincon[n]->t11 = d_alloc_1d(winsize);
      wincon[n]->t12 = d_alloc_1d(winsize);
      wincon[n]->t22 = d_alloc_1d(winsize);
      wincon[n]->q = (master->q != NULL) ? d_alloc_1d(winsize) : NULL;
      if (master->smagcode & (U1_A|U1_SA))
        wincon[n]->u1vh = d_alloc_1d(winsize);
      else
        wincon[n]->u1vh = windat->sdc;
      if (master->smagcode & (U2_A|U2_SA))
        wincon[n]->u2vh = d_alloc_1d(winsize);
      else
        wincon[n]->u2vh = windat->sdc;

      if (master->smagcode & (U1_AK|U1_SAK))
        wincon[n]->u1kh = d_alloc_1d(winsize);
      else
	wincon[n]->u1kh = windat->sdc;

      if (master->smagcode & (U2_AK|U2_SAK))
        wincon[n]->u2kh = d_alloc_1d(winsize);
      else
	wincon[n]->u2kh = windat->sdc;

      /* 2D arrays */
      winsize = window[n]->sgsizS;
      wincon[n]->oldeta = d_alloc_1d(winsize);
      wincon[n]->one = d_alloc_1d(winsize);
      for (c = 1; c <= window[n]->enonS; c++)
        wincon[n]->one[c] = 1.0;
      if (master->sigma) {
        wincon[n]->Ds = d_alloc_1d(winsize);
        wincon[n]->Hs = d_alloc_1d(winsize);
        wincon[n]->Hn1 = d_alloc_1d(winsize);
        wincon[n]->Hn2 = d_alloc_1d(winsize);
      }
      wincon[n]->Cd = d_alloc_1d(winsize);
      wincon[n]->z0 = d_alloc_1d(winsize);
      wincon[n]->topdensu1 = d_alloc_1d(winsize);
      wincon[n]->densavu1 = d_alloc_1d(winsize);
      wincon[n]->topdensu2 = d_alloc_1d(winsize);
      wincon[n]->densavu2 = d_alloc_1d(winsize);
      wincon[n]->u1c1 = d_alloc_1d(winsize);
      wincon[n]->u1c3 = d_alloc_1d(winsize);
      wincon[n]->u1c4 = d_alloc_1d(winsize);
      wincon[n]->u1c5 = d_alloc_1d(winsize);
      wincon[n]->u1c6 = d_alloc_1d(winsize);
      wincon[n]->u2c1 = d_alloc_1d(winsize);
      wincon[n]->u2c3 = d_alloc_1d(winsize);
      wincon[n]->u2c4 = d_alloc_1d(winsize);
      wincon[n]->u2c5 = d_alloc_1d(winsize);
      wincon[n]->u2c6 = d_alloc_1d(winsize);
      wincon[n]->u1adv = d_alloc_1d(winsize);
      wincon[n]->u2adv = d_alloc_1d(winsize);
      wincon[n]->u1inter = d_alloc_1d(winsize);
      wincon[n]->u2inter = d_alloc_1d(winsize);
      wincon[n]->coriolis = d_alloc_1d(winsize);
    } else {
      /* Single window. The master and the window are the same - */
      /* point the window to the master data structure arrays.  */
      /* 3D arrays */
      wincon[n]->dz = master->dz;
      wincon[n]->t11 = master->t11;
      wincon[n]->t12 = master->t12;
      wincon[n]->t22 = master->t22;
      wincon[n]->q = master->q;
      wincon[n]->u1vh = master->u1vh;
      wincon[n]->u2vh = master->u2vh;
      wincon[n]->u1kh = master->u1kh;
      wincon[n]->u2kh = master->u2kh;

      /* 2D arrays */
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
      wincon[n]->topdensu2 = master->topdensu2;
      wincon[n]->densavu2 = master->densavu2;
      wincon[n]->u1c1 = master->u1c1;
      wincon[n]->u1c3 = master->u1c3;
      wincon[n]->u1c4 = master->u1c4;
      wincon[n]->u1c5 = master->u1c5;
      wincon[n]->u1c6 = master->u1c6;
      wincon[n]->u2c1 = master->u2c1;
      wincon[n]->u2c3 = master->u2c3;
      wincon[n]->u2c4 = master->u2c4;
      wincon[n]->u2c5 = master->u2c5;
      wincon[n]->u2c6 = master->u2c6;
      wincon[n]->u1adv = master->u1adv;
      wincon[n]->u2adv = master->u2adv;
      wincon[n]->u1inter = master->u1inter;
      wincon[n]->u2inter = master->u2inter;
      wincon[n]->coriolis = master->coriolis;
    }

    /*---------------------------------------------------------------*/
    /* Initialise the wincon arrays */
    wincon[n]->g = master->g;
    wincon[n]->tz = tm_tz_offset(master->timeunit);
    wincon[n]->ambpress = master->ambpress;
    wincon[n]->trasc = master->trasc;
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
    wincon[n]->sue2 = master->sue2;
    wincon[n]->kue1 = master->kue1;
    wincon[n]->kue2 = master->kue2;
    wincon[n]->bsue1 = master->bsue1;
    wincon[n]->bsue2 = master->bsue2;
    wincon[n]->bkue1 = master->bkue1;
    wincon[n]->bkue2 = master->bkue2;
    wincon[n]->diff_scale = master->diff_scale;
    wincon[n]->smag_smooth = master->smag_smooth;
    wincon[n]->visc_method = master->visc_method;
    wincon[n]->stab = master->stab;
    wincon[n]->cfl = master->cfl;
    wincon[n]->cfl_dt = master->cfl_dt;
    wincon[n]->lnm = master->lnm;
    wincon[n]->nprof = master->nprofn;
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
    wincon[n]->means_tra = master->means_tra;
    wincon[n]->tendf = master->tendf;
    wincon[n]->trtend = master->trtend;
    wincon[n]->trsplit = master->trsplit;
    wincon[n]->trflux = master->trflux;
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
    wincon[n]->u2_f = master->u2_f;
    wincon[n]->u1av_f = master->u1av_f;
    wincon[n]->u2av_f = master->u2av_f;
    wincon[n]->save_force = master->save_force;
    wincon[n]->vinit = master->vinit;
    wincon[n]->etarlx = master->etarlx;
    wincon[n]->velrlx = master->velrlx;
    wincon[n]->alertf = master->alertf;
    wincon[n]->u1_f = master->u1_f;
    wincon[n]->u2_f = master->u2_f;
    wincon[n]->u1av_f = master->u1av_f;
    wincon[n]->u2av_f = master->u2av_f;
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
    wincon[n]->dozoom = master->dozoom;
    wincon[n]->exmapf = master->exmapf;
    wincon[n]->waves = master->waves;
    wincon[n]->orbital = master->orbital;
    wincon[n]->fillf = master->fillf;
    wincon[n]->pssinput = master->pssinput;
    wincon[n]->conserve = master->conserve;
    wincon[n]->do_closure = master->do_closure;
    wincon[n]->tmode = master->tmode;
    wincon[n]->trout = master->trout;
    wincon[n]->gint_errfcn = master->gint_errfcn;
    wincon[n]->da = master->da;
    wincon[n]->swr_type = master->swr_type;
    wincon[n]->dhwf = master->dhwf;

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
  /* Initialise tracer info */
  for (n = 1; n <= nwindows; n++) {

    if (nwindows > 1) {
      /* 3D Tracers */
      wincon[n]->trinfo_3d = (tracer_info_t *)malloc(sizeof(tracer_info_t) *
                                                   master->ntr);
      memset(wincon[n]->trinfo_3d, 0, sizeof(tracer_info_t) * master->ntr);
      for (tn = 0; tn < master->ntr; tn++) {
        tracer_copy(&wincon[n]->trinfo_3d[tn], &master->trinfo_3d[tn]);
        strcpy(wincon[n]->trinfo_3d[tn].name, master->trinfo_3d[tn].name);
      }

      /* 2D Tracers */
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

      /* Sediment tracers */
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

    /* Get the tracers to advect and diffuse for the window */
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
    /*
       for(tn=0; tn<wincon[n]->ntr; tn++) { if(master->diffuse[tn]) { c=1;
       for(cc=0; cc<wincon[n]->ntr; cc++) { if(master->advect[cc])c=0; }
       if(c)wincon[n]->ntbdy++; } } */

    /* Get the tracers which need OBC's to be set */
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

    /* Get the tracers that have a surface flux tracer set */
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

    /* Add the tracers that are diffused, and which are not already */
    /* advected.  */
    /*
       for(tn=0; tn<wincon[n]->ntr; tn++) { if(master->diffuse[tn]) { c=1;
       for(cc=0; cc<wincon[n]->ntr; cc++) { if(tn==wincon[n]->tbdy[cc])c=0;
       } if(c) { wincon[n]->tbdy[wincon[n]->ntbdy]=tn; wincon[n]->ntbdy++; }
       } } */
    /* Get the tracers which are to be decayed & diagnostic tracers */
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
  /* Initialise window data */
  for (n = 1; n <= nwindows; n++) {

    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];

      /* If a ghost cell is encountered, save it's wet neighbour in */
      /* wincon->s2 and continue.  */
      wincon[n]->s2[cc] = 0;
      if (geom->fm[c].wn != n &&
	  i == NOTVALID && j == NOTVALID && k == NOTVALID) {
	/*
        if (geom->fm[geom->xp1[c]].wn != 0)
          wincon[n]->s2[cc] = geom->fm[geom->xp1[c]].sc;
        else if (geom->fm[geom->xm1[c]].wn != 0)
          wincon[n]->s2[cc] = geom->fm[geom->xm1[c]].sc;
        else if (geom->fm[geom->yp1[c]].wn != 0)
          wincon[n]->s2[cc] = geom->fm[geom->yp1[c]].sc;
        else if (geom->fm[geom->ym1[c]].wn != 0)
          wincon[n]->s2[cc] = geom->fm[geom->ym1[c]].sc;
	*/
	wincon[n]->s2[cc] = window[n]->wgst[cc];

        continue;
      }

      /* Cell centered 3D arrays */
      if (master->q != NULL)
        wincon[n]->q[cc] = master->q[c];
      wincon[n]->dz[cc] = master->dz[c];

      /* e1 face 3D centered arrays */
      wincon[n]->u1vh[cc] = master->u1vh[c];
      wincon[n]->u1kh[cc] = master->u1kh[c];

      /* e2 face 3D centered arrays */
      wincon[n]->u2vh[cc] = master->u2vh[c];
      wincon[n]->u2kh[cc] = master->u2kh[c];

      if (cc <= window[n]->enonS) {
        /* Cell centered 2D arrays */
        wincon[n]->Cd[cc] = master->Cd[c];
        wincon[n]->z0[cc] = master->z0[c];

        /* e1 face 2D centered arrays */
        wincon[n]->u1c1[cc] = master->u1c1[c];
        wincon[n]->u1c3[cc] = master->u1c3[c];
        wincon[n]->u1c4[cc] = master->u1c4[c];
        wincon[n]->u1c5[cc] = master->u1c5[c];
        wincon[n]->u1c6[cc] = master->u1c6[c];

        /* e2 face 2D centered arrays */
        wincon[n]->u2c1[cc] = master->u2c1[c];
        wincon[n]->u2c3[cc] = master->u2c3[c];
        wincon[n]->u2c4[cc] = master->u2c4[c];
        wincon[n]->u2c5[cc] = master->u2c5[c];
        wincon[n]->u2c6[cc] = master->u2c6[c];

	/* Coriolis */
        wincon[n]->coriolis[cc] = master->coriolis[c];

	/* Fetch */
	if (master->fetch) {
	  for (tn = 0; tn < 8; tn++)
	    wincon[n]->fetch[cc][tn] = master->fetch[c][tn];
	}
      }
    }
    consts_zoom(master, geom, window[n], wincon[n]);

    /* Initialise the sigma arrays */
    if (wincon[n]->sigma) {
      for (cc = 1; cc <= window[n]->enonS; cc++) {
        c = window[n]->wsa[cc];
        wincon[n]->Ds[cc] = master->Ds[c];
        wincon[n]->Hs[cc] = master->Hs[c];
        wincon[n]->Hn1[cc] = master->Hn1[c];
        wincon[n]->Hn2[cc] = master->Hn2[c];
      }
    } else {
      wincon[n]->mdx = wincon[n]->one;
      wincon[n]->mdy = wincon[n]->one;
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
      wincon[n]->linmask_u1 = i_alloc_1d(window[n]->sgsizS);
      memset(wincon[n]->linmask_u1, 0, window[n]->sgsizS * sizeof(int));
      for (i = 0; i < window[n]->nobc; i++) {
	open_bdrys_t *open = window[n]->open[i];
	if (open->type & U1BDRY) {
	  if (open->bcond_nor & LINEAR) {
	    for(cc = 1; cc <= open->no2_e1; cc++) {
	      c = open->obc_e1[cc];
	      wincon[n]->linmask_u1[c] = 1;
	    }
	  }
	  if (open->linear_zone_nor) {
	    for(cc = 1; cc <= open->no2_e1; cc++) {
	      c = open->oi1_e1[cc];
	      for (j=0; j < open->linear_zone_nor; j++) {
		wincon[n]->linmask_u1[c] = 1;
		c = open->nmap[c];
	      }
	    }
	  }
	}
	if (open->type & U2BDRY ) {
	  if (open->bcond_tan & LINEAR) {
	    for(cc = open->no2_e1 + 1; cc <= open->to2_e1; cc++) {
	      c = open->obc_e1[cc];
	      wincon[n]->linmask_u1[c] = 1;
	    }
	  }
	  if (open->linear_zone_tan) {
	    for(cc = open->no2_e1 + 1; cc <= open->to2_e1; cc++) {
	      c = open->oi1_e1[cc];
	      for (j=0; j < open->linear_zone_tan; j++) {
		wincon[n]->linmask_u1[c] = 1;
		c = open->nmap[c];
	      }
	    }
	  }
	}
      }
    }

    /* Set the mask for u2 linear momentum cells                     */
    if (wincon[n]->dolin_u2) {
      wincon[n]->linmask_u2 = i_alloc_1d(window[n]->sgsizS);
      memset(wincon[n]->linmask_u2, 0, window[n]->sgsizS * sizeof(int));
      for (i = 0; i < window[n]->nobc; i++) {
	open_bdrys_t *open = window[n]->open[i];
	if (open->type & U2BDRY) {
	  if (open->bcond_nor & LINEAR) {
	    for(cc = 1; cc <= open->no2_e2; cc++) {
	      c = open->obc_e2[cc];
	      wincon[n]->linmask_u2[c] = 1;
	    }
	  }
	  if (open->linear_zone_nor) {
	    for(cc = 1; cc <= open->no2_e2; cc++) {
	      c = open->oi1_e2[cc];
	      for (j=0; j < open->linear_zone_nor; j++) {
		wincon[n]->linmask_u2[c] = 1;
		c = open->nmap[c];
	      }
	    }
	  }
	}
	if (open->type & U1BDRY) {
	  if (open->bcond_tan & LINEAR) {
	    for(cc = open->no2_e2 + 1; cc <= open->to2_e2; cc++) {
	      c = open->obc_e2[cc];
	      wincon[n]->linmask_u2[c] = 1;
	    }
	  }
	  if (open->linear_zone_tan) {
	    for(cc = open->no2_e2 + 1; cc <= open->to2_e2; cc++) {
	      c = open->oi1_e2[cc];
	      for (j=0; j < open->linear_zone_tan; j++) {
		wincon[n]->linmask_u2[c] = 1;
		c = open->nmap[c];
	      }
	    }
	  }
	}
      }
    }

    /* Save the amalgamated OBC conditions for each advected tracer */
    wincon[n]->obctr = i_alloc_1d(wincon[n]->ntr * sizeof(int));
    memset(wincon[n]->obctr, 0, wincon[n]->ntr * sizeof(int));
    for (tn = 0; tn < wincon[n]->ntr; tn++) {
      for (i = 0; i < window[n]->nobc; i++) {
	open_bdrys_t *open = window[n]->open[i];
	wincon[n]->obctr[tn] |= open->bcond_tra[tn];
      }
    }
    /* Save the OBC cells associated with OUTSIDE cells */
    window[n]->nobce1 = window[n]->nobce2 = 0;
    for (i = 0; i < window[n]->nobc; i++) {
      open_bdrys_t *open = window[n]->open[i];
      if (open->ocodex & R_EDGE) window[n]->nobce1 += open->no3_t;
      if (open->ocodey & F_EDGE) window[n]->nobce2 += open->no3_t;
    }
    if (window[n]->nobce1)
      window[n]->obce1 = i_alloc_1d((window[n]->nobce1 + 1) * sizeof(int));
    if (window[n]->nobce2)
      window[n]->obce2 = i_alloc_1d((window[n]->nobce2 + 1) * sizeof(int));
    window[n]->nobce1 = window[n]->nobce2 = 0;
    for (i = 0; i < window[n]->nobc; i++) {
      open_bdrys_t *open = window[n]->open[i];
      for(cc = 1; cc <= open->no3_t; cc++) {
	c = open->omap[open->obc_t[cc]];
	if (open->ocodex & R_EDGE) {
	  window[n]->nobce1 ++;
	  window[n]->obce1[window[n]->nobce1] = c;
	}
	if (open->ocodey & F_EDGE) {
	  window[n]->nobce2 ++;
	  window[n]->obce2[window[n]->nobce2] = c;
	}
      }
    }

    /* Get the maps for semi-Lagrange advection */
    if(master->trasc & LAGRANGE || master->momsc & LAGRANGE) {
      wincon[n]->osl = master->osl;
      if (wincon[n]->momsc & LAGRANGE && wincon[n]->osl > 2)
	hd_warn("LAGRANGE: Instabilities are likely to occur using %d order semi-Lagrange.\n", 
		wincon[n]->osl);
      set_lmap(window[n], wincon[n]);
    }
    if(master->trasc & FFSL) {
      wincon[n]->osl = master->osl = 0;
      set_lmap(window[n], wincon[n]);
    }

    /* Set the offset for means */
    if (master->means_os) {
      wincon[n]->means_os = master->means_os;
      for (cc = 1; cc <= window[n]->enonS; cc++)
	windat->meanc[cc] += wincon[n]->means_os;
    }

    if (DEBUG("init_w"))
      dlog("init_w", "Window %d initialised OK\n", n);
  }

  /*-----------------------------------------------------------------*/
  /* Increase friction over coarse zoom grid */
  if (master->dozoom) {
    for (n = 1; n <= nwindows; n++) {
      if (window[n]->zmfe1 > 1 && wincon[n]->u1vh0 > 0.0) {
	for (cc = 1; cc <= window[n]->n3_e1; cc++) {
	  int xm1, xp1;
	  c = window[n]->w3_e1[cc];
	  c1 = window[n]->m2d[c];	  
	  cb = xm1 = xp1 = window[n]->wsa[c];
	  if (master->diff_scale == LINEAR)
	    wincon[n]->u1vh[c] = fabs(wincon[n]->u1vh0 *
				      window[n]->h1au1[c1] / master->hmean1);
	  if (master->diff_scale == NONLIN)
	    wincon[n]->u1vh[c] = fabs(wincon[n]->u1vh0 *
				      window[n]->h1au1[c1] * 
				      window[n]->h1au1[c1] / 
				      (master->hmean1 * master->hmean1));
	  if (cc <= window[n]->b3_e1) {
	    for (i = 1; i <= window[n]->zmee1; i++) {
	      xm1 = geom->xm1[xm1];
	      xp1 = geom->xp1[xp1];
	      master->u1vh[cb] = master->u1vh[xm1] = master->u1vh[xp1] = 
		wincon[n]->u1vh[c];
	    }
	  }
	}
      }
      
      if (window[n]->zmfe2 > 1 && wincon[n]->u2vh0 > 0.0) {
	for (cc = 1; cc <= window[n]->n3_e2; cc++) {
	  int ym1, yp1;
	  c = window[n]->w3_e2[cc];
	  c1 = window[n]->m2d[c];
	  cb = ym1 = yp1 = window[n]->wsa[c];
	  if (master->diff_scale == LINEAR)
	    wincon[n]->u2vh[c] = fabs(wincon[n]->u2vh0 *
				      window[n]->h2au2[c1] / master->hmean2);
	  if (master->diff_scale == NONLIN)
	    wincon[n]->u2vh[c] = fabs(wincon[n]->u2vh0 *
				      window[n]->h2au2[c1] * 
				      window[n]->h2au2[c1] / 
				      (master->hmean2 * master->hmean2));
	  if (cc <= window[n]->b3_e2) {
	    for (j = 1; j <= window[n]->zmee2; j++) {
	      ym1 = geom->ym1[ym1];
	      yp1 = geom->yp1[yp1];
	      master->u2vh[cb] = master->u2vh[ym1] = master->u2vh[yp1] = 
		wincon[n]->u2vh[c];
	    }
	  }
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set values in ghost locations */
  for (n = 1; n <= nwindows; n++) {
    for (cc = 1; cc <= window[n]->enonS; cc++) {
      c = c1 = cc;
      cb = wincon[n]->s2[cc];
      if (cb) {
        wincon[n]->Ds[c] = wincon[n]->Ds[cb];
        /* 3D variables */
        while (c1 != window[n]->zm1[c1]) {
          wincon[n]->u1vh[c1] = wincon[n]->u1vh[cb];
          wincon[n]->u2vh[c1] = wincon[n]->u2vh[cb];
          wincon[n]->u1kh[c1] = wincon[n]->u1kh[cb];
          wincon[n]->u2kh[c1] = wincon[n]->u2kh[cb];
          c1 = window[n]->zm1[c1];
        }
      }
    }

    /* Set the flag for inner explicit maps. These maps are the      */
    /* inverse of cyclic boundaries, i.e. a part of the domain       */
    /* becomes infinitely small rather than infinitely large. These  */
    /* explicit maps are characterised by map[c] != c && map[map[c]  */
    /* == c.                                                         */
    for (cc = 1; cc <= window[n]->b2_e1; cc ++) {
      c = window[n]->w2_e1[cc];
      c1 = window[n]->xm1[c];
      if (c != c1 && c == window[n]->xm1[c1])
	wincon[n]->exmapf |= E1_INNER;
      c1 = window[n]->xp1[c];
      if (c != c1 && c == window[n]->xp1[c1])
	wincon[n]->exmapf |= E1_INNER;
    }
    for (cc = 1; cc <= window[n]->b2_e2; cc ++) {
      c = window[n]->w2_e2[cc];
      c1 = window[n]->ym1[c];
      if (c != c1 && c == window[n]->ym1[c1])
	wincon[n]->exmapf |= E2_INNER;
      c1 = window[n]->yp1[c];
      if (c != c1 && c == window[n]->yp1[c1])
	wincon[n]->exmapf |= E2_INNER;
    }

#if defined(HAVE_SEDIMENT_MODULE)
    /* Sediment temporary values */
    if (wincon[n]->do_sed) {
      for (cc = 1; cc <= window[n]->enonS; cc++) {
        c = window[n]->wsa[cc];
        wincon[n]->d1[cc] = 0.5 * (geom->thetau1[c] + geom->thetau2[c]);
      }
    }
#endif

    /* Set the pointers for horizontal mixing */
    hvisc_init(master, wincon);

    window[n]->wincon = wincon[n];
  }

  /* Set the ghost cell map                                          */
  for (n = 1; n <= nwindows; n++) {
    winsize = window[n]->sgsiz;
    memset(wincon[n]->gmap, 0, winsize * sizeof(int));
    for (cc = 1; cc <= window[n]->b3_t; cc++) {
      c = window[n]->w3_t[cc];
      c1 = window[n]->xm1[c];
      if (window[n]->xm1[c1] == c1) wincon[n]->gmap[c] |= L_EDGE;
      c1 = window[n]->xp1[c];
      if (window[n]->xp1[c1] == c1) wincon[n]->gmap[c] |= R_EDGE;
      c1 = window[n]->ym1[c];
      if (window[n]->ym1[c1] == c1) wincon[n]->gmap[c] |= B_EDGE;
      c1 = window[n]->yp1[c];
      if (window[n]->yp1[c1] == c1) wincon[n]->gmap[c] |= F_EDGE;
    }
    /* Copy to the ghost cells                                       */
    for (cc = 1; cc <= window[n]->nbpt; cc++) {
      c = window[n]->bpt[cc];
      c1 = window[n]->bin[cc];
      if (c == window[n]->xm1[c1])
	wincon[n]->gmap[c] |= (L_EDGE|U1SOLID);
      if (c == window[n]->xp1[c1])
	wincon[n]->gmap[c] |= (R_EDGE|U1SOLID);
      if (c == window[n]->ym1[c1])
	wincon[n]->gmap[c] |= (B_EDGE|U2SOLID);
      if (c == window[n]->yp1[c1])
	wincon[n]->gmap[c] |= (F_EDGE|U2SOLID);
    }

    /* Get the coordinates of blended zones */
    set_blend_zones(master, window[n], wincon[n]);
  }

  /*
    for (cc = 1; cc <= window[n]->nbpt; cc++) {
      c = window[n]->bpt[cc];
      c1 = window->bin[cc];
      wincon[n]->gmap[c] = c1;
      wincon[n]->gmap[c1] = -c;
    }

    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->bot_t[cc];
      c1 = window[n]->zm1[c];
      wincon[n]->gmap[c1] = c;
    }

    for (cc = 1; cc <= window[n]->ngsed; cc++) {
      c = window[n]->gsed_t[cc];
      wincon->gmap[c] = window->ised_t[cc];
    }
  }
  */

  /* Get the debug coordinate and window */
  if (master->dbc) {
    int wn = geom->fm[master->dbc].wn;
    c = geom->fm[master->dbc].sc;
    for (n = 1; n <= nwindows; n++) {
      wincon[n]->dbw = wn;
      wincon[n]->dbc = c;
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
    if (wincon->sigma) {
      d_free_1d(wincon->Ds);
      d_free_1d(wincon->Hs);
      d_free_1d(wincon->Hn1);
      d_free_1d(wincon->Hn2);
      d_free_1d(wincon->mdx);
      d_free_1d(wincon->mdy);
    }
    d_free_1d(wincon->rdens);
    d_free_1d(wincon->Cd);
    d_free_1d(wincon->u1c1);
    d_free_1d(wincon->u1c3);
    d_free_1d(wincon->u1c4);
    d_free_1d(wincon->u1c5);
    d_free_1d(wincon->u1c6);
    d_free_1d(wincon->u2c1);
    d_free_1d(wincon->u2c3);
    d_free_1d(wincon->u2c4);
    d_free_1d(wincon->u2c5);
    d_free_1d(wincon->u2c6);
    d_free_1d(wincon->topdensu1);
    d_free_1d(wincon->topdensu2);
    d_free_1d(wincon->densavu1);
    d_free_1d(wincon->densavu2);
    d_free_1d(wincon->u1adv);
    d_free_1d(wincon->u2adv);
    d_free_1d(wincon->u1inter);
    d_free_1d(wincon->u2inter);
    /* Don't free horizontal mixing variables if Smagorinsky         */
    /* diffusion is used: this is handled when tracers (sdc) are     */
    /* freed.                                                        */
    if (wincon->u1vh && !(wincon->smagcode & U1_SP))
      d_free_1d(wincon->u1vh);
    if (wincon->u2vh && !(wincon->smagcode & U2_SP))
      d_free_1d(wincon->u2vh);
    if (wincon->u1kh && !(wincon->smagcode & U1_SPK))
      d_free_1d(wincon->u1kh);
    if (wincon->u2kh && !(wincon->smagcode & U2_SPK))
      d_free_1d(wincon->u2kh);
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
void pre_run_setup(master_t *master,  /* Master data structure */
                   geometry_t **window, /* Window geometry */
                   window_t **windat, /* Window data */
                   win_priv_t **wincon  /* Data private to windows */
  )
{
  geometry_t *geom = master->geom;
  int nwindows = master->nwindows;
  int cc, c, lc, n, m;

  for (n = 1; n <= nwindows; n++) {
    /* Note: these must be direct calls, not via the pointers in the master-> */
    win_data_fill_3d(master, window[n], windat[n], master->nwindows);
    win_data_fill_2d(master, window[n], windat[n], master->nwindows);
    set_map_inside(window[n]);
    memcpy(wincon[n]->oldeta, windat[n]->eta, window[n]->sgsizS * sizeof(double));
    set_dz(window[n], windat[n], wincon[n]);
    get_depths(window[n], windat[n], wincon[n]);
    density_w(window[n], windat[n], wincon[n]);
    OBC_overlap(window[n], wincon[n]);
    if (!(master->vinit & NONE)) {
      win_data_fill_3d(master, window[n], windat[n], master->nwindows);
      set_dz_at_u1(window[n], windat[n], wincon[n]);
      set_dz_at_u2(window[n], windat[n], wincon[n]);
      if (master->vinit == VINIT_GEO)
	calc_geostrophic(window[n], windat[n], wincon[n]);
      set_flux_3d(window[n], windat[n], wincon[n], VEL3D);
      set_dz(window[n], windat[n], wincon[n]);

      /* Set a no-gradient over inside stagger OBCs */
      for (m = 0; m < window[n]->nobc; m++) {
	open_bdrys_t *open = window[n]->open[m];
	int *mi;
	if (open->stagger & INFACE) {
	  mi = (open->ocodex & L_EDGE) ? window[n]->xm1 : window[n]->xp1;
	  for (cc = 1; cc <= open->no3_e1; cc++) {
	    c = open->obc_e1[cc];
	    windat[n]->u1flux3d[mi[c]] = windat[n]->u1flux3d[c];
	  }
	  mi = (open->ocodey & B_EDGE) ? window[n]->ym1 : window[n]->yp1;
	  for (cc = 1; cc <= open->no3_e2; cc++) {
	    c = open->obc_e2[cc];
	    windat[n]->u2flux3d[mi[c]] = windat[n]->u2flux3d[c];
	  }
	}
      }
      win_data_empty_3d(master, window[n], windat[n], VELOCITY);
    }
    /* Set the grid angles for TIDALC                                */
    for (m = 0; m < window[n]->nobc; m++) {
      open_bdrys_t *open = window[n]->open[m];
      tidalc_setup(geom, window[n], open);
    }

      /*win_data_empty_3d(master,window[n],windat[n],TRACERS);
	win_data_refill_3d(master,window[n],windat[n],nwindows,TRACERS); */
  }

  /* Set the reef fractions                                          */
  set_reef_frac(master, window, windat, wincon);
  
  if (!(master->vinit & NONE)) {
    for (n = 1; n <= nwindows; n++) {
      /* Transfer velocities for w computation */
      win_data_refill_3d(master, window[n], windat[n], master->nwindows, VELOCITY);
      /* Vertical velocity computation */
      vel_w_update(window[n], windat[n], wincon[n]);
      vint_3d(window[n], windat[n], wincon[n]);
    }
  }

  master->is_filled = 0;
  master_fill(master, window, windat, wincon);
  master->is_filled = 0;

  for (n = 1; n <= nwindows; n++) {
    win_data_fill_3d(master, window[n], windat[n], master->nwindows);
    win_data_fill_2d(master, window[n], windat[n], master->nwindows);
    compute_ref_density(master, window[n], windat[n], wincon[n]);
    Set_lateral_BC_density_w(windat[n]->dens, window[n]->nbpt,
                             window[n]->bpt, window[n]->bin);

    windat[n]->dtf2 = master->dt / master->iratio;
    wincon[n]->calc_closure(window[n], windat[n], wincon[n]);

    /* Set a no-gradient over ghost OBCs */
    OBC_bgz_nograd(window[n]);

    /* Set the leapfrog arrays for the first iteration */
    for (cc = 1; cc <= window[n]->b3_e1; cc ++) {
      c = window[n]->w3_e1[cc];
      windat[n]->u1b[c] = windat[n]->u1[c];
    }
    for (cc = 1; cc <= window[n]->b2_e1; cc ++) {
      c = window[n]->w2_e1[cc];
      windat[n]->u1avb[c] = windat[n]->u1av[c];
    }
    for (cc = 1; cc <= window[n]->b3_e2; cc ++) {
      c = window[n]->w3_e2[cc];
      windat[n]->u2b[c] = windat[n]->u2[c];
    }
    for (cc = 1; cc <= window[n]->b2_e2; cc ++) {
      c = window[n]->w2_e2[cc];
      windat[n]->u2avb[c] = windat[n]->u2av[c];
    }
    for (cc = 1; cc <= window[n]->b2_t; cc ++) {
      c = window[n]->w2_t[cc];
      windat[n]->etab[c] = wincon[n]->oldeta[c] = windat[n]->eta[c];
    }

    /* Set the lateral boundary conditions for velocity.               */
#if !GLOB_BC
    vel2D_lbc(windat[n]->u1, window[n]->nbpte1, window[n]->nbe1,
	      window[n]->bpte1, window[n]->bine1, wincon[n]->slip);
    vel2D_lbc(windat[n]->u2, window[n]->nbpte2, window[n]->nbe2,
	      window[n]->bpte2, window[n]->bine2, wincon[n]->slip);
    vel2D_lbc(windat[n]->u1b, window[n]->nbpte1, window[n]->nbe1,
	      window[n]->bpte1, window[n]->bine1, wincon[n]->slip);
    vel2D_lbc(windat[n]->u2b, window[n]->nbpte2, window[n]->nbe2,
	      window[n]->bpte2, window[n]->bine2, wincon[n]->slip);
    vel2D_lbc(windat[n]->u1av, window[n]->nbpte1S, window[n]->nbe1S,
	      window[n]->bpte1S, window[n]->bine1S, wincon[n]->slip);
    vel2D_lbc(windat[n]->u2av, window[n]->nbpte2S, window[n]->nbe2S,
	      window[n]->bpte2S, window[n]->bine2S, wincon[n]->slip);
    set_lateral_bc_eta(windat[n]->etab, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
    set_lateral_bc_eta(wincon[n]->oldeta, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
#endif

    if (master->mode2d)
      mode2d_tracer_init(window[n], windat[n], wincon[n]);

    win_data_empty_3d(master, window[n], windat[n], VELOCITY);
    win_data_empty_3d(master, window[n], windat[n], WVEL);
    win_data_empty_2d(master, window[n], windat[n], DEPTH);

    if (master->runmode & TRANS) {
      set_map_l(window[n]);
      calc_tmass(window[n], windat[n], wincon[n], 
		 wincon[n]->tmass, wincon[n]->tsmass);
      /* Use the viscosity tensors to store streamline Courant nos. */
      windat[n]->origin = windat[n]->Vz;
      windat[n]->pc = wincon[n]->t11;
      windat[n]->qc = wincon[n]->t12;
      windat[n]->rc = wincon[n]->t22;
      wincon[n]->p1 = windat[n]->u1b;
      wincon[n]->p2 = windat[n]->u2b;
      /* Initialise the streamline origin for dumps */
      if (!(master->tmode & SP_ORIGIN)) {
	for (cc = 1; cc < window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  windat[n]->origin[c] = (double)c;
	}
      }
    }

    /* Initialise the sediments */
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

    /* Initialise the waves */
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

    /* Initialise the ecology */
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

  if (master->dozoom) {
    Set_lateral_BC_density_w(master->dens, geom->nbpt, geom->bpt,
                             geom->bin);
    global_interp(geom, master->dens, geom->nzin);
  }

  /* Transfer the backward arrays to auxiliary cells */
  /* Slave to master */
  for (n = 1; n <= nwindows; n++) {
    s2m_vel(master->u1avb,windat[n]->u1avb,
	    window[n]->s2m,window[n]->s2me1,window[n]->ns2mS);
    s2m_vel(master->u2avb,windat[n]->u2avb,
	    window[n]->s2m,window[n]->s2me2,window[n]->ns2mS); 
  }
  /* Master to slave */
  for (n = 1; n <= nwindows; n++) {
    int ce1, ce2;
    for (cc = 1; cc <= window[n]->nm2sS; cc++) {
      lc = window[n]->m2s[cc];
      ce1 = window[n]->m2se1[cc];
      ce2 = window[n]->m2se2[cc];
      windat[n]->u1avb[lc]=master->u1avb[ce1];
      windat[n]->u2avb[lc]=master->u2avb[ce2]; 
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
}

/* END pre_run_setup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets the reef fractions for the porus plate algorithm             */
/*-------------------------------------------------------------------*/
void set_reef_frac(master_t *master, 
		   geometry_t **window, 
		   window_t **windat, 
		   win_priv_t **wincon)
{
  geometry_t *geom = master->geom;
  int n, i, j, c, cs, cc, lc, wn;
  FILE *pf;
  char buf[MAXSTRLEN];
  char *fields[MAXSTRLEN * MAXNUMARGS];
  double d1;

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
/* Sets up the window structure to initialize the TIDALC angles      */
/*-------------------------------------------------------------------*/
void tidalc_setup(geometry_t *geom, geometry_t *window, open_bdrys_t *open)
{
  win_priv_t *wincon = window->wincon;
  int c, cc;

  if (open->bcond_nor & TIDALC || open->bcond_nor2d & TIDALC || 
      open->bcond_tan & TIDALC || open->bcond_tan2d & TIDALC) {
    if (window->costhu1 == NULL && window->sinthu1 == NULL) {
      window->costhu1 = d_alloc_1d(window->sgsizS);
      window->sinthu1 = d_alloc_1d(window->sgsizS);
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
	window->costhu1[cc] = geom->costhu1[c];
	window->sinthu1[cc] = geom->sinthu1[c];
      }
    }
    if (window->costhu2 == NULL && window->sinthu2 == NULL) {
      window->costhu2 = d_alloc_1d(window->sgsizS);
      window->sinthu2 = d_alloc_1d(window->sgsizS);
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
	window->costhu2[cc] = geom->costhu2[c];
	window->sinthu2[cc] = geom->sinthu2[c];
      }
    }
  }
}

/* END tidalc_setup()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets up the window structure to initialize the wave modeule       */
/*-------------------------------------------------------------------*/
void wave_setup(geometry_t *geom, geometry_t *window, int mode)
{
  win_priv_t *wincon = window->wincon;
  int c, cc;

  if (window->thetau1 == NULL) {
    window->thetau1 = d_alloc_1d(window->sgsizS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->thetau1[cc] = geom->thetau1[c];
      }
    }
  }
  if (window->thetau2 == NULL) {
    window->thetau2 = d_alloc_1d(window->sgsizS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->thetau2[cc] = geom->thetau2[c];
      }
    }
  }
  if (window->sinthcell == NULL) {
    window->sinthcell = d_alloc_1d(window->sgsizS);
    if (mode) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
        window->sinthcell[cc] = geom->sinthcell[c];
      }
    }
  }
  if (window->costhcell == NULL) {
    window->costhcell = d_alloc_1d(window->sgsizS);
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
void geom_free(master_t *master,  /* Master data structure */
               int mode         /* Code for arrays to free */
  )
{
  geometry_t *geom = master->geom;
  int n, tn, flg;

  if (mode & ALL) {
    /* Global maps */
#if !TR_CK
    free((global_map_t *)geom->fm);
#endif

    if (master->nwindows > 1) {
      master_free_nwin(master);
      i_free_1d(geom->yp1);
      i_free_1d(geom->ym1);
      i_free_1d(geom->xp1);
      i_free_1d(geom->xm1);
      i_free_1d(geom->zp1);
      i_free_1d(geom->zm1);
    }

    /* Cells to process vectors */
    i_free_1d(geom->w3_t);
    i_free_1d(geom->w3_e1);
    i_free_1d(geom->w3_e2);

    /* Sparse arrays */
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
    d_free_1d(geom->h2acell);
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
    d_free_1d(geom->u2x);
    d_free_1d(geom->u2y);
    d_free_1d(geom->dHde1);
    d_free_1d(geom->dHde2);
    d_free_1d(geom->h1au1);
    d_free_1d(geom->h2au1);
    d_free_1d(geom->botzu1);
    d_free_1d(geom->h2au2);
    d_free_1d(geom->h1au2);
    d_free_1d(geom->botzu2);
    if (geom->sm_e1)
      i_free_1d(geom->sm_e1);
    if (geom->sm_e2)
      i_free_1d(geom->sm_e2);
  }
  if (mode & (UNUSED | ALL)) {
    if (geom->mgc) i_free_1d(geom->mgc);
    /* Layer vectors */
    /*
       if(!master->waves)i_free_1d(geom->bot_t); i_free_1d(geom->bot_e1);
       i_free_1d(geom->bot_e2); */

#if defined(HAVE_SEDIMENT_MODULE)
    i_free_1d(geom->sed_t);
#endif

    /* Open boundary cell vectors */
    /*
    if (geom->nobc)
      s_free_2d(geom->owc);
    */
    return;
    if (geom->nobc)
      s_free_2d(geom->owc);
    for (n = 0; n < geom->nobc; n++) {
      /* Open boundary vectors */
      i_free_1d(geom->open[n]->oi1_t);
      i_free_1d(geom->open[n]->oi1_e1);
      i_free_1d(geom->open[n]->oi1_e2);
      i_free_1d(geom->open[n]->oi2_t);
      i_free_1d(geom->open[n]->oi2_e1);
      i_free_1d(geom->open[n]->oi2_e2);
      /* Some of the OBC vectors may be required for reading data */
      /* from file on the master; don't deallocate in this case.  */
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
void window_init(geometry_t *geom,  /* Sparse global geometry structure */
                 geometry_t **window  /* Window data structure */
  )
{
  int nwindows;                 /* Number of windows */
  int n, cc;                    /* Counters */
  int c, c1, cb;                /* Sparse locations */
  int i, j, k;                  /* Cartesian locations */
  int **s2;                     /* Dummy buffer for ghost locations */

  nwindows = geom->nwindows;
  c = 0;
  for (n = 1; n <= nwindows; n++)
    if (window[n]->enon > c)
      c = window[n]->enon;
  s2 = i_alloc_2d(c + 1, nwindows + 1);

  /*-----------------------------------------------------------------*/
  /* Initialise window data */
  for (n = 1; n <= nwindows; n++) {

    window[n]->topgrid = geom->topgrid;
    window[n]->totarea = geom->totarea;

    /* Set geog flag */
    window[n]->is_geog = geom->is_geog;
    
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];

      /* If a ghost cell is encountered, save it's wet neighbour in */
      /* wincon->s2 and continue.                                   */
      s2[n][cc] = 0;
      if (geom->fm[c].wn != n &&
	  i == NOTVALID && j == NOTVALID && k == NOTVALID) {

	/*
        if (geom->fm[geom->xp1[c]].wn != 0)
          s2[n][cc] = geom->fm[geom->xp1[c]].sc;
        else if (geom->fm[geom->xm1[c]].wn != 0)
          s2[n][cc] = geom->fm[geom->xm1[c]].sc;
        else if (geom->fm[geom->yp1[c]].wn != 0)
          s2[n][cc] = geom->fm[geom->yp1[c]].sc;
        else if (geom->fm[geom->ym1[c]].wn != 0)
          s2[n][cc] = geom->fm[geom->ym1[c]].sc;
        else if (geom->fm[geom->ym1[c]].wn != 0)
          s2[n][cc] = geom->fm[geom->ym1[c]].sc;
        else if (geom->fm[geom->yp1[geom->xp1[c]]].wn != 0)
          s2[n][cc] = geom->fm[geom->yp1[geom->xp1[c]]].sc;
        else if (geom->fm[geom->yp1[geom->xm1[c]]].wn != 0)
          s2[n][cc] = geom->fm[geom->yp1[geom->xm1[c]]].sc;
        else if (geom->fm[geom->ym1[geom->xp1[c]]].wn != 0)
          s2[n][cc] = geom->fm[geom->ym1[geom->xp1[c]]].sc;
        else if (geom->fm[geom->ym1[geom->xm1[c]]].wn != 0)
          s2[n][cc] = geom->fm[geom->ym1[geom->xm1[c]]].sc;
	*/
	/*

	c1 = ANY(c, geom->bpt, geom->nbpt);
	s2[n][cc] = geom->bin[c1];
	*/
	s2[n][cc] = window[n]->wsa[window[n]->wgst[cc]];

	if (s2[n][cc] == 0) {
	  if (geom->fm[geom->xp1[c]].wn != 0)
	    s2[n][cc] = geom->xp1[c];
	  else if (geom->fm[geom->xm1[c]].wn != 0)
	    s2[n][cc] = geom->xm1[c];
	  else if (geom->fm[geom->yp1[c]].wn != 0)
	    s2[n][cc] = geom->yp1[c];
	  else if (geom->fm[geom->ym1[c]].wn != 0)
	    s2[n][cc] = geom->ym1[c];
	  else if (geom->fm[geom->yp1[geom->xp1[c]]].wn != 0)
	    s2[n][cc] = geom->yp1[geom->xp1[c]];
	  else if (geom->fm[geom->yp1[geom->xm1[c]]].wn != 0)
	    s2[n][cc] = geom->yp1[geom->xm1[c]];
	  else if (geom->fm[geom->ym1[geom->xp1[c]]].wn != 0)
	    s2[n][cc] = geom->ym1[geom->xp1[c]];
	  else if (geom->fm[geom->ym1[geom->xm1[c]]].wn != 0)
	    s2[n][cc] = geom->ym1[geom->xm1[c]];
	}
        continue;
      }

      /* Cell centered 3D arrays */
      window[n]->gridz[cc] = geom->gridz[c];
      window[n]->cellz[cc] = geom->cellz[c];

      if (cc <= window[n]->enonS) {

        /* Cell centered 2D arrays */
        window[n]->h1acell[cc] = geom->h1acell[c];
        window[n]->h2acell[cc] = geom->h2acell[c];
        window[n]->cellarea[cc] = geom->cellarea[c];
        window[n]->botz[cc] = geom->botz[c];
        window[n]->cellx[cc] = geom->cellx[c];
        window[n]->celly[cc] = geom->celly[c];
        window[n]->dHde1[cc] = geom->dHde1[c];
        window[n]->dHde2[cc] = geom->dHde2[c];

        /* e1 face 2D centered arrays */
        window[n]->h1au1[cc] = geom->h1au1[c];
        window[n]->h2au1[cc] = geom->h2au1[c];
        window[n]->u1x[cc] = geom->u1x[c];
        window[n]->u1y[cc] = geom->u1y[c];
        window[n]->botzu1[cc] = geom->botzu1[c];

        /* e2 face 2D centered arrays */
        window[n]->h2au2[cc] = geom->h2au2[c];
        window[n]->h1au2[cc] = geom->h1au2[c];
        window[n]->u2x[cc] = geom->u2x[c];
        window[n]->u2y[cc] = geom->u2y[c];
        window[n]->botzu2[cc] = geom->botzu2[c];

        /* At grid corners */
        window[n]->gridx[cc] = geom->gridx[c];
        window[n]->gridy[cc] = geom->gridy[c];
        window[n]->botzgrid[cc] = geom->botzgrid[c];

        /* Sediments */
        if (geom->sednz > 0) {
          int kk;
          for (kk = 0; kk < geom->sednz + 1; kk++)
            window[n]->gridz_sed[kk][cc] = geom->gridz_sed[kk][c];
          for (kk = 0; kk < geom->sednz; kk++)
            window[n]->cellz_sed[kk][cc] = geom->cellz_sed[kk][c];
        }
      }
    }
    if (DEBUG("init_w"))
      dlog("init_w", "Window geometry %d initialised OK\n", n);
  }

  /*-----------------------------------------------------------------*/
  /* Set values in ghost locations */
  if (nwindows > 1) {
    if (master->dozoom & PRECOND)
      hgrid_zoom_p(geom, window);

  for (n = 1; n <= nwindows; n++) {

    /* Adjust the horizontal grid spacings in zoomed grids */
    if (!(master->dozoom & PRECOND))
      hgrid_zoom(geom, window[n]);

    for (cc = 1; cc <= window[n]->enonS; cc++) {
      c = c1 = cc;
      cb = s2[n][cc];

      if (cb) {
	/*
        window[n]->h1au1[c] = window[n]->h1au1[cb];
        window[n]->h2au2[c] = window[n]->h2au2[cb];
        window[n]->h1acell[c] = window[n]->h1acell[cb];
        window[n]->h2acell[c] = window[n]->h2acell[cb];
        window[n]->h2au1[c] = window[n]->h2au1[cb];
        window[n]->h1au2[c] = window[n]->h1au2[cb];
        window[n]->cellarea[c] = window[n]->cellarea[cb];
	*/
        window[n]->h1au1[c] = geom->h1au1[cb];
        window[n]->h2au2[c] = geom->h2au2[cb];
        window[n]->h1acell[c] = geom->h1acell[cb];
        window[n]->h2acell[c] = geom->h2acell[cb];
        window[n]->h2au1[c] = geom->h2au1[cb];
        window[n]->h1au2[c] = geom->h1au2[cb];
        window[n]->cellarea[c] = geom->cellarea[cb];
        window[n]->cellx[c] = geom->cellx[window[n]->wsa[cc]];
        window[n]->celly[c] = geom->celly[window[n]->wsa[cc]];
        window[n]->gridx[c] = geom->gridx[window[n]->wsa[cc]];
        window[n]->gridy[c] = geom->gridy[window[n]->wsa[cc]];
        window[n]->botzgrid[c] = geom->botzgrid[window[n]->wsa[c]];
      }
    }

    for (i = 0; i < window[n]->nobc; i++) {
      open_bdrys_t *open = window[n]->open[i];
      int *obc = NULL, *oi1 = NULL, nc = 0;
      if (open->ocodex & R_EDGE) {
	nc = open->no2_e1;
	obc = open->obc_e1;
	oi1 = open->oi1_e1;
      } else if (open->ocodey & F_EDGE) {
	nc = open->no2_e2;
	obc = open->obc_e2;
	oi1 = open->oi1_e2;
      }
      for(cc = 1; cc <= nc; cc++) {
	c = obc[cc];
	cb = oi1[cc];
	window[n]->h1au1[c] = window[n]->h1au1[cb];
	window[n]->h2au2[c] = window[n]->h2au2[cb];
	window[n]->h1acell[c] = window[n]->h1acell[cb];
	window[n]->h2acell[c] = window[n]->h2acell[cb];
	window[n]->h2au1[c] = window[n]->h2au1[cb];
	window[n]->h1au2[c] = window[n]->h1au2[cb];
	window[n]->cellarea[c] = window[n]->cellarea[cb];
      }
      /* Metric terms for OBC ghost cells */
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	cb = window[n]->wsa[c];
	for (j = 0; j < open->bgz; j++) {
	  c = open->omap[c];
	  cb = geom->open[open->mwn]->omap[cb];
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
	  p1 = open->nmap[c];
	  c = open->omap[c];
	  if (open->type == U1BDRY) {
	    if (open->ocodex & R_EDGE)
	      window[n]->h1au1[c] += 0.5 * (window[n]->h1acell[c] + open->rlen - window[n]->h1au1[c]);
	    else
	      window[n]->h1au1[p1] += 0.5 * (window[n]->h1acell[c] + open->rlen - window[n]->h1au1[c]);
	    window[n]->cellarea[c] = open->rlen * window[n]->h2acell[c];
	    window[n]->h1acell[c] += open->rlen;
	  }
	  if (open->type == U2BDRY) {
	    if (open->ocodey & F_EDGE)
	      window[n]->h2au2[c] += 0.5 * (window[n]->h2acell[c] + open->rlen - window[n]->h2au2[c]);
	    else
	      window[n]->h2au2[p1] += 0.5 * (window[n]->h2acell[c] + open->rlen - window[n]->h2au2[c]);
	    window[n]->cellarea[c] = open->rlen * window[n]->h1acell[c];
	    window[n]->h2acell[c] += open->rlen;
	  }
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
void win_geom_clear(geometry_t **window,  /* Window data structure */
                    int nwindows  /* Number of windows */
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
    i_free_1d(window[n]->m2se2);
    i_free_1d(window[n]->s2m);
    i_free_1d(window[n]->s2me1);
    i_free_1d(window[n]->s2me2);
    i_free_1d(window[n]->sur_t);
    i_free_1d(window[n]->nsur_t);
    i_free_1d(window[n]->bot_t);
    i_free_1d(window[n]->sur_e1);
    i_free_1d(window[n]->bot_e1);
    i_free_1d(window[n]->sur_e2);
    i_free_1d(window[n]->bot_e2);
    i_free_1d(window[n]->m2d);
    i_free_1d(window[n]->xp1);
    i_free_1d(window[n]->xm1);
    i_free_1d(window[n]->yp1);
    i_free_1d(window[n]->ym1);
    i_free_1d(window[n]->zp1);
    i_free_1d(window[n]->zm1);
    i_free_1d(window[n]->xpym1);
    i_free_1d(window[n]->xmyp1);
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
    if (window[n]->zbe1)
      i_free_1d(window[n]->zbe1);
    if (window[n]->zbe2)
      i_free_1d(window[n]->zbe2);
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

    /* Open boundary cell vectors */
    for (nn = 0; nn < window[n]->nobc; nn++) {
      i_free_1d(window[n]->open[nn]->obc_t);
      i_free_1d(window[n]->open[nn]->bot_t);
      i_free_1d(window[n]->open[nn]->obc_e1);
      i_free_1d(window[n]->open[nn]->obc_e2);
      i_free_1d(window[n]->open[nn]->oi1_t);
      i_free_1d(window[n]->open[nn]->oi1_e1);
      i_free_1d(window[n]->open[nn]->oi1_e2);
      i_free_1d(window[n]->open[nn]->oi2_t);
      i_free_1d(window[n]->open[nn]->oi2_e1);
      i_free_1d(window[n]->open[nn]->oi2_e2);
      if (window[n]->open[nn]->cyc_t)
        i_free_1d(window[n]->open[nn]->cyc_t);
      if (window[n]->open[nn]->cyc_e1)
        i_free_1d(window[n]->open[nn]->cyc_e1);
      if (window[n]->open[nn]->cyc_e2)
        i_free_1d(window[n]->open[nn]->cyc_e2);
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
  master->master_fill(master, window, windat, wincon);

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
void aux_maps(geometry_t *window, /* Processing window */
              int cc            /* Local sparse coordinate */
  )
{
  int c;                        /* Global sparse coordinate */
  int cg;                       /* Mapped global sparse counters */
  int cl;                       /* Mapped local sparse coordinate */
  int zc;

  c = window->wsa[cc];

  if (window->xp1[cc] == cc) {
    cg = geom->xp1[c];
    for (zc = 1; zc < window->zmfe1; zc++)
      cg = geom->xp1[cg];
    cl = geom->fm[cg].ac;
    if (cl)
      window->xp1[cc] = cl;
  }
  if (window->xm1[cc] == cc) {
    cg = geom->xm1[c];
    for (zc = 1; zc < window->zmfe1; zc++)
      cg = geom->xm1[cg];
    cl = geom->fm[cg].ac;
    if (cl)
      window->xm1[cc] = cl;
  }
  if (window->yp1[cc] == cc) {
    cg = geom->yp1[c];
    for (zc = 1; zc < window->zmfe2; zc++)
      cg = geom->yp1[cg];
    cl = geom->fm[cg].ac;
    if (cl)
      window->yp1[cc] = cl;
  }
  if (window->ym1[cc] == cc) {
    cg = geom->ym1[c];
    for (zc = 1; zc < window->zmfe2; zc++)
      cg = geom->ym1[cg];
    cl = geom->fm[cg].ac;
    if (cl)
      window->ym1[cc] = cl;
  }

  /* Generate ghost cell maps in directions that have not already */
  /* been defined.  */
  if (!geom->fm[c].wn) {        /* Window = 0 => ghost cell */
    if (window->xp1[cc] == cc) {  /* Ghost cell is self pointing */
      cg = geom->xp1[c];        /* Global cell to the east */
      /* for(zc=1; zc<zs; zc++)cg=geom->xp1[cg]; */
      cl = geom->fm[cg].ac;     /* Auxiliary cell to the east of cc */
      if (!cl)
        cl = geom->fm[cg].sc;   /* Wet cell if aux. cell undefined */
      if (cl)
        window->xp1[cc] = cl;
    }
    if (window->xm1[cc] == cc) {
      cg = geom->xm1[c];
      /* for(zc=1; zc<zs; zc++)cg=geom->xm1[cg]; */
      cl = geom->fm[cg].ac;
      if (!cl)
        cl = geom->fm[cg].sc;
      if (cl)
        window->xm1[cc] = cl;
    }
    if (window->yp1[cc] == cc) {
      cg = geom->yp1[c];
      /* for(zc=1; zc<zs; zc++)cg=geom->yp1[cg]; */
      cl = geom->fm[cg].ac;
      if (!cl)
        cl = geom->fm[cg].sc;
      if (cl)
        window->yp1[cc] = cl;
    }
    if (window->ym1[cc] == cc) {
      cg = geom->ym1[c];
      /* for(zc=1; zc<zs; zc++)cg=geom->ym1[cg]; */
      cl = geom->fm[cg].ac;
      if (!cl)
        cl = geom->fm[cg].sc;
      if (cl)
        window->ym1[cc] = cl;
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
void fill_multidt_aux(master_t *master, /* Master data structure */
                      geometry_t *window, /* Window structure */
                      window_t *windat  /* Window data structure */
  )
{
  int c, cc;                    /* Local sparse coordinate / counter */
  int tn;                       /* Tracer counter */

  /* Auxiliary multi-dt cells */
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
void fill_multidt_aux_2d(geometry_t *window,  /* Window structure */
                         window_t *windat,  /* Window data structure */
                         double *vel, /* Velocity array */
                         double *ae /* Multi-dt cells velocity time t+1
                                       values */
  )
{
  int c, cc;                    /* Local sparse coordinate / counter */

  /* Auxiliary multi-dt cells */
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
void get_timesteps(geometry_t **window, /* Processing window */
                   window_t **windat, /* Window data structure */
                   win_priv_t **wincon, /* Window geometry / constants */
                   int nwindows,  /* Number of windows */
                   master_t *master /* Model data structure */
  )
{
  int wn, n, *nn;               /* Window counters */
  int c, c2, cc;                /* Sparse cponters */
  double top;                   /* Surface level */
  double d1, *d2, d3, *d4;      /* Dummies */
  double maxvel = 0.5;          /* Maximum 3D velocity expected */
  double int_wave_speed = 2.0;  /* Maximum internal wave speed */
  double sf = 0.9;              /* Safety factor */

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
    /* If the time-step is already defined (non-zero) then use this */
    /* time-step for all windows.  */
    if (windat[n]->dt) {
      for (c = 1; c <= nwindows; c++)
        nn[c] = c;
      d2[n] = windat[n]->dt / windat[n]->iratio;
      d4[n] = windat[n]->dt;
    }
    /*---------------------------------------------------------------*/
    /* Get the time-steps in each window from the CFL condition */
    else {
      d2[n] = d4[n] = 1e10;

      nn[n] = n;
      for (cc = 1; cc <= window[n]->v3_t; cc++) {
        c = window[n]->w3_t[cc];
        c2 = window[n]->m2d[c];
        d1 = 1.0 / (window[n]->h1acell[c2] * window[n]->h1acell[c2]) +
          1.0 / (window[n]->h2acell[c2] * window[n]->h2acell[c2]);
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
void timeaux(geometry_t **window, /* Processing window */
             window_t **windat, /* Window data structure */
             win_priv_t **wincon, /* Window geometry / constants */
             int nwindows       /* Number of windows */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int n;                        /* Window counter */
  int gc;                       /* Global sparse coordinate */
  int ac;                       /* Auxiliary cell counter */
  int wn;                       /* Window corresponding to gc */

  /* Loop over all windows */
  for (n = 1; n <= nwindows; n++)
    window[n]->naux_t = 0;

  if (!(wincon[1]->stab & MULTI_DT))
    return;

  for (n = 1; n <= nwindows; n++) {
    window[n]->aux_t = NULL;
    ac = 0;

    /*---------------------------------------------------------------*/
    /* Count the number of auxiliary cells in window n which are */
    /* associated with a window having a longer time-step than */
    /* window n.  */
    for (cc = window[n]->v3_t + 1; cc <= window[n]->n3_t; cc++) {
      c = window[n]->w3_t[cc];  /* Local auxiliary/ghost location */
      gc = window[n]->wsa[c];   /* Global location */
      wn = geom->fm[gc].wn;     /* Window corresponding to gc */
      /* Note : ghost cells have a corresponding window = 0.  */
      if (wn && wn != n && windat[wn]->dts > windat[n]->dts)
        ac++;
    }
    if (DEBUG("init_w") && ac)
      dlog("init_w", "Window %d = %d multi time-step auxiliary cells\n", n,
           ac);

    /*---------------------------------------------------------------*/
    /* Allocate meory for long time-step auxiliary cells */
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
    /* sparse coordinates in window n and save the longer time-step */
    /* corresponding to the auxiliary cell window to array.  */
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
void prints(double *A,          /* Array to print */
            char *fname,        /* File name */
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
			 geometry_t *geom,  /* Global geometry       */
			 geometry_t *window /* Processing window     */
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
/*-------------------------------------------------------------------*/
void get_inner_exmap(geometry_t *geom, /* Processing window          */
		     geometry_t *window /* Processing window         */
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
void nan_check(geometry_t **window, /* Processing window */
	       window_t **windat, /* Window data structure */
	       win_priv_t **wincon, /* Window geometry / constants */
	       int nwindows       /* Number of windows */
	       )
{
  int n, c, cc;

  for (n = 1; n <= nwindows; n++) {
    /* eta and atmospheric pressure */
    for (cc = 1; cc <= window[n]->b2_t; cc++) {
      c = window[n]->w2_t[cc];
      if (isnan(windat[n]->eta[c]))
	hd_warn("NaN found in eta, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
      if (isnan(windat[n]->patm[c]))
	hd_warn("NaN found in patm, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
    }
    /* Wind and 2D velocity */
    for (cc = 1; cc <= window[n]->b2_e1; cc++) {
      c = window[n]->w2_e1[cc];
      if (isnan(windat[n]->u1av[c]))
	hd_warn("NaN found in u1av, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
      if (isnan(windat[n]->wind1[c]))
	hd_warn("NaN found in wind1, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
    }
    for (cc = 1; cc <= window[n]->b2_e2; cc++) {
      c = window[n]->w2_e2[cc];
      if (isnan(windat[n]->u2av[c]))
	hd_warn("NaN found in u2av, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
      if (isnan(windat[n]->wind2[c]))
	hd_warn("NaN found in wind2, window %d at (%d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c]);
    }
    /* T/S */
    for (cc = 1; cc <= window[n]->b3_t; cc++) {
      c = window[n]->w3_t[cc];
      if (isnan(windat[n]->temp[c]))
	hd_warn("NaN found in temp, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
      if (isnan(windat[n]->sal[c]))
	hd_warn("NaN found in salt, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    /* 3D velocity */
    for (cc = 1; cc <= window[n]->b3_e1; cc++) {
      c = window[n]->w3_e1[cc];
      if (isnan(windat[n]->u1[c]))
	hd_warn("NaN found in u1, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    for (cc = 1; cc <= window[n]->b3_e2; cc++) {
      c = window[n]->w3_e2[cc];
      if (isnan(windat[n]->u2[c]))
	hd_warn("NaN found in u2, window %d at (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    for (cc = window[n]->b3_e1 + 1; cc <= window[n]->a3_e1; cc++) {
      c = window[n]->w3_e1[cc];
      if (isnan(windat[n]->u1[c]))
	hd_warn("NaN found in u1, window %d at auxiliary (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
    for (cc = window[n]->b3_e2 + 1; cc <= window[n]->a3_e2; cc++) {
      c = window[n]->w3_e2[cc];
      if (isnan(windat[n]->u2[c]))
	hd_warn("NaN found in u2, window %d at auxiliary (%d %d %d)\n", n, window[n]->s2i[c], window[n]->s2j[c], window[n]->s2k[c]);
    }
  }
}

/* END nan_check()                                                   */
/*-------------------------------------------------------------------*/
