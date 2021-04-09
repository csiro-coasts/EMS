/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/tracers/tracers.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracers.c 6749 2021-03-30 00:46:08Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

#if defined(HAVE_OMP)
#include <omp.h>
#endif

/*-------------------------------------------------------------------*/
/* Routines contained in this module                                 */
/* model_step()    : Updates all prognostic variables in each window */
/* tracer_step()   : Updates tracer concentrations in each window    */
/* advect_diffuse(): Performs advection & horizontal diffusion       */
/* advect()        : Calculates the advective fluxes                 */
/* hor_diffuse()   : Calculates the horizontal diffusive fluxes      */
/* vert_diffuse_3d(): Calculates the vertical diffusion (implicit)   */
/* implicit_vdiff_tr() : Implicit vertical diffusion scheme          */
/* order1()        : Calculates first order upwind fluxes            */
/* order2()        : Calculates second order fluxes                  */
/* order4()        : Calculates fourth order fluxes                  */
/* quickest()      : Calculates QUICKEST fluxes, non-uniform grid    */
/* van_leer()      : Calculates Van Leer higher order upwind fluxes  */
/* ffsl()          : Calculates flux-form semi-lagrangian fluxes     */
/* surf_conc()     : Calculates surface tracer concentration         */
/* ultimate_filter() : Invokes the ULTIMATE limiter                  */
/* Diffusion()     : Solves the diffusion equation in each window    */
/* bdry_tracer()   : Sets open boundary conditions for each window   */
/* set_OBC_tracer()     : Invokes the OBC's                          */
/* set_lateral_BC_tr()  : Sets tracer BC's across solid boundaries   */
/* set_lateral_BC_vel() : Sets velocity BC's across solid boundaries */
/* set_dz()        : Sets the cell thickness and surface coordinates */
/* set_map_t()     : Sets spatial maps across tracer OBC's           */
/* reset_map_t()   : Resets spatial maps across tracer OBC's         */
/* perc()          : Calculates percentile values                    */
/* sorts()         : Sorts two numbers                               */
/* orders()        : Orders numbers into increasing order            */
/* eco_step()      : Performs biochem and sediment transport         */
/*-------------------------------------------------------------------*/

#define TINY  1e-20             /* Very small value                  */
#define C13   1.0 / 3.0
#define C23   2.0 / 3.0
void set_thin_tr(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_mixed_layer(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void calc_flux(geometry_t *window, window_t *windat, win_priv_t *wincon, 
	       double *fluxe1, double *fluxw, double dt, int dir1, int dir2);
void calc_age(geometry_t *window, window_t *windat, win_priv_t *wincon);
void auxiliary_routines(geometry_t *window, window_t *windat,
			win_priv_t *wincon);
void reset_dz(geometry_t *window, window_t *windat, win_priv_t *wincon);
int advect_diffuse_split(geometry_t *window, window_t *windat,
			 win_priv_t *wincon);
double check_mass(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, char *trname, int cin);
double check_tmass(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, char *trname);
void check_bounds(char *text, geometry_t *window, window_t *windat, 
		  win_priv_t *wincon);
void check_fluxes(int n, char *trname, char * text, double dtu, 
		  geometry_t *window,
		  window_t *windat, win_priv_t *wincon);
void limit_min(geometry_t *window, window_t *windat, win_priv_t *wincon);
void do_ts_relax(geometry_t *window, window_t *windat, win_priv_t *wincon);
void do_ts_increment(geometry_t *window, window_t *windat, 
		     win_priv_t *wincon);
void reset_tr_OBCflux(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      double *dtracer, double *Fx, double dt, int tn);
void swr_assimilation(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      int n);
void save_OBC_tr(geometry_t *window, window_t *windat, win_priv_t *wincon,
		 double *Fx, double *tr, int n, int mode);
void set_lateral_bdry_tr(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_bdry_cellno(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_flow(geometry_t *window, window_t *windat, win_priv_t *wincon, int mode);
double monoco(int c, int cm, int cp, double *tr, double *trmod, double sc);
void nogr(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr);
double get_upvel(double *vel, double sc, int e, int ei, int mode);
void fct_limit(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr, 
	       int ksf, int kef, int kefs, double dtu, double *osubeta, double *subeta);
double get_deriv2(geometry_t *window, window_t *windat, win_priv_t *wincon,
		  double *tr, int e, int j);
void add_cells_to_lsq(geometry_t *window, window_t *windat, win_priv_t *wincon, int ncells,
		      int *cells, int co);
void linear_limit_a(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr);
double runge_kutta(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		   double ovol, double nvol, int tc, int c, int slf);
void test_rk(int rkstage);

double rkweights[18][5] =
  { {1.00, 1.00, 1.00, 0.0, 1.00},
    {0.00, 0.00, 0.00, 0.0, 0.00},
    {1.00, 1.00, 1.00, 0.0, 0.391752226571890},
    {0.00, 0.50, 0.75, 0.0, 0.444370493651235},
    {0.00, 0.50, 0.25, 0.0, 0.555629506348765},
    {0.00, 0.50, 0.25, 0.0, 0.368410593050371},
    {0.00, 0.00, C13,  0.0, 0.620101851488403},
    {0.00, 0.00, C23,  0.0, 0.379898148511597},
    {0.00, 0.00, C23,  0.0, 0.251891774271694},
    {0.00, 0.00, 0.00, 0.0, 0.178079954393132},
    {0.00, 0.00, 0.00, 0.0, 0.821920045606868},
    {0.00, 0.00, 0.00, 0.0, 0.544974750228521},
    {0.00, 0.00, 0.00, 0.0, 0.517231671970585},
    {0.00, 0.00, 0.00, 0.0, 0.096059710526147},
    {0.00, 0.00, 0.00, 0.0, 0.063692468666290},
    {0.00, 0.00, 0.00, 0.0, 0.386708617503269},
    {0.00, 0.00, 0.00, 0.0, 0.226007483236906},
    {0.00, 0.00, 0.00, 0.0, 0.127384937332581}
  };

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Updates the tracer in each window                                 */
/*-------------------------------------------------------------------*/
void tracer_step_window(master_t *master,
                        geometry_t *window,
                        window_t *windat, 
			win_priv_t *wincon)
{

  /*-----------------------------------------------------------------*/
  /* Refill variables required for FFSL                              */
  if (wincon->trasc & FFSL)
    win_data_refill_3d(master, window, windat, master->nwindows, WVEL);

  /*-----------------------------------------------------------------*/
  /* Reset diagnostic tracers to zero if required                    */
  tr_diag_reset_w(window, windat, wincon);
  if (wincon->trtend >= 0)
    memcpy(wincon->tendency, windat->tr_wc[wincon->trtend], 
	   window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Update the tracer concentrations. Note : density and vertical   */
  /* mixing calculations are also included here.                     */
  if (!master->mode2d)
    tracer_step_3d(window, windat, wincon);
  else
    tracer_step_2d(window, windat, wincon);
  if (master->crf == RS_RESTART) return;

  /*-----------------------------------------------------------------*/
  /* Refill the master with tracer and density data from the         */
  /* window data structure.  */
  win_data_empty_3d(master, window, windat, TRACERS);

}

/* END tracer_step_window()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to solve the tracer equation over all windows             */
/*-------------------------------------------------------------------*/
void tracer_step(master_t *master,    /* Master data                 */
                 geometry_t **window, /* Window geometry             */
                 window_t **windat,   /* Window data                 */
                 win_priv_t **wincon, /* Window constants            */
                 int nwindows         /* Number windows to process   */
  )
{
  int n, nn, tn;                /* Counter                           */

  /*-----------------------------------------------------------------*/
  /* Solve the tracer equation in each window                        */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.              */
#if !GLOB_BC
    set_lateral_BC_tr(windat[n]->tr_wc, windat[n]->ntr, window[n]->nbpt,
                      window[n]->bpt, window[n]->bin);
    set_lateral_BC_vel(master->u1flux3d, window[n]->nbpte1,
                       window[n]->bpte1, window[n]->bine1);
#endif
    set_lateral_bdry_tr(window[n], windat[n], wincon[n]);

   /*---------------------------------------------------------------*/
    /* Set auxiliary cells that have been updated to the multi dt   */
    /* buffers.                                                     */
    fill_multidt_aux(master, window[n], windat[n]);
  }

  /* Invoke distributed processing step                              */
  dp_tracer_step();
  if (master->crf == RS_RESTART) return;

  /*-----------------------------------------------------------------*/
  /* Transfer mass diagnostics. This is done non-threaded, since if  */
  /* the transfers are performed within a thread, occasionally if    */
  /* multiple windows transfer to the global array at the same time, */
  /* and the global array variable sums over windows (as is the case */
  /* with total mass), then the global sum only adds mass from one   */
  /* window. This seems to be an issue with the threading libraries. */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];
    win_data_empty_3d(master, window[n], windat[n], TMASS);
  }

  /*-----------------------------------------------------------------*/
  /* Get the flushing time if required                               */
  if (master->trflsh)
    calc_flushing(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Print the region budgets if required                            */
  dump_regions(master);
}

/* END tracer_step()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialise the tracer step                                        */
/*-------------------------------------------------------------------*/
void tracer_step_init(master_t *master)
{
  int n;

  return;
  for (n = 0; n < master->ntr; n++) {
    if ((!master->advect[n]) && (master->nonlinear))
      hd_warn
        ("tracerstep: non-advecting tracer '%s' may not\nfollow water in run with non-linear eta\n",
         master->trname[n]);
  }
}

/* END tracer_step_init()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to update the tracer concentrations in each window        */
/*-------------------------------------------------------------------*/
void tracer_step_3d(geometry_t *window, /* Window geometry           */
		    window_t   *windat, /* Window data               */
		    win_priv_t *wincon  /* Window constants          */
  )
{
  double clock = dp_clock();

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across tracer OBC's, unless an  */
  /* OBC ghost zone is speified.                                     */
  set_map_t(window);
  debug_c(window, D_TS, D_INIT);

  /*-----------------------------------------------------------------*/
  /* Increment the mean counter if required                          */
  reset_means(window, windat, wincon, WIND);

  /*-----------------------------------------------------------------*/
  /* Do the advection and horizontal diffusion                       */
  TIMING_SET;
  get_bdry_cellno(window, windat, wincon);
  if (wincon->trsplit) {
    if (advect_diffuse_split(window, windat, wincon)) return;
  } else if (wincon->trasc == LAGRANGE) {
     advect_diffuse_lag(window, windat, wincon);
  } else {
    if (advect_diffuse(window, windat, wincon)) return;
  }
  TIMING_DUMP_WIN(2, "  advect_diffuse", window->wn);
  debug_c(window, D_TS, D_ADVECT);
  if (wincon->trtend >= 0)
    get_tend(window, window->w3_t, window->b3_t, windat->tr_wc[wincon->trtend],
             wincon->tendency, windat->tr_adv);

  /* Reset the maps to map to velocity cells on outward edges        */
  reset_map_t_all(window);

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion                                       */
  TIMING_SET;
  vert_diffuse_3d(window, windat, wincon);
  debug_c(window, D_TS, D_HDIFF);
  if (wincon->trtend >= 0)
    get_tend(window, window->w3_t, window->b3_t, windat->tr_wc[wincon->trtend],
             wincon->tendency, windat->tr_vdif);
  TIMING_DUMP_WIN(2, "  vert_diffuse", window->wn);

  /*-----------------------------------------------------------------*/
  /* If FFSL is used with mean velocity transport, revert arrays     */
  if (wincon->means & TRANSPORT)
    ffsl_trans_post(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the diagnostics and interfaced library routines              */
  TIMING_SET;
  auxiliary_routines(window, windat, wincon);
  TIMING_DUMP_WIN(2, "  aux_routines", window->wn);

  /*-----------------------------------------------------------------*/
  /* Set the open boundary conditions                                */
  bdry_tracer(window, windat, wincon);
  debug_c(window, D_TS, D_BDRY);

  /*-----------------------------------------------------------------*/
  /* Set the density                                                 */
  if (wincon->calc_dens)
    density_w(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the means for temperature and salinity and update the */
  /* mean counter.                                                   */
  reset_means(window, windat, wincon, RESET|TS);

  /*-----------------------------------------------------------------*/
  /* Calculate the means for tracers other than velocity             */
  reset_means(window, windat, wincon, ALL);

  /*-----------------------------------------------------------------*/
  /* Check for fatal instabilities (NaN)                             */
  if (check_unstable(window, windat, wincon, TS)) return;
  if (wincon->trtend >= 0)
    get_tend(window, window->w3_t, window->b3_t, windat->tr_wc[wincon->trtend],
             wincon->tendency, windat->tr_ncon);

  debug_c(window, D_TS, D_POST);

  windat->wclk += (dp_clock() - clock);
}

/* END tracer_step_3d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to update the tracer concentrations in each window for    */
/* the 2D mode.                                                      */
/*-------------------------------------------------------------------*/
void tracer_step_2d(geometry_t *window, /* Window geometry           */
                    window_t *windat,   /* Window data               */
                    win_priv_t *wincon  /* Window constants          */
  )
{
  int e, ee, e2;                /* Edge coordinate, counter          */
  int c, cc, c2;                /* Sparse coordinate / counter       */

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across tracer OBC's             */
  reset_map_t_all(window);

  /*-----------------------------------------------------------------*/
  /* Reset the fluxes over wet and ghost cells; substepping in the   */
  /* tracer routine may use a different timestep to windat->dt.      */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    windat->u1flux[e] /= windat->dt;
  }
  /* Set the e1 velocity                                             */
  for (ee = 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    e2 = window->m2de[e];
    windat->u1[e] = windat->u1av[e2];
    windat->u1flux3d[e] = windat->u1flux[e2];
  }

  /* Set the vertical velocity                                       */
  memset(windat->w, 0, window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Set up the cell centered dz                                     */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    wincon->dz[c] = (windat->topz[c] - window->botz[c]) * wincon->Ds[c];
  }

  /*-----------------------------------------------------------------*/
  /* Do the advection and horizontal diffusion                       */
  if (advect_diffuse_2d(window, windat, wincon)) return;

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion                                       */
  vert_diffuse_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Do the diagnostics and interfaced library routines              */
  auxiliary_routines(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the open boundary conditions                                */
  bdry_tracer(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the density                                                 */
  if (wincon->calc_dens)
    density_w(window, windat, wincon);
}

/* END tracer_step_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion                                */
/*-------------------------------------------------------------------*/
int advect_diffuse(geometry_t *window,  /* Window geometry           */
		   window_t *windat,    /* Window data               */
		   win_priv_t *wincon   /* Window constants          */
  )
{
  int nn;                       /* Tracer index                      */
  int c, cc;                    /* Sparse indices                    */
  int c2;                       /* 2D cell corresponding to 3D cell  */
  int cb;                       /* Bottom sparse coordinate (cell centre) */
  int vc;                       /* Tracer cells to process counter   */
  int vcs;                      /* Surface tracer cells to process counter */
  int vc1;                      /* As for vc excluding OBC cells     */
  int vcs1;                     /* As for vcs excluding OBC cells    */
  int ee, e, e1, e2;            /* Edge coordinate, counter          */
  int j, j1;
  int zp1;                      /* k index at k+1                    */
  int zm1;                      /* k index at k-1                    */
  int ntbdy;                    /* Number of tracers to advect       */
  int *tbdy;                    /* Tracers to advect                 */
  double top;                   /* Flux through the surface cell     */
  double surftr;                /* Surface tracer concentration      */
  double dtu;                   /* Sub-time step to use              */
  double dtm;                   /* Minimum allowable time step       */
  double trem;                  /* Time remaining                    */
  double *osubeta;              /* Old surface elevation at sub-timestep */
  double *subeta;               /* Surface elevation at sub-timestep */
  double *msubeta;              /* Mid surface elevation for FCT     */
  double d1, d2, d3, d4;        /* Dummy variables                   */
  int ksf;                      /* Surface layer for sigma / Cartesian */
  double sf = 0.8;              /* Fraction of sub-timestep to use   */
  double minval = 1e-10;        /* Minimum value for velocity        */
  int slf = 1;                  /* Set to zero on the first sub-step */
  int itermax = 20;             /* Maximum number of substeps        */
  int courf;                    /* Flag to calculate Courant numbers */
  int ii = 0, jj = 0, kk = 0;   /* (i,j,k) location of sub-step violation */
  double vm = 0.0, vs = 0.0;    /* Sub-step Courant number & grid spacing */
  char vt[4];                   /* Name of sub-step velocity component */
  double *u1, *u, *v, *w;
  double kzlim;
  int cdb = 0;
  int tc;
  double ovol, nvol;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* w1 = dtracer                                                    */
  /* w2 = Fx                                                         */
  /* w4 = Fz                                                         */
  /* w5 = Mean vertical cell spacing                                 */
  /* w6 = x direction Courant number                                 */
  /* w7 = y direction Courant number                                 */
  /* w8 = z direction Courant number                                 */
  /* w9 = vertical cell spacing between grid centers                 */
  /* d1 = osubeta                                                    */
  /* d2 = subeta                                                     */
  /* d3 = mean cell spacing in the x direction                       */
  /* d4 = mean cell spacing in the y direction                       */
  /* s1 = wet cells to process for tracers                           */

  /*-----------------------------------------------------------------*/
  /* Assign pointers                                                 */
  osubeta = wincon->d1;
  subeta = wincon->d2;
  if (wincon->trasc == FCT)
    msubeta = wincon->d4;
  else
    msubeta = osubeta;
  vc = wincon->vc;
  vc1 = wincon->vc1;
  vcs = wincon->vcs;
  vcs1 = wincon->vcs1;
  u1 = windat->u1;
  u = windat->u;
  v = windat->v;
  w = windat->w;
  if (!(wincon->alertf & NONE)) {
    memcpy(windat->salb, windat->sal, window->szc * sizeof(double));
    memcpy(windat->tempb, windat->temp, window->szc * sizeof(double));
  }
  ntbdy = wincon->ntbdy;
  tbdy = wincon->tbdy;

  /*-----------------------------------------------------------------*/
  /* Reset the fluxes for mean velocity transport                    */
  if (wincon->means & TRANSPORT) {
    if (windat->dttr == 0.0)
      return(0);
    set_lateral_BC_vel(master->u1vm, window->nbpte1,
                       window->bpte1, window->bine1);
    set_lateral_BC_vel(master->ume, window->nbpte1,
                       window->bpte1, window->bine1);
    u1 = windat->ume;
    u = windat->u1m;
    v = windat->u2m;
    w = windat->wm;
    ffsl_trans_prep(window, windat, wincon);
    if (wincon->trsplit) {
      ntbdy = wincon->ntbdys;
      tbdy = wincon->tbdys;
    }
  }

  /* Reset the river flow to a parabolic profile if required         */
  reset_flow(window, windat, wincon, 0);

  /*-----------------------------------------------------------------*/
  /* Set up advection scheme specific work arrays.                   */
  courf = 0;
  if (wincon->ultimate || wincon->trasc >= QUICKEST)
    courf = 1;
  if (wincon->trasc & (ORDER1|FCT)) {
    /* ORDER1 : Get the directional flux for upwind schemes (note :  */
    /* these schemes are monotonic and should not be used with the   */
    /* ULTIMATE filter.  */
    /* SIGMA : Don't multiply by the depth here since u1flux3d have  */
    /* already been depth multiplied.                                */
    wincon->ultimate = 0;
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      wincon->w6[e] =
        0.5 * (windat->u1flux3d[e] + fabs(windat->u1flux3d[e]));
      wincon->w7[e] =
        0.5 * (windat->u1flux3d[e] - fabs(windat->u1flux3d[e]));
    }
    courf = 0;
  }
  if (wincon->trasc == LAGRANGE) {
    /* Find the origin of the streamline & grid Courant numbers      */
    semi_lagrange_c(window, windat, wincon);
  } else if (wincon->trasc & VANLEER) {
    /* Get the mean cell spacings                                    */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[window->zm1[c]]);
    }

  } else if (wincon->trasc & FFSL) {
    itermax = 200;
    memset(wincon->nw, 0, window->sgsiz * sizeof(double));
    /* Get the grid spacing.                                         */
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      /* The layer spacing (w9) used in prep_ff_sl() is that of the  */
      /* current layer and that below. dzz uses the layer above.     */
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[window->zm1[c]]);
      zp1 = window->zp1[c];

      /* Calculate cell-centered u1 velocity. The choice of velocity */
      /* for each cell can be either an average of the face values,  */
      /* or one of the face values depending on those values - see   */
      /* Leonard et al. (1996).                                      */
      wincon->nw[c] = 0.0;
      zp1 = window->zp1[c];
      if ((w[c] >= 0.0) && (w[zp1] > 0.0)) 
        wincon->nw[c] = w[c];
      else if ((w[c] < 0.0) && (w[zp1] <= 0.0))
        wincon->nw[c] = w[zp1];
    }
    ffsl_init(window, windat, wincon);
  } else {
    /* Get the mean cell spacings                                    */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zm1 = window->zm1[c];
      wincon->w5[c] = wincon->dz[zm1] / (wincon->dz[zm1] + wincon->dz[c]);
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
    }
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      wincon->d3[e] = 0.5 * window->h1acell[e] / window->h2au1[e];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on  */
  /* the cell centers to calculate the maximum timesteps. The        */
  /* vertical sub-step for the surface layer is not included since   */
  /* the vertical advective terms are not invoked in this layer.     */
  /* Note that the FFSL advection scheme is limited by the velocity  */
  /* gradients (Lipschitz number) not the Courant number.            */
  dtm = dtu = windat->dt;
  if (wincon->means & TRANSPORT)   dtm = dtu = windat->dttr;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF | SUB_STEP_TRACER)) {

    for (cc = 1; cc <= vc; cc++) {
      c = wincon->s1[cc];       /* Wet cell to process               */
      c2 = window->m2d[c];      /* 2D cell corresponding to 3D cell  */
      zp1 = window->zp1[c];
      /*
      for (j = 1; j <= window->npe[c2] / 2; j++) {
	e1 = window->c2e[j][c];
	j1 = jo(j, window->npe[c2]);
	e2 = window->c2e[j1][c];
	d4 = 0.5 * (u1[e1] + u1[e2]);
	if (wincon->trasc == FFSL)
	  d4 = max(fabs(window->eSc[j][c2] * u1[e1] +
			window->eSc[j1][c2] * u1[e2]), 
		   wincon->u1kh[c]/ (2.0 * window->hacell[j][c2]));

	if (fabs(d4) < minval)
	  d4 = minval;
	d1 = sf * fabs(2.0 * window->hacell[j][c2] / d4);
	if (d1 < dtm) {
	  dtm = d1;
	  sprintf(vt, "u");
	  vm = d4;
	  vs = 2.0 * window->hacell[j][c2];
	  ii = c2;
	  jj = window->m2de[e1];
	  kk = window->s2k[c];
	}
      }
      */
      for (j = 1; j <= window->npe[c2]; j++) {
	e1 = window->c2e[j][c];
	e2 = window->m2de[e1];
	if (wincon->trasc & FFSL) {
	  int c1 = window->c2c[j][c];
	  d4 = (u[c] * window->costhu1[e2] + v[c] * window->sinthu1[e2]) -
  	       (u[c1] * window->costhu1[e2] + v[c1] * window->sinthu1[e2]);
	} else
	  d4 = u1[e1];

	if (fabs(d4) < minval)
	  d4 = minval;
	d1 = sf * fabs(window->h2au1[e2] / d4);
	if (d1 < dtm) {
	  dtm = d1;
	  sprintf(vt, "u");
	  vm = d4;
	  vs = window->h2au1[e2];
	  ii = c2;
	  jj = window->m2de[e1];
	  kk = window->s2k[c];
	}
      }

      /* Minimum time-step due to z velocity.                        */
      /* Note : surface tracer values are calculated on the basis of */
      /* mass conservation, hence Courant violations need not be     */
      /* considered.                                                 */
      if (wincon->trasc & FFSL)
	d4 = (cc <= vcs) ? 0.0 : fabs(w[zp1] - w[c]);
      else
	d4 = (cc <= vcs) ? 0.0 : 0.5 * (w[c] + w[zp1]);
      if (fabs(d4) < minval)
        d4 = minval;
      /* SIGMA : Multiply by depth                                   */
      d3 = (cc <= vcs) ? dtu : sf * fabs(wincon->dz[c] * wincon->Ds[c2] / d4);
      if (d3 > 0.0 && d3 < dtm) {
        dtm = d3;
        sprintf(vt, "w");
        vm = d4;
        vs = wincon->dz[c] * wincon->Ds[c2];
        ii = c;
        jj = window->m2de[e1];
        kk = window->s2k[c];
      }
    }

    if (dtm != dtu) {
      dtm = floor(dtm);
      hd_warn
        ("Sub-time stepping for tracers (%s=%e) at %8.3f days (c2=%3d e2=%3d k=%3d): h=%e dt=%5.2f\n",
         vt, vm, windat->t / 86400.0, ii, jj, kk, vs, dtm);
    }
  } else if (wincon->stab & (MULTI_DT)) {
    dtm = windat->dts;
  }
  if (dtm == 0.0 || (int)(dtu / dtm) > itermax) {
    c2 = window->m2d[ii];
    write_site(window, window->cellx[c2], window->celly[c2], "Tracer sub-step");
    hd_quit_and_dump
      ("tracer_step_3d: maximum number of sub-steps (%d) exceeded at %8.3f days (c2=%3d [%f %f] e2=%3d %3d)\n",
       itermax, windat->t / 86400.0, ii, window->cellx[c2], window->celly[c2], jj, kk);
      return(1);
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Include surface layer in the advection.                 */
  ksf = vcs + 1;

  /* Set the system up for semi-Lagrangian advection                 */
  if (wincon->trasc == LAGRANGE) {
    /* Courant numbers may be > 1                                    */
    courf = 0;
    /* No special treatment of surface layer                         */
    if (!wincon->sigma)
      ksf = 1;
    /* No ultimate limiting                                          */
    wincon->ultimate = 0;
    /* No sub-stepping (LAGRANGE is unconditionally stable)          */
    dtm = dtu = windat->dttr;
  }
  
  /*-----------------------------------------------------------------*/
  /* Initialise the surface elevation. The surface concentration     */
  /* calculated in surf_concs() must also be sub-stepped. This means */
  /* the volume of the previous sub-step must be known to get the    */
  /* concentration at the current sub-step. The elevations from the  */
  /* previous sub-steps are stored in osubeta, allowing previous     */
  /* sub-step volumes to be calculated. It is assumed that the       */
  /* velocities contributing to the fluxes into the surface layer    */
  /* are invariant over the full time-step, windat->dt.              */
  for (cc = 1; cc <= vcs; cc++) {
    c = wincon->s1[cc];         /* 3D surface wet cell to process    */
    c2 = window->m2d[c];        /* 2D cell corresponding to 3D cell  */
    osubeta[c2] = wincon->oldeta[c2];
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Get the water depth at the forward time                 */
  if (wincon->sigma) {
    for (cc = 1; cc <= vcs; cc++) {
      c = wincon->s1[cc];       /* 3D surface wet cell to process    */
      c2 = window->m2d[c];      /* 2D cell corresponding to 3D cell  */
      wincon->Hn1[c2] = windat->eta[c2] - wincon->Hs[c2];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  trem = dtu;
  while (trem > 0) {
    /* Get the sub-timestep                                          */
    dtu = dtm;
    if (trem < dtu)
      dtu = trem;
    wincon->b1 = dtu;

    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Get the courant numbers for this sub-step                   */
      if (courf) {
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,c2)
#endif
	for (ee = 1; ee <= window->n3_e1; ee++) {
	  e = window->w3_e1[ee];
	  e2 = window->m2de[e];
	  wincon->w6[e] = fmod(u1[e] * dtu / window->h2au1[e2], 1.0);
	}
	for (cc = 1; cc <= window->n3_t; cc++) {
	  c = window->w3_t[cc];
	  c2 = window->m2d[c];
	  wincon->w8[c] = fmod(w[c] * dtu / (wincon->w9[c] * wincon->Ds[c2]), 1.0);
	}
      }

      /*-------------------------------------------------------------*/
      /* Get the surface elevation for the sub-step.                 */
      /* Note : the tracer is updated in the surface layer in the    */
      /* main loop for the sigma system. Special treatment of the    */
      /* surface layer, and hence the sub-step elevation, is not     */
      /* required.                                                   */
      if (!wincon->sigma) {
	double dt = (wincon->means & TRANSPORT) ? windat->dttr : windat->dt;
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,c2)
#endif
	for (cc = 1; cc <= vcs; cc++) {
	  c = wincon->s1[cc];
	  c2 = window->m2d[c];
	  subeta[c2] =
	    wincon->oldeta[c2] + (1.0 + (dtu - trem) / dt) *
	    (windat->eta[c2] - wincon->oldeta[c2]);
	}
	if (wincon->trasc == FCT)
	  memcpy(msubeta, subeta, window->szcS * sizeof(double));
      }
      
      /*-------------------------------------------------------------*/
      /* Save the region volume fluxes                               */
      region_volume_flux_coup(window, windat, wincon, dtu, trem);
      
      /*-------------------------------------------------------------*/
      /* If the FFSL scheme is used, calculate trajectories          */
      if (wincon->trasc & FFSL)
	prep_ff_sl(window, windat, wincon, dtu);
      
      /* Check the volume continuity balance if required             */
      if (wincon->numbers1 & VOLCONT) 
	check_volcons(window, windat, wincon, osubeta, subeta, dtu);

      /*-------------------------------------------------------------*/
      /* Tracer loop                                                 */
      TIMING_SET;
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,cc,cb,c2,zp1,surftr,top)
#endif
      for (nn = 0; nn < ntbdy; nn++) {
	int n = tbdy[nn];
	double *tr = windat->tr_wc[n]; /* Tracer values              */
	tr[0] = (double)n;
	int kef, kefs, bgzf;           /* Include / exclude OBCs     */
	double *Fx = wincon->w2;       /* x tracer flux              */
	double *Fz = wincon->w4;       /* z tracer flux              */
	double *Fxh = wincon->Fxh;
	double *Fzh = wincon->Fzh;
	double *dtracer = wincon->w1;
	/* Modified tracer concentrations                            */
	double *tr_mod = wincon->tr_mod;
	double *tr_mod_x = wincon->tr_mod_x;
	double *tr_mod_z = wincon->tr_mod_z;

	/* Re-wire work arrays if parallel mode                      */
#if defined(HAVE_OMP)
	if (master->trans_num_omp > 1) {
	  int nomp = omp_get_thread_num();
	  Fx = wincon->w2n[nomp];
	  Fz = wincon->w4n[nomp];
	  dtracer = wincon->w1n[nomp];
	  tr_mod = wincon->tr_modn[nomp];
	  tr_mod_x = wincon->tr_modn_x[nomp];
	  tr_mod_z = wincon->tr_modn_z[nomp];
	}
#endif	
	
	/* Initialize                                                */
	if (windat->fluxe1) memset(windat->fluxe1, 0, window->szc * sizeof(double));
	if (windat->fluxe2) memset(windat->fluxe2, 0, window->szc * sizeof(double));
	if (windat->fluxw) memset(windat->fluxw, 0, window->szc * sizeof(double));
	if (windat->fluxkz) memset(windat->fluxkz, 0, window->szc * sizeof(double));
	for (tc = 0; tc < wincon->rkstage; tc++)
	  memset(wincon->tr_gr[tc], 0, window->szc * sizeof(double));

	/* Rugne-Kutta stage loop                                    */
	for (tc = 0; tc < wincon->rkstage; tc++) {

	  /*---------------------------------------------------------*/
	  /* Initialise the advective terms                          */
	  memset(Fx, 0, window->sze * sizeof(double));
	  memset(Fz, 0, window->szc * sizeof(double));
	  memset(dtracer, 0, window->szc * sizeof(double));
	  if (wincon->trasc & FCT) {
	    memset(Fxh, 0, window->sze * sizeof(double));
	    memset(Fzh, 0, window->szc * sizeof(double));
	  }
	  
	  /*---------------------------------------------------------*/
	  /* Set whether to include boundary cells or not            */
	  bgzf = 0;
	  kef = vc1;
	  kefs = vcs1;
	  for (cc = 0; cc < window->nobc; cc++) {
	    if (window->open[cc]->bcond_tra[n] & (TRCONC|TRFLUX|TRCONF)) {
	      bgzf = 1;
	      kef = vc;
	      kefs = vcs;
	      save_OBC_tr(window, windat, wincon, Fx, tr, n, 1);
	      break;
	    }
	  }
	  
	  /*---------------------------------------------------------*/
	  /* Set the bottom boundary condition (no-gradient)         */
	  for (cc = 1; cc <= window->b2_t; cc++) {
	    cb = window->bot_t[cc];
	    tr[window->zm1[cb]] = tr[cb];
	  }
	  /*
	  for (cc = 1; cc <= window->nbpt; cc++) {
	    int c1 = window->bpt[cc];
	    int c2 = window->bin[cc];
	    tr[c1] = tr[c2];
	  }
	  */
	  /*---------------------------------------------------------*/
	  /* Get the tracer values in multi-dt auxiliary cells. This */
	  /* is a linear interpolation between the value at the      */
	  /* start of the timestep and the end of longer timesteps.  */
	  set_multidt_t(window, windat, tr, trem, n);

	  /*---------------------------------------------------------*/
	  /* For unstructured high order schemes get the values at   */
	  /* cp and cm are reconstructed using quadratic least       */
	  /* squares.                                                */
	  if (wincon->trasc & HIORDER) {
	    for (ee = 1; ee <= window->n3_e1; ee++) {
	      e = window->w3_e1[ee];
	      wincon->trp[e] = get_quadratic_value(window, windat, wincon, tr, e, 0);
	      wincon->trm[e] = get_quadratic_value(window, windat, wincon, tr, e, 1);
	    }
	  }
	  /*
	  memcpy(wincon->tr_rk[tc], tr, window->szc * sizeof(double));
	  */
	  /*---------------------------------------------------------*/
	  /* Get the advective fluxes                                */
	  //TIMING_SET;
	  if (wincon->advect[n]) {
	    if (wincon->trasc & FCT) {
	      advect(window, windat, wincon, tr, Fx, Fz, 
		     tr_mod, tr_mod_x, tr_mod_z, dtu, ORDER1);
	      /* High order solution                                 */
	      if (wincon->trasc & ORDER3US)
		advect(window, windat, wincon, tr, Fxh, Fzh, 
		       tr_mod, tr_mod_x, tr_mod_z, dtu, ORDER3US);
	      else if (wincon->trasc & ORDER4US)
		advect(window, windat, wincon, tr, Fxh, Fzh, 
		       tr_mod, tr_mod_x, tr_mod_z, dtu, ORDER4US);
	      else
		advect(window, windat, wincon, tr, Fxh, Fzh, 
		       tr_mod, tr_mod_x, tr_mod_z, dtu, ORDER2);
	    } else
	      advect(window, windat, wincon, tr, Fx, Fz, 
		     tr_mod, tr_mod_x, tr_mod_z, dtu, wincon->trasc);
	  }
	  // TIMING_DUMP(2, "    advect");

	  /*---------------------------------------------------------*/
	  /* Save the advective fluxes on OBCs                       */
	  if (bgzf) save_OBC_tr(window, windat, wincon, Fx, tr, n, 3);

	  /*---------------------------------------------------------*/
	  /* Save the advective fluxes if required                   */
	  if (wincon->trflux == n) {
	    calc_flux(window, windat, wincon, Fx, Fz, dtu, 
		      wincon->trfd1, wincon->trfd2);
	  }

	  /*---------------------------------------------------------*/
	  /* Save the region fluxes                                  */
	  region_flux_coup(window, windat, wincon, Fx, Fz, dtu, n);
    
	  /*---------------------------------------------------------*/
	  /* Get the horizontal diffusive fluxes                     */
	  if (wincon->diffuse[n])
	    hor_diffuse(window, windat, wincon, tr, Fx);
	  
	  /*---------------------------------------------------------*/
	  /* Reset the advective fluxes on OBCs if required          */
	  if (bgzf) save_OBC_tr(window, windat, wincon, Fx, tr, n, 4);

	  /*---------------------------------------------------------*/
	  /* Store the updated horizontal advective fluxes           */
	  /* Note : land cells and cells immediately adjacent to     */
	  /* open boundaries are not updated.                        */
	  for (cc = 1; cc <= window->b3_t; cc++) {
	    c = window->w3_t[cc];
	    c2 = window->m2d[c];
	    dtracer[c] = 0.0;
	    for (j = 1; j <= window->npe[c2]; j++) {
	      e = window->c2e[j][c];
	      dtracer[c] += (window->eSc[j][c2] * Fx[e] * dtu);
	    }
	  }

	  /*---------------------------------------------------------*/
	  /* Invoke the FCT limiter if required                      */
	  if (wincon->trasc & FCT) 
	    fct_limit(window, windat, wincon, tr, ksf, kef, kefs, dtu, osubeta, msubeta);

	  /*---------------------------------------------------------*/
	  /* Reset the flux divergence at open boundary cells (if    */
	  /* required).                                              */
	  reset_tr_OBCflux(window, windat, wincon, dtracer, Fx, dtu, n);

	  /*---------------------------------------------------------*/
	  /* Evaluate the sources and sinks of tracer                */
	  ss_tracer(window, windat, wincon, n, dtracer, dtu);

	  /*---------------------------------------------------------*/
	  /* Get the updated tracer concentration                    */
	  /* Note : land cells and cells immediately adjacent to     */
	  /* open boundaries are not updated. Note, the first vcs    */
	  /* cells in the cells to process vector, wincon->s1,       */
	  /* contain surface cells and are not processed here since  */
	  /* the surface has separate treatment.                     */
	  for (cc = ksf; cc <= kef; cc++) {
	    c = wincon->s1[cc];   /* Wet cell to process             */
	    c2 = window->m2d[c];  /* 2D cell corresponding to 3D cell*/
	    zp1 = window->zp1[c]; /* Cell above cell c               */
	    wincon->tr_rk[tc][c] = tr[c];
	    /* SIGMA : Adjust tracer values for the depth            */
	    if (slf)
	      tr[c] *= wincon->Ds[c2];
	    else
	      tr[c] *= wincon->Hn1[c2];

	    wincon->tr_gr[tc][c] = -((dtracer[c] / window->cellarea[c2] +
			     (Fz[zp1] - Fz[c])) / wincon->dz[c]);
	    tr[c] = runge_kutta(window, windat, wincon, 1.0, 1.0, tc, c, slf);
	    /*tr[c] /= wincon->Hn1[c2];*/
	  }

	  /*---------------------------------------------------------*/
	  /* Set all layers above the surface to surface layer       */
	  if (wincon->trasc == LAGRANGE) {
	    if (!wincon->sigma && wincon->advect[n]) {
	      for (cc = 1; cc <= vcs1; cc++) {
		c = wincon->s1[cc];
		surftr = tr[c];
		while (c != window->zp1[c]) {
		  c = window->zp1[c];
		  tr[c] = surftr;
		}
	      }
	    }
	  }
	  
	  /*---------------------------------------------------------*/
	  /* Get the tracer concentration in the surface cells       */
	  else {
	    if (!wincon->sigma && wincon->advect[n]) {
	      /* Use 3rd order at the surface if 4th order is set.   
	      int orks = wincon->rkstage;
	      if (orks == 5) {
		wincon->rkstage = 3;
		if (tc > 2) continue;
	      }
	      */
	      for (cc = 1; cc <= kefs; cc++) {
		c = wincon->s1[cc];
		c2 = window->m2d[c];
		top = -Fz[c] * window->cellarea[c2];
		/* Get the mass at time t, sum of fluxes and volumes */
		if (wincon->rkstage == 0)
		  surftr =
		    surf_conc(window, windat, wincon, tr, msubeta[c2], 
			      subeta[c2], c, c2, cc, top, dtracer);
		else {
		  wincon->tr_rk[tc][c] =
		    surf_concrk(window, windat, wincon, tr, msubeta[c2], 
				subeta[c2], c, c2, cc, wincon->tr_gr[tc], 
				&ovol, &nvol, top, dtracer);

		  /* Get the new concentration for this stage        */
		  surftr = runge_kutta(window, windat, wincon, ovol, nvol, tc, c, slf);
		}
		tr[c] = (window->botz[c2] > 0.0 && wincon->dz[c] < wincon->hmin) ?
		  tr[c] : surftr;
	      }
	      /*if (orks == 5) wincon->rkstage = orks;*/
	    }
	    /*if(windat->days>10.022)exit(0);*/
	    if (wincon->sigma && wincon->advect[n]) {
	      /* SIGMA surface calculation : w=0 at the free surface */
	      /* Fz[zp1]=0 at this layer.                            */
	      for (cc = 1; cc < ksf; cc++) {
		c = wincon->s1[cc];
		c2 = window->m2d[c];
	    
		/* Diagnostic to check for continuity using tracer   */
		/* fluxes (tracer value must = constant).            */
		/*
		  d1=dtracer[c]/(wincon->dz[c]*tr[c]);
		  d2=-Fz[c]*window->cellarea[c2]/(wincon->dz[c]*tr[c]);
		  d3=(wincon->Hn1[c2]-wincon->Ds[c2])*window->cellarea[c2];
		  top=d1+d2+d3; */
		/* Diagnostic to check for continuity using momentum */
		/* fluxes.                                           */
		/*
		  d1=wincon->dz[c]*window->cellarea[c2]*
		  (wincon->Hn1[c2]-wincon->Ds[c2])/windat->dt;
		  d2 = 0.0;
		  for (j = 1; j <= window->npe[c2]; j++) {
		  e = window->c2e[j][c];
		  d2 += (window->eSc[j][c2] * windat->u1flux3d[e]);
		  }
		  d3=-windat->w[c]*window->cellarea[c2]; top=d1+d2+d3;
		  if(c2==2&&n==0)printf("%f %d %d %e\n",windat->t/86400,
		  window->s2k[c],c,top); */
		if (slf)
		  tr[c] *= wincon->Ds[c2];
		else
		  tr[c] *= wincon->Hn1[c2];
		tr[c] -=
		  ((dtracer[c] / window->cellarea[c2] -
		    Fz[c]) / wincon->dz[c]);
		tr[c] /= wincon->Hn1[c2];
	      }
	    }
	  }
	  
	  /* Reset boundary tracer values if required                */
	  if (bgzf) save_OBC_tr(window, windat, wincon, Fx, tr, n, 2);
	  if (cdb) printf("end %f\n",tr[cdb]);
	  /* Clip the tracer for FFSL schemes                        */
	  if (wincon->trasc & FFSL && wincon->fillf & (MONOTONIC|CLIP))
	    clip_ffsl(window, windat, wincon, tr);

	}                       /* rkstage loop end                  */
      }                         /* Tracer loop end                   */
      TIMING_DUMP(2, "    t-loop");

      /* Clip the tracer concentration if required                   */
      TIMING_SET;
      tr_bounds(window, windat, wincon);
      TIMING_DUMP(2, "    tr_bounds");
      
      /* Set the old sub-stepped elevation                           */
      if (!wincon->sigma) {
	for (cc = 1; cc <= vcs; cc++) {
	  c = wincon->s1[cc];
	  c2 = window->m2d[c];
	  osubeta[c2] = subeta[c2];
	}
      }
    }
    /* dtu!=0 loop                                                   */
    /* Decrement the time remaining                                  */
    trem -= dtu;
    slf = 0;
  }                             /* while(trem>0) loop                */

  /*-----------------------------------------------------------------*/
  /* Reset the new surface coordinate at thin layers if required     */
  if (wincon->thin_merge) {
    for (vc = 1; vc <= wincon->nkth_e1; vc++) {
      cc = (int)wincon->kth_e1[vc];
      c = window->nsur_t[cc];
      c2 = window->m2d[c];
      window->nsur_t[cc] = window->zp1[c];
      zm1 = window->zm1[window->nsur_t[cc]];
      wincon->dz[zm1] = window->gridz[window->nsur_t[cc]] - 
	max(window->gridz[zm1], window->botz[c2]);	
    }
  }

  /*-----------------------------------------------------------------*/
  /* Get the density at boundary ghost cells                         */
  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn]; 
    if (open->bgz) {
      for (ee = 1; ee <= open->no2_e1; ee++) {
	e = open->obc_e1[ee];
	c = open->obc_e2[ee];
	density_gc(window, windat, wincon, c, open->omape[ee]);
      }
    }
  }

  /* Revert river flow to original profile if required               */
  reset_flow(window, windat, wincon, 1);
  /*
  if (wincon->means & TRANSPORT)
    ffsl_trans_post(window, windat, wincon);
  */

  return(0);
}

/* END advect_diffuse()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implement the Runge-Kutta method of nth order and rkstage stages. */
/* See Shu and Osher (1998) J. Compt. Phys. 77, 439-471. Various     */
/* schemes are also described in Gottlieb (2005) J. Scientific       */
/* Comput. 25, 105-128. The 5 stage method is that of Spiteri and    */
/* Ruuth (2002) SIAM J. Numer. Anal. 40, 469-491.                    */
/*-------------------------------------------------------------------*/
double runge_kutta(geometry_t *window,   /* Window geometry          */
		   window_t *windat,     /* Window data              */
		   win_priv_t *wincon,   /* Window constants         */
		   double ovol,          /* Volume at t              */
		   double nvol,          /* Volume at t+1            */
		   int tc,               /* Stage counter            */
		   int c,                /* Mesh coordinate          */
		   int slf               /* 1st substep flag         */
		   )
{
  int c2 = window->m2d[c];
  double **ntr = wincon->tr_rk;          /* Stage tracer values      */
  double **grad = wincon->tr_gr;         /* Stage time tendency      */
  double sf = (slf) ? wincon->Ds[c2] : wincon->Hn1[c2];
  int rki = wincon->rkstage - 1;
  int rks = tc * 3;
  double tr, vol;

  /* Note: if nvol != ovol then we are solving:                      */
  /* c(n+1)vol = c(n)ovol + grad                                     */
  /* In this case vol must be computed as some intermediate value    */
  /* between olvol and nvol.                                         */
  /* Get the volume corresponding to the Runge-Kutta stage. This is  */
  /* a linear interpolation from old (t) and new (t+1) volumes with  */
  /* the linear increment added to the old volume given by the       */
  /* weight corresponding to grad[].                                 */
  vol = ovol + rkweights[rks+2][rki] * (nvol - ovol);

  /* Compute a new tracer from previous weighted values.             */
  if (wincon->rkstage == 5 && tc == 4) {
    vol = ovol + rkweights[rks+5][rki] * (nvol - ovol);
    tr = sf * rkweights[rks][rki] * ntr[2][c] + 
         sf * rkweights[rks+1][rki] * ntr[3][c] +
	      rkweights[rks+2][rki] * grad[3][c] +
	 sf * rkweights[rks+3][rki] * ntr[4][c] +
	      rkweights[rks+2][rki] * grad[4][c];
    tr /= (wincon->Hn1[c2] * vol);
  } else {
    tr = sf * rkweights[rks][rki] * ntr[0][c] + 
	 sf * rkweights[rks+1][rki] * ntr[tc][c] +
	      rkweights[rks+2][rki] * grad[tc][c];
    tr /= (wincon->Hn1[c2] * vol);
  }
  return(tr);
}

/* END runge_kutta()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Finds the volume weight for an idealised setup where initial      */
/* concentration and mass = 1, initial volume = 1, final volume = 2, */
/* final concentration = 1, final mass = 2 and flux = 1.             */
/* This should be the weight applied to volume in runge_kutta().     */
/*-------------------------------------------------------------------*/
void test_rk(int rkstage)
{
  int n;
  double im = 1.0;
  double flux = 1.0;
  double ovol = 1.0;
  double nvol = 2.0;
  int rki = rkstage - 1;
  int rks;
  double vol;

  for (n = 0; n < rkstage; n++) {
    rks = n * 3;
    vol = rkweights[rks][rki] * im + 
	 rkweights[rks+1][rki] * im +
	      rkweights[rks+2][rki] * flux;
    if (rkstage == 5 && n == 4) {
      vol = rkweights[rks][rki] * im + 
	rkweights[rks+1][rki] * im +
	rkweights[rks+2][rki] * flux +
	rkweights[rks+3][rki] * im +
	rkweights[rks+2][rki] * flux;
    }
    printf("stage%d volume weight = %17.15f\n",n, (vol-ovol)/(nvol-ovol));
  }
  hd_quit("RK test done\n");
}

/* END test_rk()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion                                */
/*-------------------------------------------------------------------*/
int advect_diffuse_2d(geometry_t *window,  /* Window geometry        */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon   /* Window constants       */
  )
{
  int nn;                       /* Tracer index                      */
  int c, cc;                    /* Sparse indices                    */
  int vcs;                      /* Surface tracer cells to process   */
  int ee, e, e1, e2;            /* Edge coordinate, counter          */
  int j, j1;
  int vcs1;                     /* As for vcs excluding OBC cells    */
  double *Fx;                   /* x tracer flux                     */
  double *Fz;                   /* z tracer flux                     */
  double *dtracer;              /* Advective terms for x,y directions*/
  double dtu;                   /* Sub-time step to use              */
  double dtm;                   /* Minimum allowable time step       */
  double trem;                  /* Time remaining                    */
  double *osubeta;              /* Old surface elevation at sub-timestep */
  double *subeta;               /* Surface elevation at sub-timestep */
  double d1, d2, d4;            /* Dummy variables                   */
  double sf = 0.8;              /* Fraction of sub-timestep to use   */
  double minval = 1e-10;        /* Minimum value for velocity        */
  int slf = 1;                  /* Set to zero on the first sub-step */
  int itermax = 20;             /* Maximum number of substeps        */
  int courf;                    /* Flag to calculate Courant numbers */
  int ii = 0, jj = 0, kk = 0;   /* (i,j,k) location of sub-step violation */
  double vm = 0.0, vs = 0.0;    /* Sub-step Courant number & grid spacing */
  char vt[4];                   /* Name of sub-step velocity component */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* w1 = dtracer                                                    */
  /* w2 = Fx                                                         */
  /* w4 = Fz                                                         */
  /* w5 = Mean vertical cell spacing                                 */
  /* w6 = x direction Courant number                                 */
  /* w7 = y direction Courant number                                 */
  /* w8 = z direction Courant number                                 */
  /* w9 = vertical cell spacing between grid centers                 */
  /* d1 = osubeta                                                    */
  /* d2 = subeta                                                     */
  /* d3 = mean cell spacing in the x direction                       */
  /* d4 = mean cell spacing in the y direction                       */
  /* s1 = wet cells to process for tracers                           */

  /*-----------------------------------------------------------------*/
  /* Assign memory for the horizontal advective terms                */
  dtracer = wincon->w1;
  Fx = wincon->w2;
  Fz = wincon->w4;
  osubeta = wincon->d1;
  subeta = wincon->d2;
  vcs = wincon->vcs;
  vcs1 = wincon->vcs1;

  /*-----------------------------------------------------------------*/
  /* Set up advection scheme specific work arrays.                   */
  if (wincon->trasc == ORDER1) {
    /* ORDER1 : Get the directional flux for upwind schemes (note :  */
    /* these schemes are monotonic and should not be used with the   */
    /* ULTIMATE filter.                                              */
    /* SIGMA : Don't multiply by the depth here since u1flux3d have  */
    /* already been depth multiplied.                                */
    wincon->ultimate = 0;
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      wincon->w6[e] =
        0.5 * (windat->u1flux3d[e] + fabs(windat->u1flux3d[e]));
      wincon->w7[e] =
        0.5 * (windat->u1flux3d[e] - fabs(windat->u1flux3d[e]));
    }
  } else if (wincon->trasc & (VANLEER|FFSL)) {
    /* Get the mean cell spacings                                    */
    memcpy(wincon->w9, wincon->dz, window->szcS * sizeof(double));
  } else {
    memcpy(wincon->w5, wincon->dz, window->szcS * sizeof(double));
    memcpy(wincon->w9, wincon->dz, window->szcS * sizeof(double));
    for (ee = 1; ee <= window->n2_e1; ee++) {
      e = window->w2_e1[ee];
      wincon->d3[e] = 0.5 * window->h1acell[e] / window->h2au1[e];
    }
  }
  courf = 0;
  if (wincon->ultimate || wincon->trasc >= QUICKEST)
    courf = 1;

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on  */
  /* the cell centers to calculate the maximum timesteps. The        */
  /* vertical sub-step for the surface layer is not included since   */
  /* the vertical advective terms are not invoked in this layer.     */
  dtm = dtu = windat->dt;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF | SUB_STEP_TRACER)) {

    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D cell corresponding to 3D cell  */

      /* Minimum time-step due to horizontal velocity                */
      for (j = 1; j <= window->npe[c] / 2; j++) {
	e1 = window->c2e[j][c];
	j1 = jo(j, window->npe[c]);
	e2 = window->c2e[j1][c];

	d4 = 0.5 * (windat->u1[e1] + windat->u1[e2]);
	if (wincon->trasc == FFSL)
	  d4 = max(fabs(window->eSc[j][c] * windat->u1[e1] +
			window->eSc[j1][c] * windat->u1[e2]), 
		   wincon->u1kh[c]/ (2.0 * window->hacell[j][c]));
	if (fabs(d4) < minval)
	  d4 = minval;
	d1 = sf * fabs(2.0 * window->hacell[j][c] / d4);
	if (d1 < dtm) {
	  dtm = d1;
	  sprintf(vt, "u");
	  vm = d4;
	  vs = 2.0 * window->hacell[j][c];
	  ii = window->s2i[c];
	  jj = window->s2j[c];
	  kk = window->s2k[c];
	}
      }
    }
    if (dtm != dtu) {
      hd_warn
        ("Sub-time stepping for tracers (%s=%6.3f) at %8.3f days (%3d %3d %3d): dt=%5.2f\n",
         vt, vm, windat->t / 86400.0, ii, jj, kk, dtm);
    }
  } else if (wincon->stab & (MULTI_DT)) {
    dtm = windat->dts;
  }
  if ((int)(dtu / dtm) > itermax) {
    hd_quit_and_dump
      ("trhadvec: maximum number of sub-steps (%d) exceeded.\n", itermax);
    return(1);
  }

  /*-----------------------------------------------------------------*/
  /* Initialise the surface elevation. The surface concentration     */
  /* calculated in surf_concs() must also be sub-stepped. This means */
  /* the volume of the previous sub-step must be known to get the    */
  /* concentration at the current sub-step. The elevations from the  */
  /* previous sub-steps are stored in osubeta, allowing previous     */
  /* sub-step volumes to be calculated. It is assumed that the       */
  /* velocities contributing to the fluxes into the surface layer    */
  /* are invariant over the full time-step, windat->dt.              */
  if (wincon->sigma)
    osubeta = wincon->one;
  else {
    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D surface wet cell to process    */
      osubeta[c] = wincon->oldeta[c] - window->botz[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Get the water depth at the forward time                 */
  if (wincon->sigma)
    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D surface wet cell to process    */
      wincon->Hn1[c] = windat->eta[c] - wincon->Hs[c];
    }

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains               */
  trem = dtu;
  while (trem > 0) {
    /* Get the sub-timestep                                          */
    dtu = dtm;
    if (trem < dtu)
      dtu = trem;

    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Get the courant numbers for this sub-step                   */
      if (courf) {
	for (ee = 1; ee <= window->n2_e1; ee++) {
	  e = window->w2_e1[ee];
	  e2 = window->m2de[e];
          wincon->w6[e] = windat->u1av[e] * dtu / window->h2au1[e2];
	}
        for (cc = 1; cc <= window->n2_t; cc++) {
          c = window->w2_t[cc];
          wincon->w8[c] = 0.0;
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the surface elevation for the sub-step.                 */
      /* Note : the tracer is updated in the surface layer in the    */
      /* main loop for the sigma system. Special treatment of the    */
      /* surface layer, and hence the sub-step elevation, is not     */
      /* required.                                                   */
      if (!wincon->sigma) {
        for (cc = 1; cc <= window->v2_t; cc++) {
          c = window->w2_t[cc];
          subeta[c] = wincon->oldeta[c] + (1.0 + (dtu - trem) / windat->dt) * 
	    (windat->eta[c] - wincon->oldeta[c]) - window->botz[c];
        }
      } else
        subeta = wincon->one;

      /*-------------------------------------------------------------*/
      /* Tracer loop                                                 */
      for (nn = 0; nn < wincon->ntbdy; nn++) {
        int n = wincon->tbdy[nn];
        double *tr = windat->tr_wc[n];  /* Tracer values             */

        /*-----------------------------------------------------------*/
        /* Initialise the advective terms                            */
        memset(Fx, 0, window->sze * sizeof(double));
        memset(dtracer, 0, window->szc * sizeof(double));

        /*-----------------------------------------------------------*/
        /* Get the tracer values in multi-dt auxiliary cells. This   */
        /* is a linear interpolation between the value at the start  */
        /* of the timestep and the end of longer timesteps.          */
        set_multidt_t(window, windat, tr, trem, n);

        /*-----------------------------------------------------------*/
        /* Get the advective fluxes                                  */
        if (wincon->advect[n])
          advect(window, windat, wincon, tr, Fx, Fz, 
		 NULL, NULL, NULL, dtu, wincon->trasc);

        /*-----------------------------------------------------------*/
        /* Get the horizontal diffusive fluxes                       */
        if (wincon->diffuse[n])
          hor_diffuse_2d(window, windat, wincon, tr, Fx);

        /*-----------------------------------------------------------*/
        /* Store the updated horizontal advective fluxes             */
        /* Note : land cells and cells immediately adjacent to open  */
        /* boundaries are not updated.                               */
        for (cc = 1; cc <= window->v2_t; cc++) {
          c = window->w2_t[cc];
	  dtracer[c] = 0.0;
	  for (j = 1; j <= window->npe[c]; j++) {
	    e = window->c2e[j][c];
	    dtracer[c] += (window->eSc[j][c] * Fx[e] * dtu);
	  }
        }

        /*-----------------------------------------------------------*/
        /* Evaluate the sources and sinks of tracer                  */
        ss_tracer(window, windat, wincon, n, dtracer, dtu);

        /*-----------------------------------------------------------*/
        /* Get the updated tracer concentration                      */
        /* Note : land cells and cells immediately adjacent to open  */
        /* boundaries are not updated. Note, the first vcs cells in  */
        /* the cells to process vector, wincon->s1, contain surface  */
        /* cells and are not processed here since the surface has    */
        /* separate treatment.                                       */
        for (cc = 1; cc <= window->v2_t; cc++) {

          c = window->w2_t[cc]; /* Wet cell to process               */

          /* SIGMA : Adjust tracer values for the depth              */
          if (slf)
            tr[c] *= wincon->Ds[c];
          else
            tr[c] *= wincon->Hn1[c];

          tr[c] = (tr[c] * osubeta[c] - dtracer[c] /
                   window->cellarea[c]) / subeta[c];

          tr[c] /= wincon->Hn1[c];
        }
      }                         /* Tracer loop end                   */

      /* Clip the tracer concentration if required                   */
      tr_bounds(window, windat, wincon);

      /* Set the old sub-stepped elevation                           */
      if (!wincon->sigma) {
        for (cc = 1; cc <= window->v2_t; cc++) {
          c = window->w2_t[cc];
          osubeta[c] = subeta[c];
        }
      }
    }
    /* dtu!=0 loop                                                   */
    /* Decrement the time remaining                                  */
    trem -= dtu;
    slf = 0;
  }                             /* while(trem>0) loop                */

  /*-----------------------------------------------------------------*/
  /* Get the mixed layer depth diagnostic if required                */
  if (!(wincon->mixlayer & NONE))
    get_mixed_layer(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the flushing diagnostic if required                   */
  if (wincon->trflsh)
    total_mass(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Evaluate the tracer decay                                       */
  for (nn = 0; nn < wincon->ntdec; nn++) {
    int n = wincon->tdec[nn];
    tracer_decay(window, windat, wincon, wincon->dectr[n],
                 windat->tr_wc[n], &wincon->trinfo_2d[n]);
  }
  return(0);
}

/* END advect_diffuse_2d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion for trsplit == 1. Here         */
/* temperature and salinity are first advected with the alternate    */
/* scheme and the remaining tracers advected using FFSL. The         */
/* arrays of tracers to advect, tbdy, are modified for two calls to  */
/* the routine advect_diffuse().                                     */
/*-------------------------------------------------------------------*/
int advect_diffuse_split(geometry_t *window,  /* Window geometry     */
			 window_t *windat,    /* Window data         */
			 win_priv_t *wincon   /* Window constants    */
  )
{
  int trasco, uo;
  double dto;

  /* Copy the original tracers to advect and diffuse vector.         */
  trasco = wincon->trasc;
  uo = wincon->ultimate;
  if (wincon->ntbdys) {
    dto = windat->dttr;
    wincon->trasc = FFSL;
    wincon->ultimate = 0;
    if (master->tratio > 1.0 && !(wincon->means & TRANSPORT)) wincon->means |= TRANSPORT;
    if (advect_diffuse(window, windat, wincon)) return(1);
    vert_diffuse_3d(window, windat, wincon);
    if (wincon->means & TRANSPORT && windat->dttr) ffsl_trans_post(window, windat, wincon);
  }

  /* Set up the tracers to advect using VANLEER                      */
  if (wincon->ntbdy) {
    if (trasco & VANLEER)
      wincon->trasc = VANLEER;
    else if (trasco & QUICKEST)
      wincon->trasc = QUICKEST;
    else if (trasco & ORDER2)
      wincon->trasc = ORDER2;
    else if (trasco & ORDER2_UW)
      wincon->trasc = ORDER2_UW;
    else if (trasco & ORDER4)
      wincon->trasc = ORDER4;
    else if (trasco & ORDER1)
      wincon->trasc = ORDER1;
    wincon->ultimate = uo;
    windat->dttr = windat->dt;
    if (master->tratio > 1.0 && wincon->means & TRANSPORT) wincon->means &= ~TRANSPORT;
    if (advect_diffuse(window, windat, wincon)) return(1);
    /* Vertical diffusion for T/S is called from tracer_step_3d()    */
  }

  /* Reset the original tracers to advect and diffuse vector         */
  wincon->trasc = trasco;
  windat->dttr = dto;
  return(0);
}

/* END advect_diffuse_split()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the advective fluxes using the appropriate   */
/* advection scheme and invoke the ULTIMATE limiter if required.     */
/*-------------------------------------------------------------------*/
void advect(geometry_t *window,               /* Window geometry     */
	    window_t *windat,                 /* Window data         */
	    win_priv_t *wincon,               /* Window constants    */
            double *tr,                       /* Tracer array        */
            double *Fx,                       /* Face direction flux */
            double *Fz,                       /* Vertical flux       */
	    double *tr_mod, 
	    double *tr_mod_x, 
	    double *tr_mod_z, 
            double dtu,                       /* Time-step           */
	    int trasc                         /* Advection scheme    */
  )
{
  int c, cc;                    /* Cell coordinate / counter         */
  int e, ee;                    /* Face coordinate, counter          */
  double *w;                    /* Vertical velocity                 */
  double trp, trm;

  if (windat->dttr == 0.0)
    return;

  /*-----------------------------------------------------------------*/
  /* Get the tracer concentration at the cell faces                  */
  if (trasc == ORDER1)
    order1(window, windat, wincon, tr, Fx, Fz, dtu);
  else if (trasc & FFSL)
    ffsl_don(window, windat, wincon, tr, Fx, Fz, dtu);
  else {
    if (trasc == ORDER2)
      order2(window, windat, wincon, tr, Fx, Fz);
    else if (trasc == ORDER4)
      order4(window, windat, wincon, tr, Fx, Fz);
    else if (trasc & QUICKEST)
      quickest(window, windat, wincon, tr, Fx, Fz);
    else if (trasc & VANLEER)
      van_leer(window, windat, wincon, tr, Fx, Fz);
    else if (trasc == ORDER2_UW)
      order2_upwind(window, windat, wincon, tr, Fx, Fz);
    else if (trasc == LAGRANGE)
      semi_lagrange(window, windat, wincon, tr);
    else if (trasc == ORDER3US)
      order3us(window, windat, wincon, tr, Fx, Fz);
    else if (trasc == ORDER4US)
      order4us(window, windat, wincon, tr, Fx, Fz);

    /*---------------------------------------------------------------*/
    /* Invoke the ULTIMATE limiter. This procedure is split into     */
    /* doing the wet cells using velocity work arrays and auxiliary  */
    /* cells using the tracer work arrays so that ULTIMATE (which    */
    /* computationally expensive) is not invoked on cell faces       */
    /* which are subsequently set to zero as a consequence of being  */
    /* associated with zero velocity.                                */
    if (wincon->ultimate) {
      int cp, cm1, cm2;
      /* Wet cell faces (x direction) in this window                 */
      for (ee = 1; ee <= window->n3_e1; ee++) {
	e = window->w3_e1[ee];
	c = window->e2c[e][0];
	cm1 = window->e2c[e][1];
	if (wincon->trasc & HIORDER && !(window->eask[e] & (W_NOBC|W_TOBC))) {
	  trm = wincon->trm[e];
	  trp = wincon->trp[e];
	} else {
	  cp = window->e2c[window->ep[e]][0];
	  cm2 = window->e2c[window->em[e]][1];
	  trm = tr[cm2];
	  trp = tr[cp];
	}
	ultimate_filter(wincon->w6, Fx, tr, e, c, trp, cm1, trm); 
      }

      /* Vertical faces                                              */
      for (cc = 1; cc <= wincon->vc; cc++) {
	c = wincon->s1[cc];
	cp = window->zp1[c];
	cm1 = window->zm1[c];
	cm2 = window->zm1[cm1];
	ultimate_filter(wincon->w8, Fz, tr, c, c, tr[cp], cm1, tr[cm2]); 
      }
    }

    /*---------------------------------------------------------------*/
    /* Get the tracer flux through the face.                         */
    /* SIGMA : Don't multiply by the depth here since u1flux3d have  */
    /* already been depth multiplied.                                */
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      Fx[e] *= windat->u1flux3d[e];
    }

    if (wincon->means & TRANSPORT)
      w = windat->wm;
    else
      w = windat->w;

    for (cc = 1; cc <= wincon->vc; cc++) {
      c = wincon->s1[cc];
      Fz[c] *= (w[c] * dtu);
    }
  }
}

/* END advect()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the first order upwind advective fluxes      */
/*-------------------------------------------------------------------*/
void order1(geometry_t *window,               /* Window geometry     */
	    window_t *windat,                 /* Window data         */
	    win_priv_t *wincon,               /* Window constants    */
            double *tr,                       /* Tracer array        */
            double *Fx,                       /* Face direction flux */
            double *Fz,                       /* Vertical flux       */
            double dtu                        /* Time-step           */
	    )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, ee;                    /* Face coordinate, counter          */
  int zm1;                      /* Sparse cell at k-1                */
  double *w;                    /* Vertical velocity                 */

  /* Horizontal fluxes                                               */
  memset(Fx, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    /* Get the fluxes                                                */
    Fx[e] = wincon->w6[e] * tr[c2] + wincon->w7[e] * tr[c1];
  }

  /* Vertical fluxes                                                 */
  memset(Fz, 0, window->szc * sizeof(double));
  if (wincon->means & TRANSPORT)
    w = windat->wm;
  else
    w = windat->w;
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zm1 = window->zm1[c];
    Fz[c] = dtu * (0.5 * (w[c] + fabs(w[c])) * tr[zm1] +
                   0.5 * (w[c] - fabs(w[c])) * tr[c]);
  }
}

/* END order1()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* second order centered flux scheme. This scheme linearly           */
/* interpolates between cell centers to get the face value and is    */
/* formulated for non-uniform grids using the mean grid spacing      */
/* arrays.                                                           */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void order2(geometry_t *window,               /* Window geometry     */
	    window_t *windat,                 /* Window data         */
	    win_priv_t *wincon,               /* Window constants    */
            double *tr,                       /* Tracer array        */
            double *Fx,                       /* Face direction flux */
            double *Fz                        /* Vertical flux       */
	    )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, e2, ee;                /* Face coordinate, counter          */
  int zm1;                      /* Sparse cell at k-1                */

  /* Horizontal fluxes                                               */
  /* This calculation sets the tracer concentration on the cell face */
  /* (using 2nd order method).                                       */
  memset(Fx, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    Fx[e] = wincon->d3[e2] * (tr[c2] - tr[c1]) + tr[c1];
  }

  /* Vertical fluxes                                                 */
  memset(Fz, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zm1 = window->zm1[c];
    Fz[c] = wincon->w5[c] * (tr[c] - tr[zm1]) + tr[zm1];
  }
}

/* END order2()                                                      */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* 3rd order centered flux scheme for unstructured grids. This       */
/* scheme linearly interpolates between cell centers to get the face */
/* value and is formulated for non-uniform grids using the mean grid */
/* spacing arrays.                                                   */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void order3us(geometry_t *window,             /* Window geometry     */
	      window_t *windat,               /* Window data         */
	      win_priv_t *wincon,             /* Window constants    */
	      double *tr,                     /* Tracer array        */
	      double *Fx,                     /* Face direction flux */
	      double *Fz                      /* Vertical flux       */
	      )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, e2, ee;                /* Face coordinate, counter          */
  int zm1, zp1, zm2;            /* Sparse cell at k-1                */
  double e2p, e2m;
  double ot = 1.0 / 12.0;;
  double beta = 0.25;
  double sgn, h2;

  /* Horizontal fluxes                                               */
  /* This calculation sets the tracer concentration on the cell face */
  /* (using 2nd order method).                                       */
  memset(Fx, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    Fx[e] = wincon->d3[e2] * (tr[c2] - tr[c1]) + tr[c1];
  }

  /* Get the second derivative contribution. Note: a no-gradient is  */
  /* set over ghost cells, so the first and second derivatives are   */
  /* zero.                                                           */
  for (ee = 1; ee <= window->v3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    sgn = (windat->u[e] >= 0.0) ? 1.0 : -1.0;
    h2 = ot * window->h2au1[e2] * window->h2au1[e2];

    e2p = get_deriv2(window, windat, wincon, tr, e, 0);
    e2m = get_deriv2(window, windat, wincon, tr, e, 1);
    Fx[e] -= h2 * (e2p + e2m);
    Fx[e] += sgn * beta * h2 * (e2p - e2m);
    /*
    e2p = (sgn == 1.0) ? e2m : e2p;
    Fx[e] -= 2.0 * h2 * e2p;
    */
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes                                                 */
  memset(Fz, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    van_leer_do(Fz, tr, windat->w, wincon->w8, c, c, zp1, zm1, zm2);
  }
}

/* END order3us()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* 4th order centered flux scheme for unstructured grids. This       */
/* scheme linearly interpolates between cell centers to get the face */
/* value and is formulated for non-uniform grids using the mean grid */
/* spacing arrays.                                                   */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void order4us(geometry_t *window,             /* Window geometry     */
	      window_t *windat,               /* Window data         */
	      win_priv_t *wincon,             /* Window constants    */
	      double *tr,                     /* Tracer array        */
	      double *Fx,                     /* Face direction flux */
	      double *Fz                      /* Vertical flux       */
	      )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, e2, ee;                /* Face coordinate, counter          */
  int zm1, zp1, zm2;            /* Sparse cell at k-1                */
  double e2p, e2m;
  double ot = 1.0 / 12.0;

  /* Horizontal fluxes                                               */
  /* This calculation sets the tracer concentration on the cell face */
  /* (using 2nd order method).                                       */
  memset(Fx, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    Fx[e] = wincon->d3[e2] * (tr[c2] - tr[c1]) + tr[c1];
  }

  /* Get the second derivative contribution. Note: a no-gradient is  */
  /* set over ghost cells, so the first and second derivatives are   */
  /* zero.                                                           */
  for (ee = 1; ee <= window->v3_e1; ee++) {
    e = window->w3_e1[ee];
    e2 = window->m2de[e];
    e2p = get_deriv2(window, windat, wincon, tr, e, 0);
    e2m = get_deriv2(window, windat, wincon, tr, e, 1);
    Fx[e] -= ot * window->h2au1[e2] * window->h2au1[e2] * (e2p + e2m);
  }

  /* Vertical fluxes                                                 */
  memset(Fz, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    order4_do(Fz, tr, c, c, zp1, zm1, zm2);
  }
}

/* END order4us()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the fourth order advective fluxes.           */
/* Formulated for uniform grids.                                     */
/*-------------------------------------------------------------------*/
void order4(geometry_t *window,               /* Window geometry     */
	    window_t *windat,                 /* Window data         */
	    win_priv_t *wincon,               /* Window constants    */
            double *tr,                       /* Tracer array        */
            double *Fx,                       /* Face direction flux */
            double *Fz                        /* Vertical flux       */
	    )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, ee;                    /* Face coordinate, counter          */
  int j, cp, cm;
  int zp1, zm1, zm2;

  /* Horizontal fluxes                                               */
  memset(Fx, 0, window->sze * sizeof(double));
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    j = window->e2e[e][0];
    cm = window->c2c[j][c2];
    j = window->e2e[e][1];
    cp = window->c2c[j][c1];

    order4_do(Fx, tr, e, c1, cp, c2, cm);
  }

  /* Vertical fluxes                                                 */
  memset(Fz, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    order4_do(Fz, tr, c, c, zp1, zm1, zm2);
  }
}

/* END order4()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* second order centered flux scheme.                                */
/*-------------------------------------------------------------------*/
void order4_do(double *F,       /* Flux array                        */
	       double *tr,      /* df_variable_t to advect           */
	       int e,           /* Edge coordinate                   */
	       int cp1,         /* Cell in front of e                */
	       int cp2,         /* Cell in front of cp1              */
	       int cm1,         /* Cell behind e                     */
	       int cm2          /* Cell behind cm1                   */
	       )
{

  F[e] = (7.0 * (tr[cm1] + tr[cp1]) - (tr[cp2] + tr[cm2])) / 12.0;
}

/* END order4_do()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the third order advective fluxes using the   */
/* QUICKEST algorithm. Formulated for non-uniform grids.             */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*        wincon->w9 = vertical grid spacing between cell centers    */
/*-------------------------------------------------------------------*/
void quickest(geometry_t *window,             /* Window geometry     */
	      window_t *windat,               /* Window data         */
	      win_priv_t *wincon,             /* Window constants    */
	      double *tr,                     /* Tracer array        */
	      double *Fx,                     /* Face direction flux */
	      double *Fz                      /* Vertical flux       */
	      )
{
  int c, cc, c1, c2, cm, cp;    /* Cell coordinate / counter         */
  int e, em, ep, es, ee;        /* Face coordinate, counter          */
  int zm1, zp1, zm2;            /* Sparse cell at k-1                */
  int j;                        /* Edge direction                    */
  double face;                  /* Concentration value at the faces  */
  double grad;                  /* Concentration difference across a face */
  double curv;                  /* Concentration curvature across a face */
  double cx, cz;                /* Courant number in the x,y,z directions */
  double db, df;                /* Mean face centered grid spacings  */
  double trp, trm;              /* Tracer values at cp and cm        */

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.                                      */
  memset(Fx, 0, window->sze * sizeof(double));
  memset(Fz, 0, window->szc * sizeof(double));
  /* Horizontal fluxes                                               */
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    em = window->em[es];
    ep = window->ep[es];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    j = window->e2e[e][0];
    cm = window->c2c[j][c2];
    j = window->e2e[e][1];
    cp = window->c2c[j][c1];

    if (wincon->trasc & HIORDER && !(window->eask[e] & (W_NOBC|W_TOBC))) {
      trp = wincon->trp[e];
      trm = wincon->trm[e];
    } else {
      trp = tr[cp];
      trm = tr[cm];
    }

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the left face.  */
    /* The face value is formulated for a non-uniform grid by        */
    /* linearly interpolating between cell centers and evaluating on */
    /* the face.                                                     */
    cx = wincon->w6[e];
    db = 2.0 / (window->h2au1[em] + window->h2au1[es]);
    df = 2.0 / (window->h2au1[ep] + window->h2au1[es]);
    face = wincon->d3[es] * (tr[c2] - tr[c1]) + tr[c1];
    grad = tr[c1] - tr[c2];
    curv = 0.0;

    /* Note: the curvature is formulated for a non-uniform grid (see */
    /* Kowalik & Murty (1993), p40) and multiplied by (h2au1[i] *    */
    /* h2au1[i-1]) to give it dimensions of [tracer].                */
    if (cx > 0.0)
      curv =
        db * (window->h2au1[em] * tr[c1] + window->h2au1[es] * trm) -
        2.0 * tr[c2];
    else if (cx < 0.0)
      curv =
        df * (window->h2au1[es] * trp +
              window->h2au1[window->m2de[ep]] * tr[c2]) - 2.0 * tr[c1];
    Fx[e] = face - 0.5 * cx * grad - (1.0 - cx * cx) * curv / 6.0;
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the lower face  */
    cz = wincon->w8[c];
    db = 2.0 / (wincon->w9[zm1] + wincon->w9[c]);
    df = 2.0 / (wincon->w9[c] + wincon->w9[zp1]);
    face = wincon->w5[c] * (tr[c] - tr[zm1]) + tr[zm1];
    grad = tr[c] - tr[zm1];
    curv = 0.0;
    if (cz > 0.0)
      curv =
        db * (wincon->w9[zm1] * tr[c] + wincon->w9[c] * tr[zm2]) -
        2.0 * tr[zm1];
    else if (cz < 0.0)
      curv =
        df * (wincon->w9[c] * tr[zp1] + wincon->w9[zp1] * tr[zm1]) -
        2.0 * tr[c];
    Fz[c] = face - 0.5 * cz * grad - (1.0 - cz * cz) * curv / 6.0;
  }
}

/* END quickest()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the third order advective fluxes using the   */
/* QUICKEST algorithm. Formulated for uniform grids.                 */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void quickest_uniform(geometry_t *window,     /* Window geometry     */
		      window_t *windat,       /* Window data         */
		      win_priv_t *wincon,     /* Window constants    */
		      double *tr,             /* Tracer array        */
		      double *Fx,             /* Face direction flux */
		      double *Fz              /* Vertical flux       */
		      )
{
  int c, cc, c1, c2, cm, cp;    /* Cell coordinate / counter         */
  int e, em, ep, es, ee;        /* Face coordinate, counter          */
  int zm1, zp1, zm2;            /* Sparse cell at k-1                */
  int j;                        /* Edge direction                    */
  double face;                  /* Concentration value at the faces  */
  double grad;                  /* Concentration difference across a face */
  double curv;                  /* Concentration curvature across a face */
  double cx, cy, cz;            /* Courant number in the x,y,z directions */

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.                                      */
  memset(Fx, 0, window->sze * sizeof(double));
  memset(Fz, 0, window->szc * sizeof(double));
  /* Horizontal fluxes                                               */
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    em = window->em[es];
    ep = window->ep[es];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    j = window->e2e[e][0];
    cm = window->c2c[j][c2];
    j = window->e2e[e][1];
    cp = window->c2c[j][c1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the left face.  */
    /* The face value is formulated for a non-uniform grid by        */
    /* linearly interpolating between cell centers and evaluating on */
    /* the face.                                                     */
    cx = wincon->w6[e];
    face = 0.5 * (tr[c1] + tr[c2]);
    grad = tr[c1] - tr[c2];
    curv = 0.0;
    if (cx > 0.0)
      curv = tr[c1] - 2.0 * tr[c2] + tr[cm];
    else if (cx < 0.0)
      curv = tr[cp] - 2.0 * tr[c1] + tr[c2];
    Fx[c] = face - 0.5 * cx * grad - (1.0 - cx * cx) * curv / 6.0;
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes                                                 */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the lower face  */
    cz = wincon->w8[c];
    face = 0.5 * (tr[c] + tr[zm1]);
    grad = tr[c] - tr[zm1];
    curv = 0.0;
    if (cz > 0.0)
      curv = tr[c] - 2.0 * tr[zm1] + tr[zm2];
    else if (cz < 0.0)
      curv = tr[zp1] - 2.0 * tr[c] + tr[zm1];
    Fz[c] = face - 0.5 * cz * grad - (1.0 - cz * cz) * curv / 6.0;
  }
}

/* END quickest_uniform()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the higher order Van Leer upwind advective   */
/* fluxes.                                                           */
/*-------------------------------------------------------------------*/
void van_leer(geometry_t *window,             /* Window geometry     */
	      window_t *windat,               /* Window data         */
	      win_priv_t *wincon,             /* Window constants    */
	      double *tr,                     /* Tracer array        */
	      double *Fx,                     /* Face direction flux */
	      double *Fz                      /* Vertical flux       */
	      )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, ee;                    /* Face coordinate, counter          */
  int j, cp, cm;
  int zp1, zm1, zm2;
  int es;

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.                                      */
  memset(Fx, 0, window->sze * sizeof(double));
  memset(Fz, 0, window->szc * sizeof(double));
  /* Horizontal fluxes                                               */
  for (ee = 1; ee <= window->n3_e1; ee++) {
    double trp, trm;
    e = window->w3_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    j = window->e2e[e][0];
    cm = window->c2c[j][c2];
    j = window->e2e[e][1];
    cp = window->c2c[j][c1];

    if (wincon->trasc & HIORDER && !(window->eask[e] & (W_NOBC|W_TOBC))) {
      trp = wincon->trp[e];
      trm = wincon->trm[e];
    } else {
      trp = tr[cp];
      trm = tr[cm];
    }

    van_leer_tr(Fx, tr, windat->u1, wincon->w6, e, tr[c1], trp, tr[c2], trm);

    /*van_leer_do(Fx, tr, windat->u1, wincon->w6, e, c1, cp, c2, cm);*/
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes                                                 */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    van_leer_do(Fz, tr, windat->w, wincon->w8, c, c, zp1, zm1, zm2);
  }
}

/* END van_leer()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using      */
/* VanLeer's method.                                                 */
/*-------------------------------------------------------------------*/
void van_leer_do(double *F,     /* Flux array                        */
                 double *tr,    /* df_variable_t to advect           */
                 double *vel,   /* Velocity at the cell center/face  */
                 double *cn,    /* Courant number for this dimension */
		 int e,         /* Edge coordinate                   */
		 int cp1,       /* Cell in front of e                */
		 int cp2,       /* Cell in front of cp1              */
		 int cm1,       /* Cell behind e                     */
		 int cm2        /* Cell behind cm1                   */
  )
{
  double df;

  if (vel[e] > 0.0) {
    df = (tr[cm1] - tr[cm2]) * (tr[cp1] - tr[cm1]);
    if (df > 0.0)
      df = 2.0 * df / (tr[cp1] - tr[cm2]);
    else
      df = 0.0;
    F[e] = tr[cm1] + 0.5 * (1.0 - cn[e]) * df;
  } else if (vel[e] < 0.0) {
    df = (tr[cp1] - tr[cm1]) * (tr[cp2] - tr[cp1]);
    
    if (df > 0.0)
      df = 2.0 * df / (tr[cp2] - tr[cm1]);
    else
      df = 0.0;
    F[e] = tr[cp1] - 0.5 * (1.0 + cn[e]) * df;
  }
}


void van_leer_tr(double *F,     /* Flux array                        */
                 double *tr,    /* df_variable_t to advect           */
                 double *vel,   /* Velocity at the cell center/face  */
                 double *cn,    /* Courant number for this dimension */
		 int e,         /* Edge coordinate                   */
		 double tp1,    /* Tracer value in front of e        */
		 double tp2,    /* Tracer value in front of tp1      */
		 double tm1,    /* Tracer value behind e             */
		 double tm2     /* Tracer value behind tm1           */
  )
{
  double df;

  if (vel[e] > 0.0) {
    df = (tm1 - tm2) * (tp1 - tm1);
    if (df > 0.0)
      df = 2.0 * df / (tp1 - tm2);
    else
      df = 0.0;
    F[e] = tm1 + 0.5 * (1.0 - cn[e]) * df;
  } else if (vel[e] < 0.0) {
    df = (tp1 - tm1) * (tp2 - tp1);
    
    if (df > 0.0)
      df = 2.0 * df / (tp2 - tm1);
    else
      df = 0.0;
    F[e] = tp1 - 0.5 * (1.0 + cn[e]) * df;
  }
}

/* END van_leer_do()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the 2nd order upwind advective fluxes. This  */
/* scheme is stable for |cn| <= 2 (Leonard, 1994, Comput. Methods    */
/* Appl. Mech. Engrg.)                                               */
/*-------------------------------------------------------------------*/
void order2_upwind(geometry_t *window,        /* Window geometry     */
		   window_t *windat,          /* Window data         */
		   win_priv_t *wincon,        /* Window constants    */
		   double *tr,                /* Tracer array        */
		   double *Fx,                /* Face direction flux */
		   double *Fz                 /* Vertical flux       */
		   )
{
  int c, cc, c1, c2;            /* Cell coordinate / counter         */
  int e, ee;                    /* Face coordinate, counter          */
  int j, cp, cm;
  int zp1, zm1, zm2;

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.                                      */
  memset(Fx, 0, window->sze * sizeof(double));
  memset(Fz, 0, window->szc * sizeof(double));
  /* Horizontal fluxes                                               */
  for (ee = 1; ee <= window->n3_e1; ee++) {
    e = window->w3_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    j = window->e2e[e][0];
    cm = window->c2c[j][c2];
    j = window->e2e[e][1];
    cp = window->c2c[j][c1];

    order2_upwind_do(Fx, tr, windat->u1, wincon->w6, e, c1, cp, c2, cm);
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes                                                 */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    order2_upwind_do(Fz, tr, windat->w, wincon->w8, c, c, zp1, zm1, zm2);
  }
}

/* END order2_upwind()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using      */
/* the second order upwind method.                                   */
/*-------------------------------------------------------------------*/
void order2_upwind_do(double *F,   /* Flux array                     */
		      double *tr,  /* df_variable_t to advect        */
		      double *vel, /* Velocity at the cell center/face  */
		      double *cn,  /* Courant number for this dimension */
		      int e,       /* Edge coordinate                */
		      int cp1,     /* Cell in frint of e             */
		      int cp2,     /* Cell in frint of ecp1          */
		      int cm1,     /* Cell behind e                  */
		      int cm2      /* Cell behind cm1                */
		      )
{
  double df;

  if (vel[e] > 0.0) {
    F[e] = 0.5 * (3.0 - cn[e]) * tr[cm1] -
	0.5 * (1.0 - cn[e]) * tr[cm2];
    } else if (vel[e] < 0.0) {
      F[e] = 0.5 * (3.0 - cn[e]) * tr[cp1] -
	0.5 * (1.0 - cn[e]) * tr[cp2];
  } else
      F[e] = 0.0;
}

/* End order2_upwind_do()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to initialize variables for the full 3D flux-form         */
/* semi-lagrangian advective scheme, based on Leonard et al. (1996). */
/* Take the cell-centered velocities in each direction, calculate a  */
/* trajectory distance and identify the cell at the head of the      */ 
/* trajectory. Repeat for the cell faces.                            */
/*-------------------------------------------------------------------*/
void prep_ff_sl(geometry_t *window, /* Window structure */
		window_t *windat,   /* Window data structure */
		win_priv_t *wincon, /* Window geometry / constants */
		double dt           /* sub-time step */
  )
{
  double *crfzf; /* Fractional factors of total trajectories for faces */
  double *crfzc; /* Fractional factors of total trajectories for cell centres */
  int *clzf;     /* Cell counter at head of trajectory, for faces    */
  int *clzc;     /* Cell counter at head of trajectory, for cell centres */
  double *nw;            /* Cell-centered velocities */
  double *w;             /* Velocities at faces */
  int c, cc, cp1;

  /* Set pointers. */
  clzf = wincon->clzf;
  clzc = wincon->clzc;
  crfzf = wincon->crfzf;
  crfzc = wincon->crfzc;
  nw = wincon->nw;
  w = (wincon->means & TRANSPORT) ? windat->wm : windat->w;

  /* Initialise */
  memset(crfzf, 0, window->sgsiz * sizeof(double));
  memset(crfzc, 0, window->sgsiz * sizeof(double));
  memset(clzf, 0, window->sgsiz * sizeof(int));
  memset(clzc, 0, window->sgsiz * sizeof(int));

  /* Calculate the cell-averaged velocities in the z direction, and  */
  /* identify the cell from which to start the backward trajectory   */
  /* (which is dependent on the velocity sense).                     */
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    cp1 = window->zp1[c];
    clzc[c] = (nw[c] < 0.0) ? window->zp1[c] : window->w3_t[cc];
  }
  
  /* Get trajectory distances and origins for the transverse terms  */
  /* based on the cell-centered velocities.                         */
  ff_sl_do_vert(nw, dt, wincon->w9, 1, window->a3_t,
       window->w3_t, window->zp1, window->zm1, clzc, crfzc);

  /* Get trajectory distance and origin for advection using */
  /* cell face velocities for each direction.               */
  /* z direction */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    clzf[c] = (w[c] > 0.0) ? window->zm1[c] : wincon->s1[cc];
  }
  ff_sl_do_vert(w, dt, wincon->dz, 1, wincon->vc, wincon->s1, window->zp1,
     window->zm1, clzf, crfzf);
}

/* END prep_ff_sl()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* First order routine to track the backwards trajectory to the      */
/* source, and note the start sparse location and fractional         */
/* component of the trajectory.                                      */
/*-------------------------------------------------------------------*/
void ff_sl_do_vert(double *vel,   /* Velocity at the cell face */
		   double dt,     /* time step */
		   double *h,     /* cell width/height at face/centre */
		   int ss,        /* Start sparse coordinate */
		   int se,        /* End sparse coordinate */
		   int *sdo,      /* Points to process array */
		   int *fmap,     /* Forward map */
		   int *bmap,     /* Backward map */
		   int *cl,        /* Index of trajectory source cell */
		   double *crf     /* Fractional component of trajectory */
  )
{
  int c, cc, cp, cpb, c2;
  double ltraj, lrem, dz;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    if (vel[c] == 0.0)
      continue;
    ltraj = fabs(vel[c]) * dt;
    lrem = ltraj;
    cp = cpb = cl[c];
    dz = fabs((h[cp] > 0) ? h[cp] : h[c]);
    while (dz < lrem) {
      lrem  -= dz;
      cp = (vel[c] < 0.0) ? fmap[cp] : bmap[cp];
      if (cpb == cp) {
	lrem = fmod(lrem, dz);
	break;
      }
      if (h[cp] > 0)
        dz = fabs(h[cp]);
      cpb = cp;
    }

    crf[c] = max(min(lrem / dz, 1.0), 0.0);
    cl[c] = cp;
 }
}

/* END ff_sl_do_vert()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using      */
/* VanLeer's method. Used with the FFSL advection scheme.            */
/*-------------------------------------------------------------------*/
void ff_sl_van_leer(double *F,  /* Flux array */
		    double *tr, /* df_variable_t to advect */
		    double *vel,/* Velocity at the cell center/face */
		    double *crf,/* Fractional component of lagrange time step */
		    int *cl,    /* Upstream coordinate for long time steps */ 
		    int ss,     /* Start sparse coordinate */
		    int se,     /* End sparse coordinate */
		    int *sdo,   /* Points to process array */
		    int *fmap,  /* Forward map */
		    int *bmap   /* Backward map */
		   )
{
  int c, cc, cp;
  int bm1, bm2;
  double df;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    if (crf[c] == 0.0)
      continue;
    cp = cl[c];
    if (vel[c] > 0.0)
      cp = fmap[cp];
    bm1 = bmap[cp];
    bm2 = bmap[bm1];

    if (vel[c] > 0.0) {
      df = (tr[bm1] - tr[bm2]) * (tr[cp] - tr[bm1]);
      if (df > 0.0)
        df = 2.0 * df / (tr[cp] - tr[bm2]);
      else
        df = 0.0;
      F[c] = tr[bm1] + 0.5 * (1.0 - crf[c]) * df;
    } else if (vel[c] < 0.0) {
      df = (tr[cp] - tr[bm1]) * (tr[fmap[cp]] - tr[cp]);
      if (df > 0.0)
        df = 2.0 * df / (tr[fmap[cp]] - tr[bm1]);
      else
        df = 0.0;
      F[c] = tr[cp] - 0.5 * (1.0 - crf[c]) * df;
    }
  }
}

/* END ff_sl_van_leer()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Calculates updated tracer concentration of surface cells          */
/*-------------------------------------------------------------------*/
double surf_conc(geometry_t *window, /* Window geometry              */
		 window_t *windat,   /* Window data                  */
		 win_priv_t *wincon, /* Window constants             */
                 double *tr,         /* Tracer array                 */
                 double osubeta, /* Old surface elevation            */
                 double subeta,  /* Surface elevation at sub-timestep */
                 int c,          /* 3D coord of lowest surface cell  */
                 int c2,         /* 2D sparse coordinate             */
                 int cc,         /* 2D sparse counter                */
                 double fcbot,   /* Tracer flux through the bottom   */
                 double *dtracer /* Advective terms x directions     */
		                 /* [tracer][m]^2.                   */
  )

{
  int cb;                /* Sparse coordinate of the bottom          */
  int co;                /* Old surface coordinate                   */
  int zp1;               /* cell above sparse coordinate c           */
  double dz;             /* Distance between layer faces             */
  double svol;           /* Volume of the surface cells              */
  double surftr;         /* Tracer concentration in the surface cell */
  double watertop;       /* z coordinate of water top in cell        */
  double waterbot;       /* z coordinate of water bottom in cell     */
  double hf, vf;         /* Diagnostics for continuity               */
  double dtr, cnt;       /* Diagnostics for continuity               */
  int dc = 0;            /* Diagnostic continuity surface coordinate */
  int dw = 1;            /* Window for diagnostic continuity         */
  double div = 0.0;
  double diff;
  double dtu = wincon->b1;
  int e, j;
  double trs = 0.0;
  /*-----------------------------------------------------------------*/
  /* Get the thickness and k index of the surface cell taking into   */
  /* account water columns that are only one cell deep. This is the  */
  /* volume between the new (updated) elevation and ks. Note that ks */
  /* is the lower vertical index of the surface at the old and       */
  /* updated time steps, i.e.                                        */
  /* ks = max(window->sur_t[cc],window->nsur_t[cc])                  */
  /* dz=windat->eta[c2]-window->gridz[c];                            */
  /* cb=window->bot_t[cc];                                           */
  cb = wincon->i1[cc];     /* Bottom coordinate for cells to process */
  co = wincon->i3[cc];     /* Old surface coodinate                  */
  /*osubeta = wincon->d1[c2];*/
  waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
  dz = subeta - waterbot;
  svol = dz * window->cellarea[c2];

  /*-----------------------------------------------------------------*/
  /* Calculate the total amount of tracer which will be present in   */
  /* the surface volume. First get the tracer due to the vertical    */
  /* flux through the bottom of this volume.                         */
  surftr = -fcbot;
  hf = vf = dtr = cnt = 0.0;
  vf = -windat->w[c] * dtu * window->cellarea[c2];

  /*-----------------------------------------------------------------*/
  /* Loop through the cells from the current layer until the layer   */
  /* below the old surface and add the tracer*volume of this cell +  */
  /* the change in tracer due to horizontal divergence. If the       */
  /* elevation has risen this is only the old surface cell. If       */
  /* elevation has dropped these are the cells from that below the   */
  /* new surface to that below the old surface. These cells were wet */
  /* at the previous time step and therefore contained tracer. The   */
  /* last cell volume calculated here is done using oldeta so that   */
  /* the change in volume correctly corresponds to the fluxes        */
  /* calculated previously.                                          */

  /* First increment the total tracer from the current layer to the  */
  /* layer below that containing the old surface. If the elevation   */
  /* has risen this loop is skipped.                                 */
  while (c != co) {
    zp1 = window->zp1[c];

    /* Get the depth of the top of the layer for the cell            */
    watertop = window->gridz[zp1];

    /* Increment the total amount of tracer with the amount in this  */
    /* cell plus the horizontal divergence through this cell.        */
    surftr +=
      (tr[c] * window->cellarea[c2] * (watertop - waterbot) - dtracer[c]);

    /* Diagnostics                                                   */
    if(window->wn == dw && c2 == dc) {
      hf = 0.0;
      for (j = 1; j <= window->npe[c2]; j++) {
        e = window->c2e[j][c];
        hf += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
      }  
      dtr += dtracer[c];
      trs += tr[c] * window->cellarea[c2] * (watertop - waterbot);
      printf("rise %d %d %f : %f %f\n",c,window->s2k[c],surftr,hf,dtr/tr[c]);
    }

    /* Get the depth of the bottom layer of the next cell up         */
    c = zp1;
    waterbot = window->gridz[c];
  }

  /* Now increment the total tracer for the layer containing the old */
  /* surface.                                                        */
  /* The old surface elevation is the upper bound for cell thickness */
  /* Note : d1=osubeta=oldeta                                        */
  watertop = osubeta;

  /* Get the depth of the bottom layer of the cell : the bottom z    */
  /* coordinate if the cell lies on the bottom.                      */
  waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
  /* Increment the total amount of tracer with the amount in this    */
  /* cell plus the horizontal divergence through this cell.          */
  surftr +=
    (tr[c] * window->cellarea[c2] * (watertop - waterbot) - dtracer[c]);

  /* Diagnostics for tracers with constant vaues. The fluxes         */
  /* calculated directly and those derived from dtracer divided      */
  /* by the tracer value should be equivalent.                       */
  if(window->wn == dw && c2 == dc) {
    hf = 0.0;
    for (j = 1; j <= window->npe[c2]; j++) {
      e = window->c2e[j][c];
      hf += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
    }
    dtr += dtracer[c];
    cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
    printf("surf_conc %f %d %d %e : %f %f : %f %f : %f\n",windat->t/86400,
	   c,window->s2k[c],cnt,hf,dtr/tr[c],vf,
	   fcbot/tr[c],surftr/svol);
  }

  /*-----------------------------------------------------------------*/
  /* Loop through the cells from above the old surface elevation to  */
  /* the top of the vertical grid. If the surface has risen then     */
  /* these cells now contain water and horizontal divergence may     */
  /* alter the total amount of tracer. If elevation has dropped then */
  /* there is no water (and hence no divergence) in these cells.     */
  if (c != c2) {
    do {
      c = window->zp1[c];
      /* Increment the horizontal divergence through this cell       */
      /*
      if(window->wn == dw && c2 == dc) {
        dtr = 0.0;
	for (j = 1; j <= window->npe[c2]; j++) {
	  e = window->c2e[j][c];
	  dtr += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
	}
	hf += dtr;
	cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
	printf("to-top %d %d %e : %f %f : %e\n",window->s2k[c],c,cnt,dtr,dtracer[c],surftr/svol);
      }
      */
      surftr -= (dtracer[c]);
    } while (c != c2);
  }

  /*-----------------------------------------------------------------*/
  /* Divide by surface volume to get concentration.                  */
  /* Note that tracers with non-zero settling velocities can go      */
  /* negative in thin water layers, so set to zero in that case.     */
  if (svol == 0.0 || (surftr < 0.0 && dz <= wincon->hmin)) {
    surftr = tr[c];
    surftr = 0.0;
  }
  else { 
   surftr /= svol;
  }
  return (surftr);
}

/* END surf_conc()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates updated tracer concentration of surface cells          */
/*-------------------------------------------------------------------*/
double surf_concrk(geometry_t *window, /* Window geometry            */
		   window_t *windat,   /* Window data                */
		   win_priv_t *wincon, /* Window constants           */
		   double *tr,     /* Tracer array                   */
		   double osubeta, /* Old surface elevation          */
		   double subeta,  /* Elevation at sub-timestep      */
		   int c,          /* Index of lowest surface cell   */
		   int c2,         /* 2D sparse coordinate           */
		   int cc,         /* 2D sparse counter              */
		   double *grad,   /* Time tendency                  */
		   double *ovol,   /* Volume at t                    */
		   double *nvol,   /* Volume at t+1                  */
		   double fcbot,   /* Tracer flux through the bottom */
		   double *dtracer /* Advective terms x directions   */
		                   /* [tracer][m]^2.                 */
  )

{
  int e, j;              /* Edge counters                            */
  int cs = c;            /* Surface coordinate                       */
  int cb;                /* Sparse coordinate of the bottom          */
  int co;                /* Old surface coordinate                   */
  int zp1;               /* cell above sparse coordinate c           */
  double dz;             /* Distance between layer faces             */
  double surftr;         /* Tracer concentration in the surface cell */
  double watertop;       /* z coordinate of water top in cell        */
  double waterbot;       /* z coordinate of water bottom in cell     */
  double vol;            /* Cell volume                              */
  double hf, vf;         /* Diagnostics for continuity               */
  double dtr, cnt;       /* Diagnostics for continuity               */
  int dc = 0;            /* Diagnostic continuity surface coordinate */
  int dw = 1;            /* Window for diagnostic continuity         */
  double div = 0.0;
  double diff;
  double dtu = wincon->b1;

  /*-----------------------------------------------------------------*/
  /* Get the thickness and k index of the surface cell taking into   */
  /* account water columns that are only one cell deep. This is the  */
  /* volume between the new (updated) elevation and ks. Note that ks */
  /* is the lower vertical index of the surface at the old and       */
  /* updated time steps, i.e.                                        */
  /* ks = max(window->sur_t[cc],window->nsur_t[cc])                  */
  /* dz=windat->eta[c2]-window->gridz[c];                            */
  /* cb=window->bot_t[cc];                                           */
  cb = wincon->i1[cc];     /* Bottom coordinate for cells to process */
  co = wincon->i3[cc];     /* Old surface coodinate                  */
  /*osubeta = wincon->d1[c2];*/
  waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
  dz = subeta - waterbot;
  *nvol = dz * window->cellarea[c2];

  /*-----------------------------------------------------------------*/
  /* Calculate the total amount of tracer which will be present in   */
  /* the surface volume. First get the tracer due to the vertical    */
  /* flux through the bottom of this volume.                         */
  surftr = *ovol = 0.0;
  grad[cs] = -fcbot;
  hf = vf = dtr = cnt = 0.0;
  vf = -windat->w[c] * dtu * window->cellarea[c2];

  /*-----------------------------------------------------------------*/
  /* Loop through the cells from the current layer until the layer   */
  /* below the old surface and add the tracer*volume of this cell +  */
  /* the change in tracer due to horizontal divergence. If the       */
  /* elevation has risen this is only the old surface cell. If       */
  /* elevation has dropped these are the cells from that below the   */
  /* new surface to that below the old surface. These cells were wet */
  /* at the previous time step and therefore contained tracer. The   */
  /* last cell volume calculated here is done using oldeta so that   */
  /* the change in volume correctly corresponds to the fluxes        */
  /* calculated previously.                                          */

  /* First increment the total tracer from the current layer to the  */
  /* layer below that containing the old surface. If the elevation   */
  /* has risen this loop is skipped.                                 */
  while (c != co) {
    zp1 = window->zp1[c];

    /* Get the depth of the top of the layer for the cell            */
    watertop = window->gridz[zp1];

    /* Increment the total amount of tracer with the amount in this  */
    /* cell plus the horizontal divergence through this cell.        */
    vol = (watertop - waterbot) * window->cellarea[c2];
    *ovol += vol;
    surftr += (tr[c] * vol);
    grad[cs] -= dtracer[c];

    /* Diagnostics                                                   */
    if(window->wn == dw && c2 == dc) {
      /*hf = 0.0;*/
      for (j = 1; j <= window->npe[c2]; j++) {
        e = window->c2e[j][c];
        hf += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
      }  
      dtr += dtracer[c];
      printf("rise %d %d %f : %f %f\n",c,window->s2k[c],surftr,hf,dtr/tr[c]);
    }

    /* Get the depth of the bottom layer of the next cell up         */
    c = zp1;
    waterbot = window->gridz[c];
  }

  /* Now increment the total tracer for the layer containing the old */
  /* surface.                                                        */
  /* The old surface elevation is the upper bound for cell thickness */
  /* Note : d1=osubeta=oldeta                                        */
  watertop = osubeta;

  /* Get the depth of the bottom layer of the cell : the bottom z    */
  /* coordinate if the cell lies on the bottom.                      */
  waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
  /* Increment the total amount of tracer with the amount in this    */
  /* cell plus the horizontal divergence through this cell.          */
  /*surftr += (tr[c] * (watertop - waterbot));*/
  vol = (watertop - waterbot) * window->cellarea[c2];
  *ovol += vol;
  surftr += tr[c] * vol;
  grad[cs] -= dtracer[c];

  /* Diagnostics for tracers with constant vaues. The fluxes         */
  /* calculated directly and those derived from dtracer divided      */
  /* by the tracer value should be equivalent.                       */
  if(window->wn == dw && c2 == dc) {
    /*hf = 0.0;*/
    for (j = 1; j <= window->npe[c2]; j++) {
      e = window->c2e[j][c];
      hf += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
    }
    dtr += dtracer[c];
    cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
    printf("surf_conc %f %d %d %e : %f %f : %f %f : %f\n",windat->t/86400,
	   c,window->s2k[c],cnt,hf,dtr/tr[c],vf,
	   fcbot/tr[c],(surftr+grad[c])/(*nvol));
  }

  /*-----------------------------------------------------------------*/
  /* Loop through the cells from above the old surface elevation to  */
  /* the top of the vertical grid. If the surface has risen then     */
  /* these cells now contain water and horizontal divergence may     */
  /* alter the total amount of tracer. If elevation has dropped then */
  /* there is no water (and hence no divergence) in these cells.     */
  if (c != c2) {
    do {
      c = window->zp1[c];
      /* Increment the horizontal divergence through this cell       */
      /*
      if(window->wn == dw && c2 == dc) {
        dtr = 0.0;
	for (j = 1; j <= window->npe[c2]; j++) {
	  e = window->c2e[j][c];
	  dtr += (window->eSc[j][c2] * windat->u1flux3d[e] * dtu);
	}
	hf += dtr;
	cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
	printf("to-top %d %d %e : %f %f : %e\n",window->s2k[c],c,cnt,dtr,dtracer[c],surftr/(*nvol));
      }
      */
      grad[cs] -= dtracer[c];
    } while (c != c2);
  }

  /*-----------------------------------------------------------------*/
  /* Divide by surface volume to get concentration.                  */
  /* Note that tracers with non-zero settling velocities can go      */
  /* negative in thin water layers, so set to zero in that case.     */
  if (*nvol == 0.0 || (surftr < 0.0 && dz <= wincon->hmin)) {
    surftr = tr[cs];
    surftr = 0.0;
  }
  /* Mass of tracer at time t (e.g. kg) and flux (e.g. kg/s) are     */
  /* returned, and converted to concentrations in runge_kutta().     */
  return (surftr);
}

/* END surf_concrk()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check for tracer range violations and print a message. */
/* Min / max violation for several special cases are clipped.        */
/*-------------------------------------------------------------------*/
void tr_bounds(geometry_t *window,        /* Window geometry         */
	       window_t *windat,          /* Window data             */
	       win_priv_t *wincon         /* Window constants        */
	       )
{
  int c, c3, cc;                /* Counters                          */
  int cs, cb;                   /* Surface and bottom coordinates    */
  int co;                       /* Old surface coordinate            */
  int n, nn;                    /* Tracer number                     */
  double surftr;                /* New surface tracer value          */
  double eps = 1e-5;            /* Minimum value for clipping        */
  double *min = wincon->mintr;
  double *max = wincon->maxtr;

  /* Average the tracer over intersecting OBCs if required           */
  for (cc = 0; cc < window->nobc; cc++) {
    open_bdrys_t *open = window->open[cc]; 
    if (open->options & OP_OBCCM) {
      for (nn = 0; nn < wincon->ntbdy; nn++) {
        n = wincon->tbdy[nn];
	average_OBC_corner(window, open, windat->tr_wc[n], 1);
      }
    }
  }
  if (wincon->mode2d) return;

  /* In thin layers the difference in diffusive or advective fluxes  */
  /* between opposite faces can be large (eg. if dzface(i) is thin   */
  /* and dzface(i+1) is not), and dividing this by the (small) cell  */
  /* volume can result in large rates of change of the tracer, which */
  /* can violate the maximum/minimum allowed concentration.          */
  /* Therefore, if the layer is thin and a max / min violation       */
  /* occurs in the surface layer only, then set the surface          */
  /* concentration to that of the layer below.                       */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = c3 = wincon->s1[cc];
    co = wincon->i3[cc];
    cb = wincon->i1[cc];
    cs = window->m2d[c];
    if (c < cb && wincon->dz[co] < wincon->hmin) {
      for (nn = 0; nn < wincon->ntbdy; nn++) {
        n = wincon->tbdy[nn];
        if (windat->tr_wc[n][c] < min[n] || windat->tr_wc[n][c] > max[n]) {
          surftr = windat->tr_wc[n][window->zm1[c]];
          windat->tr_wc[n][c3] = surftr;
          while (c3 != window->zp1[c3]) {
            c3 = window->zp1[c3];
            windat->tr_wc[n][c3] = surftr;
          }
        }
      }
    }
  }

#if defined(HAVE_OMP)
#pragma omp parallel for private(cc,c,nn,n)
#endif

  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      n = wincon->tbdy[nn];
      /* If a minimum is set to zero and the tracer undershoots      */
      /* below this then set to zero to keep positive definite.      */
      if (min[n] <= eps && windat->tr_wc[n][c] < min[n])
	windat->tr_wc[n][c] = min[n];
      
      /* Print a warning for any other min / max violations          */
      if (windat->tr_wc[n][c] < min[n] || windat->tr_wc[n][c] > max[n]) {
	hd_warn
	  ("Tracer %s outside range %3.1e to %3.1e : %e @ (%d %d %d) at %f\n",
	   wincon->trname[n], min[n], max[n], windat->tr_wc[n][c],
	   window->s2i[c], window->s2j[c], window->s2k[c], windat->t / 86400.0);
      }
#if defined(HAVE_ECOLOGY_MODULE)
      /*
	if (windat->tr_wc[n][c] > 0.0 && windat->tr_wc[n][c] < TINY) {
	hd_warn
	("Tracer %s too small : %e @ (%d %d %d) at %f\n", 
	wincon->trname[n], windat->tr_wc[n][c],
	window->s2i[c], window->s2j[c], window->s2k[c], windat->t / 86400.0);
        
	}
      */
#endif
    }
  }

  /* Set all layers above the surface to surface layer               */
  if (!wincon->sigma) {
#if defined(HAVE_OMP)
#pragma omp parallel for private(cc,c,cs,n,nn)
#endif
    for (cc = 1; cc <= wincon->vcs; cc++) {
      c = cs = wincon->s1[cc];
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	for (nn = 0; nn < wincon->ntbdy; nn++) {
	  n = wincon->tbdy[nn];
	  windat->tr_wc[n][c] = windat->tr_wc[n][cs];
	}
      }
    }
  }
}

/* END tr_bounds()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets a no gradient above the surface. This assumes the tracer     */
/* step is complete, and uses cells to process corresponding to the  */
/* updated elevation.                                                */
/*-------------------------------------------------------------------*/
void tr_set_surf(geometry_t *window,        /* Window geometry       */
		 window_t *windat,          /* Window data           */
		 win_priv_t *wincon         /* Window constants      */
		 )
{
  int c, cs, cc, n, nn;

  /* Set all layers above the surface to surface layer               */
  if (!wincon->sigma) {
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = cs = wincon->s2[cc];
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	for (n = 0; n < windat->ntr; n++) {
	  if (wincon->trinfo_3d[n].advect ||
	      wincon->trinfo_3d[n].diffuse ||
	      wincon->trinfo_3d[n].diagn)
	    windat->tr_wc[n][c] = windat->tr_wc[n][cs];
	}
      }
    }
  }
}

/* END tr_set_surf()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to print a message on tracer min / max violations         */
/* Usage : check_bounds("text", window, windat, wincon);             */
/*-------------------------------------------------------------------*/
void check_bounds(char *text,                /* Warning text         */
		  geometry_t *window,        /* Window geometry      */
		  window_t *windat,          /* Window data          */
		  win_priv_t *wincon         /* Window constants     */
		  )
{
  int c, cc;                    /* Counters                          */
  int n, nn;                    /* Tracer number                     */
  double *min = wincon->mintr;
  double *max = wincon->maxtr;

  /*-----------------------------------------------------------------*/
  /* Print min/max violations before the tracer step                 */
  for (cc = 1; cc <= wincon->vc1; cc++) {
    c = wincon->s1[cc];
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      n = wincon->tbdy[nn];
      if (windat->tr_wc[n][c] < min[n] || 
	  windat->tr_wc[n][c] > max[n] || isnan(windat->tr_wc[n][c])) {
	hd_warn
	  ("%s : Tracer %s outside range : %e @ %d(%d %d %d) at %f\n",
	   text, wincon->trname[n], windat->tr_wc[n][c], c,
	   window->s2i[c], window->s2j[c], window->s2k[c], 
	   windat->t / 86400.0);
#if defined(HAVE_ECOLOGY_MODULE)
	if (windat->tr_wc[n][c] > 0.0 && windat->tr_wc[n][c] < TINY) {
	  hd_warn
	    ("%s : Tracer %s too small : %e @ (%d %d %d) at %f\n", text,
	     wincon->trname[n], windat->tr_wc[n][c],
	     window->s2i[c], window->s2j[c], window->s2k[c], windat->t / 86400.0);
	}
#endif
      }
    }
  }
}

/* END check_bounds()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to limit minimum values to TINY                           */
/*-------------------------------------------------------------------*/
void limit_min(geometry_t *window,            /* Window geometry     */
	       window_t *windat,              /* Window data         */
	       win_priv_t *wincon             /* Window constants    */
  )
{
  int c, cc;                    /* Counters                          */
  int n, nn;                    /* Tracer number                     */

  for (nn = 0; nn < wincon->ntbdy; nn++) {
    n = wincon->tbdy[nn];
    if (n != windat->tno && n != windat->sno) {
      for (cc = 1; cc <= wincon->vc1; cc++) {
	c = wincon->s1[cc];
	/*if (windat->tr_wc[n][c] != 0.0)*/
	  windat->tr_wc[n][c] = max(windat->tr_wc[n][c], TINY);
      }
    }
  }
}

/* END limit_min()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check the bounds on flux differences                   */
/* Usage: check_fluxes(n, "MPB_N", "c", dtu, window, windat, wincon) */
/*-------------------------------------------------------------------*/
void check_fluxes(int n,                         /* Tracer number    */
		  char *trname,                  /* Tracer name      */
		  char *text,                    /* Warning text     */
		  double dtu,                    /* Time step        */
		  geometry_t *window,            /* Window geometry  */
		  window_t *windat,              /* Window data      */
		  win_priv_t *wincon             /* Window constants */
  )
{
  int c, cs, cc;                    /* Counters                      */
  int c1, c2, c3, c4;
  double *Fx = wincon->w2;
  double *Fz = wincon->w4;
  double d1, f1, f2, f3;

  if (strcmp(trname, wincon->trname[n]) == 0) {
    double min = wincon->mintr[n];
    double max = wincon->maxtr[n];
    for (cc = 1; cc <= wincon->vc1; cc++) {
      c = wincon->s1[cc];
      cs = window->m2d[c];
      c1 = window->c2e[1][c];
      c2 = window->c2e[3][c];
      c3 = window->c2e[2][c];
      c4 = window->c2e[4][c];
      d1 = dtu / (window->cellarea[cs] * wincon->dz[c]);
      f1 = (Fx[c2] - Fx[c1]) * d1;
      f2 = (Fx[c3] - Fx[c4]) * d1;
      f3 = (Fz[window->zp1[c]] - Fz[c]) / wincon->dz[c];
      if (fabs(f1) > max)
	hd_warn
	  ("%s : Tracer %s x flux divergence %e outside range : (%e - %e) (%e - %e) @ (%d %d %d) : %d at %f. dz = %f\n",
	   text, wincon->trname[n], f1, windat->tr_wc[n][c],
	   windat->tr_wc[n][window->c2c[3][c]],Fx[c2],Fx[c1],
	   window->s2i[c], window->s2j[c], window->s2k[c], c, 
	   windat->t / 86400.0, wincon->dz[c]);
      if (fabs(f2) > max)
	hd_warn
	  ("%s : Tracer %s y flux divergence %e outside range : (%e - %e) (%e - %e) @ (%d %d %d) : %d at %f. dz = %f\n",
	   text, wincon->trname[n], f2, windat->tr_wc[n][c],
	   windat->tr_wc[n][window->c2c[2][c]],Fx[c3],Fx[c4],
	   window->s2i[c], window->s2j[c], window->s2k[c], c,
	   windat->t / 86400.0, wincon->dz[c]);
      if (fabs(f3) > max)
	hd_warn
	  ("%s : Tracer %s z flux divergence %e outside range : (%e - %e) (%e - %e) @ (%d %d %d) : %d at %f. dz = %f\n",
	   text, wincon->trname[n], f3, windat->tr_wc[n][c],
	   windat->tr_wc[n][window->zp1[c]],Fz[window->zp1[c]],Fz[c],
	   window->s2i[c], window->s2j[c], window->s2k[c], c,
	   windat->t / 86400.0, wincon->dz[c]);
    }
  }
}

/* END check_fluxes()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Flux corrected transport limiting function                        */
/*-------------------------------------------------------------------*/
void fct_limit(geometry_t *window, 
	       window_t *windat, 
	       win_priv_t *wincon, 
	       double *tr,
	       int ksf,
	       int kef,
	       int kefs,
	       double dtu,
	       double *osubeta,
	       double *subeta
)
{
  int c, c1, c2, cc, cn, cb, j, zp1, zm1;
  int ee, e;
  double *Fx = wincon->w2;      /* x tracer flux */
  double *Fz = wincon->w4;      /* z tracer flux */
  double *Fxh = wincon->Fxh;
  double *Fzh = wincon->Fzh;
  double *Ax = wincon->Ax;
  double *Az = wincon->Az;
  double *dtracer = wincon->w1;
  double *pp = wincon->crfxf;
  double *pm = wincon->crfxc;
  double *rp = wincon->crfyf;
  double *rm = wincon->crfyc;
  double *wmax = wincon->crfzf;
  double *wmin = wincon->crfzc;
  double *dz = wincon->d6;
  double qm, qp, lc;
  double top, surftr, waterbot;
  int slf = 0;

  /*-----------------------------------------------------------------*/
  /* Define the antidiffusive fluxes                                 */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    Ax[e] = Fxh[e] - Fx[e];
  }
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    Az[c] = Fzh[c] - Fz[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get bounds on the antidiffusive fluxes                          */

  /*-----------------------------------------------------------------*/
  /* Get the low order solution                                      */
  memcpy(wmax, tr, window->szc * sizeof(double));
  memcpy(wmin, tr, window->szc * sizeof(double));
  for (cc = ksf; cc <= kef; cc++) {
    c = wincon->s1[cc];   /* Wet cell to process                     */
    c2 = window->m2d[c];  /* 2D cell corresponding to 3D cell        */
    zp1 = window->zp1[c]; /* Cell above cell c                       */

    /* SIGMA : Adjust tracer values for the depth                    */
    if (slf)
      tr[c] *= wincon->Ds[c2];
    else
      tr[c] *= wincon->Hn1[c2];
    tr[c] -= ((dtracer[c] / window->cellarea[c2] +
	       (Fz[zp1] - Fz[c])) / wincon->dz[c]);
    tr[c] /= wincon->Hn1[c2];
    wmax[c] = max(wmax[c], tr[c]);
    wmin[c] = min(wmin[c], tr[c]);
  }

  for (cc = 1; cc <= kefs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    top = -Fz[c] * window->cellarea[c2];
    surftr =
      surf_conc(window, windat, wincon, tr, osubeta[c2], subeta[c2], c, c2,
		cc, top, dtracer);
    /* Set all layers above the surface to surface layer             */
    tr[c] = (window->botz[c2] > 0.0 && wincon->dz[c] < wincon->hmin) ?
      tr[c] : surftr;
    wmax[c] = max(wmax[c], tr[c]);
    wmin[c] = min(wmin[c], tr[c]);
  }

  /*-----------------------------------------------------------------*/
  /* Set lateral boundary conditions                                 */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c1 = window->bpt[cc];
    c2 = window->bin[cc];
    wmax[c1] = wmax[c2];
    wmin[c1] = wmin[c2];
  }
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    wmax[window->zm1[c]] = wmax[c];
    wmin[window->zm1[c]] = wmin[c];
  }

  memcpy(Fzh, wmax, window->szc * sizeof(double));
  memcpy(Fz, wmin, window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Limit the antidiffusive fluxes in the horizontal                */
  /* Get the sum of fluxes into and out of the cell                  */
  memset(pp, 0.0, window->szc * sizeof(double));
  memset(pm, 0.0, window->szc * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    for (j = 1; j <= window->npe[c2]; j++) {
      e = window->c2e[j][c];
      cn = window->c2c[j][c];
      pp[c] += -max(0.0, -window->eSc[j][c2] * Ax[e]);
      pm[c] += max(0.0, window->eSc[j][c2] * Ax[e]);
      wmax[c] = max(wmax[c], wmax[cn]);
      wmin[c] = min(wmin[c], wmax[cn]);
    }
    /*
    wmax[c] = max(max(wmax[c], wmax[zp1]), wmax[zm1]);
    wmin[c] = min(min(wmin[c], wmin[zp1]), wmin[zm1]);
    */
  }
  /* Mid-water layers                                                */
  for (cc = ksf; cc <= kef; cc++) {
    c = wincon->s1[cc];   /* Wet cell to process                     */
    c2 = window->m2d[c];
    qp = (wmax[c] - tr[c]) * window->cellarea[c2] * wincon->dz[c];
    qm = (tr[c] - wmin[c]) * window->cellarea[c2] * wincon->dz[c];
    rp[c] = (pp[c] > 0.0) ? min(1.0, qp / pp[c]) : 0.0;
    rm[c] = (pm[c] > 0.0) ? min(1.0, qm / pm[c]) : 0.0;
  }
  /* Surface layer                                                   */
  for (cc = 1; cc <= kefs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    cb = wincon->i1[cc];
    waterbot = (c == cb) ? window->botz[c2] : window->gridz[c];
    dz[c2] = subeta[c2] - waterbot;
    qp = (wmax[c] - tr[c]) * window->cellarea[c2] * dz[c2];
    qm = (tr[c] - wmin[c]) * window->cellarea[c2] * dz[c2];
    rp[c] = (pp[c] > 0.0) ? min(1.0, qp / pp[c]) : 0.0;
    rm[c] = (pm[c] > 0.0) ? min(1.0, qm / pm[c]) : 0.0;
  }
  /* Get the limiter and update the antidiffusive fluxes             */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    lc = (Ax[e] >= 0.0) ? min(rp[c1], rm[c2]) : min(rp[c2], rm[c1]);
    Ax[e] *= lc;
  }

  /*-----------------------------------------------------------------*/
  /* Limit the antidiffusive fluxes in the vertical                  */
  /* Get the sum of fluxes into and out of the cell                  */

  memcpy(wmax, Fzh, window->szc * sizeof(double));
  memcpy(wmin, Fz, window->szc * sizeof(double));

  memset(pp, 0.0, window->szc * sizeof(double));
  memset(pm, 0.0, window->szc * sizeof(double));
  for (cc = ksf; cc <= kef; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    pp[c] += max(0.0,Az[c]) - min(0.0, Az[zp1]);
    pm[c] += max(0.0,Az[zp1]) - min(0.0, Az[c]);

    wmax[c] = max(max(wmax[c], wmax[zp1]), wmax[zm1]);
    wmin[c] = min(min(wmin[c], wmin[zp1]), wmin[zm1]);

  }
  for (cc = 1; cc <= kefs; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    pp[c] += max(0.0,Az[c]);
    pm[c] -= min(0.0, Az[c]);

    wmax[c] = max(wmax[c], wmax[zm1]);
    wmin[c] = min(wmin[c], wmin[zm1]);

  }
  for (cc = ksf; cc <= kef; cc++) {
    c = wincon->s1[cc];
    qp = (wmax[c] - tr[c]) * wincon->dz[c];
    qm = (tr[c] - wmin[c]) * wincon->dz[c];
    rp[c] = (pp[c] > 0.0) ? min(1.0, qp / pp[c]) : 0.0;
    rm[c] = (pm[c] > 0.0) ? min(1.0, qm / pm[c]) : 0.0;
  }
  for (cc = 1; cc <= kefs; cc++) {
    c = wincon->s1[cc];
    c2 = window->m2d[c];
    qp = (wmax[c] - tr[c]) * dz[c2];
    qm = (tr[c] - wmin[c]) * dz[c2];
    rp[c] = (pp[c] > 0.0) ? min(1.0, qp / pp[c]) : 0.0;
    rm[c] = (pm[c] > 0.0) ? min(1.0, qm / pm[c]) : 0.0;
  }
  /* Get the limiter                                                 */
  for (cc = ksf; cc <= kef; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    lc = (Az[zp1] >= 0.0) ? min(rp[zp1], rm[c]) : min(rp[c], rm[zp1]);
    Az[c] *= lc;
  }
  for (cc = 1; cc <= kefs; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    Az[c] *= rm[c];
  }

  /*-----------------------------------------------------------------*/
  /* Get the horizontal flux divergence                              */
  memset(dtracer, 0, window->szc * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    dtracer[c] = 0.0;
    for (j = 1; j <= window->npe[c2]; j++) {
      e = window->c2e[j][c];
      dtracer[c] += (window->eSc[j][c2] * Ax[e]);
    }
  }
  memcpy(Fz, Az, window->szc * sizeof(double));
}

/* END fct_limit()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Invokes the ULTIMATE filter                                       */
/*-------------------------------------------------------------------*/
void ultimate_filter(double *cr, /* Courant number at the cell face  */
		     double *Fx, /* Scalar flux at the cell face     */
		     double *F,  /* Scalar value                     */
		     int e,      /* Face coordinate                  */
		     int c,      /* Cell coordinate                  */
		     double trp, /* Tracer value upstream of c       */
		     int cm1,    /* Cell downstream of c             */
		     double trm  /* Tracer value downstream of cm1   */
		     )
{
  double nc;                /* Normalised center value               */
  double nu;                /* Normalised face value                 */
  double tc = 0.0;          /* tracer value at the cell center       */
  double tu = 0.0;          /* Upstream tracer value                 */
  double td = 0.0;          /* Downstream tracer value               */
  double diff;              /* Difference of upstream and downstream */

  if (cr[e] == 0.0) {
      return;
  } else {
    /* Get the normalised tracer values */
    if (cr[e] > 0.0) {
      tc = F[cm1];
      tu = trm;
      td = F[c];
    } else {
      tc = F[c];
      tu = trp;
      td = F[cm1];
    }
    diff = td - tu;
    nc = nu = 0.0;
    if (diff) {
      nc = (tc - tu) / diff;
      nu = (Fx[e] - tu) / diff;
    }

    /* Get the adjusted face values */
    if (nc >= 0.0 && nc <= 1.0) {
      if (nu > 1.0)
	Fx[e] = td;
      else if (nu < nc)
	Fx[e] = tu + nc * (td - tu);
      else if (nu > nc / fabs(cr[e]))
	Fx[e] = tu + nc * (td - tu) / fabs(cr[e]);
    } else
      Fx[e] = tu + nc * (td - tu);
  }
}

/* END ultimate_filter()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Builds the weights to get the second derivative for advection.    */
/*-------------------------------------------------------------------*/
void build_advect_weights(geometry_t *window,  /* Processing window  */
			  window_t *windat,    /* Window data        */
			  win_priv_t *wincon   /* Window constants   */
			  )
{
  int ee, e;

  if (wincon->ac0 == NULL && wincon->ac1 == NULL) {
    wincon->ac0 = d_alloc_2d(2*window->npem, window->szeS);
    wincon->ac1 = d_alloc_2d(2*window->npem, window->szeS);
  }

  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    get_advect_metrics(window, windat, wincon, e);
  }
}

/* END build_advect_weights()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Gets the metrics of cells surrounding cell c using least squares  */
/* fitting, where the solution to npe equations of a polynomial of   */
/* degree nop is obtained via singular value decomposition.          */
/* The least squares fit uses cell centred values of cells           */
/* surrounding a given cell, c, and optionally (ef=1) the area       */
/* weighted value of surrounding vertices.                           */
/* The polynomial used is:                                           */
/* t = co + c1.x + c2.y + c3.x.x + c4.x.y + c5.y.y                   */
/* This can be formally differentialed to get the second derivative, */
/* which in this case is 2.c3.                                       */
/*-------------------------------------------------------------------*/
void get_advect_metrics(geometry_t *window,   /* Processing window   */
			window_t *windat,     /* Window data         */
			win_priv_t *wincon,   /* Window constants    */
			int e)                /* Edge index          */
{
  int n, cc, c, cn, i, j, jj;
  int en, vn;
  double theta, thb, thc, a, h, npe;
  /*int nop = 6;*/    /* Number of polynomial coefficients               */
  int nos = 3;    /* Coefficient number corresponding to x * x       */
  int ef = 1;     /* Include vertex values in least squares          */
  int dotest = 0; /* Perform testing on regression and weights       */
  int npem;       /* Number of points in the least squares fit       */
  double x, y, **p, **b;
  double *f, *s, *w, *std;
  double **coeff;

  for (n = 0; n <= 1; n++) {
    c = window->e2c[e][n];
    j = window->e2e[e][n];
    npe = npem = window->npe[c];
    if (ef) npem *= 2;

    coeff = (n == 0) ? wincon->ac0 : wincon->ac1;
    p = d_alloc_2d(nop, npem);
    b = d_alloc_2d(nop, npem);
    f = d_alloc_1d(nop);
    w = d_alloc_1d(nop);
    s = d_alloc_1d(npem);
    std = d_alloc_1d(npem);

    thb = atan2(window->u1y[e] - window->celly[c],
		window->u1x[e] - window->cellx[c]);
    jj = 0;

    /* Get the metrics at the surrounding cell centres. These are    */
    /* computed according to Fig 1 in Skamarock and Gassmann (2011). */
    for (j = 1; j <= npe; j++) {
      en = window->c2e[j][c];
      cn = window->c2c[j][c];
      thc = atan2(window->u1y[en] - window->celly[c],
		  window->u1x[en] - window->cellx[c]);
      theta = thc - thb;
      x = window->h2au1[en] * cos(theta);
      y = window->h2au1[en] * sin(theta);
      p[jj][0] = 1.0;
      p[jj][1] = x;
      p[jj][2] = y;
      p[jj][3] = x * x;
      p[jj][4] = x * y;
      p[jj][5] = y * y;
      s[jj] = 1.0;
      std[jj] = 0.0;
      jj++;
    }

    /* Get the metrics at the cell vertices if required              */
    if (ef) {
      for (j = 1; j <= npe; j++) {
	vn = window->c2v[j][c];
	en = window->v2e[vn][0];
	h = sqrt(0.25 * (window->h2au1[en] * window->h2au1[en] +
			 window->h1au1[en] * window->h1au1[en]));
	thc = atan2(window->gridy[vn] - window->celly[c],
		    window->gridx[vn] - window->cellx[c]);
	theta = thc - thb;
	x = h * cos(theta);
	y = h * sin(theta);
	p[jj][0] = 1.0;
	p[jj][1] = x;
	p[jj][2] = y;
	p[jj][3] = x * x;
	p[jj][4] = x * y;
	p[jj][5] = y * y;
	s[jj] = 1.0;
	std[jj] = 0.0;
	jj++;
      }
    }

    /* Least squares fitting via singular value decomposition, where */
    /* p = U * W * V^T, then f = B * s where b = V * W^-1 * U^T.     */
    /* These weights are computed using a scalar field (s) of 1.     */
    svd_lsq_B(p, nop, npem, s, NULL, b, f);

    /* Save the coefficients corresponding to x * x                  */
    for (j = 0; j < npem; j++) {
      coeff[e][j] = b[j][nos];
    }

    /* Perform testing to check things line up                       */
    if (dotest) {
      double trp;

      printf("edge %d dir=%d\n", e, n);
      for (jj = 0; jj < nop; jj++) {
	x = 0.0;
	for (j = 0; j < npem; j++) {
	  x += b[j][jj];
	}
	printf(" coeff%d = %e : %e, curv=%e\n",jj, x, f[jj], 2.0*f[nos]);
      }
      x = window->u1x[e];
      y = window->u1y[e];
      trp = f[0] + x * f[1] + y * f[2] + x * x * f[3] +
	x * y * f[4] + y * y * f[5];
      printf(" tracer value = %f(pred) : 1.0(act)\n",trp);
    }
    d_free_2d(p);
    d_free_2d(b);
    d_free_1d(f);
    d_free_1d(w);
    d_free_1d(s);
  }
}

/* END get_advect_metrics()                                          */
/*-------------------------------------------------------------------*/

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (geod_inv_geod_fwd_sodanos(DEG2RAD(x1), DEG2RAD(y1),\
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Builds the weights to get the second derivative for advection.    */
/*-------------------------------------------------------------------*/
void build_quadratic_weights(geometry_t *window, /* Processing window*/
			     window_t *windat,   /* Window data      */
			     win_priv_t *wincon  /* Window constants */
			     )
{
  int ee, e, es, c;
  int n, m, eoe, sb, msb, mn;
  int *mask = wincon->s1;
  int nuc = 3;    /* Number of polynomial coefficients to use        */
  int mrn = 2*window->npem; /* Maximum size of mask reset array      */
  int mask_r[mrn];          /* mask reset array                      */
  
  memset(mask_r, 0, mrn * sizeof(int));
  memset(mask, 0, window->szc * sizeof(int));
  if (wincon->Bcell == NULL && wincon->nBcell == NULL) {
    msb = 0;
    /* Count the cell centres associated with each edge              */
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      sb = 0;
      /* memset(mask, 0, window->szc * sizeof(int)); */
      /* Reset mask */
      for (mn=0; mn<mrn; mn++) {
	mask[mask_r[mn]] = 0;
	mask_r[mn] = 0;
      }
      mn = 0;
      for (n = 1; n <= window->nee[es]; n++) {
	eoe = window->eSe[n][e];
	if (!eoe) continue;
	for (m = 0; m <= 1; m++) {
	  c = window->e2c[eoe][m];
	  if (!mask[c]) {
	    mask[c] = 1;
	    mask_r[mn++] = c;
	    sb++;
	  }
	}
      }
      if (mn > mrn) hd_quit("build_quadratic_weights: error in mask reset array\n");
      msb = max(msb, sb);
    }

    /* Assign the cell centres                                       */
    wincon->Bcell = i_alloc_2d(msb, window->sze);
    wincon->nBcell = i_alloc_1d(window->sze);
    wincon->B = d_alloc_3d(nuc, msb, window->sze);
    memset(wincon->nBcell, 0, window->sze * sizeof(int));
    for (ee = 1; ee <= window->n3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      /* memset(mask, 0, window->szc * sizeof(int)); */
      /* Reset mask */
      for (mn=0; mn<mrn; mn++) {
	mask[mask_r[mn]] = 0;
	mask_r[mn] = 0;
      }
      mn = 0;
      for (n = 1; n <= window->nee[es]; n++) {
	eoe = window->eSe[n][e];
	if (!eoe) continue;
	for (m = 0; m <= 1; m++) {
	  c = window->e2c[eoe][m];
	  if (!mask[c]) {
	    mask[c] = 1;
	    mask_r[mn++] = c;
	    wincon->Bcell[e][wincon->nBcell[e]++] = c;
	  }
	}
      }
      get_quadratic_metrics(window, windat, wincon, e);
    }
  }
}

/* END build_quadratic_weights()                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Gets the metrics of cells surrounding edge e using least squares  */
/* fitting, where the solution to npe equations of a polynomial of   */
/* degree nop is obtained via singular value decomposition.          */
/* The least squares fit uses cell centred values of cells           */
/* surrounding a given edge, e.                                      */
/* The polynomial used is:                                           */
/* t = co + c1.x + c2.y + c3.x.x + c4.x.y + c5.y.y                   */
/* In this case we wish to interpolate along the x axis (y=0), so    */
/* the only coefficients required are co, c1 and c3.                 */
/*-------------------------------------------------------------------*/
void get_quadratic_metrics(geometry_t *window, /* Processing window  */
			   window_t *windat,   /* Window data        */
			   win_priv_t *wincon, /* Window constants   */
			   int e)              /* Edge index         */
{
  int n, m, cc, c, cs, cn, i, j, jj;
  int eoe, en, vn;
  double theta, thb, thc, a, h;
  /*int nop = 6;*/    /* Number of polynomial coefficients               */
  int nuc = 3;    /* Number of polynomial coefficients to use        */
  int uc[3] = { 0, 1, 3 };  /* Polynomial coefficients to use        */
  int sodanos = 0; /* Compute distances using sodanos' algorithm     */
  int dotest = 0; /* Perform testing on regression and weights       */
  double d, x, y, **p, **b;
  double *f, *s, *w, *std;
  double **coeff;
  int npem = wincon->nBcell[e];
  int *cells = wincon->Bcell[e];
  int es = window->m2de[e];
  double elat, elon, clat, clon;

  p = d_alloc_2d(nop, npem);
  b = d_alloc_2d(nop, npem);
  f = d_alloc_1d(nop);
  w = d_alloc_1d(nop);
  s = d_alloc_1d(npem);
  std = d_alloc_1d(npem);

  elon = window->u1x[es];
  elat = window->u1y[es];
  thb = window->thetau1[es];
  jj = 0;

  /* Get the metrics at the surrounding cell centres                 */
  for (n = 0; n < npem; n++) {
    c = cells[n];
    cs = window->m2d[c];
    clat = window->celly[cs];
    clon = window->cellx[cs];
    thc = atan2(clat - elat, clon - elon);
    theta = thc - thb;
    if (sodanos)
      d = GEODESIC(elon, elat, clon, clat);
    else {
      x = clon - elon;
      y = clat - elat;
      d = sqrt(x * x + y * y);
    }

    x = d * cos(theta);
    y = d * sin(theta);

    p[jj][0] = 1.0;
    p[jj][1] = x;
    p[jj][2] = y;
    p[jj][3] = x * x;
    p[jj][4] = x * y;
    p[jj][5] = y * y;
    s[jj] = 1.0;
    std[jj] = 0.0;
    jj++;
  }

  /* Least squares fitting via singular value decomposition, where   */
  /* p = U * W * V^T, then f = B *s where b = V * W^-1 * U^T.        */
  /* Note: these weights are computed using a scalar field (s) of 1. */
  svd_lsq_B(p, nop, npem, s, NULL, b, f);

  /* Save the coefficients corresponding to x * x                    */
  for (i = 0; i < nuc; i++) {
    jj = uc[i];
    for (j = 0; j < npem; j++) {
      wincon->B[e][j][i] = b[j][jj];
    }
  }

  /* Perform testing to check things line up                       */
  if (dotest) {
    double trp;

    printf("edge %d\n", e);
    for (jj = 0; jj < nop; jj++) {
      x = 0.0;
      for (j = 0; j < npem; j++) {
	x += b[j][jj];
      }
      printf(" coeff%d = %e : %e\n",jj, x, f[jj]);
    }

    j = 1;
    jj = (j == window->e2e[e][0]) ? 0 : 1;
    c = window->e2c[e][jj];
    /*x = 1.5 * window->hacell[j][window->m2d[c]];*/
    x = 3.0 * window->hacell[j][window->m2d[c]];
    y = x;
    trp = f[0] + x * f[1] + y * f[2] + x * x * f[3] +
      x * y * f[4] + y * y * f[5];
    printf(" tracer value @ %f = %f(pred) : 1.0(act)\n", x, trp);
  }
  d_free_2d(p);
  d_free_2d(b);
  d_free_1d(f);
  d_free_1d(w);
  d_free_1d(s);
}

/* END get_quadratic_metrics()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
double get_quadratic_value(geometry_t *window, /* Processing window  */
			   window_t *windat,   /* Window data        */
			   win_priv_t *wincon, /* Window constants   */
			   double *tr,         /* Tracer array       */
			   int e,              /* Edge index         */
			   int j               /* Edge direction     */
			   )
{
  int i, jj, n, c, ci;
  int nuc = 3;    /* Number of polynomial coefficients to use        */
  int npem;
  double coeff[nuc];
  double x, v;
  double tmx = -HUGE;
  double tmn = HUGE;
  double sc;

  npem  = wincon->nBcell[e];
  memset(coeff, 0, nuc * sizeof(double));
  /*sc = (j == 0) ? 1.5 : -1.5;*/
  sc = (j == 0) ? 3.0 : -3.0;

  /* Get the weights for the edge                                    */
  for (n = 0; n < npem; n++) {
    c = wincon->Bcell[e][n];
    tmx = max(tmx, tr[c]);
    tmn = min(tmn, tr[c]);
    for (jj = 0; jj < nuc; jj++) {
      coeff[jj] += tr[c] * wincon->B[e][n][jj];
    }
  }

  /* Get the value                                                   */
  c = window->e2c[e][j];
  if ((ci = window->wgst[c])) return(tr[ci]);

  jj = window->e2e[e][j];
  x = sc * window->hacell[jj][window->m2d[c]];
  v = coeff[0] + coeff[1] * x + coeff[2] * x * x;
  v = min(tmx, max(tmn, v));

  return(v);
}

/* END get_quadratic_value()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Builds the weights to get linear least squares interpolation.     */
/* This is an interpolation within a volumetric cell using a linear  */
/* squares interpolation with slope limiter to be conservative. The  */
/* polynomial function is:                                           */
/* v = co + c1.x + c2.y + c3.z                                       */
/* where (x,y,z) is the 3D position in the cell. Tracer values used  */
/* in the least squares fit are from one 'ring' around the cell in   */
/* the layer of the cell, the layer above and below the cell.        */
/* No-gradient functions are imposed at the top and bottom layers.   */
/*-------------------------------------------------------------------*/
void build_linear_weights(geometry_t *window, /* Processing window   */
			  window_t *windat,   /* Window data         */
			  win_priv_t *wincon  /* Window constants    */
			  )
{
  qweights *lw;
  int cc, c, cs, c2, ci, zp1, zm1, ee;
  int n, m, eoe, sb;
  int *mask = wincon->s1;
  int *st = NULL, sz, szc, szm;
  int i, cn, cns, cg, cgs;
  double dz, **xc, **yc, **zc;
  int nuc = 4;    /* Number of polynomial coefficients               */
  int size = 3;   /* One row surrounding the centre                  */
  int df = 2;     /* Dimension of interpolation (2D or 3D)           */

  if (df == 2) wincon->trasf |= TR_LIND2;
  if (df == 3) wincon->trasf |= TR_LIND3;
  if (df == 2) nuc = 3;

  if (wincon->lw == NULL) {

    /* Note: allocate qweights (model/lib/grid/include/lsqq.h)       */
    /* rather than lweights (model/lib/grid/include/lsql.h) as       */
    /* qweights contains 6 weights rather than 3 in lweights, and we */
    /* require 4 for this lsq function.                              */
    wincon->lw = malloc(window->szc * sizeof(qweights));

    /*---------------------------------------------------------------*/
    /* Get the number of points surrounding each centre              */
    szm = 0;
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      zp1 = window->zp1[c];
      zm1 = window->zm1[c];
      lw = &wincon->lw[c];

      /* Get the stencil at location c                               */
      sz = size;
      st = stencil(window, c, &sz, ST_SIZED, 0);
      szc = sz;
      lw->ncells = 0;
      /* Ghost cells to include                                      */
      add_cells_to_lsq(window, windat, wincon, sz, st, c);

      if (df == 3) {
	/* Get the stencil at location in the layer above c          */
	if (c != zp1) {
	  sz = size;
	  st = stencil(window, zp1, &sz, ST_SIZED, 0);
	}
	add_cells_to_lsq(window, windat, wincon, sz, st, c);

	/* Get the stencil at location in the layer below c          */
	sz = szc;
	if (zm1 != window->zm1[zm1]) {
	  sz = size;
	  st = stencil(window, zm1, &sz, ST_SIZED, 0);
	}
	add_cells_to_lsq(window, windat, wincon, sz, st, c);
      }
      szm = max(szm, lw->ncells);
    }

    /* OBC ghosts cells                                              */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	lw = &wincon->lw[c];
	zp1 = window->zp1[c];
	zm1 = window->zm1[c];
	for (i = 0; i < open->bgz; i++) {
	  sz = size;
	  st = stencil(window, c, &sz, ST_SIZED, 0);
	  lw->ncells = szc = sz;
	  if (df == 3) {
	    if (c != zp1) {
	      sz = size;
	      st = stencil(window, zp1, &sz, ST_SIZED, 0);
	    }
	    lw->ncells += sz;
	    sz = szc;
	    if (zm1 != window->zm1[zm1]) {
	      sz = size;
	      st = stencil(window, zm1, &sz, ST_SIZED, 0);
	    }
	    lw->ncells += sz;
	  }
	  szm = max(szm, lw->ncells);
	  c = open->omape[ee][c];
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate                                                      */
    xc = d_alloc_2d(szm, window->szc);
    yc = d_alloc_2d(szm, window->szc);
    zc = d_alloc_2d(szm, window->szc);

    for(cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      window->cellz[c] = 0.5 * window->gridz[c];
    }

    /*---------------------------------------------------------------*/
    /* Get the points surrounding each centre                        */
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      zp1 = window->zp1[c];
      zm1 = window->zm1[c];
      lw = &wincon->lw[c];
      lw->cells = i_alloc_1d(szm);
      lw->B = d_alloc_2d(nuc, szm);

      /* Get the stencil at location c                               */
      m = 0;
      sz = size;
      st = stencil(window, c, &sz, ST_SIZED, 0);
      szc = sz;

      /* Ghost cells to include                                      */
      for (i = 0; i < sz; i++) {
	lw->cells[m] = cn = st[i];
	cns = window->m2d[cn];
	xc[c][m] = window->cellx[cns];
	yc[c][m] = window->celly[cns];
	zc[c][m] = window->cellz[cn];
	m++;
	for(n = 1; n <= window->npe[cns]; n++) {
	  cg = window->c2c[n][cn];
	  cgs = window->m2d[cg];
	  if (window->wgst[cg]) {
	    lw->cells[m] = cg;
	    xc[c][m] = window->cellx[cgs];
	    yc[c][m] = window->celly[cgs];
	    zc[c][m] = window->cellz[cn];
	    m++;
	  }
	}
      }

      if (df == 3) {
	/* Get the stencil at location in the layer above c          */
	if (c != zp1) {
	  sz = size;
	  st = stencil(window, zp1, &sz, ST_SIZED, 0);
	}
	/* Ghost cells to include                                    */
	for (i = 0; i < sz; i++) {
	  lw->cells[m] = cn = st[i];
	  cns = window->m2d[cn];
	  xc[c][m] = window->cellx[cns];
	  yc[c][m] = window->celly[cns];
	  dz = 2.0 * (window->cellz[cn] - window->gridz[cn]);
	  zc[c][m] = (c == zp1) ? window->cellz[cn] + dz : window->cellz[cn];
	  m++;
	  for(n = 1; n <= window->npe[cns]; n++) {
	    cg = window->c2c[n][cn];
	    cgs = window->m2d[cg];
	    if (window->wgst[cg]) {
	      lw->cells[m] = cg;
	      xc[c][m] = window->cellx[cgs];
	      yc[c][m] = window->celly[cgs];
	      zc[c][m] = (c == zp1) ? window->cellz[cn] + dz : window->cellz[cn];
	      m++;
	    }
	  }
	}
	
	/* Get the stencil at location in the layer below c          */
	sz = szc;
	if (zm1 != window->zm1[zm1]) {
	  sz = size;
	  st = stencil(window, zm1, &sz, ST_SIZED, 0);
	} else {
	  sz = size;
	  st = stencil(window, c, &sz, ST_SIZED, 0);
	}
	/* Ghost cells to include                                      */
	for (i = 0; i < sz; i++) {
	  lw->cells[m] = cn = st[i];
	  cns = window->m2d[cn];
	  xc[c][m] = window->cellx[cns];
	  yc[c][m] = window->celly[cns];
	  dz = 2.0 * (window->cellz[cn] - max(window->gridz[cn], window->botz[cns]));
	  zc[c][m] = (zm1 == window->zm1[zm1]) ? window->cellz[cn] - dz : window->cellz[cn];
	  m++;
	  for(n = 1; n <= window->npe[cns]; n++) {
	    cg = window->c2c[n][cn];
	    cgs = window->m2d[cg];
	    if (window->wgst[cg]) {
	      lw->cells[m] = cg;
	      xc[c][m] = window->cellx[cgs];
	      yc[c][m] = window->celly[cgs];
	      zc[c][m] = (zm1 == window->zm1[zm1]) ? window->cellz[cn] - dz : window->cellz[cn];
	      m++;
	    }
	  }
	}
      }
    }

    /* OBC ghosts                                                    */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	ci = open->obc_e2[ee];
	zp1 = window->zp1[c];
	zm1 = window->zm1[c];
	lw = &wincon->lw[c];
	lw->cells = i_alloc_1d(szm);
	lw->B = d_alloc_2d(nuc, szm);
	for (i = 0; i < open->bgz; i++) {
	  m = 0;
	  sz = size;
	  st = stencil(window, c, &sz, ST_SIZED, 0);
	  szc = sz;
	  for (i = 0; i < sz; i++) {
	    lw->cells[m] = cn = st[i];
	    cns = window->m2d[cn];
	    xc[c][m] = window->cellx[cns];
	    yc[c][m] = window->celly[cns];
	    zc[c][m] = window->cellz[ci];
	    m++;
	  }
	  if (df == 3) {
	    if (c != zp1) {
	      sz = size;
	      st = stencil(window, zp1, &sz, ST_SIZED, 0);
	    }
	    for (i = 0; i < sz; i++) {
	      lw->cells[m] = cn = st[i];
	      cns = window->m2d[cn];
	      xc[c][m] = window->cellx[cns];
	      yc[c][m] = window->celly[cns];
	      dz = 2.0 * (window->cellz[ci] - window->gridz[ci]);
	      zc[c][m] = (c == zp1) ? window->cellz[ci] + dz : window->cellz[ci];
	      m++;
	    }
	    sz = szc;
	    if (zm1 != window->zm1[zm1]) {
	      sz = size;
	      st = stencil(window, zm1, &sz, ST_SIZED, 0);
	    }
	    for (i = 0; i < sz; i++) {
	      lw->cells[m] = cn = st[i];
	      cns = window->m2d[cn];
	      xc[c][m] = window->cellx[cns];
	      yc[c][m] = window->celly[cns];
	      dz = 2.0 * (window->cellz[ci] - max(window->gridz[ci], window->botz[window->m2d[ci]]));
	      zc[c][m] = (zm1 == window->zm1[zm1]) ? window->cellz[ci] - dz : window->cellz[ci];
	      m++;
	    }
	  }
	  c = open->omape[ee][c];
	}
      }
    }

    /* Get the SVD and weights                                       */
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      get_linear_metrics(window, windat, wincon, c, xc, yc, zc);
    }
  }
}

/* END build_linear_weights()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Increments the number of cells to use in the least squares fit.   */
/*-------------------------------------------------------------------*/
void add_cells_to_lsq(geometry_t *window,     /* Processing window   */
		      window_t *windat,       /* Window data         */
		      win_priv_t *wincon,     /* Window constants    */
		      int ncells,             /* # cells in stencil  */
		      int *cells,             /* Stencil cells       */
		      int co                  /* Centre cell         */
		      )
{
  int cc, c, cs, j;
  int cn, cns;
  qweights *lw = &wincon->lw[co];

  /* Wet cells to include                                            */
  lw->ncells += ncells;

  /* Ghost cells to include                                          */
  for (cc = 0; cc < ncells; cc++) {
    c = cells[cc];
    cs = window->m2d[c];
    for(j = 1; j <= window->npe[cs]; j++) {
      cn = window->c2c[j][c];
      cns = window->m2d[cn];
      if (window->wgst[cn])
	lw->ncells++;
    }
  }
}

/* END add_cells_to_lsq()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Gets the metrics of cells surrounding centre c using least        */
/* squares fitting, where the solution to npe equations of a         */
/* polynomial of degree nuc is obtained via singular value           */
/* decomposition.                                                    */
/* The least squares fit uses cell centred values of cells           */
/* surrounding a given location, c.                                  */
/* The polynomial used is:                                           */
/* v = co + c1.x + c2.y + c3.z                                       */
/*-------------------------------------------------------------------*/
void get_linear_metrics(geometry_t *window,    /* Processing window  */
			window_t *windat,      /* Window data        */
			win_priv_t *wincon,    /* Window constants   */
			int co,                /* Centre index       */
			double **xc,
			double **yc,
			double **zc
			)
{
  qweights *lw = &wincon->lw[co];
  int n, m, cc, c, cn, i, j, jj;
  int nuc = 4;    /* Number of polynomial coefficients               */
  double x, y, z, **p, **b;
  double *f, *s, *w, *std;
  double **coeff;
  int npem = lw->ncells;
  int *cells = lw->cells;
  int cs = window->m2d[co];
  int df = (wincon->trasf & TR_LIND2) ? 2 : 3;

  if (df == 2) nuc = 3;

  p = d_alloc_2d(nuc, npem);
  b = d_alloc_2d(nuc, npem);
  f = d_alloc_1d(nuc);
  w = d_alloc_1d(nuc);
  s = d_alloc_1d(npem);
  std = d_alloc_1d(npem);

  /* Get the metrics at the surrounding cell centres                 */
  jj = 0;
  for (n = 0; n < npem; n++) {
    c = cells[n];
    x = window->cellx[cs] - xc[co][n];
    y = window->celly[cs] - yc[co][n];
    z = window->cellz[co] - zc[co][n];

    p[jj][0] = 1.0;
    p[jj][1] = x;
    p[jj][2] = y;
    if (df == 3) p[jj][3] = z;
    s[jj] = 1.0;
    std[jj] = 0.0;
    jj++;
  }

  /* Least squares fitting via singular value decomposition, where   */
  /* p = U * W * V^T, then f = B *s where b = V * W^-1 * U^T.        */
  /* Note: these weights are computed using a scalar field (s) of 1. */
  svd_lsq_B(p, nuc, npem, s, NULL, b, f);

  /* Save the coefficients                                           */
  for (jj = 0; jj < nuc; jj++) {
    for (j = 0; j < npem; j++) {
      lw->B[j][jj] = b[j][jj];
    }
  }

  /* Set the weights. These are over-written later for each tracer.  */
  memset(lw->w, 0, nuc * sizeof(double));
  for (jj = 0; jj < nuc; jj++) {
    for (j = 0; j < npem; j++) {
      lw->w[jj] += s[jj] * lw->B[j][jj];
    }
  }

  d_free_2d(p);
  d_free_2d(b);
  d_free_1d(f);
  d_free_1d(w);
  d_free_1d(s);
}

/* END get_linear_metrics()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns an interpolation from a least squares linear function.    */
/*-------------------------------------------------------------------*/
double get_linear_value(geometry_t *window, /* Processing window     */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon, /* Window constants      */
			double *tr,         /* Tracer array          */
			int co,             /* Centre index          */
			double xi,          /* x location to         */
			double yi,          /* y location to         */
			double zi           /* z location to         */
			)
{
  qweights *lw = &wincon->lw[co];
  int cs = window->m2d[co];
  double val, x, y, z;
  int df = (wincon->trasf & TR_LIND2) ? 2 : 3;

  /* Get the value                                                   */
  x = window->cellx[cs] - xi;
  y = window->celly[cs] - yi;
  z = window->cellz[co] - zi;

  if (df == 3)
    val = lw->w[0] + lw->beta * (lw->w[1] * x + lw->w[2] * y + lw->w[3] * z);
  else
    val = lw->w[0] + lw->beta * (lw->w[1] * x + lw->w[2] * y);
  /*val = lw->w[0];*/
  return(val);
}

/* END get_linear_value()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Applies a slope limiter to a linear least squares polynomial to   */
/* ensure monotinicity.                                              */
/*-------------------------------------------------------------------*/
void get_linear_limit(geometry_t *window,   /* Processing window     */
		      window_t *windat,     /* Window data           */
		      win_priv_t *wincon,   /* Window constants      */
		      double *tr            /* Tracer array          */
		      )
{
  qweights *lw;
  int nuc = 4;    /* Number of polynomial coefficients to use        */
  int cc, c, cn, cs, c1, zp1, zm1;
  int j, jj, v, vp1, vm1, vs, n;
  double vmin, vmax, cmax, cmin, beta;
  double x, y, z, d, val;
  int df = (wincon->trasf & TR_LIND2) ? 2 : 3;

  if (df == 2) nuc = 3;

  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];

    lw = &wincon->lw[c];
    lw->beta = beta = 1.0;
    lw->tmx = -HUGE;
    lw->tmn = HUGE;

    /*---------------------------------------------------------------*/
    /* Get the weights for the linear function                       */
    memset(lw->w, 0, nuc * sizeof(double));
    /* Get the weights for the edge                                  */
    for (j = 0; j < lw->ncells; j++) {
      cn = lw->cells[j];
      for (jj = 0; jj < nuc; jj++) {
	lw->w[jj] += (tr[cn] * lw->B[j][jj]);
      }
      lw->tmx = max(lw->tmx, tr[cn]);
      lw->tmn = min(lw->tmn, tr[cn]);
    }

    /*---------------------------------------------------------------*/
    /* Reset the zeroth order weight to the tracer value             */
    lw->w[0] = tr[c];

    /*---------------------------------------------------------------*/
    /* Limit the weights                                             */
    for (j = 1; j <= window->npe[cs]; j++) {
      v = window->c2v[j][c];
      vp1 = window->c2v[j][zp1];
      vm1 = window->c2v[j][zm1];
      vmin = HUGE;
      vmax = -HUGE;

      /*-------------------------------------------------------------*/
      /* Get the min/max values in layer c                           */
      if (v) {
	vs = window->m2dv[v];
	x = window->gridx[vs];
	y = window->gridy[vs];
	z = window->cellz[c];
	val = get_linear_value(window, windat, wincon, tr, c, x, y, z);

	d = val - lw->w[0];
	/* Get the min/max values in the current layer c             */
	n = 0;
	for (jj = 1; jj <= window->nvc[vs]; jj++) {
	  c1 = window->v2c[v][jj];
	  if (!c1) continue;
	  vmin = min(vmin, tr[c1]);
	  vmax = max(vmax, tr[c1]);
	  n++;
	}
	cmax = vmax;
	cmin = vmin;
	/*if (n <= 1) continue;*/

	/* Find the slope limiting factor that keeps the             */
	/* interpolated value monotonic.                             */
	if (d != 0.0 && val > vmax) {
	  beta = max(min((vmax - lw->w[0]) / d, beta), 0.0);
	}	
	if (d != 0.0 && val < vmin) {
	  beta = max(min((vmin - lw->w[0]) / d, beta), 0.0);
	}

	/*-----------------------------------------------------------*/
	/* Get interpolated values at vertices of the upper layer of */
	/* the volumetric cell and limit.                            */
	if (df == 3) {
	  z = window->cellz[c] + (window->cellz[c] - window->gridz[c]);
	  val = get_linear_value(window, windat, wincon, tr, c, x, y, z);
	  d = val - lw->w[0];
	  /* Get the min/max values in the layer above c if not the  */
	  /* surface layer.                                          */
	  if (vp1 && c != zp1) {
	    for (jj = 1; jj <= window->nvc[vs]; jj++) {
	      c1 = window->v2c[vp1][jj];
	      vmin = min(vmin, tr[c1]);
	      vmax = max(vmax, tr[c1]);
	    }
	  }

	  /* Find the slope limiting factor that keeps the           */
	  /* interpolated value monotonic.                           */
	  if (d != 0.0 && val > vmax) {
	    beta = max(min((vmax - lw->w[0]) / d, beta), 0.0);
	  }	
	  if (d != 0.0 && val < vmin) {
	    beta = max(min((vmin - lw->w[0]) / d, beta), 0.0);
	  }

	  /*---------------------------------------------------------*/
	  /* Get interpolated values at vertices of the lower layer  */
	  /* of the volumetric cell and limit.                       */
	  z = window->cellz[c] - (window->cellz[c] - max(window->gridz[c], window->botz[cs]));
	  val = get_linear_value(window, windat, wincon, tr, c, x, y, z);
	  d = val - lw->w[0];
	  /* Get the min/max values in the layer below c if not the  */
	  /* bottom layer.                                           */
	  vmax = cmax; vmin = cmin;
	  if (vm1 && zm1 != window->zm1[zm1]) {
	    for (jj = 1; jj <= window->nvc[vs]; jj++) {
	      c1 = window->v2c[vm1][jj];
	      vmin = min(vmin, tr[c1]);
	      vmax = max(vmax, tr[c1]);
	    }
	  }
	  /* Find the slope limiting factor that keeps the           */
	  /* interpolated value monotonic.                           */
	  if (d != 0.0 && val > vmax) {
	    beta = max(min((vmax - lw->w[0]) / d, beta), 0.0);
	  }
	  if (d != 0.0 && val < vmin) {
	    beta = max(min((vmin - lw->w[0]) / d, beta), 0.0);
	  }
	}
      }
    }
    lw->beta = min(max(beta, 0.0), 1.0);
    /*lw->beta = 0.0;*/
  }

  /* Streamlines should not make their way into ghost cells, but if  */
  /* they do set the leading term of the weight to the ghost value,  */
  /* and beta = 0.                                                   */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    lw = &wincon->lw[c];
    lw->w[0] = lw->tmx = lw->tmn = tr[c];
    lw->beta = 0.0;
    for (jj = 1; jj < nuc; jj++)
      lw->w[jj] = 0.0;
  }
}

/* END get_linear_limit()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Builds the weights to get 2nd order least squares interpolation.  */
/* This is an interpolation within a volumetric cell using a linear  */
/* squares interpolation with slope limiter to be conservative. The  */
/* polynomial function is:                                           *
/* v = co + c1.x + c2.y + c3.x.x + c4.x.y + c5.y.y                   */
/* where (x,y) is the 2D position in the cell. Tracer values used    */
/* in the least squares fit are from one 'ring' around the cell in   */
/* the layer of the cell.                                            */
/*-------------------------------------------------------------------*/
void build_order2_weights(geometry_t *window, /* Processing window   */
			  window_t *windat,   /* Window data         */
			  win_priv_t *wincon  /* Window constants    */
			  )
{
  qweights *lw;
  int cc, c, cs, c2, ci, ee;
  int n, m, eoe, sb;
  int *mask = wincon->s1;
  int *st = NULL, sz, szc, szm;
  int i, cn, cns, cg, cgs;
  double dz, **xc, **yc;
  int nuc = 6;    /* Number of polynomial coefficients               */
  int size = 5;   /* One row surrounding the centre                  */

  if (wincon->lw == NULL) {

    /* Note: allocate qweights (model/lib/include/lsqq.h) rather     */
    /* than lweights (model/lib/include/lsql.h) as qweights contains */
    /* 6 weights rather than 3 in lweights, and we require 6 for     */
    /* this lsq function.                                            */
    wincon->lw = malloc(window->szc * sizeof(qweights));

    /*---------------------------------------------------------------*/
    /* Get the number of points surrounding each centre              */
    szm = 0;
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      lw = &wincon->lw[c];

      /* Get the stencil at location c                               */
      sz = size;
      st = stencil(window, c, &sz, ST_SIZED, 0);
      szc = sz;
      lw->ncells = 0;
      /* Ghost cells to include                                      */
      add_cells_to_lsq(window, windat, wincon, sz, st, c);
      szm = max(szm, lw->ncells);
    }

    /* OBC ghosts cells                                              */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	lw = &wincon->lw[c];
	for (i = 0; i < open->bgz; i++) {
	  sz = size;
	  st = stencil(window, c, &sz, ST_SIZED, 0);
	  lw->ncells = szc = sz;
	  szm = max(szm, lw->ncells);
	  c = open->omape[ee][c];
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Allocate                                                      */
    xc = d_alloc_2d(szm, window->szc);
    yc = d_alloc_2d(szm, window->szc);

    /*---------------------------------------------------------------*/
    /* Get the points surrounding each centre                        */
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      lw = &wincon->lw[c];
      lw->cells = i_alloc_1d(szm);
      lw->B = d_alloc_2d(nuc, szm);

      /* Get the stencil at location c                               */
      m = 0;
      sz = size;
      st = stencil(window, c, &sz, ST_SIZED, 0);
      szc = sz;

      /* Ghost cells to include                                      */
      for (i = 0; i < sz; i++) {
	lw->cells[m] = cn = st[i];
	cns = window->m2d[cn];
	xc[c][m] = window->cellx[cns];
	yc[c][m] = window->celly[cns];
	m++;
	for(n = 1; n <= window->npe[cns]; n++) {
	  cg = window->c2c[n][cn];
	  cgs = window->m2d[cg];
	  if (window->wgst[cg]) {
	    lw->cells[m] = cg;
	    xc[c][m] = window->cellx[cgs];
	    yc[c][m] = window->celly[cgs];
	    m++;
	  }
	}
      }
    }

    /* OBC ghosts                                                    */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      int ee;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	ci = open->obc_e2[ee];
	lw = &wincon->lw[c];
	lw->cells = i_alloc_1d(szm);
	lw->B = d_alloc_2d(nuc, szm);
	for (i = 0; i < open->bgz; i++) {
	  m = 0;
	  sz = size;
	  st = stencil(window, c, &sz, ST_SIZED, 0);
	  szc = sz;
	  for (i = 0; i < sz; i++) {
	    lw->cells[m] = cn = st[i];
	    cns = window->m2d[cn];
	    xc[c][m] = window->cellx[cns];
	    yc[c][m] = window->celly[cns];
	    m++;
	  }
	  c = open->omape[ee][c];
	}
      }
    }

    /* Get the SVD and weights                                       */
    for(cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      get_order2_metrics(window, windat, wincon, c, xc, yc);
    }
  }
}

/* END build_order2_weights()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Gets the metrics of cells surrounding centre c using least        */
/* squares fitting, where the solution to npe equations of a         */
/* polynomial of degree nuc is obtained via singular value           */
/* decomposition.                                                    */
/* The least squares fit uses cell centred values of cells           */
/* surrounding a given location, c.                                  */
/* The polynomial used is:                                           */
/* v = co + c1.x + c2.y + c3.x.x + c4.x.y + c5.y.y                   */
/*-------------------------------------------------------------------*/
void get_order2_metrics(geometry_t *window,    /* Processing window  */
			window_t *windat,      /* Window data        */
			win_priv_t *wincon,    /* Window constants   */
			int co,                /* Centre index       */
			double **xc,
			double **yc
			)
{
  qweights *lw = &wincon->lw[co];
  int n, m, cc, c, cn, i, j, jj;
  int nuc = 6;    /* Number of polynomial coefficients               */
  double x, y, **p, **b;
  double *f, *s, *w, *std;
  double **coeff;
  int npem = lw->ncells;
  int *cells = lw->cells;
  int cs = window->m2d[co];

  p = d_alloc_2d(nuc, npem);
  b = d_alloc_2d(nuc, npem);
  f = d_alloc_1d(nuc);
  w = d_alloc_1d(nuc);
  s = d_alloc_1d(npem);
  std = d_alloc_1d(npem);

  /* Get the metrics at the surrounding cell centres                 */
  jj = 0;
  for (n = 0; n < npem; n++) {
    c = cells[n];
    x = window->cellx[cs] - xc[co][n];
    y = window->celly[cs] - yc[co][n];

    p[jj][0] = 1.0;
    p[jj][1] = x;
    p[jj][2] = y;
    p[jj][3] = x * x;
    p[jj][4] = x * y;
    p[jj][5] = y * y;
    s[jj] = 1.0;
    std[jj] = 0.0;
    jj++;
  }

  /* Least squares fitting via singular value decomposition, where   */
  /* p = U * W * V^T, then f = B *s where b = V * W^-1 * U^T.        */
  /* Note: these weights are computed using a scalar field (s) of 1. */
  svd_lsq_B(p, nuc, npem, s, NULL, b, f);

  /* Save the coefficients                                           */
  for (jj = 0; jj < nuc; jj++) {
    for (j = 0; j < npem; j++) {
      lw->B[j][jj] = b[j][jj];
    }
  }

  /* Set the weights. These are over-written later for each tracer.  */
  memset(lw->w, 0, nuc * sizeof(double));
  for (jj = 0; jj < nuc; jj++) {
    for (j = 0; j < npem; j++) {
      lw->w[jj] += s[jj] * lw->B[j][jj];
    }
  }
  d_free_2d(p);
  d_free_2d(b);
  d_free_1d(f);
  d_free_1d(w);
  d_free_1d(s);
}

/* END get_order2_metrics()                                          */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Returns an interpolation from a least squares second order        */
/* function.                                                         */
/*-------------------------------------------------------------------*/
double get_order2_value(geometry_t *window, /* Processing window     */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon, /* Window constants      */
			double *tr,         /* Tracer array          */
			int co,             /* Centre index          */
			double xi,          /* x location to         */
			double yi,          /* y location to         */
			double z            /* For compatibility     */
			)
{
  qweights *lw = &wincon->lw[co];
  int cs = window->m2d[co];
  double val, x, y;

  /* Get the value                                                   */
  x = window->cellx[cs] - xi;
  y = window->celly[cs] - yi;

  val = lw->w[0] + lw->beta * (lw->w[1] * x + lw->w[2] * y) +
                   lw->alpha * (lw->w[3] * x * x +
				lw->w[4] * x * y +
				lw->w[5] * y * y);
  return(val);
}

/* END get_order2_value()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Applies a slope limiter to a second order least squares           */
/* polynomial to ensure monotinicity.                                */
/*-------------------------------------------------------------------*/
void get_order2_limit(geometry_t *window,   /* Processing window     */
		      window_t *windat,     /* Window data           */
		      win_priv_t *wincon,   /* Window constants      */
		      double *tr            /* Tracer array          */
		      )
{
  qweights *lw;
  int nuc = 6;    /* Number of polynomial coefficients to use        */
  int cc, c, cn, cs, c1;
  int j, jj, v, vs, n;
  double vmin, vmax, cmax, cmin, beta, alpha;
  double x, y, d, val;

  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];

    lw = &wincon->lw[c];
    lw->beta = lw->alpha = beta = alpha = 1.0;
    lw->tmx = -HUGE;
    lw->tmn = HUGE;

    /*---------------------------------------------------------------*/
    /* Get the weights for the linear function                       */
    memset(lw->w, 0, nuc * sizeof(double));
    /* Get the weights for the edge                                  */
    for (j = 0; j < lw->ncells; j++) {
      cn = lw->cells[j];
      for (jj = 0; jj < nuc; jj++) {
	lw->w[jj] += (tr[cn] * lw->B[j][jj]);
      }
      lw->tmx = max(lw->tmx, tr[cn]);
      lw->tmn = min(lw->tmn, tr[cn]);
    }

    /*---------------------------------------------------------------*/
    /* Reset the zeroth order weight to the tracer value             */
    lw->w[0] = tr[c];

    /*---------------------------------------------------------------*/
    /* Limit the weights                                             */
    for (j = 1; j <= window->npe[cs]; j++) {
      v = window->c2v[j][c];
      vmin = HUGE;
      vmax = -HUGE;

      /*-------------------------------------------------------------*/
      /* Get the min/max values in layer c                           */
      if (v) {
	vs = window->m2dv[v];
	x = window->gridx[vs];
	y = window->gridy[vs];
	val = get_order2_value(window, windat, wincon, tr, c, x, y, 0.);
	d = val - lw->w[0];

	/* Get the min/max values in the current layer c             */
	n = 0;
	for (jj = 1; jj <= window->nvc[vs]; jj++) {
	  c1 = window->v2c[v][jj];
	  if (!c1) continue;
	  vmin = min(vmin, tr[c1]);
	  vmax = max(vmax, tr[c1]);
	  n++;
	}
	cmax = vmax;
	cmin = vmin;
	/*if (n <= 1) continue;*/

	/* Find the slope limiting factor that keeps the             */
	/* interpolated value monotonic.                             */
	if (d != 0.0 && val > vmax) {
	  lw->alpha = 0.0;
	  val = get_order2_value(window, windat, wincon, tr, c, x, y, 0.);
	  d = val - lw->w[0];
	  if (d != 0.0) {
	    beta = max(min((vmax - lw->w[0]) / d, beta), 0.0);
	    lw->alpha = 1.0;
	    val = get_order2_value(window, windat, wincon, tr, c, x, y, 0.);
	    d = val - lw->w[0];
	    if (d != 0.0) alpha = max(min((vmax - lw->w[0]) / d, alpha), 0.0);
	  }
	}	
	if (d != 0.0 && val < vmin) {
	  lw->alpha = 0.0;
	  val = get_order2_value(window, windat, wincon, tr, c, x, y, 0.);
	  d = val - lw->w[0];
	  if (d != 0.0) {
	    beta = max(min((vmin - lw->w[0]) / d, beta), 0.0);
	    lw->alpha = 1.0;
	    val = get_order2_value(window, windat, wincon, tr, c, x, y, 0.);
	    d = val - lw->w[0];
	    if (d != 0.0) alpha = max(min((vmin - lw->w[0]) / d, alpha), 0.0);
	  }
	}
      }
    }
    lw->beta = min(max(beta, 0.0), 1.0);
    lw->alpha = min(max(alpha, 0.0), 1.0);
    /*lw->beta = 1.0;*/
  }

  /* Streamlines should not make their way into ghost cells, but if  */
  /* they do set the leading term of the weight to the ghost value,  */
  /* and beta = 0.                                                   */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    lw = &wincon->lw[c];
    lw->w[0] = lw->tmx = lw->tmn = tr[c];
    lw->beta = 0.0;
    for (jj = 1; jj < nuc; jj++)
      lw->w[jj] = 0.0;
  }
}

/* END get_order2_limit()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Gets the second derivative for high order advection               */
/*-------------------------------------------------------------------*/
double get_deriv2(geometry_t *window,   /* Processing window         */
		  window_t *windat,     /* Window data               */
		  win_priv_t *wincon,   /* Window constants          */
		  double *tr,           /* Tracer array              */
		  int e,                /* Edge index                */
		  int j                 /* Edge direction (0 or 1)   */
		  )
{
  int jj, n, cn, vn;
  int es = window->m2de[e];
  int c = window->e2c[e][j];
  int cs = window->m2d[c];
  int npe = window->npe[cs];
  double *coef;
  double trv, ret = 0.0;
  int ef = 1;     /* Include vertex values in least squares          */

  coef = (j == 0) ? wincon->ac0[es] : wincon->ac1[es];
  jj = 0;
  for (n = 1; n <= npe; n++) {
    cn = window->c2c[n][c];
    ret += coef[jj++] * tr[cn];
  }
  if (ef) {
    for (n = 1; n <= npe; n++) {
      vn = window->c2v[n][c];
      trv = vertex_weighted_tr(window, windat, wincon, tr, vn);
      ret += coef[jj++] * trv;
    }
  }

  ret *= 2.0;

  return(ret);
}

/* END get_deriv2()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Computes the area weighted tracer vertex value                    */
/*-------------------------------------------------------------------*/
double vertex_weighted_tr(geometry_t *window,   /* Processing window */
			  window_t *windat,     /* Window data       */
			  win_priv_t *wincon,   /* Window constants  */
			  double *tr,           /* Tracer array      */
			  int v                 /* Vertex index      */
			  )
{
  int cc, c;
  int vs = window->m2dv[v];
  double d1 = 1.0 / window->dualarea[vs];
  double ret = 0.0;

  for (cc = 1; cc <= window->nvc[vs]; cc++) {
    c = window->v2c[v][cc];
    if (c) {
      ret += tr[c] * window->dualareap[vs][cc];

    }
  }
  ret /= window->dualarea[vs];
  return(ret);
}

/* END vertex_weighted_tr()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_t(geometry_t *window,      /* Window geometry       */
		   window_t *windat,        /* Window data           */
                   double *tr,              /* Tracer array          */
                   double trem,/* Time remaining in the sub-step loop */
                   int tn                   /* Tracer number         */
  )
{
  int cc, c;                    /* Sparse counter / coordinate */

  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    /* tr[c]=windat->tr_wc_as[tn][cc]; */

    tr[c] = (windat->dt - trem) * (windat->tr_wc_ae[tn][cc] -
                                   windat->tr_wc_as[tn][cc]) /
      window->taux_t[cc] + windat->tr_wc_as[tn][cc];

  }
}

/* END set_multidt_t()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion                                              */
/*-------------------------------------------------------------------*/
void hor_diffuse(geometry_t *window,  /* Window geometry             */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon,  /* Window constants            */
                 double *tr,          /* Tracer array                */
                 double *Fx           /* Horizontal direction flux   */
		 )
{
  int e, ee, es;                /* Sparse indices                    */
  int j, c1, c2, cp;            /* 2D cell corresponding to 3D cell  */
  double csx, csy;              /* Cross sectional areas             */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* w5 = x face cross sectional area                                */
  /* s1 = wet cells to process for tracers (calculated in advection) */

  /*-----------------------------------------------------------------*/
  /* Get the face centered diffusivity                               */
  /*
  for (ee = 1; ee <= window->a3_e1; ee++) {
    e = window->w3_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    wincon->w3[e] = 0.5 * (wincon->u1kh[c1] + wincon->u1kh[c2]);
  }
  */
  memcpy(wincon->w3, wincon->u2kh, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* For sigma, limit the horizontal diffusion to maintain           */
  /* monotinicity (Herzfeld, 1996, p136).                            */
  if (wincon->sigma) {
    double tm, dp, dm;
    for (ee = 1; ee <= window->a3_e1; ee++) {
      e = window->w3_e1[ee];
      es = window->m2de[e];
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      j = window->e2e[e][1];
      cp = window->c2c[j][c1];

      dp = tr[cp] - tr[c1];
      dm = tr[c1] - tr[c2];
      tm = max(dp, tr[c1]);
      tm = max(tm, dm);
      csx = wincon->Ds[window->m2d[cp]] * dp / window->h2au1[window->ep[es]] -
        wincon->Ds[window->m2d[c1]] * dm / window->h2au1[es];
      csx = tm * wincon->Hn1[es] / (windat->dt * csx);
      if (wincon->w3[e] > fabs(csx))
        wincon->w3[e] = fabs(csx);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the diffusive fluxes                                  */
  for (ee = 1; ee <= window->a3_e1; ee++) {
    e = window->w3_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];

    /* Get the cross sectional area of the cell faces.               */
    csx = windat->dzu1[e] * window->h1au1[es] * wincon->mdx[es];
    Fx[e] -=
      (csx * wincon->w3[e] * (tr[c1] - tr[c2]) / window->h2au1[es]);
  }
}

/* END hor_diffuse()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion                                              */
/*-------------------------------------------------------------------*/
void hor_diffuse_2d(geometry_t *window, /* Window geometry           */
		    window_t *windat,   /* Window data               */
		    win_priv_t *wincon, /* Window constants          */
                    double *tr,         /* Tracer array              */
                    double *Fx          /* Horizontal direction flux */
  )
{
  int e, ee, es;                /* Sparse indices                    */
  int j, c1, c2, cp;            /* 2D cell corresponding to 3D cell  */
  double csx, csy;              /* Cross sectional areas */

  /*-----------------------------------------------------------------*/
  /* Get the face centered diffusivity                               *
  /*
  for (ee = 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    wincon->w3[e] = 0.5 * (wincon->u1kh[c1] + wincon->u1kh[c2]);
  }
  */
  memcpy(wincon->w3, wincon->u2kh, window->sze * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Calculate the diffusive fluxes                                  */
  for (ee = 1; ee <= window->a2_e1; ee++) {
    e = window->w2_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];

    /* Get the cross sectional area of the cell faces.               */
    csx = windat->depth_e1[e] * window->h1au1[e];

    Fx[e] -=
      (csx * wincon->w3[e] * (tr[c1] - tr[c2]) / window->h1au1[e]);
  }
}

/* END hor_diffuse_2d()                                              */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Horizontal diffusion without cross products or metric terms.      */
/*-------------------------------------------------------------------*/
void hor_diffuse_simple(geometry_t *window, /* Window geometry       */
                        window_t *windat,   /* Window data           */
                        win_priv_t *wincon, /* Window constants      */
			double *tr,
			int mode
  )
{
  int c, cc, cs;              /* Sparse indices                      */
  int j, e, cp, cm;           /* 2D cell corresponding to 3D cell    */
  double *AH;                 /* e1 horizontal diffusion coefficient */
  double d1;                  /* Constant for x diffusion            */
  int *cells;                 /* Cells to process vector             */
  int vc;                     /* Size of cells[]                     */

  AH = wincon->u1kh;

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
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    cs = window->m2d[c];
    d1 = 0.0;
    for (j = 1; j < window->npe[cs]; j++) {
      cp = window->c2c[jo(j, window->npe[cs])][c];
      cm = window->c2c[j][c];
      e = window->m2de[window->c2e[j][c]];
      d1 +=
	wincon->u1kh[c] * (tr[cp] + tr[cm] - 2.0 * tr[c]) /
	(window->h1au1[e] * window->h1au1[e]);
    }
    tr[c] += windat->dt * d1;
  }
}

/* END hor_diffuse_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Vertical diffusion                                                */
/*-------------------------------------------------------------------*/
void vert_diffuse_3d(geometry_t *window, /* Window geometry          */
		     window_t *windat,   /* Window data              */
		     win_priv_t *wincon  /* Window constants         */
		     )
{
  int nn;                       /* Tracer index                      */
  int cc;                       /* Counter                           */
  int c, c2;                    /* 3D and 2D cell coordinate         */
  int zm1, zp1;                 /* Cell coordinate at k-1 and k+1    */
  int cs, cb;                   /* Surface & bottom cell coordinate  */
  int *ctp = wincon->i2;        /* Old surface cell coordinate       */
  int *cbt = wincon->i5;        /* Bottom cell coordinate            */
  int *cth = wincon->i4;        /* Thin layer locations              */
  int ntbdy;                    /* Number of tracers to advect       */
  int *tbdy;                    /* Tracers to advect                 */
  double *topflux;              /* Surface flux                      */
  double dt = windat->dt;       /* Time step for the window          */
  double dzdt;                  /* dz divided by dt                  */
  double *dzface = wincon->w7;  /* Cell thickness at the cell face   */
  double *dzcell = wincon->w9;  /* Cell centered cell thickness      */
  double *Kz = wincon->w10;     /* Vertical diffusivity              */
  double *Cm1 = wincon->w1;     /* Constant for implicit calculation */
  double *C = wincon->w2;       /* Constant for implicit calculation */
  double *Cp1 = wincon->w3;     /* Constant for implicit calculation */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* d1 = flux out of bottom layer                                   */
  /* d2 = flux out of top layer                                      */
  /* d3 = scaling = 1.0                                              */
  /* s1 = wet cells to process for tracers (calculated in advection) */
  memset(wincon->d1, 0, window->szcS * sizeof(double));
  memset(wincon->d2, 0, window->szcS * sizeof(double));
  memset(cth, 0, window->szcS * sizeof(int));
  memcpy(dzcell, wincon->dz, window->szc * sizeof(double));
  ntbdy = wincon->ntdif_v;
  tbdy = wincon->tdif_v;

  if (wincon->means & TRANSPORT) {
    if (windat->dttr == 0.0)
      return;
    memcpy(Kz, windat->Kzm, window->szc * sizeof(double));
    if (wincon->trsplit) {
      ntbdy = wincon->ntbdys;
      tbdy = wincon->tbdys;
    }
  } else
    memcpy(Kz, windat->Kz, window->szc * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Merge thin layers if required                                   */
  if (wincon->thin_merge) {
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      c = cth[cc] = ctp[cc];
      cb = wincon->i1[cc];
      if (c != cb && dzcell[c] < wincon->hmin) {
        cs = ctp[cc] = window->zm1[c];
        for (nn = 0; nn < ntbdy; nn++) {
          int n = tbdy[nn];
          double *tr = windat->tr_wc[n];
          tr[c] = tr[cs] = (tr[c] * dzcell[c] + tr[cs] * dzcell[cs]) /
            (dzcell[cs] + dzcell[c]);
        }
        dzcell[cs] += dzcell[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Calculate constants that are independent of the variable to be  */
  /* diffused in the implicit calculation (i.e. the elements of the  */
  /* tridiagonal matrix). The rhs is set up in the routine           */
  /* implicit_vdiff_tr().                                            */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = cs = ctp[cc];       /* Old surface cell coordinate           */
    c2 = window->m2d[c];    /* 2D cell location corresponding to c   */
    zm1 = window->zm1[c];

    /* Single layer case (i.e. surface lies in the bottom layer)     */
    if (zm1 == window->zm1[zm1] || !dzcell[zm1])
      continue;

    /* Get the cell thicknesses. The sparse coordinate, c, lies in   */
    /* the sediment at the end of the loop.                          */
    while (c != zm1 && dzcell[c]) {
      dzface[c] = 0.5 * (dzcell[zm1] + dzcell[c]);
      /* SIGMA : Precondition mixing coefficient and the cell for    */
      /* the sigma calculation.                                      */
      dzcell[c] *= wincon->Ds[c2];
      Kz[c] /= wincon->Ds[c2];
      c = zm1;
      zm1 = window->zm1[c];
    }

    cb = c = window->zp1[c]; /* Set cb to the bottom cell coordinate */
    zp1 = window->zp1[c];    /* Layer above the bottom               */
    cbt[cc] = cb;

    /* Set up tri-diagonal set of equations.                         */
    /* Bottom layer.                                                 */
    dzdt = dzcell[c] / dt;
    Cm1[c] = 0.0;
    Cp1[c] = -Kz[zp1] / dzface[zp1];
    C[c] = dzdt - Cp1[c];

    /* Mid-water layers */
    while (zp1 != cs) {
      c = zp1;
      zp1 = window->zp1[c];
      dzdt = dzcell[c] / dt;
      Cm1[c] = -Kz[c] / dzface[c];
      Cp1[c] = -Kz[zp1] / dzface[zp1];
      C[c] = dzdt - Cm1[c] - Cp1[c];
    }

    /* Surface layer                                                 */
    c = cs;
    dzdt = dzcell[c] / dt;
    Cm1[c] = -Kz[c] / dzface[c];
    Cp1[c] = 0.0;
    C[c] = dzdt - Cm1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Tracer loop                                                     */
  for (nn = 0; nn < ntbdy; nn++) {
    int n = tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values                 */
    double *Splus = NULL;

    /* Set the surface boundary condition                            */
    topflux = wincon->d2;
    memset(topflux, 0, window->szcS * sizeof(double));
    if (n == windat->tno) {
      if (wincon->heatflux & (ADVANCED | INVERSE | NET_HEAT | COMP_HEAT | COMP_HEAT_MOM)) {
	/* Do swr data assimilation if required */
	if (wincon->compatible & V1562) {
	  set_sbc(window, windat, wincon);
	  Splus = wincon->w8;
	  memset(Splus, 0, window->szc * sizeof(double));
	  topflux = windat->heatf;
	} else {
	  Splus = wincon->w8;
	  memset(Splus, 0, window->szc * sizeof(double));
	  calc_swr(window, windat, wincon, Splus, 0);
	  /*
    if(window->wn==6)printf("a %f %f %e\n",windat->days,Splus[2383],windat->heatf[2383]);
    if(window->wn==6)printf("b %f %f %e\n",windat->days,Splus[2360],windat->heatf[2360]);
	  */
	  topflux = windat->heatf;
	}
      }
      if (wincon->heatflux & (SURF_RELAX|AVHRR|GHRSST))
        surf_relax(window, windat, wincon);
    }
    if (n == windat->sno) {
      /* Convert to correct units (e.g. Eqn 3.3, Introduction to     */
      /* Physical Oceanography, G.L. Mellor (1996)).                 */
      if (wincon->saltflux & (ADVANCED | BULK)) {
        for (cc = 1; cc <= window->b2_t; cc++) {
          c = window->w2_t[cc];
          topflux[c] = windat->nsfd[c] * windat->sal[c];
        }
      }
    }
    /* Generic surface fluxes                                        */
    if (wincon->sflux[n] >= 0) {
      memcpy(topflux, windat->tr_wcS[wincon->sflux[n]], window->szcS * sizeof(double));
    }

    /* No surface exchanges for thin layers                          */
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      if (windat->eta[c] - window->botz[c] <= wincon->hmin) topflux[c] = 0.0;
    }

    implicit_vdiff_tr(window, windat, wincon, tr, Kz, dzcell, dzface,
                      wincon->d1, topflux, ctp, cbt, wincon->vcs, Splus,
                      NULL, wincon->one, C, Cp1, Cm1);
  }
}

/* END vert_diffuse_3d()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implicit vertical diffusion scheme                                */
/*-------------------------------------------------------------------*/
void implicit_vdiff_tr(geometry_t *window, /* Window geometry        */
		       window_t *windat,   /* Window data            */
		       win_priv_t *wincon, /* Window constants       */
                       double *var,        /* Variables to diffuse   */
                       double *Kz,      /* Mixing coefficient values */
                       double *dzcell,  /* Cell thicknesses          */
                       double *dzface,  /* Thicknesses at cell face  */
                       double *fb,      /* Flux out of bottom        */
                       double *ft,      /* Flux out of top           */
                       int *ctp,        /* Cells to process          */
                       int *cbt,        /* Bottom coordinate         */
                       int vcs,         /* Last surface cell index   */
                       double *Splus,   /* Positive part of source   */
                       double *Sminus,  /* Negative part of source   */
                       double *scale,   /* SIGMA : (depth) scaling   */
                       double *C, double *Cp1, double *Cm1)
{
  int c, k;                  /* Cell coordinate                      */
  int cs, ks;                /* Surface cell coordinate              */
  int cb, kb;                /* Bottom cell coordinate               */
  int zm1;                   /* Cell cell below c                    */
  int cc;                    /* Cell coordinate counter              */
  int c2;                    /* 2D cell corresponding to 3D location */
  double dt = windat->dt;    /* Time step for the window             */
  double dzdt;               /* dz / dt                              */
  double div;                /* Constant                             */

  /*-----------------------------------------------------------------*/
  /* Set pointers.                                                   */
  /* Note: the 3D work arrays wincon->w# could be used for the       */
  /* dummy arrays below, but execution speed is considerably faster  */
  /* when work array access is sequential in memory, hence the       */
  /* mapping to a contiguous vertical 1D work array wincon->v#.      */
  int *cth = wincon->i4;
  double *rhs = wincon->v1;
  double *sol = wincon->v2;
  double *ud = wincon->v3;
  double *B = wincon->v4;
  double *Bm1 = wincon->v5;
  double *Bp1 = wincon->v6;
  double *dz = wincon->v7;

  /* Loop therough the surface cells in this window                  */
  for (cc = 1; cc <= vcs; cc++) {

    cs = c = ctp[cc];       /* Set cs to the surface cell coordinate */
    cb = cbt[cc];           /* Bottom cell coordinate                */
    c2 = window->m2d[c];    /* 2D cell location corresponding to c   */
    zm1 = window->zm1[c];
    ks = window->s2k[cs];
    kb = window->s2k[cb];

    /*---------------------------------------------------------------*/
    /* Single layer case (i.e. surface lies in the bottom layer)     */
    if (zm1 == window->zm1[zm1] || !dzcell[zm1]) {
      var[cs] +=
        windat->dt * (fb[c2] - ft[c2]) / max(dzcell[cs], wincon->hmin);
      if (c != cth[cc]) {
        var[cth[cc]] = var[cs];
      }
      continue;
    }

    /* Map the cell arrays to local coordinates.                     */
    c = cs;
    for (k = ks; k >= kb; k--) {
      B[k] = C[c];
      Bm1[k] = Cm1[c];
      Bp1[k] = Cp1[c];
      dz[k] = dzcell[c];
      c = window->zm1[c];
    }

    /*---------------------------------------------------------------*/
    /* Set up the rhs for the system of equations.                   */
    /* Bottom layer.                                                 */
    c = cb;                   /* Set c to the bottom cell coordinate */
    dzdt = dz[kb] / dt;
    rhs[kb] = dzdt * var[cb] + fb[c2];

    /* Mid-water layers                                              */
    c = window->zp1[cb];      /* Layer above the bottom              */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / dt;
      rhs[k] = dzdt * var[c];
      c = window->zp1[c];
    }

    /* Surface layer                                                 */
    dzdt = dz[ks] / dt;
    rhs[ks] = dzdt * var[cs] - ft[c2];

    /*---------------------------------------------------------------*/
    /* Add positive part of source terms, if specified               */
    if (Splus) {
      zm1 = window->zm1[c];
      c = cb;
      for (k = kb; k <= ks; k++) {
        rhs[k] += dz[k] * Splus[c];
        c = window->zp1[c];
      }
    }

    /* Add negative part of source terms, if specified               */
    if (Sminus) {
      c = cb;
      for (k = kb; k <= ks; k++) {
        B[k] += dz[k] * Sminus[c] / var[c];
        c = window->zp1[c];
      }
    }

    /*---------------------------------------------------------------*/
    /* Solve tridiagonal system                                      */
    div = B[kb];
    sol[kb] = rhs[kb] / div;
    for (k = kb + 1; k <= ks; k++) {
      ud[k] = Bp1[k - 1] / div;
      div = B[k] - Bm1[k] * ud[k];
      if (div == 0.0) {
        hd_quit_and_dump("Tracer diffusion;implicit_vdiff_tr: zero divisor\n");
	exit(0);
      }
      sol[k] = (rhs[k] - Bm1[k] * sol[k - 1]) / div;
    }

    /*---------------------------------------------------------------*/
    /* Update the variable                                           */
    c = cs;
    var[cs] += (sol[ks] - var[cs]) * scale[c2];
    for (k = ks - 1; k >= kb; k--) {
      c = window->zm1[c];
      sol[k] -= ud[k + 1] * sol[k + 1];
      var[c] += (sol[k] - var[c]) * scale[c2];
    }

    /*---------------------------------------------------------------*/
    /* Set the concentration in thin layers                          */
    c = cth[cc];
    if (c != cs) {
      var[c] += (sol[ks] - var[cs]) * scale[c2];
    }
  }
}

/* END implicit_vdiff_tr()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Vertical diffusion for the 2d mode includes surface fluxes        */
/*-------------------------------------------------------------------*/
void vert_diffuse_2d(geometry_t *window, /* Window geometry          */
		     window_t *windat,   /* Window data              */
		     win_priv_t *wincon  /* Window constants         */
		     )
{
  int nn;                       /* Tracer index                      */
  int cc;                       /* Counter                           */
  int c, c2, cb;                /* 3D and 2D cell coordinate         */
  double *topflux;              /* Surface flux                      */
  double top, bot;              /* Surface and bottom heights        */
  double Cv = 4e3;              /* Specific heat at constant volume  */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* d1 = flux out of top layer */
  memset(wincon->d1, 0, window->szcS * sizeof(double));

  for (nn = 0; nn < wincon->ntdif_v; nn++) {
    int n = wincon->tdif_v[nn];
    double *tr = windat->tr_wc[n];

    /* Set the surface boundary condition                            */
    if (n == windat->tno) {
      if (wincon->heatflux & (ADVANCED | INVERSE | NET_HEAT | COMP_HEAT | COMP_HEAT_MOM)) {
        /* Add the shortwave to the heatflux                         */
        if (windat->swr_attn) {
          for (cc = 1; cc <= window->b2_t; cc++) {
            c = window->w2_t[cc];
            windat->heatf[c] -= (windat->swr[c] / (Cv * windat->dens[c]));
          }
        }
        topflux = windat->heatf;
      } else
        topflux = wincon->d1;
      if (wincon->heatflux & (SURF_RELAX|AVHRR))
        surf_relax(window, windat, wincon);
    } else
      topflux = wincon->d1;

    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      top = windat->eta[c];
      bot = window->botz[c];
      tr[c] +=
        windat->dt * (0.0 - topflux[c]) / max(top - bot, wincon->hmin);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the water column tracer values equal to the surface         */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values                 */

    for (cc = 1; cc <= window->b2_t; cc++) {
      c = c2 = window->w2_t[cc];
      cb = window->bot_t[cc];
      while (c != cb) {
        c = window->zm1[c];
        tr[c] = tr[c2];
      }
    }
  }
}

/* END vert_diffuse_2d()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Do the diagnostics and interfaced library routines                */
/*-------------------------------------------------------------------*/
void auxiliary_routines(geometry_t *window, /* Window geometry       */
			window_t *windat,   /* Window data           */
			win_priv_t *wincon  /* Window constants      */
			)
{
  int nn;
  int cc, c, c2, cb;
  double h1, h2, diff;

  /*-----------------------------------------------------------------*/
  /* Reset dz and the cells to process to correspond to the updated  */
  /* elevation.                                                      */
  reset_dz(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the mixed layer depth diagnostic if required                */
  if (!(wincon->mixlayer & NONE))
    get_mixed_layer(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Flux diagnostics are calculated in advect_diffuse() using the   */
  /* actual fluxes used to update tracer concentration (i.e.         */
  /* advection scheme dependent).                                    */

  /*-----------------------------------------------------------------*/
  /* Calculate the flushing diagnostic if required                   */
  if (wincon->trflsh)
    total_mass(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the age tracer diagnostic if required                 */
  if (windat->agetr)
    calc_age(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the tracer percentiles if required                    */
  if (wincon->trperc >= 0)
    perc_diag(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the steric height if required                         */
  if (wincon->lnm != 0.0)
    steric(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the vorticity if required                             */
  if (!(wincon->vorticity & NONE))
    vorticity(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the monotinicity diagnostic if required               */
  if (wincon->monon >= 0)
    calc_monotonic(window, windat, wincon);
  
  /*-----------------------------------------------------------------*/
  /* Get diagnostic numbers if required                              */
  if (!(wincon->numbers & NONE) || !(wincon->numbers1 & NONE))
    diag_numbers(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get normalized vertical profile of a tracer if required         */
  if (wincon->nprof >= 0)
    nor_vert_prof(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the DHW if required                                         */
  if (wincon->ndhw) {
    for (nn = 0; nn < wincon->ndhw; nn++) {
      if (wincon->dhwf[nn] & DHW_NOAA)
	calc_dhd(window, windat, wincon, nn);
    }
  }

#if defined(HAVE_TRACERSTATS_MODULE)
  tracerstats_prestep(window,1);
#endif

  /*-----------------------------------------------------------------*/
  /* Evaluate the tracer decay                                       */
  for (nn = 0; nn < wincon->ntdec; nn++) {
    int n = wincon->tdec[nn];
    tracer_decay(window, windat, wincon, wincon->dectr[n],
                 windat->tr_wc[n], &wincon->trinfo_3d[n]);
  }

  /*-----------------------------------------------------------------*/
  /* Do the temperature and salinity relaxation                      */
  /* DA OFF (NONE)  always                                           */
  /* DA ON  (DO_DA) reanalysis phase                                 */ 
  if (!master->da || (master->da & (NONE|DO_DA)))
    do_ts_relax(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Add temperature and salinity tracer increments                  */
  do_ts_increment(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the tracer alert diagnostics                          */
  alerts_w(window, TRACERS);

  /*-----------------------------------------------------------------*/
  /* Get the total mass in regions                                   */
  region_mass(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get total mass, heat and salt                                   */
  if (wincon->totals)
    mass_diag(window, windat, wincon);

#if defined(HAVE_TRACERSTATS_MODULE)
  tracerstats_prestep(window,2);
#endif

  /*-----------------------------------------------------------------*/
  /* Do the wave routines                                            */
#if defined(HAVE_WAVE_MODULE)
  if(wincon->do_wave)
    wave_interface_step(window);
#endif

#if defined(HAVE_TRACERSTATS_MODULE)
  tracerstats_prestep(window,6);
#endif

  /*-----------------------------------------------------------------*/
  /* Do the sediment transport routines                              */ 
#if defined(HAVE_SEDIMENT_MODULE)
  if(wincon->do_sed & LIB_DO) {
    TIMING_SET;
    sed_step(window);
    TIMING_DUMP_WIN(3,"   sed_step", window->wn);

    /* Get the total mass after sediment transport                   */
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      int n = wincon->tbdy[nn];
      region_mass_tr(window, windat, wincon, 1.0, n, RG_SED);
#if defined(HAVE_ECOLOGY_MODULE)
      if (wincon->do_eco)
	region_mass_tr(window, windat, wincon, -1.0, n, RG_ECO);
#endif
    }
  }
#endif

#if defined(HAVE_TRACERSTATS_MODULE)
  tracerstats_prestep(window,3);
#endif

  /*-----------------------------------------------------------------*/
  /* Do the ecological routines                                      */
#if defined(HAVE_ECOLOGY_MODULE)
  // Timing is in eco_step
  eco_step(window);

  /* Get the total mass after ecology                                */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    region_mass_tr(window, windat, wincon, 1.0, n, RG_ECO);
  }
#endif

#if defined(HAVE_TRACERSTATS_MODULE)
  /*
   * 0 is the default.
   * 4 and 5 are there so you can do tracerstats of tracerstats as a
   * pre- or post- respectively.
   */
  tracerstats_prestep(window,4);

  TIMING_SET;
  tracerstats_prestep(window, 0);
  TIMING_DUMP_WIN(3,"   tracerstats", window->wn);

  tracerstats_prestep(window,5);
#endif

  /* Step 0 must always be the last call to tracerstats()            */
#if defined(HAVE_TRACERSTATS_MODULE)
  /*tracerstats_prestep(window,0);*/
#endif

  if (wincon->waves & BOT_STR)
    memcpy(wincon->Cd, windat->wave_Cd, window->szcS * sizeof (double));

#if defined(HAVE_ECOLOGY_MODULE) && defined(HAVE_SEDIMENT_MODULE)
  /* Reset a no-gradient above the surface if required               */
  if (wincon->do_sed & LIB_DO || wincon->do_eco)
    tr_set_surf(window, windat, wincon);
#elif defined(HAVE_SEDIMENT_MODULE)
  /* Reset a no-gradient above the surface if required               */
  if (wincon->do_sed & LIB_DO)
    tr_set_surf(window, windat, wincon);
#endif
}

/* END auxiliary_routines()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set sub-surface tracers and density to the surface     */
/* value in the 2D mode.                                             */
/*-------------------------------------------------------------------*/
void wc_tracer_2d(geometry_t *window,       /* Window geometry       */
		  window_t *windat,         /* Window data           */
		  win_priv_t *wincon        /* Window constants      */
		  )
{
  int nn;                       /* Tracer index                      */
  int cc;                       /* Counter                           */
  int c, c2, cb;                /* 3D and 2D cell coordinate         */

  /* Density                                                         */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = c2 = window->w2_t[cc];
    cb = window->bot_t[cc];
    while (c != cb) {
      c = window->zm1[c];
      windat->dens[c] = windat->dens[c2];
    }
  }

  /* Tracers                                                         */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values                 */

    for (cc = 1; cc <= window->b2_t; cc++) {
      c = c2 = window->w2_t[cc];
      cb = window->bot_t[cc];
      while (c != cb) {
        c = window->zm1[c];
        tr[c] = tr[c2];
      }
    }
  }
}

/* END wc_tracer_2d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set all the tracer values on all open boundaries in    */
/* a given window.                                                   */
/*-------------------------------------------------------------------*/
void bdry_tracer(geometry_t *window,        /* Window geometry       */
		 window_t *windat,          /* Window data           */
		 win_priv_t *wincon         /* Window constants      */
  )
{
  int n, nn;                    /* Counters                          */
  int tn;                       /* Tracer number                     */
  open_bdrys_t **open = window->open;

  /*
  if (wincon->trasc == FFSL)
    prep_semi_lagrange_atc(window);
  */

  for (n = 0; n < window->nobc; n++) {
    for (nn = 0; nn < wincon->ntbdy; nn++) {
      tn = wincon->tbdy[nn];
      set_OBC_tr(tn, window, windat, wincon, open[n]);
    }
  }

}

/* END bdry_tracer()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set the boundary conditions on the tracers                        */
/*-------------------------------------------------------------------*/
void set_OBC_tr(int tn,             /* Tracer number                 */
		geometry_t *window, /* Window geometry               */
		window_t *windat,   /* Window data                   */
		win_priv_t *wincon, /* Window constants              */
                open_bdrys_t *open  /* OBC structure                 */
  )
{
  int c, cc;                    /* Cell coorindate / counter         */
  int c1, c2;                   /* Cell coordinates                  */
  int e, ee, es;                /* Edge coordinates / counters       */
  int bcond;                    /* Boundary condition type           */
  double *tr;                   /* Tracer array                      */
  double *newval;               /* New tracer values                 */
  int *imap = NULL;             /* Interior cell map                 */
  double *vel = NULL;           /* Velocity on the boundary          */
  double *hat = NULL;           /* Grid spacing on the boundary      */
  int rlxn = 0;                 /* Flow relaxation zone              */
  scale_details_t scale = open->sdata_t[tn];

  /* Initialise the pointers                                         */
  bcond = open->bcond_tra[tn];
  rlxn = open->relax_zone_tra[tn];
  tr = windat->tr_wc[tn];
  newval = wincon->w1;
  if (open->options & OP_TILED) return;

  /*-----------------------------------------------------------------*/
  /* Set a no gradient condition on tracer boundaries                */
  if (bcond & NOGRAD) {
    if (open->trpc[tn]) {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        c1 = open->oi1_t[cc];
        if (window->cellz[c] > open->trpc[tn])
          newval[c] = tr[c1];
      }
    } else {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        c1 = open->oi1_t[cc];
        newval[c] = tr[c1];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* No action                                                       */
  else if (bcond & NOTHIN) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] = tr[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Cyclic                                                          */
  else if (bcond & CYCLIC) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      c1 = open->cyc_t[cc];
      newval[c] = tr[c1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Clamped                                                         */
  else if (bcond & CLAMPD) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] = open->clampv[tn];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Statistical prescription                                        */
  /* perct = -2 : last value of subset (cyclic over subset)          */
  /* perct = -1 : mean value over subset                             */
  /* perct = 0 : minimum of subset                                   */
  /* perct = 100 : maximum of subset                                 */
  /* 0 < perct < 100 : percentile value of subset                    */
  /* vel_frac = velocity factor defining subset; subset limit is the */
  /* index si where (for a u1 boundary) :                            */
  /* |u1[k][j][si]-u1[k][j][i]| > |vel_frac*u1[k][j][i]|             */
  /* This boundary condition prescribes a new boundary value         */
  /* independently from the existing value. Reverts to a no gradient */
  /* condition if the subset contains <= 1 points.                   */
  else if (bcond & STATIS) {
    double perct = 20.0;
    double vel_frac = 0.1;
    double bvel, fvel;
    int e, eo, es, ci, nn, dir;
    double *vals;

    vel = windat->u1;

    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->obc_e2[ee];
      c2 = window->m2d[c];
      dir = open->ceni[ee];
      imap = (dir) ? window->em : window->ep;

      /* Get the range over which to compute the statistics. This is */
      /* taken as the cell where velocity changes by a certain       */
      /* fraction of the boundary velocity. Note : the boundary      */
      /* value is not included in the subset.                        */
      e = es = window->c2e[open->ini[ee]][c];
      bvel = vel[e];
      fvel = fabs(vel_frac * bvel);
      ci = 0;
      while (fabs(vel[e] - bvel) < fvel && e != imap[e]) {
        e = imap[e];
        ci += 1;
      }
      eo = e;
      /* Store the tracer values over the range in the array vals[]  */
      if (ci > 0) {
        vals = d_alloc_1d(ci);
        nn = 0;
	e = es;
        while (e != eo) {
	  e = imap[e];
          vals[nn] = tr[window->e2c[e][dir]];
          nn++;
        }
        /* Get the statistical value                                 */
        newval[c] = percs(vals, ci, perct);
        d_free_1d(vals);
      }
      /* Set a no gradient condition                                 */
      else {
	cc = open->e2c_e1[ee];
        c1 = open->oi1_t[cc];
        newval[c] = tr[c1];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Least squares linear interpolation                              */
  else if (bcond & LINEXT) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      int c = open->obc_t[cc];
      newval[c] = 0.0;
    }
    for (ee = 1; ee <= open->no3_e1; ee++) {
      int ei = open->obc_e2[ee];
      int cc = open->e2c_e1[ee];
      int ei0 = open->nmape[ee][ei];
      int ei1 = open->nmape[ee][ei0];
      int ei2 = open->nmape[ee][ei1];
      int ei3 = open->nmape[ee][ei2];
      double f = open->nepc[cc]; 
      newval[ei] += f * bc_leastsq(tr[ei3], tr[ei2], tr[ei1], tr[ei0]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Polynomial extrapolation (2nd order)                            */
  else if (bcond & POLEXT) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      int c = open->obc_t[cc];
      newval[c] = 0.0;
    }
    for (ee = 1; ee <= open->no3_e1; ee++) {
      int ei = open->obc_e2[ee];
      int cc = open->e2c_e1[ee];
      int ei0 = open->nmape[ee][ei];
      int ei1 = open->nmape[ee][ei0];
      int ei2 = open->nmape[ee][ei1];
      int ei3 = open->nmape[ee][ei2];
      double f = open->nepc[cc]; 
      newval[ei] += f * bc_polint(tr[ei3], tr[ei2], tr[ei1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Custom specification                                            */
  if (bcond & CUSTOM) {
    if (open->trpc[tn]) {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        if (window->cellz[c] <= open->trpc[tn])
          newval[c] =
            bdry_value_w(window, windat, wincon, open, &open->bdata_t[tn],
                         c, window->m2d[c], windat->t);
      }
    } else {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        newval[c] =
          bdry_value_w(window, windat, wincon, open, &open->bdata_t[tn], c,
                       window->m2d[c], windat->t);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* File specification                                              */
  if (bcond & FILEIN) {
    if (open->trpc[tn]) {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        if (window->cellz[c] <= open->trpc[tn])
          newval[c] = open->t_transfer[open->trm[tn]][cc];
      }
    } else {
      for (cc = 1; cc <= open->no3_t; cc++) {
        c = open->obc_t[cc];
        newval[c] = open->t_transfer[open->trm[tn]][cc];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Construct a vertical profile using surface and bottom values.   */
  /* This is constructed using an exponential profile mathched to    */
  /* an inverted exponential profile, where the matching occurs at   */
  /* the mixed layer depth having a mid-point profile value.         */
  if (bcond & PROFIL) {
    double sf = 0.1;      /* Exponential scaling factor              */
    double sf1 = 0.2;     /* Exponential scaling factor              */
    int nin = 10;         /* Number of cells inside to get mld       */
    double minml = 0.1;   /* Minimum mixed layer depth               */
    double mld;           /* Mixed layer depth                       */
    double d1, d4;        /* Dummy                                   */
    double tv;            /* Surface value                           */
    double bv;            /* Bottom value                            */
    double depth;         /* Depth in water column                   */
    double cf, cf1;       /* Scaling factors                         */
    int cs;               /* Surface cell coordinate                 */
    int cb;               /* Bottom cell coordinate                  */
    int ci;               /* Interior profile location               */
    imap = open->nmap;

    if (!windat->mixl)
      hd_quit("tracers: 'PROFIL' boundary requires 'MIX_LAYER' diagnostic set.\n");

    for (cc = 1; cc <= open->no2_t; cc++) {
      c = cs = ci = open->obc_t[cc];
      cb = open->bot_t[cc];

      /* Get the surface value                                       */
      tv = bv = newval[cs];
      /* Get the bottom value                                        */
      while (c <= cb) {
	bv = newval[c];
	c = window->zm1[c];
      }
      /* Get the mixed layer depth nin cells into the interior       */
      for (c1 = 1; c1 <= nin; c1++)
	ci = imap[ci];
      mld = max(minml, windat->mixl[ci]);

      /* Get the profile mid-point value                             */
      cf = sf * mld;
      cf1 = sf1 * mld;
      /*d1 = tv - exp(mld / cf) / exp(mld / cf - 1.0) - bv;*/
      d4 = mld / cf - log(2.0 * exp(mld / cf) / fabs(tv - bv));

      /* Calculate the profile                                       */
      c = cs;
      if(tv > bv) {
	d1 = tv - exp(mld / cf) / exp(mld / cf - d4) - bv;
	while (c <= cb) {
	  depth = window->cellz[cs] - window->cellz[c];
	  if(depth <= mld)
	    newval[c] = tv - exp(depth / cf) / exp( mld / cf - d4);
	  else
	    newval[c] = bv + d1 * exp(-(depth - mld) / (mld / cf1 - 4.0));
	  c = window->zm1[c];
	}
      }
      else {
	d1 = tv + exp(mld / cf) / exp(mld / cf - d4) - bv;
	while (c <= cb) {
	  depth = window->cellz[cs] - window->cellz[c];
	  if(depth <= mld)
	    newval[c] = tv + exp(depth / cf) / exp(mld / cf - d4);
	  else
	    newval[c] = bv + d1 * exp(-(depth - mld) / (mld / cf1 - 4.0));
	  c = window->zm1[c];
	}
      }
      /*
      while (c <= cb) {
	if(window->cellz[c] <= mld)
           newval[c] = tv - exp(window->cellz[c]/cf) / exp(mld/cf-1.0);
          else
	    newval[c] = bv + d1*exp(-(window->cellz[c]-mld) / (mld/cf1-4.0));
	c = window->zm1[c];
      }
      */
    }
  }

  /*-----------------------------------------------------------------*/
  /* Construct a vertical profile using surface and bottom values.   */
  /* This is constructed the using the actual density profile nin    */
  /* cells into the interior from the boundary. A normalized profile */
  /* is made using the surface density at this point and bottom      */
  /* density nin cells into the interior from the deepest location   */
  /* on the boundary (only surface and deepest bottom tracer values  */
  /* are supplied). The actual tracer profile is then reconstructed  */
  /* from the normalized profile.                                    */
  if (bcond & DEPROF) {
    int nin = 5;          /* Number of cells inside to get mld       */
    double tv;            /* Surface value                           */
    double bv;            /* Bottom value                            */
    double dens_t;        /* Surface density                         */
    double dens_b;        /* Bottom density                          */
    double maxt;          /* Maximum tracer value                    */
    double mint;          /* Minimum tracer value                    */
    int cs;               /* Surface cell coordinate                 */
    int cb;               /* Bottom cell coordinate                  */
    int ci;               /* Interior profile location               */
    int k;
    double *s = wincon->v1;
    imap = open->nmap;

    /* Get the profile location nin cells into the interior from the */
    /* deepest point on the boundary.                                */
    /* Surface density.                                              */
    ci = open->obc_t[open->cdeep];
    for (c1 = 1; c1 <= nin; c1++)
      ci = imap[ci];
    dens_t = windat->dens_0[ci];

    /* Bottom density.                                               */
    c = ci;
    cb = window->bot_t[window->c2cc[open->obc_t[open->cdeep]]];
    for (k = window->nz - 1; k > window->s2k[cb]; k--)
      c = window->zm1[c];
    c = min(c, window->bot_t[window->c2cc[ci]]);
    dens_b = windat->dens_0[c];

    for (cc = 1; cc <= open->no2_t; cc++) {
      c = cs = ci = open->obc_t[cc];
      cb = open->bot_t[cc];

      /* Get the profile location nin cells into the interior and    */
      /* surface density at this location.                           */
      for (c1 = 1; c1 <= nin; c1++)
	ci = imap[ci];
      dens_t = windat->dens_0[ci];

      /* Make the normalized profile. Note; if the density at any    */
      /* griven location ci is greater than the bottom density nin   */
      /* cells interior to the deepest boundary cell, then s[k] > 1. */
      /* This may progressively feed back on the boundary salinity   */
      /* and lead to maximum bound violations.                       */
      c = ci;
      k = window->s2k[ci];
      while (c != window->zm1[c]) {
	k = window->s2k[c];
	s[k] = (windat->dens_0[c] - dens_t) / (dens_b - dens_t);
	c = window->zm1[c];
      }
      for(k--; k >= 0; k--)
	s[k] = s[k+1];

      /* Get the surface and bottom tracer values                    */
      tv = newval[cs];
      bv = newval[cb];
      maxt = max(tv, bv);
      mint = min(tv, bv);

      /* Reconstruct the profile. Note; truncate the tracer value to */
      /* the maximum measured to avoid maximum bound violations.     */
      c = cs;
      while (c <= cb) {
	newval[c] = min(tv + s[window->s2k[c]] * (bv - tv), maxt);
	newval[c] = max(newval[c], mint);
	c = window->zm1[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Fit a value to the surface and bottom and a profile by          */
  /* incrementing this value by a scaled vertical density difference.*/ 
  if (bcond & DESCAL && scale.type & TRSC_DEN) {
    int zm1, zp1, cs, cb, cm;
    double dt, db, vs, ov;
    /* Fix the tracer to a mid-depth value from file and construct   */
    /* the profile upwards and downwards.                            */
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = cs = open->obc_t[cc];
      /* Find the layer corresponding to v1                          */
      while (c != window->zm1[c] && fabs(window->cellz[c]) < fabs(scale.val)) {
	c = window->zm1[c];
      }
      c = cm = window->zp1[c];
      /* Get the normalized density profile upwards                  */
      newval[c] = vs = newval[cm];
      zm1 = c;
      c = window->zp1[c];
      while (c != window->zp1[c]) {
	dt = windat->dens_0[c];
	db = windat->dens_0[zm1];
	newval[c] = vs - vs * scale.fact * 2.0 * (db - dt) / (dt + db);
	vs = newval[c];
	zm1 = c;
	c = window->zp1[c];
      }
      dt = windat->dens_0[c];
      db = windat->dens_0[zm1];
      newval[c] = vs - vs * scale.fact * 2.0 * (db - dt) / (dt + db);
      /* Get the normalized density profile downwards                */
      zp1 = cm;
      c = window->zm1[cm];
      vs = newval[cm];
      while (c != window->zm1[c]) {
	if (scale.flag & TRSC_CPY) {
	  zp1 = c;
	  c = window->zm1[c];
	  continue;
	}
	dt = windat->dens_0[zp1];
	db = windat->dens_0[c];
	ov = newval[c];
	newval[c] = (db == 0.0) ? vs : vs + vs * scale.fact * 2.0 * (db - dt) / (dt + db);
	if (scale.flag & TRSC_TRN && newval[c] > ov) newval[c] = ov;
	vs = newval[c];
	zp1 = c;
	c = window->zm1[c];
      }
    }
    /* Clip                                                          */
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] = max(newval[c], wincon->trinfo_3d[tn].valid_range_wc[0]);
      newval[c] = min(newval[c], wincon->trinfo_3d[tn].valid_range_wc[1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Scaling                                                         */
  if (scale.type == (TRSC_SUM | TRSC_NUM)) {
    double fact = scale.fact;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] += fact;
    }
  } else if (scale.type == (TRSC_SUM | TRSC_TRA)) {
    int trn = scale.ntr;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] += windat->tr_wc[trn][c];
    }
  } else if (scale.type == (TRSC_PCT | TRSC_NUM)) {
    double fact = scale.fact;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] *= fact;
    }
  } else if (scale.type == (TRSC_PCT | TRSC_TRA)) {
    int trn = scale.ntr;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      newval[c] *= windat->tr_wc[trn][c];
    }
  }
  /* Scale the ghost cells for TRCONC                                */
  if (bcond & TRCONC) {
    int m;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->ogc_t[ee];
      for (m = 1; m <= open->bgz; m++) {
	if (scale.type == (TRSC_SUM | TRSC_NUM)) {
	  double fact = scale.fact;
	  tr[c] += fact;
	} else if (scale.type == (TRSC_SUM | TRSC_TRA)) {
	  int trn = scale.ntr;
	  tr[c] += windat->tr_wc[trn][c];
	} else if (scale.type == (TRSC_PCT | TRSC_NUM)) {
	  double fact = scale.fact;
	  tr[c] *= fact;
	} else if (scale.type == (TRSC_PCT | TRSC_TRA)) {
	  int trn = scale.ntr;
	  tr[c] *= windat->tr_wc[trn][c];
	}
	c = open->omape[ee][c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set an upstream advection condition on tracer boundaries        */
  if (bcond & UPSTRM) {
    double sc;
    double f, v;
    int ei;
    vel = windat->u1;
    hat = window->h2au1;

    if (wincon->trasc == LAGRANGE) {
      /* Characteristic method for use with LAGRANGE                 */

      memset(open->dum, 0, (open->no3_t + 1) * sizeof(double));
      for (ee = 1; ee <= open->no3_e1; ee++) {
	int cg, j = open->outi[ee];
	e = open->obc_e1[ee];
	es = window->m2de[e];
	c = open->obc_e2[ee];
	c2 = window->m2d[c];
	cc = open->e2c_e1[ee];
	open->dum[cc] += (window->eSc[j][c2] * vel[e] * windat->dt / 
			  window->h2au1[es]);
	cg = open->ogc_t[ee];
	tr[cg] = newval[c];
      }
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c2 = window->m2d[c];
	f = min(open->dum[cc] * open->nepc[cc], 1.0);
	f = max(f, -1.0);
	newval[c] = (f < 0.0) ?
	  tr[c] + f * (tr[c] - newval[c]) : tr[c];
      }
    } else {
      
      /* Velocities used are the mean of boundary & interior         */
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
	ei = open->oi1_e1[ee];
	c = open->obc_e2[ee];
	cc = open->e2c_e1[ee];
	c1 = open->oi1_t[cc];
	sc = open->dir[ee];
	v = get_upvel(vel, sc, e, ei, open->upmeth);
	c2 = window->m2d[c];
	f = open->nepc[cc];
	newval[c] += f * ( - windat->dt / hat[c2]
	  * (0.5 * (v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (v - fabs(v)) * (newval[c] - tr[c])));
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Relaxation to external solutions                                */
  if (rlxn) {
    double alpha;
    int bn;
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      cc = open->e2c_e1[ee];
      c = open->oi1_t[cc];
      bn = 2;
      while (c != open->nmape[ee][c] && bn <= rlxn) {
        alpha = 1 - tanh(0.5 * (double)(bn - 1));
        tr[c] = alpha * newval[c] + (1.0 - alpha) * tr[c];
        c = open->nmape[ee][c];
        bn++;
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Update the solution */
  if (!(bcond & (NOTHIN|TRCONC|TRFLUX|TRCONF)) || open->options & OP_OWRITE) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      tr[c] = newval[c];
    }
  }
}

/* END set_OBC_tracer()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate upstream advection tracer value              */
/*-------------------------------------------------------------------*/
void upstrm(geometry_t *window,             /* Window geometry       */
	    window_t *windat,               /* Window data           */
	    win_priv_t *wincon,             /* Window constants      */
	    open_bdrys_t *open,             /* Open boundary         */
	    double *newval,                 /* New tracer value      */
	    double *tr,                     /* Tracer value          */
	    int mode)
{
  int c, cc, c1, c2;
  int e, ei, ee, es;
  double sc;
  double f, v;
  double *vel = NULL, *hat = NULL;

  vel = windat->u1;
  hat = window->h1au1;

  if (wincon->trasc == LAGRANGE) {
    /* Characteristic method for use with LAGRANGE                   */
    memset(open->dum, 0, (open->no3_t + 1) * sizeof(double));
    for (ee = 1; ee <= open->no3_e1; ee++) {
      int j = open->outi[ee];
      e = open->obc_e1[ee];
      es = window->m2de[e];
      c = open->obc_e2[ee];
      cc = open->e2c_e1[ee];
      open->dum[cc] += min(fabs(vel[e]) * windat->dt / 
			   window->h2au1[es], 1.0);
    }
    for (cc = 1; cc <= open->no3_t; cc++) {

      c = open->obc_t[cc];
      c2 = window->m2d[c];
      f = open->dum[cc] * open->nepc[cc];
      newval[c] = (f > 0.0) ?
	tr[c] - f * (tr[c] - newval[c]) : tr[c];
    }
  } else {    
    /* Tracer in boundary ghost cell is updated using face velocity  */
    if (mode & GHOST) {
      double d1;
      if (open->rlen) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = open->obc_e1[ee];
	  c = open->ogc_t[ee];
	  c1 = open->obc_e2[ee];
	  v = vel[e];
	  c2 = window->m2d[c1];
	  sc = -(double)open->dir[ee];
	  newval[c] +=  - windat->dt *
	    (0.5 * (sc * v + fabs(v)) * (newval[c] - tr[open->nmape[ee][c]]) / hat[c2]
	       + (sc * v - fabs(v)) * (tr[open->omape[ee][c]] - newval[c]) / open->rlen);
	}
      } else {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = open->obc_e1[ee];
	  c = open->ogc_t[ee];
	  c1 = open->obc_e2[ee];
	  v = vel[e];
	  c2 = window->m2d[c1];
	  sc = -(double)open->dir[ee];
	  newval[c] +=  - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (newval[c] - tr[open->nmape[ee][c]])
	       + 0.5 * (sc * v - fabs(v)) * (tr[open->omape[ee][c]] - newval[c]));
	}
      }
    } else if (mode & VANLEER) {
      set_map_t(window);
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	wincon->w6[c] = vel[c] * windat->dt / hat[window->m2d[c]];
      }
      for (cc = 1; cc <= open->no3_e1; cc++) {
	int cp;
	int e = open->obc_e1[cc];
	int c1 = window->e2c[e][0];
	int c2 = window->e2c[e][1];
	int j = window->e2e[e][0];
	int cm = window->c2c[j][c2];
	j = window->e2e[e][1];
	cp = window->c2c[j][c1];

	van_leer_do(newval, tr, windat->u1, wincon->w6, e, c1, cp, c2, cm);
      }
      reset_map_t_all(window);
    } else {
      /* Velocities used are the mean of boundary & interior           */
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
	ei = open->oi1_e1[ee];
	c = open->obc_e2[ee];
	cc = open->e2c_e1[ee];
	c1 = open->oi1_t[cc];
	sc = -(double)open->dir[ee];
	v = get_upvel(vel, sc, e, ei, mode);
	c2 = window->m2d[c];
	f = open->nepc[cc];
	newval[c] += f * ( - windat->dt / hat[c2]
	  * (0.5 * (v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (v - fabs(v)) * (newval[c] - tr[c])));
      }
    }
  }
}

/* END upstrm()                                                      */
/*-------------------------------------------------------------------*/


double get_upvel(double *vel, double sc, int e, int ei, int mode)
{
  double v;
  if (mode == CENTER)
    v = 0.5 * (vel[e] + vel[ei]);
  else if (mode == FACE)
    v = vel[e];
  else if (mode == ADAPTIVE)
    v = - sc * 0.5 * (vel[e] + vel[ei]) > 0.0 ?
      vel[e] : vel[ei];
  else
    v = vel[ei];
  return(sc * v);
}

/*-------------------------------------------------------------------*/
/* Routine to save the boundary tracer values to buffer and reset    */
/*-------------------------------------------------------------------*/
void save_OBC_tr(geometry_t *window,  /* Processing window           */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon , /* Window constants            */
		 double *Fx,          /* Flux array                  */
		 double *tr,          /* Tracer                      */
		 int tn,              /* Tracer number               */
		 int mode
		 )
{
  int n, c, cc, e, ee;

  /* Return if no tracers in this window have OBCs TRFLUX            */
  if (!mode) return;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (!(open->bcond_tra[tn] & (TRFLUX|TRCONC|TRCONF))) {
      if (mode == 1) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  open->dumn[tn][cc] = tr[c];
	}
      }
      if (mode == 2) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  tr[c] = open->dumn[tn][cc];
	}
      }
    } else {
      if (open->options & OP_NOHDIF) {
	if (mode == 3) {
	  for (ee = 1; ee <= open->no3_e1; ee++) {
	    e = open->obc_e1[ee];
	    open->dumn[tn][ee] = Fx[e];
	  }
	}
	if (mode == 4) {
	  for (ee = 1; ee <= open->no3_e1; ee++) {
	    e = open->obc_e1[ee];
	    Fx[e] = open->dumn[tn][ee];
	  }
	}
      }
    }
  }
}

/* END save_OBC_tr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the tracer flux divergence at open boundary      */
/* cells. Note: positive fluxes imply import of tracer, hence the    */
/* negative values for outward edges, where a positive value on      */
/* these edges usually implies tracer export.                        */
/* TRFLUX applies the prescribed flux directly.                      */
/* TRCONF multiplies a prescribed concentration by the volume flux,  */
/* scaling so that only inflow fluxes are modified. Note that mass   */
/* balance and UPSTRM options may be applied to this prescribed      */
/* value (in bdry_transfer_tr()).                                    */
/*-------------------------------------------------------------------*/
void reset_tr_OBCflux(geometry_t *window,  /* Processing window      */
	       window_t *windat,    /* Window data                   */
	       win_priv_t *wincon,  /* Window geometry / constants   */
	       double *dtracer,     /* Horizontal flux divergence    */
	       double *Fx,          /* Horizontal e1 fluxes          */
	       double dt,           /* Timestep                      */
	       int tn               /* Tracer number                 */
	       )
{
  int c, cc, e, ee, c2, n, m, nn, zm1, nc;

  /* Return if no tracers in this window have OBCs TRFLUX            */
  if (!(wincon->obctr[tn] & (TRFLUX|TRCONF))) return;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    double sgn;
    double *flux;
    double tf = 0.0, f = 0.0;
    int nvec, *vec;
    if (!(open->bcond_tra[tn] & (TRFLUX|TRCONF))) continue;
    /* Get the index in the tracer dummy array corresponding to tn   */
    for (nn = 0; nn < open->ntflx; nn++)
      if (open->tflx[nn] == tn) break;

    /* Save the contribution to total flux through this OBC          */
    if (wincon->obctr[tn] & TRFLUX) {
      for (ee = 1; ee <= open->no3_e1; ee++) {
	double val = open->t_transfer[open->trm[tn]][ee];
	open->dumtr[nn][ee] = val / open->ncells;
      }
    } else if (wincon->obctr[tn] & TRCONF) {
      /* open->dumtr is set in bdry_transfer_tr(), after which a     */
      /* no-gradient is set so that higher order advection schemes   */
      /* at the first interior face do not use modified boundary     */
      /* ghost values.                                               */
      /* Reset the boundary ghost values to the original prescribed  */
      /* values here.                                                */
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->ogc_t[ee];
	windat->tr_wc[tn][c] = open->dumtr[nn][ee];
	/* Save the concentration so that open->dumtr[nn][cc] can    */
	/* be reset if sub-stepping is performed.                    */
	open->dumn[tn][ee] = open->dumtr[nn][ee];
      }
      /* Use the advection scheme to modify open->dumtr              */
      /* Not yet implemented */

      /* Multipy by the flux for TRCONF (assumes the input value is  */
      /* a concentration rather than a flux).                        */
      flux = windat->u1flux3d;
      vec = open->obc_e1;
      nvec = open->no3_e1;
      for (ee = 1; ee <= nvec; ee++) {
	e = vec[ee];
	sgn = open->dir[ee];
	tf += flux[e];
	if (sgn * flux[e] > 0.0) f += flux[e];
      }
      f = (f) ? fabs(tf / f) : 0.0;
      for (ee = 1; ee <= open->no3_e1; ee++) {
	e = open->obc_e1[ee];
	sgn = open->dir[ee];
	tf = (sgn * flux[c] > 0.0) ? f : 0.0;
	cc = open->e2c_e1[ee];
	/* Get the mass flux scaled for inflow only                */
	open->dumtr[nn][ee] = open->dumn[tn][ee] * (tf * flux[e]);
      }
    }
    /* Reset the fluxes                                            */
    for (ee = 1; ee <= open->no3_e1; ee++) {
      e = open->obc_e1[ee];
      c = open->obc_e2[ee];
      sgn = open->dir[ee];
      dtracer[c] += (sgn * Fx[e] - open->dumtr[nn][ee]) * dt;
    }
    /* Reset the flux for TRCONF if sub-stepping                     */
    if (open->bcond_tra[tn] & TRCONF) {
      for (ee = 1; ee <= open->no3_e1; ee++) {
	open->dumtr[nn][ee] = open->dumn[tn][ee];
      }
    }
  }
}


/* END reset_tr_OBCflux()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set a no-flux condition on lateral boundaries for      */
/* tracers.                                                          */
/*-------------------------------------------------------------------*/
void set_lateral_BC_tr(double **tr, /* Tracer array                  */
                       int ntr,     /* Number of tracers             */
                       int sgbpt,   /* Number of boundary cells      */
                       int *bpt,    /* Boundary cell vector          */
                       int *bin     /* Interior cells to bpt         */
  )
{
  int c1, c2, cc, n, nn;            /* Counters                      */

  /* Set the boundary conditions (no flux)                           */
  for (cc = 1; cc <= sgbpt; cc++) {
    c1 = bpt[cc];
    c2 = bin[cc];
    for (n = 0; n < ntr; n++) {
      tr[n][c1] = tr[n][c2];
    }
  }
}

/* END set_lateral_BC_tr()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Set a no-gradient over OBC ghost cells                            */
/*-------------------------------------------------------------------*/
void set_lateral_bdry_tr(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int e, ee, c, cc, n, nn, tn, co;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->stagger & OUTFACE) {
      for (ee = 1; ee <= open->no3_e1; ee++) {
	c = open->obc_e2[ee];
	co = open->ogc_t[ee];
	do {
	  for (nn = 0; nn < wincon->ntbdy; nn++) {
	    tn = wincon->tbdy[nn];
	    if (!(open->bcond_tra[tn] & (TRCONC|TRCONF))) {
	      windat->tr_wc[tn][co] = windat->tr_wc[tn][c];
	    }
	  }
	  c = open->omape[ee][c];
	} while (c != open->omape[ee][c]);
      }
    }
  }
}

/* END set_lateral_OBC_tr()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set a no-gradient condition on lateral boundaries for  */
/* density. Ordinarily this is not required, since due to the grid   */
/* stagger there are always two wet cells either side of any cell    */
/* face, and hence ghost cells are never used. For nesting this may  */
/* not be the case.                                                  */
/*-------------------------------------------------------------------*/
void Set_lateral_BC_density_w(double *dens, /* Density array         */
                              int sgbpt, /* Number of boundary cells */
                              int *bpt,  /* Boundary cell vector     */
                              int *bin   /* Interior cells to bpt    */
  )
{
  int c1, c2, cc;               /* Counters                          */

  /* Set the boundary conditions (no gradient)                       */
  for (cc = 1; cc <= sgbpt; cc++) {
    c1 = bpt[cc];
    c2 = bin[cc];
    dens[c1] = dens[c2];
  }
}

/* END Set_lateral_BC_density_w()                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set a no-slip condition on lateral boundaries for      */
/* velocity.                                                         */
/*-------------------------------------------------------------------*/
void set_lateral_BC_vel(double *vel,  /* Velocity array              */
                        int sgbpt,    /* Number of boundary cells    */
                        int *bpt,     /* Boundary cell vector        */
                        int *bin      /* Interior cells to bpt       */
  )
{
  int c, cc;                    /* Counters                          */

  /* Set the boundary conditions (no slip)                           */
  for (c = 1; c <= sgbpt; c++) {
    cc = bpt[c];
    vel[cc] = 0.0;
  }
}

/* END set_lateral_BC_vel()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness and the cell coordinate of      */
/* the surface.                                                      */
/*-------------------------------------------------------------------*/
void set_dz(geometry_t *window,     /* Window geometry               */
	    window_t *windat,       /* Window data                   */
	    win_priv_t *wincon      /* Window constants              */
	    )
{
  int cc;                       /* Counter                           */
  int c, c3, c2, cs;            /* Cell coordinate                   */
  int zm1;                      /* Cell cell below cell c            */
  int cb, cbot;                 /* Cell coordinate of the bottom     */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  int vc;                       /* Tracer cells to process vector    */
  int vcs;                      /* Number of tracer cells to process */

  /*-----------------------------------------------------------------*/
  /* Sigma coordinates                                               */
  if (wincon->sigma) {
    int size = window->a2_t + 1;
    /* Get the cells to process vectors                              */
    wincon->vc = window->b3_t;
    wincon->vcs = window->v2_t;
    wincon->vc1 = window->v3_t;
    wincon->vcs1 = window->v2_t;
    memcpy(window->sur_t, window->w2_t, size * sizeof(int));
    memcpy(window->nsur_t, window->w2_t, size * sizeof(int));
    memcpy(wincon->s1, window->w3_t, (window->a3_t + 1) * sizeof(int));
    memcpy(wincon->i1, window->bot_t, size * sizeof(int));
    memcpy(wincon->i2, window->nsur_t, size * sizeof(int));
    memcpy(wincon->i3, window->sur_t, size * sizeof(int));

    for (cc = 1; cc <= window->b2_t; cc++) {
      /* Get the 2D and 3D cell coordinates of the old surface       */
      c = c3 = window->sur_t[cc];
      c2 = window->m2d[c3];
      cbot = window->bot_t[cc];
      top = windat->topz[c2];
      while (c3 != cbot) {
        bot = window->gridz[c3];
        wincon->dz[c3] = top - bot;
        top = bot;
        c3 = window->zm1[c3];
      }

      /* Set the cell thickness at the bottom                        */
      wincon->dz[cbot] = (top - window->botz[c2]);
    }
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates                                                 */
  else {

    /*---------------------------------------------------------------*/
    /* Get the vector containing the cell locations of the new       */
    /* surface (i.e. the surface after the 2D mode is complete).     */
    /* Also reset the old surface vector. This is performed over all */
    /* wet + auxiliary cells.                                        */
    for (cc = 1; cc <= window->a2_t; cc++) {

      /* Get the 2D and 3D cell coordinates                          */
      c = cs = window->w2_t[cc];
      cb = window->bot_t[cc];

      /* Copy the new surface into the old surface vector            */
      window->sur_t[cc] = window->nsur_t[cc];

      /* Get the new cell location of the surface. For elevations    */
      /* below the bottom this is set to the sediment coordinate.    */
      zm1 = window->zm1[c];
      while (c != zm1 && window->gridz[c] >= windat->eta[cs]) {
        c = zm1;
        zm1 = window->zm1[c];
      }
      window->nsur_t[cc] = c;
    }

    if (wincon->thin_merge)
      set_thin_tr(window, windat, wincon);

    /*---------------------------------------------------------------*/
    /* Get the wet cells to process vector for tracers and store in  */
    /* the buffer wincon->s1. This is performed for wet cells and    */
    /* open boundary cells only (open boundary cells are required    */
    /* to calculate vertical velocity, used in momentum advection.   */
    /* However, tracers are also updated on the boundary which uses  */
    /* T at t+1 on the boundary in the UPSTRM OBC).                  */
    /* Cells are arranged thus :                                     */
    /* 1 to wincon->vcs1 = surface layer wet cells                   */
    /* wincon->vcs1+1 to wincon->vcs = surface layer OBC cells       */
    /* wincon->vcs+1 to wincon->vc1 = sub-surface wet cells          */
    /* wincon->vc1+1 to wincon->vc = sub-surface OBC cells           */
    /* First get the cell locations of the surface : this is the     */
    /* lower of the surface elevations before and after the 2D mode. */
    vc = 1;
    for (cc = 1; cc <= window->b2_t; cc++) {
      /* The lower of the layers corresponding to current and        */
      /* previous surface elevations is required here. The cell      */
      /* coordinates in the window system (local coordinates) are    */
      /* always smallest in the surface layer, and increase towards  */
      /* the bottom. Hence the maximum of old and current layer      */
      /* locations will provide the lower level of the two.          */
      /* Store the bottom coordinate for the cells to process (1 to  */
      /* vcs) in the buffer i1. This is required in surf_conc().     */
      c = max(window->sur_t[cc], window->nsur_t[cc]);

      cs = window->m2d[c];
      top = wincon->oldeta[cs];
      top = min(windat->eta[cs], wincon->oldeta[cs]);
      bot = DRY_FRAC * wincon->hmin + window->botz[cs];
      if (top >= bot && c != window->zm1[c]) {
        wincon->s1[vc] = c;
	wincon->i1[vc] = window->bot_t[cc];
        wincon->i2[vc] = window->nsur_t[cc];
        wincon->i3[vc] = window->sur_t[cc];
        if (cc <= window->v2_t)
          wincon->vcs1 = vc;
        vc++;
      }
    }
    vcs = vc - 1;
    wincon->vcs = vcs;

    /* Loop from the layer below the surface to the bottom and get   */
    /* the cells to process.                                         */
    for (cc = 1; cc <= vcs; cc++) {
      c = window->zm1[wincon->s1[cc]];
      zm1 = window->zm1[c];
      while (c != zm1) {
        wincon->s1[vc] = c;
        vc++;
        c = zm1;
        zm1 = window->zm1[c];
      }
      if (cc == wincon->vcs1)
        wincon->vc1 = vc - 1;
    }
    vc--;
    wincon->vc = vc;

    /*---------------------------------------------------------------*/
    /* Set the cell thickness for all cells below this. Note : dz of */
    /* the surface layer is the layer thickness at the start of the  */
    /* time-step; i.e. eta-gridz if calculated before the call to    */
    /* etastep() or oldeta-gridz if calculated after etastep().      */
    for (cc = 1; cc <= window->b2_t; cc++) {
      /* Get the 2D and 3D cell coordinates of the old surface       */
      c = c3 = window->sur_t[cc];
      c2 = window->m2d[c3];
      cbot = window->bot_t[cc];

      top = wincon->oldeta[c2];
      while (c3 < cbot) {
        bot = window->gridz[c3];
        wincon->dz[c3] = top - bot;
        top = bot;
        c3 = window->zm1[c3];
      }

      /* Set the cell thickness at the bottom                        */
      wincon->dz[cbot] = top - window->botz[c2];

      /* Set all cell thickness above the old surface equal to the   */
      /* surface thickness (used in higher order vertical advection  */
      /* schemes).  */
      c3 = c;
      while (c > c2) {
        c = window->zp1[c];
        wincon->dz[c] = wincon->dz[c3];
      }
    }
  }
}

/* END set_dz()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set the velocity over thin layers                      */
/*-------------------------------------------------------------------*/
void set_thin_tr(geometry_t *window,  /* Window geometry             */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon   /* Window constants            */
		 )
{
  int cc;                       /* Counter                           */
  int c;                        /* Cell coordinate                   */
  int cs;                       /* Cell coordinate of the surface    */

  /*-----------------------------------------------------------------*/
  /* Set the thin layer vector and reset the surface layer.          */
  /* Note : if the loop is done to a2_e1 then dzu1=0 at the first    */
  /* timestep for multiple windows in auxiliary cells since it has   */
  /* not yet been copied from the master.                            */
  wincon->nkth_e1 = 1;
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->nsur_t[cc];
    cs = window->m2d[c];
    if (c != window->bot_t[cc] &&
        windat->eta[cs] - window->gridz[c] < wincon->hmin) {
      window->nsur_t[cc] = window->zm1[c];
      /* Save the cell counter of the thin layer                     */
      wincon->kth_e1[wincon->nkth_e1] = (double)cc;
      wincon->nkth_e1++;
    }
  }
  wincon->nkth_e1--;
}

/* END set_thin_tr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness and cells to process to for the */
/* layer conifuration corresponding to the updated elevation         */
/* (windat->eta and window->nsur_t). These dz and cells to process   */
/* vectors should be used in any routines called after tracers are   */
/* updated, e.g. diagnostics and interface routines (ecology,        */
/* sediments).                                                       */
/*-------------------------------------------------------------------*/
void reset_dz(geometry_t *window,  /* Window geometry                */
	      window_t *windat,    /* Window data                    */
	      win_priv_t *wincon   /* Window constants               */
	      )
{
  int cc;                       /* Counter                           */
  int c, cs, cn;                /* Cell coordinate                   */
  int zp1, zm1;                 /* Cell cell above / below cell c    */
  double top;                   /* Top face of a cell                */
  double bot;                   /* Bottom face of a cell             */
  int vc;                       /* Tracer cells to provess counter   */
  int vcs;                      /* Number of tracer cells to process */

  /*-----------------------------------------------------------------*/
  /* Set dz at the surface to account for the updated elevation.     */
  for(cc=1; cc <= window->a2_t; cc++) {
    c = max(window->sur_t[cc], window->nsur_t[cc]);
    cs = window->m2d[c];
    cn = window->nsur_t[cc];
    zp1 = window->zp1[c];
    bot = max(window->gridz[c], window->botz[cs]);
    while(c > cn) {
      top = window->gridz[zp1];
      wincon->dz[c] = top - bot;
      bot = top;
      c = zp1;
      zp1 =  window->zp1[c];
    }
    top = windat->eta[cs];
    wincon->dz[c] = (top > bot) ? top - bot : wincon->dz[c];
  }

  if (!wincon->sigma) {
    /*---------------------------------------------------------------*/
    /* Set the cells to process vector for the updated elevation.    */
    /* Cells are arranged thus :                                     */
    /* 1 to wincon->vca2 = surface layer wet cells                   */
    /* wincon->vca2+1 to wincon->vcs2 = surface layer OBC cells      */
    /* wincon->vcs2+1 to wincon->vci2 = sub-surface wet cells        */
    /* wincon->vci2+1 to wincon->vc2 = sub-surface OBC cells         */
    vc = 1;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->nsur_t[cc];
      cs = window->m2d[c];
      top = windat->eta[cs];
      bot = DRY_FRAC * wincon->hmin + window->botz[cs];
      if (top > bot && c != window->zm1[c]) {
        wincon->s2[vc] = c;
	wincon->i1[vc] = window->bot_t[cc];
        if (cc == window->v2_t)
          wincon->vca2 = vc;
        vc++;
      }
    }
    vcs = vc - 1;
    wincon->vcs2 = vcs;

    /* Loop from the layer below the surface to the bottom and get   */
    /* the cells to process.                                         */
    for (cc = 1; cc <= vcs; cc++) {
      c = window->zm1[wincon->s2[cc]];
      zm1 = window->zm1[c];
      while (c != zm1) {
        wincon->s2[vc] = c;
        vc++;
        c = zm1;
        zm1 = window->zm1[c];
      }
      if (cc == wincon->vca2)
        wincon->vci2 = vc - 1;
    }
    vc--;
    wincon->vc2 = vc;
  }

  /* Set the wet cell mask                                           */
  memset(wincon->c1, 0, window->szc * sizeof(char));
  for (cc = 1; cc <= vc; cc++) {
    c = wincon->s2[cc];
    wincon->c1[c] = 1;
  }
}

/* END reset_dz()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the correct map across open boundaries for tracers */
/*-------------------------------------------------------------------*/
void set_map_t(geometry_t *window /* Window geometry                 */
  )
{
  int n;
  int cc, cs, c;
  int ee, es, e;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & OUTFACE) {
      if (open[n]->bgz) {
	for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	  c = open[n]->obc_e2[ee];
	  open[n]->omape[ee][c] = open[n]->ogc_t[ee];
	}
      } else {
	for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	  c = open[n]->obc_e2[ee];
	  open[n]->omape[ee][c] = c;
	}
      }
    }
  }
}

/* END set_map_t()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to reset the map across tracer open boundaries            */
/*-------------------------------------------------------------------*/
void reset_map_t(geometry_t *window /* Window geometry               */
		 )
{
  int n, c, c1, ee;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodec && open[n]->stagger & OUTFACE) {
      for (ee = 1; ee <= open[n]->no3_e1; ee++) {
	c = open[n]->obc_e2[ee];
	c1 = open[n]->ogc_t[ee];
	open[n]->omape[ee][c] = c1;
      }
    }
  }
}

void reset_map_t_all(geometry_t *window) {
  int n, c, c1, ee;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    for (ee = 1; ee <= open[n]->no3_e1; ee++) {
      c = open[n]->obc_e2[ee];
      c1 = open[n]->ogc_t[ee];
      open[n]->omape[ee][c] = c1;
    }
  }
}

/* END reset_map_t()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the map across multiple ghost cell tracer open   */
/* boundaries for LAGRANGE.                                          */
/*-------------------------------------------------------------------*/
void set_map_l(geometry_t *window   /* Window geometry               */
	       )
{
  int j, n;
  int cc, ci, c, cs;
  open_bdrys_t **open = window->open;

  if (!(window->wincon->trasc & LAGRANGE)) return;

  /*-----------------------------------------------------------------*/
  /* Set the cell mappings in open boundary ghost cells            */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodec) {
      for (cc = 1; cc <= open[n]->no3_t; cc++) {
	c = open[n]->obc_t[cc];
	cs = window->m2d[c];
	j = open[n]->ocodey;
	ci = window->c2c[jp(j, window->npe[cs])][c];
	window->c2c[j][ci] = ci;
	ci = window->c2c[jm(j, window->npe[cs])][c];
	window->c2c[j][ci] = ci;
      }
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
	c = window->zm1[open[n]->bot_t[cc]];
	window->c2c[open[n]->ocodey][c] = c;
      }
    }
  }
}

/* END set_map_l()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Temperature and salinity relaxation. Performed every time-step    */
/* using data saved to tracers rtemp and rsalt.                      */
/*-------------------------------------------------------------------*/
void do_ts_relax(geometry_t *window,  /* Window geometry             */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon   /* Window constants            */
		 )
{
  int cc, c, c3;                       /* Counter                    */
  double rr;
  double sf = 1.0 / 86400.0;
  int tcf, ngf;

  rr = wincon->trinfo_3d[windat->tno].relax_rate;
  tcf = wincon->trinfo_3d[windat->tno].tctype;
  ngf = (wincon->trinfo_3d[windat->tno].flag & RLX_GRD) ? 1 : 0;
  /* Set a no-gradient in the relaxation temperature above msl. Note */
  /* that means in layers above msl can be unreliable if the free    */
  /* surface moves in and out of these layers.                       */
  if (windat->rtemp && ngf) {
    double depth = wincon->trinfo_3d[windat->tno].relax_dum;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      while (window->cellz[c] > depth && c < window->bot_t[cc])
        c = window->zm1[c];
      c3 = c;
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	windat->rtemp[c] = windat->rtemp[c3];
      }
    }
  }
  /* Do the relaxation                                               */
  if (windat->rtemp && rr) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      windat->temp[c] -= windat->dt * rr * (windat->temp[c] - windat->rtemp[c]);
    }
  }
  if (windat->rtemp && windat->temp_tc && tcf & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      rr = (windat->temp_tc[c]) ? sf / windat->temp_tc[c] : 0.0;
      windat->temp[c] -= windat->dt * (windat->temp[c] - windat->rtemp[c]) * rr;
    }
  }

  rr = wincon->trinfo_3d[windat->sno].relax_rate;
  tcf = wincon->trinfo_3d[windat->sno].tctype;
  ngf = (wincon->trinfo_3d[windat->sno].flag & RLX_GRD) ? 1 : 0;
  if (windat->rsalt && ngf) {
    double depth = wincon->trinfo_3d[windat->sno].relax_dum;
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      while (window->cellz[c] > depth && c < window->bot_t[cc])
        c = window->zm1[c];
      c3 = c;
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	windat->rsalt[c] = windat->rsalt[c3];
      }
    }
  }
  if (windat->rsalt && rr) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      windat->sal[c] -= windat->dt * rr * (windat->sal[c] - windat->rsalt[c]);
    }
  }
  if (windat->rsalt && windat->salt_tc && tcf & (RLX_ADPT|RLX_REG|RLX_OBC)) {
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      rr = (windat->salt_tc[c]) ? sf / windat->salt_tc[c] : 0.0;
      windat->sal[c] -= windat->dt * (windat->sal[c] - windat->rsalt[c]) * rr;
    }
  }

  /* Set a no-gradient above the surface                             */
  if (!wincon->sigma) {
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      int cs;
      c = cs = wincon->s2[cc];
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	windat->temp[c] = windat->temp[cs];
	windat->sal[c] = windat->sal[cs];
      }
    }
  }
}

/* END do_ts_relax()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculate the pc-th percentile. If pc < 0 calculate the mean.     */
/*-------------------------------------------------------------------*/
double percs(double *aa, int n, double pc)
{
  int m;
  double dm = 0.0;

  /* Calculate the mean value.                                       */
  if (pc == -1.0) {
    for (m = 0; m < n; m++)
      dm += aa[m];
    return (dm / (double)n);
  }
  /* Use a cyclic value over the subset                              */
  else if (pc == -2) {
    return (aa[n - 1]);
  }
  /* Calculate the percentile value.                                 */
  else {
    sorts(aa, n);
    dm = (double)(n + 1) * pc / 100.0 - 1.0;
    m = (int)floor(dm);
    if (m < 0)
      return (aa[0]);
    else if (m >= n - 1)
      return (aa[n - 1]);
    else
      return (aa[m] + (aa[m + 1] - aa[m]) * (dm - (double)m));
  }
}

/* END percs()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sorts an array into increasing order                              */
/*-------------------------------------------------------------------*/
void sorts(double *a, int n)
{
  int i, j;
  for (i = 0; i < n - 1; ++i)
    for (j = n - 1; i < j; --j)
      orders(&a[j - 1], &a[j]);
}

void orders(double *p, double *q)
{
  double t1;
  if (*p > *q) {
    t1 = *p;
    *p = *q;
    *q = t1;
  }
}

/* END sorts()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the mixed layer depth and store as a         */
/* diagnostic.                                                       */
/*-------------------------------------------------------------------*/
void get_mixed_layer(geometry_t *window,  /* Window geometry         */
		     window_t *windat,    /* Window data             */
		     win_priv_t *wincon   /* Window constants        */
		     )
{
  int cc;                   /* Counter                               */
  int c, zm1;               /* Cell coordinate                       */
  int cs;                   /* Cell coordinate of the surface        */
  int cb;                   /* Cell coordinate of the bottom         */
  int kts;                  /* Index of the top of the pycnocline    */
  int ktb;                  /* Index of the bottom of the pycnocline */
  double top_m;             /* Depth of the top of the pycnocline    */
  double bot_m;             /* Depth of the bottom of the pycnocline */
  double thr = -0.01;       /* Density gradient threshold            */

  if (wincon->mixlayer & DENS_MIX) {
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      cb = wincon->i1[cc];
      cs = wincon->i3[cc];
      zm1 = window->zm1[cs];
      c = window->m2d[cs];
      mld(window, windat, wincon, cs, cb, &top_m, &bot_m, &kts, &ktb);
      windat->mixl[c] = -top_m;
      if (cs == kts && (windat->dens_0[cs] - windat->dens_0[zm1]) /
          (window->cellz[cs] - window->cellz[zm1]) >= thr)
        windat->mixl[c] = -bot_m;
    }
  } else if (wincon->mixlayer & TKE_MIX) {
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      cb = wincon->i1[cc];
      cs = wincon->i3[cc];
      c = window->m2d[cs];
      mldk(window, windat, wincon, cs, cb, &top_m, &bot_m, &kts, &ktb);
      windat->mixl[c] = -top_m;
      if (cs == kts && cb == ktb)
        windat->mixl[c] = -bot_m;
    }
  } else if (wincon->mixlayer & TEMP_MIX) {
    thr = 0.1;
    memset(windat->mixl, 0, window->szcS * sizeof(double));
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      cb = wincon->i1[cc];
      cs = wincon->i3[cc];
      top_m = windat->temp[cs];
      for(c = cs; c <= cb; c = window->zm1[c]) {
	if(fabs(top_m - windat->temp[c]) > thr) {
	  windat->mixl[window->m2d[c]] = window->cellz[c];
	  break;
	}
      }
    }
  }
}

/* END get_mixed_layer()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the advective and vertical diffusive fluxes  */
/* of tracer n and store as a diagnostic.                            */
/*-------------------------------------------------------------------*/
void calc_flux(geometry_t *window,        /* Window geometry         */
	       window_t *windat,          /* Window data             */
	       win_priv_t *wincon,        /* Window constants        */
	       double *fluxe1,            /* Horizontal flux         */
	       double *fluxw,             /* Vertical flux           */
	       double dt,                 /* Time step               */
	       int dir1,                  /* Edge for flux 1         */
	       int dir2                   /* Edge for flux 2         */
	       )
{
  int c, cs, zm1, cc, e;           /* Cell coordinate / counter      */
  int tn = wincon->trflux;         /* Tracer to calculate fluxes for */
  double *tr;                      /* Tracer pointer                 */
  double flux;                     /* Flux value                     */

  /* double *fluxkz = wincon->w10; */

  /*-----------------------------------------------------------------*/
  /* Horizontal advective flux through e1 face                       */
  tr = windat->tr_wc[tn];

  if (windat->fluxe1) {
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	e = window->c2e[dir1][c];
	windat->fluxe1[c] = (windat->fluxe1[c] * windat->meanc[cs] + 
			     window->eSc[dir1][cs] * fluxe1[e] * dt) / 
	  (windat->meanc[cs] + dt);
	e = window->c2e[dir2][c];
	windat->fluxe2[c] = (windat->fluxe2[c] * windat->meanc[cs] + 
			     window->eSc[dir2][cs] * fluxe1[e] * dt) / 
	  (windat->meanc[cs] + dt);
      }
    } else {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	e = window->c2e[dir1][c];
	windat->fluxe1[c] += window->eSc[dir1][cs] * fluxe1[e];
	e = window->c2e[dir2][c];
	windat->fluxe2[c] += window->eSc[dir2][cs] * fluxe1[e];
      }
    }
  }
  
  /* Vertical advective flux through bottom face                     */
  if (windat->fluxw) {
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
	cs = window->m2d[c];
	flux = fluxw[c] * window->cellarea[cs];
	windat->fluxw[c] = (windat->fluxw[c] * windat->meanc[cs] + 
			    flux) / (windat->meanc[cs] + dt);
      }
    } else {
      for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      flux = fluxw[c] * window->cellarea[cs] / dt;
      windat->fluxw[c] += flux;
      }
    }
  }

  /* Vertical diffusive flux through bottom face                     */
  if (windat->fluxkz) {
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
	cs = window->m2d[c];
	zm1 = window->zm1[c];
	flux = -windat->Kz[c] * window->cellarea[cs] * 
	  (tr[c] - tr[zm1]) / wincon->dz[c] * wincon->Ds[cs];
	windat->fluxkz[c] = (windat->fluxkz[c] * windat->meanc[cs] + 
			     flux * dt) / (windat->meanc[cs] + dt);
      }
    } else {
     for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
	cs = window->m2d[c];
	zm1 = window->zm1[c];
	flux = -windat->Kz[c] * window->cellarea[cs] * 
	  (tr[c] - tr[zm1]) / wincon->dz[c] * wincon->Ds[cs];
	windat->fluxkz[c] += flux;
      }
    }
  }
}

/* END calc_flux()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the total mass in a masked region                      */
/*-------------------------------------------------------------------*/
void total_mass(geometry_t *window,      /* Window geometry          */
		window_t *windat,        /* Window data              */
		win_priv_t *wincon       /* Window constants         */
		)
{
  int c, cc, cs;

  wincon->imass = 0.0;
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs = window->m2d[c];
    if (wincon->mask[cs])
      wincon->imass += windat->fltr[c] * window->cellarea[cs] *
        wincon->dz[c] * wincon->Ds[cs];
  }
}

/* END total_mass()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to return the total mass (water column and sediment) at   */
/* a cell.                                                           */
/*-------------------------------------------------------------------*/
double check_mass(geometry_t *window,     /* Window geometry         */
		  window_t *windat,       /* Window data             */
		  win_priv_t *wincon,     /* Window constants        */
		  char *trname,           /* Tracer name             */
		  int cin                 /* Cell for balance        */
		  )
{
  int c, cc, cs, k, cis;
  int tw, ts;
  double d1, d2;

  cis = window->m2d[cin];
  tw = ts = -1;
  for (c = 0; c < windat->ntr; c++)
    if((strcmp(trname, wincon->trinfo_3d[c].name) == 0)) tw = c;
  for (c = 0; c < windat->nsed; c++)
    if((strcmp(trname, wincon->trinfo_sed[c].name) == 0)) ts = c;

  if (tw == -1 || (ts == -1 && !window->sednz))
    hd_quit("check_mass() : Can't find tracer %s\n", trname);

  d1 = d2 = 0.0;
  for(cc = 1; cc <= wincon->vc2; cc++) {
    c = wincon->s2[cc];
    cs=window->m2d[c] ;
    if(cs == cis) 
      d1 += windat->tr_wc[tw][c]*window->cellarea[cs]*wincon->dz[c];
  }
  for(cc = 1; cc <= wincon->vcs2; cc++) {
    c = wincon->s2[cc];
    cs=window->m2d[c];
    for(k=0; k<=window->sednz-1; k++) {
      if(cs == cis)
	d2 += windat->tr_sed[ts][k][cs] * window->cellarea[cs]*
	  (window->gridz_sed[k+1][cs] - window->gridz_sed[k][cs]);
    }
  }
  /*
  printf("mass : %f (%d %d) wc=%e sed=%e\n",windat->t/86400, 
	 window->s2i[cin], window->s2j[cin], d1, d2);
  */
  return(d1 + d2);
}

/* END check_mass()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to retuen the total water column mass within a domain     */
/*-------------------------------------------------------------------*/
double check_tmass(geometry_t *window,      /* Window geometry       */
		   window_t *windat,        /* Window data           */
		   win_priv_t *wincon,      /* Window constants      */
		   char *trname             /* Tracer name           */
		   )
{
  int c, cc, cs;
  int tw;
  double d1;

  tw = -1;
  for (c = 0; c < windat->ntr; c++)
    if((strcmp(trname, wincon->trinfo_3d[c].name) == 0)) tw = c;
  if (tw == -1)
    hd_quit("check_tmass() : Can't find tracer %s\n", trname);

  d1 = 0.0;
  for(cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    cs=window->m2d[c] ;
    d1 += windat->tr_wc[tw][c]*window->cellarea[cs]*wincon->dz[c];
  }
  return(d1);
}

/* END check_tmass()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to read the flushing mask, initialize the flushing tracer */
/* mass and set up the output timeseries file.                       */
/*-------------------------------------------------------------------*/
void init_flushing(master_t *master,        /* Master data structure */
                   geometry_t **window,     /* Window data           */
                   win_priv_t **wincon      /* Window constants      */
  )
{
  FILE *fp = master->prmfd;
  int c, cc, m, n, wn, lc, i, j;
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];

  /*-----------------------------------------------------------------*/
  /* Get the tracer to flush. Note: the concentration of this        */
  /* tracer will be initialized to 1 in the flushing domain and      */
  /* zero elsewhere.  */
  master->flt = NaN;
  master->flts = master->tstart;
  if (master->trflsh) {
    /* Initialize the flushing tracer */
    memset(master->fltr, 0, geom->szc * sizeof(double));
    master->mask = s_alloc_1d(geom->szcS);
    memset(master->mask, 0, geom->szcS * sizeof(short));
    sprintf(keyword, "FLUSHING_PTS");
    if (prm_read_int(fp, keyword, &n)) {
      prm_flush_line(fp);
      if (n > 0) {
	for (m = 0; m < n; ++m) {
	  /* Read i and j indices for this point                     */
	  if (fscanf(fp, "%d %d", &i, &j) != 2)
	    hd_quit("flushing: Can't read i j in points list.\n");
	  prm_flush_line(fp);
	  /* Set the mask and flushing tracer                        */
	  c = geom->map[geom->nz - 1][j][i];
	  master->mask[c] = 1;
	  while (c != geom->zm1[c]) {
	    master->fltr[c] = 1.0;
	    c = geom->zm1[c];
	  }
	}
      }
    } else if (prm_skip_to_end_of_key(fp, "FLUSHING_BLOCKS")) {
      int *iloc, *jloc;
      read_blocks(fp, "FLUSHING_BLOCKS", &n, &iloc, &jloc, NULL);
      for (m = 1; m <= n; m++) {
        c = geom->cc2s[iloc[m]];
        master->mask[c] = 1;
        while (c != geom->zm1[c]) {
          master->fltr[c] = 1.0;
          c = geom->zm1[c];
        }
      }
      i_free_1d(iloc);
    } else if (prm_read_char(fp, "FLUSHING_REGION", buf)) {
      char *files[MAXSTRLEN * MAXNUMARGS];
      int nf, nr, rgn;
      double *regionid;
      nf = parseline(buf, files, MAXNUMARGS);
      regionid = d_alloc_1d(master->geom->szc);
      nr = read_regioni(master, files[0], regionid);
      if (nr) {
	n = 0;
	for (m = 1; m < nf; m++) {
	  sscanf(files[m], "%d", &rgn);
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c]) {
	      master->mask[c] = 1;
	      while (c != geom->zm1[c]) {
		master->fltr[c] = 1.0;
		c = geom->zm1[c];
	      }
	      n++;
	    }
	  }
	}
      }
      d_free_1d(regionid);
    } else
      hd_quit("flushing: No points specified for flushing region.\n");

    /* Set the mask in the windows                                   */
    if (master->nwindows > 1) {
      window_t *windat;
      for (wn = 1; wn <= master->nwindows; wn++) {
	windat = window[wn]->windat;
        wincon[wn]->mask = s_alloc_1d(window[wn]->szcS);
        memset(wincon[wn]->mask, 0, window[wn]->szcS * sizeof(short));
        memset(windat->fltr, 0, window[wn]->szc * sizeof(double));
      }
      for (cc = 1; cc <= geom->b2_t; cc++) {
	c = geom->w2_t[cc];
	lc = geom->fm[c].sc;
	wn = geom->fm[c].wn;
	if (master->mask[c]) {
	  windat = window[wn]->windat;
	  wincon[wn]->mask[lc] = 1;
	  while (lc != window[wn]->zm1[lc]) {
	    windat->fltr[lc] = 1.0;
	    lc = window[wn]->zm1[lc];
	  }
	}
      }
    } else {
      wincon[1]->mask = master->mask;
    }

    /*---------------------------------------------------------------*/
    /* Initialise the mass                                           */
    master->imass = 0.0;
    for (wn = 1; wn <= master->nwindows; wn++) {
      total_mass(window[wn], window[wn]->windat, wincon[wn]);
      master->imass += wincon[wn]->imass;
    }

    /*---------------------------------------------------------------*/
    /* Initialize the output file                                    */
    if (tsflush.fp == NULL) {
      if (strlen(master->opath))
	sprintf(keyword, "%sflushing.ts", master->opath);
      else
	sprintf(keyword, "flushing.ts");
      strcpy(tsflush.pname, keyword);
      tsflush.fp = fopen(keyword, "w");
      tsflush.tsdt = 3600.0;
      sprintf(keyword, "FLUSHING_DT");
      prm_read_int(fp, keyword, &n);
      prm_get_time_in_secs(fp, keyword, &tsflush.tsdt);
      tsflush.master = master;
      tsflush.i = 1;
      tsflush.j = 1;
      tsflush.k = 1;
      tsflush.tsout = master->t;
      fprintf(tsflush.fp, "## COLUMNS 4\n");
      fprintf(tsflush.fp, "##\n");
      fprintf(tsflush.fp, "## COLUMN1.name  Time\n");
      fprintf(tsflush.fp, "## COLUMN1.long_name  Time\n");
      fprintf(tsflush.fp,
              "## COLUMN1.units  %s\n", master->output_tunit);
      fprintf(tsflush.fp, "## COLUMN1.missing_value -999\n");
      fprintf(tsflush.fp, "##\n");
      fprintf(tsflush.fp, "## COLUMN2.name  total_mass\n");
      fprintf(tsflush.fp,
              "## COLUMN2.long_name  Total mass in flushing region\n");
      fprintf(tsflush.fp, "## COLUMN2.units  kg\n");
      fprintf(tsflush.fp, "## COLUMN2.missing_value  0.000000\n");
      fprintf(tsflush.fp, "##\n");
      fprintf(tsflush.fp, "## COLUMN3.name  mass_ratio\n");
      fprintf(tsflush.fp, "## COLUMN3.long_name  Mass ratio\n");
      fprintf(tsflush.fp, "## COLUMN3.units  \n");
      fprintf(tsflush.fp, "## COLUMN3.missing_value  0.000000\n");
      fprintf(tsflush.fp, "##\n");
      fprintf(tsflush.fp, "## COLUMN4.name  flush_time\n");
      fprintf(tsflush.fp, "## COLUMN4.long_name  Flushing time\n");
      fprintf(tsflush.fp, "## COLUMN4.units  days\n");
      fprintf(tsflush.fp, "## COLUMN4.missing_value  0.000000\n");
      fprintf(tsflush.fp, "##\n");
    }
  }
}

/* END init_flushing()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the flushing time based on e-folding time of initially */
/* seeded tracer mass.                                               */
/*-------------------------------------------------------------------*/
void calc_flushing(master_t *master,        /* Master data structure */
		   geometry_t **window,     /* Window geometry       */
		   window_t **windat,       /* Window data           */
		   win_priv_t **wincon      /* Window constants      */
  )
{
  double e = 2.718281828;       /* Constant e                        */
  double t, mt = 0.0, fr = 0.0;
  int wn;

  /* Get the total mass in all windows                               */
  for (wn = 1; wn <= master->nwindows; wn++) {
    mt += wincon[wn]->imass;
  }
  /* If e-folding time is reached print the time                     */
  t = fabs(master->t - master->flts);
  if (mt < master->imass / e && isnan(master->flt))
    master->flt = t / 86400;
  if (master->imass)
    fr = mt / master->imass;
  if (master->t >= tsflush.tsout - DT_EPS) {
    fprintf(tsflush.fp, "%f %f %f %f\n", master->days, mt, fr, master->flt);
    tsflush.tsout += tsflush.tsdt;
    fflush(tsflush.fp);
  }
  /* Reinitialise after 2 * flushing times                           */
  if (!isnan(master->flt) && master->t > master->flts + 1.5 * 86400.0 * master->flt) {
    int c, cc, gc;
    master->flt = NaN;
    master->flts = master->t;
    master->imass = 0.0;
    for (wn = 1; wn <= master->nwindows; wn++) {
      memset(windat[wn]->fltr, 0, window[wn]->szc * sizeof(double));
      for (cc = 1; cc <= window[wn]->b2_t; cc++) {
	c = window[wn]->w2_t[cc];
	gc = window[wn]->wsa[c];
	if (master->mask[gc]) {
	  while (c != window[wn]->zm1[c]) {
	    windat[wn]->fltr[c] = 1.0;
	    master->fltr[gc] = 1.0;
	    c = window[wn]->zm1[c];
	    gc = window[wn]->wsa[c];
	  }
	}
      }
      total_mass(window[wn], window[wn]->windat, wincon[wn]);
      master->imass += wincon[wn]->imass;
    }
  }
}

/* END calc_flushing()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the age mask.                                   */
/*-------------------------------------------------------------------*/
void init_age(master_t *master,        /* Master data structure      */
	      geometry_t **window,     /* Window geometry            */
	      window_t **windat,       /* Window data                */
	      win_priv_t **wincon      /* Window constants           */
	      )
{
  geometry_t *geom = master->geom;
  int c, gc, cc, m, n, i, j;
  FILE *fp = master->prmfd;
  char *files[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  double d1, d2, top, bot;
  short *agemsk;
  int nf;

  n = 0;
  if (strlen(master->trage)) {
    nf = parseline(master->trage, files, MAXNUMARGS);
    /* Initialize the flushing tracer                                */
    agemsk = s_alloc_1d(geom->szc);
    memset(agemsk, 0, geom->szc * sizeof(short));
    /* Get the range                                                 */
    if (prm_read_char(fp, "AGE_RANGE", buf)) {
      if (sscanf(buf, "%lf %lf", &d1, &d2) == 2) {
	top = d1; bot = d2;
	if (bot > top) {
	  d1 = bot;
	  bot = top;
	  top = d1;
	}
      }
    } else {
      top = HUGE; 
      bot = -HUGE;
    }
    if (nf == 1) {
      int *iloc, *jloc;
      read_blocks(fp, "AGE_TR", &n, &iloc, &jloc, NULL);
      for (m = 1; m <= n; m++) {
	c = geom->cc2s[iloc[m]];
	while (c != geom->zm1[c] && geom->cellz[c] <= top && 
	       geom->cellz[c] >= bot) {
	  agemsk[c] = 1;
	  c = geom->zm1[c];
	}
      }
      n = 1;
      i_free_1d(iloc);
    } else {
      int nr;
      double *regionid;
      regionid = d_alloc_1d(master->geom->szc);
      nr = read_regioni(master, files[0], regionid);
      if (nr) {
	int rgn;
	for (m = 1; m < nf; m++) {
	  sscanf(files[m], "%d", &rgn);
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c]) {
	      while (c != geom->zm1[c] && geom->cellz[c] <= top && 
		     geom->cellz[c] >= bot) {
		agemsk[c] = 1;
		c = geom->zm1[c];
	      }
	    }
	  }
	}
	n = 1;
      } else
	hd_quit("age: No points specified for age region.\n");
      d_free_1d(regionid);
    }
    if (n == 1) {
      for (n = 1; n <= master->nwindows; n++) {
	wincon[n]->agemsk = s_alloc_1d(window[n]->szc);
	memset(wincon[n]->agemsk, 0, window[n]->szc);
	for (cc = 1; cc <= window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  gc = window[n]->wsa[c];
	  if (agemsk[gc])
	    wincon[n]->agemsk[c] = 1;
	}
      }
    }
    s_free_1d(agemsk);
  }
}

/* END init_age()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the age tracer increment.                              */
/*-------------------------------------------------------------------*/
void calc_age(geometry_t *window,      /* Window geometry            */
	      window_t *windat,        /* Window data                */
	      win_priv_t *wincon       /* Window constants           */
  )
{
  int c, cc;
  double s = 1.0 / 86400.0;

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    if (wincon->agemsk[c]) {
      windat->agetr[c] += windat->dttr * s;
    }
  }
}

/* END calc_age()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the tracer percentile mask.                     */
/*-------------------------------------------------------------------*/
void init_trperc(master_t *master,     /* Master data structure      */
		 geometry_t **window,  /* Window geometry            */
		 window_t **windat,    /* Window data                */
		 win_priv_t **wincon   /* Window constants           */
		 )
{
  geometry_t *geom = master->geom;
  int c, gc, cc, m, n, i, j;
  FILE *fp = master->prmfd;
  char *files[MAXSTRLEN * MAXNUMARGS], buf[MAXSTRLEN];
  double d1, d2, top, bot;
  short *percmsk;
  int nf, kl;
  int surff = 0;

  n = 0;
  if (strlen(master->trpercr)) {
    nf = parseline(master->trpercr, files, MAXNUMARGS);
    if (!nf) hd_quit("No tracer percentile region specified.\n");
    /* Initialize the flushing tracer                                */
    percmsk = s_alloc_1d(geom->szc);
    memset(percmsk, 0, geom->szc * sizeof(short));
    /* Get the range                                                 */
    if (prm_read_char(fp, "PERC_RANGE", buf)) {
      if (sscanf(buf, "%lf %lf", &d1, &d2) == 2) {
	top = d1; bot = d2;
	if (bot > top) {
	  d1 = bot;
	  bot = top;
	  top = d1;
	}
      } else if (strcmp(buf, "surf") == 0) {
	surff = 1;
      }
    } else {
      top = HUGE; 
      bot = -HUGE;
    }
    if (nf == 1) {
      int *iloc, *jloc;
      read_blocks(fp, "PERC_REGION", &n, &iloc, &jloc, NULL);
      for (m = 1; m <= n; m++) {
	c = geom->cc2s[iloc[m]];
	if (surff)
	  percmsk[c] = 1;
	else {
	  while (c != geom->zm1[c]) {
	    if (geom->cellz[c] <= top && geom->cellz[c] >= bot) percmsk[c] = 1;
	    c = geom->zm1[c];
	  }
	}
      }
      n = 1;
      i_free_1d(iloc);
    } else {
      int nr;
      double *regionid;
      regionid = d_alloc_1d(master->geom->szc);
      nr = read_regioni(master, files[0], regionid);
      if (nr) {
	int rgn;
	for (m = 1; m < nf; m++) {
	  sscanf(files[m], "%d", &rgn);
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c]) {
	      if (surff)
		percmsk[c] = 1;
	      else {
	        while (c != geom->zm1[c]) {
		  if (geom->cellz[c] <= top && geom->cellz[c] >= bot) 
		    percmsk[c] = 1;
		  c = geom->zm1[c];
		}
	      }
	    }
	  }
	}
	n = 1;
      } else
	hd_quit("tr_perc: No points specified for percentile region.\n");
      d_free_1d(regionid);
    }
  }
  if (master->trperc >= 0) {
    if (n == 1) {
      for (n = 1; n <= master->nwindows; n++) {
	wincon[n]->percmsk = s_alloc_1d(window[n]->szc);
	memset(wincon[n]->percmsk, 0, window[n]->szc * sizeof(short));
	for (cc = 1; cc <= window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  gc = window[n]->wsa[c];
	  if (percmsk[gc])
	    wincon[n]->percmsk[c] = 1;
	}
	if (surff)wincon[n]->percmsk[0] = 1;
      }
      if (percmsk) s_free_1d(percmsk);
    } else {
      for (n = 1; n <= master->nwindows; n++) {
	wincon[n]->percmsk = s_alloc_1d(window[n]->szc);
	memset(wincon[n]->percmsk, 0, window[n]->szc * sizeof(short));
	for (cc = 1; cc <= window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  wincon[n]->percmsk[c] = 1;
	}
	if (surff)wincon[n]->percmsk[0] = 1;
      }
    }
  } else
    hd_warn("No tracer percentile region specified.\n");
}

/* END init_trperc()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void mode2d_tracer_init(geometry_t *window,   /* Window geometry     */
			window_t *windat,     /* Window data         */
			win_priv_t *wincon    /* Window constants    */
			)
{
  int c, cs, n, cc;             /* Cell coordinate / counter         */
  double *tr;                   /* Tracer pointer                    */
  double depth;                 /* Water depth                       */
  int size = window->a2_t + 1;

  /* Vertically average all tracers                                  */
  for (n = 0; n < windat->ntr; n++) {
    if (wincon->trinfo_3d[n].advect) {
      memset(wincon->d1, 0, window->szcS * sizeof(double));
      tr = windat->tr_wc[n];
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	wincon->d1[cs] += (tr[c] * wincon->dz[c]);
      }
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	cs = window->m2d[c];
	depth = (windat->eta[cs] - window->botz[cs]) * wincon->Ds[cs];
	tr[c] = wincon->d1[cs] / depth;
      }
    }
  }

  /* Set the cells to process vector                                 */
  wincon->vc = window->b3_t;
  wincon->vcs = window->b2_t;
  wincon->vc1 = window->v3_t;
  wincon->vcs1 = window->v2_t;
  memcpy(window->sur_t, window->w2_t, size * sizeof(int));
  memcpy(window->nsur_t, window->w2_t, size * sizeof(int));
  memcpy(wincon->s1, window->w3_t, (window->a3_t + 1) * sizeof(int));
  memcpy(wincon->i1, window->bot_t, size * sizeof(int));
  memcpy(wincon->i2, window->nsur_t, size * sizeof(int));
  memcpy(wincon->i3, window->sur_t, size * sizeof(int));

  /* Copy the 2d cells to process into 3d arrays                     */
  memcpy(window->w3_t, window->w2_t, (window->n2_t + 1) * sizeof(int));
  window->n3_t = window->n2_t;
  window->a3_t = window->a2_t;
  window->b3_t = window->b2_t;
  window->v3_t = window->v2_t;
  memcpy(window->w3_e1, window->w2_e1, (window->n2_e1 + 1) * sizeof(int));
  window->n3_e1 = window->n2_e1;
  window->a3_e1 = window->a2_e1;
  window->b3_e1 = window->b2_e1;
  window->v3_e1 = window->v2_e1;
}

/* END mode2d_tracer_init()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise diagnostic tracers to zero after dumping to */
/* file.                                                             */
/*-------------------------------------------------------------------*/
void tr_diag_reset_w(geometry_t *window,   /* Window geometry        */
		     window_t *windat,     /* Window data            */
		     win_priv_t *wincon    /* Window constants       */
  )
{
  int n, nn;

  if (windat->df_diagn_set) {
    /* Zero diagnostic tracers                                       */
    for (nn = 0; nn < wincon->ndia; nn++) {
      n = wincon->diagn[nn];
      memset(windat->tr_wc[n], 0, window->szc * sizeof(double));
    }
    windat->df_diagn_set = 0;
  }
}

/* END tr_diag_reset_w()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise diagnostic tracers on the master to zero    */
/* after dumping to file.                                            */
/*-------------------------------------------------------------------*/
void tr_diag_reset_m(master_t *master    /* Master data structure    */
  )
{

  if (master->df_diagn_set) {
  }
  master->df_diagn_set = 0;

  /* The following code performs resets prior to hd_step. This is    */
  /* required if a tracer reset is used to create a tracer           */
  /* difference using tracerstats.                                   */
  /* No longer required : tracer resets are syncrosized : 
     see tr_reset_event().
  for (tt = 0; tt < master->nres; tt++) {
    tn = master->reset[tt];
    strcpy(schedName, "tracer_reset:");
    strcat(schedName, master->trinfo_3d[tn].name);
    tr_reset_event(sched_get_even_by_name(schedule, schedName), master->t);
  }
  for (tt = 0; tt < master->nres2d; tt++) {
    tn = master->reset2d[tt];
    strcpy(schedName, "tracer_reset2d:");
    strcat(schedName, master->trinfo_2d[tn].name);
    tr_reset2d_event(sched_get_even_by_name(schedule, schedName), master->t);
  }
  */
}

/* END tr_diag_reset_m()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* The diff(:) tracerstat fuction can be used to find the difference */
/* between a simulated variable and a variable from a previous run   */
/* read in via the tracer reset function. However, reset tracers are */
/* interpolated onto the grid using a linear or inverse weighed      */
/* scheme, which does not result in exactly the same values as the   */
/* original distribution, even on the same grid (i.e. the difference */
/* between two identical tracers is non-zero if one is read in via   */
/* the reset function). This routine will reset values to zero if    */
/* are below a threshold, so that differences reflect real values    */
/* rather than interpolation effects.                                */
/* No longer required : using cell formatted files for reset does    */
/* not perform interpolation.                                        */
/*-------------------------------------------------------------------*/
void diff_mask(geometry_t *window,   /* Window geometry              */
	       window_t *windat,     /* Window data                  */
	       win_priv_t *wincon    /* Window constants             */
	       )
{
  double thr = 0.025;       /* Difference threshold                  */
  int c, cc, tn, tm = windat->ntr;
  char buf[MAXSTRLEN]; 
  char tr1[MAXSTRLEN], tr2[MAXSTRLEN];

  for (tn = 0; tn < windat->ntr; tn++) {
    sprintf(buf, "%s", wincon->trinfo_3d[tn].tracerstat);
    if (strlen(buf)) {
      sprintf(buf,"%s",strtok(buf, "("));
      /* Check a tracer for tracerstat = 'diff'                      */
      if (strcmp(buf, "diff") == 0) {
	sprintf(tr1,"%s",strtok(NULL, ":"));
	sprintf(tr2,"%s",strtok(NULL, ")"));
	/* Check for one of the difference tracers being reset       */
	for (tm = 0; tm < windat->ntr; tm++) {
	  if ((strcmp(tr1, wincon->trinfo_3d[tm].name) == 0) &&
	      strlen(wincon->trinfo_3d[tm].reset_file)) break;
	  if ((strcmp(tr2, wincon->trinfo_3d[tm].name) == 0) &&
	      strlen(wincon->trinfo_3d[tm].reset_file)) break;
	}
	break;
      }
    }
  }
  if (tm == windat->ntr)
    return;
  /* If a tracer has a difference tracerstat and one of the          */
  /* difference tracers is being reset, then limit the difference to */
  /* the threshhold.                                                 */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    if (fabs(windat->tr_wc[tn][c]) < thr)
      windat->tr_wc[tn][c] = 0.0;
  }
}

/* END diff_mask()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the tracer increment to the relevant state         */
/* variable.                                                         */
/*-------------------------------------------------------------------*/
void do_ts_increment(geometry_t *window,   /* Window geometry        */
		     window_t *windat,     /* Window data            */
		     win_priv_t *wincon    /* Window constants       */
		     )
{
  int n, c, cc;

  for (n = 0; n < windat->ntr; n++) {
    if (wincon->trinfo_3d[n].increment & TEMP && windat->trinc[n]) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	windat->temp[c] += windat->tr_wc[n][c];
      }
      /* To update the tracer on the reset_input_dt time, then set   */
      /* windat->trinc = 0 below. To update every time-step, then    */
      /* set windat->trinc = 1.                                      */
      windat->trinc[n] = 1;
    }
    if (wincon->trinfo_3d[n].increment & SALT && windat->trinc[n]) {
      for (cc = 1; cc <= window->b3_t; cc++) {
	c = window->w3_t[cc];
	windat->sal[c] += windat->tr_wc[n][c];
      }
      windat->trinc[n] = 1;
    }
  }
}

/* END do_increment()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implicit vertical diffusion scheme for tracer var at surface      */
/* location cc.                                                      */
/*-------------------------------------------------------------------*/
void implicit_vdiff_at_cc(geometry_t *window, /* Window geometry     */
			  window_t *windat,   /* Window data         */
			  win_priv_t *wincon, /* Window constants    */
			  double *var,        /* Variable to diffuse */
			  double *Kz,      /* Mixing coefficients    */
			  double *dzcell,  /* Cell thicknesses       */
			  double *dzface,  /* Thicknesses at face    */
			  double *fb,      /* Flux out of bottom     */
			  double *ft,      /* Flux out of top        */
			  int *ctp,        /* Cells to process       */
			  int *cbt,        /* Bottom coordinates     */
			  int cc,          /* Index of surface cells */
			  double *Splus,   /* + part of source term  */
			  double *Sminus,  /* - part of source term  */
			  double *scale,   /* SIGMA : scaling        */
			  double *C, double *Cp1, double *Cm1,
			  double dt
			  )
{
  int c, k;                     /* Cell/layer coordinate             */
  int cs, ks;                   /* Surface coordinate                */
  int cb, kb;                   /* Bottom coordinate                 */
  int zm1;                      /* Cell cell below c                 */
  int c2;                       /* 2D cell corresponding to 3D cell  */
  double dzdt;                  /* dz / dt                           */
  double div;                   /* Constant                          */

  /*-----------------------------------------------------------------*/
  /* Set pointers.                                                   */
  /* Note: the 3D work arrays wincon->w# could be used for the       */
  /* dummy arrays below, but execution speed is considerably faster  */
  /* when work array access is sequential in memory, hence the       */
  /* mapping to a contiguous vertical 1D work array wincon->v#.      */
  double *rhs = wincon->v1;
  double *sol = wincon->v2;
  double *ud = wincon->v3;
  double *B = wincon->v4;
  double *Bm1 = wincon->v5;
  double *Bp1 = wincon->v6;
  double *dz = wincon->v7;

  cs = c = ctp[cc];      /* Set cs to the surface cell coordinate  */
  cb = cbt[cc];          /* Bottom cell coordinate                 */
  c2 = window->m2d[c];   /* 2D cell location corresponding to c    */
  zm1 = window->zm1[c];
  ks = window->s2k[cs];
  kb = window->s2k[cb];

  /*-----------------------------------------------------------------*/
  /* Single layer case (i.e. surface lies in the bottom layer)       */
  if (zm1 == window->zm1[zm1] || !dzcell[zm1]) {
    var[ks] +=
      windat->dt * (fb[c2] - ft[c2]) / max(dzcell[cs], wincon->hmin);
    if (c != wincon->i4[cc]) {
      k = window->s2k[wincon->i4[cc]];
      var[k] = var[ks];
    }
    return;
  }

   /* Map the cell arrays to local coordinates.                     */
  c = cs;
  for (k = ks; k >= kb; k--) {
    B[k] = C[c];
    Bm1[k] = Cm1[c];
    Bp1[k] = Cp1[c];
    dz[k] = dzcell[c];
    c = window->zm1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Set up the rhs for the system of equations.                     */
  /* Bottom layer.                                                   */
  c = cb;                     /* Set c to the bottom cell coordinate */
  dzdt = dz[kb] / dt;
  rhs[kb] = dzdt * var[kb] + fb[c2];
  
  /* Mid-water layers                                                */
  c = window->zp1[cb];        /* Layer above the bottom              */
  for (k = kb + 1; k < ks; k++) {
    dzdt = dz[k] / dt;
    rhs[k] = dzdt * var[k];
    c = window->zp1[c];
  }

  /* Surface layer                                                   */
  dzdt = dz[ks] / dt;
  rhs[ks] = dzdt * var[ks] - ft[c2];

  /*-----------------------------------------------------------------*/
  /* Add positive part of source terms, if specified                 */
  if (Splus) {
    zm1 = window->zm1[c];
    c = cb;
    for (k = kb; k <= ks; k++) {
      rhs[k] += dz[k] * Splus[c];
      c = window->zp1[c];
    }
  }

  /* Add negative part of source terms, if specified                 */
  if (Sminus) {
    c = cb;
    for (k = kb; k <= ks; k++) {
      B[k] += dz[k] * Sminus[c] / var[k];
      c = window->zp1[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Solve tridiagonal system                                        */
  div = B[kb];
  sol[kb] = rhs[kb] / div;
  for (k = kb + 1; k <= ks; k++) {
    ud[k] = Bp1[k - 1] / div;
    div = B[k] - Bm1[k] * ud[k];
    if (div == 0.0)
      hd_quit_and_dump("Tracer diffusion;implicit_vdiff_at_cc: zero divisor\n");
    sol[k] = (rhs[k] - Bm1[k] * sol[k - 1]) / div;
  }

  /*-----------------------------------------------------------------*/
  /* Update the variable                                             */
  c = cs;
  var[ks] += (sol[ks] - var[ks]) * scale[c2];
  for (k = ks - 1; k >= kb; k--) {
    c = window->zm1[c];
    sol[k] -= ud[k + 1] * sol[k + 1];
    var[k] += (sol[k] - var[k]) * scale[c2];
  }

  /*-----------------------------------------------------------------*/
  /* Set the concentration in thin layers                            */
  c = wincon->i4[cc];
  k = window->s2k[c];
  if (c != cs) {
    var[k] += (sol[ks] - var[ks]) * scale[c2];
  }
}

/* END implicit_vdiff_at_cc()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Computes the number of wet cells for each boundary                */
/*-------------------------------------------------------------------*/
void get_bdry_cellno(geometry_t *window, /* Processing window        */
		     window_t *windat,   /* Window data structure    */
		     win_priv_t *wincon  /* Window constants         */
  )
{
  int n, nn, c, cc, c2, zm1;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn]; 
    /* Currently only used with TRFLUX OBCs                          */
    if (open->ntflx) {
      open->ncells = 0.0;
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	c2 = window->m2d[c];
	zm1 = window->zm1[c];
	while (c != zm1) {
	  if (window->gridz[c] <= windat->eta[c2]) open->ncells += 1.0;
	  c = zm1;
	  zm1 = window->zm1[c];
	}
      }
    }
  }
}

/* END get_bdry_cellno()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Re-sets the river flow profile to the parabolic profile           */
/*-------------------------------------------------------------------*/
void reset_flow(geometry_t *window,    /* Processing window          */
		window_t *windat,      /* Window data                */
		win_priv_t *wincon,    /* Window constants           */
		int mode
		)
{
  int ee, e, nn;
  double *vel, d1;
  int *vec;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn]; 
    if (open->bcond_nor & CUSTOM && open->options & OP_PARFLOW && 
	open->options & (OP_GEOSTR|OP_YANKOVSKY)) {
      vel = windat->u1;
      vec = open->obc_e1;
      if (mode == 0) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = vec[ee];
	  d1 = vel[e];
	  vel[e] = open->flow[ee];
	  open->flow[ee] = d1;
	}
      }
      if (mode == 1) {
	for (ee = 1; ee <= open->no3_e1; ee++) {
	  e = vec[ee];
	  vel[e] = open->flow[ee];
	}
      }
    }
  }
}

/* END reset_flow()                                                  */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Initialises regions for swr estimation                            */
/*-------------------------------------------------------------------*/
void swr_params_init(master_t *master, geometry_t **window)
{
  parameters_t *params = master->params;
  geometry_t *geom = master->geom;
  win_priv_t *wincon;
  window_t *windat;
  int wn, nr, reg;
  int cc, c, lc, n, found;
  double *regionid, depth = 0.0;

  /*-----------------------------------------------------------------*/
  /* Initialise for no swr parameter estimation                      */
  if (!strlen(params->swr_regions)) {
    for (wn = 1; wn <= master->nwindows; wn++) {
      wincon = window[wn]->wincon;
      windat = window[wn]->windat;
      wincon->swr_next = master->t - 1;
      wincon->swr_dt = 0.0;
      wincon->nswreg = 0;
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Get the data to converge towards (GHRSST or a tracer)           */
  lc = -1;
  if (strcmp(params->swr_data, "GHRSST") == 0) {
    if (!params->ghrsst)
      hd_quit("SWR estimation must have GHRSST invoked.\n");
  } else {
    char *fields[MAXSTRLEN * MAXNUMARGS];
    n = parseline(params->swr_data, fields, MAXNUMARGS);
    if (n == 1) {           /* 2D data                               */
      found = 0;
      for (n = 0; n < master->ntrS; n++) {
	if (strcmp(fields[0], master->trinfo_2d[n].name) == 0) {
	  lc = n;
	  found = 1;
	  break;
	}
      }
      if (!found)
	hd_quit("SWR estimation must have 2D tracer %s in the tracer list.\n", fields[0]);
    }
    if (n == 2) {           /* 3D data                               */
      found = 0;
      for (n = 0; n < master->ntr; n++) {
	if (strcmp(fields[0], master->trinfo_3d[n].name) == 0) {
	  lc = n;
	  found = 1;
	  break;
	}
      }
      if (!found)
	hd_quit("SWR estimation must have 3D tracer %s in the tracer list.\n", fields[0]);
    }
    depth = atof(fields[1]);
  }
  for (wn = 1; wn <= master->nwindows; wn++) {
    wincon = window[wn]->wincon;
    wincon->swr_data = lc;
    wincon->swr_depth = -fabs(depth);
  }

  /*-----------------------------------------------------------------*/
  /* Estimate the swr parameters at every column in the grid         */
  if (strcmp(params->swr_regions, "ALL") == 0) {
    for (wn = 1; wn <= master->nwindows; wn++) {
      wincon = window[wn]->wincon;
      windat = window[wn]->windat;
      wincon->nswreg = window[wn]->b2_t;
      wincon->swr_dt = params->swreg_dt;
      wincon->swr_next = master->t;
      windat->swrc = 0.0;
      wincon->swmap = i_alloc_1d(window[wn]->sgsizS);
      wincon->swC = d_alloc_1d(window[wn]->sgsiz);
      for (cc = 1; cc <= window[wn]->b2_t; cc++) {
	c = window[wn]->w2_t[cc];
	lc = window[wn]->wsa[c];
	windat->swreg[c] = cc;
	wincon->swmap[c] = cc-1;
	if (windat->swr_attn[c] <= 0.0) {
	  windat->attn_mean[c] = windat->swr_attn[c] = 0.0;
	  master->attn_mean[lc] = master->swr_attn[lc] = 0.0;
	} else {
	  windat->attn_mean[c] = -windat->swr_attn[c];
	  master->attn_mean[lc] = -master->swr_attn[lc];
	}
	if (windat->swr_tran[c] <= 0.0) {
	  windat->tran_mean[c] = windat->swr_tran[c] = 0.0;
	  master->tran_mean[c] = master->swr_tran[c] = 0.0;
	} else {
	  windat->tran_mean[c] = -windat->swr_tran[c];
	  master->tran_mean[c] = master->swr_tran[c] = 0.0;
	}
      }
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Estimate the swr parameters at the geographic centre of the     */
  /* grid.                                                           */
  if (strcmp(params->swr_regions, "ONE") == 0) {
    for (wn = 1; wn <= master->nwindows; wn++) {
      wincon = window[wn]->wincon;
      windat = window[wn]->windat;
      wincon->nswreg = 1;
      wincon->swr_dt = params->swreg_dt;
      wincon->swr_next = master->t;
      windat->swrc = 0.0;
      wincon->swmap = i_alloc_1d(window[wn]->sgsizS);
      wincon->swC = d_alloc_1d(window[wn]->sgsiz);
      for (cc = 1; cc <= window[wn]->b2_t; cc++) {
	c = window[wn]->w2_t[cc];
	windat->swreg[c] = 0;
	wincon->swmap[c] = 0;
      }
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Estimate the swr parameters according to regional partitioning  */
  regionid = d_alloc_1d(master->geom->sgsiz);
  nr = read_regioni(master, params->swr_regions, regionid);

  if (nr) {
    /* Set the mask in the windows */
    for (wn = 1; wn <= master->nwindows; wn++) {
      wincon = window[wn]->wincon;
      memset(wincon->d1, 0, window[wn]->sgsizS * sizeof(int));
      wincon->nswreg = 0;
      wincon->swr_dt = params->swreg_dt;
      wincon->swr_next = master->t;
      windat->swrc = 0.0;
      wincon->swmap = i_alloc_1d(window[wn]->sgsizS);
      wincon->swC = d_alloc_1d(window[wn]->sgsiz);
    }
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      lc = geom->fm[c].sc;
      wn = geom->fm[c].wn;
      reg = (int)regionid[c];
      if (reg >= 0) {
	wincon = window[wn]->wincon;
	windat = window[wn]->windat;
	windat->swreg[lc] = (double)reg;
	found = 0;
	for (n = 0; n < wincon->nswreg; n++) {
	  if (reg == wincon->d1[n]) {
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  wincon->d1[wincon->nswreg] = reg;
	  wincon->nswreg++;
	}
      }
      for (n = 0; n < wincon->nswreg; n++)
	if (reg == wincon->d1[n])
	  wincon->swmap[lc] = n;
    }
  }
  d_free_1d(regionid);
}

/* END swr_params_init()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to optimize swr_attn and swr_tran from an ensemble to     */
/*  minimise SST error compared GHRSST.                              */
/*-------------------------------------------------------------------*/
double swr_params_event(geometry_t *window,
			window_t *windat,
			win_priv_t *wincon,
			int n
			)
{
  int tn, k, c, cc, cs, cb, c2, cd, ks, kb, m;
  int i1, i2, zm1, zp1;
  double *tr = windat->tr_wc[n];    /* Temperature tracer            */
  double *ghrsst;                   /* GHRSST tracer                 */
  double *topflux = windat->heatf;  /* Surface flux                  */
  double *Splus = wincon->w8;   /* swr source                        */
  int *ctp = wincon->i2;        /* Old surface sparse coordinate     */
  int *cbt = wincon->i5;        /* Bottom sparse coordinate          */
  double *dzface = wincon->w7;  /* Cell thickness at the cell face   */
  double *dzcell = wincon->w9;  /* Cell centered cell thickness      */
  double *Kz = wincon->w10;     /* Vertical diffusivity              */
  double *Cm1 = wincon->w1;     /* Constant for implicit calculation */
  double *C = wincon->w2;       /* Constant for implicit calculation */
  double *Cp1 = wincon->w3;     /* Constant for implicit calculation */
  double *temp = wincon->v1;    /* Temperature vector                */
  double attn;                  /* swr attenuation                   */
  double tran;                  /* swr transmission                  */
  double babs;                  /* swr bottom absorption             */
  double rms;                   /* Global rms error                  */
  double swr;                   /* SWR at the estimation location    */
  double *sstm;                 /* Regional mean of GHRSST           */
  double *tm;                   /* Surface temp closest to sstm      */
  double *nreg;                 /* Number of cells in each region    */
  int *ccs;                     /* Cell cc of clostest temp to sstm  */
  double *mattn, *mtran;        /* Region optimised atten and tran   */
  double dzdt;                  /* dz / dt for tridiagnol            */
  double *heatf = wincon->d2;   /* Copy of surface heatflux          */
  double *data;                 /* Target temperature array          */ 
  double tempt;                 /* Target temperature value          */
  double tmin = 0.0;            /* Minimum allowable temperature     */
  double tmax = 35.0;           /* Maximum allowable temperature     */
  double attn0 = 0.02;          /* Start attenuation for ensemble    */
  double attni = 0.05;          /* Atten increment for ensemble      */
  double attns = 10.0;          /* Attn scaling for ensemble         */
  double tran0 = 0.6;
  double trani = 0.1;
  double trans = 10.0;
  int i1s, i1e, i1i;
  int i2s, i2e, i2i;

  /* Return if next swr event isn't scheduled                        */
  if (windat->t < wincon->swr_next) return(wincon->swr_next);
  wincon->swr_next = windat->t + wincon->swr_dt;

  /*-----------------------------------------------------------------*/
  /* Allocate and ititialize                                         */
  tm = d_alloc_1d(wincon->nswreg);
  sstm = d_alloc_1d(wincon->nswreg);
  nreg = d_alloc_1d(wincon->nswreg);
  ccs = i_alloc_1d(wincon->nswreg);
  mattn = d_alloc_1d(wincon->nswreg);
  mtran = d_alloc_1d(wincon->nswreg);
  data = (wincon->swr_data >= 0) ? windat->tr_wcS[wincon->swr_data] : windat->ghrsst;
  memset(sstm, 0, wincon->nswreg * sizeof(double));
  memset(nreg, 0, wincon->nswreg * sizeof(double));
  for (m = 0; m < wincon->nswreg; m++) {
    tm[m] = HUGE;
    ccs[m] = 0;
  }

  /*-----------------------------------------------------------------*/
  /* Get the mean SST in each region                                 */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = ctp[cc];
    c2 = window->m2d[c];
    m = wincon->swmap[c2];
    if (wincon->swr_depth != 0.0) {
      while (window->gridz[c] > wincon->swr_depth && c != window->zm1[c]) {
	c = window->zm1[c];
      }
      c2 = c;
    }
    if (data[c2] >= tmin && data[c2] <= tmax) {
      sstm[m] += data[c2];
      nreg[m] += 1.0;
    }
  }
  for (m = 0; m < wincon->nswreg; m++) {
    if (nreg[m]) {
      sstm[m] /= nreg[m];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Find the cell in each region with a sst closest to the mean     */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = ctp[cc];
    c2 = window->m2d[c];
    m = wincon->swmap[c2];
    if (wincon->swr_depth != 0.0) {
      while (window->gridz[c] > wincon->swr_depth && c != window->zm1[c]) {
	c = window->zm1[c];
      }
    }
    rms = sqrt((tr[c] - sstm[m]) * (tr[c] - sstm[m]));
    if (rms < tm[m]) {
      ccs[m] = cc;
      tm[m] = rms;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set up the tridiagnol timestep dependent term                   */
  for (m = 0; m < wincon->nswreg; m++) {
    cc = ccs[m];
    cs = c = ctp[cc];
    cb = cbt[cc];
    c2 = window->m2d[cs];
    ks = window->s2k[cs];
    kb = window->s2k[cb];
    zm1 = window->zm1[c];

    /* Single layer case (i.e. surface lies in the bottom layer)     */
    if (zm1 == window->zm1[zm1] || !dzcell[zm1])
      continue;

    c = cb;
    zp1 = window->zp1[c];       /* Layer above the bottom        */

    /* Set up tri-diagonal set of equations.                         */
    /* Bottom layer.                                                 */
    dzdt = dzcell[c] / wincon->swr_dt;
    wincon->swC[c] = dzdt - Cp1[c];
    /* Mid-water layers                                              */
    while (zp1 != cs) {
      c = zp1;
      zp1 = window->zp1[c];
      dzdt = dzcell[c] / wincon->swr_dt;
      wincon->swC[c] = dzdt - Cm1[c] - Cp1[c];
    }
    /* Surface layer                                                 */
    c = cs;
    dzdt = dzcell[c] / wincon->swr_dt;
    wincon->swC[c] = dzdt - Cm1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Repopulate the mean parameters with input values for            */
  /* regionalised estimation with resets for the fixed regional      */
  /* values.                                                         */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    if (wincon->swr_type & SWR_ATTN && windat->attn_mean[c] < 0.0)
      windat->attn_mean[c] = -windat->tr_wcS[wincon->attn_tr][c];
    if (wincon->swr_type & SWR_TRAN && windat->tran_mean[c] < 0.0) 
      windat->tran_mean[c] = -windat->tr_wcS[wincon->tran_tr][c];
  }

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion over the ensemble                     */
  i1s = i2s = 0;
  i1e = i2e = 10;
  i1i = i2i = 1;
  /*
  if (wincon->swr_type & SWR_ATTN) {
    i2s = i2e = 0;
    i1e = 100;
    i1i = 1;
    trani = 0.01;
    trans = 1.0;
  }
  if (wincon->swr_type & SWR_TRAN) {
    i1s = i1e = 0;
    i2e = 100;
    i2i = 1;
    attni = 0.005;
    attns = 1.0;
  }
  */
  memset(windat->swrms, 0, window->sgsizS * sizeof(double));
  memset(Splus, 0, window->sgsiz * sizeof(double));
  memcpy(heatf, windat->heatf, window->sgsizS * sizeof(double));
  for (m = 0; m < wincon->nswreg; m++) {
    cc = ccs[m];
    if (!cc) continue;
    cs = c = cd = ctp[cc];
    cb = cbt[cc];
    c2 = window->m2d[cs];
    ks = window->s2k[cs];
    kb = window->s2k[cb];
    swr = windat->swr[c2];
    if (wincon->swr_type & SWR_ATTN && windat->attn_mean[c2] < 0.0) 
      mattn[m] = fabs(windat->attn_mean[c2]);
    if (wincon->swr_type & SWR_TRAN && windat->tran_mean[c2] < 0.0) 
      mtran[m] = fabs(windat->tran_mean[c2]);
    if (swr == 0.0) continue;
    windat->swrms[c2] = HUGE;

    /* Get the coordinate of the data                                */
    if (wincon->swr_depth != 0.0) {
      while (window->gridz[cd] > wincon->swr_depth && cd != window->zm1[cd]) {
	cd = window->zm1[cd];
      }
    } else
      cd = c2;

    tempt = data[cd];
    if (tempt < tmin || tempt > tmax) tempt = sstm[m];

    for (i1 = i1s; i1 <= i1e; i1 += i1i) {   /* Transmission range        */
      for (i2 = i2s; i2 <= i2e; i2 += i2i) { /* Attenuation range         */
	/* Set swr parameters                                        */
	attn = attn0 + (double)i2 * attni / attns;
	tran = tran0 + (double)i1 * trani / trans;
	if (tran > 1.0) tran -= 1.0;
	babs = windat->swr_babs[c2];
	/* Note: for regionalized estimation the negative values are */
	/* not subject to estimation. These are set in               */
	/* swr_params_init().                                        */
	if (wincon->swr_type & SWR_ATTN && windat->attn_mean[c2] < 0.0) 
	  attn = fabs(windat->attn_mean[c2]);
	if (wincon->swr_type & SWR_TRAN && windat->tran_mean[c2] < 0.0)
	  tran = fabs(windat->tran_mean[c2]);
	windat->swr_attn[c2] = attn;
	windat->swr_tran[c2] = tran;

	/* Set the temperature profile                               */
	c = cs;
	for (k = ks; k >= kb; k--) {
	  temp[k] = tr[c];
	  c = window->zm1[c];
	}

	/* Set the swr distribution                                  */
	calc_swr(window, windat, wincon, Splus, cc);

	/* Mix vertically                                            */
	implicit_vdiff_at_cc(window, windat, wincon, temp, Kz, dzcell, dzface,
			     wincon->d1, windat->heatf, ctp, cbt, cc, Splus,
			     NULL, wincon->one, wincon->swC, Cp1, Cm1, wincon->swr_dt);
	windat->heatf[c2] = heatf[c2];

	/* Save the swr parameters if the SST is improved            */
	rms = sqrt((data[cd] - temp[ks]) * (data[cd] - temp[ks]));
	if (rms < windat->swrms[c2]) {
	  mattn[m] = attn;
	  mtran[m] = tran;
	  windat->swrms[c2] = rms;
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Repopulate the swr parameters with optimized values             */
  k = 0;
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    if (windat->swr[c]) {
      k = 1;
      m = wincon->swmap[c];
      if (windat->attn_mean[c] >= 0)
	windat->attn_mean[c] = (windat->attn_mean[c] * windat->swrc + 
				mattn[m] * wincon->swr_dt) / 
	  (windat->swrc + wincon->swr_dt);
      if (windat->tran_mean[c] >= 0)
	windat->tran_mean[c] = (windat->tran_mean[c] * windat->swrc + 
				mtran[m] * wincon->swr_dt) / 
	  (windat->swrc + wincon->swr_dt);
      /*
      if (!(wincon->swr_type & SWR_ATTN)) windat->swr_attn[c] = windat->attn_mean[c];
      if (!(wincon->swr_type & SWR_TRAN)) windat->swr_tran[c] = windat->tran_mean[c];
      */
      windat->swr_attn[c] = fabs(windat->attn_mean[c]);
      windat->swr_tran[c] = fabs(windat->tran_mean[c]);
    }
  }
  if (k) windat->swrc += wincon->swr_dt;
  d_free_1d(tm);
  d_free_1d(sstm);
  d_free_1d(nreg);
  i_free_1d(ccs);
  d_free_1d(mattn);
  d_free_1d(mtran);
}

/* END swr_params_event()                                            */
/*-------------------------------------------------------------------*/
