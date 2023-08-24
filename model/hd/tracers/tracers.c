/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/tracers/tracers.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracers.c 6844 2021-06-30 00:15:15Z her127 $
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
/* quickest_uniform() : Calculates QUICKEST fluxes, uniform grid     */
/* van_leer()      : Calculates Van Leer higher order upwind fluxes  */
/* ff_sl()         : Calculates flux-form semi-lagrangian fluxes     */
/* surf_conc()     : Calculates surface tracer concentration         */
/* ulimate_filter() : Invokes the ULTIMATE limiter                   */
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

void set_thin_tr(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_mixed_layer(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void calc_age(geometry_t *window, window_t *windat, win_priv_t *wincon);
void auxiliary_routines(geometry_t *window, window_t *windat,
			win_priv_t *wincon);
void reset_dz(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect_diffuse_split(geometry_t *window, window_t *windat,
			  win_priv_t *wincon);
double check_mass(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, char *trname, int cin);
double check_tmass(geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, char *trname);
void check_fluxes(int n, char *trname, char * text, double dtu, 
		  geometry_t *window,
		  window_t *windat, win_priv_t *wincon);
void limit_min(geometry_t *window, window_t *windat, win_priv_t *wincon);
void do_ts_relax(geometry_t *window, window_t *windat, win_priv_t *wincon);
void do_ts_increment(geometry_t *window, window_t *windat, 
		     win_priv_t *wincon);
void reset_tr_OBCflux(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      double *dtracer, double *Fx, double *Fy, double dt, int tn,
		      int mode);
void swr_assimilation(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      int n);
void save_OBC_tr(geometry_t *window, window_t *windat, win_priv_t *wincon,
		 double *Fx, double *Fy, double *tr, int n, int mode);
void set_lateral_bdry_tr(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_bdry_cellno(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_flow(geometry_t *window, window_t *windat, win_priv_t *wincon, int mode);
double monoco(int c, int cm, int cp, double *tr, double *trmod, double sc);
void nogr(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr);

/*-------------------------------------------------------------------*/
/* Updates the tracer in each window                                 */
/*-------------------------------------------------------------------*/
void tracer_step_window(master_t *master,
                        geometry_t *window,
                        window_t *windat, win_priv_t *wincon)
{
  /*---------------------------------------------------------------*/
  /* Reset diagnostic tracers to zero if required */
  tr_diag_reset_w(window, windat, wincon);
  if (wincon->trtend >= 0)
    memcpy(wincon->tendency, windat->tr_wc[wincon->trtend], 
	   window->sgsiz * sizeof(double));

  /*---------------------------------------------------------------*/
  /* Update the tracer concentrations. Note : density and vertical */
  /* mixing calculations are also included here.  */
  if (!master->mode2d)
    tracer_step_3d(window, windat, wincon);
  else
    tracer_step_2d(window, windat, wincon);
  if (master->crf == RS_RESTART) return;

  /*--------------------------------------------------------------*/
  /* Refill the master with tracer and density data from the */
  /* window data structure.  */
  master->win_data_empty_3d(master, window, windat, TRACERS);

#if !GLOB_BC
  /* Set the latereral ghost cells with density. This is only */
  /* neccessary for zoomed windows where a step in the bathymetry */
  /* creates a ghost cell that lies within two zoom cells from the */
  /* first valid cell to process adjacent to a solid boundary.  */
  if (window->zoomf > 1)
    Set_lateral_BC_density_w(windat->dens, window->nbpt, window->bpt,
                             window->bin);
#endif

}

/* END tracer_step_window()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to solve the tracer equation over all windows             */
/*-------------------------------------------------------------------*/
void tracer_step(master_t *master,  /* Master data structure */
                 geometry_t **window, /* Processing window */
                 window_t **windat, /* Window data structure */
                 win_priv_t **wincon, /* Window geometry / constants */
                 int nwindows   /* Number of windows to process */
  )
{
  int n, nn, tn;                /* Counter */

  /*-----------------------------------------------------------------*/
  /* Solve the tracer equation in each window */
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];

    /*---------------------------------------------------------------*/
    /* Set the lateral boundary conditions for tracers.  */
#if !GLOB_BC
    set_lateral_BC_tr(windat[n]->tr_wc, windat[n]->ntr, window[n]->nbpt,
                      window[n]->bpt, window[n]->bin);
    set_lateral_BC_vel(master->u1flux3d, window[n]->nbpt,
                       window[n]->bpte1, window[n]->bine1);
    set_lateral_BC_vel(master->u2flux3d, window[n]->nbpt,
                       window[n]->bpte2, window[n]->bine2);
#endif
    set_lateral_bdry_tr(window[n], windat[n], wincon[n]);

    /*---------------------------------------------------------------*/
    /* Set auxiliary cells that have been updated to the multi dt */
    /* buffers.  */
    fill_multidt_aux(master, window[n], windat[n]);
  }

  /* Invoke distributed processing step */
  dp_tracer_step();
  if (master->crf == RS_RESTART) return;

  /*-----------------------------------------------------------------*/
  /* Transfer mass diagnostics. This is done non-threaded, since if  */
  /* the transfers are performed within a thread, occasionally if    */
  /* multiple windows transfer to the global array at the same time, */
  /* and the global array variable sums over windows (as is the case */
  /* with total mass), then the global sum only adds mass from one   */
  /* window. This seems to be an issue with the threading libraries. */
  /* xxx fix this - TMASS not supported for MPI*/
  for (nn = 1; nn <= nwindows; nn++) {
    n = wincon[1]->twin[nn];
    master->win_data_empty_3d(master, window[n], windat[n], TMASS);
  }
  if (!master->mpi_rank && master->update_master)
    master->update_master(master, windat, TMASS|TRACERS);

  /*-----------------------------------------------------------------*/
  /* Interpolate on the master for zoomed grids */
  if (!(master->dozoom == NOZOOM)) {
    global_interp(geom, master->dens, geom->nzin);
    for (tn = 0; tn < master->ntr; tn++)
      global_interp(geom, master->tr_wc[tn], geom->nzin);
  }

  /*-----------------------------------------------------------------*/
  /* Get the flushing time if required */
  if (master->trflsh)
    calc_flushing(master, window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Print the region budgets if required */
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
void tracer_step_3d(geometry_t *window, /* Processing window */
		    window_t   *windat, /* Window data structure */
		    win_priv_t *wincon  /* Window geometry / constants */
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
  /* Do the advection and horizontal diffusion */
  TIMING_SET;
  get_bdry_cellno(window, windat, wincon);
  if (wincon->trsplit)
    advect_diffuse_split(window, windat, wincon);
  else if (wincon->trasc == LAGRANGE)
     advect_diffuse_lag(window, windat, wincon);
  else {
    if (advect_diffuse(window, windat, wincon)) return;
  }
  TIMING_DUMP_WIN(2, "  advect_diffuse", window->wn);
  debug_c(window, D_TS, D_ADVECT);
  if (wincon->trtend >= 0)
    get_tend(window, window->w3_t, window->b3_t, windat->tr_wc[wincon->trtend],
             wincon->tendency, windat->tr_adv);

  /* Reset the maps to map to velocity cells on R_EDGE and F_EDGE    */
  reset_map_t_all(window);

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion */
  TIMING_SET;
  vert_diffuse_3d(window, windat, wincon);
  debug_c(window, D_TS, D_HDIFF);
  if (wincon->trtend >= 0)
    get_tend(window, window->w3_t, window->b3_t, windat->tr_wc[wincon->trtend],
             wincon->tendency, windat->tr_vdif);
  TIMING_DUMP_WIN(2, "  vert_diffuse", window->wn);

  /*-----------------------------------------------------------------*/
  /* Do the diagnostics and interfaced library routines */
  TIMING_SET;
  auxiliary_routines(window, windat, wincon);
  TIMING_DUMP_WIN(2, "  aux_routines", window->wn);

  /*-----------------------------------------------------------------*/
  /* Set the open boundary conditions */
  bdry_tracer(window, windat, wincon);
  debug_c(window, D_TS, D_BDRY);

  /*-----------------------------------------------------------------*/
  /* Set the density */
  if (wincon->calc_dens)
    density_w(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the means for temperature and salinity and update the */
  /* mean counter.                                                   */
  reset_means(window, windat, wincon, RESET|TS);

  /*-----------------------------------------------------------------*/
  /* Calculate the means for tracers other than velocity */
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
void tracer_step_2d(geometry_t *window, /* Processing window */
                    window_t *windat, /* Window data structure */
                    win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int c, cc, c2;                /* Sparse coordinate / counter */

  /*-----------------------------------------------------------------*/
  /* Set the maps to be self-mapping across tracer OBC's */
  reset_map_t_all(window);

  /*-----------------------------------------------------------------*/
  /* Reset the fluxes over wet and ghost cells; substepping in the */
  /* tracer routine may use a different timestep to windat->dt.  */
  for (c = 1; c <= window->enonS; c++) {
    windat->u1flux[c] /= windat->dt;
    windat->u2flux[c] /= windat->dt;
  }
  /* Set the e1 velocity */
  for (cc = 1; cc <= window->a2_e1; cc++) {
    c = window->w2_e1[cc];
    c2 = window->m2d[c];
    windat->u1[c] = windat->u1av[c2];
    windat->u1flux3d[c] = windat->u1flux[c2];
  }
  /* Set the e2 velocity */
  for (cc = 1; cc <= window->a2_e2; cc++) {
    c = window->w2_e2[cc];
    c2 = window->m2d[c];
    windat->u2[c] = windat->u2av[c2];
    windat->u2flux3d[c] = windat->u2flux[c2];
  }

  /* Set the vertical velocity */
  memset(windat->w, 0, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Set up the cell centered dz */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    wincon->dz[c] = (windat->topz[c] - window->botz[c]) * wincon->Ds[c];
  }

  /*-----------------------------------------------------------------*/
  /* Do the advection and horizontal diffusion */
  if (advect_diffuse_2d(window, windat, wincon)) return;

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion */
  vert_diffuse_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the open boundary conditions */
  bdry_tracer(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Set the density */
  if (wincon->calc_dens)
    density_w(window, windat, wincon);
}

/* END tracer_step_2d()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion                                */
/*-------------------------------------------------------------------*/
int advect_diffuse(geometry_t *window, /* Processing window          */
		   window_t *windat,   /* Window data structure      */
		   win_priv_t *wincon  /* Window constants           */
  )
{
  int nn;                       /* Tracer index */
  int c, cc;                    /* Sparse indices */
  int c2;                       /* 2D cell corresponding to 3D cell */
  int cb;                       /* Bottom sparse coordinate (cell centre) */
  int vc;                       /* Tracer cells to process counter */
  int vcs;                      /* Surface tracer cells to process counter */
  int vc1;                      /* As for vc excluding OBC cells */
  int vcs1;                     /* As for vcs excluding OBC cells */
  int xp1;                      /* i index at i+1 */
  int yp1;                      /* j index at j+1 */
  int zp1;                      /* k index at k+1 */
  int zm1;                      /* k index at k-1 */
  int xm1, ym1;                 /* Sparse cell at i-1, j-1 and k-1 */
  double top;                   /* Flux through the surface cell */
  double surftr;                /* Surface tracer concentration */
  double dtu;                   /* Sub-time step to use */
  double dtm;                   /* Minimum allowable time step */
  double trem;                  /* Time remaining */
  double *osubeta;              /* Old surface elevation at sub-timestep */
  double *subeta;               /* Surface elevation at sub-timestep */
  double d1, d2, d3, d4;        /* Dummy variables */
  int ksf;                      /* Surface layer for sigma / Cartesian */
  double sf = 0.8;              /* Fraction of sub-timestep to use */
  double minval = 1e-10;        /* Minimum value for velocity */
  int slf = 1;                  /* Set to zero on the first sub-step */
  int itermax = 20;             /* Maximum number of substeps */
  int courf;                    /* Flag to calculate Courant numbers */
  int ii = 0, jj = 0, kk = 0;   /* (i,j,k) location of sub-step violation */
  double vm = 0.0, vs = 0.0;    /* Sub-step Courant number & grid spacing */
  char vt[4];                   /* Name of sub-step velocity component */
  double *u1, *u2, *w;
  double kzlim;

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* w1 = dtracer */
  /* w2 = Fx */
  /* w3 = Fy */
  /* w4 = Fz */
  /* w5 = Mean vertical cell spacing */
  /* w6 = x direction Courant number */
  /* w7 = y direction Courant number */
  /* w8 = z direction Courant number */
  /* w9 = vertical cell spacing between grid centers */
  /* d1 = osubeta */
  /* d2 = subeta */
  /* d3 = mean cell spacing in the x direction */
  /* d4 = mean cell spacing in the y direction */
  /* s1 = wet cells to process for tracers */

  /*-----------------------------------------------------------------*/
  /* Assign pointers */
  osubeta = wincon->d1;
  subeta = wincon->d2;
  vc = wincon->vc;
  vc1 = wincon->vc1;
  vcs = wincon->vcs;
  vcs1 = wincon->vcs1;
  u1 = windat->u1;
  u2 = windat->u2;
  w = windat->w;
  if (!(wincon->alertf & NONE)) {
    memcpy(windat->salb, windat->sal, window->sgsiz * sizeof(double));
    memcpy(windat->tempb, windat->temp, window->sgsiz * sizeof(double));
  }
  if (wincon->tmode & SP_FFSL) itermax = 200;

  /*-----------------------------------------------------------------*/
  /* Reset the fluxes for mean velocity transport */
  if (wincon->means & TRANSPORT) {
    if (windat->dttr == 0.0)
      return(0);
    set_flux_3d(window, windat, wincon, TRANSPORT);
    set_lateral_BC_vel(master->u1flux3d, window->nbpt,
                       window->bpte1, window->bine1);
    set_lateral_BC_vel(master->u2flux3d, window->nbpt,
                       window->bpte2, window->bine2);
    set_lateral_BC_vel(master->u1m, window->nbpt,
                       window->bpte1, window->bine1);
    set_lateral_BC_vel(master->u2m, window->nbpt,
                       window->bpte2, window->bine2);
    u1 = windat->u1m;
    u2 = windat->u2m;
    w = windat->wm;
  }

  /* Reset the river flow to a parabolic profile if required         */
  reset_flow(window, windat, wincon, 0);

  /*-----------------------------------------------------------------*/
  /* Set up advection scheme specific work arrays.  */
  if (wincon->trasc == ORDER1) {
    /* ORDER1 : Get the directional flux for upwind schemes (note : */
    /* these schemes are monotonic and should not be used with the */
    /* ULTIMATE filter.  */
    /* SIGMA : Don't multiply by the depth here since u1flux3d and */
    /* u2flux3d have already been depth multiplied.  */
    wincon->ultimate = 0;
    for (cc = 1; cc <= window->n3_t; cc++) {
      c = window->w3_t[cc];
      wincon->w6[c] =
        0.5 * (windat->u1flux3d[c] + fabs(windat->u1flux3d[c]));
      wincon->w7[c] =
        0.5 * (windat->u1flux3d[c] - fabs(windat->u1flux3d[c]));
      wincon->w8[c] =
        0.5 * (windat->u2flux3d[c] + fabs(windat->u2flux3d[c]));
      wincon->w9[c] =
        0.5 * (windat->u2flux3d[c] - fabs(windat->u2flux3d[c]));
    }
  } else if (wincon->trasc == LAGRANGE) {
    /* Find the origin of the streamline & grid Courant numbers */
    semi_lagrange_c(window, windat, wincon);
  } else if (wincon->trasc == VANLEER) {
    /* Get the mean cell spacings */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[window->zm1[c]]);
    }
  } else if (wincon->trasc == FFSL) {
    memset(wincon->nu, 0, window->sgsiz * sizeof(double));
    memset(wincon->nv, 0, window->sgsiz * sizeof(double));
    memset(wincon->nw, 0, window->sgsiz * sizeof(double));
    itermax = 200;
    /* Get the grid spacing. */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[window->zm1[c]]);
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      zp1 = window->zp1[c];

      /* Calculate cell-centered u1 velocity. The choice of velocity for each cell */
      /* can be either an average of the face values, or one of the face values */
      /* depending on those values - see Leonard et al. (1996). */
      xp1 = window->xp1[c];
      if ((u1[c] >= 0.0) && (u1[xp1] > 0.0)) 
        wincon->nu[c] = u1[c];
      else if ((u1[c] < 0.0) && (u1[xp1] <= 0.0))
        wincon->nu[c] = u1[xp1];

      /* Calculate cell-centered u2 velocity */
      yp1 = window->yp1[c];
      if ((u2[c] >= 0.0) && (u2[yp1] > 0.0)) 
        wincon->nv[c] = u2[c];
      else if ((u2[c] < 0.0) && (u2[yp1] <= 0.0))
        wincon->nv[c] = u2[yp1];

      /* Calculate cell-centered w velocity */
      zp1 = window->zp1[c];
      if ((w[c] >= 0.0) && (w[zp1] > 0.0)) 
        wincon->nw[c] = w[c];
      else if ((w[c] < 0.0) && (w[zp1] <= 0.0))
        wincon->nw[c] = w[zp1];
    }
  } else {
    /* Get the mean cell spacings */
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      zm1 = window->zm1[c];
      wincon->w5[c] = wincon->dz[zm1] / (wincon->dz[zm1] + wincon->dz[c]);
      wincon->w9[c] = 0.5 * (wincon->dz[c] + wincon->dz[zm1]);
    }
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      c2 = window->m2d[c];
      xm1 = window->xm1[c2];
      ym1 = window->ym1[c2];
      wincon->d3[c2] = window->h1acell[xm1] / (window->h1acell[xm1] +
                                               window->h1acell[c2]);
      wincon->d4[c2] = window->h2acell[ym1] / (window->h2acell[ym1] +
                                               window->h2acell[c2]);
    }
  }
  courf = 0;
  if (wincon->ultimate || wincon->trasc >= QUICKEST)
    courf = 1;

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on */
  /* the cell centers to calculate the maximum timesteps. The      */
  /* vertical sub-step for the surface layer is not included since */
  /* the vertical advective terms are not invoked in this layer.  */
  /* Note that the FFSL advection scheme is limited by the velocity */
  /* gradients (Lipschitz number) not the Courant number.           */
  dtm = dtu = windat->dt;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF | SUB_STEP_TRACER)) {

    for (cc = 1; cc <= vc; cc++) {
      c = wincon->s1[cc];       /* Wet cell to process */
      c2 = window->m2d[c];      /* 2D cell corresponding to 3D cell */
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      zp1 = window->zp1[c];

      /* Minimum time-step due to x velocity */
      d4 = 0.5 * (u1[c] + u1[xp1]);
      if (wincon->trasc == FFSL)
        d4 = max(fabs(u1[xp1] - u1[c]), wincon->u1kh[c]/window->h1acell[c2]);
      if (fabs(d4) < minval)
        d4 = minval;
      d1 = sf * fabs(window->h1acell[c2] / d4);
      if (d1 < dtm) {
        dtm = d1;
        sprintf(vt, "u");
        vm = d4;
        vs = window->h1acell[c2];
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }

      /* Minimum time-step due to y velocity */
      d4 = 0.5 * (u2[c] + u2[yp1]);
      if (wincon->trasc == FFSL)
        d4 = max(fabs(u2[yp1] - u2[c]), wincon->u2kh[c]/window->h2acell[c2]);
      if (fabs(d4) < minval)
        d4 = minval;
      d2 = sf * fabs(window->h2acell[c2] / d4);
      if (d2 < dtm) {
        dtm = d2;
        sprintf(vt, "v");
        vm = d4;
        vs = window->h2acell[c2];
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }

      /* Minimum time-step due to z velocity.  */
      /* Note : surface tracer values are calculated on the basis of */
      /* mass conservation, hence Courant violations need not be */
      /* considered.  */
      d4 = (cc <= vcs) ? 0.0 : 0.5 * (w[c] + w[zp1]);
      if (wincon->trasc == FFSL)
	d4 = (cc <= vcs) ? 0.0 : fabs(w[zp1] - w[c]);
      if (fabs(d4) < minval)
        d4 = minval;
      /* SIGMA : Multiply by depth */
      d3 = (cc <= vcs) ? dtu : sf * fabs(wincon->dz[c] * wincon->Ds[c2] / d4);
      if (d3 > 0.0 && d3 < dtm) {
        dtm = d3;
        sprintf(vt, "w");
        vm = d4;
        vs = wincon->dz[c] * wincon->Ds[c2];
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }
    }
    if (dtm != dtu) {
      dtm = floor(dtm);
      hd_warn
        ("Sub-time stepping for tracers (%s=%e) at %8.3f days (%3d %3d %3d): h=%e dt=%5.2f\n",
         vt, vm, windat->t / 86400.0, ii, jj, kk, vs, dtm);
    }
  } else if (wincon->stab & (MULTI_DT)) {
    dtm = windat->dts;
  }
  if (dtm == 0.0 || (int)(dtu / dtm) > itermax) {
    hd_quit_and_dump
      ("tracer_step_3d: maximum number of sub-steps (%d) exceeded at %8.3f days (%3d %3d %3d)\n",itermax, windat->t / 86400.0, ii, jj, kk);
      return(1);
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Include surface layer in the advection.  */
  ksf = vcs + 1;

  /* Set the system up for semi-Lagrangian advection */
  if (wincon->trasc == LAGRANGE) {
    /* Courant numbers may be > 1 */
    courf = 0;
    /* No special treatment of surface layer */
    if (!wincon->sigma)
      ksf = 1;
    /* No ultimate limiting */
    wincon->ultimate = 0;
    /* No sub-stepping (LAGRANGE is unconditionally stable) */
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
    /* Get the sub-timestep */
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
        for (cc = 1; cc <= window->n3_t; cc++) {
          c = window->w3_t[cc];
          c2 = window->m2d[c];
          wincon->w6[c] = u1[c] * dtu / window->h1au1[c2];
          wincon->w7[c] = u2[c] * dtu / window->h2au2[c2];
          wincon->w8[c] = w[c] * dtu / (wincon->w9[c] * wincon->Ds[c2]);
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the surface elevation for the sub-step.                 */
      /* Note : the tracer is updated in the surface layer in the    */
      /* main loop for the sigma system. Special treatment of the    */
      /* surface layer, and hence the sub-step elevation, is not     */
      /* required.                                                   */
      if (!wincon->sigma) {
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,c2)
#endif
        for (cc = 1; cc <= vcs; cc++) {
          c = wincon->s1[cc];
          c2 = window->m2d[c];
          subeta[c2] =
            wincon->oldeta[c2] + (1.0 + (dtu - trem) / windat->dt) *
	    (windat->eta[c2] - wincon->oldeta[c2]);
        }
      }

      /*-------------------------------------------------------------*/
      /* Save the region volume fluxes                               */
      region_volume_flux_coup(window, windat, wincon, dtu, trem);

      /*-------------------------------------------------------------*/
      /* If the FFSL scheme is used, calculate trajectories          */
      if (wincon->trasc == FFSL)
        prep_ff_sl(window, windat, wincon, dtu);

      /*-------------------------------------------------------------*/
      /* Tracer loop                                                 */
      //TIMING_SET;
#if defined(HAVE_OMP)
#pragma omp parallel for private(c,cc,cb,c2,zp1,surftr,top)
#endif
      for (nn = 0; nn < wincon->ntbdy; nn++) {
        int n = wincon->tbdy[nn];
        double *tr = windat->tr_wc[n];  /* Tracer values             */
	tr[0] = (double)n;
	int kef, kefs, bgzf;          /* Include / exclude OBCs */
	double *Fx = wincon->w2;      /* x tracer flux */
	double *Fy = wincon->w3;      /* y tracer flux */
	double *Fz = wincon->w4;      /* z tracer flux */
	double *dtracer = wincon->w1;
	/* Modified tracer concentrations */
	double *tr_mod = wincon->tr_mod;
	double *tr_mod_x = wincon->tr_mod_x;
	double *tr_mod_y = wincon->tr_mod_y;
	double *tr_mod_z = wincon->tr_mod_z;

	/* Re-wire work arrays if parallel mode */
#if defined(HAVE_OMP)
	if (master->trans_num_omp > 1) {
	  int nomp = omp_get_thread_num();
	  Fx = wincon->w2n[nomp];
	  Fy = wincon->w3n[nomp];
	  Fz = wincon->w4n[nomp];
	  dtracer = wincon->w1n[nomp];
	  tr_mod = wincon->tr_modn[nomp];
	  tr_mod_x = wincon->tr_modn_x[nomp];
	  tr_mod_y = wincon->tr_modn_y[nomp];
	  tr_mod_z = wincon->tr_modn_z[nomp];
	}
#endif	

        /*-----------------------------------------------------------*/
        /* Initialise the advective terms                            */
        memset(Fx, 0, window->sgsiz * sizeof(double));
        memset(Fy, 0, window->sgsiz * sizeof(double));
        memset(Fz, 0, window->sgsiz * sizeof(double));
        memset(dtracer, 0, window->sgsiz * sizeof(double));

        /*-----------------------------------------------------------*/
	/* Set whether to include boundary cells or not              */
	bgzf = 0;
	kef = vc1;
	kefs = vcs1;
	for (cc = 0; cc < window->nobc; cc++) {
	  if (window->open[cc]->bcond_tra[n] & (TRCONC|TRFLUX|TRCONF)) {
	    bgzf = 1;
	    kef = vc;
	    kefs = vcs;
	    save_OBC_tr(window, windat, wincon, Fx, Fy, tr, n, 1);
	    break;
	  }
	}

        /*-----------------------------------------------------------*/
        /* Set the bottom boundary condition (no-gradient)           */
        for (cc = 1; cc <= window->b2_t; cc++) {
          cb = window->bot_t[cc];
          tr[window->zm1[cb]] = tr[cb];
        }
	/* Ghost cells for sediments (if UPSTRM OBCs are used)       */
	if (wincon->trasc == FFSL) {
	  for (cc = 1; cc <= window->ngsed; cc++) {
	    cb = window->gsed_t[cc];
	    tr[cb] = tr[window->zp1[cb]];
	  }
	}

        /*-----------------------------------------------------------*/
        /* Get the tracer values in multi-dt auxiliary cells. This   */
        /* is a linear interpolation between the value at the start  */
        /* of the timestep and the end of longer timesteps.          */
        set_multidt_t(window, windat, tr, trem, n);

        /*-----------------------------------------------------------*/
        /* Get the advective fluxes                                  */
	TIMING_SET;
        if (wincon->advect[n])
          advect(window, windat, wincon, tr, Fx, Fy, Fz, 
		 tr_mod, tr_mod_x, tr_mod_y, tr_mod_z, dtu);
	TIMING_DUMP(2, "    advect");

        /*-----------------------------------------------------------*/
	/* Save the advective fluxes on OBCs                         */
	if (bgzf) save_OBC_tr(window, windat, wincon, Fx, Fy, tr, n, 3);

        /*-----------------------------------------------------------*/
        /* Save the advective fluxes if required                     */
	if (wincon->trflux == n)
	  calc_flux(window, windat, wincon, Fx, Fy, Fz, dtu);

	/* Reset the OBC fluxes if required. This is done in two     */
	/* parts; the fluxes in dummy OBC arrays and fluxes Fx or Fy */
	/* are set here so that regional fluxes can be computed      */
	/* using these modified fluxes below. Horizontal mixing then */
	/* alters Fx and Fy, and dtracer is then modified by another */
	/* call to reset_tr_OBCflux() to reflect only the OBC        */
	/* specified fluxes.                                         */
	reset_tr_OBCflux(window, windat, wincon, dtracer, Fx, Fy, dtu, n, 1);

        /*-----------------------------------------------------------*/
	/* Save the region fluxes                                    */
	region_flux_coup(window, windat, wincon, Fx, Fy, Fz, dtu, n);

        /*-----------------------------------------------------------*/
        /* Get the horizontal diffusive fluxes                       */
        if (wincon->diffuse[n])
          hor_diffuse(window, windat, wincon, tr, Fx, Fy);

        /*-----------------------------------------------------------*/
	/* Reset the advective fluxes on OBCs if required            */
	if (bgzf) save_OBC_tr(window, windat, wincon, Fx, Fy, tr, n, 4);

        /*-----------------------------------------------------------*/
        /* Store the updated horizontal advective fluxes             */
        /* Note : land cells and cells immediately adjacent to open  */
        /* boundaries are not updated.                               */
        for (cc = 1; cc <= window->b3_t; cc++) {
          c = window->w3_t[cc];
          dtracer[c] = (Fx[window->xp1[c]] - Fx[c] +
                        Fy[window->yp1[c]] - Fy[c]) * dtu;	  
        }

        /*-----------------------------------------------------------*/
	/* Reset the flux divergence at open boundary cells (if      */
	/* required).                                                */
	reset_tr_OBCflux(window, windat, wincon, dtracer, Fx, Fy, dtu, n, 0);

        /*-----------------------------------------------------------*/
        /* Evaluate the sources and sinks of tracer                  */
	ss_tracer(window, windat, wincon, n, dtracer, dtu);
        if ((wincon->trasc == FFSL) && (windat->nstep == 0))
	   memset(dtracer, 0, window->sgsiz * sizeof(double));

        /*-----------------------------------------------------------*/
        /* Get the updated tracer concentration                      */
        /* Note : land cells and cells immediately adjacent to open  */
        /* boundaries are not updated. Note, the first vcs cells in  */
        /* the cells to process vector, wincon->s1, contain surface  */
        /* cells and are not processed here since the surface has    */
        /* separate treatment.                                       */
        for (cc = ksf; cc <= kef; cc++) {
          c = wincon->s1[cc];   /* Wet cell to process               */
          c2 = window->m2d[c];  /* 2D cell corresponding to 3D cell  */
          zp1 = window->zp1[c]; /* Cell above cell c                 */

          /* SIGMA : Adjust tracer values for the depth              */
          if (slf)
            tr[c] *= wincon->Ds[c2];
          else
            tr[c] *= wincon->Hn1[c2];
          tr[c] -= ((dtracer[c] / window->cellarea[c2] +
                     (Fz[zp1] - Fz[c])) / wincon->dz[c]);
          tr[c] /= wincon->Hn1[c2];
	}

        /*-----------------------------------------------------------*/
        /* Set all layers above the surface to surface layer         */
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
        /*-----------------------------------------------------------*/
        /* Get the tracer concentration in the surface cells         */
        else {
          if (!wincon->sigma && wincon->advect[n]) {
            for (cc = 1; cc <= kefs; cc++) {
              c = wincon->s1[cc];
              c2 = window->m2d[c];
              top = -Fz[c] * window->cellarea[c2];
	      surftr =
		surf_conc(window, windat, wincon, tr, subeta[c2], c, c2,
			  cc, top, dtracer);
              /* Set all layers above the surface to surface layer   */
	      tr[c] = (window->botz[c2] > 0.0 && wincon->dz[c] < wincon->hmin) ?
		tr[c] : surftr;
	      /*
              while (c != window->zp1[c]) {
                c = window->zp1[c];
                tr[c] = surftr;
              }
	      */
	    }
          }
          if (wincon->sigma && wincon->advect[n]) {
            /* SIGMA surface calculation : w=0 at the free surface   */
            /* Fz[zp1]=0 at this layer.                              */
            for (cc = 1; cc < ksf; cc++) {
              c = wincon->s1[cc];
              c2 = window->m2d[c];

              /* Diagnostic to check for continuity using tracer     */
              /* fluxes (tracer value must = constant).              */
              /*
                 d1=dtracer[c]/(wincon->dz[c]*tr[c]);
                 d2=-Fz[c]*window->cellarea[c2]/(wincon->dz[c]*tr[c]);
                 d3=(wincon->Hn1[c2]-wincon->Ds[c2])*window->cellarea[c2];
                 top=d1+d2+d3; */
              /* Diagnostic to check for continuity using momentum */
              /* fluxes.  */
              /*
                 d1=wincon->dz[c]*window->cellarea[c2]*
                 (wincon->Hn1[c2]-wincon->Ds[c2])/windat->dt;
                 d2=windat->u1flux3d[window->xp1[c]]-windat->u1flux3d[c]+
                 windat->u2flux3d[window->yp1[c]]-windat->u2flux3d[c];
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
	/* Reset boundary tracer values if required */
	if (bgzf) save_OBC_tr(window, windat, wincon, Fx, Fy, tr, n, 2);
      }                         /* Tracer loop end */
      //TIMING_DUMP(2, "    t-loop");

      /* Clip the tracer concentration if required */
      TIMING_SET;
      tr_bounds(window, windat, wincon);
      TIMING_DUMP(2, "    tr_bounds");


      /* Set the old sub-stepped elevation */
      if (!wincon->sigma) {
        for (cc = 1; cc <= vcs; cc++) {
          c = wincon->s1[cc];
          c2 = window->m2d[c];
          osubeta[c2] = subeta[c2];
        }
      }
    }
    /* dtu!=0 loop */
    /* Decrement the time remaining */
    trem -= dtu;
    slf = 0;
  }                             /* while(trem>0) loop */

  /*-----------------------------------------------------------------*/
  /* Reset the new surface coordinate at thin layers if required */
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
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	density_gc(window, windat, wincon, c, open->omap);
      }
    }
  }

  /* Revert river flow to original profile if required               */
  reset_flow(window, windat, wincon, 1);

  return(0);
}

/* END advect_diffuse()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Advection and horizontal diffusion                                */
/*-------------------------------------------------------------------*/
int advect_diffuse_2d(geometry_t *window,  /* Processing window */
		      window_t *windat,  /* Window data structure */
		      win_priv_t *wincon /* Window geometry / constants */
  )
{
  int nn;                       /* Tracer index */
  int c, cc;                    /* Sparse indices */
  int vcs;                      /* Surface tracer cells to process counter
                                 */
  int vcs1;                     /* As for vcs excluding OBC cells */
  int xp1;                      /* i index at i+1 */
  int yp1;                      /* j index at j+1 */
  int xm1, ym1;                 /* Sparse cell at i-1, j-1 and k-1 */
  double *Fx;                   /* x tracer flux */
  double *Fy;                   /* y tracer flux */
  double *Fz;                   /* z tracer flux */
  double *dtracer;              /* Advective terms for x,y directions */
  double dtu;                   /* Sub-time step to use */
  double dtm;                   /* Minimum allowable time step */
  double trem;                  /* Time remaining */
  double *osubeta;              /* Old surface elevation at sub-timestep */
  double *subeta;               /* Surface elevation at sub-timestep */
  double d1, d2, d4;            /* Dummy variables */
  double sf = 0.8;              /* Fraction of sub-timestep to use */
  double minval = 1e-10;        /* Minimum value for velocity */
  int slf = 1;                  /* Set to zero on the first sub-step */
  int itermax = 20;             /* Maximum number of substeps */
  int courf;                    /* Flag to calculate Courant numbers */
  int ii = 0, jj = 0, kk = 0;   /* (i,j,k) location of sub-step violation */
  double vm = 0.0, vs = 0.0;    /* Sub-step Courant number & grid spacing */
  char vt[4];                   /* Name of sub-step velocity component */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* w1 = dtracer */
  /* w2 = Fx */
  /* w3 = Fy */
  /* w4 = Fz */
  /* w5 = Mean vertical cell spacing */
  /* w6 = x direction Courant number */
  /* w7 = y direction Courant number */
  /* w8 = z direction Courant number */
  /* w9 = vertical cell spacing between grid centers */
  /* d1 = osubeta */
  /* d2 = subeta */
  /* d3 = mean cell spacing in the x direction */
  /* d4 = mean cell spacing in the y direction */
  /* s1 = wet cells to process for tracers */

  /*-----------------------------------------------------------------*/
  /* Assign memory for the horizontal advective terms */
  dtracer = wincon->w1;
  Fx = wincon->w2;
  Fy = wincon->w3;
  Fz = wincon->w4;
  osubeta = wincon->d1;
  subeta = wincon->d2;
  vcs = wincon->vcs;
  vcs1 = wincon->vcs1;

  /*-----------------------------------------------------------------*/
  /* Set up advection scheme specific work arrays.  */
  if (wincon->trasc == ORDER1) {
    /* ORDER1 : Get the directional flux for upwind schemes (note : */
    /* these schemes are monotonic and should not be used with the */
    /* ULTIMATE filter.  */
    /* SIGMA : Don't multiply by the depth here since u1flux3d and */
    /* u2flux3d have already been depth multiplied.  */
    wincon->ultimate = 0;
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      wincon->w6[c] = 0.5 * (windat->u1flux[c] + fabs(windat->u1flux[c]));
      wincon->w7[c] = 0.5 * (windat->u1flux[c] - fabs(windat->u1flux[c]));
      wincon->w8[c] = 0.5 * (windat->u2flux[c] + fabs(windat->u2flux[c]));
      wincon->w9[c] = 0.5 * (windat->u2flux[c] - fabs(windat->u2flux[c]));
    }
  } else if (wincon->trasc == VANLEER) {
    /* Get the mean cell spacings */
    memcpy(wincon->w9, wincon->dz, window->sgsizS * sizeof(double));
  } else if (wincon->trasc == FFSL) {
    /* Get the mean cell spacings */
    memcpy(wincon->w9, wincon->dz, window->sgsizS * sizeof(double));
  } else {
    memcpy(wincon->w5, wincon->dz, window->sgsizS * sizeof(double));
    memcpy(wincon->w9, wincon->dz, window->sgsizS * sizeof(double));
    for (cc = 1; cc <= window->n2_t; cc++) {
      c = window->w2_t[cc];
      xm1 = window->xm1[c];
      ym1 = window->ym1[c];
      wincon->d3[c] = window->h1acell[xm1] / (window->h1acell[xm1] +
                                              window->h1acell[c]);
      wincon->d4[c] = window->h2acell[ym1] / (window->h2acell[ym1] +
                                              window->h2acell[c]);
    }
  }
  courf = 0;
  if (wincon->ultimate || wincon->trasc >= QUICKEST)
    courf = 1;

  /*-----------------------------------------------------------------*/
  /* Get the minimum sub-timestep. Note: velocities are centered on */
  /* the cell centers to calculate the maximum timesteps. The */
  /* vertical sub-step for the surface layer is not included since */
  /* the vertical advective terms are not invoked in this layer.  */
  dtm = dtu = windat->dt;
  if (wincon->stab & (SUB_STEP | SUB_STEP_NOSURF | SUB_STEP_TRACER)) {

    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D cell corresponding to 3D cell */
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];

      /* Minimum time-step due to x velocity */
      d4 = 0.5 * (windat->u1av[c] + windat->u1av[xp1]);
      if (fabs(d4) < minval)
        d4 = minval;
      d1 = sf * fabs(window->h1acell[c] / d4);
      if (d1 < dtm) {
        dtm = d1;
        sprintf(vt, "u");
        vm = d4;
        vs = window->h1acell[c];
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
      }

      /* Minimum time-step due to y velocity */
      d4 = 0.5 * (windat->u2av[c] + windat->u2av[yp1]);
      if (fabs(d4) < minval)
        d4 = minval;
      d2 = sf * fabs(window->h2acell[c] / d4);
      if (d2 < dtm) {
        dtm = d2;
        sprintf(vt, "v");
        vm = d4;
        vs = window->h2acell[c];
        ii = window->s2i[c];
        jj = window->s2j[c];
        kk = window->s2k[c];
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
  /* Initialise the surface elevation. The surface concentration */
  /* calculated in surf_concs() must also be sub-stepped. This means */
  /* the volume of the previous sub-step must be known to get the */
  /* concentration at the current sub-step. The elevations from the */
  /* previous sub-steps are stored in osubeta, allowing previous */
  /* sub-step volumes to be calculated. It is assumed that the */
  /* velocities contributing to the fluxes into the surface layer */
  /* are invariant over the full time-step, windat->dt.  */
  if (wincon->sigma)
    osubeta = wincon->one;
  else {
    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D surface wet cell to process */
      osubeta[c] = wincon->oldeta[c] - window->botz[c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* SIGMA : Get the water depth at the forward time */
  if (wincon->sigma)
    for (cc = 1; cc <= window->v2_t; cc++) {
      c = window->w2_t[cc];     /* 2D surface wet cell to process */
      wincon->Hn1[c] = windat->eta[c] - wincon->Hs[c];
    }

  /*-----------------------------------------------------------------*/
  /* Loop over the time interval until no time remains */
  trem = dtu;
  while (trem > 0) {
    /* Get the sub-timestep */
    dtu = dtm;
    if (trem < dtu)
      dtu = trem;

    if (dtu > 0.0) {

      /*-------------------------------------------------------------*/
      /* Get the courant numbers for this sub-step */
      if (courf) {
        for (cc = 1; cc <= window->n2_t; cc++) {
          c = window->w2_t[cc];
          wincon->w6[c] = windat->u1av[c] * dtu / window->h1au1[c];
          wincon->w7[c] = windat->u2av[c] * dtu / window->h2au2[c];
          wincon->w8[c] = 0.0;
        }
      }

      /*-------------------------------------------------------------*/
      /* Get the surface elevation for the sub-step.  */
      /* Note : the tracer is updated in the surface layer in the */
      /* main loop for the sigma system. Special treatment of the */
      /* surface layer, and hence the sub-step elevation, is not */
      /* required.  */
      if (!wincon->sigma) {
        for (cc = 1; cc <= window->v2_t; cc++) {
          c = window->w2_t[cc];
          subeta[c] =
            wincon->oldeta[c] + (1.0 +
                                 (dtu -
                                  trem) / windat->dt) * (windat->eta[c] -
                                                         wincon->
                                                         oldeta[c]) -
            window->botz[c];
        }
      } else
        subeta = wincon->one;

      /*-------------------------------------------------------------*/
      /* Tracer loop */
      for (nn = 0; nn < wincon->ntbdy; nn++) {
        int n = wincon->tbdy[nn];
        double *tr = windat->tr_wc[n];  /* Tracer values             */

        /*-----------------------------------------------------------*/
        /* Initialise the advective terms                            */
        memset(Fx, 0, window->sgsiz * sizeof(double));
        memset(Fy, 0, window->sgsiz * sizeof(double));
        memset(Fz, 0, window->sgsiz * sizeof(double));
        memset(dtracer, 0, window->sgsiz * sizeof(double));

        /*-----------------------------------------------------------*/
        /* Get the tracer values in multi-dt auxiliary cells. This   */
        /* is a linear interpolation between the value at the start  */
        /* of the timestep and the end of longer timesteps.          */
        set_multidt_t(window, windat, tr, trem, n);

        /*-----------------------------------------------------------*/
        /* Get the advective fluxes                                  */
        if (wincon->advect[n])
          advect(window, windat, wincon, tr, Fx, Fy, Fz, 
		 NULL, NULL, NULL, NULL, dtu);

        /*-----------------------------------------------------------*/
        /* Get the horizontal diffusive fluxes                       */
        if (wincon->diffuse[n])
          hor_diffuse_2d(window, windat, wincon, tr, Fx, Fy);

        /*-----------------------------------------------------------*/
        /* Store the updated horizontal advective fluxes             */
        /* Note : land cells and cells immediately adjacent to open  */
        /* boundaries are not updated.                               */
        for (cc = 1; cc <= window->v2_t; cc++) {
          c = window->w2_t[cc];
          dtracer[c] = (Fx[window->xp1[c]] - Fx[c] +
                        Fy[window->yp1[c]] - Fy[c]) * dtu;
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
  /* Flux diagnostics are calculated in advect_diffuse() using the   */
  /* actual fluxes used to update tracer concentration (i.e.         */
  /* advection scheme dependent).                                    */
  /*
  if (wincon->trflux >= 0)
    calc_flux_old(window, windat, wincon, dtu);
  */

  /*-----------------------------------------------------------------*/
  /* Calculate the flushing diagnostic if required */
  if (wincon->trflsh)
    total_mass(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Evaluate the tracer decay */
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
/* scheme and the remaining tracers advected using LAGRANGE. The     */
/* arrays of tracers to advect, tbdy, are modified for two calls to  */
/* the routine advect_diffuse().                                     */
/*-------------------------------------------------------------------*/
void advect_diffuse_split(geometry_t *window, /* Window geometry     */
			  window_t *windat,   /* Window data         */
			  win_priv_t *wincon  /* Window constants    */
  )
{
  int *tbdy_o;
  int ntbdy_o;
  int trasco, uo;
  int nt = -1, ns = -1, n, m;
  double dto;

  /* Copy the original tracers to advect and diffuse vector          */
  ntbdy_o = wincon->ntbdy;
  tbdy_o = i_alloc_1d(ntbdy_o);
  memcpy(tbdy_o, wincon->tbdy, ntbdy_o * sizeof(int));
  for (n = 0; n < ntbdy_o; n++) {
    if (windat->sno == tbdy_o[n])ns = n;
    if (windat->tno == tbdy_o[n])nt = n;
  }
  trasco = wincon->trasc;
  uo = wincon->ultimate;
  dto = windat->dttr;

  /* Set up the tracers to advect using LAGRANGE                     */
  m = 0;
  for (n = 0; n < ntbdy_o; n++) {
    if (tbdy_o[n] != nt && tbdy_o[n] != ns) {
      wincon->tbdy[m] = tbdy_o[n];
      m++;
    }
  }
  wincon->ntbdy = m;
  wincon->trasc = LAGRANGE;
  wincon->ultimate = 0;
  advect_diffuse(window, windat, wincon);

  /* Set up the tracers to advect using VANLEER                      */
  if (ns >=0 && nt >=0) {
    wincon->ntbdy = 2;
    wincon->tbdy[0] = ns;
    wincon->tbdy[1] = nt;
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
    advect_diffuse(window, windat, wincon);
  }

  /* Reset the original tracers to advect and diffuse vector         */
  memcpy(wincon->tbdy, tbdy_o, ntbdy_o * sizeof(int));
  wincon->ntbdy = ntbdy_o;
  wincon->trasc = trasco;
  windat->dttr = dto;
  i_free_1d(tbdy_o);
}

/* END advect_diffuse_split()                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the advective fluxes using the appropriate   */
/* advection scheme and invoke the ULTIMATE limiter if required.     */
/*-------------------------------------------------------------------*/
void advect(geometry_t *window, /* Processing window */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon, /* Window geometry / constants */
            double *tr,         /* Tracer array */
            double *Fx,         /* x direction flux */
            double *Fy,         /* y direction flux */
            double *Fz,         /* z direction flux */
	    double *tr_mod, 
	    double *tr_mod_x, 
	    double *tr_mod_y, 
	    double *tr_mod_z, 
            double dtu          /* Time-step */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  double *w;                    /* Vertical velocity */

  if (windat->dttr == 0.0)
    return;

  /*-----------------------------------------------------------------*/
  /* Get the tracer concentration at the cell faces */
  if (wincon->trasc == ORDER1)
    order1(window, windat, wincon, tr, Fx, Fy, Fz, dtu);
  else if (wincon->trasc == FFSL)
    ff_sl3d(window, windat, wincon, tr, Fx, Fy, Fz, tr_mod, tr_mod_x, tr_mod_y, tr_mod_z, dtu);
  else {
    if (wincon->trasc == ORDER2)
      order2(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == ORDER4)
      order4(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == QUICKEST)
      quickest(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == QUICKEST_CO)
      quickest_uniform(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == VANLEER)
      van_leer(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == ORDER2_UW)
      order2_upwind(window, windat, wincon, tr, Fx, Fy, Fz);
    else if (wincon->trasc == LAGRANGE)
      semi_lagrange(window, windat, wincon, tr);

    /*---------------------------------------------------------------*/
    /* Invoke the ULTIMATE limiter. This procedure is split into */
    /* doing the wet cells using velocity work arrays and auxiliary */
    /* cells using the tracer work arrays so that ULTIMATE (which */
    /* computationally expensive) is not invoked on cell faces */
    /* which are subsequently set to zero as a consequence of being */
    /* associated with zero velocity.  */
    if (wincon->ultimate) {
      /* Wet cell faces (x direction) in this window */
      ulimate_filter(wincon->w6, Fx, tr, window->w3_e1, 1, window->v3_e1,
                     window->xm1, window->xp1);
      /* Auxiliary cell faces (x direction) in this window */
      ulimate_filter(wincon->w6, Fx, tr, window->w3_t, window->v3_t + 1,
                     window->a3_t, window->xm1, window->xp1);

      /* Wet cell faces (y direction) in this window */
      ulimate_filter(wincon->w7, Fy, tr, window->w3_e2, 1, window->v3_e2,
                     window->ym1, window->yp1);
      /* Auxiliary cell faces (y direction) in this window */
      ulimate_filter(wincon->w7, Fy, tr, window->w3_t, window->v3_t + 1,
                     window->a3_t, window->ym1, window->yp1);

      /* Vertical faces */
      ulimate_filter(wincon->w8, Fz, tr, wincon->s1, 1, wincon->vc,
                     window->zm1, window->zp1);
    }

    /*---------------------------------------------------------------*/
    /* Get the tracer flux through the face.  */
    /* SIGMA : Don't multiply by the depth here since u1flux3d and */
    /* u2flux3d have already been depth multiplied.  */
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
        Fx[c] *= windat->u1flux3d[c];
	Fy[c] *= windat->u2flux3d[c];
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
void order1(geometry_t *window, /* Window structure */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon, /* Window geometry / constants */
            double *tr,         /* Tracer array */
            double *Fx,         /* x tracer flux */
            double *Fy,         /* y tracer flux */
            double *Fz,         /* z tracer flux */
            double dtu          /* Time-step */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int zm1;                      /* Sparse cell at k-1 */
  double *w;                    /* Vertical velocity */

  /* Horizontal fluxes */
  for (cc = 1; cc <= window->n3_t; cc++) {
    c = window->w3_t[cc];

    /* Get the fluxes */
    Fx[c] = wincon->w6[c] * tr[window->xm1[c]] + wincon->w7[c] * tr[c];
    Fy[c] = wincon->w8[c] * tr[window->ym1[c]] + wincon->w9[c] * tr[c];
  }

  /* Vertical fluxes */
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
/* Routine to calculate the second order advective fluxes.           */
/* Formulated for non-uniform grids.                                 */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->d4 = mean cell spacing in the y direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void order2(geometry_t *window, /* Window structure */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon, /* Window geometry / constants */
            double *tr,         /* Tracer array */
            double *Fx,         /* x tracer flux */
            double *Fy,         /* y tracer flux */
            double *Fz          /* z tracer flux */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int zm1;                      /* Sparse cell at k-1 */

  /* Horizontal fluxes */
  /* x direction horizontal fluxes. This calculation sets the tracer */
  /* concentration on the cell face (using 2nd order method). Note */
  /* that some faces (e.g. faces associated with cells on southern */
  /* or western land boundary edges) will always have zero fluxes */
  /* since velocity is zero on these faces. It is inefficient to set */
  /* the tracer concentration on the face if the cell is set to zero */
  /* when multiplied by the volume flux. The flux array is first */
  /* initialised to zero so that any cells not explicitly set to a */
  /* non-zero value (e.g. ghost cells) remain zero valued. The */
  /* velocity work arrays are then used to set tracer concentrations */
  /* only at wet faces in the window. Auxiliary cell tracer */
  /* concentrations are then set using the tracer work arrays. This */
  /* method ensures that only those faces which have non-zero tracer */
  /* concentrations are actually calculated. Tracer concentrations */
  /* on the faces must be multiplied by the volume flux through the */
  /* face in the main routine to get the tracer flux. This is done */
  /* over all wet + auxiliary cells using the tracer work arrays.  */
  /* This may lead to some inefficiency by multiplying zero valued */
  /* tracer concentrations by zero on southern and western edges */
  /* (which is not severe as only one operation is involved).  */
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  /* Do the wet cell faces in this window.  */
  order2_do(Fx, tr, wincon->d3, 1, window->v3_e1, window->w3_e1,
            window->xm1, window->m2d);
  /* Do the auxiliary cell faces associated with this window. Open */
  /* boundary locations are picked up here also.  */
  order2_do(Fx, tr, wincon->d3, window->v3_t + 1, window->a3_t,
            window->w3_t, window->xm1, window->m2d);

  /* y direction horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fy, 0, window->n3_t * sizeof(double));
  /* Wet cell faces in this window */
  order2_do(Fy, tr, wincon->d4, 1, window->v3_e2, window->w3_e2,
            window->ym1, window->m2d);
  /* Auxiliary cell faces associated with this window */
  order2_do(Fy, tr, wincon->d4, window->v3_t + 1, window->a3_t,
            window->w3_t, window->ym1, window->m2d);

  /* Vertical fluxes */
  memset(Fz, 0, window->n3_t * sizeof(double));
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zm1 = window->zm1[c];
    Fz[c] = wincon->w5[c] * (tr[c] - tr[zm1]) + tr[zm1];
  }
}

/* END order2()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* second order centered flux scheme. This scheme linearly           */
/* interpolates between cell centers to get the face value and is    */
/* formulated for non-uniform grids using the mean grid spacing      */
/* arrays. Only horizontal fluxes are calculated here; the vertical  */
/* flux requires a 3D mean grid spacing array whereas only a 2D      */
/* array is used here.                                               */
/*-------------------------------------------------------------------*/
void order2_do(double *F,       /* Tracer concentration at the face */
               double *tr,      /* Tracer array */
               double *mgs,     /* Mean grid spacing */
               int ss,          /* Start of the processing vector */
               int se,          /* End of the processing vector */
               int *sdo,        /* Cells to process vector */
               int *bmap,       /* Backward spatial map */
               int *m2d         /* 3D - 2D spatial map */
  )
{
  int c, cc, c2;                /* Sparse counters / coordinates */
  int bm1;                      /* Backward sparse coordinate */

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    c2 = m2d[c];
    bm1 = bmap[c];

    /* Get the fluxes */
    F[c] = mgs[c2] * (tr[c] - tr[bm1]) + tr[bm1];
  }
}

/* END order2_do()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the fourth order advective fluxes.           */
/* Formulated for uniform grids.                                     */
/*-------------------------------------------------------------------*/
void order4(geometry_t *window, /* Window structure */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon, /* Window geometry / constants */
            double *tr,         /* Tracer array */
            double *Fx,         /* x tracer flux */
            double *Fy,         /* y tracer flux */
            double *Fz          /* z tracer flux */
  )
{

  /* x direction horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  /* Do the wet cell faces in this window.  */
  order4_do(Fx, tr, 1, window->v3_e1, window->w3_e1, window->xm1,
            window->xp1);
  /* Do the auxiliary cell faces associated with this window */
  order4_do(Fx, tr, window->v3_t + 1, window->a3_t, window->w3_t,
            window->xm1, window->xp1);

  /* y direction horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fy, 0, window->n3_t * sizeof(double));
  /* Do the wet cell faces in this window.  */
  order4_do(Fy, tr, 1, window->v3_e2, window->w3_e2, window->ym1,
            window->yp1);
  /* Do the auxiliary cell faces associated with this window */
  order4_do(Fy, tr, window->v3_t + 1, window->a3_t, window->w3_t,
            window->ym1, window->yp1);

  /* Vertical fluxes */
  memset(Fz, 0, window->n3_t * sizeof(double));
  order4_do(Fz, tr, 1, wincon->vc, wincon->s1, window->zm1, window->zp1);
}

/* END order4()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using the  */
/* second order centered flux scheme.                                */
/*-------------------------------------------------------------------*/
void order4_do(double *F,       /* Tracer concentration at the face */
               double *tr,      /* Tracer array */
               int ss,          /* Start of the processing vector */
               int se,          /* End of the processing vector */
               int *sdo,        /* Cells to process vector */
               int *bmap,       /* Backward spatial map */
               int *fmap        /* Forward spatial map */
  )
{
  int c, cc;                    /* Sparse counters / coordinates */
  int bm1, bm2;                 /* Backward sparse coordinate */

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    bm1 = bmap[c];
    bm2 = bmap[bm1];

    /* Get the fluxes */
    F[c] = (7.0 * (tr[bm1] + tr[c]) - (tr[fmap[c]] + tr[bm2])) / 12.0;
  }
}

/* END order4_do()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the third order advective fluxes using the   */
/* QUICKEST algorithm. Formulated for non-uniform grids.             */
/* Note : wincon->d3 = mean cell spacing in the x direction          */
/*        wincon->d4 = mean cell spacing in the y direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*        wincon->w9 = vertical grid spacing between cell centers    */
/*-------------------------------------------------------------------*/
void quickest(geometry_t *window, /* Window structure */
              window_t *windat, /* Window data structure */
              win_priv_t *wincon, /* Window geometry / constants */
              double *tr,       /* Tracer array */
              double *Fx,       /* x tracer flux */
              double *Fy,       /* y tracer flux */
              double *Fz        /* z tracer flux */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int c2;                       /* 2D cell corresponding to 3D cell */
  int xm1, ym1, zm1;            /* Sparse cell at i-1, j-1 and k-1 */
  int xm2, ym2, zm2;            /* Sparse cell at i-2, j-2 and k-2 */
  int xp1, yp1, zp1;            /* Sparse cell at i+1, j+1 and k+1 */
  int xm1s, xp1s;               /* 2D sparse cell at i-1, i+1 */
  int ym1s, yp1s;               /* 2D sparse cell at j-1, j+1 */
  double face;                  /* Concentration value at the faces */
  double grad;                  /* Concentration difference across a face */
  double curv;                  /* Concentration curvature across a face */
  double cx, cy, cz;            /* Courant number in the x,y,z directions */
  double db, df;                /* Mean face centered grid spacings */

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  memset(Fy, 0, window->n3_t * sizeof(double));
  memset(Fz, 0, window->n3_t * sizeof(double));
  /* Horizontal fluxes */
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    xm2 = window->xm1[xm1];
    ym2 = window->ym1[ym1];
    xm1s = window->xm1[c2];
    xp1s = window->xp1[c2];
    ym1s = window->ym1[c2];
    yp1s = window->yp1[c2];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the left face.  */
    /* The face value is formulated for a non-uniform grid by */
    /* linearly interpolating between cell centers and evaluating on */
    /* the face.  */
    cx = wincon->w6[c];
    db = 2.0 / (window->h1au1[xm1s] + window->h1au1[c2]);
    df = 2.0 / (window->h1au1[c2] + window->h1au1[xp1s]);
    face = wincon->d3[c2] * (tr[c] - tr[xm1]) + tr[xm1];
    grad = tr[c] - tr[xm1];
    curv = 0.0;

    /* Note: the curvature is formulated for a non-uniform grid (see */
    /* Kowalik & Murty (1993), p40) and multiplied by (h1au1[i] * */
    /* h1au1[i-1]) to give it dimensions of [tracer].  */
    if (cx > 0.0)
      curv =
        db * (window->h1au1[xm1s] * tr[c] + window->h1au1[c2] * tr[xm2]) -
        2.0 * tr[xm1];
    else if (cx < 0.0)
      curv =
        df * (window->h1au1[c2] * tr[xp1] +
              window->h1au1[xp1s] * tr[xm1]) - 2.0 * tr[c];
    Fx[c] = face - 0.5 * cx * grad - (1.0 - cx * cx) * curv / 6.0;

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the back face */
    cy = wincon->w7[c];
    db = 2.0 / (window->h2au2[ym1s] + window->h2au2[c2]);
    df = 2.0 / (window->h2au2[c2] + window->h2au2[yp1s]);
    face = wincon->d4[c2] * (tr[c] - tr[ym1]) + tr[ym1];
    grad = tr[c] - tr[ym1];
    curv = 0.0;
    if (cy > 0.0)
      curv =
        db * (window->h2au2[ym1s] * tr[c] + window->h2au2[c2] * tr[ym2]) -
        2.0 * tr[ym1];
    else if (cy < 0.0)
      curv =
        df * (window->h2au2[c2] * tr[yp1] +
              window->h2au2[yp1s] * tr[ym1]) - 2.0 * tr[c];
    Fy[c] = face - 0.5 * cy * grad - (1.0 - cy * cy) * curv / 6.0;
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the lower face */
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
/*        wincon->d4 = mean cell spacing in the y direction          */
/*        wincon->w5 = mean cell spacing in the z direction          */
/*-------------------------------------------------------------------*/
void quickest_uniform(geometry_t *window, /* Window structure */
                      window_t *windat, /* Window data structure */
                      win_priv_t *wincon, /* Window geometry / constants */
                      double *tr, /* Tracer array */
                      double *Fx, /* x tracer flux */
                      double *Fy, /* y tracer flux */
                      double *Fz  /* z tracer flux */
  )
{
  int c, cc;                    /* Sparse coordinate / counter */
  int xm1, ym1, zm1;            /* Sparse cell at i-1, j-1 and k-1 */
  int xm2, ym2, zm2;            /* Sparse cell at i-2, j-2 and k-2 */
  int xp1, yp1, zp1;            /* Sparse cell at i+1, j+1 and k+1 */
  double face;                  /* Concentration value at the faces */
  double grad;                  /* Concentration difference across a face */
  double curv;                  /* Concentration curvature across a face */
  double cx, cy, cz;            /* Courant number in the x,y,z directions */

  /*-----------------------------------------------------------------*/
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  memset(Fy, 0, window->n3_t * sizeof(double));
  memset(Fz, 0, window->n3_t * sizeof(double));
  /* Horizontal fluxes */
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    xm2 = window->xm1[xm1];
    ym2 = window->ym1[ym1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the left face.  */
    /* The face value is formulated for a non-uniform grid by */
    /* linearly interpolating between cell centers and evaluating on */
    /* the face.  */
    cx = wincon->w6[c];
    face = 0.5 * (tr[c] + tr[xm1]);
    grad = tr[c] - tr[xm1];
    curv = 0.0;
    if (cx > 0.0)
      curv = tr[c] - 2.0 * tr[xm1] + tr[xm2];
    else if (cx < 0.0)
      curv = tr[xp1] - 2.0 * tr[c] + tr[xm1];
    Fx[c] = face - 0.5 * cx * grad - (1.0 - cx * cx) * curv / 6.0;

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the back face */
    cy = wincon->w7[c];
    face = 0.5 * (tr[c] + tr[ym1]);
    grad = tr[c] - tr[ym1];
    curv = 0.0;
    if (cy > 0.0)
      curv = tr[c] - 2.0 * tr[ym1] + tr[ym2];
    else if (cy < 0.0)
      curv = tr[yp1] - 2.0 * tr[c] + tr[ym1];
    Fy[c] = face - 0.5 * cy * grad - (1.0 - cy * cy) * curv / 6.0;
  }

  /*-----------------------------------------------------------------*/
  /* Vertical fluxes */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    zp1 = window->zp1[c];
    zm1 = window->zm1[c];
    zm2 = window->zm1[zm1];

    /*---------------------------------------------------------------*/
    /* Get the concentration values from QUICKEST at the lower face */
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
void van_leer(geometry_t *window, /* Window structure */
              window_t *windat, /* Window data structure */
              win_priv_t *wincon, /* Window geometry / constants */
              double *tr,       /* Tracer array */
              double *Fx,       /* x tracer flux */
              double *Fy,       /* y tracer flux */
              double *Fz        /* z tracer flux */
  )
{

  /* Horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  /* x direction horizontal fluxes */
  /* Wet cell faces in this window */
  van_leer_do(Fx, tr, windat->u1, wincon->w6, 1, window->v3_e1,
              window->w3_e1, window->xp1, window->xm1);
  /* Auxiliary cell faces associated with this window */
  van_leer_do(Fx, tr, windat->u1, wincon->w6, window->v3_t + 1,
              window->a3_t, window->w3_t, window->xp1, window->xm1);

  /* y direction horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fy, 0, window->n3_t * sizeof(double));
  van_leer_do(Fy, tr, windat->u2, wincon->w7, 1, window->v3_e2,
              window->w3_e2, window->yp1, window->ym1);
  /* Auxiliary cell faces associated with this window */
  van_leer_do(Fy, tr, windat->u2, wincon->w7, window->v3_t + 1,
              window->a3_t, window->w3_t, window->yp1, window->ym1);

  /* Vertical fluxes */
  memset(Fz, 0, window->n3_t * sizeof(double));
  van_leer_do(Fz, tr, windat->w, wincon->w8, 1, wincon->vc, wincon->s1,
              window->zp1, window->zm1);

}

/* END van_leer()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using      */
/* VanLeer's method.                                                 */
/*-------------------------------------------------------------------*/
void van_leer_do(double *F,     /* Flux array */
                 double *tr,    /* df_variable_t to advect */
                 double *vel,   /* Velocity at the cell center/face */
                 double *cn,    /* Courant number for this dimension */
                 int ss,        /* Start sparse coordinate */
                 int se,        /* End sparse coordinate */
                 int *sdo,      /* Points to process array */
                 int *fmap,     /* Forward map */
                 int *bmap      /* Backward map */
  )
{
  int c, cc;
  int bm1, bm2;
  double df;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    bm1 = bmap[c];
    bm2 = bmap[bm1];

    if (vel[c] > 0.0) {
      df = (tr[bm1] - tr[bm2]) * (tr[c] - tr[bm1]);
      if (df > 0.0)
        df = 2.0 * df / (tr[c] - tr[bm2]);
      else
        df = 0.0;
      F[c] = tr[bm1] + 0.5 * (1.0 - cn[c]) * df;
    } else if (vel[c] < 0.0) {
      df = (tr[c] - tr[bm1]) * (tr[fmap[c]] - tr[c]);

      if (df > 0.0)
        df = 2.0 * df / (tr[fmap[c]] - tr[bm1]);
      else
        df = 0.0;
      F[c] = tr[c] - 0.5 * (1.0 + cn[c]) * df;
    }
  }
}

/* END van_leer_do()                                                 */
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
  double *crfxf, *crfyf, *crfzf;   /* Fractional factors of total trajectories for faces */
  double *crfxc, *crfyc, *crfzc;   /* Fractional factors of total trajectories for cell centres */
  int *clxf, *clyf, *clzf;         /* Cell counter at head of trajectory, for faces */
  int *clxc, *clyc, *clzc;         /* Cell counter at head of trajectory, for cell centres */
  double *nu, *nv, *nw;            /* Cell-centered velocities */
  double *u1, *u2, *w;             /* Velocities at faces */
  int c, cc, cp1;

  /* Set pointers. */
  clxf = wincon->clxf;
  clyf = wincon->clyf;
  clzf = wincon->clzf;
  clxc = wincon->clxc;
  clyc = wincon->clyc;
  clzc = wincon->clzc;
  crfxf = wincon->crfxf;
  crfyf = wincon->crfyf;
  crfzf = wincon->crfzf;
  crfxc = wincon->crfxc;
  crfyc = wincon->crfyc;
  crfzc = wincon->crfzc;
  nu = wincon->nu;
  nv = wincon->nv;
  nw = wincon->nw;
  u1 = windat->u1;
  u2 = windat->u2;
  w = windat->w;

  /* Initialise */
  memset(crfxf, 0, window->sgsiz * sizeof(double));
  memset(crfyf, 0, window->sgsiz * sizeof(double));
  memset(crfzf, 0, window->sgsiz * sizeof(double));
  memset(crfxc, 0, window->sgsiz * sizeof(double));
  memset(crfyc, 0, window->sgsiz * sizeof(double));
  memset(crfzc, 0, window->sgsiz * sizeof(double));
  memset(clxf, 0, window->sgsiz * sizeof(int));
  memset(clyf, 0, window->sgsiz * sizeof(int));
  memset(clzf, 0, window->sgsiz * sizeof(int));
  memset(clxc, 0, window->sgsiz * sizeof(int));
  memset(clyc, 0, window->sgsiz * sizeof(int));
  memset(clzc, 0, window->sgsiz * sizeof(int));

  /*-----------------------------------------------------------------*/
  /* Calculate the cell-averaged velocities in the e1 direction, and */
  /* identify the cell from which to start the backward trajectory   */
  /* (which is dependent on the velocity sense).                     */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cp1 = window->xp1[c];
    clxc[c] =  (nu[c] < 0.0) ? window->xp1[c] : window->w3_t[cc];
  }

  /* Calculate the cell-averaged velocities in the e2 direction, and */
  /* identify the cell from which to start the backward trajectory   */
  /* (which is dependent on the velocity sense).                     */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cp1 = window->yp1[c];
    clyc[c] =  (nv[c] < 0.0) ? window->yp1[c] : window->w3_t[cc];
  }

  /* Calculate the cell-averaged velocities in the z direction, and  */
  /* identify the cell from which to start the backward trajectory   */
  /* (which is dependent on the velocity sense).                     */
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    cp1 = window->zp1[c];
    clzc[c] = (nw[c] < 0.0) ? window->zp1[c] : window->w3_t[cc];
  }
  
  /* Get trajectory distances and origins for the transverse terms  */
  /* based on the cell-centered velocities.                         */

  ff_sl_do_hor(nu, dt, window->h1au1, window->m2d, 1, window->b3_t,
	  window->w3_t, window->xp1, window->xm1, clxc, crfxc);

  ff_sl_do_hor(nv, dt, window->h2au2, window->m2d, 1, window->b3_t,
	  window->w3_t, window->yp1, window->ym1, clyc, crfyc);

  ff_sl_do_vert(nw, dt, wincon->w9, 1, window->b3_t,
       window->w3_t, window->zp1, window->zm1, clzc, crfzc);

  /* Get trajectory distance and origin for advection using */
  /* cell face velocities for each direction.               */

  /* e1 direction */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    clxf[c] = (u1[c] > 0.0) ? window->xm1[c] : window->w3_e1[cc];
  }

  ff_sl_do_hor(u1, dt, window->h1acell, window->m2d, 1, window->b3_e1,
     window->w3_e1, window->xp1, window->xm1, clxf, crfxf);

  /* e2 direction */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    clyf[c] =  (u2[c] > 0.0) ? window->ym1[c] : window->w3_e2[cc];
  }
  
  ff_sl_do_hor(u2, dt, window->h2acell, window->m2d, 1, window->b3_e2,
     window->w3_e2, window->yp1, window->ym1, clyf, crfyf);

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
/*-------------------------------------------------------------------*/
/* Routine to calculate flux-form semi-lagrangian advective          */
/* fluxes (following Lin & Rood (1996).                              */
/*-------------------------------------------------------------------*/
void ff_sl(geometry_t *window, /* Window structure */
	   window_t *windat,   /* Window data structure */
	   win_priv_t *wincon, /* Window geometry / constants */
	   double *tr,         /* Tracer array */
	   double *Fx,         /* x tracer flux */
	   double *Fy,         /* y tracer flux */
	   double *Fz,         /* z tracer flux */
	   double dt           /* sub-time step */
  )
{
  double *crfxf, *crfyf, *crfzf;   /* Fractional factors of total trajectories for faces */
  double *crfxc, *crfyc        ;   /* Fractional factors of total trajectories for cell centres */
  int *clxf, *clyf, *clzf;         /* Cell counter at head of trajectory, for faces */
  int *clxc, *clyc;                /* Cell counter at head of trajectory, for cell centres */
  double *tr_mod;                  /* Modified tracer concentrations */
  double *nu, *nv;                 /* Cell-centered velocity */
  double *u1, *u2, *w;             /* Velocities at faces */
  double *u1vm, *u2vm;             /* Volume fluxes */
  double dx, dy, dz;
  double ltraj, lrem;
  int c, cc, cp, cp1, up1, c2;

  /* Set pointers. */
  clxf = wincon->clxf;
  clyf = wincon->clyf;
  clzf = wincon->clzf;
  clxc = wincon->clxc;
  clyc = wincon->clyc;
  crfxf = wincon->crfxf;
  crfyf = wincon->crfyf;
  crfzf = wincon->crfzf;
  crfxc = wincon->crfxc;
  crfyc = wincon->crfyc;
  tr_mod = wincon->tr_mod;
  nu = wincon->nu;
  nv = wincon->nv;
  u1 = windat->u1;
  u2 = windat->u2;
  w = windat->w;
  u1vm = windat->u1vm;
  u2vm = windat->u2vm;

  /* Initialise the arrays.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  memset(Fy, 0, window->n3_t * sizeof(double));
  memset(Fz, 0, window->n3_t * sizeof(double));

  /* Define new tracers with transverse terms */
  memcpy(tr_mod, tr, window->sgsiz * sizeof(double));

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    up1 = (nv[c] < 0.0) ? clyc[c] : window->ym1[clyc[c]];
    cp = (nv[c] < 0.0) ? window->ym1[up1] : clyc[c];
    c2 = (nv[c] < 0.0) ? window->m2d[up1] : window->m2d[cp];
    tr_mod[c] += tr[cp];
    if (crfyc[c] > 0.) {
      tr_mod[c] += (crfyc[c] * (tr[up1] - tr[cp]));
    }
    tr_mod[c] *= 0.5;
  }

  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fx, tr_mod, u1, crfxf, clxf, 1, window->b3_e1,
		 window->w3_e1, window->xp1, window->xm1);

  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    if (u1[c] == 0.0) {
      Fx[c] = 0.0;
      continue;
    }
    cp = clxf[c];
    c2 = window->m2d[cp];
    Fx[c] *= (crfxf[c] * window->h1acell[c2]) ;

    /* Integrate over the "integer" component of the trajectory */
    cp = (u1[c] > 0.0) ?  window->xm1[c] : window->w3_t[cc];
    c2 = window->m2d[cp];
    dx = window->h1acell[c2];
    ltraj = fabs(u1[c] * dt);
    lrem = ltraj;
    while (dx < lrem) {
      Fx[c] += (dx * tr_mod[cp]);
      lrem -= dx;
      cp = (u1[c] < 0) ? window->xp1[cp] : window->xm1[cp];
      c2 = window->m2d[cp];
      if (window->h1acell[c2] > 0)
	dx = window->h1acell[c2];
    }

    /* Multiply by the volume flux */
    Fx[c] *= (u1vm[c] / ltraj);
  }

  /* y direction horizontal fluxes */

  memcpy(tr_mod, tr, window->sgsiz * sizeof(double));

  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    up1 = (nu[c] < 0.0) ? clxc[c] : window->xm1[clxc[c]];
    cp = (nu[c] < 0.0) ? window->xm1[up1] : clxc[c];
    c2 = (nu[c] < 0.0) ? window->m2d[up1] : window->m2d[cp];
    tr_mod[c] += tr[cp];
    if (crfxc[c] > 0.0) {
      tr_mod[c] += (crfxc[c] * (tr[up1] - tr[cp]));
    }
    tr_mod[c] *= 0.5;
  }
		       
  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fy, tr_mod, u2, crfyf, clyf, 1, window->b3_e2,
     window->w3_e2, window->yp1, window->ym1);

  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    if (u2[c] == 0.0) {
      Fy[c] = 0.0;
      continue;
    }
    cp = clyf[c];
    c2 = window->m2d[cp];
    Fy[c] *= (crfyf[c] * window->h2acell[c2]);

    /* Integrate over the "integer" component of the trajectory */
    cp = (u2[c] > 0.0) ?  window->ym1[c] : window->w3_t[cc];
    c2 = window->m2d[cp];
    dy = window->h2acell[c2];
    ltraj = fabs(u2[c] * dt);
    lrem = ltraj;
    while (dy < lrem) {
      Fy[c] += (dy * tr_mod[cp]);
      lrem -= dy;
      cp = (u2[c] < 0) ? window->yp1[cp] : window->ym1[cp];
      c2 = window->m2d[cp];
      if (window->h2acell[c2] > 0)
	dy = window->h2acell[c2];
    }

    /* Multiply by the volume flux */
    Fy[c] *= (u2vm[c] / ltraj);
  }

  /* Vertical fluxes */

  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fz, tr, w, crfzf, clzf, 1, wincon->vc,
     wincon->s1, window->zp1, window->zm1);

  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    if (w[c] == 0.0) {
      Fz[c] = 0.0;
      continue;
    }
    cp = clzf[c];
    Fz[c] *= (crfzf[c] * wincon->dz[cp]);

    /* Integrate over the "integer" component of the trajectory */
    cp = (w[c] > 0.0) ? window->zm1[c] : wincon->s1[cc];
    dz = wincon->dz[cp];
    lrem = fabs(w[c] * dt);
    while (dz < lrem) {
      Fz[c] += (dz * tr[cp]);
      lrem -= dz;
      cp = (w[c] < 0.0) ? window->zp1[cp] : window->zm1[cp];
      if (wincon->dz[cp] > 0)
	dz = wincon->dz[cp];
    }

    /* Ensure Fz has correct sign */
    if (w[c] < 0)
      Fz[c] *= -1.0;
  }
}

/* END ff_sl()                                                        */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate flux-form advective fluxes for the full 3D   */
/* semi-lagrangian advective scheme, based on Leonard et al. (1996). */
/*-------------------------------------------------------------------*/
void ff_sl3d(geometry_t *window, /* Window structure */
	     window_t *windat,   /* Window data structure */
	     win_priv_t *wincon, /* Window geometry / constants */
	     double *tr,         /* Tracer array */
	     double *Fx,         /* x tracer flux */
	     double *Fy,         /* y tracer flux */
	     double *Fz,         /* z tracer flux */
	     double *tr_mod,     /* Modified tracer concentrations */
	     double *tr_mod_x,   /* Modified tracer concentrations */
	     double *tr_mod_y,   /* Modified tracer concentrations */
	     double *tr_mod_z,   /* Modified tracer concentrations */
	     double dt           /* sub-time step */
  )
{
  double *u1vm, *u2vm;             /* Volume fluxes */
  double *crfxf, *crfyf, *crfzf;   /* Fractional factors of total trajectories for faces */
  double *crfxc, *crfyc, *crfzc;   /* Fractional factors of total trajectories for cell centres */
  int *clxf, *clyf, *clzf;         /* Cell counter at head of trajectory, for faces */
  int *clxc, *clyc, *clzc;         /* Cell counter at head of trajectory, for cell centres */
  double *nu, *nv, *nw;            /* Cell-centered velocities */
  double *u1, *u2, *w;             /* Velocities at faces */
  double dx, dy, dz;
  double ltraj, lrem;
  double onesixth = 1.0 / 6.0;
  int c, cc, cp, cp1, up1, c2;

  /* Set pointers. */
  clxf = wincon->clxf;
  clyf = wincon->clyf;
  clzf = wincon->clzf;
  clxc = wincon->clxc;
  clyc = wincon->clyc;
  clzc = wincon->clzc;
  crfxf = wincon->crfxf;
  crfyf = wincon->crfyf;
  crfzf = wincon->crfzf;
  crfxc = wincon->crfxc;
  crfyc = wincon->crfyc;
  crfzc = wincon->crfzc;
  nu = wincon->nu;
  nv = wincon->nv;
  nw = wincon->nw;
  u1 = windat->u1;
  u2 = windat->u2;
  w = windat->w;
  u1vm = windat->u1flux3d;
  u2vm = windat->u2flux3d;

  /* Initialize the modified tracer array. */
  memset(tr_mod, 0, window->sgsiz * sizeof(double));

  /* Initialise the flux arrays.  */
  memset(Fx, 0, window->sgsiz * sizeof(double));
  memset(Fy, 0, window->sgsiz * sizeof(double));
  memset(Fz, 0, window->sgsiz * sizeof(double));

  /* Set no gradients if required */
  if (wincon->conserve & CONS_NGR) nogr(window, windat, wincon, tr); 

  /* Define a new tracer with x-direction transverse terms */
  memcpy(tr_mod_x, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    up1 = (nu[c] < 0.0) ? clxc[c] : window->xm1[clxc[c]];
    cp = (nu[c] < 0.0) ? window->xm1[up1] : clxc[c];
    tr_mod_x[c] = tr[cp];
    if (crfxc[c] > 0.) {
      tr_mod_x[c] += (crfxc[c] * (tr[up1] - tr[cp]));
    }
    /* Monotonicity constraint */
    tr_mod_x[c] = monoco(c, cp, up1, tr_mod_x, tr, 1.0);
  }
  
  /* Define a new tracer with y-direction transverse terms */
  memcpy(tr_mod_y, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    up1 = (nv[c] < 0.0) ? clyc[c] : window->ym1[clyc[c]];
    cp = (nv[c] < 0.0) ? window->ym1[up1] : clyc[c];
    tr_mod_y[c] = tr[cp];
    if (crfyc[c] > 0.) {
      tr_mod_y[c] += (crfyc[c] * (tr[up1] - tr[cp]));
    }
    /* Monotonicity constraint */
    tr_mod_y[c] = monoco(c, cp, up1, tr_mod_y, tr, 1.0);
  }

  /* Define a new tracer with z-direction transverse terms */
  memcpy(tr_mod_z, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    up1 =  (nw[c] < 0.0) ? clzc[c] : window->zm1[clzc[c]];
    cp = (nw[c] < 0.0) ?  window->zm1[up1] : clzc[c];
    tr_mod_z[c] = tr[cp];
    if (crfzc[c] > 0.) {
      tr_mod_z[c] += (crfzc[c] * (tr[up1] - tr[cp]));
    }
    /* Monotonicity constraint */
    tr_mod_z[c] = monoco(c, cp, up1, tr_mod_z, tr, 1.0);
  }

  /* Calculate the cross-coupling advective-form updates to the tracer field and
     store the modified tracer fields. */
  /* Do the x-direction terms first. */
  memcpy(tr_mod, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    tr_mod[c] += (tr[c] + tr_mod_y[c] + tr_mod_z[c]);
    up1 = (nw[c] < 0.0) ? clzc[c] : window->zm1[clzc[c]];
    cp = (nw[c] < 0.0) ? window->zm1[up1] : clzc[c];
    tr_mod[c] += (tr_mod_y[cp] + (crfzc[c] * (tr_mod_y[up1] - tr_mod_y[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_y, 5.0);*/
    up1 = (nv[c] < 0.0) ? clyc[c] : window->ym1[clyc[c]];
    cp = (nv[c] < 0.0) ? window->ym1[up1] : clyc[c];
    tr_mod[c] += (tr_mod_z[cp] + (crfyc[c] * (tr_mod_z[up1] - tr_mod_z[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_z, 6.0);*/
    tr_mod[c] *= onesixth;
  }

  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fx, tr_mod, u1, crfxf, clxf, 1, window->b3_e1,
     window->w3_e1, window->xp1, window->xm1);
  /*
  ff_sl_order1(window, Fx, tr_mod, u1, crfxf, clxf, 1, window->b3_e1,
	       window->w3_e1, window->xp1, window->xm1);
  */

  for (cc = 1; cc <= window->b3_e1; cc++) {
    c = window->w3_e1[cc];
    if (u1[c] == 0.0) {
      Fx[c] = 0.0;
      continue;
    }
    cp = clxf[c];
    c2 = window->m2d[cp];

    Fx[c] *= (crfxf[c] * window->h1acell[c2]);
    /* Integrate over the "integer" component of the trajectory */
    cp = (u1[c] > 0.0) ?  window->xm1[c] : window->w3_e1[cc];
    c2 = window->m2d[cp];
    dx = window->h1acell[c2];
    ltraj = fabs(u1[c] * dt);
    lrem = ltraj;
    while (dx < lrem) {
      Fx[c] += (dx * tr_mod[cp]);
      lrem -= dx;
      cp = (u1[c] < 0) ? window->xp1[cp] : window->xm1[cp];
      c2 = window->m2d[cp];
      if (window->h1acell[c2] > 0)
	dx = window->h1acell[c2];
    }

    /* Multiply by the volume flux */
    Fx[c] *= (u1vm[c] / ltraj);
  }

  /* Do the y-direction terms next. */
  memcpy(tr_mod, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    tr_mod[c] += (tr[c] + tr_mod_x[c] + tr_mod_z[c]);
    up1 = (nw[c] < 0.0) ? clzc[c] : window->zm1[clzc[c]];
    cp = (nw[c] < 0.0) ? window->zm1[up1] : clzc[c];
    tr_mod[c] += (tr_mod_x[cp] + (crfzc[c] * (tr_mod_x[up1] - tr_mod_x[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_x, 5.0);*/
    up1 = (nu[c] < 0.0) ? clxc[c] : window->xm1[clxc[c]];
    cp = (nu[c] < 0.0) ? window->xm1[up1] : clxc[c];
    tr_mod[c] += (tr_mod_z[cp] + (crfxc[c] * (tr_mod_z[up1] - tr_mod_z[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_z, 6.0);*/
    tr_mod[c] *= onesixth;
  }

  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fy, tr_mod, u2, crfyf, clyf, 1, window->b3_e2,
     window->w3_e2, window->yp1, window->ym1);
  /*
  ff_sl_order1(window, Fy, tr_mod, u2, crfyf, clyf, 1, window->b3_e2,
	       window->w3_e2, window->yp1, window->ym1);
  */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    c = window->w3_e2[cc];
    if (u2[c] == 0.0) {
      Fy[c] = 0.0;
      continue;
    }
    cp = clyf[c];
    c2 = window->m2d[cp];

    Fy[c] *= (crfyf[c] * window->h2acell[c2]);
    /* Integrate over the "integer" component of the trajectory */
    cp = (u2[c] > 0.0) ? window->ym1[c] : window->w3_e2[cc];
    c2 = window->m2d[cp];
    dy = window->h2acell[c2];
    ltraj = fabs(u2[c] * dt);
    lrem = ltraj;
    while (dy < lrem) {
      Fy[c] += (dy * tr_mod[cp]);
      lrem -= dy;
      cp = (u2[c] < 0) ? window->yp1[cp] : window->ym1[cp];
      c2 = window->m2d[cp];
      if (window->h2acell[c2] > 0)
	dy = window->h2acell[c2];
    }

    /* Multiply by the volume flux */
    Fy[c] *= (u2vm[c] / ltraj);
  }

  /* Do the z-direction terms next. */
  memcpy(tr_mod, tr, window->sgsiz * sizeof(double));
  for (cc = 1; cc <= window->b3_t; cc++) {
    c = window->w3_t[cc];
    tr_mod[c] += (tr[c] + tr_mod_x[c] + tr_mod_y[c]);
    up1 = (nu[c] < 0.0) ? clxc[c] : window->xm1[clxc[c]];
    cp = (nu[c] < 0.0) ? window->xm1[up1] : clxc[c];
    tr_mod[c] += (tr_mod_y[cp] + (crfxc[c] * (tr_mod_y[up1] - tr_mod_y[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_y, 5.0);*/
    up1 = (nv[c] < 0.0) ? clyc[c] : window->ym1[clyc[c]];
    cp = (nv[c] < 0.0) ? window->ym1[up1] : clyc[c];
    tr_mod[c] += (tr_mod_x[cp] + (crfyc[c] * (tr_mod_x[up1] - tr_mod_x[cp])));
    /*tr_mod[c] = monoco(c, cp, up1, tr_mod, tr_mod_x, 6.0);*/
    tr_mod[c] *= onesixth;
  }

  /* Calculate the fractional part of the volume flux, using the Van Leer 
     algorithm to get the tracer value on the face. */
  ff_sl_van_leer(Fz, tr_mod, w, crfzf, clzf, 1, wincon->vc,
     wincon->s1, window->zp1, window->zm1);
  /*
  ff_sl_order1(window, Fz, tr_mod, w, crfzf, clzf, 1, wincon->vc,
     wincon->s1, window->zp1, window->zm1);
  */
  for (cc = 1; cc <= wincon->vc; cc++) {
    c = wincon->s1[cc];
    if (w[c] == 0.0) {
      Fz[c] = 0.0;
      continue;
    }
    cp = clzf[c];

    Fz[c] *= (crfzf[c] * wincon->dz[cp]);
    /* Integrate over the "integer" component of the trajectory */
    cp = (w[c] > 0.0) ? window->zm1[c] : wincon->s1[cc];
    dz = fabs(wincon->dz[cp]);
    lrem = fabs(w[c] * dt);
    while (dz < lrem) {
      Fz[c] += (dz * tr_mod[cp]);
      lrem -= dz;
      cp = (w[c] < 0.0) ? window->zp1[cp] : window->zm1[cp];
      if (wincon->dz[cp] > 0)
	dz = wincon->dz[cp];
    }

    /* Ensure Fz has correct sign */
    if (w[c] < 0)
      Fz[c] *= -1;
  }
}

/* END ff_sl3d()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Monotonicity constraint for transverse terms in FFSL. Do not      */
/* a linear interpolation between two values be greater or less than */
/* either of those values.                                           */
/*-------------------------------------------------------------------*/
double monoco(int c, int cm, int cp, double *tr, double *trmod, double sc) {
  double v, mn, mx;

  mn = min(trmod[cm], trmod[cp]);
  mx = max(trmod[cm], trmod[cp]);
  v = max(mn * sc, tr[c]);
  v = min(mx * sc, v);

  return(v);
}

/* END monoco()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Sets no-gradient conditions for a tracer.                         */
/*-------------------------------------------------------------------*/
void nogr(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr) 
{
  int c, c2, cc, zp1;            /* Counters */

  /* Set the boundary conditions (no flux) */
  for (cc = 1; cc <= window->nbpt; cc++) {
    c = window->bpt[cc];
    c2 = window->bin[cc];
    tr[c] = tr[c2];
  }
  /* No gradient at the bottom */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->bot_t[cc];
    tr[window->zm1[c]] = tr[c];
  }
  /* No-gradient above the surface */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = c2 = wincon->s1[cc];
    zp1 = window->zm1[c];
    while (c != zp1) {
      c = zp1;
      tr[c] = tr[c2];
      zp1 = window->zp1[c];
    }
  }
}

/* END nogr()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* First order routine to track the backwards trajectory to the      */
/* source, and note the start sparse location and fractional         */
/* component of the trajectory.                                      */
/*-------------------------------------------------------------------*/
void ff_sl_do_hor(double *vel,   /* Velocity at the cell face */
                  double dt,     /* Tracer advection tine step */
                  double *h,     /* cell width/height at face/centre */
                  int *m2d,      /* 2D sparse location corresponding to c */
                  int ss,        /* Start sparse coordinate */
                  int se,        /* End sparse coordinate */
                  int *sdo,      /* Points to process array */
                  int *fmap,     /* Forward map */
                  int *bmap,     /* Backward map */
                  int *cl,       /* Index of trajectory source cell */
                  double *crf    /* Fractional component of trajectory */
  )
{
  int c, cc, c2, cp;
  double ltraj, lrem, dx;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    if (vel[c] == 0.0)
      continue;
    ltraj = fabs(vel[c]) * dt;
    lrem = ltraj;
    cp = cl[c];
    c2 = m2d[cp];
    dx = h[c2];

    while (dx < lrem) {
      lrem  -= dx;
      cp = (vel[c] < 0.0) ? fmap[cp] : bmap[cp];
      c2 = m2d[cp];
      if (h[c2] > 0)
        dx = h[c2];
    }

    crf[c] = max(min(lrem / dx, 1.0), 0.0);
    cl[c] = cp;
  }
}

/* END ff_sl_do_hor()                                                 */
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
/* Routine to get the tracer concentration on a cell face using a    */
/* first order upstream method. Used with the FFSL advection scheme. */
/*-------------------------------------------------------------------*/
void ff_sl_order1(geometry_t *window,
		  double *F,  /* Flux array */
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
  int bm1;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    if (crf[c] == 0.0)
      continue;
    cp = cl[c];
    if (vel[c] > 0.0)
      cp = fmap[cp];
    bm1 = bmap[cp];

    if (vel[c] > 0.0) {
      F[c] = tr[bm1];
    } else if (vel[c] < 0.0) {
      F[c] = tr[cp];
    }
  }
}

/* END ff_sl_order1()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to calculate the 2nd order upwind advective fluxes. This  */
/* scheme is stable for |cn| <= 2 (Leonard, 1994, Comput. Methods    */
/* Appl. Mech. Engrg.)                                               */
/*-------------------------------------------------------------------*/
void order2_upwind(geometry_t *window, /* Window structure */
		   window_t *windat, /* Window data structure */
		   win_priv_t *wincon, /* Window geometry / constants */
		   double *tr,       /* Tracer array */
		   double *Fx,       /* x tracer flux */
		   double *Fy,       /* y tracer flux */
		   double *Fz        /* z tracer flux */
  )
{

  /* Horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fx, 0, window->n3_t * sizeof(double));
  /* x direction horizontal fluxes */
  /* Wet cell faces in this window */
  order2_upwind_do(Fx, tr, windat->u1, wincon->w6, 1, window->v3_e1,
		   window->w3_e1, window->xp1, window->xm1);
  /* Auxiliary cell faces associated with this window */
  order2_upwind_do(Fx, tr, windat->u1, wincon->w6, window->v3_t + 1,
		   window->a3_t, window->w3_t, window->xp1, window->xm1);

  /* y direction horizontal fluxes */
  /* Initialise the flux array.  */
  memset(Fy, 0, window->n3_t * sizeof(double));
  order2_upwind_do(Fy, tr, windat->u2, wincon->w7, 1, window->v3_e2,
		   window->w3_e2, window->yp1, window->ym1);
  /* Auxiliary cell faces associated with this window */
  order2_upwind_do(Fy, tr, windat->u2, wincon->w7, window->v3_t + 1,
		   window->a3_t, window->w3_t, window->yp1, window->ym1);

  /* Vertical fluxes */
  memset(Fz, 0, window->n3_t * sizeof(double));
  order2_upwind_do(Fz, tr, windat->w, wincon->w8, 1, wincon->vc, wincon->s1,
		   window->zp1, window->zm1);

}

/* END order2_upwind()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the tracer concentration on a cell face using      */
/* the second order upwind method.                                   */
/*-------------------------------------------------------------------*/
void order2_upwind_do(double *F,     /* Flux array */
		      double *tr,    /* Variable to advect */
		      double *vel,   /* Velocity at the cell center/face */
		      double *cn,    /* Courant number for this dimension */
		      int ss,        /* Start sparse coordinate */
		      int se,        /* End sparse coordinate */
		      int *sdo,      /* Points to process array */
		      int *fmap,     /* Forward map */
		      int *bmap      /* Backward map */
  )
{
  int c, cc;
  int bm1, bm2;
  double df;

  for (cc = ss; cc <= se; cc++) {
    c = sdo[cc];
    bm1 = bmap[c];
    bm2 = bmap[bm1];

    if (vel[c] > 0.0) {
      F[c] = 0.5 * (3.0 - cn[c]) * tr[bm1] -
	0.5 * (1.0 - cn[c]) * tr[bm2];
    } else if (vel[c] < 0.0) {
      F[c] = 0.5 * (3.0 - cn[c]) * tr[c] -
	0.5 * (1.0 - cn[c]) * tr[fmap[c]];
    } else
      F[c] = 0.0;
  }
}

/* End order2_upwind_do()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Calculates updated tracer concentration of surface cells          */
/*-------------------------------------------------------------------*/
double surf_conc(geometry_t *window,  /* Window geometry             */
                 window_t *windat,    /* Window data structure       */
                 win_priv_t *wincon,  /* Window constants            */
                 double *tr,    /* Tracer array                      */
                 double subeta, /* Surface elevation at sub-timestep */
                 int c,         /* 3D coord of lowest surface cell   */
                 int c2,        /* 2D sparse coordinate              */
                 int cc,        /* 2D sparse counter                 */
                 double fcbot,  /* Tracer flux through the bottom    */
                 double *dtracer /* Advective terms x,y directions   */
  )
                         /* [tracer][m]^2.  */
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
  double osubeta;        /* Old surface elevation                    */
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
  osubeta = wincon->d1[c2];
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
    /* 
    if(window->wn == dw && c2 == dc) {
      hf += (windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
	     windat->u2flux3d[window->yp1[c]] -
	     windat->u2flux3d[c]) * dtu;
      dtr += dtracer[c];
      printf("rise %d %d %f : %f %f\n",c,window->s2k[c],surftr,hf,dtr/tr[c]);
    }
    */
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
    hf += (windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
	   windat->u2flux3d[window->yp1[c]] -
	   windat->u2flux3d[c]) * dtu;
    dtr += dtracer[c];
    cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
    printf("%f %d %d %e : %f %f : %f %f : %e %e\n",windat->t/86400,
	   c,window->s2k[c],cnt,hf,dtr/tr[c],vf,
	   fcbot/tr[c],surftr/svol,windat->dum3[c]);
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
	dtr = (windat->u1flux3d[window->xp1[c]] - windat->u1flux3d[c] +
	   windat->u2flux3d[window->yp1[c]] -
	   windat->u2flux3d[c]) * dtu;
	hf += dtr;
	cnt = window->cellarea[c2] * (osubeta - subeta) - (hf + vf);
	printf("to-top %d %d %e : %f %f : %e %e\n",window->s2k[c],c,cnt,dtr,dtracer[c],surftr/svol,windat->dum3[c]);
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
/* Routine to check for tracer range violations and print a message. */
/* Min / max violation for several special cases are clipped.        */
/*-------------------------------------------------------------------*/
void tr_bounds(geometry_t *window,  /* Window structure */
               window_t *windat,  /* Window data structure */
               win_priv_t *wincon /* Window geometry / constants */
  )
{
  int c, c3, cc;                /* Counters */
  int cs, cb;                   /* Surface and bottom coordinates */
  int co;                       /* Old surface coordinate */
  int n, nn;                    /* Tracer number */
  double surftr;                /* New surface tracer value */
  double eps = 1e-5;            /* Minimum value for clipping */
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

  /* In thin layers the difference in diffusive or advective fluxes */
  /* between opposite faces can be large (eg. if dzface(i) is thin */
  /* and dzface(i+1) is not), and dividing this by the (small) cell */
  /* volume can result in large rates of change of the tracer, which */
  /* can violate the maximum/minimum allowed concentration.  */
  /* Therefore, if the layer is thin and a max / min violation */
  /* occurs in the surface layer only, then set the surface */
  /* concentration to that of the layer below.  */
  for (cc = 1; cc <= wincon->vcs1; cc++) {
    c = c3 = wincon->s1[cc];
    co = wincon->i3[cc];
    cb = wincon->i1[cc];
    cs = window->m2d[c];
    /* Use eta-gridz here since dz=oldeta-gridz
       if (c < cb && windat->eta[cs] - window->gridz[co] < wincon->hmin) { */
    /* if (c < cb && (wincon->dz[co] < wincon->hmin ||
       windat->eta[cs] - window->gridz[wincon->i2[cc]] < wincon->hmin)) { */
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
      /* If a minimum is set to zero and the tracer undershoots */
      /* below this then set to zero to keep positive definite.  */
      if (min[n] <= eps && windat->tr_wc[n][c] < min[n])
        windat->tr_wc[n][c] = min[n];

      /* Print a warning for any other min / max violations */
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

  /* Set all layers above the surface to surface layer   */
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
void tr_set_surf(geometry_t *window,  /* Window structure            */
		 window_t *windat,    /* Window data structure       */
		 win_priv_t *wincon   /* Window geometry / constants */
  )
{
  int c, cs, cc, n, nn;

  /* Set all layers above the surface to surface layer               */
  if (!wincon->sigma) {
    for (cc = 1; cc <= wincon->vcs2; cc++) {
      c = cs = wincon->s2[cc];
      while (c != window->zp1[c]) {
	c = window->zp1[c];
	/*
	for (nn = 0; nn < wincon->ntbdy; nn++) {
	  n = wincon->tbdy[nn];
	*/
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
void check_bounds(char *text,                    /* Warning text     */
		  geometry_t *window,            /* Window geometry  */
		  window_t *windat,              /* Window data      */
		  win_priv_t *wincon             /* Window constants */
  )
{
  int c, cc;                    /* Counters      */
  int n, nn;                    /* Tracer number */
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
  int c, cc;                    /* Counters      */
  int n, nn;                    /* Tracer number */

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
  int c, c2, cc;                    /* Counters      */
  double *Fx = wincon->w2;
  double *Fy = wincon->w3;
  double *Fz = wincon->w4;
  double d1, f1, f2, f3;

  if (strcmp(trname, wincon->trname[n]) == 0) {
    double min = wincon->mintr[n];
    double max = wincon->maxtr[n];
    for (cc = 1; cc <= wincon->vc1; cc++) {
      c = wincon->s1[cc];
      c2 = window->m2d[c];
      d1 = dtu / (window->cellarea[c2] * wincon->dz[c]);
      f1 = (Fx[window->xp1[c]] - Fx[c]) * d1;
      f2 = (Fy[window->yp1[c]] - Fy[c]) * d1;
      f3 = (Fz[window->zp1[c]] - Fz[c]) / wincon->dz[c];
      if (fabs(f1) > max)
	hd_warn
	  ("%s : Tracer %s x flux divergence %e outside range : (%e - %e) (%e - %e) @ (%d %d %d) : %d at %f. dz = %f\n",
	   text, wincon->trname[n], f1, windat->tr_wc[n][c],
	   windat->tr_wc[n][window->xp1[c]],Fx[window->xp1[c]],Fx[c],
	   window->s2i[c], window->s2j[c], window->s2k[c], c, 
	   windat->t / 86400.0, wincon->dz[c]);
      if (fabs(f2) > max)
	hd_warn
	  ("%s : Tracer %s y flux divergence %e outside range : (%e - %e) (%e - %e) @ (%d %d %d) : %d at %f. dz = %f\n",
	   text, wincon->trname[n], f2, windat->tr_wc[n][c],
	   windat->tr_wc[n][window->yp1[c]],Fy[window->yp1[c]],Fy[c],
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
/* Invokes the ULTIMATE filter                                       */
/*-------------------------------------------------------------------*/
void ulimate_filter(double *cr, /* Courant number at the cell face */
                    double *Fx, /* Scalar flux at the cell face */
                    double *F,  /* Scalar value */
                    int *ctp,   /* Sparse cells to process */
                    int vs,     /* Start coordinate of ctp */
                    int ve,     /* End coordinate of ctp */
                    int *mapm,  /* Map to value at xm1 or ym1 or zm1 */
                    int *mapp   /* Map to value at xp1 or yp1 or zp1 */
  )
{
  double nc;                    /* Normalised center value */
  double nu;                    /* Normalised face value */
  double tc = 0.0;              /* tracer value at the cell center */
  double tu = 0.0;              /* Upstream tracer value */
  double td = 0.0;              /* Downstream tracer value */
  double diff;                  /* Difference of upstream and downstream */
  int c, cc, m1;                /* Sparse coordinates / counters */

  /* Loop through the cells to process */
  for (cc = vs; cc <= ve; cc++) {
    c = ctp[cc];
    if (cr[c] == 0.0) {
      continue;
    } else {
      /* Get the normalised tracer values */
      /* nc=nu=Fx[c]; */
      if (cr[c] > 0.0) {
        m1 = mapm[c];
        tc = F[m1];
        tu = F[mapm[m1]];
        td = F[c];
      } else {
        tc = F[c];
        tu = F[mapp[c]];
        td = F[mapm[c]];
      }

      diff = td - tu;
      nc = nu = 0.0;
      if (diff) {
        nc = (tc - tu) / diff;
        nu = (Fx[c] - tu) / diff;
      }

      /* Get the adjusted face values */
      if (nc >= 0.0 && nc <= 1.0) {
        if (nu > 1.0)
          Fx[c] = td;
        else if (nu < nc)
          Fx[c] = tu + nc * (td - tu);
        else if (nu > nc / fabs(cr[c]))
          Fx[c] = tu + nc * (td - tu) / fabs(cr[c]);
      } else
        Fx[c] = tu + nc * (td - tu);
    }
  }
}

/* END ulimate_filter()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to linearly interpolate between tracer values at time t   */
/* and t+dt for multi dt auxiliary cells.                            */
/*-------------------------------------------------------------------*/
void set_multidt_t(geometry_t *window,  /* Window structure */
                   window_t *windat,  /* Window data structure */
                   double *tr,  /* Tracer array */
                   double trem, /* Time remaining in the sub-step loop */
                   int tn       /* Tracer number */
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
void hor_diffuse(geometry_t *window,  /* Processing window           */
                 window_t *windat,    /* Window data structure       */
                 win_priv_t *wincon,  /* Window geometry / constants */
                 double *tr,          /* Tracer array                */
                 double *Fx,          /* x direction flux            */
                 double *Fy           /* y direction flux            */
  )
{
  int c, cc;                    /* Sparse indices                    */
  int c2, cp;                   /* 2D cell corresponding to 3D cell  */
  int xm1, xp1;                 /* i index at i-1                    */
  int ym1, yp1;                 /* j index at j-1                    */
  double csx, csy;              /* Cross sectional areas             */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays                                       */
  /* w5 = x face cross sectional area                                */
  /* w6 = y face cross sectional area                                */
  /* s1 = wet cells to process for tracers (calculated in advection) */

  /*-----------------------------------------------------------------*/
  /* For sigma, limit the horizontal diffusion to maintain           */
  /* monotinicity (Herzfeld, 1996, p136).                            */
  if (wincon->sigma) {
    double tm, dp, dm;
    for (cc = 1; cc <= window->a3_t; cc++) {
      c = window->w3_t[cc];
      c2 = window->m2d[c];
      xp1 = window->xp1[c];
      yp1 = window->yp1[c];
      xm1 = window->xm1[c];
      ym1 = window->ym1[c];

      cp = window->m2d[xp1];
      dp = tr[xp1] - tr[c];
      dm = tr[c] - tr[xm1];
      tm = max(dp, tr[c]);
      tm = max(tm, dm);
      csx = wincon->Ds[cp] * dp / window->h1au1[cp] -
        wincon->Ds[c2] * dm / window->h1au1[c2];
      csx = tm * wincon->Hn1[c2] / (windat->dt * csx);
      if (wincon->u1kh[c] > fabs(csx))
        wincon->u1kh[c] = fabs(csx);

      cp = window->m2d[yp1];
      dp = tr[yp1] - tr[c];
      dm = tr[c] - tr[ym1];
      tm = max(dp, tr[c]);
      tm = max(tm, dm);
      csy = wincon->Ds[cp] * dp / window->h2au2[cp] -
        wincon->Ds[c2] * dm / window->h2au2[c2];
      csy = tm * wincon->Hn2[c2] / (windat->dt * csy);
      if (wincon->u2kh[c] > fabs(csy))
        wincon->u2kh[c] = fabs(csy);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Calculate the diffusive fluxes                                  */
  for (cc = 1; cc <= window->a3_t; cc++) {
    c = window->w3_t[cc];
    c2 = window->m2d[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    /* Get the cross sectional area of the cell faces.               */
    csx = windat->dzu1[c] * window->h2au1[c2] * wincon->mdx[c2];
    csy = windat->dzu2[c] * window->h1au2[c2] * wincon->mdy[c2];
    Fx[c] -=
      (csx * wincon->u1kh[c] * (tr[c] - tr[xm1]) / window->h1au1[c2]);
    Fy[c] -=
      (csy * wincon->u2kh[c] * (tr[c] - tr[ym1]) / window->h2au2[c2]);
  }
}

/* END hor_diffuse()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Horizontal diffusion                                              */
/*-------------------------------------------------------------------*/
void hor_diffuse_2d(geometry_t *window, /* Processing window */
                    window_t *windat, /* Window data structure */
                    win_priv_t *wincon, /* Window geometry / constants */
                    double *tr, /* Tracer array */
                    double *Fx, /* x direction flux */
                    double *Fy  /* y direction flux */
  )
{
  int c, cc;                    /* Sparse indices */
  int xm1;                      /* i index at i-1 */
  int ym1;                      /* j index at j-1 */
  double csx, csy;              /* Cross sectional areas */

  /*-----------------------------------------------------------------*/
  /* Calculate the diffusive fluxes */
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->w2_t[cc];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];

    /* Get the cross sectional area of the cell faces.  */
    csx = windat->depth_e1[c] * window->h2au1[c];
    csy = windat->depth_e2[c] * window->h1au2[c];

    Fx[c] -=
      (csx * wincon->u1kh[c] * (tr[c] - tr[xm1]) / window->h1au1[c]);
    Fy[c] -=
      (csy * wincon->u2kh[c] * (tr[c] - tr[ym1]) / window->h2au2[c]);

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
  if (mode)
    AH = wincon->u1kh;
  else
    AH = wincon->u2kh;

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
    xp1 = window->xp1[c];
    yp1 = window->yp1[c];
    xm1 = window->xm1[c];
    ym1 = window->ym1[c];
    c1 =
      wincon->u1kh[c] * (tr[xp1] + tr[xm1] - 2.0 * tr[c]) /
      (window->h1au1[cs] * window->h1au1[cs]);
    c2 = wincon->u2kh[c] * (tr[yp1] + tr[ym1] - 2.0 * tr[c]) /
      (window->h2au1[cs] *  window->h2au1[cs]);
    tr[c] += windat->dt * (c1 + c2);
  }
}

/* END hor_diffuse_simple()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Vertical diffusion                                                */
/*-------------------------------------------------------------------*/
void vert_diffuse_3d(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int nn;                       /* Tracer index */
  int cc;                       /* Counter */
  int c, c2;                    /* 3D and 2D sparse coordinate */
  int zm1, zp1;                 /* Sparse coordinate at k-1 and k+1 */
  int cs, cb;                   /* Surface sparse coordinate */
  int *ctp = wincon->i2;        /* Old surface sparse coordinate */
  int *cbt = wincon->i5;        /* Bottom sparse coordinate */
  int *cth = wincon->i4;        /* Thin layer locations */
  double *topflux;              /* Surface flux */
  double dt = windat->dt;       /* Time step for the window */
  double dzdt;                  /* dz divided by dt */
  double *dzface = wincon->w7;  /* Cell thickness at the cell face */
  double *dzcell = wincon->w9;  /* Cell centered cell thickness */
  double *Kz = wincon->w10;     /* Vertical diffusivity */
  double *Cm1 = wincon->w1;     /* Constant for implicit calculation */
  double *C = wincon->w2;       /* Constant for implicit calculation */
  double *Cp1 = wincon->w3;     /* Constant for implicit calculation */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* d1 = flux out of bottom layer */
  /* d2 = flux out of top layer */
  /* d3 = scaling = 1.0 */
  /* s1 = wet cells to process for tracers (calculated in advection) */
  memset(wincon->d1, 0, window->sgsizS * sizeof(double));
  memset(wincon->d2, 0, window->sgsizS * sizeof(double));
  memset(cth, 0, window->sgsizS * sizeof(int));
  memcpy(dzcell, wincon->dz, window->sgsiz * sizeof(double));
  memcpy(Kz, windat->Kz, window->sgsiz * sizeof(double));

  /*-----------------------------------------------------------------*/
  /* Merge thin layers if required */
  if (wincon->thin_merge) {
    for (cc = 1; cc <= wincon->vcs1; cc++) {
      c = cth[cc] = ctp[cc];
      cb = wincon->i1[cc];
      if (c != cb && dzcell[c] < wincon->hmin) {
        cs = ctp[cc] = window->zm1[c];
        for (nn = 0; nn < wincon->ntdif_v; nn++) {
          int n = wincon->tdif_v[nn];
          double *tr = windat->tr_wc[n];
          tr[c] = tr[cs] = (tr[c] * dzcell[c] + tr[cs] * dzcell[cs]) /
            (dzcell[cs] + dzcell[c]);
        }
        dzcell[cs] += dzcell[c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Calculate constants that are independent of the variable to be */
  /* diffused in the implicit calculation (i.e. the elements of the */
  /* tridiagonal matrix). The rhs is set up in the routine */
  /* implicit_vdiff_tr().  */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    c = cs = ctp[cc];           /* Old surface sparse coordinate */
    c2 = window->m2d[c];        /* 2D sparse location corresponding to c */
    zm1 = window->zm1[c];

    /* Single layer case (i.e. surface lies in the bottom layer) */
    if (zm1 == window->zm1[zm1] || !dzcell[zm1])
      continue;

    /* Get the cell thicknesses. The sparse coordinate, c, lies in */
    /* the sediment at the end of the loop.  */
    while (c != zm1 && dzcell[c]) {
      dzface[c] = 0.5 * (dzcell[zm1] + dzcell[c]);
      /* SIGMA : Precondition mixing coefficient and the cell for */
      /* the sigma calculation.  */
      dzcell[c] *= wincon->Ds[c2];
      Kz[c] /= wincon->Ds[c2];
      c = zm1;
      zm1 = window->zm1[c];
    }

    cb = c = window->zp1[c];    /* Set cb to the bottom sparse coordinate */
    zp1 = window->zp1[c];       /* Layer above the bottom */
    cbt[cc] = cb;

    /* Set up tri-diagonal set of equations.  */
    /* Bottom layer.  */
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

    /* Surface layer */
    c = cs;
    dzdt = dzcell[c] / dt;
    Cm1[c] = -Kz[c] / dzface[c];
    Cp1[c] = 0.0;
    C[c] = dzdt - Cm1[c];
  }

  /*-----------------------------------------------------------------*/
  /* Tracer loop */
  for (nn = 0; nn < wincon->ntdif_v; nn++) {
    int n = wincon->tdif_v[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values */
    double *Splus = NULL;

    /* Set the surface boundary condition */
    topflux = wincon->d2;
    memset(topflux, 0, window->sgsizS * sizeof(double));
    if (n == windat->tno) {
      if (wincon->heatflux & (ADVANCED | INVERSE | NET_HEAT | COMP_HEAT | COMP_HEAT_MOM)) {
	/* Do swr data assimilation if required */
	if(wincon->trinfo_3d[n].flag & DO_SWR_INVERSE) {
	  swr_assimilation(window, windat, wincon, n);
	  topflux = windat->heatf;
	} else if (wincon->compatible & V1562) {
	  set_sbc(window, windat, wincon);
	  Splus = wincon->w8;
	  memset(Splus, 0, window->sgsiz * sizeof(double));
	  topflux = windat->heatf;
	} else {
	  Splus = wincon->w8;
	  if (wincon->nswreg) swr_params_event(window, windat, wincon, n);
	  memset(Splus, 0, window->sgsiz * sizeof(double));
	  calc_swr(window, windat, wincon, Splus, 0);
	  topflux = windat->heatf;
	}
      }
      if (wincon->heatflux & (SURF_RELAX|AVHRR|GHRSST))
        surf_relax(window, windat, wincon);
    }
    if (n == windat->sno) {
      /* Convert to correct units (e.g. Eqn 3.3, Introduction to */
      /* Physical Oceanography, G.L. Mellor (1996)).  */
      if (wincon->saltflux & (ADVANCED | BULK)) {
	/*
	for (cc = 1; cc <= vcs; cc++) {
	  c = ctp[cc];
	  cb = cbt[cc];
	  c2 = window->m2d[c];
	  if ((c != cb) && windat->eta[c2] - window->gridz[c] > wincon->hmin)
	    topflux[c2] = windat->nsfd[c2] * windat->sal[c];
	}
	*/
        for (cc = 1; cc <= window->b2_t; cc++) {
          c = window->w2_t[cc];
          topflux[c] = windat->nsfd[c] * windat->sal[c];
        }
      }
    }
    /* Generic surface fluxes */
    if (wincon->sflux[n] >= 0) {
      memcpy(topflux, windat->tr_wcS[wincon->sflux[n]], window->sgsizS * sizeof(double));
    }

    /* No surface exchanges for thin layers */
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
/* Routine to optimize swr_attn, swr_tran and swr_babs to minimise   */
/* error compared to an analysis product.                            */
/*-------------------------------------------------------------------*/
void swr_assimilation(geometry_t *window,  /* Window geometry        */
		      window_t *windat,    /* Window data            */
		      win_priv_t *wincon,  /* Window constants       */
		      int n                /* Temperature tracer no. */
		      )
{
  int tn, k, c, cc, cs, cb, c2, ks, kb;
  double *tr = windat->tr_wc[n];    /* Temperature tracer            */
  double *analysis;                 /* Analysis tracer               */
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
  double *anys = wincon->v2;    /* Analysis vector                   */
  double attn;                  /* swr attenuation                   */
  double tran;                  /* swr transmission                  */
  double babs;                  /* swr bottom absorption             */
  double rms, trms;             /* Global rms error, target rms      */
  double tf = 0.2;              /* Target rms scaling factor         */

  /* Get the analysis tracer                                         */
  for (tn = 0; tn < windat->ntr; tn++) {
    if (wincon->trinfo_3d[tn].flag & SWR_INVERSE)
      break;
  }
  analysis = windat->tr_wc[tn];

  /* Do the parameter optimisation                                   */
  /*
  rms = get_global_rms(tr, analysis, wincon->vc, wincon->s1);
  trms = tf * rms;
  */
  memset(Splus, 0, window->sgsiz * sizeof(double));
  /*
  while (rms < trms) {
    calc_swr(window, windat, wincon, Splus, 0);
    implicit_vdiff_tr(window, windat, wincon, tr, Kz, dzcell, dzface,
		      wincon->d1, topflux, ctp, cbt, wincon->vcs, Splus,
		      NULL, wincon->one, C, Cp1, Cm1);
  */
  for (cc = 1; cc <= wincon->vcs; cc++) {
    cs = c = ctp[cc];
    cb = cbt[cc];               /* Bottom sparse coordinate */
    c2 = window->m2d[cs];
    ks = window->s2k[cs];
    kb = window->s2k[cb];
    rms = get_profile_rms(tr, analysis, cs, window->zm1);
    trms = tf * rms;
    while (rms < trms) {
      /* Set swr parameters                                            */
      attn = windat->swr_attn[c2];
      tran = windat->swr_tran[c2];
      babs = windat->swr_babs[c2];
      /* Get temp and analysis vectors and initialise                */
      c = c2;
      for (k = 0; k < window->nz; k++) {
	temp[k] = 0.0;
	anys[k] = 0.0;
	Splus[c] = 0.0;
	c = window->zm1[c];
      }
      c = cs;
      for (k = ks; k >= kb; k--) {
	temp[k] = tr[c];
	anys[k] = analysis[c];
	c = window->zm1[c];
      }
      /* Compute swr source term and mix vertically                  */
      calc_swr(window, windat, wincon, Splus, cc);
      implicit_vdiff_at_cc(window, windat, wincon, temp, Kz, dzcell, dzface,
			   wincon->d1, topflux, ctp, cbt, cc, Splus,
			   NULL, wincon->one, C, Cp1, Cm1, windat->dt);

      /* Do the estimation est(temp, anys, attn, tran, babs) */
      /* Copy swr parameters to the arrays                           */
      windat->swr_attn[c2] = attn;
      windat->swr_tran[c2] = tran;
      windat->swr_babs[c2] = babs;
      /* Get the new metric                                          */
      rms = get_profile_rms(tr, analysis, cs, window->zm1);
    }
  /*rms = get_global_rms(tr, analysis, wincon->vc, ctp);*/
  }
  wincon->trinfo_3d[n].flag = NO_SWR_INVERSE;
}

/* END swr_assimilation()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Implicit vertical diffusion scheme                                */
/*-------------------------------------------------------------------*/
void implicit_vdiff_tr(geometry_t *window,  /* Processing window */
                       window_t *windat,  /* Window data structure */
                       win_priv_t *wincon,  /* Window geometry / constants
                                             */
                       double *var, /* Array of variable values to diffuse
                                     */
                       double *Kz,  /* Mixing coefficient values */
                       double *dzcell,  /* Cell thicknesses */
                       double *dzface,  /* Cell thicknesses at cell face */
                       double *fb,  /* Flux out of bottom */
                       double *ft,  /* Flux out of top */
                       int *ctp,  /* Cells to process */
                       int *cbt,  /* Bottom coordinate */
                       int vcs, /* Last index of surface cells */
                       double *Splus, /* Positive part of source term */
                       double *Sminus,  /* Negative part of source term */
                       double *scale, /* SIGMA : (depth) scaling for
                                         velocity */
                       double *C, double *Cp1, double *Cm1)
{
  int c, k;                     /* Sparse coordinate */
  int cs, ks;                   /* Surface sparse coordinate */
  int cb, kb;                   /* Bottom sparse coordinate */
  int zm1;                      /* Sparse cell below c */
  int cc;                       /* Sparse coordinate counter */
  int c2;                       /* 2D cell corresponding to 3D location */
  double dt = windat->dt;       /* Time step for the window */
  double dzdt;                  /* dz / dt */
  double div;                   /* Constant */

  /*-----------------------------------------------------------------*/
  /* Set pointers.  */
  /* Note: the 3D work arrays wincon->w# could be used for the */
  /* dummy arrays below, but execution speed is considerably faster */
  /* when work array access is sequential in memory, hence the */
  /* mapping to a contiguous vertical 1D work array wincon->v#.  */
  int *cth = wincon->i4;
  double *rhs = wincon->v1;
  double *sol = wincon->v2;
  double *ud = wincon->v3;
  double *B = wincon->v4;
  double *Bm1 = wincon->v5;
  double *Bp1 = wincon->v6;
  double *dz = wincon->v7;

  /* Loop therough the surface cells in this window */
  for (cc = 1; cc <= vcs; cc++) {

    cs = c = ctp[cc];           /* Set cs to the surface sparse coordinate
                                 */
    cb = cbt[cc];               /* Bottom sparse coordinate */
    c2 = window->m2d[c];        /* 2D sparse location corresponding to c */
    zm1 = window->zm1[c];
    ks = window->s2k[cs];
    kb = window->s2k[cb];

    /*---------------------------------------------------------------*/
    /* Single layer case (i.e. surface lies in the bottom layer) */
    if (zm1 == window->zm1[zm1] || !dzcell[zm1]) {
      var[cs] +=
        windat->dt * (fb[c2] - ft[c2]) / max(dzcell[cs], wincon->hmin);
      if (c != cth[cc]) {
        var[cth[cc]] = var[cs];
      }
      continue;
    }

    /* Map the sparse arrays to local coordinates.  */
    c = cs;
    for (k = ks; k >= kb; k--) {
      B[k] = C[c];
      Bm1[k] = Cm1[c];
      Bp1[k] = Cp1[c];
      dz[k] = dzcell[c];
      c = window->zm1[c];
    }

    /*---------------------------------------------------------------*/
    /* Set up the rhs for the system of equations.  */
    /* Bottom layer.  */
    c = cb;                     /* Set c to the bottom sparse coordinate */
    dzdt = dz[kb] / dt;
    rhs[kb] = dzdt * var[cb] + fb[c2];

    /* Mid-water layers */
    c = window->zp1[cb];        /* Layer above the bottom */
    for (k = kb + 1; k < ks; k++) {
      dzdt = dz[k] / dt;
      rhs[k] = dzdt * var[c];
      c = window->zp1[c];
    }

    /* Surface layer */
    dzdt = dz[ks] / dt;
    rhs[ks] = dzdt * var[cs] - ft[c2];

    /*---------------------------------------------------------------*/
    /* Add positive part of source terms, if specified */
    if (Splus) {
      zm1 = window->zm1[c];
      c = cb;
      for (k = kb; k <= ks; k++) {
        rhs[k] += dz[k] * Splus[c];
        c = window->zp1[c];
      }
    }

    /* Add negative part of source terms, if specified */
    if (Sminus) {
      c = cb;
      for (k = kb; k <= ks; k++) {
        B[k] += dz[k] * Sminus[c] / var[c];
        c = window->zp1[c];
      }
    }

    /*---------------------------------------------------------------*/
    /* Solve tridiagonal system */
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
    /* Update the variable */
    c = cs;
    var[cs] += (sol[ks] - var[cs]) * scale[c2];
    for (k = ks - 1; k >= kb; k--) {
      c = window->zm1[c];
      sol[k] -= ud[k + 1] * sol[k + 1];
      var[c] += (sol[k] - var[c]) * scale[c2];
    }

    /*---------------------------------------------------------------*/
    /* Set the concentration in thin layers */
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
void vert_diffuse_2d(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int nn;                       /* Tracer index */
  int cc;                       /* Counter */
  int c, c2, cb;                /* 3D and 2D sparse coordinate */
  double *topflux;              /* Surface flux */
  double top, bot;              /* Surface and bottom heights */
  double Cv = 4e3;              /* Specific heat at constant volume */

  /*-----------------------------------------------------------------*/
  /* Assignment of work arrays */
  /* d1 = flux out of top layer */
  memset(wincon->d1, 0, window->sgsizS * sizeof(double));

  for (nn = 0; nn < wincon->ntdif_v; nn++) {
    int n = wincon->tdif_v[nn];
    double *tr = windat->tr_wc[n];

    /* Set the surface boundary condition */
    if (n == windat->tno) {
      if (wincon->heatflux & (ADVANCED | INVERSE | NET_HEAT | COMP_HEAT | COMP_HEAT_MOM)) {
        /* Add the shortwave to the heatflux */
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
  /* Set the water column tracer values equal to the surface */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values */

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
void auxiliary_routines(geometry_t *window,  /* Processing window */
			window_t *windat,  /* Window data structure */
			win_priv_t *wincon /* Window geometry / constants */
  )
{
  int nn;
  int cc, c, c2, cb;
  double h1, h2, diff;

  /*-----------------------------------------------------------------*/
  /* Reset dz and the cells to process to correspond to the updated */
  /* elevation. */
  reset_dz(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get the mixed layer depth diagnostic if required */
  if (!(wincon->mixlayer & NONE))
    get_mixed_layer(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Flux diagnostics are calculated in advect_diffuse() using the   */
  /* actual fluxes used to update tracer concentration (i.e.         */
  /* advection scheme dependent).                                    */
  /*
  if (wincon->trflux >= 0)
    calc_flux_old(window, windat, wincon, dtu);
  */

  /*-----------------------------------------------------------------*/
  /* Calculate the flushing diagnostic if required */
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
  /* Calculate the steric height if required */
  if (wincon->lnm != 0.0)
    steric(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Calculate the vorticity if required */
  if (!(wincon->vorticity & NONE))
    vorticity(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Get diagnostic numbers if required                              */
  if (!(wincon->numbers & NONE))
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
  /* Evaluate the tracer decay */
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
  /* Do the ecological routines */
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

  /* Step 0 must always be the last call to tracerstats() */
#if defined(HAVE_TRACERSTATS_MODULE)
  /*tracerstats_prestep(window,0);*/
#endif

  if (wincon->waves & BOT_STR)
    memcpy(wincon->Cd, windat->wave_Cd, window->sgsizS * sizeof (double));

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
void wc_tracer_2d(geometry_t *window, /* Processing window */
                  window_t *windat, /* Window data structure */
                  win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int nn;                       /* Tracer index */
  int cc;                       /* Counter */
  int c, c2, cb;                /* 3D and 2D sparse coordinate */

  /* Density */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = c2 = window->w2_t[cc];
    cb = window->bot_t[cc];
    while (c != cb) {
      c = window->zm1[c];
      windat->dens[c] = windat->dens[c2];
    }
  }

  /* Tracers */
  for (nn = 0; nn < wincon->ntbdy; nn++) {
    int n = wincon->tbdy[nn];
    double *tr = windat->tr_wc[n];  /* Tracer values */

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
void bdry_tracer(geometry_t *window,  /* Processing window           */
                 window_t *windat,    /* Window data structure       */
                 win_priv_t *wincon   /* Window geometry / constants */
  )
{
  int n, nn;                    /* Counters                          */
  int tn;                       /* Tracer number                     */
  open_bdrys_t **open = window->open;

  if (wincon->trasc == FFSL)
    prep_semi_lagrange_atc(window);

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
                geometry_t *window, /* Processing window             */
                window_t *windat,   /* Window data structure         */
                win_priv_t *wincon, /* Window geometry / constants   */
                open_bdrys_t *open  /* OBC structure                 */
  )
{
  int c, cc;                    /* Sparse coorindate / counter       */
  int c1, c2;                   /* Sparse coordinates                */
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
    int ii, si, osi, ci, nn;
    double *vals;

    if (open->ocodex & (L_EDGE | R_EDGE)) {
      imap = open->nmap;
      vel = windat->u1;
    } else if (open->ocodey & (B_EDGE | F_EDGE)) {
      imap = open->nmap;
      vel = windat->u2;
    }

    for (cc = 1; cc <= open->no3_t; cc++) {
      c = osi = open->obc_t[cc];

      /* Get the range over which to compute the statistics. This is */
      /* taken as the cell where velocity changes by a certain       */
      /* fraction of the boundary velocity. Note : the boundary      */
      /* value is not included in the subset.                        */
      si = imap[c];
      bvel = vel[c];
      fvel = fabs(vel_frac * bvel);
      ci = 0;
      while (fabs(vel[si] - bvel) < fvel && si != imap[si]) {
        osi = si;
        si = imap[si];
        ci += 1;
      }
      si = osi;

      /* Store the tracer values over the range in the array vals[]  */
      if (ci > 0) {
        vals = d_alloc_1d(ci);
        nn = 0;
        ii = c;
        while (ii != si) {
          ii = imap[ii];
          vals[nn] = tr[ii];
          nn++;
        }
        /* Get the statistical value                                 */
        newval[c] = percs(vals, ci, perct);
        d_free_1d(vals);
      }
      /* Set a no gradient condition                                 */
      else {
        c1 = open->oi1_t[cc];
        newval[c] = tr[c1];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Least squares linear interpolation                              */
  else if (bcond & LINEXT) {
    imap = open->nmap;
    for (cc = 1; cc <= open->no3_t; cc++) {
      int ei = open->obc_t[cc];
      int ei0 = imap[ei];
      int ei1 = imap[ei0];
      int ei2 = imap[ei1];
      int ei3 = imap[ei2];
      newval[ei] = bc_leastsq(tr[ei3], tr[ei2], tr[ei1], tr[ei0]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Polynomial extrapolation (2nd order)                            */
  else if (bcond & POLEXT) {
    imap = open->nmap;
    for (cc = 1; cc <= open->no3_t; cc++) {
      int ei = open->obc_t[cc];
      int ei0 = imap[ei];
      int ei1 = imap[ei0];
      int ei2 = imap[ei1];
      int ei3 = imap[ei2];
      newval[ei] = bc_polint(tr[ei3], tr[ei2], tr[ei1]);
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
    int cs;               /* Surface sparse coordinate               */
    int cb;               /* Bottom sparse coordinate                */
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
    int cs;               /* Surface sparse coordinate               */
    int cb;               /* Bottom sparse coordinate                */
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

  /*-----------------------------------------------------------------*/
  /* Set an upstream advection condition on tracer boundaries        */
  if (bcond & UPSTRM) {
    double sc = 1.0;
    double v;
    int *cv = NULL, *civ = NULL;
    if (open->ocodex & (L_EDGE | R_EDGE)) {
      vel = windat->u1;
      hat = window->h1au1;
      cv = open->obc_e1;
      civ = open->oi1_e1;
    } else if (open->ocodey & (B_EDGE | F_EDGE)) {
      vel = windat->u2;
      hat = window->h2au2;
      cv = open->obc_e2;
      civ = open->oi1_e2;
    }
    if (open->ocodex & L_EDGE || open->ocodey & B_EDGE)
      sc = -1.0;
    if (wincon->trasc == LAGRANGE) {
      /* Characteristic method for use with LAGRANGE                 */
      double f;
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c2 = window->m2d[c];
	v = 0.5 * (vel[cv[cc]] + vel[civ[cc]]);
	f = min(fabs(v) * windat->dt / hat[c2], 1.0);
	f = 1.0;
	newval[c] = - sc * v > 0.0 ?
	  tr[c] - f * (tr[c] - newval[c]) : tr[c];
      }
    } else if (wincon->trasc == FFSL) {
      int *cl = wincon->s2;
      /* Update boundary cells with newval (e.g. FILEIN, CUSTOM)     */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c1 = cv[cc];
	tr[c] = tr[c1] = newval[c];
      }
      /* Save with the upstream value to newval                      */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	/* Get the streamline origin                                 */
	if (!cl[c])
	  streamline_atc(window, windat, wincon, c);
	/* Get the tracer value                                      */
	semi_lagrange_atc(window, windat, wincon, tr, newval, c);
      }
    } else {
      
      /* Velocities used are the mean of boundary & interior         */
      if (open->upmeth & CENTER) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  c1 = open->oi1_t[cc];
	  v = 0.5 * (vel[cv[cc]] + vel[civ[cc]]);
	  c2 = window->m2d[c];
	  newval[c] = tr[c] - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	       + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
	}
      }
      /* Velocities used are the boundary face velocities              */
      else if (open->upmeth & FACE) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  c1 = open->oi1_t[cc];
	  v = vel[cv[cc]];
	  c2 = window->m2d[c];
	  newval[c] = tr[c] - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	       + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
	}
      }
      /* Velocities used are depend on the centered velocity direction */
      else if (open->upmeth & ADAPTIVE) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  c1 = open->oi1_t[cc];
	  v = - sc * 0.5 * (vel[cv[cc]] + vel[civ[cc]]) > 0.0 ?
	    vel[cv[cc]] : vel[civ[cc]];
	  c2 = window->m2d[c];
	  newval[c] = tr[c] - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	       + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
	}
      }
      /* Velocities used are located one cell into the interior        */
      else {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->obc_t[cc];
	  c1 = open->oi1_t[cc];
	  v = vel[c1];
	  c2 = window->m2d[c];
	  newval[c] = tr[c] - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	       + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
	}
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Relaxation to external solutions                                */
  if (rlxn) {
    double alpha;
    int bn;
    imap = open->nmap;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->oi1_t[cc];
      bn = 2;
      while (c != imap[c] && bn <= rlxn) {
        alpha = 1 - tanh(0.5 * (double)(bn - 1));
        tr[c] = alpha * newval[c] + (1.0 - alpha) * tr[c];
        c = imap[c];
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

  /*-----------------------------------------------------------------*/
  /* Set the values in OBC ghost cells for TRCONC                    */
  /*
  if (open->bgz && !(bcond & TRCONC)) {
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->ogc_t[cc];
      c2 = open->obc_t[cc];
      for (c1 = 0; c1 < open->bgz; c1++) {
	tr[c] = newval[c2];
	c = open->omap[c];
      }
    }
  }
  */
  
}

/* END set_OBC_tracer()                                              */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Routine to calculate upstream advection tracer value              */
/*-------------------------------------------------------------------*/
void upstrm(geometry_t *window,  /* Processing window                */
	    window_t *windat,    /* Window data                      */
	    win_priv_t *wincon,
	    open_bdrys_t *open,
	    double *newval,
	    double *tr,
	    int mode)
{
  int c, cc, c1, c2;
  double sc = 1.0;
  double v;
  double *vel = NULL, *hat = NULL;
  int *cv = NULL, *civ = NULL;
  if (open->ocodex & (L_EDGE | R_EDGE)) {
    vel = windat->u1;
    hat = window->h1au1;
    cv = open->obc_e1;
    civ = open->oi1_e1;
  } else if (open->ocodey & (B_EDGE | F_EDGE)) {
    vel = windat->u2;
    hat = window->h2au2;
    cv = open->obc_e2;
    civ = open->oi1_e2;
  }
  if (open->ocodex & L_EDGE || open->ocodey & B_EDGE)
    sc = -1.0;
  if (wincon->trasc == LAGRANGE) {
    /* Characteristic method for use with LAGRANGE                   */
    double f;
    for (cc = 1; cc <= open->no3_t; cc++) {
      c = open->obc_t[cc];
      c2 = window->m2d[c];
      v = 0.5 * (vel[cv[cc]] + vel[civ[cc]]);
      f = min(fabs(v) * windat->dt / hat[c2], 1.0);
      f = 1.0;
      newval[c] = - sc * v > 0.0 ?
	tr[c] - f * (tr[c] - newval[c]) : tr[c];
    }
  }
  else {    
    /* Velocities used are the mean of boundary & interior           */
    if (mode & CENTER) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c1 = open->oi1_t[cc];
	v = 0.5 * (vel[cv[cc]] + vel[civ[cc]]);
	c2 = window->m2d[c];
	newval[c] = tr[c] - windat->dt / hat[c2]
	  * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
      }
    }
    /* Velocities used are the boundary face velocities              */
    else if (mode & FACE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c1 = open->oi1_t[cc];
	v = vel[cv[cc]];
	c2 = window->m2d[c];
	newval[c] = tr[c] - windat->dt / hat[c2]
	  * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
      }
    }
    /* Velocities used are depend on the centered velocity direction */
    else if (mode & ADAPTIVE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c1 = open->oi1_t[cc];
	v = - sc * 0.5 * (vel[cv[cc]] + vel[civ[cc]]) > 0.0 ?
	  vel[cv[cc]] : vel[civ[cc]];
	c2 = window->m2d[c];
	newval[c] = tr[c] - windat->dt / hat[c2]
	  * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
      }
    }
    /* Tracer in boundary ghost cell is updated using face velocity  */
    else if (mode & GHOST) {
      double d1;
      if (open->rlen) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->ogc_t[cc];
	  c1 = open->obc_t[cc];
	  v = vel[cv[cc]];
	  c2 = window->m2d[c1];
	  newval[c] +=  - windat->dt *
	    (0.5 * (sc * v + fabs(v)) * (newval[c] - tr[open->nmap[c]]) / hat[c2]
	       + (sc * v - fabs(v)) * (tr[open->omap[c]] - newval[c]) / open->rlen);
	}
      } else {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = open->ogc_t[cc];
	  c1 = open->obc_t[cc];
	  v = vel[cv[cc]];
	  c2 = window->m2d[c1];
	  newval[c] +=  - windat->dt / hat[c2]
	    * (0.5 * (sc * v + fabs(v)) * (newval[c] - tr[open->nmap[c]])
	       + 0.5 * (sc * v - fabs(v)) * (tr[open->omap[c]] - newval[c]));
	}
      }
    } else if (mode & VANLEER) {
      int *p1 = (open->type & U1BDRY) ? window->xp1 : window->yp1;
      int *m1 = (open->type & U1BDRY) ? window->xm1 : window->ym1;
      set_map_t(window);
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	wincon->w6[c] = vel[c] * windat->dt / hat[window->m2d[c]];
      }
      van_leer_do(newval, tr, vel, wincon->w6, 1, open->no3_t,
		  open->obc_t, p1, m1);
      reset_map_t_all(window);
    } else if (mode & ORDER2_UW) {
      double alp = -1.0;
      int *p1 = (open->type & U1BDRY) ? window->xp1 : window->yp1;
      int *m1 = (open->type & U1BDRY) ? window->xm1 : window->ym1;
      int p,m,pp,mm;
      set_map_t(window);
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c = m1[c];
	p = p1[c];
	pp = p1[p];
	m = m1[c];
	mm = m1[m];
	if (vel[c] < 0.0)
	  newval[c] = 0.5 * (tr[p] + tr[c]) - alp * (tr[p] - 2.0 * tr[c] + tr[m]);
	else 
	  newval[c] = 0.5 * (tr[c] + tr[m]) - alp * (tr[c] - 2.0 * tr[m] + tr[mm]);
	wincon->w6[c] = vel[c] * windat->dt / hat[window->m2d[c]];
      }
      /*
      order2_upwind_do(newval, tr, vel, wincon->w6, 1, open->no3_t,
		    open->obc_t, p1, m1);
      */
      reset_map_t_all(window);
    } else {
      /* Velocities used are located one cell into the interior        */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	c1 = open->oi1_t[cc];
	v = vel[c1];
	c2 = window->m2d[c];
	newval[c] = tr[c] - windat->dt / hat[c2]
	  * (0.5 * (sc * v + fabs(v)) * (tr[c] - tr[c1])
	     + 0.5 * (sc * v - fabs(v)) * (newval[c] - tr[c]));
      }
    }
  }
}

/* END upstrm()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to save the boundary tracer values to buffer and reset    */
/*-------------------------------------------------------------------*/
void save_OBC_tr(geometry_t *window,  /* Processing window           */
		 window_t *windat,    /* Window data                 */
		 win_priv_t *wincon , /* Window constants            */
		 double *Fx,
		 double *Fy,
		 double *tr,          /* Tracer                      */
		 int tn,              /* Tracer number               */
		 int mode
		 )
{
  int n, c, cc;

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
	  if (open->type & U1BDRY) {
	    for (cc = 1; cc <= open->no3_e1; cc++) {
	      c = open->obc_e1[cc];
	      open->dumn[tn][cc] = Fx[c];
	    }
	  } else {
	    for (cc = 1; cc <= open->no3_e2; cc++) {
	      c = open->obc_e2[cc];
	      open->dumn[tn][cc] = Fy[c];
	    }
	  }
	}
	if (mode == 4) {
	  if (open->type & U1BDRY) {
	    for (cc = 1; cc <= open->no3_e1; cc++) {
	      c = open->obc_e1[cc];
	      Fx[c] = open->dumn[tn][cc];
	    }
	  } else {
	    for (cc = 1; cc <= open->no3_e2; cc++) {
	      c = open->obc_e2[cc];
	      Fy[c] = open->dumn[tn][cc];
	    }
	  }
	}
      }
    }

    /*
    if (open->options & OP_RLEN && open->rlen && open->bcond_tra[tn] & (TRFLUX|TRCONC|TRCONF)) {
      double *length, *width;
      length = (open->type & U1BDRY) ? window->h1acell : window->h2acell;
      width = (open->type & U1BDRY) ? window->h2acell : window->h1acell;
      for (cc = 1; cc <= open->no2_t; cc++) {
	c = open->obc_t[cc];
	if (mode == 1) {
	  open->dum[c] = length[c];
	  length[c] = open->rlen;
	} else {
	  length[c] = open->dum[c];
	}
	window->cellarea[c] = length[c] * width[c];
	wincon->d3[c] = window->h1acell[window->xm1[c]] / (window->h1acell[window->xm1[c]] +
							   window->h1acell[c]);
	wincon->d3[window->xp1[c]] = window->h1acell[c] / (window->h1acell[c] +
							   window->h1acell[window->xp1[c]]);
      }
    }
    */
  }
}

/* END save_OBC_tr()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the tracer flux divergence at open boundary      */
/* cells. Note: for TRFLUX positive fluxes imply import of tracer,   */
/* hence the negative values for R_EDGE and F_EDGE, where a positive */
/* value on these edges usually implies tracer export.               */
/* For TRCONF the sign of the 3D flux dictates the direction, and    */
/* the negative sign is not required.                                */
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
	       double *Fy,          /* Horizontal e2 fluxes          */
	       double dt,           /* Timestep                      */
	       int tn,              /* Tracer number                 */
	       int mode     
	       )
{
  int c, cc, c2, n, m, nn, zm1, nc;
  double sgn = 1.0;

  /* Return if no tracers in this window have OBCs TRFLUX            */
  if (!(wincon->obctr[tn] & (TRFLUX|TRCONF))) return;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    double *flux;
    double tf = 0.0, f = 0.0;
    int nvec, *vec;
    if (!(open->bcond_tra[tn] & (TRFLUX|TRCONF))) continue;
    /* Get the index in the tracer dummy array corresponding to tn   */
    for (nn = 0; nn < open->ntflx; nn++)
      if (open->tflx[nn] == tn) break;

    /* Save the contribution to total flux through this OBC          */
    if (wincon->obctr[tn] & TRFLUX) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	double val = open->t_transfer[open->trm[tn]][cc];
	c = open->obc_t[cc];
	open->dumtr[nn][cc] = val / open->ncells;
      }
    } else if (wincon->obctr[tn] & TRCONF) {
      sgn = (open->ocodex & L_EDGE || open->ocodey & B_EDGE) ? 1.0 : -1.0;
      /* open->dumtr is set in bdry_transfer_tr(), after which a     */
      /* no-gradient is set so that higher order advection schemes   */
      /* at the first interior face do not use modified boundary     */
      /* ghost values.                                               */
      /* Reset the boundary ghost values to the original prescribed  */
      /* values here.                                                */
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->ogc_t[cc];
	windat->tr_wc[tn][c] = open->dumtr[nn][cc];
	/* Save the concentration so that open->dumtr[nn][cc] can    */
	/* be reset if sub-stepping is performed.                    */
	open->dumn[tn][cc] = open->dumtr[nn][cc];
      }
      /* Use the advection scheme to modify open->dumtr              */
      /* Not yet implemented */

      /* Multipy by the flux for TRCONF (assumes the input value is  */
      /* a concentration rather than a flux).                        */
      if (open->type & U1BDRY) {
	flux = windat->u1flux3d;
	vec = open->obc_e1;
	nvec = open->no3_e1;
      } else {
	flux = windat->u2flux3d;
	vec = open->obc_e2;
	nvec = open->no3_e2;
      }
      for (cc = 1; cc <= nvec; cc++) {
	c = vec[cc];
	tf += flux[c];
	if (sgn * flux[c] > 0.0) f += flux[c];
      }
      f = (f) ? fabs(tf / f) : 0.0;
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (open->ocodex & R_EDGE || open->ocodey & F_EDGE) c = open->omap[c];
	tf = (sgn * flux[c] > 0.0) ? f : 0.0;
	/* Get the mass flux scaled for inflow only                */
	open->dumtr[nn][cc] *= (tf * flux[c]);
      }
    }

    /* Reset the fluxes. Note: positive dumtr fluxes imply import of */
    /* tracer into the boundary cell.                                */
    if (open->ocodex & L_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (mode)
	  Fx[c] = open->dumtr[nn][cc];
	else
	  dtracer[c] += (Fx[c] - open->dumtr[nn][cc]) * dt;
      }
    }
    if (open->ocodex & R_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (mode)
	  Fx[window->xp1[c]] = -sgn * open->dumtr[nn][cc];
	else
	  dtracer[c] += (-Fx[window->xp1[c]] - sgn * open->dumtr[nn][cc]) * dt;
      }
    }
    if (open->ocodey & B_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (mode)
	  Fy[c] = open->dumtr[nn][cc];
	else
	  dtracer[c] += (Fy[c] - open->dumtr[nn][cc]) * dt;
      }
    }
    if (open->ocodey & F_EDGE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	if (mode)
	  Fy[window->yp1[c]] = -sgn * open->dumtr[nn][cc];
	else
	  dtracer[c] += (-Fy[window->yp1[c]] - sgn * open->dumtr[nn][cc]) * dt;
      }
    }
    /* Reset the flux for TRCONF if sub-stepping                     */
    if (open->bcond_tra[tn] & TRCONF) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	open->dumtr[nn][cc] = open->dumn[tn][cc];
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
void set_lateral_BC_tr(double **tr, /* Tracer array */
                       int ntr, /* Number of tracers */
                       int sgbpt, /* Number of boundary cells */
                       int *bpt,  /* Boundary cell vector */
                       int *bin /* Interior cells to bpt */
  )
{
  int c1, c2, cc, n, nn;            /* Counters */

  /* Set the boundary conditions (no flux) */
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
/* Set a no-gradient over R_EDGE and F_EDGE OBCs                     */
/*-------------------------------------------------------------------*/
void set_lateral_bdry_tr(geometry_t *window, window_t *windat, win_priv_t *wincon)
{
  int c, cc, n, nn, tn, co;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if ((open->ocodex & R_EDGE || open->ocodey & F_EDGE)
	&& open->stagger & OUTFACE) {
      for (cc = 1; cc <= open->no3_t; cc++) {
	c = open->obc_t[cc];
	co = open->omap[c];
	for (nn = 0; nn < wincon->ntbdy; nn++) {
	  tn = wincon->tbdy[nn];
	  if (!(open->bcond_tra[tn] & (TRCONC|TRCONF)))
	    windat->tr_wc[tn][co] = windat->tr_wc[tn][c];
	}
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
void Set_lateral_BC_density_w(double *dens, /* Density array */
                              int sgbpt,  /* Number of boundary cells */
                              int *bpt, /* Boundary cell vector */
                              int *bin  /* Interior cells to bpt */
  )
{
  int c1, c2, cc;               /* Counters */

  /* Set the boundary conditions (no gradient) */
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
void set_lateral_BC_vel(double *vel,  /* Velocity array */
                        int sgbpt,  /* Number of boundary cells */
                        int *bpt, /* Boundary cell vector */
                        int *bin  /* Interior cells to bpt */
  )
{
  int c, cc;                    /* Counters */

  /* Set the boundary conditions (no slip) */
  for (c = 1; c <= sgbpt; c++) {
    cc = bpt[c];
    vel[cc] = 0.0;
  }
}

/* END set_lateral_BC_vel()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to set the cell thickness and the sparse coordinate of    */
/* the surface.                                                      */
/*-------------------------------------------------------------------*/
void set_dz(geometry_t *window, /* Processing window */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon  /* Window geometry / constants */
  )
{
  int cc;                       /* Counter */
  int c, c3, c2, cs;            /* Sparse coordinate */
  int zm1;                      /* Sparse cell below cell c */
  int cb, cbot;                 /* Sparse coordinate of the bottom */
  double top;                   /* Top face of a cell */
  double bot;                   /* Bottom face of a cell */
  int vc;                       /* Tracer cells to provess counter */
  int vcs;                      /* Surface tracer cells to process counter
                                 */
  /*-----------------------------------------------------------------*/
  /* Sigma coordinates */
  if (wincon->sigma) {
    int size = window->a2_t + 1;
    /* Get the cells to process vectors */
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
      /* Get the 2D and 3D sparse coordinates of the old surface */
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

      /* Set the cell thickness at the bottom */
      wincon->dz[cbot] = (top - window->botz[c2]);
    }
  }
  /*-----------------------------------------------------------------*/
  /* 'z' coordinates */
  else {
    /*---------------------------------------------------------------*/
    /* Get the vector containing the sparse locations of the new */
    /* surface (i.e. the surface after the 2D mode is complete).  */
    /* Also reset the old surface vector. This is performed over all */
    /* wet + auxiliary cells.  */
    for (cc = 1; cc <= window->a2_t; cc++) {

      /* Get the 2D and 3D sparse coordinates */
      c = cs = window->w2_t[cc];
      cb = window->bot_t[cc];

      /* Copy the new surface into the old surface vector */
      window->sur_t[cc] = window->nsur_t[cc];

      /* Get the new sparse location of the surface. For elevations */
      /* below the bottom this is set to the sediment coordinate.  */
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
    /* Get the wet cells to process vector for tracers and store in */
    /* the buffer wincon->s1. This is performed for wet cells and */
    /* open boundary cells only (open boundary cells are required */
    /* to calculate vertical velocity, used in momentum advection.  */
    /* However, tracers are also updated on the boundary which uses */
    /* T at t+1 on the boundary in the UPSTRM OBC).  */
    /* Cells are arranged thus : */
    /* 1 to wincon->vcs1 = surface layer wet cells */
    /* wincon->vcs1+1 to wincon->vcs = surface layer OBC cells */
    /* wincon->vcs+1 to wincon->vc1 = sub-surface wet cells */
    /* wincon->vc1+1 to wincon->vc = sub-surface OBC cells */
    /* First get the sparse locations of the surface : this is the */
    /* lower of the surface elevations before and after the 2D mode. */
    vc = 1;
    for (cc = 1; cc <= window->b2_t; cc++) {
      /* The lower of the layers corresponding to current and */
      /* previous surface elevations is required here. The sparse */
      /* coordinates in the window system (local coordinates) are */
      /* always smallest in the surface layer, and increase towards */
      /* the bottom. Hence the maximum of old and current layer */
      /* locations will provide the lower level of the two.  */
      /* Store the bottom coordinate for the cells to process (1 to */
      /* vcs) in the buffer i1. This is required in surf_conc().  */
      c = max(window->sur_t[cc], window->nsur_t[cc]);
      cs = window->m2d[c];
      top = wincon->oldeta[cs];
      top = min(windat->eta[cs], wincon->oldeta[cs]);
      bot = DRY_FRAC * wincon->hmin + window->botz[cs];
      if (top > bot && c != window->zm1[c]) {
        wincon->s1[vc] = c;
	wincon->i1[vc] = window->bot_t[cc];
        wincon->i2[vc] = window->nsur_t[cc];
        wincon->i3[vc] = window->sur_t[cc];
        if (cc == window->v2_t)
          wincon->vcs1 = vc;
        vc++;
      }
    }
    vcs = vc - 1;
    wincon->vcs = vcs;

    /* Loop from the layer below the surface to the bottom and get */
    /* the cells to process.  */
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
    /* the surface layer is the layer thickness at the start of the */
    /* time-step; i.e. eta-gridz if calculated before the call to */
    /* etastep() or oldeta-gridz if calculated after etastep().  */
    for (cc = 1; cc <= window->b2_t; cc++) {
      /* Get the 2D and 3D sparse coordinates of the old surface */
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
      /* Set the cell thickness at the bottom */
      wincon->dz[cbot] = top - window->botz[c2];

      /* Set all cell thickness above the old surface equal to the */
      /* surface thickness (used in higher order vertical advection */
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
void set_thin_tr(geometry_t *window,  /* Processing window */
                 window_t *windat,  /* Window data structure */
                 win_priv_t *wincon /* Window geometry / constants */
  )
{
  int cc;                       /* Counter */
  int c;                        /* Sparse coordinate */
  int cs;                       /* Sparse coordinate of the surface */

  /*-----------------------------------------------------------------*/
  /* Set the thin layer vector and reset the surface layer.  */
  /* Note : if the loop is done to a2_e1 then dzu1=0 at the first */
  /* timestep for multiple windows in auxiliary cells since it has */
  /* not yet been copied from the master.  */
  wincon->nkth_e1 = 1;
  for (cc = 1; cc <= window->a2_t; cc++) {
    c = window->nsur_t[cc];
    cs = window->m2d[c];
    if (c != window->bot_t[cc] &&
        windat->eta[cs] - window->gridz[c] < wincon->hmin) {
      window->nsur_t[cc] = window->zm1[c];
      /* Save the sparse counter of the thin layer */
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
void reset_dz(geometry_t *window,  /* Processing window */
	      window_t *windat,  /* Window data structure */
	      win_priv_t *wincon /* Window geometry / constants */
	      )
{
  int cc;                       /* Counter */
  int c, cs, cn;                /* Sparse coordinate */
  int zp1, zm1;                 /* Sparse cell above / below cell c */
  double top;                   /* Top face of a cell */
  double bot;                   /* Bottom face of a cell */
  int vc;                       /* Tracer cells to provess counter */
  int vcs;                      /* Surface tracer cells to process counter */

  /*-----------------------------------------------------------------*/
  /* Set dz at the surface to account for the updated elevation. */
  for(cc=1; cc <= window->b2_t; cc++) {
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
    /* Set the cells to process vector for the updated elevation. */
    /* Cells are arranged thus : */
    /* 1 to wincon->vca2 = surface layer wet cells */
    /* wincon->vca2+1 to wincon->vcs2 = surface layer OBC cells */
    /* wincon->vcs2+1 to wincon->vci2 = sub-surface wet cells */
    /* wincon->vci2+1 to wincon->vc2 = sub-surface OBC cells */
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

    /* Loop from the layer below the surface to the bottom and get */
    /* the cells to process.  */
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
  memset(wincon->c1, 0, window->sgsiz * sizeof(char));
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
void set_map_t(geometry_t *window /* Window data structure */
  )
{
  int n;
  int cc, cs, c;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->stagger & OUTFACE) {
      if (open[n]->bgz) {
	for (cc = 1; cc <= open[n]->no3_t; cc++) {
	  c = open[n]->obc_t[cc];
	  open[n]->omap[c] = open[n]->ogc_t[cc];
	}
      } else {
	for (cc = 1; cc <= open[n]->no3_t; cc++) {
	  c = open[n]->obc_t[cc];
	  open[n]->omap[c] = c;
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
void reset_map_t(geometry_t *window /* Window data structure */
  )
{
  int n;
  int cc, cs, c;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodex & R_EDGE && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e1; cc++) {
	c = open[n]->obc_e1[cc];
	cs = window->xm1[c];
	window->xp1[cs] = c;
      }
    }
    if (open[n]->ocodey & F_EDGE  && open[n]->stagger & OUTFACE) {
      for (cc = 1; cc <= open[n]->no3_e2; cc++) {
	c = open[n]->obc_e2[cc];
	cs = window->ym1[c];
	window->yp1[cs] = c;
      }
    }
  }
}

void reset_map_t_all(geometry_t *window) {
  int n;
  int cc, c, *obce;
  open_bdrys_t **open = window->open;

  for (n = 0; n < window->nobc; n++) {
    if (open[n]->type & U1BDRY)
      obce = open[n]->obc_e1;
    if (open[n]->type & U2BDRY)
      obce = open[n]->obc_e2;
    for (cc = 1; cc <= open[n]->no3_t; cc++) {
      c = open[n]->obc_t[cc];
      open[n]->omap[c] = obce[cc];
    }
  }
}

/* END reset_map_t()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to reset the map across multiple ghost cell tracer open   */
/* boundaries for LAGRANGE.                                          */
/*-------------------------------------------------------------------*/
void set_map_l(geometry_t *window /* Window data structure */
  )
{
  int n;
  int cc, ci, c;
  open_bdrys_t **open = window->open;

  if (!(window->wincon->trasc & LAGRANGE)) return;

  /*-----------------------------------------------------------------*/
  /* Set the sparse mappings in open boundary ghost cells            */
  for (n = 0; n < window->nobc; n++) {
    if (open[n]->ocodex & R_EDGE) {
      for (cc = 1; cc <= open[n]->no3_t; cc++) {
	c = open[n]->obc_t[cc];
	ci = window->yp1[c];
	window->xp1[ci] = ci;
	ci = window->ym1[c];
	window->xp1[ci] = ci;
      }
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
	c = window->zm1[open[n]->bot_t[cc]];
	window->xp1[c] = c;
      }
    }
    if (open[n]->ocodey & F_EDGE) {
      for (cc = 1; cc <= open[n]->no3_t; cc++) {
	c = open[n]->obc_t[cc];
	ci = window->xp1[c];
	window->yp1[ci] = ci;
	ci = window->xm1[c];
	window->yp1[ci] = ci;
      }
      for (cc = 1; cc <= open[n]->no2_t; cc++) {
	c = window->zm1[open[n]->bot_t[cc]];
	window->yp1[c] = c;
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
void do_ts_relax(geometry_t *window,  /* Processing window */
		 window_t *windat,  /* Window data structure */
		 win_priv_t *wincon /* Window geometry / constants */
		 )
{
  int cc, c, c3;                       /* Counter */
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

  /* Calculate the mean value.  */
  if (pc == -1.0) {
    for (m = 0; m < n; m++)
      dm += aa[m];
    return (dm / (double)n);
  }
  /* use a cyclic value over the subset */
  else if (pc == -2) {
    return (aa[n - 1]);
  }
  /* Calculate the percentime value.  */
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
void get_mixed_layer(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int cc;                       /* Counter */
  int c, zm1;                   /* Sparse coordinate */
  int cs;                       /* Sparse coordinate of the surface */
  int cb;                       /* Sparse coordinate of the bottom */
  int kts;                      /* Index of the top of the pycnocline */
  int ktb;                      /* Index of the bottom of the pycnocline */
  double top_m;                 /* Depth of the top of the pycnocline */
  double bot_m;                 /* Depth of the bottom of the pycnocline */
  double thr = -0.01;           /* Density gradient threshold */

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
    memset(windat->mixl, 0, window->sgsizS * sizeof(double));
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
void calc_flux(geometry_t *window,  /* Processing window */
               window_t *windat,    /* Window data structure */
               win_priv_t *wincon,  /* Window geometry / constants */
	       double *fluxe1,
	       double *fluxe2,
	       double *fluxw,
	       double dt            /* Time step */
  )
{
  int c, cs, zm1, cc;           /* Sparse coordinate / counter */
  int tn = wincon->trflux;      /* Tracer to calculate fluxes for */
  double *tr;                   /* Tracer pointer */
  double flux;                  /* Flux value */

  /* double *fluxkz = wincon->w10; */

  /*-----------------------------------------------------------------*/
  /* Horizontal advective flux through e1 face                       */
  tr = windat->tr_wc[tn];

  if (windat->fluxe1) {
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_e1; cc++) {
	c = window->w3_e1[cc];
	cs = window->m2d[c];
	windat->fluxe1[c] = (windat->fluxe1[c] * windat->meanc[cs] + 
			     fluxe1[c] * dt) / (windat->meanc[cs] + dt);
      }
    } else
      memcpy(windat->fluxe1, fluxe1, window->sgsiz * sizeof(double));
  }

  /* Horizontal advective flux through e2 face                       */
  if (windat->fluxe2) {
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_e2; cc++) {
        c = window->w3_e2[cc];
	cs = window->m2d[c];
	windat->fluxe2[c] = (windat->fluxe2[c] * windat->meanc[cs] + 
			     fluxe2[c] * dt) / (windat->meanc[cs] + dt);
      }
    } else
      memcpy(windat->fluxe2, fluxe2, window->sgsiz * sizeof(double));
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
      windat->fluxw[c] = flux;
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
	windat->fluxkz[c] = flux;
      }
    }
  }
}


void calc_flux_old(geometry_t *window,  /* Processing window */
               window_t *windat,  /* Window data structure */
               win_priv_t *wincon, /* Window geometry / constants */
               double dtu
  )
{
  int c, cs, zm1, cc;           /* Sparse coordinate / counter */
  int tn = wincon->trflux;      /* Tracer to calculate fluxes for */
  double *tr;                   /* Tracer pointer */
  double *fluxe1 = wincon->w7;
  double *fluxe2 = wincon->w8;
  double *fluxw = wincon->w9;
  double *fluxkz = wincon->w10;

  /*-----------------------------------------------------------------*/
  /* Calculate the fluxes at each cell */
  tr = windat->tr_wc[tn];
  if (windat->fluxe1) {
    memset(fluxe1, 0, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= window->b3_e1; cc++) {
      c = window->w3_e1[cc];
      fluxe1[c] = windat->u1flux3d[c] * tr[c];
    }
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_e1; cc++) {
        c = window->w3_e1[cc];
        if (windat->meanc[window->m2d[c]] > 0)
        	windat->fluxe1[c] += (fluxe1[c] * dtu);
        /*UR if compared use sub-dt windat->fluxe1[c] += (fluxe1[c] * windat->dt); */
      }
    } else
      memcpy(windat->fluxe1, fluxe1, window->sgsiz * sizeof(double));
  }

  if (windat->fluxe2) {
    memset(fluxe2, 0, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= window->b3_e2; cc++) {
      c = window->w3_e2[cc];
      fluxe2[c] = windat->u2flux3d[c] * tr[c];
    }
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_e2; cc++) {
        c = window->w3_e2[cc];
        if (windat->meanc[window->m2d[c]] > 0)
          windat->fluxe2[c] += (fluxe2[c] * dtu);
          /*UR if compared use sub-dt windat->fluxe2[c] += (fluxe2[c] * windat->dt);*/
      }
    } else
      memcpy(windat->fluxe2, fluxe2, window->sgsiz * sizeof(double));
  }

  if (windat->fluxw) {
    memset(fluxw, 0, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      fluxw[c] = 0.5 * (windat->w[c] + windat->w[zm1]) *
        window->cellarea[cs] * tr[c];
    }
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
        if (windat->meanc[window->m2d[c]] > 0)
          windat->fluxw[c] += (fluxw[c] * dtu);
          /*UR if compared use sub-dt windat->fluxw[c] += (fluxw[c] * windat->dt);*/
      }
    } else
      memcpy(windat->fluxw, fluxw, window->sgsiz * sizeof(double));
  }
  if (windat->fluxkz) {
    memset(fluxkz, 0, window->sgsiz * sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      cs = window->m2d[c];
      zm1 = window->zm1[c];
      fluxkz[c] =
        -windat->Kz[c] * window->cellarea[cs] * (tr[c] -
                                                 tr[zm1]) / wincon->dz[c] *
        wincon->Ds[cs];
    }
    if (wincon->means & FLUX) {
      for (cc = 1; cc <= window->b3_t; cc++) {
        c = window->w3_t[cc];
        if (windat->meanc[window->m2d[c]] > 0)
          windat->fluxkz[c] += (fluxkz[c] * dtu);
          /*UR if compared use sub-dt windat->fluxkz[c] += (fluxkz[c] * windat->dt);*/
      }
    } else
      memcpy(windat->fluxkz, fluxkz, window->sgsiz * sizeof(double));
  }
}

/* END calc_flux()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the total mass in a masked region                      */
/*-------------------------------------------------------------------*/
void total_mass(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon  /* Window geometry / constants */
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
double check_mass(geometry_t *window,       /* Window geometry       */
		  window_t *windat,         /* Window data           */
		  win_priv_t *wincon,       /* Window constants      */
		  char *trname,             /* Tracer name           */
		  int cin                   /* Cell for balance      */
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
double check_tmass(geometry_t *window,       /* Window geometry       */
		  window_t *windat,         /* Window data           */
		  win_priv_t *wincon,       /* Window constants      */
		  char *trname              /* Tracer name           */
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
void init_flushing(master_t *master,  /* Master data structure */
                   geometry_t **window, /* Processing window */
                   win_priv_t **wincon  /* Window geometry / constants */
  )
{
  FILE *fp = master->prmfd;
  int c, cc, m, n, wn, lc, i, j;
  char keyword[MAXSTRLEN], buf[MAXSTRLEN];

  /*-----------------------------------------------------------------*/
  /* Get the tracer to flush. Note: the concentration of this */
  /* tracer will be initialized to 1 in the flushing domain and */
  /* zero elsewhere.  */
  master->flt = NaN;
  master->flts = master->tstart;
  if (master->trflsh) {
    /* Initialize the flushing tracer */
    memset(master->fltr, 0, geom->sgsiz * sizeof(double));
    master->mask = s_alloc_1d(geom->sgsizS);
    memset(master->mask, 0, geom->sgsizS * sizeof(short));
    sprintf(keyword, "FLUSHING_PTS");
    if (prm_read_int(fp, keyword, &n)) {
      prm_flush_line(fp);
      if (n > 0) {
	for (m = 0; m < n; ++m) {
	  /* Read i and j indices for this point */
	  if (fscanf(fp, "%d %d", &i, &j) != 2)
	    hd_quit("flushing: Can't read i j in points list.\n");
	  prm_flush_line(fp);
	  /* Set the mask and flushing tracer */
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
      read_blocks(fp, "FLUSHING_BLOCKS", &n, &iloc, &jloc);
      for (m = 1; m <= n; m++) {
        c = geom->map[geom->nz - 1][jloc[m]][iloc[m]];
        master->mask[c] = 1;
        while (c != geom->zm1[c]) {
          master->fltr[c] = 1.0;
          c = geom->zm1[c];
        }
      }
      i_free_1d(iloc);
      i_free_1d(jloc);
    } else if (prm_read_char(fp, "FLUSHING_REGION", buf)) {
      char *files[MAXSTRLEN * MAXNUMARGS];
      int nf, nr, rgn;
      double *regionid;
      nf = parseline(buf, files, MAXNUMARGS);
      regionid = d_alloc_1d(master->geom->sgsiz);
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

    /* Set the mask in the windows */
    if (master->nwindows > 1) {
      window_t *windat;
      for (wn = 1; wn <= master->nwindows; wn++) {
	windat = window[wn]->windat;
        wincon[wn]->mask = s_alloc_1d(window[wn]->sgsizS);
        memset(wincon[wn]->mask, 0, window[wn]->sgsizS * sizeof(short));
        memset(windat->fltr, 0, window[wn]->sgsiz * sizeof(double));
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
    /* Initialise the mass */
    master->imass = 0.0;
    for (wn = 1; wn <= master->nwindows; wn++) {
      total_mass(window[wn], window[wn]->windat, wincon[wn]);
      master->imass += wincon[wn]->imass;
    }

    /*---------------------------------------------------------------*/
    /* Initialize the output file */
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
  double e = 2.718281828;       /* Constant e */
  double t, mt = 0.0, fr = 0.0;
  int wn;

  /* Get the total mass in all windows */
  for (wn = 1; wn <= master->nwindows; wn++) {
    mt += wincon[wn]->imass;
  }
  /* If e-folding time is reached print the time */
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
  /* Reinitialise after 2 * flushing times */
  if (!isnan(master->flt) && master->t > master->flts + 1.5 * 86400.0 * master->flt) {
    int c, cc, gc;
    master->flt = NaN;
    master->flts = master->t;
    master->imass = 0.0;
    for (wn = 1; wn <= master->nwindows; wn++) {
      memset(windat[wn]->fltr, 0, window[wn]->sgsiz * sizeof(double));
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
void init_age(master_t *master,  /* Master data structure */
	      geometry_t **window, /* Processing window */
	      window_t **windat, /* Window data structure */
	      win_priv_t **wincon  /* Window geometry / constants */
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
    /* Initialize the flushing tracer */
    agemsk = s_alloc_1d(geom->sgsiz);
    memset(agemsk, 0, geom->sgsiz * sizeof(short));
    /* Get the range */
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
      read_blocks(fp, "AGE_TR", &n, &iloc, &jloc);
      for (m = 1; m <= n; m++) {
	c = geom->map[geom->nz - 1][jloc[m]][iloc[m]];
	while (c != geom->zm1[c] && geom->cellz[c] <= top && geom->cellz[c] >= bot) {
	  agemsk[c] = 1;
	  c = geom->zm1[c];
	}
      }
      n = 1;
      i_free_1d(iloc);
      i_free_1d(jloc);
    } else {
      int nr;
      double *regionid;
      regionid = d_alloc_1d(master->geom->sgsiz);
      nr = read_regioni(master, files[0], regionid);
      if (nr) {
	int rgn;
	for (m = 1; m < nf; m++) {
	  sscanf(files[m], "%d", &rgn);
	  for (cc = 1; cc <= geom->b2_t; cc++) {
	    c = geom->w2_t[cc];
	    if (rgn == (int)regionid[c]) {
	      while (c != geom->zm1[c] && geom->cellz[c] <= top && geom->cellz[c] >= bot) {
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
	wincon[n]->agemsk = s_alloc_1d(window[n]->sgsiz);
	memset(wincon[n]->agemsk, 0, window[n]->sgsiz);
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
void calc_age(geometry_t *window,           /* Window geometry       */
	      window_t *windat,             /* Window data           */
	      win_priv_t *wincon            /* Window constants      */
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
void init_trperc(master_t *master,    /* Master data structure */
		 geometry_t **window, /* Processing window */
		 window_t **windat,   /* Window data structure */
		 win_priv_t **wincon  /* Window geometry / constants */
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
    /* Initialize the flushing tracer */
    percmsk = s_alloc_1d(geom->sgsiz);
    memset(percmsk, 0, geom->sgsiz * sizeof(short));
    /* Get the range */
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
      read_blocks(fp, "PERC_REGION", &n, &iloc, &jloc);
      for (m = 1; m <= n; m++) {
	c = geom->map[geom->nz - 1][jloc[m]][iloc[m]];
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
      i_free_1d(jloc);
    } else {
      int nr;
      double *regionid;
      regionid = d_alloc_1d(master->geom->sgsiz);
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
		  if (geom->cellz[c] <= top && geom->cellz[c] >= bot) percmsk[c] = 1;
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
        wincon[n]->percmsk = s_alloc_1d(window[n]->sgsiz);
	memset(wincon[n]->percmsk, 0, window[n]->sgsiz * sizeof(short));
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
	wincon[n]->percmsk = s_alloc_1d(window[n]->sgsiz);
	memset(wincon[n]->percmsk, 0, window[n]->sgsiz * sizeof(short));
	for (cc = 1; cc <= window[n]->b3_t; cc++) {
	  c = window[n]->w3_t[cc];
	  wincon[n]->percmsk[c] = 1;
	}
	if (surff)wincon[n]->percmsk[0] = 1;
      }
    }
  }
}

/* END init_trperc()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void mode2d_tracer_init(geometry_t *window, /* Processing window */
                        window_t *windat, /* Window data structure */
                        win_priv_t *wincon  /* Window geometry / constants
                                             */
  )
{
  int c, cs, n, cc;             /* Sparse coordinate / counter */
  double *tr;                   /* Tracer pointer */
  double depth;                 /* Water depth */
  int size = window->a2_t + 1;

  /* Vertically average all tracers */
  for (n = 0; n < windat->ntr; n++) {
    memset(wincon->d1, 0, window->sgsizS * sizeof(double));
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

  /* Set the cells to process vector */
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

  /* Copy the 2d cells to process into 3d arrays */
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
  memcpy(window->w3_e2, window->w2_e2, (window->n2_e2 + 1) * sizeof(int));
  window->n3_e2 = window->n2_e2;
  window->a3_e2 = window->a2_e2;
  window->b3_e2 = window->b2_e2;
  window->v3_e2 = window->v2_e2;
}

/* END mode2d_tracer_init()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to initialise diagnostic tracers to zero after dumping to */
/* file.                                                             */
/*-------------------------------------------------------------------*/
void tr_diag_reset_w(geometry_t *window,  /* Processing window */
                     window_t *windat,  /* Window data structure */
                     win_priv_t *wincon /* Window geometry / constants */
  )
{
  int n, nn;

  if (windat->df_diagn_set) {
    /* Zero diagnostic tracers */
    for (nn = 0; nn < wincon->ndia; nn++) {
      n = wincon->diagn[nn];
      memset(windat->tr_wc[n], 0, window->sgsiz * sizeof(double));
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
void tr_diag_reset_m(master_t *master /* Master data structure */
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
/* No longer required : using sparse formatted files for reset does  */
/* not perform interpolation.                                        */
/*-------------------------------------------------------------------*/
void diff_mask(geometry_t *window,  /* Processing window             */
	       window_t *windat,    /* Window data structure         */
	       win_priv_t *wincon   /* Window geometry / constants   */
	       )
{
  double thr = 0.025;       /* Difference threshold                   */
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


/*-------------------------------------------------------------------*/
/* Sets the lateral boundary conditions for tangential radiation     */
/* stress. The stress adjacent to solid boundaries is placed in the  */
/* corresponding ghost cell, and replaced with the stress 2 cells    */
/* into the interior. Velocity at the boundary is computed using:    */
/* tendency(i) = dt*(stress(i+1)-stress(i-1)/2dx                     */
/*-------------------------------------------------------------------*/
void set_lateral_BC_rad(geometry_t *window, /* Processing window     */
		   window_t *windat,   /* Window data structure      */
		   win_priv_t *wincon  /* Window geometry / constants*/
	       )
{
  int c, lc, lc2, cc;
  double *Sxy = wincon->d2;
  double *Syx = wincon->d3;

  memcpy(Sxy, windat->wave_Sxy, window->sgsizS * sizeof(double));
  memcpy(Syx, windat->wave_Syx, window->sgsizS * sizeof(double));

  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    lc = window->bin[cc];
    lc2 = window->bin2[cc];

    windat->wave_Sxy[c] = 0.0;
    windat->wave_Syx[c] = 0.0;

    if (lc == window->ym1[c] && windat->wave_Sxy[c] == 0.)
      windat->wave_Sxy[c] = windat->wave_Sxy[lc2] + Sxy[lc];
    if (lc == window->yp1[c] && windat->wave_Sxy[c] == 0.)
      windat->wave_Sxy[c] = windat->wave_Sxy[lc2] - Sxy[lc];
    if (lc2) windat->wave_Sxy[lc] = windat->wave_Sxy[lc2];

    if (lc == window->xm1[c] && windat->wave_Syx[c] == 0.)
      windat->wave_Syx[c] = windat->wave_Syx[lc2] + Syx[lc];
    if (lc == window->xp1[c] && windat->wave_Syx[c] == 0.)
      windat->wave_Syx[c] = windat->wave_Syx[lc2] - Syx[lc];
    if (lc2) windat->wave_Syx[lc] = windat->wave_Syx[lc2];
  }
}

/* END set_lateral_bc_rad()                                          */


/*-------------------------------------------------------------------*/
/* Sets the lateral boundary conditions for wave induced force.      */
/* The stress adjacent to solid boundaries is placed in the          */
/* corresponding ghost cell, and replaced with the stress 2 cells    */
/* into the interior. Velocity at the boundary is computed using:    */
/* tendency(i) = dt*(stress(i+1)-stress(i-1)/2dx                     */
/*-------------------------------------------------------------------*/
void set_lateral_BC_waf(geometry_t *window, /* Processing window     */
		   window_t *windat,   /* Window data structure      */
		   win_priv_t *wincon  /* Window geometry / constants*/
	       )
{
  int c, lc, lc2, cc;
  double *Fx = wincon->d2;
  double *Fy = wincon->d3;

  memcpy(Fx, windat->wave_Fx, window->sgsizS * sizeof(double));
  memcpy(Fy, windat->wave_Fy, window->sgsizS * sizeof(double));

  for (cc = 1; cc <= window->nbptS; cc++) {
    c = window->bpt[cc];
    lc = window->bin[cc];
    lc2 = window->bin2[cc];

    windat->wave_Fx[c] = 0.0;
    windat->wave_Fy[c] = 0.0;

    if (lc == window->ym1[c] && windat->wave_Fx[c] == 0.)
      windat->wave_Fx[c] = windat->wave_Fx[lc2] + Fx[lc];
    if (lc == window->yp1[c] && windat->wave_Fx[c] == 0.)
      windat->wave_Fx[c] = windat->wave_Fx[lc2] - Fx[lc];
    if (lc2) windat->wave_Fx[lc] = windat->wave_Fx[lc2];

    if (lc == window->xm1[c] && windat->wave_Fy[c] == 0.)
      windat->wave_Fy[c] = windat->wave_Fy[lc2] + Fy[lc];
    if (lc == window->xp1[c] && windat->wave_Fy[c] == 0.)
      windat->wave_Fy[c] = windat->wave_Fy[lc2] - Fy[lc];
    if (lc2) windat->wave_Fy[lc] = windat->wave_Fy[lc2];
  }
}

/* END set_lateral_bc_waf()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to add the tracer increment to the relevant state         */
/* variable.                                                         */
/*-------------------------------------------------------------------*/
void do_ts_increment(geometry_t *window, /* Processing window     */
		     window_t *windat,   /* Window data structure      */
		     win_priv_t *wincon)
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
  int c, k;                     /* Sparse/layer coordinate           */
  int cs, ks;                   /* Surface coordinate                */
  int cb, kb;                   /* Bottom coordinate                 */
  int zm1;                      /* Sparse cell below c               */
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

  cs = c = ctp[cc];      /* Set cs to the surface sparse coordinate  */
  cb = cbt[cc];          /* Bottom sparse coordinate                 */
  c2 = window->m2d[c];   /* 2D sparse location corresponding to c    */
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

  /* Map the sparse arrays to local coordinates.                     */
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
  c = cb;                   /* Set c to the bottom sparse coordinate */
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
/* Computs the number of wet cells for each boundary                 */
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
  int cc, c, nn;
  double *vel, d1;
  int *vec;

  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn]; 
    double *vel;
    int *vec;
    if(open->type & U1BDRY) {
      vel = windat->u1;
      vec = open->obc_e1;
    } else {
      vel = windat->u2;
      vec = open->obc_e2;
    }
  }


  for (nn = 0; nn < window->nobc; nn++) {
    open_bdrys_t *open = window->open[nn]; 
    if (open->bcond_nor & CUSTOM && open->options & OP_PARFLOW && 
	open->options & (OP_GEOSTR|OP_YANKOVSKY)) {
      if(open->type & U1BDRY) {
	vel = windat->u1;
	vec = open->obc_e1;
      } else {
	vel = windat->u2;
	vec = open->obc_e2;
      }
      if (mode == 0) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = vec[cc];
	  d1 = vel[c];
	  vel[c] = open->flow[cc];
	  open->flow[cc] = d1;
	}
      }
      if (mode == 1) {
	for (cc = 1; cc <= open->no3_t; cc++) {
	  c = vec[cc];
	  vel[c] = open->flow[cc];
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
    depth = atof(fields[1]);
    }
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
	  master->tran_mean[lc] = master->swr_tran[lc] = 0.0;
	} else {
	  windat->tran_mean[c] = -windat->swr_tran[c];
	  master->tran_mean[lc] = -master->swr_tran[lc];
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
  if (wincon->attn_tr >= 0 && wincon->swr_type & SWR_ATTN) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      if (windat->attn_mean[c] < 0.0)
	windat->attn_mean[c] = -windat->tr_wcS[wincon->attn_tr][c];
    }
  }
  if (wincon->tran_tr >= 0 && wincon->swr_type & SWR_TRAN) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      c = window->w2_t[cc];
      if (windat->tran_mean[c] < 0.0) 
	windat->tran_mean[c] = -windat->tr_wcS[wincon->tran_tr][c];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Do the vertical diffusion over the ensemble                     */
  i1s = i2s = 0;
  i1e = i2e = 10;
  i1i = i2i = 1;
  /* Default range spans attn 0.02-0.07 with tran around the BGC     */
  /* estimate of 0.6.                                                */

  /* Spans the attenuation of Jerlov classes (0.03-0.127) and full   */
  /* transmission range.                                             */
  attns = trans = 1.0;
  attn0 = 0.037;
  attni = 0.01;
  tran0 = 0.0;
  trani = 0.1;
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
	/*tran = (double)i1 / 10.0;*/
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
