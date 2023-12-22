/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/slaves/transfers.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: transfers.c 7452 2023-12-13 03:45:26Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void  tr_warn(char *text, double w_m, double w_n, int c,
	      int i, int j, int k, int cg, int n, char *type, char tp);
void s2m_zoom(master_t *master, geometry_t *window, window_t *windat, 
	      int mode);
void copy_dens_profile(master_t *master, geometry_t *window, window_t *windat);

/*-------------------------------------------------------------------*/
/* win_data_fill_3d()  : Fills local 3D arrays from the master       */
/* win_data_fill2d()   : Fills local 2D arrays from the master       */
/* win_data_refill3d() : Re-fills updated 3D window data with master */
/* win_data_refill2d() : Re-fills updated 2D window data with master */
/* win_data_empty_3d()  : Empties local 3D arrays from the master    */
/* win_data_empty_2d()  : Empties local 2D arrays from the master    */
/* build_transfer_maps() : Sets up the transfer vectors              */
/* master_fill() : Fills the master with local data                  */
/* s2m_3d()      : Transfers a 3D local array to the master          */
/* s2m_2d()      : Transfers a 2D local array to the master          */
/* s2m_flux()    : Transfers a local flux array to the master        */
/* s2m_vel()     : Transfers a local velocity array to the master    */
/* s2c_2d()      : Sparse to Cartesian 2D transfers                  */
/* s2c_3d()      : Sparse to Cartesian 3D transfers                  */
/* c2s_2d()      : Cartesian to sparse 2D transfers                  */
/* c2s_3d()      : Cartesian to sparse 3D transfers                  */
/*-------------------------------------------------------------------*/


/*
 * The following slaves update routine have been split from the
 * various *_fill/*_refill routines so that they can be additionally
 * called from the TransferBuffers methods
 */
void win_data_slave_update_fill_3d(master_t *master, geometry_t *window,
				   window_t *windat, int nwindows)
{
  int c, cc, lc;                /* Local sparse coordinate / counter */
  int tn, tt;                   /* Tracer counter                    */
  int s, ce1, ce2;

  /*-----------------------------------------------------------------*/
  /* Variables to transfer for one or multiple windows               */
  windat->t = master->t;
  windat->days = master->days;
  windat->dt = master->dt;
  windat->dtf = master->dtf;
  windat->dtb = master->dtb;
  windat->dttr = master->dttr;
  windat->rampval = master->rampval;
  windat->nstep = master->nstep;
  /*windat->etarlxtc = master->etarlxtc;*/

  windat->df_diagn_set = master->df_diagn_set;
  tr_diag_reset_m(master);
  copy_dens_profile(master, window, windat);

  for (tt = 0; tt < master->ntr; tt++)
    window->wincon->trinfo_3d[tt].flag = master->trinfo_3d[tt].flag;

  /*-----------------------------------------------------------------*/
  /* Reset the windows if required                                   */
  if (master->regf & RS_RESET)
    window_reset(master, window, windat, window->wincon, RS_VH);
  if (nwindows == 1) {
    /* Transfer tracer boundary data                                 */
    bdry_transfer_tr(master, window, windat);
    /* Transfer sourcesink fluxes                                    */
    for (s = 0; s < master->npss; s++) {
      windat->wflux[s] = master->wflux[s];
      for (tt = 0; tt < master->ntr; tt++) {
        windat->tflux[tt][s] = master->tflux[tt][s];
      }
    }
    return;
  }

  /*-----------------------------------------------------------------*/
  /*-----------------------------------------------------------------*/
  /* Wet + auxiliary cells. This include variables that the master   */
  /* updates.                                                        */
  for (tt = 0; tt < master->nrlx; tt++) {
    tn = master->relax[tt];
    window->wincon->trinfo_3d[tn].relax_rate = 
      master->trinfo_3d[tn].relax_rate;
  }

  /*-----------------------------------------------------------------*/
  /* Include tracers since tracers may be altered by relaxation      */
  /* or resetting operations.                                        */
  for (cc = 1; cc <= window->b3_t; cc++) {
    lc = window->w3_t[cc];
    c = window->wsa[lc];
    for (tt = 0; tt < master->nrlx; tt++) {
      tn = master->relax[tt];
      windat->tr_wc[tn][lc] = master->tr_wc[tn][c];
    }
    for (tt = 0; tt < master->nres; tt++) {
      tn = master->reset[tt];
      windat->tr_wc[tn][lc] = master->tr_wc[tn][c];
      if (master->swr_attn == master->tr_wc[tn]) 
	windat->swr_attn[lc] = master->swr_attn[c];
    }
    if (master->rtemp) windat->rtemp[lc] = master->rtemp[c];
    if (master->rsalt) windat->rsalt[lc] = master->rsalt[c];
    if (!(master->decf & (NONE|DEC_ETA))) {
      windat->decv1[lc] = master->decv1[c];
      windat->decv2[lc] = master->decv2[c];
    }
    for (tt = 0; tt < master->ndhw; tt++) {
      if (master->dhwf[tt] & DHW_NOAA) {
	windat->dhw[tt][lc] = master->dhw[tt][c];
	windat->dhd[tt][lc] = master->dhd[tt][c];
      }
    }
  }
  for (cc = 1; cc <= window->enonS; cc++) {
    c = window->wsa[cc];
    windat->wind1[cc] = master->wind1[c];
    windat->wind2[cc] = master->wind2[c];
    windat->patm[cc] = master->patm[c];
    windat->windspeed[cc] = master->windspeed[c];
    windat->winddir[cc] = master->winddir[c];
    for (tt = 0; tt < master->nres2d; tt++) {
      tn = master->reset2d[tt];
      windat->tr_wcS[tn][cc] = master->tr_wcS[tn][c];
    }
    if (master->orbital) {
      windat->wave_ub[cc] = master->wave_ub[c];
      windat->wave_period[cc] = master->wave_period[c];
      windat->wave_dir[cc] = master->wave_dir[c];
      windat->wave_amp[cc] = master->wave_amp[c];
      windat->ustrcw[cc] = master->ustrcw[c];
    }
    if (master->waves & TAN_RAD) {
      windat->wave_Sxy[cc] = master->wave_Sxy[c];
      windat->wave_Syx[cc] = master->wave_Syx[c];
    }
    if (master->waves & WAVE_FOR) {
      windat->wave_Fx[cc] = master->wave_Fx[c];
      windat->wave_Fy[cc] = master->wave_Fy[c];
    }
    if (master->waves & STOKES) {
      windat->wave_ste1[cc] = master->wave_ste1[c];
      windat->wave_ste2[cc] = master->wave_ste2[c];
      if (windat->tau_w1) windat->tau_w1[cc] = master->tau_w1[c];
      if (windat->tau_w2) windat->tau_w2[cc] = master->tau_w2[c];
      if (windat->tau_diss1) windat->tau_diss1[cc] = master->tau_diss1[c];
      if (windat->tau_diss2) windat->tau_diss2[cc] = master->tau_diss2[c];
    }
    if (master->waves & SPECTRAL) {
      for (tt = 0; tt < master->nsfr; tt++) {
	double f = 2.0 * PI * master->freq[tt];
        windat->freq[tt] = 2.0 * (f * f)/ master->g;
      }
    }
  }
  
  /*-----------------------------------------------------------------*/
  /* Wet bulb and dew point                                          */
  if (master->sh_f & (WETBULB|DEWPOINT)) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->wetb[cc] = master->wetb[c];
    }
  }
    
  /*-----------------------------------------------------------------*/
  /* Relative humidity                                               */
  if (master->sh_f & (RELHUM|SPECHUM))
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->rh[cc] = master->rh[c];
    }

  /*-----------------------------------------------------------------*/
  /* Heatflux variables                                              */
  if (master->heatflux & NET_HEAT) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->swr[cc]   = master->swr[c];
      windat->heatf[cc] = master->heatf[c];
    }
  }
  if (master->heatflux & ADVANCED) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      if (fabs(window->wincon->albedo) <= 1.0)
	windat->swr[cc] = master->swr[c];
      windat->airtemp[cc] = master->airtemp[c];
      if (master->cloud) windat->cloud[cc] = master->cloud[c];
      if (master->heatflux & COMP_LWI)
	windat->lwri[cc] = master->lwri[c];
    }
  }
  if (master->heatflux & COMP_HEAT) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->swr[cc]  = master->swr[c];
      if (master->airtemp)
	windat->airtemp[cc] = master->airtemp[c];
      windat->lwrd[cc] = master->lwrd[c];
      windat->lwro[cc] = master->lwro[c];
      windat->shfd[cc] = master->shfd[c];
      windat->lhfd[cc] = master->lhfd[c];
    }
  }
  if (master->heatflux & (COMP_HEAT_MOM | COMP_HEAT_NONE)) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      // The if's are needed for COMP_HEAT_NONE
      if (master->swrd) windat->swrd[cc] = master->swrd[c];
      if (master->lwrd) windat->lwro[cc] = master->lwro[c];
      if (master->shfd) windat->shfd[cc] = master->shfd[c];
      if (master->lhfd) windat->lhfd[cc] = master->lhfd[c];
      // precipitation and evaporation
      if (master->precip) windat->precip[cc] = master->precip[c];
      if (master->evap)   windat->evap[cc]   = master->evap[c];
    }
  }
  if (master->heatflux & (SURF_RELAX|INVERSE)) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->hftemp[cc] = master->hftemp[c];
    }
  }
  if (master->albedo_l >= 0.0) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->light[cc] = master->light[c];
    }
  }
  /*-----------------------------------------------------------------*/
  /* Saltflux variables                                              */
  if (master->saltflux & (ADVANCED | ORIGINAL | BULK)) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->lhfd[cc] = master->lhfd[c];
      if (master->precip) windat->precip[cc] = master->precip[c];
    }
  }
  /*-----------------------------------------------------------------*/
  /* Source and sink variables                                       */
  if (master->npss) {
    for (s = 0; s < window->wincon->npss; s++) {
      pss_t *p = &window->wincon->pss[s];
      windat->wflux[s] = master->wflux[p->s2m];
      for (tt = 0; tt < master->ntr; tt++)
        windat->tflux[tt][s] = master->tflux[tt][p->s2m];
    }
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->waterss2d[cc] = master->waterss2d[c];
    }
    for (cc = 1; cc <= window->enon; cc++) {
      c = window->wsa[cc];
      windat->waterss[cc] = master->waterss[c];
    }
  }
  /*-----------------------------------------------------------------*/
  /* Surface relaxation                                              */
  if (master->etarlx & (RELAX|ALERT|BOUNDARY)) {
    relax_fill(window, master->eta_rlx, windat->eta_rlx);
    /*
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->eta_rlx[cc] = master->eta_rlx[c];
    }
    */
  }
  if (master->velrlx & RELAX)
    relax_fill(window, master->vel_rlx, windat->vel_rlx);

  /*-----------------------------------------------------------------*/
  /* Means                                                           */
  /* Mean iteration counter                                          */
  if (!(master->means & NONE)) {
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      windat->meanc[cc] = master->meanc[c];
    }
  }
  /* Mean variables; transfer only if the master was re-initialized  */
  if (master->means & RESET) {
    for (cc = 1; cc <= window->enon; cc++) {
      c = window->wsa[cc];
      for (tt = 0; tt < master->ntm_3d; tt++) {
	tn = master->tm_3d[tt];
	windat->tr_wc[tn][cc] = master->tr_wc[tn][c];
      }
    }
    for (cc = 1; cc <= window->enonS; cc++) {
      c = window->wsa[cc];
      for (tt = 0; tt < master->ntm_2d; tt++) {
	tn = master->tm_2d[tt];
	windat->tr_wcS[tn][cc] = master->tr_wcS[tn][c];
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the auxiliary multi-dt cells at the start of the timestep   */
  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    c = window->wsa[c];
    for (tn = 0; tn < windat->ntr; tn++)
      windat->tr_wc_as[tn][cc] = master->tr_wc[tn][c];
  }

  /*-----------------------------------------------------------------*/
  /* Transfer any custom data from the master to the slaves          */
  bdry_transfer_tr(master, window, windat);

  /*-----------------------------------------------------------------*/
  /* Copy the tracer increment flags to the window                   */
  memcpy(windat->trinc, master->trinc, master->ntr * sizeof(int));
  memcpy(windat->trincS, master->trincS, master->ntrS * sizeof(int));

}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master. If only one window exists the master points to the windat */
/* data structure - no need to fill.                                 */
/* 2D arrays are filled in win_data_fill_2d() which is invoked in    */
/* each iteration of the 2D mode loop.                               */
/*-------------------------------------------------------------------*/
void win_data_fill_3d(master_t *master,   /* Master data             */
                      geometry_t *window, /* Window geometry         */
                      window_t *windat,   /* Window data             */
                      int nwindows        /* Number of windows       */
  )
{
  int c, cc, lc;                /* Local sparse coordinate / counter */
  int tn, tt;                   /* Tracer counter                    */
  int s, ce1, ce2;

  /* Update from master */
  win_data_slave_update_fill_3d(master, window, windat, nwindows);

  /* No need for transfers if single window */
  if (nwindows == 1)
    return;

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require      */
  /* auxiliary cells only to be transfered to the slave since the    */
  /* wet cells already exist on the slave.                           */
  /* 3D variables.                                                   */
  for (cc = 1; cc <= window->nm2s; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    ce1 = window->m2se1[cc];
    ce2 = window->m2se2[cc];
    windat->u1[lc] = master->u1[ce1];
    windat->u2[lc] = master->u2[ce2];
    windat->w[lc] = master->w[c];
    windat->Kz[lc] = master->Kz[c];
    windat->Vz[lc] = master->Vz[c];
    windat->dens[lc] = master->dens[c];
    for (tn = 0; tn < windat->ntr; tn++) {
      windat->tr_wc[tn][lc] = master->tr_wc[tn][c];
    }
    windat->u1b[lc] = master->u1b[ce1];
    windat->u2b[lc] = master->u2b[ce2];
  }
  /* These fluxes do not need to be transferred unless diagnostic    */
  /* checks are required for continuity across windows.              */

  zflux_e1(geom,window,master->u1flux3d,windat->u1flux3d,window->m2s,
	   window->m2se1,window->nm2s);
  zflux_e2(geom,window,master->u2flux3d,windat->u2flux3d,window->m2s,
	   window->m2se2,window->nm2s); 
  /* Default update is a copy: alternative updates for zoomed grids  */
  /* are applied here.                                               */
  /*
  if (master->dozoom != NOZOOM) {
    zvel_filter(geom,window,master->u1,windat->u1,window->m2s,
		window->zmee2, window->m2se1,window->nm2s);
    zvel_filter(geom,window,master->u2,windat->u2,window->m2s,
		window->zmee1, window->m2se1,window->nm2s);
  }
  */

  /*-----------------------------------------------------------------*/
  /* 2D variables.                                                   */
  for (cc = 1; cc <= window->nm2sS; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    ce1 = window->m2se1[cc];
    ce2 = window->m2se2[cc];
    windat->wtop[lc] = master->wtop[c];
    windat->wbot[lc] = master->wbot[c];
    windat->u1bot[lc] = master->u1bot[ce1];
    windat->u2bot[lc] = master->u2bot[ce2];
    /* The updated value of eta is required at the auxiliary cells */
    /* to set dzu1 and dzu2, used in the calculation of the */
    /* advection terms.  */
    windat->eta[lc] = master->eta[c];
  }
  if (master->decf & DEC_ETA) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      lc = window->w2_t[cc];
      c = window->wsa[lc];
      windat->decv1[lc] = master->decv1[c];
      windat->decv2[lc] = master->decv2[c];
    }
  }

}

/* END win_data_fill_3d()                                            */
/*-------------------------------------------------------------------*/

void win_data_slave_update_fill_2d(master_t *master, geometry_t *window,
				   window_t *windat, int nwindows)
{
  int c, cc, lc, ce1, ce2;      /* Local sparse coordinate / counter */

  windat->dtb2 = windat->dtf2;
  windat->dtf2 = master->dt2d;
  windat->dt2d = windat->dtf2 + windat->dtb2;
  /* windat->dt2d=2.0*master->dt2d; *//* Leapfrog timestep */

  if (master->mode2d) {
    windat->dtu1 = master->dtu1;
    windat->dtu2 = master->dtu2;
  }
  if (nwindows == 1) {
    bdry_transfer_eta(master, window, windat);
    bdry_transfer_u1av(master, window, windat);
    bdry_transfer_u2av(master, window, windat);
    return;
  }

  /* These fluxes do not need to be transferred unless diagnostic */
  /* checks are required for continuity across windows.  */
  /*
  zflux_e1(geom,window,master->u1flux,windat->u1flux,window->m2s,
	   window->m2se1,window->nm2sS);
  zflux_e2(geom,window,master->u2flux,windat->u2flux,window->m2s,
	   window->m2se2,window->nm2sS); 
  */

  /*-----------------------------------------------------------------*/
  /* Set the auxiliary multi-dt cells at the start of the timestep */
  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    if (c < window->enonS) {
      c = window->wsa[c];
      windat->u1av_as[cc] = master->u1av[c];
      windat->u2av_as[cc] = master->u2av[c];
    }
  }

  /* Transfer any custom data from the master to the slaves */
  bdry_transfer_eta(master, window, windat);
  bdry_transfer_u1av(master, window, windat);
  bdry_transfer_u2av(master, window, windat);

}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master (auxiliary cells only) for the 2D mode. If only one window */
/* exists the master points to the windat data structure - no need   */
/* to fill.                                                          */
/*-------------------------------------------------------------------*/
void win_data_fill_2d(master_t *master, /* Master data structure */
                      geometry_t *window, /* Window structure */
                      window_t *windat, /* Window data structure */
                      int nwindows  /* Number of windows */
  )
{
  int c, cc, lc, ce1, ce2;      /* Local sparse coordinate / counter */

  /* Update from master */
  win_data_slave_update_fill_2d(master, window, windat, nwindows);
  
  if (nwindows == 1)
    return;

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require */
  /* auxiliary cells only to be transfered to the slave since the */
  /* wet cells already exist on the slave.  */
  for (cc = 1; cc <= window->nm2sS; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    ce1 = window->m2se1[cc];
    ce2 = window->m2se2[cc];
    windat->eta[lc] = master->eta[c];
    windat->detadt[lc] = master->detadt[c];
    windat->u1av[lc] = master->u1av[ce1];
    windat->u2av[lc] = master->u2av[ce2];
    windat->depth_e1[lc] = master->depth_e1[ce1];
    windat->depth_e2[lc] = master->depth_e2[ce2];
    windat->etab[lc] = master->etab[c];
    /* Uncomment to check master to slave transfer indicies          */
    /*
    printf("m2se1 %d window%d (%d %d)m->(%d %d)s\n",cc,window->wn,geom->s2i[ce1],
	   geom->s2j[ce1],window->s2i[lc],window->s2j[lc]);
    printf("m2se2 window%d (%d %d)m->(%d %d)s\n",window->wn,geom->s2i[ce2],
	   geom->s2j[ce2],window->s2i[lc],window->s2j[lc]);
    */

    /*
    windat->u1avb[lc]=master->u1avb[ce1];
    windat->u2avb[lc]=master->u2avb[ce2]; 
    */
  }

}

/* END win_data_fill_2d()                                            */
/*-------------------------------------------------------------------*/


void win_data_slave_update_refill_3d(master_t *master, geometry_t *window,
				     window_t *windat, int nwindows,
				     int mode)
{
  /* Set the crash recovery flags if required */
  if (master->crf & RS_WINSET) master->crf &= ~RS_WINSET;
  if (master->crf & RS_RESET) master->crf = NONE;
  
  if (MIXING) {
    /* Transfer any custom data from the master to the slaves */
    bdry_transfer_u1(master, window, windat);
    bdry_transfer_u2(master, window, windat);
  }
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master required to perform the u2av velocity update.              */
/*-------------------------------------------------------------------*/
void win_data_refill_3d(master_t *master, /* Master data structure */
                        geometry_t *window, /* Window structure */
                        window_t *windat, /* Window data structure */
                        int nwindows, /* Number of windows */
                        int mode  /* Flag to set data to empty */
  )
{
  int c, cc, lc, ce1, ce2;      /* Local sparse coordinate / counter */

  /* Update from master */
  win_data_slave_update_refill_3d(master, window, windat, nwindows, mode);

  if (nwindows == 1 && mode & MIXING)
    return;
  
  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require */
  /* auxiliary cells only to be transfered to the slave since the */
  /* wet cells already exist on the slave.  */
  /* mode = MIXING : vertical mixing coefficients and cell */
  /* thicknesses at the e1 and e2 faces.  */
  if (mode & MIXING) {
    geometry_t *geom = master->geom;
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->Kz[lc] = master->Kz[c];
      windat->Vz[lc] = master->Vz[c];
      windat->dzu1[lc] = master->dzu1[c];
      windat->dzu2[lc] = master->dzu2[c];
    }
    if (windat->tke && windat->diss) {
      for (cc = 1; cc <= window->nm2s; cc++) {
	lc = window->m2s[cc];
	c = window->wsa[lc];
	windat->tke[lc] = master->tke[c];
	windat->diss[lc] = master->diss[c];
      }
    }
    if (windat->tke && windat->omega) {
      for (cc = 1; cc <= window->nm2s; cc++) {
	lc = window->m2s[cc];
	c = window->wsa[lc];
	windat->tke[lc] = master->tke[c];
	windat->omega[lc] = master->omega[c];
      }
    }
    if (windat->Q2 && windat->Q2L) {
      for (cc = 1; cc <= window->nm2s; cc++) {
	lc = window->m2s[cc];
	c = window->wsa[lc];
	windat->Q2[lc] = master->Q2[c];
	windat->Q2L[lc] = master->Q2L[c];
      }
    }
    if (windat->sdc) {
      for (cc = 1; cc <= window->nm2s; cc++) {
	lc = window->m2s[cc];
	c = window->wsa[lc];
	windat->sdc[lc] = master->sdc[c];
      }
    }
    /* 2D */
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->sur_e1[lc] = geom->sur_e1[c];
      windat->sur_e2[lc] = geom->sur_e2[c];
    }
  }
  /* mode = VELOCITY : depth averaged adjusted velocities for the */
  /* vertical velocity equation and 3D fluxes for the tracer */
  /* equation.  */
  if (mode & VELOCITY) {
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      ce1 = window->m2se1[cc];
      ce2 = window->m2se2[cc];
      windat->u1[lc] = master->u1[ce1];
      windat->u2[lc] = master->u2[ce2];
      /* Cell thickness is required for tracer horizontal diffusion */
      windat->dzu1[lc] = master->dzu1[c];
      windat->dzu2[lc] = master->dzu2[c];
    }
    zflux_e1(geom, window, master->u1flux3d, windat->u1flux3d, window->m2s,
             window->m2se1, window->nm2s);
    zflux_e2(geom, window, master->u2flux3d, windat->u2flux3d, window->m2s,
             window->m2se2, window->nm2s);
    /* Default update is a copy: alternative updates for zoomed     */
    /* grids are applied here.                                      */
    /*
    if (master->dozoom != NOZOOM) {
      zvel_filter(geom,window,master->u1,windat->u1,window->m2s,
		  window->zmee2, window->m2se1,window->nm2s);
      zvel_filter(geom,window,master->u2,windat->u2,window->m2s,
		  window->zmee1, window->m2se1,window->nm2s);
    }
    */

    /* The updated value of eta is required at the auxiliary cells */
    /* to set wtop, used in the calculation of the advection terms.  */
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->eta[lc] = master->eta[c];
    }

    /* These fluxes do not need to be transferred unless diagnostic */
    /* checks are required for continuity across windows.  */
    /*
    zflux_e1(geom,window,master->u1flux,windat->u1flux,window->m2s,
	     window->m2se1,window->nm2sS);
    zflux_e2(geom,window,master->u2flux,windat->u2flux,window->m2s,
	     window->m2se2,window->nm2sS); 
    */
  }
  if (mode & TRACERS) {
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->dens[lc] = master->dens[c];
    }
  }
}

/* END win_data_refill_3d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master required to perform the u2av velocity update.              */
/*-------------------------------------------------------------------*/
void win_data_refill_2d(master_t *master, /* Master data structure */
                        geometry_t *window, /* Window structure */
                        window_t *windat, /* Window data structure */
                        int nwindows, /* Number of windows */
                        int mode  /* df_variable_t to transfer */
  )
{
  int cc, lc, ce1, ce2;         /* Local sparse coordinate / counter */

  /* mode = NVELOCITY : updated velocities. This is required to */
  /* maintain continuity by using the correct nuav when setting uav */
  /* in asselin() at auxiliary cells; uav is subsequently used to */
  /* compute fluxes in eta_step(). Note: if asselin() is performed */
  /* in leapfrog_update_2d() then this transfer is not neccessary but */
  /* the elevation calculation is made using unfiltered velocities.  */
  if (mode & NVELOCITY) {
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      ce1 = window->m2se1[cc];
      ce2 = window->m2se2[cc];
      windat->nu1av[lc] = master->nu1av[ce1];
      windat->nu2av[lc] = master->nu2av[ce2];
    }
    /* Transfer fine grid zoom fluxes summed over the coarse grid to */
    /* the slaves.                                                   */
    if (master->dozoom != NOZOOM) {
      zflux_e1(geom, window, master->u1flux_z, windat->u1flux_z, window->m2s,
	       window->m2se1, window->nm2sS);
      zflux_e2(geom, window, master->u2flux_z, windat->u2flux_z, window->m2s,
	       window->m2se2, window->nm2sS); 
    }
    /* Default update is a copy: alternative updates for zoomed     */
    /* grids are applied here.                                      */
    /*
    if (master->dozoom != NOZOOM) {
      zvel_filter(geom, window, master->nu1av, windat->nu1av, window->m2s,
		  window->zmee2, window->m2se1,window->nm2sS);
      zvel_filter(geom, window, master->nu2av, windat->nu2av, window->m2s,
		  window->zmee1, window->m2se1,window->nm2sS);
    }
    */
  }
  if (mode & DEPTH) {
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      ce1 = window->m2se1[cc];
      ce2 = window->m2se2[cc];
      windat->depth_e1[lc] = master->depth_e1[ce1];
      windat->depth_e2[lc] = master->depth_e2[ce2];
    }
  }
  if (mode & ETA_A) {
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      ce1 = window->wsa[lc];
      windat->eta[lc] = master->eta[ce1];
    }
  }
}

/* END win_data_refill_2d()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to empty the window data structure back into the master   */
/* (cells which are used as auxiliary cells in other windows only)   */
/* for the 3D mode. If only one window exists the master points to   */
/* the windat data structure - no need to empty.                     */
/*-------------------------------------------------------------------*/
void win_data_empty_3d(master_t *master,   /* Master data            */
                       geometry_t *window, /* Window geometry        */
                       window_t *windat,   /* Window data            */
                       int mode         /* Flag to set data to empty */
  )
{
  int c, cc, lc;                /* Local sparse coordinate / counter */
  int tn, tt;                   /* Tracer counter                    */
  geometry_t *geom = master->geom;

  /* Set the CFL diagnostics                                         */
  if (mode & CFL) {
    if (!(master->cfl & NONE)) {
      if (windat->mcfl2d < master->mcfl2d)
        master->mcfl2d = windat->mcfl2d;
      if (windat->mcfl3d < master->mcfl3d) {
        master->mcfl3d = windat->mcfl3d;
        master->cflc = windat->cflc;
      }
    }
  }
  /* Set the totals                                                  */
  if (mode & TMASS) {
    if (master->totals) {
      int nb, wn;
      master->tmass += windat->tmass;
      master->tvol += windat->tvol;
      for (tn = 0; tn < master->ntot; tn++)
	master->trtot[tn] += windat->trtot[tn];
      for (nb = 0; nb < geom->nobc; nb++) {
	wn = geom->owc[nb][window->wn];
	if (wn >= 0)
	  master->vf[nb] += windat->vf[wn];
      }
    }
    if (geom->nregions)
      region_transfer(master, window);
  }
  /* Set tracer flags                                                */
  for (tt = 0; tt < master->ntr; tt++)
    master->trinfo_3d[tt].flag = window->wincon->trinfo_3d[tt].flag;

  if(window->nwindows == 1)
    return;

  /*-----------------------------------------------------------------*/
  /* Variables that require wet cells which are used by other        */
  /* windows to be transfered to the master.                         */
  /* mode = MIXING : vertical mixing coefficients and cell           */
  /* thicknesses at the e1 and e2 faces.                             */
  if (mode & MIXING) {
    geometry_t *geom = master->geom;
    /* Transfer to the master                                        */
    s2m_vel(master->dzu1, windat->dzu1,
            window->s2m, window->s2me1, window->ns2m);
    s2m_vel(master->dzu2, windat->dzu2,
            window->s2m, window->s2me2, window->ns2m);

    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->Kz[c] = windat->Kz[lc];
      master->Vz[c] = windat->Vz[lc];
    }
    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      geom->sur_e1[c] = windat->sur_e1[lc];
      geom->sur_e2[c] = windat->sur_e2[lc];
    }
    if (windat->tke && windat->diss) {
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->tke[c] = windat->tke[lc];
	master->diss[c] = windat->diss[lc];
      }
    }
    if (windat->tke && windat->omega) {
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->tke[c] = windat->tke[lc];
	master->omega[c] = windat->omega[lc];
      }
    }
    if (windat->Q2 && windat->Q2L) {
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->Q2[c] = windat->Q2[lc];
	master->Q2L[c] = windat->Q2L[lc];
      }
    }
    if (windat->sdc) {
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->sdc[c] = windat->sdc[lc];
      }
    }
    /* 
       if(!(master->cfl&NONE)) {
       if(windat->mcfl2d<master->mcfl2d)master->mcfl2d=windat->mcfl2d;
       if(windat->mcfl3d<master->mcfl3d) { master->mcfl3d=windat->mcfl3d;
       master->cflc=windat->cflc; } } */

  }
  /* mode = VELOCITY : depth averaged adjusted velocities for the */
  /* vertical velocity equation and 3D fluxes for the tracer */
  /* equation.  */
  if (mode & VELOCITY) {
    s2m_vel(master->u1, windat->u1,
            window->s2m, window->s2me1, window->ns2m);
    s2m_flux(geom, window, master->u1flux3d, windat->u1flux3d,
             window->s2m, window->s2me1, window->ns2m, 1);
    s2m_vel(master->u2, windat->u2,
            window->s2m, window->s2me2, window->ns2m);
    s2m_flux(geom, window, master->u2flux3d, windat->u2flux3d,
             window->s2m, window->s2me2, window->ns2m, 0);

    s2m_vel(master->u1b, windat->u1b,
            window->s2m, window->s2me1, window->ns2m);
    s2m_vel(master->u2b, windat->u2b,
            window->s2m, window->s2me2, window->ns2m);

    /* Set in set_new_cells_ and used in tracer horz diffusion */
    /* Potential bug: zoomed grids do not work with this transfer */
    /*if (master->zmfe1 <= 1)*/
      s2m_vel(master->dzu1, windat->dzu1,
	      window->s2m, window->s2me1, window->ns2m);
      /*if (master->zmfe2 <= 1)*/
      s2m_vel(master->dzu2, windat->dzu2,
	      window->s2m, window->s2me2, window->ns2m);
  }
  /* mode = WVEL : 3D vertical velocity and 2D surface and bottom */
  /* vertical velocities. Also include the 3D velocity deviations */
  /* from the depth averaged velocity which was calculated in */
  /* velocity_adjust() here.  */
  if (mode & WVEL) {
    /* Transfer to the master */
    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->w[c] = windat->w[lc];
    }
    s2m_vel(master->u1bot, windat->u1bot,
            window->s2m, window->s2me1, window->ns2mS);
    s2m_vel(master->u2bot, windat->u2bot,
            window->s2m, window->s2me2, window->ns2mS);
    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->wtop[c] = windat->wtop[lc];
      master->wbot[c] = windat->wbot[lc];
    }
  }
  /* mode = TRACERS : variables calculated in the tracer routine */
  if (mode & TRACERS) {
    master->wclk[window->wn] += windat->wclk;
    memcpy(master->trinc, windat->trinc, master->ntr * sizeof(int));
    if (master->means & RESET) master->means &= ~RESET;
    /* Transfer to the master */
    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->dens[c] = windat->dens[lc];
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
    }
    /*---------------------------------------------------------------*/
    /* All wet cells (cell centered). Include relaxation and reset */
    /* tracers since tracers may be altered by updates. Include */
    /* OBC cells.  */
    for (cc = 1; cc <= window->b3_t; cc++) {
      lc = window->w3_t[cc];
      c = window->wsa[lc];
      for (tt = 0; tt < master->nrlx; tt++) {
        tn = master->relax[tt];
        master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
      for (tt = 0; tt < master->nres; tt++) {
        tn = master->reset[tt];
        master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
      }
    }
    /* Not needed now that the heatfluxes are calculated on the windows */
    /*
      if (master->heatflux & ADVANCED) {
      for (cc = 1; cc <= window->b2_t; cc++) {
      lc = window->w2_t[cc];
      c = window->wsa[lc];
      master->tr_wc[master->tno][c] = windat->tr_wc[windat->tno][lc];
      master->tr_wc[master->sno][c] = windat->tr_wc[windat->sno][lc];
      master->dens[c] = windat->dens[lc]; // already transfered above
      }
      }
    */
    if (master->waves & TAN_RAD) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	lc = window->w2_t[cc];
	c = window->wsa[lc];
	master->wave_Sxy[c] = windat->wave_Sxy[lc];
	master->wave_Syx[c] = windat->wave_Syx[lc];
      }
    }
    if (master->waves & WAVE_FOR) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	lc = window->w2_t[cc];
	c = window->wsa[lc];
	master->wave_Fx[c] = windat->wave_Fx[lc];
	master->wave_Fy[c] = windat->wave_Fy[lc];
      }
    }
    if (master->waves & STOKES) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	lc = window->w2_t[cc];
	c = window->wsa[lc];
	master->wave_ste1[c] = windat->wave_ste1[lc];
	master->wave_ste2[c] = windat->wave_ste2[lc];
	if (master->tau_w1) master->tau_w1[c] = windat->tau_w1[lc];
	if (master->tau_w2) master->tau_w2[c] = windat->tau_w2[lc];
	if (master->tau_diss1) master->tau_diss1[c] = windat->tau_diss1[lc];
	if (master->tau_diss2) master->tau_diss2[c] = windat->tau_diss2[lc];
      }
    }
    for (tt = 0; tt < master->ndhw; tt++) {
      if (master->dhwf[tt] & DHW_NOAA) {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  lc = window->w3_t[cc];
	  c = window->wsa[lc];
	  master->dhd[tt][c] = windat->dhd[tt][lc];
	}
      }
    }
    if (!(master->decf & NONE)) {
      if (master->decf == DEC_ETA) {
	for (cc = 1; cc <= window->b2_t; cc++) {
	  lc = window->w2_t[cc];
	  c = window->wsa[lc];
	  master->decv[c] = windat->eta[lc];
	}
      } else {
	for (cc = 1; cc <= window->b3_t; cc++) {
	  lc = window->w3_t[cc];
	  c = window->wsa[lc];
	  if (master->decf == DEC_U1)
	    master->decv[c] = windat->u1[lc];
	  else if (master->decf == DEC_U2)
	    master->decv[c] = windat->u2[lc];
	  else
	    master->decv[c] = windat->tr_wc[master->decn][lc];
	}
      }
    }

    /*---------------------------------------------------------------*/
    /* Particle tracking variables                                   */
    if (master->ptn) {
      win_priv_t *wincon = window->wincon;
      for (cc = 1; cc <= window->b3_t; cc++) {
        lc = window->w3_t[cc];
        c = window->wsa[lc];
        master->Kz[c] = windat->Kz[lc];
        master->dz[c] = wincon->dz[lc];
      }
      for (cc = 1; cc <= window->b2_t; cc++) {
        lc = window->w2_t[cc];
        c = window->wsa[lc];
        master->eta[c] = windat->eta[lc];
      }
      /* Get the surface sparse coordinates */
      if (!master->sigma) {
        int zm1;
        for (cc = 1; cc <= geom->b2_t; cc++) {
          c = lc = geom->w2_t[cc];
          zm1 = geom->zm1[c];
          while (c != zm1 && geom->gridz[c] > master->eta[lc]) {
            c = zm1;
            zm1 = geom->zm1[c];
          }
          geom->sur_t[cc] = c;
        }
      }
    }
    /* Mean iteration counter */
    if (!(master->means & NONE)) {
      for (cc = 1; cc <= window->enonS; cc++) {
	c = window->wsa[cc];
	master->meanc[c] = windat->meanc[cc];
      }
    }
  }
}

/* END win_data_empty_3d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to empty the window data structure back into the master   */
/* (cells which are used as auxiliary cells in other windows only)   */
/* for the 2D mode. If only one window exists the master points to   */
/* the windat data structure - no need to empty.                     */
/*-------------------------------------------------------------------*/
void win_data_empty_2d(master_t *master,  /* Master data structure */
                       geometry_t *window,  /* Window structure */
                       window_t *windat,  /* Window data structure */
                       int mode /* Flag to set data to empty */
  )
{
  int c, cc, lc;                /* Local sparse coordinate / counter */

  /*-----------------------------------------------------------------*/
  /* Cell centered 2D arrays */
  /* Variables that require wet cells which are used by other */
  /* windows to be transfered to the master.  */
  /* mode = VELOCITY : updated velocity and elevation data */
  if (mode & VELOCITY) {
    s2m_vel(master->u1av, windat->u1av,
	    window->s2m, window->s2me1, window->ns2mS);
    s2m_vel(master->u2av, windat->u2av,
	    window->s2m, window->s2me2, window->ns2mS);

    /* These velocities do not need to be transferred since the */
    /* whole uav array is copied to nuavb, auxiliary cells included, */
    /* and auxiliary cells in uav were transferred at the start of */
    /* the loop.  */
    /*
    s2m_vel(master->u1avb,windat->u1avb,
	    window->s2m,window->s2me1,window->ns2mS);
    s2m_vel(master->u2avb,windat->u2avb,
	    window->s2m,window->s2me2,window->ns2mS); 
    */
    /* These fluxes do not need to be transferred unless diagnostic */
    /* checks are required for continuity across windows.  */
    /*
    s2m_flux(geom,window,master->u1flux,windat->u1flux,
	     window->s2m,window->s2me1,window->ns2mS,1);
    s2m_flux(geom,window,master->u2flux,windat->u2flux,
	     window->s2m,window->s2me2,window->ns2mS,0); 
    */

    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->eta[c] = windat->eta[lc];
      master->detadt[c] = windat->detadt[lc];
      master->etab[c] = windat->etab[lc];
    }

    if (master->dozoom != NOZOOM) {
      s2m_zoom(master, window, windat, ETA_A);
      s2m_zoom(master, window, windat, VELOCITY);
    }

    memcpy(master->trincS, windat->trincS, master->ntrS * sizeof(int));
  }
  /* mode = DEPTH : total depth at the cell faces */
  if (mode & DEPTH) {
    s2m_vel(master->depth_e1, windat->depth_e1,
            window->s2m, window->s2me1, window->ns2mS);
    s2m_vel(master->depth_e2, windat->depth_e2,
            window->s2m, window->s2me2, window->ns2mS);
  }
  /* mode = NVELOCITY : updated velocity. Used in asselin().  */
  if (mode & NVELOCITY) {
    s2m_vel(master->nu1av, windat->nu1av,
	    window->s2m, window->s2me1, window->ns2mS);
    s2m_vel(master->nu2av, windat->nu2av,
	    window->s2m, window->s2me2, window->ns2mS);

    /* Transfer fine grid zoom fluxes to the master. These fluxes    */
    /* are summed over the corse grid when transferred to the coarse */
    /* grid auxiliary cells from the master.                         */
    if (master->dozoom != NOZOOM) {
      s2m_zoom(master, window, windat, mode);
      s2m_flux(geom,window,master->u1flux_z,windat->u1flux_z,
	       window->s2m,window->s2me1,window->ns2mS,1);
      s2m_flux(geom,window,master->u2flux_z,windat->u2flux_z,
	       window->s2m,window->s2me2,window->ns2mS,0); 
    }
  }
  if (mode & ETA_A) {
    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->eta[c] = windat->eta[lc];
    }
  }

  /* Uncomment to check slave to master transfer indicies            */
  /*
  for (cc=1; cc <= window->ns2mS; cc++) {
    c = window->s2me1[cc];
    lc = window->s2m[cc];
    printf("s2me1 window%d (%d %d)s->(%d %d)m\n",window->wn,window->s2i[lc], 
	   window->s2j[lc],geom->s2i[c],geom->s2j[c]);
    c = window->s2me2[cc];
    printf("s2me2 window%d (%d %d)s->(%d %d)m\n",window->wn,window->s2i[lc], 
	   window->s2j[lc],geom->s2i[c],geom->s2j[c]);
  }
  */
}

/* END win_data_empty_2d()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Transfers the maximum alert diagnostics from each window to the   */
/* master.                                                           */
/*-------------------------------------------------------------------*/
void master_alert_fill(master_t *master,    /* Master data           */
		       geometry_t *window,  /* Window geometry       */
		       window_t *windat     /* Window data           */
		       )
{
  int wn, nb;

  /* Summary */
  memset(master->nalert, 0, nna * sizeof(int));
  for (nb = 0; nb < nna; nb++)
    master->nalert[nb] += windat->nalert[nb];

  /* Maximum energy                                                  */
  master->em += windat->em;
  master->me += windat->me;
  master->ke += windat->ke;

  for (nb = 0; nb < geom->nobc; nb++) {
    wn = geom->owc[nb][window->wn];
    if (wn >= 0)
      master->ef[nb] += windat->ef[wn];
  }

  /* Maximum 3D velocities                                           */
  if (windat->mu1 > master->mu1) {
    master->mu1 = windat->mu1;
    master->cu1 = window->wsa[windat->cu1[0]];
  }
  if (windat->mu2 > master->mu2) {
    master->mu2 = windat->mu2;
    master->cu2 = window->wsa[windat->cu2[0]];
  }
  if (windat->mw > master->mw) {
    master->mw = windat->mw;
    master->cw = window->wsa[windat->cw];
  }
  /* Maximum 2D velocities and elevation                             */
  if (windat->meta > master->meta) {
    master->meta = windat->meta;
    master->ceta = window->wsa[windat->ceta[0]];
  }
  if (windat->mu1a  > master->mu1a) {
    master->mu1a = windat->mu1a;
    master->cu1a = window->wsa[windat->cu1a[0]];
  }
  if (windat->mu2a  > master->mu2a) {
    master->mu2a = windat->mu2a;
    master->cu2a = window->wsa[windat->cu2a[0]];
  }
  /* Maximum divergences                                             */
  if (windat->mdw > master->mdw) {
    master->mdw = windat->mdw;
    master->cdw = window->wsa[windat->cdw];
  }
  if (windat->mdeta  > master->mdeta) {
    master->mdeta = windat->mdeta;
    master->cdeta = window->wsa[windat->cdeta];
  }
  /* Minimum and maximum T/S                                         */
  if (windat->msmin  < master->msmin) {
    master->msmin = windat->msmin;
    master->csmin = window->wsa[windat->csmin];
  }
  if (windat->msmax  > master->msmax) {
    master->msmax = windat->msmax;
    master->csmax = window->wsa[windat->csmax];
  }
  if (windat->mtmin  < master->mtmin) {
    master->mtmin = windat->mtmin;
    master->ctmin = window->wsa[windat->ctmin];
  }
  if (windat->mtmax  > master->mtmax) {
    master->mtmax = windat->mtmax;
    master->ctmax = window->wsa[windat->ctmax];
  }
  /* Maximum rates of change of T/S                                  */
  if (windat->mdt  > master->mdt) {
    master->mdt = windat->mdt;
    master->cdt = window->wsa[windat->cdt];
  }
  if (windat->mds  > master->mds) {
    master->mds = windat->mds;
    master->cds = window->wsa[windat->cds];
  }
  /* CFL                                                             */
  if (windat->mcfl < master->mcfl) {
    master->mcfl = windat->mcfl;
    master->ccfl = window->wsa[windat->ccfl];
  }
  /* Shear                                                           */
  if (windat->mshear > master->mshear) {
    master->mshear = windat->mshear;
  }
  /* Mean eta mean                                                   */
  master->memean = windat->memean;
  /* Maximum tendencies                                              */
  if (windat->ma1 > master->ma1) {
    master->ma1 = windat->ma1;
    master->ca1 = window->wsa[windat->ca1];
  }
  if (windat->mh1 > master->mh1) {
    master->mh1 = windat->mh1;
    master->ch1 = window->wsa[windat->ch1];
  }
  if (windat->mv1 > master->mv1) {
    master->mv1 = windat->mv1;
    master->cv1 = window->wsa[windat->cv1];
  }
  if (windat->mb1 > master->mb1) {
    master->mb1 = windat->mb1;
    master->cb1 = window->wsa[windat->cb1];
  }
  if (windat->md1 > master->md1) {
    master->md1 = windat->md1;
    master->cd1 = window->wsa[windat->cd1];
  }
  if (windat->mc1 > master->mc1) {
    master->mc1 = windat->mc1;
    master->cc1 = window->wsa[windat->cc1];
  }
  if (windat->ma2 > master->ma2) {
    master->ma2 = windat->ma2;
    master->ca2 = window->wsa[windat->ca2];
  }
  if (windat->mh2 > master->mh2) {
    master->mh2 = windat->mh2;
    master->ch2 = window->wsa[windat->ch2];
  }
  if (windat->mv2 > master->mv2) {
    master->mv2 = windat->mv2;
    master->cv2 = window->wsa[windat->cv2];
  }
  if (windat->mb2 > master->mb2) {
    master->mb2 = windat->mb2;
    master->cb2 = window->wsa[windat->cb2];
  }
  if (windat->md2 > master->md2) {
    master->md2 = windat->md2;
    master->cd2 = window->wsa[windat->cd2];
  }
  if (windat->mc2 > master->mc2) {
    master->mc2 = windat->mc2;
    master->cc2 = window->wsa[windat->cc2];
  }
}

/* END master_alert_fill()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to generate vectors for master-slave data transfers       */
/*-------------------------------------------------------------------*/
void build_transfer_maps(geometry_t *geom,    /* Global geometry     */
			 geometry_t **window, /* Local geometry      */
                         int wn,              /* Window number       */
                         int nwindows         /* Number of windows   */
  )
{
  int c, cc, cl, ce1, ce2;      /* Sparse counters */
  int *wm;                      /* Window / ghost cell mask */
  int n;                        /* Window counter */
  int use_ghosts = 1;           /* Flag to include lateral ghosts */
  int zc, zse1, zse2;           /* Zoom counter / factor */

  wm = i_alloc_1d(geom->sgsiz);
  for (c = 1; c <= geom->sgnum; c++)
    wm[c] = geom->fm[c].wn;

  /* If lateral boundary conditions are set in the window then don't */
  /* include the ghost cells in the master-slave transfer vectors.  */
  /* The transfers will operate OK if they are included (with some */
  /* loss of speed) except if the window is zoomed where a variable */
  /* is not interpolated on the master at non-center locations, then */
  /* the lateral boundary condition is set to the non-interpolated */
  /* value and subsequently transferred back to the window. This */
  /* results in a invalid ghost value in the window.  */
#if !GLOB_BC
  use_ghosts = 0;
#endif

  /*-----------------------------------------------------------------*/
  /* Set the ghost cells in the window mask velocity ghost cells */
  /* Western edge ghost cells for e1 cells to process. These cells */
  /* are actually wet and defined to be in a window but due to the */
  /* stagger for e1 velocities at western boundaries they are ghost */
  /* cells for u1. These are identified if the west[c] is non-zero */
  /* (i.e. a wet cell in a window) and west[west[c]] is zero (i.e. a */
  /* ghost cell).  */
  for (cc = 1; cc <= geom->v3_e1; cc++) {
    c = geom->w3_e1[cc];        /* Global cell to process */
    cl = geom->xm1[c];          /* Global cell to the west of c */
    if (geom->fm[cl].wn && !(geom->fm[geom->xm1[cl]].wn))
      wm[cl] = 0;
  }
  /* Southern edge ghost cells for e2 cells to process.  */
  for (cc = 1; cc <= geom->v3_e2; cc++) {
    c = geom->w3_e2[cc];        /* Global cell to process */
    cl = geom->ym1[c];          /* Global cell to the south of c */
    if (geom->fm[cl].wn && !(geom->fm[geom->ym1[cl]].wn))
      wm[cl] = 0;
  }

  /* Count the number of auxiliary cells in window wn. Note : if */
  /* GLOB_BC=1 then ghost cells are also required to be included */
  /* as their values are set on the master.  */
  window[wn]->nm2s = window[wn]->nm2sS = 0;
  for (cc = 1; cc <= window[wn]->enon; cc++) {
    cl = window[wn]->m2d[cc];
    if (window[wn]->zoomc[cl] & (ZC | ZE1 | ZE2)) {
      c = window[wn]->wsa[cc];
      if (geom->fm[c].wn != wn) {
        /* Continue if use_ghosts=0 and fm[c].wn=0 (c=ghost cell) */
        if (!use_ghosts && !geom->fm[c].wn)
          continue;
        window[wn]->nm2s += 1;
        if (cc <= window[wn]->enonS)
          window[wn]->nm2sS += 1;
      }
    }
  }

  window[wn]->m2s = i_alloc_1d(window[wn]->nm2s + 1);

  /* Fill the master to slave transfer vector */
  window[wn]->nm2s = 1;
  for (cc = 1; cc <= window[wn]->enon; cc++) {
    cl = window[wn]->m2d[cc];
    if (window[wn]->zoomc[cl] & (ZC | ZE1 | ZE2)) {
      c = window[wn]->wsa[cc];
      if (geom->fm[c].wn != wn) {
        if (!use_ghosts && !geom->fm[c].wn)
          continue;
        window[wn]->m2s[window[wn]->nm2s] = cc;
        window[wn]->nm2s += 1;
      }
    }
  }
  window[wn]->nm2s -= 1;

  /* Get the e1 and e2 faces for master to slave transfers */
  if (window[wn]->zoom == DOZOOM) {
    window[wn]->m2se1 = i_alloc_1d(window[wn]->nm2s + 1);
    window[wn]->m2se2 = i_alloc_1d(window[wn]->nm2s + 1);
    zse1 = (int)(window[wn]->zmfe1 / 2);
    zse2 = (int)(window[wn]->zmfe2 / 2);
    for (cc = 1; cc <= window[wn]->nm2s; cc++) {
      c = window[wn]->m2s[cc];
      ce1 = ce2 = window[wn]->wsa[c];
      for (zc = 1; zc <= zse1; zc++)
        ce1 = geom->xm1[ce1];
      for (zc = 1; zc <= zse2; zc++)
        ce2 = geom->ym1[ce2];
      window[wn]->m2se1[cc] = ce1;
      window[wn]->m2se2[cc] = ce2;
    }
  } else {
    window[wn]->m2se1 = i_alloc_1d(window[wn]->nm2s + 1);
    window[wn]->m2se2 = i_alloc_1d(window[wn]->nm2s + 1);
    for (cc = 1; cc <= window[wn]->nm2s; cc++) {
      c = window[wn]->m2s[cc];
      window[wn]->m2se1[cc] = window[wn]->wsa[c];
      window[wn]->m2se2[cc] = window[wn]->wsa[c];
    }
    /* 
       window[wn]->m2se1=window[wn]->m2s;
       window[wn]->m2se2=window[wn]->m2s; */
  }

  /* Count the number of wet cells in window wn which are auxiliary */
  /* in some other window. If lateral BC's are set on the master */
  /* then the transfer vectors must include may values at the cell */
  /* interior to the ghost cell (e.g. free slip condition on */
  /* tangential velocity components).  */
  window[wn]->ns2m = window[wn]->ns2mS = 0;
  for (n = 1; n <= nwindows; n++) {
    if (n == wn)
      continue;
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      if (geom->fm[c].wn == wn) {
        cl = geom->fm[c].sc;
        window[wn]->ns2m += 1;
        if (cl <= window[wn]->enonS) {
          window[wn]->ns2mS += 1;
        }
      }
    }
  }

  if (use_ghosts) {
    for (cc = 1; cc <= geom->nbpt; cc++) {
      c = geom->bin[cc];
      if (geom->zoomc[geom->m2d[c]] & (ZC | ZE1 | ZE2) &&
          geom->fm[c].wn == wn) {
        /* if(geom->fm[c].wn == wn) { */
        cl = geom->fm[c].sc;
        window[wn]->ns2m += 1;
        if (cl <= window[wn]->enonS) {
          window[wn]->ns2mS += 1;
        }
      }
    }
  }
  window[wn]->s2m = i_alloc_1d(window[wn]->ns2m + 1);

  /* Fill the slave to master transfer vector */
  window[wn]->ns2m = window[wn]->ns2mS + 1;
  window[wn]->ns2mS = 1;
  for (n = 1; n <= nwindows; n++) {
    if (n == wn)
      continue;
    for (cc = 1; cc <= window[n]->enon; cc++) {
      c = window[n]->wsa[cc];
      if (geom->fm[c].wn == wn) {
        cl = geom->fm[c].sc;
        if (cl <= window[wn]->enonS) {
          window[wn]->s2m[window[wn]->ns2mS] = cl;
          window[wn]->ns2mS += 1;
        } else {
          window[wn]->s2m[window[wn]->ns2m] = cl;
          window[wn]->ns2m += 1;
        }
      }
    }
  }
  if (use_ghosts) {
    for (cc = 1; cc <= geom->nbpt; cc++) {
      c = geom->bin[cc];
      if (geom->zoomc[geom->m2d[c]] & (ZC | ZE1 | ZE2) &&
          geom->fm[c].wn == wn) {
        /* if(geom->fm[c].wn == wn) { */
        cl = geom->fm[c].sc;
        if (cl <= window[wn]->enonS) {
          window[wn]->s2m[window[wn]->ns2mS] = cl;
          window[wn]->ns2mS += 1;
        } else {
          window[wn]->s2m[window[wn]->ns2m] = cl;
          window[wn]->ns2m += 1;
        }
      }
    }
  }
  window[wn]->ns2m -= 1;
  window[wn]->ns2mS -= 1;

  /* Get the e1 and e2 faces for slave to master transfers */
  if (window[wn]->zoomf >= 1) {
    window[wn]->s2me1 = i_alloc_1d(window[wn]->ns2m + 1);
    window[wn]->s2me2 = i_alloc_1d(window[wn]->ns2m + 1);
    zse1 = (int)(window[wn]->zmfe1 / 2);
    zse2 = (int)(window[wn]->zmfe2 / 2);
    for (cc = 1; cc <= window[wn]->ns2m; cc++) {
      c = window[wn]->s2m[cc];
      ce1 = ce2 = window[wn]->wsa[c];
      if (window[wn]->zoom == DOZOOM) {
	for (zc = 1; zc <= zse1; zc++)
	  ce1 = geom->xm1[ce1];
	for (zc = 1; zc <= zse2; zc++)
	  ce2 = geom->ym1[ce2];
      }
      window[wn]->s2me1[cc] = ce1;
      window[wn]->s2me2[cc] = ce2;
    }
  } else {
    window[wn]->s2me1 = window[wn]->s2m;
    window[wn]->s2me2 = window[wn]->s2m;
  }

  if (DEBUG("init_w")) {
    dlog("init_w",
         "  Window %d : %d (%d 2D) master to slave transfer cells\n", wn,
         window[wn]->nm2s, window[wn]->nm2sS);
    dlog("init_w",
         "  Window %d : %d (%d 2D) slave to master transfer cells\n", wn,
         window[wn]->ns2m, window[wn]->ns2mS);
  }

  i_free_1d(wm);
}

/* END build_transfer_maps()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to get the mapping from the windows cell centered         */
/* locations to the corresponding zoomed location on the master for  */
/* e1 and e2 cells. Note that R_EDGE and F_EDGE boundary cells 'hang */
/* off' the grid and are not shifted on the zoomed master grid.      */
/*-------------------------------------------------------------------*/
void build_zoom_maps(geometry_t *geom,       /* Global geometery     */
		     geometry_t **window     /* Local geometry       */
  )
{
  int cc, c, c2, n, wn, i;

  if (geom->zoom != NOZOOM) {
    /*---------------------------------------------------------------*/
    /* e1 maps                                                       */
    geom->zse1 = i_alloc_1d(geom->sgsiz);
    for (cc = 1; cc <= geom->b3_e1; cc++) {
      c = c2 = geom->w3_e1[cc];
      wn = geom->fm[c].wn;
      for (i = 1; i <= (int)(window[wn]->zmfe1 / 2); i++)
	c2 = geom->xm1[c2];
      geom->zse1[c] = c2;
    }
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      for (cc = 1; cc <= open->no3_e1; cc++) {
	c = open->obc_e1[cc];
	if (open->ocodex & R_EDGE)
	  geom->zse1[c] = c;
      }
    }

    /*---------------------------------------------------------------*/
    /* e2 maps                                                       */
    geom->zse2 = i_alloc_1d(geom->sgsiz);
    for (cc = 1; cc <= geom->b3_e2; cc++) {
      c = c2 = geom->w3_e2[cc];
      wn = geom->fm[c].wn;
      for (i = 1; i <= (int)(window[wn]->zmfe2 / 2); i++)
	c2 = geom->ym1[c2];
      geom->zse2[c] = c2;
    }
    for (n = 0; n < geom->nobc; n++) {
      open_bdrys_t *open = geom->open[n];
      for (cc = 1; cc <= open->no3_e2; cc++) {
	c = open->obc_e2[cc];
	if (open->ocodey & F_EDGE)
	  geom->zse2[c] = c;
      }
    }
  }
}

/* END build_zoom_maps()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to populate the master with window data                   */
/*-------------------------------------------------------------------*/
void master_fill(master_t *master,  /* Master data structure */
                 geometry_t **window, /* Window geometry */
                 window_t **windat, /* Window data */
                 win_priv_t **wincon  /* Private window data */
  )
{
  int n;                        /* Window counter */

  if (master->nwindows == 1)
    return;
  if (master->is_filled)
    return;
  for (n = 1; n <= master->nwindows; n++) {
    s2m_2d(master, window[n], windat[n]);
    s2m_3d(master, window[n], windat[n], wincon[n]);
  }
  master->is_filled = 1;
}

/* END master_fill()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to populate the master with window data at time series    */
/* locations. Note: The whole column is transferred since the exact  */
/* location in the water column of the time series point may vary    */
/* according to the reference point.                                 */
/*-------------------------------------------------------------------*/
void master_fill_ts(master_t *master,     /* Master data             */
		    geometry_t **window,  /* Window geometry         */
		    window_t **windat,    /* Window data             */
		    win_priv_t **wincon   /* Window constants        */
  )
{
  /* Counters */
  int n, cc, lc, c, tn, tt; 

  if (master->nwindows == 1)
    return;
  if (master->is_filled)
    return;

  for (n = 1; n <= master->nwindows; n++) {
    if (window[n]->ns2m_ts == 0) continue;
    for (cc = 1; cc <= window[n]->ns2m_ts; cc++) {
      lc = window[n]->s2m_ts[cc];
      c = window[n]->wsa[lc];
      if (lc <= window[n]->ewetS) {
	master->eta[c] = windat[n]->eta[lc];
	master->u1av[c] = windat[n]->u1av[lc];
	master->u2av[c] = windat[n]->u2av[lc];
	for (tt = 0; tt < master->ntrS; tt++) {
	  master->tr_wcS[tt][c] = windat[n]->tr_wcS[tt][lc];
	}
      }
      master->u1[c] = windat[n]->u1[lc];
      master->u2[c] = windat[n]->u2[lc];
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	master->tr_wc[tn][c] = windat[n]->tr_wc[tn][lc];
      }
    }
  }
}

/* END master_fill_ts()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Transfers columns corresponding to a stencil surrounding the a    */
/* glider location at time t for variables to be compared to glider  */
/* observations.                                                     */
/*-------------------------------------------------------------------*/
void master_fill_glider(master_t *master,     /* Master data         */
			geometry_t **window,  /* Window geometry     */
			window_t **windat,    /* Window data         */
			win_priv_t **wincon,  /* Window constants    */
			ts_point_t *ts,
			double t
			)
{
  geometry_t *geom = master->geom;
  int wn, tn, cg, c, cs, lc;
  int m, i, found;
  int *st = NULL, ssize;    

  if (!(ts->metric & TS_GLIDER)) return;

  /* Get the cell the glider resides in                              */
  cg = get_glider_loc(master, ts, &ts->ts, t);

  /* Get a neighbourhood around the glider cell                      */
  ssize = ts->kernal;
  st = stencil(geom, cg, &ssize, 0);

  /* Transfer data to the master within the neighbourhood            */
  for (m = 0; m < ssize; m++) {
    c = st[m];                    /* Stencil coordinate              */
    cs = geom->m2d[c];            /* Surface coordinate              */
    for (wn = 1; wn <= master->nwindows; wn++) {
      if (geom->fm[c].wn == wn) { /* Check if c lies in window wn    */
	lc = geom->fm[c].sc;      /* Local coordinate                */

	c = cs;
	/* Loop down the water column                                */
	while (c != geom->zm1[c]) {
	  /* Loop over variables to be compared to glider obs        */
	  for (tn = 0; tn < ts->dnvars; tn++) {
	    if (strcmp(ts->dvars[tn], "N2") == 0) {
	      ts->data[tn][c] = windat[wn]->dens[lc];
	    } else {
	      found = 0;
	      for (i = 0; i < master->ntr; i++) {
		if (strcmp(ts->dvars[tn], master->trinfo_3d[i].name) == 0) {
		  ts->data[tn][c] = windat[wn]->tr_wc[i][c];
		  found = 1;
		}
	      }
	      if (found == 0) {
		for (i = 0; i < master->ntrS; i++) {
		  if (strcmp(ts->dvars[tn], master->trinfo_2d[i].name) == 0) {
		    ts->data[tn][geom->m2d[c]] = windat[wn]->tr_wcS[i][window[wn]->m2d[lc]];
		    found = 1;
		  }
		}
	      }
	    }
	  }
	  c = geom->zm1[c];
	}
      }
    }
  }
}

/* END master_fill_glider()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d local arrays into the master                             */
/*-------------------------------------------------------------------*/
void s2m_3d(master_t *master,   /* Model grid data structure */
            geometry_t *window, /* Window to pack */
            window_t *windat,   /* Window data structure */
            win_priv_t *wincon  /* Private window data */
  )
{
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */
  int lc;                       /* Local sparse coordinate */
  int n, tn, tt;                /* Counters */
  int zse1 = (int)(window->zmfe1 / 2);
  int zse2 = (int)(window->zmfe2 / 2);

  /* Grid centered variables */
  for (cc = 1; cc <= window->b3_t; cc++) {
    lc = window->w3_t[cc];
    c = window->wsa[lc];
    
    master->w[c] = windat->w[lc];
    master->dens[c] = windat->dens[lc];
    master->dens_0[c] = windat->dens_0[lc];
    master->Vz[c] = windat->Vz[lc];
    master->Kz[c] = windat->Kz[lc];
    master->dz[c] = wincon->dz[lc];
    master->u1vh[c] = wincon->u1vh[lc];
    master->u2vh[c] = wincon->u2vh[lc];
    for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
      tn = master->trmap_s2m_3d[tt];
      master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
    }
  }
  
  /* e1 face centered variables */
  for (cc = 1; cc <= window->b3_e1; cc++) {
    lc = window->w3_e1[cc];
    c = window->wsa[lc];
    if (zse1) c = geom->zse1[c];
    master->u1[c] = windat->u1[lc];
    master->dzu1[c] = windat->dzu1[lc];
  }
  
  /* e2 face centered variables */
  for (cc = 1; cc <= window->b3_e2; cc++) {
    lc = window->w3_e2[cc];
    c = window->wsa[lc];
    if (zse2) c = geom->zse2[c];
    master->u2[c] = windat->u2[lc];
    master->dzu2[c] = windat->dzu2[lc];
  }

  /* Fill the R_EDGE and F_EDGE cell centres if required */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    if (open->ocodex & R_EDGE) {
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	if (master->trinfo_3d[tn].type & E1VAR) {
	  for (cc = 1; cc <= open->no3_e1; cc++) {
	    lc = open->obc_e1[cc];
	    c = window->wsa[lc];
	    master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
	  }
	}
      }
    }
    if (open->ocodey & F_EDGE) {
      for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	if (master->trinfo_3d[tn].type & E2VAR) {
	  for (cc = 1; cc <= open->no3_e2; cc++) {
	    lc = open->obc_e2[cc];
	    c = window->wsa[lc];
	    master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
	  }
	}
      }
    }
  }
}

/* END s2m_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d local arrays into the master                             */
/*-------------------------------------------------------------------*/
void s2m_2d(master_t *master,   /* Model grid data structure */
            geometry_t *window, /* Window to pack */
            window_t *windat    /* Window data structure */
  )
{
  int cc;                       /* Local sparse coordinate counter */
  int c, cb;                    /* Global sparse coordinate */
  int lc;                       /* Local sparse coordinate */
  int tn, k = 0; 
  int zse1 = (int)(window->zmfe1 / 2);
  int zse2 = (int)(window->zmfe2 / 2);

  /* Grid centered variables. Note: wind1, wind2, patm and Cd are */
  /* already on the master since these variables are read onto the */
  /* master and transferred to the window at the start of the step.  */
  for (cc = 1; cc <= window->b2_t; cc++) {
    lc = window->w2_t[cc];
    c = window->wsa[lc];
    master->eta[c] = windat->eta[lc];
    master->topz[c] = windat->topz[lc];
    master->wtop[c] = windat->wtop[lc];
    for (tn = 0; tn < master->ntrS; tn++)
      master->tr_wcS[tn][c] = windat->tr_wcS[tn][lc];
    for (tn = 0; tn < windat->nsed; tn++) {
      for (k = 0; k < window->sednz; k++)
        master->tr_sed[tn][k][c] = windat->tr_sed[tn][k][lc];
    }
  }

  /* e1 face centered variables */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    lc = window->w2_e1[cc];
    cb = window->wsa[window->bot_e1[cc]];
    c = window->wsa[lc];
    if (zse1) c = geom->zse1[c];
    master->u1av[c] = windat->u1av[lc];
    /* Re-setting the master for printing purposes can interfere     */
    /* with master-slave transfers.                                  */
    /* master->u1bot[c]=windat->u1[cb]; */
    /* master->u1bot[c] = windat->u1bot[c]; */
  }
  
  /* e2 face centered variables */
  for (cc = 1; cc <= window->b2_e2; cc++) {
    lc = window->w2_e2[cc];
    cb = window->wsa[window->bot_e1[cc]];
    c = window->wsa[lc];
    if (zse2) c = geom->zse2[c];
    master->u2av[c] = windat->u2av[lc];
    /* master->u2bot[c]=windat->u2[cb]; */
    /* master->u2bot[c] = windat->u2bot[c]; */
  }
}

/* END s2m_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d local arrays into the master for refined grids. This is  */
/* done over all wet cells to process, as opped to transfer cells.   */
/*-------------------------------------------------------------------*/
void s2m_zoom(master_t *master,   /* Model grid data structure */
	      geometry_t *window, /* Window to pack */
	      window_t *windat,   /* Window data structure */
	      int mode            /* Transfer mode */
  )
{
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */
  int lc;                       /* Local sparse coordinate */
  int zse1 = window->zmee1;
  int zse2 = window->zmee2;

  if (window->zoom & PRECOND)
    zse1 = zse2 = 0;

  if (mode & NVELOCITY) {
    /* Transfer updated 2D velocity                                */
    for (cc = 1; cc <= window->b2_e1; cc++) {
      lc = window->w2_e1[cc];
      c = window->wsa[lc];
      if (zse1) c = geom->zse1[c];
      master->nu1av[c] = windat->nu1av[lc];
    }
    for (cc = 1; cc <= window->b2_e2; cc++) {
      lc = window->w2_e2[cc];
      c = window->wsa[lc];
      if (zse2) c = geom->zse2[c];
      master->nu2av[c] = windat->nu2av[lc];
    }
  } else if (mode & VELOCITY) {
    /* Transfer 2D velocity                                          */
    for (cc = 1; cc <= window->b2_e1; cc++) {
      lc = window->w2_e1[cc];
      c = window->wsa[lc];
      if (zse1) c = geom->zse1[c];
      master->u1av[c] = windat->u1av[lc];
    }
    for (cc = 1; cc <= window->b2_e2; cc++) {
      lc = window->w2_e2[cc];
      c = window->wsa[lc];
      if (zse2) c = geom->zse2[c];
      master->u2av[c] = windat->u2av[lc];
    }
  } else if (mode & ETA_A) {
    /* Transfer surface elevation                                    */
    for (cc = 1; cc <= window->b2_t; cc++) {
      lc = window->w2_t[cc];
      c = window->wsa[lc];
      master->eta[c] = windat->eta[lc];
    }
  }
}

/* END s2m_zoom()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the local fluxes to the master, scaling the   */
/* flux in the process if a zoomed flux is transferred.              */
/*-------------------------------------------------------------------*/
void s2m_flux(geometry_t *geom, /* Sparse global geometery */
              geometry_t *window, /* Processing window */
              double *Ag,       /* Global flux array */
              double *Al,       /* Local flux array */
              int *vec,         /* Cells on which Al is defined */
              int *evec,        /* Cells on which Ag is defined */
              int nvec,         /* Size of vec */
              int mode          /* 1 for e1 fluxes, 0 for e2 velocity */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */
  double *dhl, *dhg;            /* Grid spacing */

  if (window->zoom != NOZOOM) {
    /* Get the correct local and global horizontal length scale */
    if (mode) {
      dhl = window->h2au1;
      dhg = geom->h2au1;
    } else {
      dhl = window->h1au2;
      dhg = geom->h1au2;
    }
    /* Transfer the flux to the master with scaling */
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Ag[c] = Al[lc] * (dhg[geom->m2d[c]] / dhl[window->m2d[lc]]);
    }
  } else {
    /* Transfer the flux to the master */
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = evec[cc];
      Ag[c] = Al[lc];
    }
  }
}

/* END s2m_flux()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the local velocities to the master.           */
/*-------------------------------------------------------------------*/
void s2m_vel(double *Ag,        /* Global flux array */
             double *Al,        /* Local flux array */
             int *vec,          /* Cells on which A is defined */
             int *evec,         /* Cells on which Ag is defined */
             int nvec           /* Size of vec */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter */

  /* Transfer the velocity to the master */
  for (cc = 1; cc <= nvec; cc++) {
    lc = vec[cc];
    c = evec[cc];
    Ag[c] = Al[lc];
  }
}

/* END s2m_vel()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void s2c_2d(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double **ac,        /* Output cartesian array */
            int nx,             /* x dimension of 2d array */
            int ny              /* y dimension of 2d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int c;                        /* Global sparse coordinate */

  /* Initialise */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      ac[j][i] = NaN;

  /* Map the 2d sparse variable to Cartesion */
  /* Grid centered variables : vertical velocity */
  /*
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
  */
  for (c = 1; c <= geom->ewetS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    if (i == NOTVALID && j == NOTVALID)
      continue;
    if (i < nx && j < ny)
      ac[j][i] = as[c];
  }
}

/* END s2c_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void s2c_3d(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double ***ac,       /* Output cartesian array */
            int nx,             /* x dimension of 3d array */
            int ny,             /* y dimension of 3d array */
            int nz              /* z dimension of 3d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int k;                        /* Counter in the z direction */
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */

  if (dumpdata->vmap != NULL) {
    /* Linearly interpolate onto the geom->layers[] distribution */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++) {
	  ac[k][j][i] = NaN;
	  c = dumpdata->vmap[k][j][i];
	  if (c > 0) {
	    int cs = geom->m2d[c];
	    int zm1 = geom->zm1[c];
	    double d1 = geom->cellz[c] * dumpdata->botz[j][i];
	    double d2 = geom->cellz[zm1] * dumpdata->botz[j][i];
	    ac[k][j][i] = (fabs(geom->layers[k]) - d1) * 
	      (as[c] - as[zm1]) / (d1 - d2) + as[c];
	  }
	}
  } else {
    /* Initialise */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  ac[k][j][i] = 0.0;

    /* Map the 3d sparse variable to Cartesion */
    /*
      for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
    */
    for (cc = 1; cc <= geom->a3_t; cc++) {
      c = geom->wsa[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
	continue;
      if (i < nx && j < ny)
	ac[k][j][i] = as[c];
    }
  }
}

void s2c_3d_e1(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double ***ac,       /* Output cartesian array */
            int nx,             /* x dimension of 3d array */
            int ny,             /* y dimension of 3d array */
            int nz              /* z dimension of 3d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int k;                        /* Counter in the z direction */
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */

  /* Initialise */
  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
	ac[k][j][i] = 0.0;
  
  /* Map the 3d sparse variable to Cartesion */
  for (cc = 1; cc <= geom->b3_e1; cc++) {
    c = geom->w3_e1[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;
    if (i < nx && j < ny)
      ac[k][j][i] = as[c];
  }
}

void s2c_3d_e2(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double ***ac,       /* Output cartesian array */
            int nx,             /* x dimension of 3d array */
            int ny,             /* y dimension of 3d array */
            int nz              /* z dimension of 3d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int k;                        /* Counter in the z direction */
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */

  /* Initialise */
  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
	ac[k][j][i] = 0.0;
  
  /* Map the 3d sparse variable to Cartesion */
  for (cc = 1; cc <= geom->b3_e2; cc++) {
    c = geom->w3_e2[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;
    if (i < nx && j < ny)
      ac[k][j][i] = as[c];
  }
}

/* END s2c_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d sparse sediment arrays into the Cartesian array          */
/*-------------------------------------------------------------------*/
void s2c_sed(geometry_t *geom,  /* Global geometry structure */
            double **as,        /* Input sparse array */
            double ***ac,       /* Output cartesian array */
            int nx,             /* x dimension of 3d array */
            int ny,             /* y dimension of 3d array */
            int nz              /* z dimension of 3d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int k;                        /* Counter in the z direction */
  int cc;                       /* Local sparse coordinate counter */
  int c;                        /* Global sparse coordinate */

  /* Initialise */
  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        ac[k][j][i] = NaN;

  /* Map the 2d sparse variable to Cartesion */
  /* Grid centered variables : vertical velocity */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    i = geom->s2i[c];
    j = geom->s2j[c];
    for (k = 0; k < nz; k++) {
      if (i < nx && j < ny)
	ac[k][j][i] = as[nz-k-1][c];
    }
  }
}

/* END s2c_sed()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d Cartesian arrays into the sparse array                   */
/*-------------------------------------------------------------------*/
void c2s_2d(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double **ac,        /* Output cartesian array */
            int nx,             /* x dimension of 2d array */
            int ny              /* y dimension of 2d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j, k;                     /* Counter in the y direction */
  int c;                        /* Global sparse coordinate */

  for (c = 1; c <= geom->sgnumS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;
    if (i < nx && j < ny)
      as[c] = ac[j][i];
  }
  /* 
     for(i=0; i<nx; i++) for(j=0; j<ny; j++) {
     c=geom->map[geom->nz-1][j][i]; if(c>0 &&
     c<geom->sgnumS)as[c]=ac[j][i]; } */
}

/* END c2s_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void c2s_3d(geometry_t *geom,   /* Global geometry structure */
            double *as,         /* Input sparse array */
            double ***ac,       /* Output cartesian array */
            int nx,             /* x dimension of 3d array */
            int ny,             /* y dimension of 3d array */
            int nz              /* z dimension of 3d array */
  )
{
  int i;                        /* Counter in the x direction */
  int j;                        /* Counter in the y direction */
  int k;                        /* Counter in the z direction */
  int c;                        /* Global sparse coordinate */

  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID)
      continue;
    if (i < nx && j < ny)
      as[c] = ac[k][j][i];
  }
  /* 
     for(k=0; k<nz; k++) for(i=0; i<nx; i++) for(j=0; j<ny; j++) {
     c=geom->map[k][j][i]; if(c>0 && c<geom->sgnum)as[c]=ac[k][j][i]; } */
}

/* END c2s_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs the sparse array into consecutive wet cells. If the 'map'   */
/* is geom->wsa then the entire wet grid, boundaries included, is    */
/* packed.                                                           */
/* Note the the first null cell in the sparse array is removed so    */
/* that the dumped variable indies go from 0:ns-1 whereas the sparse */
/* indicies go from 1:ns with value[0]=0.                            */
/*-------------------------------------------------------------------*/
void pack_sparse(int *map, int mapsize, double *var, double *pack)
{
  int cc, c;
  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    pack[cc - 1] = var[c];
  }
}

/* END pack_sparse()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Unpacks an array of consecutive wet cells into the sparse array.  */
/* Note the the first null cell in the sparse array is padded so     */
/* that the sparse variable indies go from 1:ns with value[0]=0,     */
/* while the file variable indicies go from 0:ns-1.                  */
/*-------------------------------------------------------------------*/
void unpack_sparse(int *map, int mapsize, double *var, double *unpack)
{
  int cc, c;
  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    unpack[c] = var[cc - 1];
  }
}

/* END unpack_sparse()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Fills a relaxation structure.                                     */
/*-------------------------------------------------------------------*/
void relax_fill(geometry_t *window, relax_info_t *in, relax_info_t *out)
{
  int c, cc;

  out->rate = in->rate;
  out->dv0 = in->dv0;
  out->tc0 = in->tc0;
  out->slope = in->slope;
  out->tctype = in->tctype;
  if (out->val1) {
    for (cc = 1; cc <= out->size; cc++) {
      c = window->wsa[cc];
      out->val1[cc] = in->val1[c];
    }
  }
  if (out->val2) {
    for (cc = 1; cc <= out->size; cc++) {
      c = window->wsa[cc];
      out->val2[cc] = in->val2[c];
    }
  }
}

/* END relax_fill()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to check window transfers                                 */
/* Don't free((global_map_t *)geom->fm) in geom_free(), window.c if  */
/* this routine is used.                                             */
/*-------------------------------------------------------------------*/
void check_transfers(geometry_t *geom,    /* Global geometry         */
		     geometry_t **window, /* Window geometry         */
		     window_t **windat,   /* Window data             */
		     win_priv_t **wincon, /* Window constants        */
		     int nw,              /* Number of windows       */
		     int mode
		     )
{
  int cc, c, cg, cl, zp1;
  int n, m, tr;
  int i, j, k;
  int cellt;
  char type[MAXSTRLEN], buf[MAXSTRLEN], tp;

  if (geom->fm == NULL)
    hd_quit("Must not free((global_map_t *)geom->fm) in geom_free(), window.c\n");

  if (nw == 1)
    return;

  for (n = 1; n <= nw; n++) {
    /*---------------------------------------------------------------*/
    /* Check the fluxes after eta_step(). The fluxes used for        */
    /* elevation calculation are set in auxiliary cells and stored   */
    /* in the buffers wincon->d2, d3. The fluxes windat->u1flux,     */
    /* u2flux are used to adjust 3D velocities and are not required  */
    /* at auxiliary cells. Only fluxes in auxiliary cells at forward */
    /* and right faces are set.                                      */
    if (mode == (VEL2D|FLUX)) {
      for (cc = 1; cc <= window[n]->a2_t; cc++) {
	c = window[n]->w2_t[cc];
	cg = window[n]->wsa[c];
	i = geom->s2i[cg];
	j = geom->s2j[cg];
	k = geom->s2k[cg];
	strcpy(type,"AUX");
	if (geom->fm[cg].wn != n) {
	  m = geom->fm[cg].wn;
	  cl = geom->fm[cg].sc;
	  if (geom->fm[geom->xm1[cg]].wn == n)
	    tr_warn("u1flux",wincon[m]->d2[cl],wincon[n]->d2[c],
		    c,i,j,-1,cg,n,type,' ');
	  if (geom->fm[geom->ym1[cg]].wn == n) {
	    tr_warn("u2flux",wincon[m]->d3[cl],wincon[n]->d3[c],
		    c,i,j,-1,cg,n,type,' ');
	  }
	}
      }
    } else if (mode == (VEL3D|FLUX)) {
      for (cc = 1; cc <= window[n]->a3_t; cc++) {
	c = window[n]->w3_t[cc];
	cg = window[n]->wsa[c];
	i = geom->s2i[cg];
	j = geom->s2j[cg];
	k = geom->s2k[cg];
	strcpy(type,"AUX");
	if (geom->fm[cg].wn != n) {
	  m = geom->fm[cg].wn;
	  cl = geom->fm[cg].sc;
	  if (geom->fm[geom->xm1[cg]].wn == n)
	    tr_warn("u1flux3d",windat[m]->u1flux3d[cl],windat[n]->u1flux3d[c],
		    c,i,j,k,cg,n,type,' ');
	  if (geom->fm[geom->ym1[cg]].wn == n) {
	    tr_warn("u2flux3d",windat[m]->u2flux3d[cl],windat[n]->u2flux3d[c],
		    c,i,j,k,cg,n,type,' ');
	  }
	}
      }
    } else {
      /*---------------------------------------------------------------*/
      /* Check all variables after the window is filled                */
      for (c = 1; c <= window[n]->enon; c++) {
	cg = window[n]->wsa[c];
	i = geom->s2i[cg];
	j = geom->s2j[cg];
	k = geom->s2k[cg];
	zp1 = window[n]->wsa[window[n]->zp1[c]];
	cellt = 0;
	strcpy(type,"WET");
	if (geom->fm[cg].wn != n) {
	  cellt = 1;
	  strcpy(type,"AUX");
	}
	if (c == window[n]->zm1[c]) {
	  cellt = 2;
	  strcpy(type,"SED");
	}
	if (i == NOTVALID && j == NOTVALID && k == NOTVALID) {
	  cellt = 4;
	  strcpy(type,"GHOST");
	}	

	if (cellt == 1) {
	  m = geom->fm[cg].wn;
	  cl = geom->fm[cg].sc;
	  /*-----------------------------------------------------------*/
	  /* Check 3D variables                                        */
	  if (mode & VEL3D) {
	    /* Cell centered                                           */
	    cc = ANY(c, window[n]->w3_t, window[n]->n3_t);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v3_t)
	      tp = 'w';
	    else if (cc > window[n]->v3_t && c <= window[n]->b3_t)
	      tp = 'b';
	    else if (cc > window[n]->b3_t && c <= window[n]->a3_t)
	      tp = 'a';
	    else if (cc > window[n]->a3_t && c <= window[n]->n3_t)
	      tp = 'g';
	    if (mode & MIXING) {
	      tr_warn("Kz",windat[m]->Kz[cl],windat[n]->Kz[c],
		      c,i,j,k,cg,n,type,tp);
	      tr_warn("Vz",windat[m]->Vz[cl],windat[n]->Vz[c],
		      c,i,j,k,cg,n,type,tp);
	    } else {
	      tr_warn("w",windat[m]->w[cl],windat[n]->w[c],
		      c,i,j,k,cg,n,type,tp);
	      tr_warn("dens",windat[m]->dens[cl],windat[n]->dens[c],
		      c,i,j,k,cg,n,type,tp);
	      for (tr = 0; tr < windat[n]->ntr; tr++) {
		sprintf(buf, "3D tracer%d", tr);
		tr_warn(buf,windat[m]->tr_wc[tr][cl],windat[n]->tr_wc[tr][c],
			c,i,j,k,cg,n,type,tp);
	      }
	      /*
	      for (tr = 0; tr < windat[n]->ntrS; tr++) {
		sprintf(buf, "2D tracer%d", tr);
		tr_warn(buf,windat[m]->tr_wcS[tr][cl],windat[n]->tr_wcS[tr][c],
			c,i,j,k,cg,n,type,tp);
	      }
	      */
	      if (c <= window[n]->enonS) {
		tr_warn("eta",windat[m]->eta[cl],windat[n]->eta[c],
			c,i,j,k,cg,n,type,tp);
		tr_warn("wtop",windat[m]->wtop[cl],windat[n]->wtop[c],
			c,i,j,k,cg,n,type,tp);
		tr_warn("wbot",windat[m]->wbot[cl],windat[n]->wbot[c],
			c,i,j,k,cg,n,type,tp);
		tr_warn("cellarea",window[m]->cellarea[cl],
			window[n]->cellarea[c],c,i,j,k,cg,n,type,tp);
	      }
	    }

	    /* e1 face centered                                        */
	    cc = ANY(c, window[n]->w3_e1, window[n]->n3_e1);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v3_e1)
	      tp = 'w';
	    else if (cc > window[n]->v3_e1 && c <= window[n]->b3_e1)
	      tp = 'b';
	    else if (cc > window[n]->b3_e1 && c <= window[n]->a3_e1)
	      tp = 'a';
	    else if (cc > window[n]->a3_e1 && c <= window[n]->a3_e1)
	      tp = 'x';
	    else if (cc > window[n]->a3_e1 && c <= window[n]->n3_e1)
	      tp = 'g';
	    if (mode & MIXING) {
	      if(!wincon[m]->thin_merge || (wincon[m]->thin_merge && 
		  windat[m]->dzu1[window[m]->zp1[cl]] >= wincon[m]->hmin))
		 tr_warn("dzu1",windat[m]->dzu1[cl],windat[n]->dzu1[c],
			 c,i,j,k,cg,n,type,tp);
	    } else {
	      tr_warn("u1",windat[m]->u1[cl],windat[n]->u1[c],
		      c,i,j,k,cg,n,type,tp);
	      tr_warn("u1b",windat[m]->u1b[cl],windat[n]->u1b[c],
		      c,i,j,k,cg,n,type,tp);
	      if (c <= window[n]->enonS) {
		tr_warn("u1bot",windat[m]->u1bot[cl],windat[n]->u1bot[c],
			c,i,j,k,cg,n,type,tp);
	      }
	    }
	    
	    /* e2 face centered                                        */
	    cc = ANY(c, window[n]->w3_e2, window[n]->n3_e2);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v3_e2)
	      tp = 'w';
	    else if (cc > window[n]->v3_e2 && c <= window[n]->b3_e2)
	      tp = 'b';
	    else if (cc > window[n]->b3_e2 && c <= window[n]->a3_e2)
	      tp = 'a';
	    else if (cc > window[n]->a3_e2 && c <= window[n]->a3_e2)
	      tp = 'x';
	    else if (cc > window[n]->a3_e2 && c <= window[n]->n3_e2)
	      tp = 'g';
	    if (mode & MIXING) {
	      if(!wincon[m]->thin_merge || (wincon[m]->thin_merge && 
		  windat[m]->dzu2[window[m]->zp1[cl]] >= wincon[m]->hmin))
		 tr_warn("dzu2",windat[m]->dzu2[cl],windat[n]->dzu2[c],
			 c,i,j,k,cg,n,type,tp);
	    } else {
	      tr_warn("u2",windat[m]->u2[cl],windat[n]->u2[c],
		      c,i,j,k,cg,n,type,tp);
	      tr_warn("u2b",windat[m]->u2b[cl],windat[n]->u2b[c],
		      c,i,j,k,cg,n,type,tp);
	      if (c <= window[n]->enonS) {
		tr_warn("u2bot",windat[m]->u2bot[cl],windat[n]->u2bot[c],
			c,i,j,k,cg,n,type,tp);
	      }
	    }
	  }

	  /*-----------------------------------------------------------*/
	  /* Check 2D variables                                        */
	  if (mode & VEL2D && c <= window[n]->enonS) {
	    /* Cell centered                                           */
	    cc = ANY(c, window[n]->w2_t, window[n]->n2_t);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v2_t)
	      tp = 'w';
	    else if (cc > window[n]->v2_t && c <= window[n]->b2_t)
	      tp = 'b';
	    else if (cc > window[n]->b2_t && c <= window[n]->a2_t)
	      tp = 'a';
	    else if (cc > window[n]->a2_t && c <= window[n]->n2_t)
	      tp = 'g';
	    tr_warn("botz",window[m]->botz[cl],window[n]->botz[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h1acell",window[m]->h1acell[cl],window[n]->h1acell[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h2acell",window[m]->h2acell[cl],window[n]->h2acell[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("Cd",wincon[m]->Cd[cl],wincon[n]->Cd[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("eta",windat[m]->eta[cl],windat[n]->eta[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("etab",windat[m]->etab[cl],windat[n]->etab[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("detadt",windat[m]->detadt[cl],windat[n]->detadt[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("patm",windat[m]->patm[cl],windat[n]->patm[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("topz",windat[m]->topz[cl],windat[n]->topz[c],
		    c,i,j,-1,cg,n,type,tp);
	    
	    /* e1 face centered                                        */
	    cc = ANY(c, window[n]->w2_e1, window[n]->n2_e1);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v2_e1)
	      tp = 'w';
	    else if (cc > window[n]->v2_e1 && c <= window[n]->b2_e1)
	      tp = 'b';
	    else if (cc > window[n]->b2_e1 && c <= window[n]->a2_e1)
	      tp = 'a';
	    else if (cc > window[n]->a2_e1 && c <= window[n]->a2_e1)
	      tp = 'x';
	    else if (cc > window[n]->a2_e1 && c <= window[n]->n2_e1)
	      tp = 'g';
	    tr_warn("botzu1",window[m]->botzu1[cl],window[n]->botzu1[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h1au1",window[m]->h1au1[cl],window[n]->h1au1[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h1au2",window[m]->h1au2[cl],window[n]->h1au2[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u1av",windat[m]->u1av[cl],windat[n]->u1av[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u1avb",windat[m]->u1avb[cl],windat[n]->u1avb[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u1bot",windat[m]->u1bot[cl],windat[n]->u1bot[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("depth_e1",windat[m]->depth_e1[cl],windat[n]->depth_e1[c],
		    c,i,j,-1,cg,n,type,tp);
	    
	    /* e2 face centered                                        */
	    cc = ANY(c, window[n]->w2_e2, window[n]->n2_e2);
	    tp = ' ';
	    if (cc > 0 && cc <= window[n]->v2_e2)
	      tp = 'w';
	    else if (cc > window[n]->v2_e2 && c <= window[n]->b2_e2)
	      tp = 'b';
	    else if (cc > window[n]->b2_e2 && c <= window[n]->a2_e2)
	      tp = 'a';
	    else if (cc > window[n]->a2_e2 && c <= window[n]->a2_e2)
	      tp = 'x';
	    else if (cc > window[n]->a2_e2 && c <= window[n]->n2_e2)
	      tp = 'g';
	    tr_warn("botzu2",window[m]->botzu2[cl],window[n]->botzu2[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h2au2",window[m]->h2au1[cl],window[n]->h2au1[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("h2au1",window[m]->h2au1[cl],window[n]->h2au1[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u2av",windat[m]->u2av[cl],windat[n]->u2av[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u2avb",windat[m]->u2avb[cl],windat[n]->u2avb[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("u2bot",windat[m]->u2bot[cl],windat[n]->u2bot[c],
		    c,i,j,-1,cg,n,type,tp);
	    tr_warn("depth_e2",windat[m]->depth_e2[cl],windat[n]->depth_e2[c],
		    c,i,j,-1,cg,n,type,tp);
	  }
	}
      }
    }
  }
}

/* END check_win_2d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Resets all the arrays on the windows                              */
/*-------------------------------------------------------------------*/
void window_reset(master_t *master,      /* Master data              */
		  geometry_t *window,    /* Window geometry          */
		  window_t *windat,      /* Window data              */
		  win_priv_t *wincon,    /* Window constants         */
		  int mode
		  )
{
  int c, cc, cs, n, k, tn;

  if (mode & RS_VH) {
    for (cc = 1; cc <= window->enon; cc++) {
      c = window->wsa[cc];
      wincon->u1vh[cc] = master->u1vh[c];
      wincon->u2vh[cc] = master->u2vh[c];
    }
    master->dtb = master->dtf = master->dt = master->grid_dt;
    master->dt2d = master->grid_dt / master->iratio;
    windat->dtf2 = master->dt / master->iratio;
    return;
  }

  /* Initialise all 3D and 2D variables required by the window       */
  master->dtb = master->dtf = master->dt = master->grid_dt;
  master->dt2d = master->grid_dt / master->iratio;
  windat->dtf2 = master->dt / master->iratio;
  for (cc = 1; cc <= window->enon; cc++) {
    c = window->wsa[cc];
    for (tn = 0; tn < windat->ntr; tn++) {
      windat->tr_wc[tn][cc] = master->tr_wc[tn][c];
    }

    windat->w[cc] = master->w[c];
    windat->Kz[cc] = master->Kz[c];
    windat->Vz[cc] = master->Vz[c];
    wincon->u1vh[cc] = master->u1vh[c];
    wincon->u2vh[cc] = master->u2vh[c];
    
    if (cc <= window->enonS) {
      windat->eta[cc] = master->eta[c];
      windat->detadt[cc] = master->detadt[c];
      windat->wtop[cc] = master->wtop[c];
      windat->wbot[cc] = master->wbot[c];
      windat->patm[cc] = master->patm[c];
      windat->wind1[cc] = master->wind1[c];
      windat->wind2[cc] = master->wind2[c];
      windat->windspeed[cc] = master->windspeed[c];
      if (master->show_win)
	master->shwin[window->wsa[cc]] = windat->shwin[cc] = (double)n;
      
      if (master->ntrS) {
	for (tn = 0; tn < windat->ntrS; tn++) {
	  windat->tr_wcS[tn][cc] = master->tr_wcS[tn][c];
	}
      }
      if (master->nsed) {
	for (k = 0; k < window->sednz; k++)
	  for (tn = 0; tn < windat->nsed; tn++) {
	    windat->tr_sed[tn][k][cc] = master->tr_sed[tn][k][c];
	  }
      }
    }

    /* e1 face centered arrays */
    windat->u1[cc] = master->u1[c];
    if (cc <= window->enonS) {
      windat->u1av[cc] = master->u1av[c];
    }
    
    /* e2 face centered arrays */
    windat->u2[cc] = master->u2[c];
    if (cc <= window->enonS) {
      windat->u2av[cc] = master->u2av[c];
    }
  }

  /* Save the bottom velocity for the 2D mode */
  for (cc = 1; cc <= window->b2_e1; cc++) {
    c = window->bot_e1[cc];  /* 3D bottom coordinate */
    cs = window->m2d[c];   /* 2D coordiate corresponding to c */
    windat->u1bot[cs] = windat->u1[c] - windat->u1av[cs];
  }
  for (cc = 1; cc <= window->b2_e2; cc++) {
    c = window->bot_e2[cc];  /* 3D bottom coordinate */
    cs = window->m2d[c];   /* 2D coordiate corresponding to c */
    windat->u2bot[cs] = windat->u2[c] - windat->u2av[cs];
  }

  /* Set the leapfrog arrays for the first iteration */
  for (cc = 1; cc <= window->b3_e1; cc ++) {
    c = window->w3_e1[cc];
    windat->u1b[c] = windat->u1[c];
  }
  for (cc = 1; cc <= window->b2_e1; cc ++) {
    c = window->w2_e1[cc];
    windat->u1avb[c] = windat->u1av[c];
  }
  for (cc = 1; cc <= window->b3_e2; cc ++) {
    c = window->w3_e2[cc];
    windat->u2b[c] = windat->u2[c];
  }
  for (cc = 1; cc <= window->b2_e2; cc ++) {
    c = window->w2_e2[cc];
    windat->u2avb[c] = windat->u2av[c];
  }
  for (cc = 1; cc <= window->b2_t; cc ++) {
    c = window->w2_t[cc];
    windat->etab[c] = wincon->oldeta[c] = windat->eta[c];
  }

  set_map_inside(window);
  set_dz(window, windat, wincon);
  get_depths(window, windat, wincon);
  density_w(window, windat, wincon);
  OBC_bgz_nograd(window);
  compute_ref_density(master, window, windat, wincon);
  Set_lateral_BC_density_w(windat->dens, window->nbpt,
			   window->bpt, window->bin);
  wincon->calc_closure(window, windat, wincon);

  /* Set the lateral boundary conditions for velocity.               */
#if !GLOB_BC
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  vel2D_lbc(windat->u2, window->nbpte2, window->nbe2,
	    window->bpte2, window->bine2, wincon->slip);
  vel2D_lbc(windat->u1b, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  vel2D_lbc(windat->u2b, window->nbpte2, window->nbe2,
	    window->bpte2, window->bine2, wincon->slip);
  vel2D_lbc(windat->u1av, window->nbpte1S, window->nbe1S,
	    window->bpte1S, window->bine1S, wincon->slip);
  vel2D_lbc(windat->u2av, window->nbpte2S, window->nbe2S,
	    window->bpte2S, window->bine2S, wincon->slip);
  set_lateral_bc_eta(windat->etab, window->nbptS, window->bpt,
		     window->bin, window->bin2, 1);
  set_lateral_bc_eta(wincon->oldeta, window->nbptS, window->bpt,
		     window->bin, window->bin2, 1);
#endif

    win_data_empty_3d(master, window, windat, (VELOCITY|WVEL));
    win_data_empty_2d(master, window, windat, DEPTH);

  if (master->dozoom) {
    Set_lateral_BC_density_w(master->dens, geom->nbpt, geom->bpt,
                             geom->bin);
    global_interp(geom, master->dens, geom->nzin);
  }

  /* Transfer the backward arrays to auxiliary cells */
  /* Slave to master */
  s2m_vel(master->u1avb,windat->u1avb,
	  window->s2m,window->s2me1,window->ns2mS);
  s2m_vel(master->u2avb,windat->u2avb,
	  window->s2m,window->s2me2,window->ns2mS); 

  /* Master to slave */
  for (cc = 1; cc <= window->nm2sS; cc++) {
    int lc = window->m2s[cc];
    int ce1 = window->m2se1[cc];
    int ce2 = window->m2se2[cc];
    windat->u1avb[lc]=master->u1avb[ce1];
    windat->u2avb[lc]=master->u2avb[ce2]; 
  }
}

/* END window_reset()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Transfers the density profile at the deepest point in the domain  */
/* to the window data.                                               */
/*-------------------------------------------------------------------*/
void copy_dens_profile(master_t *master, geometry_t *window, window_t *windat)
{
  geometry_t *geom = master->geom;
  int c, cb;
  if (windat->vd == NULL) return;

  for (c = 0; c < window->nz; c++)
    windat->vd[c] = 0.0;

  /* Get the profile at the deepest point                            */
  c = geom->cdeep;
  while (c != geom->zp1[c]) {
    windat->vd[geom->s2k[c]] = master->dens_0[c];
    c = geom->zp1[c];
  }
  windat->vd[geom->s2k[c]] = master->dens_0[c];
}

/* END copy_dens_profile()                                           */
/*-------------------------------------------------------------------*/


void  tr_warn(char *text, double w_m, double w_n, int c,
	      int i, int j, int k, int cg, int n, char *type, char tp)
{
  double eps = 1e-10;
  if (fabs(w_m - w_n) > eps) {
    if (k < 0)
      hd_warn("%s mismatch at %d(%d %d)%d : %d%s %c : %8.6e vs %8.6e\n",
	      text, c, i, j, cg, n, type, tp, w_m, w_n);
    else 
      hd_warn("%s mismatch at %d(%d %d %d)%d : %d%s %c : %8.6e vs %8.6e\n",
	      text, c, i, j, k, cg, n, type, tp, w_m, w_n);
  }
}

void check_s2m(geometry_t **window, window_t **windat)
{
  int n, c, cc, lc;
  double *a = d_alloc_1d(geom->sgsiz);
  int nvec, *vecm, *vec;
  int mode = 1;

  for (n = 1; n <= master->nwindows; n++) {
    if (mode == 0) {
      nvec = window[n]->b3_t;
      vec = window[n]->w3_t;
    } else if (mode == 1) {
      nvec = window[n]->b3_e1;
      vec = window[n]->w3_e1;
    } else if (mode == 2) {
      nvec = window[n]->b3_e2;
      vec = window[n]->w3_e2;
    }

    /* Fill all wet cells from the slaves */
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = window[n]->wsa[lc];
      a[c] = windat[n]->dzu1[lc];
    }

    if (mode == 0) {
      nvec = window[n]->ns2m;
      vec = window[n]->s2m;
      vecm = window[n]->wsa;
    } else if (mode == 1) {
      nvec = window[n]->ns2m;
      vec = window[n]->s2m;
      vecm = window[n]->s2me1;
    } else if (mode == 2) {
      nvec = window[n]->ns2m;
      vec = window[n]->s2m;
      vecm = window[n]->s2me2;
    }

    /* Compare transfer cells with previously copied wet cells */
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = vecm[cc];
      if (mode==0) c = vecm[lc];
      if(a[c] != master->dzu1[c])
	printf("Mismatch : %f %d (%d %d %d)%d %f %f\n",master->t/86400,c,
	       window[n]->s2i[lc],window[n]->s2j[lc],window[n]->s2k[lc],n,
	       a[c],master->dzu1[c]);

    }
  }
  d_free_1d(a);
}


