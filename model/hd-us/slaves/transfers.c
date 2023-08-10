/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/slaves/transfers.c
 *  
 *  Description:
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: transfers.c 7377 2023-07-26 04:35:37Z her127 $
 *
 */

#include <stdlib.h>
#include <string.h>
#include "hd.h"
#include "tracer.h"

void  tr_warn(char *text, double w_m, double w_n, int c,
	      int i, int j, int k, int cg, int n, char *type, char tp);
void copy_dens_profile(master_t *master, geometry_t *window, window_t *windat);
int is_ghost(geometry_t *geom, int c);

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
/* s2m_vel()     : Transfers a local velocity array to the master    */
/* s2c_2d()      : Sparse to Cartesian 2D transfers                  */
/* s2c_3d()      : Sparse to Cartesian 3D transfers                  */
/* c2s_2d()      : Cartesian to sparse 2D transfers                  */
/* c2s_3d()      : Cartesian to sparse 3D transfers                  */
/*-------------------------------------------------------------------*/


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
  int e, ee;                    /* Edge locations                    */
  int tn, tt;                   /* Tracer counter                    */
  int s, ce1;     

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
  windat->swan_hs = master->swan_hs;

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
    windat->patm[cc] = master->patm[c];
    for (tt = 0; tt < master->nres2d; tt++) {
      tn = master->reset2d[tt];
      windat->tr_wcS[tn][cc] = master->tr_wcS[tn][c];
    }
    if (!(master->do_wave & NONE) || master->orbital) {
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
      if (windat->wave_k) windat->wave_k[cc] = master->wave_k[c];
    }
    if (master->waves & NEARSHORE) {
      if (windat->wave_Kb) windat->wave_Kb[cc] = master->wave_Kb[c];
      if (windat->wave_k) windat->wave_k[cc] = master->wave_k[c];
      if (windat->wave_fwcapx) windat->wave_fwcapx[cc] = master->wave_fwcapx[c];
      if (windat->wave_fwcapy) windat->wave_fwcapy[cc] = master->wave_fwcapy[c];
      if (windat->wave_fbrex) windat->wave_fbrex[cc] = master->wave_fbrex[c];
      if (windat->wave_fbrey) windat->wave_fbrey[cc] = master->wave_fbrey[c];
      if (windat->wave_fbotx) windat->wave_fbotx[cc] = master->wave_fbotx[c];
      if (windat->wave_fboty) windat->wave_fboty[cc] = master->wave_fboty[c];
      if (windat->wave_fsurx) windat->wave_fsurx[cc] = master->wave_fsurx[c];
      if (windat->wave_fsury) windat->wave_fsury[cc] = master->wave_fsury[c];
      if (windat->wave_wfdx) windat->wave_wfdx[cc] = master->wave_wfdx[c];
      if (windat->wave_wfdy) windat->wave_wfdy[cc] = master->wave_wfdy[c];
      if (windat->wave_wovsx) windat->wave_wovsx[cc] = master->wave_wovsx[c];
      if (windat->wave_wovsy) windat->wave_wovsy[cc] = master->wave_wovsy[c];
      if (windat->wave_frolx) windat->wave_frolx[cc] = master->wave_frolx[c];
      if (windat->wave_froly) windat->wave_froly[cc] = master->wave_froly[c];
    }
    if (master->waves & SPECTRAL) {
      for (tt = 0; tt < master->nsfr; tt++) {
	double f = 2.0 * PI * master->freq[tt];
        windat->freq[tt] = 2.0 * (f * f)/ master->g;
      }
    }
  }
  for (ee = 1; ee < window->szeS; ee++) {
    e = window->wse[ee];
    windat->wind1[ee] = master->wind1[e];
    windat->wind2[ee] = master->wind2[e];
    windat->windspeed[ee] = master->windspeed[e];
    windat->winddir[ee] = master->winddir[e];
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
  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require      */
  /* auxiliary cells only to be transfered to the slave since the    */
  /* wet cells already exist on the slave.                           */
  /* 3D variables.                                                   */
  for (cc = 1; cc <= window->nm2s; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    windat->w[lc] = master->w[c];
    /*
    if (master->means & TRANSPORT) {
      windat->u1m[lc] = master->u1m[c];
      windat->u2m[lc] = master->u2m[c];
      windat->wm[lc] = master->wm[c];
    }
    */
    windat->Kz[lc] = master->Kz[c];
    windat->Vz[lc] = master->Vz[c];
    windat->dens[lc] = master->dens[c];
    for (tn = 0; tn < windat->ntr; tn++) {
      windat->tr_wc[tn][lc] = master->tr_wc[tn][c];
    }
  }
  for (cc = 1; cc <= window->nm2se1; cc++) {
    lc = window->m2se1[cc];
    ce1 = window->wse[lc];
    windat->u1[lc] = master->u1[ce1];
    windat->u2[lc] = master->u2[ce1];
    windat->u1b[lc] = master->u1b[ce1];
    windat->u2b[lc] = master->u2b[ce1];
    windat->u1flux3d[lc] = master->u1flux3d[ce1];
    if (master->means & TRANSPORT) {
      windat->ume[lc] = master->ume[ce1];
      windat->u1vm[lc] = master->u1vm[ce1];
    }
  }

  /* These fluxes do not need to be transferred unless diagnostic    */
  /* checks are required for continuity across windows.              */
  /*
  zflux_e1(geom,window,master->u1flux3d,windat->u1flux3d,window->m2se1,
	   window->wse,window->nm2se1);
  */
  /*-----------------------------------------------------------------*/
  /* 2D variables.                                                   */
  for (cc = 1; cc <= window->nm2sS; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    ce1 = window->m2se1[cc];
    windat->wtop[lc] = master->wtop[c];
    windat->wbot[lc] = master->wbot[c];
    /* The updated value of eta is required at the auxiliary cells   */
    /* to set dzu1, used in the calculation of the advection terms.  */
    windat->eta[lc] = master->eta[c];
  }
  for (cc = 1; cc <= window->nm2se1S; cc++) {
    lc = window->m2se1[cc];
    ce1 = window->wse[lc];
    windat->u1bot[lc] = master->u1bot[ce1];
  }

  if (master->decf & DEC_ETA) {
    for (cc = 1; cc <= window->b2_t; cc++) {
      lc = window->w2_t[cc];
      c = window->wsa[lc];
      windat->decv1[lc] = master->decv1[c];
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

/* END win_data_fill_3d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master (auxiliary cells only) for the 2D mode. If only one window */
/* exists the master points to the windat data structure - no need   */
/* to fill.                                                          */
/*-------------------------------------------------------------------*/
void win_data_fill_2d(master_t *master,     /* Master data           */
                      geometry_t *window,   /* Window geometry       */
                      window_t *windat,     /* Window data           */
                      int nwindows          /* Number of windows     */
  )
{
  int c, cc, lc, ce1;           /* Local sparse coordinate / counter */

  windat->dtb2 = windat->dtf2;
  windat->dtf2 = master->dt2d;
  windat->dt2d = windat->dtf2 + windat->dtb2;

  if (master->mode2d) {
    windat->dtu1 = master->dtu1;
  }
  if (nwindows == 1) {
    bdry_transfer_eta(master, window, windat);
    bdry_transfer_u1av(master, window, windat);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Wet + auxiliary cells. This include variables that the master   */
  /* updates.                                                        */

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require      */
  /* auxiliary cells only to be transfered to the slave since the    */
  /* wet cells already exist on the slave.                           */
  for (cc = 1; cc <= window->nm2sS; cc++) {
    lc = window->m2s[cc];
    c = window->wsa[lc];
    windat->eta[lc] = master->eta[c];
    /*
    windat->uav[lc] = master->uav[c];
    windat->vav[lc] = master->vav[c];
    */
    windat->detadt[lc] = master->detadt[c];
    windat->etab[lc] = master->etab[c];
    /* Uncomment to check master to slave transfer indicies          */
    /*
    printf("m2se1 %d window%d (%d %d)m->(%d %d)s\n",cc,window->wn,geom->s2i[ce1],
	   geom->s2j[ce1],window->s2i[lc],window->s2j[lc]);
    */

  }
  for (cc = 1; cc <= window->nm2se1S; cc++) {
    lc = window->m2se1[cc];
    ce1 = window->wse[lc];
    windat->u1av[lc] = master->u1av[ce1];
    windat->u2av[lc] = master->u2av[ce1];
    windat->depth_e1[lc] = master->depth_e1[ce1];
  }
  /* These fluxes do not need to be transferred unless diagnostic    */
  /* checks are required for continuity across windows.              */
  /*
  zflux_e1(geom,window,master->u1flux,windat->u1flux,window->m2se1,
	   window->wse,window->nm2se1S);
  */

  /*-----------------------------------------------------------------*/
  /* Set the auxiliary multi-dt cells at the start of the timestep   */
  for (cc = 1; cc <= window->naux_t; cc++) {
    c = window->aux_t[cc];
    if (c < window->enonS) {
      c = window->wsa[c];
      windat->u1av_as[cc] = master->u1av[c];
    }
  }

  /* Transfer any custom data from the master to the slaves          */
  bdry_transfer_eta(master, window, windat);
  bdry_transfer_u1av(master, window, windat);

}

/* END win_data_fill_2d()                                            */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to fill the window data structure with data from the      */
/* master required to perform the u2av velocity update.              */
/*-------------------------------------------------------------------*/
void win_data_refill_3d(master_t *master,   /* Master data           */
                        geometry_t *window, /* Window geometry       */
                        window_t *windat,   /* Window data           */
                        int nwindows,       /* Number of windows     */
                        int mode            /* Data flag             */
  )
{
  int c, cc, lc, ce1;           /* Local sparse coordinate / counter */
  int e, ee, le;

  /* Set the crash recovery flags if required                        */
  if (master->crf & RS_WINSET) master->crf &= ~RS_WINSET;
  if (master->crf & RS_RESET) master->crf = NONE;

  if (nwindows == 1 && mode & MIXING) {
    bdry_transfer_u1(master, window, windat);
    return;
  }

  /*-----------------------------------------------------------------*/
  /* Auxiliary cells only. This includes variables that require      */
  /* auxiliary cells only to be transfered to the slave since the    */
  /* wet cells already exist on the slave.                           */
  /* mode = MIXING : vertical mixing coefficients and cell           */
  /* thicknesses at the e1 and e2 faces.                             */
  if (mode & MIXING) {
    geometry_t *geom = master->geom;
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->Kz[lc] = master->Kz[c];
      windat->Vz[lc] = master->Vz[c];
    }
    for (ee = 1; ee <= window->nm2se1; ee++) {
      le = window->m2se1[ee];
      e = window->wse[le];
      windat->dzu1[le] = master->dzu1[e];
    }
    for (ee = 1; ee <= window->nm2se1S; ee++) {
      le = window->m2se1[ee];
      e = window->wse[le];
      windat->sur_e1[le] = geom->sur_e1[e];
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
      /*
      if (master->means & TRANSPORT) {
	for (cc = 1; cc <= window->nm2s; cc++) {
	  lc = window->m2s[cc];
	  c = window->wsa[lc];
	  windat->u1m[lc] = master->u1m[c];
	  windat->u2m[lc] = master->u2m[c];
	  windat->wm[lc] = master->wm[c];
	}
      }
      */
    }
    if (master->waves & NEARSHORE) {
      for (cc = 1; cc <= window->nm2sS; cc++) {
	lc = window->m2s[cc];
	c = window->wsa[lc];
	windat->wave_P[lc] = master->wave_P[c];
      }
    }
    /* Transfer any custom data from the master to the slaves        */
    bdry_transfer_u1(master, window, windat);

  }
  /* mode = VELOCITY : depth averaged adjusted velocities for the    */
  /* vertical velocity equation and 3D fluxes for the tracer         */
  /* equation.                                                       */
  if (mode & VELOCITY) {
    for (ee = 1; ee <= window->nm2se1; ee++) {
      le = window->m2se1[ee];
      e = window->wse[le];
      windat->u1[le] = master->u1[e];
      /*windat->u2[le] = master->u2[e];*/
      windat->u1flux3d[le] = master->u1flux3d[e];
      if (master->means & TRANSPORT) {
	windat->ume[le] = master->ume[e];
	windat->u1vm[le] = master->u1vm[e];
      }
      /* Cell thickness is required for tracer horizontal diffusion  */
      windat->dzu1[le] = master->dzu1[e];
    }
    /* Cell centered velocity                                        */
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->u[lc] = master->u[c];
      windat->v[lc] = master->v[c];
    }

    /*
    zflux_e1(geom, window, master->u1flux3d, windat->u1flux3d, window->m2se1,
             window->wse, window->nm2se1);
    */
    /* The updated value of eta is required at the auxiliary cells   */
    /* to set wtop, used in the calculation of the advection terms.  */
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->eta[lc] = master->eta[c];
    }

    /* These fluxes do not need to be transferred unless diagnostic  */
    /* checks are required for continuity across windows.            */
    /*
    zflux_e1(geom,window,master->u1flux,windat->u1flux,window->m2se1,
	     window->wse,window->nm2se1S);
    */
  }
  if (mode & TRACERS) {
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->dens[lc] = master->dens[c];
    }
  }
  /* Refill vertical velocity and grid thickness for FFSL. This is   */
  /* required to project tracers onto a horizontal level (i.e. the   */
  /* z-direction transverse terms) in ghost cells.                   */
  if (mode & WVEL) {
    win_priv_t *wincon = window->wincon;
    for (cc = 1; cc <= window->nm2s; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->w[lc] = master->w[c];
      wincon->dz[lc] = master->dz[c];
    }
    for (cc = 1; cc <= window->nm2sS; cc++) {
      lc = window->m2s[cc];
      c = window->wsa[lc];
      windat->wtop[lc] = master->wtop[c];
      windat->wbot[lc] = master->wbot[c];
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
void win_data_refill_2d(master_t *master,   /* Master data           */
                        geometry_t *window, /* Window geometry       */
                        window_t *windat,   /* Window data           */
                        int nwindows,       /* Number of windows     */
                        int mode            /* Data flag             */
  )
{
  int cc, lc, ce1;              /* Local sparse coordinate / counter */
  int e, ee, le;

  /* mode = NVELOCITY : updated velocities. This is required to      */
  /* maintain continuity by using the correct nuav when setting uav  */
  /* in asselin() at auxiliary cells; uav is subsequently used to    */
  /* compute fluxes in eta_step(). Note: if asselin() is performed   */
  /* in leapfrog_update_2d() then this transfer is not neccessary    */
  /* but the elevation calculation is made using unfiltered          */
  /* velocities.                                                     */
  if (mode & NVELOCITY) {
    int n;
    for (ee = 1; ee <= window->nm2se1S; ee++) {
      le = window->m2se1[ee];
      e = window->wse[le];
      windat->nu1av[le] = master->nu1av[e];
    }
    /* For VERTIN OBCs the depth is reset so copy tangential OBC     */
    /* depths, so that tangential OBC auxiliaries are set in other   */
    /* windows.                                                      */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->bcond_tan2d & VERTIN) {
	for (ee = 1; ee <= open->no2_ta; ee++) {
	  le = open->obc_ta[ee];
	  e = window->wse[le];
	  windat->depth_e1[le] = master->depth_e1[e];
	}
      }
    }
  }
  if (mode & DEPTH) {
    for (ee = 1; ee <= window->nm2se1S; ee++) {
      le = window->m2se1[ee];
      e = window->wse[le];
      windat->depth_e1[le] = master->depth_e1[e];
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
                       int mode            /* Data flag              */
  )
{
  int c, cc, lc;                /* Local sparse coordinate / counter */
  int e, ee, le;                /* Edge counter                      */
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
    if (master->errornorm & N_ERRN)
      error_norms_m(master, window, windat, window->wincon);
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
            window->s2me1, window->wse, window->ns2me1);

    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->Kz[c] = windat->Kz[lc];
      master->Vz[c] = windat->Vz[lc];
    }
    for (ee = 1; ee <= window->ns2me1S; ee++) {
      le = window->s2me1[ee];
      e = window->wse[le];
      geom->sur_e1[e] = windat->sur_e1[le];
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
    if (master->waves & NEARSHORE) {
      for (cc = 1; cc <= window->ns2mS; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->wave_P[c] = windat->wave_P[lc];
      }
    }
  }

  /* mode = VELOCITY : depth averaged adjusted velocities for the    */
  /* vertical velocity equation and 3D fluxes for the tracer         */
  /* equation.                                                       */
  if (mode & VELOCITY) {
    s2m_vel(master->u1, windat->u1,
            window->s2me1, window->wse, window->ns2me1);
    /*
    s2m_vel(master->u2, windat->u2,
            window->s2me1, window->wse, window->ns2me1);
    */
    s2m_vel(master->u1flux3d, windat->u1flux3d,
            window->s2me1, window->wse, window->ns2me1);
    s2m_vel(master->u1b, windat->u1b,
            window->s2me1, window->wse, window->ns2me1);
    s2m_vel(master->u2b, windat->u2b,
            window->s2me1, window->wse, window->ns2me1);
    if (master->means & TRANSPORT) {
      s2m_vel(master->u1vm, windat->u1vm,
	      window->s2me1, window->wse, window->ns2me1);
      s2m_vel(master->ume, windat->ume,
	      window->s2me1, window->wse, window->ns2me1);
    }
    /* Set in set_new_cells_ and used in tracer horz diffusion       */
    s2m_vel(master->dzu1, windat->dzu1,
	    window->s2me1, window->wse, window->ns2me1);

    /* Cell centered velocity                                        */
    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->u[c] = windat->u[lc];
      master->v[c] = windat->v[lc];
    }

  }
  /* mode = WVEL : 3D vertical velocity and 2D surface and bottom    */
  /* vertical velocities. Also include the 3D velocity deviations    */
  /* from the depth averaged velocity which was calculated in        */
  /* velocity_adjust() here.                                         */
  if (mode & WVEL) {
    /* Transfer to the master                                        */
    for (cc = 1; cc <= window->ns2m; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->w[c] = windat->w[lc];
    }
    s2m_vel(master->u2, windat->u2,
            window->s2me1, window->wse, window->ns2me1);
    s2m_vel(master->u1bot, windat->u1bot,
            window->s2me1, window->wse, window->ns2me1S);
    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->wtop[c] = windat->wtop[lc];
      master->wbot[c] = windat->wbot[lc];
    }
    if (master->trasc & FFSL) {
      win_priv_t *wincon = window->wincon;
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->dz[c] = wincon->dz[lc];
      }
    }

    /*
    if (master->means & TRANSPORT) {
      for (cc = 1; cc <= window->ns2m; cc++) {
	lc = window->s2m[cc];
	c = window->wsa[lc];
	master->u1m[c] = windat->u1m[lc];
	master->u2m[c] = windat->u2m[lc];
	master->wm[c] = windat->wm[lc];
      }
    }
    */
  }
  /* mode = TRACERS : variables calculated in the tracer routine     */
  if (mode & TRACERS) {
    master->wclk[window->wn] += windat->wclk;
    memcpy(master->trinc, windat->trinc, master->ntr * sizeof(int));
    if (master->means & RESET) master->means &= ~RESET;
    /* Transfer to the master                                        */
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
    /* All wet cells (cell centered). Include relaxation and reset   */
    /* tracers since tracers may be altered by updates. Include      */
    /* OBC cells.                                                    */
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
      /* Cell centered velocity for Lagrange    tracking             */
      master->u[c] = windat->u[lc];
      master->v[c] = windat->v[lc];
      master->w[c] = windat->w[lc];
    }
    if (master->heatflux & ADVANCED) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	lc = window->w2_t[cc];
	c = window->wsa[lc];
	master->tr_wc[master->tno][c] = windat->tr_wc[windat->tno][lc];
	master->tr_wc[master->sno][c] = windat->tr_wc[windat->sno][lc];
	master->dens[c] = windat->dens[lc];
      }
    }
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
    /* Wave variables on the slaves are not updated; no need to      */
    /* transfer back to the master.                                  
    if (master->waves & NEARSHORE) {
      for (cc = 1; cc <= window->b2_t; cc++) {
	lc = window->w2_t[cc];
	c = window->wsa[lc];
	if (master->wave_k) master->wave_k[c] = windat->wave_k[lc];
	if (master->wave_Kb) master->wave_Kb[c] = windat->wave_Kb[lc];
	if (master->wave_P) master->wave_P[c] = windat->wave_P[lc];
	if (master->wave_fwcapx) master->wave_fwcapx[cc] = windat->wave_fwcapx[c];
	if (master->wave_fwcapy) master->wave_fwcapy[cc] = windat->wave_fwcapy[c];
	if (master->wave_fbrex) master->wave_fbrex[cc] = windat->wave_fbrex[c];
	if (master->wave_fbrey) master->wave_fbrey[cc] = windat->wave_fbrey[c];
	if (master->wave_fbotx) master->wave_fbotx[cc] = windat->wave_fbotx[c];
	if (master->wave_fboty) master->wave_fboty[cc] = windat->wave_fboty[c];
	if (master->wave_fsurx) master->wave_fsurx[cc] = windat->wave_fsurx[c];
	if (master->wave_fsury) master->wave_fsury[cc] = windat->wave_fsury[c];
	if (master->wave_wfdx) master->wave_wfdx[cc] = windat->wave_wfdx[c];
	if (master->wave_wfdy) master->wave_wfdy[cc] = windat->wave_wfdy[c];
	if (master->wave_wovsx) master->wave_wovsx[cc] = windat->wave_wovsx[c];
	if (master->wave_wovsy) master->wave_wovsy[cc] = windat->wave_wovsy[c];
	if (master->wave_frolx) master->wave_frolx[cc] = windat->wave_frolx[c];
	if (master->wave_froly) master->wave_froly[cc] = windat->wave_froly[c];
      }
    }
    */
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
      /* Get the surface sparse coordinates                          */
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
    /* Mean iteration counter                                        */
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
void win_data_empty_2d(master_t *master,    /* Master data           */
                       geometry_t *window,  /* Window geometry       */
                       window_t *windat,    /* Window data           */
                       int mode             /* Data flag             */
  )
{
  int c, cc, lc;                /* Local sparse coordinate / counter */

  /*-----------------------------------------------------------------*/
  /* Cell centered 2D arrays.                                        */
  /* Variables that require wet cells which are used by other        */
  /* windows to be transfered to the master.                         */
  /* mode = VELOCITY : updated velocity and elevation data.          */
  if (mode & VELOCITY) {
    s2m_vel(master->u1av, windat->u1av,
	    window->s2me1, window->wse, window->ns2me1S);
    s2m_vel(master->u2av, windat->u2av,
	    window->s2me1, window->wse, window->ns2me1S);

    for (cc = 1; cc <= window->ns2mS; cc++) {
      lc = window->s2m[cc];
      c = window->wsa[lc];
      master->eta[c] = windat->eta[lc];
      /*
      master->uav[c] = windat->uav[lc];
      master->vav[c] = windat->vav[lc];
      */
      master->detadt[c] = windat->detadt[lc];
      master->etab[c] = windat->etab[lc];
    }
    memcpy(master->trincS, windat->trincS, master->ntrS * sizeof(int));
  }
  /* mode = DEPTH : total depth at the cell faces                    */
  if (mode & DEPTH) {
    s2m_vel(master->depth_e1, windat->depth_e1,
            window->s2me1, window->wse, window->ns2me1S);
  }
  /* mode = NVELOCITY : updated velocity. Used in asselin().         */
  if (mode & NVELOCITY) {
    int n, ee, e, le;
    s2m_vel(master->nu1av, windat->nu1av,
	    window->s2me1, window->wse, window->ns2me1S);
    /* For VERTIN OBCs the depth is reset so copy tangential OBC     */
    /* depths, so that tangential OBC auxiliaries are set in other   */
    /* windows.                                                      */
    for (n = 0; n < window->nobc; n++) {
      open_bdrys_t *open = window->open[n];
      if (open->bcond_tan2d & VERTIN) {
	for (ee = open->no3_e1 + 1; ee <= open->to2_e1; ee++) {
	  le = open->obc_e1[ee];
	  e = window->wse[le];
	  master->depth_e1[e] = windat->depth_e1[le];
	}
      }
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

  /* Summary                                                         */
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
    master->cu1 = window->wse[windat->cu1[0]];
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
    master->cu1a = window->wse[windat->cu1a[0]];
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
    master->ca1 = window->wse[windat->ca1];
  }
  if (windat->mh1 > master->mh1) {
    master->mh1 = windat->mh1;
    master->ch1 = window->wse[windat->ch1];
  }
  if (windat->mv1 > master->mv1) {
    master->mv1 = windat->mv1;
    master->cv1 = window->wse[windat->cv1];
  }
  if (windat->mb1 > master->mb1) {
    master->mb1 = windat->mb1;
    master->cb1 = window->wse[windat->cb1];
  }
  if (windat->md1 > master->md1) {
    master->md1 = windat->md1;
    master->cd1 = window->wse[windat->cd1];
  }
  if (windat->mc1 > master->mc1) {
    master->mc1 = windat->mc1;
    master->cc1 = window->wse[windat->cc1];
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
  int c, cc, cl, ce1;           /* Sparse counters                   */
  int ee, e, ge;                /* Edge counters                     */
  int n;                        /* Window counter                    */
  int use_ghosts = 1;           /* Flag to include lateral ghosts    */
  int verbose = 0;

  /* If lateral boundary conditions are set in the window then don't */
  /* include the ghost cells in the master-slave transfer vectors.   */
  /* The transfers will operate OK if they are included (with some   */
  /* loss of speed).                                                 */
#if !GLOB_BC
  use_ghosts = 0;
#endif

  /*-----------------------------------------------------------------*/
  /* Set the ghost cells in the window mask velocity ghost cells     */
  /* Western edge ghost cells for e1 cells to process. These cells   */
  /* are actually wet and defined to be in a window but due to the   */
  /* stagger for e1 velocities at western boundaries they are ghost  */
  /* cells for u1. These are identified if the west[c] is non-zero   */
  /* (i.e. a wet cell in a window) and west[west[c]] is zero (i.e. a */
  /* ghost cell).                                                    */

  /* Count the number of auxiliary cells in window wn. Note : if     */
  /* GLOB_BC=1 then ghost cells are also required to be included     */
  /* as their values are set on the master.                          */
  window[wn]->nm2s = window[wn]->nm2sS = 0;
  for (cc = 1; cc <= window[wn]->enon; cc++) {
    cl = window[wn]->m2d[cc];
    c = window[wn]->wsa[cc];
    if (geom->fm[c].wn != wn) {
      /* Continue if use_ghosts=0 and fm[c].wn=0 (c=ghost cell)      */
      if (!use_ghosts && !geom->fm[c].wn)
	continue;

      window[wn]->nm2s += 1;
      if (cc <= window[wn]->enonS)
	window[wn]->nm2sS += 1;
    }
  }
  /* OBC ghost cells adjacent to auxiliary cells                     */
  window[wn]->m2s = i_alloc_1d(window[wn]->nm2s + 1);
 
  /* Fill the master to slave transfer vector                        */
  window[wn]->nm2s = 1;
  for (cc = 1; cc <= window[wn]->enon; cc++) {
    cl = window[wn]->m2d[cc];
    c = window[wn]->wsa[cc];
    if (geom->fm[c].wn != wn) {
      if (!use_ghosts && !geom->fm[c].wn)
	continue;
      window[wn]->m2s[window[wn]->nm2s] = cc;
      window[wn]->nm2s += 1;
    }
  }
  window[wn]->nm2s -= 1;

  /*-----------------------------------------------------------------*/
  /* Get the e1 faces for master to slave transfers                  */
  window[wn]->nm2se1 = window[wn]->nm2se1S = 0;
  for (ee = 1; ee <= window[wn]->n3_e1; ee++) {
    e = window[wn]->w3_e1[ee];
    ge = window[wn]->wse[e];
    if (geom->fm[ge].we != wn) {
      window[wn]->nm2se1 += 1;
      if (e < window[wn]->szeS) {
	window[wn]->nm2se1S += 1;
      }
    }
  }
  window[wn]->m2se1 = i_alloc_1d(window[wn]->nm2se1 + 1);

  window[wn]->nm2se1 = window[wn]->nm2se1S + 1;
  window[wn]->nm2se1S = 1;
  for (ee = 1; ee <= window[wn]->n3_e1; ee++) {
    e = window[wn]->w3_e1[ee];
    ge = window[wn]->wse[e];
    if (geom->fm[ge].we != wn) {
      if (e < window[wn]->szeS) {
	window[wn]->m2se1[window[wn]->nm2se1S] = e;
	if (verbose) printf("m2s wn=%d ee=%d e=%d(%d)\n",wn, window[wn]->nm2se1S, e, ge);
	window[wn]->nm2se1S += 1;
      } else {
	window[wn]->m2se1[window[wn]->nm2se1] = e;
	window[wn]->nm2se1 += 1;
      }
    }
  }
  window[wn]->nm2se1 -= 1;
  window[wn]->nm2se1S -= 1;

  /*-----------------------------------------------------------------*/
  /* Count the number of wet cells in window wn which are auxiliary  */
  /* in some other window. If lateral BC's are set on the master     */
  /* then the transfer vectors must include may values at the cell   */
  /* interior to the ghost cell (e.g. free slip condition on         */
  /* tangential velocity components).                                */
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
      if (geom->fm[c].wn == wn) {
        cl = geom->fm[c].sc;
        window[wn]->ns2m += 1;
        if (cl <= window[wn]->enonS) {
          window[wn]->ns2mS += 1;
        }
      }
    }
  }
  window[wn]->s2m = i_alloc_1d(window[wn]->ns2m + 1);

  /* Fill the slave to master transfer vector                        */
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
  window[wn]->ns2m -= 1;
  window[wn]->ns2mS -= 1;

  /*-----------------------------------------------------------------*/
  /* Get the e1 faces for slave to master transfers                  */
  window[wn]->ns2me1 = window[wn]->ns2me1S = 0;
  for (n = 1; n <= nwindows; n++) {
    if (n == wn)
      continue;
    for (ee = 1; ee <= window[n]->n3_e1; ee++) {
      e = window[n]->w3_e1[ee];
      ge = window[n]->wse[e];
      if (geom->fm[ge].we == wn) {
	e = geom->g2we[wn][ge];
        window[wn]->ns2me1 += 1;
        if (e < window[wn]->szeS) {
          window[wn]->ns2me1S += 1;
        }
      }
    }
  }
  window[wn]->s2me1 = i_alloc_1d(window[wn]->ns2me1 + 1);

  /* Fill the slave to master transfer vector                        */
  window[wn]->ns2me1 = window[wn]->ns2me1S + 1;
  window[wn]->ns2me1S = 1;

  for (n = 1; n <= nwindows; n++) {
    if (n == wn)
      continue;
    for (ee = 1; ee <= window[n]->n3_e1; ee++) {
      e = window[n]->w3_e1[ee];
      ge = window[n]->wse[e];
      if (geom->fm[ge].we == wn) {
	e = geom->g2we[wn][ge];
        if (e < window[wn]->szeS) {
          window[wn]->s2me1[window[wn]->ns2me1S] = e;
	  if (verbose) printf("s2m wn=%d ee=%d e=%d(%d)\n",wn, window[wn]->ns2me1S, e, ge);
          window[wn]->ns2me1S += 1;
        } else {
          window[wn]->s2me1[window[wn]->ns2me1] = e;
          window[wn]->ns2me1 += 1;
        }
      }
    }
  }
  window[wn]->ns2me1 -= 1;
  window[wn]->ns2me1S -= 1;

  if (DEBUG("init_w")) {
    dlog("init_w",
         "  Window %d : %d (%d 2D) master to slave transfer cells\n", wn,
         window[wn]->nm2s, window[wn]->nm2sS);
    dlog("init_w",
         "  Window %d : %d (%d 2D) slave to master transfer cells\n", wn,
         window[wn]->ns2m, window[wn]->ns2mS);
  }
}

/* END build_transfer_maps()                                         */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to populate the master with window data                   */
/*-------------------------------------------------------------------*/
void master_fill(master_t *master,    /* Master data                 */
                 geometry_t **window, /* Window geometry             */
                 window_t **windat,   /* Window data                 */
                 win_priv_t **wincon  /* Window constants            */
  )
{
  int n;                              /* Window counter              */

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
  /* Counters                                                        */
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
	master->uav[c] = windat[n]->uav[lc];
	master->vav[c] = windat[n]->vav[lc];
	for (tt = 0; tt < master->ntrS; tt++) {
	  if (master->do_wave & W_SWAN && master->trinfo_2d[tt].type & WAVE) continue;
	  master->tr_wcS[tt][c] = windat[n]->tr_wcS[tt][lc];
	}
      }
      master->u[c] = windat[n]->u[lc];
      master->v[c] = windat[n]->v[lc];
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

  /* Get the cell the glider resides in                              */
  cg = get_glider_loc(master, ts, &ts->ts, t);

  /* Get a neighbourhood around the glider cell                      */
  ssize = ts->kernal;
  st = stencil(geom, cg, &ssize, ST_SIZED, 0);

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
void s2m_3d(master_t *master,   /* Master data                       */
            geometry_t *window, /* Window geometry                   */
            window_t *windat,   /* Window data                       */
            win_priv_t *wincon  /* Private window data               */
  )
{
  int cc, c, lc;                /* Cell centers / counters           */
  int ee, e, le;                /* Edge centers / counters           */
  int n, tn, tt;                /* Counters                          */

  /* Cell centre variables                                           */
  for (cc = 1; cc <= window->b3_t; cc++) {
    lc = window->w3_t[cc];
    c = window->wsa[lc];
    
    master->w[c] = windat->w[lc];
    master->dens[c] = windat->dens[lc];
    master->dens_0[c] = windat->dens_0[lc];
    master->Vz[c] = windat->Vz[lc];
    master->Kz[c] = windat->Kz[lc];
    master->u[c] = windat->u[lc];
    master->v[c] = windat->v[lc];
    master->dz[c] = wincon->dz[lc];
    master->u1kh[c] = wincon->u1kh[lc];
    for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
      tn = master->trmap_s2m_3d[tt];
      master->tr_wc[tn][c] = windat->tr_wc[tn][lc];
    }
  }
  
  /* Edge centered variables                                        */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    le = window->w3_e1[ee];
    e = window->wse[le];
    master->u1[e] = windat->u1[le];
    master->u2[e] = windat->u2[le];
    master->dzu1[e] = windat->dzu1[le];
    master->u1vh[e] = wincon->u1vh[le];
  }
  if (master->means & VEL3D) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      le = window->w3_e1[ee];
      e = window->wse[le];
      master->ume[e] = windat->ume[le];
    }
  }
  if (master->means & VOLFLUX) {
    for (ee = 1; ee <= window->b3_e1; ee++) {
      le = window->w3_e1[ee];
      e = window->wse[le];
      master->u1vm[e] = windat->u1vm[le];
    }
  }
  if (master->means & VEL2D) {
    for (ee = 1; ee <= window->b2_e1; ee++) {
      le = window->w2_e1[ee];
      e = window->wse[le];
      master->uame[e] = windat->uame[le];
    }
  }

  /* Fill the R_EDGE and F_EDGE cell centres if required            */
  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    for (ee = 1; ee <= open->no3_e1; ee++) {
      c = open->obc_e2[ee];
      if (open->ceni[ee]) {
	for (tt = 0; tt < master->ntrmap_s2m_3d; tt++) {
	tn = master->trmap_s2m_3d[tt];
	if (master->trinfo_3d[tn].type & E1VAR) {
	    lc = open->ogc_t[ee];
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
void s2m_2d(master_t *master,   /* Master data                       */
            geometry_t *window, /* Window geometry                   */
            window_t *windat    /* Window data                       */
	    )
{
  int cc, c, lc;                /* Cell centers / counters           */
  int ee, e, le;                /* Edge centers / counters           */
  int cb, tn, k = 0; 

  /* Grid centered variables. Note: wind1, patm and Cd are           */
  /* already on the master since these variables are read onto the   */
  /* master and transferred to the window at the start of the step.  */
  for (cc = 1; cc <= window->b2_t; cc++) {
    lc = window->w2_t[cc];
    c = window->wsa[lc];
    master->eta[c] = windat->eta[lc];
    master->topz[c] = windat->topz[lc];
    master->wtop[c] = windat->wtop[lc];
    master->uav[c] = windat->uav[lc];
    master->vav[c] = windat->vav[lc];
    for (tn = 0; tn < master->ntrS; tn++) {
      if (master->do_wave & W_SWAN && master->trinfo_2d[tn].type & WAVE) continue;
      master->tr_wcS[tn][c] = windat->tr_wcS[tn][lc];
    }
    for (tn = 0; tn < windat->nsed; tn++) {
      for (k = 0; k < window->sednz; k++)
        master->tr_sed[tn][k][c] = windat->tr_sed[tn][k][lc];
    }
  }
  
  /* Edge centered variables                                         */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    le = window->w2_e1[ee];
    cb = window->wse[window->bot_e1[ee]];
    e = window->wse[le];
    master->u1av[e] = windat->u1av[le];
    master->u2av[e] = windat->u2av[le];
    /* Re-setting the master for printing purposes can interfere     */
    /* with master-slave transfers.                                  */
    /* master->u1bot[e]=windat->u1[cb]; */
    /* master->u1bot[e] = windat->u1bot[le]; */
  }
}

/* END s2m_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer the local velocities to the master.           */
/*-------------------------------------------------------------------*/
void s2m_vel(double *Ag,        /* Global flux array                 */
             double *Al,        /* Local flux array                  */
             int *vec,          /* Cells on which Al is defined      */
             int *evec,         /* Cells on which Ag is defined      */
             int nvec           /* Size of vec                       */
  )
{
  int c, lc, cc;                /* Sparse coordinate / counter       */

  /* Transfer the velocity to the master                             */
  for (cc = 1; cc <= nvec; cc++) {
    lc = vec[cc];
    c = evec[lc];
    Ag[c] = Al[lc];
  }
}

/* END s2m_vel()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns 1 for ghost cells                                         */
/*-------------------------------------------------------------------*/
int is_ghost(geometry_t *geom, int c)
{
  int i = geom->s2i[c];
  int j = geom->s2j[c];
  int k = geom->s2k[c];

  if (c == geom->m2d[c]) {
    if (i == NOTVALID && j == NOTVALID) return(1);
  } else {
    if (i == NOTVALID && j == NOTVALID && k == NOTVALID) return(1);
  }
  if (geom->wgst[c]) return(1);

  return(0);
}

/* END is_ghost()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void s2c_2d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double **ac,        /* Output cartesian array            */
            int nx,             /* x dimension of 2d array           */
            int ny              /* y dimension of 2d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int c;                        /* Global sparse coordinate          */

  /* Initialise                                                      */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      ac[j][i] = 0.0;

  /* Map the 2d sparse variable to Cartesion                         */
  for (c = 1; c <= geom->ewetS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    if (is_ghost(geom, c)) continue;
    if (i < nx && j < ny)
      ac[j][i] = as[c];
  }
}

/* END s2c_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d vertex arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void v2c_2d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double **ac,        /* Output cartesian array            */
            int nx,             /* x dimension of 2d array           */
            int ny              /* y dimension of 2d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int v, vv;                    /* Global sparse coordinate          */

  /* Initialise                                                      */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      ac[j][i] = 0.0;

  /* Map the 2d sparse variable to Cartesion                         */
  for (vv = 1; vv <= geom->b2_e2; vv++) {
    v = geom->w2_e2[vv];
    i = geom->v2i[v];
    j = geom->v2j[v];
    if (is_ghost(geom, geom->v2ijk[v])) continue; 

    if (i < nx && j < ny)
      ac[j][i] = as[v];
  }
}

/* END s2c_2dv()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d edge arrays into the Cartesian array                     */
/*-------------------------------------------------------------------*/
void e2c_2d(geometry_t *geom,   /* Global geometry                   */
	    double *as,         /* Input sparse array                */
	    double **ac,        /* Output cartesian array            */
	    int nx,             /* x dimension of 2d array           */
	    int ny,             /* y dimension of 2d array           */
	    int mode
	    )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int cc, c, e;                 /* Global sparse coordinate          */

  /* Initialise                                                      */
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      ac[j][i] = 0.0;

  /* Map the 2d edge variable to Cartesion centres                   */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    e = (mode) ? geom->c2e[4][c] : geom->c2e[1][c];
    i = geom->s2i[c];
    j = geom->s2j[c];
    if (is_ghost(geom, c)) continue; 
    if (i < nx && j < ny)
      ac[j][i] = as[e];
  }
}

/* END e2c_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void s2c_3d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double ***ac,       /* Output cartesian array            */
            int nx,             /* x dimension of 3d array           */
            int ny,             /* y dimension of 3d array           */
            int nz              /* z dimension of 3d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int cc;                       /* Local sparse coordinate counter   */
  int c;                        /* Global sparse coordinate          */

  if (dumpdata->vmap != NULL) {
    /* Linearly interpolate onto the geom->layers[] distribution     */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++) {
	  ac[k][j][i] = NaN;
	  c = dumpdata->vmap[k][j][i];
	  if (c > 0) {
	    int zm1 = geom->zm1[c];
	    double d1 = geom->cellz[c] * dumpdata->botz[j][i];
	    double d2 = geom->cellz[zm1] * dumpdata->botz[j][i];
	    ac[k][j][i] = (fabs(geom->layers[k]) - d1) * 
	      (as[c] - as[zm1]) / (d1 - d2) + as[c];
	  }
	}
  } else {
    /* Initialise                                                    */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  ac[k][j][i] = 0.0;

    /* Map the 3d sparse variable to Cartesion                       */
    for (cc = 1; cc <= geom->a3_t; cc++) {
      c = geom->wsa[cc];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      ac[k][j][i] = 0.0;
      if (is_ghost(geom, c)) continue; 
      if (i < nx && j < ny)
	ac[k][j][i] = as[c];
    }
  }
}

/* END s2c_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d edge arrays into the Cartesian array                     */
/*-------------------------------------------------------------------*/
void e2c_3d(geometry_t *geom,   /* Global geometry                   */
	    double *as,         /* Input sparse array                */
	    double ***ac,       /* Output cartesian array            */
	    int nx,             /* x dimension of 3d array           */
	    int ny,             /* y dimension of 3d array           */
	    int nz,             /* z dimension of 3d array           */
	    int mode            /* 0 = e1, 1 = e2                    */
	    )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int cc;                       /* Local sparse coordinate counter   */
  int c, e;                     /* Global sparse coordinate          */

  if (dumpdata->vmap != NULL) {
    /* Linearly interpolate onto the geom->layers[] distribution     */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++) {
	  ac[k][j][i] = NaN;
	  c = dumpdata->vmap[k][j][i];
	  e = (mode) ? geom->c2e[4][c] : geom->c2e[1][c];
	  if (c > 0) {
	    int zme = zm1e(geom, e);
	    int zm1 = geom->zm1[c];
	    double d1 = geom->cellz[c] * dumpdata->botz[j][i];
	    double d2 = geom->cellz[zm1] * dumpdata->botz[j][i];
	    ac[k][j][i] = (fabs(geom->layers[k]) - d1) * 
	      (as[e] - as[zme]) / (d1 - d2) + as[e];
	  }
	}
  } else {
    /* Initialise                                                    */
    for (k = 0; k < nz; k++)
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  ac[k][j][i] = 0.0;

    /* Map the 3d sparse variable to Cartesion                       */
    for (cc = 1; cc <= geom->b3_t; cc++) {
      c = geom->w3_t[cc];
      e = (mode) ? geom->c2e[4][c] : geom->c2e[1][c];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      if (is_ghost(geom, c)) continue; 
      if (i < nx && j < ny)
	ac[k][j][i] = as[e];
    }
  }
}

/* END e2c_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d sparse sediment arrays into the Cartesian array          */
/*-------------------------------------------------------------------*/
void s2c_sed(geometry_t *geom,  /* Global geometry                   */
            double **as,        /* Input sparse array                */
            double ***ac,       /* Output cartesian array            */
            int nx,             /* x dimension of 3d array           */
            int ny,             /* y dimension of 3d array           */
            int nz              /* z dimension of 3d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int cc;                       /* Local sparse coordinate counter   */
  int c;                        /* Global sparse coordinate          */

  /* Initialise                                                      */
  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        ac[k][j][i] = NaN;

  /* Map the 2d sparse variable to Cartesion                         */
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
/* Packs 2d Cartesian arrays into the cell centered array            */
/*-------------------------------------------------------------------*/
void c2s_2d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double **ac,        /* Output cartesian array            */
            int nx,             /* x dimension of 2d array           */
            int ny              /* y dimension of 2d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int c;                        /* Global sparse coordinate          */

  for (c = 1; c <= geom->sgnumS; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    if (is_ghost(geom, c)) continue; 
    if (i < nx && j < ny)
      as[c] = ac[j][i];
  }
}

/* END c2s_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d Cartesian arrays into the edge centered array            */
/*-------------------------------------------------------------------*/
void c2e_2d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double **ac,        /* Output cartesian array            */
            int nx,             /* x dimension of 2d array           */
            int ny,             /* y dimension of 2d array           */
	    int mode            /* 0 = e1, 1 = e2                    */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int e, ee, c;                 /* Global sparse coordinate          */

  for (ee = 1; ee <= geom->n2_e1; ee++) {
    e = geom->w2_e1[ee];
    if (geom->e2e[e][0] == mode) {
      c = geom->e2ijk[e];
      i = geom->s2i[c];
      j = geom->s2j[c];
      if (is_ghost(geom, c)) continue; 
      if (i < nx && j < ny)
	as[e] = ac[j][i];
    }
  }
}

/* END c2e_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 2d Cartesian arrays into the vertex centered array          */
/*-------------------------------------------------------------------*/
void c2v_2d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double **ac,        /* Output cartesian array            */
            int nx,             /* x dimension of 2d array           */
            int ny              /* y dimension of 2d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int v, vv, c;                 /* Global sparse coordinate          */

  for (vv = 1; vv <= geom->b2_e2; vv++) {
    v = geom->w2_e2[vv];
    c = geom->v2ijk[v];
    i = geom->s2i[c];
    j = geom->s2j[c];
    if (is_ghost(geom, c)) continue; 
    if (i < nx && j < ny)
      as[v] = ac[j][i];
  }
}

/* END c2v_2d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d sparse arrays into the Cartesian array                   */
/*-------------------------------------------------------------------*/
void c2s_3d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double ***ac,       /* Output cartesian array            */
            int nx,             /* x dimension of 3d array           */
            int ny,             /* y dimension of 3d array           */
            int nz              /* z dimension of 3d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int c;                        /* Global sparse coordinate          */

  for (c = 1; c <= geom->sgnum; c++) {
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];

    if (is_ghost(geom, c)) continue; 
    if (i < nx && j < ny)
      as[c] = ac[k][j][i];
  }
}

/* END c2s_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d Cartesian arrays into the edge centered array            */
/*-------------------------------------------------------------------*/
void c2e_3d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double ***ac,       /* Output cartesian array            */
            int nx,             /* x dimension of 3d array           */
            int ny,             /* y dimension of 3d array           */
            int nz,             /* z dimension of 3d array           */
	    int mode            /* 0 = e1, 1 = e2                    */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int ee, e, c;                 /* Global sparse coordinate          */

  for (ee = 1; ee <= geom->n3_e1; ee++) {
    e = geom->w3_e1[ee];
    if (geom->e2e[e][0] == mode) {
      c = geom->e2ijk[e];
      i = geom->s2i[c];
      j = geom->s2j[c];
      k = geom->s2k[c];
      if (is_ghost(geom, c)) continue; 
      if (i < nx && j < ny)
	as[c] = ac[k][j][i];
    }
  }
}

/* END c2e_3d()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Packs 3d Cartesian arrays into the vertex centered array          */
/*-------------------------------------------------------------------*/
void c2v_3d(geometry_t *geom,   /* Global geometry                   */
            double *as,         /* Input sparse array                */
            double ***ac,       /* Output cartesian array            */
            int nx,             /* x dimension of 3d array           */
            int ny,             /* y dimension of 3d array           */
            int nz              /* z dimension of 3d array           */
  )
{
  int i;                        /* Counter in the x direction        */
  int j;                        /* Counter in the y direction        */
  int k;                        /* Counter in the z direction        */
  int vv, v, c;                 /* Global sparse coordinate          */

  for (vv = 1; vv <= geom->b3_e2; vv++) {
    v = geom->w3_e2[vv];
    c = geom->v2ijk[v];
    i = geom->s2i[c];
    j = geom->s2j[c];
    k = geom->s2k[c];
    if (is_ghost(geom, c)) continue; 
    if (i < nx && j < ny)
      as[v] = ac[k][j][i];
  }
}

/* END v2e_3d()                                                      */
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
void unpack_sparse(int *map, int mapsize, double *var, double *unpack, int oset)
{
  int cc, c;
  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    unpack[cc] = var[cc-oset];
  }
}

void unpack_sparse3(int *map, int mapsize, double *var, double *unpack, int oset)
{
  int cc, c;

  for (cc = 1; cc <= mapsize; cc++) {
    c = map[cc];
    unpack[c] = var[cc - oset];
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
	/*if (i == NOTVALID && j == NOTVALID && k == NOTVALID) {*/
	if (is_ghost(geom, cg)) {
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
void window_resetm(master_t *master,      /* Master data             */
		   geometry_t **window,   /* Window geometry         */
		   window_t **windat,     /* Window data             */
		   win_priv_t **wincon,   /* Window constants        */
		   int mode
		   )
{
  int c, cc, lc, cs, n, m, k, tn;
  int ee, e, es, vv, v;
  int nwindows = master->nwindows;

  if (mode & RS_VH) {
    for (n = 1; n <= nwindows; n++) {
      for (ee = 1; ee <= window[n]->sze; ee++) {
	e = window[n]->wse[ee];
	wincon[n]->u1vh[ee] = master->u1vh[e];
      }
      master->dtb = master->dtf = master->dt = master->grid_dt;
      master->dt2d = master->grid_dt / master->iratio;
      windat[n]->dtf2 = master->dt / master->iratio;
      return;
    }
  }

  for (n = 1; n <= nwindows; n++) {
    window_reset(master, window[n], windat[n], wincon[n], RS_ALL);
    win_data_fill_3d(master, window[n], windat[n], master->nwindows);
    win_data_fill_2d(master, window[n], windat[n], master->nwindows);
    set_map_inside(window[n]);
    memcpy(wincon[n]->oldeta, windat[n]->eta, window[n]->sgsizS * sizeof(double));
    set_dz(window[n], windat[n], wincon[n]);
    get_depths(window[n], windat[n], wincon[n]);
    density_w(window[n], windat[n], wincon[n]);

    /* Set sponge zones                                              */
    set_sponge_cells(window[n]);
    set_sponge_c(window[n], wincon[n]->u1kh, windat[n]->dt);
    set_sponge_e(window[n], wincon[n]->u1vh, windat[n]->dt);
    set_sponge_e(window[n], wincon[n]->u2kh, windat[n]->dt);    

    /*---------------------------------------------------------------*/
    /* Coriolis at vertices                                          */
    for (vv = 1; vv <= window[n]->n2_e2; vv++) {
      double d1, d2;
      v = window[n]->w2_e2[vv];
      d1 = d2 = 0.0;
      for (m = 1; m <= window[n]->nvc[v]; m++) {
	c = window[n]->v2c[v][m];
	if (c && !window[n]->wgst[c]) {
	  d2 += wincon[n]->coriolis[c] * window[n]->dualareap[v][m];
	  d1 += window[n]->dualareap[v][m];
	}
      }
      if (d1) windat[n]->fv[v] = d2 / d1;
    }
  }

  /* Cell centrered velocity (for restarts)                          */
  for (n = 1; n <= nwindows; n++) {
    memset(windat[n]->u2, 0, window[n]->sze);
    memset(windat[n]->u2av, 0, window[n]->szeS);
    win_data_refill_3d(master, window[n], windat[n], master->nwindows, VELOCITY);
    vel_components_3d(window[n], windat[n], wincon[n]);
    for (ee = 1; ee <= window[n]->nm2se1S; ee++) {
      lc = window[n]->m2se1[ee];
      e = window[n]->wse[lc];
      windat[n]->u1av[lc] = master->u1av[e];
    }
    vel_components_2d(window[n], windat[n], wincon[n]);
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
    set_lateral_bc_eta(windat[n]->etab, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
    set_lateral_bc_eta(wincon[n]->oldeta, window[n]->nbptS, window[n]->bpt,
		       window[n]->bin, window[n]->bin2, 1);
#endif

    win_data_empty_3d(master, window[n], windat[n], (VELOCITY|WVEL));
    win_data_empty_2d(master, window[n], windat[n], DEPTH);
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

  obc_setup(master, window);
  init_trans(master, window, wincon);
  nan_check(window, windat, wincon, nwindows);
  pt_setup(master, window, windat, wincon);
  swr_params_init(master, window);

}

void window_reset(master_t *master,      /* Master data              */
		  geometry_t *window,    /* Window geometry          */
		  window_t *windat,      /* Window data              */
		  win_priv_t *wincon,    /* Window constants         */
		  int mode
		  )
{
  geometry_t *geom = master->geom;
  int c, cc, cs, n, k, tn;
  int ee, e, es;

  if (mode & RS_VH) {
    for (ee = 1; ee <= window->sze; ee++) {
      e = window->wse[ee];
      wincon->u1vh[ee] = master->u1vh[e];
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
    /*wincon->u1kh[cc] = master->u1kh[c];*/
    
    if (cc <= window->enonS) {
      windat->eta[cc] = master->eta[c];
      windat->detadt[cc] = master->detadt[c];
      windat->wtop[cc] = master->wtop[c];
      windat->wbot[cc] = master->wbot[c];
      windat->patm[cc] = master->patm[c];

      if (master->show_win) {
	master->shwin[window->wsa[cc]] = windat->shwin[cc] = (double)window->wn;
      }

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
  }
  /* Edge arrays                                                     */
  for (ee = 1; ee < window->sze; ee++) {
    e = window->wse[ee];
    windat->u1[ee] = master->u1[e];
    wincon->u1vh[ee] = master->u1vh[e];
    wincon->u2kh[ee] = master->u2kh[e];
  }
  for (ee = 1; ee < window->szeS; ee++) {
    e = window->wse[ee];
    windat->u1av[ee] = master->u1av[e];
    windat->wind1[ee] = master->wind1[e];
    windat->wind2[ee] = master->wind2[e];
    windat->windspeed[ee] = master->windspeed[e];
    windat->winddir[ee] = master->winddir[e];
    if (master->u1vhin) windat->u1vhin[ee] = master->u1vhin[e];
    wincon->basev[ee] = master->basev[e];
    wincon->smagv[ee] = master->smagv[e];
    wincon->basek[ee] = master->basek[e];
    wincon->smagk[ee] = master->smagk[e];
  }

  /* Save the bottom velocity for the 2D mode                        */
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->bot_e1[ee];  /* 3D bottom coordinate                 */
    es = window->m2de[e];    /* 2D coordiate corresponding to c      */
    windat->u1bot[es] = windat->u1[e] - windat->u1av[es];
  }

  /* Set the leapfrog arrays for the first iteration                 */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    windat->u1b[e] = windat->u1[e];
  }
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    windat->u1avb[e] = windat->u1av[e];
  }
  for (cc = 1; cc <= window->b2_t; cc ++) {
    c = window->w2_t[cc];
    windat->etab[c] = wincon->oldeta[c] = windat->eta[c];
  }

  set_map_inside(window);
  set_dz(window, windat, wincon);
  get_depths(window, windat, wincon);
  density_w(window, windat, wincon);

  set_sponge_cells(window);
  set_sponge_c(window, wincon->u1kh, windat->dt);
  set_sponge_e(window, wincon->u1vh, windat->dt);
  set_sponge_e(window, wincon->u2kh, windat->dt);    

  /* Cell centrered velocity (for restarts)                          */
  win_data_refill_3d(master, window, windat, master->nwindows, VELOCITY);
  vel_components_3d(window, windat, wincon);
  for (ee = 1; ee <= window->nm2se1S; ee++) {
    int lc = window->m2se1[ee];
    e = window->wse[lc];
    windat->u1av[lc] = master->u1av[e];
  }
  vel_components_2d(window, windat, wincon);

  /*-----------------------------------------------------------------*/
  /* Variable filling and initialisation                             */
  win_data_fill_3d(master, window, windat, master->nwindows);
  win_data_fill_2d(master, window, windat, master->nwindows);
  compute_ref_density(master, window, windat, wincon);
  Set_lateral_BC_density_w(windat->dens, window->nbpt,
			   window->bpt, window->bin);

  windat->dtf2 = master->dt / master->iratio;
  wincon->calc_closure(window, windat, wincon);

  /* Set a no-gradient over ghost OBCs                               */
  OBC_bgz_nograd(window);

  /* Set the leapfrog arrays for the first iteration                 */
  for (ee = 1; ee <= window->b3_e1; ee++) {
    e = window->w3_e1[ee];
    windat->u1b[e] = windat->u1[e];
    windat->u2b[e] = windat->u2[e];
  }
  for (ee = 1; ee <= window->b2_e1; ee++) {
    e = window->w2_e1[ee];
    windat->u1avb[e] = windat->u1av[e];
  }
  for (cc = 1; cc <= window->b2_t; cc ++) {
    c = window->w2_t[cc];
    windat->etab[c] = wincon->oldeta[c] = windat->eta[c];
  }
  /* Set the lateral boundary conditions for velocity.               */
#if !GLOB_BC
  vel2D_lbc(windat->u1, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  vel2D_lbc(windat->u1b, window->nbpte1, window->nbe1,
	    window->bpte1, window->bine1, wincon->slip);
  set_lateral_bc_eta(windat->etab, window->nbptS, window->bpt,
		     window->bin, window->bin2, 1);
  set_lateral_bc_eta(wincon->oldeta, window->nbptS, window->bpt,
		     window->bin, window->bin2, 1);
#endif

  win_data_empty_3d(master, window, windat, (VELOCITY|WVEL));
  win_data_empty_2d(master, window, windat, DEPTH);

  /* Transfer the backward arrays to auxiliary cells                 */
  /* Slave to master                                                 */
  s2m_vel(master->u1avb,windat->u1avb,
	  window->s2me1,window->wse,window->ns2me1S);

  /* Master to slave                                                 */
  for (cc = 1; cc <= window->nm2se1S; cc++) {
    int lc = window->m2se1[cc];
    int ce1 = window->wse[lc];
    windat->u1avb[lc]=master->u1avb[ce1];
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

    /* Fill all wet cells from the slaves                            */
    for (cc = 1; cc <= nvec; cc++) {
      lc = vec[cc];
      c = window[n]->wse[lc];
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

    /* Compare transfer cells with previously copied wet cells       */
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

/*
 * Rank 0 assumes host
 * Rank 1 assumes remote
 */
int mpi_check_multi_windows_sparse_arrays(geometry_t *geom)
{
  int ret = 0;
  
#ifdef HAVE_MPI
  int mpi_rank_other = (mpi_rank == 0 ? 1 : 0);

  /* Only works for exactly 2 processes */
  if (mpi_size != 2)
    return(0);

  /*
   * Multiple windows receives from single window
   */
  if (master->nwindows > 1) {
    /* _o denotes other */
    int v2_e1_o, v3_e1_o, v2_t_o, v3_t_o;
    int szeS_o, sze_o;
    int szcS_o, szc_o;
    int *w2_e1_o, *w3_e1_o, *w2_t_o, *w3_t_o;
    int ee, cc;
    
    /* Receive total number of 2D surface edges */
    MPI_Recv(&szeS_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->szeS != szeS_o) {
      printf("szeS mismatch : %d(multi) vs %d(single)\n", geom->szeS, szeS_o);
      ret = 1;
    } else
      printf("szeS = %d check\n", geom->szeS);

    /* Receive total number of 3D edges */
    MPI_Recv(&sze_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->sze != sze_o) {
      printf("sze mismatch : %d(multi) vs %d(single)\n", geom->sze, sze_o);
      ret = 1;
    } else
      printf("sze  = %d check\n", geom->sze);

    /* Receive total number of 2D surface faces */
    MPI_Recv(&szcS_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->szcS != szcS_o) {
      printf("szcS mismatch : %d(multi) vs %d(single)\n", geom->szcS, szcS_o);
      ret = 1;
    } else
      printf("szcS = %d check\n", geom->szcS);

    /* Receive total number of 3D faces */
    MPI_Recv(&szc_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->szc != szc_o) {
      printf("szc mismatch : %d(multi) vs %d(single)\n", geom->szc, szc_o);
      ret = 1;
    } else
      printf("szc  = %d check\n", geom->szc);

    /* Receive number of 2D edges */
    MPI_Recv(&v2_e1_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->v2_e1 != v2_e1_o) {
      printf("v2_e1 mismatch : %d(multi) vs %d(single)\n", geom->v2_e1, v2_e1_o);
      ret = 1;
    } else
      printf("v2_e1 = %d check\n", geom->v2_e1);

    /* Receive number of 3D edges */
    MPI_Recv(&v3_e1_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->v3_e1 != v3_e1_o) {
      printf("v3_e1 mismatch : %d(multi) vs %d(single)\n", geom->v3_e1, v3_e1_o);
      ret = 1;
    } else
      printf("v3_e1 = %d check\n", geom->v3_e1);

    /* Receive number of 2D faces */
    MPI_Recv(&v2_t_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->v2_t != v2_t_o) {
      printf("v2_t mismatch : %d(multi) vs %d(single)\n", geom->v2_t, v2_t_o);
      ret = 1;
    } else
      printf("v2_t = %d check\n", geom->v2_t);

    /* Receive number of 3D faces */
    MPI_Recv(&v3_t_o, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (geom->v3_t != v3_t_o) {
      printf("v3_t mismatch : %d(multi) vs %d(single)\n", geom->v3_t, v3_t_o);
      ret = 1;
    } else
      printf("v3_t = %d check\n", geom->v3_t);

    /* Send acknowledgement */
    MPI_Send(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    if (ret) return(ret);

    /* Sizes match, allocate arrays */
    w2_e1_o = i_alloc_1d(geom->szeS);
    w3_e1_o = i_alloc_1d(geom->sze);
    w2_t_o  = i_alloc_1d(geom->szcS);
    w3_t_o  = i_alloc_1d(geom->szc);

    /* Receive sparse arrays */
    MPI_Recv(&w2_e1_o[1], geom->szeS-1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&w3_e1_o[1], geom->sze-1,  MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&w2_t_o[1],  geom->szcS-1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&w3_t_o[1],  geom->szc-1,  MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);

    /* Do checks */
    /* 2D edges */
    for (ee = 1; ee < geom->szeS; ee++) {
      if (geom->w2_e1[ee] != w2_e1_o[ee]) {
	printf("w2_e1 mismatch, ee=%d, w2_e1[ee]=%d, w2_e1_o[ee]=%d\n",
	       ee, geom->w2_e1[ee], w2_e1_o[ee]);
	ret = 1;
	break;
      }
    }
    if (!ret) printf("w2_e1 check\n");
    
    /* 3D edges */
    if (!ret)
      for (ee = 1; ee < geom->sze; ee++) {
	if (geom->w3_e1[ee] != w3_e1_o[ee]) {
	  printf("w3_e1 mismatch, ee=%d, w3_e1[ee]=%d, w3_e1_o[ee]=%d\n",
		 ee, geom->w3_e1[ee], w3_e1_o[ee]);
	  ret = 1;
	  break;
	}
      }
    if (!ret) printf("w3_e1 check\n");

    /* 2D faces */
    if (!ret)
      for (cc = 1; cc < geom->szcS; cc++) {
	if (geom->w2_t[cc] != w2_t_o[cc]) {
	  printf("w2_t mismatch, cc=%d, w2_t[cc]=%d, w2_t_o[cc]=%d\n",
		 cc, geom->w2_t[cc], w2_t_o[cc]);
	  ret = 1;
	  break;
	}
      }
    if (!ret) printf("w2_t check\n");

    /* 2D faces */
    if (!ret)
      for (cc = 1; cc < geom->szc; cc++) {
	if (geom->w3_t[cc] != w3_t_o[cc]) {
	  printf("w3_t mismatch, cc=%d, w3_t[cc]=%d, w3_t_o[cc]=%d\n",
		 cc, geom->w3_t[cc], w3_t_o[cc]);
	  ret = 1;
	  break;
	}
      }
    if (!ret) printf("w3_t check\n");
    
    i_free_1d(w2_e1_o);
    i_free_1d(w3_e1_o);
    i_free_1d(w2_t_o);
    i_free_1d(w3_t_o);
    
    /* Send acknowledgement */
    MPI_Send(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    if (ret) return(ret);
    
  } else {
    /* Send sizes */
    MPI_Send(&geom->szeS,  1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->sze,   1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->szcS,  1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->szc,   1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->v2_e1, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->v3_e1, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->v2_t,  1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->v3_t,  1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);

    /* Check acknowledgement */
    MPI_Recv(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (ret) return(ret);

    /* Send sparse arrays */
    MPI_Send(&geom->w2_e1[1], geom->szeS-1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->w3_e1[1], geom->sze-1,  MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->w2_t[1],  geom->szcS-1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    MPI_Send(&geom->w3_t[1],  geom->szc-1,  MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);

    /* Check acknowledgement */
    MPI_Recv(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    if (ret) return(ret);
  }
  
#endif

  return(ret);
}

/*
 * Rank 0 assumes host
 * Rank 1 assumes remote
 */
int mpi_check_multi_windows_velocity(geometry_t *geom, master_t *master,
				     geometry_t **window, window_t **windat,
				     int flag)
{
  int ret = 0;
  
#ifdef HAVE_MPI
  static int first = 1;
  int mpi_rank_other = (mpi_rank == 0 ? 1 : 0);
  FILE *fp;
  char fname[MAXSTRLEN];
  
  /* global */
  int v_size, *w_e1, vel_sz;
  double *vel;

  /* windows */
  int win_v_size, *win_w_e1;
  double *win_vel;
  
  /* Only works for exactly 2 processes */
  if (mpi_size != 2)
    return(0);

  /* Set up sizes and arrays for this process */
  if (flag == VEL2D) {
    v_size = geom->b2_e1;
    w_e1   = geom->w2_e1;
    vel    = master->u1av; // same as windat[1]->u1av
    vel_sz = geom->szeS;
  } else
    if (flag == VEL3D) {
      /* 3D */
      v_size = geom->b3_e1;
      w_e1   = geom->w3_e1;
      vel    = master->u1; // same as windat[1]->u1
      vel_sz = geom->sze;
    } else
      if (flag == FLUX) {
	/* 3D */
	v_size = geom->b3_e1;
	w_e1   = geom->w3_e1;
	vel    = master->u1flux3d;
	vel_sz = geom->sze;
      } else {
	printf("mpi_check unsupported flag\n");
	return(1);
      }
  
  /* Construct unique file name */
  sprintf(fname, "check_vel_%d.txt", mpi_rank);
  
  /* Open file the first time around */
  if (first)
    fp = fopen(fname, "w");
  else {
    fp = fopen(fname, "a");
    first = 0;
  }
  
  fprintf(fp,"MPI Multi window check invoded\n");
  fprintf(fp,"------------------------------\n");
  fprintf(fp,"nstep = %d, t = %.5f, mpi_rank = %d, nwindows = %d, prmname = %s\n", master->nstep, master->t, mpi_rank,	 master->nwindows, prmname);
  fprintf(fp,"\n");
  
  /*
   * Multiple windows receives from single window
   */
  if (master->nwindows > 1) {
    double *vel_o;   /* _o for other */
    int sz = 0, n, ee;
    
    fprintf(fp,"Remote receives and compares\n\n");

    vel_o = d_alloc_1d(vel_sz);
    MPI_Recv(&vel_o[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);

    /* Print header */
    fprintf(fp,"es\te\tk\tle\twn\tvel_w\t\tvel_m\t\tdiff\n");

    /* Check array */
    for (n = 1; n <= master->nwindows; n++) {
      if (flag == VEL2D) {
	win_v_size = window[n]->b2_e1;
	win_w_e1   = window[n]->w2_e1;
	win_vel    = windat[n]->u1av;
      } else
	if (flag == VEL3D) {
	  win_v_size = window[n]->b3_e1;
	  win_w_e1   = window[n]->w3_e1;
	  win_vel    = windat[n]->u1;
	} else {
	  // FLUX
	  win_v_size = window[n]->b3_e1;
	  win_w_e1   = window[n]->w3_e1;
	  win_vel    = windat[n]->u1flux3d;
	}
      for (ee = 1; ee <= win_v_size; ee++) {
	int le = win_w_e1[ee];
	int e  = window[n]->wse[le];
	/* Compare results within epsilon of double precision */
	if (fabs(win_vel[le] - vel_o[e]) > DBL_EPSILON) {
	//if (fabs(win_vel[le] - vel_o[e]) > 1e-13) {
	  int es = geom->m2de[e];
	  int k  = geom->e2k[e];
	  fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%.3e\t%.3e\t%.3e\n", es, e, k, le, n,
		  win_vel[le], vel_o[e], fabs(win_vel[le]-vel_o[e]));
	  ret = 1;
	}
	sz++;
      }
    }
    fprintf(fp,"------------------------------\n\n");
    
    /* Make sure we got everything */
    if (sz != v_size) {
      fprintf(fp,"Not all edges accounted for %d vs %d\n",sz, v_size);
      ret = 1;
    }

    d_free_1d(vel_o);

    /* Send flag for a clean exit */
    MPI_Send(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    
  } else {
    fprintf(fp,"Host sends\n\n");

    /* Send velocity array itself */
    MPI_Send(&vel[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD);

    /* Receive exit status from dst */
    MPI_Recv(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
  }

  fclose(fp);

#endif

  return(ret);
}

/*
 * Rank 0 assumes host
 * Rank 1 assumes remote
 */
int mpi_check_multi_windows_Vz(geometry_t *geom, master_t *master,
			       geometry_t **window, window_t **windat)
{
  int ret = 0;
  
#ifdef HAVE_MPI
  static int first = 1;
  int mpi_rank_other = (mpi_rank == 0 ? 1 : 0);
  FILE *fp;
  char fname[MAXSTRLEN];
  
  /* global */
  int v_size, *w_t, vel_sz;
  double *vel;

  /* windows */
  int win_v_size, *win_w_t;
  double *win_vel;
  
  /* Only works for exactly 2 processes */
  if (mpi_size != 2)
    return(0);

  /* Set up sizes and arrays for this process */
  /* 3D */
  v_size = geom->v3_t;
  w_t    = geom->w3_t;
  vel    = master->Vz; // same as windat[1]->Vz
  //vel    = master->tr_wc[master->sno]; // same as windat[1]->Vz
  vel_sz = geom->szc;

  /* Construct unique file name */
  sprintf(fname, "check_Vz_%d.txt", mpi_rank);
  
  /* Open file the first time around */
  if (first)
    fp = fopen(fname, "w");
  else {
    fp = fopen(fname, "a");
    first = 0;
  }
  
  fprintf(fp,"MPI Multi window check invoded\n");
  fprintf(fp,"------------------------------\n");
  fprintf(fp,"nstep = %d, t = %.5f, mpi_rank = %d, nwindows = %d, prmname = %s\n", master->nstep, master->t, mpi_rank,
	 master->nwindows, prmname);
  fprintf(fp,"\n");

  /*
   * Multiple windows receives from single window
   */
  if (master->nwindows > 1) {
    double *vel_o;   /* _o for other */
    int sz = 0, n, cc;
    
    fprintf(fp,"Remote receives and compares\n\n");

    vel_o = d_alloc_1d(vel_sz);
    MPI_Recv(&vel_o[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    
    /* Print header */
    fprintf(fp,"cs\tc\tk\tlc\twn\tVz_w\t\tVz_m\n");
    
    /* Check array */
    for (n = 1; n <= master->nwindows; n++) {
      win_v_size = window[n]->v3_t;
      win_w_t    = window[n]->w3_t;
      win_vel    = windat[n]->Vz;
      //win_vel    = windat[n]->tr_wc[windat[n]->sno];
      
      for (cc = 1; cc <= win_v_size; cc++) {
	int lc = win_w_t[cc];
	int c  = window[n]->wsa[lc];
	/* Compare results within epsilon of double precision */
	if (fabs(win_vel[lc] - vel_o[c]) > DBL_EPSILON) {
	//if (fabs(win_vel[lc] - vel_o[c]) > 1e-8) {
	  int cs = geom->m2d[c];
	  int k  = geom->s2k[c];
	  fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%.3e\t%.3e\t%.3e\t%f\t%f\n", cs, c, k, lc, n,
		  win_vel[lc], vel_o[c], fabs(win_vel[lc]-vel_o[c]), geom->cellx[cs], geom->celly[cs]);
	  ret = 1;
	}
	sz++;
      }
    }
    fprintf(fp,"------------------------------\n\n");
    
    /* Make sure we got everything */
    if (sz != v_size) {
      fprintf(fp,"Not all edges accounted for %d vs %d\n",sz, v_size);
      ret = 1;
    }
    
    d_free_1d(vel_o);
    
    /* Send flag for a clean exit */
    MPI_Send(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    
  } else {
    fprintf(fp,"Host sends\n\n");
    
    /* Send velocity array itself */
    MPI_Send(&vel[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD);

    /* Receive exit status from dst */
    MPI_Recv(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
  }
  
  fclose(fp);
  
#endif
  
  return(ret);
}

/*
 * Rank 0 assumes host
 * Rank 1 assumes remote
 */
int mpi_check_multi_windows_tracer(geometry_t *geom, master_t *master,
				   geometry_t **window, window_t **windat)
{
  int ret = 0;
  
#ifdef HAVE_MPI
  static int first = 1;
  int mpi_rank_other = (mpi_rank == 0 ? 1 : 0);
  FILE *fp;
  char fname[MAXSTRLEN];
  
  /* global */
  int v_size, *w_t, vel_sz;
  double *vel;

  /* windows */
  int win_v_size, *win_w_t;
  double *win_vel;
  
  /* Only works for exactly 2 processes */
  if (mpi_size != 2)
    return(0);

  /* Set up sizes and arrays for this process */
  /* 3D */
  v_size = geom->v3_t;
  w_t    = geom->w3_t;
  vel    = master->tr_wc[master->sno]; // same as windat[1]->Vz
  vel_sz = geom->szc;

  /* Construct unique file name */
  sprintf(fname, "check_tr_%d.txt", mpi_rank);
  
  /* Open file the first time around */
  if (first)
    fp = fopen(fname, "w");
  else {
    fp = fopen(fname, "a");
    first = 0;
  }
  
  fprintf(fp,"MPI Multi window check invoded\n");
  fprintf(fp,"------------------------------\n");
  fprintf(fp,"nstep = %d, t = %.5f, mpi_rank = %d, nwindows = %d, prmname = %s\n", master->nstep, master->t, mpi_rank,
	 master->nwindows, prmname);
  fprintf(fp,"\n");

  /*
   * Multiple windows receives from single window
   */
  if (master->nwindows > 1) {
    double *vel_o;   /* _o for other */
    int sz = 0, n, cc;
    
    fprintf(fp,"Remote receives and compares\n\n");

    vel_o = d_alloc_1d(vel_sz);
    MPI_Recv(&vel_o[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
    
    /* Print header */
    fprintf(fp,"cs\tc\tk\tlc\twn\ttr_w\t\ttr_m\n");
    
    /* Check array */
    for (n = 1; n <= master->nwindows; n++) {
      win_v_size = window[n]->v3_t;
      win_w_t    = window[n]->w3_t;
      win_vel    = windat[n]->tr_wc[windat[n]->sno];
      
      for (cc = 1; cc <= win_v_size; cc++) {
	int lc = win_w_t[cc];
	int c  = window[n]->wsa[lc];
	/* Compare results within epsilon of double precision */
	if (fabs(win_vel[lc] - vel_o[c]) > DBL_EPSILON) {
	//if (fabs(win_vel[lc] - vel_o[c]) > 1e-8) {
	  int cs = geom->m2d[c];
	  int k  = geom->s2k[c];
	  fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%.3e\t%.3e\t%.3e\t%f\t%f\n", cs, c, k, lc, n,
		  win_vel[lc], vel_o[c], fabs(win_vel[lc]-vel_o[c]), geom->cellx[cs], geom->celly[cs]);
	  ret = 1;
	}
	sz++;
      }
    }
    fprintf(fp,"------------------------------\n\n");
    
    /* Make sure we got everything */
    if (sz != v_size) {
      fprintf(fp,"Not all edges accounted for %d vs %d\n",sz, v_size);
      ret = 1;
    }
    
    d_free_1d(vel_o);
    
    /* Send flag for a clean exit */
    MPI_Send(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD);
    
  } else {
    fprintf(fp,"Host sends\n\n");
    
    /* Send velocity array itself */
    MPI_Send(&vel[1], vel_sz-1, MPI_DOUBLE, mpi_rank_other, 0, MPI_COMM_WORLD);

    /* Receive exit status from dst */
    MPI_Recv(&ret, 1, MPI_INT, mpi_rank_other, 0, MPI_COMM_WORLD, NULL);
  }
  
  fclose(fp);
  
#endif
  
  return(ret);
}
