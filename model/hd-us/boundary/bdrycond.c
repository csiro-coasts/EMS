/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/boundary/bdrycond.c
 *  
 *  Description:
 *  Miscellaneous routines for eta boundaries.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: bdrycond.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include "hd.h"

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Set the boundary conditions on the elevation or velocity.         */
/*-------------------------------------------------------------------*/
void set_OBC(geometry_t *window,   /* Processing window              */
             window_t *windat,     /* Window data structure          */
             win_priv_t *wincon,   /* Window geometry / constants    */
             open_bdrys_t *open,   /* OBC structure                  */
             int sb,               /* Start of OBC vector            */
             int eb,               /* End of OBC vector              */
             int *obc,             /* Open boundary vector           */
             int *oi1,             /* One cell interior to obc cells */
             int *oi2,             /* Two cell interior to obc cells */
             int *cyc,             /* Cyclic sparse coordinate       */
             double *vel,          /* Variable array at time t+1     */
             double *vel_t,        /* Variable array at time t       */
             double *vel_b,        /* Variable array at time t-1     */
             int bcond,            /* Open boundary condition        */
             double dt,            /* Time step                      */
             bdry_details_t *data, /* Custom data structure          */
	     double *transfer,     /* Transfer data                  */
	     int rlxn,             /* Flow relaxation zone           */
	     int code              /* General purpose code           */
  )
{
  int e, ee, j;                    /* Edge coorindate / counter      */
  int e1, e2, e3, e4;              /* Edge coordinates               */
  int cc, c, c1, c2;               /* Cell centered coordines        */
  double *newval;                  /* New velocity values            */
  double *fval;                    /* Forced velocity values         */
  int *imap = NULL;                /* Interior cell map              */
  double *cs;                      /* Phase speed                    */
  double sf;                       /* Phase speed smoothing factor   */
  double rts;                      /* Relaxation time scale          */
  int *m2d;                        /* 3D to 2D map                   */

  /* Return if no valid cells on this boundary                       */
  if (eb < sb) return;

  /* Initialise the pointers                                         */
  newval = wincon->w5;
  fval = wincon->w6;
  memset(fval, 0, window->szm * sizeof(double));
  m2d = (code == ETASPEC) ? window->m2d : window->m2de;

  /*-----------------------------------------------------------------*/
  /* Set the pointer to phase speed to allow phase speed smoothing   */
  /* for elevation if required.                                      */
  if (open->spf && vel_t == windat->eta) {
    cs = open->sphase;
    sf = open->spf;
  } else {
    cs = open->dum;
    memset(cs, 0, (open->no3_t + 1) * sizeof(double));
    sf = 0.0;
  }

  /*-----------------------------------------------------------------*/
  /* No action taken                                                 */
  /*if (bcond & (NOTHIN|LINEAR))*/
  if (bcond == NOTHIN || bcond == (NOTHIN|TIDALH) || bcond == LINEAR ||
      bcond == (NOTHIN|TIDALC))
    return;
  if (open->adjust_flux && code == ETASPEC) {
    if (bcond == (NOTHIN|FILEIN) || bcond == (NOTHIN|TIDALH|FILEIN) || 
	bcond == (NOTHIN|TIDALC|FILEIN) || bcond == (NOTHIN|TIDEBC|FILEIN))
      return;
  }

  /*-----------------------------------------------------------------*/
  /* Set a no gradient condition                                     */
  else if (bcond & NOGRAD) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      newval[e] = vel[e1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Cyclic                                                          */
  else if (bcond & (CYCLIC | CYCLED)) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = cyc[ee];      
      newval[e] = vel[e1];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Clamped                                                         */
  else if (bcond & CLAMPD) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      newval[e] = 0.0;
    }
  }

  /*-----------------------------------------------------------------*/
  /* Gravity wave radiation                                          */
  else if (bcond & GRAVTY) {
    /* Get the distance between cells                                  */
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e2 = m2d[e];
      if (code == ETASPEC) {
	e1 = window->m2de[open->oi1_e1[ee]];
	open->dum1[e2] = window->h2au1[e1];
      } else {
	e1 = (open->ocodec) ? oi1[ee] : e;
	open->dum1[e2] = window->h2au1[m2d[e1]];
      }
    }
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = m2d[e];
      newval[e] = bc_gravity(window, windat, wincon, open->dum1[e2], dt,
			     e2, vel_t[e], vel[e1], &cs[ee], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set an Orlanski radiation boundary condition                    */
  else if (bcond & ORLANS) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      newval[e] = bc_orlanski(vel_t[e], vel_t[e1], vel_t[e2], vel[e1],
			      vel_b[e], vel_b[e1], &cs[ee], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Camerlengo & O'Brien radiation boundary condition         */
  else if (bcond & CAMOBR) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      newval[e] = bc_camerlengo(vel_t[e], vel_t[e1], vel_t[e2], vel[e1],
				vel_b[e], vel_b[e1], &cs[ee], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Miller and Thorpe radiation boundary condition            */
  else if (bcond & MILLER) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      newval[e] = bc_miller(vel_t[e], vel_t[e1], vel_t[e2], vel[e1],
                            vel_b[e], vel_b[e1], vel_b[e2],
			    &cs[ee], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set a Raymond and Kuo radiation condition                       */
  else if (bcond & RAYMND) {
    int tp, tm, tp1, tm1;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      tp = open->tmpp[e];
      tm = open->tmpm[e];
      tp1 = open->tmpp[e1];
      tm1 = open->tmpm[e1];
      newval[e] = bc_raymond(vel_t[e], vel_t[e1], vel_t[e2], vel[e1],
                             vel[e2], vel_t[tp], vel_t[tm],
                             vel_t[tp1], vel_t[tm1], &cs[ee], sf);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Least squares linear interpolation                              */
  if (bcond & LINEXT) {
    imap = open->nmap;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      e3 = imap[e2];
      e4 = imap[e3];
      newval[e] = bc_leastsq(vel[e4], vel[e3], vel[e2], vel[e1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Polynomial extrapolation (2nd order)                            */
  if (bcond & POLEXT) {
    imap = open->nmap;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = oi2[ee];
      e3 = imap[e2];
      newval[e] = bc_polint(vel[e3], vel[e2], vel[e1]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Vertically  averaged velocity  from the vertically  integrated  */
  /* velocity. Note, velocities  above the  free surface  are set to */
  /* zero.                                                           */
  if (bcond & VERTIN) {
    double *sum = wincon->d1;
    int zm1;
    double *vel = NULL;
    double *dz = NULL;
    double *bottom = NULL;
    double *depth = NULL;
    double *md = NULL;
    int *m1 = NULL;

    /*UR-TODO  - nu1 and nu2 seem to be un-initialised under pthread*/
    /* Set the pointers and initialise                               */
    memset(sum, 0, window->szeS * sizeof(double));
    vel = windat->nu1; 
    dz = windat->dzu1;
    bottom = window->botzu1;
    depth = windat->depth_e1;
    md = wincon->mdx;
    m1 = window->em;

    /* Do the vertical integration                                   */
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = window->m2de[e];
      zm1 = window->zm1e[e];
      while(e != zm1) {
	sum[e1] += vel[e] * dz[e] * md[e1];
	e = zm1;
	zm1 = window->zm1e[e];
      }
    }

    /* Get the  vertically integrated  velocity. Note: this velocity */
    /* is invariant over the 2D time-step, hence the depth used must */
    /* be  the  3D  depth calculated  from  eta at the start  of the */
    /* time-step (i.e. oldeta  =  topz). The  2D  depth (depth_e1 or */
    /* depth_e2) must  also  be set to this depth so that this depth */
    /* is used in  eta_step()  to  calculate the flux, which is used */
    /* in velocity_adjust() to adjust the 3D velocities. */
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      c1 = window->m2d[window->e2c[e][0]];
      c2 = window->m2d[window->e2c[e][1]];
      depth[e] = (max(windat->topz[c1], windat->topz[c2]) - 
		  bottom[e]);/* * md[window->m2d[e]];*/
      newval[e] = sum[e] / max(depth[e] * md[window->m2de[e]], wincon->hmin);
    }
    /* Set the depth at OBC auxiliary cells */
    for (ee = 1; ee <= open->no2_a; ee++) {
      e = open->obc_a[ee];
      c1 = window->m2d[window->e2c[e][0]];
      c2 = window->m2d[window->e2c[e][1]];
      depth[e] = (max(windat->topz[c1], windat->topz[c2]) - 
		  bottom[e]);
    }
  }

  /*-----------------------------------------------------------------*/
  /* Local solution for tangential depth averaged velocity           */
  if (bcond & LOCALT) {
    u1av_local(window, windat, wincon, open, newval, sb, eb, 1);
  }

  /*-----------------------------------------------------------------*/
  /* Local solution for normal depth averaged velocity               */
  if (bcond & LOCALN) {
    u1av_local(window, windat, wincon, open, fval, sb, eb, 0);
  }

  /*-----------------------------------------------------------------*/
  /* Local solution for elevation                                    */
  if (bcond & LOCALE) {
    eta_local(window, windat, wincon, open, fval, sb, eb);
  }

  /*-----------------------------------------------------------------*/
  /* Tidal synthesis from harmonics                                  */
  if (bcond & (TIDALH|TIDALC)) {
    double gmt = windat->t - wincon->tz;
    double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
    if (code == ETASPEC) {
      for (ee = sb; ee <= eb; ee++) {
	e = obc[ee];
	/* Note; argument to csr_tide_eval() is GMT; correct for     */
	/* timezone.                                                 */
	fval[e] += (ramp * csr_tide_eval(&open->tc, ee, gmt));
      }
    } else if (code & U1BDRY) {
      double fvx, fvy;
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	c2 = window->m2d[c];
	fvx = (ramp * csr_tide_eval(&open->tun, cc, gmt));
	fvy = (ramp * csr_tide_eval(&open->tvn, cc, gmt));
	fval[c] += (fvx * window->costhu1[c2] + fvy * window->sinthu1[c2]);
      }
    } else {
      double fvx, fvy;
      for (cc = sb; cc <= eb; cc++) {
	c = obc[cc];
	c2 = window->m2d[c];
	fvx = (ramp * csr_tide_eval(&open->tut, cc, gmt));
	fvy = (ramp * csr_tide_eval(&open->tvt, cc, gmt));
	fval[c] += (fvx * window->costhu1[c2] + fvy * window->sinthu1[c2]);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* Custom specification. Note CUSTOM is handled by the master      */
  if (bcond & CUSTOM) {
    double ramp = (wincon->rampf & CUSTOM) ? windat->rampval : 1.0;

    /* Blend geostrophic velocities with forcing over the ramp if    */
    /* required (i.e. if VELOCITY GEOSTROPHIC is ued with CUSTOM u1  */
    /* or u2 forcing).                                               */
    if (wincon->vinit & VINIT_GEO && ramp < 1.0 && (code & (U1GEN|U2GEN))) {
      double *bvel = (code & U1GEN) ? open->u1d : open->u2d;
      double fvel;
      for (ee = sb; ee <= eb; ee++) {
	e = obc[ee];
	fvel = bdry_value_w(window, windat, wincon, open, data,
			    e, ee, windat->t);
	fval[e] = (1.0 - ramp) * bvel[ee] + ramp * fvel;
      }
    } else {
      for (ee = sb; ee <= eb; ee++) {
	e = obc[ee];
	fval[e] += ramp * bdry_value_w(window, windat, 
				       wincon, open, data,
				       e, ee, windat->t);
      }
    }
  }

  /*-----------------------------------------------------------------*/
  /* File input. Note FILEN is handled by the master, including ramp */
  if (bcond & FILEIN) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      fval[e] += transfer[ee];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Geostropic 3D current                                           */
  if (bcond & CUSTOM && open->options & (OP_GEOSTR|OP_YANKOVSKY)) {
    double *dz, *md, *cf, *hat, *depth;
    double d1, d2, d3, d4, fa, fm, fb, d, sgn;
    double *sum = wincon->d1;
    sgn = 1.0;

    dz = windat->dzu1;
    md = wincon->mdx;
    cf = wincon->u1c6;
    depth = windat->depth_e1;
    if (open->ocodec) sgn = -1.0;

    /* Integrate the flow over the face (barotropic component) */
    d1 = 0.0;
    d = 0.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      d1 += fval[e] * dz[e] * window->h1au1[window->m2de[e]];
      d += dz[e];
      vel[e] = fval[e];
      open->flow[ee] = fval[e];
    }
    /* Get the baroclinic component */
    memset(sum, 0, window->szeS * sizeof(double));
    d2 = d4 = 0.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e2 = window->m2de[e];
      j = open->ocodec;
      c = window->e2c[e][j];
      c2 = window->m2d[c];
      if (window->gridz[c1] > windat->eta[c2]) continue;
      c1 = window->e2c[e][0];
      c2 = window->e2c[e][1];
      sum[e2] += wincon->g * md[e2] * (windat->dens[c1] - windat->dens[c2]) * dz[c];
      d3 = windat->dt * cf[e2] * sum[e2] / (0.5 * (windat->dens[c1] + windat->dens[c2]));
      fval[e] += d3;
      newval[e] = d3;
      /* Yankovsky (2000) OBC */
      if (open->options & OP_YANKOVSKY)
	fval[e] = vel_t[oi1[ee]];
      d2 += fval[e] * dz[e] * window->h1au1[e2];
      d4 += d3 * dz[e] * window->h1au1[e2];
    }
    /* Adjust the profile so that the total flux = barotropic component */
    d3 = 0.0;
    e2 = window->m2de[obc[sb]];
    fa = (d1 - d2) / (depth[e2] * window->h1au1[e2]);
    fm = (d2) ? fabs(d1/d2) : 1.0;
    fb = (d1) ? (d1 - d4) / d1 : 0.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e2 = window->m2de[e];
      c = open->obc_t[ee];
      c2 = window->m2d[c];
      if (window->gridz[c] > windat->eta[e2]) continue;
      if (open->options & (OP_YANKOVSKY|OP_MULTF)) {
	/* Multiplicative scaling */
	fval[e] = (d2 * sgn > 0.0) ? fval[e] * fm : fval[e] + fa;
	/*fval[e] = (d1) ? fb * vel[e] + newval[e] : newval[e];*/
      } else {
	/* Additive scaling */
	fval[e] += fa;
      }
      d3 += fval[e] * dz[e] * window->h1au1[e2];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Flather radiation (elevation components).                       */
  /* NOTHIN for eta                                                  */
  if (bcond == (FLATHE|FILEIN|NOTHIN))
    return;
  /* FLATHR + Local Solution for eta forcing                         */
  if (bcond == (FLATHE|LOCALE)) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      vel[e] = fval[ee];
    }
    return;
  }
  /* Use FILEIN on eta                                               */
  if (bcond == (FLATHE|FILEIN)) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      vel[e] = transfer[ee];
    }
    return;
  }
  /* Use FILEIN + tide on eta                                        */
  if (bcond == (FLATHE|FILEIN|TIDALH)) {
    double gmt = windat->t - wincon->tz;
    double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      vel[e] = transfer[ee] + (ramp * csr_tide_eval(&open->tc, ee, gmt));
    }
    return;
  }
  /* The (FILEIN|CUSTOM) component is used for the 2D velocity OBC,  */
  /* the radiation component is used for the eta OBC.                */
  if (bcond & FLATHE) {
    if (!open->relax_time)open->relax_time = open->relax_timei = 1.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      fval[e] = vel_t[e];
    }
  }


  /*-----------------------------------------------------------------*/
  /* Tidal constituent specification                                 */
  if (bcond & TIDEBC) {
    double ramp = (wincon->rampf & TIDEBC) ? windat->rampval : 1.0;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      fval[e] += (ramp * bc_tidal(window, windat, wincon, open->tideforce, e));
    }
  }

  /*-----------------------------------------------------------------*/
  /* LINEAR computation or NOTHIN.                                   */
  /* NOTHIN may be used with a radiation condition (NOTHIN|GRAVTY),  */
  /* where the value on the boundary is partially passive via        */
  /* radiation, or with data (NOTHIN|FILEIN), where the value on the */
  /* boundary is implicitly relaxed to data.                         */
  if (bcond & (LINEAR|NOTHIN)) {
    if (bcond & (FILEIN|CUSTOM|TIDALH|TIDALC|TIDEBC)) {
      /* Relax to data implicitly                                    */
      for (ee = sb; ee <= eb; ee++) {
        e = obc[ee];
	rts = dt / open->relax_time;
	vel[e] = fval[e] = (vel_t[e] + rts * fval[e]) / (1.0 + rts);
      }
    } else {
      for (ee = sb; ee <= eb; ee++) {
	e = obc[ee];
	fval[e] += vel[e];
      }
    }
  }

  /*
    if (open->ocodex & (L_EDGE | R_EDGE))
      hat = window->h1acell;
    else if (open->ocodey & (B_EDGE | F_EDGE))
      hat = window->h2acell;
    if (open->ocodex & L_EDGE || open->ocodey & B_EDGE) {
      sc = -1.0;
      cim = obc;
    }
    if (bcond & (FILEIN|CUSTOM|LOCALN))
      uval = fval;
    else
      uval = newval;
    for (cc = sb; cc <= eb; cc++) {
      c = obc[cc];
      c1 = oi1[cc];
      c2 = window->m2d[c1];
      ci = window->m2d[cim[cc]];
      depth = windat->eta[c2] - window->botz[c2] * wincon->Ds[c2];
      spd = (depth > 0.0) ? sqrt(wincon->g / depth) : 0.0;
      spd = (open->spf) ? open->sphase[cc] * hat[c2] / (depth * dt) : spd;
      feta = 0.0;
      if (open->bcond_ele & (FILEIN|CUSTOM))
	feta += open->transfer_eta[cc];
      if (open->bcond_ele & (TIDALH|TIDALC)) {
	double gmt = windat->t - wincon->tz;
	double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
	feta += ramp * csr_tide_eval(&open->tc, cc, gmt);
      }
      if (!(open->bcond_ele & (FILEIN|CUSTOM|TIDALH|TIDALC))) {
	ce = window->m2d[open->obc_t[cc]];
	feta += windat->eta[ce];
      }
      fval[c] = uval[c] + (sc * spd * (windat->eta[ci] - feta));
  */

  /*-----------------------------------------------------------------*/
  /* Flather radiation (velocity components)                         */
  if (bcond & FLATHR) {
    double spd;
    double depth;
    double se = 1.0;
    double feta;
    double *hat = NULL;
    double *uval;
    int cc, ci, ce, j;

    hat = window->h1acell;
    if (bcond & (FILEIN|CUSTOM|LOCALN))
      uval = fval;
    else
      uval = newval;

    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      e2 = window->m2de[e1];

      j = open->ceni[ee];
      if (!j) se = -1.0;

      ci = window->m2d[open->obc_e2[ee]];
      c = window->e2c[e][j];
      c2 = window->m2d[c];

      depth = windat->eta[c2] - window->botz[c2] * wincon->Ds[c2];
      spd = (depth > 0.0) ? sqrt(wincon->g / depth) : 0.0;
      spd = (open->spf) ? open->sphase[ee] * window->h1au1[e2] / (depth * dt) : spd;
      /* Note eta is use at the interior cell (Palma & Matano (1998) */
      /* A.1 p 1340. These authors used GRAVTY on eta and tangential */
      /* component of velocity).                                     */
      feta = 0.0;
      if (open->bcond_ele & (FILEIN|CUSTOM)) {
	cc = open->e2c_e1[ee];
	feta += open->transfer_eta[ee];
      }
      if (open->bcond_ele & (TIDALH|TIDALC)) {
	double gmt = windat->t - wincon->tz;
	double ramp = (wincon->rampf & (TIDALH|TIDALC)) ? windat->rampval : 1.0;
	feta += ramp * csr_tide_eval(&open->tc, ee, gmt);
      }
      if (!(open->bcond_ele & (FILEIN|CUSTOM|TIDALH|TIDALC))) {
	j = (j) ? 0 : 1;
	c = window->m2d[window->e2c[e][j]];
	feta += windat->eta[c];
      }
      fval[e] = uval[e] + (se * spd * (windat->eta[ci] - feta));
    }
  }

  /*-----------------------------------------------------------------*/
  /* Set the boundary value to the calculated condition.             */
  if (bcond == FILEIN || bcond == TIDEBC || 
      bcond == TIDALH || bcond == TIDALC ||
      bcond == CUSTOM || bcond == (FILEIN|TIDALH) || 
      bcond == (FILEIN|TIDALC) || bcond & FLATHR) {
    /* Active conditions - set directly to supplied data.            */
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      vel[e] = fval[e];
    }
  } else if (bcond & (FILEIN | CUSTOM | TIDEBC | TIDALH | TIDALC | 
		      LINEAR | NOTHIN)) {
    /* Partially passive conditions - relax with specified radiation */
    /* over the specified relaxation time.                           */
    /* Note : LINEAR for velocity and NOTHIN for elevation assume a  */
    /* radiation condition applied to a local solution.              */
    if (bcond & (MILLER | NOGRAD)) {
      /* Explicit fromulation.                                       */
      for (ee = sb; ee <= eb; ee++) {
        e = obc[ee];
        vel[e] =
          newval[e] + dt * (fval[e] - vel_t[e]) / open->relax_time;
      }
    } else if (bcond & (ORLANS | CAMOBR)) {
      /* Implicit formulations.                                      */
      for (ee = sb; ee <= eb; ee++) {
        e = obc[ee];
	rts = (cs[ee] > 0.0) ? open->relax_time : open->relax_timei;
        vel[e] = newval[e] + 2.0 * dt * (fval[e] - vel_b[e]) /
          (rts * (1.0 + cs[ee]));
      }
    } else if (bcond & GRAVTY) {
      /* Implicit formulation. Wave speed is always positive here    */
      /* (outgoing).                                                 */
      for (ee = sb; ee <= eb; ee++) {
        e = obc[ee];
        vel[e] = newval[e] + dt * (fval[e] - vel_t[e]) /
          (open->relax_time * (1.0 + cs[ee]));
      }
    } else if (bcond & RAYMND) {
      /* Use the adaptive relaxation (Marchesiello et al, 2001)      */
      /* where if cs>0 (outward propagation) the relax time is long  */
      /* (about 1 year) and if cs<0 no radiation is applied with     */
      /* strong nudging (relax of a few days). Note that cx = cy = 0 */
      /* for cx < 0 (Eqn 15), hence the RAYMND scheme returns the    */
      /* value vel_t[e] when cx < 0.                                 */
      for (ee = sb; ee <= eb; ee++) {
        e = obc[ee];
	rts = (cs[ee] > 0.0) ? open->relax_time : open->relax_timei;
        vel[e] = newval[e] + dt * (fval[e] - vel_t[e]) /
	                          (rts * (1.0 + cs[ee]));
      }
    }
  } else {
    /* Passive conditions - set to passive OBC.                      */
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      vel[e] = newval[e];
    }
  }

  /*-----------------------------------------------------------------*/
  /* Impose a flow relaxtion condition if required                   */
  /* Relaxation to external solutions                                */
  if (rlxn) {
    double alpha;
    int bn;
    imap = open->nmap;
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      e1 = oi1[ee];
      bn = 2;
      while (e1 != imap[e1] && bn <= rlxn) {
        alpha = 1 - tanh(0.5 * (double)(bn - 1));
	sf = (code == 0) ? transfer[ee + (bn - 1) * eb] : vel[e];
	/*
	alpha = (double)(rlxn - bn + 1) / (double)rlxn;
	alpha *= alpha;
	vel[e1] = alpha * windat->rampval * vel[e] + (1.0 - alpha) * vel[e1];
	*/
	vel[e1] = alpha * windat->rampval * sf + (1.0 - alpha) * vel[e1];
        e1 = imap[e1];
        bn++;
      }
    }
  }


  /*-----------------------------------------------------------------*/
  /* Save the phase speed if required. Note: phase speed is bounded  */
  /* by the CFL condition (0 <= obc_phase <= 1). If waves are        */
  /* in-coming then obc_phase is negative, hence bounded to zero.    */
  /* Out-going waves have obc_phase > 0.                             */
  if (windat->obc_phase && vel_t == windat->eta && !open->adjust_flux) {
    for (ee = sb; ee <= eb; ee++) {
      e = obc[ee];
      windat->obc_phase[e] = cs[ee];
    }
  }
}

/* END set_OBC()                                                     */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Routine to get the values of a given array at previous time-steps */
/* on the boundaries and store in the old data structure.            */
/*-------------------------------------------------------------------*/
void store_old_bdry_values(double *F, /* df_variable_t array */
                           open_bdrys_t *open,  /* OBC structure */
                           bdry_old_values_t *old,  /* Old values data
                                                       structure */
                           int sb,  /* Start of OBC vector */
                           int eb,  /* End of OBC vector */
                           int *obc,  /* Open boundary vector */
                           int *oi1,  /* One cell interior to obc cells */
                           int *oi2 /* Two cell interior to obc cells */
  )
{
  int co, cc, c, c1, c2;        /* Sparse counters / coordinates */
  int tp, tm, tp1, tm1;
  for (cc = sb; cc <= eb; cc++) {
    co = cc - sb + 1;
    c = obc[cc];                /* Sparse coordinate on the boundary */
    c1 = oi1[cc];               /* One cell interior to c */
    c2 = oi2[cc];               /* Two cells interior to c */
    tp = open->tmpp[c];         /* Tangential coordinates (RAYMND) */
    tm = open->tmpm[c];
    tp1 = open->nmap[tp];
    tm1 = open->nmap[tm];
    old[co].t1i = old[co].ti;   /* F at time t-1 on the boundary */
    old[co].t1i1 = old[co].ti1; /* F at time t-1 one cell interior to c */
    old[co].t1i2 = old[co].ti2; /* F at t-1 two cells interior to c */
    old[co].ti = F[c];          /* F at time t on the boundary */
    old[co].ti1 = F[c1];        /* F at time t one cell interior to c */
    old[co].ti2 = F[c2];        /* F at time t two cells interior to c */
    old[co].tjp = F[tp];        /* F, time t tangential+1 on boundary */
    old[co].tjm = F[tm];        /* F, time t tangential-1 on boundary */
    old[co].tj1p = F[tp1];      /* F, time t tangential+1 on boundary-1 */
    old[co].tj1m = F[tm1];      /* F, time t tangential-1 on boundary-1 */
  }
}

/* END store_old_bdry_values()                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate an Orlanski radiation condition              */
double bc_orlanski(double F,    /* Boundary value at time t          */
                   double Fi1,  /* Value at boundary - 1 at time t   */
                   double Fi2,  /* Value at boundary - 2 at time t   */
                   double FFi1, /* Value at boundary - 1 at time t+1 */
                   double FB,   /* Boundary value at time t-1        */
                   double FBi1, /* Value at boundary - 1 at time t-1 */
                   double *ps,  /* Phase speed                       */
		   double sf    /* Phase speed smoothing factor      */
		   )
{
  double cs;

  /* Note : Chapman (1985) writes the denominator as                 */
  /* FFi1-FBi1-2.0*Fi2 while Palma & Matano (1998) write it as       */
  /* FFi1+FBi1-2.0*FBi2. Based on consistency with the explicit      */
  /* formulation it should be : FFi1+FBi1-2.0*Fi2.                   */
  cs = FFi1 + FBi1 - 2.0 * Fi2;
  cs = (cs) ? (FBi1 - FFi1) / cs : 1.0;
  if (cs >= 1.0) cs = 1.0;
  if (cs <= 0.0) cs = 0.0;

  cs = *ps = sf * (*ps) + (1.0 - sf) * cs;
  return ((FB * (1.0 - cs) + 2.0 * cs * Fi1) / (1.0 + cs));
}


/*-------------------------------------------------------------------*/
/* Routine to calculate a Camerlengo & O'Brien radiation condition   */
double bc_camerlengo(double F,    /* Boundary value at time t        */
		     double Fi1,  /* Value at boundary - 1 at time t */
		     double Fi2,  /* Value at boundary - 2 at time t */
		     double FFi1, /* Value at boundary-1 at time t+1 */
		     double FB,   /* Boundary value at time t-1      */
		     double FBi1, /* Value at boundary-1 at time t-1 */
                     double *ps,  /* Phase speed                     */
		     double sf    /* Phase speed smoothing factor    */
		     )

{
  double cs;

  cs = FFi1 + FBi1 - 2.0 * Fi2;
  cs = (cs) ? (FBi1 - FFi1) / cs : 1.0;
  if (cs >= 1.0) cs = 1.0;
  if (cs <= 0.0) cs = 0.0;
  cs = *ps = sf * (*ps) + (1.0 - sf) * cs;

  if (cs > 0.0)
    return (Fi1);
  else
    return (FB);
}

/*-------------------------------------------------------------------*/
/* Routine to calculate Miller and Thorpe radiation condition        */
double bc_miller(double F,    /* Boundary value at time t            */
		 double Fi1,  /* Value at boundary - 1 at time t     */
		 double Fi2,  /* Value at boundary - 2 at time t     */
		 double FFi1, /* Value at boundary - 1 at time t+1   */
		 double FB,   /* Boundary value at time t-1          */
		 double FBi1, /* Value at boundary - 1 at time t-1   */
                 double FBi2, /* Value at boundary - 2 at time t-1   */
		 double *ps,  /* Phase speed                         */
		 double sf    /* Phase speed smoothing factor        */
		 )
{
  double cs, c1, c2, c3;

  c1 = (Fi2 - Fi1);
  c1 = (c1) ? (FFi1 - Fi1) / c1 : 1.0;
  c2 = (FBi1 - FB);
  c2 = (c2) ? (F - FB) / c2 : 1.0;
  c3 = (FBi2 - FBi1);
  c3 = (c3) ? (Fi1 - FBi1) / c3 : 1.0;
  cs = c1 + c2 - c3;
  /* cs=c3 */

  if (cs >= 1.0)
   cs = 1.0;
  else if (cs < 0.0)
    cs = 0.0;
  cs = *ps = sf * (*ps) + (1.0 - sf) * cs;

  return (F - cs * (F - Fi1));
}

/*-------------------------------------------------------------------*/
/* Routine to calculate the Raymond and Kuo radiation condition.     */
/* Taken from Marchesiello et al (2001), Ocean Modelling, 3, 1-20.   */
double bc_raymond(double F,    /* Value on the boundary at time t    */
                  double Fi1,  /* Value at boundary - 1 at time t    */
                  double Fi2,  /* Value at boundary - 2 at time t    */
                  double FFi1, /* Value at boundary - 1 at time t+1  */
                  double FFi2, /* Value at boundary - 2 at time t+1  */
                  double Fjp,  /* Tan. boundary value at j+1, t      */
                  double Fjm,  /* Tan. boundary value at j-1, t      */
                  double Fj1p, /* Tan. value at boundary-1 at j+1, t */
                  double Fj1m, /* Tan. value at boundary-1 at j-1, t */
                  double *ps,  /* Phase speed                        */
		  double sf    /* Phase speed smoothing factor       */
		  )

{
  double rx = 1.0, ry = 1.0, dt, dx, dy, d, FF;

  dt = FFi1 - Fi1;
  dx = FFi1 - FFi2;

  dy = dt * (Fj1p - Fj1m);
  dy = (dy > 0.0) ? Fi1 - Fj1m : Fj1p - Fi1;

  d = dx * dx + dy * dy;
  if (d) {
    /* Note : this implementation is adaptive (Marchesiello et al    */
    /* (2001), Section 4.1) such that if inward propagation occurs   */
    /* (cx = dt/dx < 0) then rx = ry = 0 and FF = F. Applied with    */
    /* partially passive conditions (e.g. FILEIN|RAYMND) then F is   */
    /* nudged to data over timescale relax_time when cx < 0, and the */
    /* radiation condition is applied with weak nudging when cx > 0. */
    /* Marchesiello et al (2001) state radiation conditions alone    */
    /* without nudging to data is insufficient to maintain           */
    /* stability.                                                    */
    /* The pahse speeds are also bounded by the CFL condition here   */
    /* (i.e. rx < 1, ry < 1).                                        */
    rx = min(max(0.0, -dt * dx / d), 1.0);    /* rx=0 if cx<0, Eqn 15*/
    ry = (rx) ? min(-dt * dy / d, 1.0) : 0.0; /* ry=0 if cx<0, Eqn 15*/
  }

  d = (ry > 0) ? F - Fjm : Fjp - F;
  FF = (F + rx * FFi1 - ry * d) / (1 + rx);

  rx = *ps = sf * (*ps) + (1.0 - sf) * rx;
  return (FF);
}

/*-------------------------------------------------------------------*/
/* Routine to calculate gravity wave radiation condition             */
double bc_gravity(geometry_t *window, window_t *windat, win_priv_t *wincon,
                  double hat,   /* Grid spacing                      */
		  double dt,    /* Time step                         */
		  int c2,       /* 2D sparse coordinate              */
		  double F,     /* Value on the boundary at time t   */
		  double FF,    /* Value at boundary - 1 at time t+1 */
                  double *ps,   /* Phase speed                       */
		  double sf     /* Phase speed smoothing factor      */
		  )
{
  double cs = sqrt(wincon->g * (windat->eta[c2] - window->botz[c2] * 
			    wincon->Ds[c2])) * dt / hat;
  cs = *ps = sf * (*ps) + (1.0 - sf) * cs;
  return (F + cs * FF) / (1 + cs);
}

/*-------------------------------------------------------------------*/
/* Least squares linear interpolation from last 4 boundary values    */
/* y=a2+a1.x, x[1]=1, y[1]=boundary-3 : x[2]=2, y[2]=boundary-2 etc. */
double bc_leastsq(double y1, double y2, double y3, double y4)
{
  double sxy = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sx2 = 0.0;
  double x[5];
  double y[5];
  double a1, a2;
  int n = 4;
  int i;

  /* n=4 */

  y[1] = y1;
  x[1] = 1.0;
  y[2] = y2;
  x[2] = 2.0;
  y[3] = y3;
  x[3] = 3.0;
  y[4] = y4;
  x[4] = 4.0;

  /* n=3 */
  /* 
     y[1]=y2; x[1]=1.0; y[2]=y3; x[2]=2.0; y[3]=y4; x[3]=3.0; */
  for (i = 1; i <= n; i++) {
    sx += x[i];
    sy += y[i];
    sxy += (x[i] * y[i]);
    sx2 += (x[i] * x[i]);
  }
  a1 = (sxy - (sx * sy / n)) / (sx2 - (sx * sx) / n);
  a2 = sy / n - a1 * sx / n;
  return (a2 + a1 * 5.0);
}

/*-------------------------------------------------------------------*/
/* Routine to use polynomial extrapolation to estimate the value of  */
/* the normal velocity one grid point beyond the boundary.           */
/* From Numerical Recipes, p102.                                     */
double bc_polint(double y1, double y2, double y3)
{
  int n = 3;                    /* Order of the approximation + 1 */
  /* n= 2 : linear extrapolation */
  /* n= 3 : 2nd order extrapolation */
  double xa[5];                 /* Input array of x values, boundary value 
                                   = n */
  double ya[5];                 /* Input array of y values, n values in
                                   from bdry */
  double x = 0;                 /* x target value, boundary+1 */
  double y;                     /* y value at boundary+1 */
  double dy;                    /* Error estimate */
  int i, m, ns;                 /* Counters */
  double den, dif, dift, ho, hp, w, c[5], d[5];

  xa[1] = 1;
  xa[2] = 2;
  xa[3] = 3;
  if (n == 2) {
    x = 3;
    ya[1] = y2;
    ya[2] = y3;
    ya[3] = y3;
  } else if (n == 3) {
    x = 4;
    ya[1] = y1;
    ya[2] = y2;
    ya[3] = y3;
  }

  ns = 1;
  dif = fabs(x - xa[1]);
  for (i = 1; i <= n; i++) {
    dift = fabs(x - xa[i]);
    if (dift < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = d[i] = ya[i];
  }
  y = ya[ns];
  ns = ns - 1;
  for (m = 1; m <= n - 1; m++) {
    for (i = 1; i <= n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      den = ho - hp;
      if (den == 0.0)
        return (y3);
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if (2 * ns < n - m)
      dy = c[ns + 1];
    else {
      dy = d[ns];
      ns -= 1;
    }
    y += dy;
  }
  return (y);
}


/*-------------------------------------------------------------------*/
/* Routine to calculate surface elevation given tidal constituents   */
/*-------------------------------------------------------------------*/
double bc_tidal(geometry_t *window, /* Processing window */
                window_t *windat, /* Window data structure */
                win_priv_t *wincon, /* Window geometry / constants */
                tide_details_t *tide, /* Tidal forcing data structure */
                int c           /* Sparse location on the boundary */
  )
{
  int t;                        /* Tidal constituent counter */
  double xdist;                 /* x distance from tide mean cell to
                                   boundary cell */
  double ydist;                 /* y distance from tide mean cell to
                                   boundary cell */
  double d1, d2, d3;            /* Dummies */
  double newval = 0.0;          /* Tidal elevation */
  int cs = window->m2d[c];

  for (t = 0; t < tide[0].ntides; t++) {
    tide_details_t t1 = tide[t];

    xdist = window->cellx[cs] - t1.ic;
    ydist = window->celly[cs] - t1.jc;

    d1 = cos(2.0 * PI * windat->t / t1.per);
    d2 = sin(2.0 * PI * windat->t / t1.per);
    d3 = t1.amp * d1 - ydist * (t1.mta * sin(t1.dta) * d1 +
                                t1.amp * t1.mtp * sin(t1.dtp) * d2) -
      xdist * (t1.mta * sin(t1.dta) * d1 +
               t1.amp * t1.mtp * cos(t1.dtp) * d2);
    newval += d3;

  }
  return (newval);
}

/* END bc_tidal()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* No gradient condition for 2D velocity (equal fluxes)              */
/*-------------------------------------------------------------------*/
void bc_nograd_2d(double *vel, double *depth, double *md, double *h,
		  double hmin, int sb, int eb, int *obc, int *oi1)
{
  int c, c1, cc;
  double mid, mid1;
  for (cc = sb; cc <= eb; cc++) {
    c = obc[cc];
    c1 = oi1[cc];
    mid = md[c];
    mid1 = md[c1];
    /* Set the no gradient condition on velocity */
    vel[c] = vel[c1] * depth[c1] * mid1 * h[c1] /
      (max(depth[c], hmin) * mid * h[c]);
  }
}

/* END bc_nograd_2d()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Local solution for u1 depth averaged velocity.                    */
/* Follows Palma and Matano (1998), Section A4 using the model       */
/* discretization.                                                   */
/*-------------------------------------------------------------------*/
void u1av_local(geometry_t *window,   /* Window geometry             */
		window_t *windat,     /* Window data                 */
		win_priv_t *wincon,   /* Window constants            */
		open_bdrys_t *open,   /* Open boundary structure     */
		double *vel,          /* Output velocity             */
		int sb,               /* Start coordinate of vec     */
		int eb,               /* End coordinate of vec       */
		int dir               /* 0=normal, 1=tangential      */
		)
{
  int n, ee, e, es, eoe;
  int c1, c2;                   /* Sparse coodinate / counter */
  double pgt = 0.0;             /* Pressure gradient term */
  double cot = 0.0;             /* Coriolis term */
  double bft;                   /* Bottom friction term */
  double sft;                   /* Surface friction term */
  double *tzp;                  /* Surface height array */
  double *depth;                /* Depth of the water column */
  double midx;                  /* Water depth at the cell face */
  double val;                   /* Dummy */
  double u2au1;                 /* u2 value at the e1 face */
  double botu1;                 /* Bottom velocity at e1 face */
  double botu2;                 /* Bottom velocity at e2 face */
  double Cdu1;                  /* Bottom drag coefficient at the e1 face */
  double ramp;                  /* Ramp value */

  depth = windat->depth_e1;
  ramp = (wincon->rampf & WIND) ? windat->rampval : 1.0;

  for (ee = sb; ee <= eb; ee++) {
    e = open->obc_e1[ee];
    es = window->m2de[e];
    c1 = window->e2c[e][0];
    c2 = window->e2c[e][1];
    midx = wincon->mdx[e];

    /*---------------------------------------------------------------*/
    /* Initialise                                                    */
    vel[e] = windat->u1avb[e];
    /* Multiply the velocity by the depth at the backward timestep   */
    /* for SIGMA calculations */
    if (wincon->sigma) windat->nu1av[e] *= mdxbs(window, windat, wincon, e);

    /*---------------------------------------------------------------*/
    /* Pressure gradient                                             */
    pgt = (dir == 0) ? 0 :
      -((windat->eta[c1] - windat->eta[c2]) * wincon->topdensu1[e] +
	(windat->patm[c1] - windat->patm[c2])) / wincon->densavu1[e];

    /*---------------------------------------------------------------*/
    /* Coriolis                                                      */
    u2au1 = 0.0;
    for (n = 1; n <= window->nee[es]; n++) {
      eoe = window->eSe[n][e];
      u2au1 += window->wAe[n][e] * windat->u1av[eoe];
    }
    cot = wincon->u1c5[e] * u2au1;

    /*---------------------------------------------------------------*/
    /* Surface friction                                              */
    sft = (ramp * windat->wind1[e] / wincon->topdensu1[e]);
				     
    /*---------------------------------------------------------------*/
    /* Bottom friction                                               */
    u2au1 = 0.0;
    botu2 = 0.0;
    for (n = 1; n <= window->nee[es]; n++) {
      eoe = window->eSe[n][e];
      u2au1 += window->wAe[n][e] * windat->u1avb[eoe];
      botu2 += window->wAe[n][e] * windat->u1bot[eoe];
    }
    botu1 = windat->u1avb[e] + windat->u1bot[e];
    botu2 = u2au1 + botu2;
    Cdu1 = 0.5 * (wincon->Cd[c1] + wincon->Cd[c2]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu1 * max(wincon->uf, val);
    /* Truncate to ensure stability */
    if (val > depth[e] * midx / windat->dt2d)
      val = depth[e] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case */
    bft = -val * botu1 / depth[e];

    /*---------------------------------------------------------------*/
    /* Calculate nu1av value                                         */
    /* SIGMA : multiply new velocity by depth (midx)                 */
    vel[e] += windat->dt2d * (midx * (pgt + cot) + sft + bft);
  }
}

/* END u1av_local()                                                  */
/*-------------------------------------------------------------------*/



/*-------------------------------------------------------------------*/
/* Local solution for elevation (normal flux divergence only).       */
/* Follows Palma and Matano (1998), Section A4 using the model       */
/* discretization.                                                   */
/*-------------------------------------------------------------------*/
void eta_local(geometry_t *window,   /* Window geometry              */
	       window_t *windat,     /* Window data                  */
	       win_priv_t *wincon,   /* Window constants             */
	       open_bdrys_t *open,   /* Open boundary structure      */
	       double *val,          /* Output velocity              */
	       int sb,               /* Start coordinate of vec      */
	       int eb                /* End coordinate of vec        */
	       )
{
  int n, c, e, es, cc, cs;      /* Sparse coodinate / counter        */
  double flux, colflux;         /* Velocity transport divergence     */

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_t[cc];
    cs = window->m2d[c];
    /* Get the fluxes at the normal boundary location, and interior  */
    /* location.                                                     */
    /* Update the elevation                                          */
    colflux = 0.0;
    for (n = 1; n <= window->npe[cs]; n++) {
      e = window->c2e[n][c];
      flux = windat->u1av[e] * windat->depth_e1[e] * window->h1au1[e] *
	wincon->mdx[e] * windat->dt2d;
      colflux += window->eSc[n][c] * flux;
    }
    val[c] = max(windat->etab[c] - colflux / window->cellarea[c],
			    window->botz[c]);
  }
}

/* END eta_local()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Averages variables on corners where 2 OBCs overlap                */
/*-------------------------------------------------------------------*/
void average_OBC_corner(geometry_t *window, open_bdrys_t *open, double *var, int mode)
{
  int n, c, cs, cc;
  int ec, cn;
  double val, nv;

  ec = (mode) ?  open->no3_t : open->no2_t;
  for (cc = 1; cc <= ec; cc++) {
    c = open->obc_t[cc];
    if (!open->olap[cc]) {
      cs = window->m2d[c];
      val = 0.0;
      nv = 0.0;
      for (n = 1; n <= window->npe[cs]; n++) {
	if ((cn = open->bcc[n][cc]) < 0) {
	  val += var[abs(cn)];
	  nv += 1.0;
	}
      }
      if (nv > 1.0)
	var[c] = val / nv;
    }
  }
}

/* END average_OBC_corner()                                          */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the flux adjustment timescale, including constant,     */
/* default, adaptive linear / exponential or temporal. The default   */
/* value may be used in adaptive methods.                            */
/*-------------------------------------------------------------------*/
double calc_flux_adjust(geometry_t *window, 
			window_t *windat, 
			win_priv_t *wincon, 
			open_bdrys_t *open, 
			double eta, 
			int c, int ci)
{
  double ts;
  fa_info_t *fas = open->fas;
  int flag = 0;

  if (fas->tctype == RLX_CONS) {
    ts = (fas->tc0 > 0) ? fas->tc0 : window->h1acell[ci] / sqrt(wincon->g * windat->depth_e1[c]);
    return(ts);
  }

  /* Set the relaxation endpoints to the default timescale if < 0    */
  if (fas->tc0 <= 0) {
    fas->tc0 = window->h1acell[ci] / sqrt(wincon->g * windat->depth_e1[c]);
    flag = 1;
  }
  if (fas->tc1 <= 0) {
    fas->tc1 = window->h1acell[ci] / sqrt(wincon->g * windat->depth_e1[c]);
    flag = 1;
  }

  /* Compute the relaxation timescale                                */
  if (fas->tctype & RLX_LINR) {
    if (flag)
      fas->slope = (fas->dv1 - fas->dv0) ? (fas->tc1 - fas->tc0) / (fas->dv1 - fas->dv0) : 0.0;
    ts = (fabs(eta) - fas->dv0) * fas->slope + fas->tc0;
    ts = min (max(ts, fas->tc0), fas->tc1);
  }
  if (fas->tctype & RLX_EXP) {
    if (flag)
      fas->slope = fas->dv1 * log(fas->tc0 / 86400.0);
    ts = exp(fas->dv0 / fabs(eta));
  }
  if (fas->tctype & RLX_TIME) {
    if (flag)
      fas->slope = (fas->dv1 - fas->dv0) ? (fas->tc1 - fas->tc0) / (fas->dv1 - fas->dv0) : 0.0;
    ts = min((windat->days - fas->dv0) * fas->slope + fas->tc0, fas->tc1);
  }
  return(ts);
}

/* END calc_flux_adjust()                                            */
/*-------------------------------------------------------------------*/
