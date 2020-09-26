/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/boundary/bdrycond.c
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
 *  $Id: bdrycond.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include "hd.h"

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
  int c, cc;                    /* Sparse coodinate / counter */
  int xm1;                      /* Sparse locations at i+1, i-1 */
  int yp1;                      /* Sparse locations at j+1, j-1 */
  int xmyp1;                    /* Sparse location at (i-1,y+1) */
  double pgt = 0.0;             /* Pressure gradient term */
  double cot = 0.0;             /* Coriolis term */
  double bft;                   /* Bottom friction term */
  double sft;                   /* Surface friction term */
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

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e1[cc];
    xm1 = window->xm1[c];
    yp1 = window->yp1[c];
    xmyp1 = window->xmyp1[c];
    midx = wincon->mdx[c];

    /*---------------------------------------------------------------*/
    /* Initialise                                                    */
    vel[c] = windat->u1avb[c];
    /* Multiply the velocity by the depth at the backward timestep   */
    /* for SIGMA calculations */
    if (wincon->sigma) windat->nu1av[c] *= mdxbs(window, windat, wincon, c);

    /*---------------------------------------------------------------*/
    /* Pressure gradient                                             */
    pgt = (dir == 0) ? 0 :
      -((windat->eta[c] - windat->eta[xm1]) * wincon->topdensu1[c] +
	(windat->patm[c] - windat->patm[xm1])) / wincon->densavu1[c];

    /*---------------------------------------------------------------*/
    /* Coriolis                                                      */
    u2au1 = 0.25 * (windat->u2av[c] + windat->u2av[xm1] +
		    windat->u2av[yp1] + windat->u2av[xmyp1]);
    cot = wincon->u1c5[c] * u2au1;

    /*---------------------------------------------------------------*/
    /* Surface friction                                              */
    sft = (ramp * windat->wind1[c] / wincon->topdensu1[c]);
				     
    /*---------------------------------------------------------------*/
    /* Bottom friction                                               */
    u2au1 = 0.25 * (windat->u2avb[c] + windat->u2avb[xm1] +
                    windat->u2avb[yp1] + windat->u2avb[xmyp1]);
    botu2 = 0.25 * (windat->u2bot[c] + windat->u2bot[xm1] +
                    windat->u2bot[yp1] + windat->u2bot[xmyp1]);
    botu1 = windat->u1avb[c] + windat->u1bot[c];
    botu2 = u2au1 + botu2;
    Cdu1 = 0.5 * (wincon->Cd[xm1] + wincon->Cd[c]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu1 * max(wincon->uf, val);
    /* Truncate to ensure stability */
    if (val > depth[c] * midx / windat->dt2d)
      val = depth[c] * midx / windat->dt2d;
    /* Note: depth=1 in the sigma case */
    bft = -val * botu1 / depth[c];

    /*---------------------------------------------------------------*/
    /* Calculate nu1av value                                         */
    /* SIGMA : multiply new velocity by depth (midx)                 */
    vel[c] += windat->dt2d * (midx * (pgt + cot) + sft + bft);
  }
}

/* END u1av_local()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Local solution for u2 depth avereged velocity.                    */
/* Follows Palma and Matano (1998), Section A4 using the model       */
/* discretization.                                                   */
/*-------------------------------------------------------------------*/
void u2av_local(geometry_t *window,   /* Window geometry             */
		window_t *windat,     /* Window data                 */
		win_priv_t *wincon,   /* Window constants            */
		open_bdrys_t *open,   /* Open boundary structure     */
		double *vel,          /* Output velocity             */
		int sb,               /* Start coordinate of vec     */
		int eb,               /* End coordinate of vec       */
		int dir               /* 0=normal, 1=tangential      */
		)
{
  int c, cc;                    /* Sparse coodinate / counter */
  int xp1;                      /* Sparse locations at i+1, i-1 */
  int ym1;                      /* Sparse locations at j+1, j-1 */
  int xpym1;                    /* Sparse location at (i-1,y+1) */
  double pgt = 0.0;             /* Pressure gradient term */
  double cot = 0.0;             /* Coriolis term */
  double bft;                   /* Bottom friction term */
  double sft;                   /* Surface friction term */
  double *depth;                /* Depth of the water column */
  double midy;                  /* Water depth at the cell face */
  double val;                   /* Dummy */
  double u1au2;                 /* u2 value at the e1 face */
  double botu1;                 /* Bottom velocity at e1 face */
  double botu2;                 /* Bottom velocity at e2 face */
  double Cdu2;                  /* Bottom drag coefficient at the e1 face */
  double ramp;                  /* Ramp value */

  depth = windat->depth_e2;
  ramp = (wincon->rampf & WIND) ? windat->rampval : 1.0;

  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_e2[cc];
    ym1 = window->ym1[c];
    xp1 = window->xp1[c];
    xpym1 = window->xpym1[c];
    midy = wincon->mdy[c];

    /*---------------------------------------------------------------*/
    /* Initialise                                                    */
    vel[c] = windat->u2avb[c];
    /* Multiply the velocity by the depth at the backward timestep   */
    /* for SIGMA calculations */
    if (wincon->sigma) windat->nu2av[c] *= mdybs(window, windat, wincon, c);

    /*---------------------------------------------------------------*/
    /* Pressure gradient                                             */
    pgt = (dir == 0) ? 0 :
      -((windat->eta[c] - windat->eta[ym1]) * wincon->topdensu2[c] +
	(windat->patm[c] - windat->patm[ym1])) / wincon->densavu2[c];

    /*---------------------------------------------------------------*/
    /* Coriolis                                                      */
    u1au2 = 0.25 * (windat->u1av[c] + windat->u1av[xp1] +
		    windat->u1av[ym1] + windat->u1av[xpym1]);
    cot = wincon->u2c5[c] * u1au2;

    /*---------------------------------------------------------------*/
    /* Surface friction                                              */
    sft = (ramp * windat->wind2[c] / wincon->topdensu2[c]);
				     
    /*---------------------------------------------------------------*/
    /* Bottom friction                                               */
    u1au2 = 0.25 * (windat->u1avb[c] + windat->u1avb[xp1] +
                    windat->u1avb[ym1] + windat->u1avb[xpym1]);
    botu1 = 0.25 * (windat->u1bot[c] + windat->u1bot[xp1] +
                    windat->u1bot[ym1] + windat->u1bot[xpym1]);
    botu2 = windat->u2avb[c] + windat->u2bot[c];
    botu1 = u1au2 + botu1;
    Cdu2 = 0.5 * (wincon->Cd[ym1] + wincon->Cd[c]);
    val = sqrt(botu1 * botu1 + botu2 * botu2);
    val = Cdu2 * max(wincon->uf, val);
    /* Truncate to ensure stability */
    if (val > depth[c] * midy / windat->dt2d)
      val = depth[c] * midy / windat->dt2d;
    /* Note: depth=1 in the sigma case */
    bft = -val * botu2 / depth[c];
    pgt=0.0;
    /*---------------------------------------------------------------*/
    /* Calculate nu1av value                                         */
    /* SIGMA : multiply new velocity by depth (midx)                 */
    vel[c] += windat->dt2d * (midy * (pgt + cot) + sft + bft);
  }
}

/* END u2av_local()                                                  */
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
  int c, c1, c2, cc;            /* Sparse coodinate / counter        */
  double colflux;               /* Velocity transport divergence     */
  double flux, fluxi;           /* Flux in normal direction          */
  double sgn;

  sgn = (open->ocodex & R_EDGE || open->ocodey & F_EDGE) ? -1.0 : 1.0;
  for (cc = sb; cc <= eb; cc++) {
    c = open->obc_t[cc];
    /* Get the fluxes at the normal boundary location, and interior  */
    /* location.                                                     */
    if (open->type == U1BDRY) {
      c1 = open->obc_e1[cc];
      c2 = open->oi1_e1[cc];
      flux = windat->u1av[c1] * windat->depth_e1[c1] * window->h2au1[c1] *
	wincon->mdx[c1] * windat->dt2d;
      fluxi = windat->u1av[c2] * windat->depth_e1[c2] * window->h2au1[c2] *
	wincon->mdx[c2] * windat->dt2d;
    }
    if (open->type == U2BDRY) {
      c1 = open->obc_e2[cc];
      c2 = open->oi1_e2[cc];
      flux = windat->u2av[c1] * windat->depth_e2[c1] * window->h1au2[c1] *
	wincon->mdy[c1] * windat->dt2d;
      fluxi = windat->u2av[c2] * windat->depth_e2[c2] * window->h1au2[c2] *
	wincon->mdy[c2] * windat->dt2d;
    }
    /* Update the elevation                                          */
    colflux = sgn * (fluxi - flux);
    val[c] = max(windat->etab[c] - colflux / window->cellarea[c],
			    window->botz[c]);
  }
}

/* END eta_local()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Averages variables on corners where 2 OBCs overlap                */
/*-------------------------------------------------------------------*/
void OBC_overlap(geometry_t *window, win_priv_t *wincon)
{
  int n, nn, c, cc;
  int *mask = wincon->i7;

  for (n = 0; n < window->nobc; n++) {
    open_bdrys_t *open = window->open[n];
    memset(mask, 0, window->sgsizS * sizeof(int));
    for (cc = 1; cc <= open->no2_t; cc++) {
      c = open->obc_t[cc];
      open->olap[cc] = -1;
      mask[c] = cc;
    }
    open->olap[0] = 0;

    for (nn = 0; nn < window->nobc; nn++) {
      open_bdrys_t *open2 = window->open[nn];
      if (nn == n) continue;
      for (cc = 1; cc <= open2->no2_t; cc++) {
	c = open2->obc_t[cc];
	if (mask[c]) {
	  open->olap[0] = 1;
	  open->olap[mask[c]] = nn;
	  if (DEBUG("init_m"))	
	    if (mask[c]) dlog("init_m", "OBC %s overlaps %s at (%d %d) : olap = %d\n", open->name, open2->name, 
			    window->s2i[c], window->s2j[c], open->olap[mask[c]]);
	  /*
	  printf("OBC %s overlaps %s at (%d %d) : olap = %d %d %d\n", open->name, open2->name, 
			    window->s2i[c], window->s2j[c], open->olap[mask[c]],nn,open2->id);
	  */
	}
      }
    }
  }
}

/* END OBC_overlap()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Averages variables on corners where 2 OBCs overlap                */
/*-------------------------------------------------------------------*/
void average_OBC_corner(geometry_t *window, open_bdrys_t *open, double *var, int mode)
{
  int n, c, cc;
  int c1, c2, c3;

  if (!open->olap[0]) return;

  for (cc = 1; cc <= open->no2_t; cc++) {
    c = open->obc_t[cc];
    n = open->olap[cc];
    if (n >= 0) {
      open_bdrys_t *open2 = window->open[n];
      c1 = open->nmap[c];
      c2 = open2->nmap[c];
      c3 = open->nmap[c2];
      var[c] = (var[c1] + var[c2] + var[c3]) / 3.0;
      if (mode) {
	while (c != window->zm1[c]) {
	  c = window->zm1[c];
	  c1 = window->zm1[c1];
	  c2 = window->zm1[c2];
	  c3 = window->zm1[c3];
	  var[c] = (var[c1] + var[c2] + var[c3]) / 3.0;
	}
      }
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
