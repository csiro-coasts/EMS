/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/boundary/riverflow.c
 *  
 *  Description:
 *  u1/u2 boundary value routine for user defined
 *  boundary behaviour
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: riverflow.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include "hd.h"
#include "tracer.h"


static void u1_flow_do(geometry_t *window, window_t *windat,
                       flow_data_t *data, int c, int cb);
static void u2_flow_do(geometry_t *window, window_t *windat,
                       flow_data_t *data, int c, int cb);
static double flow_riverprofile(double z, double hc);

/*-------------------------------------------------------------------*/
/* Master initialisation routine for u1 river forcing                */
/*-------------------------------------------------------------------*/
void bf_u1_flow_init_m(master_t *master, open_bdrys_t *open,
                       bdry_details_t *d)
{
  flow_data_t *data = (flow_data_t *)malloc(sizeof(flow_data_t));
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];

  /* Initialise */
  memset(data, 0, sizeof(flow_data_t));
  d->custdata = (void *)data;
  data->init = 0;
  data->ocode = open->ocodex;
  data->options = open->options;
  data->ncells = (double)open->no2_t; 

  /* Get the flow profile parameters */
  strcpy(buf, open->bflow);
  data->tsflow = NULL;
  if (sscanf(buf, "%lf", &data->flow) != 1) {
    data->tsflow = hd_ts_read(master, buf, 1);
    data->flowid =
      ts_get_index(data->tsflow, fv_get_varname(buf, "flow", buf2));
    if (data->flowid < 0)
      hd_quit("%s does not contain the river flow variable.\n", buf);
  }
  data->hc = open->bhc;
  data->rlen = open->rlen;
  data->init = 1;
}

/* END bf_u1_flow_init_m()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window initialisation routine for u1 river forcing                */
/*-------------------------------------------------------------------*/
void bf_u1_flow_init_w(geometry_t *window,  /* Window geometry structure */
                       open_bdrys_t *open,  /* OBC data */
                       bdry_details_t *d, /* Custom boundary data for
                                             window */
                       bdry_details_t *din  /* Custom boundary data from
                                               master */
  )
{
  flow_data_t *data = (flow_data_t *)malloc(sizeof(flow_data_t));
  flow_data_t *data_in = din->custdata;

  /* Initialise */
  memset(data, 0, sizeof(flow_data_t));
  d->custdata = (void *)data;
  data->init = 0;
  data->ocode = open->ocodex;
  data->options = open->options;
  data->ncells = (double)open->no2_t; 
  data->hc = data_in->hc;
  /* The array data->v_river is assigned to one of the work arrays */
  /* in wincon (data->v_river = wincon->w1). This cannot be done */
  /* here since the wincon structure is not initialised yet. This */
  /* assignment is made in the first call to flowbdry_w() when */
  /* data->init=0 (subsequently set data->init=1).  */

}

/* END bf_u1_flow_init_w()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calling routine from the master for river forcing. This routine   */
/* reads the flow rate from file and stores it in the custom data    */
/* structure. The window and master share the same bdry_details_t    */
/* structure, hence all windows may also access this flow rate.      */
/*-------------------------------------------------------------------*/
void bf_flow_m(geometry_t *geom, master_t *master,
               open_bdrys_t *open, bdry_details_t *d)
{
  flow_data_t *data = (flow_data_t *)d->custdata;
  if (data->tsflow != NULL)
    data->flow = ts_eval(data->tsflow, data->flowid, master->t);
}

/* END bf_flow_m()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calling routine from the window for u1 river forcing. This        */
/* calculates the river profile and returns a velocity for the u1    */
/* boundary variable. This routine uses the flow rate read from file */
/* on the master and stored in data->flow.                           */
/*-------------------------------------------------------------------*/
double bf_u1_flow_w(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, open_bdrys_t *open, double t,
                    int c, int cc, bdry_details_t *d)
{
  flow_data_t *data = (flow_data_t *)d->custdata;

  /* Point the river array to the work array on the first call */
  if (!data->init) {
    data->v_river = wincon->w1;
    data->init = 1;
  }
  /* If the boundary cell lies on the surface, then calculate the */
  /* river profile for the whole water column.  */
  if (c == window->m2d[c]) {
    u1_flow_do(window, windat, data, c, open->bot_t[cc]);
    if (windat->riverflow)
      windat->riverflow[open->obc_t[cc]] = data->flow / (double)open->no2_t;
    if (windat->riverdepth) windat->riverdepth[open->obc_t[cc]] = data->hc;
  }
  return (data->v_river[c] / (double)open->no2_t);
}

/* END bf_u1_flow_w()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the u1 river flow forcing                              */
/*-------------------------------------------------------------------*/
static void u1_flow_do(geometry_t *window,  /* Window geometry */
                       window_t *windat,  /* Window data */
                       flow_data_t *data, /* Custom data structure */
                       int c,    /* Sparse boundary coordinate */
		       int cb    /* Bottom coordinate */
  )
{
  double sum;
  double z, dz;
  double trsp;
  double flow;
  int cs, zp1;

  /* Get the depth of the outflow layer using Keulegan (1966),     */
  /* Eqn. 11.16.                                                   */
  if (data->options & OP_DYNAHC && c == (cs = window->m2d[c])) {
    double eta = windat->eta[cs];
    double botz = window->botz[cs];
    double depth = eta - botz;
    /*double rho = max(window->wincon->densavu1[cs] / window->h1au1[cs], rhof);*/
    double rho = windat->dens_0[cb];
    double rhof = (data->options & OP_IFRESH) ? 1000.0 : windat->dens_0[cs];
    /* Internal wave speed, Eqn. 11.1                              */
    double vd = sqrt(window->wincon->g * depth * 2.0 * (rho-rhof)/(rho+rhof));
    double hc = (data->hc) ? data->hc : window->botz[cs]; 
    double h = (data->options & OP_FDEPTH) ? window->botz[cs] : hc;
    double vr = data->flow / (fabs(h) * window->h2acell[cs] * data->ncells);
    double k = 1.0 / pow(2.0, 2.0/3.0);
    double f = k * pow(2.0*vr/vd, 2.0/3.0);
    f = max(min(f, 1.0), 0.0);
    /* Reference to the free surface */
    data->hc = eta - depth * f;
    /* Find the first layer deeper than hc */
    if (data->options & OP_TRUNCL) {
      while (c != window->zm1[c]) {
	if (window->gridz[c] < data->hc) {
	  data->hc = max(window->gridz[c], window->botz[cs]);
	  break;
	}
	c = window->zm1[c];
      }
    }
    data->hc = max(data->hc, window->botz[cs]);
  }

  /* Calculate river flow profile */
  trsp = 0.0;
  c = cs = window->m2d[c];
  while (c != window->zm1[c]) {
    double eta = windat->eta[cs];
    double botz = window->botz[cs];
    double z1, z0 = max(window->gridz[c], botz);
    data->v_river[c] = 0.0;
    if (eta < z0) {
      c = window->zm1[c];
      continue;
    }
    zp1 = window->zp1[c];
    z1 = (c == zp1) ? eta : min(window->gridz[zp1], eta);
    dz = (z1 - z0) / 10.0;
    sum = 0.0;
    for (z = z0 + dz / 2.0; z < z1; z += dz)
      sum += flow_riverprofile(z-eta, max(data->hc-eta, botz-eta));
    data->v_river[c] = sum / 10.0;
    trsp += data->v_river[c] * (z1 - z0) * window->h2au1[cs];
    c = window->zm1[c];
  }

  flow = data->flow;
  /* Right edge river input => flow is negative into the domain */
  if (data->ocode & R_EDGE)
    flow *= (-1.0);
  /* Scale so transport equals specified river flow */
  c = cs;
  while (c != window->zm1[c] && trsp) {
    data->v_river[c] *= flow / trsp;
    c = window->zm1[c];
  }
}

/* END u1_flow_do()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to transfer custom data from the master to the slave      */
/*-------------------------------------------------------------------*/
void bf_flow_t(master_t *master,  /* Master data */
               open_bdrys_t *open_w,  /* OBC structure in the window */
               bdry_details_t *data,  /* Custom data for this OBC */
               geometry_t *window,  /* Window geometry */
               window_t *windat /* Window data */
  )
{
  geometry_t *geom = master->geom;  /* Master geometry structure */
  open_bdrys_t *open_m = geom->open[open_w->mwn]; /* Global OBC */
  bdry_details_t *d_m;          /* Master custom data structure */
  flow_data_t *f_m, *f_w;

  f_w = data->custdata;
  d_m = &open_m->datau1;
  if (open_m->datau2.explct)
    d_m = &open_m->datau2;
  f_m = d_m->custdata;
  f_w->flow = f_m->flow;
}

/* END bf_flow_t()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Cleanup for river custom routine                                  */
/*-------------------------------------------------------------------*/
void bf_flow_free(master_t *master,      /* Master data              */
		  bdry_details_t *data   /* Custom data              */
		  )
{
  flow_data_t *flow = data->custdata;

  if (flow->tsflow)
    hd_ts_free(master, flow->tsflow);
  free(flow);
}

/* END bf_flow_clean()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Initialisation routine for u2 river forcing                       */
/*-------------------------------------------------------------------*/
void bf_u2_flow_init_m(master_t *master, open_bdrys_t *open,
                       bdry_details_t *d)
{
  flow_data_t *data = (flow_data_t *)malloc(sizeof(flow_data_t));
  char buf[MAXSTRLEN];
  char buf2[MAXSTRLEN];

  /* Initialise */
  memset(data, 0, sizeof(flow_data_t));
  d->custdata = (void *)data;
  data->init = 0;
  data->ocode = open->ocodey;
  data->options = open->options;
  data->ncells = (double)open->no2_t; 

  /* Get the flow profile parameters */
  strcpy(buf, open->bflow);
  data->tsflow = NULL;
  if (sscanf(buf, "%lf", &data->flow) != 1) {
    data->tsflow = hd_ts_read(master, buf, 1);
    data->flowid =
      ts_get_index(data->tsflow, fv_get_varname(buf, "flow", buf2));
    if (data->flowid < 0)
      hd_quit("%s does not contain the river flow variable.\n", buf);
  }
  data->hc = open->bhc;
  data->rlen = open->rlen;
  data->init = 1;
}

/* END bf_u2_flow_init_m()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Window initialisation routine for u2 river forcing                */
/*-------------------------------------------------------------------*/
void bf_u2_flow_init_w(geometry_t *window,  /* Window geometry structure */
                       open_bdrys_t *open,  /* OBC data */
                       bdry_details_t *d, /* Custom boundary data for
                                             window */
                       bdry_details_t *din  /* Custom boundary data from
                                               master */
  )
{
  flow_data_t *data = (flow_data_t *)malloc(sizeof(flow_data_t));
  flow_data_t *data_in = din->custdata;

  /* Initialise */
  memset(data, 0, sizeof(flow_data_t));
  d->custdata = (void *)data;
  data->init = 0;
  data->ocode = open->ocodey;
  data->options = open->options;
  data->ncells = (double)open->no2_t; 
  data->hc = data_in->hc;
}

/* END bf_u2_flow_init_w()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calling routine from the window for u2 river forcing. This        */
/* calculates the river profile and returns a velocity for the u1    */
/* boundary variable. This routine uses the flow rate read from file */
/* on the master and stored in data->flow.                           */
/*-------------------------------------------------------------------*/
double bf_u2_flow_w(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, open_bdrys_t *open, double t,
                    int c, int cc, bdry_details_t *d)
{
  flow_data_t *data = (flow_data_t *)d->custdata;

  /* Point the river array to the work array on the first call */
  if (!data->init) {
    data->v_river = wincon->w1;
    data->init = 1;
  }

  /* If the boundary cell lies on the surface, then calculate the */
  /* river profile for the whole water column.  */
  if (c == window->m2d[c]) {
    u2_flow_do(window, windat, data, c, open->bot_t[cc]);
    if (windat->riverflow)
      windat->riverflow[open->obc_t[cc]] = data->flow / (double)open->no2_t;
    if (windat->riverdepth) windat->riverdepth[open->obc_t[cc]] = data->hc;
  }
  return (data->v_river[c] / (double)open->no2_t);
}

/* END bf_u2_flow_w()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Calculates the u2 river flow forcing                              */
/*-------------------------------------------------------------------*/
static void u2_flow_do(geometry_t *window,  /* Window geometry */
                       window_t *windat,  /* Window data */
                       flow_data_t *data, /* Custom data structure */
                       int c,    /* Sparse boundary coordinate */
		       int cb    /* Bottom coordinate */
  )
{
  double sum;
  double z, dz;
  double trsp;
  double flow;
  int cs, zp1;

  /* Get the depth of the outflow layer using Keulegan (1966),     */
  /* Eqn. 11.16.                                                   */
  if (data->options & OP_DYNAHC && c == (cs = window->m2d[c])) {
    double eta = windat->eta[cs];
    double botz = window->botz[cs];
    double depth = eta - botz;
    /*double rho = window->wincon->densavu2[cs] / window->h2au2[cs];*/
    double rho = windat->dens_0[cb];
    double rhof = (data->options & OP_IFRESH) ? 1000.0 : windat->dens_0[cs];
    /* Internal wave speed, Eqn. 11.1                              */
    double vd = sqrt(window->wincon->g * depth * 2.0 * (rho-rhof)/(rho+rhof));
    double hc = (data->hc) ? data->hc : window->botz[cs]; 
    double h = (data->options & OP_FDEPTH) ? window->botz[cs] : hc;
    double vr = data->flow / (fabs(h) * window->h1acell[cs] * data->ncells);
    /*double vr = windat->u2[cs];*/
    double k = 1.0 / pow(2.0, 2.0/3.0);
    double f = k * pow(2.0*vr/vd, 2.0/3.0);
    f = max(min(f, 1.0), 0.0);
    /* Reference to the free surface */
    data->hc = eta - depth * f;
    /* Find the first layer deeper than hc */
    if (data->options & OP_TRUNCL) {
      while (c != window->zm1[c]) {
	if (window->gridz[c] < data->hc) {
	  data->hc = max(window->gridz[c], window->botz[cs]);
	  break;
	}
	c = window->zm1[c];
      }
    }
    data->hc = max(data->hc, window->botz[cs]);
  }
  /* Calculate river flow profile */
  trsp = 0.0;
  c = cs = window->m2d[c];

  while (c != window->zm1[c]) {
    double eta = windat->eta[cs];
    double botz = window->botz[cs];
    double z1, z0 = max(window->gridz[c], botz);
    data->v_river[c] = 0.0;
    if (eta < z0) {
      c = window->zm1[c];
      continue;
    }
    zp1 = window->zp1[c];
    z1 = (c == zp1) ? eta : min(window->gridz[zp1], eta);
    dz = (z1 - z0) / 10.0;
    sum = 0.0;
    for (z = z0 + dz / 2.0; z < z1; z += dz)
      sum += flow_riverprofile(z-eta, max(data->hc-eta, botz-eta));
    data->v_river[c] = sum / 10.0;
    trsp += data->v_river[c] * (z1 - z0) * window->h1au2[cs];
    c = window->zm1[c];
  }
  flow = data->flow;
  /* Front edge river input => flow is negative into the domain */
  if (data->ocode & F_EDGE)
    flow *= (-1.0);
  /* Scale so transport equals specified river flow */
  c = cs;
  while (c != window->zm1[c] && trsp) {
    data->v_river[c] *= flow / trsp;
    c = window->zm1[c];
  }
}

/* END u2_flow_do()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* River profile - parabolic above halocline. Value returned here is */
/* always positive.                                                  */
/*-------------------------------------------------------------------*/
double flow_riverprofile(double z, double hc)
{
  return (z > hc) ? (z - hc) * (z + hc) : 0.0;
}

/* END flow_riverprofile()                                           */
/*-------------------------------------------------------------------*/
