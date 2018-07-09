/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/master/sigma.c
 *  
 *  Description:
 *  Calculate k index for surface and bottom of
 *  water at each i,j position
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: sigma.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include "hd.h"

void check_sigma(master_t *master);
void smoothtopo(master_t *master, int n);

/*-------------------------------------------------------------------*/
/* Routine to initialize variables for sigma calculations            */
/*-------------------------------------------------------------------*/
void init_sigma(geometry_t *window, /* Window geometry               */
                window_t *windat,   /* Window data                   */
                win_priv_t *wincon  /* Window constants              */
  )
{
  int c, cp, cc, e, ee;             /* Sparse coordinates            */
  int n;
  open_bdrys_t **open = window->open;

  if (wincon->sigma) {
    for (cc = 1; cc <= window->a2_t; cc++) {
      c = window->w2_t[cc];
      wincon->Ds[c] = windat->eta[c] - wincon->Hs[c];
    }
    for (n = 0; n < window->nobc; n++) {
      for (ee = 1; ee <= open[n]->no2_e1; ee++) {
	c = open[n]->obc_e2[ee];
	e = open[n]->obc_e1[ee];
	cp = open[n]->omape[e][c];
	wincon->Ds[cp] = wincon->Ds[c];
      }
    }
    mdxs(window, wincon);
  } else
    memcpy(windat->topz, windat->eta, window->szcS * sizeof(double));
}

/* END init_sigma()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to set up the vertical sigma distribution. Use log        */
/* distributions at the bottom and top and linear spacing in the     */
/* interior.                                                         */
/*-------------------------------------------------------------------*/
void set_sigma_distrib(parameters_t *params,  /* Input parameter data
                                                 structure */
                       double bmin, /* Minimum bathymetry in grid */
                       double bmax  /* Maximum bathymetry in grid */
  )
{

  int k;
  int topl, botl, c1, c2;
  double topf;
  double cc, bb;
  double d1, d2;
  int nz;              /* Number of layers to use */
  double b = 0.0;      /* Exponential scaling factor */
  double dzlim = 0.0;  /* Maximum grid resolution */
  int explin = 0;      /* Type of profile */

  if(params->runmode & AUTO)
    explin = 1;

  /* Note : larger topf = more surface exponential layers */
  if (bmax == bmin) {
    nz = params->nz;
    topf = 0.3;
  } else if (bmax - bmin <= 50) {
    nz = params->nz = 24;
    topf = 0.175;
    b = 6;
    dzlim = 5;
  } else if (bmax - bmin > 50 && bmax - bmin <= 200) {
    nz = params->nz = 24;
    topf = 0.175;
    b = 8;
    dzlim = 50;
  } else if (bmax - bmin > 200 && bmax - bmin <= 1000) {
    nz = params->nz = 35;
    topf = 0.175;
    b = 9;
    dzlim = 100;
  } else if (bmax - bmin > 1000 && bmax - bmin <= 2500) {
    nz = params->nz = 35;
    topf = 0.2;
    b = 10;
    dzlim = 150;
  } else {
    if (params->runmode & ROAM) {
      nz = params->nz = 60;
      topf = 0.3;
      b = 15;
    } else {
      nz = params->nz = 41;
      topf = 0.25;
      b = 10;
    }
    dzlim = 200;
    explin = 1;
  }

  /* This vertical discretization assumes  an exponential profile at */
  /* at the surface downwards until the grid spacing becomes greater */
  /* than  dzlim.  A linear distribution  is then applied  until the */
  /* bottom, bmax, is  reached. Since it is  unknown how many lavels */
  /* are rquired to accomplish this, an initial guess is used to set */
  /* the  exponential  profile,  then  this  is  padded out with the */
  /* linear  distribution  until  bmax  is  reached, re-defining the */
  /* number of layers in the process. */
  if(explin) {
    double c=log(1.0/bmax)+(double)(nz)/b; /* Exponential adjustment */
    double s;            /* Linear scaling on exponential function */
    double d=1.0/exp(c); /* Slope of linear exponential scaling */
    double *nlayers;     /* Temporary layer array */
    double dz;           /* Vertical spacing */
    int n, onz=nz;       /* Inital guess at number of layers */

    nlayers = d_alloc_1d(params->nz + 1);
    /* Set up the initial exponential distribution */
    for (k = nz - 1; k >= 1; k--) {
      s = 1.0 - d * (double)(k - nz)/((double)nz - 1.0);
      nlayers[k] = 1.0/bmax - s * exp(c - (double)k / b);
    }

    /* Reset the surface layer to be 1.0 m */
    if(nlayers[nz-1] < -1.0 / bmax) {
      d1 = nlayers[nz-1] + 1.0 / bmax;
      for (k = nz - 1; k >= 1; k--)
	nlayers[k] -= d1;
    }

    /* Find the last layer with resolution less than dzlim */
    n = 0;
    for (k = nz - 1; k >= 0; k--) {
      n++;
      dz = bmax * (nlayers[k+1] - nlayers[k]);
      if(dz > dzlim  || nlayers[k] < -1) {
	break;
      }
    }

    /* Add layers until bmax is reached */
    if(k >= 0) {
      s = -nlayers[k] * bmax;
      while(s < bmax) {
	s += dzlim;
	n++;
      }
    }

    /* Define a new array, copy the surface exponential profile */
    /* from the old array, and pad out linearly until bmax is */
    /* reached. */
    nz = params->nz = n;
    params->layers = d_alloc_1d(params->nz + 1);
    n = nz - 1;
    for (k = onz - 1; k >= 0; k--) {
      params->layers[n] = nlayers[k];
      n--;
      if(nlayers[k+1] - nlayers[k] > dzlim / bmax)
	break;
      if (n < 0) 
	break;
    } 
    n++;
    for (; n >= 0; n--) {      
      params->layers[n] = params->layers[n+1] - dzlim / bmax;
    }
    d_free_1d(nlayers);
  }
  else {

    if (params->layers != NULL) d_free_1d(params->layers);
    params->layers = d_alloc_1d(params->nz + 1);

    topl = nz + 1 - (int)(topf * (nz + 1));
    botl = 2;

    c1 = (int)(topf * (nz + 1));
    c2 = nz - 1;

    bb = (double)(c2 - c1) + 4.0;
    cc = (double)(c1) - 2.0;

    d1 = 2.0 / bb / exp(0.693147 * (double)(c1 - 2));
    d2 = 2.0 / bb / exp(0.693147 * (double)(params->nz - c2));

    /* Surface exponential profile */
    for (k = nz - 1; k >= topl + 2; k--) {
      params->layers[k] = -d1 * exp(0.693147 * (double)(nz - k - 1));
    }
    /* Middle linear profile */
    for (k = topl + 1; k >= botl - 1; k--) {
      params->layers[k] = -((double)(nz - k + 1) - cc) / bb;
    }
    /* Bottom exponential profile */
    for (k = botl - 1; k >= 1; k--) {
      params->layers[k] = -(1.0 - d2 * exp(.693147 * (double)(k - 1)));
    }
  }
  params->layers[nz] = 0.0;
  params->layers[0] = -1.0;
}

/* END set_sigma_distrib()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to save the layer distribution to a diagnostic tracer.    */
/*-------------------------------------------------------------------*/
void save_sigma_layers(parameters_t *params, /* Input parameters     */
		       master_t *master      /* Master data          */
		       )
{
  geometry_t *geom = master->geom;
  int c, cs, cc, k;
  int nz = params->nz + 1;
  double *layers;

  layers = d_alloc_1d(nz);
  memcpy(layers, params->layers, nz * sizeof(double));
  set_sigma_distrib(params, 0, 0);

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    k = geom->s2k[c];
    if (master->layth)
      master->layth[c] = fabs((params->layers[k+1] - params->layers[k]) *
			      master->Ds[cs]);
  }
  
  memcpy(params->layers, layers, nz * sizeof(double));
  d_free_1d(layers);
}

/* END save_sigma_layers()                                           */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Tests the bottom gradients in the input topography file           */
/*-------------------------------------------------------------------*/
void testtopo(parameters_t *params, /* Input parameter data structure */
              master_t *master  /* Master data structure */
  )
{
  double maxgrad = 0.07;        /* Maximum allowable bottom gradient */
  int maxpass = 5;              /* Maximum number of smoothing passes */
  double mdh;                   /* Maximum bottom gradient found in file
                                   topo */
  double d1, d2;                /* Dummies */
  double dist = 0.0;            /* Grid spacing over maximum bottom
                                   gradient */
  int ic = 0, jc = 0;           /* Grid position at maximum bottom
                                   gradient */
  int nps = 0;                  /* Smoothing passes to reduce bottom
                                   gradient */
  int k;
  int c, cc;
  int xm1, ym1;

  /* Set maps self-mapping across tracer open boundaries */
  set_map_t(master->geom);

  for (k = 1; k <= maxpass; k++) {
    mdh = 0.0;
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      xm1 = geom->xm1[c];
      ym1 = geom->ym1[c];
      d1 = fabs(geom->botz[c] - geom->botz[xm1]) / geom->h1acell[c];
      if (d1 > mdh) {
        mdh = d1;
        ic = geom->s2i[c];
        jc = geom->s2j[c];
        dist = geom->h1acell[c];
      }
      d2 = fabs(geom->botz[c] - geom->botz[ym1]) / geom->h2acell[c];
      if (d2 > mdh) {
        mdh = d2;
        ic = geom->s2i[c];
        jc = geom->s2j[c];
        dist = geom->h2acell[c];
      }
    }
    if (mdh > maxgrad) {
      if (nps == 0) {
        printf
          ("\nWarning : Maximum bottom topography gradient too large:\n");
        printf("%4.0fm in %4.1fkm at point (%d,%d)\n", mdh * dist,
               dist / 1000.0, ic, jc);
      }
      nps += 1;
      smoothtopo(master, nps);
    }
  }
  if (nps != 0) {
    printf("%d smoothing passes performed.\n\n", nps);
  } else {
    smoothtopo(master, params->smooth);
  }
  /* 
     check_sigma(master); */
}

/* END testtopo()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Smooths the topography                                            */
/*-------------------------------------------------------------------*/
void smoothtopo(master_t *master, /* Model data structure */
                int n           /* Number of smoothing passes */
  )
{
  geometry_t *geom = master->geom;
  int c, cc, nn;
  double *aa;

  aa = d_alloc_1d(geom->sgsizS);
  for (nn = 1; nn <= n; nn++) {
    for (cc = 1; cc <= geom->b2_t; cc++) {
      c = geom->w2_t[cc];
      aa[c] = cvol1(master, geom->botz, c, 0);
    }
    memcpy(geom->botz, aa, geom->sgsizS * sizeof(double));
  }
  d_free_1d(aa);
}

/* END smoothtopo()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to apply a convolution smoothing filter to an element of  */
/* a 2-D array.                                                      */
/*-------------------------------------------------------------------*/
double cvol1(master_t *master,  /* Model data structure */
             double *a,         /* Array to smooth */
             int c,             /* Cell centre location */
	     int edge           /* Edge to sommoth */
  )
{
  double *kk;
  double fi, nf;
  geometry_t *geom = master->geom;
  int cc;
  int *st, sz = 3;
  int type = ST_SQ3;

  if( edge > 0) type |= ST_EDGE;

  st = stencil(geom, c, &sz, type, edge);
  /*
  kk = d_alloc_1d(sz);
  for (cc = 0; cc < sz; cc++) kk[cc] = 1.0 / (double)sz;
  */
  fi = nf = 0.0;
  for (cc = 0; cc < sz; cc++) {
    if (a[st[cc]] != 0.0) {
      fi += a[st[cc]];
      nf += 1.0;
    }
  }
  if (nf > 0.0) fi /= nf;
  i_free_1d(st);
  return (fi);
}

/* END cvol1()                                                       */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Checks that the grid satisfies hydrostatic consistency (see Haney */
/* 1991, JPO, 21, 610-619).                                          */
/*-------------------------------------------------------------------*/
void check_sigma(master_t *master /* Model data structure */
  )
{
  geometry_t *geom = master->geom;
  double d1, d2, d3 = 0.0;      /* Dummies */
  int cc, c, xm1, ym1, zp1;
  int cs;

  for (cc = 1; cc <= geom->b3_t; cc++) {
    c = geom->w3_t[cc];
    cs = geom->m2d[c];
    xm1 = geom->xm1[cs];
    ym1 = geom->ym1[cs];
    zp1 = geom->zp1[c];
    d1 = fabs(master->Ds[cs] - master->Ds[xm1]) / geom->h1acell[cs];
    d2 = fabs(geom->cellz[c] * d1 / master->Ds[cs]) * geom->h1acell[cs];
    d3 = geom->gridz[zp1] - geom->gridz[c];
    if (d2 > d3) {
      hd_warn
        ("Hydrostatic consistency violated at (x) (%d %d %d), %5.3f > %5.3f.\n",
         geom->s2i[c], geom->s2j[c], geom->s2k[c], d2, d3);
    }
    d1 = fabs(master->Ds[cs] - master->Ds[ym1]) / geom->h2acell[cs];
    d2 = fabs(geom->cellz[c] * d1 / master->Ds[cs]) * geom->h2acell[cs];
    if (d2 > d3) {
      hd_warn
        ("Hydrostatic consistency violated at (y) (%d %d %d), %5.3f > %5.3f.\n",
         geom->s2i[c], geom->s2j[c], geom->s2k[c], d2, d3);
    }
  }
}

/* END check_sigma()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Gets the mean density in each layer and stores as a reference     */
/* density which is subtracted from the actual density in the        */
/* pressure gradient calculation to reduce truncation error for the  */
/* sigma system.                                                     */
/*-------------------------------------------------------------------*/
void compute_ref_density(master_t *master,  /* Master data structure */
                         geometry_t *window,  /* Processing window */
                         window_t *windat,  /* Window data structure */
                         win_priv_t *wincon /* Window geometry / constants 
                                             */
  )
{
  int cc, c, k;
  double *rdens, *nc;
  geometry_t *geom = master->geom;

  memset(wincon->rdens, 0, window->sgsiz * sizeof(double));
  if (wincon->sigma) {
    rdens = d_alloc_1d(geom->nz * sizeof(double));
    memset(rdens, 0, geom->nz * sizeof(double));
    nc = d_alloc_1d(geom->nz * sizeof(double));
    memset(nc, 0, geom->nz * sizeof(double));
    for (cc = 1; cc <= window->b3_t; cc++) {
      c = window->w3_t[cc];
      k = window->s2k[c];
      rdens[k] += windat->dens[c];
      nc[k] += 1.0;
    }
    for (cc = 1; cc <= window->enon; cc++) {
      c = window->w3_t[cc];
      c = cc;
      k = window->s2k[c];
      if (nc[k])
        wincon->rdens[c] = rdens[k] / nc[k];
      wincon->rdens[c] = windat->dens[c];
    }
    /* No-gradient condition across the sediments */
    for (cc = 1; cc <= window->a2_t; cc++) {
      c = window->bot_t[cc];
      wincon->rdens[window->zm1[c]] = wincon->rdens[c];
    }
    d_free_1d(rdens);
    d_free_1d(nc);
  }
}

/* END compute_ref_density()                                         */
/*-------------------------------------------------------------------*/
