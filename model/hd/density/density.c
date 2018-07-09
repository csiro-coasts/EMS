/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/density/density.c
 *  
 *  Description:
 *  Calculate the model density field
 *  from salinity and temperature values
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: density.c 5841 2018-06-28 06:51:55Z riz008 $
 *
 */

#include "hd.h"


/* constants as in UNESCO report */
const double a0 = 999.842594;   /* change to -28.263737 for density
                                   anomality */
const double a1 = 6.793952e-2;
const double a2 = -9.095290e-3;
const double a3 = 1.001685e-4;
const double a4 = -1.120083e-6;
const double a5 = 6.536332e-9;

const double b0 = 8.24493e-1;
const double b1 = -4.0899e-3;
const double b2 = 7.6438e-5;
const double b3 = -8.2467e-7;
const double b4 = 5.3875e-9;

const double c0 = -5.72466e-3;
const double c1 = 1.0227e-4;
const double c2 = -1.6546e-6;

const double d0 = 4.8314e-4;

const double e0 = 19652.21;     /* change to -1930.06 for density
                                   anomality */
const double e1 = 148.4206;
const double e2 = -2.327105;
const double e3 = 1.360477e-2;
const double e4 = -5.155288e-5;

const double f0 = 54.6746;
const double f1 = -0.603459;
const double f2 = 1.09987e-2;
const double f3 = -6.1670e-5;

const double g0 = 7.944e-2;
const double g1 = 1.6483e-2;
const double g2 = -5.3009e-4;

const double h0 = 3.239908;     /* change to -0.1194975 for density
                                   anomality */
const double h1 = 1.43713e-3;
const double h2 = 1.16092e-4;
const double h3 = -5.77905e-7;

const double i0 = 2.2838e-3;
const double i1 = -1.0981e-5;
const double i2 = -1.6078e-6;

const double j_0 = 1.91075e-4;  /* renamed to avoid j0 bessel function */

const double k0 = 8.50935e-5;   /* change to 3.47718e-5 for density
                                   anomality */
const double k1 = -6.12293e-6;
const double k2 = 5.2787e-8;

const double m0 = -9.9348e-7;
const double m1 = 2.0816e-8;
const double m2 = 9.1697e-10;

INLINE double eos(double s, double t, double p)
{
  double d;                     /* final density */
  double dfw;                   /* density of pure water at p=0 */
  double dzero;                 /* density at p=0 */
  double k;
  double kzero;
  double kw;
  double a;
  double aw;
  double b;
  double bw;
  double s2;                    /* s squared */
  double s32;                   /* s to 3/2 */
  double t2;                    /* t squared */
  double t3;                    /* t cubed */
  double t4;                    /* t to the forth */
  double t5;                    /* t to the fifth */

  /* Sanity checks */
  if (s < 0.0)
    s = 0.0;
  if (t < 0.0)
    t = 0.0;                    /* Actually, very cold temperature should
                                   be handled in a compleletly different
                                   way! */

  s2 = s * s;
  s32 = sqrt(s2 * s);
  t2 = t * t;
  t3 = t2 * t;
  t4 = t3 * t;
  t5 = t4 * t;

  /* convert p to bars - note that the UNESCO code does this, but the text 
     does not make it clear that this is necessary (the UNESCO report
     talks about pressure in decibars, but their routine converts to
     bars). In our case, the input p is in N/m2 so we divide by 1e5. */
  p /= 1e5;

  dfw = a0 + a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5;
  dzero = dfw + s * (b0 + b1 * t + b2 * t2 + b3 * t3 + b4 * t4) +
    s32 * (c0 + c1 * t + c2 * t2) + d0 * s2;
  kw = e0 + e1 * t + e2 * t2 + e3 * t3 + e4 * t4;
  kzero = kw + s * (f0 + f1 * t + f2 * t2 + f3 * t3) +
    s32 * (g0 + g1 * t + g2 * t2);
  aw = h0 + h1 * t + h2 * t2 + h3 * t3;
  bw = k0 + k1 * t + k2 * t2;
  a = aw + s * (i0 + i1 * t + i2 * t2) + s32 * j_0;
  b = bw + s * (m0 + m1 * t + m2 * t2);
  k = kzero + p * (a + b * p);
  d = dzero / (1.0 - (p / k));
  return (d);
}

INLINE void eos2(double s, double t, double p, double *dens,
                 double *dens_0)
{
  double dfw;                   /* density of pure water at p=0 */
  double dzero;                 /* density at p=0 */
  double k;
  double kzero;
  double kw;
  double a;
  double aw;
  double b;
  double bw;
  double s2;                    /* s squared */
  double s32;                   /* s to 3/2 */
  double t2;                    /* t squared */
  double t3;                    /* t cubed */
  double t4;                    /* t to the forth */
  double t5;                    /* t to the fifth */

  /* Sanity checks */
  if (s < 0.0)
    s = 0.0;
  if (t < 0.0)
    t = 0.0;                    /* Actually, very cold temperature should
                                   be handled in a compleletly different
                                   way! */

  s2 = s * s;
  s32 = sqrt(s2 * s);
  t2 = t * t;
  t3 = t2 * t;
  t4 = t3 * t;
  t5 = t4 * t;

  /* convert p to bars - note that the UNESCO code does this, but the text 
     does not make it clear that this is necessary (the UNESCO report
     talks about pressure in decibars, but their routine converts to
     bars). In our case, the input p is in N/m2 so we divide by 1e5. */
  p /= 1e5;

  dfw = a0 + a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5;
  dzero = dfw + s * (b0 + b1 * t + b2 * t2 + b3 * t3 + b4 * t4) +
    s32 * (c0 + c1 * t + c2 * t2) + d0 * s2;
  kw = e0 + e1 * t + e2 * t2 + e3 * t3 + e4 * t4;
  kzero = kw + s * (f0 + f1 * t + f2 * t2 + f3 * t3) +
    s32 * (g0 + g1 * t + g2 * t2);
  aw = h0 + h1 * t + h2 * t2 + h3 * t3;
  bw = k0 + k1 * t + k2 * t2;
  a = aw + s * (i0 + i1 * t + i2 * t2) + s32 * j_0;
  b = bw + s * (m0 + m1 * t + m2 * t2);
  k = kzero + p * (a + b * p);
  *dens = dzero / (1.0 - (p / k));
  *dens_0 = dzero;
}

/**
  * Linear approximation to density equation of state. Calculated
  * from a regression on 9261 density values with s ranging from 0 to 40
  * (step 2.0), t ranging from 0 to 40 (step 2.0) and p ranging from 0 to
  * 4e6 (step 2e5)
  *
  * @param s Salinity in PSS 78  double
  * @param t Temperature in deg C double
  * @param p Pressure in N/m2 double
  * @return Density of seawater at s,t,p.
  */
INLINE double lindensity_w(double s, double t, double p)
{
  return (1002.134219 + 0.76195733 * s - 0.23438229 * t +
          0.00000044882688 * p);
}


/**
  * Quadratic approximation to density equation of state. Calculated
  * from a regression on 9261 density values with s ranging from 0 to 40
  * (step 2.0), t ranging from 0 to 40 (step 2.0) and p ranging from 0 to
  * 1000 (step 50.0).
  * 
  * @param s Salinity in PSS 78
  * @param t Temperature in deg C
  * @param p Pressure in N/m2
  * 
  * @return  Approximation of density of seawater for given s, t, p.
  */
double quaddensity_w(double s, double t, double p)
{
  double t2;
  double t3;

  t2 = t * t;
  t3 = t2 * t;
  return (1000.0441131667 + 0.799484556 * s + 0.01465623837 * t +
          0.00000049072513 * p + 6.6304848439e-5 * s * s -
          0.00530332744 * t2 - 0.00195294114 * s * t -
          3.17174874e-5 * s * t2 - 1.462783428e-9 * t * p -
          6.32128812e-10 * s * p + 1.040628257e-6 * s * t3);
}


/*-------------------------------------------------------------------*/
/* Routine to calculate the density on the slaves                    */
/*-------------------------------------------------------------------*/
void density_w(geometry_t *window,  /* Processing window */
               window_t *windat,  /* Window data structure */
               win_priv_t *wincon /* Window geometry / constants */
  )
{
  int cc;                       /* Counter */
  int c;                        /* Sparse coordinate */
  int cs;                       /* Sparse coordinate of the surface */
  int cb;                       /* Sparse coordinate of the bottom */
  double top, bot;              /* Height of the surface */
  double p;                     /* Pressure */
  double dz;
  int ndens = -wincon->calc_dens;

  /* Set pointers and initialise */
  /*
  memset(windat->dens, 0, window->sgsiz * sizeof(double));
  memset(windat->dens_0, 0, window->sgsiz * sizeof(double));
  */
  /* Overwrite density with the nominated tracer value */
  if (ndens > 0) {
    memcpy(windat->dens, windat->tr_wc[ndens], window->sgsiz * sizeof(double));
    memcpy(windat->dens_0, windat->tr_wc[ndens], window->sgsiz * sizeof(double));
    return;
  }

  /* Sanity checks */
  if (!windat->sal || !windat->temp)
    hd_quit("density_w: Tracers 'salt' and 'temp' must both be present\n");

  /* Note : cannot use the vector s1 and vcs here since newly wetted */
  /* cells are not included, hence the density will not be set in */
  /* cells where elevation has risen into the new (previously dry) */
  /* cell.  */
  for (cc = 1; cc <= window->b2_t; cc++) {
    c = window->w2_t[cc];
    cs = window->m2d[c];
    cb = window->bot_t[cc];

    /* Note : need to use undisturbed depth to calculate pressures */
    /* for density purposes, to avoid unwanted variations as eta */
    /* changes (which can lead to instability in the model).  */
    top = window->topgrid;
    p = windat->patm[cs];

    /* Loop through the water column to the bottom */
    while (c <= cb) {
      windat->dens[c] = 0.0;
      windat->dens_0[c] = 0.0;

      bot = (c == cb) ? window->botz[cs] : window->gridz[c];
      dz = top - bot;

      /* Use the pressure at the top face of the cell otherwise, if */
      /* the cell bottoms are at different depths, the density is */
      /* calculated at different effective depths, leading to */
      /* spurious horizontal gradients.  */
      eos2(windat->sal[c], windat->temp[c], p, &windat->dens[c],
           &windat->dens_0[c]);
      p += wincon->g * windat->dens[c] * dz;
      top = bot;
      c = window->zm1[c];
    }
    /* Set a no-gradient condition across the sediments */
    windat->dens[c] = windat->dens[window->zp1[c]];
  }
}

/* END density_w()                                                        */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to calculate the density on the master                    */
/*-------------------------------------------------------------------*/
void density_m(master_t *master /* Master data structure */
  )
{
  int cc;                       /* Counter */
  int c;                        /* Sparse coordinate */
  int cs;                       /* Sparse coordinate of the surface */
  int cb;                       /* Sparse coordinate of the bottom */
  double top;                   /* Height of the top of a layer */
  double bot;                   /* Height of the bottom of a layer */
  double p;                     /* Pressure */
  double dz;                    /* Cell thickness */
  int ndens = -master->calc_dens;
  geometry_t *geom = master->geom;

  /* Set pointers and initialise */
  memset(master->dens, 0, geom->sgsiz * sizeof(double));
  memset(master->dens_0, 0, geom->sgsiz * sizeof(double));

  /* Overwrite density with the nominated tracer value */
  if (ndens > 0) {
    memcpy(master->dens, master->tr_wc[ndens], geom->sgsiz * sizeof(double));
    memcpy(master->dens_0, master->tr_wc[ndens], geom->sgsiz * sizeof(double));
    return;
  }

  /* Sanity checks */
  if (!master->sal || !master->temp)
    hd_quit("density_m: Tracers 'salt' and 'temp' must both be present\n");

  /* Note : cannot use the vector s1 and vcs here since newly wetted */
  /* cells are not included, hence the density will not be set in */
  /* cells where elevation has risen into the new (previously dry) */
  /* cell.  */
  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];
    cs = geom->m2d[c];
    cb = geom->bot_t[cc];

    /* Note : need to use undisturbed depth to calculate pressures */
    /* for density purposes, to avoid unwanted variations as eta */
    /* changes (which can lead to instability in the model).  */
    bot = 0.0;
    p = 0.0;

    /* Loop through the water column to the bottom */
    while (c <= cb) {

      /* Use the pressure at the top face of the cell otherwise, if */
      /* the cell bottoms are at different depths, the density is */
      /* calculated at different effective depths, leading to */
      /* spurious horizontal gradients.  */
      eos2(master->sal[c], master->temp[c], p, &master->dens[c],
           &master->dens_0[c]);
      top = bot;
      bot = (c == cb) ? geom->botz[cs] : geom->gridz[c];
      dz = top - bot;
      p += master->g * master->dens[c] * dz;
      c = geom->zm1[c];
    }
    /* Set a no-gradient condition across the sediments */
    master->dens[c] = master->dens[geom->zp1[c]];
  }
}

/* END density_m()                                                      */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Get the density in a given wet water column                       */
/*-------------------------------------------------------------------*/
void density_c(geometry_t *window,  /* Processing window */
               window_t *windat,  /* Window data structure */
               win_priv_t *wincon, /* Window geometry / constants */
	       int cs
  )
{
  int cc;                       /* Counter */
  int c;                        /* Sparse coordinate */
  int zm1;                       /* Sparse coordinate of the bottom */
  double top, bot;              /* Height of the surface */
  double p;                     /* Pressure */
  double dz;

  /* Note : need to use undisturbed depth to calculate pressures */
  /* for density purposes, to avoid unwanted variations as eta */
  /* changes (which can lead to instability in the model).  */
  top = window->topgrid;
  p = windat->patm[cs];
  c = cs;
  zm1 = window->zm1[c];

  /* Loop through the water column to the bottom */
  while (c != zm1) {

    bot = (zm1 == window->zm1[zm1]) ? window->botz[cs] : window->gridz[c];
    dz = top - bot;

    /* Use the pressure at the top face of the cell otherwise, if */
    /* the cell bottoms are at different depths, the density is */
    /* calculated at different effective depths, leading to */
    /* spurious horizontal gradients.  */
    eos2(windat->sal[c], windat->temp[c], p, &windat->dens[c],
	 &windat->dens_0[c]);
    p += wincon->g * windat->dens[c] * dz;
    c = zm1;    top = bot;
    zm1 = window->zm1[c];
  }
  /* Set a no-gradient condition across the sediments */
  windat->dens[c] = windat->dens[window->zp1[c]];
}

/* END density_c()                                                   */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Get the density in a given ghost water column                     */
/*-------------------------------------------------------------------*/
void density_gc(geometry_t *window,  /* Processing window */
		window_t *windat,  /* Window data structure */
		win_priv_t *wincon, /* Window geometry / constants */
		int cs, 
		int *map
  )
{
  int cc;                       /* Counter */
  int c;                        /* Sparse coordinate */
  int zm1;                       /* Sparse coordinate of the bottom */
  double top, bot;              /* Height of the surface */
  double p;                     /* Pressure */
  double dz;

  /* Note : need to use undisturbed depth to calculate pressures */
  /* for density purposes, to avoid unwanted variations as eta */
  /* changes (which can lead to instability in the model).  */
  top = window->topgrid;
  p = windat->patm[cs];
  c = cs;
  zm1 = window->zm1[c];

  /* Loop through the water column to the bottom */
  while (c != zm1) {

    bot = (zm1 == window->zm1[zm1]) ? window->botz[cs] : window->gridz[c];
    dz = top - bot;

    /* Use the pressure at the top face of the cell otherwise, if */
    /* the cell bottoms are at different depths, the density is */
    /* calculated at different effective depths, leading to */
    /* spurious horizontal gradients.  */
    eos2(windat->sal[map[c]], windat->temp[map[c]], p, &windat->dens[map[c]],
	 &windat->dens_0[map[c]]);

    p += wincon->g * windat->dens[map[c]] * dz;
    c = zm1;    top = bot;
    zm1 = window->zm1[c];
  }
}

/* END density_gc()                                                  */
/*-------------------------------------------------------------------*/
