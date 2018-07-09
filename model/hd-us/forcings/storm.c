/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/forcings/storm.c
 *  
 *  Description:
 *  This file contains functions to specify moving
 *  pressure systems.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: storm.c 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#include <math.h>
#include <string.h>
#include "hd.h"

typedef struct {
  double stime;                 /* Time */
  int stype;                    /* High or low pressure system */
  double sp;                    /* Maximum pressure gradient in the high */
  double ss;                    /* Maximum radius of the high (grid
                                   spaces) */
  int si;                       /* i coordinate of the centre */
  int sj;                       /* i coordinate of the centre */
  double se;                    /* Eccentricity */
  double sr;                    /* Rotation */
} storm_t;

typedef struct {
  double dt;                    /* Time step */
  int nstorm;                   /* Number of storm tracks */
  storm_t **storm;              /* Storm track structures */
  double wind_scale;            /* Scale factor for wind speeds */
  double dlv0;                  /* Drag law v0 */
  double dlv1;                  /* Drag law v1 */
  double dlc0;                  /* Drag law c0 */
  double dlc1;                  /* Drag law c1 */
} wind_storm_t;

/* Prototypes */
void psystems(master_t *master, wind_storm_t *data, int io, int jo,
              double wmax, double dmax, double e, double ra, int mode,
              int dr);
storm_t *create_storm(void);

# define HIPR 0                 /* High pressure system */
# define LOPR 1                 /* Low pressure system */

/*-------------------------------------------------------------------*/
/* Global variables 						     */
int pgf = 6; /* pgf=0 : constant pressure gradient of pg */
             /* pgf=1 :linear, maximum at dmax/2, zero at d=0 and dmax */
           /* pgf=2 : parabolic : maximum at dmax/2, zero at d=0 and */
           /* dmax.  */
           /* pgf=3 : parabolic : maximum at dmax/2 to dmax and zero */
           /* at d=0.  */
           /* pgf=4 : parabolic : maximum at dmax/5,4dmax/5 and zero */
           /* at d=0, dmax.  */
           /* pgf=5 : parabolic : maximum at dmax/5 to dmax and zero */
           /* at d=0.  */
           /* pgf=6 : parabolic : maximum at pn to pt and zero at 0 */
           /* 0 and dmax.  */
double pnd = 300;  /* Inner maximum pressure gradient (km) for pgf = 6 */
double ptd = 1550; /* Outer maximum pressure gradient (km) for pgf = 6 */

/*-------------------------------------------------------------------*/
/* Routine to initialise the storm track                             */
/*-------------------------------------------------------------------*/
int storm_init(sched_event_t *event)
{
  int n;
  char keyword[MAXSTRLEN];
  char buf[MAXSTRLEN];

  FILE *fp = master->prmfd;
  master_t *master = (master_t *)schedGetPublicData(event);
  parameters_t *params = master->params;
  wind_storm_t *data = NULL;

  if (!master->nstorm)
    return 0;
  data = (wind_storm_t *)malloc(sizeof(wind_storm_t));
  data->nstorm = master->nstorm;

  schedSetPrivateData(event, data);
  data->dt = master->storm_dt;
  data->wind_scale = params->wind_scale;
  data->dlv0 = params->dlv0;
  data->dlv1 = params->dlv1;
  data->dlc0 = params->dlc0;
  data->dlc1 = params->dlc1;

  data->storm = (storm_t **)malloc(sizeof(storm_t *) * data->nstorm);
  prm_set_errfn(hd_quit);
  for (n = 0; n < data->nstorm; n++) {
    data->storm[n] = create_storm();
    sprintf(keyword, "ST%1d.stime", n);
    prm_read_char(fp, keyword, buf);
    tm_scale_to_secs(buf, &data->storm[n]->stime);
    sprintf(keyword, "ST%1d.stype", n);
    prm_read_char(fp, keyword, buf);
    if (strcmp(buf, "HIPR") == 0)
      data->storm[n]->stype = HIPR;
    if (strcmp(buf, "LOPR") == 0)
      data->storm[n]->stype = LOPR;
    sprintf(keyword, "ST%1d.sp", n);
    prm_read_double(fp, keyword, &data->storm[n]->sp);
    sprintf(keyword, "ST%1d.ss", n);
    prm_read_double(fp, keyword, &data->storm[n]->ss);
    /* Cannot use lat/long as input location since the position
       is often outside the grid and grid_xytoij() cannot find
       (i,j) locations of positions outside the grid.
    sprintf(keyword, "ST%1d.loc", n);
    prm_read_char(fp, keyword, buf);
    if (sscanf(buf, "%lf %lf",&x, &y) != 2)
      hd_quit_and_dump("storm_init: Can't read point list\n");
    grid_xytoij(master->xyij_tree, x, y, 
		&data->storm[n]->si, &data->storm[n]->sj);
    printf("%f %f %d %d\n",x,y,data->storm[n]->si, data->storm[n]->sj);
    */
    sprintf(keyword, "ST%1d.si", n);
    prm_read_int(fp, keyword, &data->storm[n]->si);
    sprintf(keyword, "ST%1d.sj", n);
    prm_read_int(fp, keyword, &data->storm[n]->sj);
    sprintf(keyword, "ST%1d.se", n);
    prm_read_double(fp, keyword, &data->storm[n]->se);
    sprintf(keyword, "ST%1d.sr", n);
    prm_read_double(fp, keyword, &data->storm[n]->sr);
  }
  return (1);
}

/* END storm_init()                                                  */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to move a storm systems along the route specified by the  */
/* data structure storms.                                            */
/*-------------------------------------------------------------------*/
double storm_event(sched_event_t *event, double time)
{
  master_t *master = (master_t *)schedGetPublicData(event);
  wind_storm_t *data = (wind_storm_t *)schedGetPrivateData(event);

  if (time >= (event->next_event - SEPS)) {
    int n1, n2, n = 0;
    double tm;
    int io, jo, type;
    double ecc, rot, pres, rad;
    double t1;
    double hmean = sqrt(master->hmean1 * master->hmean1 +
                        master->hmean2 * master->hmean2);
    storm_t **storm = data->storm;

    /* Get the data structures which bracket time */
    tm = storm[n]->stime;
    while (tm < time && n < data->nstorm - 1) {
      n++;
      tm = storm[n]->stime;
    }
    n2 = min(n, data->nstorm - 1);
    n1 = n - 1;

    if (n2 == 0) {
      /* Set the storm to the first enty */
      io = storm[0]->si;
      jo = storm[0]->sj;
      ecc = storm[0]->se;
      rot = storm[0]->sr;
      pres = storm[0]->sp;
      rad = storm[0]->ss * 1e3 / hmean;
      type = storm[0]->stype;
    } else if (n2 == data->nstorm) {
      /* Set the storm to the last entry */
      io = storm[data->nstorm]->si;
      jo = storm[data->nstorm]->sj;
      ecc = storm[data->nstorm]->se;
      rot = storm[data->nstorm]->sr;
      pres = storm[data->nstorm]->sp;
      rad = storm[data->nstorm]->ss * 1e3 / hmean;
      type = storm[data->nstorm]->stype;
    } else {
      /* Linearly interpolate */
      t1 =
        (time - storm[n1]->stime) / (storm[n2]->stime - storm[n1]->stime);
      io = (int)((storm[n2]->si - storm[n1]->si) * t1 + storm[n1]->si);
      jo = (int)((storm[n2]->sj - storm[n1]->sj) * t1 + storm[n1]->sj);
      ecc = (storm[n2]->se - storm[n1]->se) * t1 + storm[n1]->se;
      rot = (storm[n2]->sr - storm[n1]->sr) * t1 + storm[n1]->sr;
      pres = (storm[n2]->sp - storm[n1]->sp) * t1 + storm[n1]->sp;
      rad =
        ((storm[n2]->ss - storm[n1]->ss) * t1 +
         storm[n1]->ss) / (1e-3 * hmean);
      type = storm[n1]->stype;
    }
    if (master->wind_dt) {
      memcpy(master->wind1, master->swind1, geom->sgsizS * sizeof(double));
      memcpy(master->wind2, master->swind2, geom->sgsizS * sizeof(double));
      psystems(master, data, io, jo, pres, rad, ecc, rot, 1, type);
    } else
      psystems(master, data, io, jo, pres, rad, ecc, rot, 0, type);
    event->next_event += data->dt;
  }
  return event->next_event;
}

/* END storm_event()                                                 */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Clean up the storm events                                         */
/*-------------------------------------------------------------------*/
void storm_cleanup(sched_event_t *event, double t)
{
  wind_storm_t *data = (wind_storm_t *)schedGetPrivateData(event);
  int n;

  if (data != NULL) {
    for (n = 0; n < data->nstorm; n++)
      free(data->storm[n]);
    free(data);
  }
}

/* END storm_cleanup()                                               */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to create a wind stress field representing a high (anti-  */
/* cyclonic) or low (cyclonic) pressure systems.	 	     */
/*-------------------------------------------------------------------*/
void psystems(master_t *master, /* Master data structure */
              wind_storm_t *data, /* Storm data structure */
              int io,           /* x coordinate of centre of the high */
              int jo,           /* y coordinate of centre of the high */
              double wmax,      /* Maximum pressure gradient in the high */
              double dmax,      /* Maximum radius of the high (grid
                                   spaces) */
              double e,         /* Eccentricity (e=0 -> circular) */
              double ra,        /* Rotation angle (degrees) */
              int mode,         /* mode = 0 : reinitialize stress field */
              int dr            /* dr = 0 : high pressure system */
  )
                   /* dr = 1 : low pressure system */
{
  double d;                     /* Distance of (i,j) from centre of system 
                                 */
  double w;                     /* Wind stress at (i,j) */
  double a, ang;                /* Angle of stress vector to horizontal */
  double in;                    /* i-io */
  double jn;                    /* j-jo */
  double wu;                    /* x component of stress */
  double wv;                    /* y component of stress */
  double sf;                    /* Skew factor */
  double s;                     /* Skew factor */
  double ca;                    /* cos(ra) */
  double sa;                    /* sin(ra) */
  double ta;                    /* tan(ra) */
  double dc;                    /* 1-e*e */
  double x, y;                  /* (i,j) coordinates before rotation */
  double pi;                    /* Value of pi */
  double f;                     /* Coriolis parameter at point (i,j) */
  double dmin;                  /* Radius below which anticyclones
                                   unstable */
  double dn;                    /* Radius in m from centre of system */
  double pg = 0.0;              /* Pressure gradient */
  double rf;                    /* Rotation factor for boundary layer */
  double x1, y1, x2, y2;
  double hmean = sqrt(master->hmean1 * master->hmean1 +
                      master->hmean2 * master->hmean2);
  int pn = (int)pnd * 1e3 / hmean;
  int pt = (int)ptd * 1e3 / hmean;
  int c, cc, i, j, ee, e1, c1, c2;
  geometry_t *geom = master->geom;
  double *wind1, *wind2;

  /*-----------------------------------------------------------------*/
  /* Set parameters */
  wind1 = d_alloc_1d(geom->szc);
  wind2 = d_alloc_1d(geom->szc);
  pi = 3.141592654;
  sf = 1.0;

  if (ra >= 180.0)
    ra -= 180.0;
  if (ra < 0.0)
    ra += 180.0;
  ra *= pi / 180.0;
  ca = cos(ra);
  sa = sin(ra);
  ta = tan(ra);
  dc = 1.0 - e * e;

  for (cc = 1; cc <= geom->b2_t; cc++) {
    c = geom->w2_t[cc];

    /* Rotate the system */
    i = geom->s2i[c];
    j = geom->s2j[c];
    in = (double)(i - io);
    jn = (double)(j - jo);
    x = in * ca + jn * sa;
    y = jn * ca - in * sa;

    /* Get the distance of the major axis of the rotated system to */
    /* the origin (in grids).  */
    d = sqrt(x * x + y * y / dc);

    /* Get the wind speed via gradient wind equation */
    if (d <= dmax) {
      /* Pressure gradient */
      switch (pgf) {
      case 0:
        pg = wmax;
        break;
      case 1:
        if (d <= 0.5 * dmax)
          pg = 2.0 * wmax * d / dmax;
        else
          pg = -2.0 * wmax * (d - dmax) / dmax;
        break;
      case 2:
        pg = 4.0 * wmax * d * (dmax - d) / (dmax * dmax);
        break;
      case 3:
        if (d <= 0.5 * dmax)
          pg = 4.0 * wmax * d * (dmax - d) / (dmax * dmax);
        else
          pg = wmax;
        break;
      case 4:
        if (d < dmax / 4.0)
          pg = 16.0 * wmax * d * (dmax - d) / (3.0 * dmax * dmax);
        else if (d > 4.0 * dmax / 5.0)
          pg = 25.0 * wmax * d * (dmax - d) / (4.0 * dmax * dmax);
        else
          pg = wmax;
        break;
      case 5:
        if (d < dmax / 5.0)
          pg = 25.0 * wmax * d * (dmax - d) / (4.0 * dmax * dmax);
        else
          pg = wmax;
        break;
      case 6:
        if (d <= pn)
          pg = wmax * d * (2.0 * pn - d) / (pn * pn);
        else if (d >= pt)
          pg = wmax * d * (dmax - d) / (pt * (dmax - pt));
        else
          pg = wmax;
        break;
      }

      /* Coriolis */
      f = master->coriolis[c];
      f = fabs(f);

      /* Skew factor */
      if (sf != 1.0) {
        s = fabs(jn) * (sf - 1.0) / dmax + 1.0;
        if (jn > 0.0)
          pg *= s;
        if (jn < 0.0)
          pg /= s;
      }

      /* Distance of the rotated system to the origin (in metres) */
      dn = sqrt(x * x * geom->h1acell[c] * geom->h1acell[c] +
                y * y * geom->h2acell[c] * geom->h2acell[c] / dc);

      /* Mimimum distance for stability of system */
      dmin = 4.0 * 1.2 * pg / (f * f);
      wu = dn * f / 2.0;
      if (dn <= dmin)
        w = wu;
      else
        w = wu - sqrt(wu * wu - dn * 1.2 * pg);

      /* 
         w=10.0*d/dmax; if(jn>0.0)w*=s; if(jn<0.0)w/=s;

         w=10.0;

         if(d<=10.0)w=12.0*d*(20.0-d)/100.0; else w=12.0; */

      /* Rotate to allow for boundary layer rotation of wind */
      rf = pi * (10.0 + 1.4 * w) / 180.0;
      /* w=wspeed(w); */
      y1 = 0.0;
      windstress(&w, &y1, data->dlv0, data->dlv1, data->dlc0, data->dlc1);
      w *= data->wind_scale;

      /* Get the gradient of the rotated system */
      wu =
        dc * (in * ca * sa + jn * sa * sa) + jn * ca * ca - in * ca * sa;
      if (wu != 0) {
        a =
          (-dc * (in * ca * ca + jn * ca * sa) + jn * ca * sa -
           in * sa * sa) / wu;
        a = ang = fabs(atan(a));
      } else {
        ang = 0.0;
      }

      /* Find the coordinates (x1,y1) & (x2,y2) where slope = 0 and */
      /* slope is undefined.  */
      x1 = x2 = y1 = y2 = 0;
      if (e != 0.0) {           /* Non circular */
        if (ra != 0.0 && ra != pi / 2.0) {  /* ca != 0 and sa != 0 */
          a = (sa * sa + dc * ca * ca) / (ca * sa * (1.0 - dc));
          x1 =
            d * d * dc / (dc * (ca + sa * a) * (ca + sa * a) +
                          (ca * a - sa) * (ca * a - sa));
          x1 = sqrt(x1);
          y1 = a * x1;
          a = (ca * sa * (1.0 - dc)) / (dc * sa * sa + ca * ca);
          x2 =
            d * d * dc / (dc * (ca + sa * a) * (ca + sa * a) +
                          (ca * a - sa) * (ca * a - sa));
          x2 = sqrt(x2);
          y2 = a * x2;
        } else if (ra == pi / 2.0) {  /* ca = 0.0 */
          x1 = y2 = 0.0;
          y1 = d;
          x2 = sqrt(d * d * dc);
        } else if (ra == 0.0) { /* sa = 0.0 */
          x1 = y2 = 0.0;
          x2 = d;
          y1 = sqrt(d * d * dc);
        }
      } else {                  /* Circular : e = 0 */
        x1 = y2 = 0;
        y1 = x2 = d;
        if (ra > pi / 2.0 && ra < pi)
          y1 = -d;
      }

      /* Set the wind stress components wu & wv with correct */
      /* directions (a) : Rotation angle = 0 - 90 */
      wu = wv = 0.0;
      if (ra > 0.0 && ra < pi / 2.0) {
        if (in >= x1 && in < x2 && jn > y2 && jn <= y1) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          if (ang + rf <= pi / 2.0)
            wu *= -1.0;
        } else if (in > -x2 && in <= -x1 && jn >= -y1 && jn < -y2) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          wv *= -1.0;
          if (ang + rf > pi / 2.0)
            wu *= -1.0;
        } else if (in > -x2 && in < x1 && jn > -y2 && jn < y1 &&
                   jn > in * ta) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          wu *= -1.0;
          if (ang - rf >= 0.0)
            wv *= -1.0;
        } else if (in < x2 && in > -x1 && jn > -y1 && jn < y2) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          if (ang - rf < 0.0)
            wv *= -1.0;
        } else
          wu = wv = 0.0;
      }
      /* (b) : Rotation angle, ra = 90 - 180 */
      else if (ra > pi / 2.0 && ra < pi) {
        if (in < x2 && in >= -x1 && jn > y2 && jn <= -y1 && jn > in * ta) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          if (ang + rf <= pi / 2.0)
            wu *= -1.0;
        } else if (in <= x1 && in > -x2 && jn >= y1 && jn < -y2 &&
                   jn < in * ta) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          wv *= -1.0;
          if (ang + rf > pi / 2.0)
            wu *= -1.0;
        } else if (in < -x1 && in > -x2 && jn > -y2 && jn < -y1) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          wu *= -1.0;
          if (ang - rf >= 0.0)
            wv *= -1.0;
        } else if (in < x2 && in > x1 && jn > y1 && jn < y2) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          if (ang - rf < 0.0)
            wv *= -1.0;
        } else
          wu = wv = 0.0;
      }
      /* (c) : Rotation angle, ra = 90 or ra = 0 */
      else if (ra == pi / 2.0 || ra == 0.0) {
        if (in >= x1 && in < x2 && jn > y2 && jn <= y1) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          if (ang + rf <= pi / 2.0)
            wu *= -1.0;
        } else if (in > -x2 && in <= -x1 && jn >= -y1 && jn < -y2) {
          wu = w * fabs(cos(ang + rf));
          wv = w * sin(ang + rf);
          wv *= -1.0;
          if (ang + rf > pi / 2.0)
            wu *= -1.0;
        } else if (in > -x2 && in <= x1 && jn > -y2 && jn <= y1) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          wu *= -1.0;
          if (ang - rf >= 0.0)
            wv *= -1.0;
        } else if (in < x2 && in >= -x1 && jn >= -y1 && jn < y2) {
          wu = w * cos(ang - rf);
          wv = w * fabs(sin(ang - rf));
          if (ang - rf < 0.0)
            wv *= -1.0;
        } else
          wu = wv = 0.0;
      }
      if (jn == y2 || jn == -y2) {
        wu = w * sin(rf);
        wv = w * cos(rf);
        if (in < 0) {
          wu *= -1.0;
          wv *= -1.0;
        }
      }
      if (in == x1 || in == -x1) {
        wu = -w * cos(rf);
        wv = w * sin(rf);
        if (jn < 0) {
          wu *= -1.0;
          wv *= -1.0;
        }
        if (ra > pi / 2.0 && ra < pi && e != 0.0) {
          wu *= -1.0;
          wv *= -1.0;
        }
      }
      /* Cyclonic pressure systems */
      if (dr == 0) {
        wu *= -1.0;
        wv *= -1.0;
      }
    } else {
      wu = 0.0;
      wv = 0.0;
    }
    if (mode == 0) {
      wind1[c] = wu;
      wind2[c] = wv;
    } else {
      wind1[c] += wu;
      wind2[c] += wv;
    }
  }
  for (ee = 1; ee <= geom->b2_e1; ee++) {
    e1 = geom->w2_e1[ee];
    c1 = geom->e2c[e1][0];
    c2 = geom->e2c[e1][1];
    master->wind1[e1] = 0.5 * (wind1[c1] * geom->costhu1[e1] +
			      wind2[c1] * geom->sinthu1[e1] +
			      wind1[c2] * geom->costhu1[e1] +
			      wind2[c2] * geom->sinthu1[e1]);
  }
  d_free_1d(wind1);
  d_free_1d(wind2);
}

/* END psystems()                                                    */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Routine to allocate memory for the window data structure          */
/*-------------------------------------------------------------------*/
storm_t *create_storm(void)
{
  storm_t *storm = (storm_t *)malloc(sizeof(storm_t));
  memset(storm, 0, sizeof(storm_t));
  return storm;
}

/* END create_storm()                                                */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Returns a wind speed given a wind stress (divided by density of   */
/* water).                                                           */
/*-------------------------------------------------------------------*/
double wstress(double x         /* Wind speed */
  )
{
  double a;
  double cd = 1.4e-3;           /* Drag coefficient */
  double pa = 1.3;              /* Density of air */
  double pw = 1025.0;           /* Density of water */
  if (x != 0.0)
    a = x * sqrt(fabs(x * pw / (pa * cd))) / fabs(x);
  else
    a = 0.0;
  return (a);
}

/* END wstress()                                                     */
/*-------------------------------------------------------------------*/
