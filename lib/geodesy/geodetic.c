/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/geodesy/geodetic.c
 *
 *  \brief Spheroid calculations
 *
 *  Computes the geodetic distance between two points
 *  on a spheroid, or compute the lat/long given a point
 *  and distance.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: geodetic.c 5831 2018-06-26 23:48:06Z riz008 $
 */

#include <math.h>
#include <stdio.h>
#include "ems.h"

#define EABS 1e-4               /* Absolute error */
#define EPROP 1e-9              /* Propotional error */

/** Inverse Robbins computation for distance, and azimuth on spheroid.
  *
  * Checked against example in "The Austalian Map Grid"
  * (National Mapping Council of Australia, Special
  * Publication 7), page 10, 18 Feb. 1986 - accurate to
  * +/- .001 metre, +/- .01 sec.
  *
  * @param x1 Longitude at point 1.
  * @param y1 Latitude at point 1.
  * @param z1 Height at point 1.
  * @param x2 Longitude at point 2.
  * @param y2 Latitude at point 2.
  * @param z2 Height at point 2.
  * @param a Semi-major axis of Earth.
  * @param ecs First eccentricity squared of Earth.
  * @return Distance.
  */
double geod_inv_robbins(double x1, double y1, double z1,
                        double x2, double y2, double z2,
                        double a, double ecs)
{
  double eps = ecs / (1 - ecs);
  double sp1 = sin(y1);
  double sp2 = sin(y2);
  double f1 = a / sqrt(1 - ecs * sp1 * sp1);
  double f2 = a / sqrt(1 - ecs * sp2 * sp2);
  double cp1 = cos(y1);
  double cdl = cos(x2 - x1);
  double sdl = sin(x2 - x1);
  double t = ((1 - ecs) * sp2 + ecs * sp1 * f1 / f2) / cos(y2);
  double tau1 = cp1 * t - sp1 * cdl;
  double chi = sqrt(tau1 * tau1 + sdl * sdl);
  double sig = asin(chi / sqrt(1 + t * t));
  double caz = tau1 / chi;
  double dum = cp1 * caz;
  double h2 = dum * dum * eps;
  double dum1 = f1 * sig;
  double dums = dum1;
  double sigp = sig * sig;
  double dum3 = -sigp * h2 * (1 - h2) * 0.16666666666666667;
  double dum2 = dum1 * dum3;
  dums = dums + dum2;
  if ((fabs(dum2) > EABS) && (fabs(dum3) > EPROP)) {
    double gh;
    sigp = sigp * sig;
    gh = sp1 * dum * eps;
    dum3 = sigp * gh * (1 - 2 * h2) * 0.125;
    dum2 = dum1 * dum3;
    dums = dums + dum2;
    if ((fabs(dum2) > EABS) && (fabs(dum3) > EPROP)) {
      double g2;
      sigp = sigp * sig;
      g2 = sp1 * sp1 * eps;
      dum3 =
        sigp * (h2 * (4 - 7 * h2) -
                3 * g2 * (1 - 7 * h2)) * 8.3333333333333333e-3;
      dum2 = dum1 * dum3;
      dums = dums + dum2;
      if ((fabs(dum2) > EABS) && (fabs(dum3) > EPROP)) {
        sigp = sigp * sig;
        dum3 = -sigp * gh * 2.0833333333333333e-2;
        dum2 = dum1 * dum3;
        dums = dums + dum2;
      }
    }
  }

  /* Correction for heights above sea level ..... */
  dums = dums * (1 + (z1 + z2) / (f1 + f2));
  return sqrt(dums * dums + (z1 - z2) * (z1 - z2));
}

/** Calculates the latitude and longitude given a position and
  * geodetic distance and azimuth.
  * http://www.anzlic.orig.au/icsm/gdatm/chapter4.htm
  *
  * @param x1 Longitude at point 1.
  * @param y1 Latitude at point 1.
  * @param az1 Forward geodetic angle.
  * @param s Geodetic distance.
  * @param a Semi-major axis of Earth.
  * @param e2 First eccentricity squared of Earth.
  * @param x2 Pointer to computed longitude.
  * @param y2 Pointer to computed latitude.
  */
void geod_fwd_sodanos(double x1, double y1, double az1, double s,
                      double a, double e2, double *x2, double *y2)
{
  double term1, term2, term3, term4, term5;
  double phi0, sinphi0, cosphi0, v, omega, sinbeta2, cosbeta2;
  double e4 = e2 * e2;
  double e22 = e2 / 2;
  double e24 = e2 / 4;
  double f = (1.0 - sqrt(1.0 - e2));
  double f2 = f * f;
  double b = a * (1 - f);
  double beta1 = atan(tan(y1) * (1 - f));
  double sinbeta1 = sin(beta1);
  double sin2beta1 = sinbeta1 * sinbeta1;
  double cosbeta1 = cos(beta1);
  double saz = sin(az1);
  double caz = cos(az1);
  double cosbeta0 = cosbeta1 * saz;
  double g = cosbeta1 * caz;
  double m1 = (1 + e22 * sin2beta1) * (1 - cosbeta0 * cosbeta0);
  double phis = s / b;
  double sinphis = sin(phis);
  double cosphis = cos(phis);
  double cos2phis = cosphis * cosphis;
  double a1 =
    (1 + e22 * sin2beta1) * (sin2beta1 * cosphis + g * sinbeta1 * sinphis);

  term1 = a1 * (-e22 * sinphis);
  term2 = m1 * (-e24 * phis + e24 * sinphis * cosphis);
  term3 = a1 * a1 * ((5 * e4 / 8) * sinphis * cosphis);
  term4 =
    m1 * m1 * ((11 * e4 / 64) * phis - (13 * e4 / 64) * sinphis * cosphis -
               (e4 / 8) * phis * cos2phis +
               (5 * e4 / 32) * sinphis * cos2phis * cosphis);
  term5 =
    a1 * m1 * ((3 * e4 / 8) * sinphis + (e4 / 4) * phis * cosphis -
               (5 * e4 / 8) * sinphis * cos2phis);
  phi0 = phis + term1 + term2 + term3 + term4 + term5;
  sinphi0 = sin(phi0);
  cosphi0 = cos(phi0);

/*   az2 = atan(cosbeta0/(g*cosphi0 - sinbeta1*sinphi0)); */
  v =
    atan((sinphi0 * saz) /
         (cosbeta1 * cosphi0 - sinbeta1 * sinphi0 * caz));
  term1 = -f * phis;
  term2 = a1 * ((3 * f2 / 2) * sinphis);
  term3 = m1 * ((3 * f2 / 4) * phis - (3 * f2 / 4) * sinphis * cosphis);
  omega = (term1 + term2 + term3) * cosbeta0 + v;

  sinbeta2 = sinbeta1 * cosphi0 + g * sinphi0;
  term1 = cosbeta0;
  term2 = g * cosphi0 - sinbeta1 * sinphi0;
  cosbeta2 = sqrt(term1 * term1 + term2 * term2);

  *x2 = x1 + omega;
  *y2 = atan((sinbeta2 / cosbeta2) / (1 - f));
}


/** Compute the geodetic distance using the Sodano's algorithm
  * http://www.anzlic.orig.au/icsm/gdatm/chapter4.htm
  *
  * @param x1 Longitude at point 1.
  * @param y1 Latitude at point 1.
  * @param x2 Longitude at point 2.
  * @param y2 Latitude at point 2.
  * @param a Semi-major axis of Earth.
  * @param e2 First eccentricity squared of Earth.
  */
double geod_inv_geod_fwd_sodanos(double x1, double y1, double x2,
                                 double y2, double a, double e2)
{
  double term1, term2, term3, term4, term5, term6;
  double sinphi, cosphi, tanphi, phi, c, m;
  double f = (1.0 - sqrt(1.0 - e2));
  double f1 = f + f * f;
  double f2 = f * f / 2;
  double f3 = f1 / 2;
  double f4 = f * f / 8;
  double f5 = f * f / 16;
  double beta1 = atan(tan(y1) * (1 - f));
  double beta2 = atan(tan(y2) * (1 - f));
  double sinbeta1 = sin(beta1);
  double sinbeta2 = sin(beta2);
  double cosbeta1 = cos(beta1);
  double cosbeta2 = cos(beta2);
  double omega = x2 - x1;
  double sinomega = sin(omega);
  double cosomega = cos(omega);
  double b = a * (1 - f);

  term1 = sinomega * cosbeta2;
  term2 = sinbeta2 * cosbeta1 - sinbeta1 * cosbeta2 * cosomega;
  sinphi = sqrt(term1 * term1 + term2 * term2);
  phi = asin(sinphi);
  cosphi = cos(phi);
  tanphi = sinphi / cosphi;
  c = cosbeta1 * cosbeta2 * sinomega / sinphi;
  m = 1.0 - c * c;

  term1 = phi * (1 + f1);
  term2 = sinbeta1 * sinbeta2 * (f1 * sinphi - f2 * phi * phi / sinphi);
  term3 = m * (-phi * f3 - f3 * sinphi * cosphi + f2 * phi * phi / tanphi);
  term4 =
    (sinbeta1 * sinbeta2) * (sinbeta1 * sinbeta2) * (-f2 * sinphi *
                                                     cosphi);
  term5 =
    m * m * (f5 * phi + f5 * sinphi * cosphi - f2 * phi * phi / tanphi -
             f4 * sinphi * cosphi * cosphi * cosphi);
  term6 =
    sinbeta1 * sinbeta2 * m * (f2 * phi * phi / sinphi +
                               f2 * sinphi * cosphi * cosphi);

  return b * (term1 + term2 + term3 + term4 + term5 + term6);
}

/** Compute the geodetic forward and reverse azimuths using the
  * Sodano's algorithm
  * http://www.anzlic.orig.au/icsm/gdatm/chapter4.htm
  *
  * @param x1 Longitude at point 1.
  * @param y1 Latitude at point 1.
  * @param x2 Longitude at point 2.
  * @param y2 Latitude at point 2.
  * @param a Semi-major axis of Earth.
  * @param e2 First eccentricity squared of Earth.
  * @param raz -
  */
double geod_inv_sodanos_angles(double x1, double y1, double x2, double y2,
                               double a, double e2, double *raz)
{
  double term1, term2, term3;
  double sinphi, cosphi, tanphi, phi, c, m, v, sinv, cosv, az1, az2;
  double f = (1.0 - sqrt(1.0 - e2));
  double f1 = f + f * f;
  double f2 = f * f / 2;
  double beta1 = atan(tan(y1) * (1 - f));
  double beta2 = atan(tan(y2) * (1 - f));
  double sinbeta1 = sin(beta1);
  double sinbeta2 = sin(beta2);
  double cosbeta1 = cos(beta1);
  double cosbeta2 = cos(beta2);
  double omega = x2 - x1;
  double sinomega = sin(omega);
  double cosomega = cos(omega);

  term1 = sinomega * cosbeta2;
  term2 = sinbeta2 * cosbeta1 - sinbeta1 * cosbeta2 * cosomega;
  sinphi = sqrt(term1 * term1 + term2 * term2);
  phi = asin(sinphi);
  cosphi = cos(phi);
  tanphi = sinphi / cosphi;
  c = cosbeta1 * cosbeta2 * sinomega / sinphi;
  m = 1.0 - c * c;

  term1 = f1 * phi;
  term2 =
    sinbeta1 * sinbeta2 * (-f2 * sinphi - f * f * phi * phi / sinphi);
  term3 =
    m * ((-5 * f2 / 2) * phi + (f2 / 2) * sinphi * cosphi +
         f * f * phi * phi / tanphi);
  v = c * (term1 + term2 + term3) + omega;
  sinv = sin(v);
  cosv = cos(v);

  az1 =
    atan2((sinv * cosbeta2),
          (sinbeta2 * cosbeta1 - cosv * sinbeta1 * cosbeta2));
  az2 =
    atan2((sinv * cosbeta1),
          (sinbeta2 * cosbeta1 * cosv - sinbeta1 * cosbeta2));

  if (raz != NULL)
    *raz = az2;

  return az1;
}



/** Compute the great circle based on spherical trigonometry.
  * Sourced from http://geodesy.auslig.gov.au/gopher/bbs/georef/distance.txt
  *
  * @param x1 Longitude at point 1.
  * @param y1 Latitude at point 1.
  * @param x2 Longitude at point 2.
  * @param y2 Latitude at point 2.
  * @return great circle distance on sphere.
  */
double geod_gc_distance(double x1, double y1, double x2, double y2)
{
  return 1852.0 * 60 * (180.0 / M_PI)
    * acos(sin(y1) * sin(y2) + cos(y1) * cos(y2) * cos(x2 - x1));
}


/** Computes the auxillary latitude/longitudes given
  * equitorial latitude/longitudes and a false pole.
  *
  * @param elon Equitorial longitude (radians).
  * @param elat Equitorial latitude (radians).
  * @param plon Longitude of auxillary pole (radians).
  * @param plat Latitude of auxillary pole (radians).
  * @param alon Auxillary longitude (radians).
  * @param alat Auxillary latitude (radians).
  */
void geod_fwd_spherical_rot(double elon, double elat, double plon,
                            double plat, double *alon, double *alat)
{
  double dl = elon - plon;
  double sinplat = sin(plat);
  double cosplat = cos(plat);
  double cosdl = cos(dl);
  *alon =
    2 * M_PI - atan2(sin(dl), (cosplat * tan(elat) - sinplat * cosdl));
  *alat = asin(sin(elat) * sinplat + cos(elat) * cosplat * cosdl);
}


/** Computes the equitorial latitude/longitudes given
  * auxilary latitude/longitudes and a false pole.
  *
  * @param alon Auxilary longitude (radians).
  * @param alat Auxilary latitude (radians).
  * @param plon Longitude of auxillary pole (radians).
  * @param plat Latitude of auxillary pole (radians).
  * @param elon Equitorial longitude (radians).
  * @param elat Equitorial latitude (radians).
  */
void geod_inv_spherical_rot(double alon, double alat, double plon,
                            double plat, double *elon, double *elat)
{
  double sinplat = sin(plat);
  double cosplat = cos(plat);
  double cosalon = cos(alon);
  alon = 2 * M_PI - alon;
  *elat = asin(sinplat * sin(alat) + cosplat * cos(alat) * cosalon);
  *elon =
    atan2(sin(alon), (cosplat * tan(alat) - sinplat * cosalon)) + plon;
}

#if TEST
void main(int argc, char *argv[])
{
  ems_init();
  double x1 = 150 * M_PI / 180.0;
  double y1 = -30 * M_PI / 180.0;
  double x2 = 151 * M_PI / 180.0;
  double y2 = -80 * M_PI / 180.0;
  double a = 6378160.0;
  double f = 1 / 298.25;
  double e2 = 2 * f - f * f;
  double lat, lon;
  double plat = 90 * M_PI / 180.0;
  double plon = 0 * M_PI / 180.0;
  double tlon, tlat;
  double az2;

  double rdist = geod_inv_robbins(x1, y1, 0.0, x2, y2, 0.0, a, e2);
  double sdist = geod_inv_geod_fwd_sodanos(x1, y1, x2, y2, a, e2);
  double az1 = geod_inv_sodanos_angles(x1, y1, x2, y2, a, e2, &az2);
  double gcdist = geod_gc_distance(x1, y1, x2, y2);

  emslog(LDEBUG,"Robbins distance = %f\n", rdist);
  emslog(LDEBUG,"Sodanos distance = %f\n", sdist);
  emslog(LDEBUG,"Sodanos angles = %f %f\n", az1 * 180 / M_PI, az2 * 180 / M_PI);
  emslog(LDEBUG,"Greate circle distance = %f\n", gcdist);

  geod_fwd_sodanos(x1, y1, 315 * M_PI / 180.0, 1000000.0, a, e2, &x2, &y2);
  emslog(LDEBUG,("Sodanos computed position = %f %f\n", x2 * 180.0 / M_PI,
         y2 * 180.0 / M_PI);

#if 0
  for (lon = -160.0; lon <= 160.0; lon += 20.0) {
    for (lat = -80.0; lat <= 80.0; lat += 20.0) {
      geod_fwd_spherical_rot(lon * M_PI / 180.0, lat * M_PI / 180.0,
                             plon, plat, &tlon, &tlat);
      emslog(LDEBUG,"%g %g %g %g", lon, lat,
             tlon * 180 / M_PI, tlat * 180 / M_PI);
      geod_inv_spherical_rot(tlon, tlat, plon, plat, &tlon, &tlat);
      emslog(LDEBUG," %g %g\n", tlon * 180 / M_PI, tlat * 180 / M_PI);
    }
  }


  lon = 122.65;
  lat = -25.53;
  geod_fwd_spherical_rot(lon * M_PI / 180.0, lat * M_PI / 180.0,
                         plon, plat, &tlon, &tlat);
  emslog(LDEBUG,"%g %g %g %g", lon, lat, tlon * 180 / M_PI, tlat * 180 / M_PI);
  geod_inv_spherical_rot(tlon + 1 * M_PI / 180, tlat + 1 * M_PI / 180,
                         plon, plat, &tlon, &tlat);
  emslog(LDEBUG," %g %g\n", tlon * 180 / M_PI, tlat * 180 / M_PI);
#endif
  ems_clear();
  exit(0);
}
#endif
