/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/moon.c
 *
 *  \brief Moon angle calculations
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: moon.c 6539 2020-05-01 01:53:02Z bai155 $
 */

#include <stdlib.h>
#include <math.h>
#include "ems.h"

/*
 * The following routines were moved from model/hd-us/boundary/tide.c to here 
 * so that other libraries (eg. ecology) may also make use of it
 */

static double sgn(double a)
{
  double ret;
  if (a > 0.0) ret = 1.0;
  if (a== 0.0) ret = 0.0;
  if (a < 0.0) ret = -1.0;
  return(ret);
}

/* Four quadrant inverse tangent */
/*  a = sine of angle */
/*  b = cosine of angle */
/* y = angle (radians; 0 =< c <= 2 * pi) */
static double atan3 (double a, double b)
{
double epsilon = 0.0000000001;
double pidiv2 = 0.5 * PI;
 double y, c;

 if (fabs(a) < epsilon) {
    y = (1. - sgn(b)) * pidiv2;
    return(y);
 } else {
    c = (2. - sgn(a)) * pidiv2;
 }

 if (fabs(b) < epsilon) {
    y = c;
    return(y);
 } else
    y = c + sgn(a) * sgn(b) * (fabs(atan(a / b)) - pidiv2);
 return(y);
}

/*-------------------------------------------------------------------*/
/* lunar ephemeris                                                   */
/* input jdate = julian date                                         */
/* output                                                            */
/*  rasc  = right ascension of the moon (radians)                    */
/*          (0 <= rasc <= 2 pi)                                      */
/*  decl  = declination of the moon (radians)                        */
/*          (-pi/2 <= decl <= pi/2)                                  */
/* https://au.mathworks.com/matlabcentral/fileexchange/              */
/*       39191-low-precision-ephemeris?focused=3789856&tab=function  */
/*-------------------------------------------------------------------*/
void moonvars(double jdate, 
	      double *rasc, 
	      double *decl, 
	      double *plon, 
	      double *plat,
	      double *dist
	      )
{
double atr = PI / 648000.;
/* time arguments                                                    */
double djd = jdate - 2451545.;
double t = (djd / 36525.) + 1.;
/* fundamental trig arguments (radians)                              */
double gm = r2r(0.374897 + 0.03629164709 * djd);
double gm2 = 2. * gm;
double gm3 = 3. * gm;
double fm = r2r(0.259091 + 0.0367481952 * djd);
double fm2 = 2. * fm;
double em = r2r(0.827362 + 0.03386319198 * djd);
double em2 = 2. * em;
double em4 = 4. * em;
double gs = r2r(0.993126 + 0.0027377785 * djd);
double lv = r2r(0.505498 + 0.00445046867 * djd);
double lm = r2r(0.606434 + 0.03660110129 * djd);
double ls = r2r(0.779072 + 0.00273790931 * djd);
double rm = r2r(0.347343 - 0.00014709391 * djd);
/* geocentric, ecliptic longitude of the moon (radians)              */
 double a, b, obliq, r;
double l = 22640. * sin(gm) - 4586. * sin(gm - em2) + 2370. * sin(em2);
 l = l + 769. * sin(gm2) - 668. * sin(gs) - 412. * sin(fm2);
 l = l - 212. * sin(gm2 - em2) - 206. * sin(gm - em2 + gs);
 l = l + 192. * sin(gm + em2) + 165. * sin(em2 - gs);
 l = l + 148. * sin(gm - gs) - 125. * sin(em) - 110. * sin(gm + gs);
 l = l - 55. * sin(fm2 - em2) - 45. * sin(gm + fm2) + 40. * sin(gm - fm2);
 l = l - 38. * sin(gm - em4) + 36. * sin(gm3) - 31. * sin(gm2 - em4);
 l = l + 28. * sin(gm - em2 - gs) - 24. * sin(em2 + gs) + 19. * sin(gm - em);
 l = l + 18. * sin(em + gs) + 15. * sin(gm + em2 - gs) + 14. * sin(gm2 + em2);
 l = l + 14. * sin(em4) - 13. * sin(gm3 - em2) - 17. * sin(rm);
 l = l - 11. * sin(gm + 16. * ls - 18. * lv) + 10. * sin(gm2 - gs) + 
   9. * sin(gm - fm2 - em2);
 l = l + 9. * (cos(gm + 16. * ls - 18. * lv) - sin(gm2 - em2 + gs)) -
   8. * sin(gm + em);
 l = l + 8. * (sin(2. * (em - gs)) - sin(gm2 + gs)) - 7. * (sin(2. * gs) + 
   sin(gm - 2. * (em - gs)) - sin(rm));
 l = l - 6. * (sin(gm - fm2 + em2) + sin(fm2 + em2)) -
   4. * (sin(gm - em4 + gs) - t * cos(gm + 16. * ls - 18. * lv));
 l = l - 4. * (sin(gm2 + fm2) - t * sin(gm + 16. * ls - 18. * lv));
 l = l + 3. * (sin(gm - 3. * em) - sin(gm + em2 + gs) -
    sin(gm2 - em4 + gs) + sin(gm - 2. * gs) + sin(gm - em2 - 2. * gs));
 l = l - 2. * (sin(gm2 - em2 - gs) + sin(fm2 - em2 + gs) - sin(gm + em4));
l = l + 2. * (sin(4. * gm) + sin(em4 - gs) + sin(gm2 - em));
 *plon = lm + atr * l;
 /* geocentric, ecliptic latitude of the moon (radians)              */
 b = 18461. * sin(fm) + 1010. * sin(gm + fm) + 1000. * sin(gm - fm);
 b = b - 624. * sin(fm - em2) - 199. * sin(gm - fm - em2) -
    167. * sin(gm + fm - em2);
 b = b + 117. * sin(fm + em2) + 62. * sin(gm2 + fm) + 33. * sin(gm - fm + em2);
 b = b + 32. * sin(gm2 - fm) - 30. * sin(fm - em2 + gs) -
    16. * sin(gm2 - em2 + fm);
 b = b + 15. * sin(gm + fm + em2) + 12. * sin(fm - em2 - gs) -
   9. * sin(gm - fm - em2 + gs);
 b = b - 8. * (sin(fm + rm) - sin(fm + em2 - gs)) -
    7. * sin(gm + fm - em2 + gs);
 b = b + 7. * (sin(gm + fm - gs) - sin(gm + fm - em4));
 b = b - 6. * (sin(fm + gs) + sin(3. * fm) - sin(gm - fm - gs));
 b = b - 5. * (sin(fm + em) + sin(gm + fm + gs) + sin(gm - fm + gs) -
	       sin(fm - gs) - sin(fm - em));
 b = b + 4. * (sin(gm3 + fm) - sin(fm - em4)) - 3. * (sin(gm - fm - em4) -
						      sin(gm - 3. * fm));
 b = b - 2. * (sin(gm2 - fm - em4) + sin(3. * fm - em2) - sin(gm2 - fm + em2) -
	       sin(gm - fm + em2 - gs));
 *plat = atr * (b + 2. * (sin(gm2 - fm - em2) + sin(gm3 - fm)));
 /* obliquity of the ecliptic (radians)                              */
 obliq = atr * (84428. - 47. * t + 9. * cos(rm));
 /* geocentric distance (kilometers)                                 */
 r = 60.36298 - 3.27746 * cos(gm) - .57994 * cos(gm - em2);
 r = r - .46357 * cos(em2) - .08904 * cos(gm2) + .03865 * cos(gm2 - em2);
 r = r - .03237 * cos(em2 - gs) - .02688 * cos(gm + em2) - 
   .02358 * cos(gm - em2 + gs);
 r = r - .0203 * cos(gm - gs) + .01719 * cos(em) + .01671 * cos(gm + gs);
 r = r + .01247 * cos(gm - fm2) + .00704 * cos(gs) + .00529 * cos(em2 + gs);
 r = r - .00524 * cos(gm - em4) + .00398 * cos(gm - em2 - gs) -
   .00366 * cos(gm3);
 r = r - .00295 * cos(gm2 - em4) - .00263 * cos(em + gs) +
   .00249 * cos(gm3 - em2);
 r = r - .00221 * cos(gm + em2 - gs) + .00185 * cos(fm2 - em2) -
    .00161 * cos(2 * (em - gs));
 r = r + 0.00147 * cos(gm + fm2 - em2) - 0.00142 * cos(em4) +
   0.00139 * cos(gm2 - em2 + gs);
 *dist = 6378.14 * (r - 0.00118 * cos(gm - em4 + gs) - 0.00116 * cos(gm2 + em2) -
		    0.0011 * cos(gm2 - gs));
 /* geocentric, equatorial right ascension and declination (radians) */
 a = sin(*plon) * cos(obliq) - tan(*plat) * sin(obliq);
 b = cos(*plon);
 *rasc = atan3(a, b);
 *decl = asin(sin(*plat) * cos(obliq) + cos(*plat) * sin(obliq) * sin(*plon));
 // printf("Moonvars: jdate %e, rasc %e, decl %e \n",jdate,*rasc,*decl);
}

