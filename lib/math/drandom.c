/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/drandom.c
 *
 *  \brief Random number generators
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: drandom.c 6161 2019-03-05 04:49:33Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "ems.h"

/** Random number generator, based on the routines described
  * in Numerical Recipes. Note that although this routine
  * returns a double, ran[123]() return a float, so that the
  * distribution of random numbers can't be more fine than a
  * float can produce. Also, from test done by John Hunter,
  * ran1 is poorly distributed when plotted in 4-space and also
  * takes a long time. It should probably be avoided.
  *
  * @param min minimum of number range
  * @param max maximum of number range
  * @return Random number in the range min to max
  */
double drandom(double min, double max)
{
  int init;

  /* ran3() initialises itself automatically, so we don't ever need to do
     this */

  init = 0;
  return (min + (max - min) * ran3(&init));
}


#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/** From Numerical Recipes in C, Press et al
 * Cambridge University Press 1990.
 *
 * @param idum pointer to seed value.
 */
float ran1(int *idum)
{
  static long ix1, ix2, ix3;
  static float r[98];
  float temp;
  static int iff = 0;
  int j;

  if (*idum < 0 || iff == 0) {
    iff = 1;
    ix1 = (IC1 - (*idum)) % M1;
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix2 = ix1 % M2;
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix3 = ix1 % M3;
    for (j = 1; j <= 97; j++) {
      ix1 = (IA1 * ix1 + IC1) % M1;
      ix2 = (IA2 * ix2 + IC2) % M2;
      r[j] = (float)((ix1 + ix2 * RM2) * RM1);
    }
    *idum = 1;
  }
  ix1 = (IA1 * ix1 + IC1) % M1;
  ix2 = (IA2 * ix2 + IC2) % M2;
  ix3 = (IA3 * ix3 + IC3) % M3;
  j = 1 + ((97 * ix3) / M3);
  if (j > 97 || j < 1)
    quit("RAN1: This cannot happen.");
  temp = r[j];
  r[j] = (float)((ix1 + ix2 * RM2) * RM1);
  return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3


#include <math.h>

#define M 714025
#define IA 1366
#define IC 150889

/** From Numerical Recipes in C, Press et al
 * Cambridge University Press 1990.
 *
 * @param idum pointer to seed value.
 */
float ran2(long int *idum)
{
  static long iy, ir[98];
  static int iff = 0;
  int j;

  if (*idum < 0 || iff == 0) {
    iff = 1;
    if ((*idum = (IC - (*idum)) % M) < 0)
      *idum = -(*idum);
    for (j = 1; j <= 97; j++) {
      *idum = (IA * (*idum) + IC) % M;
      ir[j] = (*idum);
    }
    *idum = (IA * (*idum) + IC) % M;
    iy = (*idum);
  }
  j = (int)(1 + 97.0 * iy / M);
  if (j > 97 || j < 1)
    quit("RAN2: This cannot happen.");
  iy = ir[j];
  *idum = (IA * (*idum) + IC) % M;
  ir[j] = (*idum);
  return (float)iy / M;
}

#undef M
#undef IA
#undef IC

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/** From Numerical Recipes in C, Press et al
 * Cambridge University Press 1990.
 *
 * @param idum pointer to seed value.
 */
float ran3(int *idum)
{
  static int inext, inextp;
  static long ma[56];
  static int iff = 0;
  long mj, mk;
  int i, ii, k;

  if (*idum < 0 || iff == 0) {
    iff = 1;
    mj = MSEED - (*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (k = 1; k <= 4; k++)
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    inext = 0;
    inextp = 31;
    *idum = 1;
  }
  if (++inext == 56)
    inext = 1;
  if (++inextp == 56)
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ)
    mj += MBIG;
  ma[inext] = mj;
  return ((float)(mj * FAC));
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
