/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/emsmath.h
 *
 *  \brief EMS maths prototypes.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: emsmath.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef	_EMS_MATH_H
#define _EMS_MATH_H

#include <limits.h>
#include <math.h>

/* Useful defines */
#if !defined(PI)
#define PI (3.14159265358979323846)
#endif

#if !defined(M_PI)
#define M_PI (3.14159265358979323846)
#endif

#if !defined(max)
#define max(x,y) ((x)>(y) ? (x) : (y) )
#endif

#if !defined(min)
#define min(x,y) ((x)<(y) ? (x) : (y) )
#endif

#ifdef _WIN32
#include <float.h>
#if !defined(__MINGW32__)
#define isnan _isnan
#else

/* Win - MinGW  the obsolete values.h doesn't define anything */
#if defined(HAVE_LONG_LONG) && !defined(LONGLONG_MIN)
#define LONGLONG_MIN	((long long) 0x8000000000000000LL)
#define LONGLONG_MAX	((long long) 0x7FFFFFFFFFFFFFFFLL)
#endif

#ifndef DBL_MIN
#define DBL_MIN		4.94065645841246544e-324
#define FLT_MIN		((float)1.40129846432481707e-45)
#endif
#ifndef DBL_MAX
#define DBL_MAX		1.79769313486231470e+308
#define FLT_MAX		((float)3.40282346638528860e+38)
#endif

#ifndef MAXDOUBLE
#define MAXDOUBLE DBL_MAX
#endif
#ifndef MINDOUBLE
#define MINDOUBLE DBL_MIN
#endif

#endif
#endif

#include "integrator.h"

/* UR disabled in favor of nan.h */
/*
#if defined(__ICC) || defined(__ICL) || defined(__ECC) || defined(__ECL)
double setNaN(void);
#define NaN setNaN()
#else
extern double NaN;
#endif
*/

/* Prototypes */
double decay_forward(double c, double k, double dt);
double decay_centered(double c, double k, double dt);
double decay_backward(double c, double k, double dt);
double decay_exact(double c, double k, double dt);
void diffusion1d(int n, double *c, double *xc,
                 double *k, double *xk, double dt, double a);
void solvetri(double *Cm1, double *C, double *Cp1, double *rhs, double *x,
              int imin, int imax);
double wgt_tophat(double x, double scale);
double wgt_linear(double x, double scale);
double wgt_parabolic(double x, double scale);
double wgt_gaussian(double x, double scale);
double wgt_gaussian_2d(double x, double y, double scale);
double drandom(double min, double max);
float ran3(int *init);



#endif                          /* _EMS_H */
