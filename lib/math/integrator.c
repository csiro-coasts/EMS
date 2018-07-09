/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/integrator.c
 *
 *  \brief ODE solvers
 *
 *
 *  Contains a number of ODE integrators, including:
 *   - dopri8(): Solution of ODE by 7/8 order Dorman-Prince method
 *   - dopri5(): Solution of ODE by 4/5 order Dorman-Prince method
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: integrator.c 5836 2018-06-27 03:17:23Z riz008 $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ems.h"

static void stats(int nfcn, int nstep, int naccpt, int nrejct);

/** dopri8():
 *
 * Numerical solution of a system of first order ordinary differential
 * equations Y'=F(X,Y). This is an embedded Runge-Kutta method of order7(8)
 * due to Dormand & Prince (with stepsize control).
 *
 * Ported from FORTRAN77 code. E.Hairer, S.P. Norsett, G. Wanner, "Solving
 * Ordinary Differential Equations I. Nonstiff Problems", Springer-Verlag,
 * 1987.
 *
 * Excellent for intermediate precisions (1e-8 to 1e-13) and smooth functions.
 *
 * The call to `dopri8()' normally returns 1; 0 otherwise.
 * The integration results are stored in the function vector passed. The
 * procedure also stores the last valid stepsize in the initial stepsize
 * guess value passed so it could be used at the next call if necessary.
 *
 * NOTICE. The derivative arrays supplied as arguments to `calc' may contain
 * some noise, it is a duty of the derivative calculation routine to set ALL
 * of the return values.
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/* Smallest number satisfying 1.0 + UROUND > 1.0. */
#define UROUND 1.11e-16

/* Maximal number of steps. */
#define NMAX 2000

/* Coefficients. */
#define C2    (1.0 / 18.0)
#define C3    (1.0 / 12.0)
#define C4    (1.0 / 8.0)
#define C5    (5.0 / 16.0)
#define C6    (3.0 / 8.0)
#define C7    (59.0 / 400.0)
#define C8    (93.0 / 200.0)
#define C9    (5490023248.0 / 9719169821.0)
#define C10   (13.0 / 20.0)
#define C11   (1201146811.0 / 1299019798.0)
#define C12   (1.0)
#define C13   (1.0)
#define A21   (C2)
#define A31   (1.0 / 48.0)
#define A32   (1.0 / 16.0)
#define A41   (1.0 / 32.0)
#define A43   (3.0 / 32.0)
#define A51   (5.0 / 16.0)
#define A53   (-75.0 / 64.0)
#define A54   (-A53)
#define A61   (3.0 / 80.0)
#define A64   (3.0 / 16.0)
#define A65   (3.0 / 20.0)
#define A71   (29443841.0 / 614563906.0)
#define A74   (77736538.0 / 692538347.0)
#define A75   (-28693883.0 / 1125.0e+6)
#define A76   (23124283.0 / 18.0e+8)
#define A81   (16016141.0 / 946692911.0)
#define A84   (61564180.0 / 158732637.0)
#define A85   (22789713.0 / 633445777.0)
#define A86   (545815736.0 / 2771057229.0)
#define A87   (-180193667.0 / 1043307555.0)
#define A91   (39632708.0 / 573591083.0)
#define A94   (-433636366.0 / 683701615.0)
#define A95   (-421739975.0 / 2616292301.0)
#define A96   (100302831.0 / 723423059.0)
#define A97   (790204164.0 / 839813087.0)
#define A98   (800635310.0 / 3783071287.0)
#define A101  (246121993.0 / 1340847787.0)
#define A104  (-37695042795.0 / 15268766246.0)
#define A105  (-309121744.0 / 1061227803.0)
#define A106  (-12992083.0 / 490766935.0)
#define A107  (6005943493.0 / 2108947869.0)
#define A108  (393006217.0 / 1396673457.0)
#define A109  (123872331.0 / 1001029789.0)
#define A111  (-1028468189.0 / 846180014.0)
#define A114  (8478235783.0 / 508512852.0)
#define A115  (1311729495.0 / 1432422823.0)
#define A116  (-10304129995.0 / 1701304382.0)
#define A117  (-48777925059.0 / 3047939560.0)
#define A118  (15336726248.0 / 1032824649.0)
#define A119  (-45442868181.0 / 3398467696.0)
#define A1110 (3065993473.0 / 597172653.0)
#define A121  (185892177.0 / 718116043.0)
#define A124  (-3185094517.0 / 667107341.0)
#define A125  (-477755414.0 / 1098053517.0)
#define A126  (-703635378.0 / 230739211.0)
#define A127  (5731566787.0 / 1027545527.0)
#define A128  (5232866602.0 / 850066563.0)
#define A129  (-4093664535.0 / 808688257.0)
#define A1210 (3962137247.0 / 1805957418.0)
#define A1211 (65686358.0 / 487910083.0)
#define A131  (403863854.0 / 491063109.0)
#define A134  (-5068492393.0 / 434740067.0)
#define A135  (-411421997.0 / 543043805.0)
#define A136  (652783627.0 / 914296604.0)
#define A137  (11173962825.0 / 925320556.0)
#define A138  (-13158990841.0 / 6184727034.0)
#define A139  (3936647629.0 / 1978049680.0)
#define A1310 (-160528059.0 / 685178525.0)
#define A1311 (248638103.0 / 1413531060.0)
#define B1    (14005451.0 / 335480064.0)
#define B6    (-59238493.0 / 1068277825.0)
#define B7    (181606767.0 / 758867731.0)
#define B8    (561292985.0 / 797845732.0)
#define B9    (-1041891430.0 / 1371343529.0)
#define B10   (760417239.0 / 1151165299.0)
#define B11   (118820643.0 / 751138087.0)
#define B12   (-528747749.0 / 2220607170.0)
#define B13   (1.0 / 4.0)
#define BH1   (13451932.0 / 455176623.0)
#define BH6   (-808719846.0 / 976000145.0)
#define BH7   (1757004468.0 / 5645159321.0)
#define BH8   (656045339.0 / 265891186.0)
#define BH9   (-3867574721.0 / 1518517206.0)
#define BH10  (465885868.0 / 322736535.0)
#define BH11  (53011238.0 / 667516719.0)
#define BH12  (2.0 / 45.0)
#endif
int dopri8(derivfn calc, /** function computing the first derivatives */
           int n,               /** dimension of the system */
           double x,            /** initial X-value */
           double *y,           /** Y-values */
           double xend,         /** final X-value */
           double eps,          /** local tolerance */
           double hmax,         /** maximal stepsize */
           double *h0,          /** initial stepsize guess */
           intout out,          /** output procedure */
           intfinal final,      /** final procedure */
           void *custom_data    /** any custom data */
	   )
{
  /* Number of function evaluations. */
  int nfcn = 0;
  /* Number of computed steps. */
  int nstep = 0;
  /* Number of accepted steps. */
  int naccpt = 0;
  /* Number of rejected steps. Please note that the number of rejected
     steps does not increase when Dopri8 has to reduce the initial step
     size, while the number of computed steps increases. Hence, it is
     possible to have nstep > naccpt + nrejct. */
  int nrejct = 0;
  /* Direction from x to xmax. */
  double posneg = (xend > x) ? 1.0 : -1.0;
  /* If the step has been rejected. */
  int reject = 0;
  /* Work arrays. */
  double *k = malloc(n * 8 * sizeof(double));

  double *k1 = k;
  double *k2 = &k[n];
  double *k3 = &k[n * 2];
  double *k4 = &k[n * 3];
  double *k5 = &k[n * 4];
  double *k6 = &k[n * 5];
  double *k7 = &k[n * 6];
  double *y1 = &k[n * 7];

  double h = *h0;
  /* Flag -- whether the step has been trancated to match the end point */
  int trunc = 0;                /* init to eliminate gcc warning */
  double xph, err, fac, hnew;
  int i;

  memset(k, 0, n * 8 * sizeof(double));

  eps = max(eps, 13.0 * UROUND);
  h = min(max(1.0e-10, fabs(h)), hmax) * posneg;

  if (out != NULL)
    out(1, n, x, y, custom_data);

  /* Main cycle. */
  while (((x - xend) * posneg + UROUND) <= 0.0) {
    if (nstep > NMAX) {
      emslog(LERROR, "dopri8(): Could not proceed: nstep > NMAX.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      return 0;
    }
    if ((x + 0.03 * h) == x) {
      emslog(LERROR,
              "dopri8(): Could not proceed: (x + 0.03 * step) == x.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      return 0;
    }
    if (!reject) {
      if (((x + h - xend) * posneg) > 0.0) {
        h = xend - x;
        trunc = 1;
      } else
        trunc = 0;

      /* First nine stages. */
      calc(x, y, k1, custom_data);
    }
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * A21 * k1[i];

    calc(x + C2 * h, y1, k2, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);

    calc(x + C3 * h, y1, k3, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (A41 * k1[i] + A43 * k3[i]);

    calc(x + C4 * h, y1, k4, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (A51 * k1[i] + A53 * k3[i] + A54 * k4[i]);

    calc(x + C5 * h, y1, k5, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (A61 * k1[i] + A64 * k4[i] + A65 * k5[i]);

    calc(x + C6 * h, y1, k6, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * (A71 * k1[i] + A74 * k4[i] + A75 * k5[i] + A76 * k6[i]);

    calc(x + C7 * h, y1, k7, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * (A81 * k1[i] + A84 * k4[i] + A85 * k5[i] + A86 * k6[i] +
                    A87 * k7[i]);

    calc(x + C8 * h, y1, k2, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * (A91 * k1[i] + A94 * k4[i] + A95 * k5[i] + A96 * k6[i] +
                    A97 * k7[i] + A98 * k2[i]);

    calc(x + C9 * h, y1, k3, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * (A101 * k1[i] + A104 * k4[i] + A105 * k5[i] +
                    A106 * k6[i] + A107 * k7[i] + A108 * k2[i] +
                    A109 * k3[i]);

    /* Compute intermediate sums. */
    for (i = 0; i < n; ++i) {
      double y11s =
        A111 * k1[i] + A114 * k4[i] + A115 * k5[i] + A116 * k6[i] +
        A117 * k7[i] + A118 * k2[i] + A119 * k3[i];
      double y12s =
        A121 * k1[i] + A124 * k4[i] + A125 * k5[i] + A126 * k6[i] +
        A127 * k7[i] + A128 * k2[i] + A129 * k3[i];
      k4[i] =
        A131 * k1[i] + A134 * k4[i] + A135 * k5[i] + A136 * k6[i] +
        A137 * k7[i] + A138 * k2[i] + A139 * k3[i];
      k5[i] =
        B1 * k1[i] + B6 * k6[i] + B7 * k7[i] + B8 * k2[i] + B9 * k3[i];
      k6[i] =
        BH1 * k1[i] + BH6 * k6[i] + BH7 * k7[i] + BH8 * k2[i] +
        BH9 * k3[i];
      k2[i] = y11s;
      k3[i] = y12s;
    }

    /* Last 4 stages. */
    calc(x + C10 * h, y1, k7, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (k2[i] + A1110 * k7[i]);

    calc(x + C11 * h, y1, k2, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (k3[i] + A1210 * k7[i] + A1211 * k2[i]);

    xph = x + h;
    calc(xph, y1, k3, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * (k4[i] + A1310 * k7[i] + A1311 * k2[i]);

    calc(xph, y1, k4, custom_data);
    for (i = 0; i < n; ++i) {
      k5[i] =
        h * (k5[i] + B10 * k7[i] + B11 * k2[i] + B12 * k3[i] +
             B13 * k4[i]);
      k6[i] =
        k5[i] - h * (k6[i] + BH10 * k7[i] + BH11 * k2[i] + BH12 * k3[i]);
      k5[i] += y[i];
    }

    nfcn += 13;

    /* Error estimation. */
    err = 0.0;
    for (i = 0; i < n; ++i) {
      double denom =
        max(max(1.0e-6, fabs(k5[i])), max(fabs(y[i]), 2.0 * UROUND / eps));
      denom = k6[i] / denom;
      err += denom * denom;
    }
    err = sqrt(err / n);

    if (isnan(err)) {
      emslog(LERROR, "dopri8: Stopped: NaN detected\n");
      return 0;
    }

    /* New step size. We require 0.333 <= hnew / w <= 6.0 */
    fac = max(1.0 / 6.0, min(3.0, pow(err / eps, 0.125) / 0.9));
    hnew = h / fac;

    if (err < eps) {
      naccpt++;
      for (i = 0; i < n; ++i)
        y[i] = k5[i];

      x = xph;
      if (out != NULL)
        out(naccpt + 1, n, x, y, custom_data);
      if (fabs(hnew) > hmax)
        hnew = posneg * hmax;
      if (reject)
        hnew = posneg * min(fabs(hnew), fabs(h));
      reject = 0;
      h = hnew;
      /* Store back the step size if this step has not been truncated to
         match the end point or it was the first step */
      if (naccpt == 1 || !trunc)
        *h0 = hnew;
    } else {
      reject = 1;
      h = hnew;
      if (naccpt >= 1)
        nrejct++;
      nfcn--;
    }

    nstep++;
  }

  if (final != NULL)
    final(n, x, y, k1, naccpt, nrejct, nfcn, custom_data);

  free(k);

 emslog(LMETRIC,"dopri8 iterations: %d \n",nstep);

  return 1;
}

#undef NMAX
#undef C2
#undef C3
#undef C4
#undef C5
#undef C6
#undef C7
#undef C8
#undef C9
#undef C10
#undef C11
#undef C12
#undef C13
#undef A21
#undef A31
#undef A32
#undef A41
#undef A43
#undef A51
#undef A53
#undef A54
#undef A61
#undef A64
#undef A65
#undef A71
#undef A74
#undef A75
#undef A76
#undef A81
#undef A84
#undef A85
#undef A86
#undef A87
#undef A91
#undef A94
#undef A95
#undef A96
#undef A97
#undef A98
#undef A101
#undef A104
#undef A105
#undef A106
#undef A107
#undef A108
#undef A109
#undef A111
#undef A114
#undef A115
#undef A116
#undef A117
#undef A118
#undef A119
#undef A1110
#undef A121
#undef A124
#undef A125
#undef A126
#undef A127
#undef A128
#undef A129
#undef A1210
#undef A1211
#undef A131
#undef A134
#undef A135
#undef A136
#undef A137
#undef A138
#undef A139
#undef A1310
#undef A1311
#undef B1
#undef B6
#undef B7
#undef B8
#undef B9
#undef B10
#undef B11
#undef B12
#undef B13
#undef BH1
#undef BH6
#undef BH7
#undef BH8
#undef BH9
#undef BH10
#undef BH11
#undef BH12

/** dopri5():
 *
 * Numerical solution of a system of first order ordinary differential
 * equations Y'=F(X,Y). This is an embedded Runge-Kutta method of order7(8)
 * due to Dormand & Prince (with stepsize control).
 *
 * Ported from FORTRAN77 code. E.Hairer, S.P. Norsett, G. Wanner, "Solving
 * Ordinary Differential Equations I. Nonstiff Problems", Springer-Verlag,
 * 1987.
 *
 * Best for low precisions (1e-4 to 1e-7) and smooth functions.
 *
 * The call to `dopri5()' normally returns 1; 0 otherwise.
 * The integration results are stored in the function vector passed. The
 * procedure also stores the last valid stepsize in the initial stepsize
 * guess value passed so it could be used at the next call if necessary.
 *
 * NOTICE. The derivative arrays supplied as arguments to `calc' may contain
 * some noise, it is a duty of the derivative calculation routine to set ALL
 * of the return values.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Maximal number of steps. */
#define NMAX 3000
#endif

int dopri5(derivfn calc,        /** function computing the first derivatives */
           int n,               /** dimension of the system */
           double x,            /** initial X-value */
           double *y,           /** Y-values */
           double xend,         /** final X-value */
           double eps,          /** local tolerance */
           double hmax,         /** maximal stepsize */
           double *h0,          /** initial stepsize guess */
           intout out,          /** output procedure */
           intfinal final,      /** final procedure */
           void *custom_data    /** any custom data */
	   )
{
  /* Number of function evaluations. */
  int nfcn = 1;
  /* Number of computed steps. */
  int nstep = 0;
  /* Number of accepted steps. */
  int naccpt = 0;
  /* Number of rejected steps. Please note that the number of rejected
     steps does not increase when Dopri8 has to reduce the initial step
     size, while the number of computed steps increases. Hence, it is
     possible to have nstep > naccpt + nrejct. */
  int nrejct = 0;
  /* Direction from x to xmax. */
  double posneg = (xend > x) ? 1.0 : -1.0;
  /* If the step has been rejected. */
  int reject = 0;
  /* Work arrays. */
  double *k = malloc(n * 6 * sizeof(double));

  double *k1 = k;
  double *k2 = &k[n];
  double *k3 = &k[n * 2];
  double *k4 = &k[n * 3];
  double *k5 = &k[n * 4];
  double *y1 = &k[n * 5];

  double h = *h0;
  /* Flag -- whether the step has been truncated to match the end point */
  int trunc;
  double xph, err, fac, hnew;
  int i;

  memset(k, 0, n * 6 * sizeof(double));
  eps = max(eps, 13.0 * UROUND);
  h = min(max(1.0e-10, fabs(h)), hmax) * posneg;

  if (out != NULL)
    out(1, n, x, y, custom_data);
  calc(x, y, k1, custom_data);

  /* Main cycle. */
  while (((x - xend) * posneg + UROUND) <= 0.0) {
    if (nstep > NMAX) {
      emslog(LERROR, "dopri5(): Could not proceed: nstep > NMAX.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      free(k);
      return 0;
    }
    if ((x + 0.1 * h) == x) {
      emslog(LERROR,
              "dopri5(): Could not proceed: (x + 0.1 * step) == x.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      free(k);
      return 0;
    }
    if (((x + h - xend) * posneg) > 0.0) {
      h = xend - x;
      trunc = 1;
    } else
      trunc = 0;

    /* First six stages. */
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * 0.2 * k1[i];

    calc(x + 0.2 * h, y1, k2, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] = y[i] + h * ((3.0 / 40.0) * k1[i] + (9.0 / 40.0) * k2[i]);

    calc(x + 0.3 * h, y1, k3, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * ((44.0 / 45.0) * k1[i] - (56.0 / 15.0) * k2[i] +
                    (32.0 / 9.0) * k3[i]);

    calc(x + 0.8 * h, y1, k4, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * ((19372.0 / 6561.0) * k1[i] -
                    (25360.0 / 2187.0) * k2[i] +
                    (64448.0 / 6561.0) * k3[i] - (212.0 / 729.0) * k4[i]);

    calc(x + (8.0 / 9.0) * h, y1, k5, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * ((9017.0 / 3168.0) * k1[i] - (355.0 / 33.0) * k2[i] +
                    (46732.0 / 5247.0) * k3[i] + (49.0 / 176.0) * k4[i] -
                    (5103.0 / 18656.0) * k5[i]);

    xph = x + h;
    calc(xph, y1, k2, custom_data);
    for (i = 0; i < n; ++i)
      y1[i] =
        y[i] + h * ((35.0 / 384.0) * k1[i] + (500.0 / 1113.0) * k3[i] +
                    (125.0 / 192.0) * k4[i] - (2187.0 / 6784.0) * k5[i] +
                    (11.0 / 84.0) * k2[i]);

    /* Compute intermediate sum. */
    for (i = 0; i < n; ++i)
      k2[i] =
        (71.0 / 57600.0) * k1[i] - (71.0 / 16695.0) * k3[i] +
        (71.0 / 1920.0) * k4[i] - (17253.0 / 339200.0) * k5[i] +
        (22.0 / 525.0) * k2[i];

    /* Last stage. */
    calc(xph, y1, k3, custom_data);
    for (i = 0; i < n; ++i)
      k4[i] = (k2[i] - (1.0 / 40.0) * k3[i]) * h;

    nfcn += 6;

    /* Error estimation. */
    err = 0.0;
    for (i = 0; i < n; ++i) {
      double denom =
        max(max(1.0e-5, fabs(y1[i])), max(fabs(y[i]), 2.0 * UROUND / eps));
      denom = k4[i] / denom;
      err += denom * denom;
    }
    err = sqrt(err / n);

    if (isnan(err)) {
      emslog(LERROR, "dopri5: Stopped: err as NaN detected\n");
      free(k);
      return 0;
    }

    /* New step size. We require 0.2 <= hnew / w <= 10.0 */
    fac = max(0.1, min(5.0, pow(err / eps, 0.2) / 0.9));
    hnew = h / fac;

    if (err < eps) {
      naccpt++;
      for (i = 0; i < n; ++i) {
        k1[i] = k3[i];
        y[i] = y1[i];
      }
      x = xph;
      if (out != NULL)
        out(naccpt + 1, n, x, y, custom_data);
      if (fabs(hnew) > hmax)
        hnew = posneg * hmax;
      if (reject)
        hnew = posneg * min(fabs(hnew), fabs(h));
      reject = 0;
      h = hnew;
      /* Store back the step size if this step has not been truncated to
         match the end point or it was the first step */
      if (naccpt == 1 || !trunc)
        *h0 = hnew;
    } else {
      reject = 1;
      h = hnew;
      if (naccpt >= 1)
        nrejct++;
    }

    nstep++;
  }

  if (final != NULL)
    final(n, x, y, k1, naccpt, nrejct, nfcn, custom_data);

  free(k);

  emslog(LMETRIC,"dopri5 iterations: %d \n",nstep);
  return 1;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef NMAX
#define NMAX 100000
#endif

/**
 * Adaptive method of the 2nd order
*/
int adapt2(derivfn calc,
           int n,
           double x,
           double *y,
           double xend,
           double eps,
           double hmax,
           double *h0, intout out, intfinal final, void *custom_data)
{
  double posneg = (xend > x) ? 1.0 : -1.0;

  /* Number of computed steps */
  int nstep = 0;
  /* Number of function evaluations */
  int nfcn = 0;
  /* Number of accepted steps */
  int naccpt = 0;
  /* If the step has been rejected. */
  int reject = 0;
  /* Work arrays */
  double *k = malloc(sizeof(double) * n * 5); /* derivative */
  double *k1 = &k[n];           /* derivative at half-step */
  double *y1 = &k[n * 2];       /* function at half-step */
  double *y2 = &k[n * 3];       /* function after two half-steps */
  double *w = &k[n * 4];        /* function after a full step */
  /* Number of rejected steps */
  int nrejct = 0;
  /* Flag -- whether the step has been trancated to match the end point */
  int trunc = 0;

  double h = *h0;
  h = min(max(1.0e-10, fabs(h)), hmax) * posneg;
  eps = max(eps, 1.0e-10);

  if (out != NULL)
    out(1, n, x, y, custom_data);

  while ((x - xend) * posneg + UROUND <= 0) {
    double h2 = h / 2.0;
    double err = 0.0;
    double fac;
    double hnew;
    int i;

    if (nstep > NMAX) {
      emslog(LERROR, "adapt2(): Could not proceed: nstep > NMAX.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      return 0;
    }

    if ((x + 0.03 * h) == x) {
      emslog(LERROR,
              "adapt2(): Could not proceed: (x + 0.03 * step) == x.\n");
      stats(nfcn, nstep, naccpt, nrejct);
      return 0;
    }

    if (!reject) {
      if (((x + h - xend) * posneg) > 0.0) {
        h = xend - x;
        trunc = 1;
      } else
        trunc = 0;
    }

    calc(x, y, k, custom_data);

    for (i = 0; i < n; ++i) {
      y1[i] = y[i] + h2 * k[i];
      w[i] = y[i] + h * k[i];
    }

    calc(x + h2, y1, k1, custom_data);

    nfcn += 2;

    for (i = 0; i < n; ++i) {
      double denom;
      y2[i] = y1[i] + h2 * k1[i];
      denom = max(max(fabs(y[i]), fabs(y2[i])), 1.0e-3);
      denom = (y2[i] - w[i]) / denom;
      err += denom * denom;
    }

    err = sqrt(err / n);

    fac = sqrt(err / eps) / 0.9;
    if (fac < 0.1)
      fac = 0.1;
    else if (fac > 10.0)
      fac = 10.0;

    hnew = h / fac;

    if (err < eps) {
      naccpt++;
      for (i = 0; i < n; ++i)
        y[i] = 2.0 * y2[i] - w[i];
      x += h;
      if (out != NULL)
        out(naccpt + 1, n, x, y, custom_data);
      if (fabs(hnew) > hmax)
        hnew = posneg * hmax;
      if (reject)
        hnew = posneg * min(fabs(hnew), fabs(h));
      reject = 0;
      h = hnew;
      /* Store back the step size if this step has not been truncated to
         match the end point or it was the first step */
      if (naccpt == 1 || !trunc)
        *h0 = hnew;
    } else {
      reject = 1;
      h = hnew;
      if (naccpt >= 1)
        nrejct++;
      nfcn--;
    }

    nstep++;
  }

  if (final != NULL)
    final(n, x, y, k1, naccpt, nrejct, nfcn, custom_data);

  free(k);
  emslog(LMETRIC,"adapt2 iterations: %d \n",nstep);
  return 1;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef NMAX
#define NMAX 200000
#endif

/* CUT_OFF is a parameter critical for performance of adapt1().
 * A value of 0.01 was used in the PPB study.
 *
 * PPB test results (precision = 0.1):
 *
 * CUT_OFF | Av. number of integrations per ecology step
 * --------|--------------------------------------------
 * 1e-2    | 128.8
 * 1e-3    | 128.9
 * 1e-4    | 131.6
 * 1e-5    | 141.0
 * 1e-6    | 189.3 (NMAX = 268889)
 */
#define CUT_OFF 1.0e-4

extern int *essential;

/*
 * Adaptive method of the 1st order used in a PPB study
 */
int adapt1(derivfn calc, /** function computing the first derivatives */
	   int n,        /** dimension of the system */
	   double x,     /** initial X-value */
	   double *y,    /** Y-values */
	   double xend,  /** final X-value */
	   double eps,   /** local tolerance */
	   double hmax,  /** maximal stepsize */
	   double *h0,   /** not relevant here */
           intout out,   /** output function */
	   intfinal final, /** final function */
	   void *custom_data /** any custom data */)
{
  double posneg = (xend > x) ? 1.0 : -1.0;

  /* Number of computed steps */
  int nstep = 0;
  /* Work array */
  double *k = malloc(sizeof(double) * n);

  double h = xend - x;

  int index = -1;               /* index of the most unstable variable */
  int i;

  eps = max(eps, 1.0e-4);
  *h0 = xend - x;               /* just in case -- to have a sensible
                                   value there */

  while ((x - xend) * posneg + UROUND <= 0) {
    double max_del = fabs(eps / h);

    if (nstep > NMAX) {
      emslog(LERROR, "adapt1(): Could not proceed: nstep > NMAX.\n");
      emslog(LERROR, "adapt1(): y[%d] = %e, y\'[%d] = %e.\n", index,
              y[index], index, k[index]);
      stats(nstep, nstep, nstep, 0);
      return 0;
    }

    if ((x + 0.03 * h) == x) {
      emslog(LERROR,
              "adapt1(): Could not proceed: (x + 0.03 * step) == x.\n");
      emslog(LERROR, "adapt1(): y[%d] = %e, y\'[%d] = %e.\n", index,
              y[index], index, k[index]);
      stats(nstep, nstep, nstep, 0);
      return 0;
    }

    calc(x, y, k, custom_data);

    /* Find the most unstable variable and its maximum change rate */
    for (i = 0; i < n; ++i) {
      if (essential[i]) {
        double local_del = fabs(k[i] / max(fabs(y[i]), CUT_OFF));
        if (local_del > max_del) {
          max_del = local_del;
          index = i;
        }
      }
    }

    h = posneg * eps / max_del;

    if (((x + h - xend) * posneg) > 0.0)
      h = xend - x;
    else
      *h0 = h;

    for (i = 0; i < n; ++i)
      y[i] += k[i] * h;
    x += h;
    nstep++;

    if (out != NULL)
      out(nstep, n, x, y, custom_data);
  }

  if (final != NULL)
    final(n, x, y, k, nstep, 0, nstep, custom_data);

  free(k);
  emslog(LMETRIC,"adapt1 iterations: %d \n",nstep);

  return 1;
}

#undef NMAX
#undef CUT_OFF
#undef UROUND

/*
 * Euler 1st order
 */
int euler1(derivfn calc, int n, double x, double *y, double xend, double eps, double hmax, double *h0, intout out, intfinal final, void *custom_data)
{
  double* k = calloc(n, sizeof(double));
  double h = xend - x;
  int i;

  calc(x, y, k, custom_data);
  for (i = 0; i < n; ++i)
    y[i] += k[i] * h;

  if (final != NULL)
    final(n, x, y, k, 1, 0, 1, custom_data);

  free(k);
  emslog(LMETRIC,"euler1 iterations: %d \n",i);
  return 1;
}

static void stats(int nfcn, int nstep, int naccpt, int nrejct)
{
  emslog(LINFO, "  Number of function evaluations = %d.\n", nfcn);
  emslog(LINFO, "  Number of computed steps = %d.\n", nstep);
  emslog(LINFO, "  Number of accepted steps = %d.\n", naccpt);
  emslog(LINFO, "  Number of rejected steps = %d.\n", nrejct);
}
