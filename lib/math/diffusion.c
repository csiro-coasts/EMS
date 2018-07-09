/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/diffusion.c
 *
 *  \brief 1D diffusion
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: diffusion.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "ems.h"

/* Prototype for tri-diagonal system solver */
void solvetri(double *Cm1, double *C, double *Cp1, double *rhs,
              double *x, int imin, int imax);


/** Calculates 1 time step for the diffusion equation.
  *
  * \f[ \frac{\partial C}{\partial dt} = K \frac{\partial ^2 C}{\partial x^2} \f]
  *		    
  * This routine uses a weighted method using
  * a times the values at the new time step plus
  * (1-a) times the values at the current time
  * step to calculate the spatial derivatives.
  * This leads to a tri-diagonal system of
  * equations to solve for C.
  *
  * Note that consistent units must be used for
  * all these quantities. If SI units are used,
  * then the units are as follows:
  *
  * * x    - metres (m)
  * * dt   - seconds (s)
  * * k[0] - units per square metre per second (m-2s-1)
  * * k[n] - units per square metre per second (m-2s-1)
  * * k[i] - metres squared per second (m2s-1)
  *    
  * @param n number of concentration values
  * @param c array of concentration values. (n values)
  * @param xc array of x coords at c points (n values)
  * @param k array of diffusion coefficients. This array has n+1 values.
  *	 Values k[0] and k[n] are the boundary fluxes in +ve
  *	 x direction. Other values k[i] are diffusion coefficients at
  *	 xk[i], which is between points xc[i] and xc[i-1]
  * @param xk array of x coords at k points (n+1 values)
  * @param dt time step
  * @param a weighting value, in the range [0,1].
  *	  (a = 1.0) == fully implicit,
  *	  (a = 0.5) == semi implicit
  *	  (a = 0.0) == explicit
  */
void
diffusion1d(int n, double *c, double *xc, double *k,
            double *xk, double dt, double a)
{
  int i;
  double dx;
  double v;
  double *Cm1;
  double *C;
  double *Cp1;
  double *rhs;


  /* Sanity checks */
  if (n < 1)
    quit("diffusion1d: n < 1 (no points!)\n");
  if (a < 0.0 || a > 1.0)
    quit("diffusion1d: weight value must be in range [0,1]\n");

  /* Only 1 layer - simple calculation */
  if (n == 1) {
    c[0] += dt * (k[0] - k[1]) / (xk[1] - xk[0]);
    return;
  }

  /* Allocate temporary storage */
#if HAVE_ALLOCA
  Cm1 = (double *)alloca((n + 1) * sizeof(double));
  C = (double *)alloca((n + 1) * sizeof(double));
  Cp1 = (double *)alloca((n + 1) * sizeof(double));
  rhs = (double *)alloca((n + 1) * sizeof(double));
#else
  Cm1 = d_alloc_1d(n + 1);
  C = d_alloc_1d(n + 1);
  Cp1 = d_alloc_1d(n + 1);
  rhs = d_alloc_1d(n + 1);
#endif

  /* i=0 boundary */
  i = 0;
  dx = xk[i + 1] - xk[i];
  v = dt * k[i + 1] / ((xc[i + 1] - xc[i]) * dx);
  Cm1[i] = 0.0;
  C[i] = 1 + a * v;
  Cp1[i] = -a * v;
  rhs[i] =
    (1 - (1 - a) * v) * c[i] + (1 - a) * v * c[i + 1] + dt * k[i] / dx;

  /* middle points */
  for (i = 1; i < n - 1; i++) {
    double dxi = xk[i + 1] - xk[i];
    double vm = dt * k[i] / ((xc[i] - xc[i - 1]) * dxi);
    double vp = dt * k[i + 1] / ((xc[i + 1] - xc[i]) * dxi);

    Cm1[i] = -a * vm;
    C[i] = 1.0 + a * (vm + vp);
    Cp1[i] = -a * vp;
    rhs[i] = (1 - a) * vm * c[i - 1]
      + (1 - (1 - a) * (vm + vp)) * c[i]
      + (1 - a) * vp * c[i + 1];
  }

  /* i=n-1 boundary */
  i = n - 1;
  dx = xk[i + 1] - xk[i];
  v = dt * k[i] / ((xc[i] - xc[i - 1]) * dx);
  Cm1[i] = -a * v;
  C[i] = 1 + a * v;
  Cp1[i] = 0.0;
  rhs[i] =
    (1 - a) * v * c[i - 1] + (1 - (1 - a) * v) * c[i] - dt * k[i + 1] / dx;

  /* Solve the system */
  solvetri(Cm1, C, Cp1, rhs, c, 0, n - 1);

#if !HAVE_ALLOCA
  /* Free temporary storage */
  d_free_1d(Cm1);
  d_free_1d(C);
  d_free_1d(Cp1);
  d_free_1d(rhs);
#endif
}


/**
 *  Routine to solve tridiagonal system of equations
 *  Arguments:
 * 
 * * Cm1[i] - coefficient of Xi-1 in ith eqn (lower diagonal)
 * * C[i]   - coefficients of Xi in ith eqn  (diagonal)
 * * Cp1[i] - coefficients of Xi+1 in ith eqn (upper diagonal)
 * * rhs[i] - right hand side of ith eqn
 * * x[i]   - where to store solved Xi
 * * imin   - minimum index
 * * imax   - maximum index
*/
void solvetri(double *Cm1, double *C, double *Cp1, double *rhs,
              double *x, int imin, int imax)
{
  int i;
  double div;
  double *ud;

  /* check indices - note imin == imax => 1 trivial eqn */
  if (imin >= imax || imin < 0)
    quit("solvetri: bad index values\n");

  /* Allocate temporary storage */
#if HAVE_ALLOCA
  ud = (double *)alloca((imax + 1) * sizeof(double));
#else
  ud = d_alloc_1d(imax + 1);
#endif

  div = C[imin];
  if (div == 0.0)
    quit("solvetri: zero first coefficient\n");
  x[imin] = rhs[imin] / div;
  for (i = imin + 1; i <= imax; i++) {
    ud[i] = Cp1[i - 1] / div;
    div = C[i] - Cm1[i] * ud[i];
    if (div == 0.0)
      quit("solvetri: zero divisor\n");
    x[i] = (rhs[i] - Cm1[i] * x[i - 1]) / div;
  }
  for (i = imax - 1; i >= imin; i--)
    x[i] -= ud[i + 1] * x[i + 1];

#if !HAVE_ALLOCA
  /* Free temporary storage */
  d_free_1d(ud);
#endif
}
