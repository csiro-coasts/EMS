/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/weight_fn.c
 *
 *  \brief Weighting functions
 *
 *  This file contains routines to
 *  calculate various functional forms.
 *  All routines are of the form
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: weight_fn.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>
#include "ems.h"

/** Compute's a tophat weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double wgt_tophat(double x,     /* Where to evaluate function */
                  double scale  /* Length scale for function */
  )
{
  if (fabs(x) <= scale)
    return (1.0);
  else
    return (0.0);
}

/** Compute's a linear weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double wgt_linear(double x,     /* Where to evaluate function */
                  double scale  /* Length scale for function */
  )
{
  if (fabs(x) <= scale)
    return (1.0 - x / scale);
  else
    return (0.0);
}

/** Compute's a parabolic weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double wgt_parabolic(double x,  /* Where to evaluate function */
                     double scale /* Length scale for function */
  )
{
  if (fabs(x) <= scale)
    return (1.0 - (x * x) / (scale * scale));
  else
    return (0.0);
}

/** Compute's a gaussian weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double wgt_gaussian(double x,   /* Where to evaluate function */
                    double scale  /* Length scale for function */
  )
{
  return (exp(-x * x / (2 * scale * scale)));
}

/** Compute's a 2D gaussian weighting function.
  *
  * @param x where to evaluate function.
  * @param y where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double wgt_gaussian_2d(double x,
		       double y,
		       double scale
  )
{
  return (exp(-(x*x + y*y) / (4 * scale * scale)));
}
