/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/math/decay.c
 *
 *  \brief Exponential decay routines
 * 
 *   The decay equation is \f[C = C_o e^{-kt} \f] so that 
 *   \f[ \frac{dC}{dt} = -kC \f]
 *
 *  This file contains routines to implement one time step of the above
 *  equation using forward, centered or backward time steps, as well as an
 *  exact method requiring a call to exp().
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: decay.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdio.h>
#include <math.h>

/** Decay using a forward time step. Numerically, this is as follows:
  *
  * \f[C(t+dt) = C(t)*(1-k dt)\f]
  *
  * @param c concentration
  * @param k decay rate
  * @param dt decay interval
  */
double decay_forward(double c, double k, double dt)
{
  return (c * (1.0 - k * dt));
}


/** Decay using a centered time step. Numerically, this is as follows:
  * 
  *
  * \f[ C(t+dt) = \frac{2-k dt}{2+k dt} C(t) \f]
  *
  */
double decay_centered(double c, double k, double dt)
{
  return (c * (2.0 - k * dt) / (2.0 + k * dt));
}

/** Decay using a backward time step. Numerically, this is as follows:
  * 
  * \f[ C(t+dt) = \frac{1}{1+kdt} C(t) \f]
  * 
  */
double decay_backward(double c, double k, double dt)
{
  return (c / (1.0 + k * dt));
}

/** Decay using an exact time step. Numerically, this is as follows:
  * 
  * \f[ C(t+dt) = C(t) e^{-kdt} \f]
  * 
  */
double decay_exact(double c, double k, double dt)
{
  return c * exp(-k * dt);
}
