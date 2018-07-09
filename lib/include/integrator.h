/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/integrator.h
 *
 *  \brief Header file for ODE integration.
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: integrator.h 5834 2018-06-27 00:55:37Z riz008 $
 */


#if !defined(_INTEGRATOR_H)
#define _INTEGRATOR_H

typedef void (*derivfn) (double t, double *y, double *y1, void *p);

typedef void (*intout) (int nstep, int n, double x, double *y, void *p);

typedef void (*intfinal) (int n, double x, double *y, double *y1,
                          int naccpt, int nrejct, int nfcn, void *p);

typedef int (*integrator) (derivfn calc, int n, double x, double *y0,
                           double xend, double eps, double hmax, double *h,
                           intout out, intfinal final, void *p);

int dopri8(derivfn calc, int n, double x, double *y, double xend,
           double eps, double hmax, double *h0, intout out, intfinal final,
           void *custom_data);

int dopri5(derivfn calc, int n, double x, double *y, double xend,
           double eps, double hmax, double *h0, intout out, intfinal final,
           void *custom_data);

int adapt1(derivfn calc, int n, double x, double *y, double xend,
           double eps, double hmax, double *h0, intout out, intfinal final,
           void *custom_data);

int adapt2(derivfn calc, int n, double x, double *y, double xend,
           double eps, double hmax, double *h0, intout out, intfinal final,
           void *custom_data);

int euler1(derivfn calc, int n, double x, double *y, double xend,
	   double eps, double hmax, double *h0, intout out, intfinal final,
	   void *custom_data);
#endif
