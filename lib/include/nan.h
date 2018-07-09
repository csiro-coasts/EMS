/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/nan.h
 *
 *  \brief EMS definition of NaN
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: nan.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#if !defined(_NAN_H)
#define _NAN_H

#if defined(__GNUC__)

#if defined(__ia64__)

static const long long lNaN = ((unsigned long long) 1 << 63) - 1;
#define NaN (*(double*)&lNaN)
#else
  static const double NaN = 0.0 / 0.0;
#endif

#elif defined(_WIN32) &&  !defined(__MINGW32__)

static unsigned _int64 lNaN = ((unsigned _int64) 1 << 63) - 1;

#define NaN (*(double*)&lNaN)

#else

static const long long lNaN = ((unsigned long long) 1 << 63) - 1;

#define NaN (*(double*)&lNaN)

#endif

#endif

