/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  File: model/lib/grid/include/csa.h
 *  
 *  Description: A header file for csa's public functions
 *
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: csa.h 5839 2018-06-28 03:38:07Z riz008 $
 *
 */

#if !defined(_CSA_H)
#define _CSA_H

#if !defined(_CSA_STRUCT)
#define _CSA_STRUCT
typedef struct csa csa;
#endif

#include "config.h"
#include "ems.h"
#include "svd.h"
#include "minell.h"

csa* csa_create(void);
void csa_destroy(csa* a);
void csa_addpoints(csa* a, int n, point points[]);
void csa_addstd(csa* a, int n, double std[]);
void csa_calculatespline(csa* a);
void csa_approximatepoint(csa* a, point* p);
void csa_approximatepoints(csa* a, int n, point* points);
void csa_setnpmin(csa* a, int npmin);
void csa_setnpmax(csa* a, int npmax);
void csa_setk(csa* a, int k);
void csa_setnppc(csa* a, int nppc);
void csa_refindprimarycoeffs(csa* a, point *p);

#endif
