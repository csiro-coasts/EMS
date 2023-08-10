/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/render.h
 *  
 *  Description:
 *  Structures for render management
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: render.h 5873 2022-03-02 07:23:48Z her127 $
 *
 */

#ifndef _RENDER_H
#define _RENDER_H

#define ECO_MAXNUMARGS 400

#include "ems.h"
#if defined(HAVE_ECOLOGY_MODULE)
#include "einterface.h"
#include "eprocess.h"
#include "ecology_internal.h"
#include "parameter_info.h"
#endif

typedef struct {
} null_ecosed_t;

/* Hydrodynamic global management structure */
typedef struct {
  char *name;
  char *desc;
  parameters_t *param;
  open_bdrys_t *obc;
} hydro_global_t;

/* Sediment global management structure */
typedef struct {
  char *name;
  char *desc;
  char *procname;
#if defined(HAVE_SEDIMENT_MODULE)
  sed_params_t *param;
#else
  null_ecosed_t *param;
#endif
} sed_global_param_t;

/* Ecology global management structure */
typedef struct {
  char   *name;
  char *desc;
  char *trname;
  tracer_info_t *trlist;
  int    ntr;
  char *econame;
#if defined(HAVE_ECOLOGY_MODULE)
  trinfo_priv_eco_t  *eco_priv;
#else
  null_ecosed_t *eco_priv;
#endif
  int    nep;
  char *sedname;
#if defined(HAVE_SEDIMENT_MODULE)
  trinfo_priv_sed_t  *sed_priv;
#else
  null_ecosed_t *sed_priv;
#endif
  int    nsp;
} eco_global_t;

typedef struct {
  char *name;
  char *desc;
  char *wcname;
  int nwc_proc;
  char *sedname;
  int nsed_proc;
  char *epiname;
  int nepi_proc;
} eco_global_proc_t;

typedef struct {
  int index;
  char* name;
  char* desc;
  char* sym;
  char* units;
  double value[100];
  int num_values; /* size of this parameter */
  double stderr;
  char* stringvalue;
  char* ref;
} parameter_glob_t;

typedef struct {
  char *name;
  char *desc;
  char *procname;
  parameter_glob_t *param;
  int nparam;
} eco_global_param_t;

#endif
