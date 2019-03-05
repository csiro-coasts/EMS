/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd/ecology/ecology_tracer_defaults.h
 *  
 *  Description:
 *  Tracer defaults header file
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology_tracer_defaults.h 6064 2019-02-08 04:07:50Z her127 $
 *
 */

#if !defined(_ECOLOGY_TRACER_DEFAULTS_H)
#include "ems.h"
/*
 * Ecology tracers are defined using the following fields 
 */
typedef struct {
  char *name;
  char *units;
  double fill_wc;
  double fill_sed;
  int diagn;
  int advect;
  int diffuse;
  int type;
  int dissol;
  int partic;
  int inwc;
  int insed;
  int flag;
  int obc;
  int ex_type;
} eco_def_t;

/*
 * Ecology particulate tracers are defined using the following fields
 */
typedef struct {
  char *name;
  double psize;
  double b_dens;
  double i_conc;
  double f_conc;
  double svel;
} eco_def_partic_t;

/*
 * Ecology tracerstat tracer definition
 */
typedef struct {
  char   *name;
  char   *long_name;
  char   *trname;
  char   *trstat;
  char   *trdt;
} eco_def_trstat_t;

/* Ecology tracer attribute private dtat structure */
typedef struct {
  int type;              /* Private data type */
  char name[MAXSTRLEN];  /* Name of the default list */
  int obc;               /* Open boundary condition */
  int flag;              /* General purpose flag */
} trinfo_priv_eco_t;

trinfo_priv_eco_t *get_private_data_eco(tracer_info_t *tr);
void init_tracer_atts_eco(tracer_info_t *tracer, char *trname, 
			  eco_def_t *eco_def, const char *eco_def_name,
			  eco_def_partic_t *eco_def_partic);

void eco_defaults_std(tracer_info_t *tracer, char *trname, ecology *e);
void eco_defaults_est(tracer_info_t *tracer, char *trname, ecology *e);

#endif
