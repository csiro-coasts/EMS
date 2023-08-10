/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/tracer.h
 *
 *  \brief Include file for tracer routines 
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: tracer.h 7022 2022-02-28 01:10:27Z her127 $
 */

#ifndef _TRACER_H
#define _TRACER_H

#include <string.h>
#include "ems.h"

/*UR indicate to BoxModel to use the ems trace.h */
#define _DOTRACER_FIXED

#define MAXNUMTRACERS 1000

/* Check tracer positiveness */
#define TRACERTEST 1

/* Use fully implicit code for vertical diffusion of tracers */
#define TR_VD_FULLIMPLICIT

/* Use explicit code for vertical diffusion of tracers */
/* #define TR_VD_EXPLICIT */
/* Following is a typical value of critical shear stress for fine sediments 
 * consolidated for about 24 hours - see Teisson J Hydraul. Res. 29 1992.
 */
#define CRITICAL_STRESS_DEF 0.2

typedef void* (*private_data_copy_t) (void *src);

enum {
  TR_PRV_DATA_SED = 0,
  TR_PRV_DATA_ECO,
  TR_NPRV_DATA
};

typedef struct {
/* header */
  int n;                        /* index */
  int m;                        /* parameter file map index */
  char name[MAXSTRLEN];         /* name (netcdf name, etc) */
  char long_name[MAXSTRLEN];    /* long name (netcdf name, etc) */
  char std_name[MAXSTRLEN];     /* standard name (CF name, etc) */
  char units[MAXSTRLEN];        /* units */
  char groupkey[MAXSTRLEN];     /* Tag for autotracers */
  int type;                     /* type of tracer */ 
/* ranges and default values */
  double fill_value_wc;         /* watercolumn fill value */
  double valid_range_wc[2];     /* valid range in watercolumn */
  double fill_value_sed;        /* sediment fill value */
  double valid_range_sed[2];    /* valid range in sediment */
  /*UR added 5/2008 */
  double missing_value_sed;     /* sediment missing value, set to missing value or NaN if not defined  */
  double missing_value;         /* missing value, set to NaN if not defined  */



/* flags */
  char tag[MAXSTRLEN];          /* general purpose tag */
  int flag;                     /* general purpose flag */
  int diagn;                    /* flag non-zero if diagnostic; default 0.
                                   * if non-zero, 1 if "flux", other if
                                   "value". */
  int inwc;                     /* flag 1 if exists in water; default 1 */
  int insed;                    /* flag 1 if exists in sediments; default
                                   1 */
  int dissol;                   /* flag 1 if dissolved; default !partic *
                                   (if present) or 1 */
  int partic;                   /* flag 1 if particulate; default !dissol
                                   * (if present) or 0 */
  int advect;                   /* flag 1 if tracer to be advected;
                                   default 1 */
  int diffuse;                  /* flag 1 if tracer to be diffused;
                                   default 1 */
  int increment;                /* Flag for incrementing state variables */
/* candidate for exclusion */
  char decay[MAXSTRLEN];        /* decay rate */
/* candidates for exclusion */
  double ***Kz;                 /* pointer to vertical diffusivity values
                                   to * use */
  double Kz_scale;              /* scale factor for vertical diffusivity * 
                                   values */

  timeseries_t *ts;             /* time-series file */
  int var_id;                   /* variable id into TS file. */
  char tracerstat[MAXSTRLEN];   /* Tracer statistic attribute */
  char trstat_tag[MAXSTRLEN];   /* Additional tracerstat info */
  char data[MAXSTRLEN];         /* Initialisation data */

  char relax_file[MAXSTRLEN];   /* Relax input file */
  char relax_dt[MAXSTRLEN];     /* Relax time interval */
  char r_rate[MAXSTRLEN];       /* relaxation rate (s-1) or * 1 / (time
                                   constant) */
  double relax_rate;
  double relax_dum;             /* Dummy for relaxation */
  int tctype;                   /* Relaxation time constant type */

  char reset_file[MAXSTRLEN];   /* Reset input file */
  char reset_dt[MAXSTRLEN];     /* Reset time interval */

  char i_rule[MAXSTRLEN];       /* Interpolation rule */
  double scale;                 /* Scaling factor */
  char vector_name[MAXSTRLEN];  /* Name for vector tracers */
  char vector_components[MAXSTRLEN];  /* Components for vector tracers */

  /* Private data */
  void *private_data[TR_NPRV_DATA];
  private_data_copy_t private_data_copy[TR_NPRV_DATA];

  char reset_interp[MAXSTRLEN];

} tracer_info_t;

typedef struct {
  char attr[MAXLINELEN];
  char attrval[MAXLINELEN];
} tracer_att_t;


/* prototypes */
void tracer_read(FILE * fp, char *prefix, int type, errfn quitfn,
                 errfn warnfn, errfn emptyfn, int *ntr,
                 tracer_info_t *tracers[]);
void tracer_read_old(FILE * fp, char *prefix, char *keyword, errfn quitfn,
                 errfn warnfn, errfn emptyfn, int *ntr,
                 tracer_info_t *tracers[]);
void tracer_read_nc(int fid, int ndims, int dimids[], tracer_att_t *attr,
                    int *ntr, tracer_info_t *tracers[]);
void tracer_write_nc(int fid, int ntr, tracer_info_t tracers[], int nattr,
                     tracer_att_t attr[]);
void tracer_copy(tracer_info_t *tr_cpy, tracer_info_t *tr_in);
int tracer_find_index(const char *name, int ntr, tracer_info_t tracers[]);

/*UR ADDED from merger 
 */
void tracer_swap(tracer_info_t tinfo[], int i, int j);
void tracer_clear(int *ntr, tracer_info_t *tracers[]);
int tracer_read_attribute(FILE * fp, const char *prefix,const char *keyname, const char* trname, int i, const char *tag, errfn fn, char *buffer);

/*UR*/
double Ultimate(double cr, double Fx, double F, double Fm1, double Fm2,
                double Fp1);
void weights_t(double *wgt, double a, double b, double c);

#endif
