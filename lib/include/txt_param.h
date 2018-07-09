/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/txt_param.h
 *
 *  \brief Include file for parameter routines
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: txt_param.h 5834 2018-06-27 00:55:37Z riz008 $
 */

void param_fatal(int flag);
long dparam(FILE * fp, char *label, double *p, double min, double max);
long fparam(FILE * fp, char *label, float *p, double min, double max);
long lparam(FILE * fp, char *label, long *p, long min, long max);
long iparam(FILE * fp, char *label, int *p, int min, int max);
long sparam(FILE * fp, char *label, short *p, int min, int max);
long cparam(FILE * fp, char *label, char *p, char *chset);
long strparam(FILE * fp, char *label, char *p, long maxlen);
long dparamarray(FILE * fp, char *label, double **p, double min,
                 double max, long *size);
