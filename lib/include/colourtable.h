/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/colourtable.h
 *
 *  \brief Header for colour_table_t data structures
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: colourtable.h 5834 2018-06-27 00:55:37Z riz008 $
 */

typedef struct {
  int n;
  double *v;
  double *r;
  double *g;
  double *b;
} colour_table_t;

/* Prototypes */
colour_table_t *ct_read(char *fname);
void ct_get_RGB(double v, colour_table_t ct,
                double *r, double *g, double *b);
