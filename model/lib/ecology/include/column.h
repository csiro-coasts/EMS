/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/column.h
 *  
 *  Description:
 *  Column structure -- header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: column.h 6770 2021-06-04 01:49:23Z riz008 $
 *
 */

#if !defined(_COLUMN_H)

#include "declarations.h"
#include "constants.h"

struct column {
    ecology* e;                 /* ecological model */
    void* model;                /* host model */

    int b;                      /* column index */
    int topk_wc;
    int botk_wc;
    int topk_sed;
    int botk_sed;

    int ncells;                 /* number of cells */
    int n_wc;                   /* number of (water + epibenthic) cells */
    int n_sed;                  /* number of (sediment + epibenthic) cells */

    /*
     * cell array [0..n_cells-1]; the cells are listed upside down 
     */
    cell** cells;

    double** tr_wc;             /* &wctracers[b] */
    double** tr_sed;            /* &sedtracers[b] */
    double** epivar;            /* &epitracers[b] */
    double* dz_wc;
    double* dz_sed;
    double* porosity;

    /* Column based view of tracers */
    double **y;
    double *zc;
    double *dz;

    double *y_epi;
    double *y_sed0;
  
    /*
     * process communication: common column variables 
     */
    int ncv;       /* number of common variables */
    int *n_ncv;    /* length of each common variable */
    double** cv;
};

/* This structure has been created for passing flux function parameters via
 * one variable.
 */
struct intargs {
    double t;
    double* y;
    double* y1;
    void* media;
};

/** Column constructor.
 * @param e Pointer to ecology structure.
 * @param b Column (box) index.
 * @return Pointer to the column structure.
 */
column* column_create(ecology* e, int i);

/** Column destructor.
 * @param col Pointer to column
 */
void column_destroy(column* c);

/** Performs ecology step for a column.
 * @param col Pointer to column
 * @return 0 for fail, 1 for success
 */
int column_step(column* c);

#define _COLUMN_H
#endif
