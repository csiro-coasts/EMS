/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/cell.h
 *  
 *  Description:
 *  Cell structure -- header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: cell.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_CELL_H)

#include "declarations.h"

typedef enum {
    CT_NONE = 0,
    CT_WC = 1,
    CT_EPI = 2,
    CT_SED = 3
} CELLTYPE;

#define WCCHILD 0
#define SEDCHILD 1

struct cell {
    ecology* e;
    column* col;

    CELLTYPE type;
    int b;                      /* column (horizontal) index in the
                                 * ecological & host models */
    int k;                      /* vertical index in the ecological model */
    int k_wc;                   /* vertical index in the host model */
    int k_sed;                  /* vertical index in the host model */
    /*
     * cell thickness 
     */
    double dz_wc;
    double dz_sed;
    double *dz_multi_sed;
    double porosity;
    double *porosity_multi_sed;

    int nvar;                   /* number of integrated variables (tracers) */
    double* y;                  /* integration array */

    /*
     * process communication: common variables 
     */
    int ncv;
    double* cv;

    /*
     * total number of the precalc, calc and postcalc function calls 
     */
    int nsubstep;

    intargs* ia;                /* to pass flux values from integrator for
                                 * checks etc.; points to current intargs
                                 * structure if integration is under way;
                                 * NULL otherwise */

    /*
     * A CT_EPI cell must have no parent but (normally) 2 children. A CT_WC
     * or CT_SED cell must have no children; it must have a parent if and
     * only if it belongs to an epibenthic cell. 
     */
    cell* wcchild;
    cell* sedchild;
    cell* parent;

    double* h0;
};

/** Cell constructor.
 * @param col Pointer to column the cell is in
 * @param ct Cell type
 * @param k Vertical index in the column
 * @param k_wc Vertical watercolumn index in the host model
 * @param k_sed Vertical sediment index in the host model
 * @return Pointer to cell
 */
cell* cell_create(column* col, CELLTYPE ct, int k, int k_wc, int k_sed);

/** Cell destructor.
 * @param c Pointer to cell
 */
void cell_destroy(cell* c);

/** Runs precalc procedures for the cell (those to be performed before the
 * integration loop).
 * @param c Pointer to cell
 */
void cell_precalc(cell* c);

/** Integration loop for a cell.
 * @param c Pointer to cell
 * @return 0 for fail, 1 for success
 */
int cell_calc(cell* c);

/** Runs postcalc procedures for the cell (those to be performed after the
 * integration loop).
 * @param c Pointer to cell
 */
void cell_postcalc(cell* c);

/** Writes calculated tracer concentrations back to the host model.
 * The values of "flux" diagnostic tracers are written back after
 * being divided by dt to result in average flux value over the time
 * step.
 * @param c Pointer to cell
 */
void cell_writeback(cell* c);

#define _CELL_H
#endif
