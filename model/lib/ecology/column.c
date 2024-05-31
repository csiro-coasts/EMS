/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/column.c
 *  
 *  Description:
 *  Column structure -- implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: column.c 7539 2024-05-03 05:45:11Z riz008 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "externallibs.h"
#include "ecology_internal.h"
#include "einterface.h"
#include "stringtable.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "constants.h"
#include "utils.h"





/** Column constructor.
 * @param e Pointer to ecology structure.
 * @param b Column (box) index.
 * @return Pointer to the column structure.
 */
column* column_create(ecology* e, int b)
{
  column* col = calloc(1, sizeof(column));

  /*
   * k range, empty layers included 
   */
  int topk_wc = einterface_getwctopk(e->model, b);
  int botk_wc = einterface_getwcbotk(e->model, b);
  int topk_sed = einterface_getsedtopk(e->model, b);
  int botk_sed = einterface_getsedbotk(e->model, b);
  int n_wc = abs(topk_wc - botk_wc) + 1;
  int n_sed = abs(topk_sed - botk_sed) + 1;

  /*
   * k direction flags 
   */
  int up_wc = botk_wc < topk_wc;
  int up_sed = botk_sed < topk_sed;

  /*
   * bottom water / top sediment cell indices 
   */
  int ii = 0;
  int i, k, j;

  col->e = e;
  col->model = e->model;
  col->b = b;

  col->tr_wc = einterface_getwctracers(e->model, b);
  col->tr_sed = einterface_getsedtracers(e->model, b);
  col->epivar = einterface_getepivars(e->model, b);
  col->dz_wc = einterface_getwccellthicknesses(e->model, b);
  col->dz_sed = einterface_getsedcellthicknesses(e->model, b);
  col->porosity = einterface_getporosity(e->model, b);

  col->ncells = 0;
  col->n_wc = 0;
  col->n_sed = 0;

  col->cells = calloc(n_wc + n_sed, sizeof(void*));

  /*
   * What is attempted here is to automatically exclude empty water or
   * sediment cells. (By "empty" cells we count those with thickness less
   * than a specified value.) Therefore, first, closest to bottom
   * non-empty water and sediment cells are found and an epibenthic cell
   * is formed. Then, the remaining non-empty cells are processed. 
   */
  col->topk_wc = -1;
  col->botk_wc = -1;
  col->topk_sed = -1;
  col->botk_sed = -1;

  /*
   * find non-empty bottom water cell 
   */
  for (i = 0, k = botk_wc; i < n_wc; ++i, (up_wc) ? ++k : --k)
    if (col->dz_wc[k] >= THICKNESS_WC_MIN) {
      col->botk_wc = k;
      break;
    }
  /*
   * find non-empty top water cell 
   */
  for (i = 0, k = topk_wc; i < n_wc; ++i, (up_wc) ? --k : ++k)
    if (col->dz_wc[k] >= THICKNESS_WC_MIN) {
      col->topk_wc = k;
      break;
    }
  /*
   * find non-empty bottom sediment cell 
   */
  for (i = 0, k = botk_sed; i < n_sed; ++i, (up_sed) ? ++k : --k)
    if (col->dz_sed[k] >= THICKNESS_SED_MIN) {
      col->botk_sed = k;
      break;
    }

  /*
   * find non-empty top sediment cell 
   */
  for (i = 0, k = topk_sed; i < n_sed; ++i, (up_sed) ? --k : ++k)
    if (col->dz_sed[k] >= THICKNESS_SED_MIN) {
      col->topk_sed = k;
      break;
    }

  /*
   * create water cells 
   */
  if (e->npr[PT_WC] > 0) {
    if (up_wc) {
      for (k = col->topk_wc; k > col->botk_wc; --k)
	if (col->dz_wc[k] >= THICKNESS_WC_MIN) {
	  col->cells[ii] = cell_create(col, CT_WC, ii, k, -1);
	  ++ii;
	}
    } else {
      for (k = col->topk_wc; k < col->botk_wc; ++k)
	if (col->dz_wc[k] >= THICKNESS_WC_MIN) {
	  col->cells[ii] = cell_create(col, CT_WC, ii, k, -1);
	  ++ii;
	}
    }
  }

  /*
   * create epibenthic cell 
   */
  if (col->topk_wc == -1) {
    if (e->mandatory_water)
      e->quitfn("ecology: error: column %d: no non-empty water cells\n", col->b);
    else if (col->topk_sed >= 0) {
      if (e->npr[PT_SED] > 0) {
	col->cells[ii] = cell_create(col, CT_SED, ii, -1, col->topk_sed);
	++ii;
            
	/* TMP
	 */
	emstag(LTRACE,"eco:column:column_create","Have SED process, created water cell at: %d",col->b);
              
      }else
	emstag(LTRACE,"eco:column:column_create","Have NO EPI process, not cerated WC at: %d",col->b);
              
    } else
      e->quitfn("ecology: error: column %d: no non-empty cells\n", col->b);
  } else if (col->topk_sed == -1) {
    if (e->mandatory_sediment)
      e->quitfn("ecology: error: column %d: no non-empty sediment cells\n", col->b);
    else if (col->topk_wc >= 0) {
      if (e->npr[PT_WC] > 0) {
	col->cells[ii] = cell_create(col, CT_WC, ii, col->botk_wc, -1);
	++ii;
	/* TMP
	 */
	emstag(LTRACE,"eco:column:column_create","Have WC process, created water cell at: %d",col->b);
              
      }else
	emstag(LTRACE,"eco:column:column_create","Have NO WC process, not cerated WC at: %d",col->b);
          
    } else
      e->quitfn("ecology: error: column %d: no non-empty cells\n", col->b);
  } else {
    /*
     * Create epibenthic cell only if there are processes of PT_EPI
     * type; otherwhile create a water and a sediment cells. In the
     * latter case, there should not be any difference in the model
     * performance, whether there is an epibenthic cell or a
     * water+sediment cells. It just seems a bit artificial to me create
     * such a structure if there is no need for it. 
     */
    if (e->npr[PT_EPI] > 0) {
      col->cells[ii] = cell_create(col, CT_EPI, ii, col->botk_wc, col->topk_sed);
      ++ii;
    } else {
      if (e->npr[PT_WC] > 0) {
	col->cells[ii] = cell_create(col, CT_WC, ii, col->botk_wc, -1);
	++ii;
	emstag(LMETRIC,"eco:column:column_create","Cerated WC dummy cells at: %d",col->b);
    
      }
      if (e->npr[PT_SED] > 0) {
	col->cells[ii] = cell_create(col, CT_SED, ii, -1, col->topk_sed);
	++ii;
	emstag(LMETRIC,"eco:column:column_create","Cerated SED dummy cells at: %d",col->b);
      }
    }
  }

  /*
   * create sediment cells 
   */
  if (e->npr[PT_SED] > 0) {
    if (up_sed) {
      for (k = col->topk_sed - 1; k >= col->botk_sed; --k)
	if (col->dz_sed[k] >= THICKNESS_SED_MIN) {
	  col->cells[ii] = cell_create(col, CT_SED, ii, -1, k);
	  ++ii;
	}
    } else {
      for (k = col->topk_sed + 1; k <= botk_sed; ++k)
	if (col->dz_sed[k] >= THICKNESS_SED_MIN) {
	  col->cells[ii] = cell_create(col, CT_SED, ii, -1, k);
	  ++ii;
	}
    }
  }

  /* Column based tracers */
  col->y = d_alloc_2d(col->ncells, e->ntr);
  col->y_sed0 = d_alloc_1d(e->ntr);
  if (e->nepi) col->y_epi  = d_alloc_1d(e->nepi);
  col->zc = d_alloc_1d(n_wc);
  col->dz = d_alloc_1d(n_wc);
  i = 0;
  if (up_wc) {
    for (k = col->topk_wc; k >= col->botk_wc; --k) {
      col->zc[i] = einterface_getcellz(e->model, b, k);
      col->dz[i] = col->dz_wc[k];
      i++;
    }
  } else {
    for (k = col->topk_wc; k <= col->botk_wc; ++k) {
      col->zc[i] = einterface_getcellz(e->model, b, k);
      col->dz[i] = col->dz_wc[k];
      i++;
    }
  }
  
  col->ncv = e->cv_column->n;
  col->n_ncv = NULL;
  col->cv = NULL;
  if (col->ncv > 0) {
    /* Allocate column variable pointer array */
    col->cv = malloc(col->ncv * sizeof(double*));
    /* Allocate and initialise size of each column variables */
    col->n_ncv = malloc(col->ncv * sizeof(int));
    for (i=0; i<col->ncv; i++) {
      col->n_ncv[i] = stringtable_entry_get_n(e->cv_column, i);
      col->cv[i]    = malloc(col->n_ncv[i] * sizeof(double));
      for (j= 0; j < col->n_ncv[i]; ++j)
	col->cv[i][j] = NaN;  /* must be initialized by some process */
    }
  }

  return col;
}

/** Column destructor.
 * @param col Pointer to column
 */
void column_destroy(column* col)
{
  int i;

  if (col == NULL)
    return;

  /*
   * this includes writing back the new tracer values to the host model 
   */
  for (i = 0; i < col->ncells; ++i)
    cell_destroy(col->cells[i]);

  if (col->dz_wc != NULL)
    free(col->dz_wc);
  if (col->dz_sed != NULL)
    free(col->dz_sed);
  if (col->e->ntr > 0) {
    free(col->tr_wc);
    free(col->tr_sed);
  }
  if (col->e->nepi > 0)
    free(col->epivar);
  free(col->cells);
  if (col->ncv > 0) {
    /* free each pointer first */
    for(i=0; i<col->ncv; i++)
      free(col->cv[i]);
    free(col->n_ncv);
    free(col->cv);
  }
  if (col->y != NULL)
    d_free_2d(col->y);
  if (col->zc != NULL)
    d_free_1d(col->zc);
  if (col->dz != NULL)
    d_free_1d(col->dz);
  if (col->y_epi != NULL)
    d_free_1d(col->y_epi);
  if (col->y_sed0 != NULL)
    d_free_1d(col->y_sed0);
  
  free(col);
}

/** Fills local state variable values from the host model.
 * The values of diagnostic variables are set to 0.
 * @param c Pointer to cell
 * @param col Pointer to column
 */
static void cell_fill(cell* c, column* col)
{/*UR added minimum value for tracers defined in constants.h  */
  ecology* e = c->e;
  int* tdiagn = e->tracerdiagn;
  int* tflags = e->tracerflags;

  if (c->type == CT_WC) {
    double** tr_wc = &col->tr_wc[c->k_wc * e->ntr];
    double* y_wc = c->y;
    int i;

    for (i = 0; i < c->nvar; ++i)
      {
	if (tflags[i] & ECO_NORESET) {
	  /* pass through as-is */
	  y_wc[i] = *tr_wc[i];
	} else {
	  y_wc[i] = (!tdiagn[i]) ? *tr_wc[i] : 0.0;
	  if(y_wc[i] < TRACERVALUE_EPS)
	    y_wc[i] =TRACERVALUE_EPS;
	}
      }
  } else if (c->type == CT_SED) {
    double** tr_sed = &col->tr_sed[c->k_sed * e->ntr];
    double* y_sed = c->y;
    int i;

    for (i = 0; i < c->nvar; ++i)
      {
	if (tflags[i] & ECO_NORESET) {
	  /* pass through as-is */
	  y_sed[i] = *tr_sed[i];
	} else {
	  y_sed[i] = (!tdiagn[i]) ? *tr_sed[i] : 0.0;
	  if(y_sed[i] < TRACERVALUE_EPS)
	    y_sed[i] =TRACERVALUE_EPS;
	}
      }
  } else {                    /* CT_EPI */
    double** tr_wc = &col->tr_wc[c->k_wc * e->ntr];
    double** tr_sed = &col->tr_sed[c->k_sed * e->ntr];
    double* y_wc = c->y;
    double* y_sed = &c->y[e->ntr];
    double* y_epi = &c->y[e->ntr * 2];

    int* ediagn = e->epidiagn;
    int i;

    for (i = 0; i < e->ntr; ++i)
      {
	if (tflags[i] & ECO_NORESET) {
	  /* pass through as-is */
	  y_wc[i] = *tr_wc[i];
	} else {
	  y_wc[i] =  (!tdiagn[i]) ? *tr_wc[i] : 0.0;      
	  if(y_wc[i] < TRACERVALUE_EPS)
	    y_wc[i] =TRACERVALUE_EPS;
	}
      }
        
    for (i = 0; i < e->ntr; ++i)
      {
	if (tflags[i] & ECO_NORESET) {
	  /* pass through as-is */
	  y_sed[i] = *tr_sed[i];
	} else {
	  y_sed[i] = (!tdiagn[i]) ? *tr_sed[i] : 0.0;
	  if(y_sed[i] < TRACERVALUE_EPS)
	    y_sed[i] =TRACERVALUE_EPS;
	}
      }
    for (i = 0; i < e->nepi; ++i)
      {
	y_epi[i] = (!ediagn[i]) ? *col->epivar[i] : 0.0;
	if(y_epi[i] < TRACERVALUE_EPS)
	  y_epi[i] =TRACERVALUE_EPS;
      }
    /* Tack on the multi sed layers at the very end */
    if (e->use_multi_sed) {
      double* y_multi_sed  = &c->y[e->ntr * 2 + e->nepi];
      int k, ii=0;
      // Start at the second layer
      for (k = col->topk_sed-1; k >= col->botk_sed; --k) {
	double** tr_multi_sed = &col->tr_sed[k * e->ntr];
	for (i = 0; i < e->ntr; ++i, ++ii) {
	  if (tflags[i] & ECO_NORESET) {
	    /* pass through as-is */
	    y_multi_sed[ii] = *tr_multi_sed[i];
	  } else {
	    y_multi_sed[ii] = (!tdiagn[i]) ? *tr_multi_sed[i] : 0.0;
	    if (y_multi_sed[ii] < TRACERVALUE_EPS)
	      y_multi_sed[ii] =TRACERVALUE_EPS;
	  }
	}
      }
    }
  }
}

/**
 *
 */
static void column_run(column* col)
{
  ecology *e = col->e;
  int i;

  /* Loop through all the column processes */
  for (i = 0; i < e->npr[PT_COL]; ++i) {
    eprocess* p = e->processes[PT_COL][i];
    if (p->precalc == NULL)
      e_quit("column_run: only pre_calc is currently supported");
    p->precalc(p, col);
  }
}

static void column_writeback(column* col)
{
  ecology* e = col->e;
  double dt = e->dt;

  /*
   * Consists of a single stacked vector
   */
  double** y_wc   = col->y;
  double*  y_epi  = col->y_epi;
  double*  y_sed0 = col->y_sed0;
  int k, k_wc, n, ntr = e->ntr;

  for (n = 0; n < ntr; n++) {
    /* loop from the top */
    for (k = 0, k_wc = col->topk_wc; k_wc >= col->botk_wc; k_wc--, k++)
      *col->tr_wc[k_wc * ntr + n] = y_wc[n][k];

    // xxx what about fluxes that need to be divided by dt?
    
    /* Fill top of sediment */
    *col->tr_sed[col->topk_sed * ntr + n] = y_sed0[n];
  }
  
  /* Fill epi */
  for (n = 0; n < e->nepi; n++)
    *col->epivar[n] = y_epi[n];
}

/** 
 * Column based values
 */
static void column_fill(column* col)
{
  ecology* e = col->e;

  /*
   * Consists of a single stacked vector
   */
  double** y_wc   = col->y;
  double*  y_epi  = col->y_epi;
  double*  y_sed0 = col->y_sed0;
  int k, k_wc, n, ntr = e->ntr;

  for (n = 0; n < ntr; n++) {
    /* loop from the top */
    for (k = 0, k_wc = col->topk_wc; k_wc >= col->botk_wc; k_wc--, k++)
      y_wc[n][k] = *col->tr_wc[k_wc * ntr + n];

    /* Fill top of sediment */
    y_sed0[n] = *col->tr_sed[col->topk_sed * ntr + n];
  }
  
  /* Fill epi */
  for (n = 0; n < e->nepi; n++) 
    y_epi[n] = *col->epivar[n];
}

/** Performs ecology step for a column.
 * @param col Pointer to column
 * @return 0 for fail, 1 for success
 */
int column_step(column* col)
{
  int k;
  emslog(LMETRIC,"Starting Column-%d Step at %f \n",col->b,wall_time());
  /*
   * column_step -- just stepping all the cells 
   */
  for (k = 0; k < col->ncells; ++k) {
    cell* c = col->cells[k];
    assert(k == c->k);

    /*
     * column_prestep processes write directly to the state variable
     * storage in the host model, while cell processes read/write
     * from/to the local cell storage cell.y. This is why filling of
     * cell's tracer values must be done in between those two groups of
     * processes. The values of diagnostic variables are set to 0. 
     */
    cell_fill(c, col);

    /* pre integration stuff */
    cell_precalc(c); 
    /* integration */
    if (!cell_calc(c)) {
      /* Failed, bail out */
      return 0;
    }
    cell_postcalc(c); /* post integration stuff */
    /*
     * Writes calculated tracer concentrations back to the host
     * model. The values of "flux" diagnostic tracers are written back
     * after being divided by dt to result in average flux value over
     * the time step. 
     */
    cell_writeback(c);
  }
  /*UR added writeback to capture the postprocessing variables */
  /* now iterate again to write it back again after the entire column has bee calculated */
  for (k = 0; k < col->ncells; ++k) 
    {
      cell* c = col->cells[k];
      /*
       * Writes calculated tracer concentrations back to the host
       * model. The values of "flux" diagnostic tracers are written back
       * after being divided by dt to result in average flux value over
       * the time step. 
       */
      cell_writeback(c);
    }

  /* Run column based processes */
  if (col->e->npr[PT_COL] && col->topk_wc > -1) {
    column_fill(col);
    column_run(col);
    column_writeback(col);
  }
  
  return 1;
}
