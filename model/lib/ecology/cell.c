/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/cell.c
 *  
 *  Description:
 *  Cell structure -- implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: cell.c 7193 2022-09-13 06:14:06Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#if defined(__sparc)
#include <ieeefp.h>
#endif
#include <assert.h>
#include "ecology_internal.h"
// #include "einterface.h"
#include "parameter_info.h"
#include "stringtable.h"
#include "eprocess.h"
#include "cell.h"
#include "column.h"
#include "utils.h"
#include "externallibs.h"

#define NDPE_MAX 100
#define SMALL_EPS 1e-5

void ginterface_get_ij(void* model, int col, int *ij);

/** Cell constructor.
 * @param parent Parent cell
 * @param ct Created cell type
 * @return Pointer to cell
 */
static cell* cell_createchild(cell* parent, CELLTYPE ct)
{
  ecology* e = parent->e;
  cell* c = malloc(sizeof(cell));
  int i;

  assert(parent->type == CT_EPI);

  c->e = parent->e;
  c->col = parent->col;
  c->type = ct;
  c->b = parent->b;
  c->k = -1;
  c->nvar = e->ntr;
  c->ncv = e->cv_cell->n;
  c->wcchild = NULL;
  c->sedchild = NULL;
  c->parent = parent;
  c->h0 = NULL;
  c->dz_multi_sed = NULL;
  c->porosity_multi_sed = NULL;

  if (ct == CT_WC) {
    c->k_wc = parent->k_wc;
    c->k_sed = -1;
    c->dz_wc = parent->dz_wc;
    c->dz_sed = NaN;
    c->porosity = 1.0;      /* to back processes common for wc and sed */
    c->y = parent->y;
    c->cv = parent->cv;
  } else if (ct == CT_SED) {
    c->k_wc = -1;
    c->k_sed = parent->k_sed;
    c->dz_wc = NaN;
    c->dz_sed = parent->dz_sed;
    c->porosity = parent->porosity;
    c->y = &parent->y[e->ntr];      /* offset to point to sed tracers */
    if (c->ncv > 0)
      c->cv = malloc(c->ncv * sizeof(double));
    for (i = 0; i < c->ncv; ++i)
      c->cv[i] = NaN;     /* must be initialized by some process */

  } else
    e->quitfn("ecology: programming error\n");

  c->nsubstep = -1 /* do not count for a child cell */ ;
  c->ia = NULL;

  return c;
}

/** Cell constructor.
 *
 * @param col Pointer to column the cell is in
 * @param ct Cell type
 * @param k Vertical index in the column
 * @param k_wc Vertical watercolumn index in the host model
 * @param k_sed Vertical sediment index in the host model
 * @return Pointer to cell
 */
cell* cell_create(column* col, CELLTYPE ct, int k, int k_wc, int k_sed)
{
  ecology* e = col->e;
  cell* c = malloc(sizeof(cell));
  int i;

  c->e = e;
  c->col = col;
  c->type = ct;
  c->b = col->b;
  c->k = k;
  c->k_wc = k_wc;
  c->k_sed = k_sed;
  c->wcchild = NULL;
  c->sedchild = NULL;
  c->parent = NULL;
  c->dz_multi_sed = NULL;
  c->porosity_multi_sed = NULL;

  c->ncv = e->cv_cell->n;
  if (c->ncv > 0)
    c->cv = malloc(c->ncv * sizeof(double));
  for (i = 0; i < c->ncv; ++i)
    c->cv[i] = NaN;         /* must be initialized by some process */

  if (ct == CT_WC) {
    c->dz_wc = col->dz_wc[k_wc];
    c->dz_sed = NaN;
    c->porosity = 1.0;      /* to back processes common for wc and sed */
    col->n_wc++;
    c->nvar = e->ntr;
    c->y = malloc(c->nvar * sizeof(double));
    c->h0 = &e->h0[(e->nwclayers + e->nsedlayers) * c->b + c->k_wc];
  } else if (ct == CT_SED) {
    c->dz_wc = NaN;
    c->dz_sed = col->dz_sed[k_sed];
    c->porosity = col->porosity[k_sed];
    col->n_sed++;
    c->nvar = e->ntr;
    c->y = malloc(c->nvar * sizeof(double));
    c->h0 = &e->h0[(e->nwclayers + e->nsedlayers) * c->b + e->nwclayers + c->k_sed];
  } else {                    /* CT_EPI */
    c->dz_wc = col->dz_wc[k_wc];
    c->dz_sed = col->dz_sed[k_sed];
    c->porosity = col->porosity[col->topk_sed];
    col->n_wc++;
    col->n_sed++;
    /* Set up extra layers for multi sediments, if needed */
    if (e->use_multi_sed) {
      int kk;
      c->nvar = e->ntr * (e->nsedlayers + 1) + e->nepi;
      c->dz_multi_sed = malloc((e->nsedlayers) * sizeof(double));
      c->porosity_multi_sed = malloc((e->nsedlayers) * sizeof(double));
      for (kk = col->topk_sed; kk >= col->botk_sed; --kk) {
	c->dz_multi_sed[kk] = col->dz_sed[kk];
	c->porosity_multi_sed[kk] = col->porosity[kk];
      }
    } else
      c->nvar = e->ntr * 2 + e->nepi;
    c->y = malloc(c->nvar * sizeof(double));

    c->wcchild = cell_createchild(c, CT_WC);
    c->sedchild = cell_createchild(c, CT_SED);
    c->h0 = &e->h0[(e->nwclayers + e->nsedlayers) * c->b + e->nwclayers - 1];
  }
  col->ncells++;

  c->nsubstep = 0;
  c->ia = NULL;

  return c;
}

/** Writes calculated tracer concentrations back to the host model.
 ** The values of "flux" diagnostic tracers are written back after
 ** being divided by dt to result in average flux value over the time
 ** step.
 *
 * @param c Pointer to cell
 */
void cell_writeback(cell* c)
{
  ecology* e;
  column* col;
  int* tdiagn;
  int* tflags;
  double dt;

  if (c == NULL)
    return;

  col = c->col;
  e = col->e;
  tdiagn = c->e->tracerdiagn;
  tflags = c->e->tracerflags;
  dt = e->dt;

  /*
   * writing back 
   */
  if (c->type == CT_WC) {
    double** tr_wc = &col->tr_wc[c->k_wc * e->ntr];
    double* y_wc = c->y;
    int i;

    for (i = 0; i < e->ntr; ++i)
      if (!(tflags[i] & ECO_NORESET))
	*tr_wc[i] = (tdiagn[i] != 1) ? y_wc[i] : y_wc[i] / dt;
  } else if (c->type == CT_SED) {
    double** tr_sed = &col->tr_sed[c->k_sed * e->ntr];
    double* y_sed = c->y;
    int i;
    
    for (i = 0; i < e->ntr; ++i)
      if (!(tflags[i] & ECO_NORESET))
	*tr_sed[i] = (tdiagn[i] != 1) ? y_sed[i] : y_sed[i] / dt;
  } else {
    double** tr_wc = &col->tr_wc[c->k_wc * e->ntr];
    double** tr_sed = &col->tr_sed[c->k_sed * e->ntr];
    double* y_wc = c->y;
    double* y_sed = &c->y[e->ntr];
    double* y_epi = &c->y[e->ntr * 2];
    int* ediagn = e->epidiagn;
    int i;

    for (i = 0; i < e->ntr; ++i)
      if (!(tflags[i] & ECO_NORESET))
	*tr_wc[i] = (tdiagn[i] != 1) ? y_wc[i] : y_wc[i] / dt;
    for (i = 0; i < e->ntr; ++i)
      if (!(tflags[i] & ECO_NORESET))
	*tr_sed[i] = (tdiagn[i] != 1) ? y_sed[i] : y_sed[i] / dt;
    for (i = 0; i < e->nepi; ++i)
      *col->epivar[i] = (ediagn[i] != 1) ? y_epi[i] : y_epi[i] / dt;

    /* Update the deeper sediment layers, if needed */
    if (e->use_multi_sed) {
      double* y_multi_sed  = &c->y[e->ntr * 2 + e->nepi];
      int k, ii=0;
      // Start at the second layer
      for (k = col->topk_sed-1; k >= col->botk_sed; --k) {
	double** tr_multi_sed = &col->tr_sed[k * e->ntr];
	for (i = 0; i < e->ntr; ++i, ++ii) {
	  *tr_multi_sed[i] = (tdiagn[i] != 1) ? 
	                        y_multi_sed[ii] : y_multi_sed[ii] / dt;
	}
      }
    }
  }
}

/** Cell destructor.
 *
 * @param c Pointer to cell
 */
void cell_destroy(cell* c)
{
  if (c == NULL)
    return;

  if (c->wcchild != NULL) {
    cell* child = c->wcchild;

    if (child->ncv > 0 && child->cv != c->cv)
      free(child->cv);
    free(child);
  }
  if (c->sedchild != NULL) {
    cell* child = c->sedchild;

    if (child->ncv > 0 && child->cv != c->cv)
      free(child->cv);
    free(child);
  }

  if (c->dz_multi_sed)
    free(c->dz_multi_sed);
  if (c->porosity_multi_sed)
    free(c->porosity_multi_sed);

  free(c->y);
  if (c->ncv > 0)
    free(c->cv);
  free(c);
}

static EPROCESSTYPE celltype2processtype(CELLTYPE ct)
{
  if (ct == CT_WC)
    return PT_WC;
  if (ct == CT_EPI)
    return PT_EPI;
  if (ct == CT_SED)
    return PT_SED;
  e_quit("error: ecology: programming error\n");
  return 0;
}

typedef struct {
  eprocess* p;
  cell* c;
} dpe;                          /* delayed process entry */

/** Check process output for NaN or infinite values.
 * @param eprocess Process
 * @param cell : containing cell
 */
static void check_nans(eprocess* p, cell* c)
{
  ecology* e = c->col->e;
  intargs* ia = c->ia;
  int ij[2] ;
  double* y = (ia == NULL) ? c->y : ia->y;
  int i, ii;
  einterface_get_ij(c->e->model,c->b,ij);
 
  if (c->parent != NULL)
    c = c->parent;

  if (p->type == PT_WC || p->type == PT_SED) {
    for (i = 0; i < e->ntr; ++i)
      {
	if (!finite(y[i]))
	  e->quitfn("ecology: error: \"%s\": NaN or inf detected for \"%s\" in %s cell, nstep = %d, nsubstep = %d, b = %d, k = %d <%u %u>\n", p->name, e->tracernames[i], (p->type == PT_WC) ? "water" : "sediment", e->nstep, c->nsubstep, c->col->b, (p->type == PT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
      }
    if (ia == NULL)
      return;
    y = ia->y1;
    for (i = 0; i < e->ntr; ++i)
      {
	if (!finite(y[i]))
	  e->quitfn("ecology: error: \"%s\": NaN or inf detected for flux of \"%s\" in %s cell, nstep = %d, nsubstep = %d, b = %d, k = %d <%u/%u>\n", p->name, e->tracernames[i], (p->type == PT_WC) ? "water" : "sediment", e->nstep, c->nsubstep, c->col->b, (p->type == PT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
      }
  } else if (p->type == PT_EPI) {
    for (i = 0; i < e->ntr; ++i)
      if (!finite(y[i]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for \"%s\" in water, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->tracernames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr; i < e->ntr; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for \"%s\" in sediment, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->tracernames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr * 2; i < e->nepi; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for \"%s\", nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->epinames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    if (ia == NULL)
      return;
    y = ia->y1;
    for (i = 0; i < e->ntr; ++i)
      if (!finite(y[i]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for flux of \"%s\" in water, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->tracernames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr; i < e->ntr; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for flux of \"%s\" in sediment, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->tracernames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr * 2; i < e->nepi; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: \"%s\": NaN or inf detected for flux of \"%s\"\n in epibenthos, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", p->name, e->epinames[i], e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
  }
}

/** Check process output for negative values. Non-diagnostic tracers only.
 * @param p Process
 * @param media Media: either column or cell
 */
static void check_negs(cell* c, char* state)
{
  ecology* e = c->col->e;
  intargs* ia = c->ia;
  int ij[2];
  double* y = (ia == NULL) ? c->y : ia->y;
  int i, ii;
  einterface_get_ij(c->e->model,c->b,ij);

  if (c->type == CT_WC || c->type == CT_SED) {
    for (i = 0; i < e->ntr; ++i)
      if ( (y[i] + SMALL_EPS) < 0.0 && !e->tracerdiagn[i])
	e->quitfn("ecology: error: negative value of %.3g detected for \"%s\" at %s - in %s cell, nstep = %d, nsubstep = %d, b = %d, k = %d <%u/%u>", y[i],e->tracernames[i], state,(c->type == CT_WC) ? "water" : "sediment", y[i], e->nstep, c->nsubstep, c->col->b, (c->type == CT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
  } else if (c->type == CT_EPI) {
    for (i = 0; i < e->ntr; ++i)
      if ( (y[i] + SMALL_EPS) < 0.0 && !e->tracerdiagn[i])
	e->quitfn("ecology: error: negative value of %.3g detected for \"%s\" at %s - in water, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", y[i], e->tracernames[i],  state,e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr; i < e->ntr; ++i, ++ii)
      if ( (y[ii] + SMALL_EPS) < 0.0 && !e->tracerdiagn[i])
	e->quitfn("ecology: error: negative value of %.3g detected for \"%s\" at %s - in sediment, epibenthic cell, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", y[ii], e->tracernames[i],  state,e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
    for (i = 0, ii = e->ntr * 2; i < e->nepi; ++i, ++ii)
      if ( (y[ii]  + SMALL_EPS)< 0.0 && !e->epidiagn[i])
	e->quitfn("ecology: error: negative value of %.3g detected for \"%s\" at %s, nstep = %d, nsubstep = %d, b = %d <%u/%u>\n", y[ii], e->epinames[i], state,e->nstep, c->nsubstep, c->col->b,ij[0],ij[1]);
  }
}

/** Checks all tracers for infinite or NaN values. Supposed to be run prior 
 * to running any processes.
 *
 * @param c cell
 */
static void check_input_for_nans(cell* c)
{
  ecology* e = c->col->e;
  intargs* ia = c->ia;
  int ij[2] ;
  double* y = (ia == NULL) ? c->y : ia->y;
  int i, ii;
    
  einterface_get_ij(c->e->model,c->b,ij);
    
  if (c->type == CT_WC || c->type == CT_SED) {
    for (i = 0; i < e->ntr; ++i)
      if (!finite(y[i]))
	e->quitfn("ecology: error: NaN or inf detected for \"%s\" in %s cell at entering ecology, nstep = %d, b = %d, k = %d <%u %u>", e->tracernames[i], (c->type == CT_WC) ? "water" : "sediment", e->nstep, c->col->b, (c->type == CT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
  } else if (c->type == CT_EPI) {
    for (i = 0; i < e->ntr; ++i)
      if (!finite(y[i]))
	e->quitfn("ecology: error: NaN or inf detected for \"%s\" in water cell at entering ecology, nstep = %d, b = %d, k = %d <%u %u>", e->tracernames[i], e->nstep, c->col->b, (c->type == CT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
    for (i = 0, ii = e->ntr; i < e->ntr; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: NaN or inf detected for \"%s\" in sediment cell, at entering ecology, nstep = %d, b = %d, k = %d <%u %u>", e->tracernames[i], e->nstep, c->col->b, (c->type == CT_WC) ? c->k_wc : c->k_sed,ij[0],ij[1]);
    for (i = 0, ii = e->ntr * 2; i < e->nepi; ++i, ++ii)
      if (!finite(y[ii]))
	e->quitfn("ecology: error: NaN or inf detected for \"%s\" in epibenthos at entering ecology, nstep = %d, b = %d <%u %u>", e->epinames[i], e->nstep, c->col->b,ij[0],ij[1]);
  }
}

/** Runs precalc procedures for the cell (those to be performed before the
 ** integration loop).
 *
 * @param c Pointer to cell
 */
void cell_precalc(cell* c)
{
  ecology* e = c->e;
  EPROCESSTYPE pt = celltype2processtype(c->type);
  int ndpe = 0;
  dpe dpes[NDPE_MAX];
  int i;  
  char buf[MAXSTRLEN];
		
  if (e->check_nans)
    check_input_for_nans(c);

  if (c->wcchild != NULL) {
    cell* child = c->wcchild;

    if (e->check_negs)
      check_negs(child,"pre- precalc-wcchild");
    for (i = 0; i < e->npr[PT_WC]; ++i) {
      eprocess* p = e->processes[PT_WC][i];
      sprintf(buf,"post- precalc-%s-wcchild", p->name);
      if (p->precalc != NULL) {
	if (p->delayed) {
	  dpes[ndpe].p = p;
	  dpes[ndpe].c = child;
	  ndpe++;
	} else {
	  p->precalc(p, child);
	  if (e->check_nans)
	    check_nans(p, child);
	  if (e->check_negs)
	    check_negs(child,buf);
	}
      }
    }
  }

  if (e->check_negs)
    check_negs(c,"pre-precalc");
  for (i = 0; i < e->npr[pt]; ++i) {
    eprocess* p = e->processes[pt][i];
    sprintf(buf, "post-precalc-%s",p->name);
    if (p->precalc != NULL) {
      if (p->delayed) {
	dpes[ndpe].p = p;
	dpes[ndpe].c = c;
	ndpe++;
      } else {
	p->precalc(p, c);
	if (e->check_nans)
	  check_nans(p, c);
	if (e->check_negs)
	  check_negs(c,buf);
      }
    }
  }

  if (c->sedchild != NULL) {
    cell* child = c->sedchild;

    if (e->check_negs)
      check_negs(child,"pre-precalc-sedchild");
    for (i = 0; i < e->npr[PT_SED]; ++i) {
      eprocess* p = e->processes[PT_SED][i];
      sprintf(buf, "post-precalc%s-sedchild",p->name);
      if (p->precalc != NULL) {
	if (p->delayed) {
	  dpes[ndpe].p = p;
	  dpes[ndpe].c = child;
	  ndpe++;
	} else {
	  p->precalc(p, child);
	  if (e->check_nans)
	    check_nans(p, child);
	  if (e->check_negs)
	    check_negs(child,buf);
	}
      }
    }
  }

  for (i = 0; i < ndpe; ++i) {
    eprocess* p = dpes[i].p;
    sprintf(buf, "post-delayed-precalc-%s",p->name);
    p->precalc(p, dpes[i].c);
    if (e->check_nans)
      check_nans(p, dpes[i].c);
    if (e->check_negs)
      check_negs(dpes[i].c,buf);
  }

  c->nsubstep++;
}

/** Flux calculation for a watercolumn cell.
 * @param t Time (input)
 * @param y Tracer concentrations (input)
 * @param y1 Fluxes (output)
 * @param p Pointer to cell
 */
static void wcfluxes(double t, double* y, double* y1, void* pp)
{
  intargs ia;
  cell* c = (cell*) pp;
  ecology* e = c->col->e;
  int nvar = c->nvar;
  int i;

  ia.t = t;
  ia.y = y;
  ia.y1 = y1;
  ia.media = pp;

  c->ia = &ia;

  for (i = 0; i < nvar; ++i)
    y1[i] = 0.0;

  for (i = 0; i < e->npr[PT_WC]; ++i) {
    eprocess* p = e->processes[PT_WC][i];
    if (p->calc != NULL)
      p->calc(p, &ia);
  }

  c->ia = NULL;
  c->nsubstep++;
}

/** Flux calculation for a sediment cell.
 * @param t Time (input)
 * @param y Tracer concentrations (input)
 * @param y1 Fluxes (output)
 * @param p Pointer to cell
 */
static void sedfluxes(double t, double* y, double* y1, void* pp)
{
  intargs ia;
  cell* c = (cell*) pp;
  ecology* e = c->col->e;
  int nvar = c->nvar;
  int i;

  ia.t = t;
  ia.y = y;
  ia.y1 = y1;
  ia.media = pp;

  c->ia = &ia;

  for (i = 0; i < nvar; ++i)
    y1[i] = 0.0;

  for (i = 0; i < e->npr[PT_SED]; ++i) {
    eprocess* p = e->processes[PT_SED][i];
    if (p->calc != NULL)
      p->calc(p, &ia);
  }

  c->ia = NULL;
  c->nsubstep++;
}

/** Flux calculation for an epibenthic cell.
 * @param t Time (input)
 * @param y Tracer concentrations (input)
 * @param y1 Fluxes (output)
 * @param p Pointer to cell
 */
static void epifluxes(double t, double* y, double* y1, void* pp)
{
  intargs ia;
  cell* c = (cell*) pp;
  ecology* e = c->col->e;
  int nvar = c->nvar;
  int i;
		
  ia.t = t;
  ia.y = y;
  ia.y1 = y1;
  ia.media = pp;

  for (i = 0; i < nvar; ++i)
    y1[i] = 0.0;

  if (c->wcchild != NULL) {
    cell* child = c->wcchild;

    ia.media = child;
    child->ia = &ia;
    for (i = 0; i < e->npr[PT_WC]; ++i) {
      eprocess* p = e->processes[PT_WC][i];
      if (p->calc != NULL)
	p->calc(p, &ia);
    }
  }

  c->ia = &ia;
  ia.media = pp;
  ia.y = y;
  ia.y1 = y1;
  for (i = 0; i < e->npr[PT_EPI]; ++i) {
    eprocess* p = e->processes[PT_EPI][i];

    if (p->calc != NULL) {
      p->calc(p, &ia);

      if (e->check_nans)
	check_nans(p, c);

    }
  }

  if (c->sedchild != NULL) {
    cell* child = c->sedchild;

    ia.media = child;
    child->ia = &ia;
    ia.y = &y[e->ntr];
    ia.y1 = &y1[e->ntr];
    for (i = 0; i < e->npr[PT_SED]; ++i) {
      eprocess* p = e->processes[PT_SED][i];

      if (p->calc != NULL) {
	p->calc(p, &ia);

	if (e->check_nans)
	  check_nans(p, child);
      }
    }
  }

  c->ia = NULL;
  if (c->wcchild != NULL)
    c->wcchild->ia = NULL;
  if (c->sedchild != NULL)
    c->sedchild->ia = NULL;

  c->nsubstep++;
}

static void addstats(int n, double x, double* y, double* y1, int naccpt, int nrejct, int nfcn, void* p)
{
  cell* c = (cell*) p;
  CELLTYPE ct = c->type;
  integration_stats* stats = &c->e->stepstats;

#if (NCPU > 1)
  pthstuff* pth = c->e->pth;
#endif

#if (NCPU > 1)
  pthread_mutex_lock(&pth->stats_lock);
#endif
  if (ct == CT_WC) {
    stats->nwc++;
    stats->nws += naccpt + nrejct;
  } else if (ct == CT_SED) {
    stats->nsc++;
    stats->nss += naccpt + nrejct;
  } else if (ct == CT_EPI) {
    stats->nec++;
    stats->nes += naccpt + nrejct;
  }
#if (NCPU > 1)
  pthread_mutex_unlock(&pth->stats_lock);
#endif
}

/** Integration loop for a cell.
 * @param c Pointer to cell
 * @return status flag, 0 for fail, 1 for success
 */
int cell_calc(cell* c)
{
  ecology* e = c->e;
  derivfn fluxes;

  if (c->type == CT_WC)
    fluxes = wcfluxes;
  else if (c->type == CT_SED)
    fluxes = sedfluxes;
  else   /* CT_EPI */
    fluxes = epifluxes;

  /*
   * integration 
   */
  if (!e->eint(fluxes, c->nvar, e->t, c->y, e->t + e->dt, e->precision, e->dt, c->h0, NULL, addstats, c)) {
    /* Integration failed */
    int ij[2];
    int b = c->col->b;
    void *hmodel = e->model;
    char str[MAXSTRLEN];
    
    ginterface_get_ij(e->model, b, ij);
    if (c->type == CT_WC) {
      sprintf(str, "ecology: error: integration failed in water cell, nstep = %d, nsubstep = %d, b = %d, k = %d (i,j) = (%u,%u)\n", e->nstep, c->nsubstep, b, c->k_wc, ij[0],ij[1]);
      i_set_error(e->model, b, LFATAL, str);
    } else if (c->type == CT_SED) {
      sprintf(str, "ecology: error: integration failed in sediment cell, nstep = %d, nsubstep = %d, b = %d, k = %d (i,j) = (%u,%u)\n", e->nstep, c->nsubstep, b, c->k_sed, ij[0],ij[1]);
      i_set_error(e->model, b, LFATAL, str);
    } else {
      sprintf(str, "ecology: error: integration failed in epibenthic cell, nstep = %d, nsubstep = %d, b = %d (i,j) = (%u,%u)\n", e->nstep, c->nsubstep, b, ij[0],ij[1]);
      i_set_error(e->model, b, LFATAL, str);
    }
    return(0);
  }
  
  /* Integration succeeeded */
  if (e->check_negs)
    check_negs(c,"calc");
  
  return(1);
}

/** Runs postcalc procedures for the cell (those to be performed after the
 * integration loop).
 * @param c Pointer to cell
 */
void cell_postcalc(cell* c)
{
  ecology* e = c->e;
  EPROCESSTYPE pt = celltype2processtype(c->type);
  int ndpe = 0;
  dpe dpes[NDPE_MAX];
  int i;
  char buf[MAXSTRLEN];
		
  if (c->wcchild != NULL) {
    cell* child = c->wcchild;

    for (i = 0; i < e->npr[PT_WC]; ++i) {
      eprocess* p = e->processes[PT_WC][i];
      sprintf(buf, "postcalc-%s-wcchild",p->name);
      if (p->postcalc != NULL) {
	if (p->delayed) {
	  dpes[ndpe].p = p;
	  dpes[ndpe].c = child;
	  ndpe++;
	} else {
	  p->postcalc(p, child);
	  if (e->check_nans)
	    check_nans(p, child);
	  if (e->check_negs)
	    check_negs(child,buf);
	}
      }
    }
  }

  for (i = 0; i < e->npr[pt]; ++i) {
    eprocess* p = e->processes[pt][i];
    sprintf(buf, "post-%s",p->name);
    if (p->postcalc != NULL) {
      if (p->delayed) {
	dpes[ndpe].p = p;
	dpes[ndpe].c = c;
	ndpe++;
      } else {
	p->postcalc(p, c);
	if (e->check_nans)
	  check_nans(p, c);
	if (e->check_negs)
	  check_negs(c,buf);
      }
    }
  }

  if (c->sedchild != NULL) {
    cell* child = c->sedchild;

    for (i = 0; i < e->npr[PT_SED]; ++i) {
      eprocess* p = e->processes[PT_SED][i];
      sprintf(buf, "post-%s-sedchild",p->name);
      if (p->postcalc != NULL) {
	if (p->delayed) {
	  dpes[ndpe].p = p;
	  dpes[ndpe].c = child;
	  ndpe++;
	} else {
	  p->postcalc(p, child);
	  if (e->check_nans)
	    check_nans(p, child);
	  if (e->check_negs)
	    check_negs(child,buf);
	}
      }
    }
  }

  for (i = 0; i < ndpe; ++i) {
    eprocess* p = dpes[i].p;
    sprintf(buf, "delayed-%s",p->name);
    p->postcalc(p, dpes[i].c);
    if (e->check_nans)
      check_nans(p, dpes[i].c);
    if (e->check_negs)
      check_negs(dpes[i].c,buf);
  }
}
