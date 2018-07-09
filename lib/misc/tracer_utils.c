/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/tracer_utils.c
 *
 *  \brief Interface routines to tracer_info data structure
 *
 *  This file combines some common utilities for tracers used by 
 *  the different hydrodynamic models.
 *  See description of `tracer_info' in "tracers.h".
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: tracer_utils.c 5864 2018-07-02 04:09:36Z her127 $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ems.h"
#include "errfn.h"


errfn keyprm_errfn;


/** Copy the contents of one tracer_info_t structure into another.
  * @param tr_cpy 'tracer_info_t' to copy into.
  * @param tr_in 'tracer_info_t' to copy from.
  */
void tracer_copy(tracer_info_t *tr_cpy, tracer_info_t *tr_in)
{
  int i;

  strcpy(tr_cpy->name, tr_in->name);
  strcpy(tr_cpy->long_name, tr_in->long_name);
  strcpy(tr_cpy->units, tr_in->units);
  strcpy(tr_cpy->std_name, tr_in->std_name);
  tr_cpy->valid_range_wc[0] = tr_in->valid_range_wc[0];
  tr_cpy->valid_range_wc[1] = tr_in->valid_range_wc[1];
  tr_cpy->valid_range_sed[0] = tr_in->valid_range_sed[0];
  tr_cpy->valid_range_sed[1] = tr_in->valid_range_sed[1];
  tr_cpy->fill_value_sed = tr_in->fill_value_sed;
  tr_cpy->type = tr_in->type;
  tr_cpy->diagn = tr_in->diagn;
  tr_cpy->fill_value_wc = tr_in->fill_value_wc;
  tr_cpy->n = tr_in->n;
  tr_cpy->m = tr_in->m;
  tr_cpy->inwc = tr_in->inwc;
  tr_cpy->insed = tr_in->insed;
  tr_cpy->dissol = tr_in->dissol;
  tr_cpy->advect = tr_in->advect;
  tr_cpy->diffuse = tr_in->diffuse;
  strcpy(tr_cpy->decay, tr_in->decay);
  strcpy(tr_cpy->i_rule, tr_in->i_rule);
  tr_cpy->partic = tr_in->partic;
  tr_cpy->increment = tr_in->increment;
  tr_cpy->flag = tr_in->flag;
  strcpy(tr_cpy->tag, tr_in->tag);
  strcpy(tr_cpy->tracerstat, tr_in->tracerstat);
  strcpy(tr_cpy->data, tr_in->data);
  strcpy(tr_cpy->relax_file, tr_in->relax_file);
  strcpy(tr_cpy->relax_dt, tr_in->relax_dt);
  strcpy(tr_cpy->r_rate, tr_in->r_rate);
  strcpy(tr_cpy->reset_file, tr_in->reset_file);
  strcpy(tr_cpy->reset_dt, tr_in->reset_dt);
  strcpy(tr_cpy->reset_interp, tr_in->reset_interp);
  strcpy(tr_cpy->vector_name, tr_in->vector_name);
  strcpy(tr_cpy->vector_components, tr_in->vector_components);
  tr_cpy->scale = tr_in->scale;
  tr_cpy->tctype = tr_in->tctype;
  tr_cpy->relax_dum = tr_in->relax_dum;
  for (i=0; i<TR_NPRV_DATA; i++) {
    // Delete dest
    if (tr_cpy->private_data[i] != NULL) {
      free(tr_cpy->private_data[i]);
      tr_cpy->private_data[i]      = NULL;
      tr_cpy->private_data_copy[i] = NULL;
    }
    // Copy over
    if (tr_in->private_data_copy[i]) {
      tr_cpy->private_data[i] = 
	              tr_in->private_data_copy[i](tr_in->private_data[i]);
      tr_cpy->private_data_copy[i] = tr_in->private_data_copy[i];
    }
  }
}


/** Finds index of a tracer in `tracer_info_t' array by name.
 * @param name tracer name
 * @param ntr number of tracers
 * @param tracers `tracer_info_t' array */
int tracer_find_index(const char *name, int ntr, tracer_info_t tracers[])
{
  int n;
  if (ntr <= 0 || tracers == NULL)
    return -1;
  for (n = 0; n < ntr; ++n)
    if (strcasecmp(tracers[n].name, name) == 0)
      return tracers[n].n;

  return -1;
}


/** Deallocates memory of a tracer_info array
 * 
 * @param ntr the int pointer to the parent array counter
 * @param tracers the content array
 */
void tracer_clear(int *ntr, tracer_info_t *tracers[])
{
  if ((*ntr == 0) || (*tracers == NULL))
    return;
/*
 * UR diabled, since tracer_info_t - name, long_name, unit was changed to 
 * fixed length arrays
  int n;

  for (n = 0; n < *ntr; ++n) {
    tracer_info_t *tr = &(*tracers)[n];
    free(tr->name);
    free(tr->long_name);
    free(tr->units);
  }
*/
  *ntr = 0;

  free(*tracers);
  *tracers = NULL;
}


/** Swaps two `tracer_info' entries in an array.
 * @param tinfo tracer_info array within to swap
 * @param i index to swap from
 * @param j index to swap to
 */
void tracer_swap(tracer_info_t tinfo[], int i, int j)
{
  tracer_info_t tmp = tinfo[i];
  tinfo[i] = tinfo[j];
  tinfo[j] = tmp;
  /* set new indices */
  tinfo[i].n = j;
  tinfo[j].n = i;
}

/* Reads a tracer attribute from a parameter file. Here, all
 * attributes are read and returned as a string. The calling
 * code then needs to parse the string as it sees fit.
 * @param fp parameter file
 * @param prefix input prefix like: "GRID0" or NULL
 * @param suffix input suffix
 * @param i tracer index
 * @param tag attribute name
 * @param fn error function
 * @param buffer output
 * @return non-zero if and only if successful
 */
int tracer_read_attribute(FILE * fp, 
			  const char *prefix,
			  const char *keyname, 
			  const char *trname, 
			  int i, 
			  const char *tag,
			  errfn fn, char *buffer)
{
  char keybuf[MAXLINELEN];
  errfn errfn_orig = keyprm_errfn;
  int ret = 0;

  prm_set_errfn(fn);
  /* By tracer number */
  ret =
    prm_read_char(fp,
		   prm_get_key(keybuf, prefix, "%s%1.1d.%s", keyname, i,
		 	       tag), buffer);

  /* By tracer name */
  if (trname != NULL) {
    ret |=
      prm_read_char(fp,
		    prm_get_key(keybuf, prefix, "%s.%s", trname, tag), buffer);
  }
  prm_set_errfn(errfn_orig);

  return ret;
}

