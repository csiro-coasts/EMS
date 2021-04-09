/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/tracers/tracer_info.c
 *  
 *  Description:
 *  Servicing of the tracer stuff.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: tracer_info.c 6654 2020-09-08 01:46:00Z her127 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "netcdf.h"
#include "hd.h"

errfn keyprm_errfn;
int get_type(char *buf, int type);
int is_type(int type, char *stype);

/** Creates a `tracer_info' array from a PRM file.
 * @param fp parameter file
 * @param prefix input prefix like: "GRID0" or NULL
 * @param type of tracer to extract from 'type' attribute
 * @param quitfn exit function
 * @param warnfn warn function
 * @param emptyfn empty function
 * @param ntr output -- number of tracers
 * @param tracers output -- array of tracer_desc
 */
void tracer_read(FILE * fpp, char *prefix, int type, errfn quitfn,
                 errfn warnfn, errfn emptyfn, int *ntr,
                 tracer_info_t *tracers[])
{
  int i, n, m, intr, bt, ns = *ntr;
  char buf[MAXLINELEN], key[MAXLINELEN], keyname[MAXLINELEN];
  char sed_def[MAXSTRLEN], eco_def[MAXSTRLEN];
  FILE* fpt = fpp;
  int trtype = 0;
	
  *tracers = NULL;

  /* Sanity check */
  if (fpt == NULL)
    quitfn("tracer_read: NULL parameter file pointer\n");

  /* Read the sediment default function */
  sprintf(sed_def, "%c", '\0');
  prm_read_char(fpt, "SED_VARS_ATTS", sed_def);
  sprintf(eco_def, "%c", '\0');
  prm_read_char(fpt, "ECO_VARS_ATTS", eco_def);

  /* Read the number of tracers */
  /* prm_set_errfn(quitfn); */
  strcpy(keyname, "NTRACERS");

  intr = 0;
  if(!prm_read_char(fpt, "NTRACERS", buf)) {
    if (prm_read_char(fpt, "TRACERFILE", buf)) {
      fpt= fopen(buf,"r");
      if(fpt == NULL)
	hd_quit("Cannot open tracer file '%s'",buf);
      prm_read_int(fpt, prm_get_key(buf, prefix, keyname), &intr);
    }		
  } else
    prm_read_int(fpt, prm_get_key(buf, prefix, keyname), &intr);	

  if (intr < 0)
    quitfn("tracer_read: %s < 0\n", keyname);
  else if (intr == 0)
    intr = *ntr;
  if(intr == 0)
    return;

  /* Count the number of tracer of type */
  strcpy(keyname, "TRACER");
  for (m = 0; m < intr; ++m) {
    if(tracer_read_attribute(fpt, prefix, keyname, NULL, m, "type", NULL,
			     buf)) {
      if(type == WATER && ((contains_token(buf, "WATER") != NULL) ||
			   (contains_token(buf, "WC3D") != NULL) ||
			   (contains_token(buf, "WC") != NULL)))*ntr += 1;

      if(type == INTER && ((contains_token(buf, "BENTHIC") != NULL) ||
			   (contains_token(buf, "INTER") != NULL) ||
			   (contains_token(buf, "WC2D") != NULL)))*ntr += 1;

      if(type == SEDIM && ((contains_token(buf, "SEDIMENT") != NULL) ||
			   (contains_token(buf, "SED") != NULL)))*ntr += 1;
    }
    else {
#if defined(HAVE_ECOLOGY_MODULE)
      /*
       * This is an attempt to front load ecology autotracer type,
       * if it is explicitly specified in the parameter file
       */
      tracer_read_attribute(fpt, prefix, keyname, NULL, m, "name", quitfn,
			    buf);
      if (is_eco_var(buf) && strlen(eco_def)) {
	trtype = get_eco_var_type(buf, eco_def);
	
	if (type == WATER && (trtype & WATER))
	  *ntr += 1;
	if (type == SEDIM && (trtype & SEDIM))
	  *ntr += 1;
	if (type == INTER && (trtype & INTER))
	  *ntr += 1;
      } else {
	if (type == WATER)
	  *ntr += 1;
      }
#else
      if (type == WATER)
	*ntr += 1;
#endif
    }
  }

  /* Allocate memory for array of tracer_info and fill with zeros */
  {
    size_t sz;
    sz = sizeof(tracer_info_t) * *ntr;
    *tracers = (tracer_info_t *)malloc(sz);
    memset(*tracers, 0, sz);
  }

  /* Read tracer attributes for all tracers */
  n = ns - 1;
  for (m = 0; m < intr; ++m) {
    tracer_info_t *tr;

    bt = 0;
    if(tracer_read_attribute(fpt, prefix, keyname, NULL, m, "type", NULL,
			     buf)) {
      if(type == WATER && ((contains_token(buf, "WATER") != NULL) ||
			   (contains_token(buf, "WC3D") != NULL) ||
			   (contains_token(buf, "WC") != NULL)))bt = WATER;
      if(type == INTER && ((contains_token(buf, "BENTHIC") != NULL) ||
			   (contains_token(buf, "INTER") != NULL) ||
			   (contains_token(buf, "WC2D") != NULL)))bt = INTER;
      if(type == SEDIM && ((contains_token(buf, "SEDIMENT") != NULL) ||
			   (contains_token(buf, "SED") != NULL)))bt = SEDIM;
    }
    else {
#if defined(HAVE_ECOLOGY_MODULE)
      tracer_read_attribute(fpt, prefix, keyname, NULL, m, "name", quitfn,
			    buf);
      if (is_eco_var(buf) && strlen(eco_def)) {
	trtype = get_eco_var_type(buf, eco_def);
	
	if (type == WATER && (trtype & WATER))
	  bt = WATER;
	if (type == SEDIM && (trtype & SEDIM))
	  bt = SEDIM;
	if (type == INTER && (trtype & INTER))
	  bt = INTER;
      } else {
	if(type == WATER)
	  bt = WATER;
      }
#else
      if(type == WATER)
	bt = WATER;
#endif
    }
    if(bt == 0)continue;

    n++;
    tr = &(*tracers)[n];
    tr->n = n;
    tr->m = m;
    tr->type = bt;
    tr->type = get_type(buf, tr->type);
    tracer_read_attribute(fpt, prefix, keyname, NULL, m, "name", quitfn,
                          buf);
    strcpy(tr->name, buf);

    if (tracer_read_attribute (fpt, prefix, keyname, tr->name, m, "long_name", warnfn, buf))
      strcpy(tr->long_name, buf);
    else
      strcpy(tr->long_name, tr->name);

    if (tracer_read_attribute(fpt, prefix, keyname, tr->name, m, "diagn", emptyfn, buf)) {
      if(sscanf(buf, "%d", &i))
        tr->diagn = i;
      else
	      tr->diagn = is_true(buf);
    }
    else
      tr->diagn = 0;

    if (!tr->diagn) {
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "units", warnfn, buf))
        strcpy(tr->units, buf);
      else
        strcpy(tr->units, "1");
    } else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "units", emptyfn, buf))
      strcpy(tr->units, buf);
    else
      strcpy(tr->units, "1");

    /*UR 6/6/2006 changed default limits to allow negative values */
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "valid_range_wc", emptyfn, buf))
      sscanf(buf, "%lf %lf", &tr->valid_range_wc[0],
             &tr->valid_range_wc[1]);
    else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "valid_range", emptyfn, buf))
      sscanf(buf, "%lf %lf", &tr->valid_range_wc[0],
             &tr->valid_range_wc[1]);
    else {
      tr->valid_range_wc[0] = -1.0e35;
      tr->valid_range_wc[1] = 1.0e+35;
    }

    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "valid_range_sed", emptyfn, buf))
      sscanf(buf, "%lf %lf", &tr->valid_range_sed[0],
             &tr->valid_range_sed[1]);
    else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "valid_range", emptyfn, buf))
      sscanf(buf, "%lf %lf", &tr->valid_range_sed[0],
             &tr->valid_range_sed[1]);
    else {
      tr->valid_range_wc[0] = -1.0e35;
      tr->valid_range_wc[1] = 1.0e+35;
    }

    if (tr->diagn) {
      tr->fill_value_wc = 0.0;
      tr->fill_value_sed = 0.0;
      tr->inwc = 0;
      tr->insed = 0;
      tr->dissol = 0;
      tr->partic = 0;
      tr->advect = 0;
      tr->diffuse = 0;
      sprintf(tr->decay, "%c", '\0');
      /*continue;*/
    }
    tr->flag = 0;
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name,m, "tag", emptyfn, buf)) {
      int id;
      double val;
      char buf1[MAXSTRLEN];
      strcpy(tr->tag, buf);
      /* Decode the tag and set flags */
      if (decode_tag(buf, "scale_s", key)) {
	if (sscanf(key, "%lf", &val) == 1) {
	  tr->scale = val;
	  tr->flag |= SC_SC;
	} else {
	  for (id = 0; id < intr; ++id) {
	    tracer_read_attribute(fpt, prefix, keyname, NULL, id, "name", emptyfn, buf1);
	    if (strcmp(buf1, key) == 0) {
	      tr->scale = (double)(id+ns);
	      tr->flag |= SC_ST;
	      break;
	    }
	  }
	}
      }
      if (decode_tag(buf, "scale_p", key)) {
	if (sscanf(key, "%lf", &val) == 1) {
	  tr->scale = val;
	  tr->flag |= SC_PC;
	} else {
	  for (id = 0; id < intr; ++id) {
	    tracer_read_attribute(fpt, prefix, keyname, NULL, id, "name", emptyfn, buf1);
	    if (strcmp(buf1, key) == 0) {
	      tr->scale = (double)(id+ns);
	      tr->flag |= SC_PT;
	      break;
	    }
	  }
	} 
      }
      if (decode_tag(buf, "relax_nograd", key)) {
	tr->flag |= RLX_GRD;
	tr->relax_dum = atof(key);
      }
      if (decode_tag(buf, "hipass_vert", key)) {
	if (sscanf(key, "%lf", &val) == 1) {
	  tr->scale = val;
	  tr->flag |= V_HP;
	}
      }
    }
    /* Scaling */
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "scale_s", emptyfn, buf)) {
      int id, idt = 0;
      char buf1[MAXSTRLEN];
      for (id = 0; id < intr; ++id) {
	tracer_read_attribute(fpt, prefix, keyname, NULL, id, "name", emptyfn, key);
	tracer_read_attribute(fpt, prefix, keyname, NULL, id, "type", emptyfn, buf1);
	if (strcmp(buf, key) == 0) {
	  tracer_info_t *trinfo = &(*tracers)[idt+ns];
	  trinfo->scale = (double)n;
	  trinfo->flag |= SC_ST;
	  break;
	}
	if (is_type(type, buf1)) idt++;
      }
    }
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "scale_p", emptyfn, buf)) {
      int id;
      for (id = 0; id < intr; ++id) {
	tracer_read_attribute(fpt, prefix, keyname, NULL, id, "name", emptyfn, key);
	if (strcmp(buf, key) == 0) {
	  tracer_info_t *trinfo = &(*tracers)[id+ns];
	  trinfo->scale = (double)n;
	  trinfo->flag |= SC_PT;
	  break;
	}
      }
    }

    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "fill_value_wc", emptyfn, buf))
      tr->fill_value_wc = atof(buf);
    else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "fill_value", warnfn, buf))
      tr->fill_value_wc = atof(buf);
    else
      tr->fill_value_wc = 0.0;

    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "fill_value_sed", emptyfn, buf))
      tr->fill_value_sed = atof(buf);
    else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "fill_value", warnfn, buf))
      tr->fill_value_sed = atof(buf);
    else
      tr->fill_value_sed = 0.0;

    /* Search for a data file for this variable in the prm file     */
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "interp_type", warnfn, buf))
      strcpy(tr->i_rule, buf);
    else
      tr->i_rule[0] = '\0';

/*UR added 5/2008 */
    if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "missing_value", warnfn, buf))
      tr->missing_value = atof(buf);
    else
      tr->missing_value = NaN;


    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "missing_value_sed", emptyfn, buf))
      tr->missing_value_sed = atof(buf);
    else
      tr->missing_value_sed = tr->missing_value;

/*UR end added */


    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "increment", emptyfn, buf))
      tr->increment = get_inc(buf);
    else
      tr->increment = 0;

    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "inwc", emptyfn, buf))
      tr->inwc = is_true(buf);
    else
      tr->inwc = 1;
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "insed", emptyfn, buf))
      tr->insed = is_true(buf);
    else
      tr->insed = 1;
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "dissol", emptyfn, buf)) {
      tr->dissol = is_true(buf);
      tr->partic = !(is_true(buf));
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "partic", emptyfn, buf)) {
        tr->partic = is_true(buf);
        if (tr->partic && tr->dissol)
          quitfn("tracer_read: %s set both dissolved and particulate\n",
                 tr->name);
      }
    } else
      if (tracer_read_attribute
          (fpt, prefix, keyname, tr->name, m, "partic", emptyfn, buf)) {
      tr->partic = is_true(buf);
      tr->dissol = !(is_true(buf));
    } else {
      tr->dissol = 1;
      tr->partic = 0;
    }
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "advect", emptyfn, buf)) {
      char *check;
      tr->advect = strtol(buf, &check, 10);
/*UR this statement is pointless ?? */
      if (*check != 0 || tr->advect == LONG_MIN || tr->advect == LONG_MAX)
        tr->advect = is_true(buf);
    } else
      tr->advect = (tr->diagn > 0) ? 0 : 1;
    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "diffuse", emptyfn, buf))
      tr->diffuse = is_true(buf);
    else
      tr->diffuse = (tr->diagn > 0) ? 0 : 1;

    if (tracer_read_attribute
        (fpt, prefix, keyname, tr->name, m, "decay", emptyfn, buf)) {
      strcpy(tr->decay, buf);
    } else
      strcpy(tr->decay, "0.0");

    /* Read the tracer statistic. */
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "tracerstat", emptyfn, buf))
      strcpy(tr->tracerstat, buf);
    else
      sprintf(tr->tracerstat, "%c", '\0');
    /*memcpy(tr->tracerstat, 0, sizeof(tr->tracerstat));*/

    /* Read the relaxation data file. */
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "relaxation_file", emptyfn, buf))
      strcpy(tr->relax_file, buf);
    else
      sprintf(tr->relax_file, "%c", '\0');
    /*memcpy(tr->relax_file, 0, sizeof(tr->relax_file));*/
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "relaxation_input_dt", emptyfn, buf))
      strcpy(tr->relax_dt, buf);
    else
      sprintf(tr->relax_dt, "%c", '\0');
    /*memcpy(tr->relax_dt, 0, sizeof(tr->relax_dt));*/
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "relaxation_time_constant", emptyfn, buf))
      strcpy(tr->r_rate, buf);
    else
      sprintf(tr->r_rate, "%c", '\0');

    sprintf(key, "relax_%s", tr->name);
    if (prm_read_char(fpt, key, buf)) {
      char fname[MAXSTRLEN], dt[MAXSTRLEN], unit[MAXSTRLEN];
      double in_dt;
      if (sscanf(buf,"%s %s %lf %s", fname, dt, &in_dt, unit) == 4) {
	strcpy(tr->relax_file, fname);
	sprintf(buf,"%f %s\n",in_dt, unit);
	strcpy(tr->relax_dt, buf);
	strcpy(tr->r_rate, dt);
      } else
	hd_quit("relax_%s format = filename.ts dt.ts n days\n");      
    }

    /* Read sedimet defaults */
#if defined(HAVE_SEDIMENT_MODULE)
    if (is_sed_var(tr->name) 
#if defined(HAVE_ECOLOGY_MODULE)
	|| is_eco_var(tr->name)
#endif
	) {
      if (tracer_read_attribute
	  (fpt, prefix, keyname, tr->name, m, "sed_default", emptyfn, buf))
	sed_set_tracer_defaults(tr, tr->name, buf);
      else if (strlen(sed_def))
	sed_set_tracer_defaults(tr, tr->name, sed_def);
      sprintf(buf, "%s%1.1d", keyname, m);
      /* Define attributes by tracer number */
      sed_read_tr_atts(tr, fpt, buf);
      /* Define attributes by tracer name */
      sed_read_tr_atts(tr, fpt, tr->name);
    }
#endif

#if defined(HAVE_ECOLOGY_MODULE)
    if (is_eco_var(tr->name)) {
      if (tracer_read_attribute
	  (fpt, prefix, keyname, tr->name, m, "eco_default", emptyfn, buf))
	eco_set_tracer_defaults(tr, tr->name, buf, NULL);
      else if (strlen(eco_def))
	eco_set_tracer_defaults(tr, tr->name, eco_def, NULL);
      sprintf(buf, "%s%1.1d", keyname, m);
      /* Define attributes by tracer number */
      eco_read_tr_atts(tr, fpt, buf);
      /* Define attributes by tracer name */
      eco_read_tr_atts(tr, fpt, tr->name);
    }
#endif

    /* Read the initialisation data. */
    if (type == SEDIM)
      strcpy(key, "data_sed");
    else
      strcpy(key, "data");
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, key, emptyfn, buf))
      strcpy(tr->data, buf);
    else
      sprintf(tr->data, "%c", '\0');

    /* Read the reset data file. */
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "reset_file", emptyfn, buf))
      strcpy(tr->reset_file, buf);
    else
      sprintf(tr->reset_file, "%c", '\0');
    /*memcpy(tr->reset_file, 0, sizeof(tr->reset_file));*/
    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "reset_dt", emptyfn, buf))
      strcpy(tr->reset_dt, buf);
    else
      sprintf(tr->reset_dt, "%c", '\0');

    if (tracer_read_attribute
	(fpt, prefix, keyname, tr->name, m, "reset_interp", emptyfn, buf))
      strcpy(tr->reset_interp, buf);
    else
      sprintf(tr->reset_interp, "%c", '\0');
    /*memcpy(tr->reset_dt, 0, sizeof(tr->reset_dt));*/
  }
  
  if(fpt != fpp)
		fclose(fpt);
}

int get_inc(char *buf)
{
  if (strcmp(buf, "TEMP") == 0 || strcmp(buf, "temp") == 0)
    return(TEMP);
  if (strcmp(buf, "SALT") == 0 || strcmp(buf, "salt") == 0)
    return(SALT);
  if (strcmp(buf, "ETA") == 0 || strcmp(buf, "eta") == 0)
    return(ETA_M);
  if (strcmp(buf, "U1") == 0 || strcmp(buf, "u1") == 0)
    return(U1VEL);
  if (strcmp(buf, "U2") == 0 || strcmp(buf, "u2") == 0)
    return(U2VEL);
  return(0);
}

void tracer_re_read(tracer_info_t *tr, FILE *fpt, int trinfo_type)
{
  char buf[MAXSTRLEN], key[MAXSTRLEN], keyname[MAXSTRLEN];
  int i, m = tr->m;
  errfn warn = hd_silent_warn;

  strcpy(key, "TRACER");

  if (tracer_read_attribute (fpt, NULL, key, tr->name, m, "long_name", warn, buf))
    strcpy(tr->long_name, buf);

  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "diagn", warn, buf)) {
    if(sscanf(buf, "%d", &i))
      tr->diagn = i;
    else
      tr->diagn = is_true(buf);
  }

  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "units", warn, buf))
    strcpy(tr->units, buf);

  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "valid_range_wc", warn, buf))
    sscanf(buf, "%lf %lf", &tr->valid_range_wc[0], &tr->valid_range_wc[1]);             
  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "valid_range", warn, buf))
    sscanf(buf, "%lf %lf", &tr->valid_range_wc[0], &tr->valid_range_wc[1]);
  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "valid_range_sed", warn, buf))
    sscanf(buf, "%lf %lf", &tr->valid_range_sed[0], &tr->valid_range_sed[1]);

  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "fill_value_wc", warn, buf))
    tr->fill_value_wc = atof(buf);
  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "fill_value", warn, buf))          
    tr->fill_value_wc = atof(buf);
  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "fill_value_sed", warn, buf))
    tr->fill_value_sed = atof(buf);

  if (tr->type & (WATER|INTER) && trinfo_type & (WATER|INTER)) {
    if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "data", warn, buf))
      strcpy(tr->data, buf);
  }
  if (tr->type & SEDIM && trinfo_type & SEDIM) {
    if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "data_sed", warn, buf))
      strcpy(tr->data, buf);
  }

  if (tracer_read_attribute(fpt, NULL, key, tr->name, m, "tracerstat", warn, buf))
    strcpy(tr->tracerstat, buf);

  if (tracer_read_attribute(fpt, NULL, keyname, tr->name, m, "relaxation_file", warn, buf))
    strcpy(tr->relax_file, buf);
  if (tracer_read_attribute(fpt, NULL, keyname, tr->name, m, "relaxation_input_dt", warn, buf))
    strcpy(tr->relax_dt, buf);
  if (tracer_read_attribute(fpt, NULL, keyname, tr->name, m, "relaxation_time_constant", warn, buf))
    strcpy(tr->r_rate, buf);
  sprintf(keyname, "relax_%s", tr->name);
  if (prm_read_char(fpt, keyname, buf)) {
    char fname[MAXSTRLEN], dt[MAXSTRLEN], unit[MAXSTRLEN];
    double in_dt;
    if (sscanf(buf,"%s %s %lf %s", fname, dt, &in_dt, unit) == 4) {
      strcpy(tr->relax_file, fname);
      sprintf(buf,"%f %s\n",in_dt, unit);
      strcpy(tr->relax_dt, buf);
      strcpy(tr->r_rate, dt);
    }
  }

  if (tracer_read_attribute(fpt, NULL, keyname, tr->name, m, "reset_file", warn, buf))
    strcpy(tr->reset_file, buf);
  if (tracer_read_attribute(fpt, NULL, keyname, tr->name, m, "reset_dt", warn, buf))
    strcpy(tr->reset_dt, buf);

}

/** Reads the `tracer_info_t' array from a netCDF file.
 * Distinguishes tracers from other variables by:
 * -- number of dimensions being equal to `ndims' (if not 0)
 * -- dimension ids being equal `dimids[]' (if not NULL)
 * -- presence of attribute `attr->attr' (if `attr' is not NULL)
 * -- value of attribute `attr->attr' being equal to `attr->attrval'
 *    (if `attr' is not NULL and `attr->attrval' not "").
 *
 * All the variables are assumed to have attributes coinciding with
 * fields in `tracer_info_t' structure with two exceptions:
 *     _FillValueWC instead of fill_value_wc
 *     _FillValueSED instead of fill_value_sed
 * @param fid NetCDF file
 * @param ndims number of dimensions for tracers (e.g. 4 for x,y,z,t-geometry)
 * @param dimids NetCDF ids for dimension variables
 * @param attr attribute for tracers (e.g. {"tracer", "true"})
 * @param tracers output (array of `tracer_desc')
 * @param ntr output (number of tracers)
 */
void tracer_read_nc(int fid, int ndims, int dimids[], tracer_att_t *attr,
                    int *ntr, tracer_info_t *tracers[])
{
  int varids[MAXNUMTRACERS];
  int n;
  size_t sz;

  *tracers = NULL;

  /* Look for all variables with a tracer attribute */
  if (attr == NULL)
    *ntr = ncw_var_find(fid, ndims, dimids, NULL, NULL, varids);
  else if (attr->attrval[0] == 0)
    *ntr = ncw_var_find(fid, ndims, dimids, attr->attr, NULL, varids);
  else
    *ntr =
      ncw_var_find(fid, ndims, dimids, attr->attr, attr->attrval, varids);

  if (*ntr <= 0)
    /* Nothing to do */
    return;

  /* Tracers were found, so allocate memory for tracer list, and set it
     all to zero. */
  sz = sizeof(tracer_info_t) * *ntr;
  *tracers = (tracer_info_t *)malloc(sz);
  memset(*tracers, 0, sz);

  /* Read tracer attributes for all tracers */
  for (n = 0; n < *ntr; ++n) {
    tracer_info_t *tr = &(*tracers)[n];
    int varid = varids[n];
    char tmpstr[MAXLINELEN];

    tr->n = n;

    memset(tmpstr, 0, MAXLINELEN);
    nc_inq_varname(fid, varid, tmpstr);
    strcpy(tr->name, tmpstr);

    memset(tmpstr, 0, MAXLINELEN);
    nc_get_att_text(fid, varid, "long_name", tmpstr);
    strcpy(tr->long_name, tmpstr);

    memset(tmpstr, 0, MAXLINELEN);
    nc_get_att_text(fid, varid, "units", tmpstr);
    strcpy(tr->units, tmpstr);

    if (nc_get_att_int(fid, varid, "type", &tr->type) != NC_NOERR)
      tr->type = WATER;

    if (nc_get_att_double(fid, varid, "valid_range_wc", tr->valid_range_wc)
        != NC_NOERR)
      nc_get_att_double(fid, varid, "valid_range", tr->valid_range_wc);

    if (nc_get_att_double
        (fid, varid, "valid_range_sed", tr->valid_range_sed) != NC_NOERR)
      nc_get_att_double(fid, varid, "valid_range", tr->valid_range_sed);

    if (nc_get_att_int(fid, varid, "diagn", &tr->diagn) != NC_NOERR)
      tr->diagn = 0;
    if (nc_get_att_double(fid, varid, "_FillValueWC", &tr->fill_value_wc)
        != NC_NOERR)
      nc_get_att_double(fid, varid, "_FillValue", &tr->fill_value_wc);

    if (nc_get_att_double(fid, varid, "_FillValueSED", &tr->fill_value_sed)
        != NC_NOERR)
      nc_get_att_double(fid, varid, "_FillValue", &tr->fill_value_sed);

    nc_get_att_int(fid, varid, "inwc", &tr->inwc);
    nc_get_att_int(fid, varid, "insed", &tr->insed);
    nc_get_att_int(fid, varid, "dissol", &tr->dissol);
    nc_get_att_int(fid, varid, "partic", &tr->partic);
    nc_get_att_int(fid, varid, "advect", &tr->advect);
    nc_get_att_int(fid, varid, "diffuse", &tr->diffuse);
    memset(tmpstr, 0, MAXLINELEN);
    nc_get_att_text(fid, varid, "decay", tmpstr);
    if (strlen(tmpstr))
      strcpy(tr->decay, tmpstr);
    /* Get tag, if it exists */
    memset(tmpstr, 0, MAXLINELEN);
    nc_get_att_text(fid, varid, "tag", tmpstr);
    if (strlen(tmpstr))
      strcpy(tr->tag, tmpstr);
  }
}

/** Writes tracer attributes to a NetCDF file.
 * `attr' and `attrval' are common attributes added to all tracers.
 * @param fid NetCDF file
 * @param ntr number of tracers
 * @param tracers array of `tracer_desc'
 * @param attr array of additional attributes to be added to all tracers
 */
void tracer_write_nc(int fid, int ntr, tracer_info_t tracers[], int nattr,
                     tracer_att_t attr[])
{
  int n;
  char name[MAXSTRLEN];

  for (n = 0; n < ntr; n++) {
    int vid;
    tracer_info_t *tr = &tracers[n];
    strcpy(name,tr->name);
    if (!dumpdata->large && (strcmp(name, "salt") != 0) &&
	(strcmp(name, "temp") != 0)) continue;
    if (strcmp(attr[0].attr,"tracer_sed") == 0)
      strcat(name, "_sed");
    vid = ncw_var_id(fid, name);
    if (vid >= 0) {
      int i;
      if (nattr == 0) {
        nc_put_att_text(fid, vid, "tracer", strlen("true"), "true");
      }
      else
        for (i = 0; i < nattr; ++i)
          nc_put_att_text(fid, vid, attr[i].attr, strlen(attr[i].attrval),
                          attr[i].attrval);

      nc_put_att_text(fid, vid, "long_name", strlen(tr->long_name),
                      tr->long_name);
      nc_put_att_text(fid, vid, "units", strlen(tr->units), tr->units);
      nc_put_att_int(fid, vid, "type", NC_INT, 1, &tr->type);
      nc_put_att_int(fid, vid, "diagn", NC_INT, 1, &tr->diagn);
      nc_put_att_double(fid, vid, "_FillValueWC", NC_DOUBLE, 1,
                        &tr->fill_value_wc);
      nc_put_att_double(fid, vid, "valid_range_wc", NC_DOUBLE, 2,
                        tr->valid_range_wc);
      if(tr->type & INTER)continue;
      nc_put_att_double(fid, vid, "_FillValueSED", NC_DOUBLE, 1,
                        &tr->fill_value_sed);
      nc_put_att_double(fid, vid, "valid_range_sed", NC_DOUBLE, 2,
                        tr->valid_range_sed);
      nc_put_att_int(fid, vid, "inwc", NC_INT, 1, &tr->inwc);
      nc_put_att_int(fid, vid, "insed", NC_INT, 1, &tr->insed);
      nc_put_att_int(fid, vid, "dissol", NC_INT, 1, &tr->dissol);
      nc_put_att_int(fid, vid, "partic", NC_INT, 1, &tr->partic);
      nc_put_att_int(fid, vid, "advect", NC_INT, 1, &tr->advect);
      nc_put_att_int(fid, vid, "diffuse", NC_INT, 1, &tr->diffuse);
      if (strlen(tr->decay))
	nc_put_att_text(fid, vid, "decay", strlen(tr->decay), tr->decay);
      /* Add tag, if it exists */
      if (strlen(tr->tag))
	nc_put_att_text(fid, vid, "tag", strlen(tr->tag), tr->tag);
    }
  }
}

/*
 * Helper function to write out 2D tracer attributes
 * Called from below as well as sediment and ecology
 */
void tracer_write_2d(master_t *master, FILE *op, tracer_info_t *tracer, int n)
{
  parameters_t *params = master->params;
  char key[MAXSTRLEN];

  fprintf(op, "TRACER%1.1d.name            %s\n", n, tracer->name);
  fprintf(op, "TRACER%1.1d.long_name       %s\n", n, tracer->long_name);
  fprintf(op, "TRACER%1.1d.units           %s\n", n, tracer->units);
  if (tracer->fill_value_wc < 1)
    fprintf(op, "TRACER%1.1d.fill_value      %-4.2e\n", n,
	    tracer->fill_value_wc);
  else
    fprintf(op, "TRACER%1.1d.fill_value      %-6.1f\n", n,
	    tracer->fill_value_wc);
  if (n > 1)
    fprintf(op, "TRACER%1.1d.valid_range     %-4.2e %-4.2e\n", n,
	    tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
  else
    fprintf(op, "TRACER%1.1d.valid_range     %-6.1f %-6.1f\n", n,
	    tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
  fprintf(op, "TRACER%1.1d.diagn           %d\n", n, tracer->diagn);
  fprintf(op, "TRACER%1.1d.type           %s\n", n, trtypename(tracer->type, key));
  fprintf(op, "\n");
  
}

int is_type(int type, char *stype)
{
  int ret = WATER;

  if((contains_token(stype, "WATER") != NULL) ||
     (contains_token(stype, "WC3D") != NULL) ||
     (contains_token(stype, "WC") != NULL)) ret = WATER;
  if((contains_token(stype, "BENTHIC") != NULL) ||
     (contains_token(stype, "INTER") != NULL) ||
     (contains_token(stype, "WC2D") != NULL)) ret = INTER;
  if((contains_token(stype, "SEDIMENT") != NULL) ||
     (contains_token(stype, "SED") != NULL)) ret = SEDIM;
  if (ret == type) 
    return(1);
  else
    return(0);
}

int get_type(char *buf, int type)
{
  if (contains_token(buf, "HYDRO") != NULL)
    type |= HYDRO;
  if (contains_token(buf, "SEDIMENT") != NULL)
    type |= SEDIMENT;
  if (contains_token(buf, "ECOLOGY") != NULL)
    type |= ECOLOGY;
  if (contains_token(buf, "WAVE") != NULL)
    type |= WAVE;
  if (contains_token(buf, "TRACERSTAT") != NULL)
    type |= TRACERSTAT;
  if (contains_token(buf, "PROGNOSTIC") != NULL)
    type |= PROGNOSTIC;
  if (contains_token(buf, "DIAGNOSTIC") != NULL)
    type |= DIAGNOSTIC;
  if (contains_token(buf, "PARAMETER") != NULL)
    type |= PARAMETER;
  if (contains_token(buf, "FORCING") != NULL)
    type |= FORCING;
  return (type);
}

/*
 * Helper function to write out 3D tracer attributes
 * Called from below as well as sediment and ecology
 */
void tracer_write_3d(master_t *master, FILE *op, tracer_info_t *tracer, int n)
{
  parameters_t *params = master->params;
  int sn;
  char tag[MAXSTRLEN];
  char key[MAXSTRLEN];

  fprintf(op, "TRACER%1.1d.name            %s\n", n, tracer->name);
  fprintf(op, "TRACER%1.1d.long_name       %s\n", n, tracer->long_name);
  fprintf(op, "TRACER%1.1d.units           %s\n", n, tracer->units);
  if (tracer->fill_value_wc < 1)
    fprintf(op, "TRACER%1.1d.fill_value      %-4.2e\n", n,
	    tracer->fill_value_wc);
  else
    fprintf(op, "TRACER%1.1d.fill_value      %-6.1f\n", n,
	    tracer->fill_value_wc);
  if (n > 1)
    fprintf(op, "TRACER%1.1d.valid_range     %-4.2e %-4.2e\n", n,
	    tracer->valid_range_wc[0], tracer->valid_range_wc[1]);
  else
    fprintf(op, "TRACER%1.1d.valid_range     %-6.1f %-6.1f\n", n,
	    tracer->valid_range_wc[0], tracer->valid_range_wc[1]);

  if ((sn = tracer_find_index(tracer->name, master->nsed, master->trinfo_sed)) >= 0) {
    tracer_info_t *tr = &master->trinfo_sed[sn];
    if (tr->fill_value_wc < 1)
      fprintf(op, "TRACER%1.1d.fill_value_sed  %-4.2e\n", n, tr->fill_value_sed);
    else
      fprintf(op, "TRACER%1.1d.fill_value_sed  %-6.1f\n", n, tr->fill_value_sed);
    if (n > 1)
      fprintf(op, "TRACER%1.1d.valid_range_sed %-4.2e %-4.2e\n", n,
	      tr->valid_range_sed[0], tr->valid_range_sed[1]);
    else
      fprintf(op, "TRACER%1.1d.valid_range_sed %-6.1f %-6.1f\n", n,
	      tr->valid_range_sed[0], tr->valid_range_sed[1]);
  }
  fprintf(op, "TRACER%1.1d.advect          %d\n", n, tracer->advect);
  fprintf(op, "TRACER%1.1d.diffuse         %d\n", n, tracer->diffuse);
  fprintf(op, "TRACER%1.1d.diagn           %d\n", n, tracer->diagn);
  fprintf(op, "TRACER%1.1d.partic          %d\n", n, tracer->partic);
  fprintf(op, "TRACER%1.1d.dissol          %d\n", n, tracer->dissol);
  
#if defined(HAVE_SEDIMENT_MODULE)
  if (strcmp(tracer->name, "temp") == 0 && params->do_sed)
    fprintf(op, "TRACER%1.1d.type            WATER SEDIMENT\n", n);
  else if (strcmp(tracer->name, "salt") == 0 && params->do_sed)
    fprintf(op, "TRACER%1.1d.type            WATER SEDIMENT\n", n);
  else if (params->do_sed && tracer->type == (WATER|SEDIM))
    fprintf(op, "TRACER%1.1d.type            WATER SEDIMENT\n", n);
  else 
#endif
    if(tracer->type && (tracer->type & (WATER|SEDIM) || tracer->type & INTER))
      fprintf(op, "TRACER%1.1d.type           %s\n", n, trtypename(tracer->type, key));
  
  if (strlen(tracer->tracerstat))
    fprintf(op, "TRACER%1.1d.tracerstat      %s\n", n, tracer->tracerstat);
#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed)
    sed_write_tr_atts(tracer, op, n);
#endif
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco)
    eco_write_tr_atts(tracer, op, n);
#endif
  if (strlen(tracer->data))
    fprintf(op, "TRACER%1.1d.data            %s\n", n, tracer->data);
  if (params->trrlxn) {
    if (strlen(params->trrlxn[n])) {
      fprintf(op, "TRACER%1.1d.relaxation_file        %s\n",
	      n, params->trrlxn[n]);
      fprintf(op, "TRACER%1.1d.relaxation_input_dt    %s\n",
	      n, otime(params->trrlxdt[n], tag));
      fprintf(op, "TRACER%1.1d.relaxation_time_constant     %s\n",
	      n, params->trrlxtc[n]);
    }
  }
  if (params->trrest) {
    if (strlen(params->trrest[n])) {
      fprintf(op, "TRACER%1.1d.reset_file      %s\n",
	      n, params->trrest[n]);
      fprintf(op, "TRACER%1.1d.reset_dt        %s\n",
	      n, otime(params->trrestdt[n], tag));
    }
  }
  fprintf(op, "\n");
}

/*-------------------------------------------------------------------*/
/* Writes a tracer list to ascii file                                */
/*-------------------------------------------------------------------*/
int tracer_write(parameters_t *params, FILE *op)
{
  int n, tn, ntr = 0;

  /* Only print salt and temp for RECOM and PRE_MARVL */
  if (params->roammode & (A_RECOM_R1|A_RECOM_R2) || params->runmode & PRE_MARVL) {
    int num_tr = 2;
    fprintf(op, "\n# Tracers\n");
    fprintf(op, "NTRACERS                %d\n\n", num_tr);
    n = 0;
    for (tn = 0; tn < params->ntr; tn++) {
      tracer_info_t *tracer = &params->trinfo_3d[tn];
      if ((strcmp(tracer->name, "salt") == 0) || (strcmp(tracer->name, "temp") == 0)) {
	int ta = tracer->advect;
	int td = tracer->diffuse;
	if (params->runmode & PRE_MARVL) {
	  tracer->advect  = 1;
	  tracer->diffuse = 1;
	} else {
	  /* Forces read from transport files */
	  tracer->advect  = 0;
	  tracer->diffuse = 0;
	}
	tracer_write_3d(master, op, tracer, n);
	tracer->advect = ta;
	tracer->diffuse = td;
	n++;
      }
    }
    return(n);
  }

  n = 0;
  for (tn = 0; tn < params->ntrS; tn++) {
    tracer_info_t *tracer = &params->trinfo_2d[tn];
    if (strlen(tracer->name)) n++;
  }

#if defined(HAVE_SEDIMENT_MODULE)
  if (params->do_sed && strlen(params->sed_vars)) {
    n += count_sed_classes(params->sed_vars);
    n += count_sed_2d();
  }
#endif
#if defined(HAVE_ECOLOGY_MODULE)
  if (params->do_eco && strlen(params->eco_vars)) {
    n += count_eco_3d(params->eco_vars);
    n += count_eco_2d(params->eco_vars);
  }
#endif

  fprintf(op, "\n# Tracers\n");
  fprintf(op, "NTRACERS                %d\n\n", n+params->ntr-params->atr);
  for (tn = params->atr; tn < params->ntr; tn++) {
    tracer_info_t *tracer = &params->trinfo_3d[tn];
    n = tn - params->atr;
    /* Call the 3D helper function */
    tracer_write_3d(master, op, tracer, n);
    ntr++;
  }
  n = params->ntr - params->atr;
  for (tn = 0; tn < params->ntrS; tn++) {
    tracer_info_t *tracer = &params->trinfo_2d[tn];
#if defined(HAVE_ECOLOGY_MODULE)
    if (strlen(tracer->name) && !is_eco_var(tracer->name)) {
#else
    if (strlen(tracer->name)) {
#endif
      /* Call the 2D helper function */
      tracer_write_2d(master, op, tracer, n);
      n++;
      ntr++;
    }
  }
  return(ntr);
}

/* END tracer_write()                                                */
/*-------------------------------------------------------------------*/

