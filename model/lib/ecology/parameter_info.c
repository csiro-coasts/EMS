/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/parameter_info.c
 *  
 *  Description:
 *  Parameter info structure -- implementation
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: parameter_info.c 6678 2021-01-08 00:41:56Z bai155 $
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "externallibs.h"
#include "parameter_info.h"
#include "utils.h"
#include "constants.h"

/* Maximum number of values allowed per parameter */
#define MAXNUMVALUES (250)

static int get_eco_params(char name[], parameter_info *parameters[], int *nprm);
static void eco_params_auto(parameter_info *parametes[], int *nprm);
static void write_eco_params(char *name, int nparams, parameter_info parameters[]);

/** Reads entry for a given parameter attribute.
 * @param fp Parameter file
 * @param prefix Prefix, like "GRID0", or NULL
 * @param i Parameter index
 * @param tag Attribute
 * @param fn Error function
 * @param buffer Output string
 * @return 1 if successfull; 0 otherwise
 */
static int read_parameter_attribute(FILE* fp, 
				    const char* prefix,
				    char *name,
				    int i,
				    const char* tag, error_fn fn, char* buffer)
{
    char keybuf[MAXLINELEN];
    error_fn errfn_orig = prm_geterrorfn();
    int ret;

    char *key;
    if (i<0)
      key = prm_getkey(keybuf, prefix, "%s.%s", name, tag);
    else
      key = prm_getkey(keybuf, prefix, "PARAMETER%1.1d.%s", i, tag);

    prm_seterrorfn(fn);
    ret = prm_readstring(fp, key, buffer);
    prm_seterrorfn(errfn_orig);

    return ret;
}


/** Reads parameter file into parameter_info array.
 * @param biofname Name of bio param file
 * @param fp Parameter file
 * @param prefix Prefix, like "GRID0", or NULL
 * @param quitfn Quit error function
 * @param warnfn Warn error function
 * @param emptyfn Empty error function
 * @param nprm Pointer to number of parameters (output)
 * @param parameters Pointer to array of parameter_info (output)
 */
void read_parameter_info(char* biofname,
			 FILE* prmfp,
			 char* prefix, 
			 error_fn quitfn, 
			 error_fn warnfn, 
			 error_fn emptyfn, 
			 int* nprm, parameter_info** parameters)
{
    FILE* biofp = NULL;
    int n, i, sz;
    char buf[MAXLINELEN];
    char *fields[MAXSTRLEN * MAXNUMVALUES];

    prm_seterrorfn(quitfn);
    
    /*
     * Allocate max memory for array of parameter_info and fill with zeros 
     */
    sz = sizeof(parameter_info) * MAX_PARAMS;
    *parameters = (parameter_info*) malloc(sz);
    memset(*parameters, 0, sz);

    if (biofname == NULL)
      quitfn("biofname not found\n");
  
    // Try file
    biofp = fopen(biofname, "r");
      
    /* Read from file */
    if (biofp != NULL) {
      /*
       * Read the number of parameters 
       */
      prm_readint(biofp, prm_getkey(buf, prefix, "NPARAMETERS"), nprm);
      if (*nprm < 0)
        quitfn("read_parameter_info: NPARAMETERS < 0\n");
      else if (*nprm == 0)
        return;
      
      /*
       * Read parameter attributes for all parameters 
       */
      for (n = 0; n < *nprm; ++n) {
        parameter_info* prm = &(*parameters)[n];
	
        read_parameter_attribute(biofp, prefix, NULL, n, "name", quitfn, buf);
        prm->index = n;
        prm->name = strdup(buf);
        if (read_parameter_attribute(biofp, prefix, NULL, n, "desc", warnfn, buf))
	  prm->desc = strdup(buf);
        else
	  prm->desc = strdup("");
	if (read_parameter_attribute(biofp, prefix, NULL, n, "sym", warnfn, buf))
	  prm->sym = strdup(buf);
        else
	  prm->sym = strdup("");
        if (read_parameter_attribute(biofp, prefix, NULL, n, "units", warnfn, buf))
	  prm->units = strdup(buf);
        else
	  prm->units = strdup("1");
	if (read_parameter_attribute(biofp, prefix, NULL, n, "ref", warnfn, buf))
	  prm->ref = strdup(buf);
        else
	  prm->ref = strdup("Not attributed");
	
	if (read_parameter_attribute(biofp, prefix, NULL, n, "value", warnfn, buf)) {
	  prm->num_values = prm_parseline(buf, fields, MAXNUMVALUES);
	  if (prm->num_values) {
	    int  vstr  = 0; // length needed for stringvalue
	    char *pstr = NULL;
	    /* Allocate array */
	    prm->value = d_alloc_1d(prm->num_values);
	    /* And fill out values */
	    for (i=0; i<prm->num_values; i++) {
	      (void) str2double(fields[i], &prm->value[i]);
	      vstr += strlen(fields[i]) + 1;
	    }
	    pstr = prm->stringvalue = calloc(vstr+1, sizeof(char));
	    for (i=0; i<prm->num_values; i++) {
	      sprintf(pstr, "%s", fields[i]);
	      pstr += strlen(fields[i]);
	      /* Add spaces but not at the end */
	      if (i < prm->num_values-1) {
		sprintf(pstr, " ");
		pstr++;
	      }
	    }
	  }
	} else
	  prm->stringvalue = strdup("");

	if (read_parameter_attribute(biofp, prefix, NULL, n, "stderr", 
				     warnfn, buf))
	  str2double(buf, &prm->stderr);
	else
	  prm->stderr = 0.0;
      }
      fclose(biofp);
      
    } else {
      /*
       * Read defaults 
       *
       * Note: get_eco_params currently only has the one default set
       */
      if (get_eco_params(biofname, parameters, nprm))
	quitfn("read_parameter_info: Unable to allocate default eco parameters\n");

      /*
       * Now loop over and overwrite any param attributes
       * Guard against no prm-file. eg. mleco
       */
      if (prmfp != NULL) {
	for (n=0; n< *nprm; n++) {
	  parameter_info* prm = &(*parameters)[n];
	  char *name = prm->name;
	  
	  /* Description */
	  if (read_parameter_attribute(prmfp, prefix, name, -1, "desc", warnfn, buf))
	    prm->desc = strdup(buf);

	  /* Symbol */
	  if (read_parameter_attribute(prmfp, prefix, name, -1, "sym", warnfn, buf))
	    prm->sym = strdup(buf);
	  
	  /* Units */
	  if (read_parameter_attribute(prmfp, prefix, name, -1, "units", warnfn, buf))
	    prm->units = strdup(buf);
	  
	  /* Values */
	  if (read_parameter_attribute(prmfp, prefix, name, -1, "value", warnfn, buf)) {
	    prm->num_values = prm_parseline(buf, fields, MAXNUMVALUES);
	    /* De-allocate old */
	    if (prm->value)
	      d_free_1d(prm->value);
	    
	    /* Allocate array */
	    prm->value = d_alloc_1d(prm->num_values);
	    /* And fill out values */
	    for (i=0; i<prm->num_values; i++)
	      (void) str2double(fields[i], &prm->value[i]);
	    
	    if (buf != NULL)
	      prm->stringvalue = strdup(buf);
	  }
	}
	// Write out to file
	write_eco_params(biofname, *nprm, *parameters);
      }
    }
    days2seconds(*nprm, *parameters);
}


/** Clears parameter_info array.
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */
void clear_parameter_info(int* nprm, parameter_info* parameters[])
{
    int n;

    if ((*nprm == 0) || (*parameters == NULL))
        return;

    for (n = 0; n < *nprm; ++n) {
        parameter_info* pr = &(*parameters)[n];

        free(pr->name);
        free(pr->desc);
        free(pr->sym);
        free(pr->units);
        free(pr->ref);
	if (pr->stringvalue != NULL)
	  free(pr->stringvalue);
	d_free_1d(pr->value);
    }

    *nprm = 0;

    free(*parameters);
    *parameters = NULL;
}

/** Finds parameter index in a parameter_info array.
 * @param name Parameter name
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 * @return Index id successful; -1 otherwise
 */
int find_parameter_index(const char* name, int nprm, parameter_info pars[])
{
    int n;

    if (pars == NULL)
        return -1;
    for (n = 0; n < nprm; ++n)
        if (strcasecmp(pars[n].name, name) == 0)
            return pars[n].index;

    return -1;
}

/** Transfers parameters defined in d-1 to s-1.
 * @param nprm Number of parameters
 * @param parameters parameter_info array
 */
void days2seconds(int nprm, parameter_info parameters[])
{
    int i,j;

    for (i = 0; i < nprm; ++i) {
        parameter_info* p = &parameters[i];
        char* replace;

        if ((replace = strstr(p->units, "d-1")) != NULL) {
	  for (j=0; j<p->num_values; j++) {
            p->value[j] /= SEC_PER_DAY;
            strncpy(replace, "s-1", 3);
	  }
        }
    }
}


/*-------------------------------------------------------------------*/
/* Routine to read the sediment parameters                           */
/*-------------------------------------------------------------------*/
static int get_eco_params(char name[], parameter_info *parameters[], int *nprm)
{
  int i;

  struct {
    char *name;
    char *description;
    void (*init) (parameter_info *params[], int *nprm);
  } param_list[] = {
    {"standard","Standard ecology parameters",  eco_params_std},
    {"estuary", "Estuarine ecology parameters", eco_params_est},
    {"gbr4", "GBR4 ecology parameters", eco_params_gbr4},
    {"BGC2p0", "Parameters used in BGC 2.0", eco_params_bgc2p0},
    {"BGC3p1", "Parameters used in BGC 3.1", eco_params_bgc3p1},
    {"porewater", "Parameters used for porewater age", eco_params_porewater},
    {"TASSE1p0", "Parameters used for TASSE 1p0", eco_params_tasse1p0},
    {NULL, NULL, NULL}
  };
  void (*init) (parameter_info *params[], int *nprm)= NULL;

  if (strcmp(name, "auto") == 0) {
    eco_params_auto(parameters, nprm);
  } else {
    for (i = 0; (init == NULL) && param_list[i].name; ++i) {
      if (strcasecmp(name, param_list[i].name) == 0) {
	init = param_list[i].init;
      }
    }
    if (init != NULL) {
      init(parameters, nprm);
    } else
      return(1);
  }

  return(0);
}

/* END get_eco_params()                                              */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Dynamically prescribed ecology parameter specification            */
/*-------------------------------------------------------------------*/
static void eco_params_auto(parameter_info *parameters[], int *nprm)
{
  /* 
   * Prescribe values based on known system attributes 
   * TODO
   */
  e_quit("eco_params_auto: Not implemented as yet!\n");
}

/* END eco_params_auto()                                             */
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
/* Writes the ecological parameters to file                          */
/*-------------------------------------------------------------------*/
static void write_eco_params(char *name, int nparams, parameter_info parameters[])
{
  int n;
  char buf[MAXSTRLEN];
  FILE *fp;

  sprintf(buf, "bio_%s.prm", name);
  if ((fp = fopen(buf, "w")) == NULL)
    return;

  fprintf(fp, "NPARAMETERS      %d\n\n", nparams);
  for (n = 0; n < nparams; n++) {
    parameter_info* prm = &parameters[n];
    fprintf(fp, "PARAMETER%d.name    %s\n", n, prm->name);
    fprintf(fp, "PARAMETER%d.desc    %s\n", n, prm->desc);
    fprintf(fp, "PARAMETER%d.units   %s\n", n, prm->units);
    fprintf(fp, "PARAMETER%d.sym     %s\n", n, prm->sym);
    fprintf(fp, "PARAMETER%d.value   %s\n", n, prm->stringvalue);
    fprintf(fp, "PARAMETER%d.stderr  %f\n", n, prm->stderr);
    fprintf(fp, "PARAMETER%d.ref     %s\n\n", n, prm->ref);
  }

  fclose(fp);
}

/* END write_eco_params()                                            */
/*-------------------------------------------------------------------*/

