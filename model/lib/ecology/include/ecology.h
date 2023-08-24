/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/ecology.h
 *  
 *  Description:
 *  Ecology code -- header. This is supposed to be the only
 *  ecology header necessary to be included by external model
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology.h 7107 2022-04-20 11:44:57Z bai155 $
 *
 */

#if !defined(_ECOLOGY_H)
#include "ems.h"
#include "stringtable.h"
#include "ecology_version.h"

struct ecology;
typedef struct ecology ecology;

/* Ecology tracer attribute private dtat structure */
typedef struct {
  int type;              /* Private data type */
  char name[MAXSTRLEN];  /* Name of the default list */
  char trname[MAXSTRLEN]; /* Name of the tracer */
  int obc;               /* Open boundary condition */
  int flag;              /* General purpose flag */
  char optical_file[MAXSTRLEN];
  char absorp_name[MAXSTRLEN];
  char scatter_name[MAXSTRLEN];
  char backscat_name[MAXSTRLEN];
  char specresp3d_name[MAXSTRLEN];
  char abtance_name[MAXSTRLEN];
  char refltce_name[MAXSTRLEN];
  char trnmiss_name[MAXSTRLEN];
  char benreflt_name[MAXSTRLEN];
  char specresp2d_name[MAXSTRLEN];
  
} trinfo_priv_eco_t;

/** Ecology constructor.
 * @param model Pointer to host model
 * @param prmfname Parameter file name
 * @return struct ecology
 */
ecology* ecology_build(void* model, char* prmfname);

/** Ecology destructor.
 * @param e Pointer to ecology
 */
void ecology_destroy(ecology* e);

/** Performs one step for the whole model.
 * @param e Pointer to ecology
 */
void ecology_step(ecology* e, double dt);

/** Returns the number of tracers.
 * @param e Pointer to ecology
 * @return Number of tracers
 */
int ecology_getntracers(ecology* e);

/** Returns the number of epibenthic variables.
 * @param e Pointer to ecology
 * @return Number of epibenthic variables
 */
int ecology_getnepis(ecology* e);

/** Returns the name of a tracer.
 * @param e Pointer to ecology
 * @param i Tracer index
 * @return The name of a tracer
 */
char* ecology_gettracername(ecology* e, int i);

/** Returns the name of a epibenthic variable.
 * @param e Pointer to ecology
 * @param i Index of the epibenthic variable
 * @return The name of a epibenthic variable
 */
char* ecology_getepiname(ecology* e, int i);

/** Returns the value of diagnostic flag for a tracer or epibenthic variable.
 * @param name Variable name
 * @return The value of diagnostic flag; -1 if the name not found
 */
int ecology_getdiagnflag(char* name);

/** Prints integration stats for this step to specified stream.
 * @param e Pointer to ecology
 * @param f Stream to print to
 */
void ecology_printstepstats(ecology* e, FILE* f);

/** Prints integration stats for this run to specified stream.
 * @param e Pointer to ecology
 * @param f Stream to print to
 */
void ecology_printstats(ecology* e, FILE* f);

void eco_write_setup(ecology *e, const char *str, ...);

void ecology_find_rsr_waves(ecology *e);

void ecology_find_ed_waves(ecology *e);

/* Generic interface routines */

extern int i_is_geographic(void* hmodel);

#define _ECOLOGY_H
#endif
