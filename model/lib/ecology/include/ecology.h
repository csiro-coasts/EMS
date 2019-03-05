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
 *  $Id: ecology.h 6160 2019-03-05 04:35:12Z riz008 $
 *
 */

#if !defined(_ECOLOGY_H)
#include "stringtable.h"
#include "ecology_version.h"

struct ecology;
typedef struct ecology ecology;

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

#define _ECOLOGY_H
#endif
