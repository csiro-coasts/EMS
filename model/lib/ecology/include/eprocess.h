/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/eprocess.h
 *  
 *  Description:
 *  Process servicing code -- header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: eprocess.h 6770 2021-06-04 01:49:23Z riz008 $
 *
 */

#if !defined(_EPROCESS_H)

#include "declarations.h"

/* This enumeration defines major process types. The process types correspond
 * to sections in the process parameter file. They imply the following design
 * of the ecological code stepping.
 * 1. A model step consists of a number of independent column steps.
 * 2. A column step consists of:
 *    (a) Column processes
 *    (b) A number of independent cell steps
 * 3. A cell step consists of:
 *    (a) Preintegration
 *    (b) Integration
 *    (c) Postintegration
 *
 * Depending on the cell type, corresponding flux procedures are used.
 */
typedef enum {
    /*
     * generic; processes of this type can be entered as either PT_WC or
     * PT_SED 
     */
    PT_GEN = 0,
    PT_WC = 1,
    PT_SED = 2,
    PT_EPI = 3,
    PT_COL = 4,
} EPROCESSTYPE;

#define N_EPROCESS_TYPES (PT_COL + 1)

typedef struct {
    EPROCESSTYPE type;
    char* tag;
} eprocess_type_tag;

extern char* eprocess_type_tags[];

/* The process procedures.
 */
typedef void (*eprocess_initfn) (eprocess* p);
typedef void (*eprocess_calcfn) (eprocess* p, void* pp);

struct eprocess {
    ecology* ecology;
    EPROCESSTYPE type;
    char* name;
    char* fullname;             /* name + arguments */
    stringtable* prms;
    int delayed;                /* flag */

    eprocess_initfn init;
    eprocess_initfn postinit;
    eprocess_initfn destroy;
    eprocess_calcfn precalc;
    eprocess_calcfn calc;
    eprocess_calcfn postcalc;

    void* workspace;
};

/* A process entry in the global process list in "allprocesses.c".
 */
typedef struct {
    char* name;
    EPROCESSTYPE type;
    int nparams;
    int delayed;
    eprocess_initfn init;
    eprocess_initfn postinit;
    eprocess_initfn destroy;
    eprocess_calcfn precalc;
    eprocess_calcfn calc;
    eprocess_calcfn postcalc;
} eprocess_entry;

/* Global process list -- see "allprocesses.c".
 */
extern eprocess_entry eprocesslist[];
extern int NEPROCESSES;

typedef struct {
    char* name;
    int diagn;
} diagnflag_entry;

/* Global tracer list -- see "allprocesses.c".
 */
extern diagnflag_entry diagnflags[];
extern int NDIAGNFLAGS;

/** Creates and initialises process defined by the process entry string.
 * @param e Pointer to ecology
 * @param type Process type
 * @param entry Process entry string
 * @return Pointer to eprocess
 */
eprocess* eprocess_create(ecology* e, EPROCESSTYPE type, char* entry);

/** Process destructor.
 * @param p Process
 */
void eprocess_destroy(eprocess* p);

int get_eco_processes(char *name, int type, const char **procs[], int *nprocs);
void write_eco_process(ecology *e);


#define _EPROCESS_H
#endif
