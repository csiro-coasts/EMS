/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/ecology_internal.h
 *  
 *  Description:
 *  Ecology code -- "internal" header.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology_internal.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_ECOLOGY_INTERNAL_H)

#include <emslogger.h>
#include "ecology.h"
#include "einterface.h"
#include "declarations.h"
#include "integrator.h"
#include "ems.h"
#include "constants.h"
#include "bio_opt.h"

/*
 * Multi-threading is compiled in if NCPU > 1
 */
#if !defined(NCPU)
#define NCPU 1
#endif
#if (NCPU < 1)
#error NCPU < 1
#endif
#if !(NCPU == 1)
#include <pthread.h>
#include <signal.h>
#endif

#define VERBOSE_DEF 0
#define INTERNALTRACERS_DEF 1
#define MANDATORY_WATER_DEF 0
#define MANDATORY_SEDIMENT_DEF 1
#define CHECK_NANS_DEF 0
#define CHECK_NEGS_DEF 0
#define VERIFY_DIAGN_EXTERN_DEF 0
#define VERIFY_DIAGN_INTERN_DEF 1
#if (NCPU > 1)
#define MULTITHREADED_DEF 0
#endif

#ifdef HAVE_OMP
#define OMP_NUM_THREADS_DEF 1
#endif

#define ISDEBUG() (is_log_enabled(LDEBUG))

typedef struct {
    int nwc;                    /* number of watercolumn cells processed */
    int nws;                    /* total number of integration steps for wc
                                 * cells */
    int nsc;                    /* number of sediment cells processed */
    int nss;                    /* total number of integration steps for
                                 * sediment cells */
    int nec;                    /* number od epibenthic cells processed */
    int nes;                    /* total number of integration steps for
                                 * epibenthic cells */
} integration_stats;

#if (NCPU > 1)
typedef struct {
    ecology* e;
    int column_index;
    int thread_index;
} pthargs;

typedef struct {
    /*
     * Posix thread initialisation structure used to create a thread.
     */
    pthread_attr_t thread_attr;
    /*
     * Posix mutex used to lock ecostats.
     */
    pthread_mutex_t stats_lock;
    /*
     * Posix mutex used to lock thread counter.
     */
    pthread_mutex_t avail_lock;
    /*
     * Posix condition used to suspend thread requests if navail = 0. Also
     * used at waiting for all threads to finish at the end of
     * ecology_step().
     */
    pthread_cond_t have_avail;
    /*
     * Number of available threads; 0 <= navail <= NCPU
     */
    int navail;
    /*
     * Status of columns: not 0 - available, 0 - not available.
     */
    int avail[NCPU];
    /*
     * Threads to run column_step().
     */
    pthread_t threads[NCPU];
    /*
     * Information for launching a particular thread
     */
    pthargs args[NCPU];
} pthstuff;
#endif

typedef int (*findindex_fn) (stringtable* table, char* name, ecology* e);

struct ecology {
    /*
     * host model
     */
    void* model;

    int verbose;

    double tscale;

    char* prmfname;
    char* biofname;
    char* processfname;

    /*
     * Flag. If 0 (default), ecology->tracers are built using information from
     * the external model. If 1, ecology->tracers are built based on
     * requirements by the processes.
     */
    int internaltracers;
    findindex_fn find_index;
    findindex_fn try_index;

    int ntr;
    int* tracerdiagn;
    char** tracernames;         /* to get name by index */
    stringtable* tracers;       /* to get index by name */

    int nepi;
    int* epidiagn;
    char** epinames;            /* to get name by index */
    stringtable* epis;          /* to get index by name */

    int nprm;
    parameter_info* pinfo;
    stringtable* prms;

    /* tracer eco flag */
    int *tracerflags;
  
    /*
     * Processes. Kept in groups according to process type.
     */
    /*
     * Number of processes of each type. E.g., npr[PT_WC] is a number of
     * processes of PT_WC type.
     */
    int* npr;
    /*
     * Processes themselves. E.g., processes [PT_WC][0] gives (a pointer
     * to) the first process of PT_WC type.
     */
    eprocess*** processes;

    /*
     * list of common variables for all processes within a cell
     */
    stringtable* cv_cell;
    /*
     * list of common variables for all processes within a column
     */
    stringtable* cv_column;
    /*
     * list of common variables for all processes within the model
     */
    stringtable* cv_model;

    /*
     * flags; 1 by default
     */
    int mandatory_water;
    int mandatory_sediment;

    /*
     * current time
     */
    double t;
    
    /*
     * julien time in seconds
     */
    double jultime;
    
    /*
     * zenith in degres
     */
    double zenith;
    
    /*
     * current time step
     */
    double dt;
    /*
     * number of columns
     */
    int ncolumns;
    /*
     * max number of columns
     */
    int max_ncolumns;
    /*
     * number of wc layers including empty ones
     */
    int nwclayers;
    /*
     * number of sediment layers including empty ones
     */
    int nsedlayers;

    /*
     * process communication: common model variables
     */
    int ncv;
    double* cv;

    /*
     * some commonly used values
     */
    integrator eint;
    double precision;
    double* h0;                 /* initial step size */

    int check_nans;             /* flag */
    int check_negs;             /* flag */
    int verify_diagn;           /* flag */

    void (*quitfn) (const char* format, ...);

    /*
     * current step
     */
    int nstep;

    integration_stats modelstats;       /* summary of integration statistics
                                         * for the whole run */
    integration_stats stepstats;        /* integration statistics for this
                                         * step */

    /* Bio-optical properties */
    bio_opt_prop *bio_opt;

    /* This is used in RECOM/AUTOMATION */
    int pre_build;

    /* Wheter to use multiple sediment layers in epi processes*/
    int use_multi_sed;

    /* Reflectance info */
    int num_rsr_waves;
    double *rsr_waves;

    /* Ecology setup.txt file */
    FILE *eco_setup;

#if (NCPU > 1)
    int multithreaded;          /* flag */

    pthstuff* pth;
#endif

#ifdef HAVE_OMP
  int omp_num_threads;
#endif
};


/* Define a function  for the processes to verify if a certain
 * process has been initialised it might depend on
 *
 * @aaram e - the ecology struct
 * @param type - the type of process [wc, epi, sed]
 * @param fullname - the assigned name of the process (as in allprocesses.c)
 * @return flag if invoked '1' or not '0'
 */
int process_present(ecology* e, int type,char* fullname);

#define _ECOLOGY_INTERNAL_H
#endif
