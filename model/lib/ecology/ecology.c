/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/ecology.c
 *  
 *  Description:
 *  Top-level ecology functions
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: ecology.c 6224 2019-05-28 06:10:38Z riz008 $
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <emslogger.h>
#include "ems.h"
#include "ecology_internal.h"
#include "column.h"
#include "cell.h"
#include "parameter_info.h"
#include "eprocess.h"
#include "einterface.h"
#include "utils.h"
#include "ecology_version.h"
#include "externallibs.h"
#include "bio_opt.h"

#ifdef HAVE_OMP
#include "omp.h"
#endif

#define STRBUFSIZE 2048
#define DTOR 0.01745329

// defined in ecology interface
void einterface_log_error(ecology *e, void *model, int b);

/*
 * Prototypes
 */
extern prm_seterrfn_fn prm_seterrorfn;
extern prm_readstring_fn prm_readstring;
extern prm_readdouble_fn prm_readdouble;
extern prm_readint_fn prm_readint;
extern prm_getkey_fn prm_getkey;

/* Array of flag of "essential" tracers for adapt1() step control. */
int* essential = NULL;


/* Version information */
int get_ecology_major_vers(void)
{
  return(ECOLOGY_MAJOR_VERSION);
}

int get_ecology_minor_vers(void)
{
  return(ECOLOGY_MINOR_VERSION);
}

int get_ecology_patch_vers(void)
{
  return(ECOLOGY_PATCH_VERSION);
}

/** Reads processes from the process parameter file.
 * @param fname Process parameter file name
 * @param key Process group tag
 * @return Pointer to stringtable with process names
 */
static stringtable* read_process_group_from_file(char* fname, char* key)
{
    stringtable* st = NULL;
    int firstpos = -1;
    FILE *f = e_fopen(fname, "r");

    /*
     * read all process groups of `key' type
     */
    while (1) {
        char buf[STRBUFSIZE];
        char* s = NULL;
        char* end = NULL;

        if (!find_key(fname, f, key))
            return st;

        /*
         * read processes of the current group (if not already read)
         */
        while (fgets(buf, STRBUFSIZE, f) != 0 && (s = strchr(buf, '{')) == NULL);

        if (s == NULL)
            return st;

        if (firstpos == -1)
            firstpos = ftell(f);
        else if (firstpos == ftell(f))
            return st;          /* no more entries of this type */

        if (st == NULL)
            st = stringtable_create(key);

        s++;                    /* just after '{' */
        do {
            int fpos = ftell(f);

            if (!isemptyline(s)) {
                int len0 = strlen(buf);
                int len;

                /*
                 * remove junk in the head of the string
                 */
                while (isspace((int)s[0]))
                    s++;
                end = strchr(s, '}');
                if (end != NULL)
                    /*
                     * '}' found -> cut the string at '}'
                     */
                    (*end) = 0;
                /*
                 * remove junk at the tail of the string
                 */
                while ((len = strlen(s)) > 0 && isspace((int)s[len - 1]))
                    s[len - 1] = 0;
                if (strlen(s) > 0)
                    stringtable_add(st, s, -1);
                if (end != NULL)
                    /*
                     * We need to start reading from the middle of this line
                     * in this case. Move file pointer back.
                     */
                    fseek(f, fpos - (len0 - (end - buf) - 1), 0);
            }
            s = buf;
        } while (end == NULL && fgets(buf, STRBUFSIZE, f) != 0);
    }
    fclose(f);
}

/** Reads processes from the default process list
 * @param fname Process parameter file name
 * @param key Process group tag
 * @return Pointer to stringtable with process names
 */
static stringtable* read_process_group_from_defaults(char* fname, int type)
{
    stringtable* st = NULL;
    int i;
    const char **procs;
    int nprocs;

    /* Call out to the gateway function in process_defaults.c */
    if (get_eco_processes(fname, type, &procs, &nprocs))
      e_quit("ecology: Unknown default proceeses type '%s' specified\n", fname);
    
    for (i=0; i<nprocs; i++) {
      /* Create stringtable on the first iteration */
      if (st == NULL)
	st = stringtable_create(eprocess_type_tags[type]);
      
      /* Add to the stringtable */
      stringtable_add(st, (char *)procs[i], -1);
    }
    
    return(st);
}

/** Creates a list of processes based on a static list
 *
 */
static void init_default_process(ecology *e, int type,
				 int nprocs, const char *procs[])
{
  int i;
  eprocess *p;
  char buf[MAXSTRLEN];
  char *parens = "";

  /* Allocate memory for processes for this type */
  e->processes[type] = calloc(nprocs, sizeof(void*));
  
  /* Now loop over and create each one */
  for (i=0; i<nprocs; i++) {
    char *s = strstr(procs[i], "(");
    if (s == NULL)
      parens = "()";
    sprintf(buf, "%s%s", procs[i], parens);
    e->processes[type][i] = eprocess_create(e, type, buf);
  }
  e->npr[type] = nprocs;
}

/** Creates list of processes in structure `ecology'.
 * @param e Pointer to structure `ecology' to hold the processes.
 * @param f Process parameter file.
 */
void create_processes(ecology* e, char* fname)
{
    FILE* f = NULL;
    int i,type,use_defs;
    
    e->npr       = calloc(N_EPROCESS_TYPES, sizeof(int));
    e->processes = calloc(N_EPROCESS_TYPES, sizeof(void*));
    
    /*
     * Check defaults vs filenams
     */
    use_defs = 1;
    if ( fname != NULL && (f = fopen(fname, "r")) != NULL ) {
      fclose(f);
      use_defs = 0;
    }

    /*
     * Loop through and initialise all processes
     */
    for (type = PT_WC; type <= PT_EPI; ++type) {
        char* key = eprocess_type_tags[type];
        stringtable* st;

	if (use_defs)
	  st = read_process_group_from_defaults(fname, type);
	else
	  st = read_process_group_from_file(fname, key);

        if (st == NULL)
            continue;

        e->npr[type] = st->n;
        if (st->n > 0)
            e->processes[type] = calloc(st->n, sizeof(void*));

        for (i = 0; i < st->n; ++i)
            e->processes[type][i] = eprocess_create(e, type, st->se[i]->s);

        stringtable_destroy(st);
    }

    /*
     * Now loop through and run post_init
     */
    for (type = PT_WC; type <= PT_EPI; ++type) {
      for (i = 0; i < e->npr[type]; ++i) {
	eprocess* p = e->processes[type][i];
	if (p->postinit != NULL) {
	  emstag(LDEBUG,"eco:ecology:create_processes", " %s: %s: postinit:\n", eprocess_type_tags[type], p->fullname);
	  eco_write_setup(e, "Calling %s_postinit\n", p->name);	  
	  p->postinit(p);
	  eco_write_setup(e, "\n");
	}
      }
    }
}

/** Destroys existing processes.
 * @param e Pointer to structure `ecology'
 */
void destroy_processes(ecology* e)
{
    int i, j;

    for (j = 0; j < N_EPROCESS_TYPES; ++j) {
        for (i = 0; i < e->npr[j]; ++i) {
            eprocess* p = e->processes[j][i];

            eprocess_destroy(p);
        }
        if (e->npr[j] > 0)
            free(e->processes[j]);
    }
    free(e->processes);
    e->processes = NULL;
    free(e->npr);
    e->npr = 0;
}

/** Prints processes included from the process parameter file.
 * @param e Ecology
 */
static void print_processes(ecology* e)
{
    int type;

    emstag(LINFO,"eco:ecology:print_processes", " processes:");
    for (type = 0; type < N_EPROCESS_TYPES; ++type) {
        int n = e->npr[type];
        int i;

        if (n > 0) {
            emslog(LINFO, "  %s  {", eprocess_type_tags[type]);
            for (i = 0; i < n; ++i)
                emslog(LINFO, "    %s", e->processes[type][i]->fullname);
            emslog(LINFO, "  }");
        }
    }
    emslog_flush();
}

/** Prints unmapped tracers, epivariables and parameters.
 * @param e Ecology
 */
static void print_unmapped(ecology* e)
{
    int i;
    int n;
    stringtable* st;

    emslog(LDEBUG, "ecology: unused tracers:\n");
    st = e->tracers;
    n = st->n;
    for (i = 0; i < n; ++i)
        if (st->se[i]->naccess == 0)
            stringtable_printentry_fp(st, i, emslog_fp(LINFO));

    emslog(LDEBUG, "ecology: unused epivariables:\n");
    st = e->epis;
    n = st->n;
    for (i = 0; i < n; ++i)
        if (st->se[i]->naccess == 0)
            stringtable_printentry_fp(st, i,emslog_fp(LINFO));

    emslog(LDEBUG, "ecology: unused parameters:\n");
    st = e->prms;
    n = st->n;
    for (i = 0; i < n; ++i)
        if (st->se[i]->naccess == 0)
            stringtable_printentry_fp(st, i,emslog_fp(LINFO));
}

#if (NCPU > 1)
void pth_init(ecology* e)
{
    pthstuff* pth = malloc(sizeof(pthstuff));
    int i;

    e->pth = pth;
    pthread_attr_init(&pth->thread_attr);

    /*
     * The threads are run as detached threads here. This means that once
     * started, they are on their own, no formal synchronization (joining)
     * point is required. As one of the threads terminates, another thread
     * is started. The advantage of this approach is that only NCPU treads
     * for running column_step() are required for maximum effectiveness
     * (speed); above that, the threads number does not matter. An
     * alternative approach I tried is to run threads as joinable. This
     * would mean to launch them in pools and waiting for each pool to
     * terminate. It is a bit more simple from programming point of view at
     * a minimum level but is not supposed to go for maximum efficiency.
     */
    pthread_attr_setdetachstate(&pth->thread_attr, PTHREAD_CREATE_DETACHED);
    pthread_mutex_init(&pth->stats_lock, NULL);
    pthread_mutex_init(&pth->avail_lock, NULL);
    pthread_cond_init(&pth->have_avail, NULL);
    pth->navail = NCPU;
    for (i = 0; i < NCPU; ++i)
        pth->avail[i] = 1;
}
#endif

ecology* ecology_create()
{
    ecology* e = malloc(sizeof(ecology));
    memset(e, 0, sizeof(ecology));

    e->model = NULL;
    e->verbose = VERBOSE_DEF;
    e->tscale = 1.0;
    e->prmfname = NULL;
    e->biofname = NULL;
    e->processfname = NULL;
    e->ntr = 0;
    e->tracerflags = NULL;
    e->tracerdiagn = NULL;
    e->tracernames = NULL;
    e->tracers = NULL;
    e->nepi = 0;
    e->epidiagn = NULL;
    e->epinames = NULL;
    e->epis = NULL;
    e->nprm = 0;
    e->pinfo = NULL;
    e->npr = NULL;
    e->processes = NULL;
    e->prms = NULL;
    e->cv_cell = NULL;
    e->cv_column = NULL;
    e->cv_model = NULL;
    e->mandatory_water = MANDATORY_WATER_DEF;
    e->mandatory_sediment = MANDATORY_SEDIMENT_DEF;
    e->t = NaN;
    e->jultime = NaN;
    e->dt = NaN;
    e->ncolumns = -1;
    e->nwclayers = -1;
    e->nsedlayers = -1;
    e->ncv = 0;
    e->cv = NULL;
    e->eint = dopri5;
    e->precision = NaN;
    e->h0 = NULL;
    e->check_nans = CHECK_NANS_DEF;
    e->check_negs = CHECK_NEGS_DEF;
    e->verify_diagn = VERIFY_DIAGN_EXTERN_DEF;
    e->quitfn = NULL;
    e->nstep = 0;
    e->modelstats.nwc = 0;
    e->modelstats.nws = 0;
    e->modelstats.nsc = 0;
    e->modelstats.nss = 0;
    e->modelstats.nec = 0;
    e->modelstats.nes = 0;
    e->stepstats.nwc = INT_MAX;
    e->stepstats.nws = INT_MAX;
    e->stepstats.nsc = INT_MAX;
    e->stepstats.nss = INT_MAX;
    e->stepstats.nec = INT_MAX;
    e->stepstats.nes = INT_MAX;
#if (NCPU > 1)
    e->multithreaded = MULTITHREADED_DEF;
    e->pth = NULL;
#endif
#ifdef HAVE_OMP
    e->omp_num_threads = OMP_NUM_THREADS_DEF;
#endif
    e->bio_opt = NULL;
    e->pre_build = 0;
    e->use_multi_sed = 0;
    e->eco_setup = NULL;
    e->eco_osetup = NULL;
    
    return e;
}

static void verify_diagnflags(int n, char* names[], int diagn[], void (*errorfn) (const char*, ...))
{
    int i, j;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < NDIAGNFLAGS; ++j)
            if (strcmp(names[i], diagnflags[j].name) == 0)
                break;

        if (j == NDIAGNFLAGS)
        {/*UR 11/2005 changed to warning, model exit seemed inappropriate */
          emstag(LWARN,"ecology:ecology:verify_diagnflags","tracer \"%s\" not found in \"diagnflags\" list from ecology/allprocesses.c\n", names[i]);
          /*errorfn("tracer \"%s\" not found in \"diagnflags\" list from ecology/allprocesses.c\n", names[i]);*/
        }
        if (diagn[i] != diagnflags[j].diagn)
        {
          emstag(LWARN,"ecology:ecology:verify_diagnflags","tracer \"%s\" has \"diagn\" flag set to %d; according to \"diagnflags\" list from ecology/allprocesses.c it should be %d", names[i], diagn[i], diagnflags[j].diagn);
          /*errorfn("tracer \"%s\" has \"diagn\" flag set to %d; according to \"diagnflags\" list from ecology/allprocesses.c it should be %d\n", names[i], diagn[i], diagnflags[j].diagn);*/
        }
    }
}


static void ecology_initessential(ecology* e)
{
    int i;

    assert(essential == NULL);

    essential = calloc(e->ntr, sizeof(int));
    for (i = 0; i < e->ntr; ++i)
        if (einterface_gettracerdiagnflag(e->model, e->tracernames[i]) == 0)
            essential[i] = 1;
}

/*
 * Collects the RSR wavelengths from the model
 */
void ecology_find_rsr_waves(ecology *e)
{
  int n, rsr_count;
  int *rtns;

  rsr_count = einterface_get_num_rsr_tracers(e->model);
  /* Check if found any */
  if (!rsr_count) return;

  /* Allocate temp buffer */
  rtns = i_alloc_1d(rsr_count);

  /* Grab rsr tracer indicies */
  einterface_get_rsr_tracers(e->model, rtns);

  e->num_rsr_waves = rsr_count;
  e->rsr_waves = d_alloc_1d(rsr_count);

  /* Get the names and tease out the actual wavelengths */
  for (n = 0; n < rsr_count; n++) {
    char *trname = einterface_get2Dtracername(e->model, rtns[n]);
    if (strncmp(trname, "R_", 2) == 0) {
      int w2;
      if (sscanf(trname, "R_%d", &w2) != 1)
	e_quit("Unable to retrieve wavelength from '%s'\n", trname);

      /* Cache this wavelenth */
      e->rsr_waves[n] = (double)w2;

      /* Add to epi stringtable */
      // stringtable_add(e->epis, trname, -1);
    }
  }
  i_free_1d(rtns);
}

/** Ecology pre-constructor
 *
 */
void *ecology_pre_build(char *eco_vars, char *eco_defs, FILE *fp)
{
  int i;
  ecology *e = ecology_create();
  char buf[MAXLINELEN];

  /* One-off sanity check */
  einterface_check_default_tracers();

  /* Set flag */
  e->pre_build = 1;

  /* Get process name */
  if (prm_readstring(fp, "processfname", buf))
    e->processfname = strdup(buf);
  else
    e_quit("processfname not found!\n");
  
  /* Get bio parameters name */
  if (prm_readstring(fp, "biofname", buf))
    e->biofname = strdup(buf);
  else
    e_quit("biofname not found!\n");
  
  /* Pass in null for process file name */
  read_parameter_info(e->biofname, fp, NULL, quit, warn, quiet, &e->nprm, &e->pinfo);

  // Set up _or_add functions so that we don't hit the host model
  e->find_index = find_index_or_add;
  e->try_index  = try_index_or_add;

  // Set quit function, could be warn
  e->quitfn = e_quit;

  /*
   * Create all the stringtables needed to create processes 
   */
  // Params
  e->prms = stringtable_create("parameters");
  for (i = 0; i < e->nprm; ++i)
    stringtable_add(e->prms, e->pinfo[i].name, i);

  // Tracers
  e->tracers = stringtable_create("tracers");
  e->epis    = stringtable_create("epivariables");

  // Common variables
  e->cv_cell   = stringtable_create("cell variables");
  e->cv_column = stringtable_create("column variables");
  e->cv_model  = stringtable_create("model variables");

  /* 
   * Go ahead and create the processes
   */
  create_processes(e, e->processfname);

  /* Write the processes to file */
  write_eco_process(e);

  /*
   * Build up a space separated string
   */
  // 3D
  stringtable_sort(e->tracers);
  for (i=0; i<e->tracers->n; i++) {
    char *str = e->tracers->se[i]->s;
    sprintf(eco_vars, "%s ", str);
    eco_vars += strlen(str)+1;
  }
  // 2D
  stringtable_sort(e->epis);
  for (i=0; i<e->epis->n; i++) {
    char *str = e->epis->se[i]->s;
    sprintf(eco_vars, "%s ", str);
    eco_vars += strlen(str)+1;
  }

  // Terminate string after backing up space
  *(--eco_vars) = '\0';

  // We're done
  return(e);
}


/** Ecology constructor.
 * @param model pointer to host model
 * @param prmfname Parameter file name
 * @return struct ecology
 */
ecology* ecology_build(void* model, char* prmfname)
{
    ecology* e = ecology_create();
    FILE* prmfile = NULL;
    char buf[MAXLINELEN];
    int i,j;

    e->model = model;

    e->verbose = einterface_getverbosity(model);

    {
        char* tunits = einterface_gettimeunits(model);

        if (tunits != NULL)
            e->tscale = timeunits2seconds(tunits);
    }

    prm_seterrorfn(e_quit);
    e->quitfn = e_quit;

    e->prmfname = strdup(prmfname);
    emstag(LTRACE,"eco:ecology:ecology_build", "prmfname = \"%s\"\n", prmfname);

    prmfile = e_fopen(prmfname, "r");

    prm_seterrorfn(quiet);
    e->internaltracers = INTERNALTRACERS_DEF;
    prm_readint(prmfile, "internal_tracers", &e->internaltracers);
    emstag(LTRACE,"eco:ecology:ecology_build","internal_tracers = %s\n", (e->internaltracers) ? "true" : "false");
    prm_seterrorfn(e_quit);

    /*
     * [FR] I think we need better documentation with a test case for
     * the internaltracers=1 case
     */
    e->find_index = (e->internaltracers) ? find_index_or_add : find_index;
    e->try_index  = (e->internaltracers) ? try_index_or_add  : try_index;

    /*
     * get tracers and epibenthic variables from the host model if required
     */
    if (!e->internaltracers) {
        /*
         * tracers
         */
        e->ntr = einterface_getntracers(model);
        if (e->ntr > 0) {
            e->tracerdiagn = malloc(e->ntr * sizeof(int));
            e->tracerflags = malloc(e->ntr * sizeof(int));
            e->tracernames = malloc(e->ntr * sizeof(char*));
            e->tracers = stringtable_create("tracers");
        }
        for (i = 0; i < e->ntr; ++i) {
            char* name = einterface_gettracername(model, i);

            e->tracerflags[i] = einterface_get_eco_flag(model, name);
            e->tracerdiagn[i] = einterface_gettracerdiagnflag(model, name);
            e->tracernames[i] = strdup(name);
            stringtable_add(e->tracers, name, i);
        }
        emstag(LINFO,"eco:ecology:build","Found %d tracers:\n", e->ntr);
        stringtable_print_fp(e->tracers,emslog_fp(LDEBUG));

        stringtable_sort(e->tracers);

        /*
         * epibenthic variables
         */
        e->nepi = einterface_getnepis(model);
        if (e->nepi > 0) {
            e->epidiagn = malloc(e->nepi * sizeof(int));
            e->epinames = malloc(e->nepi * sizeof(char*));
            e->epis = stringtable_create("epivariables");
        }
        for (i = 0; i < e->nepi; ++i) {
            char* name = einterface_getepiname(model, i);

            e->epidiagn[i] = einterface_getepidiagnflag(model, name);
            e->epinames[i] = strdup(name);
            stringtable_add(e->epis, name, i);
        } 
        emstag(LINFO,"eco:ecology:build","Found %d epibenthic variables:\n", e->nepi);
        stringtable_print_fp(e->epis,emslog_fp(LDEBUG));

        stringtable_sort(e->epis);
    } else {
        e->tracers = stringtable_create("tracers");
        e->epis = stringtable_create("epivariables");
    }

    /*
     * parameters
     */
    // Always allocate more
    if (prm_readstring(prmfile, "biofname", buf))
      e->biofname = strdup(buf);
    emstag(LINFO,"eco:ecology:ecology_build","biofname = \"%s\"\n", buf);
    read_parameter_info(e->biofname, prmfile, NULL, quit, warn, quiet, &e->nprm, &e->pinfo);
    e->prms = stringtable_create("parameters");
    for (i = 0; i < e->nprm; ++i)
        stringtable_add(e->prms, e->pinfo[i].name, i);
    emstag(LINFO,"eco:ecology:ecology_build"," %d parameters:\n", e->prms->n);
    stringtable_print_fp(e->prms,emslog_fp(LDEBUG));

    /*
     * common variables
     */
    e->cv_cell = stringtable_create("cell variables");
    e->cv_column = stringtable_create("column variables");
    e->cv_model = stringtable_create("model variables");

    /*
     * processes
     */
    if (prm_readstring(prmfile, "processfname", buf)) {
      e->processfname = strdup(buf);
      emstag(LINFO,"eco:ecology:ecology_build","processfname = \"%s\"\n", buf);
    } else {
      e_quit("processfname not found\n");
    }
    
    /* 
     * Open ecology setup file 
     *  - Do this after the initial create_processes call to only
     *    capture the events once (see call below)
     */
    e->eco_setup = e_fopen("ecology_setup.txt", "w");
    {
      /* Open file in the outputs area */
      char *opath = einterface_get_output_path();
      if (opath) {
	sprintf(buf, "%s/ecology_setup.txt", opath);
	e->eco_osetup = e_fopen(buf, "w");
      }
    }
	
    /* 
     * Initialise number of rsr wavelengths from the model 
     *   - needs to come before bio_opt_init and slot
     *     in between finalising the epi stringtable
     */
    e->bio_opt = NULL;
    create_processes(e, e->processfname);

    if (e->verbose) {
        print_processes(e);
        print_unmapped(e);      /* tracers, epis & parameters */
    }

    if (e->internaltracers) {
        e->ntr = e->tracers->n;
        e->nepi = e->epis->n;

        emstag(LINFO,"eco:ecology:ecology_build"," %d tracers:\n", e->ntr);
        stringtable_print_fp(e->tracers,emslog_fp(LDEBUG));
        emstag(LINFO,"eco:ecology:ecology_build"," %d epibenthic variables:\n", e->nepi);
        stringtable_print_fp(e->epis,emslog_fp(LDEBUG));

        if (e->ntr > 0) {
            e->tracerdiagn = malloc(e->ntr * sizeof(int));
            e->tracerflags = malloc(e->ntr * sizeof(int));
            e->tracernames = malloc(e->ntr * sizeof(char*));
        }
        for (i = 0; i < e->ntr; ++i) {
            e->tracernames[i] = strdup(stringtable_findindexstring(e->tracers, i));
            if(!einterface_tracername_exists(model, e->tracernames[i]))
            {
            	emstag(LPANIC,"Ecology:build","Attempting to retrieve tracer: %s which does NOT exist in host model!",e->tracernames[i]);
            	e->quitfn("Missing tracer: %s !",e->tracernames[i]);
            }	
            e->tracerdiagn[i] = einterface_gettracerdiagnflag(model, e->tracernames[i]);
            e->tracerflags[i] = einterface_get_eco_flag(model, e->tracernames[i]);
        }
        if (e->nepi > 0) {
            e->epidiagn = malloc(e->nepi * sizeof(int));
            e->epinames = malloc(e->nepi * sizeof(char*));
        }
        for (i = 0; i < e->nepi; ++i) {
            e->epinames[i] = strdup(stringtable_findindexstring(e->epis, i));
            if(!einterface_tracername_exists_epi(model,e->epinames[i]))
            {
            	emstag(LPANIC,"Ecology:build","Attempting to retrieve epi tracer: %s which does NOT exist in host model!", e->epinames[i]);
            	e->quitfn("Missing epi tracer: %s !", e->epinames[i]);
            }	
            e->epidiagn[i] = einterface_getepidiagnflag(model, e->epinames[i]);
        }
        stringtable_sort(e->tracers);
        stringtable_sort(e->epis);

        /*
         * do it again (silently):
         * need to know correct value of e->ntr to produce correct mappings
	 * Also the stringtable's have been sorted (above) so need
         * correct mappings
         */
        {
            int verbose = e->verbose;

            e->verbose = 0;
            destroy_processes(e);
            create_processes(e, e->processfname);
            e->verbose = verbose;
        }
    }

    e->verify_diagn = (e->internaltracers) ? VERIFY_DIAGN_INTERN_DEF : VERIFY_DIAGN_EXTERN_DEF;
    prm_seterrorfn(quiet);
    prm_readint(prmfile, "verify_diagn", &e->verify_diagn);
    emstag(LINFO, "ecology:ecology.c:ecology_build"," verify_diagn = %s\n", (e->verify_diagn) ? "true" : "false");

    if (e->verify_diagn) {
        verify_diagnflags(e->ntr, e->tracernames, e->tracerdiagn, e->quitfn);
        verify_diagnflags(e->nepi, e->epinames, e->epidiagn, e->quitfn);
    }else
      emstag(LINFO, "ecology:ecology.c:ecology_build"," verify_diagn ignored!");

    /*
     * water and sediment flags
     */
    prm_seterrorfn(quiet);
    if (prm_readstring(prmfile, "mandatory_water", buf))
        e->mandatory_water = is_true(buf);
    emstag(LINFO, "ecology:ecology.c:ecology_build"," mandatory_water = %s\n", (e->mandatory_water) ? "true" : "false");

    if (prm_readstring(prmfile, "mandatory_sediment", buf))
        e->mandatory_sediment = is_true(buf);
    emstag(LINFO, "ecology:ecology.c:ecology_build"," mandatory_sediment = %s\n", (e->mandatory_sediment) ? "true" : "false");

    /*
     * setting some geometry stuff
     */
    /*
     * deferred to ecology_step -FR
     * e->ncolumns = einterface_getnumbercolumns(model);
     */
    e->nwclayers = einterface_getnumberwclayers(model);
    e->nsedlayers = einterface_getnumbersedlayers(model);
    e->max_ncolumns = einterface_get_max_numbercolumns(model);
    
    /*
     * integrator
     */
    if (prm_readstring(prmfile, "integrator", buf)) {
        if (strcasecmp(buf, "dopri8") == 0) {
            e->eint = dopri8;
            e->precision = 1.0e-8;
        } else if (strcasecmp(buf, "dopri5") == 0) {
            e->eint = dopri5;
            e->precision = 1.0e-5;
        } else if (strcasecmp(buf, "adapt1") == 0) {
            e->eint = adapt1;
            e->precision = 1.0e-1;
            ecology_initessential(e);
        } else if (strcasecmp(buf, "adapt2") == 0) {
            e->eint = adapt2;
            e->precision = 1.0e-2;
        } else if (strcasecmp(buf, "euler1") == 0) {
            e->eint = euler1;
            e->precision = NaN;
        } else {
            e->eint = NULL;
            e->precision = NaN;
        }
    } else {
        e->eint = dopri5;
        e->precision = 1.0e-5;
    }
    prm_readdouble(prmfile, "integration_precision", &e->precision);

    if (e->eint == dopri8)
        emslog(LDEBUG, "ecology: integrator = dopri8\n");
    else if (e->eint == dopri5)
        emslog(LDEBUG, "ecology: integrator = dopri5\n");
    else if (e->eint == adapt1)
        emslog(LDEBUG, "ecology: integrator = adapt1\n");
    else if (e->eint == adapt2)
        emslog(LDEBUG, "ecology: integrator = adapt2\n");
    else if (e->eint == euler1)
        emslog(LDEBUG, "ecology: integrator = euler1\n");
    else
        e->quitfn("ecology: error: unknown integrator\n");
    emstag(LDEBUG, "ecology:ecology.c:ecology_build"," integration precision = %.3g\n", e->precision);

    e->h0 = malloc(e->max_ncolumns * (e->nwclayers + e->nsedlayers) * sizeof(double));
    prm_readint(prmfile, "check_nans", &e->check_nans);
    emstag(LINFO, "ecology:ecology.c:ecology_build"," check_nans = %s\n", (e->check_nans) ? "true" : "false");
    prm_readint(prmfile, "check_negs", &e->check_negs);
    emstag(LINFO, "ecology:ecology.c:ecology_build"," check_negs = %s\n", (e->check_negs) ? "true" : "false");

#if (NCPU > 1)
    prm_readint(prmfile, "multithreaded", &e->multithreaded);
#endif
#ifdef HAVE_OMP
    prm_readint(prmfile, "eco_omp_num_threads", &e->omp_num_threads);
    /*
     * Check openmp consistency
     */
    if ( (e->omp_num_threads != OMP_NUM_THREADS_DEF) && 
	 (einterface_transport_mode() == 0)              ) {
      // Running ecology with OMP is currently only supported for the
      // transport mode for now
      e_quit("ecology : eco_omp_num_threads is only supported for the transport mode");
    }
#endif

    fclose(prmfile);

    /* Initialise bio optical properties, if needed */
    // e->bio_opt = NULL;
    // if (process_present(e, PT_WC, "light_spectral_wc")) {
    // e->bio_opt = bio_opt_init(e);
      // }

    e->ncv = e->cv_model->n;
    if (e->ncv > 0)
        e->cv = malloc(e->ncv * sizeof(double));
    for (i = 0; i < e->ncv; ++i)
        e->cv[i] = NaN;         /* must be initialized by some process */

    /*
     * Allocate and initialise common model variables
     */
    e->ncv = e->cv_model->n;
    if (e->ncv > 0) {
      /* Allocate column variable pointer array */
      e->cv = malloc(e->ncv * sizeof(double));
    }
    
    if (einterface_getquitfn() != NULL)
      e->quitfn = einterface_getquitfn();
    e->nstep = 0;

    /*
     * custom part
     */
    einterface_ecologyinit(model, e);

    /* Close setup file, use the runlog from now on */
    fclose(e->eco_setup);
    if (e->eco_osetup)
      fclose(e->eco_osetup);
    e->eco_setup  = NULL;
    e->eco_osetup = NULL;

#if (NCPU > 1)
    if (e->multithreaded)
        pth_init(e);
#endif

    return e;
}

/** Ecology destructor.
 * @param e pointer to ecology
 */
void ecology_destroy(ecology* e)
{
    /* int verbose = e->verbose; */
    int i;

    if (is_log_enabled(LINFO) && e->nstep > 0)
    {
      ecology_printstats(e, stderr);
      ecology_printstats(e, emslog_fp(LINFO));
    }

    if (e->prmfname != NULL)
      free(e->prmfname);
    if (e->biofname != NULL)
      free(e->biofname);
    // Always allocated
    free(e->processfname);
    if (e->ntr > 0) {
        free(e->tracerflags);
        free(e->tracerdiagn);
        for (i = 0; i < e->ntr; ++i)
            free(e->tracernames[i]);
        free(e->tracernames);
        stringtable_destroy(e->tracers);
    }
    // Free bio-optical arrays
    if (e->bio_opt != NULL) {
      bio_opt_free(e->bio_opt);
      free(e->bio_opt);
      e->bio_opt = NULL;
    }
    if (e->nepi > 0) {
        free(e->epidiagn);
        for (i = 0; i < e->nepi; ++i)
            free(e->epinames[i]);
        free(e->epinames);
        stringtable_destroy(e->epis);
    }
    clear_parameter_info(&e->nprm, &e->pinfo);
    stringtable_destroy(e->prms);
    stringtable_destroy(e->cv_cell);
    stringtable_destroy(e->cv_column);
    stringtable_destroy(e->cv_model);
    if (e->ncv > 0)
      free(e->cv);

    destroy_processes(e);
    if (essential != NULL) {
        free(essential);
        essential = NULL;
    }
#if (NCPU > 1)
    if (e->multithreaded)
        free(e->pth);
#endif
    if (e->h0 != NULL)
      free(e->h0);
    free(e);
    emstag(LINFO,"eco:ecology:ecology_destroy"," destroyed\n");
}

#if (NCPU > 1)
static void pth_finish(pthstuff* pth)
{
    pthread_mutex_lock(&pth->avail_lock);
    while (pth->navail != NCPU)
        pthread_cond_wait(&pth->have_avail, &pth->avail_lock);
    pthread_mutex_unlock(&pth->avail_lock);
}

static void pth_free(pthstuff* pth, int thread_index)
{
    pthread_mutex_lock(&pth->avail_lock);
    pth->navail++;
    pth->avail[thread_index] = 1;
    pthread_cond_signal(&pth->have_avail);
    pthread_mutex_unlock(&pth->avail_lock);
}

static void* column_step_pth(void* p)
{
    pthargs* args = (pthargs*) p;
    column* col = column_create(args->e, args->column_index);

    column_step(col);
    column_destroy(col);

    pth_free(args->e->pth, args->thread_index);

    return NULL;
}
#endif

/*
 * Perform the column steps
 *
 */
 void column_do_steps(ecology* e)
 {
   int i;
  for (i = 0; i < e->ncolumns; ++i) {
    column* col;

    if (einterface_isboundarycolumn(e->model, i))
        continue;
    /*
     * We create columns at each ecology step rather than once because
     * this results in dynamic mapping, taking account of cell wetting or
     * drying, and because this is not too expensive compared with
     * integration.
     */
    col = column_create(e, i);
    column_step(col);
    column_destroy(col);
    }

 }

/*
 * Ecology's nstep
 */
int eco_get_nstep(ecology *e)
{
  return(e->nstep);
}

/** Performs one step for the whole model.
 * @param e pointer to ecology
 */
void ecology_step(ecology* e, double dt)
{
#if (NCPU > 1)
    pthstuff* pth = e->pth;
#endif
    int i,j;

    /*
     * Now that we are using vca2, we cannot precompute the number or
     * columns until now -FR
     */
    e->ncolumns = einterface_getnumbercolumns(e->model);

    e->t = einterface_getmodeltime(e->model) * e->tscale;
    e->dt = dt * e->tscale;
    if (e->nstep == 0) {
        int ncells = e->max_ncolumns * (e->nwclayers + e->nsedlayers);

        for (i = 0; i < ncells; ++i)
            e->h0[i] = e->dt;
    }

    for (i = 0; i < e->ncv; ++i)
        e->cv[i] = NaN;
    
    e->stepstats.nwc = 0;
    e->stepstats.nws = 0;
    e->stepstats.nsc = 0;
    e->stepstats.nss = 0;
    e->stepstats.nec = 0;
    e->stepstats.nes = 0;
    
#if (NCPU == 1)
#ifdef HAVE_OMP
    omp_set_nested(1);
#pragma omp parallel for private(i) num_threads(e->omp_num_threads),schedule(dynamic,1000)
#endif
    for (i = 0; i < e->ncolumns; ++i) {
      column* col;
      /* 
       * Unlike sediments, we're not actually keying of this flag just yet
       * Note: The multi-threading issue is now fixed
       */
      // i_set_error(e->model, i, LNONE, NULL);
      if (einterface_isboundarycolumn(e->model, i))
	continue;
      
      /*
       * This is a bit of hack to log ecology error for columns that
       * already have a sediment error (implied)
       */
      if (i_get_error(e->model, i+1) != LNONE) {
	einterface_log_error(e, e->model, i);
	continue;
      }

      /*
       * We create columns at each ecology step rather than once because
       * this results in dynamic mapping, taking account of cell wetting or
       * drying, and because this is not too expensive compared with
       * the integration.
       */
      col = column_create(e, i);
      if (!column_step(col))
	einterface_log_error(e, e->model, col->b);
      
      column_destroy(col);
    }

#else
    if (e->multithreaded) {
        for (i = 0; i < e->ncolumns; ++i) {
            pthargs* args;
            int j;

            if (einterface_isboundarycolumn(e->model, i))
                continue;

            pthread_mutex_lock(&pth->avail_lock);
            if (pth->navail == 0)
                pthread_cond_wait(&pth->have_avail, &pth->avail_lock);
            pth->navail--;

            for (j = 0; j < NCPU; ++j)
                if (pth->avail[j] == 1)
                    break;

            assert(j < NCPU);
            pth->avail[j] = 0;
            pthread_mutex_unlock(&pth->avail_lock);

            args = &pth->args[j];
            args->e = e;
            args->thread_index = j;
            args->column_index = i;
            pthread_create(&pth->threads[j], &pth->thread_attr, column_step_pth, args);
        }

        pth_finish(pth);
    } else {
        for (i = 0; i < e->ncolumns; ++i) {
            column* col;

            if (einterface_isboundarycolumn(e->model, i))
                continue;
    /*
         * We create columns at each ecology step rather than once because
         * this results in dynamic mapping, taking account of cell wetting or
         * drying, and because this is not too expensive compared with
         * the integration.
         */
            col = column_create(e, i);
            column_step(col);
            column_destroy(col);
        }
    }
#endif

    e->modelstats.nwc += e->stepstats.nwc;
    e->modelstats.nws += e->stepstats.nws;
    e->modelstats.nsc += e->stepstats.nsc;
    e->modelstats.nss += e->stepstats.nss;
    e->modelstats.nec += e->stepstats.nec;
    e->modelstats.nes += e->stepstats.nes;

    e->nstep++;
}



/** Returns the number of tracers.
 * @param e Pointer to ecology
 * @return Number of tracers
 */
int ecology_getntracers(ecology* e)
{
    if (!e->internaltracers)
        e_quit("ecology_getntracers(): This request assumes that the ecological model operates in the mode when tracers are initialised from the ecological processes. However, currently the tracers are initialised from the external model. You can switch the mode by adding \"internal_tracers 1\" to \"%s\"\n", e->prmfname);

    return e->ntr;
}

/** Returns the number of epibenthic variables.
 * @param e Pointer to ecology
 * @return Number of epibenthic variables
 */
int ecology_getnepis(ecology* e)
{
    if (!e->internaltracers)
        e_quit("ecology_getnepi(): This request assumes that the ecological model operates in the mode when tracers are initialised from the ecological processes. However, currently the tracers are initialised from the external model. You can switch the mode by adding \"internal_tracers 1\" to \"%s\"\n", e->prmfname);

    return e->nepi;
}

/** Returns the name of a tracer.
 * @param e Pointer to ecology
 * @param i Tracer index
 * @return The name of a tracer
 */
char* ecology_gettracername(ecology* e, int i)
{
    if (!e->internaltracers)
        e_quit("ecology_gettracername(): This request assumes that the ecological model operates in the mode when tracers are initialised from the ecological processes. However, currently the tracers are initialised from the external model. You can switch the mode by adding \"internal_tracers 1\" to \"%s\"\n", e->prmfname);

    return e->tracernames[i];
}

/** Returns the name of a epibenthic variable.
 * @param e Pointer to ecology
 * @param i Index of the epibenthic variable
 * @return The name of a epibenthic variable
 */
char* ecology_getepiname(ecology* e, int i)
{
    if (!e->internaltracers)
        e_quit("ecology_getepiname(): This request assumes that the ecological model operates in the mode when tracers are initialised from the ecological processes. However, currently the tracers are initialised from the external model. You can switch the mode by adding \"internal_tracers 1\" to \"%s\"\n", e->prmfname);

    return e->epinames[i];
}

/** Returns the value of diagnostic flag for a tracer or epibenthic variable.
 * @param name Variable name
 * @return The value of diagnostic flag; -1 if the name not found
 */
int ecology_getdiagnflag(char* name)
{
    int i;

    for (i = 0; i < NDIAGNFLAGS; ++i)
        if (strcmp(name, diagnflags[i].name) == 0)
            return diagnflags[i].diagn;

    return -1;
}


/** Define a function  for the processes to verify if a certain
 * process has been initialised it might depend on
 *
 * @aaram e - the ecology struct
 * @param type - the type of process [wc, epi, sed]
 * @param fullname - the assigned name of the process (as in allprocesses.c)
 * @return flag if invoked '1' or not '0'
 */
int process_present(ecology* e, int type,char* fullname)
{
  int i;

  if(fullname == NULL)
  	return 0;

  emstag(LTRACE,"eco:ecology:process_present","Testing process name: %s", fullname);
  for (i = 0; i < e->npr[type]; ++i)
  {
    eprocess* p = e->processes[type][i];
    if(strncmp(p->fullname, fullname, strlen(fullname)) == 0)
      return 1;
  }
  emstag(LTRACE,"eco:ecology:process_present","Not found process name: %s", fullname);
  
  return 0;
}

/*------------------------------------------------------------------*/
/* Writes a ecology transport configuration to file                 */
/*------------------------------------------------------------------*/
void trans_write_eco(char *eco_vars, char *eco_defs,
		     double ecodt, ecology *e, FILE *fp)
{
  fprintf(fp, "############################################################################\n");
  fprintf(fp, "# Ecology\n");
  if (strlen(eco_defs))
    fprintf(fp,"\nECO_VARS_ATTS   %s\n", eco_defs);
  fprintf(fp, "DO_ECOLOGY    YES\n");
  fprintf(fp, "ECOLOGY_DT    %d secs\n", (int)ecodt); // In seconds
  fprintf(fp, "\n");
  fprintf(fp, "# Internal tracers should always 1 for SHOC\n");
  fprintf(fp, "internal_tracers %d\n", e->internaltracers);
  fprintf(fp, "\n");
  fprintf(fp, "# Name of bio parameter file - empty for defaults\n");
  fprintf(fp, "biofname\t\t%s\n", e->biofname);
  fprintf(fp, "\n");
  fprintf(fp, "# Name of process file or defaults tag\n");
  fprintf(fp, "processfname\t\t%s\n", e->processfname);
  fprintf(fp, "\n");
  fprintf(fp, "# Verify diagn flag\n");
  fprintf(fp, "verfify_diagn %d\n", e->verify_diagn);
  fprintf(fp, "\n");
  fprintf(fp, "# Whether or not water and sediment cells are required to be empty\n");
  fprintf(fp, "mandatory_water    %d\n", e->mandatory_water);
  fprintf(fp, "mandatory_sediment %d\n", e->mandatory_sediment);
  fprintf(fp, "\n");
  fprintf(fp, "# The Integrator in use and its precision\n");
  if (e->eint == dopri8)
    fprintf(fp, "integrator\t\t dopri8\n");
  else if (e->eint == dopri5)
    fprintf(fp, "integrator \t\tdopri5\n");
  else if (e->eint == adapt1)
    fprintf(fp, "integrator \t\tadapt1\n");
  else if (e->eint == adapt2)
    fprintf(fp, "integrator \t\tadapt2\n");
  else if (e->eint == euler1)
    fprintf(fp, "integrator euler1\n");
  fprintf(fp, "integration_precision \t%e\n", e->precision);
  fprintf(fp, "\n");
  fprintf(fp, "# Check for NaN/negative values in tracers\n");
  fprintf(fp, "check_nans %d\n", e->check_nans);
  fprintf(fp, "check_negs %d\n", e->check_negs);
#ifdef HAVE_OMP
  fprintf(fp, "eco_omp_num_threads %d\n", e->omp_num_threads);
#endif
  fprintf(fp, "\n");
  
  fprintf(fp, "\n############################################################################\n\n");
}


/** Calculates the solar zenith angle
 * Needs to be done only once each time ecology is called.
 *
 * No longer in use
 */
static double calc_zenith(ecology *e)
{
/* Calculate Solar Zenith Angle (zenith)
* CONVERTED FROM http://www.srrb.noaa.gov/highlights/sunrise/azel.html April 2007 by Barbara Robson
* Purpose: calculate solar position for the entered date, time and
*               location.  Results are reported in azimuth and elevation
*               (in degrees) solar zenith angle.
* Must not use daylight saving time
* Requires parameters:
*   hrsToGMT: -10 for AEST
*   refLatitude: latitude to assume in calculations
*   refLongitude: longitude to assume in calculations
*/
//  double earthRadvec, v;
    double T, eccentricity, m, sinm, sin2m, SunEqOfCentre, ObliquityCorrection, tananum, tanadenom;
    double seconds, MeanObliquityOfEcliptic, SunApparentLong, SunTrueLong, GeomMeanLongSun, tmp;
    double SunRtAscension, solarDec, eqTime, solarTimeFix, trueSolarTime,  csz;
    double exoatmElevation, te, omega, hourangle, refractionCorrection, zenith;
    double jultime = e->jultime;

		/* get mandetory parameters from parameters cache */ 
    double hrsToGMT = get_parameter_value(e, "hrsToGMT");
    double latitude = get_parameter_value(e, "refLatitude");
    double longitude = get_parameter_value(e, "refLongitude");

    if ((latitude>=-90)&&(latitude<-89.8)) {
       latitude = -89.8;
    }
    if ((latitude<=90) && (latitude>89.8)) {
       latitude = 89.8;
    }
    assert ((hrsToGMT<=12) && (hrsToGMT>=-12.5));

    T = (jultime - 2451545.0)/36525.0; // Centuries since J2000.0
    // Calculate the distance to the sun in AU
    m = DTOR*(357.52911 + T * (35999.05029 - 0.0001537 * T)); // Geometric Mean Anomaly of the Sun in radians
    sinm = sin(m);
    sin2m = sin(m+m);
    SunEqOfCentre = sinm * (1.914602 - T * (0.004817 + 0.000014 * T)) + sin2m * (0.019993 - 0.000101 * T) + sin(3*m) * 0.000289; // Sun Equation of Centre, in degrees
    eccentricity = 0.016708634 - T * (0.000042037 + 0.0000001267 * T); // eccentricity of the Earth's orbit (unitless)

//    v = m * SunEqOfCentre; // True anomaly of the sun in degrees
//    earthRadVec = (1.000001018 * (1 - eccentricity * eccentricity)) / (1 + eccentricity * cos(v*DTOR)); // Distance to the sun/sun radius vector in AU

// Calculate the right ascension of the sun in degrees
    seconds = 21.448 - T*(46.8150 + T*(0.00059 - T*(0.001813)));
    MeanObliquityOfEcliptic = 23.0 + (26.0 + (seconds/60.0))/60.0;
    omega = (125.04 - 1934.136 * T) * DTOR;
    ObliquityCorrection = (MeanObliquityOfEcliptic + 0.00256 * cos(omega))*DTOR; // Corrected obliquity of the ecliptic in radians
    GeomMeanLongSun = 280.46646 + T * (36000.76983 + 0.0003032 * T);
    while(GeomMeanLongSun > 360.0) {
    	GeomMeanLongSun = GeomMeanLongSun - 360.0;
    }
    while(GeomMeanLongSun < 0.0) {
      GeomMeanLongSun = GeomMeanLongSun + 360.0;
    }
    GeomMeanLongSun = GeomMeanLongSun * DTOR; // Geometric mean Longitude of the Sun in radians
    SunTrueLong = GeomMeanLongSun/DTOR + SunEqOfCentre; // True longitude of the sun in degrees
    SunApparentLong = (SunTrueLong - 0.00569 - 0.00478 * sin(omega))*DTOR; // Apparent longitude of the sun in radians
    tananum = (cos(ObliquityCorrection) * sin(SunApparentLong));
    tanadenom = (cos(SunApparentLong));
    SunRtAscension = atan2(tananum, tanadenom)/DTOR; // right ascension of the sun in degrees

// Calculate the declination of the sun in radians
    solarDec = asin(sin(ObliquityCorrection) * sin(SunApparentLong)); // Sun's Declination in radians
// Calculate the equation of time in minutes of time
    tmp = tan(ObliquityCorrection / 2.0);
    tmp *= tmp;
    eqTime = tmp * sin(2.0*GeomMeanLongSun) -2.0 * eccentricity * sinm + 4.0 * eccentricity * tmp * sinm * cos(2.0 * GeomMeanLongSun) - 0.5 * tmp * tmp * sin(4.0 * GeomMeanLongSun) - 1.25 * eccentricity * eccentricity * sin2m;
    eqTime = eqTime/DTOR*4.0; // in minutes of time
/* Original had these rounded off.  Not sure we need to. * /
    eqTime = (floor(100*eqTime))/100;
    solarDec = (floor(100*solarDec/DTOR))/100*DTOR;
*/
    solarTimeFix = eqTime - 4.0 * longitude + 60.0 * hrsToGMT;
    trueSolarTime = (jultime+0.5 - floor(jultime+0.5))*24*60 + solarTimeFix; // in minutes

    while (trueSolarTime>1440) trueSolarTime -= 1440;
    hourangle = trueSolarTime / 4.0 - 180.0; // in degrees
    while (hourangle <-180) hourangle += 360.0; // in degrees
    csz = sin(latitude*DTOR) * sin(solarDec) + cos(latitude*DTOR) * cos(solarDec) * cos(hourangle*DTOR);
    if (csz>1.0)
      csz = 1.0;
    else if (csz<-1.0)
      csz = -1.0;
    zenith = acos(csz)/DTOR; // uncorrected solar zenith angle in degrees
    exoatmElevation = 90 - zenith; // in degrees
    if (exoatmElevation > 85.0)
      refractionCorrection = 0.0;
    else {
      te = tan(exoatmElevation*DTOR);
      if (exoatmElevation > 5.0)
        refractionCorrection = 58.1 / te - 0.07 / (te*te*te) + 0.000086 / te*te*te*te*te;
      else if (exoatmElevation > 0.575)
        refractionCorrection = 1735.0 + exoatmElevation * (-518.2 + exoatmElevation * (103.4 + exoatmElevation * (-12.79 + exoatmElevation * 0.711) ));
      else
        refractionCorrection = -20.774 / te;
      refractionCorrection = refractionCorrection / 3600.0;
    }
    zenith = zenith - refractionCorrection; // still in degrees
    if (zenith>90.0) zenith = 90; // sun below the horizon (astronomical twilight is at 108 degrees)
    return zenith;
}
    
    
/** Prints integration stats for this step to specified stream.
 * @param e Pointer to ecology
 * @param f Stream to print to
 */
void ecology_printstats(ecology* e, FILE* f)
{
    integration_stats* stats = &e->modelstats;

    fprintf(f, "ecology: %d steps\n", e->nstep);
    fprintf(f, "ecology:       No. of cells processed | Av. No. of integration steps\n");
    fprintf(f, "ecology: -----------------------------------------------------------\n");
    fprintf(f, "ecology: wc        %10d         |     %10.1f\n", stats->nwc, (double) stats->nws / (double) stats->nwc);
    fprintf(f, "ecology: sed       %10d         |     %10.1f\n", stats->nsc, (double) stats->nss / (double) stats->nsc);
    fprintf(f, "ecology: epi       %10d         |     %10.1f\n", stats->nec, (double) stats->nes / (double) stats->nec);
}

/** Prints integration stats for this step to specified stream.
 * @param e Pointer to ecology
 * @param f Stream to print to
 */
void ecology_printstepstats(ecology* e, FILE* f)
{
    integration_stats* stats = &e->stepstats;

    fprintf(f, "ecology:       No. of cells processed | Av. No. of integration steps\n");
    fprintf(f, "ecology: -----------------------------------------------------------\n");
    fprintf(f, "ecology: wc        %10d         |     %10.1f\n", stats->nwc, (double) stats->nws / (double) stats->nwc);
    fprintf(f, "ecology: sed       %10d         |     %10.1f\n", stats->nsc, (double) stats->nss / (double) stats->nsc);
    fprintf(f, "ecology: epi       %10d         |     %10.1f\n", stats->nec, (double) stats->nes / (double) stats->nec);
}

/* Write into ecology setup file */
void eco_write_setup(ecology *e, const char *str, ...)
{
  va_list args;

  /* Guard against pre_build */
  if (e->eco_setup == NULL) return;

  va_start(args, str);
  vfprintf(e->eco_setup, str, args);
  va_end(args);

  /* Write into outputs as well */
  if (e->eco_osetup) {
    va_start(args, str);
    vfprintf(e->eco_osetup, str, args);
    va_end(args);
  }
}

/* Public function */
void eco_set_omp_num_threads(ecology *e, int n)
{
#ifdef HAVE_OMP
  e->omp_num_threads = n;
#endif
}

