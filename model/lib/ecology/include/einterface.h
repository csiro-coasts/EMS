/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/einterface.h
 *  
 *  Description:
 *  Interface functions to be provided by an external (host)
 *  model with which the ecology code is to be linked.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: einterface.h 7576 2024-05-30 03:47:52Z riz008 $
 *
 */

#if !defined(_INTERFACE_H)

/* Used for eco_flag */
#define ECO_NONE    0x0
#define ECO_NORESET 0x1

/* @param model Pointer to host model
 * @return Shallow copy of the time units string
 */
extern char* einterface_gettimeunits(void* model);

/* @param model Pointer to host model
 * @return Verbosity (0 = no log, 1 = min log, 2 = max log)
 */
extern int einterface_getverbosity(void* model);

/* @param model Pointer to host model
 * @return Number of tracers
 */
extern int einterface_getntracers(void* model);

/* @param model Pointer to host model
 * @param i Tracer index
 * @return Tracer name
 */
extern char* einterface_gettracername(void* model, int i);


/* @param model Pointer to host model
 * @param col column index
 */
extern void einterface_get_ij(void* model, int col, int* ij);


/* @param model Pointer to host model
 * @param name Tracer name
 * @return 1 if the tracer exists, 0 otherwise
 */
extern int einterface_tracername_exists(void* model, char*name);

/* @param model Pointer to host model
 * @param i Tracer index
 * @return Diagnostic flag value for the tracer
 *         (0 - non-diagnostic,
 *          1 - flux diagnostic
 *          2 - value diagnostic)
 */
extern int einterface_gettracerdiagnflag(void* model, char* name);

/* @param model Pointer to host model
 * @param i Tracer index
 * @return Particulate flag value for the tracer
 *         (0 - dissolved,
 *          1 - particulate)
 */
extern int einterface_gettracerparticflag(void* model, char* name);

/* @param model Pointer to host model
 * @return Number of epibenthic variables
 */
extern int einterface_getnepis(void* model);

/* @param model Pointer to host model
 * @param i Index of an epibenthic variable
 * @return Name of the variable
 */
extern char* einterface_getepiname(void* model, int i);

/* @param model Pointer to host model
 * @param i Index of an epibenthic variable
 * @return Diagnostic flag value for the variable
 *         (0 - non-diagnostic,
 *          1 - flux diagnostic
 *          2 - value diagnostic)
 */
extern int einterface_getepidiagnflag(void* model, char* name);

/* @param model Pointer to host model
 * @return Wrapping model time
 */
extern double einterface_getmodeltime(void* model);

/* @param model Pointer to host model
 * @return Number of watercolumn layers, including empty ones
 */
extern int einterface_getnumberwclayers(void* model);

/* @param model Pointer to host model
 * @return Number of sediment layers, including empty ones
 */
extern int einterface_getnumbersedlayers(void* model);

/* @param model Pointer to host model
 * @return number of wet columns only
 */
extern int einterface_getnumbercolumns(void* model);

/* @param model Pointer to host model
 * @return Number of columns, including boundary ones
 */
extern int einterface_get_max_numbercolumns(void* model);

/* @param model Pointer to host model
 * @param b column index
 * @return 1 if the column is a boundary one, 0 otherwhile
 */
extern int einterface_isboundarycolumn(void* model, int b);

/* @param model Pointer to host model
 * @param b column index
 * @return Index of the top non-empty watercolumn layer in the column
 */
extern int einterface_getwctopk(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Index of the bottom non-empty watercolumn layer in the column
 */
extern int einterface_getwcbotk(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Index of the top non-empty sediment layer in the column
 */
extern int einterface_getsedtopk(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Index of the bottom non-empty sediment layer in the column
 */
extern int einterface_getsedbotk(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @param k cell level, k_wc
 * @return Cell centred depth
 */
extern double einterface_getcellz(void *model, int b, int k);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of water cell thicknesses in the column
 *
 * Requires freeing memory with free() after a call.
 */
extern double* einterface_getwccellthicknesses(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of sediment cell thicknesses in the column
 *
 * Requires freeing memory with free() after a call.
 */
extern double* einterface_getsedcellthicknesses(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Short-wave radiation at water surface
 */
extern double einterface_getlighttop(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of porosity values for sediment cells in the column
 */
extern double* einterface_getporosity(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Erosion rate for the column
 */
extern double einterface_geterosionrate(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Wave-Current Friction Velocity at the bottom of the column
 */
extern double einterface_getustrcw(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of pointers to watercolumn tracer values [k * n]
 * [(n=0,k=0), (n=1,k=0), ..., (n=ntr-1,k=0), (n=0,k=1), ..., (n=ntr-1,k=nk-1)]
 *
 * Requires freeing memory with free() after a call.
 */
extern double** einterface_getwctracers(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of pointers to sediment tracer values [k * n]
 * [(n=0,k=0), (n=1,k=0), ..., (n=ntr-1,k=0), (n=0,k=1), ..., (n=ntr-1,k=nk-1)]
 *
 * Requires freeing memory with free() after a call.
 */
extern double** einterface_getsedtracers(void* model, int b);

/* @param model Pointer to host model
 * @param b Column index
 * @return Array of pointers to epibenthic variable indices
 *
 * Requires freeing memory with free() after a call.
 */
extern double** einterface_getepivars(void* model, int b);

/* @param model Pointer to host model
 */
extern void einterface_ecologyinit(void* model, void* ecology);


/**
 * retrieve the cell area for this column
 *
 * @param hmodel - the hydrodynamic host model
 * @param b - the column index
 * @return the cell area in m2
 */
extern double einterface_cellarea(void* hmodel, int b);

/**
 * Whether we're running in transport mode
 *
 */
extern int einterface_transport_mode(void);

/**
 * Get output path from host model
 *
 */
char *einterface_get_output_path(void);

/*
 * Calculates the Zenith using the library function
 */
extern double einterface_calc_zenith(void *model, double t, int b);

/** Retrieve the coordinates for this box/column
 *
 * @param model Pointer to host model
 * @return the coordinates for this box/column as int[2]
 */
/*
 extern int* einterface_getcoordinates(void* model, int b);
*/

typedef void (*quitfntype) (const char* format, ...);

/* @return quit function; must be of void(*fn)(char*, ...) type.
 */
quitfntype einterface_getquitfn();

/*
 * Get the window number from the host
 */


/** Calculates zenith given time
 * @param model host model
 * @param t time in days since 1990, beware!!!
 * @param b column number
 * @return zenith value in radians
 */
extern int einterface_get_win_num(void *model);

extern double  einterface_gettracersvel(void* model, char* name);

double einterface_botstress(void *model, int b);

double einterface_get_windspeed(void *model, int b);

void einterface_check_default_tracers(void);

char* einterface_get2Dtracername(void* model, int i);
void einterface_get_rsr_tracers(void* model, int *rtns);
int einterface_get_num_rsr_tracers(void* model);
int einterface_get_win_num(void *model);
int einterface_tracername_exists(void* model, char*name);
int einterface_tracername_exists_epi(void* model, char*name);
int einterface_get_eco_flag(void* model, char* name);

/* Generic interface */
extern void i_set_error(void* hmodel, int col, int errorf, char *text);
extern int i_get_error(void* hmodel, int col);
extern int ginterface_is_window1(void *hmodel);

/* Optical functions */
extern int einterface_is_optical(void* model, char *name);
extern int einterface_is_optical2d(void* model, char *name);
extern int einterface_get_optical_file(void* model, char *name, char *file);
extern int einterface_get_optical2d_file(void* model, char *name, char *file);
extern int einterface_get_absorp_name(void* model, char *name, char *absorp);
extern int einterface_get_scatter_name(void* model, char *name, char *scatter);
extern int einterface_get_backscat_name(void* model, char *name, char *backscat);
extern int einterface_get_abtance_name(void* model, char *name, char *abtance);
extern int einterface_get_refltce_name(void* model, char *name, char *refltce);
extern int einterface_get_trnmiss_name(void* model, char *name, char *trnmiss);
extern int einterface_get_benreflt_name(void* model, char *name, char *benreflt);
extern int einterface_get_specresp3d_name(void* model, char *name, char *specresp3d);
extern int einterface_get_specresp2d_name(void* model, char *name, char *specresp2d);
extern int einterface_get_num_ed_tracers(void* model);
extern void einterface_get_ed_tracers(void* model, int *rtns);

#define _INTERFACE_H
#endif
