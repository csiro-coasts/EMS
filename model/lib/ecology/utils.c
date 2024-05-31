/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/utils.c
 *  
 *  Description:
 *  Utility functions for ecology code.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: utils.c 7576 2024-05-30 03:47:52Z riz008 $
 *
 */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <emslogger.h>
#include "ecology_internal.h"
#include "parameter_info.h"
#include "eprocess.h"
#include "utils.h"
#include "column.h"
#include "externallibs.h"
#include "einterface.h"

#define BUFSIZE 2048

/** Wrapper to fopen(). Exits with error message if fails.
 * @param fname File name to open
 * @param mode Mode -- see fopen()
 * @return File handle
 */
FILE* e_fopen(const char* fname, const char* mode)
{
    FILE* f = NULL;

    if ((f = fopen(fname, mode)) == NULL)
        e_quit(" could not open \"%s\": %s\n", fname, strerror(errno));

    return f;
}

/*
 * THESE 2 FUNCTIONS ARE ONLY DEFINED IN DEBUG MODE
 */
#ifdef DEBUG
/** Smoothed version of max(a, 0).
 * @param a Argument
 * @return Function
 */
double e_max(double a)
{
  if (a > X_MAX)
    return (a + sqrt(a * a + EPS_MAX)) / 2.0;
  else
    return -EPS_MAX / 4.0 / a;
}

/** Smoothed version of min().
 * @param a First argument
 * @param b Second Argument
 * @return Function
 */
double e_min(double a, double b)
{
    a *= a;
    a *= a;
    b *= b;
    b *= b;
    return sqrt(sqrt(1.0 / (1.0 / a + 1.0 / b)));
}
#endif /* ifdef DEBUG */

/** Quit function for ecology.
 * @param format Quit message format
 * @param ... Quit message, arguments
 */
void e_quit(const char* format, ...)
{
    va_list args;
/*
    fflush(stdout);
*/
    /*
     * I do not put here the line:
     * fprintf(stderr, "error: ecology: ");
     * for the following reason:
     * most of error messages are handled by e->quitfn, which is initialised
     * to e_quit, but can be set from outside the ecology module. Because of
     * this, all error messages have to contain explicit module identifier.
     */
    va_start(args, format);
    /* vfprintf(stderr, format, args); */
    writelog(LFATAL,"ecology:",format,args);
    va_end(args);
    exit(1);
}

/** Warn function for ecology.
 * @param format Warning format
 * @param ... Warning message, arguments
 */
void e_warn(const char* format, ...)
{
    va_list args;

    fflush(stdout);
/*
    fprintf(stderr, "warning: ecology: ");
*/
    va_start(args, format);
    /* vfprintf(stderr, format, args); */
    writelog(LWARN,"ecology:",format,args);
    va_end(args);
}

/** Wrapper to stringtable_findstringindex(). Exits with error message
 * if the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param verbose Verbosity
 * @return Index associated with the string
 */
int find_index(stringtable* st, char* s, ecology* e)
{
    int index = stringtable_findstringindex(st, s);

    if (index < 0)
        e_quit(" %s: \"%s\" not found\n", st->name, s);
    emstag(LTRACE,"eco:utils:find_index","  %s: %s(%d)\n", st->name, s, index);
    /*
    if (verbose)
        fprintf(stderr, "  %s: %s(%d)\n", st->name, s, index);
        */
    return index;
}

/** Finds string index in a stringtable; adds a new entry if not found.
 * Sorts the stringtable after adding a new entry
 * @param st Pointer to stringtable
 * @param s String
 * @param n Size associated with this entry
 * @return Index associated with the string
 * @sa find_index_or_add
 */
int find_index_or_add_vector(stringtable* st, char* s, ecology* e, int n)
{
    int index = stringtable_findstringindex(st, s);

    if (index < 0) {
      // Not found so add
      stringtable_add(st, s, -1);
      index = st->n - 1;
      st->se[index]->naccess = 1;
      st->se[index]->n = n;
      // FR: I know sorting speeds up the lookup but the indicies that
      //     are already cached will be wrong!!
      // stringtable_sort(st);
    }
    emstag(LTRACE,"eco:utils:find_index_or_add", " %s: %s(%d)\n", st->name, s, index);
    return index;
}

/** Finds string index in a stringtable; adds a new entry if not found.
 * Sorts the stringtable after adding a new entry.
 * @param st Pointer to stringtable
 * @param s String
 * @return Index associated with the string
 * @sa find_index_or_add_vector
 */
int find_index_or_add(stringtable* st, char* s, ecology* e)
{
  return(find_index_or_add_vector(st, s, e, 1));
}

/** This is modified skipToKeyEnd() from `sjwlib'. Positions file stream just
 * after the key string.
 * @param fname Parameter file name (for error messaging)
 * @param fp File handle
 * @param key Key
 * @return 1 on success, 0 otherwhile
 */
int find_key(char* fname, FILE* fp, char* key)
{
    char buf[BUFSIZE];
    int len = strlen(key);
    int fpos;
    char* s;
    int rewound = 0;

    do {
        fpos = ftell(fp);
        s = fgets(buf, BUFSIZE, fp);
        if (s == NULL) {
	  if (!rewound) {
	    fpos = 0;
	    if (fseek(fp, fpos, 0))
	      e_quit(" \"%s\": could not rewind: %s\n", fname, strerror(errno));
	    s = fgets(buf, BUFSIZE, fp);
	    rewound = 1;
	  } else
	    return 0;
        }
	
        if (s == NULL)
	  break;

        while (isspace((int)s[0]))
	  s++;
	
	// Can't possibly match
	if (strlen(s) < len)
	  continue;

	if (isspace((int)s[len]) || /* to parse processes */ s[len] == '{')
	  s[len] = 0;
	
    } while (strcasecmp(key, s) != 0);

    if (s == NULL)
        return 0;

    if (fseek(fp, fpos + (s - buf) + len, 0))
        e_quit(" \"%s\": no characters found after \"%s\"\n", fname, key);

    return 1;
}

/** Find a process of a given type.
 * @param e Pointer to ecology
 * @param pt Process type
 * @param name Process name
 * @return Process index if successful; -1 otherwise
 */
int find_eprocess_index(ecology* e, int pt, char* name)
{
    int i;

    for (i = 0; i < e->npr[pt]; ++i)
        if (strcasecmp(e->processes[pt][i]->name, name) == 0)
            return i;
    return -1;
}

/** Gets string value for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter string
 */
char* get_parameter_stringvalue(ecology* e, char* s)
{
    int index = e->find_index(e->prms, s, e);
    char* str = e->pinfo[index].stringvalue;

    if (str[0] == 0)
        e->quitfn(" \"%s\": empty string after \"%s.value\"\n", e->prmfname, s);
    return str;
}

/** Gets value for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value
 */
double get_parameter_value(ecology* e, char* s)
{
  double *ptr = get_parameter_value_ptr(e, s);
  emstag(LDEBUG,"Ecological parameter read, %s = %e \n",s,*ptr);
  eco_write_setup(e,"Ecol. parameter forced read, %s = %e \n",s,*ptr);
  int index = find_index(e->prms, s, e);
  
  // Output parameter values to a netcdf file.

  if (ginterface_is_window1(e->model)) {
    int var_exists = 0;
    int varid;
    int ncid = e->eco_params_ncid;
    
    var_exists = nc_inq_varid(ncid,s,&varid);
    
    if (var_exists != 0){
      ncredef(ncid);
      nc_def_var(ncid,s,NC_DOUBLE,0,0,&varid);
      nc_put_att_text(ncid, varid, "name",strlen(s),s);
      char* ttz = e->pinfo[index].desc;
      nc_put_att_text(ncid, varid, "description",strlen(ttz),ttz);
      char* ttt = e->pinfo[index].units;
      nc_put_att_text(ncid, varid, "units",strlen(ttt),ttt);
      nc_put_att_text(ncid, varid, "function_call",13,"get_parameter");
      nc_enddef(ncid);
      ncw_put_var_double(ECO_PARAMS_SETUP,ncid,varid,ptr);
    }
  }
  
  return(*ptr);
}

/** Gets the value pointer for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value
 */
double *get_parameter_value_ptr(ecology* e, char* s)
{
    int index = find_index(e->prms, s, e);

    return e->pinfo[index].value;
}

/** Gets the number of values for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @return number of parameter values
 */
int get_parameter_num_values(ecology* e, char* s)
{
    int index = find_index(e->prms, s, e);

    return e->pinfo[index].num_values;
}

/** Wrapper to stringtable_findstringindex(). Unlike find_index(), goes on if
 * the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string or -1 if it wasn't found
 */
int try_index(stringtable* st, char* s, ecology* e)
{
    int index = stringtable_findstringindex(st, s);

    emstag(LTRACE,"eco:utils:try_index", "  %s: %s(%d)\n", st->name, s, index);
    return index;
}

/** Wrapper to stringtable_findstringindex(). Unlike find_index(),
 * it will add the tracer if it exists in the host model
 * and the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string or -1 if it was not found or added
 */
int try_index_or_add(stringtable* st, char* s, ecology* e)
{
    int index = stringtable_findstringindex(st, s);

    if(index < 0) {
      if (e->pre_build)
	return(-1);

      if ( (strcmp(st->name, "epivariables") == 0 &&
	    einterface_tracername_exists_epi(e->model, s)) ||
	   einterface_tracername_exists(e->model,s)) {
	stringtable_add(st, s, -1);
	index = st->n - 1;
	st->se[index]->naccess = 1;
	// FR: I know sorting speeds up the lookup but the indicies that
	//     are already cached will be wrong!!
	// stringtable_sort(st);
      } else
	emstag(LWARN,"eco:utils:try_index_or_add", " Doesn't exist:  %s: %s", st->name, s);
    }
    emstag(LDEBUG,"eco:utils:try_index_or_add", "  %s: %s(%d)", st->name, s, index);
    return index;
}

/*
 * find_index for multiple sediment layers
 */
int get_multi_sed_index(ecology *e, int k_sed, int topk_sed, char *name)
{
  int  index = e->find_index(e->tracers, name, e);
  
  /*
   * The y vector looks like wc + sed(top) + epi + sed(top-1 to bottom)
   */
  if (k_sed < topk_sed)
    return(e->ntr*2 + e->nepi + e->ntr*(topk_sed - k_sed-1) + index);
  else
    return(e->ntr + index);
}

/*
 * Wrapper to handle y values in any sediment layer
 * Make this a macro
 */
double get_multi_sed_y(column* col, double *y, int sed_index, int sed_layer)
{
  ecology *e = col->e;
  
  /* Not the first - access from the end */
  if (sed_layer < col->topk_sed) 
    return(y[e->ntr*2 + e->nepi + e->ntr*(col->topk_sed-sed_layer-1) + sed_index]);
  /* The top sed layer */
  return(y[e->ntr + sed_index]);
}

/** Wrapper to try_index(). Unlike try_index(), it returns a
 * boolean value if the string is found.
 * @param st Pointer to stringtable
 * @param s String
 * @param verbose Verbosity
 * @return 1 if the string is found, 0 otherwise
 */
int string_exists(stringtable* st, char* s, ecology* e)
{
    int index = stringtable_findstringindex(st, s);

    emstag(LTRACE,"eco:utils:string_exists", "  %s: %s(%d)\n", st->name, s, index);
    return (index > -1 );
}


/** Gets parameter value for a parameter. Unlike get_parameter_value(),
 * goes on if the parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value if found; NaN otherwise
 */
double try_parameter_value(ecology* e, char* s)
{
    int index = try_index(e->prms, s, e);
    int var_exists = 0;
    
    if (index < 0){
      eco_write_setup(e,"Ecol. parameter tried to read %s, but not in parameter file \n",s);
      return NaN;
    }else{
      eco_write_setup(e,"Ecol. parameter tried to read %s and found %e \n",s,e->pinfo[index].value[0]);
      if (ginterface_is_window1(e->model)) {
	int ncid = e->eco_params_ncid,varid;
	var_exists = nc_inq_varid(ncid,s,&varid);
	
	if (var_exists != 0){
	  ncredef(ncid);
	  nc_def_var(ncid,s,NC_DOUBLE,0,0,&varid);
	  char* tty = e->pinfo[index].name;
	  nc_put_att_text(ncid, varid, "name",strlen(tty),tty);
	  char* ttz = e->pinfo[index].desc;
	  nc_put_att_text(ncid, varid, "description",strlen(ttz),ttz);
	  char* ttt = e->pinfo[index].units;
	  nc_put_att_text(ncid, varid, "units",strlen(ttt),ttt);
	  nc_put_att_text(ncid, varid, "function_call",13,"try_parameter");
	  nc_enddef(ncid);
	  ncw_put_var_double(ECO_PARAMS_SETUP,ncid,varid,&e->pinfo[index].value[0]);
	}
      }
      return e->pinfo[index].value[0];
    }
}

/** Gets parameter value for a parameter. Unlike get_parameter_value(),
 * goes on if the parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter string value if found; NaN otherwise
 */
char* try_parameter_stringvalue(ecology* e, char* s)
{
    int index = try_index(e->prms, s, e);
    char* str;
    if(index < 0)
    	return 0;

    str = e->pinfo[index].stringvalue;
    return str;
}


/**
 * Checks whether this string contains useful information.
 * @param s String
 * @return 1 -- empty; 0 -- non-empty
 */
int isemptyline(char* s)
{
    if (s == NULL)
        return 1;
    if (s[0] == '#')
        return 1;
    while (s[0] != 0 && isspace((int)s[0]))
        s++;
    if (s[0] == 0)
        return 1;

    return 0;
}

static struct {
    char* unit;
    double mult;
} sec_conv[] = {
    {"us", 0.000001},
    {"usec", 0.000001},
    {"ms", 0.001},
    {"msec", 0.001},
    {"sec", 1.0},
    {"second", 1.0},
    {"min", 60.0},
    {"minute", 60.0},
    {"hr", 3600.0},
    {"hour", 3600.0},
    {"day", 86400.0},
    {"week", 604800.0},
    {"d", 86400.0},
    {"h", 3600.0},
    {"s", 1.0},
    {NULL, 0.0}
};

/**
 * Gets ratio of the given time units to second.
 * @param tunits Time units string
 * @return Ratio of these time units to second
 */
double timeunits2seconds(const char* tunits)
{
    int i;

    if (tunits == NULL)
        return 1.0;

    for (i = 0; sec_conv[i].unit != NULL; i++)
        if (strstr(tunits, sec_conv[i].unit) == tunits)
             return sec_conv[i].mult;

    return 1.0;
}
