/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/utils.h
 *  
 *  Description:
 *  Utility functions for ecology code -- header
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: utils.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined(_UTILS_H)
#include <stdio.h>

#include "declarations.h"
#include <time_utils.h>

/** Wrapper to fopen(). Exits with error message if fails.
 * @param fname File name to open.
 * @param mode Mode -- see fopen()
 * @return File handle.
 */
FILE* e_fopen(const char* fname, const char* mode);

/*
 * BY DEFAULT THESE ARE INLINED
 */
#define EPS_MAX 1.0e-10
#define X_MAX -100.0
#ifndef DEBUG
#define e_max(a) \
  ( (a > X_MAX) ? (a + sqrt(a * a + EPS_MAX)) / 2.0 : -EPS_MAX / 4.0 / a )
#define e_min(a, b) \
  ( sqrt(sqrt(1.0 / (1.0 / (a*a*a*a) + 1.0 / (b*b*b*b)))) )
#else
/** Smoothed version of max(a, 0).
 * @param a Argument
 * @return Function
 */
double e_max(double a);

/** Smoothed version of min().
 * @param a First argument
 * @param b Second Argument
 * @return Function
 */
double e_min(double a, double b);
#endif /* ifndef DEBUG */

/** Quit function for ecology.
 * @param format Quit message format
 * @param ... Quit message, arguments
 */
void e_quit(const char* format, ...);

/** Warn function for ecology.
 * @param format Warning format
 * @param ... Warning message, arguments
 */
void e_warn(const char* format, ...);

/** Wrapper to stringtable_findstringindex(). Exits with error message
 * if the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string
 */
int find_index(stringtable* st, char* s, ecology* e);

/** Finds string index in a stringtable; adds a new entry if not found.
 * Sorts the stringtable after adding a new entry.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string
 */
int find_index_or_add(stringtable* st, char* s, ecology* e);

/** Finds string index in a stringtable; adds a new entry if not found.
 * Sorts the stringtable after adding a new entry.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @param n Size associated with this entry
 * @return Index associated with the string
 */
int find_index_or_add_vector(stringtable* st, char* s, ecology* e, int n);

/** This is modified skipToKeyEnd() from `sjwlib'. Positions file stream just
 * after the key string.
 * @param fname Parameter file name (for error messaging)
 * @param fp File handle
 * @param key Key
 * @return 1 on success, 0 otherwhise
 */
int find_key(char* fname, FILE* fp, char* key);

/** Find a process of a given type.
 * @param e Pointer to ecology
 * @param pt Process type
 * @param name Process name
 * @return Process index if successful; -1 otherwise
 */
int find_process_index(ecology* e, int pt, char* name);

/** Gets string value for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter string
 */
char* get_parameter_stringvalue(ecology* e, char* s);

/** Gets parameter value for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value
 */
double get_parameter_value(ecology* e, char* s);

/** Gets the value pointer for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value
 */
double *get_parameter_value_ptr(ecology* e, char* s);

/** Gets the number of values for a parameter. Exits with error message if the
 * parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @return number of parameter values
 */
int get_parameter_num_values(ecology* e, char* s);

/** Wrapper to stringtable_findstringindex(). Unlike find_index(), goes on if
 * the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string or -1 if it wasn't found
 */
int try_index(stringtable* st, char* s, ecology* e);


/** Wrapper to stringtable_findstringindex(). Unlike find_index(),
 * it will add the tracer if it exists in the host model
 * and the string is not found.
 * @param st Pointer to stringtable
 * @param s String
 * @param ecology the ecology model
 * @return Index associated with the string or -1 if it was not found or added
 */
int try_index_or_add(stringtable* st, char* s, ecology* e);


/** Wrapper to try_index(). Unlike try_index(), it returns a
 * boolean value if the string is found.
 * @param st Pointer to stringtable
 * @param s String
 * @param verbose Verbosity
 * @return 1 if the string is found, 0 otherwise
 */
int string_exists(stringtable* st, char* s, ecology* e);


/** Gets parameter value for a parameter. Unlike get_parameter_value(),
 * goes on if the parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter value if found; NaN otherwise
 */
double try_parameter_value(ecology* e, char* s);


/** Gets parameter value for a parameter. Unlike get_parameter_value(),
 * goes on if the parameter is not found.
 * @param e Pointer to ecology
 * @param s Parameter name
 * @param verbose Verbosity
 * @return Parameter string value if found; NaN otherwise
 */
char* try_parameter_stringvalue(ecology* e, char* s);


/**
 * Whether this string contains useful information.
 * @param s String
 * @return 1 -- empty; 0 -- non-empty
 */
int isemptyline(char* s);

/**
 * Gets ratio of the given time units to second.
 * @param tunits Time units string
 * @return Ratio of these time units to second
 */
double timeunits2seconds(const char* tunits);

int get_multi_sed_index(ecology *e, int k_sed, int topk_sed, char *name);


#define _UTILS_H
#endif
