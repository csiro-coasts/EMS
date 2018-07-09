/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: lib/ecology/include/externallibs.h
 *  
 *  Description:
 *  External library header.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: externallibs.h 5846 2018-06-29 04:14:26Z riz008 $
 *
 */

#if !defined (_EXTERNALLIBS_H)
#define _EXTERNALLIBS_H

#include <ems.h>


/*
 * error_fn
 */
typedef void (*error_fn) (const char* format, ...);

/*
 * prm_seterrorfn()
 */
typedef void (*prm_seterrfn_fn) (error_fn f);

extern prm_seterrfn_fn prm_seterrorfn;

/*
 * prm_geterrorfn()
 */
error_fn prm_geterrorfn();

/*
 * prm_readstring()
 */
typedef int (*prm_readstring_fn) (FILE* fp, char* key, char* p);

extern prm_readstring_fn prm_readstring;

/*
 * prm_readdouble()
 */
typedef int (*prm_readdouble_fn) (FILE* fp, char* key, double* v);

extern prm_readdouble_fn prm_readdouble;

/*
 * prm_readint()
 */
typedef int (*prm_readint_fn) (FILE* fp, char* key, int* v);

extern prm_readint_fn prm_readint;

/*
 * prm_getkey()
 */
typedef char* (*prm_getkey_fn) (char* keybuf, const char* prefix, const char* suffix, ...);

extern prm_getkey_fn prm_getkey;

/*
 * prm_parseline()
 */
typedef int (*prm_parseline_fn) (char *line, char **str, int max);
extern prm_parseline_fn prm_parseline;

#endif
