/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/include/errfn.h
 *
 *  \brief Header error messages
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: errfn.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _ERRFN_H
#define _ERRFN_H

typedef void (*errfn) (const char *format, ...);


extern errfn warn;
extern errfn quit;
extern errfn quiet;


void errfn_quiet(const char *format, ...);
void errfn_quit_default(const char *format, ...);
void errfn_warn_default(const char *format, ...);

void errfn_set_quit(errfn fn);

#endif
