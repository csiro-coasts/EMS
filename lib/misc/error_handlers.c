/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file lib/misc/error_handlers.c
 *
 *  \brief Message log error handlers
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: error_handlers.c 5833 2018-06-27 00:21:35Z riz008 $
 */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "errfn.h"
#include "emslogger.h"


int writelog(int level, const char *tag, const char *s,va_list args);

/** Error message routine which does nothing.
  *
  * @param format printf-like format string.
  * @param ... ellipsis list of variable arguments appropriate to
  *        the format specified.
  */
void errfn_quiet(const char *format, ...)
{
}

/** Print an error message and exit the process.
  *
  * @param format printf-like format string.
  * @param ... ellipsis list of variable arguments appropriate to
  *        the format specified.
  */
void errfn_quit_default(const char *format, ...)
{
  va_list args;

  fflush(stdout);
  va_start(args, format);
  writelog(LFATAL,"",format,args);
  /*vfprintf(stderr, format, args);*/
  va_end(args);
  exit(1);
}

void errfn_set_quit(errfn fn)
{
  if (fn)
    quit = fn;
  else
    quit = errfn_quiet;
}

/** Print out a warning message to stderr.
  *
  * @param format pointer to the printf-like formatted string.
  * @param ... ellipsis variable arguments.
  */
void errfn_warn_default(const char *format, ...)
{
  
  va_list args;
  
  fflush(stdout);
  /*fprintf(stderr, "WARNING: "); */
  va_start(args, format);
  writelog(LWARN,"",format,args);
  /*vfprintf(stderr, format, args); */
  va_end(args);
  
  
}


void errfn_set_warn(errfn fn)
{
  if (fn)
    warn = fn;
  else
    warn = errfn_quiet;
}


errfn warn = errfn_warn_default;
errfn quit = errfn_quit_default;
errfn quiet = errfn_quiet;
