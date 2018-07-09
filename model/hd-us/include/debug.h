/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/debug.h
 *  
 *  Description:
 *  Prototypes for debug messaages.
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: debug.h 5873 2018-07-06 07:23:48Z riz008 $
 *
 */

#ifndef _DEBUG_H
#define _DEBUG_H

#include <emslogger.h>

/* debug functions and macros.
 */
void debug_init(const char *dsettings);
void debug_set_output(const char *filename);
void debug_end(void);
int debug_if(const char *dtag);
void debug_usage(void);
void dlog(const char *dtag, const char *s, ...);

#define DEBUG(dtag) ((debug && debug_if(dtag)) || is_log_enabled(LDEBUG))

#endif
