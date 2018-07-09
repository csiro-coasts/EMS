/**
 *
 *  #### ENVIRONMENTAL MODELLING SUITE (EMS)
 *
 *  \file model/hd/include/debug.h
 *
 *  \brief Prototypes for debug messaages
 *
 *  \copyright
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *
 *  $Id: emslogger.h 5834 2018-06-27 00:55:37Z riz008 $
 */

#ifndef _EMSLOG_H
#define _EMSLOG_H

#include <stdarg.h>
#include <stdio.h>



typedef int (*logfn) (int level, char *format, ...);
typedef int (*logtag) (int level, char *tag, char *format, ...);

 
/* Log levels
 */
#define LTRACE 40
#define LMAIN 0
#define LDEBUG 30
#define LWARN 20
#define LINFO 10
#define LMETRIC 50
#define LERROR -10
#define LFATAL -20
#define LPANIC -30
#define LNONE  -100

/* log functions and macros.
 */
void log_init(const char *dsettings);
void log_set_output(const char *filename);
void log_end(void);
int  log_if(const char *dtag);
void log_usage(void);
void disable_loglevel(char* tag);
void enable_loglevel(char* tag);
int  is_log_enabled(int level);
void emslog_flush();
FILE* emslog_fp(int level);
int writelog(int level, const char *tag, const char *s,va_list args);
void emslog_set_log_offset(const char str[]);

char* level2tag(int t);

extern logfn  emslog; 

extern logtag emstag;


#endif
